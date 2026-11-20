# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Parameters with provenance: each number carries its unit, status and source.

This is the pattern for **new code**. The constants of the existing modules are not
rewired to it; they are listed, read-only, in ``defaults/legacy_constants.toml`` (see
:func:`legacy_constants`) and a test keeps that list equal to what the code uses.

Statuses
    sourced   taken from a document (the source string names it)
    assumed   a placeholder, or a legacy default whose origin is not documented
    unknown   ``TBD``: asking for it raises ``MissingParameterError``
    override  changed in code or in a user file for one run

A parameter is a TOML table ``{value, unit, status, source}``. A file shipped with the
package must give all four; a user file may give a bare value, which becomes an
``override``. The ``unit`` is a label in this first version: nothing is converted.
"""

from __future__ import annotations

import copy
import hashlib
import json
import tomllib
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

from .errors import ConfigError, MissingParameterError

TBD = "TBD"
DEFAULTS_DIR = Path(__file__).with_name("defaults")


class Status(StrEnum):
    """Where a number stands (see the module docstring for what each value means)."""

    SOURCED = "sourced"
    ASSUMED = "assumed"
    UNKNOWN = "unknown"
    OVERRIDE = "override"


@dataclass(frozen=True)
class Param:
    """One parameter with its unit, status and source.

    Frozen: a changed value is a new ``Param`` with status ``override``.
    """

    section: str
    name: str
    value: Any
    unit: str
    status: Status
    source: str

    @property
    def key(self) -> str:
        """``section.name``, as used in overrides and messages."""
        return f"{self.section}.{self.name}"

    @property
    def is_unknown(self) -> bool:
        """True while the value is the ``TBD`` marker."""
        return isinstance(self.value, str) and self.value == TBD

    def require(self) -> Any:
        """Return the value, or raise :class:`~scubas.errors.MissingParameterError`.

        The error names the parameter and the reason it is unknown.
        """
        if self.is_unknown:
            raise MissingParameterError(f"{self.key} is unknown (TBD): {self.source}")
        return self.value


class ParameterSet:
    """All parameters, grouped in sections.

    Changes go through :meth:`with_overrides`, which returns a copy.
    """

    def __init__(self, data: dict[str, dict[str, Param]]) -> None:
        self._data = data

    def param(self, section: str, name: str) -> Param:
        """Return the full :class:`Param`; use :meth:`get` for just the value."""
        try:
            return self._data[section][name]
        except KeyError as e:
            raise ConfigError(f"no parameter {section}.{name}") from e

    def get(self, section: str, name: str) -> Any:
        """Return the value; raise ``MissingParameterError`` while it is ``TBD``."""
        return self.param(section, name).require()

    def has(self, section: str, name: str) -> bool:
        """Return True if the parameter exists and has a value (is not ``TBD``)."""
        try:
            return not self.param(section, name).is_unknown
        except ConfigError:
            return False

    def params(self) -> list[Param]:
        """Return every parameter, flat, in file order."""
        return [p for sec in self._data.values() for p in sec.values()]

    def sections(self) -> list[str]:
        """Return the section names, in file order."""
        return list(self._data)

    def with_overrides(self, overrides: dict[str, Any]) -> ParameterSet:
        """Return a copy with ``{"section.name": value}`` replaced (as overrides)."""
        data = copy.deepcopy(
            self._data
        )  # never mutate the original: other runs may hold it
        for key, value in overrides.items():
            section, _, name = key.partition(".")
            if not name or section not in data or name not in data[section]:
                raise ConfigError(
                    f"no parameter {key} to override (use 'section.name')"
                )
            old = data[section][name]
            data[section][name] = Param(
                section,
                name,
                value,
                old.unit,
                Status.OVERRIDE,
                f"override of: {old.source}",
            )
        return ParameterSet(data)

    def table_md(self) -> str:
        """Return every parameter as a Markdown table (the docs use it)."""
        rows = ["| key | value | unit | status | source |", "|---|---|---|---|---|"]
        rows += [
            f"| `{p.key}` | {p.value} | {p.unit} | {p.status.value} | {p.source} |"
            for p in self.params()
        ]
        return "\n".join(rows)

    def counts(self) -> dict[str, int]:
        """Return how many parameters have each status."""
        out = {s.value: 0 for s in Status}
        for p in self.params():
            out[p.status.value] += 1
        return out

    def fingerprint(self) -> str:
        """Return a short hash of every value and status (traces a result to inputs)."""
        blob = json.dumps(
            [(p.key, p.value, p.status.value) for p in self.params()],
            sort_keys=True,
            default=str,
        )
        return hashlib.sha256(blob.encode()).hexdigest()[:12]


def _parse_entry(section: str, name: str, raw: Any, *, from_user: bool) -> Param:
    """Turn one TOML entry into a :class:`Param` (shipped files: complete tables)."""
    where = f"{section}.{name}"
    if not isinstance(raw, dict):
        if not from_user:
            raise ConfigError(
                f"{where}: a shipped parameter needs value, unit, status and source"
            )
        return Param(section, name, raw, "", Status.OVERRIDE, "user file")
    if "value" not in raw:
        raise ConfigError(f"{where}: a table needs a 'value'")
    status = raw.get("status", "override" if from_user else None)
    if status is None:
        raise ConfigError(f"{where}: missing 'status'")
    if "source" not in raw and not from_user:
        raise ConfigError(f"{where}: missing 'source'")
    if "unit" not in raw and not from_user:
        raise ConfigError(
            f"{where}: missing 'unit' (use \"\" for a quantity without one)"
        )
    try:
        st = Status(status)
    except ValueError as e:
        raise ConfigError(
            f"{where}: status {status!r} not in {[s.value for s in Status]}"
        ) from e
    if st is Status.UNKNOWN and raw["value"] != TBD:
        raise ConfigError(f"{where}: status 'unknown' needs the value {TBD!r}")
    return Param(
        section,
        name,
        raw["value"],
        raw.get("unit", ""),
        st,
        raw.get("source", "user file"),
    )


def parse_parameters(data: dict[str, Any], *, from_user: bool = False) -> ParameterSet:
    """Build a :class:`ParameterSet` from a parsed TOML document."""
    out: dict[str, dict[str, Param]] = {}
    for section, entries in data.items():
        if not isinstance(entries, dict):
            raise ConfigError(f"{section}: a section must be a table of parameters")
        out[section] = {
            n: _parse_entry(section, n, raw, from_user=from_user)
            for n, raw in entries.items()
        }
    return ParameterSet(out)


def read_toml(path: str | Path) -> dict[str, Any]:
    """Parse a TOML file; a missing file raises ``ConfigError``."""
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"parameter file {p} not found")
    with p.open("rb") as f:
        return tomllib.load(f)


def load_parameters(
    defaults: str | Path,
    path: str | Path | None = None,
    overrides: dict[str, Any] | None = None,
) -> ParameterSet:
    """Load a defaults file, then an optional user file, then explicit overrides.

    Parameters
    ----------
    defaults
        A TOML file whose every parameter is a complete table.
    path
        A user TOML file merged over the defaults. A key that is not in the defaults
        is an error (a typo, not a new setting); a bare value becomes an ``override``.
    overrides
        ``{"section.name": value}`` applied last.

    Returns
    -------
    ParameterSet
        The merged parameters.
    """
    base = parse_parameters(read_toml(defaults))
    if path is not None:
        user = parse_parameters(read_toml(path), from_user=True)
        for section, entries in user._data.items():
            for name, new in entries.items():
                if name not in base._data.get(section, {}):
                    raise ConfigError(
                        f"{Path(path).name}: unknown parameter {section}.{name}"
                    )
                old = base._data[section][name]
                base._data[section][name] = Param(
                    section,
                    name,
                    new.value,
                    new.unit or old.unit,
                    new.status,
                    new.source,
                )
    return base.with_overrides(overrides) if overrides else base


def legacy_constants() -> ParameterSet:
    """Return the constants of the existing modules with their provenance.

    The list is read-only (ADR 0010).

    Returns
    -------
    ParameterSet
        ``defaults/legacy_constants.toml``. A test checks every value against the code.
    """
    return parse_parameters(read_toml(DEFAULTS_DIR / "legacy_constants.toml"))


__all__ = [
    "DEFAULTS_DIR",
    "TBD",
    "Param",
    "ParameterSet",
    "Status",
    "legacy_constants",
    "load_parameters",
    "parse_parameters",
    "read_toml",
]
