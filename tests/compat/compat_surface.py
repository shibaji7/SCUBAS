# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""Collect the public (and documented-private) surface of the ``scubas`` modules as plain JSON-able data.

The snapshot written from the 0.1.9 code is ``baseline/api_surface.json``. ``test_api_surface.py`` fails when a
baseline name disappears or when anything recorded about it changes (signature, kind, base classes, dataclass
fields, ``PROFILES`` values). A new name is allowed: additions are compatible, removals and changes are not.
"""

from __future__ import annotations

import ast
import dataclasses
import importlib
import inspect
import math
from pathlib import Path
from typing import Any

import numpy as np

MODULES = [
    "scubas",
    "scubas.cables",
    "scubas.conductivity",
    "scubas.datasets",
    "scubas.models",
    "scubas.plotlib",
    "scubas.utils",
    "scubas.sources",
    "scubas.sources.current_sources",
]
# methods that are part of the contract even though they are dunder
KEPT_DUNDERS = {"__init__", "__post_init__"}


def _norm(value: Any) -> Any:
    """A stable, JSON-able description of a default value."""
    if value is inspect.Parameter.empty:
        return "<empty>"
    if isinstance(value, np.ndarray):
        return f"ndarray{value.shape}:{_digest(value)}"
    if isinstance(value, float):
        return "inf" if math.isinf(value) else repr(value)
    if isinstance(value, (list, tuple)):
        return [_norm(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _norm(v) for k, v in value.items()}
    if callable(value):
        return f"<callable {getattr(value, '__name__', type(value).__name__)}>"
    return repr(value)


def _digest(arr: np.ndarray) -> str:
    import hashlib

    return hashlib.sha1(np.ascontiguousarray(arr, dtype=float).tobytes()).hexdigest()[
        :12
    ]


def signature_of(obj: Any) -> dict[str, Any] | None:
    """Parameter names, kinds and defaults (annotations are left out on purpose: adding type hints is compatible)."""
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        return None
    return {
        "params": [
            [p.name, p.kind.name, _norm(p.default)] for p in sig.parameters.values()
        ],
    }


def _kind(cls: type, name: str) -> str:
    raw = inspect.getattr_static(cls, name)
    for label, typ in (
        ("staticmethod", staticmethod),
        ("classmethod", classmethod),
        ("property", property),
    ):
        if isinstance(raw, typ):
            return label
    return "function" if inspect.isfunction(raw) else type(raw).__name__


def _class_surface(cls: type, node: ast.ClassDef) -> dict[str, Any]:
    out: dict[str, Any] = {"bases": [b.__name__ for b in cls.__bases__], "methods": {}}
    for item in node.body:
        if not isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        name = item.name
        if name.startswith("__") and name not in KEPT_DUNDERS:
            continue
        kind = _kind(cls, name)
        target = getattr(cls, name)
        if kind == "property":
            sig = None
        else:
            sig = signature_of(target)
        out["methods"][name] = {"kind": kind, "signature": sig}
    if dataclasses.is_dataclass(cls):
        out["dataclass"] = {
            "frozen": bool(cls.__dataclass_params__.frozen),
            "fields": [
                [
                    f.name,
                    "<empty>" if f.default is dataclasses.MISSING else _norm(f.default),
                ]
                for f in dataclasses.fields(cls)
            ],
        }
    if hasattr(cls, "__members__"):  # Enum
        out["members"] = {k: _norm(v.value) for k, v in cls.__members__.items()}
    return out


def profiles_surface() -> dict[str, Any]:
    """Every layer of every ``scubas.datasets.PROFILES`` entry (these are data users rely on)."""
    from scubas.datasets import PROFILES

    out = {}
    for key, site in vars(PROFILES).items():
        out[key] = {
            "name": site.name,
            "description": site.description,
            "layers": [
                [
                    layer.name,
                    _norm(float(layer.thickness)),
                    _norm(float(layer.conductivity)),
                ]
                for layer in site.layers
            ],
        }
    return out


def module_surface(modname: str) -> dict[str, Any]:
    """Names defined in one module: functions, classes (with methods) and top-level assignments."""
    mod = importlib.import_module(modname)
    path = Path(mod.__file__)
    tree = ast.parse(path.read_text())
    out: dict[str, Any] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out[node.name] = {
                "kind": "function",
                "signature": signature_of(getattr(mod, node.name)),
            }
        elif isinstance(node, ast.ClassDef):
            out[node.name] = {
                "kind": "class",
                **_class_surface(getattr(mod, node.name), node),
            }
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and not target.id.startswith("__"):
                    out[target.id] = {
                        "kind": "assignment",
                        "type": type(getattr(mod, target.id)).__name__,
                    }
    return out


def collect_surface() -> dict[str, Any]:
    """The whole snapshot."""
    return {
        "modules": {m: module_surface(m) for m in MODULES},
        "profiles": profiles_surface(),
    }
