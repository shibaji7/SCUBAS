# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023-2026 Shibaji Chakraborty
"""The pattern of planned capabilities: ``scubas.future`` holds only stubs that raise ``NotImplementedYetError``, and none of them also exists in the real package."""

import ast
import importlib
import pathlib

import pytest

import scubas
from scubas.errors import NotImplementedYetError, ScubasError

PKG = pathlib.Path(scubas.__file__).parent
FUTURE = PKG / "future"


def problems(future_dir: pathlib.Path, package_dir: pathlib.Path) -> list[str]:
    """What breaks the rules: a stub module that does not raise ``NotImplementedYetError``, or that already exists in the real package."""
    found = []
    for f in sorted(future_dir.rglob("*.py")):
        if f.name == "__init__.py":
            continue
        rel = f.relative_to(future_dir)
        if "NotImplementedYetError" not in f.read_text():
            found.append(
                f"{rel} is in future but does not raise NotImplementedYetError"
            )
        if (package_dir / rel).exists():
            found.append(
                f"{rel} exists in the real package as well: remove the stub or the module"
            )
    return found


def test_the_error_names_the_roadmap_item_and_is_a_not_implemented_error():
    err = NotImplementedYetError(
        "frequency-dependent cable properties", "ROADMAP: cable"
    )
    assert (
        str(err)
        == "frequency-dependent cable properties is not implemented yet (ROADMAP: cable)."
    )
    assert isinstance(err, NotImplementedError) and isinstance(err, ScubasError)
    assert "see ROADMAP.md" in str(NotImplementedYetError("x"))


def test_the_future_package_follows_the_rules_today():
    assert problems(FUTURE, PKG) == []


def test_the_checker_has_teeth(tmp_path):
    (tmp_path / "future" / "cable").mkdir(parents=True)
    (tmp_path / "real" / "cable").mkdir(parents=True)
    (tmp_path / "future" / "cable" / "working.py").write_text(
        "def run():\n    return 1\n"
    )
    (tmp_path / "future" / "cable" / "both.py").write_text(
        "from scubas.errors import NotImplementedYetError\n"
    )
    (tmp_path / "real" / "cable" / "both.py").write_text("x = 1\n")
    (tmp_path / "future" / "cable" / "fine.py").write_text(
        "from scubas.errors import NotImplementedYetError\n"
    )
    got = problems(tmp_path / "future", tmp_path / "real")
    assert any("working.py" in p and "does not raise" in p for p in got)
    assert any("both.py" in p and "exists in the real package" in p for p in got)
    assert not any("fine.py" in p for p in got)


def imported_modules(path: pathlib.Path) -> list[str]:
    """The absolute module names a file imports (relative imports are resolved against the package of the file)."""
    package = ".".join(["scubas", *path.relative_to(PKG).parent.parts])
    names = []
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            base = package.rsplit(".", node.level - 1)[0] if node.level else ""
            names += [".".join(x for x in (base, node.module or "") if x)] + [
                f"{base}.{a.name}" if node.level else f"{node.module}.{a.name}"
                for a in node.names
            ]
    return names


def test_the_future_package_has_an_introduction():
    assert ast.get_docstring(ast.parse((FUTURE / "__init__.py").read_text()))


def test_the_real_code_never_imports_a_stub():
    for f in PKG.rglob("*.py"):
        if "future" in f.relative_to(PKG).parts:
            continue
        assert not [
            n
            for n in imported_modules(f)
            if n == "scubas.future" or n.startswith("scubas.future.")
        ], f


@pytest.mark.parametrize("name", ["scubas.future", "scubas.errors"])
def test_the_new_modules_import_without_side_effects(name, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    importlib.import_module(name)
    assert list(tmp_path.iterdir()) == []
