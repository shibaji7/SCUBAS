# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Sphinx configuration for the SCUBAS documentation (the BOLTS stack: autodoc, napoleon, viewcode, intersphinx, MyST Markdown, Read the Docs theme).

Built with the ``dirhtml`` builder, so every page keeps the address it had under mkdocs (``tutorial/theory/``, ``dev/cables/``, ...).
"""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))  # the flat package is importable without installing it

project = "SCUBAS"
author = "Shibaji Chakraborty"
copyright = "2023-2026, Shibaji Chakraborty"  # noqa: A001

# The version comes from pyproject.toml, so the docs always match the release. No fallback: a missing version must fail the build.
release = tomllib.loads((_ROOT / "pyproject.toml").read_text())["project"]["version"]
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",  # pull docstrings in
    "sphinx.ext.autosummary",  # summary tables
    "sphinx.ext.napoleon",  # NumPy and Google docstring sections (the package uses NumPy style)
    "sphinx.ext.intersphinx",  # links into Python, NumPy, SciPy, pandas and matplotlib
    "sphinx.ext.mathjax",  # equations
    "sphinx.ext.viewcode",  # [source] link on every entry
    "sphinx_autodoc_typehints",  # types in the parameter descriptions
    "sphinx_copybutton",  # copy button on code blocks
    "myst_parser",  # the pages are Markdown
]

autosummary_generate = False
autodoc_default_options = {"show-inheritance": True, "member-order": "bysource"}  # members are listed per directive
autodoc_typehints = "description"
autoclass_content = "both"  # the class docstring and the __init__ docstring (the legacy classes document their parameters there)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_rtype = False
napoleon_use_param = True
default_role = "code"  # a single backtick is code, as in the Markdown pages

myst_enable_extensions = ["colon_fence", "deflist", "amsmath", "dollarmath", "attrs_inline"]
myst_heading_anchors = 3

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "requirements.txt"]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
master_doc = "index"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False, "sticky_navigation": True}
html_context = {
    "display_github": True,
    "github_user": "shibaji7",
    "github_repo": "SCUBAS",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

nitpicky = False
suppress_warnings = [
    "myst.header",  # the tutorials start at a deep heading level
    "misc.highlighting_failure",  # a tutorial code block that the lexer takes in relaxed mode (it is still highlighted)
]
html_css_files = ["css/scubas.css"]

BANNER = """:::{warning}
**Development status: alpha.** SCUBAS is research software (version 0.x). Its method has been compared with observations in published studies
(see [how to cite SCUBAS](/user/citing.md)), but results for other cables, regions or storms are not guaranteed, and the API and the numbers
may change in future versions. Check anything you rely on against the cited papers and your own data.
:::

"""


STATUS = re.compile(r"\*\*(done|planned|proposed)\*\*")


def add_banner(app, docname, source):  # noqa: ANN001, ANN201, ARG001
    """Put the development-status banner at the top of every page (``source-read``): one place for the wording, and it also covers the pages
    that are included from root files (the changelog and the roadmap). Pages that start with front matter keep it first.
    On the roadmap page the status words (``**done**``, ``**planned**``, ``**proposed**``) get a colour class: done is green."""
    text = source[0]
    if docname == "roadmap":  # the page is the root ROADMAP.md itself (an include would be expanded after this event)
        roadmap = _ROOT / "ROADMAP.md"
        app.env.note_dependency(str(roadmap))
        text = STATUS.sub(r"[**\1**]{.status-\1}", roadmap.read_text())
    if text.startswith("---\n"):
        end = text.index("\n---\n", 4) + len("\n---\n")
        source[0] = text[:end] + "\n" + BANNER + text[end:]
    else:
        source[0] = BANNER + text


def setup(app):  # noqa: ANN001, ANN201
    """Register the banner."""
    app.connect("source-read", add_banner)
    return {"parallel_read_safe": True}
