# Architecture decision records

Short records of decisions that shape the project, so the reason survives the people. Each has a status (accepted, proposed,
superseded) and the alternatives that were rejected.

| Record | Decision | Status |
|---|---|---|
| [0001](0001-legacy-numerics-are-pinned.md) | The numbers, API and file formats of 0.1.9 are pinned, not rewritten | accepted |
| [0002](0002-stay-on-mkdocs-material.md) | The documentation stays on mkdocs-material | superseded by 0007 |
| [0003](0003-flat-layout-and-runtime-dependencies.md) | Keep the flat package layout and the runtime dependencies in 0.x | accepted |
| [0004](0004-keep-loguru.md) | Keep loguru for logging | accepted |
| [0005](0005-known-defects-are-pinned.md) | Known defects are pinned and reported, never fixed silently | accepted |
| [0006](0006-license-lgpl.md) | Relicense from MIT to LGPL-3.0-or-later | accepted |
| [0007](0007-documentation-on-sphinx.md) | The documentation moves to Sphinx and MyST, built as `dirhtml` | accepted |
| [0008](0008-planned-capabilities-are-stubs.md) | A planned capability is a stub that says so | accepted |

```{toctree}
:hidden:
:maxdepth: 1

0001-legacy-numerics-are-pinned
0002-stay-on-mkdocs-material
0003-flat-layout-and-runtime-dependencies
0004-keep-loguru
0005-known-defects-are-pinned
0006-license-lgpl
0007-documentation-on-sphinx
0008-planned-capabilities-are-stubs
```
