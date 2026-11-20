# Planned capabilities: `scubas.future`

A capability that is planned but not built is a **stub that says so**, so the code, the roadmap and the documentation cannot drift apart and nobody gets a silent wrong
answer from half-written code. The pattern is the one BOLTS uses (`bolts.future`); the reason is in [ADR 0008](adr/0008-planned-capabilities-are-stubs.md).

## What a stub is

A module under `scubas/future/` that mirrors the package the capability will move to, whose public functions only raise `NotImplementedYetError`:

```python
# scubas/future/cable/frequency_dependent.py
"""Frequency-dependent cable properties: skin effect in the conductor. Planned for the roadmap item "..."; see ROADMAP.md."""

from __future__ import annotations

from typing import Any, Never

from ...errors import NotImplementedYetError


def run(*args: Any, **kwargs: Any) -> Never:
    """Not implemented yet: raises ``NotImplementedYetError`` naming the roadmap item."""
    raise NotImplementedYetError("frequency-dependent cable properties", "ROADMAP: frequency-dependent cable")


__all__ = ["run"]
```

The signature is `*args, **kwargs` on purpose: the real signature is decided when the capability is built.

## Rules

- **Only maintainers add one**, for a capability that is on the [roadmap](../roadmap.md); the stub and the roadmap row are added in the same pull request. The package is empty
  until a capability is named. Nothing is invented to fill it.
- The stub never exists in the real package at the same time (`tests/test_future.py` checks it): when the capability is implemented, the module moves out of `scubas.future` in
  the same pull request that builds it and removes its roadmap row.
- A stub is never imported by the real code and is not part of the public API of 0.1.9 (the compatibility suite pins that).
- The message names the roadmap item, so a user who hits it knows where to look.

## Where it shows up

The roadmap lists every slot with what it will do; this page and the API reference describe the pattern; `NotImplementedYetError` is a `NotImplementedError`, so code that
catches the built-in one also catches it.
