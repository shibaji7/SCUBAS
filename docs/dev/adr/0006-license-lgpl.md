# ADR 0006: relicense from MIT to LGPL-3.0-or-later

**Status:** proposed (draft pull request; waiting for the consent of every copyright holder).

**Context.** The owner decided that SCUBAS should move to LGPL 3, the license of the sibling project BOLTS. SCUBAS 0.1.1 to 0.1.9 were published
under the MIT license, with contributions credited to other people (a co-author on the documentation pages, a credited author of
`scubas/conductivity.py`).

**Decision.** From the first release after 0.1.9, SCUBAS is `LGPL-3.0-or-later`: `LICENSE` (LGPL) and `COPYING` (GPL), an SPDX header on
every Python file, the metadata, and `NOTICE.md`, which keeps the MIT notice for the earlier releases. Those stay available under MIT
(a license already granted is not withdrawn). The change merges only after every holder has agreed.

**Consequences.** The module attribute `__license__` changes in seven modules. Anyone who wants MIT terms uses a release up to 0.1.9.

**Rejected.** Keeping MIT (the owner's choice is the LGPL); a dual license (extra complexity for no stated need).
