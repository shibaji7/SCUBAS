# API Reference

:::{admonition} Evolving API Surface
:class: warning

SCUBAS API signatures and internals can change between releases as modeling and dataset workflows are refined.
:::

Structured reference for cable modeling, conductivity workflows, data containers, field-conversion models, and supporting utilities.

The SCUBAS API reference is organized by module domain:

1. Cable system and nodal analysis
2. Conductivity profiles and Earth/ocean structure
3. Datasets and profile containers
4. Ocean/electric-field transfer modeling
5. Plotting helpers
6. Shared utilities

:::{admonition} Read this first
:class: note

Start with the module summary page, then use the auto-generated
section for signatures, parameters, and return types.
:::

## API Sections

| Page | What it covers |
|---|---|
| [Cables](cables.md) | end-to-end transmission-line and nodal analysis objects |
| [Ocean/Earth conductivity](conductivity.md), [Datasets](datasets.md) | Earth/ocean profile generation and reusable site containers |
| [Ocean/Earth model](models.md), [Utilities](utils.md) | magnetic-to-electric conversion, FFT helpers, and geometry tools |
| [Plotting](plotlib.md), [Current sources](sources.md) | visualization helpers and current-source field generation |

Use the pages in this section for narrative summaries plus auto-generated API details.
