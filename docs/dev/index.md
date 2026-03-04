# API Reference

!!! warning "Evolving API Surface"
    SCUBAS API signatures and internals can change between releases as modeling and dataset workflows are refined.

<div class="hero">
  <h2>SCUBAS Developer API Guide</h2>
  <p>
    Structured reference for cable modeling, conductivity workflows, data
    containers, field-conversion models, and supporting utilities.
  </p>
</div>

The SCUBAS API reference is organized by module domain:

1. Cable system and nodal analysis
2. Conductivity profiles and Earth/ocean structure
3. Datasets and profile containers
4. Ocean/electric-field transfer modeling
5. Plotting helpers
6. Shared utilities

<span class="api-badge api-package">Package</span>
<span class="api-badge api-module">Module</span>
<span class="api-badge api-class">Class</span>
<span class="api-badge api-method">Method / Function</span>

!!! info "Read this first"
    Start with the module summary page, then use the auto-generated `mkdocstrings`
    section for signatures, parameters, and return types.

## API Sections

<div class="doc-card-grid">
  <div class="doc-card">
    <strong>Cables</strong>
    End-to-end transmission-line and nodal analysis objects.<br>
    <a href="cables/">Open Cables API</a>
  </div>
  <div class="doc-card">
    <strong>Conductivity + Datasets</strong>
    Earth/ocean profile generation and reusable site containers.<br>
    <a href="conductivity/">Conductivity</a> | <a href="datasets/">Datasets</a>
  </div>
  <div class="doc-card">
    <strong>Models + Utilities</strong>
    Magnetic-to-electric conversion, FFT helpers, and geometry tools.<br>
    <a href="models/">Models</a> | <a href="utils/">Utilities</a>
  </div>
  <div class="doc-card">
    <strong>Plotting + Sources</strong>
    Visualization helpers and current-source field generation.<br>
    <a href="plotlib/">Plotting</a> | <a href="sources/">Current Sources</a>
  </div>
</div>

Use the pages in this section for narrative summaries plus auto-generated API details.
