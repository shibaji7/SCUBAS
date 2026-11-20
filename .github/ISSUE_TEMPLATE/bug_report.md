---
name: Bug report
about: Something gives a wrong or surprising result
labels: bug
---

**What happened, and what you expected**

**Minimal way to reproduce** (a few lines of Python, the cable sections and the B-field file or its format: IAGA-2002 `.min`/`.txt` or CSV)

**Versions** (output of the line below)

```bash
python -c "import scubas, pandas, numpy, scipy, sys; from importlib.metadata import version; print(version('scubas'), pandas.__version__, numpy.__version__, scipy.__version__, sys.version.split()[0])"
```
