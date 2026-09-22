# PyFiberModes performance benchmarks

Install `.[benchmarking]` and run:

```bash
pytest benchmarks --benchmark-only
```

The suite measures stable public operations rather than private implementation
details, allowing results to be compared between commits with
`pytest-benchmark compare`.
