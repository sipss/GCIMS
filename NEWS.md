# GCIMS (development version)

- `GCIMSSample()` now validates the object right after construction, so
  passing an intensity matrix in the wrong orientation (drift time/retention
  time swapped) errors immediately with a clear message instead of
  constructing an invalid object that fails confusingly later on (#28)
- Fix the "Importing custom data formats" vignette, which built the
  intensity matrix transposed relative to what `GCIMSSample()` expects (#28)

# GCIMS 0.1.1 (2024-04-29)

- Fix undefined variable on RIP saturation detection (#21)

# GCIMS 0.1.0

- Generation of the package
- Generic functions for reading, alignment, and visualization of GC-IMS samples
