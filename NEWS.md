# GCIMS (development version)

- `getChromatogram()` and `getSpectrum()` are now S4 generics with methods
  for `GCIMSDataset`, in addition to the existing `GCIMSSample` method.
  Calling them on a dataset returns a new `GCIMSChromatogramSet`/
  `GCIMSSpectrumSet` object: one chromatogram/spectrum per sample (each kept
  on its own native axis, no interpolation across samples), together with a
  copy of `pData()`. Both new classes have a `plot()` method that can color
  by `SampleID` or by any `pData()` column via `color_by`.
- `GCIMSSample()` now validates the object right after construction, so
  passing an intensity matrix in the wrong orientation (drift time/retention
  time swapped) errors immediately with a clear message instead of
  constructing an invalid object that fails confusingly later on (#28)
- Fix the "Importing custom data formats" vignette, which built the
  intensity matrix transposed relative to what `GCIMSSample()` expects (#28)
- Remove the `Remotes: sipss/pow` field and drop `pow` from `Suggests`,
  since the `pow` repository is currently private and can't be installed
  by CI or by users from CRAN/Bioconductor. `align(method_rt = "pow")` now
  fails with a clear, actionable error (instead of breaking package
  installation for everyone) when `pow` isn't installed. Restore both
  DESCRIPTION fields if/when `pow` becomes publicly installable.

# GCIMS 0.1.1 (2024-04-29)

- Fix undefined variable on RIP saturation detection (#21)

# GCIMS 0.1.0

- Generation of the package
- Generic functions for reading, alignment, and visualization of GC-IMS samples
