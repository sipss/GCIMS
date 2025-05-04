# GCIMS 0.2.0 (2025-05-03)

- Deprecate `peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE)`. Please
  remove that argument, since the new default is better. For further details on
  the new default please check the vignette or `?GCIMS::defaultPeakDetectionCWTParams`.

# GCIMS 0.1.1 (2024-04-29)

- Fix undefined variable on RIP saturation detection (#21)

# GCIMS 0.1.0

- Generation of the package
- Generic functions for reading, alignment, and visualization of GC-IMS samples
