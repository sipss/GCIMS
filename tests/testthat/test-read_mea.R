test_that("read_mea does not crash", {
  mea_file <- system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS")
  sample <- read_mea(mea_file)
  expect_s4_class(sample, "GCIMSSample")
})

test_that("read_mea rejects files without a .mea/.mea.gz extension", {
  expect_error(read_mea("not_a_mea_file.txt"), "\\.mea or \\.mea\\.gz file")
})

test_that("read_mea uses the name of a named filename argument as the description", {
  mea_file <- system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS")

  sample <- read_mea(c(myname = mea_file))

  expect_equal(description(sample), "myname")
})

# Fixture: decompress the bundled small.mea.gz sample, split its header into
# lines, and provide helpers to patch individual "Key = value [unit]" lines
# and re-serialize a synthetic, uncompressed .mea file with the same binary
# spectrum data. Used to exercise read_mea()'s unit-validation branches,
# which the single bundled sample file can't reach on its own.
mea_header_fixture <- function() {
  mea_file <- system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS")
  con <- gzfile(mea_file, "rb")
  raw_bytes <- readBin(con, what = "raw", n = 50 * 1024 * 1024)
  close(con)
  nul_idx <- which(raw_bytes == as.raw(0))[1]
  header_text <- iconv(rawToChar(raw_bytes[1:(nul_idx - 1)]), from = "windows-1252", to = "UTF-8")
  list(
    lines = strsplit(header_text, "\n", fixed = TRUE)[[1]],
    rest = raw_bytes[nul_idx:length(raw_bytes)]
  )
}

replace_mea_header_line <- function(lines, key, new_line) {
  idx <- which(startsWith(trimws(lines), key))
  stopifnot(length(idx) == 1)
  lines[idx] <- new_line
  lines
}

write_mea_from_lines <- function(lines, rest, path = tempfile(fileext = ".mea")) {
  header_raw <- iconv(paste(lines, collapse = "\n"), from = "UTF-8", to = "windows-1252", toRaw = TRUE)[[1]]
  writeBin(c(header_raw, rest), path)
  path
}

test_that("read_mea errors clearly when Chunk sample rate isn't in kHz", {
  fx <- mea_header_fixture()
  lines <- replace_mea_header_line(fx$lines, "Chunk sample rate", "Chunk sample rate              = 150 [Hz]")
  mea_file <- write_mea_from_lines(lines, fx$rest)

  expect_error(read_mea(mea_file), "Expected Chunk sample rate to be in kHz, found Hz instead")
})

test_that("read_mea accepts Chunk trigger repetition in seconds", {
  fx <- mea_header_fixture()
  lines <- replace_mea_header_line(
    fx$lines, "Chunk trigger repetition", "Chunk trigger repetition       = 0.03 [s]"
  )
  mea_file <- write_mea_from_lines(lines, fx$rest)

  sample_s <- read_mea(mea_file)
  sample_ms <- read_mea(system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS"))

  expect_equal(rtime(sample_s), rtime(sample_ms))
})

test_that("read_mea errors clearly when Chunk trigger repetition isn't in ms or s", {
  fx <- mea_header_fixture()
  lines <- replace_mea_header_line(
    fx$lines, "Chunk trigger repetition", "Chunk trigger repetition       = 0.0005 [min]"
  )
  mea_file <- write_mea_from_lines(lines, fx$rest)

  expect_error(read_mea(mea_file), "Expected Chunk trigger repetition to be in ms, found min instead")
})

test_that("read_mea accepts a drift tube length given in mm", {
  fx <- mea_header_fixture()
  lines <- replace_mea_header_line(fx$lines, "nom Drift Tube Length", "nom Drift Tube Length          = 98 [mm]")
  mea_file <- write_mea_from_lines(lines, fx$rest)

  sample <- read_mea(mea_file)

  expect_equal(sample@drift_tube_length, 98)
})

test_that("read_mea errors clearly for an unrecognized drift tube length unit", {
  fx <- mea_header_fixture()
  lines <- replace_mea_header_line(fx$lines, "nom Drift Tube Length", "nom Drift Tube Length          = 98 [cm]")
  mea_file <- write_mea_from_lines(lines, fx$rest)

  expect_error(read_mea(mea_file), "Unit conversion cm->mm not implemented")
})

test_that("read_mea warns once and keeps the raw value for an unrecognized header key", {
  fx <- mea_header_fixture()
  # A single long padding line pushes the header past match_raw()'s 8192-byte
  # scan slice, exercising its multi-slice loop as a side effect.
  padding_value <- paste(rep("X", 9000), collapse = "")
  lines <- c(fx$lines, sprintf('PaddingKey                     = "%s"', padding_value))
  mea_file <- write_mea_from_lines(lines, fx$rest)

  expect_warning(sample <- read_mea(mea_file), "Unknown key: PaddingKey")

  expect_s4_class(sample, "GCIMSSample")
  expect_equal(dim(intensity(sample)), c(1670L, 530L))
})
