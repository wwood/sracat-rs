# Changelog

## Unreleased

### Added
- `--progress`: show progress bar for reads processed, with estimated time remaining.

### Fixed
- Performance improvements for reading and writing FIFO streams, which previously were doing too many syscalls.

## Version 0.1.0

### Added
- `--accept-singles`: stream single/orphan reads to stdout interleaved with the
  paired reads as one record-intact stream, instead of requiring a separate
  destination for them.
- `--expect-singles`: invert the default single-handling — stream single/orphan
  reads to stdout and croak if any paired spot is encountered.

## Version 0.0.3

### Added
- Added LICENSE file with MIT license text to repo.

## Version 0.0.2

### Added

### Fixed
- Write more into -1 pipe when specified, so format sniffers e.g. needletail do not block while sracat-rs tries to write to -2 pipe.

## Version 0.0.1

### Changed
- Initial release.