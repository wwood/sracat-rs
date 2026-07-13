# Changelog

## Unreleased

## Version 0.2.0

### Added
- `--progress`: show progress bar for reads processed, with estimated time remaining.
- `--eager-open-output`: open the single/orphan output up front instead of lazily
  on the first orphan read. When singles are streamed through a FIFO, the default
  lazy open means the pipe never gets a writer for a run with no orphans, so a
  reader blocking on `open(O_RDONLY)` hangs forever; this flag guarantees the
  reader sees a clean EOF.

### Fixed
- Performance improvements for reading and writing FIFO streams, which previously were doing too many syscalls.
- Intermittent crash ("double free or corruption") under `-t` > 1, caused by
  unserialized concurrent open/close of the ncbi-vdb manager; the open/close
  lifecycle is now serialized while the hot read path stays lock-free.

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