# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).



## [2.5.0] - 2024-09-20

### Added
- [GCGI-1425](https://jira.oicr.on.ca/browse/GCGI-1425) - Use MANE Select transcript for VEP annotation.

## [2.4.0] - 2024-06-25

### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add Vidarr labels to outputs (metadata changes only).

## [2.3.3] - 2024-02-12

### Added
- [GRD-730](https://jira.oicr.on.ca/browse/GRD-730) - Add `minMemory` parameter to RAM-scaling tasks.

## [2.3.2] - 2023-12-21

### Added
- [GRD-730](https://jira.oicr.on.ca/browse/GRD-730) - Add RAM scaling by chromosome size (performance update).

## [2.3.1] - 2023-09-20

### Added
- [GBS-4248](https://jira.oicr.on.ca/browse/GBS-4248) - Update mouse genome build and reference-specific modules to mm39 (from mm10).

## [2.3.0] - 2023-07-05

### Changed
- Move assembly-specific codes in Olive to WDL.

## [2.2.1] - 2023-05-25

### Added
- Add tumor and normal names as parameters.

## [2.2.0] - 2023-03-01

### Changed
- Move version 2.1.8 to version 2.2.0.

## [2.1.8] - 2023-02-14

### Added
- Add an option for VEP stats.
- Update the VCF2MAF module to 1.6.21b.

## [2.1.7] - 2022-09-22

### Changed
- Migrate 2.1.6 to Vidarr.
- Fix `getSampleNames`.

## [2.1.6] - 2022-04-13

### Changed
- Update to support VEP 105.

## [2.1.51] - 2023-03-01

### Fixed
- Small fix to version 2.1.5, fixing `getSampleNames`.

## [2.1.5] - 2022-06-28

### Changed
- Change `targetBed` file type to `String`.

## [2.1.4] - 2021-06-01

### Changed
- Migrate to Vidarr.

## [2.1.3] - 2021-02-08

### Added
- Addition of `retainInfoProvided` and `updateTagValue` options.

## [2.1.2] - 2020-11-10

### Added
- Support for non-human species.

## [2.1.1] - 2020-09-23

### Fixed
- Fix chromosome order.

## [2.1] - 2020-09-11

### Added
- Run VEP and VCF2MAF in parallel.

## [2.0.2] - 2020-07-02

### Added
- Support normal samples without `groupID`.

## [2.0.1] - 2020-06-08

### Fixed
- Fix the header for tumor-only VCF.

## [2.0] - 2020-05-07

### Added
- WDL workflow for VEP.
