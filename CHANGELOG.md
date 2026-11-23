# Changelog

Notable changes to oxiformula, by release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.11.0]

The formula engine leaves oxifit for a repository of its own, as
`oxiformula`: oxifit, oxiroot and oxiroot's dataset expressions all use it,
and it now releases on its own schedule.

### Changed

- **Renamed from `oxifit-formula`.** The API is that of oxifit-formula
  0.10.1, under the crate name `oxiformula`: replace `oxifit_formula::` with
  `oxiformula::`. Its history up to 0.10.1, with each release's notes, is in
  [oxifit's changelog](https://github.com/oxihep/oxifit/blob/v0.10.1/CHANGELOG.md).
