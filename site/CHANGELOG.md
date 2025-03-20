# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## [Unreleased]

### Added

- Added a contribution guide to document how this package is set up.
- Added a `build.yaml` workflow to test whether the package can be build sucessfully.
- Added a `test.yaml` workflow to automate tests with github actions.
- Added pre-commit hooks via the `pyproject.toml` and `.pre-commit-config.yaml` files.
- Added the `pymdownx.tabbed` extension to `mkdocs` to display tabs in the docs.
- Added this Changelog.
- Added a `codecov.yaml` file, to configure coverage reports.
- Added a utility class in `utils/check` to check for optional dependencies.
- Added a basic logger.

### Changed
- The package now uses a `pyproject.toml` file to define the build system (hatch), dependencies, etc.
- Update the `publish.yaml` workflow in `release.yaml`.
- Moved code from `hnoca` to `src/hnoca`.
- Updated dependencies in `pyproject.toml`. Extra dependencies for `maping` and `stats` are optional, to ease installation.

### Removed
- Removed the `black` badge from the `README.md`, as we did not automatically check for adherence to black code style.
- Removed the rendered docs in `docs` from git tracking, as it's not required and it yields unnecessary large git diffs when building the documentation locally.
