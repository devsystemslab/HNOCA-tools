# Changelog

Starting from version v0.2.0, all notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## [Unreleased]

## v.0.2.1 (2025-04-02)

### Fixed
- Release workflow to PyPI [#12]({{ pr_url }}12)

## v.0.2.0 (2025-04-01)

### Added

- Added a contribution guide to document how this package is set up [#8]({{ pr_url }}8).
- Added a `build.yaml` workflow to test whether the package can be build sucessfully [#8]({{ pr_url }}8).
- Added a `test.yaml` workflow to automate tests with github actions [#8]({{ pr_url }}8).
- Added pre-commit hooks via the `pyproject.toml` and `.pre-commit-config.yaml` files [#8]({{ pr_url }}8).
- Added the `pymdownx.tabbed` extension to `mkdocs` to display tabs in the docs [#8]({{ pr_url }}8).
- Added this Changelog [#8]({{ pr_url }}8).
- Added a `codecov.yaml` file, to configure coverage reports [#8]({{ pr_url }}8).
- Added a utility class in `utils/check` to check for optional dependencies [#8]({{ pr_url }}8).
- Added a basic logger [#8]({{ pr_url }}8).
- Added a check for snapseed to make sure that the marker dict has the right depth for annotate vs annotate_hierarchy [#8]({{ pr_url }}8).
- Added tests for snapseed [#8]({{ pr_url }}8).
- Set up pre-commit.ci

### Fixed
- Fixed snapseed tests by reading leiden clustering from file [#8]({{ pr_url }}8).

### Changed
- The package now uses a `pyproject.toml` file to define the build system (hatch), dependencies, etc [#8]({{ pr_url }}8).
- Update the `publish.yaml` workflow in `release.yaml` [#8]({{ pr_url }}8).
- Moved code from `hnoca` to `src/hnoca` [#8]({{ pr_url }}8).
- Updated dependencies in `pyproject.toml`. Extra dependencies for `maping` and `stats` are optional, to ease installation [#8]({{ pr_url }}8).
- Use vcs-based versioning with [hatch-vcs](https://pypi.org/project/hatch-vcs/) [#10]({{ pr_url }}10).
- Update the `build-site.yaml` workflow to run on PRs, and add package installation in the workflow using the `doc` option [#9]({{ pr_url }}9).

### Removed
- Removed the `black` badge from the `README.md`, as we did not automatically check for adherence to black code style [#8]({{ pr_url }}8).
- Removed the rendered docs in `docs` from git tracking, as it's not required and it yields unnecessary large git diffs when building the documentation locally [#8]({{ pr_url }}8).
