# Contributing guide

Scanpy provides extensive [developer documentation][scanpy developer guide], most of which applies to this project, too.
This document will not reproduce the entire content from there.
Instead, it aims at summarizing the most important information to get you started on contributing.

We assume that you are already familiar with git and with making pull requests on GitHub.
If not, please refer to the [scanpy developer guide][].

[scanpy developer guide]: https://scanpy.readthedocs.io/en/latest/dev/index.html

## Installing dev dependencies

In addition to the packages needed to _use_ this package, you need additional python packages to [run tests](#writing-tests) and [build the documentation](#docs-building).

=== "Hatch"

    The easiest way is to get familiar with [hatch environments][], with which these tasks are simply:

    ```bash
    hatch test  # defined in the [tool.hatch.envs.hatch-test] section in pyproject.toml
    hatch run docs:build  # defined in the [tool.hatch.envs.docs] section in pyproject.toml
    ```

=== "Pip"

    If you prefer managing environments manually, you can use `pip`:

    ```bash
    cd hnoca
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -e ".[dev,test,doc]"
    ```

[hatch environments]: https://hatch.pypa.io/latest/tutorials/environment/basic-usage/

## Code-style

This package uses [pre-commit][] to enforce consistent code styles.
On every commit, pre-commit checks will either automatically fix issues with the code or raise an error message.

To enable pre-commit locally, simply run

```bash
pre-commit install
```

in the root of the repository.
Pre-commit will automatically download all dependencies when it is run for the first time.

Alternatively, you can rely on the [pre-commit.ci][] service enabled on GitHub.
If you didn't run `pre-commit` before pushing changes to GitHub it will automatically commit fixes to your pull request or show an error message.

If pre-commit.ci added a commit on a branch you have been working on locally, simply use

```bash
git pull --rebase
```

to integrate the changes into yours.
While the [pre-commit.ci][] service is useful, we strongly encourage installing and running pre-commit locally first to understand its usage.

Finally, most editors have an _autoformat on save_ feature.
Consider enabling this option for [ruff][ruff-editors] and [prettier][prettier-editors].

[pre-commit]: https://pre-commit.com/
[pre-commit.ci]: https://pre-commit.ci/
[ruff-editors]: https://docs.astral.sh/ruff/integrations/
[prettier-editors]: https://prettier.io/docs/en/editors.html

## Writing tests

This package uses [pytest][] for automated testing.
Please write {doc}`scanpy:dev/testing` for every function added to the package.

Most IDEs integrate with pytest and provide a GUI to run tests.
Just point yours to one of the environments returned by

```bash
hatch env create hatch-test  # create test environments for all supported versions
hatch env find hatch-test  # list all possible test environment paths
```

Alternatively, you can run all tests from the command line by executing

=== "Hatch"

    ```bash
    hatch test  # test with the highest supported Python version
    # or
    hatch test --all  # test with all supported Python versions
    ```

=== "Pip"

    ```bash
    source .venv/bin/activate
    pytest
    ```

in the root of the repository.

[pytest]: https://docs.pytest.org/

### Continuous integration

Continuous integration will automatically run the tests on all pull requests and test
against the minimum and maximum supported Python versions.

Additionally, there's a CI job that tests against pre-releases of all dependencies (if there are any).
The purpose of this check is to detect incompatibilities of new package versions early on and
gives you time to fix the issue or reach out to the developers of the dependency before the package is released to a wider audience.

## Publishing a release

### Updating the version number

Version numbers in `pyproject.toml` are inferred automatically from git tags. Please adhere to [Semantic Versioning][semver]. In brief:

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> 1. MAJOR version when you make incompatible API changes,
> 2. MINOR version when you add functionality in a backwards compatible manner, and
> 3. PATCH version when you make backwards compatible bug fixes.
>
> Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

Once you are done, commit and push your changes and navigate to the "Releases" page of this project on GitHub.
Specify `vX.X.X` as a tag name and create a release.
For more information, see [managing GitHub releases][].
This will automatically create a git tag and trigger a GitHub workflow that creates a release on [PyPI][].

[semver]: https://semver.org/
[managing GitHub releases]: https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository
[pypi]: https://pypi.org/

## Writing Documentation

Please write documentation for new or changed features and use-cases. This project now uses [MkDocs](https://www.mkdocs.org/) with the [Material theme](https://squidfunk.github.io/mkdocs-material/) and a set of powerful plugins. Key features include:

- **Markdown-based Authoring:**
  All documentation is written in Markdown, keeping the content easy to read and maintain.

- **Automatic API Reference:**
  The [mkdocstrings](https://mkdocstrings.github.io/) plugin automatically generates API reference documentation from Python docstrings.

- **Jupyter Notebook Integration:**
  With [mknotebooks](https://mknotebooks.github.io/), you can seamlessly include and render Jupyter notebooks as part of your docs. This is ideal for tutorials and live examples.

- **Enhanced Markdown Extensions:**
  The project leverages [pymdown-extensions](https://facelessuser.github.io/pymdown-extensions/) for improved code and content rendering. In particular, we use:
    - `pymdownx.highlight`
    - `pymdownx.superfences`
    - `pymdownx.inlinehilite`
    - `pymdownx.tabbed`

- **Site Configuration:**
  The entire documentation setup—including navigation, theme settings (logo, colors, etc.), and custom CSS—is defined in the `mkdocs.yml` file.

For more insight on writing and maintaining documentation with MkDocs, please refer to the [MkDocs User Guide](https://www.mkdocs.org/user-guide/).

### Tutorials with MkDocs and Jupyter Notebooks

Our documentation integrates Jupyter notebooks using the [mknotebooks](https://mknotebooks.github.io/) plugin. Please ensure any notebooks you include are updated and re-run when necessary, so that both input and output cells remain current.

#### Building the docs locally

=== "Hatch"

    ```bash
    hatch run docs:build
    hatch run docs:open
    ```

=== "Pip"

    ```bash
    source .venv/bin/activate
    # Build the static site
    mkdocs build
    # Serve the site locally
    mkdocs serve
    ```
