# Installation

[Overview](index.md) | [Installation](installation.md) | [Commands](commands.md)

Kraut is available on PyPI as `krautils`. You can install it using `pip` or from the source repository.

## From PyPI

To install the latest stable version:

```bash
pip install krautils
```

## From Source

To install the development version from a local clone of the repository:

```bash
git clone https://github.com/quadram-institute-bioscience/kraut.git
cd kraut
pip install .
```

## Developer Installation

If you wish to contribute to Kraut or run the test suite, install the package with developer dependencies:

```bash
pip install -e ".[dev]"
```

### Running Tests

You can verify the installation by running the test suite:

```bash
pytest
```
