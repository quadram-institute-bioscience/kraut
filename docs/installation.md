# Installation

[Overview](index.md) | [Installation](installation.md) | [Commands](commands.md)


## From Bioconda

[![Conda Version](https://img.shields.io/conda/vn/bioconda/kraut)](https://bioconda.github.io/recipes/kraut/README.html)
[![Conda Downloads](https://img.shields.io/conda/d/bioconda/kraut)](https://bioconda.github.io/recipes/kraut/README.html)

```bash
mamba install -c bioconda kraut
```

## From PyPI

[![PyPI - Version](https://img.shields.io/pypi/v/krautils)](https://pypi.org/project/krautils/)

:warning: Kraut is available on PyPI as `krautils`. 

You can install it using `pip`:

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
