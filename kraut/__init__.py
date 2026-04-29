from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    tomllib = None


def _version_from_pyproject() -> str | None:
    if tomllib is None:
        return None

    pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
    try:
        with pyproject.open("rb") as handle:
            return tomllib.load(handle).get("project", {}).get("version")
    except OSError:
        return None


def _version_from_metadata() -> str:
    try:
        return version("krautils")
    except PackageNotFoundError:
        return "0+unknown"


__version__ = _version_from_pyproject() or _version_from_metadata()
