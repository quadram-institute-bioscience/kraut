from pathlib import Path

import pytest


SYNTHETIC_KRAKEN_REPORT = (
    "  5.00\t5\t5\tU\t0\tunclassified\n"
    " 95.00\t95\t0\tR\t1\troot\n"
    " 90.00\t90\t0\tD\t2\t  Bacteria\n"
    " 60.00\t60\t10\tP\t1224\t    Pseudomonadota\n"
    " 40.00\t40\t5\tG\t561\t      Escherichia\n"
    " 30.00\t30\t30\tS\t562\t        Escherichia coli\n"
    " 10.00\t10\t10\tS\t61645\t        Escherichia albertii\n"
    "  5.00\t5\t5\tD\t2157\t  Archaea\n"
)


@pytest.fixture
def synthetic_report_path(tmp_path: Path) -> Path:
    report_path = tmp_path / "synthetic.kreport"
    report_path.write_text(SYNTHETIC_KRAKEN_REPORT)
    return report_path
