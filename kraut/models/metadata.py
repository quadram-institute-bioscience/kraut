import csv
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from pandas.api.types import is_numeric_dtype

DATE_PATTERN = re.compile(
    r"("
    r"\d{4}[-/]\d{1,2}[-/]\d{1,2}"
    r"|\d{1,2}[-/]\d{1,2}[-/]\d{2,4}"
    r"|\d{4}[-/]\d{1,2}"
    r"|\d{1,2}\s+[A-Za-z]{3,9}\s+\d{4}"
    r"|[A-Za-z]{3,9}\s+\d{1,2},?\s+\d{4}"
    r")"
)


class Metadata:
    """Sample metadata parsed from a CSV or TSV file."""

    NUMERICAL = "numerical"
    CATEGORICAL = "categorical"
    DATE = "date"

    def __init__(self, data: pd.DataFrame, sample_id_column: Optional[str] = None):
        if data.empty and len(data.columns) == 0:
            raise ValueError("Metadata must contain a header row")

        self.sample_id_column = self._resolve_sample_id_column(data, sample_id_column)
        self.data = data.copy()
        self.data.columns = [str(column) for column in self.data.columns]
        self._validate_sample_ids()

        self.column_types = self._infer_column_types()
        self.numerical_columns = [
            column
            for column, column_type in self.column_types.items()
            if column_type == self.NUMERICAL
        ]
        self.categorical_columns = [
            column
            for column, column_type in self.column_types.items()
            if column_type == self.CATEGORICAL
        ]
        self.date_columns = [
            column
            for column, column_type in self.column_types.items()
            if column_type == self.DATE
        ]
        self._coerce_inferred_columns()

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        sample_id_column: Optional[str] = None,
        delimiter: Optional[str] = None,
    ) -> "Metadata":
        """Parse sample metadata from a comma- or tab-delimited file."""
        metadata_path = Path(path)
        if delimiter is None:
            delimiter = cls._detect_delimiter(metadata_path)

        data = pd.read_csv(
            metadata_path,
            sep=delimiter,
            dtype=str,
            keep_default_na=False,
        )
        return cls(data, sample_id_column=sample_id_column)

    @property
    def sample_ids(self) -> List[str]:
        return self.data[self.sample_id_column].tolist()

    def to_dataframe(self) -> pd.DataFrame:
        return self.data.copy()

    def get_sample(self, sample_id: str) -> pd.Series:
        matches = self.data[self.data[self.sample_id_column] == str(sample_id)]
        if matches.empty:
            raise KeyError(f"Sample ID not found: {sample_id}")
        return matches.iloc[0].copy()

    def type_for(self, column: str) -> str:
        if column == self.sample_id_column:
            raise ValueError("Sample ID column does not have an inferred metadata type")
        try:
            return self.column_types[column]
        except KeyError as exc:
            raise KeyError(f"Metadata column not found: {column}") from exc

    @staticmethod
    def _detect_delimiter(path: Path) -> str:
        suffix = path.suffix.lower()
        if suffix in {".tsv", ".tab"}:
            return "\t"
        if suffix == ".csv":
            return ","

        sample = path.read_text()[:4096]
        try:
            return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
        except csv.Error:
            first_line = sample.splitlines()[0] if sample.splitlines() else ""
            return "\t" if "\t" in first_line else ","

    @staticmethod
    def _resolve_sample_id_column(
        data: pd.DataFrame,
        sample_id_column: Optional[str],
    ) -> str:
        columns = [str(column) for column in data.columns]
        if not columns:
            raise ValueError("Metadata must contain at least one column")

        if sample_id_column is None:
            return columns[0]

        if sample_id_column not in columns:
            raise ValueError(f"Sample ID column not found: {sample_id_column}")
        return sample_id_column

    def _validate_sample_ids(self) -> None:
        sample_ids = self.data[self.sample_id_column].astype(str).str.strip()
        if sample_ids.eq("").any():
            raise ValueError("Sample ID column contains empty values")
        if sample_ids.duplicated().any():
            duplicates = sorted(sample_ids[sample_ids.duplicated()].unique())
            duplicate_list = ", ".join(duplicates)
            raise ValueError(
                f"Sample ID column contains duplicate values: {duplicate_list}"
            )

        self.data[self.sample_id_column] = sample_ids

    def _infer_column_types(self) -> Dict[str, str]:
        column_types = {}
        for column in self.data.columns:
            if column == self.sample_id_column:
                continue

            series = self.data[column]
            if self._is_numerical(series):
                column_types[column] = self.NUMERICAL
            elif self._is_date(series):
                column_types[column] = self.DATE
            else:
                column_types[column] = self.CATEGORICAL

        return column_types

    def _coerce_inferred_columns(self) -> None:
        for column in self.numerical_columns:
            self.data[column] = pd.to_numeric(self._blank_to_na(self.data[column]))

        for column in self.date_columns:
            self.data[column] = pd.to_datetime(
                self._blank_to_na(self.data[column]),
                errors="coerce",
            )

        for column in self.categorical_columns:
            self.data[column] = self._blank_to_na(self.data[column]).astype("string")

    @classmethod
    def _is_numerical(cls, series: pd.Series) -> bool:
        if is_numeric_dtype(series):
            return True

        values = cls._non_empty_values(series)
        if values.empty:
            return False

        parsed = pd.to_numeric(values, errors="coerce")
        return parsed.notna().all()

    @classmethod
    def _is_date(cls, series: pd.Series) -> bool:
        values = cls._non_empty_values(series)
        if values.empty:
            return False

        text_values = values.astype(str).str.strip()
        if not text_values.map(lambda value: bool(DATE_PATTERN.search(value))).all():
            return False

        parsed = pd.to_datetime(text_values, errors="coerce")
        return parsed.notna().all()

    @staticmethod
    def _blank_to_na(series: pd.Series) -> pd.Series:
        return series.replace(r"^\s*$", pd.NA, regex=True)

    @staticmethod
    def _non_empty_values(series: pd.Series) -> pd.Series:
        values = Metadata._blank_to_na(series.astype(str).str.strip())
        return values.dropna()
