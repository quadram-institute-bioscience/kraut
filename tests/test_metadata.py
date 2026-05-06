import pandas as pd
import pytest

from kraut.models.metadata import Metadata


def test_metadata_parses_csv_with_first_column_as_sample_id(tmp_path):
    metadata_file = tmp_path / "metadata.csv"
    metadata_file.write_text(
        "sample,age,group,collection_date\n"
        "S1,34,control,2024-01-02\n"
        "S2,41,treatment,2024-01-03\n"
    )

    metadata = Metadata.from_file(metadata_file)

    assert metadata.sample_id_column == "sample"
    assert metadata.sample_ids == ["S1", "S2"]
    assert metadata.numerical_columns == ["age"]
    assert metadata.categorical_columns == ["group"]
    assert metadata.date_columns == ["collection_date"]
    assert metadata.column_types == {
        "age": Metadata.NUMERICAL,
        "group": Metadata.CATEGORICAL,
        "collection_date": Metadata.DATE,
    }

    df = metadata.to_dataframe()
    assert pd.api.types.is_numeric_dtype(df["age"])
    assert pd.api.types.is_datetime64_any_dtype(df["collection_date"])


def test_metadata_accepts_named_sample_id_column_in_tsv(tmp_path):
    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text(
        "run\tsample-alias\tdepth\tbatch\n"
        "r001\talpha\t12.5\tA\n"
        "r002\tbeta\t9.0\tB\n"
    )

    metadata = Metadata.from_file(metadata_file, sample_id_column="sample-alias")

    assert metadata.sample_id_column == "sample-alias"
    assert metadata.sample_ids == ["alpha", "beta"]
    assert metadata.numerical_columns == ["depth"]
    assert set(metadata.categorical_columns) == {"run", "batch"}
    assert metadata.get_sample("alpha")["batch"] == "A"


def test_metadata_preserves_missing_values_during_type_inference(tmp_path):
    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text(
        "sample-id\treads\tphenotype\tdate\n"
        "S1\t10\tcase\t2024-01-02\n"
        "S2\t\tcontrol\t\n"
        "S3\t20\t\t2024-01-04\n"
    )

    metadata = Metadata.from_file(metadata_file, sample_id_column="sample-id")
    df = metadata.to_dataframe()

    assert metadata.numerical_columns == ["reads"]
    assert metadata.categorical_columns == ["phenotype"]
    assert metadata.date_columns == ["date"]
    assert pd.isna(df.loc[1, "reads"])
    assert pd.isna(df.loc[1, "date"])
    assert pd.isna(df.loc[2, "phenotype"])


def test_metadata_rejects_missing_sample_id_column(tmp_path):
    metadata_file = tmp_path / "metadata.csv"
    metadata_file.write_text("sample,group\nS1,A\n")

    with pytest.raises(ValueError, match="Sample ID column not found"):
        Metadata.from_file(metadata_file, sample_id_column="sample-id")


def test_metadata_rejects_duplicate_sample_ids(tmp_path):
    metadata_file = tmp_path / "metadata.csv"
    metadata_file.write_text("sample,group\nS1,A\nS1,B\n")

    with pytest.raises(ValueError, match="duplicate"):
        Metadata.from_file(metadata_file)
