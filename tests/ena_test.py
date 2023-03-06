import pytest

from viridian import ena


def test_fastq_ftp_string_to_files():
    with pytest.raises(Exception):
        ena.fastq_ftp_string_to_files(None)

    with pytest.raises(Exception):
        ena.fastq_ftp_string_to_files("")

    assert ena.fastq_ftp_string_to_files("foo.fq.gz") == {
        "unpaired": ["foo.fq.gz"],
        "paired": [],
    }

    assert ena.fastq_ftp_string_to_files("foo.fq.gz;bar.fq.gz") == {
        "unpaired": ["foo.fq.gz", "bar.fq.gz"],
        "paired": [],
    }

    assert ena.fastq_ftp_string_to_files("a_1.fastq.gz;a_2.fastq.gz") == {
        "unpaired": [],
        "paired": ["a_1.fastq.gz", "a_2.fastq.gz"],
    }
    assert ena.fastq_ftp_string_to_files(
        "foo.fq.gz;bar.fq.gz;a_1.fastq.gz;a_2.fastq.gz"
    ) == {
        "unpaired": ["foo.fq.gz", "bar.fq.gz"],
        "paired": ["a_1.fastq.gz", "a_2.fastq.gz"],
    }

    assert ena.fastq_ftp_string_to_files("a_1.fastq.gz") == {
        "unpaired": ["a_1.fastq.gz"],
        "paired": [],
    }

    assert ena.fastq_ftp_string_to_files("a_1.fastq.gz;a_2.fastq.gz;b_1.fastq.gz") == {
        "unpaired": ["b_1.fastq.gz"],
        "paired": ["a_1.fastq.gz", "a_2.fastq.gz"],
    }
