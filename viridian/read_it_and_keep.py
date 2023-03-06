import os

from viridian import utils


def run_riak(outprefix, ref_fasta, reads1, reads2=None):
    if reads2 is None:
        reads_opts = f"--tech ont --reads1 {reads1}"
    else:
        reads_opts = f"--tech illumina --reads1 {reads1} --reads2 {reads2}"

    try:
        result = utils.syscall(
            f"readItAndKeep {reads_opts} --ref_fasta {ref_fasta} -o {outprefix}"
        )
    except:
        return {}, {}, "Error running readItAndKeep"

    expect_keys = {
        "Input reads file 1",
        "Input reads file 2",
        "Kept reads 1",
        "Kept reads 2",
    }

    try:
        counts = {}
        for line in result.stdout.split("\n"):
            fields = line.split("\t")
            counts[fields[0]] = int(fields[1])
        assert set(counts.keys()) == expect_keys
    except:
        return {}, {}, "Error parsing readItAndKeep results"

    outfiles = {"reads": None, "reads_1": None, "reads_2": None}
    for key in outfiles:
        fa = f"{outprefix}.{key}.fasta.gz"
        fq = f"{outprefix}.{key}.fastq.gz"
        if os.path.exists(fa):
            assert not os.path.exists(fq)
            outfiles[key] = fa
        if os.path.exists(fq):
            assert not os.path.exists(fa)
            outfiles[key] = fq

    if reads2 is None:
        if outfiles["reads"] is None:
            return counts, outfiles, "Reads output file not found"
        outfiles["reads_1"] = outfiles["reads"]
        del outfiles["reads"]
    else:
        del outfiles["reads"]
        if outfiles["reads_1"] is None:
            return counts, outfiles, "Reads output file 1 not found"
        if outfiles["reads_2"] is None:
            return counts, outfiles, "Reads output file 2 not found"

    return counts, outfiles, None
