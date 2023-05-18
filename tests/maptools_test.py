from collections import defaultdict
import os
import random

import pyfastaq
import pytest

from viridian import maptools, utils


def test_map_reads():
    random.seed(42)
    ref_fa = "tmp.map_reads.ref.fa"
    ref_seq = pyfastaq.sequences.Fasta(
        "ref", "".join(random.choices(["A", "C", "G", "T"], k=1000))
    )
    with open(ref_fa, "w") as f:
        print(ref_seq, file=f)

    reads_1_fa = "tmp.map_reads.reads.1.fa"
    reads_2_fa = "tmp.map_reads.reads.2.fa"
    with open(reads_1_fa, "w") as f:
        print(">r1 /1", ref_seq[50:200], sep="\n", file=f)
    with open(reads_2_fa, "w") as f:
        seq = pyfastaq.sequences.Fasta("r1 /2", ref_seq[250:400])
        seq.revcomp()
        print(seq, file=f)

    bam = "tmp.map_reads.bam"
    bai = f"{bam}.bai"
    utils.syscall(f"rm -f {bam} {bai}")
    maptools.map_reads(bam, ref_fa, reads_1_fa)
    assert os.path.exists(bam)
    assert os.path.exists(bai)
    utils.syscall(f"rm {bam} {bai}")
    maptools.map_reads(bam, ref_fa, reads_1_fa, reads2=reads_2_fa)
    assert os.path.exists(bam)
    assert os.path.exists(bai)
    utils.syscall(f"rm {bam} {bai} {reads_1_fa} {reads_2_fa} {ref_fa}")


def test_pileup_indel_parse():
    assert maptools.pileup_indel_parse("A") is None
    assert maptools.pileup_indel_parse("A+2") is None
    assert maptools.pileup_indel_parse("+2AT") is None
    assert maptools.pileup_indel_parse("A+2GT") == {
        "ref": "A",
        "sign": "+",
        "length": "2",
        "alt": "GT",
    }
    assert maptools.pileup_indel_parse("C-3GAT") == {
        "ref": "C",
        "sign": "-",
        "length": "3",
        "alt": "GAT",
    }


def test_indels_from_pileup_counts():
    counts = {"A": 10, "a": 5, "C": 42}
    expect = {x: defaultdict(int) for x in ["D", "d", "I", "i"]}
    assert maptools.indels_from_pileup_counts(counts) == expect
    counts["A+1T"] = 2
    expect["I"]["T"] = 2
    assert maptools.indels_from_pileup_counts(counts) == expect
    counts["A+2GC"] = 3
    expect["I"]["GC"] = 3
    assert maptools.indels_from_pileup_counts(counts) == expect
    counts["a+2gc"] = 4
    expect["i"]["gc"] = 4
    assert maptools.indels_from_pileup_counts(counts) == expect
    counts["a-4acgt"] = 8
    expect["d"][4] = 8
    assert maptools.indels_from_pileup_counts(counts) == expect
    counts["A-3ACG"] = 16
    expect["D"][3] = 16
    assert maptools.indels_from_pileup_counts(counts) == expect


def test_pileup():
    random.seed(42)
    ref_fa = "tmp.pileup.ref.fa"
    bam = "tmp.pileup.bam"
    reads_1_fa = "tmp.pileup.reads.1.fa"
    reads_2_fa = "tmp.pileup.reads.2.fa"
    utils.syscall(f"rm -f {ref_fa}* {bam}* {reads_1_fa} {reads_2_fa}")
    ref_bases = random.choices(["A", "C", "G", "T"], k=200)
    ref_bases[20] = "A"
    ref_bases[21] = "G"
    ref_bases[22] = "T"
    ref_bases[79] = "G"
    ref_bases[80] = "A"
    ref_bases[81] = "A"
    ref_bases[130] = "A"
    ref_bases[150] = "A"
    ref_bases[160] = "T"
    ref_bases[161] = "G"
    ref_bases[162] = "C"
    read1 = ref_bases[10:100]
    read2 = ref_bases[10:21] + ref_bases[23:100]
    read3 = ref_bases[10:81] + ["GTC"] + ref_bases[81:100]
    read1_rev = pyfastaq.sequences.Fasta("r1 /2", "".join(ref_bases[90:190]))
    read1_rev.revcomp()
    read2_rev = pyfastaq.sequences.Fasta(
        "r2 /2", "".join(ref_bases[90:161] + ref_bases[162:190])
    )
    read2_rev.revcomp()
    read3_rev = pyfastaq.sequences.Fasta(
        "r3 /2", "".join(ref_bases[90:161] + ["C"] + ref_bases[162:190])
    )
    read3_rev.revcomp()

    with open(reads_1_fa, "w") as f:
        print(">r1 /1", "".join(read1), sep="\n", file=f)
        print(">r1.1 /1", "".join(read1), sep="\n", file=f)
        print(">r2 /1", "".join(read2), sep="\n", file=f)
        print(">r2.2 /1", "".join(read2), sep="\n", file=f)
        print(">r3 /1", "".join(read3), sep="\n", file=f)
    with open(reads_2_fa, "w") as f:
        print(read1_rev, file=f)
        read1_rev.id = "r1.1 /2"
        print(read1_rev, file=f)
        print(read2_rev, file=f)
        read2_rev.id = "r2.1 /2"
        print(read2_rev, file=f)
        print(read3_rev, file=f)

    ref_bases[150] = "G"
    ref_seq = pyfastaq.sequences.Fasta("ref", "".join(ref_bases))
    with open(ref_fa, "w") as f:
        print(ref_seq, file=f)

    bam = "tmp.pileup.bam"
    got_pileup = maptools.pileup(bam, ref_fa, reads_1_fa, reads2=reads_2_fa)
    expect = {i: {ref_bases[i]: 5} for i in range(10, 90)}
    expect.update(
        {i: {ref_bases[i]: 5, ref_bases[i].lower(): 5} for i in range(90, 100)}
    )
    expect.update({i: {ref_bases[i].lower(): 5} for i in range(100, 190)})
    for v in expect.values():
        v["I"] = v["i"] = v["D"] = v["d"] = 0
        v["indel"] = {x: {} for x in ["D", "d", "I", "i"]}

    expect[20].update({"A": 3, "A-2GT": 2})
    expect[21].update({"G": 3, "D": 2})
    expect[22].update({"T": 3, "D": 2})
    expect[80].update({"A": 4, "A+3GTC": 1, "I": 1})
    expect[150]["a"] = 5
    del expect[150]["g"]
    expect[160].update({"t": 3, "t-1g": 2})
    expect[161].update({"g": 2, "c": 1, "d": 2})
    expect[21]["indel"]["D"][2] = 2
    expect[80]["indel"]["I"]["GTC"] = 1
    expect[161]["indel"]["d"][1] = 2
    assert got_pileup == expect
    utils.syscall(f"rm -f {ref_fa}* {bam}* {reads_1_fa} {reads_2_fa}")
