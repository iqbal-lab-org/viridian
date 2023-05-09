from collections import defaultdict, namedtuple
import gzip
import statistics

import intervaltree

from viridian import constants

ACGT = ["A", "a", "C", "c", "G", "g", "T", "t"]
ACGT_X = [f"X_{x}" for x in ACGT]
INDEL_COLS = ["I", "i", "D", "d"]
INDEL_COLS_X = [f"X_{x}" for x in INDEL_COLS]
ACGT_INDEL = ACGT + INDEL_COLS
ACGT_INDEL_X = [f"X_{x}" for x in ACGT_INDEL]
ALL_COUNTS_COLS = ACGT + INDEL_COLS + ACGT_X + INDEL_COLS_X

TsvRow = namedtuple(
    "TsvRow",
    [
        "Ref_pos",
        "Ref_nt",
        "Cons_pos",
        "Cons_nt",
        "Masked_cons_nt",
        "Amplicon",
        "Primer",
        "Mask",
        "Total_depth",
        "Clean_depth",
        "Cons_depth",
    ]
    + ACGT
    + INDEL_COLS
    + ACGT_X
    + INDEL_COLS_X,
)


def unique_amp_coords(amplicons):
    unique_coords = []
    for i, amp in enumerate(amplicons):
        start = amp["start"] if i == 0 else amplicons[i - 1]["end"] + 1
        end = amplicons[i + 1]["start"] - 1 if i + 1 < len(amplicons) else amp["end"]
        unique_coords.append((start, end))
    return unique_coords


def get_dropped_amplicons(amplicons, ref_msa, cons_msa, max_n_percent=50):
    assert len(ref_msa) == len(cons_msa)
    no_ref_gap_cons = [cons_msa[i] for i, c in enumerate(ref_msa) if c != "-"]
    dropped = set()

    for i, (start, end) in enumerate(unique_amp_coords(amplicons)):
        n_count = no_ref_gap_cons[start:end].count("N")
        # There may be no unique part because (eg happens in ampliseq) the
        # amplicon could be covered by the amplicons to its left + right. Hence
        # check that end - star + 1 > 0.
        if end - start + 1 > 0 and 100 * (n_count / (end - start + 1)) >= max_n_percent:
            dropped.add(i)
    return dropped


def base_counts_to_pass_hets(counts, min_pc, upper_or_lower, min_single=5):
    total = 0
    for c in "A", "C", "G", "T":
        total += counts[upper_or_lower(c)]
    if total == 0:
        return -1, set()
    passed = {
        c
        for c in ["A", "C", "G", "T"]
        if 100 * counts[upper_or_lower(c)] / total >= min_pc
        and counts[upper_or_lower(c)] >= min_single
    }

    if len(passed) > 1:
        return total, passed
    elif len(passed) == 1:
        return total, passed.pop().upper()
    else:
        return -1, set()


def base_counts_to_iupac(counts, min_total=5, min_pc=10.0):
    upper_total, upper_pass = base_counts_to_pass_hets(counts, min_pc, str.upper)
    lower_total, lower_pass = base_counts_to_pass_hets(counts, min_pc, str.lower)
    if len(upper_pass) <= 1 and len(lower_pass) <= 1:
        return None
    bases_to_use = None
    if upper_total < min_total and lower_total < min_total:
        return None
    elif upper_total >= min_total:
        if lower_total < min_total:
            bases_to_use = upper_pass
        elif lower_pass == upper_pass:
            bases_to_use = lower_pass
        else:
            return None
    elif lower_total >= min_total:
        bases_to_use = lower_pass
    else:
        return None

    if len(bases_to_use) <= 1:
        return None
    else:
        return constants.IUPAC_REV["".join(sorted(bases_to_use))]


def base_counts_to_insertion(indel_counts, total_depth, min_prop):
    if total_depth == 0:
        return False
    ins_counts = {}
    for i_or_I in ["i", "I"]:
        for ins, count in indel_counts[i_or_I].items():
            ins_counts[ins.upper()] = ins_counts.get(ins.upper(), 0) + count

    return any(x for x in ins_counts.values() if x / total_depth >= min_prop)


class Qc:
    def __init__(
        self,
        amplicons,
        ref_msa,
        cons_msa,
        max_amp_n_percent=50,
        mask_min_depth=10,
        mask_min_frs=0.5,
        het_min_pc=10.0,
    ):
        self.amplicons = amplicons
        self.ref_msa = ref_msa
        self.cons_msa = cons_msa
        self.max_amp_n_percent = max_amp_n_percent
        self.mask_min_depth = mask_min_depth
        self.mask_min_frs = mask_min_frs
        self.het_min_pc = het_min_pc
        self.unique_amp_coords = unique_amp_coords(self.amplicons)
        self.dropped_amplicons = get_dropped_amplicons(
            self.amplicons,
            self.ref_msa,
            self.cons_msa,
            max_n_percent=self.max_amp_n_percent,
        )
        self._init_pileup()
        self._init_amp_and_primer_lookup()
        self.masked_cons = ""
        self.masked_cons_msa = ""
        self.masked_cons_msa_indel_as_ref = ""
        self.masked_cons_msa_indel_as_N = ""

    def _init_amp_and_primer_lookup(self):
        self.amp_tree = intervaltree.IntervalTree()
        self.primer_tree = intervaltree.IntervalTree()
        for amp in self.amplicons:
            self.amp_tree[amp["start"] : amp["end"] + 1] = amp["name"]
            for i, d in enumerate(amp["primers"]["left"]):
                self.primer_tree[d["start"] : d["end"] + 1] = amp["name"] + f"_l_{i}"
            for i, d in enumerate(amp["primers"]["right"]):
                self.primer_tree[d["start"] : d["end"] + 1] = amp["name"] + f"_r_{i}"

    def _init_pileup(self):
        self.pileup = {}
        ref_pos = -1
        cons_pos = -1
        self.ref_to_cons_coord = {}
        for ref_base, cons_base in zip(self.ref_msa, self.cons_msa):
            assert not ref_base == "-" == cons_base
            if ref_base != "-":
                ref_pos += 1
            if cons_base != "-":
                cons_pos += 1
            self.ref_to_cons_coord[ref_pos] = cons_pos
            if cons_pos not in self.pileup:
                counts = {x: 0 for x in ACGT + ACGT_X + INDEL_COLS + INDEL_COLS_X}
                indels = {x: {} for x in INDEL_COLS + INDEL_COLS_X}
                self.pileup[cons_pos] = {
                    "base": cons_base,
                    "ref_data": {},
                    "counts": counts,
                    "indels": indels,
                }
            assert ref_pos not in self.pileup[cons_pos]["ref_data"]
            self.pileup[cons_pos]["ref_data"][ref_pos] = {
                "base": ref_base,
                # "amplicons": set(),
                # "primers": set(),
            }

    def add_one_pileup(
        self, amp_index_to_add, left_primer_index, right_primer_index, pileup_dict
    ):
        amp = self.amplicons[amp_index_to_add]
        left_primer = amp["primers"]["left"][left_primer_index]
        right_primer = amp["primers"]["right"][right_primer_index]
        cons_start = self.ref_to_cons_coord[amp["start"]]
        cons_end = self.ref_to_cons_coord[amp["end"]]

        for cons_pos in range(cons_start, cons_end + 1):
            key_prefix = ""

            for ref_pos, d in self.pileup[cons_pos]["ref_data"].items():
                if (
                    left_primer["start"] <= ref_pos <= left_primer["end"]
                    or right_primer["start"] <= ref_pos <= right_primer["end"]
                ):
                    key_prefix = "X_"
                break

            new_pileup = pileup_dict.get(cons_pos, None)
            if new_pileup is None:
                continue

            for x in ACGT_INDEL:
                self.pileup[cons_pos]["counts"][f"{key_prefix}{x}"] += new_pileup.get(
                    x, 0
                )
            for d_or_i in INDEL_COLS:
                to_update = self.pileup[cons_pos]["indels"][f"{key_prefix}{d_or_i}"]
                for k, v in new_pileup["indel"][d_or_i].items():
                    to_update[k] = to_update.get(k, 0) + v

    def make_pileup(self, pileup_dict):
        for amp_index in pileup_dict:
            for (left_p_index, right_p_index), counts in pileup_dict[amp_index].items():
                self.add_one_pileup(amp_index, left_p_index, right_p_index, counts)

    def pileup_to_tsv_mask_and_count_fields(self, cons_base, pileup):
        mask = []
        base_counts = pileup["counts"]
        iupac = None
        insertion = False

        if cons_base == "-":
            mask = ["."]
            counts = ["."] * (3 + 2 * (len(ACGT) + 4))
        else:
            if cons_base == "N":
                mask.append("ASY")
                cons_depth = 0
            elif cons_base in constants.IUPAC:
                cons_depth = sum(
                    [
                        base_counts[x] + base_counts[x.lower()]
                        for x in constants.IUPAC[cons_base]
                    ]
                )
            else:
                cons_depth = base_counts[cons_base] + base_counts[cons_base.lower()]

            clean_depth = sum([base_counts[x] for x in ACGT_INDEL])
            total_depth = clean_depth + sum([base_counts[x] for x in ACGT_INDEL_X])

            if clean_depth < self.mask_min_depth:
                mask.append("DEPTH")

            if clean_depth > 0:
                if cons_depth / clean_depth < self.mask_min_frs:
                    mask.append("FRS")
                insertion = base_counts_to_insertion(
                    pileup["indels"], clean_depth, constants.MIN_INS_PROPORTION
                )

            iupac = base_counts_to_iupac(base_counts, min_pc=self.het_min_pc)
            if iupac is not None:
                mask.append("HET")

            counts = [
                total_depth,
                clean_depth,
                cons_depth,
                *[base_counts[x] for x in ALL_COUNTS_COLS],
            ]

        if len(mask) == 0:
            mask = ["PASS"]

        return mask, counts, iupac, insertion

    def make_tsv_lines_and_masked_cons(self):
        self.tsv_lines = defaultdict(list)  # ref position -> list of TsvRows
        masked_cons = []
        self.masked_positions = {}
        self.insertions = {}
        self.mask_counts = defaultdict(int)
        not_actually_masked = (["."], ["PASS"])

        for cons_pos in self.pileup:
            cons_base = self.pileup[cons_pos]["base"]
            masked = False
            mask, counts, iupac, insertion = self.pileup_to_tsv_mask_and_count_fields(
                cons_base, self.pileup[cons_pos]
            )
            if insertion:
                to_add = {
                    k: v
                    for k, v in self.pileup[cons_pos]["indels"].items()
                    if len(v) > 0
                }
                to_add.update({k: self.pileup[cons_pos]["counts"][k] for k in ACGT})
                self.insertions[cons_pos] = to_add

            if mask not in not_actually_masked:
                self.mask_counts[";".join(sorted(mask))] += 1
                ref_data = self.pileup[cons_pos]["ref_data"]
                ref_bases = "".join(ref_data[k]["base"] for k in sorted(ref_data))
                self.masked_positions[cons_pos] = {
                    "ref_start": min(ref_data.keys()),
                    "ref_end": max(ref_data.keys()),
                    "ref": ref_bases,
                    "cons": cons_base,
                    "mask": mask,
                    "counts": self.pileup[cons_pos]["counts"],
                    "indels": self.pileup[cons_pos]["indels"],
                }

            first = True

            for ref_pos, d in sorted(self.pileup[cons_pos]["ref_data"].items()):
                amps = [x[2] for x in self.amp_tree[ref_pos]]
                if len(amps) == 0:
                    amps = ["."]

                primers = [x[2] for x in self.primer_tree[ref_pos]]
                if len(primers) == 0:
                    primers = ["."]

                if cons_base != "-" and mask != ["PASS"]:
                    masked = True

                self.tsv_lines[ref_pos].append(
                    TsvRow(
                        ref_pos + 1,
                        d["base"],
                        cons_pos + 1,
                        cons_base,
                        None,  # masked cons base, we don't know it yet
                        ";".join(amps),
                        ";".join(primers),
                        ";".join(sorted(mask)) if first else ".",
                        *counts,
                    )
                )
                cons_base = "-"
                counts = ["."] * len(counts)
                first = False

            if self.pileup[cons_pos]["base"] != "-":
                if (
                    masked
                    and iupac is not None
                    and (
                        mask == ["FRS", "HET"]
                        or set(mask).isdisjoint({"ASY", "DEPTH", "FRS"})
                    )
                ):
                    masked_cons.append(iupac)
                elif masked:
                    masked_cons.append("N")
                else:
                    masked_cons.append(self.pileup[cons_pos]["base"])

        self.masked_cons = "".join(masked_cons)

    def write_qc_tsv_and_make_masked_cons_msa(self, outfile):
        msa = []
        msa_indel_ref = []
        msa_indel_N = []

        with gzip.open(outfile, "wt") as f:
            print(*TsvRow._fields, sep="\t", file=f)
            for ref_pos, lines in sorted(self.tsv_lines.items()):
                for line in lines:
                    if line.Cons_nt == "-":
                        mask_nt = "-"
                        assert line.Ref_nt != "-"
                        msa_indel_ref.append(line.Ref_nt)
                        msa_indel_N.append("N")
                    else:
                        mask_nt = self.masked_cons[line.Cons_pos - 1]
                        if line.Ref_nt != "-":
                            msa_indel_ref.append(mask_nt)
                            msa_indel_N.append(mask_nt)
                    msa.append(mask_nt)

                    print(*line[:4], mask_nt, *line[5:], sep="\t", file=f)

        # Force indels at start/end of indel-as-ref MSA to be N, not -
        assert len(msa_indel_ref) == len(msa_indel_N)
        i = 0
        while i < len(msa_indel_ref) and msa_indel_N[i] == "N":
            msa_indel_ref[i] = "N"
            i += 1
        i = len(msa_indel_ref) - 1
        while i > 0 and msa_indel_N[i] == "N":
            msa_indel_ref[i] = "N"
            i -= 1

        self.masked_cons_msa = "".join(msa)
        self.masked_cons_msa_indel_as_ref = "".join(msa_indel_ref)
        self.masked_cons_msa_indel_as_N = "".join(msa_indel_N)

    def annotated_vcf_file(
        self, vcf_header, vcf_records, outfile, sample_name="sample"
    ):
        assert vcf_header[-1].startswith("#CHROM\t")
        chrom_fields = vcf_header[-1].split("\t")
        chrom_fields[-1] = sample_name
        vcf_header[-1] = "\t".join(chrom_fields)

        new_headers = [
            '##INFO=<ID=AMP,Number=.,Type=String,Description="List of amplicon(s) overlapping the REF allele">',
            '##INFO=<ID=PRIMER,Number=.,Type=String,Description="List of primers(s) overlapping the REF allele">',
            '##INFO=<ID=CONS_POS,Number=1,Type=Integer,Description="Start position of ALT allele in consensus sequence">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total clean read depth">',
            '##FORMAT=<ID=CDP,Number=1,Type=Integer,Description="Total clean read depth that agrees with the ALT allele (ie the consensus sequence)">',
            '##FILTER=<ID=PASS,Description="Variant passed all filters">',
            '##FILTER=<ID=ASY,Description="Assembly (before QC) put an N at this position">',
            '##FILTER=<ID=DEPTH,Description="Total clean read depth (DP) is too low">',
            '##FILTER=<ID=FRS,Description="Fraction of Read Support (CDP / DP) is too low">',
        ]

        with open(outfile, "w") as f:
            print(*vcf_header[:-1], *new_headers, vcf_header[-1], sep="\n", file=f)

            for record in vcf_records:
                tsv_lines = self.tsv_lines[record.POS]
                amps = set()
                primers = set()
                record.FILTER = set()
                clean_depths = []
                cons_depths = []

                for line in tsv_lines:
                    amps.update(set(line.Amplicon.split(";")))
                    primers.update(set(line.Primer.split(";")))
                    record.FILTER.update(set(line.Mask.split(";")))
                    clean_depths.append(line.Clean_depth)
                    cons_depths.append(line.Cons_depth)

                assert len(amps) > 0
                record.INFO["AMP"] = ",".join(sorted(list(amps)))
                if primers != {"."}:
                    record.INFO["PRIMER"] = ",".join(sorted(list(primers)))
                record.INFO["CONS_POS"] = str(tsv_lines[0].Cons_pos)

                if "PASS" in record.FILTER and len(record.FILTER) > 1:
                    record.FILTER.remove("PASS")

                record.set_format_key_value(
                    "DP", str(round(statistics.mean(clean_depths)))
                )
                record.set_format_key_value(
                    "CDP", str(round(statistics.mean(cons_depths)))
                )
                print(record, file=f)
