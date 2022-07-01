import sys

from collections import namedtuple, defaultdict

BaseProfile = namedtuple(
    "BaseProfile", ["base", "in_primer", "forward_strand", "amplicon_name",],
)

Config = namedtuple("Config", ["min_frs", "min_depth"],)

default_config = Config(min_frs=0.7, min_depth=10,)


class Pileup:
    """A pileup is an array of Stats objects indexed by position in a reference
    """

    def __init__(self, refseq, msa=None, config=default_config):
        self.config = config
        self.ref = refseq
        self.seq = []

        # 1-based index translation tables
        self.ref_to_consensus = {}
        self.consensus_to_ref = {}
        if msa:
            with open(msa) as msa_fd:
                seq1 = msa_fd.readline().strip()  # consensus
                seq2 = msa_fd.readline().strip()  # reference
                # this is valid for testing when the consensus is complete
                # assert seq1.replace("-", "") == self.ref
                ref = None  # we're using 1-based coords (vcf)
                con = None

                for a, b in zip(seq1, seq2):
                    if a != "-":
                        if ref is None:
                            ref = 0
                        ref += 1
                    if b != "-":
                        if con is None:
                            con = 0
                        con += 1

                    if ref is not None:
                        self.ref_to_consensus[ref] = con
                    if con is not None:
                        self.consensus_to_ref[con] = ref
        else:
            for i in range(0, len(self.ref) + 1):
                self.ref_to_consensus[i] = i
                self.consensus_to_ref[i] = i

        for i, r in enumerate(refseq):
            if i + i in self.consensus_to_ref:
                self.seq.append(Stats(self.consensus_to_ref[i + 1], ref_base=r))
            else:
                self.seq.append(Stats(None, ref_base=r))

        # define the filters
        # filter closures take a Stats and return True on failure
        def test_amplicon_bias(s):
            passing = []
            # if len(s.amplicon_totals) > 2:
            #    for _aa in s.amplicon_totals:
            #        print(
            #            "\t\t-->",
            #            _aa.name,
            #            _aa.start,
            #            _aa.end,
            #            s.amplicon_totals[_aa],
            #            file=sys.stderr,
            #        )

            for amplicon, total in s.amplicon_totals.items():
                if total < self.config.min_depth:
                    continue
                if s.amplicon_totals[amplicon] == 0:
                    continue
                if (
                    s.refs_in_amplicons[amplicon] / s.amplicon_totals[amplicon]
                    < self.config.min_frs
                ):
                    passing.append(False)
                else:
                    passing.append(True)
            if len(passing) > 1 and not all(passing[0] == e for e in passing):
                #                print(f"{s.reference_pos}\t{passing}")
                #                for amplicon, total in s.amplicon_totals.items():
                #                    print(f"\t{amplicon.name}\t{total}\t{s.refs_in_amplicons[amplicon]} / {s.amplicon_totals[amplicon]}")

                # this is a failure
                return True
            return False

        def test_strand_bias(s):
            passing = []
            for amplicon, total in s.amplicon_totals.items():
                if total < self.config.min_depth:
                    continue
                if s.amplicon_totals[amplicon] == 0:
                    continue
                if (
                    s.refs_in_forward_strands[amplicon] / s.amplicon_totals[amplicon]
                    < self.config.min_frs
                ):
                    passing.append(False)
                else:
                    passing.append(True)
            if len(passing) > 1 and not all(passing[0] == e for e in passing):
                # True = failure
                return True
            return False

        self.filters = {
            "low_depth": (
                lambda s: s.total < self.config.min_depth,
                lambda s: f"Insufficient depth; {s.total} < {self.config.min_depth}. {s.total_reads} including primer regions.",
            ),
            "low_frs": (
                lambda s: s.refs / s.total < self.config.min_frs
                if s.total > 0
                else False,
                lambda s: f"Insufficient support of consensus base; {s.refs} / {s.total} < {self.config.min_frs}. {s.total_reads} including primer regions.",
            ),
            #            "amplicon_bias": (
            #                test_amplicon_bias,
            #                lambda s: "Per-amplicon FRS failure",
            #            ),
            #            "strand_bias": (
            #                test_strand_bias,
            #                lambda s: f"Strand-biased FRS failure"),
        }

        # initialise summary for each filter
        self.summary = defaultdict(int)
        for f in self.filters:
            self.summary[f] = 0

    def __getitem__(self, pos):
        if pos > len(self.seq):
            # be mindful of 0 vs. 1 indexing here
            print(f"position too big: {pos} {len(self.seq)}", file=sys.stderr)
            return None
        return self.seq[pos]

    def __setitem__(self, pos, profile):
        raise Exception("don't set Pileup positions")

    def update(self, pos, profile):
        self.seq[pos].update(profile)

    def __len__(self):
        return len(self.seq)

    def mask(self):
        sequence = list(self.ref)
        self.qc = {}
        for position, stats in enumerate(self.seq):
            if position >= len(sequence):
                raise Exception(
                    f"Invalid condition: mapped position {position} greater than consensus length {len(sequence)}"
                )

            self.summary["consensus_length"] += 1
            # print(position, self.consensus_to_ref[position + 1], file=sys.stderr)
            if sequence[position] == "N":
                # if a position is already masked by an upstream process skip it
                self.summary["already_masked"] += 1
                self.summary["total_masked"] += 1
                continue
            elif stats.check_for_failure(self.filters):
                for failure in stats.get_failures():
                    self.summary[failure] += 1
                self.summary["total_masked"] += 1
                sequence[position] = "N"
                self.qc[position] = stats.log
            self.qc["masking_summary"] = self.summary
        return "".join(sequence)

    def dump_tsv(self, tsv):
        fd = open(tsv, "w")
        for pos, stats in enumerate(self.seq):
            cons_pos = self.consensus_to_ref[pos + 1]
            print(f"{cons_pos}\t{pos+1}\t{stats.tsv_row()}", file=fd)
        return tsv

    def annotate_vcf(self, vcf):
        header = []
        records = []

        # TODO: assert 'chromosome' names are the same

        for line in open(vcf):
            line = line.strip()
            if line[0] == "#":
                header.append(line)
                continue
            (
                chrom,
                pos,
                mut_id,
                ref,
                alt,
                qual,
                original_filters,
                info,
                fmt,
                *r,
            ) = line.split("\t")
            pos = int(pos)

            cons_coord = self.ref_to_consensus[pos]
            stats = self.seq[cons_coord - 1]
            info_field = stats.info()
            vcf_filters = original_filters
            if stats.position_failed is None:
                print(f"Warning: attemped to evaluate N basecall", file=sys.stderr)
                continue
            if stats.position_failed:
                if original_filters == "PASS":
                    vcf_filters = ";".join(stats.get_failures())
                else:
                    vcf_filters = ";".join([original_filters, *vcf_filters])

            records.append(
                (chrom, pos, mut_id, ref, alt, qual, vcf_filters, info_field, fmt, *r)
            )

        return header, records


class Stats:
    def __init__(
        self, reference_pos, ref_base=None, cons_base=None, config=default_config,
    ):
        self.alts_in_primer = 0
        self.refs_in_primer = 0
        self.alt_bases = defaultdict(int)

        self.alts_in_amplicons = defaultdict(int)
        self.refs_in_amplicons = defaultdict(int)
        self.refs_in_forward_strands = defaultdict(int)
        self.alts_in_forward_strands = defaultdict(int)
        self.amplicon_totals = defaultdict(int)

        self.alts_forward = 0
        self.refs_forward = 0

        self.alts = 0
        self.refs = 0
        self.total = 0
        self.log = []
        self.total_reads = 0
        self.failures = None

        self.alts_matching_refs = 0
        self.reference_pos = reference_pos
        self.ref_base = ref_base
        self.cons_base = cons_base

        # self.config = config

        self.position_failed = None

    def update(self, profile, alt=None):

        if profile.in_primer:
            if profile.base != self.ref_base:
                self.alts_in_primer += 1
            else:
                self.refs_in_primer += 1
        else:
            self.alt_bases[profile.base] += 1
            if profile.amplicon_name:
                self.amplicon_totals[profile.amplicon_name] += 1

            if profile.base != self.ref_base:
                self.alts += 1
                if profile.amplicon_name:
                    self.alts_in_amplicons[profile.amplicon_name] += 1
                    if profile.forward_strand:
                        self.alts_in_forward_strands[profile.amplicon_name] += 1

                if profile.forward_strand:
                    self.alts_forward += 1

            else:
                self.refs += 1
                if profile.amplicon_name:
                    self.refs_in_amplicons[profile.amplicon_name] += 1
                    if profile.forward_strand:
                        self.refs_in_forward_strands[profile.amplicon_name] += 1

                if profile.forward_strand:
                    self.refs_forward += 1

            self.total += 1

    def check_for_failure(self, filters):
        """return whether a position should be masked

        optionally accepts a mutable reference to a dictionary that summarises masking decisions
        """

        self.position_failed = False

        for filter_name, (filter_func, msg_format) in filters.items():
            if filter_func(self):
                self.log.append(msg_format(self))
                if self.failures is None:
                    self.failures = []
                self.failures.append(filter_name)
                self.position_failed = True

        return self.position_failed

    def get_failures(self):
        """return list of failed filters for VCF FILTER field
        """
        if self.position_failed is None:
            raise Exception("pileup must be evaluated for failure first")
        return self.failures

    def info(self):
        """Output position stats as VCF INFO field
        """
        if self.position_failed is None:
            raise Exception("pileup must be evaluated for failure first")
        amplicon_totals = ",".join(
            [
                f"{str(self.refs_in_amplicons[amplicon])}/{str(self.alts_in_amplicons[amplicon])}"
                for amplicon in self.amplicon_totals
            ]
        )
        depth_symbol = ">=" if self.total >= 1000 else "="
        info = f"primer={self.refs_in_primer}/{self.alts_in_primer};total{depth_symbol}{self.total};amplicon_overlap={len(self.amplicon_totals)};amplicon_totals={amplicon_totals}"
        return info

    def tsv_row(self):
        # this will be None if the consensus is an 'N'. may not want to skip
        # if self.position_failed is None:
        #    raise Exception("pileup must be evaluated")

        amplicon_totals = ";".join(
            [
                f"{str(self.refs_in_amplicons[amplicon])}:{self.refs_in_forward_strands[amplicon]}:{str(self.alts_in_amplicons[amplicon])}"
                for amplicon in self.amplicon_totals
            ]
        )
        alts = ";".join([f"{k}:{v}" for k, v in self.alt_bases.items()])
        row = "\t".join(
            map(
                str,
                [
                    self.reference_pos,
                    self.ref_base,
                    alts,
                    self.refs_in_primer,
                    self.alts_in_primer,
                    self.total,
                    len(self.amplicon_totals),
                    amplicon_totals,
                ],
            )
        )

        return row

    def __str__(self):
        alts = " ".join([f"{alt}:{count}" for alt, count in self.alt_bases.items()])
        return f"{self.ref_base}: {alts}"


def parse_cigar(ref, query, alignment):
    """Interpret cigar string and query sequence in reference
    coords from mappy (count, op)

    Returns a list of query basecalls per reference position

    ref: AATGG
    qry: AACT-
    (1, "A")
    (2, "AC")
    (3, "T")
    (4, "-")
    (5, "-")
    etc.
    """
    positions = []
    r_pos = alignment.r_st
    cigar = alignment.cigar
    q_pos = alignment.q_st

    for count, op in cigar:
        if op == 0:
            # match/mismatch
            for _ in range(count):
                if len(query) == q_pos:
                    # done? TODO: verify that this captures the full query sequence
                    # print(f"WARNING: invalid cigar string. {len(query)}, index {q_pos}. cigar: {cigar}", file=sys.stderr)
                    break

                positions.append((r_pos, query[q_pos]))
                q_pos += 1
                r_pos += 1

        elif op == 1:
            # insertion
            # positions.append((r_pos, query[q_pos : q_pos + count + 1]))
            q_pos += count

        elif op == 2:
            # deletion
            for _ in range(count):
                positions.append((r_pos, "-"))
                r_pos += 1

        elif op == 3:
            # ref_skip
            pass

        elif op == 4:
            # soft clip
            # q_pos += count
            # may not need to be considered if q_pos offset is set
            # TODO verify
            # may need to consider where softclipping on the 3' end
            pass

        elif op == 5:
            # hard clip
            pass

        else:
            raise Exception(f"invalid cigar op {op}")

    return positions
