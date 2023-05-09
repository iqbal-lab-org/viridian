from collections import defaultdict
import csv
import logging
from operator import itemgetter
import os
import random
import statistics

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pysam

from viridian import constants, utils


LOG_PREFIX = "(scheme_id)"


class Scheme:
    def __init__(self, tsv_file=None, end_tolerance=3):
        self.left_starts = {}
        self.right_ends = {}
        self.left_dists = []
        self.right_dists = []
        self.primer_dist_hist = defaultdict(int)
        self.norm_primer_dist_hist = defaultdict(int)
        self.cumulative_primer_dists = {}
        self.norm_cumul_primer_dists = defaultdict(int)
        self.amplicons = []
        self.score = None
        self.amp_coords = None
        self.amplicon_name_indexes = {}
        self.mean_amp_length = None
        self.end_tolerance = end_tolerance
        self.last_amplicon_end = -1

        if tsv_file is not None:
            try:
                self.load_from_tsv_file(tsv_file)
            except:
                raise Exception(f"Error loading primer scheme from TSV file {tsv_file}")

            self.amp_coords = [(a["start"], a["end"]) for a in self.amplicons]
            self.amp_coords.sort()
            self._calculate_mean_amp_length()

    def load_from_tsv_file(self, tsv_file):
        with open(tsv_file) as f:
            for d in csv.DictReader(f, delimiter="\t"):
                if d["Left_or_right"] not in ["left", "right"]:
                    raise Exception(
                        f"Left_or_right column not left or right. Got: {d['Left_or_right']}"
                    )

                if d["Amplicon_name"] not in self.amplicon_name_indexes:
                    self.amplicons.append(
                        {
                            "name": d["Amplicon_name"],
                            "start": float("inf"),
                            "end": -1,
                            "primers": {"left": [], "right": []},
                        }
                    )
                    self.amplicon_name_indexes[d["Amplicon_name"]] = (
                        len(self.amplicons) - 1
                    )

                amp_index = self.amplicon_name_indexes[d["Amplicon_name"]]
                amp = self.amplicons[amp_index]
                primer_start = int(d["Position"])
                primer_end = primer_start + len(d["Sequence"]) - 1

                if d["Left_or_right"] == "left":
                    primers = amp["primers"]["left"]
                    same = [x for x in primers if x[0] == primer_start]
                    if len(same):
                        same[0][1] = max(primer_end, same[0][1])
                    else:
                        primers.append([primer_start, primer_end])
                        primers.sort()
                    self.left_starts[primer_start] = amp_index
                    amp["start"] = min(amp["start"], primer_start)
                elif d["Left_or_right"] == "right":
                    primers = amp["primers"]["right"]
                    same = [x for x in primers if x[1] == primer_end]
                    if len(same):
                        same[0][0] = min(primer_start, same[0][0])
                    else:
                        primers.append([primer_start, primer_end])
                        primers.sort()

                    self.right_ends[primer_end] = amp_index
                    amp["end"] = max(amp["end"], primer_end)
                    self.last_amplicon_end = max(amp["end"], self.last_amplicon_end)

    def _calculate_mean_amp_length(self):
        left_lengths = [
            self.amp_coords[i + 1][0] - self.amp_coords[i][0]
            for i in range(len(self.amp_coords) - 1)
        ]
        left_lengths.append(self.amp_coords[-1][1] - self.amp_coords[-1][0] + 1)
        right_lengths = [
            self.amp_coords[i][1] - self.amp_coords[i - 1][1]
            for i in range(1, len(self.amp_coords))
        ]
        right_lengths.append(self.amp_coords[0][1] - self.amp_coords[0][0] + 1)
        self.mean_amp_length = statistics.mean(left_lengths + right_lengths)

    def _init_left_distances(self, ref_length):
        starts_list = sorted(list(self.left_starts.keys()))
        self.left_dists = list(range(-starts_list[0], 0))
        for i in range(len(starts_list)):
            end = starts_list[i + 1] if i + 1 < len(starts_list) else ref_length
            self.left_dists.extend(range(end - len(self.left_dists)))

        for s in self.left_starts:
            for i in range(max(0, s - self.end_tolerance), s):
                self.left_dists[i] = 0

    def _init_right_distances(self, ref_length):
        ends_list = sorted(list(self.right_ends.keys()))
        self.right_dists = []
        for i, end in enumerate(ends_list):
            self.right_dists.extend(range(end - len(self.right_dists), -1, -1))
        self.right_dists.extend(range(-1, ends_list[-1] - ref_length, -1))

        for s in self.right_ends:
            for i in range(s, min(s + self.end_tolerance + 1, ref_length)):
                self.right_dists[i] = 0

    def init_distance_lists(self, ref_length):
        self._init_left_distances(ref_length)
        self._init_right_distances(ref_length)

    def make_primer_distance_hist(self, left_hits, right_hits):
        assert (
            len(left_hits)
            == len(right_hits)
            == len(self.left_dists)
            == len(self.right_dists)
        )
        self.primer_dist_hist = defaultdict(int)

        for i in range(len(left_hits)):
            self.primer_dist_hist[self.left_dists[i]] += left_hits[i]
            self.primer_dist_hist[self.right_dists[i]] += right_hits[i]

    def make_normalised_distance_hists(self):
        assert self.mean_amp_length is not None
        assert len(self.primer_dist_hist) > 0
        self.norm_primer_dist_hist = defaultdict(int)
        self.cumulative_primer_dists = defaultdict(int)
        total = 0

        for dist, count in self.primer_dist_hist.items():
            if dist < 0:
                c_dist = min(abs(self.mean_amp_length + dist), self.mean_amp_length)
            else:
                c_dist = min(dist, self.mean_amp_length)

            c_dist = int(100 * c_dist / self.mean_amp_length)
            self.cumulative_primer_dists[c_dist] += count
            dist = int(100 * dist / self.mean_amp_length)
            self.norm_primer_dist_hist[dist] += count
            total += count

        c_list = [self.cumulative_primer_dists.get(0, 0)]
        for i in range(1, 101):
            c_list.append(c_list[-1] + self.cumulative_primer_dists.get(i, 0))
        self.norm_cumul_primer_dists = [int(100 * x / total) for x in c_list]
        self.score = sum(
            [count - i for i, count in enumerate(self.norm_cumul_primer_dists)]
        )

    def count_primer_hits(self, left_hits, right_hits, max_dist=10):
        for amplicon in self.amplicons:
            amplicon["primer_counts"] = {"left": [], "right": []}
            left_primers = amplicon["primers"]["left"]
            right_primers = amplicon["primers"]["right"]

            for i, (start, end) in enumerate(left_primers):
                amplicon["primer_counts"]["left"].append(0)
                range_start = max(0, start - self.end_tolerance)

                if i + 1 < len(left_primers):
                    next_start = left_primers[i + 1][0] - self.end_tolerance
                    range_end = min(start + max_dist, next_start, end)
                else:
                    range_end = min(start + max_dist, len(self.left_dists), end)

                for pos in range(range_start, range_end):
                    amplicon["primer_counts"]["left"][i] += left_hits[pos]

            for i, (start, end) in enumerate(right_primers):
                amplicon["primer_counts"]["right"].append(0)
                range_end = min(end + self.end_tolerance + 1, len(self.right_dists))

                if i > 0:
                    previous_end = right_primers[i - 1][1] + self.end_tolerance + 1
                    range_start = max(end - max_dist, previous_end)
                else:
                    range_start = max(0, end - max_dist)

                for pos in range(range_start, range_end):
                    amplicon["primer_counts"]["right"][i] += right_hits[pos]

    def to_json_dict(self, debug=False):
        json_dict = {
            "amplicons": self.amplicons,
        }

        return json_dict

    def to_cleaned_amp_list_and_cylon_dict(self, min_primer_hits):
        list_out = []
        cylon_out = {}
        for amp_dict in self.amplicons:
            new_amp = {
                "name": amp_dict["name"],
                "primers": {"left": [], "right": []},
                "excluded_primers": {"left": [], "right": []},
            }
            all_primers = amp_dict["primers"]

            for l_or_r in ["left", "right"]:
                counts = amp_dict["primer_counts"][l_or_r]
                ok_indexes = {i for i, c in enumerate(counts) if c >= min_primer_hits}
                if len(ok_indexes) == 0:
                    ok_indexes = set(range(len(counts)))

                for i in range(len(counts)):
                    key = "primers" if i in ok_indexes else "excluded_primers"
                    new_amp[key][l_or_r].append(
                        {
                            "start": all_primers[l_or_r][i][0],
                            "end": all_primers[l_or_r][i][1],
                            "read_count": counts[i],
                        }
                    )

            new_amp["start"] = min([x["start"] for x in new_amp["primers"]["left"]])
            new_amp["end"] = max([x["end"] for x in new_amp["primers"]["right"]])
            new_amp["primers"]["left"].sort(key=itemgetter("start"))
            new_amp["primers"]["right"].sort(key=itemgetter("end"))
            list_out.append(new_amp)
            cylon_out[amp_dict["name"]] = {
                "start": new_amp["start"],
                "end": new_amp["end"],
                "left_primer_end": new_amp["primers"]["left"][0]["end"],
                "right_primer_start": new_amp["primers"]["right"][-1]["start"],
            }

        list_out.sort(key=itemgetter("start"))
        return list_out, cylon_out

    def simulate_reads(
        self, ref_seq, outfile, read_length=None, read_depth=50, min_read_length=50
    ):
        if read_length is not None:
            random.seed(42)

        with open(outfile, "w") as f:
            for (start, end) in self.amp_coords:
                if read_length is None:
                    print(f">{start}_{end}", file=f)
                    print(ref_seq[start : end + 1], file=f)
                else:
                    for i in range(read_depth):
                        while True:
                            middle = random.randint(
                                start + min_read_length, end - min_read_length
                            )
                            this_start = max(start, middle - int(read_length / 2))
                            this_end = min(end, this_start + read_length)
                            if this_end - this_start >= min_read_length:
                                print(
                                    f">{start}_{end}_{i}_{this_start}_{this_end-1}",
                                    file=f,
                                )
                                print(ref_seq[this_start:this_end], file=f)
                                break


def depth_increments_to_pileup_stats(depth_increments, min_depth_cutoff=20):
    depth_per_position = [depth_increments[0]]
    depth_hist = {depth_increments[0]: 1}
    x_cov_cutoffs = sorted(list({min_depth_cutoff, 1, 2, 5, 10, 15, 20, 50, 100}))
    depth_at_least_x = {x: 0 for x in x_cov_cutoffs}

    for x in depth_increments[1:-1]:
        new_depth = depth_per_position[-1] + x
        depth_per_position.append(new_depth)
        depth_hist[new_depth] = depth_hist.get(new_depth, 0) + 1
        for cov in x_cov_cutoffs:
            if new_depth < cov:
                break
            depth_at_least_x[cov] += 1

    ref_length = len(depth_per_position)
    percent_at_least = {
        k: round(100 * v / ref_length, 2) for k, v in depth_at_least_x.items()
    }

    return {
        "depth_hist": depth_hist,
        "depth_per_position": depth_per_position,
        "depth_at_least": depth_at_least_x,
        "percent_at_least_x_depth": percent_at_least,
        "mean_depth": round(statistics.mean(depth_per_position), 2),
        "mode_depth": statistics.mode(depth_per_position),
        "median_depth": statistics.median(depth_per_position),
    }


def parse_bam(filename, min_depth_cutoff=20):
    aln_file = pysam.AlignmentFile(filename)
    assert len(aln_file.lengths) == 1
    ref_length = aln_file.lengths[0]
    lefts = [0] * ref_length
    rights = [0] * ref_length
    depth_increments = [0] * (ref_length + 1)
    read_counts = {
        "total_reads": 0,
        "unpaired_reads": 0,
        "reads1": 0,
        "reads2": 0,
        "mapped": 0,
    }

    for read in aln_file:
        if read.is_secondary or read.is_supplementary:
            continue

        read_counts["total_reads"] += 1

        if read.is_paired:
            if read.is_read1:
                read_counts["reads1"] += 1
            else:
                read_counts["reads2"] += 1
        else:
            read_counts["unpaired_reads"] += 1

        if read.is_unmapped:
            continue

        read_counts["mapped"] += 1

        if read.is_paired:
            if not read.is_proper_pair:
                continue

            if read.template_length > 0:
                lefts[read.reference_start] += 1
            else:
                rights[read.reference_end - 1] += 1
        else:
            lefts[read.reference_start] += 1
            rights[read.reference_end - 1] += 1

        depth_increments[read.reference_start] += 1
        depth_increments[read.reference_end] -= 1

    pileup = depth_increments_to_pileup_stats(
        depth_increments, min_depth_cutoff=min_depth_cutoff
    )
    return ref_length, lefts, rights, pileup, read_counts


def depth_per_position_plot(
    depths, outfile, amp_coords=None, title=None, min_depth_cutoff=20
):
    x = list(range(1, len(depths) + 1))
    plt.clf()
    fig, ax = plt.subplots()
    plt.plot(
        [min(x), max(x)],
        [min_depth_cutoff, min_depth_cutoff],
        ":",
        color="black",
        linewidth=0.5,
    )
    plt.plot(x, depths, color="black", linewidth=0.75)
    y_bad = [y if y < min_depth_cutoff else None for y in depths]
    plt.plot(x, y_bad, color="red", linewidth=0.75)

    if amp_coords is not None:
        colours = ["blue", "red"]
        for i, (start, end) in enumerate(amp_coords):
            ax.axvspan(start, end, facecolor=colours[i % 2], alpha=0.15)

    if title is not None:
        ax.set_title(title)
    ax.set_xlabel("Position in genome")
    ax.set_ylabel("Read depth")
    f = plt.gcf()
    f.set_size_inches(15, 4)
    plt.savefig(outfile)
    plt.close()


def cumulative_score_plot(schemes, outfile, title=None):
    colormap = mpl.colormaps["tab10"].colors
    scheme_names = sorted(list(schemes.keys()))

    if len(scheme_names) > len(colormap):
        colormap = mpl.colormaps["tab20c"].colors
    if len(scheme_names) > len(colormap):
        logging.warn(
            f"{LOG_PREFIX} Plotting {len(scheme_names)} using {len(colormap)} colors. Have to recycle colors"
        )

    colors = {name: colormap[i % len(colormap)] for i, name in enumerate(scheme_names)}
    plt.clf()
    fig, ax = plt.subplots()
    legend_handles = []
    x = list(range(101))
    for scheme_name in scheme_names:
        colour = colors[scheme_name]
        plt.plot(x, schemes[scheme_name].norm_cumul_primer_dists, color=colour)
        legend_handles.append(mpatches.Patch(color=colour, label=scheme_name))

    plt.plot([0, 100], [0, 100], "--", color="grey")
    ax.set_xlabel("Distance to amplicon end (as % of mean amplicon length)")
    ax.set_ylabel("Cumulative percent of fragments")
    if title is not None:
        ax.set_title(title)
    ax.legend(handles=legend_handles)
    f = plt.gcf()
    f.set_size_inches(7, 7)
    plt.savefig(outfile)
    plt.close()


def get_scores_from_schemes(schemes):
    scores = {}
    best_score = None
    for scheme_name, scheme in schemes.items():
        scores[scheme_name] = scheme.score
        if scheme.score is not None:
            best_score = (
                scheme.score if best_score is None else max(best_score, scheme.score)
            )
    best_schemes = sorted(
        list({x for x in scores if scores[x] == best_score and scores[x] is not None})
    )
    return {
        "scores": scores,
        "best_schemes": best_schemes,
        "best_score": best_score,
        "best_scheme": None if len(best_schemes) == 0 else best_schemes[0],
    }


def analyse_bam(
    bam_file,
    scheme_tsvs,
    outdir,
    end_tolerance=3,
    sample_name=None,
    debug=False,
    max_primer_dist=10,
    min_depth_cutoff=20,
    min_percent_genome_cutoff=50.0,
    min_primer_hits=constants.READ_SAMPLE_MIN_PRIMER_HITS,
):
    json_dict = {
        "amplicons": None,
        "cylon_amplicons": None,
        "scheme_choice": None,
        "stats": None,
        "read_counts": None,
    }
    logging.info(
        f"{LOG_PREFIX} Analysing BAM file {bam_file}. Gathering read depths and primer matches"
    )
    ref_length, left_primer_hits, right_primer_hits, pileup, read_counts = parse_bam(
        bam_file, min_depth_cutoff=min_depth_cutoff
    )
    json_dict["stats"] = pileup
    json_dict["read_counts"] = read_counts

    covered = pileup["percent_at_least_x_depth"][min_depth_cutoff]
    logging.info(
        f"{LOG_PREFIX} {covered} percent of genome has at least {min_depth_cutoff}X read depth (require at least {min_percent_genome_cutoff} percent)"
    )
    if covered < min_percent_genome_cutoff:
        return (
            json_dict,
            f"Not enough read coverage. Got {covered}% of the genome with at least {min_depth_cutoff}X read depth, but need at least {min_percent_genome_cutoff}% of the genome with at least {min_depth_cutoff}X read depth",
        )
    if not len(left_primer_hits) == len(right_primer_hits) == ref_length:
        raise Exception(
            "Error gathering stats from BAM file {bam_file}. Cannot continue"
        )

    schemes = {}
    for scheme_name, scheme_tsv in scheme_tsvs.items():
        logging.info(f"{LOG_PREFIX} Analysing amplicon scheme {scheme_name}")
        logging.debug(f"{LOG_PREFIX} {scheme_name} Load TSV file {scheme_tsv}")
        scheme = Scheme(tsv_file=scheme_tsv, end_tolerance=end_tolerance)
        if scheme.last_amplicon_end > ref_length:
            return (
                json_dict,
                f"Scheme does not match reference genome. Scheme {scheme_name} has final amplicon ending at position {scheme.last_amplicon_end}, but the reference length is only {ref_length}",
            )
        logging.debug(f"{LOG_PREFIX} {scheme_name} initialising primer distance lists")
        scheme.init_distance_lists(ref_length)
        logging.debug(f"{LOG_PREFIX} {scheme_name} making primer distance histogram")
        scheme.make_primer_distance_hist(left_primer_hits, right_primer_hits)
        logging.debug(
            f"{LOG_PREFIX} {scheme_name} making normalised primer distance histograms etc"
        )
        scheme.make_normalised_distance_hists()
        logging.info(f"{LOG_PREFIX} {scheme_name} score: {scheme.score}")
        schemes[scheme_name] = scheme

    os.mkdir(outdir)
    score_plot = os.path.join(outdir, "score_plot.pdf")
    logging.debug(f"{LOG_PREFIX} Making cumulative score plot")
    cumulative_score_plot(schemes, score_plot, title=sample_name)
    scores = get_scores_from_schemes(schemes)
    json_dict["scheme_choice"] = scores
    best_scheme_name = scores["best_scheme"]
    best_scheme = schemes.get(best_scheme_name, None)

    per_pos_depth_plot = os.path.join(outdir, "depth_across_genome.pdf")
    logging.debug(f"{LOG_PREFIX} Making plot of read depth across genome")
    depth_per_position_plot(
        pileup["depth_per_position"],
        per_pos_depth_plot,
        amp_coords=None if best_scheme is None else best_scheme.amp_coords,
        title=sample_name,
        min_depth_cutoff=min_depth_cutoff,
    )

    if best_scheme is None:
        raise Exception("Error getting best amplicon scheme. Cannot continue")
    logging.info(f"{LOG_PREFIX} Best scheme: {best_scheme_name}")
    logging.info(f"{LOG_PREFIX} Counting primer matches for scheme {best_scheme_name}")
    best_scheme.count_primer_hits(
        left_primer_hits, right_primer_hits, max_dist=max_primer_dist
    )
    amp_list, cylon_dict = best_scheme.to_cleaned_amp_list_and_cylon_dict(
        min_primer_hits
    )
    json_dict["amplicons"] = amp_list
    json_dict["cylon_amplicons"] = cylon_dict
    if debug:
        json_dict["debug_scheme"] = best_scheme.to_json_dict()
        json_out = os.path.join(outdir, "debug.json")
        logging.debug(f"{LOG_PREFIX} Writing scheme id JSON file {json_out}")
        utils.write_json(json_out, json_dict)
    logging.info(f"{LOG_PREFIX} Finished analysing amplicon schemes vs reads")
    return json_dict, None
