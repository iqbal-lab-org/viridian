import logging
import os
import random

from viridian import amplicon_schemes, maptools, scheme_id, utils


def process_amplicon_schemes(built_in_amp_schemes=None, tsv_of_amp_schemes=None):
    if built_in_amp_schemes is None and tsv_of_amp_schemes is None:
        logging.info("No primer schemes provided. Using all built in schemes")
        built_in_amp_schemes = list(amplicon_schemes.get_built_in_schemes().keys())
    amplicon_scheme_name_to_tsv = amplicon_schemes.load_list_of_amplicon_sets(
        built_in_names_to_use=built_in_amp_schemes,
        tsv_others_to_use=tsv_of_amp_schemes,
    )
    logging.info(
        f"Processed amplicon scheme files. Amplicon scheme names: {','.join(sorted(list(amplicon_scheme_name_to_tsv.keys())))}"
    )
    return amplicon_scheme_name_to_tsv


def simulate_wgs_reads(ref_seq, outfile, read_length, step=1):
    random.seed(42)
    with open(outfile, "w") as f:
        for i in range(0, len(ref_seq) - read_length, step):
            print(f">{i}", ref_seq[i : i + read_length], sep="\n", file=f)


def sim_reads_one_scheme(
    scheme,
    scheme_name,
    amplicon_scheme_name_to_tsv,
    ref_fasta,
    ref_seq,
    outprefix,
    read_length=None,
    read_depth=50,
    wgs=False,
):
    reads_fasta = f"{outprefix}.reads.fa"
    logging.info(f"{scheme_name}: simulate reads")
    if wgs:
        assert read_length is not None
        simulate_wgs_reads(ref_seq, reads_fasta, read_length)
    else:
        scheme.simulate_reads(ref_seq, reads_fasta, read_length=read_length)

    logging.info(f"{scheme_name}: map reads")
    bam = f"{outprefix}.bam"
    maptools.map_reads(
        bam,
        ref_fasta,
        reads_fasta,
        sample_name=scheme_name,
    )

    logging.info(f"{scheme_name}: analyse mapped reads")
    results_dir = f"{outprefix}.results"
    results, error = scheme_id.analyse_bam(
        bam,
        amplicon_scheme_name_to_tsv,
        results_dir,
        sample_name=scheme_name,
        min_depth_cutoff=1,
    )
    assert error is None
    results_json = os.path.join(results_dir, "results.json")
    utils.write_json(results_json, results)
    return results


def simulate_all_schemes(
    ref_fasta,
    outdir,
    built_in_amp_schemes=None,
    tsv_of_amp_schemes=None,
    read_length=300,
):
    ref_seq = utils.load_single_seq_fasta(ref_fasta)
    logging.info(f"Loaded reference sequence from file {ref_fasta}")
    os.mkdir(outdir)
    amplicon_scheme_name_to_tsv = process_amplicon_schemes(
        built_in_amp_schemes, tsv_of_amp_schemes
    )
    all_results = {}

    for scheme_name, scheme_tsv in amplicon_scheme_name_to_tsv.items():
        logging.info(f"Processing scheme {scheme_name}")
        scheme = scheme_id.Scheme(tsv_file=scheme_tsv)

        for fragment in False, True:
            if fragment:
                outprefix = os.path.join(outdir, f"{scheme_name}.fragment")
                this_read_length = read_length
            else:
                outprefix = os.path.join(outdir, f"{scheme_name}.clean")
                this_read_length = None

            results = sim_reads_one_scheme(
                scheme,
                scheme_name,
                amplicon_scheme_name_to_tsv,
                ref_fasta,
                ref_seq,
                outprefix,
                read_length=this_read_length,
            )
            key = scheme_name + ":" + ("fragmented" if fragment else "full_length")
            all_results[key] = results["scheme_choice"]

    outprefix = os.path.join(outdir, "wgs")
    results = sim_reads_one_scheme(
        None,
        "wgs",
        amplicon_scheme_name_to_tsv,
        ref_fasta,
        ref_seq,
        outprefix,
        wgs=True,
        read_length=read_length,
    )
    all_results["wgs"] = results["scheme_choice"]
    utils.write_json(os.path.join(outdir, "summary.json"), all_results)

    schemes = sorted(list(amplicon_scheme_name_to_tsv.keys()))
    with open(os.path.join(outdir, "summary.tsv"), "w") as f:
        print(
            "Reads",
            *schemes,
            "Best_score",
            "Second_best",
            "Best_two_diff",
            "Best_two_ratio",
            sep="\t",
            file=f,
        )
        for reads, results in sorted(all_results.items()):
            scores = sorted(list(results["scores"].values()))
            print(
                reads,
                *[results["scores"][x] for x in schemes],
                scores[-1],
                scores[-2],
                scores[-1] - scores[-2],
                round(scores[-2] / scores[-1], 3),
                sep="\t",
                file=f,
            )
