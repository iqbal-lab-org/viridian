import os

from viridian import utils


def run_cylon(
    reads_dir,
    tech,
    ref_fasta,
    amp_json,
    outdir,
    debug=False,
    max_percent_n=50.0,
    max_percent_amps_fail=50.0,
):
    debug_opt = "--debug" if debug else ""
    results = {}
    command = f"cylon {debug_opt} assemble --reads_per_amp_dir {reads_dir} {tech} {ref_fasta} {amp_json} {outdir}"
    try:
        utils.syscall(command)
    except:
        return results, "Error running cylon"

    json_file = os.path.join(outdir, "run_info.json")

    try:
        results = utils.load_json(json_file)
    except:
        return results, f"Error getting assembly JSON file {json_file}"

    try:
        run_summary = results["run_summary"]
        consensus = run_summary["consensus"]
        assert consensus is not None
    except:
        return results, "No consensus sequence from assembler"

    percent_n = 100.0 * consensus.count("N") / len(consensus)
    if percent_n > max_percent_n:
        return results, f"Too many Ns in assembler consensus: {round(percent_n, 2)}%"

    expect_fasta = os.path.join(outdir, "consensus.final_assembly.fa")
    if not os.path.exists(expect_fasta):
        return results, "No FASTA file made by assembler"

    try:
        total_amps = run_summary["total_amplicons"]
        bad_amps = total_amps - run_summary["successful_amplicons"]
    except:
        return results, "Error counting good and bad amplicons from assembler log"

    if 100 * bad_amps / total_amps >= max_percent_amps_fail:
        return results, "Too many failed amplicons from assembler"

    return results, None
