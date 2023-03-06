import os

from cluster_vcf_records import vcf_file_read

from viridian import utils


def run_varifier(
    cons_fasta,
    ref_fasta,
    outdir,
    min_coord,
    max_coord,
    sanitise_gaps=False,
    indel_fix_length=None,
    debug=False,
):
    debug = "--debug" if debug else ""
    if sanitise_gaps:
        unmasked_cons_fa = os.path.join(outdir, "04.qry_sanitised_gaps.fa")
        sanitise_gaps_opt = "--sanitise_truth_gaps"
    else:
        unmasked_cons_fa = cons_fasta
        sanitise_gaps_opt = ""
    indel_fix = (
        "" if indel_fix_length is None else f"--indel_max_fix_length {indel_fix_length}"
    )
    vcf_header = None
    vcf_records = None
    sequences = {
        "unmasked_consensus": None,
        "msa_unmasked_consensus": None,
        "msa_ref": None,
    }
    command = f"varifier make_truth_vcf {debug} {sanitise_gaps_opt} {indel_fix} --global_align --global_align_min_coord {min_coord} --global_align_max_coord {max_coord} {cons_fasta} {ref_fasta} {outdir}"

    try:
        utils.syscall(command)
    except:
        return sequences, vcf_header, vcf_records, "Error running varifier"

    vcf_file = os.path.join(outdir, "04.truth.vcf")
    try:
        vcf_header, vcf_records = vcf_file_read.vcf_file_to_list(vcf_file)
    except:
        return (
            sequences,
            vcf_header,
            vcf_records,
            f"Error loading initial VCF file {vcf_file}",
        )

    msa_file = os.path.join(outdir, "04.msa")
    if not os.path.exists(msa_file):
        return sequences, vcf_header, vcf_records, f"MSA file not found {msa_file}"

    with open(msa_file) as f:
        lines = [x.rstrip() for x in f]
    if len(lines) != 2:
        return (
            sequences,
            vcf_header,
            vcf_records,
            "MSA file not in expected format (number of lines not equal to 2)",
        )

    sequences["msa_ref"] = lines[0]
    sequences["msa_unmasked_consensus"] = lines[1]

    if len(sequences["msa_ref"]) != len(sequences["msa_unmasked_consensus"]):
        return (
            sequences,
            vcf_header,
            vcf_records,
            "MSA error. Reference and consensus have different lengths",
        )

    try:
        sequences["unmasked_consensus"] = utils.load_single_seq_fasta(
            unmasked_cons_fa
        ).seq
    except:
        return (
            sequences,
            vcf_header,
            vcf_records,
            f"Error loading unmasked consensus sequence from file {unmasked_cons_fa}",
        )

    msa_cons_len = len([x for x in sequences["msa_unmasked_consensus"] if x != "-"])
    if msa_cons_len != len(sequences["unmasked_consensus"]):
        return (
            sequences,
            vcf_header,
            vcf_records,
            "Length mismatch between consensus in MSA and in FASTA",
        )

    return sequences, vcf_header, vcf_records, None
