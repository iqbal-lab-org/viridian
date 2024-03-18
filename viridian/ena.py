import logging
import os
import requests

from viridian import constants, utils


def get_ena_metadata(run_id):
    wanted_fields = ["run_accession", "instrument_platform", "fastq_ftp"]
    url = "http://www.ebi.ac.uk/ena/portal/api/filereport?"
    data = {
        "accession": run_id,
        "result": "read_run",
        "fields": ",".join(wanted_fields),
    }
    logging.info(f"Getting metadata from ENA for run {run_id}")
    try:
        r = requests.get(url, data)
    except:
        raise Exception(f"Error querying ENA to get sample from run {run_id} {r.url}")

    if r.status_code != requests.codes.ok:
        raise Exception(
            f"Error requesting data. Error code: {r.status_code}. URL:  {r.url}"
        )

    lines = r.text.rstrip().split("\n")
    if len(lines) != 2:
        lines_str = "\n".join(lines)
        raise Exception(f"Expected exactly 2 lines from ENA request. Got: {lines_str}")

    lines = [x.rstrip().split("\t") for x in lines]
    result = dict(zip(*lines))
    logging.info(f"Metadata: {result}")
    return result


def fastq_ftp_string_to_files(fastq_ftp_str):
    if fastq_ftp_str is None or len(fastq_ftp_str) == 0:
        raise Exception("No reads files reported by ENA. Cannot continue")
    all_files = fastq_ftp_str.split(";")
    if len(all_files) == 0:
        raise Exception("No reads files reported by ENA. Cannot continue")
    files = {"paired": [], "unpaired": []}
    to_match = {}
    for filename in map(os.path.basename, all_files):
        if filename.endswith("_1.fastq.gz") or filename.endswith("_2.fastq.gz"):
            prefix = filename[: -len("_1.fastq.gz")]
            if prefix in to_match:
                assert to_match[prefix] != filename
                assert len(files["paired"]) == 0
                files["paired"] = sorted([to_match[prefix], filename])
                del to_match[prefix]
            else:
                to_match[prefix] = filename
        else:
            files["unpaired"].append(filename)

    files["unpaired"].extend(to_match.values())
    return files


def wanted_reads_files(all_files, platform):
    if platform == "OXFORD_NANOPORE":
        # Just in case there were multiple files and the first two happened
        # to be called foo_1.fastq.gz, foo_2.fastq.gz, assume all reads
        # are unpaired
        return all_files["unpaired"] + all_files["paired"], None, None
    elif platform in ["ILLUMINA", "ION_TORRENT"]:
        if len(all_files["paired"]) == 2:
            return None, *all_files["paired"]
        else:
            return all_files["unpaired"] + all_files["paired"], None, None
    else:
        raise NotImplementedError(f"Unsupported sequencing platform: {platform}")


def download_run(run_id, outdir):
    metadata = get_ena_metadata(run_id)
    # possible platforms are in here:
    # https://github.com/enasequence/schema/blob/master/src/main/resources/uk/ac/ebi/ena/sra/schema/SRA.common.xsd
    # At the time of writing (20220112) these are:
    # LS454, ILLUMINA, HELICOS, ABI_SOLID, COMPLETE_GENOMICS, BGISEQ,
    # OXFORD_NANOPORE, PACBIO_SMRT, ION_TORRENT, CAPILLARY, DNBSEQ.
    # We only want to use ILLUMINA, OXFORD_NANOPORE, ION_TORRENT
    platform = metadata.get("instrument_platform", None)
    if platform is None:
        raise Exception(
            f"No instrument_platform found for ENA run ID {run_id}. Cannot continue"
        )
    if platform not in constants.ENA_PLATFORMS:
        raise Exception(
            f"Cannot use run {run_id} because it has unsupported sequencing platform: {platform}"
        )

    logging.info(f"ENA run {run_id} has instrument_platform {platform}")
    all_reads_files = fastq_ftp_string_to_files(metadata["fastq_ftp"])
    unpaired_fqs, fwd_fq, rev_fq = wanted_reads_files(all_reads_files, platform)
    os.mkdir(outdir)
    utils.syscall(f"enaDataGet -f fastq {run_id}", cwd=outdir)
    got_files = {
        os.path.join(outdir, run_id, x)
        for x in os.listdir(os.path.join(outdir, run_id))
    }
    keep_files = set()
    wanted_files = {}
    logging.info(f"All downloaded reads files: {got_files}")
    run_dir = os.path.join(outdir, run_id)

    if unpaired_fqs is None:
        wanted_files["reads1"] = os.path.join(run_dir, fwd_fq)
        wanted_files["reads2"] = os.path.join(run_dir, rev_fq)
        keep_files = set(wanted_files.values())
    else:
        assert len(unpaired_fqs) > 0
        if len(unpaired_fqs) > 1:
            to_cat = " ".join([os.path.join(run_dir, x) for x in unpaired_fqs])
            outfile = os.path.join(run_dir, "unpaired.fastq.gz")
            utils.syscall(f"cat {to_cat} > {outfile}")
            wanted_files["unpaired_reads"] = outfile
        else:
            wanted_files["unpaired_reads"] = os.path.join(run_dir, unpaired_fqs[0])

        keep_files = {wanted_files["unpaired_reads"]}

    for filename in keep_files:
        if not os.path.exists(filename):
            raise Exception(
                f"Expected to get FASTQ file {filename}, but was not downloaded"
            )

    for filename in got_files.difference(keep_files):
        logging.debug(f"Deleting unwanted reads file {filename}")
        os.unlink(filename)

    logging.debug(f"Kept files downloaded from ENA: {wanted_files}")
    return wanted_files, constants.ENA_PLATFORMS[platform], metadata
