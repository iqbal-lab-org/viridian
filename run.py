import sys
import shutil
from pathlib import Path

# import tempfile
import json

from viridian_workflow import primers, readstore, run
from viridian_workflow.subtasks import Minimap, Varifier, Viridian

if __name__ == "__main__":
    # set up the pipeline
    platform, work_dir, *fqs = sys.argv[1:]

    # load amplicon sets
    amplicon_sets = []
    for name, tsv in [
        ("artic-v3", "viridian_workflow/amplicon_scheme_data/covid-artic-v3.vwf.tsv"),
        (
            "midnight-1200",
            "viridian_workflow/amplicon_scheme_data/covid-midnight-1200.vwf.tsv",
        ),
        (
            "covid-ampliseq-v1",
            "viridian_workflow/amplicon_scheme_data/covid-ampliseq-v1.vwf.tsv",
        ),
        ("artic-v4.0", "viridian_workflow/amplicon_scheme_data/covid-artic-v4.vwf.tsv"),
    ]:
        amplicon_sets.append(primers.AmpliconSet.from_tsv(tsv, name=name))

    # work_dir = "/tmp/vwf/"
    work_dir = Path(work_dir)
    if work_dir.exists():
        exit(1)
        # print(f"work dir {work_dir} exists, clobbering", file=sys.stderr)
        # shutil.rmtree(work_dir)
    work_dir.mkdir()

    log = {"summary": {"version": "test-0.1", "status": "Interrupted"}}
    try:
        results = run.run_pipeline(work_dir, platform, fqs, amplicon_sets)
        log["results"] = results
        log["summary"]["status"] = "Success"
    except Exception as e:
        log["summary"]["status"] = {"Failure": str(e)}
        print(f"Pipeline failed with exception: {e}", file=sys.stderr)

    with open(work_dir / "log.json", "w") as json_out:
        json.dump(log, json_out, indent=4)
