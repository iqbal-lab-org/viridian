"""global config module to govern run parameters"""

import configparser
from pathlib import Path

configuration = configparser.ConfigParser()
# configuration.read('config.ini')


class config:
    def __init__(self):
        pass


varifier = config()
qcovid = config()

qcovid.min_template_match_75 = float(0.5)
qcovid.min_coverage = 10
qcovid.variant_freq = float(0.5)

primer_dir = Path(__file__).parent.resolve()

qcovid.primers = [
    str(primer_dir / "../data/artic-v3.qcovid.tsv"),
    str(primer_dir / "../data/artic-v4.qcovid.tsv"),
    str(primer_dir / "../data/midnight-1200.qcovid.tsv"),
]
