"""global config module to govern run parameters"""

import configparser

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
qcovid.primers = [
    "/build/artic-v3.qcovid.tsv",
    "/build/artic-v4.qcovid.tsv",
    "/build/midnight-1200.qcovid.tsv",
]
