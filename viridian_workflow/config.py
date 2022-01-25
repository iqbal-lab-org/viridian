"""global config module to govern run parameters"""

import configparser
from pathlib import Path

configuration = configparser.ConfigParser()
# configuration.read('config.ini')


class config:
    def __init__(self):
        pass


detect_primers = config()
self_qc = config()

detect_primers.min_template_match_75 = float(0.5)

self_qc.min_frs_threshold = 0.7
self_qc.min_depth_threshold = 10
