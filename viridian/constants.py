ALLOWED_TECH = {"illumina", "ont", "iontorrent"}

MSA_CHOICES = {"indel_as_ref", "indel_as_N"}

TECH2FRS = {
    "illumina": 0.5,
    "ont": 0.5,
    "iontorrent": 0.5,
}

CYLON_TECH = {
    "illumina": "illumina",
    "ont": "ont",
    "iontorrent": "ont",
}

ENA_PLATFORMS = {
    "ILLUMINA": "illumina",
    "OXFORD_NANOPORE": "ont",
    "ION_TORRENT": "iontorrent",
}

READ_SAMPLE_MIN_PRIMER_HITS = 25
PRIMER_END_TOLERANCE = 3
SCHEME_ID_PRIMER_WITHIN_END = 10
INDEL_FIX_LENGTH = {
    "ont": 2,
    "iontorrent": 2,
}
MIN_INS_PROPORTION = 0.1

IUPAC = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

IUPAC_REV = {"".join(sorted(v)): k for k, v in IUPAC.items()}
