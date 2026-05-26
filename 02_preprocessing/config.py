# =============================================================================
# config.py — Constants for the GEHTS preprocessing pipeline
#
# Gene panels must match 03_downstream_analysis/R/config.R exactly.
# Lookup tables map DLP barcode numbers to drug names for each dose tier.
# To use a different barcode set, edit DRUG_NAMES and the three LOOKUP_DOSE_*
# lists, then update processing.py / _get_lookup_table() if the filename
# convention for dose detection also changes.
# =============================================================================

# -----------------------------------------------------------------------------
# Gene panels (30-gene GEHTS chip)
# -----------------------------------------------------------------------------
GENES_ANABOLIC = [
    'Acan', 'Sox9', 'Col2a1', 'Matn1', 'Matn3', 'Ucma',
    'Ccnd3', 'Gadd45g', 'Pth1r', 'Gm26633', 'Col27a1',
]

GENES_CATABOLIC = [
    'Mmp3', 'Mmp13', 'Il6', 'Il17b', 'Adamts5', 'Igfbp3',
    'Ccl2', 'Cxcl5', 'Cxcl1', 'Fosl2', 'Tlr2', 'Tnfrsf1b',
]

GENES_HOUSEKEEPING = [
    'Hprt', 'Actb', 'Gapdh', 'B2m', 'Ubc', 'Ppia', 'Rpl23',
]

# -----------------------------------------------------------------------------
# DLP barcode → drug name lookup tables
# Boustrophedon (snake) plate layout; indices align with DRUG_NAMES.
# -----------------------------------------------------------------------------
DRUG_NAMES = [
    'rapamycin', 'SB431542', 'SANT-1', 'LY294002', 'SB525334',
    'KU0063794', 'ALK5 inhibitor IV', 'MK-2206 dihydrochloride',
    'SB203580', 'pamapimod', 'BMS-345541', 'JNK inhibitor V',
    'CAPE', 'XAV',
]

LOOKUP_DOSE_1  = [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14]
LOOKUP_DOSE_01 = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 12]
LOOKUP_DOSE_10 = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,  3]
