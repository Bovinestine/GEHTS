# =============================================================================
# main.py — GEHTS preprocessing pipeline entry point
#
# Usage:
#   python main.py <input_dir>           # process CSVs, return AnnData in memory
#   python main.py <input_dir> --save    # also write <input_dir>/merged.h5ad
#
# After this step, run normalization.R to apply SCTransform normalization.
# =============================================================================

from pathlib import Path
import argparse

import processing as pp


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GEHTS preprocessing: load sorted CSV gene-count files → AnnData"
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing sorted raw gene-count CSV files",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Write the merged AnnData to <input_dir>/merged.h5ad",
    )
    args = parser.parse_args()

    print(f"Input directory: {args.input_dir}")
    adata = pp.process_data(args.input_dir, save=args.save)
    print(f"Loaded {adata.n_obs} wells × {adata.n_vars} genes")
    print(f"Drug conditions: {sorted(adata.obs['drug_condition'].unique())}")


if __name__ == "__main__":
    main()
