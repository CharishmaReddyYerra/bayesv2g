"""Command-line interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from bayesv2g import __version__
from bayesv2g.config import default_config_path, load_config, project_root
from bayesv2g.data_sources import (
    fetch_open_targets_loci,
    validate_open_targets_snapshot,
    write_json,
)
from bayesv2g.finemap import credible_set, single_signal_pip
from bayesv2g.pipeline import run_pipeline


def _read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix in {".tsv", ".txt"} or path.name.endswith(".tsv.gz") else ","
    return pd.read_csv(path, sep=separator)


def _run(args: argparse.Namespace) -> int:
    paths = run_pipeline(
        config_path=args.config,
        output_directory=args.output,
        refresh_data=args.refresh,
        verify_tls=not args.insecure,
        posterior_draws=args.draws,
    )
    print(f"Report: {paths['report'].resolve()}")
    print(f"Rankings: {paths['prioritized_genes'].resolve()}")
    print(f"Manifest: {paths['manifest'].resolve()}")
    return 0


def _fetch(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    ids = [item["study_locus_id"] for item in config["demo_loci"]]
    payload = fetch_open_targets_loci(
        config["open_targets_endpoint"], ids, verify_tls=not args.insecure
    )
    validate_open_targets_snapshot(payload, ids, config["open_targets_release"])
    if args.output is None:
        configured_output = Path(
            config.get(
                "open_targets_snapshot",
                "data/snapshots/open_targets_26.09_ldl_loci.json",
            )
        ).expanduser()
        output = (
            configured_output
            if configured_output.is_absolute()
            else project_root() / configured_output
        )
    else:
        output = Path(args.output)
    write_json(payload, output)
    print(f"Wrote {output.resolve()}")
    return 0


def _fine_map(args: argparse.Namespace) -> int:
    frame = _read_table(Path(args.input))
    required = {args.beta_column, args.se_column}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Input is missing columns: {sorted(missing)}")
    frame["log_abf_rank"] = np.nan
    frame["pip"] = single_signal_pip(
        frame[args.beta_column].to_numpy(),
        frame[args.se_column].to_numpy(),
        prior_variance=args.prior_variance,
    )
    frame["log_abf_rank"] = frame["pip"].rank(method="first", ascending=False).astype(int)
    selected = set(credible_set(frame["pip"].to_numpy(), args.coverage).tolist())
    frame["in_credible_set"] = [index in selected for index in range(len(frame))]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    print(f"Wrote {output.resolve()} ({len(selected)} variants in the {args.coverage:.0%} set)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bayesv2g",
        description="Bayesian variant-to-gene prioritization with real public evidence",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="run the complete LDL demonstration")
    run.add_argument("--config", type=Path, default=default_config_path())
    run.add_argument("--output", type=Path, default=project_root() / "results" / "ldl_demo")
    run.add_argument("--draws", type=int, default=None, help="override posterior coefficient draws")
    run.add_argument("--refresh", action="store_true", help="refresh the versioned Open Targets snapshot")
    run.add_argument("--insecure", action="store_true", help="disable TLS verification (only for intercepted local networks)")
    run.set_defaults(func=_run)

    fetch = subparsers.add_parser("fetch", help="refresh the configured Open Targets loci")
    fetch.add_argument("--config", type=Path, default=default_config_path())
    fetch.add_argument(
        "--output",
        type=Path,
        default=None,
    )
    fetch.add_argument("--insecure", action="store_true")
    fetch.set_defaults(func=_fetch)

    finemap = subparsers.add_parser("fine-map", help="single-signal Wakefield-ABF fine-mapping")
    finemap.add_argument("input", type=Path)
    finemap.add_argument("--output", type=Path, required=True)
    finemap.add_argument("--beta-column", default="beta")
    finemap.add_argument("--se-column", default="standard_error")
    finemap.add_argument("--prior-variance", type=float, default=0.04)
    finemap.add_argument("--coverage", type=float, default=0.95)
    finemap.set_defaults(func=_fine_map)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        status = args.func(args)
    except (FileNotFoundError, ValueError, RuntimeError) as error:
        parser.error(str(error))
        return
    raise SystemExit(status)
