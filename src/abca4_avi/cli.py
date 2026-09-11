"""Command-line interface for ABCA4 AVI analysis."""

from __future__ import annotations

import argparse
import json

from pathlib import Path

import pandas as pd

from .atlas import annotate_table
from .clinvar import ClinVarClient
from .config import ProjectPaths, load_credentials
from .links import build_variant_url
from .pipeline import create_report, retrieve_clinvar, run_all, score_callset


def _root(value: str | None) -> Path:
    return Path(value or ".").resolve()


def _save_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _clinvar_command(args: argparse.Namespace) -> None:
    root = _root(args.root)
    credentials = load_credentials(root)
    client = ClinVarClient(
        api_key=credentials.ncbi_api_key,
        email=credentials.ncbi_email,
    )
    result = client.search(args.query, retmax=args.max_records or 0)
    _save_json(Path(args.output), result)
    print(f"Saved {result['fetched_count']} ClinVar IDs to {args.output}")


def _annotate_command(args: argparse.Namespace) -> None:
    root = _root(args.root)
    credentials = load_credentials(root)
    if not credentials.alphagenome_api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY is unavailable")
    paths = ProjectPaths.from_root(root)
    paths.create()
    annotated = annotate_table(
        Path(args.input),
        Path(args.output),
        api_key=credentials.alphagenome_api_key,
        checkpoint_path=paths.cache / "alphagenome_scores.jsonl",
        failure_path=paths.logs / "alphagenome_failures.jsonl",
        requests_per_second=args.requests_per_second,
        resume=args.resume,
        max_variants=args.max_variants,
    )
    print(f"Saved {len(annotated)} scored variants to {args.output}")


def _links_command(args: argparse.Namespace) -> None:
    if args.links_action == "variant":
        print(build_variant_url(args.variant, biosample=args.biosample))
        return
    input_path = Path(args.input)
    separator = "\t" if input_path.suffix.lower() == ".tsv" else ","
    frame = (
        pd.read_parquet(input_path)
        if input_path.suffix.lower() == ".parquet"
        else pd.read_csv(input_path, sep=separator)
    )
    frame["atlas_url"] = frame["variant"].map(
        lambda value: build_variant_url(str(value), biosample=args.biosample)
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".parquet":
        frame.to_parquet(output_path, index=False)
    else:
        output_separator = "\t" if output_path.suffix.lower() == ".tsv" else ","
        frame.to_csv(output_path, sep=output_separator, index=False)
    print(f"Saved {len(frame)} Atlas links to {output_path}")


def _analysis_command(args: argparse.Namespace) -> None:
    root = _root(args.root)
    paths = ProjectPaths.from_root(root)
    paths.create()
    if args.analysis_action == "all":
        manifest = run_all(
            root,
            resume=args.resume,
            max_records=args.max_records,
            max_variants=args.max_variants,
            bootstrap_replicates=args.bootstrap_replicates,
        )
    elif args.analysis_action == "clinvar":
        frame = retrieve_clinvar(
            paths, resume=args.resume, max_records=args.max_records
        )
        print(f"Prepared {len(frame)} unique eligible SNVs")
        return
    elif args.analysis_action == "score":
        variants = pd.read_parquet(paths.processed / "variants.parquet")
        frame = score_callset(
            paths,
            variants,
            resume=args.resume,
            max_variants=args.max_variants,
        )
        print(f"Scored {len(frame)} unique SNVs")
        return
    else:
        scored = pd.read_parquet(paths.tables / "variants_scored.parquet")
        manifest = create_report(
            paths, scored, bootstrap_replicates=args.bootstrap_replicates
        )
    print(json.dumps(manifest["record_counts"], sort_keys=True))


def _shared_root(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--root", help="Project root; defaults to the current directory"
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the unified argument parser."""
    parser = argparse.ArgumentParser(prog="abca4-avi")
    subparsers = parser.add_subparsers(dest="command", required=True)

    clinvar = subparsers.add_parser("clinvar")
    clinvar_subparsers = clinvar.add_subparsers(dest="clinvar_action", required=True)
    search = clinvar_subparsers.add_parser("search")
    search.add_argument("--query", required=True)
    search.add_argument("--output", required=True)
    search.add_argument("--max-records", type=int)
    _shared_root(search)
    search.set_defaults(handler=_clinvar_command)

    annotate = subparsers.add_parser("annotate")
    annotate.add_argument("-i", "--input", required=True)
    annotate.add_argument("-o", "--output", required=True)
    annotate.add_argument("--requests-per-second", type=float, default=5.0)
    annotate.add_argument("--max-variants", type=int)
    annotate.add_argument("--include_features", action="store_true")
    annotate.add_argument(
        "--resume", action=argparse.BooleanOptionalAction, default=True
    )
    _shared_root(annotate)
    annotate.set_defaults(handler=_annotate_command)

    links = subparsers.add_parser("links")
    links_subparsers = links.add_subparsers(dest="links_action", required=True)
    variant = links_subparsers.add_parser("variant")
    variant.add_argument("variant")
    variant.add_argument("--biosample")
    variant.set_defaults(handler=_links_command)
    table = links_subparsers.add_parser("table")
    table.add_argument("--input", required=True)
    table.add_argument("--output", required=True)
    table.add_argument("--biosample")
    table.set_defaults(handler=_links_command)

    analysis = subparsers.add_parser("analysis")
    analysis_subparsers = analysis.add_subparsers(dest="analysis_action", required=True)
    for action in ("all", "clinvar", "score", "report"):
        action_parser = analysis_subparsers.add_parser(action)
        action_parser.add_argument(
            "--resume", action=argparse.BooleanOptionalAction, default=True
        )
        action_parser.add_argument("--max-records", type=int)
        action_parser.add_argument("--max-variants", type=int)
        action_parser.add_argument("--bootstrap-replicates", type=int, default=2_000)
        _shared_root(action_parser)
        action_parser.set_defaults(handler=_analysis_command)
    return parser


def main(argv: list[str] | None = None) -> None:
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    args.handler(args)
