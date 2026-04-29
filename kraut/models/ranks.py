from __future__ import annotations

from typing import Iterator

from .kraken_data import KrakenNode, KrakenReport


CANONICAL_RANKS = ("R", "D", "K", "P", "C", "O", "F", "G", "S")
CLASSIFIED_RANKS = CANONICAL_RANKS[1:]
COUNT_MODES = ("cumulative", "exact")


def canonical_rank(rank_code: str) -> str | None:
    """Return the base taxonomic rank for a Kraken rank code."""
    rank_code = rank_code.strip().upper()
    if not rank_code:
        return None
    if rank_code == "U":
        return "U"

    rank = rank_code[0]
    if rank in CANONICAL_RANKS:
        return rank
    return None


def rank_read_counts(
    report: KrakenReport,
    mode: str = "cumulative",
    include_unclassified: bool = False,
    include_root: bool = False,
) -> dict[str, int]:
    """
    Count reads represented at each canonical rank in a Kraken report.

    In cumulative mode, counts use clade_counts and represent reads classified
    at a rank or below. In exact mode, counts use taxon_counts and represent
    reads assigned directly to nodes of that canonical rank.
    """
    mode = mode.lower()
    if mode not in COUNT_MODES:
        valid_modes = ", ".join(COUNT_MODES)
        raise ValueError(f"mode must be one of: {valid_modes}")

    counts = _empty_rank_counts(include_unclassified, include_root)
    if mode == "exact":
        for node, _ancestor_ranks in _iter_report_nodes(report):
            rank = canonical_rank(node.rank_code)
            if rank in counts:
                counts[rank] += node.taxon_counts
        return counts

    for node, ancestor_ranks in _iter_report_nodes(report):
        rank = canonical_rank(node.rank_code)
        if rank in counts and rank not in ancestor_ranks:
            counts[rank] += node.clade_counts

    return counts


def lowest_rank(report: KrakenReport, mode: str = "cumulative") -> str | None:
    """Return the deepest canonical classified rank with nonzero reads."""
    counts = rank_read_counts(report, mode=mode)
    for rank in reversed(CLASSIFIED_RANKS):
        if counts[rank] > 0:
            return rank
    return None


def _empty_rank_counts(
    include_unclassified: bool,
    include_root: bool,
) -> dict[str, int]:
    counts = {}
    if include_unclassified:
        counts["U"] = 0
    if include_root:
        counts["R"] = 0
    for rank in CLASSIFIED_RANKS:
        counts[rank] = 0
    return counts


def _iter_report_nodes(
    report: KrakenReport,
) -> Iterator[tuple[KrakenNode, tuple[str, ...]]]:
    if report.unclassified is not None:
        yield report.unclassified, ()
    if report.root is not None:
        yield from _iter_node(report.root, ())


def _iter_node(
    node: KrakenNode,
    ancestor_ranks: tuple[str, ...],
) -> Iterator[tuple[KrakenNode, tuple[str, ...]]]:
    yield node, ancestor_ranks

    rank = canonical_rank(node.rank_code)
    if rank is not None:
        ancestor_ranks = ancestor_ranks + (rank,)

    for child in node.children:
        yield from _iter_node(child, ancestor_ranks)
