"""
Shared x-axis for Figure 4A (alpha diversity) and 4B (family composition): the
IgE groups ordered by stratum, one colour per stratum, a gap between strata.

Co-author review (2026-10-22): the submitted per-group colours carried no
meaning; colour now encodes the stratum only, and the five stratum colours are
chosen to be distinct from the family colours used in Figure 4B/C
(yellow, green, blue, orange, magenta, grey).
"""

STRATA = [
    ("all sensitized", "#6a3d9a", [("any_ige", "any")]),
    ("number of sensitizations", "#8c564b", [("single", "single"), ("multiple", "multiple")]),
    ("IgE score", "#2a9d8f", [("high_score", "high"), ("medium_score", "medium"),
                              ("low_score", "low")]),
    ("IgE clusters", "#9a8c98", [(f"cluster_{c}", c) for c in "ABCDEFGHI"]),
    ("IgE-negative", "#4d4d4d", [("non_ige", "non-IgE")]),
]
GAP = 0.7   # extra x-space between strata, in box widths


def layout():
    """[(partition, label, colour, x)] in plotting order, with gaps."""
    out, x = [], 0.0
    for _, col, groups in STRATA:
        for part, lab in groups:
            out.append((part, lab, col, x))
            x += 1.0
        x += GAP
    return out
