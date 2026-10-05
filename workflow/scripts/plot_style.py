"""Register Arial for the Swap-seq figures.

Arial is not bundled with matplotlib, so asking for it in `font.sans-serif` alone falls
back to DejaVu Sans. Set SWAPSEQ_ARIAL_TTF to override the search.
"""

import os
import warnings

import matplotlib.font_manager as font_manager

ARIAL_CANDIDATES = [
    os.environ.get("SWAPSEQ_ARIAL_TTF", ""),
    "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/Library/Fonts/Arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "C:/Windows/Fonts/arial.ttf",
]


def register_arial():
    """Make Arial available to matplotlib. Warns and returns False if not found."""
    if any(font.name == "Arial" for font in font_manager.fontManager.ttflist):
        return True

    for path in ARIAL_CANDIDATES:
        if not path or not os.path.exists(path):
            continue
        font_manager.fontManager.addfont(path)
        for bold in (path.replace("Arial.ttf", "Arial_Bold.ttf"),
                     path.replace("Arial.ttf", "Arial-Bold.ttf")):
            if bold != path and os.path.exists(bold):
                font_manager.fontManager.addfont(bold)
        if any(font.name == "Arial" for font in font_manager.fontManager.ttflist):
            return True

    warnings.warn(
        "Arial was not found, so the figures will use the default sans-serif font. "
        "Set SWAPSEQ_ARIAL_TTF to an Arial .ttf to match the manuscript figures.",
        RuntimeWarning)
    return False
