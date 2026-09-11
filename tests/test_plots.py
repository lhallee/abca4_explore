from __future__ import annotations

from pathlib import Path

import pandas as pd

from PIL import Image

from abca4_avi.plots import plot_class_counts


def test_plot_is_saved_as_300_dpi_png(tmp_path: Path) -> None:
    path = tmp_path / "plot.png"

    plot_class_counts(pd.DataFrame({"clinical_class": ["P/LP", "B/LB"]}), path)

    with Image.open(path) as image:
        horizontal, vertical = image.info["dpi"]
    assert abs(horizontal - 300) < 1
    assert abs(vertical - 300) < 1
