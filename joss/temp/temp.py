from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp

WD = Path(__file__).parent
SRC_DIR = WD / "npz"
DST_DIR = WD.parent / "graphics"

for name in ["A", "A_min", "A_rec"]:
    mat = sp.csc_array(np.load(SRC_DIR / f"{name}.npz"))
    plt.figure(figsize=(6, 6), edgecolor="black", linewidth=2)
    plt.spy(mat, precision=0, ms=4, c="black")
    plt.axis("off")
    plt.savefig(DST_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close()
