import numpy as np
from resonator import Resonator
import os

SAVE_DIR = os.path.join("paper", "images")

R1s = [-0.8, 0.8, -0.8, 0.8]
R2s = [-0.8, -0.8, 1000, 0.8]
L = 0.7
names = [
    "Symmetric Resonator",
    "Convex-Concave Resonator",
    "Concave-Plane Resonator",
    "Unstable Resonator",
]

pos0s = [
    np.array([0.3, -0.1]),
    np.array([-0.1, 0.1]),
    np.array([-0.15, 0]),
    np.array([-0.1, 0.1]),
]

ns = [200, 200, 300, 10]
fig_names = [
    "symmetric",
    "cave_vex",
    "plane_cave",
    "unstable",
]


def save_one(R1, R2, L, name, pos0, n, fig_name):
    resonator = Resonator(R1, R2, L, name=name)
    fig_name = os.path.join(SAVE_DIR, fig_name) + ".png"
    resonator.propogate(
        pos0,
        n,
        show_metadata=True,
        return_fig=False,
        fig_name=fig_name,
        show_fig=False,
    )


def main():
    for i in range(len(names)):
        print(f"Saving {names[i]}.")
        save_one(
            R1s[i],
            R2s[i],
            L,
            names[i],
            pos0s[i],
            ns[i],
            fig_names[i],
        )


if __name__ == "__main__":
    main()
