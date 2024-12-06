# %%
from collections.abc import Iterable
from dataclasses import dataclass, field
import matplotlib.axes
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd
from pathlib import Path
from enum import StrEnum


# %%
class Algo(StrEnum):
    GD = "GD"
    SGD = "SGD"
    ADAM = "ADAM"


class Seed(StrEnum):
    w0 = "0"
    w1 = "0.001"
    w2 = "0.01"
    w3 = "0.1"
    w4 = "1"


class DataSet(StrEnum):
    training = "training"
    testing = "testing"


class Y(StrEnum):
    loss = "loss"
    accuracy = "accuracy"


class X(StrEnum):
    time = "time"
    epoch = "x"


# %%
data_dir = Path("~/Documents/University/Current/Diferansiyel/homework/code/out")


def load_data(algo: Algo, seed: Seed) -> pd.DataFrame:
    return pd.read_csv(data_dir / f"{algo.lower()}_{seed}.csv")


dfs = {(algo, seed): load_data(algo, seed) for algo in Algo for seed in Seed}


# %%
@dataclass
class MyPlot:
    seed: Seed
    dataset: DataSet
    y: Y
    x: X
    fig: Figure = field(init=False, default_factory=lambda: plt.figure())

    def __post_init__(self):
        ax = self.fig.add_subplot(1, 1, 1)
        for algo in Algo:
            df = dfs[(algo, self.seed)]
            ax.plot(self.x, f"{self.dataset}_{self.y}", data=df)
            ax.grid(True)
            if self.y == Y.loss:
                ax.set_yscale("log")

    def save(self):
        self.fig.savefig(f"out/{self.seed}_{self.dataset}_{self.y}_{self.x}.svg")


# %%
plots = [
    MyPlot(seed, dataset, y, x)
    for seed in Seed
    for dataset in DataSet
    for y in Y
    for x in X
]

# %% Save
for p in plots:
    p.save()

# %%
