# %%
import pandas as pd
from sklearn.manifold import TSNE
import plotly.express as px
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


# %% Load
dfs = {algo: pd.read_csv(data_dir / f"weights_{algo.lower()}.csv") for algo in Algo}

# %%
tsne = TSNE(random_state=0, perplexity=15)
gd_projection = tsne.fit_transform(dfs[Algo.GD].loc[:, "0":])

# %%
tsne = TSNE(random_state=0, perplexity=320)
sgd_projection = tsne.fit_transform(dfs[Algo.SGD].loc[:, "0":])

# %%
tsne = TSNE(random_state=0, perplexity=320)
adam_projection = tsne.fit_transform(dfs[Algo.ADAM].loc[:, "0":])

# %%
fig = px.scatter(
    dfs[Algo.GD],
    gd_projection[:, 0],
    gd_projection[:, 1],
    opacity=0.5,
    color=dfs[Algo.GD]["seed"].astype(str),
)
fig

# %%
fig = px.scatter(
    dfs[Algo.SGD],
    sgd_projection[:, 0],
    sgd_projection[:, 1],
    opacity=0.05,
    color=dfs[Algo.SGD]["seed"].astype(str),
)
fig

# %%
fig = px.scatter(
    dfs[Algo.ADAM],
    adam_projection[:, 0],
    adam_projection[:, 1],
    opacity=0.03,
    color=dfs[Algo.ADAM]["seed"].astype(str),
)
fig
