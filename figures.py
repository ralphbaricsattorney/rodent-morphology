import matplotlib.pyplot as plt
import polars as pl
import numpy as np
from scipy.stats import linregress


def make_effect_variability_scatter_plots(
    cohen_ds: pl.DataFrame, variability_ratios: pl.DataFrame, species: str
):
    model = linregress(
        np.array(cohen_ds).flatten(), np.array(variability_ratios).flatten()
    )
    m, b = model.slope, model.intercept
    plt.scatter(cohen_ds, variability_ratios)
    plt.plot(
        np.array(cohen_ds).T,
        m * np.array(cohen_ds).T + b,
        color="red",
        label="OLS best fit",
    )
    plt.xlabel("Cohen's d")
    plt.ylabel("Variability ratios")
    m = "{:.2f}".format(m)
    b = "{:.2f}".format(b)
    plt.text(
        0.95,
        0.95,
        f"m = {m}, b = {b}",
        transform=plt.gca().transAxes,
        horizontalalignment="right",
        verticalalignment="top",
    )
    plt.legend(loc="upper center")
    plt.title(f"Variability-Effect plot for {species}")
    plt.savefig(f"Variability-Effect plot for {species}.png")
    plt.close()


