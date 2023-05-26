import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from emodlib.malaria import IntrahostComponent


def run_challenge(duration):
    asexuals = np.zeros(duration)
    gametocytes = np.zeros(duration)
    fevers = np.zeros(duration)
    infects = np.zeros(duration)

    ic = IntrahostComponent.create()
    ic.challenge()

    for t in range(duration):
        ic.update(dt=1)
        asexuals[t] = ic.parasite_density
        gametocytes[t] = ic.gametocyte_density
        fevers[t] = ic.fever_temperature
        infects[t] = ic.infectiousness

    return pd.DataFrame(
        {
            "days": range(duration),
            "parasite_density": asexuals,
            "gametocyte_density": gametocytes,
            "fever_temperature": fevers,
            "infectiousness": infects,
        }
    ).set_index("days")


if __name__ == "__main__":
    IntrahostComponent.set_params()  # default params

    df = run_challenge(duration=300)
    print(df.head(10))

    fig, ax = plt.subplots(1, 1, figsize=(8, 3))

    # df.plot(ax=ax, color=dict(parasite_density="navy", gametocyte_density="darkgreen"))

    # draw parasite/gametocyte timeseries
    df.reset_index(inplace=True)
    df.plot(
        x="days",
        y="parasite_density",
        ax=ax,
        color="firebrick",
        marker="o",
        alpha=0.5,
        legend=False,
    )
    df.plot(
        x="days",
        y="gametocyte_density",
        ax=ax,
        color="darkgreen",
        marker="o",
        alpha=0.5,
        legend=False,
    )
    ax.set(ylim=(1e-3, 1e6), xlim=(0, 80))
    ax.set(yscale="log")

    # draw infectiousness measurement points
    ax.scatter(
        df.days,
        y=[4e5] * len(df),
        s=1 + 100 * df.infectiousness,
        c=100 * df.infectiousness,
        cmap="Greens",
        vmin=0,
        vmax=100,
        lw=0.5,
        edgecolors="darkgreen",
    )

    # draw fever measurement points
    ax.scatter(
        df.days,
        y=[8e5] * len(df),
        s=20 * (df.fever_temperature - 37) + 1,
        c=df.fever_temperature - 37,
        cmap="Reds",
        vmin=0,
        vmax=4,
        lw=0.5,
        edgecolors="firebrick",
    )

    fig.set_tight_layout(True)

    plt.show()
