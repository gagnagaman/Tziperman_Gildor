#!/usr/bin/env python3
"""
Run the Gildor & Tziperman (2001) Coupled 8-Box Climate Model with Sea-Ice Switch.

Usage:
    python run_simulation.py [--years 50] [--outdir figures]
"""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from seaice_model import SeaIceSwitch8BoxModel


def plot_results(res: dict, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    t = res["time_yr"]
    Teta = res["Teta"] - 273.15      # °C
    T_ocean = res["T_ocean"] - 273.15  # °C
    S_ocean = res["S_ocean"]
    fsi = res["fsi"]
    q_NP = res["q_NP_Sv"]
    q_SP = res["q_SP_Sv"]
    FWAO = res["FWAO_Sv"]

    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True)

    # 1. Atmospheric Potential Temperatures
    ax = axes[0, 0]
    ax.plot(t, Teta[:, 0], "b-", label="Box 1: South Polar (SP)")
    ax.plot(t, Teta[:, 1], "c--", label="Box 2: South Tropical (ST)")
    ax.plot(t, Teta[:, 2], "orange", linestyle="--", label="Box 3: North Tropical (NT)")
    ax.plot(t, Teta[:, 3], "r-", label="Box 4: North Polar (NP)")
    ax.set_ylabel("Atmospheric Temp (°C)", fontsize=10)
    ax.set_title("Atmospheric Potential Temperature (4 Boxes)", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    # 2. Ocean Surface Temperatures
    ax = axes[0, 1]
    ax.plot(t, T_ocean[:, 0], "b-", label="Box 1: SP Surface")
    ax.plot(t, T_ocean[:, 1], "c--", label="Box 2: ST Surface")
    ax.plot(t, T_ocean[:, 2], "orange", linestyle="--", label="Box 3: NT Surface")
    ax.plot(t, T_ocean[:, 3], "r-", label="Box 4: NP Surface")
    ax.set_ylabel("SST (°C)", fontsize=10)
    ax.set_title("Ocean Surface Temperatures", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    # 3. Sea-Ice Fractional Coverage (The Switch)
    ax = axes[1, 0]
    ax.plot(t, fsi[:, 0], "b-", linewidth=1.8, label="South Polar Sea-Ice (fsi 1)")
    ax.plot(t, fsi[:, 3], "r-", linewidth=1.8, label="North Polar Sea-Ice (fsi 4)")
    ax.set_ylabel("Sea-Ice Fraction [0..1]", fontsize=10)
    ax.set_title("Sea-Ice Switch (Polar Coverage)", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    # 4. Ocean Overturning Circulation (Sv)
    ax = axes[1, 1]
    ax.plot(t, q_NP, "r-", linewidth=1.8, label="q34 (North Overturning cell)")
    ax.plot(t, q_SP, "b-", linewidth=1.8, label="q21 (South Overturning cell)")
    ax.set_ylabel("Overturning Transport (Sv)", fontsize=10)
    ax.set_title("Ocean Overturning Circulation", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    # 5. Ocean Surface Salinities
    ax = axes[2, 0]
    ax.plot(t, S_ocean[:, 0], "b-", label="Box 1: SP Salinity")
    ax.plot(t, S_ocean[:, 1], "c--", label="Box 2: ST Salinity")
    ax.plot(t, S_ocean[:, 2], "orange", linestyle="--", label="Box 3: NT Salinity")
    ax.plot(t, S_ocean[:, 3], "r-", label="Box 4: NP Salinity")
    ax.set_xlabel("Time (years)", fontsize=10)
    ax.set_ylabel("Salinity (psu)", fontsize=10)
    ax.set_title("Ocean Surface Salinity", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    # 6. Atmospheric Freshwater Flux to Ocean (P - E)
    ax = axes[2, 1]
    ax.plot(t, FWAO[:, 0], "b-", label="FWAO 1: SP (Net precip)")
    ax.plot(t, FWAO[:, 1], "c--", label="FWAO 2: ST (Net evap)")
    ax.plot(t, FWAO[:, 2], "orange", linestyle="--", label="FWAO 3: NT (Net evap)")
    ax.plot(t, FWAO[:, 3], "r-", label="FWAO 4: NP (Net precip)")
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("Time (years)", fontsize=10)
    ax.set_ylabel("Freshwater Flux (Sv)", fontsize=10)
    ax.set_title("Atmospheric Moisture Flux into Ocean (P - E)", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    plt.suptitle("Gildor & Tziperman (2001) Coupled 8-Box Climate Model with Sea-Ice Switch", fontsize=13)
    plt.tight_layout()

    outpath = os.path.join(outdir, "gildor_tziperman_simulation.png")
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"Saved simulation plot to: {outpath}")


def main():
    parser = argparse.ArgumentParser(description="Run Gildor & Tziperman 8-Box Model")
    parser.add_argument("--years", type=float, default=50.0, help="Simulation length in years")
    parser.add_argument("--outdir", type=str, default="figures", help="Output directory for plots")
    args = parser.parse_args()

    print("=" * 65)
    print(f"Running Gildor & Tziperman (2001) 8-Box Model for {args.years:.0f} years...")
    print("=" * 65)

    model = SeaIceSwitch8BoxModel()
    res = model.run(years=args.years, output_interval_days=10.0)

    print("\nSimulation completed successfully!")
    print(f"Final Polar Atmospheric Temp NP: {res['Teta'][-1, 3] - 273.15:.2f} °C")
    print(f"Final Polar Atmospheric Temp SP: {res['Teta'][-1, 0] - 273.15:.2f} °C")
    print(f"Final Tropical Atmospheric Temp NT: {res['Teta'][-1, 2] - 273.15:.2f} °C")
    print(f"Final Sea-Ice Cover (NP Box 4): {res['fsi'][-1, 3] * 100:.1f} %")
    print(f"Final Sea-Ice Cover (SP Box 1): {res['fsi'][-1, 0] * 100:.1f} %")
    print(f"Final Northern Overturning (q34): {res['q_NP_Sv'][-1]:.2f} Sv")

    plot_results(res, args.outdir)


if __name__ == "__main__":
    main()
