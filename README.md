# Gildor & Tziperman (2001) Coupled 8-Box Model with Sea-Ice Switch (Python Implementation)

A Python implementation of the coupled ocean-atmosphere-sea ice 8-box climate model originally written in Fortran 77 by Hezi Gildor and Eli Tziperman at Harvard University:

> **Gildor, H., & Tziperman, E. (2001).**  
> *A sea ice climate switch for dramatic rapid climate transitions.*  
> **Journal of Geophysical Research: Oceans**, 106(C5), 9117–9132.  
> [DOI: 10.1029/1999JC000192](https://doi.org/10.1029/1999JC000192)  
> Source code: [Harvard Climate Dynamics Group (Eli Tziperman)](https://groups.seas.harvard.edu/climate/eli/Downloads/sea-ice-switch-8-box-model-code-hezi/)

---

## 1. Physical Components

### A. Atmosphere (4 Meridional Boxes)
- **Boxes:** South Polar (SP, 90°S–45°S), South Tropical (ST, 45°S–0°), North Tropical (NT, 0°–45°N), and North Polar (NP, 45°N–90°N).
- **Energy Balance:** Incoming shortwave solar radiation ($H_{\text{in}}$) accounts for surface albedo and cloud albedo; outgoing longwave radiation ($H_{\text{out}}$) obeys Stefan-Boltzmann with parameterised emissivity.
- **Sensible Heat Transport:** Driven down the meridional temperature gradient ($K_T$).
- **Moisture Transport & Clausius-Clapeyron:** Specific humidity $q_{\text{spec}} \propto \exp(-B/\theta)$ transports water vapor poleward ($K_{Wa}$), producing net subtropical evaporation and polar precipitation ($P - E \rightarrow F_{WAO}$).
- **Air-Sea Exchange ($H_{ao}$):** Driven by air-sea temperature difference, strongly insulated when sea ice is present.

### B. Sea-Ice Switch (Polar Boxes 1 & 4)
- **Thermodynamic Growth & Melting:** If ocean water drops below $T_{\text{ice}} = 274.0\text{ K}$ ($1^\circ\text{C}$), sea ice forms.
- **The "Climate Switch":**
  1. **Albedo Feedback:** Surface albedo jumps from open water ($\alpha_s = 0.07$) to sea-ice ($\alpha_{si} = 0.65$).
  2. **Insulation:** Ice cover suppresses air-sea heat flux and evaporation.
  3. **Hydrological Switch:** Shutting off polar evaporation dries the polar atmosphere, terminating precipitation/snowfall.

### C. Ocean (8 Boxes: 4 Surface + 4 Deep)
- **Circulation:** Hydraulic density-driven flow ($q_{21}, q_{23}, q_{34}$, etc.) following Stommel-type dynamics.
- **Seawater Density:** UNESCO (1981) non-linear equation of state (`unesco.py`).
- **Transport & Mixing:** Upwind advection, horizontal and vertical diffusion, plus convective mixing whenever static instability occurs ($\rho_{\text{surface}} > \rho_{\text{deep}}$).

---

## 2. Directory Contents

```
Rannsoknir/SeaIceSwitch8Box/
├── seaice_model.py            # Core coupled 8-box model in Python
├── unesco.py                  # UNESCO (1981) equation of state of seawater
├── run_simulation.py          # Simulation runner with multi-panel plots
├── wood_atmosphere_adapter.py # Adapter to couple this atmosphere with Wood's model
├── original_fortran/          # Complete original Fortran 77 code (from Harvard archive)
│   ├── job_cntr55
│   └── program/
│       ├── atmosphere.f, ocean.f, seaice.f, UNESCO.f, couple.f ...
└── figures/                   # Output PNG figures
```

---

## 3. Running the Model

To run a simulation and generate multi-panel diagnostic plots:

```bash
cd ~/Rannsoknir/SeaIceSwitch8Box
python3 run_simulation.py --years 50
```

This creates `figures/gildor_tziperman_simulation.png` showing:
1. Atmospheric potential temperatures ($\theta_1 \dots \theta_4$)
2. Ocean surface temperatures ($T_1 \dots T_4$)
3. Sea-ice fraction ($f_{si}$ switch in North and South polar boxes)
4. Ocean overturning transports ($q_{NP}, q_{SP}$)
5. Ocean surface salinities ($S_1 \dots S_4$)
6. Atmospheric freshwater flux into ocean ($P - E$ in Sv)

---

## 4. Coupling to the Wood et al. (2019) AMOC Model

In the Wood et al. model (`Rannsoknir/WoodAMOC`):
- North Atlantic freshwater flux $F_N$ and thermocline flux $F_T$ are prescribed constants.
- North Atlantic temperature is a simple diagnostic linear fit: $T_N = \mu q + T_0$.

With the Gildor-Tziperman atmospheric module (`wood_atmosphere_adapter.py`):
1. **Dynamic $F_N$ and $F_T$:** Poleward atmospheric moisture transport $F_{W}$ responds dynamically via Clausius-Clapeyron:
   $$F_N(T) \propto \Delta T \cdot \exp(-B/T)$$
   As the climate warms, the hydrological cycle intensifies, delivering more freshwater to the North Atlantic box.
2. **Sea-Ice Albedo/Insulation:** When the AMOC collapses and the subpolar North Atlantic cools, sea ice can trigger the "switch", reducing local precipitation and stabilizing the state.
