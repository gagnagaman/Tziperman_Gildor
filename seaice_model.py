"""
Gildor & Tziperman (2001) Coupled 8-Box Climate Model with Sea-Ice Switch.

Python translation of the coupled ocean-atmosphere-sea ice 8-box model:
    Gildor, H., & Tziperman, E. (2001).
    A sea ice climate switch for dramatic rapid climate transitions.
    Journal of Geophysical Research: Oceans, 106(C5), 9117-9132.
    https://doi.org/10.1029/1999JC000192

Architecture:
    - Atmosphere: 4 meridional boxes (Southern Polar, Southern Tropical/Midlat,
      Northern Tropical/Midlat, Northern Polar) with radiative balance, Clausius-
      Clapeyron humidity, and meridional sensible + latent heat transport.
    - Sea Ice: Dynamic sea-ice switch in polar boxes (thermodynamic growth, melt,
      ice-albedo feedback, and insulation of air-sea exchange).
    - Ocean: 8 boxes (4 surface + 4 deep) with density-driven overturning,
      upwind advection, horizontal/vertical diffusion, and convective mixing.
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
from unesco import unesco_full_density


class SeaIceSwitch8BoxModel:
    """
    Coupled Ocean-Atmosphere-Sea Ice 8-Box Model.
    """

    def __init__(self, **kwargs):
        # Time step (seconds): default 21600 s = 6 hours
        self.dt = float(kwargs.get("dt", 21600.0))

        # Box dimensions:
        # 4 zones: 1 (SP, 90S-45S), 2 (ST, 45S-EQ), 3 (NT, EQ-45N), 4 (NP, 45N-90N)
        self.L = np.array([4.15e6, 1.0e7, 1.0e7, 4.15e6], dtype=float)
        self.width = np.array([1.8e7, 1.8e7, 1.8e7, 1.8e7], dtype=float)
        self.fl = np.array([0.50, 0.20, 0.35, 0.50], dtype=float)  # Land fraction
        self.fo = 1.0 - self.fl                                    # Ocean fraction
        self.Aarea = self.L * self.width
        self.Oarea = self.fo * self.Aarea

        # Ocean depths (m): 4 surface (400 m) + 4 deep (4000 m)
        self.D = np.array([400.0, 400.0, 400.0, 400.0,
                           4000.0, 4000.0, 4000.0, 4000.0], dtype=float)
        
        # Ocean volumes (m^3)
        self.vol = np.zeros(8, dtype=float)
        for i in range(4):
            self.vol[i] = self.L[i] * self.D[i] * (self.fo[i] * self.width[i])
            self.vol[i + 4] = self.L[i] * self.D[i + 4] * (self.fo[i] * self.width[i])

        # Physical constants
        self.rhoW = float(kwargs.get("rhoW", 1028.0))       # Seawater density (kg/m^3)
        self.G = 9.8                                        # Gravity (m/s^2)
        self.sigma = 5.6701e-8                              # Stefan-Boltzmann (W/m^2/K^4)
        self.RR = 287.04                                    # Gas constant dry air (J/kg/K)
        self.CpA = 1004.0                                   # Air specific heat (J/kg/K)
        self.CpW = 4180.0                                   # Water specific heat (J/kg/K)
        self.Ps = 1.0e5                                     # Surface pressure (Pa)
        self.epsilon = 0.622                                # Rd / Rv or moisture ratio
        self.CwRhoW = 4.29e6                                # Cw * RhoW (J/m^3/K)
        self.Lice = 3.34e5                                  # Latent heat of fusion of ice (J/kg)
        self.Rhoice = 917.0                                 # Ice density (kg/m^3)
        self.rhoILice = self.Rhoice * self.Lice
        self.rhoIw = self.Rhoice / self.rhoW

        # Atmospheric parameters
        self.Gama0 = float(kwargs.get("Gama0", 1.15e-8))    # Air-sea heat exchange (s^-1)
        self.albedoL = float(kwargs.get("albedoL", 0.20))   # Land albedo
        self.albedoI = float(kwargs.get("albedoI", 0.90))   # Land-ice albedo
        self.albedoSI = float(kwargs.get("albedoSI", 0.65)) # Sea-ice albedo
        self.albedoS = float(kwargs.get("albedoS", 0.07))   # Open ocean albedo
        self.albedoC = float(kwargs.get("albedoC", 0.30))   # Cloud albedo
        self.Solar = np.array([244.67, 381.95, 381.95, 244.67], dtype=float) # (W/m^2)
        self.KT = float(kwargs.get("KT", 15.5e20))          # Atmospheric sensible heat transport
        self.KWa = float(kwargs.get("KWa", 5.4e13))         # Atmospheric moisture transport
        self.KWa2 = float(kwargs.get("KWa2", 3.3e8))        # Local polar moisture coefficient
        self.Plw0 = np.array([0.59, 0.53, 0.51, 0.73], dtype=float) # Background atmospheric emissivity
        self.RH = np.array([0.70, 0.70, 0.70, 0.70], dtype=float)   # Relative humidity
        self.ac1 = (2.0**(self.RR / self.CpA)) * self.G / (self.Ps * self.CpA)

        # Sea-ice parameters
        self.Tice = float(kwargs.get("Tice", 274.0))        # Freezing threshold (K ~ 1 °C)
        self.Di = np.array([1.5, 4.0, 4.0, 2.0], dtype=float) # Sea-ice reference thickness (m)
        self.fsi = np.zeros(4, dtype=float)                 # Sea-ice area fraction [0..1]
        self.sivol = np.zeros(4, dtype=float)               # Sea-ice volume (m^3)
        self.Idepth = np.copy(self.Di)                      # Actual ice thickness (m)
        self.fi = np.array([0.5, 0.0, 0.0, 0.15], dtype=float) # Land ice fraction (fixed or prescribed)

        # Ocean parameters
        self.r = float(kwargs.get("r", 2.6e-4))             # Ocean friction
        self.uppersection = float(kwargs.get("uppersection", 4.5e9))
        self.lowersection = float(kwargs.get("lowersection", 4.5e10))
        self.Kv = np.array([4.2e-8, 4.2e-8, 4.2e-8, 4.2e-8], dtype=float) # Vertical diffusion (m^2/s)
        self.Kh = np.array([1.5e-4, 1.5e-4, 1.5e-4, 3.5e-4, 3.5e-4, 3.5e-4], dtype=float) # Horiz diff (m^2/s)
        self.Sref = 35.0                                    # Reference salinity (psu)

        # Hydraulic transport coefficients (converting density gradients to transport in m^3/s)
        self.c21 = (self.D[4] * self.G / (self.r * self.rhoW * (self.D[0] + self.D[4]) * (self.L[0] + self.L[1]))) * self.uppersection
        self.c23 = (self.D[5] * self.G / (self.r * self.rhoW * (self.D[1] + self.D[5]) * (self.L[1] + self.L[2]))) * self.uppersection
        self.c34 = (self.D[7] * self.G / (self.r * self.rhoW * (self.D[3] + self.D[7]) * (self.L[2] + self.L[3]))) * self.uppersection
        self.c56 = self.c21 * (self.D[0] / self.D[4]) * (self.lowersection / self.uppersection)
        self.c76 = self.c23 * (self.D[1] / self.D[5]) * (self.lowersection / self.uppersection)
        self.c87 = self.c34 * (self.D[3] / self.D[7]) * (self.lowersection / self.uppersection)

        # Initial state:
        # Atmosphere potential temperature (K)
        self.Teta = np.array([273.0, 283.0, 283.0, 273.0], dtype=float)
        # Ocean temperature (K) and salinity (psu)
        self.T_ocean = np.full(8, 278.0, dtype=float)
        self.T_ocean[0] = 271.5  # South Polar surface
        self.T_ocean[1] = 295.0  # South Tropical surface
        self.T_ocean[2] = 296.0  # North Tropical surface
        self.T_ocean[3] = 272.5  # North Polar surface
        self.T_ocean[4:] = 275.0 # Deep boxes

        self.S_ocean = np.array([34.0, 35.8, 35.8, 34.3,
                                 34.6, 34.7, 34.8, 34.7], dtype=float)

        self.rho = np.zeros(8, dtype=float)
        self.update_density()

        # Diagnostics storage
        self.q_transports = {}
        self.FWAO = np.zeros(4, dtype=float) # Fresh water flux into ocean (m^3/s or Sv)
        self.Hao = np.zeros(4, dtype=float)  # Air-sea sensible heat flux (W)

    def update_density(self):
        """Update seawater density for all 8 boxes using UNESCO (1981)."""
        for i in range(8):
            t_c = self.T_ocean[i] - 273.15
            p_depth = 200.0 if i < 4 else 2000.0
            self.rho[i] = unesco_full_density(self.S_ocean[i], t_c, p_depth)

    def step_atmosphere(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute atmospheric moisture, heat transport, radiation, and air-sea flux.
        Returns:
            FWAO: Net freshwater flux into ocean for boxes 1..4 (m^3/s)
            Hao: Air-to-ocean sensible heat flux for boxes 1..4 (W)
        """
        # Specific humidity via Clausius-Clapeyron
        A_cc = 2.53e11  # Pa
        B_cc = 5.42e3   # K
        qsat = A_cc * np.exp(-B_cc / self.Teta)
        q_spec = self.RH * self.epsilon * qsat / self.Ps

        # Meridional sensible heat transport (FT) in J/s
        FT1122 = self.KT * abs(self.Teta[1] - self.Teta[0]) / (self.L[0] + self.L[1])
        FT2233 = self.KT * abs(self.Teta[2] - self.Teta[1]) / (self.L[1] + self.L[2])
        FT3344 = self.KT * abs(self.Teta[2] - self.Teta[3]) / (self.L[2] + self.L[3])

        # Directional fluxes (K/s)
        sign_12 = 1.0 if self.Teta[0] <= self.Teta[1] else -1.0
        FT21 = sign_12 * FT1122 * self.ac1 / self.Aarea[0]
        FT12 = -sign_12 * FT1122 * self.ac1 / self.Aarea[1]

        sign_23 = 1.0 if self.Teta[1] <= self.Teta[2] else -1.0
        FT32 = sign_23 * FT2233 * self.ac1 / self.Aarea[1]
        FT23 = -sign_23 * FT2233 * self.ac1 / self.Aarea[2]

        sign_34 = 1.0 if self.Teta[2] >= self.Teta[3] else -1.0
        FT43 = -sign_34 * FT3344 * self.ac1 / self.Aarea[2]
        FT34 = sign_34 * FT3344 * self.ac1 / self.Aarea[3]

        # Meridional moisture fluxes (m^3/s)
        if self.Teta[0] <= self.Teta[1]:
            FW1122 = self.KWa * (self.Teta[1] - self.Teta[0]) * q_spec[0] / (0.25 * (self.L[0] + self.L[1]))
            FW21 = 0.5 * FW1122
            FW12 = -0.5 * FW1122
        else:
            FW1122 = self.KWa * (self.Teta[0] - self.Teta[1]) * q_spec[1] / (0.25 * (self.L[0] + self.L[1]))
            FW21 = -0.5 * FW1122
            FW12 = 0.5 * FW1122

        if self.Teta[1] <= self.Teta[2]:
            FW2233 = self.KWa * (self.Teta[2] - self.Teta[1]) * q_spec[1] / (0.25 * (self.L[1] + self.L[2]))
            FW32 = FW2233
            FW23 = -FW2233
        else:
            FW2233 = self.KWa * (self.Teta[1] - self.Teta[2]) * q_spec[2] / (0.25 * (self.L[1] + self.L[2]))
            FW32 = -FW2233
            FW23 = FW2233

        if self.Teta[2] >= self.Teta[3]:
            FW3344 = self.KWa * (self.Teta[2] - self.Teta[3]) * q_spec[3] / (0.25 * (self.L[2] + self.L[3]))
            FW43 = -0.5 * FW3344
            FW34 = 0.5 * FW3344
        else:
            FW3344 = self.KWa * (self.Teta[3] - self.Teta[2]) * q_spec[2] / (0.25 * (self.L[2] + self.L[3]))
            FW43 = 0.5 * FW3344
            FW34 = -0.5 * FW3344

        # Freshwater convergence into ocean boxes (m^3/s)
        self.FWAO[0] = FW21
        self.FWAO[1] = FW12 + FW32
        self.FWAO[2] = FW23 + FW43
        self.FWAO[3] = FW34

        # Albedo and radiation in each box
        albedo = np.zeros(4)
        Hin = np.zeros(4)
        Hout = np.zeros(4)
        Hao = np.zeros(4)

        for i in range(4):
            # Surface albedo weighted by land, sea-ice, and open water
            albedo[i] = (self.fl[i] * (1.0 - self.fi[i]) * self.albedoL
                         + self.fi[i] * self.fl[i] * self.albedoI
                         + self.fsi[i] * self.fo[i] * self.albedoSI
                         + (1.0 - self.fsi[i]) * self.fo[i] * self.albedoS)
            # Incoming shortwave
            Hin[i] = self.Solar[i] * (1.0 - albedo[i]) * (1.0 - self.albedoC)
            # Outgoing longwave
            Hout[i] = self.sigma * self.Plw0[i] * (self.Teta[i]**4)

            # Air-sea sensible heat flux
            dT_ao = self.Teta[i] - self.T_ocean[i]
            # Open water heat flux
            flux_open = self.CwRhoW * self.D[i] * self.fo[i] * (1.0 - self.fsi[i]) * self.Gama0 * dT_ao
            # Insulated flux through sea ice
            ice_thick = max(self.Idepth[i], 0.1)
            flux_ice = (self.CwRhoW * self.D[i] * self.fo[i] * self.fsi[i] * 0.05 * self.Gama0 * dT_ao) / ice_thick
            Hao[i] = flux_open + flux_ice

        self.Hao = Hao

        # Potential temperature rate of change (K/s)
        dTeta_dt = np.zeros(4)
        dTeta_dt[0] = FT21 + self.ac1 * (-Hao[0] - Hout[0] + Hin[0])
        dTeta_dt[1] = FT12 + FT32 + self.ac1 * (-Hao[1] - Hout[1] + Hin[1])
        dTeta_dt[2] = FT43 + FT23 + self.ac1 * (-Hao[2] - Hout[2] + Hin[2])
        dTeta_dt[3] = FT34 + self.ac1 * (-Hao[3] - Hout[3] + Hin[3])

        # Step atmospheric potential temperature forward
        self.Teta += dTeta_dt * self.dt
        return self.FWAO, self.Hao

    def step_seaice(self):
        """
        Thermodynamic sea-ice switch in polar boxes (1 and 4).
        Updates sea-ice volume (sivol), area coverage fraction (fsi), and thickness (Idepth).
        """
        for i in [0, 3]:  # South polar (0) and North polar (3)
            T_water = self.T_ocean[i]
            if T_water < self.Tice:
                # Freezing: Heat deficit restores ocean toward Tice
                Havail = 3.8e-7 * self.CwRhoW * self.vol[i] * (self.Tice - T_water)
                new_ice_vol = (Havail * self.dt) / self.rhoILice
                self.sivol[i] += new_ice_vol

                # Compute fractional coverage or thickness
                if (self.sivol[i] / self.Di[i]) < self.Oarea[i]:
                    self.fsi[i] = self.sivol[i] / (self.Di[i] * self.Oarea[i])
                    self.Idepth[i] = self.Di[i]
                else:
                    self.fsi[i] = 1.0
                    self.Idepth[i] = self.sivol[i] / self.Oarea[i]
            else:
                # Melting: Water is warmer than freezing point
                Havail = 3.8e-7 * self.CwRhoW * self.vol[i] * (self.Tice - T_water)
                Hneed = self.sivol[i] * self.rhoILice
                if Hneed > abs(Havail) * self.dt:
                    melt_vol = (abs(Havail) * self.dt) / self.rhoILice
                    self.sivol[i] = max(self.sivol[i] - melt_vol, 0.0)
                else:
                    self.sivol[i] = 0.0

                if (self.sivol[i] / self.Di[i]) < self.Oarea[i]:
                    self.fsi[i] = self.sivol[i] / (self.Di[i] * self.Oarea[i]) if self.Oarea[i] > 0 else 0.0
                    self.Idepth[i] = self.Di[i]
                else:
                    self.fsi[i] = 1.0
                    self.Idepth[i] = self.sivol[i] / self.Oarea[i]

    def step_ocean(self):
        """
        Step the 8-box ocean forward using density-driven overturning, advection,
        horizontal/vertical diffusion, and surface heat/freshwater fluxes.
        """
        self.update_density()
        rho = self.rho

        # Meridional hydraulic transports (m^3/s):
        q21 = self.c21 * (self.D[0] * (rho[0] - rho[1]) + self.D[4] * (rho[4] - rho[5]))
        q23 = self.c23 * (self.D[1] * (rho[2] - rho[1]) + self.D[5] * (rho[6] - rho[5]))
        q34 = self.c34 * (self.D[3] * (rho[3] - rho[2]) + self.D[7] * (rho[7] - rho[6]))

        q56 = q21
        q76 = q23
        q87 = q34
        q15 = q21
        q48 = q34
        q62 = q56 + q76
        q73 = q87 - q76

        self.q_transports = {
            "q21": q21, "q23": q23, "q34": q34,
            "q62": q62, "q73": q73
        }

        # Upwind advective fluxes for Salinity and Temperature:
        def upwind(q_val, val_up, val_down):
            return q_val * val_up if q_val >= 0.0 else q_val * val_down

        sf21 = upwind(q21, self.S_ocean[1], self.S_ocean[0])
        sf15 = upwind(q15, self.S_ocean[0], self.S_ocean[4])
        sf56 = upwind(q56, self.S_ocean[4], self.S_ocean[5])
        tf21 = upwind(q21, self.T_ocean[1], self.T_ocean[0])
        tf15 = upwind(q15, self.T_ocean[0], self.T_ocean[4])
        tf56 = upwind(q56, self.T_ocean[4], self.T_ocean[5])

        sf34 = upwind(q34, self.S_ocean[2], self.S_ocean[3])
        sf48 = upwind(q48, self.S_ocean[3], self.S_ocean[7])
        sf87 = upwind(q87, self.S_ocean[7], self.S_ocean[6])
        tf34 = upwind(q34, self.T_ocean[2], self.T_ocean[3])
        tf48 = upwind(q48, self.T_ocean[3], self.T_ocean[7])
        tf87 = upwind(q87, self.T_ocean[7], self.T_ocean[6])

        sf62 = upwind(q62, self.S_ocean[5], self.S_ocean[1])
        tf62 = upwind(q62, self.T_ocean[5], self.T_ocean[1])
        sf23 = upwind(q23, self.S_ocean[1], self.S_ocean[2])
        tf23 = upwind(q23, self.T_ocean[1], self.T_ocean[2])
        sf73 = upwind(q73, self.S_ocean[6], self.S_ocean[2])
        tf73 = upwind(q73, self.T_ocean[6], self.T_ocean[2])
        sf76 = upwind(q76, self.S_ocean[6], self.S_ocean[5])
        tf76 = upwind(q76, self.T_ocean[6], self.T_ocean[5])

        # Vertical diffusion (boosted by 300x during static instability / convection)
        sdf_v = np.zeros(4)
        tdf_v = np.zeros(4)
        for i in range(4):
            deep_i = i + 4
            convect = 300.0 if rho[i] > rho[deep_i] else 1.0
            diff_coeff = self.Kv[i] * convect * self.Oarea[i]
            sdf_v[i] = diff_coeff * (self.S_ocean[i] - self.S_ocean[deep_i])
            tdf_v[i] = diff_coeff * (self.T_ocean[i] - self.T_ocean[deep_i])

        # Horizontal diffusion (m^3/s * unit)
        sdf12 = self.Kh[0] * self.uppersection * (self.S_ocean[0] - self.S_ocean[1])
        sdf23 = self.Kh[1] * self.uppersection * (self.S_ocean[1] - self.S_ocean[2])
        sdf34 = self.Kh[2] * self.uppersection * (self.S_ocean[2] - self.S_ocean[3])
        sdf56 = self.Kh[3] * self.lowersection * (self.S_ocean[4] - self.S_ocean[5])
        sdf76 = self.Kh[4] * self.lowersection * (self.S_ocean[6] - self.S_ocean[5])
        sdf78 = self.Kh[5] * self.lowersection * (self.S_ocean[6] - self.S_ocean[7])

        tdf12 = self.Kh[0] * self.uppersection * (self.T_ocean[0] - self.T_ocean[1])
        tdf23 = self.Kh[1] * self.uppersection * (self.T_ocean[1] - self.T_ocean[2])
        tdf34 = self.Kh[2] * self.uppersection * (self.T_ocean[2] - self.T_ocean[3])
        tdf56 = self.Kh[3] * self.lowersection * (self.T_ocean[4] - self.T_ocean[5])
        tdf76 = self.Kh[4] * self.lowersection * (self.T_ocean[6] - self.T_ocean[5])
        tdf78 = self.Kh[5] * self.lowersection * (self.T_ocean[6] - self.T_ocean[7])

        # Surface boundary fluxes:
        # Virtual salinity flux: -FWAO * Sref
        sflux = -1.0 * self.FWAO * self.Sref
        # Heat flux into ocean (W / (Cw * rhoW))
        tflux = self.Hao / self.CwRhoW

        # Update Salinities (psu):
        dS = np.zeros(8)
        dS[0] = (sf21 - sf15 - sdf_v[0] - sdf12 + sflux[0]) / self.vol[0]
        dS[1] = (sf62 - sf21 - sf23 + sdf12 - sdf23 - sdf_v[1] + sflux[1]) / self.vol[1]
        dS[2] = (sf23 + sf73 - sf34 + sdf23 - sdf34 - sdf_v[2] + sflux[2]) / self.vol[2]
        dS[3] = (sf34 - sf48 + sdf34 - sdf_v[3] + sflux[3]) / self.vol[3]
        dS[4] = (sf15 - sf56 + sdf_v[0] - sdf56) / self.vol[4]
        dS[5] = (sf56 + sf76 - sf62 + sdf56 + sdf_v[1] + sdf76) / self.vol[5]
        dS[6] = (sf87 - sf73 - sf76 - sdf76 + sdf_v[2] - sdf78) / self.vol[6]
        dS[7] = (sf48 - sf87 + sdf78 + sdf_v[3]) / self.vol[7]

        # Update Temperatures (K):
        dT = np.zeros(8)
        dT[0] = (tf21 - tf15 - tdf_v[0] - tdf12 + tflux[0]) / self.vol[0]
        dT[1] = (tf62 - tf21 - tf23 + tdf12 - tdf23 - tdf_v[1] + tflux[1]) / self.vol[1]
        dT[2] = (tf23 + tf73 - tf34 + tdf23 - tdf34 - tdf_v[2] + tflux[2]) / self.vol[2]
        dT[3] = (tf34 - tf48 + tdf34 - tdf_v[3] + tflux[3]) / self.vol[3]
        dT[4] = (tf15 - tf56 + tdf_v[0] - tdf56) / self.vol[4]
        dT[5] = (tf56 + tf76 - tf62 + tdf56 + tdf_v[1] + tdf76) / self.vol[5]
        dT[6] = (tf87 - tf73 - tf76 - tdf76 + tdf_v[2] - tdf78) / self.vol[6]
        dT[7] = (tf48 - tf87 + tdf78 + tdf_v[3]) / self.vol[7]

        self.S_ocean += dS * self.dt
        self.T_ocean += dT * self.dt

    def step(self):
        """Perform one coupled integration step (Atmosphere -> Sea Ice -> Ocean)."""
        self.step_atmosphere()
        self.step_seaice()
        self.step_ocean()

    def run(self, years: float = 100.0, output_interval_days: float = 10.0) -> Dict[str, np.ndarray]:
        """
        Integrate the coupled 8-box model over a period of years.
        """
        total_seconds = years * 365.25 * 86400.0
        n_steps = int(total_seconds / self.dt)
        save_every = max(1, int((output_interval_days * 86400.0) / self.dt))

        times_yr = []
        Teta_hist = []
        T_ocean_hist = []
        S_ocean_hist = []
        fsi_hist = []
        q34_hist = []
        q21_hist = []
        FWAO_hist = []

        for step_i in range(n_steps):
            self.step()

            if step_i % save_every == 0 or step_i == n_steps - 1:
                t_yr = (step_i * self.dt) / (365.25 * 86400.0)
                times_yr.append(t_yr)
                Teta_hist.append(self.Teta.copy())
                T_ocean_hist.append(self.T_ocean.copy())
                S_ocean_hist.append(self.S_ocean.copy())
                fsi_hist.append(self.fsi.copy())
                q34_hist.append(self.q_transports.get("q34", 0.0) * 1e-6) # Sv
                q21_hist.append(self.q_transports.get("q21", 0.0) * 1e-6) # Sv
                FWAO_hist.append(self.FWAO.copy() * 1e-6)                # Sv

        return {
            "time_yr": np.array(times_yr),
            "Teta": np.array(Teta_hist),
            "T_ocean": np.array(T_ocean_hist),
            "S_ocean": np.array(S_ocean_hist),
            "fsi": np.array(fsi_hist),
            "q_NP_Sv": np.array(q34_hist),
            "q_SP_Sv": np.array(q21_hist),
            "FWAO_Sv": np.array(FWAO_hist)
        }
