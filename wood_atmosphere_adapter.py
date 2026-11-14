"""
Coupling Adapter: Connecting the Gildor-Tziperman Atmospheric Component to Wood's AMOC Model.

This module shows how the atmospheric energy balance and Clausius-Clapeyron moisture
transport equations from Gildor & Tziperman (2001) can replace the fixed freshwater
fluxes and simplified temperature parameterization in Wood et al. (2019).
"""

from typing import Tuple
import sys
import os
import numpy as np

# Add WoodAMOC to path
wood_path = os.path.expanduser("~/Rannsoknir/WoodAMOC")
if wood_path not in sys.path:
    sys.path.append(wood_path)

from wood_box_model import WoodBoxModel, CALIBRATIONS
from unesco import unesco_full_density


class AtmosphericMoistureComponent:
    """
    Two-box or four-box atmospheric moisture transport model based on
    Gildor & Tziperman (2001).
    
    Computes:
        1. Specific humidity via Clausius-Clapeyron:
           q_sat(T) = A * exp(-B / T)
        2. Poleward moisture transport:
           F_W = K_Wa * (T_T - T_N) * q_spec / L
        3. Dynamic freshwater flux into North Atlantic:
           F_N(T) = F_W
           F_T(T) = -F_W (evaporation from thermocline)
    """

    def __init__(self, KWa: float = 5.4e13, RH: float = 0.70):
        self.A_cc = 2.53e11  # Pa
        self.B_cc = 5.42e3   # K
        self.Ps = 1.0e5      # Pa
        self.epsilon = 0.622
        self.RH = RH
        self.KWa = KWa
        self.L_dist = 1.0e7  # Distance between subtropical and polar centers (m)

    def specific_humidity(self, T_kelvin: float) -> float:
        """Specific humidity (kg/kg) from Clausius-Clapeyron."""
        qsat = self.A_cc * np.exp(-self.B_cc / T_kelvin)
        return self.RH * self.epsilon * qsat / self.Ps

    def moisture_flux_Sv(self, T_tropics_K: float, T_polar_K: float) -> float:
        """
        Poleward atmospheric moisture transport in Sverdrups (10^6 m^3/s).
        Increases with warming due to the exponential Clausius-Clapeyron relation.
        """
        q_v = self.specific_humidity(T_tropics_K)
        dT = max(T_tropics_K - T_polar_K, 0.0)
        # Transport in m^3/s
        F_m3s = self.KWa * dT * q_v / (0.25 * self.L_dist)
        # Convert to Sv (10^6 m^3/s)
        return F_m3s * 1e-6


class CoupledWoodAtmosphereModel:
    """
    Experimental prototype: Wood et al. (2019) 5-box ocean model coupled to
    dynamic atmospheric moisture transport from Gildor & Tziperman (2001).
    """

    def __init__(self, base_calibration: str = "FAMOUSA_1xCO2"):
        self.wood_params = CALIBRATIONS[base_calibration]
        self.wood_model = WoodBoxModel(self.wood_params)
        self.atmos = AtmosphericMoistureComponent()

        # Calibration scaling factor to match unperturbed FAMOUS FN ~ 0.375 Sv
        T_T_K = 273.15 + 20.0  # ~20 °C in Atlantic thermocline
        T_N_K = 273.15 + 4.5   # ~4.5 °C in North Atlantic
        raw_flux = self.atmos.moisture_flux_Sv(T_T_K, T_N_K)
        self.flux_scale = self.wood_params.FN / max(raw_flux, 1e-4)

    def dynamic_fluxes(self, TN_celsius: float, TT_celsius: float = 20.0) -> Tuple[float, float]:
        """
        Calculate dynamic FN and FT based on current temperatures.
        """
        T_T_K = 273.15 + TT_celsius
        T_N_K = 273.15 + max(TN_celsius, -1.8)
        fw_sv = self.atmos.moisture_flux_Sv(T_T_K, T_N_K) * self.flux_scale
        
        # FN is positive (net precipitation into N box)
        FN_dyn = fw_sv
        # FT is negative (evaporative loss from T box)
        FT_dyn = -fw_sv - 0.35  # additional subtropical export
        return FN_dyn, FT_dyn


if __name__ == "__main__":
    coupled = CoupledWoodAtmosphereModel()
    print("Coupled Wood-Atmosphere Moisture Component initialized successfully!")
    for t_n in [8.0, 4.0, 0.0, -1.5]:
        fn, ft = coupled.dynamic_fluxes(TN_celsius=t_n)
        print(f"When TN = {t_n:+4.1f} °C -> Dynamic North Atlantic Freshwater Flux FN = {fn:.4f} Sv (FT = {ft:.4f} Sv)")
