"""
Turbo expander NGL recovery plant: a simplified steady-state simulation.

Generic, educational model built with the Peng-Robinson equation of state
(thermo library). It is not a design tool and does not represent any real plant.

Flowsheet
    Feed gas (sweet, dry) -> inlet cooling (gas/gas exchanger + chiller)
    -> cold separator
        vapour -> turbo expander -> demethanizer top stage
        liquid -> JT valve      -> demethanizer top stage
    Demethanizer (two equilibrium stages: top stage + reboiler)
        top vapour   -> residue gas -> expander-driven booster + residue compressor
        bottom liquid -> NGL product (ethane plus)
"""
from dataclasses import dataclass, field
import numpy as np
from thermo import ChemicalConstantsPackage, PRMIX, CEOSLiquid, CEOSGas, FlashVL
from thermo.interaction_parameters import IPDB

COMPONENTS = ["nitrogen", "carbon dioxide", "methane", "ethane", "propane",
              "isobutane", "butane", "isopentane", "pentane", "hexane"]
LABELS = ["N2", "CO2", "C1", "C2", "C3", "iC4", "nC4", "iC5", "nC5", "C6"]

# Typical rich natural gas after sweetening and dehydration (mol %), generic values.
DEFAULT_FEED = {"N2": 1.0, "CO2": 0.5, "C1": 82.0, "C2": 8.0, "C3": 4.5,
                "iC4": 1.0, "nC4": 1.5, "iC5": 0.5, "nC5": 0.5, "C6": 0.5}

KMOL_PER_MMSCF = 1e6 / 379.49 * 0.45359237  # standard conditions 60 F, 14.696 psia

_consts, _props = ChemicalConstantsPackage.from_IDs(COMPONENTS)
_kijs = IPDB.get_ip_asymmetric_matrix("ChemSep PR", _consts.CASs, "kij")
_eos = dict(Tcs=_consts.Tcs, Pcs=_consts.Pcs, omegas=_consts.omegas, kijs=_kijs)
FLASHER = FlashVL(_consts, _props,
                  liquid=CEOSLiquid(PRMIX, _eos, HeatCapacityGases=_props.HeatCapacityGases),
                  gas=CEOSGas(PRMIX, _eos, HeatCapacityGases=_props.HeatCapacityGases))
MW = np.array(_consts.MWs)  # g/mol


@dataclass
class Stream:
    name: str
    F: float            # kmol/h
    z: np.ndarray       # mole fractions
    T: float            # K
    P: float            # Pa
    H: float = 0.0      # J/mol
    VF: float = 0.0

    @property
    def comp_flows(self):
        return self.F * self.z


def _norm(z):
    z = np.clip(np.asarray(z, float), 0, None)
    return z / z.sum()


def flash_TP(T, P, z):
    return FLASHER.flash(T=T, P=P, zs=list(_norm(z)))


def flash_PH(P, H, z):
    return FLASHER.flash(P=P, H=H, zs=list(_norm(z)))


def flash_PS(P, S, z):
    return FLASHER.flash(P=P, S=S, zs=list(_norm(z)))


def _liq(res):
    return res.liquids[0] if getattr(res, "liquids", None) else None


def split(res, F, name):
    """Return (vapour, liquid) streams from a flash result."""
    VF = res.VF if res.VF is not None else (1.0 if res.gas else 0.0)
    v = l = None
    if res.gas is not None and VF > 1e-9:
        v = Stream(name + " vapour", F * VF, np.array(res.gas.zs), res.T, res.P, res.gas.H(), 1.0)
    if _liq(res) is not None and VF < 1 - 1e-9:
        l = Stream(name + " liquid", F * (1 - VF), np.array(_liq(res).zs), res.T, res.P, _liq(res).H(), 0.0)
    return v, l


def compress(P_in, T_in, P_out, z, F, eff):
    """Adiabatic compression with isentropic efficiency. Returns (power kW, outlet T K)."""
    r_in = flash_TP(T_in, P_in, z)
    r_s = flash_PS(P_out, r_in.S(), z)
    dH = (r_s.H() - r_in.H()) / eff
    r_out = flash_PH(P_out, r_in.H() + dH, z)
    return F * dH / 3600.0, r_out.T  # kmol/h * J/mol -> kW


@dataclass
class Inputs:
    feed_mmscfd: float = 100.0
    feed_mol_pct: dict = field(default_factory=lambda: dict(DEFAULT_FEED))
    P_inlet_bar: float = 65.0
    T_inlet_C: float = 30.0
    T_cold_C: float = -35.0
    P_column_bar: float = 22.0
    expander_eff: float = 0.82
    compressor_eff: float = 0.78
    C1_to_C2_bottoms: float = 0.02     # NGL product methane spec (molar ratio)
    P_sales_bar: float = 65.0
    T_residue_warm_C: float = 25.0      # residue gas after gas/gas exchanger


def simulate(inp: Inputs):
    z0 = _norm([inp.feed_mol_pct[k] for k in LABELS])
    F0 = inp.feed_mmscfd * KMOL_PER_MMSCF / 24.0            # kmol/h
    P1, T1 = inp.P_inlet_bar * 1e5, inp.T_inlet_C + 273.15
    Pc = inp.P_column_bar * 1e5

    feed_res = flash_TP(T1, P1, z0)

    # 1. Inlet cooling to the cold separator (0.5 bar pressure drop)
    P2, T2 = P1 - 0.5e5, inp.T_cold_C + 273.15
    cold = flash_TP(T2, P2, z0)
    Q_cooling_kW = F0 * (cold.H() - feed_res.H()) / 3600.0
    sep_v, sep_l = split(cold, F0, "Cold separator")

    # 2. Turbo expander on the separator vapour
    exp_power_kW = 0.0
    exp_out = None
    if sep_v is not None:
        r_in = flash_TP(sep_v.T, sep_v.P, sep_v.z)
        r_s = flash_PS(Pc, r_in.S(), sep_v.z)
        H_out = r_in.H() - inp.expander_eff * (r_in.H() - r_s.H())
        r_out = flash_PH(Pc, H_out, sep_v.z)
        exp_power_kW = sep_v.F * (r_in.H() - H_out) / 3600.0
        exp_out = Stream("Expander outlet", sep_v.F, sep_v.z, r_out.T, Pc, H_out,
                         r_out.VF if r_out.VF is not None else 1.0)

    # 3. JT valve on the separator liquid
    jt_out = None
    if sep_l is not None:
        r_jt = flash_PH(Pc, sep_l.H, sep_l.z)
        jt_out = Stream("JT valve outlet", sep_l.F, sep_l.z, r_jt.T, Pc, sep_l.H,
                        r_jt.VF if r_jt.VF is not None else 0.0)

    # 4. Demethanizer, simplified as a cold top stage plus an ideal stripping section.
    #    Expander outlet and JT outlet mix at column pressure (adiabatic). Vapour from this
    #    cold top stage leaves as residue gas. The cold liquid is stripped: methane, nitrogen
    #    and CO2 go overhead until the bottoms meet the C1/C2 spec; ethane plus stays in the NGL.
    feeds = [s for s in (exp_out, jt_out) if s is not None]
    F_feed = sum(s.F for s in feeds)
    n_feed = sum(s.F * s.z for s in feeds)
    H_mix = sum(s.F * s.H for s in feeds) / F_feed
    top = flash_PH(Pc, H_mix, n_feed / F_feed)
    top_v, top_l = split(top, F_feed, "Demethanizer top")
    iN2, iCO2, i1, i2 = (LABELS.index(k) for k in ("N2", "CO2", "C1", "C2"))
    n_liq = top_l.comp_flows.copy() if top_l else np.zeros(len(LABELS))
    n_bot = n_liq.copy()
    n_bot[iN2] = 0.0
    n_bot[iCO2] = 0.0
    n_bot[i1] = min(n_liq[i1], inp.C1_to_C2_bottoms * n_liq[i2])
    n_strip = n_liq - n_bot
    T_top = top.T
    # Reboiler duty from an energy balance: liquid in -> bottoms at bubble point + stripped vapour at top T
    bub = FLASHER.flash(P=Pc, VF=0, zs=list(_norm(n_bot)))
    H_bot = n_bot.sum() * bub.H()
    if n_strip.sum() > 1e-9:
        H_strip = n_strip.sum() * flash_TP(T_top, Pc, n_strip).H()
    else:
        H_strip = 0.0
    H_liq_in = top_l.F * top_l.H if top_l else 0.0
    Q_reb_kW = (H_bot + H_strip - H_liq_in) / 3600.0
    T_reb = bub.T
    n_res = (top_v.comp_flows if top_v else 0) + n_strip
    residue = Stream("Residue gas", n_res.sum(), _norm(n_res), T_top, Pc)
    ngl = Stream("NGL product", n_bot.sum(), _norm(n_bot), T_reb, Pc)

    # 5. Residue gas recompression back to sales pressure
    T_warm = inp.T_residue_warm_C + 273.15
    comp_power_kW, T_disch = compress(Pc, T_warm, inp.P_sales_bar * 1e5, residue.z, residue.F, inp.compressor_eff)
    net_comp_kW = comp_power_kW - exp_power_kW   # expander drives the booster compressor

    feed_flows = F0 * z0
    recovery = {LABELS[i]: (ngl.comp_flows[i] / feed_flows[i] * 100 if feed_flows[i] > 0 else 0.0)
                for i in range(len(LABELS))}
    ngl_mass_tpd = (ngl.comp_flows * MW).sum() / 1000 * 24       # tonnes/day
    residue_mmscfd = residue.F * 24 / KMOL_PER_MMSCF

    return {
        "inputs": inp,
        "feed_kmolh": F0,
        "cold_sep_vapour_fraction": cold.VF,
        "T_expander_out_C": exp_out.T - 273.15 if exp_out else None,
        "T_jt_out_C": jt_out.T - 273.15 if jt_out else None,
        "T_reboiler_C": T_reb - 273.15,
        "T_column_top_C": T_top - 273.15,
        "expander_power_kW": exp_power_kW,
        "recompression_power_kW": comp_power_kW,
        "net_compression_kW": net_comp_kW,
        "inlet_cooling_kW": -Q_cooling_kW,
        "reboiler_duty_kW": Q_reb_kW,
        "recovery_pct": recovery,
        "ngl_tpd": ngl_mass_tpd,
        "ngl_kmolh": ngl.F,
        "ngl_comp": dict(zip(LABELS, ngl.z * 100)),
        "residue_mmscfd": residue_mmscfd,
        "residue_comp": dict(zip(LABELS, residue.z * 100)),
        "residue_discharge_T_C": T_disch - 273.15,
    }


if __name__ == "__main__":
    import time
    t = time.time()
    r = simulate(Inputs())
    print(f"solved in {time.time()-t:.1f}s")
    for k, v in r.items():
        if k in ("inputs",):
            continue
        if isinstance(v, dict):
            print(k, {a: round(b, 2) for a, b in v.items()})
        else:
            print(k, round(v, 2) if isinstance(v, float) else v)
