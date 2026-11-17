# Copied from the paper's analysis code.
"""Plain names of the fits' parameters, for figures.

Every parameter is named by the name it has in its own experiment's fit; a
combination's own parameters, ``name@experiment``, carry their experiment's
name in front, and a dataset's own copy, ``name[datasetN]``, its set's number.
"""

from __future__ import annotations

#: Each experiment's name, as a combination's parameters are labelled with it.
EXPERIMENTS = {"borexino": "Borexino", "kamland": "KamLAND", "snoplus": "SNO+", "juno": "JUNO"}
#: Names that start with a proper noun, which keep their capital after an experiment's name.
PROPER = ("Birks",)
#: Plain names of the parameters, by the name each has in its own experiment's fit.
NAMES = {"norm_crust_local_U238": "Local crust, U", "norm_crust_local_Th232": "Local crust, Th",
         "norm_crust_roc_U238": "Far crust, U", "norm_crust_roc_Th232": "Far crust, Th",
         "norm_mantle_U238": "Mantle, U", "norm_mantle_Th232": "Mantle, Th", "norm_mantle": "Mantle",
         "norm_crust_local": "Local crust", "norm_crust_roc": "Far crust",
         "norm_lithospheric_mantle": "Lithospheric mantle",
         "norm_reactor": "Reactor flux", "norm_world_reactors": "Reactors beyond 450 km",
         "norm_geo": "Geoneutrino rate",
         "osc_sin2_th12": r"$\sin^2\theta_{12}$", "osc_dm2_21": r"$\Delta m^2_{21}$",
         "sin2_th12": r"$\sin^2\theta_{12}$", "sin2_th13": r"$\sin^2\theta_{13}$",
         "dm2_21_ev2": r"$\Delta m^2_{21}$", "dm2_31_ev2": r"$\vert\Delta m^2_{31}\vert$",
         "osc_sin2_2th13": r"$\sin^2 2\theta_{13}$", "osc_dm2_32": r"$\Delta m^2_{32}$",
         "norm_li9": r"$^9$Li", "norm_li9_he8": r"$^9$Li and $^8$He", "shape_li9_he8": r"$^9$Li and $^8$He, shape",
         "norm_accidental": "Accidental coincidences", "norm_alpha_n": r"$(\alpha,n)$",
         "norm_alpha_n_gs": r"$(\alpha,n)$, ground state", "norm_alpha_n_es": r"$(\alpha,n)$, excited states",
         "norm_alpha_p": r"$\alpha$-$p$ coincidences", "norm_atmospheric": "Atmospheric neutrinos",
         "norm_bipo": r"$^{214}$Bi-$^{214}$Po", "norm_other": "Other backgrounds", "delta_e_scale": "Energy scale",
         "delta_detector_rate": "Detection rate", "delta_kb": "Birks constant", "smear_var": "Extra smearing",
         "delta_class_s": "IBD classifier efficiency", "proton_scale": "Proton energy scale",
         "class_a2": r"$(\alpha,n)$ classifier factor", "m_lightest_ev": "Lightest neutrino mass",
         "log10_tau1_over_m1": r"$\log_{10}[(\tau_1/m_1)/(\mathrm{s~eV}^{-1})]$"}


def label(name: str, combined: bool = False) -> str:
    """A parameter's plain name.

    Args:
        name: As the fit names it: ``name[datasetN]@experiment`` at most.
        combined: Whether the parameter is shown for a combination, whose own
            parameters are then labelled with their experiment.

    Returns:
        E.g. ``Mantle, U``, ``Proton energy scale (set 2)`` or
        ``KamLAND: detection rate``.

    Raises:
        KeyError: If the parameter has no plain name yet.
    """
    base, _sep, tag = name.partition("@")
    base, _bracket, dataset = base.partition("[")
    if base.startswith("reactor_shape_"):
        text = f"Reactor spectrum, component {int(base[len('reactor_shape_'):]) + 1}"
    else:
        text = NAMES[base]
    if dataset:
        text += f" (set {dataset.rstrip(']')[-1]})"
    if combined and tag:
        # The name follows the experiment's, so it starts lower case, unless it starts with a proper noun.
        proper = text.startswith(PROPER) or text[1:2].isupper()
        lower = text if proper or not text[0].isupper() else text[0].lower() + text[1:]
        text = f"{EXPERIMENTS[tag]}: {lower}"
    return text
