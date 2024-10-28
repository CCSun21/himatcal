from __future__ import annotations

from ase import units


def get_calc(label, kwargs: dict | None = None):
    """ """
    kwargs = kwargs or {}
    device = kwargs.get("device", "cuda" if label.startswith("orb") else "cpu")

    if label.startswith("orb"):  # * orb_v2, orb_d3_v2, orb_d3_xm_v2, orb_d3_sm_v2
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator

        model_path = kwargs.get("model")
        orbff = (
            getattr(pretrained, label)(device=device, weight_path=model_path)
            if model_path
            else getattr(pretrained, label)(device=device)
        )
        return ORBCalculator(orbff, device=device)

    if label == "mace_mp":
        from mace.calculators import mace_mp

        return mace_mp(
            model=kwargs.get("model", "medium"),  # * small , medium, large
            device=device,  # * cpu, cuda
            default_dtype=kwargs.get("default_type", "float32"),
            dispersion=kwargs.get("dispersion", True),
            damping=kwargs.get("damping", "bj"),
            dispersion_xc=kwargs.get("dispersion_xc", "pbe"),
            dispersion_cutoff=kwargs.get("dispersion_cutoff", 40 * units.Bohr),
        )

    if label == "mace_off":
        from mace.calculators import mace_off

        return mace_off(
            model=kwargs.get("model", "medium"),
            device=device,
            default_dtype=kwargs.get("default_dtype", "float64"),
        )

    if label == "mace_anicc":
        from mace.calculators import mace_anicc

        return mace_anicc(model_path=kwargs.get("model"), device=device)

    if label == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator

        return CHGNetCalculator(model_path=kwargs.get("model"), use_device=device)
    if label == "aimnet2":
        from himatcal.calculator.aimnet import AIMNet2ASE

        return AIMNet2ASE("aimnet2_b973c", charge=kwargs.get("charge", 0), mult=kwargs.get("mult", 1))
    return None
