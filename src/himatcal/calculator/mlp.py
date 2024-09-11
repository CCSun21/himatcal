from __future__ import annotations


def get_calc(label, kwards=None):
    if label in["orb_d3_v1","orb_v1","orb_d3_sm_v1","orb_d3_xs_v1"]:
        # * install the orb_models package from https://github.com/orbital-materials/orb-models/
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator
        device = kwards.get("device", "cuda") if kwards is not None else "cuda"
        orbff = getattr(pretrained,label)(
            device=device
        )
        return ORBCalculator(orbff, device=device)
    else:
        return None
