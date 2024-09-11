from __future__ import annotations

import sys

from mp_api.client import MPRester

from himatcal import SETTINGS

mpr = MPRester(api_key=SETTINGS.MAPI_KEY)


def get_mp_formula(formula):
    structure = mpr.materials.summary.search(
        formula=formula,
        fields=[
            "material_id",
            "formula_pretty",
            "nsites",
            "symmetry",
            "energy_above_hull",
            "band_gap",
        ],
    )
    return structure


def get_mp_id(material_id):
    structure = mpr.materials.summary.search(
        material_ids=[material_id],
        fields=[
            "material_id",
            "formula_pretty",
            "nsites",
            "symmetry",
            "energy_above_hull",
            "band_gap",
        ],
    )
    return structure


def get_mp_chemsys(chemsys):
    # chemsys = 'Li-Fe-O'
    elements = chemsys.split("-")
    structure = mpr.materials.summary.search(
        elements=elements,
        fields=[
            "material_id",
            "formula_pretty",
            "nsites",
            "symmetry",
            "energy_above_hull",
            "band_gap",
        ],
    )
    return structure


def extract_mp_prop(MPdoc):
    data = {
        "MPID": [],
        "formula": [],
        "CrystalSystem": [],
        "SpaceGroup": [],
        "nsites": [],
        "EHull": [],
        "BandGap": [],
    }
    for i in range(len(MPdoc)):
        data["MPID"].append("mp-" + str(MPdoc[i].material_id.split("-")[1]))
        data["formula"].append(MPdoc[i].formula_pretty)
        data["CrystalSystem"].append(str(MPdoc[i].symmetry.crystal_system))
        data["SpaceGroup"].append(str(MPdoc[i].symmetry.symbol))
        data["nsites"].append(MPdoc[i].nsites)
        data["EHull"].append(MPdoc[i].energy_above_hull)
        data["BandGap"].append(MPdoc[i].band_gap)
    return data


def extract_mp_data(material_id, path):
    stru = mpr.materials.summary.search(material_ids=[material_id])
    stru = stru[0].structure.to(fmt="poscar", filename=f"{path}/{material_id}.vasp")


input_type = sys.argv[1]
input_content = sys.argv[2]
# if argv[3] is not None: path = sys.argv[3]
if len(sys.argv) > 3 and sys.argv[3] is not None:
    path = sys.argv[3]

if input_type == "chemsys":
    data = get_mp_chemsys(input_content)
    extract_data = extract_mp_prop(data)
    print(extract_data)
elif input_type == "extract":
    data = extract_mp_data(input_content, path)
elif input_type == "formula":
    data = get_mp_formula(input_content)
    extract_data = extract_mp_prop(data)
    print(extract_data)
elif input_type == "id":
    data = get_mp_id(input_content)
    extract_data = extract_mp_prop(data)
    print(extract_data)
