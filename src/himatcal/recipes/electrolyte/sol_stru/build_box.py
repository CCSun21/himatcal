from __future__ import annotations

from pathlib import Path

import mbuild as mb


class ElectrolyteBuilder:
    def __init__(self, mol_path, mol_name_list: list | None = None):
        self.mol_path = Path(mol_path)
        self.EC = mb.load(str(self.mol_path / "EC.pdb"))
        self.DEC = mb.load(str(self.mol_path / "DEC.pdb"))
        self.EMC = mb.load(str(self.mol_path / "EMC.pdb"))
        self.DMC = mb.load(str(self.mol_path / "DMC.pdb"))
        self.Li = mb.load(str(self.mol_path / "Li.pdb"))
        self.PF6 = mb.load(str(self.mol_path / "PF6.pdb"))
        if mol_name_list:
            for mol_name in mol_name_list:
                setattr(self, mol_name, mb.load(str(self.mol_path / f"{mol_name}.pdb")))

    def fix_pdb(self, pdb_path, round_componds_dict, pdb_save_path, capital=True):
        # The same code as before...
        """
        将build_electrolytes生成的电解液组成的pdb文件中的RESIDUE NAME和RESIDUE INDEX修改为正确的值
        """
        mb_compound_list = list(round_componds_dict.keys())
        if "LiPF6(M)" in mb_compound_list:
            self.extracted_LiPF6(round_componds_dict, "LiPF6(M)")
        if "LiPF6" in mb_compound_list:
            self.extracted_LiPF6(round_componds_dict, "LiPF6")
        mb_compound_list = list(round_componds_dict.keys())
        # print(mb_compound_list)
        line_start = 1
        line_stop = 1
        with Path.open(Path(pdb_path), "r") as f:
            lines = f.readlines()
        for compound_index in range(len(mb_compound_list)):
            if compound_index == 0:
                line_start = 1
                line_stop = (
                    1 + round_componds_dict[mb_compound_list[compound_index]]["n_atoms"]
                )
            else:
                line_start = line_stop
                line_stop = (
                    line_start
                    + round_componds_dict[mb_compound_list[compound_index]]["n_atoms"]
                )
            # print(line_start, line_stop)
            compound = mb_compound_list[compound_index]
            resdiue_index = 1
            for j in range(line_start, line_stop):
                n_compounds = round_componds_dict[compound]["mol"]
                atom_per_compound = (
                    round_componds_dict[compound]["n_atoms"] / n_compounds
                )
                resdiue_index = int((j - line_start) / atom_per_compound) + 1
                if capital:
                    resdiue_name = mb_compound_list[compound_index]
                else:
                    resdiue_name = mb_compound_list[compound_index].lower()
                # resdiue_name = mb_compound_list[compound_index]
                lines[j] = (
                    (
                        (
                            lines[j][:17]
                            + str(" " * (3 - len(resdiue_name)) + str(resdiue_name))
                        )
                        + lines[j][20:22]
                    )
                    + str(" " * (4 - len(str(resdiue_index))) + str(resdiue_index))
                    + lines[j][26:]
                )
                # lines[j] = lines[j].replace('HETATM','ATOM  ')
        with Path.open(Path(pdb_save_path), "w") as f:
            f.writelines(lines)

    def extracted_LiPF6(self, round_componds_dict, arg1):
        MOL_LPF6 = round_componds_dict[arg1]["mol"]
        del round_componds_dict[arg1]
        round_componds_dict["PF6"] = {"mol": MOL_LPF6, "n_atoms": MOL_LPF6 * 7}
        round_componds_dict["Li"] = {"mol": MOL_LPF6, "n_atoms": MOL_LPF6 * 1}

    def build_box(
        self, box_electrolyte_composition, density, box, save_path, capital=True
    ):
        """
        使用mbuild填充盒子并修正残基名和残基序号
        """
        box_compound = [getattr(self, key) for key in box_electrolyte_composition]
        box_n_compounds = [
            box_electrolyte_composition[key]["mol"]
            for key in box_electrolyte_composition
        ]
        if density:
            box_electrolyte = mb.packing.fill_box(
                compound=box_compound,
                n_compounds=box_n_compounds,
                density=density * 1000,  # kg/m^3
                # ratio=[1,1,1],
                # box=[3,3,3]
            )
        else:
            box_electrolyte = mb.packing.fill_box(
                compound=box_compound,
                n_compounds=box_n_compounds,
                density=density *1000,# kg/m^3
                # ratio=[1,1,1],
                box=box,
            )

        box_electrolyte.save(save_path, overwrite=True)
        self.fix_pdb(save_path, box_electrolyte_composition, save_path, capital=capital)
