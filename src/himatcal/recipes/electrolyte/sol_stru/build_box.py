from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path


class ElectrolyteBuilder:
    def __init__(self, mol_path, mol_name_list: list | None = None):
        """初始化电解液构建器，加载所有需要的分子结构文件路径

        Args:
            mol_path: 分子结构文件所在的目录路径
            mol_name_list: 额外需要加载的分子名称列表
        """
        self.mol_path = Path(mol_path)
        # 存储分子结构文件的路径而不是直接加载
        self.molecule_paths = {
            "EC": str(self.mol_path / "EC.pdb"),
            "DEC": str(self.mol_path / "DEC.pdb"),
            "EMC": str(self.mol_path / "EMC.pdb"),
            "DMC": str(self.mol_path / "DMC.pdb"),
            "Li": str(self.mol_path / "Li.pdb"),
            "PF6": str(self.mol_path / "PF6.pdb"),
        }

        # 加载额外的分子结构文件
        if mol_name_list:
            for mol_name in mol_name_list:
                self.molecule_paths[mol_name] = str(self.mol_path / f"{mol_name}.pdb")

        # 存储分子的原子数量信息, 用于后续处理
        self.molecule_atom_counts = self._get_atom_counts()

    def _get_atom_counts(self) -> dict:
        """计算每个分子的原子数量

        Returns:
            包含每个分子的原子数量的字典
        """
        atom_counts = {}
        for mol_name, mol_path in self.molecule_paths.items():
            # 计算每个PDB文件中的原子数量
            with Path(mol_path).open() as f:
                lines = f.readlines()
                # 计算ATOM或HETATM行的数量
                count = sum(1 for line in lines if line.startswith(("ATOM", "HETATM")))
                atom_counts[mol_name] = count
        return atom_counts

    def fix_pdb(self, pdb_path, round_componds_dict, pdb_save_path, capital=True):
        """
        将生成的电解液组成的pdb文件中的RESIDUE NAME和RESIDUE INDEX修改为正确的值

        Args:
            pdb_path: 输入的PDB文件路径
            round_componds_dict: 包含分子信息的字典
            pdb_save_path: 保存修正后PDB文件的路径
            capital: 是否使用大写的残基名称
        """
        mb_compound_list = list(round_componds_dict.keys())
        if "LiPF6(M)" in mb_compound_list:
            self.extracted_LiPF6(round_componds_dict, "LiPF6(M)")
        if "LiPF6" in mb_compound_list:
            self.extracted_LiPF6(round_componds_dict, "LiPF6")
        mb_compound_list = list(round_componds_dict.keys())

        line_start = 1
        line_stop = 1
        with Path(pdb_path).open("r") as f:
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
        with Path(pdb_save_path).open("w") as f:
            f.writelines(lines)

    def extracted_LiPF6(self, round_componds_dict, arg1):
        """
        分离LiPF6到Li和PF6两个组分

        Args:
            round_componds_dict: 包含分子信息的字典
            arg1: LiPF6的键名
        """
        MOL_LPF6 = round_componds_dict[arg1]["mol"]
        del round_componds_dict[arg1]
        round_componds_dict["PF6"] = {"mol": MOL_LPF6, "n_atoms": MOL_LPF6 * 7}
        round_componds_dict["Li"] = {"mol": MOL_LPF6, "n_atoms": MOL_LPF6 * 1}

    def build_box(
        self,
        box_electrolyte_composition,
        density,
        box=None,
        save_path=None,
        capital=True,
    ):
        """
        使用packmol填充盒子并修正残基名和残基序号

        Args:
            box_electrolyte_composition: 电解液组成字典
            density: 目标密度(g/cm^3), 如果为None则使用提供的box尺寸
            box: 盒子尺寸[x,y,z](单位: Å)
            save_path: 保存结果的路径
            capital: 是否使用大写的残基名称

        Returns:
            保存的PDB文件路径
        """
        # 创建临时目录存储中间文件
        with tempfile.TemporaryDirectory() as temp_dir:
            # 准备packmol输入文件
            packmol_input_path = Path(temp_dir) / "packmol_input.inp"
            output_path = save_path if save_path else Path(temp_dir) / "output.pdb"

            # 计算或使用提供的盒子尺寸
            if density and not box:
                # 计算总质量
                total_mass = 0.0
                # 这里需要有分子量数据，这是一个示例实现
                # 实际应用中需要提供准确的分子量数据
                molecule_masses = {
                    "EC": 88.06,  # 示例分子量，需要替换为实际值
                    "DEC": 118.13,
                    "EMC": 104.10,
                    "DMC": 90.08,
                    "Li": 6.94,
                    "PF6": 144.96,
                    # 添加其他分子的分子量
                }

                for mol_name, info in box_electrolyte_composition.items():
                    if mol_name in molecule_masses:
                        total_mass += info["mol"] * molecule_masses[mol_name]

                # 计算体积 (cm^3)
                volume = total_mass / density
                # 转换为 Å^3
                volume_ang = volume * 1e24
                # 假设立方体盒子，计算边长
                box_length = volume_ang ** (1 / 3)
                box = [box_length, box_length, box_length]

            # 确保box是有效的
            if not box or len(box) != 3:
                raise ValueError("必须提供有效的盒子尺寸或密度")

            # 生成packmol输入文件
            with packmol_input_path.open("w") as f:
                f.write("tolerance 2.0\n")
                f.write("filetype pdb\n")
                f.write(f"output {output_path}\n\n")

                for mol_name, info in box_electrolyte_composition.items():
                    if mol_name in self.molecule_paths:
                        num_molecules = info["mol"]
                        f.write(f"structure {self.molecule_paths[mol_name]}\n")
                        f.write(f"  number {num_molecules}\n")
                        f.write(f"  inside box 0. 0. 0. {box[0]} {box[1]} {box[2]}\n")
                        f.write("end structure\n\n")

            # 运行packmol
            try:
                # 修正packmol执行方式
                with packmol_input_path.open("r") as input_file:
                    subprocess.run(
                        ["packmol"],
                        stdin=input_file,
                        check=True,
                        capture_output=True,
                        text=True,
                    )

                # 更新分子的原子数信息
                updated_composition = {}
                for mol_name, info in box_electrolyte_composition.items():
                    # 计算该类型分子的总原子数
                    if mol_name in self.molecule_atom_counts:
                        n_atoms = info["mol"] * self.molecule_atom_counts[mol_name]
                        updated_composition[mol_name] = {
                            "mol": info["mol"],
                            "n_atoms": n_atoms,
                        }

                # 修正PDB文件中的残基信息
                self.fix_pdb(
                    output_path, updated_composition, output_path, capital=capital
                )

                if save_path != output_path:
                    # 如果使用的是临时文件, 需要复制到指定路径
                    with (
                        Path(output_path).open() as src,
                        Path(save_path).open("w") as dst,
                    ):
                        dst.write(src.read())

                return save_path

            except subprocess.CalledProcessError as e:
                logger = logging.getLogger(__name__)
                logger.error(f"packmol运行出错: {e}")
                logger.error(f"标准输出: {e.stdout}")
                logger.error(f"错误输出: {e.stderr}")
                raise RuntimeError("packmol执行失败") from e
