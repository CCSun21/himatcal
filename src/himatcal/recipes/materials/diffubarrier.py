from __future__ import annotations

import logging
from pathlib import Path

import matgl
from ase import Atoms
from ase.io import read, write
from ase.mep import NEB, NEBTools
from ase.optimize import FIRE
from matgl.ext.ase import M3GNetCalculator, Relaxer
from pathos.multiprocessing import ProcessingPool as Pool
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser

from himatcal.recipes.materials.core import ase_to_pymatgen, pymatgen_to_ase

Logger = logging.getLogger(__name__)


class OgStructure:
    """
    model from OgStructure, using to generate the NEB Path in the Crystal.
    """

    def __init__(
        self,
        structure: Atoms | Structure | str = None,
    ) -> None:
        if structure is None:
            return

        if isinstance(structure, str):
            parser = CifParser.from_str(structure)
            structure = parser.get_structures()[0]
        elif isinstance(structure, Atoms):
            self.structure = ase_to_pymatgen(structure)
        elif isinstance(structure, Structure):
            self.structure = structure

        self.structure.sort()
        self.pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")

    def center(self, about_atom=None, about_point=None):
        if about_atom is not None:
            about_point = self.structure.frac_coords[about_atom]
        elif about_point is None:
            about_point = [0, 0, 0]
        sys = pymatgen_to_ase(self.structure)
        sys.center(about=about_point)
        self.structure = ase_to_pymatgen(sys)
        return self

    def epsilon(self, a, b):
        return abs(a - b) < 1e-6

    def equivalent_sites(self, i, site):
        return bool(
            self.epsilon(self.structure.frac_coords[i][0] % 1, site.frac_coords[0] % 1)
            and self.epsilon(
                self.structure.frac_coords[i][1] % 1, site.frac_coords[1] % 1
            )
            and self.epsilon(
                self.structure.frac_coords[i][2] % 1, site.frac_coords[2] % 1
            )
        )

    def _get_site_for_neighbor_site(self, neighbor):
        for i_site in range(len(self.structure)):
            if self.equivalent_sites(i_site, neighbor):
                return i_site
        logging.error(f"og:No equivalent: {neighbor}")
        return None

    def relax(self):
        relaxer = Relaxer(potential=self.pot, relax_cell=False)
        relax_results = relaxer.relax(self.structure, verbose=True, fmax=0.1)
        self.structure = relax_results["final_structure"]
        self.total_energy = relax_results["trajectory"].energies[-1]
        return self

    def generate_neb(
        self,
        moving_atom_species,
        num_images=5,
        r=4,
    ) -> None:
        structure = self.structure
        self.neb_paths = []

        def get_moving_atom_neighbors(site):
            return [
                neighbor
                for neighbor in structure.get_neighbors(site=site, r=r)
                if neighbor.specie.symbol == moving_atom_species
            ]

        def modify_structure(new_structure, i_site, i_neighbor_site):
            ogs = OgStructure(new_structure)
            new_structure = ogs.center(i_site, [0.5, 0.5, 0.5]).structure

            initial_structure = new_structure.copy()
            final_structure = new_structure.copy()
            initial_structure.remove_sites([i_site, i_neighbor_site])
            initial_structure.append(
                species=new_structure[i_site].species,
                coords=new_structure[i_site].frac_coords,
            )
            initial = pymatgen_to_ase(initial_structure)
            final_structure.remove_sites([i_site, i_neighbor_site])
            final_structure.append(
                species=new_structure[i_neighbor_site].species,
                coords=new_structure[i_neighbor_site].frac_coords,
            )
            final = pymatgen_to_ase(final_structure)
            return initial, final

        def handle_neb_path(i_site, i_neighbor_site):
            new_structure = structure.copy()
            self.neb_paths.append([i_site, i_neighbor_site])
            neb_folder = f"neb_path_{i_neighbor_site}_{i_site}"

            initial, final = modify_structure(new_structure, i_site, i_neighbor_site)

            self.images = (
                [initial] + [initial.copy() for _ in range(num_images)] + [final]
            )
            self.neb = NEB(self.images)
            Path(neb_folder).mkdir(exist_ok=True)
            self.neb.interpolate(mic=True)
            for i, img in enumerate(self.images):
                image = OgStructure(img).relax()
                image_str = image.structure.to(fmt="poscar")
                with Path.open(Path(f"{neb_folder}/{str(i).zfill(2)}.vasp"), "w") as f:
                    f.write(image_str)

        for i_site in range(len(structure)):
            if structure[i_site].specie.symbol == moving_atom_species:
                neighbors = get_moving_atom_neighbors(structure[i_site])
                logging.info(f"og:Full neighbors: {neighbors}")
                logging.info(
                    f"og:Checking site, {i_site} Surrounded by {len(neighbors)}"
                )
                for neighbor in neighbors:
                    i_neighbor_site = self._get_site_for_neighbor_site(neighbor)
                    logging.info(f"og:Checking neighbor site {i_neighbor_site}")
                    if i_neighbor_site is None:
                        logging.error("og:Really? Wrong site in neighbor list!")
                    if [i_site, i_neighbor_site] in self.neb_paths or [
                        i_neighbor_site,
                        i_site,
                    ] in self.neb_paths:
                        continue
                    handle_neb_path(i_site, i_neighbor_site)


class DiffusionBarrierCal:
    def __init__(
        self,
        atoms: Atoms,
        label: str,
        ion_type: str = "Li",
        neb_point_number: int = 5,
        nebor_cutoff: float = 4,
        parallel: int = 4,
    ):
        self.label = label
        self.atoms = atoms
        self.ion_type = ion_type
        self.neb_point_number = neb_point_number
        self.nebor_cutoff = nebor_cutoff
        self.parallel = parallel
        self._CWD = Path.cwd()

    def relax(self, atoms, relax_cell=False):
        pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        relaxer = Relaxer(potential=pot, relax_cell=relax_cell)
        relax_results = relaxer.relax(atoms, verbose=True, fmax=0.1)
        return pymatgen_to_ase(relax_results["final_structure"])

    def generate_neb_path(self):
        """
        Generate the neb path using the OgStrucure and return the path of the neb path floders.
        """
        self.structure = self.relax(atoms=self.atoms, relax_cell=True)
        write("relaxed_atoms.vasp", self.structure)
        OgStructure(structure=self.structure).generate_neb(
            moving_atom_species=self.ion_type,
            num_images=self.neb_point_number,
            r=self.nebor_cutoff,
        )
        return [item for item in self._CWD.iterdir() if item.is_dir()]

    def cal_diffusion_barrier(self, neb_path_floder: Path):
        """
        Calculate the diffusion barrier in the neb path floder
        """
        # * if Ea.txt is in cal_dir, pass
        if Path(f"{neb_path_floder}/Ea.txt").exists():
            logging.info(f"{neb_path_floder}/Ea.txt exists.")
            return

        # * read neb path data from the folder
        neb_data = []
        for i in range(self.neb_point_number + 2):
            tmp_data = read(f"{neb_path_floder}/0{i}.vasp")
            neb_data.append(tmp_data)
        neb_data[0] = self.relax(atoms=neb_data[0])
        neb_data[-1] = self.relax(atoms=neb_data[-1])

        # * neb calculation
        ## * 1. assign calculator to neb data
        images = []
        pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        calc = M3GNetCalculator(potential=pot, stress_weight=1 / 160.21766208)

        for i in neb_data:
            i.calc = calc
            images.append(i)

        ## * 2. run neb calculation
        neb = NEB(images, climb=True, allow_shared_calculator=True)
        qn = FIRE(neb, trajectory=f"{neb_path_floder}/neb.traj")
        try:
            logging.info(f"calculating image {neb_path_floder}")
            qn.run(fmax=0.15)
        except Exception:
            logging.error(f"Error: neb calculation of image {neb_path_floder} failed.")

        ## * 3. get the neb barrier
        traj = read(f"{neb_path_floder}/neb.traj@-{self.neb_point_number+2}:")
        result = []
        for j in traj:
            j.calc = calc
            result.append(j)
        logging.info(f"image {neb_path_floder} calculated.")

        nebtools = NEBTools(result)
        Ef, dE = nebtools.get_barrier()

        ## * 4. write the Ea to Ea.txt
        Path.open(Path(f"{neb_path_floder}/Ea.txt"), "w").write(f"{Ef-dE}")

    def extract_results(self, neb_path_floders: list):
        """
        extract the Ea from the Ea.txt in the neb path floders
        """
        with Path.open(Path(f"{self._CWD}/Ea_result.txt"), "w") as f:
            E_sum = []
            for i in neb_path_floders:
                try:
                    logging.info(f"Extracting 'Ea' from {i}")
                    tmp_Ea = float(
                        Path.open(Path(f"{i}/Ea.txt"), "r")
                        .readlines()[0]
                        .rstrip()
                        .split()[0]
                    )
                    logging.info(f"{i} Ea: {tmp_Ea}")
                    f.write(f"{i}: {tmp_Ea}\n")
                    E_sum.append(tmp_Ea)
                except FileNotFoundError:
                    f.write(f"{i}: Error: File not found\n")
                except ValueError:
                    f.write(f"{i}: Error: Invalid value\n")
            Ea = sum(E_sum) / len(E_sum)
            f.write(f"Ea: {Ea:.3f}")
        return Ea

    def get_diffusion_barrier(self):
        """
        Get the diffusion barrier of the atoms.
        """
        # * 1. generate neb path
        neb_path_floders = self.generate_neb_path()
        logging.info(f"neb_path_floders: {neb_path_floders}")

        # * 2. calculate diffusion barrier
        # for folder in neb_path_floders:
        #     logging.info(f"cal_diffusion_barrier {folder}")
        #     try:
        #         self.cal_diffusion_barrier(folder)
        #     except Exception:
        #         logging.error(f"Error: cal_diffusion_barrier {folder} failed.")
        pool = Pool(nodes=min(self.parallel, len(neb_path_floders)))
        pool.map(self.cal_diffusion_barrier, neb_path_floders)
        logging.info("cal_diffusion_barrier done")

        # * 3. extract results Ea
        return self.extract_results(neb_path_floders)
