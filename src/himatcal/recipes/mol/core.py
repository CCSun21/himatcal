from __future__ import annotations

import json
import logging
import re
import urllib.request
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import requests
from ase.io import write
from chemspipy import ChemSpider
from pydantic import BaseModel, field_validator
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from himatcal import SETTINGS
from himatcal.recipes.crest.core import relax as crest_relax
from himatcal.utils.rdkit.core import rdkit2ase

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms
    from rdkit.Chem.rdchem import Mol


def sanitize_mol(mol):
    """
    Sanitize and clean up a molecular structure.

    Args:
        mol (rdkit.Chem.Mol): The molecular structure to sanitize.

    Returns:
        rdkit.Chem.Mol or None: The sanitized molecular structure if successful, None otherwise.
    """
    if mol is None:
        logging.warning("Input molecule is None")
        return None

    try:
        # Process each atom before any sanitization
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol == "F":
                # Remove all bonds except one
                while atom.GetDegree() > 1:
                    bonds = atom.GetBonds()
                    mol.RemoveBond(bonds[0].GetBeginAtomIdx(), bonds[0].GetEndAtomIdx())

                # Clear any existing formal charge
                atom.SetFormalCharge(0)
                atom.UpdatePropertyCache(strict=False)

        # Initial cleanup with minimal sanitization
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ADJUSTHS)

        # Do a full sanitization without charges
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_SYMMRINGS
            | Chem.SANITIZE_KEKULIZE
            | Chem.SANITIZE_SETAROMATICITY
            | Chem.SANITIZE_CLEANUP,
        )

        # Post-process charges
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "F":
                atom.SetFormalCharge(0)  # Keep fluorine neutral for 3D embedding

            # Update atom properties
            atom.UpdatePropertyCache(strict=False)

        return mol

    except Exception as e:
        logging.error(f"Sanitization failed: {str(e)}")
        return None


def get_molecular_structure(
    molecular_cas: str,
    write_mol: bool = True,
    chemspider_api: str = SETTINGS.CHEMSPIDER_API_KEY,
):
    """
    Get molecular structure from CAS number, using chemspipy from RSC ChemSpider.

    Args:
        molecular_cas (str): The CAS number of the molecule.
        write_mol (bool): Whether to write the molecule to XYZ file. Defaults to True.
        chemspider_api (str): The ChemSpider API key.

    Returns:
        rdkit.Chem.Mol or None: The molecular structure if successful, None otherwise.

    Raises:
        ValueError: If the molecular structure cannot be retrieved or is invalid.
    """
    cs = ChemSpider(chemspider_api)
    try:
        results = cs.search(molecular_cas)
        if not results:
            return None

        c1 = results[0]
        mol_file = StringIO(c1.mol_3d)
        mol = Chem.MolFromMolBlock(
            mol_file.getvalue(),
            removeHs=False,
            sanitize=False,
        )

        mol = sanitize_mol(mol)
        if mol is None:
            return None

        mol = Chem.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
        rdForceFieldHelpers.UFFOptimizeMolecule(mol)  # Optimize using UFF

        if write_mol:
            try:
                Chem.MolToXYZFile(mol, f"{molecular_cas}.xyz")
            except Exception as e:
                logging.error(f"警告: 无法写入XYZ文件: {e!s}")

        return mol

    except Exception as e:
        logging.error(f"错误: 从ChemSpider获取结构时发生错误: {e!s}")
        return None


def get_with_retry(
    url: str, max_retries: int = 3, backoff_factor: float = 0.3
) -> requests.Response:
    """
    Makes an HTTP GET request with exponential backoff retries
    """
    session = requests.Session()
    retry_strategy = Retry(
        total=max_retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    return session.get(url, timeout=10)


def consumeApi(urlPath: str) -> str | None:
    """
    Consume an API endpoint with retry logic
    """
    try:
        response = get_with_retry(urlPath)
        return response.text if response.status_code == 200 else None
    except Exception as e:
        logging.warning(f"Failed to consume API {urlPath}: {e}")
        return None


def get_pubchem_mol(cas_id: str) -> Mol | None:
    """
    Retrieve molecular structure from PubChem using CAS ID.

    Args:
        cas_id (str): The CAS ID of the molecule

    Returns:
        Mol | None: The molecular structure if successful, None otherwise
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_id}/cids/JSON"
        response = requests.get(url)
        if response.status_code != 200:
            return None

        data = response.json()
        if "IdentifierList" not in data:
            return None

        cid = data["IdentifierList"]["CID"][0]

        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/SDF"
        sdf_response = requests.get(sdf_url)
        if sdf_response.status_code != 200:
            return None

        mol = Chem.MolFromMolBlock(sdf_response.text, removeHs=False, sanitize=False)
        return sanitize_mol(mol) if mol is not None else None
    except Exception as e:
        logging.warning("Failed to retrieve structure from PubChem: %s", e)
        return None


def cas2xyz(CAS_ID: str, relax_atoms: bool = True) -> str | None:
    """
    Converts a CAS ID into an XYZ file format representation of the corresponding molecule.

    This function first tries to convert CAS to SMILES, then retrieves molecular data from multiple sources,
    constructs a molecular structure, and optionally relaxes the atomic positions before saving the structure to an XYZ file.

    Args:
        CAS_ID (str): The Chemical Abstracts Service identifier for the desired molecule.
        relax_atoms (bool): A flag indicating whether to relax the atomic positions before saving. Defaults to True.

    Returns:
        Optional[str]: Contents of the XYZ file if successful, None otherwise

    Raises:
        ValueError: If no valid molecular structure could be obtained from any source.
    """
    mol = None
    success_source = None

    # First try getting SMILES from CAS
    try:
        logging.info("Attempting to convert CAS to SMILES...")
        if smiles := cas_to_smiles(CAS_ID):
            mol = Chem.MolFromSmiles(smiles)
            if mol := sanitize_mol(mol):
                success_source = "SMILES conversion"
                logging.info(
                    "Successfully converted CAS to SMILES and created molecule"
                )
    except Exception as e:
        logging.warning(f"Failed to convert CAS to SMILES: {e}")

    # If SMILES conversion failed, try PubChem
    if not mol:
        try:
            logging.info("Attempting to retrieve structure from PubChem...")
            mol = get_pubchem_mol(CAS_ID)
            if mol:
                success_source = "PubChem"
                logging.info("Successfully retrieved structure from PubChem")
        except Exception as e:
            logging.warning("Failed to retrieve from PubChem: %s", e)

    # Try all web sources if previous methods failed
    sources = [
        f"https://commonchemistry.cas.org/api/detail?cas_rn={CAS_ID}",
        f"https://www.chemicalbook.com/CAS/mol/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20210305/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20210111/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20180601/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20180713/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20180808/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20150408/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20200515/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20211123/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/20200331/MOL/{CAS_ID}.mol",
        f"https://www.chemicalbook.com/CAS/MOL/{CAS_ID}.mol",
    ]

    if not mol:
        for source in sources:
            try:
                if "commonchemistry" in source:
                    api_result = consumeApi(source)
                    if api_result is not None:
                        result_dict = json.loads(api_result)
                        if "inchi" not in result_dict:
                            logging.warning(
                                "No InChI found in response from %s", source
                            )
                            continue
                        mol = Chem.MolFromInchi(result_dict["inchi"])
                        if mol is None:
                            logging.warning(
                                "Failed to create molecule from InChI from %s", source
                            )
                            continue
                else:
                    request = urllib.request.Request(
                        source, headers={"User-Agent": "Mozilla/5.0"}
                    )
                    response = urllib.request.urlopen(request)
                    if response.status == 200:
                        mol_file_content = response.read().decode("utf-8")
                        mol = Chem.MolFromMolBlock(mol_file_content, removeHs=False)
                        if mol is None:
                            logging.warning(
                                "Failed to create molecule from MOL file from %s",
                                source,
                            )
                            continue

                mol = sanitize_mol(mol)
                if mol:
                    success_source = source
                    logging.info(
                        "Successfully retrieved and sanitized data from %s", source
                    )
                    break
                logging.warning("Sanitization failed for molecule from %s", source)
            except Exception as e:
                logging.warning(
                    "Failed to retrieve or process data from %s: %s", source, e
                )
                continue

    # Try ChemSpider as last resort
    if not mol:
        try:
            logging.info("Attempting to retrieve structure from ChemSpider...")
            mol = get_molecular_structure(CAS_ID)
            if mol:
                success_source = "ChemSpider"
                logging.info("Successfully retrieved structure from ChemSpider")
        except Exception as e:
            logging.warning("Failed to retrieve from ChemSpider: %s", e)

    if not mol:
        raise ValueError(
            f"Could not create valid molecular structure for CAS ID {CAS_ID} from any source. "
            f"Last successful source attempt: {success_source}"
        )

    try:
        mol = Chem.AddHs(mol, addCoords=True)
        logging.info("Successfully added hydrogens")

        # Try different 3D embedding methods
        try:
            # First try ETKDG
            params = rdDistGeom.ETKDGv3()
            params.randomSeed = 42  # For reproducibility
            rdDistGeom.EmbedMolecule(mol, params)
            logging.info("Successfully embedded molecule using ETKDG")
        except Exception as e:
            logging.warning(f"ETKDG embedding failed: {e}, trying distance geometry...")
            try:
                # Fallback to basic distance geometry
                rdDistGeom.EmbedMolecule(mol, randomSeed=42)
                logging.info("Successfully embedded molecule using distance geometry")
            except Exception as e2:
                logging.error(f"All embedding methods failed: {e2}")
                return None

        # Try UFF optimization first
        try:
            rdForceFieldHelpers.UFFOptimizeMolecule(mol, maxIters=4000)
            logging.info("Successfully optimized molecule using UFF")
        except Exception as e:
            logging.warning(f"UFF optimization failed: {e}, trying MMFF94...")
            try:
                # Fallback to MMFF94
                rdForceFieldHelpers.MMFFOptimizeMolecule(mol, maxIters=4000)
                logging.info("Successfully optimized molecule using MMFF94")
            except Exception as e2:
                logging.error(f"All force field optimizations failed: {e2}")
                return None

        if relax_atoms:
            atoms = rdkit2ase(mol)
            logging.info("Successfully converted to ASE atoms")

            atoms_relaxed = crest_relax(atoms)
            if atoms_relaxed is None:
                logging.warning("CREST relaxation failed, using unrelaxed structure")
                write(f"{CAS_ID}.xyz", atoms)
            else:
                logging.info("Successfully relaxed structure with CREST")
                filename = f"{CAS_ID}.xyz"
                write(filename, atoms_relaxed)
                with Path(filename).open() as f:
                    content = f.read()
                # 确保文件名包含在输出中
                return f"{filename}\n{content}"
        else:
            filename = f"{CAS_ID}.xyz"
            Chem.MolToXYZFile(mol, filename)
            logging.info("Successfully wrote XYZ file without relaxation")
            with Path(filename).open() as f:
                content = f.read()
            # 确保文件名包含在输出中
            return f"{filename}\n{content}"
    except Exception as e:
        raise ValueError(f"Error during molecule processing: {e}") from e


# TODO: include PubGrep ( https://github.com/grimme-lab/PubGrep )

# TODO: Methods or flow to relax the mols


def relax_mol(mol: Atoms, chg: int = 0, mult: int = 1, method: str = "crest", **kwargs):
    """
    Relax the molecular structure using the specified method.

    This function relaxes the atomic positions of the provided molecular structure using the specified method.

    Args:
        mol (Atoms): The molecular structure to relax.
        method (str): The relaxation method to use. Defaults to 'crest'.
        **kwargs: Additional keyword arguments for the relaxation method.

    Returns:
        Atoms: The relaxed molecular structure.

    Raises:
        ValueError: If the relaxation method is not supported.

    Examples:
        relax_mol(mol, method='crest')  # Relaxes the molecular structure using the CREST method.
    """
    if method == "crest":
        return crest_relax(mol, chg=chg, mult=mult, **kwargs)
    elif method == "optrs":
        from himatcal.recipes.optrs.core import relax as optrs_relax

        return optrs_relax(mol)
    elif method == "orca":
        from himatcal.recipes.orca.core import relax_job as orca_relax

        orca_result = orca_relax(mol, charge=chg, spin_multiplicity=mult, **kwargs)

        return orca_result.atoms
    else:
        raise ValueError(f"Relaxation method '{method}' is not supported")


class CASNumber(BaseModel):
    cas_number: str

    @field_validator("cas_number")
    def validate_cas_number(cls, value):
        """Validate CAS number format and length."""
        # 验证基本格式
        pattern = re.compile(r"^\d{2,7}-\d{2}-\d{1}$")
        if not re.match(pattern, value):
            raise ValueError(
                "Invalid CAS number format. It should be 2-7 digits followed by a hyphen, then 2 digits, and another hyphen followed by 1 digit."
            )

        # 验证第一部分的长度
        parts = value.split("-")
        if len(parts[0]) > 7:  # 检查第一部分是否超过7位数字
            raise ValueError("Invalid CAS number: first part cannot exceed 7 digits")

        # 验证整体长度
        total_digits = sum(len(part) for part in parts)
        if total_digits > 10:  # CAS号的数字总数不应超过10位
            raise ValueError("Invalid CAS number: total length cannot exceed 10 digits")

        return value


def cas_to_smiles(
    cas_number: str, source: Literal["auto", "pubchem", "cirpy"] = "auto"
):
    """
    convert cas to smiles using nci api, with fallback options
    """
    try:
        # * validize the cas_number
        cas_number = CASNumber(cas_number=cas_number).cas_number
    except ValueError:
        return None

    if source in ["auto", "cirpy"]:
        url = f"https://cactus.nci.nih.gov/chemical/structure/{cas_number}/smiles"
        try:
            response = requests.get(url)
            response.raise_for_status()
            smiles = response.text.strip()
            if smiles:
                return smiles
            if source == "cirpy":
                return None
        except requests.exceptions.RequestException as e:
            logging.error(f"Error fetching data for CAS {cas_number} from cirpy: {e}")
            if source == "cirpy":
                return None

    if source in ["auto", "pubchem"]:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/property/IsomericSMILES/JSON"
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
        except (
            requests.exceptions.RequestException,
            KeyError,
            IndexError,
            json.JSONDecodeError,
        ) as e:
            logging.error(f"Error fetching data from PubChem: {e}")

    return None
