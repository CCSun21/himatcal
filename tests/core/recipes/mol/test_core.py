from __future__ import annotations

"""Test suite for mol core functionality."""
import json
import os
import tempfile
from typing import Callable, Literal
from unittest.mock import Mock, patch

import pytest
from ase.atoms import Atoms
from rdkit import Chem

from himatcal.recipes.mol.core import (
    CASNumber,
    cas2xyz,
    cas_to_smiles,
    relax_mol,
    sanitize_mol,
)


def create_test_molecule(smiles: str) -> Chem.Mol:
    """Create a test molecule from SMILES."""
    return Chem.MolFromSmiles(smiles)


def create_mock_response(
    content: str | None = None, status_code: int = 200, error: Exception | None = None
) -> Mock:
    """Create a mock response for testing."""
    mock_response = Mock()
    mock_response.text = content
    mock_response.status_code = status_code
    if error:
        mock_response.raise_for_status.side_effect = error
    else:
        mock_response.raise_for_status.return_value = None
    return mock_response


def assert_valid_xyz_result(cas_id: str):
    """Assert that cas2xyz returns a valid XYZ result for the given CAS ID."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        current_dir = os.getcwd()
        try:
            os.chdir(tmp_dir)
            result = cas2xyz(cas_id, relax_atoms=False)
            assert result is not None
            assert isinstance(result, str)
            assert "xyz" in result.lower()
            # Verify the file is created in the temporary directory
            assert os.path.exists(f"{cas_id}.xyz")
        finally:
            os.chdir(current_dir)


def assert_sanitized_mol(smiles: str, validation_func: Callable[[Chem.Mol], None]):
    """Test molecule sanitization with custom validation.

    Args:
        smiles: SMILES string of the molecule to test
        validation_func: Function that takes a sanitized molecule and performs assertions
    """
    mol = create_test_molecule(smiles)
    result = sanitize_mol(mol)
    assert result is not None
    validation_func(result)


def assert_cas_to_smiles(
    cas_id: str, expected_smiles: str, source: Literal["cirpy", "pubchem"] = "cirpy"
):
    """Test CAS to SMILES conversion with expected result.

    Args:
        cas_id: CAS number to test
        expected_smiles: Expected SMILES string result
        source: Source to use for conversion
    """
    with patch("requests.get") as mock_get:
        if source == "pubchem":
            mock_response = create_mock_response(
                json.dumps(
                    {
                        "PropertyTable": {
                            "Properties": [{"IsomericSMILES": expected_smiles}]
                        }
                    }
                )
            )
            mock_response.json = lambda: {
                "PropertyTable": {"Properties": [{"IsomericSMILES": expected_smiles}]}
            }
        else:
            mock_response = create_mock_response(expected_smiles)
        mock_get.return_value = mock_response
        result = cas_to_smiles(cas_id, source=source)
        assert result == expected_smiles


class TestSanitizeMol:
    """Test suite for sanitize_mol function."""

    def test_none_input(self):
        """Test handling of None input."""
        assert sanitize_mol(None) is None

    def test_fluorine_handling(self):
        """Test handling of fluorine atoms."""

        def validate_fluorine(mol: Chem.Mol):
            f_atom = next(atom for atom in mol.GetAtoms() if atom.GetSymbol() == "F")
            assert f_atom.GetDegree() == 1
            assert f_atom.GetFormalCharge() == -1

        assert_sanitized_mol("CC(F)(F)F", lambda m: validate_fluorine(m))

    def test_charge_balancing(self):
        """Test charge balancing for ionic compounds."""

        def validate_charge(mol: Chem.Mol):
            total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            assert total_charge == 0

        assert_sanitized_mol("[Na+].[Cl-]", lambda m: validate_charge(m))


class TestCas2xyz:
    """Test suite for cas2xyz function."""

    def test_valid_cas_ethanol(self):
        """Test cas2xyz with valid CAS for ethanol."""
        assert_valid_xyz_result("64-17-5")

    def test_valid_cas_benzene(self):
        """Test cas2xyz with valid CAS for benzene."""
        assert_valid_xyz_result("71-43-2")

    def test_invalid_cas(self):
        """Test cas2xyz with invalid CAS."""
        with pytest.raises(ValueError):
            cas2xyz("invalid-cas", relax_atoms=False)

    @patch("himatcal.recipes.mol.core.cas_to_smiles")
    def test_smiles_conversion_path(self, mock_cas_to_smiles):
        """Test the SMILES conversion path in cas2xyz."""
        mock_cas_to_smiles.return_value = "CCO"  # Ethanol SMILES
        assert_valid_xyz_result("64-17-5")
        assert mock_cas_to_smiles.called


class TestCasToSmiles:
    """Test suite for cas_to_smiles function."""

    def test_valid_cas_ethanol_cirpy(self):
        """Test CAS to SMILES conversion for ethanol using cirpy."""
        assert_cas_to_smiles("64-17-5", "CCO", source="cirpy")

    def test_valid_cas_benzene_pubchem(self):
        """Test CAS to SMILES conversion for benzene using pubchem."""
        assert_cas_to_smiles("71-43-2", "c1ccccc1", source="pubchem")

    def test_invalid_cas(self):
        """Test CAS to SMILES conversion with invalid CAS."""
        with patch("requests.get") as mock_get:
            mock_get.return_value = create_mock_response(error=Exception("Invalid CAS"))
            result = cas_to_smiles("invalid-cas", source="cirpy")
            assert result is None


class TestCASNumber:
    """Test suite for CASNumber model."""

    def test_valid_cas_short(self):
        """Test validation of short valid CAS number."""
        model = CASNumber(cas_number="64-17-5")
        assert model.cas_number == "64-17-5"

    def test_valid_cas_medium(self):
        """Test validation of medium-length valid CAS number."""
        model = CASNumber(cas_number="100-41-4")
        assert model.cas_number == "100-41-4"

    def test_valid_cas_long(self):
        """Test validation of long valid CAS number."""
        model = CASNumber(cas_number="1234-56-7")
        assert model.cas_number == "1234-56-7"

    def test_invalid_cas_format(self):
        """Test validation of invalid CAS number format."""
        with pytest.raises(ValueError):
            CASNumber(cas_number="1-23-4")

    def test_invalid_cas_too_long(self):
        """Test validation of too-long CAS number."""
        with pytest.raises(ValueError):
            CASNumber(cas_number="1234567-78-9")

    def test_invalid_cas_with_letters(self):
        """Test validation of CAS number containing letters."""
        with pytest.raises(ValueError):
            CASNumber(cas_number="abc-12-3")


class TestRelaxMol:
    """Test suite for relax_mol function."""

    def test_crest_method(self):
        """Test relaxation using CREST method."""
        mock_atoms = Atoms("H2O")
        with patch("himatcal.recipes.mol.core.crest_relax") as mock_relax:
            mock_relax.return_value = mock_atoms
            result = relax_mol(mock_atoms, method="crest")
            assert result == mock_atoms
            mock_relax.assert_called_once()

    def test_invalid_method(self):
        """Test handling of invalid relaxation method."""
        mock_atoms = Atoms("H2O")
        with pytest.raises(ValueError, match=".*not supported"):
            relax_mol(mock_atoms, method="invalid_method")
