"""
Molecular Similarity Search with RDKit

This script demonstrates how to compute molecular similarity using the RDKit library.
It covers:
1. Generating Morgan (circular) fingerprints for molecules.
2. Computing Tanimoto similarity between two molecules.
3. Performing a similarity search against a database of molecules.

Author: Ajay Khanna
Date: Sep.19.2023
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def print_banner(text):
    """
    Print a banner with the given text.

    Args:
    - text (str): Text to display in the banner.
    """
    banner = f"{'=' * 20}\n{text}\n{'=' * 20}"
    print(banner)


# Call the function to print the banner
# Generate Fingerprints
def generate_fingerprint(molecule, radius=2):
    """
    Generate Morgan fingerprint for a given molecule.

    Args:
    - molecule (rdkit.Chem.rdchem.Mol): RDKit molecule object.
    - radius (int): Radius for the Morgan fingerprint.

    Returns:
    - rdkit.DataStructs.cDataStructs.ExplicitBitVect: Fingerprint of the molecule.
    """
    return AllChem.GetMorganFingerprintAsBitVect(molecule, radius)


# Compute Similarity
def compute_similarity(fp1, fp2):
    """
    Compute Tanimoto similarity between two fingerprints.

    Args:
    - fp1, fp2 (rdkit.DataStructs.cDataStructs.ExplicitBitVect): Fingerprints to compare.

    Returns:
    - float: Tanimoto similarity.
    """
    return DataStructs.TanimotoSimilarity(fp1, fp2)


# Similarity Search
def similarity_search(query_fp, database_mols, radius=2):
    """
    Perform a similarity search against a database of molecules.

    Args:
    - query_fp (rdkit.DataStructs.cDataStructs.ExplicitBitVect): Fingerprint of the query molecule.
    - database_mols (list): List of RDKit molecule objects in the database.
    - radius (int): Radius for the Morgan fingerprint.

    Returns:
    - list: List of Tanimoto similarities.
    """
    similarities = []
    for mol in database_mols:
        fp = generate_fingerprint(mol, radius)
        sim = compute_similarity(query_fp, fp)
        similarities.append(sim)
    return similarities


if __name__ == "__main__":
    # Print banner
    print_banner("Molecular Similarity Search with RDKit")
    # Sample database of molecules
    smiles_list = [
        "CCO",
        "CCN",
        "CCF",
        "CCC",
        "CCCl",
        "NCN",
        "NCC",
        "NCO",
        "NCF",
        "NCl",
        "CN",
        "CO",
        "CF",
        "Cl",
    ]
    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    # Query molecule
    query_mol = Chem.MolFromSmiles("CCBr")
    query_fp = generate_fingerprint(query_mol)

    # Search database
    results = similarity_search(query_fp, mols)

    # Print results
    for smiles, sim in zip(smiles_list, results):
        print(f"{smiles}: {sim:.2f}")
