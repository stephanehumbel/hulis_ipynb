# pipeline.py — raccourcis du labo
import molecule, huckel, huckel_analysis, chch_pubchem

def mol_from_xyz(path=None, text=None):
    syms, xyz = molecule.read_xyz(path=path, text=text)
    return molecule.Molecule(syms, xyz)

def mol_from_name(name):
    """'fulvene' → Molecule complète, prête pour Hückel."""
    m = chch_pubchem.mol_from_name(name)
    return mol_from_xyz(text=chch_pubchem.mol_to_xyz(m, name))

def full_analysis(name_or_mol):
    """Tout-en-un : molécule → Hückel → rapport. Retourne (mol, res, texte)."""
    mol = name_or_mol if isinstance(name_or_mol, molecule.Molecule) else mol_from_name(name_or_mol)
    res = huckel.solve_huckel(mol)
    texte, _ = huckel_analysis.format_report(mol, res)
    return mol, res, texte
