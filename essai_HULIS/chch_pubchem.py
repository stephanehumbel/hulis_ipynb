import urllib.request, urllib.parse
from rdkit import Chem

def mol_from_name(name):
    """Nom -> molécule RDKit 3D avec H explicites, via PubChem."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           f"{urllib.parse.quote(name)}/SDF?record_type=3d")
    try:
        sdf = urllib.request.urlopen(url).read()
    except urllib.error.HTTPError as e:
        raise ValueError(f"PubChem : '{name}' introuvable ({e.code}).") from None
    mol = Chem.MolFromMolBlock(sdf.decode())
    if mol is None:
        raise ValueError(f"SDF illisible pour '{name}'.")
    mol = Chem.AddHs(mol, addCoords=True)   # H placés sur la géométrie 3D existante
    return mol

def mol_to_xyz(mol, comment=""):
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), comment]
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol()} {p.x:.6f} {p.y:.6f} {p.z:.6f}")
    return "\n".join(lines)
