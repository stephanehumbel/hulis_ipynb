import numpy as np
from settings import ATOMTYPES,kval

BETA = -1.0   # unité d'énergie


def build_hamiltonian(mol):
    """Hπ en unités de β (β = −1 : OM liantes → e < 0, cf. docstring de huckel)."""
    if not mol._bonded:
        mol.guess_bonds()
    if not mol._typed:
        mol.assign_types()
    idx = mol.pi_idx
    n = len(idx)
    if n == 0:
        raise ValueError("Aucun atome π (npi > 0) dans la molécule.")
    reindex = {g: l for l, g in enumerate(idx)}
    H = np.zeros((n, n))
    for l, g in enumerate(idx):
        a = mol.atoms[g]
        H[l, l] = -ATOMTYPES[a.atype][1]        # −h_i  (h_C=0, h_N=0.51)
        for j in a.neighbors:
            if j in reindex:
                #H[l, reindex[j]] = -1.0          # β_ij = −k_ij ; k=1 par défaut,
                H[l, reindex[j]] = -kval(a.atype, mol.atoms[j].atype)

    return H, idx                                 # → k de settings si ton cours le définit

def print_hamiltonian(mol, H):
    """Affiche Hπ avec les symboles des atomes π en ligne/colonne."""
    idx = mol.pi_idx
    labels = [mol.atoms[g].el for g in idx]
    n = len(labels)
    w = max(len(s) for s in labels) + 8
    header = "" * (w + 2) + "".join(f"{s:>{w}}" for s in labels)
    print("Hamiltonian\n",header + "\n" + "-" * len(header))
    for i, si in enumerate(labels):
        row = f"{si:>2} |"
        for j in range(n):
            row += f"{H[i, j]:>{w}.3f}"
        print(row)

def solve_huckel(mol):
    """Diagonalisation + remplissage des OM π.

    Retourne un dict : e (énergies), c (coefficients, colonnes = OM),
    occ (occupations), Epi (énergie π totale), idx (indices des atomes π).
    """
    H, idx = build_hamiltonian(mol)
    e, c = np.linalg.eigh(H)                   # triées croissantes

    ne = int(sum(a.npi for a in mol.atoms))
    occ = np.zeros(len(e), dtype=int)
    k = 0
    while ne >= 2:
        occ[k] = 2; ne -= 2; k += 1
    if ne == 1:
        occ[k] = 1

    return {"e": e, "c": c, "occ": occ,
            "Epi": float(np.dot(e, occ)), "idx": idx}