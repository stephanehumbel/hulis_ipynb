import numpy as np

# ============================================================
#  Calculs
# ============================================================

def population_analysis(mol, res):
    idx, C, occ, e = res["idx"], res["c"], res["occ"], res["e"]

    n_e = sum(mol.atoms[g].npi for g in idx)
    P = C @ np.diag(occ) @ C.T

    charges = {}
    for l, g in enumerate(idx):
        charges[g] = mol.atoms[g].npi - P[l, l]

    # --- Ordres de liaison π : uniquement pour les paires VOISINES ---
    bond_orders = {}
    pos = {g: l for l, g in enumerate(idx)}   # index π -> position dans la base
    for l1 in range(len(idx)):
        g1 = idx[l1]
        for g2 in mol.atoms[g1].neighbors:
            if g2 > g1 and g2 in pos:          # voisin σ, une seule fois, et participant au π
                bond_orders[(g1, g2)] = P[pos[g1], pos[g2]]

    return {
        "idx": idx, "e": e, "occ": occ, "Epi": res["Epi"],
        "n_e": n_e, "P": P, "charges": charges, "bond_orders": bond_orders,
    }

# ============================================================
#  Affichages texte (retourne une chaîne, facile à mettre dans un panel)
# ============================================================

def format_report(mol, res, res_analysis=None):
    """
    Si res_analysis n'est pas fourni, lance l'analyse automatiquement.
    Retourne (texte, res_analysis).
    """
    if res_analysis is None:
        res_analysis = population_analysis(mol, res)

    lines = []
    lines.append("Énergies π (unités β) :")
    for k, (ek, ok) in enumerate(zip(res_analysis["e"], res_analysis["occ"])):
        bullet = f"● {int(ok)}e-" if ok > 0 else "○"
        lines.append(f"  OM {k+1} : ε = {-ek:+.3f} β   {bullet}")

    lines.append("\nCharges π :")
    for g, q in res_analysis["charges"].items():
        a = mol.atoms[g]
        lines.append(f"  {a.el}{g+1} : q_π = {q:+.3f}")

    lines.append("\nOrdres de liaison π :")
    for (g1, g2), p in res_analysis["bond_orders"].items():
        a1, a2 = mol.atoms[g1], mol.atoms[g2]
        lines.append(f"  {a1.el}{g1+1}–{a2.el}{g2+1} : p_π = {p:.3f}")

    lines.append(f"\nÉnergie π totale : {res_analysis['Epi']:+.3f} β")
    return "\n".join(lines), res_analysis


   