# ============================ settings.py ============================
# ---------- paramètres & constantes ----------
params = dict(zoom=1.0, orb_color='#FF0000', orb_scale=1.0, angle_step=15.0)
BOND = 1.2          # longueur de liaison imposée (unités du graphe)

# ============ Van-Catledge, J. Org. Chem. 1980, 45, 4801 ==============
# αX = α + hX·β ;  βXY = kXY·β
# Type 1 = 1 électron π ; type 2 = paire libre (2 électrons π)
# el : (n_π, h, valence π, couleur)
ATOMTYPES = {
    "C":  (1,  0.00, 3, '#333333'),
    "B":  (1, -0.45, 3, '#d9a441'),
    "N":  (1,  0.51, 2, '#3050c8'),   # N1 pyridinique : lone pair dans le plan
    "N:": (2,  1.37, 3, '#3050c8'),   # N2 pyrrole / amine : lone pair π
    "O":  (1,  0.97, 1, '#c83030'),   # O1 carbonyle
    "O:": (2,  2.09, 2, '#c83030'),   # O2 éther / furane
    "F":  (2,  2.71, 0, '#66c830'),   # substituant, 0 voisin π
    "Si": (1,  0.00, 3, '#a08060'),   # silabenzène
    "P":  (1,  0.19, 2, '#c88030'),   # P1 phospha-pyridine
    "P:": (2,  0.75, 3, '#c88030'),   # P2
    "S":  (1,  0.46, 1, '#c8c030'),   # S1 thiocarbonyle
    "S:": (2,  1.11, 2, '#c8c030'),   # S2 thioéther / thiophène
    "Cl": (2,  1.48, 0, '#30a060'),   # substituant, 0 voisin π
}

# ---- Table II Van-Catledge 1980 : demi-matrice kXY ----
# Ordre : C B N1 N2 O1 O2 F Si P1 P2 S1 S2 Cl
# N1 = N pyridine (1e π),  N2 = N pyrrole (2e)
# O1 = O carbonyle (1e),   O2 = O éther (2e), idem P1/P2, S1/S2
KHALF_KEYS = ["C", "B", "N", "N:", "O", "O:", "F", "Si", "P", "P:", "S", "S:", "Cl"]
KHALF = [
    [1.00],                                              # C
    [0.73, 0.87],                                        # B
    [1.02, 0.66, 1.09],                                  # N1
    [0.89, 0.53, 0.99, 0.98],                            # N2
    [1.06, 0.60, 1.14, 1.13, 1.26],                      # O1
    [0.66, 0.35, 0.80, 0.89, 1.02, 0.95],                # O2
    [0.52, 0.26, 0.65, 0.77, 0.92, 0.94, 1.04],          # F
    [0.75, 0.57, 0.72, 0.43, 0.65, 0.24, 0.17, 0.64],    # Si
    [0.77, 0.53, 0.78, 0.55, 0.75, 0.31, 0.21, 0.62, 0.63],          # P1
    [0.76, 0.54, 0.81, 0.64, 0.82, 0.39, 0.22, 0.52, 0.58, 0.63],     # P2
    [0.81, 0.51, 0.83, 0.68, 0.84, 0.43, 0.28, 0.61, 0.65, 0.65, 0.68],# S1
    [0.69, 0.44, 0.78, 0.73, 0.85, 0.54, 0.32, 0.40, 0.48, 0.60, 0.58, 0.63],  # S2
    [0.62, 0.41, 0.77, 0.80, 0.88, 0.70, 0.51, 0.34, 0.35, 0.55, 0.52, 0.59, 0.68],  # Cl
]
KHALF_DICT = {}
for _i, _k1 in enumerate(KHALF_KEYS):
    for _j, _k2 in enumerate(KHALF_KEYS[:_i + 1]):
        KHALF_DICT[(_k1, _k2)] = KHALF_DICT[(_k2, _k1)] = KHALF[_i][_j]

def kval(el1, el2):
    """kXY Van-Catledge, 1.0 par défaut (liaison C-C inconnue)."""
    return KHALF_DICT.get((el1, el2), 1.0)

# ============ Détection des liaisons : rayons covalents ==============
# Cordero et al., Dalton Trans. 2008, 2832-2838.
RCOV = {
    "H": 0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53,
    "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19,
    "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64,
    "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39,
    "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03,
    "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96,
    "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90,
    "Yb": 1.87, "Lu": 1.87, "Hf": 1.75, "Ta": 1.70, "W": 1.62,
    "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36,
    "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.40,
    "At": 1.50, "Rn": 1.50,
    "Fr": 2.60, "Ra": 2.21, "Ac": 2.15, "Th": 2.06, "Pa": 2.00,
    "U": 1.96, "Np": 1.90, "Pu": 1.87, "Am": 1.80, "Cm": 1.69,
}
RCOV_DEFAULT = 1.5

BOND_TOL = 1.10
MULTIPLE_CUT = 0.95

def bond_order(el1, el2, d):
    """0 = pas de liaison, 1 = simple, 2 = multiple (double/triple)."""
    r1 = RCOV.get(el1, RCOV_DEFAULT)
    r2 = RCOV.get(el2, RCOV_DEFAULT)
    rs = r1 + r2
    if d > BOND_TOL * rs:
        return 0
    if d < MULTIPLE_CUT * rs:
        return 2
    return 1