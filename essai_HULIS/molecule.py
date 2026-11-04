# ============================ molecule.py ============================
"""Molécule : géométrie, liaisons devinées, types d'atomes Van-Catledge.

Chaîne typique :
    mol = Molecule(symbols, xyz)
    mol.guess_bonds()      # liaisons par rayons covalents (H compris)
    mol.assign_types()     # types Van-Catledge depuis le nb de voisins
    e, c = huckel(mol)     # (dans huckel.py)
"""

import numpy as np
import re
from settings import (ATOMTYPES, RCOV, RCOV_DEFAULT, BOND_TOL,bond_order)

# ---------- contrôle de valence : éléments légers seulement (option b) ----------
# Un dépassement ici = quasi-certainement une erreur (lecture, géométrie
# écrasée, doublon). Les métaux ne sont PAS contrôlés (pas de valence σ
# nette : ferrocène & co doivent passer).
VALENCE_MAX_LIGHT = {
    "H": 1, "B": 4, "C": 4, "N": 4, "O": 3, "F": 1,
    "Si": 4, "P": 5, "S": 4, "Cl": 2,
    "Br": 2, "I": 2,
}

class ValenceError(Exception):
    """Atome avec trop de voisins : la molécule est probablement corrompue."""
    def __init__(self, idx, el, deg):
        self.idx, self.el, self.deg = idx, el, deg
        super().__init__(
            f"atome {idx+1} ({el}) : {deg} voisins, valence max "
            f"{VALENCE_MAX_LIGHT.get(el, '?')} — géométrie ou lecture "
            f"vraisemblablement fausse")


class Atom:
    __slots__ = ("el", "xyz", "atype", "neighbors", "orders", "npi", "flag")
    def __init__(self, el, xyz):
        self.el = el                  # symbole, ex. "C", "N"
        self.xyz = np.asarray(xyz, float)
        self.atype = None             # clé ATOMTYPES, ex. "N:" ; None = hors π
        self.neighbors = []           # indices des voisins (TOUS, H compris)
        self.orders = []              # ordre de liaison en parallèle
        self.npi = 0                  # e- π portés (0 si hors système)
        self.flag = ""                # message d'anomalie non fatale

class Molecule:
    def __init__(self, symbols, xyz):
        self.atoms = [Atom(s, x) for s, x in zip(symbols, xyz)]
        self.warnings = []
        self._pi_idx = None  
        unknown = {s for s in symbols if s not in RCOV}
        for s in sorted(unknown):
            self.warnings.append(
                f"élément '{s}' inconnu : rayon générique {RCOV_DEFAULT} Å")
        self._bonded = False
        self._typed = False

    def __iter__(self):
        return iter(self.atoms)

    # ------------------- liaisons -------------------
    def guess_bonds(self):
        """Liaison si d ≤ BOND_TOL·(r1+r2). O(n²) — largement suffisant."""
        n = len(self.atoms)
        for i in range(n):
            self.atoms[i].neighbors.clear()
            self.atoms[i].orders.clear()
        for i in range(n):
            ri = self.atoms[i].xyz
            for j in range(i + 1, n):
                d = float(np.linalg.norm(ri - self.atoms[j].xyz))
                o = bond_order(self.atoms[i].el, self.atoms[j].el, d)
                if o:
                    self.atoms[i].neighbors.append(j)
                    self.atoms[i].orders.append(o)
                    self.atoms[j].neighbors.append(i)
                    self.atoms[j].orders.append(o)
        self._bonded = True
        # coordination des éléments légers — lève ValenceError
        for i, a in enumerate(self.atoms):
            vmax = VALENCE_MAX_LIGHT.get(a.el)
            if vmax is not None and len(a.neighbors) > vmax:
                raise ValenceError(i, a.el, len(a.neighbors))
        return self

    # ------------------- types Van-Catledge -------------------
    def assign_types(self):
        """Type = entrée ATOMTYPES de l'élément avec valence π = nb TOTAL
        de voisins (H compris). Sinon : hors système π + signalement."""
        self._pi_idx = None  
        if not self._bonded:
            self.guess_bonds()
        for i, a in enumerate(self.atoms):
            deg = len(a.neighbors)               # TOUS les voisins
            cands = [k for k, v in ATOMTYPES.items()
                     if k.rstrip(':') == a.el and v[2] == deg]
            if cands:
                a.atype = cands[0]
                a.npi = ATOMTYPES[a.atype][0]
            else:
                a.atype, a.npi = None, 0
                if a.el not in ("H",):           # H hors π, c'est normal
                    a.flag = f"{a.el} à {deg} voisins : pas de type Van-Catledge, exclu du système π"
        self._typed = True
        return self

    # ------------------- accès pour Hückel / dessin -------------------
    @property
    def pi_idx(self):
        if self._pi_idx is None:
            self.assign_types()          # si pas encore fait
            self._pi_idx = [i for i, a in enumerate(self.atoms) if a.npi > 0]
        return self._pi_idx
    def pi_subgraph(self):
        """(reindex, adj) : reindex[i] = position dans le système π,
        adj = liste d'adjacence du sous-graphe π."""
        if not self._typed:
            self.assign_types()
        idx = self.pi_idx
        reindex = {g: l for l, g in enumerate(idx)}
        adj = [[] for _ in idx]
        for l, g in enumerate(idx):
            a = self.atoms[g]
            adj[l] = [reindex[j] for j in a.neighbors if j in reindex]
        return reindex, adj

    def npz(self):
        """Nombre total d'électrons π."""
        return sum(a.npi for a in self.atoms)

    def summary(self):
        if not self._typed:
            self.assign_types()
        lines = []
        for i, a in enumerate(self.atoms):
            t = a.atype if a.atype else f"— ({a.flag})" if a.flag else "—"
            lines.append(f"{i+1:>3} {a.el:<2} {len(a.neighbors)} voisins  "
                         f"type {t}  nπ={a.npi}")
        for w in self.warnings:
            lines.append("warning: " + w)
        return "\n".join(lines)


# molecule.py (partie 1 : lecture)
# ---------- lecture tolérante ----------

_LINE = re.compile(r"^\s*([A-Z][a-z]?)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")

def read_xyz(path=None, text=None):
    """Lecteur tolérant : la 1re ligne peut être un compte, un titre,
    ou déjà un atome. Renvoie (symbols, xyz) prêts pour Molecule()."""
    if path is not None:
        with open(path) as f:
            text = f.read()
    symbols, xyz = [], []
    for line in text.splitlines():
        m = _LINE.match(line)
        if not m:
            continue                          # compte, titre, blanc → skip
        g = m.groups()
        try:
            pos = (float(g[1]), float(g[2]), float(g[3]))
        except ValueError:
            continue
        sym = g[0].capitalize()               # "CL" -> "Cl"
        if sym in RCOV:                        # élément CONNU (pas ATOMTYPES !)
            symbols.append(sym)
            xyz.append(pos)
        else:
            print(f"warning read_xyz: élément '{sym}' ignoré (inconnu)")
    return symbols, xyz

