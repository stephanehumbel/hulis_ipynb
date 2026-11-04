import numpy as np

class LocalizedStructure:
    """
    Structure de Lewis : bloc-diagonale dans l'hamiltonien Hückel du parent.
    blocks : liste de blocs, chaque bloc = (indices_atomes, n_electrons)
    ex: ('liaison', [0,1], 2), ('paire Libre', [3], 2), ('sextet', [4,5,6,7,8,9], 6), ('radical', [10], 1)
    """

    def __init__(self, parent, blocks, label=None):
        self.parent = parent          # Molecule parente (possède H)
        self.blocks = blocks
        self.label = label
        self.nao = parent.nao
        self._rho = None
        self._C = None

    # ------------------------------------------------------------------ #
    #  Hamiltonien bloc-diagonal + orbitales                             #
    # ------------------------------------------------------------------ #
    def block_hamiltonian(self):
        """H du parent, avec zéros partout hors des blocs."""
        H = np.zeros((self.nao, self.nao))
        for name, idx, nelec in self.blocks:
            sub = self.parent.H[np.ix_(idx, idx)]
            H[np.ix_(idx, idx)] = sub
        return H

    def orbitals(self):
        """Diagonalisation bloc par bloc -> (C, eps, occ).
        Retourne coefficients, energies et occupations PAR ORBITALE.
        Dégénérescences -> remplissage fractionnaire (½ e⁻ par OM dégénérée,
        ou plus généralement nelec/n_degen)."""
        if self._C is not None:
            return self._C, self._eps, self._occ

        C = np.zeros((self.nao, self.nao))
        eps = np.zeros(self.nao)
        occ = np.zeros(self.nao)
        col = 0

        for name, idx, nelec in self.blocks:
            n = len(idx)
            sub = self.parent.H[np.ix_(idx, idx)]
            e, v = np.linalg.eigh(sub)
            for k in range(n):
                C[idx, col + k] = v[:, k]
                eps[col + k] = e[k]

            # --- remplissage du bloc, avec gestion des dégénérescences --- #
            # électrons restants, on descend les OM du bloc
            remaining = float(nelec)
            k = 0
            while remaining > 1e-9 and k < n:
                # dégénérescence : combien d'OM ont (presque) la même énergie ?
                deg = 1
                while k + deg < n and abs(e[k + deg] - e[k]) < 1e-8:
                    deg += 1
                # électrons à caser dans ce niveau dégénéré
                place = min(remaining, 2.0 * deg)
                occ[col + k: col + k + deg] = place / deg  # ½ e⁻ si 2e⁻ / 2 OM, etc.
                remaining -= place
                k += deg
            col += n

        self._C, self._eps, self._occ = C, eps, occ
        return C, eps, occ

    # ------------------------------------------------------------------ #
    #  Matrice densité : LE PIVOT                                         #
    # ------------------------------------------------------------------ #
    def density_matrix(self):
        """rho = sum_k  occ_k * |C_k><C_k|   (base AO orthonormée Hückel).
        Invariante par rotation des OM dégénérées -> bien définie unique."""
        if self._rho is None:
            C, eps, occ = self.orbitals()
            self._rho = (C * occ) @ C.T
        return self._rho

    def n_electrons(self):
        """Contrôle : doit redonner le nombre d'électrons du bloc."""
        return np.trace(self.density_matrix())   # AO orthonormées -> Tr(rho) = ne

    # ------------------------------------------------------------------ #
    #  Observables usuelles                                               #
    # ------------------------------------------------------------------ #
    def pi_charges(self):
        """q_mu = npi_mu - (rho)_mu_mu   (S = identité en Hückel)."""
        rho = self.density_matrix()
        return self.parent.npi - np.diag(rho)

    def bond_orders(self):
        """Indice de liaison pi: P_mu_nu = (rho)_mu_nu, mu != nu."""
        rho = self.density_matrix()
        return rho - np.diag(np.diag(rho))

    # ------------------------------------------------------------------ #
    #  Recouvrement entre structures : indice de Carbó                   #
    # ------------------------------------------------------------------ #
    def carbo(self, other):
        """Similitude normée entre les ensembles d'électrons de deux structures.
        S_IJ = Tr(rho_I rho_J) / sqrt(Tr(rho_I²) Tr(rho_J²))   dans [0, 1].
        Marche avec occupations fractionnaires et OM non orthogonales."""
        rI, rJ = self.density_matrix(), other.density_matrix()
        num = np.trace(rI @ rJ)
        den = np.sqrt(np.trace(rI @ rI) * np.trace(rJ @ rJ))
        return num / den if den > 1e-12 else 0.0

    def overlap_matrix_element(self, other):
        """Version brute (non normée): Tr(rho_I rho_J). Utile pour diag/viz."""
        return np.trace(self.density_matrix() @ other.density_matrix())


class LewisSpace:
    """Ensemble des structures d'une molécule parente."""

    def __init__(self, parent, structures):
        self.parent = parent
        self.structures = structures

    def similarity_matrix(self):
        n = len(self.structures)
        S = np.eye(n)   # Carbo(I,I) = 1
        for i in range(n):
            for j in range(i + 1, n):
                S[i, j] = S[j, i] = self.structures[i].carbo(self.structures[j])
        return S

    def weights(self):
        """Poids bruts des structures par somme des similitudes aux autres
        (premier jet : plus la structure ressemble à l'ensemble, plus elle
        est 'centrale'). À remplacer plus tard par la vraie pondération."""
        S = self.similarity_matrix()
        w = S.sum(axis=1)
        return w / w.sum()

if __name__ == "__main__":
    import numpy as np

    class FakeParent:
        nao = 4
        # butadiene-ish : Hückel simple
        H = np.array([[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]], float)
        npi = np.ones(4)

    p = FakeParent()
    # structure : deux liaisons localisées C1=C2 et C3=C4, 2 e- chacune
    s1 = LocalizedStructure(p, [("double", [0,1], 2), ("double", [2,3], 2)])
    # structure : une seule double "délocalisée" sur les 4
    s2 = LocalizedStructure(p, [("double", [1,2], 2),("radical", [0], 1),("radical", [3], 1)])

    assert abs(s1.n_electrons() - 4) < 1e-9, "perte d'électrons !"
    assert abs(s2.n_electrons() - 4) < 1e-9
    print("carbo s1-s2 :", s1.carbo(s2))
    print("charges s1 :", s1.pi_charges())
    print("indices liaison s1 :\n", s1.bond_orders())
