import numpy as np
import copy
from dataclasses import dataclass
from typing import Optional, List, Tuple
from qe_input import QEInput, Atom, AtomicSpecies

@dataclass
class DefectSpec:
    defect_type: str          # "vacancy", "substitution", "interstitial"
    site_index: Optional[int] = None       # atom index in primitive cell (0-based)
    site_species: Optional[str] = None     # species label to select site (e.g. "S", "Mo")
    new_species: str = ""     # species label for substitution/interstitial
    new_mass: float = 0.0     # mass for new species
    new_pseudo: str = ""      # pseudopotential for new species
    position: Optional[np.ndarray] = None  # for interstitial: fractional coords in supercell

    def resolve_site_index(self, prim: QEInput) -> int:
        if self.site_species is not None:
            for i, atom in enumerate(prim.atoms):
                if atom.label == self.site_species:
                    return i
            raise ValueError(
                f"No atom with species '{self.site_species}' found in primitive cell. "
                f"Available species: {[a.label for a in prim.atoms]}"
            )
        if self.site_index is not None:
            if self.site_index < 0 or self.site_index >= prim.nat:
                raise ValueError(
                    f"site_index {self.site_index} out of range [0, {prim.nat})"
                )
            return self.site_index
        raise ValueError("Either --site or --site-species must be provided.")


def build_supercell_centered(prim: QEInput, nx: int, ny: int, nz: int,
                             defect: DefectSpec) -> Tuple[QEInput, QEInput, QEInput, np.ndarray]:
    n = np.array([nx, ny, nz], dtype=float)
    defect_site_idx = defect.resolve_site_index(prim)

    # Step 1: Build unshifted supercell
    sc = copy.deepcopy(prim)

    if sc.cell_parameters is not None:
        for i in range(3):
            sc.cell_parameters[i] *= n[i]

    sc.atoms = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = np.array([ix, iy, iz], dtype=float)
                for atom in prim.atoms:
                    if prim.positions_units != 'crystal':
                        raise ValueError(
                            "Supercell generation requires ATOMIC_POSITIONS "
                            "in crystal coordinates. Got: " + prim.positions_units
                        )
                    new_pos = (atom.position + shift) / n
                    sc.atoms.append(Atom(label=atom.label, position=new_pos))
    sc.nat = len(sc.atoms)

    if sc.kpoints_mode == 'automatic' and sc.kpoints_grid:
        sc.kpoints_grid = [
            max(1, sc.kpoints_grid[0] // nx),
            max(1, sc.kpoints_grid[1] // ny),
            max(1, sc.kpoints_grid[2] // nz),
            0, 0, 0
        ]

    if sc.nbnd:
        sc.nbnd = sc.nbnd * nx * ny * nz

    # Step 2: Compute centering shift
    defect_pos_super = sc.atoms[defect_site_idx].position.copy()

    target = np.array([0.5, 0.5, 0.5])
    if nz == 1:
        target[2] = defect_pos_super[2]

    centering_shift = target - defect_pos_super

    # Step 3: Apply centering shift to build pristine supercell
    pristine = copy.deepcopy(sc)
    for atom in pristine.atoms:
        atom.position = (atom.position + centering_shift) % 1.0

    # Step 4: Build defect supercell (same shift, then introduce defect)
    defect_sc = copy.deepcopy(pristine)

    if defect.defect_type == "vacancy":
        removed = defect_sc.atoms.pop(defect_site_idx)
        defect_sc.nat -= 1
        print(f"Vacancy: removed {removed.label} at {removed.position}")

    elif defect.defect_type == "substitution":
        if not defect.new_species:
            raise ValueError("Substitution requires --new-species.")
        old_label = defect_sc.atoms[defect_site_idx].label
        defect_sc.atoms[defect_site_idx].label = defect.new_species
        print(f"Substitution: {old_label} -> {defect.new_species} "
              f"at {defect_sc.atoms[defect_site_idx].position}")
        if not any(sp.label == defect.new_species for sp in defect_sc.species):
            defect_sc.species.append(AtomicSpecies(
                label=defect.new_species,
                mass=defect.new_mass,
                pseudo_file=defect.new_pseudo
            ))

    elif defect.defect_type == "interstitial":
        if defect.position is None:
            raise ValueError("Interstitial requires --interstitial-pos.")
        if not defect.new_species:
            raise ValueError("Interstitial requires --new-species.")
        defect_sc.atoms.append(Atom(
            label=defect.new_species, position=defect.position
        ))
        defect_sc.nat += 1
        if not any(sp.label == defect.new_species for sp in defect_sc.species):
            defect_sc.species.append(AtomicSpecies(
                label=defect.new_species,
                mass=defect.new_mass,
                pseudo_file=defect.new_pseudo
            ))

    # Clean up species list
    used_labels = set(a.label for a in defect_sc.atoms)
    defect_sc.species = [sp for sp in defect_sc.species if sp.label in used_labels]
    defect_sc.ntyp = len(defect_sc.species)

    # Step 5: Build aligned primitive cell
    prim_shift = centering_shift * n

    aligned_prim = copy.deepcopy(prim)
    for atom in aligned_prim.atoms:
        atom.position = (atom.position + prim_shift) % 1.0

    print(f"Centering shift (supercell frac): {centering_shift}")
    print(f"Centering shift (primitive frac): {prim_shift}")
    print(f"Defect now at supercell center:   "
          f"{(defect_pos_super + centering_shift) % 1.0}")

    return pristine, defect_sc, aligned_prim, centering_shift


def generate_nscf_kpoints(prim_kgrid: List[int]) -> List[tuple]:
    nk1, nk2, nk3 = prim_kgrid[0], prim_kgrid[1], prim_kgrid[2]
    kpoints = []
    for i in range(nk1):
        for j in range(nk2):
            for k in range(nk3):
                kx = float(i) / nk1
                ky = float(j) / nk2
                kz = float(k) / nk3
                kpoints.append((kx, ky, kz, 1.0))
    return kpoints
