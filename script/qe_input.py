from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import numpy as np

@dataclass
class AtomicSpecies:
    label: str        # e.g. "Mo"
    mass: float       # atomic mass
    pseudo_file: str  # pseudopotential filename

@dataclass
class Atom:
    label: str                  # species label
    position: np.ndarray        # 3-vector
    if_pos: Optional[List[int]] = None  # constraint flags

@dataclass
class QEInput:
    """Full representation of a pw.x input file."""

    # &CONTROL
    calculation: str = "scf"
    prefix: str = "pwscf"
    outdir: str = "./"
    pseudo_dir: str = "./"
    verbosity: str = "high"
    restart_mode: str = "from_scratch"
    tstress: bool = True
    tprnfor: bool = True
    control_extra: dict = field(default_factory=dict)

    # &SYSTEM
    ibrav: int = 0
    A: Optional[float] = None         # lattice constant in Angstrom
    celldm: Optional[List[float]] = None
    nat: int = 0
    ntyp: int = 0
    ecutwfc: float = 60.0
    ecutrho: Optional[float] = None
    noncolin: bool = False
    lspinorb: bool = False
    nbnd: Optional[int] = None
    nspin: Optional[int] = None
    system_extra: dict = field(default_factory=dict)

    # &ELECTRONS
    diagonalization: str = "david"
    mixing_mode: str = "plain"
    mixing_beta: float = 0.7
    conv_thr: float = 1.0e-6
    electrons_extra: dict = field(default_factory=dict)

    # Cards
    species: List[AtomicSpecies] = field(default_factory=list)
    atoms: List[Atom] = field(default_factory=list)
    cell_parameters: Optional[np.ndarray] = None   # (3,3) lattice vectors
    cell_units: str = "alat"                        # alat | bohr | angstrom
    positions_units: str = "crystal"                # crystal | angstrom | alat | bohr
    kpoints_mode: str = "automatic"                 # automatic | crystal | gamma
    kpoints_grid: Optional[List[int]] = None        # [nk1, nk2, nk3, sk1, sk2, sk3]
    kpoints_list: Optional[List[Tuple[float,float,float,float]]] = None


def parse_qe_input(filepath: str) -> QEInput:
    """Parse a pw.x input file into a QEInput object."""
    inp = QEInput()
    with open(filepath, 'r') as f:
        text = f.read()

    # --- Parse namelists ---
    lines = text.split('\n')
    current_namelist = None
    namelist_content = {}

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('!'):
            continue

        # Detect namelist start
        if stripped.startswith('&'):
            current_namelist = stripped[1:].lower().split()[0]
            namelist_content[current_namelist] = {}
            continue

        # Detect namelist end
        if stripped.startswith('/'):
            current_namelist = None
            continue

        # Inside a namelist: parse key=value pairs
        if current_namelist:
            for assignment in stripped.split(','):
                assignment = assignment.strip()
                if '=' in assignment:
                    key, val = assignment.split('=', 1)
                    key = key.strip().lower()
                    val = val.strip().strip("'").strip('"')
                    namelist_content[current_namelist][key] = val

    # Map namelist values onto QEInput fields
    ctrl = namelist_content.get('control', {})
    inp.calculation = ctrl.get('calculation', 'scf')
    inp.prefix = ctrl.get('prefix', 'pwscf')
    inp.outdir = ctrl.get('outdir', './')
    inp.pseudo_dir = ctrl.get('pseudo_dir', './')
    inp.verbosity = ctrl.get('verbosity', 'low')
    inp.restart_mode = ctrl.get('restart_mode', 'from_scratch')
    inp.tstress = _parse_bool(ctrl.get('tstress', '.true.'))
    inp.tprnfor = _parse_bool(ctrl.get('tprnfor', '.true.'))

    # Store unrecognized control keys
    _known_control = {'calculation', 'prefix', 'outdir', 'pseudo_dir',
                      'verbosity', 'restart_mode', 'tstress', 'tprnfor'}
    for k, v in ctrl.items():
        if k not in _known_control:
            inp.control_extra[k] = v

    sys_ = namelist_content.get('system', {})
    inp.ibrav = int(sys_.get('ibrav', '0'))
    inp.A = _parse_float(sys_['a']) if 'a' in sys_ else None
    inp.nat = int(sys_.get('nat', '0'))
    inp.ntyp = int(sys_.get('ntyp', '0'))
    inp.ecutwfc = _parse_float(sys_.get('ecutwfc', '60.0'))
    if 'ecutrho' in sys_:
        inp.ecutrho = _parse_float(sys_['ecutrho'])
    inp.noncolin = _parse_bool(sys_.get('noncolin', '.false.'))
    inp.lspinorb = _parse_bool(sys_.get('lspinorb', '.false.'))
    if 'nbnd' in sys_:
        inp.nbnd = int(sys_['nbnd'])
    if 'nspin' in sys_:
        inp.nspin = int(sys_['nspin'])

    # Store unrecognized system keys
    _known_system = {'ibrav', 'a', 'nat', 'ntyp', 'ecutwfc', 'ecutrho',
                     'noncolin', 'lspinorb', 'nbnd', 'nspin'}
    for k, v in sys_.items():
        if k not in _known_system:
            inp.system_extra[k] = v

    elec = namelist_content.get('electrons', {})
    inp.diagonalization = elec.get('diagonalization', 'david')
    inp.mixing_mode = elec.get('mixing_mode', 'plain')
    inp.mixing_beta = _parse_float(elec.get('mixing_beta', '0.7'))
    inp.conv_thr = _parse_float(elec.get('conv_thr', '1.0e-6'))

    # Store unrecognized electrons keys
    _known_electrons = {'diagonalization', 'mixing_mode', 'mixing_beta', 'conv_thr'}
    for k, v in elec.items():
        if k not in _known_electrons:
            inp.electrons_extra[k] = v

    # --- Parse cards ---
    _parse_atomic_species(lines, inp)
    _parse_cell_parameters(lines, inp)
    _parse_atomic_positions(lines, inp)
    _parse_kpoints(lines, inp)

    return inp


def _parse_float(val: str) -> float:
    return float(val.lower().replace('d', 'e'))


def _parse_bool(val: str) -> bool:
    return val.lower().strip('.') in ('true', 't', '1', 'yes')


def _parse_atomic_species(lines, inp):
    """Parse ATOMIC_SPECIES card."""
    idx = _find_card(lines, 'ATOMIC_SPECIES')
    if idx is None:
        return
    for i in range(idx + 1, len(lines)):
        parts = lines[i].split()
        if len(parts) < 3 or _is_card_or_namelist(lines[i]):
            break
        inp.species.append(AtomicSpecies(
            label=parts[0], mass=float(parts[1]), pseudo_file=parts[2]
        ))


def _parse_cell_parameters(lines, inp):
    """Parse CELL_PARAMETERS card."""
    idx = _find_card(lines, 'CELL_PARAMETERS')
    if idx is None:
        return
    header = lines[idx].strip()
    if '{' in header:
        inp.cell_units = header.split('{')[1].split('}')[0].strip()
    elif '(' in header:
        inp.cell_units = header.split('(')[1].split(')')[0].strip()
    else:
        parts = header.split()
        if len(parts) > 1:
            inp.cell_units = parts[1].lower()

    vecs = []
    for i in range(idx + 1, idx + 4):
        vecs.append([float(x) for x in lines[i].split()])
    inp.cell_parameters = np.array(vecs)


def _parse_atomic_positions(lines, inp):
    """Parse ATOMIC_POSITIONS card."""
    idx = _find_card(lines, 'ATOMIC_POSITIONS')
    if idx is None:
        return
    header = lines[idx].strip()
    for token in ['crystal', 'angstrom', 'bohr', 'alat']:
        if token in header.lower():
            inp.positions_units = token
            break

    for i in range(idx + 1, len(lines)):
        parts = lines[i].split()
        if len(parts) < 4 or _is_card_or_namelist(lines[i]):
            break
        pos = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
        if_pos = None
        if len(parts) >= 7:
            if_pos = [int(parts[4]), int(parts[5]), int(parts[6])]
        inp.atoms.append(Atom(label=parts[0], position=pos, if_pos=if_pos))

    inp.nat = len(inp.atoms)


def _parse_kpoints(lines, inp):
    """Parse K_POINTS card."""
    idx = _find_card(lines, 'K_POINTS')
    if idx is None:
        return
    header = lines[idx].strip().lower()
    if 'automatic' in header:
        inp.kpoints_mode = 'automatic'
        parts = lines[idx + 1].split()
        inp.kpoints_grid = [int(x) for x in parts]
    elif 'gamma' in header:
        inp.kpoints_mode = 'gamma'
    elif 'crystal' in header:
        inp.kpoints_mode = 'crystal'
        nk = int(lines[idx + 1].split()[0])
        inp.kpoints_list = []
        for i in range(idx + 2, idx + 2 + nk):
            parts = lines[i].split()
            inp.kpoints_list.append(tuple(float(x) for x in parts[:4]))


def _find_card(lines, card_name):
    for i, line in enumerate(lines):
        if line.strip().upper().startswith(card_name.upper()):
            return i
    return None


def _is_card_or_namelist(line):
    s = line.strip().upper()
    cards = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS',
             'CELL_PARAMETERS', 'CONSTRAINTS', 'OCCUPATIONS',
             'ATOMIC_FORCES', 'ADDITIONAL_K_POINTS']
    return s.startswith('&') or s.startswith('/') or any(s.startswith(c) for c in cards)


def write_qe_input(inp: QEInput, filepath: str):
    """Write a QEInput object to a pw.x input file."""
    with open(filepath, 'w') as f:
        # &CONTROL
        f.write("&control\n")
        f.write(f"   calculation = '{inp.calculation}'\n")
        f.write(f"   verbosity='{inp.verbosity}'\n")
        f.write(f"   restart_mode='{inp.restart_mode}',\n")
        f.write(f"   prefix='{inp.prefix}',\n")
        if inp.tstress:
            f.write("   tstress = .true.\n")
        if inp.tprnfor:
            f.write("   tprnfor = .true.\n")
        f.write(f"   pseudo_dir = '{inp.pseudo_dir}',\n")
        f.write(f"   outdir='{inp.outdir}'\n")
        for k, v in inp.control_extra.items():
            f.write(f"   {k} = {v}\n")
        f.write("/\n")

        # &SYSTEM
        f.write("&system\n")
        f.write(f"   ibrav=  {inp.ibrav},\n")
        if inp.A is not None:
            f.write(f"   A={inp.A}\n")
        f.write(f"   nat=  {inp.nat}, ntyp= {inp.ntyp},\n")
        f.write(f"   ecutwfc ={inp.ecutwfc},\n")
        if inp.ecutrho:
            f.write(f"   ecutrho ={inp.ecutrho},\n")
        if inp.noncolin:
            f.write("   noncolin=.true.\n")
        if inp.lspinorb:
            f.write("   lspinorb=.true.\n")
        if inp.nbnd:
            f.write(f"   nbnd={inp.nbnd}\n")
        for k, v in inp.system_extra.items():
            f.write(f"   {k} = {v}\n")
        f.write("/\n\n")

        # &ELECTRONS
        f.write("&electrons\n")
        f.write(f"   diagonalization='{inp.diagonalization}'\n")
        f.write(f"   mixing_mode = '{inp.mixing_mode}'\n")
        f.write(f"   mixing_beta = {inp.mixing_beta}\n")
        f.write(f"   conv_thr =  {inp.conv_thr:.1e}\n")
        for k, v in inp.electrons_extra.items():
            f.write(f"   {k} = {v}\n")
        f.write("/\n\n\n")

        # ATOMIC_SPECIES
        f.write("ATOMIC_SPECIES\n")
        for sp in inp.species:
            f.write(f"{sp.label} {sp.mass} {sp.pseudo_file}\n")
        f.write("\n\n")

        # CELL_PARAMETERS
        if inp.cell_parameters is not None and inp.ibrav == 0:
            f.write(f"CELL_PARAMETERS {inp.cell_units}\n")
            for vec in inp.cell_parameters:
                f.write(f" {vec[0]:16.8f} {vec[1]:16.8f} {vec[2]:16.8f}\n")
            f.write("\n")

        # ATOMIC_POSITIONS
        f.write(f"ATOMIC_POSITIONS {inp.positions_units}\n")
        for atom in inp.atoms:
            f.write(f"{atom.label:14s} {atom.position[0]:18.10f}"
                    f" {atom.position[1]:18.10f} {atom.position[2]:18.10f}\n")
        f.write("\n\n")

        # K_POINTS
        if inp.kpoints_mode == 'automatic':
            f.write("K_POINTS automatic\n")
            g = inp.kpoints_grid
            f.write(f"{g[0]} {g[1]} {g[2]} {g[3]} {g[4]} {g[5]}\n")
        elif inp.kpoints_mode == 'gamma':
            f.write("K_POINTS gamma\n")
        elif inp.kpoints_mode == 'crystal':
            f.write("K_POINTS crystal\n")
            f.write(f"{len(inp.kpoints_list)}\n")
            for kp in inp.kpoints_list:
                f.write(f"  {kp[0]:.12f}   {kp[1]:.12f}   "
                        f"{kp[2]:.12f}  {kp[3]:.1f}\n")
