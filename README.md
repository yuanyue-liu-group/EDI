<p align="center">
  <img src="figs/edi.png" alt="EDI workflow" width="600">
</p>

# EDI: Electron-Defect Interaction

A Quantum ESPRESSO plugin for computing electron-defect scattering matrix elements and related quantities (e.g. carrier mobility) from first principles.

EDI uses the supercell approach to extract defect perturbation potentials, Wannier interpolation to obtain matrix elements on arbitrarily fine k-grids, and the Boltzmann transport equation (BTE) to compute defect-limited carrier mobility in both 2D and 3D materials.

## Features of EDI version 2.0

- **Electron-defect matrix elements** via the supercell difference-potential method: `M(k_i, k_f) = <psi_ki | V_defect - V_pristine | psi_kf>`, including both local and nonlocal (Kleinman-Bylander) contributions
- **Wannierization of matrix elements** from Bloch to Wannier representation `M(R, R')` for efficient interpolation to dense k-grids
- **Wannier90 integration** for interpolation of band structure (energies, velocities) and matrix elements 
- **SERTA and MRTA mobility** with irreducible Brillouin zone (IBZ) symmetry exploitation
- **Multiple delta function methods**: triangular (2D), fixed Gaussian, and adaptive Gaussian (velocity-dependent smearing)
- **Spin-orbit coupling (SOC)** support: two-component spinor wavefunctions, `dvan_so` nonlocal pseudopotentials, automatic detection from QE save files
- **MPI parallelization** with pool-based (ki, kf) pair distribution and panel broadcast for memory-efficient wavefunction sharing
- **Automated setup**: Python scripts generate all supercell and primitive cell inputs from a single `scf.in`
- **One-shot potential extraction** by extract_pot.x

Previous version: [EDI-1.0](https://github.com/yuanyue-liu-group/EDI_old)

## Workflow
<p align="center">
  <img src="figs/workflow.png" alt="EDI workflow" width="800">
</p>


## Prerequisites

- [Quantum ESPRESSO 7.5](https://www.quantum-espresso.org/) compiled with [Wannier90](https://wannier.org/)
- MPI-enabled Fortran compiler (tested with AOCC 3.1.0, GCC, Intel)
- FFTW3 library
- Python 3 with NumPy (for setup scripts only)

## Installation

If you are using Anvil, before compilation, the modules can be set up as:

```bash
module reset
module load aocc/3.1.0 openmpi/4.1.6
module load amdblis/3.0 amdlibflame/3.0 amdlibm/3.0 fftw

```


EDI is built as a plugin inside the QE source tree. Place this repository at `<QE-root>/edi-code/`:

```bash
cd <QE-root>
git clone git@github.com:yuanyue-liu-group/EDI.git edi-code
cd edi-code
make
```

This compiles QE's `pw.x` (if not already built) and then builds `edi.x` in `src/`.

To build only `edi.x` (assuming QE libraries are already compiled):

```bash
cd src
make
```
## Usage

A typical EDI workflow proceeds in 4 steps:
- Step 1: Prepare inputs for DFT calculation by QE
- Step 2: Run DFT calculations (pw.x)
- Step 3: Extract potentials and Run EDI (extract_pot.x and edi.x)
- Step 4: Post-processing


### Step 1: Prepare inputs for DFT calculation by QE 
Electron-defect matrix elements calculations need unperturbed wavefunctions from primitive cell and Kohn-Sham potentials from pristine and defective supercells. 

**Output directory structure:**
 
| File | Description | Key dependency |
|------|-------------|----------------|
| `primitive/scf.in` | Primitive cell SCF | Direct copy of the source `scf.in` |
| `primitive/nscf.in` | Primitive cell NSCF on a coarse k-grid | Same cell & pseudopotentials as `primitive/scf.in`; k-grid matches the supercell dimensions (e.g., 6×6×1) |
| `pristine_super/scf.in` | Pristine supercell SCF (no defect) | Lattice = integer multiple of primitive cell; FFT grid commensurate with primitive cell |
| `defect_super/scf.in` | Defect supercell SCF | Same lattice & FFT grid as `pristine_super/scf.in`; one site modified (vacancy / substitution / interstitial) |
| `edi/edi.in` | EDI main input | References output directories of all DFT calculations above |
| `run_all.sh` | Full pipeline job script | — |

EDI provides python scripts to help user to generate the input files. **Users can also modify the input to customize their defect, but must take care of the follow conditions**

**Important: Structural Consistency for Custom Defect Inputs**

Users who employ external tools to generate defect supercell input files must ensure full consistency between the supercell and the primitive cell. Specifically, the following conditions must be satisfied:

- **Lattice parameters.** The supercells' lattice vectors must be exact integer multiples of the primitive cell lattice vectors. Any discrepancy will corrupt the difference potential.

- **Atomic coordinates.** The positions of host atoms **before relaxation** must correspond precisely to the equivalent crystallographic sites in the pristine primitive cell upon folding back. A different origin convention, coordinate rounding, or independent relaxation of the pristine supercell will introduce spurious contributions to the defect perturbation potential.

- **FFT grid.** The FFT grid dimensions used in the supercell SCF calculation must be commensurate with those of the primitive cell. For example, if the primitive cell uses an FFT grid of $(N_1, N_2, N_3)$ and the supercells are constructed with scaling factors $(S_1, S_2, S_3)$, the supercells' FFT grid must be set to $(S_1 N_1, S_2 N_2, S_3 N_3)$. A mismatched grid will lead to inconsistent real-space sampling and erroneous difference potentials. 

Failure to satisfy any of these conditions will yield unreliable electron-defect matrix elements. Users are strongly recommended to verify consistency by comparing the lattice parameters, atomic coordinates and FFT grid of the supercells against primitive cell before proceeding with the EDI calculation.

Below is the usage of EDI python scripts:

Use `gen_supercell.py` to automatically generate all QE input files from a primitive cell `scf.in`:

```bash
# S vacancy in a 6x6x1 MoS2 supercell
python script/gen_supercell.py \
    --input scf.in \
    --nx 6 --ny 6 --nz 1 \
    --defect vacancy --site-species S \
    --nprocs 72 --output ./edi_run/
```

This generates the following directory structure:

```
edi_run/
  primitive/scf.in          # Primitive cell SCF
  primitive/nscf.in         # Primitive cell NSCF on coarse k-grid
  pristine_super/scf.in     # Pristine supercell SCF
  defect_super/scf.in       # Defect supercell SCF
  edi/edi.in                # EDI input file
  run_all.sh                # Full pipeline job script
```

Example defect types:

```bash
# Vacancy (remove one atom)
python script/gen_supercell.py --input scf.in --nx 6 --ny 6 --nz 1 \
    --defect vacancy --site-species S

# Substitution (replace S with O)
python script/gen_supercell.py --input scf.in --nx 6 --ny 6 --nz 1 \
    --defect substitution --site-species S \
    --new-species O --new-mass 15.999 --new-pseudo O_ONCV_PBE_FR-1.0.upf

# Interstitial (insert atom at specified position)
python script/gen_supercell.py --input scf.in --nx 3 --ny 3 --nz 3 \
    --defect interstitial \
    --new-species Li --new-mass 6.941 --new-pseudo Li.upf \
    --interstitial-pos 0.5 0.5 0.5
```



### Step 2: Run DFT calculations

```bash
# Primitive cell
cd primitive
mpirun -np 72 pw.x -nk 8 < scf.in  > scf.out
mpirun -np 72 pw.x -nk 8 < nscf.in > nscf.out

# Supercells
cd ../pristine_super
mpirun -np 72 pw.x < scf.in > scf.out
cd ../defect_super
mpirun -np 72 pw.x < scf.in > scf.out
```

### Step 3: Extract potentials and Run EDI

```bash
cd edi
srun -n 1 extract_pot.x -i extract_pot.in > extract_pot.out
srun -n 72 edi.x -nk 72 -i edi.in > edi.out
```

The `-nk` flag sets the number of k-point pools for MPI parallelization. For transport calculations, use `-nk` equal to the total number of MPI ranks for best performance.

### Step 4: Output files

EDI produces:

| File | Description |
|------|-------------|
| `prefix_transport.dat` | Mobility vs temperature (SERTA and MRTA, xx and yy components) |
| `prefix_inv_tau.dat` | State-resolved inverse lifetimes 1/tau(n,k) and tau(n,k) |
| `prefix_edmatw_2d.bin` | Wannier-basis matrix elements M(R, R'), reusable with `edwread = .true.` |

## Input Reference

EDI reads a single Fortran namelist `&edinput_nml` from the input file.

### General Settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `edi_prefix` | string | `'pwscf'` | Prefix of primitive cell QE calculation |
| `edi_outdir` | string | `'./'` | Output directory containing primitive cell QE save data |

### Supercell Potentials

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `potfile_d` | string | | Cube file for defect supercell potential (from `extract_pot.x`) |
| `potfile_p` | string | | Cube file for pristine supercell potential (from `extract_pot.x`) |
| `pot_align` | string | `'vacuum'` | Potential alignment: `'vacuum'` (2D), `'core'`, or `'none'` |
| `defect_center(3)` | real | 0,0,0 | Defect center in supercell fractional coords (for `pot_align='core'`) |
| `core_align_radius` | real | 2.0 | Averaging sphere radius in Bohr (for `pot_align='core'`) |

### Matrix Element Calculation

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `do_edmat` | logical | `.false.` | Compute matrix elements from supercell potentials |
| `pristine_prefix` | string | | Prefix of pristine supercell calculation (for on-the-fly loading) |
| `pristine_outdir` | string | | Output directory of pristine supercell (for on-the-fly loading) |
| `defect_prefix` | string | | Prefix of defect supercell calculation (for on-the-fly loading) |
| `defect_outdir` | string | | Output directory of defect supercell (for on-the-fly loading) |
| `edwwrite` | logical | `.true.` | Write Wannier-basis M(R) to binary file after computation |
| `edwread` | logical | `.false.` | Read M(R) from binary file, skipping recomputation |
| `band_ed` | string | | Band range for direct calculation mode (e.g., `'13-17'`) |

### Wannier90 Settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wannierize` | logical | `.false.` | Run Wannier90 minimization (`.false.` reuses existing `.chk`) |
| `nbndsub` | integer | 0 | Number of Wannier functions |
| `dis_win_min` / `dis_win_max` | real | +/-9999 | Disentanglement outer window (eV) |
| `dis_froz_min` / `dis_froz_max` | real | +/-9999 | Disentanglement frozen window (eV) |
| `num_iter` | integer | 200 | Maximum Wannier90 iterations |
| `proj(i)` | string | | Projection strings (e.g., `'Mo:d'`, `'S:p'`) |
| `wdata(i)` | string | | Additional Wannier90 `.win` file entries |
| `bands_skipped` | string | | Bands to exclude (e.g., `'exclude_bands = 1-12'`) |
| `auto_projections` | logical | `.true.` | Use automatic projections (SCDM) |

### K-grids

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `coarse_nk1` / `nk2` / `nk3` | integer | 0 | Coarse k-grid dimensions (must match the primitive cell NSCF) |
| `fine_nk1` / `nk2` / `nk3` | integer | 0 | Fine k-grid for Wannier interpolation and transport |

### Transport Settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `do_transport` | logical | `.false.` | Enable transport calculation |
| `transport_win_min` | real | -9999 | Energy window lower bound (eV) |
| `transport_win_max` | real | 9999 | Energy window upper bound (eV) |
| `carrier_conc` | real | 1e13 | Carrier concentration (cm^-2 for 2D, cm^-3 for 3D) |
| `defect_conc` | real | 1e12 | Defect concentration (cm^-2 for 2D, cm^-3 for 3D) |
| `nstemp` | integer | 1 | Number of temperature points |
| `temps(i)` | real | 300.0 | Temperature values in Kelvin |
| `delta_method` | string | `'triangular'` | Delta function: `'triangular'`, `'gaussian'`, or `'adaptive'` |
| `delta_sigma` | real | 0.1 | Fixed Gaussian broadening (eV), used when `delta_method = 'gaussian'` |

### Custom K-point Modes

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `edmat_direct_from_file` | logical | `.false.` | Compute M(ki,kf) directly for k-points listed in files |
| `filki_direct` / `filkf_direct` | string | | Input k-point files for direct mode (crystal coordinates) |
| `edmat_interp_from_file` | logical | `.false.` | Interpolate M(ki,kf) for k-points listed in files |
| `filki_interp` / `filkf_interp` | string | | Input k-point files for interpolation mode |

## Example Input

### extract_pot.in (potential extraction, run with 1 MPI rank)

```fortran
&extract_pot_nml
  defect_prefix   = 'mos2_defect'
  defect_outdir   = '../defect_super/dout/'
  pristine_prefix = 'mos2_pristine'
  pristine_outdir = '../pristine_super/dout/'
  potfile_d = 'V_d.cube'
  potfile_p = 'V_p.cube'
/
```

### edi.in (full calculation)

```fortran
&edinput_nml
  ! Primitive cell
  edi_prefix = 'mos2'
  edi_outdir = '../primitive/dout/'

  ! Pre-extracted supercell potentials (from extract_pot.x)
  potfile_d = 'V_d.cube'
  potfile_p = 'V_p.cube'
  pot_align = 'vacuum'

  ! Matrix elements
  do_edmat = .true.

  ! Wannierization
  wannierize = .true.
  nbndsub = 5
  dis_win_min = -7.0
  dis_win_max = -0.8
  dis_froz_min = -7.0
  dis_froz_max = -3.4
  num_iter = 1000
  auto_projections = .false.
  proj(1) = 'Mo:d'
  bands_skipped = 'exclude_bands = 1-12'

  ! K-grids
  coarse_nk1 = 18
  coarse_nk2 = 18
  coarse_nk3 = 1
  fine_nk1 = 300
  fine_nk2 = 300
  fine_nk3 = 1

  ! Transport
  do_transport = .true.
  transport_win_min = -4.25
  transport_win_max = -4.10
  carrier_conc = 1.0e10
  defect_conc = 1.0e12
  delta_method = 'gaussian'
  delta_sigma = 0.01
  nstemp = 3
  temps(1) = 100.0
  temps(2) = 200.0
  temps(3) = 300.0
/
```

### Reuse previously computed Wannier matrix elements

```fortran
&edinput_nml
  edi_prefix = 'mos2'
  edi_outdir = '../primitive/dout/'
  potfile_d = 'V_d.cube'
  potfile_p = 'V_p.cube'
  edwread = .true.         ! read M(R) from file
  wannierize = .false.     ! reuse existing Wannier90 checkpoint
  nbndsub = 5
  coarse_nk1 = 18, coarse_nk2 = 18, coarse_nk3 = 1
  fine_nk1 = 300, fine_nk2 = 300, fine_nk3 = 1
  do_transport = .true.
  transport_win_min = -4.25, transport_win_max = -4.10
  carrier_conc = 1.0e10, defect_conc = 1.0e12
  delta_method = 'gaussian', delta_sigma = 0.01
  nstemp = 1, temps(1) = 300.0
/
```

## Source Code Structure

All Fortran source files are in `src/` (26 files):

### Core

| File | Description |
|------|-------------|
| `edi.f90` | Main program: reads input, orchestrates Wannierization, matrix elements, transport |
| `edi_input.f90` | Namelist parsing (`&edinput_nml`), parameter broadcasting |
| `edic_mod.f90` | Shared data types (`V_file` for supercell potentials), workspace arrays |

### Potential Extraction

| File | Description |
|------|-------------|
| `extract_pot.f90` | Standalone program (`extract_pot.x`): reads QE save files and writes supercell KS potentials to cube files with full double precision. Run serially (1 MPI rank) before `edi.x`. |

### Matrix Elements

| File | Description |
|------|-------------|
| `ed_coarse.f90` | Core computation: supercell potential folding, spinor FFT, double Fourier transform, panel broadcast MPI, direct/interpolation modes, cube file reader |
| `get_betavkb.f90` | Kleinman-Bylander beta projectors for nonlocal pseudopotentials |
| `get_vloc_onthefly.f90` | On-the-fly extraction of supercell potentials (V_loc + V_xc) from QE save files (legacy, used when `potfile_d/p` not set) |

### Wannier Interpolation

| File | Description |
|------|-------------|
| `wannierize.f90` | Wannier90 driver: spread minimization, checkpoint I/O |
| `edbloch2wan.f90` | Bloch-to-Wannier transformation of electron-defect matrix elements |
| `wan2bloch.f90` | Wannier-to-Bloch interpolation of Hamiltonian and matrix elements on fine grids |
| `edi_read_hr.f90` | Reader for Wannier90 `_hr.dat` (Hamiltonian in Wannier basis) |
| `wigner_seitz.f90` | Wigner-Seitz supercell construction for real-space lattice vectors and degeneracies |
| `wann_common.f90` | Shared Wannier90 data: U matrices, lattice vectors, band mapping |

### Transport

| File | Description |
|------|-------------|
| `transport.f90` | Boltzmann transport: IBZ symmetry, SERTA/MRTA scattering rates, Fermi level bisection, mobility with spin degeneracy |
| `delta_weights.f90` | Energy-conserving delta function: triangular (2D), Gaussian (fixed/adaptive) |
| `bz_symmetry.f90` | Irreducible BZ construction and symmetry mapping |
| `ep_constants.f90` | Physical constants (Rydberg, Bohr, etc.) |

### Quantum ESPRESSO Interface

| File | Description |
|------|-------------|
| `edi_pw2wan.f90` | High-level interface between QE wavefunctions and Wannier90 |
| `pw2wan_edi.f90` | Low-level pw2wannier90 routines: AMN, MMN, eigenvalue extraction |
| `io_edi.f90` | Wavefunction I/O utilities |
| `io_var_edi.f90` | I/O variable declarations |
| `input_edi.f90` | Input module data structures |
| `kfold_edi.f90` | K-point folding and mapping for Wannier90 |
| `openfilepw_edi.f90` | QE file opening utilities |
| `low_lvl_edi.f90` | Low-level utility routines |
| `parallelism_edi.f90` | MPI pool distribution utilities |
| `global_var.f90` | Global band indices and kept-band mapping |

### Helper Scripts (`script/`)

| File | Description |
|------|-------------|
| `gen_supercell.py` | CLI tool: generates all DFT and EDI input files from a primitive cell `scf.in` |
| `supercell.py` | Library: supercell construction, defect insertion (vacancy/substitution/interstitial), coordinate transforms |
| `qe_input.py` | Library: QE input file parser and writer |
| `gen_pp_inputs.py` | Library: EDI input file and SLURM job script generator |

## Theory

EDI implements the ab initio electron-defect scattering framework [[2]](#ref2)[[3]](#ref3)[[4]](#ref4) with modifications. The central idea is to treat a point defect as a perturbation to the pristine crystal, extract the perturbation potential from supercell DFT calculations, and then evaluate scattering matrix elements between Bloch states of the pristine primitive cell. Below we summarize the key equations.


### Electron-Defect Matrix Element

The electron-defect matrix element between Bloch states $|\psi_{n\mathbf{k}}\rangle$ and $|\psi_{m\mathbf{k}'}\rangle$ of the pristine crystal is defined as:

$$M_{n \mathbf{k}, m \mathbf{k'}} = \langle \psi_{n\mathbf{k}} | \Delta V | \psi_{m\mathbf{k'}} \rangle$$

where $\Delta V$ is the defect perturbation potential, i.e., the difference between the Kohn-Sham potentials of the defect-containing and pristine supercells. Because the defect is spatially localized, $\Delta V(\mathbf{r})$ decays to zero far from the defect site, and a sufficiently large supercell captures the full perturbation.

The total Kohn-Sham potential includes local and nonlocal parts. Accordingly, the matrix element decomposes as:

$$M_{n\mathbf{k},m\mathbf{k}'} = M^{\mathrm{loc}} + M^{\mathrm{NL,defect}} - M^{\mathrm{NL,pristine}}$$

The local part involves the real-space integral of $\Delta V_{\mathrm{loc}}(\mathbf{r})$, which contains the local ionic pseudopotential, Hartree, and exchange-correlation potentials. For calculations with spin-orbit coupling (SOC), the wavefunctions are two-component spinors, and the local contribution sums over spinor components $\sigma \in \lbrace\uparrow,\downarrow\rbrace$:

$$M^{\mathrm{loc}} = \sum_{\sigma} \int \psi^{\ast}_{n\mathbf{k},\sigma}(\mathbf{r}) \Delta V_{\mathrm{loc}}(\mathbf{r}) \psi_{m\mathbf{k}',\sigma}(\mathbf{r}) d\mathbf{r}$$

In practice, the wavefunctions are expressed in a plane-wave basis and the integral is efficiently computed via FFT by folding the supercell $\Delta V$ onto the primitive-cell reciprocal grid.

The nonlocal part arises from the Kleinman-Bylander (KB) separable pseudopotentials. For each atom $I$ with projectors $|\beta^I_i\rangle$ and coupling coefficients $D^I_{ij}$, the nonlocal matrix element takes the form:

$$M^{\mathrm{NL}} = \sum_{I,i,j} \langle \psi_{n\mathbf{k}} | \beta^I_i \rangle D^I_{ij} \langle \beta^I_j | \psi_{m\mathbf{k}'} \rangle$$

The defect perturbation to the nonlocal potential is obtained by subtracting the pristine-supercell contribution from the defect-supercell contribution. Under SOC, the coupling coefficients become spin-dependent matrices $D^{I,\sigma\sigma'\!}_{ij}$ acting on the spinor indices.

Because the supercell DFT calculations for the pristine and defect systems are performed independently, an arbitrary constant offset between the two electrostatic potentials must be removed. EDI supports vacuum-level alignment (for 2D systems, where the potential plateaus in the vacuum region) and core-level alignment (averaging the potential on a sphere far from the defect center) [[2]](#ref2).


### Wannier Interpolation

Direct evaluation of $M_{n\mathbf{k},m\mathbf{k}'}$ on the dense k-grids required for converged transport calculations is computationally prohibitive, since each matrix element involves plane-wave wavefunctions and FFTs. EDI overcomes this bottleneck by Wannier interpolation [[4]](#ref4), leveraging maximally localized Wannier functions (MLWFs) constructed using Wannier90 [[5]](#ref5).

The matrix elements are first computed on a coarse $N_{\mathbf{k}}$-point grid commensurate with the supercell and then transformed to the Wannier gauge:

$$M_{ij}(\mathbf{R}, \mathbf{R}') = \frac{1}{N_{\mathbf{k}}^2} \sum_{\mathbf{k},\mathbf{k}'} e^{+i\mathbf{k}\cdot\mathbf{R}} e^{-i\mathbf{k}'\cdot\mathbf{R}'} \sum_{n,m} U^{\ast}_{ni}(\mathbf{k}) M^{(\mathrm{B})}_{nm}(\mathbf{k},\mathbf{k}') U_{mj}(\mathbf{k}')$$

where $n, m$ are band indices, $i, j$ are Wannier function indices, and $U_{ni}(\mathbf{k})$ are the unitary rotation matrices from Wannier90 (containing both the disentanglement and gauge-optimization transformations). $\mathbf{R}$, $\mathbf{R}'$ are real-space lattice vectors. Because the defect potential and the Wannier functions are both spatially localized, $M_{ij}(\mathbf{R}, \mathbf{R}')$ decays rapidly with $|\mathbf{R}|$ and $|\mathbf{R}'|$, which is the key property enabling accurate interpolation.

From the Wannier-basis matrix elements, one reconstructs $M$ at arbitrary fine k-points $(\mathbf{k}, \mathbf{k}')$ in two steps. First, a Fourier transform back to reciprocal space yields the matrix in the Wannier gauge:

$$\tilde{M}_{ij}(\mathbf{k}, \mathbf{k}') = \sum_{\mathbf{R},\mathbf{R}'} e^{-i\mathbf{k}\cdot\mathbf{R}} e^{+i\mathbf{k}'\cdot\mathbf{R}'} M_{ij}(\mathbf{R}, \mathbf{R}')$$

Then, a unitary rotation transforms back to the Bloch band basis:

$$M_{nm}(\mathbf{k}, \mathbf{k}') = \sum_{i,j} U_{ni}(\mathbf{k}) \tilde{M}_{ij}(\mathbf{k}, \mathbf{k}') U^{\ast}_{mj}(\mathbf{k}')$$

The band structure (eigenvalues and group velocities) is simultaneously interpolated from the Wannier-basis Hamiltonian $H(\mathbf{R})$ read from Wannier90's `_hr.dat` file, following the standard procedure. This two-step approach (coarse-grid calculation + Wannier interpolation) allows EDI to reach fine k-grids of $300\times300$ or denser at negligible additional cost.

<p align="center">
  <img src="figs/wannier_interp.png" alt="Wannier interpolation validation" width="500">
</p>
<p align="center"><em>Validation of Wannier-interpolated matrix elements (lines) against direct DFT calculation (dots) along the high-symmetry path Γ–K–M–Γ in monolayer MoS₂.</em></p>


### Scattering Rate

Given the matrix elements, the defect-limited scattering rate for a Bloch state $|n\mathbf{k}\rangle$ is obtained from Fermi's golden rule:

$$\frac{1}{\tau_{n\mathbf{k}}} = \frac{2\pi}{\hbar} n_{\mathrm{d}} \frac{1}{N_{\mathbf{k}}} \sum_{m,\mathbf{k}'} |M_{n\mathbf{k},m\mathbf{k}'}|^2 \delta(\varepsilon_{n\mathbf{k}} - \varepsilon_{m\mathbf{k}'})$$

where $n_{\mathrm{d}}$ is the defect concentration per unit cell and the energy-conserving delta function enforces elastic scattering. The factor $n_{\mathrm{d}}$ reflects the dilute-defect limit, where scattering events from different defect sites are treated as independent and incoherent.

For numerical evaluation of the delta function, EDI provides three methods: (i) a triangular (linear tetrahedron-like) method optimized for 2D Brillouin zones, (ii) a fixed Gaussian broadening with user-specified width $\sigma$, and (iii) an adaptive Gaussian scheme where $\sigma$ is set proportional to the energy variation $|\nabla_\mathbf{k}\varepsilon \cdot \Delta\mathbf{k}|$ across each k-mesh cell, following the approach used in EPW for electron-phonon problems.


### Carrier Mobility

The defect-limited carrier mobility tensor is computed from the linearized Boltzmann transport equation. Under the self-energy relaxation time approximation (SERTA), the mobility reads:

$$\mu_{\alpha\beta} = \frac{g_s e}{n_c} \frac{1}{N_{\mathbf{k}}} \sum_{n\mathbf{k}} v_{n\mathbf{k},\alpha} v_{n\mathbf{k},\beta} \tau_{n\mathbf{k}} \Bigl(-\frac{\partial f}{\partial \varepsilon}\Bigr)_{\varepsilon_{n\mathbf{k}}}$$

where $v_{n\mathbf{k},\alpha} = \hbar^{-1} \partial \varepsilon_{n\mathbf{k}}/\partial k_\alpha$ is the band velocity (obtained from the Wannier-interpolated Hamiltonian), $f(\varepsilon)$ is the Fermi-Dirac distribution, $n_c$ is the carrier concentration, and $g_s$ is the spin degeneracy factor ($g_s = 2$ for collinear spin-degenerate calculations, $g_s = 1$ when SOC is included). The Fermi level $\varepsilon_F$ at each temperature is determined self-consistently by bisection so that the integral of $f(\varepsilon)$ over the band structure reproduces the specified carrier concentration.

EDI also implements the momentum relaxation time approximation (MRTA), which incorporates the scattering-angle dependence through a $(1 - \cos\theta)$ factor in the scattering rate:

$$\frac{1}{\tau^{\mathrm{MRTA}}_{n\mathbf{k}}} = \frac{2\pi}{\hbar} n_{\mathrm{d}} \frac{1}{N_{\mathbf{k}}} \sum_{m,\mathbf{k}'} |M_{n\mathbf{k},m\mathbf{k}'}|^2 \delta(\varepsilon_{n\mathbf{k}} - \varepsilon_{m\mathbf{k}'}) \Bigl(1 - cos\theta\Bigr)$$

The MRTA provides an improved approximation to the full iterative BTE solution for elastic scattering processes, as forward-scattering events (small momentum transfer) contribute less to resistivity.

All transport quantities are computed on the irreducible Brillouin zone (IBZ) and symmetry-expanded, which reduces the computational cost by a factor equal to the order of the point group.

## References

<a id="ref1"></a> [1] X. XX et al., "EDI: Electron-defect interaction and defect-limited transport from first principles," (in preparation)

<a id="ref2"></a> [2] Z. Xiao et al., "Point Defect Limited Carrier Mobility in 2D Transition Metal Dichalcogenides," [ACS Nano 18 (11), 8511-8516 (2024)](https://doi.org/10.1021/acsnano.4c01033)

<a id="ref3"></a> [3] I. Lu et al., "Efficient ab initio calculations of electron-defect scattering and defect-limited carrier mobility," [Phys. Rev. Materials 3, 033804 (2019)](https://doi.org/10.1103/PhysRevMaterials.3.033804)

<a id="ref4"></a> [4] I. Lu et al., "Ab initio electron-defect interactions using Wannier functions," [npj Comput Mater 6, 17 (2020)](https://doi.org/10.1038/s41524-020-0284-y)

<a id="ref5"></a> [5] A. A. Mostofi et al., "An updated version of wannier90: A tool for obtaining maximally-localised Wannier functions," [Comput. Phys. Commun. 185, 2309 (2014)](https://doi.org/10.1016/j.cpc.2014.05.003)

<a id="ref6"></a> [6] P. Giannozzi et al., "Quantum ESPRESSO toward the exascale," [J. Chem. Phys. 152, 154105 (2020)](https://doi.org/10.1063/5.0005082)

## License

This software is distributed under the GNU General Public License v3.0. See the QE license for details.
