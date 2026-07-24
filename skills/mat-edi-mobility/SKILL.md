---
name: mat-edi-mobility
description: Compute defect-limited carrier mobility and electron-defect scattering matrix elements in 2D and 3D semiconductors from first principles with Quantum ESPRESSO and the EDI plugin.
category: materials
---

# mat-edi-mobility

## Goal
To compute the point-defect-limited carrier mobility $\mu(T)$ of a semiconductor
and the underlying electron-defect scattering matrix elements
$M(\mathbf{k}_i, \mathbf{k}_f) = \langle \psi_{\mathbf{k}_i} | \Delta V | \psi_{\mathbf{k}_f} \rangle$
from first principles, using the supercell difference-potential method of the
Quantum ESPRESSO + EDI pipeline. Here $\Delta V = V_{\text{defect}} - V_{\text{pristine}}$
is the Kohn-Sham perturbation potential of an isolated point defect, extracted as
the difference between the potentials of a defect-containing supercell and a
pristine supercell. The matrix elements are computed on a coarse k-grid,
transformed to a maximally localized Wannier basis $M(\mathbf{R}, \mathbf{R}')$,
and interpolated onto dense fine grids. The state-resolved scattering rate follows
from Fermi's golden rule,

$$\frac{1}{\tau_{n\mathbf{k}}} = \frac{2\pi}{\hbar}\, n_{\text{d}}\, \frac{1}{N_{\mathbf{k}}} \sum_{m,\mathbf{k}'} |M_{n\mathbf{k}, m\mathbf{k}'}|^2\, \delta(\varepsilon_{n\mathbf{k}} - \varepsilon_{m\mathbf{k}'}),$$

and the mobility from the linearized Boltzmann transport equation in the
self-energy (SERTA) and momentum (MRTA) relaxation-time approximations. The three
first-class outputs are the **direct** matrix elements $M$ on the coarse grid, the
**Wannier-interpolated** $M$ on a fine k-path, and the final mobility $\mu(T)$.

## Background
Defect-limited transport requires $M(\mathbf{k}_i, \mathbf{k}_f)$ on k-grids far
denser than any tractable supercell calculation can supply. EDI solves this the
way EPW solves the electron-phonon problem: compute $M$ directly on a coarse grid
commensurate with the supercell, rotate into a Wannier basis where
$M(\mathbf{R}, \mathbf{R}')$ decays rapidly with $|\mathbf{R}|$ and $|\mathbf{R}'|$,
and interpolate to $300\times300$ or denser at negligible cost. The scattering
rate scales linearly with the defect concentration $n_{\text{d}}$ (dilute,
incoherent-defect limit), so **the mobility scales as $1/n_{\text{d}}$ and the
concentration must always be reported alongside $\mu$**. For 2D materials the
pristine and defect supercells are computed independently, so their electrostatic
potentials carry an arbitrary offset that must be removed by vacuum-level
alignment (`pot_align = 'vacuum'`).

This skill computes the *defect*-limited mobility; the *phonon*-limited channel
is a separate calculation (e.g. via QE + EPW + Wannier90). At finite temperature
the two channels combine by Matthiessen's rule,
$1/\mu_{\text{total}} = 1/\mu_{\text{phonon}} + 1/\mu_{\text{defect}}$.
Defect *energetics* (formation energies, charge states), which set which defects
dominate at a given growth condition, require separate defect-thermodynamics
calculations and are not covered here.

### External binary: building Quantum ESPRESSO 7.5 + EDI
EDI is an **in-tree plugin** of Quantum ESPRESSO, not a standalone code: it links
QE's compiled static libraries and must be built against a matching QE version.

1. Obtain QE 7.5 (GitLab tag archive; no ReleasePack is required —
   `make w90` auto-fetches the pinned Wannier90 submodule even from a tarball) and
   configure it.
2. Clone EDI into the QE root and build:

   ```bash
   # Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
   cd <QE-root>
   make pw
   make w90
   git clone https://github.com/yuanyue-liu-group/EDI edi-code
   cd edi-code && make
   ```

   This produces `src/edi.x` (MPI) and `src/extract_pot.x` (serial).
3. **gfortran patch (required for GCC builds):** `ed_coarse.f90` contains lines up
   to 158 characters, so add long-line support to `edi-code/src/makefile` after the
   `include` line:

   ```make
   F90FLAGS += -ffree-line-length-none
   FFLAGS   += -ffree-line-length-none
   ```

4. **`extract_pot.x` segfault patch (required):** `extract_pot.f90` calls
   `clean_pw(.TRUE.)` at line 77, before any `read_file`, which segfaults on a
   fresh process. Comment it out (or delete that line) and rebuild:

   ```fortran
   ! CALL clean_pw(.TRUE.)   ! upstream bug: runs before any read_file, segfaults
   ```

Pseudopotentials for the worked example are PseudoDojo NC-SR v0.5 PBE **stringent**
(`https://www.pseudo-dojo.org/pseudos/nc-sr-05_pbe_stringent_upf.tgz`); the Mo and
S UPFs are used as `Mo.upf` and `S.upf`. QE 7.5 must be built with Wannier90.

## Instructions

> [!IMPORTANT]
> Quantum ESPRESSO (`pw.x`) and EDI (`edi.x`, `extract_pot.x`) are external
> compiled MPI binaries. This skill targets a source build of
> **QE 7.5 + EDI v2.0** (repo `main` @ `41fed72`) with **PseudoDojo NC-SR v0.5 PBE stringent**
> pseudopotentials. The `# Env: base-agent` annotations apply only to the Python
> parser scripts. Reference input decks for every step are in
> [`resources/inputs/`](resources/inputs/). Replace `<ranks>` with your MPI rank
> count; `mpirun -np` and `srun -n` are interchangeable. Set `pseudo_dir` in each
> deck to your UPF directory (the decks use `./pseudo/`).

The pipeline is a DAG. The two supercell SCFs (step 3) are independent of the
primitive cell and of each other and can run in parallel. Steps 5 and 6 must run
in order (6 reads the Wannier archive written by 5).

### 1. Primitive-cell SCF
Compute the ground-state charge density of the pristine primitive cell.
`assume_isolated = '2D'` truncates the out-of-plane Coulomb tail and **must** be
consistent across every calculation. Deck:
[`resources/inputs/scf_primitive.in`](resources/inputs/scf_primitive.in).

```bash
# Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
mpirun -np <ranks> pw.x -nk <pools> < scf_primitive.in > scf_primitive.out
```

Confirm `convergence has been achieved`.

### 2. Primitive-cell NSCF on the explicit coarse grid
Produce the Bloch states EDI folds into Wannier functions. The grid **must be an
explicit `K_POINTS crystal` list** (never `automatic`) and its dimensions **must
equal `coarse_nk1/nk2/nk3`** in the EDI input (12x12x1 here). Generate the list
with the helper, then splice it into the NSCF deck
([`resources/inputs/nscf_primitive.in`](resources/inputs/nscf_primitive.in)
already contains the 144-point 12x12x1 card):

```bash
# Env: base-agent
python .agents/skills/mat-edi-mobility/scripts/gen_kgrid.py 12 12 1
```

```bash
# Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
mpirun -np <ranks> pw.x -nk <pools> < nscf_primitive.in > nscf_primitive.out
```

`nbnd` (17 here) must cover the disentanglement window (the 5 kept bands 13-17).
Note: the upstream example ships an 18x18x1 NSCF that is inconsistent with its
`coarse_nk1 = 12`; this deck fixes that to 12x12x1.

### 3. Pristine and defect supercell SCFs
Compute the Kohn-Sham potentials whose difference is $\Delta V$. Both are
Gamma-only SCFs of a 6x6x1 supercell (108 atoms pristine; 107 with one S vacancy).
The supercell lattice vectors **must be exact integer multiples of the primitive
cell** and host-atom positions must fold back to the primitive sites, or $\Delta V$
is corrupted. Relax the defect supercell before this final SCF. Decks:
[`resources/inputs/scf_pristine_super.in`](resources/inputs/scf_pristine_super.in),
[`resources/inputs/scf_defect_super.in`](resources/inputs/scf_defect_super.in).

```bash
# Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
mpirun -np <ranks> pw.x < scf_pristine_super.in > scf_pristine_super.out
mpirun -np <ranks> pw.x < scf_defect_super.in > scf_defect_super.out
```

These are Gamma-only runs (`K_POINTS gamma`, a single k-point): do **not** add
`-nk` here — there is nothing to pool, unlike the primitive-cell runs.

### 4. Extract the difference potential
Run the serial `extract_pot.x` to write the pristine and defect local KS
potentials as cube files (`V_p.cube`, `V_d.cube`) on the supercell real-space
grid. Deck: [`resources/inputs/extract_pot.in`](resources/inputs/extract_pot.in).

```bash
# Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
mpirun -np 1 extract_pot.x < extract_pot.in > extract_pot.out
```

`extract_pot.x` is serial; run it with a single rank.

### 5. First EDI pass: direct matrix elements + Wannierization
Compute $M(\mathbf{k}_i, \mathbf{k}_f)$ **directly** on the coarse 12x12x1 grid,
build $M(\mathbf{R}, \mathbf{R}')$, run the Wannier90 minimization, and write the
Wannier archive `mos2_edmatw_2d.bin` that the transport pass reuses. Deck:
[`resources/inputs/edi_setup.in`](resources/inputs/edi_setup.in)
(`edwread = .false.`, `wannierize = .true.`, `do_transport = .false.`).

```bash
# Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
mpirun -np <ranks> edi.x -nk <ranks> -i edi_setup.in > edi_setup.out
```

Key namelist variables:
- `nbndsub = 5`, `proj(1) = 'Mo:d'`, `bands_skipped = 'exclude_bands = 1-12'` — the
  5-band Mo-d Wannier manifold; **material-specific**.
- `dis_win_min/max`, `dis_froz_min/max` — outer and frozen disentanglement windows
  (eV); set them around the transport bands.
- `wdata(1) = 'guiding_centres = .true.'` — required in long-vacuum 2D cells or the
  Wannier centres fold to a wrong z-image.
- `coarse_nk*` must equal the step-2 NSCF grid; `fine_nk*` sets the interpolation grid.
- `edmat_interp_from_file = .true.` with `filki_interp`/`filkf_interp` also writes
  the interpolated $M$ along the validation path (`mos2_edmat_interp.dat`), used in 6b.

Because it also computes the coarse-grid Bloch $M$, this pass writes the direct
matrix elements to `mos2_edmat_bloch.dat` (columns: `iki ikf  kix kiy kiz  kfx kfy
kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)  |M_loc|^2  |M_nl|^2`).

`-nk <ranks>` (pools = rank count) is recommended for EDI performance.

### 6. Second EDI pass: fine-grid interpolation + BTE transport
Interpolate $M$ onto the fine grid and solve the BTE for $\mu(T)$. This pass reuses
the Wannier archive (`edwread = .true.`, `wannierize = .false.`,
`do_transport = .true.`). Deck:
[`resources/inputs/edi_transport.in`](resources/inputs/edi_transport.in).

```bash
# Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
mpirun -np <ranks> edi.x -nk <ranks> -i edi_transport.in > edi_transport.out
```

Key transport variables:
- `fine_nk1 = fine_nk2 = 300` — production grid; `48` is a fast smoke test.
- `transport_win_min/max` — narrow energy window (eV) around the carrier band edge;
  only pocket states scatter. Anchor it to the band edges printed by the **SCF**
  run (`highest occupied, lowest unoccupied level`); the NSCF uses
  `calculation = 'bands'` and prints no HOMO/LUMO or Fermi summary.
- `carrier_conc` (cm^-2 in 2D) sets the Fermi level by bisection; `defect_conc`
  (cm^-2 in 2D) scales the rate. **Mobility scales as $1/$`defect_conc` — report it.**
- `delta_method = 'gaussian'`, `delta_sigma = 0.01` — energy-conserving delta; also
  `'triangular'` (2D optimized) and `'adaptive'` (velocity-dependent smearing).
- `temps(i)`, `nstemp` — temperature list.

Outputs: `mos2_transport.dat` (`T  mu_SERTA_xx  mu_MRTA_xx  mu_SERTA_yy  mu_MRTA_yy`
in cm^2/Vs) and `mos2_inv_tau.dat` (state-resolved
`ik ibnd E inv_tau_SERTA inv_tau_MRTA tau_SERTA tau_MRTA`). Parse them:

```bash
# Env: base-agent
python .agents/skills/mat-edi-mobility/scripts/parse_transport.py mos2_transport.dat --output-dir results/
python .agents/skills/mat-edi-mobility/scripts/parse_inv_tau.py mos2_inv_tau.dat --output-dir results/
```

Fine-grid convergence is cheap to sweep: because `edwread = .true.` reuses
`mos2_edmatw_2d.bin`, each rerun redoes only the interpolation + BTE (seconds to a
few minutes, e.g. ~6 s at 48x48 up to ~3 min at 300x300 for this example). The
outputs are prefix-named (`mos2_transport.dat`), so copy them to a per-grid name
between reruns. Converge `fine_nk*` before trusting $\mu$: a 48x48 smoke grid runs
21-23% low here.

### 6b. Validate: direct vs Wannier-interpolated M on a k-path
Confirm the Wannier interpolation reproduces the directly computed $M$ (this is the
interpolation's correctness gate, analogous to the DFPT-vs-EPW $|g|$ benchmark in
the electron-phonon problem). Both are evaluated on the SAME
k-path: one initial point at K ([`resources/inputs/ki.dat`](resources/inputs/ki.dat))
and a 300-point Gamma-K-M-Gamma final path
([`resources/inputs/kf.dat`](resources/inputs/kf.dat)). Deck:
[`resources/inputs/edi_direct.in`](resources/inputs/edi_direct.in), which enables
**both** `edmat_interp_from_file` and `edmat_direct_from_file` on the same path.

**Direct mode needs NSCF wavefunctions at every path k-point** — the coarse
12x12x1 NSCF only contains 8 of the 300 path points. Run it in two steps:

1. Launch the direct pass once. It detects the missing k-points, **auto-writes
   `nscf_custom.in`** listing exactly the ki + kf points it needs, and aborts.
2. Run that NSCF into a **separate** save tree (e.g. `primitive_path/dout/`; ~1.5 min
   on 32 tasks for the 301-point path here), point `edi_outdir` at it (already set to
   `../primitive_path/dout/` in the deck), and rerun the direct pass (~3m40s at
   32 ranks x 4 cores for this example — direct mode is memory-heavier; see
   Constraints).

```bash
# Requires: Quantum ESPRESSO 7.5 + EDI (external MPI build)
# 1. First attempt writes nscf_custom.in and stops:
mpirun -np <ranks> edi.x -nk <ranks> -i edi_direct.in > edi_direct.out
# 2. Run the path NSCF into primitive_path/, then rerun the direct pass:
mpirun -np <ranks> pw.x < nscf_custom.in > nscf_custom.out    # stage into primitive_path/dout/
mpirun -np <ranks> edi.x -nk <ranks> -i edi_direct.in > edi_direct.out
```

This writes `mos2_edmat_direct.dat` (direct $M$, **absolute** NSCF band indices) and
`mos2_edmat_interp.dat` (interpolated $M$, Wannier-subspace indices 1-5). Compare:

```bash
# Env: base-agent
python .agents/skills/mat-edi-mobility/scripts/parse_edmat.py mos2_edmat_interp.dat --output-dir results/
python .agents/skills/mat-edi-mobility/scripts/compare_edmat.py mos2_edmat_direct.dat mos2_edmat_interp.dat --band-offset 12 --band-sum --output-dir results/
```

`--band-offset 12` maps the two files' band conventions (direct = interp + number
of excluded bands). `--band-sum` compares the **gauge-invariant band-summed**
$\sum_{mn} |M_{mn}|^2$ per k-point — the meaningful metric, since individual band
pairs mix under the arbitrary Wannier gauge along band crossings (drop the flag
for the pairwise statistics). A small relative RMS confirms the
interpolation. Finally, benchmark $\mu(T)$ against the upstream reference with
`compare_reference.py` (see the example).

## Examples
See [`examples/mos2-s-vacancy/`](examples/mos2-s-vacancy/README.md) for the full
end-to-end run on a sulfur vacancy in monolayer MoS2, including the coarse-grid
direct $M$, the Wannier-interpolated $M$ with the direct-vs-interpolated
comparison, the mobility-vs-temperature table validated against the upstream
reference outputs, and a literature cross-check against ACS Nano **18**, 8511 (2024).

## Constraints
- **External codes**: requires an MPI build of Quantum ESPRESSO 7.5 with Wannier90
  and the EDI plugin built in-tree (`edi.x`, `extract_pot.x`); the Python parsers
  run in `base-agent`. GCC builds require `-ffree-line-length-none` (see Background).
  Version pins reflect what was targeted, not a claim that other versions are broken.
- **In-tree build only**: EDI links QE 7.5 static libraries; the QE version must
  match exactly. Our QE 7.4.1 build is not reusable.
- **Point defects only**: EDI v2.0 handles point defects (vacancies, substitutions,
  interstitials). Extended defects (surfaces, grain boundaries) are not yet supported.
- **SERTA / MRTA, not iterative BTE**: mobility is relaxation-time-based; there is no
  self-consistent iterative BTE solution.
- **Mobility scales with defect concentration**: $\mu \propto 1/$`defect_conc` in the
  dilute limit — always report the concentration alongside $\mu$.
- **2D systems** need `assume_isolated = '2D'` (vacuum-separated cell, c = 24 A here) in
  every `pw.x` run and `pot_align = 'vacuum'` in EDI, consistently.
- **Coarse-grid match**: the primitive NSCF grid must equal `coarse_nk1/nk2/nk3`, and
  the NSCF must use an explicit `K_POINTS crystal` list, or the Wannier folding is
  inconsistent.
- **Supercell consistency**: supercell lattice vectors must be exact integer multiples
  of the primitive cell, host atoms must fold to primitive sites, and the FFT grid must
  be commensurate, or $\Delta V$ is corrupted.
- **Memory scales with supercell volume x ranks-per-node**: in the setup, transport,
  and direct modes every MPI rank holds the **full** supercell real-space potentials
  ($V_d$, $V_p$, $V_{\text{colin}}$) — about 2 GB/rank for this 6x6 MoS2 example
  (240x240x300 FFT grid), and direct mode holds an extra copy. On ~2 GB/core clusters,
  under-subscribe the node: run setup/transport at half ranks
  (`--ntasks=64 --cpus-per-task=2`, `edi.x -nk 64`) and direct mode at a quarter
  (`--ntasks=32 --cpus-per-task=4`, `edi.x -nk 32`). Symptom of getting it wrong: the
  job stays `RUNNING` with the output frozen right after the `Full double-FT` banner
  while ranks are silently OOM-killed and survivors deadlock in the next MPI collective.
  CPU burn does not prove progress — check output-file mtime and live task count.
- **Material-specific inputs**: `ecutwfc`, `nbnd`, `nbndsub`, projections,
  `exclude_bands`, and the disentanglement/transport windows are set for MoS2 and must
  be adapted to your material and band manifold.

## References
- Z. Xiao, R. Guo, C. Zhang, Y. Liu, "Point Defect Limited Carrier Mobility in 2D
  Transition Metal Dichalcogenides", *ACS Nano* **18**, 8511 (2024).
  [DOI](https://doi.org/10.1021/acsnano.4c01033)
- I.-T. Lu, J.-J. Zhou, M. Bernardi, "Efficient ab initio calculations of electron-defect
  scattering and defect-limited carrier mobility", *Phys. Rev. Materials* **3**, 033804
  (2019). [DOI](https://doi.org/10.1103/PhysRevMaterials.3.033804)
- I.-T. Lu, J. Park, J.-J. Zhou, M. Bernardi, "Ab initio electron-defect interactions using
  Wannier functions", *npj Comput. Mater.* **6**, 17 (2020).
  [DOI](https://doi.org/10.1038/s41524-020-0284-y)
- A. A. Mostofi, et al., "An updated version of wannier90: A tool for obtaining
  maximally-localised Wannier functions", *Comput. Phys. Commun.* **185**, 2309 (2014).
  [DOI](https://doi.org/10.1016/j.cpc.2014.05.003)
- P. Giannozzi, et al., "Quantum ESPRESSO toward the exascale", *J. Chem. Phys.* **152**,
  154105 (2020). [DOI](https://doi.org/10.1063/5.0005082)
- The EDI code paper is in preparation; cite the repository:
  [github.com/yuanyue-liu-group/EDI](https://github.com/yuanyue-liu-group/EDI).
