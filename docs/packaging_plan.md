# Plan: Package ABSCO as a pixi project + installable Python distribution

## Context

ABSCO is currently a flat "clone-and-run-in-place" repository: top-level scripts
(`run_LBLRTM_ABSCO.py`, `ABSCO_preprocess.py`, `ABSCO_compute.py`, `build_models.py`, etc.),
a `common/` git submodule loaded with `sys.path.append('common')`, Fortran submodules
(`LBLRTM/`, `LNFL/`), small bundled data (`PT_grid/`, `VMR/`, `FSCDXS_line_params.csv`),
and a large Zenodo-downloaded line file (`AER_Line_File/`). Configuration is a hand-edited
`ABSCO_config.ini` full of **absolute machine-specific paths** injected by `build_models.py`.

The goal is to make ABSCO an **installed tool** that a user runs from an arbitrary working
directory, without editing paths by hand:

1. A **pixi-controlled project** for reproducible dev/build/runtime environments.
2. An **installable Python package** (src-layout) publishable to PyPI, shipping **prebuilt
   Fortran binaries** in platform wheels.
3. Bundled data files usable from the installed package **without the user specifying paths**,
   overridable only when using custom files.
4. All temporary/intermediate files written under the user's **current working directory** by default.
5. A **separate config-assistant CLI** that generates the `ABSCO_config.ini` consumed by the
   generation command. The user typically supplies only `wn1`, `wn2`, `outres`, `molnames`;
   the tool **suggests `lblres`** automatically.

### Decisions locked in (from interview)
- **Fortran**: not a pure-Python package — **bundle prebuilt binaries** in platform wheels (cibuildwheel).
- **Runtime artifacts** (line file, executables): live in a **user data dir** (platformdirs / `$ABSCO_DATA_DIR`), fetched/staged once, reused across working dirs.
- **CLI names**: package `absco`; commands **`absco-config`** (assistant), **`absco-generate`** (generation, replaces `run_LBLRTM_ABSCO.py`), and **`absco-init`** (end-user artifact initialization — see below).
- **Artifact initialization**: a user who installs a prebuilt binary wheel must still initialize runtime artifacts (line file, and staging binaries into the data dir). Provide a dedicated **`absco-init`** command for this, separate from the dev/build-time compile step.
- **Deps**: modernize to Python ≥3.10 with lower-bound floors; pixi lock pins exact versions.
- **Config UX**: flags with interactive fallback for missing required values.
- **Refactor**: full **src-layout** package; replace `sys.path`/`__file__`/`getcwd` path tricks with `importlib.resources` + explicit working-dir handling.
- **`common/`**: **vendor** the needed modules into `absco._common`.

### Key verified fact
`ABSCO_preprocess.py:151,183` sets `kernwidth = outres/lblres` and **hard-exits unless the
ratio is an exact power of 2** (`np.log2(kernwidth) % 1 == 0`). The notebook
`notebooks/create_absco_config-tropomi.ipynb` derives `lblres` correctly:
`power = ceil(log2(outres / 1.5e-4)); lblres = outres / 2**power` — targets `lblres ≤ 1.5e-4`
while guaranteeing a power-of-2 ratio. **This logic is consistent with the code and is what
`absco-config` will implement** (including the notebook's round-trip check that the ratio
still holds after string formatting into the `.ini`).

---

## Target structure

```
pyproject.toml            # build-system, deps, entry points, package-data
pixi.toml                 # (or [tool.pixi] in pyproject) envs + tasks
docs/packaging_plan.md    # this plan, for review
src/absco/
  __init__.py
  preprocess.py           # from ABSCO_preprocess.py
  compute.py              # from ABSCO_compute.py
  build_models.py         # dev/build-time Fortran build + Zenodo fetch (repackaged)
  paths.py                # NEW: data-dir + resource resolution (importlib.resources, platformdirs)
  config_assistant.py     # NEW: lblres suggestion + .ini writer
  cli/
    __init__.py
    generate.py           # entry point absco-generate  (was run_LBLRTM_ABSCO.py)
    config.py             # entry point absco-config
    init.py               # entry point absco-init (end-user artifact setup)
    build.py              # entry point absco-build (dev-only Fortran compile)
    read_tables.py        # was read_ABSCO_tables.py
    split_config.py / join_tables.py
  _common/                # vendored from common/ submodule (utils.py, lblTools.py, RC_utils.py, FortranFile.py)
  data/
    PT_grid/*             # AIRS_P_air.txt, build_temp_array.txt, dump_air_pres
    VMR/*.csv             # standard profiles + HDO example
    FSCDXS_line_params.csv
    ABSCO_config.template.ini   # NEW: template with placeholders, no absolute paths
```

---

## Work items

### 1. Package skeleton & build system
- Create **src-layout** package `src/absco/`; move the modules above. Fix imports
  (`import ABSCO_preprocess` → `from absco import preprocess`; drop every `sys.path.append('common')`).
- `pyproject.toml`:
  - `[build-system]` using **scikit-build-core** or **meson-python** (needed to compile/stage
    the Fortran binaries into the wheel). Pick meson-python if the LBLRTM/LNFL makefiles can be
    wrapped; otherwise a scikit-build-core custom build that invokes the existing
    `make -f make_lblrtm` / `make_lnfl` and installs the resulting exe into `absco/_bin/`.
  - `requires-python = ">=3.10"`; deps with floors: `numpy>=1.24, scipy>=1.10, pandas>=2.0,
    netCDF4>=1.6, xarray>=2023.1, platformdirs>=3, zenodo_get>=1.3`.
  - `[project.scripts]`: `absco-config = absco.cli.config:main`,
    `absco-generate = absco.cli.generate:main`, `absco-init = absco.cli.init:main`
    (end-user artifact setup), plus `absco-read`, and a dev-only `absco-build` (repackaged
    `build_models` Fortran compile).
  - Package data: include `absco/data/**` and `absco/_bin/**` via package-data / wheel install.

### 2. Path & data resolution (`paths.py`) — enables "install once, run anywhere"
- **Bundled data** (PT_grid, VMR, CSV, template ini): resolve via
  `importlib.resources.files("absco.data")`. Config keys `pfile/ptfile/vmrfile/hdofile/xs_lines`
  default to these packaged resources; a user only overrides when pointing at custom files.
- **Runtime artifacts** (line file, executables): resolve a data dir in this order —
  `$ABSCO_DATA_DIR` → `platformdirs.user_data_dir("absco")`. Executables (`lnfl_path`,
  `lbl_path`) resolve to `<data_dir>/bin/` (or the wheel-bundled `absco/_bin/` when present);
  line-file paths (`tape1_path`, `tape2_path`, `extra_params`, `xs_path`, `fscdxs`) resolve to
  `<data_dir>/AER_Line_File/...`.
- **Refactor `preprocess.py`**: replace `gitDir = os.path.dirname(__file__)` joins
  (`ABSCO_preprocess.py:48-53`) with `paths.py` resolution. When a config value is empty/omitted,
  fall back to the packaged/data-dir default instead of erroring — so a minimal user config
  (only `wn1/wn2/outres/molnames`) is valid.
- **Working-dir / temp files**: keep the existing behavior that already satisfies the
  requirement — `intdir='.'` → `os.getcwd()` (`ABSCO_preprocess.py:302-304`) and run dirs
  (`LNFL_Runs`, `LBL_Runs`, `TAPE3_dir`, `TAPE5_dir`, `outdir`) created under it. Audit the
  `os.chdir(self.gitDir)` calls in `compute.py` (lines 82, 279/328, 590/769): `self.gitDir`
  must become "the working dir captured at start" (still `os.getcwd()`), NOT a repo location —
  verify chdir round-trips land back in the user's working dir, not a package path.

### 3. Vendor `common/`
- Copy the modules actually imported (`utils.py` — `file_check`; `lblTools.py`; `RC_utils.py`;
  `FortranFile.py`) into `src/absco/_common/`. Update imports to `from absco._common import utils`.
- Remove `common` from `.gitmodules` runtime dependence (may keep submodule for dev sync, but
  package must not rely on it at runtime).

### 4. Fortran binaries & artifact initialization
Split `build_models.py` into three clearly separated concerns:
- (a) **Compile Fortran** (`make -f make_lnfl` / `make_lblrtm`, chdir into `<sub>/build`, glob the
  produced exe — `build_models.py:144-235`).
- (b) **Zenodo line-file fetch/extract** (`getLineFile`, `build_models.py:105-142`).
- (c) **Stage** compiled/bundled exes + extracted line file into the user data dir layout that
  `paths.py` expects (`<data_dir>/bin/`, `<data_dir>/AER_Line_File/...`).

Mapped to three usage paths:
- **Wheel build (CI, cibuildwheel per platform)** — runs (a) and bundles the exes into
  `absco/_bin/` inside the wheel. The line file is far too large for PyPI, so it is **never**
  in the wheel.
- **`absco-init` (end user, the requested command)** — the primary onboarding step for someone
  who installed a prebuilt binary wheel. It: verifies the wheel-bundled `absco/_bin/` binaries
  are present (and stages them into the data dir if the layout needs them there), runs (b) to
  download+extract the line file, and runs (c) to place everything under `$ABSCO_DATA_DIR` /
  `platformdirs.user_data_dir("absco")`. Idempotent: skips artifacts already present (mirrors
  the existing `checkLineFile` "already exists" logic, `build_models.py:72-103`). Flags:
  `--data-dir`, `--record <zenodo-id>` (default 3837550, from `build_models.py:258`),
  `--force`, `--lines-only`.
- **`absco-build` (developer only)** — runs (a) from a source checkout to (re)compile the
  Fortran when not using a prebuilt wheel; can then feed (c) so a from-source dev has the same
  data-dir layout an `absco-init` user would.

- **Delete** the old behavior where `configFile()` rewrites absolute paths into the ini
  (`build_models.py:186-235`) — replaced by runtime resolution in `paths.py`; `absco-init`
  populates artifacts, it does not edit any config file.

### 5. Config assistant (`config_assistant.py` + `cli/config.py`) — `absco-config`
- Inputs via flags with interactive fallback: `--outres` (required), `--wn1/--wn2`
  (or `--range`), `--molnames` (space list, lowercased), `--units` (default cm-1),
  `-o/--output` (default `./ABSCO_config.ini`), optional overrides (`--wv-vmr`, custom data files).
- **Suggest `lblres`** from `outres` using the verified formula
  (`power=ceil(log2(outres/1.5e-4)); lblres=outres/2**power`), supporting per-band arrays; then
  re-verify the power-of-2 ratio survives `.ini` string formatting (mirror notebook cells 8–9).
- Load `absco/data/ABSCO_config.template.ini`, fill `[channels]` (`wn1/wn2/lblres/outres/units`)
  and `[molecules] molnames`; leave data/exe/line-file paths **blank** so runtime resolution
  applies. Write the file and print the suggested `lblres` and RAM estimate.
- Reuse `configparser` exactly as the notebook/preprocess do so the output is round-trip valid
  for `absco-generate`.

### 6. `absco-generate` (was `run_LBLRTM_ABSCO.py`)
- Same argparse surface (`-i/--config_file` default `./ABSCO_config.ini`, `-lnfl`, `-lbl`,
  `-e2e`, `-db`, `-y`). Wrap the existing flow (`run_LBLRTM_ABSCO.py:41-113`) in a `main()`
  entry point. No logic change beyond the import/path refactors above.

### 7. pixi configuration
- `pixi.toml` (or `[tool.pixi]`): conda-forge deps + Fortran compiler (`gfortran`);
  a `pypi-dependencies` self-editable-install of `absco`.
- Tasks: `build-fortran` (compile LBLRTM/LNFL, wraps `absco-build`), `init` (wraps `absco-init`
  to fetch+stage artifacts), `install` (editable), `test`, `generate`, `config`. This gives the
  "pixi-controlled project" for dev + a documented build path feeding the wheel CI.

### 8. Docs & cleanup
- This plan lives at `docs/packaging_plan.md`.
- Update `README.md`/`AGENTS.md` install+usage sections (`pip install absco` / pixi; new commands).
- Remove stale artifacts (`ABSCO_preprocess.pyc`, `__pycache__/`); add `.gitignore` entries for
  working-dir run outputs.

---

## Verification

1. **Build & install**: `pixi run build-fortran` then `pip install -e .` (or build a wheel with
   cibuildwheel and `pip install` it in a clean env). Confirm `absco/_bin/` has `lblrtm`/`lnfl`.
2. **Artifacts (end-user path)**: in a clean env with only the installed wheel, run
   `absco-init` → confirm it downloads/extracts the line file and stages binaries into
   `$ABSCO_DATA_DIR` / user data dir with the layout `paths.py` expects. Re-run to confirm it is
   idempotent (skips existing). Dev path: `absco-build` recompiles Fortran from source and
   populates the same layout.
3. **Config assistant**: from an empty scratch dir,
   `absco-config --outres 0.01 --wn1 4166 --wn2 4358 --molnames ch4 co h2o`
   → writes `./ABSCO_config.ini` with `lblres = 7.8125e-05` (matches notebook output). Assert
   `outres/lblres` is a power of 2 by loading it through `absco.preprocess.configure(...)` with no error.
4. **Path independence**: run the above and `absco-generate` from a directory unrelated to the
   source tree; confirm (a) no absolute source paths appear in the ini, (b) `LNFL_Runs/`,
   `LBL_Runs/`, `TAPE3_dir/`, `TAPE5_dir/`, `nc_ABSCO/` are all created under **cwd**, and
   (c) bundled `PT_grid`/`VMR`/CSV are found without user paths.
5. **End-to-end smoke**: `absco-generate -e2e -y -i ABSCO_config.ini` on a narrow band + single
   molecule (e.g. h2o, small wn range) produces a netCDF; `absco-read` opens it.
6. **Custom-file override**: point `vmrfile`/`pfile` at a user file and confirm it takes precedence
   over the packaged default.
7. **Regression**: compare a small table generated post-refactor against one from the current
   in-place code for identical inputs (values match within tolerance).

## Risks / notes
- **cibuildwheel + Fortran** is the hardest piece (per-platform gfortran, manylinux toolchain,
  macOS). If wheels stall, the interim fallback is a from-source install + `absco-build` compiling
  on the user's machine, then `absco-init --lines-only` for the line file — the runtime resolution
  in `paths.py` supports both bundled-binary and compiled-locally without code change.
- **RESOLVED (item 4 follow-up): the Fortran build now succeeds end to end on linux-64.** The
  submodules were bumped to LBLRTM v12.17 and LNFL master (v3.2-30), and three things were needed:
  (1) **recursive submodule init** — both models now carry a nested `aer_rt_utils` submodule
  (LBLRTM also `cross-sections`), and `src/util_gfortran.f90` is a symlink into it; the earlier
  "No rule to make target 'util_gfortran.f90'" was an empty nested submodule, NOT a makefile bug.
  Run `git submodule update --init --recursive`. (2) the **`absco-build` PATH shim** injecting
  `-std=legacy -fallow-argument-mismatch` for the legacy AER code — with this, **the latest
  gfortran (15.2.0) builds both models cleanly**, so `pixi.toml` uses `gfortran = "*"` (no version
  pin; an 11.2.0 pin was tried and then removed once 15.2 was confirmed). (3) **netcdf-fortran**
  in the env — LBLRTM v12.17 unconditionally compiles its netCDF read module, so `compile_model`
  passes `NETCDF=yes NCI/NCL` pointing at the env's netCDF. Both ELF x86-64 executables
  (`lnfl_v3.2_linux_gnu_sgl`, `lblrtm_v12.17_linux_gnu_dbl`) build via `pixi run absco-build`,
  stage into the data dir, and resolve via `absco.paths`. The prebuilt-binary-wheel path
  (cibuildwheel) remains the intended distribution mechanism; this confirms the from-source
  `absco-build` fallback works.
- LBLRTM/LNFL license terms must permit binary redistribution — confirm before publishing wheels.
- Line file (`AER_Line_File`) is far too large for PyPI; it stays a fetched artifact regardless.
