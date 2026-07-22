# Graph Report - .  (2026-07-22)

## Corpus Check
- Large corpus: 90 files · ~543,498 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder.

## Summary
- 878 nodes · 1897 edges · 48 communities (35 shown, 13 thin omitted)
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS · INFERRED: 3 edges (avg confidence: 0.9)
- Token cost: 100,133 input · 3,500 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Optical Properties (Voigt)|Optical Properties (Voigt)]]
- [[_COMMUNITY_Optical Properties (Standard)|Optical Properties (Standard)]]
- [[_COMMUNITY_Non-LTE Calculations|Non-LTE Calculations]]
- [[_COMMUNITY_RRTMGP Flux Processing|RRTMGP Flux Processing]]
- [[_COMMUNITY_Low-Resolution RT|Low-Resolution RT]]
- [[_COMMUNITY_Atmosphere Layers|Atmosphere Layers]]
- [[_COMMUNITY_RadianceFlux Merging|Radiance/Flux Merging]]
- [[_COMMUNITY_Fortran IO Utils|Fortran I/O Utils]]
- [[_COMMUNITY_FFT Scanning (Single)|FFT Scanning (Single)]]
- [[_COMMUNITY_FFT Scanning (Double)|FFT Scanning (Double)]]
- [[_COMMUNITY_Plotting Utils|Plotting Utils]]
- [[_COMMUNITY_Continuum Models|Continuum Models]]
- [[_COMMUNITY_Post-Processing|Post-Processing]]
- [[_COMMUNITY_LBLRTM Core|LBLRTM Core]]
- [[_COMMUNITY_ABSCO Workflow|ABSCO Workflow]]
- [[_COMMUNITY_Configuration Processing|Configuration Processing]]
- [[_COMMUNITY_Table Joining|Table Joining]]
- [[_COMMUNITY_Documentation|Documentation]]
- [[_COMMUNITY_ABSCO Table Generation|ABSCO Table Generation]]
- [[_COMMUNITY_TAPE3 Generation|TAPE3 Generation]]
- [[_COMMUNITY_VMR Profiles|VMR Profiles]]
- [[_COMMUNITY_RRTMG Conversion|RRTMG Conversion]]
- [[_COMMUNITY_Build System (Root)|Build System (Root)]]
- [[_COMMUNITY_Build System (Common)|Build System (Common)]]
- [[_COMMUNITY_Solar Radiation|Solar Radiation]]
- [[_COMMUNITY_Zenodo API|Zenodo API]]
- [[_COMMUNITY_ABSCO Table Reader|ABSCO Table Reader]]
- [[_COMMUNITY_Instrument SRF|Instrument SRF]]
- [[_COMMUNITY_File Utils|File Utils]]
- [[_COMMUNITY_Join Multiple Script|Join Multiple Script]]
- [[_COMMUNITY_Earth Constants|Earth Constants]]
- [[_COMMUNITY_MT_CKD Continuum|MT_CKD Continuum]]
- [[_COMMUNITY_Mars Constants|Mars Constants]]
- [[_COMMUNITY_Test Module|Test Module]]
- [[_COMMUNITY_HDO Profiles|HDO Profiles]]
- [[_COMMUNITY_Dependency Specs|Dependency Specs]]
- [[_COMMUNITY_LBL Parameters (F90)|LBL Parameters (F90)]]
- [[_COMMUNITY_Physical Constants|Physical Constants]]
- [[_COMMUNITY_Solar Cycle|Solar Cycle]]
- [[_COMMUNITY_Multi-Config Runner|Multi-Config Runner]]
- [[_COMMUNITY_CLAUDE Instructions|CLAUDE Instructions]]
- [[_COMMUNITY_Common README|Common README]]
- [[_COMMUNITY_License|License]]
- [[_COMMUNITY_Spectral Degradation|Spectral Degradation]]

## God Nodes (most connected - your core abstractions)
1. `lblparams` - 54 edges
2. `lblparams` - 53 edges
3. `tips_2003()` - 47 edges
4. `atob()` - 47 edges
5. `tips_2003()` - 47 edges
6. `atob()` - 47 edges
7. `lblparams` - 28 edges
8. `lblparams` - 25 edges
9. `contnm()` - 24 edges
10. `lblparams` - 21 edges

## Surprising Connections (you probably didn't know these)
- `Conda Environment Specification` --semantically_similar_to--> `Python Dependencies`  [INFERRED] [semantically similar]
  environment.yml → requirements.txt
- `LBLRTM FAQ` --references--> `LBLRTM`  [EXTRACTED]
  LBLRTM/docs/FAQ_LBLRTM.pdf → LBLRTM/README.md
- `Optical Depth Calculation` --references--> `LBLRTM`  [EXTRACTED]
  README.md → LBLRTM/README.md
- `LBLRTM FAQ` --references--> `LNFL`  [EXTRACTED]
  LBLRTM/docs/FAQ_LBLRTM.pdf → LBLRTM/README.md
- `ABSCO Table Generation` --references--> `AIRS Pressure Grid`  [EXTRACTED]
  README.md → PT_grid/AIRS_P_air.txt

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **ABSCO Generation Workflow** — lblrtm_readme_lnfl, lblrtm_readme_lblrtm, readme_md_absco, readme_md_build_models_py, readme_md_run_lblrtm_absco_py [EXTRACTED 1.00]
- **LBLRTM Input File System** — lblrtm_readme_tape3, lblrtm_readme_tape5, lblrtm_readme_hitran, readme_md_aer_line_file [EXTRACTED 1.00]
- **PT Grid Configuration** — pt_grid_airs_pressures, pt_grid_temperature_array, readme_md_vmr_profiles [EXTRACTED 1.00]

## Communities (48 total, 13 thin omitted)

### Community 0 - "Optical Properties (Voigt)"
Cohesion: 0.07
Nodes (84): armstrong(), atob(), chi_fn(), cljust(), cnvfnv(), convf4(), lblparams, phys_consts (+76 more)

### Community 1 - "Optical Properties (Standard)"
Cohesion: 0.08
Nodes (85): atob(), chi_fn(), cljust(), cnvfnv(), convf4(), lblparams, phys_consts, struct_types (+77 more)

### Community 2 - "Non-LTE Calculations"
Cohesion: 0.10
Nodes (47): cnvf4q(), cnvfnq(), defnltedat(), dropspace(), lblparams, phys_consts, struct_types, getindex() (+39 more)

### Community 3 - "RRTMGP Flux Processing"
Cohesion: 0.05
Nodes (34): configSetup, lwRRTMGP, Using the netCDF template, start a netCDF that will contain the      fields from, Parse the input .ini file (inFile) and return as a dictionary for      use in th, Combine flux arrays from SW Flux Calculation output files, Compute fluxes for each RRTMGP-defined band      Keywords       broadband -- boo, Class that conforms all of LW spectra generated by RADSUM to the      RRTMGP net, Combine flux arrays from SW Flux Calculation output files (+26 more)

### Community 4 - "Low-Resolution RT"
Cohesion: 0.11
Nodes (44): ab(), abslim(), aerext(), aernsm(), aerprf(), aitk(), bs(), cirr18() (+36 more)

### Community 5 - "Atmosphere Layers"
Cohesion: 0.14
Nodes (42): alayer(), amerge(), andex(), andexd(), atmpth(), autlay(), check(), cmpalt() (+34 more)

### Community 6 - "Radiance/Flux Merging"
Cohesion: 0.19
Nodes (38): abinit(), absint(), absmrg(), absout(), adarsl(), aerf(), bbdtfn(), bbfn() (+30 more)

### Community 7 - "Fortran I/O Utils"
Cohesion: 0.08
Nodes (22): FortranFile, readReflectance(), writeReflectance(), generateHeightGrid(), generatePressureGrid(), getOD(), interP(), Read in binary LBLRTM ODint files (IMRG=1 and IOD=1, 3, or 4 in   LBLRTM specifi (+14 more)

### Community 8 - "FFT Scanning (Single)"
Cohesion: 0.14
Nodes (30): addnin(), bessi0(), bigfft(), boxcar(), ckfile(), compnin(), expan(), fftblk() (+22 more)

### Community 9 - "FFT Scanning (Double)"
Cohesion: 0.14
Nodes (30): addnin(), bessi0(), bigfft(), boxcar(), ckfile(), compnin(), expan(), fftblk() (+22 more)

### Community 10 - "Plotting Utils"
Cohesion: 0.14
Nodes (28): ax2(), axes(), axisl(), axlog(), bbscle(), dbod(), dbtr(), expt() (+20 more)

### Community 11 - "Continuum Models"
Cohesion: 0.19
Nodes (26): cld_od(), contnm(), lblparams, phys_consts, frn296(), frnco2(), herprs(), hertda() (+18 more)

### Community 12 - "Post-Processing"
Cohesion: 0.18
Nodes (25): cnvflt(), cnvrct(), cnvvrc(), cnvvrl(), convsc(), lblparams, fltmrg(), fltprt() (+17 more)

### Community 13 - "LBLRTM Core"
Cohesion: 0.20
Nodes (21): copyfl(), endfil(), lblparams, phys_consts, layer2level(), lblrtm, nwdl(), opdpth() (+13 more)

### Community 14 - "ABSCO Workflow"
Cohesion: 0.09
Nodes (20): combineVMR(), Combine attributes (ABSCO, vmrWV, and vmrO2) from multiple makeABSCO   objects g, call_CR1(), call_CR2(), check_CR2(), check_py3(), ls(), no_overwrite() (+12 more)

### Community 15 - "Configuration Processing"
Cohesion: 0.13
Nodes (12): configure, Read configuration file and set attributes of configure object     accordingly, Parse the input .ini file (inFile) and return as an object for     use in makeAB, Make sure all of the attributes in the configure object that are     necessary f, Read in user and standard atmosphere pressures (associated     with VMRs) and pe, Simple check to make sure HDO profile is consistent with      input profiles of, If the user only specifies a spectral range and no valid     molecules, try to d, Determine for which bands LNFL and LBLRTM should be run based on     the user-in (+4 more)

### Community 16 - "Table Joining"
Cohesion: 0.16
Nodes (7): AbscoFileJoiner, copy_attrs(), main(), Sort open table objects based on extent ranges, object, AbscoConfigSplitter, main()

### Community 17 - "Documentation"
Cohesion: 0.13
Nodes (22): LBLRTM FAQ, HITRAN Molecule List, LBLRTM User Instructions, LBLRTM TAPE File Convention, AER RT Utils, FASCODE, HITRAN Database, LBLRTM (+14 more)

### Community 18 - "ABSCO Table Generation"
Cohesion: 0.12
Nodes (11): makeABSCO, makeSymLinks(), For a given molecule, make a TAPE5 that can be used as input into     an LNFL ru, Loop over input files and make symbolic links for them, Run executable to generate a binary line file (TAPE3) for usage     in LBLRTM., For a given molecule, make a TAPE5 for every temperature and band     that can b, - Build TAPE5s (LNFL and LBLRTM) for each molecule of interest   - Build TAPE3 (, Inputs       mol -- string, name of molecule to process       inObj -- preproc.c (+3 more)

### Community 19 - "TAPE3 Generation"
Cohesion: 0.16
Nodes (8): genTAPE3, Make the TAPE5 LNFL specifications file for each molecule, Run LNFL for each molecule, thus generating a binary TAPE3, Assuming that line files exist for each HITRAN molecule, generate    an associat, Constructor for genTAPE3 class -- where are ASCII line files,      what are the, Make a directory if it doesn't already exist, Remove any TAPE files that were generated, Make symbolic links      Inputs       source -- string, path to file to be linke

### Community 20 - "VMR Profiles"
Cohesion: 0.13
Nodes (8): broadener, After vmrProfiles() object is constructed, run LBLATM (subroutine      of LBLRTM, Generate a small LBLRTM TAPE5 with the parameters necessary to      run LBLATM a, Read in CSV data for both the HITRAN and XS molecules and perform     a linear i, Run LBLATM with the TAPE5s generated in makeT5(), save TAPE7s in      their own, Write a CSV file that consolidates all of the LBLATM broadener      information, Calculate VMR profile on user-specified pressure grid, vmrProfiles

### Community 21 - "RRTMG Conversion"
Cohesion: 0.18
Nodes (7): Merge together the fluxes and heating rates from all profiles      into a single, Combine the flux arrays for all RFMIP experiments (each of which      should hav, Write a netCDF with the data in an rrtmg object. This is done for     each spect, Extract the RRTMG flux files for a given spectral domain      Read a single RRTM, Extract the RRTMG flux files for given spectral domain, Read a single RRTMG ASCII flux file and return the model output      in a dictio, rrtmg

### Community 22 - "Build System (Root)"
Cohesion: 0.18
Nodes (5): Retrieve line file dataset from Zenodo, extract archive, then     stage files as, Build models as needed (one at a time) and replace paths in     ABSCO_tables.py, Replace paths in ABSCO_tables.py configuration file with the paths     establish, Check if line file dataset from Zenodo needs to be downloaded, submodules

### Community 23 - "Build System (Common)"
Cohesion: 0.18
Nodes (5): Build models as needed (one at a time) and replace paths in     ABSCO_compute.py, Retrieve line file dataset from Zenodo, extract archive, then     stage files as, Replace paths in ABSCO_compute.py configuration file with the paths     establis, Check if line file dataset from Zenodo needs to be downloaded, submodules

### Community 24 - "Solar Radiation"
Cohesion: 0.33
Nodes (10): comb_solar(), phys_consts, scale_solar(), solin(), solin2(), solin_sgl(), solin_tri(), solint() (+2 more)

### Community 26 - "Zenodo API"
Cohesion: 0.25
Nodes (4): apiZenodo, `zenodo_request.py -h`, Read in secret Zenodo API Access Token, Zenodo API upload request      https://developers.zenodo.org/?python#quickstart-

### Community 27 - "ABSCO Table Reader"
Cohesion: 0.25
Nodes (4): Read in the netCDF (and time the differences in the load by the      netCDF4 and, - Time the difference between netCDF4 and xarray libraries      - find pressure,, Find array indices that correspond to the user-provided      coordinates (P, T,, testABSCO

### Community 28 - "Instrument SRF"
Cohesion: 0.38
Nodes (6): airsSRFInterpOverlap(), interpolateSRF(), Read in and return the AIRS spectral response function from HDF5    file (amongs, Interpolate the AIRS spectral response function onto the grid given    by an LBL, Designed to help save time in interpolateSRF(), which loops over    all line cen, srfAIRS()

### Community 30 - "File Utils"
Cohesion: 0.50
Nodes (4): file_check(), log(), Quick check if path exists.  Use before reading a file., Function from Andre Wehe code in    /nas/project/rc/rc2/rpernak/NH3_Retrievals/c

### Community 32 - "Earth Constants"
Cohesion: 0.50
Nodes (3): phys_consts, grav_const(), planet_consts

### Community 33 - "MT_CKD Continuum"
Cohesion: 0.67
Nodes (3): MT_CKD Continuum, MT_CKD 3.3, LBLRTM v12.10

## Knowledge Gaps
- **26 isolated node(s):** `phys_consts`, `lblparams`, `phys_consts`, `phys_consts`, `phys_consts` (+21 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **13 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `line_shrink` connect `Non-LTE Calculations` to `Optical Properties (Voigt)`, `Optical Properties (Standard)`?**
  _High betweenness centrality (0.042) - this node is a cross-community bridge._
- **Why does `configure` connect `Configuration Processing` to `ABSCO Workflow`?**
  _High betweenness centrality (0.013) - this node is a cross-community bridge._
- **Why does `convf4()` connect `Optical Properties (Voigt)` to `Non-LTE Calculations`?**
  _High betweenness centrality (0.010) - this node is a cross-community bridge._
- **What connects `Loop over input files and make symbolic links for them`, `- Build TAPE5s (LNFL and LBLRTM) for each molecule of interest   - Build TAPE3 (`, `Inputs       mol -- string, name of molecule to process       inObj -- preproc.c` to the rest of the system?**
  _123 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Optical Properties (Voigt)` be split into smaller, more focused modules?**
  _Cohesion score 0.07471264367816093 - nodes in this community are weakly interconnected._
- **Should `Optical Properties (Standard)` be split into smaller, more focused modules?**
  _Cohesion score 0.0764501470195135 - nodes in this community are weakly interconnected._
- **Should `Non-LTE Calculations` be split into smaller, more focused modules?**
  _Cohesion score 0.10180995475113122 - nodes in this community are weakly interconnected._