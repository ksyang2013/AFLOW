3.1.212 - 2018/11/12
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.212.tar.xz
        -- in --bzd command, added option to transform high-symmetry kpaths to input lattice representation with [--transform2original] (DX)
        -- in --bzd command, added option to print transformation matrices between input lattice and AFLOW standard lattice [--print_transformation_matrix] (DX)
        -- store both the coordinate system and unit cell (rigid rotation) transformations in xstructure (DX)
        -- created _kpoints class; easier to transform and print kpoint information (DX)
        (AFLOW: aflow.h, aconvasp_main.cpp, lattice.cpp, pflow.h, xatom.cpp, README_AFLOW_ACONVASP.TXT)
3.1.211 - 2018/10/19
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.211.tar.xz
        -- added Pearson coefficient to element properties (ME)
        -- fixed some constants in xscalar (ME)
        -- Rewrote AAPL to be faster, require less DFT calculations, and include four-phonon processes. (ME)
        -- Added more functionality to xvector and xmatrix for complex numbers. (ME)
        -- Added ability to update progress bar using percentages instead of indices (ME)
        -- Redesigned xtensor to be more lightweight and faster. (ME)
        -- Fixed bug in aurostd::xcombos::reset() (ME)
        -- Fixed minor typos (ME)
        -- Edited APL readme for grammar and clarity, and added new AAPL features (ME)
        -- Added APL and AAPL parameters to the aflow.rc (ME)
        -- Reformatted exception readme to be more consistent with other readmes and fixed typos (ME)
        -- updated get_datetime_formatted() and StringsAlphabetic() (CO)
        -- new chull date format avoids time stamp, too much resolution (CO)
        -- added png output option for chull and resolution option to aflowrc (CO)
        -- switched to tight tol spacegroups (vsg2) in chull (CO)
        -- centralized checking sign of distances in chull, then flipped sign of distance to hull to be positive by default (CO)
        -- funneled all points of facet through addVertex() (CO)
        -- fixed facet content tolerance scaling with increasing dimensionality, important for 6D hulls (CO)
        -- added default initializations to a bunch of variables in chull to avoid spurious warnings (CO)
        -- added polymorph case to decomposition phases and coefficients (CO)
        -- allow for eq_g_states to be calculated on the fly if not already calculated (CO)
        -- intelligently avoid printing unknown (or unset) stability criterion (CO)
        -- fixed unsorted directories in aflow_compare_structure (CO)
        -- now check for negative coefficients when balancing chemical equations, means we missed the facet (CO)
        -- correctBadDatabase() now checks for unaries too (as per apennsy) (CO)
        -- skipping ".old" prototypes (CO)
        -- fixed shidong warnings with strncat (CO)
        (AUROSTD: boot.cpp, AUROSTD/aurostd.h, main.cpp, xcombos.cpp, xcombos.h, xcomplex.cpp, xcomplex.h, xmatrix.cpp, xmatrix.h, xtensor.cpp, xtensor.h, xvector.cpp, xvector.h)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, chull.cpp, chull.h, compare_structure.cpp, contrib_shidong_auxiliary.cpp, aflowlib_web_interface.cpp, pflow_funcs.cpp, APL/aapl_cluster.cpp, APL/aapl_ifcs.cpp, APL/aapl_setup.cpp, APL/aapl_tcond.cpp, APL/apl_atomic_disp.cpp, APL/apl_dirphoncalc.cpp, APL/apl_dm_pdos_save.cpp, APL/apl_doscalc.cpp, APL/apl_group_velocity.cpp, APL/apl_gsa.cpp, APL/apl.h, APL/apl_hsqpoints.cpp, APL/apl_kphonons.cpp, APL/apl_logger.cpp, APL/apl_lrphoncalc.cpp, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/apl_supercell.cpp, APL/apl_thermalpc.cpp, APL/qha3phonons_eos.cpp, APL/qha_aflowin_creator.cpp, APL/qha_dm_pdos_save.cpp, APL/qha_eosrun.cpp, APL/qha_gruneisen.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_APL.TXT, README_AFLOW_CHULL.TXT, README_AFLOW_EXCEPTIONS.TXT)
3.1.210 - 2018/10/01
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.210.tar.xz
        -- changed 8 ANRL prototype labels (part 2) to match with Wyckoff positions in reference literature (DX)
        -- fixed prototype names for A4BC4D_tP10_123_gh_a_i_d and AB_hP6_144_a_a in aflow_anrl_list.cpp (DX) 
        -- added Strukturbericht designations to aflow_anrl_list.cpp (DX)
        -- fixed minimum enumerated Wyckoff letter determination; more robust (DX)
        -- cleaned workflow for Wyckoff functions (DX)
        -- added more debugging messages to symmetry functions (DX)
        (AFLOW: ANRL/aflow_anrl_list.cpp ANRL/aflow_anrl_A12B36CD12_cF488_210_h_3h_a_fg.cpp ANRL/aflow_anrl_A3B3C_hP14_176_h_h_c.cpp ANRL/aflow_anrl_A3B_hP8_176_h_c.cpp ANRL/aflow_anrl_AB3C_cP60_201_be_fh_g.cpp ANRL/aflow_anrl_A3B2_hP10_176_h_bc.cpp ANRL/aflow_anrl_A3BC_hP10_188_k_c_a.cpp ANRL/aflow_anrl_AB3C16_cF160_203_a_bc_eg.cpp ANRL/aflow_anrl_AB_hP4_156_ab_ab.cpp Makefile anrl.cpp symmetry_spacegroup.cpp symmetry_spacegroup.h symmetry_spacegroup_functions.cpp)
3.1.209 - 2018/08/31
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.209.tar.xz
        -- added all origin choice 2 possibilities to ITC space group list (DX)
        -- added all rhombohedral setting possibilities to ITC space group list (DX)
        -- added monoclinic unique axis choices (b or c) to ITC space group list (DX)
        -- added SYM::findRhombohedralSetting() in space group functions to distinguish between hexagonal and rhombohedral settings (DX)
        -- added other settings for Hall space group symbol (settings 1 or 2, H or R for rhombohedral, and unique axis-b or -c for monoclinic systems) (DX)
        -- added CIF reader; can read/expand CIFs with only representative Wyckoff position specified in settings 1 or 2 (H or R for rhombohedral systems and unique axis-b or -c for monoclinic systems) (DX)
        -- improved CIF writer; default functionality calculates the symmetry and prints the representative Wyckoff positions (DX)
        -- added symmetry tolerance option to CIF writer via --cif[=<tolerance>] (default: tight; other options: loose, <number>) (DX)
        -- added space group setting option to CIF writer via [--setting=1| =2] (default: 1) (DX)
        -- added [--no_symmetry] option to CIF writer, which returns CIF as space group 1 (DX)
        -- added space group setting option to --wyccar via [--setting=1| =2] (default: 1) (DX)
        -- added space group setting option to --sgdata and --edata via [--setting=1| =2] (default: 1) (DX)
        -- added [--no_scan] and [--magmom] options to --wyccar and moved function to pflow::WYCCAR() (DX)
        -- added more debugging messages to space group functions (ConventionalCell(), SpaceGroup_ITC(), etc.) (DX)
        -- changed ANRL prototype A6B2C_tP18_128_eh_d_b to A6B2C_tP18_128_eh_d_a (b vs a Wyckoff letter); consistent with reference (DX)
        -- fixed bug in minimumDistanceVector(); should return xvector<double> not double (DX)
        (AFLOW: ANRL/aflow_anrl_list.cpp, ANRL/aflow_anrl_A6B2C_tP18_128_eh_d_a.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_SYM.TXT, aflow.h, aconvasp_main.cpp, anrl.cpp, pflow.h , pflow_print.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup.h, symmetry_spacegroup_ITC_library.cpp, symmetry_spacegroup_functions.cpp, xatom.cpp)
3.1.208 - 2018/08/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.208.tar.xz
        -- new small banner style in chull (CO)
        -- perform thermo loop with static run in aflowlib_libraries (CO)
        -- added xoption to edata and sgdata functions; stores all data (DX)
        -- added more keywords to aflowlib entry (crystal_family, crystal_system, point_group_Hermann_Mauguin, point_group_Schoenflies, point_group_orbifold, point_group_type, point_group_order, point_group_structure, Bravais_lattice_lattice_type, Bravais_lattice_lattice_system, Bravais_superlattice_lattice_type, Bravais_superlattice_lattice_variation_type, Bravais_superlattice_lattice_system, Pearson_symbol_superlattice, reciprocal_geometry, reciprocal_volume_cell, reciprocal_lattice_type, reciprocal_lattice_variation_type, Wyckoff_letters, Wyckoff_multiplicities, Wyckoff_site_symmetries) (DX)
        -- new aflowlib keywords extracted from xoption and integrated into aflowlib.out/.json (DX)
        -- prepared website function to read in new keywords once database is populated; currently commented out (DX)
        -- fixed typo in sgdata; Shoenflies is now Schoenflies (DX)
        -- added functions to create Wyckoff strings in aflowlib entry; ExtractWyckoffLettersString(), ExtractWyckoffMultiplicitiesString(), and ExtractWyckoffSiteSymmetriesString() (DX)
        (AUROSTD: xoption.cpp)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, pflow.h, pflow_print.cpp, symmetry_spacegroup.h, symmetry_spacegroup_functions.cpp, aflowlib.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp)
3.1.207 - 2018/08/19
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.207.tar.xz
        -- integrated new html for entry page (JPO)
        -- fixed xaxis of bands plotter for entry page (PC)
        -- Improved distribution of APL calculations over threads. (ME)
        -- Added function to SYM that returns the minimum distance vector. (ME)
        -- Added docstring to xcombos. (ME)
        -- added xvector::normalizeSumToOne() (CO)
        -- changed decomposition reaction to atomic concentrations (CO)
        -- added fractional_compound to chull output (CO)
        -- replaced AFLOWLogicError() and AFLOWRuntimeError() with xerror() (CO)
        -- added chull plot ICSD labels mode (CO)
        -- revamped small banner setting (CO)
        -- added verbose output if skipping entries above/below half hull (CO)
        -- added --output=png --png_resolution=500 options to chull (CO)
        -- added --keep=gpl to phonons gnuplot script (CO)
        (AUROSTD: boot.cpp, AUROSTD/aurostd.h, xcombos.cpp, xscalar.cpp, xscalar.h, xvector.cpp, xvector.h)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, avasp.cpp, bader.cpp, chull.cpp, chull.h, aflow.cpp, aflow.h, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, pflow_funcs.cpp, pflow.h, pocc.cpp, symmetry.cpp, APL/aapl_tcond.cpp, APL/apl_atomic_disp.cpp, APL/apl_dm_pdos_save.cpp, APL/apl_doscalc.cpp, APL/apl_group_velocity.cpp, APL/apl.h, APL/apl_hroutines.cpp, APL/apl_pdisc.cpp, APL/qha3phonons_eos.cpp, APL/qha_dm_pdos_save.cpp, APL/qha_gruneisen.cpp, APL/scqha_eos.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_AFLOWRC.TXT, README_AFLOW_CHULL.TXT)
3.1.206 - 2018/08/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.206.tar.xz
        -- Fix JVXL (SC)
3.1.205 - 2018/07/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.205.tar.xz
        -- Improved destructor for xtensor (ME)
        -- Added function to test if xvector is a zero vector (ME)
        -- fixed parameter list for ANRL prototype (A_hR105_166_bc9h4i); z2 to x2 (DX)
  	-- more rigorous check for atoms within the new cell in GetSuperCell(); check periodic images of atoms (DX)
  	-- fixed "over-reduced" issue in primitivization routine in space group function (DX)
  	-- added roundoff to axis and SU2 matrix in JSON output of symmetry elements (DX)
  	-- fixed function for printing fractions in general Wyckoff positions (DX)
        -- created --generate_aflowin_only (CT)
        -- QHA3P and SCQHA conflict resolved (PN)
        -- reorganize QHA modes (PN)
        -- added four QHA3P options (PN)
        -- reorganize QHA, SCQHA, and QHA3P options (PN)
        -- filename changed and accordingly modified in QHA-README (PN)
        -- replaced some QHA functions with aurostd (PN)
        -- fixed slab AddAtom() function (CO+DU)
        -- added --readme=aflowrc (CO)
        -- fixed --readme=xaflow for local configuration (CO)
        -- fixed a few warnings for beta (CO)
        -- fixed AMIX/BMIX typo (CO)
        -- added defaults/options for DOD CONRAD machine (--machine=dod_conrad); also added in aflowrc (DX)
        (AUROSTD: xscalar.cpp, xscalar.h, xtensor.cpp, xvector.cpp, xvector.h)
        (AFLOW: aflowrc.cpp, avasp.cpp, chull.cpp, aflow.cpp, data.cpp, aflow.h, init.cpp, ivasp.cpp, kbin.cpp, kvasp.cpp, aflowlib_webapp_entry.js, aflowlib_web_interface.cpp, pocc_old.cpp, surface.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_functions.cpp, xatom.cpp, xclasses.cpp, ANRL_CPPS_20180710/anrl_list.cpp, APL/apl_group_velocity.cpp, APL/apl.h, APL/apl_kphonons.cpp, APL/qha3phonons_eos.cpp, APL/qha_aflowin_creator.cpp, APL/qha_energies.cpp, APL/qha_eos.cpp, APL/qha_gruneisen.cpp, APL/scqha_eos.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_AEL.TXT, README_AFLOW_AFLOWRC.TXT, README_AFLOW_AGL.TXT, README_AFLOW_POCC.TXT, README_AFLOW_QHA_SCQHA_QHA3P.TXT, README_AFLOW.TXT, README_AFLOW_XAFLOW.TXT)
3.1.204 - 2018/07/12
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.204.tar.xz
        -- updated QHA and added QHA3P and SCQHA functionality (PN)
        -- Restrucured xcombos. Added enumerations. (ME) [AUROSTD/aflow_xcombos.cpp, AUROSTD/aflow_xcombos.h]
        -- Introduced exception handlier class aurostd::xerror. (ME)
        -- Introduced xtensor class for tensors of arbitrary dimension. Other xtensor classes are obsolete now. (ME).
        -- Changed all current xtensor3 instances to the new xtensor format. (ME)
        -- Added Kronecker product to xmatrix. (ME)
        -- Option to write and use AEL data at lowest finite pressure where the material is elastically stable (CT)
        -- Option to specify a separate set of finite pressures for AEL calculations than are used for AGL post-processing (CT)
        -- Integrated workflow option for using finite pressure Poisson ratio in AGL calculations (CT)
        -- Write nominal target pressure and calculated external pressure for AEL calculations in aflow.ael.out file (CT)
        -- Fixed keyword in JSON output file (CT)
        -- Increase grid for AGL DOSCAR (CT)
        -- cleaning up webapp_bands.js and added more mouse functionality (PC)
        -- added citation information to entry page (PC)
        -- added ANRL prototypes from part 2 of library (302 prototypes) (DX)
        -- added option to print symbolic math representation of ANRL prototypes (--add_equations or --equations_only) in aims and vasp formats (DX)
        -- updated DOI for ANRL part 1 and added arXiv for part 2 to title line of each prototype (DX) 
        -- fixed bug with --vasp keyword, i.e., need to check if used with --proto command (DX) 
        -- removed unused variables in standard lattice function (DX)
        -- added applyCombo() to xcombos (CO)
        -- fixed wget *.xz issue if *.xz already exists (CO)
        -- fixed apl2agr to handle xz compression format (CO)
        -- fixed --aflow_proto empty BZ issue, occurs when structure is downloaded from online (CO)
        -- fixed reading in forces from aims.out (CO)
        -- force species input for aims structure (CO)
        -- added checks for broken API (CO)
        -- added GENERATOR to aflow.in generation (CO)
        -- added auid to chull PDF output (CO)
        -- fixed bader num_each_type issue (relax2 vs. static primitivization) (CO)
        -- fixed command line aflow_qmvasp generation (CO)
        -- fixed search for ./aflow_data issue upon initial installation (CO)
        -- PARTCAR now handles '+' in write out (CO)
        -- fixed AddAtom() to handle different occupations (CO)
        -- added combination parameters to chull (CO)
        -- fixed aflowrc load issue (remove spaces between quotes and comment) (CO)
        (AUROSTD: boot.cpp, AUROSTD/aurostd.cpp, AUROSTD/aurostd.h, main.cpp, xcombos.cpp, xcombos.h, xerror.cpp, xerror.h, xmatrix.cpp, xmatrix.h, xtensor.cpp, xtensor.h)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, anrl.cpp, avasp.cpp, bader.cpp, bader.h, chull.cpp, chull.h, aflow.cpp, data.cpp, aflow.h, init.cpp, ivasp.cpp, kbin.cpp, lattice.cpp, aflowlib_libraries.cpp, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, oaims.cpp, ovasp.cpp, pflow.h, pocc.cpp, pocc.h, surface.cpp, xatom.cpp, xclasses.cpp, xproto.cpp, ANRL_CPPS_20180710, APL/apl_atomic_disp.cpp, APL/apl_doscalc.cpp, APL/apl_group_velocity.cpp, APL/apl.h, APL/apl_hroutines.cpp, APL/apl_hsqpoints.cpp, APL/apl_kphonons.cpp, APL/apl_ltetdos.cpp, APL/apl_mpmesh.cpp, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/apl_thermalpc.cpp, APL/apl_uniform_mesh.cpp, APL/qha3phonons_eos.cpp, APL/qha_aflowin_creator.cpp, APL/qha_dm_pdos_save.cpp, APL/qha_energies.cpp, APL/qha_eoscalc.cpp, APL/qha_eos.cpp, APL/qha_gruneisen.cpp, APL/scqha_eos.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_AEL.TXT, README_AFLOW_AGL.TXT, README_AFLOW_ANRL.TXT, README_AFLOW_APL.TXT, README_AFLOW_CHULL.TXT, README_AFLOW_EXCEPTIONS.TXT, README_AFLOW_POCC.TXT, README_AFLOW_QHA_SCQHA_QHA3P.TXT, README_AFLOW.TXT, README_AFLOW_XAFLOW.TXT, README_CONTRIBS.TXT)
3.1.203 - 2018/06/15
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.203.tar.xz
        -- fixed tolerance scan issue in edata (try new tolerance from PrintSGData() on GetLatticeType() routine before changing tolerance) (DX)
        -- fixed primitivization routine in aflowSG functions (Minkowski/Niggli to fix left-handed candidate lattices and recheck moduli after Minkowski/Niggli) (DX)
        -- check all possible generator choices to match to ITC convention before changing tolerance (DX)
        -- added more rigorous check of tetragonal symmetry operations for determining conventional cell (DX)
        -- added generator information for P1 symmetry systems (DX)
        -- added check of Cartesian distance before removing fractional copies in GetPrimitive() (DX)
        (AFLOW: lattice.cpp, pflow_print.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_ITC.cpp, xatom.cpp) 
3.1.202 - 2018/06/07
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.202.tar.xz
        -- dropping EXTRA to AFLOW3_FREE/EXTRA on local machine or wget (SC)
3.1.201 - 2018/06/04
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.201.tar.xz
        -- more space saving in lib2raw (linking OUTCAR.relax instead of copying) SC
3.1.200 - 2018/05/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.200.tar.xz
        -- fixed structure rescaling issue for space group determination (affecting LIB4) (DX)
        -- print geometry file location for errors/debug (DX)
        -- added more debugging messages along with file location in symmetry functions (DX) 
        -- initialize variables for -O3 in symmetry and structure comparison routines (DX)
        -- speed increase for minimumDistance() function (DX)
        -- speed increase for primitivization routine in aflowSG functions (DX)
        -- fixed Wyckoff position typo for space group 62 (8d not 8e) (DX)
        -- fixed Wyckoff position typo for space group 89 (4i not 2i) (DX)
        (AFLOW: aflow.h, aconvasp_main.cpp, compare_structure.cpp, compare_structure_function.cpp, pflow_print.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_ITC_library.cpp, symmetry_spacegroup_functions.cpp, xatom.cpp, xproto.cpp, aflowlib_libraries.cpp)
3.1.199 - 2018/05/27
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.199.tar.xz
	      -- compressing aflow.***.json/out in LIB2RAW and linking them (SC)
        -- removing the brainy/useless stuff about file compression in aurostd (SC)
	      -- fighting sloppyness in aflow_libraries about aflow_pgroup[x][_xtal] (SC)
	      -- fixed aflow_libraries about compress, delete and link files (SC)
	      -- fixed inconsitencies in aflow_convasp_main compress (SC)
	      -- fixed as bunch of inconsistencies in compressing and XHOST.command (SC)
	      -- more inconsistencies RAW-WEB in aflow_libraries  (SC)
	      -- fixed removal useless files in LIB2RAW (SC)
	      -- added BZ2XZ engine to aurostd_main (SC)
	      -- added GZ2XZ engine to aurostd_main (SC)
	      -- added ZIP2ZIP engine to aurostd_main (SC)
3.1.198 - 2018/05/24
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.198.tar.xz
        -- fixed bader extension finder for --lib2raw (CO)
        (AFLOW: bader.cpp, bader.h, aflowlib_libraries.cpp, Makefile)
3.1.197 - 2018/05/24
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.197.tar.xz
        -- updates to entry page bands plotter (PC)
        -- added some helpful comments for getGeneralNormal(), CMdet(), and ZVAL in bader (CO)
        -- check if /www directory exists for jmol display on entry page (CO)
        (AUROSTD: xmatrix.cpp, xvector.cpp)
        (AFLOW: bader.cpp, chull.cpp, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, Makefile)
3.1.196 - 2018/05/21
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.196.tar.xz
        -- fixed entry page property line wrapping issues (PC)
        -- added button for spin-polarized band structure selection in webapp (PC)
        -- fixed precision inconsistency printing xstr.json (FK)
        -- added some -O3 compatibility (FK)
        -- citation added to aflow_aapl_pairs.cpp (ME)
        -- fixed bug in aflow_kvasp.cpp that prevented the creation of the primitive cell structure for APL, AAPL, and QHA calculations (ME)
        -- implemented combinations with repetitions with sequence taken into account (ME)
        -- moved the code to check for input and output files in APL into separate functions to make them available for AAPL, remove duplicate code, and make the code more readable (ME)
        -- removed the ENCUT and EDIFF tags from AAPL input files as they may result in lower cut-off energies and higher energy differences (ME)
        -- fixed gcc8 issue in AGL (CO)
        -- added simplex content and hypercollinearity properties to convex hull data (CO)
        -- decoupled internal links between graph2doc and withindoc (CO)
        -- fixed site error calculation in POCC for vacancies (CO)
        -- added eyes/ones xmatrix constructors (CO)
        -- added Cayley-Menger determinant to xmatrix (CO)
        -- included cstdlib in aflow_data.cpp for compilation on qrats (CO)
        (AUROSTD: boot.cpp, xcombos.cpp, xcombos.h, xmatrix.cpp, xmatrix.h, xvector.cpp)
        (AFLOW: aconvasp_main.cpp, apennsy_main.cpp, chull.cpp, chull.h, contrib_shidong_auxiliary.cpp, data.cpp, aflow.h, kvasp.cpp, aflowlib.h, aflowlib_libraries.cpp, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, pflow_funcs.cpp, pflow.h, pocc.cpp, pocc.h, pocc_old.cpp, pocc_old.h, xatom.cpp, APL/aapl_tensor.cpp, APL/apl_dirphoncalc.cpp, APL/apl.h, APL/apl_phoncalc.cpp, APL/apl_thermalpc.cpp, Makefile)
3.1.195 - 2018/05/21
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.195.tar.xz
        -- beta of release with XZ (SC)
3.1.194 - 2018/05/16
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.194.tar.xz
        -- small bug fixes for g++/gcc 7 and 8 (SC)
        -- preparing for xz compression (SC) 
        194 xatom.cpp  aurostd.h aurostd_main.cpp pthread.cpp README_SCRIPTING README_AFLOW init.cpp aflow.cpp aconvasp_main.cpp
        194c fix aconvasp_main.cpp
	194d avasp.cpp ael_elasticity.cpp ael_get_stress.cpp agl_debye.cpp agl_get_ev.cpp contrib_cormac.cpp contrib_junkai_phasediag.cpp
	194e avasp.cpp xatom.cpp  aurostd_main libraries.cpp ael_elasticity.cpp ael_get_stress.cpp agl_debye.cpp agl_get_ev.cpp contrib_cormac.cpp aconvasp_main.cpp pthread.cpp init.cpp aflow.cpp (heavy)
	194f good bye bzip2 (SC)
	194g aflow_contrib_wahyu.cpp 
	194h aflow.cpp aflow_estructure.cpp aflow_ifrozsl.cpp aflow_ivasp.cpp
	194i aflow_pthreads.cpp
	194j aflow_aflowrc.cpp
	194l BASE64 for pseudopotentials and aflow_data.cpp
	194m aflow_kbin.cpp aflow_kvasp.cpp
	194n aflowlib_web_interface.cpp aflow_matlab_funcs.cpp (EXTRA/MATLAB/plotband.m) aflow_ovasp.cpp aflow_pocc_edos.cpp
	shortened aflow_xatom.cpp working on ZVAL POMASS and removing all EXT stuff. Tests work.  Fixed even further aflow_pthread.cpp aflow_xproto_gus.cpp
3.1.193 - 2018/05/11
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.193.tar.xz
        -- extending xoptions push pop (SC)
3.1.192 - 2018/05/10
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.192.tar.xz
        -- extending xoptions push pop (SC)
3.1.191 - 2018/05/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.191.tar.xz
        -- rationalized orthogonality search in xmatrix (SC)
        aurostd_xmatrix.h aurostd_xmatrix.cpp aflow.h
3.1.190 - 2018/05/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.190.tar.xz
        -- fixed units in AGL output (CT)
        -- fixed permutation vector initialization in APL (required for compiling on DoD CONRAD machine) (DX)
        -- fixed bug in tolerance scan; was only scanning in one direction (DX) 
        -- created faster minimum cartesian distance calculator for skewed cells; fewer duplicate operations (DX) 
        -- fixed conflict between CUT_RAD and CUT_SHELL in AAPL (ME)
        -- properly added directory to bader error output (CO)
        -- fixed slow down with loadEntries() (stringElements2VectorElements()) (CO)
        -- fixed AIMS read-in issue with '#' comments (CO)
        -- added guard around BANDSDATA_JSON() for html generation (CO)
        -- fixed static compile settings (CO)
        -- fixed _sym_op.basis_map_calculated issue in ApplyAtom() for AAPL (CO)
        -- added xStream class for logging workflow updates (CO)
        -- fixed stability criterion vs. (Delta H) fonts in chull report (CO)
        -- substantial clean-up of classes in POCC in anticipation for AVASP_MakeSingleAFLOWIN integration (CO)
        -- added avasp function in preparation for AVASP_MakeSingleAFLOWIN integration (CO)
        -- added try/catches for easy debugging (CO)
        (AUROSTD: main.cpp)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, agl_debye.cpp, avasp.cpp, bader.cpp, chull.cpp, chull.h, contrib_kesong_ipocc.cpp, aflow.cpp, aflow.h, ivasp.cpp, kbin.cpp, kvasp.cpp, aflowlib_web_interface.cpp, pflow.h, pocc.cpp, pocc.h, symmetry.cpp, xatom.cpp, xclasses.cpp, APL/aapl_pairs.cpp, APL/apl_kphonons.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.189 - 2018/05/04 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.189.tar.xz
        Moved these defaults from aflow.h to aflow_aflowrc.cpp for user tuning (SC)
        #define DEFAULT_AFLOW_PRESCRIPT_OUT            string("aflow.prescript.out")  
        #define DEFAULT_AFLOW_PRESCRIPT_COMMAND        string("aflow.prescript.command")  
        #define DEFAULT_AFLOW_POSTSCRIPT_OUT           string("aflow.postscript.out")  
        #define DEFAULT_AFLOW_POSTSCRIPT_COMMAND       string("aflow.postscript.command") 
        #define DEFAULT_AFLOW_PGROUP_OUT               string("aflow.pgroup.out")
        #define DEFAULT_AFLOW_PGROUP_XTAL_OUT          string("aflow.pgroup_xtal.out")
        #define DEFAULT_AFLOW_PGROUPK_OUT              string("aflow.pgroupk.out")
        #define DEFAULT_AFLOW_PGROUPK_XTAL_OUT         string("aflow.pgroupk_xtal.out")
        #define DEFAULT_AFLOW_FGROUP_OUT               string("aflow.fgroup.out")
        #define DEFAULT_AFLOW_SGROUP_OUT               string("aflow.sgroup.out")
        #define DEFAULT_AFLOW_AGROUP_OUT               string("aflow.agroup.out")
        #define DEFAULT_AFLOW_IATOMS_OUT               string("aflow.iatoms.out")
        #define DEFAULT_AFLOW_PGROUP_JSON              string("aflow.pgroup.json")      
        #define DEFAULT_AFLOW_PGROUP_XTAL_JSON         string("aflow.pgroup_xtal.json") 
        #define DEFAULT_AFLOW_PGROUPK_JSON             string("aflow.pgroupk.json")    
        #define DEFAULT_AFLOW_PGROUPK_XTAL_JSON        string("aflow.pgroupk_xtal.json")
        #define DEFAULT_AFLOW_FGROUP_JSON              string("aflow.fgroup.json")   
        #define DEFAULT_AFLOW_SGROUP_JSON              string("aflow.sgroup.json")  
        #define DEFAULT_AFLOW_AGROUP_JSON              string("aflow.agroup.json")    
        #define DEFAULT_AFLOW_IATOMS_JSON              string("aflow.iatoms.json")   
        #define DEFAULT_AFLOW_ICAGES_OUT               string("aflow.icages.out")
        #define DEFAULT_AFLOW_SURFACE_OUT              string("aflow.surface.out")
        #define DEFAULT_AFLOW_QMVASP_OUT               string("aflow.qmvasp.out")
        #define DEFAULT_AFLOW_ERVASP_OUT               string("aflow.error.out")
        #define DEFAULT_AFLOW_IMMISCIBILITY_OUT        string("aflow.immiscibility.out")
        #define DEFAULT_AFLOW_MEMORY_OUT               string("aflow.memory.out")
        #define DEFAULT_AFLOW_FROZSL_INPUT_OUT         string("aflow.frozsl_input.out")
        #define DEFAULT_AFLOW_FROZSL_POSCAR_OUT        string("aflow.frozsl_poscar.out")
        #define DEFAULT_AFLOW_FROZSL_MODES_OUT         string("aflow.frozsl_energies.out")
        #define DEFAULT_AFLOW_FROZSL_EIGEN_OUT         string("aflow.frozsl_eigen.out")
        #define DEFAULT_AFLOW_END_OUT                  string("aflow.end.out")
3.1.188- 2018/05/03 -
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.188.tar.xz
       Fixing kvasp.cpp (SC)
3.1.187- 2018/04/26 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.187.tar.xz
        -- added missing scaling factor information to lattice for AIMS xstructure output (DX) 
        -- fixed bug in edata lattice type/variation determination; accounts for lattices that do not reflect crystal symmetry (DX) 
        -- fixed bug in edata reciprocal lattice type/variation determination; accounts for reciprocal lattices transformed from a lattice that does not reflect crystal symmetry (DX) 
        -- added directory (pwd) information to LDEBUG/ERROR messages for symmetry functions (DX)
        (AFLOW: aconvasp_main.cpp, lattice.cpp, symmetry.cpp, xatom.cpp) 
3.1.186- 2018/04/25 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.186.tar.xz
        -- fixed bug in tolerance scan; was only scanning in one direction (DX) 
        -- created faster minimum cartesian distance calculator for skewed cells; fewer duplicate operations (DX) 
        (AFLOW: symmetry_spacegroup.cpp, symmetry.cpp) 
3.1.185- 2018/04/24 -
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.185.tar.xz
       Fixing kvasp.cpp for dying jobs (SC)
3.1.184- 2018/04/23 -
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.184.tar.xz
       Fixing ivasp.cpp for AFLOW_PSEUDOPOTENTIALS.TXT AFLOW_PSEUDOPOTENTIALS_LIST.TXT (SC)
       Tuning aflow_aflowrc.cpp for HYPERTHREADING in eos,draco,cobra,hydra (SC)
3.1.183- 2018/04/20 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.183.tar.xz
        Fixing ovasp for METAGGA/isKINETIC (SC)
3.1.182- 2018/04/19 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.182.tar.xz
        AFLOW_PSEUDOPOTENTIALS.TXT AFLOW_PSEUDOPOTENTIALS_LIST.TXT (SC)
3.1.181- 2018/04/17 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.181.tar.xz
        Creating defaults in aflow_aflowrc.cpp for HYPERTHREADING in eos,draco,cobra,hydra (SC)
3.1.180- 2018/04/16 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.180.tar.xz
        Creating defaults in aflow_aflowrc.cpp for NCPUS in eos,draco,cobra,hydra (SC)
3.1.179 - 2018/04/13 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.179.tar.xz
        -- fixed PP settings for SCAN (RF)
        -- fixed warning in APIget() (CO)
        -- added relative stability criterion and latex-formatted sg's to chull properties list (CO)
        -- added control.in and geometry.in file name specification in aflowrc (CO)
        -- --generate now applies for aims output as well (CO)
        -- force aflow.in generation even for MODE=AIMS (CO)
        -- fixed k-point mismatch for non-primitive APL runs (CO)
        -- added MINATOMS_RESTRICTED tag for APL (CO)
        -- added functionality for user-defined path in phonon dispersion, specify coords/labels in aflow.in (CO)
        -- swapped out quser for qflow, but maintain backwards compatibility (CO)
        (AUROSTD: xoption.cpp, xoption.cpp, xoption.h, xoption.h, xvector.cpp, xvector.cpp, xvector.h, xvector.h)
        (AFLOW: aconvasp_main.cpp, aconvasp_main.cpp, aflowrc.cpp, aflowrc.cpp, avasp.cpp, avasp.cpp, chull.cpp, chull.cpp, chull.h, chull.h, compare_structure.cpp, compare_structure.cpp, contrib_kesong_hnfcell.cpp, contrib_kesong_hnfcell.cpp, contrib_kesong_ipocc.cpp, contrib_kesong_ipocc.cpp, contrib_kesong_pocc_basic.cpp, contrib_kesong_pocc_basic.cpp, aflow.cpp, aflow.cpp, aflow.h, aflow.h, init.cpp, init.cpp, kbin.cpp, kbin.cpp, kvasp.cpp, kvasp.cpp, aflowlib.h, aflowlib.h, aflowlib_web_interface.cpp, aflowlib_web_interface.cpp, pocc.cpp, pocc.h, poccupation_forcefield.cpp, poccupation_forcefield.cpp, symmetry.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_functions.cpp, symmetry_spacegroup_functions.cpp, symmetry_spacegroup.h, symmetry_spacegroup.h, xatom.cpp, xatom.cpp, xclasses.cpp, xclasses.cpp, APL/aapl_tcond.cpp, APL/aapl_tcond.cpp, APL/aapl_tensor.cpp, APL/aapl_tensor.cpp, APL/apl_dirphoncalc.cpp, APL/apl_dirphoncalc.cpp, APL/apl.h, APL/apl.h, APL/apl_kphonons.cpp, APL/apl_kphonons.cpp, APL/apl_lrphoncalc.cpp, APL/apl_lrphoncalc.cpp, APL/apl_pathbuilder.cpp, APL/apl_pathbuilder.cpp, APL/apl_pdisc.cpp, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/apl_phoncalc.cpp, APL/apl_supercell.cpp, APL/apl_supercell.cpp, APL/qha_eosrun.cpp, APL/qha_eosrun.cpp, Makefile, Makefile, README_AFLOW_APL.TXT, README_AFLOW_APL.TXT, README_AFLOW.TXT, README_AFLOW.TXT)
3.1.178 - 2018/04/13 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.178.tar.xz
        Moved these defaults from aflow.h to aflow_aflowrc.cpp for user tuning (SC)
        #define DEFAULT_FILE_AFLOWLIB_ENTRY_OUT        string("aflowlib.out")
        #define DEFAULT_FILE_AFLOWLIB_ENTRY_JSON       string("aflowlib.json")
        #define DEFAULT_FILE_EDATA_ORIG_OUT            string("edata.orig.out")
        #define DEFAULT_FILE_EDATA_RELAX_OUT           string("edata.relax.out")
        #define DEFAULT_FILE_EDATA_BANDS_OUT           string("edata.bands.out")
        #define DEFAULT_FILE_DATA_ORIG_OUT             string("data.orig.out")
        #define DEFAULT_FILE_DATA_RELAX_OUT            string("data.relax.out")
        #define DEFAULT_FILE_DATA_BANDS_OUT            string("data.bands.out")
        #define DEFAULT_FILE_EDATA_ORIG_JSON           string("edata.orig.json")
        #define DEFAULT_FILE_EDATA_RELAX_JSON          string("edata.relax.json")
        #define DEFAULT_FILE_EDATA_BANDS_JSON          string("edata.bands.json")
        #define DEFAULT_FILE_DATA_ORIG_JSON            string("data.orig.json")
        #define DEFAULT_FILE_DATA_RELAX_JSON           string("data.relax.json")
        #define DEFAULT_FILE_DATA_BANDS_JSON           string("data.bands.json")
        #define DEFAULT_FILE_TIME_OUT                  string("time")
        #define DEFAULT_FILE_SPACEGROUP1_OUT           string("SpaceGroup")
        #define DEFAULT_FILE_SPACEGROUP2_OUT           string("SpaceGroup2")
        #define DEFAULT_FILE_VOLDISTPARAMS_OUT         string("VOLDISTParams")
        #define DEFAULT_FILE_VOLDISTEVOLUTION_OUT      string("VOLDISTEvolution")
3.1.177 - 2018/04/06 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.177.tar.xz
        Bug fixes on autopseudootentialsl (SC)
3.1.176 - 2018/04/06 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.176.tar.xz
        Working PAW_PBE_KIN and PAW_LDA_KIN autopseudopotential (SC)
3.1.175 - 2018/04/06 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.175.tar.xz
        Working on mpcdf-cobra (SC)
3.1.174 - 2018/03/29
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.174.tar.xz    
        -- added bands app to entry page (GG)
        -- fixed scaling factor type in xstructure2json() (DX)
        -- account for Wyckoff positions in represented as fractions (hex/rhl) in minimum enumerated Wyckoff search (DX)
        -- improved minimum enumerated Wyckoff search (consider combinations of origin shifts) (DX)
        -- fixed typo in xmatrix2json function (DX)
        -- limited outlier detection to binaries only, statistics not globally favorable for ternaries and above YET (CO)
        -- added explanatory comments in aflowrc for chull settings (CO)
        -- added statistics output to chull logger (CO)
        -- fixed Greek letter issues with Helvetica fonts in chull (CO)
        -- further decoupled ChullFacet() initialization from ConvexHull() with initialize() and addVertex() (CO)
        -- added user-defined dft_type restrictions in aflowrc (CO)
        -- added LIB1 to loadEntries() default (CO)
        -- fixed stringElements2VectorElements() bugs with LIB1 colons (CO)
        -- added plotting option for iso-max latent heat (CO)
        -- added kJ/mol axis for formation enthalpy hulls (CO)
        -- added logos to hull illustrations, generally fixed header/footers (CO)
        -- added stability criterion analysis to default routine (CO)
        -- fixed sign of decomposition energy/stability criterion (CO)
        -- added equivalent ICSD structures analysis to default routine (CO)
        -- fixed --readme vs. --readme=XX (CO)
        -- fixed superfluous output in --poscar2aflowin (CO)
        -- fixed bug in aconvasp's Ewald summation function as pointed out by Wei Xie (CO)
        -- changed generic xstructure title to include cleannames in case there is pp info (CO)
        (AUROSTD: AUROSTD/aurostd.h, main.cpp, xscalar.h, xvector.cpp)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, chull.cpp, chull.h, init.cpp, aflowlib.h, aflowlib_web_interface.cpp, pflow_funcs.cpp, pflow.h, symmetry_spacegroup.cpp, webapp_bands.js, webapp_entry.js, xatom.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.173 - 2018/03/12
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.173.tar.xz    
        -- Create LIB7/LIB8/LIB9 framework (SC)
3.1.172 - 2018/02/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.172.tar.xz    
        -- Write data in JSON format in AEL/AGL (CT)
        -- Option to turn off VASP symmetry in AEL/AGL (CT)
        -- Write vibrational free energy and vibrational entropy at 300K in aflow.agl.out in AEL/AGL (CT)
        -- Write enthalpy (H = E + pV) in AEL/AGL (CT)
        -- Option to skip pressure/temperature points where no minimum Gibbs free energy is found, instead of truncating pressure/temperature range in AEL/AGL (CT)
        -- fixed tolerance scan counter for lib2raw runs (no longer static variable) (DX)
        -- edata speed increase, only calculate up to pgroup_xtal for lattice, superlattice, and reciprocal lattice (DX)
        -- fixed typo in Wyckoff position "b" for space group #160 (DX)
        -- fixed --kppra command line bug, divided by natoms twice (CO)
        -- added debug output to nanoparticle (CO)
        -- fixed --readme=chull empty string error (CO)
        -- automated plot_unstable with z_filter_cutoff in chull (CO)
        -- fixed plotting ranges with z_filter_cutoff in chull (CO)
        (AFLOW: aconvasp_main.cpp, ael_elastic_fit.cpp, ael_elasticity.cpp, ael_elasticity.h, ael_get_stress.cpp, agl_debye.cpp, agl_debye.h, agl_eqn_state.cpp, agl_get_ev.cpp, agl_hugoniot.cpp, agl_polynomial.cpp, agl_rungibbs.cpp, chull.cpp, data.cpp, aflow.h, lattice.cpp, aflowlib_libraries.cpp, pflow.h, pflow_print.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_ITC_library.cpp, xatom.cpp, Makefile, README_AFLOW_AEL.TXT, README_AFLOW_AGL.TXT, README_AFLOW_SYM.TXT)
3.1.171 - 2018/02/21
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.171.tar.xz    
        -- fixed multiple degeneracy prints in POCC (CO)
        -- added directory information to logging functions in bader/chull (CO)
        -- added --destination=|--path= flag for output of chull (CO)
        -- use full path (pwd) for working directory in AFLOW-SYM functions (DX)
        -- fixed directory flag for AFLOW-SYM functions (DX)
        -- fixed typo in site symmetry of  Wyckoff position "a" for space group #109 (DX)
        (AUROSTD: aurostd.h, main.cpp)
        (AFLOW: aconvasp_main.cpp, bader.cpp, bader.h, chull.cpp, chull.h, poccupation_edos.cpp, symmetry_spacegroup_ITC_library.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.170 - 2018/02/18
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.170.tar.xz    
        -- fixed chmod in LIB2RAW (CO)
        -- fixed DOS-extraction bug (exit() vs. return false) in LIB2RAW for POCC+AEL/AGL runs (CO)
        -- fixed joinWDelimiter() xvector<int> bug: removed erroneous delimiter at the end (CO)
        -- fixed DOSDATA_JSON() + BANDSDATA_JSON() amalgamation: removed erroneous wrapping brackets around DOS (CO)
        -- fixed misleading logging message in APL on primitivization of input (CO)
        (AUROSTD: main.cpp)
        (AFLOW: estructure.cpp, aflow.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp, poccupation_edos.cpp, APL/apl_supercell.cpp, Makefile)
3.1.169 - 2018/02/16
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.169.tar.xz    
        -- added aflow_proto functionality:  --relax_type=IONS, --module=APL, --apl_supercell=3x3x3, --no_volume_adjustment (CO)
        -- NCPUS in APL now respects parent aflow.in (default MAX, then looks at parent aflow.in, otherwise overrides with --np=XX) (CO)
        -- fixed spacegroup bug in apl post hibernation (CO)
        -- added --print and --screen_only options to chull for python integration (CO)
        -- removed extraneous ytick lines from 3D hull (CO)
        -- increased spacing between axes and ticklabels on hulls (CO)
        -- fixed bugs in entropic_temperature hull visualization (CO)
        -- fixed apl2agr/subst make commands (CO)
        -- fixed Niggli tolerance; more robust and use _ZERO_TOL_ (DX)
        -- fixed bug in Niggl step 6 (DX)
        -- fixed comparison of SU(2) to exp(theta*su(2)); parentheses issue (DX)
        -- fixed Quantum Espresso output to include lattice scaling factor (DX)
        -- fixed Quantum Espresso celldm units (Bohr) (DX)
        -- fixed Quantum Espresso alat flag to multiply cell parameters by celldm(1) or A parameter (DX)
        -- fixed tolerance issue with PrintData for lib2raw runs (DX)
        (AFLOW: aconvasp_main.cpp, avasp.cpp, chull.cpp, aflow.h, pflow_print.cpp, symmetry.cpp, xatom.cpp, APL/aapl_tcond.cpp, APL/apl_dm_pdos_save.cpp, APL/apl_doscalc.cpp, APL/apl.h, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/qha_gruneisen.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.168 - 2018/02/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.168.tar.xz
        -- fixed structure comparison bug (overwritten match)  (DX)
        (AFLOW: aflow_compare_structure.h, aflow_compare_structure_function.cpp)
3.1.167 - 2018/02/06
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.167.tar.xz
        -- speed increase for structure comparison (DX)
        -- fixed bugs for structure comparison (DX)
        -- added Wyckoff position analysis to group similar structures in directory comparison (DX)
        -- added grouped Wyckoff position information to the JSON output (DX)
        -- added logger/status for directory comparisons (DX)
        -- added functionality to read .bz2 geometry files (DX)
        -- added --no_scale option for comparisons (required for chull) (CO/DX)
        (AFLOW: aconvasp_main.cpp, chull.cpp, symmetry_spacegroup.cpp, compare_structure.cpp, compare_structure_function.cpp, compare_structure.h) (DX)
3.1.166 - 2018/02/02
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.165.tar.xz
        -- fixed broken --slab commmand (CO)
        (AFLOW: slab.cpp)
3.1.165 - 2018/01/30
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.165.tar.xz
        -- integrated .aflow.rc settings into APL/QHA/AEL/AGL (CO)
        -- added LORBIT=10 to relaxation INCAR (spinD) (CO)
        -- improved magnetic properties extraction, pull from relax first, overwrite with static (CO)
        -- added safety hull coordinates assignment (CO)
        -- new default for load entries (no load xstructures) (CO)
        (AFLOW: aconvasp_main.cpp, ael_get_stress.cpp, agl_get_ev.cpp, chull.cpp, ivasp.cpp, aflowlib_libraries.cpp, APL/apl_phoncalc.cpp, APL/qha_eosrun.cpp, Makefile) (CO)
3.1.164 - 2018/01/25
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.164.tar.xz
        -- fixed minor bugs in Quantum Espresso reader
        -- added ibrav options to Quantum Espresso reader (ibrav2lattice function)
        -- added celldm and a, b, c, cosAB, cosAC, cosBC readers for Quantum Espresso
        -- added Bohr units reader for Quantum Espresso
        (AFLOW: aconvasp_main.cpp, pflow.h, xatom.cpp) (DX)
3.1.163 - 2018/01/23
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.163.tar.xz
        -- added --dist2hull=0.25,0.25 option for chull, which provides the value of the convex hull surface at the specified coordinate/concentration (CO)
        -- added separate readme for chull (--readme=chull) (CO)
        -- customized avasp.cpp (--proto) with options from .aflow.rc, --mpi, and --np (CO)
        -- fixed QHA/AAPL aflow.in keyword bug (CO)
        -- added ANRL modifiers to directory for --proto (Ex: label:ANRL=param1,param2,...) (DX)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, avasp.cpp, chull.cpp, aflow.cpp, aflow.h, init.cpp, kbin.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT, README_AFLOW.TXT) (CO/DX)
3.1.162 - 2018/01/19
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.162.tar.xz
        -- added SU(2) complex matrix and su(2) generator coefficients representation of symmetry elements 
        -- added SU(2) and su(2) information to .out and .json symmetry files        
        -- extended functionality for complex matrices and vectors (trace, exponential, ostream, etc.)
        -- added xcomplex2json to represent complex numbers in json
        -- added pgroupk_xtal (dual of crystal point group operations) standalone function
        -- extended --proto for ANRL prototypes (to generate aflow.in)
        -- updated README_AFLOW_SYM.TXT
        (AUROSTD: aurostd.h, boot.cpp, main.cpp, xmatrix.cpp, xmatrix.h)
        (AFLOW: README_AFLOW_ACONVASP.TXT, README_AFLOW_SYM.TXT, aflow.h, aconvasp_main.cpp, avasp.cpp, symmetry.cpp, xatom.cpp) (DX)
3.1.161 - 2018/01/17
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.161.tar.xz
        -- bool xscalar and xvector gccV7 warning fixes (SC)
        -- double xscalar and xmatrix clang warning fixes (CO)
        -- uint chull gccV7 warning fix (CO)
        (AUROSTD: xscalar.h, xvector.h, xmatrix.h)
        (AFLOW: chull.cpp)
        Added README_AFLOW_SYM.TXT (init.cpp aflow.h aflow.cpp data.cpp) (DX/SC)
3.1.160 - 2018/01/13
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.160.tar.xz
        convex hull mods:
          - keep order of points (m_points) same as input from user/database
          - size will only differ if there are artificial points (they go at the end of m_points)
          - moved point banning to structurePoints() (i.e., m_coord_groups will only contain validated points)
          - only print decomposition phases with non-zero coefficients (tie-line)
          - fixed bug in initial facet.f_outside_set calculation
          - fixed bug in using hull centroid vs. hull reference to align normals
          - fixed bug in vcoords input (shiftlrows)
        (AFLOW: chull.cpp, chull.h)
3.1.159 - 2018/01/10 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.159.tar.xz
        MAJOR UPDATE - convex hull code rewritten, includes: (CO)
          - complete generalization of input (need not be thermodynamic hull, can be ANY data)
          - new major classes for naries, alloys, and coordinate (stoichiometry) groups
          - for hulls with stoichiometric coordinates, the hull and properties are calculated in steps of increasing dimensionality
          - automatic calculation of equivalent ground state structures, near-equivalent (by symmetry) ground state structures, equilibrium phases (by mixture)
          - automatic removal of equivalent ground states for calculation of stability criterion
          - automatic detection/removal of outliers (IQR/MAD)
          - moved all major plotting options to .aflow.rc
          - report any bugs to aflow.org/forum
        added xcombo class for combinations/permutations (CO)
        integrated roundoff/precision/format into stream2stream/utype2tring and variants (xvecDouble2String(),xmatDouble2String()) (CO)
        added another date format for chull (CO)
        added tol to isinteger (CO)
        added +=, -=, gcd (reduced composition), shiftlrows(), isCollinear(), getCentroid(), getGeneralAngles() (similar to Euler angles), getGeneralNormal() (n-dimensional normal), and getQuartiles()/getMAD() (outlier analysis) for xvector (CO)
        added general Householder QR decomposition (CO)
        smarter overloading of loadEntries()/loadXstructures() (CO)
        updated logger() to take directory (aflags) (CO)
        fixed Wyckoff position (a) for space group 117 (DX)
        fixed monoclinic space group determination with multiple unique axis choices (DX)
        (AUROSTD: Makefile,aurostd.cpp,aurostd.h,boot.cpp,main.cpp,xcombos.cpp,xcombos.h,xmatrix.cpp,xmatrix.h,xoption.cpp,xoption.h,xscalar.cpp,xscalar.h,xvector.cpp,xvector.h) (CO, patched SC)
	(AFLOW: Makefile,aflow.h,aconvasp_main.cpp,aflowrc.cpp,anrl.cpp,chull.cpp,chull.h,compare_structure.cpp,compare_structure.h,compare_structure_function.cpp,pflow.h,pflow_funcs.cpp,xatom.cpp) (CO)
	(AFLOWLIB: aflowlib.h,aflowlib_libraries.cpp,aflowlib_web_interface.cpp) (CO)
	(SYM: symmetry_spacegroup.cpp,symmetry_spacegroup_ITC_library.cpp) (CO)
3.1.158 - 2018/01/09 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.158.tar.xz
        Tentative distribution GNU (SC)
3.1.157 - 2017/12/19 -
	Working on mpcdf-draco (SC) 
        Replaced pgroupk with pgroupk_xtal for resolving IBZ in APL/AAPL (CO)
        Fixed pgroupk_xtal flag and wyckoff position for sg133 (DX)
        (aflowrc.cpp,APL/apl_supercell.cpp,APL/apl_mpmesh.cpp,APL/aapl_tcond.cpp,lattice.cpp,symmetry_spacegroup_ITC_library.cpp)
3.1.156 - 2017/12/18 -
        Working on mpcdf-hydra (SC)
3.1.155 - 2017/12/13 -
        updated web component of ANRL CPP files to include the prototype generator and Jmol visualization
        fixed Wyckoff position for SG #171
        improvement for monoclinic space group determination (consider alternative unique axis choices, if possible)
        added more to debug output for the space group routine
        (ANRL/,ANRL_CPPS_20171213/,aconvasp_main.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup_ITC_library.cpp,symmetry_spacegroup_functions.cpp,xatom.cpp) (DX)
      
3.1.154 - 2017/12/07 -
        added pgroupk_xtal (dual of crystal point group operations)
        added magnetic symmetry analysis (crystal + spin) to edata 
        added non-collinear magnetic symmetry analysis (crystal + spin) to 
        edata,sgdata,aflowSG,aflowSYM,fgroup,pgroup_xtal,pgroupk_xtal,agroup,and sgroup
        fixed ApplyAtom (tolerance issue)
        fixed AFLOW-SYM printing bug (--screen_only for Python environment)
        changed point group matrix comparisons (uses Uf; exact)
        added SO(3) generator expansion coefficients onto Lx, Ly, and Lz
        (README_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,ovasp.cpp,pflow.h,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup_functions.cpp,xatom.cpp,lattice.cpp,aflowlib_libraries.cpp) (DX)       
 
3.1.153 - 2017/11/23 -
        Fixing Makefile for GNU version (SC) (fixed Library_ICSD..) (SC)
	Updating directories from AFLOW2 to AFLOW3 (SC)
	Updating directories for findsym platon frozsl (SC)
	
3.1.152 - 2017/11/20 -
        Fixing Makefile for GNU version (SC) (fixed DATA_CVs..)
3.1.151 - 2017/11/17 -
        Fixing Makefile for GNU version (SC) (fixed DATA_CVs..)
3.1.150 - 2017/11/06 -
        NFS cache-cleaning HACK
        (kbin.cpp) (CO)
3.1.149 - 2017/11/06 -
        fixed --xplug
        (libraries.cpp) (SC)
3.1.148 - 2017/10/27 -
        fixed --edata (txt + json print outs)
        (aconvasp_main.cpp,pflow.h,pflow_print.cpp,aflowlib_libraries.cpp) (CO/DX)
3.1.147 - 2017/10/25 -
	added skew test for structures (sym_eps),
	incorporated full symmetry analysis in more functions,
	improved --lib2raw_local,
	improved dos/bands-2json functionality,
	actually fixed --delta_kpoints,
	made min_dist (nearest neighbor distance) function more robust/faster,
	fixed nbondxx units,
	added jsons for sym/structure/bandsdata/dosdata to lib2raw,
	incorporated scale (ReScale()) in more functions
	(aurostd.h,boot.cpp,main.cpp,xmatrix.cpp,xmatrix.h,Makefile,README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,bader.cpp,defects.cpp,estructure.cpp,init.cpp,ovasp.cpp,pflow.h,pflow_print.cpp,surface.cpp,symmetry.cpp,symmetry_spacegroup.cpp,xatom.cpp,xclasses.cpp,xproto.cpp,aflowlib.h,aflowlib_libraries.cpp) (CO/DX)
3.1.146 - 2017/10/16 -
	rewritten bandgap code,
	added extensive FHI-AIMS support (APL),
	modified loadEntries() for GFA code,
	added new xvector operations (*=,/=),
	fixed delta_kpoints function,
	fixed bands/dos-2json functions,
	added local lib2raw
	(aapl_tensor.cpp,apl.h,apl_dirphoncalc.cpp,apl_gsa.cpp,apl_kphonons.cpp,apl_lrphoncalc.cpp,apl_phoncalc.cpp,qha_eosrun.cpp,aurostd.h,boot.cpp,main.cpp,xmatrix.cpp,xvector.cpp,xvector.h,xmatrix.h,Makefile,README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,aflowrc.cpp,chull.cpp,kbin.cpp,ovasp.cpp,pflow.h,poccupation_edos.cpp,xatom.cpp,xclasses.cpp,aflowlib_libraries.cpp,apl_supercell.cpp,iaims.cpp,kaims.cpp,oaims.cpp,estructure.cpp,ivasp.cpp,aflowlib.h) (CO)
3.1.145 - 2017/10/12 -
	fixed tolerance scan range (symmetry.cpp) (DX) 
3.1.144 - 2017/10/03 -
	fixed bandgap type bug,
	fixed precision mistmatch in bands data when plotting band structure,
	changed misleading metric_tensor function name,
	fixed false positive message in structure_comparison,
	fixed APL bugs with reading aflow.in,
	avoid SYMPREC vs. PREC confusion
	(ovasp.cpp,estructure.cpp,xatom.cpp,aconvasp_main.cpp,pflow.h,aflow.h,aurostd.h,aflowlib_libraries.cpp,poccupation_edos.cpp,Makefile,ivasp.cpp,apl_phoncalc.cpp,qha_eosrun.cpp,apl.h,apl_kphonons.cpp,compare_structure.cpp,symmetry_spacegroup_functions.cpp,symmetry_spacegroup.h) (CO/DX)
3.1.143 - 2017/09/28 -
	added magnetic symmetry analysis (crystal + spin) to sgdata,
	aflowSG,aflowSYM,agroup,fgroup,pgroup_xtal,agroup,and sgroup; speed increase for space group routine,
	alter tolerance scan range
	(README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp.cpp,aconvasp_main.cpp,lattice.cpp,pflow.h,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup_functions.cpp,xatom.cpp) (DX)  
3.1.142 - 2017/09/14 -
	added mpcdf_eos machine,
	fixed proper nouns in symmetry output
	(aflow.cpp,kbin.cpp,aflowrc.cpp,init.cpp,kvasp.cpp,README_AFLOW.TXT,Makefile,pflow_print.cpp) (CO/DX)
3.1.141 - 2017/09/12 -
	fixed sgdata/wyccar issue,
	added space group functions to SYM namespace,
	more detailed warning messages
	(symmetry_spacegroup.h,symmetry_spacegroup.cpp,symmetry_spacegroup_function.cpp,symmetry_spacegroup_ITC_library.cpp,aconvasp_main.cpp,pflow_print.cpp,symmetry.cpp) (DX) 
3.1.140 - 2017/09/05 -
	added sgdata function (along with options for json output),
	edata updated (along with options for json output),
	pgroupk speed up,
	added covalent radii info,
	dos/band structure gnuplot updates
	(main.cpp,README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,lattice.cpp,pflow.h,pflow_print.cpp,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup.h,symmetry_spacegroup_functions.cpp,xatom.cpp,estructure.cpp,init.cpp) (DX/CO) 
3.1.139 - 2017/08/29 -
	fixed edata issue (agroups/iatoms were not calculated) for library runs
	(README_AFLOW_ACONVASP.TXT,aconvasp_main.cpp,lattice.cpp,xatom.cpp,aflow.h) (DX)
3.1.138 - 2017/08/25 -
	fully functional AAPL,
	Pearson symbol/lattice type speed up,
	json serializers for symmetry groups,
	tolerance scan implemented for separate symmetry groups,
	chull latex fix
	(chull.cpp,aconvasp_main.cpp,aflow.h,lattice.cpp,pflow.h,symmetry.cpp,xatom.cpp,xclasses.cpp,aapl_pairs.cpp,aapl_tcond.cpp,aapl_tensor.cpp,apl.h,apl_kphonons.cpp,aurostd.h,main.cpp) (JJPR/DX/CO)
3.1.137 - 2017/08/04 -
	incorporated aflow.Xgroup.json,
	json serializers for bands/dos,
	print to screen options for AFLOW-SYM API,
	xvector/xmatrix lrows/urows/lcols/ucols bug fix,
	APL symmetry output append to ofstream,
	removed validation for sgroup (NOT NEEDED)
	(aconvasp_main.cpp,aflow.h,pflow.h,symmetry.cpp,xatom.cpp,README_AFLOW_ACONVASP.TXT,main.cpp,aurostd.h,xvector.cpp,xmatrix.cpp,apl_supercell.cpp,Makefile,estructure.cpp) (DX/EG/CO)
3.1.136 - 2017/08/01 -
	speed up getFullSymBasis,
	fixed some Wyckoff positions,
	fixed printing in apl (spacegroup vs. lattice type)
	(apl_pdisc.cpp,apl_phoncalc.cpp,aflow.h,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup.h,symmetry_spacegroup_ITC_library.cpp,symmetry_spacegroup_functions.cpp) (DX/CO)
3.1.135 - 2017/07/25 -
	fixed xOUTCAR parsing issue with lattice vectors (vasp bug)
	(aflow.h,ovasp.cpp) (CO)
3.1.134 - 2017/07/19 -
	fixed symmetry for library runs (angle/vol tolerance),
	added directory output,
	removed relaxed data from orig structure section in entry page
	(aflowlib_web_interface.cpp,Makefile,aconvasp_main.cpp,aflow.h,lattice.cpp,aflowlib_libraries.cpp,pflow_print.cpp,symmetry.cpp,symmetry_spacegroup.cpp,xatom.cpp,xproto.cpp) (DX/GG/CO)
3.1.133 - 2017/07/19 -
	AEL/AGL: Hugoniots,
	AEL vs pressure (CT) 
3.1.132 - 2017/07/14 -
	added full_sym output to RAW
	(aflowlib_libraries.cpp) (DX/CO)
3.1.131 - 2017/07/14 -
	allow global aflowrc in /etc/aflow.conf plus bug fixes in aflowrc (SC)
3.1.130 - 2017/07/13 -
	bug fixes in aflowrc (SC)
3.1.129 - 2017/07/13 -
	bug fixes init.cpp aflowrc.cpp xvector.cpp (SC/CO)
3.1.128 - 2017/07/12 -
	more ~/.aflowrc definitions and porting to MPI-DE (SC)
3.1.127 - 2017/07/11 -
	more ~/.aflowrc aflow.h ivasp ovasp avasp kbin init VASP_POTCAR_* AFLOW_LIBRARY AFLOW_PROJECT (SC)
3.1.126 - 2017/07/10 -
	more ~/.aflowrc name_dynamics aflow.h  DEFAULT_DOS estructure.cpp poccupation_edos (SC)
3.1.125 - 2017/07/10 -
	more ~/.aflowrc aflow.h  kbin.cpp apl_phoncalc.cpp apl_kphonons.cpp (SC)
3.1.124 - 2017/07/06 -
	quaternion representation additions,
	aflowSYM bug fixes (HfV2)
	(lattice.cpp,xatom.cpp,aflow.h,symmetry.cpp,Makefile) (GG/DX/CO)
3.1.123 - 2017/07/06 -
	adding flags for METAGGA and IVDW (xclasses.cpp) (SC)
3.1.122 - 2017/07/04 -
	aflowrc bug fixes,
	clean names (SC)
3.1.121 - 2017/07/04 -
	aflowrc bug fixes,
	clean names,
	clean METAGGA_SET,
	IVDW_SET (SC)
3.1.120 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp kbin kvasp README* VASP_FORCE* PLATON/FINDSYM (SC)
3.1.119 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp kbin kvasp README* MACHINE DEPENDENT MPI* (SC)
3.1.118 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp kbin kvasp README* GENERIC MPI* (SC)
3.1.117 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp README*.TXT: AFLOW_CORE_TEMPERATURE,VASP_PREC (SC)
3.1.116 - 2017/07/01 -
	starting ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp xoptions.h.cpp avasp.cpp apennsy lattice README_AFLOW_ACONVASP.TXT (SC)
3.1.115 - 2017/06/29 -
	XVASP_INCAR_PREPARE_GENERIC aflow.h ivasp.cpp (SC)
3.1.114 - 2017/06/29 -
	added new Jmol js scripts,
	fixed DISMAG issue with APL,
	checkMEMORY() mostly removed from messaging
	(apl_dirphoncalc.cpp,apl_kphonons.cpp,aflowlib_web_interface.cpp,entry.cpp,entry.js,kvasp.cpp,kbin.cpp,init.cpp,avasp.cpp,aflow.h,README_AFLOW_POCC.TXT,README_AFLOW_APL.TXT,README_AFLOW_AGL.TXT,README_AFLOW_AEL.TXT,README_AFLOW.TXT,Makefile) (CO)
3.1.113 - 2017/06/28 -
	XVASP_INCAR_PREPARE_GENERIC aflow.h ivasp.cpp (SC)
3.1.112 - 2017/06/28 -
	added EDIFFG --ediffg
	aflow.h avasp aconvasp_main ivasp kvasp xclasses README_AFLOW_CONVASP README_AFLOW (SC)
3.1.111 - 2017/06/13 -
	SQLITE integration
	(avasp.cpp ow_chull,aflowlib_web_interface,aflowlib_libraries,kbin,apl_phoncalc,README_AFLOW_VERSIONS_HISTORY.TXT,pflow_print,pflow_funcs,main,aurostd.h,bader,Makefile,aflow,SQLITE) (CO)
3.1.110 - 2017/06/13 -
	fixed chull invalid inputs,
	removed exit(0) in loadLIBX,
	force output for mpi settings in APL,
	fixed settings in CIF writer,
	removed aflow.end.out from CompressDirectory(),
	major fixes to CHGCAR/AECCAR readers for bader analysis (CO)
3.1.109 - 2017/06/06 -
	fixed some chull issues (CO)
3.1.108 - 2017/06/06 -
	individualized QHA/AAPL workflows,
	set reliable convex hull standard to 200 entries (binaries),
	updated ill-calculated systems in the database,
	incorporated sort for loadEntries(),
	fixed DPM=AUTO in APL,
	fixed MCL systems
	(APL/aapl_tensor,apl.h,apl_kphonons,apl_lrphoncalc,apl_phoncalc,qha_eosrun,apl_dirphoncalc,Makefile,README_AFLOW.TXT,README_AFLOW_ACONVASP.TXT,README_AFLOW_AEL.TXT,README_AFLOW_AGL.TXT,README_AFLOW_
APL.TXT,README_AFLOW_POCC.TXT,README_AFLOW_VERSIONS_HISTORY.TXT,aflow,aflow.h,avasp,init,kbin,kvasp,xclasses,aconvasp_main,chull,chull.h,symmetry_spacegroup,symmetry) (CO)
3.1.107 - 2017/06/05 -
	online pubs management curtarolo/oses/toher (SC)
3.1.106 - 2017/06/02 -
	added IVDW
	(aflow.h avasp xclasses aconvasp_main ivasp kvasp) (SC)
3.1.105 - 2017/06/01 -
	added METAGGA
	(aflow.h avasp xclasses aconvasp_main ivasp kvasp) (SC)
3.1.104 - 2017/05/30 -
	set reliable convex hull standard to 200 entries (CO)
3.1.103 - 2017/05/30 -
	added pgroup_xtal propagation in supercell creation,
	also fixed ".nfs*" issues in CompressDirectory() (CO)
3.1.102 - 2017/05/29 -
	issues with symmetry in the aflowlib.out calculation (CO)
3.1.101 - 2017/05/26 -
	cosmetics and added PaddedCENTER (SC)
3.1.100 - 2017/05/26 -
	changed numbering style + cosmetics(SC)
