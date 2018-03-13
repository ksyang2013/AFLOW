3.1.173 - 2018/03/12
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.173.tgz    
        -- Create LIB7/LIB8/LIB9 framework (SC)
3.1.172 - 2018/02/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.172.tgz    
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
        (AFLOW: aconvasp_main.cpp, ael_elastic_fit.cpp, ael_elasticity.cpp, ael_elasticity.h, ael_get_stress.cpp, agl_debye.cpp, agl_debye.h, agl_eqn_state.cpp, agl_get_ev.cpp, agl_hugoniot.cpp, agl_polynomial.cpp, agl_rungibbs.cpp, chull.cpp, data.cpp, aflow.h, lattice.cpp, aflowlib_libraries.cpp, pflow.h, pflow_print.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_ITC_library.cpp, xatom.cpp, Makefile, README_AFLOW_AEL.TXT, README_AFLOW_AGL.TXT, README_AFLOW_SYMMETRY.TXT)

3.1.171 - 2018/02/21
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.171.tgz    
        -- fixed multiple degeneracy prints in POCC (CO)
        -- added directory information to logging functions in bader/chull (CO)
        -- added --destination=|--path= flag for output of chull (CO)
        -- use full path (pwd) for working directory in AFLOW-SYM functions (DX)
        -- fixed directory flag for AFLOW-SYM functions (DX)
        -- fixed typo in site symmetry of  Wyckoff position "a" for space group #109 (DX)
        (AUROSTD: aurostd.h, main.cpp)
        (AFLOW: aconvasp_main.cpp, bader.cpp, bader.h, chull.cpp, chull.h, poccupation_edos.cpp, symmetry_spacegroup_ITC_library.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)

3.1.170 - 2018/02/18
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.170.tgz    
        -- fixed chmod in LIB2RAW (CO)
        -- fixed DOS-extraction bug (exit() vs. return false) in LIB2RAW for POCC+AEL/AGL runs (CO)
        -- fixed joinWDelimiter() xvector<int> bug: removed erroneous delimiter at the end (CO)
        -- fixed DOSDATA_JSON() + BANDSDATA_JSON() amalgamation: removed erroneous wrapping brackets around DOS (CO)
        -- fixed misleading logging message in APL on primitivization of input (CO)
        (AUROSTD: main.cpp)
        (AFLOW: estructure.cpp, aflow.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp, poccupation_edos.cpp, APL/apl_supercell.cpp, Makefile)

3.1.169 - 2018/02/16
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.169.tgz    
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
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.168.tgz
        -- fixed structure comparison bug (overwritten match)  (DX)
        (AFLOW: aflow_compare_structure.h, aflow_compare_structure_function.cpp)

3.1.167 - 2018/02/06
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.167.tgz
        -- speed increase for structure comparison (DX)
        -- fixed bugs for structure comparison (DX)
        -- added Wyckoff position analysis to group similar structures in directory comparison (DX)
        -- added grouped Wyckoff position information to the JSON output (DX)
        -- added logger/status for directory comparisons (DX)
        -- added functionality to read .bz2 geometry files (DX)
        -- added --no_scale option for comparisons (required for chull) (CO/DX)
        (AFLOW: aconvasp_main.cpp, chull.cpp, symmetry_spacegroup.cpp, compare_structure.cpp, compare_structure_function.cpp, compare_structure.h) (DX)

3.1.166 - 2018/02/02
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.165.tgz
        -- fixed broken --slab commmand (CO)
        (AFLOW: slab.cpp)

3.1.165 - 2018/01/30
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.165.tgz
        -- integrated .aflow.rc settings into APL/QHA/AEL/AGL (CO)
        -- added LORBIT=10 to relaxation INCAR (spinD) (CO)
        -- improved magnetic properties extraction, pull from relax first, overwrite with static (CO)
        -- added safety hull coordinates assignment (CO)
        -- new default for load entries (no load xstructures) (CO)
        (AFLOW: aconvasp_main.cpp, ael_get_stress.cpp, agl_get_ev.cpp, chull.cpp, ivasp.cpp, aflowlib_libraries.cpp, APL/apl_phoncalc.cpp, APL/qha_eosrun.cpp, Makefile) (CO)

3.1.164 - 2018/01/25
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.164.tgz
        -- fixed minor bugs in Quantum Espresso reader
        -- added ibrav options to Quantum Espresso reader (ibrav2lattice function)
        -- added celldm and a, b, c, cosAB, cosAC, cosBC readers for Quantum Espresso
        -- added Bohr units reader for Quantum Espresso
        (AFLOW: aconvasp_main.cpp, pflow.h, xatom.cpp) (DX)

3.1.163 - 2018/01/23
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.163.tgz
        -- added --dist2hull=0.25,0.25 option for chull, which provides the value of the convex hull surface at the specified coordinate/concentration (CO)
        -- added separate readme for chull (--readme=chull) (CO)
        -- customized avasp.cpp (--proto) with options from .aflow.rc, --mpi, and --np (CO)
        -- fixed QHA/AAPL aflow.in keyword bug (CO)
        -- added ANRL modifiers to directory for --proto (Ex: label:ANRL=param1,param2,...) (DX)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, avasp.cpp, chull.cpp, aflow.cpp, aflow.h, init.cpp, kbin.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT, README_AFLOW.TXT) (CO/DX)

3.1.162 - 2018/01/19
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.162.tgz
        -- added SU(2) complex matrix and su(2) generator coefficients representation of symmetry elements 
        -- added SU(2) and su(2) information to .out and .json symmetry files        
        -- extended functionality for complex matrices and vectors (trace, exponential, ostream, etc.)
        -- added xcomplex2json to represent complex numbers in json
        -- added pgroupk_xtal (dual of crystal point group operations) standalone function
        -- extended --proto for ANRL prototypes (to generate aflow.in)
        -- updated README_AFLOW_SYMMETRY.TXT
        (AUROSTD: aurostd.h, boot.cpp, main.cpp, xmatrix.cpp, xmatrix.h)
        (AFLOW: README_AFLOW_ACONVASP.TXT, README_AFLOW_SYMMETRY.TXT, aflow.h, aconvasp_main.cpp, avasp.cpp, symmetry.cpp, xatom.cpp) (DX)

3.1.161 - 2018/01/17
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.161.tgz
        -- bool xscalar and xvector gccV7 warning fixes (SC)
        -- double xscalar and xmatrix clang warning fixes (CO)
        -- uint chull gccV7 warning fix (CO)
        (AUROSTD: xscalar.h, xvector.h, xmatrix.h)
        (AFLOW: chull.cpp)
        Added README_AFLOW_SYMMETRY.TXT (init.cpp aflow.h aflow.cpp data.cpp) (DX/SC)

3.1.160 - 2018/01/13
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.160.tgz
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
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.159.tgz
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
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.158.tgz
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
