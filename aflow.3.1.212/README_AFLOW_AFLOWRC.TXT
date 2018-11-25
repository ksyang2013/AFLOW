****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                     materials.duke.edu/AFLOW/README_AFLOW_AFLOWRC.TXT
****************************************************************************************************

AFLOWRC README:

README written by: Corey Oses (corey.oses@duke.edu)

****************************************************************************************************
GENERAL OVERVIEW

The .aflow.rc file is a general configuration file for AFLOW. 
Included in the file are machine-specific paths (VASP executables, pseudopotential files, 
database locations) and general run defaults.
Many of these settings can be overridden by tags in the aflow.in, and even by command line
flags when running AFLOW.
The general priority is as follows: command-line flags, aflow.in tags, and .aflow.rc settings.
The .aflow.rc file is located in your $HOME directory (~/.aflow.rc).
In the event that the file does not exist, AFLOW will generate a new .aflow.rc.
It is a good practice to keep the .aflow.rc file up to date with every new version.
Consider moving ~/.aflow.rc to ~/.aflow.rc_OLD and melding your system modifications.

****************************************************************************************************
LIST OF AFLOWRC SETTINGS:

AFLOWRC
  description: Indicates the version of AFLOW used to generate the .aflow.rc.
  type: string
  default: "X.X.XXX"

####################################################################################################
DEFAULT DEFINITIONS

DEFAULT_KZIP_BIN
  description: Compression format executable for files after the run has completed.
  related aflow.in tag: [AFLOW_MODE_ZIP=xz]
  type: string
  default: "xz"

DEFAULT_KZIP_EXT
  description: Compression format extension for files after the run has completed.
  type: string
  default: ".xz"

####################################################################################################
FILENAMES FOR AFLOW.ORG ANALYSIS

DEFAULT_FILE_AFLOWLIB_ENTRY_OUT
  description: File name for the entry summary properties in plain text (--lib2raw).
  type: string
  default: "aflowlib.out"

DEFAULT_FILE_AFLOWLIB_ENTRY_JSON
  description: File name for the entry summary properties in JSON format (--lib2raw).  
  type: string
  default: "aflowlib.json"

DEFAULT_FILE_EDATA_ORIG_OUT
  description: File name for the extended crystallographic data of the original structure in plain text (--edata).
  type: string
  default: "edata.orig.out"

DEFAULT_FILE_EDATA_RELAX_OUT
  description: File name for the extended crystallographic data of the relaxed structure in plain text (--edata).
  type: string
  default: "edata.relax.out"

DEFAULT_FILE_EDATA_BANDS_OUT
  description: File name for the extended crystallographic data of the primitivized relaxed structure in plain text (--edata).
  type: string
  default: "edata.bands.out"

DEFAULT_FILE_DATA_ORIG_OUT
  description: File name for the crystallographic data of the original structure in plain text (--data).
  type: string
  default: "data.orig.out"

DEFAULT_FILE_DATA_RELAX_OUT
  description: File name for the crystallographic data of the relaxed structure in plain text (--data).
  type: string
  default: "data.relax.out"

DEFAULT_FILE_DATA_BANDS_OUT
  description: File name for the crystallographic data of the primitivized relaxed structure in plain text (--data).
  type: string
  default: "data.bands.out"

DEFAULT_FILE_EDATA_ORIG_JSON
  description: File name for the extended crystallographic data of the original structure in JSON format (--edata --print=json).
  type: string
  default: "edata.orig.json"

DEFAULT_FILE_EDATA_RELAX_JSON
  description: File name for the extended crystallographic data of the relaxed structure in JSON format (--edata --print=json).
  type: string
  default: "edata.relax.json"

DEFAULT_FILE_EDATA_BANDS_JSON
  description: File name for the extended crystallographic data of the primitivized relaxed structure in JSON format (--edata --print=json).
  type: string
  default: "edata.bands.json"

DEFAULT_FILE_DATA_ORIG_JSON
  description: File name for the crystallographic data of the original structure in plain text (--data --print=json).
  type: string
  default: "data.orig.json"

DEFAULT_FILE_DATA_RELAX_JSON
  description: File name for the crystallographic data of the relaxed structure in plain text (--data --print=json).
  type: string
  default: "data.relax.json"

DEFAULT_FILE_DATA_BANDS_JSON
  description: File name for the crystallographic data of the primitivized relaxed structure in plain text (--data --print=json).
  type: string
  default: "data.bands.json"

DEFAULT_FILE_TIME_OUT
  description: File name for the entry calculation time (--lib2raw).
  type: string
  default: "time"

DEFAULT_FILE_SPACEGROUP1_OUT
  description: File name for the entry space group pre/mid/post relaxation (--lib2raw --aflowSG=loose).
  type: string
  default: "SpaceGroup"

DEFAULT_FILE_SPACEGROUP2_OUT
  description: File name for the entry space group pre/mid/post relaxation (--lib2raw --aflowSG=tight).
  type: string
  default: "SpaceGroup2"

DEFAULT_FILE_VOLDISTPARAMS_OUT
  description: File name for the entry volume/nearest-neighbor-distances post relaxation (--lib2raw).
  type: string
  default: "VOLDISTParams"

DEFAULT_FILE_VOLDISTEVOLUTION_OUT
  description: File name for the entry volume/nearest-neighbor-distances pre/mid/post relaxation (--lib2raw).
  type: string
  default: "VOLDISTEvolution"

####################################################################################################
FILENAMES FOR AFLOW OPERATION

DEFAULT_AFLOW_PRESCRIPT_OUT
  description: File name for the output of the script to be run before AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.prescript.out"

DEFAULT_AFLOW_PRESCRIPT_COMMAND
  description: Script to be run before AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.prescript.command"

DEFAULT_AFLOW_POSTSCRIPT_OUT
  description: File name for the output of the script to be run after AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.postscript.out"

DEFAULT_AFLOW_POSTSCRIPT_COMMAND
  description: Script to be run after AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.postscript.command"

DEFAULT_AFLOW_PGROUP_OUT
  description: File name for the point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroup.out"

DEFAULT_AFLOW_PGROUP_JSON
  description: File name for the point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroup.json"

DEFAULT_AFLOW_PGROUP_XTAL_OUT
  description: File name for the crystallographic point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroup_xtal.out"

DEFAULT_AFLOW_PGROUP_XTAL_JSON
  description: File name for the crystallographic point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroup_xtal.json"

DEFAULT_AFLOW_PGROUPK_OUT
  description: File name for the reciprocal point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroupk.out"

DEFAULT_AFLOW_PGROUPK_JSON
  description: File name for the reciprocal point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroupk.json"

DEFAULT_AFLOW_PGROUPK_XTAL_OUT
  description: File name for the dual of the crystallographic point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroupk_xtal.out"

DEFAULT_AFLOW_PGROUPK_XTAL_JSON
  description: File name for the dual of the crystallographic point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroupk_xtal.json"

DEFAULT_AFLOW_FGROUP_OUT
  description: File name for the factor group (coset representative) operations in plain text (--aflow-sym).
  type: string
  default: "aflow.fgroup.out"

DEFAULT_AFLOW_FGROUP_JSON
  description: File name for the factor group (coset representative) operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.fgroup.json"

DEFAULT_AFLOW_SGROUP_OUT
  description: File name for the full space group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.sgroup.out"

DEFAULT_AFLOW_SGROUP_JSON
  description: File name for the full space group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.sgroup.json"

DEFAULT_AFLOW_AGROUP_OUT
  description: File name for the site symmetry operations in plain text (--aflow-sym).
  type: string
  default: "aflow.agroup.out"

DEFAULT_AFLOW_AGROUP_JSON
  description: File name for the site symmetry operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.agroup.json"

DEFAULT_AFLOW_IATOMS_OUT
  description: File name for the inequivalent atoms calculation in plain text (--aflow-sym).
  type: string
  default: "aflow.iatoms.out"

DEFAULT_AFLOW_IATOMS_JSON
  description: File name for the inequivalent atoms calculation in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.iatoms.json"

DEFAULT_AFLOW_ICAGES_OUT
  description: File name for the interstitial defects (cages) calculation in plain text (--cages).
  type: string
  default: "aflow.icages.out"

DEFAULT_AFLOW_SURFACE_OUT
  description: File name for the surface (hkl) calculation in plain text (--hkl_search).
  type: string
  default: "aflow.surface.out"

DEFAULT_AFLOW_QMVASP_OUT
  description: File name for the VASP calculation summary in plain text (--qmvasp).
  type: string
  default: "aflow.qmvasp.out"

DEFAULT_AFLOW_ERVASP_OUT
  description: File name for the VASP calculation errors in plain text.
  type: string
  default: "aflow.error.out"

DEFAULT_AFLOW_IMMISCIBILITY_OUT
  description: File name for the miscibility determination in plain text (--miscibility).
  type: string
  default: "aflow.immiscibility.out"

DEFAULT_AFLOW_MEMORY_OUT
  description: File name for the memory check in plain text.
  type: string
  default: "aflow.memory.out"

DEFAULT_AFLOW_FROZSL_INPUT_OUT
  description: File name for the input settings of the FROZSL calculation.
  type: string
  default: "aflow.frozsl_input.out"

DEFAULT_AFLOW_FROZSL_POSCAR_OUT
  description: File name for the structure used in the FROZSL calculation.
  type: string
  default: "aflow.frozsl_poscar.out"

DEFAULT_AFLOW_FROZSL_MODES_OUT
  description: File name for the phonon eigen-energies from the FROZSL calculation.
  type: string
  default: "aflow.frozsl_energies.out"

DEFAULT_AFLOW_FROZSL_EIGEN_OUT
  description: File name for the phonon eigen-modes from the FROZSL calculation.
  type: string
  default: "aflow.frozsl_eigen.out"

DEFAULT_AFLOW_END_OUT
  description: File name for the output file indicating AFLOW is finished.
  type: string
  default: "aflow.end.out"

####################################################################################################
DEFAULT GENERIC MPI

MPI_START_DEFAULT
  description: Commands to run before running the MPI program.
  type: string
  default: ""

MPI_STOP_DEFAULT
  description: Commands to run after running the MPI program.
  type: string
  default: ""

MPI_COMMAND_DEFAULT
  description: MPI program command.
  type: string
  default: "mpirun -np"

MPI_NCPUS_DEFAULT
  description: Number of cores to run with the MPI program.
  type: integer
  default: 4

MPI_NCPUS_MAX
  description: Maximum number of cores that AFLOW can consume.
  type: integer
  default: 64

####################################################################################################
DEFAULTS BINARY

DEFAULT_VASP_GAMMA_BIN
  description: Name for serial Gamma-Point-only VASP executable.
  type: string
  default: "vasp46s_gamma"

DEFAULT_VASP_GAMMA_MPI_BIN
  description: Name for MPI Gamma-Point-only VASP executable.
  type: string
  default: "mpivasp46s_gamma"

DEFAULT_VASP_BIN
  description: Name for serial VASP executable.
  type: string
  default: "vasp46s"

DEFAULT_VASP_MPI_BIN
  description: Name for MPI VASP executable.
  type: string
  default: "mpivasp46s"

DEFAULT_VASP5_BIN
  description: Name for serial VASP 5 executable.
  type: string
  default: "vasp54s"

DEFAULT_VASP5_MPI_BIN
  description: Name for MPI VASP 5 executable.
  type: string
  default: "mpivasp54s"

DEFAULT_AIMS_BIN
  description: Name for FHI-AIMS executable.
  type: string
  default: "aims"

####################################################################################################
DEFAULTS POTCARS

DEFAULT_VASP_POTCAR_DIRECTORIES
  description: Path(s) to search for the VASP pseudopotential files, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "/common/VASP,/common/AFLOW/VASP,/home/aflow/common/AFLOW/VASP,/fslhome/fslcollab8/group/VASP,/fslhome/glh43/src/,/share/home/00470/tg457283/common/AFLOW/VASP/,/share/home/00457/tg457357/common/AFLOW/VASP/,/home/mehl/bin/AFLOW/VASP/,~/common/VASP/,~/common/AFLOW/VASP/,/nics/a/proj/aflow/common/AFLOW/VASP/,/home/users/aflow/common/VASP,/share/apps/AFLOW3/VASP,/projects/kyang-group/common/VASP,/somewhere/"

DEFAULT_VASP_POTCAR_DATE
  description: Date component of VASP pseudopotential files path, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "current"

DEFAULT_VASP_POTCAR_SUFFIX
  description: Last component of VASP pseudopotential files path, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "/POTCAR"

DEFAULT_VASP_POTCAR_DATE_POT_LDA
  description: Date component of VASP LDA pseudopotential files to be added to the --aflow_proto directory, e.g., V_pv:LDA:01Apr2000.
  type: string
  default: "01Apr2000"

DEFAULT_VASP_POTCAR_DATE_POT_GGA
  description: Date component of VASP GGA pseudopotential files to be added to the --aflow_proto directory, e.g., V_pv:LDA:01Apr2000.
  type: string
  default: "01Apr2000"

DEFAULT_VASP_POTCAR_DIR_POT_LDA
  description: Path component of VASP LDA pseudopotential files, e.g., /common/VASP/pot_LDA/current/Mn_pv/POTCAR.
  type: string
  default: "pot_LDA"

DEFAULT_VASP_POTCAR_DIR_POT_GGA
  description: Path component of VASP GGA pseudopotential files, e.g., /common/VASP/pot_GGA/current/Mn/POTCAR.
  type: string
  default: "pot_GGA"

DEFAULT_VASP_POTCAR_DIR_POT_PBE
  description: Path component of VASP PBE pseudopotential files, e.g., /common/VASP/pot_PBE/current/Mn/POTCAR.
  type: string
  default: "pot_PBE"

DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA
  description: Path component of VASP PAW LDA pseudopotential files, e.g., /common/VASP/potpaw_LDA/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_LDA"

DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA
  description: Path component of VASP PAW GGA pseudopotential files, e.g., /common/VASP/potpaw_GGA/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_GGA"

DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE
  description: Path component of VASP PAW PBE pseudopotential files, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_PBE"

DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN
  description: Path component of VASP PAW LDA pseudopotential files containing information on the kinetic energy density of the core-electrons, e.g., /common/VASP/potpaw_PBE.54/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_LDA.54"

DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN
  description: Path component of VASP PAW PBE pseudopotential files containing information on the kinetic energy density of the core-electrons, e.g., /common/VASP/potpaw_PBE.54/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_PBE.54"

####################################################################################################
DEFAULTS KPOINTS/DOS

DEFAULT_BANDS_GRID
  description: Number of points in which the energy dispersion is resolved between each leg of the high-symmetry path.
  type: integer
  default: 20

DEFAULT_BANDS_LATTICE
  description: Crystal lattice of the structure specifying the high-symmetry path.
  type: string
  choices: "AUTO", "BCC", "BCT1", "BCT2", "CUB", "FCC", "HEX", "MCL", "MCLC1", "MCLC2", "MCLC3", "MCLC4", "MCLC5", "ORC", "ORCC", "ORCF1", "ORCF2", "ORCF3", "ORCI", "RHL1", "RHL2", "TET", "TRI1a", "TRI1b", "TRI2a", "TRI2b"
  default: "AUTO"
  doi: 10.1016/j.commatsci.2010.05.010

DEFAULT_KSCHEME
  description: K-points mesh scheme.
  type: string
  choices: "M", "G"
  default: "M"

DEFAULT_KPPRA
  description: Number of k-points per reciprocal atom.
  type: integer
  default: 6000

DEFAULT_STATIC_KSCHEME
  description: K-points mesh scheme for the static (energy/DOS) run.
  type: string
  choices: "M", "G"
  default: "M"

DEFAULT_KPPRA_STATIC
  description: Number of k-points per reciprocal atom for the static (energy/DOS) run.
  type: integer
  default: 10000

DEFAULT_KPPRA_ICSD
  description: Number of k-points per reciprocal atom for structures from the ICSD catalog.
  type: integer
  default: 8000

DEFAULT_UNARY_BANDS_GRID
  description: Number of points in which the energy dispersion is resolved between each leg of the high-symmetry path for unary structures.
  type: integer
  default: 128

DEFAULT_UNARY_KPPRA
  description: Number of k-points per reciprocal atom for unary structures.
  type: integer
  default: 8000

DEFAULT_UNARY_KPPRA_STATIC
  description: Number of k-points per reciprocal atom for unary structures and the static (energy/DOS) run..
  type: integer
  default: 8000

DEFAULT_PHONONS_KSCHEME
  description: K-points mesh scheme for phonon calculations.
  type: string
  choices: "M", "G"
  default: "G"

DEFAULT_DOS_EMIN
  description: Minimum of the energy range over which the DOS is plotted.
  type: integer
  default: -10

DEFAULT_DOS_EMAX
  description: Maximum of the energy range over which the DOS is plotted.
  type: integer
  default: 10

DEFAULT_DOS_SCALE
  description: Interval for the energy axis over which the DOS is plotted.
  type: double
  default: 1.2

####################################################################################################
DEFAULTS PRECISION

DEFAULT_VASP_PREC_ENMAX_LOW
  description: Low-level (default) multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1

DEFAULT_VASP_PREC_ENMAX_MEDIUM
  description: Medium-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.3

DEFAULT_VASP_PREC_ENMAX_NORMAL
  description: Normal-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.3

DEFAULT_VASP_PREC_ENMAX_HIGH
  description: High-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.4

DEFAULT_VASP_PREC_ENMAX_ACCURATE
  description: Accurate-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.4

DEFAULT_VASP_SPIN_REMOVE_CUTOFF
  description: Magnetic moment threshold below which the spin-polarization setting will be turned off.
  type: double
  default: 0.05

DEFAULT_VASP_PREC_POTIM
  description: Sets starting value for POTIM in VASP (scaling constant for the forces). AFLOW may adjust as errors occur.
  type: double
  default: 0.5

DEFAULT_VASP_PREC_EDIFFG
  description: Sets EDIFFG for VASP. Note that a positive value defines the break condition for the ionic relaxation loop, while a negative value defines the break condition based on the forces.
  type: double
  default: -0.001

####################################################################################################
DEFAULTS OPTIONS

DEFAULT_VASP_EXTERNAL_INCAR
  description: File to load in as default VASP INCAR.
  type: string
  default: ./INCAR

DEFAULT_VASP_EXTERNAL_POSCAR
  description: File to load in as default VASP POSCAR.
  type: string
  default: ./POSCAR

DEFAULT_VASP_EXTERNAL_POTCAR
  description: File to load in as default VASP POTCAR.
  type: string
  default: ./POTCAR

DEFAULT_VASP_EXTERNAL_KPOINTS
  description: File to load in as default VASP KPOINTS.
  type: string
  default: ./KPOINTS

DEFAULT_AIMS_EXTERNAL_CONTROL
  description: File to load in as default AIMS control.in.
  type: string
  default: ./control.in

DEFAULT_AIMS_EXTERNAL_GEOM
  description: File to load in as default AIMS geom.in.
  type: string
  default: ./geometry.in

DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE
  description: Default pseudopotential type.
  type: string
  default: potpaw_PBE

DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME
  description: Default relax mode setting.
  type: string
  choices: ENERGY, FORCES, ENERGY_FORCES, FORCES_ENERGY
  default: ENERGY

DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME
  description: Default VASP PREC setting.
  type: string
  choices: LOW, MEDIUM, NORMAL, HIGH, ACCURATE
  default: ACCURATE

DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME
  description: Default VASP ALGO setting.
  type: string
  choices: NORMAL, VERYFAST, FAST, ALL, DAMPED
  default: NORMAL

DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME
  description: Default VASP METAGGA setting.
  type: string
  choices: NONE, TPSS, RTPSS, M06L, MS0, MS1 MS2, MBJ, SCAN
  default: NONE

DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME
  description: Default VASP IVDW setting.
  type: bool
  choices: 0, 1, 10, 11, 12, 2, 20, 21, 202, 4
  default: 0

DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME
  description: Default VASP integration scheme.
  type: string
  choices: DEFAULT, METAL, INSULATOR, SEMICONDUCTOR
  default: DEFAULT

DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME
  description: Default VASP Kerker mixing scheme settings. Defines AMIX/BMIX/AMIX_MAG/BMIX_MAG.
  type: string
  choices: DEFAULT (unspecified), AUTO, US, PAW
  default: DEFAULT

DEFAULT_VASP_FORCE_OPTION_SYM
  description: Default VASP ISYM option.
  type: integer
  choices: 0, 1, 2
  default: 1

DEFAULT_VASP_FORCE_OPTION_SPIN
  description: Adapts INCAR for VASP spin-polarized calculations.
  type: bool
  default: 1

DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1
  description: Turns off spin-polarized calculation after first relaxation if the magnetic moment is less than VASP_SPIN_REMOVE_CUTOFF.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2
  description: Turns off spin-polarized calculation after second relaxation if the magnetic moment is less than VASP_SPIN_REMOVE_CUTOFF.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_BADER
  description: Adapts INCAR for VASP to retrieve core charge (AECCARx) for Bader partial charge calculation.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_ELF
  description: Adapts INCAR for VASP to calculate the electron localization function.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM
  description: Adapts INCAR for VASP to initialize magnetic moment for all species. It defaults to a ferromagnetic solution.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_WAVECAR
  description: Adapts INCAR for VASP to save the WAVECAR.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_CHGCAR
  description: Adapts INCAR for VASP to save the CHGCAR.
  type: bool
  default: 1

DEFAULT_VASP_FORCE_OPTION_LSCOUPLING
  description: Adapts INCAR for VASP spin-orbit coupling calculations.
  type: bool
  default: 0

####################################################################################################
AFLOW_LIBRARY AFLOW_PROJECT

DEFAULT_AFLOW_LIBRARY_DIRECTORIES
  description: Path(s) to search for the library of materials calculations.
  type: string
  default: "/common/AFLOW/LIBS/,/home/aflow/common/AFLOW/LIBS/,/fslhome/glh43/src/,/usr/local/bin/,/fslhome/fslcollab8/group/bin/,/home/auro/work/AFLOW3/,~/common/AFLOW/LIBS/,./,/nics/a/proj/aflow/common/AFLOW/LIBS/,/home/users/aflow/common/AFLOW/LIBS,/home/junkai/PROTO_DATABASE/,/projects/kyang-group/common/LIBS,/somewhere/"

DEFAULT_AFLOW_PROJECTS_DIRECTORIES
  description: Path(s) of the various library projects.
  type: string
  default: "/common/ICSD,/common/LIB1,/common/LIB2,/common/LIB3,/common/LIB4,/common/LIB5,/common/LIB6,/common/LIB7,/common/LIB8,/common/LIB9"

####################################################################################################
DEFAULT PLATON/FINDSYM

DEFAULT_PLATON_P_EQUAL
  description: Sets PLATON EQUAL tag (incorporate atom types).
  type: bool
  default: 0

DEFAULT_PLATON_P_EXACT
  description: Sets PLATON EXACT tag (match all atoms).
  type: bool
  default: 0

DEFAULT_PLATON_P_ANG
  description: Sets PLATON ang (angle criterion) for metrical symmetry of lattice.
  type: bool
  default: 1

DEFAULT_PLATON_P_D1
  description: Sets PLATON d1 (distance criterion) for non-inversion symmetry.
  type: double
  default: 0.25

DEFAULT_PLATON_P_D2
  description: Sets PLATON d2 (distance criterion) for inversion symmetry.
  type: double
  default: 0.25

DEFAULT_PLATON_P_D3
  description: Sets PLATON d3 (distance criterion) for translation symmetry.
  type: double
  default: 0.25

DEFAULT_FINDSYM_TOL
  description: Sets FINDSYM tolerance.
  type: double
  default: 0.001

####################################################################################################
DEFAULTS GNUPLOT

DEFAULT_GNUPLOT_EPS_FONT
  description: Sets gnuplot eps terminal font.
  type: string
  default: "Helvetica"

DEFAULT_GNUPLOT_EPS_FONT_BOLD
  description: Sets gnuplot eps terminal bold font.
  type: string
  default: "Helvetica-Bold"

DEFAULT_GNUPLOT_EPS_FONT_ITALICS
  description: Sets gnuplot eps terminal italicized font.
  type: string
  default: "Helvetica-Oblique"

DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS
  description: Sets gnuplot eps terminal bold italicized font.
  type: string
  default: "Helvetica-BoldOblique"

DEFAULT_GNUPLOT_PNG_FONT
  description: Sets gnuplot png terminal font.
  type: string
  default: "Arial"

DEFAULT_GNUPLOT_PNG_FONT_BOLD
  description: Sets gnuplot png terminal bold font.
  type: string
  default: "Arial_Bold"

DEFAULT_GNUPLOT_PNG_FONT_ITALICS
  description: Sets gnuplot png terminal italicized font.
  type: string
  default: "Arial_Italic"

DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS
  description: Sets gnuplot png terminal bold italicized font.
  type: string
  default: "Arial_BoldItalic"

####################################################################################################
DEFAULTS CHULL

DEFAULT_CHULL_ALLOWED_DFT_TYPES
  description: Defines the allowed entries based on the DFT calculation type.
  type: string
  choices: "US", "GGA", "PAW_LDA", "PAW_GGA", "PAW_PBE", "GW", "HSE06"
  default: "PAW_PBE"

DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES
  description: Allows all entries independent of DFT calculation type.
  type: bool
  default: 0

DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES
  description: Defines the minimum number of entries for a reliable binary hull.
  type: integer
  default: 200

DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS
  description: Enables determination of outliers.
  type: bool
  default: 1

DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES
  description: Defines the minimum number of entries for a reliable outlier analysis. Only phases stable with respect to their end-members are considered for the outlier analysis (below the zero H_f tie-line).
  type: integer
  default: 50

DEFAULT_CHULL_OUTLIER_MULTIPLIER
  description: Defines the bounds beyond the interquartile range for which points are considered outliers.
  type: double
  default: 3.25

DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED
  description: AFLOW maintains a list of (older) prototypes known to have converged poorly. These entries are likely outliers. If this flag is on (1), then the entries are removed automatically. Turning this flag off (0) is not recommended.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_BANNER
  description: Defines banner type.
  type: integer
  choices: 0 (no banner), 1 (normal banner), 2 (small banner)
  default: 1

DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN
  description: Adds an additional column to the PDF table output with the entry compound label.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_STOICH_HEADER
  description: Toggles between a reduced composition header (default) or a stoichiometry header for the entries table.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_PLOT_UNARIES
  description: Incorporates the end-members in the convex hull illustration.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_PLOT_OFF_HULL
  description: Incorporates unstable phases in the convex hull illustration, but excludes phases unstable with respect to their end-members (above the zero H_f tie-line).
  type: integer
  choices: -1 (default: true for 2-dimensional systems, false for 3-dimensional systems), 0 (false), 1 (true).
  default: -1

DEFAULT_CHULL_LATEX_PLOT_UNSTABLE
  description: Incorporates all unstable phases in the convex hull illustration.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_FILTER_SCHEME
  description: Defines exclusion scheme for the convex hull illustration. In contrast to --neglect, this scheme is limited only to the illustration, points are still included in the analysis/report. The following strings are accepted: Z-axis (also Energy-axis) or Distance. Z-axis refers to a scheme that eliminates structures from the illustration based on their formation enthalpy. On the other hand, Distance refers to a scheme that eliminates structures from the illustration based on their distances from the hull. The criteria (value) for elimination is defined by DEFAULT_CHULL_LATEX_FILTER_VALUE.
  type: string
  choices: "z-axis", "distance"
  default: ""

DEFAULT_CHULL_LATEX_FILTER_VALUE
  description: Defines the value beyond which points are excluded per the scheme defined with DEFAULT_CHULL_LATEX_FILTER_SCHEME.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_COLOR_BAR
  description: Defines whether to show the color bar graphic. Colors can still be incorporated without the color bar graphic.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_HEAT_MAP
  description: Defines whether to color the facets with heat maps illustrating their depth (3-dimensional illustration only).
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_COLOR_GRADIENT
  description: Defines whether to incorporate a color scheme at all in the illustration. Turning this flag off will also turn off DEFAULT_CHULL_LATEX_COLOR_BAR and DEFAULT_CHULL_LATEX_HEAT_MAP.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_COLOR_MAP
  description: Defines the color map, options are presented in the PGFPLOTS Package Manual. Default is "rgb(0pt)=(0.035,0.270,0.809); rgb(63pt)=(1,0.644,0)".
  type: string
  default: ""

DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR
  description: Defines the color used for ternary compounds. Default is white.
  type: string
  default: ""

DEFAULT_CHULL_LATEX_REVERSE_AXIS
  description: Flips the axes of binary convex hulls.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW
  description: Adds a drop shadow underneath facet lines.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_LINKS
  description: Defines the links scheme. True/false, i.e., 0/1, will toggle all links on/off. 2 enables external hyperlinks only (no links to other sections of the PDF). 3 enables internal links only (no links to external pages).
  type: bool
  choices: 0 (false), 1 (true), 2 (external-only), or 3 (internal-only)
  default: 1

DEFAULT_CHULL_LATEX_LABEL_NAME
  description: Defines the labeling scheme for phases shown on the convex hull. By default, the compound labels are shown, while the prototype label can also be specified. "icsd" shows the ICSD entry number designation (lowest for multiple equivalent ground-state structures reflecting icsd_canonical_auid) if appropriate. Also acceptable: both and none.
  type: string
  choices: "compound", "prototype", "both", "icsd", "none"
  default: ""

DEFAULT_CHULL_LATEX_META_LABELS
  description: Enables verbose labels, including compound, prototype, H_f, T_S, and Delta H_f. Warning, significant overlap of labels should be expected.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_LABELS_OFF_HULL
  description: Enables labels for unstable points.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION
  description: Dictates whether compound labels on the plot are reduced. Only three values are accepted: -1 (default: true unless a binary is requested and an explicit "compounds" option is set for DEFAULT_CHULL_LATEX_LABEL_NAME), 0 (false), and 1 (true).
  type: integer
  choices: -1 (default), 0 (false), 1 (true)
  default: -1

DEFAULT_CHULL_LATEX_HELVETICA_FONT
  description: Toggles the font scheme from Computer Modern (default) to Helvetica.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_FONT_SIZE
  description: Defines the font size of the labels on the convex hull illustration. Warning, other settings may override this default.
  type: string
  choices: "tiny", "scriptsize", "footnotesize", "small", "normalsize", "large" (default), "Large", "LARGE", "huge", "Huge".
  default: ""

DEFAULT_CHULL_LATEX_ROTATE_LABELS
  description: Toggles whether labels are appropriately rotated.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_BOLD_LABELS
  description: Toggles whether labels are bolded. Only three values are accepted: -1 (default: false unless phase is a ternary), 0 (false), 1 (true).
  type: integer
  default: -1
  choices -1 (default), 0 (false), or 1 (true)

####################################################################################################
DEFAULTS CORE

AFLOW_CORE_TEMPERATURE_BEEP
  description: Threshold computer temperature (Celsius) above which the machine beeps.
  type: integer
  default: 56

AFLOW_CORE_TEMPERATURE_HALT
  description: Threshold computer temperature (Celsius) above which AFLOW stops running.
  type: integer
  default: 65

AFLOW_CORE_TEMPERATURE_REFRESH
  description: Refresh rate (seconds) for computer temperature readings.
  type: integer
  default: 5

####################################################################################################
DEFAULTS MACHINE DEPENDENT MPI

MPI_OPTIONS_DUKE_BETA_MPICH
  description: MPI options for --machine=beta.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_BETA_MPICH
  description: MPI command for --machine=beta.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_BETA_MPICH
  description: MPI executable directory for --machine=beta.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_BETA_OPENMPI
  description: MPI options for --machine=beta_openmpi.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_BETA_OPENMPI
  description: MPI command for --machine=beta_openmpi.
  type: string
  default: "/usr/bin/mpirun.openmpi -np"

MPI_BINARY_DIR_DUKE_BETA_OPENMPI
  description: MPI executable directory for --machine=beta_openmpi.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_MATERIALS
  description: MPI options for --machine=materials.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_MATERIALS
  description: MPI command for --machine=materials.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_MATERIALS
  description: MPI executable directory for --machine=materials.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_AFLOWLIB
  description: MPI options for --machine=aflowlib.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_AFLOWLIB
  description: MPI command for --machine=aflowlib.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_AFLOWLIB
  description: MPI executable directory for --machine=aflowlib.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_HABANA
  description: MPI options for --machine=habana.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_HABANA
  description: MPI command for --machine=habana.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_HABANA
  description: MPI executable directory for --machine=habana.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_QRATS_MPICH
  description: MPI options for --machine=qrats.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_QRATS_MPICH
  description: MPI command for --machine=qrats.
  type: string
  default: "/MAIN/bin/MPICH/bin/mpirun -np"

MPI_BINARY_DIR_DUKE_QRATS_MPICH
  description: MPI executable directory for --machine=qrats.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_QFLOW_OPENMPI
  description: MPI options for --machine=qflow.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_QFLOW_OPENMPI
  description: MPI command for --machine=qflow.
  type: string
  default: "/usr/bin/mpirun -n"

MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI
  description: MPI executable directory for --machine=qflow.
  type: string
  default: "/home/bin/"

MPI_OPTIONS_DOD_CONRAD
  description: MPI options for --machine=conrad.
  type: string
  default: ""

MPI_COMMAND_DOD_CONRAD
  description: MPI command for --machine=conrad.
  type: string
  default: "aprun -n"

MPI_BINARY_DIR_DOD_CONRAD
  description: MPI executable directory for --machine=conrad.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_EOS
  description: MPI options for --machine=eos.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_EOS
  description: MPI command for --machine=eos.
  type: string
  default: "/usr/bin/srun -n"

MPI_NCPUS_MPCDF_EOS
  description: MPI core limit for --machine=eos.
  type: string
  default: "32"

MPI_HYPERTHREADING_MPCDF_EOS
  description: Adjusts NCPUS for hyperthreading with --machine=eos. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "NEGLECT"

MPI_BINARY_DIR_MPCDF_EOS
  description: MPI executable directory for --machine=eos.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_DRACO
  description: MPI options for --machine=draco.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_DRACO
  description: MPI command for --machine=draco.
  type: string
  default: "/usr/bin/srun -n"

MPI_NCPUS_MPCDF_DRACO
  description: MPI core limit for --machine=draco.
  type: string
  default: "0"

MPI_HYPERTHREADING_MPCDF_DRACO
  description: Adjusts NCPUS for hyperthreading with --machine=draco. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "OFF"

MPI_BINARY_DIR_MPCDF_DRACO
  description: MPI executable directory for --machine=draco.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_COBRA
  description: MPI options for --machine=cobra.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_COBRA
  description: MPI command for --machine=cobra.
  type: string
  default: "/usr/bin/srun -n"

MPI_NCPUS_MPCDF_COBRA
  description: MPI core limit for --machine=cobra.
  type: string
  default: "0"

MPI_HYPERTHREADING_MPCDF_COBRA
  description: Adjusts NCPUS for hyperthreading with --machine=cobra. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "OFF"

MPI_BINARY_DIR_MPCDF_COBRA
  description: MPI executable directory for --machine=cobra.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_HYDRA
  description: MPI options for --machine=hydra.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_HYDRA
  description: MPI command for --machine=hydra.
  type: string
  default: "poe "

MPI_NCPUS_MPCDF_HYDRA
  description: MPI core limit for --machine=hydra.
  type: string
  default: "0"

MPI_HYPERTHREADING_MPCDF_HYDRA
  description: Adjusts NCPUS for hyperthreading with --machine=hydra. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "OFF"

MPI_BINARY_DIR_MPCDF_HYDRA
  description: MPI executable directory for --machine=hydra.
  type: string
  default: "~/bin/"

MPI_OPTIONS_TERAGRID_RANGER
  description: MPI options for --machine=teragrid_ranger.
  type: string
  default: ""

MPI_COMMAND_TERAGRID_RANGER
  description: MPI command for --machine=teragrid_ranger.
  type: string
  default: "/share/sge6.2/default/pe_scripts/ibrun"

MPI_BINARY_DIR_TERAGRID_RANGER
  description: MPI executable directory for --machine=teragrid_ranger.
  type: string
  default: "~/bin/"

MPI_OPTIONS_TERAGRID_KRAKEN
  description: MPI options for --machine=teragrid_kraken.
  type: string
  default: ""

MPI_COMMAND_TERAGRID_KRAKEN
  description: MPI command for --machine=teragrid_kraken.
  type: string
  default: "aprun -n"

MPI_BINARY_DIR_TERAGRID_KRAKEN
  description: MPI executable directory for --machine=teragrid_kraken.
  type: string
  default: "/nics/a/proj/aflow/bin/"

MPI_OPTIONS_FULTON_MARYLOU
  description: MPI options for --machine=marylou.
  type: string
  default: "export OMP_NUM_THREADS

MPI_COMMAND_FULTON_MARYLOU
  description: MPI command for --machine=marylou.
  type: string
  default: "mpiexec"

MPI_BINARY_DIR_FULTON_MARYLOU
  description: MPI executable directory for --machine=marylou.
  type: string
  default: "/fslgroup/fslg_datamining/bin/"

MPI_OPTIONS_TRINITY_PARSONS
  description: MPI options for --machine=parsons.
  type: string
  default: ""

MPI_COMMAND_TRINITY_PARSONS
  description: MPI command for --machine=parsons.
  type: string
  default: "mpirun -np"

MPI_BINARY_DIR_TRINITY_PARSONS
  description: MPI executable directory for --machine=parsons.
  type: string
  default: "/home/users/aflow/bin/"

MPI_OPTIONS_MACHINE1
  description: MPI options for --machine=host1.
  type: string
  default: ""

MPI_COMMAND_MACHINE1
  description: MPI command for --machine=host1.
  type: string
  default: "...something ..."

MPI_BINARY_DIR_MACHINE1
  description: MPI executable directory for --machine=host1.
  type: string
  default: "/somewhere/"

MPI_OPTIONS_MACHINE2
  description: MPI options for --machine=ohad.
  type: string
  default: ""

MPI_COMMAND_MACHINE2
  description: MPI command for --machine=ohad.
  type: string
  default: "stub not used"

MPI_BINARY_DIR_MACHINE2
  description: MPI executable directory for --machine=ohad.
  type: string
  default: "/home/aflow/bin/"


****************************************************************************************************

For additional information contact: Stefano Curtarolo (stefano@duke.edu) or 
                                    Corey Oses (corey.oses@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
