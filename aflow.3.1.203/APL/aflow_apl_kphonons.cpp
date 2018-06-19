// ***************************************************************************
// *                                                                         *
// *             STEFANO CURTAROLO - Duke University 2003-2018               *
// *                                                                         *
// ***************************************************************************

#include "../aflow.h"
#include "aflow_apl.h"

#define _ASTROPT_APL_OLD_ string("[AFLOW_PHONONS]") //CO 170601, ensure backwards compatibility (we ALWAYS support LEGACY)
#define _ASTROPT_APL_ string("[AFLOW_APL]") //CO 170601
#define _ASTROPT_QHA_ string("[AFLOW_QHA]") //CO 170601
#define _ASTROPT_AAPL_ string("[AFLOW_AAPL]") //CO 170601
//temporary directory for storing QH files
#define _TMPDIR_ string("ARUN.APL.QH.TMP")  //[PINKU]

bool _WITHIN_DUKE_ = false;

//CO fixing cpp version issues with auto_ptr (depreciated)
#if __cplusplus >= 201103L
template <typename T>
using auto_ptr = std::unique_ptr<T>;
#else
using std::auto_ptr;
#endif

namespace KBIN {
void VASP_RunPhonons_APL(_xvasp& xvasp,
                         string AflowIn,
                         _aflags& aflags,
                         _kflags& kflags,
                         _vflags& vflags, ofstream& messageFile) {
  _xinput xinput(xvasp);
  _xflags xflags(vflags);
  return RunPhonons_APL(xinput,AflowIn,aflags,kflags,xflags,messageFile);
}

void RunPhonons_APL(_xinput& xinput,
                    string AflowIn,
                    _aflags& aflags,
                    _kflags& kflags,
                    _xflags& xflags, 
                    ofstream& messageFile) {
  // Test
  if (!(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL)) return;

  //we make certain automatic fixes if we're within our domain, otherwise issue warning/error
  ::_WITHIN_DUKE_ = (aurostd::substring2bool(XHOST.hostname, "nietzsche") || aurostd::substring2bool(XHOST.hostname, "aflowlib") || aurostd::substring2bool(XHOST.hostname, "qrats") || aurostd::substring2bool(XHOST.hostname, "habana") || aurostd::substring2bool(XHOST.hostname, "quser"));

  //determine if we have a consistent mode defined between input and flags
  if(xinput.AFLOW_MODE_VASP){
    if(!xflags.AFLOW_MODE_VASP){
      cerr << "ERROR: KBIN::RunPhonons_APL: mismatch types between xinput(VASP) and xflags(!VASP)" << endl;
      return;
    }
  }else if(xinput.AFLOW_MODE_AIMS){
    if(!xflags.AFLOW_MODE_AIMS){
      cerr << "ERROR: KBIN::RunPhonons_APL: mismatch types between xinput(AIMS) and xflags(!AIMS)" << endl;
      return;
    }
  //}else if(xinput.AFLOW_MODE_ALIEN){  //alien doesn't have xstr, so we ignore
  }else{
    cerr << "ERROR: KBIN::RunPhonons_APL: unknown input type" << endl;
    return;
  }

  //corey
  //fix names regardless of POSCAR style, we need mass!
  if(xinput.AFLOW_MODE_VASP){pflow::fixEmptyAtomNames(xinput.xvasp.str, true);}

  // Setup our logger
  apl::Logger logger(messageFile, aflags);
  //logger.setModuleName("PHONONS");  //will rename later
  
  string _ASTROPT_; //CO 170601
  if (kflags.KBIN_PHONONS_CALCULATION_AAPL) {
    logger.setModuleName("AAPL");  //CO 170601
    _ASTROPT_ = _ASTROPT_AAPL_;    //CO 170601
  } else if (kflags.KBIN_PHONONS_CALCULATION_QHA) {
    logger.setModuleName("QHA");  //CO 170601
    _ASTROPT_ = _ASTROPT_QHA_;    //CO 170601
  } else {
    logger.setModuleName("APL");  //CO 170601
    _ASTROPT_ = _ASTROPT_APL_;    //CO 170601
  }

  logger << "RUNNING..." << apl::endl;

  // CONTROL PARAMETERS FROM _AFLOWIN_ --------------------------------------

  // General switches what to calculate
  aurostd::xoption CALCULATE_PHONON_DISPERSIONS_OPTION; CALCULATE_PHONON_DISPERSIONS_OPTION.option = false;
  aurostd::xoption CALCULATE_PHONON_DOS_OPTION; CALCULATE_PHONON_DOS_OPTION.option = false;
  aurostd::xoption CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION; CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option = false;
  aurostd::xoption CALCULATE_POLAR_CORRECTIONS_OPTION; CALCULATE_POLAR_CORRECTIONS_OPTION.option = false;
  // BEGIN JJPR
  // Anharmonic foces
  aurostd::xoption CALCULATE_APL_OPTION; CALCULATE_APL_OPTION.option = false;
  aurostd::xoption CALCULATE_KAPPA_OPTION; CALCULATE_KAPPA_OPTION.option = false;
  aurostd::xoption USER_BTE_OPTION; USER_BTE_OPTION.xscheme = "RTA"; string USER_BTE = "RTA";
  aurostd::xoption CALCULATE_ISOTOPE_OPTION; CALCULATE_ISOTOPE_OPTION.option = false;
  aurostd::xoption CALCULATE_CUMULATIVEK_OPTION; CALCULATE_CUMULATIVEK_OPTION.option = false;
  aurostd::xoption CALCULATE_BOUNDARY_OPTION; CALCULATE_BOUNDARY_OPTION.option = false;
  aurostd::xoption USER_TDISTORTION_MAGNITUDE_OPTION; USER_TDISTORTION_MAGNITUDE_OPTION.xscheme = "0.015"; double USER_TDISTORTION_MAGNITUDE = 0.015;
  aurostd::xoption USER_THERMALGRID_OPTION; USER_THERMALGRID_OPTION.xscheme = "8x8x8"; string USER_THERMALGRID = "8x8x8";
  aurostd::xoption USER_CUTOFF_DISTANCE_OPTION; USER_CUTOFF_DISTANCE_OPTION.xscheme = "5"; double USER_CUTOFF_DISTANCE = 5; //CO 180409
  aurostd::xoption USER_NANO_SIZE_OPTION; USER_NANO_SIZE_OPTION.xscheme = "13000"; double USER_NANO_SIZE = 13000;
  aurostd::xoption USER_EPS_SUM_OPTION; USER_EPS_SUM_OPTION.xscheme = "1E-5"; double USER_EPS_SUM = 1E-5;
  aurostd::xoption USER_CUTOFF_SHELL_OPTION; USER_CUTOFF_SHELL_OPTION.xscheme = "3"; int USER_CUTOFF_SHELL = 3; //CO 180409
  aurostd::xoption USER_TCT_OPTION; USER_TCT_OPTION.xscheme = "200:700:100"; double USER_TCT_TSTART = 200; double USER_TCT_TEND = 700; double USER_TCT_TSTEP = 100;
  // END JJPR
  
  //PINKU QUASI-HARMONIC START
  aurostd::xoption CALCULATE_GRUNEISEN_OPTION; CALCULATE_GRUNEISEN_OPTION.option = false;
  aurostd::xoption CALCULATE_DISPLACEMENTS_OPTION; CALCULATE_DISPLACEMENTS_OPTION.option = false;
  aurostd::xoption CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option = false;
  aurostd::xoption CALCULATE_EOS_OPTION; CALCULATE_EOS_OPTION.option = false;
  aurostd::xoption CALCULATE_EOS_SUBDIRECTORIES_OPTION; CALCULATE_EOS_SUBDIRECTORIES_OPTION.option = false;
  // different types of fitting options for EOS calculations
  // (1) BM1 => Murnaghan EOS
  // (2) BM2 => Birch-Murnaghan 3rd-order EOS
  // (3) BM3 => Birch-Murnaghan 4th-order EOS
  aurostd::xoption GP_VOL_DISTORTION_OPTION; GP_VOL_DISTORTION_OPTION.xscheme = "0.03"; double GP_VOL_DISTORTION = 0.03;
  aurostd::xoption USER_PROJECTION_DIR_OPTION; USER_PROJECTION_DIR_OPTION.xscheme = "1:1:1"; vector<double> directions(3, 0); directions[0] = 1; directions[1] = 1; directions[2] = 1;   // 3 Miller indices
  aurostd::xoption CUTOFF_FREQ_OPTION; CUTOFF_FREQ_OPTION.xscheme="1e-5"; double CUTOFF_FREQ = 1e-5;  //in amu
  aurostd::xoption EOS_VOLRANGE_DIST_OPTION; EOS_VOLRANGE_DIST_OPTION.xscheme = "-2:4:0.5"; double EOS_VOL_START = -2; double EOS_VOL_END = 4; double EOS_VOL_DISTORTION = 0.5;
  //KPPRA for static calculations, the numbers is AFLOW standard
  //STATIC_KPPRA for static calculations, the numbers are AFLOW standard
  //BANDS_GRID for static calculations
  //NEDOS number of electronic density of state
  aurostd::xoption EOS_KPOINTS_MODE_OPTION; EOS_KPOINTS_MODE_OPTION.xscheme = "8000:10000:20:50000"; string EOS_KPPRA = "8000"; string EOS_STATIC_KPPRA = "10000"; string EOS_BANDS_GRID = "20"; string NEDOS = "50000"; //PINKU make int?
  aurostd::xoption FITTING_TYPE_OPTION; FITTING_TYPE_OPTION.xscheme = "BM1"; string FITTING_TYPE = "BM1";
  //PINKU QUASI-HARMONIC END
  
  // User's control about general phonon engine
  aurostd::xoption USER_ENGINE_OPTION; USER_ENGINE_OPTION.xscheme = "DM"; string USER_ENGINE = "DM";
  aurostd::xoption AUTO_DISTORTIONS_PLUS_MINUS_OPTION; AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option = true;  //CO
  aurostd::xoption USER_DISTORTIONS_PLUS_MINUS_OPTION; USER_DISTORTIONS_PLUS_MINUS_OPTION.option = false;  //CO
  aurostd::xoption USER_DISTORTIONS_XYZ_ONLY_OPTION; USER_DISTORTIONS_XYZ_ONLY_OPTION.option = false;
  aurostd::xoption USER_DISTORTION_MAGNITUDE_OPTION; USER_DISTORTION_MAGNITUDE_OPTION.xscheme = "0.015"; double USER_DISTORTION_MAGNITUDE = 0.015;
  aurostd::xoption USER_ZEROSTATE_OPTION; USER_ZEROSTATE_OPTION.option = false;
  aurostd::xoption USER_HIBERNATE_OPTION; USER_HIBERNATE_OPTION.option = true;
  
  aurostd::xoption USER_FREQFORMAT_OPTION; USER_FREQFORMAT_OPTION.xscheme="THZ | ALLOW_NEGATIVE"; string USER_FREQFORMAT = "THZ | ALLOW_NEGATIVE";
  // User's control of supercell used for calculation
  aurostd::xoption USER_WANTS_RELAX_OPTION; USER_WANTS_RELAX_OPTION.option = false;
  bool USER_WANTS_FULL_SHELL = false;
  //  bool   USER_WANTS_FULL_ATOMS            = false;
  aurostd::xoption USER_MAXSHELL_OPTION; USER_MAXSHELL_OPTION.xscheme = "-1"; int USER_MAXSHELL = -1;  //CO set by OPTION
  aurostd::xoption USER_MINSHELL_OPTION; USER_MINSHELL_OPTION.xscheme = "6"; int USER_MINSHELL = 6;  //CO set by OPTION
  aurostd::xoption USER_MINATOMS_OPTION; USER_MINATOMS_OPTION.xscheme = "100"; int USER_MINATOMS = 100;  //CO set by OPTION
  bool USER_MINATOMS_RESTICTED_OPTION=0;  //CO 180404 - adding restricted variant to minatoms
  aurostd::xoption USER_SUPERCELL_OPTION; USER_SUPERCELL_OPTION.xscheme = ""; string USER_SUPERCELL = "";

  // User's control about calculation of dispersion curves
  // there are several ways a user can change the default path
  // a USERPATH approach is either a subset or a swapping of the default path
  // OWNPATH is a true user defined path, either by providing frac or cartesian coordinates
  // OWNPATH can be used in combination with LABELS
  // useful for matching with experimental results
  aurostd::xoption USER_DC_INITLATTICE_OPTION; USER_DC_INITLATTICE_OPTION.xscheme = ""; string USER_DC_INITLATTICE = "";  //CO set by OPTION
  aurostd::xoption USER_DC_INITSG_OPTION; USER_DC_INITSG_OPTION.xscheme = ""; int USER_DC_DCINITSG = 1;
  aurostd::xoption USER_DC_INITCOORDS_FRAC_OPTION; USER_DC_INITCOORDS_FRAC_OPTION.xscheme = ""; string USER_DC_INITCOORDS_FRAC = ""; //CO set by OPTION
  aurostd::xoption USER_DC_INITCOORDS_CART_OPTION; USER_DC_INITCOORDS_CART_OPTION.xscheme = ""; string USER_DC_INITCOORDS_CART = ""; //CO set by OPTION
  aurostd::xoption USER_DC_INITCOORDS_LABELS_OPTION; USER_DC_INITCOORDS_LABELS_OPTION.xscheme = ""; string USER_DC_INITCOORDS_LABELS = ""; //CO set by OPTION
  aurostd::xoption USER_DC_USERPATH_OPTION; USER_DC_USERPATH_OPTION.xscheme = ""; string USER_DC_USERPATH = ""; //CO set by OPTION
  aurostd::xoption USER_DC_NPOINTS_OPTION; USER_DC_NPOINTS_OPTION.xscheme = "100"; int USER_DC_NPOINTS = 100;

  // User's control about phonon DOS
  aurostd::xoption USER_DOS_MESH_OPTION; USER_DOS_MESH_OPTION.xscheme = "21x21x21"; string USER_DOS_MESH = "21x21x21";
  aurostd::xoption USER_DOS_NPOINTS_OPTION; USER_DOS_NPOINTS_OPTION.xscheme = "2000"; int USER_DOS_NPOINTS = 2000;
  aurostd::xoption USER_DOS_METHOD_OPTION; USER_DOS_METHOD_OPTION.xscheme = "LT"; string USER_DOS_METHOD = "LT";
  aurostd::xoption USER_DOS_SMEAR_OPTION; USER_DOS_SMEAR_OPTION.xscheme = "0.0"; double USER_DOS_SMEAR = 0.0;

  // User's control of thermodynamic properties calculation
  aurostd::xoption USER_TPT_OPTION; USER_TPT_OPTION.xscheme="0:2000:10"; double USER_TP_TSTART = 0; double USER_TP_TEND = 2000;  double USER_TP_TSTEP = 10;

  //if current dir DCUSERPATH is same as other PHONON sub-directories, if not make changes accordingly
  //PINKU QUASI-HARMONIC START
  apl::check_consistency_aflow_apl check_consistency(logger);
  check_consistency.getdir_name(aflags.Directory);
  check_consistency.check_consistency_in_aflow(AflowIn);
  //PINKU QUASI-HARMONIC END

  // TAR
  //  bool DO_TAR                             = false;

  // Get user's parameters from _AFLOWIN_ ///////////////////////////////////

  //for when we need double or string variable extraction
  vector<string> tokens;
  string test;
  bool override_option;
  bool found_supercell;
  try {
    // BEGIN JJPR
    // Anharmonic and thermal conductivity options

    if(kflags.KBIN_PHONONS_CALCULATION_AAPL){
    // KAPPA, e.g. KAPPA = y
      CALCULATE_KAPPA_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_AAPL; //recycle what we parsed earlier
      logger << _ASTROPT_ << "CALC is" << ( CALCULATE_KAPPA_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;  //CO 170601
      logger << "Anharmonic force constants will" << ( CALCULATE_KAPPA_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-AAPL." << apl::endl;   //CO 170601

      // Cutoff Radious (in Angs. and real space)
      USER_CUTOFF_DISTANCE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "CUT_RAD=" + "|" + _ASTROPT_APL_OLD_ + "CUT_RAD="), USER_CUTOFF_DISTANCE_OPTION.option, USER_CUTOFF_DISTANCE_OPTION.xscheme); //CO 170601
      USER_CUTOFF_DISTANCE = USER_CUTOFF_DISTANCE_OPTION.content_double; //CO 170601
      logger << (USER_CUTOFF_DISTANCE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUT_RAD=" << USER_CUTOFF_DISTANCE << "." << apl::endl;
      logger << "The cutoff to compute the anharmonic IFCs will be " << abs(USER_CUTOFF_DISTANCE) << " Angs." << apl::endl;

      // Get the users magnitute of the thershold for the sumrule
      USER_EPS_SUM_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SUMRULE=" + "|" + _ASTROPT_APL_OLD_ + "SUMRULE="), USER_EPS_SUM_OPTION.option, USER_EPS_SUM_OPTION.xscheme); //CO 170601
      USER_EPS_SUM = USER_EPS_SUM_OPTION.content_double; //CO 170601
      logger << (USER_EPS_SUM_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SUMRULE=" << USER_EPS_SUM << "." << apl::endl;
      logger << "Convergence criteria for the sumrules is set to " << abs(USER_EPS_SUM) << " eV/Angs.^3." << apl::endl;

      // Cutoff radius (Neighbors)
      USER_CUTOFF_SHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "CUT_SHELL=" + "|" + _ASTROPT_APL_OLD_ + "CUT_SHELL="), USER_CUTOFF_SHELL_OPTION.option, USER_CUTOFF_SHELL_OPTION.xscheme); //CO 170601
      USER_CUTOFF_SHELL = USER_CUTOFF_SHELL_OPTION.content_double; //CO 170601
      logger << (USER_CUTOFF_SHELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUT_SHELL=" << USER_CUTOFF_SHELL << "." << apl::endl;
      logger << "The calculation of the anharmonic IFCs will consider up to the " << USER_CUTOFF_SHELL << "th coordination shell." << apl::endl;

      // ME 180501 - If the user only specifies CUT_SHELL or CUT_RAD, unset the default values
      if (USER_CUTOFF_SHELL_OPTION.isentry && ! USER_CUTOFF_DISTANCE_OPTION.isentry){USER_CUTOFF_DISTANCE = 0.0;}
      if (! USER_CUTOFF_SHELL_OPTION.isentry && USER_CUTOFF_DISTANCE_OPTION.isentry){USER_CUTOFF_SHELL = 0;}

      // Get the users magnitude of the distortion vector (in Angs. and real space)
      USER_TDISTORTION_MAGNITUDE_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "TDMAG=" + "|" +_ASTROPT_APL_OLD_ + "TDMAG=" + "|" + _ASTROPT_ + "TDISMAG=" + "|" +_ASTROPT_APL_OLD_ + "TDISMAG="), USER_TDISTORTION_MAGNITUDE_OPTION.option, USER_TDISTORTION_MAGNITUDE_OPTION.xscheme);  //CO 170621 - TDISMAG legacy
      USER_TDISTORTION_MAGNITUDE = USER_TDISTORTION_MAGNITUDE_OPTION.content_double;
      logger << (USER_TDISTORTION_MAGNITUDE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "TDMAG=" << USER_TDISTORTION_MAGNITUDE << "." << apl::endl;
      logger << "The distortion magnitude for anharmonic IFCs will be " << USER_TDISTORTION_MAGNITUDE << " Angs." << apl::endl;

      // ISOTOPE, e.g. ISOTOPE = y
      CALCULATE_ISOTOPE_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "ISOTOPE=" + "|" +_ASTROPT_APL_OLD_ + "ISOTOPE="), CALCULATE_ISOTOPE_OPTION.option, CALCULATE_ISOTOPE_OPTION.xscheme); //CO 170601
      logger << (CALCULATE_ISOTOPE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "ISOTOPE=" << (CALCULATE_ISOTOPE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      logger << "Isotope effects " << (CALCULATE_ISOTOPE_OPTION.option ? "will" : "will NOT") << " be considered in the calculation." << apl::endl;

      //scattering at boundaries is handled either with BOUNDARY_OPTION + NANO_SIZE or CUMULATIVEK, if handled AT ALL
      // BOUNDARY, e.g. BOUNDARY = y
      CALCULATE_BOUNDARY_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "BOUNDARY=" + "|" +_ASTROPT_APL_OLD_ + "BOUNDARY="), CALCULATE_BOUNDARY_OPTION.option, CALCULATE_BOUNDARY_OPTION.xscheme); //CO 170601
      logger << (CALCULATE_BOUNDARY_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "BOUNDARY=" << (CALCULATE_BOUNDARY_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      logger << "The boundary term will " << (CALCULATE_BOUNDARY_OPTION.option ? "": "NOT ") << "be considered in the scattering time." << apl::endl;
      // CUMULATIVEK, e.g. CUMULATIVEK = y
      CALCULATE_CUMULATIVEK_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "CUMULATIVEK=" + "|" +_ASTROPT_APL_OLD_ + "CUMULATIVEK="), CALCULATE_CUMULATIVEK_OPTION.option, CALCULATE_CUMULATIVEK_OPTION.xscheme); //CO 170601
      if(!(CALCULATE_BOUNDARY_OPTION.option || CALCULATE_CUMULATIVEK_OPTION.option) ){
        logger << (CALCULATE_CUMULATIVEK_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUMULATIVEK=" << (CALCULATE_CUMULATIVEK_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
        logger << "Boundary effects will NOT be considered in calculating the lattice thermal conductivity." << apl::endl;
      } else {
        // BOUNDARY takes precendence
        override_option = false;
        if(CALCULATE_BOUNDARY_OPTION.option && CALCULATE_CUMULATIVEK_OPTION.option){CALCULATE_CUMULATIVEK_OPTION.option=false; override_option = true;}
        logger << (CALCULATE_CUMULATIVEK_OPTION.isentry || override_option ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUMULATIVEK=" << (CALCULATE_CUMULATIVEK_OPTION.option ? "ON" : "OFF") << (override_option ? " (overridden by BOUNDARY)": "") << "." << apl::endl;
        logger << "The cumulative lattice thermal conductivity will " << (CALCULATE_CUMULATIVEK_OPTION.option ? "" : "NOT") << " be calculated." << apl::endl;
      }
      // Cutoff Radious (NanoSize)
      if(CALCULATE_BOUNDARY_OPTION.option){
        USER_NANO_SIZE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "NANO_SIZE=" + "|" + _ASTROPT_APL_OLD_ + "NANO_SIZE="), USER_NANO_SIZE_OPTION.option, USER_NANO_SIZE_OPTION.xscheme); //CO 170601
        USER_NANO_SIZE = USER_NANO_SIZE_OPTION.content_double; //CO 170601
        logger << (USER_NANO_SIZE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "NANO_SIZE=" << USER_NANO_SIZE << "." << apl::endl;
        logger << "The boundary scattering will be computed for nanoparticles with " << abs(USER_NANO_SIZE) << " nm size." << apl::endl;
      }

      // BTE, e.g. BTE = RTA
      USER_BTE_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "BTE=" + "|" +_ASTROPT_APL_OLD_ + "BTE="), USER_BTE_OPTION.option, USER_BTE_OPTION.xscheme); //CO 170601
      USER_BTE = USER_BTE_OPTION.content_string;
      transform(USER_BTE.begin(), USER_BTE.end(), USER_BTE.begin(), toupper);
      logger << (USER_BTE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "BTE=" << USER_BTE << "." << apl::endl;
      if (USER_BTE == string("RTA")){logger << "Boltzmann Transport Equation will be solved using the Relaxation Time Approximation Method (RTA)." << apl::endl;}
      else if (USER_BTE == string("FULL")) {logger << "Boltzmann Transport Equation will be solved using an iterative scheme (FULL)." << apl::endl;}
      else {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"BTE. Specify as BTE=RTA or FULL.");}

      // THERMALGRID, e.g., THERMALGRID = 2x2x2
      USER_THERMALGRID_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "THERMALGRID=" + "|" + _ASTROPT_APL_OLD_ + "THERMALGRID="), USER_THERMALGRID_OPTION.option, USER_THERMALGRID_OPTION.xscheme); //CO 170601
      USER_THERMALGRID = USER_THERMALGRID_OPTION.content_string; //CO 170601
      tokens.clear();
      apl::tokenize(USER_THERMALGRID, tokens, string(" xX"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"THERMALGRID. Specify as THERMALGRID=10x10x10.");}
      logger << (USER_THERMALGRID_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "THERMALGRID=" << USER_THERMALGRID << "." << apl::endl;

      // TCT, e.g., TCT = 1000:2000:10 -> temperature from 1000 to 2000 K by step 10 K
      // no default here
      USER_TCT_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "TCT=" + "|" + _ASTROPT_APL_OLD_ + "TCT="), USER_TCT_OPTION.option, USER_TCT_OPTION.xscheme); //CO 170601
      tokens.clear();
      apl::tokenize(USER_TCT_OPTION.content_string, tokens, string(" :"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"TCT. Specify as TCT=1000:2000:10.");}
      USER_TCT_TSTART = aurostd::string2utype<double>(tokens.at(0));
      USER_TCT_TEND = aurostd::string2utype<double>(tokens.at(1));
      USER_TCT_TSTEP = aurostd::string2utype<double>(tokens.at(2));
      override_option = false;
      if (USER_TCT_TSTART == 0) {
        logger << apl::warning << "Thermal conductivity at 0K it is an easy task... infinite!! Let's try something more difficult like 50K." << apl::endl;
        USER_TCT_TSTART = 50;
        override_option = true;
      }
      logger << (USER_TCT_OPTION.isentry || override_option ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "THERMALGRID=" << USER_TCT_OPTION.content_string << "." << apl::endl;
      logger << "The thermodynamic properties will be calculated in temperature range <" << USER_TCT_TSTART << "," << USER_TCT_TEND << "> with step size " << USER_TCT_TSTEP << " K." << apl::endl;
    }

    //Anharmonic forces and thermal conductivity
    //END JJPR ANHARMONIC

    //PINKU QUASI-HARMONIC START
    if(!CALCULATE_KAPPA_OPTION.option){
      if(kflags.KBIN_PHONONS_CALCULATION_QHA){
        CALCULATE_GRUNEISEN_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA; //recycle what we parsed earlier
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen parameter will" << ( CALCULATE_GRUNEISEN_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;

        CUTOFF_FREQ_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "CUTOFF_FREQ=" + "|" + _ASTROPT_APL_OLD_ + "CUTOFF_FREQ="), CUTOFF_FREQ_OPTION.option, CUTOFF_FREQ_OPTION.xscheme); //CO 170601
        CUTOFF_FREQ = CUTOFF_FREQ_OPTION.content_double;
        //test of stupidity
        if (CUTOFF_FREQ_OPTION.isentry) { //CO 170601
          tokens.clear();
          apl::tokenize(CUTOFF_FREQ_OPTION.content_string, tokens, string(" "));
          if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"CUTOFF_FREQ. Specify as CUTOFF_FREQ=0.01.");}
        }
        logger << (CUTOFF_FREQ_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUTOFF_FREQ=" << CUTOFF_FREQ_OPTION.content_string << "." << apl::endl;

        CALCULATE_DISPLACEMENTS_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "DISPLACEMENTS=" + "|" + _ASTROPT_APL_OLD_ + "DISPLACEMENTS="), CALCULATE_DISPLACEMENTS_OPTION.option, CALCULATE_DISPLACEMENTS_OPTION.xscheme); //CO 170601
        logger << (CALCULATE_DISPLACEMENTS_OPTION.isentry ? "Setting" : "DEFAULT" ) << " " << _ASTROPT_ << "DISPLACEMENTS=" << (CALCULATE_DISPLACEMENTS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

        USER_PROJECTION_DIR_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "PROJECTION_DIR=" + "|" + _ASTROPT_APL_OLD_ + "PROJECTION_DIR="), USER_PROJECTION_DIR_OPTION.option, USER_PROJECTION_DIR_OPTION.xscheme); //CO 170601
        tokens.clear();
        apl::tokenize(USER_PROJECTION_DIR_OPTION.content_string, tokens, string(" :"));
        if (tokens.size() != 3){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"PROJECTION_DIR. Specify as PROJECTION_DIR=1:1:1.");}
        directions.clear(); directions.resize(3, 0.);
        directions[0] = aurostd::string2utype<double>(tokens.at(0)); directions[1] = aurostd::string2utype<double>(tokens.at(1)); directions[2] = aurostd::string2utype<double>(tokens.at(2));
        logger << (USER_PROJECTION_DIR_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "PROJECTION_DIR=" << USER_PROJECTION_DIR_OPTION.content_string << "." << apl::endl;
        logger << "Displacements will be calculated along direction = <" << directions[0] << "," << directions[1] << "," << directions[2] << ">." << apl::endl;

        //PINKU we may want to rename this based on new CALC tag
        CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISENSD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISENSD="), CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.xscheme); //CO 170601
        logger << (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISENSD=" << (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

        GP_VOL_DISTORTION_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "GP_VOL_DISTORTION_PERCENTAGE=" + "|" + _ASTROPT_APL_OLD_ + "GP_VOL_DISTORTION_PERCENTAGE="), GP_VOL_DISTORTION_OPTION.option, GP_VOL_DISTORTION_OPTION.xscheme); //CO 170601
        GP_VOL_DISTORTION=GP_VOL_DISTORTION_OPTION.content_double;
        //test of stupidity
        if (GP_VOL_DISTORTION_OPTION.isentry) {
          tokens.clear();
          apl::tokenize(GP_VOL_DISTORTION_OPTION.content_string, tokens, string(" "));
          if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"GP_VOL_DISTORTION_PERCENTAGE. Specify as GP_VOL_DISTORTION_PERCENTAGE=0.03.");}
        }
        logger << (GP_VOL_DISTORTION_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GP_VOL_DISTORTION_PERCENTAGE=" << GP_VOL_DISTORTION << "." << apl::endl;
        if (!CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option) {logger << "Two directories are going to be created with the following distortion magnitude = " << GP_VOL_DISTORTION << "." << apl::endl;}

        if (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option) {
          CALCULATE_EOS_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "EOSSD=" + "|" + _ASTROPT_APL_OLD_ + "EOSSD="), CALCULATE_EOS_SUBDIRECTORIES_OPTION.option, CALCULATE_EOS_SUBDIRECTORIES_OPTION.xscheme); //CO 170601
          logger << (CALCULATE_EOS_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "EOSSD=" << (CALCULATE_EOS_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
        }

        CALCULATE_EOS_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "EOS=" + "|" + _ASTROPT_APL_OLD_ + "EOS="), CALCULATE_EOS_OPTION.option, CALCULATE_EOS_OPTION.xscheme); //CO 170601
        logger << (CALCULATE_EOS_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "EOSSD=" << (CALCULATE_EOS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

        if (CALCULATE_EOS_OPTION.option) {
          EOS_VOLRANGE_DIST_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_VOLRANGE_DIST=" + "|" + _ASTROPT_APL_OLD_ + "EOS_VOLRANGE_DIST="), EOS_VOLRANGE_DIST_OPTION.option, EOS_VOLRANGE_DIST_OPTION.xscheme); //CO 170601
          tokens.clear();
          apl::tokenize(EOS_VOLRANGE_DIST_OPTION.content_string, tokens, string(" :"));
          if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_VOLRANGE_DIST. Specify as EOS_VOLRANGE_DIST=-6:6:1.");}
          EOS_VOL_START = aurostd::string2utype<double>(tokens.at(0));
          EOS_VOL_END = aurostd::string2utype<double>(tokens.at(1));
          EOS_VOL_DISTORTION = aurostd::string2utype<double>(tokens.at(2));
          logger << (EOS_VOLRANGE_DIST_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_VOLRANGE_DIST=" << EOS_VOLRANGE_DIST_OPTION.content_string << "." << apl::endl;
          logger << "The EOS properties will be calculated in distortion range <" << EOS_VOL_START << "," << EOS_VOL_END << "," << EOS_VOL_DISTORTION << "." << apl::endl;

          EOS_KPOINTS_MODE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_KPOINTS_MODE=" + "|" + _ASTROPT_APL_OLD_ + "EOS_KPOINTS_MODE="), EOS_KPOINTS_MODE_OPTION.option, EOS_KPOINTS_MODE_OPTION.xscheme); //CO 170601
          tokens.clear();
          apl::tokenize(EOS_KPOINTS_MODE_OPTION.content_string, tokens, string(" :"));
          if (tokens.size() != 4) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_KPOINTS_MODE. Specify as EOS_KPOINTS_MODE=5000:10000:20:1000.");}
          EOS_KPPRA = tokens.at(0);
          EOS_STATIC_KPPRA = tokens.at(1);
          EOS_BANDS_GRID = tokens.at(2);
          NEDOS = tokens.at(3);
          logger << (EOS_KPOINTS_MODE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_KPOINTS_MODE=" << EOS_KPOINTS_MODE_OPTION.content_string << "." << apl::endl;
          logger << "EOS_KPPRA, EOS_STATIC_KPPRA, EOS_BANDS_GRID, and NEDOS are set to <" << EOS_KPPRA << "," << EOS_STATIC_KPPRA << "," << EOS_BANDS_GRID << "," << NEDOS << ">, respectively ." << apl::endl;

          FITTING_TYPE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "FITTING_TYPE=" + "|" + _ASTROPT_APL_OLD_ + "FITTING_TYPE="), FITTING_TYPE_OPTION.option, FITTING_TYPE_OPTION.xscheme); //CO 170601
          FITTING_TYPE = FITTING_TYPE_OPTION.content_string;
          //test of stupidity
          if (FITTING_TYPE_OPTION.isentry) { //CO 170601
            tokens.clear();
            apl::tokenize(FITTING_TYPE_OPTION.content_string, tokens, string(" "));
            if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"FITTING_TYPE. Specify as FITTING_TYPE=MB2.");}
          }
          logger << (FITTING_TYPE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "FITTING_TYPE=" << FITTING_TYPE_OPTION.content_string << "." << apl::endl;
          logger << "EOS fitting type found = " << FITTING_TYPE << "." << apl::endl;
        }
      }
    }
    // PINKU QUASI-HARMONIC END

    if( !(CALCULATE_KAPPA_OPTION.option||CALCULATE_GRUNEISEN_OPTION.option) ){
      CALCULATE_APL_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_APL;
      logger << _ASTROPT_ << "CALC is" << ( CALCULATE_APL_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;  //CO 170601
      logger << "Harmonic force constants will" << ( CALCULATE_APL_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-APL." << apl::endl;   //CO 170601
    }  //CO 170601

    //BELOW HERE, all tags all shared, so we need to look for all combinations

    // HIBERNATE, e.g. HIBERNATE = y
    USER_HIBERNATE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "HIBERNATE=" + "|" + _ASTROPT_QHA_ + "HIBERNATE=" + "|" + _ASTROPT_AAPL_ + "HIBERNATE=" + "|" + _ASTROPT_APL_OLD_ + "HIBERNATE="), USER_HIBERNATE_OPTION.option, USER_HIBERNATE_OPTION.xscheme); //CO 170601
    logger << (USER_HIBERNATE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "HIBERNATE=" << (USER_HIBERNATE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The hibernate feature is switched " << (USER_HIBERNATE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

    // FREQFORMAT, e.g. FREQFORMAT = "THz | allow_negative"
    //COREY, would help if you had some sort of stupidity test here
    USER_FREQFORMAT_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "FREQFORMAT=" + "|" + _ASTROPT_QHA_ + "FREQFORMAT=" + "|" + _ASTROPT_AAPL_ + "FREQFORMAT=" + "|" + _ASTROPT_APL_OLD_ + "FREQFORMAT="), USER_FREQFORMAT_OPTION.option, USER_FREQFORMAT_OPTION.xscheme); //CO 170601
    USER_FREQFORMAT = USER_FREQFORMAT_OPTION.content_string;  //CO 170601
    transform(USER_FREQFORMAT.begin(), USER_FREQFORMAT.end(), USER_FREQFORMAT.begin(), toupper);
    logger << (USER_FREQFORMAT_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "FREQFORMAT=" << USER_FREQFORMAT << "." << apl::endl;
    logger << "The frequency will be returned in this format: " << USER_FREQFORMAT << "." << apl::endl;

    // ENGINE, e.g., ENGINE = DM or LR
    USER_ENGINE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "ENGINE=" + "|" + _ASTROPT_QHA_ + "ENGINE=" + "|" + _ASTROPT_AAPL_ + "ENGINE=" + "|" + _ASTROPT_APL_OLD_ + "ENGINE="), USER_ENGINE_OPTION.option, USER_ENGINE_OPTION.xscheme); //CO 170601
    USER_ENGINE = USER_ENGINE_OPTION.content_string;
    transform(USER_ENGINE.begin(), USER_ENGINE.end(), USER_ENGINE.begin(), toupper);
    logger << (USER_ENGINE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "ENGINE=" << USER_ENGINE << "." << apl::endl;
    if (USER_ENGINE == string("DM")) {logger << "The phonon calculator engine: Direct Method (DM)." << apl::endl;}
    else if (USER_ENGINE == string("GSA")) {
      //CO generally redirects to DM, the distinction between DM and GSA is obsolete
      logger << "The Generalized Supercell Approach (GSA) phonon calculator now directs to another engine: Direct Method (DM)." << apl::endl;
      USER_ENGINE = "DM";
    } else if (USER_ENGINE == string("LR")) {logger << "The phonon calculator engine: Linear Response (LR)." << apl::endl;}
    else {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"ENGINE. Specify as ENGINE=DM or LR.");}

    // Get the users magnitute of the distortion vector (in Angs. and real space)
    if (USER_ENGINE == string("DM")) {
      USER_DISTORTION_MAGNITUDE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DMAG=" + "|" + _ASTROPT_QHA_ + "DMAG=" + "|" + _ASTROPT_AAPL_ + "DMAG=" + "|" + _ASTROPT_APL_OLD_ + "DMAG=" + "|" + _ASTROPT_APL_ + "DISMAG=" + "|" + _ASTROPT_QHA_ + "DISMAG=" + "|" + _ASTROPT_AAPL_ + "DISMAG=" + "|" + _ASTROPT_APL_OLD_ + "DISMAG="), USER_DISTORTION_MAGNITUDE_OPTION.option, USER_DISTORTION_MAGNITUDE_OPTION.xscheme); //CO 170601, 170621 DISMAG legacy
      USER_DISTORTION_MAGNITUDE = USER_DISTORTION_MAGNITUDE_OPTION.content_double;
      logger << (USER_DISTORTION_MAGNITUDE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DMAG=" << USER_DISTORTION_MAGNITUDE << "." << apl::endl;
      logger << "The distortion magnitude will be " << USER_DISTORTION_MAGNITUDE << " Angs." << apl::endl;

      // Get flag for generation of displacements with positive/negative magnitude
      USER_DISTORTIONS_PLUS_MINUS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DPM=" + "|" + _ASTROPT_QHA_ + "DPM=" + "|" + _ASTROPT_AAPL_ + "DPM=" + "|" + _ASTROPT_APL_OLD_ + "DPM="), USER_DISTORTIONS_PLUS_MINUS_OPTION.option, USER_DISTORTIONS_PLUS_MINUS_OPTION.xscheme); //CO 170601
      AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option = !USER_DISTORTIONS_PLUS_MINUS_OPTION.option;
      if (USER_DISTORTIONS_PLUS_MINUS_OPTION.isentry){logger << (USER_DISTORTIONS_PLUS_MINUS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DPM=" << (USER_DISTORTIONS_PLUS_MINUS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;}
      if (USER_DISTORTIONS_PLUS_MINUS_OPTION.option){logger << "Each distortion will be generated with both positive and negative magnitudes." << apl::endl;}
      if (USER_DISTORTIONS_PLUS_MINUS_OPTION.isentry && !USER_DISTORTIONS_PLUS_MINUS_OPTION.option){logger << apl::warning << "Distortions will only be considered in one direction (no \"negative\" distortions) - this is NOT recommended." << apl::endl;}
      if (AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option){logger << "DEFAULT " << _ASTROPT_ << "DPM=AUTO which considers negative distortions on a per-site basis." << apl::endl;} //corey AUTO

      // Get flag for generation of displacements only along the x, y, and z axis
      USER_DISTORTIONS_XYZ_ONLY_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DXYZONLY=" + "|" + _ASTROPT_QHA_ + "DXYZONLY=" + "|" + _ASTROPT_AAPL_ + "DXYZONLY=" + "|" + _ASTROPT_APL_OLD_ + "DXYZONLY="), USER_DISTORTIONS_XYZ_ONLY_OPTION.option, USER_DISTORTIONS_XYZ_ONLY_OPTION.xscheme); //CO 170601
      logger << (USER_DISTORTIONS_XYZ_ONLY_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DXYZONLY=" << (USER_DISTORTIONS_XYZ_ONLY_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      if (USER_DISTORTIONS_XYZ_ONLY_OPTION.option){logger << "Only distortions along the x, y, and z axis will be used." << apl::endl;}
      else{logger << "Distortions will be created along lattice vectors including faces and body diagonals." << apl::endl;}

      // ZEROSTATE; One next calculation will be done with no distortion, a such obtained
      // forces will be subtracted from the all forces obtained with distortions
      USER_ZEROSTATE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "ZEROSTATE=" + "|" + _ASTROPT_QHA_ + "ZEROSTATE=" + "|" + _ASTROPT_AAPL_ + "ZEROSTATE=" + "|" + _ASTROPT_APL_OLD_ + "ZEROSTATE="), USER_ZEROSTATE_OPTION.option, USER_ZEROSTATE_OPTION.xscheme); //CO 170601
      logger << (USER_ZEROSTATE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "ZEROSTATE=" << (USER_ZEROSTATE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      if (USER_ZEROSTATE_OPTION.option){logger << "The zero state forces will be also calculated and subtracted from all forces." << apl::endl;}
    }

    // Do polar correction
    if (USER_ENGINE == string("LR") || USER_ENGINE == string("DM")) {
      CALCULATE_POLAR_CORRECTIONS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "POLAR=" + "|" + _ASTROPT_QHA_ + "POLAR=" + "|" + _ASTROPT_AAPL_ + "POLAR=" + "|" + _ASTROPT_APL_OLD_ + "POLAR="), CALCULATE_POLAR_CORRECTIONS_OPTION.option, CALCULATE_POLAR_CORRECTIONS_OPTION.xscheme); //CO 170601
      logger << (CALCULATE_POLAR_CORRECTIONS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "POLAR=" << (CALCULATE_POLAR_CORRECTIONS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      logger << "The calculation of POLAR MATERIALS corrections is switched " << (CALCULATE_POLAR_CORRECTIONS_OPTION.option? "ON": "OFF") << "." << apl::endl;
    }

    //fix vasp bin for LR or DM+POLAR
    if (USER_ENGINE == string("LR") || (USER_ENGINE == string("DM") && CALCULATE_POLAR_CORRECTIONS_OPTION.option)) {
      if(xflags.AFLOW_MODE_VASP){
      try {
        // Check the version of VASP binary
        logger << "Checking VASP version ... ";
        string vaspVersion;
        vaspVersion = apl::getVASPVersionString( (kflags.KBIN_MPI ? kflags.KBIN_MPI_BIN : kflags.KBIN_BIN ) );
        if (!vaspVersion.empty()) {
          logger << "[" << vaspVersion << "]";
          if ((vaspVersion[0] - '0') < 5) { //cool way of getting ascii value:  https://stackoverflow.com/questions/36310181/char-subtraction-in-c
            logger << apl::warning << "." << apl::endl;
            if(_WITHIN_DUKE_){
              kflags.KBIN_BIN = DEFAULT_VASP5_BIN;
              kflags.KBIN_MPI_BIN = DEFAULT_VASP5_MPI_BIN;
              logger << apl::warning << "Modifying VASP bin to " << kflags.KBIN_BIN << " (Duke machine AUTO modification)." << apl::endl;
            } else {
              throw apl::APLRuntimeError("The LR engine needs VASP5 or higher version.");
            }
          }else{logger << " OK." << apl::endl;}
        } else {
          logger << "Failed." << apl::warning << apl::endl; 
          throw apl::APLLogicError("Unexpected binary format.");
        }
      } catch (apl::APLLogicError& e) {
        logger << apl::warning << "Failed to identify the version of VASP binary." << apl::endl;
        logger << apl::warning << e.what() << apl::endl;
      }
    }
    }

    // SUPERCELL ---------------------------------------------------------

    // RELAX, Should be the primitive structure relax before any supercell is build?
    //CO looks like this has yet to be implemented, need to figure this out
    USER_WANTS_RELAX_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "RELAX=" + "|" + _ASTROPT_QHA_ + "RELAX=" + "|" + _ASTROPT_AAPL_ + "RELAX=" + "|" + _ASTROPT_APL_OLD_ + "RELAX="), USER_WANTS_RELAX_OPTION.option, USER_WANTS_RELAX_OPTION.xscheme); //CO 170601
    if(USER_WANTS_RELAX_OPTION.isentry){
      logger << "Setting " << _ASTROPT_ << "RELAX=" << (USER_WANTS_RELAX_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      if (USER_WANTS_RELAX_OPTION.option){logger << "The primitive cell is going to relax before the supercell build." << apl::endl;}
    }

    found_supercell = false;
    // SUPERCELL, e.g., SUPERCELL = 2x2x2
    USER_SUPERCELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "SUPERCELL=" + "|" + _ASTROPT_QHA_ + "SUPERCELL=" + "|" + _ASTROPT_AAPL_ + "SUPERCELL=" + "|" + _ASTROPT_APL_OLD_ + "SUPERCELL="), USER_SUPERCELL_OPTION.option, USER_SUPERCELL_OPTION.xscheme); //CO 170601
    USER_SUPERCELL = USER_SUPERCELL_OPTION.content_string;
    if(USER_SUPERCELL_OPTION.isentry){
      tokens.clear();
      apl::tokenize(USER_SUPERCELL, tokens, string(" xX"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SUPERCELL. Specify as SUPERCELL=2x2x2.");}
      // If the MAXSHELL was not specify, clear the default setting for shell restriction
      //if (!aurostd::substring2bool(AflowIn, _ASTROPT_ + "MINSHELL=", TRUE)){USER_MINSHELL = -1;}
      //CO, I think this should be MAX, not MIN
      //if(!USER_MAXSHELL_OPTION.isentry){USER_MAXSHELL = -1;}  //not really important
      logger << (USER_SUPERCELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SUPERCELL=" << USER_SUPERCELL << "." << apl::endl;
      logger << "Supercell will be built with dimensions " << USER_SUPERCELL << "." << apl::endl;
      found_supercell = true;
    }

    // Get the users minimum atoms which will be included into calculation
    // takes back seat to SUPERCELL
    if(!found_supercell){
      USER_MINATOMS_OPTION.options2entry(AflowIn, string( 
            _ASTROPT_APL_ + "MINATOMS=" + "|" + _ASTROPT_QHA_ + "MINATOMS=" + "|" + _ASTROPT_AAPL_ + "MINATOMS=" + "|" + _ASTROPT_APL_OLD_ + "MINATOMS=" + "|" +
            _ASTROPT_APL_ + "MINATOMS_RESTRICTED=" + "|" + _ASTROPT_QHA_ + "MINATOMS_RESTRICTED=" + "|" + _ASTROPT_AAPL_ + "MINATOMS_RESTRICTED=" + "|" + _ASTROPT_APL_OLD_ + "MINATOMS_RESTRICTED="  //CO 180418 - restricted means all dims are equal
            ), USER_MINATOMS_OPTION.option, USER_MINATOMS_OPTION.xscheme); //CO 170601
      USER_MINATOMS_RESTICTED_OPTION=aurostd::substring2bool(USER_MINATOMS_OPTION.keyword,"_RESTRICTED");
      USER_MINATOMS = USER_MINATOMS_OPTION.content_int;
      if(USER_MINATOMS_OPTION.isentry){
        logger << (USER_MINATOMS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MINATOMS" << (USER_MINATOMS_RESTICTED_OPTION?string("_RESTRICTED"):string("")) << "=" << USER_MINATOMS << "." << apl::endl;
        logger << "Supercell will be built with at least " << USER_MINATOMS << " atoms." << apl::endl;
        found_supercell = true;
      }
    }

    // Get the users maximum shell which will be included into calculation
    // CO, not sure how maxshell works here (F option), need to investigate further and add to README
    // also seems USER_WANTS_FULL_SHELL applies for both MAX and MIN shell settings, should one take precedence? should they be separate flags?
    if(!found_supercell){
      USER_MAXSHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "MAXSHELL=" + "|" + _ASTROPT_QHA_ + "MAXSHELL=" + "|" + _ASTROPT_AAPL_ + "MAXSHELL=" + "|" + _ASTROPT_APL_OLD_ + "MAXSHELL="), USER_MAXSHELL_OPTION.option, USER_MAXSHELL_OPTION.xscheme); //CO 170601
      test = USER_MAXSHELL_OPTION.content_string;
      if (test[test.size() - 1] == 'f' || test[test.size() - 1] == 'F') {
        USER_MAXSHELL = aurostd::string2utype<int>(test.substr(0, test.size() - 1));
        USER_WANTS_FULL_SHELL = true;
      } else {
        USER_MAXSHELL = USER_MAXSHELL_OPTION.content_int;
        USER_WANTS_FULL_SHELL = false;
      }
      if(USER_MAXSHELL_OPTION.isentry){
        logger << (USER_MAXSHELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MAXSHELL=" << USER_MAXSHELL << (USER_WANTS_FULL_SHELL ? " (FULL)" : "") << "." << apl::endl;
        logger << "Supercell will be built with at most " << USER_MAXSHELL << " shells." << apl::endl;
        found_supercell = true;
      }
    }

    // Get the users minimum shell which will be included into calculation
    // CO, not sure how minshell works here (F option), need to investigate further and add to README
    if(!found_supercell){
      USER_MINSHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "MINSHELL=" + "|" + _ASTROPT_QHA_ + "MINSHELL=" + "|" + _ASTROPT_AAPL_ + "MINSHELL=" + "|" + _ASTROPT_APL_OLD_ + "MINSHELL="), USER_MINSHELL_OPTION.option, USER_MINSHELL_OPTION.xscheme); //CO 170601
      test = USER_MINSHELL_OPTION.content_string;
      if (test[test.size() - 1] == 'f' || test[test.size() - 1] == 'F') {
        USER_MINSHELL = aurostd::string2utype<int>(test.substr(0, test.size() - 1));
        USER_WANTS_FULL_SHELL = true;
      } else {
        USER_MINSHELL = USER_MINSHELL_OPTION.content_int;
        USER_WANTS_FULL_SHELL = false;
      }
      if(USER_MINSHELL_OPTION.isentry){
        logger << (USER_MINSHELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MINSHELL=" << USER_MINSHELL << (USER_WANTS_FULL_SHELL ? " (FULL)" : "") << "." << apl::endl;
        logger << "Supercell will be built with at least " << USER_MINSHELL << " shells." << apl::endl;
        found_supercell = true;
      }
    }

    // Get the users KPPRA which will be included into calculation, get from AFLOW machinery
    if(xflags.AFLOW_MODE_VASP){
      if (xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry) {logger << "Overriding with " << _ASTROPT_ << "KPPRA=" << xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.content_uint << "." << apl::endl;}
    // Get the users KSCHEME which will be included into calculation, get from AFLOW machinery
      if (xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry) {logger << "Overriding with " << _ASTROPT_ << "KSCHEME=" << xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.content_string << "." << apl::endl;}
      if (xflags.vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.flag("EVEN")) {logger << "Overriding with " << _ASTROPT_ << "KPOINTS=EVEN" << "." << apl::endl;}
      if (xflags.vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.flag("ODD")) {logger << "Overriding with " << _ASTROPT_ << "KPOINTS=ODD" << "." << apl::endl;}
    }

    // ADDITIONAL OPTIONS ------------------------------------------------

    // DC. e.g., DC = yes
    CALCULATE_PHONON_DISPERSIONS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DC=" + "|" + _ASTROPT_QHA_ + "DC=" + "|" + _ASTROPT_AAPL_ + "DC=" + "|" + _ASTROPT_APL_OLD_ + "DC="), CALCULATE_PHONON_DISPERSIONS_OPTION.option, CALCULATE_PHONON_DISPERSIONS_OPTION.xscheme); //CO 170601
    logger << (CALCULATE_PHONON_DISPERSIONS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DC=" << (CALCULATE_PHONON_DISPERSIONS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The calculation of PHONON DISPERSION curves is switched " << (CALCULATE_PHONON_DISPERSIONS_OPTION.option? "ON" : "OFF" ) << "." << apl::endl;

    if (CALCULATE_PHONON_DISPERSIONS_OPTION.option) {
      //CO 180406 - only choose one way to define path, INITLATTICE/INITSG vs. INITCOORDS_FRAC/INITCOORDS_CART + INITCOORDS_LABELS
      bool found_user_path=false;
      
      // DCINITLATTICE, e.g.,       DCINITLATTICE = RHL
      // DCINITSSG,  e.g., DCINITSG = 166
      // DCINITCOORDSFRAC, e.g.,    DCINITCOORDSFRAC = 0,0,0;0.5,0.5,0.5
      // DCINITCOORDSCART, e.g.,    DCINITCOORDSCART = 0,0,0;0.5,0.5,0.5
      // DCINITCOORDSLABELS, e.g.,  DCINITCOORDSLABELS = G,L,M,X
      // defaults are empty strings, so don't set unless isentry
      USER_DC_INITLATTICE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITLATTICE=" + "|" + _ASTROPT_QHA_ + "DCINITLATTICE=" + "|" + _ASTROPT_AAPL_ + "DCINITLATTICE=" + "|" + _ASTROPT_APL_OLD_ + "DCINITLATTICE="), USER_DC_INITLATTICE_OPTION.option, USER_DC_INITLATTICE_OPTION.xscheme); //CO 170601
      USER_DC_INITSG_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITSG=" + "|" + _ASTROPT_QHA_ + "DCINITSG=" + "|" + _ASTROPT_AAPL_ + "DCINITSG=" + "|" + _ASTROPT_APL_OLD_ + "DCINITSG="), USER_DC_INITSG_OPTION.option, USER_DC_INITSG_OPTION.xscheme); //CO 170601
      USER_DC_INITCOORDS_FRAC_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITCOORDSFRAC=" + "|" + _ASTROPT_QHA_ + "DCINITCOORDSFRAC=" + "|" + _ASTROPT_AAPL_ + "DCINITCOORDSFRAC=" + "|" + _ASTROPT_APL_OLD_ + "DCINITCOORDSFRAC="), USER_DC_INITCOORDS_FRAC_OPTION.option, USER_DC_INITCOORDS_FRAC_OPTION.xscheme); //CO 170601
      USER_DC_INITCOORDS_CART_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITCOORDSCART=" + "|" + _ASTROPT_QHA_ + "DCINITCOORDSCART=" + "|" + _ASTROPT_AAPL_ + "DCINITCOORDSCART=" + "|" + _ASTROPT_APL_OLD_ + "DCINITCOORDSCART="), USER_DC_INITCOORDS_CART_OPTION.option, USER_DC_INITCOORDS_CART_OPTION.xscheme); //CO 170601
      USER_DC_INITCOORDS_LABELS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITCOORDSLABELS=" + "|" + _ASTROPT_QHA_ + "DCINITCOORDSLABELS=" + "|" + _ASTROPT_AAPL_ + "DCINITCOORDSLABELS=" + "|" + _ASTROPT_APL_OLD_ + "DCINITCOORDSLABELS="), USER_DC_INITCOORDS_LABELS_OPTION.option, USER_DC_INITCOORDS_LABELS_OPTION.xscheme); //CO 170601
      if(!found_user_path && USER_DC_INITLATTICE_OPTION.isentry) {
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITLATTICE = USER_DC_INITLATTICE_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCINITLATTICE=" << USER_DC_INITLATTICE << " (via DCINITLATTICE)." << apl::endl;
        found_user_path=true;
      }  //corey
      if(!found_user_path && USER_DC_INITSG_OPTION.isentry){
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_DCINITSG = USER_DC_INITSG_OPTION.content_int;
        USER_DC_INITLATTICE = LATTICE::SpaceGroup2LatticeVariation(USER_DC_DCINITSG, xinput.getXStr()); //xvasp.str);
        logger << "Setting " << _ASTROPT_ << "DCINITLATTICE=" << USER_DC_INITLATTICE << " (via DCINITSG)." << apl::endl;
        found_user_path=true;
      }  //corey
      if(!found_user_path && USER_DC_INITCOORDS_FRAC_OPTION.isentry){
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITCOORDS_FRAC = USER_DC_INITCOORDS_FRAC_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCINITCOORDSFRAC=" << USER_DC_INITCOORDS_FRAC << "." << apl::endl;
        logger << "User's q-point path will be calculated along fractional coordinates [" << USER_DC_INITCOORDS_FRAC << "]." << apl::endl;
        found_user_path=true;
      }  //corey
      if(!found_user_path && USER_DC_INITCOORDS_CART_OPTION.isentry){
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITCOORDS_CART = USER_DC_INITCOORDS_CART_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCINITCOORDSCART=" << USER_DC_INITCOORDS_CART << "." << apl::endl;
        logger << "User's q-point path will be calculated along cartesian coordinates [" << USER_DC_INITCOORDS_CART << "]." << apl::endl;
        found_user_path=true;
      }  //corey

      if(!found_user_path){
        USER_DC_INITLATTICE.clear();
        USER_DC_INITCOORDS_FRAC.clear();
        USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITCOORDS_LABELS.clear();
      }

      if(found_user_path && ( USER_DC_INITCOORDS_FRAC_OPTION.isentry || USER_DC_INITCOORDS_CART_OPTION.isentry )){
        USER_DC_INITCOORDS_LABELS = USER_DC_INITCOORDS_LABELS_OPTION.content_string;
        if(USER_DC_INITCOORDS_LABELS.empty()){throw apl::APLRuntimeError("INITCOORDS set but no corresponding labels (DCINITCOORDSLABELS) found.");}
        logger << "Setting " << _ASTROPT_ << "DCUSERPATHLABELS=" << USER_DC_INITCOORDS_LABELS << "." << apl::endl;
        logger << "User's q-point path will be labeled [" << USER_DC_INITCOORDS_LABELS << "]." << apl::endl;
      }

      // DCUSERPATH, e.g., DCUSERPATH = G-X|X-U|K-G|G-L
      // default is empty string, so don't set unless isentry
      USER_DC_USERPATH_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCUSERPATH=" + "|" + _ASTROPT_QHA_ + "DCUSERPATH=" + "|" + _ASTROPT_AAPL_ + "DCUSERPATH=" + "|" + _ASTROPT_APL_OLD_ + "DCUSERPATH="), USER_DC_USERPATH_OPTION.option, USER_DC_USERPATH_OPTION.xscheme); //CO 170601
      if(USER_DC_USERPATH_OPTION.isentry){
        USER_DC_USERPATH = USER_DC_USERPATH_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCUSERPATH=" << USER_DC_USERPATH << "." << apl::endl;
        logger << "User's q-point path will be calculated along [" << USER_DC_USERPATH << "]." << apl::endl;
      }
      
      //ASSUME SINGLE POINT OTHERWISE
      //if we initcoords, we also need to specify the path
      //if(found_user_path && ( USER_DC_INITCOORDS_FRAC_OPTION.isentry || USER_DC_INITCOORDS_CART_OPTION.isentry )){
      //  if(USER_DC_USERPATH.empty()){throw apl::APLRuntimeError("INITCOORDS set but no corresponding path (DCUSERPATH) found.");}
      //}

      // DCPOINTS, e.g., DCPOINTS = 100
      USER_DC_NPOINTS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCPOINTS=" + "|" + _ASTROPT_QHA_ + "DCPOINTS=" + "|" + _ASTROPT_AAPL_ + "DCPOINTS=" + "|" + _ASTROPT_APL_OLD_ + "DCPOINTS="), USER_DC_NPOINTS_OPTION.option, USER_DC_NPOINTS_OPTION.xscheme); //CO 170601
      USER_DC_NPOINTS = USER_DC_NPOINTS_OPTION.content_int;
      logger << (USER_DC_NPOINTS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DCPOINTS=" << USER_DC_NPOINTS << "." << apl::endl;
      logger << "Each subpath will be divided into " << USER_DC_NPOINTS << " points." << apl::endl;
    }

    // DOS, e.g., DOS = yes
    CALCULATE_PHONON_DOS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOS=" + "|" + _ASTROPT_QHA_ + "DOS=" + "|" + _ASTROPT_AAPL_ + "DOS=" + "|" + _ASTROPT_APL_OLD_ + "DOS="), CALCULATE_PHONON_DOS_OPTION.option, CALCULATE_PHONON_DOS_OPTION.xscheme); //CO 170601
    logger << (CALCULATE_PHONON_DOS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOS=" << (CALCULATE_PHONON_DOS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The calculation of PHONON DENSITY of states switched " << (CALCULATE_PHONON_DOS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

    if (CALCULATE_PHONON_DOS_OPTION.option){
      // DOSMESH, e.g., DOSMESH = 20x20x20
      USER_DOS_MESH_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSMESH=" + "|" + _ASTROPT_QHA_ + "DOSMESH=" + "|" + _ASTROPT_AAPL_ + "DOSMESH=" + "|" + _ASTROPT_APL_OLD_ + "DOSMESH="), USER_DOS_MESH_OPTION.option, USER_DOS_MESH_OPTION.xscheme); //CO 170601
      USER_DOS_MESH = USER_DOS_MESH_OPTION.content_string;
      tokens.clear();
      apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"DOSMESH. Specify as DOSMESH=20x20x20");}
      logger << (USER_DOS_MESH_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOSMESH=" << USER_DOS_MESH << "." << apl::endl;

      // DOSPOINTS, e.g., DOSPOINTS = 2000
      USER_DOS_NPOINTS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSPOINTS=" + "|" + _ASTROPT_QHA_ + "DOSPOINTS=" + "|" + _ASTROPT_AAPL_ + "DOSPOINTS=" + "|" + _ASTROPT_APL_OLD_ + "DOSPOINTS="), USER_DOS_NPOINTS_OPTION.option, USER_DOS_NPOINTS_OPTION.xscheme); //CO 170601
      USER_DOS_NPOINTS = USER_DOS_NPOINTS_OPTION.content_int;
      logger << (USER_DOS_NPOINTS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "USER_DOS_NPOINTS=" << USER_DOS_NPOINTS << "." << apl::endl;
      logger << "The phonon density of states will be calculated for " << USER_DOS_NPOINTS << " bins." << apl::endl;

      // DOSMETHOD, e.g., DOSMETHOD = LT
      USER_DOS_METHOD_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSMETHOD=" + "|" + _ASTROPT_QHA_ + "DOSMETHOD=" + "|" + _ASTROPT_AAPL_ + "DOSMETHOD=" + "|" + _ASTROPT_APL_OLD_ + "DOSMETHOD="), USER_DOS_METHOD_OPTION.option, USER_DOS_METHOD_OPTION.xscheme); //CO 170601
      USER_DOS_METHOD = USER_DOS_METHOD_OPTION.content_string;
      transform(USER_DOS_METHOD.begin(), USER_DOS_METHOD.end(), USER_DOS_METHOD.begin(), toupper);
      logger << (USER_DOS_METHOD_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOSMETHOD=" << USER_DOS_METHOD << "." << apl::endl;
      if (USER_DOS_METHOD == string("LT")) {logger << "The phonon density of states will be calculated by Linear Tetrahedron Method (LT)." << apl::endl;}
      else if (USER_DOS_METHOD == string("RS")){logger << "The phonon density of states will be calculated by Root Sampling Method (RS)." << apl::endl;}
      else{throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"DOSMETHOD. Specify as DOSMETHOD=LT.");}

      // DOSSMEAR, e.g., DOSSMEAR = 0.05
      USER_DOS_SMEAR_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSSMEAR=" + "|" + _ASTROPT_QHA_ + "DOSSMEAR=" + "|" + _ASTROPT_AAPL_ + "DOSSMEAR=" + "|" + _ASTROPT_APL_OLD_ + "DOSSMEAR="), USER_DOS_SMEAR_OPTION.option, USER_DOS_SMEAR_OPTION.xscheme); //CO 170601
      USER_DOS_SMEAR = USER_DOS_SMEAR_OPTION.content_double;
      override_option=false;
      if (USER_DOS_METHOD == string("RS")) {USER_DOS_SMEAR = 0.05; override_option=true;} // Default value, it is better with this...
      logger << ((USER_DOS_SMEAR_OPTION.isentry || override_option) ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOSSMEAR=" << USER_DOS_SMEAR << " (overridden by RS DOSMETHOD)." << apl::endl;
      if (USER_DOS_SMEAR > 1E-6){logger << "The phonon density of states will be smooth by gaussians with sigma = " << USER_DOS_SMEAR << "." << apl::endl;}
    }

    // TP, e.g., TP = yes
    CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "TP=" + "|" + _ASTROPT_QHA_ + "TP=" + "|" + _ASTROPT_AAPL_ + "TP=" + "|" + _ASTROPT_APL_OLD_ + "TP="), CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option, CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.xscheme); //CO 170601
    logger << (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "TP=" << (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The calculation of THERMODYNAMIC properties switched " << (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    if (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option && !CALCULATE_PHONON_DOS_OPTION.option){logger << apl::warning << "The thermodynamic properties may be calculated by default settings for DOS calculation." << apl::endl;}

    // TPT, e.g., TPT = 1000:2000:10 -> temperature from 1000 to 2000 K by step 10 K
    if (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option){
      USER_TPT_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "TPT=" + "|" + _ASTROPT_QHA_ + "TPT=" + "|" + _ASTROPT_AAPL_ + "TPT=" + "|" + _ASTROPT_APL_OLD_ + "TPT="), USER_TPT_OPTION.option, USER_TPT_OPTION.xscheme); //CO 170601
      tokens.clear();
      apl::tokenize(USER_TPT_OPTION.content_string, tokens, string(" :"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"TPT. Specify as TPT=1000:2000:10.");}
      USER_TP_TSTART = aurostd::string2utype<double>(tokens.at(0));
      USER_TP_TEND = aurostd::string2utype<double>(tokens.at(1));
      USER_TP_TSTEP = aurostd::string2utype<double>(tokens.at(2));
      logger << (USER_TPT_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "TPT=" << USER_TPT_OPTION.content_string << "." << apl::endl;
      logger << "The thermodynamic properties will be calculated in temperature range <" << USER_TP_TSTART << "," << USER_TP_TEND << "> with step " << USER_TP_TSTEP << " K." << apl::endl;
    }

  } catch (std::exception& e) {
    logger << apl::error << e.what() << apl::endl;
    return;
  }

  // ///////////////////////////////////////////////////////////////////////

  try {
    // Contruct the working supercell ////////////////////////////////////

    //apl::Supercell supercell(xvasp.str,logger),supercell_test(xvasp.str,logger);    //corey, slow
    apl::Supercell supercell(xinput.getXStr(),logger); //xvasp.str, logger);  //CO
    apl::Supercell supercell_test = supercell;    //CO

    //   pflow::PrintDist(xinput.getXStr(),20.0,cerr);
    if (USER_SUPERCELL.empty() && USER_MINATOMS > 0) {
      stringstream aus;
      if(USER_MINATOMS_RESTICTED_OPTION){
        for (int Ni=1; USER_SUPERCELL == ""; Ni++) {
          aus.str("");
          aus << "Ni=" << Ni
              << " "
              << "supercell=" << Ni << "x" << Ni << "x" << Ni << "  natoms=" << Ni * Ni * Ni * xinput.getXStr().atoms.size(); //xvasp.str.atoms.size();
          //	  logger << aus.getXStr() << apl::endl;
          if (Ni * Ni * Ni * ((int)xinput.getXStr().atoms.size()) > (int)USER_MINATOMS) { // xvasp.str.atoms.size()) > (int)USER_MINATOMS) {
            USER_MINATOMS = 0;
            USER_SUPERCELL = aurostd::utype2string<uint>(Ni) + "X" + aurostd::utype2string<uint>(Ni) + "X" + aurostd::utype2string<uint>(Ni);
            logger << aus.str() << apl::endl;
          }
        }
      }else{
      for (double radius = 0.01; USER_SUPERCELL == ""; radius += 0.01) {
        xvector<int> dims(3);
        dims = LatticeDimensionSphere(xinput.getXStr().lattice,radius); //xvasp.str.lattice, radius);
          aus.str("");
        aus << "Radius=" << aurostd::PaddedPOST(aurostd::utype2string<double>(radius, 3), 4)
            << " "
            << " supercell=" << dims(1) << "x" << dims(2) << "x" << dims(3) << "  natoms=" << dims(1) * dims(2) * dims(3) * xinput.getXStr().atoms.size(); //xvasp.str.atoms.size();
        //	  logger << aus.getXStr() << apl::endl;
        if (dims(1) * dims(2) * dims(3) * ((int)xinput.getXStr().atoms.size()) > (int)USER_MINATOMS) { // xvasp.str.atoms.size()) > (int)USER_MINATOMS) {
          USER_MINATOMS = 0;
          USER_SUPERCELL = aurostd::utype2string<uint>(dims(1)) + "X" + aurostd::utype2string<uint>(dims(2)) + "X" + aurostd::utype2string<uint>(dims(3));
          logger << aus.str() << apl::endl;
          }
        }
      }
    }

    //      cerr << "USER_WANTS_FULL_SHELL=" << USER_WANTS_FULL_SHELL << endl;
    //    for(int i=2;i<20;i+=2) cerr << "try " << i << ": " << supercell_test.buildSuitableForShell(i,USER_WANTS_FULL_SHELL,FALSE) << endl;

    if (USER_SUPERCELL.empty() && USER_MAXSHELL > 0) {
      logger << "a Searching for suitable cell to handle " << USER_MINSHELL << " shells..." << apl::endl;
      supercell.buildSuitableForShell(USER_MAXSHELL, USER_WANTS_FULL_SHELL, TRUE);
      supercell.setupShellRestrictions(USER_MAXSHELL);
    } else if (USER_SUPERCELL.empty() && USER_MINSHELL > 0) {
      logger << "b Searching for suitable cell to handle " << USER_MINSHELL << " shells..." << apl::endl;
      supercell.buildSuitableForShell(USER_MINSHELL, USER_WANTS_FULL_SHELL, TRUE);
    } else if (USER_SUPERCELL.find_first_of("xX") != string::npos) {
      // OK, user wants its own supercell...
      tokens.clear();
      apl::tokenize(USER_SUPERCELL, tokens, string(" xX"));
      supercell.build(aurostd::string2utype<int>(tokens.at(0)),
                      aurostd::string2utype<int>(tokens.at(1)),
                      aurostd::string2utype<int>(tokens.at(2)));
      // Did he specify also regular restriction for max shell included
      // in calculation?
      if (USER_MAXSHELL > 0)
        supercell.setupShellRestrictions(USER_MAXSHELL);
    } else {
      throw apl::APLRuntimeError("The settings for supercell construction are confusing.");
    }

    // ***** BEGIN JJPR ******
    // Supercell distances
    // Calculate inequivalent pairs and triplets for thermal conductivity purpose------------------------------------
    apl::StrPairs strPair(logger);  //CO, default, does NOTHING
    if (CALCULATE_KAPPA_OPTION.option) {
      supercell.calcRl();
      apl::StrPairs _strPair(supercell.getSupercellStructure(), supercell.getInputStructure(), supercell._pc2scMap, supercell._sc2pcMap, USER_CUTOFF_SHELL, USER_CUTOFF_DISTANCE, logger);
      strPair = _strPair;  //CO, copies everything
      USER_CUTOFF_DISTANCE = strPair.getShellDist(supercell.getSupercellStructure(), USER_CUTOFF_SHELL, USER_CUTOFF_DISTANCE);
      strPair.build(supercell.getInputStructure(), supercell.getSupercellStructure());
      strPair.Triplets();
    }
    // ******  END JJPR ******

    // Calculate phonons /////////////////////////////////////////////////

    auto_ptr<apl::PhononCalculator> phcalc;
    if (USER_ENGINE == string("DM")) {
      apl::DirectMethodPC* phcalcdm = new apl::DirectMethodPC(supercell, strPair, xinput, aflags, kflags, xflags, AflowIn, logger); //xvasp, aflags, kflags, vflags, logger);  //Modified JJPR
      phcalcdm->isPolarMaterial(CALCULATE_POLAR_CORRECTIONS_OPTION.option);                                                       // TRY POLAR [STEFANO]
      phcalcdm->setTensor(CALCULATE_KAPPA_OPTION.option);                                                                         // KAPPA JJPR
      phcalcdm->setSumRule(USER_EPS_SUM);                                                                                  // KAPPA JJPR
      //phcalcdm->setGeneratePlusMinus(USER_DISTORTIONS_PLUS_MINUS_OPTION.option); //CO auto
      phcalcdm->setGeneratePlusMinus(AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option, USER_DISTORTIONS_PLUS_MINUS_OPTION.option);  //CO auto
      phcalcdm->setGenerateOnlyXYZ(USER_DISTORTIONS_XYZ_ONLY_OPTION.option);
      phcalcdm->setDistortionMagnitude(USER_DISTORTION_MAGNITUDE);
      phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
      phcalcdm->get_special_inputs(AflowIn);  //PINKU, to include PSTRESS and LDAU_PARAMETERS in the SUPERCELL files
      phcalc.reset(phcalcdm);
    //CO generally redirects to DM, the distinction between DM and GSA is obsolete
    //} else if (USER_ENGINE == string("GSA")) {
    //  apl::GeneralizedSupercellApproach* gsa = new apl::GeneralizedSupercellApproach(supercell, strPair, xinput, aflags, kflags, xflags, logger);//xvasp, aflags, kflags, vflags, logger);  //Modified JJPR
    //  //gsa->setGeneratePlusMinus(USER_DISTORTIONS_PLUS_MINUS_OPTION.option); //CO auto
    //  gsa->setGeneratePlusMinus(AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option, USER_DISTORTIONS_PLUS_MINUS_OPTION.option);  //CO auto
    //  gsa->setGenerateOnlyXYZ(USER_DISTORTIONS_XYZ_ONLY_OPTION.option);
    //  gsa->setDistortionMagnitude(USER_DISTORTION_MAGNITUDE);
    //  gsa->setTensor(CALCULATE_KAPPA_OPTION.option);  // KAPPA JJPR
    //  gsa->setSumRule(USER_EPS_SUM);           // KAPPA JJPR
    //  //phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
    //  phcalc.reset(gsa);
    } else if (USER_ENGINE == string("LR")) {
      phcalc.reset(new apl::LinearResponsePC(supercell, strPair, xinput, aflags, kflags, xflags, AflowIn, logger)); //xvasp, aflags, kflags, vflags, logger));  //Modified JJPR
      phcalc->setTensor(CALCULATE_KAPPA_OPTION.option);                                                           // KAPPA JJPR
      phcalc->setSumRule(USER_EPS_SUM);                                                                    // KAPPA JJPR
      phcalc->isPolarMaterial(CALCULATE_POLAR_CORRECTIONS_OPTION.option);
      //phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
    } else
      throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"ENGINE. Set DM or LR only.");

    //PINKU QUASI-HARMONIC START
    // To create directories for Gruneisen and EOS calculations
    // it should be called before creation of apl.xml
    auto_ptr<apl::eos_calculate> pheos;
    if (CALCULATE_GRUNEISEN_OPTION.option) {
      pheos.reset(new apl::eos_calculate(supercell, strPair, xinput, aflags, kflags, xflags, AflowIn, logger)); //xvasp, aflags, kflags, vflags, logger));
      pheos->setGruneisen(CALCULATE_GRUNEISEN_OPTION.option);
      pheos->setGP_VOL_DISTORTION(GP_VOL_DISTORTION);
      pheos->setEOS(CALCULATE_EOS_OPTION.option);
      pheos->GP_RUN();
      if (CALCULATE_EOS_OPTION.option) {
        pheos->setTEC_VOL_DISTORTION(EOS_VOL_DISTORTION);
        pheos->setEOS_VOL_START(EOS_VOL_START);
        pheos->setEOS_VOL_END(EOS_VOL_END);
        pheos->setEOS_KPPRA(EOS_KPPRA);
        pheos->setEOS_STATIC_KPPRA(EOS_STATIC_KPPRA);
        pheos->setEOS_BANDS_GRID(EOS_BANDS_GRID);
        pheos->setNEDOS(NEDOS);
        pheos->EOS_RUN(AflowIn,EOS_VOL_START, EOS_VOL_END, EOS_VOL_DISTORTION);
      }
    }
    //PINKU QUASI-HARMONIC END

    // Run or awake
    bool isHibFileAvalaible = aurostd::EFileExist(string("apl.xml"));  //|| //CO
    //aurostd::FileExist(string("apl.xml")); //CO

    if (USER_HIBERNATE_OPTION.option && isHibFileAvalaible) {
      try {
        phcalc->awake();
      } catch (apl::APLLogicError& e) {
        logger << apl::warning << e.what() << apl::endl;
        logger << apl::warning << "Skipping awakening..." << apl::endl;
        isHibFileAvalaible = false;
      }
    }

    if (!isHibFileAvalaible) {
      phcalc->run();
      if (USER_HIBERNATE_OPTION.option)
        phcalc->hibernate();
    }
    //PINKU QUASI-HARMONIC START
    //store dynamical matrices and phonon disperson curves from all
    //APL_PHONON directories to compute Gruneisen parameter and EOSs
    if (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option) {
      apl::DM_PDOS_save store(*phcalc, logger);
      store.setdir_prefix(_TMPDIR_, GP_VOL_DISTORTION);
      string dirname = store.getdir_name(aflags.Directory);
      if (store.check_GP()) {  //Check directory is Gruneisen ?
        tokens.clear();
        apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
        //create uniform q-mesh
        store.createMPmesh(aurostd::string2utype<int>(tokens.at(0)),
                           aurostd::string2utype<int>(tokens.at(1)),
                           aurostd::string2utype<int>(tokens.at(2)),
                           phcalc->getInputCellStructure());

        apl::PhononDispersionCalculator pdisc(*phcalc, logger);
        
        // Init path according to the aflow's definition for elec. struc.
        if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_FRAC.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,false);}
        else if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_CART.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_CART,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,true);}
        else{pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);} //default!
        if(!USER_DC_USERPATH.empty()){pdisc.setPath(USER_DC_USERPATH);}   //Does user want his own path?
        std::vector<xvector<double> > qpoints = pdisc.get_qpoints();
        store.create_pdispath(qpoints);
        qpoints.clear();
        pdisc.clear();
      }
      if (CALCULATE_EOS_SUBDIRECTORIES_OPTION.option) {  //Check directory is EOS ?
        tokens.clear();
        apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
        apl::MonkhorstPackMesh qmesh(
            aurostd::string2utype<int>(tokens.at(0)),
            aurostd::string2utype<int>(tokens.at(1)),
            aurostd::string2utype<int>(tokens.at(2)),
            phcalc->getInputCellStructure(), logger);

        auto_ptr<apl::DOSCalculator> dosc;

        if (USER_DOS_METHOD == string("LT"))
          dosc.reset(new apl::LinearTetrahedronMethod(*phcalc, qmesh, logger));
        else if (USER_DOS_METHOD == string("RS"))
          dosc.reset(new apl::RootSamplingMethod(*phcalc, qmesh, logger));
        else
          throw apl::APLRuntimeError("Unknown DOS method. Check "+_ASTROPT_+"DOSMETHOD command.");
        // Calculate DOS
        dosc->calc(USER_DOS_NPOINTS, USER_DOS_SMEAR);
        if (CALCULATE_PHONON_DOS_OPTION.option) dosc->writePDOS(_TMPDIR_, dirname);
        dosc->clear();
      }
      store.clear();
      phcalc->clear();
      return;
    }
    //PINKU QUASI-HARMONIC END

    // Get the format of frequency desired by user ///////////////////////

    apl::IPCFreqFlags frequencyFormat = apl::NONE;

    if (!USER_FREQFORMAT.empty()) {
      // Convert format to machine representation
      tokens.clear();
      apl::tokenize(USER_FREQFORMAT, tokens, string(" |:;,"));
      for (uint i = 0; i < tokens.size(); i++) {
        if (tokens.at(i) == string("OMEGA")) {
          frequencyFormat |= apl::OMEGA;
          continue;
        }
        if (tokens.at(i) == string("HERTZ")) {
          frequencyFormat |= apl::HERTZ;
          continue;
        } else if (tokens.at(i) == string("THZ")) {
          frequencyFormat |= apl::THZ;
          continue;
        } else if (tokens.at(i) == string("CM-1") || tokens.at(i) == string("RECIPROCAL_CM")) {
          frequencyFormat |= apl::RECIPROCAL_CM;
          continue;
        } else if (tokens.at(i) == string("MEV")) {
          frequencyFormat |= apl::MEV;
          continue;
        }
        if (tokens.at(i) == string("ALLOW_NEGATIVE")) {
          frequencyFormat |= apl::ALLOW_NEGATIVE;
          continue;
        }
      }
      // Check if there was specified unit keyword...
      if (((frequencyFormat & ~apl::OMEGA) & ~apl::ALLOW_NEGATIVE) == apl::NONE)
        throw apl::APLLogicError("The mishmash frequency format.");
    } else {
      frequencyFormat = apl::THZ | apl::ALLOW_NEGATIVE;
    }

    //////////////////////////////////////////////////////////////////////

    if (CALCULATE_PHONON_DISPERSIONS_OPTION.option) {
      apl::PhononDispersionCalculator pdisc(*phcalc, logger);

      // Init path according to the aflow's definition for elec. struc.
      if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_FRAC.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,false);}
      else if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_CART.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_CART,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,true);}
      else{pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);} //default!
      if(!USER_DC_USERPATH.empty()){pdisc.setPath(USER_DC_USERPATH);}   //Does user want his own path?

      // Calculate frequencies on path
      pdisc.calc(frequencyFormat);

      // Write results into PDIS file
      pdisc.writePDIS();
      //PINKU QUASI-HARMONIC START
      //compute Gruneisen dispersion curve
      if (CALCULATE_GRUNEISEN_OPTION.option) {
        std::vector<xvector<double> > qpoints = pdisc.get_qpoints();
        vector<double> path = pdisc.get_path();
        vector<int> path_segment = pdisc.get_path_segment();
        apl::EOS eos(*phcalc, *pheos, logger);
        apl::QuasiHarmonicGruneisen* qh = &eos;
        if (qh->set_directories(_TMPDIR_)) {
          qh->set_frequency_cutoff(CUTOFF_FREQ);
          if (qh->cal_gp_along_path(qpoints))
            qh->write_GP_path(path, path_segment);
        }
        //clear used variables
        path.clear();
        path_segment.clear();
        qpoints.clear();
        qh->clear();
      }
      //PINKU QUASI-HARMONIC END
    }

    //////////////////////////////////////////////////////////////////////

    if (CALCULATE_PHONON_DOS_OPTION.option || CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option) {
      // Generate mesh for calculation of DOS...
      tokens.clear();
      apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
      apl::MonkhorstPackMesh qmesh(aurostd::string2utype<int>(tokens.at(0)),
                                   aurostd::string2utype<int>(tokens.at(1)),
                                   aurostd::string2utype<int>(tokens.at(2)),
                                   phcalc->getInputCellStructure(), logger);

      // Setup the DOS engine which is used also for thermodynamic properties
      auto_ptr<apl::DOSCalculator> dosc;
      if (USER_DOS_METHOD == string("LT"))
        dosc.reset(new apl::LinearTetrahedronMethod(*phcalc, qmesh, logger));
      else if (USER_DOS_METHOD == string("RS"))
        dosc.reset(new apl::RootSamplingMethod(*phcalc, qmesh, logger));
      else
        throw apl::APLRuntimeError("Unknown DOS method. Check "+_ASTROPT_+"DOSMETHOD command.");

      // Calculate DOS
      dosc->calc(USER_DOS_NPOINTS, USER_DOS_SMEAR);
      if (CALCULATE_PHONON_DOS_OPTION.option) {
        dosc->writePDOS();
      }

      // Calculate thermal properties
      if (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option) {
        if (!dosc->hasNegativeFrequencies()) {
          apl::ThermalPropertiesCalculator tpc(*dosc, logger);
          tpc.writeTHERMO(USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP);
          //PINKU QUASI-HARMONIC START
          //calculate Gruneisen
          if (CALCULATE_GRUNEISEN_OPTION.option) {
            apl::EOS eos(*phcalc, *pheos, logger);
            apl::QuasiHarmonicGruneisen* qh = &eos;
            if (qh->set_directories(_TMPDIR_)) {
              qh->set_frequency_cutoff(CUTOFF_FREQ);
              tokens.clear();
              apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
              qh->createmesh(aurostd::string2utype<int>(tokens.at(0)),
                             aurostd::string2utype<int>(tokens.at(1)),
                             aurostd::string2utype<int>(tokens.at(2)),
                             phcalc->getInputCellStructure());
              if (qh->cal_gp_in_mesh()) {
                qh->write_GP_mesh();
                qh->write_GP_FREQ_mesh();
                qh->calculate_acoustic_freqNgp();
                qh->WriteaverageGP(USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP);
                qh->Write_BrINDEXs();
                //atomic displacement calculations
                if (CALCULATE_DISPLACEMENTS_OPTION.option) {
                  qh->thermal_displacements(USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP);
                  if (directions.size() == 3)
                    qh->projected_displacement(directions, USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP);
                  else
                    throw apl::APLRuntimeError("Unknown PROJECTION_DIR format. Check "+_ASTROPT_+"PROJECTION_DIR command.");
                }
                //calculate EOS
                if (CALCULATE_EOS_OPTION.option) {
                  if (qh->setvariables()) {
                    eos.set_fitting_type(FITTING_TYPE);
                    eos.calEOS(USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP, tpc);
                  }
                  eos.clear();
                }
              }
              //else logger << apl::warning << "Some files are missing. The calculation of Gruneisen parameters  has been skipped." << apl::endl;
              qh->clear();
            }
            pheos->clear();
          }
          //PINKU QUASI-HARMONIC END
          tpc.clear();
        } else {
          logger << apl::warning << "There are negative frequencies in DOS. The calculation of thermal properties has been skipped." << apl::endl;
        }
      }

      // Clear old stuff
      dosc->clear();
      //delete dosc; //auto_ptr will do
      qmesh.clear();
    }

    //////////////////////////////////////////////////////////////////////
    // ***** BEGIN JJPR: THERMAL CONDUCTIVITY ******
    if (CALCULATE_KAPPA_OPTION.option) {
      //Constructor and store information for thermal conductivity----------------------------
      apl::ThermalConductivityCalculator tcond(*phcalc, supercell, strPair, logger);

      //Build GRID of q-points
      if (USER_THERMALGRID.find_first_of("xX") != string::npos) {
        logger << "Preparing a q-mesh of " << USER_THERMALGRID << " for thermal conductivity purpose" << apl::endl;
        tokens.clear();
        apl::tokenize(USER_THERMALGRID, tokens, string(" xX"));
        tcond.qgrid(aurostd::string2utype<int>(tokens.at(0)),
                    aurostd::string2utype<int>(tokens.at(1)),
                    aurostd::string2utype<int>(tokens.at(2)));

        logger << "Calculating freqs for thermal conductivity purpose" << apl::endl;

        //Obtain frecuencies and eigenvectors
        tokens.clear();
        apl::tokenize(USER_THERMALGRID, tokens, string(" xX"));
        tcond.calcFreq(aurostd::string2utype<int>(tokens.at(0)),
                       aurostd::string2utype<int>(tokens.at(1)),
                       aurostd::string2utype<int>(tokens.at(2)));

        //Obtain Group velocities----------------------------
        logger << "Calculating group velocity  for thermal conductivity purpose" << apl::endl;
        tcond.getGroupVelocity();

        //Compute Three phonon scattering processes----------------------------
        logger << "Calculating Scattering matrix for thermal conductivity purpose" << apl::endl;
        tcond.getScattering();
        strPair.clear();

        // Evaluate which method will be used for solve BTE-------------------------------
        bool RTA = false;
        if (USER_BTE == string("RTA")) { RTA = true; }
        tcond.Calculator(RTA, CALCULATE_ISOTOPE_OPTION.option, CALCULATE_CUMULATIVEK_OPTION.option, CALCULATE_BOUNDARY_OPTION.option, USER_NANO_SIZE, USER_TCT_TSTART, USER_TCT_TEND, USER_TCT_TSTEP);

      } else {
        throw apl::APLRuntimeError("The settings for q-points grid  construction (THERMALGRID) are confusing.");
      }
    }
    // ***** END JJPR: THERMAL CONDUCTIVITY ******

    phcalc->clear();
    // delete phcalc; //auto_ptr will do
    supercell.clear();
  } catch (apl::APLStageBreak& e) {
    logger << apl::notice << "Stopped. Waiting for finishing of required calculations..." << apl::endl;
  } catch (std::exception& e) {
    logger << apl::error << e.what() << apl::endl;
  }
}
}

//////////////////////////////////////////////////////////////////////////////

bool PHON_RunPhonons(const xstructure& _str,
                     _aflags& aflags,
                     const double& _radius,
                     const bool& osswrite, ostream& oss) {
  bool Krun = FALSE;
  oss << "not implemented" << endl;
  if (0) {  // to avoid warnings
    oss << _str << endl
        << aflags.QUIET << _radius << osswrite << endl;
  }
  return Krun;
}

// ***************************************************************************
// *                                                                         *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                         *
// ***************************************************************************
