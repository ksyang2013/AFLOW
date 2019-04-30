// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_H_
#define _AFLOW_H_

// --------------------------------------------------------------------------
// Standard stuff

#include "AUROSTD/aurostd.h"

// #define  _AFLOW_TEMP_PRESERVE_  // to preseve /tmp files for debug

#define NNN   -123456
#define GCC_VERSION (__GNUC__ * 10000  + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#define _ANRL_NOWEB_ // DX

//COMMON TOLERANCES
#define _ZERO_TOL_ 1e-10 // DX
#define _XPROTO_TOO_CLOSE_ERROR_ 0.60 // was 0.75

//XSTRUCTURE definitions
#define _AFLOW_XSTR_PRINT_PRECISION_ 14  //CO 180509
#define _AFLOW_POCC_PRECISION_ 8 //must be less than _AFLOW_XSTR_PRINT_PRECISION_, which is currently set to 14
#define _AFLOW_POCC_ZERO_TOL_ pow(10,-_AFLOW_POCC_PRECISION_) 

//moved from avasp.cpp for broader access (chull.cpp)
#define SPECIE_RAW_LIB2U string("Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")
#define SPECIE_RAW_LIB2 string("Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv")

#define SPECIE_RAW_LIB3 string("Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_sv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv")
//#define SPECIE_RAW_LIB3 string("Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")
//#define SPECIE_RAW_LIB3 string("Ag,Al,As,Au,B_h,Bi_d,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,La,Mg_pv,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Ta_pv,Te,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")

// CO 180729 - OBSOLETE - use xerror
//[OBSOLETE]// CO 180419 - global exception handling - START
//[OBSOLETE]class AFLOWRuntimeError : public std::runtime_error {
//[OBSOLETE]  public:
//[OBSOLETE]    AFLOWRuntimeError(const std::string& function,const std::string& message);
//[OBSOLETE]    AFLOWRuntimeError(const std::string& function,std::stringstream& message);
//[OBSOLETE]    string where();
//[OBSOLETE]    ~AFLOWRuntimeError() throw() {};
//[OBSOLETE]  private:
//[OBSOLETE]    string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//[OBSOLETE]};
//[OBSOLETE]class AFLOWLogicError : public std::logic_error {
//[OBSOLETE]  public:
//[OBSOLETE]    AFLOWLogicError(const std::string& function,const std::string& message);
//[OBSOLETE]    AFLOWLogicError(const std::string& function,std::stringstream& message);
//[OBSOLETE]    string where();
//[OBSOLETE]    ~AFLOWLogicError() throw() {};
//[OBSOLETE]  private:
//[OBSOLETE]    string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//[OBSOLETE]};
//[OBSOLETE]// CO 180419 - global exception handling - STOP

// --------------------------------------------------------------------------
// definitions for MULTHREADS
//#define MAX_ALLOCATABLE_PTHREADS     1024
#define MAX_ALLOCATABLE_PTHREADS     256
#define PTHREADS_DEFAULT 8
namespace AFLOW_PTHREADS {
  extern bool FLAG;         // run pthread YES/NO
  extern int MAX_PTHREADS;  // how many MAX threads I can use  default or --np
  extern int RUNNING;       // how many threads are actually running
  extern pthread_t vpthread[MAX_ALLOCATABLE_PTHREADS];  // the actual thread
  extern int viret[MAX_ALLOCATABLE_PTHREADS];          // the thread runnings
  extern bool vpthread_busy[MAX_ALLOCATABLE_PTHREADS];  // is the thread busy
}

extern string _AFLOWIN_; 
extern string _AFLOWLOCK_; 

// --------------------------------------------------------------------------
// definitions for aflow
// aflow2 default definitions
#define AFLOW_MATERIALS_SERVER_DEFAULT        string("materials.duke.edu")
#define AFLOW_WEB_SERVER_DEFAULT              string("nietzsche.mems.duke.edu")
#define AFLOWLIB_MATERIALS_SERVER             string("aflow.org")
#define AFLOWLIB_CONSORTIUM_STRING            string("AFLOW - www.aflow.org consortium")
#define _XENTRY_ string("index.php")
 
// [OBSOLETE] #define DEFAULT_KZIP_BIN              string("bzip2")           // moved to aflow_aflowrc.cpp in V3.1.194
// [OBSOLETE] #define DEFAULT_KZIP_EXT              string(".bz2")            // moved to aflow_aflowrc.cpp in V3.1.194
#define DEFAULT_KBIN_ALIEN_BIN        string("ls -las")
#define DEFAULT_KBIN_MATLAB_BIN       string("/usr/local/bin/matlab -nodesktop -nosplash -nodisplay ")
// [OBSOLETE] #define DEFAULT_KBIN_CONVERT_BIN      string("convert")
// [OBSOLETE] #define DEFAULT_KBIN_EPSTOPDF_BIN     string("epstopdf")

#define QSUB_COMMAND_DEFAULT          "qsub"
#define QSUB_PARAMS_DEFAULT           " "
  
#define KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT 3.0
#define KBIN_SYMMETRY_SGROUP_MAX_NUMBER 1000000

#define KBIN_SUBDIRECTORIES           string("ARUN.")

#define KBIN_NEIGHBOURS_MAX_NUMBER      30000
#define KBIN_NEIGHBOURS_RADIUS_DEFAULT 3.0
#define KBIN_NEIGHBOURS_DRADIUS_DEFAULT 0.1

#define ALIEN_INPUT_FILE_NAME_DEFAULT  "./input"
#define ALIEN_EXTERNAL_INPUT_DEFAULT   "../input_external"
#define ALIEN_OUTPUT_FILE_NAME_DEFAULT  "./output"

// aflow1 definitions (soon to be obsolete)
#define _MPI_NP_STRINGS_ "MPI_NP","mpi_np","-MPI_NP","-mpi_np"
#define _MPI_NCPUS_DEF_ 4
#define VASP_OPTIONS_MPI_DEFAULT         ""
#define VASPLS_BIN_POSTFIX_DEFAULT       "LS"
#define GRND_BIN_DEFAULT                 "./grnd_intel"

// --------------------------------------------------------------------------
// definitions for projects
// [OBSOLETE] #define SERVER_PROJECT_GNDSTATE       string("/common/GNDSTATE")
#define LIBRARY_NOTHING 256
extern uint LIBRARY_ICSD,LIBRARY_LIB1,LIBRARY_LIB2,LIBRARY_LIB3,LIBRARY_LIB4,LIBRARY_LIB5,LIBRARY_LIB6,LIBRARY_LIB7,LIBRARY_LIB8,LIBRARY_LIB9;  // not in order.. will be nailed by init.cpp

#define LIBRARY_ALL         100

// [OBSOLETE] #define DEFAULT_FILE_AFLOWLIB_ENTRY_OUT        string("aflowlib.out")      // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_AFLOWLIB_ENTRY_JSON       string("aflowlib.json")     // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_ORIG_OUT            string("edata.orig.out")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_RELAX_OUT           string("edata.relax.out")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_BANDS_OUT           string("edata.bands.out")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_ORIG_OUT             string("data.orig.out")     // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_RELAX_OUT            string("data.relax.out")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_BANDS_OUT            string("data.bands.out")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_ORIG_JSON           string("edata.orig.json")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_RELAX_JSON          string("edata.relax.json")  // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_BANDS_JSON          string("edata.bands.json")  // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_ORIG_JSON            string("data.orig.json")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_RELAX_JSON           string("data.relax.json")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_BANDS_JSON           string("data.bands.json")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_TIME_OUT                  string("time")              // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_SPACEGROUP1_OUT           string("SpaceGroup")        // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_SPACEGROUP2_OUT           string("SpaceGroup2")       // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_VOLDISTPARAMS_OUT         string("VOLDISTParams")     // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_VOLDISTEVOLUTION_OUT      string("VOLDISTEvolution")  // moved to aflow_aflowrc.cpp in V3.1.178 
#define _AFLOWLIB_ENTRY_SEPARATOR_       string(" | ")
#define _APENNSY_STYLE_OLD_              FALSE

// --------------------------------------------------------------------------
// definition for frozsl files
#define _FROZSL_VASPSETUP_FILE_ "./aflow.frozsl_vaspsetup_file"
// --------------------------------------------------------------------------
// definitions for WEB PHP
#define AFLOW_PHP_APOOL_REFERENCES       string("19,20,49,50,51,53,54,55,56,57,59,61,62,63,65,66,67,70,71,74,75,76,81,87,99")

// --------------------------------------------------------------------------
// Definitions
//#define DEFAULT_AFLOW_FIND_PARAMETERS "-follow"
#define DEFAULT_AFLOW_FIND_PARAMETERS_NORMAL     string("-follow")
#define DEFAULT_AFLOW_FIND_PARAMETERS_NOLEAF     string("-noleaf -follow")
#define BUFFER_MAXLEN 1024
// [OBSOLETE] #define DEFAULT_AFLOW_PRESCRIPT_OUT            string("aflow.prescript.out")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PRESCRIPT_COMMAND        string("aflow.prescript.command")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_POSTSCRIPT_OUT           string("aflow.postscript.out")  // moved to aflow_aflowrc.cpp in V3.1.189
// [OBSOLETE] #define DEFAULT_AFLOW_POSTSCRIPT_COMMAND       string("aflow.postscript.command")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_OUT               string("aflow.pgroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_XTAL_OUT          string("aflow.pgroup_xtal.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_OUT              string("aflow.pgroupk.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_XTAL_OUT         string("aflow.pgroupk_xtal.out") // moved to aflow_aflowrc.cpp in V3.1.189  // DX 12/5/17 - Added pgroupk_xtal
// [OBSOLETE] #define DEFAULT_AFLOW_FGROUP_OUT               string("aflow.fgroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_SGROUP_OUT               string("aflow.sgroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_AGROUP_OUT               string("aflow.agroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_IATOMS_OUT               string("aflow.iatoms.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_JSON              string("aflow.pgroup.json") // moved to aflow_aflowrc.cpp in V3.1.189       // DX 8/2/17 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_XTAL_JSON         string("aflow.pgroup_xtal.json") // moved to aflow_aflowrc.cpp in V3.1.189  // DX 8/2/17 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_JSON             string("aflow.pgroupk.json")  // moved to aflow_aflowrc.cpp in V3.1.189     // DX 8/2/17 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_XTAL_JSON        string("aflow.pgroupk_xtal.json") // moved to aflow_aflowrc.cpp in V3.1.189 // DX 8/2/17 - Add JSON // DX 12/5/17 - Added pgroupk_xtal
// [OBSOLETE] #define DEFAULT_AFLOW_FGROUP_JSON              string("aflow.fgroup.json") // moved to aflow_aflowrc.cpp in V3.1.189       // DX 8/2/17 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_SGROUP_JSON              string("aflow.sgroup.json")  // moved to aflow_aflowrc.cpp in V3.1.189      // DX 8/2/17 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_AGROUP_JSON              string("aflow.agroup.json")  // moved to aflow_aflowrc.cpp in V3.1.189      // DX 8/2/17 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_IATOMS_JSON              string("aflow.iatoms.json")  // moved to aflow_aflowrc.cpp in V3.1.189      // DX 8/2/17 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PHONON_FILE  string("aflow.phonons.out") // abandoned
// [OBSOLETE] #define DEFAULT_AFLOW_ICAGES_OUT               string("aflow.icages.out")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_SURFACE_OUT              string("aflow.surface.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_QMVASP_OUT               string("aflow.qmvasp.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_ERVASP_OUT               string("aflow.error.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_IMMISCIBILITY_OUT        string("aflow.immiscibility.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_MEMORY_OUT               string("aflow.memory.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_INPUT_OUT         string("aflow.frozsl_input.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_POSCAR_OUT        string("aflow.frozsl_poscar.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_MODES_OUT         string("aflow.frozsl_energies.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_EIGEN_OUT         string("aflow.frozsl_eigen.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_END_OUT                  string("aflow.end.out") // moved to aflow_aflowrc.cpp in V3.1.189 

// --------------------------------------------------------------------------
// include all prototypes for aflow
#define SWAP(a,b)      {temp=(a);(a)=(b);(b)=temp;}
#define RCYCLIC(a,b,c) {temp=(c);(b)=(a);(c)=(b);a=temp;}
#define LCYCLIC(a,b,c) {temp=(a);(a)=(b);(b)=(c);c=temp;}

#define NANOPARTICLE_RADIUS_DEFAULT   10.0
#define NANOPARTICLE_DISTANCE_DEFAULT 10.0
  using aurostd::min;
using aurostd::max;
using aurostd::mod;
using aurostd::_isodd;
using aurostd::_iseven;
using aurostd::_isfloat;
using aurostd::_iscomplex;
using aurostd::sign;
using aurostd::nint;
using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::clear;
using aurostd::xvector;
//[OBSOLETE ME180705]using aurostd::xtensor3;
//[OBSOLETE ME180705]using aurostd::xtensor4;
//[OBSOLETE ME180705]using aurostd::xtensor5;
//[OBSOLETE ME180705]using aurostd::xtensor6;
//[OBSOLETE ME180705]using aurostd::xtensor7;
//[OBSOLETE ME180705]using aurostd::xtensor8;
using aurostd::xoption;
using aurostd::xcombos;

// --------------------------------------------------------------------------
// this is a container of general global choices
class _XHOST {
 public:
  // constructor destructor                         // constructor/destructor
  _XHOST();                                         // default, just allocate
  ~_XHOST();                                        // kill everything
  // _XHOST(const _XHOST& b);                          // constructor copy
  const _XHOST& operator=(const _XHOST &b);         // copy
  // BOOT
  int PID;  // aflow.cpp
  ostringstream ostrPID; // aflow.cpp
  // machinery
  bool QUIET,TEST,DEBUG,MPI;
  bool GENERATE_AFLOWIN_ONLY; //CT 180719
  // HARDWARE/SOFTWARE
  string hostname,MachineType,Tmpfs,User,Group,Home,Shell,Progname;
  string Find_Parameters;
  bool sensors_allowed;
  // ARGUMENTS
  vector<string> argv;          // argv of line command
  // SERVERS
  string AFLOW_MATERIALS_SERVER,AFLOW_WEB_SERVER;
  long double RAM,RAM_MB,RAM_GB;
  int CPU_Cores;
  string CPU_Model;
  string CPU_MHz;
  vector<double> vTemperatureCore;
  long double Time_starting,Time_now;
  long int Date;
  string Day,Month,Year;
  string Copyright_Years; // =string("2003-YEAR_FROM_DATE");
  // MULTHREADS
  /* bool PTHREADS_FLAG;        // run pthread YES/NO */
  /* int  PTHREADS_MAX;         // how many MAX threads I can use  default or --np */
  /* int PTHREADS_RUNNING;      // how many threads are actually running */
  /* vector<pthread_t> thread;  // the actual thread */
  /* vector<int> iret;          // the thread runnings */
  /* vector<bool> thread_busy;  // is the thread busy */
  // COMMANDS
  vector<string> vcmd;  
  // RAM CHECK
  double maxmem;
  // FUNCTIONS
  string command(const string& command);
  bool is_command(const string& command);
  // AFLOW STUFF
  // vflag_aflow.flag("LOOP");
  // vflag_aflow.flag("CLEAN");
  // vflag_aflow.isflag*"XCLEAN");
  bool AFLOW_RUNDIRflag;
  bool AFLOW_MULTIflag;
  bool AFLOW_RUNXflag;
  uint AFLOW_RUNXnumber;
  // QUEQUE STUFF
  bool is_PBS;int PBS_NUM_PPN,PBS_NNODES;
  bool is_SLURM;int SLURM_CPUS_ON_NODE,SLURM_NNODES;
  bool is_MACHINE_FULTON_MARYLOU; // some flags
  // APENNST stuff
  bool APENNSY_USE_SERVER;
  bool APENNSY_USE_LIBRARY;
  bool APENNSY_SERVER_AFLOWLIB_ORG;
  // Library_CALCULATED*
  vector<uint>   vGlobal_uint;      // parameters uint
  vector<string> vGlobal_string;    // parameters as strings
  vector<vector<string> > vvGlobal_string; // parameters as vector strings
  // vector<string> vLibrary_ICSD;     // ordered by #species (needs to be allocated)
  // vector<string> vLibrary_ICSD_ALL; // line by line
  // string Library_ICSD_ALL;          // the complete library
  // vector<string> vVASP_POTCAR_DIRECTORIES;
  // vector<string> vAFLOW_LIBRARY_DIRECTORIES;
  // vector<string> vAFLOW_PROJECTS_DIRECTORIES;
  // AFLOW flags/options
  aurostd::xoption vflag_aflow;  // argv/argc options following the xoption structure
  aurostd::xoption vflag_pflow;  // argv/argc options following the xoption structure
  aurostd::xoption vflag_apennsy;  // argv/argc options following the xoption structure
  aurostd::xoption vflag_outreach;  // argv/argc options following the xoption structure
  aurostd::xoption vflag_control;  // argv/argc options following the xoption structure
  // AFLOWRC
  string aflowrc_filename;
  string aflowrc_content;
  vector<string> vaflowrc;  
  xoption adefault;            // default  xoption
  // AFLOWSYM
  bool SKEW_TEST; // DX 10/19/17
  double SKEW_TOL; // DX 10/19/17
 private:                                                //
  void free();                                           // free space
  void copy(const _XHOST& b);                            //
  void clear();                                          // free space
};

#define XHOST_vGlobal_MAX                              256
#define XHOST_Library_HTQC                             XHOST.vGlobal_string.at(0)
#define XHOST_Library_CALCULATED_ICSD_LIB              XHOST.vGlobal_string.at(1)
#define XHOST_Library_CALCULATED_ICSD_RAW              XHOST.vGlobal_string.at(2)
#define XHOST_Library_CALCULATED_LIB1_LIB              XHOST.vGlobal_string.at(3)
#define XHOST_Library_CALCULATED_LIB1_RAW              XHOST.vGlobal_string.at(4)
#define XHOST_Library_CALCULATED_LIB2_LIB              XHOST.vGlobal_string.at(5)
#define XHOST_Library_CALCULATED_LIB2_RAW              XHOST.vGlobal_string.at(6)
#define XHOST_Library_CALCULATED_LIB3_LIB              XHOST.vGlobal_string.at(7)
#define XHOST_Library_CALCULATED_LIB3_RAW              XHOST.vGlobal_string.at(8)
#define XHOST_Library_CALCULATED_LIB4_LIB              XHOST.vGlobal_string.at(9)
#define XHOST_Library_CALCULATED_LIB4_RAW              XHOST.vGlobal_string.at(10)
#define XHOST_Library_CALCULATED_LIB5_LIB              XHOST.vGlobal_string.at(11)
#define XHOST_Library_CALCULATED_LIB5_RAW              XHOST.vGlobal_string.at(12)
#define XHOST_Library_CALCULATED_LIB6_LIB              XHOST.vGlobal_string.at(13)
#define XHOST_Library_CALCULATED_LIB6_RAW              XHOST.vGlobal_string.at(14)
#define XHOST_Library_CALCULATED_LIB7_LIB              XHOST.vGlobal_string.at(15)
#define XHOST_Library_CALCULATED_LIB7_RAW              XHOST.vGlobal_string.at(16)
#define XHOST_Library_CALCULATED_LIB8_LIB              XHOST.vGlobal_string.at(17)
#define XHOST_Library_CALCULATED_LIB8_RAW              XHOST.vGlobal_string.at(18)
#define XHOST_Library_CALCULATED_LIB9_LIB              XHOST.vGlobal_string.at(19)
#define XHOST_Library_CALCULATED_LIB9_RAW              XHOST.vGlobal_string.at(20)
#define XHOST_aflowlib_icsd                            XHOST.vGlobal_string.at(21)
#define XHOST_aflowlib_lib1                            XHOST.vGlobal_string.at(22)
#define XHOST_aflowlib_lib2                            XHOST.vGlobal_string.at(23)
#define XHOST_aflowlib_lib3                            XHOST.vGlobal_string.at(24)
#define XHOST_aflowlib_lib4                            XHOST.vGlobal_string.at(25)
#define XHOST_aflowlib_lib5                            XHOST.vGlobal_string.at(26)
#define XHOST_aflowlib_lib6                            XHOST.vGlobal_string.at(27)
#define XHOST_aflowlib_lib7                            XHOST.vGlobal_string.at(28)
#define XHOST_aflowlib_lib8                            XHOST.vGlobal_string.at(29)
#define XHOST_aflowlib_lib9                            XHOST.vGlobal_string.at(30)
#define XHOST_AUID                                     XHOST.vGlobal_string.at(31)
#define XHOST_AURL                                     XHOST.vGlobal_string.at(32)
#define XHOST_LOOP                                     XHOST.vGlobal_string.at(33)
#define XHOST_Library_ICSD_ALL                         XHOST.vGlobal_string.at(34)
//  string Library_ICSD_ALL;          // the complete library

#define XHOST_vAUID                                    XHOST.vvGlobal_string.at(0)
#define XHOST_vAURL                                    XHOST.vvGlobal_string.at(1)
#define XHOST_vLOOP                                    XHOST.vvGlobal_string.at(2)
#define vAUID XHOST_vAUID
#define vAURL XHOST_vAURL
#define vLOOP XHOST_vLOOP

#define vVASP_POTCAR_DIRECTORIES                       XHOST.vvGlobal_string.at(3)
#define vAFLOW_LIBRARY_DIRECTORIES                     XHOST.vvGlobal_string.at(4)
#define vAFLOW_PROJECTS_DIRECTORIES                    XHOST.vvGlobal_string.at(5)
#define XHOST_vLibrary_ICSD                            XHOST.vvGlobal_string.at(6)
#define XHOST_vLibrary_ICSD_ALL                        XHOST.vvGlobal_string.at(7)
//  vector<string> vLibrary_ICSD;     // ordered by #species
//  vector<string> vLibrary_ICSD_ALL; // line by line

#define XHOST_README_AFLOW_LICENSE_GPL3_TXT            XHOST.vGlobal_string.at(40)
#define XHOST_README_AFLOW_TXT                         XHOST.vGlobal_string.at(41)
#define XHOST_README_AFLOW_VERSIONS_HISTORY_TXT        XHOST.vGlobal_string.at(42)
#define XHOST_README_AFLOW_PFLOW_TXT                   XHOST.vGlobal_string.at(43)
#define XHOST_README_AFLOW_APENNSY_TXT                 XHOST.vGlobal_string.at(44)
#define XHOST_README_AFLOW_SCRIPTING_TXT               XHOST.vGlobal_string.at(45)
#define XHOST_README_AFLOW_FROZSL_TXT                  XHOST.vGlobal_string.at(46)
#define XHOST_README_AFLOW_POCC_TXT                    XHOST.vGlobal_string.at(47)
#define XHOST_README_AFLOW_APL_TXT                     XHOST.vGlobal_string.at(48)
#define XHOST_README_AFLOW_QHA_SCQHA_QHA3P_TXT         XHOST.vGlobal_string.at(49)
#define XHOST_README_AFLOW_AGL_TXT                     XHOST.vGlobal_string.at(50)
#define XHOST_README_AFLOW_AEL_TXT                     XHOST.vGlobal_string.at(51)
#define XHOST_README_AFLOW_ANRL_TXT                    XHOST.vGlobal_string.at(52)
#define XHOST_README_AFLOW_COMPARE_TXT                 XHOST.vGlobal_string.at(53)
#define XHOST_README_AFLOW_SYM_TXT                     XHOST.vGlobal_string.at(54)
#define XHOST_README_AFLOW_CHULL_TXT                   XHOST.vGlobal_string.at(55)
#define XHOST_README_AFLOW_EXCEPTIONS_TXT              XHOST.vGlobal_string.at(56) //ME180705
#define XHOST_README_AFLOW_HTRESOURCES_TXT             XHOST.vGlobal_string.at(57)
#define XHOST_README_PROTO_TXT                         XHOST.vGlobal_string.at(58)
#define XHOST_README_AFLOW_XAFLOW_TXT                  XHOST.vGlobal_string.at(59)
#define XHOST_README_AFLOW_AFLOWRC_TXT                 XHOST.vGlobal_string.at(60)

#define XHOST_FINDSYM_data_space_txt                   XHOST.vGlobal_string.at(70)
#define XHOST_FINDSYM_data_wyckoff_txt                 XHOST.vGlobal_string.at(71)
#define XHOST_FROZSL_data_space_txt                    XHOST.vGlobal_string.at(72)
#define XHOST_FROZSL_data_wyckoff_txt                  XHOST.vGlobal_string.at(73)
#define XHOST_FROZSL_data_images_txt                   XHOST.vGlobal_string.at(74)
#define XHOST_FROZSL_data_irreps_txt                   XHOST.vGlobal_string.at(75)
#define XHOST_FROZSL_data_isotropy_txt                 XHOST.vGlobal_string.at(76)
#define XHOST_FROZSL_data_little_txt                   XHOST.vGlobal_string.at(77)
#define XHOST_FROZSL_symmetry2_dat                     XHOST.vGlobal_string.at(78)
#define XHOST_FROZSL_const_dat                         XHOST.vGlobal_string.at(79)
#define XHOST_FROZSL_phvaspsetup_AFLOW                 XHOST.vGlobal_string.at(80)
#define XHOST_FROZSL_phvaspsetup_POSCAR                XHOST.vGlobal_string.at(81)
#define XHOST_ElectronStoppingPower_txt                XHOST.vGlobal_string.at(82)
#define XHOST_PhotonCrossSection_txt                   XHOST.vGlobal_string.at(83)
#define XHOST_PhotonStoppingPower_txt                  XHOST.vGlobal_string.at(84)
#define XHOST_ICSD_List_txt                            XHOST.vGlobal_string.at(85)
#define XHOST_AFLOW_PSEUDOPOTENTIALS                   XHOST.vGlobal_string.at(86)
#define XHOST_AFLOW_PSEUDOPOTENTIALS_TXT               XHOST.vGlobal_string.at(87)
#define XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT          XHOST.vGlobal_string.at(88)
#define XHOST_f144468a7ccc2d3a72ba44000715efdb         XHOST.vGlobal_string.at(90)
#define XHOST_d0f1b0e47f178ae627a388d3bf65d2d2         XHOST.vGlobal_string.at(91)
#define XHOST_decf00ca3ad2fe494eea8e543e929068         XHOST.vGlobal_string.at(92)
// [OBSOLETE] #define XHOST_AFLOW_BinaryRead           XHOST.vGlobal_string.at(93)
// [OBSOLETE] #define XHOST_AFLOW_Binary_Angle_Read    XHOST.vGlobal_string.at(94)

  
// LOADENTRIES DEFAULTS
#define _AFLOW_LIB_MAX_ 10                             //LIB11 does not exist yet, modify accordingly
  
#define XHOST_LIBRARY_LIB1                             XHOST.vGlobal_uint.at(0)
#define XHOST_LIBRARY_LIB2                             XHOST.vGlobal_uint.at(1)
#define XHOST_LIBRARY_LIB3                             XHOST.vGlobal_uint.at(2)
#define XHOST_LIBRARY_LIB4                             XHOST.vGlobal_uint.at(3)
#define XHOST_LIBRARY_LIB5                             XHOST.vGlobal_uint.at(4)
#define XHOST_LIBRARY_LIB6                             XHOST.vGlobal_uint.at(5)
#define XHOST_LIBRARY_LIB7                             XHOST.vGlobal_uint.at(6)
#define XHOST_LIBRARY_LIB8                             XHOST.vGlobal_uint.at(7)
#define XHOST_LIBRARY_LIB9                             XHOST.vGlobal_uint.at(8)
#define XHOST_LIBRARY_ICSD                             XHOST.vGlobal_uint.at(9)

// max is 128
extern _XHOST XHOST; // this will be global

// ----------------------------------------------------------------------------
// aflow_aflowrc.cpp
#define _AFLOW_AFLOWRC_H_
#define _AFLOW_AFLOWRC_CPP_
#include "aflow_aflowrc.cpp"
#undef _AFLOW_AFLOWRC_CPP_
#undef _AFLOW_AFLOWRC_H_

// --------------------------------------------------------------------------

// Structures for flags and properties to share FAST !
// STRUCTURES
#define AFLOWIN_SEPARATION_LINE  string("[AFLOW] ************************************************************************************************************************** ")

// --------------------------------------------------------------------------
// general flags to run aflow
class _aflags {
 public:
  // trivial constructurs/destuctors/operators
  _aflags();                                          // default, just allocate
  ~_aflags();                                         // kill everything
  _aflags(const _aflags& b);                          // constructor copy
  const _aflags& operator=(const _aflags &b);         // copy
  void clear(void);                                   // clear
  // CONTENT
  bool QUIET;
  int  AFLOW_PTHREADS_NUMBER;                         // cant be GLOBAL as this is a local run stuff
  // particular
  string LocalDirectory;                              // where is aflow now
  string Directory;                                   // where aflow must run
  bool AFLOW_FORCE_RUN;                               // Force run also in database
  bool AFLOW_PERFORM_DIRECTORY;                       // Directory is specified (sometimes it is useful).
  bool AFLOW_PERFORM_FILE;                            // File is specified (sometimes it is useful).
  bool AFLOW_PERFORM_ORDER_SORT;                      // Sorts the _AFLOWIN_ in the list
  bool AFLOW_PERFORM_ORDER_REVERSE;                   // Reverse the _AFLOWIN_ in the list
  bool AFLOW_PERFORM_ORDER_RANDOM;                    // Randomize the _AFLOWIN_ in the list
  bool AFLOW_MODE_GENERATE;                           // TODO OVVERRIDE all _AFLOWIN_
  bool AFLOW_MODE_QSUB_MODE1;                         // TODO OVVERRIDE all _AFLOWIN_
  bool AFLOW_MODE_QSUB_MODE2;                         // TODO OVVERRIDE all _AFLOWIN_
  bool AFLOW_MODE_QSUB_MODE3;                         // TODO OVVERRIDE all _AFLOWIN_
  // general flags to operate in the directory
  bool KBIN_RUN_AFLOWIN;
  bool KBIN_GEN_GENERAL; // CO 180409
  bool KBIN_GEN_VASP_FROM_AFLOWIN;
  bool KBIN_GEN_AIMS_FROM_AFLOWIN;
  bool KBIN_GEN_AFLOWIN_FROM_VASP;
  // DX - START
  bool KBIN_GEN_SYMMETRY_OF_AFLOWIN;
  // DX - END 
  bool KBIN_DELETE_AFLOWIN;
  bool AFLOW_FORCE_MPI;     // not yet implemented
  bool AFLOW_FORCE_SERIAL;  // not yet implemented
  int  AFLOW_GLOBAL_NCPUS;         // Forced CPUS
  // Perform TASKS
  bool AFLOW_PERFORM_CLEAN;       // to clean a directory
  // host related things
  xoption AFLOW_MACHINE_GLOBAL;
  xoption AFLOW_MACHINE_LOCAL;       // flag for duke_beta_mpich
  // APENNSY
  aurostd::xoption vflag;
  string APENNSY_LATTICE_flag;                        // APENNSY flags
  string APENNSY_GNUPLOT_FONT_str;                    // APENNSY content
  string APENNSY_GNUPLOT_FONT_BOLD_str;               // APENNSY content
  string APENNSY_GNUPLOT_FONT_ITALICS_str;            // APENNSY content
  // [OBSOLETE] "APENNSY::HELP;                                   // APENNSY flags
  // [OBSOLETE] "APENNSY::VERBOSE_flag;                           // APENNSY flags
  // [OBSOLETE] "APENNSY_LIST_flag;                               // APENNSY flags
  // [OBSOLETE] "APENNSY::SERVER_AFLOWLIB_ORG_flag;               // APENNSY flags
  // [OBSOLETE] "APENNSY::LATEX_SNAPSHOT;                         // APENNSY flags
  // [OBSOLETE] "APENNSY::LATEX_OUTPUT;                           // APENNSY flags
  // [OBSOLETE] "APENNSY::LATEX_CITE;                             // APENNSY flags
  // [OBSOLETE] "APENNSY::ENTHALPY_TOT;                           // APENNSY flags
  // [OBSOLETE] "APENNSY::ENTHALPY_ATOM;                          // APENNSY flags
  // [OBSOLETE] "APENNSY::ENTHALPY_FORMATION_ATOM;                // APENNSY flags
  // [OBSOLETE] "APENNSY::LOAD_LIB2;                              // APENNSY flags
  // [OBSOLETE] "APENNSY::LOAD_LIB2U;                             // APENNSY flags
  // [OBSOLETE] "APENNSY::LOAD_LIB2PGM;                           // APENNSY flags
  // [OBSOLETE] "APENNSY::LOAD_LIB2X;                             // APENNSY flags
  // [OBSOLETE] "APENNSY::LOAD_ALLOY;                             // APENNSY flags
  // [OBSOLETE] "APENNSY::APOOL_PUBLIC;                           // APENNSY flags
  // [OBSOLETE] "APENNSY::APOOL_PRIVATE;                          // APENNSY flags
  // [OBSOLETE] "APENNSY::APOOL_TEST;                             // APENNSY flags
  // [OBSOLETE] "APENNSY::DATA;                                   // APENNSY flags
  // [OBSOLETE] "APENNSY::UNCLE;                                  // APENNSY flags
  // [OBSOLETE] "APENNSY::WEB;                                    // APENNSY flags
  // [OBSOLETE] "APENNSY::ALL;                                    // APENNSY flags
  // [OBSOLETE] "APENNSY::FCC;                                    // APENNSY flags
  // [OBSOLETE] "APENNSY::BCC;                                    // APENNSY flags
  // [OBSOLETE] "APENNSY::HCP;                                    // APENNSY flags
  // [OBSOLETE] "APENNSY::COUT;                                   // APENNSY flags
  // [OBSOLETE] "APENNSY::CERR;                                   // APENNSY flags
  // [OBSOLETE] "APENNSY::ENTHALPY_LIST"                          // APENNSY flags
  // [OBSOLETE] "APENNSY::PS_ENERGY_LIST"                         // APENNSY flags
  // [OBSOLETE] "APENNSY::CONVEX_HULL"                            // APENNSY flags
  // [OBSOLETE] "APENNSY::MATLAB"                                 // APENNSY flags
  // [OBSOLETE] "APENNSY::GNUPLOT"                                // APENNSY flags
  // [OBSOLETE] "APENNSY::SMALL_CONVEX_HULL_MATLAB"               // APENNSY flags
  // [OBSOLETE] "APENNSY::HISTOGRAM_LIST"                         // APENNSY flags
  // [OBSOLETE] "APENNSY::MATLAB_LIB"                             // APENNSY flags
  // [OBSOLETE] "APENNSY::RULES"                                  // APENNSY flags
  // [OBSOLETE] "APENNSY::STRUCTURES"                             // APENNSY flags
  // [OBSOLETE] "APENNSY::VASPIN;                                 // APENNSY flags
  // [OBSOLETE] "APENNSY::ORDER;                                  // APENNSY flags
  // [OBSOLETE] "APENNSY::INFO;                                   // APENNSY flags
  // [OBSOLETE] "APENNSY::MISCIBILITY;                            // APENNSY flags
  // [OBSOLETE] "APENNSY::MISCIBILITY_EXPERIMENTS;                // APENNSY flags
  // [OBSOLETE] "APENNSY::MISCIBILITY_MIEDEMA;                    // APENNSY flags
  // [OBSOLETE] "APENNSY::MISCIBILITY_HUMEROTHERY;                // APENNSY flags
  // [OBSOLETE] "APENNSY::MISCIBILITY_TABLE;                      // APENNSY flags
  // [OBSOLETE] "APENNSY::MISCIBILITY_STATISTICS;                 // APENNSY flags
  // [OBSOLETE] "APENNSY::STRUCTURE_VOLUMES;                      // APENNSY flags
  // [OBSOLETE] "APENNSY::REFERENCE;                              // APENNSY flags
  // [OBSOLETE] "APENNSY_PROTOCHECK;                              // APENNSY flags
  // [OBSOLETE] "APENNSY::NEGLECT_STRUCTURES;                     // APENNSY flags
  // [OBSOLETE] "APENNSY::UPDATE;                                 // APENNSY flags
  // [OBSOLETE] "APENNSY::CSWAP;                                  // APENNSY flags
  vector<string> APENNSY_NEGLECT_STRUCTURES_vstrs;    // APENNSY content
 private:                                              //
  void free();                                        // free space
  void copy(const _aflags& b);                        //
};

// --------------------------------------------------------------------------
// general flags for kbinary (all)
class _kflags {
 public:
  // trivial constructurs/destuctors/operators
  _kflags();                                          // default, just allocate
  ~_kflags();                                         // kill everything
  _kflags(const _kflags& b);                          // constructor copy
  const _kflags& operator=(const _kflags &b);         // copy
  void clear(void);                                   // clear
  // CONTENT
  // in this struct we put all the flags which are used on LOCAL DIRECTORIES in KBIN MODE
  //
  bool AFLOW_MODE_ALIEN;
  //
  bool AFLOW_MODE_MATLAB;
  bool AFLOW_MATLAB_MODE_EXPLICIT;
  bool AFLOW_MATLAB_MODE_EXPLICIT_START_STOP;
  bool AFLOW_MATLAB_MODE_IMPLICIT;
  bool AFLOW_MATLAB_MODE_EXTERNAL;
  bool AFLOW_MATLAB_FILE;
  bool AFLOW_MATLAB_FILE_FILE;
  bool AFLOW_MATLAB_FILE_COMMAND;
  //
  bool AFLOW_MODE_VASP;
  bool AFLOW_MODE_AIMS;
  //
  bool AFLOW_MODE_PRESCRIPT_EXPLICIT;
  bool AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP;
  stringstream AFLOW_MODE_PRESCRIPT;
  bool AFLOW_MODE_POSTSCRIPT_EXPLICIT;
  bool AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP;
  stringstream AFLOW_MODE_POSTSCRIPT;
  //
  bool AFLOW_MODE_EMAIL;
  // normal binary
  string KBIN_BIN;
  string KZIP_BIN;
  bool   KZIP_COMPRESS;
  // MPI binaries and flags
  bool   KBIN_MPI;
  int    KBIN_MPI_NCPUS;
  int    KBIN_MPI_NCPUS_BUFFER;
  string KBIN_MPI_START;
  string KBIN_MPI_STOP;
  string KBIN_MPI_COMMAND;
  bool   KBIN_MPI_AUTOTUNE;
  string KBIN_MPI_BIN;
  string KBIN_MPI_OPTIONS;
  // QSUB
  bool   KBIN_QSUB;
  bool   KBIN_QSUB_MODE1;
  bool   KBIN_QSUB_MODE2;
  bool   KBIN_QSUB_MODE3;
  string KBIN_QSUB_COMMAND;
  string KBIN_QSUB_PARAMS;
  bool   KBIN_QSUB_MODE_EXPLICIT;
  bool   KBIN_QSUB_MODE_EXPLICIT_START_STOP;
  bool   KBIN_QSUB_MODE_IMPLICIT;
  bool   KBIN_QSUB_FILE;
  // symmetry operation lists
  bool  KBIN_SYMMETRY_CALCULATION;
  // DX - START
  bool  KBIN_SYMMETRY_NO_SCAN;
  double KBIN_SYMMETRY_EPS;
  bool  KBIN_SYMMETRY_CALCULATE_PGROUP;       // DX 8/14/17 - Specify what to calculate/verify
  bool  KBIN_SYMMETRY_CALCULATE_PGROUPK;      // DX 8/14/17 - Specify what to calculate/verify
  bool  KBIN_SYMMETRY_CALCULATE_FGROUP;       // DX 8/14/17 - Specify what to calculate/verify
  bool  KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL;  // DX 8/14/17 - Specify what to calculate/verify
  bool  KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL; // DX 12/5/17 - Specify what to calculate/verify; Added pgroupk_xtal
  bool  KBIN_SYMMETRY_CALCULATE_IATOMS;       // DX 8/14/17 - Specify what to calculate/verify
  bool  KBIN_SYMMETRY_CALCULATE_AGROUP;       // DX 8/14/17 - Specify what to calculate/verify
  bool  KBIN_SYMMETRY_CALCULATE_SGROUP;       // DX 8/14/17 - Specify what to calculate/verify
  // DX - END
  bool  KBIN_SYMMETRY_PGROUP_WRITE;      // taken TRUE by default
  bool  KBIN_SYMMETRY_PGROUPK_WRITE;     // taken TRUE by default
  bool  KBIN_SYMMETRY_PGROUP_XTAL_WRITE; // taken TRUE by default
  bool  KBIN_SYMMETRY_PGROUPK_XTAL_WRITE;// DX 12/5/17 - Added pgroupk_xtal
  bool  KBIN_SYMMETRY_FGROUP_WRITE;      // taken TRUE by default
  bool  KBIN_SYMMETRY_SGROUP_WRITE;
  bool  KBIN_SYMMETRY_AGROUP_WRITE;      // taken TRUE by default
  bool  KBIN_SYMMETRY_IATOMS_WRITE;      // taken TRUE by default
  double KBIN_SYMMETRY_SGROUP_RADIUS;
  // neighbours operation lists
  bool  KBIN_NEIGHBOURS_CALCULATION;
  bool  KBIN_NEIGHBOURS_WRITE;
  double KBIN_NEIGHBOURS_RADIUS;
  double KBIN_NEIGHBOURS_DRADIUS;
  // pocc operation lists
  bool   KBIN_POCC;
  bool   KBIN_POCC_CALCULATION;
  // frozsl operation lists
  bool   KBIN_FROZSL;
  bool   KBIN_FROZSL_DOWNLOAD;
  bool   KBIN_FROZSL_FILE;
  string KBIN_FROZSL_FILE_NAME;
  // [OBSOLETE]  bool   KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT;
  // [OBSOLETE]  bool   KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT_START_STOP;
  // [OBSOLETE]  string KBIN_FROZSL_PRESCRIPT_STRING;
  // [OBSOLETE]  bool   KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT;
  // [OBSOLETE]  bool   KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT_START_STOP;
  // [OBSOLETE]  string KBIN_FROZSL_POSTSCRIPT_STRING;
  bool   KBIN_FROZSL_STRUCTURE_MODE_FILE;
  bool   KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP;
  string KBIN_FROZSL_STRUCTURE_STRING;
  bool   KBIN_FROZSL_DIELECTRIC_MODE_FILE;
  bool   KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP;
  bool   KBIN_FROZSL_DIELECTRIC_ZEFF;
  string KBIN_FROZSL_DIELECTRIC_STRING;
  // phonons operation lists
  bool   KBIN_PHONONS_CALCULATION_APL;
  bool   KBIN_PHONONS_CALCULATION_QHA;  // CO - 170601
  bool   KBIN_PHONONS_CALCULATION_QHA_A;    //PN180705
  bool   KBIN_PHONONS_CALCULATION_QHA_B;    //PN180705
  bool   KBIN_PHONONS_CALCULATION_QHA_C;    //PN180705
  bool   KBIN_PHONONS_CALCULATION_QHA3P;    //PN180705
  bool   KBIN_PHONONS_CALCULATION_QHA3P_A;  //PN180705
  bool   KBIN_PHONONS_CALCULATION_QHA3P_B;  //PN180705
  bool   KBIN_PHONONS_CALCULATION_QHA3P_C;  //PN180705
  bool   KBIN_PHONONS_CALCULATION_SCQHA;    //PN180705
  bool   KBIN_PHONONS_CALCULATION_SCQHA_A;  //PN180705
  bool   KBIN_PHONONS_CALCULATION_SCQHA_B;  //PN180705
  bool   KBIN_PHONONS_CALCULATION_SCQHA_C;  //PN180705
  bool   KBIN_PHONONS_CALCULATION_AAPL; // CO - 170601
  bool   KBIN_PHONONS_CALCULATION_AGL;
  bool   KBIN_PHONONS_CALCULATION_AEL;
  bool   KBIN_PHONONS_CALCULATION_FROZSL;
  string KBIN_PHONONS_CALCULATION_FROZSL_output;
  string KBIN_PHONONS_CALCULATION_FROZSL_poscars;
 private:                                             //
  void free();                                        // free space
  void copy(const _kflags& b);                        //
};

// --------------------------------------------------------------------------
// general flags for vasp mode
class xstructure; // prototype of structure, just to compile
class _vflags {
 public:
  // trivial constructurs/destuctors/operators
  _vflags();                                            // default, just allocate
  ~_vflags();                                           // kill everything
  _vflags(const _vflags& b);                            // constructor copy
  const _vflags& operator=(const _vflags &b);           // copy
  void clear(void);                                     // clear
  // CONTENT
  // in this struct we put all the flags which are used on LOCAL DIRECTORIES in VASP MODE
  int KBIN_VASP_RUN_NRELAX;
  xoption KBIN_VASP_RUN;                        // GENERATE, STATIC, KPOINTS, RELAX, RELAX_STATIC, RELAX_STATIC_BANDS, STATIC_BANDS, DIELECTRIC_STATIC, DIELECTRIC_DYNAMIC, DSCF
  xoption KBIN_VASP_REPEAT;                     // REPEAT_BANDS REPEAT_STATIC_BANDS REPEAT_DELSOL
  xoption KBIN_VASP_FORCE_OPTION_NOTUNE;        // NOTUNE
  xoption KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO;   // SYSTEM_AUTO
  xoption KBIN_VASP_FORCE_OPTION_RELAX_MODE;    // RELAX_MODE  forces/energy
  xoption KBIN_VASP_FORCE_OPTION_RELAX_TYPE;    // RELAX_TYPE  STATIC, ALL, IONS, CELL_SHAPE, CELL_VOLUME, IONS_CELL_VOLUME
  xoption KBIN_VASP_FORCE_OPTION_PREC;          // PREC 
  xoption KBIN_VASP_FORCE_OPTION_ALGO;          // ALGO 
  xoption KBIN_VASP_FORCE_OPTION_METAGGA;       // METAGGA 
  xoption KBIN_VASP_FORCE_OPTION_IVDW;          // IVDW
  xoption KBIN_VASP_FORCE_OPTION_ABMIX;         // ABMIX
  xoption KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS; // AUTO_PSEUDOPOTENTIALS
  // ENMAX_MULTIPLY
  xoption KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL;  // isentry and content_int
  // NBANDS
  xoption KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL;  // isentry and content_int
  bool KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry;
  // POTIM
  xoption KBIN_VASP_FORCE_OPTION_POTIM_EQUAL;   // isentry and content_double
  // PSTRESS
  xoption KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL; // isentry and content_double
  // EDIFFG
  xoption KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL; // isentry and content_double
  // RWIGS
  bool KBIN_VASP_FORCE_OPTION_RWIGS_STATIC;  
  xoption KBIN_VASP_FORCE_OPTION_SKIP_NOMIX;    // SKIP_NOMIX
  xoption KBIN_VASP_FORCE_OPTION_SPIN;          // SPIN 
  bool KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1;
  bool KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2;
  // xoption KBIN_VASP_FORCE_OPTION_TRISTATE;      //  SYM 
  xoption KBIN_VASP_FORCE_OPTION_BADER;         // BADER=ON | OFF | NONE
  xoption KBIN_VASP_FORCE_OPTION_ELF;           // ELF=ON | OFF | NONE
  xoption KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM;   // AUTO_MAGMOM
  xoption KBIN_VASP_FORCE_OPTION_SYM;           // SYM
  xoption KBIN_VASP_FORCE_OPTION_WAVECAR;       // WAVECAR
  xoption KBIN_VASP_FORCE_OPTION_CHGCAR;        // CHGCAR
  xoption KBIN_VASP_FORCE_OPTION_LSCOUPLING;    // LSCOUPLING
  xoption KBIN_VASP_FORCE_OPTION_LDAU0;         // LDAU0
  xoption KBIN_VASP_FORCE_OPTION_LDAU1;         // LDAU1
  xoption KBIN_VASP_FORCE_OPTION_LDAU2;         // LDAU2
  xoption KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC;// LDAU_ADIABATIC
  xoption KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF;   // LDAU_CUTOFF
  string KBIN_VASP_LDAU_SPECIES;
  string KBIN_VASP_LDAU_PARAMETERS;
  bool KBIN_VASP_LDAU_AFLOW_AUTO_flag;
  // FORCE_OPTION
  xoption KBIN_VASP_FORCE_OPTION_TYPE;          // TYPE 
  bool KBIN_VASP_FORCE_OPTION_NSW_EQUAL;
  int  KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE;

  xoption KBIN_VASP_FORCE_OPTION_IGNORE_AFIX;   // AFIX
  // xoption kopts;
  xoption KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL;   // CONVERT_UNIT_CELL
  xoption KBIN_VASP_FORCE_OPTION_VOLUME;        // EQUAL_EQUAL, MULTIPLY_EQUAL,PLUS_EQUAL
  xoption KBIN_VASP_FORCE_OPTION_KPOINTS;       // KPOINTS 
  xoption KBIN_VASP_INCAR_MODE;                 // EXPLICIT, EXPLICIT_START_STOP, IMPLICIT, EXTERNAL;
  // RELAX
  xoption KBIN_VASP_KPOINTS_MODE;               // EXPLICIT, EXPLICIT_START_STOP, IMPLICIT, EXTERNAL;
  xoption KBIN_VASP_KPOINTS_KMODE;              // isentry and content_int
  xoption KBIN_VASP_KPOINTS_KPPRA;              // isentry and content_int
  xoption KBIN_VASP_KPOINTS_KSCHEME;            // isentry and content_string
  xoption KBIN_VASP_KPOINTS_KSHIFT;             // isentry and content_string
  // STATIC
  xoption KBIN_VASP_KPOINTS_STATIC_KMODE;       // isentry and content_int
  xoption KBIN_VASP_KPOINTS_STATIC_KPPRA;       // isentry and content_int
  xoption KBIN_VASP_KPOINTS_STATIC_KSCHEME;     // isentry and content_string
  xoption KBIN_VASP_KPOINTS_STATIC_KSHIFT;      // isentry and content_string
  // PHONONS
  xoption KBIN_VASP_KPOINTS_PHONONS_KPPRA;      // isentry and content_int
  xoption KBIN_VASP_KPOINTS_PHONONS_KSCHEME;    // isentry and content_string
  xoption KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY;  // EVEN ODD
  // BANDS
  xoption KBIN_VASP_KPOINTS_BANDS_LATTICE;
  //  bool KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG;
  //  string KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE;
  bool KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG;
  bool KBIN_VASP_KPOINTS_BANDS_GRID_FLAG;
  uint KBIN_VASP_KPOINTS_BANDS_GRID_VALUE;
  bool KBIN_VASP_WRITE_KPOINTS;

  xoption KBIN_VASP_POSCAR_MODE; // EXPLICIT, EXPLICIT_START_STOP, EXPLICIT_START_STOP_POINT, IMPLICIT, EXTERNAL;
  std::vector<string> KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING;
  std::vector<xstructure> KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE;
  xoption KBIN_VASP_POTCAR_MODE;  // EXPLICIT, IMPLICIT, EXTERNAL;
  bool KBIN_VASP_INCAR_VERBOSE;   // VERBOSITY
  xoption KBIN_VASP_INCAR_FILE;   // KEYWORD, SYSTEM_AUTO, FILE, COMMAND
  xoption KBIN_VASP_KPOINTS_FILE; // KEYWORD, FILE, COMMAND
  xoption KBIN_VASP_POSCAR_FILE;  // KEYWORD, PROTOTYPE, FILE, COMMAND
  xoption KBIN_VASP_POSCAR_FILE_VOLUME; // EQUAL_EQUAL, MULTIPLY_EQUAL PLUS_EQUAL
  xoption KBIN_VASP_POTCAR_FILE; // KEYWORD, SYSTEM_AUTO, PREFIX, SUFFIX, FILE, COMMAND, WRITE
 private:                                                   //
  void free();                                              // free space
  void copy(const _vflags& b);                              //
};

// --------------------------------------------------------------------------
// general flags for aims mode
class _aimsflags {
  public:
    _aimsflags();
    ~_aimsflags();
    _aimsflags(const _aimsflags& b);
    const _aimsflags& operator=(const _aimsflags& b);
    void clear();
    // CONTENT
    // in this struct we put all the flags which are used on LOCAL DIRECTORIES in AIMS MODE
    xoption KBIN_AIMS_FORCE_OPTION_NOTUNE;
    xoption KBIN_AIMS_RUN;
    xoption KBIN_AIMS_GEOM_MODE;
    xoption KBIN_AIMS_GEOM_FILE;
    std::vector<string> KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING;
    std::vector<xstructure> KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE;
    xoption KBIN_AIMS_GEOM_FILE_VOLUME;
    xoption KBIN_AIMS_FORCE_OPTION_VOLUME;
    xoption KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL;
    xoption KBIN_AIMS_CONTROL_MODE;
    xoption KBIN_AIMS_CONTROL_FILE;
    bool KBIN_AIMS_CONTROL_VERBOSE;   // VERBOSITY
  private:
    void free();                                              // free space
    void copy(const _aimsflags& b);                           //
};

// --------------------------------------------------------------------------
// general flags for alien mode
class _alienflags {
 public:
  // trivial constructurs/destuctors/operators
  _alienflags();                                            // default, just allocate
  ~_alienflags();                                           // kill everything
  _alienflags(const _alienflags& b);                        // constructor copy
  const _alienflags& operator=(const _alienflags &b);       // copy
  void clear(void);                                         // clear
  // CONTENT
  bool KBIN_ALIEN_COMMAND_BINARY_FLAG;
  string KBIN_ALIEN_COMMAND_BINARY_VALUE;
  bool KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG;
  // in this struct we put all the flags which are used on LOCAL DIRECTORIES in ALIEN MODE
  bool KBIN_ALIEN_FORCE_OPTION_NOTUNE;
  bool KBIN_ALIEN_FORCE_OPTION_SOMETHING;                   // SOMETHING

  bool KBIN_ALIEN_INPUT_MODE_EXPLICIT;
  bool KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP;
  bool KBIN_ALIEN_INPUT_MODE_IMPLICIT;
  bool KBIN_ALIEN_INPUT_MODE_EXTERNAL;
  bool KBIN_ALIEN_INPUT_FILE;
  bool KBIN_ALIEN_INPUT_FILE_FILE_FLAG;
  string KBIN_ALIEN_INPUT_FILE_FILE_VALUE;
  bool KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG;
  string KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE;
  bool KBIN_ALIEN_INPUT_MODE_INPUT_FLAG;
  string KBIN_ALIEN_INPUT_MODE_INPUT_VALUE;
  bool KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG;
  string KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE;
 private:                                              //
  void free();                                             // free space
  void copy(const _alienflags& b);                         //
};

// --------------------------------------------------------------------------
// general container for any set of flags
class _xflags {
  public:
    _xflags();
    _xflags(_vflags& vflags);
    _xflags(_aimsflags& aimsflags);
    _xflags(_alienflags& alienflags);
    ~_xflags();
    _xflags(const _xflags& b);
    const _xflags& operator=(const _xflags& b);
    void clear();
    bool AFLOW_MODE_VASP;
    _vflags vflags;
    bool AFLOW_MODE_AIMS;
    _aimsflags aimsflags;
    bool AFLOW_MODE_ALIEN;
    _alienflags alienflags;
    void setVFlags(_vflags& vflags);
    void setAIMSFlags(_aimsflags& aimsflags);
    //add qe and others here eventually
    void setALIENFlags(_alienflags& alienflags);
  private:
    void free();                                              // free space
    void copy(const _xflags& b);                           //
};

// --------------------------------------------------------------------------
// aflow_init.cpp
namespace init {
  int GetCPUCores();
  bool InitMachine(bool INIT_VERBOSE,vector<string>& argv,vector<string>& cmds,std::ostream& outf);
  string InitLoadString(string string2load,bool=FALSE);
  string InitGlobalObject(string string2load,string="",bool=FALSE);
  string InitLibraryObject(string string2load,bool=FALSE);
  long GetRAM(void);
  uint GetTEMPs(void);
} // namespace init
uint AFLOW_getTEMP(vector<string> argv);
uint AFLOW_monitor(vector<string> argv);
double AFLOW_checkMEMORY(string="",double=102.0);
bool CheckMaterialServer(string message);
bool CheckMaterialServer(void);
string aflow_get_time_string(void);
string aflow_get_time_string_short(void);
string strPID(void);

string Message(string="");
string Message(string str1,string list2print);
string Message(const _aflags& aflags,string="",string="");
bool AFLOW_BlackList(string h);
namespace init {
  bool ErrorOption(ostream &oss,const string& options, const string& routine,vector<string> vusage);
  bool ErrorOption(ostream &oss,const string& options, const string& routine,string vusage);
}

// --------------------------------------------------------------------------
// aflow_aflowrc.cpp
namespace aflowrc {
  bool is_available(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool read(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool write_default(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool print_aflowrc(std::ostream& oss,bool AFLOWRC_VERBOSE);
} // namespace aflowrc

// --------------------------------------------------------------------------
// aflow_arguments  
uint PflowARGs(vector<string> &argv,vector<string> &cmds,aurostd::xoption &vpflow); // called inside Init::InitMachine coded in aflow_pflow_main.cpp
uint ApennsyARGs(vector<string> &argv,vector<string> &cmds,aurostd::xoption &vflag); // called inside Init::InitMachine coded in aflow_apennsy_main.cpp

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// aflow_xatom.cpp

#define _COORDS_FRACTIONAL_ 0
#define _COORDS_CARTESIAN_  1
#define _UPDATE_LATTICE_VECTORS_TO_ABCANGLES_   2
#define _UPDATE_LATTICE_ABCANGLES_TO_VECTORS_   3

#define _PGROUP_ 0        // for point group lattice
#define _PGROUPK_ 5       // for point group klattice
#define _PGROUP_XTAL_ 6    // for point group crystal
#define _PGROUPK_XTAL_ 7  // for point group kcrystal
#define _FGROUP_ 1        // for factor group
#define _SGROUP_ 2        // for space group
#define _AGROUP_ 3        // for site positions point group
#define _IATOMS_ 4        // for equivalent atoms

// --------------------------------------------------------------------------
// DX and CO - START
extern double _SYM_TOL_; // tolerance control for isequal_RHT in atom class (RHT)
// DX and CO - END

class _atom { // simple class.. nothing fancy
 public:
  // constructor destructor                              // constructor/destructor
  _atom();                                               // default, just allocate
  ~_atom();                                              // kill everything
  const _atom& operator=(const _atom &b);                // copy
  void clear();
  // content                                             // content
  xvector<double> fpos;                                  // positions are with respect to ijk lattice cell
  xvector<double> cpos;                                  // so if fpos/cpos is outside a cell, you can shift
  xvector<double> corigin;                               // origin for convasp purposes
  xvector<double> coord;                                 // general coordinate for symmetry routines (RHT)
  vector<string> fpos_equation;                         // DX 20180607 - lattice equation for atomic position 
  vector<string> cpos_equation;                         // DX 20180607 - Cartesian equation for atomic position 
  double spin;                                           // spin along z in VASP MODE
  bool spin_is_given;                                    // TRUE if spin has been set // DX 9/21/17
  xvector<double> noncoll_spin;                          // non-collinear spin                // DX 12/5/17
  bool noncoll_spin_is_given;                            // TRUE if noncoll_spin has been set // DX 12/5/17
  double mass;                                           // mass 
  int    type;                                           // with bringincell, which adjust cpos/fpos and ijk as well
  string name;                                           // the name read from the INPUT
  int info;                                              // container for misc. information  // RHT
  string cleanname;                                      // a chemical clean version of the name
  bool   name_is_given;                                  // TRUE is atom name has been given
  int    atomic_number;                                  // 0 by defauls
  int    number;                                         // atom number reference for convasp, from zero to the sky
  string sd;                                             // ?
  bool   print_RHT;                                      // a printer for coord and name (general position)   // RHT
  bool   print_cartesian;                                // print frac or cartesian
  bool   verbose;                                        // verbose in printing
  xvector<int> ijk;                                      // xvector identifier of the lattice (but you must give str)
  bool   isincell;                                       // is in cell ? (i==j==k==0 ?)
  int    basis;                                          // identifier of position in the basis, from zero to the sky
  double reference;                                      // reference/measure for ordering
  int    ireference;                                     // sort of number in the list
  // for symmetry
  int    equivalent;                                     // points to the equivalent atom in the cell (-1 if irreducible)
  bool   is_inequivalent;                                // if atom is irreducible
  uint   num_equivalents;                                // say how many they are (only for is_inequivalent)
  uint   index_iatoms;                                   // if calculated on the xstructure, the index within iatoms for the identical atoms
  // for order parameter                                 // order parameter
  int    order_parameter_value;                          // order parameter
  bool   order_parameter_atom;                           // order parameter
  // for partial occupation                              // partial occupation
  double partial_occupation_value;                       // partial occupation
  bool   partial_occupation_flag;                        // partial occupation
  int shell;                                             // neighbour shell number
  // operators/functions                                 // operator/functions
  friend ostream& operator<<(ostream &,const _atom&);    // print
  void CleanName(void);                                  // function to clean up the name
  void CleanSpin(void);                                  // function to clean up the spin from EZ vasp script
 private:                                                //
  void free();                                           // free space
};

class _atom_reference_cmp {                              // sorting through reference
 public:
  bool operator()(const _atom& atom1,const _atom& atom2) const {
    return (bool) (atom1.reference<atom2.reference);}
};
class _atom_type_cmp {                                   // sorting through type
 public:
  bool operator()(const _atom& atom1,const _atom& atom2) const {
    return (bool) (atom1.type<atom2.type);}
};

#define NUM_ELEMENTS (103+1)  // up to Uranium
extern std::vector<string> atom_symbol_vec;             // store starting from ONE
extern std::vector<string> atom_name_vec;               // store starting from ONE
extern std::vector<double> atom_mass_vec;               // store starting from ONE
extern std::vector<double> atom_volume_vec;             // store starting from ONE
extern std::vector<int> atom_valence_iupac_vec;         // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
extern std::vector<int> atom_valence_std_vec;           // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
extern std::vector<double> atom_miedema_phi_star;       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> atom_miedema_nws;            // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> atom_miedema_Vm;             // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> atom_miedema_gamma_s;        // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> atom_miedema_BVm;            // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> atom_radius_vec;             // store starting from ONE - Saxena
extern std::vector<double> atom_radius_covalent_vec;    // store starting from ONE - Codero, Covalent radii revisited, DOI: 10.1039/b801115j // DX and CO - 9/4/17
extern std::vector<double> atom_electronegativity_vec;  // store starting from ONE - Saxena
extern std::vector<string> atom_crystal_vec;            // store starting from ONE - Ashcroft Mermin
extern std::vector<double> xray_scatt_vec;              // store starting from ONE
extern std::vector<double> pettifor_scale;              // store starting from ONE - Chemical Scale Pettifor Solid State Communications 51 31-34 1984
extern std::vector<double> pearson_coefficient;         // ME 181020

void atoms_initialize(void);
uint GetAtomNumber(const string& symbol);
std::string GetAtomName(const string& symbol);
std::string GetAtomName(const uint& atnum);
std::string GetAtomSymbol(const string& symbol);
std::string GetAtomSymbol(const uint& atnum);
double GetAtomMass(const string& symbol);  // in Kg
double GetAtomMass(const uint& atnum); // in Kg
double GetAtomComptonCrossSection(const string& symbol); // barn (1 barn = 1e-28 m^2)
double GetAtomComptonCrossSection(const uint& atnum); // barn (1 barn = 1e-28 m^2)
double GetAtomPhotoelectricCrossSection(const string& symbol);  // barn (1 barn = 1e-28 m^2)
double GetAtomPhotoelectricCrossSection(const uint& atnum);  // barn (1 barn = 1e-28 m^2)
double GetAtomVolume(const string& symbol);
double GetAtomVolume(const uint& atnum);
int GetAtomValenceIupac(const string& symbol);
int GetAtomValenceIupac(const uint& atnum);
int GetAtomValenceStd(const string& symbol);
int GetAtomValenceStd(const uint& atnum);
double GetAtomRadius(const string& symbol);
double GetAtomRadius(const uint& atnum);
double GetAtomRadiusCovalent(const string& symbol); // DX and CO - 9/4/17
double GetAtomRadiusCovalent(const uint& atnum); // DX and Co - 9/4/17
double GetAtomElectronegativity(const string& symbol);
double GetAtomElectronegativity(const uint& atnum);
string GetAtomCrystal(const string& symbol);
string GetAtomCrystal(const uint& atnum);
double GetAtomPettiforScale(const string& symbol);
double GetAtomPettiforScale(const uint& atnum);
bool GetAtomPettiforScale(const vector<string>& vsymbol,vector<double>& vvalue);
bool GetAtomPettiforScale(const vector<uint>& vatnum,vector<double>& vvalue);
bool GetAtomPettiforScale(const vector<string>& vsymbol,xvector<double>& vvalue);
bool GetAtomPettiforScale(const vector<uint>& vatnum,xvector<double>& vvalue);
bool SortAtomsPettiforScale(vector<string> &vsymbols,xvector<int> &vorders,xvector<double> &vvalues);
bool SortAtomsPettiforScale(vector<string> &vsymbols,vector<int> &vorders,vector<double> &vvalues);
bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<int> &vorder);
bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<double> &vvalue);
bool SortAtomsPettiforScale(vector<string> &vsymbol);
double GetPearsonCoefficient(const string&);
double GetPearsonCoefficient(const int&);
double GetAtomXrayScatt(const string& symbol);
double GetAtomXrayScatt(const uint& atnum);
double GetCompoundAttenuationLenght(const vector<string>& species,const vector<double>& composition,const double& density);  // density in g/cm^3, return in cm
double GetCompoundAttenuationLenght(const deque<string>& _species,const deque<int>& _composition,const double& density);  // density in g/cm^3, return in cm
// DX and CO - START
bool isequalRHT(const _atom& a, const _atom& b,double=_SYM_TOL_);       // bool equality only checks 'coord' and 'name' (RHT)  // RHT
// DX and CO - END
// routines of general use
string XATOM_AlphabetizationSpecies(string speciesA,string speciesB);
string XATOM_AlphabetizationSpecies(vector<string> vspecies);
string XATOM_AlphabetizationSpecies(vector<string> vspecies,vector<double> vnumbers);
void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies,vector<double>& vnumbers);
void XATOM_AlphabetizationCompound(string& system, vector<string>& vspecies,vector<double>& vnumbers);
void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies);
void XATOM_AlphabetizationSpecies(string& system);
void XATOM_AlphabetizationCompound(string& system);
uint XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX);
uint XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX);
uint XATOM_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX);
uint XATOM_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX, vector<double> &natomsX);
// neighbour things
void GetUnitCellRep(const xvector<double>& ppos,xvector<double>& p_cell0,xvector<int>& ijk,const xmatrix<double>& lattice,const bool coord_flag);

string xstructure2json(xstructure& xstr); // DX 8/31/17 - xstructure2json
string atom2json(_atom& atom, int coord_flag, int poccupation); // DX 8/31/17 - atom2json

// --------------------------------------------------------------------------
class _sym_op {
 public:
  // constructor destructor
  _sym_op();                                                    // default, just allocate
  ~_sym_op();                                                   // kill everything
  // content
  // for _PGROUP_
  xmatrix<double>  Uc;            // 3x3                        // uniques (not irreducible) operations on positions (Uc cartesian)
  xmatrix<double>  Uf;            // 3x3                        // uniques (not irreducible) operations on indices   (Uf fractional)
  xmatrix<double>  generator;     // 3x3                        // generator A, U=exp(A*theta)
  xvector<double>  generator_coefficients;                      // generator coefficients on Lx, Ly, Lz basis // DX 12/6/17
  xmatrix<xcomplex<double> > SU2_matrix; // 2x2                 // SU(2) 2x2 complex matrix // DX 1/15/18
  xvector<xcomplex<double> > su2_coefficients;                  // su(2) coefficients on sigma_1, sigma_2, sigma_3 basis (Pauli matrices) // DX 1/15/18
  double           angle;                                       // angle axis
  xvector<double>  axis;          // 3                          // (1,2,3)=axis
  xvector<double> quaternion_vector;				//GEENA
  xmatrix<double> quaternion_matrix;				//GEENA
  string           str_type;                                    // generic type of the operation
  string           str_Hermann_Mauguin;                         // Hermann_Mauguin notation
  string           str_Schoenflies;                             // Schoenflies notation
  bool             flag_inversion;                              // flag if inversion
  vector<int>      basis_atoms_map;                             // this is the vector that tell where the basis atom gets mapped by the operation
  vector<int>      basis_types_map;                           // this is the vector that tell where the basis species gets mapped by the operation
  // DX and CO - START
  bool             basis_map_calculated;                        //have we've calculated it?
  // DX and CO - END
  bool             is_pgroup;                                   // bool is_pgroup
  // for _PGROUPXTAL_
  bool             is_pgroup_xtal;                              // bool is_pgroup_xtal
  // for _PGROUPK_
  bool             is_pgroupk;                                  // bool is_pgroupk
  // for _PGROUPK_XTAL_                  
  bool             is_pgroupk_xtal;                             // bool is_pgroupk_xtal // DX 12/5/17 - Added pgroupk_xtal
  // for _FGROUP_
  xvector<double>  ctau;          // 3                          // translation in CARTESIAN       // FACTOR GROUP only, [0,1[
  xvector<double>  ftau;          // 3                          // translation in FRACTIONAL      // FACTOR GROUP only, [0,1[
  bool             is_fgroup;                                   // bool is_fgroup
  // for _SGROUP_
  xvector<double>  ctrasl;        // 3                          // translation in CARTESIAN       // SPACE GROUP only, [integers]
  xvector<double>  ftrasl;        // 3                          // translation in FRACTIONAL      // SPACE GROUP only, [ingegers]
  bool             is_sgroup;                                   // bool is_sgroup
  // operators
  // for _AGROUP_
  bool             is_agroup;                                   // bool is_agroup     // for site operation point group
  uint             site;                                        // uint site          // site index // DX 8/3/17
  const _sym_op& operator=(const _sym_op& b);
  friend ostream& operator<<(ostream &,const _sym_op&);
 private:
  void free();
};

//DX 201801107 - add _kpoint class - START
// --------------------------------------------------------------------------
class _kpoint {
 public:
  // constructor destructor
  _kpoint();                                           // default, just allocate
  ~_kpoint();                                          // default, just allocate
  // content
  char iomode;                                         // store format (not used yet)
  xmatrix<double> klattice;                            // reciprocal lattice
  xvector<double> fpos;                                // fractional position of kpoint
  xvector<double> cpos;                                // Cartesian position of kpoint (not used yet)
  string label;                                        // kpoint label (i.e., high-symmetry point labels)
  bool is_transformed;                                 // indicates if kpoint is transformed from AFLOW standard
  const _kpoint& operator=(const _kpoint& b);          // assignment operator
  // operators/functions                               // operator/functions
  string str() const;                                  // prints "fpos ! label" (e.g., 0.0000 0.0000 0.0000 ! \\Gamma)
  void TransformKpoint(const xmatrix<double>& P);      // transforms kpoint via P matrix (k'=k*P) and klattice via Q matrix (L_recip'=Q*L_recip) (see ITC-A pg. 79)
  friend ostream& operator<<(ostream&,const _kpoint&); // ostream operator
 private:
  void free();
};
//DX 201801107 - add _kpoint class - END

// --------------------------------------------------------------------------
class wyckoffsite_ITC { //Also for wyckoff sites
 public:
  wyckoffsite_ITC(void);
  wyckoffsite_ITC(const wyckoffsite_ITC& b);
  ~wyckoffsite_ITC(void);
  // OPERATORS                                                  // --------------------------------------
  const wyckoffsite_ITC& operator=(const wyckoffsite_ITC& b);             // some operators
  friend ostream& operator<<(ostream&,const wyckoffsite_ITC&);       // ostream
  // CONTENT
  xvector<double> coord;
  string type; //chemical label etc
  string wyckoffSymbol;
 private:                                                       // ---------------------------------------
  void free();                                                  // to free everything
};

// --------------------------------------------------------------------------

#define MAX_TITLE_SIZE 512

#define IOAFLOW_AUTO   0
#define IOVASP_AUTO    1
#define IOVASP_POSCAR  2
#define IOVASP_ABCCAR  3
#define IOVASP_WYCKCAR 4
#define IOQE_AUTO      5
#define IOQE_GEOM      6
#define IOABINIT_AUTO  7
#define IOABINIT_GEOM  8
#define IOAIMS_AUTO    9
#define IOAIMS_GEOM   10
#define IOCIF         11 //DX 20180723

#define NOSG string("NNN #0")

#define _EQUIV_FPOS_EPS_    2.0e-5    // NOV 2009 Israel  used to be 1.0e-6 too small for ICSD
#define _pocc_no_sublattice_ -1
#define DEFAULT_PARTIAL_OCCUPATION_TOLERANCE 0.02

bool sortAtomsTypes(const _atom& a1,const _atom& a2);		// sort atoms by types
bool sortAtomsNames(const _atom& a1,const _atom& a2);		// sort atoms by names
bool sortAtomsDist(const _atom& a1,const _atom& a2);		// sort atoms by dist  // CO 180420

class xstructure {
 public:
  // constructors/destructors                                   // --------------------------------------
  xstructure(string="");                                        // constructor default
  xstructure(const xstructure& b);                              // constructor copy
  xstructure(istream& input,int=IOVASP_POSCAR);                 // constructor from istream
  xstructure(ifstream& input,int=IOVASP_POSCAR);                // constructor from ifstream
  xstructure(stringstream& input,int=IOVASP_POSCAR);            // constructor from stringstream
  xstructure(const string& input,int);                          // constructor from file
  xstructure(const string& url,const string& file,int=IOVASP_POSCAR); // constructor from URL
  ~xstructure();                                                // destructor
  // I/O, mutators                                              // --------------------------------------
  bool GetStoich(void);                                         // get stoich_each_type - 170724 CO
  bool FixLattices(void);                                       // Reciprocal/f2c/c2f
  void SetCoordinates(const int& mode);                         // change coordinates
  void MakeBasis(void);                                         // make basis for atoms (basis and number)
  void MakeTypes(void);                                         // refresh types based on num_each_type  // CO 180420
  void AddAtom(const _atom& atom);                              // adding an atom
  void AddAtom_POCC(const _atom& atom);                         // adding an atom FOR POCC ONLY
  void RemoveAtom(const uint& iat);                             // deleting an atom (index)
  void RemoveAtoms(vector<uint>& v_atoms_to_remove);            // deleting many atoms (indices)
  void RemoveCopies(double=1.0e-3);                             // deleting atoms too close F/C
  void RemoveFractionalCopies(double=1.0e-3);                   // deleting atoms too close F
  void RemoveCartesianCopies(double=1.0e-3);                    // deleting atoms too close C
  void AddCorners(void);                                        // for picturing purpose
  void Clear(void);                                             // clear everything
  void ClearSpecies(void);                                     // Clear all the symmetry
  void ShifOriginToAtom(const int& iat);                        // Shift the origin to atom(iat)
  void IdenticalAtoms(void);                                    // Make identical atoms
  void SwapCoordinates(const uint& i,const uint& j);            // Permute Coordinates i with j
  string SpeciesLabel(const uint& A);                           // Returns the Label of the specie A (if available)
  void SpeciesSwap(const uint& A,const uint& B);                // Permute Species A with B (safe for species C).
  bool SpeciesGetAlphabetic(void);                              // Check is species are in alphabetic order
  bool SpeciesPutAlphabetic(void);                              // Put Species in alphabetic
  string SpeciesString(void);                                   // Gives a string with the list of all the species
  uint SetSpecies(std::deque<string> vspecies);                 // Set the species
  void GetLatticeType(xstructure& sp,xstructure& sc);           // Get all lattices
  void GetLatticeType(void);                                    // Get all lattices
  void Standard_Primitive_UnitCellForm(void);                   // Reduce the Unit Cell to Standard Primitive Form
  void GetStandardPrimitive(void);                              // stub for void Standard_Primitive_UnitCellForm(void);
  void Standard_Conventional_UnitCellForm(void);                // Reduce the Unit Cell to Standard Conventional Form
  void GetStandardConventional(void);                           // stub for void Standard_Conventional_UnitCellForm(void);
  void NiggliUnitCellForm(void);                                // Reduce the Unit Cell to Niggli Form
  void MinkowskiBasisReduction(void);                           // Reduce the Basis to the max orthogonality (Minkowski)
  void LatticeReduction(void);                                  // Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
  void BringInCell(void);                                       // Bring all the atoms in the origin
  void BringInCell(double);                                     // Bring all the atoms in the origin
  void BringInCompact(void);                                    // Bring all the atoms near the origin
  void BringInWignerSeitz(void);                                // Bring all the atoms in the Wigner Seitz Cell
  void GetPrimitive(void);                                      // Make it primitive, if possible
  void GetPrimitive(double tol);                                // Make it primitive, if possible
  void GetPrimitive1(void);                                     // Make it primitive, if possible
  void GetPrimitive2(void);                                     // Make it primitive, if possible
  void GetPrimitive3(void);                                     // Make it primitive, if possible
  uint GetPrimitiveCell(void);                                  // Make it primitive, if possible. Returns 1 if routine fails (RHT)   // RHT
  double MinDist(void);                                         // get minimum interatomic distance -- CO 171024
  void ReScale(const double &in_scale);                         // Change scale but keep volume fixed
  void SetScale(const double &in_scale);                        // Change scale
  void SetVolume(const double &in_volume);                      // Change volume
  void InflateLattice(const double &coefficient);               // Inflate lattice
  void InflateVolume(const double &coefficient);                // Inflate volume
  string platon2print(bool,bool,double,double,double,double);   // Create Platon input file >=51108
  void FakeNames(void);                                         // Fix names as fakes - useful for platon
  string platon2sg(bool P_EQUAL=DEFAULT_PLATON_P_EQUAL,
		   bool P_EXACT=DEFAULT_PLATON_P_EXACT,
		   double P_ang=DEFAULT_PLATON_P_ANG,
		   double P_d1=DEFAULT_PLATON_P_D1,
		   double P_d2=DEFAULT_PLATON_P_D2,
		   double P_d3=DEFAULT_PLATON_P_D3);
  string findsym2sg(double tolerance=DEFAULT_FINDSYM_TOL);
  string findsym2execute(double tolerance=DEFAULT_FINDSYM_TOL);
  string findsym2print(double tolerance=DEFAULT_FINDSYM_TOL);
  //  string platon2sg(void);
  double GetVolume(void);                                       // Return volume
  double Volume(void);                                          // Return volume
  double GetZVAL(const vector<double>& vZVAL);                  // Given the ZVAL of each species, it returns total ZVAL of cell
  double GetPOMASS(const vector<double>& vPOMASS);              // Given the POMASS of each species, it returns total POMASS of cell
  void ClearSymmetry(void);                                     // Clear all the symmetry
  bool CalculateSymmetry(bool,double);                          // Calculate the symmetry
  bool CalculateSymmetry(void);                                 // Calculate the symmetry
  void CalculateSymmetryPointGroup(bool);                       // Calculate the symmetry
  void CalculateSymmetryPointGroup(void);                       // Calculate the symmetry
  void CalculateSymmetryFactorGroup(bool);                      // Calculate the symmetry
  void CalculateSymmetryFactorGroup(void);                      // Calculate the symmetry
  void CalculateSymmetryPointGroupCrystal(bool);                // Calculate the symmetry
  void CalculateSymmetryPointGroupCrystal(void);                // Calculate the symmetry
  void CalculateSymmetryPointGroupKlattice(bool);               // Calculate the symmetry
  void CalculateSymmetryPointGroupKlattice(void);               // Calculate the symmetry
  int  GenerateGridAtoms(int,int,int,int,int,int);              // generate grid of atoms
  int  GenerateGridAtoms(int,int,int);                          // generate grid of atoms
  int  GenerateGridAtoms(int);                                  // generate grid of atoms
  int  GenerateLIJK(double);                                    // generate lijk look up table
  // QUANTUM ESPRESSO AND ABINIT AND AIMS                       // --------------------------------------
  void buildGenericTitle(bool vasp_input=false,bool force_fix=false); // build a nice title with atoms
  void xstructure2qe(void);                                     // some wrap up IOs to convert format to QE
  void xstructure2vasp(void);                                   // some wrap up IOs to convert format to VASP
  void xstructure2abinit(void);                                 // some wrap up IOs to convert format to ABINIT
  void xstructure2aims(void);                                   // some wrap up IOs to convert format to AIMS
  //[CO 180420 - moved outside of xstructure]bool sortAtomsTypes(const _atom& a1,const _atom& a2);		// sort atoms by types
  //[CO 180420 - moved outside of xstructure]bool sortAtomsNames(const _atom& a1,const _atom& a2);		// sort atoms by names
  // OPERATORS                                                  // --------------------------------------
  const xstructure& operator=(const xstructure& b);             // some operators
  friend istream& operator>>(istream&,xstructure&);             // istream
  friend ostream& operator<<(ostream&,const xstructure&);       // ostream
  // CONTENT                                                    // --------------------------------------
  string title;                                                 // Title of the structure
  string directory;                                             // Directory where xstructure came from // DX
  string prototype;                                             // Prototype of the structure
  string info;                                                  // Info of the structure
  int iomode;                                                   // IOVASP_POSCAR/IOXXXX
  // int num_types=num_each_type.size();                        // old useless stuff
  // int num_atoms=atoms.size();                                // old useless stuff
  bool neg_scale;                                               // flag for negative scale (for printing)
  double scale;                                                 // scale (always linear A)
  bool neg_scale_second;                                        // POCC (hnf vs. tol) // CO 180409
  double scale_second;                                          // POCC hnf/stoich tol/site tol // CO 180409
  aurostd::xoption scale_third;                                 // if there is a third scale number provided, we use isentry and content_double // CO 170803 - site tol
  char coord_type[2];                                           // type of coordinates
  bool coord_flag;                                              // _COORDS_FRACTIONAL_ (0) fractional, _COORDS_CARTESIAN_ (1) cartesian.
  bool isd;                                                     // TRUE=Selective dynamics; FALSE=no selective dynamics.
  xmatrix<double> lattice;                                      // LATTICE in REAL SPACE (meters)            // vector per RAW (must trasp per algebra)
  double a,b,c,alpha,beta,gamma;                                // LATTICE in a,b,c,alpha,beta,gamma
  xmatrix<double> klattice;                                     // LATTICE in MOMENTUM SPACE (1/meters)      // vevror per RAW (must trasp per algebra)
  xvector<double> origin;                                       // origin
  xmatrix<double> f2c;                                          // transformation matrix for F vector per COLUM f2c=trasp(lattice)
  xmatrix<double> c2f;                                          // transformation matrix for C vector per ROW   c2f=inv(trasp(lattice))
  double equiv_fpos_epsilon;                                    // when they are the same DEFAULT _EQUIV_FPOS_EPS_
  std::deque<int> num_each_type;                                // WARNING: we use starting from 0
  std::deque<double> comp_each_type;                            // WARNING: we use starting from 0
  std::deque<double> stoich_each_type;                          // WARNING: we use starting from 0 - 170724
  std::deque<_atom> atoms;                                      // WARNING: we use starting from 0
  std::deque<string> species,species_pp,species_pp_type,species_pp_version; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<double> species_pp_ZVAL; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<std::deque<double> > species_pp_vLDAU;             // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<double> species_volume;                            // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<double> species_mass;                              // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  //  ----------------------------------------------------------------------------------------
  // SYMBOLIC MATH stuff
  bool symbolic_math_representation_only;                       // print symbolic math representation only //DX 20180618 
  bool constrained_symmetry_calculation;                        // append symbolic math representation for constrained symmetry calculation //DX 20180618 
  vector<vector<string> > symbolic_math_lattice;                // symbolic math representation of lattice //DX 20180618 
  uint num_parameters;                                          // number of parameters ANRL 20180618
  uint num_lattice_parameters;                                  // number of lattice parameters ANRL 20180618
  vector<string> prototype_parameter_list;                      // prototype parameter list ANRL 20180618
  vector<double> prototype_parameter_values;                    // prototype parameters values ANRL 20180618
  //  ----------------------------------------------------------------------------------------
  bool is_vasp4_poscar_format;                                  // flags for VASP4*
  bool is_vasp5_poscar_format;                                  // flags for VASP5*
  bool Niggli_calculated;                                       // flags for calculation
  bool Niggli_avoid;                                            // flags for avoiding the calculation
  bool Minkowski_calculated;                                    // flags for calculation
  bool Minkowski_avoid;                                         // flags for avoiding the calculation
  bool LatticeReduction_calculated;                             // flags for calculation
  bool LatticeReduction_avoid;                                  // flags for avoiding the calculation
  //  ----------------------------------------------------------------------------------------
  // PRINTING stuff
  string PrintSymbolicMathRepresentation(void);                 // Print symbolic math representation of structure //DX 20180618
  string PrintUNCLE(void);                                      // Print in UNCLE format
  //  ----------------------------------------------------------------------------------------
  // LATTICE stuff
  bool Standard_Lattice_calculated;                             // flags for calculation
  bool Standard_Lattice_avoid;                                  // flags for avoiding the calculation
  bool Standard_Lattice_primitive;                              // flags for calculation
  bool Standard_Lattice_conventional;                           // flags for calculation
  bool Standard_Lattice_has_failed;                             // flags for Lattice has failed ?
  string bravais_lattice_type;                                  // lattice type as a string  (14)
  string bravais_lattice_variation_type;                        // lattice type as a string wahyu mod  (with the mods of wahyu)
  string bravais_lattice_system;                                // lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
  string bravais_lattice_lattice_type;                          // lattice_lattice type as a string  (14)
  string bravais_lattice_lattice_variation_type;                // lattice_lattice type as a string wahyu mod  (with the mods of wahyu)
  string bravais_lattice_lattice_system;                        // lattice_lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
  string pearson_symbol;                                        // pearson symbol as a string
  string reciprocal_lattice_type;                               // reciprocal lattice type as a string
  string reciprocal_lattice_variation_type;                     // reciprocal lattice type as a string wahyu mod
  //string reciprocal_conventional_lattice_type;                // reciprocal lattice type as a string
  string bravais_superlattice_type;                             // super lattice type as a string (identical atoms)
  string bravais_superlattice_variation_type;                   // super lattice type as a string (identical atoms) wahyu mod
  string bravais_superlattice_system;                           // lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
  string pearson_symbol_superlattice;                           // pearson symbol of the superlattice (identical atoms)
  bool volume_changed_original2new;                             // flag for volume has changed between original and new (i.e., transformation won't work) //DX 20181105
  xmatrix<double> transform_coordinates_original2new;           // transform coordinate system from original to new; (Q in ITC notation) //DX 20181105
  xmatrix<double> transform_coordinates_new2original;           // transform coordinate system from new to original; (Q^-1 in ITC notation) //DX 20181105
  xmatrix<double> rotate_lattice_original2new;                  // rotate from original to new lattice; (P in ITC notation) //DX 20181105
  xmatrix<double> rotate_lattice_new2original;                  // rotate from new to original lattice; (P^-1 in ITC notation) //DX 20181105
  //  ----------------------------------------------------------------------------------------
  // GENERAL PURPOSE LABELS                                     // general purpose label
  uint label_uint;                                              // general purpose label_uint
  int label_int;                                                // general purpose label_int
  double label_double;                                          // general purpose label_double
  // ----------------------------------------------------------------------------------------
  // ORDER PARAMETER                                            // order parameter for xstructure
  bool order_parameter_structure;                               // order parameter for xstructure
  std::vector<uint> order_parameter_atoms;                      // indices of atoms to be shuffled
  uint order_parameter_orbit;                                   // number of equivalent configurations with the factor group
  int order_parameter_sum;                                      // sum of all the order parameters
  // ----------------------------------------------------------------------------------------
  // PARTIAL OCCUPATION                                         // partial occupation for xstructure
  bool partial_occupation_flag;                                 // flags for partial occupation TRUE/FALSE
  double partial_occupation_site_tol;                           // tolerance for partial occupation site >=0.0 <=1.0   // CO 180409
  double partial_occupation_stoich_tol;                         // tolerance for partial occupation stoich >=0.0 <=1.0 // CO 180409
  int partial_occupation_HNF;                                   // volume HNF size
  vector<int> partial_occupation_sublattice;                    // contains the information about the specie# of the sublattice in the partial occupation otherwise _pocc_no_sublattice_
  // ----------------------------------------------------------------------------------------
  // GEOMETRY ENERGETICS after the QM calculations              // --------------------------------------
  void qm_clear(void);                                          // QM create/clean all the vectors
  void qm_recycle(void);                                        // QM shift data from QM to GEOM
  void qm_load(string directory,string suffix="",int=IOVASP_POSCAR);                    // QM results load from an ab-initio calculation
  bool qm_calculated;                                           // QM calculation
  double qm_scale;                                              // QM scale (always linear A)
  xmatrix<double> qm_lattice;                                   // QM LATTICE in REAL SPACE (meters)
  xmatrix<double> qm_klattice;                                  // QM LATTICE in MOMENTUM SPACE (1/meters)     // SAVED TRASP
  xvector<double> qm_origin;                                    // QM origin
  xmatrix<double> qm_f2c;                                       // QM transformation matrix for F vector per COLUM f2c=trasp(lattice)
  xmatrix<double> qm_c2f;                                       // QM transformation matrix for C vector per ROW   c2f=inv(trasp(lattice))
  std::deque<_atom> qm_atoms;                                  // QM WARNING: we use starting from 0
  std::vector<xvector<double> > qm_forces;                      // QM FORCES calculation
  bool qm_forces_write;                                         // QM FORCES calculation
  std::vector<xvector<double> > qm_positions;                   // QM POSITIONS calculation
  bool qm_positions_write;                                      // QM POSITIONS calculation
  double qm_E_cell,qm_dE_cell,qm_H_cell,qm_PV_cell,qm_P,qm_mag_cell;                     // QM energetics PER CELL
  double qm_E_atom,qm_dE_atom,qm_H_atom,qm_PV_atom,qm_mag_atom;                     // QM energetics ATOMIC
  // ----------------------------------------------------------------------------------------
  // KPOINTS                                                    // --------------------------------------
  int kpoints_mode;                                             // mode of kpoints
  int    kpoints_k1,kpoints_k2,kpoints_k3;                      // parameters that are plug during
  double kpoints_s1,kpoints_s2,kpoints_s3;                      // parameters that are plug during
  int kpoints_kmax,kpoints_kppra;                               // load/unload and calculations
  string kpoints_kscheme;                                       // of ab-initio
  // ---------------------- SYMMETRY --------------------------------------------------------
  // A=U*B but in A and B we plug vectors as columns watch lattice is per row Uc=A*inv(B)
  // A is the lattice (vectors per colum), B is the test lattice (epr column)
  // Uc is the point_group operation which operates AFTER the vector (row)
  // as: new_vector_row=old_vector_row*Uc  and point_group is the list of all the Uc !!!
  // DX and CO - START
  double dist_nn_min; 
  // SYMMETRY TOLERANCE ----------------------------
  bool sym_eps_calculated;                                      // was it calculated automatically per symmetry operations (aflowSYM)?
  double sym_eps;                                               // universal tolerance for symmetry (dictates resolution and mapping tolerances)                     
  uint sym_eps_change_count;                                  // universal tolerance count for symmetry // DX 2/23/18 - added count to xstructure
  // DX and CO - END
  // POINT GROUP                                                // POINT GROUP LATTICE
  std::vector<_sym_op> pgroup;                                  // rotations/inversions operations
  bool pgroup_calculated;                                       // WARNING: we use starting from 0
  // POINT GROUP CRYSTAL                                        // POINT GROUP CRYSTAL
  std::vector<_sym_op> pgroup_xtal;                             // rotations/inversions operations
  bool pgroup_xtal_calculated;                                  // WARNING: we use starting from 0
  string crystal_family;                                        // crystal and point group properties
  string crystal_system;                                        // crystal and point group properties
  string point_group_crystal_class;                             // crystal and point group properties
  string point_group_Shoenflies;                                // crystal and point group properties
  string point_group_Hermann_Mauguin;                           // crystal and point group properties
  string point_group_orbifold;                                  // crystal and point group properties
  string point_group_type;                                      // crystal and point group properties
  string point_group_order;                                     // crystal and point group properties
  string point_group_structure;                                 // crystal and point group properties
  // POINT GROUP KLATTICE                                       // POINT GROUP
  std::vector<_sym_op> pgroupk;                                 // rotations/inversions operations
  bool pgroupk_calculated;                                      // WARNING: we use starting from 0
  // POINT GROUP KCRYSTAL                                       // POINT GROUP
  std::vector<_sym_op> pgroupk_xtal;                            // rotations/inversions operations
  bool pgroupk_xtal_calculated;                                 // WARNING: we use starting from 0
  // FACTOR GROUP                                               // FACTOR GROUP
  std::vector<_sym_op> fgroup;                                  // rotations/inversions + incell_translations operations
  bool fgroup_calculated;                                       // WARNING: we use starting from 0
  // SPACE GROUP                                                // SPACE GROUP
  std::vector<_sym_op> sgroup;                                  // rotations/inversions + incell        //outcell_translations operations
  bool sgroup_calculated;                                       // WARNING: we use starting from 0
  double sgroup_radius;                                         // radius of application (all ops connecting objects inside the sphere)
  xvector<int> sgroup_radius_dims;                              // dimension of the radius (in +- integers)
  // SITE POINT GROUP                                           // SITE POINT GROUP
  bool agroup_calculated;                                       // WARNING: we use starting from 0
  std::vector<std::vector<_sym_op> > agroup;                    // rotations/inversions operations on sites, has one for each atom (all)
  // INEQUIVALENTE ATOMS                                        // --------------------------------------
  bool iatoms_calculated;                                       // given the symmetry, the atoms are mapped in inequivalent
  std::vector<std::vector<int> > iatoms;                        // equivalent/inequivalent atoms lookup table
  // SPACE GROUP CALCULATION WITH PLATON/FINDSYM                // with platon >= 51108
  string spacegroup;                                            // space group symbol
  string spacegrouplabel;                                       // the number with #
  string spacegroupoption;                                      // origin, axes and so on
  int    spacegroupnumber;                                      // the number
  int    spacegroupnumberoption;                                // the option as number
  bool is_spacegroup_platon,is_spacegroup_findsym,is_spacegroup_aflow; // got spacegroup
  // SPACE GROUP ITC
  // DX and CO - START
  string crystal_system_ITC;                                    // aflow_symmetry_spacegroup.cpp (RHT)
  string point_group_ITC;                                       // aflow_symmetry_spacegroup.cpp (RHT)
  char bravais_label_ITC;                                       // aflow_symmetry_spacegroup.cpp (RHT)
  char lattice_label_ITC;                                       // aflow_symmetry_spacegroup.cpp (RHT)
  uint space_group_ITC;                                         // aflow_symmetry_spacegroup.cpp (RHT)
  string wyckoff_library_entry_ITC;                             // aflow_symmetry_spacegroup.cpp (RHT)
  int setting_ITC;                                              // aflow_symmetry_spacegroup.cpp (RHT) // DX 8/30/17 - SGDATA
  xvector<double> origin_ITC;                                   // aflow_symmetry_spacegroup.cpp (RHT) // DX 8/30/17 - SGDATA
  vector<string> general_position_ITC;                          // aflow_symmetry_spacegroup.cpp (RHT) // DX 8/30/17 - SGDATA
  // DX and CO - END
  // double volume;  USE double GetVolume(const xstructure& a);  and SetVolume //
  //  int number_of_atoms; USE (int) atoms.size();  looks you do not need.
  // DX and CO - START
  vector<string> wyccar_ITC;                                    // aflow_symmetry_spacegroup.cpp (RHT)
  xmatrix<double> standard_lattice_ITC;                         // aflow_symmetry_spacegroup.cpp (RHT)
  deque<_atom> standard_basis_ITC;                              // aflow_symmetry_spacegroup.cpp (RHT)    
  vector<wyckoffsite_ITC> wyckoff_sites_ITC;                    // aflow_symmetry_spacegroup.cpp (RHT) //(x,y,z) XX
  vector<string> wyckoff_symbols_ITC;                           // aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(void);                                    // aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol);                         // aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol, bool& no_scan);          // aflow_symmetry_spacegroup.cpp (RHT)
  //uint SpaceGroup_ITC(double& use_tol, const int& manual_it);// aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol, const int& setting);// aflow_symmetry_spacegroup.cpp (RHT) //DX 20180806
  // DX 9/5/17 - [OBSOLETE] uint SpaceGroup_ITC(double& use_tol, double& orig_tolerance, const int& manual_it,int& change_sym_count,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT)
  //uint SpaceGroup_ITC(double& use_tol,const int& manual_it,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol,const int& setting,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT) //DX 20180806
  uint SpaceGroup_ITC(double& use_tol,const int& manual_it,const int& setting,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT) //DX 20180806
  string aflow2sg(void);                                        // aflow_symmetry_spacegroup.cpp (DX)
  string aflow2sg(double& use_tol);                             // aflow_symmetry_spacegroup.cpp (DX)
  string aflow2sg(double& use_tol, const int& manual_it);       // aflow_symmetry_spacegroup.cpp (DX)
  // DX and CO - END
  // ---------------------- FROZSL ---------------------- 
  // [OBSOLETE]  string prototype_label();
  ostringstream FROZSL_output(vector<string> Kvectors);
  // ---------------------- PHONONS ---------------------------------------------------------
  // based on the Maradudin harmonic and deformation analysis
  // LIJK OBEJCTS                                               // LIJK OBEJCTS   WORKING
  bool lijk_calculated;                                         // is calculated ?
  vector<xvector<int> >    lijk_table;                          // bravais lattice look up table l <-> i,j,k
  vector<xvector<double> > lijk_fpos;                           // bravais lattice look up table fpos same as lijk_table
  vector<xvector<double> > lijk_cpos;                           // bravais lattice look up table cpos
  xvector<int> lijk_dims;                                       // dimension
  // ----------------------------------------------------------------------------------------
  // GRID ATOMS                                                 // --------------------------------------
  bool grid_atoms_calculated;                                   // GRID ATOMS from dimsL to dimsH has been calculated ?
  xvector<int> grid_atoms_dimsL;                                // dims low of i,j,k (minus is NOT included)
  xvector<int> grid_atoms_dimsH;                                // dims high of i,j,k (plus is NOT included)
  std::deque<_atom> grid_atoms;                                // WARNING: we use starting from 0
  int grid_atoms_number;                                        // how many...
  vector<int> grid_atoms_sc2pcMap;                              // 170804 CO - mapping between grid_atoms (sc) and atoms (pc)
  vector<int> grid_atoms_pc2scMap;                              // 170804 CO - mapping between grid_atoms (sc) and atoms (pc)
  // ----------------------------------------------------------------------------------------
  // NEIGHBOURS OBEJCTS EXPERIMENTAL/UNFINISHED                 // NEIGHBOURS OBEJCTS   WORKING EXPERIMENTAL
  void neighbours_clear(void);                                  // NN create/clean all the vectors
  void neighbours_calculate(void);                              // NN shift data from QM to GEOM
  bool neighbours_calculated;                                   // NN calculation
  double neighbours_radius;                                     // radius of application
  double neighbours_dradius;                                    // delta radius of shell (bins)
  vector<vector<double> > neighbours_atoms_func_r_vs_nn;        // contains function distance vs neighbours (each atom)
  vector<vector<int> > neighbours_atoms_func_num_vs_nn;         // contains function distance vs neighbours (each atom)
  vector<double> neighbours_func_r_vs_nn;                       // contains function distance vs neighbours (all atoms)
  vector<int> neighbours_func_num_vs_nn;                        // contains function number vs neighbours (all atoms)
  // xvector<int> ndims;                                        // dimension of the radius (in +- integers)
  // std::deque<_atom> ashell;                                 // all the atoms in the shell
  // std::deque<deque<_atom> > natoms;                        // vector of vectors
  // std::vector<vector<double> > rshell;                       // vector of shells
  // std::vector<vector<int> > nshell;                          // vector of density in shells
  // int nbins;                                                 // number of bins
  // NEIGHBOURS OBEJCTS OLD-ACONVASP BUT WORKS                  // NEIGHBOURS OBEJCTS 
  // GetNeighData collects all the neighbor data between rmin and rmax and stores it for each atom in a vector of atom objects in order of increasing distance.  
  void GetNeighData(const deque<_atom>& in_atom_vec,const double& rmin, const double& rmax,deque<deque<_atom> >& neigh_mat);
  // GetStrNeighData collects all the neighbor data out to some cutoff and stores it for each atom in the structure.
  void GetStrNeighData(const double cutoff,deque<deque<_atom> >& neigh_mat);
  // ----------------------------------------------------------------------------------------
  // OUTPUT/ERROR FLAGS                                         // --------------------------------------
  bool Niggli_has_failed;                                       // Niggli has failed ?
  bool Minkowski_has_failed;                                    // Minkowski has failed ?
  bool LatticeReduction_has_failed;                             // LatticeReduction has failed ?
  bool write_lattice_flag;                                      // flag for OUTPUT printing
  bool write_klattice_flag;                                     // flag for OUTPUT printing
  bool write_inequivalent_flag;                                 // flag for OUTPUT printing
  bool write_DEBUG_flag;                                        // flag for OUTPUT printing
  bool error_flag;                                              // flag TRUE is error
  string error_string;                                          // contains type of error
  // END OF CONTENT                                             //
 private:                                                       // ---------------------------------------
  void Free();                                                  // to free everything
  void Copy(const xstructure& b);                               // the flag is necessary because sometimes you need to allocate the space.
};

// CO 180420
//for stream management with objects
class xStream {
  public:
    //NECESSARY PUBLIC CLASS METHODS - START
    //constructors - START
    xStream();
    //constructors - STOP
    ~xStream();
    //NECESSARY PUBLIC CLASS METHODS - END
  protected:
    //NECESSARY private CLASS METHODS - START
    void free();
    void freeStreams();
    void freeAll();
    void copyStreams(const xStream& b);
    //NECESSARY END CLASS METHODS - END
    //logger variables
    ostream* p_oss;
    ofstream* p_FileMESSAGE;
    bool m_new_ofstream;  //for deletion later
    //general setters
    void setOFStream(ofstream& FileMESSAGE);
    void setOSS(ostream& oss);
};

// for queue
class _xqsub {
 public:
  // trivial constructurs/destuctors/operators
  _xqsub();                                                     // default, just allocate
  ~_xqsub();                                                    // kill everything
  _xqsub(const _xqsub& b);                                      // constructor copy
  const _xqsub& operator=(const _xqsub &b);                     // copy
  void clear(void);                                             // clear
  stringstream QSUB;                                            //
  stringstream QSUB_orig;                                       //
  bool         QSUB_generated;                                  //
  bool         QSUB_changed;                                    //
 private:                                                       //
  void free();                                                  // free space
  void copy(const _xqsub& b);                                   //
};

// for a vasp run
class _xvasp {
 public:
  // trivial constructurs/destuctors/operators
  _xvasp();                                                     // default, just allocate
  ~_xvasp();                                                    // kill everything
  _xvasp(const _xvasp& b);                                      // constructor copy
  const _xvasp& operator=(const _xvasp &b);                     // copy
  void clear(void);                                             // clear
  // aflow_xatom.cpp contains the code
  // CONTENT
  xstructure   str;
  string       Directory;
  string       AnalyzeLabel;
  _xqsub       xqsub;
  xoption aopts;
  // --------------------------------
  // VASP INPUT CONTENT
  // [OBSOLETE] bool         AFLOWIN_FLAG::VASP;
  stringstream POSCAR;
  stringstream POSCAR_orig;
  // [OBSOLETE] bool         POSCAR_generated;
  // [OBSOLETE] bool         POSCAR_changed;
  uint         POSCAR_index;
  stringstream INCAR;
  stringstream INCAR_orig;
  // [OBSOLETE] bool         INCAR_generated;
  // [OBSOLETE] bool         INCAR_changed;
  stringstream KPOINTS;
  stringstream KPOINTS_orig;
  // [OBSOLETE] bool         KPOINTS_generated;
  // [OBSOLETE] bool         KPOINTS_changed;
  stringstream POTCAR;
  stringstream POTCAR_orig;
  // [OBSOLETE] bool         POTCAR_generated;
  // [OBSOLETE] bool         POTCAR_changed;
  string       POTCAR_TYPE;
  bool         POTCAR_TYPE_DATE_PRINT_flag;  // to add Zr:POT_PAW:01Apr2001 to directory...
  double       POTCAR_ENMAX;
  double       POTCAR_ENMIN;
  bool         POTCAR_PAW;
  stringstream POTCAR_POTENTIALS;
  // --------------------------------
  // VASP INPUT CONTENT
  stringstream OUTCAR;  // OUTPUT
  stringstream CONTCAR; // OUTPUT
  stringstream OSZICAR; // OUTPUT
  // --------------------------------
  // QE INPUT CONTENT
  // [OBSOLETE] bool         AFLOWIN_FLAG::QE;                // FUTURE
  stringstream QE_GEOM;                        // FUTURE
  stringstream QE_GEOM_orig;                   // FUTURE
  bool         QE_GEOM_generated;              // FUTURE
  bool         QE_GEOM_changed;                // FUTURE
  uint         QE_GEOM_index;                  // FUTURE
  // --------------------------------
  // ABINIT INPUT CONTENT
  // [OBSOLETE] bool         AFLOWIN_FLAG::ABINIT;            // FUTURE
  // --------------------------------
  // AIMS INPUT CONTENT
  // [OBSOLETE] bool         AFLOWIN_FLAG::AIMS;            // FUTURE
  // --------------------------------
  int          NCPUS;
  int          NRELAX; // job to do         // -1 (static) 0(error) 1,2,3,4.... (relaxes)  -2(run kpoints)
  int          NRELAXING; // job doing (to monitor odd/even)
  // --------------------------------
  // content for AVASP generation
  bool   AVASP_aflowin_only_if_missing;                          //
  string AVASP_dirbase;
  string AVASP_libbase;
  string AVASP_label;
  string AVASP_parameters;
  double AVASP_volume_in;
  string AVASP_potential;
  bool   AVASP_alpha_fix;
  uint   AVASP_prototype_mode;
  bool   AVASP_prototype_from_library_;
  bool   AVASP_directory_from_library_;
  // [OBSOLETE] xoption AVASP_flag_AFLOW_WRITE;
  // [OBSOLETE] bool   AVASP_flag_AUTO_PSEUDOPOTENTIALS;
  // [OBSOLETE] bool   AVASP_flag_NBANDS;
  // [OBSOLETE] bool   AVASP_flag_POTIM;
  // [OBSOLETE] double AVASP_value_POTIM;
  // [OBSOLETE] bool   AVASP_flag_PSTRESS;
  // [OBSOLETE] double AVASP_value_PSTRESS;
  // [OBSOLETE] bool   AVASP_flag_EDIFFG;
  // [OBSOLETE] double AVASP_value_EDIFFG;  
  // [OBSOLETE] bool   AVASP_flag_SKIP_NOMIX;
  // [OBSOLETE] bool   AVASP_flag_WAVECAR;
  // [OBSOLETE] bool   AVASP_flag_CHGCAR;
  // [OBSOLETE] bool   AVASP_flag_SPIN;
  // [OBSOLETE] bool   AVASP_flag_SPIN_REMOVE_RELAX_1;
  // [OBSOLETE] bool   AVASP_flag_SPIN_REMOVE_RELAX_2;
  // [OBSOLETE] bool   AVASP_flag_BADER;
  // [OBSOLETE] bool   AVASP_flag_ELF;
  // [OBSOLETE] bool   AVASP_flag_LSCOUPLING;
  // [OBSOLETE] bool   AVASP_flag_AUTO_MAGMOM;
  // [OBSOLETE] bool   AVASP_flag_RELAX_FORCES;
  int    AVASP_value_NSW;
  // [OBSOLETE] bool   AVASP_flag_KPPRA;
  int    AVASP_value_KPPRA;
  string AVASP_KSCHEME;
  int    AVASP_value_KPPRA_STATIC;
  string AVASP_STATIC_KSCHEME;
  // [OBSOLETE] bool   AVASP_flag_PRECISION_flag;
  string AVASP_flag_PRECISION_scheme;
  // [OBSOLETE] bool   AVASP_flag_PRECISION_preserved;
  // [OBSOLETE] bool   AVASP_flag_ALGO_flag;
  string AVASP_flag_ALGO_scheme;
  // [OBSOLETE] bool   AVASP_flag_ALGO_preserved;
  // [OBSOLETE] string AVASP_flag_METAGGA_scheme;
  // [OBSOLETE] tring AVASP_flag_IVDW_scheme;
  // [OBSOLETE] bool   AVASP_flag_ABMIX_flag;
  string AVASP_flag_ABMIX_scheme;
  xoption AVASP_flag_TYPE;   // TYPE 
  // [OBSOLETE] bool   AVASP_flag_forceLDAU;
  // [OBSOLETE] bool   AVASP_flag_forceNOLDAU;  
  // [OBSOLETE] bool   AVASP_flag_LDAU1;
  // [OBSOLETE] bool   AVASP_flag_LDAU2;
  // [OBSOLETE] bool   AVASP_flag_LDAU_ADIABATIC;
  // [OBSOLETE] bool   AVASP_flag_LDAU_CUTOFF;
  string AVASP_LDAU_PARAMETERS_STRING;
  double AVASP_LDAU_PARAMETERS_UJSUM;
  // [OBSOLETE] xoption AVASP_flag_CONVERT_UNIT_CELL;   // CONVERT_UNIT_CELL
  // [OBSOLETE] bool   AVASP_flag_PRESERVE_VOLUME;
  // [OBSOLETE] bool   AVASP_flag_EXTRA_INCAR;
  stringstream AVASP_EXTRA_INCAR;
  bool   AVASP_flag_MPI;
  bool   AVASP_flag_RUN_RELAX;
  bool   AVASP_flag_RUN_RELAX_STATIC;
  bool   AVASP_flag_RUN_RELAX_STATIC_BANDS;
  bool   AVASP_flag_RUN_STATIC_BANDS;
  string AVASP_path_BANDS;
  uint   AVASP_value_BANDS_GRID;
  bool   AVASP_flag_RUN_STATIC;
  bool   AVASP_flag_GENERATE;
  // [OBSOLETE] xoption AVASP_flag_preserve;   // PRESERVE
  // [OBSOLETE] //  bool   AVASP_flag_preserve_POSCAR;
  // [OBSOLETE] //  bool   AVASP_flag_preserve_KPOINTS;
  // [OBSOLETE] //  bool   AVASP_flag_preserve_CHGCAR;
  // [OBSOLETE] //  bool   AVASP_flag_preserve_WAVECAR;
  // [OBSOLETE] //  bool   AVASP_flag_preserve_WAVEDER;
  // -------------------------------- FUNCTIONS
  double GetZVAL(void);
  double GetCellAtomZVAL(string mode);  // CELL ATOM
  double GetPOMASS(void);
  double GetCellAtomPOMASS(string mode); // CELL ATOM
 private:                                                       //
  void free();                                                  // free space
  void copy(const _xvasp& b);                                   //
};

//for an aims run
class _xaims {
 public:
  _xaims();
  ~_xaims();
  _xaims(const _xaims& b);
  const _xaims& operator=(const _xaims& b);
  void clear();
  uint          GEOM_index;
  xstructure    str;
  string        Directory;
  _xqsub        xqsub;
  xoption       aopts;
  int           NCPUS;
  // --------------------------------
  // AIMS INPUT CONTENT
  stringstream  CONTROL;
  stringstream  CONTROL_orig;
  bool          CONTROL_generated;
  bool          CONTROL_changed;
  string        CONTROL_FILE_NAME;
  stringstream  GEOM;
  stringstream  GEOM_orig;
  bool          GEOM_generated;
  bool          GEOM_changed;
  string        GEOM_FILE_NAME;
  string        OUTPUT_FILE_NAME;
 private:                                                       //
  void free();                                                  // free space
  void copy(const _xaims& b);                                   //
};

// for a alien run
class _xalien {
 public:
  // trivial constructurs/destuctors/operators
  _xalien();                                                    // default, just allocate
  ~_xalien();                                                   // kill everything
  _xalien(const _xalien& b);                                    // constructor copy
  const _xalien& operator=(const _xalien &b);                   // copy
  void clear(void);                                             // clear
  // aflow_xatom.cpp contains the code
  // CONTENT
  string       Directory;
  _xqsub       xqsub;
  stringstream INPUT;
  stringstream INPUT_orig;
  bool         INPUT_generated;
  bool         INPUT_changed;
  string       INPUT_FILE_NAME;
  string       OUTPUT_FILE_NAME;
  // ----------------
  int          NCPUS;
  int          NRELAX;          // -1 (static) 0(error) 1,2,3,4.... (relaxes)  -2(run kpoints)
 private:                                                       //
  void free();                                                  // free space
  void copy(const _xalien& b);                                  //
};

// for a generic run
class _xinput {
  public:
    _xinput();
    _xinput(_xvasp& xvasp);
    _xinput(_xaims& xaims);
    _xinput(_xalien& xalien);
    ~_xinput();
    _xinput(const _xinput& b);
    const _xinput& operator=(const _xinput &b);
    void clear();
    bool AFLOW_MODE_VASP;
    _xvasp xvasp;
    bool AFLOW_MODE_AIMS;
    _xaims xaims;
    bool AFLOW_MODE_ALIEN;
    _xalien xalien;
    void setXVASP(_xvasp& xvasp);
    void setXAIMS(_xaims& xaims);
    void setXALIEN(_xalien& xalien);
    xstructure& getXStr();
    string& getDirectory();
    void setXStr(const xstructure& str,bool set_all=false);
    void setDirectory(const string Directory,bool set_all=false);
  private:
    void free();
    void copy(const _xinput& b);
};

/* typedef struct { */
/*   //  _xvasp   *pxvasp; */
/*   _aflags  *paflags; */
/*   // _kflags  *pkflags; */
/*   // _vflags  *pvflags; */
/*   // string   stringA; */
/*   // string   stringB; */
/*   // int      mode; */
/*   // ofstream *pFileMESSAGE; */
/*   // bool     *pQUIET; */
/*   char  ***pargv; */
/*   bool    *pbusy; */
/* } _threaded_KBIN_params; */

//void xstructure::free();
//void xstructure::copy(const xstructure& b);
//xstructure::xstructure(string structure_title);
//xstructure::xstructure(const xstructure& b);
//xstructure::~xstructure();
//const xstructure& xstructure::operator=(const xstructure& b);
xstructure GetStructure(const int& iomode,ifstream& input);     // plug from cin
xstructure GetStructure(const int& iomode,const string& Directory); // plug from a directory
//void xstructure::SetCoordinates(const int& mode);
xstructure SetSDNumbers(const xstructure& a,const vector<string>& in_sd);
xstructure SetSDTypes(const xstructure& a,const vector<string>& in_sd);
vector<int> GetTypes(const xstructure& a);
vector<string> GetNames(const xstructure& a);
vector<string> GetCleanNames(const xstructure& a);
vector<double> GetSpins(const xstructure& a);
string GetElementName(string stringin);
string GetSpaceGroupName(int spacegroupnumber, string directory=""); //DX 20180526 - add directory
string GetSpaceGroupLabel(int spacegroupnumber);
string GetSpaceGroupSchoenflies(int spacegroupnumber, string directory=""); // DX 9/1/17 //DX 20180526 - add directory
string GetSpaceGroupHall(int spacegroupnumber, int setting=1, string directory=""); // DX 9/1/17 //DX 20180526 - add directory //DX 20180806 - added setting
string GetLaueLabel(string& point_group); // DX 9/1/17 //DX 20180526 - add directory

#define RADIANTS 0
#define DEGREES  1
#define _calculate_symmetry_default_sgroup_radius_   2.0

xmatrix<double> MetricTensor(const xstructure& a); // CO 180409
xmatrix<double> MetricTensor(const xmatrix<double>& lattice,double scale=1.0); // CO 180409
xmatrix<double> ReciprocalLattice(const xstructure& a); // CO 180409
xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice,double scale=1.0); // CO 180409
string KPPRA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const int& NK);
string KPPRA(xstructure& str,const int& _NK);
string KPPRA_DELTA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const double& DK);
string KPPRA_DELTA(xstructure& str,const double& DK);
int GetNBANDS(int electrons,int nions,int spineach,bool ispin);
double GetZVAL(const stringstream& sss,vector<double>& vZVAL);
double GetZVAL(const _xvasp& xvasp,vector<double>& vZVAL);
double GetZVAL(const string& directory,vector<double>& vZVAL);
double GetCellAtomZVAL(const stringstream& sss,vector<double>& vZVAL,const stringstream& sstr,vector<double>& sZVAL,string mode);  // sss sstr returns ZVAL cell, VAL and sZVAL
double GetCellAtomZVAL(const string& directory,vector<double>& vZVAL,vector<double>& sZVAL,string mode);  // from directory POT/POS returns total ZVAL cell, vZVAL and sZVAL
double GetPOMASS(const stringstream& sss,vector<double>& vPOMASS);
double GetPOMASS(const _xvasp& xvasp,vector<double>& vPOMASS);
double GetPOMASS(const string& directory,vector<double>& vPOMASS);
double GetCellAtomPOMASS(const stringstream& sss,vector<double>& vPOMASS,const stringstream& sstr,vector<double>& sPOMASS,string mode);  // sss sstr returns POMASS cell, VAL and sPOMASS
double GetCellAtomPOMASS(const string& directory,vector<double>& vPOMASS,vector<double>& sPOMASS,string mode);  // from directory POT/POS returns total POMASS cell, vPOMASS and sPOMASS
double GetVol(const xmatrix<double>& lat);
double det(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
double GetVol(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
double det(const double&,const double&,const double&,const double&,const double&,const double&,const double&,const double&,const double&);
//double getcos(const xvector<double>& a,const xvector<double>& b);  // removed and put in aurostd_xvector.h as cos(xvector,xvector) and sin(xvector,xvector)
xvector<double> Getabc_angles(const xmatrix<double>& lat,const int& mode);
xvector<long double> Getabc_angles(const xmatrix<long double>& lat,const int& mode);
xvector<double> Getabc_angles(const xmatrix<double>& lat,const xvector<int>& permut,const int& mode);
xvector<double> Getabc_angles(const xvector<double>& r1,const xvector<double>& r2,const xvector<double>& r3,const int& mode);
xvector<double> Getabc_angles(const xvector<double>& r1,const xvector<double>& r2,const xvector<double>& r3,const xvector<int>& permut,const int& mode);
#define _Getabc_angles Getabc_angles
//#define _Getabc_angles __NO_USE_Sortabc_angles
xvector<double> Sortabc_angles(const xmatrix<double>& lat,const int& mode);
xmatrix<double> GetClat(const xvector<double>& abc_angles);
xmatrix<double> GetClat(const double &a,const double &b,const double &c,const double &alpha,const double &beta,const double &gamma);
xstructure GetIntpolStr(xstructure strA,xstructure strB,const double& f,const string& path_flag);
double RadiusSphereLattice(const xmatrix<double>& lattice,double scale=1.0); // CO 180409
xvector<int> LatticeDimensionSphere(const xmatrix<double>& lattice,double radius,double scale=1.0); // CO 180409
xvector<int> LatticeDimensionSphere(const xstructure& str,double radius);
xvector<double> F2C(const double& scale,const xmatrix<double>& lattice,const xvector<double>& fpos);    // fpos are F components per COLUMS !
xvector<double> F2C(const xmatrix<double>& lattice,const xvector<double>& fpos);                        // fpos are F components per COLUMS !
xvector<double> C2F(const double& scale,const xmatrix<double>& lattice,const xvector<double>& cpos);    // cpos are C components per COLUMS !
xvector<double> C2F(const xmatrix<double>& lattice,const xvector<double>& cpos);                        // cpos are C components per COLUMS !
xmatrix<double> F2C(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fpos);    // fpos are F components per COLUMS !
xmatrix<double> F2C(const xmatrix<double>& lattice,const xmatrix<double>& fpos);                        // fpos are F components per COLUMS !
xmatrix<double> C2F(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cpos);    // cpos are C components per COLUMS !
xmatrix<double> C2F(const xmatrix<double>& lattice,const xmatrix<double>& cpos);                        // cpos are C components per COLUMS !
_atom F2C(const double& scale,const xmatrix<double>& lattice,const _atom& iatom);                       // atom.fpos are F components per COLUMS !
_atom F2C(const xstructure& str,const _atom& iatom);                                                    // atom.fpos are F components per COLUMS !
_atom C2F(const double& scale,const xmatrix<double>& lattice,const _atom& iatom);                       // atom.cpos are C components per COLUMS !
_atom C2F(const xmatrix<double>& lattice,const _atom& iatom);                                           // atom.cpos are C components per COLUMS !
_atom F2C(const double& scale,const xstructure& str,const _atom& iatom);                                // atom.fpos are F components per COLUMS !
_atom F2C(const xstructure& str,const _atom& iatom);                                                    // atom.fpos are F components per COLUMS !
_atom C2F(const double& scale,const xstructure& str,const _atom& iatom);                                // atom.fpos are F components per COLUMS !
_atom C2F(const xstructure& str,const _atom& iatom);                                                    // atom.cpos are C components per COLUMS !
xmatrix<double> FF2CC(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fmat);  // fmat is an operation in F coordinates
xmatrix<double> FF2CC(const xmatrix<double>& lattice,const xmatrix<double>& fmat);                      // fmat is an operation in F coordinates
xmatrix<double> CC2FF(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cmat);  // cmat is an operation in C coordinates
xmatrix<double> CC2FF(const xmatrix<double>& lattice,const xmatrix<double>& cmat);                      // cmat is an operation in C coordinates
// DX and CO - START
double BringInCell(const double& x);
double BringInCell_20161115(const double& x);
double BringInCell_20160101(const double& x);
double BringInCell(const double& x);
double BringInCell_20160101(const double& x);
xvector<double> BringInCell(const xvector<double>& v_in,double epsilon);
xvector<double> BringInCell_20161115(const xvector<double>& v_in,double epsilon);
xvector<double> BringInCell_20160101(const xvector<double>& v_in,double epsilon);
xvector<double> BringInCell(const xvector<double>& v_in);
xvector<double> BringInCell2(const xvector<double>& v_in);
xvector<double> BringInCell2_20161115(const xvector<double>& v_in);
xvector<double> BringInCell2_20160101(const xvector<double>& v_in, double tolerance);
// DX and CO - END
xstructure IdenticalAtoms(const xstructure& a);                                // Make identical atoms
//xstructure SwapSpecies(const xstructure& a,const uint& A,const uint& B);       // Permute Species A with B (safe for species C).
//xstructure SwapCoordinates(const xstructure& str,const uint& i,const uint& j); // Permute Coordinates i with j
//string SpeciesLabel(const xstructure& a,const uint& A);                        // Returns the Label of the specie A (if available)
//string SpeciesString(const xstructure& a);                                           // Gives a string with the list of all the species
bool GetNiggliCell(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q);
bool GetNiggliCell_20180213(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q); // DX 2/13/18 - new dated function
bool GetNiggliCell_20180101(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q); // DX 2/13/18 - old dated function
// standard lattice reduction and type
string GetLatticeType(xmatrix<double> lattice);
string GetLatticeType(xvector<double> data);
xstructure Standard_Primitive_UnitCellForm(const xstructure& a);
xstructure GetStandardPrimitive(const xstructure& a);
xmatrix<double> GetStandardPrimitive(xmatrix<double> lattice);
xvector<double> GetStandardPrimitive(xvector<double> data);
xstructure Standard_Conventional_UnitCellForm(const xstructure& a);
xstructure GetStandardConventional(const xstructure& a);
xmatrix<double> GetStandardConventional(xmatrix<double> lattice);
xvector<double> GetStandardConventional(xvector<double> data);
// niggli
xstructure GetNiggliStr(const xstructure& in_str);
xmatrix<double> GetNiggliStr(const xmatrix<double>& lattice);
xstructure NiggliUnitCellForm(const xstructure& a);
xmatrix<double> NiggliUnitCellForm(const xmatrix<double>& lattice);
// minkowsky
xstructure MinkowskiBasisReduction(const xstructure& a);
xmatrix<double> MinkowskiBasisReduction(const xmatrix<double>& lattice);
// optimal lattice reduction
xstructure LatticeReduction(const xstructure& a);
xmatrix<double> LatticeReduction(const xmatrix<double>& lattice);
// CO 170807 - START
deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew);
deque<_atom> foldAtomsInCell(xstructure& a, xmatrix<double>& lattice_new, bool& skew, double& tol, bool fold_in_only=false);
deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew, double& tol);
xstructure GetPrimitiveVASP(xstructure& a);
xstructure GetPrimitiveVASP(xstructure& a,double tol);
// CO 170807 - STOP
// bring cell in,compact, wigner seitz
_atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon);
// DX and CO - START
_atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon); // DX
_atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon); // DX
_atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice);
_atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice); // DX
_atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice); // DX
xstructure BringInCell(const xstructure& a,double epsilon);
xstructure BringInCell_20161115(const xstructure& a,double epsilon); // DX
xstructure BringInCell_20160101(const xstructure& a,double epsilon); // DX
xstructure BringInCell(const xstructure& a);
xstructure BringInCell_20161115(const xstructure& a); // DX
xstructure BringInCell_20160101(const xstructure& a); // DX
// DX and CO - END
xstructure BringInCompact(const xstructure& a);
xstructure BringInWignerSeitz(const xstructure& a);
// primitive stuff
xstructure GetPrimitive(const xstructure& a);
xstructure GetPrimitive(const xstructure& a,double tol);
xstructure GetPrimitive1(const xstructure& a);
xstructure GetPrimitive2(const xstructure& a);
xstructure GetPrimitive3(const xstructure& a);
bool IsTranslationFVector(const xstructure& a,const xvector<double>& ftvec);
bool IsTranslationCVector(const xstructure& a,const xvector<double>& ctvec);
// other eggs
xstructure ReScale(const xstructure& a,const double& in_scale);
xstructure SetScale(const xstructure& a,const double& in_scale);
xstructure SetVolume(const xstructure& a,const double& in_volume);
xstructure InflateLattice(const xstructure& a,const double& coefficient);
xstructure InflateVolume(const xstructure& a,const double& coefficient);
double GetVolume(const xstructure& a);
double Volume(const xstructure& a);
//DX 20180726 - START
_atom BringCloseToOrigin(_atom& atom, xmatrix<double>& f2c);
bool uniqueAtomInCell(_atom& atom, deque<_atom>& atoms);
bool alreadyInCell(_atom& atom, deque<_atom> atoms);
//DX 20180726 - END
// DX and CO - START
bool inCell(xvector<double>& pos_vec); // DX
// DX and CO - END
xstructure GetSuperCell(const xstructure& a,const xmatrix<double>& sc);
xstructure GetSuperCell(const xstructure& a,const xvector<double>& sc);
xstructure GetSuperCell(const xstructure& a,const xvector<int>& sc);
xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33);
xstructure GetSuperCell(const xstructure& a,const int& sc1,const int& sc2,const int& sc3);
//corey START
xstructure GetSuperCell(const xstructure& a,const xmatrix<double>& sc,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis);
xstructure GetSuperCell(const xstructure& a,const xvector<double>& sc,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis);
xstructure GetSuperCell(const xstructure& a,const xvector<int>& sc,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis);
xstructure GetSuperCell(const xstructure& a,const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis);
xstructure GetSuperCell(const xstructure& a,const int& sc1,const int& sc2,const int& sc3,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis);
//corey END
bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose,double radius);
bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,double radius);
bool CalculateSymmetry(xstructure& str,bool ossverbose,double radius);
bool CalculateSymmetry(xstructure& str,double radius);
bool CalculateSymmetry(xstructure& str,bool ossverbose);
bool CalculateSymmetry(xstructure& str);
void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose);
void CalculateSymmetryPointGroup(xstructure& str);
void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose);
void CalculateSymmetryPointGroupCrystal(xstructure& str);
void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose);
void CalculateSymmetryFactorGroup(xstructure& str);
void CalculateSymmetryPointGroupKlattice(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryPointGroupKlattice(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryPointGroupKlattice(xstructure& str,bool ossverbose);
void CalculateSymmetryPointGroupKlattice(xstructure& str);
xstructure Rotate(const xstructure& a,const xmatrix<double>& rm);
xstructure GetLTCell(const xmatrix<double>& lt,const xstructure& str);
xstructure GetLTFVCell(const xvector<double>& nvec,const double phi,const xstructure& str);
xstructure ShiftPos(const xstructure& a,const xvector<double>& shift,const int& flag);
xstructure ShiftCPos(const xstructure& a,const xvector<double>& shift);
xstructure ShiftFPos(const xstructure& a,const xvector<double>& shift);
double MaxStructureLattice(const xstructure& str);
double MinStructureLattice(const xstructure& str);
double AtomDist(const xstructure& str,const _atom& atom1,const _atom& atom2);
bool SameAtom(const xstructure& str,const _atom& atom1,const _atom& atom2);
bool SameAtom(const _atom& atom1,const _atom& atom2);
bool DifferentAtom(const xstructure& str,const _atom& atom1,const _atom& atom2);
vector<double> GetNBONDXX(const xstructure& a);
int GenerateGridAtoms(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2);
int GenerateGridAtoms(xstructure& str,int d1,int d2,int d3);
int GenerateGridAtoms(xstructure& str,int d);
int GenerateGridAtoms(xstructure& str,const xvector<int>& dims);
int GenerateGridAtoms(xstructure& str);
int GenerateGridAtoms(xstructure& str,const double& radius);

void l2ijk(const xstructure& str,const int &l,int &i,int &j,int &k);
void l2ijk(const xstructure& str,const int &l,xvector<int>& ijk);
xvector<int> l2ijk(const xstructure& str,const int &l);
void ijk2l(const xstructure& str,int &l,const int &i,const int &j,const int &k);
void ijk2l(const xstructure& str,int &l,const xvector<int>& ijk);
int ijk2l(const xstructure& str,const int &i,const int &j,const int &k);
int ijk2l(const xstructure& str,const xvector<int>& ijk);
xvector<double> r_lattice(const xstructure& str,const int &l);
xvector<double> r_lattice(const xstructure& str,const int &i,const int &j,const int &k);
xvector<double> r_lattice(const xstructure& str,const xvector<int>& ijk);
xstructure AQEgeom2aims(istream& input);
xstructure AQEgeom2abinit(istream& input);
xstructure AQEgeom2qe(istream& input);
xstructure AQEgeom2vasp(istream& input);

// ----------------------------------------------------------------------------
// Structure Prototypes
// aflow_xproto.cpp
#define _HTQC_PROJECT_STRING_ "HTQC Project"
#define _TERNARY_PROJECT_STRING_ "HTQC^3 Project"
#define _ICSD_STRING_ "(icsd library)"
#define _ICSD_PROJECT_STRING_ "ICSD Project"
#define _ICSD_AFLOWLIB_STRING_ "(icsd_aflowlib library)"

// aflow_xproto.cpp
namespace aflowlib {
  string PrototypeCleanLatticeString(const string& latticeIN);
}
double NearestNeighbour(const xstructure &str_in);

// for HTQC
#define STRUCTURE_MODE_NONE             0
#define STRUCTURE_MODE_RAW              1
#define STRUCTURE_MODE_ABC              2
#define STRUCTURE_MODE_WYC              3
#define STRUCTURE_MODE_ICSD             4
#define STRUCTURE_MODE_HTQC_ICSD        5
#define STRUCTURE_MODE_USE              6
#define STRUCTURE_MODE_REMOVE           7
#define STRUCTURE_MODE_SPECIES          8
#define STRUCTURE_MODE_SWAP_AB          9
#define STRUCTURE_MODE_SWAP_BC         10
#define STRUCTURE_MODE_SWAP_AC         11
#define STRUCTURE_MODE_SWAP_XY         12
#define STRUCTURE_MODE_PRIM            13
#define STRUCTURE_MODE_CONVENTIONAL    14
#define STRUCTURE_MODE_VOLUME          15
#define LIBRARY_MODE_ICSD               0
#define LIBRARY_MODE_ICSD_AFLOWLIB      1
#define LIBRARY_MODE_HTQC               2
#define LIBRARY_MODE_HTQC_ICSD          3
#define LIBRARY_MODE_HTQC_ICSD_AFLOWLIB 4
#define LIBRARY_MODE_LIB3               5
#define LIBRARY_MODE_LIB4               6
#define LIBRARY_MODE_LIB5               7
#define LIBRARY_MODE_LIB6               8
#define LIBRARY_MODE_LIB7               9
#define LIBRARY_MODE_LIB8               10
#define LIBRARY_MODE_LIB9               11
#define LIBRARY_MODE_PROTOTYPE          12
#define LIBRARY_MODE_XSTRUCTURE         13
string* LOAD_Library_ICSD(string file);

namespace aflowlib {
  struct _PROTO_PARAMS{
    string label;
    string parameters;
    deque<string> vatomX;
    deque<double> vvolumeX;
    double volume_in;
    int mode;
    bool flip_option;
  };

  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters,string atA,double volA);
  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters,string atA);
  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters);
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters,string atA,double volA);
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters,string atA);
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters);
  uint PrototypeLibrariesSpeciesNumber(const string& label);
  // xstructure PrototypeLibraries(ostream &oss,string label,string parameters,int mode=LIBRARY_MODE_HTQC);
  // xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,int mode=LIBRARY_MODE_HTQC);
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,int mode);
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,int mode);
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode);//=LIBRARY_MODE_HTQC);
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option);
  xstructure PrototypeLibraries(ostream &oss,_PROTO_PARAMS *PARAMS);

  string PrototypesHelp(void);
  string PrototypesIcsdHelp(string options);
  string CALCULATED(string options);
  string CALCULATED_ICSD_RANDOM(void);
  // aflow_xproto_gus.cpp
  xstructure PrototypeBinaryGUS(ostream &FileMESSAGE,string label);
  xstructure PrototypeBinaryGUS(ostream &FileMESSAGE,string label,string atA,string atB);
  xstructure PrototypeBinaryGUS(ostream &FileMESSAGE,string label,string atA,double volA,string atB,double volB,double vol_in);
}

extern string PrototypeBinaryGUS_Cache_Library[];

// ----------------------------------------------------------------------------
// aflow_anrl.cpp
//DX 20180710 - updated - #define DOI_ANRL " [ANRL doi: arXiv:1607.02532]"
#define DOI_ANRL " [ANRL doi: 10.1016/j.commatsci.2017.01.017 (part 1), arXiv:1806.07864 (part 2)]" //DX 20180710 - updated
#define DOI_POCC " [POCC doi: 10.1021/acs.chemmater.6b01449]"

namespace anrl {
  xstructure PrototypeANRL(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option);
  uint PrototypeANRL_LoadList(vector<string>& vproto,
			      vector<string>& vproto_label,
			      vector<uint>& vproto_nspecies,
			      vector<uint>& vproto_natoms,
			      vector<uint>& vproto_spacegroup,
			      vector<uint>& vproto_nunderscores,
			      vector<uint>& vproto_nparameters,
			      vector<string>& vproto_Pearson_symbol,
			      vector<string>& vproto_params,
			      vector<string>& vproto_Strukturbericht,
			      vector<string>& vproto_prototype,
			      vector<string>& vproto_dialect);
  bool vproto2tokens(string proto,
		     string& label,
		     uint& nspecies,
		     uint& natoms,
		     uint& spacegroup,
		     uint& nunderscores,
		     uint& nparameters,
		     string& Pearson_symbol,
		     string& params,
		     string& Strukturbericht,
		     string& prototype,
		     string& dialect);
  bool PrototypeANRL_Consistency(ostream &oss,uint vparameters_size,uint proto_nparameters,string proto_prototype,
                                 string proto_label,string proto_Strukturbericht,string proto_Pearson_symbol,
                                 uint proto_spacegroup, string proto_params, uint print_mode); //DX 20180710 - added print_mode
  xstructure rhl2hex(xstructure& str, double& a, double& c); 
}

// ----------------------------------------------------------------------------
// Various prototypes to be moved somewhere sometime
// PROTOTYPES
// uint argsprint(vector<string> argv);
// ----------------------------------------------------------------------------
// aflow.cpp
string aflow_get_time_string(void);
string aflow_get_time_string_short(void);
string strPID(void);
int AFLOW_main(vector<string> &argv);
namespace aflow {
  string License_Preamble_aflow(void);
  string Intro_aflow(string x);
  string Intro_sflow(string x);
  string Intro_HELP(string x);
  string Banner(string type);
}
int VASP_Main(vector<string> argv);
int GRND_Main(vector<string> argv);
namespace KBIN {
  int KBIN_Main(vector<string> argv);
}
string MessageTime(void);
string MessageHostTime(const _aflags& aflags);
string MessageDir(const _aflags& aflags);
string MessageDirTime(const _aflags& aflags);
string MessageDirHostTime(const _aflags& aflags);
bool AFLOW_BlackList(string hostname);

// ----------------------------------------------------------------------------
// aflow_pthreads.cpp
namespace AFLOW_PTHREADS {
  int GetTotalCPUs(void);
  bool Check_Threads(vector<string> argv,const bool& VERBOSE);
  void Clean_Threads(void);
  void No_Threads(void);
  bool Available_Free_Threads(int &fthread);
  bool Wait_Available_Free_Threads(int &fthread,const double& pthread_wait,const bool& VERBOSE);
  bool Wait_Available_Free_Threads(int &fthread,const bool& VERBOSE);
}
// interfaces
namespace KBIN {
  void RUN_Directory_PTHREADS(_aflags &aflags);
  void *_threaded_interface_RUN_Directory(void *ptr);
} // namespace KBIN
namespace aurostd { // Multithreaded add on to aurostd
  bool multithread_execute(deque<string> vcommand,int NUM_THREADS,bool VERBOSE);
  bool multithread_execute(deque<string> vcommand,int NUM_THREADS);
  bool multithread_execute(deque<string> vcommand);
  bool multithread_execute(vector<string> vcommand,int NUM_THREADS,bool VERBOSE);
  bool multithread_execute(vector<string> vcommand,int NUM_THREADS);
  bool multithread_execute(vector<string> vcommand);
} // namespace aurostd
namespace AFLOW_PTHREADS {
  bool MULTI_sh(vector<string> argv);
  bool MULTI_compress(string cmd,vector<string> argv);
  bool MULTI_zip(vector<string> argv);
  bool MULTI_bz2xz(vector<string> argv);bool MULTI_xz2bz2(vector<string> argv);
  bool MULTI_gz2xz(vector<string> argv);
}
namespace sflow {
  void KILL(string options);
  void JUST(string options,istream& input,string mode);
  void QSUB(string options);
  void QSUB(string options,string cmd);
  void QDEL(string options);
  void QDEL(string options,string cmd);
}

// ----------------------------------------------------------------------------
// aflow_kbin.cpp
//int KbinCheckInputFiles(string Directory,ofstream& FileERROR);
namespace KBIN {
  void MPI_Extract(string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  void RUN_Directory(_aflags& aflags);
  void AFLOW_RUN_Directory(const _aflags& aflags);
  void RUN_DirectoryScript(const _aflags& aflags,const string& script,const string& output);
  void CompressDirectory(const _aflags& aflags,const _kflags& kflags);
  void CompressDirectory(const _aflags& aflags);
  void Clean(const _aflags& aflags);
  void Clean(const string directory);
  void XClean(string options);
  void GenerateAflowinFromVASPDirectory(_aflags& aflags);
  void StartStopCheck(const string &AflowIn,string str1,string str2,bool &flag,bool &flagS);
  void StartStopCheck(const string &AflowIn,string str1,bool &flag,bool &flagS);
  bool Legitimate_aflowin(string aflowindir,const bool& osswrite,ostringstream& oss);
  bool Legitimate_aflowin(string aflowindir);
}

// ----------------------------------------------------------------------------
// aflow_qsub.cpp
namespace KBIN {
  bool QSUB_Extract(_xqsub& xqsub,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool QSUB_RunFinished(_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  void QSUB_WaitFinished(_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  bool QSUB_Extract_Mode1(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool QSUB_Extract_Mode2(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool QSUB_Extract_Mode3(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
}

// ----------------------------------------------------------------------------
// aflow_ialien.cpp
namespace ALIEN {
  bool Produce_INPUT(_xalien& xalien,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_alienflags &alienflags);
  bool Modify_INPUT(_xalien& xalien,ofstream &FileMESSAGE,_aflags &aflags,_alienflags &alienflags);
  bool Write_INPUT(_xalien& xalien);
  bool Produce_INPUT_FILE(_xalien& xalien,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_alienflags &alienflags);
  bool Modify_INPUT_FILE(_xalien& xalien,ofstream &FileMESSAGE,_aflags &aflags,_alienflags &alienflags);
}

// ----------------------------------------------------------------------------
// aflow_kalien.cpp
namespace ALIEN {
  _alienflags Get_Alienflags_from_AflowIN(string &AflowIn);
  bool Run_Directory(ofstream& FileERROR,_aflags& aflags,_kflags& kflags);
}

// ----------------------------------------------------------------------------
// aflow_matlab.cpp aflow_matlab_funcs.cpp
bool KBIN_MATLAB_Extract(string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
bool KBIN_MATLAB_Run(_kflags &kflags);
_kflags KBIN_MATLAB_Get_Matlabflags_from_AflowIN(string &AflowIn);
bool KBIN_MATLAB_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
string MATLAB_FUNCS_param(void);
//string MATLAB_FUNCS_plotband(string DIRECTORY,string OPTION1);
string MATLAB_FUNCS_plotband(void);

// ----------------------------------------------------------------------------
// aflow_gnuplot_plotbz.cpp
//string GNUPLOT_FUNCS_plotbz(string DIRECTORY,string OPTION1);
string GNUPLOT_FUNCS_plotbz(void);

// ----------------------------------------------------------------------------
// aflow_ifrozsl.cpp

namespace KBIN {
  void VASP_RunPhonons_FROZSL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
}

namespace FROZSL {
  bool Extract_INPUT(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool Setup_frozsl_init_input(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool Already_Calculated_Input(const string& AflowIn);
  bool WGET_INPUT(ofstream &FileMESSAGE,string AflowIn,_aflags &aflags,_kflags &kflags);
  bool WGET_OUTPUT(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool input_TO_poscar(ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  string Generate_Input_file(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool File_INPUT(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool Write(string data,string directory);
  bool Delete(string data,string directory);
}
namespace FINDSYM {
  bool Write(string data,string directory);
}

// ----------------------------------------------------------------------------
// aflow_kvasp.cpp

namespace KBIN {
  _vflags VASP_Get_Vflags_from_AflowIN(const string &AflowIn,_aflags &aflags,_kflags& kflags);
  _vflags VASP_Get_Vflags_from_AflowIN(const string &AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags& kflags);
  bool VASP_Fix_Machine_Kflags_from_AflowIN(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags);
  bool VASP_Directory(ofstream& FileERROR,_aflags& aflags,_kflags& kflags);
  void VASP_BackupOriginal(_aflags aflags);
  // [OBSOLETE] G++6 not needed  void VASP_Wait(_xvasp& xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relaxA,string relaxB,bool qmwrite,ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relaxA,bool qmwrite,ofstream &FileMESSAGE);
  bool VASP_RunFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  void WaitFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  void VASP_Error(_xvasp &xvasp,string="",string="",string="");
  void VASP_Error(_xvasp &xvasp,ofstream &FileMESSAGE,string="",string="",string="");
  string VASP_Analyze(_xvasp &xvasp,bool qmwrite);
  void VASP_CompressDirectory(_xvasp xvasp,_kflags &kflags);
  void VASP_Backup(_xvasp& xvasp,bool qmwrite,string relax);
  void VASP_CONTCAR_Save(_xvasp xvasp,string relax);
  void VASP_Recycle(_xvasp xvasp,string relax);
  void VASP_Recycle(_xvasp xvasp,int relax_number);
  void VASP_RecycleExtraFile(_xvasp xvasp,string xfile,string relax);
  void VASP_RecycleExtraFile(_xvasp xvasp,string xfile,int relax_number);
  bool VASP_CheckUnconvergedOSZICAR(string dir);
  void GetStatDiel(string& outcar, xvector<double>& eigr, xvector<double>& eigi); // CAMILO
  void GetDynaDiel(string& outcar, xvector<double>& eigr, xvector<double>& eigi); // CAMILO
}

// ----------------------------------------------------------------------------
// aflow_aims_ivasp.cpp
namespace KBIN {
  bool VASP_Produce_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp=false);
  bool VASP_Modify_INPUT(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags);
  bool VASP_Produce_and_Modify_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp=false); // CO 180418
  bool VASP_Write_INPUT(_xvasp& xvasp,_vflags &vflags);
  bool VASP_Produce_INCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Modify_INCAR(_xvasp& xvasp,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Reread_INCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  bool VASP_Produce_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Produce_POSCAR(_xvasp& xvasp);
  bool VASP_Modify_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_vflags& vflags);
  bool VASP_Reread_POSCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  bool VASP_Produce_KPOINTS(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Modify_KPOINTS(_xvasp& xvasp,ofstream& FileERROR,_aflags& aflags,_vflags& vflags);
  bool VASP_Reread_KPOINTS(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  bool VASP_Find_DATA_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar);
  bool VASP_Find_FILE_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar);
  bool VASP_Produce_POTCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Modify_POTCAR(_xvasp& xvasp,ofstream& FileERROR,_aflags& aflags,_vflags& vflags);
  bool VASP_Reread_POTCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  string VASP_PseudoPotential_CleanName(const string& specieIN);
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX);
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX);
  bool VASP_SplitAlloySpecies(string alloy_in, string &specieA, string &specieB);
  bool VASP_SplitAlloySpecies(string alloy_in, string &specieA, string &specieB, string &specieC);
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB);
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB, vector<string> &speciesC);
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX);
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX, vector<double> &natomsX);
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB);
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB, string &species_ppC);
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &species_ppsA, vector<string> &species_ppsB);
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &species_ppsA, vector<string> &species_ppsB, vector<string> &species_ppsC);
  void VASP_MPI_Autotune(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_System_Auto(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,_vflags& vflags,int number); // for steps
  void XVASP_INCAR_Static_ON(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Relax_Static_ON(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Relax_Static_Bands_ON(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_RWIGS_Static(_xvasp& xvasp,_vflags& vflags,ofstream &FileMESSAGE,bool OPERATION);
  void XVASP_INCAR_Precision(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Metagga(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Ivdw(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_ABMIX(_xvasp& xvasp,_vflags& vflags);
  int XVASP_INCAR_GetNBANDS(_xvasp& xvasp,bool ispin);
  bool XVASP_INCAR_PREPARE_GENERIC(string command,_xvasp& xvasp,_vflags& vflags,string svalue,int ivalue,double dvalue,bool bvalue);
  //  bool XVASP_INCAR_PREPARE_GENERIC(string command,_xvasp& xvasp,_kflags &kflags,_vflags& vflags,string svalue,int ivalue,double dvalue,bool bvalue);
  // ALGO, ENMAX_MULTIPLY, IMIX, IALGO, TYPE, PAW_CORRECTIONS, NBANDS, PSTRESS, EDIFFG, POTIM, SPIN, LS_COUPLING, AUTO_MAGMOM, NWS, SYM, WAVECAR, CHGCAR
  void XVASP_INCAR_SPIN_REMOVE_RELAX(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,ofstream &FileMESSAGE);
  void XVASP_KPOINTS_IBZKPT_UPDATE(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,ofstream &FileMESSAGE);
  void XVASP_INCAR_LDAU_OFF(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_LDAU_ON(_xvasp& xvasp,_vflags& vflags,uint type);
  void XVASP_INCAR_LDAU_ADIABATIC(_xvasp& xvasp,int step);
  void XVASP_INCAR_LDAU_CUTOFF(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_KPOINTS_Dielectric_SET(_xvasp& xvasp,_kflags &kflags,_vflags& vflags,string mode_dielectric);
  void XVASP_INCAR_REMOVE_ENTRY(_xvasp& xvasp,string ENTRY,string COMMENT,bool VERBOSE);
  
  bool XVASP_KPOINTS_KPOINTS(_xvasp &xvasp,ofstream &FileMESSAGE,bool VERBOSE);
  bool XVASP_KPOINTS_KPOINTS(_xvasp &xvasp);
  // bool XVASP_KPOINTS_EVEN(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_ODD(_xvasp& xvasp); TO REMOVE
  bool XVASP_KPOINTS_OPERATION(_xvasp& xvasp,string operation);
  // bool XVASP_KPOINTS_Kshift_Gamma_EVEN(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_Kshift_Gamma_ODD(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_Kscheme(_xvasp& xvasp,string kscheme);
  bool XVASP_KPOINTS_Fix_KPPRA(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE);
  bool XVASP_KPOINTS_Fix_KSHIFT(_xvasp &xvasp,_xvasp &rxvasp,bool KAUTOSHIFT,bool VERBOSE);
  bool XVASP_KPOINTS_Fix_KPOINTS(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE);
  void XVASP_string2numbers(_xvasp& xvasp);
  void XVASP_numbers2string(_xvasp& xvasp);
  void XVASP_Afix_Clean(_xvasp& xvasp,string preserve_name);
  void XVASP_Afix_ROTMAT(_xvasp& xvasp,int mode,bool verbose,_aflags &aflags,ofstream &FileMESSAGE);
  void XVASP_Afix_NBANDS(_xvasp& xvasp,int& nbands,bool VERBOSE);
  void XVASP_Afix_POTIM(_xvasp& xvasp,double& potim,bool VERBOSE);
  double XVASP_Afix_GENERIC(string mode,_xvasp& xvasp,_kflags& kflags,_vflags& vflags,double=0.0,int=0);

  string ExtractSystemName(string directory);
  double ExtractEfermiOUTCAR(string directory);
  xstructure GetMostRelaxedStructure(string directory); //CO 180627
  vector<string> ExtractAtomicSpecies(string directory);

}

// ----------------------------------------------------------------------------
// aflow_avasp.cpp
#define _AVASP_PSEUDOPOTENTIAL_AUTO_ string("AUTO")
#define _AVASP_PSEUDOPOTENTIAL_DELIMITER_ string(":")
#define _AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_ string("COMPLETE")

struct _AVASP_PROTO{
  vector<string> ucell;
  deque<int> vkppra;
  vector<double> vpressure;
  aurostd::xoption vparams;
};

bool AVASP_MakePrototype_AFLOWIN(_AVASP_PROTO *PARAMS);
bool AVASP_MakePrototypeICSD_AFLOWIN(_AVASP_PROTO *PARAMS,bool flag_AFLOW_IN_ONLY_IF_MISSING);
void AVASP_Get_LDAU_Parameters(string species,bool &LDAU,vector<string>& vLDAUspecies,
			       vector<uint>& vLDAUtype,vector<int>& vLDAUL, vector<double>& vLDAUU, vector<double> &vLDAUJ);
string AVASP_Get_PseudoPotential_PAW_PBE_KIN(string species);
string AVASP_Get_PseudoPotential_PAW_PBE(string species);
string AVASP_Get_PseudoPotential_PAW_GGA(string species);
string AVASP_Get_PseudoPotential_PAW_LDA_KIN(string species);
string AVASP_Get_PseudoPotential_PAW_LDA(string species);
string AVASP_Get_PseudoPotential_PBE(string species);
string AVASP_Get_PseudoPotential_GGA(string species);
string AVASP_Get_PseudoPotential_LDA(string species);
bool AVASP_populateXVASP(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,_xvasp& xvasp);
bool AVASP_MakeSingleAFLOWIN(_xvasp& xvaspin,stringstream &_aflowin,bool flag_WRITE,int=-1,bool flag_PRINT=TRUE);   // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN(_xvasp& xvasp_in,bool flag_WRITE,int=-1,bool flag_PRINT=TRUE);  // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN(_xvasp& xvasp_in,int=-1,bool flag_PRINT=TRUE);  // last is pthread number, if <0 then serial
bool AVASP_DefaultValuesBinary_AFLOWIN(_xvasp &xvasp);
bool AVASP_MakeSinglePOSCAR(_xvasp& xvaspin);
bool Alloys_LibraryU(vector<string> &alloy,vector<string> &pseudosA,vector<string> &pseudosB);
bool Alloys_LibraryG(vector<string> &alloy,vector<string> &pseudosA,vector<string> &pseudosB);
bool Alloys_LibraryX(vector<string> &alloy,vector<string> &pseudosA,vector<string> &pseudosB);
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// aflow_ovasp.cpp
class xOUTCAR;
class xDOSCAR;
class xEIGENVAL;
class xPOTCAR;
class xVASPRUNXML;
class xIBZKPT;
class xKPOINTS;
class xCHGCAR;
class xVASPOUT;
namespace aflowlib { class _aflowlib_entry;}

// -------------------------------------------------------------------------------------------------
class xOUTCAR {
 public:
  xOUTCAR();                                                    // default, just allocate
  ~xOUTCAR();                                                   // kill everything
  xOUTCAR(const string& fileIN,bool=TRUE);                      // constructor from filename, QUIET
  xOUTCAR(const xOUTCAR& b);                                    // constructor copy
  const xOUTCAR& operator=(const xOUTCAR &b);                   // copy
  void clear(void);                                             // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;       // the content, and lines of it
  string SYSTEM;
  int NIONS;
  double Efermi;
  bool isLSCOUPLING;
  double natoms;                                                // for aflowlib_libraries.cpp
  double energy_cell,energy_atom;                               // for aflowlib_libraries.cpp
  double enthalpy_cell,enthalpy_atom;                           // for aflowlib_libraries.cpp
  double eentropy_cell,eentropy_atom;                           // for aflowlib_libraries.cpp
  double PV_cell,PV_atom;                                       // for aflowlib_libraries.cpp
  xmatrix<double> stress;                                       // for aflowlib_libraries.cpp
  double mag_cell,mag_atom;                                     // for aflowlib_libraries.cpp
  vector<double> vmag;                                          // for aflowlib_libraries.cpp
  vector<xvector<double> > vmag_noncoll;                        // DX 12/5/17 - non-collinear
  double volume_cell,volume_atom;                               // for aflowlib_libraries.cpp
  double pressure;                                              // for aflowlib_libraries.cpp // SAME AS PSTRESS
  double pressure_residual;                                     // for aflowlib_libraries.cpp
  double Pulay_stress;                                          // for aflowlib_libraries.cpp
  vector<aurostd::xvector<double> > vforces;                    // for aflowlib_libraries.cpp
  vector<aurostd::xvector<double> > vpositions_cartesian;       // for aflowlib_libraries.cpp
  double ENCUT,EDIFF,EDIFFG,POTIM,TEIN,TEBEG,TEEND,SMASS,NPACO,APACO,PSTRESS;     // 
  int NBANDS,NKPTS,NSW,NBLOCK,KBLOCK,IBRION,NFREE,ISIF,IWAVPR,ISYM,ISPIN;   // for aflowlib_libraries.cpp
  double total_energy_change;                                   // for aflowlib_libraries.cpp
  // DOS related values
  double EMIN,EMAX,SIGMA;                                       // eV - energy-range for DOS
  int ISMEAR;                                                   // broadening in eV -4-tet -1-fermi 0-gaus
  //  Electronic relaxation
  int IALGO;              //  algorithm                         // for aflowlib_libraries.cpp
  string LDIAG;           //   sub-space diagonalisation        // for aflowlib_libraries.cpp
  int IMIX,INIMIX,MIXPRE; //     mixing-type and parameters     // for aflowlib_libraries.cpp
  double AMIX,BMIX,AMIX_MAG,BMIX_MAG,AMIN,WC; // parameters     // for aflowlib_libraries.cpp
  // Intra band minimization
  double WEIMIN,EBREAK,DEPER,TIME;  // for aflowlib_libraries.cpp
  // begin shared xPOTCAR
  double ENMAX;deque<double> vENMAX;                            // eV
  double ENMIN;deque<double> vENMIN;                            // eV
  double POMASS_sum,POMASS_min,POMASS_max;deque<double> vPOMASS;// mass
  double ZVAL_sum,ZVAL_min,ZVAL_max;deque<double> vZVAL;        // valence
  double EATOM_min,EATOM_max;deque<double> vEATOM;              // eV
  double RCORE_min,RCORE_max;deque<double> vRCORE;              // outmost cutoff radius
  double RWIGS_min,RWIGS_max;deque<double> vRWIGS;              // wigner-seitz radius (au A)
  double EAUG_min,EAUG_max;deque<double> vEAUG;                 // augmentation
  // end shared xPOTCAR
  string pp_type;
  deque<string> species,species_pp,species_pp_type,species_pp_version; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  deque<deque<double> > species_pp_vLDAU;  // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  bool isKIN;                                                   // METAGGA
  bool isMETAGGA;string METAGGA;                                // METAGGA
  string string_LDAU;                                           // for aflowlib_libraries.cpp
  uint nweights,nkpoints_irreducible;                           // kpoints reading
  vector<aurostd::xvector<double> > vkpoint_reciprocal;         // kpoints reading
  vector<aurostd::xvector<double> > vkpoint_cartesian;          // kpoints reading
  vector<double> vweights;                                      // kpoints reading
  double calculation_time;                                      // for aflowlib_libraries.cpp - calculation_time
  double calculation_memory;                                    // for aflowlib_libraries.cpp - calculation_memory
  uint calculation_cores;                                       // for aflowlib_libraries.cpp - calculation_cores
  xstructure xstr;                                              // for GetBandGap()
  vector<string> GetCorrectPositions(string line,uint expected_count);                // 170725 CO - vasp issues with lattice spacing (negative sign) 
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);          // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                      // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                    // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,uint natoms_check,bool);       // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
  // EFFECTIVE MASSES
  friend bool GetEffectiveMass(xOUTCAR& outcar, xDOSCAR& doscar, xEIGENVAL& eigenval, xstructure xstr);
  vector<int> band_index;
  vector<int> carrier_spin;
  vector<string> carrier_type;
  vector<vector<double> > extrema_cart_coord;
  vector<vector<double> > effective_mass_axes;
  vector<int> equivalent_valley;
  vector<double> effective_mass_DOS;
  vector<double> effective_mass_COND;
  vector<double> mass_elec_dos;
  vector<double> mass_hole_dos;
  vector<double> mass_elec_conduction;
  vector<double> mass_hole_conduction;
  // BAND GAPS
  bool GetXStructure();
  int isKPointLine(uint iline,xvector<double>& kpoint); //if returns 0 if not KPointLine, -1 means it gave *** for kpoint
  int isKPointLine(uint iline);                         //if returns 0 if not KPointLine, -1 means it gave *** for kpoint
  bool GetStartingKPointLines(vector<uint>& ilines);
  bool GetNextKPointLine(uint& iline);
  bool ProcessKPoint(uint iline,double EFERMI,vector<double>& b_energies,vector<double>& b_occs);
  bool GetBandEdge(vector<double>& b_energies,vector<double>& b_occs,double EFERMI,uint& iedge,double efermi_tol=AUROSTD_NAN,double energy_tol=1e-4,double occ_tol=1e-5);
  bool identicalKPoints(vector<xvector<double> >& vkpoints,uint kpt1,uint kpt2,double tol=1e-12);
  bool identicalKPoints(xvector<double>& kpoint1,xvector<double>& kpoint2,double tol=1e-12);
  bool removeDuplicateKPoints(vector<xvector<double> >& vkpoints,vector<uint>& vikpt);
  bool removeDuplicateKPoints(vector<vector<xvector<double> > >& vkpoints,vector<uint>& vikpt,vector<uint>& vispin);
  double minimumDistanceKPoints(vector<xvector<double> >& vkpoints,uint ikp1,uint ikp2);
  double minimumDistanceKPoints(xvector<double>& kpoint1,xvector<double>& kpoint2);
  struct bandEnergyOcc{
    double energy;
    double occ;
    //bool operator<(const bandEnergyOcc& other) const {return energy<other.energy;}
  };
  struct bandEnergyOccCompare{
    bandEnergyOccCompare(double _energy_tol) : energy_tol(_energy_tol) {};
    double energy_tol;
    bool operator()(const bandEnergyOcc& a,const bandEnergyOcc b) const;
  };
  bool orderBands(vector<double>& b_energies,vector<double>& b_occs,double energy_tol=1e-4);
  enum BROAD_TYPES {empty,metal,insulator};                   //bandgap types
  enum EMPTY_TYPES {empty_all,empty_partial};                 //bandgap types
  enum INSULATOR_TYPES {insulator_direct,insulator_indirect}; //bandgap types
  enum GAP_TYPES {zero_gap,non_zero_gap};                     //bandgap types
  //bool GetBandGap(void); // CO 171002 - need POSCAR for kpt_tol
  bool GetBandGap(double EFERMI=AUROSTD_NAN,double efermi_tol=AUROSTD_NAN,double energy_tol=1e-4,double occ_tol=1e-5);
  bool GetBandGap_Camilo(double kpt_tol);
  vector<double> conduction_band_min;
  double         conduction_band_min_net;
  vector<double> valence_band_max;
  double         valence_band_max_net;
  vector<double> Egap;
  double         Egap_net;
  vector<double> Egap_fit;
  double         Egap_fit_net;
  vector<string> Egap_type;
  string         Egap_type_net;
  string ERROR;
  //int number_bands,number_kpoints; // CO 171006 - camilo garbage
  //int ISPIN; // turn this into spin = 0 if ISPIN = 1 // CO 171006 - camilo garbage
  //int spin;  // CO 171006 - camilo garbage
 private:                       //
  void free();                 // free space
  void copy(const xOUTCAR& b); //
};
//-------------------------------------------------------------------------------------------------
class xDOSCAR {
 public:
  xDOSCAR();                                                    // default, just allocate
  ~xDOSCAR();                                                   // kill everything
  xDOSCAR(const string& fileIN,bool=TRUE);                      // constructor from filename QUIET
  xDOSCAR(const xDOSCAR& b);                                    // constructor copy
  const xDOSCAR& operator=(const xDOSCAR &b);                   // copy
  void clear(void);                                             // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;       // the content, and lines of it
  string title;
  uint spin;
  double Vol,POTIM;
  xvector<double> lattice;
  double temperature;
  bool RWIGS;
  double Efermi;
  double spinF;
  double energy_max;
  double energy_min;
  uint number_energies;
  double denergy;
  deque<double> venergy;                                        // venergy.at(energy_number) 
  deque<double> venergyEf;                                      // venergyEf.at(energy_number) 
  deque<deque<double> > vDOS;                                   // vDOS.at(energy_number).at(spin)
  deque<deque<double> > viDOS;                                  // viDOS.at(energy_number).at(spin)
  deque<deque<double> > vDOSs;                                  // vDOSs.at(energy_number).at(spin)
  deque<deque<double> > vDOSp;                                  // vDOSp.at(energy_number).at(spin)
  deque<deque<double> > vDOSd;                                  // vDOSd.at(energy_number).at(spin)
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                                                        //
  void free();                                                  // free space
  void copy(const xDOSCAR& b);                                  //
};
//-------------------------------------------------------------------------------------------------
class xEIGENVAL {
 public:
  xEIGENVAL();                                                  // default, just allocate
  ~xEIGENVAL();                                                 // kill everything
  xEIGENVAL(const string& fileIN,bool=TRUE);                    // constructor from filename QUIET
  xEIGENVAL(const xEIGENVAL& b);                                // constructor copy
  const xEIGENVAL& operator=(const xEIGENVAL &b);               // copy
  void clear(void);                                             // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;       // the content, and lines of it
  string title;
  uint spin;
  double Vol,POTIM;
  xvector<double> lattice;
  double temperature;
  uint number_electrons,number_kpoints,number_bands;
  deque<double> vweight;                                        // vweight.at(kpoint number)
  deque<xvector<double> > vkpoint;                              // vkpoint.at(kpoint number)[1,2,3]=xyz.
  deque<deque<deque<double> > > venergy;                        // venergy.at(kpoint number).at(band number).at(spin number)
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                                                        //
  void free();                                                  // free space
  void copy(const xEIGENVAL& b);                                //
};
//-------------------------------------------------------------------------------------------------
class xPOTCAR {
 public:
  xPOTCAR();                                                    // default, just allocate
  ~xPOTCAR();                                                   // kill everything
  xPOTCAR(const string& fileIN,bool=TRUE);                      // constructor from filename QUIET
  xPOTCAR(const xPOTCAR& b);                                    // constructor copy
  const xPOTCAR& operator=(const xPOTCAR &b);                   // copy
  void clear(void);                                             // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;       // the content, and lines of it
  string title;
  bool   POTCAR_PAW;
  string POTCAR_TYPE;
  bool   POTCAR_KINETIC;
  double ENMAX;deque<double> vENMAX;                            // eV
  double ENMIN;deque<double> vENMIN;                            // eV
  double POMASS_sum,POMASS_min,POMASS_max;deque<double> vPOMASS;// mass
  double ZVAL_sum,ZVAL_min,ZVAL_max;deque<double> vZVAL;        // valence
  double EATOM_min,EATOM_max;deque<double> vEATOM;              // eV
  double RCORE_min,RCORE_max;deque<double> vRCORE;              // outmost cutoff radius
  double RWIGS_min,RWIGS_max;deque<double> vRWIGS;              // wigner-seitz radius (au A)
  double EAUG_min,EAUG_max;deque<double> vEAUG;                 // augmentation
  string pp_type;
  deque<string> species,species_pp,species_pp_type,species_pp_version; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                                                        //
  void free();                                                  // free space
  void copy(const xPOTCAR& b);                                  //
};
// -------------------------------------------------------------------------------------------------
class xVASPRUNXML {
 public:
  xVASPRUNXML();                                                // default, just allocate
  ~xVASPRUNXML();                                               // kill everything
  xVASPRUNXML(const string& fileIN,bool=TRUE);                  // constructor from filename QUIET
  xVASPRUNXML(const xVASPRUNXML& b);                            // constructor copy
  const xVASPRUNXML& operator=(const xVASPRUNXML &b);           // copy
  void clear(void);                                             // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;       // the content, and lines of it
  double natoms;                                                // for aflowlib_libraries.cpp
  xmatrix<double> stress;                                       // for aflowlib_libraries.cpp
  vector<aurostd::xvector<double> > vkpoint;                    // for aflowlib_libraries.cpp
  vector<aurostd::xvector<double> > vweights;                   // for aflowlib_libraries.cpp
  vector<aurostd::xvector<double> > vforces;                    // for aflowlib_libraries.cpp
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                       //
  void free();                 // free space
  void copy(const xVASPRUNXML& b); //
};
// -------------------------------------------------------------------------------------------------
class xIBZKPT {
 public:
  xIBZKPT();                                                    // default, just allocate
  ~xIBZKPT();                                                   // kill everything
  xIBZKPT(const string& fileIN,bool=TRUE);                      // constructor from filename QUIET
  xIBZKPT(const xIBZKPT& b);                                    // constructor copy
  const xIBZKPT& operator=(const xIBZKPT &b);                   // copy
  void clear(void);                                             // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;       // the content, and lines of it
  uint nweights;                                                // for aflowlib_libraries.cpp
  uint nkpoints_irreducible;                                    // for aflowlib_libraries.cpp
  vector<aurostd::xvector<double> > vkpoint;                    // for aflowlib_libraries.cpp
  vector<uint> vweights;                                        // for aflowlib_libraries.cpp
  uint ntetrahedra;                                             // for aflowlib_libraries.cpp
  double wtetrahedra;                                           // for aflowlib_libraries.cpp
  vector<aurostd::xvector<int> > vtetrahedra;                   // for aflowlib_libraries.cpp
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                       //
  void free();                 // free space
  void copy(const xIBZKPT& b); //
};
// -------------------------------------------------------------------------------------------------
class xKPOINTS {
 public:
  xKPOINTS();                                                    // default, just allocate
  ~xKPOINTS();                                                    // kill everything
  xKPOINTS(const string& fileIN,bool=TRUE);                      // constructor from filename QUIET
  xKPOINTS(const xKPOINTS& b);                                   // constructor copy
  const xKPOINTS& operator=(const xKPOINTS &b);                  // copy
  void clear(void);                                              // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;        // the content, and lines of it
  string title;                                                  // first line
  int mode;                                                      // sort of mode
  string grid_type;                                              // if grid specified
  bool is_KPOINTS_NNN,is_KPOINTS_PATH;                           // control parameters
  xvector<int>    nnn_kpoints; // N*N*N                          // triplet of kpoints
  xvector<double> ooo_kpoints; // ORIGIN                         // triplet of origin
  int nkpoints;                                                  // total kpoints
  string path_mode,path;vector<string> vpath;int path_grid;      // path if any
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                       //
  void free();                 // free space
  void copy(const xKPOINTS& b); //
};
// -------------------------------------------------------------------------------------------------
class xCHGCAR {
 public:
  xCHGCAR();                                                     // default, just allocate
  ~xCHGCAR();                                                    // kill everything
  xCHGCAR(const string& fileIN,bool=TRUE);                       // constructor from filename QUIET
  xCHGCAR(const xCHGCAR& b);                                     // constructor copy
  const xCHGCAR& operator=(const xCHGCAR &b);                    // copy
  void clear(void);                                              // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;        // the content, and lines of it
  xvector<int>     grid; // N*N*N                                // triplet of grid
  vector<string>   vstring; // ORIGIN                            // string of values
  xvector<double>  vvalues; // ORIGIN                            // xvector of values
  //[OBSOLETE ME180705]xtensor3<double> tvalues; // ORIGIN       // xtensor of values
  xtensor<double> tvalues; // ORIGIN                             // xtensor of values ME180705
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);          // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                      // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                    // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                       //
  void free();                 // free space
  void copy(const xCHGCAR& b); //
};

// -------------------------------------------------------------------------------------------------
// aflow_kaims.cpp
namespace KBIN {
  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(string& AflowIn,_aflags& aflags,_kflags& kflags);
  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(string& AflowIn,ofstream& FileMESSAGE,_aflags& aflags,_kflags& kflags);
  bool AIMS_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
}
// -------------------------------------------------------------------------------------------------
// aflow_iaims.cpp
namespace KBIN {
  bool AIMS_Produce_INPUT(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Modify_INPUT(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Write_INPUT(_xaims& xaims,_aimsflags &aimsflags);
  bool AIMS_Write_CONTROL(_xaims& xaims,_aimsflags &aimsflags);
  bool AIMS_Write_GEOM(_xaims& xaims,_aimsflags &aimsflags);
  bool AIMS_Produce_CONTROL(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Modify_CONTROL(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Reread_CONTROL(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags);
  bool AIMS_Produce_GEOM(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Produce_GEOM(_xaims& xaims);
  bool AIMS_Modify_GEOM(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_aimsflags &aimsflags);
  bool AIMS_Reread_GEOM(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags);
  bool XAIMS_CONTROL_PREPARE_GENERIC(string command,_xaims& xaims,_aimsflags& aimsflags,string svalue,int ivalue,double dvalue,bool OPTION);
  void XAIMS_CONTROL_REMOVE_ENTRY(_xaims& xaims,string ENTRY,string COMMENT,bool VERBOSE);
}
// -------------------------------------------------------------------------------------------------
// aflow_oaims.cpp
class xAIMSOUT;
class xAIMSOUT {
 public:
  xAIMSOUT();                                                    // default, just allocate
  ~xAIMSOUT();                                                   // kill everything
  xAIMSOUT(const string& fileIN,bool=TRUE);                      // constructor from filename, QUIET
  xAIMSOUT(const xAIMSOUT& b);                                    // constructor copy
  const xAIMSOUT& operator=(const xAIMSOUT &b);                   // copy
  void clear(void);                                             // clear
  // CONTENT
  string content;vector<string> vcontent;string filename;       // the content, and lines of it
  vector<aurostd::xvector<double> > vforces;                    // for aflowlib_libraries.cpp
  double natoms;
  string ERROR;
  bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);          // get everything QUIET
  bool GetProperties(const string& stringIN,bool=TRUE);                      // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,bool=TRUE);                    // get everything QUIET
  bool GetPropertiesFile(const string& fileIN,uint natoms_check,bool);       // get everything QUIET
  bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
 private:                       //
  void free();                 // free space
  void copy(const xAIMSOUT& b); //
};
// -----------------------------------------------------------------------------------------------
bool PrintBandGap       (string& WorkDir, ostream &oss);
bool PrintEffectiveMass (string& WorkDir, ostream &oss);
bool PrintEigCurv       (string& WorkDir, ostream &oss);
// -----------------------------------------------------------------------------------------------
bool ParseKPOINTS(stringstream& File_Kpoints, int& GRIDS, vector<xvector<double> >& special_kpts, vector<xvector<double> >& unique_kpts, vector<int>& repeat_kpts_num);
bool AdjacencyList_KPT(vector<xvector<double> >& special_kpts,vector<xvector<double> >& unique_kpts,vector<xvector<int> >& connect_kpts,vector<int>& connect_kpts_num);
bool AdjacencyList_EIG(vector<xvector<double> >& unique_kpts,vector<xvector<int> >& connect_kpts,vector<int>& connect_kpts_num,xEIGENVAL& xeigenval,vector<xvector<double> >& unique_kpts_EIG,vector<xvector<int> >& connect_kpts_EIG,vector<xvector<double> >& vkpoint_eig);
bool RepeatsList(vector<xvector<double> >& unique_kpts_EIG,vector<int>& repeat_kpts_num,vector<xvector<double> >& vkpoint_eig,vector<xvector<int> >& repeat_kpts_EIG);
bool VertexPaths(vector<xvector<int> >& repeat_kpts_EIG,vector<xvector<int> >& connect_kpts_EIG,vector<int>& repeat_kpts_num,int& GRIDS,vector<xvector<int> >& vrtx_path);
bool RepeatedEdges(vector<xvector<int> >& vrtx_path,vector<xvector<int> >& repeat_kpts_EIG,vector<int>& repeat_kpts_num,vector<xvector<int> >& ndx_edges);
bool VertexBranches(vector<xvector<int> >& ndx_edges,vector<int>& repeat_kpts_num,vector<xvector<int> >& repeat_kpts_EIG,vector<vector<xvector<int> > >& branches);
bool PathDataStuct(xEIGENVAL& xeigenval,vector<xvector<double> >& vkpoint_eig,vector<vector<xvector<int> > >& branches,vector<vector< vector<int> > >& branches_indx,vector<vector< vector<xvector<double> > > >& branches_kpts,vector<vector< vector<vector<vector<double> > > > >& branches_bnds);
bool IBZextrema(xEIGENVAL& xeigenval, vector<xvector<double> >& vkpoint_eig, vector<vector<xvector<int> > >& branches);
void CompareDoublesChar(bool& MATCH, double& number1, double& number2);
void CompareEdges(vector<vector<xvector<int> > >& branches, vector<xvector<int> >& vertex_edges, xvector<int>& test_edge, bool& MATCH);
void NaiveCurvatures(xvector<double>& eigvec, vector<xvector<double> >& posvec, vector<double>& curvature);
double StencilLinear1D(vector<xvector<double> >& positions, xvector<double>& eigenvals);
//-------------------------------------------------------------------------------------------------
struct kEn_st {
  xvector<double> kpoint;
  double energy[2];
  int band_index;
  int band_type; // 0 -- valence band; 1 -- conduction band
};
#define _SIGMA 1.0 // default standard deviation of input data
// range of energy point to fit the ellipse curve
const double _FIT_ENERGY_RANGE = 0.026; // eV range of band
const int _FIT_POINTS_NUMBER = 8; // minimum fit points in Irreducible BZ
//range of band extremes to determine the number of bands for effective mass calculations
const double _BANDS_ENERGY_RANGE = 0.026; // eV
// used to determine cluster of points can be changed to other values
const double _BANDS_PARAMETER_MIN_RATIO = 0.2;
// factor unit
// mass is in unit of electron mass
const double _MASS_FACTOR = 3.80998; // hbar^2*10^{20}/(2.0*me*eV)
bool comparison_kEn_str_up          (const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_dn          (const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_position    (const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_band_type_up(const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_band_type_dn(const kEn_st& k1, const kEn_st& k2);
bool is_equal_position_kEn_str      (const kEn_st& k1, const kEn_st& k2);
bool near_to                        (const xvector<double> & k1, const xvector<double> & k2, const vector<double> & max_distance);
// [OBSOLETE] bool GetEffectiveMass(xOUTCAR& outcar,xDOSCAR& doscar,xEIGENVAL& eigenval,xstructure xstr,ostream& oss,const bool& osswrite);
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// aflow_estructure_dos.cpp

namespace aurostd {
  int CountWordsinString(string& input);  // put aurostd
  int CountWordsinString_web(string input); // put aurostd
}

namespace estructure {
  string PEDOS_GENERATE_GNUPLOTSCRIPT(const string&,const string&,const double&,const double&,const double&,const double&,const int&,const vector<vector<vector<double> > >&,const string&);
  bool isSpecialKPOINT(string kpoint);  // CO 170830
  string fixSpecialKPOINT_GNUPLOT(string kpoint,bool json=false);  // CO 170830
  string fixSpecialKPOINT_HTML(string kpoint);  // CO 170830
  string fixSpecialKPOINT_LATEX(string kpoint);  // CO 170830
  string fixKPOINT_GNUPLOT(string kpoint,bool json=false);  // CO 170830
  string fixKPOINT_HTML(string kpoint);  // CO 170830
  string fixKPOINT_LATEX(string kpoint);  // CO 170830
  string fixKPOINT_SPECIALONLY(string kpoint);  // CO 170830
  void PLOT_BANDDOS(string options);
  void PLOT_BAND(string options);
  void PLOT_DOS(string options);
  void PLOT_PEDOS(string options);
  void PLOT_PEDOSALL(string options);
  void PLOT_PEDOSALL_AFLOWLIB(string options,_aflags& aflags);
  void PLOT_BAND2(string options);
  // [OBSOLETE]  void PLOT_BAND3(string options);
  void PLOT_BAND_SPINSPLIT(string options);
  void PLOT_DOSWEB(string options);
  // manipulation
  string changeICSDNameGunplot(string ICSDName);
  void CombineTDOSAndTOTALPDOS(const vector<vector<double> >& TDOS, const vector<vector<double> >& TOTALPDOS, vector<vector<double> >& vvDOS);
  double GET_TDOSDATA(const string& str_dir, vector<vector<double> >& TDOS);
  double GET_TDOSDATA(stringstream& ss_dosfile, stringstream& ss_outcarfile, vector<vector<double> >& TDOS);
  double GET_TOTALPDOSDATA(const string& str_dir, vector<vector<double> >& TOTALPDOS);
  double GET_TOTALPDOSDATA(stringstream& ss_dosfile, stringstream& ss_outfile, vector<vector<double> >& TOTALPDOS);
  double GET_PDOSDATA(const string& str_dir, vector<vector<vector<double> > >& PDOS);
  double GET_PDOSDATA(stringstream& ss_dosfile, stringstream& ss_outfile, vector<vector<vector<double> > >& PDOS);
  // [OBSOLETE]  void GET_DOS_DATA(vector<string>& argv);
  bool GET_DOS_DATA(stringstream& ss_dosfile, stringstream& ss_outfile, double& Efermi, vector<vector<double> >& TDOS, vector<vector<double> >& TOTALPDOS); // CO 180216
  bool GET_DOS_DATA(const string& str_dir,   double& Efermi, vector<vector<double> >& TDOS, vector<vector<double> >& TOTALPDOS, vector<vector<vector<double> > >& PDOS);  // CO 180216
  bool GET_DOS_DATA(stringstream& ss_dosfile, stringstream& ss_outfile,   double& Efermi, vector<vector<double> >& TDOS, vector<vector<double> >& TOTALPDOS, vector<vector<vector<double> > >& PDOS); // CO 180216
  void FormatSpinofPDOS(vector<vector<vector<double> > >& vvva);

  // Functions for serializing bands data to JSON
  // Added by Eric G
  bool DOSDATA_JSON(aurostd::xoption& vpflow,ostream& oss=cout);
  bool DOSDATA_JSON(aurostd::xoption& vpflow,string directory,stringstream& json,bool wrapping_brackets=true);
  bool BANDSDATA_JSON(aurostd::xoption& vpflow,ostream& oss=cout);
  bool BANDSDATA_JSON(aurostd::xoption& vpflow,string directory,stringstream& json,bool wrapping_brackets=true);
  //uint DOSDATA_JSON(string options);
  //uint DOSDATA_JSON(string options, ostream& json);
  //uint BANDSDATA_JSON(string options);
  //uint BANDSDATA_JSON(string options, string json_dir);
  //uint BANDSDATA_JSON(string options, ostream& json);
  string linelabel2HTML(string linelabel);
  uint inequivalentAtomsJSON( vector<vector<vector<double> > >& PDOS, vector<int>& iatoms, vector<double>& numbers, vector<string>& vspecies, ostream& json);
  uint constructInequivalentAtomPDOSJSON(vector<vector<vector<double> > >& PDOS, int iatom, ostream& json); 
  // End of bands data JSON serializers

}


// ----------------------------------------------------------------------------
// aflow_poccupation_*.cpp
//  #include "aflow_pocc.h"

// aflow_poccupation_params.cpp
namespace pocc {
  bool poccInput(); //170805 CO
  
  string ReturnAtomSpecies(string atom);
  string ReturnAtomSpeciesPotential(string atom);
  string ReturnUFFParameters(string atom);
  class UFFPara {
  public:
    UFFPara(); // constructor
    ~UFFPara(); // destructor
    string symbol;
    double r1,theta0,x1,D1,zeta,Z1,Vi,Uj,Xi,hard,radius; 
    void GetUFFParameters(string);
  private:
    void Free(); //free space
  };
  string ReturnAtomProperties(string atom);
  //Atomic Properties Database
  class Atom {
  public:
    Atom();
    ~Atom();
    string name,symbol;
    int number; //atomic number
    double mass,radius,Xi; //atomic, weight radius /pauling electronegativity
    void GetAtomicProperties(string);
  private:
    void Free();
  };
} // namespace pocc

// aflow_poccupation_forcefield.cpp
namespace pocc {
  class Bond{
  public:
    Bond();
    ~Bond();
    _atom bgn,end;
    double length;
    void Set(xstructure , _atom , _atom );
    const Bond & operator=(const Bond &other);
    bool operator==(const Bond &other) const;
    bool operator!=(const Bond &other) const;
    friend ostream& operator<<(ostream&,const Bond&);
  private:
    void Free();
  };
  void SetUFFPara(_atom atomi, _atom atomj, double& R0, double& Kij, double& Xij, double& Dij);
  double CalculateBondEnergy(xstructure xstr, _atom atomi, _atom atomj);
  double CalculateNonBondEnergy(xstructure xstr, _atom atomi, _atom atomj);
  double CalculateUFFEnergy(xstructure xstr);
  void RemoveSameBond(vector<Bond>& Bonds_orig, vector<Bond>& Bonds_new);
  void ExtractBonds(const xstructure& xstr, deque<deque<_atom> >& neigh_mat_bonded, deque<deque<_atom> >& neigh_mat_nonbonded);
  void AnalyzeBonds(const xstructure& xstr, vector<Bond>& Bonds, vector<Bond>& NonBonds);
  void UFFENERGY(istream& input);
}

// ----------------------------------------------------------------------------
// aflow_mix.cpp  aflow_nomix.cpp   aflow_mix_pauling.cpp
#define MISCIBILITY_SYSTEM_NOT_STUDIED  3
#define MISCIBILITY_SYSTEM_SOLUTION     2
#define MISCIBILITY_SYSTEM_MISCIBLE     1
#define MISCIBILITY_SYSTEM_NOMIX        0
#define MISCIBILITY_SYSTEM_UNKNOWN     -1
#define MISCIBILITY_SYSTEM_CUTOFF     200
#define MIEDEMA_MIX_SLOPE 3.069              // Miedema Rule Table 1a Physica 100B (1980) 1-28

int MiscibilityCheck(int speciesA,int speciesB);                  // aflow_mix.cpp
int MiscibilityCheck(string speciesA,string speciesB);            // aflow_mix.cpp
int MiscibilityExperimentsCheck(int speciesA,int speciesB);       // aflow_mix.cpp
int MiscibilityExperimentsCheck(string speciesA,string speciesB); // aflow_mix.cpp
int MiscibilityMiedemaCheck(int speciesA,int speciesB);           // aflow_mix.cpp
int MiscibilityMiedemaCheck(string speciesA,string speciesB);     // aflow_mix.cpp
int MiscibilityMiedemaCheck(string system_in);                    // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(int speciesA,int speciesB);       // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(string speciesA,string speciesB); // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(string system_in);                // aflow_mix.cpp
int MiscibilityCheck(string system_in);                           // aflow_nomix.cpp
int MiscibilityExperimentsCheck(string system_in);                    // aflow_mix_pauling.cpp

// ----------------------------------------------------------------------------
// symmetry prototypes
// aflow_symmetry.cpp
namespace SYM {
  // DX AND CO - START
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xstructure& a); // DX AND CO
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xstructure& a, bool _incell_,bool roff); // DX AND CO
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xstructure& a,bool skew, bool _incell_,bool roff,double _eps_); // DX AND CO
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew, bool _incell_,bool roff,double _eps_); // DX AND CO
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str); // DX
  // DX and CO - END
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_);
  // DX and CO - START
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_,bool roff); // DX
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_,bool roff,bool validatePosition); // DX
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew,bool _incell_,bool roff,bool validatePosition,double eps); // DX
  _atom ApplyAtom_20161115(const _atom& atom_in,const _sym_op& symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew,bool _incell_,bool roff,bool validatePosition,double eps); // DX
  _atom ApplyAtom_20160101(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_); // DX
  // DX and CO - END
  xvector<double> ApplyCpos(const xvector<double> &cpos_in,const _sym_op &symop,const xstructure& str,bool _incell_);
  xvector<double> ApplyCpos(const xvector<double> &cpos_in,const _sym_op &symop,const xstructure& str);
  xvector<double> ApplyFpos(const xvector<double> &fpos_in,const _sym_op &symop,const xstructure& str,bool _incell_);
  xvector<double> ApplyFpos(const xvector<double> &fpos_in,const _sym_op &symop,const xstructure& str);
  xvector<int> ApplyIJK(const xvector<int> &ijk_in,const _sym_op &symop,const xstructure& str);
  int  ApplyL(const int &l_in,const _sym_op &symop,const xstructure& str);
  xstructure ApplyXstructure(const _sym_op &symop,const xstructure& str);
  xstructure ApplyXstructure(const _sym_op &symop,const xstructure& str,bool _incell_);
  // DX and CO - START
  bool AtomsEquivalent(xstructure& str,_atom& atom1,_atom& atom2); // DX
  bool AtomsEquivalent_20161115(xstructure& str,_atom& atom1,_atom& atom2); // DX
  bool AtomsEquivalent_20160101(xstructure& str,_atom& atom1,_atom& atom2); // DX
  bool AtomsEquivalent(xstructure& str, _atom& atom1, _atom& atom2,double& eps); // DX
  bool AtomsEquivalent(xstructure& str, _atom& a, _atom& b, bool& skew, double& tol); // DX
  bool AtomsEquivalent_20161115(xstructure& str,_atom& atom1,_atom& atom2,double& eps); // DX
  bool AtomsEquivalent_20160101(xstructure& str,_atom& atom1,_atom& atom2,double& eps); // DX
  bool AtomsEquivalent_Basis(xstructure& str, int atom1_indx,int atom2_indx);
  // DX and CO - END
  bool CposEquivalent(const xstructure& str,const xvector<double>& cpos1,const xvector<double>& cpos2,const double& eps);
  bool FposEquivalent(const xstructure& str,const xvector<double>& fpos1,const xvector<double>& fpos2,const double& eps);
  bool CposEquivalent(const xstructure& str,const xvector<double>& cpos1,const xvector<double>& cpos2);
  bool FposEquivalent(const xstructure& str,const xvector<double>& fpos1,const xvector<double>& fpos2);
  bool TypePointGroupOperation(const xmatrix<double>& Uc,const xmatrix<double>& Uf,string& _string,bool& _inversion,double& _angle,
			       xvector<double>& _axis,xmatrix<double>& _generator, xvector<double>& _generator_coefficients, 
                               xmatrix<xcomplex<double> >& _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double _eps_);  // calculate the symmetry inversion,type,axis,generator // DX 12/6/17 - Added generator coefficients // DX 12/7/17 - Added Uf // DX 1/17/18 - Added SU2 and su2 coefficients
  bool TypePointGroupOperationInternational(const xmatrix<double>& Uc,string& _stringHM,string& _stringSC,
					    const bool& _inversion,const double& _angle,
					    const xvector<double>& _axis,const xmatrix<double>& _generator, xvector<double>& _generator_coefficients, 
                                            xmatrix<xcomplex<double> >& _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double _eps_);  // International symbol = Hermann-Mauguin notation & Schonflies notation // DX 12/6/17 - Added generator coefficients // DX 1/17/18 - Added SU2 and su2 coefficients
  // DX and CO - START
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,
			      const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
			      const xvector<double> &ctrasl,const xvector<double> &ftrasl,
			      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group);
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,
			      const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
			      const xvector<double> &ctrasl,const xvector<double> &ftrasl,
			      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff); // DX
  uint AddSymmetryToStructure(xstructure& a,const xmatrix<double>& Uc,const xmatrix<double>& Uf,
			      const xvector<double>& ctau,const xvector<double>& ftau,const xvector<double>& ctrasl,
			      const xvector<double>& ftrasl,
			      const std::vector<int>& basis_atoms_map,const std::vector<int>& basis_types_map,bool basis_map_calculated,char group);
  uint AddSymmetryToStructure(xstructure& a,const xmatrix<double>& Uc,const xmatrix<double>& Uf,
			      const xvector<double>& ctau,const xvector<double>& ftau,const xvector<double>& ctrasl,
			      const xvector<double>& ftrasl,
			      const std::vector<int>& basis_atoms_map,const std::vector<int>& basis_types_map,bool basis_map_calculated,char group,bool roff); // DX
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
			      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group);
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
			      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff); // DX
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
			      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group); // DX
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
			      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff); // DX
  bool PointGroupsIdentical(const vector<_sym_op>& vpg1,const vector<_sym_op>& vpg2, double eps, bool is_same_lattice=false); // DX 12/7/17 - added is_same_lattice
  //GEENA START
  bool CalculateQuaternion(_sym_op& a);
  //GEENA STOP
  bool ComplexSU2Rotations(xmatrix<xcomplex<double> > & _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double& theta, xvector<double>& _axis); // DX 1/17/18 - add SU(2) and su(2) coefficients
  // DX and CO - END
  bool CalculatePointGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");      // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const xmatrix<double>& lattice, vector<_sym_op > pgroup, ofstream &FileMESSAGE,bool _write_,const bool& osswrite,ostream& oss,double _eps_);
  uint CalculatePointGroup(const xmatrix<double>& lattice, vector<_sym_op > pgroup, bool _write_,const bool& osswrite,ostream& oss,double _eps_);     // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const xmatrix<double>& lattice,double _eps_);     // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const xmatrix<double>& lattice);     // POINT GROUP      _PGROUP_
  bool CalculatePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_
  // DX and CO - START
  bool CalculatePointGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_); // DX
  bool CalculatePointGroup_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // DX
  // DX and CO - END
  bool CalculatePointGroupKlattice(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");  // POINT GROUP KLATTICE     _PGROUPK_
  bool CalculatePointGroupKCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");  // POINT GROUP KCRYSTAL     _PGROUPK_XTAL_ // DX 12/5/17 - New group: reciprocal space counterpart of pgroup_xtal
  bool TransformSymmetryFromRealToReciprocal(ofstream &FileMESSAGE, xstructure& real_space_crystal, xstructure& reciprocal_space,
                                             _aflags& aflags, const bool& osswrite, ostream& oss, string& pgroup_type); // DX 8/8/17 - New klattice routine // DX 12/5/17 - Added pgroup_type option to account for pgroupk_xtal
  bool CalculateSitePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");  // SITE POINT GROUP _AGROUP_
  // DX and CO - START
  bool CalculateSitePointGroup(ofstream &FileMESSAGE,xstructure &a,int CALCULATION_MODE,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt"); // SITE POINT GROUP _AGROUP_
  bool CalculateSitePointGroup_20160801(ofstream &FileMESSAGE,xstructure &a,int CALCULATION_MODE,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // SITE POINT GROUP _AGROUP_ // DX
  bool CalculateSitePointGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_); // SITE POINT GROUP _AGROUP_ // DX
  bool CalculateSitePointGroup_EquivalentSites(xstructure &a,double _eps_); // DX
  bool CalculateSitePointGroup_EquivalentSites(xstructure &a,bool get_full_basis,double _eps_); // DX
  // DX and CO -END
  bool CalculatePointGroupCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");     // POINT GROUP      _PGROUP_
  bool CalculatePointGroupCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_
  // DX and CO -START
  bool CalculatePointGroupCrystal_20170814(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_ // DX
  bool CalculatePointGroupCrystal_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_ // DX
  bool CalculatePointGroupCrystal_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_);      // POINT GROUP      _PGROUP_ // DX
  bool PointGroupMap(xstructure& a, string& pgname, string& operations, char group); // DX 9/6/17
  bool PointGroupLookUpTable(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format);
  // DX and CO -END
  void CalculateSitePointGroup2(xstructure &a,bool ComMidss); // for --agroup2 and --agroup2m
  // DX - START
  //xstructure and _sym_op
  bool getFullSymBasis(xstructure& a, _sym_op& symOp,bool map_types,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  bool getFullSymBasis(xstructure& a, _sym_op& symOp,bool map_types,double& tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  bool getFullSymBasis(xstructure& a, _sym_op& symOp,bool map_types,bool& skew,double& tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  bool getFullSymBasis(xstructure& a, _sym_op& symOp,bool map_types,bool& skew,double& tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  //atoms, c2f, f2c and _sym_op
  bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& lattice,xmatrix<double>& c2f, xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool& skew,double& tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& lattice,xmatrix<double>& c2f, xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool& skew,double& tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  bool getFullSymBasis_20170729(deque<_atom>& atoms,xmatrix<double>& lattice,xmatrix<double>& c2f, xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool& skew,double& tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  /*bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map);
    bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, string& str_Hermann_Mauguin, bool& skew, double& tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map);
    bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, xvector<double>& ftau, bool& skew, double& tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map);
    bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, string& str_Hermann_Mauguin, xvector<double>& ftau, bool& skew, double& tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map);*/
  bool CalculateFactorGroup_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");
  // DX - END
  bool CalculateFactorGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");     // FACTOR GROUP     _FGROUP_
  bool CalculateFactorGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // FACTOR GROUP      _FGROUP_
  // DX - START
  bool AtomsMapped(_atom& a, _atom& b, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);  
  bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, xvector<double>& out, const xmatrix<double>& c2f, const xmatrix<double>& f2c, double& tol);
  bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, xvector<double>& out, const xmatrix<double>& c2f, const xmatrix<double>& f2c, xvector<int>& ijk, bool& restriction, double& tol);
  double minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice);
  double minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice,xvector<double>& min_vec,xvector<int>& ijk);
  xvector<double> minimumCartesianVector(const xvector<double>&, const xvector<double>&,
                                         const xmatrix<double>&);  // ME 180730
  xvector<double> minimumCartesianVector(const xvector<double>&, const xvector<double>&,
                                         const xmatrix<double>&, xvector<int>&);  // ME180730
  bool PBC(xvector<double>& v_in, xvector<int>& ijk, bool& restriction);
  bool PBC(xvector<double>& v_in);
  bool AtomFPOSMatch(deque<_atom>& atom_set, _atom& atom2, uint& match_type, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);
  bool AtomFPOSMatch(_atom& atom1, _atom& atom2, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);
  bool AtomFPOSMatch(xvector<double>& atom1, xvector<double>& atom2, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);
  bool validateAtomPosition(const _atom& atom,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool& skew,double& _eps_);
  bool validateAtomPosition(const xvector<double>& cpos,const xvector<double>& fpos,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool& skew,double& _eps_);
  bool MapAtom(deque<_atom>& a_vec, _atom& b, bool map_types, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);
  bool MapAtom(vector<_atom> a_vec, _atom b, bool map_types, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);
  bool MapAtom(xvector<double>& a, xvector<double>& b, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);
  bool MapAtom(_atom& a, _atom& b, bool map_types, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol);
  bool BringInCellTolABC(xstructure& a, xvector<double> tol_abc_res);
  bool MapAtomWithBasis(vector<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol,bool fast=true);   
  bool MapAtomWithBasis(vector<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol, uint& mapped_index,bool fast=true);  
  bool MapAtomWithBasis(deque<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol,bool fast=true);
  bool MapAtomWithBasis(deque<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, double& tol, uint& mapped_index,bool fast=true);
  bool isLatticeSkewed(const xmatrix<double>& lattice, double& min_dist, double& tol);
  double minimumDistance(const xstructure& xstr);
  double minimumDistance(const deque<_atom>& atoms, const xmatrix<double>& lattice,double scale=1.0);
  double defaultTolerance(const xstructure& xstr);
  bool checkAngle(xvector<double>& v1, xvector<double>& v2, double input_angle, double& tolerance);
  bool checkAngle(xvector<double>& v1, xvector<double>& v2, double input_angle, bool& is_deg, double& tolerance);
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, double& tolerance);
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, bool& is_deg, double& tolerance);
  // DX 9/5/17 [OBSOLETE] bool change_tolerance(xstructure& xstr, double& tolerance, double& orig_tolerance, int& count , double& min_dist, bool& no_scan);
  bool change_tolerance(xstructure& xstr, double& tolerance, double& min_dist, bool& no_scan);
  deque<deque<_atom> > break_up_by_type(deque<_atom>& expanded_crystal);
  vector<vector<_atom> > break_up_by_type(vector<_atom> expanded_crystal);
  double mod_one(double d); // DX 
  _atom mod_one_atom(const _atom& atom_in); // CO
  xvector<double> mod_one_xvec(xvector<double> a); // DX
  bool CheckForIdentity(const xstructure& xstr); // DX
  bool checkSuperCellLatticePoints(xstructure& xstr, int& num_lattice_points, char& centering, uint& expand_size); // DX
  bool ComparePointGroupAndSpaceGroupString(xstructure& xstr, int& multiplicity_of_primitive, bool& derivative_structure); // DX
  bool CalculateFactorGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_);      // FACTOR GROUP      _FGROUP_
  bool CalculateSpaceGroup_20160101(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss);      // SPACE GROUP      _SGROUP_
  bool CalculateSpaceGroup_20160801(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");      // SPACE GROUP      _SGROUP_
  // DX - END
  bool CalculateSpaceGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");      // SPACE GROUP      _SGROUP_

  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt"); // EQUIVALENT ATOMS _IATOMS_
  // DX and CO - START
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ // DX
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ // DX
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ // DX
  bool CalculateInequivalentAtoms_20160801(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ // DX  
  bool CalculateInequivalentAtoms_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss); // EQUIVALENT ATOMS _IATOMS_ // DX
  // DX and CO - END
}
string AgroupSymmetryToJson(vector<vector<_sym_op> >& group, char& mode); // DX 8/3/17 - For Python wrapper
string EquivalentAtomsToJson(vector<vector<int> >& iatoms); // DX 8/3/17 - For Python wrapper
string SymmetryToJson(vector<_sym_op>& group, char& mode); // DX 8/3/17 - For Python wrapper
bool KBIN_SymmetryWrite(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,char group,const bool& osswrite,ostream& oss,const string& format="txt");
//bool KBIN_SymmetryToScreen(xstructure& a, string& format, ostream& oss); // DX 8/3/17 - For Python wrapper
bool KBIN_SymmetryToScreen(xstructure& a, string& format, ostream& oss, char mode='\0'); // DX 8/22/17 - For Python wrapper
bool KBIN_StepSymmetryPerform(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss);
// DX and CO - START
bool KBIN_StepSymmetryPerform_20161205(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss); // DX
bool KBIN_StepSymmetryPerform_20160101(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss); // DX
// DX and CO - END
vector<double> PointGroupHistogramCheck(xstructure& a);
//bool SYM_CalculatePointGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss);      // MIKNOWSKI BASIS REDUCTION

// ----------------------------------------------------------------------------
// aflow_spacegroup.cpp
// spagegroup tables and routines
extern string LibrarySPACEGROUP;
namespace spacegroup{
  class _spacegroup {
  public:
    // constructors/destructors
    _spacegroup(void);    // do nothing
    _spacegroup(const _spacegroup& b);    // do nothing
    ~_spacegroup();        // do nothing
    // OPERATORS                                                  // --------------------------------------
    const _spacegroup& operator=(const _spacegroup& b);             // some operators
    // CONTENT
    uint number;
    uint option;
    string stroption;
    string name;
    string sginfo;
    std::vector<_sym_op> fgroup;                                  // rotations/inversions + incell_translations operations
    std::vector<_sym_op> pgroup;                                  // rotations/inversions
  private:                                                       // ---------------------------------------
    void Free();                                                  // to free everything
    void Copy(const _spacegroup& b);                               // the flag is necessary because sometimes you need to allocate the space.
  };

  extern std::vector<_spacegroup> vspacegroups;
  uint SpaceGroupInitialize(void);
  bool SpaceGroupNumberStructure(xstructure &str);
  bool SpaceGroupOptionRequired(uint spacegroup);
}

// ----------------------------------------------------------------------------
// aflow_lattice.cpp
// lattice and brillouin zones
namespace LATTICE {
  bool lattice_is_working(string lat);
  string SpaceGroup2Lattice(uint sg);
  uint Lattice2SpaceGroup(string lattice,vector<uint>& vsg);
  string SpaceGroup2LatticeVariation(uint sg,const xstructure& str);
  string ConventionalLattice_SpaceGroup(uint sg,double a,double b,double c);
  string ConventionalLattice_SpaceGroup(uint sg,const xstructure& str);
  xvector<double> Getabc_angles_Conventional(const xmatrix<double>& rlattice, string lattice,const int& mode);
  bool fix_sts_sp(xstructure& str_sp,xmatrix<double> &rlattice,xmatrix<double> &plattice);
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true);
  bool Standard_Lattice_StructureDefault(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true);
  bool Standard_Lattice_StructureCoarse(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructureNormal(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructureMedium(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructurePrecise(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructureUltra(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  //bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang); // STEFANO OLD VERSION
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps);
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps,bool histogram);
  bool Standard_Lattice_Structure_20170101(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps,bool histogram);
  bool Standard_Lattice_Structure_20170718(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true);
  //bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,int mode);
  bool Bravais_Lattice_Structure(xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang); // calculate everything
  bool Bravais_Lattice_StructureDefault(xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true); // calculate everything
  // DX - START
  bool Bravais_Lattice_StructureDefault_20170401(xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true); // calculate everything
  bool Bravais_Lattice_StructureDefault_20160101(xstructure& str_in,xstructure& str_sp,xstructure& str_sc); // calculate everything
  // DX - END
  bool Lattice(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,string& bravais_lattice_type,string& bravais_lattice_variation_type,string& bravais_lattice_system,double eps,double epsang);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,double eps,double epsang);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,double eps,double epsang);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice);
  string Bravais_Lattice_System(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang);
  string Bravais_Lattice_System(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc);
  string Bravais_Lattice_System(const xmatrix<double>& lattice,double eps,double epsang);
  string Bravais_Lattice_System(const xmatrix<double>& lattice);
  string Primitive_Lattice_Type(const xstructure& str);
  string Bravais_Lattice_System(const xstructure& str);
  string Conventional_Lattice_Type(const xstructure& str);
  xstructure Standard_Primitive_Lattice_Structure(const xstructure& str);
  xstructure Standard_Conventional_Lattice_Structure(const xstructure& str);
  string Get_Primitive_Lattice_Structure(const xstructure& str);
  xmatrix<double> sc2sp(const xmatrix<double>& rlattice, string lattice,bool inverseflag);
  xmatrix<double> sp2sc(const xmatrix<double>& rlattice, string lattice,bool inverseflag);
  void BZPLOTDATA(string options,istream& poscar, int mode);
}
xvector<double> Vrotate(xvector<double> v, xvector<double> axisrot, double theta);
void CheckLatticeHistogram();
namespace LATTICE {
  // kpoints and brillouin zones
  string KPOINTS_Directions(xstructure str_in,double grid,bool &foundBZ);
  string KPOINTS_Directions(string lattice_type,xmatrix<double> sp, double _grid,int iomode,bool &foundBZ);
  string KPOINTS_Directions(string lattice_type,xmatrix<double> sp, xmatrix<double> transformation_matrix, double _grid,int iomode,bool &foundBZ); //DX 20181101
}

// ----------------------------------------------------------------------------
// neighbours prototypes
// aflow_neighbours.cpp
bool StepNeighboursPerform(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);

// ----------------------------------------------------------------------------
// surface prototypes
// aflow_surface.cpp
namespace surface {
  double PointInTriangleContribution(const xvector<double>& _point,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  double PointInRhombusContribution(const xvector<double>& _point,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& v4);
  double TriangleArea(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  bool PlaneGetABCD(double& a,double& b,double& c,double& d,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  double PlaneDistance(const xvector<double>& r,const double& a,const double& b,const double& c,const double& d);
  double PlaneDistance(const xvector<double>& r,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  xvector<double> PlaneGetProjection(const xvector<double>& r,const double& a,const double& b,const double& c,const double& d);
  xvector<double> PlaneGetProjection(const xvector<double>& r,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  xvector<double> PlaneGetHKL(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3);
  bool PlaneGetVVV(const xvector<double>& hkl,double& area,xvector<double>& v1,xvector<double>& v2,xvector<double>& v3,xvector<double>& v4,const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3);
  double GetPlaneDensityAtoms(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const int& type_at);
  double GetPlaneDensityAtoms(const xstructure& _str,const xvector<double>& hkl,const double& roughness);
  double GetPlaneDensityBBonds(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const double& bbdistance,const int& type_at1,const int& type_at2);
  double GetPlaneDensityBBonds(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const double& bbdistance);
  double GetNNeighbours(const xstructure& _str,const int& type_at1,const int& type_at2);
  double GetNNeighbours(const xstructure& _str);
  string PrintHKLSigma(int num_types,int num_types_combinations);
  string PrintHKLSigmaBB(int num_types,int num_types_combinations,const double& bbfrac,const double& bbdistance,const xmatrix<double>& bbdistances);
  bool GetSurfaceHKL(const xstructure& _str,_aflags& aflags,const xvector<double>& _hkl,vector<vector<double> >& planesreducible,vector<vector<double> >& planesirreducible,ostream& oss);
  bool GetSurfaceHKLSearch(const xstructure& _str,_aflags& aflags,const xvector<double>& iparams,vector<vector<double> >& planesreducible,vector<vector<double> >& planesirreducible,vector<vector<uint> >& planesirreducible_images,ostream& oss,const string& smode);
}

namespace slab { // ROMAN CHEPULSKYY
  xstructure MAKE_SLAB(string options, istream& cin);
  xstructure MAKE_SLAB(string options, xstructure& str_in);
  void POSCAR_reading(istream& cin);
  double VectorAbsValue(int Layer, int NinLayer,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords);
  double VectorScalarMult(int Layer1, int NinLayer1, int Layer2, int NinLayer2,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords);
  double CosAngle(int Layer1, int NinLayer1, int Layer2, int NinLayer2,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords);
  double hkl_CartCoord_Length(xvector<double>& hkl_CartCoord,const xvector<double>& hkl);
} // namespace slab

// ----------------------------------------------------------------------------
// aflow_defects.cpp
class acage {
 public:
  // constructors/destructors                                   // --------------------------------------
  acage();                                                      // constructor default
  acage(const acage& b);                                        // constructor copy
  ~acage();                                                     // destructor
  // CONTENT                                                    // --------------------------------------
  void clear(void);                                             // clear everything
  xvector<double> origin_fpos;// default 3D
  xvector<double> origin_cpos;// default 3D
  double radius;
  uint coordination_position;
  int cages_position;
  int cages_irrtype;
  deque<_atom> atoms;
 private:                                                        // ---------------------------------------
  void free();                                                  // to free everything
  void copy(const acage& b);                                    // the flag is necessary because sometimes you need to allocate the space.
};
class _isort_acage_radius {                   // sorting through reference
 public:
  bool operator()(const acage& a1, const acage& a2) const {
    return (bool) (a1.radius>a2.radius);}
};

bool GetSphereFromFourPoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& v4);
bool GetCircumCircleeFromThreePoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
bool GetCircleFromTwoPoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2);
bool FPositionInsideCell(const xvector<double>& r);
bool EmptySphere(const deque<_atom>& grid_atoms,const xvector<double>& origin_cpos,const double& radius);
bool EmptySphere(const xstructure& str,const xvector<double>& origin_cpos,const double& radius);
uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax);
uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin);
uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax);
uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin);
bool AddCageToCages(const xstructure& str,const xvector<double>& origin_cpos,const xvector<double>& origin_fpos,const double& radius,
		    const int& cage_points_type,const double& roughness,
		    vector<acage>& cages,vector<acage>& cagesX,
		    const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2,int ithread);
uint GetCages4(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages4,
               const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2);
uint GetCages3(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages3,
               const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2);
uint GetCages2(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages2,
               const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2);
bool GetCages(const xstructure& _str,_aflags& aflags,
	      vector<acage>& cagesirreducible,vector<acage>& cagesreducible,vector<acage>& cages4,
	      vector<acage>& cages3,vector<acage>& cages2,const double& _roughness,const bool& osswrite,ostream& oss);
// ----------------------------------------------------------------------------
// aflow_pocc // CO 180502
namespace KBIN {
  void VASP_RunPOCC(const _xvasp& xvasp,const string& AflowIn,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE);
}
// ----------------------------------------------------------------------------
// aflow_phonons.cpp
namespace KBIN {
  void VASP_RunPhonons_APL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  void RunPhonons_APL(_xinput &xinput,string AflowIn,_aflags &aflags,_kflags &kflags,_xflags &xflags,ofstream &FileMESSAGE);  //now it's general
  // [OBSOLETE] bool PHON_RunPhonons(const xstructure& _str,_aflags& aflags,const double& radius,const bool& osswrite,ostream& oss);
  // ----------------------------------------------------------------------------
  // aflow_agl_debye.cpp
  void VASP_RunPhonons_AGL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  // ----------------------------------------------------------------------------
  // aflow_ael_elasticity.cpp
  void VASP_RunPhonons_AEL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------
// aconvasp_aflow.cpp

namespace pflow {
  // Dane Morgan style,adjusted by Stefano Curtarolo
  // Bringing definitions inside the template helps
  // constructing the right templates.
  template<class utype>
    class matrix {
  public:
    // constructors
    matrix(void) {};
    matrix(const int m);
    matrix(const int m,const int n);
    matrix(const int m,const int n,const utype& inutype);
    matrix(const int m,const vector<utype>& inutypevec);
    ~matrix(void) {};                         // destructor
    // accessors
    void print(void);
    uint size(void) const {
      return (uint) mat.size();}
    matrix<utype> transpose(void) const;
    //  matrix<utype>::iterator begin();
    //  matrix<utype>::iterator end();
    // operator
    vector<utype>& operator[] (const int i) {
      assert(i>=0 && i<=(int) mat.size()); return mat[i];}
    const vector<utype>& operator[] (const int i) const {
      assert(i>=0 && i<=(int) mat.size()); return mat[i];};
    const matrix<utype>& operator=(const matrix<utype>& b);
    // mutators
    void push_back(const vector<utype>& inutypevec) {
      mat.push_back(inutypevec);}
    void pop_back(void) {mat.pop_back();}
    void vecvec2mat(const vector<vector<utype> >& inVV);
    void vec2mat(const vector<utype>& inV);
    void clear(void) {mat.clear();}
    void insert(const int& id,const vector<utype>& inV) {
      mat.insert(mat.begin()+id,inV);}
    void erase(const int id);
    void erase_col(const int id);
  private:
    vector<vector<utype> > mat;
  };
  template <class utype> matrix<utype>  xmatrix2matrix(const xmatrix<utype>& );
  template <class utype> xmatrix<utype> matrix2xmatrix(const matrix<utype>& );
}


xstructure PutInCell(const xstructure& a);     // Bring all atoms in the cell (will be moved to external function)
xstructure PutInCompact(const xstructure& a);  // Bring all atoms in a compact shape (will be moved to external function)
xstructure GetPrim(const xstructure& a);
bool IsTranslationFVector(const xstructure& a,const xvector<double>& ftvec);
bool IsTranslationCVector(const xstructure& a,const xvector<double>& ctvec);
xvector<double> GetMom1(const xstructure& a);          // get moment_1 position of the atoms
xstructure SetMom1(const xstructure& a,const xvector<double>& mom1_in);     // set moment_1 position of atoms
xvector<double> AtomCDisp(const _atom& at1,const _atom& at2);
double AtomDist(const xstructure& str,const _atom& atom1,const _atom& atom2);  // with structure
double AtomDist(const _atom& at1, const _atom& at2);  // without structure
xvector<double> GetCDispFromOrigin(const _atom& atom);
double GetDistFromOrigin(const _atom& atom);
//void GetUnitCellRep(const xvector<double>& ppos,xvector<double>& p_cell0,xvector<int>& ijk,const xmatrix<double>& lattice,const bool coord_flag);
_atom ConvertAtomToLat(const _atom& in_at,const xmatrix<double>& lattice);
double GetXrayScattFactor(const string& name,const double& lambda);
xmatrix<double> RecipLat(const xmatrix<double>& lat);
double Normal(const double& x,const double& mu,const double& sigma);
xstructure SetLat(const xstructure& a,const xmatrix<double>& in_lat);
xmatrix<double> GetLat(const xstructure& a);

namespace pflow {
  // [OBSOLETE] void GetNeighData(const deque<_atom>& in_atom_vec,const xstructure& in_str,const double& rmin,const double& rmax,deque<deque<_atom> >& neigh_mat);
  // [OBSOLETE] void GetStrNeighData(const xstructure& str,const double cutoff,deque<deque<_atom> >& neigh_mat);
  double GetVol(const xmatrix<double>& lat);
  double GetVol(const pflow::matrix<double>& lat);
  double GetSignedVol(const xmatrix<double>& lat);
  double GetSignedVol(const pflow::matrix<double>& lat);
  xmatrix<double> RecipLat(const xmatrix<double>& lat);
  pflow::matrix<double> RecipLat(const pflow::matrix<double>& lat);
  _atom SetCpos(const _atom& a,const vector<double>& in_cpos);
  _atom SetFpos(const _atom& a,const vector<double>& in_fpos);
  vector<double> vecF2C(const pflow::matrix<double>& lat,const vector<double>& vf);
  vector<double> vecC2F(const pflow::matrix<double>& lat,const vector<double>& vc);
  _atom SetName(const _atom& a,const string& in_name);
  _atom SetType(const _atom& a,const int in_type);
  _atom SetNum(const _atom& a,const int in_num);
  vector<int> GetTypes(const xstructure& a);
  vector<string> GetNames(const xstructure& a);
  vector<string> GetCleanNames(const xstructure& a);
  vector<double> GetSpins(const xstructure& a);
  pflow::matrix<double> GetFpos(const xstructure& str);
  pflow::matrix<double> GetCpos(const xstructure& str);
  xstructure SetNumEachType(const xstructure& a,const deque<int>& in_num_each_type);
  deque<int> GetNumEachType(const xstructure& a);
  xstructure SetLat(const xstructure& a,const pflow::matrix<double>& in_lat);
  pflow::matrix<double> GetLat(const xstructure& a);
  double GetScale(const xstructure& a);
  pflow::matrix<double> GetScaledLat(const xstructure& a);
  xstructure AddAllAtomPos(const xstructure& a,const pflow::matrix<double>& in_pos,const int in_coord_flag);
  xstructure SetAllAtomPos(const xstructure& a,const pflow::matrix<double>& in_pos,const int in_coord_flag);
  xstructure SetAllAtomNames(const xstructure& a,const vector<string>& in_names);
  xstructure SetNamesWereGiven(const xstructure& a,const vector<int>& in_names_were_given);
  xstructure SetOrigin(const xstructure& a,const vector<double>& in_origin);
  xstructure SetOrigin(const xstructure& a,const xvector<double>& in_origin);
  bool VVequal(const vector<double>& a,const vector<double>& b);
  bool VVequal(const vector<int>& a,const vector<int>& b);
  bool VVequal(const deque<double>& a,const deque<double>& b);
  bool VVequal(const deque<int>& a,const deque<int>& b);
  vector<double> SmoothFunc(const vector<double>& func,const double& sigma);
  void VVset(matrix<double>& mat,const double& value);
  void VVset(vector<vector< int> >& mat,const int& value);
  double norm(const vector<double>& v);
  double getcos(const vector<double>& a,const vector<double>& b);
  //  vector<double> Getabc_angles(const pflow::matrix<double>& lat);   // confuses namespace
  vector<double> Sort_abc_angles(const vector<double>& abc_angles);
  void Vout(const vector<double>& a,ostream& out);
  void Vout(const vector<int>& a,ostream& out);
  void Vout(const vector<string>& a,ostream& out);
  void Mout(const pflow::matrix<double>& m,ostream& out);
  void Mout(const vector<vector<double> >& m,ostream& out);
  vector<double> SVprod(const double& s,const vector<double>& b);
  vector<int> SVprod(const int& s,const vector<int>& b);
  vector<double> VVsum(const vector<double>& a,const vector<double>& b);
  vector<double> VVsum(const vector<double>& a,const vector<int>& b);
  vector<double> VVdiff(const vector<double>& a,const vector<double>& b);
  double VVprod(const vector<double>& a,const vector<double>& b);
  double VVprod(const vector<double>& a,const vector<int>& b);
  pflow::matrix<double> MMmult(const pflow::matrix<double>& a,const pflow::matrix<double>& b);
  vector<double> MVmult(const pflow::matrix<double>& A,const vector<double>& v);
  vector<double> VMmult(const vector<double>& v,const pflow::matrix<double>& A);
  vector<double> VMmult(const vector<int>& v,const pflow::matrix<double>& A);
  vector<double> VVcross(const vector<double>& a,const vector<double>& b);
  double VVdot(const vector<double>& a,const vector<double>& b);
  int GetNumAtoms(const xstructure& a);
  void SetSpline(const vector<double>& x,const vector<double>& y,const double& yp1,const double& ypn,vector<double>& y2);
  void GetSplineInt(const vector<double>& xa,const vector<double>& ya,vector<double>& y2a,const double& x,double& y);
  void PrintSpline(const vector<double>& x,const vector<double>& y,const int& npts,ostream& outf);
}

// ----------------------------------------------------------------------------
// aflow_wyckoff.cpp

xvector<double> wv(const double &x,const double &y,const double &z);
void wa(_atom& a,xstructure &str);
xstructure WyckoffPOSITIONS(uint spacegroup, xstructure strin);
xstructure WyckoffPOSITIONS(uint spacegroup, uint option, xstructure strin);

// ----------------------------------------------------------------------------
// aflow_apennsy stuff
#include "aflow_apennsy.h"

// ----------------------------------------------------------------------------
// aflowlib.h stuff
#include "aflowlib.h"

// ----------------------------------------------------------------------------
// aflowlib.h stuff
#include "aflowlib.h"

// ----------------------------------------------------------------------------
// aflowlib.h stuff
#include "aflow_pflow.h"

#endif
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

