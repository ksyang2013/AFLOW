// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

#define __XLIBS_LINK
#include "aflow.h"
#define cdebug cerr
#include <algorithm>

//bool nocase_compare(char c1,char c2) {return toupper(c1)==toupper(c2);}

//#define MaxAflowInSize 65535
//string AflowIn; //[MaxAflowInSize];

#define VRUNS_MAX_CUTOFF 32768
#define DUKE_BETANEW_DEFAULT_KILL_MEM_CUTOFF 1.50
#define DUKE_QRATS_DEFAULT_KILL_MEM_CUTOFF 1.50
#define DUKE_QUSER_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MPCDF_EOS_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MPCDF_DRACO_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MPCDF_HYDRA_DEFAULT_KILL_MEM_CUTOFF 1.50

namespace aurostd {
  // ***************************************************************************
  // Function DirectoryAlreadyInDabatase
  // ***************************************************************************
  bool DirectoryAlreadyInDatabase(string directory,bool FORCE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(FORCE) return FALSE; // neglect already in the database

    // already scanned
    if(FileExist(directory+"/ALREADY_IN_DATABASE")) return TRUE;
    if(FileExist(directory+"/ALREADY_IN_DATABASE.gz")) return TRUE;
    if(FileExist(directory+"/ALREADY_IN_DATABASE.bz2")) return TRUE;
  
    // no, then scan
    if(LDEBUG) cerr << "SEARCHING" << endl;
    bool already_in_database=FALSE;
    string tmp_directory=directory;
    aurostd::StringSubst(tmp_directory,_AFLOWIN_," ");
    aurostd::StringSubst(tmp_directory,"./"," ");
    aurostd::StringSubst(tmp_directory,"/"," ");
    vector<string> tokens;
    aurostd::string2tokens(tmp_directory,tokens);
    if(tokens.size()>=2) tmp_directory=tokens.at(tokens.size()-2)+"/"+tokens.at(tokens.size()-1);
    if(tokens.size()==1) tmp_directory=tokens.at(tokens.size()-1);

    uint library=LIBRARY_NOTHING;
    // XHOST_LIBRARY_LIB1
    if(aurostd::substring2bool(directory,"LIB1"))     library=XHOST_LIBRARY_LIB1;
    if(aurostd::substring2bool(directory,"AURO"))     library=XHOST_LIBRARY_LIB1;    // [HISTORIC]
    // XHOST_LIBRARY_LIB2
    if(aurostd::substring2bool(directory,"LIBRARYU")) library=XHOST_LIBRARY_LIB2;
    if(aurostd::substring2bool(directory,"LIBRARYX")) library=XHOST_LIBRARY_LIB2;    // [HISTORIC]
    if(aurostd::substring2bool(directory,"LIB2"))     library=XHOST_LIBRARY_LIB2;
    // XHOST_LIBRARY_ICDS
    if(aurostd::substring2bool(directory,"ICSD"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"SCINT"))    library=XHOST_LIBRARY_ICSD;    // [HISTORIC]
    if(aurostd::substring2bool(directory,"BCC"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"BCT"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"CUB"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"FCC"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"HEX"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"MCL"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"MCLC"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORC"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORCC"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORCF"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORCI"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"RHL"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"TET"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"TRI"))      library=XHOST_LIBRARY_ICSD;
    // XHOST_LIBRARY_ICSD
    if(aurostd::substring2bool(directory,"MAGNETIC")) library=XHOST_LIBRARY_LIB3;    // [HISTORIC]
    if(aurostd::substring2bool(directory,"LIB3"))     library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0001"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0002"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0001"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0004"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0005"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0006"))    library=XHOST_LIBRARY_LIB3;
    // XHOST_LIBRARY_LIB4
    if(aurostd::substring2bool(directory,"LIB4"))     library=XHOST_LIBRARY_LIB4;
    if(aurostd::substring2bool(directory,"Q0001"))    library=XHOST_LIBRARY_LIB4;
    // XHOST_LIBRARY_LIB5
    if(aurostd::substring2bool(directory,"LIB5"))     library=XHOST_LIBRARY_LIB5;
    if(aurostd::substring2bool(directory,"P0001"))    library=XHOST_LIBRARY_LIB5;
    // XHOST_LIBRARY_LIB6
    if(aurostd::substring2bool(directory,"LIB6"))     library=XHOST_LIBRARY_LIB6;
    if(aurostd::substring2bool(directory,"H0001"))    library=XHOST_LIBRARY_LIB6;
    // XHOST_LIBRARY_LIB7
    if(aurostd::substring2bool(directory,"LIB7"))     library=XHOST_LIBRARY_LIB7;
    // XHOST_LIBRARY_LIB8
    if(aurostd::substring2bool(directory,"LIB8"))     library=XHOST_LIBRARY_LIB8;
    // XHOST_LIBRARY_LIB9
    if(aurostd::substring2bool(directory,"LIB9"))     library=XHOST_LIBRARY_LIB9;

    // found something
    if(library!=LIBRARY_NOTHING) {
      tokens.clear();
      string tmp;
      vector<string> vLibrary;
      if(library==XHOST_LIBRARY_LIB1) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB1" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB1_RAW"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB2) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB2" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB2_RAW"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB3) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB3" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB3_LIB"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB4) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB4" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB4_LIB"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB5) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB5" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB5_LIB"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB6) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB6" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB6_LIB"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB7) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB7" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB7_LIB"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB8) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB8" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB8_LIB"),tokens,"\n");
      }
      if(library==XHOST_LIBRARY_LIB9) {
	if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB9" << endl;
	aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB9_LIB"),tokens,"\n");
      }
      for(uint i=0;i<tokens.size();i++)
	if(aurostd::substring2bool(tokens.at(i),"/")) {
	  tmp=tokens.at(i);
	  aurostd::StringSubst(tmp," ","");
	  vLibrary.push_back(tmp);
	}
      for(int i=vLibrary.size()-1;i>=0;i--)
	if(tmp_directory==vLibrary.at(i)) {
	  already_in_database=TRUE;
	  if(LDEBUG) cerr << vLibrary.at(i) << " FOUND .." << endl; // NEW
	}
    }
    
    if(LDEBUG) cerr << "DONE..." << endl;
    if(LDEBUG) cerr << tmp_directory << endl;

    if(already_in_database) {
      //    cerr << directory << " already in database" << endl;
      aurostd::execute("cp -f "+directory+"/"+_AFLOWIN_+" "+directory+"/ALREADY_IN_DATABASE");
      aurostd::execute(XHOST.command("bzip2")+" -f "+directory+"/"+_AFLOWIN_);
    }

    return already_in_database;
  }
}

using aurostd::DirectorySkipped;
using aurostd::DirectoryAlreadyInDatabase;
using aurostd::DirectoryUnwritable;

// int KBIN_MODE;

//#define KBIN_VOID_MODE 0             // just a shift
//#define KBIN_VASP_MODE 2             // for ab-initio VASP mode
//#define KBIN_XXXX_MODE 3             // for XXXX program
//#define KBIN_GRND_MODE 4             // for classical monte carlo

// GND MODE
#define KBIN_VASP_N_VPARS 32
#define _KBIN_VASP_SLEEP_ 2
#define _KBIN_LOOP_SLEEP_ 300
// PRIORITY
#define PRIORITY_PROBABILITY 0.2000
//#define PRIORITY_GREP_STRING string("grep -vi xxxxx ")
#define PRIORITY_GREP_STRING string("grep -vi PRIORITY ")

// ***************************************************************************
// KBIN::Legitimate_aflowin
// ***************************************************************************
namespace KBIN {
  bool Legitimate_aflowin(string _aflowindir,const bool& osswrite,ostringstream& oss) {
    string aflowindir=_aflowindir;
    aurostd::StringSubst(aflowindir,"//","/");
  
    if(aurostd::FileExist(aflowindir)) {  // file must exist
      if(!aurostd::FileEmpty(aflowindir)) { // must not be empty
	if(aurostd::substring2bool(aflowindir,_AFLOWIN_)) {  // there must be an _AFLOWIN_
	  aurostd::StringSubst(aflowindir,_AFLOWIN_,"");
	  if(!aurostd::FileExist(aflowindir+_AFLOWLOCK_)) { // it should be UNLOCKED
	    //  if(osswrite) {oss << "MMMMM  Loading Valid File Entry = " << aflowindir << MessageTime(aflags);aurostd::PrintMessageStream(oss,XHOST.QUIET);};
	    return TRUE;	
	  } else { // must be unlocked
	    if(osswrite) {oss << "MMMMM  Directory locked = " << aflowindir << Message("time") << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
	    return FALSE;
	  }
	} else { // must contain _AFLOWIN_
	  if(osswrite) {oss << "MMMMM  Not loading file without " << _AFLOWIN_ << " = " << aflowindir << Message("time") << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
	  return FALSE;
	}
      } else { // empty
	if(osswrite) {oss << "MMMMM  Not loading empty file = " << aflowindir << Message("time") << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
	return FALSE;
      }
    } else { // unexisting
      // if(osswrite) {oss << "MMMMM  Not loading unexisting file = " << aflowindir << Message("time") << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
      return FALSE;
    }
    return FALSE;
  }
}

namespace KBIN {
  bool Legitimate_aflowin(string aflowindir) { ostringstream aus; return KBIN::Legitimate_aflowin(aflowindir,FALSE,aus);};
}

// ***************************************************************************
// KBIN::Main
// ***************************************************************************
namespace KBIN {
  int KBIN_Main(vector<string> argv) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string GENERIC;
    //  string Directory;
    int i;
    ostringstream aus;
    ifstream FileAUS;
    _aflags aflags;
    aurostd::StringstreamClean(aus);
    bool _VERBOSE_=FALSE;
    int XHOST_AFLOW_RUNXnumber_multiplier=3;

    std::deque<_aflags> qaflags;

    AFLOW_PTHREADS::Clean_Threads();                                    // clean threads
    // _aflags taflags[MAX_ALLOCATABLE_PTHREADS];
    // _threaded_KBIN_params params[MAX_ALLOCATABLE_PTHREADS];
  
    // cerr << "GMODE" << endl;
    // check BlackList **************************************************
    if(AFLOW_BlackList(XHOST.hostname)) {
      aus << "MMMMM  HOSTNAME BLACKLISTED = " << XHOST.hostname << " - " << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      return 0;
    }
    // get KBIN **************************************************
    //  aus << "MMMMM  AFLOW: running " << _AFLOWIN_ << ": " << " "  << " - " << Message(aflags,"user,host,time") << endl;   // too much verbosity is annoying
    //  aurostd::PrintMessageStream(aus,XHOST.QUIET);                                                  // too much verbosity is annoying
    
    // do some running
    if(XHOST.AFLOW_RUNXflag) {
      vector<string> tokens;
      if(aurostd::args2attachedflag(argv,"--run=")) {
	aurostd::string2tokens(aurostd::args2attachedstring(argv,"--run=","1"),tokens,"=");
	XHOST.AFLOW_RUNXnumber=aurostd::string2utype<uint>(tokens.at(tokens.size()-1));
      }
      if(aurostd::args2attachedflag(argv,"-run=")) {
	aurostd::string2tokens(aurostd::args2attachedstring(argv,"-run=","1"),tokens,"=");
	XHOST.AFLOW_RUNXnumber=aurostd::string2utype<uint>(tokens.at(tokens.size()-1));
      }
      if(XHOST.AFLOW_RUNXnumber<1) XHOST.AFLOW_RUNXnumber=1;
    }
    
    if(aurostd::args2flag(argv,"-runone|--runone|-run_one|--run_one|-run1|--run1|--run=1|-run=1")) { // RUNONE COMPATIBILITY
      XHOST.AFLOW_RUNXflag=TRUE;XHOST.AFLOW_RUNXnumber=1;
    }
    
    if(XHOST.AFLOW_RUNXflag)  {
      XHOST.AFLOW_RUNDIRflag=FALSE;
      XHOST.AFLOW_MULTIflag=FALSE;
    }
    if(XHOST.AFLOW_MULTIflag)  {
      XHOST.AFLOW_RUNDIRflag=FALSE;
      XHOST.AFLOW_RUNXflag=FALSE;
    }
    
    if(XHOST.vflag_aflow.flag("LOOP")) {aus << "MMMMM  KBIN option XHOST.vflag_aflow.flag(\"LOOP\") - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_RUNDIRflag) {aus << "MMMMM  KBIN option [--run] (XHOST.AFLOW_RUNDIRflag=TRUE) - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_MULTIflag) {aus << "MMMMM  KBIN option [--run=multi] (XHOST.AFLOW_MULTIflag=TRUE) - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_RUNXflag && XHOST.AFLOW_RUNXnumber==1) {aus << "MMMMM  KBIN option [--run=1] (XHOST.AFLOW_RUNXflag=TRUE, XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_RUNXflag && XHOST.AFLOW_RUNXnumber>1) {aus << "MMMMM  KBIN option [--run=N] (XHOST.AFLOW_RUNXflag=TRUE, XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}

    aflags.KBIN_RUN_AFLOWIN=TRUE;
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN=FALSE;
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP=FALSE;
    //DX
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN=FALSE;
    //DX
    aflags.KBIN_DELETE_AFLOWIN=FALSE;

    aflags.KBIN_GEN_VASP_FROM_AFLOWIN        = aurostd::args2flag(argv,"--generate_vasp_from_aflowin|--generate");
    if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
      aflags.KBIN_RUN_AFLOWIN=FALSE;
    }
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP        = aurostd::args2flag(argv,"--generate_aflowin_from_vasp");

    //DX
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN      = aurostd::args2flag(argv,"--generate_symmetry|--generate_sym");
    //DX

    //DX
    //DX if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN || aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
    if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN || aflags.KBIN_GEN_AFLOWIN_FROM_VASP || aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN) {
      //DX
      XHOST.AFLOW_RUNXflag=TRUE;XHOST.AFLOW_RUNXnumber=1;
      XHOST.AFLOW_RUNDIRflag=TRUE;
      if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP) aflags.KBIN_RUN_AFLOWIN=FALSE;
    }
    
    aflags.KBIN_DELETE_AFLOWIN=aurostd::args2flag(argv,"--delete_aflowin");

    aflags.AFLOW_MODE_QSUB_MODE1 = aurostd::args2flag(argv,"--qsub1|-qsub1");
    aflags.AFLOW_MODE_QSUB_MODE2 = aurostd::args2flag(argv,"--qsub2|-qsub2");
    aflags.AFLOW_MODE_QSUB_MODE3 = aurostd::args2flag(argv,"--qsub3|-qsub3");

    // [OBSOLETE]  MPI=aurostd::args2flag(argv,"--MPI|--mpi");  // ABSOLUTELY otherwise the multithreads kicks
    aflags.AFLOW_FORCE_MPI = aurostd::args2flag(argv,"--MPI|--mpi");
    aflags.AFLOW_FORCE_SERIAL = aurostd::args2flag(argv,"--nompi|-nompi|--serial|-serial");
    // [OBSOLETE]  aflags.AFLOW_GLOBAL_NCPUS=aurostd::args2utype(argv,"--np",(int) 0);
    aflags.AFLOW_GLOBAL_NCPUS=aurostd::args2attachedutype<int>(argv,"--np=",(int) 0);
    // if(aflags.AFLOW_GLOBAL_NCPUS && !MPI && !aflags.AFLOW_FORCE_MPI) {}
    if(XHOST.MPI || aflags.AFLOW_FORCE_MPI) AFLOW_PTHREADS::No_Threads();

    aflags.AFLOW_MACHINE_GLOBAL.clear();
    // "MACHINE::DUKE_BETA_MPICH"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH",
				     aurostd::args2flag(argv,"--machine=beta|--machine=duke_beta|--beta|--duke_beta|--machine=beta_mpich|--machine=duke_beta_mpich"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH")) XHOST.maxmem=DUKE_BETANEW_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::DUKE_BETA_OPENMPI"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI",aurostd::args2flag(argv,"--machine=beta_openmpi|--machine=duke_beta_openmpi"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI")) XHOST.maxmem=DUKE_BETANEW_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::DUKE_QRATS_MPICH"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH",aurostd::args2flag(argv,"--machine=qrats|--machine=duke_qrats|--machine=qrats_mpich|--machine=duke_qrats_mpich"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH")) XHOST.maxmem=DUKE_QRATS_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::DUKE_QUSER_OPENMPI"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QUSER_OPENMPI",aurostd::args2flag(argv,"--machine=quser|--machine=duke_quser|--machine=quser_openmpi|--machine=duke_quser_openmpi"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QUSER_OPENMPI")) XHOST.maxmem=DUKE_QUSER_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::MPCDF_EOS_MPIIFORT"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS_MPIIFORT",aurostd::args2flag(argv,"--machine=eos|--machine=mpcdf_eos|--machine=eos_mpiifort|--machine=mpcdf_eos_mpiifort"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS_MPIIFORT")) XHOST.maxmem=MPCDF_EOS_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::MPCDF_DRACO_MPIIFORT"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO_MPIIFORT",aurostd::args2flag(argv,"--machine=draco|--machine=mpcdf_draco|--machine=draco_mpiifort|--machine=mpcdf_draco_mpiifort"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO_MPIIFORT")) XHOST.maxmem=MPCDF_DRACO_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::MPCDF_HYDRA_MPIIFORT"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA_MPIIFORT",aurostd::args2flag(argv,"--machine=hydra|--machine=mpcdf_hydra|--machine=hydra_mpiifort|--machine=mpcdf_hydra_mpiifort"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA_MPIIFORT")) XHOST.maxmem=MPCDF_HYDRA_DEFAULT_KILL_MEM_CUTOFF;
    // DUKE_MATERIALS
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_MATERIALS",aurostd::args2flag(argv,"--machine=materials|--machine=duke_materials"));
    // DUKE_AFLOWLIB
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_AFLOWLIB",aurostd::args2flag(argv,"--machine=aflowlib|--machine=duke_aflowlib"));
    // DUKE_HABANA
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_HABANA",aurostd::args2flag(argv,"--machine=habana|--machine=duke_habana"));
    // TERAGRID_RANGER
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::TERAGRID_RANGER",aurostd::args2flag(argv,"--machine=ranger|--machine=teragrid_ranger"));
    // TERAGRID_KRAKEN
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::TERAGRID_KRAKEN",aurostd::args2flag(argv,"--machine=kraken|--machine=teragrid_kraken"));
    // TRINITY_PARSONS
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::TRINITY_PARSONS",aurostd::args2flag(argv,"--machine=parsons|--machine=trinity_parsons|--machine=kelvin|--machine=trinity_kelvin|--machine=trinity"));
    // FULTON_MARYLOU
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::FULTON_MARYLOU",aurostd::args2flag(argv,"--machine=marylou|--machine=fulton_marylou"));
    // MACHINE2
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE2",aurostd::args2flag(argv,"--machine=ohad"));   
    // MACHINE1
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE1",aurostd::args2flag(argv,"--machine=host1"));
  
  
    // turn off multi pthreads on specific machines
    if(aflags.AFLOW_MACHINE_GLOBAL.flag()) {
      AFLOW_PTHREADS::No_Threads();
      AFLOW_PTHREADS::FLAG=FALSE; // safety...
      AFLOW_PTHREADS::MAX_PTHREADS=1; // safety...
      //    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
      XHOST.MPI=TRUE;
    }
    
    vector<string> vruns;
  
    aflags.AFLOW_PERFORM_CLEAN=XHOST.vflag_aflow.flag("CLEAN");// || XHOST.vflag_aflow.flag("XCLEAN"));
    // [OBSOLETE] aflags.AFLOW_PERFORM_DIRECTORY=aurostd::args2flag(argv,"--DIRECTORY|--D|--d|./");
    aflags.AFLOW_PERFORM_DIRECTORY=XHOST.vflag_control.flag("VDIR");
    // [OBSOLETE] aflags.AFLOW_PERFORM_FILE=aurostd::args2flag(argv,"--FILE|--F|--f");
    aflags.AFLOW_PERFORM_FILE=XHOST.vflag_control.flag("FILE");
    aflags.AFLOW_PERFORM_ORDER_SORT=aurostd::args2flag(argv,"--sort|-sort");                    // Sorts the _AFLOWIN_ in the list
    aflags.AFLOW_PERFORM_ORDER_REVERSE=aurostd::args2flag(argv,"--reverse|--rsort|-reverse|-rsort"); // Reverse the _AFLOWIN_ in the list
    aflags.AFLOW_PERFORM_ORDER_RANDOM=aurostd::args2flag(argv,"--random|--rnd|-random|-rnd"); // Randomize the _AFLOWIN_ in the list
    aflags.AFLOW_FORCE_RUN = aurostd::args2flag(argv,"--force|-force");

    if(aflags.AFLOW_PERFORM_DIRECTORY) {aus << "MMMMM  KBIN option PERFORM_DIRECTORY - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_FILE) {aus << "MMMMM  KBIN option PERFORM_FILE - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_ORDER_SORT) {aus << "MMMMM  KBIN option ORDER_SORT - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_ORDER_REVERSE) {aus << "MMMMM  KBIN option ORDER_REVERSE - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_ORDER_RANDOM) {aus << "MMMMM  KBIN option ORDER_RANDOM - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_FORCE_RUN) {aus << "MMMMM  KBIN option FORCE_RUN - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}

    uint maxcheck=VRUNS_MAX_CUTOFF;
    if(XHOST.AFLOW_MULTIflag) maxcheck=VRUNS_MAX_CUTOFF;
    if(XHOST.AFLOW_RUNXflag) maxcheck=XHOST.AFLOW_RUNXnumber;
    if(maxcheck==0) maxcheck=1; // safety

    // simple commands
    if(XHOST.vflag_aflow.flag("XCLEAN")) {
      KBIN::XClean(XHOST.vflag_aflow.getattachedscheme("XCLEAN"));
      return 1;
    }

    // for directory mode load them all
    if(aflags.AFLOW_PERFORM_DIRECTORY) {
      // [OBSOLETE] vruns=aurostd::args2vectorstring(argv,"--DIRECTORY|--D|--d","./");
      aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VDIR"),vruns,",");
      if(LDEBUG) { for(uint i=0;i<vruns.size();i++) cerr << "KBIN::Main: vruns.at(i)=" << vruns.at(i) << endl;}
    } else {
      if(!aflags.AFLOW_PERFORM_FILE)
	vruns.push_back(aurostd::execute2string(XHOST.command("pwd")));
      // vruns.push_back(aurostd::execute2string(XHOST.command("pwd"))+" ./");
    }
    // if file found
    if(aflags.AFLOW_PERFORM_FILE) {
      // [OBSOLETE]   string file_name=aurostd::args2string(argv,"--FILE|--F|--f","xxxx");
      string file_name=XHOST.vflag_control.getattachedscheme("FILE");
      if(!aurostd::FileExist(file_name)) {
	aus << "EEEEE  FILE_NOT_FOUND = " << file_name  << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);
      }
      if(aurostd::FileEmpty(file_name)) {
	aus << "EEEEE  FILE_EMPTY = " << file_name  << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);
      }
      aus << "MMMMM  Loading File = " << file_name << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      vector<string> vlines;
      vlines.clear();
      aurostd::file2vectorstring(file_name,vlines);
    
      aus << "MMMMM  " <<  aurostd::PaddedPOST("Legitimate VLINES = "+aurostd::utype2string(vlines.size()),40) << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      aus << "MMMMM  " <<  aurostd::PaddedPOST("         maxcheck = "+aurostd::utype2string(maxcheck),40) << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
    
      if(aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
	aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	aurostd::sort(vlines);}
      if(aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
	aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	aurostd::rsort(vlines);}
      if(aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
	aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] start (XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(aus,XHOST.QUIET);
	aurostd::random_shuffle(vlines);  // uses the std library but the seed is initialized in xrandom too
 	aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] stop (XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(aus,XHOST.QUIET);
      }
      
      for(uint i=0;(i<vlines.size() && vruns.size()<XHOST_AFLOW_RUNXnumber_multiplier*maxcheck && vruns.size()<VRUNS_MAX_CUTOFF);i++) {  // XHOST_AFLOW_RUNXnumber_multiplier times more... for safety
	if(KBIN::Legitimate_aflowin(vlines.at(i),FALSE,aus)) vruns.push_back(vlines.at(i));             // TRUE puts too much verbosity
      }
    
      aus << "MMMMM  " <<  aurostd::PaddedPOST("Legitimate VRUNS = "+aurostd::utype2string(vruns.size()),40) << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      // aus << "MMMMM  Legitimate VRUNS = " << vruns.size() << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      // cerr << vlines.size() << endl;
      // exit(0);
    }
  
    if(aflags.AFLOW_PERFORM_CLEAN && !XHOST.AFLOW_RUNDIRflag && !XHOST.AFLOW_MULTIflag && !XHOST.AFLOW_RUNXflag) {
      XHOST.AFLOW_RUNDIRflag=TRUE; // give something to clean
    }

    if(aflags.AFLOW_PERFORM_CLEAN && !aflags.AFLOW_PERFORM_DIRECTORY) {
      cerr << "AFLOW: to use --clean, you must specify one or more directories" << endl;
      exit(0);
    }
    if(aflags.AFLOW_PERFORM_CLEAN && aflags.AFLOW_PERFORM_DIRECTORY) {
      //    cout << "DEBUG CLEAN = " << vruns.size() << endl;
    }
  
    bool STOP_DEBUG=aurostd::args2flag(argv,"--STOP|--stop");
  
    // cdebug << "aflags.KBIN_RUN_AFLOWIN=" << aflags.KBIN_RUN_AFLOWIN << endl;
    // cdebug << "aflags.KBIN_GEN_VASP_FROM_AFLOWIN=" << aflags.KBIN_GEN_VASP_FROM_AFLOWIN << endl;
    // cdebug << "aflags.KBIN_GEN_AFLOWIN_FROM_VASP=" << aflags.KBIN_GEN_AFLOWIN_FROM_VASP << endl;
  
    // ------------------------------------------------------------------------------------------------------------------------------------
    // nothing to run
    if(!XHOST.AFLOW_RUNDIRflag && !XHOST.AFLOW_MULTIflag && !XHOST.AFLOW_RUNXflag && !aflags.AFLOW_PERFORM_CLEAN && !aflags.AFLOW_PERFORM_DIRECTORY) {
      aus << "MMMMM  KBIN option nothing to run [--run , --run=multi, --run=1, --run=N] - " << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //    return 0;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------
    // nothing specified : CHECK IF DAEMON
    if(XHOST.AFLOW_RUNDIRflag) { // check if daemaon aflowd
      string progname=argv.at(0);
      if(aurostd::substring2bool(progname,"aflowd")) {
	XHOST.AFLOW_MULTIflag=TRUE;
	XHOST.AFLOW_RUNXflag=FALSE;
	XHOST.vflag_aflow.flag("LOOP",TRUE); // add automatically
	aus << "MMMMM  AFLOW: running as DAEMON (aflow --kmode --multi --loop): " << " - " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(aus,XHOST.QUIET);
      }
    }
  
    // ------------------------------------------------------------------------------------------------------------------------------------
    // run specified directory: XHOST.AFLOW_RUNDIRflag
    if(XHOST.AFLOW_RUNDIRflag) {
      //    bool krun=TRUE;
      if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP0b" << endl;
      // [OBSOLETE]   vector<string> vDirectory(aurostd::args2vectorstring(argv,"--DIRECTORY|--D|--d","./"));
      vector<string> vDirectory=vruns;
      // fix the RUNS
      for(uint ii=0;ii<vDirectory.size();ii++)
	if(aurostd::substring2bool(vDirectory.at(ii),_AFLOWIN_))
	  aurostd::StringSubst(vDirectory.at(ii),_AFLOWIN_,"");
    
      // if(aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
      //   aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::sort(vDirectory);}
      // if(aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
      //   aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::rsort(vDirectory);}
      // if(aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
      //   aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::random_shuffle(vDirectory);}  // uses the std library but the seed is initialized in xrandom too
    
      // aurostd::random_shuffle(vDirectory);
      // std::random_shuffle(vDirectory.begin(),vDirectory.end());
    
      for(uint idir=0;idir<vDirectory.size();idir++) {
	bool krun=TRUE;
	aflags.Directory=vDirectory.at(idir);
	aus << "MMMMM  AFLOW: running " << _AFLOWIN_ << ", directory" << "=\""  << aflags.Directory << "\" - " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(aus,XHOST.QUIET);
	// If necessary PERFORM CLEAN
	if(krun && aflags.AFLOW_PERFORM_CLEAN) {
	  aflags.Directory=vDirectory.at(idir);
	  aurostd::StringSubst(aflags.Directory,"/"+_AFLOWLOCK_,"");
	  aurostd::StringSubst(aflags.Directory,"/OUTCAR.relax1.bz2","/OUTCAR.relax2.bz2","/OUTCAR.static.bz2","/OUTCAR.bands.bz2","");    // CLEAN UP A LITTLE
	  aurostd::StringSubst(aflags.Directory,"/EIGENVAL.relax1.bz2","/EIGENVAL.relax2.bz2","/EIGENVAL.static.bz2","/EIGENVAL.bands.bz2","");  // CLEAN UP A LITTLE
	  aurostd::StringSubst(aflags.Directory,"/OUTCAR","");  // so it is easier to search
	  aurostd::StringSubst(aflags.Directory,"/"+_AFLOWIN_,"");  // so it is easier to search
	  //  cerr << aflags.Directory << endl;
	  KBIN::Clean(aflags);
	  krun=FALSE;
	}
	// RUN
	// cerr << "STEP0b" << endl;
	if(krun) {
	  if(LDEBUG) cerr << "STEP1b" << endl;  //CO 170622 - should be debug
	  ifstream FileCHECK;string FileNameCHECK;
	  // check for directory
	  FileNameCHECK=aflags.Directory;
	  FileCHECK.open(FileNameCHECK.c_str(),std::ios::in);
	  FileCHECK.clear();FileCHECK.close();
	  if(!FileCHECK) {                                                                        // ******* Directory is non existent
	    aus << "EEEEE  DIRECTORY_NOT_FOUND = "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  } else {                                                                                // ******* Directory EXISTS
	    if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP1c" << endl;
	    if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_)) {                                               // ******* Directory is locked
	      aus << "LLLLL  DIRECTORY_LOCKED ...bzzzz... !MULTI = "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(aus,XHOST.QUIET);
	    } else {
	      if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP1d" << endl;
	      if(DirectorySkipped(aflags.Directory)) {                                            // ******* Directory is skipped
		aus << "LLLLL  DIRECTORY_SKIPPED ...bzzzz... !MULTI = "  << Message(aflags,"user,host,time") << endl;
		aurostd::PrintMessageStream(aus,XHOST.QUIET);
	      } else {        
                if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP1e" << endl;
		if(DirectoryAlreadyInDatabase(aflags.Directory,aflags.AFLOW_FORCE_RUN)) {         // ******* Directory is already in the database
		  aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ...bzzzz... !MULTI = "  << Message(aflags,"user,host,time") << endl;
		  aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ... use \"aflow --multi\" to force the calculation of this entry "  << Message(aflags,"user,host,time") << endl;
		  aurostd::PrintMessageStream(aus,XHOST.QUIET);
		} else {        
                  if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP1f" << endl;
		  if(DirectoryUnwritable(aflags.Directory)) {                                     // ******* Directory is unwritable
		    aus << "LLLLL  DIRECTORY_UNWRITABLE ...bzzzz... !MULTI = "  << Message(aflags,"user,host,time") << endl;
		    aurostd::PrintMessageStream(aus,XHOST.QUIET);
		  } else {                                                                        // ******* Directory is ok
		    if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP1g" << endl;
		    if(_VERBOSE_) aus << "LLLLL  GOOD ...bzzz... !MULTI = "  << Message(aflags,"user,host,time") << endl;
		    if(_VERBOSE_) aurostd::PrintMessageStream(aus,XHOST.QUIET);
		    if(aflags.KBIN_RUN_AFLOWIN)               FileNameCHECK=aflags.Directory+"/"+_AFLOWIN_;
		    if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)     FileNameCHECK=aflags.Directory+"/"+_AFLOWIN_;
		    if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     FileNameCHECK=aflags.Directory+"/INCAR";
		    FileCHECK.open(FileNameCHECK.c_str(),std::ios::in);
		    FileCHECK.clear();FileCHECK.close();
		    if(!FileCHECK) {                                                                    // ******* _AFLOWIN_ does not exist
		      if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP1h" << endl;
		      if(aflags.KBIN_RUN_AFLOWIN)               aus << "EEEEE  " << _AFLOWIN_ << " not found  = "  << Message(aflags,"user,host,time") << endl;
		      if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)     aus << "EEEEE  " << _AFLOWIN_ << " not found  = "  << Message(aflags,"user,host,time") << endl;
		      if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     aus << "EEEEE  INCAR not found     = "  << Message(aflags,"user,host,time") << endl;
		      aurostd::PrintMessageStream(aus,XHOST.QUIET);
		    } else {                                                                            // ******* _AFLOWIN_ exists RUN
		      if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP1i" << endl;
		      if(aflags.KBIN_RUN_AFLOWIN)  {
			KBIN::RUN_Directory(aflags);
		      }
		      // if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)    KBIN::RUN_Directory(aflags);
		      if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     KBIN::GenerateAflowinFromVASPDirectory(aflags);
		      aus << "MMMMM  AFLOW: Run Done " << " - " << Message(aflags,"user,host,time") << endl;
		      aurostd::PrintMessageStream(aus,XHOST.QUIET);
		    }
		  } // DIRECTORY WRITABLE
		} // DIRECTORY NOT IN THE DATABASE
	      } // DIRECTORY UN SKIPPED
	    } // DIRECTORY UN LOCKED
	  } // DIRECTORY FOUND
	  aus << "MMMMM  AFLOW: Done " << " - " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	}
      } // idir
      // aus << "MMMMM  AFLOW: Done " << " - " << Message(aflags,"user,host,time") << endl;
      // aurostd::PrintMessageStream(aus,XHOST.QUIET);
      if(LDEBUG) cerr << "KBIN::KBIN_Main: STEP2" << endl;
    }
  
    // ERRORS ------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------------------------------------
    // run MULTI and XHOST.AFLOW_RUNXflag (in XHOST.AFLOW_RUNXflag, runs only XHOST.AFLOW_RUNXnumber and then dies)
    // MULTI with SINGLE AND MULTI THREAD VERSION -----------------------------------------------------------------------------------------
    if(XHOST.AFLOW_MULTIflag || XHOST.AFLOW_RUNXflag) {
      uint RUN_times=0;
      bool MULTI_DEBUG=FALSE;
      if(MULTI_DEBUG) {aus << "MMMMM  AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
      if(MULTI_DEBUG) {aus << "MMMMM  AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
      bool EMPTY=FALSE,FOUND=FALSE;
      bool free_thread;
      int ithread=0;
      if(AFLOW_PTHREADS::FLAG && AFLOW_PTHREADS::MAX_PTHREADS<=1) {
	aus << "EEEEE  ERROR PTHREADS" << endl;
	aus << "MMMMM  AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << endl;
	aus << "MMMMM  AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
	aurostd::PrintMessageStream(aus,XHOST.QUIET);
	exit(0);
      }
      if(AFLOW_PTHREADS::MAX_PTHREADS>1) {
	aus << "MMMMM  AFLOW: MULTI THREAD START: phread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(aus,XHOST.QUIET);
      }
      aus << "MMMMM  AFLOW: searching subdirectories [d1] " << " - " << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      while(!EMPTY) {
	vector<string> vaflowin;
	//   cerr << "aflags.AFLOW_PERFORM_DIRECTORY=" << aflags.AFLOW_PERFORM_DIRECTORY << endl;
	// cerr << "aflags.AFLOW_PERFORM_FILE=" << aflags.AFLOW_PERFORM_FILE << endl;

	if(aflags.AFLOW_PERFORM_FILE==FALSE) { // NO FILE SPECIFIED = standard  RUN DIRECTORY
	  aus << "MMMMM  AFLOW: aflags.AFLOW_PERFORM_FILE==FALSE" << " - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	
	  string FileNameSUBDIR=aurostd::TmpFileCreate("RUN");
	  stringstream strstream;
	  aus << "rm -f " <<  FileNameSUBDIR << endl;
	  aurostd::execute(aus);  // RESET  // RESET
	
	  bool isPRIORITY=FALSE;
	  // NEW, the sorting is done internally (speed and reliability)
	  for(int ifind=0;ifind<(int)vruns.size();ifind++) {
	    isPRIORITY=(aurostd::substring2bool(vruns.at(ifind),"PRIORITY") || aurostd::substring2bool(vruns.at(ifind),"priority"));
	    if((isPRIORITY && aurostd::uniform(1.0)<=PRIORITY_PROBABILITY) || !isPRIORITY) {
	      aus << "find " << vruns.at(ifind) << " " << XHOST.Find_Parameters;
	      if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
	      // if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
	      // if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
	      // if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)  aus << " -name \"" << _AFLOWIN_ << "\" ";   // IS THIS CORRECT ?? CHECK !!!
	      if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     aus << " -name \"INCAR\" ";
	      aus << " | " << PRIORITY_GREP_STRING << " >> ";
	      aus << FileNameSUBDIR << endl;
	    }
	  }
	  // perform the command
	  aurostd::execute(aus);  // RESET  // RESET
	
	  if(_VERBOSE_) aus << "MMMMM  AFLOW: searching subdirectories [d2] " << " - " << Message(aflags,"user,host,time") << endl;
	  if(_VERBOSE_) aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  aurostd::string2tokens(aurostd::file2string(FileNameSUBDIR),vaflowin,"\n");
	  for(i=0;i<(int) vaflowin.size();i++) aurostd::StringSubst(vaflowin.at(i),_AFLOWIN_,"");
	  for(i=0;i<(int) vaflowin.size();i++) aurostd::StringSubst(vaflowin.at(i),"INCAR","");
	  if(STOP_DEBUG) for(i=0;i<(int) vaflowin.size();i++) cout << vaflowin.at(i) << endl;
	  // RANDOMIZING priority
	  if(vaflowin.size()>1) { // only if I can poll
	    // if(aurostd::substring2bool(vaflowin.at(0),"PRIORITY") || aurostd::substring2bool(vaflowin.at(0),"priority")) aurostd::random_shuffle(vaflowin);
	  }
	  // loaded up
	  aus << "rm -f " << FileNameSUBDIR << endl;
	  aurostd::execute(aus);  // RESET  // RESET
	}
      
	// FILE SPECIFIED
	if(aflags.AFLOW_PERFORM_FILE) {
	  aus << "MMMMM  AFLOW: aflags.AFLOW_PERFORM_FILE==TRUE" << " - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  vaflowin.clear();
	  for(uint i=0;(i<vruns.size() && vaflowin.size()<maxcheck);i++) {
	    if(KBIN::Legitimate_aflowin(vruns.at(i),FALSE,aus)) vaflowin.push_back(vruns.at(i)); // TRUE puts too much verbosity
	    // vaflowin.push_back(vruns.at(i)); // just load them up... they were checked before //OLD MUST RECHECH THEM as things change on the fly
	  }
	}
	// NOW TIME OF SORTING/RANDOMIZING
	if(aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
	  aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  aurostd::sort(vaflowin);}
	if(aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
	  aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  aurostd::sort(vaflowin); // PN modification 
	  aurostd::rsort(vaflowin);}
	if(aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
	  aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] - " << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  aurostd::random_shuffle(vaflowin);}  // uses the std library but the seed is initialized in xrandom too
	aus << "MMMMM  " <<  aurostd::PaddedPOST("Legitimate VAFLOWIN = "+aurostd::utype2string(vaflowin.size()),40) << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	//      aus << "MMMMM  Legitimate VAFLOWIN = " << vaflowin.size() << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);

	//	for(uint i=0;i<vaflowin.size();i++) cerr << i << " " << vaflowin.at(i) << endl; 	exit(0);
	
	// clean AFLOWIN // SAFETY
	for(uint i=0;i<vaflowin.size();i++) aurostd::StringSubst(vaflowin.at(i),_AFLOWIN_,"");  
      
	if(MULTI_DEBUG) {aus << "MMMMM  SIZE vaflowin=" << vaflowin.size() << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
	if(MULTI_DEBUG) {for(uint i=0;i<vaflowin.size();i++) {aus << "MMMMM  vaflowin.at(i)=" << vaflowin.at(i) << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}}
	// cerr << vaflowin.at(10) << " " << vaflowin.at(20) << " " << vaflowin.at(30) << endl;
	// exit(0);
	// if(aurostd::substring2bool(vaflowin.at(0),"PRIORITY") || aurostd::substring2bool(vaflowin.at(0),"priority")) aurostd::random_shuffle(vaflowin);
      
	FOUND=FALSE;
	for(uint i=0;i<vaflowin.size()&& !FOUND;i++) {
	  aflags.Directory="NULL";
	  if(aurostd::DirectoryLocked(vaflowin.at(i),_AFLOWLOCK_)) {
	    if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  LOCKED ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
	    if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
	    FOUND=FALSE;
	  } else {
	    if(DirectorySkipped(vaflowin.at(i))) {
	      if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  SKIPPED ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
	      if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
	      FOUND=FALSE;
	    } else {
	      if(DirectoryAlreadyInDatabase(vaflowin.at(i),aflags.AFLOW_FORCE_RUN)) {
		if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
		if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ... use \"aflow --multi\" to force the calculation of this entry...  MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
		if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
		FOUND=FALSE;
	      } else {
		if(DirectoryUnwritable(vaflowin.at(i))) {
		  if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  UNWRITABLE ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
		  if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
		  FOUND=FALSE;
		} else {
		  if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  GOOD ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
		  if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
		  aflags.Directory=vaflowin.at(i);
		  FOUND=TRUE;
		} // DIRECTORY WRITABLE
	      } // DIRECTORY NOT IN THE DATABASE
	    } // DIRECTORY UN SKIPPED
	  } // DIRECTORY UN LOCKED
	} // DIRECTORY FOUND
	// exiting if STOP_DEBUG
	if(STOP_DEBUG) {cout << "aflow_kbin.cpp: STOP_DEBUG" << endl; exit(0);}
      
	// FOUND SOMETHING
	if(FOUND==FALSE) {
	  EMPTY=TRUE;
	  if(AFLOW_PTHREADS::FLAG) {
	    aus << "MMMMM  AFLOW: MULTI-THREADED: FLUSHING PTHREADS - " << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(aus,XHOST.QUIET);
	    for(ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
	      if(AFLOW_PTHREADS::vpthread_busy[ithread]) {
		aus << "MMMMM  AFLOW: MULTI-THREADED: Flushing   pthread=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << Message(aflags,"user,host,time") << endl;
		aurostd::PrintMessageStream(aus,XHOST.QUIET);
		pthread_join(AFLOW_PTHREADS::vpthread[ithread],NULL);
	      }
	  }
	}
	// again another check for LOCK, because NFS (network file system might be slow in concurrent seaches
	//     if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_) && FOUND) {cerr << "AFLOW EXCEPTION on concurrent LOCK: " << aflags.Directory << endl; FOUND=FALSE;}
	if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_) && FOUND) {
	  aus << "AFLOW EXCEPTION on concurrent LOCK: " << aflags.Directory << endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  FOUND=FALSE;
	}
      
	//  ---------------------------------------------------------------------------- START RUNNING
	if(FOUND) {
	  // again another check for LOCK, because NFS (network file system might be slow in concurrent seaches
	  EMPTY=FALSE;
	  //  ---------------------------------------------------------------------------- KIN_RUN_AFLOWIN
	  if(aflags.KBIN_RUN_AFLOWIN) {
	    //	  bool PHONONS;
	    //  -------------------------------------------------------------------------- KIN_RUN_AFLOWIN multithreaded
	    // 	  bool found;
	    // 	  for(uint ii=0;ii<qaflags.size()&&!found;ii++)
	    // 	    found=(qaflags[ii].Directory==aflags.Directory);       // look in all the list of operations
	    // 	  if(found==FALSE) {                                 // new operation, generate and save it
	    // 	    qaflags.push_back(aflags);
	    // 	  }
	    //	  cerr << qaflags.size() << endl;
	    if(AFLOW_PTHREADS::FLAG) {
	      // there is something to run in aflags.
	      // wait and put in ithread there is the number of the thread
	      free_thread=AFLOW_PTHREADS::Wait_Available_Free_Threads(ithread,_VERBOSE_);        // WAIT A WHILE !!
	      if(free_thread) {
		aus << "MMMMM  AFLOW: Found subdirectory to run " << aflags.Directory<< " - " << Message(aflags,"user,host,time") << endl;
		aus << "MMMMM  AFLOW: MULTI-THREADED: Starting    pthread_free=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << Message(aflags,"user,host,time") << endl;
		aurostd::PrintMessageStream(aus,XHOST.QUIET);
		aflags.AFLOW_PTHREADS_NUMBER=ithread;
		KBIN::RUN_Directory_PTHREADS(aflags);
		RUN_times++;
		if(XHOST.AFLOW_RUNXflag) {
		  aus << "MMMMM  AFLOW: RUNFX finished running " << RUN_times <<  "/" << XHOST.AFLOW_RUNXnumber << " " << aflags.Directory<< " - " << Message(aflags,"user,host,time") << endl;
		  aurostd::PrintMessageStream(aus,XHOST.QUIET);}
		if(XHOST.AFLOW_RUNXflag && RUN_times==XHOST.AFLOW_RUNXnumber) EMPTY=TRUE; // force to end if RUXN reached
	      }
	    }
	    //  -------------------------------------------------------------------------- KIN_RUN_AFLOWIN normal
	    if(!AFLOW_PTHREADS::FLAG) {
	      KBIN::RUN_Directory(aflags);
	      RUN_times++;
	      if(XHOST.AFLOW_RUNXflag) {
		aus << "MMMMM  AFLOW: RUNFX finished running " << RUN_times <<  "/" << XHOST.AFLOW_RUNXnumber << " " << aflags.Directory<< " - " << Message(aflags,"user,host,time") << endl;
		aurostd::PrintMessageStream(aus,XHOST.QUIET);}
	      if(XHOST.AFLOW_RUNXflag && RUN_times==XHOST.AFLOW_RUNXnumber) EMPTY=TRUE; // force to end if RUXN reached
	    }
	  }
	  //  ---------------------------------------------------------------------------- KBIN_GEN_VASP_FROM_AFLOWIN normal
	  //	if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)     KBIN::RUN_Directory(argv,aflags);  // IS THIS CORRECT ?? CHECK !!!
	  //  ---------------------------------------------------------------------------- KBIN_GEN_AFLOWIN_FROM_VASP normal
	  if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
	    KBIN::GenerateAflowinFromVASPDirectory(aflags);
	  }
	}
	if(XHOST.vflag_aflow.flag("LOOP") && EMPTY && XHOST.AFLOW_RUNXflag==FALSE) {
	  EMPTY=FALSE;
	  aus << "MMMMM  AFLOW: waiting for new subdirectories: " << (int) _KBIN_LOOP_SLEEP_/60 << "mins ";
	  aus << " - " << XHOST.hostname << " - " << aflow_get_time_string() << endl;//endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	  aurostd::Sleep(_KBIN_LOOP_SLEEP_);
	}
      }
      aus << "MMMMM  AFLOW: no more subdirectories to run " << " - " << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
    }

    // ------------------------------------------------------------------------------------------------------------------------------------

    return 1;
  }
} // namespace end of MAIN


// ***************************************************************************
// KBIN::MPI_Extract
// ***************************************************************************
// This function extracts from _AFLOWIN_ the parameters for MPI run
namespace KBIN {
  void MPI_Extract(string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
    ostringstream aus;
    bool Kmpi=TRUE;
    kflags.KBIN_MPI_NCPUS=0;
    aus << "00000  [AFLOW_MODE_MPI] found in " << _AFLOWIN_ << " " << Message(aflags,"user,host,time") << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // get (integer) kflags.KBIN_MPI_NCPUS
  

    if(aflags.AFLOW_GLOBAL_NCPUS<1) {
      if(Kmpi && !aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=",TRUE) && !aflags.AFLOW_MACHINE_LOCAL.flag()) {                      // DEFAULT NO CPU SPECIFIED
	kflags.KBIN_MPI_NCPUS=MPI_NCPUS_DEFAULT;
	aus << "00000  MESSAGE MPI: NCPUS=NNNN is missing, taking NCPUS=" << kflags.KBIN_MPI_NCPUS << "  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
	Kmpi=FALSE;
      }
      if(Kmpi && (aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=MAX",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=AUTO",TRUE) || aflags.AFLOW_MACHINE_LOCAL.flag())) { // DEFAULT NCPUS=MAX
	kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;

	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;	// corey
	//XHOST.PBS_NUM_PPN;        // with DUKE_BETA force NCPUS from QUEUE
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;        // with DUKE_BETA force NCPUS from QUEUE
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN; // corey
	//HOST.PBS_NUM_PPN;        // with DUKE_QRATS force NCPUS from QUEUE
 	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;   // with DUKE_MATERIALS force NCPUS
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;   // with DUKE_AFLOWLIB force NCPUS
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;   // with DUKE_HABANA force NCPUS
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_RANGER")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NNODES ;// with TERAGRID_RANGER force NCPUS from QUEUE
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_KRAKEN")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NNODES; // with TERAGRID_KRAKEN force NCPUS
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TRINITY_PARSONS")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;         // with TRINITY_PARSONS force NCPUS
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU"))  kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;        // with FULTON_MARYLOU force NCPUS
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE2")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;             // MACHINE2 has only NCPUS
	if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE1")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;            // MACHINE1 has only NCPUS
	if(kflags.KBIN_MPI_NCPUS<1) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores; // SAFE
      
	aus << "00000  MESSAGE MPI: found NCPUS=MAX  NCPUS="<<kflags.KBIN_MPI_NCPUS<<" " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
	Kmpi=FALSE;
      }
      if(Kmpi && aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=",TRUE)) {                     // DEFAULT NCPUS=XXX
	kflags.KBIN_MPI_NCPUS=aurostd::substring2utype<int>(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=",TRUE);
	if(kflags.KBIN_MPI_NCPUS>0) {
	  aus << "00000  MESSAGE MPI: found NCPUS="<<kflags.KBIN_MPI_NCPUS<<" " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
	}
	Kmpi=FALSE;
      }
    } else {
      kflags.KBIN_MPI_NCPUS=aflags.AFLOW_GLOBAL_NCPUS;
      aus << "00000  MESSAGE MPI: NCPUS is overriden, taking NCPUS=" << kflags.KBIN_MPI_NCPUS << "  " << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
    }
  
    if(kflags.KBIN_MPI_NCPUS<1) kflags.KBIN_MPI_NCPUS=1;                                              // DEFAULT NCPUS=troubles
  
    if(kflags.KBIN_MPI_NCPUS==1) {
      kflags.KBIN_MPI=FALSE;
      aus << "00000  MESSAGE MPI: found NCPUS=1 " << Message(aflags,"user,host,time") << endl;
      aus << "00000  MESSAGE MPI: going back to SERIAL execution " << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
    }
    if(kflags.KBIN_MPI) {
      // get (string) kflags.KBIN_MPI_START
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]START=",TRUE)) {
	kflags.KBIN_MPI_START=MPI_START_DEFAULT;
	aus << "00000  MESSAGE MPI: START string is missing, taking START=\"" << kflags.KBIN_MPI_START << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
	kflags.KBIN_MPI_START=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]START=",TRUE),'"');
	aus << "00000  MESSAGE MPI: found START=\"" << kflags.KBIN_MPI_START << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      // get (string) kflags.KBIN_MPI_STOP
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]STOP=",TRUE)) {
	kflags.KBIN_MPI_STOP=MPI_STOP_DEFAULT;
	aus << "00000  MESSAGE MPI: STOP string is missing, taking STOP=\"" << kflags.KBIN_MPI_STOP << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
	kflags.KBIN_MPI_STOP=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]STOP=",TRUE),'"');
	aus << "00000  MESSAGE MPI: found STOP=\"" << kflags.KBIN_MPI_STOP << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      // get (string) kflags.KBIN_MPI_COMMAND
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]COMMAND=",TRUE)) {
	kflags.KBIN_MPI_COMMAND=MPI_COMMAND_DEFAULT;
	aus << "00000  MESSAGE MPI: COMMAND string is missing, taking COMMAND=\"" << kflags.KBIN_MPI_COMMAND << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
	kflags.KBIN_MPI_COMMAND=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]COMMAND=",TRUE),'"');
	aus << "00000  MESSAGE MPI: found COMMAND=\"" << kflags.KBIN_MPI_COMMAND << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      kflags.KBIN_MPI_AUTOTUNE=aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]AUTOTUNE",TRUE);
      if(kflags.KBIN_MPI_AUTOTUNE) {
	aus << "00000  MESSAGE MPI: found AUTOTUNE option " << Message(aflags,"user,host,time") << endl;
	aus << "00000  MESSAGE MPI: input files WILL be auto-tuned for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs "  << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
	aus << "00000  MESSAGE MPI: AUTOTUNE option NOT found " << Message(aflags,"user,host,time") << endl;
	aus << "00000  MESSAGE MPI: input files MUST be appropriate for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs "  << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      // get (string) kflags.KBIN_MPI_BIN
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]BINARY=",TRUE)) {
	kflags.KBIN_MPI_BIN=DEFAULT_VASP_MPI_BIN;
	aus << "00000  MESSAGE MPI: BINARY string is missing, taking BIN=\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
	kflags.KBIN_MPI_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]BINARY=",TRUE),'"');
	aus << "00000  MESSAGE MPI: found BINARY=\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      aus << "00000  MESSAGE MPI: Overriding BINARY=\"" << kflags.KBIN_BIN << "\" to BINARY =\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      kflags.KBIN_BIN=kflags.KBIN_MPI_BIN;
  
      // get (string) kflags.KBIN_MPI_OPTIONS
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]OPTIONS=",TRUE)) {
	kflags.KBIN_MPI_OPTIONS=VASP_OPTIONS_MPI_DEFAULT;
	aus << "00000  MESSAGE MPI: OPTIONS string is missing, taking OPTIONS=\"" << kflags.KBIN_MPI_OPTIONS << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
	kflags.KBIN_MPI_OPTIONS=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]OPTIONS=",TRUE),'"');
	aus << "00000  MESSAGE MPI: found OPTIONS=\"" << kflags.KBIN_MPI_OPTIONS << "\"  " << Message(aflags,"user,host,time") << endl;
	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
    }
  }
} // namespace KBIN


// ***************************************************************************
// KBIN::StartStopCheck
// ***************************************************************************
namespace KBIN {
  void StartStopCheck(const string &AflowIn,string str1,string str2,bool &flag,bool &flagS) {
    flag =
      aurostd::substring2bool(AflowIn,str1) || aurostd::substring2bool(AflowIn,str2) ;
    flagS= (aurostd::substring2bool(AflowIn,str1+"START") && aurostd::substring2bool(AflowIn,str1+"STOP")) ||
      (aurostd::substring2bool(AflowIn,str2+"_START") && aurostd::substring2bool(AflowIn,str2+"_STOP"));
    if(flagS) flag=FALSE;
  }
}

namespace KBIN {
  void StartStopCheck(const string &AflowIn,string str1,bool &flag,bool &flagS) {
    flag = aurostd::substring2bool(AflowIn,str1);
    flagS= aurostd::substring2bool(AflowIn,str1+"START") && aurostd::substring2bool(AflowIn,str1+"STOP");
    if(flagS) flag=FALSE;
  }
}

// ***************************************************************************
// KBIN::RUN_Directory
// ***************************************************************************
namespace KBIN {
  void RUN_Directory(_aflags& aflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    ifstream FileSUBDIR;string FileNameSUBDIR;
    FileNameSUBDIR=aflags.Directory;
    FileSUBDIR.open(FileNameSUBDIR.c_str(),std::ios::in);
    FileSUBDIR.clear();FileSUBDIR.close();
    // string::size_type sub_size1,sub_size2;
    string AflowIn,AflowInMode,subS,subS1,subS2;
    string::iterator pos;
    bool Krun=TRUE;
    //  int i;
    _kflags kflags;
  
    if(aflags.Directory.at(0)!='/' && aflags.Directory.at(0)!='.' && aflags.Directory.at(0)!=' ') aflags.Directory="./"+aflags.Directory;

    if(!FileSUBDIR) {                                                                                           // ******* Directory is non existent
      aus << "EEEEE  DIRECTORY_NOT_FOUND = "  << Message(aflags,"user,host,time") << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
    } else {                                                                                                    // ******* Directory EXISTS
      // ***************************************************************************
      // Check LOCK again
      if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_) || DirectorySkipped(aflags.Directory) || DirectoryAlreadyInDatabase(aflags.Directory,aflags.AFLOW_FORCE_RUN) || DirectoryUnwritable(aflags.Directory)) {
	// ******* Directory is locked/skipped/unwritable
	// LOCK/SKIP/UNWRITABLE exist, then RUN already RUN
	if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_)) {
	  aus << "LLLLL  LOCKED ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,"user,host,time") << endl;
	  aus << "LLLLL  LOCKED ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	}
	if(DirectorySkipped(aflags.Directory)) {
	  aus << "LLLLL  SKIPPED ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,"user,host,time") << endl;
	  aus << "LLLLL  SKIPPED ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	}
	if(DirectoryAlreadyInDatabase(aflags.Directory,aflags.AFLOW_FORCE_RUN)) {
	  aus << "LLLLL  ALREADY_IN_DATABASE ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,"user,host,time") << endl;
	  aus << "LLLLL  ALREADY_IN_DATABASE ... use \"aflow --multi\" to force the calculation of this entry "  << Message(aflags,"user,host,time") << endl;
	  aus << "LLLLL  ALREADY_IN_DATABASE ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	}
	if(DirectoryUnwritable(aflags.Directory)) {
	  aus << "LLLLL  UNWRITABLE ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,"user,host,time") << endl;
	  aus << "LLLLL  UNWRITABLE ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	}
	//     exit(1);
      } else {                                                                                                  // ******* Directory is fine
	// make a dumb lock as soon as possible -------------------------------------
	aus.clear();aus.str(std::string());
	aus << "echo \"NNNNN  KBIN LOCK ASAP for NFS concurrent jobs (aflow" << string(AFLOW_VERSION) << ")\" >> " << aflags.Directory+"/"+_AFLOWLOCK_ << endl;
	// aus << "/home/auro/bin/aflow -machine >> " << aflags.Directory+"/"+_AFLOWLOCK_ << endl;
	// aus << XHOST.command("sensors") << " >> " << aflags.Directory+"/"+_AFLOWLOCK_ << endl;
	aurostd::execute(aus);
	// now change its permission
	aurostd::ChmodFile("664",string(aflags.Directory+"/"+_AFLOWLOCK_));
	// now the lock should be done ----------------------------------------------
	ifstream FileAFLOWIN;string FileNameAFLOWIN;
	FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
	FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::in);
	if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
	  aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = "  << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(aus,XHOST.QUIET);
	} else {                                                                                                // ******* _AFLOWIN_ exists RUN    
	  // ***************************************************************************
	  // RESET LOCK
	  ofstream FileLOCK;
	  string FileNameLOCK=aflags.Directory+"/"+_AFLOWLOCK_;
	  //	FileLOCK.open(FileNameLOCK.c_str(),std::ios::out);
	  FileLOCK.open(FileNameLOCK.c_str(),std::ios::app);
	  // ***************************************************************************
	  // WRITE LOCK
	  if(0) {
	    aus <<    "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow - " << Message(aflags,"user,host,time") << endl;
	    aus <<    "MMMMM  (C) "<<XHOST.Copyright_Years<<", Stefano Curtarolo - Duke University   - " << Message(aflags,"user,host,time") << endl;
	    aus <<    "MMMMM  High-Throughput ab-initio Computing - " << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    // ***************************************************************************
	    // WRITE AFLOW VERSION
	    // aus << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow " << Message(aflags,"user,host,time") << endl;
	    // aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  }
	  // ***************************************************************************
	  // START DIRECTORY
	  aus      << "XXXXX  KBIN DIRECTORY BEGIN (aflow" << string(AFLOW_VERSION) << ")  "  << Message(aflags,"user,host,time") << endl;
	  //	aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  aus      << "XXXXX  KBIN XHOST.CPU_Model : "<<  XHOST.CPU_Model << "" << endl;// << Message(aflags,"user,host,time") << endl;
	  aus      << "XXXXX  KBIN XHOST.CPU_Cores : "<<  XHOST.CPU_Cores << "" << endl;// << Message(aflags,"user,host,time") << endl;
	  aus      << "XXXXX  KBIN XHOST.CPU_MHz   : "<<  XHOST.CPU_MHz << "" << endl;// << Message(aflags,"user,host,time") << endl;
	  aus      << "XXXXX  KBIN XHOST.RAM_GB    : "<<  XHOST.RAM_GB << "" << endl;// << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  // ***************************************************************************
	  // FLUSH & REOPEN to avoid double writing
	  FileLOCK.flush();FileLOCK.close();FileLOCK.open(FileNameLOCK.c_str(),std::ios::app);
	  // ***************************************************************************
	  // NOW Digest AFLOWIN
	  FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
	  AflowIn.clear();char c; while (FileAFLOWIN.get(c)) AflowIn+=c;               // READ _AFLOWIN_ and put into AflowIn
	  FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
	  AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
	  // ***************************************************************************
	  // FIND MPI	
	  kflags.KBIN_MPI= aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI]");  // search for MPI string
	  // ***************************************************************************
	  // FIND HOST
	  // duke_beta	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]BETA") || 
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_BETA"))   // check DUKE_BETA
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // duke_beta_openmpi	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]BETA_OPENMPI") ||    // check DUKE_BETA_OPENMPI
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_BETA_OPENMPI"))   // check DUKE_BETA_OPENMPI
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // duke_qrats	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]QRATS") || 
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_QRATS"))   // check DUKE_QRATS
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // duke_quser	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QUSER_OPENMPI") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]QUSER") || 
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_QUSER"))   // check DUKE_QUSER
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QUSER_OPENMPI")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // mpcdf_eos	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS_MPIIFORT") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]EOS") || 
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MPCDF_EOS"))   // check MPCDF_EOS
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS_MPIIFORT")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // mpcdf_draco	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO_MPIIFORT") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DRACO") || 
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MPCDF_DRACO"))   // check MPCDF_DRACO
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO_MPIIFORT")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // mpcdf_hydra	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA_MPIIFORT") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]HYDRA") || 
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MPCDF_HYDRA"))   // check MPCDF_HYDRA
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA_MPIIFORT")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // duke_materials	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_MATERIALS") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MATERIALS") ||    // check DUKE_MATERIALS
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_MATERIALS"))   // check DUKE_MATERIALS
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // duke_aflowlib	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_AFLOWLIB") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]AFLOWLIB") ||    // check DUKE_AFLOWLIB
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_AFLOWLIB"))   // check DUKE_AFLOWLIB
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // duke_habana	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_HABANA") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]HABANA") ||    // check DUKE_HABANA
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_HABANA"))   // check DUKE_HABANA
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // teragrid_ranger	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::TERAGRID_RANGER") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]RANGER") ||    // check TERAGRID_RANGER
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]TERAGRID_RANGER"))   // check TERAGRID_RANGER
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_RANGER")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // teragrid_kraken	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::TERAGRID_KRAKEN") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]KRAKEN") ||    // check TERAGRID_KRAKEN
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]TERAGRID_KRAKEN"))   // check TERAGRID_KRAKEN
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TERAGRID_KRAKEN")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // trinity_parsons	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::TRINITY_PARSONS") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]PARSONS") ||    // check TRINITY_PARSONS
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]TRINITY_PARSONS"))   // check TRINITY_PARSONS
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::TRINITY_PARSONS")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // fulton_marylou	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::FULTON_MARYLOU") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MARYLOU") ||    // check FULTON_MARYLOU
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]FULTON_MARYLOU"))   // check FULTON_MARYLOU
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // ohad	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE2") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE2") ||  // check MACHINE2
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE2"))   // check MACHINE2
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
	  if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE2")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }
	  // host1	
	  if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE1") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE1") ||  // check MACHINE1
	     aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE1"))   // check MACHINE1
	    aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
          if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE1")) {
	    aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL << " " << Message(aflags,"user,host,time","aflow_kbin.cpp") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
	  }

	  // ***************************************************************************
	  // OTHER CHECKS FOR MPI
	  // machines are done withing the VASP/ALIEN stuff, if necessary
	  if(aflags.AFLOW_FORCE_MPI) kflags.KBIN_MPI=TRUE;      // forcing
	  if(aflags.AFLOW_FORCE_SERIAL) kflags.KBIN_MPI=FALSE;  // forcing

	  kflags.KBIN_QSUB= aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]") && !aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE");  // search for QSUB string
	  kflags.KBIN_QSUB_MODE1=aflags.AFLOW_MODE_QSUB_MODE1 || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE1"); // search for QSUB string mode1
	  kflags.KBIN_QSUB_MODE2=aflags.AFLOW_MODE_QSUB_MODE2 || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE2"); // search for QSUB string mode2
	  kflags.KBIN_QSUB_MODE3=aflags.AFLOW_MODE_QSUB_MODE3 || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE3"); // search for QSUB string mode3
	  kflags.AFLOW_MODE_ALIEN=                                               // check ALIEN
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE=ALIEN]") ||             // check ALIEN
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ALIEN]") ||             // check ALIEN
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE]ALIEN");                // check ALIEN
	  kflags.AFLOW_MODE_MATLAB=                                              // check MATLAB
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE=MATLAB]") ||            // check MATLAB
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MATLAB]") ||            // check MATLAB
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE]MATLAB");               // check MATLAB
	  kflags.AFLOW_MODE_VASP=                                                // check VASP
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE=VASP]") ||              // check VASP
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE_VASP]") ||              // check VASP
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE]VASP");                 // check VASP
	  kflags.AFLOW_MODE_AIMS=                                                // check AIMS
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE=AIMS]") ||              // check AIMS
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE_AIMS]") ||              // check AIMS
	    aurostd::substring2bool(AflowIn,"[AFLOW_MODE]AIMS");                 // check AIMS
	  kflags.KBIN_SYMMETRY_CALCULATION  = aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_SYMMETRY]CALC",TRUE);
          //DX - START
	  kflags.KBIN_SYMMETRY_NO_SCAN  = aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]NO_SCAN",TRUE);
          //cerr << kflags.KBIN_SYMMETRY_EPS << endl;
          if(aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]SYM_EPS=",TRUE)){
            kflags.KBIN_SYMMETRY_EPS      = aurostd::substring2utype<double>(AflowIn,"[AFLOW_SYMMETRY]SYM_EPS=",TRUE);
          }
          //DX - END
	  // ---------------------------------------------------------
	  // parameters for AAPL - CO - 170601
    // to make backwards compatible, we need to not only look for substring, but need to see if "KAPPA=y"
    // start with AAPL first, then QHA, then APL, they are mutually exclusive
    aurostd::xoption KBIN_PHONONS_CALCULATION_AAPL;
    KBIN_PHONONS_CALCULATION_AAPL.option=false;
    KBIN_PHONONS_CALCULATION_AAPL.options2entry(AflowIn, string("[AFLOW_AAPL]KAPPA=|[AFLOW_PHONONS]KAPPA="), KBIN_PHONONS_CALCULATION_AAPL.option, KBIN_PHONONS_CALCULATION_AAPL.scheme); //CO 170601
	  KBIN_PHONONS_CALCULATION_AAPL.option |= aurostd::substring2bool(AflowIn,"[AFLOW_AAPL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AAPL]CALC",TRUE);  //legacy
	  kflags.KBIN_PHONONS_CALCULATION_AAPL  = KBIN_PHONONS_CALCULATION_AAPL.option;
	  // ---------------------------------------------------------
	  // parameters for QHA - CO - 170601
    // to make backwards compatible, we need to not only look for substring, but need to see if "GRUNEISEN=y"
    if(!kflags.KBIN_PHONONS_CALCULATION_AAPL){  //mutually exclusive
      aurostd::xoption KBIN_PHONONS_CALCULATION_QHA;
      KBIN_PHONONS_CALCULATION_QHA.option=false;
      KBIN_PHONONS_CALCULATION_QHA.options2entry(AflowIn, string("[AFLOW_QHA]GRUNEISEN=|[AFLOW_PHONONS]GRUNEISEN="), KBIN_PHONONS_CALCULATION_QHA.option, KBIN_PHONONS_CALCULATION_QHA.scheme); //CO 170601
      KBIN_PHONONS_CALCULATION_QHA.option |= aurostd::substring2bool(AflowIn,"[AFLOW_QHA]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_QHA]CALC",TRUE); //legacy
      kflags.KBIN_PHONONS_CALCULATION_QHA  = KBIN_PHONONS_CALCULATION_QHA.option;
    }
	  // ---------------------------------------------------------
	  // parameters for APL
    if(!(kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_QHA)){ //mutually exclusive
	    kflags.KBIN_PHONONS_CALCULATION_APL  = aurostd::substring2bool(AflowIn,"[AFLOW_APL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_PHONONS]CALC",TRUE);
    }
	  // ---------------------------------------------------------
	  // parameters for AGL (Debye Model)
	  kflags.KBIN_PHONONS_CALCULATION_AGL  = aurostd::substring2bool(AflowIn,"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AGL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_GIBBS]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_GIBBS]CALC",TRUE);
	  // ---------------------------------------------------------
	  // parameters for AEL (Elastic constants)
	  kflags.KBIN_PHONONS_CALCULATION_AEL  = aurostd::substring2bool(AflowIn,"[AFLOW_AEL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AEL]CALC",TRUE);
	  // ---------------------------------------------------------
	  // parameters for NEIGHBOURS
	  kflags.KBIN_NEIGHBOURS_CALCULATION  = aurostd::substring2bool(AflowIn,"[AFLOW_NEIGHBOURS]CALC",TRUE);
	  // ---------------------------------------------------------
	  // parameters for POCC CALCULATIONS, KESONG YANG
	  kflags.KBIN_POCC=FALSE;
	  kflags.KBIN_POCC_CALCULATION  = aurostd::substring2bool(AflowIn,"[AFLOW_POCC]CALC",TRUE);
	  if(kflags.KBIN_POCC_CALCULATION) {
	    aus << "00000  MESSAGE POCC_CALCULATION "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  }
	  // ---------------------------------------------------------
	  // parameters for FROZSL
	  kflags.KBIN_FROZSL=FALSE;
	  kflags.KBIN_PHONONS_CALCULATION_FROZSL  = aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]CALC",TRUE);
	  if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
	    aus << "00000  MESSAGE FROZSL_CALCULATION "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  }
	  kflags.KBIN_FROZSL_DOWNLOAD     = (aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]DOWN",TRUE) ||
					     aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]DOWNLOAD",TRUE));
	  if(kflags.KBIN_FROZSL_DOWNLOAD) {
	    aus << "00000  MESSAGE FROZSL_DOWNLOAD "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  }
	  if(kflags.KBIN_FROZSL_FILE) {
	    aus << "00000  MESSAGE FROZSL_FILE "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  }
	  kflags.KBIN_FROZSL_FILE  = aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]FILE",TRUE); // load of file somewhere else
	  if(kflags.KBIN_PHONONS_CALCULATION_FROZSL || kflags.KBIN_FROZSL_DOWNLOAD|| kflags.KBIN_FROZSL_FILE) kflags.KBIN_FROZSL=TRUE;
	  // ---------------------------------------------------------
	  // the rest of symmetry stuff is seeked inside ivasp or
	  if(kflags.AFLOW_MODE_ALIEN) {
	    kflags.AFLOW_MODE_MATLAB=FALSE;                  // fix PRIORITY
	    kflags.AFLOW_MODE_VASP=FALSE;                    // fix PRIORITY
	    kflags.KBIN_MPI=FALSE;                           // fix PRIORITY
	  }
	  if(kflags.AFLOW_MODE_MATLAB) {
	    kflags.AFLOW_MODE_VASP=FALSE;                    // fix PRIORITY
	    kflags.KBIN_MPI=FALSE;
	  }
	  if(LDEBUG) cout << "DEBUG kflags.AFLOW_MODE_ALIEN=" << kflags.AFLOW_MODE_ALIEN << endl;
	  if(LDEBUG) cout << "DEBUG kflags.AFLOW_MODE_MATLAB=" << kflags.AFLOW_MODE_MATLAB << endl;
	  if(LDEBUG) cout << "DEBUG kflags.AFLOW_MODE_VASP=" << kflags.AFLOW_MODE_VASP << endl;
	  // ***************************************************************************
	  // ZIP COMPRESS
	  // ***************************************************************************
	  kflags.KZIP_COMPRESS=TRUE;
	  aurostd::StringstreamClean(aus);
	  if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=none]") ||
	     aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=NONE]") ||
	     !aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP")) {
	    kflags.KZIP_COMPRESS=FALSE;
	    for(int i=0;i<1;i++) {
	      aus << "WWWWW  Warning no compression of output files... " << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintWarningStream(FileLOCK,aus,XHOST.QUIET);
	    }
	  } else {
	    if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP")) { // "[AFLOW_MODE_ZIP=" not found
	      kflags.KZIP_BIN=DEFAULT_KZIP_BIN;  // take default
	      aus << "00000  MESSAGE Taking DEFAULT KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    }
	    if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP]")) { // "[AFLOW_MODE_ZIP]" not found
	      kflags.KZIP_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_ZIP]");
	      aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    }
	    if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=")) { // "[AFLOW_MODE_ZIP=" found
	      kflags.KZIP_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_ZIP="),']');
	      aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    }
	  }
	  // ************************************************************************************************************************************
	  // Get the KZIP_BIN name - moved inside EACH mode
	  // ************************************************************************************************************************************
	  // LOAD PREFIX POSTFIX
	  KBIN::StartStopCheck(AflowIn,"[AFLOW_MODE_PRESCRIPT]",kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT,kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP);
	  KBIN::StartStopCheck(AflowIn,"[AFLOW_MODE_POSTSCRIPT]",kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT,kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP);
	  if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT) {  // [AFLOW_MODE_PRESCRIPT] construction
	    aus << "00000  MESSAGE Generating " << _AFLOW_PRESCRIPT_COMMAND_ << " file from " << _AFLOWIN_ << " " << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,kflags.AFLOW_MODE_PRESCRIPT,"[AFLOW_MODE_PRESCRIPT]");
	  }
	  if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_PRESCRIPT] construction
	    aus << "00000  MESSAGE Generating " << _AFLOW_PRESCRIPT_COMMAND_ << " file from START/STOP " << _AFLOWIN_ << " " << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,kflags.AFLOW_MODE_PRESCRIPT,"[AFLOW_MODE_PRESCRIPT]START","[AFLOW_MODE_PRESCRIPT]STOP");
	  }
	  if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT) {  // [AFLOW_MODE_POSTSCRIPT] construction
	    aus << "00000  MESSAGE Generating " << _AFLOW_POSTSCRIPT_COMMAND_ << " file from " << _AFLOWIN_ << " " << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,kflags.AFLOW_MODE_POSTSCRIPT,"[AFLOW_MODE_POSTSCRIPT]");
	  }
	  if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_POSTSCRIPT] construction
	    aus << "00000  MESSAGE Generating " << _AFLOW_POSTSCRIPT_COMMAND_ << " file from START/STOP " << _AFLOWIN_ << " " << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,kflags.AFLOW_MODE_POSTSCRIPT,"[AFLOW_MODE_POSTSCRIPT]START","[AFLOW_MODE_POSTSCRIPT]STOP");
	  }
	  if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_PRESCRIPT] construction
	    aurostd::stringstream2file(kflags.AFLOW_MODE_PRESCRIPT,string(aflags.Directory+"/"+_AFLOW_PRESCRIPT_COMMAND_));
	    aurostd::ChmodFile("755",string(aflags.Directory+"/"+_AFLOW_PRESCRIPT_COMMAND_));
	  }
	  if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_POSTSCRIPT] construction
	    aurostd::stringstream2file(kflags.AFLOW_MODE_POSTSCRIPT,string(aflags.Directory+"/"+_AFLOW_POSTSCRIPT_COMMAND_));
	    aurostd::ChmodFile("755",string(aflags.Directory+"/"+_AFLOW_POSTSCRIPT_COMMAND_));
	  }
	  // ************************************************************************************************************************************
	  // ALIEN MODE
	  if(kflags.AFLOW_MODE_ALIEN) {
	    AflowInMode="[AFLOW_MODE=ALIEN]";
	    aus      << "00000  MESSAGE [AFLOW_MODE=ALIEN] found in " << _AFLOWIN_ << " "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    // ***************************************************************************
	    // Getting KBIN_BIN
	    kflags.KBIN_BIN = DEFAULT_KBIN_ALIEN_BIN;  // take default
	    aurostd::StringstreamClean(aus);
	    if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
	      kflags.KBIN_BIN=DEFAULT_KBIN_ALIEN_BIN;  // take default
	      aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    }
	    if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
	      kflags.KBIN_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY]");
	      aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    }
	    if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
	      kflags.KBIN_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
	      aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    }
	    // ***************************************************************************
	    // ALIEN MODE  // must contain EMAIL perform
	    if(Krun) {
	      kflags.AFLOW_MODE_EMAIL            =
		aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
		aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL") ;
	      Krun=(Krun && ALIEN::Run_Directory(FileLOCK,aflags,kflags));
	    }
	    // ***************************************************************************
	    // COMPRESS
	    if(Krun && kflags.KZIP_COMPRESS) KBIN::CompressDirectory(aflags,kflags);
	  }
	  // ************************************************************************************************************************************
	  // MATLAB MODE
	  if(kflags.AFLOW_MODE_MATLAB) {
	    AflowInMode="[AFLOW_MODE=MATLAB]";
	    aus      << "00000  MESSAGE [AFLOW_MODE=MATLAB] found in " << _AFLOWIN_ << " "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    // ***************************************************************************
	    // MATLAB MODE  // must contain EMAIL perform
	    if(Krun) {
	      kflags.AFLOW_MODE_EMAIL            =
		aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
		aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL") ;
	      aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
	      Krun=(Krun && KBIN_MATLAB_Directory(FileLOCK,aflags,kflags));
	    }
	    // ***************************************************************************
	    // COMPRESS
	    if(Krun && kflags.KZIP_COMPRESS) KBIN::CompressDirectory(aflags,kflags);
	    Krun=FALSE;
	  }
	  // ************************************************************************************************************************************
            // AIMS MODE
            if(kflags.AFLOW_MODE_AIMS) {
              AflowInMode="[AFLOW_MODE=AIMS]";
              aus      << "00000  MESSAGE [AFLOW_MODE=AIMS] found in " << _AFLOWIN_ << " "  << Message(aflags,"user,host,time") << endl;
              aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
              aurostd::StringstreamClean(aus);
              if(1){  //no support yet
                // ***************************************************************************
                // Getting KBIN_BIN
                kflags.KBIN_BIN = DEFAULT_AIMS_BIN;  // take default  dont touch MPI as it has already been dealt by  KBIN::MPI_Extract
                if(kflags.KBIN_MPI==FALSE) { // only if no MPI is specified
                  if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
                    kflags.KBIN_BIN=DEFAULT_AIMS_BIN;  // take default
                    aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                  }
                  if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
                    kflags.KBIN_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY]");
                    aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                  }
                  if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
                    kflags.KBIN_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
                    aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
                    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                  }
                } else {
                  kflags.KBIN_BIN=kflags.KBIN_MPI_BIN;
                  aus << "00000  MESSAGE Taking KBIN_BIN=KBIN_MPI_BIN=\"" << kflags.KBIN_MPI_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
                  aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
                }
              }
              // ***************************************************************************
              // AIMS MODE  // must contain EMAIL perform
              if(Krun) {
                kflags.AFLOW_MODE_EMAIL            =
                  aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
                  aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL");
                Krun=(Krun && KBIN::AIMS_Directory(FileLOCK,aflags,kflags));
              }
              // ***************************************************************************
              // COMPRESS
              if(Krun && kflags.KZIP_COMPRESS) KBIN::CompressDirectory(aflags,kflags);
            }
            // ************************************************************************************************************************************
	  // MPI SWTICHES
	  if(kflags.KBIN_MPI) KBIN::MPI_Extract(AflowIn,FileLOCK,aflags,kflags);
	  // ************************************************************************************************************************************
	  // ************************************************************************************************************************************
	  // VASP MODE
	  if(kflags.AFLOW_MODE_VASP) {
	    AflowInMode="[AFLOW_MODE=VASP]";
	    aus      << "00000  MESSAGE [AFLOW_MODE=VASP] found in " << _AFLOWIN_ << " "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    // ***************************************************************************
	    // Getting KBIN_BIN
	    kflags.KBIN_BIN = DEFAULT_VASP_BIN;  // take default  dont touch MPI as it has already been dealt by  KBIN::MPI_Extract
	    aurostd::StringstreamClean(aus);
	    // old Get BIN
	    // 	  if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" not found
	    // 	    aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(aflags,"user,host,time") << endl;
	    // 	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    // 	    //   cerr << "take KBIN=" << kflags.KBIN_BIN << endl;
	    // 	  } else {
	    // 	    kflags.KBIN_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
	    // 	    aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(aflags,"user,host,time") << endl;
	    // 	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    // 	  }
	    if(kflags.KBIN_MPI==FALSE) { // only if no MPI is specified
	      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
		kflags.KBIN_BIN=DEFAULT_VASP_BIN;  // take default
		aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
		aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	      }
	      if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
		kflags.KBIN_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY]");
		aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
		aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	      }
	      if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
		kflags.KBIN_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
		aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
		aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	      }
	    } else {
	      kflags.KBIN_BIN=kflags.KBIN_MPI_BIN;
	      aus << "00000  MESSAGE Taking KBIN_BIN=KBIN_MPI_BIN=\"" << kflags.KBIN_MPI_BIN << "\" "  << Message(aflags,"user,host,time") << endl;
	      aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    }
	    // ***************************************************************************
	    // VASP MODE  // must contain EMAIL perform
	    if(Krun) {
	      kflags.AFLOW_MODE_EMAIL            =
		aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
		aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL");
	      Krun=(Krun && KBIN::VASP_Directory(FileLOCK,aflags,kflags));
	    }
	    // ***************************************************************************
	    // COMPRESS
	    if(Krun && kflags.KZIP_COMPRESS) KBIN::CompressDirectory(aflags,kflags);
	  }
	  // ************************************************************************************************************************************
	  // MATLAB MODE
	  if(kflags.KBIN_PHONONS_CALCULATION_FROZSL && !kflags.AFLOW_MODE_VASP) {
	    AflowInMode="[AFLOW_FROZSL]CALC";
	    aus      << "00000  MESSAGE [AFLOW_FROZSL]CALC found in " << _AFLOWIN_ << " "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	    // PRESCRIPT
	    if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)
	      KBIN::RUN_DirectoryScript(aflags,_AFLOW_PRESCRIPT_COMMAND_,_AFLOW_PRESCRIPT_FILE_);
	    // POSTSCRIPT
	    if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
	      KBIN::RUN_DirectoryScript(aflags,_AFLOW_POSTSCRIPT_COMMAND_,_AFLOW_POSTSCRIPT_FILE_);
	  }
	  // ************************************************************************************************************************************
	  // NO MODE MODE
            if(!kflags.AFLOW_MODE_VASP && !kflags.AFLOW_MODE_AIMS && !kflags.AFLOW_MODE_MATLAB && !kflags.AFLOW_MODE_ALIEN && !kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
	    aus << "EEEEE  [AFLOW_MODE=????] invalid found in     "  << Message(aflags,"user,host,time") << endl;
	    aus << "EEEEE  [AFLOW_MODE=ALIEN]        is supported "  << Message(aflags,"user,host,time") << endl;
	    aus << "EEEEE  [AFLOW_MODE=MATLAB]       is supported "  << Message(aflags,"user,host,time") << endl;
	    aus << "EEEEE  [AFLOW_MODE=VASP]         is supported "  << Message(aflags,"user,host,time") << endl;
	    aus << "EEEEE  [AFLOW_FROZSL]CALC        is supported "  << Message(aflags,"user,host,time") << endl;
	    aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
	  }
	  // ***************************************************************************
	  // FINALIZE AFLOWIN
	  // DONE, turn OFF the flag
	  kflags.AFLOW_MODE_VASP=FALSE;    
	  // ***************************************************************************
	  // FINALIZE LOCK
	  aus      << "XXXXX  KBIN DIRECTORY END (aflow" << string(AFLOW_VERSION) << ")  "  << Message(aflags,"user,host,time") << endl;
	  aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
	  // ***************************************************************************
	  // FINALIZE LOCK
	  FileLOCK.clear();FileLOCK.close();
	  // ***************************************************************************
    // NFS cache-cleaning HACK, opendir() and closedir() to invalidate cache
    // https://stackoverflow.com/questions/8311710/nfs-cache-cleaning-command
    // CO 171106
    vector<string> vfiles_dummyls;
    aurostd::DirectoryLS(aflags.Directory,vfiles_dummyls);
	  // ***************************************************************************
	  // WRITE END
	  aurostd::string2file(string(Message(aflags,"user,host,time")+"\n"),string(aflags.Directory+"/"+_AFLOW_END_FILE_));
	  // ***************************************************************************
	  // MAKE READEABLE
	  aurostd::ChmodFile("664",string(aflags.Directory+"/"+_AFLOWLOCK_));
	  aus.clear();aus.str(std::string());
	  aus << "cp " << aflags.Directory << "/" << _AFLOWLOCK_ << " " << aflags.Directory << "/" << _AFLOWLOCK_ << "." << string(AFLOW_VERSION) << endl;  
	  aurostd::execute(aus);  
	  aurostd::ChmodFile("664",string(aflags.Directory+"/*"));
	  aurostd::ChmodFile("777",string(aflags.Directory+"/*"));
	}
	FileAFLOWIN.clear();FileAFLOWIN.close();
      }
    }
  };
} // namespace

// *******************************************************************************************
namespace KBIN {
  void AFLOW_RUN_Directory(const _aflags& aflags) {
    aurostd::execute(XHOST.command("aflow")+" --run=1 --DIRECTORY="+aflags.Directory);  // run it OUTSIDE
    // this is cool as if the particula program explodes, there is still aflow running
  }
} // namespace

// *******************************************************************************************
namespace KBIN {
  void RUN_DirectoryScript(const _aflags& aflags,const string& script,const string& output) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << "cd " << aflags.Directory << endl;
    aus << DEFAULT_CHMOD_BIN << " 755 " << script << endl;
    aus << "./" << script << " >> " << output << endl;
    aurostd::execute(aus);
  }
} // namespace

// *******************************************************************************************
namespace KBIN {
  void CompressDirectory(const _aflags& aflags,const _kflags& kflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    //DX and CO - START
    //aus << "cd " << aflags.Directory << " && ";
    //aus << "ls | grep -v .bz2 | ";                                //CO, skip anything with bzip extension
    //aus << "grep -v LOCK | grep -v " << _AFLOWLOCK_ << " | ";     //CO never zip LOCK or .LOCK (agl.LOCK) newly defined LOCK
    //aus << "grep -v SKIP | ";
    //aus << "grep -v " << KBIN_SUBDIRECTORIES << " | ";
    //aus << "grep -v aflow.in | grep -v " << _AFLOWIN_ << " ";;    //CO, never zip aflow.in or _aflow.in (agl_aflow.in) or newly defined aflow.in
    vector<string> _vfiles,vfiles;
    string compressed_variant;
    aurostd::DirectoryLS(aflags.Directory,_vfiles);
    string file_path;
    for(uint i=0;i<_vfiles.size();i++){
      file_path=aflags.Directory + "/" + _vfiles[i];
      if(aurostd::IsCompressed(_vfiles[i])){continue;}  //doesn't need full path, just substring2bool for zip variants, e.g., .bz2
      if(_vfiles[i].size() && _vfiles[i][0]=='.'){continue;}  //do not try to compress any hidden files, .nfs stuff is particularly problematic to zip
      if(aurostd::substring2bool(_vfiles[i],KBIN_SUBDIRECTORIES)){continue;}
      if(aurostd::substring2bool(_vfiles[i],"LOCK")){continue;}
      if(aurostd::substring2bool(_vfiles[i],_AFLOWLOCK_)){continue;}
      if(aurostd::substring2bool(_vfiles[i],"SKIP")){continue;}
      if(aurostd::substring2bool(_vfiles[i],"aflow.in")){continue;}
      if(aurostd::substring2bool(_vfiles[i],_AFLOW_END_FILE_) || aurostd::substring2bool(_vfiles[i],"aflow.end.out")){continue;}  //CO 170613, file is special because it gets written after compression
      if(aurostd::substring2bool(_vfiles[i],_AFLOWIN_)){continue;}
      if(aurostd::CompressedFileExist(file_path,compressed_variant)){ //need full path here, also, notice the placement here, actual compressed variant would have been skipped, this is for the uncompressed variant
        //both compressed and uncompressed variants exists
        //assume compressed is from a previous run, hence obsolete
        //delete it before compressing uncompressed variant
        aurostd::RemoveFile(compressed_variant);
      }
      vfiles.push_back(_vfiles[i]);
    }
    //aurostd::string2vectorstring(aurostd::execute2string(aus),vfiles);
    //cerr << vfiles.size() << endl;
    //aurostd::StringstreamClean(aus);
    //cerr << "COREY " << aurostd::joinWDelimiter(vfiles," ") << endl;
    //exit(0);
    if(vfiles.size()){
      ostringstream aus;
      //aurostd::StringstreamClean(aus);
      aus << "cd " << aflags.Directory << " && " << endl;
      for(uint i=0;i<vfiles.size();i++){ //better than doing it all in one shot
        aus << kflags.KZIP_BIN << " " << vfiles[i] << "; " << endl;  //semi-colon is important, keeps going if it stalls on one
      }
      //aus << kflags.KZIP_BIN << " " << aurostd::joinWDelimiter(vfiles," ") << endl; //AVOID, because if one fails, the whole command stops
      aurostd::execute(aus);
      //cerr << aus.str() << endl;
    }
    //aus << kflags.KZIP_BIN << " `find . " << XHOST.Find_Parameters << " -name \"*\" | grep -v LOCK | grep -v SKIP | grep \"./\"` " << endl;
    //aus << kflags.KZIP_BIN << " `ls | grep -v LOCK | grep -v SKIP | grep -v " << KBIN_SUBDIRECTORIES << "| grep -v " << _AFLOWIN_ << " | grep -v apl.xml ` " << endl;
    //apl.xml can now be zipped
    //set up is slightly redundant (LOCK vs. agl.LOCK), but very safe
    //aus << kflags.KZIP_BIN << " `ls | grep -v .bz2 | ";                 //CO, skip anything with bzip extension
    //aus << "grep -v LOCK | grep -v " << _AFLOWLOCK_ << " | ";           //CO never zip LOCK or .LOCK (agl.LOCK) newly defined LOCK
    //aus << "grep -v SKIP | ";
    //aus << "grep -v " << KBIN_SUBDIRECTORIES << " | ";
    //aus << "grep -v aflow.in | grep -v " << _AFLOWIN_ << " ` " << endl; //CO, never zip aflow.in or _aflow.in (agl_aflow.in) or newly defined aflow.in
    //  cerr << aus.str() << endl;
    //aurostd::execute(aus);
    //DX and CO - END
  }
}

// *******************************************************************************************
namespace KBIN {
  void CompressDirectory(const _aflags& aflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    _kflags kflags;
    kflags.KZIP_BIN=XHOST.command("bzip2")+" -9q";
    KBIN::CompressDirectory(aflags,kflags);
  }
}

// *******************************************************************************************
// KBIN::Clean
// *******************************************************************************************
namespace KBIN {
  void Clean(const _aflags& aflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    KBIN::Clean(aflags.Directory);
  }
}

namespace KBIN {
  void Clean(const string _directory) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string directory=_directory;
    aurostd::StringSubst(directory,"/agl_aflow.in","");  // so it is easier to search
    aurostd::StringSubst(directory,"/ael_aflow.in","");  // so it is easier to search    
    aurostd::StringSubst(directory,"/"+_AFLOWIN_,"");  // so it is easier to search
    aurostd::StringSubst(directory,"/aflow.in","");  // so it is easier to search
    if(!aurostd::FileExist(string(directory+"/"+"NOCLEAN")) &&
       !aurostd::FileExist(string(directory+"/"+"NOCLEAR")) &&
       !aurostd::FileExist(string(directory+"/"+"noclean")) &&
       !aurostd::FileExist(string(directory+"/"+"noclear")) &&
       !aurostd::FileExist(string(directory+"/"+"Makefile")) &&
       !aurostd::FileExist(string(directory+"/"+"aflow.h")) &&
       !aurostd::FileExist(string(directory+"/"+"paper.tex")) &&
       !aurostd::FileExist(string(directory+"/"+"manuscript.tex")) &&
       !aurostd::FileExist(string(directory+"/"+"supplementary_information.tex")) &&
       !aurostd::FileExist(string(directory+"/"+"supplementary_materials.tex")) &&
       !aurostd::FileExist(string(directory+"/"+"review.tex"))) {
      if(aurostd::FileExist(string(directory+"/"+_AFLOWIN_)) ||    // normal aflow.in or specified it
	 aurostd::FileExist(string(directory+"/aflow.in")) ||      // normal aflow.in
	 aurostd::FileExist(string(directory+"/agl_aflow.in")) ||  // normal agl_aflow.in
	 aurostd::FileExist(string(directory+"/ael_aflow.in")) ) { // normal ael_aflow.in
	ostringstream aus;
	aurostd::StringstreamClean(aus);
	
	// CLEAN directory
        //DX and CO - START
        vector<string> vfiles;  //not only files, includes EVERYTHING
        string file_path;
        aurostd::DirectoryLS(directory,vfiles);
        for(uint i=0;i<vfiles.size();i++) {
          file_path=directory + "/" + vfiles.at(i);
          if(aurostd::substring2bool(vfiles.at(i),_AFLOWIN_)){continue;}
          if(aurostd::substring2bool(vfiles.at(i),"aflow.in")){continue;}
          if(aurostd::substring2bool(vfiles.at(i),"agl_aflow.in")){continue;}
          if(aurostd::substring2bool(vfiles.at(i),"ael_aflow.in")){continue;}
          if(aurostd::substring2bool(vfiles.at(i),_AFLOW_FROZSL_INPUT_FILE_)){continue;}
          if(aurostd::IsDirectory(file_path)){                 
            if(aurostd::substring2bool(vfiles.at(i),KBIN_SUBDIRECTORIES)){ // only directories we don't ignore
              aurostd::RemoveDirectory(file_path);
            }
            continue;                                                   // ignore all other directories
          }
          if(aurostd::IsFile(file_path)){                 
            if(vfiles.at(i).size() && vfiles.at(i)[0]=='.'){
              if(aurostd::substring2bool(vfiles.at(i),".nfs")){            // only hidden files we don't ignore
                aurostd::RemoveFile(file_path);
              }
              continue;                                                 // ignore all other hidden files
            }
            aurostd::RemoveFile(file_path);
          }
        }
	aurostd::execute("rm -f "+directory+"/.pam*");
	aurostd::execute("rm -f "+directory+"/*~");
	
	//DX and CO - END
	// now DECOMPRESS _AFLOWIN_
	ifstream FileCHECK;string FileNameCHECK;
	FileNameCHECK=directory+"/" + _AFLOWIN_ + ".bz2";                          // _AFLOWIN_.bz2
	FileCHECK.open(FileNameCHECK.c_str(),std::ios::in);                        // _AFLOWIN_.bz2
	FileCHECK.clear();FileCHECK.close();                                       // _AFLOWIN_.bz2
	if(FileCHECK) {                                                            // _AFLOWIN_.bz2
	  aus << XHOST.command("bzip2") << " -dq " << _AFLOWIN_ << ".bz2" << endl; // _AFLOWIN_.bz2
	  aurostd::execute(aus);                                                   // _AFLOWIN_.bz2
	}                                                                          // _AFLOWIN_.bz2
	FileNameCHECK=directory+"/"+_AFLOWIN_+".gz";                               // _AFLOWIN_.gz
	FileCHECK.open(FileNameCHECK.c_str(),std::ios::in);                        // _AFLOWIN_.gz
	FileCHECK.clear();FileCHECK.close();                                       // _AFLOWIN_.gz
	if(FileCHECK) {                                                            // _AFLOWIN_.gz
	  aus << XHOST.command("gzip") << " -d " << _AFLOWIN_ << ".gz" << endl;    // _AFLOWIN_.gz
	  aurostd::execute(aus);                                                   // _AFLOWIN_.gz
	}                                                                          // _AFLOWIN_.gz
      }
    }
  }
}

// *******************************************************************************************
// KBIN::XClean
// *******************************************************************************************
namespace KBIN {
  void XClean(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "KBIN::XClean: BEGIN" << endl;  
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=0) {
      init::ErrorOption(cout,options,"KBIN::XClean","aflow --xclean");
      exit(0);
    }
    
    vector<string> vcheck1;aurostd::string2tokens(string("OUTCAR.static,OUTCAR.relax2,OUTCAR.relax1,OUTCAR.relax1.bz2"),vcheck1,",");
    vector<string> vcheck2;aurostd::string2tokens(string("OUTCAR,OUTCAR,OUTCAR,OUTCAR.relax2.bz2"),vcheck2,",");
    
    vector<string> vfile;
    bool test=false;
    
    cout << "KBIN::XClean: checking missing " << "OUTCAR*" << " with " << _AFLOWLOCK_ << endl;  // check OUTCAR.static
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" ./ -name "+_AFLOWLOCK_),vfile);
    for(uint j=0;j<vfile.size();j++) {
      aurostd::StringSubst(vfile.at(j),_AFLOWLOCK_,"");
      if(!aurostd::FileExist(vfile.at(j)+"OUTCAR") && !aurostd::FileExist(vfile.at(j)+"OUTCAR.relax1.bz2")) {
	cout << "KBIN::XClean: cleaning=" << vfile.at(j) << endl;
	if(!test) KBIN::Clean(vfile.at(j));
      }
    }
    for(uint i=0;i<vcheck1.size();i++) {
      cout << "KBIN::XClean: checking missing " << vcheck2.at(i) << " with " << vcheck1.at(i) << endl;  // check OUTCAR.static
      aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" ./ -name "+vcheck1.at(i)),vfile);
      for(uint j=0;j<vfile.size();j++) {
	aurostd::StringSubst(vfile.at(j),vcheck1.at(i),"");
	if(!aurostd::FileExist(vfile.at(j)+vcheck2.at(i))) {
	  cout << "KBIN::XClean: cleaning=" << vfile.at(j) << endl;
	  if(!test) KBIN::Clean(vfile.at(j));
	}
      }
    }    
    if(LDEBUG) cerr << "KBIN::XClean: END" << endl;  
    // exit(0);
  }
} // namespace KBIN

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
