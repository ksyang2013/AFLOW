// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_INIT_CPP_
#define _AFLOW_INIT_CPP_
#include "aflow.h"

// ***************************************************************************

string _AFLOWIN_; 
string _AFLOWLOCK_; 

// THREADS
namespace AFLOW_PTHREADS {
  extern bool FLAG;        // run pthread YES/NO
  extern int MAX;         // how many MAX threads I can use  default or --np
  extern int RUNNING;      // how many threads are actually running
  extern pthread_t vpthread[MAX_ALLOCATABLE_PTHREADS];  // the actual thread
  extern int viret[MAX_ALLOCATABLE_PTHREADS];          // the thread runnings
  extern bool vpthread_busy[MAX_ALLOCATABLE_PTHREADS];  // is the thread busy
}

#define CPU_File     string("/proc/cpuinfo")

_XHOST XHOST;  // GLOBAL

// vector<string> vVASP_POTCAR_DIRECTORIES;
// vector<string> vAFLOW_LIBRARY_DIRECTORIES;
// vector<string> vAFLOW_PROJECTS_DIRECTORIES;

// uint LIBRARY_LIB2=LIBRARY_NOTHING;
// uint LIBRARY_ICSD=LIBRARY_NOTHING;
// uint LIBRARY_LIB3=LIBRARY_NOTHING;
// uint LIBRARY_LIB4=LIBRARY_NOTHING;
// uint LIBRARY_LIB5=LIBRARY_NOTHING;
// uint LIBRARY_LIB6=LIBRARY_NOTHING;
// uint LIBRARY_LIB7=LIBRARY_NOTHING;
// uint LIBRARY_LIB8=LIBRARY_NOTHING;
// uint LIBRARY_LIB9=LIBRARY_NOTHING;
// uint LIBRARY_LIB1=LIBRARY_NOTHING;

// ***************************************************************************
// init::InitMachine
// ***************************************************************************
namespace init {
    int GetCPUCores() { // CO 180124
        int ncpus=sysconf(_SC_NPROCESSORS_ONLN);
        if(ncpus==96) ncpus=48; // fix the hyperthreading lie
        if(ncpus<1) ncpus=1;
        return ncpus;
    }
    bool InitMachine(bool INIT_VERBOSE,vector<string>& argv,vector<string>& cmds,std::ostream& oss) {
        // DECLARATIONS
        bool LDEBUG=(FALSE || XHOST.DEBUG),found=FALSE;
        if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [BEGIN]" << endl;
        int depth_short=20,depth_long=45;
        string position;
        vector<string> tokens;
        if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;
        if(INIT_VERBOSE) oss << "* AFLOW V=" << string(AFLOW_VERSION) << " - machine information " << endl;
        if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;
        XHOST.argv.clear();for(uint i=0;i<argv.size();i++)  XHOST.argv.push_back(argv.at(i));

        // IMMEDIATELY
        XHOST.QUIET=aurostd::args2flag(argv,cmds,"--quiet|-q");
        XHOST.DEBUG=aurostd::args2flag(argv,cmds,"--debug");
        XHOST.TEST=aurostd::args2flag(argv,cmds,"--test|-test");
        XHOST.SKEW_TEST=aurostd::args2flag(argv,cmds,"--skew_test"); // DX 171025
        XHOST.MPI=aurostd::args2flag(argv,"--MPI|--mpi");
        XHOST.Tmpfs="/tmp/"; // to boot first thing to do.
        XHOST.User=aurostd::execute2string("whoami");  // AS SOON AS POSSIBLE
        XHOST.Home=aurostd::execute2string("cd && pwd");  // AS SOON AS POSSIBLE
        XHOST.GENERATE_AFLOWIN_ONLY=aurostd::args2flag(argv,cmds,"--generate_aflowin_only");  //CT 180719

        if(!aflowrc::is_available(oss,INIT_VERBOSE || XHOST.DEBUG)) aflowrc::write_default(oss,INIT_VERBOSE || XHOST.DEBUG);
        aflowrc::read(oss,INIT_VERBOSE || XHOST.DEBUG);
        XHOST.vflag_control.flag("AFLOWRC::READ",aurostd::args2flag(XHOST.argv,cmds,"--aflowrc=read|--aflowrc_read"));
        if(XHOST.vflag_control.flag("AFLOWRC::READ")) {aflowrc::print_aflowrc(oss,TRUE);exit(1);}

        // IMMEDIATELY GET PIDS
        XHOST.PID=getpid();XHOST.ostrPID.clear();XHOST.ostrPID.str(std::string());XHOST.ostrPID<<XHOST.PID;  // initialize PID ??
        XHOST.CPU_Cores=GetCPUCores();//sysconf(_SC_NPROCESSORS_ONLN);
        XHOST.CPU_Model="nan";
        XHOST.CPU_MHz="nan";
        // #ifndef _MACOSX_
        if(aurostd::FileExist(CPU_File)) {  // LINUX SYSTEMS
            // XHOST.CPU_Model
            aurostd::string2vectorstring(aurostd::execute2string("cat "+CPU_File+" | grep \"model name\""),tokens);
            if(tokens.size()>0) {
                aurostd::StringSubst(tokens.at(0),": ",":");aurostd::string2tokens(string(tokens.at(0)),tokens,":");
                if(tokens.size()>1) {
                    aurostd::StringSubst(tokens.at(1)," ","_");aurostd::StringSubst(tokens.at(1),"__","_");
                    aurostd::StringSubst(tokens.at(1),"__","_");aurostd::StringSubst(tokens.at(1),"__","_");
                    XHOST.CPU_Model=tokens.at(1);
                }
            }
            // XHOST.CPU_MHz
            aurostd::string2vectorstring(aurostd::execute2string("cat "+CPU_File+" | grep \"cpu MHz\""),tokens);
            if(tokens.size()>0) {
                aurostd::StringSubst(tokens.at(0),": ",":");aurostd::string2tokens(string(tokens.at(0)),tokens,":");
                if(tokens.size()>1) {
                    aurostd::StringSubst(tokens.at(1)," ","_");XHOST.CPU_MHz=aurostd::utype2string(ceil(aurostd::string2utype<double>(tokens.at(1))));
                }
            }
        }
        // MEMORY
        XHOST.RAM=init::GetRAM();
        XHOST.RAM_MB=XHOST.RAM/1024/1024;
        XHOST.RAM_GB=ceil(XHOST.RAM/1024/1024/1024);
        long int random_seed=aurostd::_random_initialize();
        XHOST.Tmpfs="/tmp/"; if(aurostd::FileExist("/run/shm/")) if(aurostd::DirectoryWritable("/run/shm/")) XHOST.Tmpfs="/run/shm/";
        if(INIT_VERBOSE) oss << "XHOST.Tmpfs=" << XHOST.Tmpfs << endl;
        // AFLOW_TIME
        XHOST.Day=aurostd::utype2string(aurostd::get_day());
        XHOST.Month=aurostd::utype2string(aurostd::get_month());
        XHOST.Year=aurostd::utype2string(aurostd::get_year());
        XHOST.Copyright_Years="2003-"+aurostd::utype2string(aurostd::get_year());
        XHOST.Time_starting=0.0;
        XHOST.Time_now=0.0;
        XHOST.Time_starting=aurostd::get_seconds();
        XHOST.Time_now=XHOST.Time_starting;
        // AFLOW_DATE
        XHOST.Date=aurostd::get_date();
        if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [5]" << endl;
        // some verbose
        if(INIT_VERBOSE) {
            oss << aurostd::PaddedPOST("XHOST.ostrPID = ",depth_short) << XHOST.ostrPID.str() << endl;
            oss << aurostd::PaddedPOST("DEFAULT_KZIP_BIN = ",depth_short) << DEFAULT_KZIP_BIN << endl;
            oss << aurostd::PaddedPOST("DEFAULT_KZIP_EXT = ",depth_short) << DEFAULT_KZIP_EXT << endl;
            oss << "--- MACHINE ------------ " << endl;
            oss << aurostd::PaddedPOST("CPU_Model = ",depth_short) << XHOST.CPU_Model << endl;
            oss << aurostd::PaddedPOST("CPU_MHz = ",depth_short) << XHOST.CPU_MHz << endl;
            oss << aurostd::PaddedPOST("CPU_Cores = ",depth_short) << XHOST.CPU_Cores << endl;
            oss << aurostd::PaddedPOST("RAM = ",depth_short) << (long int) XHOST.RAM << endl;
            oss << aurostd::PaddedPOST("RAM_MB = ",depth_short) << XHOST.RAM_MB << " (" << XHOST.RAM_GB << "GB)" << endl;
            oss << aurostd::PaddedPOST("random_seed = ",depth_short) << random_seed << endl;
        }

        // NAME and OS
        static struct utsname os;
        if((uname(&os)) < 0) { return FALSE;} // NULL
        XHOST.hostname=os.nodename;//aurostd::execute2string(string("hostname"));
        if(XHOST.hostname=="nietzsche") XHOST.hostname="nietzsche.mems.duke.edu";
        if(XHOST.hostname=="materials") XHOST.hostname="materials.duke.edu";
        if(XHOST.hostname=="aflowlib") XHOST.hostname="aflowlib.mems.duke.edu";
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("hostname = ",depth_short) << XHOST.hostname << endl;
        if(AFLOW_BlackList(XHOST.hostname)) {
            cout << "MMMMM  HOSTNAME BLACKLISTED = " << XHOST.hostname << endl;
            cerr << "MMMMM  HOSTNAME BLACKLISTED = " << XHOST.hostname << endl;
            exit(0);
        }

        // MACHINE TYPE
        XHOST.MachineType="linux";
        if(aurostd::substring2bool(os.sysname,"Darwin")) XHOST.MachineType="macosx";
        if(aurostd::substring2bool(os.sysname,"Linux")) XHOST.MachineType="linux";
        if(aurostd::substring2bool(os.sysname,"OSF1")) XHOST.MachineType="alpha";
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("machinetype = ",depth_short) << XHOST.MachineType << endl;
        // SERVER STUFF
        XHOST.AFLOW_MATERIALS_SERVER=AFLOW_MATERIALS_SERVER_DEFAULT;XHOST.AFLOW_WEB_SERVER=AFLOW_WEB_SERVER_DEFAULT; // DEFAULT
        if(aurostd::substring2bool(XHOST.hostname,"nietzsche")) {XHOST.AFLOW_MATERIALS_SERVER=AFLOW_MATERIALS_SERVER_DEFAULT;XHOST.AFLOW_WEB_SERVER=AFLOW_WEB_SERVER_DEFAULT;}
        if(aurostd::substring2bool(XHOST.hostname,"aflowlib")) {XHOST.AFLOW_MATERIALS_SERVER="aflowlib.mems.duke.edu";XHOST.AFLOW_WEB_SERVER="aflowlib.mems.duke.edu";}
        if(INIT_VERBOSE) {
            oss << "--- SERVER ------------------ " << endl;
            oss << aurostd::PaddedPOST("XHOST.AFLOW_MATERIALS_SERVER = ",depth_short) << XHOST.AFLOW_MATERIALS_SERVER << endl;
            oss << aurostd::PaddedPOST("XHOST.AFLOW_WEB_SERVER = ",depth_short) << XHOST.AFLOW_WEB_SERVER << endl;
        }
        // fix FIND --noleaf or not noleaf
        if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [8]" << endl;
        XHOST.Find_Parameters=DEFAULT_AFLOW_FIND_PARAMETERS_NOLEAF;
        if(XHOST.MachineType=="macosx") XHOST.Find_Parameters=DEFAULT_AFLOW_FIND_PARAMETERS_NORMAL;
        if(INIT_VERBOSE) {
            oss << "--- OS ------------------ " << endl;
            oss << aurostd::PaddedPOST("os.sysname = ",depth_short) << os.sysname << endl;
            oss << aurostd::PaddedPOST("os.release = ",depth_short) << os.release << endl;
            oss << aurostd::PaddedPOST("os.version = ",depth_short) << os.version << endl;
            oss << aurostd::PaddedPOST("os.machine = ",depth_short) << os.machine << endl;
            oss << aurostd::PaddedPOST("os.nodename = ",depth_short) << os.nodename << endl;
        }
        // USER //  uid_t uid;uid=geteuid();
        // OLD
        // struct passwd *pw; // pw=getpwuid(uid);
        // XHOST.User=string(pw->pw_name);
        // NEW
        if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [9]" << endl;
        // XHOST.User="boot";
        aurostd::StringSubst(XHOST.User,"\n","");
        // USER DONE
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("username = ",depth_short) << XHOST.User << endl;
        // GROUP
        aurostd::string2tokens(aurostd::execute2string("groups 2> /dev/null"),tokens);
        XHOST.Group="none";
        if(tokens.size()>0) XHOST.Group=tokens.at(0);
        aurostd::StringSubst(XHOST.Group,"\n","");
        // OLD
        // struct group *grp; // grp=getgrgid(pw->pw_gid);
        // XHOST.Group=string(grp->gr_name);
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("groupname = ",depth_short) << XHOST.Group << endl;
        // HOME DONE
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("home = ",depth_short) << XHOST.Home << endl;
        // SHELL
        if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [10]" << endl;
        XHOST.Shell=aurostd::execute2string(string("echo $SHELL"));
        aurostd::string2tokens(XHOST.Shell,tokens,"/");
        if(tokens.size()>0) XHOST.Shell=tokens.at(tokens.size()-1);
        aurostd::StringSubst(XHOST.Shell,"\n","");
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("shell = ",depth_short) << XHOST.Shell << endl;
        // PROGNAME
        if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [11]" << endl;
        XHOST.Progname="aflow";
        if(aurostd::substring2bool(XHOST.argv.at(0),"aconvasp") || aurostd::substring2bool(XHOST.argv.at(0),"convasp")) XHOST.Progname="aconvasp";
        if(aurostd::substring2bool(XHOST.argv.at(0),"aflow") || aurostd::substring2bool(XHOST.argv.at(0),"aflowd")) XHOST.Progname="aflow";
        if(aurostd::substring2bool(XHOST.argv.at(0),"aflow1") || aurostd::substring2bool(XHOST.argv.at(0),"aflowd1")) XHOST.Progname="aflow1";
        if(aurostd::substring2bool(XHOST.argv.at(0),"aflow2") || aurostd::substring2bool(XHOST.argv.at(0),"aflowd2")) XHOST.Progname="aflow2";
        if(aurostd::substring2bool(XHOST.argv.at(0),"apennsy") || aurostd::substring2bool(XHOST.argv.at(0),"apennsy")) XHOST.Progname="apennsy";
        if(INIT_VERBOSE) oss << "--- PROGNAME ------------ " << endl;
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("progname = ",depth_short) << XHOST.Progname << endl;

        // IP
        // ifconfig | grep inet | grep -v 127 | grep -v inet6
        // GET PROGRAMS
        if(INIT_VERBOSE) oss << "--- COMPILER ------------ " << endl;
        if(INIT_VERBOSE) oss << aurostd::PaddedPOST("G++ version = ",depth_short) << GCC_VERSION << endl;
        if(INIT_VERBOSE) oss << "--- BINARIES @ " << XHOST.hostname << " --- " << endl;
        // GET AFLOW_DATA
        found=FALSE;
        //CO 180706 - above we override Progname ./aflow with aflow, so it will never look in ./aflow
        //I am assuming this is for a good reason, but it screws things up here
        //so we need to rely on XHOST.argv.at(0) as above
        if(!found&&(XHOST.argv[0]=="aflow"||XHOST.argv[0]=="aflowd"||XHOST.argv[0]=="aconvasp"||XHOST.argv[0]=="apennsy")) {
            string aflow_data_command="aflow_data";
            //CO 180706 - note, IsCommandAvailableModify() simply overwrites aflow_data_command with true path from `which'
            if(aurostd::IsCommandAvailableModify(aflow_data_command)){XHOST.vcmd.push_back(aflow_data_command);found=TRUE;}  //CO 180703
        }
        if(!found&&(XHOST.argv[0]=="./aflow"||XHOST.argv[0]=="./aflowd"||XHOST.argv[0]=="./aconvasp"||XHOST.argv[0]=="./apennsy")) {
            if(aurostd::FileExist("./aflow_data")){XHOST.vcmd.push_back("./aflow_data");found=TRUE;}
        }
        if(!found&&(XHOST.argv[0]=="/usr/local/bin/aflow"||XHOST.argv[0]=="/usr/local/bin/aflowd"||XHOST.argv[0]=="/usr/local/bin/aconvasp"||XHOST.argv[0]=="/usr/local/bin/apennsy")) {
            if(aurostd::FileExist("/usr/local/bin/aflow_data")){XHOST.vcmd.push_back("/usr/local/bin/aflow_data");found=TRUE;}
        }
        if(!found) {
            aurostd::string2tokens(XHOST.argv.at(0),tokens,"/");  //CO 180703 - note, string2tokens without consecutive keeps beginning / with first entry
            string aflow_data="";
            for(uint i=0;i<tokens.size()-1;i++) {aflow_data+=tokens.at(i)+"/";} aflow_data+=string("aflow_data");
            if(aurostd::FileExist(aflow_data)) {
                XHOST.vcmd.push_back(aflow_data);
                for(uint i=0;i<XHOST.vcmd.size();i++)
                    if(aurostd::substring2bool(XHOST.vcmd.at(i),"aflow_data")) 
                        XHOST.vcmd.at(i)=aflow_data;
            }
        }
        // search for updates and proxies
        if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [14]" << endl;

        // if(XHOST.is_command("wget")) {aurostd::execute("wget -q http://materials.duke.edu/aflow_update/"+XHOST.User+"/"+XHOST.hostname);};

        // SOME LOADING UP
        if(INIT_VERBOSE) {
            if(XHOST.is_command("aflow_data")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"aflow_data\")=TRUE",depth_long) << "[" << XHOST.command("aflow_data") << "]" << endl;} else {oss << "XHOST.is_command(\"aflow_data\")=FALSE" << endl;}
            if(XHOST.is_command("beep")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"beep\")=TRUE",depth_long) << "[" << XHOST.command("beep") << "]" << endl;} else {oss << "XHOST.is_command(\"beep\")=FALSE" << endl;}
            if(XHOST.is_command("bkill")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"bkill\")=TRUE",depth_long) << "[" << XHOST.command("bkill") << "]" << endl;} else {oss << "XHOST.is_command(\"bkill\")=FALSE" << endl;}
            if(XHOST.is_command("bsub")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"bsub\")=TRUE",depth_long) << "[" << XHOST.command("bsub") << "]" << endl;} else {oss << "XHOST.is_command(\"bsub\")=FALSE" << endl;}
            if(XHOST.is_command("bzip2")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"bzip2\")=TRUE",depth_long) << "[" << "bzip2" << "]" << endl;} else {oss << "XHOST.is_command(\"bzip2\")=FALSE" << endl;}
            if(XHOST.is_command("compress")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"compress\")=TRUE",depth_long) << "[" << XHOST.command("compress") << "]" << endl;} else {oss << "XHOST.is_command(\"compress\")=FALSE" << endl;}
            if(XHOST.is_command("convert")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"convert\")=TRUE",depth_long) << "[" << XHOST.command("convert") << "]" << endl;} else {oss << "XHOST.is_command(\"convert\")=FALSE" << endl;}
            if(XHOST.is_command("curl")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"curl\")=TRUE",depth_long) << "[" << XHOST.command("curl") << "]" << endl;} else {oss << "XHOST.is_command(\"curl\")=FALSE" << endl;}
            if(XHOST.is_command("dvipdf")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"dvipdf\")=TRUE",depth_long) << "[" << XHOST.command("dvipdf") << "]" << endl;} else {oss << "XHOST.is_command(\"dvipdf\")=FALSE" << endl;}
            if(XHOST.is_command("dvips")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"dvips\")=TRUE",depth_long) << "[" << XHOST.command("dvips") << "]" << endl;} else {oss << "XHOST.is_command(\"dvips\")=FALSE" << endl;}
            if(XHOST.is_command("findsym")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"findsym\")=TRUE",depth_long) << "[" << XHOST.command("findsym") << "]" << endl;} else {oss << "XHOST.is_command(\"findsym\")=FALSE" << endl;}
            if(XHOST.is_command("frozsl")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"frozsl\")=TRUE",depth_long) << "[" << XHOST.command("frozsl") << "]" << endl;} else {oss << "XHOST.is_command(\"frozsl\")=FALSE" << endl;}
            if(XHOST.is_command("frozsl_init")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"frozsl_init\")=TRUE",depth_long) << "[" << XHOST.command("frozsl_init") << "]" << endl;} else {oss << "XHOST.is_command(\"frozsl_init\")=FALSE" << endl;}
            if(XHOST.is_command("gnuplot")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"gnuplot\")=TRUE",depth_long) << "[" << XHOST.command("gnuplot") << "]" << endl;} else {oss << "XHOST.is_command(\"gnuplot\")=FALSE" << endl;}
            if(XHOST.is_command("gzip")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"gzip\")=TRUE",depth_long) << "[" << "gzip" << "]" << endl;} else {oss << "XHOST.is_command(\"gzip\")=FALSE" << endl;}
            if(XHOST.is_command("halt")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"halt\")=TRUE",depth_long) << "[" << XHOST.command("halt") << "]" << endl;} else {oss << "XHOST.is_command(\"halt\")=FALSE" << endl;}
            if(XHOST.is_command("kill")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"kill\")=TRUE",depth_long) << "[" << XHOST.command("kill") << "]" << endl;} else {oss << "XHOST.is_command(\"kill\")=FALSE" << endl;}
            if(XHOST.is_command("ifconfig")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"ifconfig\")=TRUE",depth_long) << "[" << XHOST.command("ifconfig") << "]" << endl;} else {oss << "XHOST.is_command(\"ifconfig\")=FALSE" << endl;}
            if(XHOST.is_command("java")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"java\")=TRUE",depth_long) << "[" << XHOST.command("java") << "]" << endl;} else {oss << "XHOST.is_command(\"java\")=FALSE" << endl;}
            if(XHOST.is_command("jmol")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"jmol\")=TRUE",depth_long) << "[" << XHOST.command("jmol") << "]" << endl;} else {oss << "XHOST.is_command(\"jmol\")=FALSE" << endl;}
            if(XHOST.is_command("latex")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"latex\")=TRUE",depth_long) << "[" << XHOST.command("latex") << "]" << endl;} else {oss << "XHOST.is_command(\"latex\")=FALSE" << endl;}
            if(XHOST.is_command("matlab")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"matlab\")=TRUE",depth_long) << "[" << XHOST.command("matlab") << "]" << endl;} else {oss << "XHOST.is_command(\"matlab\")=FALSE" << endl;}
            if(XHOST.is_command("mpivasp46s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp46s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp46s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp46s\")=FALSE" << endl;}
            if(XHOST.is_command("mpivasp52s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp52s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp52s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp52s\")=FALSE" << endl;}
            if(XHOST.is_command("mpivasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp54s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp54s\")=FALSE" << endl;}
            if(XHOST.is_command("mpivasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp54s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp54s\")=FALSE" << endl;}
            if(XHOST.is_command("pdflatex")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"pdflatex\")=TRUE",depth_long) << "[" << XHOST.command("pdflatex") << "]" << endl;} else {oss << "XHOST.is_command(\"pdflatex\")=FALSE" << endl;}
            if(XHOST.is_command("platon")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"platon\")=TRUE",depth_long) << "[" << XHOST.command("platon") << "]" << endl;} else {oss << "XHOST.is_command(\"platon\")=FALSE" << endl;}
            if(XHOST.is_command("ps2pdf")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"ps2pdf\")=TRUE",depth_long) << "[" << XHOST.command("ps2pdf") << "]" << endl;} else {oss << "XHOST.is_command(\"ps2pdf\")=FALSE" << endl;}
            if(XHOST.is_command("pwd")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"pwd\")=TRUE",depth_long) << "[" << XHOST.command("pwd") << "]" << endl;} else {oss << "XHOST.is_command(\"pwd\")=FALSE" << endl;}
            if(XHOST.is_command("qconvex")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qconvex\")=TRUE",depth_long) << "[" << XHOST.command("qconvex") << "]" << endl;} else {oss << "XHOST.is_command(\"qconvex\")=FALSE" << endl;}
            if(XHOST.is_command("qdel")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qdel\")=TRUE",depth_long) << "[" << XHOST.command("qdel") << "]" << endl;} else {oss << "XHOST.is_command(\"qdel\")=FALSE" << endl;}
            if(XHOST.is_command("qhull")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qhull\")=TRUE",depth_long) << "[" << XHOST.command("qhull") << "]" << endl;} else {oss << "XHOST.is_command(\"qhull\")=FALSE" << endl;}
            if(XHOST.is_command("qsub")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qsub\")=TRUE",depth_long) << "[" << XHOST.command("qsub") << "]" << endl;} else {oss << "XHOST.is_command(\"qsub\")=FALSE" << endl;}
            if(XHOST.is_command("rasmol")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"rasmol\")=TRUE",depth_long) << "[" << XHOST.command("rasmol") << "]" << endl;} else {oss << "XHOST.is_command(\"rasmol\")=FALSE" << endl;}
            if(XHOST.is_command("rsync")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"rsync\")=TRUE",depth_long) << "[" << XHOST.command("rsync") << "]" << endl;} else {oss << "XHOST.is_command(\"rsync\")=FALSE" << endl;}
            if(XHOST.is_command("sbatch")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"sbatch\")=TRUE",depth_long) << "[" << XHOST.command("sbatch") << "]" << endl;} else {oss << "XHOST.is_command(\"sbatch\")=FALSE" << endl;}
            if(XHOST.is_command("scancel")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"scancel\")=TRUE",depth_long) << "[" << XHOST.command("scancel") << "]" << endl;} else {oss << "XHOST.is_command(\"scancel\")=FALSE" << endl;}
            if(XHOST.is_command("sensors")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"sensors\")=TRUE",depth_long) << "[" << XHOST.command("sensors") << "]" << endl;} else {oss << "XHOST.is_command(\"sensors\")=FALSE" << endl;}
            if(XHOST.is_command("unzip")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"unzip\")=TRUE",depth_long) << "[" << XHOST.command("unzip") << "]" << endl;} else {oss << "XHOST.is_command(\"unzip\")=FALSE" << endl;}
            if(XHOST.is_command("vasp46s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp46s\")=TRUE",depth_long) << "[" << XHOST.command("vasp46s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp46s\")=FALSE" << endl;}
            if(XHOST.is_command("vasp52s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp52s\")=TRUE",depth_long) << "[" << XHOST.command("vasp52s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp52s\")=FALSE" << endl;}
            if(XHOST.is_command("vasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp54s\")=TRUE",depth_long) << "[" << XHOST.command("vasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp54s\")=FALSE" << endl;}
            if(XHOST.is_command("vasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp54s\")=TRUE",depth_long) << "[" << XHOST.command("vasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp54s\")=FALSE" << endl;}
            if(XHOST.is_command("wget")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"wget\")=TRUE",depth_long) << "[" << XHOST.command("wget") << "]" << endl;} else {oss << "XHOST.is_command(\"wget\")=FALSE" << endl;}
            if(XHOST.is_command("zip")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"zip\")=TRUE",depth_long) << "[" << XHOST.command("zip") << "]" << endl;} else {oss << "XHOST.is_command(\"zip\")=FALSE" << endl;}
            if(XHOST.is_command("xz")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"xz\")=TRUE",depth_long) << "[" << "xz" << "]" << endl;} else {oss << "XHOST.is_command(\"xz\")=FALSE" << endl;}
            if(XHOST.is_command("xxx")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"xxx\")=TRUE",depth_long) << "[" << XHOST.command("xxx") << "]" << endl;} else {oss << "XHOST.is_command(\"xxx\")=FALSE" << endl;}
        }
        // #include "aflow_init_aus.cpp"
        // TEMPERATURE
        XHOST.sensors_allowed=TRUE;
        if(0)  if(XHOST.is_command("sensors")) {
            init::GetTEMPs();
            if(INIT_VERBOSE) oss << "XHOST.vTemperatureCore.size()=" << XHOST.vTemperatureCore.size() << endl;
            if(INIT_VERBOSE && XHOST.vTemperatureCore.size()) {
                oss << "--- TEMPERATURES --------------- " << endl;
                for(uint i=0;i<XHOST.vTemperatureCore.size();i++) {oss << (i==0?"TEMP(C)=[":"") << XHOST.vTemperatureCore.at(i) << (i<XHOST.vTemperatureCore.size()-1?",":"]\n");}
                //;if(i<XHOST.vTemperatureCore.size()-1) oss << ","; else oss << "]" << endl;}
        }
    }
    // QUEUE STUFF
    XHOST.is_PBS=FALSE;
    XHOST.PBS_NUM_PPN=aurostd::getenv2uint("PBS_NUM_PPN");
    XHOST.PBS_NNODES=aurostd::getenv2uint("PBS_NNODES");
    if(XHOST.PBS_NUM_PPN!=0 || XHOST.PBS_NNODES!=0) XHOST.is_PBS=TRUE;
    XHOST.is_SLURM=FALSE;
    XHOST.SLURM_CPUS_ON_NODE=aurostd::getenv2int("SLURM_CPUS_ON_NODE");
    XHOST.SLURM_NNODES=aurostd::getenv2int("SLURM_NNODES");          
    if(XHOST.SLURM_CPUS_ON_NODE!=0 || XHOST.SLURM_NNODES!=0) XHOST.is_SLURM=TRUE;
    if(INIT_VERBOSE) {
        oss << "--- QUEUES --------------- " << endl;
        oss << "is_PBS=" << XHOST.is_PBS << endl;
        oss << "PBS_NUM_PPN=" << XHOST.PBS_NUM_PPN << endl;
        oss << "PBS_NNODES=" << XHOST.PBS_NNODES << endl;
        oss << "is_SLURM=" << XHOST.is_SLURM << endl;
        oss << "SLURM_CPUS_ON_NODE=" << XHOST.SLURM_CPUS_ON_NODE << endl;
        oss << "SLURM_NNODES=" << XHOST.SLURM_NNODES << endl;
    }
    // maxmem
    XHOST.maxmem=aurostd::args2attachedutype<double>(XHOST.argv,"--mem=|--maxmem=",101.0);
    // DEBUG  cerr << "XHOST.maxmem=" << XHOST.maxmem << endl;exit(0);

    // MPIs
    if(INIT_VERBOSE) {
        oss << "--- DATES --------------- " << endl;
        oss << "XHOST.Date = " << XHOST.Date << endl;
        oss << "XHOST.maxmem = " << XHOST.maxmem << endl;
        oss << "--- MPI COMMANDS --- " << endl;
        oss << "MPI_OPTIONS_DUKE_BETA_MPICH=" << MPI_OPTIONS_DUKE_BETA_MPICH << "\"" << endl;
        oss << "MPI_COMMAND_DUKE_BETA_MPICH=" << MPI_COMMAND_DUKE_BETA_MPICH << "\"" << endl;
        oss << "MPI_BINARY_DIR_DUKE_BETA_MPICH=" << MPI_BINARY_DIR_DUKE_BETA_MPICH << "\"" << endl;
        oss << "MPI_OPTIONS_DUKE_BETA_OPENMPI=" << MPI_OPTIONS_DUKE_BETA_OPENMPI << "\"" << endl;
        oss << "MPI_COMMAND_DUKE_BETA_OPENMPI=" << MPI_COMMAND_DUKE_BETA_OPENMPI << "\"" << endl;
        oss << "MPI_BINARY_DIR_DUKE_BETA_OPENMPI=" << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << "\"" << endl;
        oss << "MPI_OPTIONS_DUKE_QRATS_MPICH=" << MPI_OPTIONS_DUKE_QRATS_MPICH << "\"" << endl;
        oss << "MPI_COMMAND_DUKE_QRATS_MPICH=" << MPI_COMMAND_DUKE_QRATS_MPICH << "\"" << endl;
        oss << "MPI_BINARY_DIR_DUKE_QRATS_MPICH=" << MPI_BINARY_DIR_DUKE_QRATS_MPICH << "\"" << endl;
        // CO - START
        oss << "MPI_OPTIONS_DUKE_QFLOW_OPENMPI=" << MPI_OPTIONS_DUKE_QFLOW_OPENMPI << "\"" << endl;
        oss << "MPI_COMMAND_DUKE_QFLOW_OPENMPI=" << MPI_COMMAND_DUKE_QFLOW_OPENMPI << "\"" << endl;
        oss << "MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI=" << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << "\"" << endl;
        oss << "MPI_OPTIONS_MPCDF_EOS=" << MPI_OPTIONS_MPCDF_EOS << "\"" << endl;
        oss << "MPI_COMMAND_MPCDF_EOS=" << MPI_COMMAND_MPCDF_EOS << "\"" << endl;
        oss << "MPI_NCPUS_MPCDF_EOS=" << MPI_NCPUS_MPCDF_EOS << "\"" << endl;
        oss << "MPI_HYPERTHREADING_MPCDF_EOS=" << MPI_HYPERTHREADING_MPCDF_EOS << "\"" << endl;
        oss << "MPI_BINARY_DIR_MPCDF_EOS=" << MPI_BINARY_DIR_MPCDF_EOS << "\"" << endl;
        oss << "MPI_OPTIONS_MPCDF_DRACO=" << MPI_OPTIONS_MPCDF_DRACO << "\"" << endl;
        oss << "MPI_COMMAND_MPCDF_DRACO=" << MPI_COMMAND_MPCDF_DRACO << "\"" << endl;
        oss << "MPI_NCPUS_MPCDF_DRACO=" << MPI_NCPUS_MPCDF_DRACO << "\"" << endl;
        oss << "MPI_HYPERTHREADING_MPCDF_DRACO=" << MPI_HYPERTHREADING_MPCDF_DRACO << "\"" << endl;
        oss << "MPI_BINARY_DIR_MPCDF_DRACO=" << MPI_BINARY_DIR_MPCDF_DRACO << "\"" << endl;
        oss << "MPI_OPTIONS_MPCDF_COBRA=" << MPI_OPTIONS_MPCDF_COBRA << "\"" << endl;
        oss << "MPI_COMMAND_MPCDF_COBRA=" << MPI_COMMAND_MPCDF_COBRA << "\"" << endl;
        oss << "MPI_NCPUS_MPCDF_COBRA=" << MPI_NCPUS_MPCDF_COBRA << "\"" << endl;
        oss << "MPI_HYPERTHREADING_MPCDF_COBRA=" << MPI_HYPERTHREADING_MPCDF_COBRA << "\"" << endl;
        oss << "MPI_BINARY_DIR_MPCDF_COBRA=" << MPI_BINARY_DIR_MPCDF_COBRA << "\"" << endl;
        oss << "MPI_OPTIONS_MPCDF_HYDRA=" << MPI_OPTIONS_MPCDF_HYDRA << "\"" << endl;
        oss << "MPI_COMMAND_MPCDF_HYDRA=" << MPI_COMMAND_MPCDF_HYDRA << "\"" << endl;
        oss << "MPI_NCPUS_MPCDF_HYDRA=" << MPI_NCPUS_MPCDF_HYDRA << "\"" << endl;
        oss << "MPI_HYPERTHREADING_MPCDF_HYDRA=" << MPI_HYPERTHREADING_MPCDF_HYDRA << "\"" << endl;
        oss << "MPI_BINARY_DIR_MPCDF_HYDRA=" << MPI_BINARY_DIR_MPCDF_HYDRA << "\"" << endl;
        // CO - END
        //DX 5/2/18 - DoD CONRAD - START
        oss << "MPI_OPTIONS_DOD_CONRAD=" << MPI_OPTIONS_DOD_CONRAD << "\"" << endl;
        oss << "MPI_COMMAND_DOD_CONRAD=" << MPI_COMMAND_DOD_CONRAD << "\"" << endl;
        oss << "MPI_BINARY_DIR_DOD_CONRAD=" << MPI_BINARY_DIR_DOD_CONRAD << "\"" << endl;
        //DX 5/2/18 - DoD CONRAD - END
        oss << "MPI_OPTIONS_DUKE_MATERIALS=" << MPI_OPTIONS_DUKE_MATERIALS << "\"" << endl;
        oss << "MPI_COMMAND_DUKE_MATERIALS=" << MPI_COMMAND_DUKE_MATERIALS << "\"" << endl;
        oss << "MPI_BINARY_DIR_DUKE_MATERIALS=" << MPI_BINARY_DIR_DUKE_MATERIALS << "\"" << endl;
        oss << "MPI_OPTIONS_DUKE_AFLOWLIB=" << MPI_OPTIONS_DUKE_AFLOWLIB << "\"" << endl;
        oss << "MPI_COMMAND_DUKE_AFLOWLIB=" << MPI_COMMAND_DUKE_AFLOWLIB << "\"" << endl;
        oss << "MPI_BINARY_DIR_DUKE_AFLOWLIB=" << MPI_BINARY_DIR_DUKE_AFLOWLIB << "\"" << endl;
        oss << "MPI_OPTIONS_FULTON_MARYLOU=" << MPI_OPTIONS_FULTON_MARYLOU << "\"" << endl;
        oss << "MPI_COMMAND_FULTON_MARYLOU=" << MPI_COMMAND_FULTON_MARYLOU << "\"" << endl;
        oss << "MPI_BINARY_DIR_FULTON_MARYLOU=" << MPI_BINARY_DIR_FULTON_MARYLOU << "\"" << endl;
        oss << "MPI_OPTIONS_TRINITY_PARSONS=" << MPI_OPTIONS_TRINITY_PARSONS << "\"" << endl;
        oss << "MPI_COMMAND_TRINITY_PARSONS=" << MPI_COMMAND_TRINITY_PARSONS << "\"" << endl;
        oss << "MPI_BINARY_DIR_TRINITY_PARSONS=" << MPI_BINARY_DIR_TRINITY_PARSONS << "\"" << endl;
        oss << "MPI_OPTIONS_TERAGRID_RANGER=" << MPI_OPTIONS_TERAGRID_RANGER << "\"" << endl;
        oss << "MPI_COMMAND_TERAGRID_RANGER=" << MPI_COMMAND_TERAGRID_RANGER << "\"" << endl;
        oss << "MPI_BINARY_DIR_TERAGRID_RANGER=" << MPI_BINARY_DIR_TERAGRID_RANGER << "\"" << endl;
        oss << "MPI_OPTIONS_TERAGRID_KRAKEN=" << MPI_OPTIONS_TERAGRID_KRAKEN << "\"" << endl;
        oss << "MPI_COMMAND_TERAGRID_KRAKEN=" << MPI_COMMAND_TERAGRID_KRAKEN << "\"" << endl;
        oss << "MPI_BINARY_DIR_TERAGRID_KRAKEN=" << MPI_BINARY_DIR_TERAGRID_KRAKEN << "\"" << endl;
        oss << "MPI_OPTIONS_MACHINE1=" << MPI_OPTIONS_MACHINE1 << "\"" << endl;
        oss << "MPI_COMMAND_MACHINE1=" << MPI_COMMAND_MACHINE1 << "\"" << endl;
        oss << "MPI_BINARY_DIR_MACHINE1=" << MPI_BINARY_DIR_MACHINE1 << "\"" << endl;
        oss << "MPI_OPTIONS_MACHINE2=" << MPI_OPTIONS_MACHINE2 << "\"" << endl;
        oss << "MPI_COMMAND_MACHINE2=" << MPI_COMMAND_MACHINE2 << "\"" << endl;
        oss << "MPI_BINARY_DIR_MACHINE2=" << MPI_BINARY_DIR_MACHINE2 << "\"" << endl;
        //   oss << endl;
    }
    // DO LISRS
    vector<string> vstrs;
    // DO VARIABLES
    aurostd::string2tokens(DEFAULT_VASP_POTCAR_DIRECTORIES,vVASP_POTCAR_DIRECTORIES,",");// vVASP_POTCAR_DIRECTORIES;
    // for(uint i=0;i<vVASP_POTCAR_DIRECTORIES.size();i++) oss << "vVASP_POTCAR_DIRECTORIES.at(" << i << ")=" << vVASP_POTCAR_DIRECTORIES.at(i) << endl; // exit(0);
    // LIBRARIES
    aurostd::string2tokens(DEFAULT_AFLOW_LIBRARY_DIRECTORIES,vAFLOW_LIBRARY_DIRECTORIES,",");// vAFLOW_LIBRARY_DIRECTORIES;
    //for(uint i=0;i<vAFLOW_LIBRARY_DIRECTORIES.size();i++) oss << vAFLOW_LIBRARY_DIRECTORIES.at(i) << endl;exit(0);
    // PROJECTS
    vAFLOW_PROJECTS_DIRECTORIES.clear();
    // XHOST_LIBRARY_LIB2=LIBRARY_NOTHING;
    aurostd::string2tokens(DEFAULT_AFLOW_PROJECTS_DIRECTORIES,vstrs,",");// vAFLOW_PROJECTS_DIRECTORIES;
    // DEBUG  cerr << "vAFLOW_PROJECTS_DIRECTORIES.size()=" << vAFLOW_PROJECTS_DIRECTORIES.size() << endl; 
    for(uint i=0;i<vstrs.size();i++) { // oss << vstrs.at(i) << endl;
        if(aurostd::FileExist(vstrs.at(i)) && aurostd::FileExist(vstrs.at(i)+"/LIB")) {
            if(aurostd::substring2bool(vstrs.at(i),"ICSD")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_ICSD=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB1")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB1=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB2")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB2=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB3")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB3=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB4")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB4=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB5")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB5=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB6")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB6=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB7")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB7=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB8")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB8=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
            if(aurostd::substring2bool(vstrs.at(i),"LIB9")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB9=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
        }
    }
    // DEBUG cerr << "vAFLOW_PROJECTS_DIRECTORIES.size()=" << vAFLOW_PROJECTS_DIRECTORIES.size() << endl;
    // DEBUG   for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++) cerr << "vAFLOW_PROJECTS_DIRECTORIES.at(i)=" << vAFLOW_PROJECTS_DIRECTORIES.at(i) << endl;

    if(INIT_VERBOSE) {
        oss << "--- PROJECTS @ " << XHOST.hostname << " --- " << endl;
        if(XHOST_LIBRARY_ICSD!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_ICSD=" << XHOST_LIBRARY_ICSD << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_ICSD << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD) << endl;
        if(XHOST_LIBRARY_LIB1!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB1=" << XHOST_LIBRARY_LIB1 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB1 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1) << endl;
        if(XHOST_LIBRARY_LIB2!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB2=" << XHOST_LIBRARY_LIB2 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB2 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2) << endl;
        if(XHOST_LIBRARY_LIB3!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB3=" << XHOST_LIBRARY_LIB3 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB3 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3) << endl;
        if(XHOST_LIBRARY_LIB4!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB4=" << XHOST_LIBRARY_LIB4 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB4 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4) << endl;
        if(XHOST_LIBRARY_LIB5!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB5=" << XHOST_LIBRARY_LIB5 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB5 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5) << endl;
        if(XHOST_LIBRARY_LIB6!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB6=" << XHOST_LIBRARY_LIB6 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB6 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6) << endl;
        if(XHOST_LIBRARY_LIB7!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB7=" << XHOST_LIBRARY_LIB7 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB7 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7) << endl;
        if(XHOST_LIBRARY_LIB8!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB8=" << XHOST_LIBRARY_LIB8 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB8 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8) << endl;
        if(XHOST_LIBRARY_LIB9!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB9=" << XHOST_LIBRARY_LIB9 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB9 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9) << endl;

        if(XHOST_LIBRARY_ICSD!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_ICSD_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_ICSD_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB1!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB1_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB1_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB2!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB2_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB2_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB3!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB3_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB3_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB4!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB4_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB4_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB5!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB5_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB5_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB6!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB6_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB6_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB7!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB7_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB7_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB8!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB8_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB8_LIB.size()=" << tokens.size() << endl;
        }
        if(XHOST_LIBRARY_LIB9!=LIBRARY_NOTHING) {
            aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB9_LIB"),tokens,"\n");
            oss << "Library_CALCULATED_LIB9_LIB.size()=" << tokens.size() << endl;
        }
    }

    // for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++) oss << vAFLOW_PROJECTS_DIRECTORIES.at(i) << endl;exit(0);
    // OLD aurostd::string2tokens(string(AFLOW_PROJECTS_DIRECTORIES),vAFLOW_PROJECTS_DIRECTORIES,",");// vAFLOW_PROJECTS_DIRECTORIES;
    // OLD for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++) oss << vAFLOW_PROJECTS_DIRECTORIES.at(i) << endl;exit(0);

    // check for MACHINES MARYLOU
    XHOST.is_MACHINE_FULTON_MARYLOU=FALSE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(XHOST.User=="fslcollab8" || XHOST.User=="glh43" || XHOST.User=="legoses") XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(aurostd::substring2bool(XHOST.Group,"fslcollab")) XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(aurostd::substring2bool(XHOST.Group,"fslg")) XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(aurostd::substring2bool(XHOST.Group,"glh43")) XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    // some other technique to get MARYLOU

    // check for APENNSY_USE_SERVER/AFLOWLIB
    if(XHOST.hostname=="nietzsche.mems.duke.edu" || XHOST.hostname=="materials.duke.edu" || XHOST.hostname=="aaaaflowlib.mems.duke.edu") {
        XHOST.APENNSY_USE_SERVER=TRUE;XHOST.APENNSY_USE_LIBRARY=FALSE;XHOST.APENNSY_SERVER_AFLOWLIB_ORG=FALSE;
    } else {
        XHOST.APENNSY_USE_SERVER=FALSE;XHOST.APENNSY_USE_LIBRARY=FALSE;XHOST.APENNSY_SERVER_AFLOWLIB_ORG=TRUE;
    }
    if(INIT_VERBOSE) {
        oss << "--- APENNSY_USE_*** @ " << XHOST.hostname << " --- " << endl;
        oss << "XHOST.APENNSY_USE_SERVER=" << XHOST.APENNSY_USE_SERVER << endl;
        oss << "XHOST.APENNSY_USE_LIBRARY=" << XHOST.APENNSY_USE_LIBRARY << endl;
        oss << "XHOST.APENNSY_SERVER_AFLOWLIB_ORG=" << XHOST.APENNSY_SERVER_AFLOWLIB_ORG << endl;
    }    
    // DO aflow.in and LOCK
    if(INIT_VERBOSE) oss << "--- LOADING @ _AFLOWIN_ and _AFLOWLOCK_ --- " << endl;
    _AFLOWIN_=aurostd::args2attachedstring(XHOST.argv,"--use_aflow.in=","aflow.in");
    _AFLOWLOCK_=aurostd::args2attachedstring(XHOST.argv,"--use_LOCK=","LOCK");
    if(INIT_VERBOSE) oss << "_AFLOWIN_=" << _AFLOWIN_ << endl;
    if(INIT_VERBOSE) oss << "_AFLOWLOCK_=" << _AFLOWLOCK_ << endl;

    // LOAD control stuff
    if(INIT_VERBOSE) oss << "--- LOADING @ control --- " << endl;
    XHOST.vflag_control.flag("MACHINE",aurostd::args2flag(argv,cmds,"--machine|--machine="));
    XHOST.vflag_control.flag("MULTI=SH",aurostd::args2flag(argv,cmds,"--multi=sh"));
    XHOST.vflag_control.flag("MULTI=BZIP2",aurostd::args2flag(argv,cmds,"--multi=bzip2"));
    XHOST.vflag_control.flag("MULTI=BUNZIP2",aurostd::args2flag(argv,cmds,"--multi=bunzip2"));
    XHOST.vflag_control.flag("MULTI=GZIP",aurostd::args2flag(argv,cmds,"--multi=gzip"));
    XHOST.vflag_control.flag("MULTI=GUNZIP",aurostd::args2flag(argv,cmds,"--multi=gunzip"));
    XHOST.vflag_control.flag("MULTI=XZIP",aurostd::args2flag(argv,cmds,"--multi=xz|--multi=xzip"));
    XHOST.vflag_control.flag("MULTI=XUNZIP",aurostd::args2flag(argv,cmds,"--multi=xunzip"));
    XHOST.vflag_control.flag("MULTI=BZ2XZ",aurostd::args2flag(argv,cmds,"--multi=bz2xz"));
    XHOST.vflag_control.flag("MULTI=GZ2XZ",aurostd::args2flag(argv,cmds,"--multi=gz2xz"));
    XHOST.vflag_control.flag("MULTI=ZIP",aurostd::args2flag(argv,cmds,"--multi=zip"));
    XHOST.vflag_control.flag("MONITOR",aurostd::args2flag(argv,cmds,"--monitor"));
    XHOST.vflag_control.flag("GETTEMP",aurostd::args2flag(argv,cmds,"--getTEMP|--getTEMPS|--getTEMPs|--gettemp|--gettemps"));
    XHOST.vflag_control.flag("SWITCH_AFLOW",
            aurostd::args2flag(argv,cmds,"--run|--clean|--xclean|--multi|--generate") ||
            aurostd::args2attachedflag(argv,cmds,"--run=") ||
            aurostd::args2flag(argv,cmds,"--generate_vasp_from_aflowin|--generate_aflowin_from_vasp"));
    // DX - START
    XHOST.vflag_control.flag("AFLOWIN_SYM",aurostd::args2flag(argv,cmds,"--generate_symmetry|--generate_sym")); // DX
    // DX - END
    // [OBSOLETE]    XHOST.vflag_control.flag("SWITCH_AFLOWLIB",aurostd::args2flag(argv,cmds,"--aflowlib") || aurostd::args2attachedflag(argv,cmds,"--aflowlib="));
    XHOST.vflag_control.flag("SWITCH_APENNSY1",aurostd::args2flag(argv,cmds,"--apennsy|--lib2|--lib2u|--lib2pgm|--LIB2|--LIB2U|--LIB2PGM|--libraryX|--libraryU|--libraryPGM|--alloy"));
    XHOST.vflag_control.flag("SWITCH_APENNSY2",aurostd::args2flag(argv,cmds,"--apool|--apool_private|--apool_test|--library2|-simpls|--simpls|--VASPIN|--energy|--psenergy"));

    // DIRECTORY NEEDS A SPECIAL TREATMENT
    found=FALSE;

    string dir_default="./",dir=dir_default;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--DIRECTORY|--directory|--D|--d"))) {dir=aurostd::args2string(argv,"--DIRECTORY|--directory|--D|--d",dir_default);if(INIT_VERBOSE) cerr << "--DIRECTORY " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"./"))) {dir="./";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"."))) {dir=".";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"../"))) {dir="../";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"~/"))) {dir="~/";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"~"))) {dir="~";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--DIRECTORY=|--D=|--d="))) {dir=aurostd::args2attachedstring(argv,"--DIRECTORY=|--D=|--d=",dir_default);if(INIT_VERBOSE) cerr << "--DIRECTORY=" << dir << " " << endl;}
    XHOST.vflag_control.flag("DIRECTORY",found);  // if found
    if(XHOST.vflag_control.flag("DIRECTORY")) XHOST.vflag_control.push_attached("DIRECTORY",dir);
    if(XHOST.vflag_control.flag("DIRECTORY")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"DIR\")=[" << XHOST.vflag_control.getattachedscheme("DIRECTORY") << "]" << endl; 
    // FILE NEEDS A SPECIAL TREATMENT
    string file_default="xxxx",file=file_default;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--FILE|--file|--F|--f"))) {file=aurostd::args2string(argv,"--FILE|--file|--F|--f",file_default);if(INIT_VERBOSE) cerr << "--FILE " << file << " " << endl;}
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--FILE=|--file=|--F=|--f="))) {file=aurostd::args2attachedstring(argv,"--FILE=|--file=|--F=|--f=",file_default);if(INIT_VERBOSE) cerr << "--FILE=" << file << " " << endl;}
    XHOST.vflag_control.flag("FILE",found);  // if found
    if(XHOST.vflag_control.flag("FILE")) XHOST.vflag_control.push_attached("FILE",file);
    if(XHOST.vflag_control.flag("FILE")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"FILE\")=[" << XHOST.vflag_control.getattachedscheme("FILE") << "]" << endl; 
    // VFILES NEEDS A SPECIAL TREATMENT
    string dirs_default="xxxx",dirs=dirs_default;
    vector<string> vdir;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--DIRECTORY|--directory|--D|--d"))) {dirs="";vdir=aurostd::args2vectorstring(argv,"--DIRECTORY|--directory|--D|--d",dirs_default);if(INIT_VERBOSE) cerr << "--DIRECTORY " << vdir.size() << " " << endl;}
    if(found) for(uint i=0;i<vdir.size();i++) dirs+=vdir.at(i)+(i<vdir.size()-1?",":""); // glue
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--DIRECTORY=|--D=|--d="))) {dirs=aurostd::args2attachedstring(argv,"--DIRECTORY=|--D=|--d=",dirs_default);if(INIT_VERBOSE) cerr << "--DIRECTORY=" << dirs << " " << endl;}
    XHOST.vflag_control.flag("VDIR",found);  // if found
    if(XHOST.vflag_control.flag("VDIR")) XHOST.vflag_control.push_attached("VDIR",dirs); 
    if(XHOST.vflag_control.flag("VDIR")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"VDIR\")=[" << XHOST.vflag_control.getattachedscheme("VDIR") << "]" << endl; 
    // VFILES NEEDS A SPECIAL TREATMENT
    string files_default="xxxx",files=files_default;
    vector<string> vfile;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--FILE|--file|--F|--f"))) {files="";vfile=aurostd::args2vectorstring(argv,"--FILE|--file|--F|--f",files_default);if(INIT_VERBOSE) cerr << "--FILE " << vfile.size() << " " << endl;}
    if(found) for(uint i=0;i<vfile.size();i++) files+=vfile.at(i)+(i<vfile.size()-1?",":""); // glue
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--FILE=|--file=|--F=|--f="))) {files=aurostd::args2attachedstring(argv,"--FILE=|--file=|--F=|--f=",files_default);if(INIT_VERBOSE) cerr << "--FILE=" << files << " " << endl;}
    XHOST.vflag_control.flag("VFILES",found);  // if found
    if(XHOST.vflag_control.flag("VFILES")) XHOST.vflag_control.push_attached("VFILES",files); 
    if(XHOST.vflag_control.flag("VFILES")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"VFILES\")=[" << XHOST.vflag_control.getattachedscheme("VFILES") << "]" << endl; 
    // exit(0); 

    XHOST.vflag_control.flag("AFLOW_HELP",aurostd::args2flag(argv,cmds,"-h|--help"));
    XHOST.vflag_control.flag("AFLOW_EXCEPTIONS", aurostd::args2flag(argv, cmds, "-e|--errors|--exceptions"));  // ME180531
    XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3",aurostd::args2flag(argv,cmds,"-l|--license"));
    XHOST.vflag_control.flag("README_AFLOW",aurostd::args2flag(argv,cmds,"--readme=run|--readme=aflow|--readme_aflow"));
    XHOST.vflag_control.flag("README_AFLOW_PFLOW",aurostd::args2flag(argv,cmds,"--readme=pflow|--readme=processor|--readme=aconvasp|--readme_aconvasp"));
    XHOST.vflag_control.flag("README_FROZSL",aurostd::args2flag(argv,cmds,"--readme=frozsl|--readme_frozsl"));
    XHOST.vflag_control.flag("README_APL",aurostd::args2flag(argv,cmds,"--readme=apl|--readme_apl"));
    XHOST.vflag_control.flag("README_QHA",aurostd::args2flag(argv,cmds,"--readme=qha|--readme_qha|--readme=qha3p|--readme_qha3p|--readme=scqha|--readme_scqha"));
    XHOST.vflag_control.flag("README_AAPL",aurostd::args2flag(argv,cmds,"--readme=aapl|--readme_aapl"));
    XHOST.vflag_control.flag("README_AGL",aurostd::args2flag(argv,cmds,"--readme=agl|--readme_agl"));
    XHOST.vflag_control.flag("README_AEL",aurostd::args2flag(argv,cmds,"--readme=ael|--readme_ael"));
    XHOST.vflag_control.flag("README_ANRL",aurostd::args2flag(argv,cmds,"--readme=anrl|--readme_anrl"));
    XHOST.vflag_control.flag("README_COMPARE",aurostd::args2flag(argv,cmds,"--readme=compare|--readme_compare"));
    XHOST.vflag_control.flag("README_SYMMETRY",aurostd::args2flag(argv,cmds,"--readme=symmetry|--readme_symmetry"));
    XHOST.vflag_control.flag("README_CHULL",aurostd::args2flag(argv,cmds,"--readme=chull|--readme_chull"));
    XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION",aurostd::args2flag(argv,cmds,"--readme=partial_occupation|--readme=pocc|--readme_pocc"));
    XHOST.vflag_control.flag("README_APENNSY",aurostd::args2flag(argv,cmds,"--readme=apennsy|--readme_apennsy"));
    XHOST.vflag_control.flag("README_SCRIPTING",aurostd::args2flag(argv,cmds,"--readme=scripting|--readme_scripting|--readme=script|--readme_script"));
    XHOST.vflag_control.flag("README_EXCEPTIONS", aurostd::args2flag(argv, cmds, "--readme=errors|--readme_errors|--readme=exceptions|--readme_exceptions"));  // ME 180531
    XHOST.vflag_control.flag("README_HTRESOURCES",aurostd::args2flag(argv,cmds,"--readme=htresources|--readme_htresources|--readme=resources|--readme_resources"));
    XHOST.vflag_control.flag("README_XAFLOW",aurostd::args2flag(argv,cmds,"--readme=xaflow|--readme_xaflow"));
    XHOST.vflag_control.flag("README_AFLOWRC",aurostd::args2flag(argv,cmds,"--readme=aflowrc|--readme_aflowrc"));
    if(!(XHOST.vflag_control.flag("AFLOW_HELP") || 
                XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3") ||
                XHOST.vflag_control.flag("README_AFLOW") ||
                XHOST.vflag_control.flag("README_AFLOW_PFLOW") ||
                XHOST.vflag_control.flag("README_FROZSL") ||
                XHOST.vflag_control.flag("README_APL") ||
                XHOST.vflag_control.flag("README_QHA") ||
                XHOST.vflag_control.flag("README_AAPL") ||
                XHOST.vflag_control.flag("README_AGL") ||
                XHOST.vflag_control.flag("README_AEL") ||
                XHOST.vflag_control.flag("README_ANRL") ||
                XHOST.vflag_control.flag("README_COMPARE") ||
                XHOST.vflag_control.flag("README_SYMMETRY") ||
                XHOST.vflag_control.flag("README_CHULL") ||
                XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION") ||
                XHOST.vflag_control.flag("README_APENNSY") ||
                XHOST.vflag_control.flag("README_SCRIPTING") ||
                XHOST.vflag_control.flag("README_EXCEPTIONS") ||  // ME180531
                XHOST.vflag_control.flag("README_HTRESOURCES") ||
                XHOST.vflag_control.flag("README_XAFLOW") ||
                XHOST.vflag_control.flag("README_AFLOWRC") ||
                FALSE)){  // CO 180306 - need to catch --readme such that it doesn't interfere with --readme=
                    XHOST.vflag_control.flag("AFLOW_HELP",aurostd::args2flag(argv,cmds,"--readme"));  // CO 180306
                }
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AFLOW_HELP\")=" << XHOST.vflag_control.flag("AFLOW_HELP") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOW_LICENSE_GPL3\")=" << XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOW\")=" << XHOST.vflag_control.flag("README_AFLOW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOW_PFLOW\")=" << XHOST.vflag_control.flag("README_AFLOW_PFLOW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_FROZSL\")=" << XHOST.vflag_control.flag("README_FROZSL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_APL\")=" << XHOST.vflag_control.flag("README_APL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AGL\")=" << XHOST.vflag_control.flag("README_AGL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AEL\")=" << XHOST.vflag_control.flag("README_AEL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_ANRL\")=" << XHOST.vflag_control.flag("README_ANRL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_COMPARE\")=" << XHOST.vflag_control.flag("README_COMPARE") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_SYMMETRY\")=" << XHOST.vflag_control.flag("README_SYMMETRY") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_CHULL\")=" << XHOST.vflag_control.flag("README_CHULL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_PARTIAL_OCCUPATION\")=" << XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_APENNSY\")=" << XHOST.vflag_control.flag("README_APENNSY") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_SCRIPTING\")=" << XHOST.vflag_control.flag("README_SCRIPTING") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_HTRESOURCES\")=" << XHOST.vflag_control.flag("README_HTRESOURCES") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_XAFLOW\")=" << XHOST.vflag_control.flag("README_XAFLOW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOWRC\")=" << XHOST.vflag_control.flag("README_AFLOWRC") << endl;

    // arguments
    string keep=aurostd::args2attachedstring(XHOST.argv,"--keep=","");
    XHOST.vflag_control.flag("KEEP::TEX",aurostd::substring2bool(keep,"tex") || aurostd::substring2bool(keep,"TEX"));
    XHOST.vflag_control.flag("KEEP::DVI",aurostd::substring2bool(keep,"dvi") || aurostd::substring2bool(keep,"DVI"));
    XHOST.vflag_control.flag("KEEP::TOC",aurostd::substring2bool(keep,"toc") || aurostd::substring2bool(keep,"TOC"));
    XHOST.vflag_control.flag("KEEP::EPS",aurostd::substring2bool(keep,"eps") || aurostd::substring2bool(keep,"EPS"));
    XHOST.vflag_control.flag("KEEP::PDF",aurostd::substring2bool(keep,"pdf") || aurostd::substring2bool(keep,"PDF"));
    XHOST.vflag_control.flag("KEEP::JPG",aurostd::substring2bool(keep,"jpg") || aurostd::substring2bool(keep,"JPG"));
    XHOST.vflag_control.flag("KEEP::PNG",aurostd::substring2bool(keep,"png") || aurostd::substring2bool(keep,"PNG"));
    XHOST.vflag_control.flag("KEEP::GIF",aurostd::substring2bool(keep,"gif") || aurostd::substring2bool(keep,"GIF"));
    XHOST.vflag_control.flag("KEEP::GPL",aurostd::substring2bool(keep,"gpl") || aurostd::substring2bool(keep,"GPL") || aurostd::substring2bool(keep,"gnuplot"));
    XHOST.vflag_control.flag("KEEP::MAT",aurostd::substring2bool(keep,"mat") || aurostd::substring2bool(keep,"matlab"));
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::TEX\")=" << XHOST.vflag_control.flag("KEEP::TEX") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::DVI\")=" << XHOST.vflag_control.flag("KEEP::DVI") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::TOC\")=" << XHOST.vflag_control.flag("KEEP::TOC") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::EPS\")=" << XHOST.vflag_control.flag("KEEP::EPS") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::PDF\")=" << XHOST.vflag_control.flag("KEEP::PDF") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::GPL\")=" << XHOST.vflag_control.flag("KEEP::GPL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::MAT\")=" << XHOST.vflag_control.flag("KEEP::MAT") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::JPG\")=" << XHOST.vflag_control.flag("KEEP::JPG") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::PNG\")=" << XHOST.vflag_control.flag("KEEP::PNG") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::GIF\")=" << XHOST.vflag_control.flag("KEEP::GIF") << endl;

    XHOST.vflag_control.flag("PRINT_MODE::HTML",aurostd::args2flag(XHOST.argv,cmds,"--print=html|--print_html")); 
    XHOST.vflag_control.flag("PRINT_MODE::TXT",aurostd::args2flag(XHOST.argv,cmds,"--print=txt|--print_txt")); 
    XHOST.vflag_control.flag("PRINT_MODE::JSON",aurostd::args2flag(XHOST.argv,cmds,"--print=json|--print_json")); // DX 9/7/17 - Add json
    XHOST.vflag_control.flag("PRINT_MODE::LATEX",aurostd::args2flag(XHOST.argv,cmds,"--print=latex|--print_latex"));
    XHOST.vflag_control.flag("PRINT_MODE::YEAR",aurostd::args2flag(XHOST.argv,cmds,"--print=year|--print_year"));
    XHOST.vflag_control.flag("PRINT_MODE::DOI",aurostd::args2flag(XHOST.argv,cmds,"--print=doi|--print_doi"));
    XHOST.vflag_control.flag("PRINT_MODE::EXTRA",aurostd::args2flag(argv,cmds,"--print=extra|--print=vextra_html|--print_vextra_html"));
    XHOST.vflag_control.flag("PRINT_MODE::NUMBER",aurostd::args2flag(XHOST.argv,cmds,"--print=number|--print_number"));
    XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS",aurostd::args2flag(XHOST.argv,cmds,"--print=hyperlinks|--print_hyperlinks|--print=hyperlink|--print_hyperlink"));
    XHOST.vflag_control.flag("PRINT_MODE::NOTE",aurostd::args2flag(XHOST.argv,cmds,"--print=note|--print_note|--print=notes|--print_notes"));
    XHOST.vflag_control.flag("PRINT_MODE::NEW",aurostd::args2flag(XHOST.argv,cmds,"--print=new|--print_new"));
    XHOST.vflag_control.flag("PRINT_MODE::DATE",aurostd::args2flag(XHOST.argv,cmds,"--print=date|--print_date"));
    XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT",aurostd::args2flag(argv,cmds,"--snapshot"));
    XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT",!aurostd::args2flag(argv,cmds,"--NOLATEX|--nolatex"));
    XHOST.vflag_control.flag("APENNSY::LATEX_CITE",aurostd::args2flag(argv,cmds,"--cite"));

    XHOST.vflag_control.flag("OSS::COUT",aurostd::args2flag(XHOST.argv,cmds,"--oss=cout|--oss_cout|--COUT|--cout"));
    XHOST.vflag_control.flag("OSS::CERR",aurostd::args2flag(XHOST.argv,cmds,"--oss=cerr|--oss_cerr|--CERR|--cerr"));
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\")=" << XHOST.vflag_control.flag("PRINT_MODE::HTML") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::TXT\")=" << XHOST.vflag_control.flag("PRINT_MODE::TXT") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::JSON\")=" << XHOST.vflag_control.flag("PRINT_MODE::JSON") << endl; // DX 9/7/17 - Add json
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::LATEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::LATEX") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::YEAR\")=" << XHOST.vflag_control.flag("PRINT_MODE::YEAR") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DOI\")=" << XHOST.vflag_control.flag("PRINT_MODE::DOI") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::EXTRA\")=" << XHOST.vflag_control.flag("PRINT_MODE::EXTRA") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NUMBER\")=" << XHOST.vflag_control.flag("PRINT_MODE::NUMBER") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\")=" << XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NOTE\")=" << XHOST.vflag_control.flag("PRINT_MODE::NOTE") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NEW\")=" << XHOST.vflag_control.flag("PRINT_MODE::NEW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DATE\")=" << XHOST.vflag_control.flag("PRINT_MODE::DATE") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"OSS::COUT\")=" << XHOST.vflag_control.flag("OSS::COUT") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"OSS::CERR\")=" << XHOST.vflag_control.flag("OSS::CERR") << endl;
    //  INIT_VERBOSE=FALSE;

    XHOST.vflag_control.flag("BEEP",aurostd::args2flag(XHOST.argv,cmds,"--beep")); 
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"BEEP\")=" << XHOST.vflag_control.flag("BEEP") << endl;
    XHOST.vflag_control.flag("REMOVE",aurostd::args2flag(XHOST.argv,cmds,"--remove")); 
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"REMOVE\")=" << XHOST.vflag_control.flag("REMOVE") << endl;
    XHOST.vflag_control.flag("ZIP",aurostd::args2flag(XHOST.argv,cmds,"--zip")); 
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"ZIP\")=" << XHOST.vflag_control.flag("ZIP") << endl;

    //  XHOST.vflag_control.flag("PRINT_MODE::EPS",aurostd::args2flag(XHOST.argv,cmds,"--print=eps|--print=eps"));
    XHOST.vflag_control.flag("PRINT_MODE::EPS",TRUE); // default
    XHOST.vflag_control.flag("PRINT_MODE::PDF",aurostd::args2flag(XHOST.argv,cmds,"--print=pdf|--print_pdf")); if(XHOST.vflag_control.flag("PRINT_MODE::PDF")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::GIF",aurostd::args2flag(XHOST.argv,cmds,"--print=gif|--print_gif")); if(XHOST.vflag_control.flag("PRINT_MODE::GIF")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::JPG",aurostd::args2flag(XHOST.argv,cmds,"--print=jpg|--print_jpg")); if(XHOST.vflag_control.flag("PRINT_MODE::JPG")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::PNG",aurostd::args2flag(XHOST.argv,cmds,"--print=png|--print_png")); if(XHOST.vflag_control.flag("PRINT_MODE::PNG")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::EPS\")=" << XHOST.vflag_control.flag("PRINT_MODE::EPS") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::PDF\")=" << XHOST.vflag_control.flag("PRINT_MODE::PDF") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::GIF\")=" << XHOST.vflag_control.flag("PRINT_MODE::GIF") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::JPG\")=" << XHOST.vflag_control.flag("PRINT_MODE::JPG") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::PNG\")=" << XHOST.vflag_control.flag("PRINT_MODE::PNG") << endl;

    XHOST.vflag_control.flag("XPLUG_DO_CLEAN",aurostd::args2flag(XHOST.argv,cmds,"--doclean"));
    XHOST.vflag_control.flag("XPLUG_DO_ADD",aurostd::args2flag(argv,"--add"));

    XHOST.vflag_control.flag("XPLUG_PREFIX",aurostd::args2attachedflag(argv,"--prefix="));
    if(XHOST.vflag_control.flag("XPLUG_PREFIX")) XHOST.vflag_control.push_attached("XPLUG_PREFIX",aurostd::args2attachedstring(argv,"--prefix=",""));
    XHOST.vflag_control.flag("XPLUG_NUM_ZIP",aurostd::args2attachedflag(argv,"--nzip="));
    if(XHOST.vflag_control.flag("XPLUG_NUM_ZIP")) XHOST.vflag_control.push_attached("XPLUG_NUM_ZIP",aurostd::args2attachedstring(argv,"--nzip=",""));
    if(!XHOST.vflag_control.flag("XPLUG_NUM_ZIP")) XHOST.vflag_control.push_attached("XPLUG_NUM_ZIP",aurostd::utype2string(1));
    XHOST.vflag_control.flag("XPLUG_NUM_SIZE",aurostd::args2attachedflag(argv,"--nsize="));
    if(XHOST.vflag_control.flag("XPLUG_NUM_SIZE")) XHOST.vflag_control.push_attached("XPLUG_NUM_SIZE",aurostd::args2attachedstring(argv,"--nsize=",""));
    if(!XHOST.vflag_control.flag("XPLUG_NUM_SIZE")) XHOST.vflag_control.push_attached("XPLUG_NUM_SIZE",aurostd::utype2string(128));
    XHOST.vflag_control.flag("XPLUG_NUM_THREADS",aurostd::args2attachedflag(argv,"--np="));
    XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX",aurostd::args2attachedflag(argv,"--npmax")); // CO 180124
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")) XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS",aurostd::args2attachedstring(argv,"--np=",""));
    if(!XHOST.vflag_control.flag("XPLUG_NUM_THREADS")){
        if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX")){XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS",aurostd::utype2string(XHOST.CPU_Cores));}  // CO 180124
        else{XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS",aurostd::utype2string(XHOST.CPU_Cores/2));}
    }

    // AFLOWLIB_SERVER
    XHOST.vflag_control.flag("AFLOWLIB_SERVER",aurostd::args2attachedflag(argv,"--server="));
    if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) XHOST.vflag_control.push_attached("AFLOWLIB_SERVER",aurostd::args2attachedstring(argv,"--server=",""));
    if(XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="default" || XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="aflowlib") {
        XHOST.vflag_control.pop_attached("AFLOWLIB_SERVER");
        XHOST.vflag_control.push_attached("AFLOWLIB_SERVER","aflowlib.duke.edu");
    }
    if(XHOST.vflag_control.flag("AFLOWLIB_SERVER") &&
            !(XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="aflowlib.duke.edu" || XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="materials.duke.edu")) {
        cerr << "ERROR  init::InitMachine: \"--server=\" can be only \"aflowlib.duke.edu\" or \"materials.duke.edu\"" << endl;
        exit(0);
    }
    // LOAD options
    if(INIT_VERBOSE) oss << "--- LOADING @ options --- " << endl;
    //if(INIT_VERBOSE) XHOST.DEBUG=TRUE;
    if(INIT_VERBOSE) oss << "--- LOADING @ aflow options --- " << endl;
    XHOST.vflag_aflow.flag("LOOP",aurostd::args2flag(XHOST.argv,cmds,"-loop|--loop"));
    XHOST.vflag_aflow.flag("CLEAN",aurostd::args2flag(XHOST.argv,cmds,"-c|-clean|--CLEAN|--clean"));
    XHOST.vflag_aflow.args2addattachedscheme(XHOST.argv,cmds,"XCLEAN","--xclean=","");

    XHOST.AFLOW_RUNDIRflag=aurostd::args2flag(XHOST.argv,"--run|-run");
    XHOST.AFLOW_MULTIflag=aurostd::args2flag(XHOST.argv,"--run=multi|-run=multi|--multi|-multi");	
    XHOST.AFLOW_RUNXflag=!XHOST.AFLOW_MULTIflag && (aurostd::args2attachedflag(XHOST.argv,"--run=") || aurostd::args2attachedflag(XHOST.argv,"-run="));
    XHOST.AFLOW_RUNXnumber=0;
    XHOST.vflag_pflow.clear(); 
    XHOST.vflag_apennsy.clear(); 
    XHOST.vflag_outreach.clear(); 
    if(INIT_VERBOSE) oss << "--- LOADING @ aconvasp options --- " << endl;
    PflowARGs(XHOST.argv,cmds,XHOST.vflag_pflow);
    if(INIT_VERBOSE) oss << "--- LOADING @ apennsy options --- " << endl;
    ApennsyARGs(XHOST.argv,cmds,XHOST.vflag_apennsy); 
    // LOADING OUTREACH
    if(INIT_VERBOSE) oss << "--- LOADING @ outreach options --- " << endl;   
    XHOST.vflag_control.flag("CV::PUBS",aurostd::args2flag(argv,cmds,"--cv=pubs"));
    XHOST.vflag_control.flag("CV::ITALKS",aurostd::args2flag(argv,cmds,"--cv=italks|--cv=talks|--presentations"));
    XHOST.vflag_control.flag("CV::ACADEMIC",aurostd::args2flag(argv,cmds,"--cv=academic"));
    XHOST.vflag_control.flag("CV::RESEARCH",aurostd::args2flag(argv,cmds,"--cv=research"));
    XHOST.vflag_control.flag("CV::EDUCATION",aurostd::args2flag(argv,cmds,"--cv=education"));
    XHOST.vflag_control.flag("CV::TEACHING",aurostd::args2flag(argv,cmds,"--cv=teaching"));
    XHOST.vflag_control.flag("CV::ADVISING",aurostd::args2flag(argv,cmds,"--cv=advising"));
    XHOST.vflag_control.flag("CV::AWARDS",aurostd::args2flag(argv,cmds,"--cv=awards"));
    XHOST.vflag_control.flag("CV::PRESS",aurostd::args2flag(argv,cmds,"--cv=press"));
    XHOST.vflag_control.flag("CV::PATENTS",aurostd::args2flag(argv,cmds,"--cv=patents"));
    XHOST.vflag_control.flag("CV::SERVICE_OUTSIDE",aurostd::args2flag(argv,cmds,"--cv=service_outside|--cv=serviceoutside"));
    XHOST.vflag_control.flag("CV::SERVICE_INSIDE",aurostd::args2flag(argv,cmds,"--cv=service_inside|--cv=serviceinside"));
    XHOST.vflag_control.flag("PHP::CENTER_MISSION",aurostd::args2flag(argv,cmds,"--php_center_mission|--center_mission|--center-mission"));
    XHOST.vflag_control.flag("CV::AUTHOR",aurostd::args2attachedflag(argv,"--author="));
    if(XHOST.vflag_control.flag("CV::AUTHOR")) XHOST.vflag_control.push_attached("CV::AUTHOR",aurostd::args2attachedstring(argv,"--author=",""));
    XHOST.vflag_control.flag("PHP::PUBS_ALLOY",aurostd::args2attachedflag(argv,"--php_pubs_alloy="));
    if(XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) XHOST.vflag_control.push_attached("PHP::PUBS_ALLOY",aurostd::args2attachedstring(argv,"--php_pubs_alloy=",""));
    XHOST.vflag_control.flag("PHP::PUBS_KEYWORD",aurostd::args2attachedflag(argv,"--php_pubs_keyword="));
    if(XHOST.vflag_control.flag("PHP::PUBS_KEYWORD")) XHOST.vflag_control.push_attached("PHP::PUBS_KEYWORD",aurostd::args2attachedstring(argv,"--php_pubs_keyword=",""));
    XHOST.vflag_control.flag("GRANTS",aurostd::args2flag(argv,cmds,"--grant|--grants"));
    if(XHOST.vflag_control.flag("GRANTS")) XHOST.vflag_control.push_attached("GRANTS",aurostd::args2attachedstring(argv,"--grant=|--grants=",""));
    if(LDEBUG) cout << "OUTREACH OPTIONS: xscheme=" << XHOST.vflag_control.xscheme << endl;
    if(LDEBUG) cout << "OUTREACH OPTIONS: vxscheme.size()=" << XHOST.vflag_control.vxscheme.size() << endl;
    if(LDEBUG) cout << "OUTREACH OPTIONS: vxsghost.size()=" << XHOST.vflag_control.vxsghost.size() << endl;
    if(LDEBUG) cout << "OUTREACH OPTIONS: argv.size()=" << argv.size() << endl;

    // LOADING ANRL WEB
    XHOST.vflag_control.flag("WWW",aurostd::args2flag(argv,cmds,"--www|--web|--php|-www|-web|-php"));

    // DEFAULT options
    if(INIT_VERBOSE) oss << "--- DEFAULTSs --- " << endl;
    if(INIT_VERBOSE) aflowrc::print_aflowrc(oss,INIT_VERBOSE || XHOST.DEBUG);
    //if(INIT_VERBOSE) XHOST.DEBUG=TRUE;

    // FINISHED
    if(INIT_VERBOSE) oss << endl;
    if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;
    if(INIT_VERBOSE) oss << "* AFLOW V=" << string(AFLOW_VERSION) << " - machine information " << endl;
    if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;
    if(INIT_VERBOSE) exit(0);

    return TRUE;
}
} // namespace init

// ***************************************************************************
// long init::GetRAM(void)
// ***************************************************************************
namespace init {
#ifndef _MACOSX_
#include <sys/sysinfo.h>
    long GetRAM(void) {
        long pages=sysconf(_SC_PHYS_PAGES);
        long page_size=sysconf(_SC_PAGE_SIZE);
        return pages*page_size;
    }
    long _GetRAM(void) {
        struct sysinfo s;
        if(sysinfo(&s)!=0) {cerr << "sysinfo error" << endl;exit(0);}
        return s.totalram;
    }
#endif
#ifdef _MACOSX_
#include <sys/sysctl.h>
    long GetRAM(void) {
        int mib[2]={CTL_HW,HW_MEMSIZE};
        u_int namelen=sizeof(mib)/sizeof(mib[0]);
        uint64_t size;
        size_t len=sizeof(size);
        if(sysctl(mib,namelen,&size,&len,NULL,0)<0) {cerr << "ERROR sysctl in init::GetRAM" << endl;exit(0);}
        return (long) size;
    }
#endif
} // namespace init

// ***************************************************************************
// init::InitLoadString
// ***************************************************************************
namespace init {
    string InitLoadString(string str2load,bool LVERBOSE) {
        bool LDEBUG=(FALSE || XHOST.DEBUG || LVERBOSE);
        string soliloquy="init::InitLoadString:";
        if(LDEBUG) cerr << soliloquy << " str2load=" << str2load << endl; 
        if(!XHOST.is_command("aflow_data")) {
            cerr << "AFLOW Error: " << "aflow_data" << " is not in the path... exiting.." << endl;
            exit(0);
        } 
        if(LDEBUG) cerr << "00000  MESSAGE AFLOW INIT Loading data = [" << str2load << "]";
        if(LDEBUG) cerr.flush();
        string out;
        string aflow_data_path=aurostd::args2attachedstring(XHOST.argv,"--aflow_data_path=",(string) "");
        if(aflow_data_path=="") {
            if(XHOST.hostname=="nietzsche.mems.duke.edu"&&XHOST.User=="auro"&&aurostd::FileExist(XHOST.Home+"/work/AFLOW3/aflow_data")) {  // CO, special stefano
                out=aurostd::execute2string(string(XHOST.Home+"/work/AFLOW3/aflow_data")+string(" ")+str2load);
                if(LDEBUG) cerr << soliloquy << " FOUND " << XHOST.Home << "/work/AFLOW3/aflow_data" << endl;
                if(LDEBUG) cerr << soliloquy << " out=" << out << endl; 
                if(LDEBUG) cerr << soliloquy << " str2load=" << str2load << endl; 
            } else {
                if(LDEBUG){cerr << soliloquy << " issuing command: " << XHOST.command("aflow_data") << " " << str2load << endl;}
                out=aurostd::execute2string(XHOST.command("aflow_data")+" "+str2load);
            }
        } else { // cerr << string(aflow_data_path+"/"+XHOST.command("aflow_data")) << endl;
            out=aurostd::execute2string(aflow_data_path+"/"+XHOST.command("aflow_data")+" "+str2load);
        }

        if(LDEBUG) cerr << soliloquy << "  length=" << out.length() << endl;
        if(LDEBUG) cerr.flush();
        //    if(LDEBUG) exit(0);
        return out;
    }
} // namespace init

vector<string> vvAURL,vvAUID,vvLOOP;

string vAURL_cutout(string cutout,bool LVERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || LVERBOSE);
    stringstream sss;
    if(LDEBUG) cerr << "vAURL_cutout: XHOST_AURL.size()=" << XHOST_AURL.size() << "  " << "vvAURL.size()=" << vvAURL.size() << endl;
    vector<string> tokens;
    aurostd::string2tokens(cutout,tokens,"|");
    if(!vvAURL.size()) aurostd::string2vectorstring(init::InitGlobalObject("vAURL","",LVERBOSE),vvAURL); // a bit of recursivity helps
    for(uint i=0;i<tokens.size();i++) {
        for(uint j=0;j<vvAURL.size();j++) {
            if(aurostd::substring2bool(vvAURL.at(j),tokens.at(i))) {
                aurostd::StringSubst(vvAURL.at(j),tokens.at(i),"");
                aurostd::StringSubst(vvAURL.at(j),"aflowlib.duke.edu:","");
                aurostd::StringSubst(vvAURL.at(j),"materials.duke.edu:","");
                sss << vvAURL.at(j) << endl;
            }
        }
    }
    return sss.str();
}

// ***************************************************************************
// init::InitGlobalObject
// ***************************************************************************
namespace init {
    string InitGlobalObject(string str,string grep,bool LVERBOSE) {
        bool LDEBUG=(FALSE || XHOST.DEBUG);
        // || LVERBOSE;
        string out="";
        string FileLibrary="";

        if(str=="Library_HTQC") { 
            if(XHOST_Library_HTQC.empty()) {
                return XHOST_Library_HTQC=init::InitLoadString(str,LVERBOSE);
            } else { 
                return XHOST_Library_HTQC;
            }
        } // FIX
        // FILES CALCULATED
        // README_CALCULATED_ICSD
        if(str=="Library_CALCULATED_ICSD_LIB") { if(XHOST_Library_CALCULATED_ICSD_LIB.empty()) { return XHOST_Library_CALCULATED_ICSD_LIB=vAURL_cutout("AFLOWDATA/ICSD_RAW/|AFLOWDATA/ICSD_WEB/",LVERBOSE);} else { return XHOST_Library_CALCULATED_ICSD_LIB;}} 
        if(str=="Library_CALCULATED_ICSD_RAW") { if(XHOST_Library_CALCULATED_ICSD_RAW.empty()) { return XHOST_Library_CALCULATED_ICSD_RAW=vAURL_cutout("AFLOWDATA/ICSD_RAW/|AFLOWDATA/ICSD_WEB/",LVERBOSE);} else { return XHOST_Library_CALCULATED_ICSD_RAW;}}  
        // README_CALCULATED_LIB1
        if(str=="Library_CALCULATED_LIB1_LIB") { if(XHOST_Library_CALCULATED_LIB1_LIB.empty()) { return XHOST_Library_CALCULATED_LIB1_LIB=vAURL_cutout("AFLOWDATA/LIB1_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB1_LIB;}} 
        if(str=="Library_CALCULATED_LIB1_RAW") { if(XHOST_Library_CALCULATED_LIB1_RAW.empty()) { return XHOST_Library_CALCULATED_LIB1_RAW=vAURL_cutout("AFLOWDATA/LIB1_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB1_RAW;}} 
        // README_CALCULATED_LIB2
        if(str=="Library_CALCULATED_LIB2_LIB") { if(XHOST_Library_CALCULATED_LIB2_LIB.empty()) { return XHOST_Library_CALCULATED_LIB2_LIB=vAURL_cutout("AFLOWDATA/LIB2_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB2_LIB;}} 
        if(str=="Library_CALCULATED_LIB2_RAW") { if(XHOST_Library_CALCULATED_LIB2_RAW.empty()) { return XHOST_Library_CALCULATED_LIB2_RAW=vAURL_cutout("AFLOWDATA/LIB2_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB2_RAW;}} 
        // README_CALCULATED_LIB3
        if(str=="Library_CALCULATED_LIB3_LIB") { if(XHOST_Library_CALCULATED_LIB3_LIB.empty()) { return XHOST_Library_CALCULATED_LIB3_LIB=vAURL_cutout("AFLOWDATA/LIB3_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB3_LIB;}} 
        if(str=="Library_CALCULATED_LIB3_RAW") { if(XHOST_Library_CALCULATED_LIB3_RAW.empty()) { return XHOST_Library_CALCULATED_LIB3_RAW=vAURL_cutout("AFLOWDATA/LIB3_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB3_RAW;}} 
        // README_CALCULATED_LIB4
        if(str=="Library_CALCULATED_LIB4_LIB") { if(XHOST_Library_CALCULATED_LIB4_LIB.empty()) { return XHOST_Library_CALCULATED_LIB4_LIB=vAURL_cutout("AFLOWDATA/LIB4_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB4_LIB;}} 
        if(str=="Library_CALCULATED_LIB4_RAW") { if(XHOST_Library_CALCULATED_LIB4_RAW.empty()) { return XHOST_Library_CALCULATED_LIB4_RAW=vAURL_cutout("AFLOWDATA/LIB4_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB4_RAW;}} 
        // README_CALCULATED_LIB5
        if(str=="Library_CALCULATED_LIB5_LIB") { if(XHOST_Library_CALCULATED_LIB5_LIB.empty()) { return XHOST_Library_CALCULATED_LIB5_LIB=vAURL_cutout("AFLOWDATA/LIB5_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB5_LIB;}} 
        if(str=="Library_CALCULATED_LIB5_RAW") { if(XHOST_Library_CALCULATED_LIB5_RAW.empty()) { return XHOST_Library_CALCULATED_LIB5_RAW=vAURL_cutout("AFLOWDATA/LIB5_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB5_RAW;}} 
        // README_CALCULATED_LIB6
        if(str=="Library_CALCULATED_LIB6_LIB") { if(XHOST_Library_CALCULATED_LIB6_LIB.empty()) { return XHOST_Library_CALCULATED_LIB6_LIB=vAURL_cutout("AFLOWDATA/LIB6_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB6_LIB;}} 
        if(str=="Library_CALCULATED_LIB6_RAW") { if(XHOST_Library_CALCULATED_LIB6_RAW.empty()) { return XHOST_Library_CALCULATED_LIB6_RAW=vAURL_cutout("AFLOWDATA/LIB6_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB6_RAW;}} 
        // README_CALCULATED_LIB7
        if(str=="Library_CALCULATED_LIB7_LIB") { if(XHOST_Library_CALCULATED_LIB7_LIB.empty()) { return XHOST_Library_CALCULATED_LIB7_LIB=vAURL_cutout("AFLOWDATA/LIB7_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB7_LIB;}} 
        if(str=="Library_CALCULATED_LIB7_RAW") { if(XHOST_Library_CALCULATED_LIB7_RAW.empty()) { return XHOST_Library_CALCULATED_LIB7_RAW=vAURL_cutout("AFLOWDATA/LIB7_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB7_RAW;}} 
        // README_CALCULATED_LIB8
        if(str=="Library_CALCULATED_LIB8_LIB") { if(XHOST_Library_CALCULATED_LIB8_LIB.empty()) { return XHOST_Library_CALCULATED_LIB8_LIB=vAURL_cutout("AFLOWDATA/LIB8_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB8_LIB;}} 
        if(str=="Library_CALCULATED_LIB8_RAW") { if(XHOST_Library_CALCULATED_LIB8_RAW.empty()) { return XHOST_Library_CALCULATED_LIB8_RAW=vAURL_cutout("AFLOWDATA/LIB8_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB8_RAW;}} 
        // README_CALCULATED_LIB9
        if(str=="Library_CALCULATED_LIB9_LIB") { if(XHOST_Library_CALCULATED_LIB9_LIB.empty()) { return XHOST_Library_CALCULATED_LIB9_LIB=vAURL_cutout("AFLOWDATA/LIB9_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB9_LIB;}} 
        if(str=="Library_CALCULATED_LIB9_RAW") { if(XHOST_Library_CALCULATED_LIB9_RAW.empty()) { return XHOST_Library_CALCULATED_LIB9_RAW=vAURL_cutout("AFLOWDATA/LIB9_RAW/",LVERBOSE);} else { return XHOST_Library_CALCULATED_LIB9_RAW;}} 
        // AUID AURL LOOP
        if(str=="vAUID") { if(XHOST_AUID.empty()) {XHOST_AUID=aurostd::RemoveWhiteSpaces(init::InitLoadString(str,LVERBOSE));aurostd::string2vectorstring(XHOST_AUID,XHOST_vAUID);return XHOST_AUID;} else { return XHOST_AUID;}} // 
        if(str=="vAURL") { if(XHOST_AURL.empty()) {XHOST_AURL=aurostd::RemoveWhiteSpaces(init::InitLoadString(str,LVERBOSE));aurostd::string2vectorstring(XHOST_AURL,XHOST_vAURL);return XHOST_AURL;} else { return XHOST_AURL;}} // 
        if(str=="vLOOP") { if(XHOST_LOOP.empty()) {XHOST_LOOP=aurostd::RemoveWhiteSpaces(init::InitLoadString(str,LVERBOSE));aurostd::string2vectorstring(XHOST_LOOP,XHOST_vLOOP);return XHOST_LOOP;} else { return XHOST_LOOP;}} // 
        // [OBSOLETE]   // TABLE PROTOTYPES
        // [OBSOLETE]   if(str=="AFLOW_BinaryRead") { if(XHOST_AFLOW_BinaryRead.empty()) { return XHOST_AFLOW_BinaryRead=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_BinaryRead;}} //
        // [OBSOLETE]   if(str=="AFLOW_Binary_Angle_Read") { if(XHOST_AFLOW_Binary_Angle_Read.empty()) { return XHOST_AFLOW_Binary_Angle_Read=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_Binary_Angle_Read;}} // 
        // AFLOWLIB THINGS
        // LOAD
        // if(str=="aflowlib_lib1") { if(XHOST_aflowlib_lib1.empty()) { return XHOST_aflowlib_lib1=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib1;}} //        // if(str=="aflowlib_lib2") { if(XHOST_aflowlib_lib2.empty()) { return XHOST_aflowlib_lib2=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib2;}} //
        // if(str=="aflowlib_lib3") { if(XHOST_aflowlib_lib3.empty()) { return XHOST_aflowlib_lib3=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib3;}} // 
        // if(str=="aflowlib_lib4") { if(XHOST_aflowlib_lib4.empty()) { return XHOST_aflowlib_lib4=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib4;}} // 
        // if(str=="aflowlib_lib5") { if(XHOST_aflowlib_lib5.empty()) { return XHOST_aflowlib_lib5=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib5;}} // 
        // if(str=="aflowlib_lib6") { if(XHOST_aflowlib_lib6.empty()) { return XHOST_aflowlib_lib6=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib6;}} //  
        // if(str=="aflowlib_lib7") { if(XHOST_aflowlib_lib7.empty()) { return XHOST_aflowlib_lib7=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib7;}} //  
        // if(str=="aflowlib_lib8") { if(XHOST_aflowlib_lib8.empty()) { return XHOST_aflowlib_lib8=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib8;}} //  
        // if(str=="aflowlib_lib9") { if(XHOST_aflowlib_lib9.empty()) { return XHOST_aflowlib_lib9=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib9;}} //  
        // if(str=="aflowlib_icsd") { if(XHOST_aflowlib_icsd.empty()) { return XHOST_aflowlib_icsd=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_icsd;}} // 
        // LOAD
        // STOKES THINGS
        if(str=="FINDSYM_data_space_txt") { if(XHOST_FINDSYM_data_space_txt.empty()) { return XHOST_FINDSYM_data_space_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FINDSYM_data_space_txt;}} // LOADED TXTS
        if(str=="FINDSYM_data_wyckoff_txt") { if(XHOST_FINDSYM_data_wyckoff_txt.empty()) { return XHOST_FINDSYM_data_wyckoff_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FINDSYM_data_wyckoff_txt;}} // LOADED TXTS
        if(str=="FROZSL_data_space_txt") { if(XHOST_FROZSL_data_space_txt.empty()) { return XHOST_FROZSL_data_space_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_space_txt;}} // LOADED TXTS
        if(str=="FROZSL_data_wyckoff_txt") { if(XHOST_FROZSL_data_wyckoff_txt.empty()) { return XHOST_FROZSL_data_wyckoff_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_wyckoff_txt;}} // LOADED TXTS
        if(str=="FROZSL_data_images_txt") { if(XHOST_FROZSL_data_images_txt.empty()) { return XHOST_FROZSL_data_images_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_images_txt;}} // LOADED TXTS
        if(str=="FROZSL_data_irreps_txt") { if(XHOST_FROZSL_data_irreps_txt.empty()) { return XHOST_FROZSL_data_irreps_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_irreps_txt;}} // LOADED TXTS
        if(str=="FROZSL_data_isotropy_txt") { if(XHOST_FROZSL_data_isotropy_txt.empty()) { return XHOST_FROZSL_data_isotropy_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_isotropy_txt;}} // LOADED TXTS
        if(str=="FROZSL_data_little_txt") { if(XHOST_FROZSL_data_little_txt.empty()) { return XHOST_FROZSL_data_little_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_little_txt;}} // LOADED TXTS
        if(str=="FROZSL_symmetry2_dat") { if(XHOST_FROZSL_symmetry2_dat.empty()) { return XHOST_FROZSL_symmetry2_dat=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_symmetry2_dat;}} // LOADED TXTS
        if(str=="FROZSL_const_dat") { if(XHOST_FROZSL_const_dat.empty()) { return XHOST_FROZSL_const_dat=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_const_dat;}} // LOADED TXTS
        if(str=="FROZSL_phvaspsetup_AFLOW") { if(XHOST_FROZSL_phvaspsetup_AFLOW.empty()) { return XHOST_FROZSL_phvaspsetup_AFLOW=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_phvaspsetup_AFLOW;}} // LOADED TXTS
        if(str=="FROZSL_phvaspsetup_POSCAR") { if(XHOST_FROZSL_phvaspsetup_POSCAR.empty()) { return XHOST_FROZSL_phvaspsetup_POSCAR=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_phvaspsetup_POSCAR;}} // LOADED TXTS
        // README THINGS
        if(str=="README_AFLOW_LICENSE_GPL3_TXT") { if(XHOST_README_AFLOW_LICENSE_GPL3_TXT.empty()) { return XHOST_README_AFLOW_LICENSE_GPL3_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_LICENSE_GPL3_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_TXT") { if(XHOST_README_AFLOW_TXT.empty()) { return XHOST_README_AFLOW_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_PFLOW_TXT") { if(XHOST_README_AFLOW_PFLOW_TXT.empty()) { return XHOST_README_AFLOW_PFLOW_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_PFLOW_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_APENNSY_TXT") { if(XHOST_README_AFLOW_APENNSY_TXT.empty()) { return XHOST_README_AFLOW_APENNSY_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_APENNSY_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_SCRIPTING_TXT") { if(XHOST_README_AFLOW_SCRIPTING_TXT.empty()) { return XHOST_README_AFLOW_SCRIPTING_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_SCRIPTING_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_FROZSL_TXT") { if(XHOST_README_AFLOW_FROZSL_TXT.empty()) { return XHOST_README_AFLOW_FROZSL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_FROZSL_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_POCC_TXT") { if(XHOST_README_AFLOW_POCC_TXT.empty()) { return XHOST_README_AFLOW_POCC_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_POCC_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_APL_TXT") { if(XHOST_README_AFLOW_APL_TXT.empty()) { return XHOST_README_AFLOW_APL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_APL_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_QHA_SCQHA_QHA3P_TXT") { if(XHOST_README_AFLOW_QHA_SCQHA_QHA3P_TXT.empty()) { return XHOST_README_AFLOW_QHA_SCQHA_QHA3P_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_QHA_SCQHA_QHA3P_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_AGL_TXT") { if(XHOST_README_AFLOW_AGL_TXT.empty()) { return XHOST_README_AFLOW_AGL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_AGL_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_AEL_TXT") { if(XHOST_README_AFLOW_AEL_TXT.empty()) { return XHOST_README_AFLOW_AEL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_AEL_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_ANRL_TXT") { if(XHOST_README_AFLOW_ANRL_TXT.empty()) { return XHOST_README_AFLOW_ANRL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_ANRL_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_COMPARE_TXT") { if(XHOST_README_AFLOW_COMPARE_TXT.empty()) { return XHOST_README_AFLOW_COMPARE_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_COMPARE_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_SYM_TXT") { if(XHOST_README_AFLOW_SYM_TXT.empty()) { return XHOST_README_AFLOW_SYM_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_SYM_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_CHULL_TXT") { if(XHOST_README_AFLOW_CHULL_TXT.empty()) { return XHOST_README_AFLOW_CHULL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_CHULL_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_EXCEPTIONS_TXT") {if(XHOST_README_AFLOW_EXCEPTIONS_TXT.empty()){ return XHOST_README_AFLOW_EXCEPTIONS_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_EXCEPTIONS_TXT;}}  // ME180531
        if(str=="README_AFLOW_HTRESOURCES_TXT") { if(XHOST_README_AFLOW_HTRESOURCES_TXT.empty()) { return XHOST_README_AFLOW_HTRESOURCES_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_HTRESOURCES_TXT;}} // LOADED TXTS
        if(str=="README_PROTO_TXT") { if(XHOST_README_PROTO_TXT.empty()) { return XHOST_README_PROTO_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_PROTO_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_XAFLOW_TXT") { if(XHOST_README_AFLOW_XAFLOW_TXT.empty()) { return XHOST_README_AFLOW_XAFLOW_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_XAFLOW_TXT;}} // LOADED TXTS
        if(str=="README_AFLOW_AFLOWRC_TXT") { if(XHOST_README_AFLOW_AFLOWRC_TXT.empty()) { return XHOST_README_AFLOW_AFLOWRC_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_AFLOWRC_TXT;}} // LOADED TXTS
        // SCINTILLATION THINGS
        if(str=="ElectronStoppingPower_txt") { if(XHOST_ElectronStoppingPower_txt.empty()) { return XHOST_ElectronStoppingPower_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_ElectronStoppingPower_txt;}} // LOADED TXTS
        if(str=="PhotonCrossSection_txt") { if(XHOST_PhotonCrossSection_txt.empty()) { return XHOST_PhotonCrossSection_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_PhotonCrossSection_txt;}} // LOADED TXTS
        if(str=="PhotonStoppingPower_txt") { if(XHOST_PhotonStoppingPower_txt.empty()) { return XHOST_PhotonStoppingPower_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_PhotonStoppingPower_txt;}} // LOADED TXTS
        if(str=="ICSD_List_txt") { if(XHOST_ICSD_List_txt.empty()) { return XHOST_ICSD_List_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_ICSD_List_txt;}} // LOADED TXTS
        if(str=="AFLOW_PSEUDOPOTENTIALS") { if(XHOST_AFLOW_PSEUDOPOTENTIALS.empty()) { return XHOST_AFLOW_PSEUDOPOTENTIALS=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_PSEUDOPOTENTIALS;}} // LOADED TXTS
        if(str=="AFLOW_PSEUDOPOTENTIALS_TXT") { if(XHOST_AFLOW_PSEUDOPOTENTIALS_TXT.empty()) { return XHOST_AFLOW_PSEUDOPOTENTIALS_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_PSEUDOPOTENTIALS_TXT;}} // LOADED TXTS
        if(str=="AFLOW_PSEUDOPOTENTIALS_LIST_TXT") { if(XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT.empty()) { return XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT;}} // LOADED TXTS
        if(str=="f144468a7ccc2d3a72ba44000715efdb") { if(XHOST_f144468a7ccc2d3a72ba44000715efdb.empty()) { return XHOST_f144468a7ccc2d3a72ba44000715efdb=init::InitLoadString(str,LVERBOSE);} else { return XHOST_f144468a7ccc2d3a72ba44000715efdb;}} // LOADED TXTS
        if(str=="d0f1b0e47f178ae627a388d3bf65d2d2") { if(XHOST_d0f1b0e47f178ae627a388d3bf65d2d2.empty()) { return XHOST_d0f1b0e47f178ae627a388d3bf65d2d2=init::InitLoadString(str,LVERBOSE);} else { return XHOST_d0f1b0e47f178ae627a388d3bf65d2d2;}} // LOADED TXTS
        if(str=="decf00ca3ad2fe494eea8e543e929068") { if(XHOST_decf00ca3ad2fe494eea8e543e929068.empty()) { return XHOST_decf00ca3ad2fe494eea8e543e929068=init::InitLoadString(str,LVERBOSE);} else { return XHOST_decf00ca3ad2fe494eea8e543e929068;}} // LOADED TXTS

        // SEARCH IN AFLOW_DATA AND IF NOT FROM LIBRARIES
        // pure and auro are inside aflow_data
        if(str=="Library_ICSD" || str=="aflowlib_lib1" || str=="aflowlib_lib2" || str=="aflowlib_lib3" || str=="aflowlib_lib4" || str=="aflowlib_lib5" || str=="aflowlib_lib6" || str=="aflowlib_icsd") {
            //  cerr << "(*vLibrary).length()=" << (*vLibrary).length() << endl;
            // if(LDEBUG)

            string *vLibrary; 
            if(str=="Library_ICSD")  vLibrary=&XHOST_Library_ICSD_ALL;
            if(str=="aflowlib_icsd") vLibrary=&XHOST_aflowlib_icsd;
            if(str=="aflowlib_lib1") vLibrary=&XHOST_aflowlib_lib1;
            if(str=="aflowlib_lib2") vLibrary=&XHOST_aflowlib_lib2;
            if(str=="aflowlib_lib3") vLibrary=&XHOST_aflowlib_lib3;
            if(str=="aflowlib_lib4") vLibrary=&XHOST_aflowlib_lib4;
            if(str=="aflowlib_lib5") vLibrary=&XHOST_aflowlib_lib5;
            if(str=="aflowlib_lib6") vLibrary=&XHOST_aflowlib_lib6;
            // check SHORTCUTS
            cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitGlobalObject: Loading library \"" << str << "\" from \"ram\"" << endl;
            if(!(*vLibrary).empty()) {return (*vLibrary);}
            // THEN try aflow_data
            cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitGlobalObject: Loading library \"" << str << "\" from \"aflow_data\"" << endl;
            (*vLibrary)=init::InitLoadString(str,LVERBOSE);
            if(!(*vLibrary).empty()) {return (*vLibrary);}
            // THEN try aflow libraries
            cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitGlobalObject: Loading library \"" << str << "\" from \"aflow libs\"" << endl;

            // check if available
            if((*vLibrary).empty()) {   // find and LOAD
                string str2search=str;
                aurostd::StringSubst(str2search,"Library_ICSD","aflow_library_icsd");
                (*vLibrary)="";
                for(uint j=0;j<vAFLOW_LIBRARY_DIRECTORIES.size() && (*vLibrary).empty();j++) {   // cycle through possible directories
                    FileLibrary=aurostd::CleanFileName(vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/"+str2search+".dat");
                    if(aurostd::FileExist(FileLibrary) && aurostd::FileNotEmpty(FileLibrary)) {
                        if(LDEBUG || LVERBOSE) cerr << "00000  AFLOW LIBRARY  (" << j << ")  found=" <<  FileLibrary << endl;
                        if(LDEBUG || LVERBOSE) cerr << "loading... ";
                        if(LDEBUG || LVERBOSE) cerr.flush();
                        if(grep=="") {
                            aurostd::file2string(FileLibrary,(*vLibrary));
                        } else {
                            (*vLibrary)=aurostd::execute2string("cat "+FileLibrary+" | grep -E '"+grep+"'");
                        }
                        if(LDEBUG || LVERBOSE) cerr << "length=" << (*vLibrary).size();// << " " << endl;
                        if(LDEBUG || LVERBOSE) cerr.flush();
                    }
                } // cycle through possible directories
                if((*vLibrary).empty()) {
                    cerr << "WARNING - init::InitGlobalObject: " << str << " not found! " << endl;// exit(0);
                    return "";
                }
                out=(*vLibrary);
            }
        } 

        if(out=="") {
            //    cerr << "ERROR: init::InitGlobalObject str = " << str << " not found ..." << endl; // exit(0);
        }
        return out;
    }
} // namespace init

// ***************************************************************************
// init::InitLibraryObject
// ***************************************************************************
namespace init {
    string InitLibraryObject(string str,bool LVERBOSE) {
        bool LDEBUG=FALSE;
        // Search LIBRARY
        string str2search=str;
        aurostd::StringSubst(str2search,"Library_ICSD","aflow_library_icsd.dat");
        if(str=="Library_ICSD") { 
            if(XHOST_Library_ICSD_ALL.empty()) { 
                return XHOST_Library_ICSD_ALL=init::InitLoadString(str,LVERBOSE);
            } else { 
                return XHOST_Library_ICSD_ALL;
            }
        } // FIX

        string FileLibrary,out="";
        for(uint j=0;j<vAFLOW_LIBRARY_DIRECTORIES.size() && out.empty();j++) {   // cycle through possible directories
            FileLibrary=vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/"+str2search;
            if(LDEBUG || LVERBOSE) cerr << "DDDDD  InitLibraryObject: (" << j << ")  " <<  FileLibrary << endl;
            if(aurostd::FileExist(FileLibrary) && aurostd::FileNotEmpty(FileLibrary))
                out=aurostd::file2string(FileLibrary,out);
        } // cycle through possible directories
        if(!out.empty()) {
            if(LDEBUG || LVERBOSE) cerr << "00000  MESSAGE InitLibraryObject: AFLOW LIBRARY  Found library file = [" << FileLibrary << "]" << endl;
        } else {
            cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") initLibraryObject: AFLOW_LIBRARY not found! " << endl;
            //     exit(0);
        }

        return out;
    }
} // namespace init


// ***************************************************************************
// uint init::GetTEMPs(void) // need sensors package
// ***************************************************************************
namespace init {
    pthread_mutex_t mutex_INIT_GetTEMPs=PTHREAD_MUTEX_INITIALIZER;
    uint GetTEMPs(void) {
        pthread_mutex_lock(&mutex_INIT_GetTEMPs);
        // pthread_mutex_unlock(&mutex_INIT_GetTEMPs);

        // if(aurostd::execute2string("ps aux | grep sensors | grep -v sensorsd | grep -v grep")!="") LOCAL_is_sensor=FALSE; // must postpone
        XHOST.vTemperatureCore.clear();
        if(XHOST.sensors_allowed) {
            bool LOCAL_is_sensor=XHOST.is_command("sensors");
            // test sensors
            if(LOCAL_is_sensor)
                if(aurostd::execute2utype<int>("bash -c \"sensors 2>&1 2> /dev/null\" | grep -c temp")==0)
                    LOCAL_is_sensor=FALSE;
            // now check
            if(LOCAL_is_sensor) { // need sensors package
                vector<string> vline_temp1,vline_temp2,tokens,tokens2; string sj=" ";
                aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("sensors")+" | grep temp1 | grep -v low | grep high | head -8 "),vline_temp1);
                aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("sensors")+" | grep temp2 | grep -v low | grep high | head -8 "),vline_temp2);
                for(uint i=0;i<vline_temp1.size();i++) {
                    aurostd::StringSubst(vline_temp1.at(i),"temp1"," ");
                    aurostd::StringSubst(vline_temp1.at(i),":"," ");aurostd::StringSubst(vline_temp1.at(i),"="," ");
                    aurostd::StringSubst(vline_temp1.at(i),"("," ");aurostd::StringSubst(vline_temp1.at(i),")"," ");
                    //aurostd::StringSubst(vline_temp1.at(i),""," ");aurostd::StringSubst(vline_temp1.at(i),"C"," "); // CO, PREVIOUSLY DEGREE SYMBOL, removed to get rid of warnings on mac
                    aurostd::StringSubst(vline_temp1.at(i),"\u00B0"," ");aurostd::StringSubst(vline_temp1.at(i),"C"," "); //\u00B0:  http://www.fileformat.info/info/unicode/char/b0/index.htm
                    aurostd::StringSubst(vline_temp1.at(i),"+"," ");aurostd::StringSubst(vline_temp1.at(i),"-"," ");
                    for(unsigned char j=1;j<255;j++)
                        if((j>=1 && j<=45 && j!=32) || (j==47) || (j>=58 && j<255)) {sj[0]=j;aurostd::StringSubst(vline_temp1.at(i),sj," ");}
                    //   cerr << vline_temp1.at(i) << endl;
                    aurostd::string2tokens(vline_temp1.at(i),tokens," ");
                    // cerr << vline_temp1.at(i) << endl;
                    for(uint j=0;j<tokens.size();j++) {
                        if(aurostd::substring2bool(tokens.at(j),".")) {
                            // aurostd::string2tokens(tokens.at(j),tokens2,".");
                            tokens2=tokens;
                            if(tokens2.size()>0) {XHOST.vTemperatureCore.push_back(aurostd::string2utype<double>(tokens2.at(0)));}
                            break;
                        }
                    }
                    // check for combinations
                }
                //  while(vline_temp2.size()>0) {vline_temp1.pop_back();vline_temp2.pop_back();} // remove the temp2&temp1 stuff
                if(vline_temp1.size()==9) vline_temp1.pop_back();
            } else {
                XHOST.vTemperatureCore.clear();
            }
        }
        pthread_mutex_unlock(&mutex_INIT_GetTEMPs);
        return XHOST.vTemperatureCore.size();
    }
} // namespace init

// ***************************************************************************
// AFLOW_getTEMP
// ***************************************************************************
uint AFLOW_getTEMP(vector<string> argv) {
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    bool RUNBAR=aurostd::args2flag(argv,"--runbar|--RUNBAR|--runBAR|--bar|--BAR");
    bool RUNSTAT=aurostd::args2flag(argv,"--runstat|--RUNSTAT|--runSTAT|--stat|--STAT");
    string WRITE=aurostd::args2attachedstring(argv,"--write=","");
    double maxmem=aurostd::args2attachedutype<double>(argv,"--mem=|--maxmem=",XHOST.maxmem);
    double refresh=aurostd::args2attachedutype<double>(argv,"--refresh=",AFLOW_CORE_TEMPERATURE_REFRESH);
    double warning_beep=aurostd::args2attachedutype<double>(argv,"--warning_beep=",AFLOW_CORE_TEMPERATURE_BEEP);
    double warning_halt=aurostd::args2attachedutype<double>(argv,"--warning_halt=",AFLOW_CORE_TEMPERATURE_HALT);

    if(LDEBUG) cerr << "AFLOW_getTEMP: RUNBAR=" << RUNBAR << endl;
    if(LDEBUG) cerr << "AFLOW_getTEMP: RUNSTAT=" << RUNSTAT << endl;
    if(LDEBUG) cerr << "AFLOW_getTEMP: write=" << WRITE << endl;
    if(LDEBUG) cerr << "AFLOW_getTEMP: maxmem=" << maxmem << endl;
    if(LDEBUG) cerr << "AFLOW_getTEMP: refresh=" << refresh << "   -   AFLOW_CORE_TEMPERATURE_REFRESH=" << AFLOW_CORE_TEMPERATURE_REFRESH << endl;
    if(LDEBUG) cerr << "AFLOW_getTEMP: warning_beep=" << warning_beep << "   -   AFLOW_CORE_TEMPERATURE_BEEP=" << AFLOW_CORE_TEMPERATURE_BEEP << endl;
    if(LDEBUG) cerr << "AFLOW_getTEMP: warning_halt=" << warning_halt << "   -   AFLOW_CORE_TEMPERATURE_HALT=" << AFLOW_CORE_TEMPERATURE_HALT << endl;
    if(WRITE!="") aurostd::RemoveFile(WRITE);

    while(init::GetTEMPs()) {
        stringstream oss;
        double Tmax=aurostd::max(XHOST.vTemperatureCore);
        double Tmin=aurostd::min(XHOST.vTemperatureCore);
        double Tzero=30.0;
        oss << "00000  MESSAGE " << aurostd::get_time() << " ";// << Message("host") << endl; exit(1);
        if(RUNSTAT || (!RUNSTAT && !RUNBAR)) {
            string soss="- [temp(C)=";
            for(uint i=0;i<XHOST.vTemperatureCore.size();i++) {soss+=aurostd::utype2string(XHOST.vTemperatureCore.at(i),3)+(i<XHOST.vTemperatureCore.size()-1?",":"]");}
            oss << aurostd::PaddedPOST(soss,XHOST.vTemperatureCore.size()*5+11," ");
            soss=" - ["+aurostd::utype2string(Tmin,3)+","+aurostd::utype2string(Tmax,3)+"]";
            oss << aurostd::PaddedPOST(soss,14," ") << "  beep=" << warning_beep << "  halt=" << warning_halt << " ";
        }
        if(RUNBAR) {
            for(double i=0;i<2*(Tmax-30.0);i++) oss << "*";
            oss << " " << Tmax;
        }
        if(maxmem<100.0) oss << " - [mem=" << aurostd::utype2string<double>(AFLOW_checkMEMORY("vasp",maxmem),4) << " (" << maxmem << ")]";
        if(Tmax>=warning_beep) oss << "  (beep) MAX>" << warning_beep;// << endl;
        if(Tmax>=warning_halt) oss << "  (halt) SHUTDOWN>" << warning_halt;// << endl;

        if(WRITE!="") {
            stringstream aus;vector<string> vlines;
            //      if(aurostd::FileExist(WRITE))
            {aurostd::file2vectorstring(WRITE,vlines);}
            aus<<"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">"<<endl;
            aus<<"<html> <head>"<<endl;
            aus<<"<META HTTP-EQUIV=\"expires\" CONTENT=\"0\"> <META NAME=\"robots\" CONTENT=\"none\"> <META NAME=\"robots\" CONTENT=\"noindex,nofollow\"> <META NAME=\"robots\" CONTENT=\"noarchive\">"<<endl;
            aus<<"<?php $page=$_SERVER['PHP_SELF']; $sec=\"" << int(refresh) << "\"; header(\"Refresh: $sec; url=$page\"); ?>"<<endl;
            aus << "</head> <!?php print strftime('%c'); ?> <pre>"<<endl;
            aus <<  oss.str() << endl;
            for(uint i=0;i<vlines.size();i++)
                if(i>4 && i<vlines.size()-1) aus << vlines.at(i) << endl;
            aus << "</pre> </body> </html>"<<endl;
            aurostd::stringstream2file(aus,WRITE);
        }

        if(Tmax>=warning_beep) { aurostd::execute(XHOST.command("beep")+" -l 100 -f "+aurostd::utype2string<double>(50*(Tmax-Tzero)));}
        if(Tmax>=warning_halt) {
            aurostd::execute(XHOST.command("beep")+" -f 1000");aurostd::execute(XHOST.command("beep")+" -f 1500");
            aurostd::execute(XHOST.command("halt"));aurostd::execute(XHOST.command("beep")+" -f 2000");
            aurostd::execute(XHOST.command("beep")+" -f 2500");}

        //   if(maxmem>0.0 && maxmem<100.0) oss << aurostd::utype2string<double>(AFLOW_checkMEMORY("vasp",maxmem),4);
        //   if(maxmem>0.0 && maxmem<100.0) AFLOW_checkMEMORY("vasp",maxmem);
        //    if(maxmem>0.0 && maxmem<100.0) AFLOW_checkMEMORY("aflow",maxmem);
        // if(maxmem>0.0 && maxmem<100.0) AFLOW_checkMEMORY("clamd",maxmem);

        cout << oss.str();// cerr << oss.str();
        oss.str(std::string());
        cout << endl;cout.flush();//cerr << endl;cerr.flush();

        if(!RUNSTAT && !RUNBAR) break;
        sleep(refresh*(1.0-(Tmax-Tzero)/40.0));
    }
    return 1;
}

// ***************************************************************************
// AFLOW_monitor
// ***************************************************************************
uint AFLOW_monitor(vector<string> argv) {
    cout << "MMMMM  Aflow: starting AFLOW_monitor" << endl;
    cerr << "MMMMM  Aflow: starting AFLOW_monitor" << endl;

    // double maxmem=aurostd::args2attachedutype<double>(argv,"--mem=","--maxmem=",double (95.0/XHOST.CPU_Cores));
    double maxmem=aurostd::args2attachedutype<double>(argv,"--mem=|--maxmem=",double (3.0));
    maxmem=round(10*maxmem)/10; // to get easy numbers.

    vector<string> argv_local;
    argv_local.push_back(argv.at(0));
    argv_local.push_back("--runstat");
    argv_local.push_back("--refresh=180");
    argv_local.push_back("--mem="+aurostd::utype2string<double>(maxmem));
    argv_local.push_back("--warning_beep=57");

    //  argv_local.push_back();

    return AFLOW_getTEMP(argv_local);
}

// ***************************************************************************
// CheckAFLOWLIBMaterialServer
// ***************************************************************************
bool CheckMaterialServer(void) {return CheckMaterialServer("");}
bool CheckMaterialServer(string message) {
    if(XHOST.hostname==XHOST.AFLOW_MATERIALS_SERVER) return TRUE;
    if(XHOST.hostname==XHOST.AFLOW_WEB_SERVER) return TRUE;
    if(XHOST.hostname=="habana") return TRUE;
    cerr << "AFLOW ERROR: Your machine is \"" << XHOST.hostname << "\"." << endl;
    if(message.length()>0) cerr << "AFLOW ERROR: command \"" << message << "\" can run only on \"" << XHOST.AFLOW_MATERIALS_SERVER << "\" or \"" << XHOST.AFLOW_WEB_SERVER << "\"." << endl;
    else cerr << "AFLOW ERROR: the procedure can run only on \"" << XHOST.AFLOW_MATERIALS_SERVER << "\" or \"" << XHOST.AFLOW_WEB_SERVER << "\"." << endl;
    exit(0);
    return FALSE;
}

// ***************************************************************************
// aflow_get_time_string
// ***************************************************************************
string aflow_get_time_string(void) {
#ifdef ALPHA
    ostringstream aus;
    string OUT;
    aus<<"date | sed \"s/ /_/g\" > "+XHOST.Tmpfs+"/date."<< XHOST.ostrPID.str() << " " <<endl;
    system(aus.str().c_str());
    ifstream FileAUS;
    string FileNameAUS=XHOST.Tmpfs+"/date"+XHOST.ostrPID.str();
    FileAUS.open(FileNameAUS.c_str(),std::ios::in);
    FileAUS >> OUT;
    FileAUS.clear();FileAUS.close();
    // return (char*) OUT.c_str();
    return string("NotAvailable \n");
#else
    long ltime=time(NULL);

    string date=string(ctime(&ltime));
    if(date.length()>0)
        if(date.at(date.length()-1)=='\n')
            date.erase(date.length()-1);
    return date;
#endif
}

// ***************************************************************************
// aflow_get_time_string_short
// ***************************************************************************
string aflow_get_time_string_short(void) {
    string date;
    vector<string> tokens;
    aurostd::string2tokens(aflow_get_time_string(),tokens);
    date="na";
    if(tokens.size()>4) {
        if(tokens.at(2).size()>1)
            date=tokens.at(4).substr(2,2)+tokens.at(1)+tokens.at(2);
        else
            date=tokens.at(4).substr(2,2)+tokens.at(1)+"0"+tokens.at(2);
        aurostd::StringSubst(date,"Jan","01");aurostd::StringSubst(date,"Feb","02");aurostd::StringSubst(date,"Mar","03");
        aurostd::StringSubst(date,"Apr","04");aurostd::StringSubst(date,"May","05");aurostd::StringSubst(date,"Jun","06");
        aurostd::StringSubst(date,"Jul","07");aurostd::StringSubst(date,"Aug","08");aurostd::StringSubst(date,"Sep","09");
        aurostd::StringSubst(date,"Oct","10");aurostd::StringSubst(date,"Nov","11");aurostd::StringSubst(date,"Dec","12");
        date=date.substr(0,6);
    }
    if(date.length()>0)
        if(date.at(date.length()-1)=='\n')
            date.erase(date.length()-1);
    return date;
}

// ***************************************************************************
// strPID
// ***************************************************************************
string strPID(void) {
    int PID=getpid();
    ostringstream oss;
    oss << PID;
    return (string) oss.str();
}

// ***************************************************************************
// Messages
// ***************************************************************************
double AFLOW_checkMEMORY(string progname,double memory) {
    vector<string> vps,tokens;string command;
    double maxmem=0.0;
    if(progname.empty()) aurostd::string2vectorstring(aurostd::execute2string("ps aux | grep -v \" 0.0  0.0 \" | grep "+XHOST.User),vps);
    else aurostd::string2vectorstring(aurostd::execute2string("ps aux | grep \""+progname+"\" | grep -v \" 0.0  0.0 \" | grep "+XHOST.User),vps);
    for(uint i=0;i<vps.size();i++) {
        aurostd::string2tokens(vps.at(i),tokens);
        if(tokens.size()>4) {
            if(aurostd::string2utype<double>(tokens.at(3))>maxmem) maxmem=aurostd::string2utype<double>(tokens.at(3));
            if(memory>0.0 && memory<100.0) {
                if(aurostd::string2utype<double>(tokens.at(3))>memory) {
                    command=string(XHOST.command("kill")+" -9 "+tokens.at(1));
                    aurostd::execute(command);
                    //	  cerr << endl << "AFLOW_checkMEMORY: killing(" << memory << ") = " << vps.at(i) << endl;
                    cout << endl << "AFLOW_checkMEMORY [date=" << aflow_get_time_string() << "]: kill(" << tokens.at(3) << ">" << aurostd::utype2string<double>(memory,4) << ") = [" << vps.at(i) << "]" << endl;
                }
            }
        }
    }
    return maxmem;
}

// ***************************************************************************
// Messages
// ***************************************************************************
pthread_mutex_t mutex_INIT_Message=PTHREAD_MUTEX_INITIALIZER;
string Message(string list2print) {
    pthread_mutex_lock(&mutex_INIT_Message);
    // pthread_mutex_unlock(&mutex_INIT_Message);
    stringstream oss("");
    if(aurostd::substring2bool(list2print,"user") || aurostd::substring2bool(list2print,"USER")) oss << " - [user=" << XHOST.User << "]";
    if(aurostd::substring2bool(list2print,"group") || aurostd::substring2bool(list2print,"GROUP")) oss << " - [group=" << XHOST.Group << "]";
    if(aurostd::substring2bool(list2print,"host") || aurostd::substring2bool(list2print,"HOST")) oss << " - [host=" << XHOST.hostname << "]";
    if(aurostd::substring2bool(list2print,"hostname") || aurostd::substring2bool(list2print,"HOSTNAME")) oss << " - [host=" << XHOST.hostname << "]";
    if(aurostd::substring2bool(list2print,"temperature")) if(init::GetTEMPs()) for(uint i=0;i<XHOST.vTemperatureCore.size();i++) {oss << (i==0?"- [temp(C)=":"") << XHOST.vTemperatureCore.at(i) << (i<XHOST.vTemperatureCore.size()-1?",":"]");}
    if(aurostd::substring2bool(list2print,"machine") || aurostd::substring2bool(list2print,"MACHINE")) oss << " - [host=" << XHOST.hostname << "]";
    if(list2print.empty() || aurostd::substring2bool(list2print,"time") || aurostd::substring2bool(list2print,"TIME")) oss << " - [date=" << aflow_get_time_string() << "]";
    if(aurostd::substring2bool(list2print,"date") || aurostd::substring2bool(list2print,"DATE")) oss << " - [date=" << aflow_get_time_string() << "]";
    //  if(XHOST.maxmem>0.0 && XHOST.maxmem<100.0)
    if(aurostd::substring2bool(list2print,"memory") && (XHOST.maxmem>0.0 && XHOST.maxmem<100)) oss << " - [mem=" << aurostd::utype2string<double>(AFLOW_checkMEMORY("vasp",XHOST.maxmem),4) << " (" << XHOST.maxmem << ")]"; // CO 170628 - slow otherwise!!!
    if(XHOST.vTemperatureCore.size()>0) if(max(XHOST.vTemperatureCore)>AFLOW_CORE_TEMPERATURE_BEEP) oss << " - [ERROR_TEMPERATURE=" << max(XHOST.vTemperatureCore) << ">" << AFLOW_CORE_TEMPERATURE_BEEP << "@ host=" << XHOST.hostname<< "]";
    // oss << endl;
    pthread_mutex_unlock(&mutex_INIT_Message);
    // do some killing
    //if(XHOST.maxmem>0.0 && XHOST.maxmem<100.0) AFLOW_checkMEMORY("vasp",XHOST.maxmem);  // CO 170628 - this is already run above, very slow
    // if(XHOST.maxmem>0.0 && XHOST.maxmem<100.0) AFLOW_checkMEMORY("aflow",XHOST.maxmem);
    // if(XHOST.maxmem>0.0 && XHOST.maxmem<100.0) AFLOW_checkMEMORY("clamd",XHOST.maxmem);
    return oss.str();
}

string Message(string str1,string list2print) {return string(" - "+str1+Message(list2print));}
//string Message(const _aflags& aflags) {return string(" - "+aflags.Directory + "\n");}
string Message(const _aflags& aflags) {
    string strout=" - [dir="+aflags.Directory+"]"+=Message("user,host,time");
    if(AFLOW_PTHREADS::FLAG) strout+=" - [thread="+aurostd::utype2string(aflags.AFLOW_PTHREADS_NUMBER)+"/"+aurostd::utype2string(AFLOW_PTHREADS::MAX_PTHREADS)+"]";
    return strout;}
    string Message(const _aflags& aflags,string list2print1,string list2print2) {
        stringstream strout;
        if(!list2print1.empty()) strout << " [dir=" << aflags.Directory << "]" << Message(list2print1);
        if(AFLOW_PTHREADS::FLAG) strout << " - [thread=" << aurostd::utype2string(aflags.AFLOW_PTHREADS_NUMBER) << "/" << aurostd::utype2string(AFLOW_PTHREADS::MAX_PTHREADS) << "]";
        if(!list2print2.empty()) strout << " ["  <<  list2print2 << "]";
        return strout.str();
    }

// ***************************************************************************
// AFLOW_BlackList
// ***************************************************************************
bool AFLOW_BlackList(string h) {
    // cerr << h << endl;
    // if(h=="nietzsche" || h=="nietzsche.mems.duke.edu" || h=="material.duke.edu") return TRUE;
    if(h=="blacklisted_hostname") return TRUE;
    //  if(h=="m6-11-6") return TRUE;
    return FALSE;
}

// ***************************************************************************
// init::ErrorOptions
// ***************************************************************************
namespace init {
    bool ErrorOption(ostream &oss,const string& options, const string& routine,vector<string> vusage) {
        vector<string> tokens_options;
        aurostd::string2tokens(options,tokens_options,",");

        oss << "ERROR: " << routine << ":" << endl;
        oss << "       Wrong number/type of input parameters! (" << tokens_options.size() << ")" << endl;
        string usage="       Usage: ";
        for(uint i=0;i<vusage.size();i++) {
            if(aurostd::substring2bool(vusage.at(i),"options:")) usage="              ";
            if(vusage.at(i)!="") oss << usage << vusage.at(i) << endl;
        }
        oss << "       options=[" << options << "]" << endl;
        return TRUE;
    }
    bool ErrorOption(ostream &oss,const string& options, const string& routine,string usage) {
        vector<string> vusage;
        aurostd::string2vectorstring(usage,vusage);
        return ErrorOption(oss,options,routine,vusage);
    }
}

#endif

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************
