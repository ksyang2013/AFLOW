// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// contains routines to run PTHREADS

#ifndef _AFLOW_PTHREADS_CPP
#define _AFLOW_PTHREADS_CPP
#include "aflow.h"
#include "aflow_pflow.h"

// #define  AFLOW_PTHREADS_MULTISH_PREEMPTIVE_
#define  AFLOW_PTHREADS_MULTISH_TIMESHARING_
//#define  AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_
//#define  AFLOW_PTHREADS::MULTISH_TIMESHARING_CONCURRENT_

#define _KBIN_SEEK_THREAD_SLEEP_  33
#define _KBIN_START_THREAD_SLEEP_ 2
#define _KBIN_FLUSH_THREAD_SLEEP_ 1

//#define _PTHREAD_FLUSH_TIME_ 1

using aurostd::substring2bool;

namespace AFLOW_PTHREADS {
  bool FLAG;                                     // run pthread YES/NO
  int MAX_PTHREADS;                              // how many MAX threads I can use  default or --np
  int RUNNING;                                   // how many threads are actually running
  pthread_t vpthread[MAX_ALLOCATABLE_PTHREADS];  // the actual thread
  int viret[MAX_ALLOCATABLE_PTHREADS];           // the thread runnings
  bool vpthread_busy[MAX_ALLOCATABLE_PTHREADS];  // is the thread busy
  bool MULTISH_TIMESHARING_SEQUENTIAL_=FALSE;
  bool MULTISH_TIMESHARING_CONCURRENT_=TRUE;
}

#define CPU_File     string("/proc/cpuinfo")
#define CPU_String   string("cpu MHz")

// ***************************************************************************
// AFLOW_PTHREADS::GetTotalCPUs
// ***************************************************************************
// This function returns the max number of CPUS by interrogating as
// cat /proc/cpuinfo | grep -c "cpu MHz"
// if the file is not found, it returns CPUs=1
namespace AFLOW_PTHREADS {
  int GetTotalCPUs(void) {
    int CPU_Cores=sysconf(_SC_NPROCESSORS_ONLN);
    if(CPU_Cores<1) CPU_Cores=1;
    return CPU_Cores;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Check_Threads
// **************************************************************************
// This function checks the input argv and set up the proper multithread
// parameters (SC Dec07)
namespace AFLOW_PTHREADS {
  bool Check_Threads(vector<string> argv,const bool& VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(VERBOSE) {;} // dummy load
    AFLOW_PTHREADS::FLAG=FALSE;
    AFLOW_PTHREADS::MAX_PTHREADS=1;
    bool fnp1=aurostd::args2attachedutype<bool>(argv,"--np=",FALSE);
    bool fnp=aurostd::args2attachedflag(argv,"--np=");
    bool fnpmax=aurostd::args2flag(argv,"--npmax");
    bool multish=aurostd::args2flag(argv,"--multish");
    if(!fnp && !fnpmax) {AFLOW_PTHREADS::MAX_PTHREADS=1;}
    if(fnp && !fnpmax)  {AFLOW_PTHREADS::MAX_PTHREADS=aurostd::args2attachedutype<int>(argv,"--np=",0);};
    if(!fnp && fnpmax)  {AFLOW_PTHREADS::MAX_PTHREADS=AFLOW_PTHREADS::GetTotalCPUs();}
    if(fnp && fnpmax)   {AFLOW_PTHREADS::MAX_PTHREADS=aurostd::args2attachedutype<int>(argv,"--np=",0);};
    if(multish && !fnp && !fnpmax) {AFLOW_PTHREADS::MAX_PTHREADS=AFLOW_PTHREADS::GetTotalCPUs();}

    if(AFLOW_PTHREADS::MAX_PTHREADS>1) {
      AFLOW_PTHREADS::FLAG=TRUE;
      if(AFLOW_PTHREADS::MAX_PTHREADS>MAX_ALLOCATABLE_PTHREADS) AFLOW_PTHREADS::MAX_PTHREADS=MAX_ALLOCATABLE_PTHREADS;
      if(LDEBUG) cerr << "AAAAA  AFLOW THREADED VERSION  threads=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
    }
    if(AFLOW_PTHREADS::MAX_PTHREADS<=1) {
      AFLOW_PTHREADS::FLAG=FALSE;
      AFLOW_PTHREADS::MAX_PTHREADS=1;
      if(LDEBUG) cerr << "AAAAA  AFLOW SERIAL VERSION threads=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
    }
    //  AFLOW_PTHREADS::FLAG=TRUE;
    if(LDEBUG) cerr << "AFLOW_PTHREADS::Check_Threads: fnp=" << fnp << endl;
    if(LDEBUG) cerr << "AFLOW_PTHREADS::Check_Threads: fnp1=" << fnp1 << endl;
    if(LDEBUG) cerr << "AFLOW_PTHREADS::Check_Threads: fnpmax=" << fnpmax << endl;
    if(LDEBUG) cerr << "AFLOW_PTHREADS::Check_Threads: multish=" << multish << endl;
    if(LDEBUG) cerr << "AFLOW_PTHREADS::Check_Threads: AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
    if(LDEBUG) cerr << "AFLOW_PTHREADS::Check_Threads: AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << endl;
    //  if(LDEBUG) exit(0); // for debug
    return AFLOW_PTHREADS::FLAG;
  }
} // namespace AFLOW_PTHREADS

namespace AFLOW_PTHREADS {
  bool Check_Threads_WrapperNP(vector<string> argv,uint size_to_check,const bool& VERBOSE) {
    ostringstream aus;
    AFLOW_PTHREADS::FLAG=TRUE;
    AFLOW_PTHREADS::Check_Threads(argv,VERBOSE);
    if(AFLOW_PTHREADS::MAX_PTHREADS>(int) size_to_check) {
      if(VERBOSE) aus << "AFLOW multi_XXXX (" << string(AFLOW_VERSION)  << "): WARNING (threads>commands) => threads=commands" << endl;
      AFLOW_PTHREADS::MAX_PTHREADS=(int) size_to_check;
    }
    if(AFLOW_PTHREADS::MAX_PTHREADS>=2) {if(VERBOSE) aus << "AFLOW multi_XXXX (" << string(AFLOW_VERSION)  << "): threads=" << AFLOW_PTHREADS::MAX_PTHREADS << "  -  commands=" << size_to_check << endl;}
    if(AFLOW_PTHREADS::MAX_PTHREADS<=1) {if(VERBOSE) aus << "AFLOW multi_XXXX (" << string(AFLOW_VERSION)  << "): serial  -  commands=" << size_to_check << endl;}
    //  if(VERBOSE)
    aurostd::PrintMessageStream(aus,XHOST.QUIET);

    return AFLOW_PTHREADS::FLAG;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Clean_Threads
// **************************************************************************
// This function clears the busy-ness of all the pthreads flags
namespace AFLOW_PTHREADS {
  void Clean_Threads(void) {
    for(uint ithread=0;ithread<MAX_ALLOCATABLE_PTHREADS;ithread++)         // clean threads
      AFLOW_PTHREADS::vpthread_busy[ithread]=FALSE;          // clean threads
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::No_Threads
// **************************************************************************
// This function removes threads
namespace AFLOW_PTHREADS {
  void No_Threads(void) {
    AFLOW_PTHREADS::FLAG=FALSE;
    AFLOW_PTHREADS::MAX_PTHREADS=1;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Available_Free_Threads
// **************************************************************************
// This function return TRUE and a free PTHREAD index if a free thread is
// available
namespace AFLOW_PTHREADS {
  bool Available_Free_Threads(int &fthread) {
    bool free_thread=FALSE;
    fthread=-1;
    for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
      if(AFLOW_PTHREADS::vpthread_busy[ithread]==FALSE) {
	free_thread=TRUE;
	fthread=ithread;
      }
    }
    return free_thread;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Wait_Available_Free_Threads
// **************************************************************************
// This function return TRUE and a free PTHREAD index if a free thread is
// available
namespace AFLOW_PTHREADS {
  bool Wait_Available_Free_Threads(int &fthread,const double& pthread_wait,const bool& VERBOSE) {
    ostringstream aus;
    bool free_thread=FALSE;
    while(!free_thread)  {            // waiting for free thread to start
      free_thread=AFLOW_PTHREADS::Available_Free_Threads(fthread);
      if(!free_thread) {              // do something else !
	if(VERBOSE) aus << "MMMMM  Aflow: MULTI-THREADED: Waiting for free threads: " << pthread_wait << " seconds " << endl;
	if(VERBOSE) aurostd::PrintMessageStream(aus,FALSE);
	aurostd::Sleep((int) pthread_wait);
      }
      if(free_thread) {
	if(VERBOSE) aus << "MMMMM  Aflow: MULTI-THREADED: Found free thread  fthread=" << fthread << " - " << endl;
	if(VERBOSE) aurostd::PrintMessageStream(aus,FALSE);
      }
    }
    return free_thread;
  }
} // namespace AFLOW_PTHREADS

namespace AFLOW_PTHREADS {
  bool Wait_Available_Free_Threads(int &fthread,const bool& VERBOSE) {
    return AFLOW_PTHREADS::Wait_Available_Free_Threads(fthread,_KBIN_SEEK_THREAD_SLEEP_,VERBOSE);
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// KBIN::RUN_Directory_PTHREADS
// **************************************************************************
// Interfaces for KBIN_RUN_Directory

namespace KBIN {
  typedef struct {
    _aflags *paflags;     // FOR KBIN (ALL)
    int      itbusy;      // FOR KBIN (ALL)
    bool     VERBOSE;     // FOR KBIN (ALL)
    string   command;     // FOR MULTISH PREEMPTIVE
    int      ITHREAD;     // FOR MULTISH_TIMESHARING
    int      THREADS_MAX; // FOR MULTISH_TIMESHARING
    deque<string> *dcmds; // FOR MULTISH_TIMESHARING
  } _threaded_params;
} // namespace KBIN

KBIN::_threaded_params params[MAX_ALLOCATABLE_PTHREADS];
_aflags taflags[MAX_ALLOCATABLE_PTHREADS];
pthread_mutex_t mutex_PTHREAD=PTHREAD_MUTEX_INITIALIZER;

namespace KBIN {
  void RUN_Directory_PTHREADS(_aflags &aflags) {
    int ithread=aflags.AFLOW_PTHREADS_NUMBER;
    if(ithread<0) {cerr << "ERROR KBIN::RUN_Directory_PTHREADS: ithread<0  ithread=" << ithread << endl;exit(0);};
    if(ithread>=MAX_ALLOCATABLE_PTHREADS) {cerr << "ERROR KBIN::RUN_Directory_PTHREADS: ithread>=MAX_ALLOCATABLE_PTHREADS  ithread=" << ithread << endl;exit(0);};
    if(ithread>=AFLOW_PTHREADS::MAX_PTHREADS) {cerr << "ERROR KBIN::RUN_Directory_PTHREADS: ithread>=AFLOW_PTHREADS::MAX_PTHREADS  ithread=" << ithread << endl;exit(0);};
    taflags[ithread]=aflags;
    params[ithread].paflags=&taflags[ithread];
    params[ithread].command="";//command;
    params[ithread].itbusy=ithread;
    AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,KBIN::_threaded_interface_RUN_Directory, (void*)&params[ithread]);
    aurostd::Sleep(_KBIN_START_THREAD_SLEEP_);
  }
} // namespace KBIN

namespace KBIN {
  void *_threaded_interface_RUN_Directory(void *ptr) {
    KBIN::_threaded_params* pparams;
    pparams=(KBIN::_threaded_params*) ptr;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
    AFLOW_PTHREADS::RUNNING++;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP_);
    aurostd::execute(XHOST.command("aflow")+" --run=1 --DIRECTORY="+(*pparams->paflags).Directory);  // run it OUTSIDE
    //  KBIN_RUN_Directory((*pparams->paflags)); // RUN IT INSIDE
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
    AFLOW_PTHREADS::RUNNING--;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP_);
    return NULL;
  }
} // namespace KBIN

// **************************************************************************
// FUNCTION AFLOW_PTHREADS_MULTISH_PREEMPTIVE_
// **************************************************************************
// Interfaces for AFLOW_PTHREADS_MULTISH_PREEMPTIVE_

#ifdef AFLOW_PTHREADS_MULTISH_PREEMPTIVE_
#warning "aflow_pthreads.cpp with AFLOW_PTHREADS_MULTISH_PREEMPTIVE_"

namespace KBIN {
  void *_threaded_interface_MULTIRUN_sh(void *ptr) {
    KBIN::_threaded_params* pparams;
    pparams=(KBIN::_threaded_params*) ptr;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
    AFLOW_PTHREADS::RUNNING++;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP);
    cout << pparams->itbusy << "  - " << pparams->command << endl;
    aurostd::execute(pparams->command);
    // KBIN_RUN_Directory((*pparams->paflags));
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
    AFLOW_PTHREADS::RUNNING--;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP);
    return NULL;
  }
} // namespace KBIN

namespace KBIN {
  void MULTIRUN_PTHREADS(_aflags &aflags,string command) {
    int ithread=aflags.AFLOW_PTHREADS_NUMBER;
    if(ithread<0) {cerr << "ERROR KBIN::MULTIRUN_PTHREADS: ithread<0  ithread=" << ithread << endl;exit(0);};
    if(ithread>=MAX_ALLOCATABLE_PTHREADS) {cerr << "ERROR KBIN::MULTIRUN_PTHREADS: ithread>=MAX_ALLOCATABLE_PTHREADS  ithread=" << ithread << endl;exit(0);};
    if(ithread>=AFLOW_PTHREADS::MAX_PTHREADS) {cerr << "ERROR KBIN::MULTIRUN_PTHREADS: ithread>=AFLOW_PTHREADS::MAX_PTHREADS  ithread=" << ithread << endl;exit(0);};
    taflags[ithread]=aflags;
    params[ithread].paflags=&taflags[ithread];
    params[ithread].command=command;
    params[ithread].itbusy=ithread;
    AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,KBIN::_threaded_interface_MULTIRUN_sh, (void*)&params[ithread]);
    aurostd::Sleep(_KBIN_START_THREAD_SLEEP_);
  }
} // namespace KBIN

namespace AFLOW_PTHREADS {
  bool MULTI_sh(vector<string> argv) {
    ostringstream aus;
    _aflags aflags;
    // [OBSOLETE]    string file_name=aurostd::args2string(argv,"--FILE|--F|--f","xxxx");
    string file_name=XHOST.vflag_control.getattachedscheme("FILE");
    bool free_thread;int ithread=0;
    bool VERBOSE=FALSE;

    AFLOW_PTHREADS::FLAG=TRUE;
    AFLOW_PTHREADS::MAX_PTHREADS=PTHREAD_DEFAULT; // safety...

    if(!aurostd::FileExist(file_name)) {aus << "EEEEE FILE_NOT_FOUND = " << file_name << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);}
    if( aurostd::FileEmpty(file_name)) {aus << "EEEEE FILE_EMPTY = " << file_name << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);}
    aus << "MMMMM Loading File = " << file_name << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
    vector<string> vcmds;
    vcmds.clear();
    aurostd::file2vectorstring(file_name,vcmds);
    aus << "MMMMM Loaded Lines = " << vcmds.size() << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);

    AFLOW_PTHREADS::Clean_Threads();                                  // clean threads

    for(uint i=0;i<vcmds.size();i++) {
      // loop this
      free_thread=AFLOW_PTHREADS::Wait_Available_Free_Threads(ithread,0,VERBOSE);        // WAIT A WHILE !!
      if(free_thread) {
	//  aus << "MMMMM Aflow: Found subdirectory to run " << aflags.Directory<< " - " << XHOST.hostname << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	aus << "MMMMM Aflow: MULTI-THREADED: Starting  pthread_free=" << ithread << "  pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " vline=" << i << " - " << XHOST.hostname << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
	aflags.AFLOW_PTHREADS_NUMBER=ithread;
	AFLOW_PTHREADS::vpthread_busy[ithread]=TRUE;
	KBIN::MULTIRUN_PTHREADS(aflags,vcmds.at(i));
	AFLOW_PTHREADS::vpthread_busy[ithread]=FALSE;
      }
      //    cout << free_thread << " " << ithread << endl;
    }

    aus << "MMMMM  Aflow: MULTI-THREADED: FLUSHING PTHREADS - " << XHOST.hostname << endl;
    aurostd::PrintMessageStream(aus,XHOST.QUIET);
    for(ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
      if(AFLOW_PTHREADS::vpthread_busy[ithread]==TRUE) {
	aus << "MMMMM  Aflow: MULTI-THREADED: Flushing   pthread=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << XHOST.hostname << endl;
	aurostd::PrintMessageStream(aus,XHOST.QUIET);
	pthread_join(thread[ithread],NULL);
      }
    return TRUE;
  }
} // namespace AFLOW_PTHREADS

#endif //  AFLOW_PTHREADS_MULTISH_PREEMPTIVE_

// **************************************************************************
// FUNCTION AFLOW_PTHREADS_MULTISH_TIMESHARING_
// **************************************************************************
// Interfaces for AFLOW_PTHREADS_MULTISH_TIMESHARING_

#ifdef AFLOW_PTHREADS_MULTISH_TIMESHARING_
#warning "aflow_pthreads.cpp with AFLOW_PTHREADS::MULTISH_TIMESHARING_"

namespace AFLOW_PTHREADS {
  void *_threaded_COMMANDS(void *ptr) {
    if(AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_) {
      //cerr << "AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_ AFLOW_PTHREADS::_threaded_COMMANDS" << endl;
      KBIN::_threaded_params* pparams=(KBIN::_threaded_params*) ptr;
      string command;
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
      AFLOW_PTHREADS::RUNNING++;
      if((pparams->VERBOSE)) {pthread_mutex_lock(&mutex_PTHREAD);cout << "AFLOW_PTHREADS::_threaded_COMMANDS " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;pthread_mutex_unlock(&mutex_PTHREAD);}
      for(uint ithread=(pparams->ITHREAD);ithread<(*pparams->dcmds).size();ithread+=(pparams->THREADS_MAX)) {
	if((pparams->VERBOSE)) {pthread_mutex_lock(&mutex_PTHREAD);cout <<  (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << " " << ithread << " " << (*pparams->dcmds).at(ithread) << endl;pthread_mutex_unlock(&mutex_PTHREAD);}
	command=(*pparams->dcmds).at(ithread);
	aurostd::execute(command);
	pthread_mutex_lock(&mutex_PTHREAD);cout << command << endl;cout.flush();pthread_mutex_unlock(&mutex_PTHREAD);
	// if((pparams->VERBOSE)) {cout << char('A'+(pparams->ITHREAD));cout.flush();}
      }
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
      AFLOW_PTHREADS::RUNNING--;
      // aurostd::Sleep(_KBIN_START_THREAD_SLEEP_);
      return NULL;
    }
    if(AFLOW_PTHREADS::MULTISH_TIMESHARING_CONCURRENT_) { // it bombs and I do not know why...
      bool FRONT=TRUE;
      KBIN::_threaded_params* pparams=(KBIN::_threaded_params*) ptr;
      string command;
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
      AFLOW_PTHREADS::RUNNING++;
      while((*pparams->dcmds).size()>0) {
	pthread_mutex_lock(&mutex_PTHREAD);
	if(FRONT)  {command=(*pparams->dcmds).at(0);(*pparams->dcmds).pop_front();}  // from the front
	if(!FRONT) {command=(*pparams->dcmds).at((*pparams->dcmds).size()-1);(*pparams->dcmds).pop_back();}  // from the back
	pthread_mutex_unlock(&mutex_PTHREAD);
	if((pparams->VERBOSE)) {pthread_mutex_lock(&mutex_PTHREAD);cout <<  (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << ": " << command << endl;pthread_mutex_unlock(&mutex_PTHREAD);}
	//    aurostd::Sleep((int) _KBIN_START_THREAD_SLEEP_);
	//    cerr << "[1]" << endl;
	//   cout << command << endl;cout.flush();
	aurostd::execute(command);
	// if((pparams->VERBOSE)) {cout << char('A'+(pparams->ITHREAD));cout.flush();}
      }
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
      AFLOW_PTHREADS::RUNNING--;
      // cerr << "[2]" << endl;
      return NULL;
    }
    if(!AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_ && !AFLOW_PTHREADS::MULTISH_TIMESHARING_CONCURRENT_) {
      cerr << "ERROR: AFLOW_PTHREADS::*_threaded_COMMANDS: you must specify MULTISH_TIMESHARING_SEQUENTIAL_ of MULTISH_TIMESHARING_CONCURRENT_" << endl;
      exit(0);
    }
    return NULL;
  }
} // namespace AFLOW_PTHREADS


namespace AFLOW_PTHREADS {
  bool MULTI_sh(vector<string> argv) {
    ostringstream aus;
    _aflags aflags;
    // [OBSOLETE] string file_name=aurostd::args2string(argv,"--FILE|--F|--f","xxxx");
    string file_name=XHOST.vflag_control.getattachedscheme("FILE");
    if(file_name.empty() || file_name=="--f") file_name=argv.at(argv.size()-1);
    //  bool free_thread;int ithread=0;
    bool VERBOSE=FALSE;

    if(!aurostd::FileExist(file_name)) {aus << "EEEEE  FILE_NOT_FOUND = " << file_name  << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);}
    if( aurostd::FileEmpty(file_name)) {aus << "EEEEE  FILE_EMPTY = " << file_name  << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);}
    if(VERBOSE) {aus << "MMMMM  Loading File = " << file_name << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    vector<string> vcmds;vcmds.clear();
    aurostd::file2vectorstring(file_name,vcmds);
    
    AFLOW_PTHREADS::Check_Threads_WrapperNP(argv,vcmds.size(),VERBOSE);   // check treads from NP
    aurostd::multithread_execute(vcmds,AFLOW_PTHREADS::MAX_PTHREADS,VERBOSE);
    cout << endl;
    return TRUE;
  }
} // AFLOW_PTHREADS

namespace aurostd {
  bool multithread_execute(deque<string> dcmds,int _NUM_THREADS,bool VERBOSE) {
    int NUM_THREADS=_NUM_THREADS;                                          // SAFETY
    if((int) dcmds.size()<=NUM_THREADS) NUM_THREADS=(uint) dcmds.size();   // SAFETY

    if(NUM_THREADS<=1) {                                                   // run singular
      for(uint i=0;i<dcmds.size();i++)                                     // run singular
	aurostd::execute(dcmds.at(i));                                     // run singular
    }
    if(NUM_THREADS>=2) {                                                   // multithread
      AFLOW_PTHREADS::FLAG=TRUE;AFLOW_PTHREADS::MAX_PTHREADS=NUM_THREADS;  // prepare
      if(AFLOW_PTHREADS::MAX_PTHREADS>MAX_ALLOCATABLE_PTHREADS) AFLOW_PTHREADS::MAX_PTHREADS=MAX_ALLOCATABLE_PTHREADS; // check max
      _aflags aflags;
      AFLOW_PTHREADS::Clean_Threads();                                     // clean threads
      KBIN::_threaded_params params[MAX_ALLOCATABLE_PTHREADS];             // prepare
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {  // prepare loop
	params[ithread].paflags=&aflags;                                   // prepare params
	params[ithread].ITHREAD=ithread;                                   // prepare params
	params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;                    // prepare params
	//    cerr << AFLOW_PTHREADS::MAX_PTHREADS << endl;
	params[ithread].dcmds=&dcmds;                                      // prepare params
	params[ithread].itbusy=ithread;                                    // prepare params
	params[ithread].VERBOSE=VERBOSE;                                   // prepare params
      }
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)                          // run threads
	AFLOW_PTHREADS::viret[ithread]=pthread_create(&AFLOW_PTHREADS::vpthread[ithread],NULL,AFLOW_PTHREADS::_threaded_COMMANDS,(void*)&params[ithread]); // run threads
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)                          // flush
	pthread_join(AFLOW_PTHREADS::vpthread[ithread],NULL);                                    // flush
    }
    return TRUE;
  }
} // namespace aurostd

namespace aurostd {
  bool multithread_execute(vector<string> vcmds,int _NUM_THREADS,bool VERBOSE) {
    int NUM_THREADS=_NUM_THREADS;                                          // SAFETY
    if((int) vcmds.size()<=NUM_THREADS) NUM_THREADS=(uint) vcmds.size();   // SAFETY

    deque<string> dcmds;dcmds.clear();
    for(uint i=0;i<vcmds.size();i++) dcmds.push_back(vcmds.at(i));      // copy
    return aurostd::multithread_execute(dcmds,NUM_THREADS,VERBOSE);
  }
} // namespace aurostd

#endif //  AFLOW_PTHREADS_MULTISH_TIMESHARING_

// ***************************************************************************
// MultiThread Execute vectors/deque of Strings
// ***************************************************************************
// adding something to aurostd
namespace aurostd {
  bool multithread_execute(deque<string> vcommand,int NUM_THREADS) {
    return multithread_execute(vcommand,NUM_THREADS,FALSE);
  }
  bool multithread_execute(vector<string> vcommand,int NUM_THREADS) {
    return multithread_execute(vcommand,NUM_THREADS,FALSE);
  }
  bool multithread_execute(deque<string> vcommand) {
    AFLOW_PTHREADS::MAX_PTHREADS=AFLOW_PTHREADS::GetTotalCPUs();
    return multithread_execute(vcommand,AFLOW_PTHREADS::MAX_PTHREADS,FALSE);
  }
  bool multithread_execute(vector<string> vcommand) {
    AFLOW_PTHREADS::MAX_PTHREADS=AFLOW_PTHREADS::GetTotalCPUs();
    return multithread_execute(vcommand,AFLOW_PTHREADS::MAX_PTHREADS,FALSE);
  }
} // namespace AFLOW_PTHREADS

// ***************************************************************************
// AFLOW_PTHREADS::MULTI_zip
// **************************************************************************
namespace AFLOW_PTHREADS {
  bool MULTI_zip(vector<string> argv) {
    ostringstream aus;
    _aflags aflags;
    bool VERBOSE=FALSE;
    vector<string> vdirs;vdirs.clear();
    // LOAD FILE
    // [OBSOLETE]    if(aurostd::args2flag(argv,"--FILE|--F|--f")) {
    // [OBSOLETE]      string file_name=aurostd::args2string(argv,"--FILE|--F|--f","xxxx");
    if(XHOST.vflag_control.flag("FILE")) {
      string file_name=XHOST.vflag_control.getattachedscheme("FILE");
      // if(file_name.empty() || file_name=="--f") file_name=argv.at(argv.size()-1);
      if(!aurostd::FileExist(file_name)) {aus << "EEEEE  FILE_NOT_FOUND = " << file_name  << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);}
      if( aurostd::FileEmpty(file_name)) {aus << "EEEEE  FILE_EMPTY = " << file_name  << Message(aflags,"user,host,time") << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);exit(0);}
      if(VERBOSE) {aus << "MMMMM  Loading File = " << file_name << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
      aurostd::file2vectorstring(file_name,vdirs);
    }
    // LOAD DIRECTORIES
    // [OBSOLETE] if(aurostd::args2flag(argv,"--DIRECTORY|--D|--d")) {
    // [OBSOLETE]  vdirs=aurostd::args2vectorstring(argv,"--DIRECTORY|--D|--d","./");
    // [OBSOLETE] }
    if(XHOST.vflag_control.flag("VDIR")) {
      aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VDIR"),vdirs,",");
    }
    for(uint i=0;i<vdirs.size();i++) {
      aurostd::StringSubst(vdirs.at(i),_AFLOWLOCK_,"");              // remove _AFLOWLOCK_
      aurostd::StringSubst(vdirs.at(i),_AFLOWIN_,"");                // remove _AFLOWIN_
      aurostd::StringSubst(vdirs.at(i),"OUTCAR.relax2.bz2","");      // remove OUTCAR.relax2.bz2
      aurostd::StringSubst(vdirs.at(i),"OSZICAR.relax2.bz2","");     // remove OSZICAR.relax2.bz2
      aurostd::StringSubst(vdirs.at(i),"EIGENVAL.relax2.bz2","");    // remove EIGENVAL.relax2.bz2
      aurostd::StringSubst(vdirs.at(i),"OUTCAR.static.bz2","");      // remove OUTCAR.static.bz2
      aurostd::StringSubst(vdirs.at(i),"OSZICAR.static.bz2","");     // remove OSZICAR.static.bz2
      aurostd::StringSubst(vdirs.at(i),"EIGENVAL.static.bz2","");    // remove EIGENVAL.static.bz2
      aurostd::StringSubst(vdirs.at(i),"OUTCAR.bands.bz2","");       // remove OUTCAR.bands.bz2
      aurostd::StringSubst(vdirs.at(i),"OSZICAR.bands.bz2","");      // remove OSZICAR.bands.bz2
      aurostd::StringSubst(vdirs.at(i),"EIGENVAL.bands.bz2","");     // remove EIGENVAL.bands.bz2
      //    cerr << vdirs.at(i) << endl;
    }
    //  cerr << vdirs.size() << endl; exit(0);

    int size=aurostd::args2attachedutype<int>(argv,"--size=",(int) (100));if(size<=0) size=1;
    string prefix=aurostd::args2attachedstring(argv,"--prefix=",(string) "m");
    bool flag_ADD=aurostd::args2flag(argv,"--add");

    //  cerr << prefix << endl; exit(0);

    int numzipsCE=(int) ceil(((double) vdirs.size())/((double) size));
    int numzipsFL=(int) floor(((double) vdirs.size())/((double) size));
    int numzips=0;
    if(aurostd::args2flag(argv,"--modonly")) {numzips=numzipsFL;} else {numzips=numzipsCE;}

    uint ishift=1,i=0;
    vector<string> vcommands;
    string command;
    if(!flag_ADD) while(aurostd::FileExist(prefix+"_"+aurostd::utype2string(i+ishift)+".zip")) ishift++;

    // delete POTCARs and AECCARs    
    vector<string> vremove;
    aurostd::string2tokens("POTCAR.relax1.bz2,POTCAR.relax2.bz2,POTCAR.static.bz2,POTCAR.bands.bz2,AECCAR1.static.bz2,AECCAR0.bands.bz2,AECCAR1.bands.bz2,AECCAR2.bands.bz2,core,core.bz2",vremove,",");
    for(uint j=0;j<vdirs.size();j++) {
      for(uint i=0;i<vremove.size();i++) {
	if(aurostd::FileExist(vdirs.at(j)+"/"+vremove.at(i))) {
	  aurostd::RemoveFile(vdirs.at(j)+"/"+vremove.at(i));	//	if(AFLOWLIB_VERBOSE) 
	  if(VERBOSE) cerr << "AFLOW_PTHREADS::MULTI_zip: REMOVED = " << aurostd::CleanFileName(vdirs.at(j)+"/"+vremove.at(i)) << endl;
	}
      }
    }
    
    // cerr << numzips << endl; exit(0);
    //  cerr << sysconf(_SC_ARG_MAX)  << endl;//exit(0);
    for(i=0;i<(uint) numzips;i++) {
      command="zip -9rmv "+prefix+"_"+aurostd::utype2string(i+ishift)+".zip";
      // command="zip -9rv "+prefix+"_"+aurostd::utype2string(i+ishift)+".zip";
      for(uint j=0;j<(uint) size;j++) {
	if(i*size+j < vdirs.size()) {
	  command+=" "+vdirs.at(i*size+j);
	}
      }
      command+=" | grep aflow.in";
      vcommands.push_back(command);
      //  cerr << command.size() << endl;
      //  cout << "[" << i << "] - " << command << endl;
    }

    //  for(uint i=0;i<vcommands.size();i++) aurostd::execute(vcommands.at(i));
    aurostd::multithread_execute(vcommands,PTHREADS_DEFAULT);

    return TRUE;
  }
} // namespace AFLOW_PTHREADS

// ***************************************************************************
// AFLOW_PTHREADS::MULTI_bzip2
// **************************************************************************
namespace AFLOW_PTHREADS {
  bool MULTI_bzip2(vector<string> argv) {
    // aflow --multibzip2 [--np XX | npmax | nothing] --F[ILE] file1 file2 file3 ....
    // load files
    bool VERBOSE=TRUE;
    vector<string> vf;
    // [OBSOLETE] if(aurostd::args2flag(argv,"--FILE|--F|--f")) vf=aurostd::args2vectorstring(argv,"--FILE|--F|--f","xxxx");
    if(XHOST.vflag_control.flag("VFILES")) aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VFILES"),vf,",");
    // make commands
    vector<string> vcmds;vcmds.clear();
    for(uint i=0;i<vf.size();i++) if(!aurostd::substring2bool(vf.at(i),".bz2")) vcmds.push_back(XHOST.command("bzip2")+" -9vf "+vf.at(i));
    AFLOW_PTHREADS::Check_Threads_WrapperNP(argv,vcmds.size(),VERBOSE); // check treads from NP
    // for(uint i=0;i<vcmds.size();i++) cout << vcmds.at(i) << endl;
    aurostd::multithread_execute(vcmds,AFLOW_PTHREADS::MAX_PTHREADS,VERBOSE);
    cout << endl;
    // done
    return TRUE;
  }
} // namespace AFLOW_PTHREADS

// ***************************************************************************
// AFLOW_PTHREADS::MULTI_bunzip2
// **************************************************************************
namespace AFLOW_PTHREADS {
  bool MULTI_bunzip2(vector<string> argv) {
    // aflow --multibunzip2 [--np XX | npmax | nothing] --F[ILE] file1 file2 file3 ....
    // load files
    bool VERBOSE=TRUE;
    vector<string> vf;
    // [OBSOLETE]    if(aurostd::args2flag(argv,"--FILE|--F|--f")) vf=aurostd::args2vectorstring(argv,"--FILE|--F|--f","xxxx");
    if(XHOST.vflag_control.flag("VFILES")) aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VFILES"),vf,",");
    // make commands
    vector<string> vcmds;vcmds.clear();
    for(uint i=0;i<vf.size();i++) if(aurostd::substring2bool(vf.at(i),".bz2")) vcmds.push_back(XHOST.command("bzip2")+" -9dvf "+vf.at(i));
    AFLOW_PTHREADS::Check_Threads_WrapperNP(argv,vcmds.size(),VERBOSE); // check treads from NP
    // for(uint i=0;i<vcmds.size();i++) cout << vcmds.at(i) << endl;
    aurostd::multithread_execute(vcmds,AFLOW_PTHREADS::MAX_PTHREADS,VERBOSE);
    cout << endl;
    // done
    return TRUE;
  }
} // namespace AFLOW_PTHREADS

// ***************************************************************************
// AFLOW_PTHREADS::MULTI_gzip
// **************************************************************************
namespace AFLOW_PTHREADS {
  bool MULTI_gzip(vector<string> argv) {
    // aflow --multigzip [--np XX | npmax | nothing] --F[ILE] file1 file2 file3 ....
    // load files
    bool VERBOSE=TRUE;
    vector<string> vf;
    // [OBSOLETE]    if(aurostd::args2flag(argv,"--FILE|--F|--f")) vf=aurostd::args2vectorstring(argv,"--FILE|--F|--f","xxxx");
    if(XHOST.vflag_control.flag("VFILES")) aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VFILES"),vf,",");
    // make commands
    vector<string> vcmds;vcmds.clear();
    for(uint i=0;i<vf.size();i++) if(!aurostd::substring2bool(vf.at(i),".gz")) vcmds.push_back(XHOST.command("gzip")+" -9vf "+vf.at(i));
    AFLOW_PTHREADS::Check_Threads_WrapperNP(argv,vcmds.size(),VERBOSE); // check treads from NP
    // for(uint i=0;i<vcmds.size();i++) cout << vcmds.at(i) << endl;
    aurostd::multithread_execute(vcmds,AFLOW_PTHREADS::MAX_PTHREADS,VERBOSE);
    cout << endl;
    // done
    return TRUE;
  }
} // namespace AFLOW_PTHREADS

// ***************************************************************************
// AFLOW_PTHREADS::MULTI_gunzip
// **************************************************************************
namespace AFLOW_PTHREADS {
  bool MULTI_gunzip(vector<string> argv) {
    // aflow --multigunzip [--np XX | npmax | nothing] --F[ILE] file1 file2 file3 ....
    // load files
    bool VERBOSE=TRUE;
    vector<string> vf;
    // [OBSOLETE]    if(aurostd::args2flag(argv,"--FILE|--F|--f")) vf=aurostd::args2vectorstring(argv,"--FILE|--F|--f","xxxx");
    if(XHOST.vflag_control.flag("VFILES")) aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VFILES"),vf,",");
    // make commands
    vector<string> vcmds;vcmds.clear();
    for(uint i=0;i<vf.size();i++) if(aurostd::substring2bool(vf.at(i),".gz")) vcmds.push_back(XHOST.command("gzip")+" -9dvf "+vf.at(i));
    AFLOW_PTHREADS::Check_Threads_WrapperNP(argv,vcmds.size(),VERBOSE); // check treads from NP
    // for(uint i=0;i<vcmds.size();i++) cout << vcmds.at(i) << endl;
    aurostd::multithread_execute(vcmds,AFLOW_PTHREADS::MAX_PTHREADS,VERBOSE);
    cout << endl;
    // done
    return TRUE;
  }
} // namespace AFLOW_PTHREADS


// **************************************************************************
// NOT MULTITHREAD BUT GOOD ENOUGH....

// ***************************************************************************
// sflow::KILL
// ***************************************************************************
namespace sflow {
  void KILL(string options) {
    vector<int> jobs;
    aurostd::StringCommasColumsVectorInt(options,jobs);
    vector<string> vcmds;

    stringstream aus;
    for(uint i=0;i<jobs.size();i++) {
      if(XHOST.is_command("kill")) {
	aus.clear();aus.str(std::string());
	aus << XHOST.command("kill") << " -9 " << jobs.at(i) << endl;  // command = kill
	vcmds.push_back(aus.str());
      }
    }
    for(uint i=0;i<vcmds.size();i++) {
      cout << "EXECUTING: " << vcmds.at(i);// << endl;
      aurostd::execute(vcmds.at(i));
    }
  }
} // namespace sflow

// ***************************************************************************
// sflow::JUST
// ***************************************************************************
namespace sflow {
  void JUST(string options,istream& input,string mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "sflow::JUST: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
   
    string strout="";
    ostringstream osstemp;
    osstemp << input.rdbuf();
    strout=osstemp.str();
   
    if(mode=="JUSTAFTER" || mode=="AFTER") {
      if(tokens.size()!=1 ) {
	init::ErrorOption(cout,options,"sflow::JUST","aflow --justafter=string < something");
	exit(0);
      }
      vector<string> vline;
      aurostd::string2vectorstring(strout,vline);
      bool found=FALSE;
      for(uint i=0;i<vline.size();i++) {
	if(found) cout << vline.at(i) << endl;
	if(!found) found=aurostd::substring2bool(vline.at(i),tokens.at(0));
      }
    }
    if(mode=="JUSTBEFORE" || mode=="BEFORE") {
      if(tokens.size()!=1 ) {
	init::ErrorOption(cout,options,"sflow::JUST","aflow --justbefore=string < something");
	exit(0);
      }
      if(strout.find(tokens.at(0))==string::npos) cout << strout;
      strout=strout.substr(0,strout.find(tokens.at(0)));
      strout=strout.substr(0,strout.find_last_of("\n")+1);
      cout << strout;
    }

    if(mode=="JUSTBETWEEN" || mode=="BETWEEN") {
      if(tokens.size()>2) {
	init::ErrorOption(cout,options,"sflow::JUST","aflow --justbetween=string_start[,string_stop] < something");
	exit(0);
      }
 
      string strfind_from,strfind_to,avoid1_strfind_from,avoid1_strfind_to;
      if(tokens.size()==1) {
	strfind_from="START."+tokens.at(0);
	strfind_to="STOP."+tokens.at(0);
      }
      if(tokens.size()==2) {
	strfind_from=tokens.at(0);
	strfind_to=tokens.at(1);
      }
      avoid1_strfind_from="#"+strfind_from;avoid1_strfind_to="#"+strfind_to;

      vector<string> vstranalyze;
      uint istart=0,istop=0;
      aurostd::string2vectorstring(osstemp.str(),vstranalyze);
      for(uint i=0;i<vstranalyze.size();i++) {
	if(aurostd::substring2bool(vstranalyze.at(i),strfind_from)) istart=i;
	if(aurostd::substring2bool(vstranalyze.at(i),strfind_to)) istop=i;
      }
      if(LDEBUG) cerr << "sflow::JUST: istart=" << istart << endl;
      if(LDEBUG) cerr << "sflow::JUST: istop=" << istop << endl;
      
      for(uint i=0;i<vstranalyze.size();i++)
	if(i>istart && i<istop) 
	  cout << vstranalyze.at(i) << endl;
    }
    
    if(LDEBUG) cerr << "sflow::JUST: END" << endl;
  }  
} // namespace sflow

// ***************************************************************************
// sflow::QDEL qdel bkill scancel
// ***************************************************************************
namespace sflow {
  void QDEL(string options,string cmd) {
    vector<int> jobs;
    aurostd::StringCommasColumsVectorInt(options,jobs);
    vector<string> vcmds;
    stringstream aus;
    for(uint i=0;i<jobs.size();i++) {
      aus.clear();aus.str(std::string());
      aus << cmd << " " << jobs.at(i) << endl;  // cmd = qdel OR bkill
      vcmds.push_back(aus.str());
    }
    for(uint i=0;i<vcmds.size();i++) {
      cout << "EXECUTING: " << vcmds.at(i);// << endl;
      aurostd::execute(vcmds.at(i));
    }
  }
} // namespace sflow

namespace sflow {
  void QDEL(string options) {
    vector<int> jobs;
    aurostd::StringCommasColumsVectorInt(options,jobs);
    vector<string> vcmds;
    // if(aurostd::args2flag(argv,"--scancel")) {XHOST.is_command("qdel")=FALSE;XHOST.is_command("bkill")=FALSE;} // force
    // if(aurostd::args2flag(argv,"--bkill")) {XHOST.is_command("scancel")=FALSE;XHOST.is_command("bkill")=FALSE;} // force
    // if(XHOST.is_command("qdel")) {XHOST.is_command("scancel")=FALSE;XHOST.is_command("bkill")=FALSE;} // some priority

    stringstream aus;
    for(uint i=0;i<jobs.size();i++) {
      if(XHOST.is_command("scancel")) {
	aus.clear();aus.str(std::string());
	aus << XHOST.command("scancel") << " " << jobs.at(i) << endl;  // cmd = qdel OR bkill
	vcmds.push_back(aus.str());
      } else {
	if(XHOST.is_command("qdel")) {
	  aus.clear();aus.str(std::string());
	  aus << XHOST.command("qdel") << " " << jobs.at(i) << endl;  // cmd = qdel OR bkill
	  vcmds.push_back(aus.str());
	} else {
	  if(XHOST.is_command("bkill")) {
	    aus.clear();aus.str(std::string());
	    aus << XHOST.command("bkill") << " " << jobs.at(i) << endl;  // cmd = qdel OR bkill
	    vcmds.push_back(aus.str());
	  }
	}
      }
    }
    for(uint i=0;i<vcmds.size();i++) {
      cout << "EXECUTING: " << vcmds.at(i);// << endl;
      aurostd::execute(vcmds.at(i));
    }
  }
} // namespace sflow

// ***************************************************************************
// sflow::QSUB qsub bsub sbatch
// ***************************************************************************
namespace sflow {
  void QSUB(string options,string cmd) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "sflow::QSUB: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(tokens.size()!=2) {
      init::ErrorOption(cout,options,"sflow::QSUB","aflow --qsub=N,file");
      exit(0);
    }
 
    vector<string> vcmds;
    stringstream aus;
    for(uint i=0;i<(uint) aurostd::string2utype<int>(tokens.at(0));i++) {
      aus.clear();aus.str(std::string());
      aus << cmd << " " << tokens.at(1) << endl;  // cmd = sbatch OR bsub <
      vcmds.push_back(aus.str());
    }
    for(uint i=0;i<vcmds.size();i++) {
      cout << "EXECUTING: " << vcmds.at(i);// << endl;
      aurostd::execute(vcmds.at(i));
    }
    if(LDEBUG) cerr << "sflow::QSUB: END" << endl;
  }
} // namespace sflow

namespace sflow {
  void QSUB(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "sflow::QSUB: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(tokens.size()!=2) {
      init::ErrorOption(cout,options,"sflow::QSUB","aflow --qsub=N,file");
      exit(0);
    }
 
    vector<string> vcmds;
    stringstream aus;
    for(uint i=0;i<(uint) aurostd::string2utype<int>(tokens.at(0));i++) {
      if(XHOST.is_command("sbatch")) {
	aus.clear();aus.str(std::string());
	if(XHOST.is_MACHINE_FULTON_MARYLOU) aus << XHOST.command("sbatch") << " -C beta " << tokens.at(1) << endl;  // cmd = sbatch
	else aus << XHOST.command("sbatch") << "  " << tokens.at(1) << endl;  // cmd = sbatch
	vcmds.push_back(aus.str());
      } else {
	if(XHOST.is_command("bsub")) {
	  aus.clear();aus.str(std::string());
	  aus << XHOST.command("bsub") << " <" << " " << tokens.at(1) << endl;  // cmd = bsub <
	  vcmds.push_back(aus.str());
	} else {
	  if(XHOST.is_command("qsub")) {
	    aus.clear();aus.str(std::string());
	    aus << XHOST.command("qsub") << " " << tokens.at(1) << endl;  // cmd = qsub
	    vcmds.push_back(aus.str());
	  }
	}
      }
    }
    for(uint i=0;i<vcmds.size();i++) {
      cout << "EXECUTING: " << vcmds.at(i);// << endl;
      aurostd::execute(vcmds.at(i));
    }
    if(LDEBUG) cerr << "sflow::QSUB: END" << endl;
  }
} // namespace sflow

// **************************************************************************

#endif  // _PTHREADS_IMPLEMENTATIONS_

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************
