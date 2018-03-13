// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2011               *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_kesong_ipocc.cpp
// functions written by KESONG YANG
// 2010-2011: kesong.yang@gmail.com

// This file contains the routines to prepare partial occupation input files.

#ifndef _AFLOW_CONTRIB_KESONG_IPOCC_CPP
#define _AFLOW_CONTRIB_KESONG_IPOCC_CPP

#include "aflow_contrib_kesong.h"

using aurostd::StringSubst;
using std::setfill;
const int MaxNumberPOSCAR=999;

// ***************************************************************************
// pflow::POCC_INPUT(void)
// ***************************************************************************
namespace pflow {
  void POCC_INPUT(void) {
    _aflags aflags;aflags.Directory="./";
    ofstream oss;
    string aflowin, MESSAGE="pflow::POCC_INPUT ERROR";
    aflowin=string(aflags.Directory+_AFLOWIN_);
    if(!aurostd::FileExist(aflowin)) {
      cerr << MESSAGE << ": file not found " << aflowin << endl;
      exit(1);
    }
    POCC_GENERATE_INPUT(oss,aflags);
  }
}

// ***************************************************************************
// bool POCC_GENERATE_INPUT(ofstream &FileMESSAGE,_aflags &aflags)
// ***************************************************************************
bool POCC_GENERATE_INPUT(ofstream &FileMESSAGE,_aflags &aflags) {
  ostringstream aus;
  aus << "00000  MESSAGE running POCC_GENERATE_INPUT files " << Message(aflags,"user,host,time") << endl;
  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

  string aflowin=aflags.Directory+"/"+_AFLOWIN_;

  if(!aurostd::FileExist(aflowin)) {
      FileMESSAGE << "ERROR" << ": file not found " << aflowin << endl;
      return FALSE;
  }

  string AflowIn;aurostd::file2string(aflowin,AflowIn);
  stringstream sspoccSTRUCTURE;

  // CHECK FOR INSIDE STUFF
  if(!pocc::POCC_Already_Calculated_Input(AflowIn)) {
    stringstream input_file; input_file.clear();input_file.str("");
    stringstream input_file_aus; input_file_aus.clear();input_file_aus.str("");

    aurostd::ExtractToStringstreamEXPLICIT(AflowIn,sspoccSTRUCTURE, "[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE", "[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE");

    xstructure  xstr_pocc(sspoccSTRUCTURE, IOVASP_POSCAR);
    ofstream FileMESSAGE;
    ofstream file_aflowin;
    file_aflowin.open(aflowin.c_str(),ios_base::app);

    file_aflowin << AFLOWIN_SEPARATION_LINE << endl; 
    if(!aurostd::substring2bool(AflowIn,"[AFLOW_POCC]CALC")) {
        file_aflowin << "[AFLOW_POCC]CALC"<< endl;
    }

    vector<xstructure> vecgroupxstr_sorted = Partial2Supercell(xstr_pocc);
    stringstream ssxstr_sorted;
    stringstream ss;
    ssxstr_sorted.str("");

    int Num_calculated;
    int Num_xstr=vecgroupxstr_sorted.size();
    if(Num_xstr<MaxNumberPOSCAR) {
        Num_calculated=Num_xstr;
    }
    else {
        Num_calculated=MaxNumberPOSCAR;
    }

    for(int i=0;i<Num_calculated;i++) {
      ss.str("");
      ss << "POCC_" << setfill('0') << setw(2) <<(i+1);
      ssxstr_sorted << AFLOWIN_SEPARATION_LINE<< endl;
      ssxstr_sorted << "[VASP_POSCAR_MODE_EXPLICIT]START." <<ss.str() << endl;
      ssxstr_sorted << vecgroupxstr_sorted.at(i);
      ssxstr_sorted << "[VASP_POSCAR_MODE_EXPLICIT]STOP." <<ss.str() << endl;
      ssxstr_sorted << AFLOWIN_SEPARATION_LINE<< endl;
    }

    file_aflowin << ssxstr_sorted.str();
    file_aflowin.close();
  } 
  else {
    aus << "00000  MESSAGE POCC input file already created " << Message(aflags,"user,host,time") << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    return FALSE;
  }
  return TRUE;
}

#endif     
// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2011               *
// *                                                                         *
// ***************************************************************************

