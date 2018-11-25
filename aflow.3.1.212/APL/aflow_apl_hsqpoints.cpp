// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"
#include <iterator>

namespace apl
{
  // ***************************************************************************************
  PhononHSQpoints::PhononHSQpoints(Logger& l):_logger(l)
  {
    //clear all memories before using
    clear();
  }
  // ***************************************************************************************
  PhononHSQpoints::~PhononHSQpoints()
  {
     this->clear();
  }
  // ***************************************************************************************
  void PhononHSQpoints::clear()
  {
    _qpoints.clear();
    _path.clear();
    _path_segment.clear();
    _hs_kpoints.clear();
  }
  // ***************************************************************************************
  void PhononHSQpoints::read_qpointfile()
  {
    string file=DEFAULT_APL_HSKPTS_FILE;
    _logger << "Writing " << file << apl::endl;
    string command="";
    string bz2=string(file)+string(".bz2");
    bool t1=aurostd::FileExist(file);
    bool t2=aurostd::FileExist(bz2);
    if(!t1){
      if(t2){command=string("bzip2 -d ")+string(bz2);aurostd::execute2string(command);}
      else{
        _logger << apl::error << file<<" doesn't exist" << apl::endl; exit(0);
      }
    }

    string line;
    ifstream in (file.c_str());
    if (!in.is_open()){
      _logger << apl::error << file<<" doesn't exist" << apl::endl; exit(0);
    }
    while ( getline (in,line) )
      {
        if(line==""){
        }else if(line[0]=='#'){
        vector<string> vsr;
        apl::tokenize(line, vsr, string(" "));
         xvector<double> tmp(3,1);
         tmp[1]=atof(vsr[vsr.size()-3].c_str());
         tmp[2]=atof(vsr[vsr.size()-2].c_str());
         tmp[3]=atof(vsr[vsr.size()-1].c_str());
         _hs_kpoints.push_back(tmp);
        }else{
        vector<string> vsr;
        apl::tokenize(line, vsr, string(" "));
         if(vsr.size()!=5)
         {
            _logger << apl::error << file<<" format error " << apl::endl; exit(0);
         }
         xvector<double> tmp(3,1);
         tmp[1]=atof(vsr[0].c_str());
         tmp[2]=atof(vsr[1].c_str());
         tmp[3]=atof(vsr[2].c_str());
         _qpoints.push_back(tmp);
        _path_segment.push_back(atoi(vsr[3].c_str()));
        _path.push_back(atof(vsr[4].c_str()));
        }
      }
    in.close();
 }
  // ***************************************************************************************
    vector<xvector<double> > PhononHSQpoints::get_qpoints()
    {
      return _qpoints;
    }
  // ***************************************************************************************
    vector<xvector<double> > PhononHSQpoints::get_hs_kpoints()
    {
      return _hs_kpoints;
    }
  // ***************************************************************************************
    vector<double> PhononHSQpoints::get_path()
    {
      return _path;
    }
  // ***************************************************************************************
    vector<int> PhononHSQpoints::get_path_segment()
    {
      return _path_segment;
    }
  // ***************************************************************************************
}
