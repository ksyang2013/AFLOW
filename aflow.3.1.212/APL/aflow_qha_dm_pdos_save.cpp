// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu
//main intension of this cpp file is to save dynamical matrices and phonon density of states from other sub-directories

#include "aflow_apl.h"
#include <iterator>

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

namespace apl {
  // ***************************************************************************************
  QHAsubdirectoryData::QHAsubdirectoryData( IPhononCalculator& pc, Logger& l):
    _pc(pc),_logger(l)
  {
    _gp_vol_distortion=0.0;
    _sc_vol_distortion=0.0;
    _kpoints.clear();
    _dm.clear();
    _weights.clear();
    _rlattice.clear();
    _klattice.clear();
    _running_dir="";
    _tmpdirfrefix="";
    _devnull=" >/dev/null 2>&1";
    _nBranches=_pc.getNumberOfBranches();
  }
  // ***************************************************************************************
  QHAsubdirectoryData::~QHAsubdirectoryData()
  {
    clear();
  }
  // ***************************************************************************************
  void QHAsubdirectoryData::clear()
  {
    if(_running_dir=="") return;

    _logger<<"Moving apl.xml and DYNMAT files to the dir "<<_tmpdirfrefix<< apl::endl;
    string ex=_running_dir;
    string command="";
    if(aurostd::FileExist("apl.xml")){
      command =string("mv ")+string("apl.xml")+string("  ")+
	string(_tmpdirfrefix)+string("/")+string("apl.xml.")+string(ex)+string(" ")+_devnull;
      aurostd::execute2string(command);
    }
    if(aurostd::FileExist("DYNMAT")){
      command =string("mv ")+string("DYNMAT")+string(" ")+string(_tmpdirfrefix)+string("/")+string("DYNMAT.")+string(ex)+
	string(" ")+_devnull;
      aurostd::execute2string(command);}
    //clear all variables
    _kpoints.clear();
    _dm.clear();
    _weights.clear();
    _rlattice.clear();
    _klattice.clear();
    _gp_vol_distortion=0.0;
    _sc_vol_distortion=0.0;
    _running_dir="";
    _tmpdirfrefix="";
  }
  // ***************************************************************************************
  void QHAsubdirectoryData::createMPmesh(int na, int nb, int nc,const xstructure& xs)
  {
    _logger<<"Creating uniform mesh in the dir "<<_running_dir<< apl::endl;

    xstructure xstr(xs);

    // Setup both lattices
    xstr.FixLattices();
    _rlattice = xstr.lattice;
    _klattice = ReciprocalLattice(_rlattice);

    _kpoints.clear();
    xvector<double> kpoint(3);
    for(int s = 1; s <=nc; s++)
      for(int r = 1; r <= nb; r++)
	for(int p = 1; p <= na; p++)
	  {
	    kpoint(1) = ( 2.0 * p - na - 1 ) / ( 2.0 * na );
	    kpoint(2) = ( 2.0 * r - nb - 1 ) / ( 2.0 * nb );
	    kpoint(3) = ( 2.0 * s - nc - 1 ) / ( 2.0 * nc );
	    // Transform to cartesian coordinate
	    kpoint = trasp(_klattice) * kpoint;
	    _kpoints.push_back(kpoint);
	    _weights.push_back(1.0);
	  }
    }
  // ***************************************************************************************
  void QHAsubdirectoryData::create_dm()
  {
    string ex=_running_dir;
    string dmfile_prefix= string(_tmpdirfrefix)+string("/")+string("DM_")+string("mesh");
    string dmfile       = string(dmfile_prefix)+string(".")+string(ex);
    string kfile_prefix = string(_tmpdirfrefix)+string("/")+string("k_") +string("mesh");
    string kfile        =string(kfile_prefix)+string(".")+string(ex);
    string weightfile   =string(kfile_prefix)+string("_weight")+string(".")+string(ex);
    
   _logger <<"Creating dynamical matrix "<< dmfile << apl::endl;

    //saving DM 
    _dm.clear();
    xmatrix<xcomplex<double> > tDM(_nBranches,_nBranches,1,1);
    _dm.resize(_kpoints.size(), tDM);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    // Get the number of CPUS
    int ncpus = sysconf(_SC_NPROCESSORS_ONLN);// AFLOW_MachineNCPUs;
    if(ncpus<1) ncpus=1;
//    int qpointsPerCPU = _kpoints.size() / ncpus;  OBSOLETE ME180801

    // Show info 
    if( ncpus == 1 )
      _logger.initProgressBar("Calculating dynamical matrices in k-mesh");
    else
      _logger.initProgressBar("Calculating dynamical matrices in k-mesh (" + stringify(ncpus) + " threads)");
    _logger<<"Number of kpoints:: "<< (int)_kpoints.size() <<apl::endl;

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _kpoints.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&QHAsubdirectoryData::get_dm,this,startIndex,endIndex) );
    }

/* OBSOLETE ME180801
    for(int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = icpu * qpointsPerCPU;
      endIndex = startIndex + qpointsPerCPU;
      if( ( (uint)endIndex > _kpoints.size() ) ||
          ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _kpoints.size() ) ) )
        endIndex = _kpoints.size();
      threads.push_back( new std::thread(&QHAsubdirectoryData::get_dm,this,startIndex,endIndex) );
    }
*/

    // Wait to finish all threads here!
    for(uint i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }
    // Done
    _logger.finishProgressBar();

#else
    get_dm(0, (int)_kpoints.size());
#endif
    write_dm(dmfile);
    //write_kpoints(kfile);
    //write_weight(weightfile);
    _dm.clear();
  }
  // ***************************************************************************************
  //create phonon-dispersion path
  void QHAsubdirectoryData::create_pdispath(std::vector< xvector<double> > &qpoints)
  {
    if(qpoints.size()==0){_logger<<apl::error<<"qpoints.size()==0 in create_pdispath() "<<apl::endl; exit(0);}

    _logger<<"Preparing to create dynamical matrices for high-symmetry qpoint"<<apl::endl;

    _kpoints.clear();
    _kpoints=qpoints;
    string ex=_running_dir;
    string dmfile_prefix= string(_tmpdirfrefix)+string("/")+string("DM_")+string("path");
    string dmfile       = string(dmfile_prefix)+string(".")+string(ex);

    _dm.clear();
    xmatrix<xcomplex<double> > tDM(_nBranches,_nBranches,1,1);
    _dm.resize(_kpoints.size(), tDM);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    // Get the number of CPUS
    int ncpus = sysconf(_SC_NPROCESSORS_ONLN);// AFLOW_MachineNCPUs;
    if(ncpus<1) ncpus=1;
//    int qpointsPerCPU = _kpoints.size() / ncpus;  ME180801

    // Show info 
    if( ncpus == 1 )
      _logger.initProgressBar("Calculating dynamical matrices along path");
    else
      _logger.initProgressBar("Calculating dynamical matrices along path (" + stringify(ncpus) + " threads)");
    _logger<<"Number of kpoints:: "<< (int)_kpoints.size() <<apl::endl;

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _kpoints.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&QHAsubdirectoryData::get_dm,this,startIndex,endIndex) );
    }

/* OBSOLETE ME180801
    for(int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = icpu * qpointsPerCPU;
      endIndex = startIndex + qpointsPerCPU;
      if( ( (uint)endIndex > _kpoints.size() ) ||
          ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _kpoints.size() ) ) )
        endIndex = _kpoints.size();
      threads.push_back( new std::thread(&QHAsubdirectoryData::get_dm,this,startIndex,endIndex) );
    }
*/

    // Wait to finish all threads here!
    for(uint i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }
    // Done
    _logger.finishProgressBar();
#else
    get_dm(0, (int)_kpoints.size());
#endif

    write_dm(dmfile);
  }
  // ***************************************************************************************
  void QHAsubdirectoryData::setdir_prefix(const string s)
  {
    _tmpdirfrefix=s;
  }
  // ***************************************************************************************
  void QHAsubdirectoryData::set_gp_vol_distortion(const double d)
  {
    _gp_vol_distortion=d;
  }
  // ***************************************************************************************
  void QHAsubdirectoryData::set_sc_vol_distortion(const double d)
  {
    _sc_vol_distortion=d;
  }
  // ***************************************************************************************
  //check whether running directory is Gruneisen directory
  bool
  QHAsubdirectoryData::check_GP()
  {
    string dir=_running_dir;
    string tag=NumberToString<double>(_gp_vol_distortion);
    string dir1=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("M")+string(tag);
    string dir2=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("P")+string(tag);
    if((dir==dir1)||(dir==dir2))return true;
    else return false;
  }
  // ***************************************************************************************
  //check whether running directory is SCQHA directory
  bool QHAsubdirectoryData::check_SCQHA()
  {
    string dir=_running_dir;
    string tag=NumberToString<double>(_sc_vol_distortion);
    string dir1=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("M")+string(tag);
    string dir2=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("P")+string(tag);
    if((dir==dir1)||(dir==dir2))return true;
    else return false;
  }
  // ***************************************************************************************
  string QHAsubdirectoryData::getdir_name(string path)
  {
    vector<string> vec_str=splitWdelimiter<string>(path, "/");
   _running_dir=vec_str[vec_str.size()-1];
   _logger <<"Running directory "<< _running_dir << apl::endl;

    if( !aurostd::FileExist(_tmpdirfrefix) )
      {
	string bashcommand=string("mkdir -p ")+string(_tmpdirfrefix)+string("  ")+ _devnull;
	aurostd::execute2string(bashcommand);
      }
    return _running_dir;
  }
  // ***************************************************************************************
  template<class T> 
  vector<T> QHAsubdirectoryData::splitWdelimiter(string s, string delimiter)
   {
   vector<T> lines;
   size_t pos = 0;
   std::string token;
while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    vector<T> vec=split<T>(token);
    for(uint i=0; i<vec.size(); i++)lines.push_back(vec[i]);
    s.erase(0, pos + delimiter.length());
}
    vector<T> vec1=split<T>(s) ;
    for(uint i=0; i<vec1.size(); i++)lines.push_back(vec1[i]);

return lines;
}
  // ***************************************************************************************
  //get running directory name
  vector<string>
  QHAsubdirectoryData::directory_list(const string path)
  {
    vector<string> s;
    FILE *in;
    char buff[512];
    if(!(in = popen(path.c_str(), "r"))){_logger<<apl::error<<"bash command can't not able to exacute"<<apl::endl; exit(0);}
    while(fgets(buff, sizeof(buff), in)!=NULL)
      {stringstream AFLOW;
	string dirlist;
	AFLOW << buff;
	AFLOW >> dirlist;
	s.push_back(dirlist);}
    return s;
  }
  // ***************************************************************************************
template<typename T> 
   std::vector<T> QHAsubdirectoryData::split(const std::string& line){
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
}
  // ***************************************************************************************
  //Number to String conver function
  template <typename T>
  string QHAsubdirectoryData:: NumberToString ( T Number ){
    ostringstream ss;
    ss << Number;
    return ss.str();
  }
  // ***************************************************************************************
  //this function returns dynamical matricx
  void  QHAsubdirectoryData::get_dm(int startIndex, int endIndex)
  {
    xmatrix<xcomplex<double> > tDM(_nBranches,_nBranches,1,1);
    for(int kpoint=startIndex; kpoint<endIndex; kpoint++){
      xvector<double> qpoint=_kpoints[kpoint];
      tDM = _pc.getDynamicalMatrix(qpoint);
      _dm[kpoint]=tDM;}
  }
  // ***************************************************************************************
  //write dynamica matrices to a file in a binary format to save space and reading is fast
  void QHAsubdirectoryData::write_dm(string file)
  {
    _logger <<"Writing dynamical matric to file : "<<file<<apl::endl;

    ofstream dmout;
    dmout.open(file.c_str(),ios_base::binary);
    if (!dmout.is_open()){_logger<<apl::error<<file<<" not able to open "<<apl::endl; return;}

    for(uint I=0;I<_dm.size();I++){
      for(uint J=1;J<=_nBranches;J++){
	dmout.write( (char *)(&_dm[I][J][1]), _nBranches*sizeof(xcomplex<double>) );
      }}
    dmout.clear();
    dmout.close();
  }
  // ***************************************************************************************
  //write q-point wright to a file
  void QHAsubdirectoryData::write_kpoints(string file)
  {
    ofstream kout (file.c_str());

    if (kout.is_open()){_logger<<apl::error<<file<<" not able to open"<<std::endl; return;}

    kout << std::setiosflags(std::ios::fixed|std::ios::showpoint|std::ios::right);
    kout<< setprecision(Set_QHA_Precision);
    for(uint kpoint=0; kpoint<_kpoints.size(); kpoint++){
      for(int i=1; i<=_kpoints[i].rows; i++){
	kout<<setw(Set_QHA_Precision+10)<<_kpoints[kpoint][i];}
      kout<<"\n";}
    kout.clear();
    kout.close();
  }
  // ***************************************************************************************
  //write q-point wright to a file
  void QHAsubdirectoryData::write_weight(string file)
  {
    //weight factor
    ofstream kout1(file.c_str());
    if (kout1.is_open()){_logger<<apl::error <<" QHAsubdirectoryData:: unable to create [kweights-file]  "<<apl::endl;exit(0);}

    kout1 << std::setiosflags(std::ios::fixed|std::ios::showpoint|std::ios::right);
    kout1<< setprecision(Set_QHA_Precision);
    for(uint kpoint=0; kpoint<_weights.size(); kpoint++){
      kout1<<setw(Set_QHA_Precision+10)<<_weights[kpoint]<<"\n";}
    kout1.clear();
    kout1.close();
  }
  // ***************************************************************************************
}
