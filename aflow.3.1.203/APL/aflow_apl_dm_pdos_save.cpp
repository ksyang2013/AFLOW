// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow PINKU NATH - Duke University 2014-2017             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu
//main intension of this cpp file is to save dynamical matrices and phonon density of states from other sub-directories

#include "aflow_apl.h"

//CO - START
// Some parts are written within the C++0x support in GCC, expecially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif
//CO - END

namespace apl {
// ***************************************************************************************
DM_PDOS_save::DM_PDOS_save(IPhononCalculator& pc, Logger& l) : _pc(pc), _logger(l) {
  nBranches = _pc.getNumberOfBranches();
  _kpoints.clear();
  DM.clear();
  _weights.clear();
  _rlattice.clear();
  _klattice.clear();
}
DM_PDOS_save::~DM_PDOS_save() { clear(); }
void DM_PDOS_save::clear() {
  string ex = _EOSdir;
  string command;
  if (aurostd::FileExist("apl.xml")) {
    command = string("mv ") + string("apl.xml") + string("  ") +
              string(dirfrefix) + string("/") + string("apl.xml") + string(ex) +
              string(" ") + string("2> /dev/null");
    aurostd::execute2string(command);
  }
  if (aurostd::FileExist("DYNMAT")) {
    command = string("mv ") + string("DYNMAT") + string(" ") + string(dirfrefix) + string("/") + string("DYNMAT.") + string(ex) +
              string(" ") + string("2> /dev/null");
    aurostd::execute2string(command);
  }

  //clear all variables
  _kpoints.clear();
  DM.clear();
  _weights.clear();
  _rlattice.clear();
  _klattice.clear();
}
// ***************************************************************************************
void DM_PDOS_save::createMPmesh(int na, int nb, int nc, const xstructure& xs) {
  _logger << "Creating qmesh "
          << "\n";

  xstructure xstr(xs);

  // Setup both lattices
  xstr.FixLattices();
  _rlattice = xstr.lattice;
  _klattice = ReciprocalLattice(_rlattice);

  _kpoints.clear();
  xvector<double> kpoint(3);
  for (int s = 1; s <= nc; s++)
    for (int r = 1; r <= nb; r++)
      for (int p = 1; p <= na; p++) {
        kpoint(1) = (2.0 * p - na - 1) / (2.0 * na);
        kpoint(2) = (2.0 * r - nb - 1) / (2.0 * nb);
        kpoint(3) = (2.0 * s - nc - 1) / (2.0 * nc);
        // Transform to cartesian coordinate
        kpoint = trasp(_klattice) * kpoint;
        _kpoints.push_back(kpoint);
        _weights.push_back(1.0);
      }

  string ex = _EOSdir;
  string dmfile_prefix = string(dirfrefix) + string("/") + string("DM_") + string("mesh");
  string dmfile = string(dmfile_prefix) + string(".") + string(ex);
  string kfile_prefix = string(dirfrefix) + string("/") + string("k_") + string("mesh");
  string kfile = string(kfile_prefix) + string(".") + string(ex);
  string weightfile = string(kfile_prefix) + string("_weight") + string(".") + string(ex);

  //saving DM
  DM.clear();
  xmatrix<xcomplex<double> > tDM(nBranches, nBranches, 1, 1);
  DM.resize(_kpoints.size(), tDM);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
  int qpointsPerCPU = _kpoints.size() / ncpus;

  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating dynamical matrices in k-mesh");
  else
    _logger.initProgressBar("Calculating dynamical matrices in k-mesh (" + stringify(ncpus) + " threads)");
  _logger << "Number of kpoints:: " << (int)_kpoints.size() << apl::endl;

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread*> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _kpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _kpoints.size())))
      endIndex = _kpoints.size();
    threads.push_back(new std::thread(&DM_PDOS_save::get_dm, this, startIndex, endIndex));
  }
  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
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
  DM.clear();
}
// ***************************************************************************************
//create phonon-dispersion path
void DM_PDOS_save::create_pdispath(std::vector<xvector<double> >& qpoints) {
  if (qpoints.size() == 0) {
    _logger << apl::error << "qpoints.size()==0 in create_pdispath() " << apl::endl;
    exit(0);
  }

  _logger << "Preparing to create dynamical matrices for high-symmetry qpoint" << apl::endl;

  _kpoints.clear();
  _kpoints = qpoints;
  string ex = _EOSdir;
  string dmfile_prefix = string(dirfrefix) + string("/") + string("DM_") + string("path");
  string dmfile = string(dmfile_prefix) + string(".") + string(ex);

  DM.clear();
  xmatrix<xcomplex<double> > tDM(nBranches, nBranches, 1, 1);
  DM.resize(_kpoints.size(), tDM);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
  int qpointsPerCPU = _kpoints.size() / ncpus;

  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating dynamical matrices along path");
  else
    _logger.initProgressBar("Calculating dynamical matrices along path (" + stringify(ncpus) + " threads)");
  _logger << "Number of kpoints:: " << (int)_kpoints.size() << apl::endl;

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread*> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _kpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _kpoints.size())))
      endIndex = _kpoints.size();
    threads.push_back(new std::thread(&DM_PDOS_save::get_dm, this, startIndex, endIndex));
  }
  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
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
//check whether running directory is Gruneisen directory
bool DM_PDOS_save::check_GP() {
  string dir = _EOSdir;
  string tag = NumberToString<double>(GP_VOL_DISTORTION);
  string dir1 = string(_GP_AFLOW_DIRECTORY_PREFIX_) + string("M") + string(tag);
  string dir2 = string(_GP_AFLOW_DIRECTORY_PREFIX_) + string("P") + string(tag);
  if ((dir == dir1) || (dir == dir2))
    return true;
  else
    return false;
}
// ***************************************************************************************
string
DM_PDOS_save::getdir_name(string path) {
  //listing directory names
  string bashcommand = "ls -d " + string(_GP_AFLOW_DIRECTORY_PREFIX_) + string("*") + string("  ") + string("2> /dev/null");
  vector<string> s = directory_list(bashcommand);
  string ex = "";
  for (uint i = 0; i < s.size(); i++) {
    std::size_t str1 = path.find(s[i]);
    if (str1 != std::string::npos) { ex = s[i]; }
  }
  _EOSdir = ex;

  if (!aurostd::FileExist(dirfrefix)) {
    bashcommand = string("mkdir -p ") + string(dirfrefix) + string("  ") + string("2> /dev/null");
    ;
    aurostd::execute2string(bashcommand);
  }

  _logger << "Directory name found to be to " << _EOSdir << "\n";

  return _EOSdir;
}
// ***************************************************************************************
//get running directory name
vector<string>
DM_PDOS_save::directory_list(const string path) {
  vector<string> s;
  FILE* in;
  char buff[512];
  if (!(in = popen(path.c_str(), "r"))) {
    _logger << apl::error << "bash command can't not able to exacute" << apl::endl;
    exit(0);
  }
  while (fgets(buff, sizeof(buff), in) != NULL) {
    stringstream AFLOW;
    string dirlist;
    AFLOW << buff;
    AFLOW >> dirlist;
    s.push_back(dirlist);
  }
  return s;
}
// ***************************************************************************************
//Number to String conver function
template <typename T>
string DM_PDOS_save::NumberToString(T Number) {
  ostringstream ss;
  ss << Number;
  return ss.str();
}
// ***************************************************************************************
//this function returns dynamical matricx
void DM_PDOS_save::get_dm(int startIndex, int endIndex) {
  xmatrix<xcomplex<double> > tDM(nBranches, nBranches, 1, 1);
  for (int kpoint = startIndex; kpoint < endIndex; kpoint++) {
    xvector<double> qpoint = _kpoints[kpoint];
    tDM = _pc.getDynamicMatrix(qpoint);
    DM[kpoint] = tDM;
  }
}
// ***************************************************************************************
//write dynamica matrices to a file in a binary format to save space and reading is fast
void DM_PDOS_save::write_dm(string file) {
  _logger << "writing dynamical matric to " << file << apl::endl;

  ofstream dmout;
  dmout.open(file.c_str(), ios_base::binary);
  if (!dmout.is_open()) {
    _logger << apl::error << file << " not able to open " << apl::endl;
    return;
  }

  for (uint I = 0; I < DM.size(); I++) {
    for (uint J = 1; J <= nBranches; J++) {
      dmout.write((char*)(&DM[I][J][1]), nBranches * sizeof(xcomplex<double>));
    }
  }
  dmout.clear();
  dmout.close();
}
// ***************************************************************************************
//write q-point wright to a file
void DM_PDOS_save::write_kpoints(string file) {
  //CO - START
  //ofstream kout (file.c_str());
  stringstream kout;

  //if (kout.is_open()){_logger<<apl::error<<file<<" not able to open"<<std::endl; return;}
  //CO - END

  kout << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  kout << setprecision(SetPrecision);
  for (uint kpoint = 0; kpoint < _kpoints.size(); kpoint++) {
    for (int i = 1; i <= _kpoints[i].rows; i++) {
      kout << setw(SetPrecision + 10) << _kpoints[kpoint][i];
    }
    kout << "\n";
  }
  //CO - START
  aurostd::stringstream2file(kout, file);
  if (!aurostd::FileExist(file)) {
    _logger << apl::error << file << " not able to open" << std::endl;
    return;
  }
  //kout.clear();
  //kout.close();
  //CO - END
}
// ***************************************************************************************
//write q-point wright to a file
void DM_PDOS_save::write_weight(string file) {
  //weight factor
  //CO - START
  //ofstream kout1(file.c_str());
  stringstream kout1;
  //if (kout1.is_open()){_logger<<apl::error <<" DM_PDOS_save; unable to create [kweights-file]  "<<apl::endl;exit(0);}
  //CO - START

  kout1 << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  kout1 << setprecision(SetPrecision);
  for (uint kpoint = 0; kpoint < _weights.size(); kpoint++) {
    kout1 << setw(SetPrecision + 10) << _weights[kpoint] << "\n";
  }
  //CO - START
  aurostd::stringstream2file(kout1, file);
  if (!aurostd::FileExist(file)) {
    _logger << apl::error << " DM_PDOS_save; unable to create [kweights-file]  " << apl::endl;
    exit(0);
  }
  //kout1.clear();
  //kout1.close();
  //CO - END
}
// ***************************************************************************************
//check qpoint path are same in all other sub-directories
void check_consistency_aflow_apl::check_consistency_in_aflow(string& current_running_aflowin) {
  if (_RUNNING_DIR.find(_GP_AFLOW_DIRECTORY_PREFIX_) == std::string::npos) return;
  string file_aflowin = _RUNNING_DIR + "/" + _AFLOWIN_;
  if (!aurostd::FileExist(file_aflowin)) return;
  string oldDCUSERPATH = "";
  string line;
  bool present = false;
  //ifstream read_path (file_aflowin.c_str()); //CO
  //if (!read_path.is_open()){_logger<<apl::error<<file_aflowin<<" not able to open "<<apl::endl; exit(0);} //CO

  //CO - START
  vector<string> vlines;
  uint line_count = 0;
  aurostd::efile2vectorstring(file_aflowin, vlines);
  //ifstream read_path (file_aflowin.c_str());
  //if (!read_path.is_open()){_logger<<apl::error<<file_aflowin<<" not able to open "<<apl::endl; exit(0);}
  if (!vlines.size()) {
    _logger << apl::error << file_aflowin << " not able to open " << apl::endl;
    exit(0);
  }

  //while ( getline (read_path,line) )
  while (line_count < vlines.size()) {
    line = vlines[line_count++];
    if (line == "") continue;
    if (line[0] == '#') continue;
    if (line.find("DCUSERPATH") != std::string::npos) {
      oldDCUSERPATH = line;
      present = true;
      break;
    }
  }  //read_path.close();
  if (!present) return;
  //CO - END

  string DCUSERPATH = "";
  string DIR_ = _RUNNING_DIR;
  DIR_.erase(DIR_.end() - _RUNNING_DIR.length(), DIR_.end());
  if (DIR_ == "") DIR_ = "./";
  file_aflowin = DIR_ + "/" + _AFLOWIN_;
  //CO - START
  line_count = 0;
  vlines.clear();
  //read_path.open(file_aflowin.c_str());
  aurostd::efile2vectorstring(file_aflowin, vlines);
  //if (!read_path.is_open()){ _logger<<apl::error<<file_aflowin <<" not able to open "<<apl::endl; exit(0);}
  if (!vlines.size()) {
    _logger << apl::error << file_aflowin << " not able to open " << apl::endl;
    exit(0);
  }
  //while ( getline (read_path,line) )
  while (line_count < vlines.size()) {
    line = vlines[line_count++];
    if (line == "") continue;
    if (line.find("DCUSERPATH") != std::string::npos) {
      DCUSERPATH = line;
      break;
    }
  }  //read_path.close();
  if (DCUSERPATH == "") return;

  //ofstream myaflow (tmpaflow.c_str());
  //if (!myaflow.is_open()){_logger<<apl::error <<tmpaflow <<" unable to open"<<std::endl; exit(0);}
  stringstream myaflow;

  file_aflowin = _RUNNING_DIR + "/" + _AFLOWIN_;
  //read_path.open(file_aflowin.c_str());
  line_count = 0;
  vlines.clear();
  aurostd::efile2vectorstring(file_aflowin, vlines);

  //if (!read_path.is_open()){ _logger<<apl::error << file_aflowin <<" not able to open "<<apl::endl; exit(0);}
  if (!vlines.size()) {
    _logger << apl::error << file_aflowin << " not able to open " << apl::endl;
    exit(0);
  }

  //while ( getline (read_path,line) )
  while (line_count < vlines.size()) {
    line = vlines[line_count++];
    if (line.find("DCUSERPATH") != std::string::npos) {
      myaflow << DCUSERPATH << "\n";
    } else {
      myaflow << line << "\n";
    }

  }  //read_path.close();myaflow.close();
  string tmpaflow = _RUNNING_DIR + "/tmp_aflow.in";
  aurostd::stringstream2file(myaflow, tmpaflow);
  if (aurostd::FileExist(tmpaflow)) {
    _logger << apl::error << tmpaflow << " unable to open" << std::endl;
    exit(0);
  }

  if (oldDCUSERPATH == DCUSERPATH) {
    //string command="rm "+tmpaflow+"  "+string("2> /dev/null");
    //aurostd::execute2string(command);
    aurostd::RemoveFile(tmpaflow);
    return;
  }
  file_aflowin = _RUNNING_DIR + "/" + _AFLOWIN_;
  aurostd::MoveFile(tmpaflow, file_aflowin);
  //string command="mv "+tmpaflow+"  "+file_aflowin+string(" ")+string("2> /dev/null");
  //aurostd::execute2string(command);
  aurostd::ChmodFile("777", file_aflowin);
  //command="chmod 777 "+file_aflowin+string(" ")+string("2> /dev/null");
  //aurostd::execute2string(command);
  //_logger<<"\n PINKU aflow change has been made \n";
  replace(current_running_aflowin, oldDCUSERPATH, DCUSERPATH);
  //CO - END
}
// ***************************************************************************************
//if qpoint path is not same with main aflow.in file then change with an appropriate path
bool check_consistency_aflow_apl::replace(std::string& str, const std::string& from, const std::string& to) {
  size_t start_pos = str.find(from);
  if (start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}
// ***************************************************************************************
}
