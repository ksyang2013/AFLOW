// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow PINKU NATH - Duke University 2014-2017             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"
//#include <numeric>
//#include <iostream>
//#include <vector>
//#include <functional>

#define _isnegative(a) (a < MIN_EIGEN_TRESHOLD) ? true : false

// Some parts are written within the C++0x support in GCC, expecially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

namespace apl {
// ***************************************************************************************
QuasiHarmonicGruneisen::QuasiHarmonicGruneisen(IPhononCalculator &pc,
                                               eos_calculate &runeos,
                                               Logger &l) : _pc(pc), _runeos(runeos), _logger(l) {
  clear();
  nBranches = _pc.getNumberOfBranches();
  _negative_freq_along_hsp = false;
}
QuasiHarmonicGruneisen::~QuasiHarmonicGruneisen() { this->clear(); }
void QuasiHarmonicGruneisen::clear() {
  ATOMIC_SPECIES.clear();
  ATOMIC_MASSES_AMU.clear();
  DMp.clear();
  DMm.clear();
  DM.clear();
  _eigenvectors.clear();
  _kpoints.clear();
  gp_path_test.clear();
  gp_mesh_test.clear();
  _rlattice.clear();
  _meshsize.clear();
  gpdir.clear();
  phdir.clear();
  eosdir.clear();
  gpvol.clear();
  eosvol.clear();
  _weights.clear();
  _GP_mesh.clear();
  GP_path.clear();
  _freqs_path.clear();
  _freqs_mesh.clear();
  _BrINDEXs.clear();
  _acoustic_freqs_mesh.clear();
  _acoustic_GP_mesh.clear();
}
// ***************************************************************************************
//initializing variables
bool QuasiHarmonicGruneisen::set_directories(string dir) {
  _logger << "Populationg variables to calculate Gruneisen parameters" << apl::endl;
  dirfrefix = dir;
  gpdir = _runeos.getGP_dir_names();
  gpvol = _runeos.getGP_volumes();
  phdir = _runeos.getPH_dir_names();
  eosdir = _runeos.getEOS_dir_names();
  eosvol = _runeos.getEOS_Volumes();
  ZERO_STATIC_DIR_INDEX = _runeos.getZERO_STATIC_DIR_INDEX();

  if ((gpdir.size() == 0) || (gpvol.size() == 0)) {
    _logger << apl::error << "directory sizes are zero remove LOCK and apl.xml and run again" << apl::endl;
    return false;
  }
  delta_V = 0.5 * ((gpvol[2] - gpvol[1]) + (gpvol[1] - gpvol[0]));
  V0 = gpvol[1];

  return true;
}
// ***************************************************************************************
//creating user defined q-mesh
void QuasiHarmonicGruneisen::createmesh(int na, int nb, int nc, const xstructure &xs) {
  _logger << "Preparing uniform qmesh to calculate gruneisen parameters" << apl::endl;
  xstructure xstr(xs);

  // Setup both lattices
  xstr.FixLattices();
  _rlattice = xstr.lattice;
  aurostd::xmatrix<double> _klattice = ReciprocalLattice(_rlattice);

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
  //_rlattice.clear();
  _klattice.clear();

  //populating atomic masses and species name
  ATOMIC_MASSES_AMU.clear();
  ATOMIC_SPECIES.clear();
  ATOMIC_MASSES_AMU = _pc.get_ATOMIC_MASSES_AMU();
  for (uint i = 0; i != xstr.atoms.size(); i++) {
    ATOMIC_SPECIES.push_back(xstr.atoms[i].name);
  }
  _meshsize[1] = na;
  _meshsize[2] = nb;
  _meshsize[3] = nc;
}
// ***************************************************************************************
bool QuasiHarmonicGruneisen::cal_gp_along_path(const std::vector<aurostd::xvector<double> > &kpoints) {
  _logger << "Calculating Gruneisen parameters along path" << apl::endl;
  _kpoints.clear();
  _kpoints = kpoints;
  if (_kpoints.size() == 0) {
    _logger << apl::error << "cal_gp_along_path() No qpoints to calculate" << apl::endl;
    return false;
  }
  return get_dynamicalmatrices_along_path();
}
// ***************************************************************************************
bool QuasiHarmonicGruneisen::cal_gp_in_mesh() {
  return get_dynamicalmatrices_in_mesh();
}
// ***************************************************************************************
bool QuasiHarmonicGruneisen::get_dynamicalmatrices_in_mesh() {
  _logger << "Calculating gruneisen parameters in mesh" << apl::endl;
  DMp.clear();
  DMm.clear();
  DM.clear();
  _eigenvectors.clear();
  xmatrix<xcomplex<double> > dd(nBranches, nBranches, 1, 1);
  DMp.resize(_kpoints.size(), dd);  //at +ve volume
  DMm.resize(_kpoints.size(), dd);  //at -ve volume
  DM.resize(_kpoints.size(), dd);
  _eigenvectors.resize(_kpoints.size(), dd);

  //dynamic matrix file names
  string DMmfile = string(dirfrefix) + string("/") + string("DM_") + string("mesh.") + string(gpdir[0]);
  string DMpfile = string(dirfrefix) + string("/") + string("DM_") + string("mesh.") + string(gpdir[1]);

  _logger << "Reading dynamical matrices from files" << apl::endl;
  //reading dynamical matrices
  if (!read_matrix(DMp, DMpfile)) return false;
  if (!read_matrix(DMm, DMmfile)) return false;

  _GP_mesh.clear();
  _freqs_mesh.clear();
  xvector<double> d(nBranches, 1);
  _GP_mesh.resize(_kpoints.size(), d);
  _freqs_mesh.resize(_kpoints.size(), d);

  bool gppass = false;
  gp_mesh_test.clear();
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  int qpointsPerCPU = _kpoints.size() / ncpus;
  gp_mesh_test.resize(ncpus, false);
  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating Gruneisen parameters in mesh");
  else
    _logger.initProgressBar("Calculating Gruneisen parameters in mesh (" + stringify(ncpus) + " threads)");

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread *> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _kpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _kpoints.size())))
      endIndex = _kpoints.size();
    threads.push_back(new std::thread(&QuasiHarmonicGruneisen::calculate_gp_in_mesh, this, startIndex, endIndex, icpu));
  }
  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }

  // Done
  _logger.finishProgressBar();

  uint bool_size = 0;
  for (uint id = 0; id != gp_mesh_test.size(); id++)
    bool_size += gp_mesh_test[id];
  if (bool_size == (uint)ncpus) gppass = true;
  gp_mesh_test.clear();
#else
  gp_mesh_test.resize(1, false);
  calculate_gp_in_mesh(0, (int)_kpoints.size(), 0);
  uint bool_size = 0;
  for (uint id = 0; id != gp_mesh_test.size(); id++)
    bool_size += gp_mesh_test[id];
  if (bool_size == 1) gppass = true;
#endif
  return gppass;
}
// ***************************************************************************************
//read dynamical matrices generated by other directories
bool QuasiHarmonicGruneisen::get_dynamicalmatrices_along_path() {
  DMp.clear();
  DMm.clear(), DM.clear();
  xmatrix<xcomplex<double> > dd(nBranches, nBranches, 1, 1);
  DMp.resize(_kpoints.size(), dd);  //at +ve volume
  DMm.resize(_kpoints.size(), dd);  //at -ve volume
  DM.resize(_kpoints.size(), dd);

  //dynamic matrix file names
  string DMmfile = string(dirfrefix) + string("/") + string("DM_") + string("path.") + string(gpdir[0]);
  string DMpfile = string(dirfrefix) + string("/") + string("DM_") + string("path.") + string(gpdir[1]);

  _logger << "Reading dynamical matrices for high-symmetry path" << apl::endl;
  if (!read_matrix(DMp, DMpfile)) return false;
  if (!read_matrix(DMm, DMmfile)) return false;

  GP_path.clear();
  _freqs_path.clear();
  xvector<double> d(nBranches, 1);
  GP_path.resize(_kpoints.size(), d);
  _freqs_path.resize(_kpoints.size(), d);

  bool gppass = false;
  gp_path_test.clear();

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
  int qpointsPerCPU = _kpoints.size() / ncpus;
  gp_path_test.resize(ncpus, false);
  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating Gruneisen parameters along path");
  else
    _logger.initProgressBar("Calculating Gruneisen parameters along path (" + stringify(ncpus) + " threads)");

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread *> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _kpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _kpoints.size())))
      endIndex = _kpoints.size();
    threads.push_back(new std::thread(&QuasiHarmonicGruneisen::calculate_gp_in_path, this, startIndex, endIndex, icpu));
  }
  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }

  // Done
  _logger.finishProgressBar();

  uint bool_size = 0;
  for (uint id = 0; id != gp_path_test.size(); id++)
    bool_size += gp_path_test[id];
  if (bool_size == (uint)ncpus) gppass = true;
  gp_path_test.clear();
#else
  gp_path_test.resize(1, false);
  calculate_gp_in_path(0, (int)_kpoints.size(), 0);
  uint bool_size = 0;
  for (uint id = 0; id != gp_path_test.size(); id++)
    bool_size += gp_path_test[id];
  if (bool_size == 1) gppass = true;
#endif

  return gppass;
}
// ***************************************************************************************
void QuasiHarmonicGruneisen::calculate_gp_in_mesh(int startIndex, int endIndex, int cpuid) {
  for (int In = startIndex; In < endIndex; In++) {
    xvector<double> qpoint;
    xmatrix<xcomplex<double> > DM0(nBranches, nBranches, 1, 1);  //at volume 0
    qpoint = _kpoints[In];
    DM0 = _pc.getDynamicMatrix(qpoint);
    DM[In] = DM0;

    xvector<double> eigenvalues(nBranches, 1);
    xmatrix<xcomplex<double> > eigenvectors(nBranches, nBranches, 1, 1);

    //calculate eigenvalue and eigenvectors of a dynamical matirx
    apl::aplEigensystems e;
    e.eigen_calculation(DM0, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);

    for (uint j = 1; j <= nBranches; j++) {
      if (_iszero(eigenvalues[j])) eigenvalues[j] = 0.0;
      if (eigenvalues[j] < 0) {
        if (_isnegative(eigenvalues[j])) {
          _logger << apl::warning << "gruneisen negative eigenvalue: = " << eigenvalues[j] << apl::endl;
          return;
        } else
          eigenvalues[j] = 0.00;
      }
      _freqs_mesh[In][j] = sqrt(eigenvalues[j]) * RAW2Hz;
    }  // nBranch loop end

    _eigenvectors[In] = eigenvectors;

    //calculate gruneisen
    _GP_mesh[In] = calculate_gruneisen(DM0, DMp[In], DMm[In], delta_V, V0);
  }

  gp_mesh_test[cpuid] = true;
}  //fn end
// ***************************************************************************************
void QuasiHarmonicGruneisen::calculate_gp_in_path(int startIndex, int endIndex, int cpuid) {
  for (int In = startIndex; In < endIndex; In++) {
    xvector<double> qpoint;
    xmatrix<xcomplex<double> > DM0(nBranches, nBranches, 1, 1);  //at volume 0
    qpoint = _kpoints[In];
    DM0 = _pc.getDynamicMatrix(qpoint);
    DM[In] = DM0;

    xvector<double> eigenvalues(nBranches, 1);
    xmatrix<xcomplex<double> > eigenvectors(nBranches, nBranches, 1, 1);

    for (uint j = 1; j <= nBranches; j++) {
      if (_iszero(eigenvalues[j])) eigenvalues[j] = 0.0;
      if (eigenvalues[j] < 0) {
        if (_isnegative(eigenvalues[j])) {
          _negative_freq_along_hsp = true;
          _logger << apl::warning << "gruneisen negative eigenvalue: = " << eigenvalues[j] << apl::endl;
          return;
        } else
          eigenvalues[j] = 0.00;
      }
      _freqs_path[In][j] = sqrt(eigenvalues[j]) * RAW2Hz;
    }  // nBranch loop end

    GP_path[In] = calculate_gruneisen(DM0, DMp[In], DMm[In], delta_V, V0);
  }

  gp_path_test[cpuid] = true;
}  //fn end
// ***************************************************************************************
void QuasiHarmonicGruneisen::write_GP_path(const vector<double> &path, const vector<int> &path_seg) {
  if (GP_path.size() == 0) {
    _logger << apl::error << "GP_path size is zero " << apl::endl;
    return;
  } else if (path.size() == 0) {
    _logger << apl::error << "path size is zero " << apl::endl;
    return;
  } else if (path_seg.size() == 0) {
    _logger << apl::error << "path_seg size is zero " << apl::endl;
    return;
  }

  if (_negative_freq_along_hsp) {
    _logger << apl::warning << "Negative frequencies and pass Gruneisen calculations " << apl::endl;
    return;
  }

  _logger << "Writing Gruneisen along path into file aflow.apl.gp.pdis.out ," << apl::endl;

  vector<string> hash_lines;
  if (!read_PDIS(hash_lines)) {
    _logger << apl::error << "PDIS file absent " << apl::endl;
    return;
  }

  //writing to a file
  stringstream os_gp;
  uint isize = GP_path.size();
  uint jsize = GP_path[0].rows;

  os_gp << "#Gruneisen Paramaters along high-symmetry q-points. File format is same as PDIS file"
        << "\n";
  for (uint i = 1; i != hash_lines.size(); i++)
    os_gp << hash_lines[i] << "\n";

  for (uint j = 1; j <= jsize; j++) {
    if (j == 1)
      os_gp << "#" << setw(10) << "        " << setw(10) << "q-path" << setw(10) << string("Br") + NumberToString<uint>(j);
    else
      os_gp << setw(15) << string("Br") + NumberToString<uint>(j);
  }
  os_gp << "\n";

  for (uint i = 0; i < isize; i++) {
    os_gp << setw(5) << path_seg[i];
    os_gp << setprecision(6) << std::fixed << std::showpoint
          << setw(15) << path[i];

    for (uint j = 1; j <= jsize; j++) {
      os_gp << setprecision(6) << std::fixed << std::showpoint
            << setw(15) << GP_path[i][j];
    }
    os_gp << "\n";
  }

  string outfile = "aflow.apl.gp.pdis.out";
  if (!aurostd::stringstream2file(os_gp, outfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.gp.pdis.out");
  }
  aurostd::StringstreamClean(os_gp);
}  //fn end
// ***************************************************************************************
void QuasiHarmonicGruneisen::write_GP_mesh() {
  if (_GP_mesh.size() == 0) {
    _logger << apl::error << "average GP can't be calculated since GP sizes are zero" << apl::endl;
    return;
  }
  _logger << "Writing Gruneisen parameters into file aflow.apl.gp.mesh.out," << apl::endl;

  stringstream os_gp;
  uint isize = _GP_mesh.size();
  uint jsize = _GP_mesh[0].rows;
  std::string STAR = std::string(10 * nBranches, '*');
  std::string STAR50 = std::string(50, '*');

  os_gp << "[AFLOW] " << STAR << "\n";
  os_gp << "[APL_GRUNEISEN_MESH]START"
        << "\n";

  for (uint j = 1; j <= jsize; j++) {
    if (j == 1)
      os_gp << "#" << setw(10) << string("Br") + NumberToString<uint>(j);
    else
      os_gp << setw(15) << string("Br") + NumberToString<uint>(j);
  }
  os_gp << "\n";

  for (uint i = 0; i < isize; i++) {
    for (uint j = 1; j <= jsize; j++) {
      os_gp << setprecision(6) << std::fixed << std::showpoint
            << setw(15) << _GP_mesh[i][j];
    }
    os_gp << "\n";
  }

  os_gp << "[APL_GRUNEISEN_MESH]STOP"
        << "\n";
  os_gp << "[AFLOW] " << STAR << "\n";

  string outfile = "aflow.apl.gp.mesh.out";
  if (!aurostd::stringstream2file(os_gp, outfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.gp.mesh.out");
  }
  aurostd::StringstreamClean(os_gp);

  stringstream os_qpoint;
  os_qpoint << "[AFLOW] " << STAR50 << "\n";
  os_qpoint << "[APL_QPOINT_MESH]START"
            << "\n";
  os_qpoint << "#" << aurostd::PaddedPRE("kx", 12, " ")
            << aurostd::PaddedPRE("ky", 12, " ")
            << aurostd::PaddedPRE("kz", 15, " ")
            << aurostd::PaddedPRE("weights", 20, " ") << "\n";

  for (uint i = 0; i != _kpoints.size(); i++) {
    for (int j = 1; j <= _kpoints[i].rows; j++) {
      os_qpoint << setprecision(6) << std::fixed << std::showpoint
                << setw(15) << _kpoints[i][j];
    }
    os_qpoint << setprecision(6) << std::fixed << std::showpoint
              << setw(15) << _weights[i] << "\n";
  }
  os_qpoint << "[APL_QPOINT_MESH]END"
            << "\n";
  os_qpoint << "[AFLOW] " << STAR50 << "\n";

  outfile = "aflow.apl.qpoints.mesh.out";
  if (!aurostd::stringstream2file(os_qpoint, outfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.qpoints.mesh.out");
  }
  aurostd::StringstreamClean(os_qpoint);

}  //fn end
// ***************************************************************************************
void QuasiHarmonicGruneisen::write_GP_FREQ_mesh() {
  if (_GP_mesh.size() == 0) {
    _logger << apl::error << "aflow.apl.gp.freq.mesh.out can't be calculated since GP sizes are zero" << apl::endl;
    return;
  }

  _logger << "Writing Gruneisen parameters and frequency into file aflow.apl.gp.freq.mesh.out," << apl::endl;

  stringstream os_gp;
  uint isize = _GP_mesh.size();
  uint jsize = _GP_mesh[0].rows;
  std::string STAR = std::string(5 * nBranches, '*');

  os_gp << "[AFLOW] " << STAR << "\n";
  os_gp << "[APL_GRUNEISEN_MESH]START"
        << "\n";

  os_gp << "#" << setw(15) << "Gruneisen" << setw(15) << "Freq (THz) "
        << "\n";

  for (uint i = 0; i < isize; i++) {
    for (uint j = 1; j <= jsize; j++) {
      os_gp << setprecision(6) << std::fixed << std::showpoint
            << setw(15) << _GP_mesh[i][j] << "      " << _freqs_mesh[i][j] << "\n";
    }
    os_gp << "\n";
  }

  os_gp << "[APL_GRUNEISEN_MESH]STOP"
        << "\n";
  os_gp << "[AFLOW] " << STAR << "\n";

  string outfile = "aflow.apl.gp.freq.mesh.out";
  if (!aurostd::stringstream2file(os_gp, outfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.gp.freq.mesh.out");
  }
  aurostd::StringstreamClean(os_gp);
}  //fn end
// ***************************************************************************************
//Write average Gruneisen parameter to a file
void QuasiHarmonicGruneisen::WriteaverageGP(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP) {
  if (_GP_mesh.size() == 0) {
    _logger << apl::error << "_GP_mesh size zero" << apl::endl;
    return;
  }
  if (_acoustic_GP_mesh.size() == 0) {
    _logger << apl::error << "_acoustic_GP_mesh size zero" << apl::endl;
    return;
  }

  _logger << "Writing Gruneisen into file aflow.apl.qh.out ," << apl::endl;

  std::string STAR40 = std::string(40, '*');

  stringstream os_avg;
  os_avg << "[AFLOW] " << STAR40 << "\n";
  os_avg << "[APL_GRUNEISEN_300K]START"
         << "\n";
  os_avg << "apl_average_gruneisen_300K = " << averageGP(300.0) << "\n";
  os_avg << "apl_average_acoustic_gruneisen_300K = " << AcousticAverageGP(300.0) << "\n";
  os_avg << "[APL_GRUNEISEN_300K]STOP"
         << "\n";
  os_avg << "[AFLOW] " << STAR40 << "\n";
  os_avg << "# gamma    => Gruneisen parameter"
         << "\n";
  os_avg << "# ac_gamma => Acoustic Gruneisen parameter"
         << "\n";
  os_avg << "[AFLOW] " << STAR40 << "\n";
  os_avg << "[APL_GRUNEISEN]START"
         << "\n";
  os_avg << "#" << aurostd::PaddedPRE("T(K)", 14, " ")
         << aurostd::PaddedPRE("gamma", 10, " ")
         << aurostd::PaddedPRE("ac_gamma", 18, " ") << "\n";

  int n = (int)((USER_TP_TEND - USER_TP_TSTART) / USER_TP_TSTEP);
  for (int i = 0; i <= n; i++) {
    double TEMPERATURE_IN_K = USER_TP_TSTART + i * USER_TP_TSTEP;
    double avgGP = averageGP(TEMPERATURE_IN_K);
    double ACavgGP = AcousticAverageGP(TEMPERATURE_IN_K);
    os_avg << setw(15) << setprecision(2) << std::fixed << std::showpoint << TEMPERATURE_IN_K
           << setw(15) << setprecision(8) << std::fixed << std::showpoint
           << avgGP << setw(15) << ACavgGP << "\n";
  }

  os_avg << "[APL_GRUNEISEN]STOP"
         << "\n";
  os_avg << "[AFLOW] " << STAR40 << "\n";
  string avg_out = "aflow.apl.qh.out";
  if (!aurostd::stringstream2file(os_avg, avg_out, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.qh.out");
  }
  aurostd::StringstreamClean(os_avg);
}
// ***************************************************************************************
double QuasiHarmonicGruneisen::averageGP(double temperature_in_kelvins) {
  if (_GP_mesh.size() == 0) {
    _logger << apl::error << "_GP_mesh.size()==0" << apl::endl;
    exit(0);
  }
  if (_freqs_mesh.size() == 0) { _logger << apl::error << "_freqs_mesh.size()==0" << apl::endl; }

  double hnu = 4.135667516E2;  // h*10^(12)*10^(5) = h*10^(2)
  double kB = 8.6173324;
  double beta = 1. / (kB * temperature_in_kelvins);

  double cv = 0.00;
  double gpcv = 0.00;

  for (uint i = 0; i < _freqs_mesh.size(); i++) {
    for (int j = 1; j <= _freqs_mesh[i].rows; j++) {
      //double f= _freqs_mesh[i][j]* 15.6333046177; // 1/(2*pi) * sqrt(eV/A^2* Mass) ---> 15.6333046177 sec^(-1)
      double f = _freqs_mesh[i][j];
      if (f < CUTOFF_FREQ) continue;
      double x = hnu * f * beta;
      double ex = exp(x);
      double cvij = x * x * ex * _weights[i] / ((1. - ex) * (1. - ex));
      cv = cv + cvij;
      gpcv = gpcv + _GP_mesh[i][j] * cvij;
    }
  }
  if (std::isnan(gpcv / cv)) return 0.00;
  return gpcv / cv;
}
// ***************************************************************************************
double QuasiHarmonicGruneisen::AcousticAverageGP(double temperature_in_kelvins) {
  if (_acoustic_GP_mesh.size() == 0) {
    _logger << apl::error << "_acoustic_GP_mesh.size()==0" << apl::endl;
    exit(0);
  }
  if (_acoustic_freqs_mesh.size() == 0) {
    _logger << apl::error << "_acoustic_freqs_mesh.size()==0" << apl::endl;
    exit(0);
  }

  double hnu = 4.135667516E2;  // h*10^(12)*10^(5) = h*10^(2)
  double kB = 8.6173324;
  double beta = 1. / (kB * temperature_in_kelvins);

  double cv = 0.00;
  double gpcv = 0.00;

  for (uint i = 0; i < _acoustic_freqs_mesh.size(); i++) {
    for (int j = 1; j <= _acoustic_freqs_mesh[i].rows; j++) {
      //double f= _freqs_mesh[i][j]* 15.6333046177; // 1/(2*pi) * sqrt(eV/A^2* Mass) ---> 15.6333046177 sec^(-1)
      double f = _acoustic_freqs_mesh[i][j];
      if (f < CUTOFF_FREQ) continue;
      double x = hnu * f * beta;
      double ex = exp(x);
      double cvij = x * x * ex / ((1. - ex) * (1. - ex));
      cv = cv + cvij;
      gpcv = gpcv + _acoustic_GP_mesh[i][j] * cvij;
    }
  }
  if (std::isnan(gpcv / cv)) return 0.00;
  return gpcv / cv;
}
// ***************************************************************************************
template <typename T>
std::vector<T> QuasiHarmonicGruneisen::split(const std::string &line) {
  std::istringstream is(line);
  return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
}
// ***************************************************************************************
bool QuasiHarmonicGruneisen::read_matrix(vector<xmatrix<xcomplex<double> > > &A, string file) {
  //CO - START
  //string command="";
  //string EXT=string(file)+string(".EXT");
  //bool t1=aurostd::FileExist(file);
  //bool t2=aurostd::FileExist(EXT);
  //if(!t1)
  //  {
  //if(t2){command=string("COMPRESS -d ")+string(EXT);aurostd::execute2string(command);}
  //else{_logger<<apl::error <<file<<" not found "<<apl::endl;return false;}
  //    }
  string compressed_file = "";
  bool is_compressed = aurostd::IsCompressed(file, compressed_file);
  //unavoidable here, we need to uncompress then recompress later
  if (is_compressed) {
    aurostd::UncompressFile(compressed_file);
  }
  //CO - END

  ifstream in;
  in.open(file.c_str(), ios_base::binary);
  if (!in.is_open()) {
    _logger << apl::error << file << " not able to open " << apl::endl;
    return false;
  }
  for (uint I = 0; I < A.size(); I++) {
    for (uint J = 1; J <= nBranches; J++) {
      in.read((char *)(&A[I][J][1]), nBranches * sizeof(xcomplex<double>));
    }
  }
  in.clear();
  in.close();
  //if((t1+t2)!=2){command=string("COMPRESS ")+string(file);aurostd::execute2string(command);}
  if (is_compressed) {                                //CO
    aurostd::MatchCompressed(compressed_file, file);  //CO
  }                                                   //CO
  return true;
}
// ***************************************************************************************
template <typename T>
string QuasiHarmonicGruneisen::NumberToString(T Number) {
  ostringstream ss;
  ss << Number;
  return ss.str();
}
// ***************************************************************************************
bool QuasiHarmonicGruneisen::read_PDIS(vector<string> &hash_lines) {
  string file = "PDIS";
  //CO - START
  //string command="";
  //string EXT=string(file)+string(".EXT");
  //bool t1=aurostd::FileExist(file);
  //bool t2=aurostd::FileExist(EXT);
  //if(!t1)
  //  {
  //if(t2){command=string("COMPRESS -d ")+string(EXT);aurostd::execute2string(command);}
  //else{_logger<<apl::error <<file<<" absent"<<apl::endl; return false; }
  //    }
  if (!aurostd::EFileExist(file)) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }
  vector<string> vlines;
  aurostd::efile2vectorstring(file, vlines);
  if (!vlines.size()) {
    _logger << apl::error << file << " does not exist" << apl::endl;
    exit(0);
  }

  hash_lines.clear();
  string line;
  uint line_count = 0;
  //ifstream in (file.c_str());
  //if (!in.is_open()){_logger<<apl::error <<file<<" not able to open "<<apl::endl; return false;}
  if (!vlines.size()) {
    _logger << apl::error << file << " not able to open " << apl::endl;
    return false;
  }
  //while ( getline (in,line) ){
  while (line_count < vlines.size()) {
    line = vlines[line_count++];
    if (line == "") continue;
    if (line[0] == '#')
      hash_lines.push_back(line);
    else
      break;
  }
  //in.clear();
  //in.close();
  //CO - END
  //if((t1+t2)!=2){command=string("COMPRESS ")+string(file);aurostd::execute2string(command);}
  return true;
}
// ***************************************************************************************
_CMAT_ QuasiHarmonicGruneisen::xmat2mat(const xmatrix<xcomplex<double> > &M) {
  _CMAT_ m(M.cols, _CVEC_(M.rows, _CD_(0.0, 0.0)));
  for (int i = 0; i < M.rows; i++)
    for (int j = 0; j < M.cols; j++)
      m[i][j] = _CD_(M[i + 1][j + 1].re, M[i + 1][j + 1].im);

  return m;
}
// ***************************************************************************************
_VEC_ QuasiHarmonicGruneisen::xvec2vec(const xvector<double> &V) {
  _VEC_ v(V.rows, 0.0);
  for (int i = 0; i < V.rows; i++)
    v[i] = V[i + 1];

  return v;
}
// ***************************************************************************************
_MAT_ QuasiHarmonicGruneisen::xmatd2matd(const xmatrix<double> &M) {
  _MAT_ m(M.cols, _VEC_(M.rows, 0.0));
  for (int i = 0; i < M.rows; i++)
    for (int j = 0; j < M.cols; j++)
      m[i][j] = M[i + 1][j + 1];

  return m;
}
// ***************************************************************************************
//thermal displacements calculation for a given range of temperature
//formulas can be found in the link http://atztogo.github.io/phonopy/thermal-displacement.html
void QuasiHarmonicGruneisen::thermal_displacements(double Ts, double Te, int Tinc) {
  _logger << "calculating mean square displacements " << apl::endl;
  if (ATOMIC_SPECIES.size() != ATOMIC_MASSES_AMU.size()) {
    _logger << apl::error << "error in atomic masses [remove apl.xml LOCK and run again] " << apl::endl;
    return;
  }

  uint qmesh_size = _eigenvectors.size();

  stringstream os_disp;
  os_disp << "# Cartesian atomic meansquare displacements in Angstrom uint\n";

  std::string STAR = std::string(5 * nBranches, '*');

  os_disp << "[AFLOW] " << STAR << "\n";
  os_disp << "[DISPLACEMENTS]START"
          << "\n";
  os_disp << "#" << setw(9) << "T(K)" << setw(15) << "AtomType" << setw(15) << "xdir" << setw(15) << "ydir" << setw(15) << "zdir"
          << "\n";
  os_disp << std::fixed << std::showpoint;

  for (double t = Ts; t <= Te; t = t + Tinc) {
    os_disp << setw(10) << setprecision(2) << t << "\n";
    _VEC_ disps(nBranches, 0.0);

    for (uint kvec = 0; kvec != qmesh_size; kvec++) {
      _CMAT_ eigenvectors = xmat2mat(_eigenvectors[kvec]);
      _CMAT_ vecs2(nBranches,
                   _CVEC_(nBranches, _CD_(0., 0.)));

      for (uint i = 0; i != nBranches; i++) {
        for (uint j = 0; j != nBranches; j++) {
          vecs2[i][j] = std::abs(eigenvectors[i][j]) *
                        std::abs(eigenvectors[i][j]);
        }
      }

      _MAT_ vecs2_T(nBranches,
                    _VEC_(nBranches, 0.0));

      for (uint i = 0; i != nBranches; i++) {
        for (uint j = 0; j != nBranches; j++) {
          vecs2_T[i][j] = vecs2[j][i].real();
        }
      }

      _VEC_ f = xvec2vec(_freqs_mesh[kvec]);

      for (uint i = 0; i != nBranches; i++) {
        uint cnt = 0;
        //run through all atoms and their directions
        for (uint atom_alpha = 0; atom_alpha != nBranches; atom_alpha++) {
          if (f[i] < CUTOFF_FREQ) continue;
          if ((atom_alpha % 3 == 0) && (atom_alpha != 0)) cnt++;
          double c = vecs2_T[i][atom_alpha] / (ATOMIC_MASSES_AMU[cnt] * AMU2Kg);
          disps[atom_alpha] += get_Q2(f[i], t) * c;
        }
      }
    }  //kvec loop

    os_disp << setprecision(8);
    size_t atom_cnt = 0;
    for (uint atom_alpha = 0; atom_alpha != nBranches; atom_alpha++) {
      if ((atom_alpha % 3 == 0) && atom_alpha != 0) {
        os_disp << "\n";
      }
      if (atom_alpha % 3 == 0) {
        os_disp << setw(25) << ATOMIC_SPECIES[atom_cnt];
        atom_cnt++;
      }
      os_disp << setw(15) << disps[atom_alpha] / (double)qmesh_size;
    }
    os_disp << "\n";
  }  //t loop

  os_disp << "[AFLOW] " << STAR << "\n";
  os_disp << "[DISPLACEMENTS]STOP"
          << "\n";

  string outfile = "aflow.apl.displacements.out";
  if (!aurostd::stringstream2file(os_disp, outfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.displacements.out");
  }
  aurostd::StringstreamClean(os_disp);
}
// ***************************************************************************************
//projected displacements along [hkl] direction for a given range of temperature
//formulas can be found in the link http://atztogo.github.io/phonopy/thermal-displacement.html
void QuasiHarmonicGruneisen::projected_displacement(const _VEC_ &direction, double Ts, double Te, int Tinc) {
  _logger << "calculating projected mean square displacements along [" << NumberToString<int>((int)direction[0])
          << " " << NumberToString<int>((int)direction[1]) << " " << NumberToString<int>((int)direction[2]) << "]" << apl::endl;

  uint qmesh_size = _eigenvectors.size();

  if (qmesh_size == 0) {
    _logger << apl::error << " eigenvectors are zero " << apl::endl;
    return;
  }

  if (ATOMIC_SPECIES.size() != ATOMIC_MASSES_AMU.size()) {
    _logger << apl::error << "error in atomic masses [remove apl.xml LOCK and run again]" << apl::endl;
    return;
  }

  vector<double> projector(3, 0.);
  _MAT_ lattice = xmatd2matd(_rlattice);

  //get projected vector
  for (uint i = 0; i != lattice.size(); i++) {
    vector<double> a = lattice[i];
    projector[i] = apl_inner_product(a, direction);
  }

  double norm = 0.0;
  for (uint i = 0; i != projector.size(); i++)
    norm += projector[i] * projector[i];
  norm = sqrt(norm);
  for (uint i = 0; i != projector.size(); i++)
    projector[i] /= norm;

  //make transpose of eigenvector
  vector<_CMAT_> p_eigenvectors;
  for (uint i = 0; i != qmesh_size; i++)  //sum over qpoints
  {
    //transpose of eigenvector
    _CMAT_ EVE = xmat2mat(_eigenvectors[i]);
    _CMAT_ EVE_T(nBranches, _CVEC_(nBranches, _CD_(0.0, 0.0)));

    for (uint j = 0; j != nBranches; j++)
      for (uint k = 0; k != nBranches; k++)
        EVE_T[j][k] = EVE[k][j];

    _CMAT_ tmp2d;
    for (uint j = 0; j != nBranches; j++) {
      //seperate real and imeginary parts
      _VEC_ vecsR(nBranches, 0.00);
      _VEC_ vecsI(nBranches, 0.00);
      for (uint k = 0; k < nBranches; k++) {
        vecsR[k] = EVE_T[j][k].real();
        vecsI[k] = EVE_T[j][k].imag();
      }
      _CVEC_ tmp1d;
      for (uint k = 0; k < nBranches; k++) {
        if (k % 3 == 0) {
              vector<double> tmp1(vecsR.begin()+k, vecsR.begin()+k+3);
              vector<double> tmp2(vecsI.begin()+k, vecsI.begin()+k+3);
              double real=apl_inner_product(tmp1, projector);
              double imag=apl_inner_product(tmp2, projector);
          tmp1d.push_back(_CD_(real, imag));
        }
      }
      tmp2d.push_back(tmp1d);
    }
    p_eigenvectors.push_back(tmp2d);
  }

  uint nATOMs = ATOMIC_SPECIES.size();
  _CVEC_ tmp1d(nBranches, _CD_(0., 0.));
  _CMAT_ tmp2d(nATOMs, tmp1d);
  vector<_CMAT_> p_eigenvectorsT(qmesh_size, tmp2d);

  for (uint i = 0; i != qmesh_size; i++) {
    for (uint j = 0; j != nBranches; j++) {
      for (uint k = 0; k != nATOMs; k++) {
        p_eigenvectorsT[i][k][j] = p_eigenvectors[i][j][k];
      }
    }
  }

  stringstream os_disp;
  os_disp << "# Cartesian projected atomic meansquare displacements in Angstrom uint \n";
  std::string STAR = std::string(10 * nBranches, '*');
  os_disp << "[AFLOW] " << STAR << "\n";
  os_disp << "[PROJECTED_DISPLACEMENTS]START"
          << "\n";
  os_disp << setw(10) << "# Projection direction set to [" << (int)direction[0] << " " << (int)direction[1] << " " << (int)direction[2] << "]\n";

  os_disp << "#" << setw(14) << "T(K)";
  for (uint i = 0; i != nATOMs; i++)
    os_disp << setw(15) << ATOMIC_SPECIES[i];
  os_disp << "\n";

  //temperature dependent displacements
  for (double t = Ts; t <= Te; t = t + Tinc) {
    _VEC_ disps(nATOMs, 0.00);
    for (uint kvec = 0; kvec != qmesh_size; kvec++)  //qpoint sum
    {
      _CMAT_ eigenvectors = p_eigenvectorsT[kvec];
      _CMAT_ vecs2(nATOMs, _CVEC_(nBranches, _CD_(0., 0.)));

      for (uint i = 0; i != nATOMs; i++) {
        for (uint j = 0; j != nBranches; j++) {
          vecs2[i][j] = std::abs(eigenvectors[i][j]) *
                        std::abs(eigenvectors[i][j]);
        }
      }

      //transpose of vecs2 and its complex parts are zero
      _MAT_ vecs2_T(nBranches, _VEC_(nATOMs, 0.0));

      for (uint i = 0; i != nATOMs; i++) {
        for (uint j = 0; j != nBranches; j++) {
          vecs2_T[j][i] = vecs2[i][j].real();
        }
      }

      _VEC_ f = xvec2vec(_freqs_mesh[kvec]);

      for (uint i = 0; i != nBranches; i++) {
        for (uint j = 0; j < nATOMs; j++) {
          if (f[i] < CUTOFF_FREQ) continue;
          double c = vecs2_T[i][j] / (ATOMIC_MASSES_AMU[j] * AMU2Kg);
          disps[j] += get_Q2(f[i], t) * c;
        }
      }
    }  //qpoint loop
    //PRINT
    for (uint i = 0; i != nATOMs; i++)
      disps[i] /= (double)qmesh_size;

    os_disp << setw(15) << std::setprecision(2) << t << setw(15);

    for (uint i = 0; i != nATOMs; i++) {
      os_disp << std::fixed << std::showpoint;
      os_disp << std::setprecision(8) << disps[i] << setw(15);
    }
    os_disp << "\n";
  }
  os_disp << "[AFLOW] " << STAR << "\n";
  os_disp << "[PROJECTED_DISPLACEMENTS]STOP"
          << "\n";

  string outfile = "aflow.apl.projected_displacements.out";
  if (!aurostd::stringstream2file(os_disp, outfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.projected_displacements.out");
  }
  aurostd::StringstreamClean(os_disp);
}
// ***************************************************************************************
double QuasiHarmonicGruneisen::get_population(double freq, double t) {
  if (t < 1.)
    return 0.0;
  else {
    return 1.0 / (exp((freq * 0.00413566733) / (8.61733825681e-05 * t)) - 1.);
  }
}
//   *******************************************************************************************  //
double QuasiHarmonicGruneisen::get_Q2(double freq, double t) {
  return (10.545721821978764) * (  //hbar*ev/A^2//
                                    (get_population(freq, t) + 0.5) / (freq * 2. * M_PI));
}
// ***************************************************************************************
void QuasiHarmonicGruneisen::writeDM() {
  //CO - START
  //ofstream myfile;
  stringstream myfile;
  myfile << std::fixed << std::showpoint;
  myfile << std::setprecision(8);
  //myfile.open ("mesh_dm");
  //CO - END
  myfile << "#nBranches: " << nBranches << "\n";
  myfile << "#qmesh.size: " << _kpoints.size() << "\n";
  myfile << std::endl;

  for (uint In = 0; In < _kpoints.size(); In++) {
    for (uint J = 1; J <= nBranches; J++) {
      for (uint K = 1; K <= nBranches; K++) {
        myfile << DM[In][J][K] << " ";
      }
      myfile << "\n";
    }
  }
  myfile << "\n";
  //CO - START
  string filename = "mesh_dm";
  aurostd::stringstream2file(myfile, filename);
  //myfile.clear();
  //myfile.close();
  //CO - END
}
// ***************************************************************************************
void QuasiHarmonicGruneisen::calculate_acoustic_freqNgp() {
  _logger << "Calculation acoustic frequencies and Gruneisen" << apl::endl;
  identify_acoustic_modes();
}
// ***************************************************************************************
void QuasiHarmonicGruneisen::identify_acoustic_modes() {
  //clear variables before use
  _BrINDEXs.clear();
  _acoustic_freqs_mesh.clear();
  _acoustic_GP_mesh.clear();

  std::vector<aurostd::xvector<double> > freqs_test = _freqs_mesh;
  std::vector<aurostd::xvector<double> > gp_test = _GP_mesh;

  xvector<int> d(nBranches, 1);
  _BrINDEXs.resize(_kpoints.size(), d);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  int qpointsPerCPU = _kpoints.size() / ncpus;
  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Indentifying acoustic phonon mode in q-mesh");
  else
    _logger.initProgressBar("Indentifying acoustic phonon mode in q-mesh  (" + stringify(ncpus) + " threads)");

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread *> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _kpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _kpoints.size())))
      endIndex = _kpoints.size();
    threads.push_back(new std::thread(&QuasiHarmonicGruneisen::trace_acoustic_modes_threads, this, startIndex, endIndex));
  }
  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }
  // Done
  _logger.finishProgressBar();
#else
  trace_acoustic_modes_threads(0, _kpoints.size());
#endif

  _acoustic_freqs_mesh.resize(_kpoints.size(), xvector<double>(3, 1));
  _acoustic_GP_mesh.resize(_kpoints.size(), xvector<double>(3, 1));

  for (uint i = 0; i != _kpoints.size(); i++) {
    for (int j = 1; j <= _BrINDEXs[i].rows; j++) {
      if (_BrINDEXs[i][j] == 0) {
        _acoustic_freqs_mesh[i][1] = freqs_test[i][_BrINDEXs[i][j] + 1];
        _acoustic_GP_mesh[i][1] = gp_test[i][_BrINDEXs[i][j] + 1];
      } else if (_BrINDEXs[i][j] == 1) {
        _acoustic_freqs_mesh[i][2] = freqs_test[i][_BrINDEXs[i][j] + 1];
        _acoustic_GP_mesh[i][2] = gp_test[i][_BrINDEXs[i][j] + 1];
      } else if (_BrINDEXs[i][j] == 2) {
        _acoustic_freqs_mesh[i][3] = freqs_test[i][_BrINDEXs[i][j] + 1];
        _acoustic_GP_mesh[i][3] = gp_test[i][_BrINDEXs[i][j] + 1];
      }
    }
  }
  freqs_test.clear();
  gp_test.clear();
}
// ***************************************************************************************
void QuasiHarmonicGruneisen::trace_acoustic_modes_threads(int startIndex, int endIndex) {
  int qxqy = getN(1) * getN(2);
  double dkx = abs(_kpoints[2][1] - _kpoints[1][1]);
  double dky = abs(_kpoints[2][2] - _kpoints[1][2]);
  double dkz = abs(_kpoints[2][3] - _kpoints[1][3]);
  for (int iqp = startIndex; iqp < endIndex; iqp++) {
    if (iqp < (int)(_kpoints.size() - 1)) {
      double dkx1 = abs(_kpoints[iqp + 1][1] - _kpoints[iqp][1]);
      double dky1 = abs(_kpoints[iqp + 1][2] - _kpoints[iqp][2]);
      double dkz1 = abs(_kpoints[iqp + 1][3] - _kpoints[iqp][3]);
      xvector<int> index(nBranches, 1);
      if ((iqp + 1) % qxqy == 0)  //avoid reciprocal lattice boundaries
        index = trace_acoustic_mode(_eigenvectors[iqp], _eigenvectors[iqp]);
      else if ((_isequal(dkx, dkx1) || _isequal(dky, dky1) || _isequal(dkz, dkz1)))
        index = trace_acoustic_mode(_eigenvectors[iqp], _eigenvectors[iqp + 1]);
      else
        index = trace_acoustic_mode(_eigenvectors[iqp], _eigenvectors[iqp]);
      _BrINDEXs[iqp] = index;
    } else
      for (int j = 1; j <= _BrINDEXs[iqp].rows; j++) _BrINDEXs[iqp][j] = j - 1;
  }
}
// ***************************************************************************************
void QuasiHarmonicGruneisen::Write_BrINDEXs() {
  if (_BrINDEXs.size() == 0) {
    _logger << apl::error << "_BrINDEXs size are zero" << apl::endl;
    return;
  }

  _logger << "Writing phonon branch indices  aflow.apl.branch.index.out ," << apl::endl;

  stringstream os_gp;
  uint isize = (uint)_BrINDEXs.size();
  uint jsize = (uint)_BrINDEXs[0].rows;
  std::string STAR = std::string(5 * nBranches, '*');

  os_gp << "[AFLOW] " << STAR << "\n";
  os_gp << "[APL_PHONON_BRANCH_INDEX]START"
        << "\n";

  for (uint j = 1; j <= jsize; j++) {
    if (j == 1)
      os_gp << "#" << setw(5) << string("Br") + NumberToString<uint>(j);
    else
      os_gp << setw(5) << string("Br") + NumberToString<uint>(j);
  }
  os_gp << "\n";

  for (uint i = 0; i != isize; i++) {
    for (uint j = 1; j <= jsize; j++) {
      os_gp << setw(5) << _BrINDEXs[i][j];
    }
    os_gp << "\n";
  }
  os_gp << "[APL_PHONON_BRANCH_INDEX]END"
        << "\n";
  os_gp << "[AFLOW] " << STAR << "\n";

  string outfile = "aflow.apl.branch.index.out";
  if (!aurostd::stringstream2file(os_gp, outfile, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.branch.index.out");
  }
  aurostd::StringstreamClean(os_gp);
}
// ***************************************************************************************
  double QuasiHarmonicGruneisen::apl_inner_product(const vector<double> &a, const vector<double> &b)
  {
    if(a.size()!=b.size())
        _logger<<apl::error<<"apl_inner_product() vector size must be equal" << apl::endl;

    double sum=0.00;
    for(uint i=0; i!=a.size(); i++)sum+=a[i]*b[i];
   return sum;
}
// ***************************************************************************************
}  //apl end
