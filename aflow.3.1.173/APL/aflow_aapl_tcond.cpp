// ***************************************************************************
// *                                                                         *
// *              AFlow JOSE J. PLATA   - Duke University 2017               *
// *                                                                         *
// ***************************************************************************
// tcond_contrib_jose.cpp
// written by JJPR 2013
// 2013: jose.plata@duke.edu

// These class computes thermal conductivity using anharmonic IFCs and  phonon dispersion curves.
// Usefull References:
// Appl. Phys Lett. 91, 231922 (2007)
// Phys Rev. B 80, 125203 (2009)
// Phys Rev. Lett. 109, 095901 (2012)
// ***************************************************************************

#include <map>
#include "aflow_apl.h"
// Some parts are written within the C++0x support in GCC, expecially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <sstream>
#include <iterator>

using namespace std;

namespace apl {

//Constructor////QPs////////////////////////////////////////////////////
_QPs::_QPs() {
}
//Constructor////IQPs////////////////////////////////////////////////////
_iQPs::_iQPs() {
}
//Destructor/////QPs///////////////////////////////////////////////////
_QPs::~_QPs() {
  free();
}
void _QPs::free() {
}
//Destructor/////IQPs///////////////////////////////////////////////////
_iQPs::~_iQPs() {
  free();
}
void _iQPs::free() {
}
// Copy constrcutor////QPs/////////////////////////////////////////////
const _QPs& _QPs::operator=(const _QPs& b) {  // operator=
  if (this != &b) {
    free();
  }
  return *this;
}
// Copy constrcutor////IQPs/////////////////////////////////////////////
const _iQPs& _iQPs::operator=(const _iQPs& b) {  // operator=
  if (this != &b) {
    free();
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

ThermalConductivityCalculator::ThermalConductivityCalculator(IPhononCalculator& pc, Supercell& sc, StrPairs& pa, Logger& l) : _pc(pc), _sc(sc), _pa(pa), _logger(l) {
  nBranches = _pc.getNumberOfBranches();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

ThermalConductivityCalculator::~ThermalConductivityCalculator() {
  clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::clear() {
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void ThermalConductivityCalculator::qgrid(int nx, int ny, int nz) {
  xstr = _sc.getPrimitiveStructure();

  // Outputs variables
  ostringstream aus;

  // Setup both lattices
  xstr.FixLattices();
  _rlattice = xstr.lattice;
  _klattice = ReciprocalLattice(_rlattice);

  // Generatin grid
  xvector<double> kpoint(3);
  _QPs kk;
  int counter = 0.0;
  for (int s = 1; s <= nz; s++) {
    vector<vector<int>> vec2;
    for (int r = 1; r <= ny; r++) {
      vector<int> vec1;
      for (int p = 1; p <= nx; p++) {
        kk.mapq(1) = ((s - 1.0));
        kk.mapq(2) = ((r - 1.0));
        kk.mapq(3) = ((p - 1.0));
        kpoint(1) = ((s - 1.0) / (nx));
        kpoint(2) = ((r - 1.0) / (ny));
        kpoint(3) = ((p - 1.0) / (nz));
        // Transform to cartesian coordinate
        vec1.push_back(counter);
        kk.q_dir = kpoint;
        kpoint = trasp(_klattice) * kpoint;
        kk.q = kpoint;
        kk.q_is_equivalent = FALSE;
        kk.index_ineq = counter;
        kk.index_fg = 0;
        QPs.push_back(kk);
        counter = counter + 1;
      }
      vec2.push_back(vec1);
    }
    mapQ.push_back(vec2);
  }

  // ****************** BEGIN Printing block *************************
  stringstream oss;
  oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
  oss << "[APL_THERMAL_QPOINTS]START" << std::endl;
  // ****************** END Printing block *************************

  if (!xstr.pgroupk_xtal_calculated) {
    ofstream fileDevNull("/dev/null");
    if (!fileDevNull.is_open()) {
      throw APLRuntimeError("MeshTC::MeshTC(): Cannot open output stream /dev/null.");
    }

    _aflags af;
    af.Directory = "./";
    af.QUIET = TRUE;
    xstr.LatticeReduction_avoid = TRUE;

    //SYM::CalculatePointGroupKlattice(fileDevNull, xstr, af, FALSE, FALSE, _logger.getOutputStream()); //CO 171215 pgroupk_xtal NOT pgroupk
    xstr.pgroupk_xtal_calculated = TRUE;

    fileDevNull.clear();
    fileDevNull.close();
  }

  // Checking that Uf and Uc are consistent
  for (uint i = 0; i < xstr.pgroupk_xtal.size(); i++) {
    xstr.pgroupk_xtal[i].Uf = inverse(trasp(_klattice)) * xstr.pgroupk_xtal[i].Uc * trasp(_klattice);
  }

  // Look for inequivalence
  xvector<double> newkpoint(3);
  for (uint i = 0; i < QPs.size(); i++) {
    for (uint j = i + 1; j < QPs.size(); j++) {
      for (uint symOpID = 0; symOpID < xstr.pgroupk_xtal.size(); symOpID++) {
        newkpoint = xstr.pgroupk_xtal[symOpID].Uf * QPs[i].q_dir;
        xvector<double> fdiff = newkpoint - QPs[j].q_dir;
        SYM::PBC(fdiff);
        if (aurostd::modulus(fdiff) < _AFLOW_APL_EPS_) {
          QPs[j].q_is_equivalent = TRUE;
          QPs[j].index_ineq = QPs[i].index_ineq;
          symOpID = xstr.pgroupk_xtal.size();
          j = QPs.size();
        }
      }
    }
  }

  for (uint i = 0; i < QPs.size(); i++) {
    for (uint symOpID = 0; symOpID < xstr.pgroupk_xtal.size(); symOpID++) {
      newkpoint = xstr.pgroupk_xtal[symOpID].Uf * QPs[QPs[i].index_ineq].q_dir;
      xvector<double> fdiff = newkpoint - QPs[i].q_dir;
      SYM::PBC(fdiff);
      if (aurostd::modulus(fdiff) < _AFLOW_APL_EPS_) {
        QPs[i].index_fg = symOpID;
        symOpID = xstr.pgroupk_xtal.size();
      }
    }
  }

  for (uint i = 0; i < QPs.size(); i++) {
    for (uint symOpID = 0; symOpID < xstr.pgroupk_xtal.size(); symOpID++) {
      newkpoint = xstr.pgroupk_xtal[symOpID].Uc * QPs[QPs[i].index_ineq].q;
      newkpoint = BringInCell(C2F(_klattice, newkpoint));
      xvector<double> fdiff = newkpoint - QPs[i].q_dir;
      SYM::PBC(fdiff);
      if (aurostd::modulus(fdiff) < _AFLOW_APL_EPS_) {
        QPs[i].index_fgc = symOpID;
        symOpID = xstr.pgroupk_xtal.size();
      }
    }
  }

  // ****************** BEGIN Printing block *************************
  oss << aurostd::PaddedPRE("# Index", 6, " ") << "                             " << aurostd::PaddedPRE("Q-POINTS (1/nm)", 14, " ") << std::endl;
  for (uint i = 0; i < QPs.size(); i++) {
    oss << setw(8) << i << "    ";
    oss << setw(20) << setprecision(10) << std::scientific << QPs[i].q[1];
    oss << setw(20) << setprecision(10) << std::scientific << QPs[i].q[2];
    oss << setw(20) << setprecision(10) << std::scientific << QPs[i].q[3];
    oss << std::endl;
  }
  // ****************** END  Printing block *************************

  // Change of Units [1/A] --> [1/nm]. Numerical Issues
  for (uint i = 0; i < QPs.size(); i++) {
    QPs[i].q = QPs[i].q * 10;
  }

  // ****************** BEGIN Printing block *************************
  oss << "[APL_THERMAL_QPOINTS]STOP" << std::endl;
  oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
  string ofiletemppropsname = "./aflow.apl.thermal_QPs.out";
  if (!aurostd::stringstream2file(oss, ofiletemppropsname, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_QPs.out ");
  }
  aurostd::StringstreamClean(oss);
  stringstream oss2;
  oss2 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  oss2 << "[APL_THERMAL_IRREDUCIBLE_QPOINTS]START" << std::endl;
  oss2 << aurostd::PaddedPRE("# Index", 6, " ") << "    " << aurostd::PaddedPRE("Degeneracy", 4, " ") << "                        " << aurostd::PaddedPRE("Q-POINTS (1/nm)", 14, " ") << std::endl;
  // ****************** END Printing block *************************

  //Creating inequivalent list iQPs
  for (uint i = 0; i < QPs.size(); i++) {
    _iQPs kk2;
    if (!QPs[i].q_is_equivalent) {
      kk2.q = QPs[i].q;
      kk2.mapq = QPs[i].mapq;
      for (uint j = i; j < QPs.size(); j++) {
        if (QPs[j].index_ineq == (int)i) {
          kk2.list.push_back(j);
        }
      }

      // Look for invariance symmetry operations
      xvector<double> newkpoint(3);
      for (uint symOpID = 0; symOpID < xstr.pgroupk_xtal.size(); symOpID++) {
        newkpoint = xstr.pgroupk_xtal[symOpID].Uf * QPs[i].q_dir;
        xvector<double> fdiff = newkpoint - QPs[i].q_dir;
        SYM::PBC(fdiff);
        if (aurostd::modulus(fdiff) < _AFLOW_APL_EPS_) {
          kk2.lsym.push_back(symOpID);
        }
      }
      iQPs.push_back(kk2);

      // ****************** BEGIN Printing block *************************
      oss2 << setw(6) << i << "    ";
      oss2 << setw(6) << kk2.list.size() << "    ";
      oss2 << setw(20) << setprecision(10) << std::scientific << kk2.q[1];
      oss2 << setw(20) << setprecision(10) << std::scientific << kk2.q[2];
      oss2 << setw(20) << setprecision(10) << std::scientific << kk2.q[3];
      oss2 << std::endl;
      // ****************** END Printing block *************************
    }
  }

  // ****************** BEGIN Printing block *************************
  oss2 << "[APL_THERMAL_IRREDUCIBLE_QPOINTS]STOP" << std::endl;
  oss2 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  string ofiletemppropsname2 = "./aflow.apl.thermal_iQPs.out";
  if (!aurostd::stringstream2file(oss2, ofiletemppropsname2, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_iQPs.out");
  }
  aurostd::StringstreamClean(oss2);
  // ****************** END Printing block *************************
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void ThermalConductivityCalculator::calcFreq(int nx, int ny, int nz) {
  // Get q-points for which to calculate the frequencies and store grid
  grid[1] = nx;
  grid[2] = ny;
  grid[3] = nz;
  NPQs = (double)nx * ny * nz;
  // Maybe there was some error and the list of q-points is empty, hence bye-bye...
  if (QPs.empty())
    throw apl::APLRuntimeError("There are no points for calculation.");

  // Preparing storing for eigenvectors and gradient of Dynamical Matrix
  xvector<double> pp;
  for (uint i = 0; i < QPs.size(); i++) {
    shortest.push_back(pp);
    vector<xmatrix<xcomplex<double>>> kkk;
    xmatrix<xcomplex<double>> kk(nBranches, nBranches);
    for (uint l = 0; l < 3; l++) {
      for (uint j = 1; j < nBranches + 1; j++) {
        for (uint k = 1; k < nBranches + 1; k++) {
          kk(j, k).re = 0.0;
          kk(j, k).im = 0.0;
        }
      }
      kkk.push_back(kk);
    }
    EigenVec.push_back(kk);
    dDynM.push_back(kkk);
  }

  // Use the 1st BZ image of each q point to improve the behavior of the non-analytic correction.
  double tmp1, tmp2;
  for (uint i = 0; i < QPs.size(); i++) {
    shortest[i] = QPs[i].q;
    tmp1 = aurostd::modulus(shortest[i]);
    for (int o = -2; o <= 2; o++) {
      for (int p = -2; p <= 2; p++) {
        for (int q = -2; q <= 2; q++) {
          xvector<double> pp = QPs[i].q + 10.0 * (o * _klattice(1) + p * _klattice(2) + q * _klattice(3));
          tmp2 = aurostd::modulus(pp);
          if (tmp2 < tmp1) {
            tmp1 = tmp2;
            shortest[i] = pp;
          }
        }
      }
    }
  }

// Compute frequencies for each q-point

#ifdef AFLOW_APL_MULTITHREADS_ENABLE  // If MPI....

  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
  int qpointsPerCPU = QPs.size() / ncpus;

  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating frequencies and eigenvectors for THERMALGRIDD");
  else
    _logger.initProgressBar("Calculating frequencies and eigenvectors for THERMALGRIDDD (" + stringify(ncpus) + " threads)");

  // Prepare storage
  freq.clear();
  xvector<double> zero(nBranches);
  for (uint i = 0; i < QPs.size(); i++)
    freqs.push_back(zero);

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread*> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > QPs.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < QPs.size())))
      endIndex = QPs.size();
    threads.push_back(new std::thread(&ThermalConductivityCalculator::calculateFreqMatMPI, this, startIndex, endIndex));
  }

  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }
  threads.clear();

  // Done
  _logger.finishProgressBar();

#else

  _logger.initProgressBar("Calculating frequencies and eigenvectors for THERMALGRID");
  xmatrix<xcomplex<double>> eigenVec(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double>> dummyyy(nBranches, nBranches, 1, 1);
  vector<xmatrix<xcomplex<double>>> ddynm;
  for (uint iqp = 0; iqp < 3; iqp++) {
    ddynm.push_back(eigenVec);
  }
  for (uint iqp = 0; iqp < QPs.size(); iqp++) {
    _logger.updateProgressBar(iqp, QPs.size());
    freqs.push_back(_pc.getFrequencyAAPL(QPs[iqp].q / 10.0, eigenVec, ddynm));
    EigenVec[iqp] = eigenVec;

    dDynM[iqp] = ddynm;
  }
  _logger.finishProgressBar();

#endif

  // ****************** BEGIN Printing block *************************
  stringstream oss;
  oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
  oss << "[APL_THERMAL_QPOINTS]START" << std::endl;
  oss << aurostd::PaddedPRE("# Index", 6, " ") << "                           " << aurostd::PaddedPRE("Angular frequencies (rad/ps)", 14, " ") << std::endl;
  // ****************** END Printing block *************************

  // Freq in a vector instead of xvector
  for (uint q = 0; q < QPs.size(); q++) {  //Loop over q points
    vector<double> kk;
    oss << setw(6) << q << "    ";
    for (uint br = 1; br <= nBranches; br++) {  // Loop over branches
      oss << setw(20) << setprecision(10) << std::scientific << freqs[q][br];
      double pp = freqs[q][br];
      kk.push_back(pp);
    }
    oss << std::endl;
    freq.push_back(kk);
  }

  // ****************** BEGIN Printing block *************************
  oss << "[APL_THERMAL_FREQUENCIES]STOP" << std::endl;
  oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
  string ofiletemppropsname = "./aflow.apl.thermal_omega.out";
  if (!aurostd::stringstream2file(oss, ofiletemppropsname, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_omega.out ");
  }
  // ****************** END Printing block *************************
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::calculateFreqMatMPI(int startIndex, int endIndex) {
  xmatrix<xcomplex<double>> eigenVec(nBranches, nBranches, 1, 1);
  vector<xmatrix<xcomplex<double>>> ddynm;
  for (uint iqp = 0; iqp < 3; iqp++) {
    ddynm.push_back(eigenVec);
  }
  for (int iqp = startIndex; iqp < endIndex; iqp++) {
    _logger.updateProgressBar(iqp, QPs.size());
    freqs[iqp] = _pc.getFrequencyAAPL(QPs[iqp].q / 10.0, eigenVec, ddynm);
    EigenVec[iqp] = eigenVec;
    dDynM[iqp] = ddynm;
    //std::this_thread::yield();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::getGroupVelocity() {
  xvector<double> testVec(3);

  // Preaparing matrix
  xmatrix<xcomplex<double>> eigen(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double>> ddynX(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double>> ddynY(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double>> ddynZ(nBranches, nBranches, 1, 1);

  for (uint q = 0; q < QPs.size(); q++) {  //Loop over q points

    vector<xvector<double>> branch;

    // Matrix that we need
    for (uint i = 1; i <= nBranches; i++) {
      for (uint j = 1; j <= nBranches; j++) {
        eigen[i][j] = EigenVec[q][j][i];
      }
    }

    ddynX = dDynM[q][0];
    ddynY = dDynM[q][1];
    ddynZ = dDynM[q][2];

    for (uint br = 1; br <= nBranches; br++) {  // Loop over branches
      xvector<double> gvel(3, 1);
      if (freq[q][br - 1] < _AFLOW_APL_EPS_) {
        gvel[1] = 0;
        gvel[2] = 0;
        gvel[3] = 0;
      } else {
        xvector<xcomplex<double>> _eigen(nBranches, 1);
        xvector<xcomplex<double>> _eigen_c(nBranches, 1);
        xvector<xcomplex<double>> dummyX(nBranches, 1);
        xvector<xcomplex<double>> dummyY(nBranches, 1);
        xvector<xcomplex<double>> dummyZ(nBranches, 1);

        // Prepare a column of eigenvectors

        for (uint j = 1; j <= nBranches; j++) {
          _eigen(j) = eigen(br, j);
          _eigen_c[j].re = eigen[br][j].re;
          _eigen_c[j].im = (-1) * eigen[br][j].im;
        }

        for (uint i = 1; i <= nBranches; i++) {
          for (uint j = 1; j <= nBranches; j++) {
            dummyX[i] += ddynX[i][j] * _eigen[j];
            dummyY[i] += ddynY[i][j] * _eigen[j];
            dummyZ[i] += ddynZ[i][j] * _eigen[j];
          }
        }
        xcomplex<double> xx = _eigen_c * dummyX;
        xcomplex<double> yy = _eigen_c * dummyY;
        xcomplex<double> zz = _eigen_c * dummyZ;

        // Group Velocities in Angstrom/s
        gvel[1] = (xx.re / (2 * freq[q][br - 1])) * 0.1;
        gvel[2] = (yy.re / (2 * freq[q][br - 1])) * 0.1;
        gvel[3] = (zz.re / (2 * freq[q][br - 1])) * 0.1;
      }
      branch.push_back(gvel);
    }
    groupVel.push_back(branch);
  }

  // Avoid symmtrizatiion 
  Print_GV();  // Printing group velocities
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::Print_GV() {
  stringstream oss;
  oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
  oss << "[APL_THERMAL_GROUP_VELOCITIES]START" << std::endl;
  oss << "# Frequencies and Group Velocities calculated by Aflow using Hellmann-Feynman theorem" << std::endl;
  oss << aurostd::PaddedPRE("# Index", 6, " ") << "                           " << aurostd::PaddedPRE("Group Velocoties (km/s (or nm THz))", 14, " ") << std::endl;
  for (uint q = 0; q < QPs.size(); q++) {  //Loop over q points
    oss << setw(6) << q << "    ";
    for (uint br = 1; br <= nBranches; br++) {  // Loop over branches
      for (uint i = 1; i <= 3; i++) {           // Loop over x,y,z
        oss << setw(20) << setprecision(10) << std::scientific << groupVel[q][br - 1][i];
      }
      oss << "     ";
    }
    oss << std::endl;
  }
  oss << "[APL_THERMAL_GROUP_VELOCITIES]STOP" << std::endl;
  oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
  string ofiletemppropsname = "./aflow.apl.thermal_gvel.out";
  if (!aurostd::stringstream2file(oss, ofiletemppropsname, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_gvel.out ");
  }
  aurostd::StringstreamClean(oss);
}
////////////////////////////////////////////////////////////////////////////////////////////////////


void ThermalConductivityCalculator::getScattering() {
  for (uint i = 0; i < _pa.tri_index.size(); i++) {
    //uint ipc1 = _pa.tri_index[i][1];
    //uint isc1 = _sc.pc2scMap(ipc1); // Equivalent atom in sc1
    uint isc2 = _pa.tri_index[i][2];
    uint ipc2 = _sc.sc2pcMap(isc2);  //Equivalent atom ins pc2
    uint isc3 = _pa.tri_index[i][3];
    uint ipc3 = _sc.sc2pcMap(isc3);  //Equivalent atom ins pc2
    xcomplex<double> phase1;
    xcomplex<double> phase2;
    vector<xvector<double>> kk;
    xvector<double> kk2;
    uint at2_sc = _sc.pc2scMap(ipc2);
    uint at3_sc = _sc.pc2scMap(ipc3);
    kk2 = _pa.tri_vec[i][1] - _sc.getStructure().atoms[at2_sc].cpos;
    kk.push_back(kk2);
    kk2 = _pa.tri_vec[i][2] - _sc.getStructure().atoms[at3_sc].cpos;
    kk.push_back(kk2);
    tri_vec.push_back(kk);
    kk.clear();
  }

#ifdef AFLOW_APL_MULTITHREADS_ENABLE

  // Get the number of CPUS
  int ncpus1; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus1);  //CO 180214
  if (ncpus1 < 1) ncpus1 = 1;

  // Show info
  if (ncpus1 == 1)
    _logger.initProgressBar("Calculating number of Scattering  processes");
  else
    _logger.initProgressBar("Calculating number of Scattering processes (" + stringify(ncpus1) + " threads)");

  // Prepare storage

  // Distribute the calculation
  int startIndex1;
  std::vector<std::thread*> threads1;
  for (startIndex1 = 0; startIndex1 < 2; startIndex1++) {
    threads1.push_back(new std::thread(&ThermalConductivityCalculator::calculateScatProc_MPI, this, startIndex1));
  }

  // Wait to finish all threads here!
  for (uint i = 0; i < threads1.size(); i++) {
    threads1[i]->join();
    delete threads1[i];
  }
  threads1.clear();

  // Done
  _logger.finishProgressBar();

#else

  //NEW Storage
  vector<double> tt;
  vector<double> mm;
  vector<xvector<int>> LPq;
  vector<xvector<int>> LPb;
  xvector<int> inQ;
  xvector<int> inB;
  xcomplex<double> iONE;
  iONE.re = 0.0;
  iONE.im = 1.0;

  //Looking for scattering processes
  int N_plus = 0;
  int N_minus = 0;
  _logger << "Entering Storing Plus " << apl::endl;
  for (uint iq = 0; iq < iQPs.size(); iq++) {
    uint q = iQPs[iq].list[0];
    int q1x = iQPs[iq].mapq[1];
    int q1y = iQPs[iq].mapq[2];
    int q1z = iQPs[iq].mapq[3];
    for (uint br = 0; br < nBranches; br++) {
      if (freq[q][br] != 0.0) {
        for (int q2x = 0; q2x < grid(1); q2x++) {
          for (int q2y = 0; q2y < grid(2); q2y++) {
            for (int q2z = 0; q2z < grid(3); q2z++) {
              int q2 = mapQ[q2x][q2y][q2z];

              int q3x = q1x + q2x;
              if (q3x >= grid(1)) { q3x = q3x - grid(1); }
              int q3y = q1y + q2y;
              if (q3y >= grid(2)) { q3y = q3y - grid(2); }
              int q3z = q1z + q2z;
              if (q3z >= grid(3)) { q3z = q3z - grid(3); }
              int q3 = mapQ[q3x][q3y][q3z];
              for (uint br2 = 0; br2 < nBranches; br2++) {
                // Q1+Q2=Q3
                for (uint br3 = 0; br3 < nBranches; br3++) {
                  if ((freq[q3][br3] != 0.0) && (freq[q2][br2] != 0.0)) {
                    double sigma = getSigma(q2, br2, q3, br3);
                    if (std::abs(freq[q][br] + freq[q2][br2] - freq[q3][br3]) < 2 * sigma) {
                      N_plus = N_plus + 1;
                      inQ[1] = q;
                      inQ[2] = q2;
                      inQ[3] = q3;
                      LPq.push_back(inQ);
                      inQ.clear();
                      inB[1] = br;
                      inB[2] = br2;
                      inB[3] = br3;
                      LPb.push_back(inB);
                      inB.clear();
                      tt.push_back(0.0);
                      mm.push_back(sigma);
                    }
                  }
                }
              }
            }
          }
        }
      }
      //Storing Plus
      ScatPlus.push_back(tt);
      tt.clear();
      ScatPlus_s.push_back(mm);
      mm.clear();
      ScatPlus_q.push_back(LPq);
      LPq.clear();
      ScatPlus_b.push_back(LPb);
      LPb.clear();
    }
  }
  // Q1-Q2=Q3
  for (uint iq = 0; iq < iQPs.size(); iq++) {
    uint q = iQPs[iq].list[0];
    int q1x = iQPs[iq].mapq[1];
    int q1y = iQPs[iq].mapq[2];
    int q1z = iQPs[iq].mapq[3];
    for (uint br = 0; br < nBranches; br++) {
      if (freq[q][br] != 0.0) {
        for (int q2x = 0; q2x < grid(1); q2x++) {
          for (int q2y = 0; q2y < grid(2); q2y++) {
            for (int q2z = 0; q2z < grid(3); q2z++) {
              int q2 = mapQ[q2x][q2y][q2z];
              int q3x = q1x - q2x;
              if (q3x < 0) { q3x = q3x + grid(1); }
              int q3y = q1y - q2y;
              if (q3y < 0) { q3y = q3y + grid(2); }
              int q3z = q1z - q2z;
              if (q3z < 0) { q3z = q3z + grid(3); }
              int q3 = mapQ[q3x][q3y][q3z];
              for (uint br2 = 0; br2 < nBranches; br2++) {
                for (uint br3 = 0; br3 < nBranches; br3++) {
                  if ((freq[q3][br3] != 0.0) && (freq[q2][br2] != 0.0)) {
                    double sigma = getSigma(q2, br2, q3, br3);
                    if (std::abs(freq[q][br] - freq[q2][br2] - freq[q3][br3]) < 2 * sigma) {
                      /// BEGIN SCATERING MATRIX CALCULATION ///
                      N_minus = N_minus + 1;
                      inQ[1] = q;
                      inQ[2] = q2;
                      inQ[3] = q3;
                      LPq.push_back(inQ);
                      inQ.clear();
                      inB[1] = br;
                      inB[2] = br2;
                      inB[3] = br3;
                      LPb.push_back(inB);
                      inB.clear();
                      tt.push_back(0.0);
                      mm.push_back(sigma);
                    }
                  }
                }
              }
            }
          }
        }
      }
      ScatMinus.push_back(tt);
      tt.clear();
      ScatMinus_s.push_back(mm);
      mm.clear();
      ScatMinus_q.push_back(LPq);
      LPq.clear();
      ScatMinus_b.push_back(LPb);
      LPb.clear();
    }
  }
#endif

//Computing scattering processes

#ifdef AFLOW_APL_MULTITHREADS_ENABLE

  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
  uint scat = iQPs.size() * nBranches;
  int qpointsPerCPU = scat / ncpus;

  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating Scattering processes");
  else
    _logger.initProgressBar("Calculating Scattering processes (" + stringify(ncpus) + " threads)");

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread*> threads;
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > scat) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < scat)))
      endIndex = scat;
    threads.push_back(new std::thread(&ThermalConductivityCalculator::calculateScat_MPI, this, startIndex, endIndex));
  }

  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }
  threads.clear();

  // Done
  _logger.finishProgressBar();

#else

  // Q1+Q2=Q3
  for (uint iq = 0; iq < iQPs.size(); iq++) {
    uint q = iQPs[iq].list[0];
    for (uint br = 0; br < nBranches; br++) {
      if (freq[q][br] != 0.0) {
        for (uint in = 0; in < ScatPlus[(iq * nBranches) + br].size(); in++) {  //
          int q2 = ScatPlus_q[(iq * nBranches) + br][in][2];
          int q3 = ScatPlus_q[(iq * nBranches) + br][in][3];
          int br2 = ScatPlus_b[(iq * nBranches) + br][in][2];
          int br3 = ScatPlus_b[(iq * nBranches) + br][in][3];
          xvector<double> kpoint1 = QPs[q2].q * 0.1;
          xvector<double> kpoint2 = QPs[q3].q * 0.1;
          double sigma = ScatPlus_s[(iq * nBranches) + br][in];
          /// BEGIN SCATERING MATRIX CALCULATION ///
          xcomplex<double> Vp;
          Vp.re = 0.0;
          Vp.im = 0.0;
          for (uint i = 0; i < _pa.tri_index.size(); i++) {
            uint ipc1 = _pa.tri_index[i][1];
            uint isc1 = _sc.pc2scMap(ipc1);  // Equivalent atom in sc1
            uint isc2 = _pa.tri_index[i][2];
            uint ipc2 = _sc.sc2pcMap(isc2);  //Equivalent atom ins pc2
            uint isc3 = _pa.tri_index[i][3];
            uint ipc3 = _sc.sc2pcMap(isc3);  //Equivalent atom ins pc3
            xcomplex<double> phase1;
            xcomplex<double> phase2;
            phase1 = exp(iONE * scalar_product(kpoint1, tri_vec[i][0]));
            phase2 = exp(-iONE * scalar_product(kpoint2, tri_vec[i][1]));
            xcomplex<double> prefactor = (1 / sqrt(_sc.getAtomMass(isc1) * _sc.getAtomMass(isc2) * _sc.getAtomMass(isc3))) * phase1 * phase2;
            xcomplex<double> Vp0;
            Vp0.re = 0;
            Vp0.im = 0;
            for (uint alfa = 0; alfa < 3; alfa++) {  // Loop over three directions
              for (uint beta = 0; beta < 3; beta++) {
                for (uint gam = 1; gam < 4; gam++) {
                  Vp0 = Vp0 + (_pc.getElementTensor(ipc1, isc2, isc3, alfa, beta, gam) * EigenVec[q][(ipc1 * 3) + alfa + 1][br + 1] * EigenVec[q2][(ipc2 * 3) + beta + 1][br2 + 1] * conj(EigenVec[q3][(ipc3 * 3) + gam][br3 + 1]));
                }
              }
            }
            Vp = Vp + prefactor * Vp0;
          }
          /// END SCATERING MATRIX CALCULATION ///
          ScatPlus[(iq * nBranches) + br][in] = (hbarp * NPI / 4.0) * exp(-(freq[q][br] + freq[q2][br2] - freq[q3][br3]) * (freq[q][br] + freq[q2][br2] - freq[q3][br3]) / (sigma * sigma)) / sqrt(NPI) / sigma / (freq[q][br] * freq[q2][br2] * freq[q3][br3]) * ((Vp.re * Vp.re) + (Vp.im * Vp.im));
          ScatPlus[(iq * nBranches) + br][in] = ScatPlus[(iq * nBranches) + br][in] * 5.60626442E08 / QPs.size();  //THz
        }
      }
    }
  }
  // Q1-Q2=Q3
  for (uint iq = 0; iq < iQPs.size(); iq++) {
    uint q = iQPs[iq].list[0];
    for (uint br = 0; br < nBranches; br++) {
      if (freq[q][br] != 0.0) {
        for (uint in = 0; in < ScatMinus[(iq * nBranches) + br].size(); in++) {  //
          int q2 = ScatMinus_q[(iq * nBranches) + br][in][2];
          int q3 = ScatMinus_q[(iq * nBranches) + br][in][3];
          int br2 = ScatMinus_b[(iq * nBranches) + br][in][2];
          int br3 = ScatMinus_b[(iq * nBranches) + br][in][3];
          xvector<double> kpoint1 = QPs[q2].q * 0.1;
          xvector<double> kpoint2 = QPs[q3].q * 0.1;
          double sigma = ScatMinus_s[(iq * nBranches) + br][in];
          /// BEGIN SCATERING MATRIX CALCULATION ///
          xcomplex<double> Vp;
          Vp.re = 0;
          Vp.im = 0;
          for (uint i = 0; i < _pa.tri_index.size(); i++) {
            uint ipc1 = _pa.tri_index[i][1];
            uint isc1 = _sc.pc2scMap(ipc1);  // Equivalent atom in sc1
            uint isc2 = _pa.tri_index[i][2];
            uint ipc2 = _sc.sc2pcMap(isc2);  //Equivalent atom ins pc2
            uint isc3 = _pa.tri_index[i][3];
            uint ipc3 = _sc.sc2pcMap(isc3);  //Equivalent atom ins pc2
            xcomplex<double> phase1;
            xcomplex<double> phase2;
            phase1 = exp(-iONE * scalar_product(kpoint1, tri_vec[i][0]));
            phase2 = exp(-iONE * scalar_product(kpoint2, tri_vec[i][1]));
            xcomplex<double> prefactor = (1 / sqrt(_sc.getAtomMass(isc1) * _sc.getAtomMass(isc2) * _sc.getAtomMass(isc3))) * phase1 * phase2;
            xcomplex<double> Vp0;
            Vp0.re = 0;
            Vp0.im = 0;
            for (uint alfa = 0; alfa < 3; alfa++) {  // Loop over three directions
              for (uint beta = 0; beta < 3; beta++) {
                for (uint gam = 1; gam < 4; gam++) {
                  Vp0 = Vp0 + (_pc.getElementTensor(ipc1, isc2, isc3, alfa, beta, gam) * EigenVec[q][(ipc1 * 3) + alfa + 1][br + 1] * conj(EigenVec[q2][(ipc2 * 3) + beta + 1][br2 + 1]) * conj(EigenVec[q3][(ipc3 * 3) + gam][br3 + 1]));
                }
              }
            }
            Vp = Vp + prefactor * Vp0;
          }
          /// END SCATERING MATRIX CALCULATION ///
          ScatMinus[(iq * nBranches) + br][in] = (hbarp * NPI / 4.0) * exp(-(freq[q][br] - freq[q2][br2] - freq[q3][br3]) * (freq[q][br] - freq[q2][br2] - freq[q3][br3]) / (sigma * sigma)) / sqrt(NPI) / sigma / (freq[q][br] * freq[q2][br2] * freq[q3][br3]) * ((Vp.re * Vp.re) + (Vp.im * Vp.im));
          ScatMinus[(iq * nBranches) + br][in] = ScatMinus[(iq * nBranches) + br][in] * 5.60626442E08 / QPs.size();  //THz
        }
      }
    }
  }

#endif

  _pa.clear();
}
/////////////////////////////////////////////////////////////////////////////

int ThermalConductivityCalculator::get3rdQ_P(uint q1, uint q2) {
  //Sumation of q1 and q2

  xvector<double> kk = (QPs[q1].q_dir + QPs[q2].q_dir);
  kk = BringInCell(kk);  // Obtaining fractional position
  xvector<double> qdummy(3);
  for (uint s = 0; s < QPs.size(); s++) {
    //Compare vectors
    xvector<double> fdiff = kk - QPs[s].q_dir;
    SYM::PBC(fdiff);
    if (aurostd::modulus(fdiff) < _AFLOW_APL_EPS_) {
      return s;
    }
  }

  throw APLRuntimeError("MeshTC::We dind't find q vector for mommentum conservation. How is that possible!!!!!");

  return 0;
}
/////////////////////////////////////////////////////////////////////////////

int ThermalConductivityCalculator::get3rdQ_M(uint q1, uint q2) {
  //Substraction of q1 and q2

  xvector<double> kk = (QPs[q1].q_dir - QPs[q2].q_dir);
  kk = BringInCell(kk);  // Obtaining fractional position
  xvector<double> qdummy(3);
  for (uint s = 0; s < QPs.size(); s++) {
    //Compare vectors
    xvector<double> fdiff = kk - QPs[s].q_dir;
    SYM::PBC(fdiff);
    if (aurostd::modulus(fdiff) < _AFLOW_APL_EPS_) {
      return s;
    }
  }
  throw APLRuntimeError("MeshTC::We dind't find q vector for mommentum conservation. How is that possible!!!!!");

  return 0;
}

//////////////////////////////////////////////////////////////////////////////

double ThermalConductivityCalculator::getSigma(uint q2, uint br2, uint q3, uint br3) {
  // See J. Carrete Paper (Computer Physics Communications)
  double kk1 = 0.0;
  for (uint i = 1; i < 4; i++) {
    kk1 = kk1 + pow(aurostd::scalar_product((groupVel[q2][br2] - groupVel[q3][br3]), (_klattice(i) * 10.0)) / (grid[i]), 2);
  }
  double sigmaW = sqrt(kk1 / 6.0);
  return sigmaW * 0.1;  // Scalebroad = 0.1
}

//////////////////////////////////////////////////////////////////////////////


void ThermalConductivityCalculator::Calculator(bool RTA, bool ISO, bool CUML, bool BOUND, double NANO, double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP) {
  //  Computing Fn and kappa

  _logger << "Computing kappa " << apl::endl;

  Vol = _sc._pcStructure.Volume();  //Primitive cell Volume
  double temp = USER_TP_TSTART;

  // ****************** BEGIN Printing block *************************
  stringstream oss, oss2, oss3, oss4, oss5;
  oss4 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  oss4 << "[APL_THERMAL_KAPPA]START" << std::endl;
  oss4 << aurostd::PaddedPRE("# Temp (K)", 10, " ") << "                                                                 " << aurostd::PaddedPRE("Thermal Conductivity (W/m*K)", 14, " ") << std::endl;
  oss4 << "#           " << aurostd::PaddedPRE("k(1,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(1,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(1,3)", 10, " ") << "          ";
  oss4 << aurostd::PaddedPRE("k(2,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(2,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(2,3)", 10, " ") << "          ";
  oss4 << aurostd::PaddedPRE("k(3,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(3,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(3,3)", 10, " ") << "          ";
  oss4 << std::endl;
  // ****************** END Printing block *************************

  //  Computing Wiso
  _logger << "Computing ISO " << apl::endl;
  if (ISO) { getWiso(); }
  _logger << "Computing ISO END" << apl::endl;

  while (temp <= USER_TP_TEND) {  //Temperature loop

    // ****************** BEGIN Printing block *************************
    oss2 << "[AFLOW] **************************************************************************************************************************" << std::endl;
    oss2 << "[APL_THERMAL_W= " << setw(6) << setprecision(2) << fixed << temp << " K]START" << std::endl;
    oss2 << aurostd::PaddedPRE("# Index", 6, " ") << "                           " << aurostd::PaddedPRE("Scattering rates (1/ps)", 14, " ") << std::endl;
    oss5 << "[AFLOW] **************************************************************************************************************************" << std::endl;
    oss5 << "[APL_THERMAL_Wtotal= " << setw(6) << setprecision(2) << fixed << temp << " K]START" << std::endl;
    oss5 << aurostd::PaddedPRE("# Index", 6, " ") << "                           " << aurostd::PaddedPRE("Scattering rates (1/ps)", 14, " ") << std::endl;
    oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
    oss << "[APL_THERMAL_TAU_ZERO= " << setw(6) << setprecision(2) << fixed << temp << " K]START" << std::endl;
    oss << aurostd::PaddedPRE("# Index", 6, " ") << "                           " << aurostd::PaddedPRE("Scattering time (ps)", 14, " ") << std::endl;
    // ****************** END Printing block *************************

    _logger << "Computing RTA " << apl::endl;
    getTao_RTA(temp, ISO, BOUND, NANO);  // Compute tau_zero and Wq (temperature dependent)
    _logger << "Computing RTA END" << apl::endl;

    // ****************** BEGIN Printing block *************************
    for (uint iq = 0; iq < iQPs.size(); iq++) {  //Loop over q points
      uint q = iQPs[iq].list[0];
      oss << setw(6) << q << "    ";
      oss2 << setw(6) << q << "    ";
      oss5 << setw(6) << q << "    ";
      for (uint br = 0; br < nBranches; br++) {  // Loop over branches
        oss << setw(20) << setprecision(10) << std::scientific << tau_RTA[iq][br];
        oss2 << setw(20) << setprecision(10) << std::scientific << Wq[iq][br];
        oss5 << setw(20) << setprecision(10) << std::scientific << Wtotal[iq][br];
      }
      oss << std::endl;
      oss2 << std::endl;
      oss5 << std::endl;
    }
    oss << "[APL_THERMAL_TAU_ZERO= " << setw(6) << setprecision(2) << fixed << temp << " K]STOP" << std::endl;
    oss << "[AFLOW] **************************************************************************************************************************" << std::endl;
    oss2 << "[APL_THERMAL_W= " << setw(6) << setprecision(2) << fixed << temp << " K]STOP" << std::endl;
    oss2 << "[AFLOW] **************************************************************************************************************************" << std::endl;
    oss5 << "[APL_THERMAL_Wtotal= " << setw(6) << setprecision(2) << fixed << temp << " K]STOP" << std::endl;
    oss5 << "[AFLOW] **************************************************************************************************************************" << std::endl;
    // ****************** END Printing block *************************

    _logger << "Computing iteration0 " << apl::endl;
    iteration0();  // Get Fn0 (RTA solution for BTE)
    _logger << "Computing iteration0 END " << apl::endl;
    xmatrix<double> kappa;

    // ****************** BEGIN Printing block *************************
    oss3 << "[AFLOW] **************************************************************************************************************************" << std::endl;
    oss3 << "[APL_THERMAL_KAPPA_SCF= " << setw(6) << setprecision(2) << fixed << temp << " K]START" << std::endl;
    oss3 << aurostd::PaddedPRE("# Iteration", 10, " ") << "                                                               " << aurostd::PaddedPRE("Thermal Conductivity (W/m*K)", 14, " ") << std::endl;
    oss3 << "#         " << aurostd::PaddedPRE("k(1,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(1,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(1,3)", 10, " ") << "          ";
    oss3 << aurostd::PaddedPRE("k(2,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(2,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(2,3)", 10, " ") << "          ";
    oss3 << aurostd::PaddedPRE("k(3,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(3,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(3,3)", 10, " ") << "          ";
    oss3 << std::endl;
    // ****************** END Printing block *************************

    _logger << "Computing getKappa0  " << apl::endl;
    kappa = getKappa(temp);  //Get kappa using Fn0
    _logger << "Computing getKappa0 END " << apl::endl;
    if (CUML) {
      getCumKappa(temp);
    }

    // ****************** BEGIN Printing block *************************
    oss3 << setw(6) << (int)0;
    for (uint i = 1; i < 4; i++) {
      for (uint j = 1; j < 4; j++) {
        oss3 << setw(20) << setprecision(10) << std::scientific << kappa[i][j];
      }
    }
    oss3 << std::endl;
    // ****************** END Printing block *************************

    double dmax = 1.0;
    int counter = 1;

    if (!RTA) {
      while (dmax > 0.0005) {  // SCF for Full solution BTE
        _logger << "Computing iteration   " << apl::endl;
        iteration(temp);  // Get Fn
        double kappa_old = aurostd::max(kappa);
        kappa = getKappa(temp);  // Get kappa
        double kap = aurostd::max(kappa);
        dmax = abs((kappa_old - kap) * 100 / kap);

        // ****************** BEGIN Printing block *************************
        oss3 << setw(6) << counter;
        for (uint i = 1; i < 4; i++) {
          for (uint j = 1; j < 4; j++) {
            oss3 << setw(20) << setprecision(10) << std::scientific << kappa[i][j];
          }
        }
        oss3 << std::endl;
        // ****************** END Printing block *************************

        counter = counter + 1;
      }
    }

    // ****************** BEGIN Printing block *************************
    oss4 << setw(8) << setprecision(2) << fixed << temp;
    for (uint i = 1; i < 4; i++) {
      for (uint j = 1; j < 4; j++) {
        oss4 << setw(20) << setprecision(10) << std::scientific << kappa[i][j];
      }
    }
    oss4 << std::endl;
    oss3 << "[APL_THERMAL_KAPPA= " << setw(6) << setprecision(2) << fixed << temp << " K]STOP" << std::endl;
    oss3 << "[AFLOW] **************************************************************************************************************************" << std::endl;
    // ****************** END Printing block *************************
    temp = temp + USER_TP_TSTEP;  // Increase temp
  }

  // ****************** BEGIN Printing block *************************
  oss4 << "[APL_THERMAL_KAPPA]STOP" << std::endl;
  oss4 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  string ofiletemppropsname2 = "./aflow.apl.thermal_w.out";
  if (!aurostd::stringstream2file(oss2, ofiletemppropsname2, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_w.out ");
  }
  string ofiletemppropsname5 = "./aflow.apl.thermal_w_total.out";
  if (!aurostd::stringstream2file(oss5, ofiletemppropsname5, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_wtotal.out ");
  }
  string ofiletemppropsname = "./aflow.apl.thermal_tau.out";
  if (!aurostd::stringstream2file(oss, ofiletemppropsname, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_tau.out ");
  }
  string ofiletemppropsname3 = "./aflow.apl.thermal_kappa_scf.out";
  if (!aurostd::stringstream2file(oss3, ofiletemppropsname3, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_kappa_scf.out ");
  }
  string ofiletemppropsname4 = "./aflow.apl.thermal_kappa.out";
  if (!aurostd::stringstream2file(oss4, ofiletemppropsname4, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_kappa.out ");
  }
  // ****************** END Printing block *************************

  plotKappa();
}
//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::getWiso() {
  stringstream oss2;

  for (uint iq = 0; iq < iQPs.size(); iq++) {
    vector<double> kk;
    for (uint br = 0; br < nBranches; br++) {
      kk.push_back(0.0);
    }
    Wiso.push_back(kk);
  }

  //Computing sigmas following ShengBTE recipe
  vector<vector<double>> sigma;
  for (uint q = 0; q < QPs.size(); q++) {
    vector<double> kk;
    for (uint br = 0; br < nBranches; br++) {
      kk.push_back(getSigma(q, br));
    }
    sigma.push_back(kk);
  }

  xvector<double> list(QPs.size() * nBranches);
  for (uint q = 0; q < QPs.size(); q++) {
    for (uint br = 0; br < nBranches; br++) {
      int pos = (nBranches * q) + br + 1;
      if (sigma[q][br] == 0.0) {
        list(pos) = -1E33;
      } else {
        list(pos) = log(sigma[q][br]);
      }
    }
  }

  list = aurostd::sort(list);

  for (uint q = 1; q < (QPs.size() * nBranches) + 1; q++) {
  }

  double per25 = list(25 * nBranches * QPs.size() / 100);
  double per75 = list(75 * nBranches * QPs.size() / 100);
  double delta = per75 - per25;
  double lbound = exp(per25 - 1.5 * delta);
  for (uint q = 0; q < QPs.size(); q++) {
    for (uint br = 0; br < nBranches; br++) {
      sigma[q][br] = std::max(sigma[q][br], lbound);
    }
  }

  for (uint q = 0; q < QPs.size(); q++) {
    for (uint br = 0; br < nBranches; br++) {
    }
  }

  // ****************** BEGIN Printing block *************************
  oss2 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  oss2 << "[APL_THERMAL_W_iso]START" << std::endl;
  oss2 << aurostd::PaddedPRE("# Index", 6, " ") << "                           " << aurostd::PaddedPRE("Scattering rates (1/ps)", 14, " ") << std::endl;
  // ****************** END Printing block *************************

  for (uint iq = 0; iq < iQPs.size(); iq++) {
    uint q = iQPs[iq].list[0];
    for (uint br = 0; br < nBranches; br++) {
      if (freq[q][br] != 0.0) {
        for (uint q2 = 0; q2 < QPs.size(); q2++) {
          for (uint br2 = 0; br2 < nBranches; br2++) {
            double weight = exp(-(freq[q][br] - freq[q2][br2]) * (freq[q][br] - freq[q2][br2]) / (sigma[q2][br2] * sigma[q2][br2])) / sqrt(NPI) / sigma[q2][br2];
            if (std::abs(freq[q][br] - freq[q2][br2]) < 2.5 * sigma[q2][br2]) {
              for (uint at = 0; at < xstr.atoms.size(); at++) {
                xvector<xcomplex<double>> vec1, vec2;
                vec1[1] = EigenVec[q][(at * 3) + 1][br + 1];
                vec1[2] = EigenVec[q][(at * 3) + 2][br + 1];
                vec1[3] = EigenVec[q][(at * 3) + 3][br + 1];
                vec2[1] = EigenVec[q2][(at * 3) + 1][br2 + 1];
                vec2[2] = EigenVec[q2][(at * 3) + 2][br2 + 1];
                vec2[3] = EigenVec[q2][(at * 3) + 3][br2 + 1];
                xcomplex<double> prod = vec1 * vec2;
                double prod2 = (prod.re * prod.re) + (prod.im * prod.im);
                double gfac = getMassVariance(xstr.atoms[at].atomic_number);
                Wiso[iq][br] = Wiso[iq][br] + weight * prod2 * gfac;
              }
            }
          }
        }
        Wiso[iq][br] = NPI * freq[q][br] * freq[q][br] * Wiso[iq][br] / (2.0 * QPs.size());
      }
    }
  }

  // ****************** BEGIN Printing block *************************
  for (uint iq = 0; iq < iQPs.size(); iq++) {  //Loop over q points
    uint q = iQPs[iq].list[0];
    oss2 << setw(6) << q << "    ";
    for (uint br = 0; br < nBranches; br++) {  // Loop over branches
      oss2 << setw(20) << setprecision(10) << std::scientific << Wiso[iq][br];
    }
    oss2 << std::endl;
  }
  oss2 << "[APL_THERMAL_W_iso]STOP" << std::endl;
  oss2 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  // ****************** END Printing block *************************
  string ofiletemppropsname2 = "./aflow.apl.thermal_wiso.out";
  if (!aurostd::stringstream2file(oss2, ofiletemppropsname2, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.thermal_wiso.out ");
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double ThermalConductivityCalculator::getMassVariance(uint n) {

    	double MassVariance[] = {
              0, 0.00011460743, 8.32328E-08, 0.0014588232, 0, 0.00135391428, 0.00007387218, 0.00001857771, 0.00003358805, 0, 0.00082783369, 0, 0.00073988271, 0, 0.00020046752, 0, 0.00016807795, 
            0.00058238731, 0.00003509919, 0.000164, 0.000297564, 0, 0.000286456, 9.54831E-07, 0.00013287, 1.67276E-32, 9.17912E-05, 0, 0.000430773, 0.00021086, 0.000595597, 0.000197588, 
          0.00058782, 0, 0.00046279, 0.000156277, 0.000248482, 0.000109697, 6.09969E-05, 0, 0.000342629, 0, 0.000598128, 0, 0.000406665, 1.90706E-32, 0.000309478, 8.57985E-05, 0.000271603, 
          1.24494E-05, 0.000334085, 6.60751E-05, 0.000283934, 0, 0.000267781, 0, 6.23705E-05, 4.65323E-08, 2.24956E-05, 0, 0.000231599, 0, 0.000334686, 4.32857E-05, 0.000127674, 
          0, 5.20771E-05, 2.96961E-32, 7.24618E-05, 0, 8.54557E-05, 8.27273E-07, 5.25384E-05, 3.66845E-09, 6.96679E-05, 2.70849E-05, 7.45234E-05, 2.53787E-05, 3.39206E-05, 2.08217E-32, 
          6.52519E-05, 1.99659E-05, 1.94378E-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.15611E-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	  };

  return MassVariance[n];
}

//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::getCumKappaMPI(int ind, int startIndex, int endIndex, double temp, double tmfp) {
  for (uint q = startIndex; q < (uint)endIndex; q++) {  //Loop over q points
    for (uint br = 0; br < nBranches; br++) {           // Loop over branches
      double mfp = getMFP(q, br);
      if (mfp < tmfp) {
        SUM[ind] = SUM[ind] + getKappaQ(q, br, temp);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::getTao_RTA(double temp, bool ISO, bool BOUND, double NANO) {
  for (uint iq = 0; iq < iQPs.size(); iq++) {
    vector<double> kk;
    for (uint br = 0; br < nBranches; br++) {
      kk.push_back(0.0);
    }
    Wq.push_back(kk);
    Wtotal.push_back(kk);
    tau_RTA.push_back(kk);
  }

  for (uint iq = 0; iq < iQPs.size(); iq++) {  //Loop over iQPs points
    uint q = iQPs[iq].list[0];
    for (uint br = 0; br < nBranches; br++) {  // Loop over branches
      double InvTao2 = 0.0;
      if (freq[q][br] > _AFLOW_APL_EPS_) {
        if (ScatPlus_q[(iq * nBranches) + br].size() > 0) {
          q = ScatPlus_q[(iq * nBranches) + br][0][1];
          br = ScatPlus_b[(iq * nBranches) + br][0][1];
          for (uint in = 0; in < ScatPlus_q[(iq * nBranches) + br].size(); in++) {
            int q2 = ScatPlus_q[(iq * nBranches) + br][in][2];
            int q3 = ScatPlus_q[(iq * nBranches) + br][in][3];
            int br2 = ScatPlus_b[(iq * nBranches) + br][in][2];
            int br3 = ScatPlus_b[(iq * nBranches) + br][in][3];
            double fBEprime = 1.0 / (exp(7.63823296 * freq[q2][br2] / temp) - 1);
            double fBEdprime = 1.0 / (exp(7.63823296 * freq[q3][br3] / temp) - 1);
            InvTao2 = InvTao2 + (fBEprime - fBEdprime) * (ScatPlus[(iq * nBranches) + br][in]);
          }
        }
        if (ScatMinus_q[(iq * nBranches) + br].size() > 0) {
          q = ScatMinus_q[(iq * nBranches) + br][0][1];
          br = ScatMinus_b[(iq * nBranches) + br][0][1];
          for (uint in = 0; in < ScatMinus_q[(iq * nBranches) + br].size(); in++) {
            int q2 = ScatMinus_q[(iq * nBranches) + br][in][2];
            int q3 = ScatMinus_q[(iq * nBranches) + br][in][3];
            int br2 = ScatMinus_b[(iq * nBranches) + br][in][2];
            int br3 = ScatMinus_b[(iq * nBranches) + br][in][3];
            double fBEprime = 1.0 / (exp(7.63823296 * freq[q2][br2] / temp) - 1);
            double fBEdprime = 1.0 / (exp(7.63823296 * freq[q3][br3] / temp) - 1);
            InvTao2 = InvTao2 + (fBEprime + fBEdprime + 1) * (ScatMinus[(iq * nBranches) + br][in] / 2.0);
          }
        }
        Wq[iq][br] = InvTao2;
        Wtotal[iq][br] = InvTao2;
        tau_RTA[iq][br] = 1 / Wtotal[iq][br];
        //tau_RTA[iq][br] = 100;
      } else {
        Wq[iq][br] = InvTao2;
        Wtotal[iq][br] = InvTao2;
        tau_RTA[iq][br] = 0.0;
        //tau_RTA[iq][br] = 100;
      }
    }
  }

  if (ISO) {
    for (uint iq = 0; iq < iQPs.size(); iq++) {  //Loop over iQPs points
      uint q = iQPs[iq].list[0];
      for (uint br = 0; br < nBranches; br++) {  // Loop over branches
        Wtotal[iq][br] = Wtotal[iq][br] + Wiso[iq][br];
        if (freq[q][br] > _AFLOW_APL_EPS_) {
          tau_RTA[iq][br] = 1.0 / Wtotal[iq][br];
        } else {
          tau_RTA[iq][br] = 0.0;
        }
      }
    }
  }

  if (BOUND) {
    for (uint iq = 0; iq < iQPs.size(); iq++) {  //Loop over iQPs points
      uint q = iQPs[iq].list[0];
      for (uint br = 0; br < nBranches; br++) {  // Loop over branches
        Wtotal[iq][br] = Wtotal[iq][br] + aurostd::modulus(groupVel[q][br]) / NANO;
        if (freq[q][br] > _AFLOW_APL_EPS_) {
          tau_RTA[iq][br] = 1.0 / Wtotal[iq][br];
        } else {
          tau_RTA[iq][br] = 0.0;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

xmatrix<double> ThermalConductivityCalculator::getKappa(double temp) {
  xmatrix<double> tcond;
  for (uint i = 1; i < 4; i++) {    // Loop over branches
    for (uint j = 1; j < 4; j++) {  // Loop over branches
      tcond[i][j] = 0.0;
    }
  }

  for (uint q = 1; q < QPs.size(); q++) {      //Loop over q points
    for (uint br = 0; br < nBranches; br++) {  // Loop over branches
      tcond = tcond + getKappaQ(q, br, temp);
    }
  }

  return tcond;
}

//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::getKappaMPI(int ind, int startIndex, int endIndex, double temp) {
  for (uint q = startIndex; q < (uint)endIndex; q++) {  //Loop over q points
    for (uint br = 0; br < nBranches; br++) {           // Loop over branches
      SUM[ind] = SUM[ind] + getKappaQ(q, br, temp);
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::getCumKappa(double temp) {

  // ****************** BEGIN Printing block *************************
  stringstream oss4;
  oss4 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  oss4 << "[APL_CUMULATIVE_THERMAL_KAPPA at" << temp << " K]START" << std::endl;
  oss4 << aurostd::PaddedPRE("# MFP (nm)", 18, " ") << "          " << aurostd::PaddedPRE("Thermal Conductivity (W/m*K)", 14, " ") << std::endl;
  oss4 << "#                        " << aurostd::PaddedPRE("k(1,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(1,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(1,3)", 10, " ") << "          ";
  oss4 << aurostd::PaddedPRE("k(2,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(2,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(2,3)", 10, " ") << "          ";
  oss4 << aurostd::PaddedPRE("k(3,1)", 10, " ") << "          " << aurostd::PaddedPRE("k(3,2)", 10, " ") << "          " << aurostd::PaddedPRE("k(3,3)", 10, " ") << "          ";
  oss4 << std::endl;
  // ****************** END Printing block *************************
  // Cumulative kappa

  for (uint i = 0; i < 6; i++) {     // Loop over branches
    for (uint j = 0; j < 10; j++) {  // Loop over branches
      double prefactor = i + (j * 0.1);
      double tmfp = pow(10, prefactor);
      xmatrix<double> cumk;

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
      // Get the number of CPUS
      int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
      _pc.get_NCPUS(ncpus);  //CO 180214
      if (ncpus < 1) ncpus = 1;
      int qpointsPerCPU = QPs.size() / ncpus;
      for (int icpu = 0; icpu < ncpus; icpu++) {
        SUM.push_back(cumk);
      }
      // Distribute the calculation
      int startIndex, endIndex;
      std::vector<std::thread*> threads;
      for (int icpu = 0; icpu < ncpus; icpu++) {
        startIndex = icpu * qpointsPerCPU;
        endIndex = startIndex + qpointsPerCPU;
        if (((uint)endIndex > QPs.size()) ||
            ((icpu == ncpus - 1) && ((uint)endIndex < QPs.size())))
          endIndex = QPs.size();
        threads.push_back(new std::thread(&ThermalConductivityCalculator::getCumKappaMPI, this, icpu, startIndex, endIndex, temp, tmfp));
      }

      // Wait to finish all threads here!
      for (uint i = 0; i < threads.size(); i++) {
        threads[i]->join();
        delete threads[i];
      }
      threads.clear();

      for (int icpu = 0; icpu < ncpus; icpu++) {
        cumk = cumk + SUM[icpu];
      }
      SUM.clear();
#else
      for (uint q = 0; q < QPs.size(); q++) {      //Loop over q points
        for (uint br = 0; br < nBranches; br++) {  // Loop over branches
          double mfp = getMFP(q, br);
          if (mfp < tmfp) {
            cumk = cumk + getKappaQ(q, br, temp);
          }
        }
      }
#endif
      // ****************** BEGIN Printing block *************************
      oss4 << setw(20) << setprecision(10) << tmfp;
      for (uint i = 1; i < 4; i++) {
        for (uint j = 1; j < 4; j++) {
          oss4 << setw(20) << setprecision(10) << std::scientific << cumk[i][j];
        }
      }
      oss4 << std::endl;
      // ****************** END Printing block *************************
      cumk.clear();
    }
  }

  // ****************** BEGIN Printing block *************************
  oss4 << "[APL_CUMULATIVE_THERMAL_KAPPA at" << temp << " K]STOP" << std::endl;
  oss4 << "[AFLOW] **************************************************************************************************************************" << std::endl;
  string ofiletemppropsname4 = "./aflow.apl.cum_thermal_kappa.out";
  if (!aurostd::stringstream2file(oss4, ofiletemppropsname4, "POST")) {
    throw APLRuntimeError("Cannot write aflow.apl.cum_thermal_kappa.out ");
  }
  // ****************** END Printing block *************************
  //tcondc.clear();
  // END NEW VERSION getKappa
}
//////////////////////////////////////////////////////////////////////////////

double ThermalConductivityCalculator::getSigma(uint q2, uint br2) {
  // See J. Carrete Paper (Computer Physics Communications)
  double kk1 = 0.0;
  for (uint i = 1; i < 4; i++) {
    kk1 = kk1 + pow(aurostd::scalar_product((groupVel[q2][br2]), (_klattice(i) * 10.0)) / (grid[i]), 2);
  }
  double sigmaW = sqrt(kk1 / 6.0);
  return sigmaW;  // Scalebroad = 0.1 
}

//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::calculateScatProc_MPI(int startIndex) {
  //NEW Storage
  vector<double> tt;
  vector<double> mm;
  vector<xvector<int>> LPq;
  vector<xvector<int>> LPb;
  xvector<int> inQ;
  xvector<int> inB;
  xcomplex<double> iONE;
  iONE.re = 0.0;
  iONE.im = 1.0;

  //Looking for scattering processes
  int N_plus = 0;
  int N_minus = 0;
  if (startIndex == 0) {
    for (uint iq = 0; iq < iQPs.size(); iq++) {
      uint q = iQPs[iq].list[0];
      int q1x = iQPs[iq].mapq[1];
      int q1y = iQPs[iq].mapq[2];
      int q1z = iQPs[iq].mapq[3];
      for (uint br = 0; br < nBranches; br++) {
        if (freq[q][br] != 0.0) {
          for (int q2x = 0; q2x < grid(1); q2x++) {
            for (int q2y = 0; q2y < grid(2); q2y++) {
              for (int q2z = 0; q2z < grid(3); q2z++) {
                int q2 = mapQ[q2x][q2y][q2z];
                int q3x = q1x + q2x;
                if (q3x >= grid(1)) { q3x = q3x - grid(1); }
                int q3y = q1y + q2y;
                if (q3y >= grid(2)) { q3y = q3y - grid(2); }
                int q3z = q1z + q2z;
                if (q3z >= grid(3)) { q3z = q3z - grid(3); }
                int q3 = mapQ[q3x][q3y][q3z];
                for (uint br2 = 0; br2 < nBranches; br2++) {
                  // Q1+Q2=Q3
                  for (uint br3 = 0; br3 < nBranches; br3++) {
                    if ((freq[q3][br3] != 0.0) && (freq[q2][br2] != 0.0)) {
                      double sigma = getSigma(q2, br2, q3, br3);
                      if (std::abs(freq[q][br] + freq[q2][br2] - freq[q3][br3]) < 2 * sigma) {
                        N_plus = N_plus + 1;
                        inQ[1] = q;
                        inQ[2] = q2;
                        inQ[3] = q3;
                        LPq.push_back(inQ);
                        inQ.clear();
                        inB[1] = br;
                        inB[2] = br2;
                        inB[3] = br3;
                        LPb.push_back(inB);
                        inB.clear();
                        tt.push_back(0.0);
                        mm.push_back(sigma);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        //Storing Plus
        ScatPlus.push_back(tt);
        tt.clear();
        ScatPlus_s.push_back(mm);
        mm.clear();
        ScatPlus_q.push_back(LPq);
        LPq.clear();
        ScatPlus_b.push_back(LPb);
        LPb.clear();
      }
    }
  } else {
    // Q1-Q2=Q3
    for (uint iq = 0; iq < iQPs.size(); iq++) {
      uint q = iQPs[iq].list[0];
      int q1x = iQPs[iq].mapq[1];
      int q1y = iQPs[iq].mapq[2];
      int q1z = iQPs[iq].mapq[3];
      for (uint br = 0; br < nBranches; br++) {
        if (freq[q][br] != 0.0) {
          for (int q2x = 0; q2x < grid(1); q2x++) {
            for (int q2y = 0; q2y < grid(2); q2y++) {
              for (int q2z = 0; q2z < grid(3); q2z++) {
                int q2 = mapQ[q2x][q2y][q2z];
                int q3x = q1x - q2x;
                if (q3x < 0) { q3x = q3x + grid(1); }
                int q3y = q1y - q2y;
                if (q3y < 0) { q3y = q3y + grid(2); }
                int q3z = q1z - q2z;
                if (q3z < 0) { q3z = q3z + grid(3); }
                int q3 = mapQ[q3x][q3y][q3z];
                for (uint br2 = 0; br2 < nBranches; br2++) {
                  for (uint br3 = 0; br3 < nBranches; br3++) {
                    if ((freq[q3][br3] != 0.0) && (freq[q2][br2] != 0.0)) {
                      double sigma = getSigma(q2, br2, q3, br3);
                      if (std::abs(freq[q][br] - freq[q2][br2] - freq[q3][br3]) < 2 * sigma) {
                        /// BEGIN SCATERING MATRIX CALCULATION ///
                        N_minus = N_minus + 1;
                        inQ[1] = q;
                        inQ[2] = q2;
                        inQ[3] = q3;
                        LPq.push_back(inQ);
                        inQ.clear();
                        inB[1] = br;
                        inB[2] = br2;
                        inB[3] = br3;
                        LPb.push_back(inB);
                        inB.clear();
                        tt.push_back(0.0);
                        mm.push_back(sigma);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        ScatMinus.push_back(tt);
        tt.clear();
        ScatMinus_s.push_back(mm);
        mm.clear();
        ScatMinus_q.push_back(LPq);
        LPq.clear();
        ScatMinus_b.push_back(LPb);
        LPb.clear();
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::calculateScat_MPI(int startIndex, int endIndex) {
  xcomplex<double> iONE(0.0, 1.0);

  for (int ind = startIndex; ind < endIndex; ind++) {
    _logger.updateProgressBar(ind, iQPs.size() * nBranches);
    int iq = ind / nBranches;
    int br = ind % nBranches;
    // Q1+Q2=Q3
    uint q = iQPs[iq].list[0];
    if (freq[q][br] != 0.0) {
      for (uint in = 0; in < ScatPlus[(iq * nBranches) + br].size(); in++) {  //
        int q2 = ScatPlus_q[(iq * nBranches) + br][in][2];
        int q3 = ScatPlus_q[(iq * nBranches) + br][in][3];
        int br2 = ScatPlus_b[(iq * nBranches) + br][in][2];
        int br3 = ScatPlus_b[(iq * nBranches) + br][in][3];
        xvector<double> kpoint1 = QPs[q2].q * 0.1;
        xvector<double> kpoint2 = QPs[q3].q * 0.1;
        double sigma = ScatPlus_s[(iq * nBranches) + br][in];
        /// BEGIN SCATERING MATRIX CALCULATION ///
        xcomplex<double> Vp;
        Vp.re = 0.0;
        Vp.im = 0.0;
        for (uint i = 0; i < _pa.tri_index.size(); i++) {
          uint ipc1 = _pa.tri_index[i][1];
          uint isc1 = _sc.pc2scMap(ipc1);  // Equivalent atom in sc1
          uint isc2 = _pa.tri_index[i][2];
          uint ipc2 = _sc.sc2pcMap(isc2);  //Equivalent atom ins pc2
          uint isc3 = _pa.tri_index[i][3];
          uint ipc3 = _sc.sc2pcMap(isc3);  //Equivalent atom ins pc3
          xcomplex<double> phase1;
          xcomplex<double> phase2;
          phase1 = exp(iONE * scalar_product(kpoint1, tri_vec[i][0]));
          phase2 = exp(-iONE * scalar_product(kpoint2, tri_vec[i][1]));
          xcomplex<double> prefactor = (1 / sqrt(_sc.getAtomMass(isc1) * _sc.getAtomMass(isc2) * _sc.getAtomMass(isc3))) * phase1 * phase2;
          xcomplex<double> Vp0;
          Vp0.re = 0;
          Vp0.im = 0;
          for (uint alfa = 0; alfa < 3; alfa++) {  // Loop over three directions
            for (uint beta = 0; beta < 3; beta++) {
              for (uint gam = 1; gam < 4; gam++) {
                Vp0 = Vp0 + (_pc.getElementTensor(ipc1, isc2, isc3, alfa, beta, gam) * EigenVec[q][(ipc1 * 3) + alfa + 1][br + 1] * EigenVec[q2][(ipc2 * 3) + beta + 1][br2 + 1] * conj(EigenVec[q3][(ipc3 * 3) + gam][br3 + 1]));
              }
            }
          }
          Vp = Vp + prefactor * Vp0;
        }
        /// END SCATERING MATRIX CALCULATION ///
        ScatPlus[(iq * nBranches) + br][in] = (hbarp * NPI / 4.0) * exp(-(freq[q][br] + freq[q2][br2] - freq[q3][br3]) * (freq[q][br] + freq[q2][br2] - freq[q3][br3]) / (sigma * sigma)) / sqrt(NPI) / sigma / (freq[q][br] * freq[q2][br2] * freq[q3][br3]) * ((Vp.re * Vp.re) + (Vp.im * Vp.im));
        ScatPlus[(iq * nBranches) + br][in] = ScatPlus[(iq * nBranches) + br][in] * 5.60626442E08 / QPs.size();  //THz
      }
    }
    // Q1-Q2=Q3
    if (freq[q][br] != 0.0) {
      for (uint in = 0; in < ScatMinus[(iq * nBranches) + br].size(); in++) {  //
        int q2 = ScatMinus_q[(iq * nBranches) + br][in][2];
        int q3 = ScatMinus_q[(iq * nBranches) + br][in][3];
        int br2 = ScatMinus_b[(iq * nBranches) + br][in][2];
        int br3 = ScatMinus_b[(iq * nBranches) + br][in][3];
        xvector<double> kpoint1 = QPs[q2].q * 0.1;
        xvector<double> kpoint2 = QPs[q3].q * 0.1;
        double sigma = ScatMinus_s[(iq * nBranches) + br][in];
        /// BEGIN SCATERING MATRIX CALCULATION ///
        xcomplex<double> Vp;
        Vp.re = 0;
        Vp.im = 0;
        for (uint i = 0; i < _pa.tri_index.size(); i++) {
          uint ipc1 = _pa.tri_index[i][1];
          uint isc1 = _sc.pc2scMap(ipc1);  // Equivalent atom in sc1
          uint isc2 = _pa.tri_index[i][2];
          uint ipc2 = _sc.sc2pcMap(isc2);  //Equivalent atom ins pc2
          uint isc3 = _pa.tri_index[i][3];
          uint ipc3 = _sc.sc2pcMap(isc3);  //Equivalent atom ins pc2
          xcomplex<double> phase1;
          xcomplex<double> phase2;
          phase1 = exp(-iONE * scalar_product(kpoint1, tri_vec[i][0]));
          phase2 = exp(-iONE * scalar_product(kpoint2, tri_vec[i][1]));
          xcomplex<double> prefactor = (1 / sqrt(_sc.getAtomMass(isc1) * _sc.getAtomMass(isc2) * _sc.getAtomMass(isc3))) * phase1 * phase2;
          xcomplex<double> Vp0;
          Vp0.re = 0;
          Vp0.im = 0;
          for (uint alfa = 0; alfa < 3; alfa++) {  // Loop over three directions
            for (uint beta = 0; beta < 3; beta++) {
              for (uint gam = 1; gam < 4; gam++) {
                Vp0 = Vp0 + (_pc.getElementTensor(ipc1, isc2, isc3, alfa, beta, gam) * EigenVec[q][(ipc1 * 3) + alfa + 1][br + 1] * conj(EigenVec[q2][(ipc2 * 3) + beta + 1][br2 + 1]) * conj(EigenVec[q3][(ipc3 * 3) + gam][br3 + 1]));
              }
            }
          }
          Vp = Vp + prefactor * Vp0;
        }
        /// END SCATERING MATRIX CALCULATION ///
        ScatMinus[(iq * nBranches) + br][in] = (hbarp * NPI / 4.0) * exp(-(freq[q][br] - freq[q2][br2] - freq[q3][br3]) * (freq[q][br] - freq[q2][br2] - freq[q3][br3]) / (sigma * sigma)) / sqrt(NPI) / sigma / (freq[q][br] * freq[q2][br2] * freq[q3][br3]) * ((Vp.re * Vp.re) + (Vp.im * Vp.im));
        ScatMinus[(iq * nBranches) + br][in] = ScatMinus[(iq * nBranches) + br][in] * 5.60626442E08 / QPs.size();  //THz
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::iteration0() {
  for (uint q = 0; q < QPs.size(); q++) {
    vector<xvector<double>> kk;
    for (uint br = 0; br < nBranches; br++) {
      xvector<double> kk2;
      for (uint i = 1; i < 4; i++) {
        kk2[i] = 0.0;
      }
      kk.push_back(kk2);
    }
    F_n.push_back(kk);
  }
  for (uint iq = 0; iq < iQPs.size(); iq++) {  //Loop over q points
    for (uint br = 0; br < nBranches; br++) {  // Loop over branches
      for (uint p = 0; p < iQPs[iq].list.size(); p++) {
        uint q = iQPs[iq].list[p];
        F_n[q][br] = groupVel[q][br] * freq[q][br] * tau_RTA[iq][br];
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////

void ThermalConductivityCalculator::iteration(double temp) {
  for (uint q = 0; q < QPs.size(); q++) {
    vector<xvector<double>> kk;
    for (uint br = 0; br < nBranches; br++) {
      xvector<double> kk2;
      for (uint i = 1; i < 4; i++) {
        kk2[i] = 0.0;
      }
      kk.push_back(kk2);
    }
    DeltaF.push_back(kk);
  }

  for (uint iq = 0; iq < iQPs.size(); iq++) {  //Loop over iq points
    uint q = iQPs[iq].list[0];
    for (uint br = 0; br < nBranches; br++) {  // Loop over branches
      DeltaF[q][br][1] = 0.0;
      DeltaF[q][br][2] = 0.0;
      DeltaF[q][br][3] = 0.0;
      if (ScatPlus_q[(iq * nBranches) + br].size() > 0) {
        for (uint in = 0; in < ScatPlus_q[(iq * nBranches) + br].size(); in++) {
          int q2 = ScatPlus_q[(iq * nBranches) + br][in][2];
          int q3 = ScatPlus_q[(iq * nBranches) + br][in][3];
          int br2 = ScatPlus_b[(iq * nBranches) + br][in][2];
          int br3 = ScatPlus_b[(iq * nBranches) + br][in][3];
          double fBEprime = 1.0 / (exp(7.63823296 * freq[q2][br2] / temp) - 1);
          double fBEdprime = 1.0 / (exp(7.63823296 * freq[q3][br3] / temp) - 1);
          double gamma = (fBEprime - fBEdprime) * (ScatPlus[(iq * nBranches) + br][in]);
          DeltaF[q][br] = DeltaF[q][br] + gamma * (F_n[q3][br3] - F_n[q2][br2]);
        }
      }
      if (ScatMinus_q[(iq * nBranches) + br].size() > 0) {
        for (uint in = 0; in < ScatMinus_q[(iq * nBranches) + br].size(); in++) {
          int q2 = ScatMinus_q[(iq * nBranches) + br][in][2];
          int q3 = ScatMinus_q[(iq * nBranches) + br][in][3];
          int br2 = ScatMinus_b[(iq * nBranches) + br][in][2];
          int br3 = ScatMinus_b[(iq * nBranches) + br][in][3];
          double fBEprime = 1.0 / (exp(7.63823296 * freq[q2][br2] / temp) - 1);
          double fBEdprime = 1.0 / (exp(7.63823296 * freq[q3][br3] / temp) - 1);
          double gamma = (fBEprime + fBEdprime + 1) * (ScatMinus[(iq * nBranches) + br][in] / 2.0);
          DeltaF[q][br] = DeltaF[q][br] + gamma * (F_n[q3][br3] + F_n[q2][br2]);
        }
      }

      // Self symmetrization
      xvector<double> temp(3, 1);
      for (uint kk = 0; kk < iQPs[iq].lsym.size(); kk++) {
        temp = temp + xstr.pgroupk_xtal[iQPs[iq].lsym[kk]].Uc * DeltaF[q][br];
      }
      DeltaF[q][br] = temp / (double)iQPs[iq].lsym.size();
      if (Wtotal[iq][br] > 1E-24) {
        F_n[q][br] = (1 / Wtotal[iq][br]) * (groupVel[q][br] * freq[q][br] + DeltaF[q][br]);
      } else {
        F_n[q][br][1] = 0.0;
        F_n[q][br][2] = 0.0;
        F_n[q][br][3] = 0.0;
      }
      //Symmetrize DeltaF
      for (uint qq = 1; qq < iQPs[iq].list.size(); qq++) {  //Loop over q points
        int eq = iQPs[iq].list[qq];
        int symop = QPs[eq].index_fg;
        DeltaF[eq][br] = xstr.pgroupk_xtal[symop].Uc * DeltaF[q][br];
        if (Wtotal[iq][br] > 1E-24) {
          F_n[eq][br] = (1 / Wtotal[iq][br]) * (groupVel[eq][br] * freq[eq][br] + DeltaF[eq][br]);
        } else {
          F_n[q][br][1] = 0.0;
          F_n[q][br][2] = 0.0;
          F_n[q][br][3] = 0.0;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////
void ThermalConductivityCalculator::plotKappa() {
  stringstream oss;

  oss << "#!/usr/bin/gnuplot" << std::endl;
  oss << "#" << std::endl;
  oss << "reset" << std::endl;
  oss << std::endl;
  oss << "set terminal postscript eps size 3.5,2.75 enhanced color \\" << std::endl;
  oss << "  font 'Helvetica,20' linewidth 2" << std::endl;
  oss << "set output 'kappa_tensor.eps'" << std::endl;
  oss << std::endl;
  oss << " #line styles " << std::endl;
  oss << " set style line 10 lt 1 lc rgb '#cc0000' lw 2 pt 51 ps 1.5# close red" << std::endl;
  oss << " set style line 11 lt 1 lc rgb '#0072bd' lw 2 pt 11 ps 1.5# blue" << std::endl;
  oss << " set style line 12 lt 1 lc rgb '#d95319' lw 2 pt 35 ps 1.5# orange" << std::endl;
  oss << " set style line 13 lt 1 lc rgb '#edb120' lw 2 pt 20 ps 1.5# yellow" << std::endl;
  oss << " set style line 14 lt 1 lc rgb '#7e2f8e' lw 2 pt 44 ps 1.5# purple" << std::endl;
  oss << " set style line 15 lt 1 lc rgb '#77ac30' lw 2 pt 18  ps 1.5# green" << std::endl;
  oss << " set style line 16 lt 1 lc rgb '#4dbeee' lw 2 pt 9 ps 1.5# light-blue" << std::endl;
  oss << " set style line 17 lt 1 lc rgb '#a2142f' lw 2 pt 60 ps 1.5# red" << std::endl;
  oss << " set style line 18 lt 1 lc rgb '#ff00ff' lw 2 pt 17 ps 1.5# magenta" << std::endl;
  oss << std::endl;
  oss << " # Axes " << std::endl;
  oss << " set style line 101 lc rgb '#808080' lt 1" << std::endl;
  oss << " set border 3 back ls 101" << std::endl;
  oss << " set tics nomirror out scale 0.75" << std::endl;
  oss << " set xtics font 'Helvetica,12' " << std::endl;
  oss << " set ytics font 'Helvetica,12' " << std::endl;
  oss << std::endl;
  oss << " #Grid" << std::endl;
  oss << " set style line 102 lc rgb'#808080' lt 0 lw 1" << std::endl;
  oss << " set grid back ls 102" << std::endl;
  oss << std::endl;
  oss << " # Lables" << std::endl;
  oss << " set ylabel '{/Symbol k} (W/m K)'" << std::endl;
  oss << " set xlabel 'T (K)'" << std::endl;
  oss << " set yrange [-1:]" << std::endl;
  oss << " set key at 780,0.2 below" << std::endl;
  oss << std::endl;
  oss << " plot \"aflow.apl.thermal_kappa.out\"       u 1:2  w lp ls 18 title \"k(1,1)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:3  w lp ls 12 title \"k(1,2)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:4  w lp ls 10 title \"k(1,3)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:5  w lp ls 14 title \"k(2,1)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:6  w lp ls 15 title \"k(2,2)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:7  w lp ls 16 title \"k(2,3)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:8  w lp ls 17 title \"k(3,1)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:9  w lp ls 11 title \"k(3,2)\",\\" << std::endl;
  oss << "      \"aflow.apl.thermal_kappa.out\"       u 1:10 w lp ls 13 title \"k(3,3)\" " << std::endl;

  string ofiletemppropsname = "./aflow.apl.kappa.plt";
  if (!aurostd::stringstream2file(oss, ofiletemppropsname, "WRITE")) {
    throw APLRuntimeError("Cannot write aflow.apl.kappa.plt ");
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

xmatrix<double> ThermalConductivityCalculator::getKappaQ(uint q, uint br, double temp) {
  xmatrix<double> tcond;
  if (freq[q][br] == 0.) {
    tcond = 0.;
    return tcond;
  } else {
    for (uint i = 1; i < 4; i++) {    // Loop over branches
      for (uint j = 1; j < 4; j++) {  // Loop over branches
        tcond[i][j] = groupVel[q][br][i] * F_n[q][br][j];
      }
    }
    double fB1 = 1.0 / (exp(7.63823296 * freq[q][br] / temp) - 1.0);
    tcond = fB1 * (fB1 + 1) * freq[q][br] * tcond;
    tcond = 1E24 * hbar_J * hbar_J * tcond / (kb_J * temp * temp * Vol * QPs.size());
    return tcond;
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////

double ThermalConductivityCalculator::getMFP(uint q, uint br) {
  double mfp;
  if (aurostd::modulus(groupVel[q][br]) != 0.0) {
    mfp = aurostd::scalar_product(groupVel[q][br], F_n[q][br]) / aurostd::modulus(groupVel[q][br]);
  } else {
    mfp = 1E33;
  }
  return mfp;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
}
