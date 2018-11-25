//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                  Marco Esters - Duke University 2017                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters, 2018. Based on work by Jose J. Plata (AFLOW AAPL,
// DOI: 10.1038/s41524-017-0046-7) and Jesus Carrete (ShengBTE, 
// DOI: 10.1016/j.cpc.2014.02.015).
//
// This class calculates the thermal conductivity of a material using the
// Boltzmann Transport Equation (BTE). To determine the thermal conductivity,
// this class does the following steps:
//
// 1. It generates a q-point mesh according to the user input.
// 2. The frequencies, group velocities, and eigenvectors are calculated along
//    that mesh.
// 3. It calculates the intrinsic scattering rates, i.e. the parts that do not
//    depend on the temperature, of the anharmonic contributions.
// 4. It solves the BTE and calculates the thermal conductivity.
//
// See aflow_apl.h for descriptions of the classes and their members, and for
// the structs _kcell and _qpoint.

#include "aflow_apl.h"

// Some parts are written within the C++0x support in GCC, especially std::thread,
// which is implemented in gcc 4.4 and higher. For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400
#define AFLOW_APL_MULTITHREADS_ENABLE
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

// String constants for file output and exception handling
static const string _AAPL_TCOND_ERR_PREFIX_ = "apl::TCONDCalculator::";

static const int max_iter = 100;  // Maximum number of iterations for the iterative BTE solution

// Define constants and conversion factors. See AUROSTD/aurostd_xscalar.h for more.
static const double au2THz = 9.648553873170e+02;  // eV/(A amu) -> nm * THz^2
static const double hbar = PLANCKSCONSTANTEV_hbar;  // hbar in eVs
static const double hbar_J = E_ELECTRON * 1e12 * hbar;  // hbar in J/THz;
static const double BEfactor = hbar*1e12/KBOLTZEV;  // hbar/kB in K/THz
static const double TCONDfactor = 1e24 * hbar_J * hbar_J/KBOLTZ;  // BTE constants in J*K/THz^2
static const aurostd::xcomplex<double> iONE(0.0, 1.0);  // imaginary number

using aurostd::xcombos;
using aurostd::xerror;
using std::vector;

/************************************ MPI ***********************************/

#ifdef AFLOW_APL_MULTITHREADS_ENABLE

namespace apl {

//setupMPI////////////////////////////////////////////////////////////////////
// Sets up an MPI calculation.
vector<vector<int> > setupMPI(string message, Logger& log,
                              int nproc, int& ncpus) {
  if (ncpus < 1) ncpus = 1;

  if (ncpus > 1) {
    message += " (" + stringify(ncpus) + " threads)";
  }
  log.initProgressBar(message);

  return getThreadDistribution(nproc, ncpus);
}

//finishMPI///////////////////////////////////////////////////////////////////
// Finishes the MPI progress bar and deletes the threads.
void finishMPI(vector<std::thread*>& threads, Logger& log) {
  for (uint t = 0; t < threads.size(); t++) {
    threads[t]->join();
    delete threads[t];
  }
  log.finishProgressBar();
}

}  // namespace apl
#endif

/************************** CONSTRUCTOR/DESTRUCTOR **************************/

namespace apl {

//Constructor/////////////////////////////////////////////////////////////////
// Note: Need to pass supercell because it is protected in PhononCalculator
TCONDCalculator::TCONDCalculator(PhononCalculator& pc, Supercell& sc, 
                                 Logger& l) : _pc(pc), _sc(sc), _logger(l) {
  free();
  nBranches = _pc.getNumberOfBranches();
  pcell = _sc.getInputStructure();
}

TCONDCalculator::~TCONDCalculator() {
  free();
}

void TCONDCalculator::clear() {
  free();
}

void TCONDCalculator::free() {
  eigenvectors.clear();
  freq.clear();
  gvel.clear();
  intr_trans_probs.clear();
  irred_qpoints.clear();
  irred_qpts_symops.clear();
  nBranches = 0;
  nIQPs = 0;
  nQPs = 0;
  processes.clear();
  qpoints.clear();
  qptgrid.clear();
  qptmap.clear();
  temperatures.clear();
  thermal_conductivity.clear();
}

/*********************************** SETUP ***********************************/

//setCalculationOptions///////////////////////////////////////////////////////
// Sets all options for the thermal conductivity calculations, which includes
// the corrections to the scattering rates and the temperature steps.
void TCONDCalculator::setCalculationOptions(string USER_BTE, bool isotope,
                                            bool cumulative, bool fourth_order,
                                            bool boundary, double grain_size,
                                            double temp_start, double temp_end,
                                            double temp_step) {
  if (USER_BTE == "RTA") {
    calc_options.rta_only = true;
  } else if (USER_BTE == "FULL") {
    calc_options.rta_only = false;
  } else {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "setCalculationOptions()";
    string message = "Illegal value for the flag BTE. Use RTA or FULL.";
    throw xerror(function, message, _INPUT_ILLEGAL_);
  }

  calc_options.calc_isotopes = isotope;
  calc_options.calc_boundary = boundary;
  calc_options.calc_cumulative = cumulative;
  calc_options.fourth_order = fourth_order;
  calc_options.grain_size = grain_size;
  calc_options.temp_start = temp_start;
  calc_options.temp_end = temp_end;
  calc_options.temp_step = temp_step;
}

}  // namespace apl

/********************************** QPOINTS *********************************/

namespace apl {

//buildQpoints////////////////////////////////////////////////////////////////
// Builds the q-point mesh, including the irreducible q-points.
void TCONDCalculator::buildQpoints(const xvector<int>& qgrid) {
  // Calculate all q-points
  qptgrid = qgrid;
  kcell = setupReciprocalCell();
  nQPs = qptgrid[1] * qptgrid[2] * qptgrid[3];
  qptmap.assign(qptgrid[1], vector<vector<int> > (qptgrid[2], vector<int>(qptgrid[3])));
  
  qpoints = getQpointsFromGrid();
  writeQpoints();

  // Calculate the point group of the reciprocal cell. To calculate the point
  // group, we need some dummy objects to parse into the function. These
  // objects will be removed when CalculatePointGroupKCrystal is redesigned
  // to work without output streams and flags.
  if (!pcell.pgroup_xtal_calculated) {
    ofstream FileDevNull("/dev/null");
    ofstream dummy_os(NULL);
    _aflags aflags;
    aflags.QUIET = true;
    SYM::CalculatePointGroupKCrystal(FileDevNull, pcell, aflags, 
                                     false, false, dummy_os);
    FileDevNull.clear();
    FileDevNull.close();
  }

  // Calculate irreducible q-points

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus, startIndex, endIndex;
  vector<vector<int> > thread_dist;
  _pc.get_NCPUS(ncpus);
  // 8 CPUs appears to be the ideal number. More CPUs will make the
  // process slower.
  if (ncpus >= 8) {
    ncpus = 8;
  }
  string message = "Determining irreducible q-points";
  thread_dist = setupMPI(message, _logger, nQPs, ncpus);

  // Do calculations
  vector<std::thread*> threads;
  vector<vector<vector<int> > > iqpts(ncpus);
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&TCONDCalculator::getIrreducibleQpoints,
                                      this, startIndex, endIndex, std::ref(iqpts[icpu])));
  }

  // Finish up
  finishMPI(threads, _logger);
  _logger.initProgressBar("Joining threads");
  irred_qpoints = meldIrredQpoints(iqpts);
  iqpts.clear();
  _logger.finishProgressBar();
#else
  getIrreducibleQpoints(0, nQPs, irred_qpoints);
#endif
  nIQPs = (int) irred_qpoints.size();

  // Determine invariant symmetry operations for irreducible q-points.
  // Only important for the iterative solution.
  if (!calc_options.rta_only) {
    xvector<double> fpos_trans(3);
    int q;
    double tol = _AFLOW_APL_EPS_;
    irred_qpts_symops.resize(nIQPs);
    for (int iq = 0; iq < nIQPs; iq++) {
      q = irred_qpoints[iq][0];
      vector<int> sym(1, 0);  // Identity is always invariant
      for (uint symop = 1; symop < pcell.pgroupk_xtal.size(); symop++) {
        fpos_trans = pcell.pgroupk_xtal[symop].Uf * qpoints[q].fpos;
        if (SYM::AtomFPOSMatch(fpos_trans, qpoints[q].fpos, 
                               kcell.c2f, kcell.f2c,  kcell.skewed, tol)) {
          sym.push_back(symop);
        }
      }
      irred_qpts_symops[iq] = sym;
    }
  }
  writeIrredQpoints();
}

//setupReciprocalCell/////////////////////////////////////////////////////////
// Determines the reciprocal cell.
_kcell TCONDCalculator::setupReciprocalCell() {
  _kcell rec_cell;
  rec_cell.lattice = ReciprocalLattice(pcell.lattice);
  rec_cell.f2c = trasp(rec_cell.lattice);
  rec_cell.c2f = inverse(rec_cell.f2c);

  // Determine skewedness
  double tol = _AFLOW_APL_EPS_;
  xvector<double> min_distances(3);
  for (int i = 1; i < 4; i++) {
    double n = (double) qptgrid[i];
    min_distances[i] = aurostd::modulus(rec_cell.lattice(i))/n;
  }
  double min_dist = min(min_distances);
  rec_cell.skewed = SYM::isLatticeSkewed(rec_cell.lattice, min_dist, tol);

  return rec_cell;
}

//getQpointsFromGrid//////////////////////////////////////////////////////////
// Builds the q-points from the mesh dimensions.
vector<_qpoint> TCONDCalculator::getQpointsFromGrid() {
  // Convert to cartesian coordinates
  xvector<double> qpt_fpos(3);
  _qpoint qpt;
  vector<_qpoint> qps(nQPs, qpt);
  int q = 0;
  for (int qx = 0; qx < qptgrid[1]; qx++) {
    qpt_fpos[1] = (double) qx/(double) qptgrid[1];
    for (int qy = 0; qy < qptgrid[2]; qy++) {
      qpt_fpos[2] = (double) qy/(double) qptgrid[2];
      for (int qz = 0; qz < qptgrid[3]; qz++) {
        qpt_fpos[3] = (double) qz/(double) qptgrid[3];
        qpt.fpos = qpt_fpos;
        qpt.cpos = kcell.f2c * qpt.fpos;
        qpt.symop = -1;
        xvector<int> indices(3);
        indices[1] = qx; indices[2] = qy; indices[3] = qz;
        qpt.indices = indices;
        qps[q] = qpt;
        qptmap[qx][qy][qz] = q;
        q++;
      }
    }
  }
  return qps;
}

//getIrreducibleQpoints///////////////////////////////////////////////////////
// Determines the irreducible q-points of the mesh.
void TCONDCalculator::getIrreducibleQpoints(int startIndex, int endIndex,
                                            vector<vector<int> >& iqpts) {
  xvector<double> fpos_trans(3);
  double tol = _AFLOW_APL_EPS_;
  for (int q = startIndex; q < endIndex; q++) {
    bool append = true;
    if (q > startIndex) {
      for (uint symop = 0; symop < pcell.pgroupk_xtal.size(); symop++) {
        for (uint iq = 0; iq < iqpts.size(); iq++) {
          fpos_trans = pcell.pgroupk_xtal[symop].Uf * qpoints[iqpts[iq][0]].fpos;
          if (SYM::AtomFPOSMatch(fpos_trans, qpoints[q].fpos,
                                 kcell.c2f, kcell.f2c, kcell.skewed, tol)) {
            append = false;
            qpoints[q].symop = symop;
            iqpts[iq].push_back(q);
            iq = iqpts.size();
            symop = pcell.pgroup_xtal.size();
          }
        }
      }
    }
    if (append) {
      vector<int> iqpt_init(1, q);
      iqpts.push_back(iqpt_init);
    }
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    _logger.updateProgressBar(1.0/(nQPs - 1));
#endif
  }
}

//meldIrredQpoints///////////////////////////////////////////////////////////
// Stitches the results of the MPI threads for the irreducible q-points
// together.
vector<vector<int> >
    TCONDCalculator::meldIrredQpoints(const vector<vector<vector<int> > >& iqpts) {
  vector<vector<int> > irred = iqpts[0];
  xvector<double> fpos_trans(3);
  double tol = _AFLOW_APL_EPS_;
  double nprocs = 0.0;  // for progress bar
  for (uint i = 1; i < iqpts.size(); i++) {
    nprocs += (double) iqpts[i].size();
  }
  nprocs -= 1.0;

  for (uint i = 1; i < iqpts.size(); i++) {
    for (uint j = 0; j < iqpts[i].size(); j++) {
      bool append = true;
      for (uint symop = 0; symop < pcell.pgroupk_xtal.size(); symop++) {
        for (uint iq = 0; iq < irred.size(); iq++) {
          fpos_trans = pcell.pgroupk_xtal[symop].Uf * qpoints[irred[iq][0]].fpos;
          if (SYM::AtomFPOSMatch(fpos_trans, qpoints[iqpts[i][j][0]].fpos,
                                 kcell.c2f, kcell.f2c, kcell.skewed, tol)) {
            append = false;
            for (uint q = 0; q < iqpts[i][j].size(); q++) {
              irred[iq].push_back(iqpts[i][j][q]);
              if (q == 0) {
                qpoints[iqpts[i][j][q]].symop = symop;
              } else {
                // Instead of doing all the symmetry operations again, just
                // multiply the matrices and figure out which factor group
                // it belongs to (faster).
                bool mapped = false;
                xmatrix<double> Uf;
                int oldsym = qpoints[iqpts[i][j][q]].symop;
                Uf = pcell.pgroupk_xtal[oldsym].Uf * pcell.pgroupk_xtal[symop].Uf;
                for (uint isym = 0; isym < pcell.pgroupk_xtal.size(); isym++) {
                  if (Uf == pcell.pgroupk_xtal[isym].Uf) {
                    qpoints[iqpts[i][j][q]].symop = isym;
                    isym = pcell.pgroupk_xtal.size();
                    mapped = true;
                  }
                }
                if (!mapped) {
                  string function = _AAPL_TCOND_ERR_PREFIX_ + "meldIrredQpoints()";
                  stringstream message;
                  message << "Could not find the symmetry operation to map q-point ";
                  message << qpoints[iqpts[i][j][q]].fpos;
                  throw xerror(function, message, _RUNTIME_ERROR_);
                }
              }
            }
            iq = irred.size();
            symop = pcell.pgroup_xtal.size();
          }
        }
      }
      if (append) {
        irred.push_back(iqpts[i][j]);
      }
      _logger.updateProgressBar(1.0/nprocs);
    }
  }
  return irred;
}

}  // namespace apl

/*********************** FREQUENCIES/GROUP VELOCITIES ***********************/

namespace apl {

//calculateFrequenciesGroupVelocities/////////////////////////////////////////
// Calculates the frequencies and group velocities for each q-point.
// This function is mostly overhead, the actual calculation happens in
// calculateFreqGvel.
void TCONDCalculator::calculateFrequenciesGroupVelocities() {
  // MPI variables
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus, startIndex, endIndex;
  _pc.get_NCPUS(ncpus);
  string messageMPI;
  vector<vector<int> > thread_dist;
  vector<std::thread*> threads;
#endif

  _logger << "Calculating frequencies and group velocities." << apl::endl;
  // Prepare storage
  xmatrix<xcomplex<double> > eigen(nBranches, nBranches, 1, 1);
  eigenvectors.assign(nQPs, eigen);
  freq.assign(nQPs, vector<double>(nBranches));
  xvector<double> g(3);
  gvel.assign(nQPs, vector<xvector<double> >(nBranches, g));

  // Calculate frequencies and group velocities
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  messageMPI = "Frequencies and group velocities";
  thread_dist = setupMPI(messageMPI, _logger, nQPs, ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&TCONDCalculator::calculateFreqGvel,
                                      this, startIndex, endIndex));
  }
  finishMPI(threads, _logger);
#else
  calculateFreqGvel(0, nQPs);
#endif
  writeFrequencies();
  writeGroupVelocities();
}

//calculateFreqGvel///////////////////////////////////////////////////////////
// Calculates the frequencies and group velocities using the eigenvalue
// solver implemented in apl::PhononCalculator.
void TCONDCalculator::calculateFreqGvel(int startIndex, int endIndex) {
  xmatrix<xcomplex<double> > eigen(nBranches, nBranches, 1, 1);
  vector<xmatrix<xcomplex<double> > > dDynMat(3, eigen);
  for (int q = startIndex; q < endIndex; q++) {
    // Frequency
    xvector<double> f = _pc.getFrequency(qpoints[q].cpos, apl::THZ | apl::OMEGA,
                                         eigenvectors[q], dDynMat);
    freq[q] = aurostd::xvector2vector(f);  // Convert to vector to have same indexing as gvel
    // Group velocity
    for (int br = 0; br < nBranches; br++) {
      if (freq[q][br] > _AFLOW_APL_EPS_) {
        xvector<xcomplex<double> > eigenvec = eigenvectors[q].getcol(br+1);
        xvector<xcomplex<double> > eigenvec_conj = conj(eigenvec);
        for (int i = 1; i < 4; i++) {
          xcomplex<double> integral = eigenvec_conj * (dDynMat[i-1] * eigenvec);
          gvel[q][br][i] = au2THz * integral.re/(2.0 * freq[q][br]);
        }
      } else {
        for (int i = 1; i < 4; i++) {
          gvel[q][br][i] = 0.0;
        }
      }
    }
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    _logger.updateProgressBar(1.0/(nQPs - 1));
#endif
  }
}

}  // namespace apl

/************************* TRANSITION PROBABILITIES *************************/

namespace apl {

//calculateTransitionProbabilities////////////////////////////////////////////
// Overhead for the intrinsic transition probabilities. It first determines
// which scattering processes fulfill momentum and energy conservation,
// followed by the calculation of the scattering matrix. In the end, it takes
// the scattering matrix to determine the intrinsic transition probabilities.
void TCONDCalculator::calculateTransitionProbabilities(int order) {
  if (order < 3) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "calculateScatteringRates()";
    string message = "Phonon process order needs to be three or higher.";
    throw xerror(function, message, _VALUE_RANGE_);
  } else if (order > 4) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "calculateScatteringRates()";
    string message = "Phonon process order higher than four not implemented yet.";
    throw xerror(function, message, _VALUE_RANGE_);
  }

  // Prepare processes and intr_trans_probs
  if (order == 3) {
    if (calc_options.fourth_order) {
      processes.resize(2);
      intr_trans_probs.resize(2);
    } else {
      processes.resize(1);
      intr_trans_probs.resize(1);
    }
  }

  // MPI variables
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus, startIndex, endIndex;
  _pc.get_NCPUS(ncpus);
  string messageMPI;
  vector<vector<int> > thread_dist;
  vector<std::thread*> threads;
#endif

  _logger << "Calculating scattering processes." << apl::endl;
  // Prepare variables for scattering processes
  vector<vector<vector<double> > > sgms(nIQPs, vector<vector<double> >(nBranches));
  vector<vector<vector<vector<int> > > > procs(nIQPs, vector<vector<vector<int> > >(nBranches));

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  messageMPI = "Scattering processes";
  thread_dist = setupMPI(messageMPI, _logger, nIQPs, ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&TCONDCalculator::calculateScattering, this,
                                      startIndex, endIndex, order,
                                      std::ref(sgms), std::ref(procs)));
  }
  finishMPI(threads, _logger);
#else
  calculateScattering(0, nIQPs, order, sgms, procs);
#endif

  // Flatten sgms and procs to distribute all terms evenly across threads later
  vector<double> sigmas;
  int o = order - 3;
  for (int q = 0; q < nIQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      for (uint p = 0; p < procs[q][br].size(); p++) {
        processes[o].push_back(procs[q][br][p]);
        sigmas.push_back(sgms[q][br][p]);
      }
    }
  }
  sgms.clear();
  procs.clear();

  // Calculate the real space vectors for the phases
  vector<xvector<int> > clusters;  // use xvector to index IFC tensor faster
  vector<vector<xvector<double> > > phase_vectors;
  getPhaseVectors(clusters, order, phase_vectors);

  // Calculate scatting matrices
  _logger << "Calculating scattering matrix" << apl::endl;
  xcomplex<double> xcmplx(0.0, 0.0);
  vector<xcomplex<double> > scatt_mat(processes[o].size(), xcmplx);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  messageMPI = "Scattering matrix";
  thread_dist = setupMPI(messageMPI, _logger, processes[o].size(), ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&TCONDCalculator::calculateScatteringMatrix,
                                      this, startIndex, endIndex, order,
                                      std::ref(scatt_mat), clusters, phase_vectors));
  }
  finishMPI(threads, _logger);
#else
  calculateScatteringMatrix(0, processes[o].size(), order,
                            scatt_mat, clusters, phase_vectors);
#endif

  clusters.clear();
  phase_vectors.clear();

  // Calculate scattering rates without occupation numbers
  _logger << "Calculating intrinsic transition probabilities." << apl::endl;
  intr_trans_probs[o].resize(scatt_mat.size());

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  messageMPI = "Transition Probabilities";
  thread_dist = setupMPI(messageMPI, _logger, scatt_mat.size(), ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&TCONDCalculator::calculateIntrTransProbs, this,
                                      startIndex, endIndex, order, sigmas, scatt_mat));
  }
  finishMPI(threads, _logger);
#else
  calculateIntrTransProbs(0, scatt_mat.size(), order, sigmas, scatt_mat);
#endif

  sigmas.clear();
  scatt_mat.clear();
}

//calculateScattering/////////////////////////////////////////////////////////
// Calculates all the scattering processes that conserve momentum and energy.
// The delta function is replaced with an adaptive broadening scheme.
void TCONDCalculator::calculateScattering(int startIndex, int endIndex, int order,
                                          vector<vector<vector<double> > >& sigmas,
                                          vector<vector<vector<vector<int> > > >& procs) {
  // Initialize
  int q, qi, bri, nparams;
  double freqsum, sigma;
  nparams = order - 2 + 2 * order;
  xvector<int> qsum(3), qfin(3);
  xcombos signs_combo(2, order - 2, 'C', true);
  xcombos branches_combo(nBranches, order - 1, 'E', true);
  xcombos qpts_combo(nQPs, order - 2, 'E', true);
  vector<int> signs(order - 2), branches(order - 1), qpts(order - 2), process(nparams);

  // Loop over irreducible q-points for the first q-point
  for (int iq = startIndex; iq < endIndex; iq++) {
    q = irred_qpoints[iq][0];
    for (int br = 0; br < nBranches; br++) {
      if (freq[q][br] > _AFLOW_APL_EPS_) {
        process[order - 2] = iq;
        process[order - 1] = br;
        signs_combo.reset();
        // Loop over all q-vector sign combinations
        while (signs_combo.increment()) {
          signs = signs_combo.getCombo();
          for (int i = 0; i < order - 2; i++) {
            process[i] = signs[i];
          }
          qpts_combo.reset();
          // Loop over all q-points
          while (qpts_combo.increment()) {
            qpts = qpts_combo.getCombo();
            qsum = qpoints[q].indices;
            for (int i = 0; i < order - 2; i++) {
              qi = qpts[i];
              process[order + 2 * i] = qi;
              if (signs[i] == 0) {
                qsum += qpoints[qi].indices;
              } else {
                qsum -= qpoints[qi].indices;
              }
            }
            // Last q-point must fulfill momentum conservation
            for (int i = 1; i < 4; i++) {
              if (qsum[i] < 0) {
                qfin[i] = qsum[i];
                while (qfin[i] < 0) {
                  qfin[i] += qptgrid[i];
                }
              } else {
                qfin[i] = qsum[i] % qptgrid[i];
              }
            }
            process[order + 2 * (order - 2)] = qptmap[qfin[1]][qfin[2]][qfin[3]];
            branches_combo.reset();
            // Loop over all branches
            while (branches_combo.increment()) {
              bool zerofreq = false;
              branches = branches_combo.getCombo();
              freqsum = freq[q][br];
              for (int i = 0; i < order - 1; i++) {
                bri = branches[i];
                qi = process[order + 2 * i];
                if (freq[qi][bri] > _AFLOW_APL_EPS_) {
                  process[order + 2 * i + 1] = bri;
                  if ((i < order - 2) && (signs[i] == 0)) {
                    freqsum += freq[qi][bri];
                  } else {
                    freqsum -= freq[qi][bri];
                  }
                } else {
                  zerofreq = true;
                  i = order;
                }
              }
              // Finally, test if frequencies fulfill scattering condition
              if (!zerofreq) {
                sigma = getSigma(process, order);
                if (std::abs(freqsum) < 2 * sigma) {
                  sigmas[iq][br].push_back(sigma);
                  procs[iq][br].push_back(process);
                }
              }
            }
          }
        }
      }
    }
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
      _logger.updateProgressBar(1.0/(nIQPs - 1));
#endif
  }
}

//getSigma////////////////////////////////////////////////////////////////////
// Determines the sigma value (i.e. the square root of the variance) of the
// Gaussian function for the adaptive broadening scheme.
double TCONDCalculator::getSigma(const vector<int>& process, int order) {
  double sigma = 0.0;
  int q, br;
  xvector<double> dgv(3);
  for (int i = 0; i < order - 1; i++) {
    q = process[order + 2 * i];
    br = process[order + 2 * i + 1];
    if ((i == order - 2) || (process[i] == 1)) {
      dgv -= gvel[q][br];
    } else {
      dgv += gvel[q][br];
    }
  }
  for (uint i = 1; i < 4; i++) {
    double dotprod = aurostd::scalar_product(dgv, (kcell.lattice(i) * 10.0));
    sigma += pow(dotprod/qptgrid[i], 2);
  }
  // 0.1 is the scaling factor used in ShengBTE, which works fine
  return 0.1 * (order - 2) * sqrt(sigma/6.0);
}

//getPhaseVectors/////////////////////////////////////////////////////////////
// Calculate all the real space vectors for the phase factors, and generate
// a list of clusters, including all permutations. Calculating them now saves
// time because the permutations do not need to be recalculated later with
// each iteration.
void TCONDCalculator::getPhaseVectors(vector<xvector<int> >& clusters, int order,
                                      vector<vector<xvector<double> > >& phase_vectors) {
  int o = order - 3;
  for (uint ineq = 0; ineq < _pc._clusters[o].ineq_clusters.size(); ineq++) {
    for (uint c = 0; c < _pc._clusters[o].ineq_clusters[ineq].size(); c++) {
      vector<int> atoms = _pc._clusters[o].ineq_clusters[ineq][c].atoms;
      vector<xvector<double> > vectors;
      for (uint i = 1; i < atoms.size(); i++) {
        int at_eq = _sc.sc2pcMap(atoms[i]);
        int at_eq_sc = _sc.pc2scMap(at_eq);
        xvector<double> min_vec;
        min_vec = SYM::minimumCartesianVector(_pc._clusters[o].scell.atoms[atoms[i]].cpos,
                                              _pc._clusters[o].scell.atoms[atoms[0]].cpos,
                                              _pc._clusters[o].scell.lattice);
        min_vec += _pc._clusters[o].scell.atoms[atoms[0]].cpos;
        min_vec -= _pc._clusters[o].scell.atoms[at_eq_sc].cpos;
        vectors.push_back(min_vec);
      }
      vector<int> perm_atoms(atoms.size() - 1);
      for (uint i = 0; i < perm_atoms.size(); i++) {
        perm_atoms[i] = atoms[i+1];
      }
      xcombos permut(perm_atoms);
      while (permut.increment()) {
        vector<int> permutation = permut.getCombo();
        vector<int> clust(atoms.size());
        vector<xvector<double> > vec(vectors.size());
        clust[0] = _sc.sc2pcMap(atoms[0]);  // transfer to pcell
        for (uint i = 0; i < permutation.size(); i++) {
          clust[i+1] = permutation[i];
          vec[i] = vectors[permut.m_indices[i]];
        }
        clusters.push_back(aurostd::vector2xvector(clust));
        phase_vectors.push_back(vec);
      }
    }
  }
}

//calculateScatteringMatrix///////////////////////////////////////////////////
// Calculates the scattering matrix for each valid scattering process.
void TCONDCalculator::calculateScatteringMatrix(int startIndex, int endIndex, int order,
                                                vector<xcomplex<double> >& scatt_mat,
                                                const vector<xvector<int> >& clusters,
                                                const vector<vector<xvector<double> > >& phase_vectors) {
  int o = order - 3;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  uint nprocs = processes[o].size() - 1;  // For the progress bar
#endif

  uint nclusters = clusters.size();

  // Cartesian indices for the force constant tensor
  vector<xvector<int> > cart_indices;
  xcombos cart(3, clusters[0].rows, 'E', true);
  while (cart.increment()) {
    cart_indices.push_back(aurostd::vector2xvector(cart.getCombo()));
  }
  uint ncart = cart_indices.size();
  int iq, q, br, e;

  for (int i = startIndex; i < endIndex; i++) {
    xcomplex<double> matrix(0.0, 0.0), prefactor, eigen;
    xvector<double> kpoint(3);
    vector<int> atoms;
    for (uint c = 0; c < nclusters; c++) {
      // Calculate prefactor
      double phase = 0.0;
      double mass = _sc.getAtomMass(_sc.pc2scMap(clusters[c][1]));
      for (int j = 0; j < order - 1; j++) {
        mass *= _sc.getAtomMass(clusters[c][j+2]);
        q = processes[o][i][2 * j + order];
        kpoint = qpoints[q].cpos;
        if ((j == order - 2) || (processes[o][i][j] == 1)){
          phase -= scalar_product(kpoint, phase_vectors[c][j]);
        } else {
          phase += scalar_product(kpoint, phase_vectors[c][j]);
        }
      }
      prefactor = exp(iONE * phase)/sqrt(mass);

      for (uint crt = 0; crt < ncart; crt++) {
        // Calculate product of eigenvectors
        iq = processes[o][i][order - 2];
        q = irred_qpoints[iq][0];
        e = clusters[c][1] * 3 + cart_indices[crt][1] + 1;
        br = processes[o][i][order - 1] + 1;
        eigen = eigenvectors[q][e][br];
        for (int j = 0; j < order - 1; j++) {
          q = processes[o][i][2 * j + order];
          e = _sc.sc2pcMap(clusters[c][j+2]) * 3 + cart_indices[crt][j+2] + 1;
          br = processes[o][i][2 * j + order + 1] + 1;
          if ((j == order - 2) || (processes[o][i][j] == 1)) {
            eigen *= conj(eigenvectors[q][e][br]);
          } else {
            eigen *= eigenvectors[q][e][br];
          }
        }

        // Calculate scattering matrix contribution
        double ifc = (double) _pc._anharmonicIFCs[o].force_constants(clusters[c])(cart_indices[crt]);
        matrix += prefactor * ifc * eigen;
      }
    }
    scatt_mat[i] = matrix;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    _logger.updateProgressBar(1.0/nprocs);
#endif
  }
}

//calculateIntrTransProbs/////////////////////////////////////////////////////
// Calculates the intrinsic transition probabilities for each valid process.
void TCONDCalculator::calculateIntrTransProbs(int startIndex, int endIndex, int order,
                                              const vector<double>& sigmas,
                                              const vector<xcomplex<double> >& scatt_mat) {
  int o = order - 3;
  // Conversion factor so rate unit is THz
  double prefactor = getProbabilityPrefactor(order) * pow(au2THz * 10.0, 2);
  for (int i = startIndex; i < endIndex; i++) {
    double gauss = gaussian(order, sigmas[i], processes[o][i]);
    int q, br;
    double freqprod = 1.0;
    for (int j = 0; j < order; j++) {
      q = processes[o][i][2 * j + order - 2];
      if (j == 0) {
        q = irred_qpoints[q][0];
      }
      br = processes[o][i][2 * j + order - 1];
      freqprod *= freq[q][br];
    }
    intr_trans_probs[o][i] = prefactor * gauss * magsqr(scatt_mat[i])/freqprod;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    _logger.updateProgressBar(1.0/(scatt_mat.size() - 1));
#endif
  }
}

//getProbabilityPrefactor/////////////////////////////////////////////////////
// Calculates the prefactor of the intrinsic transition probability. The
// conversion factors are chosen so that the scattering rates are in units of
// THz (or 1/ps).
double TCONDCalculator::getProbabilityPrefactor(int order) {
  double prefactor = PI/4.0;
  for (int i = 2; i < order; i++) {
    prefactor *= hbar * au2THz * 1e13/nQPs;  //convert hbar into amu A^2 THz
  }
  if (order == 4) prefactor /= 2.0;
  return prefactor;
}

//gaussian////////////////////////////////////////////////////////////////////
// Calculates the value of the Gaussian function for a given scattering
// process.
double TCONDCalculator::gaussian(int order, double sigma,
                                 const vector<int>& process) {
  int iq, q, br;
  iq = process[order - 2];
  q = irred_qpoints[iq][0];
  br = process[order - 1];
  double f = freq[q][br];
  for (int i = 0; i < order - 1; i++) {
    q = process[2 * i + order];
    br = process[2 * i + order + 1];
    if ((i == order - 2) || (process[i] == 1)) {
      f -= freq[q][br];
    } else {
      f += freq[q][br];
    }
  }
  return exp(-(f * f)/(sigma * sigma))/(sqrt(PI) * sigma);
}

}  // namespace apl

/************************************ BTE ***********************************/

namespace apl {

//calculateThermalConductivity////////////////////////////////////////////////
// Overhead to calculate the thermal conductivity tensor and to output all
// results into files. It also calculates the isotope scattering rates and the
// boundary scattering rates.
void TCONDCalculator::calculateThermalConductivity() {
  // Initialize
  double tinit, tfin, tstep;
  xmatrix<double> tcond(3, 3);
  vector<vector<double> > occ, rates_iso, rates_boundary, rates_3rd, rates_4th, rates_total;
  vector<double> iso_trans_probs;
  vector<vector<int> > iso_processes;
  stringstream stream_3rd, stream_4th, stream_total;
  tinit = calc_options.temp_start;
  tfin = calc_options.temp_end;
  tstep = calc_options.temp_step;

  // Calculate the isotope and boundary scattering rates since they do not
  // explicitly depend on temperature.
  if (calc_options.calc_isotopes) {
    _logger << "Calculating intrinsic isotope scattering rates." << apl::endl;
    rates_iso = getIsotopeRates(iso_processes, iso_trans_probs);
  }
  if (calc_options.calc_boundary) {
    _logger << "Calculating intrinsic grain boundary scattering rates." << apl::endl;
    rates_boundary = getBoundaryRates();
  }

  for (double temp = tinit; temp <= tfin; temp += tstep) {
    // Calculate the scattering rates
    temperatures.push_back(temp);
    _logger << "Starting iterations for " << temp << " K." << apl::endl;
    _logger << "Calculating relaxation times" << apl::endl;
    occ = getOccupationNumbers(temp);
    rates_total.assign(nIQPs, vector<double>(nBranches, 0.0));
    rates_3rd = getAnharmonicRates(3, occ);
    stream_3rd << tempDepRatesString("ANHARMONIC", temp, rates_3rd);
    if (calc_options.fourth_order) {
      rates_4th = getAnharmonicRates(4, occ);
      stream_4th << tempDepRatesString("ANHARMONIC_4TH", temp, rates_4th);
    }
    for (int iq = 0; iq < nIQPs; iq++) {
      for (int br = 0; br < nBranches; br++) {
        rates_total[iq][br] = rates_3rd[iq][br];
        if (calc_options.calc_isotopes) {
          rates_total[iq][br] += rates_iso[iq][br];
        }
        if (calc_options.calc_boundary) {
          rates_total[iq][br] += rates_boundary[iq][br];
        }
        if (calc_options.fourth_order) {
          rates_total[iq][br] += rates_4th[iq][br];
        }
      }
    }
    stream_total << tempDepRatesString("TOTAL", temp, rates_total);

    // Calculate the thermal conductivity tensor
    _logger << "Calculating RTA" << apl::endl;
    vector<vector<xvector<double> > > mfd = getMeanFreeDispRTA(rates_total);
    tcond = calcTCOND(temp, occ, mfd);  // RTA solution

    // Iterations for the full solution
    if (!calc_options.rta_only) {
      xmatrix<double> tcond_prev(3, 3), diff(3, 3);
      int num_iter = 1;
      double max_diff;
      _logger << "Begin SCF for the Boltzmann transport equation." << apl::endl;
      std::cout << setiosflags(std::ios::fixed | std::ios::right);
      std::cout << setw(15) << "Iteration";
      std::cout << setiosflags(std::ios::fixed | std::ios::right);
      std::cout << setw(25) << "\%-difference" << std::endl;
      do {
        tcond_prev = tcond;
        getMeanFreeDispFull(rates_total, occ, iso_processes, iso_trans_probs, mfd);
        tcond = calcTCOND(temp, occ, mfd);

        // Convergence criterion based on old AAPL. It is not ideal, but
        // better than the alternatives. Using the relative difference of all
        // tensor values would cause would overemphasize very small values and
        // using an absolute convergence criterion would add a lot more
        // unnecessary steps to materials with high thermal conductivity.
        max_diff = abs((aurostd::max(tcond_prev) - aurostd::max(tcond))/aurostd::max(tcond));
        std::cout << setiosflags(std::ios::fixed | std::ios::right);
        std::cout << setw(15) << num_iter;
        std::cout << setiosflags(std::ios::fixed | std::ios::right);
        std::cout << setw(25) << std::dec << (100*max_diff) << std::endl;
        num_iter++;
      } while ((max_diff > 5e-6) && (num_iter <= max_iter));
      if (num_iter > max_iter) {
        string function = _AAPL_TCOND_ERR_PREFIX_ + "calculateThermalConductivity";
        stringstream message;
        message << "Thermal conductivity did not converge within " << max_iter << " iterations ";
        message << "at " << temp << " K.";
        throw xerror(function, message, _RUNTIME_ERROR_);
      }
      _logger << "End SCF for the Boltzmann transport equation." << apl::endl;
    }

    thermal_conductivity.push_back(tcond);
  }

  // Write output files
  if (calc_options.calc_isotopes) {
    writeTempIndepRatesFile("isotope", DEFAULT_AAPL_ISOTOPE_FILE, rates_iso);
  }
  if (calc_options.calc_boundary) {
    writeTempIndepRatesFile("boundary", DEFAULT_AAPL_BOUNDARY_FILE, rates_boundary);
  }
  writeTempDepRatesFile(DEFAULT_AAPL_RATES_3RD_FILE, stream_3rd.str());
  if (calc_options.fourth_order) {
    writeTempDepRatesFile(DEFAULT_AAPL_RATES_4TH_FILE, stream_4th.str());
  }
  writeTempDepRatesFile(DEFAULT_AAPL_RATES_FILE, stream_total.str());
  writeThermalConductivity();
  writeThermalConductivityPlot();
}

//getIsotopeRates/////////////////////////////////////////////////////////////
// Calculates the scattering rates for the isotope correction.
vector<vector<double> > TCONDCalculator::getIsotopeRates(vector<vector<int> >& iso_processes,
                                                         vector<double>& iso_trans_probs) {
  // Calculate sigmas for Gaussian first to save time
  vector<vector<double> > sigmas = getSigmaIsotope();

  vector<int> proc(4);
  vector<vector<double> > rate(nIQPs, vector<double>(nBranches, 0.0));
  xvector<xcomplex<double> > eigen1(3), eigen2(3);
  double prefactor, pearson;
  for (int iq = 0; iq < nIQPs; iq++) {
    proc[0] = iq;
    int q = irred_qpoints[iq][0];
    for (int br = 0; br < nBranches; br++) {
      if (freq[q][br] > _AFLOW_APL_EPS_) {
        prefactor = freq[q][br] * freq[q][br] * PI/(2.0 * nQPs);
        proc[1] = br;
        for (int q2 = 0; q2 < nQPs; q2++) {
          proc[2] = q2;
          for (int br2 = 0; br2 < nBranches; br2++) {
            proc[3] = br2;
            if (std::abs(freq[q][br] - freq[q2][br2]) < 2 * sigmas[q2][br2]) {
              if (calc_options.calc_isotopes) {
                iso_processes.push_back(proc);
              }
              double gauss = gaussian(2, sigmas[q2][br2], proc);
              double rate_proc = 0;
              for (uint at = 0; at < pcell.atoms.size(); at++) {
//                int atomic_number = pcell.atoms[at].atomic_number;
//                if (pearson[atomic_number] != 0) {
                pearson = GetPearsonCoefficient(pcell.atoms[at].atomic_number);
//                pearson = pcell.atoms[at].pearson_coefficient;
                if (pearson != 0.0) {
                  for (int i = 1; i < 4; i++) {
                    eigen1(i) = conj(eigenvectors[q][(3 * at) + i][br + 1]);
                    eigen2(i) = eigenvectors[q2][(3 * at) + i][br2 + 1];
                  }
                  rate_proc += gauss * pearson * magsqr(eigen1 * eigen2);
                }
              }
              rate_proc *= prefactor;
              if (calc_options.calc_isotopes) {
                iso_trans_probs.push_back(rate_proc);
              }
              rate[iq][br] += rate_proc;
            }
          }
        }
      }
    }
  }
  return rate;
}

//getSigmaIsotope/////////////////////////////////////////////////////////////
// Calculates the sigma value for the isotope correction. The scaling factor
// is missing here because it cuts out too many scattering processes.
vector<vector<double> > TCONDCalculator::getSigmaIsotope() {
  vector<vector<double> > sigmas(nQPs, vector<double>(nBranches, 0.0));
  for (int q = 0; q < nQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      for (int i = 1; i < 4; i++) {
        double dotprod = aurostd::scalar_product(gvel[q][br], kcell.lattice(i) * 10.0);
        sigmas[q][br] += pow(dotprod/qptgrid[i], 2);
      }
      sigmas[q][br] = sqrt(sigmas[q][br]/6.0);
    }
  }

  return sigmas;
}

//getBoundaryRates////////////////////////////////////////////////////////////
// Calculates the rates of the grain boundary scattering processes.
vector<vector<double> > TCONDCalculator::getBoundaryRates() {
  vector<vector<double> > rate(nIQPs, vector<double>(nBranches, 0.0));
  for (int iq = 0; iq < nIQPs; iq++) {
    for (int br = 0; br < nBranches; br++) {
      int q = irred_qpoints[iq][0];
      rate[iq][br] = aurostd::modulus(gvel[q][br])/calc_options.grain_size;
    }
  }
  return rate;
}

//getOccupationNumbers////////////////////////////////////////////////////////
// Calculates the phonon occupation numbers for all q-points using the
// Bose-Einstein distribution
vector<vector<double> > TCONDCalculator::getOccupationNumbers(double temp) {
  vector<vector<double> > occ(nQPs, vector<double>(nBranches, 0.0));
  for (int q = 0; q < nQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      occ[q][br] = 1.0/(exp(BEfactor * freq[q][br]/temp) - 1.0);
    }
  }
  return occ;
}

//getAnharmonicRates//////////////////////////////////////////////////////////
// Calculates the scattering rates of the anharmonic phonon-phonon scattering
// processes.
vector<vector<double> > TCONDCalculator::getAnharmonicRates(int order,
                                                            const vector<vector<double> >& occ) {
  int o = order - 3;
  vector<vector<double> > rates(nIQPs, vector<double>(nBranches, 0.0));
  uint nproc = processes[o].size();
  int iq, br;
  double occ_term;
  for (uint i = 0; i < nproc; i++) {
    iq = processes[o][i][order - 2];
    br = processes[o][i][order - 1];
    occ_term = getOccupationTerm(order, processes[o][i], occ);
    rates[iq][br] += occ_term * intr_trans_probs[o][i];
  }
  return rates;
}

//getOccupationTerm///////////////////////////////////////////////////////////
// For the anharmonic scattering processes, this function calculates the term
// that contains the phonon numbers.
double TCONDCalculator::getOccupationTerm(int order, const vector<int>& process,
                                          const vector<vector<double> >& occ) {
  double term;
  if (order == 3) {
    int q2, q3, br2, br3;
    q2 = process[3]; br2 = process[4];
    q3 = process[5]; br3 = process[6];
    if (process[0] == 0) {
      term = occ[q2][br2] - occ[q3][br3];
    } else {
      term = 0.5 * (1.0 + occ[q2][br2] + occ[q3][br3]);
    }
  } else {
    int iq, q, q2, q3, q4, br, br2, br3, br4;
    iq = process[2];
    q = irred_qpoints[iq][0]; br = process[3];
    q2 = process[4]; br2 = process[5];
    q3 = process[6]; br3 = process[7];
    q4 = process[8]; br4 = process[9];
    if (process[0] + process[1] == 0) {  // ++
      term = (1 + occ[q2][br2]) * (1 + occ[q3][br3]) * occ[q4][br4]/2.0;
    } else if (process[0] + process[1] == 2) {  // --
      term = occ[q2][br2] * occ[q3][br3] * occ[q4][br4]/6.0;
    } else {  // +-
      term = (1 + occ[q2][br2]) * occ[q3][br3] * occ[q4][br4]/2.0;
    }
    term /= occ[q][br];
  }
  return term;
}

//getMeanFreeDispRTA//////////////////////////////////////////////////////////
// Calculates the mean free displacement in the relaxation time approximation.
vector<vector<xvector<double> > >
    TCONDCalculator::getMeanFreeDispRTA(const vector<vector<double> >& rates) {
  xvector<double> xvec(3);
  vector<vector<xvector<double> > > mfd(nQPs, vector<xvector<double> >(nBranches, xvec));
  for (int iq = 0; iq < nIQPs; iq++) {
    for (uint i = 0; i < irred_qpoints[iq].size(); i++) {
      int q = irred_qpoints[iq][i];
      for (int br = 0; br < nBranches; br++) {
        if (rates[iq][br] > 0.0) {
          mfd[q][br] = gvel[q][br] * freq[q][br]/rates[iq][br];
        }
      }
    }
  }
  return mfd;
}

//getMeanFreeDispFull/////////////////////////////////////////////////////////
// Calculates an iteration for the mean free displacement in the full solution
// to the BTE. It first calculates the deviations from the RTA (delta) for all
// irreducible q-points and then symmetrizes delta.
void TCONDCalculator::getMeanFreeDispFull(const vector<vector<double> >& rates_total,
                                          const vector<vector<double> >& occ,
                                          const vector<vector<int> >& iso_processes,
                                          const vector<double>& iso_trans_probs,
                                          vector<vector<xvector<double> > >& mfd) {
  xvector<double> xvec(3), correction(3);
  vector<vector<xvector<double> > > delta(nQPs, vector<xvector<double> >(nBranches, xvec));
  uint nprocs;
  int iq, q, br;

  nprocs = processes[0].size();
  for (uint i = 0; i < nprocs; i++) {
    iq = processes[0][i][1]; q = irred_qpoints[iq][0]; br = processes[0][i][2];
    correction = getMFDCorrection(3, processes[0][i], mfd, occ);
    delta[q][br] += intr_trans_probs[0][i] * correction;
  }

  if (calc_options.calc_isotopes) {
    int q2, br2;
    nprocs = iso_processes.size();
    for (uint i = 0; i < nprocs; i++) {
      q = iso_processes[i][0]; q2 = iso_processes[i][2];
      br = iso_processes[i][1]; br2 = iso_processes[i][3];
      delta [q][br] += mfd[q2][br2] * iso_trans_probs[i];
    }
  }

  if (calc_options.fourth_order) {
    nprocs = processes[1].size();
    for (uint i = 0; i < nprocs; i++) {
      iq = processes[1][i][2]; q = irred_qpoints[iq][0]; br = processes[1][i][3];
      correction = getMFDCorrection(4, processes[1][i], mfd, occ);
      delta[q][br] += intr_trans_probs[0][i] * correction;
    }
  }

  // Symmetrize delta and calculate mean free displacement
  int eq, symop;
  for (int iq = 0; iq < nIQPs; iq++) {
    // Symmetrize irreducible q-point
    q = irred_qpoints[iq][0];
    xmatrix<double> Uc(3, 3);
    for (uint i = 0; i < irred_qpts_symops[iq].size(); i++) {
      symop = irred_qpts_symops[iq][i];
      Uc += pcell.pgroupk_xtal[symop].Uc;
    }
    Uc = 1.0/irred_qpts_symops[iq].size() * Uc;
    for (int br = 0; br < nBranches; br++) {
      if (rates_total[iq][br] > 0.0) {
        delta[q][br] = Uc * delta[q][br];
        mfd[q][br] = (gvel[q][br] * freq[q][br] + delta[q][br])/rates_total[iq][br];
      } else {
        mfd[q][br] = xvec;  // Set all components to zero
      }
    }
    // Symmetrize equivalent q-points
    for (uint i = 1; i < irred_qpoints[iq].size(); i++) {
      eq = irred_qpoints[iq][i];
      symop = qpoints[irred_qpoints[iq][i]].symop;
      Uc = pcell.pgroupk_xtal[symop].Uc;
      for (int br = 0; br < nBranches; br++) {
        if (rates_total[iq][br] > 0.0) {
          mfd[eq][br] = (gvel[eq][br] * freq[eq][br] + Uc * delta[q][br])/rates_total[iq][br];
        } else {
          mfd[eq][br] = xvec;
        }
      }
    }
  }
}

//getMFDCorrection////////////////////////////////////////////////////////////
// Gets the correction to the mean free displacement.
xvector<double> TCONDCalculator::getMFDCorrection(int order, const vector<int>& process,
                                                  const vector<vector<xvector<double> > >& mfd,
                                                  const vector<vector<double> >& occ) {
  int q, br;
  q = process[3 * (order - 1) - 1];
  br = process[3 * (order - 1)];
  xvector<double> correction = mfd[q][br];
  for (int i = 0; i < order - 2; i++) {
    q = process[2 * i + order];
    br = process[2 * i + order + 1];
    if (process[i] == 0) {
      correction -= mfd[q][br];
    } else {
      correction += mfd[q][br];
    }
  }
  return getOccupationTerm(order, process, occ) * correction;
}

//calcTCOND///////////////////////////////////////////////////////////////////
// Calculates the thermal conductivity tensor.
xmatrix<double> TCONDCalculator::calcTCOND(double temp,
                                           const vector<vector<double> >& occ,
                                           const vector<vector<xvector<double> > >& mfd) {
  xmatrix<double> tcond(3, 3);
  double prefactor = 1E24 * hbar_J * hbar_J/(KBOLTZ * temp * temp * nQPs * pcell.Volume());
  for (int q = 0; q < nQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      bool include = true;
      if (freq[q][br] < _AFLOW_APL_EPS_) {
        include = false;
      } else if (calc_options.calc_cumulative) {
        // Only takes scattering processes into account that have a mean free
        // path smaller than the grain size.
        double mfpath = scalar_product(mfd[q][br], gvel[q][br])/aurostd::modulus(gvel[q][br]);
        if (mfpath > calc_options.grain_size) {
          include = false;
        }
      }
      if (include) {
        double x = occ[q][br] * (occ[q][br] + 1) * freq[q][br];
        for (int i = 1; i < 4; i++) {
          for (int j = 1; j < 4; j++) {
            double tc = gvel[q][br][i] * mfd[q][br][j];
            tcond(i, j) += x * tc;
          }
        }
      }
    }
  }
  tcond = prefactor * tcond;
  return tcond;
}

}  // namespace apl

/******************************** FILE OUTPUT *******************************/

namespace apl {

//writeQpoints////////////////////////////////////////////////////////////////
// Writes the Cartesian coordinates of each q-point into a file.
void TCONDCalculator::writeQpoints() {
  stringstream output;
  string filename = DEFAULT_AAPL_QPOINTS_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_QPOINTS]START" << std::endl;
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << setw(10) << "# Index";
  output << setw(20) << " ";
  output << "Q-points (1/Angstrom)" << std::endl;

  // Body
  for (int q = 0; q < nQPs; q++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << setw(10) << q;
    for (int i = 1; i < 4; i++) {
      output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      output << setw(20) << setprecision(10) << std::scientific << qpoints[q].cpos[i];
    }
    output << std::endl;
  }

  output << "[AAPL_QPOINTS]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  // Write to file
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeQpoints";
    string message = "Could not write q-points to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeIrredQpoints///////////////////////////////////////////////////////////
// Writes the Cartesian coordinates and the multiplicity of the irreducible
// q-points into a file.
void TCONDCalculator::writeIrredQpoints() {
  stringstream output;
  string filename = DEFAULT_AAPL_IRRQPTS_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_IRREDUCIBLE_QPOINTS]START" << std::endl;
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << setw(10) << "# Index";
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << setw(15) << "Multiplicity";
  output << setw(20) << " ";
  output << "Q-points (1/Angstrom)" << std::endl;

  // Body
  for (int iq = 0; iq < nIQPs; iq++) {
    int q = irred_qpoints[iq][0];
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << setw(10) << q;
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << setw(15) << irred_qpoints[iq].size();
    for (int i = 1; i < 4; i++) {
      output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      output << setw(20) << setprecision(10) << std::scientific << qpoints[q].cpos[i];
    }
    output << std::endl;
  }

  output << "[AAPL_IRREDUCIBLE_QPOINTS]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  // Write to file
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeIrredQpoints";
    string message = "Could not write irreducible q-points to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeFrequencies////////////////////////////////////////////////////////////
// Writes the frequencies into a file. Each row belongs to a q-point, and
// each column belongs to a phonon branch.
void TCONDCalculator::writeFrequencies() {
  stringstream output;
  string filename = DEFAULT_AAPL_FREQ_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_FREQUENCY]START" << std::endl;
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << setw(10) << "# Q-point";
  output << setw(20) << " ";
  output << "Frequencies (THz)" << std::endl;

  // Body
  for (int q = 0; q < nQPs; q++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << setw(10) << q;
    for (int br = 0; br < nBranches; br++) {
      output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      output << setw(20) << setprecision(10) << std::scientific << freq[q][br];
    }
    output << std::endl;
  }

  output << "[AAPL_FREQUENCY]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  // Write to file
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeFrequencies";
    string message = "Could not write frequencies to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeGroupVelocities////////////////////////////////////////////////////////
// Writes the group velocities into a file. Each row belongs to a q-point,
// and each column triplet belongs to a phonon branch.
void TCONDCalculator::writeGroupVelocities() {
  stringstream output;
  string filename = DEFAULT_AAPL_GVEL_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_GROUP_VELOCITY]START" << std::endl;
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << setw(10) << "# Q-point";
  output << setw(20) << " ";
  output << "Group Velocity (km/s)" << std::endl;

  // Body
  for (int q = 0; q < nQPs; q++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << setw(10) << q;
    for (int br = 0; br < nBranches; br++) {
      for (int i = 1; i < 4; i++) {
        output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        output << setw(20) << setprecision(10) << std::scientific << gvel[q][br][i];
      }
      output << setw(5) << " ";
    }
    output << std::endl;
  }

  output << "[AAPL_GROUP_VELOCITY]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  // Write to file
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeGroupVelocities";
    string message = "Could not write group velocities to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeTempIndepRatesFile/////////////////////////////////////////////////////
// Writes temperature-independent scattering rates into a file. These are
// e.g. the isopte or grain boundary scattering rates. Each row belongs to
// an irreducible q-point, and each column belongs to a phonon branch.
void TCONDCalculator::writeTempIndepRatesFile(string tag, string filename,
                                             const vector<vector<double> >& rates) {
  stringstream output;
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_SCATTERING_RATE_" << aurostd::toupper(tag) << "]START" << std::endl;
  output << setw(10) << "# Q-point";
  output << setw(20) << " ";
  output << "Scattering rates (1/ps)" << std::endl;
  for (int iq = 0; iq < nIQPs; iq++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << setw(10) << irred_qpoints[iq][0];
    for (int br = 0; br < nBranches; br++) {
      output << setw(20) << setprecision(10) << std::scientific << rates[iq][br];
    }
    output << std::endl;
  }
  output << "[AAPL_SCATTERING_RATE_" << aurostd::toupper(tag) << "]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeTempDepRateFile";
    stringstream message;
    message << "Could not write file " << filename << ".";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//tempDepRateString///////////////////////////////////////////////////////////
// Outputs the printing block of temperature-dependent scattering rates for
// one temperature step. The function writeTempDepRatesFile then outputs
// all these blocks into one output file.
string TCONDCalculator::tempDepRatesString(string tag, double temp,
                                          const vector<vector<double> >& rates) {
  stringstream ratestream;
  ratestream << AFLOWIN_SEPARATION_LINE << std::endl;
  ratestream << "[AAPL_SCATTERING_RATE_" << aurostd::toupper(tag);
  ratestream << "= " << setprecision(2) << std::fixed << temp << " K]START" << std::endl;
  ratestream << setw(10) << "# Q-point";
  ratestream << setw(20) << " ";
  ratestream << "Scattering rates (1/ps)" << std::endl;
  for (int iq = 0; iq < nIQPs; iq++) {
    ratestream << std::setiosflags(std::ios::fixed | std::ios::right);
    ratestream << setw(10) << irred_qpoints[iq][0];
    for (int br = 0; br < nBranches; br++) {
      ratestream << setw(20) << setprecision(10) << std::scientific << rates[iq][br];
    }
    ratestream << std::endl;
  }
  ratestream << "[AAPL_SCATTERING_RATE_" << aurostd::toupper(tag);
  ratestream << "= " << setprecision(2) << std::dec << temp << " K]STOP" << std::endl;
  ratestream << AFLOWIN_SEPARATION_LINE << std::endl;
  return ratestream.str();
}

//writeTempDepRatesFile///////////////////////////////////////////////////////
// Writes temperature-dependent scattering rates, e.g. anharmonic scattering
// rates, into a file. Each rows belongs to an irreducible q-point and each
// column belongs to a phonon branch.
void TCONDCalculator::writeTempDepRatesFile(string filename,
                                            const string& content) {
  aurostd::string2file(content, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeTempDepRateFile";
    stringstream message;
    message << "Could not write file " << filename << ".";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeThermalConductivity////////////////////////////////////////////////////
// Writes the thermal conductivity tensor as a function of temperature into
// a file.
void TCONDCalculator::writeThermalConductivity() {
  stringstream output;
  string filename = DEFAULT_AAPL_TCOND_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_THERMAL_CONDUCTIVITY]START" << std::endl;
  output << setw(8) << "# T (K)";
  output << setw(75) << " ";
  output << "Thermal Conductivity (W/m*K)" << std::endl;
  output << setw(8) << " ";
  string xyz [] = {"x", "y", "z"};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      output << setw(7) << " ";
      string k = "k(" + xyz[i] + "," + xyz[j] + ")";
      output << k;
      output << setw(7) << " ";
    }
  }
  output << std::endl;

  // Body
  for (uint t = 0; t < temperatures.size(); t++) {
    output << setw(8) << std::fixed << setprecision(2) << temperatures[t];
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        output << setw(2) << " ";
        output << setprecision(10) << std::scientific << thermal_conductivity[t][i][j];
        output << setw(2) << " ";
      }
    }
    output << std::endl;
  }
  output << "[AAPL_THERMAL_CONDUCTIVITY]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeThermalConductivity";
    string message = "Could not write thermal conductivities to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeThermalConductivityPlot////////////////////////////////////////////////
// Write a gnuplot script for the thermal conductivity tensor as a function
// of temperature into a file.
void TCONDCalculator::writeThermalConductivityPlot() {
  stringstream output;
  string filename = DEFAULT_AAPL_TCOND_PLOT_FILE;
  string tcondfile = DEFAULT_AAPL_TCOND_FILE;
  string font = DEFAULT_GNUPLOT_EPS_FONT;
  string font_it = DEFAULT_GNUPLOT_EPS_FONT_ITALICS;
  string font_greek = DEFAULT_GNUPLOT_GREEK_FONT_ITALICS;

  output << "#!/usr/bin/gnuplot" << std::endl << std::endl;
  output << "reset" << std::endl << std::endl;
  output << "set terminal postscript eps size 3.5,2.75 enhanced color\\" << std::endl;
  output << "  font '" << font << ",20' linewidth 2" << std::endl;
  output << "set output 'thermal_conductivity.eps'" << std::endl << std::endl;
  // The color palette was designed to be accessible for people with color
  // vision deficiencies. When changing the colors, please make sure that
  // the colors are still distinguishable for everyone!
  output << "# line styles" << std::endl;
  output << " set style line 10 lt 1 lc rgb '#000000' lw 2 pt 17 ps 1.5 # black" << std::endl;
  output << " set style line 11 lt 1 lc rgb '#004949' lw 2 pt 35 ps 1.5 # dark green" << std::endl;
  output << " set style line 12 lt 1 lc rgb '#009292' lw 2 pt 51 ps 1.5 # light green" << std::endl;
  output << " set style line 13 lt 1 lc rgb '#490092' lw 2 pt 44 ps 1.5 # purple" << std::endl;
  output << " set style line 14 lt 1 lc rgb '#b66dff' lw 2 pt 18 ps 1.5 # lavender" << std::endl;
  output << " set style line 15 lt 1 lc rgb '#6db6ff' lw 2 pt  9 ps 1.5 # light blue" << std::endl;
  output << " set style line 16 lt 1 lc rgb '#924900' lw 2 pt 60 ps 1.5 # brown" << std::endl;
  output << " set style line 17 lt 1 lc rgb '#d55e00' lw 2 pt 11 ps 1.5 # orange" << std::endl;
  output << " set style line 18 lt 1 lc rgb '#edb120' lw 2 pt 20 ps 1.5 # yellow" << std::endl;
  output << std::endl;
  output << "# Axes" << std::endl;
  output << " set style line 101 lc rgb '#000000' lt 1" << std::endl;
  output << " set border 3 back ls 101" << std::endl;
  output << " set tics nomirror out" << std::endl;
  output << " set xtics font '" << font << ",12'" << std::endl;
  output << " set ytics font '" << font << ",12'" << std::endl;
  output << std::endl;
  output << "# Grid" << std::endl;
  output << " set style line 102 lc rgb '#808080' lt 0 lw 1" << std::endl;
  output << " set grid back ls 102" << std::endl;
  output << std::endl;
  output << "# Labels" << std::endl;
  output << " set xlabel '{/" << font_it << " T} (K)'" << std::endl;
  output << " set ylabel '{/" << font_greek << " k} (W/m K)'" << std::endl;
  output << " set yrange [-1:]" << std::endl;
  output << " set key below spacing 1.1" << std::endl;
  // Spacing between key and xlabel, necessary for gnuplot version 5 and above
  output << " set bmargin 7" << std::endl;
  output << std::endl;
  output << "# Data" << std::endl;
  output << " plot \"" << tcondfile << "\" u 1:2  w lp ls 10";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " xx}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:3  w lp ls 11";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " xy}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:4  w lp ls 12";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " xz}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:5  w lp ls 13";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " yx}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:6  w lp ls 14";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " yy}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:7  w lp ls 15";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " yz}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:8  w lp ls 16";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " zx}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:9  w lp ls 17";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " zy}\",\\" << std::endl;
  output << "      \"" << tcondfile << "\" u 1:10 w lp ls 18";
  output << " title \"{/" << font_greek << " k}^{/" << font_it << " zz}\"" << std::endl;
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeThermalConductivityPlot";
    string message = "Could not write thermal conductivity plot to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow Marco Esters - Duke University 2018                *
// *                                                                         *
// ***************************************************************************
