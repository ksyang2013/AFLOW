//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                  Marco Esters - Duke University 2017                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters (ME), 2018. Based on work by Jose J. Plata (AFLOW
// AAPL, DOI: 10.1038/s41524-017-0046-7). 
//
// Routines to set up AAPL calculations.

#include "aflow_apl.h"

#define _AFLOW_AAPL_DIRECTORY_PREFIX_ string("ARUN.AAPL_C")
#define _DEBUG_AAPL_FORCES_ false
using std::vector;
using std::string;
using aurostd::xerror;

static const string _AAPL_FORCES_ERR_PREFIX_ = "apl::PhononCalculator::";

namespace apl {

//setAnharmonicOptions////////////////////////////////////////////////////////
// Sets the calculation options for the calculations of the anharmonic IFCs.
void PhononCalculator::setAnharmonicOptions(int iter, double mix, double threshold) {
  _anharmonicIFCOptions options;
  options.max_iter = iter;
  options.mixing_coefficient = mix;
  options.sumrule_threshold = threshold;
  anharmonic_IFC_options = options;
}

//buildVaspAAPL///////////////////////////////////////////////////////////////
// Creates the folders for the VASP calculations.
bool PhononCalculator::buildVaspAAPL(const ClusterSet& clst) {
  bool stagebreak = false;
  _logger << "Managing directories for ";
  if (clst.order == 3) {
    _logger << "3rd";
  } else {
    _logger << "4th";
  }
  _logger << " order IFCs." << apl::endl;

  // Check if supercell is built
  if (!_supercell.isConstructed()) {
    string function = _AAPL_FORCES_ERR_PREFIX_ + "buildVaspAAPL";
    string message = "The supercell has not been initialized yet.";
    throw xerror(function, message, _RUNTIME_INIT_);
  }

  // Determine the number of runs so the run ID in the folder name can be
  // padded with the appropriate number of zeros.
  int nruns = 0;
  for (uint ineq = 0; ineq < clst.ineq_distortions.size(); ineq++) {
    nruns += clst.ineq_distortions[ineq].distortions.size();
  }
  vector<_xinput> xinp(nruns);

  int idxRun = 0;
  for (uint ineq = 0; ineq < clst.ineq_distortions.size(); ineq++) {
    vector<int> atoms = clst.ineq_distortions[ineq].atoms;  // Declare to make code more legible
    _ineq_distortions idist = clst.ineq_distortions[ineq];  // Declare to make code more legible
    for (uint dist = 0; dist < idist.distortions.size(); dist++) {
      xinp[idxRun] = _xInput;
      vector<int> distortions = idist.distortions[dist][0];  // Declare to make code more legible

      // Set up working directory and generate distorted structure
      string runname = buildFolderNameAAPL(distortions, atoms, clst.order, idxRun, nruns);
      xinp[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname);
      applyDistortionsAAPL(xinp[idxRun], clst.distortion_vectors, distortions, atoms);

      // Check if the directory has already been run by looking for input and output files
      // If so, continue, but keep the distorted structure in the list for post-processing
      if (filesExistPhonons(xinp[idxRun])) {
        idxRun++;
        continue;
      }
      // If not, create folders and aflow.in
      _logger << "Creating " << xinp[idxRun].getDirectory() << apl::endl;
      createAflowInPhonons(xinp[idxRun], runname);
      idxRun++;
      stagebreak = true;
    }
  }
  xInputsAAPL.push_back(xinp);
  return stagebreak;
}

//buildFolderNameAAPL/////////////////////////////////////////////////////////
// Creates the name of the folder for the VASP calculation.
string PhononCalculator::buildFolderNameAAPL(const vector<int>& distortions,
                                             const vector<int>& atoms, const int& ord,
                                             const int& run, const int& nruns) {
  std::stringstream folder;
  //Prefix
  folder << "./" << _AFLOW_AAPL_DIRECTORY_PREFIX_ << ord << "_";

  // Run ID with padding
  std::stringstream runss, nrunss;
  runss << run;
  nrunss << nruns;
  for (uint i = runss.str().length(); i < nrunss.str().length(); i++) {
    folder << "0";
  }
  folder << run;

  // Atom and distortion IDs
  for (uint at = 0; at < atoms.size(); at++) {
    folder << "_A" << (at+1) << "_" << atoms[at];
    folder << "_D" << (at+1) << "_" << distortions[at];
  }
  return folder.str();
}

//applyDistortionsAAPL////////////////////////////////////////////////////////
// Applies the inequivalent distortions to the supercell structures.
void PhononCalculator::applyDistortionsAAPL(_xinput& xinp,
                                            const vector<xvector<double> >& distortion_vectors,
                                            const vector<int>& distortions,
                                            const vector<int>& atoms) {
  xinp.setXStr(_supercell.getSupercellStructureLight());  // Light copy because we don't need symmetry, etc.
  for (uint at = 0; at < atoms.size(); at++) {
    int atsc, dist_index;
    atsc = atoms[at];
    dist_index = distortions[at];
    xvector<double> dist_cart = distortion_vectors[dist_index];
    while (((at + 1) < atoms.size()) && (atoms[at] == atoms[at+1])) {
       at++;
       dist_index = distortions[at];
       dist_cart += distortion_vectors[dist_index];
    }
    // Normalize dist_cart coordinates to 1 so that distortions have the same magnitude
    for (int i = 1; i < 4; i++) {
      if (abs(dist_cart(i)) < _ZERO_TOL_) {
        dist_cart(i) = 0.0;
      } else {
        dist_cart(i) /= std::abs(dist_cart(i));
      }
    }
    dist_cart *= DISTORTION_MAGNITUDE;
    xinp.getXStr().atoms[atsc].cpos += dist_cart;
    xinp.getXStr().atoms[atsc].fpos = C2F(xinp.getXStr().lattice, xinp.getXStr().atoms[atsc].cpos);
  }
}

//calculateAnharmonicIFCs/////////////////////////////////////////////////////
// Calculates the anharmonic IFCs.
void PhononCalculator::calculateAnharmonicIFCs(ClusterSet& clst) {
  _logger << "Checking file integrity for anharmonic IFCs." << apl::endl;
  int o = clst.order - 3;
  if (!outfileFoundAnywherePhonons(xInputsAAPL[o])) {
    throw APLStageBreak();
  }
  outfileFoundEverywherePhonons(xInputsAAPL[o]);
  if (_calculateZeroStateForces) {
    subtractZeroStateForces(xInputsAAPL[o]);
  }
  AnharmonicIFCs ifcs(xInputsAAPL[o], clst, DISTORTION_MAGNITUDE,
                      anharmonic_IFC_options, _logger);
  _anharmonicIFCs.push_back(ifcs);
  // Clear runs from memory because they are no longer needed
  xInputsAAPL[o].clear();
}


void PhononCalculator::readAnharmonicIFCs(const string& filename, ClusterSet& clst) {
  _logger << "Reading anharmonic IFCs from file " << filename << "." << apl::endl;
  int o = clst.order - 3;
  AnharmonicIFCs ifcs(filename, clst, DISTORTION_MAGNITUDE,
                      anharmonic_IFC_options, _logger);
  _anharmonicIFCs.push_back(ifcs);
  // Clear runs from memory because they are no longer needed
  xInputsAAPL[o].clear();
}

} // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow Marco Esters - Duke University 2018                *
// *                                                                         *
// ***************************************************************************
