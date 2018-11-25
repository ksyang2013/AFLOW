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
// This class calculates the anharmonic interatomic force constants (IFCs). It
// reads the forces of the inequivalent distortions first, and then calculates
// the forces of the equivalent distortions. These forces are used to
// calculate the IFCs of the inequivalent clusters that are then symmetrized
// in an iterative procedure using sum rules.
//
// See aflow_apl.h for descriptions of the classes and their members, and for
// the struct _linearCombinations.

#include "aflow_apl.h"

#define _DEBUG_AAPL_IFCS_ false

using std::vector;
using aurostd::xcombos;
using aurostd::xerror;

static const string _AAPL_IFCS_ERR_PREFIX_ = "apl::AnharmonicIFCs::";

// tform represents a tensor transformation containing the index and the
// coefficients. vector<vector<int> > holds the indices, vector<double>
// the coefficients.
typedef vector<std::pair<vector<vector<int> >, vector<double> > > tform;
// v5int defined for brevity
typedef vector<vector<vector<vector<vector<int> > > > > v5int;

/************************ CONSTRUCTORS/DESTRUCTOR ***************************/

namespace apl {

//Constructors////////////////////////////////////////////////////////////////
// Default constructor
AnharmonicIFCs::AnharmonicIFCs(const vector<_xinput>& xInp,
                               ClusterSet& _clst,
                               const double& dist_mag, 
                               const _anharmonicIFCOptions& options,
                               Logger& l) : clst(_clst), _logger(l) {
  free();
  distortion_magnitude = dist_mag;
  max_iter = options.max_iter;
  order = clst.order;
  mixing_coefficient = options.mixing_coefficient;
  sumrule_threshold = options.sumrule_threshold;
  cart_indices = getCartesianIndices();

  _logger << "Reading forces for anharmonic IFCs from VASP calculations." << apl::endl;
  vector<aurostd::xtensor<double> > force_tensors = storeForces(xInp);

  _logger << "Caulating anharmonic IFCs." << apl::endl;
  vector<aurostd::xtensor<double> > ifcs_unsym(force_tensors.size());
  for (uint f = 0; f < force_tensors.size(); f++) {
    aurostd::xtensor<double> ifcs = calculateUnsymmetrizedIFCs(clst.ineq_distortions[f],
                                                               force_tensors[f]);
    ifcs_unsym[f] = ifcs;
  }
  force_tensors.clear();

  _logger << "Symmetrizing IFCs." << apl::endl;
  force_constants = symmetrizeIFCs(ifcs_unsym);
}

// From file
AnharmonicIFCs::AnharmonicIFCs(const string& filename,
                               ClusterSet& _clst,
                               const double& dist_mag,
                               const _anharmonicIFCOptions& options,
                               Logger& l) : clst(_clst), _logger(l) {
  free();
  distortion_magnitude = dist_mag;
  max_iter = options.max_iter;
  order = clst.order;
  mixing_coefficient = options.mixing_coefficient;
  sumrule_threshold = options.sumrule_threshold;
  readIFCsFromFile(filename);
}

//Copy constructor
const AnharmonicIFCs& AnharmonicIFCs::operator=(const AnharmonicIFCs& that) {
  if (this != &that) {
    _logger = that._logger;
    cart_indices = that.cart_indices;
    clst = that.clst;
    distortion_magnitude = that.distortion_magnitude;
    force_constants = that.force_constants;
    max_iter = that.max_iter;
    mixing_coefficient = that.mixing_coefficient;
    order = that.order;
    sumrule_threshold = that.sumrule_threshold;
  }
  return *this;
}

//Destructor//////////////////////////////////////////////////////////////////
AnharmonicIFCs::~AnharmonicIFCs() {
  free();
}

//free////////////////////////////////////////////////////////////////////////
// Clears all vectors and resets all values.
void AnharmonicIFCs::free() {
  cart_indices.clear();
  distortion_magnitude = 0.0;
  force_constants.clear();
  max_iter = 0;
  mixing_coefficient = 0.0;
  order = 0;
  sumrule_threshold = 0.0;
}

}  // namespace apl

/*************************** INITIAL CALCULATIONS ***************************/

namespace apl {

//getCartesianIndices/////////////////////////////////////////////////////////
// Returns a list of Cartesian indices. Since they are looped over frequently,
// it is quicker to calculate them once at the beginning.
vector<vector<int> > AnharmonicIFCs::getCartesianIndices() {
  vector<vector<int> > indices;
  xcombos crt_ind(3, order, 'E', true);
  while (crt_ind.increment()) {
    indices.push_back(crt_ind.getCombo());
  }
  return indices;
}

// BEGIN Forces
//storeForces/////////////////////////////////////////////////////////////////
// Stores the forces from the VASP calculations. Each item in the vector holds
// the force tensor for a set of distorted atoms.
vector<aurostd::xtensor<double> > AnharmonicIFCs::storeForces(const vector<_xinput>& xInp) {
  vector<aurostd::xtensor<double> > force_tensors;
  int idxRun = 0;
  for (uint id = 0; id < clst.ineq_distortions.size(); id++) {
    aurostd::xtensor<double> forces = getForces(id, idxRun, xInp);
    force_tensors.push_back(forces);
  }
  return force_tensors;
}

//getForces///////////////////////////////////////////////////////////////////
// Retrieves all forces from the calculations. Also transforms the forces
// into the forces of the equivalent distortions.
aurostd::xtensor<double> AnharmonicIFCs::getForces(int id, int& idxRun,
                                                   vector<_xinput> xInp) {
  _ineq_distortions ineq_dists = clst.ineq_distortions[id];
  int natoms = (int) clst.scell.atoms.size();
  uint ndim = ineq_dists.atoms.size() + 2;  // No. distortions x no. atoms x 3 Cart. dimensions
  vector<int> tensor_shape(ndim, 5);  // Tensor dimensions
  tensor_shape[ndim-2] = natoms - 1;
  tensor_shape[ndim-1] = 2;
  vector<int> zeros(ndim, 0);  // Makes tensor use zero-based indices
  aurostd::xtensor<double> force_tensor(tensor_shape, zeros);

  int attrans, fg;
  vector<int> indices;
  for (uint ineq = 0; ineq < ineq_dists.distortions.size(); ineq++) {
    // For the inequivalent distortion, just read the forces from VASP 
    vector<xvector<double> > qmforces = xInp[idxRun].getXStr().qm_forces;
    indices = ineq_dists.distortions[ineq][0];
    for (int at = 0; at < natoms; at++) { 
      force_tensor(indices)[at] = qmforces[at];
    }
    for (uint i = 1; i < ineq_dists.distortions[ineq].size(); i++) {
      // For each equivalent distortion, transform the forces using symmetry
      vector<xvector<double> > qmforces_trans(natoms, xvector<double>(3));
      indices = ineq_dists.distortions[ineq][i];
      fg = ineq_dists.rotations[ineq][i];
      vector<int> transformation_map = ineq_dists.transformation_maps[ineq][i];
      for (int at = 0; at < natoms; at++) {
        attrans = getTransformedAtom(transformation_map, at);
        force_tensor(indices)[at] = clst.pcell.fgroup[fg].Uc * qmforces[attrans];
      }
    }
    idxRun++;
  }
  return force_tensor;
}

//getTransformedAtom//////////////////////////////////////////////////////////
// Retrieves the atom that is obtained by the transformation in the given
// symmetry map.
int AnharmonicIFCs::getTransformedAtom(const vector<int>& symmap, const int& at) {
  for (uint i = 0; i < symmap.size(); i++) {
    if (symmap[i] == at) {
      return i;
    }
  }
  // If the for-loop runs until the end, the atom was not found
  string function = _AAPL_IFCS_ERR_PREFIX_ + "getTransformedAtom";
  stringstream message;
  message << "Could not transform atom " << at;
  throw xerror(function, message, _RUNTIME_ERROR_);
}
//END Forces

// BEGIN Calculate unsymmetrized force constants
//calculateUnsymmetrizedIFCs//////////////////////////////////////////////////
// Calculates the IFCs of the inequivalent clusters from the forces.
aurostd::xtensor<double>
    AnharmonicIFCs::calculateUnsymmetrizedIFCs(const _ineq_distortions& idist,
                                               const aurostd::xtensor<double>& forces) {
  // Initialize tensor
  vector<int> tensor_shape(order + 1, 2), zeros(order + 1, 0);
  uint natoms = clst.scell.atoms.size();
  tensor_shape[0] = natoms - 1;
  aurostd::xtensor<double> ifcs(tensor_shape, zeros);

  // Calculate force constants from the forces
  vector<int> cart_coords;
  for (uint c = 0; c < idist.clusters.size(); c++) {
    int at, cl;
    cl = idist.clusters[c];
    at = clst.ineq_clusters[cl][0].atoms[order - 1];
    for (int i = 0; i < clst.nifcs; i++) {
      ifcs[at](cart_indices[i]) = calculateIFC(forces, at, cart_indices[i], idist.atoms);
    }
  }
  return ifcs;
}

//calculateIFC////////////////////////////////////////////////////////////////
// Calculates a specific IFC from the forces using the central difference
// method. The numerator and denominator are calculated separately.
//
// The numerator is a sum of the forces for the set of atoms distorted along
// the Cartesian indices into positive and negative directions. The sign of
// the force depends on the number of negative distortions in the set.
//
// The denominator is the product of the distortion lengths times two (the
// factor of two comes from the central difference method). If the same atom
// gets distorted into the same direction multiple times, the distortion
// length needs to be adjusted so that the total length is the same as the
// distortion magnitude chosen by the user.
double AnharmonicIFCs::calculateIFC(const aurostd::xtensor<double>& forces, int at,
                                    const vector<int>& cart_ind, const vector<int>& atoms) {
  double ifc, num, denom, sign;
  // Numerator
  num = 0.0;
  uint ncart = cart_ind.size();
  vector<int> dist(ncart - 1), dist_signs(ncart - 1);
  int c = cart_ind[ncart - 1];
  xcombos bitenum(2, ncart - 1, 'E', true);  // Bit enumerator indicating the signs of the distortions
  while (bitenum.increment()) {
    for (uint i = 0; i < ncart; i++) {
      dist[i] = cart_ind[i];
    }
    dist_signs = bitenum.getCombo();
    sign = 1.0;
    for (uint d = 0; d < dist.size(); d++) {
      if (dist_signs[d] != 0) {
        dist[d] += 3;
        sign *= -1.0;
      }
    }
    num -= sign*forces(dist)[at][c]; // Subtract because negative force is needed for IFCs
  }
  
  // Denominator
  int atom, count, cart;
  denom = 1.0;
  for (uint at = 0, natoms = atoms.size(); at < natoms; at++) {
    atom = atoms[at];
    cart = cart_ind[at];
    count = 1;
    denom *= distortion_magnitude;
    while ((at + count < natoms) &&
           (atoms[at+count] == atom) &&
           (cart_ind[at+count] == cart)) {
      count++;
      denom *= distortion_magnitude;
    }
    denom *= 2.0/((double) count);
    at += count - 1;
  }

  ifc = num/denom;
  return ifc;
}
// END Calculate unsymmetrized force constants

}  // namespace apl

/****************************** SYMMETRIZATION ******************************/

namespace apl {

//symmetrizeIFCs//////////////////////////////////////////////////////////////
// Symmetrizes the IFCs using an iterative procedure to ensure 
// that the force constants sum to zero.
// 1. The IFCs of the inequivalent clusters will be symmetrized according to
//    the linear combinations found while determining ClusterSets.
// 2. The force constants will be transformed to the other clusters using the
//    symmetry of the crystal.
// 3. Determine the deviations of the IFC sums from zero.
// 4. If at least one deviation is above the chosen threshold, correct the
//    linearly dependent IFCs. If none are above the threshold or too many
//    iterations have been performed, stop the iteration procedure.
//
// Check the typedefs at the beginning of the file for tform and v5int
aurostd::xtensor<double>
    AnharmonicIFCs::symmetrizeIFCs(const vector<aurostd::xtensor<double> >& ifcs_unsym) {
  // Initialize tensors
  aurostd::xtensor<double> ifcs = initializeIFCTensor(ifcs_unsym);
  aurostd::xtensor<double> dev_from_zero = initializeDeviationsFromZero();
  aurostd::xtensor<double> abssum = dev_from_zero;
  vector<vector<int> > reduced_clusters = getReducedClusters();

  // Tensor transformations
  vector<vector<tform> > transformations(clst.ineq_clusters.size());
  v5int eq_ifcs(clst.ineq_clusters.size());
  getTensorTransformations(eq_ifcs, transformations);
  vector<vector<vector<int> > > all_clusters = getAllClusters(eq_ifcs);

  // Do iterations
  int num_iter = 0;
  double max_err;
  _logger << "Begin SCF for anharmonic force constants." << apl::endl;
  std::cout << setiosflags(std::ios::fixed | std::ios::right);
  std::cout << setw(15) << "Iteration";
  std::cout << setiosflags(std::ios::fixed | std::ios::right);
  std::cout << setw(20) << "Abs. max. error" << std::endl;
  do {
    // 1. Symmetrize using linear combinations
    applyLinCombs(ifcs);  

    // 2. Transform IFCs using symmetry and permutations
    transformIFCs(transformations, ifcs);

    // 3. Determine deviations from zero
    calcSums(reduced_clusters, ifcs, dev_from_zero, abssum);
    max_err = aurostd::max(aurostd::abs(dev_from_zero));
    std::cout << setiosflags(std::ios::fixed | std::ios::right);
    std::cout << setw(15) << num_iter;
    std::cout << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    std::cout << setw(20) << max_err << std::endl;

    // 4. Correct IFCs
    if (max_err > sumrule_threshold) {
      correctIFCs(ifcs, dev_from_zero, abssum, all_clusters, eq_ifcs);
    }
    num_iter++;
  } while ((num_iter <= max_iter) && (max_err > sumrule_threshold));
  _logger << "End SCF for anharmonic force constants." << apl::endl;
  if (num_iter > max_iter) {
    string function = _AAPL_IFCS_ERR_PREFIX_ + "symmetrizeIFCs";
    stringstream message;
    message << "Anharmonic force constants did not converge within " << max_iter << " iterations.";
    throw xerror(function, message, _RUNTIME_ERROR_);
  } else {
    return ifcs;
  }
}

// BEGIN Initializers
//initializeIFCTensor/////////////////////////////////////////////////////////
// Initializes the IFC tensor by populating the values for the inequivalent
// clusters with the unsymmetriced IFCs.
aurostd::xtensor<double>
    AnharmonicIFCs::initializeIFCTensor(const vector<aurostd::xtensor<double> >& ifcs_unsym) {
  // Initialize tensor
  uint natoms_pcell, natoms_scell;
  natoms_pcell = clst.pcell.atoms.size();
  natoms_scell = clst.scell.atoms.size();
  vector<int> tensor_shape(2*order);
  tensor_shape[0] = natoms_pcell - 1;  // The first atom has to be in the primitive cell
  for (int i = 1; i < order; i++) {
    tensor_shape[i] = natoms_scell - 1;
  }
  for (int i = order; i < 2*order; i++) {
    tensor_shape[i] = 2;
  }
  vector<int> zeros(2*order, 0); // Make tensor use zero-based indices
  aurostd::xtensor<double> tensor_init(tensor_shape, zeros);

  // Populate tensor with unsymmetrized IFCs
  for (uint ineq = 0; ineq < clst.ineq_distortions.size(); ineq++) {
    vector<int> atoms = clst.ineq_distortions[ineq].atoms;
    atoms[0] = clst.sc2pcMap[atoms[0]];  // Transform to pcell
    tensor_init(atoms) = ifcs_unsym[ineq];
  }
  return tensor_init;
} 

//initializeDeviationsFromZero////////////////////////////////////////////////
// Initializes the tensor that holds all sum rule errors. Note that all sizes
// need to be subtracted by 1 to get 0-based indexing.
aurostd::xtensor<double> AnharmonicIFCs::initializeDeviationsFromZero() {
  vector<int> tensor_shape(2*order - 1), zeros(2*order - 1);
  tensor_shape[0] = clst.pcell.atoms.size() - 1;
  for (int i = 1; i < order - 1; i++) {
    tensor_shape[i] = clst.scell.atoms.size() - 1;
  }
  for (int i = order - 1; i < 2*order - 1; i++) {
    tensor_shape[i] = 2;
  }
  aurostd::xtensor<double> deviations(tensor_shape, zeros);
  return deviations;
}

//getReducedClusters//////////////////////////////////////////////////////////
// Determines, for each set of inequivalent clusters, a uinque set of clusters
// that do not contain the last atom of the clusters. This set is important
// for the sum rules as they frequently require summations over the last atom
// within a set of inequivalent clusters.
vector<vector<int> > AnharmonicIFCs::getReducedClusters() {
  vector<vector<int> > reduced_clusters;
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    for (uint c = 0; c < clst.ineq_clusters[ineq].size(); c++) {
      vector<int> cluster(order - 1);
      cluster[0] = clst.sc2pcMap[clst.ineq_clusters[ineq][c].atoms[0]];  // transfer to pcell
      for (int i = 1; i < order - 1; i++) {
        cluster[i] = clst.ineq_clusters[ineq][c].atoms[i];
      }
      bool append = true;
      for (uint r = 0; r < reduced_clusters.size(); r++) {
        append = false;
        for (int i = 0; i < order - 1; i++) {
          if (cluster[i] != reduced_clusters[r][i]) {
            append = true;
            i = order;
          }
        }
        if (!append) {  // If append stays false, the reduced cluster is not new
          r = reduced_clusters.size();
        }
      }
      if (append) {
        reduced_clusters.push_back(cluster);
      }
    }
  }
  return reduced_clusters;
}

//getTensorTransformations////////////////////////////////////////////////////
// This algorithm does two things: it generates the tensor transformations
// for each cluster to transform the IFCs of the inequivalent clusters; and
// it generates a list of equivalent IFCs for each inequivalent cluster to
// calculate the corrections.
//
// Check the typedefs at the beginning of the file for tform and v5int
void AnharmonicIFCs::getTensorTransformations(v5int& eq_ifcs,
                                              vector<vector<tform> >& transformations) {
  typedef vector<vector<vector<vector<int> > > > v4int;
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    vector<tform> transform(clst.ineq_clusters[ineq].size() - 1);
    v4int eq(clst.nifcs, vector<vector<vector<int> > >(clst.ineq_clusters[ineq].size()));
    int ind = 0;
    for (int crt = 0; crt < clst.nifcs; crt++) {
      eq[ind][0].push_back(cart_indices[crt]);
      ind++;
    }
    for (uint c = 1; c < clst.ineq_clusters[ineq].size(); c++) {
      tform trf;
      _cluster cluster_trans = clst.ineq_clusters[ineq][c];
      int fg, perm, rw, cl, p;
      vector<int> atoms_trans = cluster_trans.atoms;
      atoms_trans[0] = clst.sc2pcMap[atoms_trans[0]];  // transfer to pcell
      fg = cluster_trans.fgroup;
      perm = cluster_trans.permutation;
      for (int itrans = 0; itrans < clst.nifcs; itrans++) {
        std::pair<vector<vector<int> >, vector<double> > t;
        int ind_orig = 0;
        for (int iorig = 0; iorig < clst.nifcs; iorig++) {
          double coeff = 1.0;
          for (int o = 0; o < order; o++) {
            rw = cart_indices[itrans][o] + 1;
            p = clst.permutations[perm][o];
            cl = cart_indices[iorig][p] + 1;
            coeff *= clst.pcell.fgroup[fg].Uc[rw][cl];
            if (abs(coeff) < _ZERO_TOL_) {
              coeff = 0.0;
              o = order;
            }
          }
          if (abs(coeff) > _ZERO_TOL_) {
            t.first.push_back(cart_indices[iorig]);
            t.second.push_back(coeff);
            eq[ind_orig][c].push_back(cart_indices[itrans]);
          }
          ind_orig++;
        }
        trf.push_back(t);
      }
      transform[c-1] = trf;
    }
    eq_ifcs[ineq] = eq;
    transformations[ineq] =  transform;
  }
}

//getAllClusters//////////////////////////////////////////////////////////////
// ineq_clusters of ClusterSets do not contain permutations, but they are
// important for the correction procedure. This function expands the clusters
// to contain all permutations. It also updates eq_ifcs.
//
// See the top of this file for the typedef of v5int.
vector<vector<vector<int> > > AnharmonicIFCs::getAllClusters(v5int& eq_ifcs) {
  vector<vector<vector<int> > > all_clusters(clst.ineq_clusters.size());
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    // First append the original clusters. This is necessary to maintain a
    // one-to-one mapping with eq_ifcs.
    for (uint c = 0; c < clst.ineq_clusters[ineq].size(); c++) {
      all_clusters[ineq].push_back(clst.ineq_clusters[ineq][c].atoms);
    }
    // Now do the permutations
    for (uint c = 0; c < clst.ineq_clusters[ineq].size(); c++) {
      vector<int> atoms_orig(all_clusters[ineq][c].size() - 1);
      vector<int> atoms_permut(all_clusters[ineq][c].size());
      atoms_permut[0] = all_clusters[ineq][c][0];
      for (uint at = 0; at < atoms_orig.size(); at++) {
        atoms_orig[at] = all_clusters[ineq][c][at+1];
      }
      xcombos perm(atoms_orig);
      ++perm;  // The first permutation is the original, so it can be skipped
      while (perm.increment()) {
        vector<int> permut = perm.getCombo();
        for (uint iperm = 0; iperm < permut.size(); iperm++) {
          atoms_permut[iperm+1] = permut[iperm];
        }
        all_clusters[ineq].push_back(atoms_permut);
        // Update eq_ifcs
        for (uint eq = 0; eq < eq_ifcs[ineq].size(); eq++) {
          vector<vector<int> > eq_indices;
          for (uint i = 0; i < eq_ifcs[ineq][eq][c].size(); i++) {
            vector<int> indices_orig = eq_ifcs[ineq][eq][c][i];
            vector<int> indices_permut(permut.size() + 1);
            indices_permut[0] = indices_orig[0];
            for (uint j = 1; j < eq_ifcs[ineq][eq][c][i].size(); j++) {
              indices_permut[j] = indices_orig[perm.m_indices[j-1] + 1];
            }
            eq_indices.push_back(indices_permut);
          }
          eq_ifcs[ineq][eq].push_back(eq_indices);
        }
      }
    }
  }
  return all_clusters;
}

// END Initializers

// BEGIN Iterations
//applyLinCombs///////////////////////////////////////////////////////////////
// Sets the lineraly dependent IFCs according to the obtained linear
// combinations.
void AnharmonicIFCs::applyLinCombs(aurostd::xtensor<double>& ifcs) {
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    vector<int> atoms = clst.ineq_clusters[ineq][0].atoms;
    atoms[0] = clst.sc2pcMap[atoms[0]];  // transfer to pcell
    _linearCombinations lcomb = clst.linear_combinations[ineq];
    for (uint d = 0; d < lcomb.dependent.size(); d++) {
      vector<int> cart_ind = cart_indices[lcomb.dependent[d]];
      ifcs(atoms)(cart_ind) = 0.0;  // reset
      for (uint lc = 0; lc < lcomb.indices[d].size(); lc++) {
        vector<int> cart_ind_indep = cart_indices[lcomb.indices[d][lc]];
        ifcs(atoms)(cart_ind) += lcomb.coefficients[d][lc] * ifcs(atoms)(cart_ind_indep);
      }
    }
  }
}

//transformIFCs///////////////////////////////////////////////////////////////
// Transforms all IFCs using symmetry. The first two loops go over all
// equivalent clusters and transform the inequivalent clusters using tensor
// transformations.
//
// See the top of this file for the typedef of tform.
void AnharmonicIFCs::transformIFCs(const vector<vector<tform> >& transformations,
                                   aurostd::xtensor<double>& ifcs) {
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    vector<int> atoms_orig = clst.ineq_clusters[ineq][0].atoms;
    applyPermutations(atoms_orig, ifcs);
    atoms_orig[0] = clst.sc2pcMap[atoms_orig[0]];  // transfer to pcell
    for (uint c = 1; c < clst.ineq_clusters[ineq].size(); c++) {
      vector<int> atoms_trans = clst.ineq_clusters[ineq][c].atoms;
      atoms_trans[0] = clst.sc2pcMap[atoms_trans[0]];  // transfer to pcell
      for (int itrans = 0; itrans < clst.nifcs; itrans++) {
        ifcs(atoms_trans)(cart_indices[itrans]) = 0.0;  // reset
        std::pair<vector<vector<int> >, vector<double> > transf = transformations[ineq][c-1][itrans];
        for (uint t = 0; t < transf.first.size(); t++) {
          vector<int> cart_indices_orig = transf.first[t];
          double coeff = transf.second[t];
          ifcs(atoms_trans)(cart_indices[itrans]) += coeff * ifcs(atoms_orig)(cart_indices_orig);
        }
      }
      atoms_trans[0] = clst.pc2scMap[atoms_trans[0]];  // transfer to scell for applyPermutations
      applyPermutations(atoms_trans, ifcs);
    }
  }
}

//applyPermutations///////////////////////////////////////////////////////////
// Applies permutation symmetry to the force constants. The first index is not
// included in the permutations to keep the first atom in the primitive cell.
void AnharmonicIFCs::applyPermutations(vector<int> atoms,
                                       aurostd::xtensor<double>& ifcs) {
  atoms[0] = clst.sc2pcMap[atoms[0]];  // transfer to pcell
  vector<int> atoms_orig(order - 1);
  for (int i = 1; i < order; i++) {
    atoms_orig[i-1] = atoms[i];
  }
  xcombos permut(atoms_orig);
  ++permut;  // The first permutation is the original, so it can be skipped
  while (permut.increment()) {
    vector<int> atoms_permut = permut.getCombo();
    for (int crt = 0; crt < clst.nifcs; crt++) {
      vector<int> indices = cart_indices[crt];
      vector<int> indices_permut(indices.size());
      indices_permut[0] = indices[0];
      for (uint i = 1; i < indices_permut.size(); i++) {
        int perm_ind = permut.m_indices[i-1];
        indices_permut[i] = indices[perm_ind + 1];
      }
      ifcs[atoms[0]](atoms_permut)(indices_permut) = ifcs(atoms)(indices);
    }
  }
}

//calcSums////////////////////////////////////////////////////////////////////
// Determines the deviations from zero and the sum of the absolute values
// for each reduced cluster. Both are used for the correction of the IFCs
// while the deviations from zero are also used as errors for the sum rules.
void AnharmonicIFCs::calcSums(const vector<vector<int> >& reduced_clusters,
                              const aurostd::xtensor<double>& ifcs, 
                              aurostd::xtensor<double>& dev_from_zero,
                              aurostd::xtensor<double>& abssum) {
  uint natoms = clst.scell.atoms.size();
  for (uint r = 0; r < reduced_clusters.size(); r++) {
    vector<int> tensor_shape(order, 2), zeros(order, 0);
    aurostd::xtensor<double> dev(tensor_shape, zeros), asum(tensor_shape, zeros);
    for (uint at = 0; at < natoms; at++) {
      dev += ifcs(reduced_clusters[r])[at];
      asum += abs(ifcs(reduced_clusters[r])[at]);
    }
    dev_from_zero(reduced_clusters[r]) = dev;
    abssum(reduced_clusters[r]) = asum;
  }
}

//correctIFCs/////////////////////////////////////////////////////////////////
// Corrects the IFCs the comply with the sum rule using weighted averages.
//
// Check the typedef at the beginning of the file for v5int.
void AnharmonicIFCs::correctIFCs(aurostd::xtensor<double>& ifcs,
                                 const aurostd::xtensor<double>& dev_from_zero, 
                                 const aurostd::xtensor<double>& abssum,
                                 const vector<vector<vector<int> > > all_clusters,
                                 const v5int& eq_ifcs) {
  for (uint ineq = 0; ineq < all_clusters.size(); ineq++) {
    uint nclusters = all_clusters[ineq].size();
    vector<int> ineq_cluster = all_clusters[ineq][0];
    ineq_cluster[0] = clst.sc2pcMap[ineq_cluster[0]];  // transfer to pcell
    // Calculate correction terms
    vector<aurostd::xtensor<double> > correction_terms(nclusters);
    for (uint c = 0; c < nclusters; c++) {
      correction_terms[c] = getCorrectionTerms(all_clusters[ineq][c],
                                               ifcs, dev_from_zero, abssum);
    }

    // Correct the linearly independent IFCs
    _linearCombinations lcomb = clst.linear_combinations[ineq];
    vector<vector<int> > indices;
    uint nindep = lcomb.independent.size();
    indices.resize(nindep);
    for (uint i = 0; i < nindep; i++) {
      indices[i] = cart_indices[lcomb.independent[i]];
    }
    for (uint i = 0; i < indices.size(); i++) {
      int neq = 0;
      double corrected_ifc = 0.0;
      for (uint c = 0; c < nclusters; c++) {
        vector<int> cluster = all_clusters[ineq][c];
        cluster[0] = clst.sc2pcMap[cluster[0]];  // transfer to pcell
        vector<vector<int> > eq;
        eq = eq_ifcs[ineq][lcomb.independent[i]][c];
        uint eqsize = eq.size();
        for (uint e = 0; e < eqsize; e++) {
          if (ifcs(cluster)(eq[e]) != 0.0) {
            neq++;
            corrected_ifc += correction_terms[c](eq[e]) * ifcs(ineq_cluster)(indices[i]) / ifcs(cluster)(eq[e]);
          }
        }
        for (uint dep = 0; dep < lcomb.indep2depMap[i].size(); dep++) {
          int dependent = lcomb.indep2depMap[i][dep];
          eq = eq_ifcs[ineq][dependent][c];
          for (uint e = 0; e < eqsize; e++) {
            if (ifcs(cluster)(eq[e]) != 0.0) {
              neq++;
              corrected_ifc += correction_terms[c](eq[e]) * ifcs(ineq_cluster)(indices[i]) / ifcs(cluster)(eq[e]);
            }
          }
        }
      }
      if (neq > 0) {
        if (mixing_coefficient == 0.0) { // faster (less tensor operations)
          ifcs(ineq_cluster)(indices[i]) = corrected_ifc/((double) neq);
        } else {
          ifcs(ineq_cluster)(indices[i]) *= mixing_coefficient;
          ifcs(ineq_cluster)(indices[i]) += (1 - mixing_coefficient) * corrected_ifc/((double) neq);
        }
      }
    }
  }
}

//getCorrectionTerms//////////////////////////////////////////////////////////
// Calculates the correction term for each IFC.
aurostd::xtensor<double>
    AnharmonicIFCs::getCorrectionTerms(vector<int> atoms,
                                       const aurostd::xtensor<double>& ifcs,
                                       const aurostd::xtensor<double>& dev_from_zero,
                                       const aurostd::xtensor<double>& abssum) {
  atoms[0] = clst.sc2pcMap[atoms[0]];  // transfer to pcell
  aurostd::xtensor<double> correction_terms = ifcs(atoms);
  uint natoms = atoms.size();
  vector<int> reduced_cluster(natoms - 1);
  for (uint at = 0; at < natoms; at++) {
    reduced_cluster[at] = atoms[at];
  }
  aurostd::xtensor<double> correction = abs(ifcs(atoms));
  for (int crt = 0; crt < clst.nifcs; crt++) {
    if (abssum(reduced_cluster)(cart_indices[crt]) != 0.0) {
      double factor = dev_from_zero(reduced_cluster)(cart_indices[crt]);
      factor /= abssum(reduced_cluster)(cart_indices[crt]);
      correction(cart_indices[crt]) *= factor;
    } else {
      correction(cart_indices[crt]) = 0.0;
    }
  }
  correction_terms -= correction;
  return correction_terms;
}

// END Iterations

}  // namespace apl

/********************************* FILE I/O *********************************/

namespace apl {

// BEGIN Write files
//writeIFCsToFile/////////////////////////////////////////////////////////////
// Writes the AnharmonicIFCs object and minimal structure information to an
// XML file.
void AnharmonicIFCs::writeIFCsToFile(const string& filename) {
  stringstream output;
  // Header
  output << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
  output << "<anharmonicifcs>" << std::endl;

  output << writeParameters();
  output << writeIFCs();
  output << "</anharmonicifcs>" << std::endl;
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_IFCS_ERR_PREFIX_ + "writeIFCsToFile";
    string message = "Could not write tensor to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeParameters/////////////////////////////////////////////////////////////
// Writes the calculation parameters and minimal structure information to the
// XML file.
string AnharmonicIFCs::writeParameters() {
  stringstream parameters;
  string tab = " ";

  // Info about calculation run
  parameters << tab << "<generator>" << std::endl;
  string time = aflow_get_time_string();
  if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
  parameters << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
  parameters << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_;
  parameters << "\" type=\"Fletcher32\">" << std::hex << getFileCheckSum("./" + _AFLOWIN_ + "");
  parameters  << "</i>" << std::endl;
  parameters << tab << "</generator>" << std::endl;

  // Distortion magnitude
  parameters << tab << "<distortion_magnitude units=\"Angstrom\" cs=\"cartesian\">";
  parameters << distortion_magnitude << "</distortion_magnitude>" << std::endl;

  //Order
  parameters << tab << "<order>" << order << "</order>" << std::endl;

  // Iteration parameters
  parameters << tab << "<iteration>" << std::endl;
  // std::dec prevents hexadecimal output
  parameters << tab << tab << "<max_iter>" << std::dec << max_iter << "</max_iter>" << std::endl;
  parameters << tab << tab << "<mixing_coefficient>";
  parameters << std::setprecision(8) << mixing_coefficient;
  parameters << "</mixing_coefficient>" << std::endl;
  parameters << tab << tab << "<sumrule_threshold>";
  parameters << std::setprecision(15) << sumrule_threshold;
  parameters << "</sumrule_threshold>" << std::endl;
  parameters << tab << "</iteration>" << std::endl;

  // Structure
  parameters << tab << "<structure units=\"Angstrom\" cs=\"fractional\">" << std::endl;
  parameters << tab << tab << "<varray name=\"pcell lattice\">" << std::endl;
  for (int i = 1; i < 4; i++) {
    parameters << tab << tab << tab << "<v>";
    for (int j = 1; j < 4; j++) {
      parameters << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      parameters << setprecision(8);
      parameters << setw(15) << clst.pcell.lattice[i][j];
    }
    parameters << "</v>" << std::endl;
  }
  parameters << tab << tab << "</varray>" << std::endl;
  parameters << tab << tab << "<varray name=\"positions\">" << std::endl;
  for (uint i = 0; i < clst.pcell.atoms.size(); i++) {
    int t = clst.pcell.atoms[i].type;
    parameters << tab << tab << tab << "<v species=\"" << clst.pcell.species[t] << "\">";
    for (int j = 1; j < 4; j++) {
      parameters << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      parameters << setprecision(8);
      parameters << setw(15) << clst.pcell.atoms[i].fpos[j];
    }
    parameters << "</v>" << std::endl;
  }
  parameters << tab << tab << "</varray>" << std::endl;
  parameters << tab << tab << "<varray name=\"supercell\">" << std::endl;
  parameters << tab << tab << tab << "<v>";
  for (int i = 1; i < 4; i++) {
    parameters << tab << clst.sc_dim[i];
  }
  parameters << "</v>" << std::endl;
  parameters << tab << tab << "</varray>" << std::endl;
  parameters << tab << "</structure>" << std::endl;
  return parameters.str();
}

//writeIFCs///////////////////////////////////////////////////////////////////
// Writes the force constants part of the XML file.
string AnharmonicIFCs::writeIFCs() {
  stringstream ifcs;
  string tab = " ";
  int precision = 15;
  int extra = 5; // first digit, decimal point, minus sign + 2 spaces
  double max_ifc = max(abs(force_constants));
  int width = (int) log10(max_ifc);
  if (width < 0) {
    width = precision + extra;
  } else {
    width += precision + extra;
  }

  ifcs << tab << "<force_constants>" << std::endl;
  uint natoms_pcell = clst.pcell.atoms.size();
  uint natoms_scell = clst.scell.atoms.size();
  for (uint a = 0; a < natoms_pcell; a++) {
    int pcat = clst.pc2scMap[a];
    xcombos ats(natoms_scell, order - 1, 'E', true);
    while (ats.increment()) {
      vector<int> atoms = ats.getCombo();
      ifcs << tab << tab << "<varray atoms=\"" << pcat;
      for (uint at = 0; at < atoms.size(); at++) {
        ifcs << " " << atoms[at];
      }
      ifcs << "\">" << std::endl;
      xcombos crt(3, order - 2, 'E', true);
      while (crt.increment()) {
        vector<int> cart_ind = crt.getCombo();
        ifcs << tab << tab << tab << "<varray slice=\"" << cart_ind[0];
        for (uint c = 1; c < cart_ind.size(); c++) {
          ifcs << " " << cart_ind[c];
        }
        ifcs << "\">" << std::endl;
        for (int i = 0; i < 3; i++) {
          ifcs << tab << tab << tab << tab << "<v>";
          for (int j =0; j < 3; j++) {
            ifcs << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
            ifcs << setprecision(precision);
            ifcs << setw(width) << force_constants[a](atoms)(cart_ind)[i][j];
          }
          ifcs << "</v>" << std::endl;
        }
        ifcs << tab << tab << tab << "</varray>" << std::endl;
      }
      ifcs << tab << tab << "</varray>" << std::endl;
    }
  }
  ifcs << tab << "</force_constants>" << std::endl;
  return ifcs.str();
}

// END Write files

// BEGIN Read files
//readIFCsFromFile////////////////////////////////////////////////////////////
// Reads an AnharmonicIFCs object from an XML file.
void AnharmonicIFCs::readIFCsFromFile(const string& filename) {
  // Open file and handle exceptions
  string function = _AAPL_IFCS_ERR_PREFIX_ + "readIFCsFromFile";
  stringstream message;

  if (!aurostd::EFileExist(filename) && !aurostd::FileExist(filename)) {
    message << "Could not open file " << filename << ". File not found.";
    throw xerror(function, message, _FILE_NOT_FOUND_);
  }
  vector<string> vlines;
  aurostd::efile2vectorstring(filename, vlines);
  if (vlines.size() == 0) {
    message << "Cannot open file " << filename << ". File empty or corrupt.";
    throw xerror(function, message, _FILE_CORRUPT_);
  }

  // Start reading
  uint line_count = 0;
  string line = vlines[line_count++];

  // Check that this is a valid xml file
  if (line.find("xml") == string::npos) {
    message << "File is not a valid xml file.";
    throw xerror(function, message, _FILE_WRONG_FORMAT_);
  }

  // Check if xml file can be used to read anharmonic IFCs
  if (checkCompatibility(line_count, vlines)) {
    force_constants = readIFCs(line_count, vlines);
  } else {
    message << "The settings in the hibernate file and the aflow.in file are incompatible.";
    throw xerror(function, message, _RUNTIME_ERROR_);
  }
}

//checkCompatibility//////////////////////////////////////////////////////////
// Checks if the hibernate XML file is compatible with the aflow.in file. If
// the checksum in the XML file is the same as the checksum of the aflow.in
// file, then the parameters are the same. If not, the function checks if the
// parameters relevant for the IFC calculation (supercell, order, calculation
// parameters) are the same. This prevents the anharmonic IFCs from being
// recalculated when only post-processing parameters are changed.
bool AnharmonicIFCs::checkCompatibility(uint& line_count, 
                                        const vector<string>& vlines) {
  string function = _AAPL_IFCS_ERR_PREFIX_ + "checkCompatibility";
  string line;
  stringstream message;
  bool compatible = true;
  int t;
  vector<string> tokens;
  uint vsize = vlines.size();

  // Compare checksum
  while (true) {
    if (line_count == vsize) {
      message << "Checksum not found in hibernate file.";
      throw xerror(function, message, _FILE_CORRUPT_);
    }
    line = vlines[line_count++];
    if (line.find("checksum") != string::npos) {
      break;
    }
  }

  t = line.find_first_of(">") + 1;
  tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
  if (strtoul(tokens[0].c_str(), NULL, 16) != getFileCheckSum("./" + _AFLOWIN_ + "")) {
    message << "The " << _AFLOWIN_ << " file has been changed from the hibernated state. ";

    tokens.clear();
    // Compare calculation parameters
    //// distortion_magnitude
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find distortion_magnitude tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("distortion_magnitude") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        double dist_mag = aurostd::string2utype<double>(tokens[0]);
        tokens.clear();
        if (dist_mag != distortion_magnitude) {
          message << "Hibernate file and aflow.in have different distortion magnitudes. ";
          compatible = false;
        }
        break;
      }
    }

    //// order
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find order tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("order") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        int ord = aurostd::string2utype<int>(tokens[0]);
        tokens.clear();
        if (ord != order) {
          message << "Hibernate file and aflow.in have different IFC order. ";
          compatible = false;
        }
        break;
      }
    }

    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find iteration tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("iteration") != string::npos) {
        break;
      }
    }
    //// max_iter
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find max_iter tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("max_iter") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        int iter = aurostd::string2utype<int>(tokens[0]);
        tokens.clear();
        if (iter != max_iter) {
          message << "Hibernate file and aflow.in have different max. number of iterations. ";
          compatible = false;
        }
        break;
      }
    }

    //// mixing_coefficient
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find mixing_coefficient tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("mixing_coefficient") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        double mix = aurostd::string2utype<bool>(tokens[0]);
        tokens.clear();
        if (mix != mixing_coefficient) {
          message << "Hibernate file and aflow.in have different mixing coefficients. ";
          compatible = false;
        }
        break;
      }
    }

    //// sumrule_threshold
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find sumrule_threshold tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("sumrule_threshold") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        double thresh = aurostd::string2utype<double>(tokens[0]);
        tokens.clear();
        if (thresh != sumrule_threshold) {
          message << "Hibernate file and aflow.in have different convergence criteria. ";
          compatible = false;
        }
        break;
      }
    }

    // Compare supercells
    //// Lattice
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find structure tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("structure") != string::npos) {
        break;
      }
    }

    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find primitive lattice vectors. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("varray name=\"pcell lattice\"") != string::npos) {
        xmatrix<double> latt(3, 3);
        for (int i = 1; i < 4; i++) {
          line = vlines[line_count++];
          if (line_count == vsize) {
            message << "pcell lattice tag is corrupt. ";
          }
          t = line.find_first_of(">") + 1;
          tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
          for (int j = 1; j < 3; j++) {
            latt(i, j) = aurostd::string2utype<double>(tokens[j - 1]);
          }
          tokens.clear();
        }
        line = vlines[line_count++];
        if (line.find("</varray>") == string::npos) {
          message << "pcell lattice tag is corrupt. ";
          compatible = false;
        } else {
          break;
        }
        if (latt != clst.pcell.lattice) {
          message << "Hibernate file and aflow.in do not have the same lattice. ";
          compatible = false;
        }
      }
    }

    //// Atomic positions
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find atomic positions. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("varray name=\"positions\"") != string::npos) {
        // First check for tag corruption and extract everything
        vector<string> species;
        vector<xvector<double> > positions;
        while (compatible) {
          if (line_count == vsize) {
            message << "positions tag is corrupt. ";
            compatible = false;
          }
          line = vlines[line_count++];
          if (line.find("species=\"") != string::npos) {
            // Extract species
            t = line.find_first_of("\"") + 1;
            tokenize(line.substr(t, line.find_last_of("\"") - t), tokens, string(" "));
            species.push_back(tokens[0]);
            tokens.clear();
            // Extract positions
            t = line.find_first_of(">") + 1;
            tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            xvector<double> fpos(3);
            for (int i = 1; i < 4; i++) {
              fpos(i) = aurostd::string2utype<double>(tokens[i-1]);
            }
            positions.push_back(fpos);
            tokens.clear();
          } else if (line.find("/varray") != string::npos) {
            break;
          } else {
            message << "positions tag is corrupt. ";
            compatible = false;
          }
        }
        // Now compare with the primitive cell from the aflow.in file
        uint nspecies = species.size();
        if (nspecies == clst.pcell.atoms.size()) {
          for (uint sp = 0; sp < nspecies; sp++) {
            int type = clst.pcell.atoms[sp].type;
            string spec = clst.pcell.species[type];
            xvector<double> pos = clst.pcell.atoms[sp].fpos;
            if (species[sp] != spec) {
              message << "The structures in the hibernate file and aflow.in ";
              message << "have different species. ";
              compatible = false;
              sp = nspecies;
            } else if (positions[sp] != pos) {
              message << "The structures in the hibernate file and aflow.in ";
              message << "have different atomic positions.";
              compatible = false;
              sp = nspecies;
            }
          }
        } else {
          message << "The structures in the hibernate file and in aflow.in ";
          message << "do not have the same number of atoms. ";
          compatible = false;
        }
        break;
      }
    }

    //// Supercell dimensions
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find supercell dimensions. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("varray name=\"supercell\"") != string::npos) {
        line = vlines[line_count++];
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        xvector<int> sc(3);
        for (int i = 1; i < 4; i++) {
          sc(i) = aurostd::string2utype<double>(tokens[i-1]);
        }
        if (sc == clst.sc_dim) {
          break;
        } else {
          message << "The supercells in the hibernate file and in aflow.in ";
          message << "have different dimensions. ";
          compatible = false;
        }
      }
    }

    if (compatible) {
      message << "The relevant settings appear to be the same, ";
      message << "so the anharmonic IFCs will not be recalculated. ";
      message << "Make sure that the changes do not impact the IFCs.";
    } else {
      message << "Anharmonic IFCs need to be recalculated.";
    }
    _logger << apl::warning << "apl::AnharmonicIFCs::readIFCsfromFile(); ";
    _logger << message.str() << apl::endl;
  }
  return compatible;
}

//readIFCs////////////////////////////////////////////////////////////////////
// Reads the IFCs from the hibernate file.
aurostd::xtensor<double> AnharmonicIFCs::readIFCs(uint& line_count,
                                                  const vector<string>& vlines) {
  string function = _AAPL_IFCS_ERR_PREFIX_ + "readIFCs";
  string line, message;
  vector<int> atoms(order), slice(order - 2);
  vector<string> tokens;
  int t;
  uint vsize = vlines.size();

  // Initialize tensor
  uint natoms_pcell, natoms_scell;
  natoms_pcell = clst.pcell.atoms.size();
  natoms_scell = clst.scell.atoms.size();
  vector<int> tensor_shape(2 * order);
  tensor_shape[0] = natoms_pcell - 1;  // The first atom has to be in the primitive cell
  for (int i = 1; i < order; i++) {
    tensor_shape[i] = natoms_scell - 1;
  }
  for (int i = order; i < 2 * order; i++) {
    tensor_shape[i] = 2;
  }
  vector<int> zeros(2 * order, 0); // Make tensor use zero-based indices
  aurostd::xtensor<double> ifcs(tensor_shape, zeros);

  // Find force_constants tag
  while (true) {
    if (line_count == vsize) {
      message = "force_constants tag not found.";
      throw xerror(function, message, _FILE_CORRUPT_);
    }
    line = vlines[line_count++];
    if (line.find("force_constants") != string::npos) {
      break;
    }
  }

  // Read IFCs
  while (line.find("/force_constants") == string::npos) {
    if (line_count == vsize) {
      message = "force_constants tag incomplete.";
      throw xerror(function, message, _FILE_CORRUPT_);
    }
    line = vlines[line_count++];
    if (line.find("atoms") != string::npos) {
      t = line.find_first_of("\"") + 1;
      tokenize(line.substr(t, line.find_last_of("\"") - t), tokens, string(" "));
      for (int at = 0; at < order; at++) {
        int atom = aurostd::string2utype<int>(tokens[at]);
        if (at == 0) {
          atom = clst.sc2pcMap[atom];
        }
        atoms[at] = atom;
      }
      tokens.clear();
    } else if (line.find("slice") != string::npos) {
      t = line.find_first_of("\"") + 1;
      tokenize(line.substr(t, line.find_last_of("\"") - t), tokens, string(" "));
      for (uint sl = 0; sl < slice.size(); sl++) {
        slice[sl] = aurostd::string2utype<int>(tokens[sl]);
      }
      tokens.clear();
      // If a slice is found, populate tensor
      for (int i = 0; i < 3; i++) {
        line = vlines[line_count++];
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        for (int j = 0; j < 3; j++) {
          double val = aurostd::string2utype<double>(tokens[j]);
          ifcs(atoms)(slice)[i][j] = val; 
        }
        tokens.clear();
      }
    }
  }

  return ifcs;
}

// END Read files

} // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow Marco Esters - Duke University 2018                *
// *                                                                         *
// ***************************************************************************
