// ***************************************************************************
// *                                                                         *
// *              AFlow JOSE J. PLATA   - Duke University 2017               *
// *                                                                         *
// ***************************************************************************
// pairs_contrib_jose.cpp
// written by JJPR 2013
// 2013: jose.plata@duke.edu
//
// Based on Jesus Carrete's thirdorder.py and ShengBTE (DOI: 10.1016/j.cpc.2014.02.015).
//
//
//   This class stores and computes the inequivalent pairs and distortions.
//   To speed up the process:
//          1) Inequivalent pairs should contain one inequivalent atom (See constructor)
//          2) Before storing forces we just need if they are inequivalent or not
//          3) Inequivalent distoritions are just computed for inequivalent pairs
//
// ***************************************************************************
//  #include <limits>

#include "aflow_apl.h"

#define _isnegative(a) (a < MIN_EIGEN_TRESHOLD) ? true : false

//#undef AFLOW_APL_MULTITHREADS_ENABLE

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

using namespace std;

namespace apl {
//Constructor////PAIR////////////////////////////////////////////////////
_pair::_pair() {
  is_inequivalent = FALSE;
}
//Constructor////IPAIR////////////////////////////////////////////////////
_ipair::_ipair() {
}
//Constructor////ITRIPLET////////////////////////////////////////////////////
_itrip::_itrip() {
}
//Destructor/////PAIR///////////////////////////////////////////////////
_pair::~_pair() {
  free();
}
void _pair::free() {
}
//Destructor/////IPAIR///////////////////////////////////////////////////
_ipair::~_ipair() {
  free();
}
void _ipair::free() {
}
//Destructor/////ITRIPLET///////////////////////////////////////////////////
_itrip::~_itrip() {
  free();
}
void _itrip::free() {
}
// Copy constrcutor////PAIR/////////////////////////////////////////////
const _pair& _pair::operator=(const _pair& b) {  // operator=
                                                 //CO - START
  if (this != &b) {
    indexA = b.indexA;
    indexB = b.indexB;
    is_inequivalent = b.is_inequivalent;
    equivalent = b.equivalent;
    tri = b.tri;
    tri_is_eq = b.tri_is_eq;
    tri_peq = b.tri_peq;
    tri_teq = b.tri_teq;
    tri_sym = b.tri_sym;
    tri_perm = b.tri_perm;
  }
  return *this;
  //CO - END
}
// Copy constrcutor////IPAIR/////////////////////////////////////////////
const _ipair& _ipair::operator=(const _ipair& b) {  // operator=
                                                    //CO - START
  if (this != &b) {
    std::copy(&b.dist[0][0], &b.dist[0][0] + 6 * 6, &dist[0][0]);
    std::copy(&b.dist_3x3[0][0], &b.dist_3x3[0][0] + 3 * 3, &dist_3x3[0][0]);
    std::copy(&b.dist_in[0][0], &b.dist_in[0][0] + 6 * 6, &dist_in[0][0]);
    eq = b.eq;
    itriplets = b.itriplets;
    itriplets_exist = b.itriplets_exist;
  }
  return *this;
  //CO - END
}
// Copy constrcutor////ITRIPLET/////////////////////////////////////////////
const _itrip& _itrip::operator=(const _itrip& b) {  // operator=
                                                    //CO - START
  if (this != &b) {
    tri_in = b.tri_in;
    tri_eq = b.tri_eq;
    p_eq = b.p_eq;
    tri_eq2 = b.tri_eq2;
    p_eq2 = b.p_eq2;
    std::copy(&b.dist_tri[0][0][0], &b.dist_tri[0][0][0] + 3 * 3 * 3, &dist_tri[0][0][0]);
    //dist_tri_eq=b.dist_tri_eq;
    pool = b.pool;
    pool2 = b.pool2;
    triplets = b.triplets;
    sym = b.sym;
    perm = b.perm;
    rot_transform_aux = b.rot_transform_aux;
    rot_independent = b.rot_independent;
    count_independent = b.count_independent;
    val_independent = b.val_independent;
    valtmp_independent = b.valtmp_independent;
    rot_transform = b.rot_transform;
    rot_transform_array = b.rot_transform_array;
  }
  return *this;
  //CO - START
}

// Constructor Triplets BTE ////////////////////////////////////
_triBTE::_triBTE() {
}
//Destructor///// Triplets BTE ///////////////////////////////////////////////////
_triBTE::~_triBTE() {
  free();
}
void _triBTE::free() {
}
//Copy constrcutor///// Triplets BTE ///////////////////////////////////////////////////
const _triBTE& _triBTE::operator=(const _triBTE& b) {  // operator=
  //CO - START
  if (this != &b) {
    v1 = b.v1;
    v2 = b.v2;
    index = b.index;
    iat = b.iat;
    ifcs = b.ifcs;
  }
  return *this;
  //CO - END
}

///Constructor //// StrPairs/////////////////////////////////////////////////////////////

StrPairs::StrPairs(const xstructure a, const xstructure pc, const vector<int>& _pc2scMap, const vector<int>& _sc2pcMap, int SHELL, double& CUTOFF, Logger& l) : _logger(l) {
  // Logger ------------------------------------------------------------------------
  _logger << "PAIRS: Creating StrPairs for thermal conductivity calculation." << apl::endl;
  // Create table--------------------------------------------------------------------

  pc2scMap = _pc2scMap;
  sc2pcMap = _sc2pcMap;
  CUTOFF = getShellDist(a, SHELL, CUTOFF);  // Defining cutoff
  cutoff = CUTOFF;                          // Storing Cutoff

  cells = a.atoms.size() / pc.atoms.size();
  Sc_atoms = a.atoms.size();
  xmatrix<double> lattice = a.lattice;
  xvector<double> at1v, at2v, at3v;

  // Creating Pairs using Cutoff
  _pair pair;
  int counter = 0;
  for (uint at1pc = 0; at1pc < pc.atoms.size(); at1pc++) {
    uint at1sc = at1pc * cells;
    at1v = a.atoms[at1sc].cpos;
    for (uint at2sc = 0; at2sc < a.atoms.size(); at2sc++) {
      at2v = a.atoms[at2sc].cpos;
      double dist = getMinDist(lattice, at1v, at2v);
      if (dist < CUTOFF) {
        test(lattice, at1v, at2v, dist);
        pair.indexA = at1pc;
        pair.indexB = at2sc;
        pair.is_inequivalent = TRUE;
        pair.equivalent = counter;
        pairs.push_back(pair);
        counter = counter + 1;
      }
    }
  }

  // Clear things ----------------------------------------------------------------
  for (uint i = 0; i < ipairs.size(); i++)
    ipairs.clear();
  //clear();
}

//Constructor // StrPairs/////////////////////////////////////////////////////////////////////////

StrPairs::StrPairs(Logger& l) : _logger(l) {
}

StrPairs::StrPairs(const StrPairs& that) : _logger(that._logger) {
  *this = that;
}

//Copy Constructor // StrPairs/////////////////////////////////////////////////////////////////////////

const StrPairs& StrPairs::operator=(const StrPairs& that) {
  //CO - START
  if (this != &that) {
    _logger = that._logger;
    a = that.a;
    pcell = that.pcell;
    _pc2scMap = that._pc2scMap;
    _sc2pcMap = that._sc2pcMap;
    pc2scMap = that.pc2scMap;
    sc2pcMap = that.sc2pcMap;
    cells = that.cells;
    cutoff = that.cutoff;
    Sc_atoms = that.Sc_atoms;
    ortoDistorsions = that.ortoDistorsions;
    pairs = that.pairs;
    ipairs = that.ipairs;
    pairsF = that.pairsF;
    ipairsF = that.ipairsF;
    itriplets = that.itriplets;
    tri_index = that.tri_index;
    tri_vec = that.tri_vec;
    permut = that.permut;
    id_equi = that.id_equi;
    rot = that.rot;
    rot2 = that.rot2;
    nonzero = that.nonzero;
  }
  return *this;
  //CO - END
}

//Destructor // StrPairs/////////////////////////////////////////////////////////////////////////

StrPairs::~StrPairs() {
  clear();
}

//Clear  ///////////////////////////////////////////////////////////////////////////

void StrPairs::clear() {
}

//Build Str /////////////////////////////////////////////////////////////////////////
void StrPairs::build(xstructure pc, xstructure b) {
  a = b;
  pcell = pc;

  // EPS: WARNING!  THINK CAREFULLY!
  double eps = 0.001;

  // Look for inequivalent pairs-------------------------------------------------------
  _logger << "PAIRS: Looking for inequivalent pairs." << apl::endl;

  _atom tatom, tatom2;
  for (uint i = 0; i < pairs.size(); i++) {
    bool update = false;
    int at1sc = pairs[i].indexA * cells;
    int at2sc = pairs[i].indexB;
    for (uint j = i + 1; j < pairs.size(); j++) {
      int at3sc = pairs[j].indexA * cells;
      int at4sc = pairs[j].indexB;
      for (uint fg = 0; fg < a.fgroup.size(); fg++) {
        tatom = SYM::ApplyAtom(a.atoms[at1sc], a.fgroup[fg], a);
        tatom2 = SYM::ApplyAtom(a.atoms[at2sc], a.fgroup[fg], a);
        if (SYM::AtomFPOSMatch(tatom.fpos, a.atoms[at3sc].fpos, a.c2f, a.f2c, update, a.sym_eps) && SYM::AtomFPOSMatch(tatom2.fpos, a.atoms[at4sc].fpos, a.c2f, a.f2c, update, a.sym_eps)) {
          pairs[j].is_inequivalent = FALSE;
          pairs[j].equivalent = pairs[i].equivalent;
          fg = a.fgroup.size();
          j = pairs.size();
        } else {
          if (SYM::AtomFPOSMatch(tatom.fpos, a.atoms[at4sc].fpos, a.c2f, a.f2c, update, a.sym_eps) && SYM::AtomFPOSMatch(tatom2.fpos, a.atoms[at3sc].fpos, a.c2f, a.f2c, update, a.sym_eps)) {
            pairs[j].is_inequivalent = FALSE;
            pairs[j].equivalent = pairs[i].equivalent;
            fg = a.fgroup.size();
            j = pairs.size();
          }
        }
      }
    }
  }

  // Crate the list of inequivalent pairs-----------------------
  _ipair ipair;
  for (uint ipa = 0; ipa < pairs.size(); ipa++) {  //That could be useful for quartic forces too...
    if (pairs.at(ipa).is_inequivalent) {
      ipair.eq.push_back(ipa);
      for (uint i = 0; i < 6; i++) {
        for (uint j = 0; j < 6; j++) {
          ipair.dist[i][j] = (i * 6) + j;
          ipair.dist_in[i][j] = TRUE;
          ipair.itriplets_exist = FALSE;
        }
      }
      for (uint ipa2 = ipa + 1; ipa2 < pairs.size(); ipa2++) {
        if (pairs[ipa2].equivalent == (int)ipa) {
          ipair.eq.push_back(ipa2);
        }
      }
      ipairs.push_back(ipair);
      ipair.eq.clear();
    }
  }

  // Evaluate inequivalent distortion just in inequivalent pairs ------------------------
  // Distortions vectors----------------------------------------------------------------
  xvector<double> testVec(3);

  // Cartesian vectors x,y,z--------------------------------------

  testVec(1) = 1.0;
  testVec(2) = 0.0;
  testVec(3) = 0.0;
  ortoDistorsions.push_back(testVec);
  testVec(1) = 0.0;
  testVec(2) = 1.0;
  testVec(3) = 0.0;
  ortoDistorsions.push_back(testVec);
  testVec(1) = 0.0;
  testVec(2) = 0.0;
  testVec(3) = 1.0;
  ortoDistorsions.push_back(testVec);

  //Negative direction too---------------------------------------

  for (uint i = 0; i < 3; i++) {
    testVec = ortoDistorsions[i] * (-1);
    ortoDistorsions.push_back(testVec);
  }

  //Initialize 6x6 matrix--------------------------------------------

  for (uint ipa = 0; ipa < ipairs.size(); ipa++) {
  }
  _logger << "PAIRS: Looking for inequivalent distortions for inequivalent pairs." << apl::endl;
  //Look for inequivalent directions just in inequivalent pairs---------------------------
  //_atom tatom, tatom2;
  for (uint ipa = 0; ipa < ipairs.size(); ipa++) {
    bool update = false;
    uint iat1 = pairs[ipairs[ipa].eq[0]].indexA * cells;
    uint iat2 = pairs[ipairs[ipa].eq[0]].indexB;
    for (uint index1 = 0; index1 < 36; index1++) {
      uint i = index1 / 6;
      uint j = index1 % 6;
      for (uint index2 = index1 + 1; index2 < 36; index2++) {
        uint k = index2 / 6;
        uint l = index2 % 6;
        for (uint fg = 0; fg < a.fgroup.size(); fg++) {
          tatom = SYM::ApplyAtom(a.atoms[iat1], a.fgroup[fg], a);
          tatom2 = SYM::ApplyAtom(a.atoms[iat2], a.fgroup[fg], a);
          if (SYM::AtomFPOSMatch(tatom.fpos, a.atoms[iat1].fpos, a.c2f, a.f2c, update, a.sym_eps) && SYM::AtomFPOSMatch(tatom2.fpos, a.atoms[iat2].fpos, a.c2f, a.f2c, update, a.sym_eps)) {
            xvector<double> vec1(3);
            xvector<double> vec2(3);
            vec1 = a.fgroup[fg].Uc * ortoDistorsions[i];
            vec2 = a.fgroup[fg].Uc * ortoDistorsions[j];
            if (iat1 == iat2) {  // Special case
              if (aurostd::modulus(ortoDistorsions[k] + ortoDistorsions[l] - vec1 - vec2) < eps) {
                ipairs[ipa].dist[k][l] = ipairs[ipa].dist[i][j];
                ipairs[ipa].dist_in[k][l] = FALSE;
                fg = a.fgroup.size();
                index2 = 36;
              }
            } else {
              if (aurostd::modulus(ortoDistorsions[k] - vec1) < eps && aurostd::modulus(ortoDistorsions[l] - vec2) < eps) {
                ipairs[ipa].dist[k][l] = ipairs[ipa].dist[i][j];
                ipairs[ipa].dist_in[k][l] = FALSE;
                fg = a.fgroup.size();
                index2 = 36;
              }
            }
          }
        }
      }
      if (iat1 == iat2 && (aurostd::modulus(ortoDistorsions[i] + ortoDistorsions[j]) < eps)) {
        ipairs[ipa].dist_in[i][j] = FALSE;
      }
    }
  }
}

//Triplets//////////////////////////////////////////////////////////////
// Include triplet in pairs class////////////////////////////////////////

void StrPairs::Triplets() {
  xmatrix<double> lattice = a.lattice;
  xvector<double> at1v, at2v, at3v;

  _logger << "TRIPLETS: Creating triplets." << apl::endl;
  //Create triplets
  bool tri_is_eq = TRUE;
  for (uint ipa = 0; ipa < pairs.size(); ipa++) {
    uint at1 = pairs[ipa].indexA * cells;
    uint at2 = pairs[ipa].indexB;
    at1v = a.atoms[at1].cpos;
    at2v = a.atoms[at2].cpos;
    uint counter = 0;
    for (uint i = 0; i < Sc_atoms; i++) {
      at3v = a.atoms[i].cpos;
      double dist1 = getMinDist(lattice, at2v, at3v);
      double dist2 = getMinDist(lattice, at3v, at1v);
      if (dist1 < cutoff && dist2 < cutoff) {
        pairs[ipa].tri.push_back(i);
        pairs[ipa].tri_is_eq.push_back(tri_is_eq);
        pairs[ipa].tri_peq.push_back(ipa);
        pairs[ipa].tri_teq.push_back(counter);
        pairs[ipa].tri_sym.push_back(0);
        pairs[ipa].tri_perm.push_back(0);
        counter = counter + 1;
      }
    }
  }

  SymUtils();

  // Comparing triplets, looking for inequivalent ones and their dependency
  int counter3 = 0;
  vector<vector<int> > list_itrip;
  vector<vector<double> > dum27(27, vector<double>(27, 0.));
  vector<int> list_itrip_p;
  for (uint kk3 = 0; kk3 < ipairs.size(); kk3++) {
    for (uint kk4 = 0; kk4 < ipairs[kk3].eq.size(); kk4++) {
      uint ipa = ipairs[kk3].eq[kk4];
      vector<int> triplet(3);
      triplet[0] = pairs[ipa].indexA * cells;
      triplet[1] = pairs[ipa].indexB;
      for (uint tri = 0; tri < pairs[ipa].tri.size(); tri++) {
        bool hola = TRUE;
        triplet[2] = pairs[ipa].tri[tri];
        counter3 += 1;
        for (int iperm = 0; iperm < 6; iperm++) {
          vector<int> triplet_perm(3);
          triplet_perm[0] = triplet[permut[iperm][0]];
          triplet_perm[1] = triplet[permut[iperm][1]];
          triplet_perm[2] = triplet[permut[iperm][2]];
          for (uint isym = 0; isym < pcell.fgroup.size(); isym++) {
            vector<int> triplet_sym(3);
            triplet_sym[0] = id_equi[isym][triplet_perm[0]];
            triplet_sym[1] = id_equi[isym][triplet_perm[1]];
            triplet_sym[2] = id_equi[isym][triplet_perm[2]];
            moveTri(triplet_sym);
            for (uint kk2 = 0; kk2 < list_itrip.size(); kk2++) {
              if (compareTri(list_itrip[kk2], triplet_sym)) {
                pairs[ipa].tri_is_eq[tri] = FALSE;
                pairs[ipa].tri_peq[tri] = list_itrip_p[kk2];
                pairs[ipa].tri_teq[tri] = list_itrip[kk2][2];
                pairs[ipa].tri_sym[tri] = isym;
                pairs[ipa].tri_sym[tri] = iperm;
                itriplets[kk2].triplets.push_back(triplet);
                itriplets[kk2].sym.push_back(isym);
                itriplets[kk2].perm.push_back(iperm);
                itriplets[kk2].rot_transform.push_back(dum27);
                itriplets[kk2].rot_transform_array.push_back(dum27);
                kk2 = list_itrip.size();
                kk2 = list_itrip.size();
                isym = pcell.fgroup.size();
                iperm = 6;
                hola = FALSE;
              }
            }
          }
        }
        if (hola) {
          list_itrip.push_back(triplet);
          list_itrip_p.push_back(ipa);
          _itrip itrip;
          itrip.triplets.push_back(triplet);
          itrip.sym.push_back(0);
          itrip.perm.push_back(0);
          itrip.rot_transform.push_back(dum27);
          itrip.rot_transform_array.push_back(dum27);
          itriplets.push_back(itrip);
        }
      }
    }
  }

  // We have [inequivalent]=O[perm]O[sym][equivalent]
  // but we want [equivalent]=O[perm]O[sym][inequivalent]
  for (uint i = 0; i < itriplets.size(); i++) {
    vector<int> triplet(3);
    triplet[0] = itriplets[i].triplets[0][0];
    triplet[1] = itriplets[i].triplets[0][1];
    triplet[2] = itriplets[i].triplets[0][2];
    for (uint j = 0; j < itriplets[i].triplets.size(); j++) {
      vector<int> triplet_eq(3);
      bool check = TRUE;
      triplet_eq[0] = itriplets[i].triplets[j][0];
      triplet_eq[1] = itriplets[i].triplets[j][1];
      triplet_eq[2] = itriplets[i].triplets[j][2];
      for (int iperm = 0; iperm < 6; iperm++) {
        vector<int> triplet_perm(3);
        triplet_perm[0] = triplet[permut[iperm][0]];
        triplet_perm[1] = triplet[permut[iperm][1]];
        triplet_perm[2] = triplet[permut[iperm][2]];
        for (uint isym = 0; isym < pcell.fgroup.size(); isym++) {
          vector<int> triplet_sym(3);
          triplet_sym[0] = id_equi[isym][triplet_perm[0]];
          triplet_sym[1] = id_equi[isym][triplet_perm[1]];
          triplet_sym[2] = id_equi[isym][triplet_perm[2]];
          moveTri(triplet_sym);
          if (compareTri(triplet_eq, triplet_sym)) {
            itriplets[i].sym[j] = isym;
            itriplets[i].perm[j] = iperm;
            iperm = 6;
            isym = pcell.fgroup.size();
            check = FALSE;
          }
        }
      }
      if (check) { cout << " Failure redoing equivalences " << std::endl; }
    }
  }

  for (uint kk = 0; kk < ipairs.size(); kk++) {
    for (uint ll = 0; ll < list_itrip_p.size(); ll++) {
      if (int(ipairs[kk].eq[0]) == list_itrip_p[ll]) {
        ipairs[kk].itriplets_exist = TRUE;
      }
    }
  }

  //Transformation array
  for (uint kk = 0; kk < itriplets.size(); kk++) {
    for (uint kk2 = 0; kk2 < itriplets[kk].triplets.size(); kk2++) {
      for (int ii = 0; ii < 27; ii++) {
        for (int jj = 0; jj < 27; jj++) {
          itriplets[kk].rot_transform[kk2][ii][jj] = rot[itriplets[kk].perm[kk2]][itriplets[kk].sym[kk2]][ii][jj];
        }
      }
    }
  }

  // Checking inequivalent components in inequivalent triplets
  vector<double> dummy(27, 0.0);
  for (uint kk = 0; kk < itriplets.size(); kk++) {
    vector<int> triplet(3);
    vector<vector<double> > coeffi;
    triplet[0] = itriplets[kk].triplets[0][0];
    triplet[1] = itriplets[kk].triplets[0][1];
    triplet[2] = itriplets[kk].triplets[0][2];
    int nnonzero = 0;
    //Permutations
    for (int iperm = 0; iperm < 6; iperm++) {
      vector<int> triplet_perm(3);
      triplet_perm[0] = triplet[permut[iperm][0]];
      triplet_perm[1] = triplet[permut[iperm][1]];
      triplet_perm[2] = triplet[permut[iperm][2]];
      //Symmetry
      for (uint isym = 0; isym < pcell.fgroup.size(); isym++) {
        vector<int> triplet_sym(3);
        triplet_sym[0] = id_equi[isym][triplet_perm[0]];
        triplet_sym[1] = id_equi[isym][triplet_perm[1]];
        triplet_sym[2] = id_equi[isym][triplet_perm[2]];
        // 1st atom must be in the original Primitive cell
        moveTri(triplet_sym);
        if (compareTri(triplet, triplet_sym)) {
          for (int indexijkprime = 0; indexijkprime < 27; indexijkprime++) {
            if (nonzero[iperm][isym][indexijkprime]) {
              coeffi.push_back(dummy);
              for (int ll = 0; ll < 27; ll++) {
                coeffi[nnonzero][ll] = rot2[iperm][isym][indexijkprime][ll];
              }
              nnonzero = nnonzero + 1;
            }
          }
        }
      }
    }
    int dim = std::max(nnonzero, 27);
    vector<vector<double> > coeffi_reduced(dim, vector<double>(27));
    for (int iaux = 0; iaux < nnonzero; iaux++) {
      for (int jaux = 0; jaux < 27; jaux++) {
        coeffi_reduced[iaux][jaux] = coeffi[iaux][jaux];
      }
    }
    gaussAAPL(kk, coeffi_reduced);
  }

  for (uint ii = 0; ii < itriplets.size(); ii++) {
    for (uint jj = 0; jj < itriplets[ii].triplets.size(); jj++) {
      for (int kk = 0; kk < 27; kk++) {
        for (uint ll = 0; ll < itriplets[ii].rot_independent.size(); ll++) {
          for (int iaux = 0; iaux < 27; iaux++) {
            itriplets[ii].rot_transform_array[jj][kk][ll] += itriplets[ii].rot_transform[jj][kk][iaux] * itriplets[ii].rot_transform_aux[iaux][ll];
          }
        }
      }
      for (int kk = 0; kk < 27; kk++) {
        for (int ll = 0; ll < 27; ll++) {
          if (abs(itriplets[ii].rot_transform_array[jj][kk][ll]) < 1e-12) {
            itriplets[ii].rot_transform_array[jj][kk][ll] = 0.;
          }
        }
      }
    }
  }

  for (uint pa = 0; pa < pairs.size(); pa++) {
    int at1pc = pairs[pa].indexA;
    int at1 = pairs[pa].indexA * cells;
    xvector<double> at1vc = a.atoms[at1].cpos;
    int at2 = pairs[pa].indexB;
    xvector<double> at2vc = a.atoms[at2].cpos;
    vector<xvector<double> > v2_eq;
    double dist = getMinDist(lattice, at1vc, at2vc);
    for (int o = -2; o <= 2; o++) {
      for (int p = -2; p <= 2; p++) {
        for (int q = -2; q <= 2; q++) {
          xvector<double> image = at2vc + o * lattice(1) + p * lattice(2) + q * lattice(3);
          if (abs(dist - aurostd::modulus(at1vc - image)) < 0.001) {
            v2_eq.push_back(image);
          }
        }
      }
    }
    for (uint i = 0; i < v2_eq.size(); i++) {
      for (uint tr = 0; tr < pairs[pa].tri.size(); tr++) {
        int at3 = pairs[pa].tri[tr];
        xvector<double> at3vc = a.atoms[at3].cpos;
        for (int o = -2; o <= 2; o++) {
          for (int p = -2; p <= 2; p++) {
            for (int q = -2; q <= 2; q++) {
              xvector<double> image2 = at3vc + o * lattice(1) + p * lattice(2) + q * lattice(3);
              if (aurostd::modulus(at1vc - image2) < cutoff && aurostd::modulus(v2_eq[i] - image2) < cutoff) {
                xvector<int> tri_kk;
                vector<xvector<double> > vec_kk;
                vec_kk.push_back(at1vc);
                vec_kk.push_back(v2_eq[i]);
                vec_kk.push_back(image2);
                tri_kk[1] = at1pc;
                tri_kk[2] = at2;
                tri_kk[3] = at3;
                tri_index.push_back(tri_kk);
                tri_vec.push_back(vec_kk);
              }
            }
          }
        }
      }
    }
    v2_eq.clear();
  }
}

////////////////////////////////////////////////////////////////////////////////////////

void StrPairs::calculateTiplets_MPI(const xstructure a, uint startIndex, uint endIndex) {
  for (uint kk = startIndex; kk < endIndex; kk++) {
    for (uint kk2 = 0; kk2 < ipairs[kk].eq.size(); kk2++) {
      uint ipa = ipairs[kk].eq[kk2];
      int at1sc = pairs[ipa].indexA * cells;
      int at2sc = pairs[ipa].indexB;
      for (uint i = 0; i < pairs[ipa].tri.size(); i++) {
        int at3sc = pairs[ipa].tri[i];
        for (uint kk3 = kk2; kk3 < ipairs[kk].eq.size(); kk3++) {
          bool EXITO = FALSE;
          uint ipa2 = ipairs[kk].eq[kk3];
          int at4sc = pairs[ipa2].indexA * cells;
          int at5sc = pairs[ipa2].indexB;
          if (ipa == ipa2) {
            for (uint j = i + 1; j < pairs[ipa2].tri.size(); j++) {
              int at6sc = pairs[ipa2].tri[j];
              for (uint fg = 0; fg < a.fgroup.size(); fg++) {
                if (MATCH3(at4sc, at5sc, at6sc, at1sc, at2sc, at3sc, fg)) {
                  pairs[ipa2].tri_is_eq[j] = FALSE;
                  pairs[ipa2].tri_peq[j] = pairs[ipa].tri_peq[i];
                  pairs[ipa2].tri_teq[j] = pairs[ipa].tri_teq[i];
                  pairs[ipa2].tri_sym[j] = fg;
                  fg = a.fgroup.size();
                  j = pairs[ipa2].tri.size();
                  EXITO = TRUE;
                }
              }
            }
          } else {
            for (uint j = 0; j < pairs[ipa2].tri.size(); j++) {
              int at6sc = pairs[ipa2].tri[j];
              for (uint fg = 0; fg < a.fgroup.size(); fg++) {
                if (MATCH3(at4sc, at5sc, at6sc, at1sc, at2sc, at3sc, fg)) {
                  pairs[ipa2].tri_is_eq[j] = FALSE;
                  pairs[ipa2].tri_peq[j] = pairs[ipa].tri_peq[i];
                  pairs[ipa2].tri_teq[j] = pairs[ipa].tri_teq[i];
                  pairs[ipa2].tri_sym[j] = fg;
                  fg = a.fgroup.size();
                  j = pairs[ipa2].tri.size();
                  EXITO = TRUE;
                }
              }
            }
          }
          if (EXITO) { kk3 = ipairs[kk].eq.size(); }
        }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////

void StrPairs::calculateTriDist_MPI(const xstructure a, uint startIndex, uint endIndex) {
  for (uint kk = startIndex; kk < endIndex; kk++) {
    uint ipa = ipairs[kk].eq[0];
    int at1sc = pairs[ipa].indexA * cells;
    int at2sc = pairs[ipa].indexB;
    for (uint t = 0; t < ipairs[kk].itriplets.size(); t++) {
      int at3sc = pairs[ipa].tri[ipairs[kk].itriplets[t].tri_in];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            for (int o = 0; o < 3; o++) {
              for (int p = 0; p < 3; p++) {
                for (int q = 0; q < 3; q++) {
                  if (((i * 9) + (j * 3) + k) < ((o * 9) + (p * 3) + q)) {
                    for (uint fg = 0; fg < a.fgroup.size(); fg++) {
                      int sign = 1;
                      if (MATCH3v(at1sc, at2sc, at3sc, at1sc, at2sc, at3sc, o, p, q, i, j, k, fg, sign)) {
                        ipairs[kk].itriplets[t].dist_tri[o][p][q] = sign * ipairs[kk].itriplets[t].dist_tri[i][j][k];
                        q = 3;
                        p = 3;
                        o = 3;
                        fg = a.fgroup.size();
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            if ((i * 9) + (j * 3) + k == ipairs[kk].itriplets[t].dist_tri[i][j][k]) { ipairs[kk].itriplets[t].dist_tri[i][j][k] = -27; }
          }
        }
      }
    }
  }
}

///GetAtom/////////////////////////////////////////////////////////////////////
// This function computed the atom of the supercell that is obtained when symmetry operator Op is applied to atom1

int StrPairs::getAtom(const xstructure& b, int atom1, int Op) {
  uint i;
  _atom tatom;
  _atom tatom2;
  xmatrix<double> mc2f = b.c2f;
  xmatrix<double> mf2c = b.f2c;
  double dsym_eps = b.sym_eps;
  tatom2 = b.atoms[atom1];
  for (i = 0; i < b.atoms.size(); i++) {
    bool update = false;
    tatom = SYM::ApplyAtom(b.atoms[i], b.fgroup[Op], b);
    if (SYM::AtomFPOSMatch(tatom2.fpos, tatom.fpos, mc2f, mf2c, update, dsym_eps)) {
      break;
    }
  }
  return i;
}

/////MATCH 2///////////////////
//// Check inequivalent pairs aplaying symmetry///////

bool StrPairs::MATCH2(int ref1, int ref2, int at1, int at2, int fg) {
  double eps = 0.00001;
  xmatrix<int> match(2, 2, 1, 1);
  int counter = 0;
  xvector<int> ref(2);
  ref[1] = ref1;
  ref[2] = ref2;
  xvector<int> at(2);
  at[1] = at1;
  at[2] = at2;
  for (uint i = 1; i < 3; i++) {
    for (uint j = 1; j < 3; j++) {
      match[i][j] = 0;
    }
  }

  _atom tatom;
  for (uint i = 1; i < 3; i++) {
    tatom = SYM::ApplyAtom(a.atoms[at[i]], a.fgroup[fg], a);
    for (uint j = 1; j < 3; j++) {
      if (match[1][j] + match[2][j] < 1) {
        if (aurostd::modulus(BringInCell(a.atoms[ref[j]].fpos - tatom.fpos)) < eps) {
          match[i][j] = 1;
          counter = counter + 1;
          j = 3;
        }
      }
    }
  }

  if (counter == 2) {
    return TRUE;
  } else {
    return FALSE;
  }
}

/////MATCH 3///////////////////
//// Check inequivalent pairs aplaying symmetry///////

bool StrPairs::MATCH3(int ref1, int ref2, int ref3, int at1, int at2, int at3, int fg) {
  double eps = 0.00001;
  xmatrix<int> match(3, 3, 1, 1);
  int counter = 0;
  xvector<int> ref(3);
  ref[1] = ref1;
  ref[2] = ref2;
  ref[3] = ref3;
  xvector<int> at(3);
  at[1] = at1;
  at[2] = at2;
  at[3] = at3;
  for (uint i = 1; i < 4; i++) {
    for (uint j = 1; j < 4; j++) {
      match[i][j] = 0;
    }
  }

  _atom tatom;
  for (uint i = 1; i < 4; i++) {
    tatom = SYM::ApplyAtom(a.atoms[at[i]], a.fgroup[fg], a);
    for (uint j = 1; j < 4; j++) {
      if (match[1][j] + match[2][j] + match[3][j] < 1) {
        if (aurostd::modulus(BringInCell(a.atoms[ref[j]].fpos - tatom.fpos)) < eps) {
          match[i][j] = 1;
          counter = counter + 1;
          j = 4;
        }
      }
    }
    if (match[i][1] + match[i][2] + match[i][3] == 0) { return FALSE; }
  }

  if (counter == 3) {
    return TRUE;
  } else {
    return FALSE;
  }
}

/////MATCH 2_2///////////////////
//// Check inequivalent pairs aplaying symmetry///////

bool StrPairs::MATCH2_2(int ref1, int ref2, int at1, int at2, int fg, xmatrix<int>& match) {
  double eps = 0.00001;
  int counter = 0;
  xvector<int> ref(2);
  ref[1] = ref1;
  ref[2] = ref2;
  xvector<int> at(2);
  at[1] = at1;
  at[2] = at2;
  for (uint i = 1; i < 3; i++) {
    for (uint j = 1; j < 3; j++) {
      match[i][j] = 0;
    }
  }

  _atom tatom;
  for (uint i = 1; i < 3; i++) {
    tatom = SYM::ApplyAtom(a.atoms[at[i]], a.fgroup[fg], a);
    for (uint j = 1; j < 3; j++) {
      if (match[1][j] + match[2][j] < 1) {
        if (aurostd::modulus(BringInCell(a.atoms[ref[j]].fpos - tatom.fpos)) < eps) {
          match[i][j] = 1;
          counter = counter + 1;
          j = 3;
        }
      }
    }
    if (match[i][1] + match[i][2] == 0) { return FALSE; }
  }

  if (counter == 2) {
    return TRUE;
  } else {
    return FALSE;
  }
}

/////MATCH 3_2///////////////////
//// Check inequivalent pairs aplaying symmetry///////

bool StrPairs::MATCH3_2(int ref1, int ref2, int ref3, int at1, int at2, int at3, int fg, xmatrix<int>& match) {
  double eps = 0.00001;
  int counter = 0;
  xvector<int> ref(3);
  ref[1] = ref1;
  ref[2] = ref2;
  ref[3] = ref3;
  xvector<int> at(3);
  at[1] = at1;
  at[2] = at2;
  at[3] = at3;
  for (uint i = 1; i < 4; i++) {
    for (uint j = 1; j < 4; j++) {
      match[i][j] = 0;
    }
  }

  _atom tatom;
  for (uint i = 1; i < 4; i++) {
    tatom = SYM::ApplyAtom(a.atoms[at[i]], a.fgroup[fg], a);
    for (uint j = 1; j < 4; j++) {
      if (match[1][j] + match[2][j] + match[3][j] < 1) {
        if (aurostd::modulus(BringInCell(a.atoms[ref[j]].fpos - tatom.fpos)) < eps) {
          match[i][j] = 1;
          counter = counter + 1;
          j = 4;
        }
      }
    }
    if (match[i][1] + match[i][2] + match[i][3] == 0) { return FALSE; }
  }

  if (counter == 3) {
    return TRUE;
  } else {
    return FALSE;
  }
}

bool StrPairs::MATCH3v2(int at1, int at2, int at3, int ref1, int ref2, int ref3, int co1, int co2, int co3, int fg, xmatrix<int> match, int& counter2) {
  xmatrix<int> match2(3, 3, 1, 1);
  double eps = 0.00001;
  int counter = 0;
  xvector<int> at, ref, co, coF;
  at[1] = at1;
  at[2] = at2;
  at[3] = at3;
  co[1] = co1;
  co[2] = co2;
  co[3] = co3;
  ref[1] = ref1;
  ref[2] = ref2;
  ref[3] = ref3;
  for (int i = 1; i < 4; i++) {
    for (uint j = 1; j < 4; j++) {
      coF[j] = coF[j] + match[i][j] * co[i];
    }
  }
  for (int i = 1; i < 4; i++) {
    for (uint j = 1; j < 4; j++) {
      match2[i][j] = 0;
    }
  }
  counter = 0;
  for (uint i = 1; i < 4; i++) {
    for (uint j = 1; j < 4; j++) {
      if (match2[1][j] + match2[2][j] + match2[3][j] < 1 && at[i] == at[j]) {
        if (aurostd::modulus(ortoDistorsions[ref[j]] - a.fgroup[fg].Uc * ortoDistorsions[coF[i]]) < eps) {
          match2[i][j] = 1;
          counter = counter + 1;
          counter2 = counter2 * 1;
          j = 4;
        } else {
          if (aurostd::modulus(ortoDistorsions[ref[j]] + a.fgroup[fg].Uc * ortoDistorsions[coF[i]]) < eps) {
            match2[i][j] = 1;
            counter = counter + 1;
            counter2 = counter2 * (-1);
            j = 4;
          }
        }
      }
    }
    if (match[i][1] + match[i][2] + match[i][3] == 0) { return FALSE; }
  }

  if (counter == 3) {
    return TRUE;
  } else {
    return FALSE;
  }
}

bool StrPairs::MATCH2v2(int at1r, int at2r, int at1, int at2, int ref1, int ref2, int co1, int co2, int fg, int& counter2) {
  double eps = 0.00001;
  xmatrix<int> match(2, 2, 1, 1);
  int counter = 0;
  xvector<int> at, atr, ref, co;
  at[1] = at1;
  at[2] = at2;
  atr[1] = at1r;
  atr[2] = at2r;
  co[1] = co1;
  co[2] = co2;
  ref[1] = ref1;
  ref[2] = ref2;

  for (int i = 1; i < 3; i++) {
    for (uint j = 1; j < 3; j++) {
      match[i][j] = 0;
    }
  }

  counter2 = 1;
  _atom tatom;
  for (uint i = 1; i < 3; i++) {
    tatom = SYM::ApplyAtom(a.atoms[at[i]], a.fgroup[fg], a);
    for (uint j = 1; j < 3; j++) {
      if (match[1][j] + match[2][j] < 1) {
        if (aurostd::modulus(BringInCell(a.atoms[atr[j]].fpos - tatom.fpos)) < eps) {
          if (aurostd::modulus(ortoDistorsions[ref[j]] - a.fgroup[fg].Uc * ortoDistorsions[co[i]]) < eps) {
            match[i][j] = 1;
            counter = counter + 1;
            counter2 = counter2 * (1);
            j = 3;
          } else {
            if (aurostd::modulus(ortoDistorsions[ref[j]] + a.fgroup[fg].Uc * ortoDistorsions[co[i]]) < eps) {
              match[i][j] = 1;
              counter = counter + 1;
              counter2 = counter2 * (-1);
              j = 3;
            }
          }
        }
      }
    }
    if (match[i][1] + match[i][2] == 0) { return FALSE; }
  }
  if (counter == 2) {
    return TRUE;
  } else {
    return FALSE;
  }
}

/////MATCH 3v///////////////////
//// Check inequivalent pairs aplaying symmetry///////

bool StrPairs::MATCH3v(int at1r, int at2r, int at3r, int at1, int at2, int at3, int ref1, int ref2, int ref3, int co1, int co2, int co3, int fg, int& counter2) {
  double eps = 0.00001;
  xmatrix<int> match(3, 3, 1, 1);
  int counter = 0;
  xvector<int> at, atr, ref, co;
  at[1] = at1;
  at[2] = at2;
  at[3] = at3;
  atr[1] = at1r;
  atr[2] = at2r;
  atr[3] = at3r;
  co[1] = co1;
  co[2] = co2;
  co[3] = co3;
  ref[1] = ref1;
  ref[2] = ref2;
  ref[3] = ref3;

  for (int i = 1; i < 4; i++) {
    for (uint j = 1; j < 4; j++) {
      match[i][j] = 0;
    }
  }

  counter2 = 1;
  _atom tatom;
  for (uint i = 1; i < 4; i++) {
    tatom = SYM::ApplyAtom(a.atoms[at[i]], a.fgroup[fg], a);
    for (uint j = 1; j < 4; j++) {
      if (match[1][j] + match[2][j] + match[3][j] < 1) {
        if (aurostd::modulus(BringInCell(a.atoms[atr[j]].fpos - tatom.fpos)) < eps) {
          if (aurostd::modulus(ortoDistorsions[ref[j]] - a.fgroup[fg].Uc * ortoDistorsions[co[i]]) < eps) {
            match[i][j] = 1;
            counter = counter + 1;
            counter2 = counter2 * (1);
            j = 4;
          } else {
            if (aurostd::modulus(ortoDistorsions[ref[j]] + a.fgroup[fg].Uc * ortoDistorsions[co[i]]) < eps) {
              match[i][j] = 1;
              counter = counter + 1;
              counter2 = counter2 * (-1);
              j = 4;
            }
          }
        }
      }
    }
    if (match[i][1] + match[i][2] + match[i][3] == 0) { return FALSE; }
  }
  if (counter == 3) {
    return TRUE;
  } else {
    return FALSE;
  }
}
/////MATCH 2v///////////////////
//// Check inequivalent pairs aplaying symmetry///////

bool StrPairs::MATCH2v(int at1, int at2, int ref1, int ref2, int co1, int co2, int fg, int& counter2) {
  double eps = 0.00001;
  xmatrix<int> match(3, 3, 1, 1);
  int counter = 0;
  xvector<int> at, ref, co, coF;
  at[1] = at1;
  at[2] = at2;
  co[1] = co1;
  co[2] = co2;
  ref[1] = ref1;
  ref[2] = ref2;

  for (int i = 1; i < 3; i++) {
    for (uint j = 1; j < 3; j++) {
      match[i][j] = 0;
    }
  }

  _atom tatom;
  for (uint i = 1; i < 3; i++) {
    tatom = SYM::ApplyAtom(a.atoms[at[i]], a.fgroup[fg], a);
    for (uint j = 1; j < 3; j++) {
      if (match[1][j] + match[2][j] < 1) {
        if (aurostd::modulus(BringInCell(a.atoms[at[j]].fpos - tatom.fpos)) < eps) {
          match[i][j] = 1;
          counter = counter + 1;
          j = 3;
        }
      }
    }
  }
  if (counter == 2) {
    for (int i = 1; i < 3; i++) {
      for (uint j = 1; j < 3; j++) {
        coF[j] = coF[j] + match[i][j] * co[i];
      }
    }
    for (int i = 1; i < 3; i++) {
      for (uint j = 1; j < 3; j++) {
        match[i][j] = 0;
      }
    }
    counter = 0;
    for (uint i = 1; i < 3; i++) {
      for (uint j = 1; j < 3; j++) {
        if (match[1][j] + match[2][j] < 1 && at[i] == at[j]) {
          if (aurostd::modulus(ortoDistorsions[ref[j]] - a.fgroup[fg].Uc * ortoDistorsions[coF[i]]) < eps) {
            match[i][j] = 1;
            counter = counter + 1;
            counter2 = counter2 * 1;
            j = 3;
          } else {
            if (aurostd::modulus(ortoDistorsions[ref[j]] + a.fgroup[fg].Uc * ortoDistorsions[coF[i]]) < eps) {
              match[i][j] = 1;
              counter = counter + 1;
              counter2 = counter2 * (-1);
              j = 3;
            }
          }
        }
      }
      if (match[i][1] + match[i][2] == 0) { return FALSE; }
    }

    if (counter == 2) {
      return TRUE;
    } else {
      return FALSE;
    }
  } else {
    return FALSE;
  }
}

///GetU/////////////////////////////////////////////////////////////////////
// This function is used to evaluate which symmetry operator link
//  pair (atom1,atom2) distorition ii with pair (atom1,atom2) distorition jj

int StrPairs::getU(const xstructure& b, int ipa, int atom1, int atom2, int ii, int jj) {
  // EPS: WARNING!
  double eps = 0.001;
  bool update = false;
  xmatrix<double> bc2f = b.c2f;
  xmatrix<double> bf2c = b.f2c;
  xvector<double> atom1v = b.atoms[atom1].fpos;
  xvector<double> atom2v = b.atoms[atom2].fpos;
  double bsym_eps = b.sym_eps;
  _atom tatom;
  _atom tatom2;
  uint k = ipairs[ipa].dist[ii][jj] / 6;
  uint l = ipairs[ipa].dist[ii][jj] % 6;
  uint fg = 0;
  for (fg = 0; fg < b.fgroup.size(); fg++) {
    tatom = SYM::ApplyAtom(b.atoms[atom1], b.fgroup[fg], b);
    tatom2 = SYM::ApplyAtom(b.atoms[atom2], b.fgroup[fg], b);
    //if(aurostd::modulus(BringInCell(b.atoms[atom1].fpos-tatom.fpos))<eps && aurostd::modulus(BringInCell(b.atoms[atom2].fpos-tatom2.fpos))<eps){
    if (SYM::AtomFPOSMatch(tatom.fpos, atom1v, bc2f, bf2c, update, bsym_eps) && SYM::AtomFPOSMatch(tatom2.fpos, atom2v, bc2f, bf2c, update, bsym_eps)) {
      xvector<double> vec1(3);
      xvector<double> vec2(3);
      vec1 = b.fgroup[fg].Uc * ortoDistorsions[k];
      vec2 = b.fgroup[fg].Uc * ortoDistorsions[l];
      if (atom1 == atom2) {
        if (aurostd::modulus(ortoDistorsions[ii] - vec1 + ortoDistorsions[jj] - vec2) < eps) {
          return fg;
        }
      } else {
        if (aurostd::modulus(ortoDistorsions[ii] - vec1) < eps && aurostd::modulus(ortoDistorsions[jj] - vec2) < eps) {
          return fg;
        }
      }
    }
  }
  _logger << error << "Problems looking for inequivalent pairs " << apl::endl;
  return fg;
}
///GetMinDist/////////////////////////////////////////////////////////////////////
// This function is used to evaluate minimun distance between two atoms
double StrPairs::getMinDist(const xmatrix<double>& lattice, xvector<double>& a, xvector<double>& b) {
  double dist = 100;
  for (int o = -2; o <= 2; o++)
    for (int p = -2; p <= 2; p++)
      for (int q = -2; q <= 2; q++)
        dist = aurostd::min(dist, aurostd::modulus(a - b + o * lattice(1) + p * lattice(2) + q * lattice(3)));
  return abs(dist);
}

///test/////////////////////////////////////////////////////////////////////
void StrPairs::test(const xmatrix<double>& lattice, xvector<double>& a, xvector<double>& b, double distM) {
  vector<double> list;
  for (int o = -2; o <= 2; o++) {
    for (int p = -2; p <= 2; p++) {
      for (int q = -2; q <= 2; q++) {
        if (abs(aurostd::modulus(a - b + o * lattice(1) + p * lattice(2) + q * lattice(3)) - distM) < 0.1) {
          list.push_back(aurostd::modulus(a - b + o * lattice(1) + p * lattice(2) + q * lattice(3)));
        }
      }
    }
  }
}

///GetMinDist/////////////////////////////////////////////////////////////////////
// This function is used to evaluate cutoof dependeng on shells

double StrPairs::getShellDist(const xstructure& a, int SHELL, double cutoff) {
  xvector<double> list(a.atoms.size(), 0);
  xvector<double> list2(a.atoms.size(), 0);
  double radious = 0.0;

  if (SHELL == 0) {
    return cutoff;
    //_logger << "Cutoff has been modified to " << cutoff <<  " Angstroms " << apl::endl;
  } else {
    for (uint i = 0; i < a.iatoms.size(); i++) {
      uint iat1 = a.iatoms[i][0];
      xvector<double> at1v = a.atoms[iat1].cpos;
      // Compute distances
      for (uint j = 0; j < a.atoms.size(); j++) {
        xvector<double> at2v = a.atoms[j].cpos;
        list[j] = abs(getMinDist(a.lattice, at1v, at2v));
      }
      // Sort
      list2 = aurostd::sort(list);
      // Check distances
      int countshell = 0;
      for (uint j = 1; j < a.atoms.size(); j++) {
        if (list2[j] > list2[j - 1] + 0.1) { countshell = countshell + 1; }
        if (countshell > SHELL) {
          radious = list2[j - 1] + 0.1;
          j = a.atoms.size();
        }
      }
      if (radious > cutoff) { cutoff = radious; }
    }
    _logger << "Cutoff has been modified to " << cutoff << " Angstroms " << apl::endl;
    return cutoff;
  }
}

///SYM_ROT/////////////////////////////////////////////////////////////////////

void StrPairs::SymUtils() {
  // Recalculating supercell dimensions
  xvector<double> sc_dim;
  for (uint ii = 1; ii < 4; ii++) {
    for (uint jj = 1; jj < 10000; jj++) {
      xvector<double> kk;
      kk = double(jj) * pcell.lattice(ii);
      if (aurostd::modulus(abs(kk - a.lattice(ii))) < _AFLOW_APL_EPS_) {
        sc_dim[ii] = jj;
        jj = 10000;
      }
    }
  }

  // Initialize id_equi
  id_equi.clear();
  id_equi.resize(pcell.fgroup.size(), vector<int>(a.atoms.size(), 0));
  for (uint isym2 = 0; isym2 < pcell.fgroup.size(); isym2++) {
    bool update = false;
    xvector<double> ftau_mod;
    ftau_mod[1] = pcell.fgroup[isym2].ftau[1] / sc_dim[1];
    ftau_mod[2] = pcell.fgroup[isym2].ftau[2] / sc_dim[2];
    ftau_mod[3] = pcell.fgroup[isym2].ftau[3] / sc_dim[3];
    for (uint iatoms = 0; iatoms < a.atoms.size(); iatoms++) {
      xvector<double> dummy;
      dummy = (pcell.fgroup[isym2].Uf * a.atoms[iatoms].fpos) + ftau_mod;
      bool check = TRUE;
      for (uint jatoms = 0; jatoms < a.atoms.size(); jatoms++) {
        if (SYM::AtomFPOSMatch(dummy, a.atoms[jatoms].fpos, a.c2f, a.f2c, update, a.sym_eps)) {
          id_equi[isym2][iatoms] = jatoms;
          jatoms = a.atoms.size();
          check = FALSE;
        }
      }
      if (check) {
        _logger << error << " SEVERE SYMMETRY ERROR. CHECK Sym Op mapping " << apl::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  //Initialize rot,rot2,nonzero
  rot.clear();
  rot2.clear();
  nonzero.clear();
  nonzero.resize(6, vector<vector<double> >(pcell.fgroup.size(), vector<double>(27, 0.)));
  rot.resize(6, vector<vector<vector<double> > >(pcell.fgroup.size(), vector<vector<double> >(27, vector<double>(27, 0.))));
  rot2.resize(6, vector<vector<vector<double> > >(pcell.fgroup.size(), vector<vector<double> >(27, vector<double>(27, 0.))));

  //Initialize permut
  permut.clear();
  //DX 5/2/18 - fix g++ version issue with vectors - START
  //DX 5/2/18 [OBSOLETE] permut.resize(6, vector<int>(3, 0));
  //DX 5/2/18 [OBSOLETE] permut[0] = {0, 1, 2};
  //DX 5/2/18 [OBSOLETE] permut[1] = {1, 0, 2};
  //DX 5/2/18 [OBSOLETE] permut[2] = {2, 1, 0};
  //DX 5/2/18 [OBSOLETE] permut[3] = {0, 2, 1};
  //DX 5/2/18 [OBSOLETE] permut[4] = {1, 2, 0};
  //DX 5/2/18 [OBSOLETE] permut[5] = {2, 0, 1};
  vector<int> tmp; 
  tmp.push_back(0); tmp.push_back(1); tmp.push_back(2);
  permut.push_back(tmp); tmp.clear();
  tmp.push_back(1); tmp.push_back(0); tmp.push_back(2);
  permut.push_back(tmp); tmp.clear();
  tmp.push_back(2); tmp.push_back(1); tmp.push_back(0);
  permut.push_back(tmp); tmp.clear();
  tmp.push_back(0); tmp.push_back(2); tmp.push_back(1);
  permut.push_back(tmp); tmp.clear();
  tmp.push_back(1); tmp.push_back(2); tmp.push_back(0);
  permut.push_back(tmp); tmp.clear();
  tmp.push_back(2); tmp.push_back(0); tmp.push_back(1);
  permut.push_back(tmp); tmp.clear();
  //DX 5/2/18 - fix g++ version issue with vectors - END

  //Rotation matrices for third derivatives and related quantities.
  for (int iperm = 0; iperm < 6; iperm++) {
    for (uint isym2 = 0; isym2 < pcell.fgroup.size(); isym2++) {
      //int isym=symor[isym2];
      for (int ip = 0; ip < 3; ip++) {
        for (int jp = 0; jp < 3; jp++) {
          for (int kp = 0; kp < 3; kp++) {
            int indexijkp = (ip * 3 + jp) * 3 + kp;
            vector<int> basis(3);
            for (int i = 0; i < 3; i++) {
              basis[0] = i;
              for (int j = 0; j < 3; j++) {
                basis[1] = j;
                for (int k = 0; k < 3; k++) {
                  basis[2] = k;
                  int indexijk = i * 9 + j * 3 + k;
                  int iper = basis[permut[iperm][0]];
                  int jper = basis[permut[iperm][1]];
                  int kper = basis[permut[iperm][2]];
                  rot[iperm][isym2][indexijkp][indexijk] = (orth(ip, iper, isym2) * orth(jp, jper, isym2) * orth(kp, kper, isym2));
                }
              }
            }
          }
        }
      }
    }
  }

  rot2 = rot;
  rot2 = rot;
  for (int iperm = 0; iperm < 6; iperm++) {
    for (uint isym = 0; isym < pcell.fgroup.size(); isym++) {
      for (int indexijkp = 0; indexijkp < 27; indexijkp++) {
        rot2[iperm][isym][indexijkp][indexijkp] -= 1.0;
        for (int indexijk = 0; indexijk < 27; indexijk++) {
          if (abs(rot2[iperm][isym][indexijkp][indexijk]) > 1E-12) {
            nonzero[iperm][isym][indexijkp] = 1.0;
          } else {
            rot2[iperm][isym][indexijkp][indexijk] = 0.0;
          }
        }
      }
    }
  }
}

void StrPairs::moveTri(vector<int>& triplet_sym) {
  vector<int> dummy(3);
  int atom1pc = sc2pcMap[triplet_sym[0]];
  int atom1sc = pc2scMap[atom1pc];
  dummy[0] = atom1sc;
  xvector<double> vec, newpos;
  vec = a.atoms[atom1sc].fpos - a.atoms[triplet_sym[0]].fpos;
  bool update = false;
  //Atom 2
  newpos = vec + a.atoms[triplet_sym[1]].fpos;
  for (uint at = 0; at < a.atoms.size(); at++) {
    if (SYM::AtomFPOSMatch(newpos, a.atoms[at].fpos, a.c2f, a.f2c, update, a.sym_eps)) {
      dummy[1] = at;
      at = a.atoms.size();
    }
  }
  //Atom 3
  newpos = vec + a.atoms[triplet_sym[2]].fpos;
  for (uint at = 0; at < a.atoms.size(); at++) {
    if (SYM::AtomFPOSMatch(newpos, a.atoms[at].fpos, a.c2f, a.f2c, update, a.sym_eps)) {
      dummy[2] = at;
      at = a.atoms.size();
    }
  }

  triplet_sym = dummy;
}

void StrPairs::gaussAAPL(int kk, vector<vector<double> > mat) {
  int i, j, k, irow;
  int row, col, ndependent, nindependent;
  double tmp;
  double _EPS_GAUSS_ = 1e-7;

  row = mat.size();
  col = mat[0].size();

  vector<int> dependent(col);
  vector<int> independent(col);
  vector<vector<double> > mat2(col, vector<double>(col, 0.0));

  irow = 0;
  ndependent = 0;
  nindependent = 0;
  for (k = 0; k < std::min(row, col); k++) {
    for (i = 0; i < row; i++) {
      if (abs(mat[i][k]) < _EPS_GAUSS_) {
        mat[i][k] = 0.;
      }
    }
    for (i = irow + 1; i < row; i++) {
      if (abs(mat[i][k]) - abs(mat[irow][k]) > _EPS_GAUSS_) {
        for (j = k; j < col; j++) {
          tmp = mat[irow][j];
          mat[irow][j] = mat[i][j];
          mat[i][j] = tmp;
        }
      }
    }
    if (abs(mat[irow][k]) > _EPS_GAUSS_) {
      dependent[ndependent] = k;
      ndependent += 1;
      for (j = col - 1; j > k; j--) {
        mat[irow][j] /= mat[irow][k];
      }
      mat[irow][k] = 1.;
      for (i = 0; i < row; i++) {
        if (i == irow) {
          continue;
        }
        for (j = col - 1; j > k; j--) {
          mat[i][j] -= mat[i][k] * mat[irow][j] / mat[irow][k];
        }
        mat[i][k] = 0.;
      }
      if (irow < row - 1) { irow += 1; }
    } else {
      independent[nindependent] = k;
      nindependent += 1;
    }
  }
  for (j = 0; j < nindependent; j++) {
    for (i = 0; i < ndependent; i++) {
      mat2[dependent[i]][j] = -mat[i][independent[j]];
    }
    mat2[independent[j]][j] = 1.;
  }

  itriplets[kk].rot_transform_aux.clear();
  itriplets[kk].rot_transform_aux.resize(27, vector<double>(27, 0.));
  itriplets[kk].rot_transform_aux = mat2;

  itriplets[kk].rot_independent.clear();
  itriplets[kk].rot_independent.resize(nindependent);
  for (j = 0; j < nindependent; j++) {
    itriplets[kk].rot_independent[j] = independent[j];
  }
}

double StrPairs::orth(int i, int j, int sym) {
  return pcell.fgroup[sym].Uc[i + 1][j + 1];
}

bool StrPairs::compareTri(vector<int> triplet, vector<int> triplet_sym) {
  for (uint i = 0; i < triplet.size(); i++) {
    if (triplet[i] != triplet_sym[i]) {
      return false;
    }
  }
  return true;
}
}

// ***************************************************************************
// *                                                                         *
// *              AFlow JOSE J. PLATA -   Duke University 2017               *
// *                                                                         *
// ***************************************************************************
