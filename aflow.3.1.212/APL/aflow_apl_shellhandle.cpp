#include "aflow_apl.h"

using namespace std;

#define _AFLOW_APL_EPS_SHELL_EQUIVALENCE_ 0.1
#define _AFLOW_APL_EPS_CPOS_EQUIVALENCE_ 0.1

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

ShellData::~ShellData() {
  for (uint i = 0; i < atoms.size(); i++)
    atoms[i].clear();
  atoms.clear();

  for (uint i = 0; i < ratoms.size(); i++)
    ratoms[i].clear();
  ratoms.clear();

  for (uint i = 0; i < index.size(); i++)
    index[i].clear();
  index.clear();
}

// ///////////////////////////////////////////////////////////////////////////

ShellData& ShellData::operator=(const ShellData& that) {
  if (this != &that) {
    occupation = that.occupation;
    occupationCapacity = that.occupationCapacity;
    isFull = that.isFull;
    radius = that.radius;
    stdevRadius = that.stdevRadius;

    atoms.clear();
    for (uint i = 0; i < that.atoms.size(); i++) {
      deque<_atom> row;
      for (uint j = 0; j < that.atoms[i].size(); j++)
        row.push_back(that.atoms[i][j]);
      atoms.push_back(row);
      row.clear();
    }

    ratoms.clear();
    for (uint i = 0; i < that.ratoms.size(); i++) {
      deque<_atom> row;
      for (uint j = 0; j < that.ratoms[i].size(); j++)
        row.push_back(that.ratoms[i][j]);
      ratoms.push_back(row);
      row.clear();
    }

    index.clear();
    for (uint i = 0; i < that.index.size(); i++)
      index.push_back(that.index[i]);
  }

  return *this;
}

// ///////////////////////////////////////////////////////////////////////////

ShellHandle::ShellHandle() {
}

ShellHandle::ShellHandle(const xstructure& xstr, int centralAtomID, int safeShellID) {
  init(xstr, centralAtomID, safeShellID);
}

// ///////////////////////////////////////////////////////////////////////////

ShellHandle::~ShellHandle() {
  this->clear();
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::clear() {
  _shells.clear();
  _centralAtomID = 0;
  _indexReductionConstant = 0.0;

  xvector<int> zero(3);
  zero.clear();
  _safeDimension = zero;

  _idSafeGeneratedShell = 0;
  _idSafeMappedShell = 0;
}

// ///////////////////////////////////////////////////////////////////////////

xvector<double> ShellHandle::getFPositionItsNearestImage(const xvector<double>& fposAtom,
                                                         const xvector<double>& fposCenter,
                                                         const xmatrix<double>& lattice) {
  double r2min = 1000000.0;
  double r2;
  xvector<double> rfmin(3), rf(3);
  xvector<double> rf0 = fposAtom - fposCenter;

  for (_AFLOW_APL_REGISTER_ int ii = 1; ii >= -1; ii--)
    for (_AFLOW_APL_REGISTER_ int jj = 1; jj >= -1; jj--)
      for (_AFLOW_APL_REGISTER_ int kk = 1; kk >= -1; kk--) {
        rf(1) = rf0(1) + (double)ii;
        rf(2) = rf0(2) + (double)jj;
        rf(3) = rf0(3) + (double)kk;
        r2 = aurostd::modulussquare(F2C(lattice, rf));
        if (r2 < r2min) {
          r2min = r2;
          rfmin = rf;
        }
      }

  return (rfmin);
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::init(const xstructure& xstr, int centralAtomID, int safeShellID) {
  _idSafeGeneratedShell = safeShellID;
  _initStructure = xstr;
  _initStructure_original = _initStructure;  //CO
  //_initStructure_atoms_original = _initStructure.atoms; //CO
  _centralAtomID = centralAtomID;

  // Generate points safely up desired shell
  int ii = 1;
  double radiusSafeShell = -1;

  while (true) {
    // Calculate shells
    calcShells(_initStructure, _centralAtomID, ii);

    // Stop if the radius does not change
    if ((((int)_shells.size() - 1) > safeShellID) &&
        (fabs(radiusSafeShell - _shells[safeShellID].radius) < _AFLOW_APL_EPS_)) {
      break;
    }
    radiusSafeShell = _shells[safeShellID].radius;

    // Increase if neccessary
    ii++;
    //cout << "++ -> " << ii << std::endl;
  }

  // Save current dimension as safe. At this setting, the order of all
  // shells is correct
  _safeDimension(1) = ii;
  _safeDimension(2) = ii;
  _safeDimension(3) = ii;
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::calcShells(const xstructure& xstr, int centralAtomID, int n) {
  // Clear old staff
  _shells.clear();

  // Count shell radius and its possible total occupations
  xvector<double> shift(3);
  for (_AFLOW_APL_REGISTER_ int ii = -n; ii <= n; ii++)
    for (_AFLOW_APL_REGISTER_ int jj = -n; jj <= n; jj++)
      for (_AFLOW_APL_REGISTER_ int kk = -n; kk <= n; kk++) {
        shift(1) = ii;
        shift(2) = jj;
        shift(3) = kk;
        shift = F2C(xstr.lattice, shift);

        for (uint i = 0; i < xstr.atoms.size(); i++) {
          double r = aurostd::modulus(xstr.atoms[i].cpos + shift - xstr.atoms[centralAtomID].cpos);

          uint j = 0;
          for (; j < _shells.size(); j++)
            if (fabs(_shells[j].radius - r) < _AFLOW_APL_EPS_SHELL_EQUIVALENCE_) break;
          if (j == _shells.size()) {
            ShellData sd;
            sd.radius = r;
            sd.occupation = 0;
            sd.occupationCapacity = 1;
            sd.isFull = false;
            deque<_atom> newAtomsList;
            newAtomsList.clear();
            sd.atoms.push_back(newAtomsList);
            _shells.push_back(sd);
          } else {
            _shells[j].occupationCapacity++;
          }
        }
      }

  // Order lists...
  for (uint i = 0; i < _shells.size() - 1; i++)
    for (uint j = i + 1; j < _shells.size(); j++) {
      if (_shells[j].radius < _shells[i].radius) {
        ShellData tsd = _shells[j];
        _shells[j] = _shells[i];
        _shells[i] = tsd;
      }
    }
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::splitBySymmetry() {
  //printReport(cout);

  xvector<double> zero(3);
  zero(1) = 0.0;
  zero(2) = 0.0;
  zero(3) = 0.0;

  // Create structure which fully fills all shells up to desired level
  _indexReductionConstant = 100000000.0;
  xmatrix<double> sclattice;
  while (true) {
    // The lattice of trial "superstructure"
    xmatrix<double> scale(3, 3);
    scale.clear();
    scale(1, 1) = (double)(2 * _safeDimension(1) + 1);
    scale(2, 2) = (double)(2 * _safeDimension(2) + 1);
    scale(3, 3) = (double)(2 * _safeDimension(3) + 1);
    sclattice = scale * _initStructure.lattice;

    // Clear old mapping
    for (uint i = 0; i < _shells.size(); i++) {
      _shells[i].occupation = 0;
      for (uint j = 0; j < _shells[i].atoms.size(); j++)
        _shells[i].atoms[j].clear();
    }

    // Fill shells...
    _atom atom;
    xvector<double> cshift(3);
    int id = 0;
    for (uint ia = 0; ia < _initStructure.iatoms.size(); ia++) {
      for (uint iia = 0; iia < _initStructure.iatoms[ia].size(); iia++) {
        // Replicate this atom by given mesh...
        //for(_AFLOW_APL_REGISTER_ int i = 0; i <= _safeDimension(1); i++)
        //for(_AFLOW_APL_REGISTER_ int j = 0; j <= _safeDimension(2); j++)
        //for(_AFLOW_APL_REGISTER_ int k = 0; k <= _safeDimension(3); k++)
        for (_AFLOW_APL_REGISTER_ int i = -_safeDimension(1); i <= _safeDimension(1); i++)
          for (_AFLOW_APL_REGISTER_ int j = -_safeDimension(2); j <= _safeDimension(2); j++)
            for (_AFLOW_APL_REGISTER_ int k = -_safeDimension(3); k <= _safeDimension(3); k++) {
              // Create position of new atoms...
              atom = _initStructure.atoms[_initStructure.iatoms[ia][iia]];
              cshift = (((double)i) * _initStructure.lattice(1) +
                        ((double)j) * _initStructure.lattice(2) +
                        ((double)k) * _initStructure.lattice(3));
              atom.cpos = atom.cpos + cshift - _initStructure.atoms[_centralAtomID].cpos;
              //atom.cpos = atom.cpos + cshift;
              atom.fpos = C2F(_initStructure.lattice, atom.cpos);
              //atom.fpos = C2F(sclattice,atom.cpos);
              atom.number = id++;

              // Which shell?
              //atom.fpos = getFPositionItsNearestImage(atom.fpos,zero,sclattice);
              //atom.cpos = F2C(sclattice,atom.fpos);

              // Get lowest index
              for (_AFLOW_APL_REGISTER_ int l = 1; l <= 3; l++) {
                if (fabs(atom.cpos(l) > _AFLOW_APL_EPS_CPOS_EQUIVALENCE_) &&
                    fabs(atom.cpos(l)) < _indexReductionConstant) {
                  _indexReductionConstant = atom.cpos(l);
                  //cout << _indexReductionConstant << std::endl;
                }
              }

              // Add to shell
              try {
                int shellID = getShell(aurostd::modulus(atom.cpos));
                addAtomToShell(shellID, atom, false);
              } catch (APLLogicError& e) {
                cout << "Atom cpos: ";
                printXVector(atom.cpos, false);
                cout << ";r = " << aurostd::modulus(atom.cpos) << std::endl;
                throw APLLogicError("apl::ShellHandle::splitBySymmetry(); No shell for this atom.");
              }
            }
      }
    }

    // If the last full shell is still lower than safe shell, increase
    // supercelll dimension and try it all again..., if satisfied break out.
    //printReport(cout);
    if (getLastFullShell() < _idSafeGeneratedShell) {
      _safeDimension(1)++;
      _safeDimension(2)++;
      _safeDimension(3)++;
      //cout << "++*-> " << _safeDimension(1) << std::endl;
      calcShells(_initStructure, _centralAtomID, _safeDimension(1));
    } else {
      break;
    }
  }

  // We will use the site point group of central atom, hence take the center
  // to this atom
  //_initStructure.ShifOriginToAtom(_centralAtomID); //CO
  //_initStructure.BringInCell(); //CO
  //COREY - START
  //we now do the same thing as with supercell, shiftorigin to central atom, then shift back original (not 0)
  //more robust
  center(_centralAtomID);  //corey
  //COREY - END

  // Calculate average shell radius
  for (uint ishell = 0; ishell < _shells.size(); ishell++) {
    const deque<_atom>& atomsAtSameShell = _shells[ishell].atoms[0];

    if (atomsAtSameShell.empty()) continue;

    // Average
    double sum = 0.0;
    for (uint i = 0; i < atomsAtSameShell.size(); i++)
      sum += aurostd::modulus(atomsAtSameShell[i].cpos);
    _shells[ishell].radius = sum / (double)atomsAtSameShell.size();

    // Standard deviation
    sum = 0.0;
    for (uint i = 0; i < atomsAtSameShell.size(); i++)
      sum += pow(aurostd::modulus(atomsAtSameShell[i].cpos) - _shells[ishell].radius, 2.0);
    _shells[ishell].stdevRadius = sqrt(sum / (double)atomsAtSameShell.size());
  }

  // Group atoms by the symmetry equivalence, hence some shells can split...
  deque<_atom> storeRow;
  vector<deque<_atom> > store;
  const vector<vector<_sym_op> >& agroup = _initStructure.agroup;

  for (uint ishell = 0; ishell < _shells.size(); ishell++) {
    const deque<_atom>& atomsAtSameShell = _shells[ishell].atoms[0];

    if (atomsAtSameShell.empty()) continue;

    for (uint i = 0; i < atomsAtSameShell.size(); i++) {
      uint k = 0;
      for (; k < store.size(); k++) {
        uint l = 0;
        for (; l < store[k].size(); l++)
          if (atomsAtSameShell[i].number == store[k][l].number) break;
        if (l != store[k].size()) break;
      }
      if (k != store.size()) continue;

      for (uint j = k; j < atomsAtSameShell.size(); j++) {
        // Is it already stored?
        uint k = 0;
        for (; k < storeRow.size(); k++)
          if (atomsAtSameShell[j].number == storeRow[k].number) break;
        if (k != storeRow.size()) continue;

        // Rotate atom and check if it coincides with the first atom of this group
        uint symOpID = 0;
        for (; symOpID < agroup[_centralAtomID].size(); symOpID++) {
          xvector<double> rotcpos = agroup[_centralAtomID][symOpID].Uc * atomsAtSameShell[j].cpos;
          xvector<double> rotfpos = C2F(sclattice, rotcpos);
          xvector<double> xmin = getFPositionItsNearestImage(rotfpos, zero, sclattice);
          rotcpos = F2C(sclattice, xmin);
          if ((fabs(fabs(rotcpos(1)) - fabs(atomsAtSameShell[i].cpos(1))) < _AFLOW_APL_EPS_CPOS_EQUIVALENCE_) &&
              (fabs(fabs(rotcpos(2)) - fabs(atomsAtSameShell[i].cpos(2))) < _AFLOW_APL_EPS_CPOS_EQUIVALENCE_) &&
              (fabs(fabs(rotcpos(3)) - fabs(atomsAtSameShell[i].cpos(3))) < _AFLOW_APL_EPS_CPOS_EQUIVALENCE_))
            break;
        }
        if (symOpID != agroup[_centralAtomID].size())
          storeRow.push_back(atomsAtSameShell[j]);
      }

      store.push_back(storeRow);
      storeRow.clear();
    }

    /*
        cout << "SHELL ID: " << ishell << std::endl;
        cout << "groups: " << store.size() << " | ";
        for(uint m = 0; m < store.size(); m++)
            cout << store[m].size() << " ";
        cout << std::endl;
        cout << "radius = " << _shells[ishell].radius << " +/- " << _shells[ishell].stdevRadius << std::endl;
        for(uint m = 0; m < store.size(); m++)
        {
            for(uint n = 0; n < store[m].size(); n++)
            {
                cout << "  - ATOM ID    : " << setw(3) << setfill('0') << store[m][n].number << " " << setfill(' ');
                printXVector(store[m][n].cpos,false);
                double r = aurostd::modulus(store[m][n].cpos);
                cout << " (" << r << ", dev = " << ( r - _shells[ishell].radius ) << ")" << std::endl;
            }
            cout << "-" << std::endl;
        }
        */

    // Clear old mapping
    for (uint i = 0; i < _shells[ishell].ratoms.size(); i++) {
      _shells[ishell].ratoms[i].clear();
      _shells[ishell].atoms[i].clear();
    }
    _shells[ishell].atoms.clear();
    _shells[ishell].ratoms.clear();
    _shells[ishell].occupation = 0;

    // Copy this configuration as reference
    for (uint i = 0; i < store.size(); i++) {
      deque<_atom> row;
      for (uint j = 0; j < store[i].size(); j++)
        row.push_back(store[i][j]);
      _shells[ishell].ratoms.push_back(row);
      // Add empty row to atoms
      row.clear();
      _shells[ishell].atoms.push_back(row);
    }
    store.clear();

    // Create index for this shell
    _shells[ishell].index.clear();
    for (uint i = 0; i < _shells[ishell].ratoms.size(); i++) {
      xvector<int> pint(3);
      if (!_shells[ishell].ratoms[i].empty()) {
        xvector<double> p = _shells[ishell].ratoms[i][0].cpos;
        p(1) = fabs(p(1));
        p(2) = fabs(p(2));
        p(3) = fabs(p(3));
        for (int j = 1; j < 3; j++)
          for (int k = j + 1; k <= 3; k++) {
            if (p(k) > p(j)) {
              double temp = p[j];
              p[j] = p[k];
              p[k] = temp;
            }
          }
        pint(1) = (int)(0.5 + (p(1) / _indexReductionConstant));
        pint(2) = (int)(0.5 + (p(2) / _indexReductionConstant));
        pint(3) = (int)(0.5 + (p(3) / _indexReductionConstant));
      } else {
        pint.clear();
      }
      _shells[ishell].index.push_back(pint);
    }
  }

  // Take the center of init structure into atom 0 (by default)
  //_initStructure.ShifOriginToAtom(0); //CO
  //_initStructure.BringInCell(); //CO
  //COREY - START
  //we now do the same thing as with supercell, shiftorigin to central atom, then shift back original (not 0)
  //more robust
  center_original();  //corey
                      //COREY - END
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::removeSplitBySymmetry() {
  for (uint ishell = 0; ishell < _shells.size(); ishell++) {
    // Clear stdev
    _shells[ishell].stdevRadius = 0.0;

    // Clear useless arrays
    for (uint j = 0; j < _shells[ishell].index.size(); j++)
      _shells[ishell].index[j].clear();
    _shells[ishell].index.clear();

    for (uint j = 0; j < _shells[ishell].ratoms.size(); j++)
      _shells[ishell].ratoms[j].clear();
    _shells[ishell].ratoms.clear();

    // Split all atoms into the first subshell (the zero index)
    for (uint j = 1; j < _shells[ishell].atoms.size(); j++) {
      for (uint k = 0; k < _shells[ishell].atoms[j].size(); k++)
        _shells[ishell].atoms[0].push_back(_shells[ishell].atoms[j][k]);
      _shells[ishell].atoms[j].clear();
    }

    // Remove rest subshells
    _shells[ishell].atoms.erase(_shells[ishell].atoms.begin() + 1,
                                _shells[ishell].atoms.end());
  }
}

// ///////////////////////////////////////////////////////////////////////////

double ShellHandle::getShellRadius(int nn) {
  if (nn >= (int)_shells.size())
    throw APLRuntimeError("ShellHandle::getShellRadius(); Index out of range.");
  return (_shells[nn].radius);
}

// ///////////////////////////////////////////////////////////////////////////

int ShellHandle::getShell(double r) {
  uint j = 0;
  for (; j < _shells.size(); j++) {
    //  3*sigma garanties there is 99.7300204% points in,..
    double sigmaMul = 3.0;
    bool found = false;

    while (sigmaMul > _AFLOW_APL_EPS_) {
      double eps = (_shells[j].stdevRadius > _AFLOW_APL_EPS_) ? sigmaMul * _shells[j].stdevRadius : _AFLOW_APL_EPS_SHELL_EQUIVALENCE_;
      //cout << j << " " <<  fabs( _shells[j].radius - r ) << " +/- " << _shells[j].stdevRadius << " " << eps << std::endl;
      //cout <<  ( fabs( _shells[j].radius - r ) < eps ) << std::endl;
      //cout <<  ( j < _shells.size()-1 ? fabs( _shells[j+1].radius - r ) > eps : true ) << std::endl;
      if ((fabs(_shells[j].radius - r) < eps) &&
          (j < _shells.size() - 1 ? fabs(_shells[j + 1].radius - r) > eps : true)) {
        found = true;
        break;
      }
      if (_shells[j].stdevRadius < _AFLOW_APL_EPS_)
        break;
      sigmaMul -= 1.0;
    }

    // If this is satisfied, we founded the shell number and break...
    if (found) break;
  }

  if (j == _shells.size())
    throw APLLogicError("apl::ShellHandle::getShell(); Problem to find shell for this radius.");

  return ((int)j);
}

// ///////////////////////////////////////////////////////////////////////////

double ShellHandle::getSafeShellRadius() {
  return (_shells[_idSafeMappedShell].radius);
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::setSafeShell(int i) {
  _idSafeMappedShell = i;
}

// ///////////////////////////////////////////////////////////////////////////

int ShellHandle::getSafeShell() {
  return _idSafeMappedShell;
}

// ///////////////////////////////////////////////////////////////////////////

int ShellHandle::getLastOccupiedShell() {
  int maxShellID = _shells.size() - 1;
  for (; maxShellID >= 0; maxShellID--)
    if (_shells[maxShellID].occupation != 0) break;

  return maxShellID;
}

// ///////////////////////////////////////////////////////////////////////////

// All shells (and subshells) below the regular shell are occupied at least by 1 atom
int ShellHandle::getLastRegularShell() {
  uint maxShellID = 0;
  for (; maxShellID < _shells.size(); maxShellID++) {
    if (_shells[maxShellID].atoms.size() > 1) {
      uint i = 0;
      for (; i < _shells[maxShellID].atoms.size(); i++)
        if (_shells[maxShellID].atoms[i].empty()) break;
      if (i != _shells[maxShellID].atoms.size()) break;
    } else if (_shells[maxShellID].occupation == 0)
      break;
  }

  return --maxShellID;
}

// ///////////////////////////////////////////////////////////////////////////

// All shells below the regular shell are fully occupied
int ShellHandle::getLastFullShell() {
  uint maxShellID = 0;
  for (; maxShellID < _shells.size(); maxShellID++)
    if (!_shells[maxShellID].isFull) break;

  return --maxShellID;
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::addAtomToShell(int nn, const _atom& atom, bool useSplittedShells) {
  if (nn >= (int)_shells.size())
    throw APLRuntimeError("apl::ShellHandle::addAtomToShell(); Index out of range.");

  // If there is more shells of the same radius try to add atom with the
  // help of reference atoms
  if (useSplittedShells && !_shells[nn].ratoms.empty() && _shells[nn].atoms.size() > 1) {
    uint i = 0;
    for (; i < _shells[nn].ratoms.size(); i++) {
      uint j = 0;
      for (; j < _shells[nn].ratoms[i].size(); j++)
        if (aurostd::modulus(atom.cpos - _shells[nn].ratoms[i][j].cpos) < _AFLOW_APL_EPS_CPOS_EQUIVALENCE_)
          break;
      if (j != _shells[nn].ratoms[i].size()) break;
    }
    if (i != _shells[nn].ratoms.size()) {
      if (_shells[nn].ratoms.size() != _shells[nn].atoms.size())
        throw APLLogicError("apl::ShellHandle::addAtomToShell(); Problem to add atom into symmetry splitted shell list. Arrays differ.");
      else
        _shells[nn].atoms[i].push_back(atom);
    } else {
      cout << "SHELL: " << nn << std::endl;
      cout << "add: ";
      printXVector(atom.cpos);
      for (i = 0; i < _shells[nn].ratoms.size(); i++) {
        for (uint j = 0; j < _shells[nn].ratoms[i].size(); j++) {
          cout << " ->  ";
          printXVector(_shells[nn].ratoms[i][j].cpos);
        }
      }
      throw APLLogicError("apl::ShellHandle::addAtomToShell(); None reference atom found.");
    }
  } else {
    _shells[nn].atoms.back().push_back(atom);
  }
  _shells[nn].occupation++;

  // Set isFull flag if ...
  if (_shells[nn].occupation == _shells[nn].occupationCapacity)
    _shells[nn].isFull = true;

  // Setup safe shell of the mapping
  for (uint j = 1; j < _shells.size(); j++) {
    if (!_shells[j].isFull) {
      _idSafeMappedShell = j - 1;
      break;
    }
  }
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::mapStructure(const xstructure& xstr, int centralAtomID, bool useSplittedShells) {
  // Clear old mapping
  for (uint i = 0; i < _shells.size(); i++) {
    _shells[i].occupation = 0;
    for (uint j = 0; j < _shells[i].atoms.size(); j++)
      _shells[i].atoms[j].clear();
  }

  //cout << "--map-begin--"<<std::endl;
  // Create new one
  for (uint i = 0; i < xstr.atoms.size(); i++) {
    xvector<double> xmin = getFPositionItsNearestImage(xstr.atoms[i].fpos,
                                                       xstr.atoms[centralAtomID].fpos,
                                                       xstr.lattice);
    xvector<double> cpos = F2C(xstr.lattice, xmin);
    //cout << i << " ";
    //printXVector(xstr.atoms[i].fpos,false);
    //cout << " ";
    //printXVector(xstr.atoms[i].cpos,false);
    //cout << " -> ";
    //printXVector(xmin,false);
    //cout << " ";
    //printXVector(cpos);
    try {
      int shellID = getShell(aurostd::modulus(cpos));
      _atom atom = xstr.atoms[i];
      atom.number = i;
      atom.fpos = xmin;
      atom.cpos = cpos;
      addAtomToShell(shellID, atom, useSplittedShells);
    } catch (APLLogicError& e) {
      cout << "Atom cpos: ";
      printXVector(cpos, false);
      cout << "; r = " << aurostd::modulus(cpos) << std::endl;
      throw APLLogicError("apl::ShellHandle::mapStructure(); There is an atom which does not correspond to any shell.");
    }
  }
  //cout << "--map-end--"<<std::endl;
}

// ///////////////////////////////////////////////////////////////////////////

int ShellHandle::getNumberOfShells() {
  return ((int)_shells.size());
}

// ///////////////////////////////////////////////////////////////////////////

int ShellHandle::getNumberOfSubshells(int ishell) {
  if (ishell > (int)_shells.size() - 1)
    throw APLRuntimeError("apl::ShellHandle::getNumberOfSubshells; The shell index is out of range.");
  return ((int)_shells[ishell].ratoms.size());
}

// ///////////////////////////////////////////////////////////////////////////

std::deque<_atom> ShellHandle::getAtomsAtSameShell(int ishell, int isubshell) {
  if (ishell > (int)_shells.size() - 1)
    throw APLRuntimeError("apl::ShellHandle::getAtomsAtSameShell(); The shell index is out of range.");
  if (isubshell > (int)_shells[ishell].atoms.size() - 1)
    throw APLRuntimeError("apl::ShellHandle::getAtomsAtSameShell(); The subshell index is out of range.");

  return (_shells[ishell].atoms[isubshell]);
}

// ///////////////////////////////////////////////////////////////////////////

const std::deque<_atom>& ShellHandle::getReferenceAtomsAtSameShell(int ishell, int isubshell) {
  if (ishell > (int)_shells.size() - 1)
    throw APLRuntimeError("apl::ShellHandle::getAtomsAtSameShell(); The shell index is out of range.");
  if (isubshell > (int)_shells[ishell].ratoms.size() - 1)
    throw APLRuntimeError("apl::ShellHandle::getAtomsAtSameShell(); The subshell index is out of range.");

  return (_shells[ishell].ratoms[isubshell]);
}

// ///////////////////////////////////////////////////////////////////////////

void ShellHandle::printReport(std::ostream& os) {
  int regShell = getLastRegularShell();
  int lastShell = getLastOccupiedShell();

  if (lastShell == -1) return;

  os << std::endl;
  os << "- To obtain the safe order of all shells up to level " << _idSafeGeneratedShell
     << ", the " << _safeDimension(1) << "x" << _safeDimension(2)
     << "x" << _safeDimension(3) << " supercell is needed." << std::endl;

  os << "- The inserted supercell is safely mapped up to shell " << _idSafeMappedShell << "." << std::endl;
  os << "- The regular shell is " << regShell << ", the last shell is " << lastShell << "." << std::endl;
  os << std::endl;
  os << "- The shell table with the center in the atom no. " << _shells[0].atoms[0][0].number << " is:" << std::endl;
  os << std::endl;
  string tab = "  ";
  os << tab << " SHELL  INDEX  OCC  MAX    Radius (Angs.)               Presented atoms               " << std::endl;
  os << tab << "--------------------------------------------------------------------------------------" << std::endl;
  for (int i = 0; i <= regShell; i++)
  //for(int i = 0; i <= ( (int)_shells.size() > _idSafeGeneratedShell ? _idSafeGeneratedShell:(int)_shells.size() ); i++)
  {
    for (uint j = 0; j < _shells[i].atoms.size(); j++) {
      os << tab << "  ";
      os << setw(2) << setfill('0') << i;
      if (_shells[i].atoms.size() > 1)
        os << (char)('a' + j) << "   ";
      else
        os << "    ";

      if (!_shells[i].index.empty()) {
        os << " " << _shells[i].index[j](1) << _shells[i].index[j](2)
           << _shells[i].index[j](3) << "   ";
      } else
        os << "       ";
      if (_shells[i].ratoms.empty())
        os << setfill('0') << setw(3) << _shells[i].occupation << "  "
           << setfill('0') << setw(3) << _shells[i].occupationCapacity;
      else
        os << setfill('0') << setw(3) << _shells[i].atoms[j].size() << "  "
           << setfill('0') << setw(3) << _shells[i].ratoms[j].size();
      if (i == _idSafeMappedShell)
        os << " S";
      else if (i == regShell)
        os << " R";
      else if (i == lastShell && lastShell != regShell)
        os << " L";
      else
        os << "  ";
      os << "     ";
      os << fixed << setfill(' ') << setw(6) << setprecision(3) << _shells[i].radius << "         ";
      for (uint k = 1; k <= _shells[i].atoms[j].size(); k++) {
        os << setfill('0') << setw(3) << _shells[i].atoms[j][k - 1].number << " ";
        if (k % 10 == 0 && k != _shells[i].atoms[j].size())
          os << std::endl
             << tab << "                                             ";
      }
      os << std::endl;
    }
  }
  os << tab << "--------------------------------------------------------------------------------------" << std::endl;
  os << std::endl
     << setfill(' ');
}

//COREY - START
void ShellHandle::center(int i) {
  //_initStructure.ShifOriginToAtom(i);
  //_initStructure.BringInCell();
  xvector<double> origin(3), frigin(3);
  origin = _initStructure_original.atoms[i].cpos;
  frigin = _initStructure_original.atoms[i].fpos;
  _initStructure.origin = origin;
  for (uint i = 0; i < _initStructure.atoms.size(); i++) {
    _initStructure.atoms[i].cpos = _initStructure_original.atoms[i].cpos - origin;
    _initStructure.atoms[i].fpos = _initStructure_original.atoms[i].fpos - frigin;
    _initStructure.atoms[i] = BringInCell(_initStructure.atoms[i], _initStructure.lattice);
  }
}

void ShellHandle::center_original(void) {
  //for(uint i=0;i<_initStructure.atoms.size();i++){
  //  _initStructure.atoms[i].fpos=_initStructure_original.atoms[i].fpos;
  //  _initStructure.atoms[i].cpos=_initStructure_original.atoms[i].cpos;
  //}
  _initStructure.origin = _initStructure_original.origin;
  for (uint i = 0; i < _initStructure.atoms.size(); i++) {
    _initStructure.atoms[i].fpos = _initStructure_original.atoms[i].fpos;
    _initStructure.atoms[i].cpos = _initStructure_original.atoms[i].cpos;
  }
}
//COREY - END

// ///////////////////////////////////////////////////////////////////////////

}  // namespace apl
