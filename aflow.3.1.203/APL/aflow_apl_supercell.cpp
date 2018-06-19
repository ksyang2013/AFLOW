// [OBSOLETE] #include <limits>

#include "aflow_apl.h"

#define _SYM_AFLOW_APL_EPS_ 0.05

//CO - START
#define ERROR_VERBOSE false
#define GETFULLSYMBASIS false  //TRUE is NOT well tested with where() and wherefrom(), also slower
#define JAHNATEK_ORIGINAL false
#define CENTER_PRIM false
#define MAP_VERBOSE false
#define COPY_XSTRUCTURE_FULL false
//CO - END

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

#if !JAHNATEK_ORIGINAL
Supercell::Supercell(const xstructure& xstr, Logger& l) : _logger(l) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  // Copy
  _inStructure = xstr;
  _inStructure.ReScale(1.0);
  _inStructure.ShifOriginToAtom(0);
  _inStructure.BringInCell();  //clean up

  if(LDEBUG){
    cerr << "APL::Supercell::Supercell: this is the structure to be analyzed" << std::endl;
    cerr << _inStructure << std::endl;
  }

  //COREY, DO NOT MODIFY THE STRUCTURE BELOW HERE, THIS INCLUDES RESCALE(), BRINGINCELL(), SHIFORIGINATOM(), etc.
  _logger << "Estimating the symmetry of structure and calculating the input structure. Please be patient." << apl::endl; //primitive cell." << apl::endl; //CO 180216 - we do NOT primitivize unless requested via [VASP_FORCE_OPTION]CONVERT_UNIT_CELL
  calculateWholeSymmetry(_inStructure);
  _pcStructure = calculatePrimitiveStructure(); //calculate and store primitive cell

  //store some mutual properties
  //need this for getInputStructureLight()
  //#if CENTER_PRIM
  _inStructure_original = _inStructure;
  LightCopy(_inStructure, _inStructure_light);
  //_inStructure_atoms_original = _inStructure.atoms;
  //#endif
  _skew = SYM::isLatticeSkewed(_inStructure.lattice, _inStructure.dist_nn_min, _inStructure.sym_eps);
  _sym_eps = _inStructure.sym_eps;
  clear();
}
#else
Supercell::Supercell(const xstructure& xstr, Logger& l) : _logger(l) {
  // Copy
  _inStructure = xstr;
  _inStructure.ReScale(1.0);
  _inStructure.ShifOriginToAtom(0);
  _inStructure.BringInCell();

  //
  _logger << "Reducing the input structure into a standard primitive form." << apl::endl;
  //_inStructure = GetPrimitive(xstr);
  _inStructure.Standard_Primitive_UnitCellForm();
  _inStructure.ReScale(1.0);
  _inStructure.ShifOriginToAtom(0);
  _inStructure.BringInCell();

  // Calculate symmetry
  _logger << "Estimating the symmetry of primitive cell." << apl::endl;
  calculateWholeSymmetry(_inStructure);

  //
  clear();
}
#endif

// ///////////////////////////////////////////////////////////////////////////

Supercell::Supercell(const Supercell& that) : _logger(that._logger) {
  *this = that;
}

// ///////////////////////////////////////////////////////////////////////////

Supercell& Supercell::operator=(const Supercell& that) {
  if (this != &that) {
    _logger = that._logger;

    _inStructure = that._inStructure;
    _inStructure_original = that._inStructure_original;  //CO
    _inStructure_light = that._inStructure_light;        //CO
    _pcStructure = that._pcStructure;                    //CO
    _scStructure = that._scStructure;                    //CO
    _scStructure_original = that._scStructure_original;  //CO
    _scStructure_light = that._scStructure_light;        //CO
    _pc2scMap.clear();
    _pc2scMap = that._pc2scMap;
    _sc2pcMap.clear();
    _sc2pcMap = that._sc2pcMap;
    //CO - START
    _skew = that._skew;
    _derivative_structure = that._derivative_structure;
    _sym_eps = that._sym_eps;
    //CO - END
    _isShellRestricted = that._isShellRestricted;
    _maxShellRadius.clear();
    _maxShellRadius = that._maxShellRadius;
    _isConstructed = that._isConstructed;  //CO
  }

  return *this;
}

// ///////////////////////////////////////////////////////////////////////////

Supercell::~Supercell() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////

void Supercell::clear() {
  _scStructure.info = "not constructed";
  _isConstructed = FALSE;
  _isShellRestricted = FALSE;
  _pc2scMap.clear();
  _sc2pcMap.clear();
  //_skew = FALSE;  //CO, same for _sc and _pc (DO NOT RESET)
  //_derivative_structure = FALSE;  //CO, same for _sc and _pc (DO NOT RESET)
  //_sym_eps = AUROSTD_NAN; //CO, same for _sc and _pc (DO NOT RESET)
  _maxShellRadius.clear();
  _maxShellID = -1;
}

// ///////////////////////////////////////////////////////////////////////////

//CO - START
void Supercell::LightCopy(const xstructure& a, xstructure& b) {
#if COPY_XSTRUCTURE_FULL
  b = a;
  b.pgroup.clear();
  b.pgroupk.clear();
  b.pgroupk_xtal.clear();
  b.fgroup.clear();
  for (uint i = 0; i < b.agroup.size(); i++) {
    b.agroup[i].clear();
  }
  b.agroup.clear();
#else
  b.Clear();
  stringstream POSCAR;
  POSCAR.str("");
  if(0){cerr << a << std::endl;}
  POSCAR << a;
  POSCAR >> b;
  //enable inequivalent flag to work
  for (uint i = 0; i < b.atoms.size(); i++) {
    b.atoms[i].equivalent = a.atoms[i].equivalent;
    b.atoms[i].is_inequivalent = a.atoms[i].is_inequivalent;
    b.atoms[i].num_equivalents = a.atoms[i].num_equivalents;
  }
  //[OBSOLETE] pseudo-potential stuff
  //for(uint i=0;i<b.species.size();i++){
  //  b.species[i]=a.species[i];
  //  b.species_pp[i]=a.species_pp[i]; //VERY IMPORTANT
  //}
  //enable inequivalent flag to work
  b.write_inequivalent_flag = a.write_inequivalent_flag;
  b.info = a.info;
  if(0){cerr << b << std::endl;}
#endif
}
//CO - END

// ///////////////////////////////////////////////////////////////////////////

#if !JAHNATEK_ORIGINAL
void Supercell::calculateWholeSymmetry(xstructure& xstr) {
  _aflags af;
  af.Directory = "./";
  af.QUIET = FALSE;

  //CO 170804 - we want to append symmetry stuff to ofstream
  _kflags kflags;
  kflags.KBIN_SYMMETRY_PGROUP_WRITE=TRUE;
  kflags.KBIN_SYMMETRY_FGROUP_WRITE=TRUE;
  kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=TRUE;
  kflags.KBIN_SYMMETRY_SGROUP_WRITE=TRUE;
  kflags.KBIN_SYMMETRY_IATOMS_WRITE=TRUE;
  kflags.KBIN_SYMMETRY_AGROUP_WRITE=TRUE;
  
  //DX CAN REMOVE string options = "";

  xstr.LatticeReduction_avoid = TRUE;
  xstr.sgroup_radius = 8.0;

  //CO 170804 - we want to append symmetry stuff to ofstream
  //if (!pflow::CalculateFullSymmetry(af, xstr)) {
  if (!pflow::PerformFullSymmetry(xstr,_logger.getOutputStream(),af,kflags,true,cout)) {
    throw APLRuntimeError("apl::Supercell::calculateWholeSymmetry(): Symmetry routine failed.");
  }
}
#else
void Supercell::calculateWholeSymmetry(xstructure& xstr) {
  ofstream fileDevNull("/dev/null");

  if (!fileDevNull.is_open()) {
    throw APLRuntimeError("apl::Supercell::calculateWholeSymmetry(): Cannot open output stream /dev/null.");
  }

  _aflags af;
  af.Directory = "./";
  // af.hostname = "";
  af.QUIET = TRUE;

  xstr.LatticeReduction_avoid = TRUE;
  xstr.sgroup_radius = 8.0;

  //if( xstr.pgroup_calculated == FALSE )
  SYM::CalculatePointGroup(fileDevNull, xstr, af, FALSE, FALSE, cout);
  //if( xstr.fgroup_calculated == FALSE )
  SYM::CalculateFactorGroup(fileDevNull, xstr, af, FALSE, FALSE, cout);
  //if( xstr.sgroup_calculated == FALSE )
  SYM::CalculateSpaceGroup(fileDevNull, xstr, af, FALSE, FALSE, cout);
  //if( xstr.iatoms_calculated == FALSE )
  SYM::CalculateInequivalentAtoms(fileDevNull, xstr, af, FALSE, FALSE, cout);
  //if( xstr.agroup_calculated == FALSE )
  SYM::CalculateSitePointGroup(fileDevNull, xstr, af, FALSE, FALSE, cout);

  fileDevNull.clear();
  fileDevNull.close();
}
#endif

// ///////////////////////////////////////////////////////////////////////////

void Supercell::reset() {
  _scStructure = _inStructure;

  _scStructure.atoms.clear();

  for (uint i = 0; i < _scStructure.num_each_type.size(); i++)
    _scStructure.num_each_type[i] = 0;

  for (uint i = 0; i < _scStructure.iatoms.size(); i++)
    _scStructure.iatoms[i].clear();

  _scStructure.agroup.clear();
  _scStructure.fgroup.clear();

  _pc2scMap.clear();
  _sc2pcMap.clear();

  _isConstructed = FALSE;
}

// ///////////////////////////////////////////////////////////////////////////

#if !JAHNATEK_ORIGINAL
void Supercell::build(int nx, int ny, int nz, bool VERBOSE) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //BEGIN JJPR
  scell(1) = nx;
  scell(2) = ny;
  scell(3) = nz;
  _derivative_structure = !(nx == ny && ny == nz);
  bool get_full_sym = GETFULLSYMBASIS;  //( GETFULLSYMBASIS || _derivative_structure ); //CO
  //END JJPR

  // Print info
  if (VERBOSE) {
    _logger << "The supercell is going to build as " << nx << " x " << ny << " x " << nz
            << " (" << (uint)(nx * ny * nz * _inStructure.atoms.size()) << " atoms)." << apl::endl;
  }

  if (VERBOSE && _derivative_structure) {
    _logger << "Derivative structure detected, be patient as we calculate symmetry of the supercell." << apl::endl;
  }
  // Create lattice of the supercell
  xmatrix<double> scale(3, 3);
  scale.clear();
  scale(1, 1) = nx;
  scale(2, 2) = ny;
  scale(3, 3) = nz;

  // Get supercell
  //_scStructure = GetSuperCell(_inStructure, scale, _sc2pcMap, _pc2scMap, TRUE, _derivative_structure);  //now gets symmetries too! no need for full_basis (just a check)
  _scStructure = GetSuperCell(_inStructure, scale, _sc2pcMap, _pc2scMap, TRUE, get_full_sym);  //now gets symmetries too! no need for full_basis (just a check)

  //  cerr << _scStructure << std::endl;
  //  for(uint i=0;i<_scStructure.agroup.size();i++){
  //    cerr << "AGROUP " << i << ": " << _scStructure.agroup[i].size() << std::endl;
  //for(uint j=0;j<a.agroup[i].size();j++){
  //}
  //  }
  /*for(uint i=0;i<_scStructure.agroup.size();i++){
    cerr << "SITE " << i << std::endl;
    for(uint j=0;j<_scStructure.agroup[i].size();j++){
      cerr << "OPERATION " << j << std::endl;
      cerr << _scStructure.agroup[i][j] << std::endl;
    //for(uint j=0;j<a.agroup[i].size();j++){
    }
  }*/
  //  exit(0);

  //for(uint i=0;i<_scStructure.agroup[0].size();i++){
  //  cerr << _scStructure.agroup[0][i] << std::endl;
  //}
  //exit(0);  //COREY REMOVE
  //_scStructure = GetSuperCell(_inStructure, scale, _sc2pcMap, _pc2scMap, TRUE, GETFULLSYMBASIS);  //now gets symmetries too! no need for full_basis (just a check)
  //_scStructure.ReScale(1.0); no longer needed, pc is already rescaled, also, try not to change structure much after calculating symmetry

  //it is TRUE that the symmetry of the supercell != symmetry of primitive cell if
  //!(nx==ny==nz), we now have a derivative structure
  //derivative structures have REDUCED symmetry, but instead of recalculating,
  //we will IGNORE failed mappings ONLY if we have a derivative structure
  //if(!( nx == ny && ny == nz)) {
  //  //SUPER slow
  //  _logger << "Supercell is not symmetric, hence we need to recalculate the whole symmetry. (very slow)" << apl::endl;
  //  calculateWholeSymmetry(_scStructure);
  //}

  // Setup output flags
  _scStructure.write_inequivalent_flag = TRUE;

  // Set the information about this construction
  _scStructure.info = "Supercell " + stringify(nx) + "x" + stringify(ny) + "x" + stringify(nz);

  // OK.
  if (VERBOSE) _logger << "Supercell successfully created." << apl::endl;
  _isConstructed = TRUE;

  _scStructure_original = _scStructure;  //COPY EVERYTHING ONCE, will be very slow
  LightCopy(_scStructure, _scStructure_light);
  //_scStructure_atoms_original = _scStructure.atoms;
  
  if(LDEBUG){
    cerr << "APL::Supercell::Supercell: this is the supercell to be analyzed" << std::endl;
    cerr << _scStructure << std::endl;
  }
}
#else
void Supercell::build(int nx, int ny, int nz, bool VERBOSE) {
  // Copy the first cell
  _scStructure = _inStructure;

  // Clear some arrays we will rebuild...
  reset();

  // Create lattice of the supercell
  xmatrix<double> scale(3, 3);
  scale.clear();
  scale(1, 1) = nx;
  scale(2, 2) = ny;
  scale(3, 3) = nz;
  _scStructure.lattice = scale * _inStructure.lattice;
  _scStructure.FixLattices();

  // Print info
  if (VERBOSE) _logger << "The supercell is going to build as " << nx << " x " << ny << " x " << nz
                       << " (" << (uint)(nx * ny * nz * _inStructure.atoms.size()) << " atoms). ";

  //
  _atom atom;
  xvector<double> cshift(3);
  for (uint ia = 0; ia < _inStructure.iatoms.size(); ia++) {
    for (uint iia = 0; iia < _inStructure.iatoms[ia].size(); iia++) {
      // Replicate this atom by given mesh...
      for (_AFLOW_APL_REGISTER_ int i = 0; i < nx; i++)
        for (_AFLOW_APL_REGISTER_ int j = 0; j < ny; j++)
          for (_AFLOW_APL_REGISTER_ int k = 0; k < nz; k++) {
            // Create position of new atoms...
            atom = _inStructure.atoms[_inStructure.iatoms[ia][iia]];
            cshift = (((double)i) * _inStructure.lattice(1) +
                      ((double)j) * _inStructure.lattice(2) +
                      ((double)k) * _inStructure.lattice(3));
            atom.cpos = atom.cpos + cshift;
            atom.fpos = C2F(_scStructure.lattice, atom.cpos);

            // Increase the number of atoms of this type...
            _scStructure.num_each_type[atom.type]++;

            // Mark this atom as equivalent or not....
            if (_scStructure.iatoms[ia].empty()) {
              atom.equivalent = _scStructure.atoms.size();
              atom.is_inequivalent = TRUE;
            } else {
              atom.equivalent = _scStructure.iatoms[ia][0];
              atom.is_inequivalent = FALSE;
            }

            // Add it to the list of all atoms...
            _scStructure.atoms.push_back(atom);

            // Add its ID number to the list of equivalent atoms of this type...
            _scStructure.iatoms[ia].push_back(_scStructure.atoms.size() - 1);

            // Add its site point group...
            _scStructure.agroup.push_back(_inStructure.agroup[_inStructure.iatoms[ia][iia]]);

            // Update our mapping arrays...
            _sc2pcMap.push_back(_inStructure.iatoms[ia][iia]);
            if (i == 0 && j == 0 && k == 0) _pc2scMap.push_back(_scStructure.atoms.size() - 1);
          }
    }
  }

  // Feed the factor group list (not efficient in this order, but we have all
  // similar operations in order just shifted...
  for (_AFLOW_APL_REGISTER_ uint l = 0; l < _inStructure.fgroup.size(); l++) {
    for (_AFLOW_APL_REGISTER_ int i = 0; i < nx; i++)
      for (_AFLOW_APL_REGISTER_ int j = 0; j < ny; j++)
        for (_AFLOW_APL_REGISTER_ int k = 0; k < nz; k++) {
          // Create position of new atoms...
          cshift = (((double)i) * _inStructure.lattice(1) +
                    ((double)j) * _inStructure.lattice(2) +
                    ((double)k) * _inStructure.lattice(3));

          // Get factor group from primitive celll
          _sym_op symOp = _inStructure.fgroup[l];

          // Add our shift(in cartesian) and transform to fractional coords
          symOp.ctau = symOp.ctau + cshift;
          symOp.ftau = C2F(_scStructure.lattice, symOp.ctau);

          // If it is out of cell, skip it
          if (symOp.ftau(1) > 1.0 - _AFLOW_APL_EPS_ ||
              symOp.ftau(2) > 1.0 - _AFLOW_APL_EPS_ ||
              symOp.ftau(3) > 1.0 - _AFLOW_APL_EPS_ ||
              symOp.ftau(1) < 0.0 - _AFLOW_APL_EPS_ ||
              symOp.ftau(2) < 0.0 - _AFLOW_APL_EPS_ ||
              symOp.ftau(3) < 0.0 - _AFLOW_APL_EPS_) continue;

          // We have to correct the Uf for each symop since we have changed the lattice...
          symOp.Uf = _scStructure.c2f * symOp.Uc * _scStructure.f2c;

          // Store it
          _scStructure.fgroup.push_back(symOp);
        }
  }

  /*
    // Old version - Feed the factor group list
    for(_AFLOW_APL_REGISTER_ uint l = 0; l < _inStructure.fgroup.size(); l++) {
    for(_AFLOW_APL_REGISTER_ int i = 0; i < nx; i++)
    for(_AFLOW_APL_REGISTER_ int j = 0; j < ny; j++)
    for(_AFLOW_APL_REGISTER_ int k = 0; k < nz; k++) {
    // Create position of new atoms...
    cshift = ( ( (double)i ) * _inStructure.lattice(1) +
    ( (double)j ) * _inStructure.lattice(2) +
    ( (double)k ) * _inStructure.lattice(3) );
    _scStructure.fgroup.push_back(_inStructure.fgroup[l]);
    (_scStructure.fgroup.back()).ctau = (_scStructure.fgroup.back()).ctau + cshift;
    (_scStructure.fgroup.back()).ftau = C2F(_scStructure.lattice,(_scStructure.fgroup.back()).ctau);
    }
    }
    */

  // Setup symmetry flags
  _scStructure.pgroup_xtal_calculated = FALSE;
  _scStructure.pgroupk_calculated = FALSE;
  _scStructure.pgroupk_xtal_calculated = FALSE;
  _scStructure.pgroup_calculated = TRUE;
  _scStructure.fgroup_calculated = TRUE;
  _scStructure.sgroup_calculated = FALSE;
  _scStructure.agroup_calculated = TRUE;

  // Setup output flags
  _scStructure.write_inequivalent_flag = TRUE;

  // Set the information about this construction
  _scStructure.info = "Supercell " + stringify(nx) + "x" + stringify(ny) + "x" + stringify(nz);

  // OK.
  if (VERBOSE) _logger << "Done." << apl::endl;
  _isConstructed = TRUE;
  //cout << _inStructure << std::endl;
  //cout << _scStructure << std::endl;

  // Bug correction for old version
  /*
      if( nx != ny || ny != nz || nx != nz ) {
      _logger << "Supercell is not symmetric, hence we need to recalculate the whole symmetry. (very slow)" << apl::endl;
      _scStructure.pgroup_calculated = FALSE;
      _scStructure.fgroup_calculated = FALSE;
      _scStructure.iatoms_calculated = FALSE;
      _scStructure.agroup_calculated = FALSE;
      calculateWholeSymmetry(_scStructure);
      }
    */
}
#endif

// ///////////////////////////////////////////////////////////////////////////

void Supercell::trimStructure(int n, const xvector<double>& a,
                              const xvector<double>& b, const xvector<double>& c, bool constructSymmetry) {
  // Clear some arrays we will rebuild...
  reset();

  // Create lattice of supercell
  _scStructure.lattice(1, 1) = a(1);
  _scStructure.lattice(2, 1) = b(1);
  _scStructure.lattice(3, 1) = c(1);
  _scStructure.lattice(1, 2) = a(2);
  _scStructure.lattice(2, 2) = b(2);
  _scStructure.lattice(3, 2) = c(2);
  _scStructure.lattice(1, 3) = a(3);
  _scStructure.lattice(2, 3) = b(3);
  _scStructure.lattice(3, 3) = c(3);
  _scStructure.FixLattices();

  // Trim from grid p x p x p, where p = 2 * n
  int p = 2 * n;
  _atom atom;
  xvector<double> cshift(3);

  for (uint ia = 0; ia < _inStructure.iatoms.size(); ia++) {
    for (uint iia = 0; iia < _inStructure.iatoms[ia].size(); iia++) {
      // Replicate this atom by given mesh...
      for (_AFLOW_APL_REGISTER_ int i = 0; i <= p; i++)
        for (_AFLOW_APL_REGISTER_ int j = 0; j <= p; j++)
          for (_AFLOW_APL_REGISTER_ int k = 0; k <= p; k++) {
            // Create position of new atoms...
            atom = _inStructure.atoms[_inStructure.iatoms[ia][iia]];
            cshift = (((double)i) * _inStructure.lattice(1) +
                      ((double)j) * _inStructure.lattice(2) +
                      ((double)k) * _inStructure.lattice(3));
            atom.cpos = atom.cpos + cshift;
            atom.fpos = C2F(_scStructure.lattice, atom.cpos);

            // Add only atoms inside the cell
            if (atom.fpos(1) > 1.0 - _AFLOW_APL_EPS_ ||
                atom.fpos(2) > 1.0 - _AFLOW_APL_EPS_ ||
                atom.fpos(3) > 1.0 - _AFLOW_APL_EPS_ ||
                atom.fpos(1) < 0.0 - _AFLOW_APL_EPS_ ||
                atom.fpos(2) < 0.0 - _AFLOW_APL_EPS_ ||
                atom.fpos(3) < 0.0 - _AFLOW_APL_EPS_) continue;

            // Increase the number of atoms of this type...
            _scStructure.num_each_type[atom.type]++;

            // Mark this atom as equivalent or not....
            if (_scStructure.iatoms[ia].empty()) {
              atom.equivalent = _scStructure.atoms.size();
              atom.is_inequivalent = TRUE;
            } else {
              atom.equivalent = _scStructure.iatoms[ia][0];
              atom.is_inequivalent = FALSE;
            }

            // Add it to the list of all atoms...
            _scStructure.atoms.push_back(atom);

            // Add its ID number to the list of equivalent atoms of this type...
            _scStructure.iatoms[ia].push_back(_scStructure.atoms.size() - 1);

            // Add its site point group...
            if (constructSymmetry) {
              _scStructure.agroup.push_back(_inStructure.agroup[_inStructure.iatoms[ia][iia]]);
              // We have to correct the Uf for each symop since we have changed the lattice...
              // Stefano formula - great help!
              vector<_sym_op>::iterator soi = _scStructure.agroup.back().begin();
              for (; soi != _scStructure.agroup.back().end(); soi++)
                soi->Uf = _scStructure.c2f * soi->Uc * _scStructure.f2c;
            }

            // Update our mapping arrays...
            _sc2pcMap.push_back(_inStructure.iatoms[ia][iia]);
            if (i == 0 && j == 0 && k == 0) _pc2scMap.push_back(_scStructure.atoms.size() - 1);
          }
    }
  }

  // Feed the factor group list (not efficient in this order, but we have all
  // similar operations in order just shifted...
  if (constructSymmetry) {
    for (_AFLOW_APL_REGISTER_ uint l = 0; l < _inStructure.fgroup.size(); l++) {
      for (_AFLOW_APL_REGISTER_ int i = 0; i < p; i++)
        for (_AFLOW_APL_REGISTER_ int j = 0; j < p; j++)
          for (_AFLOW_APL_REGISTER_ int k = 0; k < p; k++) {
            // Create position of new atoms...
            cshift = (((double)i) * _inStructure.lattice(1) +
                      ((double)j) * _inStructure.lattice(2) +
                      ((double)k) * _inStructure.lattice(3));

            _sym_op symOp = _inStructure.fgroup[l];
            symOp.ctau = symOp.ctau + cshift;
            symOp.ftau = C2F(_scStructure.lattice, symOp.ctau);

            if (symOp.ftau(1) > 1.0 - _AFLOW_APL_EPS_ ||
                symOp.ftau(2) > 1.0 - _AFLOW_APL_EPS_ ||
                symOp.ftau(3) > 1.0 - _AFLOW_APL_EPS_ ||
                symOp.ftau(1) < 0.0 - _AFLOW_APL_EPS_ ||
                symOp.ftau(2) < 0.0 - _AFLOW_APL_EPS_ ||
                symOp.ftau(3) < 0.0 - _AFLOW_APL_EPS_) continue;

            // We have to correct the Uf for each symop since we have changed the lattice...
            // Stefano formula - great help!
            symOp.Uf = _scStructure.c2f * symOp.Uc * _scStructure.f2c;

            _scStructure.fgroup.push_back(symOp);
          }
    }
  }

  // Order atoms in VASP like style
  int start = 0;
  for (_AFLOW_APL_REGISTER_ int i = 0; i < (int)_scStructure.num_each_type.size(); i++) {
    int end = start + _scStructure.num_each_type[i];
    for (_AFLOW_APL_REGISTER_ int j = start; j < end - 1; j++) {
      if (fabs(_scStructure.atoms[j].fpos(1)) < _AFLOW_APL_EPS_) _scStructure.atoms[j].fpos(1) = 0.0;
      if (fabs(_scStructure.atoms[j].fpos(2)) < _AFLOW_APL_EPS_) _scStructure.atoms[j].fpos(2) = 0.0;
      if (fabs(_scStructure.atoms[j].fpos(3)) < _AFLOW_APL_EPS_) _scStructure.atoms[j].fpos(3) = 0.0;

      for (_AFLOW_APL_REGISTER_ int k = j + 1; k < end; k++) {
        if (fabs(_scStructure.atoms[k].fpos(1)) < _AFLOW_APL_EPS_) _scStructure.atoms[k].fpos(1) = 0.0;
        if (fabs(_scStructure.atoms[k].fpos(2)) < _AFLOW_APL_EPS_) _scStructure.atoms[k].fpos(2) = 0.0;
        if (fabs(_scStructure.atoms[k].fpos(3)) < _AFLOW_APL_EPS_) _scStructure.atoms[k].fpos(3) = 0.0;

        if (_scStructure.atoms[k].fpos(1) < _scStructure.atoms[j].fpos(1)) {
          _atom ta = _scStructure.atoms[j];
          _scStructure.atoms[j] = _scStructure.atoms[k];
          _scStructure.atoms[k] = ta;

          for (_AFLOW_APL_REGISTER_ int l = 0; l < (int)_inStructure.atoms.size(); l++) {
            if (_pc2scMap[l] == j)
              _pc2scMap[l] = k;
            else if (_pc2scMap[l] == k)
              _pc2scMap[l] = j;
          }

          int ti = _sc2pcMap[j];
          _sc2pcMap[j] = _sc2pcMap[k];
          _sc2pcMap[k] = ti;
        } else if (_scStructure.atoms[k].fpos(1) == _scStructure.atoms[j].fpos(1)) {
          if (_scStructure.atoms[k].fpos(2) < _scStructure.atoms[j].fpos(2)) {
            _atom ta = _scStructure.atoms[j];
            _scStructure.atoms[j] = _scStructure.atoms[k];
            _scStructure.atoms[k] = ta;

            for (_AFLOW_APL_REGISTER_ int l = 0; l < (int)_inStructure.atoms.size(); l++) {
              if (_pc2scMap[l] == j)
                _pc2scMap[l] = k;
              else if (_pc2scMap[l] == k)
                _pc2scMap[l] = j;
            }

            int ti = _sc2pcMap[j];
            _sc2pcMap[j] = _sc2pcMap[k];
            _sc2pcMap[k] = ti;
          } else if (_scStructure.atoms[k].fpos(2) == _scStructure.atoms[j].fpos(2)) {
            if (_scStructure.atoms[k].fpos(3) < _scStructure.atoms[j].fpos(3)) {
              _atom ta = _scStructure.atoms[j];
              _scStructure.atoms[j] = _scStructure.atoms[k];
              _scStructure.atoms[k] = ta;

              for (_AFLOW_APL_REGISTER_ int l = 0; l < (int)_inStructure.atoms.size(); l++) {
                if (_pc2scMap[l] == j)
                  _pc2scMap[l] = k;
                else if (_pc2scMap[l] == k)
                  _pc2scMap[l] = j;
              }

              int ti = _sc2pcMap[j];
              _sc2pcMap[j] = _sc2pcMap[k];
              _sc2pcMap[k] = ti;
            }
          }
        }
      }
    }
    start += _scStructure.num_each_type[i];
  }

  // Setup symmetry flags
  if (constructSymmetry) {
    _scStructure.pgroup_xtal_calculated = FALSE;
    _scStructure.pgroupk_calculated = FALSE;
    _scStructure.pgroupk_xtal_calculated = FALSE;
    _scStructure.pgroup_calculated = TRUE;
    _scStructure.fgroup_calculated = TRUE;
    _scStructure.sgroup_calculated = FALSE;
    _scStructure.agroup_calculated = TRUE;
  } else {
    _scStructure.pgroup_xtal_calculated = FALSE;
    _scStructure.pgroupk_calculated = FALSE;
    _scStructure.pgroupk_xtal_calculated = FALSE;
    _scStructure.pgroup_calculated = FALSE;
    _scStructure.fgroup_calculated = FALSE;
    _scStructure.sgroup_calculated = FALSE;
    _scStructure.agroup_calculated = TRUE;
  }

  // Setup output flags
  _scStructure.write_inequivalent_flag = TRUE;

  // Set the information about this construction
  _scStructure.info = "Trimed supercell";

  // OK.
  //_logger << "Done." << apl::endl;
  _isConstructed = TRUE;
  //cout << _scStructure << std::endl;
}

// ///////////////////////////////////////////////////////////////////////////

int Supercell::buildSuitableForShell(int MIN_NN_SHELLS, bool shouldBeFullShell, bool VERBOSE) {
  // What is the dimension of the supercell ? OK, user wants to have MAX_NN_SHELLS
  // shell occupied for each nonequvalent atom. Try to find it...

  // Precompute shellhandles for each unique atom
  vector<ShellHandle> sh;
  bool useSplittedShells = true;
  for (uint i = 0; i < _inStructure.iatoms.size(); i++) {
    ShellHandle s;
    sh.push_back(s);
    sh.back().init(_inStructure, _inStructure.iatoms[i][0], MIN_NN_SHELLS);
    try {
      if (useSplittedShells)
        sh.back().splitBySymmetry();
    } catch (APLLogicError& e) {
      //COREY, we may want to kill this if we create supercells of uniform expansion (not derivative_structures)
      //come back to fix later, but for now, leave as is (NOT WELL TESTED)
      //also, see exit below (throw error), this indicates to me that we do not NEED to exit, but can proceed until next error
      //corey, kill if errors with symmetry
      //_logger << apl::error << e.what() << apl::endl;
      //_logger << apl::error << "The splitting of shells by the symmetry has failed [" << i << "]." << apl::endl;
      //throw APLRuntimeError("apl::Supercell::buildSuitableForShell(); Symmetry failed.");
      _logger << apl::error << e.what() << apl::endl;
      _logger << apl::warning << "The splitting of shells by the symmetry has failed [" << i << "]. Continuing without this..." << apl::endl;
      useSplittedShells = false;
      for (uint j = 0; j < sh.size(); j++) {
        sh[j].removeSplitBySymmetry();
      }
    }
  }
  _maxShellID = MIN_NN_SHELLS;

  _AFLOW_APL_REGISTER_ int i = 1;
  _AFLOW_APL_REGISTER_ int j = 1;
  _AFLOW_APL_REGISTER_ int k = 1;

  while (true) {
    _scStructure = GetSuperCell(_inStructure, i, j, k);
    //build(i,j,k);

    uint ia = 0;
    for (; ia < _inStructure.iatoms.size(); ia++) {
      // Get ID of origin atom from pc in the sc
      uint l = 0;
      for (; l < _scStructure.atoms.size(); l++) {
        if (aurostd::modulus(_scStructure.atoms[l].cpos -
                             _inStructure.atoms[_inStructure.iatoms[ia][0]].cpos) < _AFLOW_APL_EPS_)
          break;
      }
      if (l == _scStructure.atoms.size()) {
        throw APLRuntimeError("apl::Supercell::buildSuitableForShell(); Mapping error.");
      }

      // Map with this center
      sh[ia].mapStructure(_scStructure, l, useSplittedShells);

      // Get last shell
      int lastShell;
      if (shouldBeFullShell)
        lastShell = sh[ia].getLastFullShell();
      else
        lastShell = sh[ia].getLastRegularShell();

      if (lastShell < MIN_NN_SHELLS) break;
    }

    if (ia != _inStructure.iatoms.size()) {
      i++;
      j = k = i;
      /*
	  if( i == j && j == k )
	  i++;
	  else if( j == k )
	  j++;
	  else
	  k++;
	*/
    } else {
      break;
    }
  }

  // Build structure
  build(i, j, k, VERBOSE);

  // Print info about shells
  for (uint i = 0; i < _inStructure.iatoms.size(); i++) {
    if (VERBOSE) sh[i].printReport(cout);
    sh[i].clear();
  }
  sh.clear();

  return i * j * k * _inStructure.atoms.size();
}

// ///////////////////////////////////////////////////////////////////////////

void Supercell::setupShellRestrictions(int MAX_NN_SHELLS) {
  // Precompute shellhandles for each unique atom
  vector<ShellHandle> sh;
  for (uint i = 0; i < _inStructure.iatoms.size(); i++) {
    ShellHandle s;
    sh.push_back(s);
    sh.back().init(_inStructure, _inStructure.iatoms[i][0], MAX_NN_SHELLS);
    sh.back().mapStructure(_scStructure, _scStructure.iatoms[i][0]);
  }
  _maxShellID = MAX_NN_SHELLS;

  // Set flag to shell restriction
  _isShellRestricted = true;
  _logger << "Setting shell restrictions up to " << MAX_NN_SHELLS << "." << apl::endl;

  // Calculate the truncate radius for each atom
  _maxShellRadius.clear();
  for (uint i = 0; i < _inStructure.iatoms.size(); i++) {
    double r = sh[i].getShellRadius(MAX_NN_SHELLS);
    for (uint j = 0; j < _inStructure.iatoms[i].size(); j++) {
      _maxShellRadius.push_back(r);
    }
  }

  // Print info about shells
  for (uint i = 0; i < _inStructure.iatoms.size(); i++)
    sh[i].clear();
  sh.clear();
}

// ///////////////////////////////////////////////////////////////////////////

bool Supercell::isShellRestricted() {
  return _isShellRestricted;
}

int Supercell::getMaxShellID() {
  return _maxShellID;
}

// ///////////////////////////////////////////////////////////////////////////

xvector<double> Supercell::getFPositionItsNearestImage(const xvector<double>& fposAtom,
                                                       const xvector<double>& fposCenter,
                                                       const xmatrix<double>& lattice) {
  double r2min = numeric_limits<double>::max();
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

xvector<double> Supercell::getFPositionItsNearestImage(int atomID, int centerID) {
  return (getFPositionItsNearestImage(_scStructure.atoms[atomID].fpos,
                                      _scStructure.atoms[centerID].fpos,
                                      _scStructure.lattice));
}

// ///////////////////////////////////////////////////////////////////////////

xvector<double> Supercell::getCPositionItsNearestImage(int atomID, int centerID) {
  return (F2C(_scStructure.lattice,
              getFPositionItsNearestImage(atomID, centerID)));
}

// ///////////////////////////////////////////////////////////////////////////

bool Supercell::isConstructed() {
  return _isConstructed;
}

// ///////////////////////////////////////////////////////////////////////////

const xstructure& Supercell::getSupercellStructure() const {
  return _scStructure;
}

// ///////////////////////////////////////////////////////////////////////////

//CO, only what's necessary (no heavy symmetry copies)
const xstructure& Supercell::getSupercellStructureLight() const {
  return _scStructure_light;
  //xstructure a;
  //stringstream POSCAR;
  //POSCAR.str("");
  //POSCAR << _scStructure;
  //POSCAR >> a;
  //enable inequivalent flag to work
  //for(uint i=0;i<a.atoms.size();i++){
  //  a.atoms[i].equivalent=_scStructure.atoms[i].equivalent;
  //  a.atoms[i].is_inequivalent=_scStructure.atoms[i].is_inequivalent;
  //  a.atoms[i].num_equivalents=_scStructure.atoms[i].num_equivalents;
  //}
  //pseudo-potential stuff
  //for(uint i=0;i<a.species.size();i++){
  //  a.species[i]=_scStructure.species[i];
  //  a.species_pp[i]=_scStructure.species_pp[i]; //VERY IMPORTANT
  //}
  //enable inequivalent flag to work
  //a.write_inequivalent_flag = _scStructure.write_inequivalent_flag;
  //a.info = _scStructure.info;
  //return a;
}
// ///////////////////////////////////////////////////////////////////////////

xstructure Supercell::calculatePrimitiveStructure() const { //CO 180409
  xstructure pcStructure=_inStructure;
  pcStructure.Standard_Primitive_UnitCellForm();
  pcStructure.ReScale(1.0);
  pcStructure.ShifOriginToAtom(0);
  pcStructure.BringInCell();
  return pcStructure;
}

const xstructure& Supercell::getPrimitiveStructure() const {
  return _pcStructure;
}

const xstructure& Supercell::getInputStructure() const {
  return _inStructure;
}

// ///////////////////////////////////////////////////////////////////////////

const xstructure& Supercell::getInputStructureLight() const {
  return _inStructure_light;
  //xstructure a;
  //stringstream POSCAR;
  //POSCAR.str("");
  //POSCAR << _inStructure;
  //POSCAR >> a;
  //enable inequivalent flag to work
  //for(uint i=0;i<a.atoms.size();i++){
  //  a.atoms[i].equivalent=_inStructure.atoms[i].equivalent;
  //  a.atoms[i].is_inequivalent=_inStructure.atoms[i].is_inequivalent;
  //  a.atoms[i].num_equivalents=_inStructure.atoms[i].num_equivalents;
  //}
  //enable inequivalent flag to work
  //a.write_inequivalent_flag = _inStructure.write_inequivalent_flag;
  //a.info = _inStructure.info;
  //return a;
}

// ///////////////////////////////////////////////////////////////////////////

int Supercell::getNumberOfAtoms() {
  return _scStructure.atoms.size();
}

// ///////////////////////////////////////////////////////////////////////////

int Supercell::getNumberOfUniqueAtoms() {
  return _scStructure.iatoms.size();
}

// ///////////////////////////////////////////////////////////////////////////

int Supercell::getNumberOfUniqueAtomsOfType(int i) {
#ifndef __OPTIMIZE
  if (i >= (int)_scStructure.iatoms.size())
    throw APLRuntimeError("apl::Supercell::getNumberOfUniqueAtomsOfType: Wrong index.");
#endif
  return _scStructure.iatoms[i].size();
}

// ///////////////////////////////////////////////////////////////////////////

int Supercell::getUniqueAtomID(int i) {
#ifndef __OPTIMIZE
  if (i >= (int)_scStructure.iatoms.size())
    throw APLRuntimeError("apl::Supercell::getUniqueAtoms(): Wrong index.");
#endif
  return _scStructure.iatoms[i][0];
}

// ///////////////////////////////////////////////////////////////////////////

int Supercell::getUniqueAtomID(int i, int j) {
#ifndef __OPTIMIZE
  if (i >= (int)_scStructure.iatoms.size())
    throw APLRuntimeError("apl::Supercell::getUniqueAtoms(): Wrong index 1.");

  if (j >= (int)_scStructure.iatoms[i].size())
    throw APLRuntimeError("apl::Supercell::getUniqueAtoms(): Wrong index 2.");
#endif
  return _scStructure.iatoms[i][j];
}

// ///////////////////////////////////////////////////////////////////////////

const _atom& Supercell::getUniqueAtom(int i) {
  return _scStructure.atoms[getUniqueAtomID(i)];
}

// ///////////////////////////////////////////////////////////////////////////

#if !JAHNATEK_ORIGINAL
bool Supercell::compareFPositions(xvector<double>& v1, xvector<double>& v2) {
  //default assumption is that we are compare positions associated with symmetry rotations
  double eps = _sym_eps;
  return compareFPositions(v1, v2, eps);
}
bool Supercell::compareFPositions(xvector<double>& v1, xvector<double>& v2, double eps) {
  // Get the difference vector for SUPERCELL positions
  // if symmetry related, use eps=_sym_eps (default), otherwise eps=_AFLOW_APL_EPS_
  return SYM::AtomFPOSMatch(v1, v2, _scStructure.c2f, _scStructure.f2c, _skew, eps);
}
#else
bool Supercell::compareFPositions(const xvector<double>& v1,
                                  const xvector<double>& v2, double eps) {
  // Get the difference vector
  xvector<double> r = v1 - v2;

  // Correct and check equality
  for (_AFLOW_APL_REGISTER_ int i = 1; i <= 3; i++) {
    // Correct it for rare cases, when structure is not well relaxed or
    // there is a lot of roundoff problems, like positions like
    // [0,0,0] and [0,0,0.9993567989], but they are the equal in principle
    if (fabs(v1(i) - 1.0) < eps)
      r(i) -= 1.0;
    if (fabs(v2(i) - 1.0) < eps)
      r(i) += 1.0;

    // If this component is still nonzero -> this two possitions are not
    // the same
    if (fabs(r(i)) > eps)
      return false;
  }

  //
  return true;
}
#endif

// ///////////////////////////////////////////////////////////////////////////

/******************************************************************************/

#if !JAHNATEK_ORIGINAL

/******************************************************************************/
#if !GETFULLSYMBASIS
/******************************************************************************/
int Supercell::where(int atomID, int centerID, const _sym_op& symOp, bool translate) {
//corey
//change so that if we can retrieve from fullsymbasis,  we do so
//this functions looks at symop, and asks by applying it, which atom does atomID become?
//in fgroup, look at basis_atoms_map, return atom at index atomID
#ifndef __OPTIMIZE
  if (atomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::where(); Wrong atomID index.");

  if (centerID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::where(); Wrong centerID index.");
#endif

  // Get the center atom center...
  if (translate && symOp.is_agroup) center(centerID);

  // Transform atom...
  //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[atomID], symOp,
  //DX                                   _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
  _atom rotatedAtom;
  if (!SYM::ApplyAtomValidate(_scStructure.atoms[atomID], rotatedAtom, symOp, _scStructure, _skew, TRUE, FALSE, _sym_eps)) {
    throw APLLogicError("apl::Supercell::where(); Illegitimate mapping.");
  }

  // Find its id...
  _AFLOW_APL_REGISTER_ int l = 0;
  for (; l < (int)_scStructure.atoms.size(); l++) {
    if (compareFPositions(rotatedAtom.fpos, _scStructure.atoms[l].fpos)) {  //COREY NEW, default to symmetry tolerance
      break;
    }
  }

  if (l == (int)_scStructure.atoms.size()) {
#if ERROR_VERBOSE
    _AFLOW_APL_REGISTER_ int l = 0;
    for (; l < (int)_scStructure.atoms.size(); l++) {
      printXVector(_scStructure.atoms[atomID].fpos, false);
      cout << " -> ";
      printXVector(rotatedAtom.fpos, false);
      cout << " | ";
      printXVector(_scStructure.atoms[l].fpos, false);
      cout << " | ";
      cout << aurostd::modulus(rotatedAtom.fpos - _scStructure.atoms[l].fpos) << std::endl;
    }
#endif
    throw APLLogicError("apl::Supercell::where(); Mapping failed.");
  }

  // Move center back to zero atom...
  //if( translate && symOp.is_agroup ) center(0);
  if (translate && symOp.is_agroup) center_original();

#if MAP_VERBOSE
  bool will_translate = symOp.is_agroup && translate;
  cerr << "where: " << atomID << " " << centerID << " " << will_translate << " " << l << std::endl;
  cerr << "atomID : " << _scStructure.atoms[atomID] << std::endl;
  cerr << symOp << std::endl;
#endif

  return l;
}
/******************************************************************************/
#else
/******************************************************************************/
int Supercell::where(int atomID, int centerID, const _sym_op& symOp, bool translate) {
//corey
//change so that if we can retrieve from fullsymbasis,  we do so
//this functions looks at symop, and asks by applying it, which atom does atomID become?
//in fgroup, look at basis_atoms_map, return atom at index atomID
#ifndef __OPTIMIZE
  if (atomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::where(); Wrong atomID index.");

  if (centerID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::where(); Wrong centerID index.");
#endif
//corey
#if MAP_VERBOSE
  bool will_translate = symOp.is_agroup && translate;
  cerr << "where: " << atomID << " " << centerID << " " << will_translate << " " << symOp.basis_atoms_map.at(atomID) << std::endl;
  cerr << "atomID : " << _scStructure.atoms[atomID] << std::endl;
  cerr << symOp << std::endl;
#endif
  return symOp.basis_atoms_map.at(atomID);
}
#endif
/******************************************************************************/

#else

/******************************************************************************/
int Supercell::where(int atomID, int centerID, const _sym_op& symOp, bool translate) {
#ifndef __OPTIMIZE
  if (atomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::where(); Wrong atomID index.");

  if (centerID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::where(); Wrong centerID index.");
#endif

  // Get the center atom center...
  if (translate && symOp.is_agroup) center(centerID);

  // Transform atom...
  //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[atomID], symOp,
  //DX                                    _scStructure, TRUE);
  _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[atomID], symOp,
                                     _scStructure);

  // Find its id...
  _AFLOW_APL_REGISTER_ int l = 0;
  for (; l < (int)_scStructure.atoms.size(); l++) {
    if (compareFPositions(rotatedAtom.fpos, _scStructure.atoms[l].fpos, _SYM_AFLOW_APL_EPS_))
      break;
  }

  if (l == (int)_scStructure.atoms.size()) {
    /*
	_AFLOW_APL_REGISTER_ int l = 0;
	for(; l < (int)_scStructure.atoms.size(); l++) {
	printXVector(_scStructure.atoms[atomID].fpos,false); cout << " -> ";
	printXVector(rotatedAtom.fpos,false); cout << " | ";
	printXVector(_scStructure.atoms[l].fpos,false); cout << " | ";
	cout << aurostd::modulus( rotatedAtom.fpos - _scStructure.atoms[l].fpos ) << std::endl;
	}
      */
    throw APLLogicError("apl::Supercell::where(); Mapping failed.");
  }

  // Move center back to zero atom...
  if (translate && symOp.is_agroup) center(0);

  return l;
}
/******************************************************************************/

#endif

/******************************************************************************/

// ///////////////////////////////////////////////////////////////////////////

/******************************************************************************/

#if !JAHNATEK_ORIGINAL

/******************************************************************************/
#if !GETFULLSYMBASIS
/******************************************************************************/
int Supercell::wherefrom(int atomID, int centerID, const _sym_op& symOp, bool translate) {
//corey
//this function does the opposite (to above)
//in basis_atoms_map, return the index of the atom that is atomID
#ifndef __OPTIMIZE
  if (atomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::wherefrom(); Wrong atomID index.");

  if (centerID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::wherefrom(); Wrong centerID index.");
#endif

  // Get the center atom center...
  if (translate && symOp.is_agroup) center(centerID);

  // Find it
  int l = 0;
  _atom rotatedAtom;  //CO
  for (; l < (int)_scStructure.atoms.size(); l++) {
    //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
    if (!SYM::ApplyAtomValidate(_scStructure.atoms[l], rotatedAtom, symOp, _scStructure, _skew, TRUE, FALSE, _sym_eps)) {
      throw APLLogicError("apl::Supercell::wherefrom(); Illegitimate mapping.");
    }
    if (compareFPositions(rotatedAtom.fpos, _scStructure.atoms[atomID].fpos)) {  //COREY NEW, default to symmetry tolerance
      break;
    }
  }

  if (l == (int)_scStructure.atoms.size()) {
#if ERROR_VERBOSE
    _AFLOW_APL_REGISTER_ int l = 0;
    for (; l < (int)_scStructure.atoms.size(); l++) {
      //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
      _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE, FALSE);  //CO no roff
      printXVector(rotatedAtom.fpos, false);
      cout << " | ";
      printXVector(_scStructure.atoms[atomID].fpos, false);
      cout << " | ";
      cout << aurostd::modulus(rotatedAtom.fpos - _scStructure.atoms[atomID].fpos) << std::endl;
    }
#endif
    throw APLLogicError("apl::Supercell::wherefrom(); Mapping failed.");
  }

  // Move center back to zero atom...
  //if( translate && symOp.is_agroup ) center(0);
  if (translate && symOp.is_agroup) center_original();

#if MAP_VERBOSE
  bool will_translate = symOp.is_agroup && translate;
  cerr << "wherefrom: " << atomID << " " << centerID << " " << will_translate << " " << l << std::endl;
  cerr << "atomID : " << _scStructure.atoms[atomID] << std::endl;
  cerr << symOp << std::endl;
#endif

  return l;
}
/******************************************************************************/
#else
/******************************************************************************/
int Supercell::wherefrom(int atomID, int centerID, const _sym_op& symOp, bool translate) {
//corey
//this function does the opposite (to above)
//in basis_atoms_map, return the index of the atom that is atomID
#ifndef __OPTIMIZE
  if (atomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::wherefrom(); Wrong atomID index.");

  if (centerID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::wherefrom(); Wrong centerID index.");
#endif

  int l = 0;
  //corey
  for (; l < (int)symOp.basis_atoms_map.size(); l++) {
    if (symOp.basis_atoms_map.at(l) == atomID) {
#if MAP_VERBOSE
      bool will_translate = symOp.is_agroup && translate;
      cerr << "wherefrom: " << atomID << " " << centerID << " " << will_translate << " " << l << std::endl;
      cerr << "atomID : " << _scStructure.atoms[atomID] << std::endl;
      cerr << symOp << std::endl;
#endif
      return l;
    }
  }

#if ERROR_VERBOSE
  int l = 0;
  for (; l < (int)symOp.basis_atoms_map.size(); l++) {
    cout << "CHECKING: atomID=" << atomID << ", symOp.basis_atoms_map.at(l)=" << symOp.basis_atoms_map.at(l) << std::endl;
  }
#endif
  throw APLLogicError("apl::Supercell::wherefrom(); Mapping failed.");
}
#endif
/******************************************************************************/

#else

/******************************************************************************/
int Supercell::wherefrom(int atomID, int centerID, const _sym_op& symOp, bool translate) {
#ifndef __OPTIMIZE
  if (atomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::wherefrom(); Wrong atomID index.");

  if (centerID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::wherefrom(); Wrong centerID index.");
#endif

  // Get the center atom center...
  if (translate && symOp.is_agroup) center(centerID);

  // Find it
  int l = 0;
  for (; l < (int)_scStructure.atoms.size(); l++) {
    //DX_atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE);
    _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure);
    if (compareFPositions(rotatedAtom.fpos, _scStructure.atoms[atomID].fpos, _SYM_AFLOW_APL_EPS_))
      break;
  }

  if (l == (int)_scStructure.atoms.size()) {
    _AFLOW_APL_REGISTER_ int l = 0;
    for (; l < (int)_scStructure.atoms.size(); l++) {
      //DX _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure, TRUE);
      _atom rotatedAtom = SYM::ApplyAtom(_scStructure.atoms[l], symOp, _scStructure);
      printXVector(rotatedAtom.fpos, false);
      cout << " | ";
      printXVector(_scStructure.atoms[atomID].fpos, false);
      cout << " | ";
      cout << aurostd::modulus(rotatedAtom.fpos - _scStructure.atoms[atomID].fpos) << std::endl;
    }
    throw APLLogicError("apl::Supercell::wherefrom(); Mapping failed.");
  }

  // Move center back to zero atom...
  if (translate && symOp.is_agroup) center(0);

  return l;
}
/******************************************************************************/

#endif

/******************************************************************************/

// ///////////////////////////////////////////////////////////////////////////

/******************************************************************************/

#if !JAHNATEK_ORIGINAL

/******************************************************************************/
#if !GETFULLSYMBASIS
/******************************************************************************/
const _sym_op& Supercell::getSymOpWhichMatchAtoms(int whichAtomID, int toAtomID, int GROUP) {
//go through all fgroups, look at basis_atoms_map at index whichatomID, find toAtomID
#ifndef __OPTIMIZE
  if (whichAtomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Wrong atom1ID index.");

  if (toAtomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Wrong atom2ID index.");
#endif

  vector<_sym_op>* symPool = NULL;
  if (GROUP == _PGROUP_)
    symPool = &_scStructure.pgroup;
  else if (GROUP == _FGROUP_)
    symPool = &_scStructure.fgroup;
  else
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Unknown group type.");

  // Apply all symmetry operations on atom1 and find which one produce atom2
  uint iSymOp = 0;
  _atom newAtom;  //CO
  for (; iSymOp < symPool->size(); iSymOp++) {
    //DX _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[iSymOp],
    //DX                               _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
    if (!SYM::ApplyAtomValidate(_scStructure.atoms[whichAtomID], newAtom, (*symPool)[iSymOp], _scStructure, _skew, TRUE, FALSE, _sym_eps)) {
      throw APLLogicError("apl::Supercell::getSymOpWhichMatchAtoms(); Illegitimate mapping.");
    }
    if (compareFPositions(newAtom.fpos, _scStructure.atoms[toAtomID].fpos)) {  //COREY NEW, default to symmetry tolerance
      break;
    }
  }

  if (iSymOp == symPool->size()) {
#if ERROR_VERBOSE
    //corey some helpful output
    cout << "_scStructure.fgroup.size()=" << _scStructure.fgroup.size() << std::endl;
    cout << "symPool->size()=" << symPool->size() << std::endl;
    cout << "ATOM " << whichAtomID << ": " << _scStructure.atoms[whichAtomID].fpos << std::endl;
    cout << "ATOM " << toAtomID << ": " << _scStructure.atoms[toAtomID].fpos << std::endl;

    _AFLOW_APL_REGISTER_ uint l = 0;
    for (; l < symPool->size(); l++) {
      cout << "i=" << l << std::endl;
      cout << (*symPool)[l] << std::endl;
      printXVector(_scStructure.atoms[whichAtomID].fpos, false);
      cout << " -> ";
      //DX _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[l],
      //DX                               _scStructure, TRUE, FALSE, _derivative_structure);  //CO no roff
      _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[l],
                                     _scStructure, TRUE, FALSE);  //CO no roff
      printXVector(newAtom.fpos, false);
      cout << " | ";
      printXVector(_scStructure.atoms[toAtomID].fpos, false);
      cout << " | ";
      cout << aurostd::modulus(newAtom.fpos - _scStructure.atoms[toAtomID].fpos) << std::endl;
    }
#endif
    throw APLLogicError("apl::Supercell::getSymOpWhichMatchAtoms(); Mapping failed.");
  }

#if MAP_VERBOSE
  cerr << "getSymOp: " << whichAtomID << " " << toAtomID << " " << std::endl;
  cerr << "whichAtomID: " << _scStructure.atoms[whichAtomID] << std::endl;
  cerr << "toAtomID: " << _scStructure.atoms[toAtomID] << std::endl;
  cerr << (*symPool)[iSymOp] << std::endl;
#endif

  return (*symPool)[iSymOp];
}
/******************************************************************************/
#else
/******************************************************************************/
const _sym_op& Supercell::getSymOpWhichMatchAtoms(int whichAtomID, int toAtomID, int GROUP) {
//go through all groups, look at basis_atoms_map at index whichatomID, find toAtomID
#ifndef __OPTIMIZE
  if (whichAtomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Wrong atom1ID index.");

  if (toAtomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Wrong atom2ID index.");
#endif

  vector<_sym_op>* symPool = NULL;
  if (GROUP == _PGROUP_)
    symPool = &_scStructure.pgroup;
  else if (GROUP == _FGROUP_)
    symPool = &_scStructure.fgroup;
  else
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Unknown group type.");

  uint iSymOp = 0;
  for (; iSymOp < symPool->size(); iSymOp++) {
    if ((*symPool)[iSymOp].basis_atoms_map.at(whichAtomID) == toAtomID) {
      return (*symPool)[iSymOp];
    }
  }
  throw APLLogicError("apl::Supercell::getSymOpWhichMatchAtoms(); Mapping failed.");

#if MAP_VERBOSE
  cerr << "getSymOp: " << whichAtomID << " " << toAtomID << " " << std::endl;
  cerr << "whichAtomID: " << _scStructure.atoms[whichAtomID] << std::endl;
  cerr << "toAtomID: " << _scStructure.atoms[toAtomID] << std::endl;
  cerr << (*symPool)[iSymOp] << std::endl;
#endif

  return (*symPool)[iSymOp];
}
#endif
/******************************************************************************/

#else

/******************************************************************************/
const _sym_op& Supercell::getSymOpWhichMatchAtoms(int whichAtomID, int toAtomID, int GROUP) {
#ifndef __OPTIMIZE
  if (whichAtomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Wrong atom1ID index.");

  if (toAtomID >= (int)_scStructure.atoms.size())
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Wrong atom2ID index.");
#endif

  vector<_sym_op>* symPool = NULL;
  if (GROUP == _PGROUP_)
    symPool = &_scStructure.pgroup;
  else if (GROUP == _FGROUP_)
    symPool = &_scStructure.fgroup;
  else
    throw APLRuntimeError("apl::Supercell::getSymOpWhichMatchAtoms(); Unknown group type.");

  // Apply all symmetry operations on atom1 and find which one produce atom2
  uint iSymOp = 0;
  for (; iSymOp < symPool->size(); iSymOp++) {
    //DX _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[iSymOp],
    //DX                                _scStructure, TRUE);
    _atom newAtom = SYM::ApplyAtom(_scStructure.atoms[whichAtomID], (*symPool)[iSymOp],
                                   _scStructure);
    if (compareFPositions(newAtom.fpos, _scStructure.atoms[toAtomID].fpos, _SYM_AFLOW_APL_EPS_))
      break;
  }

  if (iSymOp == symPool->size())
    throw APLLogicError("apl::Supercell::getSymOpWhichMatchAtoms(); Mapping failed.");

  return (*symPool)[iSymOp];
}
/******************************************************************************/

#endif

/******************************************************************************/

// ///////////////////////////////////////////////////////////////////////////

int Supercell::pc2scMap(int i) {
  return _pc2scMap[i];
}

int Supercell::sc2pcMap(int i) {
  return _sc2pcMap[i];
}

// ///////////////////////////////////////////////////////////////////////////

//CO - START
void Supercell::center(int i) {
  xvector<double> origin(3), frigin(3);
#if CENTER_PRIM
  //_inStructure.ShifOriginToAtom(_sc2pcMap[i]);
  //_inStructure.BringInCell();
  origin = _inStructure_original.atoms[_sc2pcMap[i]].cpos;
  frigin = _inStructure_original.atoms[_sc2pcMap[i]].fpos;
  _inStructure.origin = origin;
  for (uint i = 0; i < _inStructure.atoms.size(); i++) {
    _inStructure.atoms[i].cpos = _inStructure_original.atoms[i].cpos - origin;
    _inStructure.atoms[i].fpos = _inStructure_original.atoms[i].fpos - frigin;
    _inStructure.atoms[i] = BringInCell(_inStructure.atoms[i], _inStructure.lattice);
  }
#endif
  //_scStructure.ShifOriginToAtom(i);
  //_scStructure.BringInCell();
  origin = _scStructure_original.atoms[i].cpos;
  frigin = _scStructure_original.atoms[i].fpos;
  _scStructure.origin = origin;
  for (uint i = 0; i < _scStructure.atoms.size(); i++) {
    _scStructure.atoms[i].cpos = _scStructure_original.atoms[i].cpos - origin;
    _scStructure.atoms[i].fpos = _scStructure_original.atoms[i].fpos - frigin;
    _scStructure.atoms[i] = BringInCell(_scStructure.atoms[i], _scStructure.lattice);
  }
}

void Supercell::center_original(void) {
//just a (fast) undo for center(atom);
//refer to ShifOriginToAtom in case more properties need to be updated
#if CENTER_PRIM
  _inStructure.origin = _inStructure_original.origin;
  for (uint i = 0; i < _inStructure.atoms.size(); i++) {
    _inStructure.atoms[i].fpos = _inStructure_original.atoms[i].fpos;
    _inStructure.atoms[i].cpos = _inStructure_original.atoms[i].cpos;
  }
#endif
  _scStructure.origin = _scStructure_original.origin;
  for (uint i = 0; i < _scStructure.atoms.size(); i++) {
    _scStructure.atoms[i].fpos = _scStructure_original.atoms[i].fpos;
    _scStructure.atoms[i].cpos = _scStructure_original.atoms[i].cpos;
  }
}

// ///////////////////////////////////////////////////////////////////////////

vector<_sym_op> Supercell::getFGROUP(void) {
  return _scStructure.fgroup;
}

// ///////////////////////////////////////////////////////////////////////////

vector<vector<_sym_op> > Supercell::getAGROUP(void) {
  return _scStructure.agroup;
}
vector<_sym_op> Supercell::getAGROUP(int i) {
  return _scStructure.agroup[i];
}

// ///////////////////////////////////////////////////////////////////////////

double Supercell::getEPS(void) {
  return _sym_eps;
}

// ///////////////////////////////////////////////////////////////////////////

bool Supercell::isDerivativeStructure(void) {
  return _derivative_structure;
}
//CO - END

// ///////////////////////////////////////////////////////////////////////////

string Supercell::getUniqueAtomSymbol(int i) {
#ifndef __OPTIMIZE
  if (i >= (int)_scStructure.iatoms.size()) {
    throw APLRuntimeError("apl::Supercell::getUniqueAtomSymbol(): Wrong index.");
  }
#endif
  return _scStructure.atoms[_scStructure.iatoms[i][0]].cleanname;
}

// ///////////////////////////////////////////////////////////////////////////

double Supercell::getUniqueAtomMass(int i) {
#ifndef __OPTIMIZE
  if (i >= (int)_scStructure.iatoms.size()) {
    throw APLRuntimeError("apl::Supercell::getUniqueAtomMass(): Wrong index.");
  }
#endif
  //return( GetAtomMass(_scStructure.atoms[_scStructure.iatoms[i][0]].cleanname) * KILOGRAM2AMU ); //JAHNATEK ORIGINAL
  //COREY - START
  double mass = GetAtomMass(_scStructure.atoms[_scStructure.iatoms[i][0]].cleanname);
  if (mass == NNN) {
    throw APLRuntimeError("apl::Supercell::getUniqueAtomMass(): Unknown atom types.");
  }
  return mass * KILOGRAM2AMU;
  //COREY - END
}

// ///////////////////////////////////////////////////////////////////////////

double Supercell::getAtomMass(int i) {
#ifndef __OPTIMIZE
  if (i >= (int)_scStructure.atoms.size()) {
    throw APLRuntimeError("apl::Supercell::getAtomMass(): Wrong index.");
  }
#endif
  return (GetAtomMass(_scStructure.atoms[i].cleanname) * KILOGRAM2AMU);
}

// ///////////////////////////////////////////////////////////////////////////

int Supercell::getAtomNumber(int i) {
#ifndef __OPTIMIZE
  if (i >= (int)_scStructure.atoms.size()) {
    throw APLRuntimeError("apl::Supercell::getAtomNumber(): Wrong index.");
  }
#endif
  return ((int)GetAtomNumber(_scStructure.atoms[i].cleanname));
}

// ///////////////////////////////////////////////////////////////////////////

bool Supercell::calcShellPhaseFactor(int atomID, int centerID, const xvector<double>& qpoint, xcomplex<double>& phase) {
  // Get the nearest image of this atom, it determine also the shell radius
  xvector<double> rf = getFPositionItsNearestImage(atomID, centerID);
  xvector<double> rc = F2C(_scStructure.lattice, rf);
  double rshell = aurostd::modulus(rc);

  //
  bool isCountable = false;
  xcomplex<double> iONE(0.0, 1.0);
  phase = 0.0;

  // Count in all atoms which sit on this shell
  if (!_isShellRestricted || (rshell <= _maxShellRadius[centerID] + _AFLOW_APL_EPS_)) {
    isCountable = true;
    int ifound = 0;
    xvector<double> pf(3), pc(3);

    for (_AFLOW_APL_REGISTER_ int ii = -1; ii <= 1; ii++)
      for (_AFLOW_APL_REGISTER_ int jj = -1; jj <= 1; jj++)
        for (_AFLOW_APL_REGISTER_ int kk = -1; kk <= 1; kk++) {
          pf(1) = rf(1) + ii;
          pf(2) = rf(2) + jj;
          pf(3) = rf(3) + kk;
          pc = F2C(_scStructure.lattice, pf);

          if (fabs(aurostd::modulus(pc) - rshell) < _AFLOW_APL_EPS_) {
            //cout << centerID+1 << " " << atomID+1 << " " << aurostd::modulus(pc) << " "; printXVector(pc);
            phase = phase + exp(iONE * scalar_product(qpoint, pc));
            ifound++;
          }
        }
    phase *= (1.0 / ifound);
  }

  //
  return isCountable;
}

// ///////////////////////////////////////////////////////////////////////////
// ***** BEGIN JJPR ***** Phase factor and derivative too

xcomplex<double> Supercell::calcShellPhaseFactor2(int atomID, int centerID, const xvector<double>& qpoint) {
  // Get the nearest image of this atom, it determine also the shell radius
  xvector<double> rf = getFPositionItsNearestImage(atomID, centerID);
  xvector<double> rc = F2C(_scStructure.lattice, rf);

  xcomplex<double> iONE(0.0, 1.0);
  xcomplex<double> phase;
  phase = exp(iONE * scalar_product(qpoint, rc));
  return phase;
}

bool Supercell::calcShellPhaseFactorAAPL(int atomID, int centerID, const xvector<double>& qpoint, xcomplex<double>& phase, xvector<xcomplex<double> >& deriv) {
  // Get the nearest image of this atom, it determine also the shell radius
  xvector<double> rf = getFPositionItsNearestImage(atomID, centerID);
  xvector<double> rc = F2C(_scStructure.lattice, rf);
  double rshell = aurostd::modulus(rc);

  //
  bool isCountable = false;
  xcomplex<double> iONE(0.0, 1.0);
  phase = 0.0;
  deriv[1].re = 0.0;
  deriv[1].im = 0.0;
  deriv[2].re = 0.0;
  deriv[2].im = 0.0;
  deriv[3].re = 0.0;
  deriv[3].im = 0.0;

  // Count in all atoms which sit on this shell
  if (!_isShellRestricted || (rshell <= _maxShellRadius[centerID] + _AFLOW_APL_EPS_)) {
    isCountable = true;
    int ifound = 0;
    xvector<double> pf(3), pc(3);

    for (_AFLOW_APL_REGISTER_ int ii = -1; ii <= 1; ii++)
      for (_AFLOW_APL_REGISTER_ int jj = -1; jj <= 1; jj++)
        for (_AFLOW_APL_REGISTER_ int kk = -1; kk <= 1; kk++) {
          pf(1) = rf(1) + ii;
          pf(2) = rf(2) + jj;
          pf(3) = rf(3) + kk;
          pc = F2C(_scStructure.lattice, pf);

          if (fabs(aurostd::modulus(pc) - rshell) < _AFLOW_APL_EPS_) {
            //cout << centerID+1 << " " << atomID+1 << " " << aurostd::modulus(pc) << " "; printXVector(pc);
            phase = phase + exp(iONE * scalar_product(qpoint, pc));
            for (uint i = 1; i < 4; i++) {
              deriv[i] = deriv[i] + exp(iONE * scalar_product(qpoint, pc)) * iONE * pc(i);
            }
            ifound++;
          }
        }
    phase *= (1.0 / ifound);
    for (uint i = 1; i < 4; i++) {
      deriv[i] *= (1.0 / ifound);
    }
  }
  //
  return isCountable;
}
// ///////////////////////////////////////////////////////////////////////////
// VECTOR NEEDED FOR ANHARMONIC GRU
void Supercell::calcRl() {
  double cells = _scStructure.atoms.size() / _inStructure.atoms.size();

  xvector<double> Rl;
  Rl(1) = 0;
  Rl(2) = 0;
  Rl(3) = 0;
  for (uint j = 0; j < _scStructure.atoms.size(); j++) {
    Rl_list.push_back(Rl);
    rl_list.push_back(Rl);
  }

  for (uint i = 0; i < (uint)cells; i++) {
    xvector<double> at1ref = _scStructure.atoms[0].cpos;
    xvector<double> at2sc = _scStructure.atoms[i].cpos;
    xmatrix<double> lattice = _scStructure.lattice;

    Rl(1) = 100;
    Rl(2) = 100;
    Rl(3) = 100;
    for (int o = -2; o <= 2; o++) {
      for (int p = -2; p <= 2; p++) {
        for (int q = -2; q <= 2; q++) {
          xvector<double> dummy = (at2sc - at1ref + o * lattice(1) + p * lattice(2) + q * lattice(3));
          if (aurostd::modulus(dummy) < aurostd::modulus(Rl)) {
            Rl = dummy;
          }
        }
      }
    }
    for (uint j = 0; j < _inStructure.atoms.size(); j++) {
      Rl_list[(j * cells) + i] = Rl;
    }
  }

  for (uint i = 0; i < _scStructure.atoms.size(); i++) {
    xvector<double> at2sc = _scStructure.atoms[i].cpos;
    xmatrix<double> lattice = _scStructure.lattice;
    Rl(1) = 100;
    Rl(2) = 100;
    Rl(3) = 100;
    for (int o = -2; o <= 2; o++) {
      for (int p = -2; p <= 2; p++) {
        for (int q = -2; q <= 2; q++) {
          xvector<double> dummy = (at2sc + o * lattice(1) + p * lattice(2) + q * lattice(3));
          if (aurostd::modulus(dummy) < aurostd::modulus(Rl)) {
            Rl = dummy;
          }
        }
      }
    }
    rl_list[i] = Rl;
  }
}
// **** END JJPR *****
// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////

}  // namespace apl
