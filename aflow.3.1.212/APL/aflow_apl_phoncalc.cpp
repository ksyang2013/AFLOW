
#include "aflow_apl.h"

#define ERROR_VERBOSE false  // CO

using namespace std;

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  PhononCalculator::PhononCalculator(Supercell& sc, vector<ClusterSet>& clst,
                                     _xinput& xinput, _aflags& aflags, _kflags& kflags,
				     _xflags& xflags, //_vflags& vflags, 
				     string& AflowIn, Logger& l)
    : _xInput(xinput), _aflowFlags(aflags), _kbinFlags(kflags), _xFlags(xflags),
      _AflowIn(AflowIn), _logger(l), _supercell(sc), _clusters(clst) {
    DISTORTION_MAGNITUDE = 0.015;
    _isGammaEwaldPrecomputed = false;
    DOtar = false;
    xInputsAAPL.clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  PhononCalculator::~PhononCalculator() {
    clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::clear() {
    for (uint i = 0; i < _uniqueDistortions.size(); i++)
      _uniqueDistortions[i].clear();
    _uniqueDistortions.clear();

    for (uint i = 0; i < _uniqueForces.size(); i++) {
      for (uint j = 0; j < _uniqueForces[i].size(); j++)
	_uniqueForces[i][j].clear();
      _uniqueForces[i].clear();
    }
    _uniqueForces.clear();

    for (uint i = 0; i < _forceConstantMatrices.size(); i++)
      _forceConstantMatrices[i].clear();
    _forceConstantMatrices.clear();

    _isGammaEwaldPrecomputed = false;
    _gammaEwaldCorr.clear();

    _bornEffectiveChargeTensor.clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::run(bool aapl_stagebreak) {
    // Check if supercell is already built
    if (!_supercell.isConstructed())
      throw APLRuntimeError("apl::PhononCalculator::run(); The supercell structure has not been initialized yet.");

    // Get all forces required for the construction of force-constant matrices
    calculateForceFields();

    // If AAPL calculations haven't run yet, stop here after
    // APL calculations have been set up
    if (aapl_stagebreak) {
      throw APLStageBreak();
	}

    // For construction of the force-constant matrices we need three
    // independent distortions. Hence, calculate the remaining distortions and
    // forces by the symmetry (if necessary).
    completeForceFields();

    // Ensure that all distortion vectors are along the cartesian directions
    projectToCartesianDirections();

    // Construct the matrix of force-constant matrices for all atoms based
    // on the force fields for the inequivalent atoms
    buildForceConstantMatrices();

    //cout << "NON-SYMMETRIZED:" << std::endl;
    //printForceConstantMatrices(cout);

    // Symmetrization of the force-constant matrices
    symmetrizeForceConstantMatrices();

    //cout << "AFTER SYMMETRIZATION:" << std::endl;
    //printForceConstantMatrices(cout);

    // Force the force-constant matrices to obey the sum-rule conditions
    correctSumRules();

    //printFCShellInfo(cout);

    // Store data into DYNMAT file format - vasp like
    writeDYNMAT();
    // Store data into FORCES file format - phon program
    //writeFORCES();

    //store masses for later uses
    store_masses();  //[PINKU]
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::completeForceFields() {
    //CO - START
    // Test of stupidity...
    if (_supercell.getEPS() == AUROSTD_NAN) {
      throw APLRuntimeError("apl::PhononCalculator::buildForceConstantMatrices(); Need to define symmetry tolerance.");
    }
    //CO - END
    // Show info
    _logger << "Calculating the missing force fields by symmetry." << apl::endl;

    // Let's go
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      // We need to have 3 linearly independent distortions
      if (_uniqueDistortions[i].size() != 3) {
	vector<xvector<double> > allDistortionsOfAtom;
	vector<xvector<double> > testForce;
	vector<vector<xvector<double> > > forcePool;
	xvector<double> testVec(3), testVec0(3);
	const vector<vector<_sym_op> >& agroup = _supercell.getAGROUP();  //CO

	// Generate next independent distortion by symmetry operations...
	uint currentSizeUniqueDistortions = _uniqueDistortions[i].size();
	for (uint idistor = 0; idistor < currentSizeUniqueDistortions; idistor++) {
	  // Apply all symmetry operations and check if it is independent
	  int inequivalentAtomID = _supercell.getUniqueAtomID(i);
	  _supercell.center(inequivalentAtomID);  //CO
	  for (uint symOpID = 0; symOpID < agroup[inequivalentAtomID].size(); symOpID++) {
	    const _sym_op& symOp = agroup[inequivalentAtomID][symOpID];

	    // Transform also all forces
	    //_supercell.center(inequivalentAtomID);  //JAHNATEK ORIGINAL
	    testForce.clear();
	    for (_AFLOW_APL_REGISTER_ int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
	      try {
		_AFLOW_APL_REGISTER_ int l = _supercell.wherefrom(k, inequivalentAtomID, symOp, FALSE);
		testForce.push_back(symOp.Uc * _uniqueForces[i][idistor][l]);
	      } catch (APLLogicError& e) {
		//corey
		//TEMPORARY CODE below by Jahnatek
		//no comments - hard to interpret what's going on
		//WILL INQUIRE SOON, exit if appropriate
		//should not happen if not a derivative structure, exit appropriately
		if (!_supercell.isDerivativeStructure()) {
		  _logger << error << "Mapping problem ? <-> " << k << "." << apl::endl;
		  throw APLLogicError("apl::PhononCalculator::completeForceFields(); Mapping failed.");
		}

		// TEMPORARY CODE =================================================================

#if ERROR_VERBOSE
		cout << "-> PROBLEM: Distortion: " << idistor << "; AtomID " << k << "; CenterAtomID: " << inequivalentAtomID
		     << "; SymOp: " << symOpID << " ( " << symOp.str_type << ") Angle:" << symOp.angle << "; Axis: ";
		printXVector(symOp.axis);
#endif

		// printXVector(_supercell.getSupercellStructure().atoms[k].cpos);
		// xvector<double> zero(3);
		// testForce.push_back( zero);
		testForce.push_back(_uniqueForces[i][idistor][k]);
		int l = 0;
		xvector<double> rotpos;  //CO
		for (; l < (int)_supercell.getNumberOfAtoms(); l++) {
		  uint symOpID2 = 0;
		  for (; symOpID2 < agroup[inequivalentAtomID].size(); symOpID2++) {
		    const _sym_op& symOp2 = agroup[inequivalentAtomID][symOpID2];
		    rotpos = symOp.Uc * inverse(symOp2.Uc) * _supercell.getSupercellStructure().atoms[l].cpos;

		    //if( aurostd::modulus( rotpos - _supercell.getSupercellStructure().atoms[k].cpos ) < _AFLOW_APL_EPS_ ) //JAHNATEK ORIGINAL
		    if (aurostd::modulus(rotpos - _supercell.getSupercellStructure().atoms[k].cpos) < _supercell.getEPS())  //CO
		      {
			break;
		      }
		  }
		  if (symOpID2 != agroup[inequivalentAtomID].size()) {
#if ERROR_VERBOSE
		    cout << l << " " << k << std::endl;
		    //printXMatrix(symOp.Uc * _supercell.getSupercellStructure().agroup[inequivalentAtomID][symOpID2].Uc);
#endif
		    testForce.push_back(symOp.Uc * inverse(agroup[inequivalentAtomID][symOpID2].Uc) * _uniqueForces[i][idistor][l]);
		    break;
		  }
		}
		if (l == (int)_supercell.getNumberOfAtoms()) {
		  throw APLLogicError("Mapping failed.2");
		}

		// TEMPORARY CODE =================================================================
	      }
	    }
	    //_supercell.center(0);  //JAHNATEK ORIGINAL

	    // Get distortion vector (it is in the cartesian form) and apply symmetry rotation
	    testVec = symOp.Uc * _uniqueDistortions[i][idistor];

	    // Ortogonalize new rotated distortion vector on all accepted distortion vectors
	    for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
	      for (_AFLOW_APL_REGISTER_ int l = 0; l < _supercell.getNumberOfAtoms(); l++) {
		testForce[l] = testForce[l] - getModeratedVectorProjection(forcePool[k][l], testVec, allDistortionsOfAtom[k]);
	      }
	      testVec = testVec - getVectorProjection(testVec, allDistortionsOfAtom[k]);
	    }

	    // If the remaining vector is non-zero length, it is new independent direction, hence store it...
	    if (aurostd::modulus(testVec) > _AFLOW_APL_EPS_) {
	      // Normalize to unit length
	      double testVectorLength = aurostd::modulus(testVec);
	      for (int l = 0; l < _supercell.getNumberOfAtoms(); l++) {
		testForce[l] = testForce[l] / testVectorLength;
	      }
	      testVec = testVec / testVectorLength;
	      allDistortionsOfAtom.push_back(testVec);

	      // We suppose the symOpID == 0 is E (Identity) operation, hence the first
	      // independent vector is already the calculated vector, hence new forces are not need to store
	      if (allDistortionsOfAtom.size() > 1) {
		// Store new distortion
		_uniqueDistortions[i].push_back(testVec);
		// Store new force field
		_uniqueForces[i].push_back(testForce);
	      }

	      // Store for next ortogonalization procedure
	      forcePool.push_back(testForce);
	    }
	    if (_uniqueDistortions[i].size() == 3) break;
	  }
	  //_supercell.center_original(); //CO
	  if (_uniqueDistortions[i].size() == 3) break;
	}
	allDistortionsOfAtom.clear();
	for (uint ii = 0; ii < forcePool.size(); ii++) forcePool[ii].clear();
	forcePool.clear();
      }

      // I hope this will never happen...
      if (_uniqueDistortions[i].size() != 3) {
	throw APLRuntimeError("apl::PhononCalculator::completeForceFields(); Can not complete force fields by symmetry.");
      }
    }
    _supercell.center_original();  //CO
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::projectToCartesianDirections() {
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      // Construct transformation matrix A
      xmatrix<double> A(3, 3), U(3, 3);
      for (uint j = 0; j < 3; j++) {
	// Ensure it is unit length
	_uniqueDistortions[i][j] = _uniqueDistortions[i][j] / aurostd::modulus(_uniqueDistortions[i][j]);

	// Copy to rows of U^T matrix
	for (uint k = 1; k <= 3; k++) {
	  U(j + 1, k) = _uniqueDistortions[i][j](k);
	}
      }
      A = inverse(U);

      // Update unique distortion vectors
      _uniqueDistortions[i][0] = trasp(A) * _uniqueDistortions[i][0];
      _uniqueDistortions[i][1] = trasp(A) * _uniqueDistortions[i][1];
      _uniqueDistortions[i][2] = trasp(A) * _uniqueDistortions[i][2];

      // Update forces
      xmatrix<double> m(3, 3);
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
	for (_AFLOW_APL_REGISTER_ int k = 0; k < 3; k++)
	  for (_AFLOW_APL_REGISTER_ int l = 1; l <= 3; l++)
	    m(k + 1, l) = _uniqueForces[i][k][j](l);
	// m = A * m * U; ??? I am not sure...
	m = A * m;
	for (_AFLOW_APL_REGISTER_ int k = 0; k < 3; k++)
	  for (_AFLOW_APL_REGISTER_ int l = 1; l <= 3; l++)
	    _uniqueForces[i][k][j](l) = m(k + 1, l);
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::buildForceConstantMatrices() {
    // Test of stupidity...
    if (!_supercell.getSupercellStructure().fgroup_calculated) {
      throw APLRuntimeError("apl::PhononCalculator::buildForceConstantMatrices(); The factor group has not been calculated yet.");
    }
    //CO - START
    if (_supercell.getEPS() == AUROSTD_NAN) {
      throw APLRuntimeError("apl::PhononCalculator::buildForceConstantMatrices(); Need to define symmetry tolerance.");
    }
    //CO - END

    // Clear old matrices
    for (uint i = 0; i < _forceConstantMatrices.size(); i++)
      _forceConstantMatrices[i].clear();
    _forceConstantMatrices.clear();

    //
    _logger << "Calculating the force constant matrices." << apl::endl;

    // We have a party. Let's fun with us...
    vector<xmatrix<double> > row;
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      // Get the number of this atom in the whole list
      int basedUniqueAtomID = _supercell.getUniqueAtomID(i);

      // This is easy. We know everything. Just construct a set of matrices.
      xmatrix<double> m(3, 3, 1, 1);
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
	for (uint k = 0; k < _uniqueDistortions[i].size(); k++) {
	  double distortionLength = aurostd::modulus(DISTORTION_MAGNITUDE * _uniqueDistortions[i][k]);
	  // FCM element = -F/d, but we will omit minus, because next force trasformations are better
	  // done without it, and in contruction of dyn. matrix we will add it to the sum
	  m(k + 1, 1) = _uniqueForces[i][k][j](1) / distortionLength;
	  m(k + 1, 2) = _uniqueForces[i][k][j](2) / distortionLength;
	  m(k + 1, 3) = _uniqueForces[i][k][j](3) / distortionLength;
	  //cout << i << " " << k << " " << j << " "; printXVector(_uniqueForces[i][k][j]);
	}
	//printXMatrix(m);
	row.push_back(m);
      }
      _forceConstantMatrices.push_back(row);
      row.clear();

      _sym_op symOp;  //CO
      // Calculate rows for next inequivalent atoms starting 1...
      for (int j = 1; j < _supercell.getNumberOfUniqueAtomsOfType(i); j++) {
	try {
	  symOp = _supercell.getSymOpWhichMatchAtoms(_supercell.getUniqueAtomID(i, j), basedUniqueAtomID, _FGROUP_);  //CO
	  //const _sym_op& symOp = _supercell.getSymOpWhichMatchAtoms(_supercell.getUniqueAtomID(i,j),basedUniqueAtomID,_FGROUP_); //JAHNATEK ORIGINAL
	  //cout << basedUniqueAtomID << " -> " << _supercell.getUniqueAtomID(i,j) << " " << symOp.str_type << " shift:"; printXVector(symOp.ftau);
	  //printXVector(_supercell.getSupercellStructure().atoms[basedUniqueAtomID].fpos);
	  //printXVector(_supercell.getSupercellStructure().atoms[_supercell.getUniqueAtomID(i,j)].fpos);
	} catch (APLLogicError& e)  //CO
	  {
	    _logger << error << "Mapping problem " << _supercell.getUniqueAtomID(i, j) << " <-> " << basedUniqueAtomID << "?" << apl::endl;
	    throw APLLogicError("apl::PhononCalculator::buildForceConstantMatrices(); Mapping failed.");
	  }

	for (_AFLOW_APL_REGISTER_ int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
	  try {
	    _AFLOW_APL_REGISTER_ int l = _supercell.where(k, _supercell.getUniqueAtomID(i, j), symOp);
	    //cout << "MAP " << k << " <-> " << l << std::endl;
	    row.push_back(inverse(symOp.Uc) * _forceConstantMatrices[basedUniqueAtomID][l] * symOp.Uc);
	  } catch (APLLogicError& e) {  //CO
	    _logger << error << "Mapping problem " << k << " <-> ?." << apl::endl;
	    throw APLLogicError("apl::PhononCalculator::buildForceConstantMatrices(); Mapping failed.");
	  }
	}
	_forceConstantMatrices.push_back(row);
	row.clear();
      }
      row.clear();
    }

    // Test of correctness
    if ((int)_forceConstantMatrices.size() != _supercell.getNumberOfAtoms()) {
      throw APLRuntimeError("apl::PhononCalculator::buildForceConstantMatrices(); Some problem with the application of factor group operations.");
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::symmetrizeForceConstantMatrices() {
    // Test of stupidity...
    if (!_supercell.getSupercellStructure().agroup_calculated) {
      throw APLRuntimeError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); The site groups have not been calculated yet.");
    }
    //CO - START
    if (_supercell.getEPS() == AUROSTD_NAN) {
      throw APLRuntimeError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); Need to define symmetry tolerance.");
    }
    //CO - END

    //
    _logger << "Symmetrizing the force constant matrices." << apl::endl;

    // Get site symmetry group
    //const vector< vector<_sym_op> >& agroup = _supercell.getSupercellStructure().agroup; //JAHNATEK ORIGINAL
    const vector<vector<_sym_op> >& agroup = _supercell.getAGROUP();  //CO

    //
    vector<xmatrix<double> > row;
    xmatrix<double> m(3, 3);
    m.clear();
    uint agroup_size;  //CO
    for (int i = 0; i < _supercell.getNumberOfAtoms(); i++) {
      if (agroup[i].size() == 0)
	throw APLRuntimeError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); Site point group operations are missing.");

      // Translate the center to this atom
      _supercell.center(i);

      //
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
	agroup_size = agroup[i].size();  //CO
	for (uint symOpID = 0; symOpID < agroup[i].size(); symOpID++) {
	  const _sym_op& symOp = agroup[i][symOpID];

	  try {
	    _AFLOW_APL_REGISTER_ int l = _supercell.where(j, i, symOp, FALSE);
	    //cout << "Mapping " << j << " <-> " << l << std::endl;
	    m = m + (inverse(symOp.Uc) * _forceConstantMatrices[i][l] * symOp.Uc);
	    //CO - START
	  } catch (APLLogicError& e) {
	    //_logger << error << "Mapping problem " << j << " <-> ?. Skipping." << apl::endl;
	    //derivative structures are expected to lose symmetry, don't bother exiting
	    if (!_supercell.isDerivativeStructure()) {
	      _logger << error << "Mapping problem " << j << " <-> ?." << apl::endl;
	      throw APLLogicError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); Mapping failed.");
	    }
	    agroup_size -= 1;  //CO, reduce agroup size
	    //CO - END
	  }
	}
	//CO - START
	//m = ( 1.0 / agroup[i].size() ) * m; //JAHNATEK ORIGINAL
	if (agroup_size) {
	  m = (1.0 / agroup_size) * m;  //CO
	}
	//CO - END
	row.push_back(m);
	m.clear();
      }
      _forceConstantMatrices[i] = row;
      row.clear();

      // Translate the center back
      //_supercell.center_original();  //CO
    }

    // Translate the center back
    //_supercell.center(0); //JAHNATEK ORIGINAL
    _supercell.center_original();  //CO
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::correctSumRules() {
    xmatrix<double> sum(3, 3), sum2(3, 3);

    for (int i = 0; i < _supercell.getNumberOfAtoms(); i++) {
      // Get SUMs
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
	if (i == j) continue;
	//if( !TRUNCATE_DYNAMICAL_MATRIX || ( rshell <= maxShellRadiusOfType[isc1] + _AFLOW_APL_EPS_ ) )
	{
	  sum = sum + _forceConstantMatrices[i][j];
	  sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
	}
      }

      //cout << "SUM RULE 1:" << std::endl;
      //printXMatrix(_forceConstantMatrices[i][i]);
      //printXMatrix(-1.0*sum);
      //cout << "SUM RULE 2:" << std::endl;
      //printXMatrix(sum);
      //printXMatrix(sum2);

      // Correct SUM2
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
	if (i == j) continue;
	_forceConstantMatrices[i][j] = 0.5 * (_forceConstantMatrices[i][j] + trasp(_forceConstantMatrices[j][i]));
	_forceConstantMatrices[j][i] = trasp(_forceConstantMatrices[i][j]);
      }

      // Get SUMs again
      sum.clear();
      sum2.clear();
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
	if (i == j) continue;
	//if( !TRUNCATE_DYNAMICAL_MATRIX || ( rshell <= maxShellRadiusOfType[isc1] + _AFLOW_APL_EPS_ ) )
	{
	  sum = sum + _forceConstantMatrices[i][j];
	  sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
	}
      }

      // Correct SUM1 to satisfied
      _forceConstantMatrices[i][i] = -sum;

      //cout << "cSUM RULE 1:" << std::endl;
      //printXMatrix(_forceConstantMatrices[i][i]);
      //printXMatrix(-1.0*sum);
      //cout << "cSUM RULE 2:" << std::endl;
      //printXMatrix(sum);
      //printXMatrix(sum2);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::printForceConstantMatrices(ostream& os) {
    int units = 1;
    double conversionFactor = 1.0;

    switch (units) {
    case (1):
      os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
      conversionFactor = 1.0;
      break;
    case (2):
      os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
      conversionFactor = 1602.17733;
      break;
    case (3):
      os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
      conversionFactor = 16.0217733;
      break;
    }
    os << std::endl;

    for (int i = 0; i < _supercell.getNumberOfAtoms(); i++)
      //for(int ii = 0; ii < _supercell.getNumberOfUniqueAtoms(); ii++)
      {
	//int i = _supercell.getUniqueAtomID(ii);
	for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
	  os << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	  os << setprecision(4);
	  os << "- MATRIX: " << i + 1 << "/" << k + 1 << " " << k + 1 << "/" << i + 1 << std::endl;
	  for (int m = 1; m <= 3; m++) {
	    for (int n = 1; n <= 3; n++)
	      os << setw(10) << (conversionFactor * _forceConstantMatrices[k][i](m, n)) << " ";
	    os << " ";
	    for (int n = 1; n <= 3; n++)
	      os << setw(10) << (conversionFactor * _forceConstantMatrices[i][k](n, m)) << " ";
	    os << std::endl;
	  }
	  os << std::endl;
	}
      }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::printFCShellInfo(ostream& os) {
    int units = 4;
    double conversionFactor = 1.0;

    switch (units) {
    case (1):
      os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
      conversionFactor = 1.0;
      break;
    case (2):
      os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
      conversionFactor = 1602.17733;
      break;
    case (3):
      os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
      conversionFactor = 16.0217733;
      break;
    case (4):
      os << "FORCE CONSTANT MATRICES in 10^3 Dyn/cm:" << std::endl;
      conversionFactor = 16.0217733;
      break;
    }
    os << std::endl;

    int maxshell = _supercell.getMaxShellID();
    if (maxshell == -1) maxshell = 25;
    std::vector<ShellHandle> sh;
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      ShellHandle s;
      sh.push_back(s);
      sh.back().init(_supercell.getInputStructure(),
		     _supercell.getInputStructure().iatoms[i][0],
		     maxshell);
      sh[i].splitBySymmetry();
      sh[i].mapStructure(_supercell.getSupercellStructure(), _supercell.getUniqueAtomID(i));
    }

    //
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      sh[i].printReport(cout);
      for (int ishell = 0; ishell <= sh[i].getLastOccupiedShell(); ishell++) {
	for (int isubshell = 0; isubshell < sh[i].getNumberOfSubshells(ishell); isubshell++) {
	  const deque<_atom>& atomsAtSameShell = sh[i].getAtomsAtSameShell(ishell, isubshell);
	  cout << "SHELL " << ishell << " " << isubshell << std::endl;

	  for (uint ai = 0; ai < atomsAtSameShell.size(); ai++) {
	    int nb = atomsAtSameShell[ai].number;
	    cout << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	    cout << setprecision(4);
	    cout << "- MATRIX: " << i << "/" << nb << " " << nb << "/" << i << std::endl;
	    //cout << "atom " << setw(3) << nb << ": "; printXVector(atomsAtSameShell[ai].cpos);
	    //cout << "atom " << setw(3) << i << ": "; printXVector(_supercell.getSupercellStructure().atoms[_supercell.getUniqueAtomID(i)].cpos);
	    for (int m = 1; m <= 3; m++) {
	      for (int n = 1; n <= 3; n++)
		cout << setw(10) << (conversionFactor * _forceConstantMatrices[nb][i](m, n)) << " ";
	      cout << " ";
	      for (int n = 1; n <= 3; n++)
		cout << setw(10) << (conversionFactor * _forceConstantMatrices[i][nb](n, m)) << " ";
	      cout << std::endl;
	    }
	    cout << std::endl;
	  }
	}
      }
    }

    // Clear
    for (uint i = 0; i < sh.size(); i++)
      sh[i].clear();
    sh.clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  // Y. Wang et.al, J. Phys.:Condens. Matter 22, 202201 (2010)
  // DOI: 10.1088/0953-8984/22/20/202201

  // ME180827 - Overloaded to calculate derivative for AAPL
  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getNonanalyticalTermWang(_q, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q,
                                                                        vector<xmatrix<xcomplex<double> > >& derivative,
                                                                        bool calc_derivative) {
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructureLight();  //CO

    // to correct the q=\Gamma as a limit
    xvector<double> q(_q);
    if (aurostd::modulus(q) < _AFLOW_APL_EPS_) {
      q(1) = _AFLOW_APL_EPS_ * 1.001;
    }

    uint pcAtomsSize = pc.atoms.size();
    uint nBranches = 3 * pcAtomsSize;

    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches);

    if (calc_derivative) {  // reset derivative
      derivative.clear();
      xmatrix<xcomplex<double> > mat(nBranches, nBranches, 1, 1);
      derivative.assign(3, mat);
    }

    // Calculation
    double fac0 = 13.605826 * 2.0 * 0.529177249;  // from a.u. to eV/A
    double volume = det(pc.lattice);
    double fac1 = 4.0 * PI / volume;
    double nbCells = det(sc.lattice) / volume;

    if (aurostd::modulus(q) > _AFLOW_APL_EPS_) {
      double dotprod = scalar_product(q, _dielectricTensor * q);
      double prefactor = fac0 * fac1/(dotprod * nbCells);
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
	for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
	  for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
	    for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
	      int typei = pc.atoms[ipc1].type;
	      int typej = pc.atoms[ipc2].type;
              double borni = (q * _bornEffectiveChargeTensor[typei])(ix);
              double bornj = (q * _bornEffectiveChargeTensor[typej])(iy);
	      dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * borni * bornj;
              if (calc_derivative) {
                for (int d = 0; d < 3; d++) {
                  xcomplex<double> coeff(0, 0);
                  coeff += borni * _bornEffectiveChargeTensor[typej](iy, d + 1);
                  coeff += bornj * _bornEffectiveChargeTensor[typei](ix, d + 1);
                  coeff -= 2 * borni * bornj * scalar_product(_dielectricTensor(d + 1), q)/dotprod;
                  derivative[d](3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * coeff;
                }
              }
	    }
          }
	    }
	}
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // X. Gonze et al., Phys. Rev. B 50, 13035 (1994)
  // X. Gonze and Ch. Lee, Phys. Rev. B 55, 10355 (1997)

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermGonze(const xvector<double> kpoint) {
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();
    //    uint nBranches = 3 * pcAtomsSize; // not needed

    if (!_isGammaEwaldPrecomputed) {
      xvector<double> zero(3);
      xmatrix<xcomplex<double> > dynamicalMatrix0(getEwaldSumDipolDipolContribution(zero, false));

      _gammaEwaldCorr.clear();
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
	xmatrix<xcomplex<double> > sum(3, 3);
	for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
	  for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
	    for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
	      sum(ix, iy) += dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy);
	}
	_gammaEwaldCorr.push_back(sum);
      }

      _isGammaEwaldPrecomputed = true;
    }

    //
    xmatrix<xcomplex<double> > dynamicalMatrix(getEwaldSumDipolDipolContribution(kpoint));

    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
	for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
	  dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= _gammaEwaldCorr[ipc1](ix, iy);
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  xmatrix<xcomplex<double> > PhononCalculator::getEwaldSumDipolDipolContribution(const xvector<double> qpoint, bool includeTerm1) {
    // Definitions
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructureLight();  //CO

    uint pcAtomsSize = pc.atoms.size();
    uint nBranches = 3 * pcAtomsSize;

    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches);

    double gmax = 14.0;
    double lambda = 1.0;
    double lambda2 = lambda * lambda;
    double lambda3 = lambda2 * lambda;
    double geg = gmax * lambda2 * 4.0;

    // Reciprocal Space
    xmatrix<double> klattice = trasp(ReciprocalLattice(pc.lattice));

    // Grid
    int n1 = (int)(sqrt(geg) / aurostd::modulus(klattice(1))) + 1;
    int n2 = (int)(sqrt(geg) / aurostd::modulus(klattice(2))) + 1;
    int n3 = (int)(sqrt(geg) / aurostd::modulus(klattice(3))) + 1;

    // Calculation
    double fac0 = 13.605826 * 2.0 * 0.529177249;  // from a.u. to eV/A
    double SQRTPI = sqrt(PI);
    double volume = det(pc.lattice);
    double fac = 4.0 * PI / volume;
    xcomplex<double> iONE(0.0, 1.0);

    // Term 1 - Reciprocal space sum

    if (includeTerm1) {
      for (_AFLOW_APL_REGISTER_ int m1 = -n1; m1 <= n1; m1++) {
	for (_AFLOW_APL_REGISTER_ int m2 = -n2; m2 <= n2; m2++) {
	  for (_AFLOW_APL_REGISTER_ int m3 = -n3; m3 <= n3; m3++) {
	    xvector<double> g = m1 * klattice(1) + m2 * klattice(2) + m3 * klattice(3) + qpoint;

	    geg = scalar_product(g, _dielectricTensor * g);

	    if (fabs(geg) > _AFLOW_APL_EPS_ && geg / lambda2 / 4.0 < gmax) {
	      double fac2 = fac * exp(-geg / lambda2 / 4.0) / geg;

	      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
		xvector<double> zag = g * _bornEffectiveChargeTensor[pc.atoms[ipc1].type];

		for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
		  xvector<double> zbg = g * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];

		  //xcomplex<double> e;
		  //(void)_supercell.calcShellPhaseFactor(ipc2,ipc1,g,e);
		  //xcomplex<double> facg = fac2 * e;
		  xcomplex<double> facg = fac2 * exp(iONE * scalar_product(g, sc.atoms[ipc2].cpos - sc.atoms[ipc1].cpos));

		  for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
		    for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
		      dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += fac0 * facg * zag(ix) * zbg(iy);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // Term 2 - Real space sum
    /*
      for(_AFLOW_APL_REGISTER_ int m1 = -n1; m1 <= n1; m1++)
      for(_AFLOW_APL_REGISTER_ int m2 = -n2; m2 <= n2; m2++)
      for(_AFLOW_APL_REGISTER_ int m3 = -n2; m3 <= n3; m3++) {
      xvector<double> rc = m1 * pc.lattice(1) + m2 * pc.lattice(2)
      + m3 * pc.lattice(3);

      //xvector<double> zero(3);
      //xvector<double> rf = _supercell.getFPositionItsNearestImage(rc,zero,pc.lattice);
      //rc = F2C(pc.lattice,rf);

      if( aurostd::modulus(rc) < _AFLOW_APL_EPS_ ) continue;

      //
      xvector<double> delta = _inverseDielectricTensor * rc;
      double D = sqrt( scalar_product(delta,rc) );

      //
      xmatrix<double> H(3,3);
      xvector<double> x = lambda * delta;
      double y = lambda * D;
      double y2 = y * y;
      double ym2 = 1.0 / y2;
      double emy2dpi = 2.0 * exp( -y2 ) / SQRTPI;
      double erfcdy = erfc(y) / y;
      double c1 = ym2 * ( 3.0 * erfcdy * ym2 + ( emy2dpi * ( 3.0 * ym2 + 2.0 ) ) );
      double c2 = ym2 * ( erfcdy + emy2dpi );
      for(_AFLOW_APL_REGISTER_ int a = 1; a <= 3; a++)
      for(_AFLOW_APL_REGISTER_ int b = 1; b <= 3; b++) {
      H(a,b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a,b) * c2;
      }

      //
      xcomplex<double> e = exp( iONE * scalar_product(qpoint,rc) );
      xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;

      //
      for(uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      xmatrix<double> zh = _bornEffectiveChargeTensor[pc.atoms[ipc1].type] * H;

      for(uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
      xmatrix<double> zhz = zh * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
 
      for(_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
      for(_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
      dynamicalMatrix(3*ipc1+ix,3*ipc2+iy) -= fac * zhz(ix,iy);
      }
      }
      }
    */

    // Term 2
    uint scAtomsSize = sc.atoms.size();
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
	uint ipc2 = _supercell.sc2pcMap(isc2);

	xvector<double> rf = _supercell.getFPositionItsNearestImage(isc2, isc1);
	xvector<double> rc = F2C(sc.lattice, rf);

	if (aurostd::modulus(rc) < _AFLOW_APL_EPS_) continue;

	//
	xvector<double> delta = _inverseDielectricTensor * rc;
	double D = sqrt(scalar_product(delta, rc));

	//
	xmatrix<double> H(3, 3);
	xvector<double> x = lambda * delta;
	double y = lambda * D;
	double y2 = y * y;
	double ym2 = 1.0 / y2;
	double emy2dpi = 2.0 * exp(-y2) / SQRTPI;
	double erfcdy = erfc(y) / y;
	double c1 = ym2 * (3.0 * erfcdy * ym2 + (emy2dpi * (3.0 * ym2 + 2.0)));
	double c2 = ym2 * (erfcdy + emy2dpi);
	for (_AFLOW_APL_REGISTER_ int a = 1; a <= 3; a++) {
	  for (_AFLOW_APL_REGISTER_ int b = 1; b <= 3; b++) {
	    H(a, b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a, b) * c2;
	  }
	}

	xmatrix<double> za = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
	xmatrix<double> zb = _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
	xmatrix<double> zhz = za * H * zb;

	//
	xcomplex<double> e;  // = exp( iONE * scalar_product(qpoint,rc) );
	(void)_supercell.calcShellPhaseFactor(isc2, isc1, qpoint, e);

	//
	xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;
	for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
	  for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
	    dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) -= fac * zhz(ix, iy);
      }
    }

    // Term 3 - Limiting contribution

    double facterm3 = fac0 * 4.0 * lambda3 * _recsqrtDielectricTensorDeterminant / (3.0 * SQRTPI);
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      xmatrix<double> z = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
      xmatrix<double> zez = z * _inverseDielectricTensor * z;

      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
	for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
	  dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= facterm3 * zez(ix, iy);
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////
  // ME180827 - Overloaded to calculate derivative for AAPL
  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getDynamicalMatrix(kpoint, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint,
                                                                vector<xmatrix<xcomplex<double> > >& dDynMat,
                                                                bool calc_derivative) {
    uint scAtomsSize = _supercell.getSupercellStructure().atoms.size();
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();

    uint nBranches = 3 * pcAtomsSize;
    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches, 1, 1);
    xmatrix<xcomplex<double> > dynamicalMatrix0(nBranches, nBranches, 1, 1);

    xcomplex<double> phase;
    double value;
    // ME 180828 - Prepare derivative calculation
    xvector<xcomplex<double> > derivative;
    double nbCells = 1.0;  // for NAC derivative
    for (int i = 1; i < 4; i++) {
      nbCells *= _supercell.scell[i];
    }
    vector<xmatrix<xcomplex<double> > > dDynMat_NAC;
    if (calc_derivative) {  // reset dDynMat
      dDynMat.clear();
      xmatrix<xcomplex<double> > mat(nBranches, nBranches, 1, 1);
      dDynMat.assign(3, mat);
    }

    // Calculate nonanalytical contribution
    xmatrix<xcomplex<double> > dynamicalMatrixNA(nBranches, nBranches, 1, 1);
    if (_isPolarMaterial)
      dynamicalMatrixNA = getNonanalyticalTermWang(kpoint, dDynMat_NAC, calc_derivative);

    // Loop over primitive cell
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
	uint ipc2 = _supercell.sc2pcMap(isc2);
        int neq;  // Important for NAC derivative
	if (_supercell.calcShellPhaseFactor(isc2, isc1, kpoint, phase,
                                            neq, derivative, calc_derivative)) {  // ME180827
	  for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
	    for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
	      value = 0.5 * (_forceConstantMatrices[isc1][isc2](ix, iy) +
			     _forceConstantMatrices[isc2][isc1](iy, ix));
	      dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) -= value * phase;
	      if (_isPolarMaterial)
		dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
	      dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) -= value;
	      if (_isPolarMaterial)
		dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy);
              if (calc_derivative) {
                for (int d = 0; d < 3; d++) {
                  dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) -= value * derivative[d+1];
                  if (_isPolarMaterial && (aurostd::modulus(kpoint) > _AFLOW_APL_EPS_)) {
                    xcomplex<double> nac = ((double) neq) * phase * dDynMat_NAC[d](3 * ipc1 + ix, 3 * ipc2 + iy);
                    dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) += nac;
	    }
	  }
	}
      }
    }
	}
      }
    }
    //printXMatrix2(dynamicalMatrix);

    // Subtract the sum of all "forces" from the central atom, this is like an automatic sum rule...
    for (uint i = 0; i < pcAtomsSize; i++) {
      for (uint j = 0; j < pcAtomsSize; j++) {
	for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
	  for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
	    dynamicalMatrix(3 * i + ix, 3 * i + iy) = dynamicalMatrix(3 * i + ix, 3 * i + iy) - dynamicalMatrix0(3 * i + ix, 3 * j + iy);
	  }
	}
      }
    }

    // Get correction for polar materials
    //if( _isPolarMaterial )
    // dynamicMatrix += getNonanalyticalTermGonze(kpoint);

    // Make it hermitian
    for (uint i = 0; i <= pcAtomsSize - 1; i++) {
      for (uint j = 0; j <= i; j++) {
	for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
	  for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
	    dynamicalMatrix(3 * i + ix, 3 * j + iy) += conj(dynamicalMatrix(3 * j + iy, 3 * i + ix));
	    dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 0.5;
	    dynamicalMatrix(3 * j + iy, 3 * i + ix) = conj(dynamicalMatrix(3 * i + ix, 3 * j + iy));
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) += conj(dDynMat[d](3 * j + iy, 3 * i + ix));
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 0.5;
                dDynMat[d](3 * j + iy, 3 * i + ix) = conj(dDynMat[d](3 * i + ix, 3 * j + iy));
              }
            }
	  }
	}
      }
    }

    // Divide by masses
    for (uint i = 0; i < pcAtomsSize; i++) {
      double mass_i = _supercell.getAtomMass(_supercell.pc2scMap(i));
      for (uint j = 0; j < pcAtomsSize; j++) {
	double mass_j = _supercell.getAtomMass(_supercell.pc2scMap(j));
	for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
	  for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
	    dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 1.0/sqrt(mass_i * mass_j);
              }
            }
	  }
	}
      }
    }

    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME180827 - Overloaded to calculate derivative and eigenvectors for AAPL
  xvector<double> PhononCalculator::getEigenvalues(const xvector<double>& kpoint) {
    const xstructure& pc = _supercell.getInputStructureLight();  //CO
    uint nBranches = 3 * pc.atoms.size();
    xmatrix<xcomplex<double> > placeholder_eigen(nBranches, nBranches, 1, 1);
    vector<xmatrix<xcomplex<double> > > placeholder_mat;
    return getEigenvalues(kpoint, placeholder_eigen, placeholder_mat, false);
  }

  xvector<double> PhononCalculator::getEigenvalues(const xvector<double>& kpoint,
                                                   xmatrix<xcomplex<double> >& eigenvectors,
                                                   vector<xmatrix<xcomplex<double> > >& dDynMat,
                                                   bool calc_derivative) {
    // Get dynamic matrix
    xmatrix<xcomplex<double> > dynamicalMatrix = getDynamicalMatrix(kpoint, dDynMat, calc_derivative);

    // Diagonalize
    xvector<double> eigenvalues(dynamicalMatrix.rows, 1);
//    xmatrix<xcomplex<double> > unitaryMatrix;  OBSOLETE ME 180827


/*
//OBSOLETE ME 180828 - zheevByJacobiRotation may not yield the correct eigenvectors.
#ifdef USE_MKL
    zheevMKL(dynamicalMatrix, eigenvalues, eigenvectors);
#else
    //tred2(dynamicalMatrix);
    zheevByJacobiRotation(dynamicalMatrix, eigenvalues, eigenvectors);
#endif
*/

// ME 180828 
    apl::aplEigensystems e;
    e.eigen_calculation(dynamicalMatrix, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);

    //
    return eigenvalues;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME180827 - Overloaded to calculate derivative and eigenvectors for AAPL
  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const IPCFreqFlags& flags) {
    const xstructure& pc = _supercell.getInputStructureLight();  //CO
    uint nBranches = 3 * pc.atoms.size();
    xmatrix<xcomplex<double> > placeholder_eigen(nBranches, nBranches, 1, 1);
    vector<xmatrix<xcomplex<double> > > placeholder_mat;
    return getFrequency(kpoint, flags, placeholder_eigen, placeholder_mat, false);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, IPCFreqFlags flags,
                                                 xmatrix<xcomplex<double> >& eigenvectors,
                                                 vector<xmatrix<xcomplex<double> > >& dDynMat,
                                                 bool calc_derivative) {
    // Compute frequency(omega) from eigenvalues [in eV/A/A/atomic_mass_unit]
    xvector<double> omega = getEigenvalues(kpoint, eigenvectors, dDynMat, calc_derivative);

    // Get value of conversion factor
    double conversionFactor = getFrequencyConversionFactor(apl::RAW | apl::OMEGA, flags);

    // Transform values to desired format
    for (_AFLOW_APL_REGISTER_ int i = omega.lrows; i <= omega.urows; i++) {
      if (omega(i) < 0) {
	if (flags & ALLOW_NEGATIVE)
	  omega(i) = -sqrt(-omega(i));
	else
	  omega(i) = 0.0;
      } else {
	omega(i) = sqrt(omega(i));
      }

      // Convert to desired units
      omega(i) *= conversionFactor;
    }

    // Return
    return (omega);
  }

  // ///////////////////////////////////////////////////////////////////////////

  double PhononCalculator::getFrequencyConversionFactor(IPCFreqFlags inFlags, IPCFreqFlags outFlags) {
    double conversionFactor = 1.0;

    // Conversion from eV/A/A/atomic_mass_unit -> something
    if (inFlags & apl::RAW) {
      if (outFlags & apl::RAW) {
	// Transform to eV/A/A/atomic_mass_unit
	conversionFactor = 1.0;
      } else if (outFlags & apl::HERTZ) {
	// Transform to s-1; sqrt(EV_TO_JOULE / (ANGSTROM_TO_METER*ANGSTROM_TO_METER) / AMU_TO_KG);
	conversionFactor = 0.98226977255434387350E14;
      } else if (outFlags & apl::THZ) {
	// Transform to THz; (in Hertz) / 1E12;
	conversionFactor = 98.226977255434387350;
      } else if (outFlags & apl::RECIPROCAL_CM) {
	// Transform to cm-1; 1/lambda(m) = freq.(s-1) / light_speed(m/s);
	conversionFactor = 1E-2 * 0.98226977255434387350E14 / 2.99792458E8;
      } else if (outFlags & apl::MEV) {
	// Transform to meV; E(eV) = h(eV.s) * freq(s-1); h[(from J.s) -> (eV.s)] = 4.1356673310E-15
	conversionFactor = 0.98226977255434387350E14 * 4.1356673310E-15 / 1E-3;
      } else {
	throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
      }
    }

    // Conversion from THz -> something
    else if (inFlags & apl::THZ) {
      if (outFlags & apl::RAW) {
	// Transform to eV/A/A/atomic_mass_unit
	conversionFactor = 1.0 / 98.226977255434387350;
      } else if (outFlags & apl::THZ) {
	conversionFactor = 1.0;
      } else if (outFlags & apl::MEV) {
	conversionFactor = 4.1356673310;
      } else {
	throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
      }
    }

    // Nothing suits?
    else
      throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");

    //
    if ((outFlags & OMEGA) && !(inFlags & OMEGA))
      conversionFactor *= 2.0 * M_PI;
    if (!(outFlags & OMEGA) && (inFlags & OMEGA))
      conversionFactor /= 2.0 * M_PI;

    //
    return (conversionFactor);
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::writeOUTPUT(_xinput& xinput) { //CO 180409
    if(!( xinput.AFLOW_MODE_VASP || xinput.AFLOW_MODE_AIMS )) { 
      throw APLRuntimeError("apl::PhononCalculator:writeOUTPUT(); Input -> aflow.in conversion unknown.");
    }

    //copying from createAFLOWIN
    _xvasp xvasp(xinput.xvasp);
    _vflags vflags(_xFlags.vflags);
    _xaims xaims(xinput.xaims);
    _aimsflags aimsflags(_xFlags.aimsflags);

    string directory=xinput.getDirectory();
    if(directory.empty()){throw APLRuntimeError("apl::PhononCalculator:writeOUTPUT(); no output directory found");}

    if(!aurostd::FileExist(directory)){aurostd::DirectoryMake(directory);}  // Create directory if it is not created
    aurostd::DirectoryChmod("777", directory);                              // CHMOD Directory 777

    stringstream outfile;

    // OK, fill it...
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[AFLOW] _ ___ _" << std::endl;
    outfile << "[AFLOW] / \\| o \\ |" << std::endl;
    outfile << "[AFLOW] | o | _/ |_" << std::endl;
    outfile << "[AFLOW] |_n_|_| |___| automatic generated file" << std::endl;
    outfile << "[AFLOW]" << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if(xinput.AFLOW_MODE_VASP){outfile << "[AFLOW_MODE=VASP]" << std::endl;}
    if(xinput.AFLOW_MODE_AIMS){outfile << "[AFLOW_MODE=AIMS]" << std::endl;}
    if(!_kbinFlags.KZIP_BIN.empty()){outfile << "[AFLOW_MODE_ZIP=" << _kbinFlags.KZIP_BIN << "]" << std::endl;}  //CO

    //CO 180130 - START
    //corey - at some point, fix alien mode for aims, for now omit!
    if(xinput.AFLOW_MODE_VASP){
      //adding aflow.rc stuff
      outfile << "[AFLOW_MODE_BINARY=";
      if(!_kbinFlags.KBIN_BIN.empty()){outfile << _kbinFlags.KBIN_BIN;}
      else{outfile << DEFAULT_VASP_BIN;}
      outfile << "]" << std::endl;
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      if(!(_kbinFlags.KBIN_MPI || XHOST.MPI)){outfile << "#";}
      outfile << "[AFLOW_MODE_MPI]" << std::endl;
      //be super cautious and avoid empty tags here
      string NCPUS_VAL; get_NCPUS(NCPUS_VAL); //CO 180214
      outfile << "[AFLOW_MODE_MPI_MODE]NCPUS=" << NCPUS_VAL << " " << std::endl;
      outfile << "[AFLOW_MODE_MPI_MODE]COMMAND =\"" << MPI_COMMAND_DEFAULT << "\" " << std::endl;
      if( _kbinFlags.KBIN_MPI_AUTOTUNE ) {outfile << "[AFLOW_MODE_MPI_MODE]AUTOTUNE " << std::endl;}
      outfile << "[AFLOW_MODE_MPI_MODE]BINARY=\"";
      if(!_kbinFlags.KBIN_MPI_BIN.empty()){outfile << _kbinFlags.KBIN_MPI_BIN;}
      else{outfile << DEFAULT_VASP_MPI_BIN;}
      outfile << "\"" << std::endl;
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      //CO 180130 - STOP

      // INCAR
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[VASP_FORCE_OPTION]WAVECAR=OFF" << std::endl;
      outfile << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl;
      outfile << "[VASP_FORCE_OPTION]PREC=PHONONS" << std::endl;  // Modified JJPR
      outfile << "[VASP_FORCE_OPTION]ALGO=NORMAL" << std::endl;

      if (vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) outfile << "[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << std::endl;
      if (vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry) outfile << "[VASP_FORCE_OPTION]ABMIX=" << vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme << std::endl;
      if (vflags.KBIN_VASP_FORCE_OPTION_TYPE.isentry) outfile << "[VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << std::endl;
      //PINKU LDAU OPTION
      if (_check_LDAU2_ON != "") {
	outfile << AFLOWIN_SEPARATION_LINE << std::endl;
	outfile << _check_LDAU2_ON << std::endl;
	outfile << _LDAU_PARAMETERS << std::endl;
	outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      }
      //PINKU LDAU OPTION
      if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry){outfile << "[VASP_FORCE_OPTION]AUTO_MAGMOM=" << (vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option ? "ON" : "OFF") << std::endl;}
    
      if (vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && vflags.KBIN_VASP_FORCE_OPTION_BADER.option){outfile << "[VASP_FORCE_OPTION]BADER=ON" << std::endl;}
      else{outfile << "[VASP_FORCE_OPTION]BADER=OFF" << std::endl;}
    
      if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && vflags.KBIN_VASP_FORCE_OPTION_ELF.option){outfile << "[VASP_FORCE_OPTION]ELF=ON" << std::endl;}
      else{outfile << "[VASP_FORCE_OPTION]ELF=OFF" << std::endl;}

      if (vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) {
	if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.option){outfile << "[VASP_FORCE_OPTION]SPIN=ON" << std::endl;}
	else{outfile << "[VASP_FORCE_OPTION]SPIN=OFF" << std::endl;}
      }

      if(vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.flag("EVEN") || vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("EVEN")){outfile << "[VASP_FORCE_OPTION]KPOINTS=EVEN" << std::endl;}
      if(vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.flag("ODD") || vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("ODD")){outfile << "[VASP_FORCE_OPTION]KPOINTS=ODD" << std::endl;}

      if(vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.isentry){outfile << "[VASP_FORCE_OPTION]IGNORE_AFIX=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.content_string << std::endl;}
      else{outfile << "[VASP_FORCE_OPTION]IGNORE_AFIX=NPARC" << std::endl;}

      // [VASP_FORCE_OPTION]KPOINTS=KEEPK
      // [VASP_FORCE_OPTION]NBANDS
      // outfile << "[VASP_INCAR_MODE_EXPLICIT]" << std::endl;
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[VASP_INCAR_MODE_EXPLICIT]START" << std::endl;
      outfile << xvasp.INCAR.str();
      if (_PSTRESS != "") outfile << _PSTRESS << std::endl;//PINKU PSTRESS OPTION
      outfile << "[VASP_INCAR_MODE_EXPLICIT]STOP" << std::endl;

      // KPOINTS
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[VASP_KPOINTS_MODE_IMPLICIT] " << std::endl;
      outfile << "[VASP_KPOINTS_FILE]KSCHEME=" << xvasp.AVASP_KSCHEME << " " << std::endl;
      outfile << "[VASP_KPOINTS_FILE]KPPRA=" << xvasp.AVASP_value_KPPRA << std::endl;
      outfile << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << xvasp.AVASP_KSCHEME << " " << std::endl;
      outfile << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << xvasp.AVASP_value_KPPRA << std::endl;

      // POTCAR
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[VASP_POTCAR_MODE_IMPLICIT] " << std::endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry){
	for(uint i = 0; i < xvasp.str.species.size(); i++){outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species.at(i) << std::endl;}
      } else {
	for(uint i = 0; i < xvasp.str.species_pp.size(); i++){outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species_pp.at(i) << std::endl;}
      }

      // POSCAR
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
      xvasp.str.is_vasp4_poscar_format = TRUE;
      xvasp.str.is_vasp5_poscar_format = FALSE;
      outfile << xvasp.str;
      outfile << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;

      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[VASP_RUN]STATIC" << std::endl;
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;

      if (xvasp.aopts.flag("AFLOWIN_FLAG::QE")) { //(AFLOWIN_QE_FLAG) {
	outfile << AFLOWIN_SEPARATION_LINE << std::endl;
	outfile << "[QE_GEOM_MODE_EXPLICIT]START " << std::endl;
	xstructure qestr(xvasp.str);
	qestr.xstructure2qe();
	outfile << qestr;
	outfile << "[QE_GEOM_MODE_EXPLICIT]STOP " << std::endl;
      }
      if (xvasp.aopts.flag("AFLOWIN_FLAG::ABINIT")) {
	outfile << AFLOWIN_SEPARATION_LINE << std::endl;
	outfile << "[ABINIT_GEOM_MODE_EXPLICIT]START " << std::endl;
	xstructure abinitstr(xvasp.str);
	abinitstr.xstructure2abinit();
	outfile << abinitstr;
	outfile << "[ABINIT_GEOM_MODE_EXPLICIT]STOP " << std::endl;
      }
      if (xvasp.aopts.flag("AFLOWIN_FLAG::AIMS")) {
	outfile << AFLOWIN_SEPARATION_LINE << std::endl;
	outfile << "[AIMS_GEOM_MODE_EXPLICIT]START " << std::endl;
	xstructure aimsstr(xvasp.str);
	aimsstr.xstructure2aims();
	outfile << aimsstr;
	outfile << "[AIMS_GEOM_MODE_EXPLICIT]STOP " << std::endl;
      }
    }
    if(xinput.AFLOW_MODE_AIMS){
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[AIMS_CONTROL_MODE_EXPLICIT]START " << std::endl;
      outfile << xaims.CONTROL.str();
      outfile << "[AIMS_CONTROL_MODE_EXPLICIT]STOP " << std::endl;
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      outfile << "[AIMS_GEOM_MODE_EXPLICIT]START " << std::endl;
      outfile << xaims.str;
      outfile << "[AIMS_GEOM_MODE_EXPLICIT]STOP " << std::endl;
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
  
      //also write out
      if(1){
	KBIN::AIMS_Write_CONTROL(xaims,aimsflags);
	xaims.GEOM.clear(); xaims.GEOM.str("");
	xaims.GEOM << xaims.str;

	string geom_filename = xaims.Directory + "/" + AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM;
	aurostd::stringstream2file(xaims.GEOM, geom_filename);
	if(!aurostd::FileExist(geom_filename)){
	  throw apl::APLRuntimeError("apl::PhononCalculator::createAIMSOUTPUT(); Cannot create [" + AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM + "] file.");
	}
	aurostd::ChmodFile("a+rw", geom_filename);
      }
    }

    //CO - START
    string filename = directory + string("/") + _AFLOWIN_;
    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)){throw apl::APLRuntimeError("apl::PhononCalculator::createAFLOWIN(); Cannot create [" + _AFLOWIN_ + "] file.");}
    aurostd::ChmodFile("a+rw", filename); // CHMOD a+rw _AFLOWIN_
    //CO - END
  }

  //////////////////////////////////////////////////////////////////////////////

  void PhononCalculator::writeDYNMAT() {
    _logger << "Writing forces into file " << DEFAULT_APL_DYNMAT_FILE << "." << apl::endl;

    //
    //CO - START
    //ofstream outfile("DYNMAT", ios_base::out);
    //if (!outfile.is_open())
    //  throw apl::APLRuntimeError("PhononCalculator::writeDYNMAT(); Cannot open output file.");
    stringstream outfile;
    //CO - END

    // 1st line
    outfile << _supercell.getNumberOfUniqueAtoms() << " ";
    outfile << _supercell.getNumberOfAtoms() << " ";
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++)
      dof += _uniqueDistortions[i].size();
    outfile << dof << std::endl;

    // 2nd line
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(3);
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      if (i != 0) outfile << " ";
      outfile << _supercell.getUniqueAtomMass(i);
    }
    outfile << std::endl;

    // forces + 1 line info about distortion
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
	// line info
	outfile << (_supercell.getUniqueAtomID(i) + 1) << " ";
	outfile << (j + 1) << " ";
	xvector<double> shift(3);
	shift = DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];
	outfile << setprecision(3);
	outfile << shift(1) << " " << shift(2) << " " << shift(3) << std::endl;
	// forces
	outfile << setprecision(6);
	for (int k = 0; k < _supercell.getNumberOfAtoms(); k++)
	  outfile << setw(15) << _uniqueForces[i][j][k](1)
		  << setw(15) << _uniqueForces[i][j][k](2)
		  << setw(15) << _uniqueForces[i][j][k](3) << std::endl;
      }
    }

    //
    //CO - START
    //outfile.clear();
    //outfile.close();

    aurostd::stringstream2file(outfile, DEFAULT_APL_DYNMAT_FILE);
    if (!aurostd::FileExist(DEFAULT_APL_DYNMAT_FILE)) {
      string function = "PhononCalculator::writeDYNMAT()";
      string message = "Cannot open output file " + DEFAULT_APL_DYNMAT_FILE + ".";
      throw aurostd::xerror(function, message, _FILE_ERROR_);
    }
//      throw apl::APLRuntimeError("PhononCalculator::writeDYNMAT(); Cannot open output file.");
    //CO - END
  }

  //////////////////////////////////////////////////////////////////////////////

  // This is the interface to phonopy code

  void PhononCalculator::writeFORCES() {
    //
    _logger << "Writing forces into file FORCES." << apl::endl;

    //
    //ifstream infile("SPOSCAR"); //CO
    xstructure ix;
    //CO - START
    string filename = "SPOSCAR";
    //if (infile.is_open()) {
    if (!aurostd::FileEmpty(filename)) {
      _logger << "Reading " << filename << apl::endl;
      stringstream SPOSCAR;
      aurostd::efile2stringstream(filename, SPOSCAR);
      SPOSCAR >> ix;
      //infile >> ix;
      //infile.close();
    } else
      ix = _supercell.getSupercellStructure();

    //
    //ofstream outfile("FORCES", ios_base::out);
    //if (!outfile.is_open())
    //  throw apl::APLRuntimeError("PhononCalculator::writeFORCES(); Cannot open output file.");
    stringstream outfile;
    //CO - END

    // 1st line
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++)
      dof += _uniqueDistortions[i].size();
    outfile << dof << std::endl;

    // forces + 1 line info about distortion
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
	// line info
	outfile << (_supercell.getUniqueAtomID(i) + 1) << " ";
	xvector<double> shift(3);
	shift = C2F(_supercell.getSupercellStructure().lattice, DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);
	outfile << setprecision(6);
	outfile << shift(1) << " " << shift(2) << " " << shift(3) << std::endl;
	// forces
	outfile << setprecision(6);
	for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
	  int l = 0;
	  for (; l < _supercell.getNumberOfAtoms(); l++)
	    if ((fabs(ix.atoms[k].cpos(1) - _supercell.getSupercellStructure().atoms[l].cpos(1)) < _AFLOW_APL_EPS_) &&
		(fabs(ix.atoms[k].cpos(2) - _supercell.getSupercellStructure().atoms[l].cpos(2)) < _AFLOW_APL_EPS_) &&
		(fabs(ix.atoms[k].cpos(3) - _supercell.getSupercellStructure().atoms[l].cpos(3)) < _AFLOW_APL_EPS_))
	      break;
	  //CO, not really mapping error, just mismatch between structure read in (ix) and current supercell structure (should be exact)
	  if (l == _supercell.getNumberOfAtoms()) {
	    cout << k << std::endl;
	    printXVector(ix.atoms[k].fpos);
	    printXVector(ix.atoms[k].cpos);
	    throw APLLogicError("apl::PhononCalculator::writeFORCES(); Mapping error.");
	  }

	  outfile << setw(15) << _uniqueForces[i][j][l](1) << " "
		  << setw(15) << _uniqueForces[i][j][l](2) << " "
		  << setw(15) << _uniqueForces[i][j][l](3) << std::endl;
	}
      }
    }

    //CO - START
    filename = "FORCES";
    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename))
      throw apl::APLRuntimeError("PhononCalculator::writeFORCES(); Cannot open output file.");

    //
    //outfile.clear();
    //outfile.close();
    //CO - END
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::writeXCrysDenForces() {
    _logger << "Writing forces into file XCrysDenForces." << apl::endl;
    //_supercell.center(0);   //JAHNATEK ORIGINAL
    _supercell.center_original();  //COREY

    stringstream outfile;  //CO
    // forces + 1 line info about distortion
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
	//string s = "FORCES_A" + stringify(_supercell.getUniqueAtomID(i)) + "D" + stringify(j) + ".xsf"; //CO
	outfile.str("");  //CO
	//ofstream outfile(s.c_str(), ios_base::out); //CO

	outfile << "CRYSTAL" << std::endl;
	outfile << "PRIMVEC 1" << std::endl;
	outfile << _supercell.getSupercellStructure().lattice << std::endl;
	outfile << "CONVEC 1" << std::endl;
	outfile << _supercell.getSupercellStructure().lattice << std::endl;
	outfile << "PRIMCOORD 1" << std::endl;
	outfile << _supercell.getNumberOfAtoms() << " 1" << std::endl;

	xvector<double> shift(3);
	shift = C2F(_supercell.getSupercellStructure().lattice, DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);

	outfile << setprecision(6);
	for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
	  outfile << _supercell.getAtomNumber(k) << " ";
	  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	  outfile << setprecision(8);
	  xvector<double> r = F2C(_supercell.getSupercellStructure().lattice,
				  _supercell.getSupercellStructure().atoms[k].fpos);
	  outfile << setw(15) << r(1) << setw(15) << r(2) << setw(15) << r(3) << " ";
	  // this is strange...
	  //outfile << setw(15) << _superCellStructure.atoms[k].cpos << " ";

	  // Scale force, it is expected in Hartree/Angs.
	  xvector<double> f = 27.212 * _uniqueForces[i][j][k];

	  outfile << setw(15) << f(1)
		  << setw(15) << f(2)
		  << setw(15) << f(3) << std::endl;
	}

	//CO - START
	string filename = "FORCES_A" + stringify(_supercell.getUniqueAtomID(i)) + "D" + stringify(j) + ".xsf";
	aurostd::stringstream2file(outfile, filename);
	if (!aurostd::FileExist(filename)) {
	  throw apl::APLRuntimeError("apl::PhononCalculator::writeXCrysDenForces(); Cannot create " + filename + " file.");
	}

	//outfile.clear();
	//outfile.close();
	//CO - END
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::hibernate() {
    _logger << "Hibernating..." << apl::endl;

    //
    //CO - START
    stringstream outfile;
    //ofstream outfile("apl.xml", ios_base::out);
    //if (!outfile.is_open())
    //  throw apl::APLRuntimeError("PhononCalculator::hibernate(); Cannot open output file.");
    //CO - END

    // XML declaration
    outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;

    // Our structure
    string tab = " ";
    outfile << "<apl>" << std::endl;

    // Info about calculation run
    outfile << tab << "<generator>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    outfile << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"Fletcher32\">"
	    << std::hex << getFileCheckSum("./" + _AFLOWIN_ + "") << "</i>" << std::endl;
    outfile << tab << "</generator>" << std::endl;

    // Unique distortions
    outfile << tab << "<distortions units=\"Angstrom\" cs=\"cartesian\">" << std::endl;

    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(8);
    outfile << tab << tab << "<i name=\"magnitude\">" << setw(15) << DISTORTION_MAGNITUDE << "</i>" << std::endl;

    outfile << tab << tab << "<varray>" << std::endl;
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      outfile << tab << tab << tab << "<varray atomID=\"" << _supercell.getUniqueAtomID(i) << "\">" << std::endl;
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
	outfile << tab << tab << tab << tab << "<v>";
	outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	outfile << setprecision(8);
	outfile << setw(15) << _uniqueDistortions[i][j](1) << " ";
	outfile << setw(15) << _uniqueDistortions[i][j](2) << " ";
	outfile << setw(15) << _uniqueDistortions[i][j](3);
	outfile << "</v>" << std::endl;
      }
      outfile << tab << tab << tab << "</varray>" << std::endl;
    }
    outfile << tab << tab << "</varray>" << std::endl;
    outfile << tab << "</distortions>" << std::endl;

    // Forces
    outfile << tab << "<forcefields units=\"eV/Angstrom\" cs=\"cartesian\">" << std::endl;
    outfile << tab << tab << "<varray>" << std::endl;
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      outfile << tab << tab << tab << "<varray atomID=\"" << _supercell.getUniqueAtomID(i) << "\">" << std::endl;
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
	outfile << tab << tab << tab << tab << "<varray distortion=\"" << j << "\">" << std::endl;
	for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
	  outfile << tab << tab << tab << tab << tab << "<v>";
	  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	  outfile << setprecision(15);
	  outfile << setw(24) << std::scientific << _uniqueForces[i][j][k](1) << " ";
	  outfile << setw(24) << std::scientific << _uniqueForces[i][j][k](2) << " ";
	  outfile << setw(24) << std::scientific << _uniqueForces[i][j][k](3);
	  outfile << "</v>" << std::endl;
	}
	outfile << tab << tab << tab << tab << "</varray>" << std::endl;
      }
      outfile << tab << tab << tab << "</varray>" << std::endl;
    }
    outfile << tab << tab << "</varray>" << std::endl;
    outfile << tab << "</forcefields>" << std::endl;

    // Force constant matrices
    outfile << tab << "<fcms units=\"eV/Angstrom^2\" cs=\"cartesian\" rows=\""
	    << _forceConstantMatrices.size() << "\" cols=\""
	    << _forceConstantMatrices[0].size() << "\">" << std::endl;

    outfile << tab << tab << "<varray>" << std::endl;
    for (uint i = 0; i < _forceConstantMatrices.size(); i++) {
      outfile << tab << tab << tab << "<varray row=\"" << i << "\">" << std::endl;
      for (uint j = 0; j < _forceConstantMatrices[i].size(); j++) {
	outfile << tab << tab << tab << tab << "<matrix row=\"" << i
		<< "\" col=\"" << j << "\">" << std::endl;
	for (int k = 1; k <= 3; k++) {
	  outfile << tab << tab << tab << tab << tab << "<v>";
	  for (int l = 1; l <= 3; l++) {
	    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	    outfile << setprecision(15);
	    outfile << setw(24) << std::scientific << _forceConstantMatrices[i][j](k, l) << " ";
	  }
	  outfile << "</v>" << std::endl;
	}
	outfile << tab << tab << tab << tab << "</matrix>" << std::endl;
      }
      outfile << tab << tab << tab << "</varray>" << std::endl;
    }
    outfile << tab << tab << "</varray>" << std::endl;
    outfile << tab << "</fcms>" << std::endl;

    //_logger << "APL-DEBUG  Only for polar materials" << apl::endl;
    // Only for polar materials
    if (_isPolarMaterial) {
      _logger << "APL-DEBUG Born effective charge tensors" << apl::endl;
      // Born effective charge tensors
      outfile << tab << "<born units=\"a.u.\" cs=\"cartesian\">" << std::endl;
      outfile << tab << tab << "<varray>" << std::endl;
      for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
	int id = _supercell.getUniqueAtomID(i);
	outfile << tab << tab << tab << "<matrix type=\"" << _supercell.getSupercellStructure().atoms[id].cleanname << "\">" << std::endl;
	for (int k = 1; k <= 3; k++) {
	  outfile << tab << tab << tab << tab << "<v>";
	  for (int l = 1; l <= 3; l++) {
	    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	    outfile << setprecision(8);
	    outfile << setw(15) << _bornEffectiveChargeTensor[i](k, l) << " ";
	  }
	  outfile << "</v>" << std::endl;
	}
	outfile << tab << tab << tab << "</matrix>" << std::endl;
      }
      outfile << tab << tab << "</varray>" << std::endl;
      outfile << tab << "</born>" << std::endl;

      // Dielectric constant matrix
      outfile << tab << "<epsilon units=\"a.u.\" cs=\"cartesian\">" << std::endl;
      outfile << tab << tab << "<matrix>" << std::endl;
      for (int k = 1; k <= 3; k++) {
	outfile << tab << tab << tab << "<v>";
	for (int l = 1; l <= 3; l++) {
	  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
	  outfile << setprecision(8);
	  outfile << setw(15) << _dielectricTensor(k, l) << " ";
	}
	outfile << "</v>" << std::endl;
      }
      outfile << tab << tab << "</matrix>" << std::endl;
      outfile << tab << "</epsilon>" << std::endl;
    }

    //
    outfile << "</apl>" << std::endl;

    //
    // outfile.clear();
    // outfile.close();

    // COREY, KBIN_ZIP will compress the whole directory, so just leave it alone
    // if (DOtar) {
    //   aurostd::stringstream2compressfile(kbinFlags.KZIP_BIN,outfile,"apl.xml");
    //   if (!aurostd::EFileExist("apl.xml"))
    //    throw apl::APLRuntimeError("PhononCalculator::hibernate(); Cannot open output apl.xml.");
    //} else {
    aurostd::stringstream2file(outfile, DEFAULT_APL_HARMIFC_FILE);
    if (!aurostd::FileExist(DEFAULT_APL_HARMIFC_FILE)) {
      string function = "PhononCalculator::hibernate()";
      string message = "Cannot open output file " + DEFAULT_APL_HARMIFC_FILE + ".";
      throw aurostd::xerror(function, message, _FILE_ERROR_);
//      throw apl::APLRuntimeError("PhononCalculator::hibernate(); Cannot open output apl.xml.");
    }
    //}

    // Compress
    //if (DOtar) aurostd::execute(string("COMPRESS -fz9 apl.xml"));
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::awake() {
    _logger << "Awakening..." << apl::endl;

    //CO, we already checked that it exists before, just open

    vector<string> vlines;                           //CO
    aurostd::efile2vectorstring(DEFAULT_APL_HARMIFC_FILE, vlines);  //CO
    // Decompress
    //bool isXMLCompressed = aurostd::FileExist(string("apl.xml.EXT")); //CO
    //if (isXMLCompressed) //CO
    //  aurostd::execute(string("EXT -d apl.xml.EXT")); //CO

    //
    //CO - START
    //ifstream infile("apl.xml", ios_base::in);
    //if (!infile.is_open())
    if (!vlines.size()) {
      string function = "PhononCalculator::awake()";
      string message = "Cannot open output file " + DEFAULT_APL_HARMIFC_FILE + ".";
      throw aurostd::xerror(function, message, _FILE_ERROR_);
//      throw apl::APLRuntimeError("apl::PhononCalculator::awake(); Cannot open input apl.xml.");
    }

    string line;
    uint line_count = 0;
    vector<string> tokens;

    // Test of xml...
    line = vlines[line_count++];
    //getline(infile, line);
    if (line.find("xml") == string::npos)
      throw APLLogicError("apl::PhononCalculator::awake(); Wrong xml file.");
    //CO - END

    // Get _AFLOWIN_ checksum and compare it to current
    while (true) {
      //getline(infile, line); //CO
      //if (infile.eof()) //CO
      if (line_count == vlines.size())  //CO
	throw APLLogicError("apl::PhononCalculator::awake(); Can not find <i name=\"checksum\" ...> tag.");
      line = vlines[line_count++];  //CO
      if (line.find("checksum") != string::npos)
	break;
    }
    int t = line.find_first_of(">") + 1;
    tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
    if (strtoul(tokens[0].c_str(), NULL, 16) != getFileCheckSum("./" + _AFLOWIN_ + ""))
      throw APLLogicError("apl::PhononCalculator::awake(); The " + _AFLOWIN_ + " file has been changed from the hibernated state.");
    tokens.clear();

    // Get force constant matrices
    while (true) {
      //getline(infile, line); //CO
      //if (infile.eof()) //CO
      if (line_count == vlines.size())  //CO
	throw APLLogicError("apl::PhononCalculator::awake(); Can not find <fcms> tag.");
      line = vlines[line_count++];  //CO
      if (line.find("fcms") != string::npos)
	break;
    }
    //CO - START
    line = vlines[line_count++];
    //getline(infile, line);
    line = vlines[line_count++];
    //getline(infile, line);
    //CO - END
    vector<xmatrix<double> > row;
    xmatrix<double> m(3, 3);
    while (true) {
      //getline(infile, line); //CO
      //if (infile.eof()) //CO
      if (line_count == vlines.size())  //CO
	throw APLLogicError("apl::PhononCalculator::awake(); Incomplete <fcms> tag.");
      line = vlines[line_count++];  //CO
      if (line.find("</varray>") != string::npos) {
	_forceConstantMatrices.push_back(row);
	row.clear();
	//getline(infile, line); //CO
	line = vlines[line_count++];  //CO
	if (line.find("</varray>") != string::npos)
	  break;
	//getline(infile, line); //CO
	line = vlines[line_count++];  //CO
      }

      for (int k = 1; k <= 3; k++) {
	//getline(infile, line); //CO
	line = vlines[line_count++];  //CO
	int t = line.find_first_of(">") + 1;
	tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
	m(k, 1) = aurostd::string2utype<double>(tokens.at(0));
	m(k, 2) = aurostd::string2utype<double>(tokens.at(1));
	m(k, 3) = aurostd::string2utype<double>(tokens.at(2));
	tokens.clear();
      }
      row.push_back(m);
      //getline(infile, line); //CO
      line = vlines[line_count++];  //CO
    }

    // Try to read born effective charges and dielectric constant
    try {  //CO
      //cerr << "APL-DEBUG Get born effective charge tensors" << std::endl; //CO
      // Get born effective charge tensors
      while (true) {
	//getline(infile, line); //CO
	//if (infile.eof()) { //CO
	if (line_count == vlines.size()) {  //CO
	  _isPolarMaterial = false;
	  throw APLLogicError("apl::PhononCalculator::awake(); Can not find <born> tag.");
	}
	line = vlines[line_count++];  //CO
	if (line.find("born") != string::npos)
	  break;
      }
      //getline(infile, line); //CO
      line = vlines[line_count++];  //CO
      while (true) {
	//getline(infile, line); //CO
	//if (infile.eof()) //CO
	if (line_count == vlines.size())  //CO
	  throw APLLogicError("apl::PhononCalculator::awake(); Incomplete <born> tag.");
	line = vlines[line_count++];  //CO
	if (line.find("</varray>") != string::npos)
	  break;
	for (int k = 1; k <= 3; k++) {
	  //getline(infile, line); //CO
	  line = vlines[line_count++];  //CO
	  int t = line.find_first_of(">") + 1;
	  tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
	  m(k, 1) = aurostd::string2utype<double>(tokens.at(0));
	  m(k, 2) = aurostd::string2utype<double>(tokens.at(1));
	  m(k, 3) = aurostd::string2utype<double>(tokens.at(2));
	  tokens.clear();
	}
	_bornEffectiveChargeTensor.push_back(m);
	//getline(infile, line); //CO
	line = vlines[line_count++];  //CO
      }

      // Get dielectric constant tensor
      while (true) {
	//getline(infile, line); //CO
	//if (infile.eof()) { //CO
	if (line_count == vlines.size()) {  //CO
	  _isPolarMaterial = false;
	  throw APLLogicError("apl::PhononCalculator::awake(); Can not find <epsilon> tag.");
	}
	line = vlines[line_count++];  //CO
	if (line.find("epsilon") != string::npos)
	  break;
      }
      //getline(infile, line); //CO
      line = vlines[line_count++];  //CO
      for (int k = 1; k <= 3; k++) {
	//getline(infile, line); //CO
	line = vlines[line_count++];  //CO
	int t = line.find_first_of(">") + 1;
	tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
	_dielectricTensor(k, 1) = aurostd::string2utype<double>(tokens.at(0));
	_dielectricTensor(k, 2) = aurostd::string2utype<double>(tokens.at(1));
	_dielectricTensor(k, 3) = aurostd::string2utype<double>(tokens.at(2));
	tokens.clear();
      }
      _inverseDielectricTensor = inverse(_dielectricTensor);
      _recsqrtDielectricTensorDeterminant = 1.0 / sqrt(determinant(_dielectricTensor));
    } catch (APLLogicError& e) {  //CO
      //_logger << apl::warning << e.what() << apl::endl;
    }

    //
    //infile.close(); //CO
    //infile.clear(); //CO

    // Compress
    //CO - START
    //if (DOtar)
    //  if (isXMLCompressed)
    //    aurostd::execute(string("COMPRESS -fz9 apl.xml"));
    //CO - END
  }

  // INTERFACE /////////////////////////////////////////////////////////////////

  const Supercell& PhononCalculator::getSupercell() { //CO 180409
    return _supercell;
  }

  const xstructure& PhononCalculator::getInputCellStructure() {
    //     cerr << _supercell.getInputStructure() << std::endl;
    return _supercell.getInputStructure();
  }

  const xstructure& PhononCalculator::getSuperCellStructure() {
    return _supercell.getSupercellStructure();
  }

  //CO - START
  double PhononCalculator::getEPS() {
    return _supercell.getEPS();
  }
  //CO - END

  uint PhononCalculator::getNumberOfBranches() {
    return (3 * _supercell.getInputStructure().atoms.size());
  }

  // ///////////////////////////////////////////////////////////////////////////
  void PhononCalculator::store_masses()  //[PINKU]
  {
    ATOMIC_MASSES_AMU.clear();
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();
    for (uint i = 0; i != pcAtomsSize; i++)
      ATOMIC_MASSES_AMU.push_back(_supercell.getAtomMass(_supercell.pc2scMap(i)));
  }
  // ///////////////////////////////////////////////////////////////////////////
  void PhononCalculator::get_special_inputs(string& AflowIn)  //[PINKU]
  {
    _check_LDAU2_ON = "";
    _LDAU_PARAMETERS = "";
    _PSTRESS = "";
    string line;
    vector<string> vlines;
    uint line_count = 0;
    aurostd::string2vectorstring(AflowIn,vlines);
    //aurostd::efile2vectorstring(_AFLOWIN_, vlines); //CO 171003
    //ifstream myfile(_AFLOWIN_.c_str());

    //CO - START
    //if (!myfile.is_open()) {
    if (!vlines.size()) {
      throw apl::APLRuntimeError("apl::PhononCalculator::get_special_inputs(); Cannot read ["+_AFLOWIN_+"] file.");
    }
    //while (getline(myfile, line)) {
    while (line_count < vlines.size()) {
      line = vlines[line_count++];
      if (line == "") continue;
      if (line[0] == '#') continue;
      if ((line.find("LDAU2=ON") != std::string::npos)) _check_LDAU2_ON = line;
      if ((line.find("LDAU_PARAMETERS") != std::string::npos)) _LDAU_PARAMETERS = line;
      if ((line.find("PSTRESS") != std::string::npos)) _PSTRESS = line;
    }
    //myfile.clear();
    //myfile.close();
    //CO - END
  }
  // ///////////////////////////////////////////////////////////////////////////

  //CO 180214 - START

  void PhononCalculator::get_NCPUS(string& ncpus) {
    ncpus="MAX";
    if(_kbinFlags.KBIN_MPI_NCPUS>0){ncpus=aurostd::utype2string(_kbinFlags.KBIN_MPI_NCPUS);}
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")){ncpus=XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS");}
  }

  void PhononCalculator::get_NCPUS(int& ncpus) {
    string ncpus_str;
    get_NCPUS(ncpus_str);
    if(ncpus_str=="MAX"){ncpus=MPI_NCPUS_MAX;return;}
    ncpus=aurostd::string2utype<int>(ncpus_str);
    if(ncpus<1){ncpus=1;}
  }

  //CO 180214 - STOP

  // BEGIN ME 180518
  //filesExistPhonons//////////////////////////////////////////////////////////////
  bool PhononCalculator::filesExistPhonons(_xinput& xinp) {
    string dir = xinp.getDirectory() + string("/");
    if (aurostd::FileExist(dir + _AFLOWIN_)) {
      return true;  //do not OVERWRITE an aflow.in
    }
    if (_kbinFlags.AFLOW_MODE_VASP){
      if(aurostd::EFileExist(dir + string("vasprun.xml.static")) ||
	 aurostd::EFileExist(dir + string("vasprun.xml"))) {
	return true;
      }
    }
    if(_kbinFlags.AFLOW_MODE_AIMS){
      if(aurostd::EFileExist(dir + string("aims.out"))) {
	return true;
      }
    }
    return false;
  }

  //createAflowInPhonons////////////////////////////////////////////////////////
  void PhononCalculator::createAflowInPhonons(_xinput& xinp, const string& runname) {
    if(xinp.AFLOW_MODE_VASP){
      // Switch off autotune, because....
      _kbinFlags.KBIN_MPI_AUTOTUNE = true;
      // Common KPOINTS settings and OVERRIDES
      xinp.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
      if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry) {
	xinp.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.content_string;
      }
      xinp.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;
      if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry) {
        xinp.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.content_uint;
      }
      // Clear old INCAR and set it as we want...
      // ME TODO: must find a better way - this removes all information on Hubbard, MAGMOM, etc.
      xinp.xvasp.INCAR.str(std::string());
      string system;
      for (uint i = 0; i < xinp.getXStr().species.size(); i++) {
	system = system + xinp.getXStr().species_pp.at(i);
      }
      system = system + "@" + runname;
      xinp.xvasp.INCAR << "SYSTEM=" << system << std::endl;
      xinp.xvasp.INCAR << "# Added by [AFLOW_APL] begin" << std::endl;
      xinp.xvasp.INCAR << "NELMIN=4         # The forces have to be well converged" << std::endl;
      xinp.xvasp.INCAR << "NELM = 120       # Many electronic steps (SC2013)" << std::endl;
      xinp.xvasp.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
      xinp.xvasp.INCAR << "# Added by [AFLOW_APL] end" << std::endl;

      // Change format of POSCAR
      if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
	  (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos))) {
	xinp.getXStr().is_vasp5_poscar_format = false;
      }
    }
    if(xinp.AFLOW_MODE_AIMS) {
      xinp.xaims.CONTROL.str(std::string());
      KBIN::AIMS_Produce_CONTROL(xinp.xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);  //DEFAULT
      KBIN::AIMS_Modify_CONTROL(xinp.xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);            //DEFAULT
    }

    // Create _AFLOWIN_
    writeOUTPUT(xinp); //CO 171009
  }
  //outfileFoundAnywherePhonons/////////////////////////////////////////////////
  bool PhononCalculator::outfileFoundAnywherePhonons(vector<_xinput>& xinps) {
    for (uint idxRun = 0; idxRun < xinps.size(); idxRun++) {
      string dir = xinps[idxRun].getDirectory();
      // If tarred and compressed directory exists...

      // COREY CHECK THIS
      deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
      string tarfilename = dir + ".tar";
      for(uint iext=0;iext<vext.size();iext++) {
	if (aurostd::FileExist(tarfilename+vext.at(iext))) {
	  aurostd::execute(string("tar -xf ") + tarfilename.at(iext));
	}
      }
      if(_kbinFlags.AFLOW_MODE_VASP) {
	if(aurostd::EFileExist(dir + string("/vasprun.xml.static")) || 
	   aurostd::EFileExist(dir + string("/vasprun.xml"))) {
	  return true;
	}
      }
      if(_kbinFlags.AFLOW_MODE_AIMS) {
	if(aurostd::EFileExist(xinps[idxRun].getDirectory() + string("/aims.out"))) {
	  return true;
	}
      }
    }
    return false;
  }


  //outfileFoundEverywherePhonons///////////////////////////////////////////////
  void PhononCalculator::outfileFoundEverywherePhonons(vector<_xinput>& xinps) {
    for (uint idxRun = 0; idxRun < xinps.size(); idxRun++) {
      string tarfilename = xinps[idxRun].getDirectory() + ".tar";
      xinps[idxRun].getXStr().qm_forces.clear();
      // Load data....
      if(_kbinFlags.AFLOW_MODE_VASP){
	string vasprunxml_file=xinps[idxRun].getDirectory() + string("/vasprun.xml.static");
	if(!aurostd::EFileExist(vasprunxml_file)) {
	  vasprunxml_file=xinps[idxRun].getDirectory() + string("/vasprun.xml");
	  if(!aurostd::EFileExist(vasprunxml_file)) {
	    _logger << apl::warning << "The vasprun.xml file in " << xinps[idxRun].getDirectory() << " directory is missing." << apl::endl;
	    throw APLRuntimeError("apl::DirectMethodPC::runVASPCalculations(); Missing data from one job.");
	  }
	}
	xVASPRUNXML vasprunxml(vasprunxml_file);
	for (uint i = 0; i < vasprunxml.vforces.size(); i++) xinps[idxRun].getXStr().qm_forces.push_back(vasprunxml.vforces.at(i));
	if (int(xinps[idxRun].getXStr().qm_forces.size()) == _supercell.getNumberOfAtoms()) { xinps[idxRun].getXStr().qm_calculated = TRUE; }
      }
      if(_kbinFlags.AFLOW_MODE_AIMS){
	if(!aurostd::EFileExist(xinps[idxRun].getDirectory() + string("/aims.out"))) {
	  _logger << apl::warning << "The aims.out file in " << xinps[idxRun].getDirectory() << " directory is missing." << apl::endl;
	  throw APLRuntimeError("apl::DirectMethodPC::runAIMSCalculations(); Missing data from one job.");
	}
	xAIMSOUT xaimsout(xinps[idxRun].getDirectory() + "/aims.out");
	for (uint i = 0; i < xaimsout.vforces.size(); i++) xinps[idxRun].getXStr().qm_forces.push_back(xaimsout.vforces.at(i));
	if (int(xinps[idxRun].getXStr().qm_forces.size()) == _supercell.getNumberOfAtoms()) { xinps[idxRun].getXStr().qm_calculated = TRUE; }
      }

      // Was it all right?
      if (!xinps[idxRun].getXStr().qm_calculated) {
	if(_kbinFlags.AFLOW_MODE_VASP){
	  _logger << apl::warning << "The vasprun.xml file in " << xinps[idxRun].getDirectory() << " is wrong." << apl::endl;
	  throw APLRuntimeError("apl::DirectMethodPC::runVASPCalculations(); Missing data from one job.");
	}
	if(_kbinFlags.AFLOW_MODE_AIMS){
	  _logger << apl::warning << "The aims.out file in " << xinps[idxRun].getDirectory() << " is wrong." << apl::endl;
	  throw APLRuntimeError("apl::DirectMethodPC::runAIMSCalculations(); Missing data from one job.");
	}
      }

      // Pack/Remove the whole directory...
      if (DOtar) {
	if (!aurostd::EFileExist(tarfilename)) {
	  aurostd::execute(string("tar -cf ") + tarfilename + " " + xinps[idxRun].getDirectory() + "/");
	  aurostd::CompressFile(tarfilename,_kbinFlags.KZIP_BIN);
	}
	if (aurostd::FileExist(tarfilename)) aurostd::execute(string("rm -rf ") + xinps[idxRun].getDirectory() + "/");
      }
    }
  }
  
  void PhononCalculator::subtractZeroStateForces(vector<_xinput>& xinps) {
    for (uint idxRun = 0; idxRun < xInputs.size() - 1; idxRun++) {
      for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
	xinps[idxRun].getXStr().qm_forces[k](1) = xinps[idxRun].getXStr().qm_forces[k](1) - xinps.back().getXStr().qm_forces[k](1);
	xinps[idxRun].getXStr().qm_forces[k](2) = xinps[idxRun].getXStr().qm_forces[k](2) - xinps.back().getXStr().qm_forces[k](2);
	xinps[idxRun].getXStr().qm_forces[k](3) = xinps[idxRun].getXStr().qm_forces[k](3) - xinps.back().getXStr().qm_forces[k](3);
      }
    }
  }
  // END ME 180518

}  // namespace apl
