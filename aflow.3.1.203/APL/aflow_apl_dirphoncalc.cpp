#include "aflow_apl.h"

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

DirectMethodPC::DirectMethodPC(Supercell& sc, StrPairs& strpair, _xinput& xinput, //_xvasp& xvasp,
                               _aflags& aflags, _kflags& kflags,
                               _xflags& xflags, //_vflags& vflags, 
                               string& AflowIn, //this is the file CONTENTS, not the file path
                               Logger& l)
    : PhononCalculator(sc, strpair, xinput, aflags, kflags, xflags, AflowIn, l) { //xvasp, aflags, kflags, vflags, l) {
  //GENERATE_PLUS_MINUS = false; //JAHNATEK ORIGINAL
  AUTO_GENERATE_PLUS_MINUS = true;   //CO
  USER_GENERATE_PLUS_MINUS = false;  //CO
  GENERATE_ONLY_XYZ = false;
  _isPolarMaterial = false;
}

// ///////////////////////////////////////////////////////////////////////////

DirectMethodPC::~DirectMethodPC() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////

void DirectMethodPC::clear() {
}

//////////////////////////////////////////////////////////////////////////////

void DirectMethodPC::calculateForceFields() {
  // Check if supercell is already built
  if (!_supercell.isConstructed()) {
    throw APLRuntimeError("apl::DirectMethodPC::calculateForceFields(); The supercell structure has not been initialized yet.");
  }

  // Determine the distortion vectors
  estimateUniqueDistortions(_supercell.getSupercellStructure(), _uniqueDistortions);

  // Print some information
  int dof = 0;
  for (uint i = 0; i < _uniqueDistortions.size(); i++)
    dof += _uniqueDistortions[i].size();
  _logger << "Found " << dof << " degree(s) of freedom." << apl::endl;
  for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
    int id = _supercell.getUniqueAtomID(i);
    for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
      _logger << "Atom [" << sf("%03d") << id << "] (" << sf("%f")
              << sw(2) << _supercell.getSupercellStructure().atoms[id].cleanname
              << ") will be distorted in direction ["
              << sf("%5.3f") << _uniqueDistortions[i][j](1) << ","
              << sf("%5.3f") << _uniqueDistortions[i][j](2) << ","
              << sf("%5.3f") << _uniqueDistortions[i][j](3) << "]." << apl::endl;
    }
  }

  // Call VASP to calculate forces
  runVASPCalculations();

  // Optional
  //if( _isPolarMaterial )
  //    removeFakeForcesOfElectrostaticField();
}

//////////////////////////////////////////////////////////////////////////////

void DirectMethodPC::estimateUniqueDistortions(const xstructure& xstr,
                                               vector<vector<xvector<double> > >& uniqueDistortions) {
  //COREY NOTES ON THIS FUNCTION
  // - this function creates symmetrically unique distortion vectors for each iatom
  // - you can have at most 3 unique (orthogonal) distortions per atom, but probably fewer considering symmetry
  // - distortions are relative to the lattice vectors
  // - we first generate different types of basic distortions (along lattice vectors, diagonal, body diagonal) in fractional coordinates
  // - convert fractional to cartesian, now distortions truly represent those along lattice vectors, etc.
  // - despite that the vector order changes with every loop, these basic distortions are not changed
  // - for each iatom,
  //   - we build three lists:
  //     - allDistortionsOfAtom (horrible name) is simply a list of symmetrically equivalent (orthogonal) distortions
  //       to each basic distortion (by rotations of the crystal) - this is true because we clear allDistortionsOfAtom with every loop,
  //       so original distortionVector is unchanged
  //     - uniqueDistortions does not play a role in this loop right now, so ignore
  //     - testVectorDim is a count of symmetrically equivalent distortions for each basic distortion
  //   - we sort the basic distortions by testVectorDim (largest first), i.e. the basic distortions that generate the highest count of
  //     equivalent distortions go first (you could consider it like the most "natural" choices for distortions based on crystal symmetry)

  // Is there a list of inequivalent atoms?
  if (!xstr.iatoms_calculated) {
    throw APLRuntimeError("apl::DirectMethodPC::estimateUniqueDistortions(); The list of the inequivalent atoms is missing.");
  }

  // Clear old stuff
  for (uint i = 0; i < uniqueDistortions.size(); i++) {
    uniqueDistortions[i].clear();
  }
  uniqueDistortions.clear();

  // Determine the distortion vectors for this system
  if (xstr.agroup_calculated) {
    // Test directions for distortion (in cartesian coordinates)
    vector<xvector<double> > testDistortions;
    xvector<double> testVec(3);

    if (GENERATE_ONLY_XYZ) {
      // Test distortion (in cartesian coordinates) along x, y, and z
      // - they have to be orthonormal
      testVec(1) = 1.0;
      testVec(2) = 0.0;
      testVec(3) = 0.0;
      testDistortions.push_back(testVec);
      testVec(1) = 0.0;
      testVec(2) = 1.0;
      testVec(3) = 0.0;
      testDistortions.push_back(testVec);
      testVec(1) = 0.0;
      testVec(2) = 0.0;
      testVec(3) = 1.0;
      testDistortions.push_back(testVec);
    } else {
      // Add lattice vectors
      testVec(1) = 1.0;
      testVec(2) = 0.0;
      testVec(3) = 0.0;
      testDistortions.push_back(testVec);
      testVec(1) = 0.0;
      testVec(2) = 1.0;
      testVec(3) = 0.0;
      testDistortions.push_back(testVec);
      testVec(1) = 0.0;
      testVec(2) = 0.0;
      testVec(3) = 1.0;
      testDistortions.push_back(testVec);

      // Add diagonal vectors
      testVec(1) = 1.0;
      testVec(2) = 1.0;
      testVec(3) = 0.0;
      testDistortions.push_back(testVec);
      testVec(1) = 1.0;
      testVec(2) = 0.0;
      testVec(3) = 1.0;
      testDistortions.push_back(testVec);
      testVec(1) = 0.0;
      testVec(2) = 1.0;
      testVec(3) = 1.0;
      testDistortions.push_back(testVec);

      // Add body diagonal vectors
      testVec(1) = 1.0;
      testVec(2) = 1.0;
      testVec(3) = 1.0;
      testDistortions.push_back(testVec);

      // Convert to cartesian representation
      for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
        testDistortions[idTestVector] = F2C(xstr.lattice, testDistortions[idTestVector]);
      }

      // Norm to one (so it will be scaled by user magnitude of distortion)
      for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
        double len = aurostd::modulus(testDistortions[idTestVector]);
        testDistortions[idTestVector](1) /= len;
        testDistortions[idTestVector](2) /= len;
        testDistortions[idTestVector](3) /= len;
      }
    }

    // Loop over inequivalent atoms
    vector<xvector<double> > uniqueDistortionsOfAtom;
    vector<xvector<double> > allDistortionsOfAtom;
    for (uint i = 0; i < xstr.iatoms.size(); i++) {
      int inequivalentAtomID = xstr.iatoms[i][0];
      //cout << "inequivalentAtomID = " << inequivalentAtomID << std::endl;
      if (xstr.agroup[inequivalentAtomID].size() == 0) {
        throw APLRuntimeError("apl::DirectMethodPC::estimateUniqueDistortions(); Site point group operations are missing.");
      }

      // Loop over test directions vectors - we count the number of unique
      // which can be obtain from the each test vector
      vector<int> testVectorDim;
      for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
        // Get the list of all equivalent distortions
        testDistortion(testDistortions[idTestVector],
                       xstr.agroup[inequivalentAtomID], allDistortionsOfAtom,
                       uniqueDistortionsOfAtom);
        testVectorDim.push_back(allDistortionsOfAtom.size());
        allDistortionsOfAtom.clear();
        uniqueDistortionsOfAtom.clear();
      }

      // Now, we order all test vectors according to their dimensionality
      // Those with the highest score go first
      for (uint j = 0; j < testDistortions.size() - 1; j++) {
        for (uint k = j + 1; k < testDistortions.size(); k++) {
          if (testVectorDim[k] > testVectorDim[j]) {
            xvector<double> temp = testDistortions[k];
            testDistortions[k] = testDistortions[j];
            testDistortions[j] = temp;
            int itemp = testVectorDim[k];
            testVectorDim[k] = testVectorDim[j];
            testVectorDim[j] = itemp;
          }
        }
      }

      // Now we are going again, but slightly different, we count all
      // generated directions together until the total count is lower than three
      for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
        testDistortion(testDistortions[idTestVector],
                       xstr.agroup[inequivalentAtomID], allDistortionsOfAtom,
                       uniqueDistortionsOfAtom);
        if (allDistortionsOfAtom.size() >= 3) break;
      }

      //cout << "XXXXX  Number of unique distortion vectors for atom ["
      //     << inequivalentAtomID << "] = " << uniqueDistortionsOfAtom.size() << std::endl;
      uniqueDistortions.push_back(uniqueDistortionsOfAtom);
      // Free useless stuff
      allDistortionsOfAtom.clear();
      uniqueDistortionsOfAtom.clear();
    }
    // Free useless stuff
    testDistortions.clear();
  } else {
    throw APLRuntimeError("apl::DirectMethodPC::estimateUniqueDistortions(); The list of the site point group operations is missing.");
  }
}

//////////////////////////////////////////////////////////////////////////////

void DirectMethodPC::testDistortion(const xvector<double>& distortionVector,
                                    const vector<_sym_op>& symPool,
                                    vector<xvector<double> >& allDistortionsOfAtom,
                                    vector<xvector<double> >& uniqueDistortionsOfAtom) {
  // Test if it is unique direction
  // Use the Gram–Schmidt method for orthogonalizing, if the final vectors
  // is non zero length -> is unique
  xvector<double> testDistortion = distortionVector;
  for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
    testDistortion = testDistortion - getVectorProjection(testDistortion, allDistortionsOfAtom[k]);
  }
  if (aurostd::modulus(testDistortion) > _AFLOW_APL_EPS_) {
    // Normalize vector
    testDistortion = testDistortion / aurostd::modulus(testDistortion);
    // Store this vector (in Cartesian form!)
    allDistortionsOfAtom.push_back(testDistortion);
    uniqueDistortionsOfAtom.push_back(testDistortion);
    //cout << "new unique distortion vector: " << testDistortion << std::endl;
  }

  // Apply all symmetry operations on test vector and generate all
  // unique directions and store them for future comparison
  for (uint iSymOp = 0; iSymOp < symPool.size(); iSymOp++) {
    testDistortion = symPool[iSymOp].Uc * distortionVector;

    for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
      testDistortion = testDistortion - getVectorProjection(testDistortion, allDistortionsOfAtom[k]);
    }
    if (aurostd::modulus(testDistortion) > _AFLOW_APL_EPS_) {
      testDistortion = testDistortion / aurostd::modulus(testDistortion);
      allDistortionsOfAtom.push_back(testDistortion);
      //cout << "new symmetry generated distortion vector: " << testDistortion << std::endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

//CO - START
bool DirectMethodPC::needMinus(uint ineq_atom_indx, uint distortion_indx) {
  //bool need_minus = true;
  uint atom_index = _supercell.getUniqueAtomID(ineq_atom_indx);
  vector<_sym_op> agroup = _supercell.getAGROUP(atom_index);
  xvector<double> rotated_distortion, distortion_sum;
  //cerr << agroup.size() << std::endl;
  //cerr << "distortion  : " << _uniqueDistortions[ineq_atom_indx][distortion_indx] << std::endl;
  for (uint i = 0; i < agroup.size(); i++) {
    rotated_distortion = agroup[i].Uc * _uniqueDistortions[ineq_atom_indx][distortion_indx];
    //cerr << "rdistortion : " << rotated_distortion << std::endl;
    if (identical(_uniqueDistortions[ineq_atom_indx][distortion_indx], -rotated_distortion, _AFLOW_APL_EPS_)) {  //just mimicking Jahnatek tolerance here
      return FALSE;
      //need_minus = false;
      //break;
    } else {
      continue;
    }
  }
  return TRUE;
  //return need_minus;
}
//CO - END

//////////////////////////////////////////////////////////////////////////////

void DirectMethodPC::runVASPCalculations() {
  // ****** BEGIN JJPR: DEFINED in aflow_apl.h because it will be used in other functions ******
  // All runs will be stored here
  //vector<_xvasp> xInputs;
  // ****** END JJPR ******

  //CO - START
  vector<bool> vgenerate_plus_minus;  //CO
  bool generate_plus_minus;           //CO
  //bool         check_minus_needed = ( AUTO_GENERATE_PLUS_MINUS && !USER_GENERATE_PLUS_MINUS && !_supercell.isDerivativeStructure() );
  bool check_minus_needed = (AUTO_GENERATE_PLUS_MINUS && !USER_GENERATE_PLUS_MINUS);
  for (uint i = 0; i < _uniqueDistortions.size(); i++) {
    for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
      vgenerate_plus_minus.push_back(true);  //assume we need plus/minus
    }
  }
  //CO - END
  // Generate calculation directories
  for (uint i = 0; i < _uniqueDistortions.size(); i++) {
    for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
      //CO - START
      if (check_minus_needed) {
        vgenerate_plus_minus[j] = needMinus(i, j);
        if (!vgenerate_plus_minus.back()) {
          _logger << "No negative distortion needed for distortion " << j << "." << apl::endl;
        }
      }
      generate_plus_minus = vgenerate_plus_minus[j];
      for (uint k = 0; k < (generate_plus_minus ? 2 : 1); k++) {
        //CO - END
        // Copy settings from common case
        xInputs.push_back(_xInput);
        int idxRun = xInputs.size() - 1;

        // Create run id name
        string runname = "./" + _AFLOW_APL_AFLOW_DIRECTORY_PREFIX_ + "APL" + stringify(idxRun) + "A" + stringify(_supercell.getUniqueAtomID(i)) + "D" + stringify(j);

        if (generate_plus_minus) {  //CO
          runname = runname + ((k == 0) ? "P" : "M");
        }

        // Setup working directory
        xInputs[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname); // = _xInput.getDirectory() + "/" + runname;

        // Generate POSCAR and apply the unique distortion for one inequvalent atom
        // This distortion vector is stored in Cartesan form, hence use C2F before applying
        //xInputs[idxRun].getXStr()                    = _supercell.getSupercellStructure(); //CO very slow
        //int idAtom                                 = xInputs[idxRun].getXStr().iatoms[i][0];
        xInputs[idxRun].setXStr(_supercell.getSupercellStructureLight()); // = _supercell.getSupercellStructureLight();  //CO faster, only what's necessary here
        int idAtom = _supercell.getUniqueAtomID(i);
        xInputs[idxRun].getXStr().atoms[idAtom].fpos = xInputs[idxRun].getXStr().atoms[idAtom].fpos + C2F(xInputs[idxRun].getXStr().lattice, ((k == 0) ? 1.0 : -1.0) * DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);
        xInputs[idxRun].getXStr().atoms[idAtom].cpos = F2C(xInputs[idxRun].getXStr().lattice,
                                                         xInputs[idxRun].getXStr().atoms[idAtom].fpos);

        // ME 180518 Moved code into separate functions, so it's accessible to AAPL as well
        //CO - START
        // If there are already LOCK or vasprun.xml.static files, it means this directory was already generated and computed,
        // hence do not touch and leave, but store this structure in the
        // list, hence it will be used in next part of code.
        if (filesExistPhonons(xInputs[idxRun])) {
          continue;
        }
        _logger << "Creating " << xInputs[idxRun].getDirectory() << apl::endl;
        //CO - START
        createAflowInPhonons(xInputs[idxRun], runname);
        //if( !AVASP_MakeSingleAFLOWIN(xInputs[idxRun]) )
        //    throw APLRuntimeError("apl::DirectMethodPC::runVASPCalculations(); Can not create [" << _AFLOWIN_ << "] in subdirectory.");
      }
    }
  }

  // Add zero state if requested
  if (_calculateZeroStateForces) {
    // Copy settings from common case
    xInputs.push_back(_xInput);
    int idxRun = xInputs.size() - 1;

    // Create run id name
    string runname = "./" + _AFLOW_APL_AFLOW_DIRECTORY_PREFIX_ + "APL" + stringify(idxRun) + "ZEROSTATE";

    // Setup working directory
    xInputs[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname);  // = _xInput.getDirectory() + "/" + runname;

    // Get structure
    //xInputs[idxRun].setXStr(_supercell.getSupercellStructure());  // = _supercell.getSupercellStructure(); //CO
    xInputs[idxRun].setXStr(_supercell.getSupercellStructureLight()); // = _supercell.getSupercellStructureLight();  //CO

    // If there is already a LOCK file, it means this directory was already generated
    // and computed, hence do not touch and leave, but store this structure in the
    // list, hence it will be used in next part of code.
    // If not, continue in this way, prepare generation of _AFLOWIN_ ...
    bool create_zero_state_output=true;
    if(aurostd::FileExist(xInputs[idxRun].getDirectory() + string("/") + _AFLOWIN_)){create_zero_state_output=false;}  //CO 180406
    if(_kbinFlags.AFLOW_MODE_VASP){
      if(aurostd::EFileExist(xInputs[idxRun].getDirectory() + string("/vasprun.xml.static")) ||
         aurostd::EFileExist(xInputs[idxRun].getDirectory() + string("/vasprun.xml"))){create_zero_state_output=false;}
    }
    if(_kbinFlags.AFLOW_MODE_AIMS){
      if(aurostd::EFileExist(xInputs[idxRun].getDirectory() + string("/aims.out"))){create_zero_state_output=false;}
    }
    //if (!aurostd::FileExist(xInputs[idxRun].getDirectory() + string("/") + _AFLOWLOCK_)) {  //CO
                                                                                           // If not, continue in this way, prepare generation of _AFLOWIN_ ...

    if(create_zero_state_output){
      _logger << "Creating " << xInputs[idxRun].getDirectory() << apl::endl;  //CO
      createAflowInPhonons(xInputs[idxRun], runname);
    }
  }

  // BEGIN STEFANO
  // Do calculation an additional calculations for polar materials
  if (_isPolarMaterial) {
    try {
      // if tarred and compressed directory exists...

      // COREY CHECK THIS
      deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
      string tarfilename =string(_AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_) + ".tar";
      for(uint iext=0;iext<vext.size();iext++) {
	if (aurostd::FileExist(tarfilename+vext.at(iext))) {
	  if(_kbinFlags.AFLOW_MODE_VASP)
	    aurostd::execute(string("tar -xf ") + tarfilename + vext.at(iext)+ " --wildcards " + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/OUTCAR*");
	  if(_kbinFlags.AFLOW_MODE_AIMS)
	    aurostd::execute(string("tar -xf ") + tarfilename + vext.at(iext)+ " --wildcards " + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/aims.out");
	}
      }

      // Calc. Born effective charge tensors and dielectric constant matrix
      runVASPCalculationsBE();
      // Parse it from OUTCAR
      if(_kbinFlags.AFLOW_MODE_VASP){readBornEffectiveChargesFromOUTCAR();}
      if(_kbinFlags.AFLOW_MODE_AIMS){readBornEffectiveChargesFromAIMSOUT();}
      // Enforce ASR (Acustic sum rules)
      symmetrizeBornEffectiveChargeTensors();
      // Parser epsilon from OUTCAR
      readDielectricTensorFromOUTCAR();
      //
      _logger << "Dielectric tensor: ";
      for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
          _logger << sf("%5.3f") << _dielectricTensor(a, b) << " ";
      _logger << apl::endl;

      // precompute
      _inverseDielectricTensor = inverse(_dielectricTensor);
      _recsqrtDielectricTensorDeterminant = 1.0 / sqrt(determinant(_dielectricTensor));

      // Pack the whole directory...
      if (DOtar) {
	// COREY CHECK THIS
        if (!aurostd::FileExist(tarfilename)) {
	  aurostd::execute(string("tar -cf ") + tarfilename + " " + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/");
	  aurostd::CompressFile(tarfilename,_kbinFlags.KZIP_BIN);
	}
 	if (aurostd::FileExist(tarfilename)) aurostd::execute(string("rm -rf ") + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/");
      }
    } catch (APLLogicError& e) {
      _logger << apl::error << e.what() << apl::endl;
      _logger << warning << "Switching the dipole-dipole correction off." << apl::endl;
      _isPolarMaterial = false;
    }
    //      cerr << "GENERATING FOR POLAR  STOP" << std::endl;exit(0);
  }
  // END STEFANO

  // Extract all forces from vasprun.xml.static ////////////////////////////////////////

  //first pass, just find if outfile is found ANYWHERE
  if(!outfileFoundAnywherePhonons(xInputs)){throw APLStageBreak();}

  //second pass, make sure it's everywhere!
  outfileFoundEverywherePhonons(xInputs);
    //CO - END

  // Remove zero state forces if possible //////////////////////////////////

  if (_calculateZeroStateForces) {
    subtractZeroStateForces(xInputs);
  }

  // Store forces //////////////////////////////////////////////////////////

  int idxRun = 0;
  for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
    vector<vector<xvector<double> > > forcesForOneAtomAndAllDistortions;
    for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
      generate_plus_minus = vgenerate_plus_minus[j];  //CO
      vector<xvector<double> > forcefield;
      xvector<double> drift(3);
      for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
        xvector<double> force(3);
        force(1) = xInputs[idxRun].getXStr().qm_forces[k](1);
        force(2) = xInputs[idxRun].getXStr().qm_forces[k](2);
        force(3) = xInputs[idxRun].getXStr().qm_forces[k](3);
        if (generate_plus_minus) {  //CO
          force(1) = 0.5 * (force(1) - xInputs[idxRun + 1].getXStr().qm_forces[k](1));
          force(2) = 0.5 * (force(2) - xInputs[idxRun + 1].getXStr().qm_forces[k](2));
          force(3) = 0.5 * (force(3) - xInputs[idxRun + 1].getXStr().qm_forces[k](3));
        }
        forcefield.push_back(force);
        drift = drift + force;
      }

      if (generate_plus_minus) {  //CO
        idxRun += 2;
      } else {
        idxRun++;
      }

      // Remove drift
      drift(1) = drift(1) / forcefield.size();
      drift(2) = drift(2) / forcefield.size();
      drift(3) = drift(3) / forcefield.size();
      for (_AFLOW_APL_REGISTER_ uint k = 0; k < forcefield.size(); k++) {
        forcefield[k] = forcefield[k] - drift;
      }

      // Store
      forcesForOneAtomAndAllDistortions.push_back(forcefield);
      forcefield.clear();
    }
    _uniqueForces.push_back(forcesForOneAtomAndAllDistortions);
    forcesForOneAtomAndAllDistortions.clear();
  }

  //******** BEGIN JJPR: Clean later to use ZERO STATE FORCES with TENSOR  *******
  // Clear useless stuff
  if (!KAPPA) {
    for (uint i = 0; i < xInputs.size(); i++) {
      xInputs[i].clear();
    }
    xInputs.clear();
  }
  //****** END JJPR *******
}

//////////////////////////////////////////////////////////////////////////////

}  // namespace apl
