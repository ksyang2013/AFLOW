// ***************************************************************************
// *                                                                         *
// *              AFlow JOSE J. PLATA   - Duke University 2017               *
// *                                                                         *
// ***************************************************************************
// tencalc_contrib_jose.cpp
// written by JJPR 2013
// 2013: jose.plata@duke.edu

// These functions are integrated in the phononcalculator developed by Michal Jahnatek.
// If TENSOR option is activated in aflow.in the anharmonic force  iconstant tensor will be calculated.
// These functions prepare inputs, store forces from static calculation and build the tensor taking into account:
//        1) Symmetry of the crystal
//        2) Permutation rules in the tensor
//        3) Zero force corrections (optional)
// ***************************************************************************

//#include <iostream>
//#include <sstream>
//#include <strings>
//#include <limits>

#include "aflow_apl.h"
//#include "aplmath.h"
//#include "shellhandle.h"

// Some parts are written within the C++0x support in GCC, expecially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

using namespace std;

namespace apl {

//VaspTensor////////////////////////////////////////////////////////////////////////////
// Prepare calculations for Tensor/////

void PhononCalculator::VaspTensor() {
  // Check if supercell is already built----------------------------------
  if (!_supercell.isConstructed()) {
    throw APLRuntimeError("apl::TensorCalculator::VaspForces(); The supercell structure has not been initialized yet.");
  }

  double cells = _supercell.getSupercellStructure().atoms.size() / _supercell.getInputStructure().atoms.size();

  xvector<double> testVec(3);
  _logger << "Managing directories for Thermal Conductivity." << apl::endl;
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

  for (uint ipa = 0; ipa < _strpair.ipairs.size(); ipa++) {
    if (_strpair.ipairs[ipa].itriplets_exist) {
      uint at1 = _strpair.pairs[_strpair.ipairs[ipa].eq[0]].indexA * cells;
      uint at2 = _strpair.pairs[_strpair.ipairs[ipa].eq[0]].indexB;
      for (uint i = 0; i < 6; i++) {
        for (uint j = 0; j < 6; j++) {
          if (_strpair.ipairs[ipa].dist_in[i][j]) {
            // Copy settings from common case
            xInputsT.push_back(_xInput);
            int idxRun = xInputsT.size() - 1;
            // Create run id name
            string runname = "./" + _AFLOW_APL_AFLOW_DIRECTORY_PREFIX_ + "AAPL_T" + stringify(idxRun) + "_A1_" + stringify(at1) + "_A2_" + stringify(at2) + "_D1_" + stringify(i) + "_D2_" + stringify(j);

            // Setup working directory
            xInputsT[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname);

            // Generate POSCAR and apply the unique distortion for one inequvalent atom
            // This distortion vector is stored in Cartesan form, hence use C2F before applying
            //CO - START
            xInputsT[idxRun].setXStr(_supercell.getSupercellStructureLight());  
            //CO - END
            xInputsT[idxRun].getXStr().atoms[at1].fpos = xInputsT[idxRun].getXStr().atoms[at1].fpos + C2F(xInputsT[idxRun].getXStr().lattice, DISTORTION_MAGNITUDE * ortoDistorsions[i]);
            xInputsT[idxRun].getXStr().atoms[at1].cpos = F2C(xInputsT[idxRun].getXStr().lattice,
                                                           xInputsT[idxRun].getXStr().atoms[at1].fpos);
            if (at1 != at2 || i != j) {
              xInputsT[idxRun].getXStr().atoms[at2].fpos = xInputsT[idxRun].getXStr().atoms[at2].fpos + C2F(xInputsT[idxRun].getXStr().lattice, DISTORTION_MAGNITUDE * ortoDistorsions[j]);
              xInputsT[idxRun].getXStr().atoms[at2].cpos = F2C(xInputsT[idxRun].getXStr().lattice,
                                                             xInputsT[idxRun].getXStr().atoms[at2].fpos);
            }

            //CO - START
            // If there are already LOCK or vasprun.xml.static files, it means this directory was already generated and computed,
            // hence do not touch and leave, but store this structure in the
            // list, hence it will be used in next part of code.
            if(aurostd::FileExist(xInputsT[idxRun].getDirectory() + string("/") + _AFLOWIN_)){continue;} //CO 180406
            if(xInputsT[idxRun].AFLOW_MODE_VASP){
            if (aurostd::FileExist(xInputsT[idxRun].getDirectory() + string("/") + _AFLOWLOCK_) ||
                aurostd::EFileExist(xInputsT[idxRun].getDirectory() + string("/OUTCAR.static")) ||
                aurostd::EFileExist(xInputsT[idxRun].getDirectory() + string("/vasprun.xml.static"))){continue;}
            }
            // If not, continue in this way, prepare generation of aflow.in ...

            _logger << "Creating " << xInputsT[idxRun].getDirectory() << apl::endl;
            //CO - END

            // Switch off autotune, because....
            _kbinFlags.KBIN_MPI_AUTOTUNE = true;

            // Common KPOINTS settings and OVERRIDES
            if(xInputsT[idxRun].AFLOW_MODE_VASP){
              xInputsT[idxRun].xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
              if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry)
                xInputsT[idxRun].xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.content_string;
              xInputsT[idxRun].xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;
              if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry)
                xInputsT[idxRun].xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.content_uint;

            // Clear old INCAR and set it as we want...
              xInputsT[idxRun].xvasp.INCAR.str(std::string());
            string system;
              for (uint i = 0; i < xInputsT[idxRun].getXStr().species.size(); i++)
                system = system + xInputsT[idxRun].getXStr().species_pp.at(i);
            system = system + "@" + runname;
              xInputsT[idxRun].xvasp.INCAR << "SYSTEM=" << system << std::endl;
              //            xInputsT[idxRun].xvasp.INCAR << std::endl;
              xInputsT[idxRun].xvasp.INCAR << "# Added by [AFLOW_PHONONS] begin" << std::endl;
              xInputsT[idxRun].xvasp.INCAR << "NELMIN=4         # The forces have to be well converged" << std::endl;
              xInputsT[idxRun].xvasp.INCAR << "NELM = 120       # Many electronic steps (SC2013)" << std::endl;
              xInputsT[idxRun].xvasp.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
//              xInputsT[idxRun].xvasp.INCAR << "ENCUT= 500   # For finer forces" << std::endl; // ME 180509
//              xInputsT[idxRun].xvasp.INCAR << "EDIFF = 1E-6   # For finer forces" << std::endl; // ME 180509
              xInputsT[idxRun].xvasp.INCAR << "# Added by [AFLOW_PHONONS] end" << std::endl;

            // Change format of POSCAR
            if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
                (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos)))
                xInputsT[idxRun].getXStr().is_vasp5_poscar_format = false;
            }

            // Create aflow.in
            writeOUTPUT(xInputsT[idxRun]);
          }
        }
      }
    }
  }
}

//StoreTensor////////////////////////////////////////////////////////////////////////////
//Store forces from vasp calculations of inequivalent pairs and distortions

void PhononCalculator::StoreForcesTensor() {
  //---------------------------------------------------------------------------------------------//
  //WARNING: Here, there is a correction in phonons code for polar meterials. Must be considered...
  //---------------------------------------------------------------------------------------------//

  //CO - START
  //Extract all forces from vasprun.xml.static --------------------------------------------------------------
  //CO - END
  _logger << "Storing forces for Lattice Thermal Conductivity." << apl::endl;

  bool found_outfile_anywhere=false;
  for (uint idxRun = 0; idxRun < xInputsT.size() && !found_outfile_anywhere; idxRun++) {
    // If tarred and compressed directory exists...

    // COREY CHECK THIS
    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
    string tarfilename = xInputsT[idxRun].getDirectory() + ".tar";
    for(uint iext=0;iext<vext.size();iext++) {
      if (aurostd::FileExist(tarfilename+vext.at(iext))) {
	aurostd::execute("tar -xf " + tarfilename+vext.at(iext));
      }
    }
    
    // If the LOCK file is missing, then it is probably corrupted run and
    // do not accept it and wait for its new run
    //CO - START
    //if (!aurostd::FileExist(xInputsT[idxRun].getDirectory() + string("/" + _AFLOWLOCK_)))
    //  throw APLStageBreak();

    if(_kbinFlags.AFLOW_MODE_VASP){
      if(aurostd::EFileExist(xInputs[idxRun].getDirectory() + string("/vasprun.xml.static"))) {
        found_outfile_anywhere=true;
        break;
      }
    }
    if(_kbinFlags.AFLOW_MODE_AIMS){
      if(aurostd::EFileExist(xInputs[idxRun].getDirectory() + string("/aims.out"))) {
        found_outfile_anywhere=true;
        break;
      }
    }
  }
  if(!found_outfile_anywhere){throw APLStageBreak();}

  //second pass, make sure it's everywhere!
  for (uint idxRun = 0; idxRun < xInputsT.size(); idxRun++) {
    string tarfilename = xInputsT[idxRun].getDirectory() + ".tar";
    xInputsT[idxRun].getXStr().qm_forces.clear();
    // Load data....
    if(_kbinFlags.AFLOW_MODE_VASP){
      if(!aurostd::EFileExist(xInputsT[idxRun].getDirectory() + string("/vasprun.xml.static"))) {
        _logger << apl::warning << "The vasprun.xml.static file in " << xInputsT[idxRun].getDirectory() << " directory is missing." << apl::endl;
        throw APLRuntimeError("apl::DirectMethodPC::runVASPCalculations(); Missing data from one job.");
      }
      xVASPRUNXML vasprunxml(xInputsT[idxRun].getDirectory() + "/vasprun.xml.static");
      for (uint i = 0; i < vasprunxml.vforces.size(); i++) xInputsT[idxRun].getXStr().qm_forces.push_back(vasprunxml.vforces.at(i));
      if (int(xInputsT[idxRun].getXStr().qm_forces.size()) == _supercell.getNumberOfAtoms()) { xInputsT[idxRun].getXStr().qm_calculated = TRUE; }
    }
    if(_kbinFlags.AFLOW_MODE_AIMS){
      if(!aurostd::EFileExist(xInputsT[idxRun].getDirectory() + string("/aims.out"))) {
        _logger << apl::warning << "The aims.out file in " << xInputsT[idxRun].getDirectory() << " directory is missing." << apl::endl;
        throw APLRuntimeError("apl::DirectMethodPC::runAIMSCalculations(); Missing data from one job.");
      }
      xAIMSOUT xaimsout(xInputsT[idxRun].getDirectory() + "/aims.out");
      for (uint i = 0; i < xaimsout.vforces.size(); i++) xInputsT[idxRun].getXStr().qm_forces.push_back(xaimsout.vforces.at(i));
      if (int(xInputsT[idxRun].getXStr().qm_forces.size()) == _supercell.getNumberOfAtoms()) { xInputsT[idxRun].getXStr().qm_calculated = TRUE; }
    }

    // Was it all right?
    if (!xInputsT[idxRun].getXStr().qm_calculated) {
      if(_kbinFlags.AFLOW_MODE_VASP){
        _logger << apl::warning << "The vasprun.xml.static file in " << xInputs[idxRun].getDirectory() << " is wrong." << apl::endl;
	throw APLRuntimeError("apl::TensorCalculator::VaspForces(); Missing data from one job.");
      }
      if(_kbinFlags.AFLOW_MODE_AIMS){
        _logger << apl::warning << "The aims.out file in " << xInputs[idxRun].getDirectory() << " is wrong." << apl::endl;
        throw APLRuntimeError("apl::TensorCalculator::AIMSForces(); Missing data from one job.");
      }
    }

    // Pack/Remove the whole directory...
    if (DOtar) {
      if (!aurostd::FileExist(tarfilename)) {
	aurostd::execute(string("tar -cf ") + tarfilename + " " + xInputsT[idxRun].getDirectory() + "/");
	aurostd::CompressFile(tarfilename,_kbinFlags.KZIP_BIN);
      }
      if (aurostd::EFileExist(tarfilename)) aurostd::execute(string("rm -rf ") + xInputsT[idxRun].getDirectory() + "/");
    }
  }

  // Remove Zero Forces...
  if (_calculateZeroStateForces) {
    for (uint idxRun = 0; idxRun < xInputsT.size(); idxRun++) {
      for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
        xInputsT[idxRun].getXStr().qm_forces[k](1) = xInputsT[idxRun].getXStr().qm_forces[k](1) - xInputs.back().getXStr().qm_forces[k](1);
        xInputsT[idxRun].getXStr().qm_forces[k](2) = xInputsT[idxRun].getXStr().qm_forces[k](2) - xInputs.back().getXStr().qm_forces[k](2);
        xInputsT[idxRun].getXStr().qm_forces[k](3) = xInputsT[idxRun].getXStr().qm_forces[k](3) - xInputs.back().getXStr().qm_forces[k](3);
      }
    }
  }

  //Create big tensor-----------------------------------

  for (uint i = 0; i < _supercell.getInputStructure().atoms.size(); i++) {
    vector<vector<vector<xvector<double> > > > pForces;
    for (uint j = 0; j < _supercell.getSupercellStructure().atoms.size(); j++) {
      vector<vector<xvector<double> > > atomsForces;
      for (uint k = 0; k < _supercell.getSupercellStructure().atoms.size(); k++) {
        vector<xvector<double> > bigmatrix;
        for (uint l = 0; l < 36; l++) {
          xvector<double> force(3);
          force(1) = 0;
          force(2) = 0;
          force(3) = 0;
          bigmatrix.push_back(force);
        }
        atomsForces.push_back(bigmatrix);
      }
      pForces.push_back(atomsForces);
    }
    bigTensor.push_back(pForces);
  }

  // Store forces inequivalent pairs------------------------------------------------------------------
  int idxRun = 0;
  for (uint ipa = 0; ipa < _strpair.ipairs.size(); ipa++) {
    if (_strpair.ipairs[ipa].itriplets_exist) {
      uint at1 = _strpair.pairs[_strpair.ipairs[ipa].eq[0]].indexA;
      uint at1sc = _supercell.pc2scMap(at1);
      uint at2 = _strpair.pairs[_strpair.ipairs[ipa].eq[0]].indexB;
      for (uint i = 0; i < 6; i++) {
        for (uint j = 0; j < 6; j++) {
          if (_strpair.ipairs[ipa].dist_in[i][j]) {
            for (int at3 = 0; at3 < _supercell.getNumberOfAtoms(); at3++) {
              bigTensor[at1][at2][at3][(i * 6) + j][1] = xInputsT[idxRun].getXStr().qm_forces[at3](1);
              bigTensor[at1][at2][at3][(i * 6) + j][2] = xInputsT[idxRun].getXStr().qm_forces[at3](2);
              bigTensor[at1][at2][at3][(i * 6) + j][3] = xInputsT[idxRun].getXStr().qm_forces[at3](3);
            }
            idxRun++;
          } else {
            if ((uint)_strpair.ipairs[ipa].dist[i][j] == (i * 6) + j) {
              for (int at3 = 0; at3 < _supercell.getNumberOfAtoms(); at3++) {
                bigTensor[at1][at2][at3][(i * 6) + j][1] = 0.0;
                bigTensor[at1][at2][at3][(i * 6) + j][2] = 0.0;
                bigTensor[at1][at2][at3][(i * 6) + j][3] = 0.0;
              }
            } else {
              int SymU = _strpair.getU(_supercell.getSupercellStructure(), ipa, at1sc, at2, i, j);
              for (int at3 = 0; at3 < _supercell.getNumberOfAtoms(); at3++) {
                int at4 = _strpair.getAtom(_supercell.getSupercellStructure(), at3, SymU);
                bigTensor[at1][at2][at3][(i * 6) + j] = _supercell.getSupercellStructure().fgroup[SymU].Uc * bigTensor[at1][at2][at4][_strpair.ipairs[ipa].dist[i][j]];
              }
            }
          }
        }
      }
    }
  }

  // Clear useless staff
  for (uint i = 0; i < xInputsT.size(); i++) {
    xInputsT[i].clear();
  }
  xInputsT.clear();
  for (uint i = 0; i < xInputs.size(); i++) {
    xInputs[i].clear();
  }
  xInputs.clear();
}

//buildTensor////////////////////////////////////////////////////////////////////////////

void PhononCalculator::buildTensor() {
  _logger << "Computing Anharmonic IFCs." << apl::endl;

  //Create small tensor-----------------------------------

  for (uint i = 0; i < _supercell.getInputStructure().atoms.size(); i++) {
    vector<vector<vector<vector<xvector<double> > > > > pForces;
    for (uint j = 0; j < _supercell.getSupercellStructure().atoms.size(); j++) {
      vector<vector<vector<xvector<double> > > > atomsForces;
      for (uint k = 0; k < _supercell.getSupercellStructure().atoms.size(); k++) {
        vector<vector<xvector<double> > > dist1;
        for (uint l = 0; l < 3; l++) {
          vector<xvector<double> > dist2;
          for (uint m = 0; m < 3; m++) {
            xvector<double> force(3);
            force(1) = 0.0;
            force(2) = 0.0;
            force(3) = 0.0;
            dist2.push_back(force);
          }
          dist1.push_back(dist2);
        }
        atomsForces.push_back(dist1);
      }
      pForces.push_back(atomsForces);
    }
    smallTensor_sumed.push_back(pForces);
    smallTensor.push_back(pForces);
  }

  // Computing Inequivalent triplets

  for (uint pp = 0; pp < _strpair.itriplets.size(); pp++) {
    _strpair.itriplets[pp].val_independent.clear();
    _strpair.itriplets[pp].val_independent.resize(_strpair.itriplets[pp].rot_independent.size(), 0.);
    _strpair.itriplets[pp].valtmp_independent.clear();
    _strpair.itriplets[pp].valtmp_independent.resize(_strpair.itriplets[pp].rot_independent.size(), 0.);
    uint at1sc = _strpair.itriplets[pp].triplets[0][0];
    uint at1pc = _supercell.sc2pcMap(at1sc);
    uint at2sc = _strpair.itriplets[pp].triplets[0][1];
    uint at3sc = _strpair.itriplets[pp].triplets[0][2];
    for (uint qq = 0; qq < _strpair.itriplets[pp].rot_independent.size(); qq++) {
      uint glob_index = _strpair.itriplets[pp].rot_independent[qq];
      uint k = glob_index % 3;
      uint aux1 = glob_index / 3;
      uint j = aux1 % 3;
      uint i = aux1 / 3;
      if (at1sc == at2sc && i == j) {  // Special case
        _strpair.itriplets[pp].val_independent[qq] = -(bigTensor[at1pc][at2sc][at3sc][(i * 6) + j][k + 1] + bigTensor[at1pc][at2sc][at3sc][(i * 6) + j + 21][k + 1]) / (DISTORTION_MAGNITUDE * DISTORTION_MAGNITUDE);
      } else {
        _strpair.itriplets[pp].val_independent[qq] = (bigTensor[at1pc][at2sc][at3sc][(i * 6) + j + 3][k + 1] + bigTensor[at1pc][at2sc][at3sc][(i * 6) + j + 18][k + 1] - bigTensor[at1pc][at2sc][at3sc][(i * 6) + j][k + 1] - bigTensor[at1pc][at2sc][at3sc][(i * 6) + j + 21][k + 1]) / (4 * DISTORTION_MAGNITUDE * DISTORTION_MAGNITUDE);
      }
    }
  }

  // Populating smallTensor

  for (uint pp = 0; pp < _strpair.itriplets.size(); pp++) {
    for (uint qq = 0; qq < _strpair.itriplets[pp].triplets.size(); qq++) {
      uint at1sc = _strpair.itriplets[pp].triplets[qq][0];
      uint at1pc = _supercell.sc2pcMap(at1sc);
      uint at2sc = _strpair.itriplets[pp].triplets[qq][1];
      uint at3sc = _strpair.itriplets[pp].triplets[qq][2];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            uint aux1 = (i * 9) + (j * 3) + k;
            for (uint rr = 0; rr < _strpair.itriplets[pp].rot_independent.size(); rr++) {
              smallTensor[at1pc][at2sc][at3sc][i][j][k + 1] += _strpair.itriplets[pp].rot_transform_array[qq][aux1][rr] * _strpair.itriplets[pp].val_independent[rr];
            }
          }
        }
      }
    }
  }

  // Counting how many times each inequivalent triplet is used

  for (uint pp = 0; pp < _strpair.itriplets.size(); pp++) {
    _strpair.itriplets[pp].count_independent.clear();
    _strpair.itriplets[pp].count_independent.resize(_strpair.itriplets[pp].rot_independent.size(), 0);
    for (uint jj = 0; jj < _strpair.itriplets[pp].rot_independent.size(); jj++) {
      for (uint qq = 0; qq < _strpair.itriplets[pp].triplets.size(); qq++) {
        for (int ii = 0; ii < 27; ii++) {
          if (abs(_strpair.itriplets[pp].rot_transform_array[qq][ii][jj]) > 1e-10) {
            _strpair.itriplets[pp].count_independent[jj] += 1;
          }
        }
      }
    }
  }
}

void PhononCalculator::MIX() {
  for (uint pp = 0; pp < _strpair.itriplets.size(); pp++) {
    for (uint jj = 0; jj < _strpair.itriplets[pp].rot_independent.size(); jj++) {
      _strpair.itriplets[pp].valtmp_independent[jj] = 0.;
      int count = 0;
      for (uint qq = 0; qq < _strpair.itriplets[pp].triplets.size(); qq++) {
        uint at1sc = _strpair.itriplets[pp].triplets[qq][0];
        uint at1pc = _supercell.sc2pcMap(at1sc);
        uint at2sc = _strpair.itriplets[pp].triplets[qq][1];
        uint at3sc = _strpair.itriplets[pp].triplets[qq][2];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              uint ii = (i * 9) + (j * 3) + k;
              if (abs(_strpair.itriplets[pp].rot_transform_array[qq][ii][jj]) > 1e-10) {
                if (smallTensor[at1pc][at2sc][at3sc][i][j][k + 1] != 0.) {
                  count = count + 1;
                  _strpair.itriplets[pp].valtmp_independent[jj] += (smallTensor_sumed[at1pc][at2sc][at3sc][i][j][k + 1] / smallTensor[at1pc][at2sc][at3sc][i][j][k + 1]) * _strpair.itriplets[pp].val_independent[jj];
                }
              }
            }
          }
        }
      }
      if (count != 0) {
        _strpair.itriplets[pp].val_independent[jj] = _strpair.itriplets[pp].valtmp_independent[jj] / count;
      }
    }
  }
}

void PhononCalculator::RESYM() {
  // Populating smallTensor
  for (uint pp = 0; pp < _strpair.itriplets.size(); pp++) {
    for (uint qq = 0; qq < _strpair.itriplets[pp].triplets.size(); qq++) {
      uint at1sc = _strpair.itriplets[pp].triplets[qq][0];
      uint at1pc = _supercell.sc2pcMap(at1sc);
      uint at2sc = _strpair.itriplets[pp].triplets[qq][1];
      uint at3sc = _strpair.itriplets[pp].triplets[qq][2];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            uint aux1 = (i * 9) + (j * 3) + k;
            smallTensor[at1pc][at2sc][at3sc][i][j][k + 1] = 0.;
            for (uint rr = 0; rr < _strpair.itriplets[pp].rot_independent.size(); rr++) {
              smallTensor[at1pc][at2sc][at3sc][i][j][k + 1] += _strpair.itriplets[pp].rot_transform_array[qq][aux1][rr] * _strpair.itriplets[pp].val_independent[rr];
            }
          }
        }
      }
    }
  }
}

double PhononCalculator::CHECK() {
  xvector<double> kk;
  double pp = 0;
  double maxR = 0;
  double counter = 0;
  for (uint ipa = 0; ipa < _strpair.pairs.size(); ipa++) {
    uint at1 = _strpair.pairs[ipa].indexA;
    uint at2 = _strpair.pairs[ipa].indexB;
    for (uint i = 0; i < 3; i++) {
      for (uint j = 0; j < 3; j++) {
        kk[1] = 0.0;
        kk[2] = 0.0;
        kk[3] = 0.0;
        for (int at3 = 0; at3 < _supercell.getNumberOfAtoms(); at3++) {
          kk = kk + smallTensor[at1][at2][at3][i][j];
        }
        counter = counter + 1;
        maxR = aurostd::max(maxR, abs(kk[1]));
        maxR = aurostd::max(maxR, abs(kk[2]));
        maxR = aurostd::max(maxR, abs(kk[3]));
        pp = pp + abs(kk[1]) + abs(kk[2]) + abs(kk[3]);
      }
    }
  }

  return maxR;
}

void PhononCalculator::CHECK2() {
  xvector<double> kk;
  double pp = 0;
  double maxR = 0;
  double counter = 0;
  for (uint ipa = 0; ipa < _strpair.pairs.size(); ipa++) {
    uint at1 = _strpair.pairs[ipa].indexA;
    uint at2 = _strpair.pairs[ipa].indexB;
    for (uint i = 0; i < 3; i++) {
      for (uint j = 0; j < 3; j++) {
        kk[1] = 0.0;
        kk[2] = 0.0;
        kk[3] = 0.0;
        for (int at3 = 0; at3 < _supercell.getNumberOfAtoms(); at3++) {
          kk = kk + smallTensor[at1][at2][at3][i][j];
        }
        counter = counter + 1;
        maxR = aurostd::max(maxR, abs(kk[1]));
        maxR = aurostd::max(maxR, abs(kk[2]));
        maxR = aurostd::max(maxR, abs(kk[3]));
        pp = pp + abs(kk[1]) + abs(kk[2]) + abs(kk[3]);
      }
    }
  }
}

void PhononCalculator::SumRuleTensor() {
  //Weight correction
  for (uint ipa = 0; ipa < _strpair.pairs.size(); ipa++) {
    uint i = _strpair.pairs[ipa].indexA;
    uint j = _strpair.pairs[ipa].indexB;
    for (uint x = 0; x < 3; x++) {
      for (uint y = 0; y < 3; y++) {
        xvector<double> kkkk;
        xvector<double> kkkk2;
        kkkk2[1] = 0.0;
        kkkk2[2] = 0.0;
        kkkk2[3] = 0.0;
        kkkk[1] = 0.0;
        kkkk[2] = 0.0;
        kkkk[3] = 0.0;
        for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
          for (uint z = 1; z < 4; z++) {
            if (abs(smallTensor[i][j][k][x][y][z]) != 0.) {
              kkkk[z] = kkkk[z] + smallTensor[i][j][k][x][y][z];
              kkkk2[z] = kkkk2[z] + abs(smallTensor[i][j][k][x][y][z]);
            }
          }
        }
        for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
          for (uint z = 1; z < 4; z++) {
            if (abs(smallTensor[i][j][k][x][y][z]) != 0.) {
              smallTensor_sumed[i][j][k][x][y][z] = smallTensor[i][j][k][x][y][z] - (abs(smallTensor[i][j][k][x][y][z]) / kkkk2[z]) * kkkk[z];
            }
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Prepare calculations for Tensor/////

double PhononCalculator::getElementTensor(uint a, uint b, uint c, uint d, uint e, uint f) {
  return smallTensor[a][b][c][d][e][f];
}

//////////////////////////////////////////////////////////////////////////////

void PhononCalculator::PrintTensor(string filename) {
  //CO - START
  //ofstream outfile(filename,ios_base::out);
  stringstream outfile;
  //if( !outfile.is_open() ) {
  //  throw apl::APLRuntimeError("Cannot open output TENSOR file.");
  //}
  //CO - END
  for (uint ipa = 0; ipa < _strpair.pairs.size(); ipa++) {
    uint at1 = _strpair.pairs[ipa].indexA;
    uint at2 = _strpair.pairs[ipa].indexB;
    for (uint t = 0; t < _strpair.pairs[ipa].tri.size(); t++) {
      uint at3 = _strpair.pairs[ipa].tri[t];
      outfile << at1 << " " << at2 << " " << at3 << std::endl;
      for (uint i = 0; i < 3; i++) {
        for (uint j = 0; j < 3; j++) {
          for (uint k = 0; k < 3; k++) {
            outfile << i << " " << j << " " << k << " " << smallTensor[at1][at2][at3][i][j][k + 1] << std::endl;
          }
        }
      }
    }
  }
  //CO - START
  //outfile.clear();
  //outfile.close();
  aurostd::stringstream2file(outfile, filename);
  if (!aurostd::FileExist(filename)) {
    throw apl::APLRuntimeError("Cannot open output TENSOR file.");
  }
  //CO - END
}

void PhononCalculator::PrintTensorBTE(string filename) {
  //CO - START
  //ofstream outfile(filename,ios_base::out);
  //if( !outfile.is_open() ) {
  //  throw apl::APLRuntimeError("Cannot open output TENSOR file.");
  //}
  stringstream outfile;
  //CO - END

  int half = _supercell.getSupercellStructure().atoms.size() / _supercell.getInputStructure().atoms.size();
  xvector<double> vec1, vec2, at1, at2, at3;

  xmatrix<double> lattice = _supercell.getSupercellStructure().lattice;

  int counter = 1;
  for (int i = 0; i < int(_supercell.getInputStructure().atoms.size()); i++) {
    uint at1sc = _supercell.pc2scMap(i);
    at1 = _supercell.getSupercellStructure().atoms[at1sc].cpos;
    for (int j = 0; j < half; j++) {
      vec1 = _supercell.getSupercellStructure().atoms[j].fpos;
      for (int r = 1; r < 4; r++) {
        if (vec1[r] > 0.5000000) {
          vec1[r] = vec1[r] - 1;
        }
      }
      vec1 = F2C(_supercell.getSupercellStructure().lattice, vec1);
      for (int j2 = 0; j2 < int(_supercell.getInputStructure().atoms.size()); j2++) {
        at2 = _supercell.getSupercellStructure().atoms[j + j2 * half].cpos;
        double dist1 = _strpair.getMinDist(lattice, at1, at2);
        if (dist1 < _strpair.cutoff) {
          for (int k = 0; k < half; k++) {
            vec2 = _supercell.getSupercellStructure().atoms[k].fpos;
            for (int r = 1; r < 4; r++) {
              if (vec2[r] > 0.5000000) {
                vec2[r] = vec2[r] - 1;
              }
            }
            vec2 = F2C(_supercell.getSupercellStructure().lattice, vec2);
            for (int k2 = 0; k2 < int(_supercell.getInputStructure().atoms.size()); k2++) {
              at3 = _supercell.getSupercellStructure().atoms[k + k2 * half].cpos;
              double dist2 = _strpair.getMinDist(lattice, at1, at3);
              double dist3 = _strpair.getMinDist(lattice, at3, at2);
              if (dist2 < _strpair.cutoff && dist3 < _strpair.cutoff) {
                outfile << counter << std::endl;
                outfile << vec1 << std::endl;
                outfile << vec2 << std::endl;
                outfile << i + 1 << "  " << j2 + 1 << "  " << k2 + 1 << std::endl;
                for (int o = 0; o < 3; o++) {
                  for (int p = 0; p < 3; p++) {
                    for (int q = 1; q < 4; q++) {
                      outfile << o + 1 << " " << p + 1 << "  " << q << "  " << std::scientific << smallTensor[i][j + j2 * half][k + k2 * half][o][p][q] << std::endl;
                    }
                  }
                }
                counter = counter + 1;
                outfile << std::endl;
              }
            }
          }
        }
      }
    }
  }
  //CO - START
  //outfile.clear();
  //outfile.close();
  aurostd::stringstream2file(outfile, filename);
  if (!aurostd::FileExist(filename)) {
    throw apl::APLRuntimeError("Cannot open output TENSOR file.");
  }
  //CO - END

  // Create supecel in PHONOPY style

  vector<xvector<double> > list_vec;
  vector<int> list_index;

  for (uint at = 0; at < _supercell.getInputStructure().atoms.size(); at++) {
    xvector<double> kk;
    kk = _supercell.getInputStructure().atoms[at].fpos;
    for (uint zz = 0; zz < 4; zz++) {
      kk[3] = (_supercell.getInputStructure().atoms[at].fpos[3] / 4.0) + double(zz) / 4.0;
      for (uint yy = 0; yy < 4; yy++) {
        kk[2] = (_supercell.getInputStructure().atoms[at].fpos[2] / 4.0) + double(yy) / 4.0;
        for (uint xx = 0; xx < 4; xx++) {
          kk[1] = (_supercell.getInputStructure().atoms[at].fpos[1] / 4.0) + double(xx) / 4.0;
          list_vec.push_back(kk);
          list_index.push_back(0);
        }
      }
    }
  }

  // Compare PHONOPY SUPERCELL and APL SUPERCELL

  for (uint at = 0; at < _supercell.getSupercellStructure().atoms.size(); at++) {
    for (uint at2 = 0; at2 < _supercell.getSupercellStructure().atoms.size(); at2++) {
      if (aurostd::modulus(_supercell.getSupercellStructure().atoms[at2].fpos - list_vec[at]) < _AFLOW_APL_EPS_) {
        list_index[at] = at2;
        at2 = _supercell.getSupercellStructure().atoms.size();
      }
    }
  }

  //CO - START
  //ofstream outfile2("FORCE_CONSTANTS_2ND",ios_base::out);
  stringstream outfile2;
  outfile2 << std::fixed << std::showpoint;
  outfile2 << std::setprecision(15);
  //if( !outfile2.is_open() ) {
  //  throw apl::APLRuntimeError("Cannot open output TENSOR file.");
  //}
  //CO - END

  uint scAtomsSize = _supercell.getSupercellStructure().atoms.size();
  //uint pcAtomsSize = _supercell.getInputStructure().atoms.size();
  outfile2 << std::fixed << std::showpoint << std::setprecision(0) << _supercell.getSupercellStructure().atoms.size() << std::endl;
  for (uint isc1 = 0; isc1 < scAtomsSize; isc1++) {
    for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
      outfile2 << std::fixed << std::showpoint << std::setprecision(0) << isc1 + 1 << "    " << isc2 + 1 << std::endl;
      for (uint ipc1 = 1; ipc1 < 4; ipc1++) {
        for (uint ipc2 = 1; ipc2 < 4; ipc2++) {
          outfile2 << std::fixed << std::showpoint << std::setprecision(15) << -_forceConstantMatrices[list_index[isc1]][list_index[isc2]][ipc1][ipc2] << " ";
        }
        outfile2 << std::endl;
      }
      outfile2 << std::endl;
    }
  }

  //CO - START
  string filename2 = "FORCE_CONSTANTS_2ND";
  aurostd::stringstream2file(outfile2, filename2);
  if (!aurostd::FileExist(filename2)) {
    throw apl::APLRuntimeError("Cannot open output TENSOR file.");
  }

  //outfile2.clear();
  //outfile2.close();
  //CO - END
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Y. Wang et.al, J. Phys.:Condens. Matter 22, 202201 (2010)
// DOI: 10.1088/0953-8984/22/20/202201

xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWangAAPL(const xvector<double> _q, vector<xmatrix<xcomplex<double> > >& kk) {
  const xstructure& sc = _supercell.getSupercellStructure();
  const xstructure& pc = _supercell.getInputStructure();

  // to correct the q=\Gamma as a limit
  xvector<double> q(_q);
  if (aurostd::modulus(q) < _AFLOW_APL_EPS_) { q(1) = _AFLOW_APL_EPS_ * 1.001; }

  uint pcAtomsSize = pc.atoms.size();
  uint nBranches = 3 * pcAtomsSize;

  xmatrix<xcomplex<double> > dynamicMatrix(nBranches, nBranches);

  // Calculation
  double fac0 = 13.605826 * 2.0 * 0.529177249;  // from a.u. to eV/A
  double volume = det(pc.lattice);
  double fac1 = 4.0 * PI / volume;
  double nbCells = det(sc.lattice) / volume;

  if (aurostd::modulus(q) > _AFLOW_APL_EPS_) {
    double fac2 = fac0 * fac1;
    double fac3 = 1.0 / scalar_product(q, _dielectricTensor * q) / nbCells;

    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++)
      for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            int typei = pc.atoms[ipc1].type;
            int typej = pc.atoms[ipc2].type;
            xcomplex<double> tmp1 = (q * _bornEffectiveChargeTensor[typei])(ix);
            xcomplex<double> tmp2 = (q * _bornEffectiveChargeTensor[typej])(iy);
            dynamicMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) = fac2 * fac3 * tmp1 * tmp2;
            kk[1](3 * ipc1 + ix, 3 * ipc2 + iy) = tmp1 * _bornEffectiveChargeTensor[typej](iy, 1) + tmp2 * _bornEffectiveChargeTensor[typei](ix, 1);
            -(2.0 * tmp1 * tmp2 * scalar_product(_dielectricTensor(1), q)) / fac3;
            kk[1](3 * ipc1 + ix, 3 * ipc2 + iy) = kk[1](3 * ipc1 + ix, 3 * ipc2 + iy) * fac2;
            kk[2](3 * ipc1 + ix, 3 * ipc2 + iy) = tmp1 * _bornEffectiveChargeTensor[typej](iy, 2) + tmp2 * _bornEffectiveChargeTensor[typei](ix, 2);
            -(2.0 * tmp1 * tmp2 * scalar_product(_dielectricTensor(2), q)) / fac3;
            kk[2](3 * ipc1 + ix, 3 * ipc2 + iy) = kk[2](3 * ipc1 + ix, 3 * ipc2 + iy) * fac2;
            kk[3](3 * ipc1 + ix, 3 * ipc2 + iy) = tmp1 * _bornEffectiveChargeTensor[typej](iy, 3) + tmp2 * _bornEffectiveChargeTensor[typei](ix, 3);
            -(2.0 * tmp1 * tmp2 * scalar_product(_dielectricTensor(3), q)) / fac3;
            kk[3](3 * ipc1 + ix, 3 * ipc2 + iy) = kk[3](3 * ipc1 + ix, 3 * ipc2 + iy) * fac2;
          }
      }
  }

  //
  return dynamicMatrix;
}

// ///////////////////////////////////////////////////////////////////////////

xmatrix<xcomplex<double> > PhononCalculator::getDynamicMatrixAAPL2(const xvector<double> kpoint, vector<xmatrix<xcomplex<double> > >& dDynM) {
  xstructure pc = _supercell.getInputStructure();
  uint sc_atoms = _supercell.getSupercellStructure().atoms.size();
  uint pcAtomsSize = pc.atoms.size();
  uint nBranches = 3 * pcAtomsSize;
  double volume = det(pc.lattice);
  xvector<int> scell = _supercell.scell;
  double tcells = scell(1) * scell(2) * scell(3);
  // Setup both lattices
  pc.FixLattices();
  xmatrix<double> _rlattice, _klattice;
  _rlattice = pc.lattice;
  _klattice = ReciprocalLattice(_rlattice);
  xmatrix<xcomplex<double> > dynamicMatrix(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dynamicMatrix_nac(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixX(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixY(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixZ(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixX_nac(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixY_nac(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixZ_nac(nBranches, nBranches, 1, 1);

  //Creating storage for total forces
  vector<vector<xmatrix<double> > > _forceConstantMatrices_total;
  vector<vector<xmatrix<double> > > _forceConstantMatrices_dipol;

  _forceConstantMatrices_total = _forceConstantMatrices;
  _forceConstantMatrices_dipol = _forceConstantMatrices;

  // Use the 1st BZ image of each q point to improve the behavior of the non-analytic correction.
  double tmp1, tmp2;
  xvector<double> shortest;
  shortest = kpoint;
  tmp1 = aurostd::modulus(shortest);
  for (int o = -2; o <= 2; o++) {
    for (int p = -2; p <= 2; p++) {
      for (int q = -2; q <= 2; q++) {
        xvector<double> pp = kpoint + 10.0 * (o * _klattice(1) + p * _klattice(2) + q * _klattice(3));
        tmp2 = aurostd::modulus(pp);
        if (tmp2 < tmp1) {
          tmp1 = tmp2;
          shortest = pp;
        }
      }
    }
  }

  double fac0 = 13.605826 * 2.0 * 0.529177249;  // from a.u. to eV/A
  double fac1 = 4.0 * PI;

  // Non analytical contribution
  if (_isPolarMaterial && aurostd::modulus(shortest) > _AFLOW_APL_EPS_) {
    double tmp3 = scalar_product(shortest, _dielectricTensor * shortest);
    for (int iatom1 = 0; iatom1 < int(pcAtomsSize); iatom1++) {
      for (int iatom2 = 0; iatom2 < int(pcAtomsSize); iatom2++) {
        for (int i = 1; i < 4; i++) {
          for (int j = 1; j < 4; j++) {
            int typei = pc.atoms[iatom1].type;
            int typej = pc.atoms[iatom2].type;
            double tmp1 = (shortest * _bornEffectiveChargeTensor[typei])(i);
            double tmp2 = (shortest * _bornEffectiveChargeTensor[typej])(j);
            dynamicMatrix_nac[(3 * iatom1) + i][(3 * iatom2) + j] = fac1 * fac0 * (tmp1 * tmp2) / (tmp3 * volume);
            dDynamicMatrixX_nac[(3 * iatom1) + i][(3 * iatom2) + j] = fac1 * fac0 * (tmp1 * _bornEffectiveChargeTensor[typej](j, 1) + tmp2 * _bornEffectiveChargeTensor[typei](i, 1) - 2 * tmp1 * tmp2 * scalar_product(_dielectricTensor(1), shortest) / tmp3) / (tmp3 * volume);
            dDynamicMatrixY_nac[(3 * iatom1) + i][(3 * iatom2) + j] = fac1 * fac0 * (tmp1 * _bornEffectiveChargeTensor[typej](j, 2) + tmp2 * _bornEffectiveChargeTensor[typei](i, 2) - 2 * tmp1 * tmp2 * scalar_product(_dielectricTensor(2), shortest) / tmp3) / (tmp3 * volume);
            dDynamicMatrixZ_nac[(3 * iatom1) + i][(3 * iatom2) + j] = fac1 * fac0 * (tmp1 * _bornEffectiveChargeTensor[typej](j, 3) + tmp2 * _bornEffectiveChargeTensor[typei](i, 3) - 2 * tmp1 * tmp2 * scalar_product(_dielectricTensor(3), shortest) / tmp3) / (tmp3 * volume);
          }
        }
      }
    }

    // Transform back to real space to obtain a correction to the short-range force constants.
    for (int at1 = 0; at1 < int(sc_atoms); at1++) {
      int at1pc = _supercell.sc2pcMap(at1);
      for (int at2 = 0; at2 < int(sc_atoms); at2++) {
        int at2pc = _supercell.sc2pcMap(at2);
        for (int i = 1; i < 4; i++) {
          for (int j = 1; j < 4; j++) {
            _forceConstantMatrices_dipol[at1][at2](i, j) = dynamicMatrix_nac[(3 * at1pc) + i][(3 * at2pc) + j].re / (tcells);
            // Force constants with long-range correction.
            _forceConstantMatrices_total[at1][at2](i, j) = _forceConstantMatrices[at1][at2](i, j) - _forceConstantMatrices_dipol[at1][at2](i, j);
          }
        }
      }
    }
  }

  //Dynamical matrix in SHENGBTE style
  xcomplex<double> iONE(0.0, 1.0);
  xcomplex<double> ztmp;
  xvector<double> rcell, r;
  xmatrix<double> lattvec = _supercell.getInputStructure().lattice;
  xmatrix<double> lattvecSc = _supercell.getSupercellStructure().lattice;
  int neq = 0;  //Warning uninizialize
  //int neq=0;
  vector<double> qr(27, 1);
  vector<xvector<double> > rr(27, 1);
  for (int iatom1 = 0; iatom1 < int(pcAtomsSize); iatom1++) {
    for (int iatom2 = 0; iatom2 < int(pcAtomsSize); iatom2++) {
      for (int ix1 = 0; ix1 < scell(1); ix1++) {
        for (int iy1 = 0; iy1 < scell(2); iy1++) {
          for (int iz1 = 0; iz1 < scell(3); iz1++) {
            int isatom1 = _supercell.pc2scMap(iatom1);
            //int isatom2; //Warnin uninizialized
            int isatom2 = 0;

            rcell = ix1 * lattvec(1) + iy1 * lattvec(2) + iz1 * lattvec(3);
            r = pc.atoms[iatom1].cpos - pc.atoms[iatom2].cpos + rcell;

            // Localize index
            xvector<double> isatom2_vec = C2F(lattvecSc, pc.atoms[iatom2].cpos - rcell);
            for (int tatoms = 0; tatoms < int(_supercell.getSupercellStructure().atoms.size()); tatoms++) {
              xvector<double> fdiff = isatom2_vec - _supercell.getSupercellStructure().atoms[tatoms].fpos;
              SYM::PBC(fdiff);
              if (aurostd::modulus(fdiff) < 1E-04) {
                isatom2 = tatoms;
                tatoms = _supercell.getSupercellStructure().atoms.size();
              }
            }

            double dmin = 10000;
            for (int ix2 = -2; ix2 < 3; ix2++) {
              for (int iy2 = -2; iy2 < 3; iy2++) {
                for (int iz2 = -2; iz2 < 3; iz2++) {
                  xvector<double> rl = ix2 * scell(1) * lattvec(1) + iy2 * scell(2) * lattvec(2) + iz2 * scell(3) * lattvec(3);
                  double Rnorm = aurostd::modulus(rl + r);
                  if (abs(Rnorm - dmin) > 1E-04) {
                    if (Rnorm < dmin) {
                      neq = 1;
                      dmin = Rnorm;
                      qr[neq] = scalar_product(kpoint, (rl + rcell));
                      rr[neq] = rl + rcell;
                    }
                  } else {
                    neq = neq + 1;
                    qr[neq] = scalar_product(kpoint, (rl + rcell));
                    rr[neq] = rl + rcell;
                  }
                }
              }
            }
            xcomplex<double> star = 0.;
            for (int ip = 1; ip < neq + 1; ip++) {
              ztmp = exp(iONE * (-qr[ip])) / (double)neq;
              star = star + ztmp;
              for (int i = 1; i < 4; i++) {
                for (int j = 1; j < 4; j++) {
                  dynamicMatrix(3 * iatom1 + i, 3 * iatom2 + j) = dynamicMatrix(3 * iatom1 + i, 3 * iatom2 + j) + ztmp * (-1.0) * _forceConstantMatrices_total[isatom1][isatom2](i, j);
                  dDynamicMatrixX(3 * iatom1 + i, 3 * iatom2 + j) = dDynamicMatrixX(3 * iatom1 + i, 3 * iatom2 + j) - iONE * rr[ip][1] * ztmp * (-1.0) * _forceConstantMatrices_total[isatom1][isatom2](i, j);
                  dDynamicMatrixY(3 * iatom1 + i, 3 * iatom2 + j) = dDynamicMatrixY(3 * iatom1 + i, 3 * iatom2 + j) - iONE * rr[ip][2] * ztmp * (-1.0) * _forceConstantMatrices_total[isatom1][isatom2](i, j);
                  dDynamicMatrixZ(3 * iatom1 + i, 3 * iatom2 + j) = dDynamicMatrixZ(3 * iatom1 + i, 3 * iatom2 + j) - iONE * rr[ip][3] * ztmp * (-1.0) * _forceConstantMatrices_total[isatom1][isatom2](i, j);
                }
              }
              if (_isPolarMaterial && aurostd::modulus(shortest) > _AFLOW_APL_EPS_) {
                for (int i = 1; i < 4; i++) {
                  for (int j = 1; j < 4; j++) {
                    dDynamicMatrixX(3 * iatom1 + i, 3 * iatom2 + j) = dDynamicMatrixX(3 * iatom1 + i, 3 * iatom2 + j) - (star * dDynamicMatrixX_nac(3 * iatom1 + i, 3 * iatom2 + j) / tcells);
                    dDynamicMatrixY(3 * iatom1 + i, 3 * iatom2 + j) = dDynamicMatrixY(3 * iatom1 + i, 3 * iatom2 + j) - (star * dDynamicMatrixY_nac(3 * iatom1 + i, 3 * iatom2 + j) / tcells);
                    dDynamicMatrixZ(3 * iatom1 + i, 3 * iatom2 + j) = dDynamicMatrixZ(3 * iatom1 + i, 3 * iatom2 + j) - (star * dDynamicMatrixZ_nac(3 * iatom1 + i, 3 * iatom2 + j) / tcells);
                  }
                }
              }
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
          dynamicMatrix(3 * i + ix, 3 * j + iy) *= 1.0 * 9648.53336213 / sqrt(mass_i * mass_j);  //THZ
          dDynamicMatrixX(3 * i + ix, 3 * j + iy) *= 1.0 * 9648.53336213 / sqrt(mass_i * mass_j);
          dDynamicMatrixY(3 * i + ix, 3 * j + iy) *= 1.0 * 9648.53336213 / sqrt(mass_i * mass_j);
          dDynamicMatrixZ(3 * i + ix, 3 * j + iy) *= 1.0 * 9648.53336213 / sqrt(mass_i * mass_j);
        }
      }
    }
  }

  // Put all together
  dDynM[0] = dDynamicMatrixX;
  dDynM[1] = dDynamicMatrixY;
  dDynM[2] = dDynamicMatrixZ;

  // Clear stuff
  _forceConstantMatrices_total.clear();
  _forceConstantMatrices_total.clear();
  return dynamicMatrix;
}
// ///////////////////////////////////////////////////////////////////////////

xmatrix<xcomplex<double> > PhononCalculator::getDynamicMatrixAAPL(const xvector<double> kpoint, vector<xmatrix<xcomplex<double> > >& dDynM) {
  uint scAtomsSize = _supercell.getSupercellStructure().atoms.size();
  uint pcAtomsSize = _supercell.getInputStructure().atoms.size();

  uint nBranches = 3 * pcAtomsSize;
  xmatrix<xcomplex<double> > dynamicMatrix(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixX(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixY(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dDynamicMatrixZ(nBranches, nBranches, 1, 1);
  xmatrix<xcomplex<double> > dynamicMatrix0(nBranches, nBranches, 1, 1);

  xcomplex<double> phase;
  xvector<xcomplex<double> > deriv;
  double value;

  // Calculate nonanalytical contribution
  xmatrix<xcomplex<double> > dynamicMatrixNA(nBranches, nBranches, 1, 1);
  vector<xmatrix<xcomplex<double> > > dDynM_NA;
  for (uint ix = 1; ix <= 3; ix++) {
    dDynM_NA.push_back(dynamicMatrixNA);
  }

  if (_isPolarMaterial)
    dynamicMatrixNA = getNonanalyticalTermWangAAPL(kpoint, dDynM_NA);

  // Loop over primitive cell
  for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
    uint isc1 = _supercell.pc2scMap(ipc1);

    for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
      uint ipc2 = _supercell.sc2pcMap(isc2);

      if (_supercell.calcShellPhaseFactorAAPL(isc2, isc1, kpoint, phase, deriv)) {
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            //value = smallMatrix[ipc1][isc2][ix-1][iy];
            value = 0.5 * (_forceConstantMatrices[isc1][isc2](ix, iy) +
                           _forceConstantMatrices[isc2][isc1](iy, ix));
            dynamicMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) -= value * phase;
            if (_isPolarMaterial) {
              dynamicMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
              dDynamicMatrixX(3 * ipc1 + ix, 3 * ipc2 + iy) += (dynamicMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) - value) * deriv[1] + dDynM_NA[1](3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
              dDynamicMatrixY(3 * ipc1 + ix, 3 * ipc2 + iy) += (dynamicMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) - value) * deriv[2] + dDynM_NA[2](3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
              dDynamicMatrixZ(3 * ipc1 + ix, 3 * ipc2 + iy) += (dynamicMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) - value) * deriv[3] + dDynM_NA[3](3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
            } else {
              dDynamicMatrixX(3 * ipc1 + ix, 3 * ipc2 + iy) += (-value) * deriv[1];
              dDynamicMatrixY(3 * ipc1 + ix, 3 * ipc2 + iy) += (-value) * deriv[2];
              dDynamicMatrixZ(3 * ipc1 + ix, 3 * ipc2 + iy) += (-value) * deriv[3];
            }
            dynamicMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) -= value;
            if (_isPolarMaterial)
              dynamicMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy);
          }
        }
      }
    }
  }

  // Subtract the sum of all "forces" from the central atom, this is like an automatic sum rule...
  for (uint i = 0; i < pcAtomsSize; i++) {
    for (uint j = 0; j < pcAtomsSize; j++) {
      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
        for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
          dynamicMatrix(3 * i + ix, 3 * i + iy) = dynamicMatrix(3 * i + ix, 3 * i + iy) - dynamicMatrix0(3 * i + ix, 3 * j + iy);
        }
      }
    }
  }

  // Make it hermitian
  for (uint i = 0; i <= pcAtomsSize - 1; i++) {
    for (uint j = 0; j <= i; j++) {
      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
        for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
          dynamicMatrix(3 * i + ix, 3 * j + iy) += conj(dynamicMatrix(3 * j + iy, 3 * i + ix));
          dynamicMatrix(3 * i + ix, 3 * j + iy) *= 0.5;
          dynamicMatrix(3 * j + iy, 3 * i + ix) = conj(dynamicMatrix(3 * i + ix, 3 * j + iy));
          dDynamicMatrixX(3 * i + ix, 3 * j + iy) += conj(dDynamicMatrixX(3 * j + iy, 3 * i + ix));
          dDynamicMatrixX(3 * i + ix, 3 * j + iy) *= 0.5;
          dDynamicMatrixX(3 * j + iy, 3 * i + ix) = conj(dDynamicMatrixX(3 * i + ix, 3 * j + iy));
          dDynamicMatrixY(3 * i + ix, 3 * j + iy) += conj(dDynamicMatrixY(3 * j + iy, 3 * i + ix));
          dDynamicMatrixY(3 * i + ix, 3 * j + iy) *= 0.5;
          dDynamicMatrixY(3 * j + iy, 3 * i + ix) = conj(dDynamicMatrixY(3 * i + ix, 3 * j + iy));
          dDynamicMatrixZ(3 * i + ix, 3 * j + iy) += conj(dDynamicMatrixZ(3 * j + iy, 3 * i + ix));
          dDynamicMatrixZ(3 * i + ix, 3 * j + iy) *= 0.5;
          dDynamicMatrixZ(3 * j + iy, 3 * i + ix) = conj(dDynamicMatrixZ(3 * i + ix, 3 * j + iy));
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
          dynamicMatrix(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
          dDynamicMatrixX(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
          dDynamicMatrixY(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
          dDynamicMatrixZ(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
        }
      }
    }
  }

  // Put all together
  dDynM[0] = dDynamicMatrixX;
  dDynM[1] = dDynamicMatrixY;
  dDynM[2] = dDynamicMatrixZ;

  //
  return dynamicMatrix;
}

// ///////////////////////////////////////////////////////////////////////////

xvector<double> PhononCalculator::getEigenvaluesAAPL(const xvector<double>& kpoint, xmatrix<xcomplex<double> >& unitaryMatrix, vector<xmatrix<xcomplex<double> > >& dDynM) {  // Modified JJPR
  xmatrix<xcomplex<double> > dynamicMatrix = getDynamicMatrixAAPL2(kpoint, dDynM);

  // Diagonalize
  xvector<double> eigenvalues(dynamicMatrix.rows, 1);

  apl::aplEigensystems e;
  e.eigen_calculation(dynamicMatrix, eigenvalues, unitaryMatrix, APL_MV_EIGEN_SORT_VAL_ASC);

  return eigenvalues;
}

// ///////////////////////////////////////////////////////////////////////////

xvector<double> PhononCalculator::getFrequencyAAPL(const xvector<double>& kpoint, xmatrix<xcomplex<double> >& unitaryMatrix, vector<xmatrix<xcomplex<double> > >& dDynM) {  // Modified by JJPR
  // Compute frequency(omega) from eigenvalues [in eV/A/A/atomic_mass_unit]
  xvector<double> omega = getEigenvaluesAAPL(kpoint, unitaryMatrix, dDynM);

  // Transform values to desired format
  for (_AFLOW_APL_REGISTER_ int i = omega.lrows; i <= omega.urows; i++) {
    if (omega(i) < 0.0) {
      omega(i) = -sqrt(-omega(i));
    } else {
      omega(i) = sqrt(omega(i));
    }
  }

  // Return
  return (omega);
}
// ///////////////////////////////////////////////////////////////////////////
}

// ***************************************************************************
// *                                                                         *
// *              AFlow JOSE J. PLATA -   Duke University 2017               *
// *                                                                         *
// ***************************************************************************
