// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                Aflow PINKU NATH - Duke University 2014-2017             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu
//this cpp file will create phonon and static directories to run VASP calculations

#include "aflow_apl.h"

namespace apl {
// ***************************************************************************************
//eos_calculate::eos_calculate(Supercell& sc, StrPairs& strpair, _xinput& xinput, //_xvasp& xvasp,
eos_calculate::eos_calculate(Supercell& sc, vector<ClusterSet>& clst,
                             _xinput& xinput, _aflags& aflags, _kflags& kflags,
                             _xflags& xflags, //_vflags& vflags,
                             string& AflowIn, StrPairs& strpair, Logger& l)
    : PhononCalculator(sc, clst, xinput, aflags, kflags, xflags, AflowIn, strpair, l) {
  GRUNEISEN = false;
  EOS = false;
  GP_dir_names.clear();
  PH_dir_names.clear();
  EOS_dir_names.clear();
  GP_volumes.clear();
  EOS_Volumes.clear();
}
eos_calculate::~eos_calculate() { this->clear(); }
void eos_calculate::clear() {}
// ***************************************************************************************
//Create directories when GRUNEISEN option is on
void eos_calculate::GP_RUN() {
  _logger << "Creating directories to calculate Gruneisen Parameter " << apl::endl;
  string APL_DIR = PhononCalculator::_xInput.getDirectory();
  vector<_xinput> GPxRuns;
  //const xstructure& pc = _supercell.getInputStructure();
  int idxRun = 0;
  double ep = 1e-6;
  for (double i = -GP_VOL_DISTORTION; i <= GP_VOL_DISTORTION; i += GP_VOL_DISTORTION) {
    GPxRuns.push_back(_xInput);
    idxRun = GPxRuns.size() - 1;
    // new scale factor calulation
    double scale = 0.00;
    GPxRuns[idxRun].setXStr(PhononCalculator::_supercell.getInputStructure());
    xmatrix<double> m(3, 3, 1, 1);
    scale = GPxRuns[idxRun].getXStr().scale;
    m = GPxRuns[idxRun].getXStr().lattice;
    GPxRuns[idxRun].getXStr().atoms = PhononCalculator::_supercell.getInputStructure().atoms;
    double det = determinant(m);
    double volume = scale * scale * scale * det;
    if (std::abs(i) <= ep) GP_volumes.push_back(volume);
    if (std::abs(i) <= ep) continue;

    string runname = "";
    string tag = NumberToString<double>(std::abs(i));
    if (std::abs(i - GP_VOL_DISTORTION) <= ep)
      runname = string(_GP_AFLOW_DIRECTORY_PREFIX_) + string("P") + string(tag);
    else
      runname = string(_GP_AFLOW_DIRECTORY_PREFIX_) + string("M") + string(tag);
    GPxRuns[idxRun].setDirectory(APL_DIR + "./" + runname);
    GP_dir_names.push_back(runname);

    double newscale = 0.0;
    newscale = volume + (i * volume) / 100.00;
    newscale = newscale / det;
    newscale = pow(newscale, 1. / 3.0);
    GPxRuns[idxRun].getXStr().scale = newscale;
    double newvolume = newscale * newscale * newscale * det;
    GP_volumes.push_back(newvolume);
    // new scale factor calulation end

    for (uint j = 0; j < GPxRuns[idxRun].getXStr().atoms.size(); j++)
      GPxRuns[idxRun].getXStr().atoms[j] = PhononCalculator::_supercell.getInputStructure().atoms[j];
    m.clear();
    if(aurostd::FileExist(GPxRuns[idxRun].getDirectory() + string("/") + string(_AFLOWIN_))){continue;} //CO 180406
    if(GPxRuns[idxRun].AFLOW_MODE_VASP){
      if (aurostd::FileExist(GPxRuns[idxRun].getDirectory() + string("/") + string(_AFLOWLOCK_)) ||
         aurostd::EFileExist(GPxRuns[idxRun].getDirectory() + string("/OUTCAR.static"))){continue;}  //CO
    }
    if(GPxRuns[idxRun].AFLOW_MODE_AIMS){
      GPxRuns[idxRun].xaims.CONTROL.str(std::string());
      KBIN::AIMS_Produce_CONTROL(GPxRuns[idxRun].xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);  //DEFAULT
      KBIN::AIMS_Modify_CONTROL(GPxRuns[idxRun].xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);            //DEFAULT
    }
    writeGPOUTPUT(GPxRuns[idxRun], true);
  }  //for loop end

  GPxRuns.clear();
}

void eos_calculate::writeGPOUTPUT(_xinput& xinput, const bool _isGP) {
  if(!( xinput.AFLOW_MODE_VASP || xinput.AFLOW_MODE_AIMS )) { 
    throw APLRuntimeError("apl::PhononCalculator:writeGPOUTPUT(); Input -> aflow.in conversion unknown.");
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

  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << setprecision(SetPrecision);

  //corey - at some point, fix alien mode for aims, for now omit!
  if(xinput.AFLOW_MODE_VASP){
  //search some key words from main aflow.in file
  string AFLOWword1 = "[VASP_POSCAR_MODE_EXPLICIT]START";
  string AFLOWword2 = "[VASP_POSCAR_MODE_EXPLICIT]STOP";
  string AFLOWword3_original = "[AFLOW_PHONONS]GRUNEISEN";  //CO 170601
  string AFLOWword3_new = "[AFLOW_QHA]CALC";                //CO 170601
  string AFLOWword4_original = "[AFLOW_PHONONS]EOS=";       //CO 170601
  string AFLOWword4_new = "[AFLOW_QHA]EOS=";                //CO 170601
  string AFLOWword5 = "[VASP_KPOINTS_FILE]KSCHEME=";
  string AFLOWword6 = "[VASP_KPOINTS_FILE]STATIC_KSCHEME=";
  string AFLOWword7 = "AFLOW_MODE_ZIP";

  //exclude all the following lines from main [aflow.in] and create a new [aflow.in]"
  vector<string> vlines;  //CO
  uint line_count = 0;    //CO
  string line;
  uint ENTRY = 0;
  aurostd::efile2vectorstring(_AFLOWIN_, vlines);  //CO
  //ifstream myfile (_AFLOWIN_.c_str()); //CO
  //if (myfile.is_open()) //CO
  if (vlines.size())  //CO
  {
    //while ( getline (myfile,line) ) //CO
    while (line_count < vlines.size())  //CO
    {
      line = vlines[line_count++];  //CO
      if (line == "") continue;
      if (line[0] == '#') continue;
      if (line.find(AFLOWword1) != std::string::npos) { ENTRY = 1; }
      if (ENTRY != 0) {
        if (line.find(AFLOWword2) != std::string::npos) { ENTRY = 0; }
      } else if ( (line.find(AFLOWword3_original) != std::string::npos) || (line.find(AFLOWword3_new) != std::string::npos) ) //CO 170601
        outfile << "[AFLOW_QHA]GRUNEISENSD=y" << std::endl; //CO 170601
      else if ( (line.find(AFLOWword4_original) != std::string::npos) || (line.find(AFLOWword4_new) != std::string::npos) ) { //CO 170601
        if (!_isGP) outfile << "[AFLOW_QHA]EOSSD=y" << std::endl; //CO 170601
      } else if (line.find(AFLOWword5) != std::string::npos) {
        EOS_KSCHEME = line;
        outfile << line << std::endl;
      } else if (line.find(AFLOWword6) != std::string::npos) {
        EOS_STATIC_KSCHEME = line;
        outfile << line << std::endl;
      } else if (line.find(AFLOWword7) != std::string::npos) {
        outfile << "#" << line << std::endl;
      } else if ((line.find("DCUSERPATH") != std::string::npos) && (!_isGP)) {
      } else if ( ((line.find("[AFLOW_APL]DC=") != std::string::npos) //CO 170601
            || (line.find("[AFLOW_QHA]DC=") != std::string::npos) //CO 170601
            || (line.find("[AFLOW_AAPL]DC=") != std::string::npos) //CO 170601
            || (line.find("[AFLOW_PHONONS]DC=") != std::string::npos)) //CO 170601
          && (!_isGP)) {
      } else if ( ((line.find("[AFLOW_APL]DOS=") != std::string::npos) || //CO 170601
            (line.find("[AFLOW_QHA]DOS=") != std::string::npos) || //CO 170601
            (line.find("[AFLOW_AAPL]DOS=") != std::string::npos) || //CO 170601
            (line.find("[AFLOW_PHONONS]DOS=") != std::string::npos)) //CO 170601
            && (_isGP)) {
      } else if ((line.find("[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM") != std::string::npos)) {
      } else
        outfile << line << std::endl;
    }
    //myfile.clear();   //CO
    //myfile.close();   //CO
  } else {
    throw apl::APLRuntimeError("apl::PhononCalculator::createGPAFLOWIN(); Cannot open [aflow.in] file.");
  }

  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
  outfile << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
  xvasp.str.is_vasp4_poscar_format = TRUE;
  xvasp.str.is_vasp5_poscar_format = FALSE;
  outfile << xvasp.str;
  outfile << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;
  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
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

  // Create file
  //CO - START
  string filename = directory + string("/") + _AFLOWIN_;
  aurostd::stringstream2file(outfile, filename);
  if (!aurostd::FileExist(filename)) {throw apl::APLRuntimeError("apl::PhononCalculator::createGPAFLOWIN; Cannot create ["+ _AFLOWIN_ +"] file.");}
  aurostd::ChmodFile("a+rw", filename); // CHMOD a+rw aflow.in
  //CO - END
}  //function end
// ***************************************************************************************
//create directories when EOS is ON
void eos_calculate::EOS_RUN(string& AflowIn,const double sd, const double ed, const double dm) {
  get_special_inputs(AflowIn);
  _logger << "Creating directories to calculate EOS " << apl::endl;
  vector<_xinput> EOSxRuns;  //static calculations
  vector<_xinput> PHxRuns;   //phonon calculations
  PHxRuns.clear();
  ZERO_STATIC_DIR_INDEX = -1;
  //const xstructure& pc = _supercell.getInputStructure();
  double ep = 1e-6;
  for (double i = sd; i <= ed; i += dm) {
    //if( (std::abs(i)<=ep)  ||  (abs(i-GP_VOL_DISTORTION)<=ep) || (std::abs(i+GP_VOL_DISTORTION)<=ep) )
    bool PH = false;
    if ((std::abs(i) <= ep)) PH = true;
    EOSxRuns.push_back(_xInput);
    PHxRuns.push_back(_xInput);

    uint idxRun = EOSxRuns.size() - 1;
    uint idPHVaspRun = 0;
    idPHVaspRun = PHxRuns.size() - 1;
    if (PH) ZERO_STATIC_DIR_INDEX = idPHVaspRun;
    string runname = "";
    string PHrunname = "";
    string tag = NumberToString<double>(std::abs(i));
    if (i < -0.001) {
      runname = string(_EOS_AFLOW_DIRECTORY_PREFIX_) + string("M") + string(tag);
      PHrunname = string(_PH_AFLOW_DIRECTORY_PREFIX_) + string("M") + string(tag);
    } else if (i > 0.001) {
      runname = string(_EOS_AFLOW_DIRECTORY_PREFIX_) + string("P") + string(tag);
      PHrunname = string(_PH_AFLOW_DIRECTORY_PREFIX_) + string("P") + string(tag);
    } else {
      runname = string(_EOS_AFLOW_DIRECTORY_PREFIX_) + string("0");
      if (!PH) PHrunname = string(_PH_AFLOW_DIRECTORY_PREFIX_) + string("0");
    }

    EOSxRuns[idxRun].setDirectory(PhononCalculator::_xInput.getDirectory() + runname);
    if (!PH) {
      PHxRuns[idPHVaspRun].setDirectory(PhononCalculator::_xInput.getDirectory() + PHrunname);
      PH_dir_names.push_back(PHrunname);
    }

    EOS_dir_names.push_back(runname);

    // new scale factor calulation
    double scale;
    EOSxRuns[idxRun].setXStr(PhononCalculator::_supercell.getInputStructure());
    if (!PH) PHxRuns[idPHVaspRun].setXStr(PhononCalculator::_supercell.getInputStructure());
    xmatrix<double> m(3, 3, 1, 1);
    scale = EOSxRuns[idxRun].getXStr().scale;
    m = EOSxRuns[idxRun].getXStr().lattice;
    double det = determinant(m);
    double volume = scale * scale * scale * det;
    double newscale = 0.00;
    newscale = volume + (i * volume) / 100.00;

    newscale = newscale / det;
    newscale = pow(newscale, 1. / 3.0);
    EOSxRuns[idxRun].getXStr().scale = newscale;
    PHxRuns[idPHVaspRun].getXStr().scale = newscale;
    for (uint j = 0; j < EOSxRuns[idxRun].getXStr().atoms.size(); j++) {
      EOSxRuns[idxRun].getXStr().atoms[j] = PhononCalculator::_supercell.getInputStructure().atoms[j];
      if (!PH) PHxRuns[idPHVaspRun].getXStr().atoms[j] = PhononCalculator::_supercell.getInputStructure().atoms[j];
    }
    EOS_Volumes.push_back(newscale * newscale * newscale * det);
    m.clear();
    // new scale factor calulation end

    //if directories doesn't exist create directories
    if(aurostd::FileExist(EOSxRuns[idxRun].getDirectory() + string("/") + string(_AFLOWIN_))){continue;} //CO 180406
    if(EOSxRuns[idxRun].AFLOW_MODE_VASP){
    if (aurostd::FileExist(EOSxRuns[idxRun].getDirectory() + string("/") + string(_AFLOWLOCK_)) ||
         aurostd::EFileExist(EOSxRuns[idxRun].getDirectory() + string("/OUTCAR.static"))){continue;}  //CO
    }
    _kbinFlags.KBIN_MPI_AUTOTUNE = true;
    if(EOSxRuns[idxRun].AFLOW_MODE_VASP){
    // Change format of POSCAR
    if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
        (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos)))
        EOSxRuns[idxRun].getXStr().is_vasp5_poscar_format = false;
    }
    if(EOSxRuns[idxRun].AFLOW_MODE_AIMS){
      EOSxRuns[idxRun].xaims.CONTROL.str(std::string());
      KBIN::AIMS_Produce_CONTROL(EOSxRuns[idxRun].xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);  //DEFAULT
      KBIN::AIMS_Modify_CONTROL(EOSxRuns[idxRun].xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);            //DEFAULT
    }

    writeEOSOUTPUT(EOSxRuns[idxRun]);

    if (!PH) {
      //if directories doesn't exist create directories
      if(aurostd::FileExist(PHxRuns[idPHVaspRun].getDirectory() + string("/") + string(_AFLOWIN_))){continue;}  //CO 180406
      if(aurostd::FileExist(PHxRuns[idPHVaspRun].getDirectory() + string("/") + string(_AFLOWLOCK_))){continue;}
      _kbinFlags.KBIN_MPI_AUTOTUNE = true;
      if(PHxRuns[idPHVaspRun].AFLOW_MODE_VASP){
      // Change format of POSCAR
      if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
          (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos)))
          PHxRuns[idPHVaspRun].getXStr().is_vasp5_poscar_format = false;
      }
      if(PHxRuns[idPHVaspRun].AFLOW_MODE_AIMS){
        PHxRuns[idPHVaspRun].xaims.CONTROL.str(std::string());
        KBIN::AIMS_Produce_CONTROL(PHxRuns[idPHVaspRun].xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);  //DEFAULT
        KBIN::AIMS_Modify_CONTROL(PHxRuns[idPHVaspRun].xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);            //DEFAULT
      }
      writeGPOUTPUT(PHxRuns[idPHVaspRun], false);
    }

  }  //for loop end
  EOSxRuns.clear();
  PHxRuns.clear();
}  //function end

void eos_calculate::writeEOSOUTPUT(const _xinput& xinput) {
  if(xinput.AFLOW_MODE_VASP){return createEOSAFLOWIN(xinput.xvasp);}
  if(xinput.AFLOW_MODE_AIMS){return createEOSAIMSOUTPUT(xinput.xaims);}
  else{
    throw APLRuntimeError("apl::eos_calculate::writeGPOUTPUT(); Input -> aflow.in conversion unknown.");
  }
}

void eos_calculate::createEOSAIMSOUTPUT(const _xaims& xaims_input) {
  _xaims xaims(xaims_input);
  _aimsflags aimsflags(_xFlags.aimsflags);
  if (!aurostd::FileExist(xaims.Directory)) {aurostd::DirectoryMake(xaims.Directory);}
  aurostd::DirectoryChmod("777", xaims.Directory); // CHMOD Directory 777
  //grab and write control.in from input
  //KBIN::AIMS_Produce_CONTROL(xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,aimsflags);
  //KBIN::AIMS_Modify_CONTROL(xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,aimsflags);
  KBIN::AIMS_Write_CONTROL(xaims,aimsflags);
  //KBIN::AIMS_Produce_INPUT(xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,aimsflags);
  //KBIN::AIMS_Modify_INPUT(xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,aimsflags);
  //KBIN::AIMS_Write_INPUT(xaims,aimsflags);
  //overwrite geometry.in
  xaims.GEOM.clear(); xaims.GEOM.str("");
  xaims.GEOM << xaims_input.str;

  string filename = xaims.Directory + "/" + AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM;
  aurostd::stringstream2file(xaims.GEOM, filename);
  if (!aurostd::FileExist(filename))
    throw apl::APLRuntimeError("apl::eos_calculate::createEOSAIMSOUTPUT(); Cannot create [" + AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM + "] file.");

  // CHMOD a+rw _AFLOWIN_
  aurostd::ChmodFile("a+rw", filename);
  //CO - END
}

// ***************************************************************************************
//create aflow.in for EOS calculations
void eos_calculate::createEOSAFLOWIN(const _xvasp& xvasp_input) {
  //const xstructure& pc = _supercell.getInputStructure();
  _xvasp xvasp(xvasp_input);
  _vflags vflags(_xFlags.vflags);
  //bool SPACES=FALSE;
  using aurostd::PaddedPOST;
  // Create directory if it is not created
  if (!aurostd::FileExist(xvasp.Directory)) aurostd::DirectoryMake(xvasp.Directory);
  // CHMOD Directory 777
  aurostd::DirectoryChmod("777", xvasp.Directory);
  bool SPACES = FALSE;
  //bool AFLOWIN_NEW=TRUE;
  // Create file
  //CO - START
  //string filename = xvasp.Directory + string("/")+string(_AFLOWIN_);
  //ofstream aflowin(filename.c_str(),ios_base::out);
  stringstream aflowin;
  aflowin << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  aflowin << setprecision(SetPrecision);
  // Test
  //if( !aflowin.is_open() )
  //  throw apl::APLRuntimeError("apl::PhononCalculator::createEOSAFLOWIN(); Cannot create [aflow.in] file.");
  // CHMOD a+rw aflow.in
  //aurostd::ChmodFile("a+rw",filename);
  //CO - END
  string system, formula;
  system = "";
  formula = "";
  for (uint i = 0; i < xvasp.str.species_pp.size(); i++) {
    system += xvasp.str.species_pp.at(i);
    formula += xvasp.str.species_pp.at(i) + aurostd::utype2string(xvasp.str.num_each_type.at(i));
  }

  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  aflowin << "[AFLOW]                                                                                     " << std::endl;
  aflowin << "[AFLOW]                     .o.        .o88o. oooo                                          " << std::endl;
  aflowin << "[AFLOW]                    .888.       888 `` `888                                          " << std::endl;
  aflowin << "[AFLOW]                   .8'888.     o888oo   888   .ooooo.  oooo oooo    ooo              " << std::endl;
  aflowin << "[AFLOW]                  .8' `888.     888     888  d88' `88b  `88. `88.  .8'               " << std::endl;
  aflowin << "[AFLOW]                 .88ooo8888.    888     888  888   888   `88..]88..8'                " << std::endl;
  aflowin << "[AFLOW]                .8'     `888.   888     888  888   888    `888'`888'                 " << std::endl;
  aflowin << "[AFLOW]               o88o     o8888o o888o   o888o `Y8bod8P'     `8'  `8'  .in             " << std::endl;
  aflowin << "[AFLOW]                                                                                     " << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  aflowin << "[AFLOW] * Stefano Curtarolo - (aflow V" << string(AFLOW_VERSION) << ") " << std::endl;
  aflowin << "[AFLOW] * D. Morgan, W. Setyawan, G. Hart, M. Jahnatek, S. Wang, O. Levy, K. Yang, J. Xue,  " << std::endl;
  aflowin << "[AFLOW] * R. Taylor, C. Calderon, C. Toher, R. Chepulskyy, K. Rasch, M. Buongiorno Nardelli " << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  //aflowin << "[AFLOW] Aflow automatically generated (aflow_avasp.cpp) " << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  // aflowin << AFLOWIN_SEPARATION_LINE << std::endl; // [AFLOW] **************************************************
  aflowin << "[AFLOW]SYSTEM=" << system << std::endl;
  if (xvasp.AVASP_prototype_mode == LIBRARY_MODE_PROTOTYPE || xvasp.AVASP_prototype_mode == LIBRARY_MODE_XSTRUCTURE) 
    aflowin << "[AFLOW]FORMULA=" << formula << std::endl;
  if (xvasp.str.species.size() == 1) aflowin << "#[AFLOW] single element calculation" << std::endl;
  aflowin << "[AFLOW_MODE=VASP]" << std::endl;
  aflowin << "[AFLOW_MODE_ZIP=" << _kbinFlags.KZIP_BIN << "]" << std::endl;  //CO
  aflowin << "[AFLOW_MODE_BINARY=";
  if(!_kbinFlags.KBIN_BIN.empty()){aflowin << _kbinFlags.KBIN_BIN;}
  else{aflowin << DEFAULT_VASP_BIN;}
  aflowin << "]" << std::endl;

  //CO 180130 - START
  //adding aflow.rc stuff
  aflowin << "[AFLOW_MODE_BINARY=";
  if(!_kbinFlags.KBIN_BIN.empty()){aflowin << _kbinFlags.KBIN_BIN;}
  else{aflowin << DEFAULT_VASP_BIN;}
  aflowin << "]" << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
  if(!(_kbinFlags.KBIN_MPI || XHOST.MPI)){aflowin << "#";}
  aflowin << "[AFLOW_MODE_MPI]" << std::endl;
  //be super cautious and avoid empty tags here
  string NCPUS_VAL; get_NCPUS(NCPUS_VAL); //CO 180214
  aflowin << "[AFLOW_MODE_MPI_MODE]NCPUS=" << NCPUS_VAL << " " << std::endl;
  aflowin << "[AFLOW_MODE_MPI_MODE]COMMAND =\"" << MPI_COMMAND_DEFAULT << "\" " << std::endl;
  if( _kbinFlags.KBIN_MPI_AUTOTUNE ) {aflowin << "[AFLOW_MODE_MPI_MODE]AUTOTUNE " << std::endl;}
  aflowin << "[AFLOW_MODE_MPI_MODE]BINARY=\"";
  if(!_kbinFlags.KBIN_MPI_BIN.empty()){aflowin << _kbinFlags.KBIN_MPI_BIN;}
  else{aflowin << DEFAULT_VASP_MPI_BIN;}
  aflowin << "\"" << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
  //CO 180130 - STOP

  // SYMMETRY WRITE
  if (!xvasp.aopts.flag("FLAGS::AVASP_SYMMMETRY=OFF")) {
    aflowin << "[AFLOW_SYMMETRY]CALC " << std::endl;
    aflowin << "#[AFLOW_SYMMETRY]SGROUP_WRITE " << std::endl;
    aflowin << "#[AFLOW_SYMMETRY]SGROUP_RADIUS=7.77 " << std::endl;
    aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  }
  // NEIGHBOURS WRITE
  if (!xvasp.aopts.flag("FLAGS::AVASP_NEIGHBOURS=OFF")) {
    aflowin << "#[AFLOW_NEIGHBOURS]CALC " << std::endl;
    aflowin << "[AFLOW_NEIGHBOURS]RADIUS=7.7 " << std::endl;
    aflowin << "[AFLOW_NEIGHBOURS]DRADIUS=0.1 " << std::endl;
    aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  }

  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
  if (SPACES) aflowin << std::endl;
  aflowin << "[VASP_RUN]STATIC" << std::endl;
  aflowin << "[VASP_FORCE_OPTION]WAVECAR=OFF" << std::endl;
  aflowin << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl;
  aflowin << "[VASP_FORCE_OPTION]PREC=ACCURATE" << std::endl;
  aflowin << "[VASP_FORCE_OPTION]ALGO=NORMAL" << std::endl;
  //SPIN
  if (vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) {
    if (vflags.KBIN_VASP_FORCE_OPTION_SPIN.option)
      aflowin << "[VASP_FORCE_OPTION]SPIN=ON" << std::endl;
    else
      aflowin << "[VASP_FORCE_OPTION]SPIN=OFF" << std::endl;
  }

  if (vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) aflowin << "[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << std::endl;
  if (vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry) aflowin << "[VASP_FORCE_OPTION]ABMIX=" << vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme << std::endl;
  if (vflags.KBIN_VASP_FORCE_OPTION_TYPE.isentry) aflowin << "[VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << std::endl;
  if (vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry) aflowin << "[VASP_FORCE_OPTION]AUTO_MAGMOM=" << (vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option ? "ON" : "OFF") << std::endl;
  if (vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) {
    if (vflags.KBIN_VASP_FORCE_OPTION_SPIN.option)
      aflowin << "[VASP_FORCE_OPTION]SPIN=ON" << std::endl;
    else
      aflowin << "[VASP_FORCE_OPTION]SPIN=OFF" << std::endl;
  }

  //aflowin << "[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM"<<std::endl;
  // TYPE WRITING
  if (xvasp.AVASP_flag_TYPE.xscheme.at(0) == 'M')
    aflowin << PaddedPOST("[VASP_FORCE_OPTION]TYPE=METAL", _aflowinpad_) << "// METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT) " << std::endl;
  if (xvasp.AVASP_flag_TYPE.xscheme.at(0) == 'I')
    aflowin << PaddedPOST("[VASP_FORCE_OPTION]TYPE=INSULATOR", _aflowinpad_) << "// METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT) " << std::endl;
  if (xvasp.AVASP_flag_TYPE.xscheme.at(0) == 'S')
    aflowin << PaddedPOST("[VASP_FORCE_OPTION]TYPE=SEMICONDUCTOR", _aflowinpad_) << "// METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT) " << std::endl;

  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
  if (vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry) {
    aflowin << "[VASP_FORCE_OPTION]LDAU1= ON" << std::endl;
    aflowin << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " << vflags.KBIN_VASP_LDAU_PARAMETERS << std::endl;
  }
  if (vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
    aflowin << "[VASP_FORCE_OPTION]LDAU2=ON " << std::endl;
    aflowin << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " << vflags.KBIN_VASP_LDAU_PARAMETERS << std::endl;
  }

  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  aflowin << "[VASP_INCAR_MODE_EXPLICIT]START" << std::endl;
  xvasp.AVASP_EXTRA_INCAR.clear();                          // WAHYU DEFAULT
  xvasp.AVASP_EXTRA_INCAR << "NELM = 120" << std::endl;     // WAHYU DEFAULT
  xvasp.AVASP_EXTRA_INCAR << "NELMIN=2" << std::endl;       // WAHYU DEFAULT
  xvasp.AVASP_EXTRA_INCAR << "LPLANE=.TRUE." << std::endl;  // WAHYU DEFAULT

  if (xvasp.str.atoms.size() <= 10)                           // cutoff for LREAL     // WAHYU DEFAULT
    xvasp.AVASP_EXTRA_INCAR << "LREAL=.FALSE." << std::endl;  // WAHYU DEFAULT
  else
    xvasp.AVASP_EXTRA_INCAR << "LREAL=Auto" << std::endl;    // WAHYU DEFAULT
  xvasp.AVASP_EXTRA_INCAR << "LSCALU=.FALSE." << std::endl;  // WAHYU DEFAULT
  // extra INCAR
  aflowin << xvasp.AVASP_EXTRA_INCAR.str();  // << endl;
  aflowin << "NEDOS=" << NEDOS << std::endl;
  aflowin << _PSTRESS << std::endl;  // WAHYU DEFAULT
  aflowin << "[VASP_INCAR_MODE_EXPLICIT]STOP" << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  if (SPACES) aflowin << std::endl;

  aflowin << "[VASP_POTCAR_MODE_IMPLICIT] " << std::endl;
  for (uint i = 0; i < xvasp.str.species.size(); i++) {
    if (xvasp.aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS") == FALSE) {
      if (xvasp.str.species.at(i) != "" && xvasp.str.species.at(i) != "X") {  // need because some times we have the "" around
        aflowin << "[VASP_POTCAR_FILE]" << xvasp.AVASP_potential << "/" << xvasp.str.species_pp.at(i) << std::endl;
      }
    } else {
      if (xvasp.str.species.at(i) != "" && xvasp.str.species.at(i) != "X")  // need because some times we have the "" around
        aflowin << "[VASP_POTCAR_FILE]" << KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(i)) << std::endl;
    }
  }

  //New Kpoins

  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
  aflowin << "[VASP_KPOINTS_MODE_IMPLICIT]" << std::endl;
  aflowin << EOS_KSCHEME << std::endl;
  aflowin << "[VASP_KPOINTS_FILE]KPPRA=" << EOS_KPPRA << std::endl;
  aflowin << EOS_STATIC_KSCHEME << std::endl;
  aflowin << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << EOS_STATIC_KPPRA << std::endl;
  aflowin << "[VASP_KPOINTS_FILE]BANDS_GRID=" << EOS_BANDS_GRID << std::endl;

  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************

  // new variable under VASP_POSCAR_MODE_EXPLICIT]START writing to new aflow.in
  aflowin << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
  xvasp.str.is_vasp4_poscar_format = TRUE;
  xvasp.str.is_vasp5_poscar_format = FALSE;
  aflowin << xvasp.str;
  aflowin << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;
  aflowin << AFLOWIN_SEPARATION_LINE << std::endl;

  //CO - START
  string filename = xvasp.Directory + string("/") + string(_AFLOWIN_);
  aurostd::stringstream2file(aflowin, filename);
  if (!aurostd::FileExist(filename)) {
    throw apl::APLRuntimeError("apl::PhononCalculator::createEOSAFLOWIN(); Cannot create [aflow.in] file.");
  }
  aurostd::ChmodFile("a+rw", filename);
  //aflowin.close();
  //aflowin.clear();
  //CO - END
}  //fn end
// ***************************************************************************************
template <typename T>
string eos_calculate::NumberToString(T Number) {
  ostringstream ss;
  ss << Number;
  return ss.str();
}
// ***************************************************************************************
void eos_calculate::get_special_inputs(string& AflowIn) {
  _PSTRESS = "";
  string line;
  vector<string> vlines;                           //CO
  uint line_count = 0;                             //CO
  aurostd::string2vectorstring(AflowIn,vlines);
  //aurostd::efile2vectorstring(_AFLOWIN_, vlines);  //CO
                                                   //  ifstream myfile (_AFLOWIN_.c_str()); //CO
  //if (!myfile.is_open()) //CO
  if (!vlines.size())  //CO
  { throw apl::APLRuntimeError("apl::PhononCalculator::get_special_inputs(); Cannot read ["+_AFLOWIN_+"] file."); }
  //while ( getline (myfile,line) ) //CO
  while (line_count < vlines.size())  //CO
  {
    line = vlines[line_count++];  //CO
    if (line == "") continue;
    if (line[0] == '#') continue;
    if ((line.find("PSTRESS") != std::string::npos)) {
      _PSTRESS = line;
      break;
    }
  }
  //myfile.clear(); //CO
  //myfile.close(); //CO
}
// ***************************************************************************************
}  //apl end
