// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2018             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu
/*
 This class creates aflow.in files for QHA, SCQHA, and QHA3P methods.
*/

#include "aflow_apl.h"

namespace apl {
  // ***************************************************************************************
  QHA_AFLOWIN_CREATOR::QHA_AFLOWIN_CREATOR(Supercell& sc, vector<ClusterSet>& clst,
                                           _xinput& xinput, //_xvasp& xvasp,
					   _aflags& aflags, _kflags& kflags,
					   _xflags& xflags, //_vflags& vflags,
					   string& AflowIn,
					   Logger& l) : PhononCalculator(sc, clst, xinput, aflags, kflags, xflags, AflowIn, l) { //xvasp, aflags, kflags, vflags, l) {
    clear();
    _log.open(_logfile.c_str());
    if (!_log.is_open())
      throw apl::APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::create_aflowin_phonon() Cannot create [aflow_qha_aflowcreate_information.out] file.");
    _log << std::setprecision(8) << std::fixed;
  }
  // ***************************************************************************************
  QHA_AFLOWIN_CREATOR::~QHA_AFLOWIN_CREATOR() { this->clear(); }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::clear()
  {
    _log.close();
    _logfile="aflow.qha.distortions.out";

    _is_gp_on=false;
    _is_gp_A_on=false;
    _is_gp_B_on=false;
    _is_gp_C_on=false;


    _is_sc_gp_on=false;
    _is_sc_gp_A_on=false;  
    _is_sc_gp_B_on=false;
    _is_sc_gp_C_on=false;

    _is_eos=false;
    _is_eos_A=false;
    _is_eos_B=false; 
    _is_eos_C=false;
    _is_edos_accurate_on=false;

    _scqha_dir_names.clear();
    _gp_dir_names.clear();
    _ph_dir_names.clear();
    _eos_dir_names.clear();
    _gp_volumes.clear();
    _scqha_volumes.clear();
    _eos_volumes.clear();
    _ph_volumes.clear();
    _zero_index=-10;;
    _pstress="";
    _EOS_VOL_START=0.0;
    _EOS_VOL_END=0.0;
    _EOS_VOL_INC=0.0;
    //_EOS_KPPRA="";
    //_EOS_BANDS_GRID="";
    _EOS_KSCHEME="G";
    _NEDOS=0;
    _EOS_STATIC_KSCHEME="G";
    _EOS_STATIC_KPPRA=0;
    _scqha_vol_distortion=0.0;
    _gp_vol_distortion=0.0;
    _EOS_VOL_START=0.0;
    _EOS_VOL_END=0.0;
    _EOS_VOL_INC=0.0;
    _lattice_index=-1;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::run_qha()
  {
    _logger<<"Arranging all qha-options "<<apl::endl;

    //options are mutually exclusive
    if((_is_gp_on) && (_is_eos)){
      _is_eos=true;
      _is_eos_A=false;
      _is_eos_B=false;
      _is_eos_C=false;
    }else if((_is_gp_A_on) && (_is_eos)){
      _is_eos=false;
      _is_eos_A=true;
      _is_eos_B=false;
      _is_eos_C=false;
    }else if((_is_gp_B_on) && (_is_eos)){
      _is_eos=false;
      _is_eos_A=false;
      _is_eos_B=true;
      _is_eos_C=false;
    }else if((_is_gp_C_on) && (_is_eos)){
      _is_eos=false;
      _is_eos_A=false;
      _is_eos_B=false;
      _is_eos_C=true;
    }


    if((_is_sc_gp_on) && (_is_eos)){
      _is_eos=true;
      _is_eos_A=false;
      _is_eos_B=false; 
      _is_eos_C=false;
    }else if((_is_sc_gp_A_on) && (_is_eos)){
      _is_eos=false;
      _is_eos_A=true;
      _is_eos_B=false; 
      _is_eos_C=false;
    }else if((_is_sc_gp_B_on) && (_is_eos)){
      _is_eos=false;
      _is_eos_A=false;
      _is_eos_B=true; 
      _is_eos_C=false;
    }else if((_is_sc_gp_C_on) && (_is_eos)){
      _is_eos=false;
      _is_eos_A=false;
      _is_eos_B=false; 
      _is_eos_C=true;
    }
 
    //error check
    if(_scqha_vol_distortion==0) 
      {
	if(_is_sc_gp_on){
	  throw APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::run_qha(): scqha_distortion==0 but _is_sc_gp is on.");
	}else if(_is_sc_gp_A_on){
	  throw APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::run_qha(): scqha_distortion==0 but _is_sc_gp_A is on.");
	}else if(_is_sc_gp_B_on){
	  throw APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::run_qha(): scqha_distortion==0 but _is_sc_gp_B is on.");
	}else if(_is_sc_gp_C_on){
	  throw APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::run_qha(): scqha_distortion==0 but _is_sc_gp_C is on.");
	}
      }
   
    //set index for anisotropic distortion
    if(_is_gp_A_on || _is_sc_gp_A_on || _is_eos_A)_lattice_index=1;
    else if(_is_gp_B_on || _is_sc_gp_B_on || _is_eos_B)_lattice_index=2;
    else if(_is_gp_C_on || _is_sc_gp_C_on || _is_eos_C)_lattice_index=3;

    //aflowin creating options
    //[phonon_option] 0->gp   || 1->sc-gp   || 2-> eos-phonon   || 3->eos-static 
    //[phonon_option] 4->gp_X || 5->sc-gp_X || 6-> eos-phonon_X || 7->eos-static-X

    if(_is_gp_on){
      create_aflowin_phonon(_gp_vol_distortion, 0);
    }else if(_is_gp_A_on || _is_gp_B_on || _is_gp_C_on){
      create_aflowin_phonon_X(_gp_vol_distortion, 4);
    }


    if(_is_sc_gp_on){
      create_aflowin_phonon(_scqha_vol_distortion, 1);
    }else if(_is_sc_gp_A_on || _is_sc_gp_B_on || _is_sc_gp_C_on){
      create_aflowin_phonon_X(_scqha_vol_distortion, 5);
    }

    if(_is_eos){
      if(_is_gp_on || _is_gp_A_on || _is_gp_B_on || _is_gp_C_on){
	create_aflowin_phonon(0, 2);
      }
      create_aflowin_phonon(0, 3);
      _log<<"Equilibrium directory index "<< _zero_index << '\n';
    }else if(_is_eos_A || _is_eos_B || _is_eos_C){
      if(_is_gp_on || _is_gp_A_on || _is_gp_B_on || _is_gp_C_on){
	create_aflowin_phonon_X(0, 6);
      }
      create_aflowin_phonon_X(0, 7);
      _log<<"Equilibrium directory index "<< _zero_index << '\n';
    }
  }
  // ***************************************************************************************
  //[phonon_option] 0->gp   || 1->sc-gp   || 2-> eos-phonon   || 3->eos-static 
  //[phonon_option] 4->gp_X || 5->sc-gp_X || 6-> eos-phonon_X || 7->eos-static-X
  void QHA_AFLOWIN_CREATOR::create_aflowin_phonon(const double distortion, const int phonon_option)
  {
    if(phonon_option==0){
      _logger<<"Creating distorted configurations to calculate QHA properties "<<apl::endl;
    }else if(phonon_option==1){
      _logger<<"Creating distorted configurations to calculate SCQHA/QHA3P properties "<<apl::endl;
    }else if(phonon_option==2){
      _logger<<"Creating distorted configurations to calculate QHA EOS "<<apl::endl;
    }else if(phonon_option==3){
      _logger<<"Creating distorted configurations to calculate QHA static energies "<<apl::endl;
    }

    if(phonon_option==0){
      _log<<"Creating QHA Gruneisen aflow.in"<<'\n';
    }else if(phonon_option==1){
      _log<<"Creating SCQHA Gruneisen aflow.in"<<'\n';
    }else if(phonon_option==2){
      _log<<"Creating EOS_phonon aflow.in"<<'\n';
    }else if(phonon_option==3){
      _log<<"Creating EOS_static aflow.in"<<'\n';
    }

    double Start=0.0, End=0.0, Inc=0.0;
    if(phonon_option==0 || phonon_option==1){
      Start=-distortion;
      End=distortion;
      Inc=distortion;
    }else if(phonon_option==2 || phonon_option==3){
      Start=_EOS_VOL_START;
      End=_EOS_VOL_END;
      Inc=_EOS_VOL_INC;
    }
    _log<<"#"<<setw(25)<<"dir_name"<<setw(15)<<"scale"<<setw(15)<<"volume"<<'\n';

    //smallest distortion value
    double smallest_distortion=100.0;
    bool smallest_distortion_found=false;
    if(phonon_option==3){
      vector<double> distortions;
      for(double i=Start; i<=End; i+=Inc){
	distortions.push_back(std::abs(i));
      }
      smallest_distortion=*std::min_element(distortions.begin(), distortions.end());
    }

    int idxRun=0;
    string APL_DIR=PhononCalculator::_xInput.getDirectory();
    vector<_xinput> vaspRuns; 
    for(double i=Start; i<=End; i+=Inc){

      vaspRuns.push_back(_xInput);
      idxRun = vaspRuns.size()-1;
      double scale=0.00;
      vaspRuns[idxRun].setXStr(PhononCalculator::_supercell.getPrimitiveStructure());
      xmatrix<double> m(3,3,1,1);
      scale=vaspRuns[idxRun].getXStr().scale;
      m= vaspRuns[idxRun].getXStr().lattice;
      double det=determinant(m);
      double volume=scale*scale*scale*det;

      if(_iszero(i)){
	if(phonon_option==0){
	  _gp_volumes.push_back(volume);
	}else if(phonon_option==1){
	  _scqha_volumes.push_back(volume);
	}
      }

      //include equilibrium configuration
      if(phonon_option==3){
	if((std::abs(i)==smallest_distortion) && (!smallest_distortion_found)){
	  create_aflowin_static_zero();
	  smallest_distortion_found=true;
	  _zero_index=idxRun;
	}
      }
 
      if(_iszero(i))continue;

      string runname="";
      if(phonon_option==0 || phonon_option==1){
	runname=get_phonon_runname(i, distortion);
      }else if(phonon_option==2){
	runname=get_phonon_runname(i);
      }else if(phonon_option==3){
	runname=get_static_runname(i);
      }

      vaspRuns[idxRun].setDirectory(APL_DIR + "./" + runname);

      // new scale factor calulation
      double newscale = 0.0;
      newscale=volume + (i*volume)/100.00;
      newscale= newscale/det;
      newscale=pow(newscale, 1./3.0);
      vaspRuns[idxRun].getXStr().scale = newscale;
      double newvolume=newscale*newscale*newscale*det;

      if(phonon_option==0){
	_gp_dir_names.push_back(runname);
	_gp_volumes.push_back(newvolume);
      }else if(phonon_option==1){
	_scqha_volumes.push_back(newvolume);
	_scqha_dir_names.push_back(runname);
      }else if(phonon_option==2){
	_ph_volumes.push_back(newvolume);
	_ph_dir_names.push_back(runname);
      }else if(phonon_option==3){
	_eos_dir_names.push_back(runname);
	_eos_volumes.push_back(newvolume);
      }

      // new scale factor calulation end
      _log<<setw(25)<<vaspRuns[idxRun].getDirectory()<<setw(15)<<newscale<<setw(15)<<newvolume<<'\n';

      vaspRuns[idxRun].getXStr().atoms = PhononCalculator::_supercell.getPrimitiveStructure().atoms;
      m.clear();
      if( aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/")+string(_AFLOWLOCK_) ) ||
	  aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/OUTCAR.static") ) ) continue;

      if(phonon_option==2){
	if (( (_is_sc_gp_on) && std::abs(i)==_scqha_vol_distortion) || ((_is_sc_gp_A_on) && std::abs(i)==_scqha_vol_distortion) ||
	    ((_is_sc_gp_B_on) && std::abs(i)==_scqha_vol_distortion) || ((_is_sc_gp_C_on) && std::abs(i)==_scqha_vol_distortion)) continue;
      }
      _logger<<"Creating "<< runname <<apl::endl;
      if(phonon_option<3){ 
	write_phonon_OUTPUT(vaspRuns[idxRun], phonon_option);
      }else if(phonon_option==3){ 
	write_static_OUTPUT(vaspRuns[idxRun]);
      }
          
    }//for loop end
    vaspRuns.clear();
  }
  // ***************************************************************************************
  //[phonon_option] 0->gp   || 1->sc-gp   || 2-> eos-phonon   || 3->eos-static 
  //[phonon_option] 4->gp_X || 5->sc-gp_X || 6-> eos-phonon_X || 7->eos-static-X
  void QHA_AFLOWIN_CREATOR::create_aflowin_phonon_X(const double distortion, const int phonon_option)
  {

    if(phonon_option==4){
      if(_is_gp_A_on){ 
	_logger<<"Creating distorted configurations to calculate Gruneisen-A Parameter "<<apl::endl;
      }else if(_is_gp_B_on){ 
	_logger<<"Creating distorted configurations to calculate Gruneisen-B Parameter "<<apl::endl;
      }else if(_is_gp_C_on){ 
	_logger<<"Creating distorted configurations to calculate Gruneisen-C Parameter "<<apl::endl;
      }
    }else if(phonon_option==5){
      if(_is_sc_gp_A_on){ 
	_logger<<"Creating distorted configurations to calculate SC-Gruneisen-A Parameter "<<apl::endl;
      }else if(_is_sc_gp_B_on){ 
	_logger<<"Creating distorted configurations to calculate SC-Gruneisen-B Parameter "<<apl::endl;
      }else if(_is_sc_gp_C_on){ 
	_logger<<"Creating distorted configurations to calculate SC-Gruneisen-C Parameter "<<apl::endl;
      }
    }else if(phonon_option==6){
      if(_is_eos_A){ 
	_logger<<"Creating distorted configurations to calculate EOS-phonon-A Parameter "<<apl::endl;
      }else if(_is_eos_B){ 
	_logger<<"Creating distorted configurations to calculate EOS-phonon-B Parameter "<<apl::endl;
      }else if(_is_eos_C){ 
	_logger<<"Creating distorted configurations to calculate EOS-phonon-C Parameter "<<apl::endl;
      }
    }else if(phonon_option==7){
      if(_is_eos_A){ 
	_logger<<"Creating distorted configurations to calculate EOS-static-A Parameter "<<apl::endl;
      }else if(_is_eos_B){ 
	_logger<<"Creating distorted configurations to calculate EOS-static-B Parameter "<<apl::endl;
      }else if(_is_eos_C){ 
	_logger<<"Creating distorted configurations to calculate EOS-static-C Parameter "<<apl::endl;
      }
    }

    if(_is_sc_gp_A_on || _is_gp_A_on){ 
      _log<<"#"<<setw(25)<<"dir_name"<<setw(15)<<"|a|"<<'\n';
    }else if(_is_sc_gp_B_on || _is_gp_B_on){ 
      _log<<"#"<<setw(25)<<"dir_name"<<setw(15)<<"|b|"<<'\n';
    } else if(_is_sc_gp_C_on || _is_gp_C_on){ 
      _log<<"#"<<setw(25)<<"dir_name"<<setw(15)<<"|c|"<<'\n';
    }

    double Start=0.0, End=0.0, Inc=0.0;
    if(phonon_option==4 || phonon_option==5){
      Start=-distortion;
      End=distortion;
      Inc=distortion;
    }else if(phonon_option==6 || phonon_option==7){
      Start=_EOS_VOL_START;
      End=_EOS_VOL_END;
      Inc=_EOS_VOL_INC;
    }

    //smallest distortion value
    double smallest_distortion=100.0;
    bool smallest_distortion_found=false;
    if(phonon_option==7){
      vector<double> distortions;
      for(double i=Start; i<=End; i+=Inc){
	distortions.push_back(std::abs(i));
      }
      smallest_distortion=*std::min_element(distortions.begin(), distortions.end());
    }

    int idxRun=0;
    string APL_DIR=PhononCalculator::_xInput.getDirectory();
    vector<_xinput> vaspRuns; vaspRuns.clear();

    for(double i=Start; i<=End; i+=Inc)
      { 
	vaspRuns.push_back(_xInput);
	idxRun = vaspRuns.size()-1;
	vaspRuns[idxRun].setXStr(PhononCalculator::_supercell.getPrimitiveStructure());
	xmatrix<double> m=vaspRuns[idxRun].getXStr().lattice;
	xvector<double> lattice_X(3,1);
	for(int j=1; j<=3; j++)lattice_X[j]=m[_lattice_index][j];
	vaspRuns[idxRun].getXStr().atoms = PhononCalculator::_supercell.getPrimitiveStructure().atoms;
	double mod_X=aurostd::modulus(lattice_X);


	if(_iszero(i)){
	  if(phonon_option==4){
	    _gp_volumes.push_back(mod_X);
	  }else if(phonon_option==5){
	    _scqha_volumes.push_back(mod_X);
	  }
	}

	//include equilibrium configuration
	if(phonon_option==7){
	  if((std::abs(i)==smallest_distortion) && (!smallest_distortion_found)){
	    create_aflowin_static_zero_X();
	    smallest_distortion_found=true;
	    _zero_index=idxRun;
	  }
	}

	if(_iszero(i))continue;

	string runname="";
	if(phonon_option==4 || phonon_option==5){
	  runname=get_phonon_runname(i, distortion);
	}else if(phonon_option==6){
	  runname=get_phonon_runname(i);
	}else if(phonon_option==7){
	  runname=get_static_runname(i);
	}

	vaspRuns[idxRun].setDirectory(APL_DIR + "./" + runname);

	//new lattice vectors
	xvector<double> lattice_X_new=lattice_X;
	for(int j=1; j<=3; j++){
	  lattice_X_new[j]= lattice_X[j]+(i*lattice_X[j])/100.00;
	}
	mod_X=aurostd::modulus(lattice_X_new);

	for(int j=1; j<=3; j++){
	  vaspRuns[idxRun].getXStr().lattice[_lattice_index][j]=lattice_X_new[j];
	}

	if(phonon_option==4){
	  _gp_dir_names.push_back(runname);
	  _gp_volumes.push_back(mod_X);
	}else if(phonon_option==5){
	  _scqha_volumes.push_back(mod_X);
	  _scqha_dir_names.push_back(runname);
	}else if(phonon_option==6){
	  _ph_volumes.push_back(mod_X);
	  _ph_dir_names.push_back(runname);
	}else if(phonon_option==7){
	  _eos_dir_names.push_back(runname);
	  _eos_volumes.push_back(mod_X);
	}

	_log<<setw(25)<<vaspRuns[idxRun].getDirectory()<<setw(15)<<mod_X<<'\n';
	vaspRuns[idxRun].getXStr().atoms = PhononCalculator::_supercell.getPrimitiveStructure().atoms;
	m.clear();

	if( aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/")+string(_AFLOWLOCK_) ) ||
	    aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/OUTCAR.static") ) ) continue;
	if(phonon_option==6){
	  if (( (_is_sc_gp_on) && std::abs(i)==_scqha_vol_distortion) || ((_is_sc_gp_A_on) && std::abs(i)==_scqha_vol_distortion) ||
	      ((_is_sc_gp_B_on) && std::abs(i)==_scqha_vol_distortion) || ((_is_sc_gp_C_on) && std::abs(i)==_scqha_vol_distortion)) continue;
	}
	_logger<<"Creating "<< runname <<apl::endl;
	if(phonon_option<7){
	  write_phonon_OUTPUT(vaspRuns[idxRun], phonon_option);
	}else if(phonon_option==7){
	  write_static_OUTPUT(vaspRuns[idxRun]);
	}
      }//for loop end
    vaspRuns.clear();
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::write_phonon_OUTPUT(const _xinput& xinput, const int phonon_option) {
    if(xinput.AFLOW_MODE_VASP){
      return write_aflowin_phonon(xinput.xvasp, phonon_option);
    }
    //if(xinput.AFLOW_MODE_AIMS){return create_aflowin_phonon_phonon(xinput.xaims);}
    else{
      throw APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::write_phonon_OUTPUT(); Input -> aflow.in conversion unknown.");
    }
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::write_static_OUTPUT(const _xinput& xinput) {
    if(xinput.AFLOW_MODE_VASP){return write_static_AFLOWIN(xinput.xvasp);}
    //if(xinput.AFLOW_MODE_AIMS){return create_aflowin_phonon_phonon(xinput.xaims);}
    else{
      throw APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::write_static_OUTPUT(); Input -> aflow.in conversion unknown.");
    }
  }
  // ***************************************************************************************
  //create aflow.in for GRUNEISEN_A phonon calculation
  //0->gp || 1->sc-gp || 2->gp_A || 3->scgp_A
  void QHA_AFLOWIN_CREATOR::write_aflowin_phonon(const _xvasp& xvasp_input, const int phonon_option)
  {
    _xvasp xvasp(xvasp_input);
    // Create directory if it is not created
    if (!aurostd::FileExist(xvasp.Directory))
      aurostd::DirectoryMake(xvasp.Directory);
    aurostd::DirectoryChmod("777", xvasp.Directory);

    stringstream outfile;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(Set_QHA_Precision);

    vector<string> vlines; vlines.clear(); 
    uint line_count = 0;    
    string line;
    uint ENTRY = 0;
    aurostd::efile2vectorstring(_AFLOWIN_, vlines);
    bool include_qha_variables=false; 
    bool eos_included=false;
    bool gp_included=false;
    bool scqha_included=false;

    if (vlines.size())  
      {
	while (line_count < vlines.size())  
	  {
	    line = vlines[line_count++]; 
	    if(line=="")continue;
	    if(line[0]=='#')continue;
	    if(line.find("AFLOW_AAPL") != std::string::npos)continue;
	    if(line.find("[VASP_POSCAR_MODE_EXPLICIT]START") != std::string::npos){ENTRY=true;}
	    if(ENTRY){if(line.find("[VASP_POSCAR_MODE_EXPLICIT]STOP") != std::string::npos){ENTRY=false;}}
	    //
	    else if(line.find("[VASP_KPOINTS_FILE]KSCHEME=") != std::string::npos){_EOS_KSCHEME=line;outfile<<line<<std::endl;}
	    else if(line.find("[VASP_KPOINTS_FILE]STATIC_KSCHEME=") != std::string::npos){_EOS_STATIC_KSCHEME=line;outfile<<line<<std::endl;}
	    //exclude SPRIM from all the subdirectory runs
	    else if( (line.find("[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM") != std::string::npos) ){}
	    //phonon entry
	    else if((line.find("[AFLOW_PHONONS]CALC") != std::string::npos) || (line.find("[AFLOW_APL]CALC") != std::string::npos)) {
	      outfile<<line<<std::endl;
            }else if(line.find("NEDOS=") != std::string::npos){
            }else if(line.find("ENGINE=") != std::string::npos){
	      outfile<<line<<std::endl;
	    }else if(line.find("SUPERCELL=") != std::string::npos){
	      outfile<<line<<std::endl;
	    }else if(line.find("MINATOMS=") != std::string::npos){
	      outfile<<line<<std::endl;
	    }else if(line.find("DISMAG=") != std::string::npos){
	      outfile<<line<<std::endl;
	    }else if(line.find("ZEROSTATE=") != std::string::npos){
	      outfile<<line<<std::endl;
	    }else if(line.find("POLAR=") != std::string::npos){
	      outfile<<line<<std::endl;
	    }else if(line.find("DC=") != std::string::npos){
	      if(phonon_option==0 || phonon_option==1 || phonon_option==4 || phonon_option==5){
		outfile<<line<<std::endl;
	      }
	    }else if(line.find("DCINITSG=") != std::string::npos){
	      outfile<<line<<std::endl;
	    }else if(line.find("DCUSERPATH=") != std::string::npos){
	      if(phonon_option==0 || phonon_option==1 || phonon_option==4 || phonon_option==5){
		outfile<<line<<std::endl;
	      }
	    }else if(line.find("DOS=") != std::string::npos){
	      if(phonon_option!=0 &&  phonon_option!=4){
		outfile<<line<<std::endl;
	      }
	    }else if(line.find("GP_DISTORTION=") != std::string::npos){
	      if(phonon_option==0 || phonon_option==4){
                outfile << setprecision(2);
		outfile<<"[AFLOW_APL]GP_DISTORTION="<<_gp_vol_distortion<<std::endl;
                outfile << setprecision(Set_QHA_Precision);
                gp_included=true;
	      }
	    }else if(line.find("SCQHA_DISTORTION=") != std::string::npos){
	      if(phonon_option==1 || phonon_option==5){
                outfile << setprecision(2);
		outfile<<"[AFLOW_APL]SCQHA_DISTORTION="<<_scqha_vol_distortion<<std::endl;
                outfile << setprecision(Set_QHA_Precision);
                scqha_included=true;
	      }
	    }else if(line.find("EOS_DISTORTION_RANGE=") != std::string::npos){
	      if(phonon_option==2 || phonon_option==6){
                outfile << setprecision(2);
		outfile<<"[AFLOW_APL]EOS_DISTORTION_RANGE="<<_EOS_VOL_START<<":"<<_EOS_VOL_END<<":"<<_EOS_VOL_INC<<std::endl;
                outfile << setprecision(Set_QHA_Precision);
                eos_included=true;
	      }
	    }else if((line.find("[AFLOW_PHONONS]") != std::string::npos) || 
		     (line.find("[AFLOW_APL]") != std::string::npos) || 
		     (line.find("[AFLOW_QHA]") != std::string::npos)){
	      if(!include_qha_variables)
		{
		  if(phonon_option==0){
		    outfile<<"[AFLOW_APL]GRUNEISEN_SD=y"<<std::endl;
		  }else if(phonon_option==1){
		    outfile<<"[AFLOW_APL]SCQHA_SD=y"<<std::endl;
		  }else if(phonon_option==2){
		    outfile<<"[AFLOW_APL]EOS_SD=y"<<std::endl;
		  }else if(phonon_option==4){
		    if(_is_gp_A_on){
		      outfile<<"[AFLOW_APL]GRUNEISEN_A_SD=y"<<std::endl;
		    }else if(_is_gp_B_on){
		      outfile<<"[AFLOW_APL]GRUNEISEN_B_SD=y"<<std::endl;
		    }else if(_is_gp_C_on){
		      outfile<<"[AFLOW_APL]GRUNEISEN_C_SD=y"<<std::endl;
		    }
		  }else if(phonon_option==5){
		    if(_is_sc_gp_A_on){
		      outfile<<"[AFLOW_APL]SCQHA_A_SD=y"<<std::endl;
		    }else if(_is_sc_gp_B_on){
		      outfile<<"[AFLOW_APL]SCQHA_B_SD=y"<<std::endl;
		    }else if(_is_sc_gp_C_on){
		      outfile<<"[AFLOW_APL]SCQHA_C_SD=y"<<std::endl;
		    }
		  }else if(phonon_option==6){
		    if(_is_sc_gp_A_on){
		      outfile<<"[AFLOW_APL]EOS_SD=y"<<std::endl;
		    }else if(_is_sc_gp_B_on){
		      outfile<<"[AFLOW_APL]EOS_SD=y"<<std::endl;
		    }else if(_is_sc_gp_C_on){
		      outfile<<"[AFLOW_APL]EOS_SD=y"<<std::endl;
		    }
		  }
		}
	      include_qha_variables=true;
	    }else{ 
	      outfile<<line<<std::endl;
	    }
	  }
      } else {
      throw apl::APLRuntimeError("apl::PhononCalculator::createGPAFLOWIN(); Cannot open [aflow.in] file.");
    }

             if(!eos_included)
            {
	      if(phonon_option==2 || phonon_option==6){
                outfile << setprecision(2);
		outfile<<"[AFLOW_APL]EOS_DISTORTION_RANGE="<<_EOS_VOL_START<<":"<<_EOS_VOL_END<<":"<<_EOS_VOL_INC<<std::endl;
                outfile << setprecision(Set_QHA_Precision);
	      }
             }
              if(!gp_included)
             {
	      if(phonon_option==0 || phonon_option==4){
                outfile << setprecision(2);
		outfile<<"[AFLOW_APL]GP_DISTORTION="<<_gp_vol_distortion<<std::endl;
                outfile << setprecision(Set_QHA_Precision);
	      }
             }
              if(!scqha_included){
	      if(phonon_option==1 || phonon_option==5){
                outfile << setprecision(2);
		outfile<<"[AFLOW_APL]SCQHA_DISTORTION="<<_scqha_vol_distortion<<std::endl;
                outfile << setprecision(Set_QHA_Precision);
	      }
    }

    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
    xvasp.str.is_vasp4_poscar_format = TRUE;
    xvasp.str.is_vasp5_poscar_format = FALSE;
    outfile << xvasp.str;
    outfile << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;

    string filename = xvasp.Directory + string("/") + string(_AFLOWIN_);
    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      throw apl::APLRuntimeError("apl::QHA_AFLOWIN_CREATOR::create_aflowin_phonon; Cannot create [aflow.in] file.");
    }
    aurostd::ChmodFile("a+rw", filename);
  }//function end
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::correcting_scqha_vol_distortion()
  {
    double min_value=std::min(std::abs(_EOS_VOL_START), _EOS_VOL_END);
    bool found=false;
    vector<double> dis; dis.clear();
    for(double i=-min_value; i<=min_value; i+=_EOS_VOL_INC)
      {
	if(_isequal(std::abs(i), _scqha_vol_distortion))
	  {
	    found=true;
	  }
	double diff=std::abs((std::abs(i)-_scqha_vol_distortion));
	if(!_iszero(i))
	  {
	    dis.push_back(diff);
	  }
      }

    if(!found)
      {
	min_value=dis[0];
	for(uint i=0; i!=dis.size(); i++)
	  {
	    if(dis[i]<min_value)
	      {
		min_value=dis[i];
	      }
	  }
	_scqha_vol_distortion+=min_value;
      }
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setGP(bool b1, bool b2, bool b3, bool b4)
  {
    _is_gp_on=b1;
    _is_gp_A_on=b2;
    _is_gp_B_on=b3;
    _is_gp_C_on=b4;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setSCGP(bool b1, bool b2, bool b3, bool b4)
  {
    _is_sc_gp_on=b1;
    _is_sc_gp_A_on=b2;
    _is_sc_gp_B_on=b3;
    _is_sc_gp_C_on=b4;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setEOS(bool b)
  {
    _is_eos=b;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setGP_VOL_DISTORTION(double d)
  {
    _gp_vol_distortion=d;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setSCGP_VOL_DISTORTION(double d)
  {
    _scqha_vol_distortion=d;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setEOS_distortion_range(double a, double b, double c)
  {
    _EOS_VOL_START =a;
    _EOS_VOL_END = b; 
    _EOS_VOL_INC = c;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setEOS_STATIC_KPPRA(int s)
  {
    _EOS_STATIC_KPPRA=s;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setEOS_NEDOS(int s)
  {
    _NEDOS=s;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::setEOS_STATIC_KSCHEME(string s)
  {
    _EOS_STATIC_KSCHEME=s;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::create_aflowin_static_zero()
  {
    int idxRun=0;
    string APL_DIR=PhononCalculator::_xInput.getDirectory();
    vector<_xinput> vaspRuns; 
    for(double i=0; i<=0; i++)
      {
	vaspRuns.push_back(_xInput);
	idxRun = vaspRuns.size()-1;
	double scale=0.00;
	vaspRuns[idxRun].setXStr(PhononCalculator::_supercell.getPrimitiveStructure());
	xmatrix<double> m(3,3,1,1);
	scale=vaspRuns[idxRun].getXStr().scale;
	m= vaspRuns[idxRun].getXStr().lattice;
	double det=determinant(m);
	double volume=scale*scale*scale*det;

	string runname=get_static_runname(i);
	vaspRuns[idxRun].setDirectory(APL_DIR + "./" + runname);
	_eos_dir_names.push_back(runname);
	_eos_volumes.push_back(volume);

	_log<<setw(25)<<vaspRuns[idxRun].getDirectory()<<setw(15)<<scale<<setw(15)<<volume<<'\n';
	// new scale factor calulation end

	vaspRuns[idxRun].getXStr().atoms = PhononCalculator::_supercell.getPrimitiveStructure().atoms;
	m.clear();
	if( aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/")+string(_AFLOWLOCK_) ) ||
	    aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/OUTCAR.static") ) ) continue;
	_logger<<"Creating "<< runname <<apl::endl;
	write_static_OUTPUT(vaspRuns[idxRun]);
      }
    vaspRuns.clear();
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::create_aflowin_static_zero_X()
  {
    int idxRun=0;
    string APL_DIR=PhononCalculator::_xInput.getDirectory();
    vector<_xinput> vaspRuns; vaspRuns.clear();

    for(double i=0; i<=0; i++)
      { 
	vaspRuns.push_back(_xInput);
	idxRun = vaspRuns.size()-1;
	vaspRuns[idxRun].setXStr(PhononCalculator::_supercell.getPrimitiveStructure());
	xmatrix<double> m=vaspRuns[idxRun].getXStr().lattice;
	xvector<double> lattice_X(3,1);
	for(int j=1; j<=3; j++)lattice_X[j]=m[_lattice_index][j];
	vaspRuns[idxRun].getXStr().atoms = PhononCalculator::_supercell.getPrimitiveStructure().atoms;
	double mod_X=aurostd::modulus(lattice_X);

	string runname=get_static_runname(i);
	vaspRuns[idxRun].setDirectory(APL_DIR + "./" + runname);

	_eos_dir_names.push_back(runname);
	_eos_volumes.push_back(mod_X);

	_log<<setw(25)<<vaspRuns[idxRun].getDirectory()<<setw(15)<<mod_X<<'\n';
	vaspRuns[idxRun].getXStr().atoms = PhononCalculator::_supercell.getPrimitiveStructure().atoms;
	m.clear();

	if( aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/")+string(_AFLOWLOCK_) ) ||
	    aurostd::FileExist( vaspRuns[idxRun].getDirectory() + string("/OUTCAR.static") ) ) continue;

	_logger<<"Creating "<< runname <<apl::endl;
	write_static_OUTPUT(vaspRuns[idxRun]);
      }
    vaspRuns.clear();
  }
  // ***************************************************************************************
  //create aflow.in for EOS calculations
  void  QHA_AFLOWIN_CREATOR::write_static_AFLOWIN(const _xvasp& xvasp_input) 
  {
    _xvasp xvasp(xvasp_input);
    _vflags vflags(_xFlags.vflags);
    using aurostd::PaddedPOST;
    if (!aurostd::FileExist(xvasp.Directory)) aurostd::DirectoryMake(xvasp.Directory);
    aurostd::DirectoryChmod("777", xvasp.Directory);
    bool SPACES = FALSE;
    stringstream aflowin;
    aflowin << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    aflowin << setprecision(Set_QHA_Precision);
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
    aflowin << "[AFLOW_MODE=VASP]" << std::endl;
    //aflowin << "[AFLOW_MODE_ZIP=bzip]" << std::endl; //CO
    aflowin << "[AFLOW_MODE_ZIP=" << _kbinFlags.KZIP_BIN << "]" << std::endl;  //CO
    if (xvasp.str.species.size() == 1)
      aflowin << "#[AFLOW] single element calculation" << std::endl;
    aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
    if (_kbinFlags.KBIN_MPI) {
      aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
      // michal outfile << "#[AFLOW_MODE_BINARY=" << _kbinFlags.KBIN_BIN << "]" << std::endl;
      aflowin << "[AFLOW_MODE_BINARY=" << _kbinFlags.KBIN_BIN << "]" << std::endl;
      aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
      aflowin << "[AFLOW_MODE_MPI]" << std::endl;
      if (_kbinFlags.KBIN_MPI_AUTOTUNE) {
	aflowin << "[AFLOW_MODE_MPI_MODE]AUTOTUNE" << std::endl;
      } else {
	aflowin << "[AFLOW_MODE_MPI_MODE]NCPUS=MAX" << std::endl;
      }
      aflowin << "[AFLOW_MODE_MPI_MODE]BINARY=" << _kbinFlags.KBIN_MPI_BIN << std::endl;
    } else {
      aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
      aflowin << "[AFLOW_MODE_BINARY=" << _kbinFlags.KBIN_BIN << "]" << std::endl;
      aflowin << AFLOWIN_SEPARATION_LINE << std::endl;
      // michal outfile << "#[AFLOW_MODE_MPI_MODE]BINARY=\"mpi" << _kbinFlags.KBIN_BIN << "\"" << std::endl;
      aflowin << "[AFLOW_MODE_MPI_MODE]BINARY=\"mpi" << _kbinFlags.KBIN_BIN << "\"" << std::endl;
      aflowin << "[AFLOW_MODE_MPI_MODE]NCPUS=MAX" << std::endl;
      aflowin << "[AFLOW_MODE_MPI_MODE]COMMAND=\"mpirun -np\" " << std::endl;
      aflowin << "[AFLOW_MODE_MPI_MODE]AUTOTUNE" << std::endl;
    }
    aflowin << AFLOWIN_SEPARATION_LINE << std::endl;  // [AFLOW] **************************************************
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
    aflowin << "NEDOS=" << _NEDOS << std::endl;
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
    aflowin << _EOS_KSCHEME << std::endl;
    aflowin << "[VASP_KPOINTS_FILE]KPPRA=" << _EOS_STATIC_KPPRA  << std::endl;
    aflowin << _EOS_STATIC_KSCHEME << std::endl;
    aflowin << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << _EOS_STATIC_KPPRA << std::endl;
    //aflowin << "[VASP_KPOINTS_FILE]BANDS_GRID=" << _EOS_BANDS_GRID << std::endl;

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
  }//fn end
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::get_pstress()
  {
    _PSTRESS = "";
    string line;
    vector<string> vlines;                           //CO
    uint line_count = 0;                             //CO
    aurostd::string2vectorstring(_AFLOWIN_,vlines);
    if (!vlines.size())  //CO
      { throw apl::APLRuntimeError("apl::PhononCalculator::get_special_inputs(); Cannot read ["+_AFLOWIN_+"] file."); }
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
  }
  // ***************************************************************************************
  vector<string> QHA_AFLOWIN_CREATOR::get_scqha_dir_names()
  {
    return _scqha_dir_names;
  }
  // ***************************************************************************************
  vector<string> QHA_AFLOWIN_CREATOR::get_gp_dir_names()
  {
    return _gp_dir_names;
  }
  // ***************************************************************************************
  vector<string> QHA_AFLOWIN_CREATOR::get_ph_dir_names()
  {
    return _ph_dir_names;
  }
  // ***************************************************************************************
  vector<string> QHA_AFLOWIN_CREATOR::get_eos_dir_names()
  {
    return _eos_dir_names;
  }
  // ***************************************************************************************
  vector<double> QHA_AFLOWIN_CREATOR::get_gp_volumes()
  {
    return _gp_volumes;
  }
  // ***************************************************************************************
  vector<double> QHA_AFLOWIN_CREATOR::get_scqha_volumes()
  {
    return _scqha_volumes;
  }
  // ***************************************************************************************
  vector<double> QHA_AFLOWIN_CREATOR::get_eos_volumes()
  {
    return _eos_volumes;
  }
  // ***************************************************************************************
  int QHA_AFLOWIN_CREATOR::get_zero_index()
  {
    return _zero_index;
  }
  // ***************************************************************************************
  double QHA_AFLOWIN_CREATOR::get_scqha_vol_distortion()
  {
    return _scqha_vol_distortion;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::set_edos_accurate(bool b)
  {
    _is_edos_accurate_on=b;
  }
  // ***************************************************************************************
  void QHA_AFLOWIN_CREATOR::close_log()
  {
    _log.close();
  }
  // ***************************************************************************************
  template <typename T>
  string QHA_AFLOWIN_CREATOR::NumToStr ( T Number )
  {
    ostringstream ss;
    ss << Number;
    return ss.str();
  }
  // ***************************************************************************************
  string QHA_AFLOWIN_CREATOR::get_phonon_runname(const double i, const double distortion)
  {
    string runname="";
    string tag=NumToStr<double>(std::abs(i));
    if( _iszero(std::abs(i-distortion)) )
      {
	runname=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("P")+string(tag);
      }else{
      runname=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("M")+string(tag);
    }
    return runname;
  }
  // ***************************************************************************************
  string QHA_AFLOWIN_CREATOR::get_phonon_runname(const double i)
  {
    string runname="";
    string tag=NumToStr<double>(std::abs(i));
    if(i<0.0){
      runname=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("M")+string(tag);
    }else if(i>0.0){
      runname=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("P")+string(tag);
    }
    else{
      runname=string(_AFLOW_QHA_PHONS_DIRECTORY_PREFIX_)+string("0");
    }
    return runname;
  }
  // ***************************************************************************************
  string QHA_AFLOWIN_CREATOR::get_static_runname(const double i)
  {
    string runname="";
    string tag=NumToStr<double>(std::abs(i));
    if(i<0.0){
      runname=string(_EOS_AFLOW_DIRECTORY_PREFIX_)+string("M")+string(tag);
    }else if(i>0.0){
      runname=string(_EOS_AFLOW_DIRECTORY_PREFIX_)+string("P")+string(tag);
    }
    else{
      runname=string(_EOS_AFLOW_DIRECTORY_PREFIX_)+string("0");
    }
    return runname;
  }
  // ***************************************************************************************
}//apl end
