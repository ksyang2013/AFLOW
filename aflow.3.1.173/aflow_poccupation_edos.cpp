// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// functions written by KESONG YANG
// 2010-2011: kesong.yang@gmail.com

#ifndef _AFLOW_POCCUPATION_EDOS_CPP_
#define _AFLOW_POCCUPATION_EDOS_CPP_

#include "aflow_contrib_kesong.h"
// const double KBOLTZEV = 8.617332478E-5; //(ev.K^(-1))

// ***************************************************************************
// void pocc::POCC_SetOptions(string options, string& directory, double& T, double& DOS_Emin, double& DOS_Emax, double& DOSSCALE)
// ***************************************************************************
namespace pocc {
  void POCC_SetOptions(string options, string& directory, double& T, double& DOS_Emin, double& DOS_Emax, double& DOSSCALE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: BEGIN" << endl;
    }
    directory="./"; DOS_Emin=DEFAULT_DOS_EMIN/2.0; DOS_Emax=DEFAULT_DOS_EMAX/2.0; DOSSCALE=DEFAULT_DOS_SCALE; // some defaults  
    T = 300.0;  //300K room temperature default
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: options=[" << options << "]" << endl;
    }
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: tokens.size()=" << tokens.size() << endl;
    }
    if(LDEBUG) {
        for(uint i=0;i<tokens.size();i++) cerr << "pocc::POCC_SetOptions: tokens.at(i)=" << tokens.at(i) << endl;
    }
    
    //Usage: program 0; --pocc_dos 1; directory 2; temperature 3, DOS_Emin 4; DOS_Emax 5; DOSSCALE 6
    
    if(tokens.size()>=1) {
        directory = tokens.at(0); 
    }
    if(tokens.size()>=2) {
        T = aurostd::string2utype<double>(tokens.at(1)); 
    }
    if(tokens.size()>=3) {
        DOS_Emin = aurostd::string2utype<double>(tokens.at(2)); 
    }
    if(tokens.size()>=4) {
        DOS_Emax = aurostd::string2utype<double>(tokens.at(3));
    }
    if(tokens.size()>=5) {
        DOSSCALE = aurostd::string2utype<double>(tokens.at(4));
    }
    
    // [OBSOLETE] if(directory=="") directory="./";
    // [OBSOLETE] if(T<0.01) T=0.01;
    // [OBSOLETE] if(DOS_Emin<=-999.9) DOS_Emin=DEFAULT_DOS_EMIN;
    // [OBSOLETE] if(DOS_Emax<=-999.9) DOS_Emax=DEFAULT_DOS_EMAX;
    // [OBSOLETE] if(DOSSCALE<=-999.9) DOSSCALE=DEFAULT_DOS_SCALE;
    
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: directory=[" << directory << "]" << endl;
    }
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: T=" << T << endl;
    }
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: DOS_Emin=" << DOS_Emin << endl;
    }
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: DOS_Emax=" << DOS_Emax << endl;
    }
    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: DOSSCALE=" << DOSSCALE << endl;
    }

    if(LDEBUG) {
        cerr << "pocc::POCC_SetOptions: END" << endl;
    }
  }
} // namespace pocc

// ***************************************************************************
// pocc::POCC_DOS
// ***************************************************************************
namespace pocc {
  void POCC_DOS(ostream &oss,string options) {
    string directory;
    double T, DOS_Emin, DOS_Emax, DOSSCALE; 
    pocc::POCC_SetOptions(options,directory,T,DOS_Emin,DOS_Emax,DOSSCALE);
    oss << "pocc::POCC_DOS is working on " << directory  << endl;
    oss << "Temperature is " << T << "K" << endl;
    pocc::POCC_GENERATE_OUTPUT(directory,T,DOS_Emin,DOS_Emax,DOSSCALE);
  }
} // namespace pocc

// ***************************************************************************
// pocc::POCC_Already_Calculated_Input(const string& str_AflowIn)
// ***************************************************************************
namespace pocc {
  bool POCC_Already_Calculated_Input(const string& str_AflowIn) {
    return aurostd::substring2bool(str_AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]START.POCC");
  }
} // namespace pocc

// ***************************************************************************
// void pocc::POCC_CalculateDistribution(vector<double>& vdelta_toten, vector<double>& vprob, const double& T, const string& NameDist, const vector<int>& vDE)
// ***************************************************************************
namespace pocc {
  void POCC_CalculateDistribution(vector<double>& vdelta_toten, vector<double>& vprob, const double& T, const string& NameDist, const vector<int>& vDE) {
    //Check 
    if(vdelta_toten.size()!=vDE.size()) {
      cerr << "Error!!! Check the number of POSCARs, Degeneracy!" << endl;
    }
    vector<double> vFi; vprob.clear();
    double Fi, DEi;
    for (uint i=0; i<vdelta_toten.size();i++) {
      double delta_toten = vdelta_toten.at(i);
      DEi = vDE.at(i);
      if(NameDist=="B")  {
          Fi = DEi*exp((-1.0)*delta_toten/(KBOLTZEV*T));    // Boltzmann distribution
      }
      if(NameDist=="FD") {
          Fi = DEi*1.0/(exp(delta_toten/(KBOLTZEV*T))+1.0); // Fermi-Dirac distribution
      }
      if(NameDist=="BE") {
          Fi = DEi*1.0/(exp(delta_toten/(KBOLTZEV*T))-1.0); // Bose-Einstein distribution
      }
      vFi.push_back(Fi);
    }
    double sum_vFi = aurostd::sum<double> (vFi);
    for (uint i=0; i<vFi.size();i++) {
      double probi = vFi.at(i)/sum_vFi;
      vprob.push_back(probi);
    }
  }
} // namespace pocc

// ***************************************************************************
// vector<vector<double> > pocc::SpinSplitDOS(const vector<vector<double> >& vva)
// ***************************************************************************
namespace pocc {
  vector<vector<double> > SpinFlipDOS(const vector<vector<double> >& vva) {
    vector<vector<double> > vvb; vvb.clear();
    if(vva.at(0).size()%2==0) {cerr << "It is not spin-polarized! Aborting! " << endl; exit(1);}
    for (uint i=0; i<vva.size();i++) {
      vector<double> vtmp; vtmp.clear();
      vtmp.push_back(vva.at(i).at(0));
      for (uint j=1; j<vva.at(0).size();j=j+2) {
	vtmp.push_back(vva.at(i).at(j+1)*(-1));   //exchange columns 1 and 2 spin_dn becomes spin_up
	vtmp.push_back(vva.at(i).at(j)*(-1));   //spin_up becomes spin_dn
      }
      vvb.push_back(vtmp);
    }
    return vvb;
  }
} // namespace pocc

// ***************************************************************************
// vector<vector<double> > pocc::SpinSplitDOS(const vector<vector<double> >& vva)
// ***************************************************************************
namespace pocc {
  vector<vector<double> > SpinSplitDOS(const vector<vector<double> >& vva) {
    vector<vector<double> > vvb;
    //if(vva.at(0).size()%2==1) {cerr << "It is already spin-polarized! Aborting! " << endl; exit(1);}
    for (uint i=0; i<vva.size();i++) {
      vector<double> vtmp; vtmp.clear();
      vtmp.push_back(vva.at(i).at(0));
      for (uint j=1; j<vva.at(i).size();j++) {
	double tmp = 0.5*vva.at(i).at(j);  
	vtmp.push_back(tmp); vtmp.push_back((-1)*tmp);
      }
      vvb.push_back(vtmp);
    }
    return vvb;
  }
} // namespace pocc


// ***************************************************************************
// vector<vector<double> > pocc::POCC_Formalise(SPIN_FLAG, vweight, TOTALPDOS)
// ***************************************************************************
namespace pocc {
  vector<vector<double> > POCC_Formalise(const bool& SPIN_FLAG, const double& weight, const double& mag, const vector<vector<double> >& TDOS) {
    vector<vector<double> > vva = TDOS;
    vector<vector<double> > vvb;
    //weight
    for (uint i=0; i<vva.size();i++) {
      for (uint j=1; j<vva.at(i).size();j++) {
	vva.at(i).at(j) *= weight;
      }
    }
    //spin
    if(SPIN_FLAG) {
      if(TDOS.at(0).size()==3 || TDOS.at(0).size()==4) { //TDOS or TOTALPDOS, E, s, p, d
	vvb = pocc::SpinSplitDOS(vva);
      } else if(TDOS.at(0).size()==5 && abs(mag) < 1E-6) { //must be TOTALPDOS: E, s, p, d, f, no-spin-polarized
	vvb = pocc::SpinSplitDOS(vva);
      } else if((TDOS.at(0).size()==5) && ((-1)*mag > 1E-3 )) { //if magnetic moment is smaller than -0.001 and size is 5,  must be TDOS, 
	vvb = SpinFlipDOS(vva);
      } else if((TDOS.at(0).size()==7) && ((-1)*mag > 1E-3 )) { //if magnetic moment is smaller than -0.001, size is 7, must be PDOS
	vvb = SpinFlipDOS(vva);
      } else if((TDOS.at(0).size()==9) && ((-1)*mag > 1E-3 )) { //if magnetic moment is smaller than -0.001, size is 7, must be PDOS
	vvb = SpinFlipDOS(vva);
      } else {
	vvb = vva;
      }
    } else {
      vvb = vva;
    }
    return vvb;
  }
} // namespace pocc

// ***************************************************************************
//void ExtracAllPOSCARSFromAflowin(vector<xstructure>& vxstr, const string& str_aflowin)
// ***************************************************************************
void ExtracAllPOSCARSFromAflowin(vector<xstructure>& vxstr, const string& str_aflowin) {
  vxstr.clear();
  vector<string> vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING;
  aurostd::substring2strings(str_aflowin,vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING,"[VASP_POSCAR_MODE_EXPLICIT]START.");
  // load up the structures
  for(uint i=0;i<vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
    string START="[VASP_POSCAR_MODE_EXPLICIT]START";
    string STOP="[VASP_POSCAR_MODE_EXPLICIT]STOP";
    START="[VASP_POSCAR_MODE_EXPLICIT]START."+vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i);
    STOP="[VASP_POSCAR_MODE_EXPLICIT]STOP."+vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i);
    stringstream POSCAR;POSCAR.clear();POSCAR.str(std::string());
    if(aurostd::substring2bool(str_aflowin,START) && aurostd::substring2bool(str_aflowin,STOP)) {
      aurostd::ExtractToStringstreamEXPLICIT(str_aflowin,POSCAR,START,STOP);
    }
    vxstr.push_back(xstructure(POSCAR,IOVASP_AUTO));
  }
}

// ***************************************************************************
// void GetDegeneracyFromVPOSCAR(const vector<xstructure>& vxstr, vector<int>& vDE)
// ***************************************************************************
void GetDegeneracyFromVPOSCAR(const vector<xstructure>& vxstr, vector<int>& vDE) {
  for (uint i=0; i<vxstr.size();i++) {
    string title = vxstr.at(i).title;
    if(!aurostd::substring2bool(title, "DG=")) {
      cerr << "Error!!! There are no degeneracy data! Please regenerate your " << _AFLOWIN_ << " file!" << endl;
      exit(1);
    }
    //CO 180220 - multiply occupied sites will yield titles with more than one DG= value
    vector<string> vtitle,vtitle2;
    aurostd::string2tokensAdd(title,vtitle," ");
    int DEI=1;
    for(uint j=0;j<vtitle.size();j++){
      if(!aurostd::substring2bool(vtitle[j],"DG=")){continue;}
      aurostd::string2tokens(vtitle[j],vtitle2,"=");
      if(vtitle2.size()!=2){
        cerr << "Error! Bad degeneracy value read! Please regenerate your " << _AFLOWIN_ << " file!" << endl;
        exit(1);
      }
      DEI*=aurostd::string2utype<int>(vtitle2[1]);
    }
    //[OBSOLETE CO 180220]string last_part_title = vtitle.at(vtitle.size()-1);
    //[OBSOLETE CO 180220]vector<string> vtitle2;
    //[OBSOLETE CO 180220]aurostd::string2tokensAdd(last_part_title,vtitle2,"=");
    //[OBSOLETE CO 180220]int DEI;
    //[OBSOLETE CO 180220]DEI=aurostd::string2utype<int>(vtitle2.at(vtitle2.size()-1));
    vDE.push_back(DEI);
  }
}

// ***************************************************************************
// void pocc::POCC_GENERATE_DOSDATA(const string& directory, const double& T, vector<vector<double> >& DOS, vector<double>& POCC_Efermi, double& POCC_vmag, vector<double>& vprob)
// ***************************************************************************
namespace pocc {
  void POCC_GENERATE_DOSDATA(const string& directory, const double& T, vector<vector<double> >& TDOS_ONLY, vector<vector<double> >& PDOS_ONLY, vector<double>& POCC_Efermi, double& POCC_mag, double& Egap_net, vector<double>& Egap, vector<double>& vprob) {
    //Warnning: DOS is absolute value, no shift, and the output DOS, it is format is Energy, s, p, d, f, TDOS, TDOS_sum
    
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) {
        cerr << "pocc::POCC_GENERATE_DOSDATA: BEGIN" << endl;
    }
    if(LDEBUG) {
        cerr << "pocc::POCC_GENERATE_DOSDATA: directory=[" << directory << "]" << endl;
    }
 
    string aflowin,MESSAGE="pocc::POCC_DOS ERROR";
    aflowin=string(directory +"/"+_AFLOWIN_);
    if(!aurostd::FileExist(aflowin)) {cerr << MESSAGE << ": file not found " << aflowin << endl; exit(1);}
    string str_AflowIn; aurostd::file2string(aflowin, str_AflowIn);
    if(pocc::POCC_Already_Calculated_Input(str_AflowIn)) {
      //Fix Degeneracy, Reading all POSCARs from aflowin
      vector<xstructure> vxstr;
      ExtracAllPOSCARSFromAflowin(vxstr,str_AflowIn);
      vector<int> vDE; vDE.clear();
      GetDegeneracyFromVPOSCAR(vxstr, vDE);

      //Fixing Degeneracy can be very easy if we begin from below, but I do not want to repeat calculations, so I fixt it above
      vector<string> vrun, vdoscar_files, voutcar_files;
      string command, aflow_pocc_out;
      command = "grep \"\\[VASP_POSCAR_MODE_EXPLICIT\\]START\\.\" " + aflowin;
      aflow_pocc_out = aurostd::execute2string(command);
      aurostd::StringSubst(aflow_pocc_out,"[VASP_POSCAR_MODE_EXPLICIT]START.",""); //replace string
      aurostd::string2vectorstring(aflow_pocc_out,vrun);

      //Check files
      bool FLAG_ALLFILES_EXIST = true;
      ostringstream aus;
      ofstream FileMESSAGE;
      //double kpt_tol;
      for(uint i=0;i<vrun.size();i++) {
	string outcar_file = aurostd::CleanFileName(directory + "/ARUN."+vrun.at(i)+"/OUTCAR.static.bz2");
	string doscar_file = aurostd::CleanFileName(directory + "/ARUN."+vrun.at(i)+"/DOSCAR.static.bz2");
	if(aurostd::FileExist(outcar_file)) {
	  if(LDEBUG) {
          aus << "00000  MESSAGE POCC OUTCAR file OK: " << outcar_file << " " << endl;
      }
	  if(LDEBUG) {
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
	} else {
	  aus << "ERROR" << ": file not found " << outcar_file << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  FLAG_ALLFILES_EXIST = false;
	}
	if(aurostd::FileExist(doscar_file)) {
	  if(LDEBUG) {
          aus << "00000  MESSAGE POCC DOSCAR file OK: " << doscar_file << " " << endl;
      }
	  if(LDEBUG) {
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
	} else {
	  aus << "ERROR" << ": file not found " << doscar_file << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
	  FLAG_ALLFILES_EXIST = false;
	}
	voutcar_files.push_back(outcar_file);  // worse
	vdoscar_files.push_back(doscar_file);
      }
      if(!FLAG_ALLFILES_EXIST) {
          exit(1);
      }

      //ions, total energy, spin flag, magnetic moment
      vector<int> vions; vector<double> vtoten_per_atom; vector<double> vmag; vector<double> vEgap_net;
      vector<vector<double> > vEgap;
      bool POCC_SPIN_FLAG = false;  //default non-spin-polarized
      for (uint i=0; i<vrun.size();i++) {
	if(LDEBUG) {
        cerr << "pocc::POCC_GENERATE_DOSDATA: vrun.at(" << i << ")=" << vrun.at(i) << endl;
    }
  //xstructure xstr=vxstr[i]; //CO 171002 - grab xstr's too
  //CO 171002 - using tolerance from symmetry calc - START
  //if(xstr.CalculateSymmetry()){kpt_tol=xstr.sym_eps;}
  //else{kpt_tol=SYM::defaultTolerance(xstr);}
  //CO 171002 - using tolerance from symmetry calc - STOP

	xOUTCAR outcar_aus;
	if(!outcar_aus.GetPropertiesFile(voutcar_files.at(i))){
	  aus << "ERROR" << ": OUTCAR.static reading error " << outcar_aus.ERROR << endl;
	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    exit(1);
  }
  double EFERMI=outcar_aus.Efermi;
	outcar_aus.GetBandGap();
	vions.push_back(outcar_aus.NIONS);
	vtoten_per_atom.push_back(outcar_aus.enthalpy_atom);
	vmag.push_back(outcar_aus.mag_cell);
	string outcar_bands=voutcar_files.at(i);aurostd::StringSubst(outcar_bands,"static","bands");
	if(aurostd::FileExist(outcar_bands)) {
	  if(LDEBUG) {
          cerr << "pocc::POCC_GENERATE_DOSDATA: loading outcar=" << outcar_bands << endl;
      }
	  outcar_aus.GetPropertiesFile(outcar_bands);
	  outcar_aus.GetBandGap(EFERMI);
	} 
	vEgap_net.push_back(outcar_aus.Egap_net);
	//	cerr << "CAMILO ??? outcar_aus.Egap.size()=" << outcar_aus.Egap.size() << endl;
	vEgap.push_back(outcar_aus.Egap);

    if(LDEBUG||1) {
        cerr << "pocc::POCC_GENERATE_DOSDATA: vrun.at(" << i << ")=" << vrun.at(i) << "; NIONS=" << vions.at(i) << "; enthalpy_atom=" << vtoten_per_atom.at(i) << "; mag_cell=" << vmag.at(i) << "; Egap_net=" << vEgap_net.at(i) << endl;    // corey - should read out enthalpy/atom from STATIC, not BANDS
    }

//	if(LDEBUG||1) {
//	cerr << "pocc::POCC_GENERATE_DOSDATA: vrun.at(" << i << ")=" << vrun.at(i) << "; NIONS=" << outcar_aus.NIONS << "; enthalpy_atom=" << outcar_aus.enthalpy_atom << "; mag_cell=" << outcar_aus.mag_cell << "; Egap_net=" << outcar_aus.Egap_net << endl;   //corey
//	}
      }

      int min_IONS = min(vions); 
      vector<double> vweight; //different ions in each structure
      for (uint i=0; i<vrun.size();i++) {   // will use this for DOS as it is extensive
	vweight.push_back(double(1.0*min_IONS/vions.at(i))); //1.0 make int to double, otherwise it is zero!
      }
  
      double min_toten = min(vtoten_per_atom);
      vector<double> vdelta_toten;
      for (uint i=0; i<vtoten_per_atom.size();i++) {
	vdelta_toten.push_back(vtoten_per_atom.at(i) - min_toten);
      }
      
      // Calculate the probability boltzmann distribution
      vprob.clear(); 
      pocc::POCC_CalculateDistribution(vdelta_toten, vprob, T, "B", vDE); //Boltzmann

      if(0) { 
	cerr << "check vprob: ";
	for (uint i=0; i<vprob.size();i++) {
        cerr << vprob.at(i) << " ";
    }
	cerr << " = " << aurostd::sum(vprob) << endl;// exit(0);
      }

      // vector<double> vmag_abs; 
      //  for (uint i=0; i<vmag.size();i++) {
      //  vmag_abs.push_back(abs(vmag.at(i))); //ferromagnetic alignment
      //  }

      POCC_mag=0.0;
      for (uint i=0;i<vmag.size();i++) {
	POCC_mag+=abs(vmag.at(i))*vprob.at(i); //Calculate average absolute magnetic moment
      }
      
      Egap_net=0.0;
      for (uint i=0;i<vEgap_net.size();i++) {
	Egap_net+=vEgap_net.at(i)*vprob.at(i); //Calculate average Egap_net
      }

      Egap.clear();
      for(uint j=0;j<vEgap.at(0).size();j++) {
	Egap.push_back(0);
	for (uint i=0;i<vEgap.size();i++) {
      //Calculate average Egap
	  Egap.at(j)+=vEgap.at(0).at(j)*vprob.at(i);    //corey - otherwise, it crashes if you have Egap up AND Egap down (magnetized system)
      //Egap.at(j)+=vEgap.at(i).at(j)*vprob.at(i);  //corey
      }
      }

      // Get TDOS & TOTALPDOS DATA
      vector<vector<vector<double> > > POCC_TDOS, POCC_TOTALPDOS;
      for (uint i=0; i<vdoscar_files.size();i++) {
	string doscar_file = vdoscar_files.at(i); stringstream ss_doscar; aurostd::efile2stringstream(doscar_file, ss_doscar); 
	string outcar_file = voutcar_files.at(i); stringstream ss_outcar; aurostd::efile2stringstream(outcar_file, ss_outcar);
	double Efermi; vector<vector<double> > TDOS, TOTALPDOS;
  //CO 180218 - let's not mess around with kesong's functions too much
  //for now, assume POCC runs have standard DOSCAR.static (PDOS in it)
  //if no PDOS, exit
  //I will fix later
  if(!(estructure::GET_DOS_DATA(ss_doscar, ss_outcar, Efermi, TDOS, TOTALPDOS) && TOTALPDOS.size()>0)){
    cerr << "ERROR: DOSCAR extraction failed, perhaps there is no PDOS, needed for POCC" << endl;
    exit(1);
  }
  //format TDOS, if spin and non-spin coexist, then format them into spin
	vector<vector<double> > TDOSf, TOTALPDOSf;
	TDOSf = pocc::POCC_Formalise(POCC_SPIN_FLAG, vweight.at(i), vmag.at(i), TDOS); 
	TOTALPDOSf = pocc::POCC_Formalise(POCC_SPIN_FLAG, vweight.at(i), vmag.at(i), TOTALPDOS);
	POCC_Efermi.push_back(Efermi);
	POCC_TDOS.push_back(aurostd::ShiftFirstColumn(TDOSf, -1*Efermi)); //conver it into absolute value
	POCC_TOTALPDOS.push_back(aurostd::ShiftFirstColumn(TOTALPDOSf, -1*Efermi));
      }

      vector<vector<double> > POCC_TDOS_normalized = aurostd::NormalizeAndSum3DVector(POCC_TDOS, vprob);  //normalize TDOS
      vector<vector<double> > POCC_TOTALPDOS_normalized = aurostd::NormalizeAndSum3DVector(POCC_TOTALPDOS, vprob); //normalize pdos
      TDOS_ONLY = POCC_TDOS_normalized;
      PDOS_ONLY = POCC_TOTALPDOS_normalized;
    }
  }
} // namespace pocc

// ***************************************************************************
// void pocc::POCC_COMBINE_TDOS_PDOS_ONEDOS(const vector<vector<double> >& TDOS, const vector<vector<double> >& PDOS, vector<vector<double> >& DOS)
// ***************************************************************************
namespace pocc {
  void POCC_COMBINE_TDOS_PDOS_ONEDOS(const vector<vector<double> >& TDOS, const vector<vector<double> >& PDOS, vector<vector<double> >& DOS) {
    if(TDOS.size()!=PDOS.size()) {cerr << " TDOS and PDOS have different size! Aborting! " << endl; exit(1);}
    vector<double> vtmp;
    if(TDOS.at(0).size()==3) { //non-spin
      for (uint i=0; i<TDOS.size();i++) {
	vtmp.clear();
	vtmp.push_back(TDOS.at(i).at(0)); //Energy
	for (uint j=1; j<PDOS.at(i).size();j++) {
        vtmp.push_back(PDOS.at(i).at(j)); //s, p, d
    }
	vtmp.push_back(TDOS.at(i).at(1)); //TDOS
	DOS.push_back(vtmp);
      }
    }
    else if(TDOS.at(0).size()==5) { //spin
      for (uint i=0; i<TDOS.size();i++) {
	vtmp.clear();
	vtmp.push_back(TDOS.at(i).at(0)); //Energy
	for (uint j=1; j<PDOS.at(i).size();j++) {
        vtmp.push_back(PDOS.at(i).at(j)); //s, p ,d, f
    }
	vtmp.push_back(TDOS.at(i).at(1)); //TDOS up
	vtmp.push_back(TDOS.at(i).at(2)); //TDOS dn
	DOS.push_back(vtmp);
      }
    }
    else {
      cerr << "ERROR!!!!!!!!" << endl;
    }
  }
} // namespace pocc

// ***************************************************************************
// string pocc::POCC_GENERATE_GNUPLOTSCRIPT(vector<vector<double> >& DOS, const string& SystemName, const string& dosdatafile, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE, const double& DOSMAX)
// ***************************************************************************
namespace pocc {
  string POCC_GENERATE_GNUPLOTSCRIPT(vector<vector<double> >& DOS, const string& SystemName, const string& dosdatafile, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE, const double& DOSMAX) {
    stringstream  ss_gnu; ss_gnu.str(std::string());
    ss_gnu << "#Generated by AFLOW (Kesong Yang [kesong.yang@gmail.com], 2011, Duke)" << endl;
    ss_gnu << "set term postscript eps enhanced color font \"Times-Roman, 40\" size 18, 10.125" << endl;
    ss_gnu << "set output " << "\"" << SystemName <<"_DOS.eps" << "\"" << endl;
    ss_gnu << "set title \"" << SystemName << "\""<< endl;
    ss_gnu << endl;

    ss_gnu << "#DOS PLOT" << endl;
    ss_gnu << "set xtics 2" << endl;
    ss_gnu << "set ytics" << endl;
    ss_gnu << "set xrange [" << DOS_Emin << ":" << DOS_Emax << "]" << endl;
    if(DOS.at(0).size()==10||(DOS.at(0).size()==12)) {
      if(DOSMAX*DOSSCALE!=0) {
	ss_gnu << "set yrange [" << DOSMAX*DOSSCALE*(-1) <<":"  << DOSMAX*DOSSCALE << "]" << endl;
      } else {
	ss_gnu << "set yrange [0:2]" << endl;
      }
    }
    if(DOS.at(0).size()==5||(DOS.at(0).size()==6)) {
      if(DOSMAX*DOSSCALE!=0) {
	ss_gnu << "set yrange [0:" << DOSMAX*DOSSCALE << "]" << endl;
      }
    }
    ss_gnu << endl;
    ss_gnu << "set label '" << AFLOWLIB_CONSORTIUM_STRING << "' at screen 0.70, 0.02 font \"Times-Roman, 32\"" << endl;
    ss_gnu << "set xlabel 'Energy (eV)' offset graph 0.00" << endl;
    ss_gnu << "set ylabel 'Density of States (States/eV)' offset graph 0.00" << endl;
    ss_gnu << "set label 'E_F' at 0.05, graph 0.95" << endl;
    ss_gnu << "set label 'T = " << aurostd::utype2string(T) << "K' " << "at 0.3, graph 0.95" << endl;
    ss_gnu << "set arrow from 0, 0 to first 0, graph 1 nohead lt 3 lw 1.5" << endl;
    ss_gnu << "set arrow from 0, 0 to first 0, graph 0 nohead lt 3 lw 1.5" << endl;
    ss_gnu << endl;
    ss_gnu << "plot[][] \\" << endl;
    if(DOS.at(0).size()==5) {
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt -1 lw 2 title 'Total'" << endl;
    }
    if(DOS.at(0).size()==6) {
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt  5 lw 6 title 'f', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:6 w l lt -1 lw 2 title 'Total'" << endl;
    }
    if(DOS.at(0).size()==9) { // Only works for s, p, d and f orbitals; Spin-polarized
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  2 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt  3 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:6 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:7 w l lt  8 lw 6 notitle, \\"  << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:8 w l lt -1 lw 2 title 'Total', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:9 w l lt -1 lw 2 notitle" << endl;
    }
    if(DOS.at(0).size()==11) { // Only works for s, p, d and f orbitals; Spin-polarized
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  2 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt  3 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:6 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:7 w l lt  8 lw 6 notitle, \\"  << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:8 w l lt  5 lw 6 title 'f', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:9 w l lt  5 lw 6 notitle, \\"  << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:10 w l lt -1 lw 2 title 'Total', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:11 w l lt -1 lw 2 notitle"<< endl;
    }
    ss_gnu << endl;
    return ss_gnu.str();
  }
} // namespace pocc


// ***************************************************************************
// void pocc::POCC_GENERATE_OUTPUT(const string& directory, const double& T, const double& DOS_Emax, double& DOS_Emin, const double& DOSSCALE)
// ***************************************************************************
namespace pocc {
  void POCC_GENERATE_OUTPUT(const string& directory, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE) {
    //produce DOS data
    if(!XHOST.is_command("gnuplot")) {cerr << "AFLOW V" << string(AFLOW_VERSION) << " - pocc::POCC_GENERATE_OUTPU. ERROR gnuplot is necessary." << endl;exit(1);}; 
    if(!XHOST.is_command("convert")) {cerr << "AFLOW V" << string(AFLOW_VERSION) << " - pocc::POCC_GENERATE_OUTPU. ERROR convert is necessary." << endl;exit(1);}; 
    vector<vector<double> > TDOS_ONLY, PDOS_ONLY, DOS;  vector<double> vEfermi,Egap,vprob; double mag,Egap_net;
    POCC_GENERATE_DOSDATA(directory,T,TDOS_ONLY,PDOS_ONLY,vEfermi,mag,Egap_net,Egap,vprob);
    POCC_COMBINE_TDOS_PDOS_ONEDOS(TDOS_ONLY,PDOS_ONLY,DOS);
    DOS = aurostd::ShiftFirstColumn(DOS, max(vEfermi)); //shift to Efermi
    double DOSMAX = aurostd::FindMaxIn2DvectorExcept1stColumn(DOS, DOS_Emin, DOS_Emax);
    
    //write 2D vector into files
    string str_dos = aurostd::vector2string(DOS);
    string dosdatafile = "DOSDATA_" + aurostd::utype2string(T) + "K";
    aurostd::string2file(str_dos, dosdatafile);

    //generate gnuplot script
    string POCC_NAME = "untitled_" + aurostd::utype2string(T) + "K";
    string gnuplotscript = "GNUPLOT_POCC_" + POCC_NAME + "_DOS.gp";
    string str_gnu = POCC_GENERATE_GNUPLOTSCRIPT(DOS, POCC_NAME, dosdatafile, T, DOS_Emin, DOS_Emax, DOSSCALE, DOSMAX);
    aurostd::string2file(str_gnu, gnuplotscript);

    //run system command to plot
    aurostd::execute(XHOST.command("gnuplot")+" " + gnuplotscript);
    aurostd::execute(XHOST.command("convert")+" -background white ./" + POCC_NAME + "_DOS.eps ./" + POCC_NAME + "_DOS.png");

    if(!XHOST.vflag_control.flag("KEEP::GPL")) {
      aurostd::execute("rm  -f " + gnuplotscript);
      //     aurostd::execute("rm  -f " + dosdatafile); // always keep it
    }
    if(!XHOST.vflag_control.flag("KEEP::EPS")) {
      aurostd::execute("rm  -f " + POCC_NAME + "_DOS.eps");
    } 
  }
} // namespace pocc

// ***************************************************************************
// pocc::POCC_BANDGAP(string options) {
// ***************************************************************************
namespace pocc {
  void POCC_BANDGAP(string options) {
    // modified: Camilo Calderon 04 AUG 2014
    string directory;
    double Temperature, DOS_Emin, DOS_Emax, DOSSCALE; 
    pocc::POCC_SetOptions(options, directory, Temperature, DOS_Emin, DOS_Emax, DOSSCALE);
    vector<vector<double> > DOS, PDOS;  
    vector<double> vEfermi,vprob, Egap; double mag = 0.0, Egap_net=0.0;
    POCC_GENERATE_DOSDATA(directory,Temperature,DOS,PDOS,vEfermi,mag,Egap_net,Egap,vprob);
    // [NON_NECESSARY] DOS = aurostd::ShiftFirstColumn(DOS, max(vEfermi)); //shift to Efermi
    //  cout << "Egap.size()=" << Egap.size() << endl;
    if(Egap.size() == 0) {  // CAMILO, can Egap have size = 0 ?? this is the <xOUTCAR.Egap> with vprob
      cout << "Temperature: " << Temperature 
           << "  Egap_net:  " << Egap_net << endl ;
    }
    if(Egap.size() == 1) {
      cout << "Temperature: " << Temperature 
           << "  Egap up:  " << Egap.at(0) 
           << "  Egap_net:  " << Egap_net << endl ;
    }
    if(Egap.size() == 2) {
      cout << "Temperature: " << Temperature 
           << "  Egap up:  " << Egap.at(0) 
           << "  Egap dn:  " << Egap.at(1) 
           << "  Egap net: " << Egap_net << endl ;
    }
  }
}

// ***************************************************************************
// pocc::POCC_MAG(string options)
// ***************************************************************************
namespace pocc {
  void POCC_MAG(string options) {
    string directory;
    double T, DOS_Emin, DOS_Emax, DOSSCALE; 
    pocc::POCC_SetOptions(options, directory, T, DOS_Emin, DOS_Emax, DOSSCALE);
    vector<vector<double> > DOS, PDOS;  vector<double> vEfermi, vprob,Egap; double mag = 0.0,Egap_net=0.0;
    POCC_GENERATE_DOSDATA(directory,T,DOS,PDOS,vEfermi,mag,Egap_net,Egap,vprob);
    cout << "magnetic moment: "  << mag << endl;
  }
}

#endif //  _AFLOW_POCCUPATION_EDOS_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
