// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2013               *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_kesong.cpp
// functions written by KESONG YANG
// 2010-2013: kesong.yang@gmail.com

#ifndef _AFLOW_CONTRIB_KESONG_MULTIENUM_CPP_
#define _AFLOW_CONTRIB_KESONG_MULTIENUM_CPP_

#include "aflow_contrib_kesong.h"

// ***************************************************************************
// pocc::POSCAR2ENUM(istream& input)
// ***************************************************************************
namespace pocc {
  void POSCAR2ENUM(istream& input) {

    xstructure xstr_in; //, xstr_out;
    xstr_in.Clear();
    xstr_in=xstructure(input, IOVASP_POSCAR);
    xstr_in.ReScale(1.00000);

    stringstream oss;
    oss.str("");
    ofstream FileMESSAGE;
    FileMESSAGE.open("LOG.ENUM.INPUT");
    _aflags aflags;

    ostringstream aus;
    aus <<"0000 MESSAGE    Printing input POSCAR " << Message(aflags, "user, host, time");
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus << xstr_in;
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);

    pocc::POSCAR2ENUM(xstr_in, oss, FileMESSAGE, aflags);
    aurostd::stringstream2file(oss,"struct_enum.in");
    FileMESSAGE.close();
  }
}

// ***************************************************************************
// pocc::POSCAR2ENUM(xstructure &a, stringstream &oss, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
namespace pocc {
  void POSCAR2ENUM(xstructure &a, stringstream &oss, ofstream &FileMESSAGE, _aflags &aflags) {
    xstructure xstr_in=a;

    string str_species_ori = xstr_in.SpeciesString();
    vector<string> vxstr_species_ori;
    aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");

    int nHNF = InitializeXstr(xstr_in, vxstr_species_ori, FileMESSAGE, aflags);
    int nMin = nHNF;
    int nMax = nHNF;

    //OptimizeXstr(xstr_in, FileMESSAGE, aflags);

    oss.precision(6);
    oss.setf(ios_base::fixed,ios_base::floatfield);
    oss << xstr_in.title << endl;
    oss << "bulk" << endl;
    oss <<xstr_in.lattice(1,1)<<"  "<<xstr_in.lattice(1,2) <<"  "<<xstr_in.lattice(1,3)<<endl;
    oss <<xstr_in.lattice(2,1)<<"  "<<xstr_in.lattice(2,2) <<"  "<<xstr_in.lattice(2,3)<<endl;
    oss <<xstr_in.lattice(3,1)<<"  "<<xstr_in.lattice(3,2) <<"  "<<xstr_in.lattice(3,3)<<endl;

    vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr_in);
    int nDFull=NormalisedNumber.size();
    int Ntype;
    if(CheckVacancy(xstr_in)) {Ntype=xstr_in.num_each_type.size()+1;}
    else {Ntype=xstr_in.num_each_type.size();}

    oss << Ntype << " -nary case" << endl;
    oss << nDFull << "  #Number of points in the multilattice" << endl;

    vector<vector<int> > LabelXstr=CalculateLableXstructure(xstr_in);
    for (unsigned int i=0; i<NormalisedNumber.size();i++) {
      int j=NormalisedNumber.at(i).at(0);
      oss << xstr_in.atoms.at(j).cpos[1] << "   " << xstr_in.atoms.at(j).cpos[2] << "   " << xstr_in.atoms.at(j).cpos[3] << "   ";
      for (unsigned int k=0; k<LabelXstr.at(i).size();k++) { 
	oss << LabelXstr.at(i).at(k);
	if(LabelXstr.at(i).size()-k-1) {
	  oss << "/";
	}
      }
      oss << endl;
    }

    vector<vector<int> > cRangeFrac=CalculatecRange(xstr_in);
    int NumCell=cRangeFrac.at(0).at(2)/nDFull;
    if(NumCell>nMax) {nMax=NumCell;}
    if(!CheckPartialOccupation(xstr_in)) {nMax=1;}
    oss << nMin << "  " << nMax << "  # Starting and ending cell sizes for search" << endl;
    oss << "1E-6   # Epsilon (finite precision parameter)" << endl;
    oss << "full list of labelings" << endl;
    oss << "# Concentration ranges" << endl;
    for (int i=0; i<Ntype;i++) {
      oss << cRangeFrac.at(i).at(0) << " " << cRangeFrac.at(i).at(1) << " " << cRangeFrac.at(i).at(2) << endl;
    }

    ostringstream aus;
    aus.str("");
    aus <<"0000 MESSAGE    Printing \'struct_enum.in\' file " << Message(aflags, "user, host, time");
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus << oss.str();
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  }
}

// ***************************************************************************
// pocc::POSCAR2ENUM(xstructure &a)
// ***************************************************************************
namespace pocc {
  string POSCAR2ENUM(xstructure &a) {
    stringstream oss;
    oss.str("");
    xstructure xstr_in=a;
    string str_species_ori = xstr_in.SpeciesString();
    vector<string> vxstr_species_ori;
    aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");

    ofstream FileMESSAGE;
    _aflags aflags;
    int nHNF = InitializeXstr(xstr_in, vxstr_species_ori, FileMESSAGE, aflags);
    int nMin = nHNF;
    int nMax = nHNF;

    //OptimizeXstr(xstr_in, FileMESSAGE, aflags);

    oss.precision(6);
    oss.setf(ios_base::fixed,ios_base::floatfield);
    oss << xstr_in.title << endl;
    oss << "bulk" << endl;
    oss <<xstr_in.lattice(1,1)<<"  "<<xstr_in.lattice(1,2) <<"  "<<xstr_in.lattice(1,3)<<endl;
    oss <<xstr_in.lattice(2,1)<<"  "<<xstr_in.lattice(2,2) <<"  "<<xstr_in.lattice(2,3)<<endl;
    oss <<xstr_in.lattice(3,1)<<"  "<<xstr_in.lattice(3,2) <<"  "<<xstr_in.lattice(3,3)<<endl;

    vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr_in);
    int nDFull=NormalisedNumber.size();
    int Ntype;
    if(CheckVacancy(xstr_in)) {Ntype=xstr_in.num_each_type.size()+1;}
    else {Ntype=xstr_in.num_each_type.size();}

    oss << Ntype << " -nary case" << endl;
    oss << nDFull << "  #Number of points in the multilattice" << endl;

    vector<vector<int> > LabelXstr=CalculateLableXstructure(xstr_in);
    for (unsigned int i=0; i<NormalisedNumber.size();i++) {
      int j=NormalisedNumber.at(i).at(0);
      oss << xstr_in.atoms.at(j).cpos[1] << "   " << xstr_in.atoms.at(j).cpos[2] << "   " << xstr_in.atoms.at(j).cpos[3] << "   ";
      for (unsigned int k=0; k<LabelXstr.at(i).size();k++) { 
	oss << LabelXstr.at(i).at(k);
	if(LabelXstr.at(i).size()-k-1) {
	  oss << "/";
	}
      }
      oss << endl;
    }

    vector<vector<int> > cRangeFrac=CalculatecRange(xstr_in);
    int NumCell=cRangeFrac.at(0).at(2)/nDFull;
    if(NumCell>nMax) {nMax=NumCell;}
    if(!CheckPartialOccupation(xstr_in)) {nMax=1;}
    oss << nMin << "  " << nMax << "  # Starting and ending cell sizes for search" << endl;
    oss << "1E-6   # Epsilon (finite precision parameter)" << endl;
    oss << "full list of labelings" << endl;
    oss << "# Concentration ranges" << endl;
    for (int i=0; i<Ntype;i++) {
      oss << cRangeFrac.at(i).at(0) << " " << cRangeFrac.at(i).at(1) << " " << cRangeFrac.at(i).at(2) << endl;
    }
    return oss.str();
  }
}


// ***************************************************************************
// pocc::MultienumPrintAllXstr(istream& input)
// ***************************************************************************
namespace pocc {
  bool MultienumPrintAllXstr(istream& input) {
    xstructure xstr(input,IOVASP_AUTO);
    xstr.ReScale(1.00000);

    ofstream FileMESSAGE;
    FileMESSAGE.open("LOG.ENUM");
    _aflags aflags;
    ostringstream aus;
    aus << "0000 MESSAGE    Printing input POSCAR " << Message(aflags, "user, host, time");
    aus << AFLOWIN_SEPARATION_LINE<< endl;
    aus << xstr;
    aus << AFLOWIN_SEPARATION_LINE<< endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);

    //OptimizeXstr(xstr, FileMESSAGE, aflags);

    stringstream ss;
    stringstream oss;
    vector<xstructure> groupxstr;
#ifdef _AFLOW_GUS_POCC_
    groupxstr=pocc::MultienumGenerateXstr(xstr, FileMESSAGE, aflags);
#endif
    for (uint i=0; i<groupxstr.size();i++) {
      ss.str("");
      ss << setfill('0') << setw(6) <<(i+1);
      oss << AFLOWIN_SEPARATION_LINE<< endl;
      oss << "[VASP_POSCAR_MODE_EXPLICIT]START." <<ss.str() << endl;
      oss << groupxstr.at(i);
      oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP." <<ss.str() << endl;
      oss << AFLOWIN_SEPARATION_LINE<< endl;
    }
    aus << "0000 MESSAGE    Printing derivate POSCARs " << Message(aflags, "user, host, time");
    aus << oss.str();
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    FileMESSAGE.close();
    return TRUE;
  }
}

// ***************************************************************************
// vector<xstructure> pocc::MultienumGenerateXstr(xstructure& xstr, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
extern "C" {
#ifdef _AFLOW_GUS_POCC_
  void __aflow_call_enum_MOD_aflow_pass_parameter(double *parLV, int *nDFull,double *dFull,int *rdFull,int *cdFull,int *k, int *nMin, int *nMax, char *pLatTyp, const double *eps, bool *full, int *labelFull, int *rlabelFull, int *clabelFull, int *digitFull, int *ndigitFull, int *equivalencies, int *nequivalencies,bool *conc_check,int *cRange, int *rcRange, int *ccRange);
  void __aflow_call_makestr_MOD_makestr(int *strNi, int *strNf);
#endif
}

#ifdef _AFLOW_GUS_POCC_
namespace pocc {
  vector<xstructure> MultienumGenerateXstr(xstructure& xstr, ofstream &FileMESSAGE, _aflags &aflags) {
    const double eps=1E-6;
    vector<string> vectitle;
    string title;
    aurostd::string2tokens(xstr.title, vectitle, " ");
    title=vectitle.at(0);
    
    string str_species_ori = xstr.SpeciesString();
    vector<string> vxstr_species_ori;
    aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");
  
    int nHNF = InitializeXstr(xstr, vxstr_species_ori,FileMESSAGE, aflags);
    int nMin = nHNF;
    int nMax = nHNF;


    ostringstream aus;
    aus << "0000 MESSAGE    Print multienum input file " << Message(aflags, "user, host, time");
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aus << pocc::POSCAR2ENUM(xstr);
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    double parLV[3*3];
    parLV[0]=xstr.lattice(1,1); parLV[1]=xstr.lattice(1,2); parLV[2]=xstr.lattice(1,3);
    parLV[3]=xstr.lattice(2,1); parLV[4]=xstr.lattice(2,2); parLV[5]=xstr.lattice(2,3);
    parLV[6]=xstr.lattice(3,1); parLV[7]=xstr.lattice(3,2); parLV[8]=xstr.lattice(3,3);

    vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
    int nDFull=NormalisedNumber.size();

    int Ntype;
    bool VacancyFLAG=CheckVacancy(xstr);
    if(VacancyFLAG) {Ntype=xstr.num_each_type.size()+1;} //Ntype means the number of types of atoms
    else {Ntype=xstr.num_each_type.size();}
    int k=Ntype;

    char pLatTyp[1];
    pLatTyp[0]='B';

    bool full=TRUE;
    bool conc_check=TRUE;

    int nequivalencies=nDFull;
    int equivalencies[nDFull];
    for (int i=0; i<nDFull; i++) {
      equivalencies[i]=i+1;
    }

    int rdFull=3;                  //rdFull means the row number of dFull
    int cdFull=nDFull;             //cdFull means the column number of dFull
    double dFull[rdFull*cdFull];   //dFull is the data array storing the atom sites
    for (unsigned int i=0; i<NormalisedNumber.size();i++) {
      int j=NormalisedNumber.at(i).at(0);
      dFull[3*i+0]=xstr.atoms.at(j).cpos[1];
      dFull[3*i+1]=xstr.atoms.at(j).cpos[2];
      dFull[3*i+2]=xstr.atoms.at(j).cpos[3];
    }

    unsigned int MaxNum=Ntype;
    vector<vector<int> > LabelXstr=CalculateLableXstructure(xstr);
    vector<vector<int> > NewLabelXstr; //This label fills the data array
    vector<int> LabelTmp;
    for (unsigned int i=0; i<LabelXstr.size();i++) {
      LabelTmp.clear();
      for (unsigned int j=0; j<LabelXstr.at(i).size();j++) {
	LabelTmp.push_back(LabelXstr.at(i).at(j));
      }
      if(LabelXstr.at(i).size()<MaxNum) {
	for (unsigned int k=0; k<(MaxNum-LabelXstr.at(i).size());k++) {
	  LabelTmp.push_back(-1);
	}
      }
      NewLabelXstr.push_back(LabelTmp);
    }
    //-----------------------------------------------------------------------
    //Convert the multi dimensional array into one dimensional
    vector<int> OneDimensionLabelXstr;
    for (unsigned int i=0; i<NewLabelXstr.size();i++) {
      for (unsigned int j=0; j<NewLabelXstr.at(i).size();j++) {
	OneDimensionLabelXstr.push_back(NewLabelXstr.at(i).at(j));
      }
    }
    //-----------------------------------------------------------------------
    int clabelFull=nDFull;    
    int rlabelFull=MaxNum;         //rlabelFull means number of columns, different from c++
    int LabelSize=clabelFull*rlabelFull;
    int labelFull[LabelSize];
    for (int i=0;i<LabelSize;i++) {
      labelFull[i]=OneDimensionLabelXstr.at(i);
    }

    int ndigitFull=nDFull;
    int digitFull[nDFull];
    for(int i=0; i<nDFull; i++) {
      digitFull[i]=LabelXstr.at(i).size();
    }

    //Note the format of cRange is different from that of dFull, if you take a look at its fortran code
    vector<vector<int> > cRangeFrac=CalculatecRange(xstr);
    int ccRange=3;
    int rcRange=Ntype;
    int cRange[rcRange*ccRange];
    int cRange2[rcRange*ccRange];
    for (int i=0; i<Ntype;i++) {
      cRange2[i*3+0]=cRangeFrac.at(i).at(0);
      cRange2[i*3+1]=cRangeFrac.at(i).at(1);
      cRange2[i*3+2]=cRangeFrac.at(i).at(2);
    }
    for( int i=0; i<3; i++) {
      for (int j=0; j<Ntype; j++) {
	cRange[i*Ntype+j]=cRange2[j*3+i];
      }
    }
    int NumCell=cRangeFrac.at(0).at(2)/nDFull;
    if(NumCell>nMax) {nMax=NumCell;}
    if(!CheckPartialOccupation(xstr)) {nMax=1;}
    char *cur_dir_name = getcwd(NULL, 0);
    string tmpdir=aurostd::TmpDirectoryCreate("PartialOccupation");
    char new_dir_name[1024];
    strcpy(new_dir_name, tmpdir.c_str());
    chdir(new_dir_name);

    __aflow_call_enum_MOD_aflow_pass_parameter(parLV, &nDFull,dFull,&rdFull,&cdFull,&k, &nMin, &nMax, pLatTyp,&eps, &full, labelFull,&rlabelFull,&clabelFull,digitFull,&ndigitFull,equivalencies,&nequivalencies,&conc_check,cRange,&rcRange,&ccRange);

    //debug
    //--------------------------------------------------------------------------------
    //Call makestr
    stringstream str_out;
    string filename_str_out="struct_enum.out";
    if(!aurostd::FileExist(filename_str_out)) {
      cerr << "File \"struct_enum.out\" does not exist! " << endl;
      exit(1);
    }
    string line,lastline;
    file2stringstream(filename_str_out, str_out);
    while(!str_out.eof()) {
      line.clear();
      getline(str_out,line);
      if(line.size()>0) {
	lastline=line;
      }
    }
    int strNi;
    int strNf;
    vector<string> vec_str_out;
    aurostd::string2tokens(lastline, vec_str_out, " ");
    if(is_number(vec_str_out.at(0))) {
      strNi=1;
      strNf=aurostd::string2utype<int>(vec_str_out.at(0));
    }
    else {
      cerr << "Error! Multienum did not generate a valid struct_enum.out file" << endl;
      exit(1);
    }

    ////output
    aus << "0000 MESSAGE    Print multienum output file " << Message(aflags, "user, host, time");
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    aus << str_out.str();
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    //output

    __aflow_call_makestr_MOD_makestr(&strNi, &strNf);

    stringstream ss;
    xstructure xstr_cleaned;
    vector<xstructure> groupxstr, groupxstr_names;
    for (int i=1; i<=strNf;i++) {
      ss.str("");
      ss << setfill('0') << setw(6) << i;
      string strnumber=ss.str();

      xstructure xstr_partocc="xstr"+strnumber;
      string vaspname="vasp."+strnumber;
      xstr_partocc=xstructure(vaspname, IOVASP_POSCAR);
      string_replaceAll(xstr_partocc.title,"MULTIENUM", title); //Replace the name of POSCAR using the actual title
      //if xstr has vacancy, then we need to clean them
      if(VacancyFLAG) {
	xstr_cleaned=CleanVacancy(xstr_partocc);
	groupxstr.push_back(xstr_cleaned);
      }
      else {
	groupxstr.push_back(xstr_partocc);
      }
    }
    //To make sure safely remove the directory, we divide it into two steps
    //Clean all the files in the temporary directory
    stringstream ss_cmd;
    ss_cmd << "rm -rf " << tmpdir << endl;
    aurostd::execute(ss_cmd);
    
    chdir(cur_dir_name);
    return groupxstr;
  }
}
#endif

// ***************************************************************************
// pocc::MULTIENUM(vector<string> argv,istream& input)
// ***************************************************************************
// pocc::MultienumFilterXstr(vector<xstructure> groupxstr, vector<string> AtomSpecies) 
namespace pocc {
  bool MultienumPrintSortedXstr(istream& input) {
    xstructure xstr; 
    xstr.Clear();
    xstr=xstructure(input, IOVASP_POSCAR);

    vector<string> AtomSpecies;
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");

    if(!xstr.atoms.at(0).name_is_given) {
      cerr << "Warnning!!! Atom Species are not assigned!" << endl;
      exit(1);
    }

    ofstream FileMESSAGE;
    FileMESSAGE.open("LOG.ENUM.SORTED");
    _aflags aflags;
    ostringstream aus;
    OptimizeXstr(xstr, FileMESSAGE, aflags);
    
    vector<xstructure> groupxstr;
#ifdef _AFLOW_GUS_POCC_
    groupxstr=pocc::MultienumGenerateXstr(xstr, FileMESSAGE, aflags);
#endif

    vector<xstructure> sortedgroupxstr;
    //sortedgroupxstr=pocc::SortGroupXstr(groupxstr, AtomSpecies);

    vector<xstructure> groupxstr_names = AssignNameXstr(groupxstr, AtomSpecies);
    sortedgroupxstr=pocc::SortGroupXstrUFFEnergy(groupxstr_names);

    //vector<xstructure> vxstr_final = sortedgroupxstr;
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(sortedgroupxstr, FileMESSAGE, aflags);
    //Format xstructure
    vector<xstructure> vxstr_final_alphabetic;
    for (uint i=0; i<vxstr_final.size();i++) {
      xstructure xstr_tmp = vxstr_final.at(i);
      vector<string> vxstr_tmp;
      aurostd::string2tokens(xstr_tmp.SpeciesString(), vxstr_tmp, " ");
      xstructure xstr2 = AssignNameXstr(xstr_tmp, vxstr_tmp);
      xstr2.SpeciesPutAlphabetic();
      vxstr_final_alphabetic.push_back(xstr2);
    }

    stringstream ss;
    stringstream oss;
    for (uint i=0; i<vxstr_final_alphabetic.size();i++) {
      ss.str("");
      ss << setfill('0') << setw(6) <<(i+1);
      oss << AFLOWIN_SEPARATION_LINE<< endl;
      oss << "[VASP_POSCAR_MODE_EXPLICIT]START." <<ss.str() << endl;
      oss << vxstr_final_alphabetic.at(i);
      oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP." <<ss.str() << endl;
      oss << AFLOWIN_SEPARATION_LINE<< endl;
    }
    aus << "0000 MESSAGE    Printing sorted derivate POSCARs " << Message(aflags, "user, host, time");
    aus << oss.str();
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    /*
    //This test proves that hfncell program can exactly give same number supercells with multienum
    //if they are different, this is due to the energy tollerance in gulp.
    cout << "total number " << vxstr_final.size() << endl;
    */
    FileMESSAGE.close();
    
    return TRUE;
  }
}
// ***************************************************************************
#endif
// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG - Duke University 2010-2013              *
// *                                                                         *
// ***************************************************************************

