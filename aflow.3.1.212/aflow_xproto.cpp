// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2008-2015
// this is one of the most complicate parts of aflow
// fix Ag1Te3_ICSD_37186

///home/auro/work/AFLOW3/aflow --potential=potpaw_GGA --ldau2 --potential_complete --aflow_proto=ICSD_52482.A,ICSD_52483.A,ICSD_76031.A,ICSD_652630.A,ICSD_652632.A,ICSD_652633.A,ICSD_652635.A,ICSD_652637.A,ICSD_652640.A,ICSD_108746.A,ICSD_246657.A:Sm_2,Sm_3

// /home/auro/work/AFLOW3/aflow --potential=potpaw_GGA  --potential_complete --aflow_proto=A3,ICSD_94429.A,ICSD_108026.A,ICSD_165132.A,ICSD_240995.A,ICSD_22300,A,ICSD_240995.A,ICSD_14288.A,ICSD_43431.A:B,B_s,B_h


#ifndef _AFLOW_XPROTO_CPP
#define _AFLOW_XPROTO_CPP
#include "aflow.h"
#include "aflow_xproto_library_default.cpp"

#define _EPS_ 0.02
//#define _XPROTO_TOO_CLOSE_ERROR_ 0.60 // was 0.75 // CO 171023 - see aflow.h
#define AFLOWLIB_SERVER_DEFAULT string("aflowlib.duke.edu")

using aurostd::isdifferent;
using aurostd::isequal;
using spacegroup::SpaceGroupOptionRequired;
using std::ostream;
using std::vector;
using std::deque;
using aurostd::ExtractToStringEXPLICIT;

// [OBSOLETE] std::string Library_ICSD;
// [OBSOLETE] std::vector<std::string> vLibrary_ICSD_ALL;

string* LOAD_Library_ICSD(string file);

namespace aurostd {
  template<class utype> void swap(utype &a,utype &b) {utype temp=a;a=b;b=temp;}
}

// ***************************************************************************
// Function extra operator << for vector
template<class utype>                            // operator <<  vector<>
std::ostream& operator<< (std::ostream& buf,const std::vector<utype>& x) {
  for(uint i=0;i<x.size();i++) {
    buf << x.at(i) << " ";
  }
  return buf;
}
// ***************************************************************************
// Function extra operator << for deque
template<class utype>                            // operator <<  deque<>
std::ostream& operator<< (std::ostream& buf,const std::deque<utype>& x) {
  for(uint i=0;i<x.size();i++) {
    buf << x.at(i) << " ";
  }
  return buf;
}

// ***************************************************************************
// ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL
// ***************************************************************************

// ***************************************************************************
namespace aflowlib {
  string PrototypeCleanLatticeString(const string& latticeIN) {
    string lattice=latticeIN;
    lattice=aurostd::RemoveSubStringFirst(lattice,"fcc");
    lattice=aurostd::RemoveSubStringFirst(lattice,"FCC");
    lattice=aurostd::RemoveSubStringFirst(lattice,"bcc");
    lattice=aurostd::RemoveSubStringFirst(lattice,"BCC");
    lattice=aurostd::RemoveSubStringFirst(lattice,"hcp");
    lattice=aurostd::RemoveSubStringFirst(lattice,"HCP");
    lattice=aurostd::RemoveSubStringFirst(lattice,"sc");
    lattice=aurostd::RemoveSubStringFirst(lattice,"SC");
    lattice=aurostd::RemoveSubStringFirst(lattice,"f");
    lattice=aurostd::RemoveSubStringFirst(lattice,"F");
    lattice=aurostd::RemoveSubStringFirst(lattice,"b");
    lattice=aurostd::RemoveSubStringFirst(lattice,"B");
    lattice=aurostd::RemoveSubStringFirst(lattice,"h");
    lattice=aurostd::RemoveSubStringFirst(lattice,"H");
    lattice=aurostd::RemoveSubStringFirst(lattice,"c");
    lattice=aurostd::RemoveSubStringFirst(lattice,"C");
    lattice=aurostd::RemoveSubStringFirst(lattice,"p");
    lattice=aurostd::RemoveSubStringFirst(lattice,"P");
    lattice=aurostd::RemoveSubStringFirst(lattice,"s");
    lattice=aurostd::RemoveSubStringFirst(lattice,"S");
    lattice=aurostd::RemoveSubStringFirst(lattice," ");
    return lattice;
  }
} // namespace aflowlib

// ***************************************************************************
// PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE
namespace aflowlib {
  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters,string atomA) {
    return aflowlib::PrototypePureHTQC(FileMESSAGE,label,parameters,atomA);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters) {
    return aflowlib::PrototypePureHTQC(FileMESSAGE,label,parameters);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters,string atomA,double volumeA) {
    return aflowlib::PrototypePureHTQC(FileMESSAGE,label,parameters,atomA,volumeA);
  }
} // namespace aflowlib

// ***************************************************************************
// PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE PURE

namespace aflowlib {
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters,string atomA) {
    double atomvolumeA;
    atomvolumeA=GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomA));
    return aflowlib::PrototypePureHTQC(FileMESSAGE,label,parameters,atomA,atomvolumeA);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters) {
    return aflowlib::PrototypePureHTQC(FileMESSAGE,label,parameters,"A",1.0);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters,string atomA,double volumeA) {
    xstructure str("");str.lattice.clear();

    // *********************************************************************
    // FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC FCC
    // *********************************************************************
    if(label=="FCC" || label=="fcc" || label=="1" || label=="2" || label=="A1") {
      deque<string> atomX;atomX.push_back(atomA);atomX.push_back("B");
      deque<double> volumeX;volumeX.push_back(volumeA);volumeX.push_back(0.0);
      return aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);}
    // *********************************************************************
    // BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC BCC
    // *********************************************************************
    if(label=="BCC" || label=="bcc" || label=="58" || label=="59" || label=="A2") {
      deque<string> atomX;atomX.push_back(atomA);atomX.push_back("B");
      deque<double> volumeX;volumeX.push_back(volumeA);volumeX.push_back(0.0);
      return aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);}
    // *********************************************************************
    // HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP HCP
    // *********************************************************************
    if(label=="HCP" || label=="hcp" || label=="115" || label=="117" || label=="A3") {
      deque<string> atomX;atomX.push_back(atomA);atomX.push_back("B");
      deque<double> volumeX;volumeX.push_back(volumeA);volumeX.push_back(0.0);
      return aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);}
    // *********************************************************************
    // DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA DIA
    // *********************************************************************
    if(label=="DIA" || label=="dia" || label=="301" || label=="302" || label=="A4") {
      deque<string> atomX;atomX.push_back(atomA);atomX.push_back("B");
      deque<double> volumeX;volumeX.push_back(volumeA);volumeX.push_back(0.0);
      return aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);}
    // *********************************************************************
    // INDIUM INDIUM INDIUM INDIUM INDIUM INDIUM INDIUM INDIUM INDIUM IND
    // *********************************************************************
    if(label=="Indium" || label=="In" || label=="303" || label=="304" || label=="A6") {
      deque<string> atomX;atomX.push_back(atomA);atomX.push_back("B");
      deque<double> volumeX;volumeX.push_back(volumeA);volumeX.push_back(0.0);
      return aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);}
    // *********************************************************************
    // betaSn betaSn betaSn betaSn betaSn betaSn betaSn betaSn betaSn
    // *********************************************************************
    if(label=="betaSn" || label=="Sn" || label=="305" || label=="306" || label=="A5") {
      deque<string> atomX;atomX.push_back(atomA);atomX.push_back("B");
      deque<double> volumeX;volumeX.push_back(volumeA);volumeX.push_back(0.0);
      return aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);}
    // *********************************************************************
    // alphaAs alphaAs alphaAs alphaAs alphaAs alphaAs alphaAs alphaAs
    // *********************************************************************
    if(0)  if(label=="alphaAs" || label=="As" || label=="307" || label=="308" || label=="A7") {
	deque<string> atomX;atomX.push_back(atomA);atomX.push_back("B");
	deque<double> volumeX;volumeX.push_back(volumeA);volumeX.push_back(0.0);
	return aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);}
    // *********************************************************************
    return str; // empty
    // *********************************************************************
  }
} // namespace aflowlib

// ***************************************************************************
// ***************************************************************************
namespace aflowlib {
  uint PrototypeLibrariesSpeciesNumber(const string& label) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesSpeciesNumber: label=" << label << endl;
    // search for _ICSD_XXX or ICSD_XXX
    if(aurostd::substring2bool(label,"ICSD_")) { // found ICSD
      vector<string> vline,tokens;
      string _label=label;
      if(aurostd::substring2bool(_label,".")) {
	aurostd::string2tokens(_label,tokens,".");
	_label="";
	for(uint i=0;i<tokens.size()-1;i++)
          _label+=tokens.at(i)+(i<tokens.size()-2?".":"");  //CO 180705 - simply removing .ABC and rebuilding without last .
      } //   cerr << _label << endl;exit(0);
      aurostd::string2vectorstring(string(init::InitGlobalObject("ICSD_List_txt")),vline);
      for(uint iline=0;iline<vline.size();iline++) {
	if(aurostd::substring2bool(vline.at(iline),_label)) {
	  aurostd::string2tokens(vline.at(iline),tokens," ");
	  for(uint itoken=0;itoken<tokens.size();itoken++)
	    if(tokens.at(itoken)==_label) {
	      aurostd::StringSubst(tokens.back(),"(","");
	      aurostd::StringSubst(tokens.back(),")","");
	      return aurostd::string2utype<uint>(tokens.back());
	    }
	}
      }
    }
    //CO 180705 - START
    if(aurostd::substring2bool(label,"_ICSD_")) { // if we get here, then it was not found in our current ICSD_List_txt, could be old, so derive based in compond input
      vector<string> parts;
      aurostd::string2tokens(label,parts,"_");
      if(parts.size()==3){  //compound_ICSD_XXXXXX
        string compound=parts[0];
        vector<string> valloy;
        KBIN::VASP_SplitAlloySpecies(compound,valloy);
        if(valloy.size()>0){return valloy.size();}
      }
    }
    //CO 180705 - START
    if(LDEBUG) { cerr << "uint aflowlib::PrototypeLibrariesSpeciesNumber(const string& label)" << endl; }
    // check ternaries 
    if(label[0]=='T' || label[0]=='t') {return (uint) 3;}
    // check quaternaries
    if(label[0]=='Q' || label[0]=='q') {return (uint) 4;}
    // check unaries
    if(label=="FCC" || label=="fcc" || label=="A1") {return (uint) 1;}
    if(label=="BCC" || label=="bcc" || label=="A2") {return (uint) 1;}
    if(label=="HCP" || label=="hcp" || label=="A3") {return (uint) 1;}
    if(label=="HCP" || label=="hcp" || label=="A3") {return (uint) 1;}
    if(label=="A1" || label=="A2" || label=="A3" || label=="A4" || label=="A5") {return (uint) 1;}
    if(label=="A6" || label=="A7" || label=="A8" || label=="A9" || label=="A10") {return (uint) 1;}
    if(label=="A11" || label=="A12" || label=="A13" || label=="A14" || label=="A15") {return (uint) 1;}
    // check ANRL
    vector<string> tokens;
    aurostd::string2tokens(label,tokens,"_");
    if(tokens.size()>=4) return (uint) tokens.size()-3;
    
    // otherwise binary
    return (uint) 2;
  }
} // namespace aflowlib

// ***************************************************************************
// ***************************************************************************
namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,int mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,int mode)" << endl; }
    deque<string> atomX;
    deque<double> volumeX;
    uint nspeciesHTQC=aflowlib::PrototypeLibrariesSpeciesNumber(label);
    for(uint i=0;i<nspeciesHTQC;i++) {
      string string_tmp="A";string_tmp[0]+=i;
      atomX.push_back(string_tmp);
      volumeX.push_back(1.0);
    }
    return aflowlib::PrototypeLibraries(oss,label,parameters,atomX,volumeX,-2.0,mode,(bool) FALSE);
  }
} // namespace aflowlib

namespace aflowlib { 
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX,int mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX)" << endl; }
    deque<double> volumeX;
    for(uint i=0;i<atomX.size();i++)
      volumeX.push_back(GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomX.at(i))));
    return aflowlib::PrototypeLibraries(oss,label,parameters,atomX,volumeX,-3.0,mode,(bool) FALSE);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX,deque<double> &volumeX,int mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX,deque<double> &volumeX)" << endl; }
    return aflowlib::PrototypeLibraries(oss,label,parameters,atomX,volumeX,-4.0,mode,(bool) FALSE);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX,deque<double> &volumeX,double volume_in,int mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX,deque<double> &volumeX,double volume_in)" << endl; }
    // try GUS
    if(label.length()>1 && !aurostd::substring2bool(label,"_ICSD_") &&
       (label[0]=='F' || label[0]=='f' || label[0]=='B' || label[0]=='b' || label[0]=='H' || label[0]=='h' || label[0]=='S' || label[0]=='s') &&
       (label[1]!='c') && (label[1]!='C')) {
      return aflowlib::PrototypeBinaryGUS(oss,label,atomX.at(0),volumeX.at(0),atomX.at(1),volumeX.at(1),volume_in);   // found FCC,BCC,HCP
    } else {
      return aflowlib::PrototypeLibraries(oss,label,parameters,atomX,volumeX,volume_in,mode,(bool) FALSE);  // try HTQC
    }
  }
} // namespace aflowlib

//double NearestNeighbour(const xstructure &str_in);

// ***************************************************************************
namespace aflowlib {
  void PrototypeFixTET(ostream &oss,xstructure &str,string optionsTET) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string _options=optionsTET;
    //  aurostd::StringSubst(_options,"%","");  cerr << aurostd::string2utype<double>(_options) << endl;
    double fraction=aurostd::string2utype<double>(optionsTET);if(aurostd::substring2bool(optionsTET,"%")) fraction=fraction/100+1;
    if(LDEBUG) { oss << "DEBUG: (aflowlib::PrototypeFixTET) fraction=" << fraction << endl; }
    if(LDEBUG) { oss << "DEBUG: (aflowlib::PrototypeFixTET) optionsTET=" << optionsTET << endl; }
    //  cerr << str.title << endl;exit(0);
    xmatrix<double> mlt(3,3);
    mlt(1,1)=1.0;mlt(2,2)=1.0;mlt(3,3)=fraction;
    if(abs(fraction-1.0)>0.001) {
      str.GetStandardConventional();
      str=GetLTCell(mlt,str);
      // cerr << str.lattice << endl;
      // str.lattice(3,1)=fraction*str.lattice(3,1);
      // str.lattice(3,2)=fraction*str.lattice(3,2);
      // str.lattice(3,3)=fraction*str.lattice(3,3);
      // cerr << str.lattice << endl;
      // DO WE NEED ? str.GetStandardPrimitive();
      // str.SetScale(1.0);
    }
  }
} // namespace aflowlib

// ***************************************************************************
// the mother of all the prototypes
namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option) { // COMPLETE ONE
    // XHOST.DEBUG=TRUE;
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option)" << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: label=" << label << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: parameters=" << parameters << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: vatomX.size()=" << vatomX.size() << endl;}
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: vatomX ="; for(uint i=0;i<vatomX.size();i++) {cerr << " " << vatomX.at(i);} cerr << endl;}
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: vvolumeX.size()=" << vvolumeX.size() << endl;}
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: vvolumeX ="; for(uint i=0;i<vvolumeX.size();i++) {cerr << " " << vvolumeX.at(i);} cerr << endl;}
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: volume_in=" << volume_in << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: mode=" << mode << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: flip_option=" << flip_option << endl; }

    // check for ANRL
    vector<string> vlabel_ANRL;
    if(aurostd::string2tokens(label,vlabel_ANRL,"_")>=4) {
      if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: ANRL=" << 1 << endl; }
//      return anrl::PrototypeANRL(oss,label,parameters,vatomX,vvolumeX,volume_in,mode,flip_option); KESONG
    }
    // done
    
    deque<string> vatomX_backup(vatomX);
    deque<double> vvolumeX_backup(vvolumeX);

    xstructure str("");str.lattice.clear();
    string optionsTET;bool isTET=FALSE;

    if(label.length()>1 && !aurostd::substring2bool(label,"_ICSD_") &&
       (label[0]=='F' || label[0]=='f' || label[0]=='B' || label[0]=='b' || label[0]=='H' || label[0]=='h' || label[0]=='S' || label[0]=='s') &&
       (label[1]!='c') && (label[1]!='C')) {
      str=aflowlib::PrototypeBinaryGUS(oss,label,vatomX.at(0),vvolumeX.at(0),vatomX.at(1),vvolumeX.at(1),volume_in);   // found FCC,BCC,HCP
      return str;
    }

    ostream *voss;if(LDEBUG || XHOST.DEBUG) voss=&cerr; else voss=&oss;
    voss=&cerr;
 
    vector<uint> natomsX;
    string structure_name="";
    str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
    bool found=FALSE;
    std::string label_swap="",index_swap="",label_use="",label_remove="",label_species="",string_tmp;
    string *_pstrstream;
    vector<string> tokens,ttokens,database,speciesX;
    vector<double> vnatomsX;
    stringstream aus;
    ostringstream oaus;
    uint kstructure=0,kprototype=0,kinfo=0,kmode=0,klattice=0,katoms=0,kstop=0;
    uint iat,natoms=0,nspecies=0;
    double volume=0.0;
    int mode_load=0,mode_prim=0,mode_conventional=0,mode_volume=0,mode_remove=0,mode_species=0;
    if(LDEBUG) { *voss << "aflowlib::PrototypeLibraries: XHOST_Library_HTQC.length()=" << XHOST_Library_HTQC.length() << endl; }
    init::InitGlobalObject("Library_HTQC");

    // *********************************************************************
    speciesX.clear();natomsX.clear();
    str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();;str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
    nspecies=0;
    _pstrstream=&XHOST_Library_HTQC; // default is HTQC so it does not complain
    // *********************************************************************
    // PARSE and STORE the file
    // LDEBUG=TRUE;
    // ICSD
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [aflowlib::PrototypeLibraries] [0]" << endl; }
    //  *voss << "MODE=" << mode << endl; exit(0);
    found=FALSE; // start the search
    // *********************************************************************
    if(mode==LIBRARY_MODE_ICSD) {
      // recognize which LIBRARY of ICSD
      if(LDEBUG) {
	if(mode==LIBRARY_MODE_ICSD) {
	  *voss << "DEBUG: (aflowlib::PrototypeLibraries) mode==LIBRARY_MODE_ICSD" << endl;
	}
      }
      
      string icsd_label=label;
      if(aurostd::substring2bool(icsd_label,"_ICSD_")) {
	aurostd::string2tokens(icsd_label,tokens,"_");
	icsd_label=tokens.at(0);
      }
      nspecies=KBIN::VASP_SplitAlloySpecies(aurostd::RemoveSubStringFirst(icsd_label,"C-"),speciesX,vnatomsX);

      //    for(uint i=0;i<speciesX.size();i++) { if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) " << nspecies << " " << speciesX.at(i) << " " << vnatomsX.at(i) << endl; } } // some debug always helps exit(0);
      if(aurostd::substring2bool(LibraryDEFAULT,label)) _pstrstream=&LibraryDEFAULT;
      else {
	if(nspecies<1||nspecies>10) {
	  _pstrstream=LOAD_Library_ICSD("Library_ICSD2");//default
	} else {
	  _pstrstream=LOAD_Library_ICSD("Library_ICSD"+aurostd::utype2string<int>(nspecies));
	}
      }
      if(_pstrstream==NULL) {
	return PrototypeLibraries(oss,label,parameters,vatomX,vvolumeX,volume_in,LIBRARY_MODE_ICSD_AFLOWLIB,flip_option); 
	exit(0);
      }
	
      vatomX.clear();vvolumeX.clear();
      for(uint i=0;i<nspecies;i++) {
	str.species.push_back(speciesX.at(i));
	str.species_pp.push_back(speciesX.at(i));
	str.species_pp_type.push_back("");
	str.species_pp_version.push_back("");
	str.species_pp_ZVAL.push_back(0.0);
	str.species_pp_vLDAU.push_back(deque<double>());
	str.species_volume.push_back(0.0);
	str.species_mass.push_back(0.0);
	vatomX.push_back(speciesX.at(i));
	vvolumeX.push_back(0.0);
	natomsX.push_back(0);
	// str.species.push_back(speciesX.at(i));vatomX.push_back("");vvolumeX.push_back(0.0);natomsX.push_back(0);
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [1] nspecies=" << nspecies << endl; }
    }
    // *********************************************************************
    xstructure str_LIBRARY_MODE_HTQC_ICSD_AFLOWLIB;
    string label_prefix_ICSD_AFLOWLIB,label_postfix_ICSD_AFLOWLIB;
    string label_icsd_ICSD_AFLOWLIB=label;
    if(mode==LIBRARY_MODE_ICSD_AFLOWLIB || mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) {
      if(!XHOST.vflag_control.flag("AFLOWLIB_SERVER") && !XHOST.vflag_control.flag("AFLOWLIB_SERVER"))
        XHOST.vflag_control.push_attached("AFLOWLIB_SERVER",AFLOWLIB_SERVER_DEFAULT); // some default
      // recognize __AFLOWLIB
      if(LDEBUG) { if(mode==LIBRARY_MODE_ICSD_AFLOWLIB) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) mode==LIBRARY_MODE_ICSD_AFLOWLIB" << endl; } }
      if(LDEBUG) { if(mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB" << endl; } }
      vector<string> vlattice;aurostd::string2tokens("BCC,BCT,CUB,FCC,HEX,MCL,MCLC,ORC,ORCC,ORCF,ORCI,RHL,TET,TRI",vlattice,",");
      aflowlib::_aflowlib_entry data;
      vector<aflowlib::_aflowlib_entry> vdata;
      string aurl="";
      *voss << "aflowlib::PrototypeLibraries: Contacting \"" <<  XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER") << "\"" << endl;
      vector<string> tokens,tokensl,tokensd;
      aurostd::string2tokens(label,tokensl,"_");
      if(tokensl.size()>0) label_icsd_ICSD_AFLOWLIB=tokensl.at(tokensl.size()-1);
      aurostd::string2tokens(label_icsd_ICSD_AFLOWLIB,tokensl,".");
      if(tokensl.size()>0) label_icsd_ICSD_AFLOWLIB=tokensl.at(0);
      if(tokensl.size()>1) label_postfix_ICSD_AFLOWLIB=tokensl.at(1);
      
      
      for(uint i=0;i<vlattice.size()&&aurl.empty();i++) {
     	data.clear();
	*voss << "aflowlib::PrototypeLibraries: Scanning " << vlattice.at(i) << " ";
	data.url2aflowlib(XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")+":AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/?format=text",*voss,FALSE);vdata.push_back(data);
	if(LDEBUG) { *voss << "aflowlib::PrototypeLibraries: AFLOWLIB " << vlattice.at(i) << "=" << data.vaflowlib_entries.size() << endl; }
	for(uint j=0;j<data.vaflowlib_entries.size()&&aurl.empty();j++) {
	  if(mode==LIBRARY_MODE_ICSD_AFLOWLIB) {
	    if(data.vaflowlib_entries.at(j)==label) {
	      aurl=XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")+":AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/"+data.vaflowlib_entries.at(j);
	      aurostd::StringSubst(aurl,AFLOWLIB_SERVER_DEFAULT,XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER"));
	      found=TRUE;
	      *voss << "found aurl=" << aurl;
	    }
	  }
	  if(mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) {
	    aurostd::string2tokens(data.vaflowlib_entries.at(j),tokensd,"_");
	    if(tokensd.size()>0 && !label_icsd_ICSD_AFLOWLIB.empty()) {
	      if(label_icsd_ICSD_AFLOWLIB==tokensd.at(tokensd.size()-1)) {
		aurl=XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")+":AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/"+data.vaflowlib_entries.at(j);
		aurostd::StringSubst(aurl,AFLOWLIB_SERVER_DEFAULT,XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER"));
		found=TRUE;
	        *voss << "found aurl=" << aurl;
		label_prefix_ICSD_AFLOWLIB=data.vaflowlib_entries.at(j);
		aurostd::string2tokens(data.vaflowlib_entries.at(j),tokens,"_");
		speciesX.clear();
		nspecies=KBIN::VASP_SplitAlloySpecies(aurostd::RemoveSubStringFirst(tokens.at(0),"C-"),speciesX,vnatomsX);
	      }
	    }
	  }
	}
	if(aurl.empty()) *voss << "[" << data.vaflowlib_entries.size() << " entries]";
	*voss << endl;
      }
      if(LDEBUG) { *voss << "aflowlib::PrototypeLibraries: AFLOWLIB aurl=" << aurl << endl; }
      if(found) {
	stringstream stream;
	xstructure str;
	if(XHOST.vflag_pflow.flag("PROTO::VASP")) {xstructure str_vasp(aurl,"CONTCAR.relax.vasp",IOVASP_AUTO);str=str_vasp;str.title+=" [ICSD from "+XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")+"]" ;}
	if(XHOST.vflag_pflow.flag("PROTO::QE")) {xstructure str_qe(aurl,"CONTCAR.relax.qe",IOQE_AUTO);str=str_qe;str.title+=" [ICSD from "+XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")+"]" ;}
	if(XHOST.vflag_pflow.flag("PROTO::ABINIT")) {stream << aurl << endl;aurostd::url2stringstream(aurl+"/CONTCAR.relax.abinit",stream);cout << stream.str();exit(0);}  // abinit constructor not prepared, yet
	// if(XHOST.vflag_pflow.flag("PROTO::ABINIT")) {xstructure str_abinit(aurl,"CONTCAR.relax.abinit",IOABINIT_AUTO);return str_abinit;}
	if(XHOST.vflag_pflow.flag("PROTO::AIMS")) {stream << aurl << endl;aurostd::url2stringstream(aurl+"/CONTCAR.relax.aims",stream);cout << stream.str();exit(0);}  // aims constructor not prepared, yet
	// if(XHOST.vflag_pflow.flag("PROTO::AIMS")) {xstructure str_aims(aurl,"CONTCAR.relax.aims",IOAIMS_AUTO);return str_aims;}
  //CO 180622 - aflow_proto does not declare PROTO::VASP default for some reason
  if(str.atoms.size()==0){xstructure str_vasp(aurl,"CONTCAR.relax.vasp",IOVASP_AUTO);str=str_vasp;str.title+=" [ICSD from "+XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")+"]" ;} //catch all
	if(mode==LIBRARY_MODE_ICSD_AFLOWLIB) return str;
	if(mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) {
	  str_LIBRARY_MODE_HTQC_ICSD_AFLOWLIB=str;
	}
      }
    }
    // *********************************************************************
    // HTQC and ICSD
    if(mode==LIBRARY_MODE_HTQC_ICSD || mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) {
      //   *voss << "MODE ICSD" << endl;
      string line_ICSD,label_compound;
      string label_prefix=label;
      string label_postfix="",label_icsd="";

      if(LDEBUG) { if(mode==LIBRARY_MODE_HTQC_ICSD) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) mode==LIBRARY_MODE_HTQC_ICSD" << endl; } }
      if(LDEBUG) { if(mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB" << endl; } }

      if(mode==LIBRARY_MODE_HTQC_ICSD) {
	if(aurostd::substring2bool(label_prefix,".")) {
	  aurostd::string2tokens(label_prefix,tokens,".");
	  if(tokens.size()>0) label_prefix=tokens.at(0);
	  if(tokens.size()>0) label_icsd=tokens.at(0);
	  if(tokens.size()>0) label_postfix=tokens.at(1);
	}
	
	if(aurostd::substring2bool(label_prefix,"ICSD_") && !aurostd::substring2bool(label_prefix,"_ICSD_")) {  // LOAD THEM ALL
	  bool found_HTQC_ICSD=FALSE;
	  for(uint n=1;n<=10&&!found_HTQC_ICSD;n++) {
	    _pstrstream=LOAD_Library_ICSD("Library_ICSD"+aurostd::utype2string<uint>(n));
	    if(_pstrstream==NULL) {
	      return PrototypeLibraries(oss,label,parameters,vatomX,vvolumeX,volume_in,LIBRARY_MODE_HTQC_ICSD_AFLOWLIB,flip_option); 
	      exit(0);
	    }
	    aurostd::string2tokens(*_pstrstream,database,"\n");
	    for(uint i=0;i<database.size()&&!found_HTQC_ICSD;i++)
	      if(aurostd::substring2bool(database[i],"STRUCTURE")) {
		aurostd::string2tokens(database[i],tokens);
		for(uint j=0;j<tokens.size();j++)
		  if(label_prefix==tokens.at(j)) {
		    found_HTQC_ICSD=TRUE;//   *voss << database[i] << endl;
		    line_ICSD=database[i+1];
		    for(uint k=0;k<tokens.size();k++)
		      if(aurostd::substring2bool(tokens.at(k),"_ICSD_")) {
			label_prefix=tokens.at(k);
			label_compound=tokens.at(1);
			//    *voss << label_compound << endl;
		      }
		  }
	      }
	  }
	}
  	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) label_prefix=" << label_prefix << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) label_icsd=" << label_icsd << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) label_postfix=" << label_postfix << endl; }
 	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) line_ICSD=" << line_ICSD << endl; }
      }
      
      if(mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) {
	label_prefix=label_prefix_ICSD_AFLOWLIB;
	label_icsd="ICSD_"+label_icsd_ICSD_AFLOWLIB;
	label_postfix=label_postfix_ICSD_AFLOWLIB;
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) label_prefix=" << label_prefix << endl; }
        if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) label_icsd=" << label_icsd << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) label_postfix=" << label_postfix << endl; }
      }
      
      if(label_prefix!="" && (line_ICSD!="" || mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB)) {
	found=TRUE;
	deque<string> vatomX_tmp;
	deque<double> vvolumeX_tmp;
	deque<string> abc;
	if(mode==LIBRARY_MODE_HTQC_ICSD)          {str=aflowlib::PrototypeLibraries(oss,label_prefix,parameters,vatomX_tmp,vvolumeX_tmp,volume_in,LIBRARY_MODE_ICSD,flip_option);}
	if(mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB) {str=str_LIBRARY_MODE_HTQC_ICSD_AFLOWLIB;vatomX_tmp=vatomX;vvolumeX_tmp=vvolumeX;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str=" << str << endl; }
	
	// need to fix this because I could get stuff from WYCKOFIZATION
	
	//  *voss << str << endl;exit(0);
	if(vatomX_tmp.size()!=label_postfix.size()) {
	  *voss << "ERROR - aflowlib::PrototypeLibraries:  your label needs to contain the right species number/order, i.e. " << label_icsd << ".ABC..." << endl;
	  *voss << "ERROR - aflowlib::PrototypeLibraries:  vatomX.size()!=label_postfix.size() (" << vatomX_tmp.size() << "!=" << label_postfix.size() << ")" << endl;
	  exit(0);}
	// create the mapping

	vector<int> mapij;
	for(uint i=0;i<label_postfix.size();i++) {
	  mapij.push_back(-1);
	  string digit=" ";
	  for(uint j=0;j<26;j++) {
	    digit[0]=65+j; // 65 is "A"
	    if(label_postfix.at(i)==digit[0])
	      mapij.at(mapij.size()-1)=j;
	  }
	}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) mapij.size()=" << mapij.size() << endl; for(uint i=0;i<mapij.size();i++) *voss << "mapij.at(" << i << ")=" << mapij.at(i) << endl;}

	
	abc.clear();
	if(vatomX.size()==0 && vvolumeX.size()==0) for(uint i=0;i<vatomX_tmp.size();i++) {string digit=" ";digit[0]=65+mapij.at(i);abc.push_back(digit);}
	if(vatomX.size()==str.species.size() && (vvolumeX.size()==0 || vvolumeX.size()==str.species.size())) for(uint i=0;i<vatomX_tmp.size();i++)  {abc.push_back(vatomX.at(mapij.at(i)));}
	//      *voss << label_prefix << endl;

	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) vatomX.size()=" << vatomX.size() << endl; for(uint i=0;i<vatomX.size();i++) *voss << "vatomX.at(" << i << ")=" << vatomX.at(i) << endl;}
	
	if(vatomX.size()!=vatomX_tmp.size()) {
	  *voss << "EXCEPTION FOUND - aflowlib::PrototypeLibraries:  vatomX.size()!=vatomX_tmp.size() (" << vatomX.size() << "!=" << vatomX_tmp.size() << ") ... FIXING" << endl;
	  vatomX_tmp.clear(); for(uint i=0;i<vatomX.size();i++) vatomX_tmp.push_back(vatomX.at(i));
	}
	// /*
	// if(vvolumeX.size()!=vvolumeX_tmp.size()) {
	// 	*voss << "EXCEPTION FOUND - aflowlib::PrototypeLibraries:  vvolumeX.size()!=vvolumeX_tmp.size() (" << vvolumeX.size() << "!=" << vvolumeX_tmp.size() << ") ... FIXING" << endl;
	// 	//	vvolumeX_tmp.clear(); for(uint i=0;i<vvolumeX.size();i++) vvolumeX_tmp.push_back(vvolumeX.at(i));
	// }
	// */
	if(vatomX.size()>0 && vatomX.size()!=vatomX_tmp.size()) {
	  *voss << "ERROR - aflowlib::PrototypeLibraries:  vatomX.size()!=vatomX_tmp.size() (" << vatomX.size() << "!=" << vatomX_tmp.size() << ")" << endl;
	  exit(0);
	}

	// fil A,B,C ..
	if(vatomX.size()==0 && vvolumeX.size()==0) { // fill A B C
	  str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
	  for(uint i=0;i<vatomX_tmp.size()&&i<abc.size();i++) {
	    str.species.push_back(abc.at(i));str.species_pp.push_back(abc.at(i));str.species_pp_type.push_back("");str.species_pp_version.push_back("");str.species_pp_ZVAL.push_back(0.0);
	    str.species_pp_vLDAU.push_back(deque<double>());str.species_volume.push_back(1.0);str.species_mass.push_back(0.0);
	  }
	}
	// fil real names/volumes
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) vatomX.size()=" << vatomX.size() << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.species.size()=" << str.species.size() << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) vvolumeX.size()=" << vvolumeX.size() << endl; }
 	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.GetVolume()=" << str.GetVolume() << endl; }

	if(vatomX.size()==str.species.size() && (vvolumeX.size()==0 || vvolumeX.size()==str.species.size())) { // fill A B C
	  //	*voss << "IN" << endl;
	  str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
	  for(uint i=0;i<vatomX.size()&&i<abc.size();i++) {
	    str.species.push_back(abc.at(i));str.species_pp.push_back(abc.at(i));str.species_pp_type.push_back("");str.species_pp_version.push_back("");str.species_pp_ZVAL.push_back(0.0);
	    str.species_pp_vLDAU.push_back(deque<double>());str.species_mass.push_back(0.0);
	    str.species_volume.push_back(GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i))));
	    //	  *voss << "GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i)))=" << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i))) << endl;
	  }
	}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) vatomX.size()=" << vatomX.size() << endl; for(uint i=0;i<vatomX.size();i++) *voss << "vatomX.at(" << i << ")=" << vatomX.at(i) << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(vatomX.at(i))) << endl;}
	//	exit(0);
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) abc.size()=" << abc.size() << endl; for(uint i=0;i<abc.size();i++) *voss << "abc.at(" << i << ")=" << abc.at(i) << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i))) << endl;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) volumeX.size()=" << vvolumeX.size() << endl; for(uint i=0;i<vvolumeX.size();i++) *voss << "vvolumeX.at(" << i << ")=" << vvolumeX.at(i) << endl;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.species.size()=" << str.species.size() << endl; for(uint i=0;i<str.species.size();i++) *voss << "str.species.at(" << i << ")=" << str.species.at(i) << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(str.species.at(i))) << endl;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.atoms.size()=" << str.atoms.size() << endl; for(uint i=0;i<str.atoms.size();i++) *voss << "str.atoms.at(" << i << ").type=" << str.atoms.at(i).type << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(str.atoms.at(i).name)) << endl;}
	//

	double volume=0.0;
	for(uint iat=0;iat<str.atoms.size();iat++) {
	  if((uint) str.atoms.at(iat).type<str.species.size()) {
	    str.atoms.at(iat).name=str.species.at(str.atoms.at(iat).type);
	    str.atoms.at(iat).name_is_given=TRUE;
	    volume+=str.species_volume.at(str.atoms.at(iat).type);
	    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) volume=" << volume << "  (" << str.atoms.at(iat).type << ")" << endl; }
	  }
	}
	if(abs(volume)<1.0e-6) { // wrap up to fix it
	  *voss << "EXCEPTION FOUND - aflowlib::PrototypeLibraries:  volume==0 ... FIXING" << endl;
	  for(uint iat=0;iat<vvolumeX.size();iat++) volume+=vvolumeX.at(iat);
	}

	//  if(vatomX.size()==0 && vvolumeX.size()==0)
	{  // fix volume in WYC
	  str.scale=std::pow((double) (abs(volume)/det(str.lattice)),(double) 1.0/3.0);
	  str.neg_scale=TRUE;
	}
	// strip unstrip atoms
	std::deque<_atom> atoms;
	atoms=str.atoms;
	while(str.atoms.size()>0) str.RemoveAtom(0);
	for(uint i=0;i<atoms.size();i++) str.AddAtom(atoms.at(i));
	// make alphabetic
	str.SpeciesPutAlphabetic();
	// fix title
	//     aurostd::StringSubst(str.title,label_icsd,label_icsd+"."+label_postfix);
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(1): " << str.title << endl; }
	aurostd::StringSubst(str.title,label_icsd,"");aurostd::StringSubst(str.title,label_compound+"_","");aurostd::StringSubst(str.title,"["+label_compound+"]","");
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(2): " << str.title << endl; }
	aurostd::StringSubst(str.title,label_compound,"");aurostd::StringSubst(label_compound,"1","");aurostd::StringSubst(str.title,label_compound,"");
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(3): " << str.title << endl; }
	aurostd::StringSubst(str.title,"()","");aurostd::StringSubst(str.title,"(icsd library)","");aurostd::StringSubst(str.title,"  "," ");
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(4): " << str.title << endl; }
	string aus=" "+str.title;aurostd::StringSubst(str.title,"  "," ");
	str.title="";
	for(uint i=0;i<vatomX.size();i++) str.title+=vatomX.at(i);
	str.title+="/"+label_icsd+"."+label_postfix+aus;
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(5): " << str.title << endl; }
	//      *voss << str.title << endl;exit(0);
	// done
	return str;
      }
    }

    // *********************************************************************
    // HTQC
    // *voss << "MODE=" << mode << endl; exit(0);
    if(mode==LIBRARY_MODE_HTQC) { // mode =LIBRARY_MODE_HTQC
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) mode==LIBRARY_MODE_HTQC" << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) vatomX.size()=" << vatomX.size() << endl; }
      // check if binary ternary or so on..
      uint nspeciesHTQC=aflowlib::PrototypeLibrariesSpeciesNumber(label);
      if(nspeciesHTQC>vatomX.size()) {
	if(LDEBUG) { *voss << "gotta fix if: nspeciesHTQC=" << nspeciesHTQC << " > vatomX.size()=" << vatomX.size() << endl; }
	for(uint i=vatomX.size();i<nspeciesHTQC;i++) {
	  string_tmp="A";string_tmp[0]+=i;
	  vatomX.push_back(string_tmp);
	  vvolumeX.push_back(vvolumeX.at(0));
	}
      }
      // fixed now build up
      for(uint i=0;i<vatomX.size();i++) {
	if(LDEBUG) { *voss << "DEBUG Building up speciesX and natomsX. We have vatomX.at(i)=" << vatomX.at(i) << endl; }
	string_tmp="A";string_tmp[0]+=i; // do ABC... in ascii
	speciesX.push_back(string_tmp);
	// str.species.push_back(vatomX.at(i));
	natomsX.push_back(0); // start from NOTHING
      }
      _pstrstream=&XHOST_Library_HTQC; // default is HTQC
    }
    database.clear();aurostd::string2tokens(*_pstrstream,database,"\n");
    // FIND STRUCTURE kstructure
    string label_library="";
    kstructure=0;kstop=0;

    for(uint i=0;i<database.size();i++)
      if(aurostd::substring2bool(database[i],"STRUCTURE"))
	if(aurostd::substring2bool(database[i]," "+label+" ")) { // label or " "+label+" "
	  aurostd::string2tokens(database[i],tokens);
	  for(uint j=0;j<tokens.size();j++)
	    if(tokens[j]==label) {found=TRUE;kstructure=i;label_library=tokens[1];}
	}

    // LOOKING FOR TET MODIFICATIONS
    if(!found) {
      if(aurostd::substring2bool(label,"-tet")) {
	isTET=TRUE;
	if(isTET && LDEBUG) *voss << "TET modification" << endl;
	aurostd::string2tokens(label,tokens,"-tet");
	string labelAUS=tokens.at(0);
	optionsTET=label;aurostd::StringSubst(optionsTET,tokens.at(0),"");
	if(aurostd::substring2bool(label,"%")) {aurostd::string2tokens(label,tokens,"%");labelAUS+=tokens.at(1);aurostd::StringSubst(optionsTET,tokens.at(1),"");}
	aurostd::StringSubst(optionsTET,"-tet","");
	if(LDEBUG) { *voss << "LOADING PROTO=" << labelAUS << endl; }
	if(LDEBUG) { *voss << "optionsTET=" << optionsTET << endl; }
	for(uint i=0;i<database.size();i++)
	  if(aurostd::substring2bool(database[i],"STRUCTURE"))
	    if(aurostd::substring2bool(database[i]," "+labelAUS+" ")) { // label or " "+label+" "
	      aurostd::string2tokens(database[i],tokens);
	      for(uint j=0;j<tokens.size();j++)
		if(tokens[j]==labelAUS) {found=TRUE;kstructure=i;label_library=label;}   // label_library=tokens[1];
	    }
      }
    }

    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) FOUND" << endl; }
    
    // *********************************************************************
    // IF NOT FOUND
    if(!found) {
      if(mode==LIBRARY_MODE_HTQC)
	*voss << "ERROR - aflowlib::PrototypeLibraries: Structure not present in the HTQC database " << _HTQC_PROJECT_STRING_ << "  (" << label << ")" << endl;
      if(mode==LIBRARY_MODE_ICSD_AFLOWLIB || mode==LIBRARY_MODE_HTQC_ICSD_AFLOWLIB)
	*voss << "ERROR - aflowlib::PrototypeLibraries: Structure not present in the ICSD_AFLOWLIB database " << _ICSD_AFLOWLIB_STRING_ << "  (" << label << ")" << endl;
      if(mode==LIBRARY_MODE_ICSD || mode==LIBRARY_MODE_HTQC_ICSD) {
	*voss << "WARNING - aflowlib::PrototypeLibraries: Structure not present in the ICSD database " << _ICSD_STRING_ << "  (" << label << ")" << endl;
  string server=XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER"); //CO 180622 - might be empty if not specified
  if(server.empty()){server=AFLOWLIB_SERVER_DEFAULT;}
	*voss << "WARNING - aflowlib::PrototypeLibraries: Switching to " << server << " download" << endl;
	if(mode==LIBRARY_MODE_ICSD) return PrototypeLibraries(oss,label,parameters,vatomX,vvolumeX,volume_in,LIBRARY_MODE_ICSD_AFLOWLIB,flip_option);
	if(mode==LIBRARY_MODE_HTQC_ICSD) return PrototypeLibraries(oss,label,parameters,vatomX,vvolumeX,volume_in,LIBRARY_MODE_HTQC_ICSD_AFLOWLIB,flip_option);
	
	// *voss  << "XHOST_vLibrary_ICSD.at(1).length()=" << XHOST_vLibrary_ICSD.at(1).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(2).length()=" << XHOST_vLibrary_ICSD.at(2).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(3).length()=" << XHOST_vLibrary_ICSD.at(3).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(4).length()=" << XHOST_vLibrary_ICSD.at(4).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(5).length()=" << XHOST_vLibrary_ICSD.at(5).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(6).length()=" << XHOST_vLibrary_ICSD.at(6).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(7).length()=" << XHOST_vLibrary_ICSD.at(7).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(8).length()=" << XHOST_vLibrary_ICSD.at(8).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(9).length()=" << XHOST_vLibrary_ICSD.at(9).length() << endl;
	// *voss  << "XHOST_vLibrary_ICSD.at(10).length()=" << XHOST_vLibrary_ICSD.at(10).length() << endl;
      }
      exit(0);
    }

    
    // *********************************************************************
    // FIND STRUCTURE kstop
    kstop=0;
    for(uint i=kstructure+1;i<database.size() && !kstop;i++) {
      if(aurostd::substring2bool(database[i],"STRUCTURE")) {
	kstop=i;
      }
    }

    if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kstructure=" << kstructure << endl; }
    if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kstop=" << kstop << endl; }
    // FIND STRUCTURE katoms
    klattice=kstructure;katoms=kstructure; // but must be fixed !
    for(uint i=kstructure;i<kstop;i++) {
      //  *voss << database[i] << endl;
      if(aurostd::substring2bool(database[i],"PROTOTYPE")) kprototype=i;
      if(aurostd::substring2bool(database[i],"INFO")) kinfo=i;
      if(aurostd::substring2bool(database[i],"MODE")) kmode=i;
    }
    if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kprototype=" << kprototype << endl; }
    if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kinfo=" << kinfo << endl; }
    if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kmode=" << kmode << endl; }
    // FIX PROTOTYPE
    aurostd::string2tokens(database[kprototype],tokens);
    for(uint j=1;j<tokens.size();j++) str.prototype+=tokens[j]+" ";
    // FIX INFO
    aurostd::string2tokens(database[kinfo],tokens);
    for(uint j=1;j<tokens.size();j++) str.info+=tokens[j]+" ";
    str.info+=_HTQC_PROJECT_STRING_;
    // FIX MODE
    if(LDEBUG) { *voss << database[kmode] << endl; }
    aurostd::string2tokens(database[kmode],tokens," ");
    mode_prim=STRUCTURE_MODE_NONE;
    mode_conventional=STRUCTURE_MODE_NONE;
    mode_volume=STRUCTURE_MODE_NONE;
    mode_load=STRUCTURE_MODE_NONE;
    mode_remove=STRUCTURE_MODE_NONE;
    for(uint j=1;j<tokens.size();j++) {
      aurostd::string2tokens(tokens.at(j),ttokens,"=");
      if(LDEBUG) { *voss << tokens.at(j) << endl; }
      if(tokens.at(j)=="VOLUME" || tokens.at(j)=="VOL") {
	mode_volume=STRUCTURE_MODE_VOLUME;}
      if(tokens.at(j)=="PRIM") {
	mode_prim=STRUCTURE_MODE_PRIM;}
      if(tokens.at(j)=="CONVENTIONAL" || tokens.at(j)=="CONV") {
	mode_conventional=STRUCTURE_MODE_CONVENTIONAL;}
      if(tokens.at(j)=="RAW") {
	mode_load=STRUCTURE_MODE_RAW;}
      if(tokens.at(j)=="ABC") {
	mode_load=STRUCTURE_MODE_ABC;}
      if(tokens.at(j)=="WYC") {
	mode_load=STRUCTURE_MODE_WYC;}
      if(tokens.at(j)=="WYCK_ICSD" || tokens.at(j)=="WYC_ICSD" || tokens.at(j)=="ICSD") {
	mode_load=STRUCTURE_MODE_ICSD;mode_volume=STRUCTURE_MODE_VOLUME;}
      if(aurostd::substring2bool(tokens.at(j),"USE=")) {
	mode_load=STRUCTURE_MODE_USE;aurostd::string2tokens(tokens.at(j),ttokens,"=");label_use=ttokens.at(1);}
      if(tokens.at(j)=="REMOVE") {
	mode_remove=STRUCTURE_MODE_REMOVE;label_remove=tokens.at(j+1);}
      if(aurostd::substring2bool(tokens.at(j),"SPECIES=")) {
	mode_species=STRUCTURE_MODE_SPECIES;aurostd::string2tokens(tokens.at(j),ttokens,"=");label_species=ttokens.at(1);}
      if(aurostd::substring2bool(tokens.at(j),"SWAP_AB=") || aurostd::substring2bool(tokens.at(j),"SWAP_BA=")) {
	mode_load=STRUCTURE_MODE_SWAP_AB;aurostd::string2tokens(tokens.at(j),ttokens,"=");label_swap=ttokens.at(1);}
      if(aurostd::substring2bool(tokens.at(j),"SWAP(") && aurostd::substring2bool(tokens.at(j),")=")) {
	mode_load=STRUCTURE_MODE_SWAP_XY;
	aurostd::string2tokens(tokens.at(j),ttokens,"=");
	index_swap=ttokens.at(0);aurostd::StringSubst(index_swap,"SWAP(","");aurostd::StringSubst(index_swap,")","");
	label_swap=ttokens.at(1);
      }
    }

    if(mode_prim==STRUCTURE_MODE_PRIM) { if(LDEBUG) { *voss << "(mode_prim==STRUCTURE_MODE_PRIM)" << endl; } }
    if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL) { if(LDEBUG) { *voss << "(mode_conventional==STRUCTURE_MODE_CONVENTIONAL)" << endl; } }
    if(mode_load==STRUCTURE_MODE_NONE) {
      *voss << "invalid mode" << endl;
      *voss << "rerun with --debug" << endl;
      exit(0);
    }

    if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) mode=" << mode_load << endl; }

    // *voss << "HERE" << endl;exit(0);

    // FIX TITLE
    if(mode_load!=STRUCTURE_MODE_ICSD) {
      str.title="";
      for(uint i=0;i<vatomX.size();i++) str.title+=vatomX.at(i);
      str.title+="/"+label_library+" - "+"("+label+")";
      str.title+=" - "+str.prototype+" "+_HTQC_PROJECT_STRING_+" ";
    }
    if(mode_load==STRUCTURE_MODE_ICSD) {
      str.title=label_library+" - "+"("+label+")";
      str.title+=" - "+str.prototype+" "+_ICSD_STRING_;
    }

    // --------------------------------------------------------------------------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_USE && mode_species!=STRUCTURE_MODE_SPECIES) {
      // LDEBUG=TRUE;
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] (mode_load==STRUCTURE_MODE_USE)" << endl; }
      string aus_title=str.title;
      str=aflowlib::PrototypeLibraries(*voss,label_use,parameters,vatomX,vvolumeX,volume_in,mode);
      // if(LDEBUG) { {oss << str << endl;exit(0);}
      str.title=aus_title+" (use of "+label_use+")";
      if(LDEBUG) { *voss << str << endl; }
      str.FixLattices();
      natoms=str.atoms.size();
      nspecies=str.species.size();
      natomsX.clear();
      for(uint i=0;i<str.species.size();i++) natomsX.push_back(0);
      for(uint i=0;i<str.atoms.size();i++) natomsX.at(str.atoms.at(i).type)++;
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.num_each_type.size()=" << str.num_each_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.copmp_each_type.size()=" << str.comp_each_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp.size()=" << str.species_pp.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_type.size()=" << str.species_pp_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_version.size()=" << str.species_pp_version.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_ZVAL.size()=" << str.species_pp_ZVAL.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_vLDAU.size()=" << str.species_pp_vLDAU.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_volume.size()=" << str.species_volume.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_mass.size()=" << str.species_mass.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.num_each_type.size()=" << str.num_each_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] natoms=" << natoms << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] nspecies=" << nspecies << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] natomsX.size()=" << natomsX.size() << endl; }
      if(LDEBUG) {
	for(uint i=0;i<str.species.size();i++) {
	  *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species.at(" << i << ")=" << str.species.at(i)
		<< " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
		<< " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
		<< " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
		<< " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
		<< " str.species_pp_vLDAU.at(" << i << ").size()=" << str.species_pp_vLDAU.at(i).size()
		<< " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
		<< " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
		<< " str.num_each_type.at(" << i << ")=" << str.num_each_type.at(i)
		<< " str.comp_each_type.at(" << i << ")=" << str.comp_each_type.at(i)
		<< " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;
	}
      }
      // natoms=str.atoms.size();
      // nspecies=str.species.size();
      // if(isTET) { aflowlib::PrototypeFixTET(*voss,str,optionsTET); }
      // return str;
    }
    
    // -------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_SWAP_AB) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP] (mode_load==STRUCTURE_SWAP_AB)" << endl; }
      xstructure str_tmp;
      deque<string> vatomX_tmp;
      deque<double> vvolumeX_tmp;
      vatomX_tmp.push_back(vatomX.at(1));vvolumeX_tmp.push_back(vvolumeX.at(1));
      vatomX_tmp.push_back(vatomX.at(0));vvolumeX_tmp.push_back(vvolumeX.at(0));

      str_tmp=aflowlib::PrototypeLibraries(*voss,label_swap,parameters,vatomX_tmp,vvolumeX_tmp,volume_in,mode);
      // FIX THE SPECIES
      if(str_tmp.species.size() != str_tmp.num_each_type.size()) {
	*voss << "ERROR: aflow_xproto.cpp  in the SWAP_AB MODE" << endl; exit(0);
      }

      nspecies=str_tmp.num_each_type.size();
      str.scale=str_tmp.scale;
      str.lattice=str_tmp.lattice;
      str.num_each_type.clear();
      str.comp_each_type.clear();
      str.neg_scale=str_tmp.neg_scale;
      str.title+=" (swap of "+label_swap+")";
      // str.species.clear(); for(uint i=0;i<nspecies;i++) str.species.push_back(""); // create enough species
      // str.species_pp.clear(); for(uint i=0;i<nspecies;i++) str.species_pp.push_back(""); // create enough species_pp
      // str.species_pp_type.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_type.push_back(""); // create enough species_pp_type
      // str.species_pp_version.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_version.push_back(""); // create enough species_pp_version
      // str.species_pp_ZVAL.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_ZVAL.push_back(0.0); // create enough species_pp_ZVAL
      // str.species_pp_vLDAU.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_vLDAU.push_back(deque<double>()); // create enough species_pp_vLDAU
      // str.species_volume.clear(); for(uint i=0;i<nspecies;i++) str.species_volume.push_back(0.0); // create enough species_volume
      // str.species_mass.clear(); for(uint i=0;i<nspecies;i++) str.species_mass.push_back(0.0); // create enough species_mass
      
      if(nspecies>1)
	for(iat=0;iat<str_tmp.atoms.size();iat++)
	  if(str_tmp.atoms.at(iat).type==1)
	    str.AddAtom(str_tmp.atoms.at(iat));
      for(iat=0;iat<str_tmp.atoms.size();iat++)
	if(str_tmp.atoms.at(iat).type==0)
	  str.AddAtom(str_tmp.atoms.at(iat));
      if(nspecies>2)
	for(iat=0;iat<str_tmp.atoms.size();iat++)
	  if(str_tmp.atoms.at(iat).type==2)
	    str.AddAtom(str_tmp.atoms.at(iat));
      
      if(nspecies==1) {
	str.num_each_type.clear();str.num_each_type.push_back(str_tmp.num_each_type.at(0));
	str.comp_each_type.clear();str.comp_each_type.push_back(str_tmp.comp_each_type.at(0));
	str.species.clear();str.species.push_back(str_tmp.species.at(0));
	str.species_pp.clear();str.species_pp.push_back(str_tmp.species_pp.at(0));
	str.species_pp_type.clear();str.species_pp_type.push_back("");
	str.species_pp_version.clear();str.species_pp_version.push_back("");
	str.species_pp_ZVAL.clear();str.species_pp_ZVAL.push_back(0.0);
	str.species_pp_vLDAU.clear();str.species_pp_vLDAU.push_back(deque<double>());
	str.species_volume.clear();str.species_volume.push_back(str_tmp.species_volume.at(0));
	str.species_mass.clear();str.species_mass.push_back(str_tmp.species_mass.at(0));
      }
      
      if(nspecies>1) {
	str.num_each_type.clear();str.num_each_type.push_back(str_tmp.num_each_type.at(1));str.num_each_type.push_back(str_tmp.num_each_type.at(0));
	str.comp_each_type.clear();str.comp_each_type.push_back(str_tmp.comp_each_type.at(1));str.comp_each_type.push_back(str_tmp.comp_each_type.at(0));
	str.species.clear();str.species.push_back(str_tmp.species.at(1));str.species.push_back(str_tmp.species.at(0));
	str.species_pp.clear();str.species_pp.push_back(str_tmp.species_pp.at(1));str.species_pp.push_back(str_tmp.species_pp.at(0));
	str.species_pp_type.clear();str.species_pp_type.push_back("");str.species_pp_type.push_back("");
	str.species_pp_version.clear();str.species_pp_version.push_back("");str.species_pp_version.push_back("");
	str.species_pp_ZVAL.clear();str.species_pp_ZVAL.push_back(0.0);str.species_pp_ZVAL.push_back(0.0);
	str.species_pp_vLDAU.clear();str.species_pp_vLDAU.push_back(deque<double>());str.species_pp_vLDAU.push_back(deque<double>());
	str.species_volume.clear();str.species_volume.push_back(str_tmp.species_volume.at(1));str.species_volume.push_back(str_tmp.species_volume.at(0));
	str.species_mass.clear();str.species_mass.push_back(str_tmp.species_mass.at(1));str.species_mass.push_back(str_tmp.species_mass.at(0));

	if(nspecies>2) {
	  for(uint isp=2;isp<nspecies;isp++) {
	    str.num_each_type.push_back(str_tmp.num_each_type.at(isp));
	    str.comp_each_type.push_back(str_tmp.comp_each_type.at(isp));
	    str.species.push_back(str_tmp.species.at(isp));
	    str.species_pp.push_back(str_tmp.species_pp.at(isp));
	    str.species_pp_type.push_back("");
	    str.species_pp_version.push_back("");
	    str.species_pp_ZVAL.push_back(0.0);
	    str.species_pp_vLDAU.push_back(deque<double>());
	    str.species_volume.push_back(str_tmp.species_volume.at(isp));
	    str.species_mass.push_back(str_tmp.species_mass.at(isp));
	  }
	}
      }

      str.species_pp=str.species;
      // patch for species size to avoid pushing back inside
      if(str.num_each_type.size()<str.species.size()) {
	while(str.num_each_type.size()<str.species.size()) {
	  oss << "DELETING str.species.at(str.species.size()-1)=" << str.species.at(str.species.size()-1) << endl;
	  str.species.pop_back();
	  str.species_pp.pop_back();
	  str.species_pp_type.pop_back();
	  str.species_pp_version.pop_back();
	  str.species_pp_ZVAL.pop_back();
	  str.species_pp_vLDAU.pop_back();
	  str.species_volume.pop_back();
	  str.species_mass.pop_back();
	}
      }
      if(isTET) aflowlib::PrototypeFixTET(*voss,str,optionsTET);
      return str;
    }
    
    // -------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_SWAP_XY) {
      // LDEBUG=TRUE;
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_SWAP_XY)" << endl; }
      // xstructure str;
      deque<string> vatomX_tmp(vatomX);
      deque<double> vvolumeX_tmp(vvolumeX);
      aurostd::StringSubst(index_swap,"SWAP(","");aurostd::StringSubst(index_swap,")",""); // SAFETY
      aurostd::string2tokens(index_swap,tokens,",");
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 2 " << endl; }
      if(tokens.size()<2) {oss << "ERROR: xproto.cpp tokens.size()>=2 in STRUCTURE_SWAP_XY" << endl;exit(0);}
      vector<uint> ispecies(2);
      for(uint i=0;i<tokens.size();i++) {
	if(tokens.at(i)=="A" || tokens.at(i)=="a" || tokens.at(i)=="0") ispecies.at(i)=0;
	if(tokens.at(i)=="B" || tokens.at(i)=="b" || tokens.at(i)=="1") ispecies.at(i)=1;
	if(tokens.at(i)=="C" || tokens.at(i)=="c" || tokens.at(i)=="2") ispecies.at(i)=2;
	if(tokens.at(i)=="D" || tokens.at(i)=="d" || tokens.at(i)=="3") ispecies.at(i)=3;
	if(tokens.at(i)=="E" || tokens.at(i)=="e" || tokens.at(i)=="4") ispecies.at(i)=4;
	if(tokens.at(i)=="F" || tokens.at(i)=="f" || tokens.at(i)=="5") ispecies.at(i)=5;
	if(tokens.at(i)=="G" || tokens.at(i)=="g" || tokens.at(i)=="6") ispecies.at(i)=6;
	if(tokens.at(i)=="H" || tokens.at(i)=="h" || tokens.at(i)=="7") ispecies.at(i)=7;
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 3 " << endl; }
      if(LDEBUG) { *voss << "ispecies.size()=" << ispecies.size() << endl; }
      
      if(ispecies.at(0)!=ispecies.at(1)) {
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 4 " << endl; }
	// SWAPOUT
	xstructure str_tmp;
	aurostd::swap(vatomX_tmp.at(ispecies.at(0)),vatomX_tmp.at(ispecies.at(1)));
	aurostd::swap(vvolumeX_tmp.at(ispecies.at(0)),vvolumeX_tmp.at(ispecies.at(1)));
	str_tmp=aflowlib::PrototypeLibraries(*voss,label_swap,parameters,vatomX_tmp,vvolumeX_tmp,volume_in,mode);
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 5 " << endl; }
	
	// FIX THE SPECIES
	if(str_tmp.species.size() != str_tmp.num_each_type.size()) {
	  oss << "ERROR: aflow_xproto.cpp  in the SWAP_XY MODE" << endl;
	  exit(0);
	}
	
	nspecies=str_tmp.num_each_type.size();
	str.scale=str_tmp.scale;
	str.lattice=str_tmp.lattice;
	str.num_each_type.clear();
	str.comp_each_type.clear();
	str.neg_scale=str_tmp.neg_scale;
	str.title+=" (swap("+index_swap+") of "+label_swap+")";
	// str.num_each_type=str_tmp.num_each_type; // copy deque num_each_type
	// str.comp_each_type=str_tmp.comp_each_type; // copy deque comp_each_type
	// str.species=str_tmp.species; // copy deque species
	// str.species_pp=str_tmp.species_pp; // copy deque species_pp
	// str.species_pp_type=str_tmp.species_pp_type; // copy deque species_pp_type
	// str.species_pp_version=str_tmp.species_pp_version; // copy deque species_pp_version
	// str.species_pp_ZVAL=str_tmp.species_pp_ZVAL; // copy deque species_pp_ZVAL
	// str.species_pp_vLDAU=str_tmp.species_pp_vLDAU.size(); // copy deque species_pp_vLDAU
	// str.species_volume=str_tmp.species_volume; // copy deque species_volume
	// str.species_mass=str_tmp.species_mass; // copy deque species_mass
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.num_each_type=" << str_tmp.num_each_type << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.comp_each_type=" << str_tmp.comp_each_type << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species=" << str_tmp.species << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp=" << str_tmp.species_pp << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_type=" << str_tmp.species_pp_type << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_version=" << str_tmp.species_pp_version << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_ZVAL=" << str_tmp.species_pp_ZVAL << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_vLDAU=" << str_tmp.species_pp_vLDAU.size() << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_volume=" << str_tmp.species_volume << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_mass=" << str_tmp.species_mass << endl; }
	// nspecies=str.num_each_type.size();
	// *voss << nspecies << endl;
	// *voss << str_tmp.atoms.at(0).type << str_tmp.atoms.at(1).type << str_tmp.atoms.at(2).type << endl;
	// *voss << vatomX.at(ispecies.at(0)) <<  vatomX.at(ispecies.at(1)) << endl;
	// *voss << vatomX_tmp.at(ispecies.at(0)) <<  vatomX_tmp.at(ispecies.at(1)) << endl;
	if(nspecies>1) {
	  for(uint isp=0;isp<vatomX.size();isp++)
	    for(uint iat=0;iat<str_tmp.atoms.size();iat++)
	      if(str_tmp.atoms.at(iat).name==vatomX.at(isp)) {
		if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] plugging=" <<  str_tmp.atoms.at(iat).name << endl; }
		str.AddAtom(str_tmp.atoms.at(iat));
	      }
	}
	// aurostd::swap(str.num_each_type.at(ispecies.at(0)),str.num_each_type.at(ispecies.at(1)));  // taken care by AddAtom
	// aurostd::swap(str.comp_each_type.at(ispecies.at(0)),str.comp_each_type.at(ispecies.at(1)));  // taken care by AddAtom
	// aurostd::swap(str.species.at(ispecies.at(0)),str.species.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp.at(ispecies.at(0)),str.species_pp.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_type.at(ispecies.at(0)),str.species_pp_type.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_version.at(ispecies.at(0)),str.species_pp_version.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_ZVAL.at(ispecies.at(0)),str.species_pp_ZVAL.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_vLDAU.at(ispecies.at(0)),str.species_pp_vLDAU.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_volume.at(ispecies.at(0)),str.species_volume.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_mass.at(ispecies.at(0)),str.species_mass.at(ispecies.at(1))); // taken care by AddAtom

	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.num_each_type=" << str.num_each_type << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.comp_each_type=" << str.comp_each_type << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species=" << str.species << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp=" << str.species_pp << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_type=" << str.species_pp_type << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_version=" << str.species_pp_version << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_ZVAL=" << str.species_pp_ZVAL << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_vLDAU=" << str.species_pp_vLDAU.size() << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_volume=" << str.species_volume << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_mass=" << str.species_mass << endl; }

      }

      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] X " << endl; }
      // str.SpeciesPutAlphabetic();
      if(isTET) aflowlib::PrototypeFixTET(*voss,str,optionsTET);
      return str;
    }
    // -------------------------------------------------------------------
    // LATTICE
    if(mode_load==STRUCTURE_MODE_RAW) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_MODE_RAW)" << endl; }
      klattice=kmode+1;
      katoms=klattice+3;
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) klattice=" << klattice << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) katoms=" << katoms << endl; }
      // GET LATTICE
      for(uint i=1;i<=3;i++) {
	aus.clear();aus.str(aurostd::RemoveRounding(database[klattice+i-1]));
	for(uint j=1;j<=3;j++)
	  aus >> str.lattice(i,j);
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.lattice=" << endl << str.lattice << endl; }
      str.FixLattices(); // abc alpha beta gamma calculations are inside
    }
    if(mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      //   LDEBUG=TRUE;
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_MODE_WYC)" << endl; }
      klattice=kmode+1;
      katoms=klattice+1;
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) klattice=" << klattice << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) katoms=" << katoms << endl; }
      aus.clear();aus.str(aurostd::RemoveRounding(database[klattice])+" 0 "); // put the phony origin choice if asked.
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) " << aus.str() << endl; }
      aus >> str.a >> str.b >> str.c;
      aus >> str.alpha >> str.beta >> str.gamma;
      aus >> str.spacegroupnumber;
      aus >> str.spacegroupnumberoption;
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.alpha=" << str.alpha << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.beta=" << str.beta << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.gamma=" << str.gamma << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.a=" << str.a << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.b=" << str.b << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.c=" << str.c << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=" << str.spacegroupnumber << "." << str.spacegroupnumberoption << endl; }
      if(SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) {
	//     *voss << "SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE" << endl;
	if(flip_option==FALSE) {
	  // *voss << "SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE && flip_option==FALSE" << endl;
	  // *voss << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
	  // *voss << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;
	  // some ICSD hacks
	  if(str.spacegroupnumber==3 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==4 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==5 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==6 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==8 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==9 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==10 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==11 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==12 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==13 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==14 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==15 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==48 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==50 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  // if(str.spacegroupnumber==59 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK (sometimes it is 2, though)
	  if(str.spacegroupnumber==62 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=2; } // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==68 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==70 && str.spacegroupnumberoption==0) {str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==85 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==86 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==88 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==125 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==126 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==129 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=2; } // ICSD seems std OK
	  if(str.spacegroupnumber==130 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==133 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==134 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==137 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==138 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==141 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==142 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==146 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  // if(str.spacegroupnumber==148 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK
	  if(str.spacegroupnumber==148 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=2; } // ICSD seems std OK
	  if(str.spacegroupnumber==155 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; } // ICSD seems std OK
	  if(str.spacegroupnumber==160 && str.spacegroupnumberoption==0) {
	    if(isequal(str.alpha,str.beta,0.1) && isequal(str.beta,str.gamma,0.1)) { str.spacegroupnumberoption=1; } // for ICSD in RHL format try 1
	    if(isequal(str.alpha,str.beta,0.1) && isequal(str.beta,str.gamma,0.1)) { str.spacegroupnumberoption=1; } // for ICSD in RHL format try 1
	    if(isequal(str.alpha,120.0,0.1) || isequal(str.beta,120.0,0.1) || isequal(str.gamma,120.0,0.1)) { str.spacegroupnumberoption=2; } // for ICSD in HEX format try 2
	  } //912 267 9998
	  // if(str.spacegroupnumber==161 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=2; } // ICSD seems std OK
	  if(str.spacegroupnumber==166 && str.spacegroupnumberoption==0) { //str.spacegroupnumberoption=2;  // ICSD seems std OK
	    // for 166 you can switch to 2 if there is a 120 angle, otherwise you should keep 1
	    if(isequal(str.alpha,str.beta,0.1) && isequal(str.beta,str.gamma,0.1)) { str.spacegroupnumberoption=1; }// for ICSD in RHL format try 1
	    if(isequal(str.alpha,120.0,0.1) || isequal(str.beta,120.0,0.1) || isequal(str.gamma,120.0,0.1)) { str.spacegroupnumberoption=2; } // for ICSD in HEX format try 2
	  }
	  // if(str.spacegroupnumber==167 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=2; } // ICSD seems std OK
	  if(str.spacegroupnumber==201 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; } // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==203 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; } // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==222 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==224 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  // if(str.spacegroupnumber==227 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK (sometimes it is 2, though)
	  if(str.spacegroupnumber==227 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==228 && str.spacegroupnumberoption==0) { str.spacegroupnumberoption=1; }  // ICSD seems std OK CHECK
	  // *voss << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
	  // *voss << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;
	}
      }
      if(flip_option==TRUE && SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) { str.spacegroupnumberoption=3-str.spacegroupnumberoption; }
      if(SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE && str.spacegroupnumberoption!=1 && str.spacegroupnumberoption!=2) { str.spacegroupnumberoption=2; }
      //  if(SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) { *voss << "SPACEGROUP.OPTION=" << str.spacegroupnumber << "." << str.spacegroupnumberoption << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=" << str.spacegroupnumber << "." << str.spacegroupnumberoption << endl; }
      if(SpaceGroupOptionRequired(str.spacegroupnumber)==FALSE) { str.spacegroupnumberoption=0; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=" << str.spacegroupnumber << "." << str.spacegroupnumberoption << endl; }
      if(str.spacegroupnumberoption==3) { str.spacegroupnumberoption=2; }
      //   *voss << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
      // *voss << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;
      
      str.spacegrouplabel=GetSpaceGroupLabel(str.spacegroupnumber);
      str.spacegroup=GetSpaceGroupName(str.spacegroupnumber,str.directory); //DX 20180526 - add directory
      str.lattice=GetClat(str.a,str.b,str.c,str.alpha,str.beta,str.gamma);
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) " << str.a << " " <<  str.b << " " <<  str.c << " " <<  str.alpha << " " <<  str.beta << " " <<  str.gamma
			 << " " <<  str.spacegroupnumber << " " <<  str.spacegroupnumberoption << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.lattice=" << str.lattice << endl; }
      str.FixLattices();
      if(mode_load==STRUCTURE_MODE_ICSD) {
	// load species
	aurostd::string2tokens(database[kstructure],tokens);
	KBIN::VASP_SplitAlloySpecies(aurostd::RemoveSubStringFirst(tokens[1],"C-"),speciesX);
      }
    }
    //*voss << str << endl;exit(0);
    if(mode_load==STRUCTURE_MODE_ABC) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_MODE_ABC)" << endl; }
      klattice=kmode+1;
      katoms=klattice+1;
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) klattice=" << klattice << endl; }
      if(LDEBUG) { *voss <<"DEBUG: (aflowlib::PrototypeLibraries) katoms=" << katoms << endl; }
      aus.clear();aus.str(aurostd::RemoveRounding(database[klattice]));
      aurostd::string2tokens(database[klattice],tokens);
      // *voss << tokens.size() << endl;
      aus >> str.a >> str.b >> str.c;
      aus >> str.alpha >> str.beta >> str.gamma;
      // if(tokens.size()>=7) aus >> str.spacegroupnumber;
      if(tokens.size()>=8)  aus >> str.spacegroupnumberoption;
      str.lattice=GetClat(str.a,str.b,str.c,str.alpha,str.beta,str.gamma);
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.lattice=" << str.lattice << endl;
      str.FixLattices();
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupoption=[" << str.spacegroupoption << "]" << endl;}
    // -------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_RAW || mode_load==STRUCTURE_MODE_ABC || mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      // LDEBUG=TRUE;   // FIX ATOMS
      uint nspecies_old=nspecies;
      natoms=0;nspecies=0;
      // for(uint i=0;i<speciesX.size();i++) *voss << speciesX.at(i) << endl;
      // for(uint i=0;i<vvolumeX.size();i++) *voss << vvolumeX.at(i) << endl;
	
      // check species and volumes from library
      vector<string> lvatomX;
      for(iat=0;iat<(uint) kstop-katoms;iat++) {
	aurostd::string2tokens(database[iat+katoms],tokens);
	bool found=FALSE;
	for(uint i=0;i<lvatomX.size()&&!found;i++) if(tokens[3]==lvatomX.at(i)) found=TRUE;
	if(!found) lvatomX.push_back(tokens[3]);
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.size()=" << lvatomX.size() << endl; }
      if(LDEBUG) { for(uint i=0;i<lvatomX.size();i++) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.at(" << i << ")=" << lvatomX.at(i) << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] speciesX.size()=" << speciesX.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.size()=" << str.species.size() << endl; }
      if(lvatomX.size()>speciesX.size()) {
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] FIX nspecies" << endl; }
	for(uint i=speciesX.size();i<lvatomX.size();i++) {
	  vatomX.push_back(lvatomX.at(i));
	  speciesX.push_back(lvatomX.at(i));
	  vvolumeX.push_back(1.0);
	  natomsX.push_back(0);
	}
	nspecies_old=vatomX.size();
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.size()=" << lvatomX.size() << endl; }
      if(LDEBUG) { for(uint i=0;i<lvatomX.size();i++) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.at(" << i << ")=" << lvatomX.at(i) << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] speciesX.size()=" << speciesX.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.size()=" << str.species.size() << endl; }

      // fixed species
      for(iat=0;iat<(uint) kstop-katoms;iat++) {
	aurostd::string2tokens(database[iat+katoms],tokens);
	// *voss << iat << " " << tokens[3] << endl;
	natoms++;
	for(uint i=0;i<speciesX.size();i++)
	  if(tokens[3]==speciesX.at(i)) natomsX.at(i)++;
      }
      for(uint i=0;i<speciesX.size();i++)
	if(natomsX.at(i)>0) {str.num_each_type.push_back(natomsX.at(i));str.comp_each_type.push_back((double) natomsX.at(i));nspecies++;}
      if(mode_load==STRUCTURE_MODE_ICSD) {
	if(nspecies!=nspecies_old) {
	  oss << "ERROR (aflow_xproto.cpp): label=" << label;
	  oss << "  nspecies!=nspecies_old (" << nspecies << "," << nspecies_old << ") *************" << endl;
	  exit(0);
	}
      }
      // now I have nspecies
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] NSPECIES_OLD=" << nspecies_old << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] NSPECIES=" << nspecies << endl; }

      str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
      for(uint i=0;i<nspecies;i++) {
	str.species.push_back(vatomX.at(i));    // only the required ones
	str.species_pp.push_back(vatomX.at(i)); // only the required ones
	str.species_pp_type.push_back(""); // only the required ones
	str.species_pp_version.push_back(""); // only the required ones
	str.species_pp_ZVAL.push_back(0.0); // only the required ones
	str.species_pp_vLDAU.push_back(deque<double>()); // only the required ones
	str.species_volume.push_back(0.0);     // only the required ones
	str.species_mass.push_back(0.0);     // only the required ones
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] natoms=" << natoms << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] nspecies=" << nspecies << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) {
	for(uint i=0;i<nspecies;i++) {
	  *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.at(" << i << ")=" << str.species.at(i)
		<< " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
		<< " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
		<< " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
		<< " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
		<< " str.species_pp_vLDAU.at(" << i << ")=" << str.species_pp_vLDAU.at(i).size()
		<< " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
		<< " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
		<< " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;
	}
      }
      // Plug ATOMS
      volume=0.0;
      // perform atoms A,B,C,D,E... ---------------------------------

      for(uint ispecies=1;ispecies<=nspecies;ispecies++) { // pick atom X
	if(speciesX.size()>=ispecies) {
	  for(iat=0;iat<natoms;iat++) {
	    aurostd::string2tokens(database[iat+katoms],tokens);
	    if(tokens[3]==speciesX.at(ispecies-1)) { // atoms X
	      _atom atom;
	      atom.CleanName();atom.CleanSpin();
	      aus.clear();aus.str(database[iat+katoms]);
	      aus >> atom.fpos(1) >> atom.fpos(2) >> atom.fpos(3);
	      atom.cpos=F2C(str.lattice,atom.fpos);
	      atom.name=str.species.at(ispecies-1);atom.type=0;
	      atom.CleanName(); // to fix the real name
	      for(uint icheck=0;icheck<ispecies-1;icheck++)
		if(str.species.at(icheck)!="") atom.type++; // X is 0,1,2,3 to ispecies-1
	      atom.name_is_given=TRUE;str.atoms.push_back(atom);
	    }
	  }
	}
      }
      // ------ --------------------------------------------
      for(uint i=0;i<nspecies;i++) {
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] vvolumeX.at(" << i << ")=" << vvolumeX.at(i) << endl; }
      }
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupoption=[" << str.spacegroupoption << "]" << endl;}
    // -------------------------------------------------------------------
    if(mode_remove==STRUCTURE_MODE_REMOVE) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] BEFORE REMOVE ********************* " << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natoms=" << natoms << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] nspecies=" << nspecies << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp.size()=" << str.species_pp.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_type.size()=" << str.species_pp_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_version.size()=" << str.species_pp_version.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_ZVAL.size()=" << str.species_pp_ZVAL.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_vLDAU.size()=" << str.species_pp_vLDAU.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natomsX.size()=" << natomsX.size() << endl; }
      if(LDEBUG) {
	for(uint i=0;i<nspecies;i++) {
	  *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.at(" << i << ")=" << str.species.at(i)
		<< " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
		<< " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
		<< " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
		<< " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
		<< " str.species_pp_vLDAU.at(" << i << ")=" << str.species_pp_vLDAU.at(i).size()
		<< " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
		<< " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
		<< " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;
	}
      }
      tokens.clear();
      aurostd::string2tokens(label_remove,tokens,",");
      vector<uint> ratoms;
      for(iat=0;iat<tokens.size();iat++)
	ratoms.push_back(aurostd::string2utype<uint>(tokens.at(iat)));
      std::sort(ratoms.begin(), ratoms.end());
      for(int iat=ratoms.size()-1;iat>=0;iat--)
	str.RemoveAtom(ratoms.at(iat));
      // cout << ratoms.at(iat) << endl;
      // for(uint j=1;j<tokens.size();j++) str.prototype+=tokens[j]+" ";

      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] AFTER REMOVE ********************* " << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natoms=" << natoms << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] nspecies=" << nspecies << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp.size()=" << str.species_pp.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_type.size()=" << str.species_pp_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_version.size()=" << str.species_pp_version.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_ZVAL.size()=" << str.species_pp_ZVAL.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_vLDAU.size()=" << str.species_pp_vLDAU.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natomsX.size()=" << natomsX.size() << endl; }
      if(LDEBUG) {
	for(uint i=0;i<nspecies;i++) {
	  *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.at(" << i << ")=" << str.species.at(i)
		<< " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
		<< " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
		<< " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
		<< " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
		<< " str.species_pp_vLDAU.at(" << i << ")=" << str.species_pp_vLDAU.at(i).size()
		<< " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
		<< " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
		<< " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;
	  // exit(0);
	}
      }
    }

    // --------------------------------------------------------------------------------------------------------------------------------------
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupoption=[" << str.spacegroupoption << "]" << endl;}

    if(mode_load==STRUCTURE_MODE_USE && mode_species==STRUCTURE_MODE_SPECIES) {
      //   LDEBUG=TRUE;

      // this is USE_SPECIES
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES] (mode_load==STRUCTURE_MODE_USE) (mode_species==STRUCTURE_MODE_SPECIES) " << endl;}
      string aus_title=str.title;
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vatomX.size()=" << vatomX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vatomX="; for(uint i=0;i<vatomX.size();i++) *voss << vatomX.at(i) << " ";*voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vvolumeX.size()=" << vvolumeX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vvolumeX="; for(uint i=0;i<vvolumeX.size();i++) *voss << vvolumeX.at(i) << " ";*voss  << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  label_use=" << label_use << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  volume_in=" << volume_in << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  mode=" << mode << endl;}
      str=aflowlib::PrototypeLibraries(*voss,label_use,parameters,vatomX,vvolumeX,volume_in,2);//_mode); // STEFANO HACK
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vatomX.size()=" << vatomX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vatomX="; for(uint i=0;i<vatomX.size();i++) *voss << vatomX.at(i) << " ";*voss  << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vvolumeX.size()=" << vvolumeX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vvolumeX="; for(uint i=0;i<vvolumeX.size();i++) *voss << vvolumeX.at(i) << " ";*voss  << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  label_use=" << label_use << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  volume_in=" << volume_in << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  mode=" << mode << endl;}
      //    exit(0);
      //  if(LDEBUG) {oss << str << endl;exit(0);}
      str.title=aus_title+" (use of "+label_use+")";
      if(LDEBUG) *voss << str << endl;
      str.FixLattices();
      natoms=str.atoms.size();
      nspecies=str.species.size();
      natomsX.clear();
      for(uint i=0;i<str.species.size();i++) natomsX.push_back(0);
      for(uint i=0;i<str.atoms.size();i++) natomsX.at(str.atoms.at(i).type)++;

      // this is SPECIES

      // *voss << str.scale << endl; exit(0);

      tokens.clear();
      aurostd::string2tokens(label_species,tokens,",");
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  label_species=" << label_species << endl;}
      if(tokens.size()!=str.num_each_type.size()) {oss << "tokens.size()!=str.num_each_type.size()" << endl;exit(0);}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  tokens.size()=" << tokens.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.num_each_type.size()=" << str.num_each_type.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.comp_each_type.size()=" << str.comp_each_type.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  speciesX.size()=" << speciesX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species.size()=" << str.species.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_volume.size()=" << str.species_volume.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_mass.size()=" << str.species_mass.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vvolumeX.size()=" << vvolumeX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vatomX.size()=" << vatomX.size() << endl;}

      if(LDEBUG) {for(uint i=0;i<vvolumeX.size();i++) *voss << vvolumeX.at(i) << " "; *voss << endl;}
      deque<double> vol(vvolumeX);
      for(uint i=0;i<str.num_each_type.size();i++) {
	uint j=(int) tokens.at(i)[0]-(int)'A';
	tokens.at(i)=str.species.at(j);
	vol.at(i)=vvolumeX.at(j);
      }
      vvolumeX=vol;
      for(uint i=0;i<vvolumeX.size();i++) str.species_volume.at(i)=vol.at(i);
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  volume="; for(uint i=0;i<vvolumeX.size();i++) *voss << vvolumeX.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  tokens="; for(uint i=0;i<tokens.size();i++) *voss << tokens.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species="; for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_volume="; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_mass="; for(uint i=0;i<str.species_mass.size();i++) *voss << str.species_mass.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vatomX="; for(uint i=0;i<vatomX.size();i++) *voss << vatomX.at(i) << " "; *voss << endl;}
      uint iat=0;
      for(uint i=0;i<str.num_each_type.size();i++) {
	// *voss << str.species.at(i) << endl;
	str.species.at(i)=tokens.at(i);
	for(uint j=0;j<(uint) str.num_each_type.at(i);j++) {
	  str.atoms.at(iat).type=i;
	  str.atoms.at(iat).name=str.species.at(i);
	  str.atoms.at(iat).CleanName(); // cool little name
	  iat++;
	}
      }
      // for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;
      // JUNKAI
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      str.SpeciesPutAlphabetic();
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      //   exit(0);
      // for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;
      for(uint i=0;i<str.species_volume.size();i++) vvolumeX.at(i)=str.species_volume.at(i);
      // make species together
      // if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
      for(uint i=1;i<str.num_each_type.size();i++) {
	if(str.species.at(i)==str.species.at(i-1)) {
	  str.species.erase(str.species.begin()+i);
	  str.num_each_type.at(i-1)+=str.num_each_type.at(i);str.num_each_type.erase(str.num_each_type.begin()+i);
	  str.comp_each_type.at(i-1)+=str.comp_each_type.at(i);str.comp_each_type.erase(str.comp_each_type.begin()+i);
	  vvolumeX.erase(vvolumeX.begin()+i);
	  str.species_volume.erase(str.species_volume.begin()+i);
	  str.species_mass.erase(str.species_mass.begin()+i);
	  // vatomX.erase(vatomX.begin()+i);
	  tokens.erase(tokens.begin()+i);
	  i=0;
	}
      }
      // fix types
      uint k=0;
      for(uint i=0;i<str.num_each_type.size();i++)
	for(uint j=0;j<(uint) str.num_each_type.at(i);j++)
	  str.atoms.at(k++).type=i; // done
      nspecies=str.num_each_type.size();
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  [4] nspecies=" << nspecies << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  species.size()=" << str.species.size() << endl; }
      if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  tokens.size()=" << tokens.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.num_each_type.size()=" << str.num_each_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.comp_each_type.size()=" << str.comp_each_type.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  speciesX.size()=" << speciesX.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_volume.size()=" << str.species_volume.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_mass.size()=" << str.species_mass.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vvolumeX.size()=" << vvolumeX.size() << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  vatomX.size()=" << vatomX.size() << endl; }
      if(LDEBUG) {for(uint i=0;i<vvolumeX.size();i++) *voss << vvolumeX.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<tokens.size();i++) *voss << tokens.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.species_mass.size();i++) *voss << str.species_mass.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<vatomX.size();i++) *voss << vatomX.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
      if(LDEBUG) { *voss << str << endl; }
      //   exit(0);
      //   str.SpeciesPutAlphabetic();
    }
    // -------------------------------------------------------------------
    // FIX ALL ATOMS
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) FIX ALL ATOMS (if necessary)" << endl; }
    bool DEBUG_ICSD=TRUE;
    xstructure strw=str;
    if(mode_load==STRUCTURE_MODE_ICSD) DEBUG_ICSD=TRUE;

    if(mode_load==STRUCTURE_MODE_RAW) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] RAW str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] label = " << label << " calculating WyckoffPOSITIONS()" << endl; }
      
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] str.atoms.size()=[" << str.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] vatomX.size()=[" << vatomX.size()  << "]" << endl;}
      if(LDEBUG) {for(uint i=0;i<vatomX.size();i++) *voss << vatomX.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] vvolumeX.size()=[" << vvolumeX.size()  << "]" << endl;}
      if(LDEBUG) {for(uint i=0;i<vvolumeX.size();i++) *voss << vvolumeX.at(i) << " "; *voss << endl;}
    }
    
    if(mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] WYC/ICSD str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] WYC/ICSD strw.species.size()=" << strw.species.size() << endl; }
      if(LDEBUG) {for(uint i=0;i<strw.species.size();i++) *voss << strw.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] label = " << label << " calculating WyckoffPOSITIONS()" << endl; }

      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] str.spacegroupnumberoption=[" << str.spacegroupnumberoption << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] strw.spacegroupnumber=[" << strw.spacegroupnumber << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] strw.spacegroupnumberoption=[" << strw.spacegroupnumberoption << "]" << endl;}

      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] BEFORE WyckoffPOSITIONS" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] str.atoms.size()=[" << str.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] strw.atoms.size()=[" << strw.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] vatomX.size()=[" << vatomX.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] vvolumeX.size()=[" << vvolumeX.size()  << "]" << endl;}
      str=WyckoffPOSITIONS(strw.spacegroupnumber,strw.spacegroupnumberoption,strw);
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] AFTER WyckoffPOSITIONS" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] str.atoms.size()=[" << str.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] strw.atoms.size()=[" << strw.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] vatomX.size()=[" << vatomX.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] vvolumeX.size()=[" << vvolumeX.size()  << "]" << endl;}

      //   *voss << str << endl;exit(0);
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_ALL_ATOMS] WYC/ICSD str.species.size()=" << str.species.size() << endl; }
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      
      string sgstring=(SpaceGroupOptionRequired(str.spacegroupnumber)?aurostd::utype2string(str.spacegroupnumber)+"."+aurostd::utype2string(str.spacegroupnumberoption):aurostd::utype2string(str.spacegroupnumber));
      // if(mode_load==STRUCTURE_MODE_WYC) {
      // 	label_library="";
      // 	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) vatomX.size()=" << vatomX.size() << endl; }
      // 	for(uint i=0;i<str.species.size();i++) label_library+=vatomX.at(i);
      // 	label_library+="/";	//	vatomX.at(0)+vatomX.at(1)+"/"+label_library;
      // }
      
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(6): " << str.title << endl; }

      if(mode_load==STRUCTURE_MODE_ICSD) {
	str.title=label_library+" "+"#"+sgstring+" - "+"("+label+")";  // SG# AS Wahyu wants
	str.title+=" - "+str.prototype+" "+_ICSD_STRING_;                                                // SG# as Wahyu wants
	str.title=str.title+" (WICKOFF "+sgstring+" "+str.spacegrouplabel+")";
      }
      if(mode_load==STRUCTURE_MODE_WYC) {
	str.title=str.title+" (WICKOFF "+sgstring+" "+str.spacegrouplabel+")";
      }
      if(mode_load==STRUCTURE_MODE_ICSD) {
	for(uint i=0;i<str.atoms.size();i++) {
	  str.atoms.at(i).name=speciesX.at(str.atoms.at(i).type);
	  str.atoms.at(i).CleanName();
	}
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(7): " << str.title << endl; }
      //  exit(0);

      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) WYC/ICSD str.species.size()=" << str.species.size() << endl; }
      if(0 && DEBUG_ICSD) {
	*voss << endl << "calculating space group" << endl;
	int spacegroup_test=str.spacegroupnumber;
	str.platon2sg();
	*voss << "  sgs(" << spacegroup_test << "," << str.spacegroupnumber << ")" << endl;
	if(spacegroup_test!=str.spacegroupnumber) {
	  oss << "ERROR (aflow_xproto.cpp): label=" << label;
	  oss << " " << spacegroup_test << " " << str.spacegroupnumber << " CHANGE space group option in compound " << label << "   *************" << endl;
	  exit(0);
	}
      }
      if(mode==LIBRARY_MODE_ICSD) {
	if(NearestNeighbour(str)<_XPROTO_TOO_CLOSE_ERROR_ && flip_option==FALSE && SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) {
	  // *voss << "AFLOW WARNING (aflow_xproto.cpp): label=" << label << " WRONG NNdist too close =" << NearestNeighbour(str) << "   *************"  << endl;
	  *voss << "AFLOW WARNING (aflow_xproto.cpp): Too close NNdist(" << NearestNeighbour(str) << "), Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption << " not enough, try flip_option " << endl;
	  xstructure str_flipped;
	  str_flipped=aflowlib::PrototypeLibraries(*voss,label,parameters,vatomX,vvolumeX,volume_in,mode,TRUE);
	  str_flipped.title+=" (sg.opt flipped)";
	  str_flipped.species_pp=str_flipped.species;
	  *voss << "AFLOW WARNING (aflow_xproto.cpp): Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption <<" flipped to option=" << str_flipped.spacegroupnumberoption << endl;
	  if(isTET) aflowlib::PrototypeFixTET(*voss,str_flipped,optionsTET);
	  return str_flipped;
	}
      }
      if(mode==LIBRARY_MODE_ICSD) {
	if(DEBUG_ICSD) {
	  double natoms_label=0.0,natoms_icsd=0.0;
	  for(uint i=0;i<str.num_each_type.size();i++) {
	    natoms_icsd+=str.num_each_type.at(i);
	    natoms_label+=vnatomsX.at(i);
	  }
	  for(uint i=0;i<str.num_each_type.size();i++) {
	    //   *voss << str.num_each_type.at(i) << " " << natoms_icsd << " " << vnatomsX.at(i) << " " << natoms_icsd <<  endl;
	    if(abs((double) (str.num_each_type.at(i)/natoms_icsd-vnatomsX.at(i)/natoms_label))>0.01) {
	      if(flip_option==FALSE && SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) {
		*voss << "AFLOW WARNING (aflow_xproto.cpp): Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption << " not enough, try flip_option " << endl;
		xstructure str_flipped;
		str_flipped=aflowlib::PrototypeLibraries(*voss,label,parameters,vatomX,vvolumeX,volume_in,mode,TRUE);
		str_flipped.title+=" (sg.opt flipped)";
		str_flipped.species_pp=str_flipped.species;
		*voss << "AFLOW WARNING (aflow_xproto.cpp): Spacegroup=" << str.spacegroupnumber << " found option=" << str_flipped.spacegroupnumberoption << endl;
		if(isTET) aflowlib::PrototypeFixTET(*voss,str_flipped,optionsTET);
		return str_flipped;
	      }
	      // *voss << endl;
	      *voss << "ERROR (aflow_xproto.cpp): label=" << label
		    << " WRONG concentrations species=" << str.species.at(i)
		    << " conc_prototype=" << vnatomsX.at(i)/natoms_label << " "
		    << " conc_icsd=" << str.num_each_type.at(i)/natoms_icsd << " : CHECK ALSO THE SPACEGROUP" << "   *************"  << endl;
	      if(LDEBUG) { *voss << str << endl; }
	      exit(0);
	    }
	  }
	}
      }
    }
    // -------------------------------------------------------------------
    // FIX SCALE
    // *voss << str << endl;
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] FIX SCALE (if necessary)" << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] nspecies=" << nspecies << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.atoms.size()=" << str.atoms.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] species.size()=" << str.species.size() << endl; }
    if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] tokens.size()=" << tokens.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.num_each_type.size()=" << str.num_each_type.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.comp_each_type.size()=" << str.comp_each_type.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] speciesX.size()=" << speciesX.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.species.size()=" << str.species.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.species_volume.size()=" << str.species_volume.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.species_mass.size()=" << str.species_mass.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] vvolumeX.size()=" << vvolumeX.size() << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] vatomX.size()=" << vatomX.size() << endl; }

    str.species_volume.clear();
    for(uint i=0;i<str.species.size();i++)
      str.species_volume.push_back(0.0);

    if(mode_volume==STRUCTURE_MODE_NONE) {
      if(volume_in>0.0) {
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume choice 1 (volume_in=" << volume_in << ")" << endl; }
	volume=str.atoms.size()*volume_in;
      } else {
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume choice 2: vvolumeX.size()=" << vvolumeX.size() << endl; }
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume choice 2: str.species_volume.size()="<<str.species_volume.size()<< ": "; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}
	volume=0.0;
	for(iat=0;iat<str.atoms.size();iat++) {
	  if(LDEBUG) { *voss << iat << " " << str.atoms.at(iat).name << " " << nspecies << " " << str.atoms.size() << endl; }
	  for(uint i=0;i<nspecies;i++)
	    if(str.atoms.at(iat).name==str.species.at(i)) {
	      volume+=vvolumeX.at(i);
	      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) volume atom (" << iat << ")=" << volume << endl; }
	    }

	  for(uint i=0;i<str.species.size();i++)
	    str.species_volume.at(i)=vvolumeX.at(i);
	}
      }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume=" << volume << endl; }
      // fix volume in WYC or ABC
      str.scale=std::pow((double) (abs(volume)/det(str.lattice)),(double) 1.0/3.0);
      str.neg_scale=TRUE;
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] structure_mode_volume" << endl; }
    if(mode_volume==STRUCTURE_MODE_VOLUME) {
      str.scale=1.0;
      str.ReScale(1.0);
      str.neg_scale=FALSE;
    }

    // check for nndist

    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] nndist" << endl; }
    if(NearestNeighbour(str)<_XPROTO_TOO_CLOSE_ERROR_) {
      *voss << "ERROR (aflow_xproto.cpp): label=" << label << " WRONG NNdist too close =" << NearestNeighbour(str) << " Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption <<"   *************"  << endl;
      *voss << str << endl;
      //   str.SetCoordinates(_COORDS_CARTESIAN_);
      //  *voss << str << endl;

      exit(0);
    }

    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1" << endl; }

    double eps=-1.0,epsang=-1.0; // for lattice calculation. initialize negative
    if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL || mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      //    LDEBUG=TRUE;
      /* OLD
      // get btravais_lattice_type before it is too hard
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) calling LATTICE::SpaceGroup2Lattice from xproto.cpp" << endl;
      // str.bravais_lattice_type=LATTICE::SpaceGroup2Lattice(str.spacegroupnumber);
      // str.bravais_lattice_variation_type=LATTICE::SpaceGroup2LatticeVariation(str.spacegroupnumber,str);
      // str.bravais_lattice_system=LATTICE::Lattice_System_SpaceGroup(str.spacegroupnumber);
      */
      // NEW
      //    *voss << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
      //    *voss << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;

      xstructure str_in(str),str_sp,str_sc;
      // DX - START
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1 - Symmetry tolerance scan (DX and COREY) " << endl; }
      LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
      /* DX
      //    LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.1 - precise" << endl; }
      LATTICE::Standard_Lattice_StructurePrecise(str_in,str_sp,str_sc);
      // LATTICE::Standard_Lattice_StructureMedium(str_in,str_sp,str_sc);
      if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.2 - try again - medium+" << endl; }
      int ss=0; eps=0.005;epsang=0.05; LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang,ss,_EPS_);} // medium++
      if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.3 - try again - medium" << endl; }
      int ss=0; eps=0.01;epsang=0.1; LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang,ss,_EPS_);} // medium++
      if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.4 - try again - medium-default" << endl; }
      LATTICE::Standard_Lattice_StructureMedium(str_in,str_sp,str_sc);}
      if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.5 - try again - loose" << endl; }
      LATTICE::Standard_Lattice_StructureDefault(str,str_sp,str_sc);}
      */ 
      // DX - END

      str.bravais_lattice_type=str_sp.bravais_lattice_type;
      str.bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
      str.bravais_lattice_system=str_sp.bravais_lattice_system;
      if(LDEBUG) { *voss << str.bravais_lattice_type << endl; }
      if(LDEBUG) { *voss << str_sp.bravais_lattice_type << endl; }
      if(LDEBUG) { *voss << str_sc.bravais_lattice_type << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.99" << endl; }
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 2" << endl; }
    if(mode_load==STRUCTURE_MODE_ABC) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 2a" << endl; }
      str.GetLatticeType();
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] str.bravais_lattice_type=" << str.bravais_lattice_type << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] str.bravais_lattice_variation_type=" << str.bravais_lattice_variation_type << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 2b" << endl; }
      if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL) { str.GetStandardPrimitive(); }
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 3" << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] str.species_volume.size()="<<str.species_volume.size()<< ": "; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}
      
    // now if you do the PRIM you get screwed because it is very hard to reconstruct the lattice type
    if(mode_prim==STRUCTURE_MODE_PRIM) {
      // LDEBUG=TRUE;
      // *voss << "mode_prim=" << mode_prim << endl;
      // str.neg_scale=FALSE;
      // str.ReScale(1.0);
      // str.SetVolume(str.atoms.size());
      bool fixed=FALSE;
      if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL || mode_load==STRUCTURE_MODE_ICSD) {
	//   *voss << "mode_conventional=" << mode_conventional << endl;
	fixed=TRUE;
	xstructure str_sp,str_sc;
	str.FixLattices();
	str_sp=str;str_sc=str; // so it creates species and stuff.
	if(eps<0.0 && epsang <0.0) {  // still undefined
	  // eps=0.1;epsang=1.0; // coarse
	  // eps=0.05;epsang=0.5;// normal
	  // eps=0.01;epsang=0.1; // medium//
	  eps=0.005;epsang=0.05; // medium++
	  // eps=0.002;epsang=0.02; // precise
	  if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): redefining  eps=" << eps << "  epsang=" << epsang << " " << endl; }
	} else {
	  if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): inheriting  eps=" << eps << "  epsang=" << epsang << " " << endl; }
	}
	str_sp.bravais_lattice_type="X";
        // DX - START
	// [[JUNKAI OBSOLETE]]  uint step=0;
	//  while (str_sp.bravais_lattice_type!=str.bravais_lattice_type && str_sp.bravais_lattice_type!="UNKNOWN" && step<10) {
        LATTICE::Standard_Lattice_StructureDefault(str,str_sp,str_sc); // DX
        // [[JUNKAI OBSOLETE]]
	// [[JUNKAI OBSOLETE]]while (str_sp.bravais_lattice_type!=str.bravais_lattice_type && str_sp.bravais_lattice_type!="UNKNOWN" && step<10) {
	// [[JUNKAI OBSOLETE]]  step++;
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): step=" << step << "  eps=" << eps << "  epsang=" << epsang << " "; }
	// [[JUNKAI OBSOLETE]]  // LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang); // STEFANO OLD VERSION
	// [[JUNKAI OBSOLETE]]  int ss=0;
	// [[JUNKAI OBSOLETE]]  LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang,ss,_EPS_);
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) { *voss << str.bravais_lattice_type << endl; }
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) { *voss << str_sp.bravais_lattice_type << endl; }
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) { *voss << str_sc.bravais_lattice_type << endl; }
	// [[JUNKAI OBSOLETE]]  // *voss << str << endl; *voss << str_sp << endl; *voss << str_sc << endl;exit(0);
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) { *voss << str_sp.bravais_lattice_variation_type << endl; }
	// [[JUNKAI OBSOLETE]]  // *voss << str_sp.bravais_lattice_type << " " << eps << " " << epsang << " " << endl;
	// [[JUNKAI OBSOLETE]]  // if(str_sp.bravais_lattice_type!=str.bravais_lattice_type) {  // try different option str.bravais_lattice_type known only for ICSD
	// [[JUNKAI OBSOLETE]]  if(str_sp.bravais_lattice_type!=str.bravais_lattice_type) {
	// [[JUNKAI OBSOLETE]]    xstructure stropt=strw,stropt_sp=str,stropt_sc=str;
	// [[JUNKAI OBSOLETE]]    stropt.spacegroupnumberoption=3-stropt.spacegroupnumberoption;
	// [[JUNKAI OBSOLETE]]    if(stropt.spacegroupnumberoption==3) stropt.spacegroupnumberoption=2;
	// [[JUNKAI OBSOLETE]]    stropt=WyckoffPOSITIONS(stropt.spacegroupnumber,stropt.spacegroupnumberoption,stropt);
	// [[JUNKAI OBSOLETE]]    stropt.title=label_library+" "+"#"+aurostd::utype2string(stropt.spacegroupnumber)+" - "+"("+label+")";  // SG# AS Wahyu wants
	// [[JUNKAI OBSOLETE]]    stropt.title+=" - "+stropt.prototype+" "+_ICSD_STRING_;                                                // SG# as Wahyu wants
	// [[JUNKAI OBSOLETE]]    stropt.title=stropt.title+" (WICKOFF "+stropt.spacegroup+" "+stropt.spacegrouplabel+")";
	// [[JUNKAI OBSOLETE]]    for(uint i=0;i<stropt.atoms.size();i++) {stropt.atoms.at(i).name=speciesX.at(stropt.atoms.at(i).type);stropt.atoms.at(i).CleanName();}
	// [[JUNKAI OBSOLETE]]    // stropt=LATTICE::Conventional_Lattice_Structure(stropt); // put it in the right order...
	// [[JUNKAI OBSOLETE]]    stropt.bravais_lattice_type=LATTICE::SpaceGroup2Lattice(stropt.spacegroupnumber);
	// [[JUNKAI OBSOLETE]]    stropt.bravais_lattice_variation_type=LATTICE::SpaceGroup2LatticeVariation(stropt.spacegroupnumber,stropt);
	// [[JUNKAI OBSOLETE]]    stropt.bravais_lattice_system=stropt.bravais_lattice_type; // FIX
	// [[JUNKAI OBSOLETE]]    stropt.FixLattices();
	// [[JUNKAI OBSOLETE]]    stropt_sp=stropt;stropt_sc=stropt;stropt_sp.bravais_lattice_type="X";  // so it creates species and stuff.
	// [[JUNKAI OBSOLETE]]    step++;
	// [[JUNKAI OBSOLETE]]    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): step=" << step << "  eps=" << eps << "  epsang=" << epsang << " "; }
	// [[JUNKAI OBSOLETE]]    // LATTICE::Standard_Lattice_Structure(stropt,stropt_sp,stropt_sc,eps,epsang); // STEFANO OLD VERSION
	// [[JUNKAI OBSOLETE]]    ss=0; // JUNKAI
	// [[JUNKAI OBSOLETE]]    LATTICE::Standard_Lattice_Structure(stropt,stropt_sp,stropt_sc,eps,epsang,ss,_EPS_); //JUNKAI
	// [[JUNKAI OBSOLETE]]    if(LDEBUG) { *voss << str_sp.bravais_lattice_variation_type << endl; }
	// [[JUNKAI OBSOLETE]]    if(stropt_sp.bravais_lattice_type==stropt.bravais_lattice_type) { // gotta a wrong OPTION
	// [[JUNKAI OBSOLETE]]      *voss << "WARNING: (aflow_xproto.cpp): label=" << label << "  wrong option in lattice, put option=" << stropt.spacegroupnumberoption << endl;
	// [[JUNKAI OBSOLETE]]      str=stropt;str_sp=stropt_sp;str_sc=stropt_sc;
	// [[JUNKAI OBSOLETE]]    }
	// [[JUNKAI OBSOLETE]]  }
	// [[JUNKAI OBSOLETE]]  eps=eps/3.0;epsang=epsang/3.0;
	// [[JUNKAI OBSOLETE]]}
        // DX - END

	if(str_sp.bravais_lattice_type!=str.bravais_lattice_type) {
	  oss << "ERROR (aflow_xproto.cpp): label=" << label << "  " << str_sp.bravais_lattice_type << " " << str.bravais_lattice_type << " "
	      << "FOUND Lattice differs from the ICSD one." << endl;
	  str_sp.error_flag=TRUE;         // flag TRUE is error
	  str_sp.error_string="ERROR: Found lattice="+str_sp.bravais_lattice_type+"  SGlattice="+str.bravais_lattice_type+"  ";
	}
	// *voss << str_sp.bravais_lattice_type << endl;
	// exit(0);
	str_sp.FixLattices();str_sp.BringInCell();
	str_sc.FixLattices();str_sc.BringInCell();
	str=str_sp; // standard primitive
      } else {
	//     *voss << "PRIM not CONV" << endl;
	uint natoms=str.atoms.size();
	str.GetPrimitive(0.005);
	if(natoms<str.atoms.size()) { fixed=TRUE; }
      }
      if(!fixed) { str.GetPrimitive(0.05); } // now you are totally screwed..... but you save plenty of ab-initio time. Wisdom of stefano.
    }
    // done
    str.BringInCell();
    str.FixLattices();
    str.MakeBasis();
    str.species_pp=str.species;
    // done
    if(LDEBUG) {
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.num_each_type.size()="<<str.num_each_type.size()<<": "; for(uint i=0;i<str.num_each_type.size();i++) *voss << str.num_each_type.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.comp_each_type.size()="<<str.comp_each_type.size()<<": "; for(uint i=0;i<str.comp_each_type.size();i++) *voss << str.comp_each_type.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species.size()="<<str.species.size()<< ": "; for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp.size()="<<str.species_pp.size()<< ": "; for(uint i=0;i<str.species_pp.size();i++) *voss << str.species_pp.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_type.size()="<<str.species_pp_type.size()<< ": "; for(uint i=0;i<str.species_pp_type.size();i++) *voss << str.species_pp_type.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_version.size()="<<str.species_pp_version.size()<< ": "; for(uint i=0;i<str.species_pp_version.size();i++) *voss << str.species_pp_version.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_ZVAL.size()="<<str.species_pp_ZVAL.size()<< ": "; for(uint i=0;i<str.species_pp_ZVAL.size();i++) *voss << str.species_pp_ZVAL.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_vLDAU.size()="<<str.species_pp_vLDAU.size()<< ": "; for(uint i=0;i<str.species_pp_vLDAU.size();i++) *voss << str.species_pp_vLDAU.at(i).size() << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_volume.size()="<<str.species_volume.size()<< ": "; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_mass.size()="<<str.species_mass.size()<< ": "; for(uint i=0;i<str.species_mass.size();i++) *voss << str.species_mass.at(i) << " "; *voss << endl;
      *voss << str << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [X] READY TO RETURN" << endl;
    }

    if(isTET) { aflowlib::PrototypeFixTET(*voss,str,optionsTET); }
    return str;
  }
} // namespace aflowlib

// ***************************************************************************

double NearestNeighbour(const xstructure &str_in) {
  return SYM::minimumDistance(str_in);
  //[CO 171024 OBSOLETE]xstructure str(str_in);
  //[CO 171024 OBSOLETE]// if(LDEBUG) { cerr << "NearestNeighbour 1" << endl; }
  //[CO 171024 OBSOLETE]// if(LDEBUG) { cerr << str.scale << endl; }
  //[CO 171024 OBSOLETE]str.ReScale(1.0);
  //[CO 171024 OBSOLETE]xvector<int> ndims(3);
  //[CO 171024 OBSOLETE]// str.neighbours_radius=RadiusSphereLattice(str.lattice);
  //[CO 171024 OBSOLETE]// str.neighbours_radius=max(modulus(str.lattice(1)),modulus(str.lattice(2)),modulus(str.lattice(3)));
  //[CO 171024 OBSOLETE]// ndims=LatticeDimensionSphere(str.lattice,str.neighbours_radius);
  //[CO 171024 OBSOLETE]ndims[1]=ndims[2]=ndims[3]=1;
  //[CO 171024 OBSOLETE]deque<_atom> vatoms;
  //[CO 171024 OBSOLETE]_atom atom;
  //[CO 171024 OBSOLETE]// if(LDEBUG) { cerr << "NearestNeighbour 2" << endl; }
  //[CO 171024 OBSOLETE]for(int i=-ndims[1];i<=ndims[1];i++) {
  //[CO 171024 OBSOLETE]  for(int j=-ndims[2];j<=ndims[2];j++) {
  //[CO 171024 OBSOLETE]    for(int k=-ndims[3];k<=ndims[3];k++) {
	//[CO 171024 OBSOLETE]for(uint iat=0;iat<str.atoms.size();iat++) {
	//[CO 171024 OBSOLETE]  atom=str.atoms.at(iat);
	//[CO 171024 OBSOLETE]  atom.cpos=atom.cpos+i*str.lattice(1)+j*str.lattice(2)+k*str.lattice(3);
	//[CO 171024 OBSOLETE]  vatoms.push_back(atom);
	//[CO 171024 OBSOLETE]}
  //[CO 171024 OBSOLETE]    }
  //[CO 171024 OBSOLETE]  }
  //[CO 171024 OBSOLETE]}
  //[CO 171024 OBSOLETE]double nndist=RadiusSphereLattice(str.lattice);
  //[CO 171024 OBSOLETE]for(uint i=0;i<vatoms.size();i++) {
  //[CO 171024 OBSOLETE]  for(uint j=0;j<vatoms.size();j++) {
  //[CO 171024 OBSOLETE]    if(i!=j)
	//[CO 171024 OBSOLETE]if(modulus(vatoms.at(i).cpos-vatoms.at(j).cpos) < nndist) nndist=modulus(vatoms.at(i).cpos-vatoms.at(j).cpos);
  //[CO 171024 OBSOLETE]  }
  //[CO 171024 OBSOLETE]}
  //[CO 171024 OBSOLETE]return nndist;
}

// ***************************************************************************
// ***************************************************************************
string* LOAD_Library_ICSD(string file) {  // LOAD ONE BY ONE
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(LDEBUG) { cerr << "LOAD_Library_ICSD: BEGIN" << endl; }
  
  vector<string> vn;
  vn.push_back("0");vn.push_back("1");vn.push_back("2");vn.push_back("3");vn.push_back("4");
  vn.push_back("5");vn.push_back("6");vn.push_back("7");vn.push_back("8");vn.push_back("9");
  vn.push_back("10");
  
  if(1) {
    for(uint i=1;i<vn.size();i++) {
      if(file==string("Library_ICSD"+vn.at(i)) && XHOST_vLibrary_ICSD.at(i).length()>0) {
	return &XHOST_vLibrary_ICSD.at(i);
      }
    }
  }
  if(LDEBUG) { cerr << "file=" << file << endl; }
  if(0 &&LDEBUG) {
    cerr << "XHOST_vLibrary_ICSD_ALL.size()=" << XHOST_vLibrary_ICSD_ALL.size() << endl;
    for(uint i=1;i<vn.size();i++)
      cerr << "XHOST_vLibrary_ICSD.at(" << i << ").length()=" << XHOST_vLibrary_ICSD.at(i).length() << endl;
  }
  //  DEBUG=TRUE;
  
  // *********************************************************************
  // *********************************************************************
  if(XHOST_vLibrary_ICSD_ALL.size()==0) { // LOAD from aflow_data && aflow libs
    XHOST_Library_ICSD_ALL=init::InitGlobalObject("Library_ICSD");
    aurostd::string2vectorstring(XHOST_Library_ICSD_ALL,XHOST_vLibrary_ICSD_ALL);
  }
  // *********************************************************************
  // *********************************************************************
  // now somehow I do have Library_ICSD
  // START EXTRACTING
  uint jstart=0,jstop=0;
  if(1) { //new
    for(uint i=1;i<vn.size();i++) {
      if(file==string("Library_ICSD"+vn.at(i))) { // Library_ICSDX
	for(uint j=0;j<XHOST_vLibrary_ICSD_ALL.size();j++) {
	  if(aurostd::substring2bool(XHOST_vLibrary_ICSD_ALL.at(j),string("[README_LIBRARY_ICSD"+vn.at(i)+".TXT]START"))) jstart=j;
	  if(jstart>0) if(aurostd::substring2bool(XHOST_vLibrary_ICSD_ALL.at(j),string("[README_LIBRARY_ICSD"+vn.at(i)+".TXT]STOP"))) {jstop=j;j=XHOST_vLibrary_ICSD_ALL.size();}
	}
      }
    }
  }

  if(LDEBUG) { cerr << jstart+1 << " " << jstop-1 << endl; }
  if(1) {
    for(uint i=1;i<vn.size();i++)
      if(file==string("Library_ICSD"+vn.at(i))) {
	XHOST_vLibrary_ICSD.at(i)="";  // a clean one
	for(uint j=jstart+1;j<jstop;j++)
	  if(XHOST_vLibrary_ICSD_ALL.at(j).at(0)!='/')
	    XHOST_vLibrary_ICSD.at(i)+=XHOST_vLibrary_ICSD_ALL.at(j)+" \n";
	if(LDEBUG) { cerr << file << "=" << XHOST_vLibrary_ICSD.at(i).length() << endl; }
	return &XHOST_vLibrary_ICSD.at(i);
      }  // used at(0)
    //   if(LDEBUG) { cerr << file << "=" << libaus.length() << endl; }
  }

  cerr << "aflowlib::PrototypeLibraries: AFLOW_LIBRARY not found! " << endl;
  return NULL;
}

// ***************************************************************************
// aflowlib::PrototypesHelp
// ***************************************************************************
namespace aflowlib {
  string PrototypesHelp(void) {
    stringstream strstream;
    // intro(strstream);
    strstream << endl;
    //  strstream << aflow::Banner("BANNER_TINY") << endl;
    strstream << aflow::Banner("BANNER_BIG") << endl;
    strstream << init::InitGlobalObject("README_PROTO_TXT");
    strstream << " " << endl;
    // exit(1);
    return strstream.str();
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::PrototypesIcsdHelp
// ***************************************************************************
namespace aflowlib {
  string PrototypesIcsdHelp(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::PrototypesIcsdHelp: BEGIN" << endl; }
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()>1) {
      init::ErrorOption(cout,options,"aflowlib::CALCULATED","aflow --protos_icsd[=N]");
      exit(0);
    }
    
    string *_pstrstream;
    vector<string> database;
    stringstream strstream;
    vector<uint> vnspecies;
    strstream << " " << endl;

    if(tokens.size()==0) strstream << aflow::Banner("BANNER_BIG") << endl;

    // LOAD THEM UP
    if(tokens.size()==0) {
      aurostd::string2tokens("1,2,3,4,5,6,7,8,9",vnspecies,",");
    } else {
      aurostd::string2tokens(tokens.at(0),vnspecies,",");      
    }

    // TEST **********************************
    if(!XHOST.vflag_control.flag("AFLOWLIB_SERVER")) _pstrstream=LOAD_Library_ICSD("Library_ICSD2");
    // cerr << "HERE *_pstrstream=\"" <<*_pstrstream << "\""<< endl;exit(0);
    if(!XHOST.vflag_control.flag("AFLOWLIB_SERVER") && *_pstrstream!="") { // FROM THE MACHINE **********************************
      for(uint i=0;i<vnspecies.size();i++)
	if(XHOST_vLibrary_ICSD.at(vnspecies.at(i)).length()==0) 
	  LOAD_Library_ICSD("Library_ICSD"+aurostd::utype2string(vnspecies.at(i)));
      
      for(uint i=0;i<vnspecies.size();i++) {
	_pstrstream=&LibraryDEFAULT;
	if(vnspecies.at(i)>=1 && vnspecies.at(i)<=10) {
	  _pstrstream=&XHOST_vLibrary_ICSD.at(vnspecies.at(i));
	}
	
	if(vnspecies.at(i)<=10 || vnspecies.at(i)>=11) {
	  if(vnspecies.at(i)!=11) { strstream << vnspecies.at(i) << "-COMPONENTS " << endl;}
	  else { strstream << "DEFAULT-COMPONENTS " << endl; }
	  database.clear();aurostd::string2tokens(*_pstrstream,database,"\n");
	  // strstream << vnspecies.at(i) << " " << database.size() << endl;
	  for(uint j=0;j<database.size();j++) {
	    if(aurostd::substring2bool(database.at(j),"STRUCTURE") && !aurostd::substring2bool(database.at(j),"PROTO_END")) {
	      if(vnspecies.at(i)!=11) { strstream << " " << database.at(j+1).substr(10) << "     (" << vnspecies.at(i) << ")" << endl; }
	      else { strstream << " " << database.at(j+1).substr(10) << "(DEFAULT)" << endl; }
	    }
	  }
	  strstream << "" << endl;
	}
      }
    }
    //    cerr << "HERE" << endl;exit(0);
    //   if((!XHOST.vflag_control.flag("AFLOWLIB_SERVER") && _pstrstream==NULL) || XHOST.vflag_control.flag("AFLOWLIB_SERVER")) { // FROM THE server **********************************
    if((!XHOST.vflag_control.flag("AFLOWLIB_SERVER") && *_pstrstream=="") || XHOST.vflag_control.flag("AFLOWLIB_SERVER")) { // FROM THE server **********************************
      if(!XHOST.vflag_control.flag("AFLOWLIB_SERVER") && !XHOST.vflag_control.flag("AFLOWLIB_SERVER"))
	XHOST.vflag_control.push_attached("AFLOWLIB_SERVER",AFLOWLIB_SERVER_DEFAULT); // some default
      cerr << "aflowlib::PrototypeLibraries: Contacting \"" <<  XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER") << "\"" << endl;
      vector<string> vlattice,speciesX,tokens;aurostd::string2tokens("BCC,BCT,CUB,FCC,HEX,MCL,MCLC,ORC,ORCC,ORCF,ORCI,RHL,TET,TRI",vlattice,",");
      // vector<string> vlattice,speciesX,tokens;aurostd::string2tokens("BCC,BCT,CUB",vlattice,",");
      vector<double> vnatomsX;
      vector<vector<string> > voutput;for(uint i=0;i<10;i++) voutput.push_back(*(new vector<string>(0)));
      aflowlib::_aflowlib_entry data;
      vector<aflowlib::_aflowlib_entry> vdata;
      
      for(uint i=0;i<vlattice.size();i++) {
	data.clear();
	cerr << "aflowlib::PrototypeLibraries: Scanning .. " << vlattice.at(i) << " ";
	data.url2aflowlib(XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")+":AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/?format=text",cout,FALSE);vdata.push_back(data);
	for(uint j=0;j<data.vaflowlib_entries.size();j++) {
	  aurostd::string2tokens(data.vaflowlib_entries.at(j),tokens,"_");
	  KBIN::VASP_SplitAlloySpecies(tokens.at(0),speciesX,vnatomsX);
	  if(speciesX.size()<=voutput.size()) {
	    voutput.at(speciesX.size()-1).push_back(tokens.at(0)+" ["+tokens.at(0)+"] "+vlattice.at(i)+" "+data.vaflowlib_entries.at(j)+" ICSD_"+tokens.at(tokens.size()-1)+"       ("+aurostd::utype2string(speciesX.size())+")");
	  }
	}
	cerr << "[" << data.vaflowlib_entries.size() << " entries]" << endl;
      }
      for(uint j=0;j<voutput.size();j++) {
	aurostd::sort(voutput.at(j));
      }
      for(uint i=0;i<vnspecies.size();i++) {
	for(uint j=0;j<voutput.size();j++) {
	  if(j+1==vnspecies.at(i) || vnspecies.at(i)==0) {
	    if(voutput.at(j).size()>0) {
	      strstream << j+1 << "-COMPONENTS " << endl;
	      for(uint k=0;k<voutput.at(j).size();k++) {
		strstream << voutput.at(j).at(k) << endl;
	      }
	      strstream << endl;
	    }
	  }
	}
      }
    }
    if(LDEBUG) { cerr << "aflowlib::PrototypesIcsdHelp: END" << endl; }
    return strstream.str();

  }
} // namespace aflowlib
// B1H3In1Na1O10P2 [B1H3In1Na1O10P2] mS72 C12/c1 NaIn(B(P2O7)(OH)3) 15 B1H3In1Na1O10P2_ICSD_409585 ICSD_409585       (6)

// ***************************************************************************
// aflowlib::CALCULATED
// ***************************************************************************
namespace aflowlib {
  string CALCULATED(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::CALCULATED: BEGIN" << endl; }
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()>1) {
      init::ErrorOption(cout,options,"aflowlib::CALCULATED","aflow --calculated[[=]all|icsd|lib1|lib2|lib3|lib4|lib5|lib6|lib7|lib8|lib9]");
      exit(0);
    }

    string library;
    if(tokens.size()==0) library="all";
    if(tokens.size()==1) library=tokens.at(0);
    library=aurostd::tolower(library);
    
    if(LDEBUG) { cerr << "aflowlib::CALCULATED: library=" << library << endl; }

    stringstream out;
    // LIBRARY_ICSD or LIBRARY_LIBN
    if(aurostd::substring2bool(library,"icsd") || 
       aurostd::substring2bool(library,"lib1") || 
       aurostd::substring2bool(library,"lib2") || 
       aurostd::substring2bool(library,"lib3") || 
       aurostd::substring2bool(library,"lib4") || 
       aurostd::substring2bool(library,"lib5") || 
       aurostd::substring2bool(library,"lib6") || 
       aurostd::substring2bool(library,"lib7") || 
       aurostd::substring2bool(library,"lib8") || 
       aurostd::substring2bool(library,"lib9")) {
      vector<string> ltokens,vLibrary;
      if(aurostd::substring2bool(library,"icsd")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_ICSD_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib1")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB1_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib2")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB2_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib3")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB3_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib4")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB4_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib5")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB5_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib6")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB6_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib7")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB7_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib8")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB8_LIB"),ltokens,"\n");
      if(aurostd::substring2bool(library,"lib9")) aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB9_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary.push_back(ltokens.at(i));
      if(XHOST.vflag_control.flag("REMOVE")) {
	for(int i=vLibrary.size()-1;i>=0;i--)
	  out << "rm -rfv ./" << vLibrary.at(i)<< endl;

      } else {
	if(XHOST.vflag_control.flag("ZIP")) {
	  for(int i=vLibrary.size()-1;i>=0;i--)
	    out << "zip -9rmv already_done.zip ./" << vLibrary.at(i)<< endl;
	
	} else {
	  for(int i=vLibrary.size()-1;i>=0;i--)
	    //    out << i+1 << " " << vLibrary.at(i)<< endl; // OLD
	    out << i+1 << "/" << vLibrary.size() << "   " << vLibrary.at(i)<< endl; // NEW
	}
      }
    }
    // out << "XHOST.Progname=" << XHOST.Progname << endl;
    // out << "XHOST.command(\"aflow_data\")=" << XHOST.command("aflow_data") << endl;
    // LIBRARY_ALL
    if(aurostd::substring2bool(library,"all")) {
      out << "AFLOW V"<<string(AFLOW_VERSION)<< "  "  << TODAY << endl;
      vector<string> ltokens,vLibrary_CALCULATED_ICSD,
	vLibrary_CALCULATED_LIB1,
	vLibrary_CALCULATED_LIB2,
	vLibrary_CALCULATED_LIB3,
	vLibrary_CALCULATED_LIB4,
	vLibrary_CALCULATED_LIB5,
	vLibrary_CALCULATED_LIB6,
	vLibrary_CALCULATED_LIB7,
	vLibrary_CALCULATED_LIB8,
	vLibrary_CALCULATED_LIB9;
      
      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_ICSD_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_ICSD.push_back(ltokens.at(i));
      out << "Library_ICSD_CALCULATED = " << vLibrary_CALCULATED_ICSD.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB1_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB1.push_back(ltokens.at(i));
      out << "Library_LIB1_CALCULATED = " << vLibrary_CALCULATED_LIB1.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB2_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB2.push_back(ltokens.at(i));
      out << "Library_LIB2_CALCULATED = " << vLibrary_CALCULATED_LIB2.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB3_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB3.push_back(ltokens.at(i));
      out << "Library_LIB3_CALCULATED = " << vLibrary_CALCULATED_LIB3.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB4_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB4.push_back(ltokens.at(i));
      out << "Library_LIB4_CALCULATED = " << vLibrary_CALCULATED_LIB4.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB5_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB5.push_back(ltokens.at(i));
      out << "Library_LIB5_CALCULATED = " << vLibrary_CALCULATED_LIB5.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB6_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB6.push_back(ltokens.at(i));
      out << "Library_LIB6_CALCULATED = " << vLibrary_CALCULATED_LIB6.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB7_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB7.push_back(ltokens.at(i));
      out << "Library_LIB7_CALCULATED = " << vLibrary_CALCULATED_LIB7.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB8_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB8.push_back(ltokens.at(i));
      out << "Library_LIB8_CALCULATED = " << vLibrary_CALCULATED_LIB8.size() << endl;

      aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_LIB9_LIB"),ltokens,"\n");
      for(uint i=0;i<ltokens.size();i++)
	if(aurostd::substring2bool(ltokens.at(i),"/"))
	  vLibrary_CALCULATED_LIB9.push_back(ltokens.at(i));
      out << "Library_LIB9_CALCULATED = " << vLibrary_CALCULATED_LIB9.size() << endl;

      out << "AFLOW V" 
	  << string(AFLOW_VERSION)
	  << " [date="  <<  aurostd::get_datetime() << "]" 
	  << " [built="  << TODAY << "]" 
	  << " Calcs=" 
	  << vLibrary_CALCULATED_ICSD.size() << "," 
	  << vLibrary_CALCULATED_LIB1.size() << "," 
	  << vLibrary_CALCULATED_LIB2.size() << "," 
	  << vLibrary_CALCULATED_LIB3.size() << "," 
	  << vLibrary_CALCULATED_LIB4.size() << "," 
	  << vLibrary_CALCULATED_LIB5.size() << "," 
	  << vLibrary_CALCULATED_LIB6.size() << "," 
	  << vLibrary_CALCULATED_LIB7.size() << "," 
	  << vLibrary_CALCULATED_LIB8.size() << "," 
	  << vLibrary_CALCULATED_LIB9.size() << "=" << 
	vLibrary_CALCULATED_ICSD.size() +
	vLibrary_CALCULATED_LIB1.size() +
	vLibrary_CALCULATED_LIB2.size() +
	vLibrary_CALCULATED_LIB3.size() +
	vLibrary_CALCULATED_LIB4.size() +
	vLibrary_CALCULATED_LIB5.size() +
	vLibrary_CALCULATED_LIB6.size() +
	vLibrary_CALCULATED_LIB7.size() +
	vLibrary_CALCULATED_LIB8.size() +
	vLibrary_CALCULATED_LIB9.size() << endl;

    }
    if(LDEBUG) { cerr << "aflowlib::CALCULATED: END" << endl; }
    return out.str();
  }
}  // namespace aflowlib

// ***************************************************************************
// Aflowlib_CALCULATED_ICSD_RANDOM
// ***************************************************************************
namespace aflowlib {
  string CALCULATED_ICSD_RANDOM(void) {
    stringstream out;
    vector<string> tokens,vLibrary;
    aurostd::string2tokens(init::InitGlobalObject("Library_CALCULATED_ICSD_LIB"),tokens,"\n");
    for(uint i=0;i<tokens.size();i++)
      if(aurostd::substring2bool(tokens.at(i),"/"))
	vLibrary.push_back(tokens.at(i));

    uint index=(uint) floor(vLibrary.size()*aurostd::ran0());
    aurostd::string2tokens(vLibrary.at(index),tokens,"/");
    out << vLibrary.at(index) << "/" << tokens.at(1) << endl;

    return out.str();
  }
}  // namespace aflowlib


// ***************************************************************************

// ***************************************************************************

#endif // _AFLOW_XPROTO_CPP
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************



// ***************************************************************************// ***************************************************************************
// ***************************************************************************// ***************************************************************************
// ***************************************************************************// ***************************************************************************
// the mother of all the prototypes
namespace aflowlib {
  xstructure _____PrototypeLibraries(ostream &oss,_PROTO_PARAMS *PARAMS) { // COMPLETE ONE
    //  DEBUG=TRUE;
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries(ostream &oss,_PROTO_PARAMS *PARAMS)" << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: PARAMS->label=" << PARAMS->label << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: PARAMS->volume_in=" << PARAMS->volume_in << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: PARAMS->mode=" << PARAMS->mode << endl; }
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: PARAMS->flip_option=" << PARAMS->flip_option << endl; }
   
    deque<string> vatomX_backup(PARAMS->vatomX);
    deque<double> vvolumeX_backup(PARAMS->vvolumeX);

    xstructure str("");str.lattice.clear();
    string optionsTET;bool isTET=FALSE;

    if(PARAMS->label.length()>1 && !aurostd::substring2bool(PARAMS->label,"_ICSD_") &&
       (PARAMS->label[0]=='F' || PARAMS->label[0]=='f' || PARAMS->label[0]=='B' || PARAMS->label[0]=='b' || PARAMS->label[0]=='H' || PARAMS->label[0]=='h' || PARAMS->label[0]=='S' || PARAMS->label[0]=='s') &&
       (PARAMS->label[1]!='c') && (PARAMS->label[1]!='C')) {
      str=aflowlib::PrototypeBinaryGUS(oss,PARAMS->label,PARAMS->vatomX.at(0),PARAMS->vvolumeX.at(0),PARAMS->vatomX.at(1),PARAMS->vvolumeX.at(1),PARAMS->volume_in);   // found FCC,BCC,HCP
      return str;
    }

    ostream *voss;if(LDEBUG || XHOST.DEBUG) voss=&cerr; else voss=&oss;

    vector<uint> natomsX;
    string structure_name="";
    str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
    bool found=FALSE;
    std::string label_swap="",index_swap="",label_use="",label_remove="",label_species="",string_tmp;
    string *_pstrstream;
    vector<string> tokens,ttokens,database,speciesX;
    vector<double> vnatomsX;
    stringstream aus;
    ostringstream oaus;
    uint kstructure=0,kprototype=0,kinfo=0,kmode=0,klattice=0,katoms=0,kstop=0;
    uint iat,natoms=0,nspecies=0;
    double volume=0.0;
    int mode_load=0,mode_prim=0,mode_conventional=0,mode_volume=0,mode_remove=0,mode_species=0;
    if(LDEBUG) { cerr << "aflowlib::PrototypeLibraries: XHOST_Library_HTQC.length()=" << XHOST_Library_HTQC.length() << endl; }
    init::InitGlobalObject("Library_HTQC");

    // *********************************************************************
    speciesX.clear();natomsX.clear();
    str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
    nspecies=0;
    _pstrstream=&XHOST_Library_HTQC; // default is HTQC so it does not complain
    // *********************************************************************
    // PARSE and STORE the file
    // LDEBUG=TRUE;
    // ICSD
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [aflowlib::PrototypeLibraries] [0]" << endl;
    //  *voss << "PARAMS->MODE=" << PARAMS->mode << endl; exit(0);
    found=FALSE; // start the search
    // *********************************************************************
    if(PARAMS->mode==LIBRARY_MODE_ICSD) {
      // recognize which LIBRARY of ICSD
      if(LDEBUG) if(PARAMS->mode==LIBRARY_MODE_ICSD) *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->mode==LIBRARY_MODE_ICSD" << endl;

      string icsd_label=PARAMS->label;
      if(aurostd::substring2bool(icsd_label,"_ICSD_")) {
	aurostd::string2tokens(icsd_label,tokens,"_");
	icsd_label=tokens.at(0);
      }
      nspecies=KBIN::VASP_SplitAlloySpecies(aurostd::RemoveSubStringFirst(icsd_label,"C-"),speciesX,vnatomsX);

      //    for(uint i=0;i<speciesX.size();i++) if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) " << nspecies << " " << speciesX.at(i) << " " << vnatomsX.at(i) << endl; // some debug always helps exit(0);
      if(aurostd::substring2bool(LibraryDEFAULT,PARAMS->label)) _pstrstream=&LibraryDEFAULT;
      else {
	if(nspecies<1||nspecies>10) {
	  _pstrstream=LOAD_Library_ICSD("Library_ICSD2");//default
	} else {
	  _pstrstream=LOAD_Library_ICSD("Library_ICSD"+aurostd::utype2string<int>(nspecies));
	}
      }
      PARAMS->vatomX.clear();PARAMS->vvolumeX.clear();
      for(uint i=0;i<nspecies;i++) {
	str.species.push_back(speciesX.at(i));
	str.species_pp.push_back(speciesX.at(i));
	str.species_pp_type.push_back("");
	str.species_pp_version.push_back("");
	str.species_pp_ZVAL.push_back(0.0);
	str.species_pp_vLDAU.push_back(deque<double>());
	str.species_volume.push_back(0.0);
	str.species_mass.push_back(0.0);
	PARAMS->vatomX.push_back(speciesX.at(i));
	PARAMS->vvolumeX.push_back(0.0);
	natomsX.push_back(0);
	// str.species.push_back(speciesX.at(i));PARAMS->vatomX.push_back("");PARAMS->vvolumeX.push_back(0.0);natomsX.push_back(0);
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [1] nspecies=" << nspecies << endl;
    }
    // *********************************************************************
    // HTQC and ICSD
    if(PARAMS->mode==LIBRARY_MODE_HTQC_ICSD) {
      //   cerr << "MODE ICSD" << endl;
      string line_ICSD,label_compound;
      string label_prefix=PARAMS->label;
      string label_postfix="",label_icsd="";

      if(aurostd::substring2bool(label_prefix,".")) {
	aurostd::string2tokens(label_prefix,tokens,".");
	if(tokens.size()>0) label_prefix=tokens.at(0);
	if(tokens.size()>0) label_icsd=tokens.at(0);
	if(tokens.size()>0) label_postfix=tokens.at(1);
      }
      // cerr << label_prefix << endl;
      // cerr << label_postfix << endl;

      if(LDEBUG) if(PARAMS->mode==LIBRARY_MODE_HTQC_ICSD) *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->mode==LIBRARY_MODE_HTQC_ICSD" << endl;
      if(aurostd::substring2bool(label_prefix,"ICSD_") && !aurostd::substring2bool(label_prefix,"_ICSD_")) {  // LOAD THEM ALL
	bool found_HTQC_ICSD=FALSE;
	for(uint n=1;n<=10&&!found_HTQC_ICSD;n++) {
	  aurostd::string2tokens(*LOAD_Library_ICSD("Library_ICSD"+aurostd::utype2string<uint>(n)),database,"\n");
	  for(uint i=0;i<database.size()&&!found_HTQC_ICSD;i++)
	    if(aurostd::substring2bool(database[i],"STRUCTURE")) {
	      aurostd::string2tokens(database[i],tokens);
	      for(uint j=0;j<tokens.size();j++)
		if(label_prefix==tokens.at(j)) {
		  found_HTQC_ICSD=TRUE;//   cerr << database[i] << endl;
		  line_ICSD=database[i+1];
		  for(uint k=0;k<tokens.size();k++)
		    if(aurostd::substring2bool(tokens.at(k),"_ICSD_")) {
		      label_prefix=tokens.at(k);
		      label_compound=tokens.at(1);
		      //    cerr << label_compound << endl;
		    }
		}
	    }
	}
      }

      if(label_prefix!="" && line_ICSD!="") {
	found=TRUE;
	deque<string> vatomX_tmp;
	deque<double> vvolumeX_tmp;
	deque<string> abc;
	str=aflowlib::PrototypeLibraries(oss,label_prefix,PARAMS->parameters,vatomX_tmp,vvolumeX_tmp,PARAMS->volume_in,LIBRARY_MODE_ICSD,PARAMS->flip_option);

	// need to fix this because I could get stuff from WYCKOFIZATION

	//  cerr << str << endl;exit(0);
	if(vatomX_tmp.size()!=label_postfix.size()) {
	  *voss << "ERROR - aflowlib::PrototypeLibraries:  your label needs to contain the right species number/order, i.e. " << label_icsd << ".ABC..." << endl;
	  *voss << "ERROR - aflowlib::PrototypeLibraries:  PARAMS->vatomX.size()!=label_postfix.size() (" << vatomX_tmp.size() << "!=" << label_postfix.size() << ")" << endl;
	  exit(0);}
	// create the mapping

	vector<int> mapij;
	for(uint i=0;i<label_postfix.size();i++) {
	  mapij.push_back(-1);
	  string digit=" ";
	  for(uint j=0;j<26;j++) {
	    digit[0]=65+j; // 65 is "A"
	    if(label_postfix.at(i)==digit[0])
	      mapij.at(mapij.size()-1)=j;
	  }
	}
	// cerr << "mapij.size()=" << mapij.size() << endl; for(uint i=0;i<mapij.size();i++) cerr << "mapij.at(" << i << ")=" << mapij.at(i) << endl;
	// exit(0);
	abc.clear();
	if(PARAMS->vatomX.size()==0 && PARAMS->vvolumeX.size()==0) for(uint i=0;i<vatomX_tmp.size();i++) {string digit=" ";digit[0]=65+mapij.at(i);abc.push_back(digit);}
	if(PARAMS->vatomX.size()==str.species.size() && (PARAMS->vvolumeX.size()==0 || PARAMS->vvolumeX.size()==str.species.size())) for(uint i=0;i<vatomX_tmp.size();i++)  {abc.push_back(PARAMS->vatomX.at(mapij.at(i)));}
	//      cerr << label_prefix << endl;


	if(PARAMS->vatomX.size()!=vatomX_tmp.size()) {
	  *voss << "EXCEPTION FOUND - aflowlib::PrototypeLibraries:  PARAMS->vatomX.size()!=vatomX_tmp.size() (" << PARAMS->vatomX.size() << "!=" << vatomX_tmp.size() << ") ... FIXING" << endl;
	  vatomX_tmp.clear(); for(uint i=0;i<PARAMS->vatomX.size();i++) vatomX_tmp.push_back(PARAMS->vatomX.at(i));
	}
	// /*
	// if(PARAMS->vvolumeX.size()!=PARAMS->vvolumeX_tmp.size()) {
	// 	*voss << "EXCEPTION FOUND - aflowlib::PrototypeLibraries:  PARAMS->vvolumeX.size()!=PARAMS->vvolumeX_tmp.size() (" << PARAMS->vvolumeX.size() << "!=" << PARAMS->vvolumeX_tmp.size() << ") ... FIXING" << endl;
	// 	//	vvolumeX_tmp.clear(); for(uint i=0;i<PARAMS->vvolumeX.size();i++) vvolumeX_tmp.push_back(PARAMS->vvolumeX.at(i));
	// }
	// */
	if(PARAMS->vatomX.size()>0 && PARAMS->vatomX.size()!=vatomX_tmp.size()) {
	  *voss << "ERROR - aflowlib::PrototypeLibraries:  PARAMS->vatomX.size()!=vatomX_tmp.size() (" << PARAMS->vatomX.size() << "!=" << vatomX_tmp.size() << ")" << endl;
	  exit(0);
	}

	// fil A,B,C ..
	if(PARAMS->vatomX.size()==0 && PARAMS->vvolumeX.size()==0) { // fill A B C
	  str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
	  for(uint i=0;i<vatomX_tmp.size()&&i<abc.size();i++) {
	    str.species.push_back(abc.at(i));str.species_pp.push_back(abc.at(i));str.species_pp_type.push_back("");str.species_pp_version.push_back("");str.species_pp_ZVAL.push_back(0.0);
	    str.species_pp_vLDAU.push_back(deque<double>());str.species_volume.push_back(1.0);str.species_mass.push_back(0.0);
	  }
	}
	// fil real names/volumes
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.species.size()=" << str.species.size() << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vvolumeX.size()=" << PARAMS->vvolumeX.size() << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.GetVolume()=" << str.GetVolume() << endl;

	if(PARAMS->vatomX.size()==str.species.size() && (PARAMS->vvolumeX.size()==0 || PARAMS->vvolumeX.size()==str.species.size())) { // fill A B C
	  //	cerr << "IN" << endl;
	  str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
	  for(uint i=0;i<PARAMS->vatomX.size()&&i<abc.size();i++) {
	    str.species.push_back(abc.at(i));str.species_pp.push_back(abc.at(i));str.species_pp_type.push_back("");str.species_pp_version.push_back("");str.species_pp_ZVAL.push_back(0.0);
	    str.species_pp_vLDAU.push_back(deque<double>());str.species_mass.push_back(0.0);
	    str.species_volume.push_back(GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i))));
	    //	  cerr << "GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i)))=" << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i))) << endl;
	  }
	}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl; for(uint i=0;i<PARAMS->vatomX.size();i++) cerr << "PARAMS->vatomX.at(" << i << ")=" << PARAMS->vatomX.at(i) << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(PARAMS->vatomX.at(i))) << endl;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) abc.size()=" << abc.size() << endl; for(uint i=0;i<abc.size();i++) cerr << "abc.at(" << i << ")=" << abc.at(i) << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(abc.at(i))) << endl;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) volumeX.size()=" << PARAMS->vvolumeX.size() << endl; for(uint i=0;i<PARAMS->vvolumeX.size();i++) cerr << "PARAMS->vvolumeX.at(" << i << ")=" << PARAMS->vvolumeX.at(i) << endl;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.species.size()=" << str.species.size() << endl; for(uint i=0;i<str.species.size();i++) cerr << "str.species.at(" << i << ")=" << str.species.at(i) << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(str.species.at(i))) << endl;}
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.atoms.size()=" << str.atoms.size() << endl; for(uint i=0;i<str.atoms.size();i++) cerr << "str.atoms.at(" << i << ").type=" << str.atoms.at(i).type << " " << GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(str.atoms.at(i).name)) << endl;}
	//

	double volume=0.0;
	for(uint iat=0;iat<str.atoms.size();iat++) {
	  if((uint) str.atoms.at(iat).type<str.species.size()) {
	    str.atoms.at(iat).name=str.species.at(str.atoms.at(iat).type);
	    str.atoms.at(iat).name_is_given=TRUE;
	    volume+=str.species_volume.at(str.atoms.at(iat).type);
	    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) volume=" << volume << "  (" << str.atoms.at(iat).type << ")" << endl;
	  }
	}
	if(abs(volume)<1.0e-6) { // wrap up to fix it
	  *voss << "EXCEPTION FOUND - aflowlib::PrototypeLibraries:  volume==0 ... FIXING" << endl;
	  for(uint iat=0;iat<PARAMS->vvolumeX.size();iat++) volume+=PARAMS->vvolumeX.at(iat);
	}

	//  if(PARAMS->vatomX.size()==0 && PARAMS->vvolumeX.size()==0)
	{  // fix volume in WYC
	  str.scale=std::pow((double) (abs(volume)/det(str.lattice)),(double) 1.0/3.0);
	  str.neg_scale=TRUE;
	}
	// strip unstrip atoms
	std::deque<_atom> atoms;
	atoms=str.atoms;
	while(str.atoms.size()>0) str.RemoveAtom(0);
	for(uint i=0;i<atoms.size();i++) str.AddAtom(atoms.at(i));
	// make alphabetic
	str.SpeciesPutAlphabetic();
	// fix title
	//     aurostd::StringSubst(str.title,label_icsd,label_icsd+"."+label_postfix);
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(1): " << str.title << endl;
	aurostd::StringSubst(str.title,label_icsd,"");aurostd::StringSubst(str.title,label_compound+"_","");aurostd::StringSubst(str.title,"["+label_compound+"]","");
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(2): " << str.title << endl;
	aurostd::StringSubst(str.title,label_compound,"");aurostd::StringSubst(label_compound,"1","");aurostd::StringSubst(str.title,label_compound,"");
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(3): " << str.title << endl;
	aurostd::StringSubst(str.title,"()","");aurostd::StringSubst(str.title,"(icsd library)","");aurostd::StringSubst(str.title,"  "," ");
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(4): " << str.title << endl;
	string aus=" "+str.title;aurostd::StringSubst(str.title,"  "," ");
	str.title="";
	for(uint i=0;i<PARAMS->vatomX.size();i++) str.title+=PARAMS->vatomX.at(i);
	str.title+="/"+label_icsd+"."+label_postfix+aus;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(5): " << str.title << endl;
	//      cerr << str.title << endl;exit(0);
	// done
	return str;
      }
    }

    // *********************************************************************
    // HTQC
    // *voss << "PARAMS->MODE=" << PARAMS->mode << endl; exit(0);
    if(PARAMS->mode==LIBRARY_MODE_HTQC) { // mode =LIBRARY_MODE_HTQC
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->mode==LIBRARY_MODE_HTQC" << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl;
      // check if binary ternary or so on..
      uint nspeciesHTQC=aflowlib::PrototypeLibrariesSpeciesNumber(PARAMS->label);
      if(nspeciesHTQC>PARAMS->vatomX.size()) {
	if(LDEBUG) *voss << "gotta fix if: nspeciesHTQC=" << nspeciesHTQC << " > PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl;
	for(uint i=PARAMS->vatomX.size();i<nspeciesHTQC;i++) {
	  string_tmp="A";string_tmp[0]+=i;
	  PARAMS->vatomX.push_back(string_tmp);
	  PARAMS->vvolumeX.push_back(PARAMS->vvolumeX.at(0));
	}
      }
      // fixed now build up
      for(uint i=0;i<PARAMS->vatomX.size();i++) {
	if(LDEBUG) *voss << "DEBUG Building up speciesX and natomsX. We have PARAMS->vatomX.at(i)=" << PARAMS->vatomX.at(i) << endl;
	string_tmp="A";string_tmp[0]+=i; // do ABC... in ascii
	speciesX.push_back(string_tmp);
	// str.species.push_back(PARAMS->vatomX.at(i));
	natomsX.push_back(0); // start from NOTHING
      }
      _pstrstream=&XHOST_Library_HTQC; // default is HTQC
    }
    database.clear();aurostd::string2tokens(*_pstrstream,database,"\n");
    // FIND STRUCTURE kstructure
    string label_library="";
    kstructure=0;kstop=0;

    for(uint i=0;i<database.size();i++)
      if(aurostd::substring2bool(database[i],"STRUCTURE"))
	if(aurostd::substring2bool(database[i]," "+PARAMS->label+" ")) { // PARAMS->label or " "+PARAMS->label+" "
	  aurostd::string2tokens(database[i],tokens);
	  for(uint j=0;j<tokens.size();j++)
	    if(tokens[j]==PARAMS->label) {found=TRUE;kstructure=i;label_library=tokens[1];}
	}

    // LOOKING FOR TET MODIFICATIONS
    if(!found) {
      if(aurostd::substring2bool(PARAMS->label,"-tet")) {
	isTET=TRUE;
	if(isTET && LDEBUG) cerr << "TET modification" << endl;
	aurostd::string2tokens(PARAMS->label,tokens,"-tet");
	string labelAUS=tokens.at(0);
	optionsTET=PARAMS->label;aurostd::StringSubst(optionsTET,tokens.at(0),"");
	if(aurostd::substring2bool(PARAMS->label,"%")) {aurostd::string2tokens(PARAMS->label,tokens,"%");labelAUS+=tokens.at(1);aurostd::StringSubst(optionsTET,tokens.at(1),"");}
	aurostd::StringSubst(optionsTET,"-tet","");
	if(LDEBUG) *voss << "LOADING PROTO=" << labelAUS << endl;
	if(LDEBUG) *voss << "optionsTET=" << optionsTET << endl;
	for(uint i=0;i<database.size();i++)
	  if(aurostd::substring2bool(database[i],"STRUCTURE"))
	    if(aurostd::substring2bool(database[i]," "+labelAUS+" ")) { // PARAMS->label or " "+PARAMS->label+" "
	      aurostd::string2tokens(database[i],tokens);
	      for(uint j=0;j<tokens.size();j++)
		if(tokens[j]==labelAUS) {found=TRUE;kstructure=i;label_library=PARAMS->label;}   // label_library=tokens[1];
	    }
      }
    }

    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) FOUND" << endl;
    // IF NOT FOUND
    if(!found) {
      if(PARAMS->mode==LIBRARY_MODE_HTQC)
	*voss << "ERROR - aflowlib::PrototypeLibraries:  Structure not present in the HTQC database " << _HTQC_PROJECT_STRING_ << "  (" << PARAMS->label << ")" << endl;
      if(PARAMS->mode==LIBRARY_MODE_ICSD || PARAMS->mode==LIBRARY_MODE_HTQC_ICSD) {
	*voss << "ERROR - aflowlib::PrototypeLibraries:  Structure not present in the ICSD database " << _ICSD_STRING_ << "  (" << PARAMS->label << ")" << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(1).length()=" << XHOST_vLibrary_ICSD.at(1).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(2).length()=" << XHOST_vLibrary_ICSD.at(2).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(3).length()=" << XHOST_vLibrary_ICSD.at(3).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(4).length()=" << XHOST_vLibrary_ICSD.at(4).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(5).length()=" << XHOST_vLibrary_ICSD.at(5).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(6).length()=" << XHOST_vLibrary_ICSD.at(6).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(7).length()=" << XHOST_vLibrary_ICSD.at(7).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(8).length()=" << XHOST_vLibrary_ICSD.at(8).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(9).length()=" << XHOST_vLibrary_ICSD.at(9).length() << endl;
	*voss  << "XHOST_vLibrary_ICSD.at(10).length()=" << XHOST_vLibrary_ICSD.at(10).length() << endl;
      }
      exit(0);
    }
    // FIND STRUCTURE kstop
    kstop=0;
    for(uint i=kstructure+1;i<database.size() && !kstop;i++) {
      if(aurostd::substring2bool(database[i],"STRUCTURE")) {
	kstop=i;
      }
    }

    if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kstructure=" << kstructure << endl;
    if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kstop=" << kstop << endl;
    // FIND STRUCTURE katoms
    klattice=kstructure;katoms=kstructure; // but must be fixed !
    for(uint i=kstructure;i<kstop;i++) {
      //  cerr << database[i] << endl;
      if(aurostd::substring2bool(database[i],"PROTOTYPE")) kprototype=i;
      if(aurostd::substring2bool(database[i],"INFO")) kinfo=i;
      if(aurostd::substring2bool(database[i],"PARAMS->MODE")) kmode=i;
    }
    if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kprototype=" << kprototype << endl;
    if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kinfo=" << kinfo << endl;
    if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) kmode=" << kmode << endl;
    // FIX PROTOTYPE
    aurostd::string2tokens(database[kprototype],tokens);
    for(uint j=1;j<tokens.size();j++) str.prototype+=tokens[j]+" ";
    // FIX INFO
    aurostd::string2tokens(database[kinfo],tokens);
    for(uint j=1;j<tokens.size();j++) str.info+=tokens[j]+" ";
    str.info+=_HTQC_PROJECT_STRING_;
    // FIX MODE
    if(LDEBUG) *voss << database[kmode] << endl;
    aurostd::string2tokens(database[kmode],tokens," ");
    mode_prim=STRUCTURE_MODE_NONE;
    mode_conventional=STRUCTURE_MODE_NONE;
    mode_volume=STRUCTURE_MODE_NONE;
    mode_load=STRUCTURE_MODE_NONE;
    mode_remove=STRUCTURE_MODE_NONE;
    for(uint j=1;j<tokens.size();j++) {
      aurostd::string2tokens(tokens.at(j),ttokens,"=");
      if(LDEBUG) *voss << tokens.at(j) << endl;
      if(tokens.at(j)=="VOLUME" || tokens.at(j)=="VOL") {
	mode_volume=STRUCTURE_MODE_VOLUME;}
      if(tokens.at(j)=="PRIM") {
	mode_prim=STRUCTURE_MODE_PRIM;}
      if(tokens.at(j)=="CONVENTIONAL" || tokens.at(j)=="CONV") {
	mode_conventional=STRUCTURE_MODE_CONVENTIONAL;}
      if(tokens.at(j)=="RAW") {
	mode_load=STRUCTURE_MODE_RAW;}
      if(tokens.at(j)=="ABC") {
	mode_load=STRUCTURE_MODE_ABC;}
      if(tokens.at(j)=="WYC") {
	mode_load=STRUCTURE_MODE_WYC;}
      if(tokens.at(j)=="WYCK_ICSD" || tokens.at(j)=="WYC_ICSD" || tokens.at(j)=="ICSD") {
	mode_load=STRUCTURE_MODE_ICSD;mode_volume=STRUCTURE_MODE_VOLUME;}
      if(aurostd::substring2bool(tokens.at(j),"USE=")) {
	mode_load=STRUCTURE_MODE_USE;aurostd::string2tokens(tokens.at(j),ttokens,"=");label_use=ttokens.at(1);}
      if(tokens.at(j)=="REMOVE") {
	mode_remove=STRUCTURE_MODE_REMOVE;label_remove=tokens.at(j+1);}
      if(aurostd::substring2bool(tokens.at(j),"SPECIES=")) {
	mode_species=STRUCTURE_MODE_SPECIES;aurostd::string2tokens(tokens.at(j),ttokens,"=");label_species=ttokens.at(1);}
      if(aurostd::substring2bool(tokens.at(j),"SWAP_AB=") || aurostd::substring2bool(tokens.at(j),"SWAP_BA=")) {
	mode_load=STRUCTURE_MODE_SWAP_AB;aurostd::string2tokens(tokens.at(j),ttokens,"=");label_swap=ttokens.at(1);}
      if(aurostd::substring2bool(tokens.at(j),"SWAP(") && aurostd::substring2bool(tokens.at(j),")=")) {
	mode_load=STRUCTURE_MODE_SWAP_XY;
	aurostd::string2tokens(tokens.at(j),ttokens,"=");
	index_swap=ttokens.at(0);aurostd::StringSubst(index_swap,"SWAP(","");aurostd::StringSubst(index_swap,")","");
	label_swap=ttokens.at(1);
      }
    }

    if(mode_prim==STRUCTURE_MODE_PRIM) if(LDEBUG) *voss << "(mode_prim==STRUCTURE_MODE_PRIM)" << endl;
    if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL) if(LDEBUG) *voss << "(mode_conventional==STRUCTURE_MODE_CONVENTIONAL)" << endl;
    if(mode_load==STRUCTURE_MODE_NONE) {
      *voss << "invalid mode" << endl;
      *voss << "rerun with --debug" << endl;
      exit(0);
    }

    if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) mode=" << mode_load << endl;

    // *voss << "HERE" << endl;exit(0);

    // FIX TITLE
    if(mode_load!=STRUCTURE_MODE_ICSD) {
      str.title="";
      for(uint i=0;i<PARAMS->vatomX.size();i++) str.title+=PARAMS->vatomX.at(i);
      str.title+="/"+label_library+" - "+"("+PARAMS->label+")";
      str.title+=" - "+str.prototype+" "+_HTQC_PROJECT_STRING_+" ";
    }
    if(mode_load==STRUCTURE_MODE_ICSD) {
      str.title=label_library+" - "+"("+PARAMS->label+")";
      str.title+=" - "+str.prototype+" "+_ICSD_STRING_;
    }

    // --------------------------------------------------------------------------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_USE && mode_species!=STRUCTURE_MODE_SPECIES) {
      // LDEBUG=TRUE;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] (mode_load==STRUCTURE_MODE_USE)" << endl;
      string aus_title=str.title;
      str=aflowlib::PrototypeLibraries(*voss,label_use,PARAMS->parameters,PARAMS->vatomX,PARAMS->vvolumeX,PARAMS->volume_in,PARAMS->mode);
      // if(LDEBUG) {oss << str << endl;exit(0);}
      str.title=aus_title+" (use of "+label_use+")";
      if(LDEBUG) *voss << str << endl;
      str.FixLattices();
      natoms=str.atoms.size();
      nspecies=str.species.size();
      natomsX.clear();
      for(uint i=0;i<str.species.size();i++) natomsX.push_back(0);
      for(uint i=0;i<str.atoms.size();i++) natomsX.at(str.atoms.at(i).type)++;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.num_each_type.size()=" << str.num_each_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.copmp_each_type.size()=" << str.comp_each_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species.size()=" << str.species.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp.size()=" << str.species_pp.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_type.size()=" << str.species_pp_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_version.size()=" << str.species_pp_version.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_ZVAL.size()=" << str.species_pp_ZVAL.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_pp_vLDAU.size()=" << str.species_pp_vLDAU.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_volume.size()=" << str.species_volume.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species_mass.size()=" << str.species_mass.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.num_each_type.size()=" << str.num_each_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] natoms=" << natoms << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] nspecies=" << nspecies << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] natomsX.size()=" << natomsX.size() << endl;
      if(LDEBUG) for(uint i=0;i<str.species.size();i++)
		   *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE] str.species.at(" << i << ")=" << str.species.at(i)
			 << " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
			 << " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
			 << " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
			 << " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
			 << " str.species_pp_vLDAU.at(" << i << ").size()=" << str.species_pp_vLDAU.at(i).size()
			 << " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
			 << " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
			 << " str.num_each_type.at(" << i << ")=" << str.num_each_type.at(i)
			 << " str.comp_each_type.at(" << i << ")=" << str.comp_each_type.at(i)
			 << " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;
      // natoms=str.atoms.size();
      // nspecies=str.species.size();
      // if(isTET) aflowlib::PrototypeFixTET(*voss,str,optionsTET);
      // return str;
    }

    // -------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_SWAP_AB) {
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP] (mode_load==STRUCTURE_SWAP_AB)" << endl;
      xstructure str_tmp;
      deque<string> vatomX_tmp;
      deque<double> vvolumeX_tmp;
      vatomX_tmp.push_back(PARAMS->vatomX.at(1));vvolumeX_tmp.push_back(PARAMS->vvolumeX.at(1));
      vatomX_tmp.push_back(PARAMS->vatomX.at(0));vvolumeX_tmp.push_back(PARAMS->vvolumeX.at(0));

      str_tmp=aflowlib::PrototypeLibraries(*voss,label_swap,PARAMS->parameters,vatomX_tmp,vvolumeX_tmp,PARAMS->volume_in,PARAMS->mode);
      // FIX THE SPECIES
      if(str_tmp.species.size() != str_tmp.num_each_type.size()) {
	*voss << "ERROR: aflow_xproto.cpp  in the SWAP_AB MODE" << endl; exit(0);
      }

      nspecies=str_tmp.num_each_type.size();
      str.scale=str_tmp.scale;
      str.lattice=str_tmp.lattice;
      str.num_each_type.clear();
      str.comp_each_type.clear();
      str.neg_scale=str_tmp.neg_scale;
      str.title+=" (swap of "+label_swap+")";
      // str.species.clear(); for(uint i=0;i<nspecies;i++) str.species.push_back(""); // create enough species
      // str.species_pp.clear(); for(uint i=0;i<nspecies;i++) str.species_pp.push_back(""); // create enough species_pp
      // str.species_pp_type.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_type.push_back(""); // create enough species_pp_type
      // str.species_pp_version.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_version.push_back(""); // create enough species_pp_version
      // str.species_pp_ZVAL.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_ZVAL.push_back(0.0); // create enough species_pp_ZVAL
      // str.species_pp_vLDAU.clear(); for(uint i=0;i<nspecies;i++) str.species_pp_vLDAU.push_back(deque<double>()); // create enough species_pp_vLDAU
      // str.species_volume.clear(); for(uint i=0;i<nspecies;i++) str.species_volume.push_back(0.0); // create enough species_volume
      // str.species_mass.clear(); for(uint i=0;i<nspecies;i++) str.species_mass.push_back(0.0); // create enough species_mass

      if(nspecies>1)
	for(iat=0;iat<str_tmp.atoms.size();iat++)
	  if(str_tmp.atoms.at(iat).type==1) str.AddAtom(str_tmp.atoms.at(iat));
      for(iat=0;iat<str_tmp.atoms.size();iat++)
	if(str_tmp.atoms.at(iat).type==0) str.AddAtom(str_tmp.atoms.at(iat));
      if(nspecies>2)
	for(iat=0;iat<str_tmp.atoms.size();iat++)
	  if(str_tmp.atoms.at(iat).type==2) str.AddAtom(str_tmp.atoms.at(iat));

      if(nspecies==1) {
	str.num_each_type.clear();str.num_each_type.push_back(str_tmp.num_each_type.at(0));
	str.comp_each_type.clear();str.comp_each_type.push_back(str_tmp.comp_each_type.at(0));
	str.species.clear();str.species.push_back(str_tmp.species.at(0));
	str.species_pp.clear();str.species_pp.push_back(str_tmp.species_pp.at(0));
	str.species_pp_type.clear();str.species_pp_type.push_back("");
	str.species_pp_version.clear();str.species_pp_version.push_back("");
	str.species_pp_ZVAL.clear();str.species_pp_ZVAL.push_back(0.0);
	str.species_pp_vLDAU.clear();str.species_pp_vLDAU.push_back(deque<double>());
	str.species_volume.clear();str.species_volume.push_back(str_tmp.species_volume.at(0));
	str.species_mass.clear();str.species_mass.push_back(str_tmp.species_mass.at(0));
      }

      if(nspecies>1) {
	str.num_each_type.clear();str.num_each_type.push_back(str_tmp.num_each_type.at(1));str.num_each_type.push_back(str_tmp.num_each_type.at(0));
	str.comp_each_type.clear();str.comp_each_type.push_back(str_tmp.comp_each_type.at(1));str.comp_each_type.push_back(str_tmp.comp_each_type.at(0));
	str.species.clear();str.species.push_back(str_tmp.species.at(1));str.species.push_back(str_tmp.species.at(0));
	str.species_pp.clear();str.species_pp.push_back(str_tmp.species_pp.at(1));str.species_pp.push_back(str_tmp.species_pp.at(0));
	str.species_pp_type.clear();str.species_pp_type.push_back("");str.species_pp_type.push_back("");
	str.species_pp_version.clear();str.species_pp_version.push_back("");str.species_pp_version.push_back("");
	str.species_pp_ZVAL.clear();str.species_pp_ZVAL.push_back(0.0);str.species_pp_ZVAL.push_back(0.0);
	str.species_pp_vLDAU.clear();str.species_pp_vLDAU.push_back(deque<double>());str.species_pp_vLDAU.push_back(deque<double>());
	str.species_volume.clear();str.species_volume.push_back(str_tmp.species_volume.at(1));str.species_volume.push_back(str_tmp.species_volume.at(0));
	str.species_mass.clear();str.species_mass.push_back(str_tmp.species_mass.at(1));str.species_mass.push_back(str_tmp.species_mass.at(0));

	if(nspecies>2) {
	  for(uint isp=2;isp<nspecies;isp++) {
	    str.num_each_type.push_back(str_tmp.num_each_type.at(isp));
	    str.comp_each_type.push_back(str_tmp.comp_each_type.at(isp));
	    str.species.push_back(str_tmp.species.at(isp));
	    str.species_pp.push_back(str_tmp.species_pp.at(isp));
	    str.species_pp_type.push_back("");
	    str.species_pp_version.push_back("");
	    str.species_pp_ZVAL.push_back(0.0);
	    str.species_pp_vLDAU.push_back(deque<double>());
	    str.species_volume.push_back(str_tmp.species_volume.at(isp));
	    str.species_mass.push_back(str_tmp.species_mass.at(isp));
	  }
	}
      }

      str.species_pp=str.species;
      // patch for species size to avoid pushing back inside
      if(str.num_each_type.size()<str.species.size()) {
	while(str.num_each_type.size()<str.species.size()) {
	  oss << "DELETING str.species.at(str.species.size()-1)=" << str.species.at(str.species.size()-1) << endl;
	  str.species.pop_back();
	  str.species_pp.pop_back();
	  str.species_pp_type.pop_back();
	  str.species_pp_version.pop_back();
	  str.species_pp_ZVAL.pop_back();
	  str.species_pp_vLDAU.pop_back();
	  str.species_volume.pop_back();
	  str.species_mass.pop_back();
	}
      }
      if(isTET) aflowlib::PrototypeFixTET(*voss,str,optionsTET);
      return str;
    }

    // -------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_SWAP_XY) {
      // LDEBUG=TRUE;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_SWAP_XY)" << endl;
      // xstructure str;
      deque<string> vatomX_tmp(PARAMS->vatomX);
      deque<double> vvolumeX_tmp(PARAMS->vvolumeX);
      aurostd::StringSubst(index_swap,"SWAP(","");aurostd::StringSubst(index_swap,")",""); // SAFETY
      aurostd::string2tokens(index_swap,tokens,",");
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 2 " << endl;
      if(tokens.size()<2) {oss << "ERROR: xproto.cpp tokens.size()>=2 in STRUCTURE_SWAP_XY" << endl;exit(0);}
      vector<uint> ispecies(2);
      for(uint i=0;i<tokens.size();i++) {
	if(tokens.at(i)=="A" || tokens.at(i)=="a" || tokens.at(i)=="0") ispecies.at(i)=0;
	if(tokens.at(i)=="B" || tokens.at(i)=="b" || tokens.at(i)=="1") ispecies.at(i)=1;
	if(tokens.at(i)=="C" || tokens.at(i)=="c" || tokens.at(i)=="2") ispecies.at(i)=2;
	if(tokens.at(i)=="D" || tokens.at(i)=="d" || tokens.at(i)=="3") ispecies.at(i)=3;
	if(tokens.at(i)=="E" || tokens.at(i)=="e" || tokens.at(i)=="4") ispecies.at(i)=4;
	if(tokens.at(i)=="F" || tokens.at(i)=="f" || tokens.at(i)=="5") ispecies.at(i)=5;
	if(tokens.at(i)=="G" || tokens.at(i)=="g" || tokens.at(i)=="6") ispecies.at(i)=6;
	if(tokens.at(i)=="H" || tokens.at(i)=="h" || tokens.at(i)=="7") ispecies.at(i)=7;
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 3 " << endl;
      if(LDEBUG) *voss << "ispecies.size()=" << ispecies.size() << endl;

      if(ispecies.at(0)!=ispecies.at(1)) {
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 4 " << endl;
	// SWAPOUT
	xstructure str_tmp;
	aurostd::swap(vatomX_tmp.at(ispecies.at(0)),vatomX_tmp.at(ispecies.at(1)));
	aurostd::swap(vvolumeX_tmp.at(ispecies.at(0)),vvolumeX_tmp.at(ispecies.at(1)));
	str_tmp=aflowlib::PrototypeLibraries(*voss,label_swap,PARAMS->parameters,vatomX_tmp,vvolumeX_tmp,PARAMS->volume_in,PARAMS->mode);
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] 5 " << endl;

	// FIX THE SPECIES
	if(str_tmp.species.size() != str_tmp.num_each_type.size()) {
	  oss << "ERROR: aflow_xproto.cpp  in the SWAP_XY MODE" << endl; exit(0);
	}

	nspecies=str_tmp.num_each_type.size();
	str.scale=str_tmp.scale;
	str.lattice=str_tmp.lattice;
	str.num_each_type.clear();
	str.comp_each_type.clear();
	str.neg_scale=str_tmp.neg_scale;
	str.title+=" (swap("+index_swap+") of "+label_swap+")";
	// str.num_each_type=str_tmp.num_each_type; // copy deque num_each_type
	// str.comp_each_type=str_tmp.comp_each_type; // copy deque comp_each_type
	// str.species=str_tmp.species; // copy deque species
	// str.species_pp=str_tmp.species_pp; // copy deque species_pp
	// str.species_pp_type=str_tmp.species_pp_type; // copy deque species_pp_type
	// str.species_pp_version=str_tmp.species_pp_version; // copy deque species_pp_version
	// str.species_pp_ZVAL=str_tmp.species_pp_ZVAL; // copy deque species_pp_ZVAL
	// str.species_pp_vLDAU=str_tmp.species_pp_vLDAU.size(); // copy deque species_pp_vLDAU
	// str.species_volume=str_tmp.species_volume; // copy deque species_volume
	// str.species_mass=str_tmp.species_mass; // copy deque species_mass
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.num_each_type=" << str_tmp.num_each_type << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.comp_each_type=" << str_tmp.comp_each_type << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species=" << str_tmp.species << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp=" << str_tmp.species_pp << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_type=" << str_tmp.species_pp_type << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_version=" << str_tmp.species_pp_version << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_ZVAL=" << str_tmp.species_pp_ZVAL << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_pp_vLDAU=" << str_tmp.species_pp_vLDAU.size() << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_volume=" << str_tmp.species_volume << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str_tmp.species_mass=" << str_tmp.species_mass << endl;
	// nspecies=str.num_each_type.size();
	// *voss << nspecies << endl;
	// *voss << str_tmp.atoms.at(0).type << str_tmp.atoms.at(1).type << str_tmp.atoms.at(2).type << endl;
	// *voss << PARAMS->vatomX.at(ispecies.at(0)) <<  PARAMS->vatomX.at(ispecies.at(1)) << endl;
	// *voss << PARAMS->vatomX_tmp.at(ispecies.at(0)) <<  PARAMS->vatomX_tmp.at(ispecies.at(1)) << endl;
	if(nspecies>1) {
	  for(uint isp=0;isp<PARAMS->vatomX.size();isp++)
	    for(uint iat=0;iat<str_tmp.atoms.size();iat++)
	      if(str_tmp.atoms.at(iat).name==PARAMS->vatomX.at(isp)) {
		if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] plugging=" <<  str_tmp.atoms.at(iat).name << endl;
		str.AddAtom(str_tmp.atoms.at(iat));
	      }
	}
	// aurostd::swap(str.num_each_type.at(ispecies.at(0)),str.num_each_type.at(ispecies.at(1)));  // taken care by AddAtom
	// aurostd::swap(str.comp_each_type.at(ispecies.at(0)),str.comp_each_type.at(ispecies.at(1)));  // taken care by AddAtom
	// aurostd::swap(str.species.at(ispecies.at(0)),str.species.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp.at(ispecies.at(0)),str.species_pp.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_type.at(ispecies.at(0)),str.species_pp_type.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_version.at(ispecies.at(0)),str.species_pp_version.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_ZVAL.at(ispecies.at(0)),str.species_pp_ZVAL.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_pp_vLDAU.at(ispecies.at(0)),str.species_pp_vLDAU.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_volume.at(ispecies.at(0)),str.species_volume.at(ispecies.at(1))); // taken care by AddAtom
	// aurostd::swap(str.species_mass.at(ispecies.at(0)),str.species_mass.at(ispecies.at(1))); // taken care by AddAtom

	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.num_each_type=" << str.num_each_type << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.comp_each_type=" << str.comp_each_type << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species=" << str.species << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp=" << str.species_pp << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_type=" << str.species_pp_type << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_version=" << str.species_pp_version << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_ZVAL=" << str.species_pp_ZVAL << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_pp_vLDAU=" << str.species_pp_vLDAU.size() << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_volume=" << str.species_volume << endl;
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] str.species_mass=" << str.species_mass << endl;

      }

      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [SWAP(,)] X " << endl;
      // str.SpeciesPutAlphabetic();
      if(isTET) aflowlib::PrototypeFixTET(*voss,str,optionsTET);
      return str;
    }
    // -------------------------------------------------------------------
    // LATTICE
    if(mode_load==STRUCTURE_MODE_RAW) {
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_MODE_RAW)" << endl;
      klattice=kmode+1;
      katoms=klattice+3;
      if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) klattice=" << klattice << endl;
      if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) katoms=" << katoms << endl;
      // GET LATTICE
      for(uint i=1;i<=3;i++) {
	aus.clear();aus.str(aurostd::RemoveRounding(database[klattice+i-1]));
	for(uint j=1;j<=3;j++)
	  aus >> str.lattice(i,j);
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.lattice=" << endl << str.lattice << endl;
      str.FixLattices(); // abc alpha beta gamma calculations are inside
    }
    if(mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_MODE_WYC)" << endl;
      klattice=kmode+1;
      katoms=klattice+1;
      if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) klattice=" << klattice << endl;
      if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) katoms=" << katoms << endl;
      aus.clear();aus.str(aurostd::RemoveRounding(database[klattice])+" 0 "); // put the phony origin choice if asked.
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) " << aus.str() << endl;
      aus >> str.a >> str.b >> str.c;
      aus >> str.alpha >> str.beta >> str.gamma;
      aus >> str.spacegroupnumber;
      aus >> str.spacegroupnumberoption;
      if(SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) {
	//     cerr << "SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE" << endl;
	if(PARAMS->flip_option==FALSE) {
	  // cerr << "SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE && PARAMS->flip_option==FALSE" << endl;
	  // cerr << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
	  // cerr << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;
	  // some ICSD hacks
	  if(str.spacegroupnumber==3 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==4 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==5 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==6 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==8 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==9 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==10 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==11 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==12 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==13 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==14 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==15 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==48 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==50 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  // if(str.spacegroupnumber==59 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK (sometimes it is 2, though)
	  if(str.spacegroupnumber==62 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=2;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==68 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==70 && str.spacegroupnumberoption==0) {str.spacegroupnumberoption=1;}  // ICSD seems std OK
	  if(str.spacegroupnumber==85 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==86 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==88 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==125 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==126 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==129 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=2;  // ICSD seems std OK
	  if(str.spacegroupnumber==130 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==133 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==134 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==137 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==138 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==141 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==142 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==146 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  // if(str.spacegroupnumber==148 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==148 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=2;  // ICSD seems std OK
	  if(str.spacegroupnumber==155 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK
	  if(str.spacegroupnumber==160 && str.spacegroupnumberoption==0) {
	    if(isequal(str.alpha,str.beta,0.1) && isequal(str.beta,str.gamma,0.1)) str.spacegroupnumberoption=1; // for ICSD in RHL format try 1
	    if(isequal(str.alpha,str.beta,0.1) && isequal(str.beta,str.gamma,0.1)) str.spacegroupnumberoption=1; // for ICSD in RHL format try 1
	    if(isequal(str.alpha,120.0,0.1) || isequal(str.beta,120.0,0.1) || isequal(str.gamma,120.0,0.1)) str.spacegroupnumberoption=2; // for ICSD in HEX format try 2
	  } //912 267 9998
	  // if(str.spacegroupnumber==161 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=2;  // ICSD seems std OK
	  if(str.spacegroupnumber==166 && str.spacegroupnumberoption==0) { //str.spacegroupnumberoption=2;  // ICSD seems std OK
	    // for 166 you can switch to 2 if there is a 120 angle, otherwise you should keep 1
	    if(isequal(str.alpha,str.beta,0.1) && isequal(str.beta,str.gamma,0.1)) str.spacegroupnumberoption=1; // for ICSD in RHL format try 1
	    if(isequal(str.alpha,120.0,0.1) || isequal(str.beta,120.0,0.1) || isequal(str.gamma,120.0,0.1)) str.spacegroupnumberoption=2; // for ICSD in HEX format try 2
	  }
	  // if(str.spacegroupnumber==167 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=2;  // ICSD seems std OK
	  if(str.spacegroupnumber==201 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==203 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==222 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==224 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  // if(str.spacegroupnumber==227 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK (sometimes it is 2, though)
	  if(str.spacegroupnumber==227 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  if(str.spacegroupnumber==228 && str.spacegroupnumberoption==0) str.spacegroupnumberoption=1;  // ICSD seems std OK CHECK
	  // cerr << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
	  // cerr << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;
	}
      }
      if(PARAMS->flip_option==TRUE && SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) str.spacegroupnumberoption=3-str.spacegroupnumberoption;
      if(SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE && str.spacegroupnumberoption!=1 && str.spacegroupnumberoption!=2) str.spacegroupnumberoption=2;
      //  if(SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) cerr << "SPACEGROUP.OPTION=" << str.spacegroupnumber << "." << str.spacegroupnumberoption << endl;
      if(SpaceGroupOptionRequired(str.spacegroupnumber)==FALSE) str.spacegroupnumberoption=0;
      if(str.spacegroupnumberoption==3) str.spacegroupnumberoption=2;
      //   cerr << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
      // cerr << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;

      str.spacegrouplabel=GetSpaceGroupLabel(str.spacegroupnumber);
      str.spacegroup=GetSpaceGroupName(str.spacegroupnumber,str.directory);//DX 20180526 - add directory
      str.lattice=GetClat(str.a,str.b,str.c,str.alpha,str.beta,str.gamma);
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) " << str.a << " " <<  str.b << " " <<  str.c << " " <<  str.alpha << " " <<  str.beta << " " <<  str.gamma
		       << " " <<  str.spacegroupnumber << " " <<  str.spacegroupnumberoption << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.lattice=" << str.lattice << endl;
      str.FixLattices();
      if(mode_load==STRUCTURE_MODE_ICSD) {
	// load species
	aurostd::string2tokens(database[kstructure],tokens);
	KBIN::VASP_SplitAlloySpecies(aurostd::RemoveSubStringFirst(tokens[1],"C-"),speciesX);
      }
    }
    //cerr << str << endl;exit(0);
    if(mode_load==STRUCTURE_MODE_ABC) {
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (mode_load==STRUCTURE_MODE_ABC)" << endl;
      klattice=kmode+1;
      katoms=klattice+1;
      if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) klattice=" << klattice << endl;
      if(LDEBUG) *voss <<"DEBUG: (aflowlib::PrototypeLibraries) katoms=" << katoms << endl;
      aus.clear();aus.str(aurostd::RemoveRounding(database[klattice]));
      aurostd::string2tokens(database[klattice],tokens);
      // *voss << tokens.size() << endl;
      aus >> str.a >> str.b >> str.c;
      aus >> str.alpha >> str.beta >> str.gamma;
      // if(tokens.size()>=7) aus >> str.spacegroupnumber;
      if(tokens.size()>=8)  aus >> str.spacegroupnumberoption;
      str.lattice=GetClat(str.a,str.b,str.c,str.alpha,str.beta,str.gamma);
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.lattice=" << str.lattice << endl;
      str.FixLattices();
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupoption=[" << str.spacegroupoption << "]" << endl;}
    // -------------------------------------------------------------------
    if(mode_load==STRUCTURE_MODE_RAW || mode_load==STRUCTURE_MODE_ABC || mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      // LDEBUG=TRUE;   // FIX ATOMS
      uint nspecies_old=nspecies;
      natoms=0;nspecies=0;
      // for(uint i=0;i<speciesX.size();i++) *voss << speciesX.at(i) << endl;
      // for(uint i=0;i<PARAMS->vvolumeX.size();i++) *voss << PARAMS->vvolumeX.at(i) << endl;

      // check species and volumes from library
      vector<string> lvatomX;
      for(iat=0;iat<(uint) kstop-katoms;iat++) {
	aurostd::string2tokens(database[iat+katoms],tokens);
	bool found=FALSE;
	for(uint i=0;i<lvatomX.size()&&!found;i++) if(tokens[3]==lvatomX.at(i)) found=TRUE;
	if(!found) lvatomX.push_back(tokens[3]);
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.size()=" << lvatomX.size() << endl;
      if(LDEBUG) {for(uint i=0;i<lvatomX.size();i++) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.at(" << i << ")=" << lvatomX.at(i) << endl;}
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] speciesX.size()=" << speciesX.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.size()=" << str.species.size() << endl;
      if(lvatomX.size()>speciesX.size()) {
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] FIX nspecies" << endl;
	for(uint i=speciesX.size();i<lvatomX.size();i++) {
	  PARAMS->vatomX.push_back(lvatomX.at(i));
	  speciesX.push_back(lvatomX.at(i));
	  PARAMS->vvolumeX.push_back(1.0);
	  natomsX.push_back(0);
	}
	nspecies_old=PARAMS->vatomX.size();
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.size()=" << lvatomX.size() << endl;
      if(LDEBUG) {for(uint i=0;i<lvatomX.size();i++) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] lvatomX.at(" << i << ")=" << lvatomX.at(i) << endl;}
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] speciesX.size()=" << speciesX.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.size()=" << str.species.size() << endl;

      // fixed species
      for(iat=0;iat<(uint) kstop-katoms;iat++) {
	aurostd::string2tokens(database[iat+katoms],tokens);
	// *voss << iat << " " << tokens[3] << endl;
	natoms++;
	for(uint i=0;i<speciesX.size();i++)
	  if(tokens[3]==speciesX.at(i)) natomsX.at(i)++;
      }
      for(uint i=0;i<speciesX.size();i++)
	if(natomsX.at(i)>0) {str.num_each_type.push_back(natomsX.at(i));str.comp_each_type.push_back((double) natomsX.at(i));nspecies++;}
      if(mode_load==STRUCTURE_MODE_ICSD) {
	if(nspecies!=nspecies_old) {
	  oss << "ERROR (aflow_xproto.cpp): PARAMS->label=" << PARAMS->label;
	  oss << "  nspecies!=nspecies_old (" << nspecies << "," << nspecies_old << ") *************" << endl;
	  exit(0);
	}
      }
      // now I have nspecies
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] NSPECIES_OLD=" << nspecies_old << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] NSPECIES=" << nspecies << endl;

      str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear();
      for(uint i=0;i<nspecies;i++) {
	str.species.push_back(PARAMS->vatomX.at(i));    // only the required ones
	str.species_pp.push_back(PARAMS->vatomX.at(i)); // only the required ones
	str.species_pp_type.push_back(""); // only the required ones
	str.species_pp_version.push_back(""); // only the required ones
	str.species_pp_ZVAL.push_back(0.0); // only the required ones
	str.species_pp_vLDAU.push_back(deque<double>()); // only the required ones
	str.species_volume.push_back(0.0);     // only the required ones
	str.species_mass.push_back(0.0);     // only the required ones
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] natoms=" << natoms << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] nspecies=" << nspecies << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.size()=" << str.species.size() << endl;
      if(LDEBUG) for(uint i=0;i<nspecies;i++)
		   *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] str.species.at(" << i << ")=" << str.species.at(i)
			 << " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
			 << " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
			 << " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
			 << " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
			 << " str.species_pp_vLDAU.at(" << i << ")=" << str.species_pp_vLDAU.at(i).size()
			 << " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
			 << " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
			 << " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;

      // Plug ATOMS
      volume=0.0;
      // perform atoms A,B,C,D,E... ---------------------------------

      for(uint ispecies=1;ispecies<=nspecies;ispecies++) { // pick atom X
	if(speciesX.size()>=ispecies) {
	  for(iat=0;iat<natoms;iat++) {
	    aurostd::string2tokens(database[iat+katoms],tokens);
	    if(tokens[3]==speciesX.at(ispecies-1)) { // atoms X
	      _atom atom;
	      atom.CleanName();atom.CleanSpin();
	      aus.clear();aus.str(database[iat+katoms]);
	      aus >> atom.fpos(1) >> atom.fpos(2) >> atom.fpos(3);
	      atom.cpos=F2C(str.lattice,atom.fpos);
	      atom.name=str.species.at(ispecies-1);atom.type=0;
	      atom.CleanName(); // to fix the real name
	      for(uint icheck=0;icheck<ispecies-1;icheck++)
		if(str.species.at(icheck)!="") atom.type++; // X is 0,1,2,3 to ispecies-1
	      atom.name_is_given=TRUE;str.atoms.push_back(atom);
	    }
	  }
	}
      }
      // ------ --------------------------------------------
      for(uint i=0;i<nspecies;i++)
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [RAWI] PARAMS->vvolumeX.at(" << i << ")=" << PARAMS->vvolumeX.at(i) << endl;
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupoption=[" << str.spacegroupoption << "]" << endl;}
    // -------------------------------------------------------------------
    if(mode_remove==STRUCTURE_MODE_REMOVE) {
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] BEFORE REMOVE ********************* " << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natoms=" << natoms << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] nspecies=" << nspecies << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.size()=" << str.species.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp.size()=" << str.species_pp.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_type.size()=" << str.species_pp_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_version.size()=" << str.species_pp_version.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_ZVAL.size()=" << str.species_pp_ZVAL.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_vLDAU.size()=" << str.species_pp_vLDAU.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natomsX.size()=" << natomsX.size() << endl;
      if(LDEBUG) for(uint i=0;i<nspecies;i++)
		   *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.at(" << i << ")=" << str.species.at(i)
			 << " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
			 << " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
			 << " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
			 << " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
			 << " str.species_pp_vLDAU.at(" << i << ")=" << str.species_pp_vLDAU.at(i).size()
			 << " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
			 << " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
			 << " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;

      tokens.clear();
      aurostd::string2tokens(label_remove,tokens,",");
      vector<uint> ratoms;
      for(iat=0;iat<tokens.size();iat++)
	ratoms.push_back(aurostd::string2utype<uint>(tokens.at(iat)));
      std::sort(ratoms.begin(), ratoms.end());
      for(int iat=ratoms.size()-1;iat>=0;iat--)
	str.RemoveAtom(ratoms.at(iat));
      // cout << ratoms.at(iat) << endl;
      // for(uint j=1;j<tokens.size();j++) str.prototype+=tokens[j]+" ";

      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] AFTER REMOVE ********************* " << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natoms=" << natoms << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] nspecies=" << nspecies << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.size()=" << str.species.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp.size()=" << str.species_pp.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_type.size()=" << str.species_pp_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_version.size()=" << str.species_pp_version.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_ZVAL.size()=" << str.species_pp_ZVAL.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species_pp_vLDAU.size()=" << str.species_pp_vLDAU.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] natomsX.size()=" << natomsX.size() << endl;
      if(LDEBUG) for(uint i=0;i<nspecies;i++)
		   *voss << "DEBUG: (aflowlib::PrototypeLibraries) [3] str.species.at(" << i << ")=" << str.species.at(i)
			 << " str.species_pp.at(" << i << ")=" << str.species_pp.at(i)
			 << " str.species_pp_type.at(" << i << ")=" << str.species_pp_type.at(i)
			 << " str.species_pp_version.at(" << i << ")=" << str.species_pp_version.at(i)
			 << " str.species_pp_ZVAL.at(" << i << ")=" << str.species_pp_ZVAL.at(i)
			 << " str.species_pp_vLDAU.at(" << i << ")=" << str.species_pp_vLDAU.at(i).size()
			 << " str.species_volume.at(" << i << ")=" << str.species_volume.at(i)
			 << " str.species_mass.at(" << i << ")=" << str.species_mass.at(i)
			 << " natomsX.at(" << i << ")=" << natomsX.at(i) << endl;
      // exit(0);
    }

    // --------------------------------------------------------------------------------------------------------------------------------------
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupoption=[" << str.spacegroupoption << "]" << endl;}

    if(mode_load==STRUCTURE_MODE_USE && mode_species==STRUCTURE_MODE_SPECIES) {
      //   LDEBUG=TRUE;

      // this is USE_SPECIES
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES] (mode_load==STRUCTURE_MODE_USE) (mode_species==STRUCTURE_MODE_SPECIES) " << endl;
      string aus_title=str.title;
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vatomX="; for(uint i=0;i<PARAMS->vatomX.size();i++) *voss << PARAMS->vatomX.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vvolumeX.size()=" << PARAMS->vvolumeX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vvolumeX="; for(uint i=0;i<PARAMS->vvolumeX.size();i++) *voss << PARAMS->vvolumeX.at(i) << " "; *voss << endl;}

      str=aflowlib::PrototypeLibraries(*voss,label_use,PARAMS->parameters,PARAMS->vatomX,PARAMS->vvolumeX,PARAMS->volume_in,PARAMS->mode);
      //  if(LDEBUG) {oss << str << endl;exit(0);}
      str.title=aus_title+" (use of "+label_use+")";
      if(LDEBUG) *voss << str << endl;
      str.FixLattices();
      natoms=str.atoms.size();
      nspecies=str.species.size();
      natomsX.clear();
      for(uint i=0;i<str.species.size();i++) natomsX.push_back(0);
      for(uint i=0;i<str.atoms.size();i++) natomsX.at(str.atoms.at(i).type)++;

      // this is SPECIES

      // *voss << str.scale << endl; exit(0);

      tokens.clear();
      aurostd::string2tokens(label_species,tokens,",");
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  label_species=" << label_species << endl;}
      if(tokens.size()!=str.num_each_type.size()) {oss << "tokens.size()!=str.num_each_type.size()" << endl;exit(0);}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  tokens.size()=" << tokens.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.num_each_type.size()=" << str.num_each_type.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.comp_each_type.size()=" << str.comp_each_type.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  speciesX.size()=" << speciesX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species.size()=" << str.species.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_volume.size()=" << str.species_volume.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_mass.size()=" << str.species_mass.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vvolumeX.size()=" << PARAMS->vvolumeX.size() << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl;}
      if(LDEBUG) {for(uint i=0;i<PARAMS->vvolumeX.size();i++) *voss << PARAMS->vvolumeX.at(i) << " "; *voss << endl;}
      deque<double> vol(PARAMS->vvolumeX);
      for(uint i=0;i<str.num_each_type.size();i++) {
	uint j=(int) tokens.at(i)[0]-(int)'A';
	tokens.at(i)=str.species.at(j);
	vol.at(i)=PARAMS->vvolumeX.at(j);
      }
      PARAMS->vvolumeX=vol;
      for(uint i=0;i<PARAMS->vvolumeX.size();i++) str.species_volume.at(i)=vol.at(i);
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  volume="; for(uint i=0;i<PARAMS->vvolumeX.size();i++) *voss << PARAMS->vvolumeX.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  tokens="; for(uint i=0;i<tokens.size();i++) *voss << tokens.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species="; for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_volume="; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_mass="; for(uint i=0;i<str.species_mass.size();i++) *voss << str.species_mass.at(i) << " "; *voss << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vatomX="; for(uint i=0;i<PARAMS->vatomX.size();i++) *voss << PARAMS->vatomX.at(i) << " "; *voss << endl;}
      uint iat=0;
      for(uint i=0;i<str.num_each_type.size();i++) {
	// *voss << str.species.at(i) << endl;
	str.species.at(i)=tokens.at(i);
	for(uint j=0;j<(uint) str.num_each_type.at(i);j++) {
	  str.atoms.at(iat).type=i;
	  str.atoms.at(iat).name=str.species.at(i);
	  str.atoms.at(iat).CleanName(); // cool little name
	  iat++;
	}
      }
      // for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;
      // JUNKAI
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      str.SpeciesPutAlphabetic();
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      //   exit(0);
      // for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;
      for(uint i=0;i<str.species_volume.size();i++) PARAMS->vvolumeX.at(i)=str.species_volume.at(i);
      // make species together
      // if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
      for(uint i=1;i<str.num_each_type.size();i++) {
	if(str.species.at(i)==str.species.at(i-1)) {
	  str.species.erase(str.species.begin()+i);
	  str.num_each_type.at(i-1)+=str.num_each_type.at(i);str.num_each_type.erase(str.num_each_type.begin()+i);
	  str.comp_each_type.at(i-1)+=str.comp_each_type.at(i);str.comp_each_type.erase(str.comp_each_type.begin()+i);
	  PARAMS->vvolumeX.erase(PARAMS->vvolumeX.begin()+i);
	  str.species_volume.erase(str.species_volume.begin()+i);
	  str.species_mass.erase(str.species_mass.begin()+i);
	  // PARAMS->vatomX.erase(PARAMS->vatomX.begin()+i);
	  tokens.erase(tokens.begin()+i);
	  i=0;
	}
      }
      // fix types
      uint k=0;
      for(uint i=0;i<str.num_each_type.size();i++)
	for(uint j=0;j<(uint) str.num_each_type.at(i);j++)
	  str.atoms.at(k++).type=i; // done
      nspecies=str.num_each_type.size();
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  [4] nspecies=" << nspecies << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  species.size()=" << str.species.size() << endl;
      if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  tokens.size()=" << tokens.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.num_each_type.size()=" << str.num_each_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.comp_each_type.size()=" << str.comp_each_type.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  speciesX.size()=" << speciesX.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species.size()=" << str.species.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_volume.size()=" << str.species_volume.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  str.species_mass.size()=" << str.species_mass.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vvolumeX.size()=" << PARAMS->vvolumeX.size() << endl;
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [USE_SPECIES]  PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl;
      if(LDEBUG) {for(uint i=0;i<PARAMS->vvolumeX.size();i++) *voss << PARAMS->vvolumeX.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<tokens.size();i++) *voss << tokens.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.species_mass.size();i++) *voss << str.species_mass.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<PARAMS->vatomX.size();i++) *voss << PARAMS->vatomX.at(i) << " "; *voss << endl;}
      if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
      if(LDEBUG) *voss << str << endl;
      //   exit(0);
      //   str.SpeciesPutAlphabetic();
    }
    // -------------------------------------------------------------------
    // FIX ALL ATOMS
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) FIX ALL ATOMS (if necessary)" << endl;
    bool DEBUG_ICSD=TRUE;
    xstructure strw=str;
    if(mode_load==STRUCTURE_MODE_ICSD) DEBUG_ICSD=TRUE;
    if(mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) WYC/ICSD str.species.size()=" << str.species.size() << endl;
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) WYC/ICSD strw.species.size()=" << strw.species.size() << endl;
      if(LDEBUG) {for(uint i=0;i<strw.species.size();i++) *voss << strw.species.at(i) << " "; *voss << endl;}
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->label = " << PARAMS->label << " calculating WyckoffPOSITIONS()" << endl;

      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumber=[" << str.spacegroupnumber << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.spacegroupnumberoption=[" << str.spacegroupnumberoption << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) strw.spacegroupnumber=[" << strw.spacegroupnumber << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) strw.spacegroupnumberoption=[" << strw.spacegroupnumberoption << "]" << endl;}

      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) BEFORE WyckoffPOSITIONS" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.atoms.size()=[" << str.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) strw.atoms.size()=[" << strw.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vatomX.size()=[" << PARAMS->vatomX.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vvolumeX.size()=[" << PARAMS->vvolumeX.size()  << "]" << endl;}
      str=WyckoffPOSITIONS(strw.spacegroupnumber,strw.spacegroupnumberoption,strw);
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) AFTER WyckoffPOSITIONS" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) str.atoms.size()=[" << str.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) strw.atoms.size()=[" << strw.atoms.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vatomX.size()=[" << PARAMS->vatomX.size()  << "]" << endl;}
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) PARAMS->vvolumeX.size()=[" << PARAMS->vvolumeX.size()  << "]" << endl;}

      //   cerr << str << endl;exit(0);
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) WYC/ICSD str.species.size()=" << str.species.size() << endl;
      if(LDEBUG) {for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;}
      string sgstring=(SpaceGroupOptionRequired(str.spacegroupnumber)?aurostd::utype2string(str.spacegroupnumber)+"."+aurostd::utype2string(str.spacegroupnumberoption):aurostd::utype2string(str.spacegroupnumber));
      if(mode_load==STRUCTURE_MODE_WYC) {
	label_library="";
	PARAMS->vatomX.at(0)+PARAMS->vatomX.at(1)+"/"+label_library;
      }


      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(6): " << str.title << endl;

      if(mode_load==STRUCTURE_MODE_ICSD) {
	str.title=label_library+" "+"#"+sgstring+" - "+"("+PARAMS->label+")";  // SG# AS Wahyu wants
	str.title+=" - "+str.prototype+" "+_ICSD_STRING_;                                                // SG# as Wahyu wants
	str.title=str.title+" (WICKOFF "+sgstring+" "+str.spacegrouplabel+")";
      }
      if(mode_load==STRUCTURE_MODE_WYC) {
	str.title=str.title+" (WICKOFF "+sgstring+" "+str.spacegrouplabel+")";
      }
      if(mode_load==STRUCTURE_MODE_ICSD) {
	for(uint i=0;i<str.atoms.size();i++) {
	  str.atoms.at(i).name=speciesX.at(str.atoms.at(i).type);
	  str.atoms.at(i).CleanName();
	}
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) TITLE(7): " << str.title << endl;
      //  exit(0);

      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) WYC/ICSD str.species.size()=" << str.species.size() << endl;
      if(0 && DEBUG_ICSD) {
	*voss << endl << "calculating space group" << endl;
	int spacegroup_test=str.spacegroupnumber;
	str.platon2sg();
	*voss << "  sgs(" << spacegroup_test << "," << str.spacegroupnumber << ")" << endl;
	if(spacegroup_test!=str.spacegroupnumber) {
	  oss << "ERROR (aflow_xproto.cpp): PARAMS->label=" << PARAMS->label;
	  oss << " " << spacegroup_test << " " << str.spacegroupnumber << " CHANGE space group option in compound " << PARAMS->label << "   *************" << endl;
	  exit(0);
	}
      }
      if(PARAMS->mode==LIBRARY_MODE_ICSD) {
	if(NearestNeighbour(str)<_XPROTO_TOO_CLOSE_ERROR_ && PARAMS->flip_option==FALSE && SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) {
	  // *voss << "AFLOW WARNING (aflow_xproto.cpp): PARAMS->label=" << PARAMS->label << " WRONG NNdist too close =" << NearestNeighbour(str) << "   *************"  << endl;
	  *voss << "AFLOW WARNING (aflow_xproto.cpp): Too close NNdist(" << NearestNeighbour(str) << "), Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption << " not enough, try PARAMS->flip_option " << endl;
	  xstructure str_flipped;
	  str_flipped=aflowlib::PrototypeLibraries(*voss,PARAMS->label,PARAMS->parameters,PARAMS->vatomX,PARAMS->vvolumeX,PARAMS->volume_in,PARAMS->mode,TRUE);
	  str_flipped.title+=" (sg.opt flipped)";
	  str_flipped.species_pp=str_flipped.species;
	  *voss << "AFLOW WARNING (aflow_xproto.cpp): Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption <<" flipped to option=" << str_flipped.spacegroupnumberoption << endl;
	  if(isTET) aflowlib::PrototypeFixTET(*voss,str_flipped,optionsTET);
	  return str_flipped;
	}
      }
      if(PARAMS->mode==LIBRARY_MODE_ICSD) {
	if(DEBUG_ICSD) {
	  double natoms_label=0.0,natoms_icsd=0.0;
	  for(uint i=0;i<str.num_each_type.size();i++) {
	    natoms_icsd+=str.num_each_type.at(i);
	    natoms_label+=vnatomsX.at(i);
	  }
	  for(uint i=0;i<str.num_each_type.size();i++) {
	    //   cerr << str.num_each_type.at(i) << " " << natoms_icsd << " " << vnatomsX.at(i) << " " << natoms_icsd <<  endl;
	    if(abs((double) (str.num_each_type.at(i)/natoms_icsd-vnatomsX.at(i)/natoms_label))>0.01) {
	      if(PARAMS->flip_option==FALSE && SpaceGroupOptionRequired(str.spacegroupnumber)==TRUE) {
		*voss << "AFLOW WARNING (aflow_xproto.cpp): Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption << " not enough, try PARAMS->flip_option " << endl;
		xstructure str_flipped;
		str_flipped=aflowlib::PrototypeLibraries(*voss,PARAMS->label,PARAMS->parameters,PARAMS->vatomX,PARAMS->vvolumeX,PARAMS->volume_in,PARAMS->mode,TRUE);
		str_flipped.title+=" (sg.opt flipped)";
		str_flipped.species_pp=str_flipped.species;
		*voss << "AFLOW WARNING (aflow_xproto.cpp): Spacegroup=" << str.spacegroupnumber << " found option=" << str_flipped.spacegroupnumberoption << endl;
		if(isTET) aflowlib::PrototypeFixTET(*voss,str_flipped,optionsTET);
		return str_flipped;
	      }
	      // *voss << endl;
	      *voss << "ERROR (aflow_xproto.cpp): PARAMS->label=" << PARAMS->label
		    << " WRONG concentrations species=" << str.species.at(i)
		    << " conc_prototype=" << vnatomsX.at(i)/natoms_label << " "
		    << " conc_icsd=" << str.num_each_type.at(i)/natoms_icsd << " : CHECK ALSO THE SPACEGROUP" << "   *************"  << endl;
	      if(LDEBUG) *voss << str << endl;
	      exit(0);
	    }
	  }
	}
      }
    }
    // -------------------------------------------------------------------
    // FIX SCALE
    // *voss << str << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] FIX SCALE (if necessary)" << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] nspecies=" << nspecies << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.atoms.size()=" << str.atoms.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] species.size()=" << str.species.size() << endl;
    if(LDEBUG) {for(uint i=0;i<str.atoms.size();i++) *voss << str.atoms.at(i).type << " "; *voss << endl;}
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] tokens.size()=" << tokens.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.num_each_type.size()=" << str.num_each_type.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.comp_each_type.size()=" << str.comp_each_type.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] speciesX.size()=" << speciesX.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.species.size()=" << str.species.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.species_volume.size()=" << str.species_volume.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] str.species_mass.size()=" << str.species_mass.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] PARAMS->vvolumeX.size()=" << PARAMS->vvolumeX.size() << endl;
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [FIX_SCALE] PARAMS->vatomX.size()=" << PARAMS->vatomX.size() << endl;

    str.species_volume.clear();
    for(uint i=0;i<str.species.size();i++)
      str.species_volume.push_back(0.0);

    if(mode_volume==STRUCTURE_MODE_NONE) {
      if(PARAMS->volume_in>0.0) {
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume choice 1 (PARAMS->volume_in=" << PARAMS->volume_in << ")" << endl;
	volume=str.atoms.size()*PARAMS->volume_in;
      } else {
	if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume choice 2: PARAMS->vvolumeX.size()=" << PARAMS->vvolumeX.size() << endl;
	if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume choice 2: str.species_volume.size()="<<str.species_volume.size()<< ": "; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}
	volume=0.0;
	for(iat=0;iat<str.atoms.size();iat++) {
	  if(LDEBUG) *voss << iat << " " << str.atoms.at(iat).name << " " << nspecies << " " << str.atoms.size() << endl;
	  for(uint i=0;i<nspecies;i++)
	    if(str.atoms.at(iat).name==str.species.at(i)) {
	      volume+=PARAMS->vvolumeX.at(i);
	      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) volume atom (" << iat << ")=" << volume << endl;
	    }

	  for(uint i=0;i<str.species.size();i++)
	    str.species_volume.at(i)=PARAMS->vvolumeX.at(i);
	}
      }
      if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [6] volume=" << volume << endl;
      // fix volume in WYC or ABC
      str.scale=std::pow((double) (abs(volume)/det(str.lattice)),(double) 1.0/3.0);
      str.neg_scale=TRUE;
    }
    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] structure_mode_volume" << endl;
    if(mode_volume==STRUCTURE_MODE_VOLUME) {
      str.scale=1.0;
      str.ReScale(1.0);
      str.neg_scale=FALSE;
    }

    // check for nndist

    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] nndist" << endl;
    if(NearestNeighbour(str)<_XPROTO_TOO_CLOSE_ERROR_) {
      *voss << "ERROR (aflow_xproto.cpp): PARAMS->label=" << PARAMS->label << " WRONG NNdist too close =" << NearestNeighbour(str) << " Spacegroup=" << str.spacegroupnumber << " option=" << str.spacegroupnumberoption <<"   *************"  << endl;
      *voss << str << endl;
      //   str.SetCoordinates(_COORDS_CARTESIAN_);
      //  *voss << str << endl;

      exit(0);
    }

    if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1" << endl;

    double eps=-1.0,epsang=-1.0; // for lattice calculation. initialize negative
    if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL || mode_load==STRUCTURE_MODE_WYC || mode_load==STRUCTURE_MODE_ICSD) {
      //    LDEBUG=TRUE;
      /* OLD
      // get btravais_lattice_type before it is too hard
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) calling LATTICE::SpaceGroup2Lattice from xproto.cpp" << endl;
      // str.bravais_lattice_type=LATTICE::SpaceGroup2Lattice(str.spacegroupnumber);
      // str.bravais_lattice_variation_type=LATTICE::SpaceGroup2LatticeVariation(str.spacegroupnumber,str);
      // str.bravais_lattice_system=LATTICE::Lattice_System_SpaceGroup(str.spacegroupnumber);
      */
      // NEW
      //    cerr << "str.spacegroupnumber=" << str.spacegroupnumber << endl;
      //    cerr << "str.spacegroupnumberoption=" << str.spacegroupnumberoption << endl;

      // DX - START
      xstructure str_in(str),str_sp,str_sc;
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1 - Symmetry tolerance scan (DX and COREY) " << endl; }
      LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);      
 
      /* DX
         xstructure str_in(str),str_sp,str_sc;
         //    LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
         if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.1 - precise" << endl;
         LATTICE::Standard_Lattice_StructurePrecise(str_in,str_sp,str_sc);
         // LATTICE::Standard_Lattice_StructureMedium(str_in,str_sp,str_sc);
         if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
         if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.2 - try again - medium+" << endl;
         int ss=0; eps=0.005;epsang=0.05; LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang,ss,_EPS_);} // medium++
         if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
         if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.3 - try again - medium" << endl;
         int ss=0; eps=0.01;epsang=0.1; LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang,ss,_EPS_);} // medium++
         if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
         if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.4 - try again - medium-default" << endl;
         LATTICE::Standard_Lattice_StructureMedium(str_in,str_sp,str_sc);}
         if(str_sp.bravais_lattice_type.empty() || str_sp.bravais_lattice_type=="UNKNOWN") {
         if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.5 - try again - loose" << endl;
         LATTICE::Standard_Lattice_StructureDefault(str,str_sp,str_sc);}
      */
      // DX - END

      str.bravais_lattice_type=str_sp.bravais_lattice_type;
      str.bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
      str.bravais_lattice_system=str_sp.bravais_lattice_system;
      if(LDEBUG) { cerr << str.bravais_lattice_type << endl; }
      if(LDEBUG) { cerr << str_sp.bravais_lattice_type << endl; }
      if(LDEBUG) { cerr << str_sc.bravais_lattice_type << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 1.99" << endl; }
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 2" << endl; }
    if(mode_load==STRUCTURE_MODE_ABC) {
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 2a" << endl; }
      str.GetLatticeType();
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] str.bravais_lattice_type=" << str.bravais_lattice_type << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] str.bravais_lattice_variation_type=" << str.bravais_lattice_variation_type << endl; }
      if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 2b" << endl; }
      if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL) str.GetStandardPrimitive();
    }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] symmetry 3" << endl; }
    if(LDEBUG) { *voss << "DEBUG: (aflowlib::PrototypeLibraries) [7] str.species_volume.size()="<<str.species_volume.size()<< ": "; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;}

    // now if you do the PRIM you get screwed because it is very hard to reconstruct the lattice type
    if(mode_prim==STRUCTURE_MODE_PRIM) {
      // LDEBUG=TRUE;
      // cerr << "mode_prim=" << mode_prim << endl;
      // str.neg_scale=FALSE;
      // str.ReScale(1.0);
      // str.SetVolume(str.atoms.size());
      bool fixed=FALSE;
      if(mode_conventional==STRUCTURE_MODE_CONVENTIONAL || mode_load==STRUCTURE_MODE_ICSD) {
	//   cerr << "mode_conventional=" << mode_conventional << endl;
	fixed=TRUE;
	xstructure str_sp,str_sc;
	str.FixLattices();
	str_sp=str;str_sc=str; // so it creates species and stuff.
	if(eps<0.0 && epsang <0.0) {  // still undefined
	  // eps=0.1;epsang=1.0; // coarse
	  // eps=0.05;epsang=0.5;// normal
	  // eps=0.01;epsang=0.1; // medium//
	  eps=0.005;epsang=0.05; // medium++
	  // eps=0.002;epsang=0.02; // precise
	  if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): redefining  eps=" << eps << "  epsang=" << epsang << " " << endl;
	} else {
	  if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): inheriting  eps=" << eps << "  epsang=" << epsang << " " << endl;
	}
	str_sp.bravais_lattice_type="X";
        // DX - START
	// [[JUNKAI OBSOLETE]] uint step=0;
        LATTICE::Standard_Lattice_StructureDefault(str,str_sp,str_sc); // DX
        // [[JUNKAI OBSOLETE]]
	// [[JUNKAI OBSOLETE]]//  while (str_sp.bravais_lattice_type!=str.bravais_lattice_type && str_sp.bravais_lattice_type!="UNKNOWN" && step<10) {
	// [[JUNKAI OBSOLETE]]while (str_sp.bravais_lattice_type!=str.bravais_lattice_type && str_sp.bravais_lattice_type!="UNKNOWN" && step<10) {
	// [[JUNKAI OBSOLETE]]  step++;if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): step=" << step << "  eps=" << eps << "  epsang=" << epsang << " ";
	// [[JUNKAI OBSOLETE]]  // LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang); // STEFANO OLD VERSION
	// [[JUNKAI OBSOLETE]]  int ss=0;
	// [[JUNKAI OBSOLETE]]  LATTICE::Standard_Lattice_Structure(str,str_sp,str_sc,eps,epsang,ss,_EPS_);
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) *voss << str.bravais_lattice_type << endl;
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) *voss << str_sp.bravais_lattice_type << endl;
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) *voss << str_sc.bravais_lattice_type << endl;
	// [[JUNKAI OBSOLETE]]  // cerr << str << endl; cerr << str_sp << endl; cerr << str_sc << endl;exit(0);
	// [[JUNKAI OBSOLETE]]  if(LDEBUG) *voss << str_sp.bravais_lattice_variation_type << endl;
	// [[JUNKAI OBSOLETE]]  // *voss << str_sp.bravais_lattice_type << " " << eps << " " << epsang << " " << endl;
	// [[JUNKAI OBSOLETE]]  // if(str_sp.bravais_lattice_type!=str.bravais_lattice_type) {  // try different option str.bravais_lattice_type known only for ICSD
	// [[JUNKAI OBSOLETE]]  if(str_sp.bravais_lattice_type!=str.bravais_lattice_type) {
	// [[JUNKAI OBSOLETE]]    xstructure stropt=strw,stropt_sp=str,stropt_sc=str;
	// [[JUNKAI OBSOLETE]]    stropt.spacegroupnumberoption=3-stropt.spacegroupnumberoption;
	// [[JUNKAI OBSOLETE]]    if(stropt.spacegroupnumberoption==3) stropt.spacegroupnumberoption=2;
	// [[JUNKAI OBSOLETE]]    stropt=WyckoffPOSITIONS(stropt.spacegroupnumber,stropt.spacegroupnumberoption,stropt);
	// [[JUNKAI OBSOLETE]]    stropt.title=label_library+" "+"#"+aurostd::utype2string(stropt.spacegroupnumber)+" - "+"("+PARAMS->label+")";  // SG# AS Wahyu wants
	// [[JUNKAI OBSOLETE]]    stropt.title+=" - "+stropt.prototype+" "+_ICSD_STRING_;                                                // SG# as Wahyu wants
	// [[JUNKAI OBSOLETE]]    stropt.title=stropt.title+" (WICKOFF "+stropt.spacegroup+" "+stropt.spacegrouplabel+")";
	// [[JUNKAI OBSOLETE]]    for(uint i=0;i<stropt.atoms.size();i++) {stropt.atoms.at(i).name=speciesX.at(stropt.atoms.at(i).type);stropt.atoms.at(i).CleanName();}
	// [[JUNKAI OBSOLETE]]    // stropt=LATTICE::Conventional_Lattice_Structure(stropt); // put it in the right order...
	// [[JUNKAI OBSOLETE]]    stropt.bravais_lattice_type=LATTICE::SpaceGroup2Lattice(stropt.spacegroupnumber);
	// [[JUNKAI OBSOLETE]]    stropt.bravais_lattice_variation_type=LATTICE::SpaceGroup2LatticeVariation(stropt.spacegroupnumber,stropt);
	// [[JUNKAI OBSOLETE]]    stropt.bravais_lattice_system=stropt.bravais_lattice_type; // FIX
	// [[JUNKAI OBSOLETE]]    stropt.FixLattices();
	// [[JUNKAI OBSOLETE]]    stropt_sp=stropt;stropt_sc=stropt;stropt_sp.bravais_lattice_type="X";  // so it creates species and stuff.
	// [[JUNKAI OBSOLETE]]    step++;if(LDEBUG) *voss << "DEBUG: (aflowlib::PrototypeLibraries) (aflow_xproto.cpp): step=" << step << "  eps=" << eps << "  epsang=" << epsang << " ";
	// [[JUNKAI OBSOLETE]]    // LATTICE::Standard_Lattice_Structure(stropt,stropt_sp,stropt_sc,eps,epsang); // STEFANO OLD VERSION
	// [[JUNKAI OBSOLETE]]    ss=0; // JUNKAI
	// [[JUNKAI OBSOLETE]]    LATTICE::Standard_Lattice_Structure(stropt,stropt_sp,stropt_sc,eps,epsang,ss,_EPS_); //JUNKAI
	// [[JUNKAI OBSOLETE]]    if(LDEBUG) *voss << str_sp.bravais_lattice_variation_type << endl;
	// [[JUNKAI OBSOLETE]]    if(stropt_sp.bravais_lattice_type==stropt.bravais_lattice_type) { // gotta a wrong OPTION
	// [[JUNKAI OBSOLETE]]      *voss << "WARNING: (aflow_xproto.cpp): PARAMS->label=" << PARAMS->label << "  wrong option in lattice, put option=" << stropt.spacegroupnumberoption << endl;
	// [[JUNKAI OBSOLETE]]      str=stropt;str_sp=stropt_sp;str_sc=stropt_sc;
	// [[JUNKAI OBSOLETE]]    }
	// [[JUNKAI OBSOLETE]]  }
	// [[JUNKAI OBSOLETE]]  eps=eps/3.0;epsang=epsang/3.0;
	// [[JUNKAI OBSOLETE]]}

	if(str_sp.bravais_lattice_type!=str.bravais_lattice_type) {
	  oss << "ERROR (aflow_xproto.cpp): PARAMS->label=" << PARAMS->label << "  " << str_sp.bravais_lattice_type << " " << str.bravais_lattice_type << " "
	      << "FOUND Lattice differs from the ICSD one." << endl;
	  str_sp.error_flag=TRUE;         // flag TRUE is error
	  str_sp.error_string="ERROR: Found lattice="+str_sp.bravais_lattice_type+"  SGlattice="+str.bravais_lattice_type+"  ";
	}
	// *voss << str_sp.bravais_lattice_type << endl;
	// exit(0);
	str_sp.FixLattices();str_sp.BringInCell();
	str_sc.FixLattices();str_sc.BringInCell();
	str=str_sp; // standard primitive
      } else {
	//     cerr << "PRIM not CONV" << endl;
	uint natoms=str.atoms.size();
	str.GetPrimitive(0.005);
	if(natoms<str.atoms.size()) fixed=TRUE;
      }
      if(!fixed) str.GetPrimitive(0.05); // now you are totally screwed..... but you save plenty of ab-initio time. Wisdom of stefano.
    }
    // done
    str.BringInCell();
    str.FixLattices();
    str.MakeBasis();
    str.species_pp=str.species;
    // done
    if(LDEBUG) {
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.num_each_type.size()="<<str.num_each_type.size()<<": "; for(uint i=0;i<str.num_each_type.size();i++) *voss << str.num_each_type.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.comp_each_type.size()="<<str.comp_each_type.size()<<": "; for(uint i=0;i<str.comp_each_type.size();i++) *voss << str.comp_each_type.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species.size()="<<str.species.size()<< ": "; for(uint i=0;i<str.species.size();i++) *voss << str.species.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp.size()="<<str.species_pp.size()<< ": "; for(uint i=0;i<str.species_pp.size();i++) *voss << str.species_pp.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_type.size()="<<str.species_pp_type.size()<< ": "; for(uint i=0;i<str.species_pp_type.size();i++) *voss << str.species_pp_type.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_version.size()="<<str.species_pp_version.size()<< ": "; for(uint i=0;i<str.species_pp_version.size();i++) *voss << str.species_pp_version.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_ZVAL.size()="<<str.species_pp_ZVAL.size()<< ": "; for(uint i=0;i<str.species_pp_ZVAL.size();i++) *voss << str.species_pp_ZVAL.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_pp_vLDAU.size()="<<str.species_pp_vLDAU.size()<< ": "; for(uint i=0;i<str.species_pp_vLDAU.size();i++) *voss << str.species_pp_vLDAU.at(i).size() << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_volume.size()="<<str.species_volume.size()<< ": "; for(uint i=0;i<str.species_volume.size();i++) *voss << str.species_volume.at(i) << " "; *voss << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [8] str.species_mass.size()="<<str.species_mass.size()<< ": "; for(uint i=0;i<str.species_mass.size();i++) *voss << str.species_mass.at(i) << " "; *voss << endl;
      *voss << str << endl;
      *voss << "DEBUG: (aflowlib::PrototypeLibraries) [X] READY TO RETURN" << endl;
    }

    if(isTET) aflowlib::PrototypeFixTET(*voss,str,optionsTET);
    return str;
  }
} // namespace aflowlib
