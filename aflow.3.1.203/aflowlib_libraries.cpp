// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AFLOWLIB_LIBRARIES_CPP_
#define _AFLOWLIB_LIBRARIES_CPP_

#include "aflow.h"
#include "aflowlib.h"
#include "aflow_pflow.h"
#include "aflow_bader.h"

using std::vector;
using std::deque;
using std::string;

vector<string> vLibrary_ALL;vector<vector<string> > vLibrary_ALL_tokens;
vector<string> vLibrary_ICSD;vector<vector<string> > vLibrary_ICSD_tokens;
vector<string> vLibrary_LIB1;vector<vector<string> > vLibrary_LIB1_tokens;
vector<string> vLibrary_LIB2;vector<vector<string> > vLibrary_LIB2_tokens;
vector<string> vLibrary_LIB3;vector<vector<string> > vLibrary_LIB3_tokens;
vector<string> vLibrary_LIB4;vector<vector<string> > vLibrary_LIB4_tokens;
vector<string> vLibrary_LIB5;vector<vector<string> > vLibrary_LIB5_tokens;
vector<string> vLibrary_LIB6;vector<vector<string> > vLibrary_LIB6_tokens;
vector<string> vLibrary_LIB7;vector<vector<string> > vLibrary_LIB7_tokens;
vector<string> vLibrary_LIB8;vector<vector<string> > vLibrary_LIB8_tokens;
vector<string> vLibrary_LIB9;vector<vector<string> > vLibrary_LIB9_tokens;

#define SPECIE_RAW_LIB3 string("Ag,Al,As,Au,B_h,Bi_d,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,La,Mg_pv,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Ta_pv,Te,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")

bool AFLOWLIB_VERBOSE=TRUE; // FALSE;
#define _EPSILON_COMPOSITION_ 0.001

#define USE_AFLOW_SG
//#define USE_PLATON_SG

#define AFLOW_MAX_ARGV 1024

// ******************************************************************************************************************************************************
using aurostd::FileExist;

/*
  double EnthalpyReference(string pseudopotential,string type,deque<double> LDAU, ostream& oss);
  bool EnthalpyReference(string pp,string type,deque<double> LDAU,string& gs,double& enthalpy_atom,
  double& volume_atom,double& spin_atom, ostream& oss);
  bool EnthalpyReferenceAvailable(string pp,string type,deque<double> LDAU);
*/

double EnthalpyReference(string pp_version,deque<double> LDAU, ostream& oss);
bool EnthalpyReference(string pp_version,deque<double> LDAU,string& gs,double& enthalpy_atom,double& volume_atom,double& spin_atom, ostream& oss);
bool EnthalpyReferenceAvailable(string pp_version,deque<double> LDAU);

// ***************************************************************************
// aflowlib::TokenPresentAFLOWLIB
// ***************************************************************************
namespace aflowlib {
  bool TokenPresentAFLOWLIB(string line,string query) { return aurostd::substring2bool(line,query,TRUE); }
}
// ***************************************************************************
// aflowlib::TokenExtractAFLOWLIB
// ***************************************************************************
namespace aflowlib {
  string TokenExtractAFLOWLIB(string line,string query) {
    if(!TokenPresentAFLOWLIB(line,query)) return "";
    vector<string> tokens,tk;aurostd::string2tokens(line,tokens,"|");
    string qquery=query;aurostd::StringSubst(qquery,"=","");aurostd::StringSubst(qquery," ","");
    for(uint i=0;i<tokens.size();i++) {
      aurostd::StringSubst(tokens.at(i)," ","");
      aurostd::string2tokens(tokens.at(i),tk,"=");
      if(tk.size()==2) { //   cerr << tk.at(0) << endl;
        if(tk.at(0)==qquery) { return tk.at(1); }
        // return aurostd::substring2string(tokens.at(i),query,TRUE);
      }
    }
    return "";
  }
}

namespace aflowlib {
  bool TokenExtractAFLOWLIB(string line,string query,string &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,string &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=TokenExtractAFLOWLIB(line,query);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,int &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,int &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<int>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,uint &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,uint &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<uint>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,float &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,float &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<float>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,double &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,double &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<double>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<string> &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<string> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),value,",");return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<int> &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<int> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectorint(vvalue);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<uint> &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<uint> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectoruint(vvalue);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<float> &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<float> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectorfloat(vvalue);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<double> &value) {
    // cerr << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<double> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectordouble(vvalue);return TokenPresentAFLOWLIB(line,query); }
}

// ******************************************************************************************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(vector<string> vspecies,                         // IN  [0,nspecies[ the species Ag,Cd,...   nspecies=number of these items    nspecies=naries
			vector<string>& vspecies_pp,                     // IN  [0,nspecies[ the pseudopotentials Ag_pv, Cd_sv
			vector<vector<string> >& vList,                  // OUT [0,naries[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
			vector<vector<vector<string> > > &vList_tokens,  // OUT [0,naries[*[0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
			vector<vector<vector<string> > > &vList_species, // OUT [0,naries[*[0,vList.size()[*nspecies  returns the species present
			vector<double> &vList_Hmin,                      // OUT [0,nspecies[ returns the min enthalpy for reference
			vector<string> &vList_Pmin,                      // OUT [0,nspecies[ returns the prototype for reference
			vector<uint> &vList_Imin,                        // OUT [0,nspecies[ returns the line index of vList.at(ispecies), in which we have the min enthalpy for reference
			vector<vector<vector<double> > > &vList_concs,   // OUT [0,naries[*[0,vList.size()[*[0.nspecies[ the concentrations AxAyCz... where x+y+z=1 and it contains also ZEROS so that 0 0.25 0.75 is allowed
			vector<vector<double> > &vList_Hf) {             // OUT [0,naries[*[0,vList.size()[ returns the formation enthalpy of the list
    bool LVERBOSE=FALSE;
    bool LIBVERBOSE=TRUE;
    if(LVERBOSE) cerr << "GREP_Species_ALL: START" << endl;
    if(vspecies.size()==1) LOAD_Library_ALL("no_aflowlib_lib3,no_aflowlib_lib4,no_aflowlib_lib5,no_aflowlib_lib6,no_aflowlib_lib7,no_aflowlib_lib8,no_aflowlib_lib9,no_aflowlib_lib2",LIBVERBOSE);
    else LOAD_Library_ALL(LIBVERBOSE); // might not have been loaded, yet
    // clean stuff

    // start
    vector<string> tokens,tokens_species,tokens_species_pp;
    vector<uint> tokens_composition;
    uint natoms,nspecies=vspecies.size(),naries=vspecies.size();
    string species;
    if(nspecies==0) return 0;
    if(vspecies_pp.size()==0 || vspecies.size()!=vspecies_pp.size())  { // generate some default ones
      vspecies_pp.clear(); for(uint i=0;i<nspecies;i++)vspecies_pp.push_back(AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)));
    }
    for(uint i=0;i<nspecies;i++)
      if(vspecies_pp.at(i).empty()) vspecies_pp.at(i)=AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)); // in the event that I pass an empty string

    // generate space
    vector<string> dummy_vstring;
    vList.clear(); for(uint i=0;i<naries;i++) vList.push_back(dummy_vstring);                  // create space  [0,naries[*[0,vList.size()[  <strings>
    vector<vector<string> > dummy_vvstring;
    vList_tokens.clear(); for(uint i=0;i<naries;i++) vList_tokens.push_back(dummy_vvstring);   // create space  [0,naries[*[0,vList.size()[*[0,size_tokens[  <strings>
    vList_Hmin.clear(); for(uint i=0;i<nspecies;i++) vList_Hmin.push_back(1E6);                // create space  [0,nspecies[ <double>
    vList_Pmin.clear(); for(uint i=0;i<nspecies;i++) vList_Pmin.push_back("");                 // create space  [0,nspecies[ <string>
    vList_Imin.clear(); for(uint i=0;i<nspecies;i++) vList_Imin.push_back(0);                  // create space  [0,nspecies[ <uint>
    vector<double> dummy_vdouble;
    vList_Hf.clear(); for(uint i=0;i<naries;i++) vList_Hf.push_back(dummy_vdouble);            // create space  [0,naries[*[0,vList.size()[  <double>
    vector<vector<double> > dummy_vvdouble;
    vList_concs.clear(); for(uint i=0;i<naries;i++) vList_concs.push_back(dummy_vvdouble);     // create space  [0,naries[*[0,vList.size()[*[0.nspecies[   <double>
    vList_species.clear(); for(uint i=0;i<naries;i++) vList_species.push_back(dummy_vvstring); // create space  [0,naries[*[0,vList.size()[*[0.nspecies[   <string>

    uint num_species=0;string species_pp;double enthalpy_atom,enthalpy;
    // for(uint ispecies=1;ispecies<=nspecies;ispecies++) {
    for(uint iall=0;iall<vLibrary_ALL.size();iall++) {
      bool found_some_species=FALSE;
      for(uint ispecies=0;((ispecies<nspecies) && (!found_some_species));ispecies++)
        if(aurostd::substring2bool(vLibrary_ALL.at(iall),vspecies_pp.at(ispecies),TRUE)) found_some_species=TRUE;
      if(found_some_species) {
        TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"nspecies=",num_species); // search for num_species
        for(uint ispecies=1;ispecies<=nspecies;ispecies++) {
          bool found_nspecies=FALSE,found_species_pp=FALSE;
          if(num_species==ispecies) found_nspecies=TRUE;
          if(found_nspecies && num_species==1) {
            if(!TokenPresentAFLOWLIB(vLibrary_ALL.at(iall),"LIB1")) found_nspecies=FALSE; // patch to use only pure
            if(found_nspecies) {
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Rh_pv" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Ni_pv" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
              // if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="La" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A1") found_nspecies=FALSE;  // patches for incorrect relaxation
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A7") found_nspecies=FALSE;  // patches for incorrect relaxation
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A8") found_nspecies=FALSE;  // patches for incorrect relaxation
              // fix La Mn and Hg
            }
          }
          // search for species_pp
          if(found_nspecies) {
            TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=",species_pp);
            TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=",tokens_species);
            uint species_pp_matches=0;
            for(uint k1=0;k1<vspecies_pp.size();k1++)
              for(uint k2=0;k2<tokens_species.size();k2++)
                if(vspecies_pp.at(k1)==tokens_species.at(k2)) species_pp_matches++;
            if(species_pp_matches==ispecies) found_species_pp=TRUE;
          }
          if(found_species_pp) {
            // cerr << ispecies << " " << vLibrary_ALL.at(iall) << endl;
            string prototype="";
            TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=",prototype);
            if(prototype=="64" || prototype=="65" || prototype=="549" || prototype=="550") {  // check for broken calcs
              if(LVERBOSE) cerr << "GREP_Species_ALL: Removing prototype=" << prototype << endl;
            } else {
              vList.at(ispecies-1).push_back(vLibrary_ALL.at(iall));       // add lines
              aurostd::string2tokens(vLibrary_ALL.at(iall),tokens,"|");    // add tokens
              vList_tokens.at(ispecies-1).push_back(tokens);               // add tokens
              if(LVERBOSE) if(!mod((int) vList.at(ispecies-5).at(vList.at(ispecies-1).size()-1).size(),1)) cerr << vLibrary_ALL_tokens.at(iall).at(0) << " " << species_pp << endl;  // for DEBUG printout
            }
          }
        }
      }
    }
    if(LVERBOSE) for(uint ispecies=0;ispecies<nspecies;ispecies++) cerr << ispecies << " " << vList.at(ispecies).size() << " " << vList_tokens.at(ispecies).size() << endl;
    // LIST PREPARED
    if(LVERBOSE) cerr << "GREP_Species_ALL: LIST PREPARED" << endl;
    // creating reference enthalpies
    if(LVERBOSE) cerr << "GREP_Species_ALL: creating reference enthalpies" << endl;
    for(uint ispecies=0;ispecies<nspecies;ispecies++) {
      if(LVERBOSE) cerr << "GREP_Species_ALL: getting [" << vspecies.at(ispecies) << "]" << endl;
      for(uint j=0;j<vList.at(0).size();j++) { // only pure
        if(LVERBOSE) cerr << j << " " << vList_tokens.at(0).at(j).at(0) << endl;
        TokenExtractAFLOWLIB(vList.at(0).at(j),"species=",species);
        TokenExtractAFLOWLIB(vList.at(0).at(j),"enthalpy_atom=",enthalpy_atom);
        if(species==vspecies.at(ispecies)) {
          if(LVERBOSE) cerr << j << " " << vList_tokens.at(0).at(j).at(0) << endl;
          if(enthalpy_atom<vList_Hmin.at(ispecies)) {
            if(LVERBOSE) cerr << "GREP_Species_ALL: ***** FOUND " << j << " " << vList_tokens.at(0).at(j).at(0) << endl;
            vList_Hmin.at(ispecies)=enthalpy_atom;
            vList_Pmin.at(ispecies)=TokenExtractAFLOWLIB(vList.at(0).at(j),"prototype=");
            vList_Imin.at(ispecies)=j;
          }
        }
      }
    }

    // JUNKAI THIS IS THE WAY TO READ vList_Hmin.at(ispecies) vList_Pmin.at(ispecies) vList_Imin.at(ispecies) right use of vList_tokens.at(0,1,2 A,AB,ABC).at(vList_Imin.at(ispecies)).at(TOKENS INDEX THAT YOU WANT)
    // HOWEVER YOU CAN EXTRACT INFORMATION WITH TokenExtractAFLOWLIB(vList..at(j),"species=",species);
    for(uint ispecies=0;ispecies<nspecies;ispecies++) {
      if(LVERBOSE) cerr << vspecies.at(ispecies) << " " << vList_Hmin.at(ispecies) << " " << vList_Pmin.at(ispecies) << " " << vList_Imin.at(ispecies) << " " << vList_tokens.at(0).at(vList_Imin.at(ispecies)).at(0) << endl;
    }
    for(uint iary=0;iary<naries;iary++) {
      if(LVERBOSE) cerr << vList_tokens.at(iary).size() << endl;
    }

    // creating formation enthalpies
    for(uint iaries=0;iaries<naries;iaries++) {
      for(uint j=0;j<vList.at(iaries).size();j++) { // only pure
        vector<string> vspecies_local;
        vector<double> vconcs_local;
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"natoms=",natoms); // extract natoms
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"enthalpy_cell=",enthalpy); // extract enthalpy
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"species=",tokens_species);
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"species_pp=",tokens_species_pp);
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"composition=",tokens_composition);

        if(tokens_species.size()!=tokens_composition.size()) { cerr << "GREP_Species_ALL: tokens_species.size()!=tokens_composition.size()" << endl;exit(0); }
        // cerr << "GREP_Species_ALL: E=" << enthalpy << " ";
        for(uint k=0;k<nspecies;k++) // order by species, no entries
          for(uint j=0;j<tokens_species.size();j++)
            // if(tokens_species.at(j)==vspecies.at(k)) { // found specie k associated with position k
            if(tokens_species_pp.at(j)==vspecies_pp.at(k)) { // found specie_pp k associated with position k
              enthalpy=enthalpy-(double) vList_Hmin.at(k)*tokens_composition.at(j); // remove formation
              vspecies_local.push_back(vspecies_pp.at(k));
              vconcs_local.push_back((double) tokens_composition.at(j)/(double) natoms);
            }
        enthalpy=enthalpy/(double) natoms; // normalize per atom
        if(aurostd::abs(enthalpy)<0.0001) enthalpy=0.0; // <0.1meV equal zero
        // cerr << enthalpy << endl;

        vList_Hf.at(iaries).push_back(enthalpy);              // plug scalar enthalpy_formation_atom
        vList_species.at(iaries).push_back(vspecies_local); // plug VECTOR species
        vList_concs.at(iaries).push_back(vconcs_local);     // plug VECTOR concentrations
      }
    }

    for(uint iaries=0;iaries<naries;iaries++) {
      for(uint jentry=0;jentry<vList.at(iaries).size();jentry++) { // only pure
	if(LVERBOSE) cerr << vList_tokens.at(iaries).at(jentry).at(0) << "  " << vList_tokens.at(iaries).at(jentry).at(1) << "   Hf = " << vList_Hf.at(iaries).at(jentry) << "  Compound = ";
	for(uint k=0;k<vList_species.at(iaries).at(jentry).size();k++) if(LVERBOSE) cerr << vList_species.at(iaries).at(jentry).at(k) << "_" << vList_concs.at(iaries).at(jentry).at(k) << " ";
	if(LVERBOSE) cerr << endl;
      }
    }

    // RETURN
    uint numentries=0;
    for(uint iaries=0;iaries<naries;iaries++) numentries+=vList.at(iaries).size();
    return numentries;
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(string species,                                   // IN   the species Ag,...   nspecies=number of these items
			string& species_pp,                               // IN   the pseudopotentials Ag_pv,
			vector<string> & vList,                           // OUT  [0,nspecies[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
			vector<vector<string> > &vList_tokens,            // OUT [0,nspecies[*[0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
			double& List_Hmin,                                // OUT  returns the min enthalpy for reference
			string& List_Pmin,                                // OUT  returns the prototype for reference
			uint& List_Imin) {                                // OUT  returns the line index of vList, in which we have the min enthalpy for reference
    vector<string> vspecies;vspecies.push_back(species);
    vector<string> vspecies_pp;vspecies_pp.push_back(species_pp);

    vector<vector<string> > _vList;
    vector<vector<vector<string> > > _vList_tokens;
    vector<vector<vector<string> > > vList_species;
    vector<double> vList_Hmin;
    vector<string> vList_Pmin;
    vector<uint> vList_Imin;
    vector<vector<vector<double> > > vList_concs;
    vector<vector<double> > vList_Hf;

    uint out=GREP_Species_ALL(vspecies,vspecies_pp,_vList,_vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);
    species_pp=vspecies_pp.at(0);
    List_Hmin=vList_Hmin.at(0);  // they are for ONE SPECIE =
    List_Pmin=vList_Pmin.at(0);  // they are for ONE SPECIE
    List_Imin=vList_Imin.at(0);  // they are for ONE SPECIE
    for(uint i=0;i<_vList.at(0).size();i++) vList.push_back(_vList.at(0).at(i));
    for(uint i=0;i<_vList_tokens.at(0).size();i++) vList_tokens.push_back(_vList_tokens.at(0).at(i));

    return out;
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList,vector<vector<vector<string> > > &vList_tokens) {
    vector<vector<vector<string> > > vList_species;
    vector<vector<vector<double> > > vList_concs;
    vector<vector<double> > vList_Hf;
    vector<double> vList_Hmin;
    vector<string> vList_Pmin;
    vector<uint> vList_Imin;

    return GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList) {
    bool LVERBOSE=FALSE;
    bool LIBVERBOSE=TRUE;
    if(LVERBOSE) cerr << "GREP_Species_ALL: START" << endl;
    if(vspecies.size()==1) LOAD_Library_ALL("no_aflowlib_lib3,no_aflowlib_lib4,no_aflowlib_lib5,no_aflowlib_lib6,no_aflowlib_lib7,no_aflowlib_lib8,no_aflowlib_lib9,no_aflowlib_lib2",LIBVERBOSE);
    else LOAD_Library_ALL(LIBVERBOSE); // might not have been loaded, yet
    // clean stuff

    // start //
    vector<string> tokens,tokens_species,tokens_species_pp;
    uint nspecies=vspecies.size(),naries=vspecies.size();
    string species;
    if(nspecies==0) return 0;
    if(vspecies_pp.size()==0 || vspecies.size()!=vspecies_pp.size())  { // generate some default ones
      vspecies_pp.clear();
      for(uint i=0;i<nspecies;i++) {
	vspecies_pp.push_back(AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)));
      }
    }
    for(uint i=0;i<nspecies;i++) {
      if(vspecies_pp.at(i).empty()) {
	vspecies_pp.at(i)=AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)); // in the event that I pass an empty string
      }
    }
    // generate space
    vector<string> dummy_vstring;
    vList.clear(); for(uint i=0;i<naries;i++) vList.push_back(dummy_vstring);                  // create space  [0,naries[*[0,vList.size()[  <strings>

    cerr << "START GREP" << endl;

    uint num_species=0;string species_pp;
    // for(uint ispecies=1;ispecies<=nspecies;ispecies++) {
    string vLibrary_ALL_at_iall;
    for(uint iall=0;iall<vLibrary_ALL.size();iall++) {
      vLibrary_ALL_at_iall=vLibrary_ALL.at(iall);
      bool found_some_species=FALSE;
      for(uint ispecies=0;((ispecies<nspecies) && (!found_some_species));ispecies++)
	if(aurostd::substring2bool(vLibrary_ALL_at_iall,vspecies_pp.at(ispecies),TRUE)) found_some_species=TRUE;
      if(found_some_species) {
	for(uint ispecies=1;ispecies<=nspecies;ispecies++) {
	  if(aurostd::substring2bool(vLibrary_ALL_at_iall,string("nspecies="+aurostd::utype2string(ispecies)))) {
	    TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"nspecies=",num_species); // search for num_species
	    // cerr << " [nspecies=" << aurostd::utype2string(ispecies) << "] " << endl;
	    bool found_nspecies=FALSE,found_species_pp=FALSE;
	    if(num_species==ispecies) found_nspecies=TRUE;
	    if(found_nspecies && num_species==1) {
	      if(!TokenPresentAFLOWLIB(vLibrary_ALL_at_iall,"LIB1")) found_nspecies=FALSE; // patch to use only pure
	      if(found_nspecies) {
		if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Rh_pv" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
		if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Ni_pv" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
		// if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="La" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A1") found_nspecies=FALSE;  // patches for incorrect relaxation
		if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A7") found_nspecies=FALSE;  // patches for incorrect relaxation
		if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A8") found_nspecies=FALSE;  // patches for incorrect relaxation
		// fix La Mn and Hg
	      }
	    }
	    // search for species_pp
	    if(found_nspecies) {
	      TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=",species_pp);
	      TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=",tokens_species);
	      uint species_pp_matches=0;
	      for(uint k1=0;k1<vspecies_pp.size();k1++)
		for(uint k2=0;k2<tokens_species.size();k2++)
		  if(vspecies_pp.at(k1)==tokens_species.at(k2)) species_pp_matches++;
	      if(species_pp_matches==ispecies) found_species_pp=TRUE;
	    }
	    if(found_species_pp) {
	      // cerr << ispecies << " " << vLibrary_ALL_at_iall << endl;
	      string prototype="";
	      TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=",prototype);
	      if(prototype=="64" || prototype=="65" || prototype=="549" || prototype=="550") {  // check for broken calcs
		if(LVERBOSE) cerr << "GREP_Species_ALL: Removing prototype=" << prototype << endl;
	      } else {
		vList.at(ispecies-1).push_back(vLibrary_ALL_at_iall);       // add lines
		aurostd::string2tokens(vLibrary_ALL_at_iall,tokens,"|");    // add tokens
	      }
	    }
	  }
	}
      }
    }

    cerr << "END GREP" << endl;
    // RETURN
    uint numentries=0;
    for(uint iaries=0;iaries<naries;iaries++) numentries+=vList.at(iaries).size();
    return numentries;
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(string species,                                   // IN   the species Ag,...   nspecies=number of these items
			string& species_pp,                               // IN   the pseudopotentials Ag_pv,
			double& List_Hmin,                                // OUT  returns the min enthalpy for reference
			string& List_Pmin) {                              // OUT  returns the prototype for reference
    vector<string> vList;
    vector<vector<string> > vList_tokens;
    uint List_Imin;
    return GREP_Species_ALL(species,species_pp,vList,vList_tokens,List_Hmin,List_Pmin,List_Imin);
  }
}

// ***************************************************************************
// aflowlib::LIBS_EFormation
// ***************************************************************************
namespace aflowlib {
  bool LIBS_EFormation(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::ALIBRARIES: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    bool LVERBOSE=FALSE;
    aflowlib::LOAD_Library_LIBRARY("pure","",LVERBOSE);

    string species,species_pp;
    // chose somehow the species
    if(tokens.size()>=1) { species=tokens.at(0); } else { species="Nb"; }

    //  cerr << "vLibrary_ALL.size()=" << vLibrary_ALL.size() << " vLibrary_ALL_tokens.size()=" << vLibrary_ALL_tokens.size() << endl;

    vector<string> vspecies,vspecies_pp;
    aurostd::string2tokens(SPECIE_RAW_LIB3,vspecies,",");
    for(uint i=0;i<vspecies.size();i++) {
      vspecies.at(i)=KBIN::VASP_PseudoPotential_CleanName(vspecies.at(i));
      vspecies_pp.push_back("");//AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)));

      vector<string> vList;
      vector<vector<string> > vList_tokens;
      double List_Hmin=0.0;
      string List_Pmin="";
      uint List_Imin=0;

      species_pp=AVASP_Get_PseudoPotential_PAW_PBE(species);

      GREP_Species_ALL(vspecies.at(i),          // IN    the species Ag,...   nspecies=number of these items
		       vspecies_pp.at(i),       // IN    the pseudopotentials Ag_pv,
		       vList,            // OUT   [0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
		       vList_tokens,     // OUT   [0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
		       List_Hmin,        // OUT   returns the min enthalpy for reference
		       List_Pmin,        // OUT   returns the prototype for reference
		       List_Imin);
      cerr << vspecies.at(i) << " " << vspecies_pp.at(i) << " " << List_Hmin << " " << List_Pmin << " " << List_Imin << " " << vList_tokens.at(List_Imin).at(0) << endl;
    }
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::ALIBRARIES
// ***************************************************************************
namespace aflowlib {
  bool ALIBRARIES(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::ALIBRARIES: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(tokens.size()>10) {
      init::ErrorOption(cout,options,"aflowlib::ALIBRARIES","aflow --qhull[=Nb[,Pt[,Rh[,...]]]]");
      exit(0);
    }

    bool LVERBOSE=TRUE;

    LOAD_Library_ALL(LVERBOSE);
    //  LOAD_Library_ALL();
    cerr << "vLibrary_ALL.size()=" << vLibrary_ALL.size() << " vLibrary_ALL_tokens.size()=" << vLibrary_ALL_tokens.size() << endl;

    vector<string> vspecies,vspecies_pp;
    vector<vector<string> > vList;
    vector<vector<vector<string> > > vList_tokens;
    vector<double> vList_Hmin;
    vector<string> vList_Pmin;
    vector<uint> vList_Imin;
    vector<vector<vector<string> > > vList_species;
    vector<vector<vector<double> > > vList_concs;
    vector<vector<double> > vList_Hf;

    // uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);  // complete
    // uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imi);

    // [OBSOLETE] double List_Hmin=1E6;
    string List_Pmin="";
    // [OBSOLETE] uint List_Imin=0;

    // chose somehow the species
    if(tokens.size()==0) {
      vspecies.push_back("Nb");vspecies.push_back("Pt");vspecies.push_back("Rh");  // DO IN ALPHABETIC ORDER OTHERWISE COMPLAIN
    } else {
      for(uint i=0;i<tokens.size();i++)
	vspecies.push_back(tokens.at(i));
    }

    //  uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens);
    uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);
    // uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens);
    //  uint nentries=GREP_Species_ALL(vspecies.at(0),"",vList,vList_tokens,List_Hmin,List_Pmin,List_Imin);
    for(uint iaries=0;iaries<vList.size();iaries++) {
      //  cerr << iaries << " " << vList.at(iaries).size() << endl;
      //    for(uint j=0;j<vList.at(iaries).size();j++) cerr << vList.at(iaries).at(j) << endl;
    }

    for(uint iaries=0;iaries<vList.size();iaries++) {
      if(LVERBOSE) cerr << "***************[" << iaries+1 << "]*************" << endl;
      for(uint jentry=0;jentry<vList.at(iaries).size();jentry++) { // only pure
	if(LVERBOSE) cerr << vList_tokens.at(iaries).at(jentry).at(0) << "  " << vList_tokens.at(iaries).at(jentry).at(1) << "   Hf = " << vList_Hf.at(iaries).at(jentry) << "  Compound = ";
	for(uint k=0;k<vList_species.at(iaries).at(jentry).size();k++) if(LVERBOSE) cerr << vList_species.at(iaries).at(jentry).at(k) << "_" << vList_concs.at(iaries).at(jentry).at(k) << " ";
	if(LVERBOSE) cerr << endl;
      }
    }

    if(LVERBOSE) cerr << nentries << endl;
    cerr << "vList.at(iaries).size()="; for(uint iaries=0;iaries<vList.size();iaries++) cerr << vList.at(iaries).size() << " "; cerr << endl;
    cerr << "vList_tokens.at(iaries).size()="; for(uint iaries=0;iaries<vList_tokens.size();iaries++) cerr << vList_tokens.at(iaries).size() << " "; cerr << endl;
    // cerr << "List_Hmin=" << List_Hmin << endl;
    // cerr << "List_Pmin=" << List_Pmin << endl;
    // cerr << "List_Imin=" << List_Imin << endl;

    // exit(0);
    if(LDEBUG) cerr << "aflowlib::ALIBRARIES: END" << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LOAD_Library_LIBRARY
// ***************************************************************************
namespace aflowlib {
  uint LOAD_Library_LIBRARY(string file,string grep,bool LVERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // LDEBUG=TRUE;
    string FileLibrary="";
    vector<string> tokens;
    bool found=FALSE;
    vector<string>* pvLibrary=NULL; vector<vector<string> >* pvLibrary_tokens=NULL;
    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    if(!found && (file=="aflowlib_lib1" ||
		  file=="aflowlib_lib2" ||
		  file=="aflowlib_lib3" ||
		  file=="aflowlib_lib4" ||
		  file=="aflowlib_lib5" ||
		  file=="aflowlib_lib6" ||
		  file=="aflowlib_lib7" ||
		  file=="aflowlib_lib8" ||
		  file=="aflowlib_lib9" ||
		  file=="aflowlib_icsd" )) { // from aflow_data
      if(file=="aflowlib_icsd") { // ICSD IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_ICSD;pvLibrary_tokens=&vLibrary_ICSD_tokens;
	init::InitGlobalObject("aflowlib_icsd",grep,LVERBOSE);
	if(XHOST_aflowlib_icsd.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_icsd,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib1") { // LIB1 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB1;pvLibrary_tokens=&vLibrary_LIB1_tokens;
	init::InitGlobalObject("aflowlib_lib1",grep,LVERBOSE);
	if(XHOST_aflowlib_lib1.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib1,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib2") { // LIB2 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB2;pvLibrary_tokens=&vLibrary_LIB2_tokens;
	init::InitGlobalObject("aflowlib_lib2",grep,LVERBOSE);
	if(XHOST_aflowlib_lib2.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib2,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib3") { // LIB3 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB3;pvLibrary_tokens=&vLibrary_LIB3_tokens;
	init::InitGlobalObject("aflowlib_lib3",grep,LVERBOSE);
	if(XHOST_aflowlib_lib3.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib3,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib4") { // LIB4 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB4;pvLibrary_tokens=&vLibrary_LIB4_tokens;
	init::InitGlobalObject("aflowlib_lib4",grep,LVERBOSE);
	if(XHOST_aflowlib_lib4.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib4,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib5") { // LIB5 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB5;pvLibrary_tokens=&vLibrary_LIB5_tokens;
	init::InitGlobalObject("aflowlib_lib5",grep,LVERBOSE);
	if(XHOST_aflowlib_lib5.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib5,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib6") { // LIB6 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB6;pvLibrary_tokens=&vLibrary_LIB6_tokens;
	init::InitGlobalObject("aflowlib_lib6",grep,LVERBOSE);
	if(XHOST_aflowlib_lib6.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib6,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib7") { // LIB7 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB7;pvLibrary_tokens=&vLibrary_LIB7_tokens;
	init::InitGlobalObject("aflowlib_lib7",grep,LVERBOSE);
	if(XHOST_aflowlib_lib7.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib7,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib8") { // LIB8 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB8;pvLibrary_tokens=&vLibrary_LIB8_tokens;
	init::InitGlobalObject("aflowlib_lib8",grep,LVERBOSE);
	if(XHOST_aflowlib_lib8.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib8,(*pvLibrary));
	found=TRUE; }
      if(file=="aflowlib_lib9") { // LIB9 IS INSIDE init::InitGlobalObject
	pvLibrary=&vLibrary_LIB9;pvLibrary_tokens=&vLibrary_LIB9_tokens;
	init::InitGlobalObject("aflowlib_lib9",grep,LVERBOSE);
	if(XHOST_aflowlib_lib9.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib9,(*pvLibrary));
	found=TRUE; }
      // cerr << "(*pvLibrary_tokens).size()=" << (uint) (*pvLibrary_tokens).size() << endl;
      // cerr << "(*pvLibrary).size()=" << (*pvLibrary).size() << endl;
      if(found==TRUE) {
	if(LDEBUG || LVERBOSE) cerr << " ... size=" << (*pvLibrary).size() << " ..  ";// << endl;
      } else {
	cerr << "WWWWW  AFLOW_LIBRARY not found! " << endl;exit(0);
	return FALSE;
      }
      if((*pvLibrary_tokens).size()==0) {
	if(LDEBUG || LVERBOSE) cerr << "processing...   ";
	for(uint i=0;i<(*pvLibrary).size();i++) {
	  tokens.clear();
	  aurostd::string2tokens((*pvLibrary).at(i),tokens,"|");
	  (*pvLibrary_tokens).push_back(tokens);
	}
	if(LDEBUG || LVERBOSE) cerr << "done." << endl;
      }
      // ----------------------------------------------------------------------
    }
    return (*pvLibrary).size();
  }
}

// ***************************************************************************
// aflowlib::LOAD_Library_ALL
// ***************************************************************************
namespace aflowlib {
  uint LOAD_Library_ALL(string options,string grep,bool LVERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(!aurostd::substring2bool(options,"no_aflowlib_icsd",TRUE) && vLibrary_ICSD.size()==0) LOAD_Library_LIBRARY("aflowlib_icsd",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib1",TRUE) && vLibrary_LIB1.size()==0) LOAD_Library_LIBRARY("aflowlib_lib1",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib2",TRUE) && vLibrary_LIB2.size()==0) LOAD_Library_LIBRARY("aflowlib_lib2",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib3",TRUE) && vLibrary_LIB3.size()==0) LOAD_Library_LIBRARY("aflowlib_lib3",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib4",TRUE) && vLibrary_LIB4.size()==0) LOAD_Library_LIBRARY("aflowlib_lib4",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib5",TRUE) && vLibrary_LIB5.size()==0) LOAD_Library_LIBRARY("aflowlib_lib5",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib6",TRUE) && vLibrary_LIB6.size()==0) LOAD_Library_LIBRARY("aflowlib_lib6",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib7",TRUE) && vLibrary_LIB7.size()==0) LOAD_Library_LIBRARY("aflowlib_lib7",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib8",TRUE) && vLibrary_LIB8.size()==0) LOAD_Library_LIBRARY("aflowlib_lib8",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib9",TRUE) && vLibrary_LIB9.size()==0) LOAD_Library_LIBRARY("aflowlib_lib9",grep,LVERBOSE);
    if(LDEBUG) cerr << "vLibrary_ICSD.size()=" << vLibrary_ICSD.size() << " vLibrary_ICSD_tokens.size()=" << vLibrary_ICSD_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB1.size()=" << vLibrary_LIB1.size() << " vLibrary_LIB1_tokens.size()=" << vLibrary_LIB1_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB2.size()=" << vLibrary_LIB2.size() << " vLibrary_LIB2_tokens.size()=" << vLibrary_LIB2_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB3.size()=" << vLibrary_LIB3.size() << " vLibrary_LIB3_tokens.size()=" << vLibrary_LIB3_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB4.size()=" << vLibrary_LIB4.size() << " vLibrary_LIB4_tokens.size()=" << vLibrary_LIB4_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB5.size()=" << vLibrary_LIB5.size() << " vLibrary_LIB5_tokens.size()=" << vLibrary_LIB5_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB6.size()=" << vLibrary_LIB6.size() << " vLibrary_LIB6_tokens.size()=" << vLibrary_LIB6_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB7.size()=" << vLibrary_LIB7.size() << " vLibrary_LIB7_tokens.size()=" << vLibrary_LIB7_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB8.size()=" << vLibrary_LIB8.size() << " vLibrary_LIB8_tokens.size()=" << vLibrary_LIB8_tokens.size() << endl;
    if(LDEBUG) cerr << "vLibrary_LIB9.size()=" << vLibrary_LIB9.size() << " vLibrary_LIB9_tokens.size()=" << vLibrary_LIB9_tokens.size() << endl;
    // now make vLibrary_ALL
    if(vLibrary_ALL.size()==0) {
      vLibrary_ALL.clear();
      vector<string> tokens;
      if(LDEBUG || LVERBOSE) cerr << "00000  MESSAGE AFLOW LIBRARY  ALL  ";
      if(LDEBUG || LVERBOSE) cerr << "loading... ";
      for(uint i=0;i<vLibrary_ICSD.size();i++) vLibrary_ALL.push_back(vLibrary_ICSD.at(i));
      for(uint i=0;i<vLibrary_LIB1.size();i++) vLibrary_ALL.push_back(vLibrary_LIB1.at(i));
      for(uint i=0;i<vLibrary_LIB2.size();i++) vLibrary_ALL.push_back(vLibrary_LIB2.at(i));
      for(uint i=0;i<vLibrary_LIB3.size();i++) vLibrary_ALL.push_back(vLibrary_LIB3.at(i));
      for(uint i=0;i<vLibrary_LIB4.size();i++) vLibrary_ALL.push_back(vLibrary_LIB4.at(i));
      for(uint i=0;i<vLibrary_LIB5.size();i++) vLibrary_ALL.push_back(vLibrary_LIB5.at(i));
      for(uint i=0;i<vLibrary_LIB6.size();i++) vLibrary_ALL.push_back(vLibrary_LIB6.at(i));
      for(uint i=0;i<vLibrary_LIB7.size();i++) vLibrary_ALL.push_back(vLibrary_LIB7.at(i));
      for(uint i=0;i<vLibrary_LIB8.size();i++) vLibrary_ALL.push_back(vLibrary_LIB8.at(i));
      for(uint i=0;i<vLibrary_LIB9.size();i++) vLibrary_ALL.push_back(vLibrary_LIB9.at(i));
      if(LDEBUG || LVERBOSE) cerr << "processing...   ";
      for(uint i=0;i<vLibrary_ALL.size();i++) {
	tokens.clear();
	aurostd::string2tokens(vLibrary_ALL.at(i),tokens,"|");
	vLibrary_ALL_tokens.push_back(tokens);
      }
      if(LDEBUG || LVERBOSE) cerr << "done." << endl;
    }
    if(LDEBUG) cerr << "vLibrary_ALL.size()=" << vLibrary_ALL.size() << " vLibrary_ALL_tokens.size()=" << vLibrary_ALL_tokens.size() << endl;
    return vLibrary_ALL.size();
  }
}

namespace aflowlib {
  uint LOAD_Library_ALL(string options,bool LVERBOSE) {
    return LOAD_Library_ALL(options,"",LVERBOSE);
  }
}

namespace aflowlib {
  uint LOAD_Library_ALL(bool LVERBOSE) {
    return LOAD_Library_ALL("","",LVERBOSE);
  }
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

// ***************************************************************************
// aflowlib::XLIB2RAW_CheckProjectFromDirectory
// ***************************************************************************
namespace aflowlib {
  string LIB2RAW_CheckProjectFromDirectory(string directory) {
    CheckMaterialServer("aflowlib::LIB2RAW_CheckProjectFromDirectory"); // must be in AFLOW_MATERIALS_SERVER
    // find from PWD
    string PROJECT_LIBRARY="NOTHING",directory_pwd=aurostd::execute2string("pwd");aurostd::StringSubst(directory_pwd,"\n","");
    if(directory_pwd=="/common/GNDSTATE") directory_pwd="/common/LIB2"; // [HISTORIC]
    aurostd::StringSubst(directory_pwd,"common/SCINT","common/ICSD"); // [HISTORIC]
    aurostd::StringSubst(directory_pwd,"common/ELPASOLITES","common/AURO"); // [HISTORIC]
    aurostd::StringSubst(directory_pwd,"LIB2/RAW","LIB2/LIB"); // [HISTORIC]

    // if not found by PWD switch to directory
    for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++)
      if(aurostd::substring2bool(directory,vAFLOW_PROJECTS_DIRECTORIES.at(i))) PROJECT_LIBRARY=vAFLOW_PROJECTS_DIRECTORIES.at(i);
    if(PROJECT_LIBRARY!="NOTHING") {
      //     cerr << "aflowlib::LIB2RAW_CheckProjectFromDirectory: FOUND from directory: " << PROJECT_LIBRARY << endl;
      return PROJECT_LIBRARY;
    }

    // cerr << directory_pwd << endl;
    for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++) {
      //   cerr << vAFLOW_PROJECTS_DIRECTORIES.at(i) << endl;
      if(aurostd::substring2bool(directory_pwd,vAFLOW_PROJECTS_DIRECTORIES.at(i))) PROJECT_LIBRARY=vAFLOW_PROJECTS_DIRECTORIES.at(i);
    }
    if(PROJECT_LIBRARY!="NOTHING") {
      //  cerr << "aflowlib::LIB2RAW_CheckProjectFromDirectory: FOUND from pwd: " << PROJECT_LIBRARY << endl;
      return PROJECT_LIBRARY;
    }

    if(PROJECT_LIBRARY=="NOTHING") {
      cerr << "aflowlib::LIB2RAW_CheckProjectFromDirectory: Nothing found from pwd or directory ... exiting... [directory=" << directory << "]   [pwd=" << directory_pwd << "]" << endl;
      exit(0);
    }
    return "/tmp";
  }
}

// ***************************************************************************
// aflowlib::XFIX_LIBRARY_ALL
// ***************************************************************************
namespace aflowlib {
  void XFIX_LIBRARY_ALL(string LIBRARY_IN,vector<string> argv) {
    CheckMaterialServer("aflowlib::XFIX_LIBRARY_ALL"); // must be in AFLOW_MATERIALS_SERVER
    stringstream aus_exec;
    string system=argv.at(2),structure=argv.at(3);   // if 3 inputs need to fix with fewer inputs
    string systemstructure=system+"/"+structure; // if 3 inputs need to fix with fewer inputs
    cout << "Fixing " << LIBRARY_IN << "/LIB/" << systemstructure << endl;
    aurostd::RemoveFile(string("_aflow_xfix_LIBRARY_"));
    aus_exec << "#!/bin/csh" << endl;
    aus_exec << "# Cleaning up " << LIBRARY_IN << " system/structure " << endl;
    aus_exec << "if !(-e " << LIBRARY_IN << "/FIX/" << system << ") then" << endl;
    aus_exec << "   mkdir " << LIBRARY_IN << "/FIX/" << system << endl;
    aus_exec << "endif" << endl;
    aus_exec << "mv " << LIBRARY_IN << "/LIB/" << systemstructure << " " << LIBRARY_IN << "/FIX/" << system << " " << endl;
    aus_exec << "mv " << LIBRARY_IN << "/RAW/" << systemstructure << " /tmp/" << system << " " << endl;
    aus_exec << "rm -rf " << LIBRARY_IN << "/RAW/" << systemstructure << endl;
    aurostd::stringstream2file(aus_exec,string("_aflow_xfix_LIBRARY_"));
    aurostd::ChmodFile("755",string("_aflow_xfix_LIBRARY_"));
    // aurostd::execute((char*) "cat ./_aflow_xfix_LIBRARY_");
    aurostd::execute((char*) "./_aflow_xfix_LIBRARY_");
    aurostd::RemoveFile(string("_aflow_xfix_LIBRARY_"));
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_ALL
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_ALL(string options,bool flag_FORCE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_ALL: BEGIN" << endl;
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_ALL: options=" << options << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_ALL: tokens.size()=" << tokens.size() << endl;
    if(tokens.size()>2) {
      init::ErrorOption(cout,options,"aflowlib::LIB2RAW_ALL",aurostd::liststring2string("aflow --lib2raw=directory","aflow --lib2raw=all[,dir]"));
      exit(0);
    }
    if(tokens.size()>=1) {
      if(tokens.at(0)!="all")  {
	init::ErrorOption(cout,options,"aflowlib::LIB2RAW_ALL",aurostd::liststring2string("aflow --lib2raw=directory","aflow --lib2raw=all[,dir]"));
	exit(0);
      }
    }

    CheckMaterialServer("aflowlib::LIB2RAW_ALL"); // must be in AFLOW_MATERIALS_SERVER
    string PROJECT_LIBRARY;
    if(tokens.size()==2) PROJECT_LIBRARY=aflowlib::LIB2RAW_CheckProjectFromDirectory(tokens.at(1));
    else PROJECT_LIBRARY=aflowlib::LIB2RAW_CheckProjectFromDirectory(aurostd::execute2string("pwd"));
    cerr << "aflowlib::LIB2RAW_ALL FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;

    int multi_sh_value=XHOST.CPU_Cores;
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")) multi_sh_value=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"));

    cerr << "multi_sh_value=" << multi_sh_value << endl;
    //    exit(0);
    deque<string> dcmds;

    if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8) ||
       PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)) {
      string command="find "+PROJECT_LIBRARY+"/LIB/ -name \"" + _AFLOWIN_ + "\" | sort ";
      string directory_list=aurostd::execute2string(command);
      vector<string> tokens;
      aurostd::string2tokens(directory_list,tokens,"\n");
      for(uint i=0;i<tokens.size();i++) {
	aurostd::StringSubst(tokens.at(i),"/"+_AFLOWIN_,"");
	aurostd::StringSubst(tokens.at(i),"/"+_AFLOWLOCK_,"");
	aurostd::StringSubst(tokens.at(i),"/core","");
	string cmd="aflow";
	if(XHOST.vflag_control.flag("BEEP")) cmd+=" --beep";
	if(flag_FORCE) cmd+=" --force";
	cmd+=" --lib2raw="+ tokens.at(i);
	dcmds.push_back(cmd);
      };
      if(multi_sh_value==0 || multi_sh_value==1) {
	long double delta_seconds,eta_seconds,reference_seconds=aurostd::get_seconds();
	vector<long double> vdelta_seconds;
	for(uint i=0;i<dcmds.size();i++) {
	  //  cout << delta_seconds << " " << dcmds.at(i) << endl;
	  aflowlib::LIB2RAW(tokens.at(i),flag_FORCE);
	  delta_seconds=aurostd::get_delta_seconds(reference_seconds);
	  if(delta_seconds>1.0) {
	    vdelta_seconds.push_back(delta_seconds);
	    eta_seconds=aurostd::mean(vdelta_seconds)*(dcmds.size()-vdelta_seconds.size());
	    cout << "aflowlib::LIB2RAW_ALL: [STEP]"
		 << "  DONE= " << vdelta_seconds.size() << " / " << dcmds.size()-vdelta_seconds.size() << " (" << 100*vdelta_seconds.size()/dcmds.size()  << ")"
		 << "  iSEC=" << vdelta_seconds.at(vdelta_seconds.size()-1)
		 << "  aSEC=" << aurostd::mean(vdelta_seconds)
		 << "  TO_DO=" << dcmds.size()-vdelta_seconds.size()
	      //	     << "  ETA(secs)=" << eta_seconds
	      //	     << "  ETA(mins)=" << eta_seconds/(60.0)
	      //	     << "  ETA(hours)=" << eta_seconds/(60*60)
		 << "  ETA(days)=" << eta_seconds/(60*60*24)
	      //	     << "  ETA(weeks)=" << eta_seconds/(60*60*24*7)
	      //	     << "  ETA(years)=" << eta_seconds/(60*60*24*365)
		 << endl;
	  }
	}
      } else {
	// DO MULTI
	aurostd::multithread_execute(dcmds,multi_sh_value,FALSE);
      }
      return TRUE;
    }
    cerr << "aflowlib::LIB2RAW_ALL Project Not Found" << endl;exit(0);
    return FALSE;
  }
}

// ***************************************************************************
// aflowlib::GetSpeciesDirectory
// ***************************************************************************
namespace aflowlib {
  uint GetSpeciesDirectory(string directory,vector<string>& vspecies) {
    vspecies.clear();vector<string> vs,tokens;
    stringstream oss;
    string temp_file=aurostd::TmpFileCreate("getspecies");
    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
    deque<string> vcmd; aurostd::string2tokens("cat,bzcat,xzcat,gzcat",vcmd,",");
    if(vext.size()!=vcmd.size()) { cerr << "ERROR - aflowlib::GetSpeciesDirectory: vext.size()!=vcmd.size()" << endl;exit(0); }

    for(uint iext=0;iext<vext.size();iext++) { // check _AFLOWIN_.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/" + _AFLOWIN_ + vext.at(iext))) {
	aurostd::string2vectorstring(aurostd::execute2string(vcmd.at(iext)+" "+directory+"/" + _AFLOWIN_ + vext.at(iext)+" | grep VASP_POTCAR_FILE"),vspecies);
	for(uint i=0;i<vspecies.size();i++) aurostd::StringSubst(vspecies.at(i),"[VASP_POTCAR_FILE]","");
	for(uint i=0;i<vspecies.size();i++) KBIN::VASP_PseudoPotential_CleanName(vspecies.at(i));
      }
    }
    for(uint iext=0;iext<vext.size();iext++) { // check POSCAR.orig.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/POSCAR.orig"+vext.at(iext))) {
	oss.str(aurostd::execute2string(vcmd.at(iext)+" "+directory+"/POSCAR.orig"+vext.at(iext)));
	xstructure xstr(oss,IOVASP_POSCAR);
	if(xstr.species.size()>0) if(xstr.species.at(0)!="") for(uint i=0;i<xstr.species.size();i++) vspecies.push_back(xstr.species.at(i)); // dont change order
      }
    }
    for(uint iext=0;iext<vext.size();iext++) {  // check POSCAR.relax1.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/POSCAR.relax1"+vext.at(iext))) {
	oss.str(aurostd::execute2string(vcmd.at(iext)+" "+directory+"/POSCAR.relax1"+vext.at(iext)));
	xstructure xstr(oss,IOVASP_POSCAR);
	if(xstr.species.size()>0) if(xstr.species.at(0)!="") for(uint i=0;i<xstr.species.size();i++) vspecies.push_back(xstr.species.at(i)); // dont change order
      }
    }
    for(uint iext=0;iext<vext.size();iext++) { // check POSCAR.bands.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/POSCAR.bands"+vext.at(iext))) {
	oss.str(aurostd::execute2string(vcmd.at(iext)+" "+directory+"/POSCAR.bands"+vext.at(iext)));
	xstructure xstr(oss,IOVASP_POSCAR);
	if(xstr.species.size()>0) if(xstr.species.at(0)!="") for(uint i=0;i<xstr.species.size();i++) vspecies.push_back(xstr.species.at(i)); // dont change order
      }
    }
    for(uint iext=0;iext<vext.size();iext++) { // check OUTCAR.relax1.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/OUTCAR.relax1"+vext.at(iext))) {
	aurostd::string2vectorstring(aurostd::execute2string(vcmd.at(iext)+" "+directory+"/OUTCAR.relax1"+vext.at(iext)+" | grep TITEL"),vs);
	for(uint i=0;i<vs.size();i++) { aurostd::string2tokens(vs.at(i),tokens," ");vspecies.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(3))); }
      }
    }
    for(uint iext=0;iext<vext.size();iext++) { // check OUTCAR.static.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/OUTCAR.static"+vext.at(iext))) {
	aurostd::string2vectorstring(aurostd::execute2string(vcmd.at(iext)+" "+directory+"/OUTCAR.static"+vext.at(iext)+" | grep TITEL"),vs);
	for(uint i=0;i<vs.size();i++) { aurostd::string2tokens(vs.at(i),tokens," ");vspecies.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(3))); }
      }
    }
    //  for(uint i=0;i<vspecies.size();i++) cerr << i << " - " << vspecies.at(i) << endl; exit(0);
    return vspecies.size();
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW(string options,bool flag_FORCE,bool LOCAL) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW: BEGIN" << endl;
    if(LDEBUG) cerr << "aflowlib::LIB2RAW: options=" << options << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(LDEBUG) cerr << "aflowlib::LIB2RAW: tokens.size()=" << tokens.size() << endl;
    if(tokens.size()==0 || tokens.size()>2) {
      init::ErrorOption(cout,options,"aflowlib::LIB2RAW",aurostd::liststring2string("aflow --lib2raw=directory","aflow --lib2raw=all[,dir]"));
      exit(0);
    }

    if(tokens.size()>=1) {
      if(tokens.at(0)=="all") {
	XHOST.sensors_allowed=FALSE;
	XHOST.vflag_pflow.flag("MULTI=SH",FALSE);
	return LIB2RAW_ALL(options,flag_FORCE);
	XHOST.sensors_allowed=TRUE;
      }
    }

    string directory_LIB,directory_RAW,directory_WEB;
    bool flag_WEB=FALSE;
    bool flag_files_LIB=FALSE,flag_files_RAW=FALSE,flag_files_WEB=FALSE;
    string PROJECT_LIBRARY;
    if(LOCAL) {
      flag_FORCE=true;
      string directory=aurostd::CleanFileName(options);
      aurostd::StringSubst(directory,"./","");
      if(directory=="." || directory.empty()) { directory=aurostd::execute2string("pwd"); }
      flag_WEB=FALSE;
      flag_files_LIB=FALSE,flag_files_RAW=FALSE,flag_files_WEB=FALSE;
      directory_LIB=directory;
      directory_RAW=directory+"/RAW";
      directory_WEB=directory+"/WEB";
      PROJECT_LIBRARY=directory_LIB;
    } else {    //normal run
      //  cout << "aflowlib::LIB2RAW: AFLOW (" << VERSION << ")" << endl;
      CheckMaterialServer("aflowlib::LIB2RAW"); // must be in AFLOW_MATERIALS_SERVER
      string directory=aurostd::CleanFileName(options);
      if(directory.at(directory.size()-1)=='/')  directory=directory.substr(0,directory.size()-1);
      aurostd::StringSubst(directory,"/"+_AFLOWIN_,"");
      aurostd::StringSubst(directory,"/"+_AFLOWLOCK_,"");
      aurostd::StringSubst(directory,"/core","");
      aurostd::StringSubst(directory,"common/SCINT","common/ICSD");   // [HISTORIC]
      aurostd::StringSubst(directory,"common/ELPASOLITES","common/AURO");   // [HISTORIC]
      aurostd::StringSubst(directory,"LIBRARYX/RAW","LIBRARYX/LIB");   // [HISTORIC]
      aurostd::StringSubst(directory,"LIB2/RAW","LIB2/LIB");

      PROJECT_LIBRARY=aflowlib::LIB2RAW_CheckProjectFromDirectory(directory);
      // cout << "aflowlib::LIB2RAW: FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
      cout << "aflowlib::LIB2RAW: directory=" << directory << endl;

      //  bool flag_LIB=FALSE;
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: scan libraries BEGIN [1]" << endl;
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: PROJECT_LIBRARY=" << PROJECT_LIBRARY << endl;
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB1=" << XHOST_LIBRARY_LIB1 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)) { flag_WEB=FALSE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB2=" << XHOST_LIBRARY_LIB2 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)) { flag_WEB=FALSE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB3=" << XHOST_LIBRARY_LIB3 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB4=" << XHOST_LIBRARY_LIB4 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB5=" << XHOST_LIBRARY_LIB5 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB6=" << XHOST_LIBRARY_LIB6 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB7=" << XHOST_LIBRARY_LIB7 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB8=" << XHOST_LIBRARY_LIB8 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_LIB9=" << XHOST_LIBRARY_LIB9 << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: XHOST_LIBRARY_ICSD=" << XHOST_LIBRARY_ICSD << endl;
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)) { flag_WEB=TRUE;flag_files_WEB=TRUE; }
      if(LDEBUG) cerr << "aflowlib::LIB2RAW: scan libraries END [2]"  << endl;

      stringstream aus_exec;
      //   vector<string> tokens;
      aurostd::string2tokens(directory,tokens,"/");
      directory_LIB.clear();directory_RAW.clear();
      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)) {
	for(uint i=tokens.size()-1;i>0;i--)
	  if(aurostd::substring2bool(tokens.at(i),"_ICSD_")) {
	    if(i>=1) {
	      directory_LIB=tokens.at(i-1)+"/"+tokens.at(i);
	      //    directory_WEB=tokens.at(i-1);  removed to make things consistent
	      directory_WEB=directory;
	    }
	  }
	if(directory_LIB.length()==0) {
	  cout << "aflowlib::LIB2RAW: FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
	  cout << "aflowlib::LIB2RAW: you must specify the directory including the whole lattice type" << endl;
	  cout << " such as  aflow --lib2raw=FCC/La1Se1_ICSD_27104  " << endl;
	  exit(0);
	}
      }

      // strip the directory_LIB of everything else
      directory_LIB=directory; // somewhere to start
      if(aurostd::substring2bool(directory_LIB,"LIB/")) directory_LIB=aurostd::substring2string(directory_LIB,"LIB/",FALSE);
      directory_RAW=directory_LIB;
      directory_WEB=directory_LIB;
      directory_LIB=aurostd::CleanFileName(PROJECT_LIBRARY+"/LIB/"+directory_LIB);
      directory_RAW=aurostd::CleanFileName(PROJECT_LIBRARY+"/RAW/"+directory_RAW);aurostd::StringSubst(directory_RAW,"RAW/LIB","RAW");
      directory_WEB=aurostd::CleanFileName(PROJECT_LIBRARY+"/WEB/"+directory_WEB);aurostd::StringSubst(directory_WEB,"WEB/LIB","WEB");
      if(flag_WEB==FALSE) directory_WEB=aurostd::CleanFileName(PROJECT_LIBRARY+"/WEB/");

      if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1) ||
	 PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)) {
	//      directory_WEB=directory_RAW;
      }
      cout << "aflowlib::LIB2RAW: PROJECT_LIBRARY=" << PROJECT_LIBRARY << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9) << endl;
      // cout << "aflowlib::LIB2RAW: vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD) << endl;
    }
    cout << "aflowlib::LIB2RAW: directory_LIB=" << directory_LIB << endl;
    cout << "aflowlib::LIB2RAW: directory_RAW=" << directory_RAW << endl;
    cout << "aflowlib::LIB2RAW: directory_WEB=" << directory_WEB << endl;
    
    bool perform_LOCK=TRUE,perform_BANDS=FALSE,perform_BADER=FALSE,perform_MAGNETIC=FALSE,perform_THERMODYNAMICS=FALSE;
    bool perform_AGL=FALSE,perform_AEL=FALSE;

    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
    for(uint iext=0;iext<vext.size();iext++) {
      if(aurostd::FileExist(directory_LIB+"/OUTCAR.relax1"+vext.at(iext)) ||
	 aurostd::FileExist(directory_LIB+"/OUTCAR.relax2"+vext.at(iext)) ||
	 aurostd::FileExist(directory_LIB+"/OUTCAR.relax3"+vext.at(iext))) perform_THERMODYNAMICS=TRUE;
      if(aurostd::FileExist(directory_LIB+"/OUTCAR.static"+vext.at(iext))) perform_MAGNETIC=TRUE;
      if(aurostd::FileExist(directory_LIB+"/OUTCAR.bands"+vext.at(iext))) perform_BANDS=TRUE;
      if(aurostd::FileExist(directory_LIB+"/AECCAR0.static"+vext.at(iext)) && aurostd::FileExist(directory_LIB+"/AECCAR2.static"+vext.at(iext))) perform_BADER=TRUE;
      if(aurostd::FileExist(directory_LIB+"/aflow.agl.out"+vext.at(iext))) perform_AGL=TRUE;
      if(aurostd::FileExist(directory_LIB+"/aflow.ael.out"+vext.at(iext))) perform_AEL=TRUE;
    }
 
    //    directory_RAW="/home/junkai/RAW";
    // override for the web
    //    if(PROJECT_LIBRARY==vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)) directory_WEB=aurostd::CleanFileName(PROJECT_AFLOWLIB_WEB);

    // cerr << directory_WEB << endl; cerr << directory_LIB << endl; cerr << directory_RAW << endl; exit(0);
    if(!aurostd::FileExist(directory_LIB+"/"+_AFLOWIN_)) {
      cout << "aflowlib::LIB2RAW: FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
      cout << "aflowlib::LIB2RAW: directory does not exist: " << directory_LIB << endl;
      return FALSE;
      // exit(0);
    }
    if(flag_FORCE==FALSE) {
      if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {  // directory_RAW+"/"+_AFLOWIN_)
	_aflags aflags;
 	aurostd::ZIP2ZIP(directory_LIB,"bz2","xz",FALSE); // PATCH FOR REFRESH (to be removed)
	aurostd::ZIP2ZIP(directory_LIB,"gz","xz",FALSE); // PATCH FOR REFRESH (to be removed)
	cout << "aflowlib::LIB2RAW: ALREADY CALCULATED = " << directory_RAW << "   END_DATE - [v=" << string(AFLOW_VERSION) << "] -" << Message(aflags,"time") << endl;
	return FALSE;
      }
      if(perform_BANDS) {
	if(aurostd::FileExist(directory_RAW+"/EIGENVAL.bands") || aurostd::EFileExist(directory_RAW+"/EIGENVAL.bands")) {
	  // return FALSE;
	  cout << "aflowlib::LIB2RAW: directory is skip because of BANDS: " << directory_RAW << endl;
	  return FALSE;
	  //  exit(0);
	}
      }
    }
    // directory_LIB exists and directory_RAW does not exist can move on:
    { aurostd::DirectoryMake(directory_RAW);aurostd::execute("rm -f /"+directory_RAW+"/*"); }
    if(!aurostd::IsDirectory(directory_RAW)) {
      cout << "aflowlib::LIB2RAW: directory is skip because cannot create directory_RAW: " << directory_RAW << endl;
      return FALSE;
    }
    if(flag_WEB) {
      aurostd::DirectoryMake(directory_WEB);aurostd::execute("rm -f /"+directory_WEB+"/*");
      if(!aurostd::IsDirectory(directory_WEB)) {
	cout << "aflowlib::LIB2RAW: directory is skip because cannot create directory_WEB: " << directory_WEB << endl;
	return FALSE;
      }
    }
    cout << "aflowlib::LIB2RAW: FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
    if((perform_THERMODYNAMICS || perform_BANDS ||  perform_MAGNETIC)) {
      _aflags aflags;
      cout << "aflowlib::LIB2RAW: dir=" << directory_LIB << "   BEGIN_DATE = " << Message(aflags,"user,host,time") << endl;
      aurostd::ZIP2ZIP(directory_LIB,"bz2","xz");
      aurostd::ZIP2ZIP(directory_LIB,"gz","xz");
         
      vector<string> vfile;   // the needed files
      aflowlib::_aflowlib_entry aflowlib_data;
      vector<string> vspecies;   // the species
      GetSpeciesDirectory(directory_LIB,vspecies);
      for(uint i=0;i<vspecies.size();i++) { aflowlib_data.species+=vspecies.at(i);if(i<vspecies.size()-1) aflowlib_data.species+=","; }    // plug vspecies into aflowlib_data

      if(LOCAL) {
	aflowlib_data.aurl=aflowlib_data.auid=directory_LIB; //dummy
      } else {
	// build aflowlib_data.aurl
	aflowlib_data.aurl="aflowlib.duke.edu:"+directory_LIB;
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIBRARYX")) { aurostd::StringSubst(aflowlib_data.aurl,"common/GNDSTATE/LIBRARYX/LIB","AFLOWDATA/LIB2_RAW");aflowlib_data.catalog="LIBRARYX"; } // [HISTORIC]
	if(aurostd::substring2bool(aflowlib_data.aurl,"ICSD")) { aurostd::StringSubst(aflowlib_data.aurl,"common/ICSD/LIB","AFLOWDATA/ICSD_WEB");aflowlib_data.catalog="ICSD"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB1")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB1/LIB","AFLOWDATA/LIB1_RAW");aflowlib_data.catalog="LIB1"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB2")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB2/LIB","AFLOWDATA/LIB2_RAW");aflowlib_data.catalog="LIB2"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB3")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB3/LIB","AFLOWDATA/LIB3_RAW");aflowlib_data.catalog="LIB3"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB4")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB4/LIB","AFLOWDATA/LIB4_RAW");aflowlib_data.catalog="LIB4"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB5")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB5/LIB","AFLOWDATA/LIB5_RAW");aflowlib_data.catalog="LIB5"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB6")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB6/LIB","AFLOWDATA/LIB6_RAW");aflowlib_data.catalog="LIB6"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB7")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB7/LIB","AFLOWDATA/LIB7_RAW");aflowlib_data.catalog="LIB7"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB8")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB8/LIB","AFLOWDATA/LIB8_RAW");aflowlib_data.catalog="LIB8"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"LIB9")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB9/LIB","AFLOWDATA/LIB9_RAW");aflowlib_data.catalog="LIB9"; }
	if(aurostd::substring2bool(aflowlib_data.aurl,"AURO")) { aurostd::StringSubst(aflowlib_data.aurl,"common/AURO/LIB","AFLOWDATA/AURO_RAW");aflowlib_data.catalog="AURO"; }
	aurostd::StringSubst(aflowlib_data.aurl,":/AFLOWDATA",":AFLOWDATA");
	cout << "aflowlib::LIB2RAW: AURL = " << aurostd::PaddedPOST(aflowlib_data.aurl,60) << endl;//"   " << directory_LIB << endl; 
	// build aflowlib_data.auid
	if(LDEBUG) cerr << "aflowlib::LIB2RAW: [AUID=0] directory_LIB=" << directory_LIB << endl;
	aflowlib_data.auid=aflowlib::directory2auid(directory_LIB,aflowlib_data.aurl);
	//  if(AFLOWLIB_VERBOSE)
	cout << "aflowlib::LIB2RAW: AUID = " << aurostd::PaddedPOST(aflowlib_data.auid,60) << endl;//"   " << directory_LIB << endl;
	if(LDEBUG) cerr << "aflowlib::LIB2RAW: [AUID=2]" << endl;
	cout << "aflowlib::LIB2RAW: CATALOG = " << aurostd::PaddedPOST(aflowlib_data.catalog,60) << endl;//"   " << directory_LIB << endl;
	if(LDEBUG) cerr << "aflowlib::LIB2RAW: [AUID=3]" << endl;
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // do the THERMNODYNAMYCS
      if(perform_THERMODYNAMICS) {
	cout << "aflowlib::LIB2RAW: THERMODYNAMIC LOOP ---------------------------------------------------------------------------------" << endl;//perform_BANDS=FALSE;
	aflowlib::LIB2RAW_Loop_Thermodynamics(directory_LIB,directory_RAW,vfile,aflowlib_data,"aflowlib::LIB2RAW (thermodynamics):",LOCAL);	// identifier inside
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // do the BANDS
      if(perform_BANDS) {
	cout << "aflowlib::LIB2RAW: BANDS LOOP ---------------------------------------------------------------------------------" << endl;
	aflowlib::LIB2RAW_Loop_Bands(directory_LIB,directory_RAW,vfile,aflowlib_data,"aflowlib::LIB2RAW (bands):");
	// MOVE/LINK PICS data
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // do the MAGNETIC
      if((perform_MAGNETIC || perform_BANDS)) { // JUNKAI
	cout << "aflowlib::LIB2RAW: MAGNETIC LOOP ---------------------------------------------------------------------------------" << endl;
	aflowlib::LIB2RAW_Loop_Magnetic(directory_LIB,directory_RAW,vfile,aflowlib_data,"aflowlib::LIB2RAW (magnetic):");	
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // do the BADER
      //      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB6")) perform_BADER=FALSE; // hack
      if(perform_BADER) {
	cout << "aflowlib::LIB2RAW: BADER LOOP ---------------------------------------------------------------------------------" << endl;
	aflowlib::LIB2RAW_Loop_Bader(directory_LIB,directory_RAW,vfile,aflowlib_data,"aflowlib::LIB2RAW (bader):");
	// MOVE/LINK PICS data
	// [OBSOLETE] [MOVED DOWN] if(flag_WEB) {
	// [OBSOLETE] [MOVED DOWN] aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/*_abader.out")+" "+directory_WEB);      // LINK
	// [OBSOLETE] [MOVED DOWN] aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/*jvxl")+" "+directory_WEB);            // LINK
	// [OBSOLETE] [MOVED DOWN] }
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // do the AGL
      if(perform_AGL) {
	cout << "aflowlib::LIB2RAW: AGL LOOP ---------------------------------------------------------------------------------" << endl;
	aflowlib::LIB2RAW_Loop_AGL(directory_LIB,directory_RAW,vfile,aflowlib_data,"aflowlib::LIB2RAW (agl):");
	if(flag_WEB) {
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/aflow.agl.out")+" "+directory_WEB);    // LINK
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/AGL.out")+" "+directory_WEB);    // LINK
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/AGL_energies_temperature.out")+" "+directory_WEB);    // LINK
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/AGL_thermal_properties_temperature.out")+" "+directory_WEB);    // LINK
	}
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // do the AEL
      if(perform_AEL) {
	cout << "aflowlib::LIB2RAW: AEL LOOP ---------------------------------------------------------------------------------" << endl;
	aflowlib::LIB2RAW_Loop_AEL(directory_LIB,directory_RAW,vfile,aflowlib_data,"aflowlib::LIB2RAW (ael):");
	if(flag_WEB) {
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/aflow.ael.out")+" "+directory_WEB);    // LINK
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/AEL_Elastic_constants.out")+" "+directory_WEB);    // LINK
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/AEL_Compliance_tensor.out")+" "+directory_WEB);    // LINK
	}
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // do the LOCK
      if(perform_LOCK) {
	cout << "aflowlib::LIB2RAW: LOCK LOOP ---------------------------------------------------------------------------------" << endl;
	aflowlib::LIB2RAW_Loop_LOCK(directory_LIB,directory_RAW,vfile,aflowlib_data,"aflowlib::LIB2RAW (LOCK):");	
      }
      // ---------------------------------------------------------------------------------------------------------------------------------
      // write DOS + BANDS JSON // CO 171025
      string system_name=KBIN::ExtractSystemName(directory_LIB);
      if((aurostd::FileExist(directory_LIB+"/DOSCAR.static") || aurostd::EFileExist(directory_LIB+"/DOSCAR.static")) &&
	 (aurostd::FileExist(directory_LIB+"/POSCAR.static") || aurostd::EFileExist(directory_LIB+"/POSCAR.static"))) {
	stringstream _json_file,json_file;
	aurostd::xoption vpflow;  //dummy
	if(estructure::DOSDATA_JSON(vpflow,directory_LIB,_json_file,true)) {
	  json_file << _json_file.str() << endl; _json_file.str("");
	  cout << "aflowlib::LIB2RAW: compressing file: " << string(directory_RAW+"/"+system_name+"_dosdata.json") << endl; cout.flush();
	  aurostd::stringstream2compressfile(DEFAULT_KZIP_BIN,json_file,directory_RAW+"/"+system_name+"_dosdata.json");
	  json_file.str("");
	  if(aurostd::EFileExist(directory_LIB+"/EIGENVAL.bands") && aurostd::EFileExist(directory_LIB+"/KPOINTS.bands")) {
	    if(estructure::BANDSDATA_JSON(vpflow,directory_LIB,_json_file,true)) {
	      json_file << _json_file.str() << endl; _json_file.str("");
	      cout << "aflowlib::LIB2RAW: compressing file: " << string(directory_RAW+"/"+system_name+"_bandsdata.json") << endl; cout.flush();
	      aurostd::stringstream2compressfile(DEFAULT_KZIP_BIN,json_file,directory_RAW+"/"+system_name+"_bandsdata.json");
	      json_file.str("");
	    }
	  }
	}
      }
      
      // ---------------------------------------------------------------------------------------------------------------------------------
      // DO THE COMPRESSING
      //      cout << "COMPRESSING" << endl;

      // generic for compressing and linking
      deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
      deque<string> vtype; aurostd::string2tokens(".orig,.relax,.relax1,.relax2,.relax3,.relax4,.static,.bands",vtype,",");
      deque<string> vout; aurostd::string2tokens(".out,.json",vout,",");

      for(uint iext=0;iext<vext.size();iext++) {
	aurostd::RemoveFile(directory_RAW+"/*.eps"+vext.at(iext));
      }
       
      // COREY DESTROYED THE LOGICS OF THIS ROUTINE AND I HAD TO REWIRE...
      // [OBSOLETE] for(uint iext=0;iext<vext.size();iext++) {
      // [OBSOLETE] 	//  aurostd::RemoveFile(directory_RAW+"/core*");
      // [OBSOLETE] 	  for(uint ifile=0;ifile<vfile.size();ifile++) {
      // [OBSOLETE] 	  cout << vfile.at(ifile) << endl;
      // [OBSOLETE] 	  if(!aurostd::substring2bool(vfile.at(ifile),vext.at(iext)) && aurostd::FileExist(directory_RAW+"/"+vfile.at(ifile))) {
      // [OBSOLETE] 	    if(LOCAL) { //CO 171025
      // [OBSOLETE] 	      if(aurostd::EFileExist(directory_LIB+"/"+vfile.at(ifile))) { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); } //it's all sitting in the directory above, wasteful
      // [OBSOLETE] 	      if(vfile.at(ifile)=="KPOINTS.relax") { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); }
      // [OBSOLETE] 	      if(vfile.at(ifile)=="EIGENVAL.bands.old") { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); }
      // [OBSOLETE] 	      if(vfile.at(ifile)=="OUTCAR.relax") { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); }
      // [OBSOLETE] 	    } 
      // [OBSOLETE] 	} // vfile
      // [OBSOLETE]  } // iext
      
      // FILES to remove
      deque<string> vfile2remove; aurostd::string2tokens("KPOINTS.bands.old,EIGENVAL.bands.old,OUTCAR.relax1,OUTCAR.bands",vfile2remove,",");
      vfile2remove.push_back("aflow.pgroupk_xtal.out"); // comes from nowere (DAVID) DX
      for(uint iremove=0;iremove<vfile2remove.size();iremove++) {
 	if(aurostd::FileExist(directory_RAW+"/"+vfile2remove.at(iremove))) { // need to be present
	  //if(LDEBUG)
	  cout << "aflowlib::LIB2RAW: removing file: " << string(directory_RAW+"/"+vfile2remove.at(iremove)) << endl;
	  aurostd::RemoveFile(directory_RAW+"/"+vfile2remove.at(iremove));
	} // FileExist
      } // iremove
      // FILES to compress if not compressed already (linked), in such case they will be deleted.
      deque<string> vfile2compress0; aurostd::string2tokens("OUTCAR.relax",vfile2compress0,",");
      for(uint icompress=0;icompress<vfile2compress0.size();icompress++) {
	if(aurostd::FileExist(directory_RAW+"/"+vfile2compress0.at(icompress)+"."+DEFAULT_KZIP_BIN)) { // test if there is one already compressed (possibly linked)
	  if(1||LDEBUG) { cout << "aflowlib::LIB2RAW: found compressed file: " << string(directory_RAW+"/"+vfile2compress0.at(icompress)+"."+DEFAULT_KZIP_BIN) << endl;  cout.flush(); }
	  if(1||LDEBUG) { cout << "aflowlib::LIB2RAW: removing file: " << string(directory_RAW+"/"+vfile2compress0.at(icompress)) << endl;  cout.flush(); }
	  aurostd::RemoveFile(directory_RAW+"/"+vfile2compress0.at(icompress));
	}
	if(aurostd::FileExist(directory_RAW+"/"+vfile2compress0.at(icompress))) { // need to be present
	  if(1||LDEBUG) { cout << "aflowlib::LIB2RAW: compressing file: " << string(directory_RAW+"/"+vfile2compress0.at(icompress)) << endl;  cout.flush(); }
	  aurostd::CompressFile(directory_RAW+"/"+vfile2compress0.at(icompress),DEFAULT_KZIP_BIN);
	} // FILES to compress
      } // icompress
      if(perform_BADER) {
	if(LDEBUG) { cout << "aflowlib::LIB2RAW: compressing file: " << string(directory_RAW+"/"+"*.jvxl") << endl;  cout.flush(); }
	aurostd::CompressFile(directory_RAW+"/"+"*.jvxl",DEFAULT_KZIP_BIN);
      }
      // FILES to compress
      deque<string> vfile2compress1; aurostd::string2tokens("aflow.pgroup,aflow.pgroup_xtal,aflow.pgroupk,aflow.pgroupk_xtal,aflow.fgroup,aflow.iatoms,aflow.agroup",vfile2compress1,",");
      for(uint ilink=0;ilink<vfile2compress1.size();ilink++) {
	for(uint iout=0;iout<vout.size();iout++) {
	  for(uint itype=0;itype<vtype.size();itype++) {	  
	    if(aurostd::FileExist(directory_RAW+"/"+vfile2compress1.at(ilink)+vtype.at(itype)+vout.at(iout))) {
	      if(1||LDEBUG) { cout << "aflowlib::LIB2RAW: compressing file (" << ilink << ","<< iout << "," << itype << "): " << string(directory_RAW+"/"+vfile2compress1.at(ilink)+vtype.at(itype)+vout.at(iout)) << endl; cout.flush(); }
	      aurostd::CompressFile(directory_RAW+"/"+vfile2compress1.at(ilink)+vtype.at(itype)+vout.at(iout),DEFAULT_KZIP_BIN);
	    } // FileExist
	  } // itype
	} // iout
      } // ilink
       
      // FILES to link LIB to RAW
      for(uint ifile=0;ifile<vfile.size();ifile++) {
	if(aurostd::FileExist(directory_RAW+"/"+vfile.at(ifile))) { // need to be present
	  deque<string> vfile2link0; aurostd::string2tokens("DOSCAR.static,OUTCAR.static,CHGCAR.static,AECCAR0.static,AECCAR2.static,EIGENVAL.bands,OSZICAR.bands",vfile2link0,",");
	  for(uint ilink=0;ilink<vfile2link0.size();ilink++) {
	    if(vfile.at(ifile)==vfile2link0.at(ilink)) {
	      // cout << "aflowlib::LIB2RAW: linking file RAW->LIB: " << vfile2link0.at(ilink) << endl;	      
	      if(aurostd::FileExist(directory_LIB+"/"+vfile.at(ifile)) || aurostd::EFileExist(directory_LIB+"/"+vfile.at(ifile))) { // need to be present in LIB also		
		aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); // remove RAW original
		if(aurostd::FileExist(directory_LIB+"/"+vfile.at(ifile))) { 
		  cout << "aflowlib::LIB2RAW: linking file RAW->LIB: " << string(directory_LIB+"/"+vfile.at(ifile)) << endl; cout.flush();   
		  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_LIB+"/"+vfile.at(ifile))+" "+directory_RAW); // link LIB to RAW (save space)
		}
		for(uint iext=0;iext<vext.size();iext++) {
		  if(aurostd::FileExist(directory_LIB+"/"+vfile.at(ifile)+vext.at(iext))) {
		    cout << "aflowlib::LIB2RAW: linking file RAW->LIB: " << string(directory_LIB+"/"+vfile.at(ifile)+vext.at(iext)) << endl;  cout.flush();   	      
		    aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_LIB+"/"+vfile.at(ifile)+vext.at(iext))+" "+directory_RAW); // link LIB to RAW (save space)
		  }
		}
	      }
	    }
	  } // files to link LIB to RAW
	} // File Exist
      } // ifile

      // FILES to leave as is
      for(uint ifile=0;ifile<vfile.size();ifile++) {
	if(aurostd::FileExist(directory_RAW+"/"+vfile.at(ifile))) { // need to be present
	  // [NO COMPRESS] if(vfile.at(ifile)=="KPOINTS.relax") aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/"+vfile.at(ifile)); // seems to work, although I do not like if(SC-0914)
	  // [NO COMPRESS] if(vfile.at(ifile)=="KPOINTS.static") aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/"+vfile.at(ifile)); // seems to work, although I do not like if(SC-0914)
	  // [NO COMPRESS] if(vfile.at(ifile)=="KPOINTS.bands") aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/"+vfile.at(ifile)); // seems to work, although I do not like if(SC-0914)
	  // [NO COMPRESS] if(vfile.at(ifile)=="INCAR.relax") aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/"+vfile.at(ifile)); // seems to work, although I do not like if(SC-0914)
	  // [NO COMPRESS] if(vfile.at(ifile)=="INCAR.static") aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/"+vfile.at(ifile)); // seems to work, although I do not like if(SC-0914)
	  // [NO COMPRESS] if(vfile.at(ifile)=="INCAR.bands") aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/"+vfile.at(ifile)); // seems to work, although I do not like if(SC-0914)
	  // [NO COMPRESS] if(vfile.at(ifile)=="OUTCAR.relax") aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/"+vfile.at(ifile));  // seems to work, although I do not like if(SC-0512)
	} // File Exist
      } // ifile
      
      // DO THE FINISH LINK/COPY FOR WEB
      
      if(perform_BANDS || 1) {
	cout << "aflowlib::LIB2RAW: linking stuff flag_WEB=" << flag_WEB << endl;
	// MOVE/LINK PICS data

	if(flag_WEB) {
	  string system_name=KBIN::ExtractSystemName(directory_LIB);
	  cout << "aflowlib::LIB2RAW: linking SYSTEM=" << system_name << endl;
	  if(aurostd::FileExist(directory_RAW+"/"+system_name+".png"))
	    aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/*png")+" "+directory_WEB);            // LINK
	  if(aurostd::FileExist(directory_RAW+"/"+system_name+".cif"))
	    aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/*cif")+" "+directory_WEB);            // LINK

	  aurostd::execute("ln -sf /www/AFLOWDATA/api_index.php "+directory_RAW+"/index.php");                      // LINK
	  aurostd::execute("ln -sf /www/AFLOWDATA/api_index.php "+directory_WEB+"/index.php");                      // LINK
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)+" "+directory_WEB);    // LINK
	  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON)+" "+directory_WEB);    // LINK
	  
	  deque<string> vfile2link1; aurostd::string2tokens("aflow.pgroup,aflow.pgroup_xtal,aflow.pgroupk,aflow.pgroupk_xtal,aflow.fgroup,aflow.iatoms,aflow.agroup,edata,data",vfile2link1,",");
	  for(uint ilink=0;ilink<vfile2link1.size();ilink++) {
	    for(uint iout=0;iout<vout.size();iout++) {
	      for(uint itype=0;itype<vtype.size();itype++) {	  
		if(aurostd::FileExist(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout))) { // no compression
		  if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)) << endl; cout.flush(); }	      	      
		  aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout))+" "+directory_WEB);
		}  // FileExist
		for(uint iext=0;iext<vext.size();iext++) {
		  if(aurostd::FileExist(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)+vext.at(iext))) { // with compression
		    if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)+vext.at(iext)) << endl; cout.flush(); }
		    aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)+vext.at(iext))+" "+directory_WEB);
		  } // FileExist
		} // iext
	      } //itype
	    } // iout
	  } // ilink
	  // PREVIOUSLY IN NETFLIX deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
	  for(uint iout=0;iout<vout.size();iout++) {
	    if(aurostd::FileExist(directory_RAW+"/"+system_name+"_structure_relax"+vout.at(iout))) {
	      if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_structure_relax"+vout.at(iout)) << endl; cout.flush(); }
	      aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+system_name+"_structure_relax"+vout.at(iout))+" "+directory_WEB);  // CO 171024
	    }  // FileExist
	  } // iout
	  for(uint iout=0;iout<vout.size();iout++) {
	    if(aurostd::FileExist(directory_RAW+"/"+system_name+"_structure_relax1"+vout.at(iout))) {
	      if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_structure_relax1"+vout.at(iout)) << endl; cout.flush(); }
	      aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+system_name+"_structure_relax1"+vout.at(iout))+" "+directory_WEB);  // CO 171024
	    } // FileExist
	  } // iout
	  for(uint iout=0;iout<vout.size();iout++) {
	    if(aurostd::FileExist(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)))  {  // NO EXTENSION
	      if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)) << endl; cout.flush(); }
	      aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout))+" "+directory_WEB);       // CO 171024
	    }
	    for(uint iext=0;iext<vext.size();iext++) {
	      if(aurostd::FileExist(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)+vext.at(iext)))  {
		if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)+vext.at(iext)) << endl; cout.flush(); }
		aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)+vext.at(iext))+" "+directory_WEB);       // CO 171024
	      } // FileExist
	    } // iext
	  } // iout
	  for(uint iout=0;iout<vout.size();iout++) {
	    if(aurostd::FileExist(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout))) { // NO EXTENSION
	      if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)) << endl; cout.flush(); }
	      aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout))+" "+directory_WEB);     // CO 171024
	    } // FileExist
	    for(uint iext=0;iext<vext.size();iext++) {
	      if(aurostd::FileExist(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)+vext.at(iext))) {
		if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)+vext.at(iext)) << endl; cout.flush(); }
		aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)+vext.at(iext))+" "+directory_WEB);     // CO 171024
	      } // FileExist
	    } // iext
	  } // iout
	  deque<string> vfile2link2; aurostd::string2tokens("EIGENVAL.bands,DOSCAR.static,OUTCAR.static,CONTCAR.relax,CONTCAR.relax1,POSCAR.bands,CONTCAR.relax.vasp,CONTCAR.relax.qe,CONTCAR.relax.abinit,CONTCAR.relax.aims,KPOINTS.relax,KPOINTS.static,KPOINTS.bands,INCAR.bands,AECCAR0.static,AECCAR2.static,CHGCAR.static",vfile2link2,",");
	  for(uint ilink=0;ilink<vfile2link2.size();ilink++) {
	    if(aurostd::FileExist(directory_RAW+"/"+vfile2link2.at(ilink))) { // NO EXTENSION
	      if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link2.at(ilink)) << endl; cout.flush(); }
	      aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+vfile2link2.at(ilink))+" "+directory_WEB);         // LINK
	    }
	    for(uint iext=0;iext<vext.size();iext++) {
	      if(aurostd::FileExist(directory_RAW+"/"+vfile2link2.at(ilink)+vext.at(iext))) {
		if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link2.at(ilink)+vext.at(iext)) << endl; cout.flush(); }
		aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+vfile2link2.at(ilink)+vext.at(iext))+" "+directory_WEB);         // LINK
	      }
	    } // iext
	  } // ilink
	  
	  if(perform_BADER) {
	    if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+"*jvxl") << endl; cout.flush(); }
	    aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+"*jvxl*")+" "+directory_WEB);         // LINK
	    if(LDEBUG) { cout << "aflowlib::LIB2RAW: linking file WEB->RAW: " << string(directory_RAW+"/"+"*_abader.out") << endl; cout.flush(); }
	    aurostd::execute("ln -sf "+aurostd::CleanFileName(directory_RAW+"/"+"*_abader.out*")+" "+directory_WEB);         // LINK
	  }
	} // flag_WEB
      }
      // DONE
      // write files if necessary
      vector<string> vdirectory;
      // do the directories
      if(flag_files_LIB) {
	aurostd::DirectoryLS(directory_LIB,vdirectory);
	for(uint i=0;i<vdirectory.size();i++)
	  aflowlib_data.vfiles.push_back(vdirectory.at(i));
      }
      if(flag_files_RAW) {
	aurostd::DirectoryLS(directory_RAW,vdirectory);
	for(uint i=0;i<vdirectory.size();i++)
	  aflowlib_data.vfiles.push_back(vdirectory.at(i));
      }
      if(flag_files_WEB) {
	aurostd::DirectoryLS(directory_WEB,vdirectory);
	for(uint i=0;i<vdirectory.size();i++)
	  aflowlib_data.vfiles.push_back(vdirectory.at(i));
      }
      // DO THE FINAL WRITING

      //     cout << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ": " << aflowlib_data.aflowlib2file(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
      //      aurostd::execute("ln -sf ../../"+_XENTRY_+" "+directory_RAW+"/"+_XENTRY_);
      if(!LOCAL) { // CO 171025
	aurostd::execute("ln -sf /www/AFLOWDATA/api_index.php "+directory_RAW+"/"+_XENTRY_);
      }
      // write aflowlib.out
      cout << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ": " << aflowlib_data.aflowlib2file(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"out");
      //      cout << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ": " << aflowlib_data.aflowlib2file(directory_WEB+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
      // write aflowlib.json
      cout << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << ": " << aflowlib_data.aflowlib2file(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON,"json");
      //      cout << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << ": " << aflowlib_data.aflowlib2file(directory_WEB+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON);
      if(flag_WEB) {
	if(!LOCAL) { // CO 171025
	  aurostd::execute("ln -sf /www/AFLOWDATA/api_index.php "+directory_WEB+"/"+_XENTRY_);
	}
      }
      // DONE
      cout << "aflowlib::LIB2RAW: dir=" << directory_LIB << "   END_DATE - [v=" << string(AFLOW_VERSION) << "] -" << Message(aflags,"time") << endl;
      if(XHOST.vflag_control.flag("BEEP")) aurostd::beep(aurostd::min(6000,aurostd::abs(int(3*aflowlib_data.aflowlib2string().length()-2000))),50);
    }
    // COMPRESS
    /*
    //  aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/POSCAR.bands");
    //  aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/plotbz.sh");
    */
    // DELETE STUFF

    // CHANGE PERMISSIONS
    // changing order of permission editing, if LOCAL, then order matters, otherwise NOT REALLY
    // files first, since we do /* (just in case directory_RAW is inside directory_LIB)
    // then directories
    // FILES
 
    // CO 180216 - more robust for any type of directory/file setup
    aurostd::execute("chmod 755 `find "+directory_LIB+" -type d`");
    aurostd::execute("chmod 644 `find "+directory_LIB+" -type f`");
    aurostd::execute("chmod 755 `find "+directory_RAW+" -type d`");
    aurostd::execute("chmod 644 `find "+directory_RAW+" -type f`");
    if(CHMODWEB) {
      aurostd::execute("chmod 755 `find "+directory_WEB+" -type d`");
      aurostd::execute("chmod 644 `find "+directory_WEB+" -type f`");
    }

    // ALTERNATIVE but bombs xchmod {}
    // aurostd::execute("find "+directory_LIB+" -type d ! -perm 755 -print -exec xchmod 755 {} \\;");
    // aurostd::execute("find "+directory_LIB+" -type f ! -perm 644 -print -exec xchmod 644 {} \\;");
    // aurostd::execute("find "+directory_RAW+" -type d ! -perm 755 -print -exec xchmod 755 {} \\;");
    // aurostd::execute("find "+directory_RAW+" -type f ! -perm 644 -print -exec xchmod 644 {} \\;");
    // if(CHMODWEB) aurostd::execute("find "+directory_WEB+" -type d ! -perm 755 -print -exec xchmod 755 {} \\;");
    // if(CHMODWEB) aurostd::execute("find "+directory_WEB+" -type f ! -perm 644 -print -exec xchmod 644 {} \\;");
    
    // [OBSOLETE CO 180216] aurostd::ChmodFile("755",directory_LIB);
    // [OBSOLETE CO 180216] aurostd::ChmodFile("644",directory_LIB+"/*");
    // [OBSOLETE CO 180216] aurostd::ChmodFile("755",directory_LIB+"/ARUN*");
    // [OBSOLETE CO 180216] aurostd::ChmodFile("755",directory_RAW);
    // [OBSOLETE CO 180216] aurostd::ChmodFile("644",directory_RAW+"/*");
    // [OBSOLETE CO 180216] aurostd::ChmodFile("755",directory_RAW+"/ARUN*");
    // [OBSOLETE CO 180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("755",directory_WEB); }
    // [OBSOLETE CO 180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("644",directory_WEB+"/*.png"); }
    // [OBSOLETE CO 180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("644",directory_WEB+"/*.html"); }
    // [OBSOLETE CO 180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("644",directory_WEB+"/*.jpg"); }
    // [OBSOLETE CO 180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("755",directory_WEB+"/ARUN*"); }
    // NOW WRITE DOWN THE FILE FOR THE LIBRARY

    // done
    //    cout << "FIX aflowlib::LIB2RAW: [1]" << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_FileNeeded
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_FileNeeded(string directory_LIB,string fileLIB,string directory_RAW,string fileRAW,vector<string> &vfile,string MESSAGE) {
    //  bool LDEBUG=(FALSE || XHOST.DEBUG);

    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
    deque<string> vcmd; aurostd::string2tokens("bzip2,xz,gzip",vcmd,",");
    if(vext.size()!=vcmd.size()) { cerr << "ERROR - aflowlib::LIB2RAW_FileNeeded: vext.size()!=vcmd.size()" << endl;exit(0); }

    string file_LIB=directory_LIB+"/"+fileLIB;
    string file_RAW=directory_RAW+"/"+fileRAW;
    string file_LIB_nocompress=directory_LIB+"/"+fileLIB;
    for(uint iext=0;iext<vext.size();iext++) aurostd::StringSubst(file_LIB_nocompress,vext.at(iext),"");
    string file_RAW_nocompress=directory_RAW+"/"+fileRAW;
    for(uint iext=0;iext<vext.size();iext++) aurostd::StringSubst(file_RAW_nocompress,vext.at(iext),"");
    if(aurostd::FileExist(file_RAW)) return TRUE;   // already there
    if(aurostd::FileExist(file_RAW_nocompress)) return TRUE; // already there

    if(!aurostd::FileExist(file_LIB) && !aurostd::FileExist(file_LIB_nocompress) && !aurostd::EFileExist(file_LIB_nocompress)) {
      if(!aurostd::FileExist(file_LIB)) {
	cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_FileNeeded: file not found " << file_LIB << endl;exit(0); }
      if(!aurostd::FileExist(file_LIB_nocompress)) {
	cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_FileNeeded: file not found " << file_LIB_nocompress << endl;exit(0); }
      if(!aurostd::EFileExist(file_LIB_nocompress)) {
	cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_FileNeeded: file not found " << file_LIB_nocompress << ".EXT" << endl;exit(0); }
    }
    if(aurostd::FileExist(file_LIB)) aurostd::CopyFile(file_LIB,file_RAW);
    if(aurostd::FileExist(file_LIB_nocompress)) aurostd::CopyFile(file_LIB_nocompress,file_RAW_nocompress);
    for(uint iext=0;iext<vext.size();iext++)
      if(aurostd::FileExist(file_LIB_nocompress+vext.at(iext)))
	aurostd::CopyFile(file_LIB_nocompress+vext.at(iext),file_RAW_nocompress+vext.at(iext));

    for(uint iext=0;iext<vext.size();iext++) {
      if(aurostd::FileExist(file_RAW) && aurostd::substring2bool(file_RAW,vext.at(iext))) aurostd::execute(vcmd.at(iext)+" -dqf "+file_RAW);
      if(aurostd::FileExist(file_RAW_nocompress+vext.at(iext))) aurostd::execute(vcmd.at(iext)+" -dqf "+file_RAW_nocompress+vext.at(iext));
    }
    string file2add=fileRAW;
    for(uint iext=0;iext<vext.size();iext++) aurostd::StringSubst(file2add,vext.at(iext),"");
    vfile.push_back(file2add);
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Bands
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_Bands(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // LDEBUG=TRUE;
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Bands [1]" << endl;
    // Stefano Curtarolo 2009-2010-2011-2012
    vector<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");

    cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Bands: begin " << directory_LIB << endl;
    cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Bands: species = " << vspecies.size() << endl;
    data.vloop.push_back("bands");

    stringstream command;command.clear();command.str(std::string());
    stringstream aus_exec;
    // directories must exist already
    bool flag_DATA_BANDS_=FALSE;
    bool flag_use_MATLAB=FALSE,flag_use_GNUPLOT=!flag_use_MATLAB;  // KESONG
    // bool flag_use_MATLAB=TRUE,flag_use_GNUPLOT=!flag_use_MATLAB;   // WAHYU
    // [OBSOLETE]  bool flag_ORIG=FALSE;

    aurostd::StringSubst(directory_LIB,"ELPASOLITES","AURO"); // PATCH
    aurostd::StringSubst(directory_LIB,"SCINT","ICSD"); // PATCH

    // copy _AFLOWIN_ LOCK DOSCAR.static.EXT EIGENVAL.bands.EXT KPOINTS.bands.EXT POSCAR.bands.EXT
    string file_LIB,file_RAW;
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWIN_,directory_RAW,_AFLOWIN_,vfile,MESSAGE);  // _AFLOWIN_
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // LOCK
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"DOSCAR.static",directory_RAW,"DOSCAR.static",vfile,MESSAGE);  // DOSCAR.static
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.static",directory_RAW,"OUTCAR.static",vfile,MESSAGE);  // OUTCAR.static
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OSZICAR.static",directory_RAW,"OSZICAR.static",vfile,MESSAGE);  // OSZICAR.static
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OSZICAR.bands",directory_RAW,"OSZICAR.bands",vfile,MESSAGE);  // OSZICAR.bands
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"EIGENVAL.bands",directory_RAW,"EIGENVAL.bands",vfile,MESSAGE);  // EIGENVAL.bands
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"KPOINTS.static",directory_RAW,"KPOINTS.static",vfile,MESSAGE);  // KPOINTS.static  // not needed but good for show off SC 0914
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"KPOINTS.bands",directory_RAW,"KPOINTS.bands",vfile,MESSAGE);  // KPOINTS.bands
    //  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"INCAR.static",directory_RAW,"INCAR.static",vfile,MESSAGE);  // INCAR.static  // not needed but good for show off SC 0914
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"INCAR.bands",directory_RAW,"INCAR.bands",vfile,MESSAGE);  // INCAR.bands  // not needed but good for show off SC 0914
    //  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1",directory_RAW,"POSCAR.relax1",vfile,MESSAGE);  // POSCAR.relax1, Get Atomic Species Name


    if(TRUE || flag_DATA_BANDS_) { // POSCAR.bands.EXT
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.bands",directory_RAW,"POSCAR.bands",vfile,MESSAGE);  // POSCAR.bands
    }

    xKPOINTS kpoints_static;
    kpoints_static.GetPropertiesFile(directory_RAW+"/KPOINTS.static");
    data.kpoints_nnn_static=kpoints_static.nnn_kpoints;
    data.kpoints+=";"+aurostd::utype2string(kpoints_static.nnn_kpoints[1])+","+aurostd::utype2string(kpoints_static.nnn_kpoints[2])+","+aurostd::utype2string(kpoints_static.nnn_kpoints[3]);
    xKPOINTS kpoints_bands;
    kpoints_bands.GetPropertiesFile(directory_RAW+"/KPOINTS.bands");
    //get pairs
    data.kpoints_pairs.clear();
    if(kpoints_bands.vpath.size()%2==0) {  //if even
      for(uint i=0;i<kpoints_bands.vpath.size();i+=2) {
	data.kpoints_pairs.push_back(kpoints_bands.vpath.at(i)+"-"+kpoints_bands.vpath.at(i+1));
      }
    }
    data.kpoints_bands_path_grid=kpoints_bands.path_grid;
    data.kpoints+=";"+kpoints_bands.path+";"+aurostd::utype2string(kpoints_bands.path_grid);
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " KPOINTS = " << data.kpoints << endl;

    if(flag_use_MATLAB) { // MATLAB STUFF  OLD WAHYU-STEFANO
      // PERFORM THE MATLAB STEP
      cout << MESSAGE << " MATLAB start: " << directory_RAW << endl;
      //WRITE plotbz.sh
      stringstream gnuplot_plotbz;
      gnuplot_plotbz.clear();gnuplot_plotbz.str(std::string());
      // [OBSOLETE]    gnuplot_plotbz << GNUPLOT_FUNCS_plotbz("./","") << endl;
      gnuplot_plotbz << GNUPLOT_FUNCS_plotbz() << endl;
      aurostd::stringstream2file(gnuplot_plotbz,string(directory_RAW+"/plotbz.sh"));

      // WRITE PARAM.M
      stringstream matlab_param;
      matlab_param.clear();matlab_param.str(std::string());
      matlab_param << MATLAB_FUNCS_param() << endl;
      // aurostd::stringstream2file(matlab_param,string(directory_RAW+"/param.m")); // it seems that param is not needed anymore
      // WRITE PLOTBAND.M
      stringstream matlab_plotband;
      matlab_plotband.clear();matlab_plotband.str(std::string());
      // [OBSOLETE]   matlab_plotband << MATLAB_FUNCS_plotband("./","normal") << endl; // normal OR log
      matlab_plotband << MATLAB_FUNCS_plotband() << endl; // normal OR log
      matlab_plotband << "exit;"  << endl;
      aurostd::stringstream2file(matlab_plotband,string(directory_RAW+"/plotband.m"));

      // NEW STUFF
      command.clear();command.str(std::string());
      command << "cd " << directory_RAW << endl;
      command << "mv KPOINTS.bands KPOINTS.bands.old" << endl; vfile.push_back("KPOINTS.bands.old"); // so it is compressed
      command << "mv EIGENVAL.bands EIGENVAL.bands.old" << endl; vfile.push_back("EIGENVAL.bands.old"); // so it is compressed
      aurostd::execute(command);
      command.clear();command.str(std::string());
      _aflags aflags;
      aflags.Directory=directory_RAW;

      if(pflow::FIXBANDS(aflags,"POSCAR.bands,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands,EIGENVAL.bands")==FALSE) {
	cout << "ERROR_RERUN " << directory_LIB << endl;
	return FALSE;
      }

      // EXECUTE PLOTBZ.SH using ksh
      command.clear();command.str(std::string());
      command << "cd " << directory_RAW << endl;
      command << "ksh plotbz.sh" << endl;
      aurostd::execute(command);

      // EXECUTE MATLAB
      command.clear();command.str(std::string());
      command << "cd " << directory_RAW << endl;
      command << "export DISPLAY=:0.0" << endl;
      aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
      command << DEFAULT_KBIN_MATLAB_BIN << " -r " << string("plotband") << endl;
      aurostd::execute(command);
    }

    if(flag_use_GNUPLOT) { // GNUPLOT STUFF NEW KESONG-STEFANO
      // KESONG WRITE THE CODE HERE
      cout << MESSAGE << " GNUPLOT start: " << directory_RAW << endl;
      // WRITE plotbz.sh
      stringstream gnuplot_plotbz;
      gnuplot_plotbz.clear();gnuplot_plotbz.str(std::string());
      // [OBSOLETE]    gnuplot_plotbz << GNUPLOT_FUNCS_plotbz("./","") << endl;
      gnuplot_plotbz << GNUPLOT_FUNCS_plotbz() << endl;
      aurostd::stringstream2file(gnuplot_plotbz,string(directory_RAW+"/plotbz.sh"));

      // NEW STUFF
      command.clear();command.str(std::string());
      command << "cd " << directory_RAW << endl;
      command << "mv KPOINTS.bands KPOINTS.bands.old" << endl; vfile.push_back("KPOINTS.bands.old"); // so it is compressed
      command << "mv EIGENVAL.bands EIGENVAL.bands.old" << endl; vfile.push_back("EIGENVAL.bands.old"); // so it is compressed
      aurostd::execute(command);
      command.clear();command.str(std::string());
      _aflags aflags;
      aflags.Directory=directory_RAW;
      if(pflow::FIXBANDS(aflags,"POSCAR.bands,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands,EIGENVAL.bands")==FALSE) {
	cout << "ERROR_RERUN " << directory_LIB << endl;
	return FALSE;
      }

      // EXECUTE PLOTBZ.SH using ksh
      command.clear();command.str(std::string());
      command << "cd " << directory_RAW << endl;
      command << "ksh plotbz.sh" << endl;
      aurostd::execute(command);

      // EXECUTE PLOTBAND
      char work_dir[1024];
      string cdir; //, wdir;
      cdir = getcwd(work_dir, 1024);  //Get the working directory

      char raw_dir[1024];
      strcpy(raw_dir, directory_RAW.c_str());
      chdir(raw_dir);               //Change into the RAW direcotry

      // [OBSOLETE]  vector<string> directory;
      // [OBSOLETE]  directory.push_back(" ");
      // [OBSOLETE]  directory.push_back(" ");
      // [OBSOLETE]  directory.push_back("./");
      // [OBSOLETE]  estructure::PLOT_BANDDOS(directory);
      // [OBSOLETE]  estructure::PLOT_PEDOSALL_AFLOWLIB(directory, aflags);
      estructure::PLOT_BANDDOS("./");
      estructure::PLOT_PEDOSALL_AFLOWLIB("./", aflags);

      chdir(work_dir);  //Go to the working direcotry
    }

    // Kesong adds it
    command << "cd " << directory_RAW << endl;
    command << "rm -f *pdf *jpg " << endl;
    aurostd::execute(command);

    // DONE
    cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Bands: end " << directory_LIB << endl;
    return TRUE;
  }
}
// ***************************************************************************
// aflowlib::LIB2RAW_Loop_DATA
// ***************************************************************************
// namespace aflowlib {
//   bool LIB2RAW_Loop_DATA(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE) {
//     bool LDEBUG=(FALSE || XHOST.DEBUG);
//     // LDEBUG=TRUE;
//     if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_DATAs [1]" << endl;
//     // Stefano Curtarolo 2009-2010-2011-2012
//     vector<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");

//     cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: begin " << directory_LIB << endl;
//     cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: species = " << vspecies.size() << endl;
//     data.vloop.push_back("data");

//     uint relax_max=10;
//     bool flag_EDATA_ORIG_=FALSE,flag_EDATA_RELAX_=FALSE;
//     bool flag_DATA_ORIG_=FALSE,flag_DATA_RELAX_=FALSE,flag_DATA_BANDS_=FALSE;


//     aurostd::StringSubst(directory_LIB,"ELPASOLITES","AURO"); // PATCH
//     aurostd::StringSubst(directory_LIB,"SCINT","ICSD"); // PATCH
//     if(aurostd::substring2bool(directory_LIB,"ICSD"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB1"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB2"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB3"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB4"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB5"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB6"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB7"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB8"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB9"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"AURO"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIBRARYX")) { flag_EDATA_ORIG_=TRUE; } // [HISTORIC]

//     string file_LIB,file_RAW;
//     aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWIN_,directory_RAW,_AFLOWIN_,vfile,MESSAGE);  // _AFLOWIN_
//     aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // LOCK

//     if(flag_DATA_ORIG_ || flag_EDATA_ORIG_) {  // POSCAR.orig.EXT
//       // if(flag_ORIG==FALSE) {
//       bool found=FALSE;
//       file_LIB=directory_LIB+"/POSCAR.orig"+EXT;
//       if(!found && (found=aurostd::FileExist(file_LIB))) {
// 	cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: building POSCAR.orig from POSCAR.orig"+EXT << endl;
// 	aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.orig",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
//       }
//       file_LIB=directory_LIB+"/POSCAR.relax1"+EXT;
//       if(!found && (found=aurostd::FileExist(file_LIB))) {
// 	cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: building POSCAR.orig from POSCAR.relax1"+EXT << endl;
// 	aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
//       }
//       if(!found) {
// 	found=TRUE;
// 	cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: building POSCAR.orig from " << _AFLOWIN_ << "" << endl;
// 	aurostd::execute("cat "+directory_RAW+"/" + _AFLOWIN_ + " | aflow --justbetween=\"[VASP_POSCAR_MODE_EXPLICIT]START\",\"[VASP_POSCAR_MODE_EXPLICIT]STOP\" > "+directory_RAW+"/POSCAR.orig");
// 	//  ExtractToStringEXPLICIT(Library_ICSD,Library_ICSD0,"[README_LIBRARY_ICSD1.TXT]START","[README_LIBRARY_ICSD1.TXT]STOP");
//       }
//     }

//     if(flag_DATA_RELAX_ || flag_EDATA_RELAX_) {  // CONTCAR.relax.EXT
//       for(uint i=1;i<=relax_max;i++) {
// 	file_LIB=directory_LIB+"/CONTCAR.relax"+aurostd::utype2string(i)+EXT;file_RAW=directory_RAW+"/CONTCAR.relax"+EXT;
// 	if(aurostd::FileExist(file_LIB)) { aurostd::CopyFile(file_LIB,file_RAW);aurostd::execute(DEFAULT_KZIP_BIN+" -9fd "+file_RAW);vfile.push_back("CONTCAR.relax"); }
//       }
//       if(aurostd::FileExist(file_RAW)) { cout << MESSAGE << ": ERROR - aflowlib::LIB2RAW_Loop_DATAs:[1] - file not prepared " << file_LIB << endl;exit(0); }
//     }

//     if(!flag_DATA_BANDS_) {
//       file_LIB=directory_LIB+"/POSCAR.bands"+EXT;
//       if(aurostd::FileExist(file_LIB)) { aurostd::CopyFile(file_LIB,file_RAW);aurostd::execute(DEFAULT_KZIP_BIN+" -9fd "+file_RAW);vfile.push_back("POSCAR.bands"); }
//     }

//     xstructure str,str_sp,str_sc;
//     vector<xstructure> vcif;

//     // PERFORM EDATA STEP
//     if(flag_EDATA_ORIG_ || flag_EDATA_RELAX_) cout << MESSAGE << " EDATA start: " << directory_RAW << endl;
//     if(flag_EDATA_ORIG_) { // ORIG
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
// 	cout << MESSAGE << " EDATA doing orig (POSCAR.orig): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"EDATA"); // 1=EDATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT);
// 	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
//       }
//     }
//     if(flag_EDATA_RELAX_) { // RELAX
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
// 	cout << MESSAGE << " EDATA doing relax (CONTCAR.relax): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"EDATA"); // EDATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT);
// 	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
//       }
//     }

//     //  if(flag_EDATA_BANDS_)
//     { // BANDS IF AVAILABLE DO IT
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
// 	cout << MESSAGE << " EDATA doing bands (POSCAR.bands): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"EDATA"); // EDATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT);
// 	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
//       }
//     }

//     // PERFORM DATA STEP
//     if(flag_DATA_ORIG_ || flag_DATA_RELAX_ || flag_DATA_BANDS_) cout << MESSAGE << " DATA start: " << directory_RAW << endl;
//     if(flag_DATA_ORIG_) { // ORIG
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
// 	cout << MESSAGE << " DATA doing orig (POSCAR.orig): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"DATA"); // DATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT);
//       }
//     }
//     if(flag_DATA_RELAX_) { // RELAX
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
// 	cout << MESSAGE << " DATA doing relax (CONTCAR.relax): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"DATA"); // DATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT);
//       }
//     }
//     if(flag_DATA_BANDS_) { // BANDS
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
// 	cout << MESSAGE << " DATA doing bands (POSCAR.bands): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"DATA"); // DATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT);
//       }
//     }

//     // NOW DO THE CIFS
//     aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CONTCAR.relax2",directory_RAW,"CONTCAR.relax2",vfile,MESSAGE);  // POSCAR.bands

//     //      cout << MESSAGE << " CIF creation: " << directory_LIB << endl;
//     if(vcif.size()==0) vcif.push_back(xstructure(directory_RAW+"/POSCAR.bands",IOVASP_AUTO));
//     xvector<double> nvec(3);nvec(1)=1;nvec(2)=1;nvec(3)=1;
//     double angle=45;

//     for (uint j=0;j<vcif.size();j++) {
//       vcif.at(j)=GetLTFVCell(nvec,angle,vcif.at(j));
//       if(j==0) cout << MESSAGE << " CIF creation: data.spacegroup_relax=" << data.spacegroup_relax << endl;
//       if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal" << endl;
//       if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim" << endl;
//       if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv" << endl;
//       stringstream oss;
//       //    for(uint i=0;i<vspecies.size();i++) cerr << vspecies.at(i) << endl;
//       vcif.at(j).species.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species.push_back(vspecies.at(i));
//       vcif.at(j).species_pp.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp.push_back(vspecies.at(i));
//       vcif.at(j).species_pp_type.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_type.push_back("");
//       vcif.at(j).species_pp_version.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_version.push_back("");
//       vcif.at(j).species_pp_ZVAL.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_ZVAL.push_back(0.0);
//       vcif.at(j).species_pp_vLDAU.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_vLDAU.push_back(deque<double>());
//       for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).name=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
//       for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).cleanname=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
//       pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
//       if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+".cif");
//       if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim.cif");
//       if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv.cif");
//       vcif.at(j).AddCorners();
//       oss.clear();oss.str("");
//       pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
//       if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal_corner" << endl;
//       if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim_corner" << endl;
//       if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv_corner" << endl;
//       if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_corner.cif");
//       if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim_corner.cif");
//       if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv_corner.cif");

//       // [OBSOLETE] aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+".cif");
//     }

//     // DONE
//     cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATA: end " << directory_LIB << endl;
//     return TRUE;
//   }
// }


// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Thermodynamics
// ***************************************************************************
namespace aflowlib {
  bool ExtractOUT_from_VASP_OUTCAR(string _file,const double& data_natoms,xOUTCAR& outcar) {
    bool flag=outcar.GetPropertiesFile(_file);
    if(aurostd::abs(data_natoms-(double) outcar.natoms)>0.1) { cerr << "ERROR ExtractOUT_from_VASP_OUTCAR: data_natoms(" << data_natoms << ")!= (int) outcar.natoms(" << outcar.natoms << ") ..." << endl;exit(0); }
    return flag;
  }
}

namespace aflowlib {
  bool AddFileNameBeforeExtension(string _file,string addendum,string& out_file) { // CO 171025
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string file=aurostd::CleanFileName(_file);
    out_file=file;
    if(file[file.size()-1]=='/' || file[file.size()-1]=='\\') { return false; } //not doing directories
    //if(aurostd::IsDirectory(file)) { return false; }
    //if(!aurostd::FileExist(file)) { return false; }
    //grab actually file
    string path="",rfile="";
    std::size_t pos=file.find_last_of("/\\");
    if(pos!=string::npos) { path=file.substr(0,pos); } //keep path empty if there's no directory info (just relative file name)
    rfile=file.substr(pos+1);
    if(LDEBUG) {
      cerr << "aflowlib::AddFileNameBeforeExtension:: path=" << path << endl;
      cerr << "aflowlib::AddFileNameBeforeExtension:: file=" << rfile << endl;
    }
    string possible_extensions;
    possible_extensions="out,txt,json,png,pdf,jpg,jpeg";  //normal file types first
    possible_extensions+=",tar,tbz,bz2,zip,gz";           //zip extensions last
    vector<string> vpossible_extensions;
    aurostd::string2tokens(possible_extensions,vpossible_extensions,",");
    //now split file by "."
    vector<string> parts;
    aurostd::string2tokens(rfile,parts,".");
    bool found=false;
    string rfile_new="";
    for(uint i=0;i<parts.size();i++) {
      for(uint j=0;j<vpossible_extensions.size()&&!found;j++) {
	if(parts[i]==vpossible_extensions[j]) {
	  found=true;
	  rfile_new+=(rfile_new.empty()?"":".")+addendum;
	}
      }
      rfile_new+=(rfile_new.empty()?"":".")+parts[i];
    }
    if(!found) { rfile_new+=(rfile_new.empty()?"":".")+addendum; }
    out_file=aurostd::CleanFileName((path.empty()?"":path)+(pos==string::npos?"":string(file.substr(pos,1)))+rfile_new);
    return true;
  }
}

namespace aflowlib {
  bool LIB2RAW_Loop_Thermodynamics(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE,bool LOCAL) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [1]" << endl;
    // ZIP-AGNOSTIC
    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
    deque<string> vcmd; aurostd::string2tokens("bzip2,xz,gzip",vcmd,",");
    if(vext.size()!=vcmd.size()) { cerr << "ERROR - aflowlib::LIB2RAW_Loop_Thermodynamics: vext.size()!=vcmd.size()" << endl;exit(0); }
    // CO and DX START 170713 - adding symmetry output to RAW
    _aflags aflags;
    aflags.Directory=directory_RAW;
    ofstream FileMESSAGE; //dummy ofstream, not really used
    stringstream message; //dummy stringstream, can output to cout, but not used right now
    string system_name=KBIN::ExtractSystemName(directory_LIB);
    // CO and DX STOP 170713 - adding symmetry output to RAW
    //    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics data.aurl=" << data.aurl << endl;
    // aurostd::StringSubst(directory_LIB,"/./","/");
    // aurostd::StringSubst(directory_RAW,"/./","/");
    vector<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Thermodynamics - begin " << directory_LIB << endl;
    if(LDEBUG) cerr << "directory_LIB=\"" << directory_LIB << "\"" << endl;
    if(LDEBUG) cerr << "directory_RAW=\"" << directory_RAW << "\"" << endl;
    if(directory_LIB.at(directory_LIB.size()-1)=='/')  directory_LIB=directory_LIB.substr(0,directory_LIB.size()-1);
    if(directory_RAW.at(directory_RAW.size()-1)=='/')  directory_RAW=directory_RAW.substr(0,directory_RAW.size()-1);
    if(LDEBUG) cerr << "directory_LIB=\"" << directory_LIB << "\"" << endl;
    if(LDEBUG) cerr << "directory_RAW=\"" << directory_RAW << "\"" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Thermodynamics - species = " << vspecies.size() << endl;
    uint relax_max=10;
    stringstream command;command.clear();command.str(std::string());
    stringstream aus_exec;
    // directories must exist already
    bool flag_EDATA_ORIG_=FALSE,flag_EDATA_RELAX_=FALSE,flag_EDATA_BANDS_=FALSE;
    bool flag_DATA_ORIG_=FALSE,flag_DATA_RELAX_=FALSE,flag_DATA_BANDS_=FALSE;
    bool flag_TIMING=FALSE;
    bool flag_ENERGY1=FALSE;
    bool flag_SG1=FALSE;
    bool flag_SG2=FALSE;
    bool flag_VOLDISTPARAMS=FALSE;
    bool flag_VOLDISTEVOLUTION=FALSE;
    bool flag_ICSD=FALSE,flag_MAGNETIC=FALSE,flag_LIB1=FALSE,flag_LIB2=FALSE;
    bool flag_ERROR=FALSE;

    
    string fileA_LIB,fileA_RAW,fileE_LIB,fileE_RAW,fileX_LIB,fileX_RAW,fileK_LIB,fileK_RAW,fileI_LIB,fileI_RAW,fileJ_LIB,fileJ_RAW,cmd;
    string FileName_OUTCAR_relax="";
    
    if(LOCAL) {
      flag_EDATA_ORIG_=flag_EDATA_RELAX_=flag_EDATA_BANDS_=TRUE;
      flag_DATA_ORIG_=flag_DATA_RELAX_=flag_DATA_BANDS_=TRUE;
      flag_TIMING=TRUE;
      flag_ENERGY1=FALSE;
      flag_SG1=TRUE;
      flag_SG2=TRUE;
      flag_VOLDISTPARAMS=TRUE;
      flag_VOLDISTEVOLUTION=TRUE;
      flag_ICSD=flag_LIB1=flag_LIB2=FALSE;
      flag_MAGNETIC=TRUE;
      flag_ERROR=FALSE;
    } else {
      aurostd::StringSubst(directory_LIB,"ELPASOLITES","AURO"); // PATCH
      aurostd::StringSubst(directory_LIB,"SCINT","ICSD"); // PATCH
      if(aurostd::substring2bool(directory_LIB,"ICSD"))     { flag_ICSD=TRUE;    flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB1"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB2"))     { flag_LIB2=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB3"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB4"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB5"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB6"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB7"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB8"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB9"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIBRARYX")) { flag_LIB2=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; } // [HISTORIC]
    }

    // check for flag_EDATA_BANDS_
    if(flag_EDATA_ORIG_ && flag_EDATA_RELAX_) {
      for(uint iext=0;iext<vext.size();iext++) {
	fileA_LIB=directory_LIB+"/POSCAR.bands"+vext.at(iext);
	fileA_RAW=directory_RAW+"/POSCAR.bands"+vext.at(iext);
	if(aurostd::FileExist(fileA_LIB)) {
	  flag_EDATA_BANDS_=TRUE;
	  aurostd::CopyFile(fileA_LIB,fileA_RAW);
	  aurostd::execute(vcmd.at(iext)+" -dqf "+fileA_RAW);
	  vfile.push_back("POSCAR.bands");
	}
      }
    }
    // check for flag_DATA_BANDS_
    if(flag_DATA_ORIG_ && flag_DATA_RELAX_) {
      for(uint iext=0;iext<vext.size();iext++) {
	fileA_LIB=directory_LIB+"/POSCAR.bands"+vext.at(iext);
	fileA_RAW=directory_RAW+"/POSCAR.bands"+vext.at(iext);
	if(aurostd::FileExist(fileA_LIB)) {
	  flag_DATA_BANDS_=TRUE;
	  aurostd::CopyFile(fileA_LIB,fileA_RAW);
	  aurostd::execute(vcmd.at(iext)+" -dqf "+fileA_RAW);
	  vfile.push_back("POSCAR.bands");
	}
      }
    }
    if(flag_ICSD) { ; } // dummy load
    if(flag_MAGNETIC) { ; } // dummy load
    if(flag_LIB1) { ; } // dummy load

    xstructure str_orig,str_relax1,str_relax;

    // DX - START
    // DX 20180526 [OBSOLETE] str_orig.directory = str_relax1.directory = str_relax.directory = aflags.Directory;
    // DX - END

    vector<string> tokens;

    double data1_energy_cell=0.0,data1_energy_atom=0.0;
    string data_sg1_pre="",data_sg1_mid="",data_sg1_post="",data_sg2_pre="",data_sg2_mid="",data_sg2_post="";
    string data_v_atom="",data_ucelld="";
    xvector<double> data_abcabc;
    vector<double> data_vcomposition;
    vector<xvector<double> > data_vforces;                      // QM FORCES calculation
    vector<xvector<double> > data_vpositions_cartesian;         // QM POSITIONS_CARTESIAN calculation

    xOUTCAR outcar;
    // xDOSCAR doscar;
    // xEIGENVAL eigenval;
    xKPOINTS kpoints;

    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [2]" << endl;
    // copy _AFLOWIN_ LOCK
    // _AFLOWIN_
    data.vloop.push_back("thermodynamics");
    // star

    aurostd::string2tokens(directory_LIB,tokens,"/");
    data.prototype=tokens.at(tokens.size()-1);
    aurostd::StringSubst(data.prototype,":LDAU2","");
    aurostd::StringSubst(data.prototype,"\n","");aurostd::StringSubst(data.prototype," ","");aurostd::StringSubst(data.prototype," ","");

    // FILES
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWIN_,directory_RAW,_AFLOWIN_,vfile,MESSAGE);  // _AFLOWIN_
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // LOCK

    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [3]" << endl;
    if(TRUE || flag_DATA_ORIG_ || flag_EDATA_ORIG_ || flag_SG1 || flag_SG2) {  // POSCAR.orig.EXT
      deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
      // if(flag_ORIG==FALSE) {
      bool found=FALSE;
      for(uint iext=0;iext<vext.size();iext++) { //      cerr << "BUILDING POSCAR.orig from POSCAR.orig.EXT" << endl;
	if(!found && aurostd::FileExist(directory_LIB+"/POSCAR.orig"+vext.at(iext))) {
	  found=TRUE;
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.orig",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
	}
      }
      for(uint iext=0;iext<vext.size();iext++) { //      cerr << "BUILDING POSCAR.orig from POSCAR.relax1.EXT" << endl;
	if(!found && aurostd::FileExist(directory_LIB+"/POSCAR.relax1"+vext.at(iext))) {
	  found=TRUE;
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
	}
      }
      if(!found) { found=TRUE;//      cerr << "BUILDING POSCAR.orig from " << _AFLOWIN_ << "" << endl;
	aurostd::execute(string("cat ")+directory_RAW+"/"+_AFLOWIN_+" | aflow --justbetween \"[VASP_POSCAR_MODE_EXPLICIT]START\" \"[VASP_POSCAR_MODE_EXPLICIT]STOP\" > "+directory_RAW+"/POSCAR.orig");
	//    ExtractToStringEXPLICIT(Library_ICSD,Library_ICSD0,"[README_LIBRARY_ICSD1.TXT]START","[README_LIBRARY_ICSD1.TXT]STOP");
      }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [4]" << endl;
    if(TRUE || flag_DATA_RELAX_ || flag_EDATA_RELAX_ || TRUE || flag_SG1 || flag_SG2) {  // CONTCAR.relax.EXT
      for(uint iext=0;iext<vext.size();iext++) {
	for(uint i=1;i<=relax_max;i++) {
	  fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax"+aurostd::utype2string(i)+vext.at(iext));
	  fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+vext.at(iext));
	  fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax"+aurostd::utype2string(i)+vext.at(iext));
	  fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+vext.at(iext));
	  fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.relax"+aurostd::utype2string(i)+vext.at(iext));
	  fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+vext.at(iext));
	  fileI_LIB=aurostd::CleanFileName(directory_LIB+"/INCAR.relax"+aurostd::utype2string(i)+vext.at(iext));
	  fileI_RAW=aurostd::CleanFileName(directory_RAW+"/INCAR.relax"+vext.at(iext));
	  if(aurostd::FileExist(fileX_LIB)) {
	    aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileX_RAW);vfile.push_back("CONTCAR.relax");
	  }
	  if(aurostd::FileExist(fileE_LIB)) {
	    aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileE_RAW);vfile.push_back("OUTCAR.relax");
	    FileName_OUTCAR_relax=fileE_LIB;
	  }
	  // if(aurostd::FileExist(fileK_LIB)) {
	  //   aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileK_RAW);vfile.push_back("KPOINTS.relax");
	  // }
	  // if(aurostd::FileExist(fileI_LIB)) {
	  //   aurostd::CopyFile(fileI_LIB,fileI_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileI_RAW);vfile.push_back("INCAR.relax");
	  // }
	}
      }
      for(uint iext=0;iext<vext.size();iext++) {
	if(!aurostd::FileExist(directory_RAW+"/CONTCAR.relax") && aurostd::FileExist(directory_LIB+"/CONTCAR.static"+vext.at(iext)) &&
	   !aurostd::FileExist(directory_RAW+"/OUTCAR.relax") && aurostd::FileExist(directory_LIB+"/OUTCAR.static"+vext.at(iext))) {
	  fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.static"+vext.at(iext));
	  fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+vext.at(iext));
	  fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.static"+vext.at(iext));
	  fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+vext.at(iext));
	  //  if(AFLOWLIB_VERBOSE)
	  cout << MESSAGE << " WARNING - PATCHING CONTCAR.relax with CONTCAR.static " << fileX_LIB << endl;
	  //	  if(AFLOWLIB_VERBOSE)
	  cout << MESSAGE << " WARNING - PATCHING OUTCAR.relax with OUTCAR.static " << fileE_LIB << endl;
	  if(aurostd::FileExist(fileX_LIB)) {
	    aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileX_RAW);vfile.push_back("CONTCAR.relax");
	  }
	  if(aurostd::FileExist(fileE_LIB)) {
	    aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileE_RAW);vfile.push_back("OUTCAR.relax");
	    FileName_OUTCAR_relax=fileE_LIB;
	  }
	}
      }
      if(aurostd::FileExist(fileX_RAW)) { cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_Loop_Thermodynamics:[1] - file not prepared " << fileX_LIB << endl;exit(0); }
      if(aurostd::FileExist(fileE_RAW)) { cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_Loop_Thermodynamics:[2] - file not prepared " << fileE_LIB << endl;exit(0); }
      if(aurostd::FileExist(fileK_RAW)) { cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_Loop_Thermodynamics:[3] - file not prepared " << fileK_LIB << endl;exit(0); }
      if(aurostd::FileExist(fileI_RAW)) { cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_Loop_Thermodynamics:[3] - file not prepared " << fileI_LIB << endl;exit(0); }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [5]" << endl;
    if(flag_SG1 || flag_SG2) {  // CONTCAR.relax1
      for(uint iext=0;iext<vext.size();iext++) {
	fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax1"+vext.at(iext));
	fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax1"+vext.at(iext));
	fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax1"+vext.at(iext));
	fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax1"+vext.at(iext));
	if(aurostd::FileExist(fileX_LIB)) {
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CONTCAR.relax1",directory_RAW,"CONTCAR.relax1",vfile,MESSAGE);
	}
	if(aurostd::FileExist(fileE_LIB)) {
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.relax1",directory_RAW,"OUTCAR.relax1",vfile,MESSAGE);
	}
	if(!aurostd::FileExist(directory_RAW+"/CONTCAR.relax1") && aurostd::FileExist(directory_LIB+"/CONTCAR.static"+vext.at(iext)) &&
	   !aurostd::FileExist(directory_RAW+"/OUTCAR.relax1") && aurostd::FileExist(directory_LIB+"/OUTCAR.static"+vext.at(iext))) {
	  fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.static"+vext.at(iext));
	  fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax1"+vext.at(iext));
	  fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.static"+vext.at(iext));
	  fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax1"+vext.at(iext));
	  //  if(AFLOWLIB_VERBOSE)
	  cout << MESSAGE << " WARNING - PATCHING CONTCAR.relax1 with CONTCAR.static " << fileX_LIB << endl;
	  //	  if(AFLOWLIB_VERBOSE)
	  cout << MESSAGE << " WARNING - PATCHING OUTCAR.relax1 with OUTCAR.static " << fileE_LIB << endl;
	  if(aurostd::FileExist(fileX_LIB)) {
	    aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileX_RAW);vfile.push_back("CONTCAR.relax1");
	  }
	  if(aurostd::FileExist(fileE_LIB)) {
	    aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileE_RAW);vfile.push_back("OUTCAR.relax1");
	  }
	}
	if(aurostd::FileExist(fileX_RAW)) {
	  cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_Loop_Thermodynamics:[4] - file not prepared " << fileX_LIB << endl;exit(0); }
	if(aurostd::FileExist(fileE_RAW)) {
	  cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_Loop_Thermodynamics:[5] - file not prepared " << fileE_LIB << endl;exit(0); }
      }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [6]" << endl;
    if(flag_DATA_RELAX_ || flag_EDATA_RELAX_ || TRUE) {  // OSZICAR.relax.EXT
      for(uint iext=0;iext<vext.size();iext++) {
	// trying to get an OUTCAR
	if(aurostd::FileExist(directory_LIB+"/OUTCAR.relax1"+vext.at(iext)))
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.relax1",directory_RAW,"OUTCAR.relax1",vfile,MESSAGE);  // _AFLOWIN_
	if(aurostd::FileExist(directory_LIB+"/CONTCAR.relax1"+vext.at(iext)))
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CONTCAR.relax1",directory_RAW,"CONTCAR.relax1",vfile,MESSAGE);  // _AFLOWIN_
	// if(aurostd::FileExist(directory_LIB+"/KPOINTS.relax1"+vext.at(iext)))
	// aflowlib::LIB2RAW_FileNeeded(directory_LIB,"KPOINTS.relax1",directory_RAW,"KPOINTS.relax1",vfile,MESSAGE);  // _AFLOWIN_
	fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.static"+vext.at(iext));
	fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+vext.at(iext));
	fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.static"+vext.at(iext));
	fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+vext.at(iext));
	fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.static"+vext.at(iext));
	fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+vext.at(iext));
	/* STEFANO STATIC=>RELAX
	  if(aurostd::FileExist(fileE_LIB) && aurostd::FileExist(fileX_LIB) && aurostd::FileExist(fileK_LIB)) {
	  // a bug in old aflow the fileK_LIB was directory_LIB+"/KPOINTS.relax"+aurostd::utype2string(i)+vext.at(iext) with i=3 or more, so it was never picked
	  // and this part was never done.
	  aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileE_RAW);vfile.push_back("OUTCAR.relax");
	  FileName_OUTCAR_relax=fileE_LIB;
	  // aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileX_RAW);vfile.push_back("CONTCAR.relax");
	  aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileK_RAW);vfile.push_back("KPOINTS.relax");
	  } else
	*/ {
	  fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax2"+vext.at(iext));
	  fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+vext.at(iext));
	  fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax2"+vext.at(iext));
	  fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+vext.at(iext));
	  fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.relax2"+vext.at(iext));
	  fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+vext.at(iext));
	  fileI_LIB=aurostd::CleanFileName(directory_LIB+"/INCAR.relax2"+vext.at(iext));
	  fileI_RAW=aurostd::CleanFileName(directory_RAW+"/INCAR.relax"+vext.at(iext));
	  if(aurostd::FileExist(fileE_LIB) && aurostd::FileExist(fileX_LIB) && aurostd::FileExist(fileK_LIB) && aurostd::FileExist(fileI_LIB)) {
	    aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileE_RAW);vfile.push_back("OUTCAR.relax");
	    FileName_OUTCAR_relax=fileE_LIB;
	    aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileX_RAW);vfile.push_back("CONTCAR.relax");
	    aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileK_RAW);vfile.push_back("KPOINTS.relax");
	    // nocopy aurostd::CopyFile(fileI_LIB,fileI_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileI_RAW);vfile.push_back("INCAR.relax");
	  } else {
	    fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax1"+vext.at(iext));
	    fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+vext.at(iext));
	    fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax1"+vext.at(iext));
	    fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+vext.at(iext));
	    fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.relax1"+vext.at(iext));
	    fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+vext.at(iext));
	    fileI_LIB=aurostd::CleanFileName(directory_LIB+"/INCAR.relax1"+vext.at(iext));
	    fileI_RAW=aurostd::CleanFileName(directory_RAW+"/INCAR.relax"+vext.at(iext));
	    if(aurostd::FileExist(fileE_LIB) && aurostd::FileExist(fileX_LIB) && aurostd::FileExist(fileK_LIB) && aurostd::FileExist(fileI_LIB)) {
	      aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileE_RAW);vfile.push_back("OUTCAR.relax");
	      FileName_OUTCAR_relax=fileE_LIB;
	      aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileX_RAW);vfile.push_back("CONTCAR.relax");
	      aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileK_RAW);vfile.push_back("KPOINTS.relax");
	      // nocopy aurostd::CopyFile(fileI_LIB,fileI_RAW);aurostd::execute(vcmd.at(iext)+" -dqf "+fileI_RAW);vfile.push_back("INCAR.relax");
	    }
	  }
	}
      }
    }

    // get code
    data.code="nan";
    aurostd::string2tokens(aurostd::execute2string("cat "+directory_RAW+"/OUTCAR.relax"+" | grep vasp | head -1"),tokens," ");
    if(tokens.size()>1) data.code=tokens.at(0);
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CODE = " << data.code << endl;

    // create structures
    // CO 171021 - and write xstr_json
    stringstream xstr_js;
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [7]" << endl;
    if(str_orig.num_each_type.size()==0 && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
      xstructure _str_orig(directory_RAW+"/POSCAR.orig",IOVASP_AUTO);
      str_orig=_str_orig;
      str_orig.ReScale(1.0);
      /*xstr_js.str("");xstr_js << xstructure2json(str_orig);aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_orig.json");*/
    } // CO 171025
    if(str_relax.num_each_type.size()==0 && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
      xstructure _str_relax(directory_RAW+"/CONTCAR.relax",IOVASP_AUTO);
      str_relax=_str_relax;
      str_relax.ReScale(1.0);
      xstr_js.str("");
      xstr_js << xstructure2json(str_orig);
      aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_relax.json");
    } // CO 171025
    if(str_relax.num_each_type.size()==0 && aurostd::FileExist(directory_RAW+"/CONTCAR.static")) {
      xstructure _str_relax(directory_RAW+"/CONTCAR.static",IOVASP_AUTO);
      str_relax=_str_relax;
      str_relax.ReScale(1.0);
      xstr_js.str("");
      xstr_js << xstructure2json(str_orig);
      aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_relax.json");
    } // CO 171025
    if(aurostd::FileExist(directory_RAW+"/CONTCAR.relax1")) {
      xstructure _str_relax1(directory_RAW+"/CONTCAR.relax1",IOVASP_AUTO);
      str_relax1=_str_relax1;
      str_relax1.ReScale(1.0);
      xstr_js.str("");
      xstr_js << xstructure2json(str_orig);
      aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_relax1.json");
    } // CO 171025
    // do the extractions
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [8]" << endl;

    // LOAD STRUCTURES
    data.nspecies=str_relax.num_each_type.size();
    data.natoms=str_relax.atoms.size();
    data.volume_cell=str_relax.GetVolume();
    data.volume_atom=str_relax.GetVolume()/(double) data.natoms;
    data_abcabc=Getabc_angles(str_relax.lattice,DEGREES);
    data.vgeometry.clear(); for(uint i=1;i<=6;i++) data.vgeometry.push_back(data_abcabc(i));

    //corey, get fpos now, cpos comes from outcar later (either way works)
    vector<string> fpos_strings;
    vector<string> fpos_strings_combined;
    data.vpositions_fractional.clear();
    for(uint i=0;i<str_relax.atoms.size();i++) {
      data.vpositions_fractional.push_back(str_relax.atoms.at(i).fpos);
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " POSITIONS_FRACTIONAL = " << data.vpositions_fractional.at(i)[1] << "," << data.vpositions_fractional.at(i)[2] << "," << data.vpositions_fractional.at(i)[3] << "," << endl;
      //prepare for string variant
      for(uint j=1;j<(uint)str_relax.atoms.at(i).fpos.rows+1;j++) {
	fpos_strings.push_back(aurostd::utype2string(str_relax.atoms.at(i).fpos[j],8));
      }
      fpos_strings_combined.push_back(aurostd::joinWDelimiter(fpos_strings,","));
      fpos_strings.clear();
    }
    data.positions_fractional=aurostd::joinWDelimiter(fpos_strings_combined,";");

    data.geometry="";
    if(data.vgeometry.size()) {
      for(uint i=0;i<data.vgeometry.size();i++) {
	data.geometry+=aurostd::utype2string(data.vgeometry.at(i),7)+(i<data.vgeometry.size()-1?",":"");
      }
    }

    data_vcomposition.clear();
    for(uint i=0;i<str_relax.num_each_type.size();i++) data_vcomposition.push_back(double(str_relax.num_each_type.at(i))/double(str_relax.atoms.size()));
    data.vstoichiometry=data_vcomposition;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NSPECIES = " << data.nspecies << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NATOMS = " << data.natoms << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VOLUME (A^3) = " << data.volume_cell << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VOLUME_ATOM (A^3) = " << data.volume_atom << "   " << directory_LIB << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " GEOMETRY (A,A,A,deg,deg,deg) = " << data.geometry << endl;
    // deque<int> num_each_type;              // WARNING: we use starting from 0
    // std::deque<_atom> atoms;               // WARNING: we use starting from 0
    // deque<string> species,species_pp;      // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    // deque<double> species_volume;          // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    // deque<double> species_mass;            // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    // cerr << str_relax.atoms.size() << " " << str_relax.atoms.at(0) << endl;
    // cerr << str_relax.num_each_type.size() << " " << str_relax.num_each_type.at(0) << endl;
    // cerr << str_relax.species.size() << " " << str_relax.species.size() << endl;
    // cerr << str_relax.species_pp.size() << " " << str_relax.species_pp.size() << endl;
    // cerr << str_relax.species_pp_type.size() << " " << str_relax.species_pp_type.size() << endl;
    // cerr << str_relax.species_pp_version.size() << " " << str_relax.species_pp_version.size() << endl;
    // cerr << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_pp_ZVAL.size() << endl;
    // cerr << str_relax.species_pp_vLDAU.size() << " " << str_relax.species_pp_vLDAU.size() << endl;
    // cerr << str_relax.species_volume.size() << " " << str_relax.species_volume.size() << endl;
    // cerr << str_relax.species_mass.size() << " " << str_relax.species_mass.size() << endl;

    // LOAD ENERGY DATA1
    if(flag_ENERGY1) {
      outcar.GetPropertiesFile(directory_RAW+"/OUTCAR.relax1",data.natoms,TRUE);
      //   ExtractDataOSZICAR(directory_RAW+"/OSZICAR.relax1",data.natoms,data1_dE,data1_dEN,data1_mag,data1_mag_atom);
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY1 total E0 (eV) = " << (data1_energy_cell=outcar.energy_cell) << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY1 per atom E0/N (eV) = " << (data1_energy_atom=outcar.energy_atom) << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY1 total E0 (eV) = " << outcar.enthalpy_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY1 per atom E0/N (eV) = " << outcar.enthalpy_atom << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN1 mag (\\mu) = " << outcar.mag_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN1 per atom mag/N (\\mu) = " << outcar.mag_atom << endl;
    }
    // LOAD ENERGY DATA
    outcar.GetPropertiesFile(directory_RAW+"/OUTCAR.relax",data.natoms,TRUE);
    kpoints.GetPropertiesFile(directory_RAW+"/KPOINTS.relax",TRUE);

    data.pressure=outcar.pressure;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PRESSURE (kB) = " << data.pressure << endl;
    data.vstress_tensor.clear();
    data.vstress_tensor.push_back(outcar.stress(1,1));data.vstress_tensor.push_back(outcar.stress(1,2));data.vstress_tensor.push_back(outcar.stress(1,3));
    data.vstress_tensor.push_back(outcar.stress(2,1));data.vstress_tensor.push_back(outcar.stress(2,2));data.vstress_tensor.push_back(outcar.stress(2,3));
    data.vstress_tensor.push_back(outcar.stress(3,1));data.vstress_tensor.push_back(outcar.stress(3,2));data.vstress_tensor.push_back(outcar.stress(3,3));
    data.stress_tensor="";
    if(data.vstress_tensor.size())
      for(uint i=0;i<data.vstress_tensor.size();i++)
	data.stress_tensor+=aurostd::utype2string(data.vstress_tensor.at(i),7)+(i<data.vstress_tensor.size()-1?",":"");
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " STRESS_TENSOR (kB) = " << data.stress_tensor << endl;
    data.pressure_residual=outcar.pressure_residual;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PRESSURE_RESIDUAL (kB) = " << data.pressure_residual << endl;
    data.Pulay_stress=outcar.Pulay_stress;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PULAY_STRESS (kB) = " << data.Pulay_stress << endl;
    data.energy_cell=outcar.energy_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY total E0 (eV) = " << data.energy_cell << endl;
    data.energy_atom=outcar.energy_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY per atom E0/N (eV) = " << data.energy_atom << endl;
    data.enthalpy_cell=outcar.enthalpy_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY total E0 (eV) = " << data.enthalpy_cell << endl;
    data.enthalpy_atom=outcar.enthalpy_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY per atom E0/N (eV) = " << data.enthalpy_atom << "   " << directory_LIB << endl;
    data.eentropy_cell=outcar.eentropy_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " E-ENTROPY total E0 (eV) = " << data.eentropy_cell << endl;
    data.eentropy_atom=outcar.eentropy_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " E-ENTROPY per atom E0/N (eV) = " << data.eentropy_atom << endl;
    data.PV_cell=outcar.PV_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PV total E0 (eV) = " << data.PV_cell << endl;
    data.PV_atom=outcar.PV_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PV per atom E0/N (eV) = " << data.PV_atom << endl;
    data.spin_cell=outcar.mag_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN mag (\\mu) = " << data.spin_cell << endl;
    data.spin_atom=outcar.mag_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN per atom mag/N (\\mu) = " << data.spin_atom << "   " << directory_LIB << endl;
    // CO 180130 - START
    //moving FROM magnetic loop so we keep spin/cell, spin/atom, and spinD all together
    data.spinD="";
    data.vspinD.clear();
    if(outcar.vmag.size()) {
      for(uint i=0;i<(uint) outcar.vmag.size();i++) {
	data.spinD+=aurostd::utype2string<double>(outcar.vmag.at(i),5)+(i<outcar.vmag.size()-1?",":"");
	data.vspinD.push_back(outcar.vmag.at(i));
      }
    } else {
      for(uint i=0;i<outcar.natoms;i++) {  //use outcar.natoms as there can be a primitivization between relax and static
	data.spinD+=aurostd::utype2string<double>(0)+(i<outcar.natoms-1?",":"");
	data.vspinD.push_back(0.0);
      }
    }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIND (\\mu) = " << data.spinD << "   " << directory_LIB << endl;
    // CO 180130 - STOP
    data.energy_cutoff=outcar.ENCUT;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY_CUTOFF (eV) = " << data.energy_cutoff << endl;
    data.delta_electronic_energy_convergence=outcar.total_energy_change;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DELTA_ELECTRONIC_ENERGY_CONVERGENCE (eV) = " << data.delta_electronic_energy_convergence << endl; // CORMAC
    data.delta_electronic_energy_threshold=outcar.EDIFF;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DELTA_ELECTRONIC_ENERGY_THRESHOLD (eV) = " << data.delta_electronic_energy_threshold << endl; // CORMAC
    data.nkpoints=kpoints.nkpoints;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NKPOINTS (from KPOINTS) = " << data.nkpoints << endl;
    data.kpoints_nnn_relax=kpoints.nnn_kpoints;
    data.kpoints=aurostd::utype2string(kpoints.nnn_kpoints[1])+","+aurostd::utype2string(kpoints.nnn_kpoints[2])+","+aurostd::utype2string(kpoints.nnn_kpoints[3]);
    data.nkpoints_irreducible=outcar.nkpoints_irreducible;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NKPOINTS_IRREDUCIBLE (from OUTCAR) = " << data.nkpoints_irreducible << endl;
    data.kppra=data.natoms*kpoints.nkpoints;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NKPPRA = " << data.kppra << endl;
    data_vforces.clear(); for(uint i=0;i<outcar.vforces.size();i++) { data_vforces.push_back(outcar.vforces.at(i)); }
    data.vforces=data_vforces;
    data_vpositions_cartesian.clear(); for(uint i=0;i<outcar.vpositions_cartesian.size();i++)  data_vpositions_cartesian.push_back(outcar.vpositions_cartesian.at(i));
    data.vpositions_cartesian=data_vpositions_cartesian;
    uint precfp=8;
    data.forces="";
    data.positions_cartesian="";
    for(uint i=0;i<(uint) data.natoms;i++) {
      data.positions_cartesian+=aurostd::utype2string(data_vpositions_cartesian.at(i)[1],precfp)+","+aurostd::utype2string(data_vpositions_cartesian.at(i)[2],precfp)+","+aurostd::utype2string(data_vpositions_cartesian.at(i)[3],precfp);
      if(i<data.natoms-1) data.positions_cartesian+=";";
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " POSITIONS_CARTESIAN = " << data_vpositions_cartesian.at(i)[1] << "," << data_vpositions_cartesian.at(i)[2] << "," << data_vpositions_cartesian.at(i)[3] << "," << endl;
    }
    for(uint i=0;i<(uint) data.natoms;i++) {
      data.forces+=aurostd::utype2string(data_vforces.at(i)[1],precfp)+","+aurostd::utype2string(data_vforces.at(i)[2],precfp)+","+aurostd::utype2string(data_vforces.at(i)[3],precfp);
      if(i<data.natoms-1) data.forces+=";";
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " FORCES(eV/Angst) = " << data_vforces.at(i)[1] << "," << data_vforces.at(i)[2] << "," << data_vforces.at(i)[3] << "," << endl;
    }

    // LOAD SPECIES
    //  cerr << data.nspecies << " " << str_relax.species.size() << " " << str_relax.species_pp.size() << " " << str_relax.species_pp_type.size() << " " << str_relax.species_pp_version.size() << " " << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_volume.size() << " " << str_relax.species_mass.size() << endl;

    str_relax.species.clear();str_relax.species_pp.clear();str_relax.species_pp_type.clear();str_relax.species_pp_version.clear();str_relax.species_pp_ZVAL.clear();
    str_relax.species_pp_vLDAU.clear();str_relax.species_volume.clear();str_relax.species_mass.clear();

    // try OUTCARs
    data.vspecies.clear();
    data.vspecies_pp.clear();
    data.vspecies_pp_version.clear();
    data.vspecies_pp_ZVAL.clear();

    for(uint itry=1;itry<=5&&str_relax.species.size()==0;itry++) {  // ONLY OUTCAR POTCARS ARE OBSOLETE
      deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
      for(uint iext=0;iext<vext.size();iext++) {
	string stry="";
	if(itry==1) stry="OUTCAR.relax1"+vext.at(iext);
	if(itry==2) stry="OUTCAR.relax2"+vext.at(iext);
	if(itry==3) stry="OUTCAR.relax3"+vext.at(iext);
	if(itry==4) stry="OUTCAR.static"+vext.at(iext);
	if(itry==5) stry="OUTCAR.bands"+vext.at(iext);
	fileE_LIB=directory_LIB+"/"+stry;
	if(aurostd::FileExist(fileE_LIB)) {
	  // [OBSOLETE]	stringstream stream_outcar(aurostd::execute2string("bzcat"+" "+fileE_LIB));
	  // [OBSOLETE]	outcar.GetProperties(stream_outcar);
	  outcar.GetPropertiesFile(fileE_LIB);
	  // DEBUG	for(uint i=0;i<outcar.species.size();i++) cerr << "outcar.species.at(i)=" << outcar.species.at(i) << endl;
	  str_relax.species.clear(); for(uint i=0;i<outcar.species.size();i++) { str_relax.species.push_back(outcar.species.at(i));data.vspecies.push_back(outcar.species.at(i)); } // for aflowlib_libraries.cpp
	  str_relax.species_pp.clear(); for(uint i=0;i<outcar.species_pp.size();i++) { str_relax.species_pp.push_back(outcar.species_pp.at(i));data.vspecies_pp.push_back(outcar.species_pp.at(i)); } // for aflowlib_libraries.cpp
	  str_relax.species_pp_type.clear(); for(uint i=0;i<outcar.species_pp_type.size();i++) { str_relax.species_pp_type.push_back(outcar.species_pp_type.at(i));/*data.vspecies_pp_type.push_back(outcar.vspecies_pp_type.at(i));*/} // for aflowlib_libraries.cpp
	  str_relax.species_pp_version.clear(); for(uint i=0;i<outcar.species_pp_version.size();i++) { str_relax.species_pp_version.push_back(outcar.species_pp_version.at(i));data.vspecies_pp_version.push_back(outcar.species_pp_version.at(i)); } // for aflowlib_libraries.cpp
	  str_relax.species_pp_ZVAL.clear(); for(uint i=0;i<outcar.vZVAL.size();i++) { str_relax.species_pp_ZVAL.push_back(outcar.vZVAL.at(i));data.vspecies_pp_ZVAL.push_back(outcar.vZVAL.at(i)); } // for aflowlib_libraries.cpp
	  data.dft_type=outcar.pp_type;
	  data.vdft_type.push_back(outcar.pp_type);  // CO, this is technically a vector (RESTAPI paper)
	  str_relax.species_pp_vLDAU.clear(); for(uint i=0;i<outcar.species_pp_vLDAU.size();i++) str_relax.species_pp_vLDAU.push_back(outcar.species_pp_vLDAU.at(i));  // for aflowlib_libraries.cpp
	  data.ldau_TLUJ=outcar.string_LDAU;
	  if(AFLOWLIB_VERBOSE && data.ldau_TLUJ.size()) cout << MESSAGE << " LDAU_string=" << data.ldau_TLUJ << endl;
	}
      }
    }

    if(str_relax.species.size()==0) { cerr << MESSAGE << " ERROR - OUTCAR/POTCAR not FOUND in " << directory_LIB << endl;exit(0); }

    for(uint j=0;j<str_relax.species.size();j++) {
      str_relax.species_volume.push_back(GetAtomVolume(str_relax.species.at(j)));
      str_relax.species_mass.push_back(GetAtomMass(str_relax.species.at(j)));
    }

    // cerr << data.nspecies << " " << str_relax.species.size() << " " << str_relax.species_pp.size() << " " << str_relax.species_pp_type.size() << " " << str_relax.species_pp_version.size() << " " << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_pp_vLDAU.size() << " " << str_relax.species_volume.size() << " " << str_relax.species_mass.size() << endl;
    if(data.nspecies!=str_relax.species.size()) { cerr << MESSAGE << " ERROR - data.nspecies[" << data.nspecies << "]!=str_relax.species.size()[" << str_relax.species.size() << "]" << endl << str_relax << endl;;exit(0); }
    if(data.nspecies!=str_relax.species_pp.size()) { cerr << MESSAGE << " ERROR - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp.size()[" << str_relax.species_pp.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_pp_type.size()) { cerr << MESSAGE << " ERROR - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp_type.size()[" << str_relax.species_pp_type.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_pp_version.size()) { cerr << MESSAGE << " ERROR - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp_version.size()[" << str_relax.species_pp_version.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_pp_ZVAL.size()) { cerr << MESSAGE << " ERROR - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp_ZVAL.size()[" << str_relax.species_pp_ZVAL.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_volume.size()) { cerr << MESSAGE << " ERROR - data.nspecies[" << data.nspecies << "]!=str_relax.species_volume.size()[" << str_relax.species_volume.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_mass.size()) { cerr << MESSAGE << " ERROR - data.nspecies[" << data.nspecies << "]!=str_relax.species_mass.size()[" << str_relax.species_mass.size() << "]" << endl;exit(0); }

    data.compound="";
    data.composition="";
    data.vcomposition.clear();
    data.density=0.0;
    data.stoichiometry="";
    //data.stoich=""; // CO 171026 - we do this here now, mostly obsolete, keep for legacy
    data.species="";
    data.species_pp="";
    data.species_pp_version="";
    data.species_pp_ZVAL="";
    data.valence_cell_iupac=0.0;
    data.valence_cell_std=0.0;
    for(uint i=0;i<(uint) data.nspecies;i++) {
      data.compound+=str_relax.species.at(i)+aurostd::utype2string(str_relax.num_each_type.at(i)); // if(i<data.nspecies-1) data.compound+=",";
      data.composition+=aurostd::utype2string(str_relax.num_each_type.at(i)); if(i<data.nspecies-1) data.composition+=",";
      data.vcomposition.push_back(str_relax.num_each_type.at(i));
      data.density+=(double) str_relax.num_each_type.at(i)*GetAtomMass(str_relax.species.at(i));
      data.stoichiometry+=aurostd::utype2string(data_vcomposition.at(i),9); if(i<data.nspecies-1) data.stoichiometry+=",";
      //data.stoich+=aurostd::utype2string(data_vcomposition.at(i),4); if(i<data.nspecies-1) data.stoichiometry+="   "; //mimic old BS format, 4 digits of accuracy and 3 spaces between, stoich=0.5000   0.1667   0.3333
      data.species+=str_relax.species.at(i);if(i<data.nspecies-1) data.species+=",";
      data.species_pp+=str_relax.species_pp.at(i);if(i<data.nspecies-1) data.species_pp+=",";
      // [UNUSED]    data.species_pp_type+=str_relax.species_pp_type.at(i);if(i<data.nspecies-1) data.species_pp_type+=",";
      data.species_pp_version+=str_relax.species_pp_version.at(i);if(i<data.nspecies-1) data.species_pp_version+=",";
      data.species_pp_ZVAL+=aurostd::utype2string(str_relax.species_pp_ZVAL.at(i));if(i<data.nspecies-1) data.species_pp_ZVAL+=",";
      //  cerr << "SIMPLE=" << str_relax.species.at(i) << endl;
      data.valence_cell_iupac+=str_relax.num_each_type.at(i)*GetAtomValenceIupac(str_relax.species.at(i));
      data.valence_cell_std+=str_relax.num_each_type.at(i)*GetAtomValenceStd(str_relax.species.at(i));
    }

    // density
    data.density/=data.volume_cell;
    data.density*=1000.0; // grams instead of kilos
    data.density*=1e8*1e8*1e8; // cm^3 instead of A^3
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DENSITY (grams/cm^3) = " << data.density << endl;
    // scintillation_attenuation_length
    data.scintillation_attenuation_length=0.0;
    data.scintillation_attenuation_length=GetCompoundAttenuationLenght(str_relax.species,str_relax.num_each_type,(double) data.density);
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SCINTILLATION_ATTENUATION_LENGTH (cm) = " << data.scintillation_attenuation_length << endl;

    // [UNUSED] if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_type = " << data.species_pp_type << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL dft_type=" << data.dft_type << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_version = " << data.species_pp_version << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_ZVAL = " << data.species_pp_ZVAL << endl;
    // reference
    bool isLDAUcalc=FALSE;
    if(str_relax.species_pp_vLDAU.at(0).size()>0) isLDAUcalc=TRUE;

    bool FORMATION_CALC=FALSE;
    if(isLDAUcalc==FALSE) FORMATION_CALC=TRUE;

    if(FORMATION_CALC==TRUE) {
      for(uint i=0;i<(uint) data.nspecies;i++) {
	//	if(EnthalpyReferenceAvailable(str_relax.species_pp.at(i),str_relax.species_pp_type.at(i),str_relax.species_pp_vLDAU.at(i))==FALSE) {  // OLD
	if(EnthalpyReferenceAvailable(str_relax.species_pp_version.at(i),str_relax.species_pp_vLDAU.at(i))==FALSE) {
	  FORMATION_CALC=FALSE; // SCAN AVAILABILITY
	  if(AFLOWLIB_VERBOSE) cout << MESSAGE << " REFERENCE NOT AVAILABLE SPECIE_pp=" << str_relax.species_pp_version.at(i) << "   Href=nothing" << endl;
	}
      }
    }

    if(FORMATION_CALC==TRUE) { // no LDAU yet
      if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [FCALC=1]" << endl;
      vector<double> venthalpy_atom_ref;
      for(uint i=0;i<(uint) data.nspecies;i++) {
	//      string pseudopotential,string type,vector<double> LDAU
	double enthalpy_atom_ref=data.enthalpy_atom; // if there is 1 then there is only one
	//	if(data.nspecies>1)
	//	enthalpy_atom_ref=EnthalpyReference(str_relax.species_pp.at(i),str_relax.species_pp_type.at(i),str_relax.species_pp_vLDAU.at(i),cerr);  // OLD
	enthalpy_atom_ref=EnthalpyReference(str_relax.species_pp_version.at(i),str_relax.species_pp_vLDAU.at(i),cerr);
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " REFERENCE SPECIE_pp_version=" << str_relax.species_pp_version.at(i) << "  Href=" << enthalpy_atom_ref << endl;
	venthalpy_atom_ref.push_back(enthalpy_atom_ref);
      }

      if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [FCALC=2]" << endl;
      if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [FCALC=2] data.nspecies=" << data.nspecies << endl;
      if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [FCALC=2] venthalpy_atom_ref.size()=" << venthalpy_atom_ref.size() << endl;

      // calculation of REF
      //     cerr << data.enthalpy << endl;
      data.enthalpy_formation_cell=data.enthalpy_cell;
      data.enthalpy_formation_atom=data.enthalpy_atom;

      for(uint i=0;i<(uint) data.nspecies;i++) data.enthalpy_formation_cell=data.enthalpy_formation_cell-(double(venthalpy_atom_ref.at(i)*str_relax.num_each_type.at(i)));
      //   for(uint i=0;i<(uint) data.nspecies;i++) data.enthalpy_formation_atom=data.enthalpy_formation_atom-venthalpy_atom_ref.at(i)*double(str_relax.num_each_type.at(i))/double(str_relax.atoms.size());
      data.enthalpy_formation_atom=data.enthalpy_formation_cell/double(str_relax.atoms.size());

      if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [FCALC=3]" << endl;

      data.entropic_temperature=0;
      if(data_vcomposition.size()>1) {
	for(uint i=0;i<(uint) data_vcomposition.size();i++)
	  if(data_vcomposition.at(i)>_EPSILON_COMPOSITION_ && data_vcomposition.at(i)<1-_EPSILON_COMPOSITION_)
	    data.entropic_temperature+=data_vcomposition.at(i)*logl(data_vcomposition.at(i));
	data.entropic_temperature=data.enthalpy_formation_atom/(data.entropic_temperature*KBOLTZEV);
      }
      // cerr << data.enthalpy_formation_cell << endl << data.enthalpy_formation_cell/str_relax.atoms.size() << endl << data.enthalpy_formation_atom << endl;
      if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [FCALC=4]" << endl;
    }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VALENCE_IUPAC = " << data.valence_cell_iupac << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VALENCE_STD = " << data.valence_cell_std << endl;

    //   aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energyd=" << data_dE;
    //   aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energyd_atom=" << data_dEN;

    if(FORMATION_CALC==TRUE) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION total E0 (eV) = " << data.enthalpy_formation_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION per atom E0/N (eV) = " << data.enthalpy_formation_atom << "   " << directory_LIB << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTROPIC_TEMPERATURE (eV) = " << data.entropic_temperature*KBOLTZEV << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTROPIC_TEMPERATURE (K) = " << data.entropic_temperature << "   " << directory_LIB << endl;
    }
    // [FIX]  if(flag_ENERGY1) aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy1_cell=" << data1_energy_cell;
    // [FIX] if(flag_ENERGY1) aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy1_atom=" << data1_energy_atom;
    // DONE WITH THERMO

    // do the TIMING
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [9]" << endl;
    data.calculation_cores=1;
    data.calculation_time=0.0;
    data.calculation_memory=0.0;
    if(flag_TIMING) {  // OUTCAR.relax.EXT
      for(uint i=1;i<=relax_max+2;i++) {
	deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
	for(uint iext=0;iext<vext.size();iext++) {
	  fileE_LIB=directory_LIB+"/OUTCAR.relax"+aurostd::utype2string(i)+vext.at(iext);
	  if(i==relax_max+1) fileE_LIB=directory_LIB+"/OUTCAR.static"+vext.at(iext);  // do the static
	  if(i==relax_max+2) fileE_LIB=directory_LIB+"/OUTCAR.bands"+vext.at(iext);  // do the bands
	  // cerr << fileE_LIB << endl;
	  if(aurostd::FileExist(fileE_LIB)) {
	    xOUTCAR outcar_tmp; // cerr << fileE_LIB << endl;
	    outcar_tmp.GetPropertiesFile(fileE_LIB); // OK
	    data.calculation_cores=aurostd::max((int) data.calculation_cores,(int) outcar_tmp.calculation_cores);
	    //	data.calculation_time+=double(data.calculation_cores)*outcar_tmp.calculation_time;
	    data.calculation_time+=outcar_tmp.calculation_time;  // will multiply after
	    data.calculation_memory=aurostd::max(data.calculation_memory,outcar_tmp.calculation_memory);
	    xOUTCAR outcar;
	    if(i==relax_max+2) { outcar.GetPropertiesFile(directory_LIB+"/OUTCAR.bands"+vext.at(iext)); } // cerr << "xOUTCAR.Efermi=" << outcar.Efermi << endl; }
	    xEIGENVAL eigenval;
	    if(i==relax_max+2) { eigenval.GetPropertiesFile(directory_LIB+"/EIGENVAL.bands"+vext.at(iext)); }
	    xDOSCAR doscar;
	    if(i==relax_max+2) { doscar.GetPropertiesFile(directory_LIB+"/DOSCAR.bands"+vext.at(iext)); } // cerr << "xDOSCAR.Efermi=" << doscar.Efermi << " " << outcar.Efermi-doscar.Efermi<< endl; }
	  }
	}
      }

      if(_APENNSY_STYLE_OLD_) aurostd::string2file(aurostd::utype2string((long(100*data.calculation_time))/100)+"\n",directory_RAW+"/"+DEFAULT_FILE_TIME_OUT);
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CALCULATION time (sec) = " << data.calculation_time << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CALCULATION mem (MB) = " << data.calculation_memory << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CALCULATION cores = " << data.calculation_cores << endl;
      //	exit(0);
    }
    // do the SG
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [10]" << endl;
    if(flag_SG1 || flag_SG2) {  // POSCAR.orig CONTCAR.relax1 CONTCAR.relax
#ifdef USE_PLATON_SG
      string space_group_calculator_sg1_cmd=" | aflow --platonSG=3.0,0.5,0.5,0.5";
      string space_group_calculator_sg2_cmd=" | aflow --platonSG=1.5,0.25,0.25,0.25";
#endif
#ifdef USE_AFLOW_SG
      // DX and CO - START
      // DX [OBSOLETE] string space_group_calculator_sg1_cmd=" | aflow --aflowSG=0.001";
      // DX [OBSOLETE] string space_group_calculator_sg2_cmd=" | aflow --aflowSG=0.00075";
      string space_group_calculator_sg1_cmd=" | aflow --aflowSG=loose";
      string space_group_calculator_sg2_cmd=" | aflow --aflowSG=tight";
      // DX and CO - END
#endif
      string ilattice_cmd=" | aflow --ilattice 2.0";
      if(flag_SG1) {
	stringstream ssfile;
	cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg1_cmd << endl;
	ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg1_cmd << endl;
	ssfile << directory_RAW << endl;
	// I run outside so a segfault+core does not block the aflow
	// INSIDE no more inside to avoid BOMBING
	// str_orig.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_pre=str_orig.spacegroup;     //  aflow --platonSG 3.0 0.5 0.5 0.5
	// str_relax1.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_mid=str_relax1.spacegroup; //  aflow --platonSG 3.0 0.5 0.5 0.5
	// str_relax.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_post=str_relax.spacegroup;  //  aflow --platonSG 3.0 0.5 0.5 0.5
	// OUTSIDE
	// sg1_pre
	data_sg1_pre=aurostd::execute2string("cat "+directory_RAW+"/POSCAR.orig"+space_group_calculator_sg1_cmd);
	aurostd::string2tokens(data_sg1_pre,tokens,"#");if(tokens.size()!=2) { data_sg1_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_pre=NOSG; } else { aurostd::StringSubst(data_sg1_pre,"\n",""); }}
	if(data_sg1_pre==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_pre" << endl;data_sg1_pre=aurostd::execute2string("cat "+directory_RAW+"/POSCAR.orig"+ilattice_cmd+space_group_calculator_sg1_cmd); }
	aurostd::string2tokens(data_sg1_pre,tokens,"#");if(tokens.size()!=2) { data_sg1_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_pre=NOSG; } else { aurostd::StringSubst(data_sg1_pre,"\n",""); }}
	if(aurostd::substring2bool(data_sg1_pre,"SymbolnotKnown")) data_sg1_pre=NOSG;  // give up
	// sg1_mid
	data_sg1_mid=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax1"+space_group_calculator_sg1_cmd);
	aurostd::string2tokens(data_sg1_mid,tokens,"#");if(tokens.size()!=2) { data_sg1_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_mid=NOSG; } else { aurostd::StringSubst(data_sg1_mid,"\n",""); }}
	if(data_sg1_mid==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_mid" << endl;data_sg1_mid=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax1"+ilattice_cmd+space_group_calculator_sg1_cmd); }
	aurostd::string2tokens(data_sg1_mid,tokens,"#");if(tokens.size()!=2) { data_sg1_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_mid=NOSG; } else { aurostd::StringSubst(data_sg1_mid,"\n",""); }}
	if(aurostd::substring2bool(data_sg1_mid,"SymbolnotKnown")) data_sg1_mid=NOSG;  // give up
	// sg1_mid
	data_sg1_post=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax"+space_group_calculator_sg1_cmd);
	aurostd::string2tokens(data_sg1_post,tokens,"#");if(tokens.size()!=2) { data_sg1_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_post=NOSG; } else { aurostd::StringSubst(data_sg1_post,"\n",""); }}
	if(data_sg1_post==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_post" << endl;data_sg1_post=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax"+ilattice_cmd+space_group_calculator_sg1_cmd); }
	aurostd::string2tokens(data_sg1_post,tokens,"#");if(tokens.size()!=2) { data_sg1_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_post=NOSG; } else { aurostd::StringSubst(data_sg1_post,"\n",""); }}
	if(aurostd::substring2bool(data_sg1_post,"SymbolnotKnown")) data_sg1_post=NOSG;  // give up
	// DONE
	ssfile << "PRE  " << data_sg1_pre << endl;
	ssfile << "MID  " << data_sg1_mid << endl;
	ssfile << "POST " << data_sg1_post << endl;
	//      aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP1_OUT);
	data.sg=data_sg1_pre+","+data_sg1_mid+","+data_sg1_post;
	data.vsg.clear();data.vsg.push_back(data_sg1_pre);data.vsg.push_back(data_sg1_mid);data.vsg.push_back(data_sg1_post); // CO 171202
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP1 = " << data.sg << endl;
      }
      if(flag_SG2) {
	stringstream ssfile;
	cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg2_cmd << endl;
	ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg2_cmd << endl;
	ssfile << directory_RAW << endl;
	// I run outside so a segfault+core does not block the aflow
	// INSIDE
	// str_orig.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_pre=str_orig.spacegroup;     //  aflow --platonSG 1.5 0.25 0.25 0.25
	// str_relax1.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_mid=str_relax1.spacegroup; //  aflow --platonSG 1.5 0.25 0.25 0.25
	// str_relax.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_post=str_relax.spacegroup;  //  aflow --platonSG 1.5 0.25 0.25 0.25
	// OUTSIDE
	// sg2_pre
	data_sg2_pre=aurostd::execute2string("cat "+directory_RAW+"/POSCAR.orig"+space_group_calculator_sg2_cmd);
	aurostd::string2tokens(data_sg2_pre,tokens,"#");if(tokens.size()!=2) { data_sg2_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_pre=NOSG; } else { aurostd::StringSubst(data_sg2_pre,"\n",""); }}
	if(data_sg2_pre==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_pre" << endl;data_sg2_pre=aurostd::execute2string("cat "+directory_RAW+"/POSCAR.orig"+ilattice_cmd+space_group_calculator_sg2_cmd); }
	aurostd::string2tokens(data_sg2_pre,tokens,"#");if(tokens.size()!=2) { data_sg2_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_pre=NOSG; } else { aurostd::StringSubst(data_sg2_pre,"\n",""); }}
	if(aurostd::substring2bool(data_sg2_pre,"SymbolnotKnown")) data_sg2_pre=NOSG;  // give up
	// sg2_mid
	data_sg2_mid=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax1"+space_group_calculator_sg2_cmd);
	aurostd::string2tokens(data_sg2_mid,tokens,"#");if(tokens.size()!=2) { data_sg2_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_mid=NOSG; } else { aurostd::StringSubst(data_sg2_mid,"\n",""); }}
	if(data_sg2_mid==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_mid" << endl;data_sg2_mid=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax1"+ilattice_cmd+space_group_calculator_sg2_cmd); }
	aurostd::string2tokens(data_sg2_mid,tokens,"#");if(tokens.size()!=2) { data_sg2_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_mid=NOSG; } else { aurostd::StringSubst(data_sg2_mid,"\n",""); }}
	if(aurostd::substring2bool(data_sg2_mid,"SymbolnotKnown")) data_sg2_mid=NOSG;  // give up
	// sg2_mid
	data_sg2_post=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax"+space_group_calculator_sg2_cmd);
	aurostd::string2tokens(data_sg2_post,tokens,"#");if(tokens.size()!=2) { data_sg2_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_post=NOSG; } else { aurostd::StringSubst(data_sg2_post,"\n",""); }}
	if(data_sg2_post==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_post" << endl;data_sg2_post=aurostd::execute2string("cat "+directory_RAW+"/CONTCAR.relax"+ilattice_cmd+space_group_calculator_sg2_cmd); }
	aurostd::string2tokens(data_sg2_post,tokens,"#");if(tokens.size()!=2) { data_sg2_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_post=NOSG; } else { aurostd::StringSubst(data_sg2_post,"\n",""); }}
	if(aurostd::substring2bool(data_sg2_post,"SymbolnotKnown")) data_sg2_post=NOSG;  // give up
	// DONE
	ssfile << "PRE  " << data_sg2_pre << endl;
	ssfile << "MID  " << data_sg2_mid << endl;
	ssfile << "POST " << data_sg2_post << endl;
	// aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP2_OUT);
	data.sg2=data_sg2_pre+","+data_sg2_mid+","+data_sg2_post;
	data.vsg2.clear();data.vsg2.push_back(data_sg2_pre);data.vsg2.push_back(data_sg2_mid);data.vsg2.push_back(data_sg2_post); // CO 171202
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP2 = " << data.sg2 << endl;

	aurostd::string2tokens(data_sg2_pre,tokens,"#");tokens.at(tokens.size()-1);
	data.spacegroup_orig="0"; if(tokens.size()>1) { data.spacegroup_orig=tokens.at(1); }
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP_ORIG = " << data.spacegroup_orig << endl;

	aurostd::string2tokens(data_sg2_post,tokens,"#");tokens.at(tokens.size()-1);
	data.spacegroup_relax="0"; if(tokens.size()>1) { data.spacegroup_relax=tokens.at(1); }
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP_RELAX = " << data.spacegroup_relax << endl;

      }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [11]" << endl;
    // VOLDISTParams
    if(flag_VOLDISTPARAMS) {  // CONTCAR.relax
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [VOLDISTParams]" << endl;
      stringstream ssfile;
      ssfile << "Volume Distance Bonds Analyzer: RELAX" << endl;
      ssfile << directory_RAW << endl;
      data.vnbondxx=GetNBONDXX(str_relax);  // CO 171024
      string data1=pflow::PrintData1(str_relax,-1.0);
      vector<string> vdata1;aurostd::string2vectorstring(data1,vdata1);
      for(uint i=0;i<vdata1.size();i++) {
	if(aurostd::substring2bool(vdata1.at(i),"Stoich Str1")) {
	  data.stoich=aurostd::substring2string(vdata1.at(i),"Stoich Str1: ",FALSE);
	  ssfile << "STOICH " << data.stoich << endl;
	}
	if(aurostd::substring2bool(vdata1.at(i),"Vol Str1")) {
	  data_v_atom=aurostd::substring2string(vdata1.at(i),"Vol Str1: ",FALSE);
	  ssfile << "V_ATOM " << data_v_atom << endl;
	}
	if(aurostd::substring2bool(vdata1.at(i),"Cell Str1")) {
	  data_ucelld=aurostd::substring2string(vdata1.at(i),"Cell Str1: ",FALSE);
	  ssfile << "UCELLD " << data_ucelld << endl;
	}
      }
      // aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "v_atom=" << data_v_atom;
      // aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "ucelld=" << data_ucelld;
      vector<double> vnbondxx_OLD;
      uint _nspecies=data.nspecies;
      for(uint isp1=0;isp1<_nspecies;isp1++)
	for(uint isp2=isp1;isp2<_nspecies;isp2++)
	  for(uint i=0;i<vdata1.size();i++) {
	    string string2find,string2search="Pairs between types: "+aurostd::utype2string(isp1)+" "+aurostd::utype2string(isp2);
	    tokens.clear();tokens.push_back("");
	    if(aurostd::substring2bool(vdata1.at(i),string2search)) {
	      string2find=aurostd::substring2string(vdata1.at(i),string2search,FALSE);
	      aurostd::string2tokens(string2find,tokens);
	      vnbondxx_OLD.push_back(aurostd::string2utype<double>(tokens.at(0))); // will do better but for now it is OK
	      //[CO 171024 OBSOLETE]data.vnbondxx.push_back(aurostd::string2utype<double>(tokens.at(0))); // will do better but for now it is OK
	    }
	  }

      if(data.vnbondxx.size())
	for(uint i=0;i<(uint) data.vnbondxx.size();i++)
	  data.nbondxx+=aurostd::utype2string<double>(data.vnbondxx.at(i),5)+(i<data.vnbondxx.size()-1?",":"");
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NBONDXX_OLD = " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vnbondxx_OLD,6),',') << endl; // CO 171025
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NBONDXX = " << data.nbondxx << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " MIN_DIST = " << SYM::minimumDistance(str_relax) << endl; // CO 171025

      if(aurostd::abs((double) data.vnbondxx.size()-_nspecies*(_nspecies+1.0)/2.0)>0.1)
	{ cerr << MESSAGE << " incompatible data.vnbondxx.size()[" << data.vnbondxx.size() << "]!=_nspecies*(_nspecies+1)/2)[" << (_nspecies*(_nspecies+1.0)/2.0) << "]" << endl;exit(0); }
      if(flag_LIB2==TRUE && _nspecies==1) { _nspecies++;data.vnbondxx.push_back(0.0);data.vnbondxx.push_back(0.0); } // FIX for purity control and alien generation
      for(uint isp1=0,inbondxx=0;isp1<_nspecies;isp1++)
	for(uint isp2=isp1;isp2<_nspecies;isp2++)
	  ssfile << "BOND_" << char('A'+isp1) << char('A'+isp2) << " " << data.vnbondxx.at(inbondxx++) << " [norm V_ATOM^0.33]" << endl;
      // aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_VOLDISTPARAMS_OUT);
    }
    // PRINT FORCES/POSITIONS_CARTESIAN

    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [12]" << endl;
    // VOLDISTEvolution
    if(flag_VOLDISTEVOLUTION) {  // CONTCAR.relax
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [VOLDISTEvolution]" << endl;
      stringstream ssfile;
      xstructure str;
      bool VOLDISTEvolution_flag=FALSE;//TRUE;
      for(int istep=1;istep<=3;istep++) {
	if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [AUID=] istep=" << istep << endl;
	if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [start istep] " << istep << endl;
	if(istep==1) str=str_orig;
	if(istep==2) str=str_relax1;
	if(istep==3) str=str_relax;
	if(istep==1) {
	  if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [Volume Distance Bonds Analyzer: ORIGINAL]" << endl;
	  ssfile << "Volume Distance Bonds Analyzer: ORIGINAL" << endl;
	  stringstream sss;vector<string> vline;
	  pflow::PrintData(xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO),sss,"DATA"); // DATA
	  aurostd::string2vectorstring(sss.str(),vline);
	  for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"real space volume")) ssfile << vline.at(iline) << endl;
	  for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"Real space a b c alpha beta gamma")) ssfile << vline.at(iline) << endl;
	}
	if(istep==2) {
	  if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [Volume Distance Bonds Analyzer: RELAX1]" << endl;
	  ssfile << "Volume Distance Bonds Analyzer: RELAX1" << endl;
	  stringstream sss;vector<string> vline;
	  pflow::PrintData(xstructure(directory_RAW+"/CONTCAR.relax1",IOAFLOW_AUTO),sss,"DATA"); // DATA
	  aurostd::string2vectorstring(sss.str(),vline);
	  for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"real space volume")) ssfile << vline.at(iline) << endl;
	  for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"Real space a b c alpha beta gamma")) ssfile << vline.at(iline) << endl;
	}
	if(istep==3) {
	  if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [Volume Distance Bonds Analyzer: RELAX]" << endl;
	  ssfile << "Volume Distance Bonds Analyzer: RELAX" << endl;
	  stringstream sss;vector<string> vline;
	  pflow::PrintData(xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO),sss,"DATA"); // DATA
	  aurostd::string2vectorstring(sss.str(),vline);
	  for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"real space volume")) ssfile << vline.at(iline) << endl;
	  for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"Real space a b c alpha beta gamma")) ssfile << vline.at(iline) << endl;
	}
	if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [stop istep] " << istep << endl;
	if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << str << endl;
	if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << ssfile.str() << endl;

	vector<string> aus_data_nbondxx;
	string data1=pflow::PrintData1(str,-1.0);
	vector<string> vdata1;aurostd::string2vectorstring(data1,vdata1);
	uint _nspecies=data.nspecies;
	if(VOLDISTEvolution_flag) cerr << _nspecies << " " << data1.size() << " " << vdata1.size() << endl << data1 << endl;
	for(uint isp1=0;isp1<_nspecies;isp1++)
	  for(uint isp2=isp1;isp2<_nspecies;isp2++)
	    for(uint i=0;i<vdata1.size();i++) {
	      if(VOLDISTEvolution_flag) cerr << isp1 << " " << isp2 << " " << i << endl;
	      string string2find,string2search="Pairs between types: "+aurostd::utype2string(isp1)+" "+aurostd::utype2string(isp2);
	      tokens.clear();tokens.push_back("");
	      if(aurostd::substring2bool(vdata1.at(i),string2search)) {
		string2find=aurostd::substring2string(vdata1.at(i),string2search,FALSE);
		aurostd::string2tokens(string2find,tokens);
		aus_data_nbondxx.push_back(tokens.at(0));
	      }
	    }
	if(aurostd::abs((double) aus_data_nbondxx.size()-_nspecies*(_nspecies+1.0)/2.0)>0.1)
	  { cerr << MESSAGE << " incompatible aus_data_nbondxx.size()[" << aus_data_nbondxx.size() << "]!=_nspecies*(_nspecies+1)/2)[" << (_nspecies*(_nspecies+1.0)/2.0) << "]" << endl;exit(0); }
	if(flag_LIB2==TRUE && _nspecies==1) { _nspecies++;aus_data_nbondxx.push_back("");aus_data_nbondxx.push_back(""); } // FIX for purity control and alien generation
	for(uint isp1=0,inbondxx=0;isp1<_nspecies;isp1++)
	  for(uint isp2=isp1;isp2<_nspecies;isp2++)
	    ssfile << "BOND_" << char('A'+isp1) << char('A'+isp2) << " " << aus_data_nbondxx.at(inbondxx++) << " [norm V_ATOM^0.33]" << endl;
      }
      //   aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_VOLDISTEVOLUTION_OUT);
    } // END VOLDISTEvolution
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [13]" << endl;
    // check for ERRORS
    if(flag_ENERGY1) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [flag_ENERGY1]" << endl;
      if(aurostd::abs(data.energy_atom-data1_energy_atom)>ENERGY_ATOM_ERROR_meV/1000.0)
	if(data.energy_atom<0.0 && data1_energy_atom>0.0) flag_ERROR=TRUE;
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [14]" << endl;
    // done now WRITE aflowlib.out
    if(flag_ERROR) data.error_status="ERROR";

    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [15]" << endl;

    // PERFORM EDATA DATA AND CIF STEPS -------------------------------------------------------------------------
    xstructure str,str_sp,str_sc;
    vector<xstructure> vcif;

    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [16]" << endl;
    // PERFORM EDATA STEP
    if(flag_EDATA_ORIG_ || flag_EDATA_RELAX_) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA start: " << directory_RAW << endl;
    if(flag_EDATA_ORIG_) { // ORIG
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing orig (POSCAR.orig) text format: " << directory_RAW << endl;
	// [OBSOLETE] aurostd::execute("cd "+directory_RAW+" && cat POSCAR.orig | aflow --edata > "+DEFAULT_FILE_EDATA_ORIG_OUT);
	str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
	// DX - START
	//DX 20180526 [OBSOLETE] str.directory = str_sp.directory = str_sc.directory = aflags.Directory;
	// DX - END
	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
	xstructure str_sym=str; // CO171027 // DX 2/26/18 - set equal str
	pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"EDATA","txt",false); // 1=EDATA // CO 171025 // CO171027
	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT);
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing orig (POSCAR.orig) json format: " << directory_RAW << endl;
	stringstream jjj; // CO 171025
	pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"EDATA","json",true); // 1=EDATA, already_calculated!  // CO 171025 // CO171027
	aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_JSON); // CO 171025
	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
	// CO and DX START 170713 - adding symmetry output to RAW
	string new_sym_file;
	//txt variants
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,new_sym_file);
	} // DX 12/7/17 - added pgroupk_xtal
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,new_sym_file);
	} // CO 171025
	//json variants
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,new_sym_file);
	} // DX 12/7/17 - added pgroupk_xtal
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,"orig",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,new_sym_file);
	} // CO 171025
	//cout << message;
	// CO and DX STOP 170713 - adding symmetry output to RAW
	// now extract info
	vector<string> vline_edata;
	aurostd::string2vectorstring(sss.str(),vline_edata);
	for(uint iline=0;iline<vline_edata.size();iline++) {
	  if(data.Bravais_lattice_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.Bravais_lattice_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.lattice_variation_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice Variation")) { // Bravais_Lattice_Variation
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.lattice_variation_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.lattice_system_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice System")) { // Lattice_System
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.lattice_system_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.Pearson_symbol_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata.at(iline),"Superlattice")) { // Pearson
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.Pearson_symbol_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	}
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_orig << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_orig << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_orig << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_orig << endl;
      }
    }
    if(flag_EDATA_RELAX_) { // RELAX
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing relax (CONTCAR.relax) text format: " << directory_RAW << endl;
	// [OBSOLETE] aurostd::execute("cd "+directory_RAW+" && cat CONTCAR.relax | aflow --edata > "+DEFAULT_FILE_EDATA_RELAX_OUT);
	str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
	xstructure str_sym=str; // CO171027 // DX 2/26/18 - set equal str
	pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"EDATA","txt",false); // EDATA // CO 171025 // CO171027
	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT);
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing relax (CONTCAR.relax) json format: " << directory_RAW << endl;
	stringstream jjj; // CO 171025
	pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"EDATA","json",true); // EDATA already_calculated! // CO 171025 // CO171027
	aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_JSON); // CO 171025
	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
	// CO and DX START 170713 - adding symmetry output to RAW
	string new_sym_file;
	//txt variants
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,new_sym_file);
	} // DX 12/7/17 - added pgroupk_xtal
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,new_sym_file);
	} // CO 171025
	//json variants
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,new_sym_file);
	} // DX 12/7/17 - added pgroupk_xtal
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,"relax",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,new_sym_file);
	} // CO 171025
	//cout << message;
	// CO and DX STOP 170713 - adding symmetry output to RAW
	// now extract info
	vector<string> vline_edata;
	aurostd::string2vectorstring(sss.str(),vline_edata);
	for(uint iline=0;iline<vline_edata.size();iline++) {
	  if(data.Bravais_lattice_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.Bravais_lattice_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.lattice_variation_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice Variation")) { // Bravais_Lattice_Variation
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.lattice_variation_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.lattice_system_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice System")) { // Lattice_System
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.lattice_system_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.Pearson_symbol_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata.at(iline),"Superlattice")) { // Pearson
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.Pearson_symbol_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	}
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_relax << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_relax << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_relax << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_relax << endl;
      }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [17]" << endl;
    if(flag_EDATA_BANDS_) { // BANDS
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing bands (POSCAR.bands) text format: " << directory_RAW << endl;
	// [OBSOLETE] aurostd::execute("cd "+directory_RAW+" && cat POSCAR.bands | aflow --edata > "+DEFAULT_FILE_EDATA_BANDS_OUT);
	str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
	xstructure str_sym=str; // CO171027 // DX 2/26/18 - set equal str
	pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"EDATA","txt",false); // EDATA // CO 171025 // CO171027
	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT);
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing bands (POSCAR.bands) json format: " << directory_RAW << endl;
	stringstream jjj; // CO 171025
	pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"EDATA","json",true); // EDATA already_calculated! // CO 171025 // CO171027
	aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_JSON); // CO 171025
	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
	// CO and DX START 170713 - adding symmetry output to RAW
	string new_sym_file;
	//txt variants
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,new_sym_file);
	} // DX 12/7/17 - added pgroupk_xtal
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"txt");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,new_sym_file);
	} // CO 171025
	//json variants
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,new_sym_file);
	} // DX 12/7/17 - added pgroupk_xtal
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,new_sym_file);
	} // CO 171025
	KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"json");
	if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON)) {
	  AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,"bands",new_sym_file);
	  aurostd::RenameFile(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,new_sym_file);
	} // CO 171025
	//cout << message;
	// CO and DX STOP 170713 - adding symmetry output to RAW
	// now extract info
	vector<string> vline_edata;
	aurostd::string2vectorstring(sss.str(),vline_edata);
	for(uint iline=0;iline<vline_edata.size();iline++) {
	  if(data.Bravais_lattice_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.Bravais_lattice_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.lattice_variation_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice Variation")) { // Bravais_Lattice_Variation
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.lattice_variation_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.lattice_system_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice System")) { // Lattice_System
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.lattice_system_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	  if(data.Pearson_symbol_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata.at(iline),"Superlattice")) { // Pearson
	    aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
	    data.Pearson_symbol_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
	}
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_relax << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_relax << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_relax << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_relax << endl;
      }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [18]" << endl;
    // PERFORM DATA STEP
    if(flag_DATA_ORIG_ || flag_DATA_RELAX_) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA start: " << directory_RAW << endl;
    if(flag_DATA_ORIG_) { // ORIG
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing orig (POSCAR.orig) text format: " << directory_RAW << endl;
	str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
	xstructure str_sym=str; // CO171027 // DX 2/26/18 - set equal str
	pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"DATA","txt",false); // DATA // CO 171025 // CO171027
	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT);
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing orig (POSCAR.orig) json format: " << directory_RAW << endl;
	stringstream jjj; // CO 171025
	pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"DATA","json",true); // DATA already_calculated! // CO 171025 // CO171027
	aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_JSON); // CO 171025
      }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [19]" << endl;
    if(flag_DATA_RELAX_) { // RELAX
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (CONTCAR.relax) text format: " << directory_RAW << endl;
	str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
	xstructure str_sym=str; // CO171027 // DX 2/26/18 - set equal str
	pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"DATA","txt",false); // DATA // CO 171025 // CO171027
	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT);
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (CONTCAR.relax) json format: " << directory_RAW << endl;
	stringstream jjj; // CO 171025
	pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"DATA","json",true); // DATA already_calculated! // CO 171025 // CO171027
	aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_JSON); // CO 171025
      }
    }
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [20]" << endl;
    if(flag_DATA_BANDS_) { // BANDS
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (POSCAR.bands) text format: " << directory_RAW << endl;
	str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.Clear();str_sc.Clear();
	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
	xstructure str_sym=str; // CO171027 // DX 2/26/18 - set equal str
	pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"DATA","txt",false); // DATA // CO 171025 // CO171027
	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT);
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (POSCAR.bands) json format: " << directory_RAW << endl;
	stringstream jjj; // CO 171025
	pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"DATA","json",true); // DATA already_calculated! // CO 171025 // CO171027
	aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_JSON); // CO 171025
      }
    }

    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Thermodynamics [21]" << endl;
    // CREATING CIFS
    //      cout << MESSAGE << " CIF creation: " << directory_LIB << endl;
    if(vcif.size()==0 && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) vcif.push_back(xstructure(directory_RAW+"/POSCAR.bands",IOVASP_AUTO));
    xvector<double> nvec(3);nvec(1)=1;nvec(2)=1;nvec(3)=1;
    double angle=45;

    for (uint j=0;j<vcif.size();j++) {
      vcif.at(j)=GetLTFVCell(nvec,angle,vcif.at(j));
      if(j==0) cout << MESSAGE << " CIF creation: data.spacegroup_relax=" << data.spacegroup_relax << endl;
      if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal" << endl;
      if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim" << endl;
      if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv" << endl;
      stringstream oss;
      //    for(uint i=0;i<vspecies.size();i++) cerr << vspecies.at(i) << endl;
      vcif.at(j).species.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species.push_back(vspecies.at(i));
      vcif.at(j).species_pp.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp.push_back(vspecies.at(i));
      vcif.at(j).species_pp_type.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_type.push_back("");
      vcif.at(j).species_pp_version.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_version.push_back("");
      vcif.at(j).species_pp_ZVAL.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_ZVAL.push_back(0.0);
      vcif.at(j).species_pp_vLDAU.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_vLDAU.push_back(deque<double>());
      for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).name=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
      for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).cleanname=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
      pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
      if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+".cif");
      if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim.cif");
      if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv.cif");
      vcif.at(j).AddCorners();
      oss.clear();oss.str("");
      pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
      if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal_corner" << endl;
      if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim_corner" << endl;
      if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv_corner" << endl;
      if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_corner.cif");
      if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim_corner.cif");
      if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv_corner.cif");
    }

    // CREATING GEOMETRY
    if(aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " LOAD SPECIES=" << data.species << endl;
      xstructure xstr(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);
      deque<string> vspecies;
      aurostd::string2tokens(data.species,vspecies,",");
      xstr.SetSpecies(vspecies);
      stringstream sss_vasp;sss_vasp << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for VASP" << endl;
      aurostd::stringstream2file(sss_vasp,directory_RAW+"/CONTCAR.relax.vasp");
      //	cerr << xstr << endl;
      xstr.xstructure2qe();
      stringstream sss_qe;sss_qe << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for QE" << endl;
      aurostd::stringstream2file(sss_qe,directory_RAW+"/CONTCAR.relax.qe");
      //	cerr << xstr << endl;
      xstr.xstructure2abinit();
      stringstream sss_abinit;sss_abinit << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for ABINIT" << endl;
      aurostd::stringstream2file(sss_abinit,directory_RAW+"/CONTCAR.relax.abinit");
      //	cerr << xstr << endl;
      xstr.xstructure2aims();
      stringstream sss_aims;sss_aims << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for AIMS" << endl;
      aurostd::stringstream2file(sss_aims,directory_RAW+"/CONTCAR.relax.aims");
      //	cerr << xstr << endl;
    }

    if(flag_LIB2==TRUE) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PATCHING for LIB2: " << directory_RAW << endl;
      if(_APENNSY_STYLE_OLD_) {
	for(uint iext=0;iext<vext.size();iext++) {
	  fileA_LIB=directory_LIB+"/aflow.agroup.out"+vext.at(iext);
	  fileA_RAW=directory_RAW+"/aflow.agroup.out"+vext.at(iext);
	  if(aurostd::FileExist(fileA_LIB)) {
	    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.agroup.out"+vext.at(iext),directory_RAW,"aflow.agroup.out"+vext.at(iext),vfile,MESSAGE);
	  }
	}
      }
      if(_APENNSY_STYLE_OLD_) {
	for(uint iext=0;iext<vext.size();iext++) {
	  fileA_LIB=directory_LIB+"/aflow.fgroup.out"+vext.at(iext);
	  fileA_RAW=directory_RAW+"/aflow.fgroup.out"+vext.at(iext);
	  if(aurostd::FileExist(fileA_LIB)) {
	    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.fgroup.out"+vext.at(iext),directory_RAW,"aflow.fgroup.out"+vext.at(iext),vfile,MESSAGE);
	  }
	}
      }
      if(_APENNSY_STYLE_OLD_) {
	for(uint iext=0;iext<vext.size();iext++) {
	  fileA_LIB=directory_LIB+"/aflow.pgroup.out"+vext.at(iext);
	  fileA_RAW=directory_RAW+"/aflow.pgroup.out"+vext.at(iext);
	  if(aurostd::FileExist(fileA_LIB)) {
	    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.pgroup.out"+vext.at(iext),directory_RAW,"aflow.pgroup.out"+vext.at(iext),vfile,MESSAGE);
	  }
	}
      }
      if(_APENNSY_STYLE_OLD_) {
	for(uint iext=0;iext<vext.size();iext++) {
	  fileA_LIB=directory_LIB+"/aflow.iatoms.out"+vext.at(iext);
	  fileA_RAW=directory_RAW+"/aflow.iatoms.out"+vext.at(iext);
	  if(aurostd::FileExist(fileA_LIB)) {
	    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.iatoms.out"+vext.at(iext),directory_RAW,"aflow.iatoms.out"+vext.at(iext),vfile,MESSAGE);
	  }
	}
      }
      for(uint iext=0;iext<vext.size();iext++) {
	fileA_LIB=directory_LIB+"/aflow.qmvasp.out"+vext.at(iext);
	fileA_RAW=directory_RAW+"/aflow.qmvasp.out"+vext.at(iext);
	if(aurostd::FileExist(fileA_LIB)) {
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.qmvasp.out"+vext.at(iext),directory_RAW,"aflow.qmvasp.out"+vext.at(iext),vfile,MESSAGE);
	}
      }
      // for(uint iext=0;iext<vext.size();iext++) {
      // no vasp.out.relax1  fileA_LIB=directory_LIB+"/vasp.out.relax1"+vext.at(iext);fileA_RAW=directory_RAW+"/vasp.out.relax1"+vext.at(iext);if(aurostd::FileExist(fileA_LIB)) { aflowlib::LIB2RAW_FileNeeded(directory_LIB,"vasp.out.relax1"+vext.at(iext),directory_RAW,"vasp.out.relax1"+vext.at(iext),vfile,MESSAGE); }
      // no vasp.out.relax2  fileA_LIB=directory_LIB+"/vasp.out.relax2"+vext.at(iext);fileA_RAW=directory_RAW+"/vasp.out.relax2"+vext.at(iext);if(aurostd::FileExist(fileA_LIB)) { aflowlib::LIB2RAW_FileNeeded(directory_LIB,"vasp.out.relax2"+vext.at(iext),directory_RAW,"vasp.out.relax2"+vext.at(iext),vfile,MESSAGE); }
      // }
      if(_APENNSY_STYLE_OLD_) {
	for(uint iext=0;iext<vext.size();iext++) {
	  fileX_LIB=directory_LIB+"/POSCAR.relax1"+vext.at(iext);
	  fileX_RAW=directory_RAW+"/POSCAR.relax1"+vext.at(iext);
	  if(aurostd::FileExist(fileX_LIB)) {
	    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1"+vext.at(iext),directory_RAW,"POSCAR.relax1"+vext.at(iext),vfile,MESSAGE);
	  }
	}
      }
      for(uint iext=0;iext<vext.size();iext++) {
	fileX_LIB=directory_LIB+"/POSCAR.relax2"+vext.at(iext);
	fileX_RAW=directory_RAW+"/POSCAR.relax2"+vext.at(iext);
	if(aurostd::FileExist(fileX_LIB)) {
	  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax2"+vext.at(iext),directory_RAW,"POSCAR.relax2"+vext.at(iext),vfile,MESSAGE);
	}
      }
      aurostd::CopyFile(directory_RAW+"/CONTCAR.relax",directory_RAW+"/CONTCAR.relax2");vfile.push_back("CONTCAR.relax2");
      // aurostd::CopyFile(directory_LIB+"/OSZICAR.relax1"+vext.at(iext),directory_RAW+"/OSZICAR.relax1"+vext.at(iext));vfile.push_back("OSZICAR.relax1");
      // aurostd::CopyFile(directory_LIB+"/OSZICAR.relax2"+vext.at(iext),directory_RAW+"/OSZICAR.relax2"+vext.at(iext));vfile.push_back("OSZICAR.relax2");
      if(_APENNSY_STYLE_OLD_) {
	aurostd::CopyFile(directory_RAW+"/OSZICAR.relax",directory_RAW+"/OSZICAR.relax2");vfile.push_back("OSZICAR.relax2");
      }
      if(_APENNSY_STYLE_OLD_) {
	stringstream aus_exec;
	aus_exec << "cd " << directory_RAW << endl;
	aus_exec << "bin2ascii ./OSZICAR*" << endl << "bin2ascii ./OSZICAR*" << endl;
	aus_exec << "subst \"\\n\\n\" \"\\n\" ./OSZICAR*" << endl << "subst \"\\n\\n\" \"\\n\" ./OSZICAR*" << endl;
	aus_exec << "rm -f ./*~" << endl;
	aurostd::execute(aus_exec);
      }
    }
    // COMPRESS
    // aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/OSZICAR.relax");

    // LINKING
    if(FileName_OUTCAR_relax!="") {
      string FROM=FileName_OUTCAR_relax;
      string TO=FileName_OUTCAR_relax;
      aurostd::StringSubst(TO,"LIB/","RAW/");
      aurostd::StringSubst(TO,"relax1","relax");
      aurostd::StringSubst(TO,"relax2","relax");
      aurostd::StringSubst(TO,"relax3","relax");
      aurostd::StringSubst(TO,"static","relax");
      aurostd::StringSubst(TO,"bands","relax");
      cout << MESSAGE << " linking " << FROM << "->" << TO << endl;
      aurostd::execute("ln -sf "+FROM+" "+TO);
      //      exit(0);
    }
    
    // DONE
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Thermodynamics - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Magnetic
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_Magnetic(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE) {
    // CO 180130 - note that here we extract STATIC/BANDS properties only
    //some spin properties (spin/cell, spin/atom, spinD) can be extracted from relax2
    //we do this in the thermo loop and DO NOT attempt to redo here
    //spinF is nonzero for magnetic metals, so we need to determine Egap here as well (need BANDS)
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Magnetic [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Magnetic - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=magnetic";
    data.vloop.push_back("magnetic");
    xOUTCAR outcar_bands,outcar_static;
    xDOSCAR doscar;

    double MAG_EPS=1e-6;

    //try to grab STATIC properties first, Efermi and spin/cell, spin/atom, spinD
    //NB, spin properties from relax2 are okay (lower kppra), we grab it in thermo loop, so we don't attempt to here
    double EFERMI=AUROSTD_NAN;
    outcar_static.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"OUTCAR.static") || aurostd::EFileExist(directory_LIB+"/"+"OUTCAR.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.static",directory_RAW,"OUTCAR.static",vfile,MESSAGE);  // OUTCAR.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"OUTCAR.static") << endl;
      if(outcar_static.GetPropertiesFile(directory_RAW+"/"+"OUTCAR.static")) {
	EFERMI=outcar_static.Efermi;
	data.spin_cell=outcar_static.mag_cell;
	data.spin_atom=outcar_static.mag_atom;
	data.spinD="";
	data.vspinD.clear();
	if(outcar_static.vmag.size()) {
	  for(uint i=0;i<(uint) outcar_static.vmag.size();i++) {
	    data.spinD+=aurostd::utype2string<double>(outcar_static.vmag.at(i),5)+(i<outcar_static.vmag.size()-1?",":"");
	    data.vspinD.push_back(outcar_static.vmag.at(i));
	  }
	} else {
	  for(uint i=0;i<outcar_static.natoms;i++) {  //use outcar_static.natoms as there can be a primitivization between relax and static
	    data.spinD+=aurostd::utype2string<double>(0)+(i<outcar_static.natoms-1?",":"");
	    data.vspinD.push_back(0.0);
	  }
	}
      } else { cout << MESSAGE << " ERROR OUTCAR.static properties cannot be extracted: " << outcar_static.ERROR << endl; }
    } else { cout << MESSAGE << " MISSING OUTCAR.static" << endl; }
    if(EFERMI==AUROSTD_NAN) { cout << MESSAGE << " unable to load OUTCAR.static, using Efermi from bands" << endl; }

    //ideally we grab both bands and static
    //we want static because we trust this Efermi (self-consistent)
    //however, we fall back on Efermi from bands

    // BANDS CALCULATION
    data.spinF=0.0;// DEFAULT
    outcar_bands.clear();
    doscar.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"OUTCAR.bands") || aurostd::EFileExist(directory_LIB+"/"+"OUTCAR.bands")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.bands",directory_RAW,"OUTCAR.bands",vfile,MESSAGE);  // OUTCAR.bands
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"OUTCAR.bands") << endl;
      if(outcar_bands.GetPropertiesFile(directory_RAW+"/"+"OUTCAR.bands")) {
	outcar_bands.GetBandGap(EFERMI);
	data.Egap=outcar_bands.Egap_net;
	data.Egap_fit=outcar_bands.Egap_fit_net;
	data.Egap_type=outcar_bands.Egap_type_net;
	if(aurostd::substring2bool(data.Egap_type,"metal")) data.Egap=0.0;      //half-metal
	if(aurostd::substring2bool(data.Egap_type,"metal")) data.Egap_fit=0.0;  //half-metal
   
	// SPIN POLARIZATION AT FERMI LEVEL
	if(data.Egap<MAG_EPS && aurostd::abs(data.spin_cell)>MAG_EPS) { //must be metal and magnetic
	  if(doscar.content=="" && (aurostd::FileExist(directory_LIB+"/"+"DOSCAR.static") || aurostd::EFileExist(directory_LIB+"/"+"DOSCAR.static"))) {
	    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"DOSCAR.static",directory_RAW,"DOSCAR.static",vfile,MESSAGE);  // DOSCAR.static
	    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"DOSCAR.static") << endl;
	    doscar.GetPropertiesFile(directory_RAW+"/"+"DOSCAR.static");
	  }
	  if(doscar.content=="" && (aurostd::FileExist(directory_LIB+"/"+"DOSCAR.relax2") || aurostd::EFileExist(directory_LIB+"/"+"DOSCAR.relax2"))) {
	    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"DOSCAR.relax2",directory_RAW,"DOSCAR.relax",vfile,MESSAGE);  // DOSCAR.relax2
	    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"DOSCAR.relax") << endl;
	    doscar.GetPropertiesFile(directory_RAW+"/"+"DOSCAR.relax");
	  }
	  if(!doscar.content.empty()) {
	    data.spinF=doscar.spinF;
	  } else {
	    cout << MESSAGE << " MISSING DOSCAR.static and DOSCAR.relax[2]" << endl;
	  }
	}
      } else { cout << MESSAGE << " ERROR OUTCAR.static properties cannot be extracted: " << outcar_static.ERROR << endl; }
    } else { cout << MESSAGE << " MISSING OUTCAR.static" << endl; }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spin/cell = "         <<   data.spin_cell << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spin/atom = "         <<   data.spin_atom << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spinD     = "         <<   data.spinD << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spinF     = "         << ((data.spinF   !=AUROSTD_NAN)?aurostd::utype2string(data.spinF,    5):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " Egap (eV)     = "     << ((data.Egap    !=AUROSTD_NAN)?aurostd::utype2string(data.Egap,     5):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " Egap_fit (eV) = " << ((data.Egap_fit!=AUROSTD_NAN)?aurostd::utype2string(data.Egap_fit, 5):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " Egap_type     = "     << ((data.Egap_type.size())?data.Egap_type:"unavailable") << endl;

    // done
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Magnetic - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Bader  // COREY WORK HERE
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_Bader(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_Bader [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Bader - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=bader";
    data.vloop.push_back("bader");

    // [OBSOLETE]   xOUTCAR outcar;
    // [OBSOLETE]aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.static",directory_RAW,"OUTCAR.static",vfile,MESSAGE);  // OUTCAR.static
    // [OBSOLETE]if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"OUTCAR.static") << endl;
    // [OBSOLETE]outcar.GetPropertiesFile(directory_RAW+"/"+"OUTCAR.static");
    if(aurostd::FileExist(directory_LIB+"/"+"CHGCAR.static") || aurostd::EFileExist(directory_LIB+"/"+"CHGCAR.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CHGCAR.static",directory_RAW,"CHGCAR.static",vfile,MESSAGE);  // CHGCAR.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"CHGCAR.static") << endl;
    } else {
      cout << MESSAGE << " MISSING CHGCAR.static" << endl;
      return FALSE;
    }
    if(aurostd::FileExist(directory_LIB+"/"+"AECCAR0.static") || aurostd::EFileExist(directory_LIB+"/"+"AECCAR0.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AECCAR0.static",directory_RAW,"AECCAR0.static",vfile,MESSAGE);  // AECCAR0.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"AECCAR0.static") << endl;
    } else {
      cout << MESSAGE << " MISSING AECCAR0.static" << endl;
      return FALSE;
    }
    if(aurostd::FileExist(directory_LIB+"/"+"AECCAR2.static") || aurostd::EFileExist(directory_LIB+"/"+"AECCAR2.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AECCAR2.static",directory_RAW,"AECCAR2.static",vfile,MESSAGE);  // AECCAR2.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"AECCAR2.static") << endl;
    } else {
      cout << MESSAGE << " MISSING AECCAR2.static" << endl;
      return FALSE;
    }

    deque<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");
    for(uint i=0;i<vspecies.size();i++)
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << vspecies.at(i) << endl;

    vector<double> vspecies_pp_ZVAL=data.vspecies_pp_ZVAL;//;aurostd::string2tokens(data.species_pp_ZVAL,vspecies_pp_ZVAL,",");
    for(uint i=0;i<vspecies_pp_ZVAL.size();i++)
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << vspecies_pp_ZVAL.at(i) << endl;

    deque<int> num_each_type;
    aurostd::string2tokens<int>(data.composition,num_each_type,",");

    vector<double> vZVAL;
    for(uint i=0;i<vspecies_pp_ZVAL.size();i++) {
      for(uint j=0;j<(uint)num_each_type.at(i);j++) {
	vZVAL.push_back(vspecies_pp_ZVAL.at(i));
	//vZVAL.push_back(aurostd::string2utype<double>(vspecies_pp_ZVAL.at(i)));
      }
    }

    //flags
    aurostd::xoption bader_flags;
    // DX and CO - START
    bader_flags.flag("BADER::AFLOWLIB_LIBRARY",TRUE); //net charge file
    bader_flags.flag("BADER::JVXL_ALL_SPECIES",TRUE); //for flag automation, no need to specify cutoffs,downsample_ratios here
    bader_flags.flag("BADER::JVXL_CYCLIC",TRUE);      //CYCLIC MODE
    //bader_flags.flag("BADER::JVXL_CYCLIC",FALSE);   //SETS MODE
    // DX and CO - END

    vector<double> cutoffs;
    vector<int> downsample_ratios;
    //cutoffs.push_back(0.1);
    cutoffs.push_back(0.2);
    cutoffs.push_back(0.3);
    cutoffs.push_back(0.4);
    cutoffs.push_back(0.5);
    //cutoffs.push_back(0.75);
    //downsample_ratios.push_back(1);
    downsample_ratios.push_back(2);   //for all
    /*if(data.natoms<=7) {              //for 0.5 and 0.75, downsample_ratio depends on natoms
      downsample_ratios.push_back(1);
      //downsample_ratios.push_back(1);
      } else {                          //8 or more
      downsample_ratios.push_back(2);
      //downsample_ratios.push_back(2);
      }
      //downsample_ratios.push_back(1);*/
    string bader_options;
    ostringstream oss;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " beginning BADER calculation--please be patient" << endl;
    bader_functions::Flags2BaderCommands(bader_flags,bader_options,oss);
    // DX and CO - START
    bader_functions::BaderCalc(bader_flags,bader_options,data.prototype,vspecies,num_each_type,vZVAL,cutoffs,downsample_ratios,directory_RAW,oss);
    // DX and CO - END
    if(AFLOWLIB_VERBOSE) cout << oss.str();

    vector<string> vline,tokens;
    stringstream abader_ss;
    string abader_out=data.prototype+"_abader.out";

    abader_ss.clear();
    if(aurostd::FileExist(directory_RAW+"/"+abader_out)) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+abader_out) << endl;
      aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(directory_RAW+"/"+abader_out),abader_ss,"[BADER_RESULTS]START","[BADER_RESULTS]STOP");
      aurostd::stream2vectorstring(abader_ss,vline);
      for (uint i=0;i<vline.size();i++) {
	aurostd::StringSubst(vline.at(i),"="," ");
	aurostd::string2tokens(vline.at(i),tokens," ");
	if(tokens.size()>=2) {
	  if(tokens.at(0)=="bader_net_charges") data.bader_net_charges=tokens.at(1);
	  if(tokens.at(0)=="bader_atomic_volumes") data.bader_atomic_volumes=tokens.at(1);
	}
      }
      aurostd::string2tokens<double>(data.bader_net_charges,data.vbader_net_charges,",");         //be careful, sets precision to 20
      aurostd::string2tokens<double>(data.bader_atomic_volumes,data.vbader_atomic_volumes,",");   //be careful, sets precision to 20
      if(data.vbader_net_charges.size()==0) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " data.vbader_net_charges was not extracted correctly" << endl;
	return FALSE;
      }
      if(data.vbader_net_charges.size()!=data.vbader_atomic_volumes.size()) {
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " length of data.vbader_net_charges does not match length of data.vbader_atomic_volumes" << endl;
	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " there was a problem with the bader data extraction" << endl;
	return FALSE;
      }
      stringstream num_prec;
      if(AFLOWLIB_VERBOSE) {
	for(uint i=0;i<data.vbader_net_charges.size();i++) {
	  num_prec << std::fixed << setprecision(4) << data.vbader_net_charges.at(i);
	  cout << MESSAGE << " bader net charge for atom " << i+1 << " = " << num_prec.str() << endl;
	  num_prec.str("");
	}
	for(uint i=0;i<data.vbader_atomic_volumes.size();i++) {
	  num_prec << std::fixed << setprecision(4) << data.vbader_atomic_volumes.at(i);
	  cout << MESSAGE << " bader atomic volume for atom " << i+1 << " = " << num_prec.str() << endl;
	  num_prec.str("");
	}
      }
    } else {
      return FALSE;
    }
    // done
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << "aflowlib::LIB2RAW_Loop_Bader - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_AGL  // CORMAC
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_AGL(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_AGL [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AGL - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=agl";
    data.vloop.push_back("agl");

    vector<string> vline,tokens;
    stringstream aflow_agl_out;

    // AFLOW AGL
    aflow_agl_out.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"aflow.agl.out") || aurostd::EFileExist(directory_LIB+"/"+"aflow.agl.out")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.agl.out",directory_RAW,"aflow.agl.out",vfile,MESSAGE);  // aflow.agl.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL.out",directory_RAW,"AGL.out",vfile,MESSAGE);  // AGL.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_energies_temperature.out",directory_RAW,"AGL_energies_temperature.out",vfile,MESSAGE);  // AGL.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_thermal_properties_temperature.out",directory_RAW,"AGL_thermal_properties_temperature.out",vfile,MESSAGE);  // AGL.out
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"aflow.agl.out") << endl;
      aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(directory_RAW+"/"+"aflow.agl.out"),aflow_agl_out,"[AGL_RESULTS]START","[AGL_RESULTS]STOP");
      aurostd::stream2vectorstring(aflow_agl_out,vline);
      for (uint i=0;i<vline.size();i++) {
	aurostd::StringSubst(vline.at(i),"="," ");
	aurostd::string2tokens(vline.at(i),tokens," ");
	if(tokens.size()>=2) {
	  if(tokens.at(0)=="agl_thermal_conductivity_300K") data.agl_thermal_conductivity_300K=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_debye") data.agl_debye=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_acoustic_debye") data.agl_acoustic_debye=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_gruneisen") data.agl_gruneisen=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_heat_capacity_Cv_300K") data.agl_heat_capacity_Cv_300K=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_heat_capacity_Cp_300K") data.agl_heat_capacity_Cp_300K=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_thermal_expansion_300K") data.agl_thermal_expansion_300K=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_bulk_modulus_static_300K") data.agl_bulk_modulus_static_300K=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="agl_bulk_modulus_isothermal_300K") data.agl_bulk_modulus_isothermal_300K=aurostd::string2utype<double>(tokens.at(1));
	}
      }
    } else {
      return FALSE;
    }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_thermal_conductivity_300K (W/m*K) = " << ((data.agl_thermal_conductivity_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_thermal_conductivity_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_debye (K) = " << ((data.agl_debye!=AUROSTD_NAN)?aurostd::utype2string(data.agl_debye,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_acoustic_debye (K) = " << ((data.agl_acoustic_debye!=AUROSTD_NAN)?aurostd::utype2string(data.agl_acoustic_debye,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_gruneisen = " << ((data.agl_gruneisen!=AUROSTD_NAN)?aurostd::utype2string(data.agl_gruneisen,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_heat_capacity_Cv_300K (kB/cell) = " << ((data.agl_heat_capacity_Cv_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_heat_capacity_Cv_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_heat_capacity_Cp_300K (kB/cell) = " << ((data.agl_heat_capacity_Cp_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_heat_capacity_Cp_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_thermal_expansion_300K (1/K) = " << ((data.agl_thermal_expansion_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_thermal_expansion_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_bulk_modulus_static_300K (GPa) = " << ((data.agl_bulk_modulus_static_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_bulk_modulus_static_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_bulk_modulus_isothermal_300K (GPa) = " << ((data.agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_bulk_modulus_isothermal_300K,10):"unavailable") << endl;
    // done
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AGL - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_AEL  // CORMAC
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_AEL(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_AEL [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AEL - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=ael";
    data.vloop.push_back("ael");

    vector<string> vline,tokens;
    stringstream aflow_ael_out;

    // AFLOW AEL
    aflow_ael_out.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"aflow.ael.out") || aurostd::EFileExist(directory_LIB+"/"+"aflow.ael.out")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.ael.out",directory_RAW,"aflow.ael.out",vfile,MESSAGE);  // aflow.ael.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AEL_Elastic_constants.out",directory_RAW,"AEL_Elastic_constants.out",vfile,MESSAGE);  // AEL_Elastic_constants.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AEL_Compliance_tensor.out",directory_RAW,"AEL_Compliance_tensor.out",vfile,MESSAGE);  // AEL_Compliance_tensor.out
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"aflow.ael.out") << endl;
      aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(directory_RAW+"/"+"aflow.ael.out"),aflow_ael_out,"[AEL_RESULTS]START","[AEL_RESULTS]STOP");
      aurostd::stream2vectorstring(aflow_ael_out,vline);
      for (uint i=0;i<vline.size();i++) {
	aurostd::StringSubst(vline.at(i),"="," ");
	aurostd::string2tokens(vline.at(i),tokens," ");
	if(tokens.size()>=2) {
	  if(tokens.at(0)=="ael_poisson_ratio") data.ael_poisson_ratio=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="ael_bulk_modulus_voigt" || tokens.at(0)=="ael_bulk_modulus_voight") data.ael_bulk_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="ael_bulk_modulus_reuss") data.ael_bulk_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="ael_shear_modulus_voigt" || tokens.at(0)=="ael_shear_modulus_voigth") data.ael_shear_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="ael_shear_modulus_reuss") data.ael_shear_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="ael_bulk_modulus_vrh") data.ael_bulk_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="ael_shear_modulus_vrh") data.ael_shear_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
	  if(tokens.at(0)=="ael_elastic_anistropy") data.ael_elastic_anistropy=aurostd::string2utype<double>(tokens.at(1));
	}
      }
    } else {
      return FALSE;
    }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_poisson_ratio = " << ((data.ael_poisson_ratio!=AUROSTD_NAN)?aurostd::utype2string(data.ael_poisson_ratio,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_bulk_modulus_voigt (GPa) = " << ((data.ael_bulk_modulus_voigt!=AUROSTD_NAN)?aurostd::utype2string(data.ael_bulk_modulus_voigt,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_bulk_modulus_reuss (GPa) = " << ((data.ael_bulk_modulus_reuss!=AUROSTD_NAN)?aurostd::utype2string(data.ael_bulk_modulus_reuss,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_shear_modulus_voigt (GPa) = " << ((data.ael_shear_modulus_voigt!=AUROSTD_NAN)?aurostd::utype2string(data.ael_shear_modulus_voigt,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_shear_modulus_reuss (GPa) = " << ((data.ael_shear_modulus_reuss!=AUROSTD_NAN)?aurostd::utype2string(data.ael_shear_modulus_reuss,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_bulk_modulus_vrh (GPa) = " << ((data.ael_bulk_modulus_vrh!=AUROSTD_NAN)?aurostd::utype2string(data.ael_bulk_modulus_vrh,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_shear_modulus_vrh (GPa) = " << ((data.ael_shear_modulus_vrh!=AUROSTD_NAN)?aurostd::utype2string(data.ael_shear_modulus_vrh,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_elastic_anistropy = " << ((data.ael_elastic_anistropy!=AUROSTD_NAN)?aurostd::utype2string(data.ael_elastic_anistropy,10):"unavailable") << endl;
    // done
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AEL - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_LOCK
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_LOCK(string& directory_LIB,string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,string MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_LOCK [1]" << endl;
    // Stefano Curtarolo 2009-2010-2011-2012
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_LOCK - begin " << directory_LIB << endl;
    //  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=LOCK";
    data.vloop.push_back("lock");

    vector<string> vlock,tokens;
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // OUTCAR.static
    aurostd::file2vectorstring(directory_RAW+"/"+_AFLOWLOCK_,vlock) ;
    _XHOST aus_XHOST;
    // ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&data.aflow_version.empty();iline++)
      if(aurostd::substring2bool(vlock.at(iline),"NFS") && aurostd::substring2bool(vlock.at(iline),"(") && aurostd::substring2bool(vlock.at(iline),")")) {
	aurostd::string2tokens(vlock.at(iline),tokens);
	data.aflow_version=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1));
	aurostd::StringSubst(data.aflow_version,"(","");aurostd::StringSubst(data.aflow_version,")","");
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflow_version = " << ((data.aflow_version.size())?data.aflow_version:"unavailable") << endl;
    // XHOST.CPU_Model ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.CPU_Model.empty();iline++)
      if(aurostd::substring2bool(vlock.at(iline),"XHOST.CPU_Model") && aurostd::substring2bool(vlock.at(iline),":")) {
	aurostd::string2tokens(vlock.at(iline),tokens,":");
	aus_XHOST.CPU_Model=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1));
	data.node_CPU_Model=aus_XHOST.CPU_Model;
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.CPU_Model = " << ((aus_XHOST.CPU_Model.size())?aus_XHOST.CPU_Model:"unavailable") << endl;
    // XHOST.CPU_Cores ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.CPU_Cores==0;iline++)
      if(aurostd::substring2bool(vlock.at(iline),"XHOST.CPU_Cores") && aurostd::substring2bool(vlock.at(iline),":")) {
	aurostd::string2tokens(vlock.at(iline),tokens,":");
	aus_XHOST.CPU_Cores=aurostd::string2utype<int>(aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)));
	aus_XHOST.CPU_Cores=aurostd::string2utype<int>(aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)));
	data.node_CPU_Cores=aus_XHOST.CPU_Cores;
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.CPU_Cores = " << ((aus_XHOST.CPU_Cores)?aurostd::utype2string<int>(aus_XHOST.CPU_Cores):"unavailable") << endl;
    // XHOST.CPU_MHz ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.CPU_MHz.empty();iline++)
      if(aurostd::substring2bool(vlock.at(iline),"XHOST.CPU_MHz") && aurostd::substring2bool(vlock.at(iline),":")) {
	aurostd::string2tokens(vlock.at(iline),tokens,":");
	aus_XHOST.CPU_MHz=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1));
	data.node_CPU_MHz=ceil(aurostd::string2utype<double>(aus_XHOST.CPU_MHz));
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.CPU_MHz = " << ((aus_XHOST.CPU_MHz.size())?aus_XHOST.CPU_MHz:"unavailable") << endl;
    // XHOST.RAM_GB ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.RAM_GB<0.001;iline++)
      if(aurostd::substring2bool(vlock.at(iline),"XHOST.RAM_GB") && aurostd::substring2bool(vlock.at(iline),":")) {
	aurostd::string2tokens(vlock.at(iline),tokens,":");
	aus_XHOST.RAM_GB=ceil(aurostd::string2utype<double>(aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1))));
	data.node_RAM_GB=aus_XHOST.RAM_GB;
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.RAM_GB = " << ((aus_XHOST.RAM_GB)?aurostd::utype2string<double>(aus_XHOST.RAM_GB):"unavailable") << endl;
    // REMOVING  ---------------------------------------------------------------
    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
    vector<string> vremove,vremovej;
    aurostd::string2tokens("aflow.in~,POTCAR.relax1,POTCAR.relax2,POTCAR.relax3,POTCAR.static,POTCAR.bands,AECCAR1.static,AECCAR0.bands,AECCAR1.bands,AECCAR2.bands,AECCAR0.relax1,AECCAR1.relax1,AECCAR2.relax1,AECCAR0.relax2,AECCAR1.relax2,AECCAR2.relax2,",vremove,",");
    for(uint iext=0;iext<vext.size();iext++) {
      for(uint iremove=0;iremove<vremove.size();iremove++) {
	vremovej.clear();
	aurostd::string2vectorstring(aurostd::execute2string("find "+directory_LIB+" -name \""+vremove.at(iremove)+vext.at(iext)+"\""),vremovej);
	for(uint j=0;j<vremovej.size();j++) {
	  if(aurostd::FileExist(vremovej.at(j))) {
	    aurostd::RemoveFile(vremovej.at(j));
	    if(AFLOWLIB_VERBOSE) 
	      cout << MESSAGE << " REMOVED = " << aurostd::CleanFileName(vremovej.at(j)) << endl;
	  }
	}
      } // iremove
    } // iext

    vector<string> vremoveERROR;
    aurostd::string2tokens("aflow.error.rotmat,aflow.error.nbands,aflow.error.symprec,aflow.error.read_kpoints_rd_sym,aflow.error.ibzkpt,aflow.error.gamma_shift,aflow.error.mpich11,aflow.error.mpich139,aflow.error.nkxyz_ikptd,aflow.error.eddrmm,aflow.error.lreal,aflow.error.exccor,aflow.error.brmix,aflow.error.dav,aflow.error.edddav,aflow.error.efield_pead,aflow.error.zpotrf,aflow.error.zpotrf_potim,aflow.error.natoms,aflow.error.psmaxn,aflow.error.npar,aflow.error.nparc,aflow.error.nparn,aflow.error.npar_remove,aflow.error.csloshing,aflow.error.dentet",vremoveERROR,",");
    for(uint iext=0;iext<vext.size();iext++) {
      for(uint iremove=0;iremove<vremoveERROR.size();iremove++) {
	if(aurostd::FileExist(directory_LIB+"/"+vremoveERROR.at(iremove)+vext.at(iext))) {
	  aurostd::RemoveFile(directory_LIB+"/"+vremoveERROR.at(iremove)+vext.at(iext));
	}
      } // iremove
    } // iext
    
    // [OBSOLETE]  for(uint i=0;i<vremove.size();i++) {
    // [OBSOLETE]  if(aurostd::FileExist(directory_LIB+"/"+vremove.at(iremove)+vext.at(iext))) {
    // [OBSOLETE]	aurostd::RemoveFile(directory_LIB+"/"+vremove.at(iremove)+vext.at(iext));
    // [OBSOLETE]	if(AFLOWLIB_VERBOSE) cout << MESSAGE << " REMOVED = " << aurostd::CleanFileName(directory_LIB+"/"+vremove.at(iremove)+vext.at(iext)) << endl;
    // [OBSOLETE]   }
    // [OBSOLETE]  }
    // done   ---------------------------------------------------------------
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_LOCK - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::XPLUG
// ***************************************************************************
namespace aflowlib {
  bool XPLUG_Directory_ok(const string &directory,const uint &i,const uint &j,const int& ithread,stringstream &oss);
  void *_threaded_interface_XPLUG_Directory(void *ptr);

  typedef struct {
    int      itbusy;          // FOR XPLUG (ALL)
    bool     VERBOSE;         // FOR XPLUG (ALL)
    //  ostringstream  oss;   // printout
    deque<bool>   *vok;       // directory to check
    deque<string> *vdirs;     // FOR
    deque<uint>   *vrun;      // index to run
    int      ITHREAD;         // FOR
    int      THREADS_MAX;     // FOR
  } _threaded_XPLUG_params;

  //_threaded_XPLUG_params params[MAX_ALLOCATABLE_PTHREADS];
  pthread_mutex_t mutex_LIBRARIES=PTHREAD_MUTEX_INITIALIZER;

  //#define _PTHREAD_FLUSH_TIME_ 1
#define _XPLUG_FLUSH_THREAD_SLEEP_ 1

  void *_threaded_interface_XPLUG_Directory(void *ptr) {
    bool CONCURRENT=TRUE;
    _threaded_XPLUG_params* pparams=(_threaded_XPLUG_params*) ptr;
    string directory; uint ith;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
    AFLOW_PTHREADS::RUNNING++;
    if((pparams->VERBOSE)) { pthread_mutex_lock(&mutex_LIBRARIES);cout << "_threaded_COMMANDS " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;pthread_mutex_unlock(&mutex_LIBRARIES); }
    if(CONCURRENT==FALSE) { // SPLITS tasks in MOD(THREAD_MAX);
      for(ith=(pparams->ITHREAD);ith<(*pparams->vdirs).size();ith+=(pparams->THREADS_MAX)) {
	if((pparams->VERBOSE)) { pthread_mutex_lock(&mutex_LIBRARIES);cout <<  (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << " " << ith << " " << (*pparams->vdirs).at(ith) << endl;pthread_mutex_unlock(&mutex_LIBRARIES); }
	stringstream oss;
	(*pparams->vok).at(ith)=XPLUG_Directory_ok((*pparams->vdirs).at(ith),ith,(*pparams->vdirs).size(),ith,oss); // do this
	pthread_mutex_lock(&mutex_LIBRARIES);
	cerr << oss.str();cerr.flush();
	pthread_mutex_unlock(&mutex_LIBRARIES);
      }
    } else { // RUNS in a queue
      while((*pparams->vrun).size()>0) {
	bool FRONT=TRUE;
	pthread_mutex_lock(&mutex_LIBRARIES);
	if(FRONT)  { ith=(*pparams->vrun).at(0);(*pparams->vrun).pop_front(); } // from the front
	if(!FRONT) { ith=(*pparams->vrun).at((*pparams->vrun).size()-1);(*pparams->vrun).pop_back(); } // from the back
	pthread_mutex_unlock(&mutex_LIBRARIES);
	if((pparams->VERBOSE)) { pthread_mutex_lock(&mutex_LIBRARIES);cout <<  (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << ": " << ith << endl;pthread_mutex_unlock(&mutex_LIBRARIES); }
	stringstream oss;
	(*pparams->vok).at(ith)=XPLUG_Directory_ok((*pparams->vdirs).at(ith),ith,(*pparams->vdirs).size(),ith,oss); // do this
	pthread_mutex_lock(&mutex_LIBRARIES);
	cerr << oss.str();cerr.flush();
	pthread_mutex_unlock(&mutex_LIBRARIES);
      }
    }
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
    AFLOW_PTHREADS::RUNNING--;
    return NULL;
  }
}

// DO ONE
namespace aflowlib {
  bool XPLUG_Directory_ok(const string &directory,const uint &i,const uint &j,const int& ithread,stringstream &oss) {
    string dir=directory;
    bool ok=FALSE,print=TRUE;
    int answer=0;
    stringstream obb;

    // clean
    aurostd::RemoveSubString(dir,"/"+_AFLOWIN_);
    aurostd::RemoveSubString(dir,"/aflow.end.out");

    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
    deque<string> vcmd; aurostd::string2tokens("bzcat,xzcat,zcat",vcmd,",");
    deque<string> vrelax; aurostd::string2tokens(".relax1,.relax2,.relax3,.static,.bands",vrelax,",");
    deque<string> vrem;aurostd::string2tokens("OUTCAR,OSZICAR,CONTCAR,EIGENVAL",vrem,",");

    for(uint irem=0;irem<vrem.size();irem++) {
      for(uint iext=0;iext<vext.size();iext++) {
	for(uint irelax=0;irelax<vrelax.size();irelax++) {
	  aurostd::RemoveSubString(dir,"/"+vrem.at(irem)+vrelax.at(irelax)+vext.at(iext));
	}
      }
    }

    if(!aurostd::FileExist(dir+"/"+_AFLOWLOCK_)) {
      return FALSE;
    }
    if(!aurostd::EFileExist(dir+"/OUTCAR.relax1") && 
       !aurostd::EFileExist(dir+"/OUTCAR.relax2") &&
       !aurostd::EFileExist(dir+"/OUTCAR.relax3") &&
       !aurostd::EFileExist(dir+"/OUTCAR.static") &&
       !aurostd::EFileExist(dir+"/OUTCAR.bands")) {
      return FALSE;
    }

    obb<<"[" << i+1 << "/" << j << "] ";

    if(ithread<0)  obb<<"DIR=" << dir << " ";
    if(ithread>=0) obb<<"DIR(" << ithread << ")="<< dir << " ";
    // CHECK ALL FILES
    print=TRUE;
    ok=TRUE;
    // TEST INCOMPLETE
    if(ok) { obb<<".";if(!aurostd::FileExist(dir+"/"+_AFLOWLOCK_)) { ok=FALSE;obb<<" no=LOCK"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/EIGENVAL")) { ok=FALSE;obb<<" yes=EIGENVAL"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/vasp.out")) { ok=FALSE;obb<<" yes=vasp.out"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/AECCAR0")) { ok=FALSE;obb<<" yes=AECCAR0"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/AECCAR1")) { ok=FALSE;obb<<" yes=AECCAR1"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/AECCAR2")) { ok=FALSE;obb<<" yes=AECCAR2"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/CHGCAR")) { ok=FALSE;obb<<" yes=CHGCAR"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/DOSCAR")) { ok=FALSE;obb<<" yes=DOSCAR"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/OUTCAR")) { ok=FALSE;obb<<" yes=OUTCAR"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/OSZICAR")) { ok=FALSE;obb<<" yes=OSZICAR"; }}
    if(ok) { obb<<".";if(aurostd::FileExist(dir+"/POTCAR")) { ok=FALSE;obb<<" yes=POTCAR"; }}
    for(uint irelax=0;irelax<vrelax.size();irelax++) {
      if(ok) { obb<<".";if(aurostd::FileExist(dir+"/OUTCAR"+vrelax.at(irelax))) { ok=FALSE;obb<<" yes=OUTCAR"+vrelax.at(irelax); }}
      if(ok) { obb<<".";if(aurostd::FileExist(dir+"/OSZICAR"+vrelax.at(irelax))) { ok=FALSE;obb<<" yes=OSZICAR"+vrelax.at(irelax); }}
      if(ok) { obb<<".";if(aurostd::FileExist(dir+"/POTCAR"+vrelax.at(irelax))) { ok=FALSE;obb<<" yes=POTCAR"+vrelax.at(irelax); }}
    }
    // TEST RELAX1 RELAX2 RELAX3
    for(uint irelax=0;irelax<vrelax.size();irelax++) {
      if(ok && aurostd::EFileExist(dir+"/OUTCAR"+vrelax.at(irelax))) { // relax1 relax2 relax3
	if(ok) { obb<<".";if(!aurostd::EFileExist(dir+"/OUTCAR"+vrelax.at(irelax))) { ok=FALSE;print=FALSE;obb<<" no=OUTCAR"+vrelax.at(irelax)+".EXT"; }}
	if(ok) { obb<<".";if(!aurostd::EFileExist(dir+"/OSZICAR"+vrelax.at(irelax))) { ok=FALSE;obb<<" no=OSZICAR"+vrelax.at(irelax)+".EXT"; }}
	if(ok) { obb<<".";if(!aurostd::EFileExist(dir+"/vasp.out"+vrelax.at(irelax))) { ok=FALSE;obb<<" no=vasp.out"+vrelax.at(irelax)+".EXT"; }}
	if(ok) { obb<<".";if(!aurostd::EFileExist(dir+"/EIGENVAL"+vrelax.at(irelax))) { ok=FALSE;obb<<" no=EIGENVAL"+vrelax.at(irelax)+".EXT"; }}
	if(vrelax.at(irelax)==".static")
	  if(ok) { obb<<".";if(!aurostd::EFileExist(dir+"/DOSCAR"+vrelax.at(irelax))) { ok=FALSE;obb<<" no=DOSCAR"+vrelax.at(irelax)+".EXT"; }}
	// OUTCAR"+vrelax.at(irelax)+".EXT
	if(ok) { obb<<".";  // check Answer 4 or 5 in OUTCAR.RELAX.EXT
	  for(uint iext=0;iext<vext.size();iext++) {
	    if(aurostd::FileExist(dir+"/OUTCAR"+vrelax.at(irelax)+vext.at(iext))) {
	      answer=aurostd::execute2utype<int>(vcmd.at(iext)+" "+dir+"/OUTCAR"+vrelax.at(irelax)+vext.at(iext)+" | grep -c \"(sec)\" ");
	      if(answer!=4 && answer!=5) { ok=FALSE;obb<<" error(" << answer << ")=OUTCAR"+vrelax.at(irelax)+".EXT"; }}
	  }
	}
	if(ok) { obb<<".";
	  for(uint iext=0;iext<vext.size();iext++) {
	    if(aurostd::FileExist(dir+"/vasp.out"+vrelax.at(irelax)+vext.at(iext))) {
	      answer=aurostd::execute2utype<int>(vcmd.at(iext)+" "+dir+"/vasp.out"+vrelax.at(irelax)+vext.at(iext)+" | grep -c \"The distance between some ions is very small\" ");
	      if(answer!=0) { ok=FALSE;obb<<" ions=vasp.out"+vrelax.at(irelax)+".EXT"; }}
	  }
	}
      }
    }
    // DONE
    if(ok==TRUE) obb<<" good";
    if(ok==FALSE) obb<<" bad";
    if(aurostd::FileExist(dir+"/OUTCAR.relax2.gz")) obb<<"      gz";
    if(aurostd::FileExist(dir+"/OUTCAR.relax2.bz2")) obb<<"      bz2";
    if(aurostd::FileExist(dir+"/OUTCAR.relax2.xz")) obb<<"      xz";
    obb<<endl;
    if(print==TRUE) oss << obb.str();
    oss.flush();
    return ok;
  }
}

namespace aflowlib {
  bool XPLUG(vector<string> argv) {
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    int NUM_THREADS=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"));
    int NUM_ZIP=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_ZIP"));
    int NUM_SIZE=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_SIZE"));
    bool FLAG_DO_CLEAN=XHOST.vflag_control.flag("XPLUG_DO_CLEAN");
    bool FLAG_DO_ADD=XHOST.vflag_control.flag("XPLUG_DO_ADD");
    string PREFIX=XHOST.vflag_control.getattachedscheme("XPLUG_PREFIX");

    if(NUM_THREADS<1) NUM_THREADS=XHOST.CPU_Cores/2;
    if(NUM_ZIP<1) NUM_ZIP=1;
    if(NUM_SIZE<1) NUM_SIZE=128;

    if(LDEBUG) cerr << "aflowlib::XPLUG: NUM_THREADS=" << NUM_THREADS << endl;
    if(LDEBUG) cerr << "aflowlib::XPLUG: NUM_ZIP=" << NUM_ZIP << endl;
    if(LDEBUG) cerr << "aflowlib::XPLUG: NUM_SIZE=" << NUM_SIZE << endl;
    if(LDEBUG) cerr << "aflowlib::XPLUG: FLAG_DO_CLEAN=" << FLAG_DO_CLEAN << endl;
    if(LDEBUG) cerr << "aflowlib::XPLUG: FLAG_DO_ADD=" << FLAG_DO_ADD << endl;
    if(LDEBUG) cerr << "aflowlib::XPLUG: PREFIX=" << PREFIX << endl;
    //  exit(0);
    // if(LDEBUG) cerr << "aflowlib::XPLUG = " << XHOST.hostname << endl;
    deque<string> vdirs,vzips,vcleans;
    deque<bool> vok;
    deque<uint> vrun;

    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");
    deque<string> vcmd; aurostd::string2tokens("bzcat,xzcat,zcat",vcmd,",");
    deque<string> vrelax; aurostd::string2tokens(".relax1,.relax2,.relax3,.static,.bands",vrelax,",");
    deque<string> vrem;aurostd::string2tokens("OUTCAR,OSZICAR,CONTCAR,EIGENVAL",vrem,",");

    for(uint i=0;i<argv.size();i++) {
      if(aurostd::substring2bool(argv.at(i),_AFLOWIN_))
	vdirs.push_back(aurostd::RemoveSubString(argv.at(i),"/"+_AFLOWIN_));
      if(aurostd::substring2bool(argv.at(i),"aflow.end.out"))
	vdirs.push_back(aurostd::RemoveSubString(argv.at(i),"/aflow.end.out"));
      for(uint irem=0;irem<vrem.size();irem++) {
	for(uint iext=0;iext<vext.size();iext++) {
	  for(uint irelax=0;irelax<vrelax.size();irelax++) {
	    if(aurostd::substring2bool(argv.at(i),vrem.at(irem)+vrelax.at(irelax)+vext.at(iext)))
	      vdirs.push_back(aurostd::RemoveSubString(argv.at(i),"/"+vrem.at(irem)+vrelax.at(irelax)+vext.at(iext)));
	  }
	}
      }
    }

    std::sort(vdirs.begin(),vdirs.end());
    if(LDEBUG) cerr << "aflowlib::XPLUG: vdirs.size()=" << vdirs.size() << endl;
    //  if(LDEBUG) cerr << "aflowlib::XPLUG: [1]" << endl;
    for(uint i=0;i<vdirs.size();i++) { vok.push_back(TRUE);vrun.push_back(i); }
    //  if(LDEBUG) cerr << "aflowlib::XPLUG: [2]" << endl;

    if((int) vdirs.size()<=NUM_THREADS) NUM_THREADS=(uint) vdirs.size();        // SAFETY
    //  if(LDEBUG) cerr << "aflowlib::XPLUG: [3]" << endl;

    if(NUM_THREADS<=1) {                                                        // run singular
      //  if(LDEBUG) cerr << "aflowlib::XPLUG: [3b]" << endl;
      for(uint i=0;i<vdirs.size();i++) {
	stringstream oss;
	vok.at(i)=XPLUG_Directory_ok(vdirs.at(i),i,vdirs.size(),-1,oss);
	cerr << oss.str();cerr.flush();
      }
    }
    //  if(LDEBUG) cerr << "aflowlib::XPLUG: [4]" << endl;
    if(NUM_THREADS>=2) {                                                       // multithread
      AFLOW_PTHREADS::FLAG=TRUE;AFLOW_PTHREADS::MAX_PTHREADS=NUM_THREADS;      // prepare
      if(AFLOW_PTHREADS::MAX_PTHREADS>MAX_ALLOCATABLE_PTHREADS) AFLOW_PTHREADS::MAX_PTHREADS=MAX_ALLOCATABLE_PTHREADS; // check max
      //  if(LDEBUG) cerr << "aflowlib::XPLUG: [5]" << endl;
      AFLOW_PTHREADS::Clean_Threads();                                         // clean threads
      //    if(LDEBUG) cerr << "aflowlib::XPLUG: [6]" << endl;
      _threaded_XPLUG_params params[MAX_ALLOCATABLE_PTHREADS];                 // prepare
      //    if(LDEBUG) cerr << "aflowlib::XPLUG: [7]" << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {      // prepare loop
	params[ithread].ITHREAD=ithread;                                       // prepare params
	params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;              // prepare params
	params[ithread].vok=&vok;                                              // prepare params
	params[ithread].vdirs=&vdirs;                                          // prepare params
	params[ithread].vrun=&vrun;                                            // prepare params
	params[ithread].itbusy=ithread;                                        // prepare params
	params[ithread].VERBOSE=FALSE;                                         // prepare params
      }
      //    if(LDEBUG) cerr << "aflowlib::XPLUG: [8]" << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) AFLOW_PTHREADS::viret[ithread]=pthread_create(&AFLOW_PTHREADS::vpthread[ithread],NULL,_threaded_interface_XPLUG_Directory,(void*)&params[ithread]); // run threads
      //    if(LDEBUG) cerr << "aflowlib::XPLUG: [9]" << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) pthread_join(AFLOW_PTHREADS::vpthread[ithread],NULL); // flush threads
      //  if(LDEBUG) cerr << "aflowlib::XPLUG: [10]" << endl;
    }

    // has OK. now do the counting
    for(uint i=0;i<vdirs.size();i++) {
      if(vok.at(i)==TRUE) vzips.push_back(vdirs.at(i));
      if(vok.at(i)==FALSE) vcleans.push_back(vdirs.at(i));
    }
    if(LDEBUG) cerr << "aflowlib::XPLUG: vdirs.size()=" << vdirs.size() << endl;
    if(LDEBUG) cerr << "aflowlib::XPLUG: vzips.size()=" << vzips.size() << endl;
    if(LDEBUG) cerr << "aflowlib::XPLUG: vcleans.size()=" << vcleans.size() << endl;


    if(vzips.size()>0) {
      stringstream command;
      if(aurostd::substring2bool(XHOST.hostname,"m7int0") || aurostd::substring2bool(XHOST.hostname,"m6int0")) XHOST.hostname="marylou";
      aurostd::RemoveSubString(XHOST.hostname,".egr.duke.edu");
      aurostd::RemoveSubString(XHOST.hostname,".mems.duke.edu");
      aurostd::RemoveSubString(XHOST.hostname,".pratt.duke.edu");
      aurostd::RemoveSubString(XHOST.hostname,".duke.edu");
      for(uint i=0;i<vzips.size();i+=AFLOW_MAX_ARGV) {
	command << "aflow --multi=zip " << (FLAG_DO_ADD?"--add ":"") << "--np=" << NUM_ZIP
		<< " --prefix=update_" << (PREFIX!=""?string(PREFIX+"_"):string("")) << aurostd::get_date()
		<< "-" << aurostd::get_hour() << aurostd::get_min() << aurostd::get_sec() << "_" << i/AFLOW_MAX_ARGV+1 << "_" << XHOST.hostname
		<< " --size=" << NUM_SIZE << " --DIRECTORY ";
	for(uint j=i;j<i+AFLOW_MAX_ARGV && j<vzips.size();j++)
	  command << " " << vzips.at(j);
	command << " " << endl << endl;
      }
      cerr << command.str() << endl;
      aurostd::execute(command);
    }

    /*
      if(vzips.size()>0) {
      stringstream command;
      if(aurostd::substring2bool(XHOST.hostname,"m7int0") || aurostd::substring2bool(XHOST.hostname,"m6int0")) XHOST.hostname="marylou";
      aurostd::RemoveSubString(XHOST.hostname,".egr.duke.edu");
      aurostd::RemoveSubString(XHOST.hostname,".mems.duke.edu");
      aurostd::RemoveSubString(XHOST.hostname,".pratt.duke.edu");
      aurostd::RemoveSubString(XHOST.hostname,".duke.edu");
      command << "aflow --multi=zip " << (FLAG_DO_ADD?"--add ":"") << "--np=" << NUM_ZIP
      << " --prefix=update_" << (PREFIX!=""?string(PREFIX+"_"):string("")) << aurostd::get_date()
      << "-" << aurostd::get_hour() << aurostd::get_min() << aurostd::get_sec() << "_" << XHOST.hostname
      << " --size=" << NUM_SIZE << " --DIRECTORY ";
      for(uint j=0;j<vzips.size()&&j<AFLOW_MAX_ARGV;j++)
      command << " " << vzips.at(j);
      cerr << command.str() << endl;
      aurostd::execute(command);
      }
    */
    if(FLAG_DO_CLEAN && vcleans.size()>0) {
      for(uint i=0;i<vcleans.size();i++) {
	// cerr << "Cleaning=" << vcleans.at(i) << endl;
	KBIN::Clean(vcleans.at(i));
      }
    }
    return FALSE;
  }
}

// ***************************************************************************
// ***************************************************************************

bool EnthalpyReference(string pp_version,deque<double> vLDAU,string& gs,double& enthalpy_atom,double& volume_atom,double& spin_atom, ostream& oss) {
  gs="",enthalpy_atom=999999,volume_atom=999999,spin_atom=999999;
  string ppp=pp_version;
  // oss << ppp << endl;
  // oss << type << endl;
  // oss << vLDAU.size() << endl;
  vector<string> tokens;aurostd::string2tokens(ppp,tokens,":");
  if(tokens.size()==0) { cerr << "ERROR (EnthalpyReference), [E0]  pp_version=" << pp_version << endl; exit(0); }
  if(tokens.size()==1) { cerr << "ERROR (EnthalpyReference), [E1]  pp_version=" << pp_version << endl; exit(0); }
  if(tokens.size()==2) { cerr << "ERROR (EnthalpyReference), [E2]  pp_version=" << pp_version << endl; exit(0); }
  string type=tokens.at(1);
  string date=tokens.at(2);
  bool flag_LDAU=FALSE;
  // OVERRIDE  if(vLDAU.size()>0) flag_LDAU=TRUE;

  //  if(type=="US") { oss << "ERROR (EnthalpyReference): US_LDA/GGA not supported yet. pp_version=" << ppp << endl;return FALSE; }  // LDA and GGA
  if(type=="LDA") {
    if(!flag_LDAU) { //  NO LDAU
      if(ppp=="Zr:LDA:01Apr2000") { gs="A1";enthalpy_atom=-9.33851;volume_atom=21.116;spin_atom=0.0;return TRUE; }
    }
  }
  if(type=="GGA") {
    if(!flag_LDAU) { //  NO LDAU
      if(ppp=="Pt:GGA:01Apr2000") { gs="A1";enthalpy_atom=-6.01482;volume_atom=15.849;spin_atom=0.0;return TRUE; }
    }
  }
  if(type=="PAW" || type=="PAW_LDA") {
    if(!flag_LDAU) { // PAW_LDA NO LDAU
      if(ppp=="Ag:PAW_LDA:17Apr2000") { gs="A1";enthalpy_atom=-3.74531;volume_atom=16.0612;spin_atom=0.0;return TRUE; }
      if(ppp=="Al:PAW_LDA:17Apr2000") { gs="A1";enthalpy_atom=-4.19212;volume_atom=15.8018;spin_atom=0.0;return TRUE; }
      if(ppp=="La_s:PAW_LDA:17Apr2000") { gs="A1";enthalpy_atom=-5.5382;volume_atom=32.3194;spin_atom=0.0;return TRUE; }
    }
  }
  if(type=="PAW_GGA") {
    if(!flag_LDAU) { // GGA_PBE NO LDAU
      if(ppp=="Ac_s:PAW_GGA:11Apr2000") { gs="A1";enthalpy_atom=-4.04579;volume_atom=44.9982;spin_atom=0.0;return TRUE; }
      if(ppp=="Ag:PAW_GGA:18Jul2000") { gs="A1";enthalpy_atom=-2.72686;volume_atom=17.7851;spin_atom=0.0;return TRUE; }
      if(ppp=="Al:PAW_GGA:05Jan2001") { gs="A1";enthalpy_atom=-3.6961;volume_atom=16.5601;spin_atom=0.0;return TRUE; }
      if(ppp=="B:PAW_GGA:18Jul2000") { gs="ICSD_240995";enthalpy_atom=0;volume_atom=7.76493;spin_atom=0.0;return TRUE; }  // FIX
      if(ppp=="B_s:PAW_GGA:21Dec2000") { gs="ICSD_240995";enthalpy_atom=-6.52874;volume_atom=7.79803;spin_atom=0.0;return TRUE; } // FIX
      if(ppp=="B_h:PAW_GGA:18Jul2000") { gs="ICSD_108026";enthalpy_atom=-6.177995603333333;volume_atom=6.717450;spin_atom=0.0;return TRUE; }  // FIX
      if(ppp=="Ge_h:PAW_RPBE:09Apr2002") { gs="A4";enthalpy_atom=-4.53028;volume_atom=23.6635;spin_atom=0;return TRUE; }
      if(ppp=="Ir:PAW_GGA:04May1998") { gs="A1";enthalpy_atom=-8.79333;volume_atom=14.5807;spin_atom=0.0;return TRUE; }
      if(ppp=="La_s:PAW_GGA:17Apr2000") { gs="A1";enthalpy_atom=-4.86344;volume_atom=36.6893;spin_atom=0.0;return TRUE; }
      if(ppp=="Li_sv:PAW_GGA:23Jan2001") { gs="A2*";enthalpy_atom=-1.89772;volume_atom=20.3318;spin_atom=0.0;return TRUE; }
      if(ppp=="Mg_pv:PAW_GGA:10Feb1998") { gs="A3";enthalpy_atom=-1.4779;volume_atom=22.8092;spin_atom=0.0;return TRUE; }
      if(ppp=="Na_sv:PAW_GGA:28Sep2000") { gs="A7*";enthalpy_atom=-1.31443;volume_atom=35.0106;spin_atom=0.0;return TRUE; }
      if(ppp=="Rh_pv:PAW_GGA:17Apr2000") { gs="A1";enthalpy_atom=-7.13614;volume_atom=14.1992;spin_atom=0.0;return TRUE; }
      if(ppp=="Si:PAW_GGA:05Jan2001") { gs="A4";enthalpy_atom=-5.43028;volume_atom=20.4209;spin_atom=0.0;return TRUE; }
      if(ppp=="Sm_2:PAW_GGA:07Jan2002") { gs="A2";enthalpy_atom=0.0;volume_atom=0.0;spin_atom=0.0;return TRUE; } // FIX
      if(ppp=="Sm_3:PAW_GGA:11May2000") { gs="ICSD_246657";enthalpy_atom=-4.621400;volume_atom=33.447633;spin_atom=0.0;return TRUE; } // FIX

      if(ppp=="Sm_3:PAW_GGA:11May2000" && 0) { gs="ICSD_652637";enthalpy_atom=-4.64136;volume_atom=33.5075;spin_atom=0.0;return TRUE; } // IT HAS LDAU
    }
  }
  if(type=="PAW_PBE" || type=="PAW_RPBE") {
    if(!flag_LDAU) { // PAW_PBE NO LDAU
      if(ppp=="Ac:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-4.09443;volume_atom=45.4098;spin_atom=0.0;return TRUE; }
      if(ppp=="Ac_s:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-4.02653;volume_atom=44.7612;spin_atom=0.0;return TRUE; }
      if(ppp=="Ag:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-2.82769;volume_atom=17.9065;spin_atom=0.0;return TRUE; }
      if(ppp=="Al:PAW_PBE:04Jan2001") { gs="A1";enthalpy_atom=-3.74357;volume_atom=16.4694;spin_atom=0.0;return TRUE; }
      if(ppp=="Al_h:PAW_PBE:08Apr2002") { gs="A1";enthalpy_atom=-3.80008;volume_atom=14.8772;spin_atom=0.0;return TRUE; }
      if(ppp=="As:PAW_PBE:06Sep2000") { gs="A7";enthalpy_atom=-4.65233;volume_atom=22.6706;spin_atom=0.0;return TRUE; }
      if(ppp=="Au:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-3.27222;volume_atom=18.0933;spin_atom=0.0;return TRUE; }
      // OLD if(ppp=="B:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-5.96676;volume_atom=9.11282;spin_atom=0.0;return TRUE; }
      // OLD if(ppp=="B_s:PAW_PBE:22Jan2003") { gs="A3";enthalpy_atom=-5.9736;volume_atom=8.8771;spin_atom=0.0;return TRUE; }
      // OLD  if(ppp=="B_h:PAW_PBE:07Sep2000") { gs="A3";enthalpy_atom=-5.98337;volume_atom=9.39294;spin_atom=0.0;return TRUE; }
      if(ppp=="B:PAW_PBE:06Sep2000") { gs="ICSD_94429";enthalpy_atom=-6.6775;volume_atom=7.25055;spin_atom=0.0;return TRUE; }  // FIX
      if(ppp=="B_s:PAW_PBE:22Jan2003") { gs="ICSD_240995";enthalpy_atom=-6.52874;volume_atom=7.79803;spin_atom=0.0;return TRUE; }  // FIX
      if(ppp=="B_h:PAW_PBE:07Sep2000") { gs="";enthalpy_atom=0;volume_atom=0;spin_atom=0.0;return TRUE; }  // FIX
      if(ppp=="Ba_sv:PAW_PBE:06Sep2000") { gs="A7";enthalpy_atom=-1.92399;volume_atom=63.1684;spin_atom=0.0;return TRUE; }
      if(ppp=="Be:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-3.75362;volume_atom=7.92077;spin_atom=0.0;return TRUE; }
      if(ppp=="Be_sv:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-3.7407;volume_atom=7.92019;spin_atom=0.0;return TRUE; }
      if(ppp=="Bi:PAW_PBE:08Apr2002") { gs="A7";enthalpy_atom=-3.87275;volume_atom=36.8656;spin_atom=0.0;return TRUE; }
      if(ppp=="Bi_d:PAW_PBE:06Sep2000") { gs="A7";enthalpy_atom=-4.03743;volume_atom=36.2786;spin_atom=0.0;return TRUE; }
      if(ppp=="Br:PAW_PBE:06Sep2000") { gs="A11";enthalpy_atom=-1.5898;volume_atom=45.5907;spin_atom=0.0;return TRUE; }
      if(ppp=="C:PAW_PBE:08Apr2002") { gs="A9";enthalpy_atom=-9.22165;volume_atom=10.4453;spin_atom=0.0;return TRUE; }
      // if(ppp="C_s:PAW_PBE:06Sep2000") { gs="A9";enthalpy_atom=;volume_atom=;spin_atom=0.0;return TRUE; }
      // if(ppp="C_h:PAW_PBE:20Dec2001") { gs="A9";enthalpy_atom=;volume_atom=;spin_atom=0.0;return TRUE; }
      if(ppp=="Ca_pv:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-1.97683;volume_atom=41.7973;spin_atom=0.0;return TRUE; }
      if(ppp=="Ca_sv:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-2.00151;volume_atom=42.2101;spin_atom=0.0;return TRUE; }
      if(ppp=="Cd:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-0.905693;volume_atom=22.4532;spin_atom=0.0;return TRUE; }
      // if(ppp=="Ce") { gs="ICSD_2284-mS4";enthalpy_atom=-5.93013;volume_atom=26.0697;spin_atom=0.0;return TRUE; }
      if(ppp=="Ce") { gs="A1";enthalpy_atom=-5.92998;volume_atom=26.0579;spin_atom=0.0;return TRUE; }
      // if(ppp=="Cl_h:PAW_PBE:08Apr2002") { gs="A11";enthalpy_atom=-1.8156;volume_atom=37.3299;spin_atom=0.0;return TRUE; } WAITING
      if(ppp=="Cl:PAW_PBE:17Jan2003") { gs="A11";enthalpy_atom=-1.8156;volume_atom=37.3299;spin_atom=0.0;return TRUE; }
      if(ppp=="Co:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-7.10872;volume_atom=10.8414;spin_atom=1.59707;return TRUE; }
      if(ppp=="Cr:PAW_PBE:06Sep2000") { gs="A2";enthalpy_atom=-9.51277;volume_atom=11.3963;spin_atom=0.0;return TRUE; }
      if(ppp=="Cr_pv:PAW_PBE:07Sep2000") { gs="A2";enthalpy_atom=-9.62922;volume_atom=11.5175;spin_atom=0.0;return TRUE; }
      if(ppp=="Cu:PAW_PBE:05Jan2001") { gs="A1";enthalpy_atom=-3.71935;volume_atom=11.9654;spin_atom=0.0;return TRUE; }
      if(ppp=="Cu_pv:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-4.09817;volume_atom=11.8644;spin_atom=0.0;return TRUE; }
      if(ppp=="Dy_3:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-4.58786;volume_atom=31.857;spin_atom=0.0;return TRUE; }
      if(ppp=="Fe:PAW_PBE:06Sep2000") { gs="A2";enthalpy_atom=-8.3111;volume_atom=11.332;spin_atom=2.19866;return TRUE; }
      if(ppp=="Fe_pv:PAW_PBE:06Sep2000") { gs="A2";enthalpy_atom=-8.45475;volume_atom=11.3277;spin_atom=2.19889;return TRUE; }
      if(ppp=="Ga_h:PAW_PBE:09Apr2002") { gs="A11";enthalpy_atom=-2.88121;volume_atom=18.6097;spin_atom=0.0;return TRUE; }
      if(ppp=="Ge_h:PAW_RPBE:09Apr2002") { gs="A4";enthalpy_atom=-4.53028;volume_atom=23.6635;spin_atom=0;return TRUE; } // PAW_GGA
      if(ppp=="Ge:PAW_PBE:05Jan2001") { gs="A4";enthalpy_atom=-4.2768;volume_atom=19.4258;spin_atom=0.0;return TRUE; }
      if(ppp=="Ge_h:PAW_PBE:09Apr2002") { gs="A4";enthalpy_atom=-4.50403;volume_atom=23.6965;spin_atom=0.0;return TRUE; }
      if(ppp=="Ge_d:PAW_PBE:06Sep2000") { gs="A4";enthalpy_atom=-4.62207;volume_atom=23.7871;spin_atom=0.0;return TRUE; }
      if(ppp=="Ge:PAW_PBE:05Jan2001") { gs="A4";enthalpy_atom=-4.49258;volume_atom=24.1781;spin_atom=0.0;return TRUE; }
      if(ppp=="Hf_pv:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-9.9527;volume_atom=22.4349;spin_atom=0.0;return TRUE; }
      if(ppp=="Hg:PAW_PBE:06Sep2000") { gs="A10(A3)";enthalpy_atom=-0.299947;volume_atom=29.356;spin_atom=0.0;return TRUE; }
      if(ppp=="Ho_3:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-4.5683;volume_atom=31.3589;spin_atom=0.0;return TRUE; }
      if(ppp=="In_d:PAW_PBE:06Sep2000") { gs="A6";enthalpy_atom=-2.72095;volume_atom=27.1064;spin_atom=0.0;return TRUE; }
      if(ppp=="Ir:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-8.85703;volume_atom=14.526;spin_atom=0.0;return TRUE; }
      if(ppp=="K_pv:PAW_PBE:17Jan2003") { gs="A2";enthalpy_atom=-1.026832;volume_atom=72.4893;spin_atom=0.0;return TRUE; }
      if(ppp=="K_sv:PAW_PBE:06Sep2000") { gs="A2";enthalpy_atom=-1.09689;volume_atom=73.4415;spin_atom=0.0;return TRUE; }
      if(ppp=="La:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-4.91809;volume_atom=36.5244;spin_atom=0.0;return TRUE; }
      if(ppp=="Li:PAW_PBE:17Jan2003") { gs="A2*";enthalpy_atom=-1.897466;volume_atom=20.3189;spin_atom=0.0;return TRUE; }
      if(ppp=="Li_sv:PAW_PBE:23Jan2001") { gs="A2*";enthalpy_atom=-1.90503;volume_atom=20.2771;spin_atom=0.0;return TRUE; }
      if(ppp=="Mg:PAW_PBE:05Jan2001") { gs="A3";enthalpy_atom=-1.54144;volume_atom=22.8463;spin_atom=0.0;return TRUE; }
      if(ppp=="Mg_pv:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-1.594;volume_atom=22.7983;spin_atom=0.0;return TRUE; }
      if(ppp=="Mn:PAW_PBE:06Sep2000") { gs="A12";enthalpy_atom=-9.02786;volume_atom=10.7275;spin_atom=0.0;return TRUE; }
      if(ppp=="Mn_pv:PAW_PBE:07Sep2000") { gs="A12";enthalpy_atom=-9.1535;volume_atom=10.6973;spin_atom=0.0;return TRUE; }
      if(ppp=="Mo:PAW_PBE:08Apr2002") { gs="A2";enthalpy_atom=-10.945857;volume_atom=15.5844;spin_atom=0.0;return TRUE; }
      if(ppp=="Mo_pv:PAW_PBE:08Apr2002") { gs="A2";enthalpy_atom=-10.8434;volume_atom=15.8555;spin_atom=0.0;return TRUE; }
      if(ppp=="Na_pv:PAW_PBE:05Jan2001") { gs="A2";enthalpy_atom=-1.31096;volume_atom=36.2541;spin_atom=0.0;return TRUE; }
      if(ppp=="Na_sv:PAW_PBE:28Sep2000") { gs="A7*";enthalpy_atom=-1.31537;volume_atom=35.606;spin_atom=0.0;return TRUE; }
      if(ppp=="Nb_sv:PAW_PBE:17Jan2003") { gs="A2";enthalpy_atom=-10.2246;volume_atom=18.0525;spin_atom=0.0;return TRUE; }
      if(ppp=="Ni_pv:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-5.7786;volume_atom=10.8136;spin_atom=0.628888;return TRUE; }
      if(ppp=="O:PAW_PBE:08Apr2002") { gs="A8";enthalpy_atom=-4.50622;volume_atom=12.5315;spin_atom=0.0;return TRUE; }
      if(ppp=="Os_pv:PAW_PBE:20Jan2003") { gs="A3";enthalpy_atom=-11.2191;volume_atom=14.2979;spin_atom=0.0;return TRUE; }
      if(ppp=="P:PAW_PBE:17Jan2003") { gs="A7";enthalpy_atom=-5.32414;volume_atom=16.0927;spin_atom=0.0;return TRUE; }
      if(ppp=="Pb_d:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-3.70507;volume_atom=31.6585;spin_atom=0.0;return TRUE; }
      if(ppp=="Pd_pv:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-5.38277;volume_atom=15.3962;spin_atom=0.0;return TRUE; }
      if(ppp=="Pt:PAW_PBE:05Jan2001") { gs="A1";enthalpy_atom=-6.05482;volume_atom=15.6527;spin_atom=0.0;return TRUE; }
      if(ppp=="Rb_sv:PAW_PBE:06Sep2000") { gs="A2* ";enthalpy_atom=-0.962354;volume_atom=90.2922;spin_atom=0.0;return TRUE; }
      if(ppp=="Re_pv:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-12.4322;volume_atom=14.9985;spin_atom=0.0;return TRUE; }
      if(ppp=="Rh_pv:PAW_PBE:06Sep2000") { gs="A1";enthalpy_atom=-7.34057;volume_atom=14.1854;spin_atom=0.0;return TRUE; }
      if(ppp=="Ru_pv:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-9.27196;volume_atom=13.8821;spin_atom=0.0;return TRUE; }
      if(ppp=="Sb:PAW_PBE:06Sep2000") { gs="A7";enthalpy_atom=-3.88932;volume_atom=27.1685;spin_atom=0.0;return TRUE; }
      if(ppp=="Sc_sv:PAW_PBE:07Sep2000") { gs="A3";enthalpy_atom=-6.33209;volume_atom=24.4214;spin_atom=0.0;return TRUE; }
      if(ppp=="Se:PAW_PBE:06Sep2000") { gs="A8";enthalpy_atom=-3.48283;volume_atom=29.6441;spin_atom=0.0;return TRUE; }
      if(ppp=="Si:PAW_PBE:05Jan2001") { gs="A4";enthalpy_atom=-5.42373;volume_atom=20.4311;spin_atom=0.0;return TRUE; }
      if(ppp=="Si_h:PAW_PBE:08Apr2002") { gs="A4";enthalpy_atom=-5.44183;volume_atom=20.3614;spin_atom=0.0;return TRUE; }
      if(ppp=="Sm_3:PAW_PBE:07Sep2000") { gs="FIX";enthalpy_atom=0;volume_atom=0;spin_atom=0.0;return TRUE; } // FIX
      if(ppp=="Sm_3:PAW_PBE:07Sep2000" && 0) { gs="A1";enthalpy_atom=-4.7062;volume_atom=33.8339;spin_atom=0.0;return TRUE; } // IT HAS LDAU
      if(ppp=="Sm:PAW_PBE:08Apr2002") { gs="FIX";enthalpy_atom=0;volume_atom=0;spin_atom=0.0;return TRUE; } // FIX
      if(ppp=="Sn:PAW_PBE:08Apr2002") { gs="A5";enthalpy_atom=-3.79537;volume_atom=28.341;spin_atom=0;return TRUE; }
      if(ppp=="Sn_d:PAW_PBE:06Sep2000") { gs="A5";enthalpy_atom=-3.96419;volume_atom=28.1306;spin_atom=0;return TRUE; }
      if(ppp=="Sr_sv:PAW_PBE:07Sep2000") { gs="A1";enthalpy_atom=-1.68357;volume_atom=54.6082;spin_atom=0.0358379;return TRUE; }
      if(ppp=="Ta_pv:PAW_PBE:07Sep2000") { gs="A2";enthalpy_atom=-11.8492;volume_atom=18.2778;spin_atom=0.0;return TRUE; }
      if(ppp=="Tc_pv:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-10.3594;volume_atom=14.5515;spin_atom=0.0;return TRUE; }
      if(ppp=="Te:PAW_PBE:08Apr2002") { gs="A8";enthalpy_atom=-3.14134;volume_atom=34.8509;spin_atom=0.0;return TRUE; }
      if(ppp=="Ti_sv:PAW_PBE:07Sep2000") { gs="A3";enthalpy_atom=-7.93818;volume_atom=17.2475;spin_atom=0.0;return TRUE; }
      if(ppp=="Tl_d:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-2.36244;volume_atom=30.7601;spin_atom=0.0;return TRUE; }
      if(ppp=="V_sv:PAW_PBE:07Sep2000") { gs="A2";enthalpy_atom=-9.1149;volume_atom=13.3971;spin_atom=0.0;return TRUE; }
      if(ppp=="W_pv:PAW_PBE:06Sep2000") { gs="A2";enthalpy_atom=-12.9534;volume_atom=16.1649;spin_atom=0.0;return TRUE; }
      if(ppp=="Y_sv:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-6.46301;volume_atom=32.7578;spin_atom=0.0;return TRUE; }
      if(ppp=="Zn:PAW_PBE:06Sep2000") { gs="A3";enthalpy_atom=-1.26559;volume_atom=15.1693;spin_atom=0.0;return TRUE; }
      if(ppp=="Zr_sv:PAW_PBE:07Sep2000") { gs="A3";enthalpy_atom=-8.54331;volume_atom=23.4268;spin_atom=0.0;return TRUE; }
    }
  }
  if(type=="PAW_PBE_KIN") {
    if(!flag_LDAU) { // PAW_PBE_KIN NO LDAU
      if(ppp=="Ac:PAW_PBE_KIN:06Sep2000") { gs="A1";enthalpy_atom=-4.09443;volume_atom=45.4098;spin_atom=0.0;return TRUE; } // GOTTA FIX IT
      if(ppp=="Ac_s:PAW_PBE_KIN:06Sep2000") { gs="A1";enthalpy_atom=-4.02653;volume_atom=44.7612;spin_atom=0.0;return TRUE; } // GOTTA FIX IT
    }
  }
  if(type=="PAW_LDA_KIN") {
    if(!flag_LDAU) { // PAW_LDA_KIN NO LDAU
      if(ppp=="Ac:PAW_LDA_KIN:06Sep2000") { gs="A1";enthalpy_atom=-4.09443;volume_atom=45.4098;spin_atom=0.0;return TRUE; } // GOTTA FIX IT
      if(ppp=="Ac_s:PAW_LDA_KIN:06Sep2000") { gs="A1";enthalpy_atom=-4.02653;volume_atom=44.7612;spin_atom=0.0;return TRUE; } // GOTTA FIX IT
    }
  }
  oss << "ERROR (EnthalpyReference):  pp_version=" << pp_version  << " vLDAU.size()=" << vLDAU.size() << endl;// exit(0);
  return FALSE;
};

double EnthalpyReference(string pp_version,deque<double> vLDAU, ostream& oss) {
  string gs="";
  double enthalpy_atom=0.0,volume_atom=0.0,spin_atom=0.0;
  bool found=EnthalpyReference(pp_version,vLDAU,gs,enthalpy_atom,volume_atom,spin_atom,oss);
  if(found) return enthalpy_atom;
  oss << "ERROR (EnthalpyReference):  pp_version=" << pp_version  << " vLDAU.size()=" << vLDAU.size() << endl;// exit(0);
  return 0.0;
}

bool EnthalpyReferenceAvailable(string pp_version,deque<double> vLDAU) {
  string gs="";
  double enthalpy_atom=0.0,volume_atom=0.0,spin_atom=0.0;
  bool found=EnthalpyReference(pp_version,vLDAU,gs,enthalpy_atom,volume_atom,spin_atom,cout);
  if(found) return TRUE;
  return FALSE;
}


// ***************************************************************************
// void aflowlib::vaspfile2stringstream(string& str_dir, const string& FILE, stringstream& ss_vaspfile)
// ***************************************************************************
namespace aflowlib {
  bool VaspFileExist(const string& str_dir, const string& FILE) {
    bool RUN_FLAG=FALSE;
    if(!RUN_FLAG && aurostd::FileExist(str_dir+"/"+FILE)) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::EFileExist(str_dir+"/"+FILE)) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::FileExist(str_dir+"/"+FILE+".static")) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::EFileExist(str_dir+"/"+FILE+".static")) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::FileExist(str_dir+"/"+FILE+".bands")) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::EFileExist(str_dir+"/"+FILE+".bands")) RUN_FLAG=TRUE;
    if(!RUN_FLAG) {
      cerr<< FILE+" or "+FILE+".static/bands or "+FILE+".static/bands.EXT not found in the directory!"<<endl;
    }
    return RUN_FLAG;
  }
}

// ***************************************************************************
// void aflowlib::vaspfile2stringstream(string& str_dir, const string& FILE, stringstream& ss_vaspfile)
// ***************************************************************************
namespace aflowlib {
  uint _OLD_vaspfile2stringstream(const string& str_dir, const string& FILE, stringstream& sss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
    vector<string> vcat; aurostd::string2tokens("cat,bzcat,xzcat,gzcat",vcat,",");

    if(LDEBUG) cerr << "vaspfile2stringstream: BEGIN" << endl;
    //if XXXX.bands exists, you may also use this function by setting FILE=vaspfile.bands
    sss.clear(); sss.str(std::string(""));
    bool gfound=FALSE;
    if(!gfound && (FILE=="EIGENVAL" || FILE=="KPOINTS")) {
      gfound=TRUE;
      if(LDEBUG) cerr << "vaspfile2stringstream: 1=" << FILE << endl;
      for(uint iext=0;iext<vext.size();iext++) {
	if(aurostd::FileExist(str_dir+"/"+FILE+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".bands"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".bands"+vext.at(iext));return sss.str().length(); }
      }
      cerr<< FILE+" or "+FILE+".bands or "+FILE+".bands.EXT not found in the directory, aborting!"<<endl;
      exit(1);
    }
    if(!gfound && (FILE=="DOSCAR")) {
      gfound=TRUE;
      if(LDEBUG) cerr << "vaspfile2stringstream: 2=" << FILE << endl;
      for(uint iext=0;iext<vext.size();iext++) {
	if(aurostd::FileExist(str_dir+"/"+FILE+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".static"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".static"+vext.at(iext));return sss.str().length(); }
      }
      cerr<< FILE+" or "+FILE+".static or "+FILE+".static.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(!gfound && (FILE=="POSCAR")) {
      gfound=TRUE;
      if(LDEBUG) cerr << "vaspfile2stringstream: 1=" << FILE << endl;
      for(uint iext=0;iext<vext.size();iext++) {
	if(aurostd::FileExist(str_dir+"/"+FILE+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".bands"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".bands"+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".static"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".static"+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".relax1"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".relax1"+vext.at(iext));return sss.str().length(); }
      }
      cerr<< FILE+" or "+FILE+".bands/static/relax or "+FILE+".bands./static/relax.EXT not found in the directory, aborting!"<<endl;
      exit(1);
    }
    if(!gfound) {
      if(LDEBUG) cerr << "vaspfile2stringstream: 3=" << FILE << endl;
      for(uint iext=0;iext<vext.size();iext++) {
	if(aurostd::FileExist(str_dir+"/"+FILE+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".static"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".static"+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".relax1"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".relax1"+vext.at(iext));return sss.str().length(); }
	if(aurostd::FileExist(str_dir+"/"+FILE+".bands"+vext.at(iext))) { sss << aurostd::execute2string(vcat.at(iext)+" "+ str_dir+"/"+FILE+".bands"+vext.at(iext));return sss.str().length(); }
      }
      cerr<< FILE+" or "+FILE+".static/relax1/bands or "+FILE+".static/relax1/bands.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(LDEBUG) cerr << "vaspfile2stringstream: END" << endl;
    return sss.str().length();
  }
}


// ***************************************************************************
// void aflowlib::vaspfile2stringstream
// ***************************************************************************
namespace aflowlib {
  string vaspfile2stringstream(const string& str_dir, const string& FILE, stringstream& sss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: BEGIN" << endl;
    //if XXXX.bands exists, you may also use this function by setting FILE=vaspfile.bands
    sss.clear(); sss.str(std::string(""));
    bool gfound=FALSE;
    if(!gfound && (FILE=="EIGENVAL" || FILE=="KPOINTS")) {
      gfound=TRUE;
      if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.bands,.bands.bz2,.bands.gz,.bands.xz";
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front(""); // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las "+str_dir);
      for(uint i=0;i<vtype.size();i++) {
	if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
	  aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
	  if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FOUND FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	  return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
	}
      }
      cerr<< FILE+" or "+FILE+".bands or "+FILE+".bands.EXT not found in the directory, aborting!"<<endl;
      exit(1);
    }
    if(!gfound && (FILE=="DOSCAR")) {
      gfound=TRUE;
      if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.static,.static.bz2,.static.gz,.static.xz";
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front("");  // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las "+str_dir);
      for(uint i=0;i<vtype.size();i++) {
	if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
	  aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
	  if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FOUND FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	  return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
	}
      }
      cerr<< FILE+" or "+FILE+".static or "+FILE+".static.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(!gfound && (FILE=="POSCAR")) {
      gfound=TRUE;
      if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.bands,.bands.bz2,.bands.gz,.bands.xz,.static,.static.bz2,.static.gz,.static.xz,.relax1,.relax1.bz2,.relax1.gz,.relax1.xz";
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front(""); // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las "+str_dir);
      for(uint i=0;i<vtype.size();i++) {
	if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
	  aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
	  if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FOUND FILE=" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	  return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
	}
      }
      cerr<< FILE+" or "+FILE+".bands/static/relax or "+FILE+".bands./static/relax.EXT not found in the directory, aborting!"<<endl;
      exit(1);
    }
    if(!gfound) {
      if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.static,.static.bz2,.static.gz,.static.xz,.relax1,.relax1.bz2,.relax1.gz,.relax1.xz,.bands,.bands.bz2,.bands.gz,.bands.xz";
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front(""); // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las "+str_dir);
      for(uint i=0;i<vtype.size();i++) {
	if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
	  aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
	  if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FOUND FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
	  return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
	}
      }
      cerr<< FILE+" or "+FILE+".static/relax1/bands or "+FILE+".static/relax1/bands.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(LDEBUG) cerr << "vaspfile2stringstream: END" << endl;
    return string("");
  }
}

namespace aflowlib {
  string vaspfile2stringstream(const string& str_dir, const string& FILE) {
    stringstream sss;
    return vaspfile2stringstream(str_dir,FILE,sss);
  }
}

#endif //  _AFLOWLIB_LIBRARIES_CPP_
// ***************************************************************************
