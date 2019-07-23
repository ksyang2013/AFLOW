// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_PFLOW_H_
#define _AFLOW_PFLOW_H_

#include "aflow.h"

// pflow prototypes
// aflow_pflow.cpp
#define kwarning        "[WARNING] aflow: "
#define kintro          " "

// LOGGER MODES
static const char _LOGGER_ERROR_ = 'E';
static const char _LOGGER_WARNING_ = 'W';
static const char _LOGGER_COMPLETE_ = 'C';
static const char _LOGGER_OPTION_ = 'O';
static const char _LOGGER_RAW_ = 'R';
static const char _LOGGER_MESSAGE_ = 'M';

//[MOVED to aflow.h]// LOADENTRIES DEFAULTS
//[MOVED to aflow.h]static const uint _AFLOW_LIB_MAX_ = 10;  //LIB11 does not exist yet, modify accordingly

void DEVELOP(vector<string> argv);

// aflow_pflow_main.cpp
namespace pflow {
  int main(vector<string> &argv,vector<string> &cmds);
}

namespace pflow {
  bool CheckCommands(vector<string>,const vector<string>& cmds);
  void CheckIntegritiy();
  string Intro_pflow(string x);
}
string AFLOW_PrototypesIcsdReadme(void);
namespace pflow {
  xstructure XXX(istream& input,const int& iatom);
  bool XXX(vector<string>,istream& input);
}
namespace pflow {
  xstructure ABCCAR(istream& input);
  void ACE(istream& input);
  bool AddSpinToXstructure(xstructure& a, vector<double>& vmag); // DX 9/27/17 - Magnetic symmetry
  bool AddSpinToXstructure(xstructure& a, vector<xvector<double> >& vmag_noncoll); // DX 12/5/17 - Magnetic symmetry (non-collinear)
  void AGROUP(_aflags &aflags,istream& input);
  void AGROUP2(istream& input);
  void AGROUP2m(istream& input);
  xstructure ALPHABETIC(istream& input);
  string ALPHACompound(string options);
  string ALPHASpecies(string options);
  string AFLOWIN(istream& input);
  void ANGLES(string options,istream& input);
  string ATOMSMAX(string options,istream& input);
  void BANDS(string options,istream& input);
  void BANDGAP(aurostd::xoption& vpflow,ostream& oss); // CAMILO  // CO 171006
  void BANDSTRUCTURE(_aflags &aflags);
  string BZDirectionsLATTICE(string options);
  //DX 20181102 [OBSOLETE] string BZDirectionsSTRUCTURE(istream& input);
  string BZDirectionsSTRUCTURE(istream& input, aurostd::xoption& vpflow); //DX 20181102 - add options
  void CAGES(_aflags &aflags,string options,istream& input);
  // DX and CO - START
  bool PerformFullSymmetry(xstructure& a);
  bool PerformFullSymmetry(xstructure& a,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss, string format="txt");
  bool PerformFullSymmetry(xstructure& a,double& tolerance,bool no_scan,bool force_perform,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss, string format="txt");
  void defaultKFlags4SymWrite(_kflags& kflags,bool write=true);
  void defaultKFlags4SymCalc(_kflags& kflags,bool calc=true);
  bool CalculateFullSymmetry(istream& input,aurostd::xoption& vpflow,ostream& oss);
  bool CalculateFullSymmetry(_aflags &aflags, _kflags& kflags, xstructure& a,aurostd::xoption& vpflow,bool osswrite,ostream& oss);
  bool fixEmptyAtomNames(xstructure& xstr,bool force_fix=false);  //force_fix=true if you want to override what is already in species
  // DX and CO - END
  xstructure CART(istream& input);
  xstructure CORNERS(istream& input);
  void ChangeSuffix(string options);
  string CHGDIFF(aurostd::xoption vpflow);
  bool CHGDIFF(const string& chgcar1_file,const string& chgcar2_file, const string& output_File,ostream&oss);
  // DX and CO - START
  void CHGINT(vector<string>);
  // DX and CO - END
  string CHGSUM(aurostd::xoption vpflow);
  bool CHGSUM(const string& chgcar_in1,const string& chgcar_in2,ostream& oss);
  bool CHGSUM(string& species_header,const string& chgcar_in1,const string& chgcar_in2,const string& output_file,ostream& oss);
  bool CHGSUM(const string& chgcar_in1,const string& chgcar_in2,const string& output_file,ostream& oss);
  bool CHGSUM(const vector<string>& chgcar_files,ostream& oss);
  bool CHGSUM(const vector<string>& chgcar_files,const string& output,ostream& oss);
  bool CHGSUM(string& species_header,const vector<string>& chgcar_files,ostream& oss);
  bool CHGSUM(string& species_header,const vector<string>& chgcar_files,const string& output_file,ostream& oss);
  //DX 20180806 [OBSOLETE] void CIF(istream& input);
  void CIF(istream& input,aurostd::xoption& vpflow);
  void CLAT(string options);
  void CLEAN(vector<string>);
  void CLEANALL(istream& input);
  void CMPSTR(vector<string>);
  void COMPARE(string options);
  string compareStructures(aurostd::xoption& vpflow); //DAVID
  string compareStructureDirectory(aurostd::xoption& vpflow); //DAVID
  // DX 9/1/17 [OBSOLETE] void DATA(string smode,istream& input);
  bool DATA(string smode, istream& input, aurostd::xoption& vpflow, ostream& oss); // DX 9/1/17 - SGDATA + JSON
  void DATA1(string options,istream& input);
  void DATA2(istream& input);
  void DEBYE(string options);
  void DISP(string options,istream& input);
  void DIST(string options,istream& input);
  void DYNADIEL(vector<string>& argv); // CAMILO
  void EDOS(vector<string>);
  void EFFMASS(vector<string>& argv, ostream& oss); // CAMILO
  void EIGCURV(string options, ostream& oss); // CAMILO
  // DX 8/18/17 [OBSOLETE] xstructure EQUIVALENT(_aflags &aflags,istream& input);
  string EQUIVALENT(_aflags &aflags,istream& input, aurostd::xoption& vpflow);
  void EWALD(string options,istream& input);
  string EXTRACT_xcar(_aflags &aflags,vector<string>,string,string);
  string EXTRACT_Symmetry(_aflags &aflags,vector<string>);
  void FGROUP(_aflags &aflags,istream& input);
  bool FIXBANDS(_aflags &aflags,string opts);
  // DX 9/26/17 [OBSOLETE] void FINDSYM(string options,uint mode,istream& input);
  void FINDSYM(aurostd::xoption& vpflow,uint mode,istream& input);
  xstructure FRAC(istream& input);
  string FROZSL_VASPSETUP(vector<string> argv,int mode);
  string FROZSL_ANALYZE(istream& input);
  string FROZSL_INPUT(void);
  string FROZSL_OUTPUT(void);
  bool GetCollinearMagneticInfo(int& num_atoms, string& magmom_info, vector<double>& vmag); // DX 9/27/17 - Magnetic symmetry
  bool GetNonCollinearMagneticInfo(int& num_atoms, string& magmom_info, vector<xvector<double> >& vmag_noncoll); // DX 12/5/17 - Magnetic symmetry non-collinear
  void GULP(istream& input);
  void HKL(string options,_aflags &aflags,istream& input);
  void HKLSearch(string options,_aflags &aflags,istream& input,const string& smode);
  string HNF(vector<string> argv,istream& input,ostream& oss);
  string HNFTOL(vector<string> argv,istream& input,ostream& oss);
  void ICSD_2WYCK(istream& input,bool SOF);
  void ICSD(vector<string> argv, istream& input);
  xstructure IDENTICAL(istream& input);
  xstructure INCELL(istream& input);
  xstructure INCOMPACT(istream& input);
  void INTPOL(string options);
  xstructure INWS(istream& input);
  void JMOL(string options,istream& input);
  void KBAND(vector<string>);
  xstructure INFLATE_LATTICE(string options,istream& input);
  xstructure INFLATE_VOLUME(string options,istream& input);
  void KPATH(istream& input,double grid,bool WWW); 
  xstructure KPOINTS(string options,istream& input,ostream& oss);
  xstructure KPOINTS_DELTA(aurostd::xoption& vpflow, istream& input, ostream& oss);
  void JOINSTRLIST(vector<string>);
  void MAKESTRLIST(vector<string>);
  xstructure LATTICEREDUCTION(istream& input);
  string LATTICE_TYPE(istream& input);
  string LATTICE_LATTICE_TYPE(istream& input);
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //START - all relevent functions for loading entries here
  //Added by Corey Oses - May 2017
  //load entries is heavily overloaded, mostly to accommodate entries separated as
  //vector<vector<vector<> > > entries (unaries vs. binaries, then species-specific, good for convex hull),
  //vector<vector<> > entries (unaries vs. binaries), OR
  //vector<> entries (all together)
  ////////////////////////////////////////////////////////////////////////////////
  string arity_string(uint arity,bool capital=false,bool plural=false); // CO 180329
  // loadEntries
  bool loadEntries(vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadEntries
  bool loadEntries(vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadEntries
  bool loadEntries(vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadEntries(vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  bool loadFromCommon(aurostd::xoption& vpflow);
  ////////////////////////////////////////////////////////////////////////////////
  //load and merging LIBX
  bool loadAndMergeLIBX(string combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss=cout);
  bool loadAndMergeLIBX(string _combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss=cout);
  bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, string combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss=cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, string _combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ostream& oss=cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, ofstream& FileMESSAGE, ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadLIBX string elements
  bool loadLIBX(string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  // loadLIBX vector elements
  bool loadLIBX(string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadLIBX string elements, organized by -naries
  bool loadLIBX(string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  // loadLIBX_nested vector elements
  bool loadLIBX(string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ostream& oss=cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry> >& entries, ofstream& FileMESSAGE, ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  //merge vector entries lists
  bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, vector<vector<vector<aflowlib::_aflowlib_entry> > >& new_entries);
  bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, vector<vector<aflowlib::_aflowlib_entry> >& new_entries, bool assume_same_type = false);
  bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, vector<aflowlib::_aflowlib_entry>& new_entries, bool assume_same_type = false);
  bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, aflowlib::_aflowlib_entry& new_entries);
  bool mergeEntries(vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries, aflowlib::_aflowlib_entry& new_entries, int& match1, int& match2);
  bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries, vector<vector<aflowlib::_aflowlib_entry> >& new_entries, bool assume_same_type = false, bool sort_by_species = true);
  bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries, vector<aflowlib::_aflowlib_entry>& new_entries, bool assume_same_type = false, bool sort_by_species = true);
  bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries, aflowlib::_aflowlib_entry& new_entry, int& match, bool sort_by_species = true);
  bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries, aflowlib::_aflowlib_entry& new_entry, bool sort_by_species = true);
  bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries, vector<aflowlib::_aflowlib_entry>& new_entries);
  bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries, aflowlib::_aflowlib_entry& new_entry);
  //start combining
  bool mergeEntries(vector<vector<aflowlib::_aflowlib_entry> >& naries, vector<vector<vector<aflowlib::_aflowlib_entry> > >& new_entries, bool sort_by_species = true);
  bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries, vector<vector<vector<aflowlib::_aflowlib_entry> > >& new_entries);
  bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries, vector<aflowlib::_aflowlib_entry>& new_entries);
  //trivial
  bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries, vector<vector<vector<aflowlib::_aflowlib_entry> > >& new_entries);
  bool mergeEntries(vector<aflowlib::_aflowlib_entry>& naries, vector<vector<aflowlib::_aflowlib_entry> >& new_entries);
  ////////////////////////////////////////////////////////////////////////////////
  //get elemental combinations (recursively)
  //[OBSOLETE CO 180528]void getCombination(const vector<string>& velements, vector<string>& combination, vector<vector<string> >& combinations, uint offset, uint nary);
  vector<vector<string> > elementalCombinations(const vector<string>& velements, uint nary);
  ////////////////////////////////////////////////////////////////////////////////
  //easy way to think about it:  do compounds belong to the hull?
  bool compoundsBelong(const vector<string>& velements, const string& input, ostream& oss=cout, bool sort_input=false, bool clean=true);
  bool compoundsBelong(const vector<string>& velements, const string& input, ofstream& FileMESSAGE, ostream& oss=cout, bool sort_input=false, bool clean=true);
  bool compoundsBelong(const vector<string>& velements, const vector<string>& elements, ostream& oss=cout, bool sort_elements=false);
  bool compoundsBelong(const vector<string>& velements, const vector<string>& elements, ofstream& FileMESSAGE, ostream& oss=cout, bool sort_elements=false);
  ////////////////////////////////////////////////////////////////////////////////
  // loads xstructures
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, ostream& oss=cout, bool relaxed_only=true, string path="", bool is_url_path=false);
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, ofstream& FileMESSAGE, ostream& oss=cout, bool relaxed_only=true, string path="", bool is_url_path=false);
  ////////////////////////////////////////////////////////////////////////////////
  // returns UNSORTED vector<string> from string
  vector<string> stringElements2VectorElements(const string& input, ostream& oss=cout, bool clean=true);
  vector<string> stringElements2VectorElements(const string& input, ofstream& FileMESSAGE, ostream& oss=cout, bool clean=true);
  ////////////////////////////////////////////////////////////////////////////////
  // functions for making input alphabetic
  // PdMn -> MnPd, does it by CAPITAL letters
  string makeAlphabeticString(const string& input, ostream& oss=cout);
  string makeAlphabeticString(const string& input, ofstream& FileMESSAGE, ostream& oss=cout);
  vector<string> makeAlphabeticVector(const string& input, ostream& oss=cout);
  vector<string> makeAlphabeticVector(const string& input, ofstream& FileMESSAGE, ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  // sets default flags
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, ostream& oss=cout, string input="A", bool entry_output=true, bool silent=false);
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& oss=cout, string input="A", bool entry_output=true, bool silent=false);
  ////////////////////////////////////////////////////////////////////////////////
  bool prototypeMatch(string proto_database, string proto_search); //smarter than == for prototype matches, deals with 549 vs. 549.bis
  ////////////////////////////////////////////////////////////////////////////////
  //Added by Corey Oses - May 2017
  //END - all relevent functions for loading entries here
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // START - added by Corey Oses - May 2017
  // effectively logs EVERYTHING, deals with cout and logger
  void logger(const string& function_name, stringstream& message, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, stringstream& message, ofstream& FileMESSAGE, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, stringstream& message, const string& directory, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, stringstream& message, const string& directory, ofstream& FileMESSAGE, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, stringstream& message, const _aflags& aflags, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, stringstream& message, const _aflags& aflags, ofstream& FileMESSAGE, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, const string& _message, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, const string& _message, ofstream& FileMESSAGE, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, const string& _message, const string& directory, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, const string& _message, const _aflags& aflags, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // overload
  void logger(const string& function_name, const string& _message, const string& directory, ofstream& FileMESSAGE, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // main function
  void logger(const string& function_name, const string& _message, const _aflags& aflags, ofstream& FileMESSAGE, ostream& oss=cout, const char& type=_LOGGER_MESSAGE_, bool silent=false, const string& message_metadata="user, host, time");  // main function
  // END - added by Corey Oses - May 2017
  xstructure LTCELL(string options,istream& input);
  // [OBSOLETE]  xstructure LTCELLFV(string options,istream& input);
  void MagneticParameters(string _directory, ostream& oss);
  // [OBSOLETE]  xstructure MILLER(string options,istream& input);
  xstructure MINKOWSKIBASISREDUCTION(istream& input);
  string MISCIBILITY(vector<string> argv);
  void MOM(istream& input);
  void MSI(istream& input);
  uint NATOMS(istream& input);
  xmatrix<double> GetDistMatrix(const xstructure& a); // CO 171025
  string NBONDXX(istream& input,bool aflowlib_legacy_format=false); // CO 171025
  string NBONDXX(const xstructure& a,bool aflowlib_legacy_format=false); // CO 171025
  xstructure NAMES(vector<string>,istream& input);
  xstructure NANOPARTICLE(istream& input,const xvector<double>& iparams);
  xstructure NIGGLI(istream& input);
  void NDATA(istream& input);
  double NNDIST(istream& input);
  xstructure NOORDERPARAMETER(istream& input);
  xstructure NOSD(istream& input);
  xstructure NUMNAMES(vector<string>,istream& input);
  uint NSPECIES(istream& input);
  bool OPARAMETER(vector<string>,istream& input);
  bool SHIRLEY(vector<string>,istream& input);
  bool SHIRLEY2(vector<string>,istream& input);
  string PEARSON_SYMBOL(istream& input);
  bool POCCUPATION(vector<string>,istream& input);
  void PDB(istream& input);
  void PDOS(vector<string>);
  void PHONONS(_aflags &aflags,istream& input,const double& radius);
  void PGROUP(_aflags &aflags,istream& input);
  void PGROUPXTAL(_aflags &aflags,istream& input);
  void PGROUPK(_aflags &aflags,istream& input);
  void PLANEDENS(vector<string>);
  // [OBSOLETE] string PLATON(vector<string>,istream& input);
  string PLATON(string options,istream& input);
  // DX 9/26/17 [OBSOLETE] string SG(string options,istream& input,string mode,string print);
  string SG(aurostd::xoption& vpflow,istream& input,string mode,string print);
  // [OBSOLETE]  string SG(string mode,string print,vector<string>,istream& input);
  void STATDIEL(vector<string>& argv); // CAMILO
  bool SYMMETRY_GROUPS(_aflags &aflags,istream& input, aurostd::xoption& vpflow, ostream& oss); // DX 8/18/17 - Add no_scan option to all symmetry Xgroups
  void POCC(vector<string>);
  string POSCAR2AFLOWIN(istream& input);
  void POSCAR2WYCKOFF(istream& input);
  bool QMVASP(aurostd::xoption& vpflow);  //vector<string> argv); //CO 180703
  xstructure POSCAR(istream& input);
  xmatrix<double> QE_ibrav2lattice(const int& ibrav, const xvector<double>& parameters, const bool& isabc); // DX 1/23/18 - added more robust QE reader

}
bool RequestedAlphabeticLabeling(string &label);
bool AlphabetizePrototypeLabelSpecies(deque<string> &species,deque<string> &species_pp,deque<double> &volumes,deque<double> &masses,string &label);
bool AlphabetizePrototypeLabelSpecies(deque<string> &species,deque<string> &species_pp,string &label);
bool AlphabetizePrototypeLabelSpecies(deque<string> &species,deque<double> &volumes,string &label);
bool AlphabetizePrototypeLabelSpecies(deque<string> &species,string &label);
string AlphabetizePrototypeLabelSpeciesArgv(vector<string> &argv);
namespace pflow {
  xstructure PROTO_LIBRARIES(aurostd::xoption vpflow);
  bool PROTO_AFLOW(aurostd::xoption vpflow,bool flag_REVERSE);  // too many options
  // [OBSOLETE] xstructure PROTO_HTQC_PURE(vector<string>);
  // [OBSOLETE] xstructure PROTO_GUS_PURE(vector<string>);
  // [OBSOLETE] bool PROTO_GUS_CPP(vector<string>);
  bool PROTO_ICSD_AFLOWIN(vector<string> &argv);
  xstructure PRIM(istream& input,uint mode);
  void RASMOL(string options,istream& input);
  void RBANAL(vector<string>);
  void RBDIST(vector<string>);
  xstructure RMATOM(istream& input,const int& iatom);
  xstructure RMCOPIES(istream& input);
  void RAYTRACE(vector<string>);
  xstructure SCALE(string options,istream& input);
  void RDF(string options,istream& input);
  void RDFCMP(string options);
  void RSM(vector<string>, istream& input);
  xstructure SD(vector<string>,istream& input);
  xstructure SETCM(istream& input,const xvector<double>& cm);
  xstructure SETORIGIN(istream& input,const xvector<double>& origin);
  xstructure SETORIGIN(istream& input,const int& natom);
  void SEWALD(vector<string>,istream& input);
  void SG(istream& input);
  bool SGDATA(istream& input, aurostd::xoption& vpflow, ostream& oss); // DX 8/31/17 - SGDATA
  void SGROUP(_aflags &aflags,istream& input,double radius);
  void SHELL(string options,istream& input);
  string SPECIES(istream& input);
  xstructure SHIFT(string options,istream& input);
  void SPLINE(vector<string>);
  void SUMPDOS(vector<string>);
  xstructure SUPERCELL(string options,istream& input);
  void SUPERCELLSTRLIST(string options);
  xstructure xstrSWAP(vector<string>, istream& input);
  xstructure VOLUME(string options, istream& input);
  string WYCCAR(aurostd::xoption& vpflow,istream& input); //DX 20180807 - added wyccar to pflow
  xstructure WYCKOFF(vector<string>,istream& input);
  void XRAY(string options,istream& input);
  void XYZ(string options,istream& input);
  void XYZINSPHERE(istream& input,double radius);
  void XYZWS(istream& input);
  void ZVAL(string options);
}

// aflow_pflow_print.cpp
namespace pflow {
  void PrintACE(const xstructure&,ostream& oss);
  void PrintAngles(xstructure str,const double& cutoff,ostream& oss);
  class projdata;
  void PrintBands(const pflow::projdata& pd);
  bool PrintCHGCAR(const xstructure& str,const stringstream& chgcar_header,const vector<int>& ngrid,const vector<int>& format_dim,const vector<double>& chg_tot,const vector<double>& chg_diff,const string& output_name,ostream& oss);
  void PrintChgInt(vector<pflow::matrix<double> >& rad_chg_int,pflow::matrix<double>& vor_chg_int,ostream& oss);  
  void PrintCIF(ostream& oss,const xstructure&,int=1,int=1); //DX 20180806 - added setting default
  void PrintClat(const xvector<double>& data,ostream& oss);
  void PrintCmpStr(const xstructure& str1,const xstructure& str2,const double& rcut,ostream& oss);  
  void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc,ostream& oss,string mode, const string& format="txt",bool already_calculated=false); // CO171027
  void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss,string smode, double tolerance, bool no_scan, const int& sg_setting=1, const string& format="txt",bool already_calculated=false); // CO171027
  void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss,string smode, aurostd::xoption& vpflow, const string& format="txt",bool already_calculated=false); //DX 20180823
  void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss_final,string smode, aurostd::xoption& vpflow, double tolerance, bool no_scan, const int& sg_setting=1, const string& format="txt",bool already_calculated=false); //DX 20180822
  // DX 9/1/17 [OBSOLETE] void PrintData(const xstructure&,ostream& oss,string mode,const string& format="txt",bool already_calculated=false);
  void PrintData(const xstructure& str,ostream& oss,string smode,double tolerance, bool no_scan, const int& sg_setting=1, const string& format="txt");
  void PrintData(const xstructure& str,ostream& oss,string smode,const string& format="txt");
  void PrintData1(const xstructure& str1,const double& rcut,ostream& oss);
  string PrintData1(const xstructure& str1,const double& rcut);
  void PrintData2(const xstructure&,ostream& oss);
  void PrintDisplacements(xstructure str,const double cutoff,ostream& oss);
  void PrintDistances(xstructure str,const double cutoff,ostream& oss);
  void PrintEwald(const xstructure& in_str,double& epoint,double& ereal,double& erecip,double& eewald,double& eta,const double& SUMTOL,ostream& oss);
  void PrintGulp(const xstructure&,ostream& oss);
  bool PrintSGData(xstructure& str_sg, ostream& oss, bool standalone=true, const string& format="txt",bool already_calculated=false); // DX 8/30/17 - SGDATA
  bool PrintSGData(xstructure& str_sg, double& tolerance, ostream& oss, bool no_scan, const int& setting=1, bool standalone=true, const string& format="txt",bool already_calculated=false); // DX 2/26/18 - added & to tolerance
  bool PrintSGData(xstructure& str_sg, double& tolerance, ostream& oss_final, aurostd::xoption& vpflow, bool no_scan, const int& sg_setting=1, bool standalone=true, const string& format="txt",bool already_calculated=false); //DX 20180822
}
void PrintKmesh(const xmatrix<double>& kmesh,ostream& oss);    // HERE
void PrintImages(xstructure strA,xstructure strB,const int& ni,const string& path_flag);
void PrintMSI(const xstructure&,ostream& oss);
void PrintNdata(const xstructure&,ostream& oss);
//void PrintNeatProj(projdata& pd);
void PrintPDB(const xstructure&,ostream& oss);
void platon2print(xstructure,bool P_EQUAL,bool P_EXACT,double P_ang,double P_d1,double P_d2,double P_d3,ostream& sout);
void PrintRDF(const xstructure& str,const double& rmax,const int& nbins,const int& smooth_width,const pflow::matrix<double>& rdf_all,
	      pflow::matrix<double>& rdfsh_all,pflow::matrix<double>& rdfsh_loc,ostream& oss);
void PrintRDFCmp(const xstructure& str_A,const xstructure& str_B,const double& rmax,const int nbins,
		 const double& smooth_width,const int nsh,const pflow::matrix<double>& rdfsh_all_A,
		 const pflow::matrix<double>& rdfsh_all_B,const vector<int>& best_match,
		 const pflow::matrix<double>& rms_mat,ostream& oss);
void PrintRSM(const xstructure&,ostream& oss);
void PrintShell(const xstructure& str,const int& ns,const double& rmin,const double& rmax,const string& sname,const int lin_dens,ostream& oss);
void PrintXray(const xstructure& str,const double& l,ostream& oss);
void PrintXYZ(const xstructure& a,const xvector<int>& n,ostream& oss);
void PrintXYZws(const xstructure& a,ostream& oss);
void PrintXYZInSphere(const xstructure& a,const double& radius,ostream& oss);

// aflow_pflow_funcs.cpp
double DebyeWallerFactor(const double& theta,const double& lambda,const double& temp,const double& debye_temp,const double& mass);
xvector<double> balanceChemicalEquation(const vector<xvector<double> >& _lhs,const vector<xvector<double> >& _rhs,
    bool normalize,double tol); //CO 180817
xvector<double> balanceChemicalEquation(const xmatrix<double>& _composition_matrix,bool normalize,double tol);
void ParseChemFormula(string& ChemFormula,vector<string>& ChemName,vector<float>& ChemConc);
void ParseChemFormulaIndividual(uint nchar,string& ChemFormula,string& AtomSymbol,float& AtomConc);

namespace pflow {
  void GetXray(const xstructure& str,vector<double>& dist,vector<double>& sf,const double& lambda,
	       vector<double>& scatt_fact,vector<double>& mass,vector<double>& twoB_vec);
  void GetRDF(xstructure str,const double& rmax,const int& nbins,matrix<double>& rdf_all);
  void GetRDFShells(const xstructure& str,const double& rmax,const int& nbins,const int& smooth_width,
		    const pflow::matrix<double>& rdf,matrix<double>& rdfsh,matrix<double>& rdfsh_loc);
  double RdfSh_RMS(const int iaA,const int iaB,const int nsh_max,const int nt,
		   const pflow::matrix<double>& rdfsh_all_A,const pflow::matrix<double>& rdfsh_all_B);
  void CmpRDFShells(const xstructure& str_A,const xstructure& str_B,const pflow::matrix<double>& rdfsh_all_A,
		    const pflow::matrix<double>& rdfsh_all_B,const int nsh,vector<int>& best_match,
		    pflow::matrix<double>& rms_mat);
  pflow::matrix<double> GetSmoothRDF(const pflow::matrix<double>& rdf,const double& sigma);
  void CmpStrDist(xstructure str1,xstructure str2,const double& cutoff,
		  pflow::matrix<double>& dist1,pflow::matrix<double>& dist2,
		  pflow::matrix<double>& dist_diff,matrix<double>& dist_diff_n);
}

// aflow_pflow.cpp
int SignNoZero(const double& x);
int Nint(const double& x);
int Sign(const double& x);

// ---------------------------------------------------------------------------
// PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDO
namespace pflow {
  class pdosdata
  {
  public:
    // Constructors
    pdosdata(); // default
    void PrintParams(ostream& oss, const vector<string>& Ltotnames);
    void PrintPDOS(ostream& oss, const int& sp);
    // variables
    string PDOSinfile;
    pflow::matrix<int> pdos_at;
    pflow::matrix<int> pdos_k;
    pflow::matrix<int> pdos_b;
    pflow::matrix<int> pdos_lm;
    pflow::matrix<double> pdos;
    double emin;
    double emax;
    double efermi;
    double smooth_sigma;
    int spin;
    int nlm;
    int natoms;
    int print_params;
    int nbins;
  };
}

// ---------------------------------------------------------------------------
// RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY
namespace pflow {
  class rtparams
  {
  public:
    void free();
    void copy(const rtparams& b);
    // Constructors
    rtparams(); // default
    rtparams(const rtparams& b); // default
    //Operators
    const rtparams& operator=(const rtparams& b);
    // Accessors
    void Print() const;
    // variables
    double resx;
    double resy;
    vector<double> viewdir;
    int viewdir_s;
    vector<double> updir;
    int updir_s;
    double zoom;
    double aspectratio;
    double antialiasing;
    double raydepth;
    vector<double> center;
    vector<double> center_guide;
    int center_s;
    vector<double> background;
    pflow::matrix<double> lightcenter;
    vector<double> lightrad;
    pflow::matrix<double> lightcolor;
    // Sphere texture variables (ambient, diffuse, specular, opacity)
    pflow::matrix<double> sphtex_tex;
    vector<double> sphtex_tex_def;
    pflow::matrix<double> sphtex_color;
    vector<double> sphtex_color_def;
    pflow::matrix<double> sphtex_phong;
    vector<double> sphtex_phong_def;
    vector<string> sphtex_names;
    vector<double> sph_rad;
    // Plane variables
    int plane;
    int plane_s;
    vector<double> plane_center;
    vector<double> plane_normal;
    vector<double> plane_color;
    vector<double> planetex_tex;
    vector<double> plane_center_def;
    vector<double> plane_normal_def;
    vector<double> plane_color_def;
    vector<double> planetex_tex_def;
    int plane_center_s;
    int plane_normal_s;
    int plane_color_s;
    int planetex_tex_s;
    double sph_rad_def;
    string shading;
    string outfile;
    pflow::matrix<double> sc;
    int sc_s;
    int calc_type;
    vector<string> input_files;
    int first_set;
    string insert_file;
    vector<double> rotation;
    vector<double> struct_origin;
    int struct_origin_s;
  };

  void SetRTParams(xstructure& str, pflow::rtparams& rtp);
  vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(ifstream& outcar_inf, ifstream& xdatcar_inf);
  void GetDatFromOutcar(vector<pflow::matrix<double> >& lat_vec, deque<int>& num_each_type, ifstream& outcar_inf);
  void GetDatFromXdatcar(vector<pflow::matrix<double> >& fpos_vec, ifstream& xdatcar_inf);
  vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(ifstream& outcar_inf, ifstream& xdatcar_inf);
  void PrintStrVec(const vector<xstructure>& str_vec, ostream& outf);
  void ReadInRTParams(ifstream& rtinfile, pflow::rtparams& rtp);
  void ReadInStrVec(vector<xstructure>& str_vec, ifstream& strlist_inf);
  void JoinStrVec(vector<xstructure> str_vec_1,vector<xstructure> str_vec_2,vector<xstructure>& str_vec_tot);
  void SetStrFromRTParams(xstructure& str, pflow::rtparams& rtp);
  void SuperCellStrVec(vector<xstructure>& str_vec, const pflow::matrix<double>& sc);
  void UpDateRTParams(pflow::rtparams& rtp, const int& istr, int nstr);
  void SetRTParams(xstructure& str, pflow::rtparams& rtp);
  void GetRTDatFile(xstructure str, const pflow::rtparams& rtp, ostringstream& rtdat_file);
  string PrintRTDatFile(ostringstream& rtdat_file, const pflow::rtparams& rt_params);
  string CreateRTtgaFile(const string& datfile, const pflow::rtparams& rt_params);
  string CreateRTjpgFile(const string& tgafile, const pflow::rtparams& rt_params);  
  void GetRTencFile(const pflow::rtparams& rtp, const int nim,ostringstream& os);
  string PrintRTencFile(const pflow::rtparams& rt_params, ostringstream& rtenc_file);
  string CreateRTmpgFile(const pflow::rtparams& rt_params, const string& encfile);
  void RayTraceManager(vector<string>);
  pflow::matrix<double> GetRotationMatrix(const vector<double>& angles);
  void RotateStrVec(vector<xstructure>& str_vec, const vector<double>& rot);

}

// ---------------------------------------------------------------------------
// PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PRO

namespace pflow {
  class projdata
  {
  public:
    // Constructors
    projdata(); // default
    void Print(ostream& outf);
    // variables
    int nl_max; // 4 for s,p,d,f orbitals
    int nlm_max; // 16 for s,p,d,f orbitals
    int nlmtot_max; // 20 for s,p,d,f orbitals + p,d,f,all totals
    int nl; // 3 for spd
    int nlm; // 9 for spd
    int nlmtot; // 9+Psum+Dsum+Allsum=12 (for spd)
    int nkpts;
    int nbands;
    int nions;
    int ntypes;
    vector<int> num_each_type;
    pflow::matrix<double> wfermi_u;
    pflow::matrix<double> wfermi_d;
    vector<double> wkpt;
    vector<pflow::matrix<pflow::matrix<std::complex<double> > > > pdat_u;
    vector<pflow::matrix<pflow::matrix<std::complex<double> > > > pdat_d;
    pflow::matrix<pflow::matrix<double> > occ_vs_ion_kpt_bnd_lm_u;
    pflow::matrix<pflow::matrix<double> > occ_vs_ion_kpt_bnd_lm_d;
    pflow::matrix<pflow::matrix<double> > occ_vs_ion_kpt_bnd_l_u;
    pflow::matrix<pflow::matrix<double> > occ_vs_ion_kpt_bnd_l_d;
    pflow::matrix<pflow::matrix<double> > occ_vs_ion_kpt_bnd_lmtot_u;
    pflow::matrix<pflow::matrix<double> > occ_vs_ion_kpt_bnd_lmtot_d;
    vector<pflow::matrix<double> > occ_vs_ion_kpt_lm_u;
    vector<pflow::matrix<double> > occ_vs_ion_kpt_lm_d;
    vector<pflow::matrix<double> > occ_vs_ion_kpt_l_u;
    vector<pflow::matrix<double> > occ_vs_ion_kpt_l_d;
    vector<pflow::matrix<double> > occ_vs_ion_bnd_lm_u;
    vector<pflow::matrix<double> > occ_vs_ion_bnd_lm_d;
    vector<pflow::matrix<double> > occ_vs_ion_bnd_l_u;
    vector<pflow::matrix<double> > occ_vs_ion_bnd_l_d;
    pflow::matrix<double> occ_vs_ion_lm_u;
    pflow::matrix<double> occ_vs_ion_lm_d;
    pflow::matrix<double> occ_vs_ion_l_u;
    pflow::matrix<double> occ_vs_ion_l_d;
    pflow::matrix<double> ener_k_b_u;
    pflow::matrix<double> ener_k_b_d;
    int sp;
    int rspin;
    pflow::matrix<double> kpts;
    vector<string> LMnames;
    vector<string> Lnames;
    vector<string> LLMnames;
    string PROOUTinfile;
    pflow::matrix<double> lat;
  };
}

// ---------------------------------------------------------------------------
// PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJF

// in aflow_pflow_funcs.cpp
namespace pflow {
  std::complex<double> ProcessProjection(const std::complex<double>& proj);
  void ReadInProj(projdata& pd);
  void CalcNeatProj(projdata& pd, int only_occ);
  void ReadInPDOSData(const projdata& prd, pdosdata& pdd);
  void CalcPDOS(const projdata& prd, pdosdata& pdd);
  void SmoothPDOS(const projdata& prd, pdosdata& pdd);
  void AtomCntError(const string& tok, const int tokcnt, const int atom_cnt);
}

// in aflow_pflow_print.cpp
void PrintNeatProj(pflow::projdata& pd, ostream& outf);

// ---------------------------------------------------------------------------
// SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUN

// in aflow_pflow_funcs.cpp
namespace pflow {
  void ReadSumDOSParams(ifstream& infile, pflow::pdosdata& pdd);
  void ReadInPDOSData(matrix<pflow::matrix<double> >& allpdos, pflow::pdosdata& pdd,ifstream& infile);
  void SumPDOS(const pflow::matrix<pflow::matrix<double> >& allpdos, pflow::pdosdata& pdd);
}

// in pflow_print
void PrintSumPDOS(pflow::pdosdata& pdd, ostream& outf);

// ---------------------------------------------------------------------------
// RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBF

// in aflow_pflow_funcs.cpp
namespace pflow {
  double TotalAtomDist(xstructure str, xstructure str00, const string& path_flag);
  vector<string> GetRBDir(const int& nim);
  vector<double> GetRBEner(const int& nim);
  vector<xstructure> GetRBStruct(const int& nim);
  vector<double> GetRBDistCum(const vector<xstructure>& str_vec, const string& path_flag);
  vector<double> GetRBDistFromStrI(const vector<xstructure>& str_vec,const xstructure& strI,const string& path_flag);
  void RBPoscarDisp(const xstructure& str1in, const xstructure& str2in,xstructure& diffstr, double& totdist, pflow::matrix<double>& cm,const string& path_flag);
}

// in aflow_pflow_print.cpp
void PrintRBAnal(const int& nim, const string& path_flag, ostream& outf);
void PrintRBPoscarDisp(const xstructure& diffstr, double& totdist, pflow::matrix<double>& cm, const string& path_flag, ostream& outf);

// ---------------------------------------------------------------------------
// CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUN

namespace pflow {
  class pd_params {
  public:
    string type;
    double scale;
    pflow::matrix<double> pts;
    pflow::matrix<double> dpts;
    int Nx,Ny;
    string orig_loc;
    string ortho;
    void Print(ostream& outf) const;
  };
  
  bool ReadCHGCAR(xstructure& str,stringstream& chgcar_header, vector<int>& ngrid, vector<int>& format_dim, vector<double>& chg_tot,
                  vector<double>& chg_diff, stringstream& chgcar_ss,ostream& oss);
  bool ReadChg(xstructure& str,vector<int>& ngrid, vector<double>& chg_tot,
               vector<double>& chg_diff, istream& chgfile);
  void GetChgInt(vector<pflow::matrix<double> >& rad_chg_int, pflow::matrix<double>& vor_chg_int,
		 xstructure& str,vector<int>& ngrid,vector<double>& chg_tot, vector<double>& chg_diff);
  void ReadPlaneDensParams(const xstructure& str, pd_params& pdp, istream& infile);
  void GetPlaneDens(const pd_params& pdp, vector<double>& dens2d_tot, vector<double>& dens2d_diff,
		    const xstructure& str, const vector<int>& ngrid,
		    const vector<double>& chg_tot, const vector<double>& chg_diff);
  void PrintPlaneDens(const pd_params& pdp, const vector<double>& dens2d_tot,
		      const vector<double>& dens2d_diff, const xstructure& str);
}

// ---------------------------------------------------------------------------
// EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWA
namespace pflow {
  void Ewald(const xstructure& in_str,double& epoint,double& ereal,
	     double& erecip,double& eewald,double& eta,const double& SUMTOL);
  double GetEta(const int& natoms,const double& vol);
  double GetPointEner(const double& rteta,const vector<double>& atchg,const double& vol);
  double GetRecipEner(const double& eta,const vector<double>& atchg,const double& vol,
		      const pflow::matrix<double>& rlat,const pflow::matrix<double>& fpos,const double& SUMTOL);
  double GetRealEner(const double& eta,const vector<double>& atchg,const double& vol,
		     const pflow::matrix<double>& lat,const pflow::matrix<double>& fpos,const double& SUMTOL);
  double GetScreenedESEner(void);
  double ScreenedESEner(const xstructure& in_str,const double& Ks,const double& SUMTOL);
}

// Output help information of an option
void helpIndividualOption(vector<string> & argv);

// ---------------------------------------------------------------------------
// FORMER WAHYU.H

void AConvaspBandgap(vector<string>& bandsdir);
void AConvaspBandgaps(istream& bandsdir,ostream& oss);
void AConvaspBandgaps(istream& bandsdir,ostringstream& oss);
void AConvaspBandgapFromDOS(istream& doscar);
void AConvaspBandgapListFromDOS(istream& doscar);
namespace pflow {
  void ICSD(vector<string> argv, istream& input);
  void ICSD_CheckRaw(vector<string> argv);
  void ICSD_2POSCAR(istream& input);
  void ICSD_2PROTO(istream& input);
  void ICSD_2WYCK(istream& input,bool SOF);
  void ICSD_ListMetals();
}
float GetBandGap_WAHYU(stringstream& straus,float Efermi,char& gaptype);
float GetBandgapFromDOS(ifstream& doscar);
float GetBandgapFromDOS(istream& doscar);
vector<string> GetMetalsList(bool v);
bool IsMetal(const string s);
void ParseChemicalFormula(string Formula,vector<string>& Zname,vector<float>& Zconc);
string RemoveCharFromString(const string s, char c);
string RemoveStringFromTo(const string s, char cstart, char cstop);
int StringCrop(string s,vector<string>& vstring);
string StringCropAt(const string s,int icrop);
vector<float> SortFloat(vector<float> v, int mode);
xvector<double> cross(const xvector<double> a, const xvector<double> b);


// ---------------------------------------------------------------------------
// FORMER RICHARD.H
namespace pflow {
  double GetAtomicPlaneDist(string options,istream& input);
  double frac2dbl(string str);
  bool havechar(string str_in, char c);
  int whereischar(string str, char c);
  void cleanupstring(string & str);
} // namespace pflow

//from kesong's old files
namespace pflow {
  void BZMAX(istream& input);
}

#endif
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

