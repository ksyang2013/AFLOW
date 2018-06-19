// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AFLOWLIB_H_
#define _AFLOWLIB_H_

#include "aflow.h"
//[OBSOLETE] [KESONG] #include "aflow_contrib_kesong.h" //CO 180515

using std::vector;
using std::string;

#define NOSG string("NNN #0")
#define CHMODWEB FALSE
#define INF 1E9
#define ENERGY_ATOM_ERROR_meV 50
#define PRESSURE_ZERO_ENTHALPY_ENERGY 1e-6

extern vector<string> vLibrary_LIB2;   // need to remove somwhoe
extern vector<string> vLibrary_LIB3;
extern vector<string> vLibrary_LIB1;
extern vector<string> vLibrary_ICSD;
extern vector<string> vLibrary_AURO;
extern vector<vector<string> > vLibrary_LIB2_tokens;
extern vector<vector<string> > vLibrary_LIB3_tokens;
extern vector<vector<string> > vLibrary_LIB1_tokens;
extern vector<vector<string> > vLibrary_ICSD_tokens;
extern vector<vector<string> > vLibrary_AURO_tokens;

// namespace for the web
// ***************************************************************************
namespace aflowlib {
  class _aflowlib_entry {
  public:
    // constructor destructor                                 // constructor/destructor
    _aflowlib_entry();                                         // default, just allocate
    ~_aflowlib_entry();                                        // kill everything
    _aflowlib_entry(const _aflowlib_entry& b);                  // constructor copy
    const _aflowlib_entry& operator=(const _aflowlib_entry &b); // copy
    // CONTROL
    string entry;vector<string> ventry;                       // ventry split by "|"
    string auid;                                              // AFLOW UNIQUE IDENTIFIER
    string aurl;vector<string> vaurl;                         // AFLOW RESEARCH LOCATOR and TOKENS
    string keywords;vector<string> vkeywords;                 // keywords inside
    string aflowlib_date,aflowlib_version;                    // version/creation
    string aflowlib_entries;vector<string> vaflowlib_entries; // this contains the subdirectories that can be associated
    int aflowlib_entries_number;                              // their number
    string aflow_version;                                     // version/creation
    string catalog;                                           // ICSD,LIB2, etc.
    string data_api,data_source,data_language;                // version/source/language
    string error_status;                                            // ERROR ??
    string author;vector<string> vauthor; 
    int calculation_cores;double calculation_memory,calculation_time;
    string corresponding;vector<string> vcorresponding; 
    string loop;vector<string> vloop;                         // postprocessing
    int node_CPU_Cores;                                       // computer
    double node_CPU_MHz;                                      // computer
    string node_CPU_Model;                                    // computer
    double node_RAM_GB;                                       // computer
    // materials
    string Bravais_lattice_orig,Bravais_lattice_relax;        // structure
    string code;
    string composition;vector<double> vcomposition;
    string compound;
    double density;
    string dft_type;vector<string> vdft_type;
    double eentropy_cell,eentropy_atom;
    double Egap,Egap_fit;
    string Egap_type;  
    double energy_cell,energy_atom,energy_atom_relax1;
    double energy_cutoff;
    double delta_electronic_energy_convergence;
    double delta_electronic_energy_threshold;
    uint nkpoints,nkpoints_irreducible,kppra;
    string kpoints;
    xvector<int> kpoints_nnn_relax,kpoints_nnn_static;
    vector<string> kpoints_pairs;
    double kpoints_bands_path_grid;
    double enthalpy_cell,enthalpy_atom;
    double enthalpy_formation_cell,enthalpy_formation_atom;
    double entropic_temperature;
    string files;vector<string> vfiles;
    string files_LIB;vector<string> vfiles_LIB;
    string files_RAW;vector<string> vfiles_RAW;
    string files_WEB;vector<string> vfiles_WEB;
    string forces;vector<xvector<double> > vforces;
    string geometry;vector<double> vgeometry; // a,b,c and unit_cell_angles (b,c) (a,c) (a,b)
    string lattice_system_orig,lattice_variation_orig,lattice_system_relax,lattice_variation_relax;
    string ldau_TLUJ;
    uint natoms;
    string nbondxx;vector<double> vnbondxx;
    uint nspecies;
    string Pearson_symbol_orig,Pearson_symbol_relax;
    string positions_cartesian;vector<xvector<double> > vpositions_cartesian;
    string positions_fractional;vector<xvector<double> > vpositions_fractional;
    double pressure;           // the true applied pressure (PSTRESS)
    string stress_tensor;vector<double> vstress_tensor; // (1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3)
    double pressure_residual;  // the leftover pressure due to convergence 
    double Pulay_stress;       // the leftover pressure due to incomplete basis set
    string prototype;
    double PV_cell,PV_atom;
    double scintillation_attenuation_length;
    string sg,sg2;vector<string> vsg,vsg2; // CO 180101
    string spacegroup_orig,spacegroup_relax;
    string species;vector<string> vspecies;
    string species_pp;vector<string> vspecies_pp;
    string species_pp_version;vector<string> vspecies_pp_version;
    string species_pp_ZVAL;vector<double> vspecies_pp_ZVAL;
    double spin_cell,spin_atom;
    string spinD;vector<double> vspinD;
    string spinD_magmom_orig;vector<double> vspinD_magmom_orig;
    double spinF;
    string sponsor;vector<string> vsponsor; 
    string stoichiometry;vector<double> vstoichiometry;
    double valence_cell_std,valence_cell_iupac;
    double volume_cell,volume_atom;
    // AGL/AEL
    double agl_thermal_conductivity_300K;//  (W/m*K)
    double agl_debye;//  (K)
    double agl_acoustic_debye;//  (K)
    double agl_gruneisen;// 
    double agl_heat_capacity_Cv_300K;//  (kB/cell)
    double agl_heat_capacity_Cp_300K;//  (kB/cell)
    double agl_thermal_expansion_300K;//  (1/K)
    double agl_bulk_modulus_static_300K;//  (GPa)
    double agl_bulk_modulus_isothermal_300K;//  (GPa)
    double ael_poisson_ratio ;//
    double ael_bulk_modulus_voigt;//  (GPa)
    double ael_bulk_modulus_reuss;//  (GPa)
    double ael_shear_modulus_voigt;//  (GPa)
    double ael_shear_modulus_reuss;//  (GPa)
    double ael_bulk_modulus_vrh;//  (GPa)
    double ael_shear_modulus_vrh;//  (GPa)
    double ael_elastic_anistropy;//
    // BADER
    string bader_net_charges;vector<double> vbader_net_charges;//electrons
    string bader_atomic_volumes;vector<double> vbader_atomic_volumes;//Angst^3
    // legacy
    string server;vector<string> vserver;vector<vector<string> > vserverdir; 
    string icsd;
    string stoich;vector<double> vstoich;
    // apennsy
    string structure_name;                            // name for apennsy
    string structure_description;                     // description for apennsy
    double distance_gnd;                              // distance_gnd
    double distance_tie;                              // distance_tie
    bool pureA,pureB;                                 // pureA,pureB
    bool fcc,bcc,hcp;                                 // options for lattices
    double stoich_a,stoich_b;                         // stoich_b,stoich_b
    double bond_aa,bond_ab,bond_bb;                   // bond_xx // BOND_XX [norm V_ATOM^0.33]
    vector<uint> vNsgroup;                            // vNsgroups
    vector<string> vsgroup;                           // vsgroups
    vector<xstructure> vstr;                          // vstructures
    bool FixDescription(void);                        // fix description names
    void GetSGROUP(string aflowlibentry);             // disassemble SG
    // functions
    void clear();                                              // free space
    uint Load(const stringstream& stream,ostream& oss);        // load from stringstream it std is cout
    uint Load(const string& entry,ostream& oss);               // load from string it std is cout
    uint file2aflowlib(const string& file,ostream& oss);       // load from file
    uint url2aflowlib(const string& url,ostream& oss,bool=TRUE); // load from the web (VERBOSE)
    string aflowlib2string(string="out");                      //
    string aflowlib2file(string file,string="out");            //
    void correctBadDatabase(bool verbose=true,ostream& oss=cout);                             // CO 171202 - apennsy fixes
    void correctBadDatabase(ofstream& FileMESSAGE,bool verbose=true,ostream& oss=cout);       // CO 171202 - apennsy fixes
    bool ignoreBadDatabase() const;                                                           // CO 171202 - apennsy fixes
    bool ignoreBadDatabase(string& reason) const;                                             // CO 171202 - apennsy fixes
    string getPathAURL(ostream& oss=cout, bool load_from_common=false);                   // converts entry.aurl to url/path (common)
    string getPathAURL(ofstream& FileMESSAGE, ostream& oss, bool load_from_common=false); // converts entry.aurl to url/path (common)
  private:                                                     //
    void free();                                               // free space
    void copy(const _aflowlib_entry& b);                       //
   };
}

// ***************************************************************************
// AFLUX STUFF
namespace aflowlib {
  class APIget {
    private:
      struct sockaddr_in client;
      int sock;
      int PORT;// = 80; // CO 180401
      string Summons;
      string API_Path;
      string Domain;

      bool establish();
    public:
      APIget( string a_Summons="", string a_API_Path="/search/API/?", string a_Domain="aflowlib.duke.edu" ): Summons(a_Summons), API_Path(a_API_Path), Domain(a_Domain) {};
      void reset( string a_Summons="#", string a_API_Path="", string a_Domain="" );
      friend ostream& operator<<( ostream& output, APIget& a );
  };
}

// ***************************************************************************
namespace aflowlib {
  string directory2auid(string directory,string aurl);
  uint auid2present(string="");
  bool AflowlibLocator(const string& in,string& out,const string& mode);
  string AflowlibLocator(string options,string mode);
  uint WEB_Aflowlib_Entry_PHP(string options,ostream& oss);
}

// ***************************************************************************
// to create/analyze the tokens and load the stuff up for qhull
namespace aflowlib {
  bool TokenPresentAFLOWLIB(string line,string query);
  string TokenExtractAFLOWLIB(string line,string query);
  bool TokenExtractAFLOWLIB(string line,string query,string &value);
  bool TokenExtractAFLOWLIB(string line,string query,int &value);
  bool TokenExtractAFLOWLIB(string line,string query,uint &value);
  bool TokenExtractAFLOWLIB(string line,string query,float &value);
  bool TokenExtractAFLOWLIB(string line,string query,double &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<string> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<int> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<uint> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<float> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<double> &value);
}

namespace aflowlib {
  uint LOAD_Library_LIBRARY(string file,string="",bool=FALSE);
  uint LOAD_Library_ALL(string options,string="",bool=FALSE);
  uint LOAD_Library_ALL(string options,bool=FALSE);
  uint LOAD_Library_ALL(bool=FALSE);
}

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
			vector<vector<double> > &vList_Ef);              // OUT [0,naries[*[0,vList.size()[ returns the formation energy of the list
  uint GREP_Species_ALL(string species,                                  // IN  the species Ag,...   nspecies=number of these items
			string& species_pp,                              // IN  the pseudopotentials Ag_pv,
			vector<string>& vList,                           // OUT [0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
			vector<vector<string> > &vList_tokens,           // OUT [0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
			double& List_Hmin,                               // OUT returns the min enthalpy for reference
			string& List_Pmin,                               // OUT returns the prototype for reference
			uint& List_Imin);                                // OUT returns the line index of vList, in which we have the min enthalpy for reference
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList,vector<vector<vector<string> > > &vList_tokens);
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList);
  uint GREP_Species_ALL(string species,                                  // IN  the species Ag,...   nspecies=number of these items
			string& species_pp,                              // IN  the pseudopotentials Ag_pv,
			double& List_Hmin,                               // OUT returns the min enthalpy for reference
			string& List_Pmin);                              // OUT returns the prototype for reference
}

namespace aflowlib {
  bool ALIBRARIES(string options);
  bool LIBS_EFormation(string options);
}

// ***************************************************************************
// LIBS to RAWS of each entry
namespace aflowlib {
  uint GetSpeciesDirectory(string directory,vector<string>& vspecies);
}
namespace aflowlib {
  void XFIX_LIBRARY_ALL(string LIBRARY_IN,vector<string>);
  // void LIB2RAW_LIBRARY_ALL(string LIBRARY_IN);
  string LIB2RAW_CheckProjectFromDirectory(string directory);
  bool LIB2RAW_ALL(string options,bool overwrite);
  bool LIB2RAW_FileNeeded(string directory_LIB,string fileLIB,string directory_RAW,string fileRAW,vector<string> &vfiles,string MESSAGE);
  // [OBSOLETE] bool LIB2RAW(vector<string> argv,bool overwrite);
  bool LIB2RAW(string options,bool overwrite,bool LOCAL=false);
  bool XPLUG(vector<string> argv);
  bool AddFileNameBeforeExtension(string _file,string addendum,string& out_file); // CO 171025
  bool LIB2RAW_Loop_Thermodynamics(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,string MESSAGE,bool LOCAL=false);
  // [OBSOLETE]  bool LIB2RAW_Loop_DATA(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry& data,string MESSAGE);
  bool LIB2RAW_Loop_Bands(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,string MESSAGE);
  bool LIB2RAW_Loop_Magnetic(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,string MESSAGE);
  bool LIB2RAW_Loop_Bader(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,string MESSAGE);
  bool LIB2RAW_Loop_AGL(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,string MESSAGE);
  bool LIB2RAW_Loop_AEL(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,string MESSAGE);
  bool LIB2RAW_Loop_LOCK(string& directory_LIB,string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry& data,string MESSAGE);
}

namespace aflowlib {
  bool VaspFileExist(const string& str_dir, const string& FILE);
  string vaspfile2stringstream(const string& str_dir, const string& FILE, stringstream& streamFILE);
  string vaspfile2stringstream(const string& str_dir, const string& FILE);
}

// ***************************************************************************
#define HTRESOURCE_MODE_NONE  0
#define HTRESOURCE_MODE_PHP_AUTHOR   4
#define HTRESOURCE_MODE_PHP_THRUST   5
#define HTRESOURCE_MODE_PHP_ALLOY    6

// ***************************************************************************
// _OUTREACH CLASS
class _outreach {
 public:
  // constructors/destructors
  _outreach(void);    // do nothing
  _outreach(const _outreach& b);    // do nothing
  ~_outreach();        // do nothing
  // OPERATORS                                              // --------------------------------------
  const _outreach& operator=(const _outreach& b);             // some operators
  void clear(void);                                         // clear
  // CONTENT
  // [OBSOLETE] string print_mode; // "TXT","LATEX","HTML"
  // [OBSOLETE] inside  XHOST.vflag_control.flag("PRINT_MODE::TXT");
  // [OBSOLETE] inside  XHOST.vflag_control.flag("PRINT_MODE::LATEX");
  // [OBSOLETE] inside  XHOST.vflag_control.flag("PRINT_MODE::HTML");
  // [OBSOLETE] bool print_doi; inside  XHOST.vflag_control.flag("PRINT_MODE::DOI");
  // [OBSOLETE] bool print_pdf; inside  XHOST.vflag_control.flag("PRINT_MODE::PDF");
  // [OBSOLETE] bool print_wnumber; inside  XHOST.vflag_control.flag("PRINT_MODE::NUMBER");
  uint wnumber;
  bool newflag;
  uint year;
  vector<string> vauthor;
  string title;
  string journal,link,arxiv,supplementary,bibtex;
  string place,date;
  string type;   // ARTICLE PRESENTATION_TALK PRESENTATION_SEMINAR PRESENTATION_COLLOQUIUM PRESENTATION_KEYNOTE PRESENTATION_PLENARY PRESENTATION_TUTORIAL PRESENTATION_CONTRIBUTED PRESENTATION_POSTER
  bool _isinvited;       // YES
  string host;           // who invited
  string abstract;       // if available and in LaTeX
  string pdf;
  string doi;
  vector<string> vextra_html,vextra_latex;
  vector<string> vkeyword,vsponsor,valloy;
  // operators/functions                                    // operator/functions
  friend ostream& operator<<(ostream &,const _outreach&);   // print
 private:                                                   // ---------------------------------------
  void free();                                              // to free everything
  void copy(const _outreach& b);                            // the flag is necessary because sometimes you need to allocate the space.
};

uint voutreach_load(vector<_outreach>& voutreach,string what2print);
void voutreach_print(uint _mode,ostream& oss,string what2print);
void voutreach_print_everything(ostream& oss,const vector<string>& vitems,string msg1,string msg2,string sectionlabel);

// sort

class _sort_outreach_outreach_year_ {
 public:
  bool operator()(const _outreach& a1, const _outreach& a2) const {
    return (bool) (a1.year<a2.year);}}; // sorting through reference
class _rsort_outreach_outreach_year_ {
 public:
  bool operator()(const _outreach& a1, const _outreach& a2) const {
    return (bool) (a1.year>a2.year);}}; // sorting through reference
class _sort_outreach_outreach_wnumber_ {
 public:
  bool operator()(const _outreach& a1, const _outreach& a2) const {
    return (bool) (a1.wnumber<a2.wnumber);}}; // sorting through reference
class _rsort_outreach_outreach_wnumber_ {
 public:
  bool operator()(const _outreach& a1, const _outreach& a2) const {
    return (bool) (a1.wnumber>a2.wnumber);}}; // sorting through reference

uint voutreach_sort_year(vector<_outreach>& voutreach);
uint voutreach_rsort_year(vector<_outreach>& voutreach);
uint voutreach_sort_wnumber(vector<_outreach>& voutreach);
uint voutreach_rsort_wnumber(vector<_outreach>& voutreach);

uint voutreach_remove_duplicate(vector<_outreach>& voutreach);

// automatic load up
bool ProcessPhpLatexCv(void);
// ***************************************************************************
void center_print(uint mode,ostream& oss);

// ***************************************************************************
// references search
void SystemReferences(const string& system_in,vector<string>& vrefs,bool=FALSE);  // if true then put et_al
void SystemReferences(const string& system_in,vector<uint>& vwnumber); 
bool SystemInSystems(const string& system,const string& systems);
bool SystemInSystems(const string& system,const vector<string>& vsystems);
bool AlloyAlphabeticLIBRARY(const string& system);

// ***************************************************************************
// OLD STUFF FOR SECURITY
// bool ProcessSecurityOptions(vector<string> argv,vector<string> cmds);
// void Aflowlib_AUROHOUSE(void);
// void Aflowlib_AUROHOUSE(int deltat);
// void Aflowlib_AUROHOUSE(int deltat,bool CURATOR);
// ***************************************************************************

#endif //  _AFLOWLIB_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
