// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2018             *
// *                 Kesong Yang - Duke University 2010-2011                 *
// *                  Corey Oses - Duke University 2013-2018                 *
// *                                                                         *
// ***************************************************************************
// aflow_pocc.h and aflow_pocc*.cpp*
//
// completely revised approach to KESONG YANG's original implementation
// no recursion needed --- too memory intensive
// issues with UFF (structure comparison), multiply occupied sites, and 
// vacancies are all addressed robustly here
// 2013-2017: corey.oses@duke.edu
// 2010-2011: kesong.yang@gmail.com (LEGACY)

#ifndef _AFLOW_POCC_H_
#define _AFLOW_POCC_H_

//precision defines tol and paddings
//const int _AFLOW_POCC_PRECISION_ = 8;    //moved to aflow.h
//tolerances
//const double _AFLOW_POCC_ZERO_TOL_ = pow(10,-_AFLOW_POCC_PRECISION_);  //moved to aflow.h

//uff param modes
const uint BOND_MODE    = 0;
const uint NONBOND_MODE = 1;

//files
const string UNIQUE_DERIVATIVE_STRUCTURES_FILE="aflow_pocc_structures_unique.out";
const string ALL_DERIVATIVE_STRUCTURES_FILE="aflow_pocc_structures_all.out";
const string AFLOW_POCC_TAG="[AFLOW_POCC]";

//elements
//to make POCC chemistry-independent, calculate UFF based on standard set of elements
//we should not need to have more than 30 species, but if you want to add to list
//I recommend adding IN ORDER so that elements become more different left to right (like periodic table)
const string STD_ELEMENTS_LIST="Sc Ti V Cr Mn Fe Co Ni Cu Zn "
                               "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
                               "La Hf Ta W Re Os Ir Pt Au Hg ";

namespace pocc {
  bool structuresGenerated();
  xstructure extractPARTCAR(const string& AflowIn);
  
  void generateStructures(const string& AflowIn, ostream& oss=cout);
  void generateStructures(const string& AflowIn, ofstream& FileMESSAGE, ostream& oss=cout);
  void generateStructures(const _xvasp& xvasp,const string& AflowIn,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags, ofstream& FileMESSAGE, ostream& oss=cout);

	void updateProgressBar(unsigned long long int current, unsigned long long int end, ostream& oss=cout);
  vector<string> getElementsList();

  class POccSiteConfiguration;    //specify later
 
  //bool sortAtoms(const _atom& a1,const _atom& a2);
  //bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,vector<uint>& v_config_order,vector<int>& v_config_iterators,vector<vector<int> >& v_types_config);
  bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,vector<int>& v_config_iterators,vector<vector<int> >& v_types_config);
  bool getNextSiteConfiguration(vector<uint>& v_config_order,vector<vector<int> >& v_types_config);
  bool getNextSiteConfiguration(vector<int>& site_config);
  void getConfigOrder(vector<vector<int> >& v_types_config,vector<uint>& v_config_order);

  //simple structure for sorting by vacancy count
  struct SiteVacancyCount {
    uint site;
    uint vacancy_count;
    bool operator<(const SiteVacancyCount& other) const {return vacancy_count<other.vacancy_count;}
  };

  struct StructureConfiguration {
    vector<POccSiteConfiguration> site_configs;
    double max_stoich_error;
    double max_site_error;
    //xvector<double> site_errors;
  };

  struct POccSite {
    uint site;                                        //reflexive
    bool partial_occupation_flag;
    vector<uint> v_occupants;
    vector<uint> v_types;
  };

  struct POccGroup{                                   //groupings of occupants on a site that have same pocc values (optimize together)
    uint site;                                        //reflexive
    double partial_occupation_value;
    vector<uint> v_occupants;
    vector<uint> v_types;
    bool operator<(const POccGroup& other) const {return partial_occupation_value<other.partial_occupation_value;}
  };

  struct POccReducedDerivativeStructure{
    unsigned long long int hnf_index;
    unsigned long long int site_config_index;
    double energy;
    bool operator<(const POccReducedDerivativeStructure& other) const {return energy<other.energy;}
    const POccReducedDerivativeStructure& operator=(const POccReducedDerivativeStructure& b){ //fix me later
      if(this!=&b){
        hnf_index=b.hnf_index;
        site_config_index=b.site_config_index;
        energy=b.energy;
      }
      return *this;
    }
  };
  
  struct POccReducedDerivativeStructuresSet{
    const POccReducedDerivativeStructuresSet& operator=(const POccReducedDerivativeStructuresSet& b){
      prds_set.clear(); prds_set=b.prds_set;
      return *this;
    }
    bool operator<(const POccReducedDerivativeStructuresSet& other) const {return getEnergy()<other.getEnergy();}
    vector<POccReducedDerivativeStructure> prds_set;
    //figure out what to do if vector is not filled, I think it should always be filled
    unsigned long long int getDegeneracy() const {return prds_set.size();}
    const POccReducedDerivativeStructure& getPRDS() const {return prds_set[0];}
    double getHNFIndex() const {return getPRDS().hnf_index;} //fix me eventually to do robust matrix comparison
    double getEnergy() const {return getPRDS().energy;}
    //const double& getEnergy(){if(getDegeneracy()==0){return AUROSTD_MAX_DOUBLE;} return getPRDS().energy;}
  };

  struct POccDerivativeStructure{
    unsigned long long int hnf_index;
    xmatrix<double> hnf_mat;
    unsigned long long int site_config_index;
    vector<vector<int> > v_types_config;
    double energy;
    unsigned long long int degeneracy;
    const POccDerivativeStructure& operator=(const POccDerivativeStructure& b){
      hnf_index=b.hnf_index;
      hnf_mat=b.hnf_mat;
      site_config_index=b.site_config_index;
      for(uint i=0;i<v_types_config.size();i++){v_types_config[i].clear();} v_types_config.clear(); v_types_config=b.v_types_config;
      energy=b.energy;
      return *this;
    }
    bool operator<(const POccDerivativeStructure& other) const {return energy<other.energy;}
    //void reset(){
    //  hnf_mat.clear();
    //  for(uint i=0;i<v_types_config.size();i++){v_types_config[i].clear();} v_types_config.clear();
    //  energy=0.0;
    //  degeneracy=0;
    //}
  };

  struct UFFParamAtom{
    string symbol;
    double r1;        //bond distance
    double theta0;
    double x1;        //nonbond distance
    double D1;        //nonbond energy
    double zeta;      //scale
    double Z1;        //effective charge
    double Vi;
    double Uj;
    double ChiI;       //electronegativity
    double hard;
    double radius;
  };

  struct UFFParamBond{
    //UFFParamAtom* uffp1=new UFFParamAtom(); //safe
    //UFFParamAtom* uffp2=new UFFParamAtom(); //safe
    //double distij;
    double ren;
    double R0;
    double Kij;
    double Xij;
    double Dij;
    double delta;
    double X6;
    double X12;
    void calculate(UFFParamAtom& uffp1,UFFParamAtom& uffp2,double distij){
      //uffp1=_uffp1; uffp2=_uffp2; distij=_distij;
      //cerr << "SYMBOL " << _uffp1.symbol << " " << uffp1.r1 << "  " << _uffp2.symbol << " " << uffp2.r1 << endl;
      //SKIP equation 3 - zero for single bond (rbo)
      //equation 4
      ren = uffp1.r1 * uffp2.r1 * (pow( sqrt(uffp1.ChiI) - sqrt(uffp2.ChiI),2.0)) / (uffp1.ChiI*uffp1.r1 + uffp2.ChiI*uffp2.r1);
      //equation 2
      //NOTE: See http://towhee.sourceforge.net/forcefields/uff.html: there is a typo in the published paper
      R0 = uffp1.r1 + uffp2.r1 - ren;
      //equation 6
      Kij = 664.12 * uffp1.Z1 * uffp2.Z1 / (R0 * R0 * R0);
      Xij = sqrt(uffp1.x1 * uffp2.x1);
      Dij = sqrt(uffp1.D1 * uffp2.D1);
      delta = distij - R0;
      //cerr << "x1 " << uffp1.x1 << endl;
      //cerr << "x2 " << uffp2.x1 << endl;
      //cerr << "Xij " << Xij << endl;
      //cerr << "distij " << distij << endl;
      //cerr << "Xij/distij " << Xij/distij << endl;
      X6 = pow(Xij/distij,6);
      X12 = pow(Xij/distij,12);
      //cerr << "Dij " << Dij << endl;
      //cerr << "X12 " << X12<< endl;
      //cerr << "X6 " << X6<< endl;
    }
  };

  class POccSiteConfiguration {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccSiteConfiguration();
      POccSiteConfiguration(int _site,int _i_hnf,vector<POccGroup>& _v_pocc_groups);
      POccSiteConfiguration(const POccSiteConfiguration& b);
      //constructors - STOP

      ~POccSiteConfiguration();
      const POccSiteConfiguration& operator=(const POccSiteConfiguration& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END
      
      //debug
      vector<int> types_configuration_debug;                  //atom types, vacancy is -1
      //debug
   
      int site;                                         //reflexive
      int i_hnf;                                        //reflexive
      bool partial_occupation_flag;
      //any vector or xvector is over pocc_groups
      vector<POccGroup> v_pocc_groups;                  //reflexive
      xvector<int> xv_occupation_count_input;      //pre multiplication, from xstr_pocc
      xvector<int> xv_occupation_multiple;              //increment with each POccGroup
      xvector<int> xv_occupation_count_supercell;           //post multiplication with multiple
      xvector<double> xv_partial_occupation_value;
      xvector<double> xv_site_error;
      //sum of occupation_count_total and vacancy_count yields i_hnf (total sites)
      int occupation_count_total;
      int vacancy_count;
      double max_site_error;
      //double error_total;
      
      void prepareNoPOccConfig();
      void preparePOccConfig();
      int getNextOccupationMultiple(int i_hnf,xvector<int>& xv_occupation_multiple);
      int calculateOccupationCountTotal(xvector<int>& xv_next_occupation_multiple);
      void updateOccupationCounts(int _i_hnf, xvector<int> & xv_next_occupation_multiple);
      void calculateError();
      //double getErrorTotal() const;
      bool isPartiallyOccupied() const;
      vector<int> getStartingTypesConfiguration() const;

      //const vector<POccGroup>& getVPOccGroups() const;  //fix me
      //const vector<int>& getTypesConfiguration() const;
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccSiteConfiguration& b);
      //NECESSARY END CLASS METHODS - END

  };

  class POccStructureTemplate {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccStructureTemplate();
      //constructors - STOP
      ~POccStructureTemplate();
      //NECESSARY PUBLIC CLASS METHODS - END
      
      xstructure xstr_pocc;                   //input from PARTCAR
      xvector<double> stoich_each_type;       //converting deque<double> to xvector<double>
      xstructure xstr_nopocc;                 //will contain symmetry objects (_sym_op, pgroups most important here)
      vector<uint> types2pc_map;              //list of atom indices where types2pc_map(0) is 1st type 0 atom, types2pc_map(1) is 1st type 1 atom, etc.
      vector<UFFParamAtom> types2uffparams_map;
      
      void setPOccStructure(const xstructure& xstr_pocc);
      void setNonPOccStructure(const xstructure& xstr_nonpocc);
      void setTypes2PcMap(const vector<uint>& types2pc_map);
      void setTypes2UFFParamsMap(const vector<UFFParamAtom>& types2uffparams_map);
    protected:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccStructureTemplate& b);
      //NECESSARY PRIVATE CLASS METHODS - END
  };

  class POccStructure: public POccStructureTemplate,xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccStructure(ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      POccStructure(ofstream& FileMESSAGE,ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccStructure(const POccStructure& b);
      //POccStructure(xmatrix<double>*& _hnf_mat,xstructure*& _xstr_nopocc);
      //constructors - STOP
      ~POccStructure();
      const POccStructure& operator=(const POccStructure& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END
      
      bool m_initialized;
      _aflags m_aflags;                         //standard aflow flags
      _kflags m_kflags;                         //standard aflow flags
      xstructure xstr_sym;                    //will contain symmetry objects (_sym_op, pgroups most important here)
      int n_hnf;
      vector<POccSite> v_sites;                   //groupings of atoms that are on the same site non-vacant sites only, relative to xstr_nopocc
      int pocc_atoms_total;
      vector<vector<POccGroup> > vv_pocc_groups;
      vector<StructureConfiguration> v_str_configs;
      
      bool initialize(ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      bool initialize(ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);


      void setAFlags(const _aflags& Aflags);                      //standard _aflags
      void setKFlags(const _kflags& Kflags);                      //standard _kflags
      //void setOFStream(ofstream& FileMESSAGE);  //set output stream
      //void setOSS(ostream& oss);                //set output stream
      void preparePOccStructure();

      const xmatrix<double>& getLattice() const;
      const vector<_sym_op>& getPGroup() const;
      const StructureConfiguration& getXStrCountConfiguration(uint i) const;
      //const vector<vector<POccGroup> >& getVVPOccGroups() const;
      
      bool calculateNHNF();                   //get n_hnf
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccStructure& b);
      //NECESSARY PRIVATE CLASS METHODS - END

      //table stuff
      //set some nice printing precisions and paddings, mostly definitions
      uint hnf_table_general_precision;
      uint hnf_table_iteration_padding;
      uint hnf_table_error_padding;
      uint hnf_table_column_padding;
      string header_max_stoich_error;
      string header_max_site_error;
      
      string hnfTableHeader();
      bool getHNFTableOutput(int i_hnf,double& stoich_error,double& site_error);
      string hnfTableLineOutput(int i_hnf,int str_config);
      //void setHNFTablePadding(int _AFLOW_POCC_PRECISION_);
      void partitionPOccSites();              //get v_sites
      xvector<double> calculateStoichEachType(vector<vector<int> >& v_types_config);
      //double calculateStoichDiff(deque<double>& s1,deque<double>& s2);
      //bool updatePOccValues();                //update pocc_values and comp_each_type
      void calculateNonPOccStructure();             //convert pocc xstructure to non-pocc
      void calculateSymNonPOccStructure();          //calculate symmetry of non-pocc structure
      
      //void getSiteCountConfigurations(int i_hnf,double& stoich_error);
      void getSiteCountConfigurations(int i_hnf);
      void getOptimizedSiteCountConfigurations(int site,int i_hnf,vector<POccSiteConfiguration>& v_site_configs);
      bool getUFFParameters(const string& atom,string& params);              //return entry from kesong's messy table
  };

  class POccSuperStructure: public POccStructureTemplate {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccSuperStructure();
      POccSuperStructure(const POccSuperStructure& b);
      //POccSuperStructure(xmatrix<double>*& _hnf_mat,xstructure*& _xstr_nopocc);
      //constructors - STOP
      ~POccSuperStructure();
      const POccSuperStructure& operator=(const POccSuperStructure& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END
      
      //underlying data structures
      //xstructure xstr_bonding;
      xmatrix<double> hnf_mat;
      //POccStructure* p_str;
      double exploration_radius;
      //vector<vector<int> >* v_types_config;
      xmatrix<double> distance_matrix;                    //references xstr_nopocc
      vector<double> v_dist_nn;                           //references xstr_nopocc
      xstructure xstr_ss;       //superstructure
      vector<int> sc2pc_map;
      vector<int> pc2sc_map;
      
      void initialize(POccStructure& p_str,double _exploration_radius);
      void getCluster(xmatrix<double>& _hnf_mat);
      void setBonds(vector<vector<int> >& v_types_config);
      //void replaceTypes(xstructure& supercell, vector<vector<int> >& _v_types_config, vector<uint>& v_vacancies);
      //vector<uint> getVacancies(xstructure& supercell, vector<vector<int> >& v_types_config, bool replace_types=false);
      void replaceRandomSites(xstructure& supercell,const vector<vector<int> >& v_types_config);
      vector<uint> getVacancies(const vector<vector<int> >& v_types_config);
      void copyAtomAttributes(const _atom& atom_qualities,_atom& atom_position);
      xstructure createDerivativeStructure(POccDerivativeStructure& pds,int n_hnf=0,unsigned long long int hnf_count=0,unsigned long long int types_config_permutations_count=0,bool clean_structure=false,bool primitivize=false);
      double getEnergy();
      //void getSpecificBonding();
      bool isVacancy(vector<uint>& v_vacancies,uint atom);
      double getBondEnergy(xstructure& xstr,vector<vector<uint> >& v_bonded_atom_indices,uint MODE);
      //void removeVacancies(xstructure& xstr,vector<uint> v_vacancies);
      void rebuildStructure(xstructure& xstr,vector<uint>& v_vacancies);
      //void removeVacancies(xstructure& xstr,vector<uint>& v_vacancies);
      void removeVacancies(deque<_atom>& atoms,vector<uint>& v_vacancies);
      void rebuildStructure(xstructure& xstr);
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccSuperStructure& b);
      //NECESSARY PRIVATE CLASS METHODS - END

      //for bonding, we need to create a super-superstructure (cluster)
      //deque<_atom> cluster;                                   //cluster of atoms within radius
      //vector<int> ssc2sc_map;
      //vector<int> sc2ssc_map;
      xstructure xstr_cluster;                                //cluster of atoms within radius
      vector<vector<uint> > v_bonded_atom_indices;            //references xstr_cluster
      vector<vector<uint> > v_nonbonded_atom_indices;         //references xstr_cluster

      bool has_vacancies;                                     //are there vacancies present in v_types_config_bonding?
      bool bonding_set;                                       //have we already found bonding for this configuration?
      vector<vector<int> > v_types_config_bonding;            //the config for which we determined bonding
      vector<uint> v_vacancies_bonding;

      double energy;

      //bool hasVacancies(vector<vector<int> >& v_types_config);
      //bool hasVacancies();
      uint NNDistancesMapPC(uint atom);
      uint NNDistancesMapSC(uint atom);
      void calculateNNDistances(xstructure& xstr,vector<uint>& v_vacancies);
  };

  /*class POccDerivativeStructure {
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccDerivativeStructure& b);
      //NECESSARY PRIVATE CLASS METHODS - END

      //underlying data structures
      uint hnf_index;                       //refers to v_hnf
      vector<vector<int> > v_occupants;     //first index refers to v_pocc_sites, second index is site in supercell, value is occupant and refers to xstr_pocc - NOTE that it's an int and not uint, if negative, it's a vacancy
      double energy;                        //uff
      unsigned long long int degerancy;     //number of duplicate structures this structure represents
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccDerivativeStructure();
      POccDerivativeStructure(const POccDerivativeStructure& b);
      //constructors - STOP
      ~POccDerivativeStructure();
      const POccDerivativeStructure& operator=(const POccDerivativeStructure& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END
  };*/

  class POccCalculator : public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccCalculator(ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      POccCalculator(ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const POccCalculator& b);
      //constructors - STOP
      ~POccCalculator();
      const POccCalculator& operator=(const POccCalculator& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END

      //keep going, need to decide what's public/private
      //then figure out poccstructure/pocccalculator

      //inputs
      bool m_initialized;
      _aflags m_aflags;                         //standard aflow flags
      _kflags m_kflags;                         //standard aflow flags
      POccStructure p_str;
      POccSuperStructure ps_str;
      unsigned long long int hnf_count;
      unsigned long long int types_config_permutations_count;
      unsigned long long int total_permutations_count;
      std::list<POccReducedDerivativeStructuresSet> l_derivative_structures_sets;
      //std::list<POccDerivativeStructure> l_derivative_structures;
      //vector<xmatrix<double> > v_hnf;         //vector of unique hnf's (lattices for supercell)
      //vector<POccDerivativeStructure> v_sc;   //vector of supercells (bare info - reduction of memory)
      //standard flags - ALL options will be handled via xoptions
      aurostd::xoption enumerator_mode;       //how do we determine duplicates - UFF, SNF, ...
      
      bool initialize(ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      bool initialize(ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);

      //external methods
      void setPOccStructure(const xstructure& xstr_pocc);
      void setAFlags(const _aflags& Aflags);                      //standard _aflags
      void setKFlags(const _kflags& Kflags);                      //standard _kflags
      //void CopyStreams(ofstream& ofs1,ofstream& ofs2);      //copy two streams
      //void CopyStreams(ostream& os1,ostream& os2);          //copy two streams
      //void setOFStream(ofstream& _FileMESSAGE);             //set output stream
      //void setOSS(ostream& _oss);                           //set output stream
      //double getBondEnergy(xstructure& xstr,vector<uint>& v_vacancies,xmatrix<double>& distance_matrix,vector<vector<uint> >& v_bonded_atom_indices,uint MODE=BOND_MODE);
      bool areEquivalentStructuresByUFF(std::list<POccReducedDerivativeStructuresSet>::iterator it, const POccReducedDerivativeStructure& prds) const;
      void add2DerivativeStructuresList(const POccReducedDerivativeStructure& prds,std::list<POccReducedDerivativeStructuresSet>::iterator i_start,std::list<POccReducedDerivativeStructuresSet>::iterator i_end);
      void add2DerivativeStructuresList(const POccReducedDerivativeStructure& prds);
      //void add2DerivativeStructuresList(POccDerivativeStructure& pds,std::list<POccDerivativeStructure>::iterator i_start,std::list<POccDerivativeStructure>::iterator i_end);
      //void add2DerivativeStructuresList(POccDerivativeStructure& pds);
      POccDerivativeStructure getParticularPOccDerivativeStructure(const POccReducedDerivativeStructure& prds);
      bool areEquivalentStructures(const POccReducedDerivativeStructure& prds_a,const POccReducedDerivativeStructure& prds_b);
      bool areEquivalentStructures(const xstructure& a,const xstructure& b);
      unsigned long long int runRobustStructureComparison(std::list<POccReducedDerivativeStructuresSet>::iterator it);
      bool calculate();
      string getARUNString(unsigned long long int i);
      xstructure getUniqueDerivativeStructure(unsigned long long int i);
      vector<xstructure> getUniqueDerivativeStructures();
      unsigned long long int getUniqueDerivativeStructuresCount() const;
      //bool printUniqueDerviativeStructures();
      //void resetMaxSLRadius();
      void resetHNFMatrices();
      void resetSiteConfigurations();
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccCalculator& b);
      //NECESSARY END CLASS METHODS - END
      
      //hnf matrices
      int a_start, b_start, c_start;
      int d_start, e_start, f_start;
      vector<xmatrix<double> > v_unique_superlattices;
      //double max_superlattice_radius;
      //configurations
      xmatrix<double> hnf_mat;
      vector<vector<int> > v_types_config;
      //vector<int> v_config_iterators;
      uint config_iterator;
      vector<uint> v_config_order;

      //useful internal methods
      void getTotalPermutationsCount();
      //vector<_sym_op> getPGroups();           //fetch pgroups of xstr_nopocc
      //uint getHNFCount();                     //get count of hnf matrices, refers to v_hnf
      //xmatrix<double> getHNFMatrix(uint i);   //fetch specific hnf matrix, refers to v_hnf
      //_atom getAtom(uint i);                  //grab specific atom, refers to xstr_pocc
      
      //determines site occupancy and vacancy count, given n_hnf
      bool getNextHNFMatrix();                //calculate all unique hnf's
      void getConfigOrder();
      bool getNextSiteConfiguration();
      bool getNextSiteConfiguration(vector<vector<int> >& v_site_config);
      double getEnergy();
      //void RemoveAtoms(xstructure& xstr,vector<uint>& v_atoms_to_remove);
      //void getBonding(xmatrix<double>& hnf_mat,xmatrix<double>& _distance_matrix,vector<vector<uint> >& v_bonded_atom_indices,vector<vector<uint> >& v_nonbonded_atom_indices);
  };

}

#endif  // _AFLOW_POCC_H_

// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2018             *
// *                 Kesong Yang - Duke University 2010-2011                 *
// *                  Corey Oses - Duke University 2013-2018                 *
// *                                                                         *
// ***************************************************************************
