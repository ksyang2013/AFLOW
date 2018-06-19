// ***************************************************************************
// *                                                                         *
// *              AFlow STEFANO CURTAROLO  Duke University 2003-2017         *
// *              AFlow COREY OSES  Duke University 2013-2017                *
// *                                                                         *
// ***************************************************************************
// aflow_pocc.cpp
// completely revised approach to KESONG YANG's original implementation
// no recursion needed --- too memory intensive
// issues with UFF (structure comparison), multiply occupied sites, and 
// vacancies are all addressed robustly here
// 2013-2017: corey.oses@duke.edu
// 2010-2011: kesong.yang@gmail.com (LEGACY)

// This file contains the routines to prepare partial occupation input files.

#ifndef _AFLOW_POCC_CPP_
#define _AFLOW_POCC_CPP_

#include "aflow.h"
#include "aflow_pocc.h"
#include "aflow_compare_structure.h"

//tolerances
const double ZERO_TOL = pow(10,-PRECISION);

//make defaults in AFLOW_RC
const double BONDING_THRESHOLD = 0.5; //angstroms
const double ENERGY_TOLERANCE = 1e-6; //eV
const double ENERGY_RADIUS = 10; //angstroms
//const bool STRICT_STOICH_EQUIVALENCE = true;
const bool FIX_POSCAR_TITLES = false;
const bool PERFORM_ROBUST_STRUCTURE_COMPARISON = false;
const bool WRITE_OUT_ALL_STRUCTURES = false;
//const int CLUSTER_EXPLORATION_MULTIPLE = 1;
//const bool REDUCE_VACANCIES = true;

//some constants
const int BAR_WIDTH = 70;
const int A_START = 1, C_START = 1, F_START = 1;
const int B_START = 0, D_START = 0, E_START = 0;

//conversion for uff params
const double KCAL_2_EV  = 4.336443203200000E-002; // 1(kcal/mol) = 4.33644E-2 eV

bool PRIMITIVIZE=false;

//testing
bool PRINT_STRUCT_COMP_PRANAB=false;
bool COMPARE_WITH_KESONG=false;
bool SET_KESONG_STANDARD_DIST=(true||COMPARE_WITH_KESONG);
bool ENUMERATE_ALL_HNF=false;

namespace pocc {
  //temporary function to replace kesong's code, but should never be called
  //does NOT handle aflow.in generation well at all
  bool poccInput() {
    string soliloquy="pocc::poccInput():";
    stringstream message;
    ostream& oss=cout;
    ofstream FileMESSAGE; FileMESSAGE.open("aflow_pocc.log");
    _aflags aflags;aflags.Directory="./";
    
    string AflowIn_file=string(aflags.Directory+"/"+_AFLOWIN_);
    if(!aurostd::FileExist(AflowIn_file)) {
      message << "Input file does not exist: " << AflowIn_file;
      pflow::logger(soliloquy,message,FileMESSAGE,oss,_LOGGER_ERROR_);
      return false;
    }
    message << "Using input file: " << AflowIn_file;
    pflow::logger(soliloquy,message,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    string AflowIn;
    aurostd::file2string(AflowIn_file,AflowIn);
    AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
    
    try{generateStructures(AflowIn,FileMESSAGE,oss);}
    catch(AFLOWRuntimeError& re){
      pflow::logger(re.where(), re.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      FileMESSAGE.close();
      return false;
    }
    catch(AFLOWLogicError& le){
      pflow::logger(le.where(), le.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      FileMESSAGE.close();
      return false;
    }
    FileMESSAGE.close();
    return true;
  }
} // namespace pocc

namespace KBIN {
	void VASP_RunPOCC(const _xvasp& xvasp,const string& AflowIn,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE) {
    string soliloquy="KBIN::VASP_RunPOCC():";
    stringstream message;
    ostream& oss=cout;
    
    //at some point, program it in that you should check if UNIQUE_DERIVATIVE_STRUCTURES_FILE exists
    //if so, don't calculate, grab structures from there
    if(!pocc::structuresGenerated()){
      try{pocc::generateStructures(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);}
      catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);}
      catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);}
      return;
    }

    //continue with post-processing here
    //do try/catch and rerun generateStructures()
  
  }
} // namespace pocc

namespace pocc {
  bool structuresGenerated(){return aurostd::FileNotEmpty(UNIQUE_DERIVATIVE_STRUCTURES_FILE);}
  xstructure extractPARTCAR(const string& AflowIn){
    stringstream ss_pocc_structure;
    aurostd::ExtractToStringstreamEXPLICIT(AflowIn,ss_pocc_structure, "[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE", "[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE");
    xstructure xstr_pocc(ss_pocc_structure, IOVASP_POSCAR);
    return xstr_pocc;
  }
  void generateStructures(const string& AflowIn, ostream& oss) {ofstream FileMESSAGE;return generateStructures(AflowIn,FileMESSAGE,oss);}
  void generateStructures(const string& AflowIn, ofstream& FileMESSAGE, ostream& oss) {
    _xvasp xvasp;_aflags aflags;_kflags kflags;
    aflags.Directory="./";
    _vflags vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags);
    return generateStructures(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void generateStructures(const _xvasp& in_xvasp,const string& AflowIn,const _aflags& in_aflags,const _kflags& in_kflags,const _vflags& in_vflags, ofstream& FileMESSAGE, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="pocc::generateStructures():";
    stringstream message;
    _xvasp xvasp=in_xvasp;_aflags aflags=in_aflags;_kflags kflags=in_kflags;_vflags vflags=in_vflags; //make copies

    xstructure xstr_pocc=pocc::extractPARTCAR(AflowIn);
    pocc::POccCalculator pcalc(xstr_pocc,aflags,kflags,FileMESSAGE,oss);
    if(!pcalc.m_initialized){throw AFLOWLogicError(soliloquy,"POccCalculator failed to initialized");}
    if(!pcalc.calculate()){throw AFLOWLogicError(soliloquy,"POccCalculator failed to calculate");}
    
    message << "Creating list of unique derivative supercells.";
    pflow::logger(soliloquy,message,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    
    stringstream unique_derivative_structures_ss;
    stringstream new_AflowIn_ss;
    bool original_quiet_flag=XHOST.QUIET;


    for(unsigned long long int i=0;i<pcalc.getUniqueDerivativeStructuresCount();i++) {

      //populate UNIQUE_DERIVATIVE_STRUCTURES_FILE
      unique_derivative_structures_ss << AFLOWIN_SEPARATION_LINE<< endl;
      unique_derivative_structures_ss << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl; // ." << ss_pocc_count.str() << endl;
      unique_derivative_structures_ss << pcalc.getUniqueDerivativeStructure(i);
      unique_derivative_structures_ss << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl; // ." << ss_pocc_count.str() << endl;
      unique_derivative_structures_ss << AFLOWIN_SEPARATION_LINE<< endl;
    
      string arun_directory=pcalc.getARUNString(i);

      xvasp.str.Clear();
      xvasp.str=pcalc.getUniqueDerivativeStructure(i);
      xvasp.AVASP_prototype_from_library_=false;
      xvasp.AVASP_prototype_mode=LIBRARY_MODE_XSTRUCTURE;
      xvasp.Directory=arun_directory;
      
      XHOST.QUIET=true;
      if(!KBIN::VASP_Produce_and_Modify_INPUT(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags,true)){  //load from xvasp.str
        throw AFLOWLogicError(soliloquy,"Cannot produce/modify input");
      }
      if(LDEBUG){
        cerr << soliloquy << " printing xstructure" << endl;
        cerr << xvasp.str << endl;
        cerr << soliloquy << " printing POSCAR" << endl;
        cerr << xvasp.POSCAR.str() << endl;
        cerr << soliloquy << " xvasp.str.species.size()=" << xvasp.str.species.size() << endl;
      }
      AVASP_MakeSingleAFLOWIN(xvasp,new_AflowIn_ss,false,-1,false);  //don't write/print and hence don't pthread
      XHOST.QUIET=original_quiet_flag;

      message << "Creating " << arun_directory;
      pflow::logger(soliloquy,message,FileMESSAGE,oss,_LOGGER_MESSAGE_);

      if(!aurostd::FileExist(arun_directory)){
        aurostd::DirectoryMake(arun_directory);
        aurostd::stringstream2file(new_AflowIn_ss,string(arun_directory+"/"+_AFLOWIN_));
      }else{
        message << arun_directory << " already exists, skipping!";
        pflow::logger(soliloquy,message,FileMESSAGE,oss,_LOGGER_WARNING_);
      }
      //updateProgressBar(i,pcalc.getUniqueDerivativeStructuresCount()-1,oss);
    }

    aurostd::stringstream2file(unique_derivative_structures_ss,UNIQUE_DERIVATIVE_STRUCTURES_FILE);
  }
} // namespace pocc

namespace pocc {
  vector<string> getElementsList() {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="pocc::getElementsList():";
    uint skip_every=2;
    vector<string> elements;
    vector<string> VSTD_ELEMENTS_LIST;
    aurostd::string2tokens(STD_ELEMENTS_LIST,VSTD_ELEMENTS_LIST," ");
    for(uint i=0;i<skip_every+1;i++){
      for(uint j=i;j<VSTD_ELEMENTS_LIST.size();j+=skip_every+1){
        elements.push_back(VSTD_ELEMENTS_LIST[j]);
      }
    }
    if(LDEBUG){
      cerr << soliloquy << " ";
      for(uint i=0;i<elements.size();i++){
        cerr << elements[i] << " ";
      }
      cerr << endl;
    }
    return elements;
  }
} // namespace pocc

namespace pocc {
//--------------------------------------------------------------------------------
// class POccCalculator
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// constructor
//--------------------------------------------------------------------------------
POccCalculator::POccCalculator(ostream& oss) : xStream() {initialize(oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,ostream& oss) : xStream() {initialize(xstr_pocc,oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) : xStream() {initialize(xstr_pocc,kflags,oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,kflags,oss);}
POccCalculator::POccCalculator(ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(FileMESSAGE,oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,FileMESSAGE,oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,FileMESSAGE,oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,kflags,FileMESSAGE,oss);}
POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,kflags,FileMESSAGE,oss);}
POccCalculator::POccCalculator(const POccCalculator& b) {copy(b);} // copy PUBLIC

POccCalculator::~POccCalculator() {freeAll();}

const POccCalculator& POccCalculator::operator=(const POccCalculator& b) {  // operator= PUBLIC
  if(this!=&b) {freeAll();copy(b);}
  return *this;
}

void POccCalculator::clear() {POccCalculator a;copy(a);}
void POccCalculator::free() {
  m_initialized=false;
  m_aflags.clear();
  m_kflags.clear();
  p_str.clear();
  ps_str.clear();
  resetHNFMatrices();
  resetSiteConfigurations();
  hnf_count=0;
  types_config_permutations_count=0;
  total_permutations_count=0;
  //for(std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();it!=l_derivative_structures_sets.end();++it){(*it).clear();}
  l_derivative_structures_sets.clear();
  enumerator_mode.clear();
}

void POccCalculator::copy(const POccCalculator& b) { // copy PRIVATE
  xStream::copyStreams(b);
  m_initialized=b.m_initialized;
  m_aflags=b.m_aflags;
  m_kflags=b.m_kflags;
  p_str=b.p_str;
  ps_str=b.ps_str;
  hnf_mat=b.hnf_mat;
  a_start=b.a_start;b_start=b.b_start;c_start=b.c_start;
  d_start=b.d_start;e_start=b.e_start;f_start=b.f_start;
  v_unique_superlattices=b.v_unique_superlattices;
  v_types_config=b.v_types_config;
  config_iterator=b.config_iterator;
  v_config_order=b.v_config_order;
  hnf_count=b.hnf_count;
  types_config_permutations_count=b.types_config_permutations_count;
  total_permutations_count=b.total_permutations_count;
  //for(std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();it!=l_derivative_structures_sets.end();++it){(*it).clear();}
  l_derivative_structures_sets.clear();
  for(std::list<POccReducedDerivativeStructuresSet>::const_iterator it=b.l_derivative_structures_sets.begin();it!=b.l_derivative_structures_sets.end();++it){l_derivative_structures_sets.push_back(*it);}
  enumerator_mode=b.enumerator_mode;
}

bool POccCalculator::initialize(ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,aflags,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,kflags,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,aflags,kflags,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccCalculator::initialize(ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    m_initialized=false;  //no point
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setAFlags(aflags);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setKFlags(kflags);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setAFlags(aflags);
    setKFlags(kflags);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}

void POccCalculator::setPOccStructure(const xstructure& xstr_pocc) {
  p_str.setPOccStructure(xstr_pocc);
  p_str.preparePOccStructure();
}

void POccCalculator::setAFlags(const _aflags& Aflags) {
  m_aflags=Aflags;
  p_str.setAFlags(m_aflags);
}

void POccCalculator::setKFlags(const _kflags& Kflags) {
  m_kflags=Kflags;
  p_str.setKFlags(m_kflags);
}

//bool POccCalculator::updatePOccValues() {
//  //This function is only used to update partial occupation value i.e., comp_each_type
//  xstr_pocc.comp_each_type.clear();
//  for(uint i=0;i<xstr_pocc.num_each_type.size();i++) {xstr_pocc.comp_each_type.push_back(0.0);}
//
//  for(uint i=0;i<xstr_pocc.atoms.size();i++) {
//    xstr_pocc.atoms[i].partial_occupation_flag = !aurostd::isequal(xstr_pocc.atoms[i].partial_occupation_value,1.0,ZERO_TOL);
//    xstr_pocc.atoms[i].partial_occupation_value = ( xstr_pocc.atoms[i].partial_occupation_flag ? xstr_pocc.atoms[i].partial_occupation_value : 1.0 );
//    xstr_pocc.comp_each_type[xstr_pocc.atoms[i].type] += xstr_pocc.atoms[i].partial_occupation_value;
//  }
//
//  //corey add check for if all sites are fully occupied
//  return true;
//}

//void POccCalculator::RemoveAtoms(xstructure& xstr,vector<uint>& v_atoms_to_remove){
//  std::sort(v_atoms_to_remove.rbegin(),v_atoms_to_remove.rend()); //NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
//  for(uint atom=0;atom<v_atoms_to_remove.size();atom++){xstr.RemoveAtom(v_atoms_to_remove[atom]);}
//}

//void POccCalculator::resetMaxSLRadius() {
//  max_superlattice_radius=0.0;
//}

void POccCalculator::resetHNFMatrices() {
  hnf_mat.clear();
  //a_start=1;c_start=1;f_start=1;
  //b_start=0;d_start=0;e_start=0;
  a_start=A_START;c_start=C_START;f_start=F_START;
  b_start=B_START;d_start=D_START;e_start=E_START;
  v_unique_superlattices.clear();
}

bool POccCalculator::getNextHNFMatrix(){
  bool LDEBUG=(FALSE || ENUMERATE_ALL_HNF || XHOST.DEBUG);
  string soliloquy="POccCalculator::getNextHNFMatrix()";
  
  xmatrix<double> _hnf_mat(3,3),duplicate_mat(3,3);
  xmatrix<double> superlattice(3,3), rotated_superlattice(3,3);
  bool duplicate_H;
  int n_hnf=p_str.n_hnf;
  const xmatrix<double>& lattice=p_str.getLattice();
  const vector<_sym_op>& pgroup=p_str.getPGroup();

  //this allows us to reset completely everytime (SAFE but inefficient)
  //int a_start=1;int c_start=1;int f_start=1;
  //int b_start=0;int d_start=0;int e_start=0;

  // PHYSICAL REVIEW B 77, 224115 2008
  // Algorithm for generating derivative structureslgorithm for generating derivative structures
  // Gus L. W. Hart and Rodney W. Forcade
  // IMPORTANT DISTINCTION FROM PAPER, our lattice is in row matrix, not columns
  // let's work in row space to avoid so many transpose operations
  // also, slight optimization for both speed and memory:
  // we leverage a scheme that remembers where it was last, and picks up here in the nested for-loop
  // it relies on saving a to a_start, b to b_start, etc. when it finds a unique hnf_mat
  // then we start the loops at these values only for the initial loop through
  // starting_config ensures that we reset the loops back to A_START, B_START, etc. after the initial loop
  // that way, we explore all possibilities efficiently without having to save the hnf's (only corresponding superlattices)
  // this iterative scheme helps save some memory, which we need elsewhere in the code
  // starting_config ensures we don't ALSO lose calculation iterations in recalculating hnf possibilities we already explored
  // REMEMBER when comparing with kesong's hnf's, his are TRANSPOSED!
	bool starting_config;
  for(int a=a_start;a<=n_hnf;a++){
		starting_config=(a==a_start);
    for(int c=( starting_config ? c_start : C_START );c<=n_hnf/a;c++){ //for(int c=( a==a_start ? c_start : C_START );c<=n_hnf/a;c++){ //for(int c=c_start;c<=n_hnf/a;c++){
			starting_config&=(c==c_start);
      for(int f=( starting_config ? f_start : F_START );f<=n_hnf/a/c;f++){ //for(int f=( (a==a_start && c==c_start) ? f_start : F_START );f<=n_hnf/a/c;f++){ //for(int f=f_start;f<=n_hnf/a/c;f++){
        if(a*c*f==n_hnf){
          //found new viable diagonal set, now enumerate based on off-diagonals
					starting_config&=(f==f_start);
          for(int b=( starting_config ? b_start : B_START );b<c;b++){ //for(int b=( (a==a_start && c==c_start && f==f_start) ? b_start : B_START );b<c;b++){ //for(int b=b_start;b<c;b++){
						starting_config&=(b==b_start);
            for(int d=( starting_config ? d_start : D_START );d<f;d++){ //for(int d=( (a==a_start && c==c_start && f==f_start && b==b_start) ? d_start : D_START );d<f;d++){ //for(int d=d_start;d<f;d++){
							starting_config&=(d==d_start);
              for(int e=( starting_config ? e_start : E_START );e<f;e++){ //for(int e=( (a==a_start && c==c_start && f==f_start && b==b_start && d==d_start) ? e_start : E_START );e<f;e++){ //for(int e=e_start;e<f;e++){
                if(LDEBUG){
								  cerr << " a=" << a << "(max=n_hnf=" << n_hnf << ")";
								  cerr << " c=" << c << "(max=n_hnf/a=" << n_hnf/a << ")";
								  cerr << " f=" << f << "(max=n_hnf/a/c=" << n_hnf/a/c << ")";
								  cerr << " b=" << b << "(max=c=" << c-1 << ")";
								  cerr << " d=" << d << "(max=f=" << f-1 << ")";
								  cerr << " e=" << e << "(max=f=" << f-1 << ")";
								  cerr << endl;
                }
                
                _hnf_mat(1,1)=a;_hnf_mat(1,2)=b;_hnf_mat(1,3)=d;
                _hnf_mat(2,1)=0;_hnf_mat(2,2)=c;_hnf_mat(2,3)=e;
                _hnf_mat(3,1)=0;_hnf_mat(3,2)=0;_hnf_mat(3,3)=f;
                
                if(LDEBUG){
                  cerr << "HNF MAT " << endl;
                  cerr << _hnf_mat << endl;
                }

                superlattice=_hnf_mat*lattice;

                if(LDEBUG){
                  cerr << "SUPERLATTICE " << endl;
                  cerr << superlattice << endl;
                }

                duplicate_H=false;
                for(uint i=0;i<v_unique_superlattices.size() && !duplicate_H;i++){
                  if(LDEBUG){
                    cerr << "SUPERLATTICES OLD " << endl;
                    cerr << v_unique_superlattices[i] << endl;
                  }
                  for(uint pg=0;pg<pgroup.size() && !duplicate_H;pg++){
                    if(LDEBUG){
                      cerr << "POINT GROUP " << endl;
                      cerr << pgroup[pg].Uc << endl;
                    }
                    //original matrix in column space: Bi = R * Bj * H, Gus derives H = Bj^-1 * R^-1 * Bi
                    //in row space:  
                    //(Bi)^T = (R * Bj * H)^T
                    //Bi^T = H^T * Bj^T * R^T
                    //Bi^T * (R^T)^-1 * (Bj^T)^-1 = H^T
                    duplicate_mat = (superlattice * aurostd::inverse(pgroup[pg].Uc) * aurostd::inverse(v_unique_superlattices[i]));
                    duplicate_H = aurostd::isinteger(duplicate_mat);
                    if(LDEBUG){
                      cerr << "DUPLICATE MATRIX" << endl;
                      cerr << duplicate_mat << endl;
                    }
                    //duplicate_H=false;
                  }
                }
                //cerr << "n_hnf " << n_hnf << endl;
                //cerr << "lattice " << endl;
                //cerr << lattice << endl;
                //cerr << "hnf_mat" << endl;
                //cerr << _hnf_mat << endl;
                //cerr << "superlattice " << endl;
                //cerr << superlattice << endl;
                ////duplicate_H=true;
                if(ENUMERATE_ALL_HNF){
                  if(!duplicate_H){ //still get FOUND output
                    cerr << "----------------------" << endl;
                    cerr << "FOUND " << endl;
                    cerr << _hnf_mat << endl;
                    cerr << "----------------------" << endl;
                    duplicate_H=true; //so we go through all iterations
                  }
                }
                if(!duplicate_H){
                  if(LDEBUG){
                    cerr << "----------------------" << endl;
                    cerr << "FOUND " << endl;
                    cerr << _hnf_mat << endl;
                    cerr << "----------------------" << endl;
                  }
                  //found new matrix, set new starting conditions and return
                  a_start=a;b_start=b;c_start=c;
                  d_start=d;e_start=e;f_start=f;
                  hnf_mat=_hnf_mat;
                  v_unique_superlattices.push_back(superlattice);
                  //cerr << hnf_mat << endl;
                  //cerr << endl;
                  //max_superlattice_radius=max(max_superlattice_radius,RadiusSphereLattice(superlattice));

                //cerr << "n_hnf " << n_hnf << endl;
                //cerr << "lattice " << endl;
                //cerr << lattice << endl;
                //cerr << "hnf_mat" << endl;
                //cerr << _hnf_mat << endl;
                //cerr << "superlattice " << endl;
                //cerr << superlattice << endl;
                //cerr << "max_radius " << max_superlattice_radius << endl;
                //cerr << "SEE " << max_superlattice_radius << endl;

									//cerr << endl;
                  return true;
                //}else{
                //  cerr << "DUPLICATE " << endl;
                //cerr << "n_hnf " << n_hnf << endl;
                //cerr << "lattice " << endl;
                //cerr << lattice << endl;
                //cerr << "hnf_mat" << endl;
                //cerr << _hnf_mat << endl;
                //cerr << "superlattice " << endl;
                //cerr << superlattice << endl;
                  
                }
              }
            }
          }
        }
      }
    }
  }

  return false;
}

void POccCalculator::resetSiteConfigurations() {
  v_types_config.clear();
  config_iterator=0;
  //v_config_iterators.clear();
  v_config_order.clear();
}

void POccCalculator::getConfigOrder(){
  void getConfigOrder(vector<vector<int> >& v_types_config,vector<uint>& v_config_order);
  return getConfigOrder(v_types_config,v_config_order);
}

void getConfigOrder(vector<vector<int> >& v_types_config,vector<uint>& v_config_order){
  vector<SiteVacancyCount> v_svc;
  for(uint site=0;site<v_types_config.size();site++){
    v_svc.push_back(SiteVacancyCount());
    v_svc.back().site=site;
    v_svc.back().vacancy_count=0;
    for(uint pos=0;pos<v_types_config[site].size();pos++){
      if(v_types_config[site][pos]==-1){v_svc.back().vacancy_count++;}
    }
  }
  std::sort(v_svc.begin(),v_svc.end()); //most vacancies should go last
  v_config_order.clear();
  for(uint i=0;i<v_svc.size();i++){v_config_order.push_back(v_svc[i].site);}
}

bool POccCalculator::getNextSiteConfiguration() {
  //const vector<vector<POccSiteConfiguration> >& vv_count_configs=p_str.v_str_configs;
  //bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,
  //    vector<uint>& v_config_order,
  //    vector<int>& v_config_iterators,
  //    vector<vector<int> >& v_types_config);
  //return getNextSiteConfiguration(p_str.v_str_configs,v_count_order,v_config_iterators,v_types_config);
  if(v_types_config.size()==0){
    config_iterator=0;
    //v_config_iterators.clear();
    const vector<StructureConfiguration>& v_str_configs=p_str.v_str_configs;
    //cerr << "SEE " << v_str_configs.size() << endl;
    for(uint site=0;site<v_str_configs[config_iterator].site_configs.size();site++){
      //v_config_iterators.push_back(0);
      v_types_config.push_back(v_str_configs[config_iterator].site_configs[site].getStartingTypesConfiguration());
    }
    getConfigOrder();
    //cerr << "WOW2 " << v_types_config.size() << endl;
    //exit(0);
    return true;
  }
  bool getNextSiteConfiguration(vector<uint>& v_config_order,vector<vector<int> >& v_types_config);
  if(!getNextSiteConfiguration(v_config_order,v_types_config)){
    if(config_iterator>=p_str.v_str_configs.size()-1){return false;}
    config_iterator++;
    const vector<StructureConfiguration>& v_str_configs=p_str.v_str_configs;
    for(uint site=0;site<v_types_config.size();site++){v_types_config[site].clear();} v_types_config.clear();
    for(uint site=0;site<v_str_configs[config_iterator].site_configs.size();site++){
      //v_config_iterators.push_back(0);
      v_types_config.push_back(v_str_configs[config_iterator].site_configs[site].getStartingTypesConfiguration());
    }
    getConfigOrder();
    return true;
  }
  return true;
}

bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,
    vector<int>& v_config_iterators, 
    vector<vector<int> >& v_types_config) {
  if(v_types_config.size()==0){
    v_config_iterators.clear();
    //const vector<vector<POccSiteConfiguration> >& vv_count_configs=p_str.v_str_configs;
    //cerr << "SEE " << vv_count_configs.size() << endl;
    for(uint site=0;site<vv_count_configs.size();site++){
      v_config_iterators.push_back(0);
      v_types_config.push_back(vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration());
    }
    //getConfigOrder(v_types_config,v_config_order);
    //cerr << "WOW2 " << v_types_config.size() << endl;
    //exit(0);
    return true;
  }
  for(uint site=0;site<v_config_iterators.size();site++){
    if(v_config_iterators[site]>=(int)vv_count_configs[site].size()-1){  //greater than not necessary but safe
      if(site==v_config_iterators.size()-1){break;} //stop condition
      v_config_iterators[site]=0; 
      v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();
      continue;
    }
    v_config_iterators[site]++;
    v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();
    return true;
  }
  return false;
}

//bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,
//    vector<uint>& v_config_order,
//    vector<int>& v_config_iterators, 
//    vector<vector<int> >& v_types_config) {
//  //starting condition START - v_types_config is empty!
//  if(v_types_config.size()==0){
//    v_config_iterators.clear();
//    //const vector<vector<POccSiteConfiguration> >& vv_count_configs=p_str.v_str_configs;
//    //cerr << "SEE " << vv_count_configs.size() << endl;
//    for(uint site=0;site<vv_count_configs.size();site++){
//      v_config_iterators.push_back(0);
//      v_types_config.push_back(vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration());
//    }
//    getConfigOrder(v_types_config,v_config_order);
//    //cerr << "WOW2 " << v_types_config.size() << endl;
//    //exit(0);
//    return true;
//  }
//  //starting condition STOP - v_types_config is empty!
//  //otherwise, similar next bitstring generator
//  if(!getNextSiteConfiguration(v_config_order,v_types_config)){
//    //const vector<vector<POccSiteConfiguration> >& vv_count_configs=p_str.v_str_configs;
//    //if(1){
//    //for(uint site=0;site<vv_count_configs.size();site++){
//    //  if(v_config_iterators[site]>=(int)vv_count_configs[site].size()-1){
//    //    if(site==vv_count_configs.size()-1){break;} //stop condition
//    //    v_config_iterators[site]=0;
//    //    v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
//    //    continue;
//    //  }
//    //  v_config_iterators[site]++;
//    //  v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
//    //  getConfigOrder(v_types_config,v_config_order);
//    //  return true;
//    //}
//    //}
//    //if(0){
//    uint site;
//    for(uint index=0;index<v_config_order.size();index++){
//      site=v_config_order[index];
//      //cerr << "index = " << index << " site=" << site << " v_config_iterators[site]=" << v_config_iterators[site] << endl;
//      if(v_config_iterators[site]>=(int)vv_count_configs[site].size()-1){
//        //if(site==vv_count_configs.size()-1){break;} //stop condition
//        if(index==vv_count_configs.size()-1){break;} //stop condition
//        v_config_iterators[site]=0;
//        v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
//        continue;
//      }
//      //cerr << "HERE" << endl;
//      v_config_iterators[site]++;
//      v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
//      getConfigOrder(v_types_config,v_config_order);
//      return true;
//    }
//    //}
//    return false; //stop condition SAFE
//  }
//  return true;
//}

bool POccCalculator::getNextSiteConfiguration(vector<vector<int> >& v_types_config) {
  bool getNextSiteConfiguration(vector<uint>& v_config_order,vector<vector<int> >& v_types_config);
  return getNextSiteConfiguration(v_config_order,v_types_config);
}

bool getNextSiteConfiguration(vector<uint>& v_config_order,
    vector<vector<int> >& v_types_config) {
  //bool getNextSiteConfiguration(vector<int>& types_config);
  //if(1){
  //for(uint site=0;site<v_types_config.size();site++){
  //  if(!getNextSiteConfiguration(v_types_config[site])){
  //    //reset site first
  //    //we know that the last permutation is the reverse sorted of original vector
  //    //we can either a) completely sort or b) simply reverse (slightly faster)
  //    //std::sort(v_types_config[site].begin(),v_types_config[site].end());
  //    std::reverse(v_types_config[site].begin(),v_types_config[site].end());
  //    if(site==v_types_config.size()-1){return false;} //stop condition
  //    continue;
  //  }
  //  return true;
  //}
  //return false; //stop condition SAFE
  //}
  //here we optimized ordering to avoid fewest number of bonding recalculations
  //if(0){
  uint site;
  for(uint index=0;index<v_config_order.size();index++){
    site=v_config_order[index];
    if(!getNextSiteConfiguration(v_types_config[site])){
      //reset site first
      //we know that the last permutation is the reverse sorted of original vector
      //we can either a) completely sort or b) simply reverse (slightly faster)
      //std::sort(v_types_config[site].begin(),v_types_config[site].end());
      std::reverse(v_types_config[site].begin(),v_types_config[site].end());
      //if(site==v_types_config.size()-1){return false;} //stop condition
      if(index==v_types_config.size()-1){return false;} //stop condition
      continue;
    }
    return true;
  }
  //}
  return false; //stop condition SAFE
}

bool getNextSiteConfiguration(vector<int>& types_config) {
  //Shen, MK. BIT (1962) 2: 228. doi:10.1007/BF01940170
  int _i=-1;
  int _j=-1;
  for(int i=1;i<(int)types_config.size();i++){if(types_config[i-1]<types_config[i]&&(i>_i)){_i=i;}}
  if(_i==-1){return false;} //stop condition
  for(int j=0;j<(int)types_config.size();j++){if(types_config[_i-1]<types_config[j]&&(j>_j)){_j=j;}}
  //cerr << "i=" << _i << "  j=" << _j << " ";
  std::swap(types_config[_i-1],types_config[_j]);
  for(int i=0;i<((int)types_config.size()-_i+1)/2;i++){std::swap(types_config[_i+i],types_config[types_config.size()-i-1]);}
  return true;
}

//simultaneously define and calculate
double POccCalculator::getEnergy() {ps_str.setBonds(v_types_config);return ps_str.getEnergy();}

void POccCalculator::getTotalPermutationsCount(){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy = "POccCalculator::getTotalPermutationsCount():";
  stringstream message;
  
  message << "Getting total number of derivative structure permutations";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  
  //unsigned long long int hnf_count=0;
  //unsigned long long int types_config_permutations_count=0;   //per hnf matrix

  //starting criteria for HNF matrices
  //int a_start=1;int c_start=1;int f_start=1;
  //int b_start=0;int d_start=0;int e_start=0;
  xmatrix<double> hnf_mat;                            //really xmatrix of int's, but we rule for int * double in xmatrix, no big deal
  //vector<xmatrix<double> > v_unique_superlattices;    //only store locally, as we need to make comparisons
  
  //START KESONG
  //cerr << "START KESONG " << endl;
  //CalculateHNF(p_str.xstr_pocc,p_str.n_hnf);

  //cerr << endl;
  //cerr << "START COREY " << endl;
  //resetMaxSLRadius();
  hnf_count=0;
  resetHNFMatrices();
  while(getNextHNFMatrix()){hnf_count++;}
	//cerr << "max_radius " << max_superlattice_radius << endl;
  //for(uint i=0;i<v_unique_superlattices.size();i++){
  //  //cerr << LatticeDimensionSphere(v_unique_superlattices[i],max_superlattice_radius/24.0) << endl;
  //  xvector<int> dims = LatticeDimensionSphere(v_unique_superlattices[i],10.0);
  //  cerr << "dims = " << dims << endl;
  //  xmatrix<double> scell(3,3);
  //  scell(1,1)=dims(1);scell(2,2)=dims(2);scell(3,3)=dims(3);
  //  xmatrix<double> newlattice = scell * v_unique_superlattices[i];
  //  cerr << "radius = " << RadiusSphereLattice(newlattice) << endl;
  //}
  //if(ENUMERATE_ALL_HNF){exit(0);}  //for debugging purposes - EXCEPTION
  message << "Total count of unique HNF matrices = " << hnf_count;
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  //exit(0);

  //starting criteria for site combinations
  //vector<vector<int> > v_types_config;
  //vector<int> v_config_iterators; //, v_site_iterators;
  
  //while(getNextSiteConfiguration(v_types_config,v_config_iterators,v_site_iterators)){
  
  // COREY MAKE THIS A TRUE CALCULATOR, simply go through each site, add between configs, 
  //multiple across sites
  types_config_permutations_count=0;   //per hnf matrix
  unsigned long long int str_config_permutations_count;
  unsigned long long int config_permutations_count=0;         //per hnf matrix
  total_permutations_count=0;  //per hnf matrix
  vector<int> current_config;
  const vector<StructureConfiguration>& v_str_configs=p_str.v_str_configs;
  bool getNextSiteConfiguration(vector<int>& types_config);
  for(uint str_config=0;str_config<v_str_configs.size();str_config++){
    str_config_permutations_count=1;
    for(uint site=0;site<v_str_configs[str_config].site_configs.size();site++){
      //for(uint config=0;config<v_str_configs[site].size();config++){
      config_permutations_count=0;
      current_config=v_str_configs[str_config].site_configs[site].getStartingTypesConfiguration();
      //if(config>0){config_permutations_count++;}
      config_permutations_count++;	//for starting config
      while(getNextSiteConfiguration(current_config)){config_permutations_count++;}
      str_config_permutations_count*=config_permutations_count;
      //}
      //str_config_permutations_count*=config_permutations_count;
    }
    types_config_permutations_count+=str_config_permutations_count;
  }
  //cerr << types_config_permutations_count << endl;
  //exit(0);
  message << "Total count of unique types-configuration permutations = " << types_config_permutations_count;
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  
  if(LDEBUG){
    cerr << soliloquy << " Checking unique types-configuration permutation count" << endl;
    types_config_permutations_count=0;
    resetSiteConfigurations();
    while(POccCalculator::getNextSiteConfiguration()){
      types_config_permutations_count++; 
      //cerr << "WOW1 " << v_types_config.size() << endl;
      //for(uint i=0;i<v_types_config.size();i++){
      //  cerr << "LOOK i=" << i << " " << v_types_config[i].size() << endl;
      //}
      //exit(0);
      cerr << soliloquy << " Permutation " << types_config_permutations_count << "  "; 
      for(uint i=0;i<v_types_config.size();i++){
        for(uint j=0;j<v_types_config[i].size();j++){
          cerr << v_types_config[i][j] << " ";
        } 
        cerr << "    ";
      }
      cerr << endl;
    }
    cerr << soliloquy << "Total count of unique types-configuration permutations = " << types_config_permutations_count << endl;
  }

  total_permutations_count=hnf_count*types_config_permutations_count;
  
  message << "Total number of derivative structure permutations = " <<total_permutations_count;
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

}

void POccSuperStructure::copyAtomAttributes(const _atom& atom_qualities,_atom& atom_position){
  //atom_position will be updated with qualities from atom_qualities
  //atom_qualities will not be updated
  xvector<double> fpos,cpos;
  fpos=atom_position.fpos; cpos=atom_position.cpos;
  atom_position=atom_qualities;
  atom_position.fpos=fpos;atom_position.cpos=cpos;
  /*
  _atom atom_clean;
  //grab data from atom_position
  atom_clean.fpos=atom_position.fpos;
  atom_clean.cpos=atom_position.cpos;
  //grab data from atom_qualities
  atom_clean.type=atom_qualities.type;
  atom_clean.spin=atom_qualities.spin;
  atom_clean.mass=atom_qualities.mass;
  atom_clean.name=atom_qualities.name;
  atom_clean.cleanname=atom_qualities.cleanname;
  atom_clean.atomic_number=atom_qualities.atomic_number;
  atom_clean.name_is_given=atom_qualities.name_is_given;
  //end properties transfer
  //fully replace
  atom_position=atom_clean;
  */
}

//void POccCalculator::getUFFParamBond(uint type1,uint type2,UFFParamBond& uffb){
//}

/*double POccCalculator::getBondEnergy(xstructure& xstr,vector<uint>& v_vacancies,
    xmatrix<double>& distance_matrix,vector<vector<uint> >& v_bonded_atom_indices,
    uint MODE){
 
  double energy=0.0;
  bool found_vacancy;
  uint atom1,atom2,type1,type2;
  UFFParamBond uffb;
  
  for(uint i=0;i<v_bonded_atom_indices.size();i++){
    found_vacancy=false;
    for(uint vi=0;vi<v_vacancies.size() && !found_vacancy;vi++){
      found_vacancy=(v_bonded_atom_indices[i][0]==v_vacancies[vi] ||
          v_bonded_atom_indices[i][1]==v_vacancies[vi]);
    }
    if(found_vacancy){cerr<<"VACANCYYYY"<<endl;continue;}
    atom1=v_bonded_atom_indices[i][0]; atom2=v_bonded_atom_indices[i][1];
    type1=xstr.atoms[atom1].type; type2=xstr.atoms[atom2].type;
    //uffb.uffp1=&types2uffparams_map[type1]; uffb.uffp2=&types2uffparams_map[type2];
    //uffb.distij=distance_matrix(atom1,atom2);
    //simultaneously initialize and calculate uffb
    uffb.calculate(types2uffparams_map[type1],types2uffparams_map[type2],distance_matrix(atom1,atom2));
    if(MODE==BOND_MODE){energy += 0.5 * uffb.Kij * uffb.delta * uffb.delta;}
    else if(MODE==NONBOND_MODE){energy += uffb.Dij * (uffb.X12 - 2.0 * uffb.X6);}
    else{exit(1);}
  }
  energy*=KCAL_2_EV;
  return energy;
}*/

/*xstructure POccSuperStructure::getXStructure(vector<vector<int> >& _v_types_config){
  v_types_config=&_v_types_config;
  
  const xstructure& xstr_nopocc = p_str.xstr_nopocc;
  const vector<Site>& v_sites = p_str.v_sites;
  vector<int> sc2pc_map, pc2sc_map;
  xstructure supercell=GetSuperCell(xstr_nopocc,hnf_mat,sc2pc_map,pc2sc_map,false,false);
  
  uint starting_supercell_atom_index,pocc_atom_index,supercell_atom_index;
  vector<uint> v_atoms_to_remove;
  for(uint site=0;site<v_types_config.size();site++){
    if(v_sites[site].partial_occupation_flag){  //slight optimization
      //test of stupidity
      if(v_sites[site].v_occupants.size()<2){cerr << "WHOOOOOOPS" << endl;exit(1);}
      //find index of atom first in v_occupants list, that's the one that remained
      starting_supercell_atom_index=pc2sc_map[v_sites[site].v_occupants[0]];
      for(uint i=0;i<v_types_config[site].size();i++){
        //switch these atoms
        //supercell.atoms[starting_supercell_atom_index+i]
        //xstr_pocc.atoms[types2pc_map[v_types_config[site][i]]]
        //pocc_atom_index=types2pc_map[v_types_config[site][i]];
        pocc_atom_index=p_str.types2pc_map[v_types_config[site][i]];
        supercell_atom_index=starting_supercell_atom_index+i;
        const _atom& pocc_atom = p_str.xstr_pocc.atoms[pocc_atom_index];
        if(v_types_config[site][i]>=0){copyAtomAttributes(pocc_atom,supercell.atoms[supercell_atom_index]);}
        else{v_atoms_to_remove.push_back(supercell_atom_index);}
      }
    }
  }
  //std::sort(v_atoms_to_remove.rbegin(),v_atoms_to_remove.rend()); //NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
  //for(uint atom=0;atom<v_atoms_to_remove.size();atom++){supercell.RemoveAtom(v_atoms_to_remove[atom]);}

  //calculate UFF energy here before you rearrange the structure
  double energy=0.0;
  //energy += getBondEnergy(supercell,v_atoms_to_remove,distance_matrix,v_bonded_atom_indices,BOND_MODE);
  //energy += getBondEnergy(supercell,v_atoms_to_remove,distance_matrix,v_nonbonded_atom_indices,NONBOND_MODE);

  cerr << energy << endl;



  //rearrange the structure below here

  //RemoveAtoms(supercell,v_atoms_to_remove);
  supercell.RemoveAtoms(v_atoms_to_remove);
  //supercell.SpeciesPutAlphabetic();

  //cerr << pocc::CalculateUFFEnergy(supercell) << endl;
  //cerr << supercell << endl;
  //exit(0);
  return supercell;
}*/

bool POccStructure::getUFFParameters(const string& atom,string& params){
  params.clear();
  if (atom=="H")       {params="H   0.354   180.000 2.886   0.044   12.000  0.712   0.000   0.000   4.528   6.945   0.371   ";return true;}   
  else if (atom=="He") {params="He  0.849   90.000  2.362   0.056   15.240  0.098   0.000   0.000   9.660   14.920  1.300   ";return true;}   
  else if (atom=="Li") {params="Li  1.336   180.000 2.451   0.025   12.000  1.026   0.000   2.000   3.006   2.386   1.557   ";return true;}   
  else if (atom=="Be") {params="Be  1.074   109.470 2.745   0.085   12.000  1.565   0.000   2.000   4.877   4.443   1.240   ";return true;}  
  else if (atom=="B")  {params="B   0.838   109.470 4.083   0.180   12.052  1.755   0.000   2.000   5.110   4.750   0.822   ";return true;}  
  else if (atom=="C")  {params="C   0.706   180.000 3.851   0.105   12.730  1.912   0.000   2.000   5.343   5.063   0.759   ";return true;}  
  else if (atom=="N")  {params="N   0.656   180.000 3.660   0.069   13.407  2.544   0.000   2.000   6.899   5.880   0.715   ";return true;}  
  else if (atom=="O")  {params="O   0.639   180.000 3.500   0.060   14.085  2.300   0.000   2.000   8.741   6.682   0.669   ";return true;}  
  else if (atom=="F")  {params="F   0.668   180.000 3.364   0.050   14.762  1.735   0.000   2.000   10.874  7.474   0.706   ";return true;}  
  else if (atom=="Ne") {params="Ne  0.920   90.000  3.243   0.042   15.440  0.194   0.000   2.000   11.040  10.550  1.768   ";return true;}  
  else if (atom=="Na") {params="Na  1.539   180.000 2.983   0.030   12.000  1.081   0.000   1.250   2.843   2.296   2.085   ";return true;}  
  else if (atom=="Mg") {params="Mg  1.421   109.470 3.021   0.111   12.000  1.787   0.000   1.250   3.951   3.693   1.500   ";return true;}  
  else if (atom=="Al") {params="Al  1.244   109.470 4.499   0.505   11.278  1.792   0.000   1.250   4.060   3.590   1.201   ";return true;}  
  else if (atom=="Si") {params="Si  1.117   109.470 4.295   0.402   12.175  2.323   1.225   1.250   4.168   3.487   1.176   ";return true;}  
  else if (atom=="P")  {params="P   1.101   93.800  4.147   0.305   13.072  2.863   2.400   1.250   5.463   4.000   1.102   ";return true;}  
  else if (atom=="S")  {params="S   1.064   92.1000 4.035   0.274   13.969  2.703   0.000   1.250   6.928   4.486   1.047   ";return true;}  
  else if (atom=="Cl") {params="Cl  1.044   180.000 3.947   0.227   14.866  2.348   0.000   1.250   8.564   4.946   0.994   ";return true;}  
  else if (atom=="Ar") {params="Ar  1.032   90.000  3.868   0.185   15.763  0.300   0.000   1.250   9.465   6.355   2.108   ";return true;}  
  else if (atom=="K")  {params="K   1.953   180.000 3.812   0.035   12.000  1.165   0.000   0.700   2.421   1.920   2.586   ";return true;}  
  else if (atom=="Ca") {params="Ca  1.761   90.000  3.399   0.238   12.000  2.141   0.000   0.700   3.231   2.880   2.000   ";return true;}  
  else if (atom=="Sc") {params="Sc  1.513   109.470 3.295   0.019   12.000  2.592   0.000   0.700   3.395   3.080   1.750   ";return true;}  
  else if (atom=="Ti") {params="Ti  1.412   109.470 3.175   0.017   12.000  2.659   0.000   0.700   3.470   3.380   1.607   ";return true;}  
  else if (atom=="V")  {params="V   1.402   109.470 3.144   0.016   12.000  2.679   0.000   0.700   3.650   3.410   1.470   ";return true;}  
  else if (atom=="Cr") {params="Cr  1.345   90.000  3.023   0.015   12.000  2.463   0.000   0.700   3.415   3.865   1.402   ";return true;}  
  else if (atom=="Mn") {params="Mn  1.382   90.000  2.961   0.013   12.000  2.430   0.000   0.700   3.325   4.105   1.533   ";return true;}  
  else if (atom=="Fe") {params="Fe  1.270   109.470 2.912   0.013   12.000  2.430   0.000   0.700   3.760   4.140   1.393   ";return true;}  
  else if (atom=="Co") {params="Co  1.241   90.000  2.872   0.014   12.000  2.430   0.000   0.700   4.105   4.175   1.406   ";return true;}  
  else if (atom=="Ni") {params="Ni  1.164   90.000  2.834   0.015   12.000  2.430   0.000   0.700   4.465   4.205   1.398   ";return true;}  
  else if (atom=="Cu") {params="Cu  1.302   109.470 3.495   0.005   12.000  1.756   0.000   0.700   4.200   4.220   1.434   ";return true;}  
  else if (atom=="Zn") {params="Zn  1.193   109.470 2.763   0.124   12.000  1.308   0.000   0.700   5.106   4.285   1.400   ";return true;}  
  else if (atom=="Ga") {params="Ga  1.260   109.470 4.383   0.415   11.000  1.821   0.000   0.700   3.641   3.160   1.211   ";return true;}  
  else if (atom=="Ge") {params="Ge  1.197   109.470 4.280   0.379   12.000  2.789   0.701   0.700   4.051   3.438   1.189   ";return true;}  
  else if (atom=="As") {params="As  1.211   92.100  4.230   0.309   13.000  2.864   1.500   0.700   5.188   3.809   1.204   ";return true;}  
  else if (atom=="Se") {params="Se  1.190   90.600  4.205   0.291   14.000  2.764   0.335   0.700   6.428   4.131   1.224   ";return true;}  
  else if (atom=="Br") {params="Br  1.192   180.000 4.189   0.251   15.000  2.519   0.000   0.700   7.790   4.425   1.141   ";return true;}  
  else if (atom=="Kr") {params="Kr  1.147   90.000  4.141   0.220   16.000  0.452   0.000   0.700   8.505   5.715   2.270   ";return true;}  
  else if (atom=="Rb") {params="Rb  2.260   180.000 4.114   0.040   12.000  1.592   0.000   0.200   2.331   1.846   2.770   ";return true;}  
  else if (atom=="Sr") {params="Sr  2.052   90.000  3.641   0.235   12.000  2.449   0.000   0.200   3.024   2.440   2.415   ";return true;}  
  else if (atom=="Y")  {params="Y   1.698   109.470 3.345   0.072   12.000  3.257   0.000   0.200   3.830   2.810   1.998   ";return true;}  
  else if (atom=="Zr") {params="Zr  1.564   109.470 3.124   0.069   12.000  3.667   0.000   0.200   3.400   3.550   1.758   ";return true;}  
  else if (atom=="Nb") {params="Nb  1.473   109.470 3.165   0.059   12.000  3.618   0.000   0.200   3.550   3.380   1.603   ";return true;}  
  else if (atom=="Mo") {params="Mo  1.467   90.000  3.052   0.056   12.000  3.400   0.000   0.200   3.465   3.755   1.530   ";return true;}  
  else if (atom=="Tc") {params="Tc  1.322   90.000  2.998   0.048   12.000  3.400   0.000   0.200   3.290   3.990   1.500   ";return true;}  
  else if (atom=="Ru") {params="Ru  1.478   90.000  2.963   0.056   12.000  3.400   0.000   0.200   3.575   4.015   1.500   ";return true;}  
  else if (atom=="Rh") {params="Rh  1.332   90.000  2.929   0.053   12.000  3.500   0.000   0.200   3.975   4.005   1.509   ";return true;}  
  else if (atom=="Pd") {params="Pd  1.338   90.000  2.899   0.048   12.000  3.210   0.000   0.200   4.320   4.000   1.544   ";return true;}  
  else if (atom=="Ag") {params="Ag  1.386   180.000 3.148   0.036   12.000  1.956   0.000   0.200   4.436   3.134   1.622   ";return true;}  
  else if (atom=="Cd") {params="Cd  1.403   109.470 2.848   0.228   12.000  1.650   0.000   0.200   5.034   3.957   1.600   ";return true;}  
  else if (atom=="In") {params="In  1.459   109.470 4.463   0.599   11.000  2.070   0.000   0.200   3.506   2.896   1.404   ";return true;}  
  else if (atom=="Sn") {params="Sn  1.398   109.470 4.392   0.567   12.000  2.961   0.199   0.200   3.987   3.124   1.354   ";return true;}  
  else if (atom=="Sb") {params="Sb  1.407   91.600  4.420   0.449   13.000  2.704   1.100   0.200   4.899   3.342   1.404   ";return true;}  
  else if (atom=="Te") {params="Te  1.386   90.250  4.470   0.398   14.000  2.882   0.300   0.200   5.816   3.526   1.380   ";return true;}  
  else if (atom=="I")  {params="I   1.382   180.000 4.500   0.339   15.000  2.650   0.000   0.200   6.822   3.762   1.333   ";return true;}  
  else if (atom=="Xe") {params="Xe  1.267   90.000  4.404   0.332   12.000  0.556   0.000   0.200   7.595   4.975   2.459   ";return true;}  
  else if (atom=="Cs") {params="Cs  2.570   180.000 4.517   0.045   12.000  1.573   0.000   0.100   2.183   1.711   2.984   ";return true;}  
  else if (atom=="Ba") {params="Ba  2.277   90.000  3.703   0.364   12.000  2.727   0.000   0.100   2.814   2.396   2.442   ";return true;}  
  else if (atom=="La") {params="La  1.943   109.470 3.522   0.017   12.000  3.300   0.000   0.100   2.836   2.742   2.071   ";return true;}  
  else if (atom=="Ce") {params="Ce  1.841   90.000  3.556   0.013   12.000  3.300   0.000   0.100   2.774   2.692   1.925   ";return true;}  
  else if (atom=="Pr") {params="Pr  1.823   90.000  3.606   0.010   12.000  3.300   0.000   0.100   2.858   2.564   2.007   ";return true;}  
  else if (atom=="Nd") {params="Nd  1.816   90.000  3.575   0.010   12.000  3.300   0.000   0.100   2.869   2.621   2.007   ";return true;}  
  else if (atom=="Pm") {params="Pm  1.801   90.000  3.547   0.009   12.000  3.300   0.000   0.100   2.881   2.673   2.000   ";return true;}  
  else if (atom=="Sm") {params="Sm  1.780   90.000  3.520   0.008   12.000  3.300   0.000   0.100   2.912   2.720   1.978   ";return true;}  
  else if (atom=="Eu") {params="Eu  1.771   90.000  3.493   0.008   12.000  3.300   0.000   0.100   2.879   2.788   2.227   ";return true;}  
  else if (atom=="Gd") {params="Gd  1.735   90.000  3.368   0.009   12.000  3.300   0.000   0.100   3.167   2.975   1.968   ";return true;}  
  else if (atom=="Tb") {params="Tb  1.732   90.000  3.451   0.007   12.000  3.300   0.000   0.100   3.018   2.834   1.954   ";return true;}  
  else if (atom=="Dy") {params="Dy  1.710   90.000  3.428   0.007   12.000  3.300   0.000   0.100   3.056   2.872   1.934   ";return true;}  
  else if (atom=="Ho") {params="Ho  1.696   90.000  3.409   0.007   12.000  3.416   0.000   0.100   3.127   2.891   1.925   ";return true;}  
  else if (atom=="Er") {params="Er  1.673   90.000  3.391   0.007   12.000  3.300   0.000   0.100   3.187   2.915   1.915   ";return true;}  
  else if (atom=="Tm") {params="Tm  1.660   90.000  3.374   0.006   12.000  3.300   0.000   0.100   3.251   2.933   2.000   ";return true;}  
  else if (atom=="Yb") {params="Yb  1.637   90.000  3.355   0.228   12.000  2.618   0.000   0.100   3.289   2.965   2.158   ";return true;}  
  else if (atom=="Lu") {params="Lu  1.671   90.000  3.640   0.041   12.000  3.271   0.000   0.100   2.963   2.463   1.896   ";return true;}  
  else if (atom=="Hf") {params="Hf  1.611   109.470 3.141   0.072   12.000  3.921   0.000   0.100   3.700   3.400   1.759   ";return true;}  
  else if (atom=="Ta") {params="Ta  1.511   109.470 3.170   0.081   12.000  4.075   0.000   0.100   5.100   2.850   1.605   ";return true;}  
  else if (atom=="W")  {params="W   1.392   90.000  3.069   0.067   12.000  3.700   0.000   0.100   4.630   3.310   1.538   ";return true;}  
  else if (atom=="Re") {params="Re  1.372   90.000  2.954   0.066   12.000  3.700   0.000   0.100   3.960   3.920   1.600   ";return true;}  
  else if (atom=="Os") {params="Os  1.372   90.000  3.120   0.037   12.000  3.700   0.000   0.100   5.140   3.630   1.700   ";return true;}  
  else if (atom=="Ir") {params="Ir  1.371   90.000  2.840   0.073   12.000  3.731   0.000   0.100   5.000   4.000   1.866   ";return true;}  
  else if (atom=="Pt") {params="Pt  1.364   90.000  2.754   0.080   12.000  3.382   0.000   0.100   4.790   4.430   1.557   ";return true;}  
  else if (atom=="Au") {params="Au  1.262   90.000  3.293   0.039   12.000  2.625   0.000   0.100   4.894   2.586   1.618   ";return true;}  
  else if (atom=="Hg") {params="Hg  1.340   180.000 2.705   0.385   12.000  1.750   0.000   0.100   6.270   4.160   1.600   ";return true;}  
  else if (atom=="Tl") {params="Tl  1.518   120.000 4.347   0.680   11.000  2.068   0.000   0.100   3.200   2.900   1.530   ";return true;}  
  else if (atom=="Pb") {params="Pb  1.459   109.470 4.297   0.663   12.000  2.846   0.100   0.100   3.900   3.530   1.444   ";return true;}  
  else if (atom=="Bi") {params="Bi  1.512   90.000  4.370   0.518   13.000  2.470   1.000   0.100   4.690   3.740   1.514   ";return true;}  
  else if (atom=="Po") {params="Po  1.500   90.000  4.709   0.325   14.000  2.330   0.300   0.100   4.210   4.210   1.480   ";return true;}  
  else if (atom=="At") {params="At  1.545   180.000 4.750   0.284   15.000  2.240   0.000   0.100   4.750   4.750   1.470   ";return true;}  
  else if (atom=="Rn") {params="Rn  1.420   90.000  4.765   0.248   16.000  0.583   0.000   0.100   5.370   5.370   2.200   ";return true;}  
  else if (atom=="Fr") {params="Fr  2.880   180.000 4.900   0.050   12.000  1.847   0.000   0.000   2.000   2.000   2.300   ";return true;}  
  else if (atom=="Ra") {params="Ra  2.512   90.000  3.677   0.404   12.000  2.920   0.000   0.000   2.843   2.434   2.200   ";return true;}  
  else if (atom=="Ac") {params="Ac  1.983   90.000  3.478   0.033   12.000  3.900   0.000   0.000   2.835   2.835   2.108   ";return true;}  
  else if (atom=="Th") {params="Th  1.721   90.000  3.396   0.026   12.000  4.202   0.000   0.000   3.175   2.905   2.018   ";return true;}  
  else if (atom=="Pa") {params="Pa  1.711   90.000  3.424   0.022   12.000  3.900   0.000   0.000   2.985   2.905   1.800   ";return true;}  
  else if (atom=="U")  {params="U   1.684   90.000  3.395   0.022   12.000  3.900   0.000   0.000   3.341   2.853   1.713   ";return true;}  
  else if (atom=="Np") {params="Np  1.666   90.000  3.424   0.019   12.000  3.900   0.000   0.000   3.549   2.717   1.800   ";return true;}  
  else if (atom=="Pu") {params="Pu  1.657   90.000  3.424   0.016   12.000  3.900   0.000   0.000   3.243   2.819   1.840   ";return true;}  
  else if (atom=="Am") {params="Am  1.660   90.000  3.381   0.014   12.000  3.900   0.000   0.000   2.990   3.004   1.942   ";return true;}  
  else if (atom=="Cm") {params="Cm  1.801   90.000  3.326   0.013   12.000  3.900   0.000   0.000   2.832   3.190   1.900   ";return true;}  
  else if (atom=="Bk") {params="Bk  1.761   90.000  3.339   0.013   12.000  3.900   0.000   0.000   3.194   3.036   1.900   ";return true;}  
  else if (atom=="Cf") {params="Cf  1.750   90.000  3.313   0.013   12.000  3.900   0.000   0.000   3.197   3.101   1.900   ";return true;}  
  else if (atom=="Es") {params="Es  1.724   90.000  3.299   0.012   12.000  3.900   0.000   0.000   3.333   3.089   1.900   ";return true;}  
  else if (atom=="Fm") {params="Fm  1.712   90.000  3.286   0.012   12.000  3.900   0.000   0.000   3.400   3.100   1.900   ";return true;}  
  else if (atom=="Md") {params="Md  1.689   90.000  3.274   0.011   12.000  3.900   0.000   0.000   3.470   3.110   1.900   ";return true;}  
  else if (atom=="No") {params="No  1.679   90.000  3.248   0.011   12.000  3.900   0.000   0.000   3.475   3.175   1.900   ";return true;}  
  else if (atom=="Lw") {params="Lw  1.698   90.000  3.236   0.011   12.000  3.900   0.000   0.000   3.500   3.200   1.900   ";return true;}  
  else {return false;}                                                                                                      
}

//void POccCalculator::getBonding(xmatrix<double>& hnf_mat,xmatrix<double>& _distance_matrix,vector<vector<uint> >& v_bonded_atom_indices,
//    vector<vector<uint> >& v_nonbonded_atom_indices) {
  //vector<int> sc2pc_map, pc2sc_map;
  //xstructure supercell=GetSuperCell(xstr_nopocc,hnf_mat,sc2pc_map,pc2sc_map,false,false);

  //[OBSOLETE]xvector<double> min_vec; xvector<int> ijk;  //dummy
  //[OBSOLETE]xmatrix<double> distance_matrix(supercell.atoms.size()-1,supercell.atoms.size()-1,0,0);
  //[OBSOLETE]vector<double> v_nn_dists;
  //[OBSOLETE]for(uint atom1=0;atom1<supercell.atoms.size();atom1++){v_nn_dists.push_back(std::numeric_limits<double>::max());} //initialize with big double
  //[OBSOLETE]for(uint atom1=0;atom1<supercell.atoms.size();atom1++){
  //[OBSOLETE]  for(uint atom2=atom1+1;atom2<supercell.atoms.size();atom2++){
  //[OBSOLETE]    distance_matrix(atom1,atom2)=distance_matrix(atom2,atom1)=SYM::minimumCartesianDistance(supercell.atoms[atom1].cpos,supercell.atoms[atom2].cpos,supercell.lattice,min_vec,ijk);
  //[OBSOLETE]    if(distance_matrix(atom1,atom2)<v_nn_dists[atom1]){v_nn_dists[atom1]=distance_matrix(atom1,atom2);}
  //[OBSOLETE]  }
  //[OBSOLETE]}
  //[OBSOLETE]_distance_matrix=distance_matrix;

  //[OBSOLETE]for(uint atom1=0;atom1<supercell.atoms.size();atom1++){
  //[OBSOLETE]  for(uint atom2=atom1+1;atom2<supercell.atoms.size();atom2++){
  //[OBSOLETE]    if(abs(distance_matrix(atom1,atom2)-v_nn_dists[atom1])<0.5){  //kesong standard for bonding, keep for now
  //[OBSOLETE]      v_bonded_atom_indices.push_back(vector<uint>(0));
  //[OBSOLETE]      v_bonded_atom_indices.back().push_back(atom1);
  //[OBSOLETE]      v_bonded_atom_indices.back().push_back(atom2);
  //[OBSOLETE]      //cerr << "BOND   " << atom1 << " " << atom2 << endl;
  //[OBSOLETE]    }else{
  //[OBSOLETE]      v_nonbonded_atom_indices.push_back(vector<uint>(0));
  //[OBSOLETE]      v_nonbonded_atom_indices.back().push_back(atom1);
  //[OBSOLETE]      v_nonbonded_atom_indices.back().push_back(atom2);
  //[OBSOLETE]      //cerr << "NOBOND " << atom1 << " " << atom2 << endl;
  //[OBSOLETE]    }
  //[OBSOLETE]  }
  //[OBSOLETE]}

  //[OBSOLETE]//cerr << distance_matrix << endl;
  //[OBSOLETE]//exit(0);
  

//}

bool POccCalculator::areEquivalentStructuresByUFF(std::list<POccReducedDerivativeStructuresSet>::iterator it, const POccReducedDerivativeStructure& prds) const {
  bool energy_equal=aurostd::isequal((*it).getEnergy(),prds.energy,ENERGY_TOLERANCE);
  bool hnf_equal=((*it).getHNFIndex()==prds.hnf_index);
  return energy_equal && hnf_equal;
}

//NEW
void POccCalculator::add2DerivativeStructuresList(const POccReducedDerivativeStructure& prds,
    std::list<POccReducedDerivativeStructuresSet>::iterator i_start,
    std::list<POccReducedDerivativeStructuresSet>::iterator i_end){
//  cerr << "FULL LIST START" << endl;
//  for(std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();it!=l_derivative_structures_sets.end();++it){
//    cerr << std::fixed << std::setprecision(9) << (*it).getEnergy() << " ";
//  }
//  cerr << endl;
//  cerr << "FULL LIST END" << endl;
//  cerr << "CURRENT " << prds.energy << endl;
//  cerr << "LOOK start " << std::distance(l_derivative_structures_sets.begin(),i_start) << " " << (*i_start).getEnergy() << endl;
//  cerr << "LOOK end " << std::distance(l_derivative_structures_sets.begin(),i_end) << " " << (*i_end).getEnergy() << endl;
  if(i_start==i_end){ //std::distance(i_start,i_end)==0){
//    cerr << "start==stop" << endl;
//
    if(i_start==l_derivative_structures_sets.end()){--i_start;}

    if(areEquivalentStructuresByUFF(i_start,prds)){//aurostd::isequal(prds.energy,(*i_start).getEnergy(),ENERGY_TOLERANCE)){
//      cerr << "found degeneracy with " << std::distance(l_derivative_structures_sets.begin(),i_start) << endl;
      (*i_start).prds_set.push_back(prds);
      //(*i_start).degeneracy++;
      return;
    }

    if(prds.energy>(*i_start).getEnergy()){++i_start;}
    //  if(i_start!=l_derivative_structures_sets.end()){++i_start;}
    //}
//    cerr << "Adding to " << std::distance(l_derivative_structures_sets.begin(),i_start) << endl;
    POccReducedDerivativeStructuresSet prdss; prdss.prds_set.push_back(prds);
    l_derivative_structures_sets.insert(i_start,prdss);
    //l_derivative_structures_sets.insert(i_start,prds);
    return;
  }
//  cerr << "start!=stop" << endl;
  std::list<POccReducedDerivativeStructuresSet>::iterator i_middle=i_start;
  //std::advance(i_middle,std::distance(l_derivative_structures_sets.begin(),i_start));
  std::advance(i_middle,std::distance(i_start,i_end)/2);
//  cerr << "Looking at middle=" << std::distance(l_derivative_structures_sets.begin(),i_middle) << " " << (*i_middle).getEnergy() << endl;
  //cerr << "SEE HERE " << i_middle << endl;
  if(areEquivalentStructuresByUFF(i_middle,prds)){ //aurostd::isequal(prds.energy,(*i_middle).getEnergy(),ENERGY_TOLERANCE)){
//    cerr << "found degeneracy with " << std::distance(l_derivative_structures_sets.begin(),i_middle) << endl;
    (*i_middle).prds_set.push_back(prds);
    //(*i_middle).degeneracy++;
    return;
  }
  if(prds.energy<(*i_middle).getEnergy()){
//    cerr << "less than middle" << endl;
    return add2DerivativeStructuresList(prds,i_start,i_middle);
  }
//  cerr << "greater than middle" << endl;
  return add2DerivativeStructuresList(prds,++i_middle,i_end);
}

//ORIGINAL
//void POccCalculator::add2DerivativeStructuresList(POccDerivativeStructure& pds,
//    std::list<POccDerivativeStructure>::iterator i_start,
//    std::list<POccDerivativeStructure>::iterator i_end){
////  cerr << "FULL LIST START" << endl;
////  for(std::list<POccDerivativeStructure>::iterator it=l_derivative_structures.begin();it!=l_derivative_structures.end();++it){
////    cerr << std::fixed << std::setprecision(9) << (*it).energy << " ";
////  }
////  cerr << endl;
////  cerr << "FULL LIST END" << endl;
////  cerr << "CURRENT " << pds.energy << endl;
////  cerr << "LOOK start " << std::distance(l_derivative_structures.begin(),i_start) << " " << (*i_start).energy << endl;
////  cerr << "LOOK end " << std::distance(l_derivative_structures.begin(),i_end) << " " << (*i_end).energy << endl;
//  if(i_start==i_end){ //std::distance(i_start,i_end)==0){
////    cerr << "start==stop" << endl;
//    if(aurostd::isequal(pds.energy,(*i_start).energy,ENERGY_TOLERANCE)){
////      cerr << "found degeneracy with " << std::distance(l_derivative_structures.begin(),i_start) << endl;
//      (*i_start).degeneracy++;
//      return;
//    }
//    if(pds.energy>(*i_start).energy){
//      if(i_start!=l_derivative_structures.end()){++i_start;}
//    }
////    cerr << "Adding to " << std::distance(l_derivative_structures.begin(),i_start) << endl;
//    l_derivative_structures.insert(i_start,pds);
//    return;
//  }
////  cerr << "start!=stop" << endl;
//  std::list<POccDerivativeStructure>::iterator i_middle=i_start;
//  //std::advance(i_middle,std::distance(l_derivative_structures.begin(),i_start));
//  std::advance(i_middle,std::distance(i_start,i_end)/2);
////  cerr << "Looking at middle=" << std::distance(l_derivative_structures.begin(),i_middle) << " " << (*i_middle).energy << endl;
//  //cerr << "SEE HERE " << i_middle << endl;
//  if(aurostd::isequal(pds.energy,(*i_middle).energy,ENERGY_TOLERANCE)){
////    cerr << "found degeneracy with " << std::distance(l_derivative_structures.begin(),i_middle) << endl;
//    (*i_middle).degeneracy++;
//    return;
//  }
//  if(pds.energy<(*i_middle).energy){
////    cerr << "less than middle" << endl;
//    return add2DerivativeStructuresList(pds,i_start,i_middle);
//  }
////  cerr << "greater than middle" << endl;
//  return add2DerivativeStructuresList(pds,++i_middle,i_end);
//}

void POccCalculator::add2DerivativeStructuresList(const POccReducedDerivativeStructure& prds){
//  cerr << "SIZE " << l_derivative_structures_sets.size() << endl;
  if(l_derivative_structures_sets.size()==0){
    POccReducedDerivativeStructuresSet prdss; prdss.prds_set.push_back(prds);
    l_derivative_structures_sets.push_back(prdss);
    //l_derivative_structures.push_back(pds);
    return;
  }
  return add2DerivativeStructuresList(prds,l_derivative_structures_sets.begin(),l_derivative_structures_sets.end());
  //return add2DerivativeStructuresList(pds,l_derivative_structures.begin(),l_derivative_structures.end());
}

void updateProgressBar(unsigned long long int current, unsigned long long int end, ostream& oss){
  double progress = (double)current/(double)end;
  int pos = BAR_WIDTH * progress;

  oss << "[";
  for (int i = 0; i < BAR_WIDTH; ++i) {
    if (i < pos) oss << "=";
    else if (i == pos) oss << ">";
  	else oss << " ";
  }
	oss << "] " << int(progress * 100.0) << " %\r";
	oss.flush();
	if(current==end){ oss << endl; }
}

POccDerivativeStructure POccCalculator::getParticularPOccDerivativeStructure(const POccReducedDerivativeStructure& prds){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccCalculator::getParticularPOccDerivativeStructure():";
  if(LDEBUG){
    cerr << soliloquy << " prds.hnf_index=" << prds.hnf_index << endl;
    cerr << soliloquy << " prds.site_config_index=" << prds.site_config_index << endl;
    cerr << soliloquy << " prds.energy=" << prds.energy << endl;
  }
  unsigned long long int hnf_index=0;
  unsigned long long int site_config_index=0;
  //add an initialization flag for ps_str eventually, and check here
  POccDerivativeStructure pds;
  pds.hnf_index=prds.hnf_index;
  pds.site_config_index=prds.site_config_index;
  pds.energy=prds.energy;
  resetHNFMatrices();
  resetSiteConfigurations();
  while(getNextHNFMatrix()){  //sets hnf_mat
    if(hnf_index!=pds.hnf_index){hnf_index++;continue;}
    pds.hnf_mat=hnf_mat;  //abstract this away, confusing where it comes from
    while(getNextSiteConfiguration()){  //sets v_types_config
      if(site_config_index!=pds.site_config_index){site_config_index++;continue;}
      pds.v_types_config=v_types_config;  //abstract this away, confusing where it comes from
      return pds;
    }
  }
  throw AFLOWLogicError(soliloquy,"Invalid hnf site configuration indices");
}

bool POccCalculator::areEquivalentStructures(const POccReducedDerivativeStructure& prds_a,const POccReducedDerivativeStructure& prds_b) {
  POccDerivativeStructure pds_a,pds_b;
  pds_a=getParticularPOccDerivativeStructure(prds_a);
  pds_b=getParticularPOccDerivativeStructure(prds_b);

  xstructure a=ps_str.createDerivativeStructure(pds_a,p_str.n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
  xstructure b=ps_str.createDerivativeStructure(pds_b,p_str.n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
  return areEquivalentStructures(a,b);
}

bool POccCalculator::areEquivalentStructures(const xstructure& a, const xstructure& b) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccCalculator::areEquivalentStructures():";
  if(LDEBUG){
    cerr << soliloquy << " comparing structure:" << endl;
    cerr << a;
    cerr << "========================================== vs. ==========================================" << endl;
    cerr << b;
  }
  //xstructure aa("/home/corey/work/work/pranab/POCC/new_tests/new_sets_test_new_pocc/test/sc_test/POSCAR1",IOAFLOW_AUTO);
  //xstructure bb("/home/corey/work/work/pranab/POCC/new_tests/new_sets_test_new_pocc/test/sc_test/POSCAR2",IOAFLOW_AUTO);
  //cerr << aa << endl;
  //cerr << bb << endl;
  bool are_equivalent=compare::aflowCompareStructure(a,b,true,false,true); //match species and use fast match, but not scale volume, two structures with different volumes (pressures) are different! // DX 1/23/18 - added fast_match = true
  //cerr << are_equivalent << endl;
  //exit(1);
  if(LDEBUG){cerr << soliloquy << " structures are " << (are_equivalent?"":"NOT ") << "equivalent" << endl;}
  return are_equivalent;
}

unsigned long long int POccCalculator::runRobustStructureComparison(std::list<POccReducedDerivativeStructuresSet>::iterator it_prdss){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccCalculator::runRobustStructureComparison():";
  stringstream message;

  const vector<POccReducedDerivativeStructure>& vprds=(*it_prdss).prds_set;
  vector<vector<uint> > unique_structure_bins;
  unique_structure_bins.push_back(vector<uint>(0));
  unique_structure_bins.back().push_back(0);  //initialize with first structure

  //when you parallelize this, make a new structure that holds bool, equivalent or not
  //and then run through remaining set of structures

  unsigned long long int i_prdss=std::distance(l_derivative_structures_sets.begin(),it_prdss);

  message << "Performing robust structure comparison for structure group[" << i_prdss+1 << "]";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); 

  bool are_equivalent;
  xstructure a,b;
  uint starting_index=1;
  updateProgressBar(0,vprds.size()-starting_index,*p_oss);
  for(uint i=starting_index;i<vprds.size();i++){
    POccDerivativeStructure pds_b=getParticularPOccDerivativeStructure(vprds[i]); //fix me
    b=ps_str.createDerivativeStructure(pds_b,p_str.n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
    for(uint j=0;j<unique_structure_bins.size();j++){
      if(unique_structure_bins[j].size()==0){
        throw AFLOWLogicError(soliloquy,"No structure groups found");
      }
      if(LDEBUG){cerr << soliloquy << " comparing structure[" << unique_structure_bins[j][0]+1 << "] with structure[" << i+1 << "] (max=" << vprds.size() << ")" << endl;}
      POccDerivativeStructure pds_a=getParticularPOccDerivativeStructure(vprds[unique_structure_bins[j][0]]); //fix me
      a=ps_str.createDerivativeStructure(pds_a,p_str.n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
      are_equivalent=areEquivalentStructures(a,b);//vprds[unique_structure_bins[j][0]],vprds[i]);
      if(are_equivalent){unique_structure_bins[j].push_back(i);break;}
      //if(!are_equivalent){cerr << "AHHHHHHHHHHHHHHHHHH" << endl;}
      if(j==unique_structure_bins.size()-1){  //we've reached the end and no match, need to make a new bin
        unique_structure_bins.push_back(vector<uint>(0));
        unique_structure_bins.back().push_back(i);
      }
    }
    updateProgressBar(i,vprds.size(),*p_oss);
  }

  if(unique_structure_bins.size()==1){return unique_structure_bins.size();}  //they are all equivalent as dictated by UFF
  //otherwise, we need to rearrange structures

  if(LDEBUG){}
  message << "Splitting structure group[" << std::distance(l_derivative_structures_sets.begin(),it_prdss)+1;
  message << "] into " << unique_structure_bins.size() << " groups as determined by the robust structure comparison";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_); 

  std::list<POccReducedDerivativeStructuresSet>::iterator it=it_prdss;
  vector<uint> vi_prds_to_remove;
  for(uint i=1;i<unique_structure_bins.size();i++){
    ++it;
    POccReducedDerivativeStructuresSet prdss;
    for(uint j=0;j<unique_structure_bins[i].size();j++){
      prdss.prds_set.push_back(vprds[unique_structure_bins[i][j]]);
      vi_prds_to_remove.push_back(unique_structure_bins[i][j]);
    }
    if(prdss.prds_set.size()){l_derivative_structures_sets.insert(it,prdss);}
  }

  std::sort(vi_prds_to_remove.rbegin(),vi_prds_to_remove.rend()); //sort in reverse for removal

  for(uint i=0;i<vi_prds_to_remove.size();i++){(*it_prdss).prds_set.erase((*it_prdss).prds_set.begin()+vi_prds_to_remove[i]);}  //remove!

  return unique_structure_bins.size();
}

bool POccCalculator::calculate(){
  string soliloquy = "POccCalculator::calculate():";
  stringstream message;

  if(!p_str.calculateNHNF()){return false;}               //get n_hnf

  //starting criteria for HNF matrices
  //xmatrix<double> hnf_mat;                            //really xmatrix of int's, but we rule for int * double in xmatrix, no big deal
  //vector<xmatrix<double> > v_unique_superlattices;    //only store locally, as we need to make comparisons

  //starting criteria for site combinations
  //vector<vector<int> > v_types_config;
  //vector<int> v_config_iterators, v_site_iterators;
  
  //debug
  /*for(uint i=0;i<vv_count_configs.size();i++){vv_count_configs[i].clear();}
  vv_count_configs.clear();
  vv_count_configs.push_back(vector<POccSiteConfiguration>(0));
  vv_count_configs.push_back(vector<POccSiteConfiguration>(0));
  POccSiteConfiguration test;
  test.types_configuration_debug.push_back(-1);
  test.types_configuration_debug.push_back(-1);
  test.types_configuration_debug.push_back(0);
  test.types_configuration_debug.push_back(0);
  vv_count_configs[0].push_back(test);
  test.types_configuration_debug.clear();
  test.types_configuration_debug.push_back(1);
  test.types_configuration_debug.push_back(1);
  test.types_configuration_debug.push_back(1);
  test.types_configuration_debug.push_back(2);
  vv_count_configs[0].push_back(test);
  test.types_configuration_debug.clear();
  test.types_configuration_debug.push_back(3);
  test.types_configuration_debug.push_back(3);
  test.types_configuration_debug.push_back(4);
  test.types_configuration_debug.push_back(4);
  vv_count_configs[1].push_back(test);
  test.types_configuration_debug.clear();
  test.types_configuration_debug.push_back(5);
  test.types_configuration_debug.push_back(6);
  test.types_configuration_debug.push_back(7);
  test.types_configuration_debug.push_back(8);
  vv_count_configs[1].push_back(test);
  for(uint i=0;i<vv_count_configs.size();i++){
    for(uint j=0;j<vv_count_configs[i].size();j++){
    cerr << "site " << i << " config " << j <<  endl;
      for(uint k=0;k<vv_count_configs[i][j].types_configuration_debug.size();k++){
        cerr << vv_count_configs[i][j].types_configuration_debug[k] << " ";
      }
      cerr << endl;
    }
  } */
  //exit(0);
  
  getTotalPermutationsCount();

  //POccSuperStructure ps_str;
  //const xstructure& xstr_nopocc = p_str.xstr_nopocc;

  //xmatrix<double> distance_matrix;
  //vector<vector<uint> > v_bonded_atom_indices, v_nonbonded_atom_indices; 

  //BIG REVELATION: nearest neighbor distances change with vacancies! 
  //this is a problem, because we need to recalculate bonding for every new site configuration
  //originally, we iterated by hnf_mat first, then site configurations
  //because they are completely independent, we can reverse it
  //iterate each new site configuration first, then hnf_mat's
  //there are some obvious speed improvements to consider
  //if no vacancies, only perform bonding analysis once
  //if no vacancies, we can extrapolate nn distances from primitive cell
  //need to set two booleans, one for if vacancies, and if bonding set
  //if vacancies, bonding set won't matter, as we redo everytime
  //remember, wait to do minkowski reduction until the end
  //it does reduce the size of the cluster, but it also changes the number of atoms (BAD for v_types_config)
  //the cluster only needs to be calculated once, from there we can determine each new bonding
  //make functions for getting nn distances and then determining bonding, which simply returns if no vacancies and already set

  message << "Calculating unique derivative supercells. Please be patient";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

  double energy_radius=RadiusSphereLattice(p_str.getLattice())*1.5; //figure this out corey, only works if energy radius = 10
  if(SET_KESONG_STANDARD_DIST){energy_radius=ENERGY_RADIUS;}

  //energy_radius=12;

  ps_str.initialize(p_str,energy_radius);
  
  unsigned long long int hnf_index=0;
  unsigned long long int site_config_index;
  unsigned long long int current_iteration=0;
  
  //NEW
  POccReducedDerivativeStructure prds;
  resetHNFMatrices();
	updateProgressBar(current_iteration,total_permutations_count,*p_oss);
  while(getNextHNFMatrix()){
    ps_str.getCluster(hnf_mat);
    prds.hnf_index=hnf_index;
    resetSiteConfigurations();
    site_config_index=0;
    while(getNextSiteConfiguration()){
      prds.site_config_index=site_config_index;
      prds.energy=getEnergy();
      add2DerivativeStructuresList(prds);
      site_config_index++;
			updateProgressBar(++current_iteration,total_permutations_count,*p_oss);
    }
    hnf_index++;
  }
  
  //ORIGINAL
  //POccDerivativeStructure pds;
  //pds.reset();
  //resetHNFMatrices();
  ////cerr << "HERE2" << endl;
  //while(getNextHNFMatrix()){
  //  //cerr << "HERE3" << endl;
  //  ps_str.getCluster(hnf_mat);
  //  pds.hnf_index=hnf_index;
  //  //cerr << "HERE4" << endl;
  //  pds.hnf_mat=hnf_mat;
  //  resetSiteConfigurations();
  //  site_config_index=0;
  //  //cerr << "HERE5" << endl;
  //  while(getNextSiteConfiguration()){
  //    //cerr << "HERE6" << endl;
  //    pds.site_config_index=site_config_index;
  //    //cerr << pds.site_config_index << endl;
  //    pds.v_types_config=v_types_config;
  //    //cerr << "HERE7" << endl;
  //    getBonds();
  //    //cerr << "HERE8" << endl;
  //    pds.energy=getEnergy();
  //    pds.degeneracy=1;
  //    //cerr << "HERE9" << endl;
  //    add2DerivativeStructuresList(pds);
  //    //cerr << "DONE ADDING " << endl;
  //    site_config_index++;
	//		updateProgressBar(++current_iteration,total_permutations_count,*p_oss);
  //    //cerr << "HERE10" << endl;
  //  }
  //  hnf_index++;
  //}
  
  //tests of stupidity - START
  if(!l_derivative_structures_sets.size()){
    message << "No derivative structures detected by UFF. Please run calculate() to determine unique derivative structures";
    throw AFLOWLogicError(soliloquy,message);
    //pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
    //cerr << "probably broken pointers" << endl;
    //exit(1);
  }
  

  if(PERFORM_ROBUST_STRUCTURE_COMPARISON){
    message << "Calculated " << l_derivative_structures_sets.size() << " unique derivative structure groups via UFF, now performing robust structure comparison";
    pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);//_LOGGER_COMPLETE_);
    
    unsigned long long int count_new_groups;
    for(std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();it!=l_derivative_structures_sets.end();++it){
      POccReducedDerivativeStructuresSet& prdss=(*it);
      std::sort(prdss.prds_set.begin(),prdss.prds_set.end()); //sort by energies (should all be about the same, ensures splits via struct compar stay sorted by energy
      count_new_groups=runRobustStructureComparison(it);  //l_derivative_structures_sets may change size, get count of the increase
      std::advance(it,count_new_groups-1);  //adjust iterator appropriately
    }
  }
  
  message << "Calculated " << l_derivative_structures_sets.size() << " unique derivative supercells";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_COMPLETE_);
  
  unsigned long long int total_degeneracy=0;
  
  stringstream message_prec;
  for(std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();it!=l_derivative_structures_sets.end();++it){
    //cerr << "UNIQUE " << std::distance(l_derivative_structures_sets.begin(),it) << endl;
    message_prec << std::fixed << std::setprecision(15) << "Unique structure "  << std::distance(l_derivative_structures_sets.begin(),it)+1 << ": energy=" << (*it).getEnergy() << ", " << "degeneracy=" << (*it).getDegeneracy();
    pflow::logger(soliloquy,message_prec,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    //cerr << endl;
    total_degeneracy+=(*it).getDegeneracy();
  }
  
  if(total_permutations_count!=total_degeneracy){
    throw AFLOWLogicError(soliloquy,"Unexpected degeneracy count (does not match expected total permutations count)");
  }

  if(WRITE_OUT_ALL_STRUCTURES){
    message << "Writing out " << ALL_DERIVATIVE_STRUCTURES_FILE << ". Please be patient.";
    pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    stringstream all_derivative_structures_ss;
    POccDerivativeStructure pds;
    for(std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();it!=l_derivative_structures_sets.end();++it){
      const POccReducedDerivativeStructuresSet& prdss=(*it);
      //all_derivative_structures_ss << "---------------------------------------------------------------------------------------------" << endl;
      all_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
      all_derivative_structures_ss << AFLOW_POCC_TAG << "STRUCTURE GROUP #" << std::distance(l_derivative_structures_sets.begin(),it)+1 << endl;
      all_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
      for(uint i=0;i<prdss.prds_set.size();i++){
        pds=getParticularPOccDerivativeStructure(prdss.prds_set[i]);
        pds.degeneracy=1;
        all_derivative_structures_ss << AFLOW_POCC_TAG << "UFF_ENERGY=" << std::fixed << std::setprecision(15) << prdss.prds_set[i].energy << endl;
        all_derivative_structures_ss.unsetf(std::ios_base::floatfield);
        all_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
        all_derivative_structures_ss << ps_str.createDerivativeStructure(pds,p_str.n_hnf,hnf_count,types_config_permutations_count,true,PRIMITIVIZE); // << endl;
        all_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
      }
      all_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
    }
    aurostd::stringstream2file(all_derivative_structures_ss,ALL_DERIVATIVE_STRUCTURES_FILE);
  }

  //tests of stupidity - END

  return true;
}

unsigned long long int POccCalculator::getUniqueDerivativeStructuresCount() const {
  return l_derivative_structures_sets.size();
}

int getZeroPadding(unsigned long long int num){return int(log10(num))+1;}

string POccCalculator::getARUNString(unsigned long long int i) {
  stringstream ARUN;

  std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();
  std::advance(it,i);


  ARUN << "ARUN.POCC_";
  ARUN << std::setfill('0') << std::setw(getZeroPadding(l_derivative_structures_sets.size())) << i+1 << "_";  //+1 so we start at 1, not 0 (count)
  ARUN << "H";
  ARUN << std::setfill('0') << std::setw(getZeroPadding(hnf_count-1)) << (*it).getPRDS().hnf_index; //-1 because this is treated as a HASH and it starts at 0
  ARUN << "C";
  ARUN << std::setfill('0') << std::setw(getZeroPadding(types_config_permutations_count-1)) << (*it).getPRDS().site_config_index; //-1 because this is treated as a HASH and it starts at 0
  
  //string ARUN;
  //ARUN="ARUN.POCC"+aurostd::utype2string(i);
  //ARUN+="H"+aurostd::utype2string((*it).getPRDS().hnf_index);//+1);
  //ARUN+="C"+aurostd::utype2string((*it).getPRDS().site_config_index);//+1);
  //return ARUN;

  return ARUN.str();
}

xstructure POccCalculator::getUniqueDerivativeStructure(unsigned long long int i) {
  string soliloquy = "POccCalculator::getUniqueDerivativeStructure():";
  stringstream message;

  if(i>l_derivative_structures_sets.size()-1){
    throw AFLOWLogicError(soliloquy,"Invalid derivative structure index");
  }

  std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();
  std::advance(it,i);

  POccDerivativeStructure pds=getParticularPOccDerivativeStructure((*it).getPRDS());
  pds.degeneracy=(*it).getDegeneracy();
	
  return ps_str.createDerivativeStructure(pds,p_str.n_hnf,hnf_count,types_config_permutations_count,true,PRIMITIVIZE);
}

vector<xstructure> POccCalculator::getUniqueDerivativeStructures() {
  string soliloquy = "POccCalculator::getUniqueDerivativeStructures():";
  stringstream message;

  vector<xstructure> v_xstr;
  if(!l_derivative_structures_sets.size()){
    message << "No derivative structures detected. Please run calculate() to determine unique derivative structures";
    pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
    return v_xstr;
  }

  unsigned long long int current_iteration=0;
  POccDerivativeStructure pds;
	updateProgressBar(current_iteration,l_derivative_structures_sets.size()-1,*p_oss);
  for(std::list<POccReducedDerivativeStructuresSet>::iterator it=l_derivative_structures_sets.begin();it!=l_derivative_structures_sets.end();++it){
    pds=getParticularPOccDerivativeStructure((*it).getPRDS());
    pds.degeneracy=(*it).getDegeneracy();
    v_xstr.push_back(ps_str.createDerivativeStructure(pds,p_str.n_hnf,hnf_count,types_config_permutations_count,true,PRIMITIVIZE));
		updateProgressBar(++current_iteration,l_derivative_structures_sets.size()-1,*p_oss);
    //cout << AFLOWIN_SEPARATION_LINE << endl;
    //cout << ps_str.createDerivativeStructure((*it),true);
    //cout << AFLOWIN_SEPARATION_LINE << endl;
  }
  
  message << "Full list generated";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

  return v_xstr;

//  resetHNFMatrices();
//  while(getNextHNFMatrix()){
  ////unsigned long long int total_permutations_count=0;
  //  //cerr << "NEW HNF" << endl;
  //  //cerr << hnf_mat << endl;
  //  
  //  v_types_config.clear();
  //  
    
    //getBonding(hnf_mat,distance_matrix,v_bonded_atom_indices,v_nonbonded_atom_indices);
    //ps_str.getBonding(p_str,hnf_mat);
//    getBonds();

//    resetSiteConfigurations();
//    while(getNextSiteConfiguration()){
      //getSpecificBonding();
//      getEnergy();
      //getXStructure();
  //    //total_permutations_count++;
  //    //cerr << "Possibility " << total_permutations_count << "    ";
  //    //for(uint i=0;i<v_types_config.size();i++){
  //    //  for(uint j=0;j<v_types_config[i].size();j++){
  //    //    cerr << v_types_config[i][j] << " ";
  //    //  }
  //    //  cerr << "    ";
  //    //}
  //    //cerr << endl;
//    }
    //exit(0);
  //
  //  //cerr << "end hnf" << endl << endl;
//  }



  //ostream& oss=cout;
  //message << "Total number of derivative structure possibilities = " << total_permutations_count;
  //pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

  //exit(0);

  //return true;
}

} // namespace pocc

namespace pocc {
//--------------------------------------------------------------------------------
// class POccSiteConfiguration (nested in POccCalculator)
//--------------------------------------------------------------------------------
POccSiteConfiguration::POccSiteConfiguration() {
  free();
  //v_pocc_groups= new vector<POccGroup>(0);    //safe
}

POccSiteConfiguration::POccSiteConfiguration(int _site,int _i_hnf,vector<POccGroup>& _v_pocc_groups) {
  free();
  site=_site;
  i_hnf=_i_hnf;
  v_pocc_groups=_v_pocc_groups;
}

POccSiteConfiguration::POccSiteConfiguration(const POccSiteConfiguration& b){
  copy(b);
}

POccSiteConfiguration::~POccSiteConfiguration() {
  free();
}

void POccSiteConfiguration::free() {
  //v_pocc_groups.clear();
  site=0;
  i_hnf=0;
  v_pocc_groups.clear();
  //cerr << "WOO1" << endl;
  //delete(v_pocc_groups);
  //cerr << "WOO2" << endl;
  //vector<POccGroup> tmp; v_pocc_groups=tmp;
  //delete v_pocc_groups;
  partial_occupation_flag=false;
  xv_occupation_count_input.clear();
  xv_occupation_multiple.clear();
  xv_occupation_count_supercell.clear();
  xv_partial_occupation_value.clear();
  xv_site_error.clear();
  occupation_count_total=0;
  vacancy_count=0;
  max_site_error=0.0;
  //error_total=0.0;
  //types_configuration.clear();
}

void POccSiteConfiguration::copy(const POccSiteConfiguration& b) {
  //v_pocc_groups=b.v_pocc_groups;
  site=b.site;
  i_hnf=b.i_hnf;
  v_pocc_groups=b.v_pocc_groups;
  partial_occupation_flag=b.partial_occupation_flag;
  xv_occupation_count_input=b.xv_occupation_count_input;
  xv_occupation_multiple=b.xv_occupation_multiple;
  xv_occupation_count_supercell=b.xv_occupation_count_supercell;
  xv_partial_occupation_value=b.xv_partial_occupation_value;
  xv_site_error=b.xv_site_error;
  occupation_count_total=b.occupation_count_total;
  vacancy_count=b.vacancy_count;
  //error_total=b.error_total;
  max_site_error=b.max_site_error;
  types_configuration_debug=b.types_configuration_debug;
}

const POccSiteConfiguration& POccSiteConfiguration::operator=(const POccSiteConfiguration& b) { // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

void POccSiteConfiguration::clear() { // clear PRIVATE
  POccSiteConfiguration _temp;
  copy(_temp);
}

void POccSiteConfiguration::prepareNoPOccConfig() {
  //////////////////////////////////////////////////////////////////////////////
  //clear all fields, but save that which is declared in constructor
  int _site=site, _i_hnf=i_hnf;
  vector<POccGroup> _v_pocc_groups=v_pocc_groups;
  clear();
  site=_site; i_hnf=_i_hnf;
  v_pocc_groups=_v_pocc_groups;
  //////////////////////////////////////////////////////////////////////////////
  //v_pocc_groups.push_back(vv_pocc_groups[site][0]);
  partial_occupation_flag = false;
  xvector<int> xvi_dummy(0,0); 
  xvector<double> xvd_dummy(0,0); 
  xvi_dummy(0)=1; xv_occupation_count_input=xvi_dummy;
  xvi_dummy(0)=i_hnf; xv_occupation_multiple=xvi_dummy; xv_occupation_count_supercell=xvi_dummy;
  occupation_count_total=i_hnf;
  xv_site_error=xvd_dummy;
  xv_partial_occupation_value=xvd_dummy;
}

void POccSiteConfiguration::preparePOccConfig() {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="preparePOccConfig():";
  //////////////////////////////////////////////////////////////////////////////
  //clear all fields, but save that which is declared in constructor
  int _site=site, _i_hnf=i_hnf;
  vector<POccGroup> _v_pocc_groups=v_pocc_groups;
  clear();
  site=_site; i_hnf=_i_hnf;
  v_pocc_groups=_v_pocc_groups;
  //////////////////////////////////////////////////////////////////////////////
  partial_occupation_flag = true;
  int pocc_groups_count=v_pocc_groups.size(); //vv_pocc_groups[site].size();
  xvector<int> xvi_dummy(pocc_groups_count-1,0);
  xvector<double> xvd_dummy(pocc_groups_count-1,0);
  //cerr << pocc_groups_count << endl;
  //cerr << xvd_dummy.rows << endl;
  //exit(0);
  xv_occupation_count_input=xvi_dummy;
  xv_occupation_multiple=xvi_dummy;
  xv_occupation_count_supercell=xvi_dummy;
  xv_partial_occupation_value=xvd_dummy;
  xv_site_error=xvd_dummy;
  //create xv of counts of _pocc_group occupants
  //for(int i=0;i<xv_occupation_count_input.rows;i++){xv_occupation_count_input(i)=(*v_pocc_groups)[i].v_occupants.size();cerr << "xv_occupation_count_input(" << i << ")=" << xv_occupation_count_input(i) << endl;}
  for(int i=0;i<xv_occupation_count_input.rows;i++){xv_occupation_count_input[i]=v_pocc_groups[i].v_occupants.size();}
  //debug
  if(LDEBUG){
    cerr << soliloquy << " site=" << site << ", i_hnf=" << i_hnf << endl;
    for(uint i=0;i<(uint)pocc_groups_count;i++){
      cerr << soliloquy << " pocc_group=" << i+1 << "/" << v_pocc_groups.size() << ": ";
      cerr << "v_occupants.size()=" << v_pocc_groups[i].v_occupants.size() << " ?= ";
      cerr << "xv_occupation_count_input=" << xv_occupation_count_input[i] << " ";
      cerr << endl;
    }
  }
  //debug
}

int POccSiteConfiguration::getNextOccupationMultiple(int i_hnf, xvector<int>& xv_next_occupation_multiple) {
  //this is a "dummy" bit-string constructor, except canonical bit-string constructors rely on recursion to enumerate ALL possibilities
  //in general, this is not memory-safe
  //so we enumerator ONE-BY-ONE, with a stop condition in place
  //also, most bit-string constructors start with the last digit, and iterate to the first digit (this is how binary numbers are read)
  //we do the opposite, iterate first digit first
  //it conveniently returns index (starts at 0) of _pocc_group impacted
  //i_hnf=_i_hnf;
  //cerr << "HERE  " << xv_next_occupation_multiple << endl;
  for(int i=0;i<xv_next_occupation_multiple.rows;i++){
    if(xv_next_occupation_multiple(i)>=i_hnf){  //greater than not necessary, but safe
      if(i==xv_next_occupation_multiple.rows-1){break;}  //STOP condition, safe
      xv_next_occupation_multiple(i)=0; continue;
    }
    xv_next_occupation_multiple(i)++;
    //cerr << xv_next_occupation_multiple << endl;
    return i;
  }
  return -1;  //STOP condition, safe
}

int POccSiteConfiguration::calculateOccupationCountTotal(xvector<int>& xv_next_occupation_multiple){
  return aurostd::scalar_product(xv_next_occupation_multiple,xv_occupation_count_input);
}

void POccSiteConfiguration::updateOccupationCounts(int _i_hnf, xvector<int> & xv_next_occupation_multiple){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccSiteConfiguration::updateOccupationCounts():";
  i_hnf=_i_hnf;
  xv_occupation_multiple=xv_next_occupation_multiple;
  occupation_count_total=calculateOccupationCountTotal(xv_occupation_multiple); //aurostd::scalar_product(xv_occupation_multiple,xv_occupation_count_input);
  vacancy_count=i_hnf-occupation_count_total;
  if(LDEBUG){
    cerr << soliloquy << " xv_occupation_multiple=" << xv_occupation_multiple << endl;
    cerr << soliloquy << " xv_occupation_count_input" << xv_occupation_count_input << endl;
    cerr << soliloquy << " i_hnf=" << i_hnf << endl;
    cerr << soliloquy << " occupation_count_total=" << occupation_count_total << endl;
    cerr << soliloquy << " vacancy_count=" << vacancy_count << endl;
  }
}

void POccSiteConfiguration::calculateError() {
  string soliloquy="POccSiteConfiguration::calculateError():";
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //now check max_site_error
  max_site_error=0.0;
  for(uint pocc_group=0;pocc_group<v_pocc_groups.size();pocc_group++){
    xv_occupation_count_supercell[pocc_group]=xv_occupation_multiple[pocc_group]*xv_occupation_count_input[pocc_group];
    xv_partial_occupation_value[pocc_group] = (xv_occupation_multiple[pocc_group]/(double)i_hnf); //0.0;
    if(LDEBUG){
      cerr << soliloquy << " v_pocc_groups[" << pocc_group << "].partial_occupation_value=" << v_pocc_groups[pocc_group].partial_occupation_value << endl;
      cerr << soliloquy << " xv_occupation_count_supercell[" << pocc_group << "]=" << xv_occupation_count_supercell[pocc_group] << endl;
      cerr << soliloquy << " occupation_count_total=" << occupation_count_total << endl;
      cerr << soliloquy << " xv_occupation_count_input[" << pocc_group << "]=" << xv_occupation_count_input[pocc_group] << endl;
    }
    //err = (double)xv_occupation_count_input[pocc_group] * vv_pocc_groups[site][pocc_group].partial_occupation_value;
    xv_site_error[pocc_group] = v_pocc_groups[pocc_group].partial_occupation_value;
    if(LDEBUG){
      cerr << soliloquy << " xv_site_error[pocc_group] = v_pocc_groups[pocc_group].partial_occupation_value = " << v_pocc_groups[pocc_group].partial_occupation_value << endl;
    }
    //if(occupation_count_total > 0){
    //  xv_partial_occupation_value[pocc_group] = ((double)xv_occupation_count_supercell[pocc_group] / (double)i_hnf); //occupation_count_total);
    //  if(LDEBUG){cerr << soliloquy << " xv_partial_occupation_value[pocc_group]  = ((double)xv_occupation_count_supercell[pocc_group] / (double)i_hnf) = " << xv_partial_occupation_value[pocc_group] << endl;}
    //  xv_partial_occupation_value[pocc_group] /= (double)xv_occupation_count_input[pocc_group];  //NORMALIZE to 1.0
    //  if(LDEBUG){cerr << soliloquy << " xv_partial_occupation_value[pocc_group] /= (double)xv_occupation_count_input[pocc_group] = " << xv_partial_occupation_value[pocc_group] << endl;}
    xv_site_error[pocc_group] -= xv_partial_occupation_value[pocc_group];
    if(LDEBUG){cerr << soliloquy << " xv_site_error[pocc_group] -= xv_partial_occupation_value[pocc_group] = " << xv_site_error[pocc_group] << endl;}
    //}
    //err/=(double)xv_occupation_count_input[pocc_group];  //NORMALIZE to 1.0
    //xv_site_error[pocc_group] = abs(xv_site_error[pocc_group]); //stability
    //error_total += xv_site_error[pocc_group];
    //cerr << "GROUP " << i << endl;
    ////cerr << "POCC " << vv_pocc_groups[site][pocc_group].partial_occupation_value << endl;
    //cerr << "MULTIPLE " << xv_occupation_multiple[pocc_group] << endl;
    //cerr << "COUNT " << xv_occupation_count_input[pocc_group] << endl;
    //cerr << "TOTAL " << xv_occupation_count_supercell[pocc_group] << endl;
    //cerr << "ERR1 " << xv_site_error[pocc_group] << endl;
    ////cerr << "ERR " << err << endl;
    //cerr << endl;
  }
  xv_site_error=aurostd::abs(xv_site_error);  //could be negative
  max_site_error = aurostd::max(xv_site_error);
  //cerr << "ERROR " << max_site_error << endl;
  //cerr << endl << endl;
}

//double POccSiteConfiguration::getErrorTotal() const {
//  return error_total;
//}

bool POccSiteConfiguration::isPartiallyOccupied() const {
  return partial_occupation_flag;
}

vector<int> POccSiteConfiguration::getStartingTypesConfiguration() const {
  //return types_configuration_debug; //debug
  vector<int> starting_config(occupation_count_total+vacancy_count);
  //cerr << "WOOOOT " << starting_config.size() << endl;
  //cerr << "occupation_count_total=" << occupation_count_total << endl;
  //cerr << "vacancy_count=" << vacancy_count << endl;
  int site_count=0;
  for(int i=0;i<vacancy_count;i++){starting_config[i]=-1;} site_count+=vacancy_count;
  //for(int i=0;i<vacancy_count;i++){starting_config[site_count]=-1; site_count++;}
  for(uint group=0;group<v_pocc_groups.size();group++){
    for(uint occ=0;occ<v_pocc_groups[group].v_types.size();occ++){
      for(int i=0;i<xv_occupation_multiple[group];i++){
        starting_config[site_count]=v_pocc_groups[group].v_types[occ]; site_count++;
      }
    }
  }
  std::sort(starting_config.begin(),starting_config.end());
  //cerr << "site " << site << ": "; 
  //for(uint i=0;i<starting_config.size();i++){
  //  cerr << starting_config[i] << " ";
  //}
  //cerr << endl;
  //types_configuration=starting_configuration;
  return starting_config;
}

//const vector<POccGroup>& POccSiteConfiguration::getVPOccGroups() const {
//  return v_pocc_groups;
//}

//const vector<int>& POccSiteConfiguration::getTypesConfiguration() const {
//  return types_configuration;
//}

} // namespace pocc

namespace pocc {
//--------------------------------------------------------------------------------
// class POccStructureTemplate
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// constructor
//--------------------------------------------------------------------------------
POccStructureTemplate::POccStructureTemplate() {} //{free();}
POccStructureTemplate::~POccStructureTemplate() {free();}

void POccStructureTemplate::free() {
  xstr_pocc.Clear();
  xstr_nopocc.Clear();
  types2pc_map.clear();
  types2uffparams_map.clear();
}

void POccStructureTemplate::copy(const POccStructureTemplate& b){ // copy PRIVATE
  xstr_pocc=b.xstr_pocc;
  xstr_nopocc=b.xstr_nopocc;
  types2pc_map.clear(); for(uint i=0;i<b.types2pc_map.size();i++){types2pc_map.push_back(b.types2pc_map[i]);}
  types2uffparams_map.clear(); for(uint i=0;i<b.types2uffparams_map.size();i++){types2uffparams_map.push_back(b.types2uffparams_map[i]);}
}

void POccStructureTemplate::setPOccStructure(const xstructure& _xstr_pocc){
  xstr_pocc=_xstr_pocc;
  xvector<double> _stoich_each_type(xstr_pocc.stoich_each_type.size()-1,0);stoich_each_type=_stoich_each_type;
  for(uint i=0;i<xstr_pocc.stoich_each_type.size();i++){stoich_each_type[i]=xstr_pocc.stoich_each_type[i];}
}

void POccStructureTemplate::setNonPOccStructure(const xstructure& _xstr_nopocc){xstr_nopocc=_xstr_nopocc;}
void POccStructureTemplate::setTypes2PcMap(const vector<uint>& _types2pc_map){types2pc_map=_types2pc_map;}
void POccStructureTemplate::setTypes2UFFParamsMap(const vector<UFFParamAtom>& _types2uffparams_map){types2uffparams_map=_types2uffparams_map;}
} // namespace pocc

namespace pocc {
//--------------------------------------------------------------------------------
// class POccStructure
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// constructor
//--------------------------------------------------------------------------------
POccStructure::POccStructure(ostream& oss) : xStream() {initialize(oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,ostream& oss) : xStream() {initialize(xstr_pocc,oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) : xStream() {initialize(xstr_pocc,kflags,oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,kflags,oss);}
POccStructure::POccStructure(ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(FileMESSAGE,oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,FileMESSAGE,oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,FileMESSAGE,oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,kflags,FileMESSAGE,oss);}
POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream() {initialize(xstr_pocc,aflags,kflags,FileMESSAGE,oss);}
POccStructure::POccStructure(const POccStructure& b) {copy(b);} // copy PUBLIC

POccStructure::~POccStructure() {freeAll();}

const POccStructure& POccStructure::operator=(const POccStructure& b) {  // operator= PUBLIC
  if(this!=&b) {freeAll();copy(b);}
  return *this;
}

void POccStructure::clear() {POccStructure a;copy(a);}  // clear PRIVATE
void POccStructure::free() {
  POccStructureTemplate::free();
  m_initialized=false;
  m_aflags.clear();
  m_kflags.clear();
  xstr_sym.Clear();
  n_hnf=0;
  v_sites.clear();
  pocc_atoms_total=0;
  for(uint i=0;i<vv_pocc_groups.size();i++){vv_pocc_groups[i].clear();} vv_pocc_groups.clear();
  v_str_configs.clear();
  
  //fix this!
  hnf_table_general_precision=PRECISION;
  hnf_table_iteration_padding=4;
  hnf_table_error_padding=PRECISION+2;
  hnf_table_column_padding=2*hnf_table_general_precision+2*hnf_table_iteration_padding+15;
  header_max_stoich_error="stoich_error";
  header_max_site_error="site_error";
}

void POccStructure::copy(const POccStructure& b) { // copy PRIVATE
  POccStructureTemplate::copy(b);
  xStream::copyStreams(b);
  m_initialized=b.m_initialized;
  m_aflags=b.m_aflags;
  m_kflags=b.m_kflags;
  xstr_sym=b.xstr_sym;
  n_hnf=b.n_hnf;
  v_sites=b.v_sites;
  pocc_atoms_total=b.pocc_atoms_total;
  vv_pocc_groups.clear(); for(uint i=0;i<vv_pocc_groups.size();i++){vv_pocc_groups[i].clear();}
  for(uint i=0;i<b.vv_pocc_groups.size();i++){vv_pocc_groups.push_back(vector<POccGroup>(0));for(uint j=0;j<b.vv_pocc_groups[i].size();j++){vv_pocc_groups[i].push_back(b.vv_pocc_groups[i][j]);}}
  v_str_configs.clear(); for(uint i=0;i<b.v_str_configs.size();i++){v_str_configs.push_back(b.v_str_configs[i]);}

  hnf_table_general_precision=b.hnf_table_general_precision;
  hnf_table_iteration_padding=b.hnf_table_iteration_padding;
  hnf_table_error_padding=b.hnf_table_error_padding;
  hnf_table_column_padding=b.hnf_table_column_padding;
  header_max_stoich_error=b.header_max_stoich_error;
  header_max_site_error=b.header_max_site_error;
}

bool POccStructure::initialize(ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,aflags,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,kflags,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) {
  freeAll();
  ofstream* _p_FileMESSAGE=new ofstream();m_new_ofstream=true;
  initialize(xstr_pocc,aflags,kflags,*_p_FileMESSAGE,oss);
  m_new_ofstream=true;  //override
  return m_initialized;
}
bool POccStructure::initialize(ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    m_initialized=false;  //no point
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setAFlags(aflags);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setKFlags(kflags);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}
bool POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
  free();
  try{
    setOFStream(FileMESSAGE); m_new_ofstream=false;
    setOSS(oss);
    setAFlags(aflags);
    setKFlags(kflags);
    setPOccStructure(xstr_pocc);
    m_initialized=true;
  }
  catch(AFLOWRuntimeError& re){pflow::logger(re.where(), re.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  catch(AFLOWLogicError& le){pflow::logger(le.where(), le.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  return m_initialized;
}

void POccStructure::setAFlags(const _aflags& Aflags) {m_aflags=Aflags;}
void POccStructure::setKFlags(const _kflags& Kflags) {m_kflags=Kflags;}

void POccStructure::preparePOccStructure() {
  string soliloquy="POccStructure::preparePOccStructure():";
  stringstream message;
  
  if(!xstr_pocc.partial_occupation_flag){
    message << "No partially occupied sites detected in this structure";
    throw AFLOWLogicError(soliloquy,message);
  }
  
  //fix up structure
  pflow::fixEmptyAtomNames(xstr_pocc, false); //true);  //aflow likes pp info
  xstr_pocc.ReScale(1.0);
  xstr_pocc.ShifOriginToAtom(0);
  xstr_pocc.BringInCell();

  //get types2pc_map
  types2pc_map.clear();
  uint atom_index=0;
  for(uint i=0;i<xstr_pocc.num_each_type.size();i++){
    types2pc_map.push_back(atom_index);
    atom_index+=xstr_pocc.num_each_type[i];
  }
  //WRONG IMPOSSIBLE - get types2pc_map, assumes xstr_pocc can have types out of order
  /*for(uint i=0;i<xstr_pocc.species.size();i++){types2pc_map.push_back(0);}
  bool found_type;
  for(int i=0;i<(int)types2pc_map.size();i++){
    found_type=false;
    for(uint j=0;j<xstr_pocc.atoms.size() && !found_type;j++){
      if(xstr_pocc.atoms[j].type==i){
        types2pc_map[i]=j;
        //remove pocc stuff - START
        //types2pc_map[i].partial_occupation_value=1.0;
        //types2pc_map[i].partial_occupation_flag=false;
        //remove pocc stuff - END
        found_type=true;
      }
    }
  }*/

  //quick test of stupidity - START
  //see if types are all unique
  for(uint i=0;i<types2pc_map.size();i++){
    for(uint j=i+1;j<types2pc_map.size();j++){
      if(xstr_pocc.atoms[types2pc_map[i]].name==xstr_pocc.atoms[types2pc_map[j]].name){
        if(xstr_pocc.atoms[types2pc_map[i]].name==""){
          message << "Bad input - please assign names to all inputs. Results vary with different atom types.";
        }else{
          message << "Bad input - atoms " << types2pc_map[i]+1 << "/" << xstr_pocc.atoms.size() << " ";
          message << "and " << types2pc_map[j]+1 << "/" << xstr_pocc.atoms.size() << " ";
          message << "show to have the same name (" << xstr_pocc.atoms[types2pc_map[i]].name << "). ";
          message << "Please reorganize your input by type.";
        }
        throw AFLOWLogicError(soliloquy,message);
      }
    }
  }
  //check that code can handle arbitrary reassignment of random elements
  vector<string> elements=getElementsList();
  bool use_mixed_set_elements=true;
  if(xstr_pocc.num_each_type.size()>elements.size()){
    use_mixed_set_elements=false;
    message << "Unexpectedly large count of species in structure (nspecies=" << xstr_pocc.num_each_type.size();
    message <<  ">" << elements.size()<< "), resorting to exact UFF energy determination";
    pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_);
  }
  //quick test of stupidity - END
  
  //exit(0);

  //DO UFF STUFF IF IT'S THE RIGHT MODE
  //MAKE DIFFERENT MODES HERE
  string element,uff_parameters_string;
  vector<string> uff_params;
  UFFParamAtom uffp;
  uint elements_index=0;
  for(uint i=0;i<types2pc_map.size();i++){
    if(COMPARE_WITH_KESONG||use_mixed_set_elements==false){element=xstr_pocc.atoms[types2pc_map[i]].cleanname;} //exact UFF energies
    else{element=elements[elements_index++];}                                                                   //reassigned UFF energies for orthogonality of species properties
    if(!getUFFParameters(element,uff_parameters_string)){throw AFLOWLogicError(soliloquy,"Unable to fetch UFF parameters (requested "+element+")");}

    aurostd::string2tokens(uff_parameters_string,uff_params," ");
    types2uffparams_map.push_back(uffp);
    types2uffparams_map.back().symbol=uff_params[0];
    //cerr << "i=" << i << " " << uff_params[0] << endl;
    types2uffparams_map.back().r1=aurostd::string2utype<double>(uff_params[1]);
    types2uffparams_map.back().theta0=aurostd::string2utype<double>(uff_params[2]);
    types2uffparams_map.back().x1=aurostd::string2utype<double>(uff_params[3]);
    types2uffparams_map.back().D1=aurostd::string2utype<double>(uff_params[4]);
    types2uffparams_map.back().zeta=aurostd::string2utype<double>(uff_params[5]);
    types2uffparams_map.back().Z1=aurostd::string2utype<double>(uff_params[6]);
    types2uffparams_map.back().Vi=aurostd::string2utype<double>(uff_params[7]);
    types2uffparams_map.back().Uj=aurostd::string2utype<double>(uff_params[8]);
    types2uffparams_map.back().ChiI=aurostd::string2utype<double>(uff_params[9]);
    types2uffparams_map.back().hard=aurostd::string2utype<double>(uff_params[10]);
    types2uffparams_map.back().radius=aurostd::string2utype<double>(uff_params[11]);
  }
  //exit(0);

  partitionPOccSites(); //get v_site

  //get non-pocc structure and calculate its symmetry (pgroups)
  calculateNonPOccStructure();
  calculateSymNonPOccStructure();

}

const xmatrix<double>& POccStructure::getLattice() const {return xstr_sym.lattice;}
const vector<_sym_op>& POccStructure::getPGroup() const {return xstr_sym.pgroup;}

void POccStructure::partitionPOccSites() {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccStructure::partitionPOccSites():";
  //go through atoms in PARTCAR, grab site_indices
  
  //reset
  v_sites.clear();
  for(uint i=0;i<vv_pocc_groups.size();i++){vv_pocc_groups[i].clear();}
  vv_pocc_groups.clear();

  bool new_site;
  POccGroup dummy_pg;

  //occupy v_sites first
  for(uint i=0;i<xstr_pocc.atoms.size();i++){
    new_site=true;
    for(uint j=0;j<v_sites.size() && new_site;j++){
      if(v_sites[j].v_occupants.size() && aurostd::modulus(xstr_pocc.atoms[i].cpos-xstr_pocc.atoms[v_sites[j].v_occupants[0]].cpos)<ZERO_TOL){
        v_sites[j].v_occupants.push_back(i);
        v_sites[j].v_types.push_back(xstr_pocc.atoms[i].type);
        new_site=false;
      }
    }
    if(new_site){
      v_sites.push_back(POccSite());
      v_sites.back().site=v_sites.size()-1;
      v_sites.back().v_occupants.push_back(i);
      v_sites.back().v_types.push_back(xstr_pocc.atoms[i].type);
      v_sites.back().partial_occupation_flag=!aurostd::isequal(xstr_pocc.atoms[i].partial_occupation_value,1.0,ZERO_TOL); //if the only other occupant is a vacancy
    }
  }
  pocc_atoms_total=0;
  for(uint i=0;i<v_sites.size();i++){
    v_sites[i].partial_occupation_flag=bool(v_sites[i].partial_occupation_flag||v_sites[i].v_occupants.size()>1); //if the only other occupant is a vacancy
    if(v_sites[i].partial_occupation_flag){pocc_atoms_total+=v_sites[i].v_occupants.size();}
  } //is it a pocc site?

  //occupy vv_pocc_groups next
  for(uint site=0;site<v_sites.size();site++){
    vv_pocc_groups.push_back(vector<POccGroup>(0));
    for(uint occ=0;occ<v_sites[site].v_occupants.size();occ++){
      new_site=true;
      for(uint k=0;k<vv_pocc_groups[site].size() && new_site;k++){
        if(aurostd::isequal(xstr_pocc.atoms[v_sites[site].v_occupants[occ]].partial_occupation_value,\
              vv_pocc_groups[site][k].partial_occupation_value,ZERO_TOL)){
          //cerr << xstr_pocc.atoms[v_sites[site].v_occupants[occ]].partial_occupation_value << " vs. " << vv_pocc_groups[site][k].partial_occupation_value << endl;
          vv_pocc_groups[site][k].v_occupants.push_back(v_sites[site].v_occupants[occ]);
          vv_pocc_groups[site][k].v_types.push_back(v_sites[site].v_types[occ]);
          new_site=false;
        }
      }
      if(new_site){
        vv_pocc_groups[site].push_back(dummy_pg);
        vv_pocc_groups[site].back().site=site;
        vv_pocc_groups[site].back().partial_occupation_value=xstr_pocc.atoms[v_sites[site].v_occupants[occ]].partial_occupation_value;
        vv_pocc_groups[site].back().v_occupants.push_back(v_sites[site].v_occupants[occ]);
        vv_pocc_groups[site].back().v_types.push_back(v_sites[site].v_types[occ]);
      }
    }
    //std::sort(vv_pocc_groups[site].rbegin(),vv_pocc_groups[site].rend()); //NOTE the r, reverse sort, bigger poccs first
  }

  if(LDEBUG){
    for(uint i=0;i<vv_pocc_groups.size();i++){
      cerr << soliloquy << " v_sites[site=" << i << "]: ";
      cerr << "partial_occupation_flag=" << v_sites[i].partial_occupation_flag << endl;
      for(uint j=0;j<vv_pocc_groups[i].size();j++){
        cerr << soliloquy << " vv_pocc_groups[site=" << i << " group=" << j << "]: ";
        cerr << "partial_occupation_value=" << vv_pocc_groups[i][j].partial_occupation_value << ", ";
        cerr << "v_occupants.size()=" << vv_pocc_groups[i][j].v_occupants.size() << " ";
        cerr << endl;
      }
    }
  }
}

void POccStructure::getOptimizedSiteCountConfigurations(
    int site,
    int i_hnf,
    vector<POccSiteConfiguration>& v_site_configs) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccStructure::getOptimizedSiteCountConfigurations():";
  stringstream message;
  //reset
  v_site_configs.clear();
  
  //trivial site_config, no pocc
  POccSiteConfiguration site_config(site,i_hnf,vv_pocc_groups[site]); //sc_dummy;

  //test of stupidity
  if(vv_pocc_groups[site].size()<1){
    message << "Invalid count of POCC site groups (vv_pocc_groups[" << site << "].size()=" << vv_pocc_groups[site].size() << ")";
    throw AFLOWLogicError(soliloquy,message);
  }
  
  //if it's truly non-pocc, then there's only one vv_pocc_groups[site]; 
  //otherwise, there's more than 1 and it doesn't matter
  site_config.prepareNoPOccConfig();
  
  if(!v_sites[site].partial_occupation_flag){v_site_configs.push_back(site_config);return;}  //no pocc, trivial
  
  //otherwise, pocc!
  site_config.preparePOccConfig();
  //int pocc_groups_count=vv_pocc_groups[site].size();
  //already initialized to 0, xvector<uint> is ill-defined (in AUROSTD), so just use int, it's okay
  //aurostd::xvector<int> xv_occupation_count_input(pocc_groups_count);
  //aurostd::xvector<int> xv_occupation_multiple(pocc_groups_count), xv_occupation_multiple_last(pocc_groups_count);


  //this is the stop condition
  //for(int i=0;i<site_config.xv_occupation_multiple_last.rows;i++){site_config.xv_occupation_multiple_last(i)=(int)i_hnf;}

  //cerr << "START, n=" << pocc_groups_count << ",k=" << i_hnf << endl;
  int pocc_group_iterator=0;    //current group impacted
  int pocc_group_to_skip=-1;    //slight optimizer, if occupation_number hits i_hnf early, skip remaining possibilities
  bool skipping;                //needed for resetting pocc_group_to_skip
  
  double eps=1.0;
  xvector<int> xv_next_occupation_multiple=site_config.xv_occupation_multiple;
  //int minimum_vacancy_count=i_hnf;
  while(pocc_group_iterator!=-1){
    if(LDEBUG){
      cerr << soliloquy << " looking at xv_next_occupation_multiple=" << xv_next_occupation_multiple << endl;
      cerr << soliloquy << " pocc_group_to_skip=" << pocc_group_to_skip << " ?= pocc_group_iterator=" << pocc_group_iterator << endl;
    }
    if(pocc_group_to_skip!=pocc_group_iterator){  //slight optimizer
      if(LDEBUG){
        cerr << soliloquy << " exploring xv_next_occupation_multiple=" << xv_next_occupation_multiple << endl;
      }
  //    cerr << "calculateOccupationCountTotal(xv_next_occupation_multiple)=" << site_config.calculateOccupationCountTotal(xv_next_occupation_multiple) << endl;
      
      if(site_config.calculateOccupationCountTotal(xv_next_occupation_multiple)>i_hnf){ //start ignoring
        if(LDEBUG){
          cerr << soliloquy << " will skip multiples from now on" << endl;
          cerr << soliloquy << " site_config.calculateOccupationCountTotal(xv_next_occupation_multiple)=" << site_config.calculateOccupationCountTotal(xv_next_occupation_multiple) << " > " << i_hnf << "=i_hnf" << endl;
        }
        pocc_group_to_skip=pocc_group_iterator; skipping=true; continue;
      }
      
      site_config.updateOccupationCounts(i_hnf,xv_next_occupation_multiple);
      
      //found viable configuration
      site_config.calculateError();

  //    cerr << "eps=" << eps << endl;
  //    cerr << "site_config.max_site_error=" << site_config.max_site_error << endl;
      if(aurostd::isequal(site_config.max_site_error,eps,ZERO_TOL) || site_config.max_site_error<eps){
        if(!aurostd::isequal(site_config.max_site_error,eps,ZERO_TOL)){
          eps=site_config.max_site_error;
          v_site_configs.clear();
        }
  //      cerr << "adding one!!!!!!1" << endl;
        v_site_configs.push_back(site_config);
        //minimum_vacancy_count=min(minimum_vacancy_count,site_config.vacancy_count);
      }
    }
    //getNextOccupationMultiple() is basically a bitstring iterator:  000, 001, 010, ...
    
    //it will increment xv_occupation_multiple and return the index of the digit that was modified
    //a return value of -1 means it's exhausted all possibilities 
    pocc_group_iterator = site_config.getNextOccupationMultiple(i_hnf,xv_next_occupation_multiple);
    if(skipping && pocc_group_iterator!=pocc_group_to_skip){pocc_group_to_skip=-1; skipping=false;}
  //  cerr << "i_hnf=" << i_hnf << endl;
  //  cerr << "xv_next_occupation_multiple=" << xv_next_occupation_multiple << endl;
  //  cerr << "pocc_group_iterator=" << pocc_group_iterator << endl;
  //  cerr << endl;
  }
//  for(uint i=0;i<v_site_configs.size();i++){
//    cerr << "site config " << i << endl;
//    cerr << "v_site_configs[i].xv_occupation_count_input = " << v_site_configs[i].xv_occupation_count_input << endl;
//    cerr << "v_site_configs[i].xv_occupation_multiple = " << v_site_configs[i].xv_occupation_multiple << endl;
//    cerr << "v_site_configs[i].xv_occupation_count_supercell = " << v_site_configs[i].xv_occupation_count_supercell << endl;
//    cerr << "v_site_configs[i].occupation_count_total = " << v_site_configs[i].occupation_count_total << endl;
//    cerr << "v_site_configs[i].vacancy_count = " << v_site_configs[i].vacancy_count << endl;
//    cerr << "v_site_configs[i].max_site_error = " << v_site_configs[i].max_site_error << endl;
//  }
//  cerr << "END" << endl;
  //
  //if(REDUCE_VACANCIES){
  //  if(v_site_configs.size()){
  //    for(uint config=v_site_configs.size()-1;config-- > 0; ){  //will stop at 0
  //      if(v_site_configs[config].vacancy_count>minimum_vacancy_count){
  //        v_site_configs.erase(v_site_configs.begin()+config);
  //      }
  //    }
  //  }
  //}

}

xvector<double> POccStructure::calculateStoichEachType(vector<vector<int> >& v_types_config) {
  xvector<double> stoich_each_type(types2pc_map.size()-1,0);
  //for(uint type=0;type<types2pc_map.size();type++){stoich_each_type.push_back(0.0);}
  double total_atoms=0.0;
  for(uint site=0;site<v_types_config.size();site++){
    for(uint atom=0;atom<v_types_config[site].size();atom++){
      if(v_types_config[site][atom]!=-1){stoich_each_type[v_types_config[site][atom]]+=1.0;total_atoms+=1.0;}
    }
  }
  if(total_atoms>0.0){
    for(int type=0;type<stoich_each_type.rows;type++){stoich_each_type[type]/=total_atoms;}
  }
  return stoich_each_type;
}

//double POccStructure::calculateStoichDiff(deque<double>& s1,deque<double>& s2){
//  if(s1.size()!=s2.size()){exit(1);}
//  double error=0.0;
//  for(uint i=0;i<s1.size();i++){
//    error+=abs(s1[i]-s2[i]);
//  }
//  return error;
//}

//void POccStructure::getSiteCountConfigurations(int i_hnf, double& min_stoich_error) {
void POccStructure::getSiteCountConfigurations(int i_hnf) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccStructure::getSiteCountConfigurations():";
  //return these values, do not necessarily store as "global" values

  //reset
  //for(uint str_config=0;str_config<v_str_configs.size();str_config++){
  //  v_str_configs[str_config].clear();
  //}
  v_str_configs.clear();
  //v_site_occupancy.clear();
  //v_site_effective_pocc.clear();
  //v_site_error.clear();
  //v_vacancy_count.clear();

  //only useful here, but super important!
  //look at case i_hnf = 5, and a site with 8 occupants, evenly occupied (easy)
  //the original algorithm (look at each pocc atom individually), 
  //would want to fill this site with one atom (+8 atoms total, not even space!)
  //so we group by occupation as well, and allowed occupations are when all equally occupants can fill together

  //create space
  //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {v_site_error.push_back(1.0);}
  //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {v_site_effective_pocc.push_back(0.0);}
  //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {v_site_occupancy.push_back(0);}
  //for(uint i=0;i<v_sites.size();i++){v_vacancy_count.push_back(i_hnf);}
    
  //perform calculation
  //cerr << "HERE4" << endl;
  vector<vector<POccSiteConfiguration> > vv_count_configs;  //first layer site, second layer config (opposite from v_str_configs)
  vector<POccSiteConfiguration> v_site_configs;
  for(uint site=0;site<vv_pocc_groups.size();site++){
    //cerr << "HERE6" << endl;
    //cerr << "site=" << site <<endl;
    //cerr << "i_hnf=" << i_hnf <<endl;
    //cerr << "v_site_configs.size()=" << v_site_configs.size() <<endl;
    getOptimizedSiteCountConfigurations(site,i_hnf,v_site_configs);
  //  cerr << "site = " << site << " count: " << v_site_configs.size() << endl;
    vv_count_configs.push_back(v_site_configs);
    //cerr << "HERE8" << endl;
    //cerr << "site = " << site << endl;
  }
  if(LDEBUG){
  for(uint site=0;site<vv_count_configs.size();site++){
    cerr << soliloquy << " site " << site << endl;
    for(uint config=0;config<vv_count_configs[site].size();config++){
      cerr << soliloquy << " config " << config << endl;    //
      cerr << soliloquy << " vv_count_configs[" << site << "][" << config << "].xv_site_error=" << vv_count_configs[site][config].xv_site_error << endl;
      cerr << soliloquy << " vv_count_configs[" << site << "][" << config << "].max_site_error=" << vv_count_configs[site][config].max_site_error << endl;
      cerr << soliloquy << " vv_count_configs[" << site << "][" << config << "].xv_occupation_count_input=" << vv_count_configs[site][config].xv_occupation_count_input << endl;
      cerr << soliloquy << " vv_count_configs[" << site << "][" << config << "].xv_occupation_multiple=" << vv_count_configs[site][config].xv_occupation_multiple << endl;
      cerr << soliloquy << " vv_count_configs[" << site << "][" << config << "].xv_partial_occupation_value=" << vv_count_configs[site][config].xv_partial_occupation_value << endl;
    }
    cerr << endl << endl;
  }
  }
  //exit(0);
  //we need to fundamentally change scheme here to satisfy STRICT_STOICH_EQUIVALENCE
  //we cannot simply permute through each config on each site
  //only certain configs work together, we must define these
  //vector<uint> v_config_order;    //dummy, don't really care about order just yet
  vector<int> v_config_iterators;
  vector<vector<int> > v_types_config;
  xvector<double> stoich_each_type_new;
  double stoich_error;
  double eps=1.0; //(double)vv_pocc_groups.size();   //number of sites
  //while(getNextSiteConfiguration(vv_count_configs,v_config_order,v_config_iterators,v_types_config)){
  //cerr << "COREY START" << endl;
  std::vector<double>::iterator it_max_stoich;
  while(getNextSiteConfiguration(vv_count_configs,v_config_iterators,v_types_config)){
    stoich_each_type_new=calculateStoichEachType(v_types_config);
    stoich_error=aurostd::max(aurostd::abs(stoich_each_type-stoich_each_type_new));
    //stoich_error=calculateStoichDiff(stoich_each_type,xstr_pocc.stoich_each_type);
    //cerr << "TYPES" << endl;
    //for(uint i=0;i<v_types_config.size();i++){
    //  for(uint j=0;j<v_types_config[i].size();j++){
    //    cerr << v_types_config[i][j] << " ";
    //  }
    //  cerr << endl;
    //}
    //cerr << "STOICH" << endl;
    //for(int i=0;i<stoich_each_type_new.rows; i++){
    //  cerr << stoich_each_type_new[i] << endl;
    //}
    //cerr << "STOICH_SET" << endl;
    //for(int i=0;i<stoich_each_type.rows; i++){
    //  cerr << stoich_each_type[i] << endl;
    //}
    //cerr << "ERROR " << stoich_error << " eps" << eps  << endl;
    
    if(aurostd::isequal(stoich_error,eps,ZERO_TOL) || stoich_error<eps){
      //cerr << "FOUND 1" << endl;
      if(!aurostd::isequal(stoich_error,eps,ZERO_TOL)){
        eps=stoich_error;
        //for(uint str_config=0;str_config<v_str_configs.size();str_config++){v_str_configs[str_config].clear();} 
        v_str_configs.clear();
      }
      v_str_configs.push_back(StructureConfiguration());
      v_str_configs.back().max_site_error=0.0;
      for(uint site=0;site<v_config_iterators.size();site++){
        v_str_configs.back().site_configs.push_back(vv_count_configs[site][v_config_iterators[site]]);
        v_str_configs.back().max_site_error=max(v_str_configs.back().max_site_error,v_str_configs.back().site_configs.back().max_site_error);
      }
      //cerr << "SEE " << v_str_configs.back().max_site_error << endl;
      v_str_configs.back().max_stoich_error=stoich_error;
    }
  }
  //cerr << "COREY Stop" << endl;

  //vector<int> starting_config;
  //for(uint str_config=0;str_config<v_str_configs.size();str_config++){
  //  //cerr << "config " << str_config << endl;
  //  for(uint site=0;site<v_str_configs[str_config].size();site++){
  //    //cerr << "site " << site << ":  ";
  //    starting_config=v_str_configs[str_config][site].getStartingTypesConfiguration();
  //    //for(uint i=0;i<starting_config.size();i++){
  //    //  cerr << starting_config[i] << " ";
  //    //}
  //    //cerr << endl;
  //  }
  //}
  //cerr << "HERE" << endl;
  //exit(0);
  
  //perform calculation
  //double eps;
  //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {
  //  v_site_error[iatom]=1.0;
  //  v_site_occupancy[iatom]=0;
  //  for(uint j=0;j<=i_hnf;j++) {
  //    eps=abs((double) j/i_hnf-xstr_pocc.atoms[iatom].partial_occupation_value);
  //    if(eps<v_site_error[iatom]) {
  //      v_site_error[iatom]=eps;
  //      v_site_occupancy[iatom]=j;
  //      v_site_effective_pocc[iatom]=(double) v_site_occupancy[iatom]/i_hnf;
  //    }
  //  }
  //}

  //min_stoich_error=eps;
}

void POccStructure::calculateNonPOccStructure(){

  xstr_nopocc=xstr_pocc;
  //clear structure pocc properties, set defaults
  xstr_nopocc.scale_second=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
  xstr_nopocc.neg_scale_second=FALSE;
  xstr_nopocc.partial_occupation_flag=FALSE;
  xstr_nopocc.partial_occupation_site_tol=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
  xstr_nopocc.partial_occupation_HNF=0;
  xstr_nopocc.partial_occupation_sublattice.clear();

  //make pocc sites non-pocc sites by removing all but one occupant
  //gather indices of atoms to remove, then remove backwards
  vector<uint> v_atoms_to_remove;
  for(uint site=0;site<v_sites.size();site++){
    if(v_sites[site].partial_occupation_flag){
      //test of stupidity
      //[COREY COME BACK, add another reasonable check that includes vacancies]if(v_sites[site].v_occupants.size()<2){return false;}
      //set the 0th index to non-pocc status
      xstr_nopocc.atoms[v_sites[site].v_occupants[0]].partial_occupation_value=1.0;
      xstr_nopocc.atoms[v_sites[site].v_occupants[0]].partial_occupation_flag=false;
      //add others to list to be removed
      for(uint atom=1;atom<v_sites[site].v_occupants.size();atom++){
        v_atoms_to_remove.push_back(v_sites[site].v_occupants[atom]);
      }
    }
  }

  //RemoveAtoms(xstr_nopocc,v_atoms_to_remove);
  xstr_nopocc.RemoveAtoms(v_atoms_to_remove);
  //std::sort(v_atoms_to_remove.rbegin(),v_atoms_to_remove.rend()); //NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
  //for(uint atom=0;atom<v_atoms_to_remove.size();atom++){xstr_nopocc.RemoveAtom(v_atoms_to_remove[atom]);}

  //make copy of symmetry calculation
  xstr_sym=xstr_nopocc;

}

void POccStructure::calculateSymNonPOccStructure() {
  string soliloquy="POccStructure::calculateSymNonPOccStructure():";
  // CO - default should always be to turn everything on, modify input
  //get ROBUST symmetry determination
  pflow::defaultKFlags4SymWrite(m_kflags,true);
  pflow::defaultKFlags4SymCalc(m_kflags,true);
  bool sym_calculated=pflow::PerformFullSymmetry(xstr_sym,*p_FileMESSAGE,m_aflags,m_kflags,true,*p_oss);
  if(sym_calculated==false){throw AFLOWLogicError(soliloquy,"AFLOW-SYM failed to calculate the symmetry of the clean (non-pocc'd) structure");}
}

string POccStructure::hnfTableHeader(){

  stringstream ss_header;

  //HEADER - START
  //fix output to 
  //a) accommodate for fixed hnf (see top)
  //b) when two configs show up, make it apparent!
  //ss_header << aurostd::PaddedPRE(aurostd::utype2string(0),hnf_table_iteration_padding) << "  " << string(hnf_table_error_padding,'-') << " ";
  ss_header << aurostd::PaddedPRE("n",hnf_table_iteration_padding) << "  " << string(hnf_table_error_padding,'-') << " ";
  uint pocc_atom=1;   //start from 1
  //for(uint atom=0;atom<xstr_pocc.atoms.size();atom++) {
  for(uint site=0;site<v_sites.size();site++){
    if(v_sites[site].partial_occupation_flag){
      for(uint occ=0;occ<v_sites[site].v_occupants.size();occ++){
        //output - START
        ss_header << "| ";
        ss_header << aurostd::PaddedCENTER("pocc_atom = " + 
            aurostd::utype2string(pocc_atom) + "/" + 
            aurostd::utype2string(pocc_atoms_total),hnf_table_column_padding+2);  //need to +2, bug in PaddedCENTER()
        ss_header << " " ;  // CO 170629
        //output - END
        pocc_atom++;
      }
    }
  }
  //ss_header << " | " << "total_error" << endl;
  ss_header << " | " << header_max_stoich_error << " | " << header_max_site_error  << endl;

  return ss_header.str();
}

string POccStructure::hnfTableLineOutput(int i_hnf,int str_config){

  string soliloquy="POccStructure::hnfTableLineOutput():";
  stringstream message;
  //TABLE - START
  stringstream line_output;
  stringstream ausf;ausf.precision(hnf_table_general_precision);ausf.setf(std::ios::fixed,std::ios::floatfield);
  //string multiple_configs_partiion="-------------------------------------";
    
  double stoich_error,site_error=1.0/i_hnf;
  POccSiteConfiguration site_config;

  //for(uint str_config=0;str_config<v_str_configs.size();str_config++){
  //output - START
  line_output << aurostd::PaddedPRE(aurostd::utype2string(i_hnf),hnf_table_iteration_padding);
  line_output << "  ";
  ausf.str(""); ausf << site_error;
  line_output << aurostd::PaddedPOST(ausf.str(),hnf_table_error_padding);
  line_output << " ";
  //output - END
  
  stoich_error=v_str_configs[str_config].max_stoich_error;
  site_error=v_str_configs[str_config].max_site_error;
  if(v_str_configs.size()>1){
    //cerr << "COREY CHECK IT OUT" << endl;
    //for(uint str_config=0;str_config<v_str_configs.size();str_config++){
    //  cerr << "config " << str_config << endl;
    //  for(uint site=0;site<v_str_configs[str_config].size();site++){
    //    cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_site_error=" << v_str_configs[str_config].site_configs[site].xv_site_error << endl;
    //    cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].max_site_error=" << v_str_configs[str_config].site_configs[site].max_site_error << endl;
    //    cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_occupation_count_input=" << v_str_configs[str_config].site_configs[site].xv_occupation_count_input << endl;
    //    cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_occupation_multiple=" << v_str_configs[str_config].site_configs[site].xv_occupation_multiple << endl;
    //    cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_partial_occupation_value=" << v_str_configs[str_config].site_configs[site].xv_partial_occupation_value << endl;
    //    cerr << endl << endl;
    //  }
    //  //exit(0);
    //}
  }  //i am interested in finding an example!

  //cerr << "HERE4" << endl;

  //site_error=0.0;
  //print out results
  for(uint site=0;site<v_str_configs[str_config].site_configs.size();site++){
    site_config = v_str_configs[str_config].site_configs[site]; //error output is the same across all configs
    //site_error+=site_config.max_site_error;
    //site_error=max(site_error,site_config.max_site_error);
    if(site_config.isPartiallyOccupied()){
      for(int pocc_group=0;pocc_group<site_config.xv_occupation_count_input.rows;pocc_group++){
        for(int pocc_atom=0;pocc_atom<site_config.xv_occupation_count_input[pocc_group];pocc_atom++){
          //IMPORTANT CHECK
          if(!aurostd::isequal(site_config.xv_occupation_multiple[pocc_group]/((double)i_hnf),
                site_config.xv_partial_occupation_value[pocc_group],_ZERO_TOL_)){
            message << "site_config.xv_occupation_multiple[pocc_group]/((double)i_hnf)=" << site_config.xv_occupation_multiple[pocc_group]/((double)i_hnf);
            message << " != ";
            message << site_config.xv_partial_occupation_value[pocc_group] << "=site_config.xv_partial_occupation_value[pocc_group]";
            throw AFLOWLogicError(soliloquy,message);
          }
          //output - START
          ausf.str("");
          ausf << aurostd::PaddedPRE(aurostd::utype2string(site_config.xv_occupation_multiple[pocc_group]),hnf_table_iteration_padding);
          ausf << "/";
          ausf << aurostd::PaddedPOST(aurostd::utype2string(i_hnf),hnf_table_iteration_padding);
          ausf << " ";
          ausf << site_config.xv_partial_occupation_value[pocc_group];
          ausf << " ";
          ausf << "(err = " << site_config.xv_site_error[pocc_group] << ")";
          line_output << "| " << aurostd::PaddedPOST(ausf.str(),hnf_table_column_padding) << " ";
          //output - END
        }
      }
    }
  }
  //cerr << "stoich_error=" << stoich_error << " site_error=" << site_error << endl;
  //error=(stoich_error+site_error)/2.0;
  //line_output << " | " << ( site_error < ZERO_TOL ? 0.0 : stoich_error ) << endl;
  //line_output << " | " << ( site_error < ZERO_TOL ? 0.0 : site_error ) << endl;
  ausf.str(""); ausf << stoich_error;
  line_output << " | " << aurostd::PaddedPOST(ausf.str(),header_max_stoich_error.size()); //hnf_table_error_padding);
  ausf.str(""); ausf << site_error;
  line_output << " | " << aurostd::PaddedPOST(ausf.str(),header_max_stoich_error.size()); //hnf_table_error_padding);
  line_output << endl;
  //}

  return line_output.str();

}

//void POccStructure::setHNFTablePadding(int PRECISION){
//  //set some nice printing precisions and paddings
//  hnf_table_general_precision=PRECISION;
//  hnf_table_iteration_padding = 4;
//  hnf_table_error_padding=PRECISION+2;
//  hnf_table_column_padding=2*hnf_table_general_precision+2*hnf_table_iteration_padding+15;
//  //header_max_stoich_error="stoich_error";
//  //header_max_site_error="site_error";
//}

bool POccStructure::getHNFTableOutput(int i_hnf,double& stoich_error,double& site_error) {
  
  string soliloquy="POccStructure::getHNFTableOutput()";
  stringstream multiple_configs_ss;
  
  if(!v_str_configs.size()){cerr << "WHOOPS" << endl;return false;}
  
  stoich_error=v_str_configs[0].max_stoich_error; //same for all
  site_error=v_str_configs[0].max_site_error;     //same for all

  //if(v_str_configs.size()>1){message_raw << multiple_configs_partiion << endl;}
  if(v_str_configs.size()>1){
    multiple_configs_ss.str("");
    multiple_configs_ss << "*** MULTIPLE CONFIGURATIONS POSSIBLE FOR HNF=" << i_hnf << " - START ***" << endl;
    pflow::logger(soliloquy,aurostd::PaddedCENTER(multiple_configs_ss.str(),AFLOWIN_SEPARATION_LINE.size()+2),*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
  }

  for(uint str_config=0;str_config<v_str_configs.size();str_config++){
    if(!v_str_configs[str_config].site_configs.size()){cerr << "WHOOPS2" << endl;return false;}
    pflow::logger(soliloquy,hnfTableLineOutput(i_hnf,str_config),*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
  }

  //if(v_str_configs.size()>1){message_raw << multiple_configs_partiion << endl;}
  if(v_str_configs.size()>1){
    multiple_configs_ss.str("");
    multiple_configs_ss << "*** MULTIPLE CONFIGURATIONS POSSIBLE FOR HNF=" << i_hnf << " - END ***" << endl;
    pflow::logger(soliloquy,aurostd::PaddedCENTER(multiple_configs_ss.str(),AFLOWIN_SEPARATION_LINE.size()+2),*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
  }
  return true;
}

//tol -> HNF optimizer
bool POccStructure::calculateNHNF(){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccStructure::calculateNHNF():";

  stringstream message,message_raw;
  message_raw.precision(hnf_table_general_precision);
  message << "Fetching HNF value";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  
  //setHNFTablePadding(PRECISION);
  
  //vector<vector<POccSiteConfiguration> > v_str_configs;

  double stoich_error,site_error;
  stoich_error=site_error=1.0;
  if(LDEBUG){
    cerr << soliloquy << " POcc structure:" << endl;
    cerr << xstr_pocc;
  }
  if (xstr_pocc.partial_occupation_HNF!=0) {
    n_hnf=xstr_pocc.partial_occupation_HNF;
    message << "Using input HNF value = " << n_hnf;
    pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    //getSiteCountConfigurations(n_hnf,stoich_error);
    getSiteCountConfigurations(n_hnf);

    message_raw << AFLOWIN_SEPARATION_LINE << endl; pflow::logger(soliloquy,message_raw,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    pflow::logger(soliloquy,hnfTableHeader(),*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    if(!getHNFTableOutput(n_hnf,stoich_error,site_error)){return false;}
    message_raw << AFLOWIN_SEPARATION_LINE << endl; pflow::logger(soliloquy,message_raw,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);

    return true;
  } //if HNF exists, no need to optimize pocc values


  //retrieve input/default tolerance
  const double& stoich_tolerance=xstr_pocc.partial_occupation_stoich_tol;
  const double& site_tolerance=xstr_pocc.partial_occupation_site_tol;
  //stoich_tolerance=site_tolerance=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
  //if(xstr_pocc.partial_occupation_site_tol>0) {stoich_tolerance=xstr_pocc.partial_occupation_stoich_tol;}
  //if(xstr_pocc.scale_third.isentry) {site_tolerance=xstr_pocc.partial_occupation_site_tol;}
  //cerr << "stoich_tolerance=" << stoich_tolerance << " site_tolerance=" << site_tolerance << endl;
  //exit(0);

  //output - START
  message << "Stoichiometry tolerance = " << stoich_tolerance;
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  message << "Site tolerance = " << site_tolerance;
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  message << "Optimizing HNF value";
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

  //output - END

  int i_hnf=0;

  // COREY seems everything is working up to here
  //error needs to be defined not PER site but for the whole structure
  //make header stoich, not sites
  //split header and output, show header and output for fixed hnf
  //

  //pflow::logger(_LOGGER_RAW_,soliloquy,hnfTableHeader(),*p_oss,*p_FileMESSAGE);
  //HEADER - STOP
  
  //stringstream aus;aus.precision(hnf_table_general_precision);//aus.setf(std::ios::floatfield);
  //cerr << "HERE1" << endl;
  //cerr << "HERE2" << endl;
  //site_error=1.0;
  
  message_raw << AFLOWIN_SEPARATION_LINE  << endl; pflow::logger(soliloquy,message_raw,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
  pflow::logger(soliloquy,hnfTableHeader(),*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
  while(stoich_error>stoich_tolerance || site_error>site_tolerance) {
    i_hnf++;
    //site_error=1.0/i_hnf;
    //cerr << "HERE3" << endl;
    //main calculation
    //getSiteCountConfigurations(i_hnf,stoich_error);
    getSiteCountConfigurations(i_hnf);
    //cerr << "HERE4" << endl;
    if(!getHNFTableOutput(i_hnf,stoich_error,site_error)){return false;}
  }
  message_raw << AFLOWIN_SEPARATION_LINE << endl; 
  pflow::logger(soliloquy,message_raw,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
  //TABLE - STOP

  //DO NOT UPDATE YET, we can have multiple configurations, so we do one a time later
  //set values in xstr_pocc
  /*for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {
    if(!aurostd::isequal(xstr_pocc.atoms[iatom].partial_occupation_value,1.0,ZERO_TOL)) {
      xstr_pocc.atoms[iatom].partial_occupation_value = v_site_effective_pocc[iatom];
    }
  }*/

  //if(!updatePOccValues()){return false;} //update pocc_values and comp_each_type
  n_hnf=i_hnf;

  message << "Optimized HNF value = " << n_hnf;
  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  // COREY FIX ME
  //for(uint site=0;site<vv_count_configs.size();site++){
  //  message << "Site " << aurostd::PaddedPRE(aurostd::utype2string(site+1),hnf_table_iteration_padding);
  //  message << "/" << aurostd::PaddedPOST(aurostd::utype2string(site+1),hnf_table_iteration_padding) << ": number of count configurations = " << vv_count_configs[site].size();
  //  pflow::logger(soliloquy,message,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
  //}
  //exit(0);
  return true;
}

} //namespace pocc

namespace pocc {
//--------------------------------------------------------------------------------
// class POccSuperStructure (nested in POccCalculator)
//--------------------------------------------------------------------------------
POccSuperStructure::POccSuperStructure() {
  //------------------------------------------------------------------------------
  // constructor
  free();
  //hnf_mat=new xmatrix<double>();              //safe
  //p_str=new POccStructure();                  //safe
  //v_types_config=new vector<vector<int> >(0); //safe
}

POccSuperStructure::~POccSuperStructure() {
  free();
}

void POccSuperStructure::free() {
  POccStructureTemplate::free();
  hnf_mat.clear();
  xstr_ss.Clear();
  xstr_cluster.Clear();
  //cluster.clear();
  for(uint i=0;i<v_bonded_atom_indices.size();i++){v_bonded_atom_indices[i].clear();} v_bonded_atom_indices.clear();
  for(uint i=0;i<v_nonbonded_atom_indices.size();i++){v_nonbonded_atom_indices[i].clear();} v_nonbonded_atom_indices.clear();
}

void POccSuperStructure::copy(const POccSuperStructure& b) { // copy PRIVATE
  hnf_mat=b.hnf_mat;
  xstr_pocc=b.xstr_pocc;
  xstr_nopocc=b.xstr_nopocc;
  types2pc_map=b.types2pc_map;
  types2uffparams_map=b.types2uffparams_map;
  //p_str=b.p_str;
  //v_types_config=b.v_types_config;
  //xstr_nopocc=b.xstr_nopocc;
  xstr_ss=b.xstr_ss;
  xstr_cluster=b.xstr_cluster;
  //cluster=b.cluster;
  v_bonded_atom_indices=b.v_bonded_atom_indices;
  v_nonbonded_atom_indices=b.v_nonbonded_atom_indices;
}

const POccSuperStructure& POccSuperStructure::operator=(const POccSuperStructure& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

POccSuperStructure::POccSuperStructure(const POccSuperStructure& b) { // copy PUBLIC
  copy(b);
}

void POccSuperStructure::clear() {  // clear PRIVATE
  POccSuperStructure _temp;
  copy(_temp);
}

void POccSuperStructure::calculateNNDistances(xstructure& xstr,vector<uint>& v_vacancies){
  xmatrix<double> _distance_matrix(xstr.atoms.size()-1,xstr.atoms.size()-1,0,0); distance_matrix=_distance_matrix;

  //get distance matrix first, then find nearest-neighbor distances
  xvector<double> min_vec; xvector<int> ijk;  //dummy
  //xmatrix<double> lattice(p_str->getLattice()); //hack, need copy to avoid const issues
  
  //get distance matrix first
  for(uint atom1=0;atom1<xstr.atoms.size();atom1++){
    if(isVacancy(v_vacancies,atom1)){continue;}
    for(uint atom2=atom1+1;atom2<xstr.atoms.size();atom2++){
      if(isVacancy(v_vacancies,atom2)){continue;}
      //distance_matrix(atom1,atom2)=distance_matrix(atom2,atom1)=SYM::minimumCartesianDistance(xstr.atoms[atom1].cpos,xstr.atoms[atom2].cpos,xstr.lattice,min_vec,ijk);
      distance_matrix(atom1,atom2)=distance_matrix(atom2,atom1)=SYM::minimumCartesianDistance(xstr.atoms[atom1].cpos,xstr.atoms[atom2].cpos,xstr.lattice,min_vec,ijk);
    }
  }

  //now get n-n distances
  v_dist_nn.clear();
  for(uint i=0;i<xstr.atoms.size();i++){v_dist_nn.push_back(std::numeric_limits<double>::max());} //initialize with big number
  for(uint atom1=0;atom1<xstr.atoms.size();atom1++){
    if(isVacancy(v_vacancies,atom1)){continue;}
    for(uint atom2=0;atom2<xstr.atoms.size();atom2++){
      if(isVacancy(v_vacancies,atom2)){continue;}
      if(atom1!=atom2 && distance_matrix(atom1,atom2)<v_dist_nn[atom1]){v_dist_nn[atom1]=distance_matrix(atom1,atom2);}
    }
  }
}

void POccSuperStructure::initialize(POccStructure& p_str,double _exploration_radius){
  setPOccStructure(p_str.xstr_pocc);
  setNonPOccStructure(p_str.xstr_nopocc);
  setTypes2PcMap(p_str.types2pc_map);
  setTypes2UFFParamsMap(p_str.types2uffparams_map);

  exploration_radius=_exploration_radius;
  //xstr_nopocc=&_xstr_nopocc;
  //const xstructure& xstr_nopocc = p_str->getNonPOccStructure();
}
  
void POccSuperStructure::getCluster(xmatrix<double>& _hnf_mat){
  hnf_mat=_hnf_mat;
  sc2pc_map.clear(); pc2sc_map.clear();
  xstr_ss=GetSuperCell(xstr_nopocc,hnf_mat,sc2pc_map,pc2sc_map,false,false);

  //double radius=10;
  //xvector<int> dims(3);
  //double radius = ENERGY_RADIUS/2.0;  //LatticeDimensionSphere radius fits within 2*dim
  //double radius = exploration_radius/2.0;  //LatticeDimensionSphere radius fits within 2*dim
  //cerr << "exploration radius = " << exploration_radius << endl;
  //cerr << RadiusSphereLattice(xstr_ss.lattice) << endl;
  //safety, no need, gives 1x1x1 as long as radius >= 0.0
  //if(RadiusSphereLattice(xstr_ss.lattice)<radius){dims(1)=1;dims(2)=1;dims(3)=3;}
  //else{dims=LatticeDimensionSphere(xstr_ss.lattice,radius);}
  //cerr << xstr_ss.lattice << endl;
  //dims=LatticeDimensionSphere(xstr_ss.lattice,exploration_radius);
  //cerr << "dims = " << dims << endl;

  //xmatrix<double> slattice(3,3);
  //slattice(1,1)=dims(1);slattice(2,2)=dims(2);slattice(3,3)=dims(3);
  //cerr << "START CLUSTER" << endl;
  //cluster=GetCluster(xstr_ss,radius,ssc2sc_map,sc2ssc_map); //GetSuperCell(xstr_ss,slattice,ssc2sc_map,sc2ssc_map,false,false);
  //xstr_cluster=GetSuperCell(xstr_ss,slattice,ssc2sc_map,sc2ssc_map,false,false);
  
  //[OBSOLETE]ssc2sc_map.clear(); sc2ssc_map.clear();
  //[OBSOLETE]cluster=GetCluster(xstr_ss,2.0*exploration_radius,ssc2sc_map,sc2ssc_map);
  xstr_cluster=xstr_ss;
  GenerateGridAtoms(xstr_cluster,exploration_radius);
  
  bonding_set=false;
  //cerr << "ssc2sc_map" << endl;
  //int ko=-1;
  //for(uint i=0;i<xstr_cluster.grid_atoms_sc2pcMap.size();i++){
  //  if(ko!=(int)xstr_cluster.grid_atoms_sc2pcMap[i]){
  //  cerr << "i=" << i << " " << xstr_cluster.grid_atoms_sc2pcMap[i] << endl;
  //  ko=(int)xstr_cluster.grid_atoms_sc2pcMap[i];
  //  }
  //}
  //cerr << xstr_cluster.grid_atoms_number << endl;
  //cerr << "STOP CLUSTER" << endl;
  //cerr << "CLUSTER " << endl;
  
  //xstr_bonding=xstr_nopocc; //will be updated if vacancies present
}

//bool POccSuperStructure::hasVacancies(vector<vector<int> >& v_types_config){
//bool POccSuperStructure::hasVacancies(){
//  return v_vacancies.size()>0;
  //for(uint i=0;i<v_types_config.size();i++){
  //  for(uint j=0;j<v_types_config[i].size();j++){
  //    if(v_types_config[i][j]==-1){return true;}
  //  }
  //}
  //return false;
//}

//nice function pointers so we don't need to rewrite so much code, and for loops are fast
uint POccSuperStructure::NNDistancesMapPC(uint atom){return sc2pc_map[xstr_cluster.grid_atoms_sc2pcMap[atom]];}
uint POccSuperStructure::NNDistancesMapSC(uint atom){return xstr_cluster.grid_atoms_sc2pcMap[atom];}

//shortcut for recalculating bonds - only if bonds have been calculated already
//take in current v_types_config
//we will generate a list of vacancies
//then we compare with previously determined vacancies
//if they match, do not recalculate
void POccSuperStructure::setBonds(vector<vector<int> >& v_types_config){

  //cout << "v_types_config_bonding OLD ";
  //for(uint i=0;i<v_types_config_bonding.size();i++){for(uint j=0;j<v_types_config_bonding[i].size();j++){cout << v_types_config_bonding[i][j] << " ";} }
  //cout << endl;
  //cout << "v_types_config         NEW ";
  //for(uint i=0;i<v_types_config.size();i++){for(uint j=0;j<v_types_config[i].size();j++){cout << v_types_config[i][j] << " ";} }
  //cout << endl;


  v_types_config_bonding=v_types_config;
  //vector<uint> v_vacancies=getVacancies(xstr_ss, v_types_config, v_vacancies, false);  //do not modify xstr_ss (no replace_types, just get vacancies)
  vector<uint> v_vacancies=getVacancies(v_types_config);  //get vacancies //do not modify xstr_ss (no replace_types, just get vacancies)
  has_vacancies=(v_vacancies.size()>0);



  //cout << "v_vacancies_bonding OLD ";
  //for(uint i=0;i<v_vacancies_bonding.size();i++){cout << v_vacancies_bonding[i] << " ";}
  //cout << endl;
  //cout << "v_vacancies         NEW ";
  //for(uint i=0;i<v_vacancies.size();i++){cout << v_vacancies[i] << " ";}
  //cout << endl;
 


  xstructure* ptr_xstr_bonding=&xstr_nopocc;  //slight optimization, if no vacancies, use primitive cell n-n distances
  uint (POccSuperStructure::*NNDistancesMap)(uint)=NULL; //nice function pointers so we don't need to rewrite so much code, and for loops are fast
  NNDistancesMap=&POccSuperStructure::NNDistancesMapPC;
 
  //cerr << "CHERE4" << endl;
  if(has_vacancies){
    if(bonding_set){bonding_set=(v_vacancies==v_vacancies_bonding);}  //bonding_set==false when cluster is generated, check to recalculate only if has been calculated before
    ptr_xstr_bonding=&xstr_ss;
    NNDistancesMap=&POccSuperStructure::NNDistancesMapSC;
    //xstr_bonding=xstr_ss; //copies it all over
    //removeVacancies(xstr_bonding,v_vacancies);
  //}else{
  }

  //cerr << "CHERE5" << endl;
  xstructure& xstr_bonding = *ptr_xstr_bonding;
  //cout << distance_matrix << endl;

  if(!bonding_set){
    //set v_types_config_bonding and v_vacancies_bonding!!!!!
    //cout << "HERE" << endl;
    //cerr << "CHERE6" << endl;
    v_vacancies_bonding=v_vacancies;
    bonding_set=true;
    calculateNNDistances(xstr_bonding,v_vacancies_bonding); //v_vacancies_bonding is empty if xstr_bonding==xstr_nopocc, so this is okay
    //cerr << "CHERE7" << endl;

    //cout << distance_matrix << endl;

    //now get bonding indices for super-superstructure
    double distij;
    for(uint i=0;i<v_bonded_atom_indices.size();i++){v_bonded_atom_indices[i].clear();} v_bonded_atom_indices.clear();
    for(uint i=0;i<v_nonbonded_atom_indices.size();i++){v_nonbonded_atom_indices[i].clear();} v_nonbonded_atom_indices.clear();
    //cerr << "CHERE8" << endl;
    //if(0){
    //for(uint atom1=0;atom1<xstr_cluster.atoms.size();atom1++){
    //  for(uint atom2=atom1+1;atom2<xstr_cluster.atoms.size();atom2++){
    //    distij=AtomDist(xstr_cluster.atoms[atom1],xstr_cluster.atoms[atom2]);
    //    //if(distij<=radius){
    //    if(abs(distij-v_dist_nn[sc2pc_map[xstr_cluster.grid_atoms_sc2pcMap[atom1]]])<BONDING_THRESHOLD){
    //      v_bonded_atom_indices.push_back(vector<uint>(0));
    //      v_bonded_atom_indices.back().push_back(atom1);
    //      v_bonded_atom_indices.back().push_back(atom2);
    //      //cerr << "BONDING    " << atom1 << " " << atom2 << endl;
    //    }else{
    //      v_nonbonded_atom_indices.push_back(vector<uint>(0));
    //      v_nonbonded_atom_indices.back().push_back(atom1);
    //      v_nonbonded_atom_indices.back().push_back(atom2);
    //      //cerr << "NONBONDING " << atom1 << " " << atom2 << endl;
    //    }
    //    //}
    //  }
    //}
    //}
    //cerr << "BONDS" << endl;
    uint atom1;
    //uint count;
    for(uint sc_atom1=0;sc_atom1<xstr_ss.atoms.size();sc_atom1++){
      //count=0;
      if(isVacancy(v_vacancies_bonding,sc_atom1)){continue;}
      atom1=xstr_cluster.grid_atoms_pc2scMap[sc_atom1];
      //for(uint atom2=0;atom2<xstr_cluster.atoms.size();atom2++){
      for(uint atom2=0;atom2<(uint)xstr_cluster.grid_atoms_number;atom2++){
        if(atom1!=atom2){
        //distij=AtomDist(xstr_cluster.atoms[atom1],xstr_cluster.atoms[atom2]);
        distij=AtomDist(xstr_cluster.grid_atoms[atom1],xstr_cluster.grid_atoms[atom2]);
        if(distij<=exploration_radius){
        //  cerr << "OK cluster atom 1 " << atom1 << " " << xstr_cluster.grid_atoms[atom1].cpos << endl;
        //  cerr << "OK cluster atom 2 " << atom2 << " " << xstr_cluster.grid_atoms[atom2].cpos << endl;
        //if(abs(distij-v_dist_nn[sc2pc_map[xstr_cluster.grid_atoms_sc2pcMap[atom1]]])<=BONDING_THRESHOLD){
        //if(abs(distij-v_dist_nn[ ( has_vacancies ? xstr_cluster.grid_atoms_sc2pcMap[atom1] : sc2pc_map[xstr_cluster.grid_atoms_sc2pcMap[atom1]]  )        ])<=BONDING_THRESHOLD){
        if(abs(distij-v_dist_nn[(this->*NNDistancesMap)(atom1)])<=BONDING_THRESHOLD){
          v_bonded_atom_indices.push_back(vector<uint>(0));
          v_bonded_atom_indices.back().push_back(atom1);
          v_bonded_atom_indices.back().push_back(atom2);
          //cerr << "BOND " << atom1 << " " << atom2 << endl;
          //count++;
          //cerr << "BONDING    " << atom1 << " " << atom2 << endl;
        }else{
          v_nonbonded_atom_indices.push_back(vector<uint>(0));
          v_nonbonded_atom_indices.back().push_back(atom1);
          v_nonbonded_atom_indices.back().push_back(atom2);
          //cerr << "NONBONDING " << atom1 << " " << atom2 << endl;
        }
        }
        }
      }
      //cerr << "_atom1=" << _atom1 << " " << count << endl;
    }
    //exit(0);

    //for(uint i=0;i<v_dist_nn.size();i++){cerr << "i=" << i << " " << v_dist_nn[i] << endl;}
    //cerr << distance_matrix << endl;
    //exit(0);

    //cerr << *xstr_nopocc << endl;
    //cerr << xstr_ss << endl;
    //cerr << dims << endl;
    //cerr << xstr_cluster.atoms << endl;
    //exit(0);
    //cerr << "DONE" << endl;
  }
}

//xstructure POccSuperStructure::replaceTypes(vector<vector<int> >& _v_types_config, vector<uint>& v_vacancies) {
//void POccSuperStructure::replaceTypes(xstructure& supercell, vector<vector<int> >& v_types_config, vector<uint>& v_vacancies) {
//void POccSuperStructure::getVacancies(xstructure& supercell, vector<uint>& v_vacancies,bool replace_types) {

void POccSuperStructure::replaceRandomSites(xstructure& supercell,const vector<vector<int> >& v_types_config){
  uint starting_supercell_atom_index,pocc_atom_index,supercell_atom_index;
  int type;
  for(uint site=0;site<v_types_config.size();site++){
    starting_supercell_atom_index=pc2sc_map[site];
    //cerr << "starting_supercell_atom_index=" << starting_supercell_atom_index << endl;
    for(uint i=0;i<v_types_config[site].size();i++){
      type=v_types_config[site][i];
      //switch these atoms
      //supercell.atoms[starting_supercell_atom_index+i]
      //xstr_pocc.atoms[types2pc_map[(*v_types_config)[site][i]]]
      //pocc_atom_index=types2pc_map[(*v_types_config)[site][i]];
      supercell_atom_index=starting_supercell_atom_index+i;
      if(type>=0){
        //if(replace_types){
        pocc_atom_index=types2pc_map[type]; //p_str->types2pcMap(type);
        const _atom& pocc_atom = xstr_pocc.atoms[pocc_atom_index]; //xstr_pocc.atoms[pocc_atom_index]; //p_str->xstr_pocc.atoms[pocc_atom_index];
        //cerr << pocc_atom << endl;
        //cerr << supercell_atom_index << endl;
        //cerr << supercell.atoms.size() << endl;
        copyAtomAttributes(pocc_atom,supercell.atoms[supercell_atom_index]);
        supercell.atoms[supercell_atom_index].type=type;
        supercell.atoms[supercell_atom_index].partial_occupation_flag=false;
        supercell.atoms[supercell_atom_index].partial_occupation_value=1.0;
        //supercell.num_each_type[type]++;
        //supercell.comp_each_type[type]+=supercell.atoms[supercell_atom_index].partial_occupation_value;
        //cerr << "supercell.atoms[" << supercell_atom_index << "].type=" << supercell.atoms[supercell_atom_index].type << endl;
        //}
      }
    }
  }
}

vector<uint> POccSuperStructure::getVacancies(const vector<vector<int> >& v_types_config) {
  //v_vacancies.clear();
  vector<uint> v_vacancies;
  
  uint starting_supercell_atom_index,/*pocc_atom_index,*/supercell_atom_index;
  int type;
  for(uint site=0;site<v_types_config.size();site++){
    //if(v_sites[site].partial_occupation_flag){  //slight optimization
    //  //test of stupidity
    //  if(v_sites[site].v_occupants.size()<2){cerr << "WHOOOOOOPS" << endl;exit(1);}
    //find index of atom first in v_occupants list, that's the one that remained
    //cerr << "v_sites[" << site << "].v_occupants[0]=" << v_sites[site].v_occupants[0] << endl;
    //cerr << "pc2sc_map.size()=" << pc2sc_map.size() << endl;
    //cerr << "pc2sc_map[0]=" << pc2sc_map[0] << endl; 
    //cerr << "pc2sc_map[1]=" << pc2sc_map[1] << endl; 
    starting_supercell_atom_index=pc2sc_map[site];//v_sites[site].v_occupants[0]];
    //cerr << "starting_supercell_atom_index=" << starting_supercell_atom_index << endl;
    //for(uint i=0;i<(*v_types_config)[site].size();i++){
    for(uint i=0;i<v_types_config[site].size();i++){
      type=v_types_config[site][i];
      //switch these atoms
      //supercell.atoms[starting_supercell_atom_index+i]
      //xstr_pocc.atoms[types2pc_map[(*v_types_config)[site][i]]]
      //pocc_atom_index=types2pc_map[(*v_types_config)[site][i]];
      supercell_atom_index=starting_supercell_atom_index+i;
      if(type<0){v_vacancies.push_back(supercell_atom_index);}
    }
  }

  std::sort(v_vacancies.rbegin(),v_vacancies.rend());
  return v_vacancies;
}

bool POccSuperStructure::isVacancy(vector<uint>& v_vacancies, uint atom){
  for(uint vi=0;vi<v_vacancies.size();vi++){if(atom==v_vacancies[vi]){return true;}}
  return false;
}

double POccSuperStructure::getBondEnergy(xstructure& xstr,vector<vector<uint> >& v_bonded_atom_indices,uint MODE) {

  string soliloquy="POccSuperStructure::getBondEnergy()";
  stringstream message;
  //cerr << xstr << endl;
  double energy=0.0;
  //bool found_vacancy;
  double distij;
  uint atom1,atom2,type1,type2;
  UFFParamAtom uffai, uffaj;
  UFFParamBond uffb;
  uint count=0;
  for(uint i=0;i<v_bonded_atom_indices.size();i++){
    //cerr << "HERE1" << endl;
    
    atom1=xstr_cluster.grid_atoms_sc2pcMap[v_bonded_atom_indices[i][0]];
    if(isVacancy(v_vacancies_bonding,atom1)){continue;}
    
    atom2=xstr_cluster.grid_atoms_sc2pcMap[v_bonded_atom_indices[i][1]];
    if(isVacancy(v_vacancies_bonding,atom2)){continue;}
    
    //cerr << "HERE2" << endl;
    //found_vacancy=false;
    //skip anything with a vacancy - START
    //for(uint vi=0;vi<v_vacancies.size() && !found_vacancy;vi++){if(atom1==v_vacancies[vi] || atom2==v_vacancies[vi]){found_vacancy=true;}}
    //if(found_vacancy){
      //cerr << endl;
      //cerr << "SKIPPING!" << endl;
      //cerr << "atom1 " << v_bonded_atom_indices[i][0] << "(" << atom1 << ")  " << "(cluster " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]].cpos << ") (supercell " << xstr.atoms[atom1].cpos << ")" << endl;
      //cerr << "atom3 " << v_bonded_atom_indices[i][1] << "(" << atom2 << ")  " << "(cluster " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]].cpos << ") (supercell " << xstr.atoms[atom2].cpos << ")" << endl;
      //cerr << "SKIPPING!" << endl;
      //cerr << endl;
    //  continue;
    //}
    count++;
    //skip anything with a vacancy - END
    //cerr << "atom1=" << atom1 << endl;
    //cerr << "atom2=" << atom2 << endl;
    //cerr << xstr.atoms.size() << endl;
    type1=xstr.atoms[atom1].type;
    type2=xstr.atoms[atom2].type;
    //cerr << "v_bonded_atom_indices[i][0]=" << v_bonded_atom_indices[i][0] << endl;
    //cerr << "v_bonded_atom_indices[i][1]=" << v_bonded_atom_indices[i][1] << endl;
    //cerr << xstr_cluster.grid_atoms_number << endl;
    //distij=AtomDist(xstr_cluster.atoms[v_bonded_atom_indices[i][0]],xstr_cluster.atoms[v_bonded_atom_indices[i][1]]);
    distij=AtomDist(xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]],xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]]);
    if(distij>exploration_radius){
      message << "Attempting to explore atoms outside of the search radius ";
      message << "(";
      message << "cluster atom 1 " << v_bonded_atom_indices[i][0] << " " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]].cpos;
      message << " vs. ";
      message << "cluster atom 2 " << v_bonded_atom_indices[i][1] << " " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]].cpos;
      message << ")";
      throw AFLOWLogicError(soliloquy,message);
    }
    //cerr << "OK " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]].cpos << " " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]]  << " " ;
    //cerr << "DIST " << distij <<endl;
    //cerr << "HERE4" << endl;
    //cerr << "sss atom1=" << v_bonded_atom_indices[i][0] << " type=" << xstr_cluster.atoms[v_bonded_atom_indices[i][0]].type << " ";
    //cerr << "sss atom2=" << v_bonded_atom_indices[i][1] << " type=" << xstr_cluster.atoms[v_bonded_atom_indices[i][1]].type << endl;
    //cerr << "ss  atom1=" << atom1 << " type=" << type1 << " ";
    //cerr << "ss  atom2=" << atom2 << " type=" << type2 << endl;
    //cerr << "distij=" << distij << " " << (MODE==BOND_MODE ? "BOND" : "NONBOND") << endl;
    uffai=types2uffparams_map[type1];
    uffaj=types2uffparams_map[type2];
    uffb.calculate(uffai,uffaj,distij);
    if(MODE==BOND_MODE){energy += 0.5 * uffb.Kij * uffb.delta * uffb.delta;
      //cerr << "atomi " << xstr.atoms[atom1] << endl;
      //cerr << "atomj " << xstr.atoms[atom2] << endl;
      //cerr << "distij " << distij << endl;
      //cerr << "uffb.Kij" << uffb.Kij << endl;
      //cerr << "uffb.delta" << uffb.delta << endl;
      //cerr << "uffb.ren " << uffb.ren << endl;
      //cerr << "uffb.R0 " << uffb.R0 << endl;
      //cerr << "energy " << KCAL_2_EV * 0.5 * uffb.Kij * uffb.delta * uffb.delta << endl;
    }        //spring potential
    else if(MODE==NONBOND_MODE){energy += uffb.Dij * (uffb.X12 - 2.0 * uffb.X6);
      //cerr << "type1 " << type1 << endl;
      //cerr << "type2 " << type2 << endl;
      //cerr << "p_str->types2uffparams_map[type1] " << p_str->types2uffparams_map[type1] << endl;
      //cerr << "p_str->types2uffparams_map[type2] " << p_str->types2uffparams_map[type2] << endl;
      //cerr << "atomi " << xstr.atoms[atom1] << endl;
      //cerr << "atomj " << xstr.atoms[atom2] << endl;
      //cerr << "atomi.type " << xstr.atoms[atom2].type << endl;
      //cerr << "atomj.type " << xstr.atoms[atom2].type << endl;
      //cerr << "distij " << distij << endl;
      //cerr << uffb.Dij * (uffb.X12 - 2.0 * uffb.X6) << endl;
    }   //lennard-jones potential
    else{
      message << "Invalid energy mode (MODE=" << MODE << ")";
      throw AFLOWLogicError(soliloquy,message);
    }
  }
  //if(MODE==BOND_MODE) cerr << "BONDING LENGTH C " << count << endl;
  energy*=KCAL_2_EV;
  return energy;
}

//bool sortAtoms(const _atom& a1,const _atom& a2) {
  //if(a1.type!=a2.type){return a1.type<a2.type;}
  //double dist1=aurostd::modulus(a1.cpos);
  //double dist2=aurostd::modulus(a2.cpos);
  //return dist1<dist2;
//}

//void POccSuperStructure::removeVacancies(xstructure& xstr,vector<uint> v_vacancies) {
//void POccSuperStructure::removeVacancies(xstructure& xstr,vector<uint>& v_vacancies) {
void POccSuperStructure::removeVacancies(deque<_atom>& atoms,vector<uint>& v_vacancies) {
  //only used in rebuildStructure(), so don't worry about doing a full RemoveAtom()
  if(v_vacancies.size()){
    std::sort(v_vacancies.rbegin(),v_vacancies.rend());
    for(uint i=0;i<v_vacancies.size();i++){atoms.erase(atoms.begin()+v_vacancies[i]);}
  }
}

void POccSuperStructure::rebuildStructure(xstructure& xstr,vector<uint>& v_vacancies) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccSuperStructure::rebuildStructure():";
  stringstream message;
  
  //remove vacancies from atoms first
  //removeVacancies(xstr,v_vacancies);
  //we have the right atoms, we just need to rearrange and get right counts
  //std::sort(xstr.atoms.begin(),xstr.atoms.end(),sortAtomsDist); //pocc::sortAtoms);
  //if(LDEBUG){cerr << soliloquy << " new atom count=" << xstr.atoms.size() << endl;}

  deque<_atom>& atoms=xstr.atoms;
  removeVacancies(atoms,v_vacancies);
  std::stable_sort(atoms.begin(),atoms.end(),sortAtomsDist); //pocc::sortAtoms);
  if(LDEBUG){cerr << soliloquy << " new atom count=" << atoms.size() << endl;}

  //I have two options here, I can:
  //1) make a full copy of atoms and use AddAtom(), followed by copying over all pp info, or 
  //2) I can simply copy all the stuff over from xstr_pocc

  //get vector of applicable types
  vector<uint> new_types; //vacancies are -1, and they have been removed
  vector<uint> new_types_atoms;
  vector<uint> new_species;
  bool found;
  for(uint atom=0;atom<atoms.size();atom++){
    found=false;
    for(uint type=0;type<new_types.size()&&!found;type++){
      if(atoms[atom].type==(int)new_types[type]){found=true;}
    }
    if(!found){
      new_types.push_back(atoms[atom].type);
      new_types_atoms.push_back(atom);
    }
  }

  //make sure they are in order
  for(uint type=1;type<new_types.size();type++){
    if(new_types[type-1]>new_types[type]){throw AFLOWLogicError(soliloquy,"Error in atom types reassignment, atoms are out of order");}
  }

  //load in to new clean structure
  xstructure xstr_clean;
  xstr_clean.lattice=xstr.lattice;
  xstr_clean.scale=xstr.scale;
  xstr_clean.FixLattices();

  for(uint atom=0;atom<atoms.size();atom++){xstr_clean.AddAtom(atoms[atom]);}

  //we need to copy the follow species info over, which is not copied with AddAtom(), space is created for them though
  //species_pp_type
  //species_pp_version
  //species_pp_ZVAL
  //species_pp_vLDAU

  uint atom,type;
  for(uint it=0;it<new_types.size();it++){
    type=new_types[it];
    atom=new_types_atoms[it];
    if(LDEBUG){
      cerr << soliloquy << " type[" << it << "]=" << type << endl;
      cerr << soliloquy << " atom[" << it << "]=" << atom << endl;
    }
    //the if statement ensures we don't seg fault
    if((uint)atoms[atom].type<xstr_pocc.species_pp_type.size()){xstr_clean.species_pp_type[it]=xstr_pocc.species_pp_type[type];}  //species_pp_type
    if((uint)atoms[atom].type<xstr_pocc.species_pp_version.size()){xstr_clean.species_pp_version[it]=xstr_pocc.species_pp_version[type];}  //species_pp_version
    if((uint)atoms[atom].type<xstr_pocc.species_pp_ZVAL.size()){xstr_clean.species_pp_ZVAL[it]=xstr_pocc.species_pp_ZVAL[type];}  //species_pp_ZVAL
    if((uint)atoms[atom].type<xstr_pocc.species_pp_vLDAU.size()){xstr_clean.species_pp_vLDAU[it]=xstr_pocc.species_pp_vLDAU[type];}  //species_pp_vLDAU
  }

  if(0){  //this works, but it's not stable, what if we add more species properties to xstructure? better to start fresh and load in what you need!
    xstr.ClearSpecies();
    for(uint type=0;type<new_types.size();type++){
      xstr.num_each_type.push_back(0);
      xstr.comp_each_type.push_back(0.0);
    }

    uint last_type=0;
    for(uint atom=0;atom<xstr.atoms.size();atom++){
      found=false;
      for(uint type=last_type;type<new_types.size()&&!found;type++){
        if(xstr.atoms[atom].type==(int)new_types[type]){
          xstr.num_each_type[type]++;
          xstr.comp_each_type[type]+=xstr.atoms[atom].partial_occupation_value;
          last_type=type;
          found=true;
        }
      }
    }
    
    //clear what you can with ClearSpecies()
    //add here what's needed for species specification in POSCAR
    for(uint type=0;type<new_types.size();type++){
      if(type>xstr_pocc.species.size()-1){
        message << "Bad index with xstr.species[i=" << type << ">" << xstr_pocc.species.size()-1 << "=species.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species.push_back(xstr_pocc.species[type]);
      
      if(type>xstr_pocc.species_pp.size()-1){
        message << "Bad index with xstr.species_pp[i=" << type << ">" << xstr_pocc.species_pp.size()-1 << "=species_pp.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species_pp.push_back(xstr_pocc.species_pp[type]);
      
      if(type>xstr_pocc.species_pp_type.size()-1){
        message << "Bad index with xstr.species_pp_type[i=" << type << ">" << xstr_pocc.species_pp_type.size()-1 << "=species_pp_type.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species_pp_type.push_back(xstr_pocc.species_pp_type[type]);
      
      if(type>xstr_pocc.species_pp_version.size()-1){
        message << "Bad index with xstr.species_pp_version[i=" << type << ">" << xstr_pocc.species_pp_version.size()-1 << "=species_pp_version.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species_pp_version.push_back(xstr_pocc.species_pp_version[type]);
      
      if(type>xstr_pocc.species_pp_ZVAL.size()-1){
        message << "Bad index with xstr.species_pp_ZVAL[i=" << type << ">" << xstr_pocc.species_pp_ZVAL.size()-1 << "=species_pp_ZVAL.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species_pp_ZVAL.push_back(xstr_pocc.species_pp_ZVAL[type]);
      
      if(type>xstr_pocc.species_pp_vLDAU.size()-1){
        message << "Bad index with xstr.species_pp_vLDAU[i=" << type << ">" << xstr_pocc.species_pp_vLDAU.size()-1 << "=species_pp_vLDAU.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species_pp_vLDAU.push_back(xstr_pocc.species_pp_vLDAU[type]);
      
      if(type>xstr_pocc.species_volume.size()-1){
        message << "Bad index with xstr.species_volume[i=" << type << ">" << xstr_pocc.species_volume.size()-1 << "=species_volume.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species_volume.push_back(xstr_pocc.species_volume[type]);
      
      if(type>xstr_pocc.species_mass.size()-1){
        message << "Bad index with xstr.species_mass[i=" << type << ">" << xstr_pocc.species_mass.size()-1 << "=species_mass.size()-1], likely not properly copied over";
        throw AFLOWLogicError(soliloquy,message);
      }
      xstr.species_mass.push_back(xstr_pocc.species_mass[type]);
    }

    //fix volume and mass if they are missing
    if(xstr_pocc.species.size()){
      //volume
      double volume_sum=0.0;
      for(uint species=0;species<xstr_pocc.species_volume.size();species++){volume_sum+=xstr_pocc.species_volume[species];}
      if(volume_sum<_ZERO_TOL_){
        xstr_pocc.species_volume.clear();
        for(uint species=0;species<xstr_pocc.species.size();species++){
          xstr_pocc.species_volume.push_back(GetAtomVolume(xstr_pocc.species[species]));
        }
      }
      //mass
      double mass_sum=0.0;
      for(uint species=0;species<xstr_pocc.species_mass.size();species++){mass_sum+=xstr_pocc.species_mass[species];}
      if(mass_sum<_ZERO_TOL_){
        xstr_pocc.species_mass.clear();
        for(uint species=0;species<xstr_pocc.species.size();species++){
          xstr_pocc.species_mass.push_back(GetAtomMass(xstr_pocc.species[species]));
        }
      }
    }
    //leave for end - fix types
    last_type=0;
    for(uint atom=0;atom<xstr.atoms.size();atom++){
      found=false;
      for(uint type=last_type;type<new_types.size()&&!found;type++){
        if(xstr.atoms[atom].type==(int)new_types[type]){
          xstr.atoms[atom].type=type;
          found=true;
        }
      }
    }
  }

  if(0){  //this assumes that ALL atom types are created... not true, depends on tol/hnf/vacancies
    for(uint i=0;i<types2pc_map.size();i++){
      xstr.num_each_type.push_back(0);
      xstr.comp_each_type.push_back(0.0);
    }

    for(uint atom=0;atom<xstr.atoms.size();atom++){
      xstr.num_each_type[xstr.atoms[atom].type]++;
      xstr.comp_each_type[xstr.atoms[atom].type]+=xstr.atoms[atom].partial_occupation_value;
    }

    //assume we can just copy from xstr_pocc, which should be true unless you screwed something up royally with input
    xstr.species=xstr_pocc.species;
    xstr.species_pp=xstr_pocc.species_pp;
    xstr.species_pp_type=xstr_pocc.species_pp_type;
    xstr.species_pp_version=xstr_pocc.species_pp_version;
    xstr.species_pp_ZVAL=xstr_pocc.species_pp_ZVAL;
    xstr.species_pp_vLDAU=xstr_pocc.species_pp_vLDAU;
    xstr.species_volume=xstr_pocc.species_volume;
    xstr.species_mass=xstr_pocc.species_mass;
    xstr.order_parameter_atoms=xstr_pocc.order_parameter_atoms;
  }
  
  xstr_clean.MakeBasis();
  xstr_clean.MakeTypes();
  
  //change title
  xstr_clean.title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title);
  if(FIX_POSCAR_TITLES || xstr_clean.title.empty()){xstr_clean.buildGenericTitle(true,false);}

  //full copy over
  xstr=xstr_clean;

  //cerr << "FINAL" << endl;
  //cerr << xstr << endl;
}

xstructure POccSuperStructure::createDerivativeStructure(POccDerivativeStructure& pds,int n_hnf,
    unsigned long long int hnf_count,
    unsigned long long int types_config_permutations_count,
    bool clean_structure,bool primitivize){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="POccSuperStructure::createDerivativeStructure():";
  sc2pc_map.clear(); pc2sc_map.clear();
  //vector<int> _sc2pc_map, _pc2sc_map; //make local variants, can't since getVacancies() and rebuild() depend on it, fix later
  
  if(LDEBUG){cerr << soliloquy << " creating superstructure" << endl;}
  xstructure supercell=GetSuperCell(xstr_nopocc,pds.hnf_mat,sc2pc_map,pc2sc_map,false,false);

  if(LDEBUG){cerr << soliloquy << " configuring vacancies (if any)" << endl;}
  //vector<uint> v_vacancies=getVacancies(supercell, pds.v_types_config, true);  //replace types
  vector<uint> v_vacancies=getVacancies(pds.v_types_config);  //get vacancies
  replaceRandomSites(supercell,pds.v_types_config); //replace random sites (ignore vacancies for now, deal in rebuildStructure ONLY)

  if(LDEBUG){cerr << soliloquy << " rebuilding structure" << endl;}
  rebuildStructure(supercell,v_vacancies);    //remove vacancies and clean all xstructure properties

  if(LDEBUG){cerr << soliloquy << " structure built!" << endl;}
  
  if(clean_structure){
    //supercell.MinkowskiBasisReduction();
    //cerr << "PRE PRIM" << endl;
    //cerr << supercell << endl;
    //cerr << "HERE5 " << endl;
    //cerr << "SUPERCELL 1 " << endl;
    //cerr << supercell << endl;
    //xstructure _supercell=GetStandardPrimitive(supercell); 
    if(primitivize){supercell.GetStandardPrimitive();}
    //supercell=GetPrimitiveVASP(supercell);
    //xstructure _supercell=GetPrimitiveVASP(supercell); //GetStandardPrimitive(supercell); 
    //cerr << "SUPERCELL 2" << endl;
    //cerr << _supercell << endl;
    //cerr << "HERE6 " << endl;
    //supercell=_supercell;
    //cerr << "SUPERCELL 3 " << endl;
    //cerr << supercell << endl;
    //cerr << "HERE7 " << endl;
    //cerr << "POST PRIM" << endl;
    //cerr << supercell << endl;
    supercell.ReScale(1.0);
    supercell.ShifOriginToAtom(0);
    supercell.BringInCell();
    //cerr << "HERE8 " << endl;
  }
  //cerr << "HERE9" << endl;
  //[HNF(4)=7/7= 1 0 0; 1 2 0; 1 0 2]
  if(LDEBUG){cerr << soliloquy << " rewriting title" << endl;}
  supercell.title+=" HNF(n="+aurostd::utype2string(n_hnf)+",#"+aurostd::utype2string(pds.hnf_index+1)+"/"+aurostd::utype2string(hnf_count)+")";
  supercell.title+="=["+aurostd::utype2string(pds.hnf_mat(1,1))+" "+aurostd::utype2string(pds.hnf_mat(1,2))+" "+aurostd::utype2string(pds.hnf_mat(1,3))+";";
  supercell.title+= " "+aurostd::utype2string(pds.hnf_mat(2,1))+" "+aurostd::utype2string(pds.hnf_mat(2,2))+" "+aurostd::utype2string(pds.hnf_mat(2,3))+";";
  supercell.title+= " "+aurostd::utype2string(pds.hnf_mat(3,1))+" "+aurostd::utype2string(pds.hnf_mat(3,2))+" "+aurostd::utype2string(pds.hnf_mat(3,3))+"]";
  supercell.title+=" site_config(#"+aurostd::utype2string(pds.site_config_index+1)+"/"+aurostd::utype2string(types_config_permutations_count)+")";
  supercell.title+="=[";
  for(uint site=0;site<pds.v_types_config.size();site++){
    //if(!pds.v_types_config[site].size()){continue;} //should never happen, but be safe //keep off to spot bugs
    for(uint atom=0;atom<pds.v_types_config[site].size();atom++){
      supercell.title+=(site==0&&atom==0?"":" ")+aurostd::utype2string(pds.v_types_config[site][atom]);
    }
    supercell.title+=((site!=pds.v_types_config.size()-1)?";":"");
  }
  supercell.title+="]";
  supercell.title+=" DG="+aurostd::utype2string(pds.degeneracy);
  if(LDEBUG){
    cerr << soliloquy << " new structure" << endl;
    cerr << supercell;
  }
  return supercell;
}

double POccSuperStructure::getEnergy() {

  xstructure supercell(xstr_ss);
  //vector<uint> v_vacancies;
  //replaceTypes(supercell, v_types_config_bonding, v_vacancies);
  //getVacancies(supercell, v_types_config_bonding, true);  //replace types
  replaceRandomSites(supercell,v_types_config_bonding); //replace random sites (ignore vacancies for now, deal in rebuildStructure ONLY)
  //for(uint i=0;i<v_vacancies_bonding.size();i++){
  //  cerr << "vacancy = " << v_vacancies_bonding[i] << endl;
  //}
  //cerr << supercell << endl;

  energy=0.0;
  energy+=getBondEnergy(supercell,v_bonded_atom_indices,BOND_MODE);
  //cerr << "BONDING ENERGY C " << energy << endl;
  energy+=getBondEnergy(supercell,v_nonbonded_atom_indices,NONBOND_MODE);


  //cerr << "WOO" << endl;
  //if(v_vacancies.size()){supercell.RemoveAtoms(v_vacancies);}
  //exit(0);

  //const vector<Site>& v_sites = p_str->v_sites;
  
  //rebuildStructure(supercell,v_vacancies_bonding);
  
  
  if(COMPARE_WITH_KESONG){
    cout << "COREY  " << std::fixed << std::setprecision(15) <<  energy << endl;
    if(SET_KESONG_STANDARD_DIST){
      rebuildStructure(supercell,v_vacancies_bonding); 
      double energy_k=pocc::CalculateUFFEnergy(supercell); 
      cout << "KESONG " << std::fixed << std::setprecision(15) << energy_k << endl; 
      cout << "DIFF " << std::fixed << std::setprecision(15) << abs(energy-energy_k) << endl;
      cout << supercell << endl;
      //exit(0);
    }
  }

  //exit(0);
  //cerr << "AAA" << endl;
  return energy;
}

} // namespace pocc

/*
namespace pocc {
//--------------------------------------------------------------------------------
// class POccDerivativeStructure (nested in POccCalculator)
//--------------------------------------------------------------------------------
POccDerivativeStructure::POccDerivativeStructure() {
  //------------------------------------------------------------------------------
  // constructor
  free();
}        

POccDerivativeStructure::~POccDerivativeStructure() {
  free();
}

void POccDerivativeStructure::free() {
  //------------------------------------------------------------------------------
  for(uint i=0;i<v_occupants.size();i++){v_occupants[i].clear();} v_occupants.clear();
}

void POccDerivativeStructure::copy(const POccDerivativeStructure& b) { // copy PRIVATE
  hnf_index=b.hnf_index;
  v_occupants=b.v_occupants;
  energy=b.energy;
  degerancy=b.degerancy;
}

const POccDerivativeStructure& POccDerivativeStructure::operator=(const POccDerivativeStructure& b) {  // operator= PUBLIC
  if(this!=&b) {free();copy(b);}
  return *this;
}

POccDerivativeStructure::POccDerivativeStructure(const POccDerivativeStructure& b) { // copy PUBLIC
  copy(b);
}

void POccDerivativeStructure::clear() {  // clear PRIVATE
  POccDerivativeStructure _temp;
  copy(_temp);
}

} // namespace pocc
*/

#endif  // _AFLOW_POCC_CPP_

// ***************************************************************************
// *                                                                         *
// *              AFlow STEFANO CURTAROLO  Duke University 2003-2017         *
// *              AFlow COREY OSES  Duke University 2013-2017                *
// *                                                                         *
// ***************************************************************************
