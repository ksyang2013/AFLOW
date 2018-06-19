// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// MAKEFILE FOR AFLOW_APL
// Written by Michal Jahnatek

#ifndef _AFLOW_APL_H_
#define _AFLOW_APL_H_

//#define USE_MKL 1

// Almost generally used precision in the whole apl, however, somewhere it is
// hard coded based on the tests and it work much better...
#define _AFLOW_APL_EPS_ 1E-6
extern bool _WITHIN_DUKE_;  //will define it immediately in kphonons

// Aflow core libraries
#include "../aflow.h"
#include "../AUROSTD/aurostd.h"
//The functions in the header file "aflow_qha_operations.h" include various vector and matrix operations.//
//These functions are used to calculate eigenvalues and eigenvectors of complex symmetric Hermitian matrices and .//
//Also perform other operations involved in calculating quasiharmonic properties.//
//These functions have been built based on vectorview concepts//
#include "aflow_qha_operations.h"
//#define _AFLOW_APL_REGISTER_ register   // register ?
#define _AFLOW_APL_REGISTER_

// Create the version of GCC, we will uset it for multithread parts of code,
// to check if the current compiling version of gcc is able to compile the
// std::thead features
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

// Basic objects ...
// ***************************************************************************
// "aplexcept.h"
// in aurostd.h // [OBSOLETE]
#include <stdexcept>
namespace apl {
//
class APLRuntimeError : public std::runtime_error {
 public:
  APLRuntimeError(const std::string& s) : std::runtime_error(s) {}
};
//
class APLLogicError : public std::logic_error {
 public:
  APLLogicError(const std::string& s) : std::logic_error(s) {}
};
//
class APLStageBreak : public std::exception {
 public:
  APLStageBreak() {}
};
}

// ***************************************************************************
// "logger.h"
namespace apl {
class Logger;
// Templates used for the definition of the one parameter manipulation
// functions of Logger's streams
template <class T>
class LMANIP;
template <class T>
Logger& operator<<(Logger&, const LMANIP<T>&);
template <class T>
class LMANIP {
 public:
  LMANIP(Logger& (*a)(Logger&, T), T v) {
    _action = a;
    _value = v;
  }
  friend Logger& operator<<<>(Logger&, const LMANIP<T>&);

 private:
  Logger& (*_action)(Logger&, T);
  T _value;
};
template <class T>
Logger& operator<<(Logger& s, const LMANIP<T>& m) {
  return (*m._action)(s, m._value);
}
// Function without parameter...
typedef Logger& (*MyStreamManipulator)(Logger&);
// this is the type of std::cout
typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
// this is the function signature of std::endl
typedef CoutType& (*StandardEndLine)(CoutType&);
//
class Logger {
 private:
  ofstream* _os;
  std::string _barCode;
  std::string _typeofMessage;
  std::string _moduleName;
  _aflags _aflowFlags;
  stringstream _ss;
  int _progressBarLastPercent;
  bool _isQuiet;

 public:
  Logger(std::ofstream&, const _aflags&);
  Logger(Logger&);
  ~Logger();
  ofstream& getOutputStream();
  void initProgressBar(const char*);
  void initProgressBar(const string&);
  void updateProgressBar(int, int);
  void finishProgressBar();
  void setTypeOfMessage(const string&);
  void setBarCode(const string&);
  void setBarCode(const char&);
  void setModuleName(const string&);
  void setQuietMode(bool);
  Logger& operator=(Logger&);
  Logger& operator<<(const string& s);
  Logger& operator<<(const int& s);
  Logger& operator<<(const uint& s);
  Logger& operator<<(const double& s);
  // Manipulation functions without parameter
  Logger& operator<<(StandardEndLine);
  Logger& operator<<(MyStreamManipulator);
  void endl();
  friend Logger& endl(Logger&);
  friend Logger& message(Logger&);
  friend Logger& notice(Logger&);
  friend Logger& warning(Logger&);
  friend Logger& error(Logger&);
  friend Logger& quiet(Logger&);
  // Manipulation functions with one parameter
  friend Logger& setwidth(Logger&, int);
  friend Logger& setformat(Logger&, const char*);
};
// Friend functions without parameters
Logger& endl(Logger&);
Logger& message(Logger&);
Logger& notice(Logger&);
Logger& warning(Logger&);
Logger& error(Logger&);
Logger& quiet(Logger&);
// Friend functions with one parameter
Logger& setwidth(Logger&, int);
LMANIP<int> sw(int);
Logger& setformat(Logger&, const char*);
LMANIP<const char*> sf(const char*);
}  // namespace APL

// ***************************************************************************
// "hroutines.h"
// in aurostd.h // [OBSOLETE]
#include <typeinfo>

namespace apl {
template <typename T>
inline std::string stringify(const T& x) {
  std::ostringstream o;
  if (!(o << x)) {
    throw APLRuntimeError(std::string("stringify(") + typeid(x).name() + ")");
  }
  return o.str();
}
void tokenize(const string&, vector<string>&, string);
string getVASPVersionString(const string&);
unsigned int getFileCheckSum(const string&);
void printXVector(const xvector<double>&, bool = true);
void printXVector(const xvector<xcomplex<double> >&);
void printXMatrix(const xmatrix<double>&);
void printXMatrix(const xmatrix<xcomplex<double> >&);
void printXMatrix2(const xmatrix<xcomplex<double> >&);
}  // namespace apl

// ***************************************************************************
// ***************************************************************************
// BEGIN JJPR: PAIRS AND TRIPLETS
// ***************************************************************************
#include <vector>

namespace apl {
class _pair {  //Pair of atoms: Needed for Thermal Conductivity
 public:
  _pair();                                 // default, just allocate
  ~_pair();                                // kill everything
  const _pair& operator=(const _pair& b);  // copy
  //Actuales
  uint indexA;           // ATOM INDEX (Primitive Cell)
  uint indexB;           // ATOM INDEX (Supercell)
  bool is_inequivalent;  // TRUE = INEQUIVAENT
  int equivalent;        // Index inequivalent pair
  // Arrays for TRIPLETS
  std::vector<uint> tri;        // ATOM INDEX
  std::vector<bool> tri_is_eq;  // TRUE = INEQUIVALENT
  std::vector<uint> tri_peq;    // Index Inequivalent pair
  std::vector<uint> tri_teq;    // Index Inequivalent triplet
  std::vector<int> tri_sym;     // Index Sym Operator
  std::vector<int> tri_perm;    // Index Permutation Operator

 private:
  void free();  // free space
};

class _itrip {  //Inequivalent triplets  of atoms: Needed for Thermal Conductivity
 public:
  _itrip();                                  // default, just allocate
  ~_itrip();                                 // kill everything
  const _itrip& operator=(const _itrip& b);  // copy
  uint tri_in;                               // Index triplet in pairs
  std::vector<uint> tri_eq;                  // Array of index for equivalent triplets (First values is inequivalent one)
  std::vector<uint> p_eq;                    // Array of index for equivalent pairs (First values is inequivalent one)
  std::vector<uint> tri_eq2;                 // Array of index for equivalent triplets (First values is inequivalent one)
  std::vector<uint> p_eq2;                   // Array of index for equivalent pairs (First values is inequivalent one)
  int dist_tri[3][3][3];
  //std::vector< vector< vector < vector<int> > > >    dist_tri_eq; // 3x3x3 Distortions tensor for inequivalent triplet
  std::vector<vector<double> > pool;   // Pool for different purposes
  std::vector<vector<double> > pool2;  // Pool for different purposes

  std::vector<vector<int> > triplets;                         //triplets
  std::vector<int> sym;                                       //symmetry operator
  std::vector<int> perm;                                      //permutation operator
  std::vector<vector<double> > rot_transform_aux;             //transform rot
  std::vector<int> rot_independent;                           //transform rot
  std::vector<int> count_independent;                         //transform rot
  std::vector<double> val_independent;                        //transform rot
  std::vector<double> valtmp_independent;                     //transform rot
  std::vector<vector<vector<double> > > rot_transform;        //transform rot
  std::vector<vector<vector<double> > > rot_transform_array;  //transform rot

 private:
  void free();  // free space
};
class _ipair {  // Inequivalent pairs of atoms: Needed for Thermal Conductivity
 public:
  _ipair();                                  // default, just allocate
  ~_ipair();                                 // kill everything
  const _ipair& operator=(const _ipair& b);  // copy
  int dist[6][6];                            // Distortions in pairs x/y/z/-x/-y/-z
  int dist_3x3[3][3];                        // Distortions in pairs x/y/z/-x/-y/-z
  bool dist_in[6][6];                        // TRUE = INEQUIVALENT DISTORTION
  std::vector<uint> eq;                      // Arrays with equivalent pairs
  std::vector<_itrip> itriplets;             // WARNING: we use starting from 0
  bool itriplets_exist;                      // TRUE = There are INEQUIVALENT TRIPLETS
 private:
  void free();  // free space
};

class StrPairs {
 private:
  Logger& _logger;

 public:
  StrPairs(const xstructure, const xstructure, const vector<int>&, const vector<int>&, int, double&, Logger&);  // default, just allocate
  StrPairs(const StrPairs&);                                                                                    // default, just allocate
  StrPairs(Logger&);                                                                                            // VOID
  ~StrPairs();                                                                                                  // kill everything
  void clear();
  const StrPairs& operator=(const StrPairs&);  // copy

  xstructure a;
  xstructure pcell;
  vector<int> _pc2scMap;
  vector<int> _sc2pcMap;
  vector<int> pc2scMap;
  vector<int> sc2pcMap;
  double cells;                              // Number of unit cells in the supercell
  double cutoff;                             // Cutoff Radious
  uint Sc_atoms;                             // Number of Atoms in supercell
  vector<xvector<double> > ortoDistorsions;  // Distorsions y cartesian coordinates
  std::vector<_pair> pairs;                  // WARNING: we use starting from 0
  std::vector<_ipair> ipairs;                // equivalent/inequivalent atoms lookup table
  std::vector<_pair> pairsF;                 // WARNING: we use starting from 0
  std::vector<_ipair> ipairsF;               // equivalent/inequivalent atoms lookup table
  std::vector<_itrip> itriplets;             // Triplets list sort by inequivalents

  vector<xvector<int> > tri_index;
  vector<vector<xvector<double> > > tri_vec;

  void build(xstructure pc, xstructure a);
  void Triplets();
  void calculateTiplets_MPI(xstructure ai, uint, uint);
  void calculateTriDist_MPI(xstructure ai, uint, uint);
  int getU(const xstructure& b, int ipa, int atom1, int atom2, int vec1, int vec2);
  int getAtom(const xstructure& b, int atom1, int Op);
  double getShellDist(const xstructure& b, int SHELL, double cutoff);
  double getMinDist(const xmatrix<double>& lattice, xvector<double>& a, xvector<double>& b);
  void test(const xmatrix<double>& lattice, xvector<double>& a, xvector<double>& b, double);
  bool MATCH2(int, int, int, int, int);
  bool MATCH3(int, int, int, int, int, int, int);
  bool MATCH3_2(int, int, int, int, int, int, int, xmatrix<int>&);
  bool MATCH2_2(int, int, int, int, int, xmatrix<int>&);
  bool MATCH3v(int, int, int, int, int, int, int, int, int, int, int, int, int, int&);
  bool MATCH2v(int, int, int, int, int, int, int, int&);
  bool MATCH3v2(int, int, int, int, int, int, int, int, int, int, xmatrix<int>, int&);
  bool MATCH2v2(int, int, int, int, int, int, int, int, int, int&);
  void SymUtils();
  void gaussAAPL(int, vector<vector<double> > mat);
  void moveTri(vector<int>&);
  double orth(int, int, int);
  bool compareTri(vector<int>, vector<int>);
  vector<vector<int> > permut;
  vector<vector<int> > id_equi;
  vector<vector<vector<vector<double> > > > rot, rot2;
  vector<vector<vector<double> > > nonzero;
  void get3x3();
};

class _triBTE {  //Pair of atoms: Needed for Thermal Conductivity
 public:
  _triBTE();                                   // default, just allocate
  ~_triBTE();                                  // kill everything
  const _triBTE& operator=(const _triBTE& b);  // copy

  xvector<double> v1;
  xvector<double> v2;
  xvector<int> index;
  xvector<int> iat;
  vector<xmatrix<double> > ifcs;

 private:
  void free();  // free space
};

}  //namespace apl

// ***************************************************************************
// END JJPR: PAIRS AND TRIPLETS
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// Linear and nonlinear fitting functions, this functions compute fitting parameters correctly
//although chi-square matrix is not computing correctly [PINKU]
#define FIT_LIMIT 0.001
namespace apl {
class aflowFITTING {
 public:
  aflowFITTING() {}
  ~aflowFITTING() { clear(); }
  void clear() {}

 private:
  double chisq_quadratic;
  double chisq_birch_murnaghan;
  uint Iteration_birch_murnaghan;
  double alamda_birch_murnaghan;
  vector<double> Uncertainties_birch_murnaghan;
  //user defined fitting functions
  void funcs(double x, xvector<double>& afunc);
  void birch_murnaghan_function(double x,
                                const xvector<double> a,
                                double* y,
                                xvector<double>& dyda);

  //linear leatsquare fitting function
  template <class utype>
  bool lfit(xvector<utype>& x,
            xvector<utype>& y,
            xvector<utype>& sig,
            xvector<utype>& a,
            xvector<int>& ia,
            xmatrix<utype>& covar,
            utype& chisq,
            void (aflowFITTING::*funcs)(utype, xvector<utype>&));

  //nonlinear leatsquare fitting function
  bool mrqmin(xvector<double> x,
              xvector<double> y,
              xvector<double>& sig,
              xvector<double>& a,
              xvector<int>& ia,
              xmatrix<double>& covar,
              xmatrix<double>& alpha,
              double& chisq,
              void (aflowFITTING::*birch_murnaghan_function)(double, xvector<double>, double*, xvector<double>&),
              double* alamda);

  template <class utype>
  void covsrt(xmatrix<utype>& covar, xvector<int>& ia, int& mfit);
  template <class utype>
  bool gaussj(xmatrix<utype>& a, int& n, xmatrix<utype>& b, int m);

  void mrqcof(xvector<double> x,
              xvector<double> y,
              xvector<double>& sig,
              xvector<double>& a,
              xvector<int>& ia,
              xmatrix<double>& alpha,
              xvector<double>& beta,
              double& chisq,
              void (aflowFITTING::*birch_murnaghan_function)(double, xvector<double>, double*, xvector<double>&));

 public:
  bool
  quadraticfit(xvector<double> energy, xvector<double> volume, xvector<double>& params);
  bool
  birch_murnaghan_fitting(xvector<double> energy, xvector<double> volume, xvector<double> guess, xvector<double>& params);
  double get_chisq_quadratic() { return chisq_quadratic; }
  double get_chisq_birch_murnaghan() { return chisq_birch_murnaghan; }
  uint getIteration_birch_murnaghan() { return Iteration_birch_murnaghan; }
  double getalamda_birch_murnaghan() { return alamda_birch_murnaghan; }
  vector<double> getUncertainties_birch_murnaghan() { return Uncertainties_birch_murnaghan; }
};
}
// ***************************************************************************
namespace apl {
#ifndef EIGEN_H
#define EIGEN_H
//Functions in this class calculate eigenvalues and eigenvectors of complex symmetrix
//nxn matrices and can sort them according to following options
//eval and evec SORTING OPTIONS
// [ 1. [APL_MV_EIGEN_SORT_VAL_ASC]  => ascending order   ]
// [ 2. [APL_MV_EIGEN_SORT_VAL_DESC] => descending order ]
// [ 3. [APL_MV_EIGEN_SORT_ABS_ASC]  => absolute ascending order ]
// [ 4. [APL_MV_EIGEN_SORT_ABS_DESC] => absolute descending order]
//Compute Eigenvalues and Eigenvectors for (nxn) Complex Hermitian matrices
class aplEigensystems : public MVops {
 public:
  aplEigensystems() {}
  ~aplEigensystems() {}
  void clear() { } //this->clear(); }
  void eigen_calculation(const aurostd::xmatrix<xcomplex<double> >& M, aurostd::xvector<double>& eval, aurostd::xmatrix<xcomplex<double> >& evec);
  void eigen_calculation(const xmatrix<xcomplex<double> >& M, xvector<double>& eval, xmatrix<xcomplex<double> >& evec, apl_eigen_sort_t t);
};
#endif
}
// ***************************************************************************
namespace apl {
//The functions in this class are used to perform linear and nonlinear fittings//
//User defined functions can be used to perform fitting//
//Both Levenberg-Marquardt and unscaled Levenberg-Marquardt algorithms can be performed
//either with keyword "APL_multifit_fdfsolver_lmsder" or "APL_multifit_fdfsolver_lmder" //
#ifndef FIT_H
#define FIT_H
#define spread 1e-6               //accuracy of data fitting
#define allowed_fit_error 100.00  // if spread is small then allowed_fit_error is relatively take high value
#define Absolute_Error 1E-6
#define Relative_Error 1E-6
#define APL_multifit_fdfsolver_lmsder  //This is a robust and efficient version of the Levenberg-Marquardt algorithm
#undef APL_multifit_fdfsolver_lmder    //This is an unscaled version of the Levenberg-Marquardt algorithm

class md_lsquares : public MVops {
  // multi-parameter linear regression
  // multidimensional nonlinear least-squares fitting
 public:
  md_lsquares() {}
  ~md_lsquares() {}

  //user define data sets
  vector<double> Xdata;  //usder defined input data
  vector<double> Ydata;  //usder defined input data
  //guess values calculated automatically
  vector<double> guess;  // guess values for nonlinear fitting, It will calculate automatically.
  //error managments
  bool data_read_error;   //return true if Xdata.size and Ydata.size are same
  int nl_success_status;  //return error message interms of int numbers
  string nl_err_msg;      //return type of error mesages

  //lfit output
  double luncertanity_V0;  //uncertainties linear fit(Quadratic function) in Volume
  double luncertanity_E0;  //uncertainties linear fit(Quadratic function) in Energies
  double luncertanity_B0;  //uncertainties linear fit(Quadratic function) in Bulk Modulus
  double lchisq;           //linear fit chisquare
  double leqmV0;           //Equilibrium Volume from linear fit (Quadratic function)
  double leqmE0;           //Equilibrium Energy from linear fit (Quadratic function)
  double leqmB0;           //Equilibrium Bulk Modulus from linear fit (Quadratic function)
  //nlfit output
  double uncertanity_V0;   //uncertainties non linear fit(Birch Murnaghan Function) in Volume
  double uncertanity_E0;   //uncertainties non linear fit(Birch Murnaghan Function) in Energy
  double uncertanity_B0;   //uncertainties non linear fit(Birch Murnaghan Function) in Bulk Modulus
  double uncertanity_Bp;   //uncertainties non linear fit(Birch Murnaghan Function) in Pressure derivatives Bulk Modulus
  double uncertanity_Bpp;  //uncertainties non linear fit(Birch Murnaghan Function) in Pressure derivatives Bp
  double chisq_dof;        //chi square per degrees of freedom
  double nleqmV0;          //Equilibrium Volume from nonlinear fit (Birch Murnaghan Function)
  double nleqmE0;          //Equilibrium Energy from nonlinear fit (Birch Murnaghan Function)
  double nleqmB0;          //Equilibrium Bulk Modulus from nonlinear fit (Birch Murnaghan Function)
  double nleqmBp;          //Equilibrium Pressure derivatives Bulk Modulus from nonlinear fit (Birch Murnaghan Function)
  double nleqmBpp;         //Equilibrium Pressure derivatives Bp from nonlinear fit (Birch Murnaghan Function)
  string fdfsolver_name;   //return nonlinear method name
  //fitting functions
  void cubic_polynomial_fit();                                            //cubic linear function fitting
  void birch_murnaghan_fit();                                             //Birch-Murnaghan nonlinear function fitting
  void birch_murnaghan_4th_order_fit(const xvector<double>& user_guess);  //Birch-Murnaghan 4th order nonlinear function fitting
  void birch_murnaghan_3rd_order_fit(const xvector<double>& user_guess);  //Birch-Murnaghan 3rd order nonlinear function fitting
  void clear();
};
#endif
}
// ***************************************************************************
// ***************************************************************************
// Engine core for phonon calculation
// ***************************************************************************
// "supercell.h"
namespace apl {
class Supercell {
 private:
  Logger& _logger;
  xstructure _inStructure;
  xstructure _inStructure_original;  //CO
  xstructure _inStructure_light;     //CO, does not include HEAVY symmetry stuff
  //deque<_atom> _inStructure_atoms_original; //CO
  xstructure _pcStructure;           //CO 180406 - for the path
  xstructure _scStructure;
  xstructure _scStructure_original;  //CO
  xstructure _scStructure_light;     //CO, does not include HEAVY symmetry stuff
  //deque<_atom> _scStructure_atoms_original; //CO
  //CO - START
  bool _skew;                  //SYM::isLatticeSkewed(), same for pc and sc
  bool _derivative_structure;  //vs. simple expanded_lattice, derivative structure's lattice has LESS symmetry, so be careful ApplyAtom()'ing
  double _sym_eps;             //same for pc and sc
  //CO - END
  bool _isShellRestricted;
  int _maxShellID;
  vector<double> _maxShellRadius;
  bool _isConstructed;

 private:
  void calculateWholeSymmetry(xstructure&);
  xstructure calculatePrimitiveStructure() const;

 public:
  Supercell(const xstructure&, Logger&);
  Supercell(const Supercell&);
  ~Supercell();
  void LightCopy(const xstructure& a, xstructure& b);
  void clear();
  Supercell& operator=(const Supercell&);
  bool isConstructed();
  void reset();
  void build(int, int, int, bool = TRUE);
  void trimStructure(int, const xvector<double>&,
                     const xvector<double>&, const xvector<double>&,
                     bool = true);
  int buildSuitableForShell(int, bool, bool VERBOSE);
  void setupShellRestrictions(int);
  bool isShellRestricted();
  int getMaxShellID();
  int getNumberOfAtoms();
  int getNumberOfUniqueAtoms();
  int getNumberOfUniqueAtomsOfType(int);
  int getUniqueAtomID(int);
  int getUniqueAtomID(int, int);
  const _atom& getUniqueAtom(int);
  string getUniqueAtomSymbol(int);
  double getUniqueAtomMass(int);
  double getAtomMass(int);
  int getAtomNumber(int);
  const xstructure& getSupercellStructure() const;
  const xstructure& getSupercellStructureLight() const;
  const xstructure& getPrimitiveStructure() const;
  const xstructure& getInputStructure() const;
  const xstructure& getInputStructureLight() const;
  int where(int, int, const _sym_op&, bool = TRUE);
  int wherefrom(int, int, const _sym_op&, bool = TRUE);
  const _sym_op& getSymOpWhichMatchAtoms(int, int, int);
  xvector<double> getFPositionItsNearestImage(const xvector<double>&,
                                              const xvector<double>&,
                                              const xmatrix<double>&);
  xvector<double> getFPositionItsNearestImage(int, int);
  xvector<double> getCPositionItsNearestImage(int, int);
  bool compareFPositions(xvector<double>&, xvector<double>&);          //CO
  bool compareFPositions(xvector<double>&, xvector<double>&, double);  //CO
  bool calcShellPhaseFactor(int, int, const xvector<double>&, xcomplex<double>&);
  int pc2scMap(int);
  int sc2pcMap(int);
  void center(int);
  //CO - START
  void center_original(void);
  //corey
  vector<vector<_sym_op> > getAGROUP(void);
  vector<_sym_op> getFGROUP(void);
  vector<_sym_op> getAGROUP(int);
  bool isDerivativeStructure();
  double getEPS(void);
  // CO - END
  // **** BEGIN JJPR *****
  xvector<int> scell;
  xcomplex<double> calcShellPhaseFactor2(int, int, const xvector<double>&);
  bool calcShellPhaseFactorAAPL(int, int, const xvector<double>&, xcomplex<double>&, xvector<xcomplex<double> >&);  //JJPR modification
  vector<xvector<double> > Rl_list;
  vector<xvector<double> > rl_list;
  void calcRl();
  friend class ThermalConductivityCalculator;
  friend class StrPairs;
  vector<int> _pc2scMap;
  vector<int> _sc2pcMap;
  // **** END  JJPR *****
};
}  // namespace apl

// ***************************************************************************
// "ipc.h"

namespace apl {
enum IPCFreqFlags {
  NONE = 0L,
  ALLOW_NEGATIVE = 1L << 1,
  OMEGA = 1L << 2,
  RAW = 1L << 3,  // eV/A/A/atomic_mass_unit
  HERTZ = 1L << 4,
  THZ = 1L << 5,
  RECIPROCAL_CM = 1L << 6,
  MEV = 1L << 7
};
inline IPCFreqFlags operator&(const IPCFreqFlags& __a, const IPCFreqFlags& __b) {
  return IPCFreqFlags(static_cast<int>(__a) & static_cast<int>(__b));
}
inline IPCFreqFlags operator|(const IPCFreqFlags& __a, const IPCFreqFlags& __b) {
  return IPCFreqFlags(static_cast<int>(__a) | static_cast<int>(__b));
}
inline IPCFreqFlags operator|=(IPCFreqFlags& __a, const IPCFreqFlags& __b) {
  return (__a = (__a | __b));
}
// //////////////////////////////////////////////////////////////////////////

class IPhononCalculator {
 public:
  virtual ~IPhononCalculator() {}
  virtual xvector<double> getFrequency(const xvector<double>&, IPCFreqFlags) = 0;
  virtual double getEPS() = 0;  //CO
  // **** BEGIN JJPR ******
  virtual xvector<double> getFrequencyAAPL(const xvector<double>&, xmatrix<xcomplex<double> >&, vector<xmatrix<xcomplex<double> > >&) = 0;  // Modified by JJPR
  virtual double getElementTensor(uint, uint, uint, uint, uint, uint) = 0;                                                                  // Modified by JJPR
  // **** END JJPR ****
  virtual double getFrequencyConversionFactor(IPCFreqFlags, IPCFreqFlags) = 0;
  virtual const Supercell& getSupercell() = 0;
  virtual const xstructure& getInputCellStructure() = 0;
  virtual const xstructure& getSuperCellStructure() = 0;
  virtual uint getNumberOfBranches() = 0;
  // **** BEGIN PINKU ******
  virtual xmatrix<xcomplex<double> > getDynamicMatrix(const xvector<double>) = 0;
  virtual vector<double> get_ATOMIC_MASSES_AMU() = 0;
  virtual void get_special_inputs(string& AflowIn) = 0;
  virtual void get_NCPUS(string& ncpus) = 0;  //CO 180214
  virtual void get_NCPUS(int& ncpus) = 0;     //CO 180214
  // **** END   PINKU ******
};
}  // namespace apl

// ***************************************************************************
// "phoncalc.h"
#define _AFLOW_APL_AFLOW_DIRECTORY_PREFIX_ string("ARUN.")
#define _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ string(_AFLOW_APL_AFLOW_DIRECTORY_PREFIX_ + "APL0LRBE")
#define _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_ string(_AFLOW_APL_AFLOW_DIRECTORY_PREFIX_ + "APL1LRFF")

namespace apl {
class PhononCalculator : virtual public IPhononCalculator {
 protected:
  // USER PARAMETERS
  double DISTORTION_MAGNITUDE;
  bool KAPPA;        // JJPR: TENSOR OPTION
  double THRESHOLD;  // JJPR: SUMRULE THERSHOLD
  // Aflow's stuff required for running some routines
  _xinput& _xInput; //_xvasp& _vaspRun;
  _aflags& _aflowFlags;
  _kflags& _kbinFlags;
  _xflags& _xFlags; //_vflags& _vaspFlags;
  string& _AflowIn;

  // Our standard and common output system
  Logger& _logger;
  // To tar or not to tar
  bool DOtar;
  // The computational model, supercell (now used by both implementations)
  Supercell& _supercell;
  //******BEGIN JJPR***********
  // Class for pairs info
  StrPairs& _strpair;
  // All runs from Tensor will be stored here
  vector<_xinput> xInputsT; //vector<_xvasp> vaspRunsT;
  vector<_xinput> xInputs; //vector<_xvasp> vaspRuns;
  //Forces for all atoms with 36 distortions
  vector<vector<vector<vector<xvector<double> > > > > bigTensor;
  //*******END JJPR************
  // For each inequivalent atom, there is a set of unique distortions
  vector<vector<xvector<double> > > _uniqueDistortions;
  // For each inequivalent atom and unique distortion, there is a field
  // of forces (for each atom of the supercell)
  vector<vector<vector<xvector<double> > > > _uniqueForces;
  // For each atom of supercell, there is a full force field
  vector<vector<xmatrix<double> > > _forceConstantMatrices;
  // Calculate forces at no distortion - since for some structure
  // (not well relaxed, or with other problems) these forces have to be
  // known and substracted from the calculated forces with distortion
  bool _calculateZeroStateForces;

  vector<double> ATOMIC_MASSES_AMU;  //[PINKU]
  string _check_LDAU2_ON;            //PINKU
  string _LDAU_PARAMETERS;           //PINKU
  string _PSTRESS;                   //PINKU

  // Stuff for polar materials
  bool _isPolarMaterial;
  // For each atom there is a matrix 3x3 of Born effective charge
  vector<xmatrix<double> > _bornEffectiveChargeTensor;
  // Dielectric tensor
  xmatrix<double> _dielectricTensor;
  // Precomputed values used in non-analytical term
  xmatrix<double> _inverseDielectricTensor;
  double _recsqrtDielectricTensorDeterminant;
  // Precomputed Ewald sum at Gamma point
  bool _isGammaEwaldPrecomputed;
  vector<xmatrix<xcomplex<double> > > _gammaEwaldCorr;

 private:
  virtual void calculateForceFields() {}
  // *********BEGIN JJPR***********
  void VaspTensor();
  void StoreForcesTensor();
  void buildTensor();
  void SumRuleTensor();
  double CHECK();
  void CHECK2();
  void MIX();
  void RESYM();
  void RESYMM();
  void PrintTensor(string);
  void ReadTensor();
  void PrintTensorBTE(string);
  void get2ndF();
  double CHECKM();
  void SumRuleM();
  xmatrix<xcomplex<double> > getDynamicMatrixAAPL(const xvector<double>, vector<xmatrix<xcomplex<double> > >&);          // Modified by JJPR
  xmatrix<xcomplex<double> > getDynamicMatrixAAPL2(const xvector<double>, vector<xmatrix<xcomplex<double> > >&);         // Modified by JJPR
  xmatrix<xcomplex<double> > getNonanalyticalTermWangAAPL(const xvector<double>, vector<xmatrix<xcomplex<double> > >&);  // Modified JJPR
  // ***********END JJPR*********
  void completeForceFields();
  void projectToCartesianDirections();
  void buildForceConstantMatrices();
  void symmetrizeForceConstantMatrices();
  void correctSumRules();
  void printForceConstantMatrices(ostream&);
  void printFCShellInfo(ostream&);
  xmatrix<xcomplex<double> > getDynamicMatrix(const xvector<double>);
  xmatrix<xcomplex<double> > getNonanalyticalTermWang(const xvector<double>);
  xmatrix<xcomplex<double> > getNonanalyticalTermGonze(const xvector<double>);
  xmatrix<xcomplex<double> > getEwaldSumDipolDipolContribution(const xvector<double>, bool = true);
  vector<double> get_ATOMIC_MASSES_AMU() { return ATOMIC_MASSES_AMU; }  //[PINKU]
  void store_masses();                                                  //[PINKU]
 protected:
  void writeOUTPUT(_xinput&);
  //void createAIMSOUTPUT(const _xaims&); //CO 180406 - obsolete
  //void createAFLOWIN(const _xvasp&);    //CO 180406 - obsolete

 public:
  PhononCalculator(Supercell&, StrPairs&, _xinput&, _aflags&, _kflags&, _xflags&, string&, Logger&); //_xvasp&, _aflags&, _kflags&, _vflags&, Logger&);  // Modified by JJPR
  virtual ~PhononCalculator();
  void clear();
  void run();
  xvector<double> getEigenvalues(const xvector<double>&);
  xvector<double> getEigenvaluesAAPL(const xvector<double>&, xmatrix<xcomplex<double> >&, vector<xmatrix<xcomplex<double> > >&);  // Modified by JJPR
  void isPolarMaterial(bool b) { _isPolarMaterial = b; }
  void setDistortionMagnitude(double f) { DISTORTION_MAGNITUDE = f; }
  void setCalculateZeroStateForces(bool b) { _calculateZeroStateForces = b; }
  void setTensor(bool b) { KAPPA = b; }         // JJPR
  void setSumRule(double f) { THRESHOLD = f; }  // JJPR
  void hibernate();
  void awake();
  void writeFORCES();
  void writeDYNMAT();
  void writeXCrysDenForces();
  //*******BEGIN JJPR*****
  vector<xvector<double> > ortoDistorsions;
  vector<double> Ineq_Forc_In;
  vector<double> Ineq_Forc_Fin;
  //Forces for all atoms with 9 distortions
  vector<vector<vector<vector<vector<xvector<double> > > > > > smallTensor;
  vector<vector<vector<vector<vector<xvector<double> > > > > > smallTensor_sumed;
  vector<vector<vector<vector<vector<xvector<double> > > > > > smallTensorOLD;
  vector<vector<vector<xvector<double> > > > smallMatrix;
  //******END JJPR*******
  // Interface
  xvector<double> getFrequencyAAPL(const xvector<double>&, xmatrix<xcomplex<double> >&, vector<xmatrix<xcomplex<double> > >&);  //Modified by JJPR
  xvector<double> getFrequency(const xvector<double>&, IPCFreqFlags);
  double getFrequencyConversionFactor(IPCFreqFlags, IPCFreqFlags);
  const Supercell& getSupercell();
  const xstructure& getInputCellStructure();
  const xstructure& getSuperCellStructure();
  double getEPS();  //CO
  uint getNumberOfBranches();
  /* friend void runVASPCalculationsBE(apl::PhononCalculator*); */
  /* friend void readBornEffectiveChargesFromOUTCAR(apl::PhononCalculator *pcalculator); */
  /* friend void symmetrizeBornEffectiveChargeTensors(apl::PhononCalculator *pcalculator); */
  /* friend void readDielectricTensorFromOUTCAR(apl::PhononCalculator *pcalculator); */
  void runVASPCalculationsBE(void);
  void readBornEffectiveChargesFromAIMSOUT(void);
  void readBornEffectiveChargesFromOUTCAR(void);
  void symmetrizeBornEffectiveChargeTensors(void);
  void readDielectricTensorFromAIMSOUT(void);
  void readDielectricTensorFromOUTCAR(void);
  double getElementTensor(uint, uint, uint, uint, uint, uint);  //JJPR
  void get_special_inputs(string& AflowIn);                                    //PINKU
  void get_NCPUS(string& ncpus);  //CO 180214
  void get_NCPUS(int& ncpus);     //CO 180214
  // BEGIN ME 180518
  bool filesExistPhonons(_xinput&);
  void createAflowInPhonons(_xinput&, const string&);
  bool outfileFoundAnywherePhonons(vector<_xinput>&);
  void outfileFoundEverywherePhonons(vector<_xinput>&);
  void subtractZeroStateForces(vector<_xinput>&);
  // END ME 180518
};
}  // namespace apl

// ***************************************************************************
// "dirphoncalc.h"
namespace apl {
class DirectMethodPC : public PhononCalculator {
 protected:
  //bool GENERATE_PLUS_MINUS;  //JAHNATEK ORIGINAL
  //CO - START
  bool AUTO_GENERATE_PLUS_MINUS;
  bool USER_GENERATE_PLUS_MINUS;
  //CO - END
  bool GENERATE_ONLY_XYZ;

 protected:
  void estimateUniqueDistortions(const xstructure&,
                                 vector<vector<xvector<double> > >&);
  void testDistortion(const xvector<double>&, const vector<_sym_op>&,
                      vector<xvector<double> >&,
                      vector<xvector<double> >&);
  bool needMinus(uint inq_atom_indx, uint distortion_indx);  //CO
  void runVASPCalculations();
  // BORN
  /* void runVASPCalculationsBE(); */
  /* void readBornEffectiveChargesFromOUTCAR(); */
  /* void symmetrizeBornEffectiveChargeTensors(); */
  /* void readDielectricTensorFromOUTCAR(); */

 public:
  DirectMethodPC(Supercell&, StrPairs&, _xinput&, _aflags&, _kflags&, _xflags&, string&, Logger&); //_xvasp&, _aflags&, _kflags&, _vflags&, Logger&);  //Modified by JJPR
  ~DirectMethodPC();
  void clear();
  void calculateForceFields();
  // Easy access to global parameters
  //void setGeneratePlusMinus(bool b) { GENERATE_PLUS_MINUS = b; } //JAHNATEK ORIGINAL
  void setGeneratePlusMinus(bool _auto_, bool _user_) {
    AUTO_GENERATE_PLUS_MINUS = _auto_;
    USER_GENERATE_PLUS_MINUS = _user_;
  }  //CO
  void setGenerateOnlyXYZ(bool b) { GENERATE_ONLY_XYZ = b; }
};
}  // namespace apl

//PINKU - START
// ***************************************************************************
#define _GP_AFLOW_DIRECTORY_PREFIX_ string("ARUN.APL_PHONON_")
#define _EOS_AFLOW_DIRECTORY_PREFIX_ string("ARUN.APL_STATIC_")
#define _PH_AFLOW_DIRECTORY_PREFIX_ string("ARUN.APL_PHONON_")
#define _aflowinpad_ 54
#define SetPrecision 15
//Functions in this class are used to calculate equation of state properties
namespace apl {
class eos_calculate : public PhononCalculator {
 private:
  vector<string> GP_dir_names;
  vector<string> PH_dir_names;
  vector<string> EOS_dir_names;
  vector<double> GP_volumes;
  vector<double> EOS_Volumes;
  int ZERO_STATIC_DIR_INDEX;
  string _PSTRESS;

 public:
  eos_calculate(Supercell& sc, StrPairs&, _xinput& xinput, //_xvasp& xvasp,
                _aflags& aflags, _kflags& kflags,
                _xflags& xflags, //_vflags& vflags, 
                string&,
                Logger& l);
  ~eos_calculate();
  void clear();
  void GP_RUN();                                           //VASP run for Gruneisen parameter
  void EOS_RUN(string&,const double, const double, const double);  //VASP run for EOS
 private:
  bool GRUNEISEN, EOS;
  double GP_VOL_DISTORTION, EOS_VOL_DISTORTION;
  int EOS_VOL_START, EOS_VOL_END;
  string EOS_KPPRA;
  string EOS_STATIC_KPPRA;
  string EOS_BANDS_GRID;
  string NEDOS;
  string EOS_KSCHEME;
  string EOS_STATIC_KSCHEME;
  void writeGPOUTPUT(_xinput& xinput, const bool _isGP);
  //void createGPAIMSOUTPUT(const _xaims& xaims_input);                //CO 180406 - obsolete
  //void createGPAFLOWIN(const _xvasp& xvasp_input, const bool _isGP); //CO 180406 - obsolete
  void writeEOSOUTPUT(const _xinput& xinput);
  void createEOSAIMSOUTPUT(const _xaims& xaims_input);
  void createEOSAFLOWIN(const _xvasp& xvasp_input);
  void get_special_inputs(string& AflowIn);

 public:
  void setGruneisen(bool b) { GRUNEISEN = b; }
  void setEOS(bool b) { EOS = b; }
  void setGP_VOL_DISTORTION(double f) { GP_VOL_DISTORTION = f; }
  void setTEC_VOL_DISTORTION(double f) { EOS_VOL_DISTORTION = f; }
  void setEOS_VOL_START(int Int) { EOS_VOL_START = Int; }
  void setEOS_VOL_END(int Int) { EOS_VOL_END = Int; }
  void setEOS_KPPRA(std::string s) { EOS_KPPRA = s; }
  void setEOS_STATIC_KPPRA(std::string s) { EOS_STATIC_KPPRA = s; }
  void setEOS_BANDS_GRID(std::string s) { EOS_BANDS_GRID = s; }
  void setNEDOS(std::string s) { NEDOS = s; }
  template <typename T>
  string NumberToString(T Number);
  vector<string> getGP_dir_names() { return GP_dir_names; }
  vector<string> getPH_dir_names() { return PH_dir_names; }
  vector<string> getEOS_dir_names() { return EOS_dir_names; }
  vector<double> getGP_volumes() { return GP_volumes; }
  vector<double> getEOS_Volumes() { return EOS_Volumes; }
  int getZERO_STATIC_DIR_INDEX() { return ZERO_STATIC_DIR_INDEX; }
};
}  // namespace apl
//PINKU - END

// ***************************************************************************
// "lrphoncalc.h"
namespace apl {
class LinearResponsePC : public PhononCalculator {
 private:
  /* void runVASPCalculationsBE(); */
  void runVASPCalculationsFF();
  /* void readBornEffectiveChargesFromOUTCAR(); */
  /* void symmetrizeBornEffectiveChargeTensors(); */
  /* void readDielectricTensorFromOUTCAR(); */
  void readForceFieldsFromDYNMAT();

 public:
  LinearResponsePC(Supercell&, StrPairs&, _xinput&, _aflags&, _kflags&, _xflags&, string&, Logger&); //_xvasp&, _aflags&, _kflags&, _vflags&, Logger&);  //Modified by JJPR
  ~LinearResponsePC();
  void clear();
  void calculateForceFields();
};
}  // namespace apl

// ***************************************************************************
// "gsa.h"

//CO generally redirects to DM, the distinction between DM and GSA is obsolete
namespace apl {
class GeneralizedSupercellApproach : public DirectMethodPC {
 public:
  GeneralizedSupercellApproach(Supercell&, StrPairs&, _xinput&, //_xvasp&, 
      _aflags&, _kflags&, _xflags&, //_vflags&, 
      string&,
      Logger&);  //Modified by JJPR
  ~GeneralizedSupercellApproach();
  void clear();
};
}  // namespace apl

// ***************************************************************************

//PINKU - START
#define AFLOW_APL_VASP_USE_LEPSILON
#undef AFLOW_APL_VASP_USE_LCALCEPS  // HAS some problem [PINKU]
//PINKU - END

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// Supplementary classes for calculation of dispersion curves and density of states

// ***************************************************************************
// "pathbuilder.h"
namespace apl {
class PathBuilder {
 public:
  enum StoreEnumType { RECIPROCAL_LATTICE,
                       CARTESIAN_LATTICE };
  enum ModeEnumType { SINGLE_POINT_MODE,
                      COUPLE_POINT_MODE };

 private:
  std::vector<aurostd::xvector<double> > _path;
  std::vector<aurostd::xvector<double> > _points;
  std::vector<std::string> _labels;
  aurostd::xmatrix<double> reciprocalLattice;
  aurostd::xmatrix<double> cartesianLattice;
  int _pointsVectorDimension;
  int _pointsVectorStartingIndex;
  uint _nPointsPerSubPath;
  ModeEnumType _mode;
  StoreEnumType _store;

 private:
  void buildPath();
  void tokenize(const string&, vector<string>&, string);

 public:
  PathBuilder();
  PathBuilder(ModeEnumType);
  ~PathBuilder();
  void clear();
  void addPoint(const std::string& l, int dim, ...);
  void addPoint(const std::string&, const aurostd::xvector<double>&);
  void transform(const aurostd::xmatrix<double>&);
  void pointsAreDefinedFor(const xstructure&, StoreEnumType);
  void transformPointsFor(const xstructure&, StoreEnumType);
  void defineCustomPoints(const string&,const string&,const Supercell&,bool CARESTIAN_COORDS=false);
  void takeAflowElectronicPath(const string&,const Supercell&);//, const xstructure&, const xstructure&);
  void setMode(ModeEnumType);
  void setStore(StoreEnumType);
  void setDensity(int);
  int getDensity();
  uint getPathSize();
  uint getPointSize();
  aurostd::xvector<double> getPoint(uint);
  uint getPointIndexOnPath(uint);
  std::string getPointLabel(uint);
  std::vector<aurostd::xvector<double> > getPath();
  std::vector<aurostd::xvector<double> > getPath(ModeEnumType, const string&);
  double getPathLength();
  double getPathLength(uint);
};
}  // namespace apl

// ***************************************************************************
// "pdisc.h"
namespace apl {
class PhononDispersionCalculator {
 private:
  IPhononCalculator& _pc;
  Logger& _logger;
  PathBuilder _pb;
  std::vector<xvector<double> > _qpoints;
  std::vector<xvector<double> > _freqs;
  IPCFreqFlags _frequencyFormat;
  vector<double> path;       //[PINKU]
  vector<int> path_segment;  //[PINKU]
 private:
  void calculateInOneThread(int, int);

 public:
  PhononDispersionCalculator(IPhononCalculator&, Logger&);
  ~PhononDispersionCalculator();
  void clear();
  void initPathCoords(const string&,const string&,int,bool=false);  //CO 180406
  void initPathLattice(const string&, int);
  void setPath(const string&);
  void calc(const IPCFreqFlags);
  void writePDIS();
  bool isExactQPoint(const xvector<double>&, const xmatrix<double>&);
  std::vector<xvector<double> > get_qpoints() { return _qpoints; }  //[PINKU]
  std::vector<double> get_path() { return path; }                   //[PINKU]
  std::vector<int> get_path_segment() { return path_segment; }      //[PINKU]
};
}  // namespace apl

// ***************************************************************************
// "irpg.h"
namespace apl {
class IReciprocalPointGrid {
 public:
  virtual ~IReciprocalPointGrid() {}
  virtual int getN(int) = 0;
  virtual aurostd::xmatrix<double> getReciprocalLattice() = 0;
  virtual std::vector<aurostd::xvector<double> > getPoints() = 0;
  virtual std::vector<double> getWeights() = 0;
  virtual aurostd::xtensor3<int> getAllToIrrPointMap() = 0;
};
}  // end namespace apl

// ***************************************************************************
// "mpmesh.h"
namespace apl {
class MonkhorstPackMesh : virtual public IReciprocalPointGrid {
 private:
  Logger& _logger;
  aurostd::xvector<int> _n;
  std::vector<aurostd::xvector<double> > _kpoints;
  std::vector<double> _weights;
  aurostd::xvector<double> _shift;
  aurostd::xmatrix<double> _rlattice;
  aurostd::xmatrix<double> _klattice;
  aurostd::xtensor3<int> _allToIrrPointMap;
  bool _isIrreducible;

 private:
  void generateAllGridPoints();
  void makeIrreducible(const std::vector<_sym_op>&);

 public:
  MonkhorstPackMesh(int, int, int, const xstructure&, Logger&);
  ~MonkhorstPackMesh();
  void clear();
  void setDensity(int);
  void setDensity(int, int, int);
  // Interface
  int getN(int);
  aurostd::xmatrix<double> getReciprocalLattice();
  std::vector<aurostd::xvector<double> > getPoints();
  std::vector<double> getWeights();
  aurostd::xtensor3<int> getAllToIrrPointMap();
};
}  // namespace apl

// ***************************************************************************
// "idosc.h"
namespace apl {
class IDOSCalculator {
 public:
  virtual std::vector<double> getBins() = 0;
  virtual std::vector<double> getDOS() = 0;
  virtual bool hasNegativeFrequencies() = 0;

 public:
  virtual ~IDOSCalculator() {}
};
}  // namespace apl

// ***************************************************************************
// "doscalc.h"
namespace apl {
#define MIN_FREQ_TRESHOLD -0.1  //in AMU
#define RAW2Hz 15.6333046177
#define AMU2Kg 1.66053904
#define MIN_EIGEN_TRESHOLD -1e-2  // eigenvalue treshold in AMU
class DOSCalculator : virtual public IDOSCalculator {
 protected:
  IPhononCalculator& _pc;
  IReciprocalPointGrid& _rg;
  Logger& _logger;
  std::vector<aurostd::xvector<double> > _qpoints;
  std::vector<double> _qweights;
  std::vector<aurostd::xvector<double> > _freqs;
  double _minFreq;
  double _maxFreq;
  double _stepDOS;
  double _halfStepDOS;
  std::vector<double> _bins;
  std::vector<double> _dos;
  //CO - START
 //private:
  void calculateInOneThread(int, int);
  //CO - END
 protected:
  void calculateFrequencies();
  void smearWithGaussian(vector<double>&, double, double);

 public:
  DOSCalculator(IPhononCalculator&, IReciprocalPointGrid&, Logger&);
  virtual ~DOSCalculator();
  //
  virtual void rawCalc(int) {}
  void calc(int);
  void calc(int, double);
  void clear();
  void writePDOS();
  void writePDOS(string, string);  //[PINKU]
  // Interface IDOSCalculator
  std::vector<double> getBins();
  std::vector<double> getDOS();
  bool hasNegativeFrequencies();
};
}  // namespace apl

// ***************************************************************************
// "rsmdos.h"
namespace apl {
class RootSamplingMethod : public DOSCalculator {
 public:
  RootSamplingMethod(IPhononCalculator&, IReciprocalPointGrid&, Logger&);
  ~RootSamplingMethod();
  void rawCalc(int);
};
}  // namespace apl

// ***************************************************************************
// "ltetdos.h"
namespace apl {
class LinearTetrahedronMethod : public DOSCalculator {
 private:
  double _weightVolumeOfEachTetrahedron;
  std::vector<std::vector<int> > _irrTetrahedraList;
  std::vector<int> _irrTetrahedraWeightList;

 private:
  void generateTetrahedras();

 public:
  LinearTetrahedronMethod(IPhononCalculator&, IReciprocalPointGrid&, Logger&);
  ~LinearTetrahedronMethod();
  void clear();
  void rawCalc(int);
};
}  // end namespace apl

// ***************************************************************************
// thermalpc.h
namespace apl {
enum ThermalPropertiesUnits { eV,
                              meV,
                              ueV,
                              eVK,
                              meVK,
                              ueVK,
                              kB };

class ThermalPropertiesCalculator {
 private:
  IDOSCalculator& _dosc;  //CO
  Logger& _logger;
  std::vector<double> _bins;
  std::vector<double> _dos;
  double _stepDOS;
  double _zeroPointVibrationEnergy_meV;
  bool _isCalcZeroPointVibrationEnergy_meV;

 private:
  double getScalingFactor(ThermalPropertiesUnits units);

 public:
  ThermalPropertiesCalculator(IDOSCalculator&, Logger&);
  ~ThermalPropertiesCalculator();
  void clear();
  void writeTHERMO(double, double, double);
  double getZeroPointVibrationEnergy(ThermalPropertiesUnits);
  double getInternalEnergy(double, ThermalPropertiesUnits);
  double getVibrationalFreeEnergy(double, ThermalPropertiesUnits);
  double getVibrationalEntropy(double, ThermalPropertiesUnits);
  double getIsochoricSpecificHeat(double, ThermalPropertiesUnits);
};

}  // namespace apl
   // ***************************************************************************
   // BEGIN JJPR: LATTICE THERMAL CONDUCTIVITY
   // ***************************************************************************
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include "aflow_apl.h"

//Constants
#define _kB 8.61734315E-05                 // Boltzmann Constant [eV/K]
#define kb_J 0.13806488E-22                // Boltzmann Constant [J/K]
#define NPI 3.14159265358979               // Pi
#define hbar 6.58211928E-16                // hbar [eVs]
#define hbarp 1.05457172647                // hbar Carrete JTHs* 1E-22
#define hbar_J 0.10545717E-21              // hbar Carrete JTHs* 1E-22
#define Thershold_freq 0.0001              //
#define FACTORR 0.98226977255434387350E14  //
#define ToPico 0.98226977255434387350E02

namespace apl {
class _QPs {  // List of q points
 public:
  _QPs();                                // allocate
  ~_QPs();                               // kill everything
  const _QPs& operator=(const _QPs& b);  // copy
  aurostd::xvector<double> q;
  aurostd::xvector<double> q_dir;
  aurostd::xvector<int> mapq;
  bool q_is_equivalent;
  int index_ineq;
  int index_fg;
  int index_fgc;

 private:
  void free();  // free space
};

class _iQPs {  // List of q points
 public:
  _iQPs();                                 // allocate
  ~_iQPs();                                // kill everything
  const _iQPs& operator=(const _iQPs& b);  // copy
  aurostd::xvector<double> q;
  aurostd::xvector<int> mapq;
  std::vector<uint> list;
  std::vector<uint> lsym;

 private:
  void free();  // free space
};

class ThermalConductivityCalculator {
 private:
  IPhononCalculator& _pc;
  Supercell& _sc;
  StrPairs& _pa;
  Logger& _logger;

  // Stored constants for scattering
  double Npc;
  uint nBranches;
  double alfaScat;

  // Group velocities
  vector<vector<xvector<double> > > groupVel;
  vector<vector<xvector<double> > > groupVel2;
  vector<xvector<double> > phon_vel;

  // Gamma centered grid
  double NPQs;
  xvector<int> grid;

  // Contributions to Q in MPI version
  std::vector<xvector<double> > qThreads;

  // Frequencies
  vector<vector<double> > freq;
  vector<xvector<double> > freqs;
  vector<double> phon_freq;

  //Thermal Conductivity
  vector<xvector<double> > ThermalCond;

  vector<vector<xvector<double> > > timeF;

  // USER OPTIONS
  IPCFreqFlags _frequencyFormat;
  bool _iso;

  // Eigenvectors and Derivative of the Dynamical Matrix
  std::vector<xmatrix<xcomplex<double> > > EigenVec;
  std::vector<xmatrix<xcomplex<double> > > _DMatrix;
  std::vector<vector<xmatrix<xcomplex<double> > > > dDynM;
  vector<xcomplex<double> > phon_eig;
  vector<xvector<double> > shortest;

  // Vectors for  phase factors
  vector<vector<xvector<double> > > tri_vec;

 public:
  ThermalConductivityCalculator(IPhononCalculator&, Supercell&, StrPairs&, Logger&);
  ~ThermalConductivityCalculator();

  // QPOINTS
  std::vector<_QPs> QPs;
  std::vector<_iQPs> iQPs;
  vector<vector<vector<int> > > mapQ;
  void qgrid(int, int, int);
  int get3rdQ_P(uint, uint);
  int get3rdQ_M(uint, uint);
  aurostd::xmatrix<double> _klattice;
  aurostd::xmatrix<double> _rlattice;
  xstructure xstr;

  // Rl
  vector<xvector<double> > Rl_list;
  void Rl();

  //Scaterring and elastic processes
  vector<vector<double> > ScatPlus;
  vector<vector<double> > ScatPlus_s;
  vector<vector<xvector<int> > > ScatPlus_q;
  vector<vector<xvector<int> > > ScatPlus_b;
  vector<vector<double> > ScatMinus;
  vector<vector<double> > ScatMinus_s;
  vector<vector<xvector<int> > > ScatMinus_q;
  vector<vector<xvector<int> > > ScatMinus_b;
  vector<vector<vector<vector<double> > > > IsoM;
  vector<vector<double> > Wq;
  vector<vector<double> > Wtotal;
  vector<vector<double> > Wiso;
  vector<vector<double> > tau_RTA;
  vector<vector<xvector<double> > > F_n;
  vector<vector<xvector<double> > > DeltaF;

  void clear();

  // Frequency function
  void calcFreq(int, int, int);
  void calculateFreqMatMPI(int, int);

  // Group velocity function
  void getdDynM(uint);
  void getGroupVelocity();
  void Print_GV();

  // Scaterring processes
  void getScattering();
  void calculateScat_MPI(int, int);
  void calculateScatProc_MPI(int);
  double getSplus(uint, uint, uint, uint, uint, uint, double);
  double getSminus(uint, uint, uint, uint, uint, uint, double);
  double getSigma(uint, uint, uint, uint);
  double getSigma(uint, uint);
  double getDeltaplus(uint, uint, uint, uint, uint, uint, double);
  double getDeltaminus(uint, uint, uint, uint, uint, uint, double);
  double getSiso(uint, uint, uint, uint, double);
  double getMassVariance(uint n);
  double getSigmaE(uint, uint);
  double getDeltaE(uint, uint, uint, uint, double);
  double getScatMatrixP(int, int, int, int, int, int);
  double getScatMatrixM(int, int, int, int, int, int);
  double Vol;

  void Calculator(bool, bool, bool, bool, double, double, double, double);
  void iteration0();
  void iteration(double);
  //double getTao_RTA( uint, uint, double);
  void getTao_RTA(double, bool, bool, double);
  xmatrix<double> getKappa(double);
  void getKappaMPI(int, int, int, double);
  xmatrix<double> getKappaQ(uint, uint, double);
  double getMFP(uint, uint);
  void getCumKappa(double);
  void getCumKappaMPI(int, int, int, double, double);
  vector<xmatrix<double> > SUM;
  void plotKappa();
  void calculateRTA_MPI(int, int, double);
  double getSpecificHeat(uint, uint, double);
  double getBoseDist(uint, uint, double);
  double getQ(uint, uint, uint, double);

  //Write TCOND file
  void writeThermalCond(bool, bool, double, double, double);

  // Scattering processes taking intoo account the temperature
  //double getWplus(uint, uint,uint, uint,uint, uint, double);
  //double getWminus(uint, uint,uint, uint,uint, uint, double);
  void getWiso();

  double getEigenProduct2(uint, uint, uint, uint, uint);
  double getMassProd(int, int, int);
  xcomplex<double> geteiRq(int, int);
  // PHONOPY
  template <typename T>
  std::vector<T> split(const std::string& line);
  void phonopy();
};
}
// ***************************************************************************
// END JJPR: LATTICE THERMAL CONDUCTIVITY
// ***************************************************************************

// ***************************************************************************
#include <iterator>
//Functions in this class are used to calculate Gruneisen Paramater related properties//
namespace apl {
class QuasiHarmonicGruneisen : public MVops {
 private:
  vector<string> ATOMIC_SPECIES;
  vector<double> ATOMIC_MASSES_AMU;
  vector<xmatrix<xcomplex<double> > > DMp;  //at +ve volume
  vector<xmatrix<xcomplex<double> > > DMm;  //at -ve volume
  vector<xmatrix<xcomplex<double> > > DM;
  vector<xmatrix<xcomplex<double> > > _eigenvectors;  //eigenvectors at equilibrium volume
  uint nBranches;
  bool _negative_freq_along_hsp;
  std::vector<aurostd::xvector<double> > _kpoints;
  vector<bool> gp_path_test;
  vector<bool> gp_mesh_test;
  aurostd::xmatrix<double> _rlattice;
  aurostd::xvector<int> _meshsize;
  double apl_inner_product(const vector<double> &a, const vector<double> &b);

 protected:
  IPhononCalculator& _pc;
  eos_calculate& _runeos;
  Logger& _logger;

  string dirfrefix;
  vector<string> gpdir;
  vector<string> phdir;
  vector<string> eosdir;
  vector<double> gpvol;
  vector<double> eosvol;
  std::vector<double> _weights;
  vector<xvector<double> > _GP_mesh;     // gruneisen parameter in qmesh
  vector<xvector<double> > GP_path;      // gruneisen parameter along path
  vector<xvector<double> > _freqs_path;  //frequecies along high symmetry path in THz
  vector<xvector<double> > _freqs_mesh;  //fequecies in qmesh in THz
  double delta_V;
  double V0;
  template <typename T>
  std::vector<T> split(const std::string& line);
  int ZERO_STATIC_DIR_INDEX;
  double CUTOFF_FREQ;

  //acoustic properties related variables
  std::vector<aurostd::xvector<int> > _BrINDEXs;
  std::vector<aurostd::xvector<double> > _acoustic_freqs_mesh;
  std::vector<aurostd::xvector<double> > _acoustic_GP_mesh;

 public:
  QuasiHarmonicGruneisen(IPhononCalculator&, eos_calculate&, Logger&);
  ~QuasiHarmonicGruneisen();
  void clear();
  bool set_directories(string dir);
  void set_frequency_cutoff(double d) { CUTOFF_FREQ = d; }
  void createmesh(int na, int nb, int nc, const xstructure&);
  bool cal_gp_along_path(const std::vector<aurostd::xvector<double> >& kpoints);
  bool cal_gp_in_mesh();
  void write_GP_path(const vector<double>& path, const vector<int>& path_seg);
  void write_GP_mesh();
  void WriteaverageGP(double, double, double);
  double averageGP(double);
  double AcousticAverageGP(double);
  void thermal_displacements(double Ts, double Te, int Tinc);
  virtual bool setvariables() { return false; }
  _CMAT_ xmat2mat(const xmatrix<xcomplex<double> >& M);
  _VEC_ xvec2vec(const xvector<double>& V);
  _MAT_ xmatd2matd(const xmatrix<double>& M);
  void projected_displacement(const _VEC_& direction, double Ts, double Te, int Tinc);
  void writeDM();
  void write_GP_FREQ_mesh();
  void calculate_acoustic_freqNgp();
  void Write_BrINDEXs();
  int getN(int i) { return _meshsize[i]; }

 private:
  bool read_PDIS(vector<string>& hash_lines);
  bool get_dynamicalmatrices_in_mesh();
  bool get_dynamicalmatrices_along_path();
  bool read_matrix(vector<xmatrix<xcomplex<double> > >& A, string file);
  void calculate_gp_in_path(int startIndex, int endIndex, int cpuid);
  void calculate_gp_in_mesh(int startIndex, int endIndex, int cpuid);
  template <typename T>
  string NumberToString(T Number);
  double get_Q2(double freq, double t);
  double get_population(double freq, double t);
  void identify_acoustic_modes();
  void trace_acoustic_modes_threads(int startIndex, int endIndex);
};
}
// ***************************************************************************
namespace apl {
//Functions in this class are used to calculate full equation of states//
class EOS : public QuasiHarmonicGruneisen {
 public:
  EOS(IPhononCalculator&, eos_calculate&, Logger&);
  ~EOS();
  void clear();
  bool setvariables();
  //calculate EOS with user defined temperatures
  void calEOS(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP, ThermalPropertiesCalculator& e);
  void set_fitting_type(string s) { FITTING_TYPE = s; }

 private:
  //error managments
  bool data_read_error;
  int nl_success_status;
  string nl_err_msg;
  double CUTOFF_FREQ;
  //lfit output
  double luncertanity_V0;  //uncertainties in V0 from linear least squares fit
  double luncertanity_E0;  //uncertainties in E0 from linear least squares fit
  double luncertanity_B0;  //uncertainties in B0 from linear least squares fit
  double lchisq;           //chisquare from linear least squares fit
  double leqmV0;           //Equilibrium V0 from linear least squares fit
  double leqmE0;           //Equilibrium E0 from linear least squares fit
  double leqmB0;           //Equilibrium B0 from linear least squares fit
  //nlfit output
  double uncertanity_V0;   //uncertainties in V0 from nonlinear least squares fit
  double uncertanity_E0;   //uncertainties in E0 from nonlinear least squares fit
  double uncertanity_B0;   //uncertainties in B0 from nonlinear least squares fit
  double uncertanity_Bp;   //uncertainties in Bp from nonlinear least squares fit
  double uncertanity_Bpp;  //uncertainties in Bpp from nonlinear least squares fit
  double chisq_dof;        //chi square per degrees of freedom
  double nleqmV0;          //Equilibrium V0 from linear least squares fit
  double nleqmE0;          //Equilibrium E0 from linear least squares fit
  double nleqmB0;          //Equilibrium B0 from linear least squares fit
  double nleqmBp;          //Equilibrium Bp from linear least squares fit
  double nleqmBpp;         //Equilibrium Bpp from linear least squares fit
  string fdfsolver_name;
  string FITTING_TYPE;
  double Feqm, Beqm, Veqm, Bp, Bpp;  //equilibrium properties

  vector<string> phdir;
  vector<string> eosdir;
  string dirfrefix;
  vector<double> E0K;                     //0K energies
  vector<double> pV;                      //pV energies
  vector<double> MagCell;                 //Magnetization per cell
  vector<double> zpe;                     //zero point energy
  vector<double> fermi;                   //fermi energy
  vector<double> eosvol;                  //eos volumes
  bool _ismagnetic;                       //check magnetic material
  vector<vector<vector<double> > > pDOS;  //phonon DOS
  vector<vector<vector<double> > > eDOS;  //electronic DOS
  int ZERO_STATIC_DIR_INDEX;
  vector<uint> NEGATIVE_FREQ;
  double chisq_quadratic;
  double chisq_birch_murnaghan;
  uint Iteration_birch_murnaghan;
  double alamda_birch_murnaghan;
  vector<double> Uncertainties_birch_murnaghan;
  void calculateE0K();
  bool calculatePDOS();
  double getE0K(string file, double& pv, double& magcell);
  template <class T>
  vector<vector<T> >
  readfile(string file, const vector<uint>& readCOL);
  void getZeroPointVibrationEnergy();
  void geteDOS();
  template <class T>
  vector<vector<T> > geteDOS(T& fermi, string file);
  double VibrationEnergy(double temperature_in_kelvins, uint dir_index);
  double ElectronicEnergy(double temperature_in_kelvins, uint dir_index);
  double fermi_dirac_distribution(double Ei, double fermi_energy, double temperature_in_kelvins);
  bool EVfit(const xvector<double>& E, const xvector<double>& V);
  double getIsochoricSpecificHeat(double temperature_in_kelvins, uint dir_index);
  void resize_EOSnPH();
  void md_lsquares_call(const xvector<double>& xdata, const xvector<double>& ydata);
  void initialize_output_variables();
  bool more_refinement(const xvector<double>& E, const xvector<double>& V, xvector<double>& guess, xvector<double>& out);

 public:
  template <typename T>
  std::vector<T> split(const std::string& line);
};
}
// ***************************************************************************
namespace apl {
//this class store dynamical matrices and phonon dispersion files from other directories
class DM_PDOS_save {
 private:
  IPhononCalculator& _pc;
  Logger& _logger;
  string _EOSdir;
  string dirfrefix;
  uint nBranches;
  double GP_VOL_DISTORTION;
  std::vector<aurostd::xvector<double> > _kpoints;
  vector<xmatrix<xcomplex<double> > > DM;
  std::vector<double> _weights;
  aurostd::xmatrix<double> _rlattice;
  aurostd::xmatrix<double> _klattice;

 public:
  DM_PDOS_save(IPhononCalculator&, Logger&);
  ~DM_PDOS_save();
  void clear();
  void createMPmesh(int, int, int, const xstructure&);
  void create_pdispath(std::vector<xvector<double> >& qpoints);
  void setdir_prefix(std::string s, double d) {
    dirfrefix = s;
    GP_VOL_DISTORTION = d;
  }
  void set_GP_VOL_DISTORTION(double d) { GP_VOL_DISTORTION = d; }
  bool check_GP();
  void PDOSsave();
  void DMsave();
  string
  getdir_name(string path);
  std::vector<aurostd::xvector<double> >
  get_kpoints() { return _kpoints; }
  std::vector<double>
  get_weights() { return _weights; }
  vector<string>
  directory_list(const string path);

 private:
  template <typename T>
  string NumberToString(T Number);
  void
  get_dm(int startIndex, int endIndex);
  void write_dm(string file);
  void write_kpoints(string file);
  void write_weight(string file);
};
}
// ***************************************************************************
namespace apl {
//Functions in this class check the inconsistencies of GRUNEISEN subdirectories
//and make changes according to the main aflow.in
class check_consistency_aflow_apl {
 private:
  string _RUNNING_DIR;
  Logger& _logger;

 public:
  check_consistency_aflow_apl(Logger& l) : _logger(l) {}
  ~check_consistency_aflow_apl() {}
  void clear() { } //this->clear(); }
  void getdir_name(string s) { _RUNNING_DIR = s; }
  void check_consistency_in_aflow(string& current_running_aflowin);
  bool replace(std::string& str, const std::string& from, const std::string& to);
};
}
// ***************************************************************************
// shellhandle.h

namespace apl {

struct ShellData {
  int occupation;
  int occupationCapacity;
  bool isFull;
  double radius;
  double stdevRadius;
  vector<xvector<int> > index;
  vector<deque<_atom> > atoms;
  vector<deque<_atom> > ratoms;

  ~ShellData();
  ShellData& operator=(const ShellData&);
};

class ShellHandle {
 private:
  int _idSafeGeneratedShell;
  int _idSafeMappedShell;

  int _centralAtomID;
  double _indexReductionConstant;
  xstructure _initStructure;
  xstructure _initStructure_original;  //corey, does not include HEAVY symmetry stuff
  //deque<_atom> _initStructure_atoms_original; //corey

  xvector<int> _safeDimension;
  vector<ShellData> _shells;

 private:
  xvector<double> getFPositionItsNearestImage(const xvector<double>&,
                                              const xvector<double>&,
                                              const xmatrix<double>&);

 public:
  ShellHandle();
  ShellHandle(const xstructure&, int, int);
  ~ShellHandle();

  void clear();
  void init(const xstructure&, int, int);

  double getShellRadius(int);
  int getShell(double);

  double getSafeShellRadius();
  void setSafeShell(int);
  int getSafeShell();

  void calcShells(const xstructure&, int, int);
  void splitBySymmetry();
  void removeSplitBySymmetry();
  void addAtomToShell(int, const _atom&, bool = true);
  void mapStructure(const xstructure&, int, bool = true);

  int getLastOccupiedShell();
  int getLastRegularShell();
  int getLastFullShell();

  int getNumberOfShells();
  int getNumberOfSubshells(int);
  std::deque<_atom> getAtomsAtSameShell(int, int = 0);
  const std::deque<_atom>& getReferenceAtomsAtSameShell(int, int = 0);

  double getIndexReductionConstant() { return _indexReductionConstant; }

  void printReport(ostream&);

  //COREY added here
  void center(int);
  void center_original(void);
  //COREY added here
};

}  // end namespace apl

// ***************************************************************************
// aplmath.h

namespace aurostd {
// Calculate vector projection b on a
template <class utype>
xvector<utype> getVectorProjection(const xvector<utype>& b, const xvector<utype>& a) {
  return (a * (utype)(scalar_product(a, b) / scalar_product(a, a)));
}
// Calculate vector projection c on a like b on a
template <class utype>
xvector<utype> getModeratedVectorProjection(const xvector<utype>& c, const xvector<utype>& b, const xvector<utype>& a) {
  return (c * (utype)(scalar_product(a, b) / scalar_product(a, a)));
}
// Calculate vector convolution
template <class utype>
xvector<utype> getVectorConvolution(const xvector<utype>& a, const xvector<utype>& b) {
  xvector<utype> v(a.rows);
  for (int i = v.lrows; i <= v.urows; i++) {
    v(i) = a(i) * b(i);
  }
  return (v);
}
}

namespace apl {
void tred2(xmatrix<xcomplex<double> >&);
void zheevByJacobiRotation(xmatrix<xcomplex<double> >&, xvector<double>&, xmatrix<xcomplex<double> >&);
#ifdef USE_MKL
void zheevMKL(xmatrix<xcomplex<double> >&, xvector<double>&, xmatrix<xcomplex<double> >&);
#endif
}

// ***************************************************************************
// cursor.h

#define cursor_moveyx(y, x) printf("\033[%d;%dH", y, x) /*Move cursor to position y,x (rows, columns) with (1,1) as origin*/
#define cursor_moveup(y) printf("\033[%dA", y)          /*Move cursor up y*/
#define cursor_movedown(y) printf("\033[%dB", y)        /*Move cursor down y*/
#define cursor_moveright(x) printf("\033[%dC", x)       /*Move cursor right x*/
#define cursor_moveleft(x) printf("\033[%dD", x)        /*Move cursor left x*/
#define cursor_store() printf("\033[s")                 /*Store current cursor position and color*/
#define cursor_restore() printf("\033[u")               /*Restore cursor position and color from cursor_store()*/
#define cursor_clear() printf("\033[2J")                /*Clear screen and leave cursor where is*/
#define cursor_clearline() printf("\033[K")             /*Clear to end of line and leave cursor where is*/
#define cursor_fore_black() printf("\033[30m")          /*Change foreground color to black*/
#define cursor_fore_red() printf("\033[31m")            /*Change foreground color to red*/
#define cursor_fore_green() printf("\033[32m")          /*Change foreground color to green*/
#define cursor_fore_orange() printf("\033[33m")         /*Change foreground color to orange*/
#define cursor_fore_blue() printf("\033[34m")           /*Change foreground color to blue*/
#define cursor_fore_magenta() printf("\033[35m")        /*Change foreground color to magenta*/
#define cursor_fore_cyan() printf("\033[36m")           /*Change foreground color to cyan*/
#define cursor_fore_yellow() printf("\033[33m\033[1m")
#define cursor_fore_white() printf("\033[37m")    /*Change foreground color to white*/
#define cursor_back_black() printf("\033[40m")    /*Change background color to black*/
#define cursor_back_red() printf("\033[41m")      /*Change background color to red*/
#define cursor_back_green() printf("\033[42m")    /*Change background color to green*/
#define cursor_back_orange() printf("\033[43m")   /*Change background color to orange*/
#define cursor_back_blue() printf("\033[44m")     /*Change background color to blue*/
#define cursor_back_magenta() printf("\033[45m")  /*Change background color to magenta*/
#define cursor_back_cyan() printf("\033[46m")     /*Change background color to cyan*/
#define cursor_back_white() printf("\033[47m")    /*Change background color to white*/
#define cursor_attr_none() printf("\033[0m")      /*Turn off all cursor attributes*/
#define cursor_attr_bold() printf("\033[1m")      /*Make test bold*/
#define cursor_attr_underline() printf("\033[4m") /*Underline text*/
#define cursor_attr_blink() printf("\033[5m")     /*Supposed to make text blink, usually bolds it instead*/
#define cursor_attr_reverse() printf("\033[7m")   /*Swap background and foreground colors*/

// ***************************************************************************
// xtensor.hpp

namespace apl {

//////////////////////////////////////////////////////////////////////////////
template <typename T, unsigned int TENSOR_ORDER>
class xtensor {
  T* _data;
  int _lindex;
  int _hindex;
  unsigned int _indexSize;
  unsigned long long _arraySize;
  std::vector<unsigned long long> _precomputedOffsets;

 private:
  unsigned long long getOffset(const std::vector<unsigned int>&);

 public:
  xtensor(int, int);
  ~xtensor();
  void zero();
  void fill(const T&);
  T& operator()(int, ...);
};

//////////////////////////////////////////////////////////////////////////////
template <typename T, unsigned int TENSOR_ORDER>
xtensor<T, TENSOR_ORDER>::xtensor(int index1, int index2) {
  _hindex = index1 > index2 ? index1 : index2;
  _lindex = index1 < index2 ? index1 : index2;

  _indexSize = _hindex - _lindex + 1;

  // Allocate data
  unsigned long long arraySizeBefore;
  _arraySize = 1ULL;
  _precomputedOffsets.push_back(_arraySize);
  for (_AFLOW_APL_REGISTER_ unsigned int i = 0; i < TENSOR_ORDER; i++) {
    arraySizeBefore = _arraySize;
    _arraySize *= _indexSize;
    // Detect overflow
    if (arraySizeBefore > _arraySize)
      throw APLRuntimeError("apl::xtensor<T>::xtensor(); The setting is producing an array which can not by handled by this implementation.");
    _precomputedOffsets.push_back(_arraySize);
  }

  try {
    // FIX: Problem; new(std::size_t), where size_t is hardware/platform specific,
    // for 32 bit systems it is uint = 2^32, if our array is bigger than this,
    // there is a overflow -> use more pages of the same block _data[PAGE][2^32]???,
    // or go to the 64 bits systems -> 2^64 (unsigned long long)
    if (_arraySize > (1ULL << (8 * sizeof(std::size_t) - 1)))  // dont go overboard  STEFANO (-1)
      throw APLRuntimeError("apl::xtensor<T>::xtensor(); Problem to allocate a required array. Hardware specific problem.");
    _data = new T[(std::size_t)_arraySize];
  } catch (std::bad_alloc& e) {
    throw APLRuntimeError("apl::xtensor<T>::xtensor(); Bad allocation.");
  }

  // Reverse precomputed for better manipulation and shift by 1
  std::vector<unsigned long long> temp(_precomputedOffsets.rbegin() + 1, _precomputedOffsets.rend());
  _precomputedOffsets = temp;
  temp.clear();
}

//////////////////////////////////////////////////////////////////////////////
template <typename T, unsigned int TENSOR_ORDER>
xtensor<T, TENSOR_ORDER>::~xtensor() {
  _precomputedOffsets.clear();
  delete[] _data;
  _data = NULL;
  _arraySize = 0ULL;
}

//////////////////////////////////////////////////////////////////////////////
template <typename T, unsigned int TENSOR_ORDER>
void xtensor<T, TENSOR_ORDER>::zero() {
  // Take care! The last argument is of std::size_t, hence platform specific
  memset(_data, 0, _arraySize * sizeof(T));
}

//////////////////////////////////////////////////////////////////////////////
template <typename T, unsigned int TENSOR_ORDER>
void xtensor<T, TENSOR_ORDER>::fill(const T& val) {
  for (unsigned long long i = 0ULL; i < _arraySize; i++)
    _data[i] = val;
}

//////////////////////////////////////////////////////////////////////////////
template <typename T, unsigned int TENSOR_ORDER>
unsigned long long xtensor<T, TENSOR_ORDER>::getOffset(const std::vector<unsigned int>& idxs) {
  unsigned long long offset = idxs.back();
  for (_AFLOW_APL_REGISTER_ int i = idxs.size() - 2; i >= 0; i--)
    offset += (unsigned long long)(idxs[i] * _precomputedOffsets[i]);
  return offset;
}

//////////////////////////////////////////////////////////////////////////////
template <typename T, unsigned int TENSOR_ORDER>
T& xtensor<T, TENSOR_ORDER>::operator()(int idx1, ...) {
  // Get indices and transform to form 0...max
  va_list arguments;
  std::vector<unsigned int> idxs;

  // The 1st index
  idxs.push_back((unsigned int)(idx1 - _lindex));
  //    if( idxs.back() < 0 || idxs.back() > _indexSize ) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");
  if (idxs.back() > _indexSize) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");  // unsigned long long cant be negative

  // The rest of indices
  va_start(arguments, idx1);
  for (_AFLOW_APL_REGISTER_ unsigned int i = 0; i < TENSOR_ORDER - 1; i++) {
    idxs.push_back((unsigned int)(va_arg(arguments, int) - _lindex));
    //   if( idxs.back() < 0 || idxs.back() > _indexSize ) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");
    if (idxs.back() > _indexSize) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");  // unsigned long long cant be negative
  }
  va_end(arguments);

  // Calculate offset
  unsigned long long offset = getOffset(idxs);
  idxs.clear();

  // Return value
  return _data[offset];
}
//////////////////////////////////////////////////////////////////////////////
}
// ***************************************************************************

#endif  // _AFLOW_APL_H_

// ***************************************************************************
