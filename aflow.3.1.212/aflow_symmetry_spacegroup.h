// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2016           *
// *                                                                         *
// ***************************************************************************
// Written by Richard H. Taylor
// UPDATED by David Hicks
// d.hicks@duke.edu

#ifndef _AFLOW_SYMMETRY_SPACEGROUP_H_
#define _AFLOW_SYMMETRY_SPACEGROUP_H_
#include <cstring>
#include <map>
#include "aflow.h"
const double Pi_r = 3.141592653589793;

// ******************************************************************************
// Tolerance Functions
// ******************************************************************************
namespace SYM{
  void SetTolerance(const double& starting_tol);
}

// ******************************************************************************
// Class Declarations
// ******************************************************************************
//xstructure
class xstructure;  //forward define xstructure for compilation.

vector<int> AllCombination41(int num, int total_num, int index);
vector<int> AllCombination42(int num, int total_num, vector<int>& str_in);
unsigned long int CombinationNr(int num, int total_num);

namespace SYM {
  //stringdouble (Used for Wyckoff position algebra)
  class stringdouble {
    friend ostream& operator<<(ostream& output, const stringdouble& a);

   public:
    string type;
    double distance;
    xvector<double> coord;
    stringdouble() {
      type = "XX";
      distance = 0;
    }
    ~stringdouble(){};
  };

  //wyckoff site
  class wyckoffsite  //Also for wyckoff sites
  {
    //friend bool operator==(const wyckoffsite& a, const wyckoffsite& b);
    friend ostream& operator<<(ostream& output, const wyckoffsite& a);

   public:
    wyckoffsite() {
      wyckoffSymbol = " ";
    };
    ~wyckoffsite(){};
    xvector<double> coord;
    string type;  //chemical label etc
    string wyckoffSymbol;
  };

  //atom class manipulation function
  deque<int> arrange_atoms(deque<_atom>& atoms);

  //Screw (Screw Operation)
  class Screw {
    friend xvector<double> operator*(Screw& S, const xvector<double>& P);
    friend xvector<double> operator*(const xvector<double>& P, Screw& S);
    friend ostream& operator<<(ostream& output, const Screw& S);

   private:
    xmatrix<double> A;  //operator matrix
    void get_A();
    string linestring;          //parametric line
    xvector<double> T;          //translation component of screw
    double angle;               //angle of rotation
    xvector<double> one_point;  //a point on the axis (in get_A)
    xvector<double> direction_vector;

   public:
    Screw(){};
    ~Screw(){};
    void get_screw_direct(xvector<double> axis_dir, xvector<double> point, double order);
    void get_line_string(string str);
    void get_trans_comp(xvector<double> trans);
    void get_angle(double theta);
    xvector<double> return_direction();
    xvector<double> return_point();
    double return_order();
  };

  //Glide (Glide Operation)
  class Glide {  //for a plane defined on the plane ax+by+cz=d with translation
		 //(t1,t2,t3) the glide operation works as follows:
		 // A.P + d*a + |a.P - d|*a + (t1,t2,t3)^T
   private:
    bool HEX;
    bool DIRECT;
    string planestring;           //parametric plane
    xmatrix<double> A;            //reflecting matrix
    xvector<double> a;            //normal vector or fixed point (if HEX)
    xvector<double> T;            //translation vector
    xvector<double> plane_point;  //point in plane for get_glide_direct
    double d;                     // d in ax+by+cz+d
    void get_A(xvector<double> n);

   public:
    Glide() {
      HEX = false;
      DIRECT = false;
      d = 0.0;
    };
    ~Glide(){};
    //overload operator * for Glide
    // cartesian
    void get_glide_direct(xvector<double> n, xvector<double> p);
    void get_glide_direct(xvector<double> n, xvector<double> p, xvector<double> trans);
    xvector<double> return_point();      //returns point on plane
    xvector<double> return_direction();  //returns normal to plane
    friend xvector<double> operator*(Glide& G, const xvector<double>& P);
    friend xvector<double> operator*(const xvector<double>& P, Glide& G);
    friend ostream& operator<<(ostream& output, const Glide& S);
  };

  //Translation
  class Translation {
   private:
    xvector<double> translation_vec;

   public:
    Translation(){};
    ~Translation(){};
    void get_translation(string ITCstring);
    friend xvector<double> operator*(Translation& T, const xvector<double>& P);
    friend xvector<double> operator*(const xvector<double>& P, Translation& T);
    friend xvector<double> operator+(Translation& T, const xvector<double>& P);
    friend xvector<double> operator+(const xvector<double>& P, Translation& T);
  };

  //Inversion
  class Inversion {
    friend xvector<double> operator*(Inversion& I, const xvector<double>& P);
    friend xvector<double> operator*(const xvector<double>& P, Inversion& I);

   private:
   public:
    Inversion(){};
    ~Inversion(){};
    //overload operator * for Glide
    xvector<double> inversion_point;
    void get_inversion(string ITCstring);
  };
} //namespace SYM

//End Classes
// ******************************************************************************

// ******************************************************************************
// Structure Declarations
// ******************************************************************************

namespace SYM {
  //Symfolder (Stores all symmetry elements/operators)
  struct symfolder {
    vector<Screw> twofold_ops;
    vector<Screw> rot_ops;
    vector<Glide> mirror_ops;
    vector<xvector<double> > twofold_lattice_vectors;
    vector<xvector<double> > rot_lattice_vectors;
    vector<xvector<double> > mirror_lattice_vectors;
    bool commensurate;
    string crystalsystem;
    string latticesystem;
    vector<int> lattice_pgroups;                // DX NEW
    vector<xmatrix<double> > lattice_sym_mats;  // DX NEW
    vector<xmatrix<double> > crystal_sym_mats;  // DX NEW
    vector<int> insym;
    vector<xvector<double> > translations;
    vector<vector<int> > all_atom_maps;
    vector<vector<int> > all_type_maps;
    vector<xvector<double> > centeringops;
    vector<string> pointgroupops;
  };

  //Symop
  struct symop {
    friend ostream& operator<<(ostream& output, const symop& a);
    string symbol;
    xvector<double> direction;
    xvector<double> screwglide;  //will also store inversion pnts for rotoinversion
    xvector<double> shift;

    void clear();
  };

  //Eqatoms
  struct eqatoms {
    string label;
    vector<vector<double> > atoms;
  };

  //sdouble (String and double)
  struct sdouble {
    char chr;
    double dbl;
    sdouble() {
      chr = '\0';
      dbl = 0;
    }
  };

  //enum_alpha (enumerate alphabet: for Wyckoff letters)
  struct enum_alpha {
    string letter;
    int number;
  };

  //End Structures
  // ******************************************************************************

  // ******************************************************************************
  // Template Declarations
  // ******************************************************************************

  template <class d>
  d infnorm(vector<d> vec) {
    if (vec.size() == 0) {
      cerr << "Empty vector!:infnorm" << endl;
      exit(1);
    }
    d largest = vec[0];
    for (uint i = 0; i < vec.size(); i++) {
      if (vec[i] > largest) {
	largest = vec[i];
      }
    }
    return largest;
  }

  template <class d>
  int infnormindex(vector<d> vec) {
    if (vec.size() == 0) {
      cerr << "Empty vector!:infnormindex" << endl;
      exit(1);
    }
    d largest = vec[0];
    int out;
    for (uint i = 0; i < vec.size(); i++) {
      if (vec[i] > largest) {
	//largest = vec[i];
	out = i;
      }
    }
    return out;
  }

  template <class d>
  d infnorm(xvector<d> xvec) {
    if (xvec.urows == 0) {
      cerr << "Empty vector!:infnorm xvec" << endl;
      exit(1);
    }
    d largest = aurostd::abs(xvec[1]);
    for (int i = 0; i < xvec.urows; i++) {
      if (aurostd::abs(xvec[i + 1]) > largest) {
	largest = aurostd::abs(xvec[i + 1]);
      }
    }
    return largest;
  }

  template <typename T>
  void remove_duplicates_rt(vector<T>& vec) {
    //sort(vec.begin(), vec.end());//sorts ascending order
    vec.erase(unique(vec.begin(), vec.end()), vec.end());
    //sort(vec.begin(), vec.end(), std::greater<int>());//sorts in descending order
  }

  template <class dmmy>
  bool allsame(vector<dmmy> vec) {
    bool all = true;
    for (uint i = 0; i < vec.size(); i++) {
      if (vec[i] != vec[0])
	all = false;
    }
    return all;
  }

  template <class dmmy>
  bool invec(vector<dmmy> vec, dmmy n) {
    bool contains = false;
    for (uint i = 0; i < vec.size(); i++) {
      if (vec[i] == n)
	contains = true;
    }
    return contains;
  }

  template <class dmmy>
  int invec_index(vector<dmmy> vec, dmmy n) {
    int out;
    for (uint i = 0; i < vec.size(); i++) {
      if (vec[i] == n)
	out = i;
    }
    return out;
  }

  int invecindex(vector<_atom> vec, _atom n);

  template <class dmmy>
  bool inveconce(vector<dmmy> vec, dmmy n) {
    bool contains = false;
    int count = 0;
    for (uint i = 0; i < vec.size(); i++) {
      if (vec[i] == n) {
	count++;
      }
    }
    //if(count >= 1){
    if (count == 1) {
      contains = true;
    }
    return contains;
  }
} //namespace SYM

//End Templates
// ******************************************************************************

// ******************************************************************************
// Function Declarations
// ******************************************************************************
//MAIN FUNCITONS
namespace SYM {
  string OrthoDefect(istream& cin);
  xstructure SpaceGroup(istream& cin);
  void rgcd(vector<string> num);
  // DX void AtomicEnvironment(istream & cin, vector<string> av);
  string ReverseSpaceGroup(vector<string> num);
  //END MAIN FUNCTIONS

  //FUNCTION THAT TAKES WYCCAR FROM XSTRUCTURE AND PUTS IT IN OSTREAM
  void printWyccar(ofstream& FileMESSAGE, xstructure& str, const bool& osswrite, ostream& oss);

  //LINEAR ALGEBRA FUNCTIONS
  bool solve_overdetermined_system(vector<xvector<double> >& LHS, vector<double>& RHS, xvector<double>& SOL, xmatrix<double>& lattice, double& min_dist);
  void ReducedRowEchelonForm(xmatrix<double>& M);
  bool find_solution_UD(xmatrix<double> M, xvector<double>& SOL);  //check if underdetermined system has a solution, and get a particular solution.
  bool checkLinearSystem(vector<xvector<double> >& LHS, vector<double>& RHS, xmatrix<double>& lattice);
  //WYCKOFF FUNCITONS
  vector<double> system_solve(vector<double> numeric, vector<vector<sdouble> > variable);
  vector<vector<double> > wyckoff_solve(vector<vector<vector<sdouble> > > win, vector<eqatoms> pin);
  vector<vector<double> > wyckoff_solve_2(vector<vector<vector<sdouble> > > win, eqatoms pin);
  vector<vector<double> > wyckoff_sites(vector<vector<vector<sdouble> > > win, vector<_atom> atomgroup);

  //SPACE GROUP LIBRARY FUNCTIONS
  bool initsgs(string axis_cell);
  bool initsymops();
  bool initsymmats();
  bool initglides();
  bool initgenerators(string axis_cell);
  vector<int> generatorindex(int spacegroup);
  vector<xvector<double> > ReturnITCGenShift(int spacegroupnum, string axis_cell);

  string point_group_library_search(vector<string> symmetryoperations, string crystalsystem, int centeringops);
  string crystal_system_library_search(vector<string> symmetryoperations);
  vector<int> spacegroups_CS_LS(string crystalsystem, char latticesystem);
  char discern_sym_op(string ITCstring);
  int* PointGroup_SpaceGroup(string pgroup);
  vector<int> PointGroup_SpaceGroup(string pgroup, char cl);
  xmatrix<double> spacegroup_to_lattice(int spacegroup);

  vector<int> get_multiplicities(string sg);
  vector<string> get_symmetry_symbols(string sg);

  vector<string> get_wyckoff_equation(string spaceg, int mult);  // DX 8/30/17
  vector<vector<vector<string> > > get_wyckoff_pos(string spaceg, int mult);
  vector<vector<vector<string> > > get_wyckoff_pos(string spaceg, int mult, bool getcentering);
  vector<string> get_minimum_enumerated_Wyckoff_letters(string spacegroupstring, vector<int>& multiplicities, vector<string> site_symmetries);
  int enumerate_wyckoff_letter(string& wyckoff_letter);
  vector<int> enumerate_wyckoff_letters(vector<string>& wyckoff_letters); //DX 20180927
  void get_all_wyckoff_for_site_symmetry(string spaceg, int mult, string site_symmetry, vector<vector<string> >& all_positions);
  xvector<double> Wyckoff_position_string2xvector(string& string_position);
  vector<xvector<double> > get_possible_origin_shifts(string spacegroupstring, int multiplicity, string site_symmetry);
  void get_certain_wyckoff_pos(string spaceg, int mult, string site_symmetry, vector<string>& site_symmetries, vector<string>& letters, vector<string>& positions);
  void getGeneralWyckoffMultiplicityAndPosition(uint space_group_number, string& space_group_setting, int& general_wyckoff_multiplicity, vector<string>& general_wyckoff_position);
  vector<string> findGeneralWyckoffPosition(string& spacegroupstring, int& general_wyckoff_multiplicity);
  bool shiftWyckoffPositions(deque<deque<_atom> >& equivalent_atoms_shifted, xvector<double>& previous_shift, xvector<double>& new_shift);
  bool findWyckoffPositions(xstructure& CCell, deque<_atom>& atomicbasis, vector<vector<vector<string> > >& tmpvvvstring,
			    deque<deque<_atom> >& equivalent_atoms, deque<deque<_atom> >& equivalent_atoms_shifted,
			    bool& foundspacegroup, string& spacegroupstring, bool& orig_origin_shift, xvector<double>& OriginShift,
			    vector<int>& wyckoffmult, vector<string>& wyckoffsymbols, vector<wyckoffsite_ITC>& wyckoffVariables,
			    deque<_atom>& wyckoffPositionsVector, vector<string>& wyckoffSymbols, ostringstream& woss,
			    bool& obverse_force_transformed);

  void print_wyckoff_pos(vector<vector<vector<string> > > wyckoff_positions);
  vector<vector<vector<vector<sdouble> > > > convert_wyckoff_pos_sd(vector<vector<vector<string> > > wyckoff_positions);
  void convert_wyckoff_pos(vector<vector<vector<string> > > wyckoff_positions);
  vector<vector<string> > get_centering(string spaceg);

  string ExtractWyckoffAttributesString(vector<string>& wyccar_ITC, uint attribute_index); //DX 201780823 
  string ExtractWyckoffLettersString(vector<string>& wyccar_ITC); //DX 201780823
  string ExtractWyckoffMultiplicitiesString(vector<string>& wyccar_ITC); //DX 201780823
  string ExtractWyckoffSiteSymmetriesString(vector<string>& wyccar_ITC); //DX 201780823

  //TOPOLOGY FUNCTIONS
  vector<xvector<double> > find_vectors_inplane(const vector<xvector<double> >& big_expanded, const xvector<double>& perp_to_vec);

  //check if three distances can define a screw translation
  //shortest vectors. "m" specifies how many you want (e.g., 2 --> the smallest 2)
  vector<int> shortest_vec(vector<xvector<double> > vecs, int m);
  //Closest point in terms of lattice:
  //[OBSOLETE] xvector<double> closest_point(xvector<double> norm, xmatrix<double> Linv, xmatrix<double> L);
  xmatrix<double> get_mod_angle_tensor(vector<xvector<double> >& points);  //Points is a vector of points (ordered pair, triplet, etc.);
  //Check if two points are equivalent under the lattice L
  bool points_equivalent(xmatrix<double>& c2f, xmatrix<double>& f2c, xvector<double> P1, xvector<double> P2, xvector<double>& lattice_vector, double& radius, bool& skew);
  //Check if the candidate operation is a symmetry of the crystal C
  bool symmetry_axis_equivalent(xmatrix<double> L, xmatrix<double> Linv, xvector<double> P1, xvector<double> N1, xvector<double> P2, xvector<double> N2);
  bool screw_equivalent(Screw S1, Screw S2);
  bool mirror_plane_equivalent(xvector<double> normal);
  xvector<double> next_point_on_line(xvector<double> P, xmatrix<double> L);
  void add_3d_point(vector<xvector<double> >& points, double x, double y, double z);
  char discern_rot_sym(xmatrix<double> m);
  bool allsame(vector<double> v);
  bool is_lattice_point(xmatrix<double> L, xvector<double> point, xvector<double>& lattice_vector, double& radius, bool& skew);
  bool in_cell(xmatrix<double> Linv, xvector<double> point);
  bool in_cell(xvector<double> P);                     //For use
  bool in_cell(xvector<double> P, xvector<double> N);  //Check if in cell modulo the direction of N (normal)
  double distance_between_points(const xvector<double>& a, const xvector<double>& b);

  double get_angle(xvector<double> a, xvector<double> b, string c);
  xvector<double> get_random_point_in_plane(string plane);
  xvector<double> get_point_in_line(string line, double param);
  xvector<double> random_point();

  //STRING FUNCTIONS
  bool blank(string str_in);
  bool containschar(string str_in);
  bool iselem(string str_in);
  bool havechar(string str_in, char c);  //see template function "invec"
  bool intinvec(vector<int> vint, int n);
  int whereischar(string str, char c);
  char whichchar(string str_in);
  double whichnum(string str_in);
  double frac2dbl(string str);  //expand to cover case when input is e.g., ".5"
  string dbl2frac(double a, bool sign_prefix=true);
  void multiply(vector<string> A, vector<string> B);
  void xstring(ostream& output, xmatrix<double> a);
  void cleanupstring(string& str);  //eliminates blank spaces before and after string
  vector<string> splitstring(string str, char c);  //c is delimeter
  vector<string> splitstringbyspaces(string str);
  vector<string> splitstring(string str);
  vector<string> splitstring_2(string str);
  vector<sdouble> simplify(string str);
  xvector<double> get_triplet(string str);

  // Used to check atomic basis
  bool GCD_conventional_atomic_basis(deque<_atom>& conventional_basis_atoms, deque<deque<_atom> >& prim_split_atom_types, int& prim_GCD);
  //[CO 180409 - moved to xatom]deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew);
  //[CO 180409 - moved to xatom]deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew, double& tol);
  bool MapAtomsInNewCell(_atom& a, _atom& b, xmatrix<double>& c2f_orig, xmatrix<double>& f2c_new, bool& skew, double& tol);
  bool MapAtomsInNewCell(xvector<double>& a, xvector<double>& b, xmatrix<double>& c2f_orig, xmatrix<double>& f2c_new, bool& skew, double& tol);
  deque<deque<_atom> > groupSymmetryEquivalentAtoms(deque<_atom>& atoms, xmatrix<double>& lattice, vector<xmatrix<double> >& sym_ops,
						    vector<xvector<double> >& translations, double& min_dist);
  deque<deque<_atom> > shiftSymmetryEquivalentAtoms(deque<deque<_atom> >& equivalent_atoms, xmatrix<double>& lattice, xvector<double>& translation, double& min_dist);

  // ******************************************************************************
  // RSTD Namespace Functions
  // ******************************************************************************
//  namespace rstd {
  typedef std::map<int, xvector<double> > hash;
  xvector<double> CrossPro(const xvector<double>& a, const xvector<double>& b);
  double DotPro(xvector<double> a, xvector<double> b);
//  double modulus(xvector<double> a);
  double modulus(vector<double> a);
//  //Linear Algebra Functions:
  void swap_rows(xmatrix<double>& M, int a, int b);
//  }

  xmatrix<double> concatenate(vector<xvector<double> >& V);
  xmatrix<double> concatenate(vector<xmatrix<double> >& V);

  void normalize(xvector<double>& v);
  void normalize(vector<double>& v);
  xmatrix<double> xvec2xmat(xvector<double> a, xvector<double> b, xvector<double> c);
  xmatrix<double> xvec2xmat(vector<xvector<double> > V);
  xmatrix<double> xvec2xmat(vector<xvector<double> > V, vector<double> R);

  xvector<double> extract_row(xmatrix<double> a, int row);
  xvector<double> dir(double a, double b, double c);
  xmatrix<double> a2m3x3(double* array);

  double get_random_double(double min, double max);

  bool vec_compare(vector<double> a, vector<double> b);
  bool vec_compare(xvector<double> a, xvector<double> b);

  //SYMMETRY OPERATIONS FUNCTIONS
  symfolder check_ccell(xstructure& xstr);
  vector<xvector<double> > expand_space_group_on_point(int sg, xvector<double> point);
  vector<xvector<double> > grid(double t);
  xvector<double> find_inversion_point(double tol, vector<eqatoms> poscar_atoms);
  vector<xvector<double> > symmetry_directions(char lattice_type);

  vector<Glide> mirror_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew);
  vector<Screw> triplet_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew);
  vector<Screw> twofold_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew);
  vector<xvector<double> > getLatticeVectorsFromOriginalMirrorOperations(vector<Glide>& old_mirrors, vector<Glide>& new_mirrors,
									 vector<xvector<double> >& lattice_vectors, bool& all_matched);
  vector<xvector<double> > getLatticeVectorsFromOriginalRotationOperations(vector<Screw>& old_rotations_twofold,
									   vector<Screw>& old_rotations_higher, vector<Screw>& new_rotations,
									   vector<xvector<double> >& twofold_lattice_vectors,
									   vector<xvector<double> >& rot_lattice_vectors, bool& all_matched);

  // COMBINATORICS FUNCTIONS
  void reduce_atom_deques(deque<_atom>& expanded, xmatrix<double>& lattice, double& min_dist);
//NOT IN SYM  vector<int> AllCombination41(int num, int total_num, int index);
//NOT IN SYM  vector<int> AllCombination42(int num, int total_num, vector<int>& str_in);
//NOT IN SYM  unsigned long int CombinationNr(int num, int total_num);
  vector<vector<int> > permute(int n);

  //LATTICE FUNCTIONS
  double length_lattice_vectors(xmatrix<double> lattice);
  double orthogonality_defect(xmatrix<double> xmat);

  //NUMBER THEORY FUNCTIONS
  long long int gcd(long long int u, long long int v);                             //Euclid's
  unsigned long long int gcd(unsigned long long int u, unsigned long long int v);  //Dijkstra's GCD Algorithm
  int gcdD(int u, int v);
  //[OBSOLETE]long long int cast2int(double d, long long int prec);
  //modulo reduce
  double mod_one(double d);
  double smallest_gt_min(double min, vector<double> vec);
  int smallest_gt_min_index(double min, int not_index1, int not_index2, vector<double> vec);

  //VISUALIZATION FUNCTIONS
  stringstream* latex_plot_lattice(xmatrix<double> L, string color);
  string plot_lattice(vector<string> files);

  //ATOM CLASS ELEMENT MANIPULATION FUNCTIONS
  //Do (1-atm.coord) for the atomic basis, for the column associated with the lattice vector row
  void minus_one(deque<_atom>& atomicBasis, int lvec);
  //when lattice vectors are permuted, you must swap columns in basis (in direct)
  void swap_columns(deque<_atom>& atomicBasis, int col1, int col2);
  void rearrange_columns(deque<_atom>& atomicBasis, int c1, int c2, int c3);
} //namespace SYM
/*
//A structure to store "decomposition" of points in plane--plane locations and distances from planes:
struct Proj {
  vector<_atom> inplane_locations;
  vector<double> distances_from_plane;
  xvector<double> plane_normal;
};
*/
void xb();  //print a break
template <class d>
void print(vector<d> vec) {
  for (uint i = 0; i < vec.size(); i++) {
    cout << vec[i] << endl;
  }
}
void xb();  //print a break
void print(vector<int> vec);
void print(const vector<xvector<double> >& points);
void print(const vector<vector<double> >& points);
void print(const vector<xvector<double> >& points, xmatrix<double> T);
void print(const vector<xmatrix<double> >& mats);
void print(const vector<_atom>& atoms);
void print(const deque<_atom>& atoms);

//Function to eliminate duplicate projections:
// DX TEST void eliminate_duplicates(Proj& P);

namespace SYM {
  //Operations on class-atoms
  vector<int> count_types(deque<_atom>& vatom);

  //EXPAND LATTICE/CRYSTAL FUNCTIONS
  vector<xvector<double> > expand_lattice_positive_only(int& a, int& b, int& c, xmatrix<double>& L);
  vector<xvector<double> > expand_lattice(int& a, int& b, int& c, xmatrix<double>& L);
  vector<xvector<double> > expand_cell(xmatrix<double>& L);
  deque<_atom> add_basis(vector<xvector<double> >& expanded_lattice_points, xmatrix<double>& L, xstructure& xstr);

  // CONVENTIONAL LATTICE VECTOR FUNCTIONS
  void orient(xstructure& xstr, bool update_atom_positions=true);
  xstructure ConventionalCell(xstructure& xstr, int& IT, int& cell_choice, bool& last_orientation, string& crystalsystem_prev,
			      xstructure& CrystOut_prev, vector<xmatrix<double> >& candidate_lattice_vectors_prev,
			      vector<char>& candidate_lattice_chars_prev, symfolder& checkops, bool& lattice_reformed,
			      vector<int>& lattice_pgroups, vector<xmatrix<double> >& lattice_sym_mats,
			      vector<xmatrix<double> >& crystal_sym_mats, bool& symmetry_found);
  bool findCubicLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec,
			vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars);
  bool findTrigonalLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec, vector<xvector<double> >& big_expanded,
			   vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars);
  bool findTetragonalLattice(vector<xvector<double> >& rot_lattice_vectors, vector<xvector<double> >& twofold_lattice_vectors,
			     vector<Screw>& rot_ops_vec, vector<Screw>& twofold_ops_vec, vector<xvector<double> >& big_expanded,
			     vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars);
  bool findMonoclinicLattice(vector<xvector<double> >& mirror_lattice_vectors, vector<xvector<double> >& twofold_lattice_vectors,
			     vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors,
			     vector<char>& candidate_lattice_chars, int& cell_choice); //DX 20180816 - added cell_choice
  bool findTriclinicLattice(xmatrix<double>& lattice, vector<xmatrix<double> >& candidate_lattice_vectors,
			    vector<char>& candidate_lattice_chars);
  bool findOrthorhombicLattice(vector<xvector<double> >& twofold_lattice_vectors, vector<xvector<double> >& mirror_lattice_vectors,
			       vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars);
  bool findRhombohedralLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec,
			       vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars);
  bool findRhombohedralSetting(vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors,
                               vector<char>& candidate_lattice_chars);
  bool determineLatticeCentering(vector<xvector<double> >& bravais_basis, int& bravais_count, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, vector<xvector<double> >& big_expanded, string& crystalsystem, vector<char>& candidate_lattice_chars);
  string getPearsonSymbol(char& centering, char& lattice_char, deque<_atom> atoms);
  bool getAtomGCD(deque<_atom>& atomic_basis, deque<deque<_atom> >& split_atom_types, int& GCD);
  deque<_atom> updateAtomPositions(deque<_atom>& atoms, Screw& S, xmatrix<double>& lattice);

  //RHOMBOHEDRAL OBVERSE/REVERSE FUNCTIONS
  bool isObverseSetting(xstructure& xstr, double& tolerance);
  bool isObverseSetting(xmatrix<double>& lattice, deque<_atom>& atomic_basis_, double& dist_nn_min, double& tolerance);
  bool isObverseSetting(xmatrix<double>& lattice, deque<deque<_atom> >& equivalent_atoms, double& dist_nn_min, double& tolerance);
  bool transformToObverse(xmatrix<double>& lattice, deque<_atom>& atoms);
  bool transformToObverse(xmatrix<double>& lattice, deque<deque<_atom> >& equivalent_atoms);

  //LATTICE VECTOR CHECK FUNCTIONS
  bool latticeVectorsSame(xvector<double>& a, xvector<double>& b, xvector<double>& c, double& tol);
  bool orientVectorsPositiveQuadrant(vector<xvector<double> >& lattice_vectors, double& tol);
  bool orientVectorsPositiveQuadrant(xvector<double>& vec, double& tol);
  bool alignLatticeWithXYZ(xvector<double>& a, xvector<double>& b, xvector<double>& c, double& tol);
  bool anyNullVectors(vector<xvector<double> >& vecs, double& tol);
  bool nullVector(xvector<double>& vec, double& tol);
} //namespace SYM

// external variables
namespace SYM {
  extern vector<xvector<double> > glideplanes;
  extern vector<xvector<double> > glidetrans;
  extern vector<string> glidesymbols;
  extern vector<xmatrix<double> > sym_mats;
  extern vector<string> symbol;
  extern vector<string> dirparam;
  extern hash sym_mats_direction;
  extern hash sym_mats_direction_hex;
  extern vector<xmatrix<double> > sym_mats_hex;
  extern vector<string> symbol_hex;
  extern vector<string> dirparam_hex;
  extern vector<xvector<double> > glideplanes_hex;
  extern vector<xvector<double> > glidetrans_hex;
  extern vector<string> glidesymbols_hex;
  extern vector<int> index_cubic;
  extern vector<int> index_hex;
  extern vector<int> index_rhom;
  extern vector<int> index_tetr;
  extern vector<int> index_ortho;
  extern vector<int> index_mono_b;
  extern vector<int> index_mono_c;
  extern vector<int> index_tric;
  extern vector<vector<symop> > generators;
  extern vector<int> sgindex;
  extern vector<xmatrix<double> > sym_mats;
  extern vector<string> symbol;
  extern hash sym_mats_direction;
  extern vector<string> gl_sgs;
}


#endif
