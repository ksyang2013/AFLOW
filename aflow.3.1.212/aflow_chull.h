// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *           Aflow COREY OSES - Duke University 2013-2018                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu
// Previous versions also written by Eric Perim and Eric Gossett

#ifndef _AFLOW_CHULL_H_
#define _AFLOW_CHULL_H_

// UNITS
const char _std_ = 'S';  // standard aflow/vasp units
const char _m_ = 'm';    // convert to milli-

// FORMATS
const char _apool_ = 'a';  // apool
const char _json_ = 'j';   // standard json
const char _pdf_ = 'p';    // pdf
const char _txt_ = 't';    // plain text
const char _web_ = 'w';    // web json
const char _latex_ = 'l';    // latex
const char _gnuplot_ = 'g';  // gnuplot

// REDUCTION MODES
const char _frac_ = 'f';  //fractional
const char _gcd_ = 'g';   //gcd
const char _none_ = 'n';  //none

// DEFAULTS
const int CHULL_PRECISION = 8;                          //must be less than _precision_ in aflow_xatom.cpp, which is currently set to 14
const int FULL_PRECISION = 15;                          //max printing precision
const int COEF_PRECISION = 4;
const double ZERO_TOL = pow(10,-CHULL_PRECISION);       //lower bound for absolute resolution of floats, significant differences among floats should be well above this threshold
const double ROUNDOFF_TOL = pow(10,-CHULL_PRECISION+2); //make less strigent so we don't get 1e-6
const double ZERO_FULL_TOL = pow(10,-FULL_PRECISION);
const double ZERO_COEF_TOL = pow(10,-COEF_PRECISION);
const double ENERGY_TOL = 0.015;                        //eV, CO NOTES - structures within this thresold may be equivalent, I've seen as large as 5meV, keep at 15 to be safe
const int ZERO_RANGE_TOL = 1;
//[CO 180316 - moved to aflowrc]const uint BINARY_ENTRIES_THRESHOLD = 200;

// CO 180419 - moved to AFLOWRuntimeError and AFLOWLogicError
//namespace chull {
//  class CHullRuntimeError : public std::runtime_error {
//    public:
//      CHullRuntimeError(const std::string& function,const std::string& message);
//      CHullRuntimeError(const std::string& function,std::stringstream& message);
//      string where();
//      ~CHullRuntimeError() throw() {};
//    private:
//      string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//  };
//  class CHullLogicError : public std::logic_error {
//    public:
//      CHullLogicError(const std::string& function,const std::string& message);
//      CHullLogicError(const std::string& function,std::stringstream& message);
//      string where();
//      ~CHullLogicError() throw() {};
//    private:
//      string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//  };
//} // namespace chull

namespace chull {
  class ChullPoint; //forward declaration
  class ConvexHull; //forward declaration
  bool convexHull(const aurostd::xoption& vpflow);
  ////////////////////////////////////////////////////////////////////////////////
  // gets path to redirect output
  string getPath(bool add_backslash=true);
  string getPath(const aurostd::xoption& vpflow, ostream& oss=cout, bool silent=true); // CO 180220
  string getPath(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& oss=cout, bool silent=true);  // CO 180220
  string getPath(string _path, ostream& oss=cout, bool silent=true); // CO 180220
  string getPath(string _path, ofstream& FileMESSAGE, ostream& oss=cout, bool silent=true);  // CO 180220
  ////////////////////////////////////////////////////////////////////////////////
  //logs which flags are on
  void flagCheck(aurostd::xoption& vpflow, const vector<string>& velements, ostream& oss=cout, bool silent=false);
  void flagCheck(aurostd::xoption& vpflow, const vector<string>& velements, ofstream& FileMESSAGE, ostream& oss=cout, bool silent=false);
  ////////////////////////////////////////////////////////////////////////////////
  //stability criterion calculation
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const vector<string>& velements,const string& auid,double& scriterion,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const vector<string>& velements,const string& auid,double& scriterion,ofstream& FileMESSAGE,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ofstream& FileMESSAGE,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,double& scriterion,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,double& scriterion,ofstream& FileMESSAGE,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const string& auid,double& scriterion,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const string& auid,double& scriterion,ofstream& FileMESSAGE,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,vector<double>& vscriterion,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,vector<double>& vscriterion,ofstream& FileMESSAGE,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ostream& oss=cout);
  //[OBSOLETE now inside ConvexHull]bool calculateStabilityCriterion(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<string>& vauid,vector<double>& vscriterion,ofstream& FileMESSAGE,ostream& oss=cout);
  ////////////////////////////////////////////////////////////////////////////////
  // returns value in desired units
  double convertUnits(double value, char units=_std_);
  double H_f_atom(const ChullPoint& point, char units=_std_);
  double H_f_atom(const aflowlib::_aflowlib_entry& entry, char units=_std_);
  double T_S(const ChullPoint& point);
  double T_S(const aflowlib::_aflowlib_entry& entry);
  double isoMaxLatentHeat(const ChullPoint& point, double x, char units=_std_);
  double isoMaxLatentHeat(const aflowlib::_aflowlib_entry& entry, double x, char units=_std_);
  ////////////////////////////////////////////////////////////////////////////////
  int roundDouble(double doub, int multiple, bool up);
  ////////////////////////////////////////////////////////////////////////////////
  bool greaterEqualZero(double val);
  bool lessEqualZero(double val);
  bool notPositive(double val,bool soft_cutoff,double tol=ZERO_TOL);
  bool notNegative(double val,bool soft_cutoff,double tol=ZERO_TOL);
  bool zeroWithinTol(double val,double tol=ZERO_TOL);
  bool nonZeroWithinTol(double val,double tol=ZERO_TOL);
  bool subspaceBelongs(const xvector<int>& space,const xvector<int>& subspace);
  bool correctSignVerticalDistance(double dist_2_hull,bool should_be_positive);
} // namespace chull

// CO 180420 - moved to xStream (aflow.h)
//namespace chull {
//  class ChullClassTemplate {
//    public:
//      //NECESSARY PUBLIC CLASS METHODS - START
//      //constructors - START
//      ChullClassTemplate();
//      //constructors - STOP
//      ~ChullClassTemplate();
//      //NECESSARY PUBLIC CLASS METHODS - END
//
//      ////flags
//      //aurostd::xoption m_cflags;
//      //_aflags m_aflags;                  //used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags
//
//    protected:
//      //NECESSARY private CLASS METHODS - START
//      void free();
//      void freeStreams();
//      void freeAll();
//      //NECESSARY END CLASS METHODS - END
//      
//      //logger variables
//      ostream* p_oss;
//      ofstream* p_FileMESSAGE;
//      bool m_new_ofstream;  //for deletion later
//      
//      //general setters
//      //void setDefaultCFlags();
//      //void setCFlags(const aurostd::xoption& vpflow);
//      void setOFStream(ofstream& FileMESSAGE);
//      void setOSS(ostream& oss);
//      //void setDirectory();
//  };
//} // namespace chull

namespace chull {
  class ChullPoint : public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      ChullPoint(ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      ChullPoint(const xvector<double>& coord,ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      ChullPoint(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ostream& oss=cout,bool formation_energy_coord=true);
      ChullPoint(ofstream& FileMESSAGE,ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      ChullPoint(const xvector<double>& coord,ofstream& FileMESSAGE,ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      ChullPoint(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ofstream& FileMESSAGE,ostream& oss=cout,bool formation_energy_coord=true);
      ChullPoint(const ChullPoint& b);
      //constructors - STOP
      ~ChullPoint();
      const ChullPoint& operator=(const ChullPoint& other);
      void HullCopy(const ChullPoint& b);  //copies b without entry
      bool operator<(const ChullPoint& other) const;
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP
      
      //general attributes
      bool m_initialized;
      xvector<double> m_coords;           //most general hull coordinates
      bool m_has_stoich_coords;
      aflowlib::_aflowlib_entry m_entry;
      bool m_has_entry;
      bool m_formation_energy_coord;
      bool m_is_artificial;
      
      //for organization of points
      uint m_i_nary;    //stoich_coords only
      uint m_i_alloy;   //stoich_coords_only
      uint m_i_coord_group;
      
      //stoich_coords only!
      xvector<double> s_coords;         //stoich_coords, m_coords[0:rows-1]+hidden_dimension
      xvector<double> c_coords;         //composition coords, similar to stoich_coords, but with integers (except for POCC)
      xvector<int> m_elements_present;  //1 if s_coords[i]>ZERO_TOL, zero otherwise, that way, nary=sum(m_elements_present)

      //calculate per hull
      bool m_is_on_hull;  //one max per coordgroup
      bool m_is_g_state;  //one max per coordgroup, must have m_entry (artificialMap())
      bool m_is_equivalent_g_state; //can be many, includes original g_state
      bool m_is_sym_equivalent_g_state; //can be many, includes original g_state
      double m_dist_2_hull; //warning, this is not TRUE dist to hull (facet), this is vertical distance
      //[OBSOLETE - reduce by _frac_ always! so use coord_group values]xvector<double> m_decomp_coefs; //un-reduced coefficients here, reduced exist at the m_coord_groups level
      double m_stability_criterion;
      
      //flags - MOVED TO xStream
      //aurostd::xoption m_cflags;
      //_aflags m_aflags;                  //used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags

      //initialization methods
      bool initialize(ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      bool initialize(const xvector<double>& coord,ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      bool initialize(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ostream& oss=cout,bool formation_energy_coord=true);
      bool initialize(ofstream& FileMESSAGE,ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      bool initialize(const xvector<double>& coord,ofstream& FileMESSAGE,ostream& oss=cout,bool has_stoich_coords=false,bool formation_energy_coord=false,bool is_artificial=false);
      bool initialize(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,ofstream& FileMESSAGE,ostream& oss=cout,bool formation_energy_coord=true);
      
      //getters
      bool isWithinHalfHull(bool lower_hull=true) const;
      bool isGState() const;
      xvector<double> getStoichiometricCoords() const; //get stoichiometric coordinates (sans energetic coordinate)
      xvector<double> getTruncatedReducedCoords(const xvector<int>& elements_present,char reduce_mode=_frac_) const;
      xvector<double> getTruncatedCoords(const xvector<double>& coords,const xvector<int>& elements_present) const; //truncated arbitrary coords
      xvector<double> getTruncatedSCoords(const xvector<int>& elements_present) const; //truncate stoichiometry
      xvector<double> getTruncatedCCoords(const xvector<int>& elements_present,bool reduce=true) const; //similar to truncated stoichiometry, but in integer form (if not POCC)
      xvector<double> getReducedCCoords() const;  //reduce by gcd()
      uint loadXstructures(bool relaxed_only);
      bool getMostRelaxedXstructure(xstructure& xstr) const;
      double getLastCoord() const;
      uint getDim() const;
      bool isUnary() const;
      double getFormationEnthalpy() const;
      double getEntropicTemperature() const;
      const vector<string>& getVSG() const;
      const string& getSG() const;
      double getDist2Hull(char units=_std_) const;
      double getStabilityCriterion(char units=_std_) const;
      double getRelativeStabilityCriterion() const;
      
      //setters
      void setHullCoords();
      void setHullCoords(const xvector<double>& coords);
      void setHullCoords(const xvector<int>& elements_present);
     
      //general methods
      bool entryIdentical(const aflowlib::_aflowlib_entry& other) const;

      friend class ChullFacet;  //ConvexHull defines everything!
      friend class ConvexHull;  //ConvexHull defines everything!
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const ChullPoint& b);
      //NECESSARY END CLASS METHODS - END
      
      //logger variables - MOVED TO xStream
      //ostream* p_oss;
      //ofstream* p_FileMESSAGE;
      //bool m_new_ofstream;  //for deletion later
      
      //general setters - MOVED TO xStream
      //void setOFStream(ofstream& FileMESSAGE);
      //void setOSS(ostream& oss);
      
      //temporary variables that come with every new calculate() command
      xvector<double> h_coords;  //dummy coord that changes per convex hull iteration, mere projection

      //initialization methods
      void initializeCoords(const xvector<double>& coord,bool formation_energy_coord=false);
      void initializeCoords(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,bool formation_energy_coord=true);
      void addEntry(const aflowlib::_aflowlib_entry& entry);
      void setGenCoords(const xvector<double>& coord,bool formation_energy_coord=false);
      void setGenCoords(const vector<string>& velements,const aflowlib::_aflowlib_entry& entry,bool formation_energy_coord=true);

      //hull specific methods
      vector<uint> getRelevantIndices(const xvector<int>& elements_present) const;
      void setStoichCoords();
      void cleanPointForHullCalc();  //clean hull state of point
      void cleanPointForHullTransfer();
  };
} // namespace chull

namespace chull {
  //allows chullfacet to be (mostly) INDEPENDENT of convex hull
  //there is a HUGE advantage to considering points by index, rather than coordinate
  //comparisons/groupings are easier
  //make sure that such comparisons are being made per facets of the same hull, never among different hulls!
  class FacetPoint {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      FacetPoint();
      FacetPoint(const ChullPoint& point,uint index,bool full_copy=true);
      FacetPoint(const FacetPoint& b);
      //constructors - STOP
      ~FacetPoint();
      const FacetPoint& operator=(const FacetPoint& other);
      bool operator<(const FacetPoint& other) const;
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP

      //attributes
      bool m_initialized;
      uint ch_index;       //belongs to chull
      ChullPoint ch_point; //take only what you need from chullpoint, don't copy the whole thing (has entry which is large)
      
      //initializer
      void initialize(const ChullPoint& point,uint index,bool full_copy=true);
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const FacetPoint& b);
      //NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

//nice sorting for points we know have energy vs. stoich coords
namespace chull {
  struct sortThermoPoints{
    sortThermoPoints(bool sort_stoich_ascending=true,bool sort_energy_ascending=true) : 
      m_sort_stoich_ascending(sort_stoich_ascending),m_sort_energy_ascending(sort_energy_ascending) {};
    bool m_sort_stoich_ascending; //good for sorting points in facets (lower vs. upper hemisphere), if facet in lower hemisphere, then sort ascending order (CLOCKWISE, graphing)
    bool m_sort_energy_ascending; //good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
    bool operator() (const FacetPoint& fpi,const FacetPoint& fpj) const;
    bool operator() (const ChullPoint& ci,const ChullPoint& cj) const;
  };
} // namespace chull

namespace chull {
  class ChullFacet: public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      ChullFacet(ostream& oss=cout);
      ChullFacet(ofstream& FileMESSAGE, ostream& oss=cout);
      ChullFacet(const ChullFacet& b);
      //constructors - STOP
      ~ChullFacet();
      const ChullFacet& operator=(const ChullFacet& other);
      bool operator<(const ChullFacet& other) const;
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP

      //general attributes
      bool m_initialized;
      vector<FacetPoint> m_vertices;
      uint m_dim;
      bool m_has_stoich_coords;
      bool m_formation_energy_coord;
      double m_content; //length of line, area of triangle, vol of tetrahedron...
      vector<xvector<double> > m_directive_vectors;
      xvector<double> m_normal;
      double m_offset;
      xvector<double> m_facet_centroid;
      xvector<double> m_hull_reference;
      bool m_hypercollinear;
      bool m_is_vertical;
      bool m_is_artificial;
      bool m_in_lower_hemisphere;                //this determines how we sort wrt stoich (descending in lower_hemisphere)
      vector<ChullFacet> m_ridges;

      //flags - MOVED TO xStream
      //aurostd::xoption m_cflags;
      //_aflags m_aflags;                  //used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags
      
      //getters
      vector<uint> getCHIndices() const;

      //setters
      void addVertex(const ChullPoint& point,uint index=AUROSTD_MAX_UINT);
      void addVertex(const FacetPoint& fp);
      void initialize(const xvector<double>& hull_centroid,uint hull_dim,bool check_validity=true); //needs some hull data (context) to align correctly

      //methods using facets to build hull
      bool shareRidge(const ChullFacet& other) const;
      bool isPointOnFacet(const FacetPoint& fp) const;
      bool isPointOnFacet(uint i_point) const;
      bool isPointOutside(const FacetPoint& fp) const;
      bool isPointOutside(const ChullPoint& point) const;
      double getSignedPointPlaneDistance(const ChullPoint& point) const;
      double getSignedPointPlaneDistance(const xvector<double>& point) const;
      double getSignedVerticalDistanceToZero(const ChullPoint& point) const;
      double getSignedVerticalDistanceToZero(const xvector<double>& point) const;
      double getSignedVerticalDistance(const ChullPoint& point) const;
      double getSignedVerticalDistance(const xvector<double>& point) const;
      
      friend class ConvexHull;  //ConvexHull defines everything!
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const ChullFacet& b);
      //NECESSARY END CLASS METHODS - END

      //logger variables - MOVED TO xStream
      //ostream* p_oss;
      //ofstream* p_FileMESSAGE;
      //bool m_new_ofstream;  //for deletion later
      
      //temporary variables that come with every new calculate() command
      bool f_visited;                         //temporary per calculation(), has it been visited?
      vector<FacetPoint> f_outside_set;
      FacetPoint f_furthest_point;
      vector<uint> f_neighbors;

      //general setters
      void create(ostream& oss=cout);
      void create(ofstream& FileMESSAGE,ostream& oss=cout);
      //MOVED TO xStream
      //void setCFlags(const aurostd::xoption& vpflow);
      //void setOFStream(ofstream& FileMESSAGE);
      //void setOSS(ostream& oss);
      
      //validation methods - we split up as we build facets in pieces, and not all work
      //sometimes we find a (d-1)flat facet, so we remove
      //we want to return false vs. throw exception appropriately
      bool hasValidPoints(string& error);
      void setContent();
      void setDirectiveVectors(bool check_validity=true);
      bool pointsMatchDirectiveVectors(string& error);
      bool hasValidDirectiveVectors(string& error);
      bool hasCollinearVectors(bool check_validity=true);
      bool isValid(string& error);

      //initialization methods, again split up because we build in pieces
      void setNormal(bool check_validity=true);
      void setOffset();
      void setCentroid();
      void setVertical();
      void setArtificial();
      void alignNormalInward();
      void setHemisphere();
      void setFurthestPoint();
      void setRidges();
      
      //hull specific
      void cleanFacet();  //clean state of facet
  };
} // namespace chull

//Just another group of ChullPoints
namespace chull {
  class CoordGroup {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      CoordGroup();
      CoordGroup(const xvector<double>& coord,bool has_stoich_coords);
      CoordGroup(const CoordGroup& b);
      //constructors - STOP
      ~CoordGroup();
      const CoordGroup& operator=(const CoordGroup& other);
      bool operator<(const CoordGroup& other) const;
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP
      
      //initializer
      void initialize(const xvector<double>& coord,bool has_stoich_coords);

      //getters
      xvector<int> getElementsPresent() const;
      uint getDim() const;
      
      //attributes
      bool m_initialized;
      xvector<double> m_coords;               // ChullPoint.getStoichiometricCoords()
      vector<uint> m_points;                  // points to ChullPoints
      bool m_has_stoich_coords;
      bool m_has_artificial_unary;
      bool m_is_on_hull;                      // is ground state
      uint m_hull_member;                     // points to ChullPoints
      uint m_ref_state;                       // first in order (artificial map)
      vector<uint> m_candidate_hull_points;   // also points to ChullPoints, refers to extremes of CoordGroup (half_hull's only has extrema)
      
      //for organization of points
      uint m_i_nary;    //stoich_coords only
      uint m_i_alloy;   //stoich_coords_only

      //the follow properties are ONLY found if stoich_coords + half_hull + artificial_points
      //assumes ONE gstate per coordgroup
      uint m_nearest_facet;         // which facet is directly below/above me?
      double m_nearest_distance;    // how close is the nearest facet?
      vector<uint> m_decomp_phases;
      xvector<double> m_decomp_coefs;
      vector<vector<uint> > m_equilibrium_phases;
      bool m_calculated_equivalent_g_states;
      vector<uint> m_equivalent_g_states; //structure comparison
      vector<uint> m_sym_equivalent_g_states; //structure comparison
      double m_stability_criterion; //g-states only
      bool m_icsd_g_state;          //whether icsd exists among equivalent states
      uint m_i_canonical_icsd;      //canonical icsd entry (lowest number)

      friend class ConvexHull;  //ConvexHull defines everything!
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const CoordGroup& b);
      //NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

//simple, nice organization of points
//usually, we further organize points into CoordGroups EXCEPT for hull_points_unary
//hull_points_unary === single CoordGroup, trivial
namespace chull {
  class Alloy {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      Alloy();
      Alloy(const xvector<int>& elements_present);
      Alloy(const Alloy& b);
      //constructors - STOP
      ~Alloy();
      const Alloy& operator=(const Alloy& other);
      bool operator<(const Alloy& other) const;
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP
      
      //initializer
      void initialize(const xvector<int>& elements_present);
      bool belongs2Hull(const xvector<int>& hull_elements_present) const;  //checks if alloy is a member of the hull (via hull_elements_present)

      //attributes
      bool m_initialized;
      xvector<int> m_elements_present;
      uint m_dim;                                       //projected hull dimensionality
      vector<uint> m_coord_groups;
      vector<uint> m_facets;
      
      friend class ConvexHull;  //ConvexHull defines everything!
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const Alloy& b);
      //NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

//simple, nice organization of points
namespace chull {
  class Nary {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      Nary();
      Nary(uint nary);
      Nary(const Nary& b);
      //constructors - STOP
      ~Nary();
      const Nary& operator=(const Nary& other);
      bool operator<(const Nary& other) const;
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP
      
      //initializer
      void initialize(uint nary);

      //attributes
      bool m_initialized;
      uint nary;
      vector<Alloy> m_alloys;
      
      friend class ConvexHull;  //ConvexHull defines everything!
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const Nary& b);
      //NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

namespace chull {
  class ConvexHull: public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      ConvexHull(ostream& _oss=cout);
      ConvexHull(string alloy,ostream& _oss=cout);
      ConvexHull(const vector<string>& velements,ostream& _oss=cout);
      ConvexHull(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& _oss=cout);
      ConvexHull(const vector<xvector<double> >& vcoords,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const vector<ChullPoint>& vpoints,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(string alloy,ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const aurostd::xoption& vpflow,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,string alloy,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const aurostd::xoption& vpflow,ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,string alloy,ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& _oss=cout);
      ConvexHull(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      ConvexHull(const ConvexHull& b);
      //constructors - STOP
      ~ConvexHull();
      const ConvexHull& operator=(const ConvexHull& other);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP
    
      //attributes
      bool m_initialized;
      vector<string> m_velements;
      vector<uint> m_icsd_entries;      //save these for checks for experimental validation
      vector<ChullPoint> m_points;
      vector<Nary> m_naries;             //one really nice container for indicies for points, naries, then alloys, then CoordGroups, then points
      vector<CoordGroup> m_coord_groups; //index to m_points
      uint m_dim;                        //full dimensionality
      bool m_half_hull;                  //huge speed up
      bool m_lower_hull;                 //true == lower half only, false == upper half only
      bool m_has_stoich_coords;          //[0,1]
      bool m_add_artificial_unaries;     //force unaries to 0
      bool m_thermo_hull;                //enthalpy_formation/entropic_temperature with stoich_coords and artificial unaries
      bool m_formation_energy_hull;
      vector<ChullFacet> m_facets;
      vector<uint> m_i_facets;           //contains indices of most recently calculated facets, points to m_facets
      bool m_sort_energy_ascending;      //true for lower half hulls
      
      //flags
      aurostd::xoption m_cflags;
      _aflags m_aflags;                  //used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags
      
      //initialization methods
      bool initialize(ostream& _oss=cout);
      bool initialize(string alloy,ostream& _oss=cout);
      bool initialize(const vector<string>& velements,ostream& _oss=cout);
      bool initialize(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& _oss=cout);
      bool initialize(const vector<xvector<double> >& vcoords,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const vector<ChullPoint>& vpoints,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(string alloy,ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const aurostd::xoption& vpflow,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,string alloy,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,const vector<string>& velements,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const aurostd::xoption& vpflow,ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,string alloy,ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries,ofstream& FileMESSAGE,ostream& _oss=cout);
      bool initialize(const aurostd::xoption& vpflow,const vector<xvector<double> >& vcoords,ofstream& FileMESSAGE,ostream& _oss=cout,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool initialize(const aurostd::xoption& vpflow,const vector<ChullPoint>& vpoints,const vector<string>& velements,ofstream& FileMESSAGE,ostream& _oss=cout,bool formation_energy_hull=false,bool add_artificial_unaries=false);

      //initialize points ONLY
      void initializePoints(string alloy);
      void initializePoints(const vector<string>& velements);
      void initializePoints(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries);
      void initializePoints(const vector<xvector<double> >& vcoords,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      void initializePoints(const vector<ChullPoint>& vpoints,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      void initializePoints(const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_energy_hull=false,bool add_artificial_unaries=false);

      //getters
      uint getDim() const;
      uint getEntriesCount(bool only_within_half_hull=true) const;
      uint getEntriesCount(uint i_nary,bool only_within_half_hull=true) const;
      uint getEntriesCount(uint i_nary,uint i_alloy,bool only_within_half_hull=true) const;
      vector<vector<uint> > getHullSizes(bool only_within_half_hull=true) const;
      uint getGStateCount() const;
      uint getGStateCount(uint i_nary) const;
      vector<uint> getHullPoints(bool sort_stoich_ascending=true) const;
      vector<uint> getGStates(bool include_unaries=true,bool sort_stoich_ascending=true) const;
      uint getUnaryGState(uint i_alloy) const;
      bool isViablePoint(uint i_point) const;
      bool isViableGState(uint g_state) const;
      bool findPoint(const string& auid,uint& i_point) const;
      bool findPoint(const xvector<double>& coords,uint& i_point) const;
      bool getNariesIndex(uint i_point,uint& i_nary,uint& i_alloy,uint& i_coord_group,bool redo=false) const;
      bool getNariesIndex(const ChullPoint& point,uint& i_nary,uint& i_alloy,uint& i_coord_group,bool redo=false) const;
      bool getCoordGroupIndex(uint i_point,uint& i_coord_group,bool redo=false) const;
      bool getCoordGroupIndex(const ChullPoint& point,uint& i_coord_group,bool redo=false) const;
      bool getCoordGroupIndex(const xvector<double>& r_coords,uint& i_coord_group) const;
      bool getAlloyIndex(const ChullPoint& point,uint& i_nary,uint& i_alloy,bool redo=false) const;
      bool getAlloyIndex(const CoordGroup& cg,uint& i_nary,uint& i_alloy,bool redo=false) const;
      bool getAlloyIndex(const xvector<int>& elements_present,uint& i_nary,uint& i_alloy) const;
      uint artificialMap(uint i_point) const;
      uint getNearestFacetVertically(const vector<uint>& i_facets,const ChullPoint& point) const;
      uint getNearestFacetVertically(const vector<uint>& i_facets,const xvector<double>& point) const;
      double getSignedVerticalDistanceWithinCoordGroup(uint i_coord_group,uint i_point) const;
      double getSignedVerticalDistanceWithinCoordGroup(uint i_coord_group,const ChullPoint& point) const;
      double getDistanceToHull(uint i_point,bool redo=false) const;
      double getDistanceToHull(const ChullPoint& point,bool redo=false) const;
      vector<double> getDistancesToHull(const vector<string>& vauid,bool redo=false) const;
      vector<uint> extractDecompositionPhases(const ChullFacet& facet) const;
      vector<uint> getDecompositionPhases(uint i_point) const;
      vector<uint> getDecompositionPhases(const ChullPoint& point) const;
      xvector<double> getDecompositionCoefficients(uint i_point,char reduce_mode=_frac_) const;
      xvector<double> getDecompositionCoefficients(const ChullPoint& point,char reduce_mode=_frac_) const;
      xvector<double> getDecompositionCoefficients(uint i_point,const vector<uint>& decomp_phases,char reduce_mode=_frac_) const;
      xvector<double> getDecompositionCoefficients(const ChullPoint& point,const vector<uint>& decomp_phases,char reduce_mode=_frac_) const;
      vector<uint> getAdjacentFacets(uint hull_member,bool ignore_hypercollinear=true,bool ignore_vertical=true,bool ignore_artificial=true) const;
      vector<vector<uint> > getEquilibriumPhases(uint hull_member) const;
      vector<uint> getEquivalentGStates(uint g_state) const;
      vector<uint> getSymEquivalentGStates(uint g_state) const;
      double getStabilityCriterion(const string& cauid) const;
      double getStabilityCriterion(uint cpoint) const;
      vector<double> getStabilityCriterion(const vector<string>& vcauid) const;
      vector<double> getStabilityCriterion(const vector<uint>& vcpoint) const;

      //writer
      bool write(char mode=_pdf_) const;
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const ConvexHull& b);
      //NECESSARY END CLASS METHODS - END
      
      //logger variables - MOVED TO xStream
      //ostream* p_oss;
      //ofstream* p_FileMESSAGE;
      //bool m_new_ofstream;  //for deletion later
      
      //aflowrc definitions we don't want to redefine over and over again, just once at initialization
      bool m_allow_all_formation_energies;
      vector<string> m_allowed_dft_types;
      
      //temporary variables that come with every new calculate() command
      uint h_dim;                                       //projected hull dimensionality
      xvector<int> m_elements_present;                  //stoich_coords only
      vector<uint> h_points;                            //current hull points, subset of m_points
      xvector<double> h_centroid;                       //centroid of ALL these points
      xvector<double> h_reference;                      //centroid of INITIAL points on hull (more stable reference)
      vector<ChullFacet> h_facets;                      //current hull facets, real facets are stored in Alloys (m_naries)
      vector<uint> h_visible_facets;
      vector<ChullFacet> h_horizon_ridges;

      //general setters
      void setDefaultCFlags();
      void setCFlags(const aurostd::xoption& vpflow);
      void setDirectory();
      // MOVED TO xStream
      //void setOFStream(ofstream& FileMESSAGE);
      //void setOSS(ostream& oss);
      
      //wrapper for try/catch's
      bool createHull(string alloy);
      bool createHull(const vector<string>& velements);
      bool createHull(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries);
      bool createHull(const vector<xvector<double> >& vcoords,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool createHull(const vector<ChullPoint>& vpoints,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      bool createHull(const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_energy_hull=false,bool add_artificial_unaries=false);

      //methods associated with points
      bool entryValid(const aflowlib::_aflowlib_entry& entry,bool ignore_bad_database=true) const;
      bool entryValid(const aflowlib::_aflowlib_entry& entry,string& reason,bool ignore_bad_database=true) const;
      bool entryValid(const aflowlib::_aflowlib_entry& entry,string& reason,char& LOGGER_TYPE,bool ignore_bad_database=true) const;
      bool entryUnique(const vector<uint>& unique_entries,const aflowlib::_aflowlib_entry& entry2,string& canonical_auid) const;
      void addArtificialUnaries(uint dim);
      void addArtificialUnaries();
      void loadPoints(string alloy);
      void loadPoints(const vector<string>& velements);
      void loadPoints(const vector<string>& velements,const vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries);
      void loadPoints(const vector<xvector<double> >& vcoords,bool has_stoich_coords=false,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      void loadPoints(const vector<ChullPoint>& vpoints,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      void loadPoints(const vector<ChullPoint>& vpoints,const vector<string>& velements,bool formation_energy_hull=false,bool add_artificial_unaries=false);
      void calculateOutlierThreshold(const vector<double>& _energies,double& upper_threshold,double& lower_threshold);
      void calculateOutlierThreshold(const xvector<double>& energies,double& upper_threshold,double& lower_threshold);
      vector<uint> calculateOutliers(const vector<uint>& points_to_consider);
      vector<uint> getOutliers();
      vector<uint> getOutliers(const xvector<int>& elements_present);
      vector<uint> findArtificialPoints(uint i_coord_group);
      uint findArtificialUnary(uint i_coord_group);
      void organizeHullPoints(uint i_coord_group);
      void organizeHullPoints();
      void initializeNaries();
      void structurePoints();
      vector<string> alloyToElements(const ChullPoint& point) const;
      vector<string> alloyToElements(uint i_nary,uint i_alloy) const;
      void checkStructurePoints();
      
      //other methods associated with calculating hull
      void addPointToFacet(ChullFacet& facet,uint i_point);
      void initializeFacet(ChullFacet& facet,bool check_validity=true);
      bool isPointOutsideFacet(ChullFacet& facet,uint i_point);
      uint getExtremePoint(uint dim);
      uint getExtremePoint(uint dim,const vector<FacetPoint>& points_to_avoid);
      void setCentroid();
      vector<FacetPoint> getInitialExtremePoints();
      void setNeighbors();
      void createInitializeSimplex();
      void setVisibleFacets(uint i_facet);
      void setHorizonRidges();
      uint createNewFacets(FacetPoint furthest_point);   //returns true if new facets were, in fact created (not necessary at each iteration)
      void updateOutsideSet(uint new_facet_count);
      void deleteVisibleFacets();
      void removeDuplicateHullPoints();
      void calculateFacets();
      const xvector<int>& getElementsPresent(uint i_nary,uint i_alloy) const;
      void setElementsPresent(uint i_nary,uint i_alloy);
      void addRelevantUnariesToHullCalculation(uint i_nary,uint i_alloy);
      void addRelevantUnariesToHullCalculation(xvector<int>& elements_present);
      void addLowerDimensionPointsToHullCalculation(uint i_nary_max);
      void addPointToHullCalculation(uint i_point);
      void addPointToHullCalculation(uint i_point,xvector<int>& elements_present);
      void preparePointsForHullCalculation(uint i_nary,uint i_alloy);
      void preparePointsForHullCalculation();
      const vector<uint>& getRelevantFacets(uint i_nary,uint i_alloy) const;
      void setHullMembers();
      void setHullMembers(uint i_nary,uint i_alloy);
      void setHullMembers(const vector<uint>& i_facets);
      void setNearestFacet(uint i_nary,uint i_alloy,uint i_coord_group);
      bool hullDistanceExtractionRequired() const;
      bool thermoPropertiesExtractionRequired() const;
      bool thermoPostProcessingExtractionRequired() const;
      void setDistancesToHull(uint i_nary,uint i_alloy);
      void setDistancesToHull(uint i_nary,uint i_alloy,uint i_coord_group);
      void setDecompositionPhases(uint i_nary,uint i_alloy,uint i_coord_group);
      void setDecompositionCoefficients(uint i_nary,uint i_alloy,uint i_coord_group);
      void setOffHullProperties(uint i_nary,uint i_alloy);
      void setEquilibriumPhases(uint i_nary,uint i_alloy,uint i_coord_group);
      bool energiesDiffer(uint i_point1,uint i_point2,bool strict=true) const;
      bool spacegroupsDiffer(uint i_point1,uint i_point2,bool strict=true) const;
      bool structuresEquivalent(uint i_point1,uint i_point2) const;
      bool isICSD(uint i_point) const;
      void setEquivalentGStates(uint i_nary,uint i_alloy,uint i_coord_group);
      void setSymEquivalentGStates(uint i_nary,uint i_alloy,uint i_coord_group);
      void setOnHullProperties(uint i_nary,uint i_alloy);
      void storeHullData(uint i_nary,uint i_alloy);
      void storeHullData();
      void extractThermodynamicProperties(uint i_nary,uint i_alloy);
      void thermodynamicsPostProcessing();
      void calculate();
      void setStabilityCriterion();
      void cleanHull(); //clean state of hull
      
      //writer functions
      string prettyPrintCompound(const ChullPoint& point,char reduce_mode=_gcd_,bool exclude1=true,char mode=_latex_) const;
      string prettyPrintCompound(const aflowlib::_aflowlib_entry& entry,char reduce_mode=_gcd_,bool exclude1=true,char mode=_latex_) const;
      string prettyPrintCompound(const vector<string>& vspecies,const vector<double>& vcomposition,char reduce_mode=_gcd_,bool exclude1=true,char mode=_latex_) const;
      string prettyPrintCompound(const vector<string>& vspecies,const xvector<double>& vcomposition,char reduce_mode=_gcd_,bool exclude1=true,char mode=_latex_) const;
      string getICSDNumber(uint i_point,bool remove_suffix=true) const;
      string getICSDNumber(const ChullPoint& point,bool remove_suffix=true) const;
      string getICSDNumber(const aflowlib::_aflowlib_entry& entry,bool remove_suffix=true) const;
      string prettyPrintPrototype(const ChullPoint& point, bool double_back_slash,bool icsd_label_skim=false) const;
      string prettyPrintPrototype(const aflowlib::_aflowlib_entry& entry, bool double_back_slash,bool icsd_label_skim=false) const;
      string fixStringLatex(const string& input, bool double_back_slash,bool symmetry_string) const;
      string getPlotHeaderPDF(char function_mode,const string& column_header,bool points_color_gradient=DEFAULT_CHULL_LATEX_COLOR_GRADIENT) const;
      string getPlotPointContentPDF(const ChullPoint& point,bool zero_end_point=true,bool zero_dist_2_hull=false) const;
      string getNodeCoordPosition(const ChullPoint& point) const;
      string getNodeCoordPosition(const aflowlib::_aflowlib_entry& entry,const xvector<double>& coord) const;
      string nodeCreator(stringstream& option, stringstream& position, stringstream& content) const;
      string nodeCreator(const string& option, const string& position, const string& content) const;
      bool unwantedFacetLine(uint vi,uint vj,bool check_border=true) const;
      bool unwantedFacetLine(uint vi,uint vj,vector<vector<uint> >& facet_lines,bool check_border=true) const;
      string getPointsPropertyHeaderList(char mode) const;
      string getDelta(bool helvetica_font) const;
      string getSnapshotTableHeader(string headers,bool designate_HEADER=false) const;
      bool addInternalHyperlinks(bool internal_links_graph2report=true,bool internal_links_withinreport=true) const;
      bool addExternalHyperlinks() const;
      double getRoundToValue(double point_range) const;
      double getYTickDistance(double y_range,int approx_num_ticks,double round_to_value) const;
      vector<string> grabAcceptableLatexColors(bool replace_pranab_standard=true,bool allow_dvips_colors=true,uint count=10) const;
      vector<string> grabAcceptableLatexColors(const string& banned_colors_str,bool replace_pranab_standard=true,bool allow_dvips_colors=true,uint count=10) const;
      vector<string> grabAcceptableLatexColors(const vector<string>& banned_colors,bool replace_pranab_standard=true,bool allow_dvips_colors=true,uint count=10) const;
      aurostd::xoption resolvePlotLabelSettings() const;
      void writePDF() const;
      string getPlainTextHeader() const;
      string getJSONHeader() const;
      string grabCHPointProperty(const ChullPoint& point,const string& property,char mode=_txt_) const;
      string grabCHFacetProperty(const ChullFacet& facet,const string& property,char mode=_txt_) const;
      vector<vector<string> > getPointsData(const string& properties_str,vector<string>& headers,char mode=_txt_) const;
      vector<vector<vector<vector<string> > > > getFacetsData(const string& properties_str,vector<string>& headers,char mode=_txt_) const;
      void getPlainTextColumnSizes(const vector<string>& headers,const vector<vector<string> >& ventry,vector<uint>& sizes) const;
      void getPlainTextColumnSizesPoints(const vector<string>& headers,const vector<vector<string> >& ventries,vector<uint>& sizes) const;
      void getPlainTextColumnSizesFacets(const vector<string>& headers,const vector<vector<vector<vector<string> > > >& ventries,vector<uint>& sizes) const;
      string getPlainTextTable(const vector<string>& headers,const vector<vector<string> >& ventries,const vector<uint>& sizes) const;
      string getJSONTable(const vector<string>& headers,const vector<vector<string> >& ventries) const;
      void writeText(char mode=_txt_) const;
      void writeWebApp() const;
      void writeAPool() const;

      //sorting structures
      struct sortWithinCoordGroup {
        sortWithinCoordGroup(const vector<ChullPoint>& vpoints,bool sort_energy_ascending=true) : 
          m_points(vpoints),m_sort_energy_ascending(sort_energy_ascending) {};
        const vector<ChullPoint>& m_points;
        bool m_sort_energy_ascending; //good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
        bool operator() (uint ci,uint cj);
      };
      struct sortCHullPoints{
        sortCHullPoints(const vector<ChullPoint>& vpoints,bool sort_stoich_ascending=true,bool sort_energy_ascending=true) : 
          m_points(vpoints),m_sort_stoich_ascending(sort_stoich_ascending),m_sort_energy_ascending(sort_energy_ascending) {};
        const vector<ChullPoint>& m_points;
        bool m_sort_stoich_ascending; //good for sorting points in facets (lower vs. upper hemisphere), if facet in lower hemisphere, then sort ascending order (CLOCKWISE, graphing)
        bool m_sort_energy_ascending; //good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
        bool operator() (uint i,uint j) const;
      };
      struct sortFacetsByPoints {
        sortFacetsByPoints(const vector<ChullPoint>& vpoints,bool auto_sort_stoich=true,bool sort_stoich_ascending=true,bool auto_sort_energy=true,bool sort_energy_ascending=true) : 
          m_points(vpoints),m_auto_sort_stoich(auto_sort_stoich),m_sort_stoich_ascending(sort_stoich_ascending),
          m_auto_sort_energy(auto_sort_energy),m_sort_energy_ascending(sort_energy_ascending){};
        const vector<ChullPoint>& m_points;
        bool m_auto_sort_stoich;
        bool m_sort_stoich_ascending; //good for sorting points in facets (lower vs. upper hemisphere), if facet in lower hemisphere, then sort ascending order (CLOCKWISE, graphing)
        bool m_auto_sort_energy;
        bool m_sort_energy_ascending; //good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
        bool operator() (const ChullFacet& fi,const ChullFacet& fj) const;
      };
  };
} // namespace chull

#endif  // _AFLOW_CHULL_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *           Aflow COREY OSES - Duke University 2013-2018                  *
// *                                                                         *
// ***************************************************************************
