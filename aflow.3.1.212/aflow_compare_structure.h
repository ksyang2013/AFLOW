// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

#include<iostream>
#include "aflow_pflow.h"
#include "aflow.h"
#include "math.h"
#include<vector>
#include<string>
//#include <stdatomic.h> // Needed to communicate between threads
//#include<atomic>

using aurostd::isequal;
using aurostd::string2utype;
using aurostd::FileExist;

#ifndef __AFLOW_COMPARE_STRUCTURE_H_
#define __AFLOW_COMPARE_STRUCTURE_H_

#define JSON_MODE 0

// ===== GroupedWyckoffPosition Class ===== //
class GroupedWyckoffPosition{
  public:
    GroupedWyckoffPosition();
    ~GroupedWyckoffPosition();
    friend ostream& operator<<(ostream& oss, const GroupedWyckoffPosition& GroupedWyckoffPosition);
    const GroupedWyckoffPosition& operator=(const GroupedWyckoffPosition& b);
    GroupedWyckoffPosition(const GroupedWyckoffPosition& b);
    uint type;
    string element;
    vector<string> site_symmetries;
    vector<uint> multiplicities;
  private:
    void Free();
    void Copy(const GroupedWyckoffPosition& b);
};

// ===== StructurePrototype Class ===== //
class StructurePrototype{
  public:
    StructurePrototype();
    ~StructurePrototype();
    friend ostream& operator<<(ostream& oss, const StructurePrototype& StructurePrototype);
    const StructurePrototype& operator=(const StructurePrototype& b);
    StructurePrototype(const StructurePrototype& b);
    int iomode;
    string master_structure_name;
    xstructure master_structure;
    int number_types;
    vector<string> elements;
    vector<uint> stoichiometry;
    uint number_of_atoms;
    vector<string> unique_permutations;
    string pearson;
    uint space_group;
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    vector<string> wyckoff_site_symmetry;
    vector<int> wyckoff_multiplicity;
    vector<int> wyckoff_letter;
    vector<string> proto_structures_names;
    vector<xstructure> proto_structures;
    vector<string> family_structures_names;
    vector<xstructure> family_structures;
    vector<double> misfits;
    vector<double> family_misfits;
  private:
    void Free();
    void Copy(const StructurePrototype& b);
};



namespace compare{

  // ===== Main functions ===== //
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
                               const bool& same_species, const bool& scale_volume, const bool& fast_match, ostream& oss, double& final_misfit); //Main function
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species); //Overloaded, returns true (match), false (no match)
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species, const bool& scale_volume, const bool& fast_match); 
  double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species, const bool& fast_match); //Overloaded, returns misfit value
//  string CompareStructures(aurostd::xoption& vpflow);
//  string CompareStructureDirectory(aurostd::xoption& vpflow);

  // ===== Compare Directory Functions ===== //
  vector<uint> getStoichiometry(const xstructure& xstr, const bool& same_species);
  vector<string> getElements(xstructure& xstr);
  vector<uint> gcdStoich(const deque<int>& numbers);
  void calculateSymmetry(xstructure& xstr, vector<string>& vpearsons, vector<uint>& vsgroups,
                         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions);
  bool groupWyckoffPositions(xstructure& xstr, vector<GroupedWyckoffPosition>& grouped_positions);
  bool matchableWyckoffPositions(vector<GroupedWyckoffPosition>& temp_grouped_Wyckoffs,
                                 vector<GroupedWyckoffPosition>& master_grouped_Wyckoffs,
                                 const bool& same_species);

  void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes,
                                 const vector<xstructure>& vxstrs, const bool& same_species,
                                 const vector< vector<string> >& vvelements,
                                 vector< vector<uint> >& vstoichs, vector<string>& vpearsons,
                                 vector<uint>& vsgroups,
                                 vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
                                 const string& directory, const vector<string>& vfiles);

  vector<StructurePrototype> runComparisonScheme(uint& num_proc, vector<StructurePrototype>& comparison_schemes, const bool& same_species, 
                                                 const bool& scale_volume, const bool& fast_match, ostringstream& oss);
  vector<std::pair<uint,uint> > calculateDivisors(const int& number);
  bool checkNumberOfGroupings(vector<StructurePrototype>& comparison_schemes, uint number);
  int numberMismatches(const vector<StructurePrototype> comparison_schemes);
  void appendStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
                                 vector<StructurePrototype>& final_prototypes);
  void checkPrototypes(const uint& num_proc, const bool& same_species, vector<StructurePrototype>& final_prototypes);
  void prepareTextOutput(ostream& ss_out, const bool& same_species, const vector<StructurePrototype>& final_prototypes);

  bool SVD(xmatrix<double>& A);
  bool groupSameRatios(vector<int>& stoich, vector<int>& unique_stoich, vector<vector<int> >& type_index);
  //vector<vector<int> > generatePermutations(uint& num_elements, vector<int>& indices);




  // ===== All other functions ===== //
  bool matchableSpecies(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species);
  xstructure GetSuperCell3x3x3(const xstructure& aa, const xmatrix<double> &supercell);
  bool sameSpecies(const xstructure& x1, const xstructure& x2, const bool& display); 
  void rescaleStructure(xstructure& x1, xstructure& x2);
  void atomicNumberDensity(xstructure& xstr1, xstructure& xstr2);
  void fakeAtomsName(xstructure& x1);
  void printParameters(xstructure& xstr, ostream& oss);
  string leastFrequentAtom(const xstructure& xstr);
  vector<string> leastFrequentAtom2(const xstructure& xstr);
  bool checkTolerance(xvector<double> d1, xmatrix<double> Q2);
  bool checkABCTolerance(xvector<double> d1, xvector<double> d2);
  bool checkAngleTolerance(xvector<double> d1, xvector<double> d2);
  xvector<double> centroid_with_PBC(const xstructure& xstr);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, const xmatrix<double>& lattice);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, vector<double>& weights,
                                    const xmatrix<double>& lattice);
  bool findMatch(const xstructure& xstr1, const xstructure& PROTO,vector<uint>& im1, vector<uint>& im2, vector<double>& min_dists, const int& type_match);
  void clusterize(const xstructure& xstr1, const vector<uint>& im1, vector<string>& TYPE1,
                  vector<uint>& QTA1, vector<vector<uint> >& I1);
  bool sameAtomType(const xstructure& xstr1, const xstructure& xstr2, const vector<uint>& im1,
                    const vector<uint>& im2, const int& type_match);
  bool cleanMatch(const vector<uint>& im1);
  void cellDiagonal(xstructure& xstr, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  void cellDiagonal(xmatrix<double>& lattice, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  double latticeDeviation(const vector<double>& diag_sum1,const vector<double>& diag_sum2, 
                          const vector<double>& diag_diff1,const vector<double>& diag_diff2); 
  vector<double> computeNearestNeighbors(xstructure& xstr);
  double shortestDistance(const xstructure& xstr, const uint& k);
  void coordinateDeviation(const xstructure& xstr1, const xstructure& xstr2,
                      const vector<double>& all_nn1, const vector<double>& all_nn_proto,
                      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, vector<double>& min_dists,
                      double& cd, double& fail_figure);
  double computeMisfit(const double& dev, const double& dis, const double& fail);
  void printMatch(const vector<uint>& indexMatch1, const vector<uint>& indexMatch2,
                  const xstructure& PROTO, const xstructure& xstr1, ostream& oss);
  xvector<double> bringCoordinateInCell(xvector<double>& coord);
  double checkLatticeDeviation(double& xstr1_vol, xmatrix<double>& q2,vector<double>& D1,vector<double>& F1);
  bool quadrupletPeriodic(const xmatrix<double>& quad, const xstructure& lfa_supercell, const int& i, 
                       const int& j, const int& k, const int& w);
  bool vectorPeriodic(const xvector<double>& vec, const xstructure& lfa_supercell, const int& i, 
                       const int& j);
  void threadGeneration(const uint& num_proc,xmatrix<double>& q1, xstructure& xstr2, 
                        vector<xstructure> &vprotos, xstructure &xstr1, const int& type_match, 
                        const bool& fast_match, double& minMis, ostream& oss);
  void quadrupletSearch(const xmatrix<double>& q1, const xstructure& xstr_LFA_supercell,
                        vector<xvector<double> >& lattice_vecs, vector<vector<uint> >& ij_index);
//  void quadrupletSearch(std::atomic_bool &misfit_in_threshold_found, const string& lfa,
//                        const xmatrix<double>& q1,
//                        const xstructure& xstr_LFA_supercell,
//                        vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2,
//                        const int& type_match, double& possible_minMis, ostream& oss,
//                        vector<xvector<double> >& lattice_vecs,
//                        vector<vector<uint> >& ij_index); // 1 bracket
  bool buildSimilarLattices(vector<xvector<double> >& translation_vectors, xmatrix<double>& q1, double& xstr1_vol, double& abs_det_q1, xvector<double>& abc_angles_q1, vector<xmatrix<double> >& lattices, vector<xmatrix<double> >& clattices, vector<double>& latt_devs, const bool& fast_match);
  bool structureSearch(const string& lfa,
                        const vector<double>& all_nn1,
                        const xstructure& xstr,
                        vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2,
                        const int& type_match, double& possible_minMis,
                        vector<xmatrix<double> >& lattices,
                        vector<xmatrix<double> >& clattices,
                        vector<double>& latt_devs,
                        const bool& fast_match);

} // end namespace

#endif 

// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

