// ***************************************************************************
// 			    AFLOW Compare Structure 
//		David Hicks (d.hicks@duke.edu) and Carlo de Santo
// ***************************************************************************
#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include<exception>
#include<algorithm>
#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_compare_structure.h"
//#include "aflow_contrib_richard_wyckoff.h"

// ***************************************************************************
// AFLOW COMMANDS 
// ***************************************************************************

// ======== Single Comparison ======== //

  // aflow --compare_material=POSCAR_1,POSCAR_2
  //       Description: Compares POSCAR_1 to POSCAR_2 if the ratios are 
  //                    commensurate and types of atoms are the same 
  //                    (i.e. same material). 

  // aflow --compare_structure=POSCAR_1,POSCAR_2
  //       Description: Compares POSCAR_1 to POSCAR_2 (no requirement on type 
  //                    of atoms, just stoichiometry ratios). 

// ======== Directory Comparison ======== //

  // aflow --compare_material_directory|--compare_material_dir
  //       Description: Determines the unique materials (same atomic species) 
  //                     within a given directory. Returns a JSON and TXT file 
  //                     with the results. Default Directory: "."

  // aflow --compare_structure_directory|--compare_structure_dir
  //       Description: Determines the unique structure prototypes within a 
  //                    given directory. Returns a JSON and TXT file with the 
  //                    results. Default Directory: "."

  // OPTIONS:
  //       -D "PATH":   User can specify a specific directory to compare. 
  //                    Output will be placed there also.

// ***************************************************************************
// OPTIONS FOR ALL AFLOW COMMANDS
// ***************************************************************************
// --np=xx          : Number of processors. Algorithm is thread-friendly. 
//                    (Default: 8 processors)

// ***************************************************************************
// OVERVIEW OF ALGORITHM
// ***************************************************************************
// This algorithm takes two crystal structures and compares them to one another
// and determines their level of similarity as described by H. Burzlaff 
// (see his paper for more information: Acta Cryst., A53, 217-224 (1997)).

// Steps:
//   1) Scale volumes of two structures to be commensurate
//   2) Determine the least frequently occuring atom (LFA)
//   3) Shift the least frequently occuring atom (LFA) for 
//      both structures; these will be used to indicate our lattice
//   4) Create a 3x3x3 supercell of structure2
//   5) Search for all possible quadruplets of LFA atoms in the 
//      supercell and see if we can match it with structure1 
//      using Burzlaff's criteria
//   6) Once possible quadruplet/lattice is found, check contents
//      of lattice (i.e. atoms). Check we can have a one-to-one 
//      mapping
//   7) Of the best match (smallest misfit (mis)):
//      If mis <= 0.1:
//        Structures similar. Print out new representation of 
//        structure2 and the figure of misfit
//      else if 0.1 < mis <=0.2:  
//        Structures in the same family (possible symmetric
//        group-subgroup relation). Print out new representation 
//        of structure2 and the figure of misfit
//      else mis>0.2:
//        Structures are not the same. Print "No match"

// ***************************************************************************

// ***************************************************************************
// pflow::CompareStructures - Prepares comparison from command line input 
// ***************************************************************************
namespace pflow {
  string compareStructures(aurostd::xoption& vpflow){ 
    bool LDEBUG=(false || XHOST.DEBUG);
    ostringstream oss;
   
    string usage_material="aflow --compare_material=POSCAR1,POSCAR2";
    string usage_structure="aflow --compare_structure=POSCAR1,POSCAR2";
    //[THREADS]string options_single="--np=xx (default 8),--print";
    string options_single="--print"; //fast
 
    // ========== PARSE FLAGS ========== //
    // ===== FLAG: NUMBER OF PROCESSORS ===== //
    //[THREADS]int num_proc=8; //Defalut=8
    uint num_proc=1; //Defalut=1
    if(vpflow.flag("COMPARE_STRUCTURE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_STRUCTURE::NP"));
    }

    // ===== FLAG: TYPE OF COMPARISON ===== //
    bool same_species=false;
    // Material comparisons (find duplicates)
    if(vpflow.flag("COMPARE_MATERIAL")){
      same_species=true;
    }
    // Structure comparisons (find structure prototypes)
    else if (vpflow.flag("COMPARE_STRUCTURE")){ 
      same_species=false;
    }

    // ===== FLAG: STRUCTURES BEING COMPARED ===== //
    vector<string> vinput;
    if(vpflow.flag("COMPARE_STRUCTURE::STRUCTURES_1")) {
      vinput.push_back(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURES_1"));
    }
    if(vpflow.flag("COMPARE_STRUCTURE::STRUCTURES_2")) {
      vinput.push_back(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURES_2"));
    }
 
    // ===== FLAG: PRINT MAPPINGS ===== //
    bool print=false;
    if(vpflow.flag("COMPARE_STRUCTURE::PRINT")) {
      print=true;
    }
    
    // ===== FLAG: FAST MATCH ===== //
    bool fast_match=true; //false
    if(vpflow.flag("COMPARE_STRUCTURE::FAST")) {
      fast_match=true;
    }
    
    // ===== FLAG: SCALE VOLUME ===== //
    bool scale_volume=true;
    if(vpflow.flag("COMPARE_STRUCTURE::NO_SCALE_VOLUME")) {
      scale_volume=false;
    }

    // ========== Check structures ========== //
    double final_misfit=-1.0;
    if(LDEBUG) cerr << "pflow::compareStructures: begin" << endl;
    if(LDEBUG) cerr << "pflow::compareStructures: vinput.size()=" << vinput.size() << endl;
    // === Ensure 2 structures === //
    if(vinput.size()!=2) {      cerr << "pflow::compareStructures: ERROR vinput.size()!=2" << endl;
      if(vpflow.flag("COMPARE_MATERIAL")) {
        init::ErrorOption(cout,options_single,"pflow::compareStructures()",usage_material);
      }
      else if(vpflow.flag("COMPARE_STRUCTURE")) {
        init::ErrorOption(cout,options_single,"pflow::compareStructures()",usage_structure);
      }
      return oss.str();
    }

    if(LDEBUG) cerr << "pflow::compareStructures: loading vinput.at(0)=" << vinput.at(0) << endl;

    // === Check structure 1 === //
    if(!aurostd::FileExist(vinput.at(0))){
      oss << "pflow::compareStructures: ERROR file vinput.at(0)=" << vinput.at(0) << "  not found" << endl;
      return oss.str();
    }
    stringstream sss1;
    aurostd::efile2stringstream(vinput.at(0),sss1);
    xstructure xstr1(sss1);  

    // === Check structure 2 === //
    if(LDEBUG) cerr << "pflow::compareStructures: loading vinput.at(1)=" << vinput.at(1) << endl;
    if(!aurostd::FileExist(vinput.at(1))){
      cerr << "pflow::compareStructures: ERROR file vinput.at(1)=" << vinput.at(1) << "  not found" << endl;
      return oss.str();
    }
    stringstream sss2;
    aurostd::efile2stringstream(vinput.at(1),sss2);
    xstructure xstr2(sss2);  
  

    // ===== Call main comparison function ===== //
    compare::aflowCompareStructure(num_proc,xstr1,xstr2,same_species, scale_volume, fast_match, oss,final_misfit);
    if(print==true){
      return oss.str();
    }
    else{
      oss.str("");
      oss.clear();
      oss << final_misfit << " : ";
      if(final_misfit <=0.1 && (final_misfit+1.0)> 1e-3){
        oss << "MATCH" << endl;
      }
      else if(final_misfit > 0.1 && final_misfit <= 0.2){
        oss << "SAME FAMILY" << endl;
      }
      else if(final_misfit > 0.2 || (final_misfit+1.0) < 1e-3){ 
        oss << "NOT A MATCH" << endl;
      }
      return oss.str();
    }
  }
}


// ***************************************************************************
// pflow::COMPARE_STRUCTURE_DIRECTORY
// ***************************************************************************
namespace pflow {
  string compareStructureDirectory(aurostd::xoption& vpflow){

    string function_name = "pflow::compareStructureDirectory()";
    // This function takes compares all the structures within a given directory.
  
    bool LDEBUG=(false || XHOST.DEBUG);
    ostringstream oss;
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;
   
    // ========== PARSE FLAGS ========== //
    // FLAG: DIRECTORY
    string directory=vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY");
    if(!aurostd::FileExist(directory)) {
      message << "Unable to locate directory: " << directory << ". Exiting." << endl;
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
      //oss << "Unable to locate directory: " << directory << endl;
      //oss << "Exiting." << endl;
      return oss.str();
    }
    message << "Comparison directory: " << directory;
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // FLAG: NUMBER OF PROCESSORS
    //[THREADS]int num_proc=8; //Defalut=8
    uint num_proc=1; //Defalut=1
    if(vpflow.flag("COMPARE_STRUCTURE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_STRUCTURE::NP"));
      message << "OPTIONS: Using multiple threads; np = " << num_proc << ".";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
  
    // FLAG: TYPE OF COMPARISON
    bool same_species=false;
    if(vpflow.flag("COMPARE_MATERIAL_DIRECTORY")){  // Material comparisons (find duplicates)
      same_species=true;
      message << "OPTIONS: Performing material type comparisons (comparing alike atomic species).";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else if (vpflow.flag("COMPARE_STRUCTURE_DIRECTORY")){ // Structure comparisons (find structure prototypes)
      same_species=false;
      message << "OPTIONS: Performing structure type comparisons (any atomic species).";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
 
    
    // ===== FLAG: FAST MATCH ===== //
    bool fast_match=true; //false
    if(vpflow.flag("COMPARE_STRUCTURE::FAST")) {
      fast_match=true;
      message << "OPTIONS: Using fast match; if any match exists, return without exploring all options (a better match may exist)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ===== FLAG: SCALE VOLUME ===== //
    bool scale_volume=true;
    if(vpflow.flag("COMPARE_STRUCTURE::NO_SCALE_VOLUME")) {
      scale_volume=false;
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ========== Determine contents in directory ========== //
    vector<string> vfiles;
    vector<xstructure> vxstrs;
    vector< vector<uint> > vstoichs;
    vector< vector<string> > vvelements;
    aurostd::DirectoryLS(directory, vfiles);
    for(uint i=0; i<vfiles.size(); i++){
      //if(LDEBUG){cerr << "compare:: " << i << "/" << vfiles.size() << " " << vfiles[i] << endl;}
      // Convert files to xstructure objects
      // NOTE: Can only have geometry files in directory. Otherwise, xstructure will break and exits automatically.  
      // Only neglecting files which contain ".json" and ".out", since this routine writes these files and since
      // multiple runs may be done in a directory.
      if(vfiles[i].find("comparison_output.json") != std::string::npos || 
         vfiles[i].find("comparison_output.out") != std::string::npos || 
         vfiles[i].find("nohup.out") != std::string::npos){
        vfiles.erase(vfiles.begin()+i);
        i--;
      }
      else{
        stringstream sss1;
        aurostd::efile2stringstream(directory+"/"+vfiles[i],sss1);
        xstructure xstr1(sss1);
        vstoichs.push_back(compare::getStoichiometry(xstr1,same_species));
        vvelements.push_back(compare::getElements(xstr1));     
        vxstrs.push_back(xstr1);
        if(LDEBUG){
          cerr << "pflow::compareStructureDirectory() Found structure: " << directory+"/"+vfiles[i] << endl;
        }
      }
    }

    message << "Total number of structures to compare: " << vxstrs.size();
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ========== Calculate Symmetry of Structures ========== //
    message << "Calculating the symmetry of the structures.";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<string> vpearsons;
    vector<uint> vsgroups;
    vector<vector<GroupedWyckoffPosition> > vgrouped_Wyckoff_positions;
    for(uint i=0; i<vxstrs.size(); i++){
      compare::calculateSymmetry(vxstrs[i],vpearsons,vsgroups,vgrouped_Wyckoff_positions);
    }   
    message << "Symmetries calculated.";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
  
    // ========== Initialize StructurePrototype Objects ========== //
    vector< vector<vector<xstructure> > > vvvxstrs;
    vector< vector<uint> > vvunique_sgroups;
    vector<string> vunique_pearsons;
    vector<StructurePrototype> comparison_schemes;
  

    message << "Grouping sets of comparisons by stoichiometry and symmetry.";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // === Organize into objects based on stoichiometry and symmetry (Pearson and space group)
    compare::createStructurePrototypes(comparison_schemes, vxstrs, same_species, 
                                       vvelements, vstoichs, vpearsons, vsgroups, 
                                       vgrouped_Wyckoff_positions, directory, vfiles);
   

    message << "Running comparisons...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = compare::runComparisonScheme(num_proc, comparison_schemes, same_species, scale_volume, fast_match, oss); 
    
    if(final_prototypes.size()==0){
      return oss.str();
    }
    comparison_schemes.clear();
 
    message << "Number of unique prototypes: " << final_prototypes.size() << " (out of " << vxstrs.size() << " structures).";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

   
    // ========== Prepare JSON output ========== //
    stringstream ss_json;
    ss_json << "[" << endl;
    for(uint j=0; j<final_prototypes.size(); j++){
       ss_json << final_prototypes[j];
       if(j!=final_prototypes.size()-1){
         ss_json << "," << endl;
       }
    }
    ss_json << endl << "]" << endl;
  
    // ========== Prepare TEXT (.out) output ========== //
    stringstream ss_out;
    compare::prepareTextOutput(ss_out, same_species, final_prototypes);
  
    if(same_species==true){
      aurostd::stringstream2file(ss_json,directory+"/material_comparison_output.json");
      aurostd::stringstream2file(ss_out,directory+"/material_comparison_output.out");
      message << "RESULTS: See " << directory << "/material_comparison_output.out" << " or " << directory << "/material_comparison_output.json" << " for list of unique/duplicate materials.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else if(same_species==false){
      aurostd::stringstream2file(ss_json,directory+"/structure_comparison_output.json");
      aurostd::stringstream2file(ss_out,directory+"/structure_comparison_output.out");
      message << "RESULTS: See " << directory << "/structure_comparison_output.out" << " or " << directory << "/structure_comparison_output.json" << " for list of unique/duplicate structures.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
  
    return oss.str();
  
  
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species) {
    ostringstream oss;
    uint num_proc=1;
    double final_misfit=-1;
    bool scale_volume=true; //default is true
    bool fast_match=false; //default is false
    return aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, fast_match, oss, final_misfit);
  }
}

namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species, const bool& scale_volume, const bool& fast_match) {
    ostringstream oss;
    uint num_proc = 1;
    double final_misfit = -1;
    return aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, fast_match, oss, final_misfit);
  }
}


// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species) {
    ostringstream oss;
    uint num_proc=1;
    double final_misfit=-1;
    bool scale_volume=true; //default is true
    bool fast_match=false; //default is false
    aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, fast_match, oss, final_misfit);
    return final_misfit;
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
                             const bool &same_species, const bool& scale_volume, const bool& fast_match, 
                             ostream& oss, double& final_misfit) {

    // This is the main function. This compares two crystal structure together 
    // and determines their level of similarity based on the idea discussed 
    // in H. Burzlaff's paper (Acta Cryst., A53, 217-224 (1997)).

    bool LDEBUG=(false || XHOST.DEBUG);

    oss << "==================================================================================" << endl;

    xstructure xstr_base, xstr_test;
    // Need to expand larger system to ensure we find a commensurate unit cell between each of the structures
    if(xstr1.atoms.size()>xstr2.atoms.size()){
      xstr_base = xstr2;
      xstr_test = xstr1;
      if(LDEBUG){
        cerr << "compare::aflowCompareStructure: WARNING: Swapping order of xstructure 1 and 2 since, 1 is larger than the other." << endl;
      }
    }
    else{
      xstr_base = xstr1;
      xstr_test = xstr2;
    }

    xstr_base.FixLattices();
    xstr_test.FixLattices();
    xstr_base.ReScale(1.0);
    xstr_test.ReScale(1.0);

    // Below is not necessary, algorithm can handle supercells/conventional/prim
    //xstr_base=GetStandardPrimitive(xstr_base);
    //xstr_test=GetStandardPrimitive(xstr_test);
    //xstr_base.NiggliUnitCellForm();
    //xstr_test.NiggliUnitCellForm();

    int type_match=0;
    bool criteria_met = false;
    if(same_species == true){
      type_match=2;
      // If atoms in poscar not labeled in either POSCAR; assign fake names
      if(xstr_base.atoms.at(0).name == "" || xstr_test.atoms.at(0).name == ""){ 
        if(LDEBUG){cerr << "compare:: " << "Atoms not labeled ... Assigning fake names." << endl;}
        fakeAtomsName(xstr_base);
        fakeAtomsName(xstr_test);
      }
    }
    else if(same_species == false){
      type_match=1;
    }
    if(matchableSpecies(xstr_base,xstr_test,same_species)==true){
      criteria_met = true;
    }
    //cerr << "type_match: " << type_match << endl; 
    if(criteria_met == true){
      oss << "=========================================================" << endl; 

      oss << "STRUCTURE 1: " << endl;  
      oss << xstr_base << endl;
      //cerr << xstr_base << endl;

      oss << "=========================================================" << endl;

      oss << "STRUCTURE 2: " << endl;
      oss << xstr_test << endl;	
      //cerr << xstr_test << endl;

      oss << "=========================================================" << endl;

      //---------------------------------------------------------------------------//
      // DESCRIPTION OF OPTIONS

      // type_match determines how the atoms should be matched
      //  1: Assigns fake names to atoms (allows for structural comparison regardless of type of atom)
      //  2: Uses the names given in POSCAR (allows for structural comparison of material; type of atom necessary)

      if(same_species==true){
        type_match=2;
      }
      if(same_species==false){
        type_match=1;
      }	

      //====//VARIABLES
      uint i=0;
      xvector<double> origin;
      xstructure proto;           
      vector<xstructure> vprotos,vprotos_tmp;
      vector<vector<uint> > IM1, IM2;
      vector<vector<double> > vmin_dists;
      vector<uint> im1, im2;
      vector<string> PAIR1, PAIR2;
      double minMis=1;

      //====//PREPARATION
      oss<<"-------------------------------------------------------"<<endl;

      if(LDEBUG){cerr << "compare:: " << "Scale structures."<<endl;} 
      //THIS CONTAINS AN IF-STATEMENT...SHOULD I ALWAYS RESCALE THE STRUCTURES?
      rescaleStructure(xstr_base,xstr_test);

      if(scale_volume==true){atomicNumberDensity(xstr_base, xstr_test);}

      //----- Option: Fake Atoms name
      if(type_match==1){
        fakeAtomsName(xstr_base);
        fakeAtomsName(xstr_test);
      }

      // OBSOLETE THIS PRINTS OUT XSTRUCTURES WITH ATOM ZERO SHIFTED TO ORIGIN...
      // OBSOLETE oss<<"========================================================="<<endl;
      // OBSOLETE oss << xstr_base << endl;
      // OBSOLETE oss<<"========================================================="<<endl;
      // OBSOLETE oss << xstr_test << endl;		

      printParameters(xstr_base,oss);
      printParameters(xstr_test,oss);


      oss << "========================================================="<<endl;    
      oss << "QUADRUPLETS METHOD" << endl;

      xmatrix<double> q_base=xstr_base.lattice;; 
      
      if(LDEBUG){cerr << "compare:: " << "WAIT... Computing quadruplets..."<<endl;} 
      // ===== Creates the threads for checking quadruplets (lattices) ===== //
      threadGeneration(num_proc,q_base,xstr_test,vprotos,xstr_base,type_match,fast_match,minMis,oss);

      if(LDEBUG){cerr << "compare:: " << "Total # of possible matching representations: " << vprotos.size() << endl;}	

      //*************//======//FIND MATCH: Fake(1) or Pre-defined(2)

      // This first match finder is based on the best fitting between each atoms. This means that 
      // the routine looks for the closest atoms in order to find the match.
      // Can happen that one atom is the best matching for more than one atom in the second structure;
      // in this case the match is cancelled (cleanMatch)

      for(i=0; i<vprotos.size(); i++){
        //cerr << "xstr_base " << xstr_base << endl;
        //cerr << "vprotos[i] " << vprotos[i] << endl;
        //cerr << "orig: " << endl;
        //for(uint j=0;j<im1.size();j++){
        //  cerr << im1[j] << " == " << im2[j] << endl; 
        //}
        im1.clear(); im2.clear();
        //cerr << "after: " << endl;
        vector<double> min_dists;
        //findMatch(xstr_base,vprotos.at(i),im1,im2);
        //cerr << "find new match" << endl;
        findMatch(xstr_base,vprotos.at(i),im1,im2,min_dists,type_match);
        //cerr << "im1.size(): " << im1.size() << endl;
        //for(uint j=0;j<im1.size();j++){
        //  cerr << im1[j] << " == " << im2[j] << " (" << min_dists[j] << ")" << endl; 
        //}
        if(cleanMatch(im2)==false && cleanMatch(im1)==false){
          //cerr << "cleanMatch" << endl;
          vprotos_tmp.push_back(vprotos.at(i));
          IM1.push_back(im1);
          IM2.push_back(im2);
          vmin_dists.push_back(min_dists);
        }
      }
      if(LDEBUG){cerr << "compare:: " << "Number of matching representations: "<< vprotos_tmp.size() << endl;}

      vprotos.clear();
      vprotos=vprotos_tmp;
      vprotos_tmp.clear();

      if(vprotos.size()!=0){

        vector<vector<uint> > auxstr_base,auxstr_test;

        auxstr_base=IM1;
        auxstr_test=IM2;	
        IM1.clear();
        IM2.clear();

        // sameAtomType allow to check that the atoms matched are of the same type.
        // (can happen that certain transformation find matches between atoms of different type)

        for(i=0; i<vprotos.size(); i++){
          //cerr << "Same atom" << endl;
          //if(sameAtomType(xstr_base,vprotos.at(i),auxstr_base.at(i),auxstr_test.at(i),type_match)==true){
          //  cerr << "IN Same atom" << endl;
            vprotos_tmp.push_back(vprotos.at(i));
            IM1.push_back(auxstr_base.at(i));
            IM2.push_back(auxstr_test.at(i));
          //}
        }
        if(LDEBUG){cerr << "compare:: " << "Number of valid matches with the same type: " << vprotos_tmp.size() << endl;}

        vprotos.clear();
        vprotos=vprotos_tmp;
        vprotos_tmp.clear();

      }

      //---------------------

      if(vprotos.size()!=0){
        //======//MISFIT
        //cerr << "MIIIIIIIIIIIIIIIIIIIIIIINNNNNNN" << endl;
        //  This last part follows the computation of the figure of Misfit described by Burzlaff in his paper: 
        //  Burzlaff H., Malinovsky Y. (1996), "A Procedure for the Classification of Non-Organic Crystal structures."

        xstructure xstr_base_tmp = xstr_base;
        vector<double> diag_sum1,diag_sum2,diag_diff1,diag_diff2;
        double scale, lattdev;
        vector<double> vLattDevs;
        vector<double> vCoordDevs, vfails;
        double coorddev=1e9, fail_figure=1e9;
        double mis=1e9;
        vector<double> misfits;
        double min=1e9;
        xstructure better;

        for(i=0; i<vprotos.size(); i++){
          xstructure proto = vprotos[i];
          diag_sum1.clear();	diag_sum2.clear();
          diag_diff1.clear();	diag_diff2.clear();

          scale=xstr_base.Volume()/vprotos.at(i).Volume();
          scale=pow(scale,0.3333);

          cellDiagonal(xstr_base,diag_sum1,diag_diff1,1);
          cellDiagonal(vprotos.at(i),diag_sum2,diag_diff2,scale);

          lattdev=latticeDeviation(diag_sum1,diag_sum2,diag_diff1,diag_diff2);
          //cerr << "lattdev: " << lattdev << endl;
          vLattDevs.push_back(lattdev);

          vector<double> all_nn1 = computeNearestNeighbors(xstr_base); 
          vector<double> all_nn_proto = computeNearestNeighbors(proto);
          coordinateDeviation(xstr_base_tmp,proto,all_nn1,all_nn_proto,IM1.at(i),IM2.at(i),vmin_dists[i],coorddev,fail_figure);
          // DX ORIG coordinateDeviation(xstr_base,vprotos.at(i),IM1.at(i),IM2.at(i),coorddev,fail_figure);
          //cerr << "coorddev: " << coorddev << endl;
          //cerr << "fail_figure: " << fail_figure << endl;

          vCoordDevs.push_back(coorddev);
          vfails.push_back(fail_figure);

          mis=computeMisfit(lattdev,coorddev,fail_figure);
          min=mis;
	  if(LDEBUG){cerr << "compare:: " << "misfit: " << mis << "   (lattice deviation: " << lattdev << "  coordinate displacement: " 
                          << coorddev << "  figure of fail: " << fail_figure << ")" << endl;}
          misfits.push_back(mis);
          vprotos_tmp.push_back(vprotos.at(i));

          /*
          //If between 0.1 and 0.2 will try shifting method to see if we can obtain a figure of misfit under 0.1
          if(minMis <= 0.2 && minMis > 0.1){	
            if(fail_figure!=0 && std::isnan(misfits.at(i))==false){
              // This part is fundamental because it allows us to correct 
              // a trivial simplification done during the transformation: 
              // When the structure has been rotated, then it is shifted 
              // to each atom and brought in the cell to look for matchings. 
              // The atom shifted to the origin will coincide perfectly with 
              // the atom in the origin for the reference structure. 
              // However, this can lead to a matching failure between other 
              // pairs of atoms in a position different from the origin 
              // (see definition of failure on the paper). This routine 
              // aims to take the structure from the atom in the origin and 
              // move the entire structures little by little around this 
              // position. The rigid translation of all the atoms allow us 
              // to compute many different figures of misfit with the 
              // possibility that some failures disapper returning in the 
              // allowed tolerances.	
              if(LDEBUG){cerr << "compare:: " << "Attempting shift of structure since the minimum misfit is just above the similarity threshold..." << endl;}
              proto=vprotos.at(i);	
              for(j=0; j<proto.atoms.size(); j++){
                if(proto.atoms.at(j).cpos==origin){
                  delta=0.01*shortestDistance(proto,j);		
                  inc=0.2*delta;
                }
              }
              //cerr << "delta: " << delta << endl;
              for(double j=-delta; j<=delta; j=j+inc){
                //cerr << "j: " << j << endl;
                for(double k=-delta; k<=delta; k=k+inc){
                  for(double w=-delta; w<=delta; w=w+inc){
                    proto=vprotos.at(i);
                    for(uint iat=0; iat<proto.atoms.size(); iat++){
                      proto.atoms.at(iat).cpos(1)+=j;
                      proto.atoms.at(iat).cpos(2)+=k;
                      proto.atoms.at(iat).cpos(3)+=w;
                    }	
                    	
                    vector<double> min_dists; 
                    findMatch(xstr_base,proto,im1,im2,min_dists,type_match);
                    //coordinateDeviation(xstr_base,proto,IM1.at(i),IM2.at(i),coorddev,fail_figure);
                    coordinateDeviation(xstr_base,proto,im1,im2,coorddev,fail_figure);
                    //cerr << "lattdev: " << lattdev << endl;
                    //cerr << "coorddev: " << coorddev << endl;
                    mis_tmp=computeMisfit(lattdev,coorddev,fail_figure);
                    //cerr << "mis_tmp: " << mis_tmp << endl;
		    if(mis_tmp<min){
                      min=mis_tmp;
                      coorddev_tmp=coorddev;
                      fail_figure_tmp=fail_figure;
                      better=proto;
                    }
                  }
                }
              }		
              if(min<mis){
                vprotos_tmp.at(i)=better;
                vCoordDevs.at(i)=coorddev_tmp;
                vfails.at(i)=fail_figure_tmp;
                misfits.at(i)=min;
              }
            }
          }
          */
        }

        vprotos.clear();
        vprotos=vprotos_tmp;
        vprotos_tmp.clear();

        if(LDEBUG){cerr << "compare:: " << "Number of Misfits Computed:	"<<vprotos.size()<<endl;}

        // ====== RESULTS ====== //
        if(vprotos.size()!=0){

          uint min_index = 0; //DX 5/14/18 - added initialization
          int flag=0;

          for(i=0; i<misfits.size(); i++){
            if(flag==0){
              min=misfits.at(i);
              min_index=i;
              flag=1;
            }
            else{
              if(misfits.at(i)<min){
                min=misfits.at(i);
                min_index=i;
              }   
            }   
          }	   

          // ==== PRINT RESULTS ==== //
          oss << endl <<"**************************** RESULTS ****************************"<<endl;
          final_misfit=misfits.at(min_index);
          if(misfits.at(min_index)<0.1){
            oss << endl <<"MISFIT" <<":			" << misfits.at(min_index)<<"  STRUCTURES ARE COMPATIBLE" << endl;
            oss <<"----------------------------------------------------"<<endl;
            oss << "Figure of Deviation:	"<< vLattDevs.at(min_index) << endl;
            oss << "Figure of Displacement:	"<<vCoordDevs.at(min_index) << endl;
            oss << "Figure of Failure:	"<<vfails.at(min_index) << endl;
            oss <<"----------------------------------------------------"<<endl;
            printMatch(IM1.at(min_index),IM2.at(min_index),vprotos.at(min_index),xstr_base,oss);
            oss <<"----------------------------------------------------"<<endl;
            oss << "FINAL - REFERENCE STRUCTURE: " << endl;	
            oss << xstr_base << endl;
            oss <<"----------------------------------------------------"<<endl;
            oss << "FINAL - MAPPED STRUCTURE: " << endl;
            oss << vprotos.at(min_index);
          }
          else{
            if(misfits.at(min_index)<0.2){
              oss << endl <<"MISFIT" <<":                    " << misfits.at(min_index)<<"  STRUCTURES ARE IN THE SAME FAMILY" << endl;
            }
            else{
              oss << endl <<"MISFIT" <<":			" << misfits.at(min_index)<<"  STRUCTURES ARE INCOMPATIBLE (No match found)" << endl;
            }
            oss <<"----------------------------------------------------"<<endl;
            oss << "Figure of Deviation:	"<< vLattDevs.at(min_index) << endl;
            oss << "Figure of Displacement:	"<<vCoordDevs.at(min_index) << endl;
            oss << "Figure of Failure:	"<<vfails.at(min_index) << endl;
            oss <<"----------------------------------------------------"<<endl;
            printMatch(IM1.at(min_index),IM2.at(min_index),vprotos.at(min_index),xstr_base,oss);
            oss <<"----------------------------------------------------"<<endl;
            oss << "FINAL - REFERENCE STRUCTURE: " << endl;
            oss << xstr_base << endl;
            oss <<"----------------------------------------------------"<<endl;
            oss << "FINAL - MAPPED STRUCTURE: " << endl;
            oss << vprotos.at(min_index) << endl;
          } 

          //---------------------------------------------------------------------------//
        }
        else{
          oss << "[ERROR]: No match found!" << endl;
        }
      }
      else{ 
        oss << "[ERROR]: No match found!" << endl;
      }
      oss << endl << "*********************  THE END - FINE  **********************" << endl << endl;
      if(final_misfit<0.1 && !((final_misfit+1.0)<1e-3)){
        return true;
      }
      else{
        return false;
      }
    } //end of criteria_met
    return false;
  }
} //end of compare namespace

// ***************************************************************************
//				 END  -  FINE
// 			    AFLOW Compare Structure 
//		David Hicks (d.hicks@duke.edu) and Carlo de Santo
// ***************************************************************************
