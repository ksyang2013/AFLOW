// ***************************************************************************
// 			AFLOW Compare Structure - Functions
//		David Hicks (d.hicks@duke.edu) and Carlo de Santo
// ***************************************************************************
#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_compare_structure.h"
//#include <thread>
//#include <atomic>
#include "aflow_symmetry_spacegroup.h"
//#include <stdatomic.h> // Needed to communicate between threads

// ***************************************************************************
// Prototype Class 
// ***************************************************************************
// ===== Constructor ===== //
StructurePrototype::StructurePrototype(){ 
  iomode=JSON_MODE;
  master_structure_name="";
  master_structure.Clear();
  number_types=0;
  elements.clear();
  stoichiometry.clear();
  number_of_atoms=0;
  unique_permutations.clear();
  pearson="";
  space_group=0;
  grouped_Wyckoff_positions.clear();
  wyckoff_site_symmetry.clear();
  wyckoff_multiplicity.clear();
  wyckoff_letter.clear();
  proto_structures_names.clear();
  proto_structures.clear();
  family_structures_names.clear();
  family_structures.clear();
  misfits.clear();
  family_misfits.clear();
}

// ===== Free  ===== //
void StructurePrototype::Free(){
}

// ===== Destructor ===== //
StructurePrototype::~StructurePrototype(){ 
  master_structure.Clear();
  elements.clear();
  stoichiometry.clear();
  unique_permutations.clear();
  grouped_Wyckoff_positions.clear();
  wyckoff_site_symmetry.clear();
  wyckoff_multiplicity.clear();
  wyckoff_letter.clear();
  proto_structures_names.clear();
  proto_structures.clear();
  family_structures_names.clear();
  family_structures.clear();
  misfits.clear();
  family_misfits.clear();
  Free();
}

// ===== Copy Constructor ===== //
StructurePrototype::StructurePrototype(const StructurePrototype& b){
  Copy(b);
}

// ===== Copy Constructor Function ===== //
void StructurePrototype::Copy(const StructurePrototype& b) {
  if(this != &b){
    iomode=b.iomode;
    master_structure_name=b.master_structure_name; 
    master_structure=b.master_structure; 
    number_types=b.number_types;
    elements=b.elements;
    stoichiometry=b.stoichiometry;
    number_of_atoms=b.number_of_atoms;
    unique_permutations=b.unique_permutations;
    pearson=b.pearson;
    space_group=b.space_group;
    grouped_Wyckoff_positions=b.grouped_Wyckoff_positions;
    wyckoff_site_symmetry=b.wyckoff_site_symmetry;
    wyckoff_multiplicity=b.wyckoff_multiplicity;
    wyckoff_letter=b.wyckoff_letter;
    proto_structures_names=b.proto_structures_names;
    proto_structures=b.proto_structures;
    family_structures_names=b.family_structures_names;
    family_structures=b.family_structures;
    misfits=b.misfits;
    family_misfits=b.family_misfits;
  }
}

// ===== Assignment Operator (operator=) ===== //
const StructurePrototype& StructurePrototype::operator=(const StructurePrototype& b){
  if(this!=&b){
    Free();
    Copy(b);
  }
  return *this;
}

// ===== Output Handler ===== //
ostream& operator<<(ostream& oss, const StructurePrototype& StructurePrototype){
  //if(StructurePrototype.iomode!=JSON_MODE){ //A safeguard until we construct more output schemes.
  //  StructurePrototype.iomode=JSON_MODE;
  //}
  if(StructurePrototype.iomode==JSON_MODE){
    oss << "{" << endl;
    oss << "  \"master_structure\": \"" << StructurePrototype.master_structure_name << "\"," << endl;
    oss << "  \"number_types\": \"" << StructurePrototype.number_types << "\"," << endl;
    oss << "  \"elements\": [";
    for(uint i=0; i<StructurePrototype.elements.size(); i++){
      if(i==0){
        oss << endl << "    \"" << StructurePrototype.elements[i] << "\"";
      }
      else{
        oss << "," << endl << "    \"" << StructurePrototype.elements[i] << "\"";
      }
    }
    oss << endl << "  ]," << endl;
    oss << "  \"stoichiometry\": \"";
    for(uint i=0; i<StructurePrototype.stoichiometry.size(); i++){
      if(i==0){
        oss << StructurePrototype.stoichiometry[i];
      }
      else{
        oss << ":" << StructurePrototype.stoichiometry[i];
      }
    }
    oss << "\"," << endl;
    oss << "  \"number_of_atoms\": \"" << StructurePrototype.number_of_atoms << "\"," << endl;
    //oss << "  \"unique_permutations\": [";
    //for(uint i=0; i<StructurePrototype.unique_permutations.size(); i++){
    //  if(i==0){
    //    oss << endl << "    \"" << StructurePrototype.unique_permutations[i] << "\"";
    //  }
    //  else{
    //    oss << "," << endl << "    \"" << StructurePrototype.unique_permutations[i] << "\"";
    //  }
    //}
    //if(StructurePrototype.unique_permutations.size() != 0){  
    //  oss << endl << "  ]," << endl;
    //}
    //else{
    //  oss << "]," << endl;
    //}
    //oss << "  \"pearson\": \"" << StructurePrototype.pearson << "\"," << endl;
    oss << "  \"space_group\": \"" << StructurePrototype.space_group << "\"," << endl;
    oss << "  \"grouped_Wyckoff_positions\": ["; // << endl;
    for(uint i=0; i<StructurePrototype.grouped_Wyckoff_positions.size(); i++){
      if(i==0){
        oss << endl << "    " << StructurePrototype.grouped_Wyckoff_positions[i];
      }
      else{
        oss << "," << endl << "    " << StructurePrototype.grouped_Wyckoff_positions[i];
      }
    }
    if(StructurePrototype.grouped_Wyckoff_positions.size() != 0){  
      oss << "  ]," << endl;
    }
    else{
      oss << "]," << endl;
    }
    //oss << "  \"wyckoff_positions\": ["; // << endl;
    //for(uint i=0; i<StructurePrototype.wyckoff_site_symmetry.size(); i++){
    //  if(i==0){
    //    oss << StructurePrototype.wyckoff_multiplicity[i] << " " << StructurePrototype.wyckoff_site_symmetry[i] <<  StructurePrototype.wyckoff_letter[i];
    //  }
    //  else{
    //    oss << "; " << StructurePrototype.wyckoff_multiplicity[i] << " " << StructurePrototype.wyckoff_site_symmetry[i] <<  StructurePrototype.wyckoff_letter[i];
    //  }
    //}
    //if(StructurePrototype.wyckoff_site_symmetry.size() != 0){  
    //  oss << "  ]," << endl;
    //}
    //else{
    //  oss << "]," << endl;
    //}
    oss << "  \"proto_structures\": ["; // << endl;
    for(uint i=0; i<StructurePrototype.proto_structures_names.size(); i++){
      if(i==0){
        oss << endl << "    \"" << StructurePrototype.proto_structures_names[i] << "\"";
      }
      else{
        oss << "," << endl << "    \"" << StructurePrototype.proto_structures_names[i] << "\"";
      }
    }
    if(StructurePrototype.proto_structures_names.size() != 0){  
      oss << endl << "  ]," << endl;
    }
    else{
      oss << "]," << endl;
    }
    oss << "  \"misfits\": [";// << endl;
    for(uint i=0; i<StructurePrototype.misfits.size(); i++){
      if(i==0){
        oss << endl << "    \"" << StructurePrototype.misfits[i] << "\"";
      }
      else{
        oss << "," << endl << "    \"" << StructurePrototype.misfits[i] << "\"";
      }
    }
    if(StructurePrototype.misfits.size() != 0){
      oss << endl << "  ]" << endl;
    }
    else{
      oss << "]" << endl;
    }
    //oss << "  \"family_structures\": ["; // << endl;
    //for(uint i=0; i<StructurePrototype.family_structures_names.size(); i++){
    //  if(i==0){
    //    oss << endl << "    \"" << StructurePrototype.family_structures_names[i] << "\"";
    //  }
    //  else{
    //    oss << "," << endl << "    \"" << StructurePrototype.family_structures_names[i] << "\"";
    //  }
    //}
    //if(StructurePrototype.family_structures_names.size() != 0){  
    //  oss << endl << "  ]," << endl;
    //}
    //else{
    //  oss << "]," << endl;
    //}
    //oss << "  \"family_misfits\": [";// << endl;
    //for(uint i=0; i<StructurePrototype.family_misfits.size(); i++){
    //  if(i==0){
    //    oss << endl << "    \"" << StructurePrototype.family_misfits[i] << "\"";
    //  }
    //  else{
    //    oss << "," << endl << "    \"" << StructurePrototype.family_misfits[i] << "\"";
    //  }
    //}
    //if(StructurePrototype.family_misfits.size() != 0){
    //  oss << endl << "  ]" << endl;
    //}
    //else{
    //  oss << "]" << endl;
    //}
    oss << "}";
  }
  return oss;
}

// ***************************************************************************
// END:: Prototype Class
// ***************************************************************************

// ***************************************************************************
// START:: GroupedWyckoffPosition Class
// ***************************************************************************
// ===== Constructor ===== //
GroupedWyckoffPosition::GroupedWyckoffPosition(){
  type = 0;
  element = "";
  site_symmetries.clear();
  multiplicities.clear();
}

// ===== Free  ===== //
void GroupedWyckoffPosition::Free(){
}

// ===== Destructor  ===== //
GroupedWyckoffPosition::~GroupedWyckoffPosition(){ 
  site_symmetries.clear();
  multiplicities.clear();
}

// ===== Copy Constructor ===== //
GroupedWyckoffPosition::GroupedWyckoffPosition(const GroupedWyckoffPosition& b){
  Copy(b);
}

// ===== Copy Constructor Function ===== //
void GroupedWyckoffPosition::Copy(const GroupedWyckoffPosition& b) {
  if(this != &b){
    type=b.type;
    element=b.element;
    site_symmetries=b.site_symmetries;
    multiplicities=b.multiplicities;
  }
}

// ===== Assignment Operator (operator=) ===== //
const GroupedWyckoffPosition& GroupedWyckoffPosition::operator=(const GroupedWyckoffPosition& b){
  if(this!=&b){
    Free();
    Copy(b);
  }
  return *this;
}

// ===== Output Handler ===== //
ostream& operator<<(ostream& oss, const GroupedWyckoffPosition& GroupedWyckoffPosition){
  //if(StructurePrototype.iomode!=JSON_MODE){ //A safeguard until we construct more output schemes.
  //  StructurePrototype.iomode=JSON_MODE;
  //}
  oss << "{" << endl;
  //oss << "  \"type\": \"" << GroupedWyckoffPosition.type << "\"," << endl;
  oss << "  \"element\": \"" << GroupedWyckoffPosition.element << "\"," << endl;
  oss << "  \"site_symmetries\": [";
  for(uint i=0; i<GroupedWyckoffPosition.site_symmetries.size(); i++){
    if(i==0){
      oss << endl << "    \"" << GroupedWyckoffPosition.site_symmetries[i] << "\"";
    }
    else{
      oss << "," << endl << "    \"" << GroupedWyckoffPosition.site_symmetries[i] << "\"";
    }
  }
  if(GroupedWyckoffPosition.site_symmetries.size() != 0){  
    oss << endl << "  ]," << endl;
  }
  else{
    oss << "]," << endl;
  }
  oss << "  \"multiplicities\": [";
  for(uint i=0; i<GroupedWyckoffPosition.multiplicities.size(); i++){
    if(i==0){
      oss << endl << "    " << GroupedWyckoffPosition.multiplicities[i] << "";
    }
    else{
      oss << "," << endl << "    " << GroupedWyckoffPosition.multiplicities[i] << "";
    }
  }
  if(GroupedWyckoffPosition.multiplicities.size() != 0){  
    oss << endl << "  ]" << endl;
  }
  else{
    oss << "]" << endl;
  }
  oss << "}" << endl;
  return oss;
}

// ***************************************************************************
// END:: GroupedWyckoffPosition Class
// ***************************************************************************

// ***************************************************************************


// ***************************************************************************
// *                                                                         *
// *                             FUNCTIONS                                   *
// *                                                                         *
// ***************************************************************************

// ***************************************************************************
// groupSameRatios
// ***************************************************************************
namespace compare{
  bool groupSameRatios(vector<int>& stoich, vector<int>& unique_stoich, vector<vector<int> >& type_index){
    for(uint i=0;i<stoich.size();i++){
      bool stoich_stored = false;      
      for(uint j=0;j<unique_stoich.size();j++){     
        if(stoich[i] == unique_stoich[j]){
          stoich_stored = true;
          type_index[j].push_back(i);
          break;
        }
      }
      if(!stoich_stored){
        unique_stoich.push_back(stoich[i]);
        vector<int> tmp; tmp.push_back(i);
        type_index.push_back(tmp);
      }
    }
    return true;
  }
}


// ***************************************************************************
//  Get Stoichiometry - Obtain stoichiometries from xstructure
// ***************************************************************************
namespace compare{
  vector<uint> getStoichiometry(const xstructure& xstr, const bool& same_species){

    // Obtains the least common multiple representation of the stoichiometry.

    vector<uint> stoich;
    if(xstr.species.size()==1){
       stoich.push_back(1);
    }
    else{
      stoich=gcdStoich(xstr.num_each_type);
    }
    // If a structure prototype comparison (not material type), ensure 
    // stoichiometries are in numerical order for comparison.  Else, 
    // leave in position indicating atomic species count.
    if(same_species==false){
      for(uint i=0; i<stoich.size(); i++){
	std::sort(stoich.begin(),stoich.end());
      }
    }
    return stoich;
  }
}

// ***************************************************************************
//  Get Elements 
// ***************************************************************************
namespace compare{
  vector<string> getElements(xstructure& xstr){
    // Obtains the elements in the structure.

    bool LDEBUG=(false || XHOST.DEBUG);
    vector<string> velements;
    // If atoms in poscar not labeled in either POSCAR; assign fake names
    if (xstr.atoms[0].name == ""){
      if(LDEBUG){cerr << "compare:: " << "WARNING!!!!! Atoms not labeled ... Assigning Fake names" << endl;}
      fakeAtomsName(xstr);
    }
    string prev_element="";
    for(uint i=0; i<xstr.atoms.size(); i++){
      if(xstr.atoms[i].name != prev_element){
	velements.push_back(xstr.atoms[i].name);
	prev_element=xstr.atoms[i].name;
      }
    }
    return velements;
  }
}

// ***************************************************************************
// gcdStoich - Euler's Greatest Common Divisor Algorithm
// ***************************************************************************
namespace compare{
  vector<uint> gcdStoich(const deque<int>& numbers){

    // This is Euler's Greates Common Divisor Algorithm.  It is used to determine 
    // the least common multiple representation for the stoichiometry.

    int global_GCD = 0; //DX 5/14/18 - added initialization
    int GCD = 0; //DX 5/14/18 - added initialization
    vector<uint> reduced_numbers;
    // Find min number first
    int min=0;
    for(uint i=0; i<numbers.size(); i++){
      if(i==0){
	min=numbers[i];
      }
      else{
	if(numbers[i]<min){
	   min=numbers[i];
	}
      }
    }
    bool found_GCD=true;
    for(uint i=0; i<numbers.size(); i++){
      if(numbers[i]%min != 0){
	found_GCD=false;
	break;
      }
    }
    if(found_GCD==true){
      global_GCD=min;
    }
    else if(found_GCD==false){
      int remainder=1000;
      int divisor=min;
      for(uint i=0; i<numbers.size(); i++){
	int num=numbers[i];
	while(remainder != 0){
	  remainder=(num%divisor);
	  if(remainder==0){
	    GCD=divisor;
	  }
	  else{
	    num=divisor;
	    divisor=remainder;
	  }
	}
	divisor=GCD;
	if(i==0){
	  global_GCD=GCD;
	}
	else if(GCD < global_GCD){
	   global_GCD=GCD;
	}
	remainder=100;
      }
    }
    for(uint i=0; i<numbers.size(); i++){
       reduced_numbers.push_back((uint)(numbers[i]/global_GCD));
       if(numbers[i]%global_GCD){
	  cerr << "compare::ERROR - Error in GCD procedure. Exiting. (David Hicks: d.hicks@duke.edu)" << endl;
	  exit(1);
       }
    }
    return reduced_numbers;
  }
}

// ***************************************************************************
// calculateSymmetry - Calculate Symmetry
// ***************************************************************************
namespace compare{
  void calculateSymmetry(xstructure& xstr, vector<string>& vpearsons, vector<uint>& vsgroups, 
                         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions){

    // Calculates the Pearson symbol and space group of the structure

    //xstr.GetLatticeType(); //A little slow. Consider a different method -> (SpaceGroup_ITC) finds this
    //vpearsons.push_back(xstr.pearson_symbol);
    vpearsons.push_back("xX");
    vsgroups.push_back(xstr.SpaceGroup_ITC());
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions; 
    groupWyckoffPositions(xstr, grouped_Wyckoff_positions);
    vgrouped_Wyckoff_positions.push_back(grouped_Wyckoff_positions);
  }
}

// ***************************************************************************
// groupWyckoffPositions
// ***************************************************************************
namespace compare{
  bool groupWyckoffPositions(xstructure& xstr, vector<GroupedWyckoffPosition>& grouped_positions){
    for(uint i=0;i<xstr.wyckoff_sites_ITC.size();i++){
      vector<string> tokens;
      aurostd::string2tokens(xstr.wyckoff_sites_ITC[i].wyckoffSymbol,tokens," ");  
      uint multiplicity = aurostd::string2utype<uint>(tokens[0]);     
      string letter = aurostd::string2utype<string>(tokens[1]);     
      string site_symmetry = aurostd::string2utype<string>(tokens[2]);     

      bool element_found = false;
      uint element_index = 0;
      for(uint j=0;j<grouped_positions.size();j++){
        if(xstr.wyckoff_sites_ITC[i].type == grouped_positions[j].element){
          element_found = true;
          element_index = j;
          break;
        }
      }
      if(element_found == false){
        GroupedWyckoffPosition tmp;
        tmp.element = xstr.wyckoff_sites_ITC[i].type;   
        tmp.site_symmetries.push_back(site_symmetry);          
        tmp.multiplicities.push_back(multiplicity);          
        grouped_positions.push_back(tmp);
      }
      else{
        grouped_positions[element_index].site_symmetries.push_back(site_symmetry);
        grouped_positions[element_index].multiplicities.push_back(multiplicity);
      }
    }
    //cerr << "xstr: " << xstr << endl;
    //for(uint j=0;j<grouped_positions.size();j++){
    //  cerr << "grouped wyckoffs: " << grouped_positions[j] << endl;
    //}
    return true;
  }
}

// ***************************************************************************
// matchableWyckoffPositions
// ***************************************************************************
namespace compare{
  bool matchableWyckoffPositions(vector<GroupedWyckoffPosition>& temp_grouped_Wyckoffs,
                                 vector<GroupedWyckoffPosition>& master_grouped_Wyckoffs, 
                                 const bool& same_species){
    
    // quick check: are the number of Wyckoff positions the same; cannot match otherwise
    if(temp_grouped_Wyckoffs.size() != master_grouped_Wyckoffs.size()){
      return false;
    }

    vector<vector<bool> > found_matches;
    for(uint i=0;i<temp_grouped_Wyckoffs.size();i++){
      vector<bool> tmp;
      for(uint m=0;m<temp_grouped_Wyckoffs[i].multiplicities.size();m++){
        tmp.push_back(false);
      }
      found_matches.push_back(tmp);
    }

    for(uint i=0;i<temp_grouped_Wyckoffs.size();i++){
      for(uint j=0;j<master_grouped_Wyckoffs.size();j++){
        if(same_species && temp_grouped_Wyckoffs[i].element == master_grouped_Wyckoffs[j].element &&
           temp_grouped_Wyckoffs[i].multiplicities.size() == master_grouped_Wyckoffs[j].multiplicities.size()){
          uint match_counts = 0;
          for(uint m=0;m<temp_grouped_Wyckoffs[i].multiplicities.size();m++){
            for(uint n=0;n<master_grouped_Wyckoffs[j].multiplicities.size();n++){
              if(temp_grouped_Wyckoffs[i].multiplicities[m] == master_grouped_Wyckoffs[j].multiplicities[n] &&
                 temp_grouped_Wyckoffs[i].site_symmetries[m] == master_grouped_Wyckoffs[j].site_symmetries[n]){
                found_matches[i][m] = true;
                match_counts++;
              }
            }
          }
        }
        else if(!same_species && temp_grouped_Wyckoffs[i].multiplicities.size() == master_grouped_Wyckoffs[j].multiplicities.size()){
          uint match_counts = 0;
          for(uint m=0;m<temp_grouped_Wyckoffs[i].multiplicities.size();m++){
            for(uint n=0;n<master_grouped_Wyckoffs[j].multiplicities.size();n++){
              if(temp_grouped_Wyckoffs[i].multiplicities[m] == master_grouped_Wyckoffs[j].multiplicities[n] &&
                 temp_grouped_Wyckoffs[i].site_symmetries[m] == master_grouped_Wyckoffs[j].site_symmetries[n] &&
                 !found_matches[i][m]){ // and not already matched
                found_matches[i][m] = true;
                match_counts++;
              }
            }
          } 
          // if any match, all need to match; otherwise the Wyckoff positions are not matchable
          if(match_counts>0 && match_counts != temp_grouped_Wyckoffs[i].multiplicities.size()){ 
            vector<bool> tmp; for(uint m=0;m<temp_grouped_Wyckoffs[i].multiplicities.size();m++){tmp.push_back(false);}
            found_matches[i] = tmp;
          }
        }
      }
    }

    for(uint i=0;i<found_matches.size();i++){
      for(uint j=0;j<found_matches[i].size();j++){
        if(found_matches[i][j] == false){
          //cerr << "could not match!!!: " << i << " " << j << endl;
          return false;
        }
      }
    }
    return true;
  }
}

// ***************************************************************************
// createStructurePrototypes - Group structures by Pearson symbol, then space group
// ***************************************************************************
namespace compare{
  void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
				 const vector<xstructure>& vxstrs, const bool& same_species, 
				 const vector< vector<string> >& vvelements,
				 vector< vector<uint> >& vstoichs, vector<string>& vpearsons, 
				 vector<uint>& vsgroups, 
                                 vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
                                 const string& directory, const vector<string>& vfiles){

    // Populates the structure information into the StructurePrototype object.
    // It groups structure based on their stoichiometry and pearson symbol and 
    // space group. A "master" structure is chosen and will be compared to the 
    // possible "duplicates". The misfit values are set to -1.0 until compared.

    // First, separate by stoichiometry
    for(uint i=0;i<vstoichs.size(); i++){
      bool scheme_created=false;
      if(i==0){
        StructurePrototype tmp;
        tmp.master_structure_name=directory+"/"+vfiles[i];
        tmp.number_types=vxstrs[i].num_each_type.size();
        tmp.elements=vvelements[i];
        tmp.stoichiometry=vstoichs[i];
        tmp.number_of_atoms=vxstrs[i].atoms.size();
        tmp.pearson=vpearsons[i];
        tmp.space_group=vsgroups[i];
        tmp.grouped_Wyckoff_positions=vgrouped_Wyckoff_positions[i];
        tmp.master_structure=vxstrs[i];
        comparison_schemes.push_back(tmp);
      }
      else{
        for(uint j=0; j<comparison_schemes.size(); j++){
          bool same_material_stoich=false;
          ostringstream tmp;
          tmp.clear();
          if(same_species==true && 
             matchableSpecies(vxstrs[i],comparison_schemes[j].master_structure,same_species)==true){
            same_material_stoich=true;
          }
          else if(same_species==false){
            same_material_stoich=true;
          }
          if(same_material_stoich==true && 
             vstoichs[i] == comparison_schemes[j].stoichiometry && 
             vpearsons[i] == comparison_schemes[j].pearson && 
             vsgroups[i] == comparison_schemes[j].space_group &&
             matchableWyckoffPositions(vgrouped_Wyckoff_positions[i], comparison_schemes[j].grouped_Wyckoff_positions,same_species)){
            if(same_species==false){
              for(uint e=0;e<vvelements[i].size();e++){
                bool already_in=false;
                for(uint f=0;f<comparison_schemes[j].elements.size();f++){
                  if(vvelements[i][e]==comparison_schemes[j].elements[f]){
                    already_in=true;
                    break;
                  }
                }
                if(already_in==false){
                  comparison_schemes[j].elements.push_back(vvelements[i][e]);
                }
              }
            }
            comparison_schemes[j].proto_structures_names.push_back(directory+"/"+vfiles[i]);
            comparison_schemes[j].proto_structures.push_back(vxstrs[i]);
            comparison_schemes[j].misfits.push_back(-1.0);
            scheme_created=true;
            break;
          }
        }
        if(scheme_created==false){
          StructurePrototype tmp;
          tmp.master_structure_name=directory+"/"+vfiles[i];
          tmp.number_types=vxstrs[i].num_each_type.size();
          tmp.elements=vvelements[i];
          tmp.stoichiometry=vstoichs[i];
          tmp.number_of_atoms=vxstrs[i].atoms.size();
          tmp.pearson=vpearsons[i];
          tmp.space_group=vsgroups[i];
          tmp.grouped_Wyckoff_positions=vgrouped_Wyckoff_positions[i];
          tmp.master_structure=vxstrs[i];
          comparison_schemes.push_back(tmp);
        }
      }
    }
  }
}


// ***************************************************************************
// runComparisonScheme: Runs comparisons automatically 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> runComparisonScheme(uint& num_proc, vector<StructurePrototype>& comparison_schemes, const bool& same_species, const bool& scale_volume, const bool& fast_match, ostringstream& oss){

    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;
    string function_name = "compare::runComparisonScheme()";

    // ========== Initial Comparisons ========== //
    // Compare within each comparison_scheme object
    for(uint i=0; i<comparison_schemes.size(); i++){
      for(uint j=0; j<comparison_schemes[i].proto_structures.size(); j++){
        ostringstream tmp_oss;
        tmp_oss.clear();
        double final_misfit=-1.0;
        //bool scale_volume=true; //default is true
        //bool fast_match=false; //default is false
        // Call the main comparison function
        compare::aflowCompareStructure(num_proc, comparison_schemes[i].master_structure,
                                       comparison_schemes[i].proto_structures[j],
                                       same_species, scale_volume, fast_match, tmp_oss, final_misfit);
        // Store the figure of misfit
        comparison_schemes[i].misfits[j]=final_misfit;
      }
    }

    // === Count the number of mismatches (i.e. mis > 0.1)
    int num_mismatches_orig=compare::numberMismatches(comparison_schemes);
    int num_mismatches=num_mismatches_orig;
    vector<StructurePrototype> final_prototypes;

    if(num_mismatches==0){
      compare::appendStructurePrototypes(comparison_schemes, final_prototypes);
    }
    
    // ========== Continue comparison until all strucutures are matched or all comparisons schemes exhaused ========== //
    while(num_mismatches!=0){
      // Clean comparison_schemes object; store final prototypes
      compare::appendStructurePrototypes(comparison_schemes, final_prototypes);

      for(uint i=0; i<comparison_schemes.size(); i++){
        for(uint j=0; j<comparison_schemes[i].proto_structures.size(); j++){
          ostringstream tmp_oss;
          tmp_oss.clear();
          double final_misfit=-1.0;
          compare::aflowCompareStructure(num_proc, comparison_schemes[i].master_structure,
                                         comparison_schemes[i].proto_structures[j],
                                         same_species, scale_volume, fast_match, tmp_oss, final_misfit);
          comparison_schemes[i].misfits[j]=final_misfit;
        }
      }
      // Update number of mismatches
      num_mismatches_orig=num_mismatches;
      num_mismatches=compare::numberMismatches(comparison_schemes);

      if(num_mismatches > 0){
        message << "Number of unmatched structures: " << num_mismatches << ". Continuing comparisons ...";
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      }


      // Ensure while loop is not uncontrolled
      if(num_mismatches>=num_mismatches_orig){
         oss << "compare::ERROR - Number of mismatches is increasing...impossible (comparision framework flawed). "
              << "Contact David Hicks (d.hicks@duke.edu)." << endl;
         final_prototypes.clear();
         return final_prototypes;
      }
    }

    // ========== Append new prototypes ========== //
    final_prototypes.insert(final_prototypes.end(),comparison_schemes.begin(),comparison_schemes.end());
    return final_prototypes;
  }
}

// ***************************************************************************
// calculateDivisors 
// ***************************************************************************
namespace compare{
  vector<std::pair<uint,uint> > calculateDivisors(const uint& number){
    vector<std::pair<uint,uint> > divisor_pairs;
    std::pair<uint,uint> tmp_pair; 
    for(uint i=1; i<=number; i++){
      if(number%i==0){
        tmp_pair.first = i; tmp_pair.second = number/i;
        divisor_pairs.push_back(tmp_pair);
      }
    }
    return divisor_pairs;   
  }
}

// ***************************************************************************
// calculateDivisors 
// ***************************************************************************
namespace compare{
  bool checkNumberOfGroupings(vector<StructurePrototype>& comparison_schemes, uint number){
    vector<std::pair<uint,uint> > divisor_pairs = calculateDivisors(number);
    uint divisor_index = 0;
    for(uint i=0;i<divisor_pairs.size();i++){
      //cerr << "divisor_pairs: " << divisor_pairs[i].first << ", " << divisor_pairs[i].second << endl;
      //cerr << "comparison_schemes.size(): " << comparison_schemes.size() << endl;
      if(comparison_schemes.size() == divisor_pairs[i].first){
        divisor_index = i;
        break;
      }
    }
    uint num_consistent = 0;
    for(uint j=0; j<comparison_schemes.size(); j++){
      //cerr << "comparisons_schemes[j]: " << comparison_schemes[j] << endl;
      //cerr << comparison_schemes[j].proto_structures_names.size()+1 << " vs " << divisor_pairs[divisor_index].second << endl;
      if(comparison_schemes[j].proto_structures_names.size()+1 == divisor_pairs[divisor_index].second){ //+1 to include master
        num_consistent++;
      }
    }
    if(num_consistent != comparison_schemes.size()){
      return false;
    }
    return true;
  }
}



// ***************************************************************************
// NumberMismaches - Count the number of non-matches
// ***************************************************************************
namespace compare{
  int numberMismatches(const vector<StructurePrototype> comparison_schemes){
    int num_mismatches=0;
    for(uint i=0; i<comparison_schemes.size(); i++){
      for(uint j=0; j<comparison_schemes[i].misfits.size(); j++){
        if(comparison_schemes[i].misfits[j] > 0.1 || comparison_schemes[i].misfits[j] == -1.0 ){
          num_mismatches+=1;
        }
      }
    }
    return num_mismatches;
  }
}

// ***************************************************************************
// appendStructurePrototypes - Create new structure prototypes after comparisons
// ***************************************************************************
namespace compare{
  void appendStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
				 vector<StructurePrototype>& final_prototypes){

    // This cleans the StrucuturePrototype objects be removing all the mismatches.
    // Then, it takes the mismatches and makes them into new StructurePrototype objects
    // to be compared.

    ostringstream oss;
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;
    string function_name = "compare::appendStructurePrototypes()";

    vector<StructurePrototype> tmp_list;
    for(uint i=0; i<comparison_schemes.size(); i++){
      bool first_mismatch=true;
      for(uint j=0; j<comparison_schemes[i].misfits.size(); j++){
        if(comparison_schemes[i].misfits[j] > 0.1 || comparison_schemes[i].misfits[j] == -1.0 ){
          // First, store any family prototype information
          if(comparison_schemes[i].misfits[j] > 0.1 && comparison_schemes[i].misfits[j] <= 0.2){
            comparison_schemes[i].family_structures_names.push_back(comparison_schemes[i].proto_structures_names[j]);
            comparison_schemes[i].family_misfits.push_back(comparison_schemes[i].misfits[j]);
          }
          // Take first mismatch and make as the master structure in the new object
          if(first_mismatch==true){
            StructurePrototype tmp;
            tmp.master_structure_name=comparison_schemes[i].proto_structures_names[j];
            tmp.number_types=comparison_schemes[i].proto_structures[j].num_each_type.size();
            tmp.elements=comparison_schemes[i].elements;
            tmp.stoichiometry=comparison_schemes[i].stoichiometry;
            tmp.number_of_atoms=comparison_schemes[i].proto_structures[j].atoms.size();
            tmp.pearson=comparison_schemes[i].pearson;
            tmp.space_group=comparison_schemes[i].space_group;
            tmp.master_structure=comparison_schemes[i].proto_structures[j];
            tmp_list.push_back(tmp);
            comparison_schemes[i].proto_structures_names.erase(comparison_schemes[i].proto_structures_names.begin()+j);
            comparison_schemes[i].proto_structures.erase(comparison_schemes[i].proto_structures.begin()+j);
            comparison_schemes[i].misfits.erase(comparison_schemes[i].misfits.begin()+j);
            j--;
            first_mismatch=false;
          }
          // If not the first mismatch, add as a proto structure in the new object
          else if(first_mismatch==false){
            tmp_list.back().proto_structures_names.push_back(comparison_schemes[i].proto_structures_names[j]);
            tmp_list.back().proto_structures.push_back(comparison_schemes[i].proto_structures[j]);
            tmp_list.back().misfits.push_back(-1.0);
            comparison_schemes[i].proto_structures_names.erase(comparison_schemes[i].proto_structures_names.begin()+j);//
            comparison_schemes[i].proto_structures.erase(comparison_schemes[i].proto_structures.begin()+j);//
            comparison_schemes[i].misfits.erase(comparison_schemes[i].misfits.begin()+j);//
            j--;
          }
        }
      }
      // Store finished (already compared) schemes in final_prototypes
      message << "Identified unique prototype: " << comparison_schemes[i].master_structure_name << endl;
      if(comparison_schemes[i].proto_structures_names.size()==0){
        message << "   No duplicates. " << endl;
      }
      else{
        message << "   " << setw(80) << std::left << "List of duplicates"
                << setw(15) << std::left << "misfit value" << endl;
        message << "   " << setw(80) << std::left 
                << "-----------------------------------------------------------------------------------------------" << endl;
        for(uint d=0;d<comparison_schemes[i].proto_structures_names.size();d++){
          message << "   " << setw(80) << std::left << comparison_schemes[i].proto_structures_names[d]
                  << setw(15) << std::left << comparison_schemes[i].misfits[d] << endl;
        }
      }
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      
      final_prototypes.push_back(comparison_schemes[i]);
      comparison_schemes.erase(comparison_schemes.begin()+i);
      i--;
    }
    // Store newly generated schemes (not compared yet) into comparison_schemes
    comparison_schemes=tmp_list;
  }
}

// ***************************************************************************
// checkPrototypes - Ensure prototypes of different SG are compared
// ***************************************************************************
namespace compare{
  void checkPrototypes(const uint& num_proc, const bool& same_species, 
		       vector<StructurePrototype>& final_prototypes){

    // Checks to see if prototypes of different space groups are similar. 
    // If they are, then we combine the StructurePrototype objects into one.
    // When combining, we keep the "master" prototype as the one with a higher 
    // symmetry (i.e. higher space group). 

    for(uint i=0;i<final_prototypes.size();i++){
      vector<int> store_indices;
      vector<double> store_misfits;
      int min_index=-1;
      double min_misfit=-1.0;
      for(uint j=i;j<final_prototypes.size();j++){
        ostringstream tmp_oss;
        tmp_oss.clear();
        if(
           // If same_species==true
           (same_species==true && 
            matchableSpecies(final_prototypes[i].master_structure,final_prototypes[j].master_structure,same_species)==true &&
            final_prototypes[i].stoichiometry==final_prototypes[j].stoichiometry && 
            final_prototypes[i].pearson==final_prototypes[j].pearson &&
            final_prototypes[i].space_group!=final_prototypes[j].space_group) ||
           // If same_species==false
           (same_species==false && 
            final_prototypes[i].stoichiometry==final_prototypes[j].stoichiometry &&
            final_prototypes[i].pearson==final_prototypes[j].pearson && 
            final_prototypes[i].space_group!=final_prototypes[j].space_group)
          ){
          double final_misfit=-1.0;
          bool scale_volume=true; //default is true
          bool fast_match=false; //default is false
          aflowCompareStructure(num_proc, final_prototypes[i].master_structure, 
					 final_prototypes[j].master_structure, same_species, 
					 scale_volume, fast_match, tmp_oss, final_misfit);
          if(final_misfit < min_misfit){
             min_misfit=final_misfit;
             min_index=j;
          }
        }
      }
      // If one prototype is similar to another, add to one with higher space group
      if(min_misfit!=-1.0){
        int sg_ind=-1;
        uint other_ind=-1;
        if(final_prototypes[i].space_group > final_prototypes[min_index].space_group){
          sg_ind=i;
          other_ind=min_index;
        }
        else{
          sg_ind=min_index;
          other_ind=i;
        }
        // Transfer info to prototype with higher space group
        final_prototypes[sg_ind].proto_structures_names.push_back(final_prototypes[other_ind].master_structure_name);
        final_prototypes[sg_ind].proto_structures_names.insert(final_prototypes[sg_ind].proto_structures_names.end(),
				final_prototypes[other_ind].proto_structures_names.begin(),
				final_prototypes[other_ind].proto_structures_names.end());
        final_prototypes[sg_ind].proto_structures.push_back(final_prototypes[other_ind].master_structure);
        final_prototypes[sg_ind].proto_structures.insert(final_prototypes[sg_ind].proto_structures.end(),
				final_prototypes[other_ind].proto_structures.begin(),
				final_prototypes[other_ind].proto_structures.end());
        // Delete the prototype with the lower space group
        final_prototypes.erase(final_prototypes.begin()+other_ind);
        // If the index deleted was less than the initial loop (i), then need to reduce iterator
        if(other_ind<=i){
          i--;
        }
      }
    }
  }
}

// ***************************************************************************
// prepareTextOutput - Displays results for .txt file
// ***************************************************************************
namespace compare{
  void prepareTextOutput(ostream& ss_out, const bool& same_species, 
			 const vector<StructurePrototype>& final_prototypes){
   
    // Displays comparison information in an TXT file.

    for(uint j=0; j<final_prototypes.size(); j++){
      ss_out << "==============================================================" << endl;
      ss_out << "# ";
      if(same_species==true){
        for(uint k=0;k<final_prototypes[j].elements.size();k++){
          ss_out << final_prototypes[j].elements[k] << final_prototypes[j].stoichiometry[k];
        }
        // TEST ss_out << "   " << final_prototypes[j].pearson
        ss_out << "   " << final_prototypes[j].space_group 
               << "   " << "(" << final_prototypes[j].master_structure_name << ")" << endl;
      }
      else if(same_species==false){
        for(uint k=0;k<final_prototypes[j].stoichiometry.size();k++){
          if(k==0){
            ss_out << final_prototypes[j].stoichiometry[k];
          }
          else{
            ss_out << ":" << final_prototypes[j].stoichiometry[k];
          }
        }
        //ss_out << "   " << final_prototypes[j].pearson
        ss_out << "   " << final_prototypes[j].space_group 
               << "   " << "(" << final_prototypes[j].master_structure_name << ")" << endl;
      }
      ss_out << "==============================================================" << endl;
      for(uint k=0;k<final_prototypes[j].proto_structures_names.size();k++){
	       ss_out << "   " << setw(80) << std::left << final_prototypes[j].proto_structures_names[k] 
	              << setw(15) << std::left << final_prototypes[j].misfits[k] << endl;
      }
    }
  }
}

// ***************************************************************************
// Matchable species
// ***************************************************************************
namespace compare{
  bool matchableSpecies(const xstructure& xstr1, const xstructure& xstr2, 
			const bool& same_species){

    // Determine if it is possible to match up species based on the number of atoms
    // (i.e reduced stoichiometries equal)

    bool LDEBUG=(false || XHOST.DEBUG);
    vector<uint> stoich1;
    vector<uint> stoich2;
    bool matchable=true;
    if(xstr1.species.size()==xstr2.species.size()){
      if(xstr1.species.size()==1){
     	  stoich1.push_back(1); stoich2.push_back(1);
      }
      else{
        stoich1=gcdStoich(xstr1.num_each_type);
        stoich2=gcdStoich(xstr2.num_each_type);
      }
      uint matches=0;
      // Check if we can match to same species (atoms and stoichs)
      if(same_species==true){
        bool commensurate=false;
        for(uint i=0; i<stoich1.size(); i++){
          for(uint j=0; j<stoich2.size(); j++){
            if(stoich1[i]==stoich2[j] && xstr1.species[i]==xstr2.species[j]){
              //cerr << "matching: " << stoich1[i] << "==" << stoich2[j] << " && " << xstr1.species[i] << "==" << xstr2.species[j] << endl;
              matches++;
              commensurate=true;
              break;
            }
	  }
          if(commensurate==false){
            matchable=false;
            break;
          }
	      }
      }
      // Check if we can match stoichs only
      else{
        for(uint i=0; i<stoich1.size(); i++){
          std::sort(stoich1.begin(),stoich1.end());
          std::sort(stoich2.begin(),stoich2.end());
        }
        for(uint i=0; i<stoich1.size(); i++){
          if(stoich1[i]!=stoich2[i]){
            matchable=false;
            break;
          }
          else{
            matches++;
          }
        }
      }
      if(matchable==true && matches==stoich1.size()){
        //cerr << "match found" << endl;
	      return true;
      }
      else{
	      return false;
      }
    }
    else{
      if(LDEBUG){cerr << "compare:: " << "NUMBER OF TYPES OF ATOMIC SPECIES:   xstr1:  " << xstr1.num_each_type.size()
	                    << "           xstr2:  " << xstr2.num_each_type.size() << endl;}
      if(LDEBUG){cerr << "compare:: " << "NUMBER OF TYPES OF ATOMIC SPECIES IS NOT THE SAME...QUITTING..." << endl;}
      return false;
    }
  }
}

// ***************************************************************************
// Same Species 	
// ***************************************************************************
namespace compare{
  bool sameSpecies(const xstructure& xstr1, const xstructure& xstr2, const bool& display){

    bool LDEBUG=(false || XHOST.DEBUG);

    // Check number of types
    if(xstr1.num_each_type.size() != xstr2.num_each_type.size()){
      // Display counts
      if(display==true){
        if(LDEBUG){
          cerr << "compare::NUMBER OF TYPES OF ATOMIC SPECIES:   xstr1:  " << xstr1.num_each_type.size() 
               << " 		xstr2:  " << xstr2.num_each_type.size() << endl;
          cerr << "compare::NUMBER OF TYPES OF ATOMIC SPECIES IS NOT THE SAME..." << endl;
        }
        return false;
      }
    }

    // Check counts and species
    for(uint i=0;i<xstr1.num_each_type.size();i++){
      bool matched = false;
      for(uint j=0;j<xstr2.num_each_type.size();j++){
        // DX IS SPECIES CHECK TOO STRICT? if(xstr1.num_each_type[i] == xstr2.num_each_type[j] &&
        // DX IS SPECIES CHECK TOO STRICT?   xstr1.species[i] == xstr2.species[j]){
        if(xstr1.num_each_type[i] == xstr2.num_each_type[j]){
          matched = true;
          break;
        }
      }
      if(matched == false){
        if(display==true){ 
          if(LDEBUG){
            cerr << "compare::WARNING:: TYPE OF ATOMIC SPECIES OR NUMBER PER TYPE ARE NOT THE SAME..." << endl;  
          }
        }
        return false;
      }
    }
    if(display==true){
      if(LDEBUG){
	      cerr << "compare::NUMBER AND TYPE OF ATOMIC SPECIES ARE THE SAME...PROCEEDING..." << endl;
      }
    }
    return true;
  }
}

// ***************************************************************************
// Rescale Structure
// ***************************************************************************
namespace compare{
  void rescaleStructure(xstructure& xstr1, xstructure& xstr2){

    // If the scale factor is different, the two structures are rescaled to 1.00 

    if(abs(xstr1.scale-xstr2.scale)>0.001){
       xstr1.ReScale(1.0);
       xstr2.ReScale(1.0);
       xstr1.FixLattices();
       xstr2.FixLattices();
    }
  }
}

// ***************************************************************************
// Atomic Number Desnity
// ***************************************************************************
namespace compare{
  void atomicNumberDensity(xstructure& xstr1, xstructure& xstr2) { 

    // In order to compare structure with different volumes
    // we rescale the cell so that the volume divided by
    // the number of atoms is the same.

    double scale;
    scale=(xstr1.Volume()/xstr1.atoms.size())/(xstr2.Volume()/xstr2.atoms.size());
    xstr2.InflateVolume(scale);
  }
}

// ***************************************************************************
// Fake Atoms Name
// ***************************************************************************
namespace compare{
  void fakeAtomsName(xstructure& xstr){

    // Assign a fake letter to each atom type. In case of materials with more 
    // than 26 species it is necessary to add more characters to this string

    string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int iat=0;

    for(uint i=0; i<xstr.num_each_type.size(); i++){
      xstr.species[i]=letters[i];
      for(int j=0; j<xstr.num_each_type[i]; j++){
        xstr.atoms[iat].name=letters[i];
	iat++;
      }
    }
  }
}

// ***************************************************************************
// Print Structure Parameters
// ***************************************************************************
namespace compare{
  void printParameters(xstructure& xstr, ostream& oss) { 
    oss << "========================================================" << endl;
    oss << xstr.title << endl;
    oss << "Parameters:" << endl;
    oss << "a:     " << xstr.a << endl;
    oss << "b:     " << xstr.b << endl;
    oss << "c:     " << xstr.c << endl;
    oss << "alpha: " << xstr.alpha << endl;
    oss << "beta:  " << xstr.beta << endl;
    oss << "gamma: " << xstr.gamma << endl;
    oss << "volume:" << xstr.GetVolume() << endl;
  }
}

// ***************************************************************************
// Least Frequent Atom
// ***************************************************************************
namespace compare{
  string leastFrequentAtom(const xstructure& xstr) {

    // Least frequent atom can be exploited in case we want to reduce the number 
    // of quadruplets to compute the rotational matrix by considering only the 
    // LFA species. It differs from other leastfrequentatom2 function because it 
    // finds the first LFA atom in xstructure 1. We do not need to find more the 
    // one LFA species in xstructure 1 because we will find all the LFAs for 
    // xstructure 2 and if the structures are a match we are guaranteed to find 
    // a one-to-one correspondence from one of the LFAs in xstructure 2.

    int flag=0, leastFrequentAtomCount=0;
    string leastFrequentAtomType;


    for(uint i=0; i<xstr.num_each_type.size(); i++){
      if(flag==0){
        // Number of occurrences of the LFA
        leastFrequentAtomCount = xstr.num_each_type[i];   
        // LFA Species 
        leastFrequentAtomType = xstr.species[i];         
        flag=1;
      }
      else{
        if(leastFrequentAtomCount>xstr.num_each_type[i]){
          leastFrequentAtomCount = xstr.num_each_type[i];
          leastFrequentAtomType = xstr.species[i];
        }   
      }   
    }   
    return leastFrequentAtomType;
  }
}

// ***************************************************************************
// Least Frequent Atom 2
// ***************************************************************************
namespace compare{
  vector<string> leastFrequentAtom2(const xstructure& xstr) {

    // This least frequent atom function finds all possible least frequent atoms 
    // for an xstructure and stores them in a vector. All of these LFAs are used 
    // in the quadruplet search.
    // We may not need to search over multiple LFAs during the quadruplet search. 
    // If a match is not found for one LFA, it won't be found for another since 
    // we need to map all atoms in one structure to the other structure. We will 
    // leave this implementation in for now, but may speed up the quadruplet 
    // search if we consider only one LFA.

    int flag=0, leastFrequentAtomCount=0;
    vector<string> leastFrequentAtomType;

    for(uint i=0; i<xstr.num_each_type.size(); i++){
      if(flag==0){
        // Number of occurrences of the LFA
        leastFrequentAtomCount = xstr.num_each_type[i];   
        // LFA Species 
        leastFrequentAtomType.push_back(xstr.species[i]);  
        flag=1;
      }
      else{
        if(leastFrequentAtomCount>xstr.num_each_type[i]){
          leastFrequentAtomCount = xstr.num_each_type[i];
          leastFrequentAtomType.clear();
          leastFrequentAtomType.push_back(xstr.species[i]);
        }
        if(leastFrequentAtomCount==xstr.num_each_type[i] && leastFrequentAtomType[0]!=xstr.species[i]){ 
          // Added the statement after '&&' (above); ensures no double counting from previous if statement
          leastFrequentAtomType.push_back(xstr.species[i]);
        }
      }
    }
    return leastFrequentAtomType;
  }
}

// ***************************************************************************
// Check Tolerances
// ***************************************************************************
namespace compare{
  bool checkTolerance(xvector<double> d1, xvector<double> d2){

    // When we look for 2 corresponding reference frames, we have to check that
    // the length of the 3 vectors and the angles between them are equal or within
    // a given tolerance (in this case 10% of the reference value has been set).

    double tol_length=0.3, tol_angle=0.3;

    if( abs(d1(1)-d2(1)) < tol_length*abs(d1(1)) &&
        abs(d1(2)-d2(2)) < tol_length*abs(d1(2)) &&
        abs(d1(3)-d2(3)) < tol_length*abs(d1(3)) &&
        abs(d1(4)-d2(4)) < tol_angle*abs(d1(4)) &&
        abs(d1(5)-d2(5)) < tol_angle*abs(d1(5)) &&
        abs(d1(6)-d2(6)) < tol_angle*abs(d1(6))
      ){
      return false;
    }
    else{	
      return true;
    }
  }
}

// ***************************************************************************
// Check ABC Tolerances
// ***************************************************************************
namespace compare{
  bool checkABCTolerance(xvector<double> d1, xvector<double> d2){

    // Similar to checkTolerance, but it only looks at the length of the lattice 
    // vectors (for screening structures; makes faster)

    double tol_length=0.3;

    for(uint i=1;i<4;i++){
      for(uint j=1;j<4;j++){
	if(j!=i){
	  for(uint k=1;k<4;k++){
	    if(k!=i && k!=j){
	      if(abs(d1(i)-d2(1)) < tol_length*abs(d1(i)) &&
		 abs(d1(j)-d2(2)) < tol_length*abs(d1(j)) &&
		 abs(d1(k)-d2(3)) < tol_length*abs(d1(k))
		){
		return false;
	      }
	    }
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// Check Angle Tolerances
// ***************************************************************************
namespace compare{
  bool checkAngleTolerance(xvector<double> d1, xvector<double> d2){

    // Similar to checkTolerance, but it only looks at the angles betweeen the 
    // lattice vectors (for screening structures; makes faster)

    double tol_angle=0.3;

    for(uint i=4;i<7;i++){
      for(uint j=4;j<7;j++){
	if(j!=i){
	  for(uint k=4;k<7;k++){
	    if(k!=i && k!=j){	
	      if(abs(d1(i)-d2(4)) < tol_angle*abs(d1(4)) &&
		 abs(d1(j)-d2(5)) < tol_angle*abs(d1(5)) &&
		 abs(d1(k)-d2(6)) < tol_angle*abs(d1(6))
		){
		return false;
	      }
	    }
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// Find centroid for system with periodic boundary conditions
// ***************************************************************************
namespace compare{
  xvector<double> centroid_with_PBC(const xstructure& xstr){
    // Calculate the "best" centroid in a system with periodic boundary conditions.
    // This is based on the algorithm proposed in: 
    // https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    // Used to find the best origin/centroid for a crystal structure.
    vector<xvector<double> > coordinates;
    for(uint i=0;i<xstr.atoms.size();i++){
      coordinates.push_back(xstr.atoms[i].cpos); //or cpos
    }
    return centroid_with_PBC(coordinates,xstr.lattice);
  }
}

// ***************************************************************************
// Find centroid for system with periodic boundary conditions
// ***************************************************************************
namespace compare{
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, const xmatrix<double>& lattice){
    // Calculate the "best" centroid in a system with periodic boundary conditions.
    // This is based on the algorithm proposed in: 
    // https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    // If there are no weights (geometric center), then the weights are set to 1
    vector<double> weights;
    for(uint i=0;i<coordinates.size();i++){
      weights.push_back(1.0);
    }
    return centroid_with_PBC(coordinates,weights,lattice);
  }
}

// ***************************************************************************
// Find centroid for system with periodic boundary conditions
// ***************************************************************************
namespace compare{
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, vector<double>& weights, 
				    const xmatrix<double>& lattice){
    // Calculate the "best" centroid in a system with periodic boundary conditions.
    // This is based on the algorithm proposed in: 
    // https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions

    xvector<double> centroid;
    for(uint i=1;i<4;i++){
      double zi_avg = 0.0;
      double zeta_avg = 0.0;
      double theta_avg =0.0;
      for(uint j=0;j<coordinates.size();j++){
        double theta = coordinates[j][i]*(2.0*Pi_r)/(aurostd::modulus(lattice(i)));
	double zi = std::cos(theta)*weights[j];
	double zeta = std::sin(theta)*weights[j];
	zi_avg += zi/coordinates.size();
	zeta_avg += zeta/coordinates.size();
      }
      theta_avg = std::atan2(-zeta_avg,-zi_avg);
      centroid(i) = theta_avg*(aurostd::modulus(lattice(i))/(2.0*Pi_r));
    }
    return centroid;
  } 
}


// ***************************************************************************
// Find Matches
// ***************************************************************************
namespace compare{
  bool findMatch(const xstructure& xstr1, const xstructure& PROTO, 
		    vector<uint>& im1, vector<uint>& im2, vector<double>& min_dists, 
                    const int& type_match) {

    // In order to find the best matchings the routine computes 
    // the difference between one atom and all the others, 
    // building a matrix of differences of coordinates.
    // Then, it checks which atoms have the best matching
    // with another atom in the second structure.

    // A | 1    A1, A2     I can check which is the best matching
    // B | 2 -> B1, B2  -> for the atom 1 and 2 in the structure
    // C |      C1, C2     with A,B,C,D 
    // D |      D1, D2    

    uint j=0,k=0;
    int i1=0,i2=0;                                  //Indices corresponding atoms

    vector<double> vdiffs;                      //Difference btwn atoms coords
    vector<std::pair<xvector<double>,xvector<double> > > min_positions;                      //Store sets of Cartesian coords which minimize distance
    vector<vector<double> > all_vdiffs;         //For all the atoms

    vector<uint> im1_tmp;
    vector<uint> im2_tmp;
    vector<string> im1_name;
    vector<string> im2_name;
    im1.clear();
    im2.clear();
    vdiffs.clear();
    all_vdiffs.clear();

    xmatrix<double> klattice=PROTO.lattice;

    double tmp = 1e9;
    uint i1_real=0;
    uint i2_real=0;
    string i1_name = "";
    string i2_name = "";

    xvector<double> best_centroid1 = centroid_with_PBC(xstr1); 
    xvector<double> best_centroid2 = centroid_with_PBC(PROTO); 

    for(j=0;j<xstr1.atoms.size();j++){
      std::pair<xvector<double>,xvector<double> > tmp_pair;
      xvector<double> tmp_xvec = xstr1.atoms[j].cpos;
      tmp_pair.first = tmp_xvec;
      vdiffs.clear();
      double match_dist=1e9;
      for(k=0;k<PROTO.atoms.size();k++){
        double dist=1e9;
        xvector<double> min_xvec;
        xvector<int> abc;
        // Need to find the min distance; thus check distance between neighboring cells to find true minimum.
        for(int a=-2;a<=2;a++){
          for(int b=-2;b<=2;b++){
            for(int c=-2;c<=2;c++){
              tmp_xvec = PROTO.atoms[k].cpos+a*klattice(1)+b*klattice(2)+c*klattice(3);
              tmp = aurostd::modulus(xstr1.atoms[j].cpos-tmp_xvec);
              if(tmp < dist){
                abc(1)=a;
                abc(2)=b;
                abc(3)=c;
                i1 = j;
                i2 = k;
                //cerr << tmp << " vs " << dist << " - " << j << ", " << k << endl;
                dist = tmp;
                min_xvec = tmp_xvec;
              }
              // DXdist=aurostd::min(dist,aurostd::modulus(PROTO.atoms[j].cpos-xstr1.atoms[k].cpos+a*klattice(1)+b*klattice(2)+c*klattice(3)));
            }
          }
        }
        //cerr << "match_dist: " << match_dist << endl;
        if(dist<match_dist){
          i1_real=i1;
          i2_real=i2;
          match_dist = dist;
          i1_name = xstr1.atoms[i1_real].name;
          i2_name = PROTO.atoms[i2_real].name;
          tmp_pair.second = min_xvec;
        }
        vdiffs.push_back(dist);
      }
      // Check if same species match
      if(type_match == 2){ // same species
        if(i1_name != i2_name){
          //cerr << "WARNING: Matching species are not the same type, throwing out match (same species comparison)" << endl;
          return false;
        }
      }
      // Check for one-to-one mappings 
      for(uint i=0;i<im1_name.size();i++){
        // Check if i1 index has multiple mappings
        if(i1_real == im1[i]){
          //cerr << "WARNING: Used the same index for matching in i1! (" << i1_real << " == " << im1[i] << ")"<< endl;
          //cerr << match_dist << " vs " << min_dists[i] << endl;
          return false;
        }
        // Check if i2 index has multiple mappings
        if(i2_real == im2[i]){
          //cerr << "WARNING: Used the same index for matching in i2! (" << i2_real << " == " << im2[i] << ")"<< endl;
          //cerr << match_dist << " vs " << min_dists[i] << endl;
          return false;
        }
        // Check if types are not consistently mapped to a single type
        if(i1_name == im1_name[i]){
          if(i2_name != im2_name[i]){
            //cerr << "WARNING: Matching one type of atom to more than one type! (" << i1_name << " == " << i2_name << " | " << im1_name[i] << " == " << im2_name[i] << ")" <<  endl;
            //for(uint j=0;j<im1_name.size();j++){
            //  cerr << im1[j] << " == " << im2[j] << " | " << xstr1.atoms[im1[j]].cpos << " == " << PROTO.atoms[im2[j]].cpos << " (" << min_dists[i] << ") | " << im1_name[j] << " == " << im2_name[j] << endl;
            //}
            //cerr << i1_real << " == " << i2_real << " | " << xstr1.atoms[i1_real].cpos << " == " << PROTO.atoms[i2_real].cpos << " (" << match_dist << ") | " << i1_name << " == " << i2_name << endl;            
            return false;
          }
        }
      }
      im1_tmp.push_back(i1_real);
      im2_tmp.push_back(i2_real);
      im1.push_back(i1_real);
      im2.push_back(i2_real);
      im1_name.push_back(i1_name);
      im2_name.push_back(i2_name);
      min_dists.push_back(match_dist);
      all_vdiffs.push_back(vdiffs);
      min_positions.push_back(tmp_pair);
    }
    return true;
  }
}

// ***************************************************************************
// clusterize
// ***************************************************************************
namespace compare{
  void clusterize(const xstructure& xstr1, const vector<uint>& im1, vector<string>& TYPE1, 
		  vector<uint>& QTA1, vector<vector<uint> >& I1){

    // This function builds clusters/vectors of atoms of the same type. It is necessary when
    // we want to check correspondences between atoms of the same type or between atoms
    //of a specific species.

    uint i,j;
    int flag;
    vector<uint> i1;

    i1.clear();

    for(i=0; i<im1.size(); i++){
      i1.clear();
      if(i==0){
	TYPE1.push_back(xstr1.atoms[im1[i]].name);
	QTA1.push_back(1);
	i1.push_back(im1[i]);
	I1.push_back(i1);
      }
      else{
	flag=0;
	for(j=0; j<TYPE1.size(); j++){
	  if(xstr1.atoms[im1[i]].name==TYPE1[j]){
	    QTA1[j]++;
	    i1=I1[j];
	    i1.push_back(im1[i]);
	    I1[j]=i1;
	    flag=1;
	    break;
	  }
	}
	if(flag==0){
	  TYPE1.push_back(xstr1.atoms[im1[i]].name);
	  QTA1.push_back(1);
	  i1.push_back(im1[i]);
	  I1.push_back(i1);
	}
      }
    }
  }
}

// ***************************************************************************
// Same Atom Type
// ***************************************************************************
namespace compare{
  bool sameAtomType(const xstructure& xstr1, const xstructure& xstr2, const vector<uint>& im1, 
		    const vector<uint>& im2, const int& type_match){

    // Bool function; checks when 2 atoms matched are OK 

    uint i,j,k,w,z;
    vector<string> TYPE1, TYPE2;
    vector<uint> QTA1,QTA2;
    vector<vector<uint> > I1,I2;

    vector<int> flag,checkType;

    clusterize(xstr1,im1,TYPE1,QTA1,I1);
    clusterize(xstr2,im2,TYPE2,QTA2,I2);

    checkType.clear();

    for(i=0; i<TYPE1.size(); i++){
      for(j=0; j<TYPE2.size(); j++){
	if(QTA1[i]==QTA2[j]){	
	  flag.clear();
	  for(w=0; w<I1[i].size(); w++){
	    for(z=0; z<I2[j].size(); z++){
	      for(k=0; k<im1.size(); k++){
		//Means we want the same atomic species to be matched up.
		if(type_match==2){           
		  if(I1[i][w]==im1[k] && I2[j][z]==im2[k] && 
		     xstr1.atoms[im1[k]].name == xstr2.atoms[im2[k]].name
		    ) 
		    flag.push_back(1);
		}
		else{
		  if(I1[i][w]==im1[k] && I2[j][z]==im2[k]) flag.push_back(1);
		}
	      }
	    }
	  }
	  if(flag.size()==QTA1[i])
	    checkType.push_back(1);	
	}
      }
    }
    if(checkType.size()==TYPE1.size()) return true;
    else return false;
  }
}

// ***************************************************************************
// Clean Match
// ***************************************************************************
namespace compare{
  bool cleanMatch(const vector<uint>& im1) {
    
    // The result of the findMatch function is a pair of vectors containing the indices of the matched
    // atoms of the structures. This function allows to delate the matchings where the same
    // index appears more the one time.

    uint i,j;

    for(i=0; i<im1.size(); i++){
      for(j=0; j<im1.size(); j++){
	if(i!=j){
	  if(im1[i]==im1[j]){
	    return true;
          }
        }
      }	
    }
    return false;
  }
}

// ***************************************************************************
// Cell Diagonal
// ***************************************************************************
namespace compare{
  void cellDiagonal(xstructure& xstr, vector<double>& diag_sum, 
		    vector<double>& diag_diff, const double& scale) {
    xmatrix<double> lattice = xstr.lattice;
    cellDiagonal(lattice, diag_sum, diag_diff, scale);
  }
}
namespace compare{
  void cellDiagonal(xmatrix<double>& lattice, vector<double>& diag_sum, 
		    vector<double>& diag_diff, const double& scale) {

    // The cell diagonals are represented by the shortest and longest diagonals
    // (in the case of cubic lattice the 2 diagonals are equal)
    // They are obtained as the vectorial sum and difference of the lattice
    // basis vectors.

    xvector<double> origin;

    diag_sum.push_back(abs(distance(origin,(lattice(1)+lattice(2))*scale)));
    diag_sum.push_back(abs(distance(origin,(lattice(1)+lattice(3))*scale)));
    diag_sum.push_back(abs(distance(origin,(lattice(2)+lattice(3))*scale)));
    diag_diff.push_back(abs(distance(lattice(1),lattice(2))*scale));
    diag_diff.push_back(abs(distance(lattice(1),lattice(3))*scale));
    diag_diff.push_back(abs(distance(lattice(2),lattice(3))*scale));
  }
}

// ***************************************************************************
// Lattice Deviation
// ***************************************************************************
namespace compare{
  double latticeDeviation(const vector<double>& diag_sum1,const vector<double>& diag_sum2, 
			  const vector<double>& diag_diff1,const vector<double>& diag_diff2) {

    // Lattice deviation is computed as the deviation of the 2 diagonals of each face
    // normalized on the sum of the diagonals of the faces of the reference structure.
    // The images of the diagonals of the mapped structure must be rescaled such that
    // the volume of the image of its unit cell of equals the volume of the unit
    // cell of the reference one.

    uint i;
    double d;
    vector<double> dev;

    for(i=0; i<3; i++) 
      dev.push_back((abs(diag_sum1[i]-diag_sum2[i])+abs(diag_diff1[i]-diag_diff2[i]))/(diag_sum2[i]+diag_diff2[i]));

    d=1;
    for(i=0;i<dev.size();i++) 
      d=d*(1-dev[i]);
    d=1-d;

    return d;
  }
}

// ***************************************************************************
// Compute nearest neighbors 
// ***************************************************************************
namespace compare{
  vector<double> computeNearestNeighbors(xstructure& xstr){
    vector<double> all_nn_distances;
    double nn = 1e9;
    for(uint i=0;i<xstr.atoms.size();i++){
      nn = shortestDistance(xstr,i);
      all_nn_distances.push_back(nn);
    }
    return all_nn_distances;
  }
}

// ***************************************************************************
// Shortest Distance from one atom
// ***************************************************************************
namespace compare{
  double shortestDistance(const xstructure& xstr, const uint& k) {

    // Find the distance with the closest distance

    //double min;
    double min_dist=1e9;
    xmatrix<double> klattice = xstr.lattice; //NEW
    for(uint ii=0; ii<xstr.atoms.size(); ii++){
      if(ii!=k){
	for(int a=-2;a<=2;a++){
	  for(int b=-2;b<=2;b++){
	    for(int c=-2;c<=2;c++){
              double tmp = aurostd::modulus(xstr.atoms[k].cpos-xstr.atoms[ii].cpos+a*klattice(1)+b*klattice(2)+c*klattice(3));
	      min_dist=aurostd::min(min_dist,tmp);
	    }
	  }
	}      
      }
    }
    return min_dist;
  }
}

// ***************************************************************************
// Coordinates Deviation
// ***************************************************************************
namespace compare{
  void coordinateDeviation(const xstructure& xstr1, const xstructure& xstr2, 
                      const vector<double>& all_nn1, const vector<double>& all_nn_proto,
		      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, vector<double>& min_dists,
		      double& cd, double& fail_figure) {

    // To compute the coordinates deviation we look at each pair of atoms from
    // the reference and mapped structure

    uint j;
    double num=0, den=0, nfail=0;
    double dd, nn1, nn2; //dd=delta distance, nn=nearest neighbour 
    int  fail1, fail2;
    xmatrix<double> klattice = xstr1.lattice;
    for(j=0; j<indexMatch1.size(); j++){
      //nn1=shortestDistance(xstr1,indexMatch1[j]);
      //nn2=shortestDistance(xstr2, indexMatch2[j]);
      nn1 = all_nn1[indexMatch1[j]];
      nn2 = all_nn_proto[indexMatch2[j]];
      dd = min_dists[j];
//      dd=1e9;
//      // Need to find the min distance; thus check distance between neighboring cells to find true minimum.
//      for(int a=-2;a<=2;a++){
//	for(int b=-2;b<=2;b++){
//	  for(int c=-2;c<=2;c++){
//	    dd=aurostd::min(dd,modulus(xstr1.atoms.at(indexMatch1[j]).cpos-xstr2.atoms.at(indexMatch2[j]).cpos+a*klattice(1)+b*klattice(2)+c*klattice(3)));
//	  }
//	}
//      }

      if(dd<=0.5*nn1) fail1=0;
      if(dd>0.5*nn1) fail1=1;
      if(dd<=0.5*nn2) fail2=0;
      if(dd>0.5*nn2) fail2=1;

      if(fail1==0){
	num=num+dd;
	den=den+nn1;
      }   
      if(fail2==0){
	num=num+dd;
	den=den+nn2;
      }   
      if(fail1==1) nfail++;
      if(fail2==1) nfail++;
    }   

    if(den==0) cd=1;
    else cd=num/den;
    
    //Consider unmatched atoms
    int flag=0;
    for(uint i=0; i<xstr1.atoms.size();i++){
      flag=0;
      for(uint k=0; k<indexMatch1.size();k++){
	if(i==indexMatch1[k]){
	  flag=1;
	}
      }
      if(flag==0){	//Meaning this atom does not have a match; increase the figure of failure
	nfail++;
      }
    }
    for(uint i=0; i<xstr2.atoms.size();i++){
      flag=0;
      for(uint k=0; k<indexMatch2.size();k++){
	if(i==indexMatch2[k]){
	  flag=1;
	}
      }
      if(flag==0){        //Meaning this atom does not have a match; increase the figure of failure
	nfail++;
      }
    }
    
    fail_figure=(nfail/(xstr1.atoms.size()+xstr2.atoms.size()));
  }   
}

// ***************************************************************************
// Compute Misfit
// ***************************************************************************
namespace compare{
  double computeMisfit(const double& dev, const double& dis, const double& fail) {

    // Combines differences between all aspects of crystal structure into a 
    // figure of misfit. (See Burzlaff)

    double mis;

    mis=1-((1-dev)*(1-dis)*(1-fail));
    return mis;
  }   
}

// ***************************************************************************
// Print Matching Between Atoms
// ***************************************************************************
namespace compare{
  void printMatch(const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, 
		  const xstructure& PROTO, const xstructure& xstr1, ostream& oss) {

    // With this function we print the atoms matched in the previous function
    // whose indices are stored in the indexMatch vector.
    // This allows us to call them directly.
    // Lastly, the translational term, specific of each mapping, is printed.

    uint i,j;

    oss << "              Reference                               Mapped"<<endl;
    for(j=0; j< indexMatch1.size(); j++){
      oss << indexMatch1[j]<<"-"<<indexMatch2[j]<<"    ";
      oss << xstr1.atoms[indexMatch1[j]].cpos << " " << xstr1.atoms[indexMatch1[j]].name;
      oss << "       ";
      oss << PROTO.atoms[indexMatch2[j]].cpos << " " << PROTO.atoms[indexMatch2[j]].name << endl;
    }

    int flag=0;

    oss <<"----------------------------------------------------"<<endl;
    oss << "Missing Atoms in Reference structure:"<< endl;
    for(j=0; j<xstr1.atoms.size(); j++){
      flag=0;
      for(i=0; i<indexMatch1.size(); i++){
	if(j==indexMatch1[i])
	  flag=1;
      }
      if(flag==0){
	oss << "# "<< j << "   " << xstr1.atoms[j].cpos << "   " << xstr1.atoms[j].name << endl;
      }
    }

    oss << "Missing Atoms in Mapped structure:" << endl;
    for(j=0; j<PROTO.atoms.size(); j++){
      flag=0;
      for(i=0; i<indexMatch2.size(); i++){
	if(indexMatch2[i]==j) flag=1;
      }
      if(flag==0) oss << "# "<< j << "   " << PROTO.atoms[j].cpos << "   " << PROTO.atoms[j].name << endl;
    }
  }
}

// ***************************************************************************
// Bring Coordinate in the cell (Similar to xstructure.BringInCell())
// ***************************************************************************
namespace compare{
  xvector<double> bringCoordinateInCell(xvector<double>& coord){

    // This function brings a coordinate back in the unit cell.  
    // It is not an xstructure attribute, so we cannot use xstructure.BringInCell().

    double tol=1e-6;
    for(uint i=1;i<4;i++){
      if(coord(i)<-tol){
	coord(i)+=1.0;
      }
      else if(coord(i)>=1.0-tol){
	coord(i)-=1.0;
      }
    }
    return coord;
  }
}

// ***************************************************************************
// Determine if vector is Periodic
// ***************************************************************************
namespace compare{
  bool vectorPeriodic(const xvector<double>& vec, const xstructure& lfa_supercell, 
		       const int& i, const int& j){

    // Once we have a possible quadruplet (lattice), we need to make sure that this 
    // choice of the primitive cell preserves the periodicity o the lattice. 
    // Therefore, we check that each of the quadruplet atoms maps onto another atom 
    // in the supercell. Helpful analogy: Lattice periodicty vs crystal periodicity. 
    // The quadruplets form the lattice and in this function we check for lattice
    // periodicity. The misfit criteria checks the crystal periodicity.

    double tolerance = 0.01; // Hundredth of an angstrom
    deque<_atom> atoms = lfa_supercell.atoms;
    xmatrix<double> lattice = lfa_supercell.lattice;
    xmatrix<double> f2c = trasp(lattice);
    xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = false;

    vector<int> ind(2); ind[0]=i, ind[1]=j;
    xvector<double> tmp;

    uint count=0;
    deque<_atom> transformed;
    deque<uint> index_to_check;

    // ===== Check if applying the symmetry element along with internal translation maps to another atom ===== //
    for(uint d=0;d<atoms.size();d++){
      _atom tmp;
      tmp.type = atoms[d].type;
      tmp.cpos = atoms[d].cpos+vec;
      tmp.fpos = C2F(lattice,tmp.cpos);
      if(SYM::MapAtom(atoms,tmp,true,c2f,f2c,skew,tolerance)){
        transformed.push_back(tmp);
        index_to_check.push_back(d);
        count++;
      }
      else{
        return false;
      }
    }
    if(count == atoms.size()){
      return true;
    }
    return false;
  }
}

// ***************************************************************************
// Determine if LFA Quadruplet is Periodic
// ***************************************************************************
namespace compare{
  bool quadrupletPeriodic(const xmatrix<double>& quad, const xstructure& lfa_supercell, 
		       const int& i, const int& j, const int& k, const int& w){

    // Once we have a possible quadruplet (lattice), we need to make sure that this 
    // choice of the primitive cell preserves the periodicity o the lattice. 
    // Therefore, we check that each of the quadruplet atoms maps onto another atom 
    // in the supercell. Helpful analogy: Lattice periodicty vs crystal periodicity. 
    // The quadruplets form the lattice and in this function we check for lattice
    // periodicity. The misfit criteria checks the crystal periodicity.

    vector<int> ind(4); ind[0]=i, ind[1]=j; ind[2]=k; ind[3]=w;
    vector<xvector<double> > latt_vecs(3); 
    latt_vecs[0]=quad(1); latt_vecs[1]=quad(2); latt_vecs[2]=quad(3);
    xvector<double> tmp;

    for(uint b=0; b<ind.size();b++){
      for(uint c=0;c<latt_vecs.size();c++){
	tmp=lfa_supercell.atoms.at(ind[b]).cpos+latt_vecs[c];
	bool match_found=false;
	for(uint d=0;d<lfa_supercell.atoms.size();d++){
	  if(abs(lfa_supercell.atoms.at(d).cpos(1)-tmp(1))<0.01 && 
	     abs(lfa_supercell.atoms.at(d).cpos(2)-tmp(2))<0.01 && 
	     abs(lfa_supercell.atoms.at(d).cpos(3)-tmp(3))<0.01){ // Less than hundredth of Angstrom?
	    match_found=true;
	    break;
	  }
	}
	if(match_found==false){
	  xvector<double> tmp_frac=C2F(lfa_supercell.lattice,tmp);
	  tmp_frac=bringCoordinateInCell(tmp_frac);
	  xstructure lfa_supercell_tmp = lfa_supercell;
	  for(uint f=0;f<lfa_supercell_tmp.atoms.size();f++){
	    lfa_supercell_tmp.atoms.at(f).fpos=C2F(lfa_supercell_tmp.lattice,lfa_supercell_tmp.atoms.at(f).cpos);
	    if(abs(lfa_supercell_tmp.atoms.at(f).fpos(1)-tmp_frac(1))<0.01 && 
	       abs(lfa_supercell_tmp.atoms.at(f).fpos(2)-tmp_frac(2))<0.01 && 
	       abs(lfa_supercell_tmp.atoms.at(f).fpos(3)-tmp_frac(3))<0.01){
	      match_found=true;
	      break;
	    }
	  }
	  if(match_found==false){
	    return false;
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// Thread Generation (For parallel processing of quadruplets)
// ***************************************************************************
namespace compare{
  void threadGeneration(const uint& num_proc,xmatrix<double>& q1, xstructure& xstr2, 
			vector<xstructure> &vprotos, xstructure &xstr1, const int& type_match, 
			const bool& fast_match, double& minMis, ostream& oss){ 

    // This function creates the supercell of the second structure and begins 
    // the quadruplets search within a supercell. Due to the costly nature of 
    // this algorithm, the quadruplet search is parallelized. The splitting 
    // of computation of done here

    bool LDEBUG=(false || XHOST.DEBUG);

    xstructure xstr1_tmp = xstr1;

    //cerr << "LFA" << endl;
    vector<string> LFA_str1=leastFrequentAtom2(xstr1);
    vector<string> LFA_str2=leastFrequentAtom2(xstr2);
    string lfa, lfa_str1;

    //cerr << "SUPERCELL" << endl;
    xstructure xstr=GetSuperCell(xstr2,3,0,0,0,3,0,0,0,3);

    uint y=0;
    uint x=0;

    // DX TEST
    // Consider all LFAs
    for(y=0;y<LFA_str2.size();y++){
    for(x=0;x<LFA_str1.size();x++){
    // DX TEST
      lfa_str1=LFA_str1[x];
      lfa=LFA_str2[y];
      if(type_match == 2 && lfa_str1 != lfa){
        continue;
      }
      oss << "===> LFA: "<<lfa<<endl;

      if(LDEBUG){
        cerr << "===> LFA_1: " << lfa_str1 <<endl;
        cerr << "===> LFA: "<<lfa<<endl;
      }
      
      // DX NEW so we don't need to do this in an inner loop
      for(uint a=0;a<xstr.atoms.size();a++){
        if(xstr.atoms[a].name == lfa){
          xstr.ShifOriginToAtom(a);
          break;
        }
      }
      // DX NEW

      //cerr << "xstr1 centroid: " << endl;

      //cerr << "SHIFT" << endl;
      // NEED TO SHIFT origin of xstr1_tmp to one of the LFA (this was missing before and caused ICSD_102428.BCA, and CBA to not match, but they should
      for(uint i=0;i<xstr1_tmp.atoms.size();i++){
        if(xstr1_tmp.atoms[i].name==lfa_str1){
          xstr1_tmp.ShifOriginToAtom(i);
          xstr1_tmp.BringInCell(1e-10);
          break;
        }
      }

      // NEED TO SHIFT origin of xstr2 to one of the LFA
      for(uint i=0;i<xstr2.atoms.size();i++){
        if(xstr2.atoms[i].name==lfa){
          xstr2.ShifOriginToAtom(i);
          xstr2.BringInCell(1e-10);
          break;
        }
      }

      // When checking the quadruplets/lattice, we only need to generate a 
      // LFA supercell (i.e. take out the other atoms).  This greatly reduces 
      // the time of computation (don't need to scan through unnecessary atoms) 
      xstructure xstr_LFA_only=xstr2;
      uint num_atoms=xstr_LFA_only.atoms.size();
      for(uint i=0;i<num_atoms;i++){
	if(xstr_LFA_only.atoms[i].name!=lfa){
	  xstr_LFA_only.RemoveAtom(i);
	  num_atoms--;
	  i--;
	}
      }
      //cerr << "LFA SUPERELL" << endl;
      xstructure xstr_LFA_supercell=GetSuperCell(xstr_LFA_only,3,0,0,0,3,0,0,0,3);

      // Determines the number of LFAs in the supercell.
      int num_LFAs=-1; //-1 as default value 
      for(uint q=0; q<xstr.num_each_type.size();q++){
	if(xstr.species[q] == lfa){ 
	  num_LFAs= xstr.num_each_type[q];
	  if(LDEBUG){cerr << "compare:: " << "Number of LFAs in supercell: " << xstr.species[q] << "= " << num_LFAs << endl;}
	}
      }

      // === THREAD PREPARATION FOR PARALLEL PROCESSING OF QUADRUPLET SEARCH === //

      //[THREADS]// DECLARATION OF ATOMIC BOOL SECTION: Allows the threads to communicate. 
      //[THREADS]// This is useful for stopping the threads if the misfit falls below
      //[THREADS]// the compatible misfit criterion (mis<0.1) in any of the threads. 
      //[THREADS]std::atomic_bool misfit_in_threshold_found (false);

      //[NONTHREADS]bool misfit_in_threshold_found=false;
   
      vector<xstructure> xstr1_for_thread;
      vector<double> possible_minMis;
      vector<vector<xstructure> > vvprotos;
      //[THREADS]vector<std::thread> threads;

      // compute xstr1 information once only (perhaps we can use this in the directory scheme) 
      // and really only calculate once; but not that expensive, may be more expensive to store --memory!)
      vector<double> D1,F1;
      cellDiagonal(xstr1_tmp,D1,F1,1);
      xstr1_tmp.lattice=GetClat(xstr1_tmp.a,xstr1_tmp.b,xstr1_tmp.c,xstr1_tmp.alpha,xstr1_tmp.beta,xstr1_tmp.gamma);
      for(uint iat=0; iat<xstr1_tmp.atoms.size(); iat++){
        xstr1_tmp.atoms[iat].cpos=F2C(xstr1_tmp.lattice,xstr1_tmp.atoms[iat].fpos);
      }
      vector<double> all_nn1 = computeNearestNeighbors(xstr1_tmp);
      // DX NEW
      for(uint n=0; n<num_proc; n++){
	vector<xstructure> vprotos_tmp;
	vvprotos.push_back(vprotos_tmp);
	xstr1_for_thread.push_back(xstr1_tmp);
	possible_minMis.push_back(1.0);
      }
      
      vector<xmatrix<double> > lattices;
      vector<xmatrix<double> > clattices;
      vector<vector<uint> > ij_index;
      vector<double> latt_devs;

      //[THREADS]vector<std::thread> threads0;
      vector<xvector<double> > translation_vectors;
      //cerr << "FINDING TRANSLATION VECTORS: " << endl;
      //[THREADS]for(uint n=0; n<1; n++){
      //[THREADS]	//cerr << "here: " << n <<  endl;
      //[THREADS]	threads0.push_back(std::thread(quadrupletSearch,q1,xstr_LFA_supercell,std::ref(translation_vectors),std::ref(ij_index)));
      //[THREADS]}        
      //[THREADS]for(uint t=0;t<threads0.size();t++){
      //[THREADS]	threads0[t].join();
      //[THREADS]}
      //NON-THREADED
      quadrupletSearch(q1,xstr_LFA_supercell,translation_vectors,ij_index);

      double abs_det_q1=abs(det(q1));
      xvector<double> abc_angles_q1=Getabc_angles(q1,DEGREES);

      buildSimilarLattices(translation_vectors, q1, abs_det_q1, abs_det_q1, abc_angles_q1, lattices, clattices, latt_devs, fast_match);
      if(LDEBUG){cerr << "pflow::threadGeneration: Number of lattices to compare: " << lattices.size() << endl;}

      if(lattices.size()>0){
	//[THREADS]vector<std::thread> threads1;
	vector<vector<xmatrix<double> > > lattices_split;
	vector<vector<xmatrix<double> > > clattices_split;
	vector<vector<double> > latt_devs_split;

	//[THREADS]uint num_per_thread = lattices.size()/num_proc;
	//[THREADS]uint residual = lattices.size()%num_proc;
	uint num_per_thread = lattices.size()/1; //NON-THREADED
	uint residual = lattices.size()%1; //NON-THREADED
	bool accounted_for_residual=false;
	if(residual!=0){
	  num_per_thread+=1;
	}
	uint count = 0;
	uint thread_count = 0;
	vector<xmatrix<double> > latt_tmp, clatt_tmp;
	vector<double> tmp_dev;
	for(uint l=0; l<lattices.size(); l++){
	  latt_tmp.push_back(lattices[l]); clatt_tmp.push_back(clattices[l]); tmp_dev.push_back(latt_devs[l]);
	  count+=1;
	  if(count == num_per_thread && thread_count<num_proc-1){
	    thread_count+=1;
	    lattices_split.push_back(latt_tmp);
	    clattices_split.push_back(clatt_tmp);
	    latt_devs_split.push_back(tmp_dev);
	    latt_tmp.clear(); clatt_tmp.clear(); tmp_dev.clear();
	    count = 0;
	  }
	  else if(thread_count==num_proc-1 && l==lattices.size()-1){
	    thread_count+=1;
	    lattices_split.push_back(latt_tmp);
	    clattices_split.push_back(clatt_tmp);
	    latt_devs_split.push_back(tmp_dev);
	    latt_tmp.clear(); clatt_tmp.clear(); tmp_dev.clear();
	    count = 0;
	  }
	  if(!accounted_for_residual && residual!=0 && thread_count==residual){
	    accounted_for_residual=true;
	    num_per_thread=num_per_thread-1;
	  }
	}
	//[THREADS]uint recovered=0;
	//[THREADS]//Need the following safety in case the number of threads is greater than the number of lattices to test
	//[THREADS]uint num_of_threads=0;
	//[THREADS]if(lattices_split.size()>=num_proc){
	//[THREADS]  num_of_threads=num_proc;
	//[THREADS]}
	//[THREADS]else if(lattices_split.size()<num_proc){
	//[THREADS]  num_of_threads=lattices_split.size();
	//[THREADS]}
	//[THREADS]for(uint n=0; n<num_of_threads; n++){
	//[THREADS]  for(uint h=0;h<lattices_split[n].size();h++){
	//[THREADS]    recovered+=1;
	//[THREADS]    //cerr << "recovered: " << recovered << " - " << lattices_split[n][h] << endl;
	//[THREADS]  }
	//[THREADS]} 
	//[THREADS]if(recovered != lattices.size()){
	//[THREADS]  cerr << "The splitting of jobs failed...not all were accounted for: " << recovered << " != " << lattices.size() << endl;
	//[THREADS]  exit(1);
	//[THREADS]}
        
	//[THREADS]for(uint n=0; n<num_of_threads; n++){
	//[THREADS]  threads1.push_back(std::thread(structureSearch,lfa,all_nn1,xstr,
	//[THREADS]		    std::ref(vvprotos[n]),std::ref(xstr1_for_thread[n]),xstr2,type_match,std::ref(possible_minMis[n]),
	//[THREADS]		    std::ref(lattices_split[n]),std::ref(clattices_split[n]),std::ref(latt_devs_split[n]),
        //[THREADS]                    fast_match));
	//[THREADS]}         
	//[THREADS]for(uint t=0;t<threads1.size();t++){
	//[THREADS]  threads1[t].join();
	//[THREADS]}
        uint n=0;
	structureSearch(lfa,all_nn1,xstr,vvprotos[n],xstr1_for_thread[n],xstr2,type_match,possible_minMis[n],
	   	        lattices_split[n],clattices_split[n],latt_devs_split[n],fast_match);
      }
      //cerr << "========== possible_minMis.size(): " << possible_minMis.size() << endl;
      for(uint p=0;p<possible_minMis.size();p++){
	if(p==0 && y==0 && x==0){ // DX 2/8/17 - need to add x==0 ortherwise matches can be overwritten
	  minMis=possible_minMis[p];
	  xstr1=xstr1_for_thread[p];
	  vprotos=vvprotos[p];
	}
	else{
	  if(possible_minMis[p]<=minMis){
	    minMis=possible_minMis[p];
	    xstr1=xstr1_for_thread[p];
	    vprotos=vvprotos[p];
	  }
	}
        //cerr << "minMis: " << minMis << endl;
      }
      //break if(minMis<=0.1) break;
    } 
    } 
  }  
}

// ***************************************************************************
// Quadruplet Search
// ***************************************************************************
namespace compare{
  void quadrupletSearch(const xmatrix<double>& q1, const xstructure& xstr_LFA_supercell, 
			vector<xvector<double> >& lattice_vecs, vector<vector<uint> >& ij_index){

    // This function scans through the possible quadruplets (sets of 4 LFA atoms) i
    // to find a lattice which is commensurate with the reference structure (xstr1). 
    // This function is parallelized since it is the time-limiting function.

    bool LDEBUG=(false || XHOST.DEBUG);

    double min_q1_a = aurostd::modulus(q1(1))-aurostd::modulus(q1(1))*0.1;
    double max_q1_a = aurostd::modulus(q1(1))+aurostd::modulus(q1(1))*0.1;
    double min_q1_b = aurostd::modulus(q1(2))-aurostd::modulus(q1(2))*0.1;
    double max_q1_b = aurostd::modulus(q1(2))+aurostd::modulus(q1(2))*0.1;
    double min_q1_c = aurostd::modulus(q1(3))-aurostd::modulus(q1(3))*0.1;
    double max_q1_c = aurostd::modulus(q1(3))+aurostd::modulus(q1(3))*0.1;

    if(LDEBUG){
      cerr << "compare::quadrupletSearch: Lattice parameters: " << aurostd::modulus(q1(1)) << ", " << aurostd::modulus(q1(2)) << ", " << aurostd::modulus(q1(3)) << endl;
      cerr << "compare::quadrupletSearch: Modulus search range for lattice vector a: " << min_q1_a << " - " << max_q1_a << endl;
      cerr << "compare::quadrupletSearch: Modulus search range for lattice vector b: " << min_q1_b << " - " << max_q1_b << endl;
      cerr << "compare::quadrupletSearch: Modulus search range for lattice vector c: " << min_q1_c << " - " << max_q1_c << endl;
    }

    xvector<double> tmp_vec;
    double tmp_mod = 0.0;
    
    // Search all possible vectors with modulus comparable to one of the lattice vectors 
    for(uint i=0; i<xstr_LFA_supercell.atoms.size(); i++){
      for(uint j=i+1; j<xstr_LFA_supercell.atoms.size(); j++){ //upper triangular
	tmp_vec = xstr_LFA_supercell.atoms[j].cpos-xstr_LFA_supercell.atoms[i].cpos;
        tmp_mod = aurostd::modulus(tmp_vec);
        if((tmp_mod <= max_q1_a && tmp_mod >= min_q1_a) || 
           (tmp_mod <= max_q1_b && tmp_mod >= min_q1_b) || 
           (tmp_mod <= max_q1_c && tmp_mod >= min_q1_c)){ 
          bool vec_stored = false;
          for(uint p=0;p<lattice_vecs.size();p++){
            if(identical(lattice_vecs[p],tmp_vec,1e-10)){
              vec_stored = true;
              break;
            }
          }       
          if(vec_stored == false){
            lattice_vecs.push_back(tmp_vec);
            // Store indices of atoms comprising the vector
            vector<uint> ij;
            ij.push_back(i); ij.push_back(j);
            ij_index.push_back(ij); 
            // Store negative (may not be needed)
            //lattice_vecs.push_back(-tmp_vec);
            //vector<uint> ji;
            //ji.push_back(j); ji.push_back(i);
            //ij_index.push_back(ji); 
          }
        }
      }
    }

    if(LDEBUG){
      cerr << "compare::quadrupletSearch: Number of potential lattice vectors: " << lattice_vecs.size() << endl;
    } 
    // Removing non-periodic lattice vectors
    vector<xvector<double> > lattice_vecs_periodic;
    for(uint i=0;i<lattice_vecs.size();i++){
      if(vectorPeriodic(lattice_vecs[i],xstr_LFA_supercell,ij_index[i][0],ij_index[i][1])){
        lattice_vecs_periodic.push_back(lattice_vecs[i]);
        lattice_vecs_periodic.push_back(-lattice_vecs[i]);
      }
    }
    lattice_vecs = lattice_vecs_periodic;
    if(LDEBUG){
      cerr << "compare::quadrupletSearch: Number of lattice vectors (preserves periodicity): " << lattice_vecs.size() << endl;
      for(uint i=0;i<lattice_vecs.size();i++){
        cerr << "compare::quadrupletSearch: lattice vector " << i << ": " << lattice_vecs[i] << " (" << aurostd::modulus(lattice_vecs[i]) << ")" << endl; 
      }
    } 
  }
}

// ***************************************************************************
// Build All Lattices
// ***************************************************************************
namespace compare{
  bool buildSimilarLattices(vector<xvector<double> >& translation_vectors, xmatrix<double>& q1, double& xstr1_vol, double& abs_det_q1, 
                            xvector<double>& abc_angles_q1, vector<xmatrix<double> >& lattices, vector<xmatrix<double> >& clattices, 
                            vector<double>& latt_devs, const bool& fast_match){

    bool LDEBUG=(false || XHOST.DEBUG);

    vector<double> D1,F1;
    cellDiagonal(q1,D1,F1,1);

    double tol_vol=0.1;
    double det_tol=tol_vol*abs_det_q1;

    double tol_a=abc_angles_q1[1]*0.3;
    double tol_b=abc_angles_q1[2]*0.3;
    double tol_c=abc_angles_q1[3]*0.3;

    if(LDEBUG){
      cerr << "compare::buildSimilarLattices: Tolerance for a (Angstroms): " << tol_a << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for b (Angstroms): " << tol_b << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for c (Angstroms): " << tol_c << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for alpha (degrees): " << abc_angles_q1[4]*0.3 << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for beta (degrees): " << abc_angles_q1[5]*0.3 << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for gamma (degrees): " << abc_angles_q1[6]*0.3 << endl;
    }

    xmatrix<double> tmp_lattice(3,3);
    xmatrix<double> tmp_clatt(3,3);

    int store=0;
    vector<double> translations_mod;
    for(uint i=0;i<translation_vectors.size();i++){
      translations_mod.push_back(aurostd::modulus(translation_vectors[i]));
    }
    if(LDEBUG){
      cerr << "buildSimilarLattices:: Number of lattice vectors: " << translation_vectors.size() << endl;
    }

    // Build all possible unit cells with combinations of lattice vectors (order matters, hence not upper triangular)
    for(uint i=0;i<translation_vectors.size();i++){
      if(abs(translations_mod[i]-abc_angles_q1[1])<tol_a){ //check a
	for(uint j=0;j<translation_vectors.size();j++){
	  if(j!=i){
	    if(abs(translations_mod[j]-abc_angles_q1[2])<tol_b){ // check b
	      for(uint k=0;k<translation_vectors.size();k++){
		if(k!=i && k!=j){
		  if(abs(translations_mod[k]-abc_angles_q1[3])<tol_c){ //check c
		    tmp_lattice = SYM::xvec2xmat(translation_vectors[i],translation_vectors[j],translation_vectors[k]);         
		    if(abs(abs_det_q1-abs(det(tmp_lattice))) < det_tol){ //check determinant/volume
		      xvector<double> abc_angles_q2=Getabc_angles(tmp_lattice,DEGREES);
		      if(checkTolerance(abc_angles_q1,abc_angles_q2)==false){
			tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
			bool unique = true;
			for(uint t=0;t<lattices.size();t++){
			  if(identical(tmp_lattice,lattices[t],1e-10)){
			    // DX TEST - CANNOT DO: Eliminates potential matches for(uint t=0;t<clattices.size();t++){
			    // DX TEST - CANNOT DO: Eliminates potential matches if(identical(tmp_clatt,clattices[t],1e-10)){
			    unique=false;
			    break;
			  }
			}
			double tmp_latt_dev = checkLatticeDeviation(xstr1_vol,tmp_lattice,D1,F1);
                        // Time-saver, keep lattices that have a deviation smaller than Burzlaff's same-family requirement
                        // otherwise, there is no possible way that it could match with anything or be in the same-family
			if(unique && fast_match && tmp_latt_dev <= 0.1){ //fast match doesn't care about finding same family information
			  lattices.push_back(tmp_lattice); // stores original original orientation
			  clattices.push_back(tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
			  latt_devs.push_back(tmp_latt_dev);
			  store++;
                        }
			else if(unique && !fast_match && tmp_latt_dev <= 0.2){
			  lattices.push_back(tmp_lattice); // stores original original orientation
			  clattices.push_back(tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
			  latt_devs.push_back(tmp_latt_dev);
			  store++;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// checkLatticeDeviation
// ***************************************************************************
namespace compare{
  double checkLatticeDeviation(double& xstr1_vol, xmatrix<double>& q2,vector<double>& D1,vector<double>& F1){
    double scale=xstr1_vol/(aurostd::abs(aurostd::det(q2)));
    scale=pow(scale,0.3333);
    vector<double> D2,F2;
    cellDiagonal(q2,D2,F2,scale);
    double latt_dev=latticeDeviation(D1,D2,F1,F2);
    return latt_dev;
  }
}

// ***************************************************************************
// Internal structure
// ***************************************************************************
namespace compare{
  bool structureSearch(const string& lfa, 
                        const vector<double>& all_nn1, 
			const xstructure& xstr, 
			vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2, 
			const int& type_match, double& possible_minMis,
			vector<xmatrix<double> >& lattices,
			vector<xmatrix<double> >& clattices, 
			vector<double>& latt_devs, 
			const bool& fast_match){ 


    bool LDEBUG=(false || XHOST.DEBUG);

    double mis=1;  
    xstructure proto;
    int flag=0;
    xstructure xstr2_tmp = xstr2;

    for(uint p=0;p<lattices.size();p++){
      if(LDEBUG){
        cerr << "compare::structureSearch: Trying lattice " << p << endl;
      }
      proto=xstr;
      proto.lattice=lattices[p];

      // Transform
      for(uint iat=0;iat<proto.atoms.size();iat++){
	proto.atoms[iat].fpos=C2F(proto.lattice,proto.atoms[iat].cpos);
      }
      proto.lattice=clattices[p];
	for(uint iat=0;iat<proto.atoms.size();iat++){
	proto.atoms[iat].cpos=F2C(proto.lattice,proto.atoms[iat].fpos);
      }
      xstructure proto_new;
      proto_new.title=proto.title;
      proto_new.lattice=clattices[p];

      // DX NEW - START =======================
      xmatrix<double> f2c = trasp(proto.lattice);
      xmatrix<double> c2f = aurostd::inverse(trasp(proto.lattice));
      bool skew = false;
      double tol=0.01;
      deque<_atom> new_basis;
      for(uint j=0;j<proto.atoms.size();j++){
	if(new_basis.size()==0){
	  proto.atoms[j].fpos = BringInCell(proto.atoms[j].fpos,1e-10);
	  proto.atoms[j].cpos = f2c*proto.atoms[j].fpos;
	  new_basis.push_back(proto.atoms[j]);
	  //proto_new.AddAtom(proto.atoms[j]);
	}
	else{
	  bool duplicate_lattice_point=false;
	  for(uint a=0; a<new_basis.size(); a++){
	    xvector<double> tmp = BringInCell(proto.atoms[j].fpos,1e-10);
	    if(SYM::MapAtom(new_basis[a].fpos,tmp,c2f,f2c,skew,tol)){
	      duplicate_lattice_point=true;
	      break;
	    }
	  }
	  if(duplicate_lattice_point==false){
	    proto.atoms[j].fpos = BringInCell(proto.atoms[j].fpos,1e-10);
	    proto.atoms[j].cpos = f2c*proto.atoms[j].fpos;
	    new_basis.push_back(proto.atoms[j]);
	    //proto_new.AddAtom(proto.atoms[j]);
	  }
	}
      }
      proto_new.atoms = new_basis;
      proto_new.BringInCell(1e-10); 
      proto_new.FixLattices();
      proto_new.SpeciesPutAlphabetic();
      deque<int> sizes = SYM::arrange_atoms(new_basis);
      proto_new = pflow::SetNumEachType(proto_new, sizes);
      proto = proto_new;
      if(sameSpecies(proto,xstr1,false)){
	vector<double> all_nn_proto;
	bool all_nn_calculated = false;
	for(uint iat=0; iat<proto.atoms.size();iat++){
	  if(proto.atoms[iat].name==lfa){
	    proto.ShifOriginToAtom(iat);
	    proto.BringInCell(1e-10);
	    vector<uint> im1, im2;
	    vector<double> min_dists;
	    if(findMatch(xstr1,proto,im1,im2,min_dists,type_match)){;
	      double cd, f;
	      // Only calculate the NN for the proto if we found suitable matches.  
	      // Only calculate once, nothing changes between shifts to origin (affine)
	      if(!all_nn_calculated){
                all_nn_proto = computeNearestNeighbors(proto);
		if(LDEBUG){
		  cerr << "compare::structureSearch: Nearest neighbors:" << endl;
		  for(uint a=0;a<all_nn_proto.size();a++){
		    cerr << "compare::structureSearch: Nearest neighbor distance from " << a << " atom: " << all_nn_proto[a] << endl;
		  }
		}
		all_nn_calculated = true;
              }
	      coordinateDeviation(xstr1,proto,all_nn1,all_nn_proto,im1,im2,min_dists,cd,f);
	      mis=computeMisfit(latt_devs[p],cd,f);
	      //if(LDEBUG){
	      //  cerr << "mis,latt_dev,cd,f: " << mis << ", " <<latt_devs[p] << ", " << cd << ", " << f <<  endl;
	      //}
	      if(flag==0){
	        flag=1;
		//cerr << "storing: " << proto << endl;
		vprotos.push_back(proto);
		possible_minMis=mis;
	      }
	      else{
	        if(mis<possible_minMis){
		  vprotos.clear();
		  possible_minMis=mis;
		  //cerr << "storing: " << proto << endl;
		  vprotos.push_back(proto); //to here
		}
	      }
              // If we want to simply find a match and not find the best match, we can exit early
	      if(mis<0.1 && fast_match) {
                return true;
	        //DEBUGGING
		//cerr <<"Winning combo: "<<i<<","<<j<<","<<k<<","<<w<<endl;
		//cerr << "proto.lattice: " << proto.lattice << endl;
		//cerr << "lattice(1): " << modulus(proto.lattice(1)) << endl;
		//cerr << "lattice(2): " << modulus(proto.lattice(2)) << endl;
		//cerr << "lattice(3): " << modulus(proto.lattice(3)) << endl;
	      }
	    }
	  }
	}
      }// end of if protos.size()...
    }
    return true;
  }
}

//---------------------------------------------------------------

// ***************************************************************************
//                                 END - FINE
// 			AFLOW Compare Structure - Functions
//		David Hicks (d.hicks@duke.edu) and Carlo de Santo
// ***************************************************************************
