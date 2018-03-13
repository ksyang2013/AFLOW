// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2013               *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_kesong_hnfcell.cpp
// functions written by KESONG YANG
// 2010-2013: kesong.yang@gmail.com

#ifndef _AFLOW_CONTRIB_KESONG_HNFCELL_CPP_
#define _AFLOW_CONTRIB_KESONG_HNFCELL_CPP_

#include "aflow_contrib_kesong.h"
#include "aflow.h"

bool CALCULATE_ENERGY_UFF = true;

// ***************************************************************************
// void combine(vector<int> &range, vector<int> &cur, vector<vector<int> > &final_result, int start, int depth)
// ***************************************************************************
// subroutine of combine
void combine(vector<int> &range, vector<int> &cur, vector<vector<int> > &final_result, int start, int depth) {
  if(depth == (int) cur.size()) {
    vector<int> tmp;
    for (vector<int>::iterator i = cur.begin(); i != cur.end(); ++i) {
      tmp.push_back(*i);
    }   
    final_result.push_back(tmp);
    return;
  }   

  for (int i = start; i != (int) range.size(); ++i) {
    cur[depth] = range[i];
    combine(range, cur, final_result, i + 1, depth + 1); 
  }   
}

// ***************************************************************************
// void combine(vector<int> &range, vector<vector<int> > &final_result, int n)
// ***************************************************************************
void combine(vector<int> &range, vector<vector<int> > &final_result, int n) {
  //Combine function, choose all the combinations of n numbers from range
  if(n > (int) range.size()) {
    cerr << "Error in combination! " << endl;
    exit(0);
  }   
  vector<int> result(n);
  combine(range, result, final_result, 0, 0); 
}

// ***************************************************************************
// vector<xmatrix<double> > CalculateHNF(xstructure a, int n)
// ***************************************************************************
vector<xmatrix<double> > CalculateHNF(xstructure str, int n) {
  //xstructure str=a; //CO, already copies
  str.BringInCell();

  // GET SYMMETRY
  str.FixLattices();
  xstructure non_pocc_xstr=NormalizeXstructure(str);
  non_pocc_xstr.CalculateSymmetryPointGroup();
  //str.CalculateSymmetryPointGroup();  //CO, cannot calculate symmetry of POCC structure as atoms are on top of each other
  // str.CalculateSymmetryFactorGroup();

  // PHYSICAL REVIEW B 77, 224115 2008
  // Algorithm for generating derivative structures
  // Gus L. W. Hart and Rodney W. Forcade
  xmatrix<double> HNF(3,3);
  vector<xmatrix<double> > vHNF;
  xmatrix<double> A(3,3),B(3,3),Bi(3,3),Bj(3,3),R(3,3),H(3,3);A=trasp(str.lattice);
  vector<xmatrix<double> > vB;    // contains lattices of the potential xstructures
  vector<xmatrix<double> > sHNF;  // contains only the supercells HNF that are uniques
  long long int i=0;
  bool found=FALSE;
  for(int a=1;a<=n;a++) {
    for(int c=1;c<=n/a;c++) {
      for(int f=1;f<=n/a/c;f++) {
	if(a*c*f==n) {
	  for(int b=0;b<c;b++) {
	    for(int d=0;d<f;d++) {
	      for(int e=0;e<f;e++) {
		i++;
		HNF(1,1)=a;HNF(1,2)=0;HNF(1,3)=0;
		HNF(2,1)=b;HNF(2,2)=c;HNF(2,3)=0;
		HNF(3,1)=d;HNF(3,2)=e;HNF(3,3)=f;
		vHNF.push_back(HNF);
	      }
        }
      }
    }
      }
    }
  }

  for(i=0;i<(int) vHNF.size();i++) {
    Bi=A*vHNF.at(i);
    found=FALSE;
    for(uint istr=0;istr<vB.size()&&found==FALSE;istr++)  {                              // cycle through the unique vBs
      Bj=vB.at(istr);
      for(uint pgroup=0;pgroup<non_pocc_xstr.pgroup.size()&&found==FALSE;pgroup++) {                  // cycle through the pgroup of str
	R=trasp(non_pocc_xstr.pgroup.at(pgroup).Uc); 
	H=aurostd::inverse(Bj)*aurostd::inverse(R)*Bi;
	if(aurostd::isinteger(H)) {found=TRUE;}
      }
    }
    if(found==FALSE) { // not found, then plug
      vB.push_back(Bi);
      sHNF.push_back(vHNF.at(i));
    }
  }
  //debug 
  return sHNF;
}

// ***************************************************************************
// xstructure XstrSubstitute(xstructure &a, int n, string b)
// ***************************************************************************
xstructure XstrSubstitute(xstructure &a, int n, string b) {
  xstructure xstrb=a;

  _atom atom_replaced=xstrb.atoms.at(n); //changed the name of the substituted atom
  atom_replaced.name=b;

  xstrb.RemoveAtom(n);    //Remove this atom
  xstrb.AddAtom_POCC(atom_replaced); //Add the substituted atom, note the added atom always appear in the last 
  return xstrb;
}

// ***************************************************************************
// bool sort_function_struct_num_name (atom_number_name i, atom_number_name j)
// ***************************************************************************
bool sort_function_struct_num_name (atom_number_name i, atom_number_name j) {
  return (i.number>j.number);
}

// ***************************************************************************
// vector<atom_number_name> SortMax2Min(vector<atom_number_name>& a)
// ***************************************************************************
vector<atom_number_name> SortMax2Min(vector<atom_number_name>& a) {
  sort(a.begin(), a.end(), sort_function_struct_num_name);
  return a;
}

// ***************************************************************************
// bool myfunction (int i,int j) 
// ***************************************************************************
bool myfunction_int (int i,int j) 
{return (i>j);}

// ***************************************************************************
// vector<int> SortMax2Min(vector<int> a)
// ***************************************************************************
vector<int> SortMax2Min(vector<int>& a) {
  //sort from max to min
  sort(a.begin(), a.end(), myfunction_int);
  return a;
}

// ***************************************************************************
// bool myfunction_double (double i,double j) 
// ***************************************************************************
bool myfunction_double (double i,double j) {
  return (i>j);
}

// ***************************************************************************
// vector<double> SortMax2Min(vector<double>& a)
// ***************************************************************************
vector<double> SortMax2Min(vector<double>& a) {
  //sort from max to min
  sort(a.begin(), a.end(), myfunction_double);
  return a;
}

// ***************************************************************************
// xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, string b)
// ***************************************************************************
xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, string b) {
  xstructure xstr_tmp=a;
  vector<int> vec_n_sorted = SortMax2Min(vec_n);

  //Firstly, we add all the atoms we need
  for (uint i=0; i<vec_n_sorted.size();i++) {
    int k= vec_n_sorted.at(i);
    _atom atom_replaced = xstr_tmp.atoms.at(k);
    atom_replaced.name = b;
    xstr_tmp.AddAtom_POCC(atom_replaced);
  }
  //Now we remove all the atoms that we do not need
  for (uint i=0; i<vec_n_sorted.size();i++) {
    int k= vec_n_sorted.at(i);
    xstr_tmp.RemoveAtom(k);
  }
  int ntype = xstr_tmp.num_each_type.size();
  xstr_tmp.SpeciesSwap(0,ntype-1);   //Force to re-assign
  xstr_tmp.SpeciesPutAlphabetic();
  xstr_tmp.MakeBasis();
  return xstr_tmp;
}

// ***************************************************************************
// xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, vector<string> b)
// ***************************************************************************
xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, vector<string> b) {
  if(vec_n.size()!=b.size()) { 
    cerr <<"Error! vec_n.size()!=b.size() !" << endl;
    exit(1);
  }
  vector<atom_number_name> vec_num_name;
  atom_number_name tmp;
  for (uint i=0; i<vec_n.size();i++) {
    tmp.number = vec_n.at(i);
    tmp.name = b.at(i);
    vec_num_name.push_back(tmp);
  }
  vector<atom_number_name> vec_num_name_new = SortMax2Min(vec_num_name);
  vector<int> vec_num_sorted;
  vector<string> vec_str_sorted;
  for (uint i=0; i<vec_num_name_new.size();i++) {
    vec_num_sorted.push_back(vec_num_name_new.at(i).number);
    vec_str_sorted.push_back(vec_num_name_new.at(i).name);
  }

  xstructure xstr_tmp=a;
  //Firstly, we add all the atoms we need
  for (uint i=0; i<vec_num_sorted.size();i++) {
    int k= vec_num_sorted.at(i);
    _atom atom_replaced = xstr_tmp.atoms.at(k);
    atom_replaced.name = vec_str_sorted.at(i);
    xstr_tmp.AddAtom_POCC(atom_replaced);
  }
  //Now we remove all the atoms that we do not need
  for (uint i=0; i<vec_num_sorted.size();i++) {
    int k= vec_num_sorted.at(i);
    xstr_tmp.RemoveAtom(k);
  }
  xstr_tmp.SpeciesPutAlphabetic();
  int ntype = xstr_tmp.num_each_type.size();
  xstr_tmp.SpeciesSwap(0,ntype-1);   //Force to re-assign
  xstr_tmp.SpeciesPutAlphabetic();
  xstr_tmp.MakeBasis();

  //Clean Vacancies
  int natoms = xstr_tmp.atoms.size();
  if(xstr_tmp.atoms.at(natoms-1).name == name_vacancy) {
    return CleanVacancy(xstr_tmp);
  }
  return xstr_tmp;
}

// ***************************************************************************
// xstructure NormalizeXstructure(xstructure a)
// ***************************************************************************
xstructure NormalizeXstructure(xstructure a) {
  xstructure b=a;
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(a);

  vector<int> removed_number; 
  for(uint i=0; i<NormalisedNumber.size(); i++) {
    if(NormalisedNumber.at(i).size() >= 2) {
      for (uint j=1; j<NormalisedNumber.at(i).size();j++) { //Keep the first one
	int k=NormalisedNumber.at(i).at(j);
	removed_number.push_back(k);
      }
    }
  }
  vector<int> sorted_removed_number = SortMax2Min(removed_number);
  for (uint i=0; i<sorted_removed_number.size();i++) {
    int k=sorted_removed_number.at(i);
    b.RemoveAtom(k);
  }
  return b;
}

// ***************************************************************************
// vector<vector<int> > CalculateXstrNumberSupercell(xstructure a, int n)
// ***************************************************************************
vector<vector<int> > CalculateXstrNumberSupercell(xstructure a, int n) {
  //This function is used to produce the coresponding number in the supercell for the special atom in xstructure a
  vector<vector<int> > number_supercell;
  vector<int> number_tmp;

  for (int i=0; i< (int) a.atoms.size(); i++) {
    number_tmp.clear();
    for (int j=i*n; j<(i+1)*n; j++) {
      number_tmp.push_back(j);
    }
    number_supercell.push_back(number_tmp);
  }
  return number_supercell;
}

// ***************************************************************************
// int hnf_double2int(double a)
// ***************************************************************************
int hnf_double2int(double a) {   
  double adiff=a-int(a);
  if(adiff >=0 && adiff < 0.5 ) {
      return int(a);
  }
  else {
      return int(a)+1;
  }
}

// ***************************************************************************
// bool CheckDoubleOccupied<xstructure xstr_orig) {
// ***************************************************************************
bool CheckDoubleOccupied(xstructure xstr_orig) {
  //DoubleOccupied means two atoms occupy one position or one atom and one vacancy occupy one position
  bool FLAG_double_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==1 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_double_occupied = true;}
    if(num_norm_xstr.at(i).size()==2 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_double_occupied = true;}
  }
  return FLAG_double_occupied;
}

// ***************************************************************************
// bool CheckTripleOccupied(xstructure xstr_orig)
// ***************************************************************************
bool CheckTripleOccupied(xstructure xstr_orig) {
  //TripleOccupied means three atoms occupy one position or two atom and one vacancy occupy one position
  bool FLAG_triple_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==2 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_triple_occupied = true;}
    if(num_norm_xstr.at(i).size()==3 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_triple_occupied = true;}
  }
  return FLAG_triple_occupied;
}

// ***************************************************************************
// bool CheckFourfoldOccupied(xstructure xstr_orig)
// ***************************************************************************
bool CheckFourfoldOccupied(xstructure xstr_orig) {
  //FourfoldOccupied means three atoms occupy one position or two atom and one vacancy occupy one position
  bool FLAG_fourfold_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==3 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fourfold_occupied = true;}
    if(num_norm_xstr.at(i).size()==4 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fourfold_occupied = true;}
  }
  return FLAG_fourfold_occupied;
}

// ***************************************************************************
// bool CheckFivefoldOccupied(xstructure xstr_orig)
// ***************************************************************************
bool CheckFivefoldOccupied(xstructure xstr_orig) {
  //FourfoldOccupied means three atoms occupy one position or two atom and one vacancy occupy one position
  bool FLAG_fivefold_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==4 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fivefold_occupied = true;}
    if(num_norm_xstr.at(i).size()==5 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fivefold_occupied = true;}
  }
  return FLAG_fivefold_occupied;
}


// ***************************************************************************
// void CombineAll(const vector<vector<vector<int> > > &allVecs, size_t vecIndex, vector<vector<int> > &intSoFar, vector<vector<vector<int> > > &results)
// ***************************************************************************
void CombineAll(const vector<vector<vector<int> > > &allVecs, size_t vecIndex, vector<vector<int> > &intSoFar, vector<vector<vector<int> > > &results) {
  //function to combine all the vectors (choose one from each one)
  //Call method: CombineAll(d, 0, e, f); d initial vector; e temporary vector, f final result
  if(vecIndex >= allVecs.size())
    {   
      vector<vector<int> > str_tmp;
      for(uint i=0;i<intSoFar.size();i++) { 
	str_tmp.push_back(intSoFar.at(i));
      }   
      results.push_back(str_tmp);
      return;
    }   
  for (size_t i=0; i<allVecs.at(vecIndex).size(); i++) {
    vector<vector<int> > tmp=intSoFar;
    tmp.push_back(allVecs.at(vecIndex).at(i));
    CombineAll(allVecs, vecIndex+1, tmp, results);
  }   
}

// ***************************************************************************
// vector<string> CalculateSecondAtomicNameSupercell(xstructure xstr_orig, int n)
// ***************************************************************************
vector<string> CalculateSecondAtomicNameSupercell(xstructure xstr_orig, int n) {
  vector<string> AtomicNameSupercell;
  //    if(CheckDoubleOccupied(xstr_orig)) {
  xstructure xstr_norm = NormalizeXstructure(xstr_orig); //Format POSCAR
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig); //Format normalized POSCAR's number using 2D vector
  vector<vector<int> > vec_num_supercell=CalculateXstrNumberSupercell(xstr_norm,n); //Produce supercell's number using 2D vector

  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if(xstr_norm.atoms.at(i).partial_occupation_flag) {
      //if this is partial occupation
      if(num_norm_xstr.at(i).size()==1) {
	//This means this site only contains vacancy
	for (uint j=0;j<vec_num_supercell.at(i).size();j++) {
	  AtomicNameSupercell.push_back(name_vacancy);
	}
      }
      else{
	//This means site site contains another element
	int int_tmp = num_norm_xstr.at(i).at(1); //Read the atomic number in xstr_orig
	string tmp_second_element_name = xstr_orig.atoms.at(int_tmp).name; //Read the name
	for (uint j=0;j<vec_num_supercell.at(i).size();j++) {
	  AtomicNameSupercell.push_back(tmp_second_element_name);
	}
      }
    }
    else {
      for (uint j=0;j<vec_num_supercell.at(i).size();j++) {
	string tmp_name=xstr_norm.atoms.at(i).name;
	AtomicNameSupercell.push_back(tmp_name);
      }
    }
  }
  return AtomicNameSupercell;
}


// ***************************************************************************
// vector<vector<int> > GenerateSupercellAtomNumber(xstructure xstr_orig, int n)
// ***************************************************************************
vector<vector<int> > GenerateSupercellAtomNumber(xstructure xstr_orig, int n) {
  //This function generates numbers of all the atoms which will be replaced
  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<vector<int> > > vec_two_num_vacancy;

  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig); //format like 0 1 2; 3

  xstructure xstr_norm = NormalizeXstructure(xstr_orig);
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, n); //format like 0 1 2 3; 4 5 6 7
  //    vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_orig, n); //format like 0 1 2 3; 4 5 6 7
  //    We just need to store the atoms' number using normalized xstructure

  //if(CheckDoubleOccupied(xstr_orig)) {
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if(xstr_norm.atoms.at(i).partial_occupation_flag) {
      int num_vacancy= hnf_double2int((1-xstr_norm.atoms.at(i).partial_occupation_value)*n);
      vector<int> vec_num_part=atom_num_supercell.at(i);
      vector<vector<int> > num_com_vacancy;
      num_com_vacancy.clear();
      combine(vec_num_part, num_com_vacancy, num_vacancy); 
      //generate all the combinations and store them into 2D vectors num_com_vacancy
      vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
      //push the 2D vector into 3D vector according to the atoms' number
    }
  }

  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  //choose only one combination (vector) from each atom's site, and recombine them as one entire combination (2D vector) for one structure
  //All the possible combinations are stored as 3D vector final_num_result
  vector<vector<int> > formated_final_num;

  //Extract the actural atom's number from the 3D vector, and store them into 2D vector formated_final_num
  vector<int> tmp; 
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
	tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }
  return formated_final_num;
}


// ***************************************************************************
// vector<xstructure> AssignPartialOccupation(xstructure xstr_supercell, xstructure xstr_orig, int n)
// ***************************************************************************
vector<xstructure> AssignPartialOccupation(xstructure xstr_supercell, xstructure xstr_orig, int n) {
  //This function use two of the most important function GenerateSupercellAtomNumber(xstr_orig, n) & CalculateSecondAtomicNameSupercell(xstr_orig, n) 

  vector<xstructure> xstr_final_list_supercell;

  //Produce atom number which needs to be replaced
  vector<vector<int> > tmp_atom_num_supercell = GenerateSupercellAtomNumber(xstr_orig, n);
  //DEBUG
  //Produce coressponding atomic names which will replace original atom
  vector<string> name_str=CalculateSecondAtomicNameSupercell(xstr_orig, n);

  for (uint i=0; i<tmp_atom_num_supercell.size();i++) {
    xstructure xstr_tmp=xstr_supercell;
    vector<int> atom_number = tmp_atom_num_supercell.at(i);
    vector<int> atom_number_sorted = SortMax2Min(tmp_atom_num_supercell.at(i));

    //Add all the needed atoms 
    for (uint i=0; i<atom_number_sorted.size();i++) {
      int k=atom_number_sorted.at(i);
      _atom atom_replaced = xstr_tmp.atoms.at(k);
      atom_replaced.name = name_str.at(k);
      xstr_tmp.AddAtom_POCC(atom_replaced);
    }

    //Remove all the additional atoms
    for (uint i=0; i<atom_number_sorted.size();i++) {
      int k=atom_number_sorted.at(i);
      xstr_tmp.RemoveAtom(k);
    }

    //Force to re-assign the structure 
    int ntype = xstr_tmp.num_each_type.size();
    xstr_tmp.SpeciesSwap(0,ntype-1);   //Force to re-assign
    xstr_tmp.SpeciesPutAlphabetic();
    xstr_tmp.MakeBasis();

    int LastAtom = xstr_tmp.atoms.size()-1;
    if(xstr_tmp.atoms.at(LastAtom).name == name_vacancy) { //this is safe 
      xstr_final_list_supercell.push_back(CleanVacancy(xstr_tmp));
    }
    else {
      xstr_final_list_supercell.push_back(xstr_tmp);
    }
  }
  return xstr_final_list_supercell;
}

// ***************************************************************************
// int CalculateSupercellSize(xstructure xstr_orig)
// ***************************************************************************
int CalculateSupercellSize(xstructure xstr_orig) {
  int n;
  if(!xstr_orig.partial_occupation_HNF) {
    cerr << "partial_occupation_HNF is not found!" << endl;
    cerr << "Calculating supercell size using occupation value ..." << endl; 
    vector<int> vec_denom;
    for (uint i=0; i<xstr_orig.atoms.size();i++) {
      double partial_value_tmp = xstr_orig.atoms.at(i).partial_occupation_value;
      str_num_data tmp;
      tmp = double2str_num_data(partial_value_tmp);
      vec_denom.push_back(tmp.denom);
    }
    n = CalculateLcmofVector(vec_denom);
  }
  else {
    n = xstr_orig.partial_occupation_HNF;
  }
  return n;
}

// ***************************************************************************
// xstructure AssignAlphabeticNameXstr(xstructure& xstr)
// ***************************************************************************
xstructure AssignAlphabeticNameXstr(xstructure& xstr) {
  vector<string> names;
  names.push_back("As"); names.push_back("Ba"); names.push_back("Cl"); names.push_back("Dy");
  names.push_back("Eu"); names.push_back("Fe"); names.push_back("Ge"); names.push_back("He");

  xstructure xstr_new = xstr;
  int iatom=0;
  for (uint itype=0; itype<xstr.num_each_type.size();itype++) {
    string species = string(names.at(xstr_new.atoms.at(iatom).type));
    xstr_new.species.at(itype)=species;
    for (int j=0;j<xstr.num_each_type.at(itype);j++) {
      xstr_new.atoms.at(iatom).name=species;
      xstr_new.atoms.at(iatom).CleanName();
      xstr_new.atoms.at(iatom).CleanSpin();
      xstr_new.atoms.at(iatom).name_is_given=TRUE;
      iatom++;
    }
  }
  return xstr_new;
}

// ***************************************************************************
// xstructure CleanNameXstr(xstructure& xstr)
// ***************************************************************************
xstructure CleanNameXstr(xstructure& xstr) {
  xstructure xstr_new = xstr;
  for (uint i=0; i<xstr_new.atoms.size();i++) {
    xstr_new.atoms.at(i).name_is_given=false;
  }
  return xstr_new;
}

// ***************************************************************************
// xstructure AssignNameXstr(xstructure& xstr, vector<string> names)
// ***************************************************************************
xstructure AssignNameXstr(xstructure& xstr, vector<string> names) {
  //Assign names to xstr
  xstructure xstr_new = xstr;
  int iatom=0;
  for (uint itype=0; itype<xstr.num_each_type.size();itype++) {
    string species = string(names.at(xstr_new.atoms.at(iatom).type));
    xstr_new.species.at(itype)=species;
    for (int j=0;j<xstr.num_each_type.at(itype);j++) {
      xstr_new.atoms.at(iatom).name=species;
      xstr_new.atoms.at(iatom).CleanName();
      xstr_new.atoms.at(iatom).CleanSpin();
      xstr_new.atoms.at(iatom).name_is_given=TRUE;
      iatom++;
    }
  }
  return xstr_new;
}


// ***************************************************************************
// vector<xstructure> AssignNameXstr(vector<xstructure>& vxstr, vector<string> names)
// ***************************************************************************
vector<xstructure> AssignNameXstr(vector<xstructure>& vxstr, vector<string> names) {
  vector<xstructure> vxstr_new;
  xstructure xstr_tmp;
  for (uint i=0; i<vxstr.size(); i++) {
    xstr_tmp = AssignNameXstr(vxstr.at(i), names);
    vxstr_new.push_back(xstr_tmp);
  }
  return vxstr_new;
}


// ***************************************************************************
// pocc::HNFCELL
// ***************************************************************************
namespace pocc {
  void HNFCELL(istream& input) {
    xstructure xstr_ori(input,IOVASP_AUTO);
    vector<xstructure> vxstr = Partial2Supercell(xstr_ori);
  }
}

// ***************************************************************************
// vector<xstructure> Partial2Supercell(xstructure xstr)
// ***************************************************************************
vector<xstructure> Partial2Supercell(xstructure xstr_ori) {
  xstructure xstr = AssignAlphabeticNameXstr(xstr_ori);
  string str_species_ori;
  if(xstr_ori.atoms.at(0).name_is_given) {str_species_ori = xstr_ori.SpeciesString();}
  else {str_species_ori = "As Ba Cl Dy Eu Fe Ga He In K";}
  vector<string> vxstr_species_ori;
  aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");

  ofstream FileMESSAGE;
  FileMESSAGE.open("LOG.POCC");
  _aflags aflags;

  int nHNF = InitializeXstr(xstr, vxstr_species_ori, FileMESSAGE, aflags);
  vector<vector<xstructure> > vv_xstr;
  if(CheckFivefoldOccupied(xstr)) { //fivefold or quadruple
    vv_xstr = Partial2Xstr_Fivefold_Occupied(xstr, nHNF, FileMESSAGE, aflags);
  }
  else if(CheckFourfoldOccupied(xstr)) { //fourfold or quadruple
    vv_xstr = Partial2Xstr_Fourfold_Occupied(xstr, nHNF, FileMESSAGE, aflags);
  }
  else if(CheckTripleOccupied(xstr)) { //Triple
    vv_xstr = Partial2Xstr_Triple_Occupied(xstr, nHNF, FileMESSAGE, aflags);
  }
  else if(CheckDoubleOccupied(xstr)) {  //Double
    vv_xstr = Partial2Xstr(xstr, nHNF, FileMESSAGE, aflags);
  }

  vector<xstructure> vxstr;
  for (uint i=0; i<vv_xstr.size();i++) {
    for (uint j=0; j<vv_xstr.at(i).size();j++) {
      xstructure xstr_tmp = vv_xstr.at(i).at(j);
      xstructure xstr_tmp_name = AssignNameXstr(xstr_tmp, vxstr_species_ori);
      xstr_tmp_name.SpeciesPutAlphabetic();
      vxstr.push_back(xstr_tmp_name);
    }
  }
  vector<xstructure> vxstr_sorted = pocc::SortGroupXstrUFFEnergy(vxstr);

  ostringstream aus;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aus << "0000 MESSAGE    Printing Derivative Supercells ... " << Message(aflags, "user, host, time") << endl;
  for (uint i=0; i<vxstr_sorted.size();i++) {
    aus << AFLOWIN_SEPARATION_LINE << endl;
    aus << "[VASP_POSCAR_MODE_EXPLICIT]START.HNFCELL_" << i+1 << endl;
    aus << vxstr_sorted.at(i);
    aus << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNFCELL_" << i+1 << endl;
    aus << AFLOWIN_SEPARATION_LINE << endl;
    aus.flush();
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus.str("");
  }
  FileMESSAGE.close();
  return vxstr_sorted;
}

// ***************************************************************************
// vector<xstructure> CalculateInitialSupercell(xstructure xstr, int n, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<xstructure> CalculateInitialSupercell(xstructure xstr, int n, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET=FALSE;
  ostringstream aus;
  aus << "0000 MESSAGE    Calculating HNF index ... " << Message(aflags, "user, host, time") << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  vector<xmatrix<double> > sHNF=CalculateHNF(xstr, n);
  vector<xstructure> vxstr_init;
  xstructure b_tmp;
  aus << "0000 MESSAGE    Generating supercells ... " << Message(aflags, "user, host, time") << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  stringstream ss;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  for (int i=0; i<(int) sHNF.size();i++) {
    b_tmp=GetSuperCell(xstr,trasp(sHNF.at(i)));//when generates structure, using trasp
    b_tmp.title+="  [HNF("+aurostd::utype2string(n)+")="+aurostd::utype2string(i+1)+"/"+aurostd::utype2string(sHNF.size())+"=";
    string str_sHNF = " " 
      + aurostd::utype2string(sHNF.at(i)(1,1))+" " 
      + aurostd::utype2string(sHNF.at(i)(1,2))+" "
      + aurostd::utype2string(sHNF.at(i)(1,3))+"; "
      + aurostd::utype2string(sHNF.at(i)(2,1))+" "
      + aurostd::utype2string(sHNF.at(i)(2,2))+" "
      + aurostd::utype2string(sHNF.at(i)(2,3))+"; "
      + aurostd::utype2string(sHNF.at(i)(3,1))+" "
      + aurostd::utype2string(sHNF.at(i)(3,2))+" "
      + aurostd::utype2string(sHNF.at(i)(3,3))+"]";
    b_tmp.title+=str_sHNF;
    b_tmp.partial_occupation_flag=FALSE; //Not partial occupation anymore
    vxstr_init.push_back(b_tmp);
    aus << "[VASP_POSCAR_MODE_EXPLICIT]START.HNF_" << n << "_" << i+1 << "_" << sHNF.size() << endl;
    aus << CleanNameXstr(b_tmp);
    aus << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNF_" << n << "_" << i+1 << "_" << sHNF.size() << endl;
    aus << AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus.str("");
  }
  aus << "0000 MESSAGE    Total number of initial supercells is " << vxstr_init.size() <<" " << Message(aflags, "user, host, time") << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");

  return vxstr_init;
}

// ***************************************************************************
// int InitializeXstr(xstructure &xstr, vector<string> vxstr_species_ori, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
int InitializeXstr(xstructure &xstr, vector<string> vxstr_species_ori, ofstream &FileMESSAGE, _aflags &aflags) {
  //This part is based on "pflow::HNFTOL" 
  ostringstream oss;
  xstr.BringInCell();
  //xstr.CalculateSymmetryPointGroup(); //CO, cannot calculate symmetry of a structure with overlapping atoms
  xstr.ReScale(1.0);

  uint digits=6, digits1=4, digits2=2*digits+11;

  if(xstr.partial_occupation_HNF) {return xstr.partial_occupation_HNF;} //if HNF exists, no need to optimize pocc values
  double tolerance=DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
  if(xstr.partial_occupation_tol>0) {tolerance=xstr.partial_occupation_tol;}
  oss << "[AFLOW] hnf_tolerance=" << tolerance << endl;
  oss << "Printing Input PARTCAR ... " << Message(aflags, "user, host, time") << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
  oss << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl;
  oss << AssignNameXstr(xstr, vxstr_species_ori);// << endl;
  oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
  oss << "Optimizing HNF number ... " << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
  oss.precision(digits);
  double error=1.0,eps;
  int HNFi=0;
  vector<double> error_iatom;  vector<double> effective_pocc_iatom;  vector<uint> M_iatom;
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {error_iatom.push_back(1.0);}
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {effective_pocc_iatom.push_back(0.0);}
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {M_iatom.push_back(0);}

  oss << aurostd::PaddedPRE(aurostd::utype2string(0),digits1) << "  " << "--------" << " ";
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
    if(xstr.atoms.at(iatom).partial_occupation_value<1.0) {
      oss << "| "+aurostd::PaddedPOST("iatom="+aurostd::utype2string(iatom)+"/"+aurostd::utype2string(xstr.atoms.size()),digits2) << " " ;
    }
  }
  oss << " | " << "error" << endl;
  while(error>tolerance) {
    HNFi++;
    error=1.0/HNFi;
    oss << aurostd::PaddedPRE(aurostd::utype2string(HNFi),digits1) << "  " << error << " ";//endl;
    for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
      error_iatom.at(iatom)=1.0;
      M_iatom.at(iatom)=0;
      for(int j=0;j<=HNFi;j++) {
	eps=abs((double) j/HNFi-xstr.atoms.at(iatom).partial_occupation_value);
	if(eps<error_iatom.at(iatom)) {
	  error_iatom.at(iatom)=eps;
	  M_iatom.at(iatom)=(uint) j;
	  effective_pocc_iatom.at(iatom)=(double) M_iatom.at(iatom)/HNFi;
	}
      }
    }

    for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
      if(xstr.atoms.at(iatom).partial_occupation_value<1.0) {
	stringstream aus;aus.precision(digits);aus.setf(std::ios::fixed,std::ios::floatfield);
	aus << effective_pocc_iatom.at(iatom) << "," << error_iatom.at(iatom) << "," << M_iatom.at(iatom) << "/" << HNFi;
	oss << "| " << aurostd::PaddedPOST(aus.str(),digits2) << " " ;
      }
    }
    error=aurostd::max(error_iatom);
    oss << " | " << error << endl;
  }
  oss << AFLOWIN_SEPARATION_LINE << endl;

  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
    if(xstr.atoms.at(iatom).partial_occupation_value<1.0) {xstr.atoms.at(iatom).partial_occupation_value = effective_pocc_iatom.at(iatom) ;}
  }
  UpdateXstr_comp_each_type(xstr); //update partial occupation
  oss << "Printing optimized input PARTCAR ..." << Message(aflags, "user, host, time") << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;
  oss << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl;
  oss << AssignNameXstr(xstr, vxstr_species_ori);// << endl;
  oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;

  aurostd::PrintMessageStream(FileMESSAGE, oss,XHOST.QUIET);
  int nHNF=HNFi;
  return nHNF;
}

// ***************************************************************************
//vector<vector<xstructure> > Partial2Xstr_Fivefold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr_Fivefold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET= TRUE;
  //convert fivefold into fourfold firstly
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
  xstructure xstr_norm = NormalizeXstructure(xstr);
  xstructure xstr_tri = xstr;
  for (uint i=0;i<xstr_norm.atoms.size();i++) {
    if(NormalisedNumber.at(i).size()==5) {
      int k=NormalisedNumber.at(i).at(4);  //begin from 0, so the 5th atom should be at(4)
      xstr_tri.RemoveAtom(k);
    }
  }
  //Call fourfold routine
  vector<vector<xstructure> > vec_vec_xstr = Partial2Xstr_Fourfold_Occupied(xstr_tri, nHNF, FileMESSAGE, aflags);
  ostringstream aus;

  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, nHNF);
  vector<string>  element_name; //replacing original name
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if((NormalisedNumber.at(i).size()==5) || 
       ((NormalisedNumber.at(i).size()==4)&& CheckVacancyOnOnesite(xstr,i)))
      { //tripled occupied
	int k1= NormalisedNumber.at(i).at(0);
	int k2= NormalisedNumber.at(i).at(1);
	int k3= NormalisedNumber.at(i).at(2);
	int k4= NormalisedNumber.at(i).at(3);
	double sum_first_four_pocc_value = xstr.atoms.at(k1).partial_occupation_value
	  + xstr.atoms.at(k2).partial_occupation_value 
	  + xstr.atoms.at(k3).partial_occupation_value
	  + xstr.atoms.at(k4).partial_occupation_value;
	int num_vacancy= hnf_double2int((1-sum_first_four_pocc_value)*nHNF);
	for (int k=0; k<num_vacancy; k++) {
	  if(NormalisedNumber.at(i).size()==5) { //If it contains four elements
	    int k5=NormalisedNumber.at(i).at(4); //get the the name of the 5th element
	    element_name.push_back(xstr.atoms.at(k5).name);
	  }
	  else { //if it contains three elements and one vacancy
	    element_name.push_back(name_vacancy);
	  }
	}
	vector<int> vec_num_part=atom_num_supercell.at(i);
	//get rid of the numbers of the first three atom
	//the 4th and 5th atoms will be thought as one entry, so we need to remove the first three
	float first_three_pocc_value = xstr.atoms.at(k1).partial_occupation_value 
	  + xstr.atoms.at(k2).partial_occupation_value
	  + xstr.atoms.at(k3).partial_occupation_value;
	int num_first_three = hnf_double2int(first_three_pocc_value*nHNF);
	vector<int> vec_num_part_new; 
	for (uint i=num_first_three; i<vec_num_part.size();i++) {
	  vec_num_part_new.push_back(vec_num_part.at(i));
	}
	vector<vector<int> > num_com_vacancy; num_com_vacancy.clear();
	combine(vec_num_part_new, num_com_vacancy, num_vacancy);
	vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
      }
  }
  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  vector<vector<int> > formated_final_num;

  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  vector<int> tmp;
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
	tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }


  vector<vector<xstructure> > vec_vec_xstr_fivefold;
  int total=0;
  for (uint k=0; k<vec_vec_xstr.size();k++) {  //2-D vector
    vector<xstructure> vec_xstr=vec_vec_xstr.at(k);
    vector<xstructure> vec_new_xstr;
    for (uint i=0; i<vec_xstr.size();i++) {
      xstructure xstr_tmp = vec_xstr.at(i);
      for (uint j=0; j<formated_final_num.size();j++) {
	total++;
	vector<int> vec_int = formated_final_num.at(j);
	xstructure new_xstr = XstrSubstitute(xstr_tmp, vec_int, element_name);
	vec_new_xstr.push_back(new_xstr);
      }
    }
    vec_vec_xstr_fivefold.push_back(vec_new_xstr);
  }

  vector<vector<xstructure> > vec_vec_xstr_final;
  for (uint i=0; i<vec_vec_xstr_fivefold.size();i++) {
    vector<xstructure> vec_new_xstr = vec_vec_xstr_fivefold.at(i);
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(vec_new_xstr, FileMESSAGE, aflags); 
    vec_vec_xstr_final.push_back(vxstr_final);
  }

  return vec_vec_xstr_final;
}

// ***************************************************************************
// vector<vector<xstructure> > Partial2Xstr_Fourfold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr_Fourfold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET= TRUE;
  //convert fourfold into triple
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
  xstructure xstr_norm = NormalizeXstructure(xstr);
  xstructure xstr_tri = xstr;
  for (uint i=0;i<xstr_norm.atoms.size();i++) {
    if(NormalisedNumber.at(i).size()==4) { //if the 4th element is vacancy, then it is fake triple one
      int k=NormalisedNumber.at(i).at(3);
      xstr_tri.RemoveAtom(k);
    }
  }

  vector<vector<xstructure> > vec_vec_xstr = Partial2Xstr_Triple_Occupied(xstr_tri, nHNF, FileMESSAGE, aflags);
  ostringstream aus;

  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, nHNF);
  vector<string>  element_name; //replacing original name
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if((NormalisedNumber.at(i).size()==4) || 
       ((NormalisedNumber.at(i).size()==3)&& CheckVacancyOnOnesite(xstr,i)))
      { //tripled occupied
	int k1= NormalisedNumber.at(i).at(0);
	int k2= NormalisedNumber.at(i).at(1);
	int k3= NormalisedNumber.at(i).at(2);
	double sum_first_three_pocc_value = xstr.atoms.at(k1).partial_occupation_value
	  + xstr.atoms.at(k2).partial_occupation_value + xstr.atoms.at(k3).partial_occupation_value;
	int num_vacancy= hnf_double2int((1-sum_first_three_pocc_value)*nHNF);
	for (int k=0; k<num_vacancy; k++) {
	  if(NormalisedNumber.at(i).size()==4) { //If it contains four elements
	    int k4=NormalisedNumber.at(i).at(3); //get the the name of  4th element
	    element_name.push_back(xstr.atoms.at(k4).name);
	  }
	  else { //if it contains three elements and one vacancy
	    element_name.push_back(name_vacancy);
	  }
	}
	vector<int> vec_num_part=atom_num_supercell.at(i);
	//get rid of the numbers of the first two atom
	float first_two_pocc_value = xstr.atoms.at(k1).partial_occupation_value + xstr.atoms.at(k2).partial_occupation_value;
	int num_first_two = hnf_double2int(first_two_pocc_value*nHNF);
	vector<int> vec_num_part_new; 
	for (uint i=num_first_two; i<vec_num_part.size();i++) {
	  vec_num_part_new.push_back(vec_num_part.at(i));
	}
	vector<vector<int> > num_com_vacancy; num_com_vacancy.clear();
	combine(vec_num_part_new, num_com_vacancy, num_vacancy);
	vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
      }
  }
  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  vector<vector<int> > formated_final_num;

  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  vector<int> tmp;
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
	tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }


  vector<vector<xstructure> > vec_vec_xstr_fourfold;
  int total=0;
  for (uint k=0; k<vec_vec_xstr.size();k++) {  //2-D vector
    vector<xstructure> vec_xstr=vec_vec_xstr.at(k);
    vector<xstructure> vec_new_xstr;
    for (uint i=0; i<vec_xstr.size();i++) {
      xstructure xstr_tmp = vec_xstr.at(i);
      for (uint j=0; j<formated_final_num.size();j++) {
	total++;
	vector<int> vec_int = formated_final_num.at(j);
	xstructure new_xstr = XstrSubstitute(xstr_tmp, vec_int, element_name);
	vec_new_xstr.push_back(new_xstr);
      }
    }
    vec_vec_xstr_fourfold.push_back(vec_new_xstr);
  }

  vector<vector<xstructure> > vec_vec_xstr_final;
  for (uint i=0; i<vec_vec_xstr_fourfold.size();i++) {
    vector<xstructure> vec_new_xstr = vec_vec_xstr_fourfold.at(i);
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(vec_new_xstr, FileMESSAGE, aflags); 
    vec_vec_xstr_final.push_back(vxstr_final);
  }
  return vec_vec_xstr_final;
}

// ***************************************************************************
// vector<vector<xstructure> > Partial2Xstr_Triple_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr_Triple_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET=TRUE;
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
  xstructure xstr_norm = NormalizeXstructure(xstr);
  //convert triple into double, in fact, it is not very necessary here
  xstructure xstr_tri = xstr;
  for (uint i=0;i<xstr_norm.atoms.size();i++) {
    if(NormalisedNumber.at(i).size()==3) { //if the 3th element is vacancy, then it is fake triple one
      int k=NormalisedNumber.at(i).at(2);
      xstr_tri.RemoveAtom(k);
    }
  }

  vector<vector<xstructure> > vec_vec_xstr = Partial2Xstr(xstr_tri, nHNF, FileMESSAGE, aflags);
  ostringstream aus;

  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, nHNF);
  vector<string>  element_name; //replacing original name
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if((NormalisedNumber.at(i).size()==3) || 
       ((NormalisedNumber.at(i).size()==2)&& CheckVacancyOnOnesite(xstr,i)))
      { //tripled occupied
	int k1= NormalisedNumber.at(i).at(0);
	int k2= NormalisedNumber.at(i).at(1);
	double sum_first_two_pocc_value = xstr.atoms.at(k1).partial_occupation_value
	  + xstr.atoms.at(k2).partial_occupation_value;
	int num_vacancy= hnf_double2int((1-sum_first_two_pocc_value)*nHNF);
	for (int k=0; k<num_vacancy; k++) {
	  if(NormalisedNumber.at(i).size()==3) { //If it contains three elements
	    int k3=NormalisedNumber.at(i).at(2);
	    element_name.push_back(xstr.atoms.at(k3).name);
	  }
	  else { //if it contains two elements and one vacancy
	    element_name.push_back(name_vacancy);
	  }
	}
	vector<int> vec_num_part=atom_num_supercell.at(i);
	//get rid of the numbers of the first atom
	int num_first_one = hnf_double2int(xstr.atoms.at(k1).partial_occupation_value*nHNF);
	vector<int> vec_num_part_new; 
	for (uint i=num_first_one; i<vec_num_part.size();i++) {
	  vec_num_part_new.push_back(vec_num_part.at(i));
	}
	vector<vector<int> > num_com_vacancy; num_com_vacancy.clear();
	combine(vec_num_part_new, num_com_vacancy, num_vacancy);
	vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
      }
  }
  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  vector<vector<int> > formated_final_num;

  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  vector<int> tmp;
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
	tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }


  vector<vector<xstructure> > vec_vec_xstr_triple;
  int total=0;
  for (uint k=0; k<vec_vec_xstr.size();k++) {  //2-D vector
    vector<xstructure> vec_xstr=vec_vec_xstr.at(k);
    vector<xstructure> vec_new_xstr;
    for (uint i=0; i<vec_xstr.size();i++) {
      xstructure xstr_tmp = vec_xstr.at(i);
      for (uint j=0; j<formated_final_num.size();j++) {
	total++;
	vector<int> vec_int = formated_final_num.at(j);
	xstructure new_xstr = XstrSubstitute(xstr_tmp, vec_int, element_name);
	vec_new_xstr.push_back(new_xstr);
      }
    }
    vec_vec_xstr_triple.push_back(vec_new_xstr);
  }

  vector<vector<xstructure> > vec_vec_xstr_final;
  for (uint i=0; i<vec_vec_xstr_triple.size();i++) {
    vector<xstructure> vec_new_xstr = vec_vec_xstr_triple.at(i);
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(vec_new_xstr, FileMESSAGE, aflags); 
    vec_vec_xstr_final.push_back(vxstr_final);
  }
  return vec_vec_xstr_final;
}


// ***************************************************************************
// vector<vector<xstructure> > Partial2Xstr(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  //doubly occupied
  xstr.ReScale(1.0);
  xstr.iomode=IOVASP_POSCAR;

  XHOST.QUIET= FALSE;
  ostringstream aus;
  aus << "0000 MESSAGE    POCC running " << Message(aflags, "user, host, time") << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  xstructure xstr_c = NormalizeXstructure(xstr);  

  int n = nHNF;
  vector<xstructure> vxstr_init =  CalculateInitialSupercell(xstr_c, n, FileMESSAGE, aflags);
  aus << "0000 MESSAGE    Generating derivate supercells ..." << Message(aflags, "user, host, time") << endl;
  aus << "0000 MESSAGE    Please be patient ..." << Message(aflags, "user, host, time") << endl;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);

  vector<xstructure> vxstr_mid; vxstr_mid.clear(); 
  vector<vector<xstructure> > vv_xstr; vv_xstr.clear();
  for (uint i=0; i<vxstr_init.size();i++) {
    xstructure xstr_tmp =vxstr_init.at(i);
    vxstr_mid=AssignPartialOccupation(xstr_tmp, xstr, n);
    vv_xstr.push_back(vxstr_mid);
  }

  vector<xstructure> vxstr_mid2; vxstr_mid2.clear(); 
  vector<vector<xstructure> > vv_xstr2; vv_xstr2.clear();
  int total=0;
  for (uint i=0; i<vv_xstr.size();i++) {
    for (uint j=0; j<vv_xstr.at(i).size();j++) {
      total++;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      aus << "[VASP_POSCAR_MODE_EXPLICIT]START.HNF_" << i+1 << "_" << j+1 << "_" << total << endl;
      //string str_title = " HNF_" + aurostd::double2string(i+1) + "_" + aurostd::double2string(j+1) + "_" + aurostd::double2string(total);
      //vv_xstr.at(i).at(j).title+=str_title;
      aus << vv_xstr.at(i).at(j);
      vxstr_mid2.push_back(vv_xstr.at(i).at(j));
      aus << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNF_" << i+1 << "_" << j+1 << "_" << total << endl;
      //aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
      aus.str("");
    }
    vv_xstr2.push_back(vxstr_mid2);
    vxstr_mid2.clear();
  }

  vector<vector<xstructure> > vv_xstr_final; vv_xstr_final.clear();
  for (uint i=0;i<vv_xstr2.size();i++) {
    vector<xstructure> vxstr_tmp = vv_xstr2.at(i);
    vector<xstructure> vxstr_nq2 = RemoveEquivalentXstr(vxstr_tmp, FileMESSAGE, aflags); //
    vv_xstr_final.push_back(vxstr_nq2);
  } 
  return vv_xstr_final;
}

// ***************************************************************************
// vector<xstructure> CalculatePrimitiveCell(vector<xstructure> &vec_xstr)
// ***************************************************************************
vector<xstructure> CalculatePrimitiveCell(vector<xstructure> &vec_xstr) {
  cerr << "Calculating Pritimive cells" << endl;
  vector<xstructure> vec_xstr_sp;
  vec_xstr_sp.clear();
  xstructure xstr_tmp_sp, a;
  for (uint i=0; i<vec_xstr.size();i++) {
    a = vec_xstr.at(i);
    xstr_tmp_sp = GetStandardPrimitive(a);
    vec_xstr_sp.push_back(xstr_tmp_sp);
    cerr << i+1 << " is finished!" << endl;
  }
  return vec_xstr_sp;
}
// ***************************************************************************
// bool comparison_atom_fpos(_atom atom1, _atom atom2)
// ***************************************************************************
bool comparison_atom_fpos(_atom atom1, _atom atom2) {
  double _EQUAL_DOUBLE = 1E-6;
  if( abs(atom1.fpos[1] - atom2.fpos[1]) > _EQUAL_DOUBLE ) {
    return (atom1.fpos[1] < atom2.fpos[1]);
  } else if( abs(atom1.fpos[2] - atom2.fpos[2]) > _EQUAL_DOUBLE ) {
    return (atom1.fpos[2] < atom2.fpos[2]);
  } else if( abs(atom1.fpos[3] - atom2.fpos[3]) > _EQUAL_DOUBLE ) {
    return (atom1.fpos[3] < atom2.fpos[3]);
  } else {
    return true;
  }
}

// ***************************************************************************
// pflow::Sort_atom_fpos(xstructure &xstr)
// ***************************************************************************
namespace pflow {
  void Sort_atom_fpos(xstructure &xstr) {
    xstr.BringInCell();
    xstr.SetCoordinates(_COORDS_FRACTIONAL_);
    int begin = 0;
    for (uint i=0; i <xstr.num_each_type.size(); i++) {
      int j= xstr.num_each_type.at(i);
      sort(xstr.atoms.begin()+begin, xstr.atoms.begin()+begin+j, comparison_atom_fpos);
      begin += j;
    }
  }
}

// ***************************************************************************
// bool comparison_atom_cpos(_atom atom1, _atom atom2)
// ***************************************************************************
bool comparison_atom_cpos(_atom atom1, _atom atom2) {
  double _EQUAL_DOUBLE = 1E-6;
  if( abs(atom1.cpos[1] - atom2.cpos[1]) > _EQUAL_DOUBLE ) {
    return (atom1.cpos[1] < atom2.cpos[1]);
  } else if( abs(atom1.cpos[2] - atom2.cpos[2]) > _EQUAL_DOUBLE ) {
    return (atom1.cpos[2] < atom2.cpos[2]);
  } else if( abs(atom1.cpos[3] - atom2.cpos[3]) > _EQUAL_DOUBLE ) {
    return (atom1.cpos[3] < atom2.cpos[3]);
  } else {
    return true;
  }
}

// ***************************************************************************
// pflow::Sort_atom_cpos(xstructure &xstr)
// ***************************************************************************
namespace pflow {
  void Sort_atom_cpos(xstructure &xstr) {
    xstr.BringInCell();
    xstr.SetCoordinates(_COORDS_CARTESIAN_);
    int begin = 0;
    for (uint i=0; i <xstr.num_each_type.size(); i++) {
      int j= xstr.num_each_type.at(i);
      sort(xstr.atoms.begin()+begin, xstr.atoms.begin()+begin+j, comparison_atom_cpos);
      begin += j;
    }
  }
}

// ***************************************************************************
// bool LatticeCompare(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
bool LatticeCompare(xstructure xstr1, xstructure xstr2) {
  bool RUN_FLAG = false;
  double tmp1=0, tmp2=0, sum=0, epsilon=1E-4;
  tmp1 = abs(xstr1.a-xstr2.a) + abs(xstr1.b-xstr2.b) + abs(xstr1.c-xstr2.c);
  tmp2 = abs(xstr1.alpha-xstr2.alpha) + abs(xstr1.beta-xstr2.beta) + abs(xstr1.gamma-xstr2.gamma);
  sum = tmp1 +tmp2;  //Make sure sum is absolute value
  if(sum <epsilon) {RUN_FLAG = true;}
  return RUN_FLAG;
}

// ***************************************************************************
// bool CoordCompare(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
bool CoordCompare(xstructure xstr1, xstructure xstr2) {
  bool RUN_FLAG = false;
  double sum=0, epsilon=1E-4;
  if(!LatticeCompare(xstr1, xstr2)) {return RUN_FLAG;}
  else  {
    if(xstr1.atoms.size()==xstr2.atoms.size()) {
      for (uint i=0;i<xstr1.atoms.size();i++) {
	double b1=0, b2=0, b3=0;
	b1= xstr1.atoms.at(i).fpos[1] - xstr2.atoms.at(i).fpos[1];
	b2= xstr1.atoms.at(i).fpos[2] - xstr2.atoms.at(i).fpos[2];
	b3= xstr1.atoms.at(i).fpos[3] - xstr2.atoms.at(i).fpos[3];
	sum += (b1*b1+b2*b2+b3*b3);
      }
      if(sum <epsilon) {RUN_FLAG = true;}
    }
  }
  return RUN_FLAG;
}

// ***************************************************************************
// vector<xstr_atom_number> CalculateXstrNumEachType(xstructure xstr)
// ***************************************************************************
vector<xstr_atom_number> CalculateXstrNumEachType(xstructure xstr) {
  vector<xstr_atom_number> vec_xstr_number;
  int begin=0;
  for (uint i=0; i<xstr.num_each_type.size();i++) {
    int k=xstr.num_each_type.at(i);
    xstr_atom_number atom_number_tmp;
    atom_number_tmp.number = k;
    for (int j=0; j<k; j++) {
      atom_number_tmp.vec_atom_number.push_back(begin);
      begin++;
    }
    vec_xstr_number.push_back(atom_number_tmp);
  }
  return vec_xstr_number;
}

bool sort_atom_number(xstr_atom_number str1, xstr_atom_number str2) {
  return(str1.number<str2.number);
}

// ***************************************************************************
// vector<strno_energy> xstructure2strno_energy(vector<xstructure> vec_xstr)
// ***************************************************************************
vector<xstr_energy> xstructure2xstr_energy(vector<xstructure> vec_xstr) {
  vector<xstr_energy> vxstr_energy;
  xstr_energy xstr_energy_tmp;
  for (uint i=0; i<vec_xstr.size();i++) {
    //xstructure xstr_tmp = vec_xstr.at(i);
    //xstr_energy_tmp.xstr = xstr_tmp;
    xstr_energy_tmp.xstr = vec_xstr.at(i);
    double TotalEnergyTmp;
    if(CALCULATE_ENERGY_UFF) {  
      //TotalEnergyTmp = pocc::CalculateUFFEnergy(xstr_tmp);
      TotalEnergyTmp = pocc::CalculateUFFEnergy(vec_xstr.at(i));
    }
    else {
      //TotalEnergyTmp = pocc::CalculateEnergyUsingGulp(xstr_tmp);
      TotalEnergyTmp = pocc::CalculateEnergyUsingGulp(vec_xstr.at(i));
    }
    xstr_energy_tmp.energy = TotalEnergyTmp;
    vxstr_energy.push_back(xstr_energy_tmp);
  }
  return vxstr_energy;
}

// ***************************************************************************
// vector<xstructure> RemoveEquivalentXstr(vector<xstructure> vec_xstr, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<xstructure> RemoveEquivalentXstr(vector<xstructure> vec_xstr, ofstream &FileMESSAGE, _aflags &aflags) {
  vector<xstr_energy> vec_xstr_energy = xstructure2xstr_energy(vec_xstr);
  const double epsilon_etot = 1E-6;
  vector<xstr_energy> vec_xstr_energy_final;
  vec_xstr_energy_final.push_back(vec_xstr_energy.at(0));
  for (uint i=1; i<vec_xstr_energy.size();i++) {
    uint j;
    for (j=0; j<vec_xstr_energy_final.size();j++) {
      if(abs(vec_xstr_energy.at(i).energy-vec_xstr_energy_final.at(j).energy)<epsilon_etot) { //if they have same energy
	break;
      }
    }
    if(j==vec_xstr_energy_final.size()) {
      vec_xstr_energy_final.push_back(vec_xstr_energy.at(i));
    }
  }

  //cerr << "KESONG TEST" << endl;
  vector<int> vecDG; vecDG.clear(); //Degenracy
  for (uint i=0;i<vec_xstr_energy_final.size();i++) {
    int DGi =0;
    for (uint j=0; j<vec_xstr_energy.size();j++) {
      if(abs(vec_xstr_energy_final.at(i).energy- vec_xstr_energy.at(j).energy)<epsilon_etot) {DGi++;}
    }
    vecDG.push_back(DGi);
  }

  int number = vec_xstr_energy.size() - vec_xstr_energy_final.size();
  ostringstream aus;
  aus << "0000 MESSAGE    " << number << " equivalent structures are removed!" << Message(aflags, "user, host, time") << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  vector<xstructure> vec_xstr_final;
  for (uint i=0;i<vec_xstr_energy_final.size();i++) {
    xstructure xstr_tmp = vec_xstr_energy_final.at(i).xstr;
    string str_title = " DG=" + aurostd::utype2string(vecDG.at(i));
    xstr_tmp.title += str_title;
    vec_xstr_final.push_back(xstr_tmp);  //xstr_energy, structure, which contains xstr, and energy
  }
  return vec_xstr_final;
}

// ***************************************************************************
// pocc::DIFF(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF_UFFENERGY(xstructure xstr1, xstructure xstr2) {
    const double  energy_eps_tol_ =  1E-10;
    double e1 = pocc::CalculateUFFEnergy(xstr1);
    double e2 = pocc::CalculateUFFEnergy(xstr2);
    double e_delta = abs(e1-e2);
    if(e_delta < energy_eps_tol_) {return true;}
    else {return false;}
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF(xstructure xstr1, xstructure xstr2) {
    bool RUN_FLAG = false;

    if(pocc::DIFF_LEVEL1(xstr1, xstr2)) { RUN_FLAG = true; return RUN_FLAG;}
    //Highest lelvel
    if(pocc::DIFF_LEVEL4(xstr1, xstr2)) { RUN_FLAG = true; return RUN_FLAG;}
    if(pocc::DIFF_LEVEL4(xstr2, xstr1)) { RUN_FLAG = true; return RUN_FLAG;}
    return RUN_FLAG;
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF_LEVEL1(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF_LEVEL1(xstructure xstr1, xstructure xstr2) {
    //pocc::DIFF_Level: 
    //Only compare coordinates, assume names_is_given
    bool RUN_FLAG = false;
    xstr1.SpeciesPutAlphabetic();
    xstr2.SpeciesPutAlphabetic();
    xstr1.ShifOriginToAtom(0);
    pflow::Sort_atom_cpos(xstr1);
    for (int i=0; i<xstr2.num_each_type.at(0);i++) {
      xstr2.ShifOriginToAtom(i);
      xstr2.BringInCell();
      pflow::Sort_atom_cpos(xstr2);
      if(CoordCompare(xstr1, xstr2)) {
	RUN_FLAG = true;
	break;
      }
    }
    return RUN_FLAG;
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF_LEVEL4(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF_LEVEL4(xstructure xstr1, xstructure xstr2) {
    //pocc::DIFF_Level: 
    //This is the highest level to compare two POSCARs
    bool RUN_FLAG = false;
    xstructure xstr1_sp, xstr2_sp;
    xstr1_sp = GetStandardPrimitive(xstr1);
    xstr2_sp = GetStandardPrimitive(xstr2);
    xstr1_sp.ShifOriginToAtom(0);
    xstr1_sp.BringInCell();
    pflow::Sort_atom_cpos(xstr1_sp);

    vector<xstructure> vec_rot = pocc::GenerateRotatedXstr(xstr2_sp);
    for (uint i=0; i<vec_rot.size();i++) {
      xstructure xstr2_sp_rot_tmp = vec_rot.at(i);
      for (uint j=0; j<xstr2_sp_rot_tmp.atoms.size();j++) {
	xstr2_sp_rot_tmp.ShifOriginToAtom(j);
	xstr2_sp.BringInCell();
	pflow::Sort_atom_cpos(xstr2_sp_rot_tmp);
	if(CoordCompare(xstr1_sp, xstr2_sp_rot_tmp)) {
	  RUN_FLAG = true;
	  break;
	}
      }
    }
    return RUN_FLAG;
  }
} // namespace pocc

// ***************************************************************************
// vector<xstructure> pocc::GenerateRotatedXstr(xstructure xstr)
// ***************************************************************************
namespace pocc {
  vector<xstructure> GenerateRotatedXstr(xstructure xstr) {
    //only including rotations, because origin shift is alreay included in pocc::DIFF
    vector<xstructure> vec_xstr_rotated;
    xstructure xstr_tmp=xstr;
    
    xstr.CalculateSymmetryPointGroup();
    
    for(uint i=0; i<xstr.pgroup.size(); i++) {
      for (uint j=0; j<xstr.atoms.size();j++) {
	xstr_tmp.atoms.at(j).fpos = xstr.pgroup.at(i).Uf*xstr.atoms.at(j).fpos;
      }
      vec_xstr_rotated.push_back(xstr_tmp); 
    }
    return vec_xstr_rotated;
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF(string options)
// ***************************************************************************
namespace pocc {
  void DIFF(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) {cerr << "pocc::DIFF: BEGIN" << endl;}
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    
    if(LDEBUG) {cerr << "pocc::DIFF: options=[" << options << "]" << endl;}
    if(LDEBUG) {cerr << "pocc::DIFF: tokens.size()=" << tokens.size() << endl;}
    if(LDEBUG) {for(uint i=0;i<tokens.size();i++) cerr << "pocc::DIFF: tokens.at(i)=" << tokens.at(i) << endl;}
    
    if(tokens.size()<2) {
      init::ErrorOption(cout,options,"pocc::DIFF","aflow --diff=POSCAR1,POSCAR2");
      exit(0);
    }
    if(!aurostd::FileExist(tokens.at(0))) {cerr << tokens.at(0) << " does not exist!" <<endl; exit(1);}
    if(!aurostd::FileExist(tokens.at(1))) {cerr << tokens.at(1) << " does not exist!" <<endl; exit(1);}
    xstructure xstr1, xstr2;
    xstr1 = xstructure(tokens.at(0), IOAFLOW_AUTO);
    xstr2 = xstructure(tokens.at(1), IOAFLOW_AUTO);
    if(tokens.at(0)==tokens.at(1)) {
      cerr << "Do not try to fool me" << endl;
      cerr << "But I will still continue for test!" << endl;
    }
    cerr << "Now running void pocc::DIFF(string options) " << endl;
    
    if(!pocc::DIFF_UFFENERGY(xstr1, xstr2)) {
      cerr << "Not equivalent!" << endl;
    }
    else if(pocc::DIFF(xstr1, xstr2)) {
      cerr << "Equivalent!" << endl;
    }
    else {
      cerr << "Not equivalent!" << endl;
    }
  }
} // namespace pocc

// **************************************************************************
// xstructure::AddAtom_POCC
// **************************************************************************
// This adds an atom to the structure.
// When this code was written, AddAtom added atoms to the end of the structure, 
// majority of the code intact, this variant of AddAtom was added. 

void xstructure::AddAtom_POCC(const _atom& atom) {
  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  _atom btom;btom=atom;

  // check that this atom is not already present
  xvector<double> a1(3),a2(3),a3(3),aijk(3);
  a1=lattice(1);a2=lattice(2);a3=lattice(3);

  bool FOUND_POSITION=FALSE;
  for(uint iat=0;iat<atoms.size()&&FOUND_POSITION==FALSE;iat++) {
    if(atoms.at(iat).type==atom.type && atoms.at(iat).name==atom.name) {
      for(int i=-1;i<=1&&FOUND_POSITION==FALSE;i++) {
	for(int j=-1;j<=1&&FOUND_POSITION==FALSE;j++) {
	  for(int k=-1;k<=1&&FOUND_POSITION==FALSE;k++) {
	    aijk[1]=i;aijk[2]=j;aijk[3]=k;
	    //	if(aurostd::modulus(atoms.at(iat).cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atom.cpos))<0.1) FOUND_POSITION=TRUE;
	    if(aurostd::modulus(atoms.at(iat).fpos-(aijk+atom.fpos))<0.01) {FOUND_POSITION=TRUE;}
	  }
    }
      }
    }
  }
  if(FOUND_POSITION==TRUE) {return;} // found no need to add it further

  // now found that it does not exist check type
  //  cerr << "AddAtom new atom" << endl;
  bool FOUND_SPECIES=FALSE;
  uint species_position=0;
  for(uint isp=0;isp<species.size()&&FOUND_SPECIES==FALSE;isp++) {
    if(atom.name==species.at(isp)) {FOUND_SPECIES=TRUE;species_position=isp;}
  }
  
  if(FOUND_SPECIES==FALSE) {
    if(LDEBUG) {cerr << "AddAtom new_specie=" << atom.name << endl;}
    num_each_type.push_back(1);
    comp_each_type.push_back(atom.partial_occupation_value);
    species.push_back(atom.name); // cerr << "AddAtom=" << atom.name << endl;
    species_pp.push_back(atom.name); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_type.push_back(""); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_version.push_back(""); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_ZVAL.push_back(0.0); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_vLDAU.push_back(deque<double>()); // cerr << "AddAtom=" << atom.name << endl;
    species_volume.push_back(GetAtomVolume(atom.name)); // cerr << "AddAtom=" << atom.name << endl;
    species_mass.push_back(GetAtomMass(atom.name)); // cerr << "AddAtom=" << atom.name << endl;
  } else {
    // cerr <<  num_each_type.size() << " " <<  btom.type << endl;
    // cerr <<  comp_each_type.size() << " " <<  btom.type << endl;
    if(LDEBUG) {cerr << "AddAtom increasing species_position " << species_position << endl;}
    num_each_type.at(species_position)++;
    comp_each_type.at(species_position)+=atom.partial_occupation_value;
  }
  if(btom.name_is_given) {
    btom.CleanName();
    btom.CleanSpin();
  }

  // OLD STYLE
  atoms.push_back(btom);  MakeBasis(); return;

  /*// NEW STYLE
  bool found=FALSE;
  if(0)  for(uint iat=0;iat<atoms.size()&&!found;iat++) {
    if(iat<atoms.size()-1) {
      if(atoms.at(iat).type==btom.type && atoms.at(iat+1).type!=btom.type) {
	//	if(LDEBUG)
	cerr << "HERE1 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
	atoms.insert(iat+atoms.begin()+1,btom); // potential problem  with CAST
	found=TRUE;
      }
    }
  }
  if(1) {
    std::deque<_atom>::iterator it=atoms.begin();
    for(uint iat=0;iat<atoms.size()&&!found;iat++,it++) {
      if(iat<atoms.size()-1) {
	//	cerr << "HERE0 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;
	if(atoms.at(iat).type==btom.type && atoms.at(iat+1).type!=btom.type) {
	  //	if(LDEBUG)
	  //	  cerr << "HERE1 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
	  atoms.insert(it+1,btom);  // it is iterator, fine for insert.
	  found=TRUE;
	}
      }
    }
  }
  // if never found add at the end
  if(!found) atoms.push_back(btom);
  // 
  //  atoms.push_back(btom);
  // done
  MakeBasis(); // need to update NUMBER and BASIS */
}

#endif
// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG - Duke University 2010-2013              *
// *                                                                         *
// ***************************************************************************

