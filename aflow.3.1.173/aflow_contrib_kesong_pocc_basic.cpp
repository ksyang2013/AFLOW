// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2013               *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_kesong.cpp
// functions written by KESONG YANG
// 2010-2013: kesong.yang@gmail.com

#ifndef _AFLOW_CONTRIB_KESONG_MULTIENUM_CPP_
#define _AFLOW_CONTRIB_KESONG_MULTIENUM_CPP_

#include "aflow_contrib_kesong.h"



const double Tol_DecimalPart=1E-6;
const double Tol_FracDeciDiff=1E-6;
//This is the optimized tolerance
const int MaxDenominator=999999;

// [OBSOLETE] [WARNING_UNUSED] const int nCellMin=1;
// [OBSOLETE] [WARNING_UNUSED] const int nCellMax=3;

// ***************************************************************************
// void remove_duplicates(vector<T>& vec)
// ***************************************************************************
template <typename T>
void remove_duplicates(vector<T>& vec) {
  sort(vec.begin(), vec.end());
  vec.erase(unique(vec.begin(), vec.end()), vec.end());
}

str_num_data OptimizePoccValue(double dvalue, double tol) {
  str_num_data optimized_dvalue;

  if(tol<1E-6) {
      tol = 0.001;
  }

  int Nmax=100;
  vector<double> vec_double;
  for (int i=1; i<=Nmax; i++) {
    for (int j=i+1; j<=Nmax;j++) {
      double tmp= (1.0*i)/j;
      vec_double.push_back(tmp);
    }
  }
  vec_double.push_back(1.0);
  sort(vec_double.begin(), vec_double.end());
  remove_duplicates(vec_double);

  vector<double> double_range;
  for (uint i=0; i<vec_double.size();i++) {
    double num_diff = abs(vec_double.at(i)-dvalue);
    if(num_diff<=tol) {
      double_range.push_back(vec_double.at(i));
    }
  }

  vector<int> frac_tmp;
  vector<str_num_data> vec_str_num_data;
  str_num_data tmp_str_num_data;
  for (uint i=0; i<double_range.size();i++) {
    double db_tmp =  double_range.at(i);
    tmp_str_num_data = double2str_num_data(db_tmp);
    vec_str_num_data.push_back(tmp_str_num_data);
  }

  sort(vec_str_num_data.begin(), vec_str_num_data.end(), sortdenom); //sorting according to denom
  optimized_dvalue = vec_str_num_data.at(0); 


  return optimized_dvalue;
}

// ***************************************************************************
// void UpdateXstructure(xstructure &b)
// ***************************************************************************
void UpdateXstr_comp_each_type(xstructure &b) {
  //This function is only used to update partial occupation value i.e., comp_each_type
  vector<double> poccaus; // partial occupation local host
  for (uint i=0; i<b.atoms.size();i++) {
    double tmp=b.atoms.at(i).partial_occupation_value;
    poccaus.push_back(tmp);
  }

  b.comp_each_type.clear();
  for(uint i=0;i<b.num_each_type.size();i++) {
    b.comp_each_type.push_back(0.0);
  }

  for(uint i=0;i<b.atoms.size();i++) {
    if(poccaus.at(i)<1.0-1e-5) {
      b.atoms.at(i).partial_occupation_flag=TRUE;
      b.atoms.at(i).partial_occupation_value=poccaus.at(i);
      b.comp_each_type.at(b.atoms.at(i).type)+=b.atoms.at(i).partial_occupation_value;
    } else {
      b.atoms.at(i).partial_occupation_flag=FALSE;
      b.atoms.at(i).partial_occupation_value=1.0;
      b.comp_each_type.at(b.atoms.at(i).type)+=b.atoms.at(i).partial_occupation_value;
    }
  }
}


// ***************************************************************************
// void UpdateXstr(xstructure &xstr_orig, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
void UpdateXstr(xstructure &xstr_orig, ofstream &FileMESSAGE, _aflags &aflags) {
  double epsilon = 1E-2;
  //Update paritial occupation value
  if(xstr_orig.partial_occupation_flag) {
    ostringstream aus;
    aus << "0000 MESSAGE    \"partial_occupation_tol\" is " << xstr_orig.partial_occupation_tol << " " << Message(aflags, "user, host, time");
    if(xstr_orig.partial_occupation_tol>1E-6) {
        epsilon = xstr_orig.partial_occupation_tol;
    }
    else {
      aus << "0000 MESSAGE    \"partial_occupation_tol\" is not set " << Message(aflags, "user, host, time");
      aus << "0000 MESSAGE    Default value (1E-2) of \"partial_occupation_tol\" is set " << Message(aflags, "user, host, time");
    }
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    double partial_value_tmp;
    for (uint i=0; i<xstr_orig.atoms.size();i++) {
      partial_value_tmp = xstr_orig.atoms.at(i).partial_occupation_value;
      str_num_data tmp;
      tmp = OptimizePoccValue(partial_value_tmp, epsilon);
      xstr_orig.atoms.at(i).partial_occupation_value = tmp.decimal;
      UpdateXstr_comp_each_type(xstr_orig);
    }
  }
}

// ***************************************************************************
// Function pocc::HNFCELL
// ***************************************************************************
void OptimizeXstr(xstructure &a, ofstream &FileMESSAGE, _aflags &aflags) {
  //This function optimize partial values of xstructure
  ostringstream aus;
  aus <<"0000 MESSAGE    Optimizing POSCAR " << Message(aflags, "user, host, time");
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  UpdateXstr(a, FileMESSAGE, aflags);
  aus <<"0000 MESSAGE    Printing optimized POSCAR " << Message(aflags, "user, host, time");
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aus << a;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
}

// ***************************************************************************
// void string_replaceAll(std::string& str, const std::string& from, const std::string& to) {
// ***************************************************************************
void string_replaceAll(std::string& str, const std::string& from, const std::string& to) 
{
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

// ***************************************************************************
//  bool is_number(const std::string& s)
// ***************************************************************************
bool is_number(const std::string& s) {
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) {
      ++it;
  }
  return !s.empty() && it == s.end();
}

// ***************************************************************************
//  double CalculateDistanceXstructure(xstructure& xstr, int i, int j) {
// ***************************************************************************
double CalculateDistanceXstructure(xstructure& xstr, int i, int j) {
  double d1=xstr.atoms.at(i).cpos(1)-xstr.atoms.at(j).cpos(1);
  double d2=xstr.atoms.at(i).cpos(2)-xstr.atoms.at(j).cpos(2);
  double d3=xstr.atoms.at(i).cpos(3)-xstr.atoms.at(j).cpos(3);
  double d=sqrt(d1*d1+d2*d2+d3*d3);
  return d;
}

// ***************************************************************************
// vector<vector<int> > CalculateNumberMultiOccupied(xstructure& xstr)
// ***************************************************************************
vector<vector<int> > NormalisedNumberXstructure(xstructure& xstr) {
  //Format Xstructure, if two different atoms occupy the same site, then they will
  //be put in the same vector
  vector<vector<int> > a; 
  bool insert_flag;
  vector<int> temp;
  temp.push_back(0); //Put the first number of xstr
  a.push_back(temp);
  for (unsigned int i=1; i<xstr.atoms.size();i++) {
    int k1=i;
    insert_flag = false;
    for (unsigned int j=0; j<a.size();j++) {
      int k2=a.at(j).at(0);
      if(CheckOneSite(xstr,k1,k2)) {
	insert_flag = true;
	a.at(j).push_back(k1);
	break;
      }
    }
    if(! insert_flag) {
      temp.clear();
      temp.push_back(k1);
      a.push_back(temp);
    }
  }
  return a;
}

// ***************************************************************************
// bool CheckPartialOccupation(xstructure& xstr)
// ***************************************************************************
vector<vector<int> > CalculateLableXstructure(xstructure& xstr) {
  vector<vector<int> > a=NormalisedNumberXstructure(xstr);
  vector<vector<int> > NormalisedLable;
  vector<int> LabelTmp;
  for (unsigned int i=0; i<a.size();i++) {
    LabelTmp.clear();
    for (unsigned int j=0; j<a.at(i).size();j++) {
      int k=a.at(i).at(j);
      LabelTmp.push_back(xstr.atoms.at(k).type);
    }
    if(CheckVacancyOnOnesite(xstr,i)) {
      int NumVacancy=xstr.num_each_type.size();
      LabelTmp.push_back(NumVacancy);
    }
    //LabelTmp.push_back(-1); //Add "-1" to terminate
    NormalisedLable.push_back(LabelTmp);
  }
  return NormalisedLable;
}
// ***************************************************************************
// bool CheckPartialOccupation(xstructure& xstr)
// ***************************************************************************
bool CheckPartialOccupation(xstructure& xstr) {
  //if it is a partially occupied, then return true, else return false
  bool RunFLAG=false;
  if(xstr.partial_occupation_flag) {
    for (unsigned int i=0; i<xstr.atoms.size(); i++) {
      if(xstr.atoms.at(i).partial_occupation_flag) {
          RunFLAG=true;
      }
    }
  }
  return RunFLAG;
}

// ***************************************************************************
// bool CheckOneSite(xstructure& xstr, int i, int j)
// ***************************************************************************
bool CheckOneSite(xstructure& xstr, int i, int j) {
  //if xstr.atoms.at(i) and xstr.atoms.at(j) belong to the same site, return true, else return false
  double epsilon=1E-3;
  bool RunFLAG=false; 
  if(CalculateDistanceXstructure(xstr, i, j)<epsilon) {
      RunFLAG=true;
  }
  return RunFLAG;
}

// ***************************************************************************
// bool CheckMultiOccupied(xstructure& xstr)
// ***************************************************************************
bool CheckMultiOccupied(xstructure& xstr) {
  //if it is double occupied, then return true, else return false
  return(NormalisedNumberXstructure(xstr).size()<xstr.atoms.size()); 
}

// ***************************************************************************
// double CalculatePartialValueOfVacancy(xstructure& xstr, unsigned int k)
// ***************************************************************************
double CalculatePartialValueOfVacancy(xstructure& xstr, unsigned int k) {
  double SumTmp=0.0;
  vector<vector<int> > a=NormalisedNumberXstructure(xstr);
  if(k>=a.size()) {
    cerr << "Error! There must be wrong with your input! I can not go on!" << endl;
    exit(1);
  }
  if(a.at(k).size()==1) {
    int n_xstr_orig=a.at(k).at(0);
    SumTmp=xstr.atoms.at(n_xstr_orig).partial_occupation_value;
  } else{
    SumTmp=0.0;
    for (unsigned int i=0; i<a.at(k).size();i++) {
      int j=a.at(k).at(i);
      SumTmp+=xstr.atoms.at(j).partial_occupation_value;
    }
  }
  return(1.0-SumTmp);
}
// ***************************************************************************
// bool CheckVacancyOnOnesite(xstructure& xstr, int k)
// ***************************************************************************
bool CheckVacancyOnOnesite(xstructure& xstr, unsigned int k) {
  //Note that this k is the number is the normalised xstructure, not orignal xstructure
  bool RunFLAG=false;
  vector<vector<int> > a=NormalisedNumberXstructure(xstr);
  if(k>=a.size()) {
    cerr << "Error! Errors are found in your input file!" << endl;
    exit(1);
  }

  if(a.at(k).size()==1) {
    int n_xstr_orig=a.at(k).at(0);
    if(xstr.atoms.at(n_xstr_orig).partial_occupation_value<1.0) {
        RunFLAG=true;
    }
  }  else{
    double SumTmp=0.0;
    for (unsigned int i=0; i<a.at(k).size();i++) {
      int j=a.at(k).at(i);
      SumTmp+=xstr.atoms.at(j).partial_occupation_value;
    }
    if(SumTmp<1.0) {
        RunFLAG=true;
    }
  }
  return RunFLAG;
}
// ***************************************************************************
//bool CleanVacancy(xstructure& xstr) 
// ***************************************************************************
xstructure CleanVacancy(xstructure& xstr) 
{
  //This routine assumes that the vacancy is at the last type of the atoms
  xstructure CleanedXstr;
  int Nions=xstr.atoms.size();
  int Ntype=xstr.num_each_type.size();
  int NumLastType=xstr.num_each_type.at(Ntype-1);
  for (int i=(Nions-1); i>=(Nions-NumLastType);i--) {
    xstr.RemoveAtom(i);
  } 
  CleanedXstr=xstructure(xstr);
  return CleanedXstr;
}

// ***************************************************************************
//bool CheckVacancy(xstructure& xstr) 
// ***************************************************************************
bool CheckVacancy(xstructure& xstr) 
{
  //if xstr has a vacancy, then return true, else return false
  bool RunFLAG=false;
  if(CheckPartialOccupation(xstr)) {
    if(!CheckMultiOccupied(xstr)) {
        RunFLAG=true;
    }
    else {
      vector<vector<int> > a=NormalisedNumberXstructure(xstr);
      for (uint i=0;i<a.size();i++) {
	for (uint j=0; j<a.at(i).size();j++) {
	  cout << a.at(i).at(j) << " ";
	}
	cout << endl;
      }

      vector<double> SumOccupation;
      double SumTmp;
      for (unsigned int i=0; i<a.size();i++) {
	SumTmp=0.0;
	for (unsigned int j=0; j<a.at(i).size(); j++) {
	  int k=a.at(i).at(j);
	  SumTmp+=xstr.atoms.at(k).partial_occupation_value;
	}
	SumOccupation.push_back(SumTmp);
      }
      for(unsigned int i=0; i<SumOccupation.size();i++) { 
	if((1-SumOccupation.at(i))>1E-4) {
        RunFLAG=true;
    }
	else if(SumOccupation.at(i)>1) {
	  cerr << "Error! Please check your occupation value carefully  !!! " << endl;
	  cerr << "The total occupation on some site exceeds 1.0 !!!" << endl;
	  cerr << "I refuse to run this silly job!" << endl;
	  exit(1);
	}
      }
    }
  }
  return RunFLAG;
}

// ***************************************************************************
// vector<vector<int> > CalculatecRange(xstructure& xstr)
// ***************************************************************************
vector<vector<int> > CalculatecRange(xstructure& xstr) {
  int Ntype;
  if(CheckVacancy(xstr)) {
      Ntype=xstr.num_each_type.size()+1;
  }
  else {
      Ntype=xstr.num_each_type.size();
  }
  int nDFull=NormalisedNumberXstructure(xstr).size(); 
    
  vector<vector<int> > XstrcRange;
  vector<int> cRangeTmp;
  if(CheckPartialOccupation(xstr)) {
    double ctmp=0.0;
    vector<double>  ConUniqueSite;
    for (int i=0; i<(Ntype-1);i++) {
      ctmp+=xstr.comp_each_type.at(i);
      ConUniqueSite.push_back(xstr.comp_each_type.at(i)/nDFull);
    }
    ConUniqueSite.push_back((nDFull-ctmp)/nDFull);

    vector<int> FracTmp; //Midlle Step
    vector<int> ConFrac; //Get the fraction for each unique site
    vector<vector<int> > AllConFrac; //Combine all the fraction into vector
    for (unsigned int i=0; i<ConUniqueSite.size();i++) {
      FracTmp=Decimal2Fraction(ConUniqueSite.at(i));
      ConFrac=GetFraction(FracTmp);
      AllConFrac.push_back(ConFrac);
    }

    vector<int> denominatortmp;
    for (uint i=0; i<AllConFrac.size();i++) {
      denominatortmp.push_back(AllConFrac.at(i).at(1));
    }
    int Xstrlcm=CalculateLcmofVector(denominatortmp);

    int tmp;
    for (uint i=0; i<AllConFrac.size();i++) {
      cRangeTmp.clear();
      tmp=Xstrlcm/AllConFrac.at(i).at(1)*AllConFrac.at(i).at(0);
      cRangeTmp.push_back(tmp);
      cRangeTmp.push_back(tmp);
      cRangeTmp.push_back(Xstrlcm);
      XstrcRange.push_back(cRangeTmp);
    }
    /*The second algorithm
      int sum_num=0;
      for (uint i=0; i<XstrcRange.size();i++) {
      sum_num+=XstrcRange.at(i).at(0);
      }
      cRangeTmp.clear();
      cRangeTmp.push_back(XstrcRange.at(0).at(2)-sum_num);
      cRangeTmp.push_back(XstrcRange.at(0).at(2)-sum_num);
      cRangeTmp.push_back(XstrcRange.at(0).at(2));
      XstrcRange.push_back(cRangeTmp);
    */

    int sumtmp=0;
    for (uint i=0; i<XstrcRange.size();i++) {
      sumtmp+=XstrcRange.at(i).at(0); 
    }
        
    if(sumtmp!=XstrcRange.at(0).at(2)) { //sumtmp2 does not contain vacancy; sumtmp contains vacancy
      int sumtmp2=0;
      int k=XstrcRange.size();
      for (int i=0; i<(k-1);i++) {
	sumtmp2+=XstrcRange.at(i).at(0);
      }
      XstrcRange.at(k-1).at(0)=XstrcRange.at(0).at(2)-sumtmp2;
      XstrcRange.at(k-1).at(1)=XstrcRange.at(0).at(2)-sumtmp2;
    }
  }
  else {
    for (unsigned int i=0; i<xstr.num_each_type.size();i++) {
      cRangeTmp.clear();
      int k=xstr.num_each_type.at(i);
      int Nions=xstr.atoms.size();
      cRangeTmp.push_back(k);
      cRangeTmp.push_back(k);
      cRangeTmp.push_back(Nions);
      XstrcRange.push_back(cRangeTmp);
    }
  }
  return XstrcRange;
}

// ***************************************************************************
// int CalculateLcmofVector(vector<int> Num)
// ***************************************************************************
int CalculateLcmofVector(vector<int> Num) {
  vector<int> vlcm;
  int lcmtmp=0;
  if(Num.size()<2) {
    cerr << "Error ! Input values are less than 2! Its Least Common Multiple cannot be generated!" << endl;
    exit(1);
  }
  lcmtmp=lcm(Num.at(0),Num.at(1));
  vlcm.push_back(lcmtmp);
  for (unsigned int i=1; i<(Num.size()-1);i++) {
    lcmtmp=lcm(Num.at(i),Num.at(i+1));
    vlcm.push_back(lcmtmp);
  }
  lcmtmp = max(vlcm);
  return lcmtmp;
}

// ***************************************************************************
// int lcm (int a, int b)
// ***************************************************************************
int lcm (int a, int b) {
  return a/gcd(a,b)*b;
}
// ***************************************************************************
// int gcd(int a, int b)
// ***************************************************************************
int gcd(int a, int b) {
  if(a%b==0) {
      return b;
  }
  else {
      return gcd(b, a%b);
  }
}
// ***************************************************************************
// vector<int> Decimal2Fraction(double a)
// ***************************************************************************
vector<int> Decimal2Fraction(double a) {
  //This function converts decimal number to fraction number (not real),
  //You need to call GetFraction to get the real fraction
  vector<int> IntList;
  vector<int> FracTmp;
  int b;
  double DecNum=a; 
  //Debug
  double FracTmpDouble;
  double NumDiff=1;
  if(abs(a-int(a))<Tol_DecimalPart) {
    IntList.push_back(0);
    IntList.push_back(int(a));
  }
  else {
    while(abs(a)>Tol_DecimalPart && NumDiff>Tol_FracDeciDiff) {
      b=int(1.0/a);
      if(b>MaxDenominator) {
	break;
      }
      IntList.push_back(b);
      a=1.0/a-b;
      FracTmp=GetFraction(IntList);
      FracTmpDouble=double (FracTmp.at(0))/FracTmp.at(1);
      NumDiff=abs(DecNum-FracTmpDouble);
    }
  }
  return IntList;
}

// ***************************************************************************
// vector<int> GetFraction(vector<int> IntList)
// ***************************************************************************
vector<int> GetFraction(vector<int> IntList) {
  //This function read the vector generated by Decimal2Fraction(), and generate real fraction
  vector<int> Fraction;
  int LenghIntList=IntList.size();
  int tmp;
  int fac=1;
  int Denom;
  Denom=IntList.at(LenghIntList-1);
  for (int i=LenghIntList-2; i>=0; i--) {
    tmp=Denom;
    Denom=fac+Denom*IntList.at(i);
    fac=tmp;
  }
  Fraction.push_back(fac);
  Fraction.push_back(Denom);
  return Fraction;
}

// ***************************************************************************
// vector<int> double2fraction(double a)
// ***************************************************************************
vector<int> double2fraction(double a) {
  vector<int> IntList = Decimal2Fraction(a);
  vector<int> Fraction = GetFraction(IntList);
  return Fraction;
}

// ***************************************************************************
// bool sortdenom (str_num_data str_i, str_num_data str_j)
// ***************************************************************************
bool sortdenom (str_num_data str_i, str_num_data str_j) {
  return (str_i.denom<str_j.denom); //Ascending sort according to denom 
}


// ***************************************************************************
// str_num_data double2str_num_data(double a)
// ***************************************************************************
str_num_data double2str_num_data(double a) { 
  str_num_data result;
  vector<int> fraction = double2fraction(a);

  result.decimal = a;
  result.fac = fraction.at(0);
  result.denom = fraction.at(1);
  return result;
}

// ***************************************************************************
// pocc::POSCAR2GULP(istream& input)
// ***************************************************************************
namespace pocc {
  void POSCAR2GULP(istream& input) {
    xstructure xstr; 
    xstr.Clear();
    xstr=xstructure(input, IOVASP_POSCAR);
    xstr.ReScale(1.0);
    string GulpInput;
    GulpInput.clear();

    vector<string> AtomSpecies;
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");
    if(!xstr.atoms.at(0).name_is_given) {
      cerr << "Warnning!!! Atom Species are not assigned! They are defaulted to H." << endl;
      AtomSpecies.clear();
      for (uint i=0; i<xstr.num_each_type.size();i++) {
	AtomSpecies.push_back("H");
      }
    }

    GulpInput=pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    cout << GulpInput ;
  }
}

// ***************************************************************************
// pocc::POSCAR2GulpInput(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  string POSCAR2GulpInput(xstructure& xstr) {
    vector<string> AtomSpecies;
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");
    if(!xstr.atoms.at(0).name_is_given) {
      cerr << "Warnning!!! Atom Species are not assigned! They are defaulted to H." << endl;
      AtomSpecies.clear();
      for (uint i=0; i<xstr.num_each_type.size();i++) {
	AtomSpecies.push_back("H");
      }
    }

    string GulpInput = pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    return GulpInput;
  }
}

// ***************************************************************************
// pocc::POSCAR2GulpInput(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  string POSCAR2GulpInput(xstructure& xstr, vector<string> AtomSpecies) {
    if(xstr.num_each_type.size()!=AtomSpecies.size()) {
      cerr << "Error! Please check your names of each atom type!" << endl;
      exit(1);
    }
    //Assign names to xstructure 
    int iatom=0;
    for(uint itype=0; itype<xstr.num_each_type.size();itype++) {
      for (int j=0; j<xstr.num_each_type.at(itype);j++) {
	xstr.atoms.at(iatom).name=AtomSpecies.at(itype);
	xstr.atoms.at(iatom).name_is_given=TRUE;
	iatom++;
      }
    }

    stringstream  oss;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    oss << "single dist comp conp" << endl;
    oss << "cell" << endl;
    oss << xstr.a << '\t' << xstr.b << '\t'  << xstr.c << '\t' << xstr.alpha << '\t' << xstr.beta << '\t' << xstr.gamma << endl;
    oss << "fractional" << endl;
    oss.setf(ios_base::fixed,ios_base::floatfield);
    for (uint i=0; i < xstr.atoms.size(); i++) {
      oss <<xstr.atoms.at(i).name << '\t' << xstr.atoms.at(i).fpos(1) << '\t' << xstr.atoms.at(i).fpos(2) << '\t' << xstr.atoms.at(i).fpos(3) << endl;
    }
    oss << "#" << endl;
    oss << "species    " << xstr.num_each_type.size() << endl;

    string SpeciesTmp;
    for (uint i=0; i<AtomSpecies.size();i++) {
      SpeciesTmp.clear();
      SpeciesTmp=AtomSpecies.at(i);
      oss << pocc::ReturnAtomSpecies(SpeciesTmp) << endl;
    }
    oss << "#" << endl;
    oss << "uff1 kcal" << endl;
    for (uint i=0; i<AtomSpecies.size();i++) {
      SpeciesTmp.clear();
      SpeciesTmp=AtomSpecies.at(i);
      oss << pocc::ReturnAtomSpeciesPotential(SpeciesTmp) << endl;
    }
    return oss.str();
  }
}

// ***************************************************************************
// double pocc::CalculateEnergyUsingGulp(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  double CalculateEnergyUsingGulp(xstructure & xstr, vector<string> AtomSpecies) {
    double TotalEnergy=0.0;
    string XstrGulpInput;
    string GulpOut;
    string GulpTmp=aurostd::TmpFileCreate("GulpTmp");

    if(AtomSpecies.at(0).compare("unkown")==0) {
      cerr << "Warnning!!! Atom Species are not assigned!" << endl;
      exit(1);
    }

    XstrGulpInput=pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    aurostd::string2file(XstrGulpInput,GulpTmp);

    string RunGulp="cat "+ GulpTmp +" | gulp";
    GulpOut=aurostd::execute2string(RunGulp);
    aurostd::execute("rm " + GulpTmp);

    stringstream ss_gulp_out;
    ss_gulp_out << GulpOut;

    string anchor_word_Total_energy="Total lattice energy";
    string line;
    vector<string> vgulp_output_line, vgulp_output;
    while(getline(ss_gulp_out, line)) {
      if(line.find(anchor_word_Total_energy) !=string::npos) {
	vgulp_output.push_back(line);
      }
    }
    aurostd::string2tokens(vgulp_output.at(0), vgulp_output_line, " ");
    TotalEnergy=aurostd::string2utype<double>(vgulp_output_line.at(4));
    return TotalEnergy;
  }
} // namespace pocc

// ***************************************************************************
// double pocc::CalculateEnergyUsingGulp(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  double CalculateEnergyUsingGulp(xstructure& xstr) {
    double TotalEnergy=0.0;
    string XstrGulpInput;
    string GulpOut;
    string GulpTmp=aurostd::TmpFileCreate("GulpTmp");
    vector<string> AtomSpecies;
    AtomSpecies.clear();
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");

    XstrGulpInput=pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    aurostd::string2file(XstrGulpInput,GulpTmp);

    string RunGulp="cat "+ GulpTmp +" | gulp";
    GulpOut=aurostd::execute2string(RunGulp);
    aurostd::execute("rm " + GulpTmp);

    stringstream ss_gulp_out;
    ss_gulp_out << GulpOut;

    string anchor_word_Total_energy="Total lattice energy";
    string line;
    vector<string> vgulp_output_line, vgulp_output;
    while(getline(ss_gulp_out, line)) {
      if(line.find(anchor_word_Total_energy) !=string::npos) {
	vgulp_output.push_back(line);
      }
    }
    aurostd::string2tokens(vgulp_output.at(0), vgulp_output_line, " ");
    TotalEnergy=aurostd::string2utype<double>(vgulp_output_line.at(4));
    return TotalEnergy;
  }
} // namespace pocc


// ***************************************************************************
// double pocc::CalculateEnergyUsingGulp(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  vector<double> CalculateEnergyUsingGulp(vector<xstructure>& vxstr) {
    vector<double> Etot_vxstr; 
    for (uint i=0; i<vxstr.size();i++) {
      double dbtmp=0.0;
      xstructure xstr_tmp;
      xstr_tmp = vxstr.at(i);
      dbtmp=pocc::CalculateEnergyUsingGulp(xstr_tmp);
      Etot_vxstr.push_back(dbtmp);
    }
    return  Etot_vxstr;
  }
} // namespace pocc

void PrintGulpEnergies(vector<xstructure>& vxstr) {
  vector<double> etot_vxstr = pocc::CalculateEnergyUsingGulp(vxstr);
  for (uint i=0;i<etot_vxstr.size();i++) {
    ostringstream oss;
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(12);
    oss << "No. " << i+1 << " " << setw(15) << etot_vxstr.at(i) << endl;
    cout << oss.str();
  }
}

// ***************************************************************************
// bool sortenergy (strno_energy str_i, strno_energy str_j)
// ***************************************************************************
bool sortenergy (strno_energy str_i, strno_energy str_j) {
  return (str_i.energy<str_j.energy);
}

// ***************************************************************************
// vector<int> pocc::SortGroupXstr(vector<xstructure> groupxstr, vector<string> AtomSpecies)
// ***************************************************************************
namespace pocc {
  vector<xstructure> SortGroupXstr(vector<xstructure> groupxstr, vector<string> AtomSpecies) {
    vector<strno_energy> vec_group_strno_energy;
    vector<int> number;

    cerr.precision(8);
    cerr.setf(ios_base::fixed,ios_base::floatfield);

    double TotalEnergyTmp;
    for (uint i=0; i<groupxstr.size();i++) {
      strno_energy strno_energy_tmp;
      TotalEnergyTmp=pocc::CalculateEnergyUsingGulp(groupxstr.at(i),AtomSpecies);
      cerr << "Estimating Total Energy: Str # " <<  setfill('0') << setw(6) << i+1 << " : " << TotalEnergyTmp << " eV" << endl;
      strno_energy_tmp.energy=TotalEnergyTmp;
      strno_energy_tmp.number=i;
      vec_group_strno_energy.push_back(strno_energy_tmp);
    }
    sort (vec_group_strno_energy.begin(), vec_group_strno_energy.end(), sortenergy);

    cerr << "Sorting energy: " << endl;
    for (uint i=0; i<vec_group_strno_energy.size();i++) {
      cerr << "No. " << setfill('0') << setw(6) << i+1 << "   Str # " << setfill('0') << setw(6) << vec_group_strno_energy.at(i).number << " : " << vec_group_strno_energy.at(i).energy << " eV" << endl; 
      number.push_back(vec_group_strno_energy.at(i).number);
    }

    vector<xstructure> SortedGroupXstr;

    for (uint i=0; i<number.size();i++) {
      int k=number.at(i);
      SortedGroupXstr.push_back(groupxstr.at(k));
    }
    return SortedGroupXstr;
  }
} // namespace pocc

// ***************************************************************************
// vector<int> pocc::SortGroupXstrUFFEnergy(vector<xstructure> groupxstr, vector<string> AtomSpecies)
// ***************************************************************************
namespace pocc {
  vector<xstructure> SortGroupXstrUFFEnergy(vector<xstructure> groupxstr) {
    vector<strno_energy> vec_group_strno_energy;
    vector<int> number;

    cout.precision(8);
    cout.setf(ios_base::fixed,ios_base::floatfield);

    double TotalEnergyTmp;
    for (uint i=0; i<groupxstr.size();i++) {
      strno_energy strno_energy_tmp;
      TotalEnergyTmp=pocc::CalculateUFFEnergy(groupxstr.at(i));
      cout << "Estimating Total Energy: Str # " <<  setfill('0') << setw(6) << i+1 << " : " << TotalEnergyTmp << " eV" << endl;
      strno_energy_tmp.energy=TotalEnergyTmp;
      strno_energy_tmp.number=i;
      vec_group_strno_energy.push_back(strno_energy_tmp);
    }
    sort (vec_group_strno_energy.begin(), vec_group_strno_energy.end(), sortenergy);

    cout << "Sorting energy: " << endl;
    for (uint i=0; i<vec_group_strno_energy.size();i++) {
      cout << "No. " << setfill('0') << setw(6) << i+1 << "   Str # " << setfill('0') << setw(6) << vec_group_strno_energy.at(i).number << " : " << vec_group_strno_energy.at(i).energy << " eV" << endl; 
      number.push_back(vec_group_strno_energy.at(i).number);
    }

    vector<xstructure> SortedGroupXstr;

    for (uint i=0; i<number.size();i++) {
      int k=number.at(i);
      SortedGroupXstr.push_back(groupxstr.at(k));
    }
    return SortedGroupXstr;
  }
} // namespace pocc

#endif
// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG - Duke University 2010-2013              *
// *                                                                         *
// ***************************************************************************

