// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_APENNSY_3D_CPP
#define _AFLOW_APENNSY_3D_CPP

#include "aflow.h"
#include <map>

using std::ostringstream;
using std::string;
using std::ifstream;
using std::vector;
using std::cerr;
using aurostd::StringCommasColumsVectorInt;
using std::map;

#define PIC_POINTSIZE  string("3.0")
//#define PIC_SIZE     string("18,12")
#define PIC_SIZE       string("15,12")
#define PIC_FONTSIZE_DEFAULT 30
#define PIC_RATIO      0.866// sqrt(3.0)/2.0
bool TERDATA_COMPLETE=TRUE;// FALSE; // JESUS
#define epsilon 1E-6

//**************************************************************
bool checkBIG(vector<vector<double> >& b) {
  vector<double> tp;
  vector<vector<double> > a;
  for(uint i=0;i<b.size();i++) {
    tp.push_back(b.at(i).at(0));
    tp.push_back(b.at(i).at(1));
    a.push_back(tp);
    tp.clear();
  }			 
  sort(a.begin(),a.end());
  if(aurostd::abs(a.at(0).at(0))<epsilon && aurostd::abs(a.at(0).at(1))<epsilon &&
     aurostd::abs(a.at(1).at(0)-0.5)<epsilon && aurostd::abs(a.at(1).at(1)-0.866)<epsilon && 
     aurostd::abs(a.at(2).at(0)-1)<epsilon && aurostd::abs(a.at(2).at(1))<epsilon) return TRUE;
  return FALSE;
}

bool checkfinished(vector<double>& FX, vector<double>& FY,double x,double y) {
  if(FX.size()==0) {
    return FALSE;
  }
  for(uint i=0;i<FX.size();i++) {
    if(aurostd::abs(x-FX.at(i))<epsilon && aurostd::abs(y-FY.at(i))<epsilon) return TRUE;
  }
  return FALSE;
}

double dotproforter(const vector<double>& a,const vector<double>& b) {
  double c=0;
  for(uint i=0;i<a.size();i++) {
    c+=a.at(i)*b.at(i);
    //    cerr << a.at(i)*b.at(i) << endl;
  }
  //*****DEBUG*****
  if(0) {
    for(uint i=0;i<a.size();i++) cerr << a.at(i) << "  ";
    cerr << endl;
    for(uint i=0;i<b.size();i++) cerr << b.at(i) << "  ";
    cerr << endl;
  }
  return c;
}

void NumReduction(vector<double>& a) {
  vector<string> data;
  double divider;
  string primelist="2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61";
  aurostd::string2tokens(primelist,data,",");
  bool nextflag=FALSE;
  for(uint i=0;i<data.size();) {
    divider=aurostd::string2utype<double>(data.at(i));
    if(nextflag) {
      i++;
      nextflag=FALSE;
      continue;
    }
    bool breakflag=TRUE;
    for(uint j=0;j<a.size() && breakflag;j++) {
      if(int(a.at(j))%int(divider)==0) continue;
      else breakflag=FALSE;
    }
    if(breakflag) {
      for(uint j=0;j<a.size();j++) a.at(j)=a.at(j)/divider;
    }
    else nextflag=TRUE;
  }
}

void PointsCalculatefor3D(vector<double>& inp,vector<double>& pha, vector<double>& cat) {
  double t, a,b,c,e,AB,BC,CA,X,Y;
  a=inp.at(0);b=inp.at(1);c=inp.at(2);e=inp.at(3);t=a+b+c;
  AB=a/t;BC=b/t;CA=c/t;
  pha.push_back(AB);pha.push_back(CA);pha.push_back(BC);pha.push_back(e);
  X=CA+BC/2;Y=BC*sqrt(3)/2;
  cat.push_back(X);cat.push_back(Y);cat.push_back(e);
}

void FindGNDState(vector<string>& a, double& b) {
  double enthalpy_min=100000.0,enthalpy;
  for(uint i=0;i<a.size();i++) {
    vector<string> SB;
    aurostd::string2tokens(a.at(i),SB,"|");
    for(uint j=0;j<SB.size();j++) {
      if(aurostd::substring2bool(SB.at(j),"enthalpy_atom=")) {
	enthalpy=aurostd::substring2utype<double>(SB.at(j),"enthalpy_atom=");
	//	cerr << enthalpy_min << endl;
	if(enthalpy_min > enthalpy) {
	  enthalpy_min=enthalpy;
	}
      }
    }
  }
  b=enthalpy_min;
}
 
bool CheckFacet(vector<vector<double> >& a, vector<double>& b) {
  xvector<double> v1(3),v2(3),v3(3),tp(3);
  v3(1)=0;v3(2)=0;v3(3)=1;
  for(uint i=1;i<3;i++) {
    v1(i)=a.at(0).at(i-1)-a.at(1).at(i-1);
    v2(i)=a.at(1).at(i-1)-a.at(2).at(i-1);
  }
  //  v1(4)=b.at(0)-b.at(1);   v2(4)=b.at(1)-b.at(2);
  tp=aurostd::vector_product(v1,v2);
  if(aurostd::abs(aurostd::scalar_product(tp,v3)) < 0.001) return TRUE;
  return FALSE;
}
 
double normc(xvector<double>& a) {return sqrt(a(1)*a(1)+a(2)*a(2));}
double dpc(xvector<double>& a,xvector<double>& b) { return (a(1)*b(1)+a(2)*b(2)); }

void TriangleForm(vector<string>& PURE_SPE, xvector<double>& axis, string& name) {
  vector<string> Species;
  vector<double> Num;
  int t=0;
  xvector<double> A(3),B(3),C(3);

  A(1)=1.0;A(2)=0;A(3)=0;
  B(1)=0;B(2)=1.0;B(3)=0;
  C(1)=0;C(2)=0;C(3)=1.0;

  XATOM_SplitAlloySpecies(name,Species,Num);
  for(uint x=0;x<Num.size();x++) t+=Num.at(x);
  for(uint x=0;x<Num.size();x++) Num.at(x)=Num.at(x)/t;
  if(Species.size()==1) {
    if(Species.at(0)==PURE_SPE.at(0)) axis=A;
    if(Species.at(0)==PURE_SPE.at(1)) axis=B;
    if(Species.at(0)==PURE_SPE.at(2)) axis=C;
  }
  if(Species.size()==2) {
    if(Species.at(0)==PURE_SPE.at(0)) {
      if(Species.at(1)==PURE_SPE.at(1)) axis=(Num.at(0)*A+Num.at(1)*B);
      if(Species.at(1)==PURE_SPE.at(2)) axis=(Num.at(0)*A+Num.at(1)*C);
    }
    if(Species.at(0)==PURE_SPE.at(1) && Species.at(1)==PURE_SPE.at(2)) axis=(Num.at(0)*B+Num.at(1)*C);
  }
  if(Species.size()==3) axis=(Num.at(0)*A+Num.at(1)*B+Num.at(2)*C);
  // cerr << axis(1) << " " << axis(2) << " " << axis(3) << endl;
  if(Species.size()==0) {
    cerr << "ERROR: NO ELEMENT FOUND IN DATABASE!" << endl;
    exit(0);
  }
}

//**************************PHASE POINTS*****************************

class PhasePoints{
public:
  PhasePoints();
  PhasePoints(string& a, double& b, int& c, xvector<double>& d,double& f, string& g,string& h);
  PhasePoints(string& a, double& b, int& c, xvector<double>& d, vector<double>& e,double& f, string& g,bool flag);
  string CompoundName,protos,Dirs;
  int ID,DIM;
  double DisttoPlane,Ts;
  vector<double> Cart, Phase;
  bool sortingflag;
};

PhasePoints::PhasePoints() {
  ID=0;
  DIM=0;
  protos="";
  Ts=0.0;
  DisttoPlane=0.0;
  Dirs="";
  CompoundName="NONE";
  ID=-1;
}

PhasePoints::PhasePoints(string& a,double& b, int& c, xvector<double>& d,double& f, string& g, string& h) {
  vector<string> Species;
  vector<double> tp,Num;
  CompoundName=a;
  Dirs=h;
  ID=c;
  protos=g;
  Ts=f;
  XATOM_SplitAlloySpecies(a,Species,tp);
  DIM=Species.size();
  for(uint i=1;i<4;i++) Num.push_back(d(i));
  Num.push_back(b);
  PointsCalculatefor3D(Num,Phase,Cart);
  DisttoPlane=-1.0;
  sortingflag=FALSE;
}

PhasePoints::PhasePoints(string& a,double& b, int& c, xvector<double>& d, vector<double>& e,double& f,string& g,bool flag) {
  vector<string> Species;
  vector<double> tp,Num;
  CompoundName=a;
  ID=c;
  Ts=f;
  protos=g;
  XATOM_SplitAlloySpecies(a,Species,tp);
  for(uint i=1;i<4;i++) Num.push_back(d(i));
  //  cerr << b << "  ";
  if(flag) b=b-Num.at(0)*e.at(0)-Num.at(1)*e.at(1)-Num.at(2)*e.at(2);
  //  cerr << b << endl;
  DIM=Species.size();
  Num.push_back(b);
  PointsCalculatefor3D(Num,Phase,Cart);
  DisttoPlane=-1.0;
  sortingflag=FALSE;
}

//*************************CONVEX SURFACE************************
class Surface{
public:
  Surface();
  Surface(vector<PhasePoints>& a, int& b, int c);
  bool INSIDEpoint(PhasePoints& a);
  double CalculateDIST(PhasePoints& b);
  vector<PhasePoints> corners;
  xvector<double> normalv;
  int ID,Dim;
};

Surface::Surface() {
  ID=-1;
}

Surface::Surface(vector<PhasePoints>& a, int& b, int c) {
  for(uint i=0;i<a.size();i++) corners.push_back(a.at(i));
  Dim=c;
  xvector<double> a12(3),a23(3);
  a12(1)=corners.at(0).Cart.at(0)-corners.at(1).Cart.at(0);a12(2)=corners.at(0).Cart.at(1)-corners.at(1).Cart.at(1);a12(3)=corners.at(0).Cart.at(2)-corners.at(1).Cart.at(2);
  a23(1)=corners.at(1).Cart.at(0)-corners.at(2).Cart.at(0);a23(2)=corners.at(1).Cart.at(1)-corners.at(2).Cart.at(1);a23(3)=corners.at(1).Cart.at(2)-corners.at(2).Cart.at(2);
  normalv(1)=a12(2)*a23(3)-a23(2)*a12(3);  normalv(2)=a12(3)*a23(1)-a23(3)*a12(1);  normalv(3)=a12(1)*a23(2)-a23(1)*a12(2);
  //  cerr << normalv << endl;
  ID=b;
}

bool Surface::INSIDEpoint(PhasePoints& a) {
  xvector<double> v1(2),v2(2),v3(2);
  double a1,a2,a3;
  v1(1)=a.Cart.at(0)-corners.at(0).Cart.at(0);v1(2)=a.Cart.at(1)-corners.at(0).Cart.at(1);
  v2(1)=a.Cart.at(0)-corners.at(1).Cart.at(0);v2(2)=a.Cart.at(1)-corners.at(1).Cart.at(1);
  v3(1)=a.Cart.at(0)-corners.at(2).Cart.at(0);v3(2)=a.Cart.at(1)-corners.at(2).Cart.at(1);
  a1=acos(dpc(v1,v2)/normc(v1)/normc(v2));a2=acos(dpc(v2,v3)/normc(v2)/normc(v3));a3=acos(dpc(v1,v3)/normc(v1)/normc(v3));
  if(aurostd::abs((a1+a2+a3)/PI*180-360)< 0.01) return TRUE;
  return FALSE;
}

double Surface::CalculateDIST(PhasePoints& b) {
  double pinplane;
  if(aurostd::abs(normalv(1))<epsilon && aurostd::abs(normalv(2))<epsilon && aurostd::abs(normalv(3))<epsilon) return (corners.at(0).Cart.at(2)-b.Cart.at(2));
  pinplane=(normalv(1)*(b.Cart.at(0)-corners.at(0).Cart.at(0))+normalv(2)*(b.Cart.at(1)-corners.at(0).Cart.at(1)))/normalv(3)+corners.at(0).Cart.at(2);
  return (b.Cart.at(2)-pinplane);
}

//**************************************************CLASS END**********************************

vector<double> Faxiangliang(vector<PhasePoints>& a) {
  vector<double> b,ab,bc;
  for(uint i=0;i<a.at(0).Cart.size();i++) {
    ab.push_back(a.at(0).Cart.at(i)-a.at(1).Cart.at(i));
    bc.push_back(a.at(1).Cart.at(i)-a.at(2).Cart.at(i));
  }
  //****DEBUG*******
  if(0) {
    for(uint i=0;i<a.size();i++) {
      for(uint j=0;j<a.at(i).Cart.size();j++) {
 	cerr << "(" << i << "," << j << ") " << a.at(i).Cart.at(j) << "  ";
      }
      cerr << endl;
    }
    cerr << endl;
    for(uint i=0;i<bc.size();i++) cerr << bc.at(i) << "  "; 
    cerr << endl;
    for(uint i=0;i<ab.size();i++) cerr << ab.at(i) << "  "; 
    cerr << endl << endl;
  }
  //*****CORSSPRODUCT*****
  b.push_back(ab.at(1)*bc.at(2)-ab.at(2)*bc.at(1));
  b.push_back(ab.at(2)*bc.at(0)-ab.at(0)*bc.at(2));
  b.push_back(ab.at(0)*bc.at(1)-ab.at(1)*bc.at(0));
  return b;
}

double GetFormation(string& add, string& title) {  //map<string,double>& pure_list,string& title) {
  ostringstream oss,ossclean;
  ifstream inputfile;
  string line;
  double enthalpy_atom, FORMATION,ANUM;
  vector<string> spe;
  vector<double> nmv; 
  //cerr << add+"OUTCAR.static.bz2" << endl;
  //cerr << aurostd::FileExist(add+"/OUTCAR.static.bz2") << endl;;
  if(aurostd::FileExist(add+"/OUTCAR.static.bz2")) oss << XHOST.command("bzcat") << " " << add << "/OUTCAR.static.bz2 > phase_tp.dat";
  else{
    if(aurostd::FileExist(add+"/OUTCAR.static")) oss << "cat " << add << "/OUTCAR.static > phase_tp.dat";
    else{
      cerr << "ERROR: NO OUTCAR IN DIR! " << add << " WILL NOT BE INVOLVED!";
      return 1000.0;
    }
  }
  aurostd::execute(oss);
  inputfile.open("phase_tp.dat");
  while(!inputfile.eof()) {
    vector<string> data;
    getline(inputfile,line);
    //    cerr << line << endl;
    aurostd::string2tokens(line,data," ");
    //cerr << line << endl;
    if(data.size() > 5 && data.at(0)=="free") enthalpy_atom=aurostd::string2utype<double>(data.at(4));
  }
  inputfile.close();
  ossclean << "rm -f phase_tp.dat";
  aurostd::execute(ossclean);
  XATOM_SplitAlloySpecies(title,spe,nmv);
  FORMATION=enthalpy_atom;
  ANUM=0;
  for(uint i=0;i<nmv.size();i++) {
    //    FORMATION-=nmv.at(i)*pure_list[spe.at(i)];
    ANUM+=nmv.at(i);
  }
  FORMATION=FORMATION/ANUM;
  //  cerr << "ENTHALPY_ATOM" << FORMATION << endl;
  return FORMATION;
}

//**************************************************************

void INPUTDATAFORTERPHASE_NEW(vector<string>& vspecies, vector<PhasePoints>& POINTS,int& TOTAL, vector<PhasePoints>& TERS, PhasePoints& maxP);
void INPUTDATAFORTERPHASE_NEW(vector<string>& vspecies, vector<PhasePoints>& POINTS,int& TOTAL, vector<PhasePoints>& TERS) {
  PhasePoints maxP;
  INPUTDATAFORTERPHASE_NEW(vspecies,POINTS,TOTAL,TERS,maxP);
}

void INPUTDATAFORTERPHASE_NEW(vector<string>& vspecies, vector<PhasePoints>& POINTS,int& TOTAL, vector<PhasePoints>& TERS, PhasePoints& maxP) {
  ostringstream oss;
  int uID=0;
  bool LVERBOSE=TRUE; 
  string grep="";
  for(uint i=0;i<vspecies.size();i++) grep=grep+vspecies.at(i)+((i<vspecies.size()-1)?"|":"");
  aflowlib::LOAD_Library_ALL("",grep,LVERBOSE);
  vector<string> vspecies_pp;
  vector<vector<string> > vList,RAW_list;
  
  map<string,double>::iterator it;
  
  // for(uint i=0;i<vspecies.size();i++) cerr << "vspecies.at(i)=" << vspecies.at(i) << endl;
  // for(uint i=0;i<vspecies_pp.size();i++) cerr << "vspecies_pp.at(i)=" << vspecies_pp.at(i) << endl;
  if(vspecies.size()!=3) {
    cerr << "ERROR: INPUTDATAFORTERPHASE_NEW is made for ternaries..... exiting.." << endl; exit(0);
  }
  aflowlib::GREP_Species_ALL(vspecies,vspecies_pp,vList);
  TOTAL=int(vList.at(0).size()+vList.at(1).size()+vList.at(2).size());
  cout << "00000  MESSAGE - INPUTDATAFORTERPHASE_NEW: 1ary vList.at(0).size()=" << vList.at(0).size() << endl;
  cout << "00000  MESSAGE - INPUTDATAFORTERPHASE_NEW: 2ary vList.at(1).size()=" << vList.at(1).size() << endl;
  cout << "00000  MESSAGE - INPUTDATAFORTERPHASE_NEW: 3ary vList.at(2).size()=" << vList.at(2).size() << endl;
  // for(uint i=0;i<vList.at(1).size();i++) cout << vList.at(1).at(i) << endl;  exit(0);

  //********************PURE*************************
  string SEP,name,Protos,DIRS,LDAU;
  double enthalpy=0.0,TS=0.0;
  PhasePoints tmp;
  vector<PhasePoints> Min;
  //   vector<double> Min_e;
  vector<int> ID;
  for(uint i=0;i<vspecies.size();i++) {
    PhasePoints ta;
    name=vspecies[i]+"1";
    TS=0.0;
    Protos="";
    DIRS="";
    xvector<double> axis(3);
    TriangleForm(vspecies,axis,name);
    ta=PhasePoints(name,enthalpy,uID,axis,TS,Protos,DIRS);
    uID++;
    POINTS.push_back(ta);
  }
  
  for(uint i=0;i<vList.at(1).size();i++) {
    aflowlib::_aflowlib_entry aflowlib_tmp;
    aflowlib_tmp.Load(vList.at(1).at(i),cerr);
    name=aflowlib_tmp.compound;
    // cout << name << endl;
    DIRS=aflowlib_tmp.aurl;
    if(aurostd::substring2bool(vList.at(1).at(i),"enthalpy_formation_atom=")) {
      enthalpy=aflowlib_tmp.enthalpy_formation_atom; //flag=FALSE;
    } else {
      if(aurostd::substring2bool(vList.at(1).at(i),"ldau_TLUJ=")) {
	cerr << "WARNING: IGNORE THE COMPOUND (LDAU) FOR " << name << " DIR=[" << DIRS << "]" << endl;
	continue;
      } else{
	cerr << "WARNING: IGNORE THE COMPOUND (NO ENTHALPY_FORMATION) FOR  " << name << " DIR=[" << DIRS << "]" << endl;
	continue;
      }
      // enthalpy=aurostd::string2utype<double>(aflowlib::TokenExtractAFLOWLIB(vList.at(1).at(i),"enthalpy_atom="));
      // flag=TRUE;
    }
    Protos=aflowlib_tmp.prototype; 
    TS=aflowlib_tmp.entropic_temperature;
    xvector<double> axis(3);
    TriangleForm(vspecies,axis,name);
    PhasePoints tppt;
    tppt=PhasePoints(name,enthalpy,uID,axis,TS,Protos,DIRS);
    if(tppt.Cart.at(2)<0) {
      uID++;
      POINTS.push_back(tppt);
      if(tppt.Ts > maxP.Ts) maxP=tppt;
    }  else continue;
  }
  
  //******************************DIST ADDED********************
  //******************************GET MINIMUM ENTHALPY************
  
  vector<PhasePoints> Sorting;
  for(uint i=0;i<vList.at(2).size();i++) {
    aflowlib::_aflowlib_entry aflowlib_tmp;
    aflowlib_tmp.Load(vList.at(2).at(i),cerr);
    name=aflowlib_tmp.compound;
    Protos=aflowlib_tmp.prototype;
    TS=aflowlib_tmp.entropic_temperature;
    //   if(aurostd::substring2bool(Protos,"_ICSD_")) {name=Protos;aurostd::StringSubst(name,"_ICSD_","=");}
    
    cerr << "name=" << name << "  proto=" << Protos << "  ";
    // cerr << aflowlib::TokenExtractAFLOWLIB(vList.at(2).at(i),"enthalpy_atom=") << "  ";
    DIRS=aflowlib_tmp.aurl;
    if(aurostd::substring2bool(vList.at(2).at(i),"enthalpy_formation_atom=")) {
      enthalpy=aflowlib_tmp.enthalpy_formation_atom;
      //flag=FALSE;
    } else {
      // calculate your way
      if(aurostd::substring2bool(vList.at(2).at(i),"ldau_TLUJ=")) {
	cerr << "WARNING: IGNORE THE COMPOUND (LDAU) FOR " << name << " DIR=[" << DIRS << "]" << endl;
	continue;
      } else {
	cerr << "WARNING: IGNORE THE COMPOUND (NO ENTHALPY_FORMATION) FOR  " << name << " DIR=[" << DIRS << "]" << endl;
	continue;
      }
      // flag=TRUE;
      // enthalpy=aurostd::string2utype<double>(aflowlib::TokenExtractAFLOWLIB(vList.at(2).at(i),"enthalpy_atom="));
    }
    cerr << "  enthalpy_formation_atom=" << enthalpy;
    cerr << "  pearson=" << aflowlib_tmp.Pearson_symbol_relax;
    cerr << endl;

    xvector<double> axis(3);
    TriangleForm(vspecies,axis,name);
    PhasePoints tppt;
    tppt=PhasePoints(name,enthalpy,uID,axis,TS,Protos,DIRS);
    //cerr << tppt.Cart.at(2) << " ";
    if(tppt.Cart.at(2) < 0) { // stefano
      uID++;
      Sorting.push_back(tppt);
      if(tppt.Ts > maxP.Ts) maxP=tppt;
    } else continue;  // stefano
    //cerr << tppt.CompoundName << endl;
  }
  
  for(uint i=0;i<Sorting.size();i++) {
    TERS.push_back(Sorting[i]);
    double MinSe=Sorting.at(i).Cart.at(2);
    PhasePoints Mintp=Sorting.at(i);
    if(! Sorting.at(i).sortingflag) Sorting.at(i).sortingflag=TRUE;
    else continue;
    for(uint j=i+1;j<Sorting.size();j++) {
      if(! Sorting.at(j).sortingflag) {
 	if(aurostd::abs(Mintp.Cart.at(0)-Sorting.at(j).Cart.at(0))<0.01) {
 	  if(aurostd::abs(Mintp.Cart.at(1)-Sorting.at(j).Cart.at(1))<0.01) {
 	    Sorting.at(j).sortingflag=TRUE;
 	    if(MinSe <  Sorting.at(i).Cart.at(2)) {
	      MinSe=Sorting.at(j).Cart.at(2);
	      Mintp=Sorting.at(j);
 	    } 
 	    else continue;
 	  }
 	  else continue;
 	}
 	else continue;
      }
      else continue;
    }
    POINTS.push_back(Mintp); 
  }
}

void TERNARY_CONVEX(vector<PhasePoints>& POINTS, vector<vector<vector<double> > >& LINE2D,vector<vector<double> >& ENTHALPY) {
  //  sort (POINTS.begin(),POINTS.end(),SortPoints);
  
  //**************************INPUT FOR QHULL***********************
  //ostringstream outfile;
  string qhull_in=aurostd::TmpFileCreate("TERNARY_qhull_in");
  string qhull_out=aurostd::TmpFileCreate("TERNARY_qhull_out");
  stringstream oss,oss1,oss2,oss3;
  oss << "3 PhaseDiagram " << POINTS.size() << endl << POINTS.size() << endl;
  for(uint i=0;i<POINTS.size();i++) {
    if(POINTS.at(i).Cart.at(2) > 0) continue;
    for(uint j=0;j<POINTS.at(i).Cart.size();j++) {
      oss << POINTS.at(i).Cart.at(j) << " ";
    }
    oss << endl;
  }
  aurostd::stringstream2file(oss,qhull_in);
  oss.clear();
  oss1 << "qconvex G < " << qhull_in << " > " << qhull_out;
  aurostd::execute(oss1);
  oss1.clear();
  oss2 << "rm -f " << qhull_in;
  aurostd::execute(oss2);
  oss2.clear();

  //***********************READ QHULL RESULTS**********************
 
  stringstream infile;
  string line,c;
  bool readflag=FALSE;
  vector<double> conten,enthalpy;
  vector<vector<double> > CONTEN;
  aurostd::file2stringstream(qhull_out,infile);
  while(!infile.eof()) {
    getline(infile,line);
    if(line.substr(0,5)=="{ OFF") {
      readflag=TRUE;
      continue;
    }

    if(readflag) {
      vector<string> tmp;
      aurostd::string2tokens(line,tmp);
      if(tmp.size()==3) {
	conten.push_back(aurostd::string2utype<double>(tmp.at(0)));
	conten.push_back(aurostd::string2utype<double>(tmp.at(1)));
	CONTEN.push_back(conten);
	enthalpy.push_back(aurostd::string2utype<double>(tmp.at(2)));
	conten.clear();
      } else{
	readflag=FALSE;
	LINE2D.push_back(CONTEN);
	ENTHALPY.push_back(enthalpy);
	enthalpy.clear();
	conten.clear();
	CONTEN.clear();
      }
    }
  }
  oss3 << "rm -f " << qhull_out << endl;;
  aurostd::execute(oss3);
  oss3.clear();
}

void GetListName(vector<vector<vector<double> > >& LINE2D, vector<vector<double> >& ENTHALPY,vector<PhasePoints>& POINTS,ostream& oss, PhasePoints& maxP) {
  stringstream outdat;
  vector<double> X,Y,FX,FY,CX,CY;
  map<string,PhasePoints> lists;
  //for(uint i=0;i<POINTS.size();i++) {
  //  lists[POINTS[i].CompoundName]=-1;
  //}
  for(uint i=0;i<LINE2D.size();i++) {
    if(!CheckFacet(LINE2D.at(i),ENTHALPY.at(i))) {
      vector<PhasePoints> tplist;
      if(!checkBIG(LINE2D.at(i))) {
	for(uint j=0;j<LINE2D.at(i).size();j++) {
	  for(unsigned int k=0;k<POINTS.size();k++) {
	    //	if(arr[k]==FALSE) continue;
	    if(aurostd::abs(ENTHALPY.at(i).at(j)-POINTS.at(k).Cart.at(2))<0.01) {
	      if(aurostd::abs(LINE2D.at(i).at(j).at(0)-POINTS.at(k).Cart.at(0))<0.01) {
		if(aurostd::abs(LINE2D.at(i).at(j).at(1)-POINTS.at(k).Cart.at(1))<0.01) {
		  if(POINTS[k].DIM==3) {
		    // arr[k]=FALSE;
		    // cerr << "TEST: " << POINTS[k].CompoundName << " " << POINTS[k].Ts << " " << lists[POINTS[k].CompoundName] << endl;
		    if((lists[POINTS[k].CompoundName]).Ts < POINTS[k].Ts) lists[POINTS[k].CompoundName]=POINTS[k];
		    break;
		  }
		}
	      }
	    }
	  }
	}
      }
    }  
  }
  map<string,PhasePoints>::iterator it;
  for(it=lists.begin();it!=lists.end();it++) {
    if((it->second).Ts > 0) 
      oss << aurostd::PaddedPRE(it->first,15) << aurostd::PaddedPRE((it->second).Cart[2],15) << aurostd::PaddedPRE((it->second).Ts,10) 
	  << aurostd::PaddedPRE((it->second).protos,15) << aurostd::PaddedPRE(maxP.CompoundName+" "+aurostd::utype2string(maxP.Ts),30)
	  << aurostd::PaddedPRE((it->second).Dirs,70) << endl;
  }

}

void DigramPlot3ary(vector<vector<vector<double> > >& LINE2D, vector<vector<double> >& ENTHALPY,const vector<string>& vspecies,vector<PhasePoints>& POINTS,int& TOTAL, vector<PhasePoints>& TERS,int PIC_FONTSIZE) {
  if(vspecies.size()!=3) {
    cerr << "ERROR: DigramPlot3ary is made for ternaries..... exiting.." << endl; exit(0);
  }

  stringstream outdat;
  vector<double> X,Y,FX,FY,CX,CY;
  //  vector<uint> CX,CY;
  //outdat.open("tmp.dat");
  string gnu_dat=aurostd::TmpFileCreate("gnu_dat");

  //  for(uint i=0;i<POINTS.size();i++) cout << POINTS.at(i).CompoundName << endl;
  outdat << "0 0" << endl << "0.5 0.866" << endl << endl;  
  outdat << "1 0" << endl << "0.5 0.866" << endl << endl;
  outdat << "1 0" << endl << "0 0 " << endl << endl;
  for(uint i=0;i<LINE2D.size();i++) {
    if(!CheckFacet(LINE2D.at(i),ENTHALPY.at(i))) {
      for(uint j=0;j<LINE2D.at(i).size();j++) {
	if(aurostd::abs(ENTHALPY.at(i).at(j))>0.001) {
	  if(!(aurostd::abs(LINE2D.at(i).at(j).at(0))<epsilon && aurostd::abs(LINE2D.at(i).at(j).at(1))<epsilon)) {
	    if(!(aurostd::abs(LINE2D.at(i).at(j).at(0)-0.5)<epsilon && aurostd::abs(LINE2D.at(i).at(j).at(1)-0.866)<epsilon)) {
	      if(!(aurostd::abs(LINE2D.at(i).at(j).at(0)-1)<epsilon && aurostd::abs(LINE2D.at(i).at(j).at(1))<epsilon)) {
		X.push_back(LINE2D.at(i).at(j).at(0));
		Y.push_back(LINE2D.at(i).at(j).at(1));
		CX.push_back(i);
		CY.push_back(j);
	      }
	    }
	  }
	}
      }
    }
  }
  

  for(uint i=0;i<LINE2D.size();i++) {
    if(!CheckFacet(LINE2D.at(i),ENTHALPY.at(i))) {
      for(uint j=0;j<LINE2D.at(i).size();j++) {
	if(j==LINE2D.at(i).size()-1) {
	  //if(CheckAngle(LINE2D.at(i).at(0),LINE2D.at(i).at(j))) {
	  outdat << LINE2D.at(i).at(j).at(0) << " " << LINE2D.at(i).at(j).at(1) << endl << LINE2D.at(i).at(0).at(0) << " " << LINE2D.at(i).at(0).at(1) << endl << endl;
	} else {
	  //	    if(CheckAngle(LINE2D.at(i).at(j),LINE2D.at(i).at(j+1))) {
	  outdat << LINE2D.at(i).at(j).at(0) << " " << LINE2D.at(i).at(j).at(1) << endl << LINE2D.at(i).at(j+1).at(0) << " " << LINE2D.at(i).at(j+1).at(1) << endl << endl;
	  //}
	}
      }
    }
  }
  
  aurostd::stringstream2file(outdat,gnu_dat);
 
  //**********************CONSTRUCT DATA POINTS*****************************************
  vector<Surface> SURFACE;
  Surface tpsur;
  int recnum=0;

  for(uint i=0;i<LINE2D.size();i++) {
    if(CheckFacet(LINE2D.at(i),ENTHALPY.at(i))) {continue;}
    vector<PhasePoints> tplist;
    if(checkBIG(LINE2D.at(i))) continue;
    for(uint j=0;j<LINE2D.at(i).size();j++) {
      for(unsigned int k=0;k<POINTS.size();k++) {
	if(aurostd::abs(ENTHALPY.at(i).at(j)-POINTS.at(k).Cart.at(2))<0.01) {
	  //cerr << LINE2D.at(i).at(j).at(0) << " " << LINE2D.at(i).at(j).at(1) << " " << ENTHALPY.at(i).at(j) << "  " << POINTS.at(k).Cart.at(0) << " " << POINTS.at(k).Cart.at(1) << " " << POINTS.at(k).Cart.at(2) << endl;
	  if(aurostd::abs(LINE2D.at(i).at(j).at(0)-POINTS.at(k).Cart.at(0))<0.01) {
	    if(aurostd::abs(LINE2D.at(i).at(j).at(1)-POINTS.at(k).Cart.at(1))<0.01) {
	      tplist.push_back(POINTS.at(k));
	      break;
	    }
	  }
	}
      }
    }
    //****DEBUG***
    if(0) {
      cerr << "LINE: ";
      cerr << LINE2D.at(i).size() << endl;
      for(uint sa=0;sa<LINE2D.at(i).size();sa++) { cerr << LINE2D.at(i).at(sa).at(0) << " " << LINE2D.at(i).at(sa).at(1) << " " << ENTHALPY.at(i).at(sa) << endl;}
      cerr << "TPLI: " << tplist.size() << endl;;
      for(uint sa=0;sa<tplist.size();sa++) cerr << tplist.at(sa).Cart.at(0) << " " << tplist.at(sa).Cart.at(1) << " " << tplist.at(sa).Cart.at(2) << endl;
    }

    //**********
    //    cerr << tplist.size() << endl;
    if(tplist.size()==3) {tpsur=Surface(tplist,recnum,3);recnum++;SURFACE.push_back(tpsur);}
    else cerr << "ERROR: POINTS NOT COMPLETE!";
  }
  //  cerr << SURFACE.size() << endl;;
  for(uint i=0;i<TERS.size();i++) {
    int found_surface=-1;
    for(uint j=0;j<SURFACE.size();j++) {
      //   cerr << "SURFACE.size()=" << SURFACE.size() << endl;
      if(SURFACE.at(j).INSIDEpoint(TERS.at(i))) {
	TERS.at(i).DisttoPlane=SURFACE.at(j).CalculateDIST(TERS.at(i));
	//	cerr << "POINT: " << TERS.at(i).Phase.at(0) << " " << TERS.at(i).Phase.at(1) << " " << TERS.at(i).Phase.at(2) << " " << TERS.at(i).DisttoPlane << endl;
	//****DEBUG****
	if(0) {
	  cerr << "CORNERS:" << endl;
	  for(uint s=0;s<SURFACE.at(j).corners.size();s++) {
	    cerr << SURFACE.at(j).corners.at(s).Cart.at(0) << " " << SURFACE.at(j).corners.at(s).Cart.at(1) << " " << SURFACE.at(j).corners.at(s).Cart.at(2) << endl;
	  }
	  cerr << "POINT:" << endl;
	  cerr << TERS.at(i).Cart.at(0) << " " << TERS.at(i).Cart.at(1) << " " << TERS.at(i).Cart.at(2) << endl;
	  cerr << endl;
	}
	if(aurostd::abs(TERS.at(i).DisttoPlane)< 0.001) TERS.at(i).DisttoPlane=0.0;  // roundoff 1meV
	//cerr << SURFACE.at(j).normalv << endl;
	//cerr << "DIST: " << TERS.at(i).DisttoPlane << endl;
	found_surface=j;
	break;
      }
    }
    if(0) { // stefano
      cerr << "POINT (" << i << "/" << TERS.size() << "): " 
	   << TERS.at(i).Phase.at(0) << " " << TERS.at(i).Phase.at(1) << " " << TERS.at(i).Phase.at(2) << " "
	   << found_surface << " " << TERS.at(i).DisttoPlane << endl;
    }
  }

  //*********************************END THIS PART**************************************

  cerr << "00000  MESSAGE WARNING: COMPOUNDS WITH POSITIVE FORMATION ENTHALPIES ARE REMOVED FROM CONVEX HULL!" << endl;
  stringstream outfile;
  vector<int> rec;
  int NUM=1,s;
  for(uint i=0;i<POINTS.size();i++) rec.push_back(1);
  string filename="phasediagram_"+vspecies.at(0)+vspecies.at(1)+vspecies.at(2)+"."+aurostd::utype2string(TOTAL);
  string convex=aurostd::TmpFileCreate("3Dconvexgp");
  outfile << "set terminal postscript eps color enhanced size " << PIC_SIZE << " font " << aurostd::utype2string<int>(PIC_FONTSIZE) << endl;
  outfile << "unset border" << endl;
  outfile << "unset xtics" << endl;
  outfile << "unset ytics" << endl;
  outfile << "unset key" << endl; 
  outfile << "set xrange [-0.02:1.02]" << endl;
  outfile << "set yrange [-0.02:0.92]" << endl;
  if(TERDATA_COMPLETE) outfile << "set label " << NUM++ << " \"" << vspecies.at(0) << "\\nE_{F}= 0.0 eV\" at -0.015,-0.015" << endl;
  if(TERDATA_COMPLETE) outfile << "set label " << NUM++ << " \"" << vspecies.at(1) << "\\nE_{F}= 0.0 eV\" at 0.49,0.88" << endl;
  if(TERDATA_COMPLETE) outfile << "set label " << NUM++ << " \"" << vspecies.at(2) << "\\nE_{F}= 0.0 eV\" at 1.0,-0.015" << endl;
  if(!TERDATA_COMPLETE) outfile << "set label " << NUM++ << " \"" << vspecies.at(0) << "\\n\" at -0.015,-0.015" << endl;  // JESUS
  if(!TERDATA_COMPLETE) outfile << "set label " << NUM++ << " \"" << vspecies.at(1) << "\\n\" at 0.49,0.88" << endl; // JESUS
  if(!TERDATA_COMPLETE) outfile << "set label " << NUM++ << " \"" << vspecies.at(2) << "\\n\" at 1.0,-0.015" << endl;  // JESUS
  outfile << "set label " << NUM++ << " \"aflow.org\\n[" << TODAY << "]\\nAFLOW V" << VERSION << " calcs=" << TOTAL << "\\nstefano.curtarolo\" at 0.05,0.9" << endl;  // JESUS

  //*************************************
  bool check_Heusler=false;
  bool check_Half_Heusler=false;
  bool check_metastable_Heusler=false;
  bool check_metastable_Half_Heusler=false;

  bool tflag=TRUE;
  for(uint i=0;i<X.size();i++) {
    tflag=TRUE;
    if(aurostd::abs(X.at(i))<0.001) X.at(i)=0;
    if(aurostd::abs(Y.at(i))<0.001) Y.at(i)=0;
    if(checkfinished(FX,FY,X.at(i),Y.at(i))) tflag=FALSE;
    //cerr << X.at(i) << " " << Y.at(i) << " " << tflag << endl;
    for(uint j=0;j<POINTS.size()&&tflag;j++) {
      if(aurostd::abs(ENTHALPY.at(CX.at(i)).at(CY.at(i))-POINTS.at(j).Cart.at(2))<0.01) {
	if(aurostd::abs(X.at(i)-POINTS.at(j).Cart.at(0))<0.01) {
	  if(aurostd::abs(Y.at(i)-POINTS.at(j).Cart.at(1))<0.01) {
	    if(rec.at(j)==1) {
	      rec.at(j)=0;
	      FX.push_back(X.at(i));
	      FY.push_back(Y.at(i));
	      s=(int)j;
	    }
	  }
	}
      }
    }
    if(tflag&&s!=-1) {
      vector<string> Species;
      vector<double> Num;
      XATOM_SplitAlloySpecies(POINTS.at(s).CompoundName,Species,Num);
      //**********DEBUG***********
      if(0) {
	cerr << POINTS.at(s).CompoundName << endl;
	for(uint r=0;r<Species.size();r++) {
	  cerr << Species.at(r) << "  " << Num.at(r) << endl;
	}
	cerr << POINTS.at(s).Cart.at(0) << " " << POINTS.at(s).Cart.at(1) << " " << POINTS.at(s).Cart.at(2) << endl;
	cerr << X.at(i) << " " << Y.at(i) << " " << ENTHALPY.at(CX.at(i)).at(CY.at(i)) << endl;
	cerr << aurostd::abs(ENTHALPY.at(CX.at(i)).at(CY.at(i))-POINTS.at(s).Cart.at(2)) << endl;
	cerr << endl;
      }
      outfile << "set label " << NUM++ << " \"";
      NumReduction(Num);
      for(uint g=0;g<Species.size();g++) {
	outfile << Species.at(g);
	if(aurostd::abs(Num.at(g)-1.0)<epsilon) continue;
	else outfile << "_{" << Num.at(g) << "}";
      }
      if(TERDATA_COMPLETE) outfile << "\\nE_{F}= " << ENTHALPY.at(CX.at(i)).at(CY.at(i)) << " eV " << "\\nT_{s}= " << POINTS.at(s).Ts << "K\" at " << (X.at(i)-0.02) << "," << (Y.at(i)+0.02) << endl; 
      if(!TERDATA_COMPLETE) outfile << "\" at " << (X.at(i)-0.02) << "," << (Y.at(i)+0.02) << endl; // JESUS was Y+0.02
      if(Species.size()==3) {
	if(aurostd::isequal(Num.at(0),1.0) && aurostd::isequal(Num.at(1),1.0) && aurostd::isequal(Num.at(2),1.0)) check_Half_Heusler=TRUE;
	if(aurostd::isequal(Num.at(0),2.0) && aurostd::isequal(Num.at(1),1.0) && aurostd::isequal(Num.at(2),1.0)) check_Heusler=TRUE;
	if(aurostd::isequal(Num.at(0),1.0) && aurostd::isequal(Num.at(1),2.0) && aurostd::isequal(Num.at(2),1.0)) check_Heusler=TRUE;
	if(aurostd::isequal(Num.at(0),1.0) && aurostd::isequal(Num.at(1),1.0) && aurostd::isequal(Num.at(2),2.0)) check_Heusler=TRUE;
      }
      s=-1;
    }
  }

  outfile.setf(ios_base::fixed,ios_base::floatfield);
  map<string,int> existlists;
  for(unsigned int i=0;i<TERS.size();i++) {
    if(TERS.at(i).DisttoPlane>0) {
      if(existlists[TERS[i].CompoundName]==1) continue;
      else existlists[TERS[i].CompoundName]=1;
      vector<string> Species;
      vector<double> Num;
      outfile << "set label " << NUM++ << " \"";
      XATOM_SplitAlloySpecies(TERS.at(i).CompoundName,Species,Num);
      NumReduction(Num);
      for(uint g=0;g<Species.size();g++) {
	outfile << Species.at(g);
	if(aurostd::abs(Num.at(g)-1.0)<epsilon) continue;
	else outfile << "_{" << int(Num.at(g)) << "}";
      }
      if(TERDATA_COMPLETE) outfile << "\\nE_{F}= " << std::setprecision(4) << TERS.at(i).Cart.at(2) << " eV " << "\\nDis=" << std::setprecision(4) << TERS.at(i).DisttoPlane << " eV " << "\\nT_{s}= " << TERS.at(i).Ts << "K\" at " << (TERS.at(i).Cart.at(0)-0.02) << "," << (TERS.at(i).Cart.at(1)+0.02) << endl;
      if(!TERDATA_COMPLETE) outfile << "\" at " << (TERS.at(i).Cart.at(0)-0.02) << "," << (TERS.at(i).Cart.at(1)+0.02) << endl;
      
      if(Species.size()==3) {
	if(aurostd::isequal(Num.at(0),1.0) && aurostd::isequal(Num.at(1),1.0) && aurostd::isequal(Num.at(2),1.0)) check_metastable_Half_Heusler=TRUE;
	if(aurostd::isequal(Num.at(0),2.0) && aurostd::isequal(Num.at(1),1.0) && aurostd::isequal(Num.at(2),1.0)) check_metastable_Heusler=TRUE;
	if(aurostd::isequal(Num.at(0),1.0) && aurostd::isequal(Num.at(1),2.0) && aurostd::isequal(Num.at(2),1.0)) check_metastable_Heusler=TRUE;
	if(aurostd::isequal(Num.at(0),1.0) && aurostd::isequal(Num.at(1),1.0) && aurostd::isequal(Num.at(2),2.0)) check_metastable_Heusler=TRUE;
      }
    }
  }

  cerr << "check_Half_Heusler=" << check_Half_Heusler << endl;
  if(check_Half_Heusler) outfile << "set label " << NUM++ << " \"check=HalfHeusler\" at 0.05,0.76" << endl;  // JESUS
  cerr << "check_Heusler=" << check_Heusler << endl;
  if(check_Heusler) outfile << "set label " << NUM++ << " \"check=Heusler\" at 0.05,0.73" << endl;  // JESUS
  cerr << "check_metastable_Half_Heusler=" << check_metastable_Half_Heusler << endl;
  if(check_metastable_Half_Heusler) outfile << "set label " << NUM++ << " \"check=metaHalfHeusler\" at 0.05,0.70" << endl;  // JESUS
  cerr << "check_metastable_Heusler=" << check_metastable_Heusler << endl;
  if(check_metastable_Heusler) outfile << "set label " << NUM++ << " \"check=metaHeusler\" at 0.05,0.67" << endl;  // JESUS

  outfile << "set output '" << filename << ".eps" << endl;
  outfile << "set pointsize " << PIC_POINTSIZE << endl;
  outfile << "plot \'" << gnu_dat << "\' w lp lc -1";
  for(unsigned int i=0;i<TERS.size();i++) {
    if(TERS.at(i).DisttoPlane>0) {
      outfile << ", \\" << endl << "\"<echo \'" << TERS.at(i).Cart.at(0) << " " << TERS.at(i).Cart.at(1) << "\'\" w p lc 1 pt 2";
    }
  }
  outfile << endl;
  
  aurostd::stringstream2file(outfile,convex);
  ostringstream aus;
  aus << XHOST.command("gnuplot") << " " << convex << endl;
  aus << XHOST.command("epstopdf") << " " << filename << ".eps" << endl;
  if(XHOST.vflag_control.flag("PRINT_MODE::GIF")) aus << XHOST.command("convert") << " " << filename << ".eps" << " " << filename << ".gif" << endl;
  if(XHOST.vflag_control.flag("PRINT_MODE::JPG")) aus << XHOST.command("convert") << " " << filename << ".eps" << " " << filename << ".jpg" << endl;
  if(XHOST.vflag_control.flag("PRINT_MODE::PNG")) aus << XHOST.command("convert") << " " << filename << ".eps" << " " << filename << ".png" << endl;
  aurostd::execute(aus);

  // [OBSOLETE]  cerr << "XHOST.vflag_control.flag(\"KEEP::GPL\")=" << XHOST.vflag_control.flag("KEEP::GPL") << endl;
  // [OBSOLETE]  cerr << "XHOST.vflag_control.flag(\"KEEP::EPS\")=" << XHOST.vflag_control.flag("KEEP::EPS") << endl;
  if(!XHOST.vflag_control.flag("KEEP::GPL")) {
    aurostd::RemoveFile(convex);
  } else {
    aurostd::execute("mv -f "+convex+" "+filename+"_convex ");
  }
  if(!XHOST.vflag_control.flag("KEEP::GPL")) { 
    aurostd::RemoveFile(gnu_dat);
  } else {
    aurostd::execute("mv -f "+gnu_dat+" "+filename+"_gnu_dat");
  }
  if(!XHOST.vflag_control.flag("KEEP::EPS")) 
    aurostd::execute("rm -f "+filename+".eps");
}

//**********************************FINAL FUNCTION**********************
void GENERATESTABLELIST(vector<string>& argv) {
  vector<PhasePoints> POINTS,TERS;
  vector<vector<vector<double> > > LINE2D;
  vector<vector<double> > ENTHALPY;
  ostringstream oss;
  PhasePoints maxP;
  int TOTAL;
  ifstream infile;
  infile.open(argv.at(2).c_str());
  vector<string> tmp;
  vector<vector<string> > result;
  string line;
  while(!infile.eof()) {
    getline(infile,line);
    aurostd::string2tokens(line,tmp," ");
    result.push_back(tmp);
    tmp.clear();
  }
  infile.close();
  
  oss << aurostd::PaddedPRE("SYSTEM",10) << aurostd::PaddedPRE("COMPOUND",15) << aurostd::PaddedPRE("ENTHALPY(eV)",15) 
      << aurostd::PaddedPRE("TS(K)",10) << aurostd::PaddedPRE("PROTO",15) << aurostd::PaddedPRE("MAX_TS",30) << aurostd::PaddedPRE("DIRS",70) << endl;
  for(uint i=0;i<result.size();i++) {
    oss << aurostd::PaddedPRE(result[i][2]+result[i][3]+result[i][4],10);
    INPUTDATAFORTERPHASE_NEW(result[i],POINTS,TOTAL,TERS,maxP);
    TERNARY_CONVEX(POINTS,LINE2D,ENTHALPY);
    GetListName(LINE2D,ENTHALPY,POINTS,oss,maxP);
    POINTS.clear();TERS.clear();LINE2D.clear();ENTHALPY.clear();
  }
  ofstream outfile;
  outfile.open("exist_lists.txt");
  outfile << oss.str();
  outfile.close();
}

void INPUTDATAFORTERPHASE(vector<string>& argv) {
  // bool LDEBUG=(FALSE || XHOST.DEBUG);
  vector<PhasePoints> POINTS,TERS;
  vector<vector<vector<double> > > LINE2D;
  vector<vector<double> > ENTHALPY;
  int TOTAL;

  vector<string> vspecies;
  if(!aurostd::getproto_itemized_vector_string_from_input(argv,"--terdata",vspecies)) {}; 
  aurostd::sort(vspecies);
  if(vspecies.size()!=3) {
    cerr << "ERROR: INPUTDATAFORTERPHASE is made for ternaries..... exiting.." << endl; exit(0);
  }

  int PIC_FONTSIZE=aurostd::args2attachedutype<int>(argv,"--font=|--fonts=",PIC_FONTSIZE_DEFAULT);

  INPUTDATAFORTERPHASE_NEW(vspecies,POINTS,TOTAL,TERS);
  TERNARY_CONVEX(POINTS,LINE2D,ENTHALPY);
  DigramPlot3ary(LINE2D,ENTHALPY,vspecies,POINTS,TOTAL,TERS,PIC_FONTSIZE);  

}

#endif //  _AFLOW_APENNSY_3D_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
