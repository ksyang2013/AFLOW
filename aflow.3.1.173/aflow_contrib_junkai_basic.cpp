// ***************************************************************************
// *                                                                         *
// *           Aflow JUNKAI XUE - Duke University 2003-2018                  *
// *                                                                         *
// ***************************************************************************
// Written by Junkai: 2010-2011
// junkai.xue@duke.edu

#include<iostream>
#include<cmath>
#include "aflow.h"
#include<vector>

using aurostd::det;
using aurostd::execute;



void SORT(vector<int>& a) {
  int tmp;
  for(uint i=0;i<a.size();i++) {
    for(uint j=i+1;j<a.size();j++) {
      if(a.at(i)<a.at(j)) {
        tmp=a.at(i);
        a.at(i)=a.at(j);
        a.at(j)=tmp;
      }
      else continue;
    }
  }
}


void SORT(vector<double>& a) {
  double tmp;
  for(uint i=0;i<a.size();i++) {
    for(uint j=i+1;j<a.size();j++) {
      if(a.at(i)>a.at(j)) {
        tmp=a.at(i);
        a.at(i)=a.at(j);
        a.at(j)=tmp;
      }
      else continue;
    }
  }
}



void  SORTDecrease(vector<uint>& a) {
  uint tmp;
  for(uint i=0;i<a.size();i++) {
    for(uint j=i+1;j<a.size();j++) {
      if(a.at(i)>a.at(j)) {
        tmp=a.at(i);
        a.at(i)=a.at(j);
        a.at(j)=tmp;
      }
      else continue;
    }
  }
}



double NORM(const xvector<double>& v) {
  //if(v[1]==0&&v[2]==0&v[3]==0) cerr<<"0 comes!"<<endl;                                                                                                                                                            
   return sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
}



double COS(const xvector<double>& v1, const xvector<double>& v2) {
  double dotpro,normpro;
  dotpro=v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
  normpro=NORM(v1)*NORM(v2);
  if(normpro==0) {cerr<<"ERROR:Norm product is Zero!"<<endl;exit(0);}
  else return(dotpro/normpro);
}



double DotPro(const xvector<double>& v1, const xvector<double>& v2) {
  double dotpro;
  dotpro=v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
  return dotpro;
}



void CrossPro(const xvector<double>& a, const xvector<double>& b, xvector<double>& resu) {
  resu[1]=a[2]*b[3]-b[2]*a[3];
  resu[2]=a[3]*b[1]-b[3]*a[1];
  resu[3]=a[1]*b[2]-b[1]*a[2];
}
void POSCAR_OutputFun(xstructure& a) {
  int r=0;
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.precision(10);
  cout<<a.title<<endl;
  cout<<-a.Volume()<<endl;
  //for(int h=1;h<4;h++) cout<<a.lattice(h)<<endl;
  for(int h=1;h<4;h++) {
for(int hh=1; hh<4; hh++) {
cout<<a.lattice(h,hh)<<" ";
}
cout << endl;
}
  for(int g=0;g<(int)a.num_each_type.size();g++) cout<<a.num_each_type.at(g)<<" ";
  cout<<endl<<"Direct("<<a.atoms.size()<<")"<<endl;
  //cerr<<a.num_each_type.size()<<endl;
  for(int g=0;g<(int)a.num_each_type.size();g++) {
    //    cerr<<a.num_each_type.at(g)<<endl;
    for(int h=0;h<(int)a.num_each_type.at(g);h++) {
      if(a.atoms.at(r).name!="") {
for(int hh=1; hh<4; hh++) {
cout<< a.atoms.at(r).fpos[hh] << " ";
}
cout <<" "<<a.atoms.at(r).name<<endl;
//cout<<a.atoms.at(r).fpos<<" "<<a.atoms.at(r).name<<endl;
}
      else {
for(int hh=1; hh<4; hh++) {
cout<< a.atoms.at(r).fpos[hh] << " ";
}
cout <<" "<<char('A'+g)<<endl;
}
      r++;
    }
  }
  cout<<endl;
}

void POSCAR_Reread(xstructure& a) {
  a.Standard_Primitive_UnitCellForm();
  stringstream oss;
  int r=0;
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  oss<<a.title<<endl;
  oss<<-a.Volume()<<endl;
  //for(int h=1;h<4;h++) oss<<a.lattice(h)<<endl;                                                                                                                                                                  
  for(int h=1;h<4;h++) {
    for(int hh=1; hh<4; hh++) {
      oss<<a.lattice(h,hh)<<" ";
    }
    oss << endl;
  }
  for(int g=0;g<(int)a.num_each_type.size();g++) oss<<a.num_each_type.at(g)<<" ";
  oss<<endl<<"Direct("<<a.atoms.size()<<")"<<endl;
  //cerr<<a.num_each_type.size()<<endl;                                                                                                                                                                              
  for(int g=0;g<(int)a.num_each_type.size();g++) {
    //    cerr<<a.num_each_type.at(g)<<endl;                                                                                                                                                                        
    for(int h=0;h<(int)a.num_each_type.at(g);h++) {
      if(a.atoms.at(r).name!="") {
	for(int hh=1; hh<4; hh++) {
	  oss<< a.atoms.at(r).fpos[hh] << " ";
	}
	oss <<" "<<a.atoms.at(r).name<<endl;
	//oss<<a.atoms.at(r).fpos<<" "<<a.atoms.at(r).name<<endl;                                                                                                                                                          
      }
      else {
	for(int hh=1; hh<4; hh++) {
	  oss<< a.atoms.at(r).fpos[hh] << " ";
	}
	oss <<" "<<char('A'+g)<<endl;
      }
      r++;
    }
  }
  oss<<endl;
  xstructure b(oss,IOAFLOW_AUTO);
  a=b;
}


string POSCAR_Output(xstructure& a) {
  int r=0;
  ostringstream OUT;
  OUT.setf(std::ios::fixed,std::ios::floatfield);
  OUT.precision(10);
  OUT<<a.title<<endl;
  OUT<<-a.Volume()<<endl;
  for(int h=1;h<4;h++) {
    for(int hh=1; hh<4; hh++) {
      OUT<<a.lattice(h,hh)<<" ";
    }
    OUT << endl;
  }
  for(int g=0;g<(int)a.num_each_type.size();g++) OUT<<a.num_each_type.at(g)<<" ";
  OUT<<endl<<"Direct("<<a.atoms.size()<<")"<<endl;
  for(int g=0;g<(int)a.num_each_type.size();g++) {
    for(int h=0;h<(int)a.num_each_type.at(g);h++) {
      if(a.atoms.at(r).name!="") {
        for(int hh=1; hh<4; hh++) {
          OUT<< a.atoms.at(r).fpos[hh] << " ";
        }
        OUT <<" "<<a.atoms.at(r).name<<endl;
      }
      else {
        for(int hh=1; hh<4; hh++) {
          OUT<< a.atoms.at(r).fpos[hh] << " ";
        }
        OUT <<" "<<char('A'+g)<<endl;
      }
      r++;
    }
  }
  OUT<<endl;
  return OUT.str();
}


bool LinearIndependence(const xmatrix<double>& a) {
  if(abs(det(a))<0.1) return false;
  return true;
}

bool VectSort(const xvector<double>& a,const  xvector<double>& b) {
  return(NORM(a)<NORM(b));
}

bool AtomSort(const _atom& a, const _atom& b) {
  return (NORM(a.fpos)<NORM(b.fpos));
}

void SpeciesSeparate(string& a, vector<string>& spe, vector<double>& nv) {
  string tpname="",tpnum="",rec;
  bool flag=false,pushflag=false;
  for(uint i=0;i<a.length();i++) {
    if(!isdigit(a[i])) {
      if(flag) {
	rec=tpname;
	tpname="";
	pushflag=true;
	flag=false;
      }
      tpname+=a[i];
    }
    else{
      tpnum+=a[i];
      flag=true;
    }
      
    if(pushflag) {
      spe.push_back(rec);
      nv.push_back(aurostd::string2utype<double>(tpnum));
      tpnum="";
      pushflag=false;
    }
    //    cerr<<i<<" "<<a.length()<<endl;
    if(i==a.length()-1) {
      // cerr<<"HERE??"<<endl;
      spe.push_back(tpname);
      nv.push_back(aurostd::string2utype<double>(tpnum));
    }
  }
  
}

bool DsortF(double i, double j) {return (i>j);}
