// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2008           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011
// Positions taken from the Bilbao Crystallographic Database

#ifndef _WYCKOFF_CPP_
#define _WYCKOFF_CPP_

#include "aflow.h"

// all vector./matrix must be loaded before

void SpaceGroupOptionRequired(uint &spacegroup, uint &option) {
  cerr << "AFLOW WARNING: Wyckoff Spacegroup " << spacegroup << " requires option 1 or 2 (" << option << ")" << endl;
  // cout << "AFLOW WARNING: Wyckoff Spacegroup " << spacegroup << " requires option 1 or 2 (" << option << ")" << endl;
  //  cerr << "AFLOW WARNING: check it out and double-check concentrations and space-group" << endl;
  // cerr << "EXITING..." << endl;
  //  exit(0);SpaceGroupOptionRequired(spacegroup,option);
  if(option==0 || option>3) {
    option=1;
    cerr << "AFLOW WARNING: Wyckoff Spacegroup " << spacegroup << " taking option=" << option << " (let`s hope it is the right one, check the concentrations and space-group)" << endl;
    //   cout << "AFLOW WARNING: Taking option=" << option << " (let`s hope it is the right one, check the concentrations and space-group)" << endl;
  }
  if(option==3) {
    option=2;
    cerr << "AFLOW WARNING: Wyckoff Spacegroup " << spacegroup << " taking option=" << option << " (let`s hope it is the right one, check the concentrations and space-group)" << endl;
    //  cout << "AFLOW WARNING: Taking option=" << option << " (let`s hope it is the right one, check the concentrations and space-group)" << endl;
  }
}

bool SpaceGroupOptionRequired(uint sg) {
  if(sg==3   || sg==4   || sg==5   || sg==6   || sg==7   ||   sg==8) return TRUE;
  if(sg==9   || sg==10  || sg==11  || sg==12  || sg==13  ||  sg==14) return TRUE;
  if(sg==15  || sg==48  || sg==50  || sg==59  || sg==68  ||  sg==70) return TRUE;
  if(sg==85  || sg==86  || sg==88  || sg==125 || sg==126 || sg==129) return TRUE;
  if(sg==130 || sg==133 || sg==134 || sg==137 || sg==138 || sg==141) return TRUE;
  if(sg==142 || sg==146 || sg==148 || sg==155 || sg==160 || sg==161) return TRUE;
  if(sg==166 || sg==167 || sg==201 || sg==203 || sg==222 || sg==224) return TRUE;
  if(sg==227 || sg==228) return TRUE;
  return FALSE;
}


// ----------------------------------------------------------------------------
xvector<double> wv(const double &x,const double &y,const double &z) {
  xvector<double> out(3);
  out[1]=x;out[2]=y;out[3]=z;
  return out;
}

// ----------------------------------------------------------------------------
void wa(_atom& a,xstructure &str) {
  a=BringInCell(a,str.lattice);
  a.cpos=F2C(str.lattice,a.fpos);
  str.AddAtom(a);  
}

// ----------------------------------------------------------------------------
xstructure WyckoffPOSITIONS(uint spacegroup, xstructure strin) {
  return  WyckoffPOSITIONS(spacegroup,(uint)1, strin);
}

xstructure WyckoffPOSITIONS(uint spacegroup_in, uint option_in, xstructure strin) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  
  if(LDEBUG) cerr << "WyckoffPOSITIONS [0]" << endl;
  
  xvector<double> o(3);
  xstructure str(strin);
  double x=0,y=0,z=0;
  _atom a;
  uint spacegroup=spacegroup_in;
  uint option=option_in;

  while (str.atoms.size()>0) { str.RemoveAtom(0); } // strip all atoms
  // str.species.clear();str.species_pp.clear();str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();species_pp_vLDAU.clear();str.species_volume.clear();str.species_mass.clear(); // patch for RemoveAtom

  if(LDEBUG) cerr << "WyckoffPOSITIONS [1]" << endl;
  
  for(uint i=0;i<strin.atoms.size();i++) {
    x=strin.atoms.at(i).fpos[1];
    y=strin.atoms.at(i).fpos[2];
    z=strin.atoms.at(i).fpos[3];
    a=strin.atoms.at(i);

//     //(0,0,0)(1./2,1./2,0)
//     for(uint j=1;j<=2;j++) {
//       if(j==1) o=wv(0,0,0);
//       if(j==2) o=wv(1./2,1./2,0);
//     }
//     //(0,0,0)(0,1./2,1./2)
//     for(uint j=1;j<=2;j++) {
//       if(j==1) o=wv(0,0,0);
//       if(j==2) o=wv(0,1./2,1./2);
//     }
//     //(0,0,0)(1./2,1./2,1./2)
//     for(uint j=1;j<=2;j++) {
//       if(j==1) o=wv(0,0,0);
//       if(j==2) o=wv(1./2,1./2,1./2);
//     }
//     //(0,0,0)(2./3,1./3,1./3)(1./3,2./3,2./3)
//     for(uint j=1;j<=3;j++) {
//       if(j==1) o=wv(0,0,0);
//       if(j==2) o=wv(2./3,1./3,1./3);
//       if(j==3) o=wv(1./3,2./3,2./3);
//     }
//     //(0,0,0)(0,1./2,1./2)(1./2,0,1./2)(1./2,1./2,0)
//     //(0,0,0)(1./2,0,1./2)(0,1./2,1./2)(1./2,1./2,0)
//     for(uint j=1;j<=4;j++) {
//       if(j==1) o=wv(0,0,0);
//       if(j==2) o=wv(0,1./2,1./2);
//       if(j==3) o=wv(1./2,0,1./2);
//       if(j==4) o=wv(1./2,1./2,0);
//     }
    
    // for all spacegroups
    o=wv(0,0,0);

    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    if(option!=1 && option!=2) {
      if(spacegroup::SpaceGroupOptionRequired(spacegroup)) SpaceGroupOptionRequired(spacegroup,option);
      //    cerr << spacegroup << " " << spacegroup::SpaceGroupOptionRequired(spacegroup) << endl;
     // if(SpaceGroupOptionRequired(spacegroup)) SpaceGroupOptionRequired(spacegroup,option);
    }
    if(spacegroup::SpaceGroupOptionRequired(spacegroup)) str.spacegroupnumberoption=option;    // plug option inside   
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 1  P1 #1
    if(spacegroup==1) {
      str.spacegroup="P1";
      str.spacegrouplabel="#1";
      str.spacegroupnumber=1;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 2 P-1 #2
    if(spacegroup==2) {
      str.spacegroup="P-1";
      str.spacegrouplabel="#2";
      str.spacegroupnumber=2;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 3  P2 #3
    if(spacegroup==3 && option==1) {
      str.spacegroup="P2";
      str.spacegrouplabel="#3";
      str.spacegroupnumber=3;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
    }
    if(spacegroup==3 && option==2) {
      str.spacegroup="P2";
      str.spacegrouplabel="#3";
      str.spacegroupnumber=3;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 4  P2_{1} #4
    if(spacegroup==4 && option==1) {
      str.spacegroup="P2_{1}";
      str.spacegrouplabel="#4";
      str.spacegroupnumber=4;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
    }
    if(spacegroup==4 && option==2) {
      str.spacegroup="P2_{1}";
      str.spacegrouplabel="#4";
      str.spacegroupnumber=4;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 5  C2 #5
    if(spacegroup==5 && option==1) {
      str.spacegroup="C2";
      str.spacegrouplabel="#5";
      str.spacegroupnumber=5;
      str.spacegroupoption="unique axis b";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
      }
    }
    if(spacegroup==5 && option==2) {
      str.spacegroup="C2";
      str.spacegrouplabel="#5";
      str.spacegroupnumber=5;
      str.spacegroupoption="unique axis c";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 6  Pm #6
    if(spacegroup==6 && option==1) {
      str.spacegroup="Pm";
      str.spacegrouplabel="#6";
      str.spacegroupnumber=6;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
    }
    if(spacegroup==6 && option==2) {
      str.spacegroup="Pm";
      str.spacegrouplabel="#6";
      str.spacegroupnumber=6;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 7  Pc #7
    if(spacegroup==7 && option==1) {
      str.spacegroup="Pc";
      str.spacegrouplabel="#7";
      str.spacegroupnumber=7;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
    }
    if(spacegroup==7 && option==2) {
      str.spacegroup="Pc";
      str.spacegrouplabel="#7";
      str.spacegroupnumber=7;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 8  Cm #8
    if(spacegroup==8 && option==1) {
      str.spacegroup="Cm";
      str.spacegrouplabel="#8";
      str.spacegroupnumber=8;
      str.spacegroupoption="unique axis b";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
      }
    }
    if(spacegroup==8 && option==2) {
      str.spacegroup="Cm";
      str.spacegrouplabel="#8";
      str.spacegroupnumber=8;
      str.spacegroupoption="unique axis c";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 9  Cc #9
    if(spacegroup==9 && option==1) {
      str.spacegroup="Cc";
      str.spacegrouplabel="#9";
      str.spacegroupnumber=9;
      str.spacegroupoption="unique axis b";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      }
    }
    if(spacegroup==9 && option==2) {
      str.spacegroup="Cc";
      str.spacegrouplabel="#9";
      str.spacegroupnumber=9;
      str.spacegroupoption="unique axis c";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 10  P2./m #10
    if(spacegroup==10 && option==1) {
      str.spacegroup="P2./m";
      str.spacegrouplabel="#10";
      str.spacegroupnumber=10;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
    }
    if(spacegroup==10 && option==2) {
      str.spacegroup="P2./m";
      str.spacegrouplabel="#10";
      str.spacegroupnumber=10;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 11  P2_{1}./m #11
    if(spacegroup==11 && option==1) {
      str.spacegroup="P2_{1}./m";
      str.spacegrouplabel="#11";
      str.spacegroupnumber=11;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
    }
    if(spacegroup==11 && option==2) {
      str.spacegroup="P2_{1}./m";
      str.spacegrouplabel="#11";
      str.spacegroupnumber=11;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 12  C2./m #12
    if(spacegroup==12 && option==1) {
      str.spacegroup="C2./m";
      str.spacegrouplabel="#12";
      str.spacegroupnumber=12;
      str.spacegroupoption="unique axis b";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
      }
    }
    if(spacegroup==12 && option==2) {
      str.spacegroup="C2./m";
      str.spacegrouplabel="#12";
      str.spacegroupnumber=12;
      str.spacegroupoption="unique axis c";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 13  P2./c #13
    if(spacegroup==13 && option==1) {
      str.spacegroup="P2./c";
      str.spacegrouplabel="#13";
      str.spacegroupnumber=13;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
    }
    if(spacegroup==13 && option==2) {
      str.spacegroup="P2./c";
      str.spacegrouplabel="#13";
      str.spacegroupnumber=13;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 14  P2_{1}./c #14
    if(spacegroup==14 && option==1) {
      str.spacegroup="P2_{1}./c";
      str.spacegrouplabel="#14";
      str.spacegroupnumber=14;
      str.spacegroupoption="unique axis b";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
    }
    if(spacegroup==14 && option==2) {
      str.spacegroup="P2_{1}./c";
      str.spacegrouplabel="#14";
      str.spacegroupnumber=14;
      str.spacegroupoption="unique axis c";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 15  C2./c #15
    if(spacegroup==15 && option==1) {
      str.spacegroup="C2./c";
      str.spacegrouplabel="#15";
      str.spacegroupnumber=15;
      str.spacegroupoption="unique axis b";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      }
    }
    if(spacegroup==15 && option==2) {
      str.spacegroup="C2./c";
      str.spacegrouplabel="#15";
      str.spacegroupnumber=15;
      str.spacegroupoption="unique axis c";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 16  P222 #16
    if(spacegroup==16) {
      str.spacegroup="P222";
      str.spacegrouplabel="#16";
      str.spacegroupnumber=16;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 17  P222_{1} #17
    if(spacegroup==17) {
      str.spacegroup="P222_{1}";
      str.spacegrouplabel="#17";
      str.spacegroupnumber=17;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 18  P2_{1}2_{1}2 #18
    if(spacegroup==18) {
      str.spacegroup="P2_{1}2_{1}2";
      str.spacegrouplabel="#18";
      str.spacegroupnumber=18;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 19  P2_{1}2_{1}2_{1} #19
    if(spacegroup==19) {
      str.spacegroup="P2_{1}2_{1}2_{1}";
      str.spacegrouplabel="#19";
      str.spacegroupnumber=19;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 20  C222_{1} #20
    if(spacegroup==20) {
      str.spacegroup="C222_{1}";
      str.spacegrouplabel="#20";
      str.spacegroupnumber=20;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 21  C222 #21
    if(spacegroup==21) {
      str.spacegroup="C222";
      str.spacegrouplabel="#21";
      str.spacegroupnumber=21;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 22  F222 #22
    if(spacegroup==22) {
      str.spacegroup="F222";
      str.spacegrouplabel="#22";
      str.spacegroupnumber=22;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 23  I222 #23
    if(spacegroup==23) {
      str.spacegroup="I222";
      str.spacegrouplabel="#23";
      str.spacegroupnumber=23;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 24  I2_{1}2_{1}2_{1} #24
    if(spacegroup==24) {
      str.spacegroup="I2_{1}2_{1}2_{1}";
      str.spacegrouplabel="#24";
      str.spacegroupnumber=24;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 25  Pmm2 #25
    if(spacegroup==25) {
      str.spacegroup="Pmm2";
      str.spacegrouplabel="#25";
      str.spacegroupnumber=25;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 26  Pmc2_{1} #26
    if(spacegroup==26) {
      str.spacegroup="Pmc2_{1}";
      str.spacegrouplabel="#26";
      str.spacegroupnumber=26;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 27  Pcc2 #27
    if(spacegroup==27) {
      str.spacegroup="Pcc2";
      str.spacegrouplabel="#27";
      str.spacegroupnumber=27;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 28  Pma2 #28
    if(spacegroup==28) {
      str.spacegroup="Pma2";
      str.spacegrouplabel="#28";
      str.spacegroupnumber=28;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 29  Pca2_{1} #29
    if(spacegroup==29) {
      str.spacegroup="Pca2_{1}";
      str.spacegrouplabel="#29";
      str.spacegroupnumber=29;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 30  Pnc2 #30
    if(spacegroup==30) {
      str.spacegroup="Pnc2";
      str.spacegrouplabel="#30";
      str.spacegroupnumber=30;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 31  Pmn2_{1} #31
    if(spacegroup==31) {
      str.spacegroup="Pmn2_{1}";
      str.spacegrouplabel="#31";
      str.spacegroupnumber=31;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 32  Pba2 #32
    if(spacegroup==32) {
      str.spacegroup="Pba2";
      str.spacegrouplabel="#32";
      str.spacegroupnumber=32;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 33  Pna2_{1} #33
    if(spacegroup==33) {
      str.spacegroup="Pna2_{1}";
      str.spacegrouplabel="#33";
      str.spacegroupnumber=33;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 34  Pnn2 #34
    if(spacegroup==34) {
      str.spacegroup="Pnn2";
      str.spacegrouplabel="#34";
      str.spacegroupnumber=34;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 35  Cmm2 #35
    if(spacegroup==35) {
      str.spacegroup="Cmm2";
      str.spacegrouplabel="#35";
      str.spacegroupnumber=35;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 36  Cmc2_{1} #36
    if(spacegroup==36) {
      str.spacegroup="Cmc2_{1}";
      str.spacegrouplabel="#36";
      str.spacegroupnumber=36;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 37  Ccc2 #37
    if(spacegroup==37) {
      str.spacegroup="Ccc2";
      str.spacegrouplabel="#37";
      str.spacegroupnumber=37;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 38  Amm2 #38
    if(spacegroup==38) {
      str.spacegroup="Amm2";
      str.spacegrouplabel="#38";
      str.spacegroupnumber=38;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 39  Aem2 #39
    if(spacegroup==39) {
      str.spacegroup="Aem2";
      str.spacegrouplabel="#39";
      str.spacegroupnumber=39;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 40  Ama2 #40
    if(spacegroup==40) {
      str.spacegroup="Ama2";
      str.spacegrouplabel="#40";
      str.spacegroupnumber=40;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 41  Aea2 #41
    if(spacegroup==41) {
      str.spacegroup="Aea2";
      str.spacegrouplabel="#41";
      str.spacegroupnumber=41;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 42  Fmm2 #42
    if(spacegroup==42) {
      str.spacegroup="Fmm2";
      str.spacegrouplabel="#42";
      str.spacegroupnumber=42;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 43  Fdd2 #43
    if(spacegroup==43) {
      str.spacegroup="Fdd2";
      str.spacegrouplabel="#43";
      str.spacegroupnumber=43;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x+1./4,-y+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,y+1./4,z+1./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 44  Imm2 #44
    if(spacegroup==44) {
      str.spacegroup="Imm2";
      str.spacegrouplabel="#44";
      str.spacegroupnumber=44;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 45  Iba2 #45
    if(spacegroup==45) {
      str.spacegroup="Iba2";
      str.spacegrouplabel="#45";
      str.spacegroupnumber=45;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 46  Ima2 #46
    if(spacegroup==46) {
      str.spacegroup="Ima2";
      str.spacegrouplabel="#46";
      str.spacegroupnumber=46;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 47  Pmmm #47
    if(spacegroup==47) {
      str.spacegroup="Pmmm";
      str.spacegrouplabel="#47";
      str.spacegroupnumber=47;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 48  Pnnn #48
    if(spacegroup==48 && option==1) {
      str.spacegroup="Pnnn";
      str.spacegrouplabel="#48";
      str.spacegroupnumber=48;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
    }
    if(spacegroup==48 && option==2) {
      str.spacegroup="Pnnn";
      str.spacegrouplabel="#48";
      str.spacegroupnumber=48;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 49  Pccm #49
    if(spacegroup==49) {
      str.spacegroup="Pccm";
      str.spacegrouplabel="#49";
      str.spacegroupnumber=49;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 50  Pban #50
    if(spacegroup==50 && option==1) {
      str.spacegroup="Pban";
      str.spacegrouplabel="#50";
      str.spacegroupnumber=50;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
    }
    if(spacegroup==50 && option==2) {
      str.spacegroup="Pban";
      str.spacegrouplabel="#50";
      str.spacegroupnumber=50;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 51  Pmma #51
    if(spacegroup==51) {
      str.spacegroup="Pmma";
      str.spacegrouplabel="#51";
      str.spacegroupnumber=51;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 52  Pnna #52
    if(spacegroup==52) {
      str.spacegroup="Pnna";
      str.spacegrouplabel="#52";
      str.spacegroupnumber=52;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 53  Pmna #53
    if(spacegroup==53) {
      str.spacegroup="Pmna";
      str.spacegrouplabel="#53";
      str.spacegroupnumber=53;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 54  Pcca #54
    if(spacegroup==54) {
      str.spacegroup="Pcca";
      str.spacegrouplabel="#54";
      str.spacegroupnumber=54;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 55  Pbam #55
    if(spacegroup==55) {
      str.spacegroup="Pbam";
      str.spacegrouplabel="#55";
      str.spacegroupnumber=55;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 56  Pccn #56
    if(spacegroup==56) {
      str.spacegroup="Pccn";
      str.spacegrouplabel="#56";
      str.spacegroupnumber=56;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 57  Pbcm #57
    if(spacegroup==57) {
      str.spacegroup="Pbcm";
      str.spacegrouplabel="#57";
      str.spacegroupnumber=57;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 58  Pnnm #58
    if(spacegroup==58) {
      str.spacegroup="Pnnm";
      str.spacegrouplabel="#58";
      str.spacegroupnumber=58;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 59  Pmmn #59
    if(spacegroup==59 && option==1) {
      str.spacegroup="Pmmn";
      str.spacegrouplabel="#59";
      str.spacegroupnumber=59;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
    }
    if(spacegroup==59 && option==2) {
      str.spacegroup="Pmmn";
      str.spacegrouplabel="#59";
      str.spacegroupnumber=59;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 60  Pbcn #60
    if(spacegroup==60) {
      str.spacegroup="Pbcn";
      str.spacegrouplabel="#60";
      str.spacegroupnumber=60;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 61  Pbca #61
    if(spacegroup==61) {
      str.spacegroup="Pbca";
      str.spacegrouplabel="#61";
      str.spacegroupnumber=61;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 62  Pnma #62
    if(spacegroup==62) {
      str.spacegroup="Pnma";
      str.spacegrouplabel="#62";
      str.spacegroupnumber=62;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 63  Cmcm #63
    if(spacegroup==63) {
      str.spacegroup="Cmcm";
      str.spacegrouplabel="#63";
      str.spacegroupnumber=63;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 64  Cmce #64
    if(spacegroup==64) {
      str.spacegroup="Cmce";
      str.spacegrouplabel="#64";
      str.spacegroupnumber=64;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 65  Cmmm #65
    if(spacegroup==65) {
      str.spacegroup="Cmmm";
      str.spacegrouplabel="#65";
      str.spacegroupnumber=65;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 66  Cccm #66
    if(spacegroup==66) {
      str.spacegroup="Cccm";
      str.spacegrouplabel="#66";
      str.spacegroupnumber=66;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 67  Cmme #67
    if(spacegroup==67) {
      str.spacegroup="Cmme";
      str.spacegrouplabel="#67";
      str.spacegroupnumber=67;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 68  Ccce #68
    if(spacegroup==68 && option==1) {
      str.spacegroup="Ccce";
      str.spacegrouplabel="#68";
      str.spacegroupnumber=68;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
      }
    }
    if(spacegroup==68 && option==2) {
      str.spacegroup="Ccce";
      str.spacegrouplabel="#68";
      str.spacegroupnumber=68;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 69  Fmmm #69
    if(spacegroup==69) {
      str.spacegroup="Fmmm";
      str.spacegrouplabel="#69";
      str.spacegroupnumber=69;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 70  Fddd #70
    if(spacegroup==70 && option==1) {
      str.spacegroup="Fddd";
      str.spacegrouplabel="#70";
      str.spacegroupnumber=70;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x+1./4,-y+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,y+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-y+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,y+1./4,z+1./4);wa(a,str);
      }
    }
    if(spacegroup==70 && option==2) {
      str.spacegroup="Fddd";
      str.spacegrouplabel="#70";
      str.spacegroupnumber=70;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+3./4,-y+3./4,z);wa(a,str);
	a.fpos=o+wv(-x+3./4,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(x,-y+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./4,y+1./4,-z);wa(a,str);
	a.fpos=o+wv(x+1./4,-y,z+1./4);wa(a,str);
	a.fpos=o+wv(-x,y+1./4,z+1./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 71  Immm #71
    if(spacegroup==71) {
      str.spacegroup="Immm";
      str.spacegrouplabel="#71";
      str.spacegroupnumber=71;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 72  Ibam #72
    if(spacegroup==72) {
      str.spacegroup="Ibam";
      str.spacegrouplabel="#72";
      str.spacegroupnumber=72;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 73  Ibca #73
    if(spacegroup==73) {
      str.spacegroup="Ibca";
      str.spacegrouplabel="#73";
      str.spacegroupnumber=73;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 74  Imma #74
    if(spacegroup==74) {
      str.spacegroup="Imma";
      str.spacegrouplabel="#74";
      str.spacegroupnumber=74;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 75  P4 #75
    if(spacegroup==75) {
      str.spacegroup="P4";
      str.spacegrouplabel="#75";
      str.spacegroupnumber=75;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 76  P4_{1} #76
    if(spacegroup==76) {
      str.spacegroup="P4_{1}";
      str.spacegrouplabel="#76";
      str.spacegroupnumber=76;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./4);wa(a,str);
      a.fpos=o+wv(y,-x,z+3./4);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 77  P4_{2} #77
    if(spacegroup==77) {
      str.spacegroup="P4_{2}";
      str.spacegrouplabel="#77";
      str.spacegroupnumber=77;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 78  P4_{3} #78
    if(spacegroup==78) {
      str.spacegroup="P4_{3}";
      str.spacegrouplabel="#78";
      str.spacegroupnumber=78;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,z+3./4);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./4);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 79  I4 #79
    if(spacegroup==79) {
      str.spacegroup="I4";
      str.spacegrouplabel="#79";
      str.spacegroupnumber=79;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 80  I4_{1} #80
    if(spacegroup==80) {
      str.spacegroup="I4_{1}";
      str.spacegrouplabel="#80";
      str.spacegroupnumber=80;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 81  P-4 #81
    if(spacegroup==81) {
      str.spacegroup="P-4";
      str.spacegrouplabel="#81";
      str.spacegroupnumber=81;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 82  I-4 #82
    if(spacegroup==82) {
      str.spacegroup="I-4";
      str.spacegrouplabel="#82";
      str.spacegroupnumber=82;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 83  P4./m #83
    if(spacegroup==83) {
      str.spacegroup="P4./m";
      str.spacegrouplabel="#83";
      str.spacegroupnumber=83;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 84  P4_{2}./m #84
    if(spacegroup==84) {
      str.spacegroup="P4_{2}./m";
      str.spacegrouplabel="#84";
      str.spacegroupnumber=84;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 85  P4./n #85
    if(spacegroup==85 && option==1) {
      str.spacegroup="P4./n";
      str.spacegrouplabel="#85";
      str.spacegroupnumber=85;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
    }
    if(spacegroup==85 && option==2) {
      str.spacegroup="P4./n";
      str.spacegrouplabel="#85";
      str.spacegroupnumber=85;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 86  P4_{2}./n #86
    if(spacegroup==86 && option==1) {
      str.spacegroup="P4_{2}./n";
      str.spacegrouplabel="#86";
      str.spacegroupnumber=86;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
    }
    if(spacegroup==86 && option==2) {
      str.spacegroup="P4_{2}./n";
      str.spacegrouplabel="#86";
      str.spacegroupnumber=86;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 87  I4./m #87
    if(spacegroup==87) {
      str.spacegroup="I4./m";
      str.spacegrouplabel="#87";
      str.spacegroupnumber=87;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 88  I4_{1}./a #88
    if(spacegroup==88 && option==1) {
      str.spacegroup="I4_{1}./a";
      str.spacegrouplabel="#88";
      str.spacegroupnumber=88;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
      }
    }
    if(spacegroup==88 && option==2) {
      str.spacegroup="I4_{1}./a";
      str.spacegrouplabel="#88";
      str.spacegroupnumber=88;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(y+3./4,-x+3./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,x+1./4,-z+1./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 89  P422 #89
    if(spacegroup==89) {
      str.spacegroup="P422";
      str.spacegrouplabel="#89";
      str.spacegroupnumber=89;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 90  P42_{1}2 #90
    if(spacegroup==90) {
      str.spacegroup="P42_{1}2";
      str.spacegrouplabel="#90";
      str.spacegroupnumber=90;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 91  P4_{1}22 #91
    if(spacegroup==91) {
      str.spacegroup="P4_{1}22";
      str.spacegrouplabel="#91";
      str.spacegroupnumber=91;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./4);wa(a,str);
      a.fpos=o+wv(y,-x,z+3./4);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z+3./4);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./4);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 92  P4_{1}2_{1}2 #92
    if(spacegroup==92) {
      str.spacegroup="P4_{1}2_{1}2";
      str.spacegrouplabel="#92";
      str.spacegroupnumber=92;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./4);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+3./4);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./4);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+3./4);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 93  P4_{2}22 #93
    if(spacegroup==93) {
      str.spacegroup="P4_{2}22";
      str.spacegrouplabel="#93";
      str.spacegroupnumber=93;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 94  P4_{2}2_{1}2 #94
    if(spacegroup==94) {
      str.spacegroup="P4_{2}2_{1}2";
      str.spacegrouplabel="#94";
      str.spacegroupnumber=94;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 95  P4_{3}22 #95
    if(spacegroup==95) {
      str.spacegroup="P4_{3}22";
      str.spacegrouplabel="#95";
      str.spacegroupnumber=95;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,z+3./4);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./4);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./4);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+3./4);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 96  P4_{3}2_{1}2 #96
    if(spacegroup==96) {
      str.spacegroup="P4_{3}2_{1}2";
      str.spacegrouplabel="#96";
      str.spacegroupnumber=96;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+3./4);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./4);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+3./4);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./4);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 97  I422 #97
    if(spacegroup==97) {
      str.spacegroup="I422";
      str.spacegrouplabel="#97";
      str.spacegroupnumber=97;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 98  I4_{1}22 #98
    if(spacegroup==98) {
      str.spacegroup="I4_{1}22";
      str.spacegrouplabel="#98";
      str.spacegroupnumber=98;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,-z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 99  P4mm #99
    if(spacegroup==99) {
      str.spacegroup="P4mm";
      str.spacegrouplabel="#99";
      str.spacegroupnumber=99;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 100  P4bm #100
    if(spacegroup==100) {
      str.spacegroup="P4bm";
      str.spacegrouplabel="#100";
      str.spacegroupnumber=100;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 101  P4_{2}cm #101
    if(spacegroup==101) {
      str.spacegroup="P4_{2}cm";
      str.spacegrouplabel="#101";
      str.spacegroupnumber=101;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 102  P4_{2}nm #102
    if(spacegroup==102) {
      str.spacegroup="P4_{2}nm";
      str.spacegrouplabel="#102";
      str.spacegroupnumber=102;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 103  P4cc #103
    if(spacegroup==103) {
      str.spacegroup="P4cc";
      str.spacegrouplabel="#103";
      str.spacegroupnumber=103;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 104  P4nc #104
    if(spacegroup==104) {
      str.spacegroup="P4nc";
      str.spacegrouplabel="#104";
      str.spacegroupnumber=104;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 105  P4_{2}mc #105
    if(spacegroup==105) {
      str.spacegroup="P4_{2}mc";
      str.spacegrouplabel="#105";
      str.spacegroupnumber=105;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 106  P4_{2}bc #106
    if(spacegroup==106) {
      str.spacegroup="P4_{2}bc";
      str.spacegrouplabel="#106";
      str.spacegroupnumber=106;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 107  I4mm #107
    if(spacegroup==107) {
      str.spacegroup="I4mm";
      str.spacegrouplabel="#107";
      str.spacegroupnumber=107;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 108  I4cm #108
    if(spacegroup==108) {
      str.spacegroup="I4cm";
      str.spacegrouplabel="#108";
      str.spacegroupnumber=108;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
	a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 109  I4_{1}md #109
    if(spacegroup==109) {
      str.spacegroup="I4_{1}md";
      str.spacegrouplabel="#109";
      str.spacegroupnumber=109;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,x,z+3./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 110  I4_{1}cd #110
    if(spacegroup==110) {
      str.spacegroup="I4_{1}cd";
      str.spacegrouplabel="#110";
      str.spacegroupnumber=110;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
	a.fpos=o+wv(-y,-x+1./2,z+3./4);wa(a,str);
	a.fpos=o+wv(y+1./2,x,z+1./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 111  P-42m #111
    if(spacegroup==111) {
      str.spacegroup="P-42m";
      str.spacegrouplabel="#111";
      str.spacegroupnumber=111;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 112  P-42c #112
    if(spacegroup==112) {
      str.spacegroup="P-42c";
      str.spacegrouplabel="#112";
      str.spacegroupnumber=112;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 113  P-42_{1}m #113
    if(spacegroup==113) {
      str.spacegroup="P-42_{1}m";
      str.spacegrouplabel="#113";
      str.spacegroupnumber=113;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 114  P-42_{1}c #114
    if(spacegroup==114) {
      str.spacegroup="P-42_{1}c";
      str.spacegrouplabel="#114";
      str.spacegroupnumber=114;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 115  P-4m2 #115
    if(spacegroup==115) {
      str.spacegroup="P-4m2";
      str.spacegrouplabel="#115";
      str.spacegroupnumber=115;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 116  P-4c2 #116
    if(spacegroup==116) {
      str.spacegroup="P-4c2";
      str.spacegrouplabel="#116";
      str.spacegroupnumber=116;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 117  P-4b2 #117
    if(spacegroup==117) {
      str.spacegroup="P-4b2";
      str.spacegrouplabel="#117";
      str.spacegroupnumber=117;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 118  P-4n2 #118
    if(spacegroup==118) {
      str.spacegroup="P-4n2";
      str.spacegrouplabel="#118";
      str.spacegroupnumber=118;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 119  I-4m2 #119
    if(spacegroup==119) {
      str.spacegroup="I-4m2";
      str.spacegrouplabel="#119";
      str.spacegroupnumber=119;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 120  I-4c2 #120
    if(spacegroup==120) {
      str.spacegroup="I-4c2";
      str.spacegrouplabel="#120";
      str.spacegroupnumber=120;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
	a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 121  I-42m #121
    if(spacegroup==121) {
      str.spacegroup="I-42m";
      str.spacegrouplabel="#121";
      str.spacegroupnumber=121;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 122  I-42d #122
    if(spacegroup==122) {
      str.spacegroup="I-42d";
      str.spacegrouplabel="#122";
      str.spacegroupnumber=122;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(y+1./2,x,z+3./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 123  P4./mmm #123
    if(spacegroup==123) {
      str.spacegroup="P4./mmm";
      str.spacegrouplabel="#123";
      str.spacegroupnumber=123;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 124  P4./mcc #124
    if(spacegroup==124) {
      str.spacegroup="P4./mcc";
      str.spacegrouplabel="#124";
      str.spacegroupnumber=124;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 125  P4./nbm #125
    if(spacegroup==125 && option==1) {
      str.spacegroup="P4./nbm";
      str.spacegrouplabel="#125";
      str.spacegroupnumber=125;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);  
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    if(spacegroup==125 && option==2) {
      str.spacegroup="P4./nbm";
      str.spacegrouplabel="#125";
      str.spacegroupnumber=125;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 126  P4./nnc #126
    if(spacegroup==126 && option==1) {
      str.spacegroup="P4./nnc";
      str.spacegrouplabel="#126";
      str.spacegroupnumber=126;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    if(spacegroup==126 && option==2) {
      str.spacegroup="P4./nnc";
      str.spacegrouplabel="#126";
      str.spacegroupnumber=126;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 127  P4./mbm #127
    if(spacegroup==127) {
      str.spacegroup="P4./mbm";
      str.spacegrouplabel="#127";
      str.spacegroupnumber=127;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 128  P4./mnc #128
    if(spacegroup==128) {
      str.spacegroup="P4./mnc";
      str.spacegrouplabel="#128";
      str.spacegroupnumber=128;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 129  P4./nmm #129
    if(spacegroup==129 && option==1) {
      str.spacegroup="P4./nmm";
      str.spacegrouplabel="#129";
      str.spacegroupnumber=129;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    if(spacegroup==129 && option==2) {
      str.spacegroup="P4./nmm";
      str.spacegrouplabel="#129";
      str.spacegroupnumber=129;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 130  P4./ncc #130
    if(spacegroup==130 && option==1) {
      str.spacegroup="P4./ncc";
      str.spacegrouplabel="#130";
      str.spacegroupnumber=130;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);  
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    if(spacegroup==130 && option==2) {
      str.spacegroup="P4./ncc";
      str.spacegrouplabel="#130";
      str.spacegroupnumber=130;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 131  P4_{2}./mmc #131
    if(spacegroup==131) {
      str.spacegroup="P4_{2}./mmc";
      str.spacegrouplabel="#131";
      str.spacegroupnumber=131;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 132  P4_{2}./mcm #132
    if(spacegroup==132) {
      str.spacegroup="P4_{2}./mcm";
      str.spacegrouplabel="#132";
      str.spacegroupnumber=132;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 133  P4_{2}./nbc #133
    if(spacegroup==133 && option==1) {
      str.spacegroup="P4_{2}./nbc";
      str.spacegrouplabel="#133";
      str.spacegroupnumber=133;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);  
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    if(spacegroup==133 && option==2) {
      str.spacegroup="P4_{2}./nbc";
      str.spacegrouplabel="#133";
      str.spacegroupnumber=133;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);  
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 134  P4_{2}./nnm #134
    if(spacegroup==134 && option==1) {
      str.spacegroup="P4_{2}./nnm";
      str.spacegrouplabel="#134";
      str.spacegroupnumber=134;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);  
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    if(spacegroup==134 && option==2) {
      str.spacegroup="P4_{2}./nnm";
      str.spacegrouplabel="#134";
      str.spacegroupnumber=134;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 135  P4_{2}./mbc #135
    if(spacegroup==135) {
      str.spacegroup="P4_{2}./mbc";
      str.spacegrouplabel="#135";
      str.spacegroupnumber=135;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 136  P4_{2}./mnm #136
    if(spacegroup==136) {
      str.spacegroup="P4_{2}./mnm";
      str.spacegrouplabel="#136";
      str.spacegroupnumber=136;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 137  P4_{2}./nmc #137
    if(spacegroup==137 && option==1) {
      str.spacegroup="P4_{2}./nmc";
      str.spacegrouplabel="#137";
      str.spacegroupnumber=137;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
    }
    if(spacegroup==137 && option==2) {
      str.spacegroup="P4_{2}./nmc";
      str.spacegrouplabel="#137";
      str.spacegroupnumber=137;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 138  P4_{2}./ncm #138
    if(spacegroup==138 && option==1) {
      str.spacegroup="P4_{2}./ncm";
      str.spacegrouplabel="#138";
      str.spacegroupnumber=138;
      str.spacegroupoption="origin choice 1";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z);wa(a,str);
    }
    if(spacegroup==138 && option==2) {
      str.spacegroup="P4_{2}./ncm";
      str.spacegrouplabel="#138";
      str.spacegroupnumber=138;
      str.spacegroupoption="origin choice 2";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 139  I4./mmm #139
    if(spacegroup==139) {
      str.spacegroup="I4./mmm";
      str.spacegrouplabel="#139";
      str.spacegroupnumber=139;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 140  I4./mcm #140
    if(spacegroup==140) {
      str.spacegroup="I4./mcm";
      str.spacegrouplabel="#140";
      str.spacegroupnumber=140;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
	a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 141  I4_{1}./amd #141
    if(spacegroup==141 && option==1) {
      str.spacegroup="I4_{1}./amd";
      str.spacegrouplabel="#141";
      str.spacegroupnumber=141;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,-z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(y,x+1./2,z+1./4);wa(a,str);
      }
    }
    if(spacegroup==141 && option==2) {
      str.spacegroup="I4_{1}./amd";
      str.spacegrouplabel="#141";
      str.spacegroupnumber=141;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(y+1./4,x+3./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-y+3./4,-x+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(y+3./4,x+3./4,z+1./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 142  I4_{1}./acd #142
    if(spacegroup==142 && option==1) {
      str.spacegroup="I4_{1}./acd";
      str.spacegrouplabel="#142";
      str.spacegroupnumber=142;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,-z+1./4);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,-z+3./4);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,z);wa(a,str);
	a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x,z+1./4);wa(a,str);
	a.fpos=o+wv(y,x+1./2,z+3./4);wa(a,str);
      }
    }
    if(spacegroup==142 && option==2) {
      str.spacegroup="I4_{1}./acd";
      str.spacegrouplabel="#142";
      str.spacegroupnumber=142;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./4,x+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+3./4,-x+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(y+3./4,x+3./4,z+3./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 143  P3 #143
    if(spacegroup==143) {
      str.spacegroup="P3";
      str.spacegrouplabel="#143";
      str.spacegroupnumber=143;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 144  P3_{1} #144
    if(spacegroup==144) {
      str.spacegroup="P3_{1}";
      str.spacegrouplabel="#144";
      str.spacegroupnumber=144;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 145  P3_{2} #145
    if(spacegroup==145) {
      str.spacegroup="P3_{2}";
      str.spacegrouplabel="#145";
      str.spacegroupnumber=145;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 146  R3 #146
    if(spacegroup==146 && option==1) {
      str.spacegroup="R3";
      str.spacegrouplabel="#146";
      str.spacegroupnumber=146;
      str.spacegroupoption="rhombohedral axes";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
    }
    if(spacegroup==146 && option==2) {
      str.spacegroup="R3";
      str.spacegrouplabel="#146";
      str.spacegroupnumber=146;
      str.spacegroupoption="hexagonal axes";
      for(uint j=1;j<=3;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(2./3,1./3,1./3);
	if(j==3) o=wv(1./3,2./3,2./3);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-y,x-y,z);wa(a,str);
	a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 147  P-3 #147
    if(spacegroup==147) {
      str.spacegroup="P-3";
      str.spacegrouplabel="#147";
      str.spacegroupnumber=147;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 148  R-3 #148
    if(spacegroup==148 && option==1) {
      str.spacegroup="R-3";
      str.spacegrouplabel="#148";
      str.spacegroupnumber=148;
      str.spacegroupoption="rhombohedral axes";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
    }
    if(spacegroup==148 && option==2) {
      str.spacegroup="R-3";
      str.spacegrouplabel="#148";
      str.spacegroupnumber=148;
      str.spacegroupoption="hexagonal axes";
      for(uint j=1;j<=3;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(2./3,1./3,1./3);
	if(j==3) o=wv(1./3,2./3,2./3);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-y,x-y,z);wa(a,str);
	a.fpos=o+wv(-x+y,-x,z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(y,-x+y,-z);wa(a,str);
	a.fpos=o+wv(x-y,x,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 149  P312 #149
    if(spacegroup==149) {
      str.spacegroup="P312";
      str.spacegrouplabel="#149";
      str.spacegroupnumber=149;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 150  P321 #150
    if(spacegroup==150) {
      str.spacegroup="P321";
      str.spacegrouplabel="#150";
      str.spacegroupnumber=150;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 151  P3_{1}12 #151
    if(spacegroup==151) {
      str.spacegroup="P3_{1}12";
      str.spacegrouplabel="#151";
      str.spacegroupnumber=151;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+1./3);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 152  P3_{1}21 #152
    if(spacegroup==152) {
      str.spacegroup="P3_{1}21";
      str.spacegrouplabel="#152";
      str.spacegroupnumber=152;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z+2./3);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+1./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 153  P3_{2}12 #153
    if(spacegroup==153) {
      str.spacegroup="P3_{2}12";
      str.spacegrouplabel="#153";
      str.spacegroupnumber=153;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+2./3);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 154  P3_{2}21 #154
    if(spacegroup==154) {
      str.spacegroup="P3_{2}21";
      str.spacegrouplabel="#154";
      str.spacegroupnumber=154;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z+1./3);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+2./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 155  R32 #155
    if(spacegroup==155 && option==1) {
      str.spacegroup="R32";
      str.spacegrouplabel="#155";
      str.spacegroupnumber=155;
      str.spacegroupoption="rhombohedral axes";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-z,-y,-x);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-z,-y);wa(a,str);
    }
    if(spacegroup==155 && option==2) {
      str.spacegroup="R32";
      str.spacegrouplabel="#155";
      str.spacegroupnumber=155;
      str.spacegroupoption="hexagonal axes";
      for(uint j=1;j<=3;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(2./3,1./3,1./3);
	if(j==3) o=wv(1./3,2./3,2./3);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-y,x-y,z);wa(a,str);
	a.fpos=o+wv(-x+y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(x-y,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 156  P3m1 #156
    if(spacegroup==156) {
      str.spacegroup="P3m1";
      str.spacegrouplabel="#156";
      str.spacegroupnumber=156;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+y,y,z);wa(a,str);
      a.fpos=o+wv(x,x-y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 157  P31m #157
    if(spacegroup==157) {
      str.spacegroup="P31m";
      str.spacegrouplabel="#157";
      str.spacegroupnumber=157;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x-y,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 158  P3c1 #158
    if(spacegroup==158) {
      str.spacegroup="P3c1";
      str.spacegrouplabel="#158";
      str.spacegroupnumber=158;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 159  P31c #159
    if(spacegroup==159) {
      str.spacegroup="P31c";
      str.spacegrouplabel="#159";
      str.spacegroupnumber=159;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 160  R3m #160
    if(spacegroup==160 && option==1) {
      str.spacegroup="R3m";
      str.spacegrouplabel="#160";
      str.spacegroupnumber=160;
      str.spacegroupoption="rhombohedral axes";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(z,y,x);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x,z,y);wa(a,str);
    }
    if(spacegroup==160 && option==2) {
      str.spacegroup="R3m";
      str.spacegrouplabel="#160";
      str.spacegroupnumber=160;
      str.spacegroupoption="hexagonal axes";
      for(uint j=1;j<=3;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(2./3,1./3,1./3);
	if(j==3) o=wv(1./3,2./3,2./3);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-y,x-y,z);wa(a,str);
	a.fpos=o+wv(-x+y,-x,z);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(-x+y,y,z);wa(a,str);
	a.fpos=o+wv(x,x-y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 161  R3c #161
    if(spacegroup==161 && option==1) {
      str.spacegroup="R3c";
      str.spacegrouplabel="#161";
      str.spacegroupnumber=161;
      str.spacegroupoption="rhombohedral axes";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
    }
    if(spacegroup==161 && option==2) {
      str.spacegroup="R3c";
      str.spacegrouplabel="#161";
      str.spacegroupnumber=161;
      str.spacegroupoption="hexagonal axes";
      for(uint j=1;j<=3;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(2./3,1./3,1./3);
	if(j==3) o=wv(1./3,2./3,2./3);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-y,x-y,z);wa(a,str);
	a.fpos=o+wv(-x+y,-x,z);wa(a,str);
	a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
	a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 162  P-31m #162
    if(spacegroup==162) {
      str.spacegroup="P-31m";
      str.spacegrouplabel="#162";
      str.spacegroupnumber=162;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x-y,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 163  P-31c #163
    if(spacegroup==163) {
      str.spacegroup="P-31c";
      str.spacegrouplabel="#163";
      str.spacegroupnumber=163;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 164  P-3m1 #164
    if(spacegroup==164) {
      str.spacegroup="P-3m1";
      str.spacegrouplabel="#164";
      str.spacegroupnumber=164;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+y,y,z);wa(a,str);
      a.fpos=o+wv(x,x-y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 165  P-3c1 #165
    if(spacegroup==165) {
      str.spacegroup="P-3c1";
      str.spacegrouplabel="#165";
      str.spacegroupnumber=165;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 166  R-3m #166
    if(spacegroup==166 && option==1) {
      str.spacegroup="R-3m";
      str.spacegrouplabel="#166";
      str.spacegroupnumber=166;
      str.spacegroupoption="rhombohedral axes";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-z,-y,-x);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-z,-y);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(z,y,x);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x,z,y);wa(a,str);
    }
    if(spacegroup==166 && option==2) {
      str.spacegroup="R-3m";
      str.spacegrouplabel="#166";
      str.spacegroupnumber=166;
      str.spacegroupoption="hexagonal axes";
      for(uint j=1;j<=3;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(2./3,1./3,1./3);
	if(j==3) o=wv(1./3,2./3,2./3);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-y,x-y,z);wa(a,str);
	a.fpos=o+wv(-x+y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(x-y,-y,-z);wa(a,str);
	a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(y,-x+y,-z);wa(a,str);
	a.fpos=o+wv(x-y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(-x+y,y,z);wa(a,str);
	a.fpos=o+wv(x,x-y,z);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 167  R-3c #167
    if(spacegroup==167 && option==1) {
      str.spacegroup="R-3c";
      str.spacegrouplabel="#167";
      str.spacegroupnumber=167;
      str.spacegroupoption="rhombohedral axes";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
    }
    if(spacegroup==167 && option==2) {
      str.spacegroup="R-3c";
      str.spacegrouplabel="#167";
      str.spacegroupnumber=167;
      str.spacegroupoption="hexagonal axes";
      for(uint j=1;j<=3;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(2./3,1./3,1./3);
	if(j==3) o=wv(1./3,2./3,2./3);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-y,x-y,z);wa(a,str);
	a.fpos=o+wv(-x+y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
	a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(y,-x+y,-z);wa(a,str);
	a.fpos=o+wv(x-y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
	a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 168  P6 #168
    if(spacegroup==168) {
      str.spacegroup="P6";
      str.spacegrouplabel="#168";
      str.spacegroupnumber=168;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z);wa(a,str);
      a.fpos=o+wv(x-y,x,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 169  P6_{1} #169
    if(spacegroup==169) {
      str.spacegroup="P6_{1}";
      str.spacegrouplabel="#169";
      str.spacegroupnumber=169;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+5./6);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./6);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 170  P6_{5} #170
    if(spacegroup==170) {
      str.spacegroup="P6_{5}";
      str.spacegrouplabel="#170";
      str.spacegroupnumber=170;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./6);wa(a,str);
      a.fpos=o+wv(x-y,x,z+5./6);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 171  P6_{2} #171
    if(spacegroup==171) {
      str.spacegroup="P6_{2}";
      str.spacegrouplabel="#171";
      str.spacegroupnumber=171;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+2./3);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 172  P6_{4} #172
    if(spacegroup==172) {
      str.spacegroup="P6_{4}";
      str.spacegrouplabel="#172";
      str.spacegroupnumber=172;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./3);wa(a,str);
      a.fpos=o+wv(x-y,x,z+2./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 173  P6_{3} #173
    if(spacegroup==173) {
      str.spacegroup="P6_{3}";
      str.spacegrouplabel="#173";
      str.spacegroupnumber=173;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 174  P-6 #174
    if(spacegroup==174) {
      str.spacegroup="P-6";
      str.spacegrouplabel="#174";
      str.spacegroupnumber=174;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 175  P6./m #175
    if(spacegroup==175) {
      str.spacegroup="P6./m";
      str.spacegrouplabel="#175";
      str.spacegroupnumber=175;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z);wa(a,str);
      a.fpos=o+wv(x-y,x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 176  P6_{3}./m #176
    if(spacegroup==176) {
      str.spacegroup="P6_{3}./m";
      str.spacegrouplabel="#176";
      str.spacegroupnumber=176;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 177  P622 #177
    if(spacegroup==177) {
      str.spacegroup="P622";
      str.spacegrouplabel="#177";
      str.spacegroupnumber=177;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z);wa(a,str);
      a.fpos=o+wv(x-y,x,z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 178  P6_{1}22 #178
    if(spacegroup==178) {
      str.spacegroup="P6_{1}22";
      str.spacegrouplabel="#178";
      str.spacegroupnumber=178;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+5./6);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./6);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./3);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+2./3);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+5./6);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+1./6);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 179  P6_{5}22 #179
    if(spacegroup==179) {
      str.spacegroup="P6_{5}22";
      str.spacegrouplabel="#179";
      str.spacegroupnumber=179;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./6);wa(a,str);
      a.fpos=o+wv(x-y,x,z+5./6);wa(a,str);
      a.fpos=o+wv(y,x,-z+2./3);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+1./3);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./6);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+5./6);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 180  P6_{2}22 #180
    if(spacegroup==180) {
      str.spacegroup="P6_{2}22";
      str.spacegrouplabel="#180";
      str.spacegroupnumber=180;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+1./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+2./3);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./3);wa(a,str);
      a.fpos=o+wv(y,x,-z+2./3);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+1./3);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+2./3);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+1./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 181  P6_{4}22 #181
    if(spacegroup==181) {
      str.spacegroup="P6_{4}22";
      str.spacegrouplabel="#181";
      str.spacegroupnumber=181;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z+2./3);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./3);wa(a,str);
      a.fpos=o+wv(x-y,x,z+2./3);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./3);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+2./3);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./3);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+2./3);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 182  P6_{3}22 #182
    if(spacegroup==182) {
      str.spacegroup="P6_{3}22";
      str.spacegrouplabel="#182";
      str.spacegroupnumber=182;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 183  P6mm #183
    if(spacegroup==183) {
      str.spacegroup="P6mm";
      str.spacegrouplabel="#183";
      str.spacegroupnumber=183;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z);wa(a,str);
      a.fpos=o+wv(x-y,x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+y,y,z);wa(a,str);
      a.fpos=o+wv(x,x-y,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x-y,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 184  P6cc #184
    if(spacegroup==184) {
      str.spacegroup="P6cc";
      str.spacegrouplabel="#184";
      str.spacegroupnumber=184;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z);wa(a,str);
      a.fpos=o+wv(x-y,x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 185  P6_{3}cm #185
    if(spacegroup==185) {
      str.spacegroup="P6_{3}cm";
      str.spacegrouplabel="#185";
      str.spacegroupnumber=185;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x-y,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 186  P6_{3}mc #186
    if(spacegroup==186) {
      str.spacegroup="P6_{3}mc";
      str.spacegrouplabel="#186";
      str.spacegroupnumber=186;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+y,y,z);wa(a,str);
      a.fpos=o+wv(x,x-y,z);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 187  P-6m2 #187
    if(spacegroup==187) {
      str.spacegroup="P-6m2";
      str.spacegrouplabel="#187";
      str.spacegroupnumber=187;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+y,y,z);wa(a,str);
      a.fpos=o+wv(x,x-y,z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 188  P-6c2 #188
    if(spacegroup==188) {
      str.spacegroup="P-6c2";
      str.spacegrouplabel="#188";
      str.spacegroupnumber=188;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 189  P-62m #189
    if(spacegroup==189) {
      str.spacegroup="P-62m";
      str.spacegrouplabel="#189";
      str.spacegroupnumber=189;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x-y,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 190  P-62c #190
    if(spacegroup==190) {
      str.spacegroup="P-62c";
      str.spacegrouplabel="#190";
      str.spacegroupnumber=190;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 191  P6./mmm #191
    if(spacegroup==191) {
      str.spacegroup="P6./mmm";
      str.spacegrouplabel="#191";
      str.spacegroupnumber=191;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z);wa(a,str);
      a.fpos=o+wv(x-y,x,z);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+y,y,z);wa(a,str);
      a.fpos=o+wv(x,x-y,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x-y,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 192  P6./mcc #192
    if(spacegroup==192) {
      str.spacegroup="P6./mcc";
      str.spacegrouplabel="#192";
      str.spacegroupnumber=192;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(y,-x+y,z);wa(a,str);
      a.fpos=o+wv(x-y,x,z);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 193  P6_{3}./mcm #193
    if(spacegroup==193) {
      str.spacegroup="P6_{3}./mcm";
      str.spacegrouplabel="#193";
      str.spacegroupnumber=193;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z);wa(a,str);
      a.fpos=o+wv(x,x-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(x-y,-y,z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 194  P6_{3}./mmc #194
    if(spacegroup==194) {
      str.spacegroup="P6_{3}./mmc";
      str.spacegrouplabel="#194";
      str.spacegroupnumber=194;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-y,x-y,z);wa(a,str);
      a.fpos=o+wv(-x+y,-x,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+y,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(x-y,-y,-z);wa(a,str);
      a.fpos=o+wv(-x,-x+y,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(y,-x+y,-z);wa(a,str);
      a.fpos=o+wv(x-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x-y,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+y,-x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(-x+y,y,z);wa(a,str);
      a.fpos=o+wv(x,x-y,z);wa(a,str);
      a.fpos=o+wv(y,x,z+1./2);wa(a,str);
      a.fpos=o+wv(x-y,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,-x+y,z+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 195  P23 #195
    if(spacegroup==195) {
      str.spacegroup="P23";
      str.spacegrouplabel="#195";
      str.spacegroupnumber=195;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 196  F23 #196
    if(spacegroup==196) {
      str.spacegroup="F23";
      str.spacegrouplabel="#196";
      str.spacegroupnumber=196;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 197  I23 #197
    if(spacegroup==197) {
      str.spacegroup="I23";
      str.spacegrouplabel="#197";
      str.spacegroupnumber=197;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 198  P2_{1}3 #198
    if(spacegroup==198) {
      str.spacegroup="P2_{1}3";
      str.spacegrouplabel="#198";
      str.spacegroupnumber=198;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
      a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 199  I2_{1}3 #199
    if(spacegroup==199) {
      str.spacegroup="I2_{1}3";
      str.spacegrouplabel="#199";
      str.spacegroupnumber=199;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
	a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
	a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 200  Pm-3 #200
    if(spacegroup==200) {
      str.spacegroup="Pm-3";
      str.spacegrouplabel="#200";
      str.spacegroupnumber=200;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,x,y);wa(a,str);
      a.fpos=o+wv(z,x,-y);wa(a,str);
      a.fpos=o+wv(z,-x,y);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,x);wa(a,str);
      a.fpos=o+wv(-y,z,x);wa(a,str);
      a.fpos=o+wv(y,z,-x);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 201  Pn-3 #201
    if(spacegroup==201 && option==1) {
      str.spacegroup="Pn-3";
      str.spacegrouplabel="#201";
      str.spacegroupnumber=201;
      str.spacegroupoption="origin choice 1";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,z+1./2,-x+1./2);wa(a,str);
    }
    if(spacegroup==201 && option==2) {
      str.spacegroup="Pn-3";
      str.spacegrouplabel="#201";
      str.spacegroupnumber=201;
      str.spacegroupoption="origin choice 2";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x+1./2,y);wa(a,str);
      a.fpos=o+wv(-z+1./2,x,-y+1./2);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y+1./2,z,-x+1./2);wa(a,str);
      a.fpos=o+wv(y,-z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z+1./2,x);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,x+1./2,-y);wa(a,str);
      a.fpos=o+wv(z+1./2,-x,y+1./2);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(y+1./2,-z,x+1./2);wa(a,str);
      a.fpos=o+wv(-y,z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,z+1./2,-x);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 202  Fm-3 #202
    if(spacegroup==202) {
      str.spacegroup="Fm-3";
      str.spacegrouplabel="#202";
      str.spacegroupnumber=202;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,x,y);wa(a,str);
	a.fpos=o+wv(z,x,-y);wa(a,str);
	a.fpos=o+wv(z,-x,y);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,x);wa(a,str);
	a.fpos=o+wv(-y,z,x);wa(a,str);
	a.fpos=o+wv(y,z,-x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 203  Fd-3 #203
    if(spacegroup==203 && option==1) {
      str.spacegroup="Fd-3";
      str.spacegrouplabel="#203";
      str.spacegroupnumber=203;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(-x+1./4,-y+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,y+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-y+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,y+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-x+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,x+1./4,y+1./4);wa(a,str);
	a.fpos=o+wv(z+1./4,x+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(z+1./4,-x+1./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-z+1./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-z+1./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,z+1./4,x+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,z+1./4,-x+1./4);wa(a,str);
      }
    }
    if(spacegroup==203 && option==2) {
      str.spacegroup="Fd-3";
      str.spacegrouplabel="#203";
      str.spacegroupnumber=203;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+3./4,-y+3./4,z);wa(a,str);
	a.fpos=o+wv(-x+3./4,y,-z+3./4);wa(a,str);
	a.fpos=o+wv(x,-y+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x+3./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,-x+3./4,y);wa(a,str);
	a.fpos=o+wv(-z+3./4,x,-y+3./4);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y+3./4,z,-x+3./4);wa(a,str);
	a.fpos=o+wv(y,-z+3./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,-z+3./4,x);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./4,y+1./4,-z);wa(a,str);
	a.fpos=o+wv(x+1./4,-y,z+1./4);wa(a,str);
	a.fpos=o+wv(-x,y+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,x+1./4,y+1./4);wa(a,str);
	a.fpos=o+wv(z+1./4,x+1./4,-y);wa(a,str);
	a.fpos=o+wv(z+1./4,-x,y+1./4);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y+1./4,-z,x+1./4);wa(a,str);
	a.fpos=o+wv(-y,z+1./4,x+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,z+1./4,-x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 204  Im-3 #204
    if(spacegroup==204) {
      str.spacegroup="Im-3";
      str.spacegrouplabel="#204";
      str.spacegroupnumber=204;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,x,y);wa(a,str);
	a.fpos=o+wv(z,x,-y);wa(a,str);
	a.fpos=o+wv(z,-x,y);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,x);wa(a,str);
	a.fpos=o+wv(-y,z,x);wa(a,str);
	a.fpos=o+wv(y,z,-x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 205  Pa-3 #205
    if(spacegroup==205) {
      str.spacegroup="Pa-3";
      str.spacegrouplabel="#205";
      str.spacegroupnumber=205;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
      a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z+1./2,x+1./2,y);wa(a,str);
      a.fpos=o+wv(z+1./2,x,-y+1./2);wa(a,str);
      a.fpos=o+wv(z,-x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(y,-z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,z+1./2,x);wa(a,str);
      a.fpos=o+wv(y+1./2,z,-x+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 206  Ia-3 #206
    if(spacegroup==206) {
      str.spacegroup="Ia-3";
      str.spacegrouplabel="#206";
      str.spacegroupnumber=206;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
	a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
	a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,x+1./2,y);wa(a,str);
	a.fpos=o+wv(z+1./2,x,-y+1./2);wa(a,str);
	a.fpos=o+wv(z,-x+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y,-z+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,z+1./2,x);wa(a,str);
	a.fpos=o+wv(y+1./2,z,-x+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 207  P432 #207
    if(spacegroup==207) {
      str.spacegroup="P432";
      str.spacegrouplabel="#207";
      str.spacegroupnumber=207;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(x,z,-y);wa(a,str);
      a.fpos=o+wv(-x,z,y);wa(a,str);
      a.fpos=o+wv(-x,-z,-y);wa(a,str);
      a.fpos=o+wv(x,-z,y);wa(a,str);
      a.fpos=o+wv(z,y,-x);wa(a,str);
      a.fpos=o+wv(z,-y,x);wa(a,str);
      a.fpos=o+wv(-z,y,x);wa(a,str);
      a.fpos=o+wv(-z,-y,-x);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 208  P4_{2}32 #208
    if(spacegroup==208) {
      str.spacegroup="P4_{2}32";
      str.spacegrouplabel="#208";
      str.spacegroupnumber=208;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 209  F432 #209
    if(spacegroup==209) {
      str.spacegroup="F432";
      str.spacegrouplabel="#209";
      str.spacegroupnumber=209;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(x,z,-y);wa(a,str);
	a.fpos=o+wv(-x,z,y);wa(a,str);
	a.fpos=o+wv(-x,-z,-y);wa(a,str);
	a.fpos=o+wv(x,-z,y);wa(a,str);
	a.fpos=o+wv(z,y,-x);wa(a,str);
	a.fpos=o+wv(z,-y,x);wa(a,str);
	a.fpos=o+wv(-z,y,x);wa(a,str);
	a.fpos=o+wv(-z,-y,-x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 210  F4_{1}32 #210
    if(spacegroup==210) {
      str.spacegroup="F4_{1}32";
      str.spacegrouplabel="#210";
      str.spacegroupnumber=210;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x,-y+1./2);wa(a,str);
	a.fpos=o+wv(-z,-x+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,x+1./2,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y+1./2,z+1./2,-x);wa(a,str);
	a.fpos=o+wv(y+1./2,-z,-x+1./2);wa(a,str);
	a.fpos=o+wv(-y,-z+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+3./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(x+3./4,z+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,z+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-z+3./4,y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,y+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(z+1./4,-y+3./4,x+3./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,y+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 211  I432 #211
    if(spacegroup==211) {
      str.spacegroup="I432";
      str.spacegrouplabel="#211";
      str.spacegroupnumber=211;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(x,z,-y);wa(a,str);
	a.fpos=o+wv(-x,z,y);wa(a,str);
	a.fpos=o+wv(-x,-z,-y);wa(a,str);
	a.fpos=o+wv(x,-z,y);wa(a,str);
	a.fpos=o+wv(z,y,-x);wa(a,str);
	a.fpos=o+wv(z,-y,x);wa(a,str);
	a.fpos=o+wv(-z,y,x);wa(a,str);
	a.fpos=o+wv(-z,-y,-x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 212  P4_{3}32 #212
    if(spacegroup==212) {
      str.spacegroup="P4_{3}32";
      str.spacegrouplabel="#212";
      str.spacegroupnumber=212;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
      a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./4,x+3./4,-z+3./4);wa(a,str);
      a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
      a.fpos=o+wv(y+3./4,-x+3./4,z+1./4);wa(a,str);
      a.fpos=o+wv(-y+3./4,x+1./4,z+3./4);wa(a,str);
      a.fpos=o+wv(x+1./4,z+3./4,-y+3./4);wa(a,str);
      a.fpos=o+wv(-x+3./4,z+1./4,y+3./4);wa(a,str);
      a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
      a.fpos=o+wv(x+3./4,-z+3./4,y+1./4);wa(a,str);
      a.fpos=o+wv(z+1./4,y+3./4,-x+3./4);wa(a,str);
      a.fpos=o+wv(z+3./4,-y+3./4,x+1./4);wa(a,str);
      a.fpos=o+wv(-z+3./4,y+1./4,x+3./4);wa(a,str);
      a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 213  P4_{1}32 #213
    if(spacegroup==213) {
      str.spacegroup="P4_{1}32";
      str.spacegrouplabel="#213";
      str.spacegroupnumber=213;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
      a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
      a.fpos=o+wv(y+3./4,x+1./4,-z+1./4);wa(a,str);
      a.fpos=o+wv(-y+3./4,-x+3./4,-z+3./4);wa(a,str);
      a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
      a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
      a.fpos=o+wv(x+3./4,z+1./4,-y+1./4);wa(a,str);
      a.fpos=o+wv(-x+1./4,z+3./4,y+1./4);wa(a,str);
      a.fpos=o+wv(-x+3./4,-z+3./4,-y+3./4);wa(a,str);
      a.fpos=o+wv(x+1./4,-z+1./4,y+3./4);wa(a,str);
      a.fpos=o+wv(z+3./4,y+1./4,-x+1./4);wa(a,str);
      a.fpos=o+wv(z+1./4,-y+1./4,x+3./4);wa(a,str);
      a.fpos=o+wv(-z+1./4,y+3./4,x+1./4);wa(a,str);
      a.fpos=o+wv(-z+3./4,-y+3./4,-x+3./4);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 214  I4_{1}32 #214
    if(spacegroup==214) {
      str.spacegroup="I4_{1}32";
      str.spacegrouplabel="#214";
      str.spacegroupnumber=214;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
	a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
	a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,x+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,-x+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(x+3./4,z+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,z+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,-z+3./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-z+1./4,y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,y+1./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(z+1./4,-y+1./4,x+3./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,y+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,-y+3./4,-x+3./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 215  P-43m #215
    if(spacegroup==215) {
      str.spacegroup="P-43m";
      str.spacegrouplabel="#215";
      str.spacegroupnumber=215;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(x,z,y);wa(a,str);
      a.fpos=o+wv(-x,z,-y);wa(a,str);
      a.fpos=o+wv(-x,-z,y);wa(a,str);
      a.fpos=o+wv(x,-z,-y);wa(a,str);
      a.fpos=o+wv(z,y,x);wa(a,str);
      a.fpos=o+wv(z,-y,-x);wa(a,str);
      a.fpos=o+wv(-z,y,-x);wa(a,str);
      a.fpos=o+wv(-z,-y,x);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 216  F-43m #216
    if(spacegroup==216) {
      str.spacegroup="F-43m";
      str.spacegrouplabel="#216";
      str.spacegroupnumber=216;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(x,z,y);wa(a,str);
	a.fpos=o+wv(-x,z,-y);wa(a,str);
	a.fpos=o+wv(-x,-z,y);wa(a,str);
	a.fpos=o+wv(x,-z,-y);wa(a,str);
	a.fpos=o+wv(z,y,x);wa(a,str);
	a.fpos=o+wv(z,-y,-x);wa(a,str);
	a.fpos=o+wv(-z,y,-x);wa(a,str);
	a.fpos=o+wv(-z,-y,x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 217  I-43m #217
    if(spacegroup==217) {
      str.spacegroup="I-43m";
      str.spacegrouplabel="#217";
      str.spacegroupnumber=217;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(x,z,y);wa(a,str);
	a.fpos=o+wv(-x,z,-y);wa(a,str);
	a.fpos=o+wv(-x,-z,y);wa(a,str);
	a.fpos=o+wv(x,-z,-y);wa(a,str);
	a.fpos=o+wv(z,y,x);wa(a,str);
	a.fpos=o+wv(z,-y,-x);wa(a,str);
	a.fpos=o+wv(-z,y,-x);wa(a,str);
	a.fpos=o+wv(-z,-y,x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 218  P-43n #218
    if(spacegroup==218) {
      str.spacegroup="P-43n";
      str.spacegrouplabel="#218";
      str.spacegroupnumber=218;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 219  F-43c #219
    if(spacegroup==219) {
      str.spacegroup="F-43c";
      str.spacegrouplabel="#219";
      str.spacegroupnumber=219;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 220  I-43d #220
    if(spacegroup==220) {
      str.spacegroup="I-43d";
      str.spacegrouplabel="#220";
      str.spacegroupnumber=220;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
	a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
	a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
	a.fpos=o+wv(y+1./4,x+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+3./4,z+3./4);wa(a,str);
	a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,z+1./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,z+3./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-z+3./4,y+3./4);wa(a,str);
	a.fpos=o+wv(x+3./4,-z+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(z+1./4,y+1./4,x+1./4);wa(a,str);
	a.fpos=o+wv(z+3./4,-y+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,y+3./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-y+3./4,x+3./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 221  Pm-3m #221
    if(spacegroup==221) {
      str.spacegroup="Pm-3m";
      str.spacegrouplabel="#221";
      str.spacegroupnumber=221;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(x,z,-y);wa(a,str);
      a.fpos=o+wv(-x,z,y);wa(a,str);
      a.fpos=o+wv(-x,-z,-y);wa(a,str);
      a.fpos=o+wv(x,-z,y);wa(a,str);
      a.fpos=o+wv(z,y,-x);wa(a,str);
      a.fpos=o+wv(z,-y,x);wa(a,str);
      a.fpos=o+wv(-z,y,x);wa(a,str);
      a.fpos=o+wv(-z,-y,-x);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,x,y);wa(a,str);
      a.fpos=o+wv(z,x,-y);wa(a,str);
      a.fpos=o+wv(z,-x,y);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,x);wa(a,str);
      a.fpos=o+wv(-y,z,x);wa(a,str);
      a.fpos=o+wv(y,z,-x);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-z,y);wa(a,str);
      a.fpos=o+wv(x,-z,-y);wa(a,str);
      a.fpos=o+wv(x,z,y);wa(a,str);
      a.fpos=o+wv(-x,z,-y);wa(a,str);
      a.fpos=o+wv(-z,-y,x);wa(a,str);
      a.fpos=o+wv(-z,y,-x);wa(a,str);
      a.fpos=o+wv(z,-y,-x);wa(a,str);
      a.fpos=o+wv(z,y,x);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 222  Pn-3n #222
    if(spacegroup==222 && option==1) {
      str.spacegroup="Pn-3n";
      str.spacegrouplabel="#222";
      str.spacegroupnumber=222;
      str.spacegroupoption="origin choice 1";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y,x,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(y,-x,z);wa(a,str);
      a.fpos=o+wv(-y,x,z);wa(a,str);
      a.fpos=o+wv(x,z,-y);wa(a,str);
      a.fpos=o+wv(-x,z,y);wa(a,str);
      a.fpos=o+wv(-x,-z,-y);wa(a,str);
      a.fpos=o+wv(x,-z,y);wa(a,str);
      a.fpos=o+wv(z,y,-x);wa(a,str);
      a.fpos=o+wv(z,-y,x);wa(a,str);
      a.fpos=o+wv(-z,y,x);wa(a,str);
      a.fpos=o+wv(-z,-y,-x);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
    }
    if(spacegroup==222 && option==2) {
      str.spacegroup="Pn-3n";
      str.spacegrouplabel="#222";
      str.spacegroupnumber=222;
      str.spacegroupoption="origin choice 2";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x+1./2,y);wa(a,str);
      a.fpos=o+wv(-z+1./2,x,-y+1./2);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y+1./2,z,-x+1./2);wa(a,str);
      a.fpos=o+wv(y,-z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z+1./2,x);wa(a,str);
      a.fpos=o+wv(y,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,z);wa(a,str);
      a.fpos=o+wv(x,z,-y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,z,y);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(x,-z+1./2,y);wa(a,str);
      a.fpos=o+wv(z,y,-x+1./2);wa(a,str);
      a.fpos=o+wv(z,-y+1./2,x);wa(a,str);
      a.fpos=o+wv(-z+1./2,y,x);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,x+1./2,-y);wa(a,str);
      a.fpos=o+wv(z+1./2,-x,y+1./2);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(y+1./2,-z,x+1./2);wa(a,str);
      a.fpos=o+wv(-y,z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,z+1./2,-x);wa(a,str);
      a.fpos=o+wv(-y,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-z,y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-z,-y);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x,z+1./2,-y);wa(a,str);
      a.fpos=o+wv(-z,-y,x+1./2);wa(a,str);
      a.fpos=o+wv(-z,y+1./2,-x);wa(a,str);
      a.fpos=o+wv(z+1./2,-y,-x);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 223  Pm-3n #223
    if(spacegroup==223) {
      str.spacegroup="Pm-3n";
      str.spacegrouplabel="#223";
      str.spacegroupnumber=223;
      str.spacegroupoption="";
      o=wv(0,0,0);
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,z);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,x,y);wa(a,str);
      a.fpos=o+wv(z,x,-y);wa(a,str);
      a.fpos=o+wv(z,-x,y);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,x);wa(a,str);
      a.fpos=o+wv(-y,z,x);wa(a,str);
      a.fpos=o+wv(y,z,-x);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 224  Pn-3m #224
    if(spacegroup==224 && option==1) {
      str.spacegroup="Pn-3m";
      str.spacegrouplabel="#224";
      str.spacegroupnumber=224;
      str.spacegroupoption="origin choice 1";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x,-y,z);wa(a,str);
      a.fpos=o+wv(-x,y,-z);wa(a,str);
      a.fpos=o+wv(x,-y,-z);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,-x,y);wa(a,str);
      a.fpos=o+wv(-z,x,-y);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y,z,-x);wa(a,str);
      a.fpos=o+wv(y,-z,-x);wa(a,str);
      a.fpos=o+wv(-y,-z,x);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,-y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,-x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,-z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-y,-x,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(-y,x,-z);wa(a,str);
      a.fpos=o+wv(y,-x,-z);wa(a,str);
      a.fpos=o+wv(-x,-z,y);wa(a,str);
      a.fpos=o+wv(x,-z,-y);wa(a,str);
      a.fpos=o+wv(x,z,y);wa(a,str);
      a.fpos=o+wv(-x,z,-y);wa(a,str);
      a.fpos=o+wv(-z,-y,x);wa(a,str);
      a.fpos=o+wv(-z,y,-x);wa(a,str);
      a.fpos=o+wv(z,-y,-x);wa(a,str);
      a.fpos=o+wv(z,y,x);wa(a,str);
    }
    if(spacegroup==224 && option==2) {
      str.spacegroup="Pn-3m";
      str.spacegrouplabel="#224";
      str.spacegroupnumber=224;
      str.spacegroupoption="origin choice 2";
      a.fpos=o+wv(x,y,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,-y+1./2,z);wa(a,str);
      a.fpos=o+wv(-x+1./2,y,-z+1./2);wa(a,str);
      a.fpos=o+wv(x,-y+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(z,x,y);wa(a,str);
      a.fpos=o+wv(z,-x+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-x+1./2,y);wa(a,str);
      a.fpos=o+wv(-z+1./2,x,-y+1./2);wa(a,str);
      a.fpos=o+wv(y,z,x);wa(a,str);
      a.fpos=o+wv(-y+1./2,z,-x+1./2);wa(a,str);
      a.fpos=o+wv(y,-z+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(-y+1./2,-z+1./2,x);wa(a,str);
      a.fpos=o+wv(y+1./2,x+1./2,-z);wa(a,str);
      a.fpos=o+wv(-y,-x,-z);wa(a,str);
      a.fpos=o+wv(y+1./2,-x,z+1./2);wa(a,str);
      a.fpos=o+wv(-y,x+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(x+1./2,z+1./2,-y);wa(a,str);
      a.fpos=o+wv(-x,z+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(-x,-z,-y);wa(a,str);
      a.fpos=o+wv(x+1./2,-z,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,y+1./2,-x);wa(a,str);
      a.fpos=o+wv(z+1./2,-y,x+1./2);wa(a,str);
      a.fpos=o+wv(-z,y+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(-z,-y,-x);wa(a,str);
      a.fpos=o+wv(-x,-y,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,y+1./2,-z);wa(a,str);
      a.fpos=o+wv(x+1./2,-y,z+1./2);wa(a,str);
      a.fpos=o+wv(-x,y+1./2,z+1./2);wa(a,str);
      a.fpos=o+wv(-z,-x,-y);wa(a,str);
      a.fpos=o+wv(-z,x+1./2,y+1./2);wa(a,str);
      a.fpos=o+wv(z+1./2,x+1./2,-y);wa(a,str);
      a.fpos=o+wv(z+1./2,-x,y+1./2);wa(a,str);
      a.fpos=o+wv(-y,-z,-x);wa(a,str);
      a.fpos=o+wv(y+1./2,-z,x+1./2);wa(a,str);
      a.fpos=o+wv(-y,z+1./2,x+1./2);wa(a,str);
      a.fpos=o+wv(y+1./2,z+1./2,-x);wa(a,str);
      a.fpos=o+wv(-y+1./2,-x+1./2,z);wa(a,str);
      a.fpos=o+wv(y,x,z);wa(a,str);
      a.fpos=o+wv(-y+1./2,x,-z+1./2);wa(a,str);
      a.fpos=o+wv(y,-x+1./2,-z+1./2);wa(a,str);
      a.fpos=o+wv(-x+1./2,-z+1./2,y);wa(a,str);
      a.fpos=o+wv(x,-z+1./2,-y+1./2);wa(a,str);
      a.fpos=o+wv(x,z,y);wa(a,str);
      a.fpos=o+wv(-x+1./2,z,-y+1./2);wa(a,str);
      a.fpos=o+wv(-z+1./2,-y+1./2,x);wa(a,str);
      a.fpos=o+wv(-z+1./2,y,-x+1./2);wa(a,str);
      a.fpos=o+wv(z,-y+1./2,-x+1./2);wa(a,str);
      a.fpos=o+wv(z,y,x);wa(a,str);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 225  Fm-3m #225
    if(spacegroup==225) {
      str.spacegroup="Fm-3m";
      str.spacegrouplabel="#225";
      str.spacegroupnumber=225;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(x,z,-y);wa(a,str);
	a.fpos=o+wv(-x,z,y);wa(a,str);
	a.fpos=o+wv(-x,-z,-y);wa(a,str);
	a.fpos=o+wv(x,-z,y);wa(a,str);
	a.fpos=o+wv(z,y,-x);wa(a,str);
	a.fpos=o+wv(z,-y,x);wa(a,str);
	a.fpos=o+wv(-z,y,x);wa(a,str);
	a.fpos=o+wv(-z,-y,-x);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,x,y);wa(a,str);
	a.fpos=o+wv(z,x,-y);wa(a,str);
	a.fpos=o+wv(z,-x,y);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,x);wa(a,str);
	a.fpos=o+wv(-y,z,x);wa(a,str);
	a.fpos=o+wv(y,z,-x);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-x,-z,y);wa(a,str);
	a.fpos=o+wv(x,-z,-y);wa(a,str);
	a.fpos=o+wv(x,z,y);wa(a,str);
	a.fpos=o+wv(-x,z,-y);wa(a,str);
	a.fpos=o+wv(-z,-y,x);wa(a,str);
	a.fpos=o+wv(-z,y,-x);wa(a,str);
	a.fpos=o+wv(z,-y,-x);wa(a,str);
	a.fpos=o+wv(z,y,x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 226  Fm-3c #226
    if(spacegroup==226) {
      str.spacegroup="Fm-3c";
      str.spacegrouplabel="#226";
      str.spacegroupnumber=226;
      str.spacegroupoption="";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,-y+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,y+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,x,y);wa(a,str);
	a.fpos=o+wv(z,x,-y);wa(a,str);
	a.fpos=o+wv(z,-x,y);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,x);wa(a,str);
	a.fpos=o+wv(-y,z,x);wa(a,str);
	a.fpos=o+wv(y,z,-x);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,-z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,-y+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,-y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 227  Fd-3m #227
    if(spacegroup==227 && option==1) {
      str.spacegroup="Fd-3m";
      str.spacegrouplabel="#227";
      str.spacegroupnumber=227;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x,-y+1./2);wa(a,str);
	a.fpos=o+wv(-z,-x+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,x+1./2,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y+1./2,z+1./2,-x);wa(a,str);
	a.fpos=o+wv(y+1./2,-z,-x+1./2);wa(a,str);
	a.fpos=o+wv(-y,-z+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+3./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(x+3./4,z+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,z+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-z+3./4,y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,y+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(z+1./4,-y+3./4,x+3./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,y+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-y+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,y+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(x+3./4,-y+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,y+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-x+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,x+1./4,y+3./4);wa(a,str);
	a.fpos=o+wv(z+1./4,x+3./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,-x+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-z+1./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(y+3./4,-z+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,z+1./4,x+3./4);wa(a,str);
	a.fpos=o+wv(y+1./4,z+3./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x,z+1./2);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
	a.fpos=o+wv(-y,x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-x+1./2,-z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-z,y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-z+1./2,-y);wa(a,str);
	a.fpos=o+wv(x,z,y);wa(a,str);
	a.fpos=o+wv(-x,z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,-y,x+1./2);wa(a,str);
	a.fpos=o+wv(-z,y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,-y+1./2,-x);wa(a,str);
	a.fpos=o+wv(z,y,x);wa(a,str);
      }
    }
    if(spacegroup==227 && option==2) {
      //   cerr << "aflow_wyckoff.cpp: spacegroup=" << spacegroup << " option=" << option << endl;
      str.spacegroup="Fd-3m";
      str.spacegrouplabel="#227";
      str.spacegroupnumber=227;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+3./4,-y+1./4,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./4,y+1./2,-z+3./4);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+3./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x+3./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,-x+1./4,y+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./4,x+1./2,-y+3./4);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y+1./4,z+1./2,-x+3./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-z+3./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,-z+1./4,x+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,x+1./4,-z+1./2);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+1./2,z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./2,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(x+3./4,z+1./4,-y+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,z+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x,-z,-y);wa(a,str);
	a.fpos=o+wv(x+1./4,-z+1./2,y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,y+1./4,-x+1./2);wa(a,str);
	a.fpos=o+wv(z+1./4,-y+1./2,x+3./4);wa(a,str);
	a.fpos=o+wv(-z+1./2,y+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-z,-y,-x);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./4,y+3./4,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+3./4,-y+1./2,z+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,x+1./4,y+3./4);wa(a,str);
	a.fpos=o+wv(z+1./4,x+3./4,-y+1./2);wa(a,str);
	a.fpos=o+wv(z+3./4,-x+1./2,y+1./4);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y+3./4,-z+1./2,x+1./4);wa(a,str);
	a.fpos=o+wv(-y+1./2,z+1./4,x+3./4);wa(a,str);
	a.fpos=o+wv(y+1./4,z+3./4,-x+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+3./4,z+1./2);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+1./2,-z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-z+3./4,y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-z+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(x,z,y);wa(a,str);
	a.fpos=o+wv(-x+3./4,z+1./2,-y+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-y+3./4,x+1./2);wa(a,str);
	a.fpos=o+wv(-z+3./4,y+1./2,-x+1./4);wa(a,str);
	a.fpos=o+wv(z+1./2,-y+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(z,y,x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 228  Fd-3c #228
    if(spacegroup==228 && option==1) {
      str.spacegroup="Fd-3c";
      str.spacegrouplabel="#228";
      str.spacegroupnumber=228;
      str.spacegroupoption="origin choice 1";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,-y,-z+1./2);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x,-y+1./2);wa(a,str);
	a.fpos=o+wv(-z,-x+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-z+1./2,x+1./2,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y+1./2,z+1./2,-x);wa(a,str);
	a.fpos=o+wv(y+1./2,-z,-x+1./2);wa(a,str);
	a.fpos=o+wv(-y,-z+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+3./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(x+3./4,z+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,z+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-z+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-z+3./4,y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,y+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(z+1./4,-y+3./4,x+3./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,y+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-y+1./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,-y+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(x+3./4,y+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-y+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,y+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,-x+3./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,x+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(z+3./4,x+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(z+1./4,-x+1./4,y+3./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,-z+3./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-z+1./4,x+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,z+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(y+3./4,z+1./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(-y,-x+1./2,z);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,x,-z);wa(a,str);
	a.fpos=o+wv(y,-x,-z+1./2);wa(a,str);
	a.fpos=o+wv(-x,-z+1./2,y);wa(a,str);
	a.fpos=o+wv(x,-z,-y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,z,-y);wa(a,str);
	a.fpos=o+wv(-z,-y+1./2,x);wa(a,str);
	a.fpos=o+wv(-z+1./2,y,-x);wa(a,str);
	a.fpos=o+wv(z,-y,-x+1./2);wa(a,str);
	a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
      }
    }
    if(spacegroup==228 && option==2) {
      str.spacegroup="Fd-3c";
      str.spacegrouplabel="#228";
      str.spacegroupnumber=228;
      str.spacegroupoption="origin choice 2";
      for(uint j=1;j<=4;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(0,1./2,1./2);
	if(j==3) o=wv(1./2,0,1./2);
	if(j==4) o=wv(1./2,1./2,0);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./4,-y+3./4,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+3./4,y+1./2,-z+1./4);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-x+3./4,y+1./2);wa(a,str);
	a.fpos=o+wv(-z+3./4,x+1./2,-y+1./4);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y+3./4,z+1./2,-x+1./4);wa(a,str);
	a.fpos=o+wv(y+1./2,-z+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,-z+3./4,x+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,x+1./4,-z);wa(a,str);
	a.fpos=o+wv(-y+1./2,-x+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(y+1./4,-x,z+3./4);wa(a,str);
	a.fpos=o+wv(-y,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(x+3./4,z+1./4,-y);wa(a,str);
	a.fpos=o+wv(-x,z+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,-z+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(x+1./4,-z,y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,y+1./4,-x);wa(a,str);
	a.fpos=o+wv(z+1./4,-y,x+3./4);wa(a,str);
	a.fpos=o+wv(-z,y+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./2,-y+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+3./4,y+1./4,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./4,-y+1./2,z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,x+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(z+3./4,x+1./4,-y+1./2);wa(a,str);
	a.fpos=o+wv(z+1./4,-x+1./2,y+3./4);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y+1./4,-z+1./2,x+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./2,z+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(y+3./4,z+1./4,-x+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+3./4,z);wa(a,str);
	a.fpos=o+wv(y+1./2,x+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-y+3./4,x,-z+1./4);wa(a,str);
	a.fpos=o+wv(y,-x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-z+3./4,y);wa(a,str);
	a.fpos=o+wv(x,-z+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(x+1./2,z+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-x+3./4,z,-y+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-y+3./4,x);wa(a,str);
	a.fpos=o+wv(-z+3./4,y,-x+1./4);wa(a,str);
	a.fpos=o+wv(z,-y+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(z+1./2,y+1./2,x+1./2);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 229  Im-3m #229
    if(spacegroup==229) {
      str.spacegroup="Im-3m";
      str.spacegrouplabel="#229";
      str.spacegroupnumber=229;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,-x,y);wa(a,str);
	a.fpos=o+wv(-z,x,-y);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,-x);wa(a,str);
	a.fpos=o+wv(-y,-z,x);wa(a,str);
	a.fpos=o+wv(y,x,-z);wa(a,str);
	a.fpos=o+wv(-y,-x,-z);wa(a,str);
	a.fpos=o+wv(y,-x,z);wa(a,str);
	a.fpos=o+wv(-y,x,z);wa(a,str);
	a.fpos=o+wv(x,z,-y);wa(a,str);
	a.fpos=o+wv(-x,z,y);wa(a,str);
	a.fpos=o+wv(-x,-z,-y);wa(a,str);
	a.fpos=o+wv(x,-z,y);wa(a,str);
	a.fpos=o+wv(z,y,-x);wa(a,str);
	a.fpos=o+wv(z,-y,x);wa(a,str);
	a.fpos=o+wv(-z,y,x);wa(a,str);
	a.fpos=o+wv(-z,-y,-x);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x,y,-z);wa(a,str);
	a.fpos=o+wv(x,-y,z);wa(a,str);
	a.fpos=o+wv(-x,y,z);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z,x,y);wa(a,str);
	a.fpos=o+wv(z,x,-y);wa(a,str);
	a.fpos=o+wv(z,-x,y);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y,-z,x);wa(a,str);
	a.fpos=o+wv(-y,z,x);wa(a,str);
	a.fpos=o+wv(y,z,-x);wa(a,str);
	a.fpos=o+wv(-y,-x,z);wa(a,str);
	a.fpos=o+wv(y,x,z);wa(a,str);
	a.fpos=o+wv(-y,x,-z);wa(a,str);
	a.fpos=o+wv(y,-x,-z);wa(a,str);
	a.fpos=o+wv(-x,-z,y);wa(a,str);
	a.fpos=o+wv(x,-z,-y);wa(a,str);
	a.fpos=o+wv(x,z,y);wa(a,str);
	a.fpos=o+wv(-x,z,-y);wa(a,str);
	a.fpos=o+wv(-z,-y,x);wa(a,str);
	a.fpos=o+wv(-z,y,-x);wa(a,str);
	a.fpos=o+wv(z,-y,-x);wa(a,str);
	a.fpos=o+wv(z,y,x);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // 230  Ia-3d #230
    if(spacegroup==230) {
      str.spacegroup="Ia-3d";
      str.spacegrouplabel="#230";
      str.spacegroupnumber=230;
      str.spacegroupoption="";
      for(uint j=1;j<=2;j++) {
	if(j==1) o=wv(0,0,0);
	if(j==2) o=wv(1./2,1./2,1./2);
	a.fpos=o+wv(x,y,z);wa(a,str);
	a.fpos=o+wv(-x+1./2,-y,z+1./2);wa(a,str);
	a.fpos=o+wv(-x,y+1./2,-z+1./2);wa(a,str);
	a.fpos=o+wv(x+1./2,-y+1./2,-z);wa(a,str);
	a.fpos=o+wv(z,x,y);wa(a,str);
	a.fpos=o+wv(z+1./2,-x+1./2,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,-x,y+1./2);wa(a,str);
	a.fpos=o+wv(-z,x+1./2,-y+1./2);wa(a,str);
	a.fpos=o+wv(y,z,x);wa(a,str);
	a.fpos=o+wv(-y,z+1./2,-x+1./2);wa(a,str);
	a.fpos=o+wv(y+1./2,-z+1./2,-x);wa(a,str);
	a.fpos=o+wv(-y+1./2,-z,x+1./2);wa(a,str);
	a.fpos=o+wv(y+3./4,x+1./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,-x+3./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(y+1./4,-x+1./4,z+3./4);wa(a,str);
	a.fpos=o+wv(-y+1./4,x+3./4,z+1./4);wa(a,str);
	a.fpos=o+wv(x+3./4,z+1./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,z+3./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,-z+3./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(x+1./4,-z+1./4,y+3./4);wa(a,str);
	a.fpos=o+wv(z+3./4,y+1./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(z+1./4,-y+1./4,x+3./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,y+3./4,x+1./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,-y+3./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(-x,-y,-z);wa(a,str);
	a.fpos=o+wv(x+1./2,y,-z+1./2);wa(a,str);
	a.fpos=o+wv(x,-y+1./2,z+1./2);wa(a,str);
	a.fpos=o+wv(-x+1./2,y+1./2,z);wa(a,str);
	a.fpos=o+wv(-z,-x,-y);wa(a,str);
	a.fpos=o+wv(-z+1./2,x+1./2,y);wa(a,str);
	a.fpos=o+wv(z+1./2,x,-y+1./2);wa(a,str);
	a.fpos=o+wv(z,-x+1./2,y+1./2);wa(a,str);
	a.fpos=o+wv(-y,-z,-x);wa(a,str);
	a.fpos=o+wv(y,-z+1./2,x+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./2,z+1./2,x);wa(a,str);
	a.fpos=o+wv(y+1./2,z,-x+1./2);wa(a,str);
	a.fpos=o+wv(-y+1./4,-x+3./4,z+3./4);wa(a,str);
	a.fpos=o+wv(y+1./4,x+1./4,z+1./4);wa(a,str);
	a.fpos=o+wv(-y+3./4,x+3./4,-z+1./4);wa(a,str);
	a.fpos=o+wv(y+3./4,-x+1./4,-z+3./4);wa(a,str);
	a.fpos=o+wv(-x+1./4,-z+3./4,y+3./4);wa(a,str);
	a.fpos=o+wv(x+3./4,-z+1./4,-y+3./4);wa(a,str);
	a.fpos=o+wv(x+1./4,z+1./4,y+1./4);wa(a,str);
	a.fpos=o+wv(-x+3./4,z+3./4,-y+1./4);wa(a,str);
	a.fpos=o+wv(-z+1./4,-y+3./4,x+3./4);wa(a,str);
	a.fpos=o+wv(-z+3./4,y+3./4,-x+1./4);wa(a,str);
	a.fpos=o+wv(z+3./4,-y+1./4,-x+3./4);wa(a,str);
	a.fpos=o+wv(z+1./4,y+1./4,x+1./4);wa(a,str);
      }
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
  }
  // cycl atoms
  // cerr << "DONE atoms=" << str.atoms.size() << endl;
  //DONE
  str.MakeBasis();
  return str;
}

//----------------------------------------------------------------------------

#endif//_WYCKOFF_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2008           *
// *                                                                         *
// ***************************************************************************


