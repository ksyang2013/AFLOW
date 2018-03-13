// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2011               *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_kesong.cpp
// functions written by KESONG YANG
// 2010-2011: kesong.yang@gmail.com

#ifndef _AFLOW_CONTRIB_KESONG_MAIN_CPP_
#define _AFLOW_CONTRIB_KESONG_MAIN_CPP_

#include "aflow_contrib_kesong.h"

// ***************************************************************************
// pflow::BZMAX
// ***************************************************************************
namespace pflow {
  void BZMAX(istream & input) {
    //Usage: aflow --BZmax < POSCAR

    xstructure str_in(input,IOAFLOW_AUTO);
    xstructure str_sp,str_sc;
    LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc); //This functions format the str_in, and stores data in str_sp and str_sc
    double grid=16.0;
    string lattice, stmp, line, kpath_string;
    stringstream KPATH, strline;
    bool foundBZ;
    xvector<double> b1(3); xvector<double> b2(3); xvector<double> b3(3); //reciprocal lattice
    xmatrix<double> klattice;    //Declare a matrix to store reciprocal lattice

    lattice=str_sp.bravais_lattice_variation_type;
    kpath_string=LATTICE::KPOINTS_Directions(lattice, str_sp.lattice, grid,str_sp.iomode,foundBZ);
    KPATH.str(kpath_string);  //converting kpath from string type to stringstream
    klattice=str_sp.klattice;

    b1[1]=str_sp.klattice(1,1); b1[2]=str_sp.klattice(1,2); b1[3]= str_sp.klattice(1,3);
    b2[1]=str_sp.klattice(2,1); b2[2]=str_sp.klattice(2,2); b2[3]= str_sp.klattice(2,3);
    b3[1]=str_sp.klattice(3,1); b3[2]=str_sp.klattice(3,2); b3[3]= str_sp.klattice(3,3);

    //////////////READING KPOINTS//////////////////////////////////////////////////////////
    for (int i=0; i<4; i++) {getline(KPATH, stmp);} //Read the first 4 lines
    vector<vector<double> > kpoints;
    vector<string> kpointslabel;
    int count=0, j=0; //count is the number of rows of kpoints
    while (getline(KPATH, line)) {
      if(aurostd::CountWordsinString(line)>=3) {
	vector<double> kpt_tmp;
	double a1, a2, a3;
	string a5;
	strline.clear();
	strline.str(line);
	strline >> a1 >> a2 >> a3 >> stmp >> a5;
	kpt_tmp.push_back(a1);
	kpt_tmp.push_back(a2);
	kpt_tmp.push_back(a3);
	kpoints.push_back(kpt_tmp);
	kpointslabel.push_back(a5);
	j++;
      }
      count =j;
    }

    //////////////READING KPOINTS//////////////////////////////////////////////////////////
    //Converting kpoints into Cartesian
    int NKPOINTS=count;
    vector<vector<double> > kpcart;
    vector<double> klinedistance;
    //vector<double> klinedistance(NKPOINTS);
    for (int i=0; i<NKPOINTS; i++) {
      vector<double> kpcart_tmp;
      double c1=0, c2=0, c3=0, distance;
      c1=kpoints[i][0]*b1[1]+kpoints[i][1]*b2[1]+kpoints[i][2]*b3[1];
      c2=kpoints[i][0]*b1[2]+kpoints[i][1]*b2[2]+kpoints[i][2]*b3[2];
      c3=kpoints[i][0]*b1[3]+kpoints[i][1]*b2[3]+kpoints[i][2]*b3[3];
      kpcart_tmp.push_back(c1);
      kpcart_tmp.push_back(c2);
      kpcart_tmp.push_back(c3);
      kpcart.push_back(kpcart_tmp);
      distance=sqrt(c1*c1+c2*c2+c3*c3);
      klinedistance.push_back(distance);
    }
    //Sorting
    for (int i=0; i<(NKPOINTS-1); i++) {
      for (int j=i+1; j<NKPOINTS; j++) {
	if(klinedistance[j] > klinedistance[i]) {
	  //Sorting kline distance
	  double tmp;
	  tmp=klinedistance[i];
	  klinedistance[i]=klinedistance[j];
	  klinedistance[j]=tmp;
	  //Meanwhile change the order of the kpoints label
	  string kplabel_tmp;
	  kplabel_tmp=kpointslabel[i];
	  kpointslabel[i]=kpointslabel[j];
	  kpointslabel[j]=kplabel_tmp;
	  //Meanwhile change the order of the kpoints label
	  double k0, k1, k2;
	  k0=kpoints[i][0];
	  k1=kpoints[i][1];
	  k2=kpoints[i][3];
	  kpoints[i][0]=kpoints[j][0];
	  kpoints[i][1]=kpoints[j][1];
	  kpoints[i][2]=kpoints[j][2];
	  kpoints[j][0]=k0;
	  kpoints[j][1]=k1;
	  kpoints[j][2]=k2;
	}
      }
    }
    //Formating output
    cout.precision(8);
    cout << "Bravais Lattice:   " << lattice << endl;
    cout << "              KPOINTS              DISTANCE   TYPE" << endl;
    cout << kpoints[0][0] <<" ";
    cout << kpoints[0][1] <<" ";
    cout << kpoints[0][2]<<"  ";
    cout << klinedistance[0] <<"   ";
    cout << kpointslabel[0] << endl;

    for (int i=1; i<(NKPOINTS-1); i++) {
      if(kpointslabel[i].compare(kpointslabel[i-1])!=0) {
	cout << kpoints[i][0] <<" ";
	cout << kpoints[i][1] <<" ";
	cout << kpoints[i][2]<<"  ";
	cout << klinedistance[i] <<"   ";
	cout << kpointslabel[i] << endl;
      }
    }
  }
}

#endif
// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG - Duke University 2010-2011              *
// *                                                                         *
// ***************************************************************************

