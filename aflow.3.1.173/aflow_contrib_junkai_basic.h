// ***************************************************************************
// *                                                                         *
// *           Aflow JUNKAI XUE - Duke University 2003-2018                  *
// *                                                                         *
// ***************************************************************************
// Written by Junkai: 2010-2011
// junkai.xue@duke.edu

#include<iostream>


#ifndef _JUNKAI_BASIC_H_
#define _JUNKAI_BASIC_H_


// BASIC
void SORT(vector<int>& a);
void SORT(vector<double>& a);
void SORTDecrease(vector<uint>& a);
double NORM(const xvector<double>& v);
double COS(const xvector<double>& v1, const xvector<double>& v2);
double DotPro(const xvector<double>& v1, const xvector<double>& v2);
void CrossPro(const xvector<double>& a, const xvector<double>& b, xvector<double>& resu);
void  POSCAR_OutputFun(xstructure& a);
void  POSCAR_Reread(xstructure& a);
string POSCAR_Output(xstructure& a);
bool LinearIndependence(const xmatrix<double>& a);
bool VectSort(const xvector<double>& a, const xvector<double>& b);
bool AtomSort(const _atom& a, const _atom& b);
bool DsortF(double i, double j);
void SpeciesSeparate(string& a, vector<string>& spe, vector<double>& n);

// PROTOCHECK
void CheckProtoTypeLow(void);//(vector<string>& argv);
void CheckProtoTypeMedian(void);
void CheckProtoTypeHigh(void);
void CheckProtoType_Angle(void);
void CheckProtoTypeMutiLow(void);
void CheckProtoTypeMutiMedian(void);
void CheckProtoTypeMutiHigh(void);
void CheckProtoTypeMutiLow(ostream& oss,vector<string>& result);
void CheckProtoTypeMutiMedian(ostream& oss,vector<string>& result);
void CheckProtoTypeMutiHigh(ostream& oss,vector<string>& result);
void CheckProtoType(xstructure a, int LEV, ostream& oss,vector<string>& result);
void CheckProtoType_Angle(xstructure a, ostream& oss,vector<string>& result);
void CheckProtoTypeLow(ostream& oss,vector<string>& result);//(vector<string>& argv);
void CheckProtoTypeMedian(ostream& oss,vector<string>& result);
void CheckProtoTypeHigh(ostream& oss,vector<string>& result);

// CPROTO
void ProtoTypeInitial(vector<string>& argv);
void ProtoClassification(int LEV);
void ProtoClassificationHigh(void);
void ProtoClassificationMedian(void);
void ProtoClassificationLow(void);
void ProtoClassificationCoarse(void);
void ProtoClassificationLow(vector<string>& argv);
void PROTOOUTFUN(vector<vector<xstructure> >& PROTO, vector<string>& argv, int& size, string& TYPE);  
void ProtoTypeCompare( vector<string> argv);
void Databasegenerate(void);//(vector<string> argv);
bool CoordCompare(const xvector<double>& v1, const xvector<double>& v2,int& LEV);
bool CoordCompare(const xvector<double>& v,const _atom& atom2, int& LEV);
bool CoordCompare(const _atom& atom1, const deque<_atom>& V, int& LEV, vector<int>& label2);
bool CoordCompare(const deque<_atom>& V1, const deque<_atom>& V2, int& LEV);
bool CoordCompare(const xstructure& a, const xstructure& b, int& LEV);
void AtomicEnvironment(const deque<deque<_atom> >& neigh_mat, vector<int>& shellnum);
bool AETCompare(vector<int>& shellnum1, vector<int>& shellnum2);
bool PRIMCHANGEPRO(const xstructure& a);
// [OBSOLETE] string AFLOW_BinaryRead(void);
// [OBSOLETE] string AFLOW_Binary_Angle_Read(void);
void ProtoClassification_ANGLE(void);
bool LatticeCompare(xmatrix<double>& a, xmatrix<double>& b);

class AETS{
public:
  AETS(void);
  AETS(const vector<double>& b);
  AETS(const AETS& b);
  //Parameter                                                                                                                                                                                                        
  double RelDis[7];
  double AtomNum[7];
  double MaxG;
  int MaxS;
  int NumPoly;
  int totalnum;
  //Function                                                                                                                                                                                                        
  void GetMaxGap(void);
  void GetMaxShell(void);
};

// PRIMITIVETRANS
void Primitive_Transform(xstructure& a, bool& Tflag);
void Primitive_Transform(void);
// MAGNET
void GetMagneticParameters(vector<string>& argv,ostream& oss);
void GetMagPara_heusler(vector<string>& argv, ostream& oss);
// PHASEDIAGRAM
void INPUTDATAFORTERPHASE(vector<string>& argv);
void GENERATESTABLELIST(vector<string>& argv);

#endif
