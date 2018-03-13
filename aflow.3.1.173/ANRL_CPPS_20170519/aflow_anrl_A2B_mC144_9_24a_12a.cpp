// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2016
// FILE "ANRL/aflow_anrl_A2B_mC144_9_24a_12a.cpp"

#ifndef _AFLOW_ANRL_A2B_mC144_9_24a_12a_CPP
#define _AFLOW_ANRL_A2B_mC144_9_24a_12a_CPP
#include "../aflow.h"

namespace anrl {
  uint WebANRL_A2B_mC144_9_24a_12a(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A2B_mC144_9_24a_12a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG) {
    // system A2B_mC144_9_24a_12a

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A2B_mC144_9_24a_12a(web,LDEBUG); // PLUG WEB STUFF
      #ifdef _ANRL_NOWEB_
      cout << "no web" << endl;
      #else
      cout << web.str() << endl;
      #endif
      exit(0);
    }

    vector<double> vparameters;
    aurostd::string2tokens(parameters,vparameters,",");

    uint nspecies,natoms,spacegroup,nunderscores,nparameters;
    string label,Pearson_symbol,params,Strukturbericht,prototype,dialect;

    anrl::vproto2tokens(proto_line,label,nspecies,natoms,spacegroup,nunderscores,nparameters,Pearson_symbol,params,Strukturbericht,prototype,dialect);

    if(!anrl::PrototypeANRL_Consistency(oss,vparameters.size(),nparameters,prototype,label,
                 Strukturbericht,Pearson_symbol,spacegroup, params)) { exit(0);}    

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: vparameters.size()=" << vparameters.size() << endl;}

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    _atom atom;

    uint i=0;
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: c=" << c << " (c/a=" << covera << ")" << endl;}
    double beta=vparameters.at(i++);               if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: beta=" << beta << endl;}
    
    double x1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x1=" << x1 << endl;}
    double y1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y1=" << y1 << endl;}
    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z1=" << z1 << endl;}
    double x2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x2=" << x2 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y2=" << y2 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z2=" << z2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x3=" << x3 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y3=" << y3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x11=" << x11 << endl;}
    double y11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y11=" << y11 << endl;}
    double z11=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x12=" << x12 << endl;}
    double y12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y12=" << y12 << endl;}
    double z12=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z12=" << z12 << endl;}
    double x13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x13=" << x13 << endl;}
    double y13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y13=" << y13 << endl;}
    double z13=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z13=" << z13 << endl;}
    double x14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x14=" << x14 << endl;}
    double y14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y14=" << y14 << endl;}
    double z14=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z14=" << z14 << endl;}
    double x15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x15=" << x15 << endl;}
    double y15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y15=" << y15 << endl;}
    double z15=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z15=" << z15 << endl;}
    double x16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x16=" << x16 << endl;}
    double y16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y16=" << y16 << endl;}
    double z16=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z16=" << z16 << endl;}
    double x17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x17=" << x17 << endl;}
    double y17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y17=" << y17 << endl;}
    double z17=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z17=" << z17 << endl;}
    double x18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x18=" << x18 << endl;}
    double y18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y18=" << y18 << endl;}
    double z18=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z18=" << z18 << endl;}
    double x19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x19=" << x19 << endl;}
    double y19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y19=" << y19 << endl;}
    double z19=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z19=" << z19 << endl;}
    double x20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x20=" << x20 << endl;}
    double y20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y20=" << y20 << endl;}
    double z20=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z20=" << z20 << endl;}
    double x21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x21=" << x21 << endl;}
    double y21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y21=" << y21 << endl;}
    double z21=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z21=" << z21 << endl;}
    double x22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x22=" << x22 << endl;}
    double y22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y22=" << y22 << endl;}
    double z22=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z22=" << z22 << endl;}
    double x23=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x23=" << x23 << endl;}
    double y23=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y23=" << y23 << endl;}
    double z23=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z23=" << z23 << endl;}
    double x24=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x24=" << x24 << endl;}
    double y24=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y24=" << y24 << endl;}
    double z24=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z24=" << z24 << endl;}
    double x25=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x25=" << x25 << endl;}
    double y25=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y25=" << y25 << endl;}
    double z25=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z25=" << z25 << endl;}
    double x26=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x26=" << x26 << endl;}
    double y26=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y26=" << y26 << endl;}
    double z26=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z26=" << z26 << endl;}
    double x27=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x27=" << x27 << endl;}
    double y27=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y27=" << y27 << endl;}
    double z27=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z27=" << z27 << endl;}
    double x28=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x28=" << x28 << endl;}
    double y28=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y28=" << y28 << endl;}
    double z28=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z28=" << z28 << endl;}
    double x29=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x29=" << x29 << endl;}
    double y29=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y29=" << y29 << endl;}
    double z29=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z29=" << z29 << endl;}
    double x30=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x30=" << x30 << endl;}
    double y30=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y30=" << y30 << endl;}
    double z30=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z30=" << z30 << endl;}
    double x31=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x31=" << x31 << endl;}
    double y31=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y31=" << y31 << endl;}
    double z31=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z31=" << z31 << endl;}
    double x32=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x32=" << x32 << endl;}
    double y32=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y32=" << y32 << endl;}
    double z32=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z32=" << z32 << endl;}
    double x33=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x33=" << x33 << endl;}
    double y33=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y33=" << y33 << endl;}
    double z33=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z33=" << z33 << endl;}
    double x34=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x34=" << x34 << endl;}
    double y34=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y34=" << y34 << endl;}
    double z34=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z34=" << z34 << endl;}
    double x35=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x35=" << x35 << endl;}
    double y35=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y35=" << y35 << endl;}
    double z35=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z35=" << z35 << endl;}
    double x36=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: x36=" << x36 << endl;}
    double y36=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: y36=" << y36 << endl;}
    double z36=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: z36=" << z36 << endl;}
        
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: cos(beta)=" << cos(deg2rad*beta)  << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A2B_mC144_9_24a_12a: sin(beta)=" << sin(deg2rad*beta)  << endl;}
        
    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG#="+aurostd::utype2string(spacegroup)+DOI_ANRL;
    str.scale=1.0;

    a1=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn;
    a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
    a3=c*cos(deg2rad*beta)*xn+c*sin(deg2rad*beta)*zn;
    
    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    str.FixLattices(); // Reciprocal/f2c/c2f
    
    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=(x1-y1);atom.fpos(2)=(x1+y1);atom.fpos(3)=z1;        // atom B1
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1
    
    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=(x1+y1);atom.fpos(2)=(x1-y1);atom.fpos(3)=((1.0/2.0)+z1);// atom B2
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2
    
    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=(x2-y2);atom.fpos(2)=(x2+y2);atom.fpos(3)=z2;        // atom B3
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3
    
    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=(x2+y2);atom.fpos(2)=(x2-y2);atom.fpos(3)=((1.0/2.0)+z2);// atom B4
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4
    
    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=(x3-y3);atom.fpos(2)=(x3+y3);atom.fpos(3)=z3;        // atom B5
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5
    
    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=(x3+y3);atom.fpos(2)=(x3-y3);atom.fpos(3)=((1.0/2.0)+z3);// atom B6
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6
    
    atom.name="A"; atom.type=0;                                       // atom B7
    atom.fpos(1)=(x4-y4);atom.fpos(2)=(x4+y4);atom.fpos(3)=z4;        // atom B7
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7
    
    atom.name="A"; atom.type=0;                                       // atom B8
    atom.fpos(1)=(x4+y4);atom.fpos(2)=(x4-y4);atom.fpos(3)=((1.0/2.0)+z4);// atom B8
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8
    
    atom.name="A"; atom.type=0;                                       // atom B9
    atom.fpos(1)=(x5-y5);atom.fpos(2)=(x5+y5);atom.fpos(3)=z5;        // atom B9
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9
    
    atom.name="A"; atom.type=0;                                       // atom B10
    atom.fpos(1)=(x5+y5);atom.fpos(2)=(x5-y5);atom.fpos(3)=((1.0/2.0)+z5);// atom B10
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10
    
    atom.name="A"; atom.type=0;                                       // atom B11
    atom.fpos(1)=(x6-y6);atom.fpos(2)=(x6+y6);atom.fpos(3)=z6;        // atom B11
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11
    
    atom.name="A"; atom.type=0;                                       // atom B12
    atom.fpos(1)=(x6+y6);atom.fpos(2)=(x6-y6);atom.fpos(3)=((1.0/2.0)+z6);// atom B12
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12
    
    atom.name="A"; atom.type=0;                                       // atom B13
    atom.fpos(1)=(x7-y7);atom.fpos(2)=(x7+y7);atom.fpos(3)=z7;        // atom B13
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13
    
    atom.name="A"; atom.type=0;                                       // atom B14
    atom.fpos(1)=(x7+y7);atom.fpos(2)=(x7-y7);atom.fpos(3)=((1.0/2.0)+z7);// atom B14
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14
    
    atom.name="A"; atom.type=0;                                       // atom B15
    atom.fpos(1)=(x8-y8);atom.fpos(2)=(x8+y8);atom.fpos(3)=z8;        // atom B15
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15
    
    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=(x8+y8);atom.fpos(2)=(x8-y8);atom.fpos(3)=((1.0/2.0)+z8);// atom B16
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16
    
    atom.name="A"; atom.type=0;                                       // atom B17
    atom.fpos(1)=(x9-y9);atom.fpos(2)=(x9+y9);atom.fpos(3)=z9;        // atom B17
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17
    
    atom.name="A"; atom.type=0;                                       // atom B18
    atom.fpos(1)=(x9+y9);atom.fpos(2)=(x9-y9);atom.fpos(3)=((1.0/2.0)+z9);// atom B18
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18
    
    atom.name="A"; atom.type=0;                                       // atom B19
    atom.fpos(1)=(x10-y10);atom.fpos(2)=(x10+y10);atom.fpos(3)=z10;   // atom B19
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19
    
    atom.name="A"; atom.type=0;                                       // atom B20
    atom.fpos(1)=(x10+y10);atom.fpos(2)=(x10-y10);atom.fpos(3)=((1.0/2.0)+z10);// atom B20
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20
    
    atom.name="A"; atom.type=0;                                       // atom B21
    atom.fpos(1)=(x11-y11);atom.fpos(2)=(x11+y11);atom.fpos(3)=z11;   // atom B21
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21
    
    atom.name="A"; atom.type=0;                                       // atom B22
    atom.fpos(1)=(x11+y11);atom.fpos(2)=(x11-y11);atom.fpos(3)=((1.0/2.0)+z11);// atom B22
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22
    
    atom.name="A"; atom.type=0;                                       // atom B23
    atom.fpos(1)=(x12-y12);atom.fpos(2)=(x12+y12);atom.fpos(3)=z12;   // atom B23
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23
    
    atom.name="A"; atom.type=0;                                       // atom B24
    atom.fpos(1)=(x12+y12);atom.fpos(2)=(x12-y12);atom.fpos(3)=((1.0/2.0)+z12);// atom B24
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24
    
    atom.name="A"; atom.type=0;                                       // atom B25
    atom.fpos(1)=(x13-y13);atom.fpos(2)=(x13+y13);atom.fpos(3)=z13;   // atom B25
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25
    
    atom.name="A"; atom.type=0;                                       // atom B26
    atom.fpos(1)=(x13+y13);atom.fpos(2)=(x13-y13);atom.fpos(3)=((1.0/2.0)+z13);// atom B26
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26
    
    atom.name="A"; atom.type=0;                                       // atom B27
    atom.fpos(1)=(x14-y14);atom.fpos(2)=(x14+y14);atom.fpos(3)=z14;   // atom B27
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27
    
    atom.name="A"; atom.type=0;                                       // atom B28
    atom.fpos(1)=(x14+y14);atom.fpos(2)=(x14-y14);atom.fpos(3)=((1.0/2.0)+z14);// atom B28
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28
    
    atom.name="A"; atom.type=0;                                       // atom B29
    atom.fpos(1)=(x15-y15);atom.fpos(2)=(x15+y15);atom.fpos(3)=z15;   // atom B29
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29
    
    atom.name="A"; atom.type=0;                                       // atom B30
    atom.fpos(1)=(x15+y15);atom.fpos(2)=(x15-y15);atom.fpos(3)=((1.0/2.0)+z15);// atom B30
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30
    
    atom.name="A"; atom.type=0;                                       // atom B31
    atom.fpos(1)=(x16-y16);atom.fpos(2)=(x16+y16);atom.fpos(3)=z16;   // atom B31
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31
    
    atom.name="A"; atom.type=0;                                       // atom B32
    atom.fpos(1)=(x16+y16);atom.fpos(2)=(x16-y16);atom.fpos(3)=((1.0/2.0)+z16);// atom B32
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32
    
    atom.name="A"; atom.type=0;                                       // atom B33
    atom.fpos(1)=(x17-y17);atom.fpos(2)=(x17+y17);atom.fpos(3)=z17;   // atom B33
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33
    
    atom.name="A"; atom.type=0;                                       // atom B34
    atom.fpos(1)=(x17+y17);atom.fpos(2)=(x17-y17);atom.fpos(3)=((1.0/2.0)+z17);// atom B34
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34
    
    atom.name="A"; atom.type=0;                                       // atom B35
    atom.fpos(1)=(x18-y18);atom.fpos(2)=(x18+y18);atom.fpos(3)=z18;   // atom B35
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35
    
    atom.name="A"; atom.type=0;                                       // atom B36
    atom.fpos(1)=(x18+y18);atom.fpos(2)=(x18-y18);atom.fpos(3)=((1.0/2.0)+z18);// atom B36
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36
    
    atom.name="A"; atom.type=0;                                       // atom B37
    atom.fpos(1)=(x19-y19);atom.fpos(2)=(x19+y19);atom.fpos(3)=z19;   // atom B37
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B37 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B37
    
    atom.name="A"; atom.type=0;                                       // atom B38
    atom.fpos(1)=(x19+y19);atom.fpos(2)=(x19-y19);atom.fpos(3)=((1.0/2.0)+z19);// atom B38
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B38 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B38
    
    atom.name="A"; atom.type=0;                                       // atom B39
    atom.fpos(1)=(x20-y20);atom.fpos(2)=(x20+y20);atom.fpos(3)=z20;   // atom B39
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B39 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B39
    
    atom.name="A"; atom.type=0;                                       // atom B40
    atom.fpos(1)=(x20+y20);atom.fpos(2)=(x20-y20);atom.fpos(3)=((1.0/2.0)+z20);// atom B40
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B40 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B40
    
    atom.name="A"; atom.type=0;                                       // atom B41
    atom.fpos(1)=(x21-y21);atom.fpos(2)=(x21+y21);atom.fpos(3)=z21;   // atom B41
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B41 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B41
    
    atom.name="A"; atom.type=0;                                       // atom B42
    atom.fpos(1)=(x21+y21);atom.fpos(2)=(x21-y21);atom.fpos(3)=((1.0/2.0)+z21);// atom B42
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B42 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B42
    
    atom.name="A"; atom.type=0;                                       // atom B43
    atom.fpos(1)=(x22-y22);atom.fpos(2)=(x22+y22);atom.fpos(3)=z22;   // atom B43
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B43 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B43
    
    atom.name="A"; atom.type=0;                                       // atom B44
    atom.fpos(1)=(x22+y22);atom.fpos(2)=(x22-y22);atom.fpos(3)=((1.0/2.0)+z22);// atom B44
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B44 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B44
    
    atom.name="A"; atom.type=0;                                       // atom B45
    atom.fpos(1)=(x23-y23);atom.fpos(2)=(x23+y23);atom.fpos(3)=z23;   // atom B45
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B45 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B45
    
    atom.name="A"; atom.type=0;                                       // atom B46
    atom.fpos(1)=(x23+y23);atom.fpos(2)=(x23-y23);atom.fpos(3)=((1.0/2.0)+z23);// atom B46
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B46 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B46
    
    atom.name="A"; atom.type=0;                                       // atom B47
    atom.fpos(1)=(x24-y24);atom.fpos(2)=(x24+y24);atom.fpos(3)=z24;   // atom B47
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B47 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B47
    
    atom.name="A"; atom.type=0;                                       // atom B48
    atom.fpos(1)=(x24+y24);atom.fpos(2)=(x24-y24);atom.fpos(3)=((1.0/2.0)+z24);// atom B48
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B48 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B48
    
    atom.name="B"; atom.type=1;                                       // atom B49
    atom.fpos(1)=(x25-y25);atom.fpos(2)=(x25+y25);atom.fpos(3)=z25;   // atom B49
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B49 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B49
    
    atom.name="B"; atom.type=1;                                       // atom B50
    atom.fpos(1)=(x25+y25);atom.fpos(2)=(x25-y25);atom.fpos(3)=((1.0/2.0)+z25);// atom B50
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B50 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B50
    
    atom.name="B"; atom.type=1;                                       // atom B51
    atom.fpos(1)=(x26-y26);atom.fpos(2)=(x26+y26);atom.fpos(3)=z26;   // atom B51
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B51 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B51
    
    atom.name="B"; atom.type=1;                                       // atom B52
    atom.fpos(1)=(x26+y26);atom.fpos(2)=(x26-y26);atom.fpos(3)=((1.0/2.0)+z26);// atom B52
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B52 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B52
    
    atom.name="B"; atom.type=1;                                       // atom B53
    atom.fpos(1)=(x27-y27);atom.fpos(2)=(x27+y27);atom.fpos(3)=z27;   // atom B53
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B53 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B53
    
    atom.name="B"; atom.type=1;                                       // atom B54
    atom.fpos(1)=(x27+y27);atom.fpos(2)=(x27-y27);atom.fpos(3)=((1.0/2.0)+z27);// atom B54
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B54 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B54
    
    atom.name="B"; atom.type=1;                                       // atom B55
    atom.fpos(1)=(x28-y28);atom.fpos(2)=(x28+y28);atom.fpos(3)=z28;   // atom B55
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B55 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B55
    
    atom.name="B"; atom.type=1;                                       // atom B56
    atom.fpos(1)=(x28+y28);atom.fpos(2)=(x28-y28);atom.fpos(3)=((1.0/2.0)+z28);// atom B56
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B56 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B56
    
    atom.name="B"; atom.type=1;                                       // atom B57
    atom.fpos(1)=(x29-y29);atom.fpos(2)=(x29+y29);atom.fpos(3)=z29;   // atom B57
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B57 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B57
    
    atom.name="B"; atom.type=1;                                       // atom B58
    atom.fpos(1)=(x29+y29);atom.fpos(2)=(x29-y29);atom.fpos(3)=((1.0/2.0)+z29);// atom B58
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B58 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B58
    
    atom.name="B"; atom.type=1;                                       // atom B59
    atom.fpos(1)=(x30-y30);atom.fpos(2)=(x30+y30);atom.fpos(3)=z30;   // atom B59
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B59 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B59
    
    atom.name="B"; atom.type=1;                                       // atom B60
    atom.fpos(1)=(x30+y30);atom.fpos(2)=(x30-y30);atom.fpos(3)=((1.0/2.0)+z30);// atom B60
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B60 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B60
    
    atom.name="B"; atom.type=1;                                       // atom B61
    atom.fpos(1)=(x31-y31);atom.fpos(2)=(x31+y31);atom.fpos(3)=z31;   // atom B61
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B61 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B61
    
    atom.name="B"; atom.type=1;                                       // atom B62
    atom.fpos(1)=(x31+y31);atom.fpos(2)=(x31-y31);atom.fpos(3)=((1.0/2.0)+z31);// atom B62
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B62 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B62
    
    atom.name="B"; atom.type=1;                                       // atom B63
    atom.fpos(1)=(x32-y32);atom.fpos(2)=(x32+y32);atom.fpos(3)=z32;   // atom B63
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B63 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B63
    
    atom.name="B"; atom.type=1;                                       // atom B64
    atom.fpos(1)=(x32+y32);atom.fpos(2)=(x32-y32);atom.fpos(3)=((1.0/2.0)+z32);// atom B64
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B64 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B64
    
    atom.name="B"; atom.type=1;                                       // atom B65
    atom.fpos(1)=(x33-y33);atom.fpos(2)=(x33+y33);atom.fpos(3)=z33;   // atom B65
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B65 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B65
    
    atom.name="B"; atom.type=1;                                       // atom B66
    atom.fpos(1)=(x33+y33);atom.fpos(2)=(x33-y33);atom.fpos(3)=((1.0/2.0)+z33);// atom B66
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B66 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B66
    
    atom.name="B"; atom.type=1;                                       // atom B67
    atom.fpos(1)=(x34-y34);atom.fpos(2)=(x34+y34);atom.fpos(3)=z34;   // atom B67
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B67 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B67
    
    atom.name="B"; atom.type=1;                                       // atom B68
    atom.fpos(1)=(x34+y34);atom.fpos(2)=(x34-y34);atom.fpos(3)=((1.0/2.0)+z34);// atom B68
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B68 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B68
    
    atom.name="B"; atom.type=1;                                       // atom B69
    atom.fpos(1)=(x35-y35);atom.fpos(2)=(x35+y35);atom.fpos(3)=z35;   // atom B69
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B69 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B69
    
    atom.name="B"; atom.type=1;                                       // atom B70
    atom.fpos(1)=(x35+y35);atom.fpos(2)=(x35-y35);atom.fpos(3)=((1.0/2.0)+z35);// atom B70
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B70 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B70
    
    atom.name="B"; atom.type=1;                                       // atom B71
    atom.fpos(1)=(x36-y36);atom.fpos(2)=(x36+y36);atom.fpos(3)=z36;   // atom B71
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B71 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B71
    
    atom.name="B"; atom.type=1;                                       // atom B72
    atom.fpos(1)=(x36+y36);atom.fpos(2)=(x36-y36);atom.fpos(3)=((1.0/2.0)+z36);// atom B72
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B72 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B72
    

    return str.atoms.size();  
  }
} // namespace anrl::

namespace anrl {
  uint WebANRL_A2B_mC144_9_24a_12a(stringstream& web,bool LDEBUG) {
    #ifndef _ANRL_NOWEB_
    web << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << endl;
    web << "<head profile=\"http://www.w3.org/1999/xhtml/vocab\">" << endl;
    web << "<title>" << endl;
    web << "" << endl;
    web << "</title>" << endl;
    web << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\" />" << endl;
    web << "<script type=\"text/x-mathjax-config\">" << endl;
    web << "MathJax.Hub.Config({" << endl;
    web << "TeX: {" << endl;
    web << "noErrors: {disabled: true}, // Show error messages" << endl;
    web << "MAXBUFFER: 50*1024, // Set size of buffer in bytes" << endl;
    web << "}," << endl;
    web << "tex2jax: {" << endl;
    web << "inlineMath: [ [\'$\',\'$\'], [\'\\\\(\',\'\\\\)\'] ]" << endl;
    web << "}," << endl;
    web << "});" << endl;
    web << "</script>" << endl;
    web << "<script type=\"text/javascript\"src=\"https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script>" << endl;
    web << "<link rel=\"stylesheet\" href=\"css/bootstrap-grid.min.css\" />" << endl;
    web << "<link rel=\"stylesheet\" href=\"css/style.css\" />" << endl;
    web << "<meta charset=\"utf-8\" content=\"charset\" />" << endl;
    web << "</head>" << endl;
    web << "<body>" << endl;
    web << "<nav class=\"nav\">" << endl;
    web << "<div class=\"container\">" << endl;
    web << "<ul class=\"nav-links\">" << endl;
    web << "<li class=\"dropdown\">" << endl;
    web << "<div class=\"dropbtn\">Lattice Info</div>" << endl;
    web << "<div class=\"dropdown-content\">" << endl;
    web << "<a href=\"triclinic_lattice.html\">Triclinic</a>" << endl;
    web << "<a href=\"monoclinic_lattice.html\">Monoclinic</a>" << endl;
    web << "<a href=\"orthorhombic_lattice.html\">Orthorhombic</a>" << endl;
    web << "<a href=\"tetragonal_lattice.html\">Tetragonal</a>" << endl;
    web << "<a href=\"trigonal_lattice.html\">Trigonal</a>" << endl;
    web << "<a href=\"hexagonal_lattice.html\">Hexagonal</a>" << endl;
    web << "<a href=\"cubic_lattice.html\">Cubic</a>" << endl;
    web << "</div>" << endl;
    web << "</li>" << endl;
    web << "<li> <a href=\"./chemical_symbols.html\">Chemical Symbol</a> </li>" << endl;
    web << "<li> <a href=\"./strukturberichts.html\">Strukurbericht</a> </li>" << endl;
    web << "<li> <a href=\"./pearson_symbols.html\">Pearson Symbol</a> </li>" << endl;
    web << "<li> <a href=\"./space_groups.html\">Space Group</a> </li>" << endl;
    web << "</ul>" << endl;
    web << "</div>" << endl;
    web << "</nav>" << endl;
    web << "<div class=\"header\">" << endl;
    web << "<div class=\"container header_container\">" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div id=\"logo\" class=\"col-md-4\">" << endl;
    web << "<a href=\"./index.html\">" << endl;
   web << "<svg class=\"logo\" xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 189.08 253.13\"><defs><style>.cls-1{fill:#4e8150;}.cls-2{fill:#659b65;}.cls-3,.cls-6{fill:#f3f2f2;}.cls-4{fill:#eaeaea;}.cls-5{fill:#ffba54;}.cls-6{stroke:#efefef;}.cls-10,.cls-13,.cls-15,.cls-16,.cls-18,.cls-19,.cls-20,.cls-21,.cls-6,.cls-7,.cls-9{stroke-miterlimit:10;}.cls-13,.cls-6{stroke-width:0.5px;}.cls-13,.cls-15,.cls-19,.cls-20,.cls-21,.cls-7{fill:none;}.cls-7{stroke:#000;}.cls-19,.cls-20,.cls-21,.cls-7{stroke-width:0.25px;}.cls-8{fill:#db3b19;}.cls-10,.cls-9{fill:#3edbc0;}.cls-18,.cls-9{stroke:#aedfac;}.cls-10{stroke:#db3b19;}.cls-11{fill:#667ba3;}.cls-12{fill:#9fbdbe;}.cls-13{stroke:#fff;}.cls-14,.cls-18{fill:#aedfac;}.cls-15{stroke:#92b58f;}.cls-16,.cls-24{fill:#55a4e5;}.cls-16{stroke:#c75c14;}.cls-17{fill:#92b58f;}.cls-19{stroke:#55a4e5;}.cls-20{stroke:#dc4a3d;}.cls-21{stroke:#4ece8e;}.cls-22,.cls-23,.cls-24{font-size:5px;font-family:MyriadPro-Regular, Myriad Pro;}.cls-22{fill:#dc4a3d;}.cls-23{fill:#4ece8e;}.cls-25{fill:#fff;}</style></defs><title>encyclopedia_teal_natural</title><g id=\"book\"><path class=\"cls-1\" d=\"M20,0H189.08a0,0,0,0,1,0,0V237.18a0,0,0,0,1,0,0H25a25,25,0,0,1-25-25V20A20,20,0,0,1,20,0Z\"/><rect class=\"cls-2\" x=\"24.16\" width=\"164.92\" height=\"183.88\"/><path class=\"cls-3\" d=\"M25.5,190.76H189.08a0,0,0,0,1,0,0v40a0,0,0,0,1,0,0H26.2a20,20,0,0,1-20-20v-.7a19.3,19.3,0,0,1,19.3-19.3Z\"/><path class=\"cls-4\" d=\"M6.2,210.51H189.08a0,0,0,0,1,0,0v20.25a0,0,0,0,1,0,0H26.2a20,20,0,0,1-20-20v-.25a0,0,0,0,1,0,0Z\"/><path class=\"cls-5\" d=\"M140.2,210.51v42.62l14.5-15.84L169,252.85V210.51Z\"/></g><g id=\"background\"><rect class=\"cls-6\" x=\"33.91\" y=\"59.71\" width=\"146.95\" height=\"115.01\"/></g><g id=\"atoms\"><line class=\"cls-7\" x1=\"118.49\" y1=\"148.26\" x2=\"74.33\" y2=\"137.33\"/><line class=\"cls-7\" x1=\"74.21\" y1=\"87.2\" x2=\"74.33\" y2=\"137.33\"/><circle class=\"cls-8\" cx=\"74.33\" cy=\"137.33\" r=\"2.92\" transform=\"translate(-75.33 92.78) rotate(-45)\"/><line class=\"cls-9\" x1=\"73.16\" y1=\"126.25\" x2=\"82.87\" y2=\"119.42\"/><line class=\"cls-10\" x1=\"63.85\" y1=\"132.8\" x2=\"73.16\" y2=\"126.25\"/><polygon class=\"cls-11\" points=\"42.33 152.39 58.27 119.7 64.4 132.51 42.33 152.39\"/><polyline class=\"cls-11\" points=\"80.46 150.2 64.4 132.51 58.27 119.7\"/><polyline class=\"cls-12\" points=\"80.46 150.2 42.33 152.39 64.4 132.51\"/><polyline class=\"cls-12\" points=\"96.27 117.64 58.27 119.7 80.33 99.95\"/><polyline class=\"cls-11\" points=\"96.27 117.64 80.33 99.95 74.21 87.2\"/><polygon class=\"cls-11\" points=\"58.27 119.7 74.21 87.2 80.33 99.95 58.27 119.7\"/><line class=\"cls-13\" x1=\"74.21\" y1=\"87.2\" x2=\"80.33\" y2=\"99.95\"/><line class=\"cls-13\" x1=\"80.33\" y1=\"99.95\" x2=\"58.27\" y2=\"119.7\"/><line class=\"cls-9\" x1=\"81.81\" y1=\"109.83\" x2=\"83.38\" y2=\"118.95\"/><line class=\"cls-10\" x1=\"80.33\" y1=\"99.95\" x2=\"81.81\" y2=\"109.83\"/><line class=\"cls-13\" x1=\"96.27\" y1=\"117.64\" x2=\"80.33\" y2=\"99.95\"/><line class=\"cls-7\" x1=\"42.33\" y1=\"102.26\" x2=\"86.65\" y2=\"112.76\"/><line class=\"cls-7\" x1=\"74.21\" y1=\"87.2\" x2=\"42.33\" y2=\"102.26\"/><line class=\"cls-7\" x1=\"118.49\" y1=\"97.89\" x2=\"74.21\" y2=\"87.2\"/><line class=\"cls-7\" x1=\"86.65\" y1=\"112.76\" x2=\"118.49\" y2=\"97.89\"/><line class=\"cls-7\" x1=\"118.49\" y1=\"148.26\" x2=\"118.49\" y2=\"97.89\"/><line class=\"cls-7\" x1=\"42.33\" y1=\"102.26\" x2=\"42.33\" y2=\"152.39\"/><line class=\"cls-7\" x1=\"86.65\" y1=\"163.01\" x2=\"42.33\" y2=\"152.39\"/><line class=\"cls-7\" x1=\"118.49\" y1=\"148.26\" x2=\"86.65\" y2=\"163.01\"/><polygon class=\"cls-11\" points=\"80.33 150.14 96.27 117.64 102.58 130.45 80.33 150.14\"/><line class=\"cls-9\" x1=\"93.19\" y1=\"124.83\" x2=\"83.38\" y2=\"118.95\"/><line class=\"cls-10\" x1=\"102.58\" y1=\"130.45\" x2=\"93.19\" y2=\"124.83\"/><polyline class=\"cls-12\" points=\"118.49 148.26 80.33 150.14 102.58 130.45\"/><polyline class=\"cls-11\" points=\"118.49 148.26 102.58 130.45 96.27 117.64\"/><circle class=\"cls-8\" cx=\"42.33\" cy=\"102.26\" r=\"2.92\" transform=\"translate(-59.91 59.89) rotate(-45)\"/><circle class=\"cls-8\" cx=\"74.21\" cy=\"87.2\" r=\"2.92\" transform=\"translate(-39.93 78.01) rotate(-45)\"/><circle class=\"cls-8\" cx=\"118.49\" cy=\"97.89\" r=\"2.92\" transform=\"translate(-34.51 112.45) rotate(-45)\"/><line class=\"cls-13\" x1=\"58.27\" y1=\"119.7\" x2=\"64.4\" y2=\"132.51\"/><line class=\"cls-13\" x1=\"42.33\" y1=\"152.39\" x2=\"64.4\" y2=\"132.51\"/><line class=\"cls-13\" x1=\"80.46\" y1=\"150.2\" x2=\"64.4\" y2=\"132.51\"/><line class=\"cls-13\" x1=\"96.27\" y1=\"117.64\" x2=\"102.58\" y2=\"130.45\"/><line class=\"cls-13\" x1=\"118.49\" y1=\"148.26\" x2=\"102.58\" y2=\"130.45\"/><line class=\"cls-13\" x1=\"80.46\" y1=\"150.2\" x2=\"102.58\" y2=\"130.45\"/><circle class=\"cls-8\" cx=\"58.27\" cy=\"119.7\" r=\"2.92\" transform=\"translate(-67.57 76.26) rotate(-45)\"/><circle class=\"cls-8\" cx=\"96.27\" cy=\"117.64\" r=\"2.92\" transform=\"translate(-54.99 102.53) rotate(-45)\"/><circle class=\"cls-8\" cx=\"80.46\" cy=\"150.2\" r=\"2.92\" transform=\"translate(-82.64 100.89) rotate(-45)\"/><circle class=\"cls-8\" cx=\"118.49\" cy=\"148.01\" r=\"2.92\" transform=\"translate(-69.96 127.14) rotate(-45)\"/><circle class=\"cls-8\" cx=\"42.33\" cy=\"152.39\" r=\"2.92\" transform=\"translate(-95.36 74.57) rotate(-45)\"/><circle class=\"cls-14\" cx=\"83.38\" cy=\"118.95\" r=\"2.92\" transform=\"translate(-45.78 185.43) rotate(-82.21)\"/><circle class=\"cls-8\" cx=\"80.33\" cy=\"99.95\" r=\"2.92\" transform=\"translate(-47.15 86.08) rotate(-45)\"/><circle class=\"cls-8\" cx=\"64.4\" cy=\"132.51\" r=\"2.92\" transform=\"translate(-74.84 84.35) rotate(-45)\"/><circle class=\"cls-8\" cx=\"102.58\" cy=\"130.45\" r=\"2.92\" transform=\"translate(-62.2 110.75) rotate(-45)\"/><line class=\"cls-7\" x1=\"86.65\" y1=\"112.76\" x2=\"86.65\" y2=\"163.01\"/><circle class=\"cls-8\" cx=\"86.65\" cy=\"163.01\" r=\"2.92\" transform=\"translate(-89.89 109.01) rotate(-45)\"/><line class=\"cls-15\" x1=\"83.38\" y1=\"118.95\" x2=\"84.95\" y2=\"115.98\"/><line class=\"cls-16\" x1=\"84.95\" y1=\"115.98\" x2=\"86.63\" y2=\"112.86\"/><circle class=\"cls-17\" cx=\"83.4\" cy=\"118.91\" r=\"0.5\" transform=\"translate(-60.73 136.82) rotate(-62.04)\"/><circle class=\"cls-8\" cx=\"86.65\" cy=\"112.76\" r=\"2.92\" transform=\"translate(-54.36 94.3) rotate(-45)\"/><line class=\"cls-18\" x1=\"143.6\" y1=\"95.3\" x2=\"153.32\" y2=\"88.47\"/><line class=\"cls-10\" x1=\"134.3\" y1=\"101.85\" x2=\"143.6\" y2=\"95.3\"/><line class=\"cls-18\" x1=\"152.26\" y1=\"78.88\" x2=\"153.82\" y2=\"88\"/><line class=\"cls-10\" x1=\"150.78\" y1=\"69\" x2=\"152.26\" y2=\"78.88\"/><line class=\"cls-18\" x1=\"163.64\" y1=\"93.88\" x2=\"153.82\" y2=\"88\"/><line class=\"cls-10\" x1=\"173.03\" y1=\"99.5\" x2=\"163.64\" y2=\"93.88\"/><circle class=\"cls-14\" cx=\"153.82\" cy=\"88\" r=\"2.92\" transform=\"translate(45.78 228.47) rotate(-82.21)\"/><line class=\"cls-15\" x1=\"153.82\" y1=\"88\" x2=\"155.4\" y2=\"85.03\"/><line class=\"cls-16\" x1=\"155.4\" y1=\"85.03\" x2=\"157.08\" y2=\"81.91\"/><circle class=\"cls-17\" cx=\"153.85\" cy=\"87.96\" r=\"0.5\" transform=\"translate(4.02 182.6) rotate(-62.04)\"/><line class=\"cls-7\" x1=\"150.78\" y1=\"69\" x2=\"134.84\" y2=\"101.56\"/><line class=\"cls-7\" x1=\"173.03\" y1=\"99.5\" x2=\"134.84\" y2=\"101.56\"/><line class=\"cls-7\" x1=\"150.78\" y1=\"69\" x2=\"173.03\" y2=\"99.5\"/><line class=\"cls-7\" x1=\"157.08\" y1=\"81.91\" x2=\"173.03\" y2=\"99.5\"/><line class=\"cls-7\" x1=\"150.78\" y1=\"69\" x2=\"157.08\" y2=\"81.91\"/><circle class=\"cls-8\" cx=\"157.09\" cy=\"81.81\" r=\"2.92\" transform=\"translate(-11.84 135.05) rotate(-45)\"/><circle class=\"cls-8\" cx=\"150.78\" cy=\"69\" r=\"2.92\" transform=\"translate(-4.63 126.83) rotate(-45)\"/><circle class=\"cls-8\" cx=\"173.03\" cy=\"99.5\" r=\"2.92\" transform=\"translate(-19.68 151.5) rotate(-45)\"/><circle class=\"cls-8\" cx=\"134.84\" cy=\"101.56\" r=\"2.92\" transform=\"translate(-32.32 125.1) rotate(-45)\"/></g><g id=\"axis\"><line class=\"cls-19\" x1=\"167.74\" y1=\"157.98\" x2=\"150.97\" y2=\"153.83\"/><line class=\"cls-20\" x1=\"150.98\" y1=\"136.52\" x2=\"150.98\" y2=\"153.97\"/><line class=\"cls-21\" x1=\"138.82\" y1=\"159.55\" x2=\"150.97\" y2=\"153.83\"/><text class=\"cls-22\" transform=\"translate(148.42 136.26)\">a1</text><text class=\"cls-23\" transform=\"translate(133.84 162.4)\">a2</text><text class=\"cls-24\" transform=\"translate(168.7 160.56)\">a3</text></g><g id=\"AFLOW\"><path class=\"cls-25\" d=\"M55.26,42.28V37.45H38.44v4.83H33.91V25.38c0-4.24,2.81-7.24,6.79-7.24H53c3.91,0,6.79,3.07,6.79,7.24v16.9Zm0-16.9A2.15,2.15,0,0,0,53,23H40.7a2.12,2.12,0,0,0-2.26,2.41v7.24H55.26Z\"/><path class=\"cls-25\" d=\"M66.58,32.62v9.66H62.06V25.38c0-4.24,2.81-7.24,6.79-7.24H87.94V23H68.85a2.12,2.12,0,0,0-2.26,2.41v2.41H87.94v4.83Z\"/><path class=\"cls-25\" d=\"M97,42.28c-4,0-6.79-3-6.79-7.24V18.14h4.53V35A2.12,2.12,0,0,0,97,37.45h19.09v4.83Z\"/><path class=\"cls-25\" d=\"M137.44,42.28H125.14c-4,0-6.79-3-6.79-7.24V25.38a7,7,0,0,1,6.79-7.24h12.29c3.88,0,6.79,3.1,6.79,7.24V35C144.23,39.35,141.48,42.28,137.44,42.28Zm2.26-16.9A2.18,2.18,0,0,0,137.44,23H125.14a2.29,2.29,0,0,0-2.26,2.41V35a2.12,2.12,0,0,0,2.26,2.41h12.29c1.55,0,2.26-.76,2.26-2.41Z\"/><path class=\"cls-25\" d=\"M174.48,39.76a3.34,3.34,0,0,1-6.47,0L164.19,27.1a.56.56,0,0,0-1,0l-3.82,12.66a3.34,3.34,0,0,1-6.47,0L146.5,18.14H151l4.56,15.17a.56.56,0,0,0,1,0l3.82-12.66a3.34,3.34,0,0,1,6.47,0l3.82,12.66a.56.56,0,0,0,1,0l4.56-15.17h4.53Z\"/></g></svg>" << endl;
    web << "</a>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"col-md-8\">" << endl;
    web << "<h1 class=\"page-title\"> Library of Crystallographic Prototypes</h1>" << endl;
    web << "<p class=\"page-subtitle\">" << endl;
    web << "<p class=\"page-cite\">M. J. Mehl, et al. (<a href=\"http://doi.org/doi=10.1016/j.commatsci.2017.01.017\" style=\"color:D0E5FF\">doi=10.1016/j.commatsci.2017.01.017</a>)</p>" << endl;
    web << "" << endl;
    web << "</p>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"container content-container\">" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<h2>" << endl;
    web << "Monoclinic (Cc) Low Tridymite (SiO<sub>2</sub>) Structure:   A2B_mC144_9_24a_12a" << endl;
    web << "</h2>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<a href=\"./PICS/A2B_mC144_9_24a_12a_composite_full.png\"><img src=\"./PICS/A2B_mC144_9_24a_12a_composite_full.png\" align=\"middle\" alt=\"Picture of Structure; Click for Big Picture\" class=\"img-res\"/></a>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<table class=\"tables\">" << endl;
    web << "<tbody>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Prototype</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>SiO<sub>2</sub></td" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong><span style=\"font-weight:400\">AFLOW</span> prototype label</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>A2B_mC144_9_24a_12a</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><em><strong>Strukturbericht designation</strong></em></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>None</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Pearson symbol</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>mC144</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Space group number</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>9</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Space group symbol</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>Cc</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong><span style=\"font-weight:400\">AFLOW</span> prototype command</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td><span class=\"monospaced\">aflow --proto=A2B_mC144_9_24a_12a <br>--params=</span>$a,b/a,c/a,\\beta,x_{1},y_{1},z_{1},x_{2},y_{2},z_{2},x_{3},y_{3},z_{3},x_{4},y_{4},z_{4},x_{5},y_{5},z_{5},x_{6},y_{6},z_{6},     \\\\x_{7},y_{7},z_{7},x_{8},y_{8},z_{8},x_{9},y_{9},z_{9},x_{10},y_{10},z_{10},x_{11},y_{11},z_{11},x_{12},y_{12},z_{12},x_{13},y_{13},z_{13},x_{14},     \\\\y_{14},z_{14},x_{15},y_{15},z_{15},x_{16},y_{16},z_{16},x_{17},y_{17},z_{17},x_{18},y_{18},z_{18},x_{19},y_{19},z_{19},x_{20},y_{20},z_{20},     \\\\x_{21},y_{21},z_{21},x_{22},y_{22},z_{22},x_{23},y_{23},z_{23},x_{24},y_{24},z_{24},x_{25},y_{25},z_{25},x_{26},y_{26},z_{26},x_{27},y_{27},     \\\\z_{27},x_{28},y_{28},z_{28},x_{29},y_{29},z_{29},x_{30},y_{30},z_{30},x_{31},y_{31},z_{31},x_{32},y_{32},z_{32},x_{33},y_{33},z_{33},x_{34},     \\\\y_{34},z_{34},x_{35},y_{35},z_{35},x_{36},y_{36},z_{36}$</td>" << endl;
    web << "</tr>" << endl;
    web << "</tbody>" << endl;
    web << "</table>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<a href=\"./PICS/A2B_mC144_9_24a_12a_composite_full.png\">View the structure from several different perspectives</a>" << endl;
    web << "<br>" << endl;
    web << "<a href=\"./PICS/A2B_mC144_9_24a_12a_prim.png\">View the primitive and conventional cell</a>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<h3><a href=\"./monoclinic_lattice.html#lattice3\">Base-centered Monoclinic</a> primitive vectors: </h3>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-6\" style=\"padding-top: 60px;\">" << endl;
    web << "\\[" << endl;
    web << "\\begin{array}{ccc}" << endl;
    web << "\\mathbf{a}_1&=&\\frac12\\,a\\,\\mathbf{\\hat{x}}-\\frac12\\,b\\,\\mathbf{\\hat{y}}\\\\" << endl;
    web << "\\mathbf{a}_2&=&\\frac12\\,a\\,\\mathbf{\\hat{x}}+\\frac12\\,b\\,\\mathbf{\\hat{y}}\\\\" << endl;
    web << "\\mathbf{a}_3&=&c\\cos\\beta\\,\\mathbf{\\hat{x}}+c\\sin\\beta\\,\\mathbf{\\hat{z}}" << endl;
    web << "\\end{array}" << endl;
    web << "\\]" << endl;
    web << "</div>" << endl;
    web << "<div class=\"col-md-6\">" << endl;
    web << "<a href=\"./PICS/A2B_mC144_9_24a_12a_prim.png\">" << endl;
    web << "<img src=\"./PICS/A2B_mC144_9_24a_12a_prim.png\" align=\"middle\" style=\"width:60%\" alt=\"Picture of Structure; Click for Big Picture\" />" << endl;
    web << "</a>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12 overflow-x-scroll\">" << endl;
    web << "<p>" << endl;
    web << "<strong>Basis vectors:</strong>" << endl;
    web << "</p>" << endl;
    web << "\\[" << endl;
    web << "\\begin{array}{ccccccc}" << endl;
    web << "&&\\mbox{Lattice Coordinates}&&\\mbox{Cartesian Coordinates}&\\mbox{Wyckoff Position}&\\mbox{Atom Type}\\\\" << endl;
    web << "\\mathbf{B_{1}}&=&\\left(x_{1}-y_{1}\\right)\\,\\mathbf{a_1}+\\left(x_{1}+y_{1}\\right)\\,\\mathbf{a_2}+z_{1}\\,\\mathbf{a_3}&=&\\left(x_{1}\\,a+z_{1}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{1}\\,b\\,\\mathbf{\\hat{y}}+z_{1}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O I}\\\\" << endl;
    web << "\\mathbf{B_{2}}&=&\\left(x_{1}+y_{1}\\right)\\,\\mathbf{a_1}+\\left(x_{1}-y_{1}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{1}\\right)\\,\\mathbf{a_3}&=&\\left(x_{1}\\,a+\\left(\\frac12+z_{1}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{1}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{1}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O I}\\\\" << endl;
    web << "\\mathbf{B_{3}}&=&\\left(x_{2}-y_{2}\\right)\\,\\mathbf{a_1}+\\left(x_{2}+y_{2}\\right)\\,\\mathbf{a_2}+z_{2}\\,\\mathbf{a_3}&=&\\left(x_{2}\\,a+z_{2}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{2}\\,b\\,\\mathbf{\\hat{y}}+z_{2}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O II}\\\\" << endl;
    web << "\\mathbf{B_{4}}&=&\\left(x_{2}+y_{2}\\right)\\,\\mathbf{a_1}+\\left(x_{2}-y_{2}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{2}\\right)\\,\\mathbf{a_3}&=&\\left(x_{2}\\,a+\\left(\\frac12+z_{2}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{2}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{2}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O II}\\\\" << endl;
    web << "\\mathbf{B_{5}}&=&\\left(x_{3}-y_{3}\\right)\\,\\mathbf{a_1}+\\left(x_{3}+y_{3}\\right)\\,\\mathbf{a_2}+z_{3}\\,\\mathbf{a_3}&=&\\left(x_{3}\\,a+z_{3}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{3}\\,b\\,\\mathbf{\\hat{y}}+z_{3}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O III}\\\\" << endl;
    web << "\\mathbf{B_{6}}&=&\\left(x_{3}+y_{3}\\right)\\,\\mathbf{a_1}+\\left(x_{3}-y_{3}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{3}\\right)\\,\\mathbf{a_3}&=&\\left(x_{3}\\,a+\\left(\\frac12+z_{3}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{3}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{3}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O III}\\\\" << endl;
    web << "\\mathbf{B_{7}}&=&\\left(x_{4}-y_{4}\\right)\\,\\mathbf{a_1}+\\left(x_{4}+y_{4}\\right)\\,\\mathbf{a_2}+z_{4}\\,\\mathbf{a_3}&=&\\left(x_{4}\\,a+z_{4}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{4}\\,b\\,\\mathbf{\\hat{y}}+z_{4}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O IV}\\\\" << endl;
    web << "\\mathbf{B_{8}}&=&\\left(x_{4}+y_{4}\\right)\\,\\mathbf{a_1}+\\left(x_{4}-y_{4}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{4}\\right)\\,\\mathbf{a_3}&=&\\left(x_{4}\\,a+\\left(\\frac12+z_{4}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{4}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{4}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O IV}\\\\" << endl;
    web << "\\mathbf{B_{9}}&=&\\left(x_{5}-y_{5}\\right)\\,\\mathbf{a_1}+\\left(x_{5}+y_{5}\\right)\\,\\mathbf{a_2}+z_{5}\\,\\mathbf{a_3}&=&\\left(x_{5}\\,a+z_{5}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{5}\\,b\\,\\mathbf{\\hat{y}}+z_{5}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O V}\\\\" << endl;
    web << "\\mathbf{B_{10}}&=&\\left(x_{5}+y_{5}\\right)\\,\\mathbf{a_1}+\\left(x_{5}-y_{5}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{5}\\right)\\,\\mathbf{a_3}&=&\\left(x_{5}\\,a+\\left(\\frac12+z_{5}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{5}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{5}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O V}\\\\" << endl;
    web << "\\mathbf{B_{11}}&=&\\left(x_{6}-y_{6}\\right)\\,\\mathbf{a_1}+\\left(x_{6}+y_{6}\\right)\\,\\mathbf{a_2}+z_{6}\\,\\mathbf{a_3}&=&\\left(x_{6}\\,a+z_{6}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{6}\\,b\\,\\mathbf{\\hat{y}}+z_{6}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O VI}\\\\" << endl;
    web << "\\mathbf{B_{12}}&=&\\left(x_{6}+y_{6}\\right)\\,\\mathbf{a_1}+\\left(x_{6}-y_{6}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{6}\\right)\\,\\mathbf{a_3}&=&\\left(x_{6}\\,a+\\left(\\frac12+z_{6}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{6}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{6}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O VI}\\\\" << endl;
    web << "\\mathbf{B_{13}}&=&\\left(x_{7}-y_{7}\\right)\\,\\mathbf{a_1}+\\left(x_{7}+y_{7}\\right)\\,\\mathbf{a_2}+z_{7}\\,\\mathbf{a_3}&=&\\left(x_{7}\\,a+z_{7}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{7}\\,b\\,\\mathbf{\\hat{y}}+z_{7}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O VII}\\\\" << endl;
    web << "\\mathbf{B_{14}}&=&\\left(x_{7}+y_{7}\\right)\\,\\mathbf{a_1}+\\left(x_{7}-y_{7}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{7}\\right)\\,\\mathbf{a_3}&=&\\left(x_{7}\\,a+\\left(\\frac12+z_{7}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{7}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{7}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O VII}\\\\" << endl;
    web << "\\mathbf{B_{15}}&=&\\left(x_{8}-y_{8}\\right)\\,\\mathbf{a_1}+\\left(x_{8}+y_{8}\\right)\\,\\mathbf{a_2}+z_{8}\\,\\mathbf{a_3}&=&\\left(x_{8}\\,a+z_{8}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{8}\\,b\\,\\mathbf{\\hat{y}}+z_{8}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O VIII}\\\\" << endl;
    web << "\\mathbf{B_{16}}&=&\\left(x_{8}+y_{8}\\right)\\,\\mathbf{a_1}+\\left(x_{8}-y_{8}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{8}\\right)\\,\\mathbf{a_3}&=&\\left(x_{8}\\,a+\\left(\\frac12+z_{8}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{8}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{8}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O VIII}\\\\" << endl;
    web << "\\mathbf{B_{17}}&=&\\left(x_{9}-y_{9}\\right)\\,\\mathbf{a_1}+\\left(x_{9}+y_{9}\\right)\\,\\mathbf{a_2}+z_{9}\\,\\mathbf{a_3}&=&\\left(x_{9}\\,a+z_{9}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{9}\\,b\\,\\mathbf{\\hat{y}}+z_{9}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O IX}\\\\" << endl;
    web << "\\mathbf{B_{18}}&=&\\left(x_{9}+y_{9}\\right)\\,\\mathbf{a_1}+\\left(x_{9}-y_{9}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{9}\\right)\\,\\mathbf{a_3}&=&\\left(x_{9}\\,a+\\left(\\frac12+z_{9}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{9}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{9}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O IX}\\\\" << endl;
    web << "\\mathbf{B_{19}}&=&\\left(x_{10}-y_{10}\\right)\\,\\mathbf{a_1}+\\left(x_{10}+y_{10}\\right)\\,\\mathbf{a_2}+z_{10}\\,\\mathbf{a_3}&=&\\left(x_{10}\\,a+z_{10}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{10}\\,b\\,\\mathbf{\\hat{y}}+z_{10}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O X}\\\\" << endl;
    web << "\\mathbf{B_{20}}&=&\\left(x_{10}+y_{10}\\right)\\,\\mathbf{a_1}+\\left(x_{10}-y_{10}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{10}\\right)\\,\\mathbf{a_3}&=&\\left(x_{10}\\,a+\\left(\\frac12+z_{10}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{10}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{10}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O X}\\\\" << endl;
    web << "\\mathbf{B_{21}}&=&\\left(x_{11}-y_{11}\\right)\\,\\mathbf{a_1}+\\left(x_{11}+y_{11}\\right)\\,\\mathbf{a_2}+z_{11}\\,\\mathbf{a_3}&=&\\left(x_{11}\\,a+z_{11}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{11}\\,b\\,\\mathbf{\\hat{y}}+z_{11}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XI}\\\\" << endl;
    web << "\\mathbf{B_{22}}&=&\\left(x_{11}+y_{11}\\right)\\,\\mathbf{a_1}+\\left(x_{11}-y_{11}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{11}\\right)\\,\\mathbf{a_3}&=&\\left(x_{11}\\,a+\\left(\\frac12+z_{11}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{11}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{11}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XI}\\\\" << endl;
    web << "\\mathbf{B_{23}}&=&\\left(x_{12}-y_{12}\\right)\\,\\mathbf{a_1}+\\left(x_{12}+y_{12}\\right)\\,\\mathbf{a_2}+z_{12}\\,\\mathbf{a_3}&=&\\left(x_{12}\\,a+z_{12}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{12}\\,b\\,\\mathbf{\\hat{y}}+z_{12}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XII}\\\\" << endl;
    web << "\\mathbf{B_{24}}&=&\\left(x_{12}+y_{12}\\right)\\,\\mathbf{a_1}+\\left(x_{12}-y_{12}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{12}\\right)\\,\\mathbf{a_3}&=&\\left(x_{12}\\,a+\\left(\\frac12+z_{12}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{12}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{12}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XII}\\\\" << endl;
    web << "\\mathbf{B_{25}}&=&\\left(x_{13}-y_{13}\\right)\\,\\mathbf{a_1}+\\left(x_{13}+y_{13}\\right)\\,\\mathbf{a_2}+z_{13}\\,\\mathbf{a_3}&=&\\left(x_{13}\\,a+z_{13}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{13}\\,b\\,\\mathbf{\\hat{y}}+z_{13}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XIII}\\\\" << endl;
    web << "\\mathbf{B_{26}}&=&\\left(x_{13}+y_{13}\\right)\\,\\mathbf{a_1}+\\left(x_{13}-y_{13}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{13}\\right)\\,\\mathbf{a_3}&=&\\left(x_{13}\\,a+\\left(\\frac12+z_{13}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{13}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{13}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XIII}\\\\" << endl;
    web << "\\mathbf{B_{27}}&=&\\left(x_{14}-y_{14}\\right)\\,\\mathbf{a_1}+\\left(x_{14}+y_{14}\\right)\\,\\mathbf{a_2}+z_{14}\\,\\mathbf{a_3}&=&\\left(x_{14}\\,a+z_{14}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{14}\\,b\\,\\mathbf{\\hat{y}}+z_{14}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XIV}\\\\" << endl;
    web << "\\mathbf{B_{28}}&=&\\left(x_{14}+y_{14}\\right)\\,\\mathbf{a_1}+\\left(x_{14}-y_{14}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{14}\\right)\\,\\mathbf{a_3}&=&\\left(x_{14}\\,a+\\left(\\frac12+z_{14}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{14}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{14}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XIV}\\\\" << endl;
    web << "\\mathbf{B_{29}}&=&\\left(x_{15}-y_{15}\\right)\\,\\mathbf{a_1}+\\left(x_{15}+y_{15}\\right)\\,\\mathbf{a_2}+z_{15}\\,\\mathbf{a_3}&=&\\left(x_{15}\\,a+z_{15}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{15}\\,b\\,\\mathbf{\\hat{y}}+z_{15}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XV}\\\\" << endl;
    web << "\\mathbf{B_{30}}&=&\\left(x_{15}+y_{15}\\right)\\,\\mathbf{a_1}+\\left(x_{15}-y_{15}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{15}\\right)\\,\\mathbf{a_3}&=&\\left(x_{15}\\,a+\\left(\\frac12+z_{15}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{15}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{15}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XV}\\\\" << endl;
    web << "\\mathbf{B_{31}}&=&\\left(x_{16}-y_{16}\\right)\\,\\mathbf{a_1}+\\left(x_{16}+y_{16}\\right)\\,\\mathbf{a_2}+z_{16}\\,\\mathbf{a_3}&=&\\left(x_{16}\\,a+z_{16}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{16}\\,b\\,\\mathbf{\\hat{y}}+z_{16}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XVI}\\\\" << endl;
    web << "\\mathbf{B_{32}}&=&\\left(x_{16}+y_{16}\\right)\\,\\mathbf{a_1}+\\left(x_{16}-y_{16}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{16}\\right)\\,\\mathbf{a_3}&=&\\left(x_{16}\\,a+\\left(\\frac12+z_{16}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{16}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{16}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XVI}\\\\" << endl;
    web << "\\mathbf{B_{33}}&=&\\left(x_{17}-y_{17}\\right)\\,\\mathbf{a_1}+\\left(x_{17}+y_{17}\\right)\\,\\mathbf{a_2}+z_{17}\\,\\mathbf{a_3}&=&\\left(x_{17}\\,a+z_{17}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{17}\\,b\\,\\mathbf{\\hat{y}}+z_{17}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XVII}\\\\" << endl;
    web << "\\mathbf{B_{34}}&=&\\left(x_{17}+y_{17}\\right)\\,\\mathbf{a_1}+\\left(x_{17}-y_{17}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{17}\\right)\\,\\mathbf{a_3}&=&\\left(x_{17}\\,a+\\left(\\frac12+z_{17}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{17}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{17}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XVII}\\\\" << endl;
    web << "\\mathbf{B_{35}}&=&\\left(x_{18}-y_{18}\\right)\\,\\mathbf{a_1}+\\left(x_{18}+y_{18}\\right)\\,\\mathbf{a_2}+z_{18}\\,\\mathbf{a_3}&=&\\left(x_{18}\\,a+z_{18}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{18}\\,b\\,\\mathbf{\\hat{y}}+z_{18}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XVIII}\\\\" << endl;
    web << "\\mathbf{B_{36}}&=&\\left(x_{18}+y_{18}\\right)\\,\\mathbf{a_1}+\\left(x_{18}-y_{18}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{18}\\right)\\,\\mathbf{a_3}&=&\\left(x_{18}\\,a+\\left(\\frac12+z_{18}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{18}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{18}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XVIII}\\\\" << endl;
    web << "\\mathbf{B_{37}}&=&\\left(x_{19}-y_{19}\\right)\\,\\mathbf{a_1}+\\left(x_{19}+y_{19}\\right)\\,\\mathbf{a_2}+z_{19}\\,\\mathbf{a_3}&=&\\left(x_{19}\\,a+z_{19}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{19}\\,b\\,\\mathbf{\\hat{y}}+z_{19}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XIX}\\\\" << endl;
    web << "\\mathbf{B_{38}}&=&\\left(x_{19}+y_{19}\\right)\\,\\mathbf{a_1}+\\left(x_{19}-y_{19}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{19}\\right)\\,\\mathbf{a_3}&=&\\left(x_{19}\\,a+\\left(\\frac12+z_{19}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{19}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{19}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XIX}\\\\" << endl;
    web << "\\mathbf{B_{39}}&=&\\left(x_{20}-y_{20}\\right)\\,\\mathbf{a_1}+\\left(x_{20}+y_{20}\\right)\\,\\mathbf{a_2}+z_{20}\\,\\mathbf{a_3}&=&\\left(x_{20}\\,a+z_{20}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{20}\\,b\\,\\mathbf{\\hat{y}}+z_{20}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XX}\\\\" << endl;
    web << "\\mathbf{B_{40}}&=&\\left(x_{20}+y_{20}\\right)\\,\\mathbf{a_1}+\\left(x_{20}-y_{20}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{20}\\right)\\,\\mathbf{a_3}&=&\\left(x_{20}\\,a+\\left(\\frac12+z_{20}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{20}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{20}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XX}\\\\" << endl;
    web << "\\mathbf{B_{41}}&=&\\left(x_{21}-y_{21}\\right)\\,\\mathbf{a_1}+\\left(x_{21}+y_{21}\\right)\\,\\mathbf{a_2}+z_{21}\\,\\mathbf{a_3}&=&\\left(x_{21}\\,a+z_{21}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{21}\\,b\\,\\mathbf{\\hat{y}}+z_{21}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXI}\\\\" << endl;
    web << "\\mathbf{B_{42}}&=&\\left(x_{21}+y_{21}\\right)\\,\\mathbf{a_1}+\\left(x_{21}-y_{21}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{21}\\right)\\,\\mathbf{a_3}&=&\\left(x_{21}\\,a+\\left(\\frac12+z_{21}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{21}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{21}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXI}\\\\" << endl;
    web << "\\mathbf{B_{43}}&=&\\left(x_{22}-y_{22}\\right)\\,\\mathbf{a_1}+\\left(x_{22}+y_{22}\\right)\\,\\mathbf{a_2}+z_{22}\\,\\mathbf{a_3}&=&\\left(x_{22}\\,a+z_{22}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{22}\\,b\\,\\mathbf{\\hat{y}}+z_{22}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXII}\\\\" << endl;
    web << "\\mathbf{B_{44}}&=&\\left(x_{22}+y_{22}\\right)\\,\\mathbf{a_1}+\\left(x_{22}-y_{22}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{22}\\right)\\,\\mathbf{a_3}&=&\\left(x_{22}\\,a+\\left(\\frac12+z_{22}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{22}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{22}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXII}\\\\" << endl;
    web << "\\mathbf{B_{45}}&=&\\left(x_{23}-y_{23}\\right)\\,\\mathbf{a_1}+\\left(x_{23}+y_{23}\\right)\\,\\mathbf{a_2}+z_{23}\\,\\mathbf{a_3}&=&\\left(x_{23}\\,a+z_{23}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{23}\\,b\\,\\mathbf{\\hat{y}}+z_{23}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXIII}\\\\" << endl;
    web << "\\mathbf{B_{46}}&=&\\left(x_{23}+y_{23}\\right)\\,\\mathbf{a_1}+\\left(x_{23}-y_{23}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{23}\\right)\\,\\mathbf{a_3}&=&\\left(x_{23}\\,a+\\left(\\frac12+z_{23}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{23}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{23}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXIII}\\\\" << endl;
    web << "\\mathbf{B_{47}}&=&\\left(x_{24}-y_{24}\\right)\\,\\mathbf{a_1}+\\left(x_{24}+y_{24}\\right)\\,\\mathbf{a_2}+z_{24}\\,\\mathbf{a_3}&=&\\left(x_{24}\\,a+z_{24}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{24}\\,b\\,\\mathbf{\\hat{y}}+z_{24}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXIV}\\\\" << endl;
    web << "\\mathbf{B_{48}}&=&\\left(x_{24}+y_{24}\\right)\\,\\mathbf{a_1}+\\left(x_{24}-y_{24}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{24}\\right)\\,\\mathbf{a_3}&=&\\left(x_{24}\\,a+\\left(\\frac12+z_{24}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{24}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{24}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{O XXIV}\\\\" << endl;
    web << "\\mathbf{B_{49}}&=&\\left(x_{25}-y_{25}\\right)\\,\\mathbf{a_1}+\\left(x_{25}+y_{25}\\right)\\,\\mathbf{a_2}+z_{25}\\,\\mathbf{a_3}&=&\\left(x_{25}\\,a+z_{25}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{25}\\,b\\,\\mathbf{\\hat{y}}+z_{25}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si I}\\\\" << endl;
    web << "\\mathbf{B_{50}}&=&\\left(x_{25}+y_{25}\\right)\\,\\mathbf{a_1}+\\left(x_{25}-y_{25}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{25}\\right)\\,\\mathbf{a_3}&=&\\left(x_{25}\\,a+\\left(\\frac12+z_{25}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{25}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{25}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si I}\\\\" << endl;
    web << "\\mathbf{B_{51}}&=&\\left(x_{26}-y_{26}\\right)\\,\\mathbf{a_1}+\\left(x_{26}+y_{26}\\right)\\,\\mathbf{a_2}+z_{26}\\,\\mathbf{a_3}&=&\\left(x_{26}\\,a+z_{26}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{26}\\,b\\,\\mathbf{\\hat{y}}+z_{26}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si II}\\\\" << endl;
    web << "\\mathbf{B_{52}}&=&\\left(x_{26}+y_{26}\\right)\\,\\mathbf{a_1}+\\left(x_{26}-y_{26}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{26}\\right)\\,\\mathbf{a_3}&=&\\left(x_{26}\\,a+\\left(\\frac12+z_{26}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{26}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{26}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si II}\\\\" << endl;
    web << "\\mathbf{B_{53}}&=&\\left(x_{27}-y_{27}\\right)\\,\\mathbf{a_1}+\\left(x_{27}+y_{27}\\right)\\,\\mathbf{a_2}+z_{27}\\,\\mathbf{a_3}&=&\\left(x_{27}\\,a+z_{27}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{27}\\,b\\,\\mathbf{\\hat{y}}+z_{27}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si III}\\\\" << endl;
    web << "\\mathbf{B_{54}}&=&\\left(x_{27}+y_{27}\\right)\\,\\mathbf{a_1}+\\left(x_{27}-y_{27}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{27}\\right)\\,\\mathbf{a_3}&=&\\left(x_{27}\\,a+\\left(\\frac12+z_{27}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{27}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{27}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si III}\\\\" << endl;
    web << "\\mathbf{B_{55}}&=&\\left(x_{28}-y_{28}\\right)\\,\\mathbf{a_1}+\\left(x_{28}+y_{28}\\right)\\,\\mathbf{a_2}+z_{28}\\,\\mathbf{a_3}&=&\\left(x_{28}\\,a+z_{28}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{28}\\,b\\,\\mathbf{\\hat{y}}+z_{28}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si IV}\\\\" << endl;
    web << "\\mathbf{B_{56}}&=&\\left(x_{28}+y_{28}\\right)\\,\\mathbf{a_1}+\\left(x_{28}-y_{28}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{28}\\right)\\,\\mathbf{a_3}&=&\\left(x_{28}\\,a+\\left(\\frac12+z_{28}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{28}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{28}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si IV}\\\\" << endl;
    web << "\\mathbf{B_{57}}&=&\\left(x_{29}-y_{29}\\right)\\,\\mathbf{a_1}+\\left(x_{29}+y_{29}\\right)\\,\\mathbf{a_2}+z_{29}\\,\\mathbf{a_3}&=&\\left(x_{29}\\,a+z_{29}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{29}\\,b\\,\\mathbf{\\hat{y}}+z_{29}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si V}\\\\" << endl;
    web << "\\mathbf{B_{58}}&=&\\left(x_{29}+y_{29}\\right)\\,\\mathbf{a_1}+\\left(x_{29}-y_{29}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{29}\\right)\\,\\mathbf{a_3}&=&\\left(x_{29}\\,a+\\left(\\frac12+z_{29}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{29}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{29}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si V}\\\\" << endl;
    web << "\\mathbf{B_{59}}&=&\\left(x_{30}-y_{30}\\right)\\,\\mathbf{a_1}+\\left(x_{30}+y_{30}\\right)\\,\\mathbf{a_2}+z_{30}\\,\\mathbf{a_3}&=&\\left(x_{30}\\,a+z_{30}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{30}\\,b\\,\\mathbf{\\hat{y}}+z_{30}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si VI}\\\\" << endl;
    web << "\\mathbf{B_{60}}&=&\\left(x_{30}+y_{30}\\right)\\,\\mathbf{a_1}+\\left(x_{30}-y_{30}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{30}\\right)\\,\\mathbf{a_3}&=&\\left(x_{30}\\,a+\\left(\\frac12+z_{30}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{30}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{30}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si VI}\\\\" << endl;
    web << "\\mathbf{B_{61}}&=&\\left(x_{31}-y_{31}\\right)\\,\\mathbf{a_1}+\\left(x_{31}+y_{31}\\right)\\,\\mathbf{a_2}+z_{31}\\,\\mathbf{a_3}&=&\\left(x_{31}\\,a+z_{31}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{31}\\,b\\,\\mathbf{\\hat{y}}+z_{31}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si VII}\\\\" << endl;
    web << "\\mathbf{B_{62}}&=&\\left(x_{31}+y_{31}\\right)\\,\\mathbf{a_1}+\\left(x_{31}-y_{31}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{31}\\right)\\,\\mathbf{a_3}&=&\\left(x_{31}\\,a+\\left(\\frac12+z_{31}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{31}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{31}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si VII}\\\\" << endl;
    web << "\\mathbf{B_{63}}&=&\\left(x_{32}-y_{32}\\right)\\,\\mathbf{a_1}+\\left(x_{32}+y_{32}\\right)\\,\\mathbf{a_2}+z_{32}\\,\\mathbf{a_3}&=&\\left(x_{32}\\,a+z_{32}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{32}\\,b\\,\\mathbf{\\hat{y}}+z_{32}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si VIII}\\\\" << endl;
    web << "\\mathbf{B_{64}}&=&\\left(x_{32}+y_{32}\\right)\\,\\mathbf{a_1}+\\left(x_{32}-y_{32}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{32}\\right)\\,\\mathbf{a_3}&=&\\left(x_{32}\\,a+\\left(\\frac12+z_{32}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{32}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{32}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si VIII}\\\\" << endl;
    web << "\\mathbf{B_{65}}&=&\\left(x_{33}-y_{33}\\right)\\,\\mathbf{a_1}+\\left(x_{33}+y_{33}\\right)\\,\\mathbf{a_2}+z_{33}\\,\\mathbf{a_3}&=&\\left(x_{33}\\,a+z_{33}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{33}\\,b\\,\\mathbf{\\hat{y}}+z_{33}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si IX}\\\\" << endl;
    web << "\\mathbf{B_{66}}&=&\\left(x_{33}+y_{33}\\right)\\,\\mathbf{a_1}+\\left(x_{33}-y_{33}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{33}\\right)\\,\\mathbf{a_3}&=&\\left(x_{33}\\,a+\\left(\\frac12+z_{33}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{33}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{33}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si IX}\\\\" << endl;
    web << "\\mathbf{B_{67}}&=&\\left(x_{34}-y_{34}\\right)\\,\\mathbf{a_1}+\\left(x_{34}+y_{34}\\right)\\,\\mathbf{a_2}+z_{34}\\,\\mathbf{a_3}&=&\\left(x_{34}\\,a+z_{34}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{34}\\,b\\,\\mathbf{\\hat{y}}+z_{34}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si X}\\\\" << endl;
    web << "\\mathbf{B_{68}}&=&\\left(x_{34}+y_{34}\\right)\\,\\mathbf{a_1}+\\left(x_{34}-y_{34}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{34}\\right)\\,\\mathbf{a_3}&=&\\left(x_{34}\\,a+\\left(\\frac12+z_{34}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{34}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{34}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si X}\\\\" << endl;
    web << "\\mathbf{B_{69}}&=&\\left(x_{35}-y_{35}\\right)\\,\\mathbf{a_1}+\\left(x_{35}+y_{35}\\right)\\,\\mathbf{a_2}+z_{35}\\,\\mathbf{a_3}&=&\\left(x_{35}\\,a+z_{35}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{35}\\,b\\,\\mathbf{\\hat{y}}+z_{35}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si XI}\\\\" << endl;
    web << "\\mathbf{B_{70}}&=&\\left(x_{35}+y_{35}\\right)\\,\\mathbf{a_1}+\\left(x_{35}-y_{35}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{35}\\right)\\,\\mathbf{a_3}&=&\\left(x_{35}\\,a+\\left(\\frac12+z_{35}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{35}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{35}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si XI}\\\\" << endl;
    web << "\\mathbf{B_{71}}&=&\\left(x_{36}-y_{36}\\right)\\,\\mathbf{a_1}+\\left(x_{36}+y_{36}\\right)\\,\\mathbf{a_2}+z_{36}\\,\\mathbf{a_3}&=&\\left(x_{36}\\,a+z_{36}\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}+y_{36}\\,b\\,\\mathbf{\\hat{y}}+z_{36}\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si XII}\\\\" << endl;
    web << "\\mathbf{B_{72}}&=&\\left(x_{36}+y_{36}\\right)\\,\\mathbf{a_1}+\\left(x_{36}-y_{36}\\right)\\,\\mathbf{a_2}+\\left(\\frac12+z_{36}\\right)\\,\\mathbf{a_3}&=&\\left(x_{36}\\,a+\\left(\\frac12+z_{36}\\right)\\,c\\,\\cos\\beta\\right)\\,\\mathbf{\\hat{x}}-y_{36}\\,b\\,\\mathbf{\\hat{y}}+\\left(\\frac12+z_{36}\\right)\\,c\\,\\sin\\beta\\,\\mathbf{\\hat{z}}&\\left(4a\\right)&\\mbox{Si XII}\\\\" << endl;
    web << "\\end{array}" << endl;
    web << "\\]" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<h3>References</h3>" << endl;
    web << "<ul>" << endl;
    web << "<li> W. A. Dollase and W. H. Baur, <em>The superstructure of meteoritic low tridymite solved by computer simulation</em>, Am. Mineral. <strong>61</strong>, 971&ndash;978 (1976).</li>" << endl;
    web << "</ul>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<h3>Geometry files</h3>" << endl;
    web << "<ul>" << endl;
    web << "<li><a href=\"./CIF/A2B_mC144_9_24a_12a.cif\">CIF</a></li>" << endl;
    web << "<li><a href=\"./POSCAR/POSCAR.A2B_mC144_9_24a_12a\">POSCAR</a></li>" << endl;
    web << "</ul>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<!-- /#page -->" << endl;
    web << "<footer class=\"footer\">" << endl;
    web << "Citation Information:" << endl;
    web << "<ul>" << endl;
    web << "M. J. Mehl, D. Hicks, C. Toher, O. Levy, R. M. Hanson, G. L. W. Hart, and S. Curtarolo," << endl;
    web << "The AFLOW Library of Crystallographic Prototypes, in press Comp. Mat. Sci. (2017). " << endl;
    web << "(<a href=\"http://doi.org/doi=10.1016/j.commatsci.2017.01.017\" style=\"color:D0E5FF\">doi=10.1016/j.commatsci.2017.01.017</a>)" << endl;
    web << "</ul>" << endl;
    web << "<p class=\"center\"><a href=\"./index.html\">Return to the AFLOW Library of Crystallographic Prototypes Home Page</a></p>" << endl;
    web << "</footer>" << endl;
    web << "</body>" << endl;
    #endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A2B_mC144_9_24a_12a: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************