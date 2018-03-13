// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2016
// FILE "ANRL/aflow_arnl_A_tP32_138_j.cpp"

#ifndef _AFLOW_ANRL_A_tP32_138_j_CPP
#define _AFLOW_ANRL_A_tP32_138_j_CPP
#include "../aflow.h"


namespace aflowlib {
  uint PrototypeLibrariesANRL_A_tP32_138_j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG) {
    // system A_tP32_138_j

    if(XHOST.vflag_control.flag("WWW")) {
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
    string label,Pearson_symbol,params,Strukturbericht,prototype;

    aflowlib::vproto_line2tokens(proto_line,label,nspecies,natoms,spacegroup,nunderscores,nparameters,Pearson_symbol,params,Strukturbericht,prototype);

    if(!aflowlib::PrototypeLibrariesANRL_Consistency(oss,vparameters.size(),nparameters,prototype,label,
                 Strukturbericht,Pearson_symbol,spacegroup, params)) exit(0);    

    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: FOUND" << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: label=" << label << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: nspecies=" << nspecies << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: natoms=" << natoms << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: spacegroup=" << spacegroup << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: nunderscores=" << nunderscores << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: nparameters=" <<  nparameters << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: Pearson_symbol=" << Pearson_symbol << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: params=" << params << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: Strukturbericht=" << Strukturbericht << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: prototype=" << prototype << endl;
    if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: vparameters.size()=" << vparameters.size() << endl;

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    _atom atom;

    uint i=0;
    double a=vparameters.at(i++);                  if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: a=" << a << endl;
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: c=" << c << " (c/a=" << covera << ")" << endl;
    
    double x1=vparameters.at(i++);                 if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: x1=" << x1 << endl;
    double y1=vparameters.at(i++);                 if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: y1=" << y1 << endl;
    double z1=vparameters.at(i++);                 if(LDEBUG) cerr << "aflowlib::PrototypeLibrariesANRL_A_tP32_138_j: z1=" << z1 << endl;
        
    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG#="+aurostd::utype2string(spacegroup)+DOI_ANRL;
    str.scale=1.0;

    a1=a*xn;
    a2=a*yn;
    a3=c*zn;
    
    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    str.FixLattices(); // Reciprocal/f2c/c2f
    
    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=x1;atom.fpos(2)=y1;atom.fpos(3)=z1;                     // atom B1
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1
    
    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=((1.0/2.0)-x1);atom.fpos(2)=((1.0/2.0)-y1);atom.fpos(3)=z1;                     // atom B2
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2
    
    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=((1.0/2.0)-y1);atom.fpos(2)=x1;atom.fpos(3)=((1.0/2.0)+z1);                     // atom B3
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3
    
    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=y1;atom.fpos(2)=((1.0/2.0)-x1);atom.fpos(3)=((1.0/2.0)+z1);                     // atom B4
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4
    
    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=-x1;atom.fpos(2)=((1.0/2.0)+y1);atom.fpos(3)=((1.0/2.0)-z1);                     // atom B5
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5
    
    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=((1.0/2.0)+x1);atom.fpos(2)=-y1;atom.fpos(3)=((1.0/2.0)-z1);                     // atom B6
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6
    
    atom.name="A"; atom.type=0;                                       // atom B7
    atom.fpos(1)=((1.0/2.0)+y1);atom.fpos(2)=((1.0/2.0)+x1);atom.fpos(3)=-z1;                     // atom B7
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7
    
    atom.name="A"; atom.type=0;                                       // atom B8
    atom.fpos(1)=-y1;atom.fpos(2)=-x1;atom.fpos(3)=-z1;                     // atom B8
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8
    
    atom.name="A"; atom.type=0;                                       // atom B9
    atom.fpos(1)=-x1;atom.fpos(2)=-y1;atom.fpos(3)=-z1;                     // atom B9
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9
    
    atom.name="A"; atom.type=0;                                       // atom B10
    atom.fpos(1)=((1.0/2.0)+x1);atom.fpos(2)=((1.0/2.0)+y1);atom.fpos(3)=-z1;                     // atom B10
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10
    
    atom.name="A"; atom.type=0;                                       // atom B11
    atom.fpos(1)=((1.0/2.0)+y1);atom.fpos(2)=-x1;atom.fpos(3)=((1.0/2.0)-z1);                     // atom B11
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11
    
    atom.name="A"; atom.type=0;                                       // atom B12
    atom.fpos(1)=-y1;atom.fpos(2)=((1.0/2.0)+x1);atom.fpos(3)=((1.0/2.0)-z1);                     // atom B12
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12
    
    atom.name="A"; atom.type=0;                                       // atom B13
    atom.fpos(1)=x1;atom.fpos(2)=((1.0/2.0)-y1);atom.fpos(3)=((1.0/2.0)+z1);                     // atom B13
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13
    
    atom.name="A"; atom.type=0;                                       // atom B14
    atom.fpos(1)=((1.0/2.0)-x1);atom.fpos(2)=y1;atom.fpos(3)=((1.0/2.0)+z1);                     // atom B14
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14
    
    atom.name="A"; atom.type=0;                                       // atom B15
    atom.fpos(1)=((1.0/2.0)-y1);atom.fpos(2)=((1.0/2.0)-x1);atom.fpos(3)=z1;                     // atom B15
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15
    
    atom.name="A"; atom.type=0;                                       // atom B16
    atom.fpos(1)=y1;atom.fpos(2)=x1;atom.fpos(3)=z1;                     // atom B16
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16
    

    return str.atoms.size();  
  }
} // namespace aflowlib::

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
