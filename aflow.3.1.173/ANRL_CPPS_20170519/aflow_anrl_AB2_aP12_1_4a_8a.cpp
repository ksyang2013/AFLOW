// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2016
// FILE "ANRL/aflow_anrl_AB2_aP12_1_4a_8a.cpp"

#ifndef _AFLOW_ANRL_AB2_aP12_1_4a_8a_CPP
#define _AFLOW_ANRL_AB2_aP12_1_4a_8a_CPP
#include "../aflow.h"

namespace anrl {
  uint WebANRL_AB2_aP12_1_4a_8a(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_AB2_aP12_1_4a_8a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG) {
    // system AB2_aP12_1_4a_8a

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_AB2_aP12_1_4a_8a(web,LDEBUG); // PLUG WEB STUFF
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
  
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: vparameters.size()=" << vparameters.size() << endl;}
  
    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);
  
    _atom atom;
  
    uint i=0;
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: c=" << c << " (c/a=" << covera << ")" << endl;}
    double alpha=vparameters.at(i++);              if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: alpha=" << alpha << endl;}
    double beta=vparameters.at(i++);               if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: beta=" << beta << endl;}
    double gamma=vparameters.at(i++);              if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: gamma=" << gamma << endl;}
  
    double x1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x1=" << x1 << endl;}
    double y1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y1=" << y1 << endl;}
    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z1=" << z1 << endl;}
    double x2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x2=" << x2 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y2=" << y2 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z2=" << z2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x3=" << x3 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y3=" << y3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z3=" << z3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z4=" << z4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y5=" << y5 << endl;}
    double z5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z5=" << z5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y6=" << y6 << endl;}
    double z6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z6=" << z6 << endl;}
    double x7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x7=" << x7 << endl;}
    double y7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y7=" << y7 << endl;}
    double z7=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z7=" << z7 << endl;}
    double x8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x8=" << x8 << endl;}
    double y8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y8=" << y8 << endl;}
    double z8=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z8=" << z8 << endl;}
    double x9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x9=" << x9 << endl;}
    double y9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y9=" << y9 << endl;}
    double z9=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z9=" << z9 << endl;}
    double x10=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x10=" << x10 << endl;}
    double y10=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y10=" << y10 << endl;}
    double z10=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z10=" << z10 << endl;}
    double x11=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x11=" << x11 << endl;}
    double y11=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y11=" << y11 << endl;}
    double z11=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z11=" << z11 << endl;}
    double x12=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: x12=" << x12 << endl;}
    double y12=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: y12=" << y12 << endl;}
    double z12=vparameters.at(i++);                if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: z12=" << z12 << endl;}
    
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: cos(alpha)=" << cos(deg2rad*alpha) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: sin(alpha)=" << sin(deg2rad*alpha) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: cos(beta)=" << cos(deg2rad*beta) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: sin(beta)=" << sin(deg2rad*beta) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: cos(gamma)=" << cos(deg2rad*gamma) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2_aP12_1_4a_8a: sin(gamma)=" << sin(deg2rad*gamma) << endl;}
    
    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG#="+aurostd::utype2string(spacegroup)+DOI_ANRL;
    str.scale=1.0;

    double cx=c*cos(deg2rad*beta);
    double cy=c*(cos(deg2rad*alpha)-cos(deg2rad*beta)*cos(deg2rad*gamma))/sin(deg2rad*gamma);
    double cz=sqrt(pow(c,2.0)-pow(cx,2.0)-pow(cy,2.0));

    a1=a*xn;
    a2=b*cos(deg2rad*gamma)*xn+b*sin(deg2rad*gamma)*yn;
    a3=cx*xn+cy*yn+cz*zn;
  
    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    str.FixLattices(); // Reciprocal/f2c/c2f
  
    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=x1;atom.fpos(2)=y1;atom.fpos(3)=z1;                  // atom B1
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1
  
    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=x2;atom.fpos(2)=y2;atom.fpos(3)=z2;                  // atom B2
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2
  
    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=x3;atom.fpos(2)=y3;atom.fpos(3)=z3;                  // atom B3
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3
  
    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=x4;atom.fpos(2)=y4;atom.fpos(3)=z4;                  // atom B4
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4
  
    atom.name="B"; atom.type=1;                                       // atom B5
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=z5;                  // atom B5
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5
  
    atom.name="B"; atom.type=1;                                       // atom B6
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=z6;                  // atom B6
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6
  
    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=x7;atom.fpos(2)=y7;atom.fpos(3)=z7;                  // atom B7
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7
  
    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=x8;atom.fpos(2)=y8;atom.fpos(3)=z8;                  // atom B8
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8
  
    atom.name="B"; atom.type=1;                                       // atom B9
    atom.fpos(1)=x9;atom.fpos(2)=y9;atom.fpos(3)=z9;                  // atom B9
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9
  
    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=x10;atom.fpos(2)=y10;atom.fpos(3)=z10;               // atom B10
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10
  
    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=x11;atom.fpos(2)=y11;atom.fpos(3)=z11;               // atom B11
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11
  
    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=x12;atom.fpos(2)=y12;atom.fpos(3)=z12;               // atom B12
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12
  

    return str.atoms.size();  
  }
} // namespace anrl::

namespace anrl {
  uint WebANRL_AB2_aP12_1_4a_8a(stringstream& web,bool LDEBUG) {
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
    web << "FeS<sub>2</sub> (P1) Structure: AB2_aP12_1_4a_8a" << endl;
    web << "</h2>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<a href=\"./PICS/AB2_aP12_1_4a_8a_composite_full.png\"><img src=\"./PICS/AB2_aP12_1_4a_8a_composite_full.png\" align=\"middle\" alt=\"Picture of Structure; Click for Big Picture\" class=\"img-res\"/></a>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<table class=\"tables\">" << endl;
    web << "<tbody>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Prototype</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>FeS<sub>2</sub></td" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong><span style=\"font-weight:400\">AFLOW</span> prototype label</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>AB2_aP12_1_4a_8a</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><em><strong>Strukturbericht designation</strong></em></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>None</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Pearson symbol</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>aP12</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Space group number</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>1</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong>Space group symbol</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td>P1</td>" << endl;
    web << "</tr>" << endl;
    web << "<tr>" << endl;
    web << "<td><strong><span style=\"font-weight:400\">AFLOW</span> prototype command</strong></td>" << endl;
    web << "<td>:</td>" << endl;
    web << "<td><span class=\"monospaced\">aflow --proto=AB2_aP12_1_4a_8a <br>--params=</span>$a,b/a,c/a,\\alpha,\\beta,\\gamma,x_{1},y_{1},z_{1},x_{2},y_{2},z_{2},x_{3},y_{3},z_{3},x_{4},y_{4},z_{4},x_{5},y_{5},z_{5},x_{6},     \\\\y_{6},z_{6},x_{7},y_{7},z_{7},x_{8},y_{8},z_{8},x_{9},y_{9},z_{9},x_{10},y_{10},z_{10},x_{11},y_{11},z_{11},x_{12},y_{12},z_{12}$</td>" << endl;
    web << "</tr>" << endl;
    web << "</tbody>" << endl;
    web << "</table>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<a href=\"./PICS/AB2_aP12_1_4a_8a_composite_full.png\">View the structure from several different perspectives</a>" << endl;
    web << "<br>" << endl;
    web << "<a href=\"./PICS/AB2_aP12_1_4a_8a_prim.png\">View the primitive and conventional cell</a>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<ul>" << endl;
    web << "<li>This structure is just a slightly distorted version of <a href=\"AB2_cP12_205_a_c.html\">pyrite (C2)</a>, with no rotational symmetry whatsoever. </li>" << endl;
    web << "</ul>" << endl;
    web << "</li>" << endl;
    web << "</ul>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<h3><a href=\"./triclinic_lattice.html#lattice1\">Triclinic</a> primitive vectors: </h3>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-6\" style=\"padding-top: 60px;\">" << endl;
    web << "\\[" << endl;
    web << "\\begin{array}{ccc}" << endl;
    web << "\\mathbf{a}_1&=&a\\mathbf{\\hat{x}}\\\\" << endl;
    web << "\\mathbf{a}_2&=&b\\cos\\gamma\\,\\mathbf{\\hat{x}}+b\\sin\\gamma\\,\\mathbf{\\hat{y}}\\\\" << endl;
    web << "\\mathbf{a}_3&=&c_x\\mathbf{\\hat{x}}+c_y\\,\\mathbf{\\hat{y}}+c_z\\,\\mathbf{\\hat{z}}" << endl;
    web << "\\end{array}" << endl;
    web << "\\]" << endl;
    web << "</div>" << endl;
    web << "<div class=\"col-md-6\">" << endl;
    web << "<a href=\"./PICS/AB2_aP12_1_4a_8a_prim.png\">" << endl;
    web << "<img src=\"./PICS/AB2_aP12_1_4a_8a_prim.png\" align=\"middle\" style=\"width:60%\" alt=\"Picture of Structure; Click for Big Picture\" />" << endl;
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
    web << "\\mathbf{B_{1}}&=&x_{1}\\,\\mathbf{a_1}+y_{1}\\,\\mathbf{a_2}+z_{1}\\,\\mathbf{a_3}&=&\\left(x_{1}\\,a+y_{1}\\,b\\,\\cos\\gamma\\,+z_{1}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{1}\\,b\\,\\sin\\gamma+z_{1}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{1}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{Fe I}\\\\" << endl;
    web << "\\mathbf{B_{2}}&=&x_{2}\\,\\mathbf{a_1}+y_{2}\\,\\mathbf{a_2}+z_{2}\\,\\mathbf{a_3}&=&\\left(x_{2}\\,a+y_{2}\\,b\\,\\cos\\gamma\\,+z_{2}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{2}\\,b\\,\\sin\\gamma+z_{2}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{2}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{Fe II}\\\\" << endl;
    web << "\\mathbf{B_{3}}&=&x_{3}\\,\\mathbf{a_1}+y_{3}\\,\\mathbf{a_2}+z_{3}\\,\\mathbf{a_3}&=&\\left(x_{3}\\,a+y_{3}\\,b\\,\\cos\\gamma\\,+z_{3}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{3}\\,b\\,\\sin\\gamma+z_{3}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{3}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{Fe III}\\\\" << endl;
    web << "\\mathbf{B_{4}}&=&x_{4}\\,\\mathbf{a_1}+y_{4}\\,\\mathbf{a_2}+z_{4}\\,\\mathbf{a_3}&=&\\left(x_{4}\\,a+y_{4}\\,b\\,\\cos\\gamma\\,+z_{4}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{4}\\,b\\,\\sin\\gamma+z_{4}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{4}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{Fe IV}\\\\" << endl;
    web << "\\mathbf{B_{5}}&=&x_{5}\\,\\mathbf{a_1}+y_{5}\\,\\mathbf{a_2}+z_{5}\\,\\mathbf{a_3}&=&\\left(x_{5}\\,a+y_{5}\\,b\\,\\cos\\gamma\\,+z_{5}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{5}\\,b\\,\\sin\\gamma+z_{5}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{5}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S I}\\\\" << endl;
    web << "\\mathbf{B_{6}}&=&x_{6}\\,\\mathbf{a_1}+y_{6}\\,\\mathbf{a_2}+z_{6}\\,\\mathbf{a_3}&=&\\left(x_{6}\\,a+y_{6}\\,b\\,\\cos\\gamma\\,+z_{6}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{6}\\,b\\,\\sin\\gamma+z_{6}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{6}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S II}\\\\" << endl;
    web << "\\mathbf{B_{7}}&=&x_{7}\\,\\mathbf{a_1}+y_{7}\\,\\mathbf{a_2}+z_{7}\\,\\mathbf{a_3}&=&\\left(x_{7}\\,a+y_{7}\\,b\\,\\cos\\gamma\\,+z_{7}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{7}\\,b\\,\\sin\\gamma+z_{7}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{7}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S III}\\\\" << endl;
    web << "\\mathbf{B_{8}}&=&x_{8}\\,\\mathbf{a_1}+y_{8}\\,\\mathbf{a_2}+z_{8}\\,\\mathbf{a_3}&=&\\left(x_{8}\\,a+y_{8}\\,b\\,\\cos\\gamma\\,+z_{8}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{8}\\,b\\,\\sin\\gamma+z_{8}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{8}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S IV}\\\\" << endl;
    web << "\\mathbf{B_{9}}&=&x_{9}\\,\\mathbf{a_1}+y_{9}\\,\\mathbf{a_2}+z_{9}\\,\\mathbf{a_3}&=&\\left(x_{9}\\,a+y_{9}\\,b\\,\\cos\\gamma\\,+z_{9}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{9}\\,b\\,\\sin\\gamma+z_{9}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{9}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S V}\\\\" << endl;
    web << "\\mathbf{B_{10}}&=&x_{10}\\,\\mathbf{a_1}+y_{10}\\,\\mathbf{a_2}+z_{10}\\,\\mathbf{a_3}&=&\\left(x_{10}\\,a+y_{10}\\,b\\,\\cos\\gamma\\,+z_{10}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{10}\\,b\\,\\sin\\gamma+z_{10}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{10}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S VI}\\\\" << endl;
    web << "\\mathbf{B_{11}}&=&x_{11}\\,\\mathbf{a_1}+y_{11}\\,\\mathbf{a_2}+z_{11}\\,\\mathbf{a_3}&=&\\left(x_{11}\\,a+y_{11}\\,b\\,\\cos\\gamma\\,+z_{11}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{11}\\,b\\,\\sin\\gamma+z_{11}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{11}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S VII}\\\\" << endl;
    web << "\\mathbf{B_{12}}&=&x_{12}\\,\\mathbf{a_1}+y_{12}\\,\\mathbf{a_2}+z_{12}\\,\\mathbf{a_3}&=&\\left(x_{12}\\,a+y_{12}\\,b\\,\\cos\\gamma\\,+z_{12}\\,c_x\\right)\\,\\mathbf{\\hat{x}}+\\left(y_{12}\\,b\\,\\sin\\gamma+z_{12}\\,c_y\\right)\\,\\mathbf{\\hat{y}}+z_{12}\\,c_z\\,\\mathbf{\\hat{z}}&\\left(1a\\right)&\\mbox{S VIII}\\\\" << endl;
    web << "\\end{array}" << endl;
    web << "\\]" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<h3>References</h3>" << endl;
    web << "<ul>" << endl;
    web << "<li> P. Bayliss, <em>Crystal structure refinement of a weakly anisotropic pyrite</em>, Am. Mineral. <strong>62</strong>, 1168&ndash;1172 (1977).</li>" << endl;
    web << "</ul>" << endl;
    web << "</div>" << endl;
    web << "</div>" << endl;
    web << "<hr />" << endl;
    web << "<div class=\"row\">" << endl;
    web << "<div class=\"col-md-12\">" << endl;
    web << "<h3>Geometry files</h3>" << endl;
    web << "<ul>" << endl;
    web << "<li><a href=\"./CIF/AB2_aP12_1_4a_8a.cif\">CIF</a></li>" << endl;
    web << "<li><a href=\"./POSCAR/POSCAR.AB2_aP12_1_4a_8a\">POSCAR</a></li>" << endl;
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

    if(LDEBUG) {cerr << "anrl:: WebANRL_AB2_aP12_1_4a_8a: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************