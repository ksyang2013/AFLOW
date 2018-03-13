// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo, Wahyu Setyawan - 2007-2011 Duke

#ifndef _AFLOW_LATTICE_CPP
#define _AFLOW_LATTICE_CPP

#include "aflow.h"

#define _RHL_HEX_SC_   // if you want RHL conventional to be the HEX
#define _EPS_ 0.02

// ***************************************************************************
namespace LATTICE {
  bool lattice_is_working(string lat) {
    if(lat=="CUB" || lat=="cP") return TRUE;
    if(lat=="BCC" || lat=="cI") return TRUE;
    if(lat=="FCC" || lat=="cF") return TRUE;
    if(lat=="TET" || lat=="tP") return TRUE;
    if(lat=="BCT" || lat=="tI") return TRUE;
    if(lat=="ORC" || lat=="oP") return TRUE;
    if(lat=="ORCF" || lat=="oF") return TRUE;
    if(lat=="ORCI" || lat=="oI") return TRUE;
    if(lat=="ORCC" || lat=="oS") return TRUE;
    if(lat=="HEX" || lat=="hP") return TRUE;
    if(lat=="RHL" || lat=="hR") return TRUE;
    if(lat=="MCL" || lat=="mP") return TRUE;
    if(lat=="MCLC" || lat=="mS") return TRUE;
    if(lat=="TRI" || lat=="aP") return TRUE;

    cerr << "ERROR: BZ for " << lat << " is not ready" << endl;
    exit(0); // you do not want to produce or run stuff that is not ready so you do not clog the computers
    return FALSE;
  }
}

// ***************************************************************************
namespace LATTICE {
  bool SpaceGroup2Lattice(uint sg,string &lattice_type,string& lattice_system) {
    /* FROM SAXENA BOOK
       1 TRI  (triclinic)
       2 MCL  (monoclinic)
       3 MCLC (C-centered monoclinic)
       4 ORC  (orthorhombic)
       5 ORCC (C-centered orthorhombic)
       6 ORCF (face-centered orthorhombic)
       7 ORCI (body-centered orthorhombic)
       8 TET  (tetragonal)
       9 BCT  (body-centered tetragonal)
       10 RHL  (rhombohedral/trigonal)
       11 HEX  (hexagonal)
       12 CUB  (simple cubic)
       13 FCC  (face-centered cubic)
       14 BCC  (body-centered cubic)
    */
    lattice_type="";lattice_system="";
    // tri
    if(sg==1 || sg==2) lattice_type="TRI";
    // mcl
    if(sg==3 || sg==4 || sg==6 ||sg==7 || sg==10 || sg==11 || sg==13 || sg==14) lattice_type="MCL";
    if(sg==5 || sg==8 || sg==9 || sg==12 || sg==15) lattice_type="MCLC";
    // orc
    if((sg>=16 && sg <=19) || (sg>=25 && sg<=34) || (sg>=47 && sg<=62)) lattice_type="ORC";
    if(sg==20 || sg==21 || (sg>=35 && sg<=41) || (sg>=63 && sg<=68)) lattice_type="ORCC";
    if(sg==22 || sg==42 || sg==43 ||sg==69 || sg==70) lattice_type="ORCF";
    if(sg==23 || sg==24 || (sg>=44 && sg<=46) || (sg>=71 && sg<=74)) lattice_type="ORCI";
    // tet
    if((sg>=75 && sg<=78) || sg==81 || (sg>=83 && sg<=86) || (sg>=89 && sg<=96) ||
       (sg>=99 && sg<=106) || (sg>=111 && sg<=118) || (sg>=123 && sg<=138)) lattice_type="TET";
    if(sg==79 || sg==80 || sg==82 || sg==87 || sg==88 || sg==97 || sg==98 ||
       (sg>=107 && sg<=110) || (sg>=119 && sg<=122) || (sg>=139 && sg<=142)) lattice_type="BCT";
    // hex
    if(sg>=143 && sg<=167) {
      if(sg==146 || sg==148 || sg==155 || sg==160 || sg==161 || sg==166 || sg==167) lattice_type="RHL";
      else lattice_type="HEX";
    }
    // http://en.wikipedia.org/wiki/Space_group
    if(sg>=168 && sg<=194) lattice_type="HEX";
    // cub
    if(sg==195 || sg==198 || sg==200 || sg==201 || sg==205 || sg==207 || sg==208 ||
       sg==212 || sg==213 || sg==215 || sg==218 || (sg>=221 && sg<=224)) lattice_type="CUB";
    if(sg==196 || sg==202 || sg==203 || sg==209 || sg==210 || sg==216 || sg==219 || (sg>=225 && sg<=228)) lattice_type="FCC";
    if(sg==197 || sg==199 || sg==204 || sg==206 || sg==211 || sg==214 || sg==217 || sg==220 || sg==229 || sg==230) lattice_type="BCC";
    // system
    if(lattice_type=="TRI") {lattice_system="TRI";return TRUE;}
    if(lattice_type=="MCL" || lattice_type=="MCLC") {lattice_system="MCL";return TRUE;}
    if(lattice_type=="ORC" || lattice_type=="ORCC" || lattice_type=="ORCF" || lattice_type=="ORCI") {lattice_system="ORC";return TRUE;}
    if(lattice_type=="TET" || lattice_type=="BCT") {lattice_system="TET";return TRUE;}
    if(lattice_type=="HEX" || lattice_type=="RHL") {lattice_system="HEX";return TRUE;}
    if(lattice_type=="CUB" || lattice_type=="FCC" || lattice_type=="BCC") {lattice_system="CUB";return TRUE;}

    lattice_type="error";lattice_system="error";
    return FALSE;
  }
}

namespace LATTICE {
  string SpaceGroup2Lattice(uint sg) {
    string lattice_type,lattice_system;
    LATTICE::SpaceGroup2Lattice(sg,lattice_type,lattice_system);
    return lattice_type;
  }
}

namespace LATTICE {
  string Lattice_System_SpaceGroup(uint sg) {
    string lattice_type,lattice_system;
    LATTICE::SpaceGroup2Lattice(sg,lattice_type,lattice_system);
    return lattice_system;
  }
}

// ***************************************************************************
namespace LATTICE {
  uint Lattice2SpaceGroup(string lattice,vector<uint>& vsgn) {
    vsgn.clear();
    string lattice_type,lattice_system;
    for(uint sg=1;sg<=230;sg++) {
      LATTICE::SpaceGroup2Lattice(sg,lattice_type,lattice_system);
      if(lattice_type=="error") {cerr << "LATTICE::Lattice2SpaceGroup error" << endl; exit(0);}
      if(lattice_type==lattice) vsgn.push_back(sg);
    }
    return vsgn.size();  // if zero, error
  }
}

// ***************************************************************************
namespace LATTICE {
  string SpaceGroup2LatticeVariation(uint sg,const xstructure& str_in) {
    /*
      This function return the bravais_lattice_variation_type according to
      spacegroup number sg and if necessary xstructure str. xstructure str is
      needed when there is more than one variation of the lattice, in this case,
      sg will be used as additional check of the lattice type.
    
      List of bravais_lattice_variation_type:
      1. TRI order: kalpha,kbeta,kgamma  > 90 (kgamma<kalpha, kgamma<kbeta)
      or kalpha,kbeta,kgamma  < 90 (kgamma>kalpha, kgamma>kbeta)
      special case when kgamma=90
      "TRI1a" kalpha>90 kbeta>90 kgamma>90
      "TRI1b" kalpha<90 kbeta<90 kgamma<90
      "TRI2a" kalpha>90 kbeta>90 kgamma=90
      "TRI2b" kalpha<90 kbeta<90 kgamma=90
      2. "MCL" unique (order b<=c)
      3. MCLC (order alpha<90)
      "MCLC1"  kgamma>90
      "MCLC2"  kgamma=90
      "MCLC3"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 < 1
      "MCLC4"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 = 1
      "MCLC5"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 > 1
      4. "ORC" unique (order a<b<c)
      5. "ORCC" unique (order a<b)
      6. ORCF (order a<b<c)
      "ORCF1" "ORCF_invb2+invc2<inva2"  for 1/a^2 > 1/b^2 + 1/c^2
      "ORCF2" "ORCF_inva2<invb2+invc2"  for 1/a^2 < 1/b^2 + 1/c^2
      "ORCF3"                           for 1/a^2 = 1/b^2 + 1/c^2
      7. "ORCI" unique (order a<b<c)
      8. "TET" unique (order a a c)
      9. BCT (order a a c)
      "BCT1" "BCT_c<a" for c<a
      "BCT2" "BCT_c>a" for c>a
      10. "RHL1" alpha<90
      "RHL2" alpha>90
      11. "HEX" unique (order 60 90 90)
      12. "CUB" unique
      13. "FCC" unique (order 60 60 60)
      14. "BCC" unique
    */
    
    bool found=false;
    string lattice="";
    xstructure str_sp, str_sc;
    // 1. TRI
    if(sg==1 || sg==2) {
      found=LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
      if(found) {
	lattice=str_sp.bravais_lattice_variation_type;
	if(aurostd::substring2bool(lattice,"TRI") || aurostd::substring2bool(lattice,"aP")) return lattice;
	else {
	  return "TRI";
	  cerr<<"WARNING LATTICE::Standard_Lattice_Structure found "<<lattice<<" instead of TRI for sg "<<sg<<endl;
	}
      }
      else return "TRI";
    }
    // 2. MCL
    if(sg==3 || sg==4 || sg==6 ||sg==7 || sg==10 || sg==11 || sg==13 || sg==14) return "MCL";
    // 3. MCLC
    if(sg==5 || sg==8 || sg==9 || sg==12 || sg==15) {
      found=LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
      if(found) {
	lattice=str_sp.bravais_lattice_variation_type;
	if(aurostd::substring2bool(lattice,"MCLC") || aurostd::substring2bool(lattice,"mS")) return lattice;
	else {
	  return "MCLC";
	  cerr<<"WARNING LATTICE::Standard_Lattice_Structure found "<<lattice<<" instead of MCLC for sg "<<sg<<endl;
	}
      }
      else return "MCLC";
    }
    // 4. ORC
    if((sg>=16 && sg <=19) || (sg>=25 && sg<=34) || (sg>=47 && sg<=62)) return "ORC";
    // 5. ORCC
    if(sg==20 || sg==21 || (sg>=35 && sg<=41) || (sg>=63 && sg<=68)) {
      return "ORCC";
    }
    // 6. ORCF
    if(sg==22 || sg==42 || sg==43 ||sg==69 || sg==70) {
      found=LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
      if(found) {
	lattice=str_sp.bravais_lattice_variation_type;
	if(aurostd::substring2bool(lattice,"ORCF") || aurostd::substring2bool(lattice,"oF")) return lattice;
	else {
	  return "ORCF";
	  cerr<<"WARNING LATTICE::Standard_Lattice_Structure found "<<lattice<<" instead of ORCF for sg "<<sg<<endl;
	}
      }
      else return "ORCF";
    }
    // 7. ORCI
    if(sg==23 || sg==24 || (sg>=44 && sg<=46) || (sg>=71 && sg<=74)) return "ORCI";
    // 8. TET
    if((sg>=75 && sg<=78) || sg==81 || (sg>=83 && sg<=86) || (sg>=89 && sg<=96) || (sg>=99 && sg<=106) || (sg>=111 && sg<=118) || (sg>=123 && sg<=138)) return "TET";
    // 9. BCT
    if(sg==79 || sg==80 || sg==82 || sg==87 || sg==88 || sg==97 || sg==98 || (sg>=107 && sg<=110) || (sg>=119 && sg<=122) || (sg>=139 && sg<=142)) {
      found=LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
      if(found) {
	lattice=str_sp.bravais_lattice_variation_type;
	if(aurostd::substring2bool(lattice,"BCT") || aurostd::substring2bool(lattice,"tI")) return lattice;
	else {
	  return "BCT";
	  cerr<<"WARNING LATTICE::Standard_Lattice_Structure found "<<lattice<<" instead of BCT for sg "<<sg<<endl;
	}
      }
      else return "BCT";
    }
    // 10. Trigonal, can be rhombohedral or hexagonal
    if(sg>=143 && sg<=167) {
      if(sg==146 || sg==148 || sg==155 || sg==160 || sg==161 || sg==166 || sg==167) {
	found=LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
	if(found) {
	  lattice=str_sp.bravais_lattice_variation_type;
	  if(aurostd::substring2bool(lattice,"RHL") || aurostd::substring2bool(lattice,"hR")) return lattice;
	  else {
	    return "RHL";
	    cerr<<"WARNING LATTICE::Standard_Lattice_Structure found "<<lattice<<" instead of RHL for sg "<<sg<<endl;
	  }
	}
	else return "RHL";
      }
      else return "HEX";
    }
    // 11. HEX
    if(sg>=168 && sg<=194) return "HEX";
    // 12. CUB
    if(sg==195 || sg==198 || sg==200 || sg==201 || sg==205 || sg==207 || sg==208 || sg==212 || sg==213 || sg==215 || sg==218 || (sg>=221 && sg<=224)) return "CUB";
    // 13. FCC
    if(sg==196 || sg==202 || sg==203 || sg==209 || sg==210 || sg==216 || sg==219 || (sg>=225 && sg<=228)) return "FCC";
    // 14. BCC
    if(sg==197 || sg==199 || sg==204 || sg==206 || sg==211 || sg==214 || sg==217 || sg==220 || sg==229 || sg==230) return "BCC";

    return "error";
  }
}

// ***************************************************************************
namespace LATTICE {
  string ConventionalLattice_SpaceGroup(uint sg,double a,double b,double c) {
    double eps=0.05;
    /*
    // 1. "TRI"
    // 2. "MCL"
    // 3. "MCLC"
    // 4. "ORC"
    // 5. "ORCC"
    // 6. "ORCF1" "ORCF_invb2+invc2<inva2"
    //    "ORCF2" "ORCF_inva2<invb2+invc2"
    //    "ORCF3" "ORCF_inva2=invb2+invc2"
    // 7. "ORCI"
    // 8. "TET"
    // 9. "BCT1"  "BCT_c<a"
    //    "BCT2"  "BCT_a<c"
    // 10. "RHL"
    // 11. "HEX"
    // 12. "CUB"
    // 13. "FCC"
    // 14. "BCC"
    */
  
    // 1. TRI
    // there are 4 variations depending on kalpha,kbeta,kgamma.
    // Since we don't have the lattice vectors as input parameters of this function,
    // we will return just "TRI" for now.
    if(sg==1 || sg==2) return "TRI";
    // 2. MCL
    if(sg==3 || sg==4 || sg==6 ||sg==7 || sg==10 || sg==11 || sg==13 || sg==14) return "MCL";
    // 3. MCLC
    if(sg==5 || sg==8 || sg==9 || sg==12 || sg==15) {
      //There are 5 variations of MCLC depending on alpha, kgamma, and some testphi
      //However, this function needs lattice vectors to determine alpha,kgamma,testphi
      //and to return the proper variation, for now it will return "MCLC".
      return "MCLC";
      cerr << "LATTICE::ConventionalLattice_SpaceGroup: error in MCLCC" << endl;
      cerr << a << " " << b << " " << c << " " << endl << endl;
      exit(0);
    }
    // 4. ORC
    if((sg>=16 && sg <=19) || (sg>=25 && sg<=34) || (sg>=47 && sg<=62)) return "ORC";
    // 5. ORCC
    if(sg==20 || sg==21 || (sg>=35 && sg<=41) || (sg>=63 && sg<=68)) {
      return "ORCC";
      cerr << "LATTICE::ConventionalLattice_SpaceGroup: error in ORCC" << endl;
      cerr << a << " " << b << " " << c << " " << endl << endl;
      exit(0);
    }
    // 6. ORCF
    if(sg==22 || sg==42 || sg==43 ||sg==69 || sg==70) {
      double mismatch=pow(1.0/a,2.0)/(pow(1.0/b,2.0)+pow(1.0/c,2.0));
      // cerr << mismatch << " " << eps/5.0 << endl;  // DEBUG purposes
      if(mismatch>1+eps/5.0) return "ORCF1";
      if(mismatch<1-eps/5.0) return "ORCF2";
      if(aurostd::isequal(mismatch,1.0,eps/5.0)) return "ORCF3";
      cerr << "LATTICE::ConventionalLattice_SpaceGroup: error in ORCF" << endl;
      cerr << a << " " << b << " " << c << " " << endl << endl;
      exit(0);
    }
    // 7. ORCI
    if(sg==23 || sg==24 || (sg>=44 && sg<=46) || (sg>=71 && sg<=74)) return "ORCI";
    // 8. TET
    if((sg>=75 && sg<=78) || sg==81 || (sg>=83 && sg<=86) || (sg>=89 && sg<=96) || (sg>=99 && sg<=106) || (sg>=111 && sg<=118) || (sg>=123 && sg<=138)) return "TET";
    // 9. BCT
    if(sg==79 || sg==80 || sg==82 || sg==87 || sg==88 || sg==97 || sg==98 || (sg>=107 && sg<=110) || (sg>=119 && sg<=122) || (sg>=139 && sg<=142)) {
      if(c<=a) return "BCT1";
      if(a<c)  return "BCT2";
      // test if a=b=c
      if(aurostd::isequal(a,b,eps) && aurostd::isequal(a,c,eps)) return "BCT2";
      cerr << "LATTICE::ConventionalLattice_SpaceGroup: error in BCT" << endl;
      cerr << a << " " << b << " " << c << " " << endl << endl;
      exit(0);
    }
    // 10. Trigonal, can be rhombohedral or hexagonal
    if(sg>=143 && sg<=167) {
      if(sg==146 || sg==148 || sg==155 || sg==160 || sg==161 || sg==166 || sg==167) return "RHL";
      else return "HEX";
    }
    // 11. HEX
    if(sg>=168 && sg<=194) return "HEX";
    // 12. CUB
    if(sg==195 || sg==198 || sg==200 || sg==201 || sg==205 || sg==207 || sg==208 || sg==212 || sg==213 || sg==215 || sg==218 || (sg>=221 && sg<=224)) return "CUB";
    // 13. FCC
    if(sg==196 || sg==202 || sg==203 || sg==209 || sg==210 || sg==216 || sg==219 || (sg>=225 && sg<=228)) return "FCC";
    // 14. BCC
    if(sg==197 || sg==199 || sg==204 || sg==206 || sg==211 || sg==214 || sg==217 || sg==220 || sg==229 || sg==230) return "BCC";

    return "error";
  }
}

// ***************************************************************************
namespace LATTICE {
  string ConventionalLattice_SpaceGroup(uint sg,const xstructure& str) {
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES); // always
    return LATTICE::ConventionalLattice_SpaceGroup(sg,data(1),data(2),data(3));
  }
}

// ***************************************************************************
// stantard conventional vs standard primitive operations
namespace LATTICE {
  xmatrix<double> sc2sp(const xmatrix<double>& rlattice, string lattice,bool inverseflag) {
    // inverseflag==0   take rlattice as standard conventional and output standard primitive
    // inverseflag==1   take rlattice as standard primitive and output standard conventional
    //
    // lattice: bravais_lattice_type or bravais_lattice_variation_type
    if(inverseflag==0) inverseflag=1;
    else inverseflag=0;
    return(LATTICE::sp2sc(rlattice,lattice,inverseflag));
  }
}

// ***************************************************************************
// stantard conventional vs standard primitive operations
namespace LATTICE {
  xmatrix<double> sp2sc(const xmatrix<double>& rlattice, string lattice,bool inverseflag) {
    // inverseflag==0   take rlattice as standard primitive and output standard conventional
    // inverseflag==1   take rlattice as standard conventional and output standard primitive
    //
    // lattice: bravais_lattice_type or bravais_lattice_variation_type
    
    bool found;
    xmatrix<double> clattice(3,3), identity(3,3), M(3,3);
    
    identity[1][1]=1.0;identity[1][2]=0.0;identity[1][3]=0.0;
    identity[2][1]=0.0;identity[2][2]=1.0;identity[2][3]=0.0;
    identity[3][1]=0.0;identity[3][2]=0.0;identity[3][3]=1.0;
    
    found=false;
    if(lattice=="cub" || lattice=="CUB" || lattice=="cP") {// unique choice
      found=true; M=identity;
    }
    if(lattice=="bcc" || lattice=="BCC" || lattice=="cI") {// unique choice
      found=true;
      xmatrix<double> bcc(3,3);
      bcc[1][1]=0.0;bcc[1][2]=1.0;bcc[1][3]=1.0;
      bcc[2][1]=1.0;bcc[2][2]=0.0;bcc[2][3]=1.0;
      bcc[3][1]=1.0;bcc[3][2]=1.0;bcc[3][3]=0.0;
      M=bcc;
    }
    if(lattice=="fcc" || lattice=="FCC" || lattice=="cF") {// unique choice
      found=true;
      xmatrix<double> fcc(3,3);
      fcc[1][1]=-1.0;fcc[1][2]= 1.0;fcc[1][3]= 1.0;
      fcc[2][1]= 1.0;fcc[2][2]=-1.0;fcc[2][3]= 1.0;
      fcc[3][1]= 1.0;fcc[3][2]= 1.0;fcc[3][3]=-1.0;
      M=fcc;
    }
    if(lattice=="tet" || lattice=="TET" || lattice=="tP") {
      found=true; M=identity;
    }
    if(aurostd::substring2bool(lattice,"bct") || aurostd::substring2bool(lattice,"BCT")
       || aurostd::substring2bool(lattice,"tI")) {
      found=true;
      xmatrix<double> bct(3,3);
      bct[1][1]=0.0;bct[1][2]=1.0;bct[1][3]=1.0;
      bct[2][1]=1.0;bct[2][2]=0.0;bct[2][3]=1.0;
      bct[3][1]=1.0;bct[3][2]=1.0;bct[3][3]=0.0;
      M=bct;
    }  
    if(lattice=="orc" || lattice=="ORC" || lattice=="oP") {// unique choice
      found=true; M=identity;
    }  
    if(aurostd::substring2bool(lattice,"orcf") || aurostd::substring2bool(lattice,"ORCF")
       || aurostd::substring2bool(lattice,"oF")) {
      found=true;
      xmatrix<double> orcf(3,3);
      orcf[1][1]=-1.0;orcf[1][2]= 1.0;orcf[1][3]= 1.0;
      orcf[2][1]= 1.0;orcf[2][2]=-1.0;orcf[2][3]= 1.0;
      orcf[3][1]= 1.0;orcf[3][2]= 1.0;orcf[3][3]=-1.0;
      M=orcf;
    }  
    if(aurostd::substring2bool(lattice,"orci") || aurostd::substring2bool(lattice,"ORCI")
       || aurostd::substring2bool(lattice,"oI")) {
      found=true;
      xmatrix<double> orci(3,3);
      orci[1][1]=0.0;orci[1][2]=1.0;orci[1][3]=1.0;
      orci[2][1]=1.0;orci[2][2]=0.0;orci[2][3]=1.0;
      orci[3][1]=1.0;orci[3][2]=1.0;orci[3][3]=0.0;
      M=orci;
    }  
    if(aurostd::substring2bool(lattice,"orcc") || aurostd::substring2bool(lattice,"ORCC")
       || aurostd::substring2bool(lattice,"oS") || aurostd::substring2bool(lattice,"oC")) {
      found=true;
      xmatrix<double> orcc(3,3);
      orcc[1][1]=1.0;orcc[1][2]=1.0;orcc[1][3]=0.0;
      orcc[2][1]=-1.0;orcc[2][2]=1.0;orcc[2][3]=0.0;
      orcc[3][1]=0.0;orcc[3][2]=0.0;orcc[3][3]=1.0;
      M=orcc;
    }
    if(lattice=="hex" || lattice=="HEX" || lattice=="hP") {
      found=true; M=identity;
    }
    if(aurostd::substring2bool(lattice,"rhl") || aurostd::substring2bool(lattice,"RHL")
       || aurostd::substring2bool(lattice,"hR")) {
      found=true; M=identity;
    }
    if(lattice=="mcl" || lattice=="MCL" || lattice=="mP") {// unique choice
      found=true; M=identity;
    }
    if(aurostd::substring2bool(lattice,"mclc") || aurostd::substring2bool(lattice,"MCLC")
       || aurostd::substring2bool(lattice,"mS") || aurostd::substring2bool(lattice,"mC")) {
      found=true;
      xmatrix<double> mclc(3,3);
      mclc[1][1]=1.0;mclc[1][2]=-1.0;mclc[1][3]=0.0;
      mclc[2][1]=1.0;mclc[2][2]=1.0;mclc[2][3]=0.0;
      mclc[3][1]=0.0;mclc[3][2]=0.0;mclc[3][3]=1.0;
      M=mclc;
    }
    if(aurostd::substring2bool(lattice,"tri") || aurostd::substring2bool(lattice,"TRI") || aurostd::substring2bool(lattice,"aP")) {// unique choice
      found=true; M=identity;
    }
    if(found) {
      if(inverseflag==0) clattice=M*rlattice;
      else clattice=inverse(M)*rlattice;
      return clattice;
    }
    else {
      cerr<<"ERROR: lattice "<<lattice<<" not found in function LATTICE::sp2sc."<<endl; abort();
    }
  }
}

// ***************************************************************************
namespace LATTICE {
  xvector<double> Getabc_angles_Conventional(const xmatrix<double>& rlattice, string lattice,const int& mode) {
    xmatrix<double> clattice(3,3);
    clattice=rlattice;
    // lattices
    xmatrix<double> identity(3,3);
    identity[1][1]=1.0;identity[1][2]=0.0;identity[1][3]=0.0;
    identity[2][1]=0.0;identity[2][2]=1.0;identity[2][3]=0.0;
    identity[3][1]=0.0;identity[3][2]=0.0;identity[3][3]=1.0;
  
    if(lattice=="cub" || lattice=="CUB") {// unique choice
      xmatrix<double> cub(3,3);cub=identity;
      clattice=cub*rlattice;
    }
    if(lattice=="bcc" || lattice=="BCC") {// unique choice
      xmatrix<double> bcc(3,3);
      bcc[1][1]=0.0;bcc[1][2]=1.0;bcc[1][3]=1.0;
      bcc[2][1]=1.0;bcc[2][2]=0.0;bcc[2][3]=1.0;
      bcc[3][1]=1.0;bcc[3][2]=1.0;bcc[3][3]=0.0;
      clattice=bcc*rlattice;
    }
    if(lattice=="fcc" || lattice=="FCC") {// unique choice
      xmatrix<double> fcc(3,3);
      fcc[1][1]=-1.0;fcc[1][2]= 1.0;fcc[1][3]= 1.0;
      fcc[2][1]= 1.0;fcc[2][2]=-1.0;fcc[2][3]= 1.0;
      fcc[3][1]= 1.0;fcc[3][2]= 1.0;fcc[3][3]=-1.0;
      clattice=fcc*rlattice;
    }
    if(aurostd::substring2bool(lattice,"tet") || aurostd::substring2bool(lattice,"TET")) {
      xmatrix<double> tet(3,3);tet=identity;
      clattice=tet*rlattice;
    }
    if(aurostd::substring2bool(lattice,"bct") || aurostd::substring2bool(lattice,"BCT")) {
      xmatrix<double> bct(3,3);
      bct[1][1]=0.0;bct[1][2]=1.0;bct[1][3]=1.0;
      bct[2][1]=1.0;bct[2][2]=0.0;bct[2][3]=1.0;
      bct[3][1]=1.0;bct[3][2]=1.0;bct[3][3]=0.0;
      clattice=bct*rlattice;
    }
  
    if(aurostd::substring2bool(lattice,"rhl") || aurostd::substring2bool(lattice,"RHL")) {
      // xmatrix<double> rhl(3,3);
      // rhl[1][1]=1.0;rhl[1][2]=-1.0;rhl[1][3]= 0.0;
      // rhl[2][1]=0.0;rhl[2][2]= 1.0;rhl[2][3]=-1.0;
      // rhl[3][1]=1.0;rhl[3][2]= 1.0;rhl[3][3]= 1.0;
      // clattice=rhl*rlattice;
      xmatrix<double> rhl(3,3);rhl=identity;
      clattice=rhl*rlattice;
    }
  
    if(lattice=="hex" || lattice=="HEX") {// unique choice
      xmatrix<double> hex(3,3);hex=identity;
      clattice=hex*rlattice;
    }
  
    if(lattice=="orc" || lattice=="ORC") {// unique choice
      xmatrix<double> orc(3,3);orc=identity;
      clattice=orc*rlattice;
    }
  
    if(aurostd::substring2bool(lattice,"orcf") || aurostd::substring2bool(lattice,"ORCF")) {
      xmatrix<double> orcf(3,3);
      orcf[1][1]=-1.0;orcf[1][2]= 1.0;orcf[1][3]= 1.0;
      orcf[2][1]= 1.0;orcf[2][2]=-1.0;orcf[2][3]= 1.0;
      orcf[3][1]= 1.0;orcf[3][2]= 1.0;orcf[3][3]=-1.0;
      clattice=orcf*rlattice;
    }
  
    if(aurostd::substring2bool(lattice,"orci") || aurostd::substring2bool(lattice,"ORCI")) {
      xmatrix<double> orci(3,3);
      orci[1][1]=0.0;orci[1][2]=1.0;orci[1][3]=1.0;
      orci[2][1]=1.0;orci[2][2]=0.0;orci[2][3]=1.0;
      orci[3][1]=1.0;orci[3][2]=1.0;orci[3][3]=0.0;
      clattice=orci*rlattice;
    }
  
    if(aurostd::substring2bool(lattice,"orcc") || aurostd::substring2bool(lattice,"ORCC")) {
      xmatrix<double> orcc(3,3);
      orcc[1][1]=1.0;orcc[1][2]=1.0;orcc[1][3]=0.0;
      orcc[2][1]=-1.0;orcc[2][2]=1.0;orcc[2][3]=0.0;
      orcc[3][1]=0.0;orcc[3][2]=0.0;orcc[3][3]=1.0;
      clattice=orcc*rlattice;
    }
  
    if(lattice=="mcl" || lattice=="MCL") {// unique choice
      xmatrix<double> mcl(3,3);mcl=identity;
      clattice=mcl*rlattice;
    }

    if(aurostd::substring2bool(lattice,"mclc") || aurostd::substring2bool(lattice,"MCLC")) {
      xmatrix<double> mclc(3,3);
      mclc[1][1]=1.0;mclc[1][2]=-1.0;mclc[1][3]=0.0;
      mclc[2][1]=1.0;mclc[2][2]=1.0;mclc[2][3]=0.0;
      mclc[3][1]=0.0;mclc[3][2]=0.0;mclc[3][3]=1.0;
      clattice=mclc*rlattice;
    }
  
    if(lattice=="tri" || lattice=="TRI") {// unique choice
      xmatrix<double> tri(3,3);tri=identity;
      clattice=tri*rlattice;
    }
    // done with the lattices
  
    xvector<double> data(6);
    data(1)=modulus(clattice(1));
    data(2)=modulus(clattice(2));
    data(3)=modulus(clattice(3));
    data(4)=angle(clattice(2),clattice(3));
    data(5)=angle(clattice(3),clattice(1));
    data(6)=angle(clattice(1),clattice(2));
    if(mode==DEGREES) {
      data(4)*=rad2deg;
      data(5)*=rad2deg;
      data(6)*=rad2deg;
    }
    return data;
  }
}

// ***************************************************************************
namespace LATTICE {
  bool fix_sts_sp(xstructure& str_sp,xmatrix<double> &rlattice,xmatrix<double> &plattice) {
    // this routine sets plattice as the str_sp.lattice, and fixes the atoms
    // plattice is a rotated version of rlattice.
    str_sp.lattice=rlattice;
    str_sp.FixLattices();
    // cpos are identical in the old and new lattice (the new rlattice is just redefined)
    for(uint i=0;i<str_sp.atoms.size();i++)
      str_sp.atoms.at(i).fpos=C2F(str_sp.lattice,str_sp.atoms.at(i).cpos); // DONE, now it has new fpos with the same cpos !
    str_sp.BringInCell();
    // now the new fpos are correct.  make new plattice which is a rotation
    str_sp.lattice=plattice;
    str_sp.FixLattices();
    // fpos are correct with respect to the new basis (it has the right angles, lenghts) so get the cpos
    for(uint i=0;i<str_sp.atoms.size();i++)
      str_sp.atoms.at(i).cpos=F2C(str_sp.lattice,str_sp.atoms.at(i).fpos); // DONE, now it has new cpos with the same fpos !
    str_sp.BringInCell();
    
    return TRUE;
  }
}

// ***************************************************************************
// DX - START
namespace LATTICE {
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc, bool full_sym) {
    return LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc,full_sym);
  }
}
/* OLD DX
// ***************************************************************************
// WORKING WORKING WORKING WORKING WORKING WORKING WORKING WORKING WORKING WORKING
namespace LATTICE {
bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc) {
return LATTICE::Standard_Lattice_StructureNormal(str_in,str_sp,str_sc);}
}
*/
// DX - END

//**************************JUNKAI EDITED START****************************
namespace LATTICE {
  bool Standard_Lattice_StructureDefault(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym) {
    //DX OBSOLETEreturn LATTICE::Standard_Lattice_StructureMedium(str_in,str_sp,str_sc);}
    return LATTICE::Standard_Lattice_Structure_20170718(str_in,str_sp,str_sc,full_sym);} //DX
/*
  bool Standard_Lattice_StructureCoarse(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc) {
    int time=0;
    return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.05,0.5,time,_EPS_);}
  //  return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.1,1.0,time,_EPS_);}
  bool Standard_Lattice_StructureNormal(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc) {
    int time=0;
    return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.05,0.5,time,_EPS_);}
  bool Standard_Lattice_StructureMedium(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc) {
    int time=0;
    return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.01,0.1,time,_EPS_);}
  bool Standard_Lattice_StructurePrecise(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc) {
    int time=0;
    return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.002,0.02,time,_EPS_);}
  bool Standard_Lattice_StructureUltra(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc) {
    int time=0;
    return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.0004,0.004,time,_EPS_);} // epsabc //epsang
*/
}
//**************************JUNKAI EDITED END****************************

namespace LATTICE {
  //bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps,bool histogram) ;
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps) {
    return LATTICE::Standard_Lattice_Structure_20170101(str_in,str_sp,str_sc,eps,epsang,time,symeps,false);
  }
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps,bool histogram) {
    return LATTICE::Standard_Lattice_Structure_20170101(str_in,str_sp,str_sc,eps,epsang,time,symeps,histogram);
  }
}

namespace LATTICE {
  // transformation to standardization
  bool Standard_Lattice_Structure_20170718(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym) {//**************************JUNKAI EDITED*****************
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    double eps = 0.001;
    if(str_in.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
      eps=str_in.sym_eps;
    }

    bool is_deg = true; //DX
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: BEGIN" << endl;
    // cerr << "eps=" << eps << " " << "epsang=" << epsang << endl;
    //  LDEBUG=TRUE;
    // starting from a str_in (whatever lattices) this routine returns a standard primitive (str_sp)
    // and a standard conventional (str_sc) following the ideas of stefano and wahyu
    if(str_in.Standard_Lattice_avoid==TRUE) return FALSE; // if you do not want to calculate
    // preparation
    str_sp=str_in; // copy it
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [1]" << endl;
    str_sp.GetPrimitive(0.005);
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [2]" << endl;
    str_sp.FixLattices();

    double str_sp_volume=str_sp.Volume();    // backups
    bool str_sp_neg_scale=str_sc.neg_scale;  // backups

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [3]" << endl;
    
    str_sp.ReScale(1.0);     // get it off the way, I might save and plug it back but it is not necessary now
    str_sp.neg_scale=FALSE;  // get it off the way, I might save and plug it back but it is not necessary now
    {str_sp.MinkowskiBasisReduction();} // shorten the vectors as much as possible and as perpendicular as possible

    str_sc=str_in; // copy it
    str_sc.FixLattices();
    str_sc.ReScale(1.0);     // get it off the way, I might save and plug it back but it is not necessary now
    str_sc.neg_scale=FALSE;  // get it off the way, I might save and plug it back but it is not necessary now
    {str_sc.MinkowskiBasisReduction();} // shorten the vectors as much as possible and as perpendicular as possible
    str_sc.ClearSymmetry();
  
    string crystal_system="none";
    bool SYS_VERBOSE=FALSE;
    //DX - START
    //DX 8/14/17 - if(str_sp.pgroup_calculated==FALSE || str_sp.fgroup_calculated==FALSE || str_sp.pgroup_xtal_calculated==FALSE) {
    if(str_sp.crystal_system=="") {
      ofstream FileDevNull("/dev/null");
      _aflags aflags;
      _kflags kflags;
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE=FALSE; //DX 12/7/17 - was mising before
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE=FALSE; //DX 12/7/17 - added pgroupk_xtal
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
      aflags.Directory="./";aflags.QUIET=FALSE;
      str_sp.LatticeReduction_avoid=TRUE;
      if(LDEBUG){ cerr << "LATTICE::Standard_Lattice_Structure: [4b]" << endl;}
      //pflow::PerformFullSymmetry(str_sp,FileDevNull,aflags,kflags,SYS_VERBOSE,cout);
      string format="out";
      if(str_sp.pgroup_xtal_calculated==TRUE) { //DX 8/14/17 - Speed increase, if already calculated just look use look-up table
        SYM::PointGroupLookUpTable(FileDevNull,str_sp,aflags,kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE,SYS_VERBOSE,cout,format); //DX 8/14/17
      }
      else if(str_sp.pgroup_xtal_calculated==FALSE) {
        //DX Calculate up to pgroup_xtal only
        if(full_sym){
          kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;
          kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;
          kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;
          kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=TRUE;
          kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=TRUE; //DX 12/7/17 - added pgroupk_xtal
          kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=TRUE;
          kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=TRUE;
          kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=TRUE;
        }
        else{
        kflags.KBIN_SYMMETRY_CALCULATE_PGROUP=TRUE;
        kflags.KBIN_SYMMETRY_CALCULATE_FGROUP=TRUE;
        kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL=TRUE;
        kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK=FALSE;
        kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL=FALSE; //DX 12/7/17 - added pgroupk_xtal //DX 12/18/17 - missing the "K" in "PGROUPK_XTAL"
        kflags.KBIN_SYMMETRY_CALCULATE_IATOMS=FALSE;
        kflags.KBIN_SYMMETRY_CALCULATE_AGROUP=FALSE;
        kflags.KBIN_SYMMETRY_CALCULATE_SGROUP=FALSE;
        }
      pflow::PerformFullSymmetry(str_sp,FileDevNull,aflags,kflags,SYS_VERBOSE,cout);
    }
    }
    crystal_system=str_sp.crystal_system;

    eps = str_sp.sym_eps;

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4]" << endl;

    xvector<int> dims(3);    // to scan
    dims=LatticeDimensionSphere(str_sp.lattice,RadiusSphereLattice(str_sp.lattice)/1.0);
    //  dims(1)=4;dims(2)=4;dims(3)=4;
    xmatrix<double> rlattice(3,3),plattice(3,3),clattice(3,3),klattice(3,3),t(3,3);
    xvector<double> aus(3),a1(3),a2(3),a3(3),rdata(6),kdata(6);
    double volume,vaus,nndist,cutoff,a,b,c,ka,kb,kc,alpha,beta,gamma,ac,bc,cc,minabc,maxabc,kalpha,kbeta,kgamma;
    // generate all the possible transformations
    bool found=FALSE;
    vector<xmatrix<double> > vrlattice1,vrlattice1_aaa;
    vector<xvector<double> > vvectors;
    double amin=min(aurostd::modulus(str_sp.lattice(1)),aurostd::modulus(str_sp.lattice(2)),aurostd::modulus(str_sp.lattice(3))); //DX
    volume=det(str_sp.lattice);
    double eps_volume = eps*amin*amin*amin; //DX
    vrlattice1.push_back(str_sp.lattice);
    
    nndist=1.0e6;
    cutoff=1.0e6;

    // find NN distance
    for(int i=-dims[1];i<=dims[1];i++)
      for(int j=-dims[2];j<=dims[2];j++)
	for(int k=-dims[3];k<=dims[3];k++) {
	  aus=i*str_sp.lattice(1)+j*str_sp.lattice(2)+k*str_sp.lattice(3);
	  if(modulus(aus)>eps && modulus(aus)<=nndist) nndist=modulus(aus);
	}

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4a]" << endl;


    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4c]" << endl;

    if(LDEBUG) {
      cerr << "str_sp.atoms.size()=" << str_sp.atoms.size() << endl;
      cerr << "str_sc.atoms.size()=" << str_sc.atoms.size() << endl;
      cerr << "str_sp.pgroup.size()=" << str_sp.pgroup.size() << endl;
      cerr << "str_sp.fgroup.size()=" << str_sp.fgroup.size() << endl;
      cerr << "str_sp.pgroup_xtal.size()=" << str_sp.pgroup_xtal.size() << endl;
      cerr << "str_sp.crystal_family=" << str_sp.crystal_family << endl;
      cerr << "str_sp.crystal_system=" << str_sp.crystal_system << endl;
      cerr << "crystal_system=" << crystal_system << endl;
      cerr << "radius_sp=" << RadiusSphereLattice(str_sp.lattice) << endl;
      cerr << "radius_sc=" << RadiusSphereLattice(str_sc.lattice) << endl;
      //  exit(0);
    }
  
    cutoff=0.0;
    rdata=Getabc_angles(str_sc.lattice,DEGREES);
    if(crystal_system=="cubic") {cutoff=sqrt(2)*max(rdata[1],rdata[2],rdata[3]);}
    if(crystal_system=="hexagonal") {cutoff=max(rdata[1],rdata[2],rdata[3]);} // the primitives always get the max of the lenghts
    //  if(crystal_system=="tetragonal") {cutoff=(2.0/sqrt(3.0))*max(rdata[1],rdata[2],rdata[3]);}
    if(crystal_system=="tetragonal") {cutoff=1.05*(2.0/sqrt(3.0))*RadiusSphereLattice(str_sc.lattice);}
    if(crystal_system=="trigonal") {cutoff=max(rdata[1],rdata[2],rdata[3]);} // the primitives always get the max of the lenghts
    //  if(crystal_system=="orthorhombic") {cutoff=sqrt(2)*max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;} //stefano recipe
    if(crystal_system=="orthorhombic") {cutoff=1.05*(2.0/sqrt(3.0))*RadiusSphereLattice(str_sc.lattice);}
    if(crystal_system=="monoclinic") {cutoff=sqrt(2)*max(rdata[1],rdata[2],rdata[3]);} //sqrt(2) catches the centered face
    if(crystal_system=="triclinic") {cutoff=max(rdata[1],rdata[2],rdata[3]);} // the primitives always get the max of the lenghts
    if(crystal_system=="triclinic") {cutoff=max(rdata[2]+rdata[3],rdata[1]+rdata[3],rdata[1]+rdata[2]);} // the primitives always get the max of the lenghts
    //  if(crystal_system=="triclinic") {cutoff=rdata[1]+rdata[2]+rdata[3];} // the primitives always get the max of the lenghts
  
    // cerr << crystal_system << endl;//exit(0);

    for(int i=-dims[1];i<=dims[1];i+=1)
      for(int j=-dims[2];j<=dims[2];j+=1)
	for(int k=-dims[3];k<=dims[3];k+=1) {
	  aus=i*str_sp.lattice(1)+j*str_sp.lattice(2)+k*str_sp.lattice(3);
	  if(modulus(aus)<=cutoff+0.01) vvectors.push_back(aus);   // cubics are all in nn distances
	}
  
    // most of the time is spent inside this routine
    for(uint i=0;i<vvectors.size();i++) { for(uint ii=1;ii<=3;ii++) rlattice[1][ii]=vvectors[i][ii];
      for(uint j=0;j<vvectors.size();j++)  { for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=vvectors[j][ii];
	for(uint k=0;k<vvectors.size();k++)  { for(uint ii=1;ii<=3;ii++) rlattice[3][ii]=vvectors[k][ii];
	  // vaus=det(vvectors[i],vvectors[j],vvectors[k]);                                                 // SLOW
	  vaus=rlattice[1][1]*rlattice[2][2]*rlattice[3][3]+rlattice[1][2]*rlattice[2][3]*rlattice[3][1]+        // FAST
	    rlattice[1][3]*rlattice[2][1]*rlattice[3][2]-rlattice[1][3]*rlattice[2][2]*rlattice[3][1]-        // FAST
	    rlattice[1][2]*rlattice[2][1]*rlattice[3][3]-rlattice[1][1]*rlattice[2][3]*rlattice[3][2];        // FAST
	  if(abs(vaus)>eps_volume){ //DX 9/4/17 - Otherwise, Getabc_angles may fail if vaus=0.0
	  if(abs(vaus-1.0*volume)<eps_volume || abs(vaus-2.0*volume)<eps_volume || abs(vaus-4.0*volume)<eps_volume) { //DX
	    rdata=Getabc_angles(rlattice,DEGREES);a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    //DX if(abs(vaus-volume)<eps) {
	    if(abs(vaus-volume)<eps_volume) { //DX
	      vrlattice1.push_back(rlattice);
	        if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps)){
		vrlattice1_aaa.push_back(rlattice);
                }
	      } //DX 9/4/17
	    }
	  }
	}
      }
    }
    if(LDEBUG || 0) {
      cerr << "DEBUG dims=" << dims << endl;
      cerr << "DEBUG vvectors.size()=" << vvectors.size() << endl;
      cerr << "DEBUG vrlattice1.size()=" << vrlattice1.size() << "   vrlattice1_aaa.size()=" << vrlattice1_aaa.size() << " " << endl;
      cerr << "DEBUG str_sp.pgroup.size()=" << str_sp.pgroup.size() << endl;
      cerr << "DEBUG str_sc.pgroup.size()=" << str_sc.pgroup.size() << endl;
      cerr << "DEBUG str_sp.fgroup.size()=" << str_sp.fgroup.size() << endl;
      cerr << "DEBUG str_sc.fgroup.size()=" << str_sc.fgroup.size() << endl;
      cerr << "DEBUG crystal_system=" << crystal_system << endl;
    }
  

    // ***************************************************************************
    // start scan
    found=FALSE;

    // ***************************************************************************************************
    bool VERBOSE_PROGRESS=FALSE; // TRUE;
    //  bool VERBOSE_PROGRESS=TRUE;
    // try different solutions
    for (uint choice=0;choice<=2&&!found;choice++) {
      if(VERBOSE_PROGRESS) if(choice==0) cerr << "LATTICE::Standard_Lattice_Structure: DEFINITION PRISTINE" << endl;
      if(VERBOSE_PROGRESS) if(choice==1) cerr << "LATTICE::Standard_Lattice_Structure: DEFINITION RELAX1" << endl; // locura
      if(VERBOSE_PROGRESS) if(choice==2) cerr << "LATTICE::Standard_Lattice_Structure: DEFINITION RELAX2" << endl; // locura
      if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: crystal_system=" << crystal_system << endl;
      
      // ***************************************************************************************************
      // CUBIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="cubic")) {
	// --------------------------------------------------------------------------
	// CUB FROM PRIMITIVE
	// CUBIC SYSTEMS looking for lattices a a a 90 90 90
	if(VERBOSE_PROGRESS) cerr << "CUBIC SYSTEMS looking for lattices a a a 90 90 90" << endl;
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                           // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if((choice==0 && (aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps)))) { // just to help not necessary
	    // CUB: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 a]: a a a 90 90 90  primitive
	    // CUB: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 a]: a a a 90 90 90  conventional
	    //DX if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps) &&
	    //DX   aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps) && //DX
	       SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX 
	      found=TRUE;
	      // do the STANDARD PRIMITIVE
	      str_sp.bravais_lattice_type="CUB";
	      str_sp.bravais_lattice_variation_type="CUB";
	      str_sp.bravais_lattice_system="cubic";
	      str_sp.pearson_symbol="cP"+aurostd::utype2string(str_sp.atoms.size());
	      ac=a; // make the standard primitive cubic
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive cubic
	      plattice[2][1]=0.0;plattice[2][2]= ac;plattice[2][3]=0.0; // make the standard primitive cubic
	      plattice[3][1]=0.0;plattice[3][2]=0.0;plattice[3][3]= ac; // make the standard primitive cubic
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // CUB standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE CUB
	// --------------------------------------------------------------------------
	// BCC FROM PRIMITIVE
	// BCC SYSTEMS looking for lattices a a a 109.471 109.471 109.471
	if(VERBOSE_PROGRESS)  cerr << "BCC SYSTEMS looking for lattices a a a 109.471 109.471 109.471" << endl;
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                             // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if((choice==0 && (aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps)))) { // just to help not necessary
	    // BCC:  a1=[-a/2 a/2 a/2] a2=[a/2 -a/2 a/2]/ a3=[a/2 a/2 -a/2]     a a a 109.471 109.471 109.471 primitive
	    // BCC:  [a 0 0] [0 a 0] [0 0 a] with 0,0,0 and 1/2 1/2 1/2 basis
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps) &&
               SYM::checkAngle(b,c,alpha,109.471,is_deg,eps) && SYM::checkAngle(a,c,beta,109.471,is_deg,eps) && SYM::checkAngle(a,b,gamma,109.471,is_deg,eps)) { //DX
	       //aurostd::isequal(alpha,109.471,epsang) && aurostd::isequal(beta,109.471,epsang) && aurostd::isequal(gamma,109.471,epsang)) {
	      found=TRUE;
	      str_sp.bravais_lattice_type="BCC";
	      str_sp.bravais_lattice_variation_type="BCC";
	      str_sp.bravais_lattice_system="cubic";
	      str_sp.pearson_symbol="cI"+aurostd::utype2string(2*str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=2.0*a/sqrt(3.0); // make the standard primtive bcc
	      plattice[1][1]=-ac/2.0;plattice[1][2]= ac/2.0;plattice[1][3]= ac/2.0; // make the standard primtive bcc
	      plattice[2][1]= ac/2.0;plattice[2][2]=-ac/2.0;plattice[2][3]= ac/2.0; // make the standard primtive bcc
	      plattice[3][1]= ac/2.0;plattice[3][2]= ac/2.0;plattice[3][3]=-ac/2.0; // make the standard primtive bcc
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // make the STANDARD CONVENTIONAL
	      // C=supercell*P => supercell = C*inv(P)
	      str_sc=GetSuperCell(str_sp,0,1,1,1,0,1,1,1,0); // very cool (dual of bcc is fcc)
	      str_sc.bravais_lattice_type="CUB";
	      str_sc.bravais_lattice_variation_type="CUB";
	      str_sc.bravais_lattice_system="cubic";
	      str_sc.pearson_symbol="cP"+aurostd::utype2string(str_sc.atoms.size());
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE BCC
	// --------------------------------------------------------------------------
	// FCC FROM PRIMITIVE
	// FCC SYSTEMS looking for lattices a a a 60.0 60.0 60.0
	if(VERBOSE_PROGRESS)  cerr << "FCC SYSTEMS looking for lattices a a a 60.0 60.0 60.0" << endl;
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                           // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if((choice==0 && (aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps)))) { // just to help not necessary
	    rdata=Getabc_angles(rlattice,DEGREES);a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    // FCC:  a1=[0 a/2 a/2]  a2=[a/2 0 a/2]  a3=[a/2 a/2 0] a a a 60.0 60.0 60.0
	    // FCC:  [a 0 0] [0 a 0] [0 0 a] with 0,0,0 and 0 1/2 1/2 + 1/2 0 1/2 + 1/2 1/2 0 basis
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps) &&
               SYM::checkAngle(b,c,alpha,60.0,is_deg,eps) && SYM::checkAngle(a,c,beta,60.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,60.0,is_deg,eps)) { //DX
	       //DX aurostd::isequal(alpha,60.0,epsang) && aurostd::isequal(beta,60.0,epsang) && aurostd::isequal(gamma,60.0,epsang)) {
	      found=TRUE;
	      str_sp.bravais_lattice_type="FCC";
	      str_sp.bravais_lattice_variation_type="FCC";
	      str_sp.bravais_lattice_system="cubic";
	      str_sp.pearson_symbol="cF"+aurostd::utype2string(4*str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=2.0*a/sqrt(3.0); // make the standard primtive fcc
	      plattice[1][1]=   0.0;plattice[1][2]=ac/2.0;plattice[1][3]=ac/2.0; // make the standard primtive fcc
	      plattice[2][1]=ac/2.0;plattice[2][2]=   0.0;plattice[2][3]=ac/2.0; // make the standard primtive fcc
	      plattice[3][1]=ac/2.0;plattice[3][2]=ac/2.0;plattice[3][3]=   0.0; // make the standard primtive fcc
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // make the STANDARD CONVENTIONAL
	      // C=supercell*P => supercell = C*inv(P)
	      str_sc=GetSuperCell(str_sp,-1,1,1,1,-1,1,1,1,-1); // very cool (dual of fcc is bcc)
	      str_sc.bravais_lattice_type="CUB";
	      str_sc.bravais_lattice_variation_type="CUB";
	      str_sc.bravais_lattice_system="cubic";
	      str_sc.pearson_symbol="cP"+aurostd::utype2string(str_sc.atoms.size());
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE FCC
	// --------------------------------------------------------------------------
      } // DONE CUBIC CRYSTAL SYSTEMS
      // ***************************************************************************************************
  
      // ***************************************************************************************************
      // TETRAGONAL CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="tetragonal")) {
	// cerr << "TETRAGONAL CRYSTAL SYSTEM SEARCH " << vrlattice1.size() << " " << eps << endl;
	// TET SYSTEM
	if(VERBOSE_PROGRESS) cerr << "TET SYSTEMS looking for lattices a a c 90.0 90.0 90.0" << endl;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  //   cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << endl;
	  // DX if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
          if(SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	    // cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << endl;
	    if((choice==0 && ((aurostd::isequal(a,b,eps)&&aurostd::isdifferent(b,c,eps))||(aurostd::isequal(b,c,eps)&&aurostd::isdifferent(c,a,eps))||(aurostd::isequal(c,a,eps)&&aurostd::isdifferent(a,b,eps)))) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps)&&aurostd::isequal(b,c,eps))))) {
	      // if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
    
	      // TET: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 c]: a a c 90 90 90  primitive
	      // TET: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 c]: a a c 90 90 90  conventional
	      // ordering vectors to a a c
	      if(aurostd::isequal(a,b,eps)) {;}                                // a a c, ok
	      if(aurostd::isequal(a,c,eps)) {swap_rows(rlattice,2,3);for(uint ii=1;ii<=3;ii++) rlattice[1][ii]=-rlattice[1][ii];}  // a c a, swap 2,3
	      if(aurostd::isequal(b,c,eps)) {swap_rows(rlattice,1,3);for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=-rlattice[2][ii];}  // c a a, swap 1,3  
	      rdata=Getabc_angles(rlattice,DEGREES);
	      a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6]; // DONE
	      // TET a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 c]: a a c 90 90 90
	      found=TRUE;
	      str_sp.bravais_lattice_type="TET";
	      str_sp.bravais_lattice_variation_type="TET";
	      str_sp.bravais_lattice_system="tetragonal";
	      str_sp.pearson_symbol="tP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;cc=c; // make the standard primitive tet
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive tet
	      plattice[2][1]=0.0;plattice[2][2]= ac;plattice[2][3]=0.0; // make the standard primitive tet
	      plattice[3][1]=0.0;plattice[3][2]=0.0;plattice[3][3]= cc; // make the standard primitive tet
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // do the
	      // TET standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE TET
	// --------------------------------------------------------------------------
	// BCT FROM PRIMITIVE (TYPE 2, SMART)
	// looking for rlattice that is a transformed version of a1=[-a/2 a/2 c/2] a2=[a/2 -a/2 c/2] a3=[a/2 a/2 -c/2]
	if(VERBOSE_PROGRESS) cerr << "BCT FROM PRIMITIVE looking for rlattice that is a transformed version of a1=[-a/2 a/2 c/2] a2=[a/2 -a/2 c/2] a3=[a/2 a/2 -c/2]" << endl;
	xmatrix<double> bct(3,3);
	bct[1][1]=0.0;bct[1][2]=1.0;bct[1][3]=1.0;
	bct[2][1]=1.0;bct[2][2]=0.0;bct[2][3]=1.0;
	bct[3][1]=1.0;bct[3][2]=1.0;bct[3][3]=0.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(bct*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of BCT then bct*rlattice==conventional lattice of BCT
	  // STD_CONV  [a 0 0,  0 a 0 , 0 0 c ]
	  // STD_PRIM: a1=[-a/2 a/2 c/2] a2=[a/2 -a/2 c/2] a3=[a/2 a/2 -c/2]
	  //if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
          if(SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) {
	    if((choice==0 && ((aurostd::isequal(a,b,eps)&&aurostd::isdifferent(b,c,eps)))) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps)&&aurostd::isequal(b,c,eps))))) {
	      ac=a;cc=c;
	      //DX if(abs(ac*ac*cc-2.0*volume)<eps) {
	      if(abs(ac*ac*cc-2.0*volume)<eps_volume) { //DX
		found=TRUE;
		str_sp.bravais_lattice_type="BCT";
		str_sp.bravais_lattice_system="tetragonal";
		str_sp.pearson_symbol="tI"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]=-ac/2.0;plattice[1][2]= ac/2.0;plattice[1][3]= cc/2.0; // make the standard primitive bct
		plattice[2][1]= ac/2.0;plattice[2][2]=-ac/2.0;plattice[2][3]= cc/2.0; // make the standard primitive bct
		plattice[3][1]= ac/2.0;plattice[3][2]= ac/2.0;plattice[3][3]=-cc/2.0; // make the standard primitive bct
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
		if(cc<ac) {str_sp.bravais_lattice_variation_type="BCT1";}
		//DX - 9/7/17 - Need to include a=c if(ac<cc) {str_sp.bravais_lattice_variation_type="BCT2";}
		if(ac<=cc) {str_sp.bravais_lattice_variation_type="BCT2";} //DX 9/7/17 - Add a=c case; pick BCT2 since less high sym points 
		// do the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,bct);
		str_sc.bravais_lattice_type="TET";
		str_sc.bravais_lattice_variation_type="TET";
		str_sc.bravais_lattice_system="tetragonal";
		str_sc.pearson_symbol="tP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    } // DONE BCT
	  }
	}
      } // DONE TETRAGONAL CRYSTAL SYSTEMS
      // **************************************************************************************************
  
      // ***************************************************************************************************
      // TRIGONAL CRYSTAL SYSTEM                          
      if(found==FALSE && (crystal_system=="all" || crystal_system=="trigonal")) {  // trigonal can be rhl or hex
	if(VERBOSE_PROGRESS) cerr << "RHL SYSTEMS looking for lattices a a a alpha alpha alpha" << endl;
	//DX double rhl_ratio=1.0; // you might need to reduce this to reproduce tiny distortions... insted of increasing eps, you reduce it
	xmatrix<double> rhl(3,3);
	rhl[1][1]=1.0;rhl[1][2]=-1.0;rhl[1][3]= 0.0;
	rhl[2][1]=0.0;rhl[2][2]= 1.0;rhl[2][3]=-1.0;
	rhl[3][1]=1.0;rhl[3][2]= 1.0;rhl[3][3]= 1.0;
	// RHL SYSTEMS looking for lattices a a a alpha alpha alpha
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                           // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  //  cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << endl;
	  if(choice==0) { // only pristine
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps)) {  // aaa
	      //  cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << " " << epsang << endl;
	      //DX if(aurostd::isequal(alpha,beta,epsang) && aurostd::isequal(beta,gamma,epsang)) {  // alpha alpha alpha
              if(SYM::checkAngle(b,c,alpha,beta,is_deg,eps) && SYM::checkAngle(a,c,beta,gamma,is_deg,eps)) { //DX
		//  cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << " " << epsang << endl;
		//DX if(aurostd::isdifferent(alpha,90.0,epsang/rhl_ratio) && // no cubic
		//DX    aurostd::isdifferent(alpha,60.0,epsang/rhl_ratio) && aurostd::isdifferent(alpha,180.0-60.0,epsang/rhl_ratio) && // no fcc (60 and 120)
		//DX   aurostd::isdifferent(alpha,109.471,epsang/rhl_ratio) && aurostd::isdifferent(alpha,180.0-109.471,epsang/rhl_ratio)) {  // no bcc stuff
		//DX 9/6/17 [OBSOLETE] if(!SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && // no cubic
		//DX 9/6/17 [OBSOLETE]    !SYM::checkAngle(b,c,alpha,60.0,is_deg,eps) && !SYM::checkAngle(b,c,alpha,120.0,is_deg,eps) && // no fcc (60 and 120)
		//DX 9/6/17 [OBSOLETE]    !SYM::checkAngle(b,c,alpha,109.471,is_deg,eps) && !SYM::checkAngle(b,c,alpha,180.0-109.471,is_deg,eps)) {  // no bcc stuff
		  //   cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << " " << epsang << endl;
		  // make standard primitive
		  // RHL: a1=[a*cos(alpha/2) -a*sin(alpha/2) 0]
		  // RHL: a2=[a*cos(alpha/2)  a*sin(alpha/2) 0]
		  // RHL: a3=[a*cos(alpha)/cos(alpha/2) 0 a*sqrt(1-cos(alpha)*cos(alpha)/cos(alpha/2)/cos(alpha/2))]
		  found=TRUE;
		  str_sp.bravais_lattice_type="RHL";
		  str_sp.bravais_lattice_system="rhombohedral";
		  str_sp.pearson_symbol="hR"+aurostd::utype2string(str_sp.atoms.size()); // FIX
		  if(alpha<90.0) {
		    str_sp.bravais_lattice_variation_type="RHL1";
		  } else {
		    str_sp.bravais_lattice_variation_type="RHL2";
		  }
		  alpha=alpha*pi/180.0; // change to RADIANTS !
		  plattice[1][1]=a*cos(alpha/2.0);plattice[1][2]=-a*sin(alpha/2.0);plattice[1][3]=0.0; // make the standard primitive rhl
		  plattice[2][1]=a*cos(alpha/2.0);plattice[2][2]= a*sin(alpha/2.0);plattice[2][3]=0.0; // make the standard primitive rhl
		  plattice[3][1]=a*cos(alpha)/cos(alpha/2.0);plattice[3][2]=0.0;
		  plattice[3][3]=a*sqrt(1.0-cos(alpha)*cos(alpha)/cos(alpha/2.0)/cos(alpha/2.0));      // make the standard primitive rhl
                  //DX 9/29/17 - check if sqrt results in nan - START 
                  if(std::isnan(plattice[3][3])){
                    found=FALSE;
                    continue;
                  }
                  //DX 9/29/17 - check if sqrt results in nan - END
		  LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
		  // make the STANDARD CONVENTIONAL
#ifndef _RHL_HEX_SC_
		  // RHL standard_conventional=standard_primitive
		  str_sc=str_sp;
#else
		  // C=supercell*P => supercell = C*inv(P)
		  str_sc=GetSuperCell(str_sp,rhl); // very cool (make it hexagonal)
		  str_sc.bravais_lattice_type="HEX";
		  str_sc.bravais_lattice_variation_type="HEX";
		  str_sc.bravais_lattice_system="hexagonal";
		  str_sc.pearson_symbol="hP"+aurostd::utype2string(str_sc.atoms.size());
#endif
		  if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
		//DX 9/6/17 [OBSOLETE]}
	      }
	    }
	  }
	} // DONE RHL  
	if(LDEBUG) cerr << "TRIGONAL CRYSTAL SYSTEM found=" << found << endl;
      } // DONE TRIGONAL CRYSTAL SYSTEMS
      // ***************************************************************************************************

      // ***************************************************************************************************
      // HEXAGONAL CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="hexagonal" || crystal_system=="trigonal")) {
	if(VERBOSE_PROGRESS) cerr << "HEX SYSTEMS looking for lattices a a c 90 90 120" << endl;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // HEX a1=[a/2 -a/2*sqrt(3) 0]  a2=[a/2 a/2*sqrt(3) 0]  a3=[0 0 c]   a a c 90 90 120 primitive
	  // HEX a1=[a/2 -a/2*sqrt(3) 0]  a2=[a/2 a/2*sqrt(3) 0]  a3=[0 0 c]   a a c 90 90 120 conventional
	  if(choice==0) { // only pristine
	    //DX if((aurostd::isequal(a,b,eps) && aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang)) &&
	    //DX    (aurostd::isequal(gamma,120.0,epsang) || aurostd::isequal(gamma,60.0,epsang))) {
	    //DX   if(aurostd::isequal(gamma,60.0,epsang)) {
	    if((aurostd::isequal(a,b,eps) && SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps)) &&
	       (SYM::checkAngle(a,b,gamma,120.0,is_deg,eps) || SYM::checkAngle(a,b,gamma,60.0,is_deg,eps))) {
	      if(SYM::checkAngle(a,b,gamma,60.0,is_deg,eps)) {
		for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=rlattice[2][ii]-rlattice[1][ii];  //  a1 (a2-a1) a3 are 90 90 120
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	      }
	      found=TRUE;
	      str_sp.bravais_lattice_type="HEX";
	      str_sp.bravais_lattice_variation_type="HEX";
	      str_sp.bravais_lattice_system="hexagonal";
	      str_sp.pearson_symbol="hP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;cc=c; // make the conventional hex
	      plattice[1][1]=ac/2.0;plattice[1][2]=-ac/2.0*sqrt(3.0);plattice[1][3]=0.0; // make the conventional hex
	      plattice[2][1]=ac/2.0;plattice[2][2]= ac/2.0*sqrt(3.0);plattice[2][3]=0.0; // make the conventional hex
	      plattice[3][1]=   0.0;plattice[3][2]=              0.0;plattice[3][3]= cc; // make the conventional hex
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // make standard conventional lattice
	      // HEX standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE HEX
	if(LDEBUG) cerr << "HEXAGONAL CRYSTAL SYSTEM found=" << found << endl;
      } // DONE HEXAGONAL CRYSTAL SYSTEMS
      // ***************************************************************************************************
    
      // ***************************************************************************************************
      // ORTHORHOMBIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="orthorhombic")) {
	if(VERBOSE_PROGRESS) cerr << "ORC SYSTEM looking for lattices a<=b<=c 90 90 90" << endl;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  //DX if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	  if(SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	    if((choice==0 && ((aurostd::isdifferent(a,b,eps) && aurostd::isdifferent(b,c,eps) && aurostd::isdifferent(c,a,eps)))) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps) || aurostd::isequal(b,c,eps) || aurostd::isequal(c,a,eps))))) {
	      // ORC: a1=[a 0 0]  a2=[0 b 0]  a3=[0 0 c]: a<b<c 90 90 90  primitive
	      // ORC: a1=[a 0 0]  a2=[0 b 0]  a3=[0 0 c]: a<b<c 90 90 90  conventional
	      // ordering vectors to get a<b<c
	      rdata=Getabc_angles(rlattice,DEGREES);
	      a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      if(aurostd::isequal(a,minabc,eps)) {;}                           // min * * , ok
	      if(aurostd::isequal(b,minabc,eps)) {                             // * min * , swap 1,2
		swap_rows(rlattice,1,2);if(det(rlattice)<0.0) {for(uint ii=1;ii<=3;ii++) rlattice[3][ii]=-rlattice[3][ii];}
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      }
	      if(aurostd::isequal(c,minabc,eps) || aurostd::isequal(a,maxabc,eps)) { // * * min , swap 1,3    OR    max * * , swap 1,3
		swap_rows(rlattice,1,3);if(det(rlattice)<0.0) {for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=-rlattice[2][ii];}
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      }
	      if(aurostd::isequal(b,maxabc,eps)) {                             // * max * , swap 2,3
		swap_rows(rlattice,3,2);if(det(rlattice)<0.0) {for(uint ii=1;ii<=3;ii++) rlattice[1][ii]=-rlattice[1][ii];}
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      }
	      if(aurostd::isequal(c,maxabc,eps)) {;}                           // * * max , ok
	      // done
	      a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6]; // DONE
	      // ORC a1=[a 0 0]  a2=[0 b 0]  a3=[0 0 c]: a<b<c 90 90 90
	      found=TRUE;
	      str_sp.bravais_lattice_type="ORC";
	      str_sp.bravais_lattice_variation_type="ORC";
	      str_sp.bravais_lattice_system="orthorhombic";
	      str_sp.pearson_symbol="oP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;bc=b;cc=c; // make the standard primitive orc
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive orc
	      plattice[2][1]=0.0;plattice[2][2]= bc;plattice[2][3]=0.0; // make the standard primitive orc
	      plattice[3][1]=0.0;plattice[3][2]=0.0;plattice[3][3]= cc; // make the standard primitive orc
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // ORC standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE ORC
	// --------------------------------------------------------------------------
	// ORCF from PRIMITIVE
	if(VERBOSE_PROGRESS) cerr << "ORCF SYSTEM looking for orcf*rlattice==conventional lattice of ORCF" << endl;
	xmatrix<double> orcf(3,3);
	orcf[1][1]=-1.0;orcf[1][2]=1.0;orcf[1][3]=1.0;
	orcf[2][1]=1.0;orcf[2][2]=-1.0;orcf[2][3]=1.0;
	orcf[3][1]=1.0;orcf[3][2]=1.0;orcf[3][3]=-1.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(orcf*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of ORCF then orcf*rlattice==conventional lattice of ORCF !
	  // standard conventional:
	  // [ac 0 0,  0 bc 0 , 0 0 cc ]
	  // standard primitive:
	  // a1=[0 bc/2 cc/2]
	  // a2=[ac/2 0 cc/2]
	  // a3=[ac/2 bc/2 0]
	  //DX if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
          if(SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	    if((choice==0 && (a<b && b<c)) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps) && b<c) || (a<b && aurostd::isequal(b,c,eps)))) ||
	       (choice==2 && ((aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps))))) {
	      ac=a;bc=b;cc=c;
	      //DX if(abs(ac*bc*cc-4.0*volume)<eps) {
	      if(abs(ac*bc*cc-4.0*volume)<eps_volume) { //DX
		found=TRUE;
		str_sp.bravais_lattice_type="ORCF";
		str_sp.bravais_lattice_system="orthorhombic";
		str_sp.pearson_symbol="oF"+aurostd::utype2string(4*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		double mismatch=pow(1.0/ac,2.0)/(pow(1.0/bc,2.0)+pow(1.0/cc,2.0));
		//  cerr << mismatch << " " << eps/5.0 << endl;  // DEBUG purposes
		if(mismatch>1+eps/5.0) {str_sp.bravais_lattice_variation_type="ORCF1";}
		if(mismatch<1-eps/5.0) {str_sp.bravais_lattice_variation_type="ORCF2";}
		if(aurostd::isequal(mismatch,1.0,eps/5.0)) {str_sp.bravais_lattice_variation_type="ORCF3";}
		plattice[1][1]=   0.0;plattice[1][2]=bc/2.0;plattice[1][3]=cc/2.0; // make the standard primitive orcf
		plattice[2][1]=ac/2.0;plattice[2][2]=   0.0;plattice[2][3]=cc/2.0; // make the standard primitive orcf
		plattice[3][1]=ac/2.0;plattice[3][2]=bc/2.0;plattice[3][3]=   0.0; // make the standard primitive orcf
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // plattice is a rotated version of rlattice
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,orcf); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="ORC";
		str_sc.bravais_lattice_variation_type="ORC";
		str_sc.bravais_lattice_system="orthorhombic";
		str_sc.pearson_symbol="oP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  }
	}
	// --------------------------------------------------------------------------
	// ORCI FROM PRIMITIVE SMART
	if(VERBOSE_PROGRESS) cerr << "ORCF SYSTEM looking for orci*rlattice==conventional lattice of ORCI" << endl;
	xmatrix<double> orci(3,3);
	orci[1][1]=0.0;orci[1][2]=1.0;orci[1][3]=1.0;
	orci[2][1]=1.0;orci[2][2]=0.0;orci[2][3]=1.0;
	orci[3][1]=1.0;orci[3][2]=1.0;orci[3][3]=0.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(orci*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of ORCI then orci*rlattice==conventional lattice of ORCI !
	  // ORCI
	  // Standard Conventional:
	  // [ac 0 0,  0 bc 0 , 0 0 cc ]
	  // Standard Primitive:
	  // a1=[-a/2 b/2 c/2]
	  // a2=[a/2 -b/2 c/2]
	  // a3=[a/2 b/2 -c/2]
	  //DX if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
          if(SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	    if((choice==0 && (a<b && b<c)) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps) && b<c) || (a<b && aurostd::isequal(b,c,eps)))) ||
	       (choice==2 && ((aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps))))) {
	      ac=a;bc=b;cc=c;
	      //DX if(abs(ac*bc*cc-2.0*volume)<eps) {
	      if(abs(ac*bc*cc-2.0*volume)<eps_volume) { //DX
		found=TRUE;
		str_sp.bravais_lattice_type="ORCI";
		str_sp.bravais_lattice_variation_type="ORCI";
		str_sp.bravais_lattice_system="orthorhombic";
		str_sp.pearson_symbol="oI"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]=-ac/2.0;plattice[1][2]= bc/2.0;plattice[1][3]= cc/2.0; // make the standard primitive orci
		plattice[2][1]= ac/2.0;plattice[2][2]=-bc/2.0;plattice[2][3]= cc/2.0; // make the standard primitive orci
		plattice[3][1]= ac/2.0;plattice[3][2]= bc/2.0;plattice[3][3]=-cc/2.0; // make the standard primitive orci
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,orci); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="ORC";
		str_sc.bravais_lattice_variation_type="ORC";
		str_sc.bravais_lattice_system="orthorhombic";
		str_sc.pearson_symbol="oP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  } // DONE ORCI
	}
	// --------------------------------------------------------------------------
	// ORCC FROM PRIMITIVE SMART
	if(VERBOSE_PROGRESS) cerr << "ORCC SYSTEM looking for orcc*rlattice==conventional lattice of ORCC" << endl;
	xmatrix<double> orcc(3,3);
	orcc[1][1]= 1.0;orcc[1][2]=1.0;orcc[1][3]=0.0;
	orcc[2][1]=-1.0;orcc[2][2]=1.0;orcc[2][3]=0.0;
	orcc[3][1]= 0.0;orcc[3][2]=0.0;orcc[3][3]=1.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(orcc*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of ORCC then orcc*rlattice==conventional lattice of ORCC !
	  // ORCC (a<b)
	  // [ac 0 0,  0 bc 0 , 0 0 cc ]      Conventional direct lattice:
	  // Direct Lattice (VASP and Bilbao):
	  // a1=[a/2 -b/2 0]
	  // a2=[a/2 b/2 0]
	  // a3=[0 0 c]
	  // pick only the case when a<b
	  //DX if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
          if(SYM::checkAngle(b,c,alpha,90.0,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	    if((choice==0 && (a<b)) || (choice==1 && (aurostd::isequal(a,b,eps)))) {
	      ac=a;bc=b;cc=c;
	      //DX if(abs(ac*bc*cc-2.0*volume)<eps) {
	      if(abs(ac*bc*cc-2.0*volume)<eps_volume) { //DX
		found=TRUE;
		str_sp.bravais_lattice_type="ORCC";
		str_sp.bravais_lattice_variation_type="ORCC";
		str_sp.bravais_lattice_system="orthorhombic";
		str_sp.pearson_symbol="oS"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]=ac/2.0;plattice[1][2]=-bc/2.0;plattice[1][3]=0.0; // make the standard primitive orcc
		plattice[2][1]=ac/2.0;plattice[2][2]= bc/2.0;plattice[2][3]=0.0; // make the standard primitive orcc
		plattice[3][1]=   0.0;plattice[3][2]=    0.0;plattice[3][3]= cc; // make the standard primitive orcc
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,orcc); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="ORC";
		str_sc.bravais_lattice_variation_type="ORC";
		str_sc.bravais_lattice_system="orthorhombic";
		str_sc.pearson_symbol="oP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  } // DONE ORCC
	}
      }
      // ***************************************************************************************************

      // ***************************************************************************************************
      // MONOCLINIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="monoclinic")) {
	// MCL SYSTEM
	// STANDARD: a b<=c alpha<90 90 90  (a b c can be whatever so wikipedia must be relaxed... but put alpha as close as possible to 90)
	// with abc ordering
	if(VERBOSE_PROGRESS) cerr << "MCL SYSTEM a b<=c alpha<~90 90 90" << endl;
	double MCL_alpha_max=0.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  //DX if(alpha<90.0 && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	  if(alpha<90.0 && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	    if((choice==0 && b<c)|| (choice==1 && aurostd::isequal(b,c,eps))) {  // && not necessary aurostd::isdifferent(a,b,eps) && aurostd::isdifferent(b,c,eps) && aurostd::isdifferent(c,a,eps)
	      if(alpha>MCL_alpha_max) MCL_alpha_max=alpha;
	    }
	  }
	}
	// got alpha_max
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  //DX if(aurostd::isequal(alpha,MCL_alpha_max,epsang)  && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
          if(SYM::checkAngle(b,c,alpha,MCL_alpha_max,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	    if((choice==0 && b<c)|| (choice==1 && aurostd::isequal(b,c,eps))) {  // && not necessary aurostd::isdifferent(a,b,eps) && aurostd::isdifferent(b,c,eps) && aurostd::isdifferent(c,a,eps)
	      // MCL: a1=[a 0 0]  a2=[0 b 0]  a3=[0 c*cos(alpha) c*sin(alpha)]: a b<c alpha<90 90 90  primitive
	      // MCL: a1=[a 0 0]  a2=[0 b 0]  a3=[0 c*cos(alpha) c*sin(alpha)]: a b<c alpha<90 90 90  conventional
	      found=TRUE;
	      str_sp.bravais_lattice_type="MCL";
	      str_sp.bravais_lattice_variation_type="MCL";
	      str_sp.bravais_lattice_system="monoclinic";
	      str_sp.pearson_symbol="mP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;bc=b;cc=c; // make the standard primitive mcl
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive mcl
	      plattice[2][1]=0.0;plattice[2][2]= bc;plattice[2][3]=0.0; // make the standard primitive mcl
	      plattice[3][1]=0.0;plattice[3][2]=cc*cos(alpha/180.0*pi);plattice[3][3]=cc*sin(alpha/180.0*pi); // make the standard primitive mcl
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // do the
	      // MCL standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE MCL
	// --------------------------------------------------------------------------
	// MCLC from PRIMITIVE
	// MCLC (a b c can be whatever so wikipedia must be relaxed... but put alpha as close as possible to 90)
	if(VERBOSE_PROGRESS) cerr << "MCLC SYSTEM a b c alpha<~90 90 90" << endl;
	xmatrix<double> mclc(3,3);
	mclc[1][1]=1.0;mclc[1][2]=-1.0;mclc[1][3]=0.0;
	mclc[2][1]=1.0;mclc[2][2]=1.0;mclc[2][3]=0.0;
	mclc[3][1]=0.0;mclc[3][2]=0.0;mclc[3][3]=1.0;
	if(choice==0) { // only pristine
	  double MCLC_alpha_max=0.0,testphi,alpharad;
	  for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	    rlattice=vrlattice1.at(i);
	    rdata=Getabc_angles(mclc*rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    //DX if(alpha<90.0 && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
            if(alpha<90.0 && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	      ac=a;bc=b;cc=c;
	      //DX if(abs(det(mclc*rlattice)-2.0*volume)<eps) {
	      if(abs(det(mclc*rlattice)-2.0*volume)<eps_volume) { //DX
		if(alpha>MCLC_alpha_max) MCLC_alpha_max=alpha;
	      }
	    }
	  }
	  // got alpha_max
	  for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	    rlattice=vrlattice1.at(i);
	    rdata=Getabc_angles(mclc*rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    // if rlattice==primitive lattice of MCLC then mclc*rlattice==conventional lattice of MCLC !
	    // MCLC a b c alpha<90.0 90 90
	    // [ac 0 0] [0 bc 0] [0 cc*cos(alpha) cc*sin(alpha) ]    Conventional direct lattice:
	    // Direct Lattice (VASP and Bilbao):
	    // a1=[ a/2 b/2 0] a2=[-a/2 b/2 0] a3=[0 c*cos(alpha) c*sin(alpha)]
	    //DX if(aurostd::isequal(alpha,MCLC_alpha_max,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
            if(SYM::checkAngle(b,c,alpha,MCLC_alpha_max,is_deg,eps) && SYM::checkAngle(a,c,beta,90.0,is_deg,eps) && SYM::checkAngle(a,b,gamma,90.0,is_deg,eps)) { //DX
	      ac=a;bc=b;cc=c;
	      //DX if(abs(det(mclc*rlattice)-2.0*volume)<eps) {
	      if(abs(det(mclc*rlattice)-2.0*volume)<eps_volume) { //DX
		str_sp.bravais_lattice_type="MCLC";
		str_sp.bravais_lattice_system="monoclinic";
		str_sp.pearson_symbol="mS"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]= ac/2.0;plattice[1][2]=bc/2.0;plattice[1][3]=0.0; // make the standard primitive mclc
		plattice[2][1]=-ac/2.0;plattice[2][2]=bc/2.0;plattice[2][3]=0.0; // make the standard primitive mclc
		plattice[3][1]=  0.0;
		plattice[3][2]=cc*cos(alpha*pi/180.0);
		plattice[3][3]=cc*sin(alpha*pi/180.0);                           // make the standard primitive mclc
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);                    // project out of rlattice and in plattice
		str_sp.FixLattices();
		klattice=str_sp.klattice;
		kdata=Getabc_angles(klattice,DEGREES);kalpha=kdata[4];kbeta=kdata[5];kgamma=kdata[6];
		//DX if(aurostd::isequal(kgamma,90.0,epsang)) {
		if(SYM::checkAngle(kdata[1],kdata[2],kgamma,90.0,is_deg,eps)) { //DX
		  found=TRUE; str_sp.bravais_lattice_variation_type="MCLC2";
		} else {
		  if(kgamma>90.0) {
		    found=TRUE; str_sp.bravais_lattice_variation_type="MCLC1";
		  } else {
		    // kgamma<90
		    alpharad=alpha*PI/180.0;
		    testphi=0.75-0.25*b*b/a/a-0.25*b/tan(alpharad)*(1/c/sin(alpharad)-1/b/tan(alpharad));
		    if(aurostd::isequal(testphi,0.5,1e-2)) { //DX: This is a normalized quantity, so we do not need to use checkAngle, should be approximately 0.5
		      found=TRUE; str_sp.bravais_lattice_variation_type="MCLC4";
		    } else {
		      if(testphi>0.5) {
			found=TRUE; str_sp.bravais_lattice_variation_type="MCLC3";
		      } else {
			found=TRUE; str_sp.bravais_lattice_variation_type="MCLC5";
		      }
		    }
		  }
		}
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,mclc); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="MCL";
		str_sc.bravais_lattice_variation_type="MCL";
		str_sc.bravais_lattice_system="monoclinic";
		str_sc.pearson_symbol="mP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  }
	}
      } // DONE MCLC
      // ***************************************************************************
      // TRICLINIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="triclinic"))  {
	found=FALSE;
	double weight_TRI1a=1e6,weight_TRI2a=1e6,weight_TRI1b=1e6,weight_TRI2b=1e6;
	xmatrix<double> rrlattice_TRI1a(3,3),rrlattice_TRI2a(3,3),rrlattice_TRI1b(3,3),rrlattice_TRI2b(3,3);
	if(choice==0) { // only pristine
	  for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++)  {
	    rlattice=vrlattice1.at(i);
	    rdata=Getabc_angles(rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    klattice=ReciprocalLattice(rlattice);
	    kdata=Getabc_angles(klattice,DEGREES);
	    ka=kdata[1];kb=kdata[2];kc=kdata[3];kalpha=kdata[4];kbeta=kdata[5];kgamma=kdata[6];
	    // find which better fits
	    if(kalpha>90.0 && kbeta>90.0 && kgamma>90.0) {
	      if((kalpha-90.0)+(kbeta-90.0)+(kgamma-90.0)<weight_TRI1a) {
		weight_TRI1a=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI1a=rlattice; } }
	    //DX if(kalpha>90.0 && kbeta>90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
	    if(kalpha>90.0 && kbeta>90.0 && SYM::checkAngle(ka,kb,kgamma,90.0,is_deg,eps)) { //DX
	      if((kalpha-90.0)+(kbeta-90.0)<weight_TRI2a) {
		weight_TRI2a=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI2a=rlattice; } }
	    if(kalpha<90.0 && kbeta<90.0 && kgamma<90.0) {
	      if((90.0-kalpha)+(90.0-kbeta)+(90.0-kgamma)<weight_TRI1b) {
		weight_TRI1b=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI1b=rlattice; } }
	    //DX if(kalpha<90.0 && kbeta<90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
	    if(kalpha<90.0 && kbeta<90.0 && SYM::checkAngle(ka,kb,kgamma,90.0,is_deg,eps)) {
	      if((90.0-kalpha)+(90.0-kbeta)<weight_TRI2b) {
		weight_TRI2b=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI2b=rlattice; } }
	  }
	  // now it found it
	  {
	    // lets pick the nicest TRI 1/2 a/b
	    a=b=c=alpha=beta=gamma=ka=kb=kc=kalpha=kbeta=kgamma=0.0;
	    bool found_TRI=FALSE;
	    if(weight_TRI1a<=min(weight_TRI2a,weight_TRI1b,weight_TRI2b) && found_TRI==FALSE) { rlattice=rrlattice_TRI1a; found_TRI=TRUE; }
	    if(weight_TRI2a<=min(weight_TRI1b,weight_TRI2b,weight_TRI1a) && found_TRI==FALSE) { rlattice=rrlattice_TRI2a; found_TRI=TRUE; }
	    if(weight_TRI1b<=min(weight_TRI2b,weight_TRI1a,weight_TRI2a) && found_TRI==FALSE) { rlattice=rrlattice_TRI1b; found_TRI=TRUE; }
	    if(weight_TRI2b<=min(weight_TRI1a,weight_TRI2a,weight_TRI1b) && found_TRI==FALSE) { rlattice=rrlattice_TRI2b; found_TRI=TRUE; }
	    rdata=Getabc_angles(rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //2-17-16 COses, Karol Jarolimek discovered B1Cd1Li1O3_ICSD_200615 was calculated to be TRI1B, not TRI1A like shown in
            //Comp. Mat. Sci. 49, 299 (2010) Fig 47.  Turns out we forgot to set ka...kgamma, so the following lines were added.
            //Interestingly, this system has weight_TRIa==weight_TRIb, so the check above (TRI1a before TRI1b) is very important.
            klattice=ReciprocalLattice(rlattice);
            kdata=Getabc_angles(klattice,DEGREES);
            ka=kdata[1];kb=kdata[2];kc=kdata[3];kalpha=kdata[4];kbeta=kdata[5];kgamma=kdata[6];	    
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // picked, put the label
	    if(kalpha>90.0 && kbeta>90.0 && kgamma>90.0) {
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI1a";
	    }
	    //DX if(kalpha>90.0 && kbeta>90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
            if(kalpha>90.0 && kbeta>90.0 && SYM::checkAngle(ka,kb,kgamma,90.0,is_deg,eps)) { //DX
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI2a";
	    }
	    if(kalpha<90.0 && kbeta<90.0 && kgamma<90.0) {
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI1b";
	    }
	    //DX if(kalpha<90.0 && kbeta<90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
            if(kalpha<90.0 && kbeta<90.0 && SYM::checkAngle(ka,kb,kgamma,90.0,is_deg,eps)) { //DX
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI2b";
	    }
	    // found=false;
	    // cerr << kalpha << " " << kbeta << " " << kgamma << " " << found << endl;
	    if(found) {
	      str_sp.bravais_lattice_system="triclinic";
	      str_sp.pearson_symbol="aP"+aurostd::utype2string(str_sp.atoms.size());
	      // make standard primitive
	      double alf=alpha*PI/180.0, bet=beta*PI/180.0, gam=gamma*PI/180.0;
	      plattice[1][1]=a; plattice[1][2]=0.0; plattice[1][3]=0.0;
	      plattice[2][1]=b*cos(gam); plattice[2][2]=b*sin(gam); plattice[2][3]=0.0;
	      plattice[3][1]=c*cos(bet);
	      plattice[3][2]=c*(cos(alf)-cos(bet)*cos(gam))/sin(gam);
	      plattice[3][3]=1.0-cos(alf)*cos(alf)-cos(bet)*cos(bet)-cos(gam)*cos(gam)+2*cos(alf)*cos(bet)*cos(gam);
	      plattice[3][3]=c*sqrt(plattice[3][3])/sin(gam);
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
	      str_sp.FixLattices();
	      // TRI standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  } //for vrlattice1
	} // choice
      }// TRIclinic DONE  
      // ***************************************************************************************************
      if(VERBOSE_PROGRESS) if(found) cerr << "str_sp.bravais_lattice_type=" << str_sp.bravais_lattice_type << endl;
      if(VERBOSE_PROGRESS) if(found) cerr << "str_sp.bravais_lattice_system=" << str_sp.bravais_lattice_system << endl;
      if(VERBOSE_PROGRESS) if(found) cerr << "str_sp.pearson_symbol=" << str_sp.pearson_symbol << endl;
    }
    // cerr << str_sp.bravais_lattice_type << endl; exit(0);
  
    if(found==FALSE) str_sp.bravais_lattice_type="UNKNOWN";
    if(found==TRUE) {
      str_sp.title=str_sp.title+"  ["+str_sp.bravais_lattice_type+","+str_sp.bravais_lattice_variation_type+","+str_sp.pearson_symbol+"]"+" (STD_PRIM doi:10.1016/j.commatsci.2010.05.010)";
    }
    if(found==FALSE) str_sc.bravais_lattice_type="UNKNOWN";
    if(found==TRUE) {
      str_sc.title=str_sc.title+"  ["+str_sc.bravais_lattice_type+","+str_sc.bravais_lattice_variation_type+","+str_sc.pearson_symbol+"]"+" (STD_CONV doi:10.1016/j.commatsci.2010.05.010)";
    }
    // DONE

    if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: X1 found=" << found << endl;
    
    if(found==FALSE) {
      str_sc.Standard_Lattice_calculated=TRUE;str_sc.Standard_Lattice_avoid=FALSE;
      str_sc.Standard_Lattice_primitive=FALSE;str_sc.Standard_Lattice_conventional=FALSE;
      str_sc.Standard_Lattice_has_failed=TRUE;
      str_sp.Standard_Lattice_calculated=TRUE;str_sp.Standard_Lattice_avoid=FALSE;
      str_sp.Standard_Lattice_primitive=FALSE;str_sp.Standard_Lattice_conventional=FALSE;
      str_sp.Standard_Lattice_has_failed=TRUE;
      str_sp=str_in;
      str_sp.GetPrimitive();
      str_sp.MinkowskiBasisReduction(); // shorten the vectors as much as possible and as perpendicular as possible
      str_sc=str_in; // copy it
      str_sc.MinkowskiBasisReduction(); // shorten the vectors as much as possible and as perpendicular as possible
    
      //    if(str_in.title!="NO_RECURSION")
/*
      if(0)      {
	cerr << "REMAKING LATTICE" << endl;
	xstructure str_in2(str_in),str_sp2(str_sp),str_sc2(str_sc);
	while(str_sp2.atoms.size()>1) str_sp2.RemoveAtom(0);
	while(str_sc2.atoms.size()>1) str_sc2.RemoveAtom(0);
	str_in2.title="NO_RECURSION";
	LATTICE::Standard_Lattice_Structure(str_in2,str_sp2,str_sc2,eps,epsang,time,symeps,histogram);
	cerr << str_sc2.bravais_lattice_type << endl;
	cerr << str_sc2.bravais_lattice_system << endl;
	cerr << str_sc2.pearson_symbol << endl;
	exit(0);
      }
*/
    }
    if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: X2 found=" << found << endl;

    if(found==TRUE) {
      str_sc.Standard_Lattice_calculated=TRUE;str_sc.Standard_Lattice_avoid=FALSE;
      str_sc.Standard_Lattice_primitive=FALSE;str_sc.Standard_Lattice_conventional=TRUE;
      str_sc.Standard_Lattice_has_failed=FALSE;
      str_sp.Standard_Lattice_calculated=TRUE;str_sp.Standard_Lattice_avoid=FALSE;
      str_sp.Standard_Lattice_primitive=TRUE;str_sp.Standard_Lattice_conventional=FALSE;
      str_sp.Standard_Lattice_has_failed=FALSE;
      str_sp.FixLattices();str_sc.FixLattices();
      // get reciprocal lattice
      // cerr << "STEP" << endl;
/*
      if(0 && str_in.title!="NO_RECURSION") {
	// RECIPROCAL
	xstructure str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc;
	str_reciprocal_in.lattice=str_sp.klattice;str_reciprocal_in.FixLattices();
	str_reciprocal_in.title="NO_RECURSION";
	_atom atom;str_reciprocal_in.AddAtom(atom);
	// LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang); //STEFANO OLD VERSION
	int ss=0; //JUNKAI
	LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang,ss,_EPS_);//JUNKAI
	str_sp.reciprocal_lattice_type=str_reciprocal_sp.bravais_lattice_type;
	str_sp.reciprocal_lattice_variation_type=str_reciprocal_sp.bravais_lattice_variation_type;//wahyu mod
	//str_sp.reciprocal_conventional_lattice_type=str_reciprocal_sp.bravais_lattice_system;
	str_sc.reciprocal_lattice_type=str_reciprocal_sp.bravais_lattice_type;
	str_sc.reciprocal_lattice_variation_type=str_reciprocal_sc.bravais_lattice_variation_type;//wahyu mod
	//str_sc.reciprocal_conventional_lattice_type=str_reciprocal_sp.bravais_lattice_system;
	// SUPERLATTICE
	xstructure str_superlattice_in,str_superlattice_sp,str_superlattice_sc;
	str_superlattice_in=str_sp;
	str_superlattice_in.title="NO_RECURSION";
	str_superlattice_in.IdenticalAtoms();
	str_superlattice_in.Minkowski_calculated=FALSE;
	//LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,eps,epsang); //STEFANO OLD VERSION
	ss=0; //JUNKAI
	LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,eps,epsang,ss,_EPS_);//JUNKAI
	str_sp.bravais_superlattice_type=str_superlattice_sp.bravais_lattice_type;
	str_sp.bravais_superlattice_variation_type=str_superlattice_sp.bravais_lattice_variation_type;
	str_sp.bravais_superlattice_system=str_superlattice_sp.bravais_lattice_system;
	str_sp.pearson_symbol_superlattice=str_superlattice_sp.pearson_symbol;
	str_sc.bravais_superlattice_type=str_superlattice_sp.bravais_lattice_type;
	str_sc.bravais_superlattice_variation_type=str_superlattice_sp.bravais_lattice_variation_type;
	str_sc.bravais_superlattice_system=str_superlattice_sp.bravais_lattice_system;
	str_sc.pearson_symbol_superlattice=str_superlattice_sp.pearson_symbol;
      }
*/
      str_sp.SetVolume(str_sp_volume);
      str_sc.SetVolume(str_sp_volume*((double) str_sc.atoms.size()/str_sp.atoms.size()));
      str_sp.neg_scale=str_sp_neg_scale;  // reload from backup
      str_sc.neg_scale=str_sp_neg_scale;  // reload from backup
    }

    if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: X3 found=" << found << endl;

    // last checks
    /* STEFANO OLD VERSION
       if(abs(NearestNeighbour(str_in)-NearestNeighbour(str_sp))>0.01) {
       cerr << "ERROR:  LATTICE::Standard_Lattice_Structure" << endl;
       cerr << "ERROR:  str_in.title=" << str_in.title << endl;
       cerr << "ERROR:  nn(str_in)=" << NearestNeighbour(str_in) << endl;
       cerr << "ERROR:  nn(str_sp)=" << NearestNeighbour(str_sp) << endl;
       cerr << "ERROR:  nn(str_sc)=" << NearestNeighbour(str_sc) << endl;
       // exit(0);
       }
    */  
    // some coding for test
    //   if(0) {
    //     _aflags aflags;APENNSY_Parameters params;params.LoadStructuresHQTC(aflags);
    //     for(uint i=1;i<params.structures[2].size();i++) {
    //       cout << "echo " << aurostd::RemoveCharacter(params.structures[2].at(i).name,'/') << endl;
    //       cout << "./aflow --proto " << aurostd::RemoveCharacter(params.structures[2].at(i).name,'/') << " Fe Co | ./aflow --sprim >> /dev/null " << endl;
    //     }
    //   }

    //  cerr << "OK=" << str_sp.bravais_lattice_type << endl;
    return found;
  }
}



namespace LATTICE {
  // transformation to standardization
  bool Standard_Lattice_Structure_20170101(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps,bool histogram) {//**************************JUNKAI EDITED*****************
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: BEGIN" << endl;
    //bool LATTICE::Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang) {   //STEFANO OLD VERSION
    // cerr << "eps=" << eps << " " << "epsang=" << epsang << endl;
    //  LDEBUG=TRUE;
    // starting from a str_in (whatever lattices) this routine returns a standard primitive (str_sp)
    // and a standard conventional (str_sc) following the ideas of stefano and wahyu
    if(str_in.Standard_Lattice_avoid==TRUE) return FALSE; // if you do not want to calculate
    // preparation
    str_sp=str_in; // copy it
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [1]" << endl;
    str_sp.GetPrimitive(0.005);
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [2]" << endl;
    str_sp.FixLattices();

    double str_sp_volume=str_sp.Volume();    // backups
    bool str_sp_neg_scale=str_sc.neg_scale;  // backups

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [3]" << endl;
    
    str_sp.ReScale(1.0);     // get it off the way, I might save and plug it back but it is not necessary now
    str_sp.neg_scale=FALSE;  // get it off the way, I might save and plug it back but it is not necessary now
    //   if(aurostd::isequal(str_sp.a,str_sp.b,eps) && aurostd::isequal(str_sp.b,str_sp.c,eps) && aurostd::isequal(str_sp.c,str_sp.a,eps)) {
    //     cerr << "LATTICE::Standard_Lattice_Structure: skipping str_sp.MinkowskiBasisReduction()" << endl;
    //   } else
    {str_sp.MinkowskiBasisReduction();} // shorten the vectors as much as possible and as perpendicular as possible

    str_sc=str_in; // copy it
    str_sc.FixLattices();
    str_sc.ReScale(1.0);     // get it off the way, I might save and plug it back but it is not necessary now
    str_sc.neg_scale=FALSE;  // get it off the way, I might save and plug it back but it is not necessary now
    //   if(aurostd::isequal(str_sc.a,str_sc.b,eps) && aurostd::isequal(str_sc.b,str_sc.c,eps) && aurostd::isequal(str_sc.c,str_sc.a,eps)) {
    //     cerr << "LATTICE::Standard_Lattice_Structure: skipping str_sc.MinkowskiBasisReduction()" << endl;
    //   } else
    {str_sc.MinkowskiBasisReduction();} // shorten the vectors as much as possible and as perpendicular as possible
    str_sc.ClearSymmetry();
  
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4]" << endl;

    xvector<int> dims(3);    // to scan
    dims=LatticeDimensionSphere(str_sp.lattice,RadiusSphereLattice(str_sp.lattice)/1.0);
    //  dims(1)=4;dims(2)=4;dims(3)=4;
    xmatrix<double> rlattice(3,3),plattice(3,3),clattice(3,3),klattice(3,3),t(3,3);
    xvector<double> aus(3),a1(3),a2(3),a3(3),rdata(6),kdata(6);
    double volume,vaus,nndist,cutoff,a,b,c,ka,kb,kc,alpha,beta,gamma,ac,bc,cc,minabc,maxabc,kalpha,kbeta,kgamma;
    // generate all the possible transformations
    bool found=FALSE;
    vector<xmatrix<double> > vrlattice1,vrlattice1_aaa;
    vector<xvector<double> > vvectors;
    volume=det(str_sp.lattice);
    vrlattice1.push_back(str_sp.lattice);

    nndist=1.0e6;
    cutoff=1.0e6;
    // find NN distance
    for(int i=-dims[1];i<=dims[1];i++)
      for(int j=-dims[2];j<=dims[2];j++)
	for(int k=-dims[3];k<=dims[3];k++) {
	  aus=i*str_sp.lattice(1)+j*str_sp.lattice(2)+k*str_sp.lattice(3);
	  if(modulus(aus)>eps && modulus(aus)<=nndist) nndist=modulus(aus);
	}

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4a]" << endl;

    string crystal_system="none";
    bool lattice_found=FALSE;
    //  str_sp.CalculateSymmetryPointGroup(FALSE);
    bool SYS_VERBOSE=FALSE;
    /* STEFANO OLD VERSION
       if(str_sp.pgroup_calculated==FALSE) str_sp.CalculateSymmetryPointGroup(SYS_VERBOSE); // point group gets the highest symmetry within a crystal_system
       if(str_sp.fgroup_calculated==FALSE) str_sp.CalculateSymmetryFactorGroup(SYS_VERBOSE);
       if(str_sp.pgroup_xtal_calculated==FALSE) str_sp.CalculateSymmetryPointGroupCrystal(SYS_VERBOSE);
    */

    //DX - START
    if(str_sp.pgroup_calculated==FALSE || str_sp.fgroup_calculated==FALSE || str_sp.pgroup_xtal_calculated==FALSE) {
      ofstream FileDevNull("/dev/null");
      _aflags aflags;
      _kflags kflags;
      kflags.KBIN_SYMMETRY_PGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE=FALSE;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE=FALSE;
      aflags.Directory="./";aflags.QUIET=FALSE;
      str_sp.LatticeReduction_avoid=TRUE;
      if(LDEBUG){ cerr << "LATTICE::Standard_Lattice_Structure: [4b]" << endl;}
      pflow::PerformFullSymmetry(str_sp,FileDevNull,aflags,kflags,SYS_VERBOSE,cout);
    }

    //**************************JUNKAI EDITED START****************************
    /* OLD (DX)
    // JUNKAI CHECK THE POSCAR OF KESONG it has 1.21% strain and finds it cubic but is it tet ?
    //  cerr << "[1]" << endl;
    if(str_sp.pgroup_calculated==FALSE) {
    //   cerr << "str_sp.pgroup_calculated" << endl;
    ofstream FileDevNull("/dev/null");
    _aflags aflags;
    aflags.Directory="./";aflags.QUIET=TRUE;
    str_sp.LatticeReduction_avoid=TRUE;
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4a1]" << endl;
    //DX SYM::CalculatePointGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps/2.0);
    SYM::CalculatePointGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps); //DX
    }

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4b]" << endl;
    //  cerr << "[2]" << endl;
    if(str_sp.fgroup_calculated==FALSE) {
    //  cerr << "str_sp.fgroup_calculated" << endl;
    ofstream FileDevNull("/dev/null");
    _aflags aflags;
    aflags.Directory="./";aflags.QUIET=TRUE;
    str_sp.LatticeReduction_avoid=TRUE;
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4b1]" << endl;
    //DX SYM::CalculatePointGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps/2.0);
    SYM::CalculatePointGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps); //DX
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4b2]" << endl;
    //DX SYM::CalculateFactorGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps/2.0);
    SYM::CalculateFactorGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps); //DX
    }

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4c]" << endl;
    //  cerr << "[3]" << endl;
    if(str_sp.pgroup_xtal_calculated==FALSE) {
    //  cerr << "str_sp.pgroup_xtal_calculated" << endl;
    //   cerr << symeps << endl;
    ofstream FileDevNull("/dev/null");
    _aflags aflags;
    aflags.Directory="./";aflags.QUIET=TRUE;
    str_sp.LatticeReduction_avoid=TRUE;
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4c1]" << endl;
    //DX if(str_sp.pgroup_calculated==FALSE) SYM::CalculatePointGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps/2.0);
    if(str_sp.pgroup_calculated==FALSE) SYM::CalculatePointGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps); //DX
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4c2]" << endl;
    //DX if(str_sp.fgroup_calculated==FALSE) SYM::CalculateFactorGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps/2.0);
    if(str_sp.fgroup_calculated==FALSE) SYM::CalculateFactorGroup(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps); //DX
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4c3]" << endl;
    //DX SYM::CalculatePointGroupCrystal(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps/2.0);
    SYM::CalculatePointGroupCrystal(FileDevNull,str_sp,aflags,FALSE,SYS_VERBOSE,cout,symeps); //DX
    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4c4]" << endl;
    }

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4d]" << endl;
    //  cerr << "[4]" << endl;
    */ //OLD (DX)
    //DX - END
    //**************************JUNKAI EDITED END****************************

    // http://en.wikipedia.org/wiki/Crystal_system  table from bottom to top
    //  if(str_sp.pgroup.size()==24 && str_sp.crystal_system=="") str_sp.crystal_system="hexagonal";
    crystal_system=str_sp.crystal_system;
    // cerr << crystal_system << endl;exit(0);

    /* OLD
       if(str_sp.pgroup.size()==48) {crystal_system="cubic";lattice_found=TRUE;};
       // if(str_sp.pgroup.size()==24) {crystal_system="hexagonal";lattice_found=TRUE;}; // hexagonal and trigonal must suffer
       //  if(str_sp.pgroup.size()==12) {crystal_system="trigonal";lattice_found=TRUE;}; // hexagonal and trigonal must suffer
       if(str_sp.pgroup.size()==16) {crystal_system="tetragonal";lattice_found=TRUE;};
       if(str_sp.pgroup.size()== 8) {crystal_system="orthorhombic";lattice_found=TRUE;};
       if(str_sp.pgroup.size()== 4) {crystal_system="monoclinic";lattice_found=TRUE;};
       if(str_sp.pgroup.size()== 2) {crystal_system="triclinic";lattice_found=TRUE;};
       if(lattice_found==FALSE) { // if not found need more stuff
       //    cerr << "DEBUG: LATTICE::Standard_Lattice_Structure crystal_system not from point group lattice, need factor group" << endl;
       crystal_system=str_sp.crystal_family; // always good
       // if(LDEBUG) cout << str_sp.crystal_family << endl;
       // if(LDEBUG) cout << str_sp.crystal_system << endl;
       if(str_sp.crystal_family=="hexagonal" && str_sp.crystal_system=="hexagonal") crystal_system="hexagonal"; // 1 sixfold axis of rotation
       if(str_sp.crystal_family=="hexagonal" && str_sp.crystal_system=="trigonal")  crystal_system="trigonal"; // 1 threefold axis of rotation
       //  if(LDEBUG) cerr << "DEBUG: crystal_system = " << crystal_system << endl;
       // Since do not have a basis (equal to 1 atom in position 000)
       // then the max symmetry is given by the max order in the symmetry number in the class of
       // page http://en.wikipedia.org/wiki/crystal_system
       // pick the centrosymmetric with highest order
       }
    */

    if(LDEBUG) cerr << "LATTICE::Standard_Lattice_Structure: [4c]" << endl;

    if(LDEBUG) {
      cerr << "str_sp.atoms.size()=" << str_sp.atoms.size() << endl;
      cerr << "str_sc.atoms.size()=" << str_sc.atoms.size() << endl;
      cerr << "str_sp.pgroup.size()=" << str_sp.pgroup.size() << endl;
      cerr << "str_sp.fgroup.size()=" << str_sp.fgroup.size() << endl;
      cerr << "str_sp.pgroup_xtal.size()=" << str_sp.pgroup_xtal.size() << endl;
      cerr << "str_sp.crystal_family=" << str_sp.crystal_family << endl;
      cerr << "str_sp.crystal_system=" << str_sp.crystal_system << endl;
      cerr << "crystal_system=" << crystal_system << endl;
      cerr << "radius_sp=" << RadiusSphereLattice(str_sp.lattice) << endl;
      cerr << "radius_sc=" << RadiusSphereLattice(str_sc.lattice) << endl;
      //  exit(0);
    }
  
    cutoff=0.0;
    rdata=Getabc_angles(str_sc.lattice,DEGREES);
    if(crystal_system=="cubic") {cutoff=sqrt(2)*max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;}
    if(crystal_system=="hexagonal") {cutoff=max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;} // the primitives always get the max of the lenghts
    //  if(crystal_system=="tetragonal") {cutoff=(2.0/sqrt(3.0))*max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;}
    if(crystal_system=="tetragonal") {cutoff=1.05*(2.0/sqrt(3.0))*RadiusSphereLattice(str_sc.lattice);lattice_found=TRUE;}
    if(crystal_system=="trigonal") {cutoff=max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;} // the primitives always get the max of the lenghts
    //  if(crystal_system=="orthorhombic") {cutoff=sqrt(2)*max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;} //stefano recipe
    if(crystal_system=="orthorhombic") {cutoff=1.05*(2.0/sqrt(3.0))*RadiusSphereLattice(str_sc.lattice);lattice_found=TRUE;}
    if(crystal_system=="monoclinic") {cutoff=sqrt(2)*max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;} //sqrt(2) catches the centered face
    if(crystal_system=="triclinic") {cutoff=max(rdata[1],rdata[2],rdata[3]);lattice_found=TRUE;} // the primitives always get the max of the lenghts
    if(crystal_system=="triclinic") {cutoff=max(rdata[2]+rdata[3],rdata[1]+rdata[3],rdata[1]+rdata[2]);lattice_found=TRUE;} // the primitives always get the max of the lenghts
    //  if(crystal_system=="triclinic") {cutoff=rdata[1]+rdata[2]+rdata[3];lattice_found=TRUE;} // the primitives always get the max of the lenghts
  
    // cerr << crystal_system << endl;//exit(0);
    //  LDEBUG=TRUE;
    /* STEFANO OLD VERSION
       if(lattice_found==FALSE || crystal_system=="none" || cutoff<=0.01) {
       cerr << "ERROR: (LATTICE::Standard_Lattice_Structure) NO crystal_system FOUND !" << endl;
       cerr << " str_sp.pgroup.size()=" << str_sp.pgroup.size() << endl;
       cerr << " str_sp.fgroup.size()=" << str_sp.fgroup.size() << endl;
       cerr << " str_sp.pgroup_xtal.size()=" << str_sp.pgroup_xtal.size() << endl;
       exit(0);
       }
    */

    //**************************JUNKAI EDITED START****************************
    if((lattice_found==FALSE || crystal_system=="none" || cutoff<=0.01)&&!histogram) {
      vector<double> tolerance;
      int ij;
      double aeps;
      tolerance=PointGroupHistogramCheck(str_sp);
      ij=(int)tolerance.size()-time-1;
      time++;
      aeps=tolerance.at(ij);
      // cerr << "RETRY tolerance.size=" << tolerance.size() << "  aeps=" << aeps << endl;
      return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.05,0.5,time,aeps);
      //   return LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.05,0.5,time,aeps);
    }

    //**************************JUNKAI EDITED END****************************  

    // crystal_system="all";cutoff=sqrt(2)*max(rdata[1],rdata[2],rdata[3]);

    for(int i=-dims[1];i<=dims[1];i+=1)
      for(int j=-dims[2];j<=dims[2];j+=1)
	for(int k=-dims[3];k<=dims[3];k+=1) {
	  aus=i*str_sp.lattice(1)+j*str_sp.lattice(2)+k*str_sp.lattice(3);
	  if(modulus(aus)<=cutoff+0.01) vvectors.push_back(aus);   // cubics are all in nn distances
	}
  
    // most of the time is spent inside this routine
    for(uint i=0;i<vvectors.size();i++) { for(uint ii=1;ii<=3;ii++) rlattice[1][ii]=vvectors[i][ii];
      for(uint j=0;j<vvectors.size();j++)  { for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=vvectors[j][ii];
	for(uint k=0;k<vvectors.size();k++)  { for(uint ii=1;ii<=3;ii++) rlattice[3][ii]=vvectors[k][ii];
	  // vaus=det(vvectors[i],vvectors[j],vvectors[k]);                                                 // SLOW
	  vaus=rlattice[1][1]*rlattice[2][2]*rlattice[3][3]+rlattice[1][2]*rlattice[2][3]*rlattice[3][1]+        // FAST
	    rlattice[1][3]*rlattice[2][1]*rlattice[3][2]-rlattice[1][3]*rlattice[2][2]*rlattice[3][1]-        // FAST
	    rlattice[1][2]*rlattice[2][1]*rlattice[3][3]-rlattice[1][1]*rlattice[2][3]*rlattice[3][2];        // FAST
	  if(abs(vaus-1.0*volume)<eps || abs(vaus-2.0*volume)<eps || abs(vaus-4.0*volume)<eps) {
	    rdata=Getabc_angles(rlattice,DEGREES);a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    if(abs(vaus-volume)<eps) {
	      vrlattice1.push_back(rlattice);
	      if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps))
		vrlattice1_aaa.push_back(rlattice);
	    }
	  }
	}
      }
    }
  
    if(LDEBUG || 0) {
      cerr << "DEBUG dims=" << dims << endl;
      cerr << "DEBUG vvectors.size()=" << vvectors.size() << endl;
      cerr << "DEBUG vrlattice1.size()=" << vrlattice1.size() << "   vrlattice1_aaa.size()=" << vrlattice1_aaa.size() << " " << endl;
      cerr << "DEBUG str_sp.pgroup.size()=" << str_sp.pgroup.size() << endl;
      cerr << "DEBUG str_sc.pgroup.size()=" << str_sc.pgroup.size() << endl;
      cerr << "DEBUG str_sp.fgroup.size()=" << str_sp.fgroup.size() << endl;
      cerr << "DEBUG str_sc.fgroup.size()=" << str_sc.fgroup.size() << endl;
      cerr << "DEBUG crystal_system=" << crystal_system << endl;
    }
  

    // ***************************************************************************
    // start scan
    found=FALSE;

    // ***************************************************************************************************
    bool VERBOSE_PROGRESS=FALSE; // TRUE;
    //  bool VERBOSE_PROGRESS=TRUE;
    // try different solutions
    for (uint choice=0;choice<=2&&!found;choice++) {
      if(VERBOSE_PROGRESS) if(choice==0) cerr << "LATTICE::Standard_Lattice_Structure: DEFINITION PRISTINE" << endl;
      if(VERBOSE_PROGRESS) if(choice==1) cerr << "LATTICE::Standard_Lattice_Structure: DEFINITION RELAX1" << endl; // locura
      if(VERBOSE_PROGRESS) if(choice==2) cerr << "LATTICE::Standard_Lattice_Structure: DEFINITION RELAX2" << endl; // locura
      if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: crystal_system=" << crystal_system << endl;
      
      // ***************************************************************************************************
      // CUBIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="cubic")) {
	// --------------------------------------------------------------------------
	// CUB FROM PRIMITIVE
	// CUBIC SYSTEMS looking for lattices a a a 90 90 90
	if(VERBOSE_PROGRESS) cerr << "CUBIC SYSTEMS looking for lattices a a a 90 90 90" << endl;
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                           // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if((choice==0 && (aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps)))) { // just to help not necessary
	    // CUB: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 a]: a a a 90 90 90  primitive
	    // CUB: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 a]: a a a 90 90 90  conventional
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps) &&
	       aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	      found=TRUE;
	      // do the STANDARD PRIMITIVE
	      str_sp.bravais_lattice_type="CUB";
	      str_sp.bravais_lattice_variation_type="CUB";
	      str_sp.bravais_lattice_system="cubic";
	      str_sp.pearson_symbol="cP"+aurostd::utype2string(str_sp.atoms.size());
	      ac=a; // make the standard primitive cubic
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive cubic
	      plattice[2][1]=0.0;plattice[2][2]= ac;plattice[2][3]=0.0; // make the standard primitive cubic
	      plattice[3][1]=0.0;plattice[3][2]=0.0;plattice[3][3]= ac; // make the standard primitive cubic
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // CUB standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE CUB
	// --------------------------------------------------------------------------
	// BCC FROM PRIMITIVE
	// BCC SYSTEMS looking for lattices a a a 109.471 109.471 109.471
	if(VERBOSE_PROGRESS)  cerr << "BCC SYSTEMS looking for lattices a a a 109.471 109.471 109.471" << endl;
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                             // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if((choice==0 && (aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps)))) { // just to help not necessary
	    // BCC:  a1=[-a/2 a/2 a/2] a2=[a/2 -a/2 a/2]/ a3=[a/2 a/2 -a/2]     a a a 109.471 109.471 109.471 primitive
	    // BCC:  [a 0 0] [0 a 0] [0 0 a] with 0,0,0 and 1/2 1/2 1/2 basis
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps) &&
	       aurostd::isequal(alpha,109.471,epsang) && aurostd::isequal(beta,109.471,epsang) && aurostd::isequal(gamma,109.471,epsang)) {
	      found=TRUE;
	      str_sp.bravais_lattice_type="BCC";
	      str_sp.bravais_lattice_variation_type="BCC";
	      str_sp.bravais_lattice_system="cubic";
	      str_sp.pearson_symbol="cI"+aurostd::utype2string(2*str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=2.0*a/sqrt(3.0); // make the standard primtive bcc
	      plattice[1][1]=-ac/2.0;plattice[1][2]= ac/2.0;plattice[1][3]= ac/2.0; // make the standard primtive bcc
	      plattice[2][1]= ac/2.0;plattice[2][2]=-ac/2.0;plattice[2][3]= ac/2.0; // make the standard primtive bcc
	      plattice[3][1]= ac/2.0;plattice[3][2]= ac/2.0;plattice[3][3]=-ac/2.0; // make the standard primtive bcc
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // make the STANDARD CONVENTIONAL
	      // C=supercell*P => supercell = C*inv(P)
	      str_sc=GetSuperCell(str_sp,0,1,1,1,0,1,1,1,0); // very cool (dual of bcc is fcc)
	      str_sc.bravais_lattice_type="CUB";
	      str_sc.bravais_lattice_variation_type="CUB";
	      str_sc.bravais_lattice_system="cubic";
	      str_sc.pearson_symbol="cP"+aurostd::utype2string(str_sc.atoms.size());
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE BCC
	// --------------------------------------------------------------------------
	// FCC FROM PRIMITIVE
	// FCC SYSTEMS looking for lattices a a a 60.0 60.0 60.0
	if(VERBOSE_PROGRESS)  cerr << "FCC SYSTEMS looking for lattices a a a 60.0 60.0 60.0" << endl;
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                           // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if((choice==0 && (aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps)))) { // just to help not necessary
	    rdata=Getabc_angles(rlattice,DEGREES);a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    // FCC:  a1=[0 a/2 a/2]  a2=[a/2 0 a/2]  a3=[a/2 a/2 0] a a a 60.0 60.0 60.0
	    // FCC:  [a 0 0] [0 a 0] [0 0 a] with 0,0,0 and 0 1/2 1/2 + 1/2 0 1/2 + 1/2 1/2 0 basis
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps) && aurostd::isequal(c,a,eps) &&
	       aurostd::isequal(alpha,60.0,epsang) && aurostd::isequal(beta,60.0,epsang) && aurostd::isequal(gamma,60.0,epsang)) {
	      found=TRUE;
	      str_sp.bravais_lattice_type="FCC";
	      str_sp.bravais_lattice_variation_type="FCC";
	      str_sp.bravais_lattice_system="cubic";
	      str_sp.pearson_symbol="cF"+aurostd::utype2string(4*str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=2.0*a/sqrt(3.0); // make the standard primtive fcc
	      plattice[1][1]=   0.0;plattice[1][2]=ac/2.0;plattice[1][3]=ac/2.0; // make the standard primtive fcc
	      plattice[2][1]=ac/2.0;plattice[2][2]=   0.0;plattice[2][3]=ac/2.0; // make the standard primtive fcc
	      plattice[3][1]=ac/2.0;plattice[3][2]=ac/2.0;plattice[3][3]=   0.0; // make the standard primtive fcc
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // make the STANDARD CONVENTIONAL
	      // C=supercell*P => supercell = C*inv(P)
	      str_sc=GetSuperCell(str_sp,-1,1,1,1,-1,1,1,1,-1); // very cool (dual of fcc is bcc)
	      str_sc.bravais_lattice_type="CUB";
	      str_sc.bravais_lattice_variation_type="CUB";
	      str_sc.bravais_lattice_system="cubic";
	      str_sc.pearson_symbol="cP"+aurostd::utype2string(str_sc.atoms.size());
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE FCC
	// --------------------------------------------------------------------------
      } // DONE CUBIC CRYSTAL SYSTEMS
      // ***************************************************************************************************
  
      // ***************************************************************************************************
      // TETRAGONAL CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="tetragonal")) {
	// cerr << "TETRAGONAL CRYSTAL SYSTEM SEARCH " << vrlattice1.size() << " " << eps << endl;
	// TET SYSTEM
	if(VERBOSE_PROGRESS) cerr << "TET SYSTEMS looking for lattices a a c 90.0 90.0 90.0" << endl;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  //   cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << endl;
	  if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    // cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << endl;
	    if((choice==0 && ((aurostd::isequal(a,b,eps)&&aurostd::isdifferent(b,c,eps))||(aurostd::isequal(b,c,eps)&&aurostd::isdifferent(c,a,eps))||(aurostd::isequal(c,a,eps)&&aurostd::isdifferent(a,b,eps)))) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps)&&aurostd::isequal(b,c,eps))))) {
	      // if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
    
	      // TET: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 c]: a a c 90 90 90  primitive
	      // TET: a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 c]: a a c 90 90 90  conventional
	      // ordering vectors to a a c
	      if(aurostd::isequal(a,b,eps)) {;}                                // a a c, ok
	      if(aurostd::isequal(a,c,eps)) {swap_rows(rlattice,2,3);for(uint ii=1;ii<=3;ii++) rlattice[1][ii]=-rlattice[1][ii];}  // a c a, swap 2,3
	      if(aurostd::isequal(b,c,eps)) {swap_rows(rlattice,1,3);for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=-rlattice[2][ii];}  // c a a, swap 1,3  
	      rdata=Getabc_angles(rlattice,DEGREES);
	      a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6]; // DONE
	      // TET a1=[a 0 0]  a2=[0 a 0]  a3=[0 0 c]: a a c 90 90 90
	      found=TRUE;
	      str_sp.bravais_lattice_type="TET";
	      str_sp.bravais_lattice_variation_type="TET";
	      str_sp.bravais_lattice_system="tetragonal";
	      str_sp.pearson_symbol="tP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;cc=c; // make the standard primitive tet
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive tet
	      plattice[2][1]=0.0;plattice[2][2]= ac;plattice[2][3]=0.0; // make the standard primitive tet
	      plattice[3][1]=0.0;plattice[3][2]=0.0;plattice[3][3]= cc; // make the standard primitive tet
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // do the
	      // TET standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE TET
	// --------------------------------------------------------------------------
	// BCT FROM PRIMITIVE (TYPE 2, SMART)
	// looking for rlattice that is a transformed version of a1=[-a/2 a/2 c/2] a2=[a/2 -a/2 c/2] a3=[a/2 a/2 -c/2]
	if(VERBOSE_PROGRESS) cerr << "BCT FROM PRIMITIVE looking for rlattice that is a transformed version of a1=[-a/2 a/2 c/2] a2=[a/2 -a/2 c/2] a3=[a/2 a/2 -c/2]" << endl;
	xmatrix<double> bct(3,3);
	bct[1][1]=0.0;bct[1][2]=1.0;bct[1][3]=1.0;
	bct[2][1]=1.0;bct[2][2]=0.0;bct[2][3]=1.0;
	bct[3][1]=1.0;bct[3][2]=1.0;bct[3][3]=0.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(bct*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of BCT then bct*rlattice==conventional lattice of BCT
	  // STD_CONV  [a 0 0,  0 a 0 , 0 0 c ]
	  // STD_PRIM: a1=[-a/2 a/2 c/2] a2=[a/2 -a/2 c/2] a3=[a/2 a/2 -c/2]
	  if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if((choice==0 && ((aurostd::isequal(a,b,eps)&&aurostd::isdifferent(b,c,eps)))) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps)&&aurostd::isequal(b,c,eps))))) {
	      ac=a;cc=c;
	      if(abs(ac*ac*cc-2.0*volume)<eps) {
		found=TRUE;
		str_sp.bravais_lattice_type="BCT";
		str_sp.bravais_lattice_system="tetragonal";
		str_sp.pearson_symbol="tI"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]=-ac/2.0;plattice[1][2]= ac/2.0;plattice[1][3]= cc/2.0; // make the standard primitive bct
		plattice[2][1]= ac/2.0;plattice[2][2]=-ac/2.0;plattice[2][3]= cc/2.0; // make the standard primitive bct
		plattice[3][1]= ac/2.0;plattice[3][2]= ac/2.0;plattice[3][3]=-cc/2.0; // make the standard primitive bct
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
		if(cc<ac) {str_sp.bravais_lattice_variation_type="BCT1";}
		if(ac<cc) {str_sp.bravais_lattice_variation_type="BCT2";}
		// do the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,bct);
		str_sc.bravais_lattice_type="TET";
		str_sc.bravais_lattice_variation_type="TET";
		str_sc.bravais_lattice_system="tetragonal";
		str_sc.pearson_symbol="tP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    } // DONE BCT
	  }
	}
      } // DONE TETRAGONAL CRYSTAL SYSTEMS
      // **************************************************************************************************
  
      // ***************************************************************************************************
      // TRIGONAL CRYSTAL SYSTEM                          
      if(found==FALSE && (crystal_system=="all" || crystal_system=="trigonal")) {  // trigonal can be rhl or hex
	if(VERBOSE_PROGRESS) cerr << "RHL SYSTEMS looking for lattices a a a alpha alpha alpha" << endl;
	double rhl_ratio=1.0; // you might need to reduce this to reproduce tiny distortions... insted of increasing eps, you reduce it
	xmatrix<double> rhl(3,3);
	rhl[1][1]=1.0;rhl[1][2]=-1.0;rhl[1][3]= 0.0;
	rhl[2][1]=0.0;rhl[2][2]= 1.0;rhl[2][3]=-1.0;
	rhl[3][1]=1.0;rhl[3][2]= 1.0;rhl[3][3]= 1.0;
	// RHL SYSTEMS looking for lattices a a a alpha alpha alpha
	for(uint i=0;i<vrlattice1_aaa.size()&&found==FALSE;i++) {  // FASTER with _aaa
	  rlattice=vrlattice1_aaa.at(i);                           // FASTER with _aaa
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  //  cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << endl;
	  if(choice==0) { // only pristine
	    if(aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps)) {  // aaa
	      //  cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << " " << epsang << endl;
	      if(aurostd::isequal(alpha,beta,epsang) && aurostd::isequal(beta,gamma,epsang)) {  // alpha alpha alpha
		//  cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << " " << epsang << endl;
		if(aurostd::isdifferent(alpha,90.0,epsang/rhl_ratio) && // no cubic
		   aurostd::isdifferent(alpha,60.0,epsang/rhl_ratio) && aurostd::isdifferent(alpha,180.0-60.0,epsang/rhl_ratio) && // no fcc (60 and 120)
		   aurostd::isdifferent(alpha,109.471,epsang/rhl_ratio) && aurostd::isdifferent(alpha,180.0-109.471,epsang/rhl_ratio)) {  // no bcc stuff
		  //   cerr << a << " " << b << " " << c << " " << alpha << " " << beta << " " << gamma << " " << eps << " " << epsang << endl;
		  // make standard primitive
		  // RHL: a1=[a*cos(alpha/2) -a*sin(alpha/2) 0]
		  // RHL: a2=[a*cos(alpha/2)  a*sin(alpha/2) 0]
		  // RHL: a3=[a*cos(alpha)/cos(alpha/2) 0 a*sqrt(1-cos(alpha)*cos(alpha)/cos(alpha/2)/cos(alpha/2))]
		  found=TRUE;
		  str_sp.bravais_lattice_type="RHL";
		  str_sp.bravais_lattice_system="rhombohedral";
		  str_sp.pearson_symbol="hR"+aurostd::utype2string(str_sp.atoms.size()); // FIX
		  if(alpha<90.0) {
		    str_sp.bravais_lattice_variation_type="RHL1";
		  } else {
		    str_sp.bravais_lattice_variation_type="RHL2";
		  }
		  alpha=alpha*pi/180.0; // change to RADIANTS !
		  plattice[1][1]=a*cos(alpha/2.0);plattice[1][2]=-a*sin(alpha/2.0);plattice[1][3]=0.0; // make the standard primitive rhl
		  plattice[2][1]=a*cos(alpha/2.0);plattice[2][2]= a*sin(alpha/2.0);plattice[2][3]=0.0; // make the standard primitive rhl
		  plattice[3][1]=a*cos(alpha)/cos(alpha/2.0);plattice[3][2]=0.0;
		  plattice[3][3]=a*sqrt(1.0-cos(alpha)*cos(alpha)/cos(alpha/2.0)/cos(alpha/2.0));      // make the standard primitive rhl
		  LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
		  // make the STANDARD CONVENTIONAL
#ifndef _RHL_HEX_SC_
		  // RHL standard_conventional=standard_primitive
		  str_sc=str_sp;
#else
		  // C=supercell*P => supercell = C*inv(P)
		  str_sc=GetSuperCell(str_sp,rhl); // very cool (make it hexagonal)
		  str_sc.bravais_lattice_type="HEX";
		  str_sc.bravais_lattice_variation_type="HEX";
		  str_sc.bravais_lattice_system="hexagonal";
		  str_sc.pearson_symbol="hP"+aurostd::utype2string(str_sc.atoms.size());
#endif
		  if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
		}
	      }
	    }
	  }
	} // DONE RHL  
	if(LDEBUG) cerr << "TRIGONAL CRYSTAL SYSTEM found=" << found << endl;
      } // DONE TRIGONAL CRYSTAL SYSTEMS
      // ***************************************************************************************************

      // ***************************************************************************************************
      // HEXAGONAL CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="hexagonal" || crystal_system=="trigonal")) {
	if(VERBOSE_PROGRESS) cerr << "HEX SYSTEMS looking for lattices a a c 90 90 120" << endl;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // HEX a1=[a/2 -a/2*sqrt(3) 0]  a2=[a/2 a/2*sqrt(3) 0]  a3=[0 0 c]   a a c 90 90 120 primitive
	  // HEX a1=[a/2 -a/2*sqrt(3) 0]  a2=[a/2 a/2*sqrt(3) 0]  a3=[0 0 c]   a a c 90 90 120 conventional
	  if(choice==0) { // only pristine
	    if((aurostd::isequal(a,b,eps) && aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang)) &&
	       (aurostd::isequal(gamma,120.0,epsang) || aurostd::isequal(gamma,60.0,epsang))) {
	      if(aurostd::isequal(gamma,60.0,epsang)) {
		for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=rlattice[2][ii]-rlattice[1][ii];  //  a1 (a2-a1) a3 are 90 90 120
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	      }
	      found=TRUE;
	      str_sp.bravais_lattice_type="HEX";
	      str_sp.bravais_lattice_variation_type="HEX";
	      str_sp.bravais_lattice_system="hexagonal";
	      str_sp.pearson_symbol="hP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;cc=c; // make the conventional hex
	      plattice[1][1]=ac/2.0;plattice[1][2]=-ac/2.0*sqrt(3.0);plattice[1][3]=0.0; // make the conventional hex
	      plattice[2][1]=ac/2.0;plattice[2][2]= ac/2.0*sqrt(3.0);plattice[2][3]=0.0; // make the conventional hex
	      plattice[3][1]=   0.0;plattice[3][2]=              0.0;plattice[3][3]= cc; // make the conventional hex
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // make standard conventional lattice
	      // HEX standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE HEX
	if(LDEBUG) cerr << "HEXAGONAL CRYSTAL SYSTEM found=" << found << endl;
      } // DONE HEXAGONAL CRYSTAL SYSTEMS
      // ***************************************************************************************************
    
      // ***************************************************************************************************
      // ORTHORHOMBIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="orthorhombic")) {
	if(VERBOSE_PROGRESS) cerr << "ORC SYSTEM looking for lattices a<=b<=c 90 90 90" << endl;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if((choice==0 && ((aurostd::isdifferent(a,b,eps) && aurostd::isdifferent(b,c,eps) && aurostd::isdifferent(c,a,eps)))) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps) || aurostd::isequal(b,c,eps) || aurostd::isequal(c,a,eps))))) {
	      // ORC: a1=[a 0 0]  a2=[0 b 0]  a3=[0 0 c]: a<b<c 90 90 90  primitive
	      // ORC: a1=[a 0 0]  a2=[0 b 0]  a3=[0 0 c]: a<b<c 90 90 90  conventional
	      // ordering vectors to get a<b<c
	      rdata=Getabc_angles(rlattice,DEGREES);
	      a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      if(aurostd::isequal(a,minabc,eps)) {;}                           // min * * , ok
	      if(aurostd::isequal(b,minabc,eps)) {                             // * min * , swap 1,2
		swap_rows(rlattice,1,2);if(det(rlattice)<0.0) {for(uint ii=1;ii<=3;ii++) rlattice[3][ii]=-rlattice[3][ii];}
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      }
	      if(aurostd::isequal(c,minabc,eps) || aurostd::isequal(a,maxabc,eps)) { // * * min , swap 1,3    OR    max * * , swap 1,3
		swap_rows(rlattice,1,3);if(det(rlattice)<0.0) {for(uint ii=1;ii<=3;ii++) rlattice[2][ii]=-rlattice[2][ii];}
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      }
	      if(aurostd::isequal(b,maxabc,eps)) {                             // * max * , swap 2,3
		swap_rows(rlattice,3,2);if(det(rlattice)<0.0) {for(uint ii=1;ii<=3;ii++) rlattice[1][ii]=-rlattice[1][ii];}
		rdata=Getabc_angles(rlattice,DEGREES);
		a=rdata[1];b=rdata[2];c=rdata[3];minabc=min(a,b,c);maxabc=max(a,b,c);
	      }
	      if(aurostd::isequal(c,maxabc,eps)) {;}                           // * * max , ok
	      // done
	      a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6]; // DONE
	      // ORC a1=[a 0 0]  a2=[0 b 0]  a3=[0 0 c]: a<b<c 90 90 90
	      found=TRUE;
	      str_sp.bravais_lattice_type="ORC";
	      str_sp.bravais_lattice_variation_type="ORC";
	      str_sp.bravais_lattice_system="orthorhombic";
	      str_sp.pearson_symbol="oP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;bc=b;cc=c; // make the standard primitive orc
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive orc
	      plattice[2][1]=0.0;plattice[2][2]= bc;plattice[2][3]=0.0; // make the standard primitive orc
	      plattice[3][1]=0.0;plattice[3][2]=0.0;plattice[3][3]= cc; // make the standard primitive orc
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // ORC standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE ORC
	// --------------------------------------------------------------------------
	// ORCF from PRIMITIVE
	if(VERBOSE_PROGRESS) cerr << "ORCF SYSTEM looking for orcf*rlattice==conventional lattice of ORCF" << endl;
	xmatrix<double> orcf(3,3);
	orcf[1][1]=-1.0;orcf[1][2]=1.0;orcf[1][3]=1.0;
	orcf[2][1]=1.0;orcf[2][2]=-1.0;orcf[2][3]=1.0;
	orcf[3][1]=1.0;orcf[3][2]=1.0;orcf[3][3]=-1.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(orcf*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of ORCF then orcf*rlattice==conventional lattice of ORCF !
	  // standard conventional:
	  // [ac 0 0,  0 bc 0 , 0 0 cc ]
	  // standard primitive:
	  // a1=[0 bc/2 cc/2]
	  // a2=[ac/2 0 cc/2]
	  // a3=[ac/2 bc/2 0]
	  if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if((choice==0 && (a<b && b<c)) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps) && b<c) || (a<b && aurostd::isequal(b,c,eps)))) ||
	       (choice==2 && ((aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps))))) {
	      ac=a;bc=b;cc=c;
	      if(abs(ac*bc*cc-4.0*volume)<eps) {
		found=TRUE;
		str_sp.bravais_lattice_type="ORCF";
		str_sp.bravais_lattice_system="orthorhombic";
		str_sp.pearson_symbol="oF"+aurostd::utype2string(4*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		double mismatch=pow(1.0/ac,2.0)/(pow(1.0/bc,2.0)+pow(1.0/cc,2.0));
		//  cerr << mismatch << " " << eps/5.0 << endl;  // DEBUG purposes
		if(mismatch>1+eps/5.0) {str_sp.bravais_lattice_variation_type="ORCF1";}
		if(mismatch<1-eps/5.0) {str_sp.bravais_lattice_variation_type="ORCF2";}
		if(aurostd::isequal(mismatch,1.0,eps/5.0)) {str_sp.bravais_lattice_variation_type="ORCF3";}
		plattice[1][1]=   0.0;plattice[1][2]=bc/2.0;plattice[1][3]=cc/2.0; // make the standard primitive orcf
		plattice[2][1]=ac/2.0;plattice[2][2]=   0.0;plattice[2][3]=cc/2.0; // make the standard primitive orcf
		plattice[3][1]=ac/2.0;plattice[3][2]=bc/2.0;plattice[3][3]=   0.0; // make the standard primitive orcf
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // plattice is a rotated version of rlattice
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,orcf); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="ORC";
		str_sc.bravais_lattice_variation_type="ORC";
		str_sc.bravais_lattice_system="orthorhombic";
		str_sc.pearson_symbol="oP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  }
	}
	// --------------------------------------------------------------------------
	// ORCI FROM PRIMITIVE SMART
	if(VERBOSE_PROGRESS) cerr << "ORCF SYSTEM looking for orci*rlattice==conventional lattice of ORCI" << endl;
	xmatrix<double> orci(3,3);
	orci[1][1]=0.0;orci[1][2]=1.0;orci[1][3]=1.0;
	orci[2][1]=1.0;orci[2][2]=0.0;orci[2][3]=1.0;
	orci[3][1]=1.0;orci[3][2]=1.0;orci[3][3]=0.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(orci*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of ORCI then orci*rlattice==conventional lattice of ORCI !
	  // ORCI
	  // Standard Conventional:
	  // [ac 0 0,  0 bc 0 , 0 0 cc ]
	  // Standard Primitive:
	  // a1=[-a/2 b/2 c/2]
	  // a2=[a/2 -b/2 c/2]
	  // a3=[a/2 b/2 -c/2]
	  if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if((choice==0 && (a<b && b<c)) ||
	       (choice==1 && ((aurostd::isequal(a,b,eps) && b<c) || (a<b && aurostd::isequal(b,c,eps)))) ||
	       (choice==2 && ((aurostd::isequal(a,b,eps) && aurostd::isequal(b,c,eps))))) {
	      ac=a;bc=b;cc=c;
	      if(abs(ac*bc*cc-2.0*volume)<eps) {
		found=TRUE;
		str_sp.bravais_lattice_type="ORCI";
		str_sp.bravais_lattice_variation_type="ORCI";
		str_sp.bravais_lattice_system="orthorhombic";
		str_sp.pearson_symbol="oI"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]=-ac/2.0;plattice[1][2]= bc/2.0;plattice[1][3]= cc/2.0; // make the standard primitive orci
		plattice[2][1]= ac/2.0;plattice[2][2]=-bc/2.0;plattice[2][3]= cc/2.0; // make the standard primitive orci
		plattice[3][1]= ac/2.0;plattice[3][2]= bc/2.0;plattice[3][3]=-cc/2.0; // make the standard primitive orci
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,orci); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="ORC";
		str_sc.bravais_lattice_variation_type="ORC";
		str_sc.bravais_lattice_system="orthorhombic";
		str_sc.pearson_symbol="oP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  } // DONE ORCI
	}
	// --------------------------------------------------------------------------
	// ORCC FROM PRIMITIVE SMART
	if(VERBOSE_PROGRESS) cerr << "ORCC SYSTEM looking for orcc*rlattice==conventional lattice of ORCC" << endl;
	xmatrix<double> orcc(3,3);
	orcc[1][1]= 1.0;orcc[1][2]=1.0;orcc[1][3]=0.0;
	orcc[2][1]=-1.0;orcc[2][2]=1.0;orcc[2][3]=0.0;
	orcc[3][1]= 0.0;orcc[3][2]=0.0;orcc[3][3]=1.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(orcc*rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  // if rlattice==primitive lattice of ORCC then orcc*rlattice==conventional lattice of ORCC !
	  // ORCC (a<b)
	  // [ac 0 0,  0 bc 0 , 0 0 cc ]      Conventional direct lattice:
	  // Direct Lattice (VASP and Bilbao):
	  // a1=[a/2 -b/2 0]
	  // a2=[a/2 b/2 0]
	  // a3=[0 0 c]
	  // pick only the case when a<b
	  if(aurostd::isequal(alpha,90.0,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if((choice==0 && (a<b)) || (choice==1 && (aurostd::isequal(a,b,eps)))) {
	      ac=a;bc=b;cc=c;
	      if(abs(ac*bc*cc-2.0*volume)<eps) {
		found=TRUE;
		str_sp.bravais_lattice_type="ORCC";
		str_sp.bravais_lattice_variation_type="ORCC";
		str_sp.bravais_lattice_system="orthorhombic";
		str_sp.pearson_symbol="oS"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]=ac/2.0;plattice[1][2]=-bc/2.0;plattice[1][3]=0.0; // make the standard primitive orcc
		plattice[2][1]=ac/2.0;plattice[2][2]= bc/2.0;plattice[2][3]=0.0; // make the standard primitive orcc
		plattice[3][1]=   0.0;plattice[3][2]=    0.0;plattice[3][3]= cc; // make the standard primitive orcc
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,orcc); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="ORC";
		str_sc.bravais_lattice_variation_type="ORC";
		str_sc.bravais_lattice_system="orthorhombic";
		str_sc.pearson_symbol="oP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  } // DONE ORCC
	}
      }
      // ***************************************************************************************************

      // ***************************************************************************************************
      // MONOCLINIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="monoclinic")) {
	// MCL SYSTEM
	// STANDARD: a b<=c alpha<90 90 90  (a b c can be whatever so wikipedia must be relaxed... but put alpha as close as possible to 90)
	// with abc ordering
	if(VERBOSE_PROGRESS) cerr << "MCL SYSTEM a b<=c alpha<~90 90 90" << endl;
	double MCL_alpha_max=0.0;
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if(alpha<90.0 && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if((choice==0 && b<c)|| (choice==1 && aurostd::isequal(b,c,eps))) {  // && not necessary aurostd::isdifferent(a,b,eps) && aurostd::isdifferent(b,c,eps) && aurostd::isdifferent(c,a,eps)
	      if(alpha>MCL_alpha_max) MCL_alpha_max=alpha;
	    }
	  }
	}
	// got alpha_max
	for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	  rlattice=vrlattice1.at(i);
	  rdata=Getabc_angles(rlattice,DEGREES);
	  a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	  if(aurostd::isequal(alpha,MCL_alpha_max,epsang)  && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	    if((choice==0 && b<c)|| (choice==1 && aurostd::isequal(b,c,eps))) {  // && not necessary aurostd::isdifferent(a,b,eps) && aurostd::isdifferent(b,c,eps) && aurostd::isdifferent(c,a,eps)
	      // MCL: a1=[a 0 0]  a2=[0 b 0]  a3=[0 c*cos(alpha) c*sin(alpha)]: a b<c alpha<90 90 90  primitive
	      // MCL: a1=[a 0 0]  a2=[0 b 0]  a3=[0 c*cos(alpha) c*sin(alpha)]: a b<c alpha<90 90 90  conventional
	      found=TRUE;
	      str_sp.bravais_lattice_type="MCL";
	      str_sp.bravais_lattice_variation_type="MCL";
	      str_sp.bravais_lattice_system="monoclinic";
	      str_sp.pearson_symbol="mP"+aurostd::utype2string(str_sp.atoms.size());
	      str_sp.lattice=rlattice;
	      ac=a;bc=b;cc=c; // make the standard primitive mcl
	      plattice[1][1]= ac;plattice[1][2]=0.0;plattice[1][3]=0.0; // make the standard primitive mcl
	      plattice[2][1]=0.0;plattice[2][2]= bc;plattice[2][3]=0.0; // make the standard primitive mcl
	      plattice[3][1]=0.0;plattice[3][2]=cc*cos(alpha/180.0*pi);plattice[3][3]=cc*sin(alpha/180.0*pi); // make the standard primitive mcl
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);
	      // do the
	      // MCL standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  }
	} // DONE MCL
	// --------------------------------------------------------------------------
	// MCLC from PRIMITIVE
	// MCLC (a b c can be whatever so wikipedia must be relaxed... but put alpha as close as possible to 90)
	if(VERBOSE_PROGRESS) cerr << "MCLC SYSTEM a b c alpha<~90 90 90" << endl;
	xmatrix<double> mclc(3,3);
	mclc[1][1]=1.0;mclc[1][2]=-1.0;mclc[1][3]=0.0;
	mclc[2][1]=1.0;mclc[2][2]=1.0;mclc[2][3]=0.0;
	mclc[3][1]=0.0;mclc[3][2]=0.0;mclc[3][3]=1.0;
	if(choice==0) { // only pristine
	  double MCLC_alpha_max=0.0,testphi,alpharad;
	  for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	    rlattice=vrlattice1.at(i);
	    rdata=Getabc_angles(mclc*rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    if(alpha<90.0 && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	      ac=a;bc=b;cc=c;
	      if(abs(det(mclc*rlattice)-2.0*volume)<eps) {
		if(alpha>MCLC_alpha_max) MCLC_alpha_max=alpha;
	      }
	    }
	  }
	  // got alpha_max
	  for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++) {
	    rlattice=vrlattice1.at(i);
	    rdata=Getabc_angles(mclc*rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    // if rlattice==primitive lattice of MCLC then mclc*rlattice==conventional lattice of MCLC !
	    // MCLC a b c alpha<90.0 90 90
	    // [ac 0 0] [0 bc 0] [0 cc*cos(alpha) cc*sin(alpha) ]    Conventional direct lattice:
	    // Direct Lattice (VASP and Bilbao):
	    // a1=[ a/2 b/2 0] a2=[-a/2 b/2 0] a3=[0 c*cos(alpha) c*sin(alpha)]
	    if(aurostd::isequal(alpha,MCLC_alpha_max,epsang) && aurostd::isequal(beta,90.0,epsang) && aurostd::isequal(gamma,90.0,epsang)) {
	      ac=a;bc=b;cc=c;
	      if(abs(det(mclc*rlattice)-2.0*volume)<eps) {
		str_sp.bravais_lattice_type="MCLC";
		str_sp.bravais_lattice_system="monoclinic";
		str_sp.pearson_symbol="mS"+aurostd::utype2string(2*str_sp.atoms.size());
		str_sp.lattice=rlattice;
		plattice[1][1]= ac/2.0;plattice[1][2]=bc/2.0;plattice[1][3]=0.0; // make the standard primitive mclc
		plattice[2][1]=-ac/2.0;plattice[2][2]=bc/2.0;plattice[2][3]=0.0; // make the standard primitive mclc
		plattice[3][1]=  0.0;
		plattice[3][2]=cc*cos(alpha*pi/180.0);
		plattice[3][3]=cc*sin(alpha*pi/180.0);                           // make the standard primitive mclc
		LATTICE::fix_sts_sp(str_sp,rlattice,plattice);                    // project out of rlattice and in plattice
		str_sp.FixLattices();
		klattice=str_sp.klattice;
		kdata=Getabc_angles(klattice,DEGREES);kalpha=kdata[4];kbeta=kdata[5];kgamma=kdata[6];
		if(aurostd::isequal(kgamma,90.0,epsang)) {
		  found=TRUE; str_sp.bravais_lattice_variation_type="MCLC2";
		} else {
		  if(kgamma>90.0) {
		    found=TRUE; str_sp.bravais_lattice_variation_type="MCLC1";
		  } else {
		    // kgamma<90
		    alpharad=alpha*PI/180.0;
		    testphi=0.75-0.25*b*b/a/a-0.25*b/tan(alpharad)*(1/c/sin(alpharad)-1/b/tan(alpharad));
		    if(aurostd::isequal(testphi,0.5,eps)) {
		      found=TRUE; str_sp.bravais_lattice_variation_type="MCLC4";
		    } else {
		      if(testphi>0.5) {
			found=TRUE; str_sp.bravais_lattice_variation_type="MCLC3";
		      } else {
			found=TRUE; str_sp.bravais_lattice_variation_type="MCLC5";
		      }
		    }
		  }
		}
		// make the STANDARD CONVENTIONAL
		// C=supercell*P => supercell = C*inv(P)
		str_sc=GetSuperCell(str_sp,mclc); // very cool (dual of fcc is bcc)
		str_sc.bravais_lattice_type="MCL";
		str_sc.bravais_lattice_variation_type="MCL";
		str_sc.bravais_lattice_system="monoclinic";
		str_sc.pearson_symbol="mP"+aurostd::utype2string(str_sc.atoms.size());
		if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	      }
	    }
	  }
	}
      } // DONE MCLC
      // ***************************************************************************
      // TRICLINIC CRYSTAL SYSTEM
      if(found==FALSE && (crystal_system=="all" || crystal_system=="triclinic"))  {
	found=FALSE;
	double weight_TRI1a=1e6,weight_TRI2a=1e6,weight_TRI1b=1e6,weight_TRI2b=1e6;
	xmatrix<double> rrlattice_TRI1a(3,3),rrlattice_TRI2a(3,3),rrlattice_TRI1b(3,3),rrlattice_TRI2b(3,3);
	if(choice==0) { // only pristine
	  for(uint i=0;i<vrlattice1.size()&&found==FALSE;i++)  {
	    rlattice=vrlattice1.at(i);
	    rdata=Getabc_angles(rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
	    klattice=ReciprocalLattice(rlattice);
	    kdata=Getabc_angles(klattice,DEGREES);
	    ka=kdata[1];kb=kdata[2];kc=kdata[3];kalpha=kdata[4];kbeta=kdata[5];kgamma=kdata[6];
	    // find which better fits
	    if(kalpha>90.0 && kbeta>90.0 && kgamma>90.0) {
	      if((kalpha-90.0)+(kbeta-90.0)+(kgamma-90.0)<weight_TRI1a) {
		weight_TRI1a=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI1a=rlattice; } }
	    if(kalpha>90.0 && kbeta>90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
	      if((kalpha-90.0)+(kbeta-90.0)<weight_TRI2a) {
		weight_TRI2a=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI2a=rlattice; } }
	    if(kalpha<90.0 && kbeta<90.0 && kgamma<90.0) {
	      if((90.0-kalpha)+(90.0-kbeta)+(90.0-kgamma)<weight_TRI1b) {
		weight_TRI1b=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI1b=rlattice; } }
	    if(kalpha<90.0 && kbeta<90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
	      if((90.0-kalpha)+(90.0-kbeta)<weight_TRI2b) {
		weight_TRI2b=abs(kalpha-90.0)+abs(kbeta-90.0)+abs(kgamma-90.0);rrlattice_TRI2b=rlattice; } }
	  }
	  // now it found it
	  {
	    // lets pick the nicest TRI 1/2 a/b
	    a=b=c=alpha=beta=gamma=ka=kb=kc=kalpha=kbeta=kgamma=0.0;
	    bool found_TRI=FALSE;
	    if(weight_TRI1a<=min(weight_TRI2a,weight_TRI1b,weight_TRI2b) && found_TRI==FALSE) { rlattice=rrlattice_TRI1a; found_TRI=TRUE; }
	    if(weight_TRI2a<=min(weight_TRI1b,weight_TRI2b,weight_TRI1a) && found_TRI==FALSE) { rlattice=rrlattice_TRI2a; found_TRI=TRUE; }
	    if(weight_TRI1b<=min(weight_TRI2b,weight_TRI1a,weight_TRI2a) && found_TRI==FALSE) { rlattice=rrlattice_TRI1b; found_TRI=TRUE; }
	    if(weight_TRI2b<=min(weight_TRI1a,weight_TRI2a,weight_TRI1b) && found_TRI==FALSE) { rlattice=rrlattice_TRI2b; found_TRI=TRUE; }
	    rdata=Getabc_angles(rlattice,DEGREES);
	    a=rdata[1];b=rdata[2];c=rdata[3];alpha=rdata[4];beta=rdata[5];gamma=rdata[6];
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //2-17-16 COses, Karol Jarolimek discovered B1Cd1Li1O3_ICSD_200615 was calculated to be TRI1B, not TRI1A like shown in
            //Comp. Mat. Sci. 49, 299 (2010) Fig 47.  Turns out we forgot to set ka...kgamma, so the following lines were added.
            //Interestingly, this system has weight_TRIa==weight_TRIb, so the check above (TRI1a before TRI1b) is very important.
            klattice=ReciprocalLattice(rlattice);
            kdata=Getabc_angles(klattice,DEGREES);
            ka=kdata[1];kb=kdata[2];kc=kdata[3];kalpha=kdata[4];kbeta=kdata[5];kgamma=kdata[6];	    
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // picked, put the label
	    if(kalpha>90.0 && kbeta>90.0 && kgamma>90.0) {
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI1a";
	    }
	    if(kalpha>90.0 && kbeta>90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI2a";
	    }
	    if(kalpha<90.0 && kbeta<90.0 && kgamma<90.0) {
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI1b";
	    }
	    if(kalpha<90.0 && kbeta<90.0 && aurostd::isequal(kgamma,90.0,epsang)) {
	      found=TRUE;str_sp.bravais_lattice_type="TRI";str_sp.bravais_lattice_variation_type="TRI2b";
	    }
	    // found=false;
	    // cerr << kalpha << " " << kbeta << " " << kgamma << " " << found << endl;
	    if(found) {
	      str_sp.bravais_lattice_system="triclinic";
	      str_sp.pearson_symbol="aP"+aurostd::utype2string(str_sp.atoms.size());
	      // make standard primitive
	      double alf=alpha*PI/180.0, bet=beta*PI/180.0, gam=gamma*PI/180.0;
	      plattice[1][1]=a; plattice[1][2]=0.0; plattice[1][3]=0.0;
	      plattice[2][1]=b*cos(gam); plattice[2][2]=b*sin(gam); plattice[2][3]=0.0;
	      plattice[3][1]=c*cos(bet);
	      plattice[3][2]=c*(cos(alf)-cos(bet)*cos(gam))/sin(gam);
	      plattice[3][3]=1.0-cos(alf)*cos(alf)-cos(bet)*cos(bet)-cos(gam)*cos(gam)+2*cos(alf)*cos(bet)*cos(gam);
	      plattice[3][3]=c*sqrt(plattice[3][3])/sin(gam);
	      LATTICE::fix_sts_sp(str_sp,rlattice,plattice);  // project out of rlattice and in plattice
	      str_sp.FixLattices();
	      // TRI standard_conventional=standard_primitive
	      str_sc=str_sp;
	      if(VERBOSE_PROGRESS) cerr << str_sp.bravais_lattice_type << " found" << endl;
	    }
	  } //for vrlattice1
	} // choice
      }// TRIclinic DONE  
      // ***************************************************************************************************
      if(VERBOSE_PROGRESS) if(found) cerr << "str_sp.bravais_lattice_type=" << str_sp.bravais_lattice_type << endl;
      if(VERBOSE_PROGRESS) if(found) cerr << "str_sp.bravais_lattice_system=" << str_sp.bravais_lattice_system << endl;
      if(VERBOSE_PROGRESS) if(found) cerr << "str_sp.pearson_symbol=" << str_sp.pearson_symbol << endl;
    }
    // cerr << str_sp.bravais_lattice_type << endl; exit(0);
  
    if(found==FALSE) str_sp.bravais_lattice_type="UNKNOWN";
    if(found==TRUE) {
      str_sp.title=str_sp.title+"  ["+str_sp.bravais_lattice_type+","+str_sp.bravais_lattice_variation_type+","+str_sp.pearson_symbol+"]"+" (STD_PRIM doi:10.1016/j.commatsci.2010.05.010)";
    }
    if(found==FALSE) str_sc.bravais_lattice_type="UNKNOWN";
    if(found==TRUE) {
      str_sc.title=str_sc.title+"  ["+str_sc.bravais_lattice_type+","+str_sc.bravais_lattice_variation_type+","+str_sc.pearson_symbol+"]"+" (STD_CONV doi:10.1016/j.commatsci.2010.05.010)";
    }
    // DONE

    if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: X1 found=" << found << endl;
    
    if(found==FALSE) {
      str_sc.Standard_Lattice_calculated=TRUE;str_sc.Standard_Lattice_avoid=FALSE;
      str_sc.Standard_Lattice_primitive=FALSE;str_sc.Standard_Lattice_conventional=FALSE;
      str_sc.Standard_Lattice_has_failed=TRUE;
      str_sp.Standard_Lattice_calculated=TRUE;str_sp.Standard_Lattice_avoid=FALSE;
      str_sp.Standard_Lattice_primitive=FALSE;str_sp.Standard_Lattice_conventional=FALSE;
      str_sp.Standard_Lattice_has_failed=TRUE;
      str_sp=str_in;
      str_sp.GetPrimitive();
      str_sp.MinkowskiBasisReduction(); // shorten the vectors as much as possible and as perpendicular as possible
      str_sc=str_in; // copy it
      str_sc.MinkowskiBasisReduction(); // shorten the vectors as much as possible and as perpendicular as possible
    
      //    if(str_in.title!="NO_RECURSION")
      if(0)      {
	cerr << "REMAKING LATTICE" << endl;
	xstructure str_in2(str_in),str_sp2(str_sp),str_sc2(str_sc);
	while(str_sp2.atoms.size()>1) str_sp2.RemoveAtom(0);
	while(str_sc2.atoms.size()>1) str_sc2.RemoveAtom(0);
	str_in2.title="NO_RECURSION";
	LATTICE::Standard_Lattice_Structure(str_in2,str_sp2,str_sc2,eps,epsang,time,symeps,histogram);
	cerr << str_sc2.bravais_lattice_type << endl;
	cerr << str_sc2.bravais_lattice_system << endl;
	cerr << str_sc2.pearson_symbol << endl;
	exit(0);
      }
    }

    if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: X2 found=" << found << endl;

    if(found==TRUE) {
      str_sc.Standard_Lattice_calculated=TRUE;str_sc.Standard_Lattice_avoid=FALSE;
      str_sc.Standard_Lattice_primitive=FALSE;str_sc.Standard_Lattice_conventional=TRUE;
      str_sc.Standard_Lattice_has_failed=FALSE;
      str_sp.Standard_Lattice_calculated=TRUE;str_sp.Standard_Lattice_avoid=FALSE;
      str_sp.Standard_Lattice_primitive=TRUE;str_sp.Standard_Lattice_conventional=FALSE;
      str_sp.Standard_Lattice_has_failed=FALSE;
      str_sp.FixLattices();str_sc.FixLattices();
      // get reciprocal lattice
      // cerr << "STEP" << endl;
      if(0 && str_in.title!="NO_RECURSION") {
	// RECIPROCAL
	xstructure str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc;
	str_reciprocal_in.lattice=str_sp.klattice;str_reciprocal_in.FixLattices();
	str_reciprocal_in.title="NO_RECURSION";
	_atom atom;str_reciprocal_in.AddAtom(atom);
	// LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang); //STEFANO OLD VERSION
	int ss=0; //JUNKAI
	LATTICE::Standard_Lattice_Structure(str_reciprocal_in,str_reciprocal_sp,str_reciprocal_sc,eps,epsang,ss,_EPS_);//JUNKAI
	str_sp.reciprocal_lattice_type=str_reciprocal_sp.bravais_lattice_type;
	str_sp.reciprocal_lattice_variation_type=str_reciprocal_sp.bravais_lattice_variation_type;//wahyu mod
	//str_sp.reciprocal_conventional_lattice_type=str_reciprocal_sp.bravais_lattice_system;
	str_sc.reciprocal_lattice_type=str_reciprocal_sp.bravais_lattice_type;
	str_sc.reciprocal_lattice_variation_type=str_reciprocal_sc.bravais_lattice_variation_type;//wahyu mod
	//str_sc.reciprocal_conventional_lattice_type=str_reciprocal_sp.bravais_lattice_system;
	// SUPERLATTICE
	xstructure str_superlattice_in,str_superlattice_sp,str_superlattice_sc;
	str_superlattice_in=str_sp;
	str_superlattice_in.title="NO_RECURSION";
	str_superlattice_in.IdenticalAtoms();
	str_superlattice_in.Minkowski_calculated=FALSE;
	//LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,eps,epsang); //STEFANO OLD VERSION
	ss=0; //JUNKAI
	LATTICE::Standard_Lattice_Structure(str_superlattice_in,str_superlattice_sp,str_superlattice_sc,eps,epsang,ss,_EPS_);//JUNKAI
	str_sp.bravais_superlattice_type=str_superlattice_sp.bravais_lattice_type;
	str_sp.bravais_superlattice_variation_type=str_superlattice_sp.bravais_lattice_variation_type;
	str_sp.bravais_superlattice_system=str_superlattice_sp.bravais_lattice_system;
	str_sp.pearson_symbol_superlattice=str_superlattice_sp.pearson_symbol;
	str_sc.bravais_superlattice_type=str_superlattice_sp.bravais_lattice_type;
	str_sc.bravais_superlattice_variation_type=str_superlattice_sp.bravais_lattice_variation_type;
	str_sc.bravais_superlattice_system=str_superlattice_sp.bravais_lattice_system;
	str_sc.pearson_symbol_superlattice=str_superlattice_sp.pearson_symbol;
      }

      str_sp.SetVolume(str_sp_volume);
      str_sc.SetVolume(str_sp_volume*((double) str_sc.atoms.size()/str_sp.atoms.size()));
      str_sp.neg_scale=str_sp_neg_scale;  // reload from backup
      str_sc.neg_scale=str_sp_neg_scale;  // reload from backup
    }

    if(VERBOSE_PROGRESS) cerr << "LATTICE::Standard_Lattice_Structure: X3 found=" << found << endl;

    // last checks
    /* STEFANO OLD VERSION
       if(abs(NearestNeighbour(str_in)-NearestNeighbour(str_sp))>0.01) {
       cerr << "ERROR:  LATTICE::Standard_Lattice_Structure" << endl;
       cerr << "ERROR:  str_in.title=" << str_in.title << endl;
       cerr << "ERROR:  nn(str_in)=" << NearestNeighbour(str_in) << endl;
       cerr << "ERROR:  nn(str_sp)=" << NearestNeighbour(str_sp) << endl;
       cerr << "ERROR:  nn(str_sc)=" << NearestNeighbour(str_sc) << endl;
       // exit(0);
       }
    */  
    // some coding for test
    //   if(0) {
    //     _aflags aflags;APENNSY_Parameters params;params.LoadStructuresHQTC(aflags);
    //     for(uint i=1;i<params.structures[2].size();i++) {
    //       cout << "echo " << aurostd::RemoveCharacter(params.structures[2].at(i).name,'/') << endl;
    //       cout << "./aflow --proto " << aurostd::RemoveCharacter(params.structures[2].at(i).name,'/') << " Fe Co | ./aflow --sprim >> /dev/null " << endl;
    //     }
    //   }

    //  cerr << "OK=" << str_sp.bravais_lattice_type << endl;
    return found;
  }
}

void CheckLatticeHistogram(xstructure& a) {
  xstructure str_in(a),str_sp,str_sc;
  int ss=0;
  vector<double> tolerance;
  tolerance=PointGroupHistogramCheck(str_in);
  for(uint i=0;i<tolerance.size();i++) {
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,0.05,0.5,ss,tolerance.at(i),true);
    cout<<"TOLERANCE: "<<tolerance.at(i);
    cout<<"  LATTICETYPE: "<<str_sp.bravais_lattice_type<<endl;
  }
}

void CheckLatticeHistogram() {
  xstructure a(cin,IOAFLOW_AUTO);
  CheckLatticeHistogram(a);
}

namespace LATTICE {
  string Primitive_Lattice_Type(const xstructure& str) {
    xstructure str_in(str),str_sp,str_sc;
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc);
    return str_sp.bravais_lattice_type;
  }
}
/*
  namespace LATTICE {
  string Conventional_Lattice_Type(const xstructure& str) {
  xstructure str_in(str),str_sp,str_sc;
  LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc);
  return str_sp.bravais_lattice_system;
  }
  }
*/

namespace LATTICE {
  string Lattice_Variation_Type(const xstructure& str) {
    xstructure str_in(str),str_sp,str_sc;
    return str_sp.bravais_lattice_variation_type;
  }
}

namespace LATTICE {
  xstructure Standard_Primitive_Lattice_Structure(const xstructure& str) {
    xstructure str_in(str),str_sp,str_sc;
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc);
    return str_sp;
  }
}

namespace LATTICE {
  xstructure Standard_Conventional_Lattice_Structure(const xstructure& str) {
    xstructure str_in(str),str_sp,str_sc;
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc);
    return str_sc;
  }
}

// ***************************************************************************
namespace LATTICE {
  bool Bravais_Lattice_Structure(xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang) {
    //  LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang);// STEFANO OLD VERSION
    int ss=0; //JUNKAI                                                                                                                                                                                                                                                
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang,ss,_EPS_);//JUNKAI
    str_in.bravais_lattice_type=str_sp.bravais_lattice_type;
    str_in.bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
    str_in.bravais_lattice_system=str_sp.bravais_lattice_system;
    xstructure _str_in,_str_sp,_str_sc;
    _str_in.lattice=str_in.lattice; _str_in.scale=1.0;
    _atom atom; atom.cpos.clear();atom.fpos.clear();atom.type=0; _str_in.AddAtom(atom);
    //  LATTICE::Standard_Lattice_Structure(_str_in,str_sp,str_sc,eps,epsang); //STEFANO OLD VERSION
    ss=0; //JUNKAI                                        
    LATTICE::Standard_Lattice_Structure(_str_in,str_sp,str_sc,eps,epsang,ss,_EPS_);//JUNKAI
    str_in.bravais_lattice_lattice_type=_str_sp.bravais_lattice_type;
    str_in.bravais_lattice_lattice_variation_type=_str_sp.bravais_lattice_variation_type;
    str_in.bravais_lattice_lattice_system=_str_sp.bravais_lattice_system;
    return TRUE;
  }
}

//DX - START 
namespace LATTICE {
  bool Bravais_Lattice_StructureDefault(xstructure& str_in,xstructure& str_sp,xstructure& str_sc, bool full_sym) {
    return Bravais_Lattice_StructureDefault_20170401(str_in,str_sp,str_sc,full_sym);
  }
}
namespace LATTICE {
  bool Bravais_Lattice_StructureDefault_20170401(xstructure& str_in,xstructure& str_sp,xstructure& str_sc, bool full_sym) {
    bool same_eps = false;
    bool no_scan = false;
    bool ignore_checks = false;
    uint count=0;
    while(same_eps == false && count++ < 100){
      if(ignore_checks==true){
        same_eps = true; //force single while loop, no check
      }
      if(!LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc,full_sym)){
        cerr << "LATTICE::WARNING: Could not find crystal lattice type." << endl;
        if(!SYM::change_tolerance(str_sp,str_sp.sym_eps,str_sp.dist_nn_min,no_scan)){
          cerr << "LATTICE::WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << "[dir=" << str_sp.directory << "]" << endl;
          ignore_checks = true;
        }
        str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = str_sp.sym_eps;
        str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = str_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
        continue;
      }
      str_in.bravais_lattice_type=str_sp.bravais_lattice_type;
      str_in.bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
      str_in.bravais_lattice_system=str_sp.bravais_lattice_system;
      str_sp.bravais_lattice_type=str_sp.bravais_lattice_type;
      str_sp.bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
      str_sp.bravais_lattice_system=str_sp.bravais_lattice_system;
      xstructure _str_in,_str_sp,_str_sc;
      _str_in.lattice=str_in.lattice; _str_in.scale=1.0;
      _str_in.sym_eps=_str_sp.sym_eps=_str_sc.sym_eps=str_sp.sym_eps; //DX
      _str_in.sym_eps_calculated=_str_sp.sym_eps_calculated=_str_sc.sym_eps_calculated=str_sp.sym_eps_calculated; //DX
      _str_in.sym_eps_change_count=_str_sp.sym_eps_change_count=_str_sc.sym_eps_change_count=str_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
      _atom atom; atom.cpos.clear();atom.fpos.clear();atom.type=0; _str_in.AddAtom(atom);
      //DX 8/14/17 - START - Use real pgroup to calculate pgroupk and then set pgrouk from str_sp to the pgroup and pgroup_xtal of str_reciprocal_in
      //DX 8/14/17 The pgroup and pgroup_xtal are the same for the str_reciprocal structure because there is only one atom at the origin
      //DX 8/14/17 (i.e. lattice and crystal symmetry are the same for the reciprocal space crystal) 
      _str_in.pgroup = _str_sp.pgroup = _str_sc.pgroup = str_sp.pgroup;
      _str_in.pgroup_xtal = _str_sp.pgroup_xtal = _str_sc.pgroup_xtal = str_sp.pgroup;
      _str_in.pgroup_calculated = _str_sp.pgroup_calculated = _str_sc.pgroup_calculated = str_sp.pgroup_calculated;
      _str_in.pgroup_xtal_calculated = _str_sp.pgroup_xtal_calculated = _str_sc.pgroup_xtal_calculated = str_sp.pgroup_xtal_calculated;
      //DX 8/14/17 - END
      //DX 2/26/18 [OBSOLETE] if(!LATTICE::Standard_Lattice_StructureDefault(_str_in,_str_sp,_str_sc,full_sym)){;
      if(!LATTICE::Standard_Lattice_StructureDefault(_str_in,_str_sp,_str_sc,false)){; //DX 2/26/18 - do not need to do full sym on lattice
        cerr << "LATTICE::WARNING: Could not find lattice type." << endl;
        if(!SYM::change_tolerance(str_sp,str_sp.sym_eps,str_sp.dist_nn_min,no_scan)){
          cerr << "LATTICE::WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << "[dir=" << str_sp.directory << "]" << endl;
        }
        str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = str_sp.sym_eps;
        str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = str_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
        continue;
      }
      str_in.bravais_lattice_lattice_type=_str_sp.bravais_lattice_type;
      str_in.bravais_lattice_lattice_variation_type=_str_sp.bravais_lattice_variation_type;
      str_in.bravais_lattice_lattice_system=_str_sp.bravais_lattice_system;
      str_sp.bravais_lattice_lattice_type=_str_sp.bravais_lattice_type;
      str_sp.bravais_lattice_lattice_variation_type=_str_sp.bravais_lattice_variation_type;
      str_sp.bravais_lattice_lattice_system=_str_sp.bravais_lattice_system;
      if(str_sp.sym_eps == _str_sp.sym_eps){
        same_eps = true;
      }
      else{ 
        str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = _str_sp.sym_eps;
        str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = _str_sp.sym_eps_change_count; //DX 2/22/18 - added sym_eps change count
      }
    }
    if(count==100){
      cerr << "ERROR in Bravais_Lattice_StructureDefault(): Unable to find reliable sym_eps." << endl;
      return FALSE;
    }
    return TRUE;
  }
}
//DX - END

//DX - START
namespace LATTICE {
  bool Bravais_Lattice_StructureDefault_20160101(xstructure& str_in,xstructure& str_sp,xstructure& str_sc) {
    LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
    str_in.bravais_lattice_type=str_sp.bravais_lattice_type;
    str_in.bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
    str_in.bravais_lattice_system=str_sp.bravais_lattice_system;
    str_sp.bravais_lattice_type=str_sp.bravais_lattice_type;
    str_sp.bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
    str_sp.bravais_lattice_system=str_sp.bravais_lattice_system;
    xstructure _str_in,_str_sp,_str_sc;
    _str_in.lattice=str_in.lattice; _str_in.scale=1.0;
    _atom atom; atom.cpos.clear();atom.fpos.clear();atom.type=0; _str_in.AddAtom(atom);
    LATTICE::Standard_Lattice_StructureDefault(_str_in,_str_sp,_str_sc);
    str_in.bravais_lattice_lattice_type=_str_sp.bravais_lattice_type;
    str_in.bravais_lattice_lattice_variation_type=_str_sp.bravais_lattice_variation_type;
    str_in.bravais_lattice_lattice_system=_str_sp.bravais_lattice_system;
    str_sp.bravais_lattice_lattice_type=_str_sp.bravais_lattice_type;
    str_sp.bravais_lattice_lattice_variation_type=_str_sp.bravais_lattice_variation_type;
    str_sp.bravais_lattice_lattice_system=_str_sp.bravais_lattice_system;
    return TRUE;
  }
}
//DX - END

namespace LATTICE {
  bool Lattice(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,string& bravais_lattice_type,string& bravais_lattice_variation_type,string& bravais_lattice_system,double eps,double epsang) {
    xstructure str_in,str_sp,str_sc;
    str_in.lattice=lattice; str_in.scale=1.0;
    _atom atom; atom.cpos.clear();atom.fpos.clear();atom.type=0; str_in.AddAtom(atom);
    //  LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang); // STEFANO OLD VERSION
    int ss=0; //JUNKAI                                                                                                                                                                                                                                                
    LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang,ss,_EPS_); //JUNKAI
    lattice_sp=str_sp.lattice;lattice_sc=str_sc.lattice;
    bravais_lattice_type=str_sp.bravais_lattice_type;
    bravais_lattice_variation_type=str_sp.bravais_lattice_variation_type;
    bravais_lattice_system=str_sp.bravais_lattice_system;
    return TRUE;
  }
}

namespace LATTICE {
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,eps,epsang);
    return bravais_lattice_type;}
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,0.01,0.1);
    return bravais_lattice_type;}
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,double eps,double epsang) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    xmatrix<double> lattice_sp(3,3),lattice_sc(3,3);
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,eps,epsang);
    return bravais_lattice_type;}
  string Bravais_Lattice_Type(const xmatrix<double>& lattice) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    xmatrix<double> lattice_sp(3,3),lattice_sc(3,3);
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,0.01,0.1);
    return bravais_lattice_type;}
}

namespace LATTICE {
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,eps,epsang);
    return bravais_lattice_variation_type;}
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,0.01,0.1);
    return bravais_lattice_variation_type;}
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,double eps,double epsang) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    xmatrix<double> lattice_sp(3,3),lattice_sc(3,3);
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,eps,epsang);
    return bravais_lattice_variation_type;}
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    xmatrix<double> lattice_sp(3,3),lattice_sc(3,3);
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,0.01,0.1);
    return bravais_lattice_variation_type;}
}

namespace LATTICE {
  string Bravais_Lattice_System(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,eps,epsang);
    return bravais_lattice_system;}
  string Bravais_Lattice_System(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,0.01,0.1);
    return bravais_lattice_system;}
  string Bravais_Lattice_System(const xmatrix<double>& lattice,double eps,double epsang) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    xmatrix<double> lattice_sp(3,3),lattice_sc(3,3);
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,eps,epsang);
    return bravais_lattice_system;}
  string Bravais_Lattice_System(const xmatrix<double>& lattice) {
    string bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system;
    xmatrix<double> lattice_sp(3,3),lattice_sc(3,3);
    Lattice(lattice,lattice_sp,lattice_sc,bravais_lattice_type,bravais_lattice_variation_type,bravais_lattice_system,0.01,0.1);
    return bravais_lattice_system;}
}

// ***************************************************************************
// RECIPROCAL LATTICE FUNCTIONS

// ***************************************************************************
// Vrotate
// ***************************************************************************
xvector<double> Vrotate(xvector<double> V, xvector<double> Vaxis, double zzeta) {
  /*
    rotate x,y,z coordinate in V about Vaxis by zzeta(rad)
  */

  double theta,phi,sinphi,cosphi; //the spherical angles of vector Vaxis
  xvector<double> newV(3);
  xmatrix<double> MzMy(3,3),Mz(3,3),invMzMy(3,3),Mrot(3,3),tmpM(3,3);
  
  theta=acos(Vaxis(3));
  if(theta<1e-6) phi=0;
  else {
    cosphi=Vaxis(1)/sin(theta);
    sinphi=Vaxis(2)/sin(theta);
    if(cosphi>=0) phi=asin(sinphi);
    else {
      if(sinphi>=0) phi=acos(cosphi);
      else phi=-acos(cosphi);
    }
  }

  MzMy(1,1)=cos(phi)*cos(theta); MzMy(1,2)=-sin(phi); MzMy(1,3)=cos(phi)*sin(theta);
  MzMy(2,1)=sin(phi)*cos(theta); MzMy(2,2)= cos(phi); MzMy(2,3)=sin(phi)*sin(theta);
  MzMy(3,1)=-sin(theta);         MzMy(3,2)=0;         MzMy(3,3)=cos(theta);

  invMzMy=trasp(MzMy); //rotation matrix is unitary (transpose=inverse)

  Mz(1,1)=cos(zzeta); Mz(1,2)=-sin(zzeta); Mz(1,3)=0;
  Mz(2,1)=sin(zzeta); Mz(2,2)=cos(zzeta);  Mz(2,3)=0;
  Mz(3,1)=0;         Mz(3,2)=0;          Mz(3,3)=1;

  tmpM=MzMy*Mz;
  Mrot=tmpM*invMzMy;
  newV=Mrot*V;

  return newV;
}

// ***************************************************************************
// lattice::BZPLOTDATA
// ***************************************************************************
namespace LATTICE {
  void BZPLOTDATA(string options,istream& poscar, int mode) {
    /*
      for mode=0
      Synopsis:
      aflow --bzplotdata < POSCAR > plotbz.dat
      It reads POSCAR from stdin and generates data for brillouin zone and kpath plotting.
      The data is outputted to stdout (cout). bravais_lattice_variation_type of POSCAR will
      be calculated and used to get kpath

      for mode=10
      Synopsis:
      aflow --bzplotdatauseKPOINTS=KPOINTS < POSCAR > plotbz.dat
      It reads POSCAR from stdin and generates data for brillouin zone and kpath plotting.
      The data is outputted to stdout (cout). first word of the first line in KPOINTS file
      will be used as lattice_type to get kpath

      for mode=1
      Synopsis:
      aflow --bzplot < POSCAR
      It reads POSCAR from stdin and makes brillouin zone image file "bzplot.eps" and "bzplot.png".
      The data for plotting is saved in bzplot.dat and the script to make the plot is saved in plotbz.sh.
      plotbz.sh can be modified by hand to make further adjustment of the figure.

      for mode=11
      Synopsis:
      aflow --bzplotuseKPOINTS=KPOINTS < POSCAR
      analogous to --bzplotdatauseKPOINTS but for --bzplot.
    */

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "lattice::BZPLOTDATA: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(LDEBUG) cerr << "lattice::BZPLOTDATA: options=" << options << endl;
    if(LDEBUG) cerr << "lattice::BZPLOTDATA: tokens.size()=" << tokens.size() << endl;
    if(LDEBUG) cerr << "lattice::BZPLOTDATA: mode=" << mode << endl;
    
    if(mode==0 && tokens.size()!=0) {
      init::ErrorOption(cout,options,"lattice::BZPLOTDATA","aflow --bzplotdata < POSCAR > plotbz.dat");
      exit(0);
    }
    if(mode==10 && tokens.size()!=1) {
      init::ErrorOption(cout,options,"lattice::BZPLOTDATA","aflow --bzplotdatauseKPOINTS=KPOINTS < POSCAR > plotbz.dat");
      exit(0);
    }
    if(mode==1 && tokens.size()!=0) {
      init::ErrorOption(cout,options,"lattice::BZPLOTDATA","aflow --bzplot < POSCAR");
      exit(0);
    }
    if(mode==11 && tokens.size()!=1) {
      init::ErrorOption(cout,options,"lattice::BZPLOTDATA","aflow --bzplotuseKPOINTS=KPOINTS < POSCAR");
      exit(0);
    }

    string filename="";
    if(mode==10) filename=tokens.at(0);
    if(mode==11) filename=tokens.at(0);
   
#define bzsetf std::setw(8)<<std::fixed<<setprecision(4)

    bool found=false;
    int i,j;
    double a,b,c,alpha,beta,gamma,pLtheta;
    xstructure str_in(poscar,IOVASP_AUTO),str_sp;
    if(mode<10) {
      str_sp=GetStandardPrimitive(str_in);
    } else {
      str_sp=GetStandardPrimitive(str_in);
      mode=mode-10;
      vector<string> tokens;
      if(aurostd::FileExist(filename)) {
	aurostd::efile2vectorstring(filename,tokens);
      } else {
	cerr<<"ERROR: lattice::BZPLOTDATA: "<<filename<<" can not be opened. aborted.";abort();
      }
      if(tokens.size()>0) str_sp.bravais_lattice_variation_type=tokens.at(0);
      aurostd::string2tokens(string(tokens.at(0)),tokens," ");
      if(tokens.size()>0) str_sp.bravais_lattice_variation_type=tokens.at(0);
      if(tokens.size()>1) {
	if(tokens.at(0)=="KPOINTS:") 
	  str_sp.bravais_lattice_variation_type=tokens.at(1);
      }
    }
    cerr << "LATTICE TYPE="<<str_sp.bravais_lattice_variation_type<<endl;
    xvector<double> data(6),a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),pL(3),
      pLaxis(3),pz(3),ptmp(3),tmpb1(3),tmpb2(3),tmpb3(3);
    xmatrix<double> klattice(3,3);
    string LattVar="";

    LattVar=str_sp.bravais_lattice_variation_type;
    data=LATTICE::Getabc_angles_Conventional(str_sp.lattice,LattVar,0);//angles in radiant
    a=data(1); b=data(2); c=data(3);  alpha=data(4); beta=data(5); gamma=data(6);
    if(beta) {;} // dummy load
    if(gamma) {;} // dummy load

    klattice=ReciprocalLattice(str_sp.lattice);
    for(i=1;i<4;i++) {
      a1(i)=str_sp.lattice(1,i);
      a2(i)=str_sp.lattice(2,i);
      a3(i)=str_sp.lattice(3,i);
      b1(i)=klattice(1,i);
      b2(i)=klattice(2,i);
      b3(i)=klattice(3,i);
    }

    pL(1)=0;pL(2)=0;pL(3)=1.0;
    int Nbz=0,Nirbz=0,Ncon=0,Nconback=0,Nircon=0,Nirconback=0;
    double eta=0.0,zzeta=0.0,nu=0.0,delta=0.0,lambda=0.0,mu=0.0,omega=0.0,psi=0.0,phi=0.0,rho=0.0,tau=0.0,theta=0.0;
    double xrotview=10,zrotview=115;//rotation view for gnuplot's set view command in degrees
    double elview=10,azview=120;//elevation and azimuth view for matlab in degrees
    double b1arrow=0.25,b2arrow=0.25,b3arrow=0.25;//lenght of recip. vec. arrow in percent
    double b1text=0.25,b2text=0.25,b3text=0.25;//b1 b2 b3 label of recip. vec. arrow in percent
    xmatrix<double> bz(50,3),irbz(50,3);
    xmatrix<int> con(2,100),conback(2,100),ircon(2,50),irconback(2,50);
    vector<string> irbzlab(50+1);
    string title="",glatt="",gpath="";

    //----------- CUB  --------------------------
    if(LattVar=="CUB" || LattVar=="cub") {
      found=true;
      glatt="CUB";
      xrotview=80; zrotview=115;
      azview=120; elview=10;
      b1arrow=0.6; b2arrow=0.4; b3arrow=0.3;
      b1text=0.75; b2text=0.42; b3text=0.33;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="M"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="R"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="X"; Nirbz=i;
      gpath="{/Symbol G}-X-M-{/Symbol G}-R-X|M-R";
      i=1;
      ircon(1,i)=1; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=4;  i++;
      ircon(1,i)=2; ircon(2,i)=3;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0; Nbz=i;
      //connection index
      i=1;
      con(1,i)=2; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=8;  i++;
      con(1,i)=8; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=2;  i++;
      con(1,i)=10; con(2,i)=3;  i++;
      con(1,i)=3; con(2,i)=7;  i++;
      con(1,i)=5; con(2,i)=3;  i++;
      con(1,i)=7; con(2,i)=6;  i++;
      con(1,i)=5; con(2,i)=4;  Ncon=i;
      i=1;
      conback(1,i)=9; conback(2,i)=5;  i++;
      conback(1,i)=9; conback(2,i)=8;  i++;
      conback(1,i)=7; conback(2,i)=9;  Nconback=i;
    }
    //----------- FCC  --------------------------
    if(LattVar=="FCC" || LattVar=="fcc") {
      found=true;
      glatt="FCC";
      xrotview=80; zrotview=115;
      azview=120; elview=10;
      b1arrow=0.2; b2arrow=0.2; b3arrow=0.2;
      b1text=0.22; b2text=0.28; b3text=0.23;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=3.0/8; irbz(i,2)=3.0/8; irbz(i,3)=0.75;  irbzlab[i]="K"; i++;
      irbz(i,1)=0.5;   irbz(i,2)=0.5;   irbz(i,3)=0.5;   irbzlab[i]="L"; i++;
      irbz(i,1)=5.0/8; irbz(i,2)=0.25;  irbz(i,3)=5.0/8; irbzlab[i]="U"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.25;  irbz(i,3)=0.75;  irbzlab[i]="W"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0;     irbz(i,3)=0.5;   irbzlab[i]="X"; Nirbz=i;
      gpath="{/Symbol G}-X-W-K-{/Symbol G}-L-U-W-L-K|U-X";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=2;  i++;
      ircon(1,i)=4; ircon(2,i)=6;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=3/8.0; bz(i,2)=3/8.0; bz(i,3)=0.75; i++;
      bz(i,1)=0.5; bz(i,2)=0.25; bz(i,3)=0.75; i++;
      bz(i,1)=5/8.0; bz(i,2)=0.25; bz(i,3)=5/8.0; i++;
      //right 4
      bz(i,1)=0.75; bz(i,2)=0.25; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-0.25; bz(i,3)=0.25; i++;
      bz(i,1)=0.25; bz(i,2)=-0.25; bz(i,3)=0.5; i++;
      //front 4
      bz(i,1)=0.25; bz(i,2)=0.5; bz(i,3)=0.75; i++;
      bz(i,1)=0.25; bz(i,2)=0.75; bz(i,3)=0.5; i++;
      bz(i,1)=-0.25; bz(i,2)=0.5; bz(i,3)=0.25; i++;
      bz(i,1)=-0.25; bz(i,2)=0.25; bz(i,3)=0.5; i++;
      //top 4
      bz(i,1)=0.5; bz(i,2)=0.75; bz(i,3)=0.25; i++;
      bz(i,1)=0.75; bz(i,2)=0.5; bz(i,3)=0.25; i++;
      bz(i,1)=0.5; bz(i,2)=0.25; bz(i,3)=-0.25; i++;
      bz(i,1)=0.25; bz(i,2)=0.5; bz(i,3)=-0.25; i++;
      //bottom 4
      bz(i,1)=-0.5; bz(i,2)=-0.25; bz(i,3)=0.25; i++;
      bz(i,1)=-0.25; bz(i,2)=-0.5; bz(i,3)=0.25; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.75; bz(i,3)=-0.25; i++;
      bz(i,1)=-0.75; bz(i,2)=-0.5; bz(i,3)=-0.25; i++;
      //left 4
      bz(i,1)=-0.75; bz(i,2)=-0.25; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.25; bz(i,3)=-0.75; i++;
      bz(i,1)=-0.25; bz(i,2)=0.25; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.25; bz(i,3)=-0.25; i++;
      //back 4
      bz(i,1)=-0.25; bz(i,2)=-0.5; bz(i,3)=-0.75; i++;
      bz(i,1)=-0.25; bz(i,2)=-0.75; bz(i,3)=-0.5; i++;
      bz(i,1)=0.25; bz(i,2)=-0.5; bz(i,3)=-0.25; i++;
      bz(i,1)=0.25; bz(i,2)=-0.25; bz(i,3)=-0.5; Nbz=i;
      //connection index
      i=1;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=7; con(2,i)=3;  i++;
      con(1,i)=2; con(2,i)=8;  i++;
      con(1,i)=8; con(2,i)=9;  i++;
      con(1,i)=9; con(2,i)=10;  i++;
      con(1,i)=7; con(2,i)=17;  i++;
      con(1,i)=8; con(2,i)=11;  i++;
      con(1,i)=9; con(2,i)=12;  i++;
      con(1,i)=5; con(2,i)=13;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=13; con(2,i)=14;  i++;
      con(1,i)=14; con(2,i)=15;  i++;
      con(1,i)=15; con(2,i)=12;  i++;
      con(1,i)=11; con(2,i)=16;  i++;
      con(1,i)=16; con(2,i)=17;  i++;
      con(1,i)=19; con(2,i)=16;  i++;
      con(1,i)=10; con(2,i)=23;  i++;
      con(1,i)=23; con(2,i)=20;  i++;
      con(1,i)=22; con(2,i)=23;  i++;
      con(1,i)=22; con(2,i)=15;  i++;
      con(1,i)=19; con(2,i)=20;  Ncon=i;
      i=1;
      conback(1,i)=6; conback(2,i)=26;  i++;
      conback(1,i)=17; conback(2,i)=18;  i++;
      conback(1,i)=18; conback(2,i)=25;  i++;
      conback(1,i)=18; conback(2,i)=19;  i++;
      conback(1,i)=14; conback(2,i)=27;  i++;
      conback(1,i)=20; conback(2,i)=21;  i++;
      conback(1,i)=21; conback(2,i)=22;  i++;
      conback(1,i)=21; conback(2,i)=24;  i++;
      conback(1,i)=24; conback(2,i)=25;  i++;
      conback(1,i)=25; conback(2,i)=26;  i++;
      conback(1,i)=26; conback(2,i)=27;  i++;
      conback(1,i)=27; conback(2,i)=24;  Nconback=i;
    }
    //----------- BCC  --------------------------
    if(LattVar=="BCC" || LattVar=="bcc") {
      found=true;
      glatt="BCC";
      xrotview=84; zrotview=118;
      azview=115; elview=5;
      b1arrow=0.4; b2arrow=0.5; b3arrow=1.0;
      b1text=0.42; b2text=0.56; b3text=1.05;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0.5; irbzlab[i]="H"; i++;
      irbz(i,1)=0.25; irbz(i,2)=0.25; irbz(i,3)=0.25; irbzlab[i]="P"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="N"; Nirbz=i;
      gpath="{/Symbol G}-H-N-{/Symbol G}-P-H|P-N";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=2;  i++;
      ircon(1,i)=3; ircon(2,i)=4;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //right 4
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.25; bz(i,2)=0.25; bz(i,3)=0.25; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.25; bz(i,2)=-0.25; bz(i,3)=0.75; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      //front 4
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.75; bz(i,2)=0.25; bz(i,3)=0.25; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.25; bz(i,2)=0.75; bz(i,3)=-0.25; i++;
      //back 4
      bz(i,1)=0.75; bz(i,2)=-0.25; bz(i,3)=-0.25; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=0.25; bz(i,2)=-0.75; bz(i,3)=0.25; i++;
      //left 4
      bz(i,1)=-0.25; bz(i,2)=-0.25; bz(i,3)=-0.25; i++;
      bz(i,1)=0.25; bz(i,2)=0.25; bz(i,3)=-0.75; Nbz=i;
      //connection index
      i=1;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=5; con(2,i)=2;  i++;
      con(1,i)=4; con(2,i)=8;  i++;
      con(1,i)=8; con(2,i)=9;  i++;
      con(1,i)=7; con(2,i)=5;  i++;
      con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=3; con(2,i)=6;  i++;
      con(1,i)=9; con(2,i)=10;  i++;
      con(1,i)=4; con(2,i)=10;  i++;
      con(1,i)=2; con(2,i)=11;  i++;
      con(1,i)=2; con(2,i)=13;  i++;
      con(1,i)=7; con(2,i)=13;  i++;
      con(1,i)=6; con(2,i)=10;  i++;
      con(1,i)=6; con(2,i)=11;  Ncon=i;
      i=1;
      conback(1,i)=7; conback(2,i)=14;  i++;
      conback(1,i)=9; conback(2,i)=14;  i++;
      conback(1,i)=9; conback(2,i)=15;  i++;
      conback(1,i)=6; conback(2,i)=15;  i++;
      conback(1,i)=11; conback(2,i)=12;  i++;
      conback(1,i)=13; conback(2,i)=12;  i++;
      conback(1,i)=12; conback(2,i)=14;  i++;
      conback(1,i)=12; conback(2,i)=15;  Nconback=i;
    }
    //----------- TET  --------------------------
    if(LattVar=="TET" || LattVar=="tet") {
      found=true;
      glatt="TET";
      xrotview=80; zrotview=115;
      azview=110; elview=7;
      //b1arrow=0.6; b2arrow=0.35; b3arrow=0.4;
      //b1text=0.75; b2text=0.36; b3text=0.44;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="A"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="M"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="R"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-X-M-{/Symbol G}-Z-R-A-Z|X-R|M-A";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=6;  i++;
      ircon(1,i)=5; ircon(2,i)=4;  i++;
      ircon(1,i)=3; ircon(2,i)=2;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0.5; Nbz=i;
      //connection index
      i=1;
      con(1,i)=3; con(2,i)=7;  i++;
      con(1,i)=5; con(2,i)=3;  i++;
      con(1,i)=7; con(2,i)=6;  i++;
      con(1,i)=5; con(2,i)=4;  i++;
      con(1,i)=6; con(2,i)=8;  i++;
      con(1,i)=8; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=2;  i++;
      con(1,i)=3; con(2,i)=10;  i++;
      con(1,i)=6; con(2,i)=11;  Ncon=i;
      i=1;
      conback(1,i)=9; conback(2,i)=5;  i++;
      conback(1,i)=9; conback(2,i)=7;  i++;
      conback(1,i)=9; conback(2,i)=8;  Nconback=i;
    }
    //----------- BCT1 --------------------------
    if(LattVar=="BCT1" || LattVar=="bct1") {
      found=true;
      glatt="BCT_1";
      xrotview=72; zrotview=100;
      azview=105; elview=15;
      //b1arrow=0.25; b2arrow=0.4; b3arrow=0.5;
      //b1text=0.26; b2text=0.42; b3text=0.52;
      eta=0.25*(1.0+c*c/a/a);
      i=1;
      irbz(i,1)=0;    irbz(i,2)=0;       irbz(i,3)=0;    irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=-0.5; irbz(i,2)=0.5;     irbz(i,3)=0.5;  irbzlab[i]="M"; i++;
      irbz(i,1)=0;    irbz(i,2)=0.5;     irbz(i,3)=0;    irbzlab[i]="N"; i++;
      irbz(i,1)=0.25; irbz(i,2)=0.25;    irbz(i,3)=0.25; irbzlab[i]="P"; i++;
      irbz(i,1)=0;    irbz(i,2)=0;       irbz(i,3)=0.5;  irbzlab[i]="X"; i++;
      irbz(i,1)=eta;  irbz(i,2)=eta;     irbz(i,3)=-eta; irbzlab[i]="Z"; i++;
      irbz(i,1)=-eta; irbz(i,2)=1.0-eta; irbz(i,3)=eta;  irbzlab[i]="Z_1"; Nirbz=i;
      gpath="{/Symbol G}-X-M-{/Symbol G}-Z-P-N-Z_1-M|X-P";
      //    1-5-2-1-6-4-3-7-2  5-4
      i=1;
      ircon(1,i)=1; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=7;  i++;
      ircon(1,i)=7; ircon(2,i)=2;  i++;
      ircon(1,i)=5; ircon(2,i)=4;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=3; irconback(2,i)=6;  i++;
      irconback(1,i)=4; irconback(2,i)=7;  Nirconback=i;
      //bz
      //front 6
      i=1;
      bz(i,1)=0;       bz(i,2)=0;       bz(i,3)=0; i++;
      bz(i,1)=-0.25;   bz(i,2)=-0.25;   bz(i,3)=0.75; i++;
      bz(i,1)=eta;     bz(i,2)=eta-1.0; bz(i,3)=1.0-eta; i++;
      bz(i,1)=1.0-eta; bz(i,2)=-eta;    bz(i,3)=eta; i++;
      bz(i,1)=0.25;    bz(i,2)=0.25;    bz(i,3)=0.25; i++;
      bz(i,1)=-eta;    bz(i,2)=1.0-eta; bz(i,3)=eta; i++;
      bz(i,1)=-0.5;    bz(i,2)=0.5;     bz(i,3)=0.5; i++;
      bz(i,1)=eta-1.0; bz(i,2)=eta;     bz(i,3)=1.0-eta; Nbz=i;
      //back 6
      for(j=1;j<4;j++) {
	bz(Nbz+1,j)=-bz(2,j);
	bz(Nbz+2,j)=-bz(3,j);
	bz(Nbz+3,j)=-bz(4,j);
	bz(Nbz+4,j)=-bz(5,j);
	bz(Nbz+5,j)=-bz(6,j);
	bz(Nbz+6,j)=-bz(8,j);
      }
      //middle 6
      i=Nbz+7;
      bz(15,1)=-0.25; bz(15,2)=0.75;  bz(15,3)=-0.25; i++;
      bz(16,1)=eta;   bz(16,2)=eta;   bz(16,3)=-eta; i++;
      bz(17,1)=0.75;  bz(17,2)=-0.25; bz(17,3)=-0.25; Nbz=i;
      for(j=1;j<4;j++) {
	bz(Nbz+1,j)=-bz(15,j);
	bz(Nbz+2,j)=-bz(16,j);
	bz(Nbz+3,j)=-bz(17,j);
      }
      Nbz=Nbz+3;
      //connection index
      i=1;
      con(1,i)=2;   con(2,i)=3;  i++;
      con(1,i)=3;   con(2,i)=4;  i++;
      con(1,i)=4;   con(2,i)=5;  i++;
      con(1,i)=7;   con(2,i)=8;  i++;
      con(1,i)=8;   con(2,i)=2;  i++;
      con(1,i)=9;   con(2,i)=10;  i++;
      con(1,i)=10;  con(2,i)=11;  i++;
      con(1,i)=6;  con(2,i)=15;  i++;
      con(1,i)=15; con(2,i)=10;  i++;
      con(1,i)=11; con(2,i)=20;  i++;
      con(1,i)=20; con(2,i)=8;  i++;
      con(1,i)=17; con(2,i)=4;  i++;
      con(1,i)=16; con(2,i)=9;  i++;
      con(1,i)=15; con(2,i)=16;  i++;
      con(1,i)=16; con(2,i)=17;  i++;
      con(1,i)=19; con(2,i)=20;  i++;
      con(1,i)=2;  con(2,i)=19;  Ncon=i;
      i=1;
      conback(1,i)=11; conback(2,i)=12;  i++;
      conback(1,i)=12; conback(2,i)=13;  i++;
      conback(1,i)=13; conback(2,i)=14;  i++;
      conback(1,i)=14; conback(2,i)=9;  i++;
      conback(1,i)=3;  conback(2,i)=18;  i++;
      conback(1,i)=18; conback(2,i)=13;  i++;
      conback(1,i)=14; conback(2,i)=17;  i++;
      conback(1,i)=18; conback(2,i)=19;  i++;
      conback(1,i)=19; conback(2,i)=12;  Nconback=i;
    }
    //----------- BCT2  --------------------------
    if(LattVar=="BCT2" || LattVar=="bct2") {
      found=true;
      glatt="BCT_2";
      xrotview=80; zrotview=115;
      azview=112; elview=12;
      //b1arrow=0.32; b2arrow=0.5; b3arrow=0.55;
      //b1text=0.33; b2text=0.65; b3text=0.57;
      eta=0.25*(1+a*a/c/c); zzeta=0.5*a*a/c/c;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="N"; i++;
      irbz(i,1)=0.25; irbz(i,2)=0.25; irbz(i,3)=0.25; irbzlab[i]="P"; i++;
      irbz(i,1)=-eta; irbz(i,2)=eta; irbz(i,3)=eta; irbzlab[i]="{/Symbol S}"; i++;
      irbz(i,1)=eta; irbz(i,2)=1-eta; irbz(i,3)=-eta; irbzlab[i]="{/Symbol S}_1"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="X"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=zzeta; irbz(i,3)=0.5; irbzlab[i]="Y"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=-zzeta; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=-0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-X-Y-{/Symbol S}-{/Symbol G}-Z-{/Symbol S}_1-N-P-Y_1-Z|X-P";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=7;  i++;
      ircon(1,i)=7; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=9;  i++;
      ircon(1,i)=9; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=8;  i++;
      ircon(1,i)=8; ircon(2,i)=9;  i++;
      ircon(1,i)=6; ircon(2,i)=3;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=2; irconback(2,i)=4;  i++;
      irconback(1,i)=3; irconback(2,i)=7;  i++;
      irconback(1,i)=5; irconback(2,i)=8;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //front 4
      bz(i,1)=-0.25; bz(i,2)=-0.25; bz(i,3)=0.75; i++;
      bz(i,1)=zzeta; bz(i,2)=-zzeta; bz(i,3)=0.5; i++;
      bz(i,1)=0.25; bz(i,2)=0.25; bz(i,3)=0.25; i++;
      bz(i,1)=-zzeta; bz(i,2)=zzeta; bz(i,3)=0.5; i++;
      //back 4
      bz(i,1)=0.25; bz(i,2)=0.25; bz(i,3)=-0.75; i++;
      bz(i,1)=-zzeta; bz(i,2)=zzeta; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.25; bz(i,2)=-0.25; bz(i,3)=-0.25; i++;
      bz(i,1)=zzeta; bz(i,2)=-zzeta; bz(i,3)=-0.5; i++;
      //left 4
      bz(i,1)=-0.75; bz(i,2)=0.25; bz(i,3)=0.25; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=zzeta; i++;
      bz(i,1)=-0.25; bz(i,2)=0.75; bz(i,3)=-0.25; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=-zzeta; i++;
      //right 4
      bz(i,1)=0.75; bz(i,2)=-0.25; bz(i,3)=-0.25; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=-zzeta; i++;
      bz(i,1)=0.25; bz(i,2)=-0.75; bz(i,3)=0.25; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=zzeta; i++;
      //top 5
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=-zzeta; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=zzeta; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=zzeta-1.0; i++;
      bz(i,1)=zzeta; bz(i,2)=1.0-zzeta; bz(i,3)=-0.5; i++;
      bz(i,1)=eta; bz(i,2)=1.0-eta; bz(i,3)=-eta; i++;
      //bottom
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=zzeta; i++;
      bz(i,1)=-1.0+zzeta; bz(i,2)=-zzeta; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=-zzeta+1.0; i++;
      bz(i,1)=-zzeta; bz(i,2)=-1.0+zzeta; bz(i,3)=0.5; i++;
      bz(i,1)=-eta; bz(i,2)=eta; bz(i,3)=eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=2; con(2,i)=3;  i++;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=2; con(2,i)=5;  i++;
      con(1,i)=2; con(2,i)=25;  i++;
      con(1,i)=3; con(2,i)=17;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      con(1,i)=13; con(2,i)=10;  i++;
      con(1,i)=16; con(2,i)=17;  i++;
      con(1,i)=17; con(2,i)=14;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=19; con(2,i)=20;  i++;
      con(1,i)=10; con(2,i)=24;  i++;
      con(1,i)=16; con(2,i)=26;  i++;
      con(1,i)=11; con(2,i)=27;  i++;
      con(1,i)=14; con(2,i)=19;  i++;
      con(1,i)=12; con(2,i)=21;  i++;
      con(1,i)=20; con(2,i)=21;  i++;
      con(1,i)=21; con(2,i)=22;  i++;
      con(1,i)=24; con(2,i)=25;  i++;
      con(1,i)=25; con(2,i)=26;  Ncon=i;
      i=1;
      conback(1,i)=6; conback(2,i)=7;  i++;
      conback(1,i)=7; conback(2,i)=8;  i++;
      conback(1,i)=8; conback(2,i)=9;  i++;
      conback(1,i)=9; conback(2,i)=6;  i++;
      conback(1,i)=14; conback(2,i)=15;  i++;
      conback(1,i)=15; conback(2,i)=16;  i++;
      conback(1,i)=23; conback(2,i)=24;  i++;
      conback(1,i)=26; conback(2,i)=23;  i++;
      conback(1,i)=23; conback(2,i)=8;  i++;
      conback(1,i)=15; conback(2,i)=9;  i++;
      conback(1,i)=13; conback(2,i)=7;  i++;
      conback(1,i)=20; conback(2,i)=6;  Nconback=i;
    }
    //----------- ORC --------------------------
    if(LattVar=="ORC" || LattVar=="orc") {
      found=true;
      glatt="ORC";
      xrotview=80; zrotview=115;
      azview=115; elview=10;
      //b1arrow=0.4; b2arrow=0.4; b3arrow=0.4;
      //b1text=0.55; b2text=0.42; b3text=0.43;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="R"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="S"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="T"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="U"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Y"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-X-S-Y-{/Symbol G}-Z-U-R-T-Z|Y-T|U-X|S-R";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=7;  i++;
      ircon(1,i)=7; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=8;  i++;
      ircon(1,i)=8; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=8;  i++;
      ircon(1,i)=7; ircon(2,i)=4;  i++;
      ircon(1,i)=5; ircon(2,i)=6;  i++;
      ircon(1,i)=3; ircon(2,i)=2;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0.5; Nbz=i;
      //connection index
      i=1;
      con(1,i)=12; con(2,i)=6;  i++;
      con(1,i)=11; con(2,i)=4;  i++;
      con(1,i)=10; con(2,i)=3;  i++;
      con(1,i)=6; con(2,i)=8;  i++;
      con(1,i)=8; con(2,i)=4;  i++;
      con(1,i)=3; con(2,i)=7;  i++;
      con(1,i)=5; con(2,i)=3;  i++;
      con(1,i)=7; con(2,i)=6;  i++;
      con(1,i)=5; con(2,i)=4;  Ncon=i;
      i=1;
      conback(1,i)=7; conback(2,i)=9;  i++;
      conback(1,i)=8; conback(2,i)=9;  i++;
      conback(1,i)=5; conback(2,i)=9;  Nconback=i;
    }
    //----------- ORCC --------------------------
    if(LattVar=="ORCC" || LattVar=="orcc") {
      found=true;
      glatt="ORCC";
      xrotview=80; zrotview=120;
      azview=110; elview=10;
      //b1arrow=0.3; b2arrow=1; b3arrow=0.6;
      //b1text=0.4; b2text=1.05; b3text=0.65;
      zzeta=0.25*(1+a*a/b/b);
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=zzeta; irbz(i,3)=0.5; irbzlab[i]="A"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=1.0-zzeta; irbz(i,3)=0.5; irbzlab[i]="A_1"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="R"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="S"; i++;
      irbz(i,1)=-0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="T"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=zzeta; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=1.0-zzeta; irbz(i,3)=0; irbzlab[i]="X_1"; i++;
      irbz(i,1)=-0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Y"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-X-S-R-A-Z-{/Symbol G}-Y-X_1-A_1-T-Y|Z-T";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=7;  i++;
      ircon(1,i)=7; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=2;  i++;
      ircon(1,i)=1; ircon(2,i)=9;  i++;
      ircon(1,i)=9; ircon(2,i)=8;  i++;
      ircon(1,i)=8; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=9;  i++;
      ircon(1,i)=10; ircon(2,i)=2;  i++;
      ircon(1,i)=10; ircon(2,i)=6;  i++;
      ircon(1,i)=10; ircon(2,i)=1;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      i=1;
      irconback(1,i)=2; irconback(2,i)=7;  i++;
      irconback(1,i)=3; irconback(2,i)=4;  i++;
      irconback(1,i)=8; irconback(2,i)=5;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //top
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=zzeta-1.0; bz(i,2)=zzeta; bz(i,3)=0.5; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta-1.0; bz(i,3)=0.5; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=-zzeta; bz(i,3)=0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=0.5; i++;
      //bottom
      bz(i,1)=-zzeta; bz(i,2)=1.0-zzeta; bz(i,3)=-0.5; i++;
      bz(i,1)=zzeta-1.0; bz(i,2)=zzeta; bz(i,3)=-0.5; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=-0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta-1.0; bz(i,3)=-0.5; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=-zzeta; bz(i,3)=-0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=-0.5; i++;
      //middle
      bz(i,1)=-zzeta; bz(i,2)=1.0-zzeta; bz(i,3)=0; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=0; Nbz=i;
      //connection index
      i=1;
      con(1,i)=2; con(2,i)=3;  i++;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=8; con(2,i)=9;  i++;
      con(1,i)=9; con(2,i)=3;  i++;
      con(1,i)=8; con(2,i)=13;  i++;
      con(1,i)=8; con(2,i)=14;  i++;
      con(1,i)=6; con(2,i)=12;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      con(1,i)=13; con(2,i)=15;  Ncon=i;
      i=1;
      conback(1,i)=9; conback(2,i)=10;  i++;
      conback(1,i)=4; conback(2,i)=10;  i++;
      conback(1,i)=5; conback(2,i)=11;  i++;
      conback(1,i)=10; conback(2,i)=11;  i++;
      conback(1,i)=11; conback(2,i)=12;  Nconback=i;
    }
    //----------- ORCI --------------------------
    if(LattVar=="ORCI" || LattVar=="orci") {
      found=true;
      glatt="ORCI";
      xrotview=80; zrotview=115;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=.25; b3text=0.25;
      zzeta=0.25*(1+a*a/c/c); eta=0.25*(1+b*b/c/c);
      delta=0.25*(b*b-a*a)/c/c; mu=0.25*(a*a+b*b)/c/c;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=-mu; irbz(i,2)=mu; irbz(i,3)=0.5-delta; irbzlab[i]="L"; i++;
      irbz(i,1)=mu; irbz(i,2)=-mu; irbz(i,3)=0.5+delta; irbzlab[i]="L_1"; i++;
      irbz(i,1)=0.5-delta; irbz(i,2)=0.5+delta; irbz(i,3)=-mu; irbzlab[i]="L_2"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="R"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="S"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="T"; i++;
      irbz(i,1)=0.25; irbz(i,2)=0.25; irbz(i,3)=0.25; irbzlab[i]="W"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=zzeta; irbz(i,3)=zzeta; irbzlab[i]="X"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=1.0-zzeta; irbz(i,3)=-zzeta; irbzlab[i]="X_1"; i++;
      irbz(i,1)=eta; irbz(i,2)=-eta; irbz(i,3)=eta; irbzlab[i]="Y"; i++;
      irbz(i,1)=1.0-eta; irbz(i,2)=eta; irbz(i,3)=-eta; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=-0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-X-L-T-W-R-X_1-Z-{/Symbol G}-Y-S-W|L_1-Y|Y_1-Z";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=9;  i++;
      ircon(1,i)=9; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=7;  i++;
      ircon(1,i)=7; ircon(2,i)=8;  i++;
      ircon(1,i)=8; ircon(2,i)=5;  i++;
      ircon(1,i)=6; ircon(2,i)=8;  i++;
      ircon(1,i)=5; ircon(2,i)=10;  i++;
      ircon(1,i)=1; ircon(2,i)=13;  i++;
      ircon(1,i)=1; ircon(2,i)=11;  i++;
      ircon(1,i)=6; ircon(2,i)=11;  i++;
      ircon(1,i)=3; ircon(2,i)=11;  i++;
      ircon(1,i)=10; ircon(2,i)=13;  i++;
      ircon(1,i)=12; ircon(2,i)=13;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=5; irconback(2,i)=9;  i++;
      irconback(1,i)=10; irconback(2,i)=4;  i++;
      irconback(1,i)=12; irconback(2,i)=4;  i++;
      irconback(1,i)=2; irconback(2,i)=8;  i++;
      irconback(1,i)=4; irconback(2,i)=8;  i++;
      irconback(1,i)=7; irconback(2,i)=3;  i++;
      irconback(1,i)=8; irconback(2,i)=3;  i++;
      irconback(1,i)=6; irconback(2,i)=12;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //top
      bz(i,1)=1.0-eta; bz(i,2)=eta; bz(i,3)=-eta; i++;
      bz(i,1)=1.0-mu; bz(i,2)=mu; bz(i,3)=-delta-0.5; i++;
      bz(i,1)=0.5+delta; bz(i,2)=0.5-delta; bz(i,3)=mu-1.0; i++;
      bz(i,1)=mu; bz(i,2)=1.0-mu; bz(i,3)=delta-0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=1.0-zzeta; bz(i,3)=-zzeta; i++;
      //bot
      bz(i,1)=delta-0.5; bz(i,2)=-delta-0.5; bz(i,3)=mu; i++;
      bz(i,1)=-1.0+mu; bz(i,2)=-mu; bz(i,3)=delta+0.5; i++;
      bz(i,1)=-0.5-delta; bz(i,2)=-0.5+delta; bz(i,3)=-mu+1.0; i++;
      bz(i,1)=-mu; bz(i,2)=-1.0+mu; bz(i,3)=-delta+0.5; i++;
      //front
      bz(i,1)=-0.25; bz(i,2)=-0.25; bz(i,3)=0.75; i++;
      bz(i,1)=mu; bz(i,2)=-mu; bz(i,3)=0.5+delta; i++;
      bz(i,1)=0.25; bz(i,2)=0.25; bz(i,3)=0.25; i++;
      bz(i,1)=-mu; bz(i,2)=mu; bz(i,3)=0.5-delta; i++;
      //J
      bz(i,1)=0.5-delta; bz(i,2)=0.5+delta; bz(i,3)=-mu; i++;
      //back
      bz(i,1)=0.25; bz(i,2)=0.25; bz(i,3)=-0.75; i++;
      bz(i,1)=-mu; bz(i,2)=mu; bz(i,3)=-0.5-delta; i++;
      bz(i,1)=-0.25; bz(i,2)=-0.25; bz(i,3)=-0.25; i++;
      bz(i,1)=mu; bz(i,2)=-mu; bz(i,3)=-0.5+delta; i++;
      //left
      bz(i,1)=-0.75; bz(i,2)=0.25; bz(i,3)=0.25; i++;
      bz(i,1)=delta-0.5; bz(i,2)=0.5-delta; bz(i,3)=mu; i++;
      bz(i,1)=-0.25; bz(i,2)=0.75; bz(i,3)=-0.25; i++;
      bz(i,1)=-0.5-delta; bz(i,2)=0.5+delta; bz(i,3)=-mu; i++;
      //sigma
      bz(i,1)=-zzeta; bz(i,2)=zzeta; bz(i,3)=zzeta; i++;
      //right
      bz(i,1)=0.75; bz(i,2)=-0.25; bz(i,3)=-0.25; i++;
      bz(i,1)=-delta+0.5; bz(i,2)=-0.5+delta; bz(i,3)=-mu; i++;
      bz(i,1)=0.25; bz(i,2)=-0.75; bz(i,3)=0.25; i++;
      bz(i,1)=0.5+delta; bz(i,2)=-0.5-delta; bz(i,3)=mu; i++;
      //Y
      bz(i,1)=eta; bz(i,2)=-eta; bz(i,3)=eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=2; con(2,i)=3;  i++;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=8; con(2,i)=9;  i++;
      con(1,i)=9; con(2,i)=11;  i++;
      con(1,i)=9; con(2,i)=10;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      //    con(1,i)=12; con(2,i)=13;  i++;
      //con(1,i)=13; con(2,i)=14;  i++;
      con(1,i)=14; con(2,i)=11;  i++;
      //con(1,i)=13; con(2,i)=15;  i++;
      con(1,i)=20; con(2,i)=21;  i++;
      con(1,i)=21; con(2,i)=22;  i++;
      con(1,i)=22; con(2,i)=23;  i++;
      con(1,i)=23; con(2,i)=20;  i++;
      con(1,i)=21; con(2,i)=24;  i++;
      con(1,i)=27; con(2,i)=28;  i++;
      con(1,i)=28; con(2,i)=25;  i++;
      con(1,i)=28; con(2,i)=29;  i++;
      con(1,i)=5; con(2,i)=22;  i++;
      con(1,i)=8; con(2,i)=20;  i++;
      con(1,i)=3; con(2,i)=25;  i++;
      con(1,i)=10; con(2,i)=27;  Ncon=i;
      i=1;
      conback(1,i)=8; conback(2,i)=7;  i++;
      conback(1,i)=7; conback(2,i)=10;  i++;
      conback(1,i)=7; conback(2,i)=18;  i++;
      conback(1,i)=4; conback(2,i)=16;  i++;
      conback(1,i)=16; conback(2,i)=17;  i++;
      conback(1,i)=17; conback(2,i)=18;  i++;
      conback(1,i)=18; conback(2,i)=19;  i++;
      conback(1,i)=19; conback(2,i)=16;  i++;
      conback(1,i)=17; conback(2,i)=23;  i++;
      conback(1,i)=25; conback(2,i)=26;  i++;
      conback(1,i)=26; conback(2,i)=27;  i++;
      conback(1,i)=19; conback(2,i)=26;  Nconback=i;
    }
    //----------- ORCF1 --------------------------
    if(LattVar=="ORCF1" || LattVar=="orcf1") {
      found=true;
      glatt="ORCF_1";
      xrotview=82; zrotview=132;
      azview=125; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      zzeta=0.25*(1.0+a*a/b/b-a*a/c/c); eta=0.25*(1.0+a*a/b/b+a*a/c/c);
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5+zzeta; irbz(i,3)=zzeta; irbzlab[i]="A"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5-zzeta; irbz(i,3)=1.0-zzeta; irbzlab[i]="A_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="L"; i++;
      irbz(i,1)=1; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="T"; i++;
      irbz(i,1)=0; irbz(i,2)=eta; irbz(i,3)=eta; irbzlab[i]="X"; i++;
      irbz(i,1)=1; irbz(i,2)=1.0-eta; irbz(i,3)=1.0-eta; irbzlab[i]="X_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Y"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-T-Z-{/Symbol G}-X-A_1-Y|T-X_1|X-A-Z|L-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=8;  i++;
      ircon(1,i)=8; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=9;  i++;
      ircon(1,i)=9; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=8;  i++;
      ircon(1,i)=5; ircon(2,i)=7;  i++;
      ircon(1,i)=6; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=9;  i++;
      ircon(1,i)=4; ircon(2,i)=1;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=2; irconback(2,i)=7;  i++;
      irconback(1,i)=3; irconback(2,i)=7;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //top
      bz(i,1)=1; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=1; bz(i,2)=eta; bz(i,3)=eta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5-zzeta; bz(i,3)=-zzeta; i++;
      bz(i,1)=0; bz(i,2)=eta; bz(i,3)=eta-1.0; i++;
      bz(i,1)=0; bz(i,2)=1.0-eta; bz(i,3)=-eta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5+zzeta; bz(i,3)=zzeta; i++;
      //bottom
      bz(i,1)=-1; bz(i,2)=eta-1.0; bz(i,3)=eta-1.0; i++;
      bz(i,1)=-1; bz(i,2)=-eta; bz(i,3)=-eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5-zzeta); bz(i,3)=zzeta; i++;
      bz(i,1)=0; bz(i,2)=-eta; bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=0; bz(i,2)=-(1.0-eta); bz(i,3)=eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5+zzeta); bz(i,3)=-zzeta; i++;
      //middle
      bz(i,1)=0.5; bz(i,2)=0.5-zzeta; bz(i,3)=1.0-zzeta; i++;
      bz(i,1)=0.5; bz(i,2)=zzeta-0.5; bz(i,3)=zzeta; i++;
      bz(i,1)=0; bz(i,2)=-eta; bz(i,3)=-eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5-zzeta); bz(i,3)=-(1.0-zzeta); i++;
      bz(i,1)=-0.5; bz(i,2)=-(zzeta-0.5); bz(i,3)=-zzeta; i++;
      bz(i,1)=0; bz(i,2)=eta; bz(i,3)=eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=2; con(2,i)=3;  i++;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=9; con(2,i)=10;  i++;
      con(1,i)=9; con(2,i)=18;  i++;
      con(1,i)=6; con(2,i)=18;  i++;
      con(1,i)=3; con(2,i)=15;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      con(1,i)=11; con(2,i)=14;  i++;
      con(1,i)=12; con(2,i)=15;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=10; con(2,i)=19;  Ncon=i;
      i=1;
      conback(1,i)=8; conback(2,i)=9;  i++;
      conback(1,i)=4; conback(2,i)=16;  i++;
      conback(1,i)=8; conback(2,i)=13;  i++;
      conback(1,i)=8; conback(2,i)=17;  i++;
      conback(1,i)=5; conback(2,i)=17;  i++;
      conback(1,i)=12; conback(2,i)=13;  i++;
      conback(1,i)=15; conback(2,i)=16;  i++;
      conback(1,i)=16; conback(2,i)=17;  i++;
      conback(1,i)=13; conback(2,i)=16;  Nconback=i;
    }
    //----------- ORCF2 --------------------------
    if(LattVar=="ORCF2" || LattVar=="orcf2") {
      found=true;
      glatt="ORCF_2";
      xrotview=80; zrotview=115;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      eta=0.25*(1.0+a*a/b/b-a*a/c/c); phi=0.25*(1.0+c*c/b/b-c*c/a/a);
      delta=0.25*(1.0+b*b/a/a-b*b/c/c);
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5-eta; irbz(i,3)=1.0-eta; irbzlab[i]="C"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5+eta; irbz(i,3)=eta; irbzlab[i]="C_1"; i++;
      irbz(i,1)=0.5-delta; irbz(i,2)=0.5; irbz(i,3)=1.0-delta; irbzlab[i]="D"; i++;
      irbz(i,1)=0.5+delta; irbz(i,2)=0.5; irbz(i,3)=delta; irbzlab[i]="D_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="L"; i++;
      irbz(i,1)=1.0-phi; irbz(i,2)=0.5-phi; irbz(i,3)=0.5; irbzlab[i]="H"; i++;
      irbz(i,1)=phi; irbz(i,2)=0.5+phi; irbz(i,3)=0.5; irbzlab[i]="H_1"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="X"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Y"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-C-D-X-{/Symbol G}-Z-D_1-H-C|C_1-Z|X-H_1|H-Y|L-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=2; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=9;  i++;
      ircon(1,i)=9; ircon(2,i)=1;  i++;
      ircon(1,i)=5; ircon(2,i)=7;  i++;
      ircon(1,i)=7; ircon(2,i)=2;  i++;
      ircon(1,i)=9; ircon(2,i)=8;  i++;
      ircon(1,i)=6; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=10;  i++;
      ircon(1,i)=2; ircon(2,i)=10;  i++;
      ircon(1,i)=1; ircon(2,i)=11;  i++;
      ircon(1,i)=5; ircon(2,i)=11;  i++;
      ircon(1,i)=3; ircon(2,i)=11;  i++;
      ircon(1,i)=7; ircon(2,i)=10;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=3; irconback(2,i)=8;  i++;
      irconback(1,i)=4; irconback(2,i)=8;  i++;
      irconback(1,i)=3; irconback(2,i)=5;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //top
      bz(i,1)=0.5+delta; bz(i,2)=0.5; bz(i,3)=delta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5-eta; bz(i,3)=-eta; i++;
      bz(i,1)=0.5-delta; bz(i,2)=0.5; bz(i,3)=-delta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5+eta; bz(i,3)=eta; i++;
      //bottom
      bz(i,1)=-(0.5+delta); bz(i,2)=-0.5; bz(i,3)=-delta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5-eta); bz(i,3)=eta; i++;
      bz(i,1)=-(0.5-delta); bz(i,2)=-0.5; bz(i,3)=delta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5+eta); bz(i,3)=-eta; i++;
      //front
      bz(i,1)=-phi; bz(i,2)=0.5-phi; bz(i,3)=0.5; i++;
      bz(i,1)=0.5-delta; bz(i,2)=0.5; bz(i,3)=1.0-delta; i++;
      bz(i,1)=phi; bz(i,2)=phi+0.5; bz(i,3)=0.5; i++;
      bz(i,1)=delta-0.5; bz(i,2)=0.5; bz(i,3)=delta; i++;
      //back
      bz(i,1)=phi; bz(i,2)=-(0.5-phi); bz(i,3)=-0.5; i++;
      bz(i,1)=-(0.5-delta); bz(i,2)=-0.5; bz(i,3)=-(1.0-delta); i++;
      bz(i,1)=-phi; bz(i,2)=-(phi+0.5); bz(i,3)=-0.5; i++;
      bz(i,1)=-(delta-0.5); bz(i,2)=-0.5; bz(i,3)=-delta; i++;
      //left
      bz(i,1)=phi-1.0; bz(i,2)=phi-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=eta-0.5; bz(i,3)=eta-1.0; i++;
      bz(i,1)=-phi; bz(i,2)=0.5-phi; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5-eta; bz(i,3)=-eta; i++;
      //right
      bz(i,1)=-(phi-1.0); bz(i,2)=-(phi-0.5); bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-(eta-0.5); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=phi; bz(i,2)=-(0.5-phi); bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-(0.5-eta); bz(i,3)=eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=2; con(2,i)=3;  i++;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=6; con(2,i)=18;  i++;
      con(1,i)=7; con(2,i)=10;  i++;
      con(1,i)=4; con(2,i)=20;  i++;
      con(1,i)=8; con(2,i)=24;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      con(1,i)=13; con(2,i)=10;  i++;
      con(1,i)=21; con(2,i)=18;  i++;
      con(1,i)=21; con(2,i)=13;  i++;
      con(1,i)=20; con(2,i)=21;  i++;
      con(1,i)=23; con(2,i)=24;  i++;
      con(1,i)=24; con(2,i)=25;  i++;
      con(1,i)=25; con(2,i)=22;  Ncon=i;
      i=1;
      conback(1,i)=8; conback(2,i)=9;  i++;
      conback(1,i)=9; conback(2,i)=6;  i++;
      conback(1,i)=3; conback(2,i)=14;  i++;
      conback(1,i)=9; conback(2,i)=16;  i++;
      conback(1,i)=14; conback(2,i)=15;  i++;
      conback(1,i)=15; conback(2,i)=16;  i++;
      conback(1,i)=16; conback(2,i)=17;  i++;
      conback(1,i)=17; conback(2,i)=14;  i++;
      conback(1,i)=18; conback(2,i)=19;  i++;
      conback(1,i)=15; conback(2,i)=19;  i++;
      conback(1,i)=19; conback(2,i)=20;  i++;
      conback(1,i)=17; conback(2,i)=25;  Nconback=i;
    }
    //----------- ORCF3 --------------------------
    if(LattVar=="ORCF3" || LattVar=="orcf3") {
      found=true;
      glatt="ORCF_3";
      xrotview=80; zrotview=120;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      zzeta=0.25*(1.0+a*a/b/b-a*a/c/c); eta=0.25*(1.0+a*a/b/b+a*a/c/c);
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5+zzeta; irbz(i,3)=zzeta; irbzlab[i]="A"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5-zzeta; irbz(i,3)=1.0-zzeta; irbzlab[i]="A_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="L"; i++;
      irbz(i,1)=1; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="T"; i++;
      irbz(i,1)=0; irbz(i,2)=eta; irbz(i,3)=eta; irbzlab[i]="X"; i++;
      irbz(i,1)=1; irbz(i,2)=1.0-eta; irbz(i,3)=1.0-eta; irbzlab[i]=""; i++;//"X_1"; i++
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Y"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-T-Z-{/Symbol G}-X-A_1-Y|X-A-Z|L-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=8;  i++;
      ircon(1,i)=8; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=9;  i++;
      ircon(1,i)=9; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=8;  i++;
      ircon(1,i)=6; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=9;  i++;
      ircon(1,i)=4; ircon(2,i)=1;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=2; irconback(2,i)=5;  i++;
      irconback(1,i)=3; irconback(2,i)=5;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //top
      bz(i,1)=1; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=1; bz(i,2)=eta; bz(i,3)=eta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5-zzeta; bz(i,3)=-zzeta; i++;
      bz(i,1)=0; bz(i,2)=eta; bz(i,3)=eta-1.0; i++;
      bz(i,1)=0; bz(i,2)=1.0-eta; bz(i,3)=-eta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5+zzeta; bz(i,3)=zzeta; i++;
      //bottom
      bz(i,1)=-1; bz(i,2)=eta-1.0; bz(i,3)=eta-1.0; i++;
      bz(i,1)=-1; bz(i,2)=-eta; bz(i,3)=-eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5-zzeta); bz(i,3)=zzeta; i++;
      bz(i,1)=0; bz(i,2)=-eta; bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=0; bz(i,2)=-(1.0-eta); bz(i,3)=eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5+zzeta); bz(i,3)=-zzeta; i++;
      //middle
      bz(i,1)=0.5; bz(i,2)=0.5-zzeta; bz(i,3)=1.0-zzeta; i++;
      bz(i,1)=0.5; bz(i,2)=zzeta-0.5; bz(i,3)=zzeta; i++;
      bz(i,1)=0; bz(i,2)=-eta; bz(i,3)=-eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(0.5-zzeta); bz(i,3)=-(1.0-zzeta); i++;
      bz(i,1)=-0.5; bz(i,2)=-(zzeta-0.5); bz(i,3)=-zzeta; i++;
      bz(i,1)=0; bz(i,2)=eta; bz(i,3)=eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=9; con(2,i)=10;  i++;
      con(1,i)=3; con(2,i)=15;  i++;
      con(1,i)=9; con(2,i)=18;  i++;
      con(1,i)=6; con(2,i)=18;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      con(1,i)=11; con(2,i)=14;  i++;
      con(1,i)=12; con(2,i)=15;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=10; con(2,i)=19;  Ncon=i;
      i=1;
      conback(1,i)=2; conback(2,i)=3;  i++;
      conback(1,i)=8; conback(2,i)=9;  i++;
      conback(1,i)=8; conback(2,i)=13;  i++;
      conback(1,i)=4; conback(2,i)=16;  i++;
      conback(1,i)=8; conback(2,i)=17;  i++;
      conback(1,i)=5; conback(2,i)=17;  i++;
      conback(1,i)=12; conback(2,i)=13;  i++;
      conback(1,i)=15; conback(2,i)=16;  i++;
      conback(1,i)=16; conback(2,i)=17;  i++;
      conback(1,i)=13; conback(2,i)=16;  Nconback=i;
    }
    //----------- HEX --------------------------
    if(LattVar=="HEX" || LattVar=="hex") {
      found=true;
      glatt="HEX";
      xrotview=80; zrotview=100;
      azview=100; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="A"; i++;
      irbz(i,1)=1/3.0; irbz(i,2)=1/3.0; irbz(i,3)=0.5; irbzlab[i]="H"; i++;
      irbz(i,1)=1/3.0; irbz(i,2)=1/3.0; irbz(i,3)=0; irbzlab[i]="K"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="L"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="M"; Nirbz=i;
      gpath="{/Symbol G}-M-K-{/Symbol G}-A-L-H|L-M|K-H";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=2;  i++;
      ircon(1,i)=2; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=2;  i++;
      ircon(1,i)=5; ircon(2,i)=6;  i++;
      ircon(1,i)=4; ircon(2,i)=3;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      //top
      bz(i,1)=0.3333; bz(i,2)=0.3333; bz(i,3)=0.5; i++;
      bz(i,1)=-0.3333; bz(i,2)=0.6667; bz(i,3)=0.5; i++;
      bz(i,1)=-0.6667; bz(i,2)=0.3333; bz(i,3)=0.5; i++;
      bz(i,1)=-0.3333; bz(i,2)=-0.3333; bz(i,3)=0.5; i++;
      bz(i,1)=0.3333; bz(i,2)=-0.6667; bz(i,3)=0.5; i++;
      bz(i,1)=0.6667; bz(i,2)=-0.3333; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0.5; i++;
      //bottom
      bz(i,1)=-0.3333; bz(i,2)=0.6667; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.6667; bz(i,2)=0.3333; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.3333; bz(i,2)=-0.3333; bz(i,3)=-0.5; i++;
      bz(i,1)=0.3333; bz(i,2)=-0.6667; bz(i,3)=-0.5; i++;
      bz(i,1)=0.6667; bz(i,2)=-0.3333; bz(i,3)=-0.5; i++;
      bz(i,1)=0.3333; bz(i,2)=0.3333; bz(i,3)=-0.5; i++;
      bz(i,1)=0.3333; bz(i,2)=0.3333; bz(i,3)=0; Nbz=i;
      //connection index
      i=1;
      con(1,i)=2; con(2,i)=3;  i++;
      con(1,i)=3; con(2,i)=4;  i++;
      con(1,i)=4; con(2,i)=5;  i++;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=3; con(2,i)=9;  i++;
      con(1,i)=9; con(2,i)=10;  i++;
      con(1,i)=4; con(2,i)=10;  i++;
      con(1,i)=7; con(2,i)=13;  i++;
      con(1,i)=9; con(2,i)=14;  i++;
      con(1,i)=13; con(2,i)=14;  i++;
      con(1,i)=14; con(2,i)=15;  Ncon=i;
      i=1;
      conback(1,i)=5; conback(2,i)=11;  i++;
      conback(1,i)=6; conback(2,i)=12;  i++;
      conback(1,i)=10; conback(2,i)=11;  i++;
      conback(1,i)=11; conback(2,i)=12;  i++;
      conback(1,i)=12; conback(2,i)=13;  Nconback=i;
    }
    //----------- RHL1 --------------------------
    if(LattVar=="RHL1" || LattVar=="rhl1") {
      double pb1theta;
      found=true;
      glatt="RHL_1";
      xrotview=75; zrotview=125;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;

      //-----rotating b1 b2 b3 so that
      // pL --> kz
      // b1 --> on kxkz-plane
      pL(1)=0.5; pL(2)=0.5; pL(3)=0.5;
      pL=pL*klattice;
      pLtheta=acos(pL(3)/modulus(pL));//=acos(scalar_product(pL,pz)/modulus(pL))
      ptmp(1)=pL(2); ptmp(2)=-pL(1); ptmp(3)=0; //=cross(pL,pz);
      pLaxis=ptmp/modulus(ptmp);
      tmpb1=Vrotate(b1,pLaxis,pLtheta);
      tmpb2=Vrotate(b2,pLaxis,pLtheta);
      tmpb3=Vrotate(b3,pLaxis,pLtheta);
      pb1theta=atan(abs(tmpb1(2)/tmpb1(1)));
      if(tmpb1(1)<0 && tmpb1(2)>0) pb1theta=pi-pb1theta;
      if(tmpb1(1)<0 && tmpb1(2)<0) pb1theta=pi+pb1theta;
      if(tmpb1(1)>0 && tmpb1(2)<0) pb1theta=-pb1theta;
      pLaxis(1)=0; pLaxis(2)=0; pLaxis(3)=-1;
      b1=Vrotate(tmpb1,pLaxis,pb1theta);
      b2=Vrotate(tmpb2,pLaxis,pb1theta);
      b3=Vrotate(tmpb3,pLaxis,pb1theta);
      for(j=1;j<4;j++) {
	klattice(1,j)=b1(j);
	klattice(2,j)=b2(j);
	klattice(3,j)=b3(j);
      }

      eta=(1.0+4*cos(alpha))/(2.0+4*cos(alpha));
      nu=0.75-0.5*eta;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=eta; irbz(i,2)=0.5; irbz(i,3)=1.0-eta; irbzlab[i]="B"; i++;
      irbz(i,1)=0.5; irbz(i,2)=1.0-eta; irbz(i,3)=eta-1.0; irbzlab[i]="B_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="F"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="L"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=-0.5; irbzlab[i]="L_1"; i++;
      irbz(i,1)=eta; irbz(i,2)=nu; irbz(i,3)=nu; irbzlab[i]="P"; i++;
      irbz(i,1)=1.0-nu; irbz(i,2)=1.0-nu; irbz(i,3)=1.0-eta; irbzlab[i]="P_1"; i++;
      irbz(i,1)=nu; irbz(i,2)=nu; irbz(i,3)=eta-1.0; irbzlab[i]="P_2"; i++;
      irbz(i,1)=1.0-nu; irbz(i,2)=nu; irbz(i,3)=0; irbzlab[i]="Q"; i++;
      irbz(i,1)=nu; irbz(i,2)=0; irbz(i,3)=-nu; irbzlab[i]="X"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-L-B_1|B-Z-{/Symbol G}-X|Q-F-P_1-Z|L-P";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=5;  i++;
      ircon(1,i)=5; ircon(2,i)=3;  i++;
      ircon(1,i)=5; ircon(2,i)=7;  i++;
      ircon(1,i)=4; ircon(2,i)=8;  i++;
      ircon(1,i)=2; ircon(2,i)=1;  i++;
      ircon(1,i)=2; ircon(2,i)=12;  i++;
      ircon(1,i)=1; ircon(2,i)=11;  i++;
      ircon(1,i)=4; ircon(2,i)=10;  i++;
      ircon(1,i)=8; ircon(2,i)=12;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=6; irconback(2,i)=9;  i++;
      irconback(1,i)=9; irconback(2,i)=4;  i++;
      irconback(1,i)=6; irconback(2,i)=1;  i++;
      irconback(1,i)=3; irconback(2,i)=9;  i++;
      irconback(1,i)=3; irconback(2,i)=2;  i++;
      irconback(1,i)=2; irconback(2,i)=7;  i++;
      irconback(1,i)=2; irconback(2,i)=8;  i++;
      irconback(1,i)=3; irconback(2,i)=11;  i++;
      irconback(1,i)=5; irconback(2,i)=11;  i++;
      irconback(1,i)=6; irconback(2,i)=11;  i++;
      irconback(1,i)=7; irconback(2,i)=12;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=nu; bz(i,2)=0; bz(i,3)=-nu; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=1.0-eta; bz(i,3)=eta-1.0; i++;
      bz(i,1)=nu; bz(i,2)=nu; bz(i,3)=eta-1.0; i++;
      bz(i,1)=1.0-nu; bz(i,2)=nu; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=1.0-nu; bz(i,2)=1.0-nu; bz(i,3)=1.0-eta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=eta; bz(i,2)=nu; bz(i,3)=nu; i++;
      bz(i,1)=eta; bz(i,2)=0.5; bz(i,3)=1.0-eta; i++;
      bz(i,1)=0.5; bz(i,2)=eta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=1.0-eta; bz(i,2)=eta; bz(i,3)=0.5; i++;
      bz(i,1)=1.0-eta; bz(i,2)=0.5; bz(i,3)=eta; i++;
      bz(i,1)=0.5; bz(i,2)=1.0-eta; bz(i,3)=eta; i++;
      bz(i,1)=eta; bz(i,2)=1.0-eta; bz(i,3)=0.5; i++;
      bz(i,1)=1.0-eta; bz(i,2)=0.5; bz(i,3)=eta-1.0; i++;
      bz(i,1)=eta-1.0; bz(i,2)=1.0-eta; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=1.0-eta; bz(i,3)=eta-1.0; i++;
      bz(i,1)=eta-1.0; bz(i,2)=0.5; bz(i,3)=1.0-eta; i++;
      bz(i,1)=eta-1.0; bz(i,2)=1.0-eta; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-(1.0-eta); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(1.0-eta); bz(i,2)=-0.5; bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(eta-1.0); bz(i,2)=-(1.0-eta); bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=-(1.0-eta); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(eta-1.0); bz(i,2)=-0.5; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-(eta-1.0); bz(i,2)=-(1.0-eta); bz(i,3)=-0.5; i++;
      bz(i,1)=-(1.0-eta); bz(i,2)=-0.5; bz(i,3)=-eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-(1.0-eta); bz(i,3)=-eta; i++;
      bz(i,1)=-eta; bz(i,2)=-(1.0-eta); bz(i,3)=-0.5; i++;
      bz(i,1)=-eta; bz(i,2)=-0.5; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-0.5; bz(i,2)=-eta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-(1.0-eta); bz(i,2)=-eta; bz(i,3)=-0.5; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=30;  i++;
      con(1,i)=8; con(2,i)=20;  i++;
      con(1,i)=11; con(2,i)=15;  i++;
      con(1,i)=16; con(2,i)=17;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=19; con(2,i)=13;  i++;
      con(1,i)=15; con(2,i)=16;  i++;
      con(1,i)=17; con(2,i)=18;  i++;
      con(1,i)=15; con(2,i)=20;  i++;
      con(1,i)=16; con(2,i)=23;  i++;
      con(1,i)=18; con(2,i)=27;  i++;
      con(1,i)=19; con(2,i)=28;  i++;
      con(1,i)=21; con(2,i)=22;  i++;
      con(1,i)=20; con(2,i)=21;  i++;
      con(1,i)=22; con(2,i)=23;  i++;
      con(1,i)=27; con(2,i)=28;  i++;
      con(1,i)=28; con(2,i)=29;  i++;
      con(1,i)=29; con(2,i)=36;  i++;
      con(1,i)=29; con(2,i)=30;  i++;
      con(1,i)=21; con(2,i)=32;  i++;
      con(1,i)=32; con(2,i)=31;  i++;
      con(1,i)=36; con(2,i)=31;  i++;
      con(1,i)=31; con(2,i)=30;  Ncon=i;
      i=1;
      conback(1,i)=17; conback(2,i)=24;  i++;
      conback(1,i)=22; conback(2,i)=33;  i++;
      conback(1,i)=25; conback(2,i)=24;  i++;
      conback(1,i)=24; conback(2,i)=23;  i++;
      conback(1,i)=25; conback(2,i)=26;  i++;
      conback(1,i)=26; conback(2,i)=27;  i++;
      conback(1,i)=26; conback(2,i)=35;  i++;
      conback(1,i)=22; conback(2,i)=33;  i++;
      conback(1,i)=25; conback(2,i)=34;  i++;
      conback(1,i)=36; conback(2,i)=35;  i++;
      conback(1,i)=35; conback(2,i)=34;  i++;
      conback(1,i)=34; conback(2,i)=33;  i++;
      conback(1,i)=33; conback(2,i)=32;  Nconback=i;
    }
    //----------- RHL2 --------------------------
    if(LattVar=="RHL2" || LattVar=="rhl2") {
      double pb1theta;
      found=true;
      glatt="RHL_2";
      xrotview=85; zrotview=40;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;

      //-----rotating b1 b2 b3 so that
      // pL --> kz
      // b1 --> on kxkz-plane
      pL(1)=0.5; pL(2)=0.5; pL(3)=0.5;
      pL=pL*klattice;
      pLtheta=acos(pL(3)/modulus(pL));//=acos(scalar_product(pL,pz)/modulus(pL))
      ptmp(1)=pL(2); ptmp(2)=-pL(1); ptmp(3)=0; //=cross(pL,pz);
      pLaxis=ptmp/modulus(ptmp);
      tmpb1=Vrotate(b1,pLaxis,pLtheta);
      tmpb2=Vrotate(b2,pLaxis,pLtheta);
      tmpb3=Vrotate(b3,pLaxis,pLtheta);
      pb1theta=atan(abs(tmpb1(2)/tmpb1(1)));
      if(tmpb1(1)<0 && tmpb1(2)>0) pb1theta=pi-pb1theta;
      if(tmpb1(1)<0 && tmpb1(2)<0) pb1theta=pi+pb1theta;
      if(tmpb1(1)>0 && tmpb1(2)<0) pb1theta=-pb1theta;
      pLaxis(1)=0; pLaxis(2)=0; pLaxis(3)=-1;
      b1=Vrotate(tmpb1,pLaxis,pb1theta);
      b2=Vrotate(tmpb2,pLaxis,pb1theta);
      b3=Vrotate(tmpb3,pLaxis,pb1theta);
      for(j=1;j<4;j++) {
	klattice(1,j)=b1(j);
	klattice(2,j)=b2(j);
	klattice(3,j)=b3(j);
      }

      eta=1.0/(2*tan(alpha/2.0)*tan(alpha/2.0));
      nu=0.75-0.5*eta;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="F"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="L"; i++;
      irbz(i,1)=1.0-nu; irbz(i,2)=-nu; irbz(i,3)=1.0-nu; irbzlab[i]="P"; i++;
      irbz(i,1)=nu; irbz(i,2)=nu-1.0; irbz(i,3)=nu-1.0; irbzlab[i]="P_1"; i++;
      irbz(i,1)=eta; irbz(i,2)=eta; irbz(i,3)=eta; irbzlab[i]="Q"; i++;
      irbz(i,1)=1.0-eta; irbz(i,2)=-eta; irbz(i,3)=-eta; irbzlab[i]="Q_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-P-Z-Q-{/Symbol G}-F-P_1-Q_1-L-Z";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=4;  i++;
      ircon(1,i)=4; ircon(2,i)=8;  i++;
      ircon(1,i)=5; ircon(2,i)=7;  i++;
      ircon(1,i)=8; ircon(2,i)=6;  i++;
      ircon(1,i)=6; ircon(2,i)=1;  i++;
      ircon(1,i)=2; ircon(2,i)=5;  i++;
      ircon(1,i)=1; ircon(2,i)=2;  i++;
      ircon(1,i)=7; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=8;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=1; irconback(2,i)=5;  i++;
      irconback(1,i)=4; irconback(2,i)=2;  i++;
      irconback(1,i)=7; irconback(2,i)=8;  i++;
      irconback(1,i)=3; irconback(2,i)=6;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=1.0-nu; bz(i,2)=-nu; bz(i,3)=1.0-nu; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0; i++;
      bz(i,1)=nu; bz(i,2)=nu-1.0; bz(i,3)=nu-1.0; i++;
      bz(i,1)=1.0-eta; bz(i,2)=-eta; bz(i,3)=-eta; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=eta; bz(i,2)=eta; bz(i,3)=eta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-(1.0-eta); bz(i,2)=eta; bz(i,3)=eta; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=-eta; bz(i,2)=-eta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=eta; bz(i,2)=eta-1.0; bz(i,3)=eta; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=-eta; bz(i,2)=-eta; bz(i,3)=-eta; i++;
      bz(i,1)=-0.5; bz(i,2)=0.5; bz(i,3)=-0.5; i++;
      bz(i,1)=eta; bz(i,2)=eta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-eta; bz(i,2)=-(eta-1.0); bz(i,3)=-eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=9;  i++;
      con(1,i)=8; con(2,i)=11;  i++;
      con(1,i)=9; con(2,i)=15;  i++;
      con(1,i)=5; con(2,i)=16;  i++;
      con(1,i)=7; con(2,i)=17;  i++;
      con(1,i)=18; con(2,i)=17;  i++;
      con(1,i)=17; con(2,i)=16;  i++;
      con(1,i)=16; con(2,i)=14;  i++;
      con(1,i)=14; con(2,i)=18;  i++;
      con(1,i)=11; con(2,i)=10;  i++;
      con(1,i)=10; con(2,i)=12;  i++;
      con(1,i)=15; con(2,i)=12;  i++;
      con(1,i)=14; con(2,i)=15;  i++;
      con(1,i)=17; con(2,i)=20;  i++;
      con(1,i)=11; con(2,i)=20;  Ncon=i;
      i=1;
      conback(1,i)=12; conback(2,i)=13;  i++;
      conback(1,i)=13; conback(2,i)=14;  i++;
      conback(1,i)=13; conback(2,i)=19;  i++;
      conback(1,i)=19; conback(2,i)=18;  i++;
      conback(1,i)=19; conback(2,i)=20;  i++;
      conback(1,i)=19; conback(2,i)=21;  i++;
      conback(1,i)=11; conback(2,i)=21;  i++;
      conback(1,i)=12; conback(2,i)=21;  Nconback=i;
    }
    //----------- MCL --------------------------
    if(LattVar=="MCL" || LattVar=="mcl") {
      double pb2theta;
      found=true;
      glatt="MCL";
      xrotview=80; zrotview=115;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;

      //-----rotating b1 b2 b3 so that
      // pL --> kz
      // b2 --> on kxkz-plane
      pL(1)=0.5; pL(2)=0; pL(3)=0;
      pL=pL*klattice;
      pLtheta=acos(pL(3)/modulus(pL));//=acos(scalar_product(pL,pz)/modulus(pL))
      ptmp(1)=pL(2); ptmp(2)=-pL(1); ptmp(3)=0; //=cross(pL,pz);
      pLaxis=ptmp/modulus(ptmp);
      tmpb1=Vrotate(b1,pLaxis,pLtheta);
      tmpb2=Vrotate(b2,pLaxis,pLtheta);
      tmpb3=Vrotate(b3,pLaxis,pLtheta);
      pb2theta=atan(abs(tmpb2(2)/tmpb2(1)));
      if(tmpb2(1)<0 && tmpb2(2)>0) pb2theta=pi-pb2theta;
      if(tmpb2(1)<0 && tmpb2(2)<0) pb2theta=pi+pb2theta;
      if(tmpb2(1)>0 && tmpb2(2)<0) pb2theta=-pb2theta;
      pLaxis(1)=0; pLaxis(2)=0; pLaxis(3)=-1;
      b1=Vrotate(tmpb1,pLaxis,pb2theta);
      b2=Vrotate(tmpb2,pLaxis,pb2theta);
      b3=Vrotate(tmpb3,pLaxis,pb2theta);
      for(j=1;j<4;j++) {
	klattice(1,j)=b1(j);
	klattice(2,j)=b2(j);
	klattice(3,j)=b3(j);
      }

      eta=(1.0-b/c*cos(alpha))/(2*sin(alpha)*sin(alpha));
      nu=0.5-eta*c/b*cos(alpha);
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="A"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="C"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="D"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=-0.5; irbzlab[i]="D_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="E"; i++;
      irbz(i,1)=0; irbz(i,2)=eta; irbz(i,3)=1.0-nu; irbzlab[i]="H"; i++;
      irbz(i,1)=0; irbz(i,2)=1.0-eta; irbz(i,3)=nu; irbzlab[i]="H_1"; i++;
      irbz(i,1)=0; irbz(i,2)=eta; irbz(i,3)=-nu; irbzlab[i]="H_2"; i++;
      irbz(i,1)=0.5; irbz(i,2)=eta; irbz(i,3)=1.0-nu; irbzlab[i]="M"; i++;
      irbz(i,1)=0.5; irbz(i,2)=1.0-eta; irbz(i,3)=nu; irbzlab[i]="M_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=eta; irbz(i,3)=-nu; irbzlab[i]="M_2"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Y"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=-0.5; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-H-C-E-M_1-A-X-{/Symbol G}-Z-D-M|Z-A|D-Y|X-H_1";
      //old: gpath="{/Symbol G}-Y-H-C-E-M_1-A-X-H_1|M-D-Z|Y-D";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=13;  i++;
      ircon(1,i)=1; ircon(2,i)=16;  i++;
      ircon(1,i)=2; ircon(2,i)=16;  i++;
      ircon(1,i)=7; ircon(2,i)=3;  i++;
      ircon(1,i)=3; ircon(2,i)=6;  i++;
      ircon(1,i)=1; ircon(2,i)=14;  i++;
      ircon(1,i)=7; ircon(2,i)=14;  i++;
      ircon(1,i)=6; ircon(2,i)=11;  i++;
      ircon(1,i)=2; ircon(2,i)=11;  i++;
      ircon(1,i)=2; ircon(2,i)=13;  i++;
      ircon(1,i)=8; ircon(2,i)=13;  i++;
      ircon(1,i)=4; ircon(2,i)=10;  i++;
      ircon(1,i)=4; ircon(2,i)=16;  i++;
      ircon(1,i)=4; ircon(2,i)=14;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=3; irconback(2,i)=8;  i++;
      irconback(1,i)=7; irconback(2,i)=10;  i++;
      irconback(1,i)=6; irconback(2,i)=10;  i++;
      irconback(1,i)=8; irconback(2,i)=11;  i++;
      irconback(1,i)=9; irconback(2,i)=13;  i++;
      irconback(1,i)=9; irconback(2,i)=15;  i++;
      irconback(1,i)=1; irconback(2,i)=15;  i++;
      irconback(1,i)=2; irconback(2,i)=12;  i++;
      irconback(1,i)=5; irconback(2,i)=12;  i++;
      irconback(1,i)=5; irconback(2,i)=16;  i++;
      irconback(1,i)=5; irconback(2,i)=15;  i++;
      irconback(1,i)=9; irconback(2,i)=12;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++; //Z
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++; //X
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=-0.5; i++;
      bz(i,1)=0; bz(i,2)=eta; bz(i,3)=-nu; i++;
      bz(i,1)=0; bz(i,2)=1.0-eta; bz(i,3)=nu; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0; bz(i,2)=eta; bz(i,3)=1.0-nu; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=-0.5; i++;
      bz(i,1)=0.5; bz(i,2)=eta; bz(i,3)=-nu; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0; i++; //A
      bz(i,1)=0.5; bz(i,2)=1.0-eta; bz(i,3)=nu; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0.5; i++;
      bz(i,1)=0.5; bz(i,2)=eta; bz(i,3)=1.0-nu; i++;
      bz(i,1)=0.5; bz(i,2)=-eta; bz(i,3)=nu; i++;
      bz(i,1)=0.5; bz(i,2)=eta-1.0; bz(i,3)=-nu; i++;
      bz(i,1)=0.5; bz(i,2)=-eta; bz(i,3)=nu-1.0; i++;
      bz(i,1)=-0.5; bz(i,2)=-eta; bz(i,3)=-(1.0-nu); i++;
      bz(i,1)=-0.5; bz(i,2)=eta; bz(i,3)=-nu; i++;
      bz(i,1)=-0.5; bz(i,2)=-(eta-1.0); bz(i,3)=nu; i++;
      bz(i,1)=-0.5; bz(i,2)=eta; bz(i,3)=-(nu-1.0); i++;
      bz(i,1)=-0.5; bz(i,2)=-eta; bz(i,3)=nu; i++;
      bz(i,1)=-0.5; bz(i,2)=-(1.0-eta); bz(i,3)=-nu; Nbz=i;
      //connection index
      i=1;
      con(1,i)=6; con(2,i)=21;  i++;
      con(1,i)=7; con(2,i)=22;  i++;
      con(1,i)=9; con(2,i)=23;  i++;
      con(1,i)=19; con(2,i)=11;  i++;
      con(1,i)=10; con(2,i)=17;  i++;
      con(1,i)=17; con(2,i)=18;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=17; con(2,i)=24;  i++;
      con(1,i)=21; con(2,i)=22;  i++;
      con(1,i)=22; con(2,i)=23;  i++;
      con(1,i)=23; con(2,i)=24;  Ncon=i;
      i=1;
      conback(1,i)=18; conback(2,i)=25;  i++;
      conback(1,i)=19; conback(2,i)=20;  i++;
      conback(1,i)=21; conback(2,i)=20;  i++;
      conback(1,i)=20; conback(2,i)=25;  i++;
      conback(1,i)=25; conback(2,i)=24;  Nconback=i;
    }
    //----------- MCLC1 --------------------------
    if(LattVar=="MCLC1" || LattVar=="mclc1") {
      found=true;
      glatt="MCLC_1";
      xrotview=75; zrotview=170;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      zzeta=(2.0-b/c*cos(alpha))/(4*sin(alpha)*sin(alpha));
      eta=0.5+2*zzeta*c/b*cos(alpha);
      psi=0.75-a*a/(4*b*b*sin(alpha)*sin(alpha));
      phi=psi+(0.75-psi)*b/c*cos(alpha);
      mu=psi+0.25*(2-b/c/cos(alpha))/tan(alpha)/tan(alpha);
      delta=1.0-0.5*b/tan(alpha)*(1/c/sin(alpha)-2/b/tan(alpha))-mu;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="N"; i++;
      irbz(i,1)=0; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="N_1"; i++;
      irbz(i,1)=1.0-zzeta; irbz(i,2)=1.0-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="F"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=zzeta; irbz(i,3)=eta; irbzlab[i]="F_1"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="F_2"; i++;
      //irbz(i,1)=1.0-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="F_3"; i++;
      irbz(i,1)=phi; irbz(i,2)=1.0-phi; irbz(i,3)=0.5; irbzlab[i]="I"; i++;
      irbz(i,1)=1.0-phi; irbz(i,2)=phi-1.0; irbz(i,3)=0.5; irbzlab[i]="I_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="L"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="M"; i++;
      irbz(i,1)=1.0-psi; irbz(i,2)=psi-1.0; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=psi; irbz(i,2)=1.0-psi; irbz(i,3)=0; irbzlab[i]="X_1"; i++;
      irbz(i,1)=psi-1.0; irbz(i,2)=-psi; irbz(i,3)=0; irbzlab[i]="X_2"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Y"; i++;
      irbz(i,1)=-0.5; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i; i++;
      irbz(i,1)=1.0-mu; irbz(i,2)=-delta; irbz(i,3)=1.0-eta; irbzlab[i]=""; i++;
      irbz(i,1)=1.0-delta; irbz(i,2)=1.0-mu; irbz(i,3)=1.0-eta; irbzlab[i]=""; i++;
      irbz(i,1)=mu; irbz(i,2)=delta; irbz(i,3)=eta; irbzlab[i]=""; i++;
      irbz(i,1)=-delta; irbz(i,2)=-mu; irbz(i,3)=1.0-eta; irbzlab[i]=""; i++;
      gpath="{/Symbol G}-Y-F-L-I|I_1-Z-{/Symbol G}-X|X_1-Y|M-{/Symbol G}-N|Z-F_1";
      //old gpath="{/Symbol G}-Y-F-L-I|I_1-Z-F_1|Y-X_1|X-{/Symbol G}-N|M-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=2;  i++;
      ircon(1,i)=1; ircon(2,i)=14;  i++;
      ircon(1,i)=1; ircon(2,i)=16;  i++;
      ircon(1,i)=4; ircon(2,i)=14;  i++;
      ircon(1,i)=4; ircon(2,i)=9;  i++;
      ircon(1,i)=7; ircon(2,i)=9;  i++;
      ircon(1,i)=8; ircon(2,i)=16;  i++;
      ircon(1,i)=5; ircon(2,i)=16;  i++;
      ircon(1,i)=12; ircon(2,i)=14;  i++;
      ircon(1,i)=1; ircon(2,i)=11;  i++;
      ircon(1,i)=1; ircon(2,i)=10;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=5; irconback(2,i)=9;  i++;
      irconback(1,i)=1; irconback(2,i)=15;  i++;
      irconback(1,i)=4; irconback(2,i)=18;  i++;
      irconback(1,i)=7; irconback(2,i)=18;  i++;
      irconback(1,i)=7; irconback(2,i)=19;  i++;
      irconback(1,i)=5; irconback(2,i)=19;  i++;
      irconback(1,i)=7; irconback(2,i)=17;  i++;
      irconback(1,i)=8; irconback(2,i)=19;  i++;
      irconback(1,i)=8; irconback(2,i)=17;  i++;
      irconback(1,i)=6; irconback(2,i)=15;  i++;
      irconback(1,i)=6; irconback(2,i)=16;  i++;
      //irconback(1,i)=5; irconback(2,i)=20;  i++;
      irconback(1,i)=6; irconback(2,i)=20;  i++;
      irconback(1,i)=8; irconback(2,i)=20;  i++;
      irconback(1,i)=11; irconback(2,i)=12;  i++;
      irconback(1,i)=11; irconback(2,i)=13;  i++;
      irconback(1,i)=13; irconback(2,i)=15;  i++;
      irconback(1,i)=12; irconback(2,i)=18;  i++;
      irconback(1,i)=11; irconback(2,i)=17;  i++;
      //irconback(1,i)=16; irconback(2,i)=15;  i++;
      irconback(1,i)=13; irconback(2,i)=20;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++; //Z
      bz(i,1)=1.0-phi; bz(i,2)=phi-1.0; bz(i,3)=-0.5; i++;
      bz(i,1)=mu; bz(i,2)=delta; bz(i,3)=eta-1.0; i++;
      bz(i,1)=delta; bz(i,2)=mu; bz(i,3)=eta-1.0; i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=-(phi-1.0); bz(i,3)=-0.5; i++;
      bz(i,1)=-mu; bz(i,2)=-delta; bz(i,3)=-eta; i++;
      bz(i,1)=-delta; bz(i,2)=-mu; bz(i,3)=-eta; i++;
      bz(i,1)=phi-1.0; bz(i,2)=-phi; bz(i,3)=-0.5; i++;
      bz(i,1)=delta; bz(i,2)=mu-1.0; bz(i,3)=eta-1.0; i++;
      bz(i,1)=1.0-psi; bz(i,2)=psi-1.0; bz(i,3)=0; i++;
      bz(i,1)=1.0-mu; bz(i,2)=-delta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=phi; bz(i,2)=1.0-phi; bz(i,3)=0.5; i++;
      bz(i,1)=1.0-delta; bz(i,2)=1.0-mu; bz(i,3)=1.0-eta; i++;
      bz(i,1)=psi; bz(i,2)=1.0-psi; bz(i,3)=0; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=1.0-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=1.0-mu; bz(i,2)=1.0-delta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=-(phi-1.0); bz(i,2)=phi; bz(i,3)=0.5; i++;
      bz(i,1)=-delta; bz(i,2)=-(mu-1.0); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(1.0-mu); bz(i,2)=delta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-phi; bz(i,2)=-(1.0-phi); bz(i,3)=-0.5; i++;
      bz(i,1)=-(1.0-delta); bz(i,2)=-(1.0-mu); bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-mu; bz(i,2)=-delta; bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=-(phi-1.0); bz(i,3)=0.5; i++;
      bz(i,1)=delta; bz(i,2)=mu; bz(i,3)=eta; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=eta; i++;
      bz(i,1)=mu; bz(i,2)=delta; bz(i,3)=eta; i++;
      bz(i,1)=1.0-phi; bz(i,2)=phi-1.0; bz(i,3)=0.5; i++;
      bz(i,1)=-delta; bz(i,2)=-mu; bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=psi-1.0; bz(i,2)=-psi; bz(i,3)=0; i++;
      bz(i,1)=-(1.0-mu); bz(i,2)=-(1.0-delta); bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=6; con(2,i)=17;  i++;
      con(1,i)=5; con(2,i)=12;  i++;
      con(1,i)=7; con(2,i)=19;  i++;
      con(1,i)=8; con(2,i)=22;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      //con(1,i)=13; con(2,i)=14;  i++;
      //con(1,i)=14; con(2,i)=15;  i++;
      //con(1,i)=15; con(2,i)=16;  i++;
      //con(1,i)=16; con(2,i)=17;  i++;
      //con(1,i)=16; con(2,i)=18;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=19; con(2,i)=20;  i++;
      //con(1,i)=15; con(2,i)=29;  i++;
      //con(1,i)=14; con(2,i)=30;  i++;
      con(1,i)=20; con(2,i)=21;  i++;
      con(1,i)=21; con(2,i)=22;  i++;
      //con(1,i)=29; con(2,i)=28;  i++;
      con(1,i)=28; con(2,i)=27;  i++;
      con(1,i)=27; con(2,i)=20;  i++;
      con(1,i)=27; con(2,i)=26;  i++;
      con(1,i)=26; con(2,i)=21;  i++;
      con(1,i)=26; con(2,i)=25;  i++;
      //con(1,i)=29; con(2,i)=30;  i++;
      con(1,i)=25; con(2,i)=34;  Ncon=i;
      //con(1,i)=34; con(2,i)=31;  i++;
      //con(1,i)=31; con(2,i)=30;  Ncon=i;
      i=1;
      conback(1,i)=9; conback(2,i)=8;  i++;
      conback(1,i)=5; conback(2,i)=10;  i++;
      conback(1,i)=9; conback(2,i)=10;  i++;
      conback(1,i)=9; conback(2,i)=23;  i++;
      conback(1,i)=12; conback(2,i)=11;  i++;
      conback(1,i)=10; conback(2,i)=11;  i++;
      conback(1,i)=11; conback(2,i)=33;  i++;
      conback(1,i)=22; conback(2,i)=23;  i++;
      conback(1,i)=23; conback(2,i)=24;  i++;
      conback(1,i)=24; conback(2,i)=25;  i++;
      conback(1,i)=24; conback(2,i)=33;  i++;
      conback(1,i)=33; conback(2,i)=32;  Nconback=i;
    }
    //----------- MCLC2 --------------------------
    if(LattVar=="MCLC2" || LattVar=="mclc2") {
      found=true;
      glatt="MCLC_2";
      xrotview=82; zrotview=165;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      zzeta=(2.0-b/c*cos(alpha))/(4*sin(alpha)*sin(alpha));
      eta=0.5+2*zzeta*c/b*cos(alpha);
      psi=0.75-a*a/(4*b*b*sin(alpha)*sin(alpha));
      phi=psi+(0.75-psi)*b/c*cos(alpha);
      mu=psi+0.25*(2-b/c/cos(alpha))/tan(alpha)/tan(alpha);
      delta=1.0-0.5*b/tan(alpha)*(1/c/sin(alpha)-2/b/tan(alpha))-mu;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="N"; i++;
      irbz(i,1)=0; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="N_1"; i++;
      irbz(i,1)=1.0-zzeta; irbz(i,2)=1.0-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="F"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=zzeta; irbz(i,3)=eta; irbzlab[i]="F_1"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="F_2"; i++;
      irbz(i,1)=1.0-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="F_3"; i++;
      irbz(i,1)=phi; irbz(i,2)=1.0-phi; irbz(i,3)=0.5; irbzlab[i]="I"; i++;
      irbz(i,1)=1.0-phi; irbz(i,2)=phi-1.0; irbz(i,3)=0.5; irbzlab[i]="I_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="L"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="M"; i++;
      irbz(i,1)=1.0-psi; irbz(i,2)=psi-1.0; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=psi; irbz(i,2)=1.0-psi; irbz(i,3)=0; irbzlab[i]=""; i++; //"X_1"; i++;
      irbz(i,1)=psi-1.0; irbz(i,2)=-psi; irbz(i,3)=0; irbzlab[i]=""; i++; //"X_2"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Y"; i++;
      irbz(i,1)=-0.5; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-F-L-I|I_1-Z-{/Symbol G}-M|N-{/Symbol G}|Z-F_1";
      //old gpath="{/Symbol G}-Y-F-L-I|I_1-Z-F_1|N-{/Symbol G}-M";
      //connection table for kpath points
      i=1;
      ircon(1,i)=1; ircon(2,i)=2;  i++;
      ircon(1,i)=1; ircon(2,i)=15;  i++;
      ircon(1,i)=1; ircon(2,i)=17;  i++;
      ircon(1,i)=4; ircon(2,i)=15;  i++;
      ircon(1,i)=4; ircon(2,i)=10;  i++;
      ircon(1,i)=8; ircon(2,i)=10;  i++;
      ircon(1,i)=5; ircon(2,i)=17;  i++;
      ircon(1,i)=9; ircon(2,i)=17;  i++;
      ircon(1,i)=1; ircon(2,i)=2;  i++;
      ircon(1,i)=1; ircon(2,i)=11;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=4; irconback(2,i)=8;  i++;
      irconback(1,i)=5; irconback(2,i)=8;  i++;
      irconback(1,i)=8; irconback(2,i)=7;  i++;
      irconback(1,i)=7; irconback(2,i)=9;  i++;
      irconback(1,i)=9; irconback(2,i)=5;  i++;
      irconback(1,i)=6; irconback(2,i)=9;  i++;
      irconback(1,i)=1; irconback(2,i)=16;  i++;
      irconback(1,i)=7; irconback(2,i)=12;  i++;
      irconback(1,i)=5; irconback(2,i)=10;  i++;
      irconback(1,i)=6; irconback(2,i)=17;  i++;
      irconback(1,i)=6; irconback(2,i)=16;  i++;
      irconback(1,i)=15; irconback(2,i)=12;  i++;
      irconback(1,i)=12; irconback(2,i)=16;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=1.0-phi; bz(i,2)=phi-1.0; bz(i,3)=-0.5; i++;
      bz(i,1)=mu; bz(i,2)=delta; bz(i,3)=eta-1.0; i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=-(phi-1.0); bz(i,3)=-0.5; i++;
      bz(i,1)=-mu; bz(i,2)=-delta; bz(i,3)=-eta; i++;
      bz(i,1)=0.5; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=1.0-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=-(phi-1.0); bz(i,2)=phi; bz(i,3)=0.5; i++;
      bz(i,1)=-delta; bz(i,2)=-(mu-1.0); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(1.0-mu); bz(i,2)=delta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-phi; bz(i,2)=-(1.0-phi); bz(i,3)=-0.5; i++;
      bz(i,1)=-(1.0-zzeta); bz(i,2)=-(1.0-zzeta); bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=phi-1.0; bz(i,2)=-phi; bz(i,3)=-0.5; i++;
      bz(i,1)=delta; bz(i,2)=mu-1.0; bz(i,3)=eta-1.0; i++;
      bz(i,1)=1.0-psi; bz(i,2)=psi-1.0; bz(i,3)=0; i++;
      bz(i,1)=1.0-mu; bz(i,2)=-delta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=phi; bz(i,2)=1.0-phi; bz(i,3)=0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=eta; i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=-(phi-1.0); bz(i,3)=0.5; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=1.0-phi; bz(i,2)=phi-1.0; bz(i,3)=0.5; i++;
      bz(i,1)=-0.5; bz(i,2)=-0.5; bz(i,3)=0; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=6; con(2,i)=9;  i++;
      //con(1,i)=9; con(2,i)=10;  i++;
      con(1,i)=7; con(2,i)=13;  i++;
      con(1,i)=5; con(2,i)=17;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      con(1,i)=17; con(2,i)=18;  i++;
      //con(1,i)=18; con(2,i)=19;  i++;
      //con(1,i)=19; con(2,i)=20;  i++;
      //con(1,i)=10; con(2,i)=20;  i++;
      con(1,i)=11; con(2,i)=21;  i++;
      con(1,i)=12; con(2,i)=22;  i++;
      //con(1,i)=19; con(2,i)=24;  i++;
      con(1,i)=21; con(2,i)=22;  i++;
      con(1,i)=23; con(2,i)=22;  Ncon=i;
      //con(1,i)=20; con(2,i)=21;  i++;
      //con(1,i)=24; con(2,i)=21;  i++;
      //con(1,i)=24; con(2,i)=23;  Ncon=i;
      i=1;
      conback(1,i)=5; conback(2,i)=8;  i++;
      conback(1,i)=8; conback(2,i)=7;  i++;
      conback(1,i)=8; conback(2,i)=16;  i++;
      conback(1,i)=8; conback(2,i)=14;  i++;
      conback(1,i)=16; conback(2,i)=17;  i++;
      conback(1,i)=14; conback(2,i)=13;  i++;
      conback(1,i)=14; conback(2,i)=15;  i++;
      conback(1,i)=15; conback(2,i)=16;  i++;
      conback(1,i)=15; conback(2,i)=25;  Nconback=i;
    }
    //----------- MCLC3 --------------------------
    if(LattVar=="MCLC3" || LattVar=="mclc3") {
      found=true;
      glatt="MCLC_3";
      xrotview=75; zrotview=160;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      mu=0.25*(1+b*b/a/a);
      delta=b*c*cos(alpha)/(2*a*a);
      zzeta=mu-0.25+(1.0-b/c*cos(alpha))/(4*sin(alpha)*sin(alpha));
      eta=0.5+2*zzeta*c/b*cos(alpha);
      phi=1.0+zzeta-2*mu;
      psi=eta-2.0*delta;

      mu = 0.25+0.25*b*b/a/a;
      delta = (2*mu-0.5)*c/b*cos(alpha);
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=1.0-phi; irbz(i,2)=1.0-phi; irbz(i,3)=1.0-psi; irbzlab[i]="F"; i++;
      irbz(i,1)=phi; irbz(i,2)=phi-1.0; irbz(i,3)=psi; irbzlab[i]="F_1"; i++;
      irbz(i,1)=1.0-phi; irbz(i,2)=-phi; irbz(i,3)=1.0-psi; irbzlab[i]="F_2"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=zzeta; irbz(i,3)=eta; irbzlab[i]="H"; i++;
      irbz(i,1)=1.0-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="H_1"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="H_2"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0.5; irbzlab[i]="I"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="M"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="N"; i++;
      irbz(i,1)=0; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="N_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=mu; irbz(i,2)=mu; irbz(i,3)=delta; irbzlab[i]="Y"; i++;
      irbz(i,1)=1.0-mu; irbz(i,2)=-mu; irbz(i,3)=-delta; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=-mu; irbz(i,2)=-mu; irbz(i,3)=-delta; irbzlab[i]="Y_2"; i++;
      irbz(i,1)=mu; irbz(i,2)=mu-1.0; irbz(i,3)=delta; irbzlab[i]="Y_3"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-F-H-Z-I-X-{/Symbol G}-Z|M-{/Symbol G}-N|X-Y_1-H_1|I-F_1";
      //old: gpath="{/Symbol G}-Y-F-H-Z-I-F_1|H_1-Y_1-X-{/Symbol G}-N|M-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=2; ircon(2,i)=5;  i++;
      ircon(1,i)=1; ircon(2,i)=13;  i++;
      ircon(1,i)=1; ircon(2,i)=17;  i++;
      ircon(1,i)=2; ircon(2,i)=13;  i++;
      ircon(1,i)=5; ircon(2,i)=17;  i++;
      ircon(1,i)=8; ircon(2,i)=12;  i++;
      ircon(1,i)=8; ircon(2,i)=17;  i++;
      ircon(1,i)=8; ircon(2,i)=3;  i++;
      ircon(1,i)=6; ircon(2,i)=14;  i++;
      ircon(1,i)=1; ircon(2,i)=12;  i++;
      ircon(1,i)=1; ircon(2,i)=10;  i++;
      ircon(1,i)=1; ircon(2,i)=9;  i++;
      ircon(1,i)=14; ircon(2,i)=12;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=4; irconback(2,i)=8;  i++;
      irconback(1,i)=4; irconback(2,i)=7;  i++;
      irconback(1,i)=2; irconback(2,i)=6;  i++;
      irconback(1,i)=6; irconback(2,i)=3;  i++;
      irconback(1,i)=3; irconback(2,i)=5;  i++;
      irconback(1,i)=1; irconback(2,i)=15;  i++;
      irconback(1,i)=4; irconback(2,i)=16;  i++;
      irconback(1,i)=7; irconback(2,i)=17;  i++;
      irconback(1,i)=7; irconback(2,i)=15;  i++;
      irconback(1,i)=13; irconback(2,i)=14;  i++;
      irconback(1,i)=12; irconback(2,i)=16;  i++;
      irconback(1,i)=16; irconback(2,i)=15;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=phi; bz(i,2)=phi-1.0; bz(i,3)=psi-1.0; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=phi; bz(i,3)=-(1.0-psi); i++;
      bz(i,1)=-phi; bz(i,2)=-(phi-1.0); bz(i,3)=-psi; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=-eta; i++;	
      bz(i,1)=1.0-phi; bz(i,2)=-phi; bz(i,3)=-psi; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta-1.0; bz(i,3)=eta-1.0; i++;
      bz(i,1)=mu; bz(i,2)=mu-1.0; bz(i,3)=delta; i++;
      bz(i,1)=1.0-phi; bz(i,2)=-phi; bz(i,3)=1.0-psi; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=phi; bz(i,2)=phi-1.0; bz(i,3)=psi; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=1.0-mu; bz(i,2)=-mu; bz(i,3)=-delta; i++;
      bz(i,1)=1.0-phi; bz(i,2)=1.0-phi; bz(i,3)=1.0-psi; i++;
      bz(i,1)=mu; bz(i,2)=mu; bz(i,3)=delta; i++;
      bz(i,1)=-zzeta; bz(i,2)=-(zzeta-1.0); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=phi; bz(i,3)=psi; i++;
      bz(i,1)=-phi; bz(i,2)=-(phi-1.0); bz(i,3)=-(psi-1.0); i++;
      bz(i,1)=-(1.0-zzeta); bz(i,2)=zzeta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=-(1.0-phi); bz(i,3)=-(1.0-psi); i++;
      bz(i,1)=-mu; bz(i,2)=-mu; bz(i,3)=-delta; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      //con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=5; con(2,i)=10;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=5; con(2,i)=17;  i++;
      con(1,i)=6; con(2,i)=19;  i++;
      con(1,i)=7; con(2,i)=20;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      //con(1,i)=12; con(2,i)=13;  i++;
      //con(1,i)=13; con(2,i)=14;  i++;
      //con(1,i)=14; con(2,i)=15;  i++;
      //con(1,i)=15; con(2,i)=16;  i++;
      //con(1,i)=16; con(2,i)=17;  i++;
      //con(1,i)=16; con(2,i)=18;  i++;
      //con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=18; con(2,i)=20;  i++;
      //con(1,i)=15; con(2,i)=27;  i++;
      //con(1,i)=18; con(2,i)=27;  i++;
      //con(1,i)=13; con(2,i)=26;  i++;
      con(1,i)=22; con(2,i)=21;  i++;
      con(1,i)=21; con(2,i)=20;  i++;
      con(1,i)=27; con(2,i)=21;  i++;
      con(1,i)=22; con(2,i)=26;  Ncon=i;
      //con(1,i)=26; con(2,i)=25;  Ncon=i;
      i=1;
      conback(1,i)=10; conback(2,i)=9;  i++;
      conback(1,i)=9; conback(2,i)=8;  i++;
      conback(1,i)=8; conback(2,i)=7;  i++;
      conback(1,i)=8; conback(2,i)=23;  i++;
      conback(1,i)=9; conback(2,i)=24;  i++;
      conback(1,i)=11; conback(2,i)=24;  i++;
      conback(1,i)=24; conback(2,i)=23;  i++;
      conback(1,i)=23; conback(2,i)=22;  i++;
      conback(1,i)=24; conback(2,i)=25;  Nconback=i;
    }
    //----------- MCLC4 --------------------------
    if(LattVar=="MCLC4" || LattVar=="mclc4") {
      found=true;
      glatt="MCLC_4";
      xrotview=80; zrotview=160;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      mu=0.25*(1+b*b/a/a);
      delta=b*c*cos(alpha)/(2*a*a);
      zzeta=mu-0.25+(1.0-b/c*cos(alpha))/(4*sin(alpha)*sin(alpha));
      eta=0.5+2*zzeta*c/b*cos(alpha);
      phi=1.0+zzeta-2*mu;
      psi=eta-2.0*delta;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=1.0-phi; irbz(i,2)=1.0-phi; irbz(i,3)=1.0-psi; irbzlab[i]="F"; i++;
      irbz(i,1)=phi; irbz(i,2)=phi-1.0; irbz(i,3)=psi; irbzlab[i]=""; i++; //"F_1"; i++;
      irbz(i,1)=1.0-phi; irbz(i,2)=-phi; irbz(i,3)=1.0-psi; irbzlab[i]=""; i++; //"F_2"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=zzeta; irbz(i,3)=eta; irbzlab[i]="H"; i++;
      irbz(i,1)=1.0-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="H_1"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="H_2"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0.5; irbzlab[i]="I"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="M"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="N"; i++;
      irbz(i,1)=0; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="N_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=mu; irbz(i,2)=mu; irbz(i,3)=delta; irbzlab[i]="Y"; i++;
      irbz(i,1)=1.0-mu; irbz(i,2)=-mu; irbz(i,3)=-delta; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=-mu; irbz(i,2)=-mu; irbz(i,3)=-delta; irbzlab[i]="Y_2"; i++;
      irbz(i,1)=mu; irbz(i,2)=mu-1.0; irbz(i,3)=delta; irbzlab[i]="Y_3"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-F-H-Z-I-X-{/Symbol G}-Z|M-{/Symbol G}-N|X-Y_1-H_1";
      //old: gpath="{/Symbol G}-Y-F-H-Z-I|H_1-Y_1-X-{/Symbol G}-N|M-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=2; ircon(2,i)=5;  i++;
      ircon(1,i)=1; ircon(2,i)=13;  i++;
      ircon(1,i)=1; ircon(2,i)=17;  i++;
      ircon(1,i)=2; ircon(2,i)=13;  i++;
      ircon(1,i)=5; ircon(2,i)=17;  i++;
      ircon(1,i)=8; ircon(2,i)=12;  i++;
      ircon(1,i)=8; ircon(2,i)=17;  i++;
      ircon(1,i)=6; ircon(2,i)=14;  i++;
      ircon(1,i)=1; ircon(2,i)=12;  i++;
      ircon(1,i)=1; ircon(2,i)=10;  i++;
      ircon(1,i)=1; ircon(2,i)=9;  i++;
      ircon(1,i)=14; ircon(2,i)=12;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=2; irconback(2,i)=6;  i++;
      irconback(1,i)=6; irconback(2,i)=8;  i++;
      irconback(1,i)=8; irconback(2,i)=5;  i++;
      irconback(1,i)=8; irconback(2,i)=7;  i++;
      irconback(1,i)=1; irconback(2,i)=15;  i++;
      irconback(1,i)=8; irconback(2,i)=16;  i++;
      irconback(1,i)=7; irconback(2,i)=17;  i++;
      irconback(1,i)=7; irconback(2,i)=15;  i++;
      irconback(1,i)=13; irconback(2,i)=14;  i++;
      irconback(1,i)=12; irconback(2,i)=16;  i++;
      irconback(1,i)=16; irconback(2,i)=15;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=phi; bz(i,2)=phi-1.0; bz(i,3)=psi-1.0; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=phi; bz(i,3)=-(1.0-psi); i++;
      bz(i,1)=-phi; bz(i,2)=-(phi-1.0); bz(i,3)=-psi; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=-eta; i++;
      bz(i,1)=1.0-phi; bz(i,2)=-phi; bz(i,3)=-psi; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta-1.0; bz(i,3)=eta-1.0; i++;
      bz(i,1)=mu; bz(i,2)=mu-1.0; bz(i,3)=delta; i++;
      bz(i,1)=1.0-phi; bz(i,2)=-phi; bz(i,3)=1.0-psi; i++;
      bz(i,1)=0.5; bz(i,2)=-0.5; bz(i,3)=0.5; i++;
      bz(i,1)=phi; bz(i,2)=phi-1.0; bz(i,3)=psi; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=1.0-mu; bz(i,2)=-mu; bz(i,3)=-delta; i++;
      bz(i,1)=1.0-phi; bz(i,2)=1.0-phi; bz(i,3)=1.0-psi; i++;
      bz(i,1)=mu; bz(i,2)=mu; bz(i,3)=delta; i++;
      bz(i,1)=-zzeta; bz(i,2)=-(zzeta-1.0); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=phi; bz(i,3)=psi; i++;
      bz(i,1)=-phi; bz(i,2)=-(phi-1.0); bz(i,3)=-(psi-1.0); i++;
      bz(i,1)=-(1.0-zzeta); bz(i,2)=zzeta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-(1.0-phi); bz(i,2)=-(1.0-phi); bz(i,3)=-(1.0-psi); i++;
      bz(i,1)=-mu; bz(i,2)=-mu; bz(i,3)=-delta; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=eta; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      //con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=5; con(2,i)=10;  i++;
      con(1,i)=5; con(2,i)=17;  i++;
      con(1,i)=6; con(2,i)=19;  i++;
      con(1,i)=7; con(2,i)=20;  i++;
      con(1,i)=10; con(2,i)=11;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      //con(1,i)=12; con(2,i)=13;  i++;
      //con(1,i)=13; con(2,i)=14;  i++;
      //con(1,i)=14; con(2,i)=15;  i++;
      //con(1,i)=15; con(2,i)=16;  i++;
      //con(1,i)=16; con(2,i)=17;  i++;
      //con(1,i)=16; con(2,i)=18;  i++;
      //con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=18; con(2,i)=20;  i++;
      //con(1,i)=18; con(2,i)=27;  i++;
      //con(1,i)=13; con(2,i)=26;  i++;
      //con(1,i)=15; con(2,i)=27;  i++;
      con(1,i)=22; con(2,i)=21;  i++;
      con(1,i)=21; con(2,i)=20;  i++;
      con(1,i)=27; con(2,i)=21;  i++;
      con(1,i)=22; con(2,i)=26;  Ncon=i;
      //con(1,i)=26; con(2,i)=25;  Ncon=i;
      i=1;
      conback(1,i)=9; conback(2,i)=8;  i++;
      conback(1,i)=8; conback(2,i)=7;  i++;
      conback(1,i)=9; conback(2,i)=10;  i++;
      conback(1,i)=8; conback(2,i)=23;  i++;
      conback(1,i)=9; conback(2,i)=24;  i++;
      conback(1,i)=11; conback(2,i)=24;  i++;
      conback(1,i)=24; conback(2,i)=23;  i++;
      conback(1,i)=23; conback(2,i)=22;  i++;
      conback(1,i)=24; conback(2,i)=25;  Nconback=i;
    }
    //----------- MCLC5 --------------------------
    if(LattVar=="MCLC5" || LattVar=="mclc5") {
      found=true;
      glatt="MCLC_5";
      xrotview=80; zrotview=155;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
      zzeta=0.25*(b*b/a/a+(1.0-b/c*cos(alpha))/sin(alpha)/sin(alpha));
      eta=0.5+2*zzeta*c/b*cos(alpha);
      mu=0.5*eta+0.25*b*b/a/a-0.5*b*c/a/a*cos(alpha);
      nu=2.0*mu-zzeta;
      omega=c*(4*nu-1.0-b*b/a/a*sin(alpha)*sin(alpha))/(2*b*cos(alpha));
      delta=zzeta*c/b*cos(alpha)+0.5*omega-0.25;
      rho=1.0-zzeta*a*a/b/b;
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=nu; irbz(i,2)=nu; irbz(i,3)=omega; irbzlab[i]="F"; i++;
      irbz(i,1)=1.0-nu; irbz(i,2)=1.0-nu; irbz(i,3)=1.0-omega; irbzlab[i]="F_1"; i++;
      irbz(i,1)=nu; irbz(i,2)=nu-1.0; irbz(i,3)=omega; irbzlab[i]="F_2"; i++;
      irbz(i,1)=zzeta; irbz(i,2)=zzeta; irbz(i,3)=eta; irbzlab[i]="H"; i++;
      irbz(i,1)=1.0-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="H_1"; i++;
      irbz(i,1)=-zzeta; irbz(i,2)=-zzeta; irbz(i,3)=1.0-eta; irbzlab[i]="H_2"; i++;
      irbz(i,1)=rho; irbz(i,2)=1.0-rho; irbz(i,3)=0.5; irbzlab[i]="I"; i++;
      irbz(i,1)=1.0-rho; irbz(i,2)=rho-1.0; irbz(i,3)=0.5; irbzlab[i]="I_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="L"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="M"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="N"; i++;
      irbz(i,1)=0; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="N_1"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=mu; irbz(i,2)=mu; irbz(i,3)=delta; irbzlab[i]="Y"; i++;
      irbz(i,1)=1.0-mu; irbz(i,2)=-mu; irbz(i,3)=-delta; irbzlab[i]="Y_1"; i++;
      irbz(i,1)=-mu; irbz(i,2)=-mu; irbz(i,3)=-delta; irbzlab[i]="Y_2"; i++;
      irbz(i,1)=mu; irbz(i,2)=mu-1.0; irbz(i,3)=delta; irbzlab[i]="Y_3"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="{/Symbol G}-Y-F-L-I|I_1-Z-{/Symbol G}-X-Y_1-H_1|H-F_1|F_2-X|M-{/Symbol G}-N|H-Z";
      //connection table for kpath points
      i=1;
      ircon(1,i)=5; ircon(2,i)=3;  i++;
      ircon(1,i)=1; ircon(2,i)=15;  i++;
      ircon(1,i)=2; ircon(2,i)=15;  i++;
      ircon(1,i)=2; ircon(2,i)=10;  i++;
      ircon(1,i)=8; ircon(2,i)=10;  i++;
      ircon(1,i)=9; ircon(2,i)=19;  i++;
      ircon(1,i)=5; ircon(2,i)=19;  i++;
      ircon(1,i)=6; ircon(2,i)=16;  i++;
      ircon(1,i)=1; ircon(2,i)=14;  i++;
      ircon(1,i)=1; ircon(2,i)=12;  i++;
      ircon(1,i)=1; ircon(2,i)=11;  i++;
      ircon(1,i)=1; ircon(2,i)=19;  i++;
      ircon(1,i)=4; ircon(2,i)=14;  i++;
      ircon(1,i)=16; ircon(2,i)=14;  Nircon=i;
      //segments in irbz that are non kpath
      i=1;
      irconback(1,i)=9; irconback(2,i)=7;  i++;
      irconback(1,i)=5; irconback(2,i)=9;  i++;
      irconback(1,i)=2; irconback(2,i)=8;  i++;
      irconback(1,i)=8; irconback(2,i)=3;  i++;
      irconback(1,i)=8; irconback(2,i)=6;  i++;
      irconback(1,i)=6; irconback(2,i)=4;  i++;
      irconback(1,i)=4; irconback(2,i)=9;  i++;
      irconback(1,i)=3; irconback(2,i)=10;  i++;
      irconback(1,i)=4; irconback(2,i)=18;  i++;
      irconback(1,i)=1; irconback(2,i)=17;  i++;
      //irconback(1,i)=9; irconback(2,i)=17;  i++;
      irconback(1,i)=7; irconback(2,i)=19;  i++;
      irconback(1,i)=7; irconback(2,i)=17;  i++;
      irconback(1,i)=14; irconback(2,i)=18;  i++;
      irconback(1,i)=15; irconback(2,i)=16;  i++;
      irconback(1,i)=18; irconback(2,i)=17;  Nirconback=i;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=1.0-rho; bz(i,2)=rho-1.0; bz(i,3)=-0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=eta-1.0; i++;
      bz(i,1)=rho-1.0; bz(i,2)=1.0-rho; bz(i,3)=-0.5; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=-eta; i++;
      bz(i,1)=nu-1.0; bz(i,2)=nu-1.0; bz(i,3)=omega-1.0; i++;
      bz(i,1)=rho-1.0; bz(i,2)=-rho; bz(i,3)=-0.5; i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta-1.0; bz(i,3)=eta-1.0; i++;
      bz(i,1)=1.0-nu; bz(i,2)=-nu; bz(i,3)=-omega; i++;
      bz(i,1)=1.0-mu; bz(i,2)=-mu; bz(i,3)=-delta; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=rho; bz(i,2)=1.0-rho; bz(i,3)=0.5; i++;
      bz(i,1)=nu; bz(i,2)=nu; bz(i,3)=omega; i++;
      bz(i,1)=mu; bz(i,2)=mu; bz(i,3)=delta; i++;
      bz(i,1)=-(rho-1.0); bz(i,2)=rho; bz(i,3)=0.5; i++;
      bz(i,1)=-zzeta; bz(i,2)=-(zzeta-1.0); bz(i,3)=-(eta-1.0); i++;
      bz(i,1)=-nu; bz(i,2)=-(nu-1.0); bz(i,3)=-omega; i++;
      bz(i,1)=-(1.0-zzeta); bz(i,2)=zzeta; bz(i,3)=-(1.0-eta); i++;
      bz(i,1)=-rho; bz(i,2)=-(1.0-rho); bz(i,3)=-0.5; i++;
      bz(i,1)=-nu; bz(i,2)=-nu; bz(i,3)=-omega; i++;
      bz(i,1)=-mu; bz(i,2)=-mu; bz(i,3)=-delta; i++;
      bz(i,1)=-zzeta; bz(i,2)=-zzeta; bz(i,3)=1.0-eta; i++;
      bz(i,1)=1.0-rho; bz(i,2)=rho-1.0; bz(i,3)=0.5; i++;
      bz(i,1)=nu; bz(i,2)=nu-1.0; bz(i,3)=omega; i++;
      bz(i,1)=mu; bz(i,2)=mu-1.0; bz(i,3)=delta; i++;
      bz(i,1)=-(nu-1.0); bz(i,2)=-(nu-1.0); bz(i,3)=-(omega-1.0); i++;
      bz(i,1)=zzeta; bz(i,2)=zzeta; bz(i,3)=eta; i++;
      bz(i,1)=-(1.0-rho); bz(i,2)=-(rho-1.0); bz(i,3)=0.5; i++;
      bz(i,1)=-(1.0-nu); bz(i,2)=nu; bz(i,3)=omega; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=5; con(2,i)=12;  i++;
      con(1,i)=6; con(2,i)=17;  i++;
      con(1,i)=7; con(2,i)=20;  i++;
      con(1,i)=11; con(2,i)=12;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      //con(1,i)=13; con(2,i)=14;  i++;
      //con(1,i)=14; con(2,i)=15;  i++;
      //con(1,i)=15; con(2,i)=16;  i++;
      //con(1,i)=16; con(2,i)=17;  i++;
      con(1,i)=16; con(2,i)=18;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=19; con(2,i)=20;  i++;
      con(1,i)=11; con(2,i)=28;  i++;
      //con(1,i)=14; con(2,i)=27;  i++;
      //con(1,i)=15; con(2,i)=29;  i++;
      con(1,i)=18; con(2,i)=29;  i++;
      con(1,i)=19; con(2,i)=32;  i++;
      //con(1,i)=28; con(2,i)=27;  i++;
      //con(1,i)=24; con(2,i)=25;  i++;
      //con(1,i)=25; con(2,i)=26;  i++;
      //con(1,i)=26; con(2,i)=27;  i++;
      con(1,i)=19; con(2,i)=32;  i++;
      //con(1,i)=26; con(2,i)=30;  i++;
      //con(1,i)=29; con(2,i)=30;  i++;
      con(1,i)=25; con(2,i)=31;  i++;
      con(1,i)=30; con(2,i)=31;  i++;
      con(1,i)=31; con(2,i)=32;  Ncon=i;
      i=1;
      conback(1,i)=5; conback(2,i)=8;  i++;
      conback(1,i)=8; conback(2,i)=7;  i++;
      conback(1,i)=8; conback(2,i)=9;  i++;
      conback(1,i)=9; conback(2,i)=10;  i++;
      conback(1,i)=10; conback(2,i)=11;  i++;
      conback(1,i)=10; conback(2,i)=23;  i++;
      conback(1,i)=9; conback(2,i)=22;  i++;
      conback(1,i)=23; conback(2,i)=24;  i++;
      conback(1,i)=23; conback(2,i)=22;  i++;
      conback(1,i)=22; conback(2,i)=21;  i++;
      conback(1,i)=21; conback(2,i)=20;  i++;
      conback(1,i)=21; conback(2,i)=32;  i++;
      conback(1,i)=23; conback(2,i)=24;  Nconback=i;
    }
    //----------- TRI1a or TRI2a --------------------------
    if(LattVar=="TRI1A" || LattVar=="TRI1a" || LattVar=="tri1a" || LattVar=="TRI2A" || LattVar=="TRI2a" || LattVar=="tri2a") {	//corey and david
      found=true;
      if(LattVar=="TRI1A" || LattVar=="TRI1a" || LattVar=="tri1a")    glatt="TRI_1_a";	//corey and david
      if(LattVar=="TRI2A" || LattVar=="TRI2a" || LattVar=="tri2a")    glatt="TRI_2_a";	//corey and david
      xrotview=75; zrotview=140;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;
    
      zzeta=1.0+(b1(3)+b2(3))*(b1(3)+b2(3)+b3(3))/b1(1)/b1(1);
      zzeta=zzeta+ (b1(2)+b2(2))/b1(1)/b1(1)*(b1(2)-b2(3)*(b2(3)+b3(3))/b2(2));
      zzeta=zzeta*0.5;
      phi=0.5+0.5*b2(3)*(b2(3)+b3(3))/b2(2)/b2(2)-zzeta*b1(2)/b2(2);
      psi=0.5+zzeta*b1(3)/b3(3)+phi*b2(3)/b3(3);
      rho=1.0+b1(2)*(b1(2)+b2(2))/b1(1)/b1(1);
      rho=rho+b1(2)/b1(1)/b1(1)*(b1(3)*(b1(3)+b3(3))/b1(2)+b2(3)*(b2(3)+b3(3))/b2(2));
      rho=rho*0.5;
      tau=0.5+0.5*b2(3)*(b2(3)+b3(3))/b2(2)/b2(2) +rho*b1(2)/b2(2);
      nu=0.5+(rho*b1(3)-tau*b2(3))/b3(3);
      theta=1.0 -b1(2)*(b1(2)+b2(2))/b1(1)/b1(1);
      theta=theta+(b1(3)*(b1(3)+b3(3))*(b1(2)+b2(2))
		   -b1(2)*(b1(3)+b2(3))*(b1(3)+b2(3)+b3(3)))/b1(1)/b1(1)/b2(2);
      theta=theta*0.5;
      lambda=0.5+b1(2)/b2(2)*(1-theta)+0.5/b2(2)/b2(2)*((b1(3)+b2(3))*(b1(3)+b2(3)+b3(3))-b1(3)*(b1(3)+b3(3)));
      omega=0.5+(theta*b1(3)+lambda*b2(3))/b3(3);
      eta=0.5*(zzeta-rho);
      mu=0.5+0.5*b2(3)*(b2(3)+b3(3))/b2(2)/b2(2)-eta*b1(2)/b2(2);
      delta=0.5+(eta*b1(3)+mu*b2(3))/b3(3);
    
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="L"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="M"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="N"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0.5; irbz(i,3)=0.5; irbzlab[i]="R"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=0; irbz(i,2)=0.5; irbz(i,3)=0; irbzlab[i]="Y"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="X-{/Symbol G}-Y|L-{/Symbol G}-Z|N-{/Symbol G}-M|R-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=6; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=7;  i++;
      ircon(1,i)=2; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=8;  i++;
      ircon(1,i)=4; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=3;  i++;
      ircon(1,i)=5; ircon(2,i)=1;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=rho; bz(i,2)=-tau; bz(i,3)=-nu; i++;
      bz(i,1)=theta; bz(i,2)=lambda; bz(i,3)=-omega; i++;
      bz(i,1)=zzeta; bz(i,2)=phi; bz(i,3)=-psi; i++;
      bz(i,1)=-rho; bz(i,2)=tau; bz(i,3)=nu-1.0; i++;
      bz(i,1)=-theta; bz(i,2)=-lambda; bz(i,3)=omega-1.0; i++;
      bz(i,1)=-zzeta; bz(i,2)=-phi; bz(i,3)=psi-1.0; i++;
      bz(i,1)=-rho; bz(i,2)=tau-1.0; bz(i,3)=nu-1.0; i++;
      bz(i,1)=zzeta; bz(i,2)=phi-1.0; bz(i,3)=-psi; i++;
      bz(i,1)=1.0-theta; bz(i,2)=-lambda; bz(i,3)=omega; i++;
      bz(i,1)=1.0-rho; bz(i,2)=tau; bz(i,3)=nu; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=1.0-phi; bz(i,3)=psi; i++;
      bz(i,1)=1.0-theta; bz(i,2)=1.0-lambda; bz(i,3)=omega; i++;
      bz(i,1)=rho; bz(i,2)=-(tau-1.0); bz(i,3)=-(nu-1.0); i++;
      bz(i,1)=-zzeta; bz(i,2)=-(phi-1.0); bz(i,3)=psi; i++;
      bz(i,1)=-(1.0-theta); bz(i,2)=lambda; bz(i,3)=-omega; i++;
      bz(i,1)=-(1.0-rho); bz(i,2)=-tau; bz(i,3)=-nu; i++;
      bz(i,1)=-(1.0-zzeta); bz(i,2)=-(1.0-phi); bz(i,3)=-psi; i++;
      bz(i,1)=-theta; bz(i,2)=-lambda; bz(i,3)=omega; i++;
      bz(i,1)=-rho; bz(i,2)=tau; bz(i,3)=nu; i++;
      bz(i,1)=zzeta; bz(i,2)=phi; bz(i,3)=-(psi-1.0); i++;
      bz(i,1)=theta; bz(i,2)=lambda; bz(i,3)=-(omega-1.0); i++;
      bz(i,1)=rho; bz(i,2)=-tau; bz(i,3)=-(nu-1.0); i++;
      bz(i,1)=-zzeta; bz(i,2)=-phi; bz(i,3)=psi; i++;
      bz(i,1)=-(1.0-theta); bz(i,2)=-(1.0-lambda); bz(i,3)=-omega; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=5; con(2,i)=12;  i++;
      con(1,i)=6; con(2,i)=15;  i++;
      con(1,i)=7; con(2,i)=16;  i++;
      con(1,i)=8; con(2,i)=19;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      con(1,i)=13; con(2,i)=14;  i++;
      con(1,i)=14; con(2,i)=15;  i++;
      con(1,i)=15; con(2,i)=16;  i++;
      con(1,i)=16; con(2,i)=17;  i++;
      con(1,i)=17; con(2,i)=18;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=18; con(2,i)=23;  i++;
      con(1,i)=17; con(2,i)=24;  i++;
      con(1,i)=14; con(2,i)=25;  i++;
      con(1,i)=13; con(2,i)=26;  i++;
      con(1,i)=22; con(2,i)=23;  i++;
      con(1,i)=23; con(2,i)=24;  i++;
      con(1,i)=24; con(2,i)=25;  i++;
      con(1,i)=25; con(2,i)=26;  i++;
      con(1,i)=26; con(2,i)=27;  i++;
      con(1,i)=27; con(2,i)=22;  i++;
      con(1,i)=22; con(2,i)=23;  Ncon=i;
      i=1;
      conback(1,i)=9; conback(2,i)=8;  i++;
      conback(1,i)=5; conback(2,i)=10;  i++;
      conback(1,i)=9; conback(2,i)=10;  i++;
      conback(1,i)=9; conback(2,i)=20;  i++;
      conback(1,i)=11; conback(2,i)=12;  i++;
      conback(1,i)=10; conback(2,i)=11;  i++;
      conback(1,i)=11; conback(2,i)=28;  i++;
      conback(1,i)=19; conback(2,i)=20;  i++;
      conback(1,i)=20; conback(2,i)=21;  i++;
      conback(1,i)=21; conback(2,i)=22;  i++;
      conback(1,i)=21; conback(2,i)=28;  i++;
      conback(1,i)=28; conback(2,i)=27;  Nconback=i;
    }
    //----------- TRI1b or TRI2b --------------------------
    if( LattVar=="TRI1B" || LattVar=="TRI1b" || LattVar=="tri1b" || LattVar=="TRI2B" || LattVar=="TRI2b" || LattVar=="tri2b") {	//corey and david
      found=true;
      if(LattVar=="TRI1B" || LattVar=="TRI1b" || LattVar=="tri1b")    glatt="TRI_1_b";	//corey and david
      if(LattVar=="TRI2B" || LattVar=="TRI2b" || LattVar=="tri2b")    glatt="TRI_2_b";	//corey and david
      xrotview=75; zrotview=55;
      azview=120; elview=10;
      b1arrow=0.25; b2arrow=0.25; b3arrow=0.25;
      b1text=0.25; b2text=0.25; b3text=0.25;

      double aa,bb,cc;
      rho= 1.0+b1(2)*(b1(2)+b2(2))/b1(1)/b1(1)-b1(3)*(b3(3)-b1(3))/b1(1)/b1(1);
      rho=rho-b2(3)*(b3(3)-b2(3))*b1(2)/b1(1)/b1(1)/b2(2);
      rho=0.5*rho;
      tau=0.5+rho*b1(2)/b2(2)-0.5*b2(3)*(b3(3)-b2(3))/b2(2)/b2(2);
      nu=0.5+(rho*b1(3)-tau*b2(3))/b3(3);
      zzeta= 1.0+b1(3)*(b3(3)-b1(3))/b1(1)/b1(1)-b2(3)*(b3(3)-b2(3))*b1(2)/b1(1)/b1(1)/b2(2);
      zzeta=zzeta+b1(2)*(b2(2)-b1(2))/b1(1)/b1(1);
      zzeta=0.5*zzeta;
      phi=0.5+(zzeta-1.0)*b1(2)/b2(2)-0.5*b2(3)*(b3(3)-b2(3))/b2(2)/b2(2);
      psi=0.5+((zzeta-1.0)*b1(3)-phi*b2(3))/b3(3);
      eta=0.5*(1.0+zzeta-rho);
      mu=0.5*(tau-phi);
      delta=0.5*(1.0+psi-nu);
      aa=b1(3)/(b1(3)*b2(2)-b1(2)*b2(3));
      bb=(b3(3)-b1(3)-b2(3))/(b2(2)*(b3(3)-b1(3)-b2(3))+b2(3)*(b1(2)+b2(2)));
      cc=0.5*aa*(b1(3)-b2(3)+(b1(2)*b1(2)+b1(1)*b1(1))/b1(3)-b2(2)*b2(2)/b2(3));
      cc=cc+0.5*bb*(b3(3)-b1(3)+b2(2)*b2(2)/b2(3)+((b1(2)+b2(2))*(b1(2)+b2(2))+b1(1)*b1(1))/(b3(3)-b1(3)-b2(3)));
      theta=cc/b1(1)/(aa/b1(3)+bb/(b3(3)-b1(3)-b2(3)));
      theta=theta/b1(1);
      lambda=0.5*(b1(3)-b2(3))-b1(1)*b1(1)*(theta-0.5)/b1(3)+0.5*b1(2)*b1(2)/b1(3)-0.5*b2(2)*b2(2)/b2(3);
      lambda=1.0+theta*b1(2)/b2(2)+lambda*aa*b2(3)/b2(2);
      omega=0.5*(b1(3)-b2(3))-b1(1)*b1(1)*(theta-0.5)/b1(3)+0.5*b1(2)*b1(2)/b1(3)-0.5*b2(2)*b2(2)/b2(3);
      omega= -0.5*b2(3)-b2(2)*(omega*aa*b2(3)+0.5*b2(2))/b2(3)+theta*b1(3)+(1.0-lambda)*b2(3);
      omega=omega/b3(3);

      //rotating klattice so that pL --> kz
      //after those theta omega lambda etc have been calculated, we now can rotate
      //this is because those parameters are expressed in formula that exploit
      //the special coordinate of the original klattice
      pL(1)=-0.5; pL(2)=0; pL(3)=0.5;
      pL=pL*klattice;
      pLtheta=acos(pL(3)/modulus(pL));//=acos(scalar_product(pL,pz)/modulus(pL))
      ptmp(1)=pL(2); ptmp(2)=-pL(1); ptmp(3)=0; //=cross(pL,pz);
      pLaxis=ptmp/modulus(ptmp);
      tmpb1=Vrotate(b1,pLaxis,pLtheta);
      tmpb2=Vrotate(b2,pLaxis,pLtheta);
      tmpb3=Vrotate(b3,pLaxis,pLtheta);
      b1=tmpb1; b2=tmpb2; b3=tmpb3;
      for(j=1;j<4;j++) {
	klattice(1,j)=b1(j);
	klattice(2,j)=b2(j);
	klattice(3,j)=b3(j);
      }
    
      i=1;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="{/Symbol G}"; i++;
      irbz(i,1)=0.5; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="L"; i++;
      irbz(i,1)=0; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="M"; i++;
      irbz(i,1)=-0.5; irbz(i,2)=-0.5; irbz(i,3)=0.5; irbzlab[i]="N"; i++;
      irbz(i,1)=0; irbz(i,2)=-0.5; irbz(i,3)=0.5; irbzlab[i]="R"; i++;
      irbz(i,1)=0; irbz(i,2)=-0.5; irbz(i,3)=0; irbzlab[i]="X"; i++;
      irbz(i,1)=0.5; irbz(i,2)=0; irbz(i,3)=0; irbzlab[i]="Y"; i++;
      irbz(i,1)=-0.5; irbz(i,2)=0; irbz(i,3)=0.5; irbzlab[i]="Z"; Nirbz=i;
      gpath="X-{/Symbol G}-Y|L-{/Symbol G}-Z|N-{/Symbol G}-M|R-{/Symbol G}";
      //connection table for kpath points
      i=1;
      ircon(1,i)=6; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=7;  i++;
      ircon(1,i)=2; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=8;  i++;
      ircon(1,i)=4; ircon(2,i)=1;  i++;
      ircon(1,i)=1; ircon(2,i)=3;  i++;
      ircon(1,i)=5; ircon(2,i)=1;  Nircon=i;
      //segments in irbz that are non kpath
      Nirconback=0;
      //-- bz --
      i=1;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0.5; bz(i,2)=0; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0.5; bz(i,3)=0; i++;
      bz(i,1)=0; bz(i,2)=0; bz(i,3)=0.5; i++;
      bz(i,1)=rho; bz(i,2)=-tau; bz(i,3)=-nu; i++;
      bz(i,1)=theta; bz(i,2)=-lambda; bz(i,3)=-omega; i++;
      bz(i,1)=zzeta; bz(i,2)=-phi; bz(i,3)=-psi; i++;
      bz(i,1)=-(rho-1.0); bz(i,2)=tau; bz(i,3)=-(1.0-nu); i++;
      bz(i,1)=-(theta-1.0); bz(i,2)=lambda; bz(i,3)=-(1.0-omega); i++;
      bz(i,1)=-(zzeta-1.0); bz(i,2)=phi; bz(i,3)=-(1.0-psi); i++;
      bz(i,1)=-rho; bz(i,2)=tau; bz(i,3)=-(1.0-nu); i++;
      bz(i,1)=zzeta-1.0; bz(i,2)=-phi; bz(i,3)=-psi; i++;
      bz(i,1)=-theta; bz(i,2)=lambda-1.0; bz(i,3)=omega; i++;
      bz(i,1)=-rho; bz(i,2)=tau-1.0; bz(i,3)=nu; i++;
      bz(i,1)=1.0-zzeta; bz(i,2)=phi-1.0; bz(i,3)=psi; i++;
      bz(i,1)=1.0-theta; bz(i,2)=lambda-1.0; bz(i,3)=omega; i++;
      bz(i,1)=rho; bz(i,2)=-tau; bz(i,3)=1.0-nu; i++;
      bz(i,1)=-(zzeta-1.0); bz(i,2)=phi; bz(i,3)=psi; i++;
      bz(i,1)=theta; bz(i,2)=-(lambda-1.0); bz(i,3)=-omega; i++;
      bz(i,1)=rho; bz(i,2)=-(tau-1.0); bz(i,3)=-nu; i++;
      bz(i,1)=-(1.0-zzeta); bz(i,2)=-(phi-1.0); bz(i,3)=-psi; i++;
      bz(i,1)=-theta; bz(i,2)=lambda; bz(i,3)=omega; i++;
      bz(i,1)=-rho; bz(i,2)=tau; bz(i,3)=nu; i++;
      bz(i,1)=zzeta-1.0; bz(i,2)=-phi; bz(i,3)=1.0-psi; i++;
      bz(i,1)=theta-1.0; bz(i,2)=-lambda; bz(i,3)=1.0-omega; i++;
      bz(i,1)=rho-1.0; bz(i,2)=-tau; bz(i,3)=1.0-nu; i++;
      bz(i,1)=-zzeta; bz(i,2)=phi; bz(i,3)=psi; i++;
      bz(i,1)=-(1.0-theta); bz(i,2)=-(lambda-1.0); bz(i,3)=-omega; Nbz=i;
      //connection index
      i=1;
      con(1,i)=5; con(2,i)=6;  i++;
      con(1,i)=6; con(2,i)=7;  i++;
      con(1,i)=7; con(2,i)=8;  i++;
      con(1,i)=5; con(2,i)=12;  i++;
      con(1,i)=6; con(2,i)=15;  i++;
      con(1,i)=7; con(2,i)=16;  i++;
      con(1,i)=8; con(2,i)=19;  i++;
      con(1,i)=12; con(2,i)=13;  i++;
      con(1,i)=13; con(2,i)=14;  i++;
      con(1,i)=14; con(2,i)=15;  i++;
      con(1,i)=15; con(2,i)=16;  i++;
      con(1,i)=16; con(2,i)=17;  i++;
      con(1,i)=17; con(2,i)=18;  i++;
      con(1,i)=18; con(2,i)=19;  i++;
      con(1,i)=18; con(2,i)=23;  i++;
      con(1,i)=17; con(2,i)=24;  i++;
      con(1,i)=14; con(2,i)=25;  i++;
      con(1,i)=13; con(2,i)=26;  i++;
      con(1,i)=22; con(2,i)=23;  i++;
      con(1,i)=23; con(2,i)=24;  i++;
      con(1,i)=24; con(2,i)=25;  i++;
      con(1,i)=25; con(2,i)=26;  i++;
      con(1,i)=26; con(2,i)=27;  i++;
      con(1,i)=27; con(2,i)=22;  i++;
      con(1,i)=22; con(2,i)=23;  Ncon=i;
      i=1;
      conback(1,i)=9; conback(2,i)=8;  i++;
      conback(1,i)=5; conback(2,i)=10;  i++;
      conback(1,i)=9; conback(2,i)=10;  i++;
      conback(1,i)=9; conback(2,i)=20;  i++;
      conback(1,i)=11; conback(2,i)=12;  i++;
      conback(1,i)=10; conback(2,i)=11;  i++;
      conback(1,i)=19; conback(2,i)=20;  i++;
      conback(1,i)=11; conback(2,i)=28;  i++;
      conback(1,i)=20; conback(2,i)=21;  i++;
      conback(1,i)=21; conback(2,i)=22;  i++;
      conback(1,i)=21; conback(2,i)=28;  i++;
      conback(1,i)=28; conback(2,i)=27;  Nconback=i;
    }
  
    if(!found) glatt="unknown";
  
    xmatrix<double> cbz,cirbz;
    cbz=bz*klattice;
    cirbz=irbz*klattice;
  
    //output
    if(mode==0) {
      cout<<glatt<<endl
	  <<gpath<<endl
	  <<bzsetf<<xrotview<<" "<<bzsetf<<zrotview<<"  xrot and zrot view in degrees (gnuplot)"<<endl
	  <<bzsetf<<azview<<" "<<bzsetf<<elview<<"  azimuth and elevation view in degrees (matlab)"<<endl
	  <<bzsetf<<b1arrow<<" "<<bzsetf<<b2arrow<<" "<<bzsetf<<b3arrow
	  <<"  lenght of recip. vec. arrow in fractions of recip vec norm"<<endl
	  <<bzsetf<<b1text<<" "<<bzsetf<<b2text<<" "<<bzsetf<<b3text
	  <<"  arrow label position in fractions of recip. vec. norm"<<endl
	  <<"Reciprocal lattice vectors:"<<endl
	  <<klattice<<endl;
      cout<<Nirbz<<" #kpoints in kpath (cartesian coordinates)"<<endl;
      for(i=1;i<Nirbz+1;i++) {
	cout<<bzsetf<<cirbz(i,1)<<" "
	    <<bzsetf<<cirbz(i,2)<<" "
	    <<bzsetf<<cirbz(i,3)<<"  "
	    <<irbzlab[i]<<endl;
      }
      cout<<"Each line represents a line segment p1 to p2, "<<endl
	  <<"with format p1x p1y p1z p2x-p1x p2y-p1y p2z-p1z."<<endl;
      cout<<Nircon<<" #segments for kpath (cartesian coordinates)"<<endl;
      for(i=1;i<Nircon+1;i++) {
	cout<<bzsetf<<cirbz(ircon(1,i),1)<<" "
	    <<bzsetf<<cirbz(ircon(1,i),2)<<" "
	    <<bzsetf<<cirbz(ircon(1,i),3)<<" "
	    <<bzsetf<<cirbz(ircon(2,i),1)-cirbz(ircon(1,i),1)<<" "
	    <<bzsetf<<cirbz(ircon(2,i),2)-cirbz(ircon(1,i),2)<<" "
	    <<bzsetf<<cirbz(ircon(2,i),3)-cirbz(ircon(1,i),3)<<endl;
      }
      cout<<Nirconback<<" #segments for irbz non kpath (cartesian coordinates)"<<endl;
      for(i=1;i<Nirconback+1;i++) {
	cout<<bzsetf<<cirbz(irconback(1,i),1)<<" "
	    <<bzsetf<<cirbz(irconback(1,i),2)<<" "
	    <<bzsetf<<cirbz(irconback(1,i),3)<<" "
	    <<bzsetf<<cirbz(irconback(2,i),1)-cirbz(irconback(1,i),1)<<" "
	    <<bzsetf<<cirbz(irconback(2,i),2)-cirbz(irconback(1,i),2)<<" "
	    <<bzsetf<<cirbz(irconback(2,i),3)-cirbz(irconback(1,i),3)<<endl;
      }
      cout<<Ncon<<" #segments for front bz (cartesian coordinates)"<<endl;
      for(i=1;i<Ncon+1;i++) {
	cout<<bzsetf<<cbz(con(1,i),1)<<" "
	    <<bzsetf<<cbz(con(1,i),2)<<" "
	    <<bzsetf<<cbz(con(1,i),3)<<" "
	    <<bzsetf<<cbz(con(2,i),1)-cbz(con(1,i),1)<<" "
	    <<bzsetf<<cbz(con(2,i),2)-cbz(con(1,i),2)<<" "
	    <<bzsetf<<cbz(con(2,i),3)-cbz(con(1,i),3)<<endl;
      }
      cout<<Nconback<<" #segments for back bz (cartesian coordinates)"<<endl;
      for(i=1;i<Nconback+1;i++) {
	cout<<bzsetf<<cbz(conback(1,i),1)<<" "
	    <<bzsetf<<cbz(conback(1,i),2)<<" "
	    <<bzsetf<<cbz(conback(1,i),3)<<" "
	    <<bzsetf<<cbz(conback(2,i),1)-cbz(conback(1,i),1)<<" "
	    <<bzsetf<<cbz(conback(2,i),2)-cbz(conback(1,i),2)<<" "
	    <<bzsetf<<cbz(conback(2,i),3)-cbz(conback(1,i),3)<<endl;
      }
      cout<<Nirbz<<" #kpoints in kpath (fractional coordinates)"<<endl;
      for(i=1;i<Nirbz+1;i++) {
	cout<<bzsetf<<irbz(i,1)<<" "
	    <<bzsetf<<irbz(i,2)<<" "
	    <<bzsetf<<irbz(i,3)<<"  "
	    <<irbzlab[i]<<endl;
      }
      cout<<"Each line represents a line segment p1 to p2, "<<endl
	  <<"with format p1(1) p1(2) p1(3) to p2(1) p2(2) p2(3)."<<endl;
      cout<<Nircon<<" #segments for kpath (fractional coordinates)"<<endl;
      for(i=1;i<Nircon+1;i++) {
	cout<<bzsetf<<irbz(ircon(1,i),1)<<" "
	    <<bzsetf<<irbz(ircon(1,i),2)<<" "
	    <<bzsetf<<irbz(ircon(1,i),3)<<"  to "
	    <<bzsetf<<irbz(ircon(2,i),1)<<" "
	    <<bzsetf<<irbz(ircon(2,i),2)<<" "
	    <<bzsetf<<irbz(ircon(2,i),3)<<endl;
      }
      cout<<Nirconback<<" #segments for irbz non kpath (fractional coordinates)"<<endl;
      for(i=1;i<Nirconback+1;i++) {
	cout<<bzsetf<<irbz(irconback(1,i),1)<<" "
	    <<bzsetf<<irbz(irconback(1,i),2)<<" "
	    <<bzsetf<<irbz(irconback(1,i),3)<<"  to "
	    <<bzsetf<<irbz(irconback(2,i),1)<<" "
	    <<bzsetf<<irbz(irconback(2,i),2)<<" "
	    <<bzsetf<<irbz(irconback(2,i),3)<<endl;
      }
      cout<<Ncon<<" #segments for front bz"<<endl;
      for(i=1;i<Ncon+1;i++) {
	cout<<bzsetf<<bz(con(1,i),1)<<" "
	    <<bzsetf<<bz(con(1,i),2)<<" "
	    <<bzsetf<<bz(con(1,i),3)<<"  to "
	    <<bzsetf<<bz(con(2,i),1)<<" "
	    <<bzsetf<<bz(con(2,i),2)<<" "
	    <<bzsetf<<bz(con(2,i),3)<<endl;
      }
      cout<<Nconback<<" #segments for back bz (fractional coordinates)"<<endl;
      for(i=1;i<Nconback+1;i++) {
	cout<<bzsetf<<bz(conback(1,i),1)<<" "
	    <<bzsetf<<bz(conback(1,i),2)<<" "
	    <<bzsetf<<bz(conback(1,i),3)<<"  to "
	    <<bzsetf<<bz(conback(2,i),1)<<" "
	    <<bzsetf<<bz(conback(2,i),2)<<" "
	    <<bzsetf<<bz(conback(2,i),3)<<endl;
      }
    }
    if(mode==1) {
      ofstream fout;
      fout.open("bzplot.dat");
      cerr<<"Generating plotbz.dat file"<<endl;
      fout<<glatt<<endl
	  <<gpath<<endl
	  <<bzsetf<<xrotview<<" "<<bzsetf<<zrotview<<"  xrot and zrot view in degrees (gnuplot)"<<endl
	  <<bzsetf<<azview<<" "<<bzsetf<<elview<<"  azimuth and elevation view in degrees (matlab)"<<endl
	  <<bzsetf<<b1arrow<<" "<<bzsetf<<b2arrow<<" "<<bzsetf<<b3arrow
	  <<"  lenght of recip. vec. arrow in fractions of recip vec norm"<<endl
	  <<bzsetf<<b1text<<" "<<bzsetf<<b2text<<" "<<bzsetf<<b3text
	  <<"  arrow label position in fractions of recip. vec. norm"<<endl
	  <<"Reciprocal lattice vectors:"<<endl
	  <<klattice<<endl;
      fout<<Nirbz<<" #kpoints in kpath (cartesian coordinates)"<<endl;
      for(i=1;i<Nirbz+1;i++) {
	fout<<bzsetf<<cirbz(i,1)<<" "
	    <<bzsetf<<cirbz(i,2)<<" "
	    <<bzsetf<<cirbz(i,3)<<"  "
	    <<irbzlab[i]<<endl;
      }
      fout<<"Each line represents a line segment p1 to p2, "<<endl
	  <<"with format p1x p1y p1z p2x-p1x p2y-p1y p2z-p1z."<<endl;
      fout<<Nircon<<" #segments for kpath (cartesian coordinates)"<<endl;
      for(i=1;i<Nircon+1;i++) {
	fout<<bzsetf<<cirbz(ircon(1,i),1)<<" "
	    <<bzsetf<<cirbz(ircon(1,i),2)<<" "
	    <<bzsetf<<cirbz(ircon(1,i),3)<<" "
	    <<bzsetf<<cirbz(ircon(2,i),1)-cirbz(ircon(1,i),1)<<" "
	    <<bzsetf<<cirbz(ircon(2,i),2)-cirbz(ircon(1,i),2)<<" "
	    <<bzsetf<<cirbz(ircon(2,i),3)-cirbz(ircon(1,i),3)<<endl;
      }
      fout<<Nirconback<<" #segments for irbz non kpath (cartesian coordinates)"<<endl;
      for(i=1;i<Nirconback+1;i++) {
	fout<<bzsetf<<cirbz(irconback(1,i),1)<<" "
	    <<bzsetf<<cirbz(irconback(1,i),2)<<" "
	    <<bzsetf<<cirbz(irconback(1,i),3)<<" "
	    <<bzsetf<<cirbz(irconback(2,i),1)-cirbz(irconback(1,i),1)<<" "
	    <<bzsetf<<cirbz(irconback(2,i),2)-cirbz(irconback(1,i),2)<<" "
	    <<bzsetf<<cirbz(irconback(2,i),3)-cirbz(irconback(1,i),3)<<endl;
      }
      fout<<Ncon<<" #segments for front bz (cartesian coordinates)"<<endl;
      for(i=1;i<Ncon+1;i++) {
	fout<<bzsetf<<cbz(con(1,i),1)<<" "
	    <<bzsetf<<cbz(con(1,i),2)<<" "
	    <<bzsetf<<cbz(con(1,i),3)<<" "
	    <<bzsetf<<cbz(con(2,i),1)-cbz(con(1,i),1)<<" "
	    <<bzsetf<<cbz(con(2,i),2)-cbz(con(1,i),2)<<" "
	    <<bzsetf<<cbz(con(2,i),3)-cbz(con(1,i),3)<<endl;
      }
      fout<<Nconback<<" #segments for back bz (cartesian coordinates)"<<endl;
      for(i=1;i<Nconback+1;i++) {
	fout<<bzsetf<<cbz(conback(1,i),1)<<" "
	    <<bzsetf<<cbz(conback(1,i),2)<<" "
	    <<bzsetf<<cbz(conback(1,i),3)<<" "
	    <<bzsetf<<cbz(conback(2,i),1)-cbz(conback(1,i),1)<<" "
	    <<bzsetf<<cbz(conback(2,i),2)-cbz(conback(1,i),2)<<" "
	    <<bzsetf<<cbz(conback(2,i),3)-cbz(conback(1,i),3)<<endl;
      }
      fout<<Nirbz<<" #kpoints in kpath (fractional coordinates)"<<endl;
      for(i=1;i<Nirbz+1;i++) {
	fout<<bzsetf<<irbz(i,1)<<" "
	    <<bzsetf<<irbz(i,2)<<" "
	    <<bzsetf<<irbz(i,3)<<"  "
	    <<irbzlab[i]<<endl;
      }
      fout<<"Each line represents a line segment p1 to p2, "<<endl
	  <<"with format p1(1) p1(2) p1(3) to p2(1) p2(2) p2(3)."<<endl;
      fout<<Nircon<<" #segments for kpath (fractional coordinates)"<<endl;
      for(i=1;i<Nircon+1;i++) {
	fout<<bzsetf<<irbz(ircon(1,i),1)<<" "
	    <<bzsetf<<irbz(ircon(1,i),2)<<" "
	    <<bzsetf<<irbz(ircon(1,i),3)<<"  to "
	    <<bzsetf<<irbz(ircon(2,i),1)<<" "
	    <<bzsetf<<irbz(ircon(2,i),2)<<" "
	    <<bzsetf<<irbz(ircon(2,i),3)<<endl;
      }
      fout<<Nirconback<<" #segments for irbz non kpath (fractional coordinates)"<<endl;
      for(i=1;i<Nirconback+1;i++) {
	fout<<bzsetf<<irbz(irconback(1,i),1)<<" "
	    <<bzsetf<<irbz(irconback(1,i),2)<<" "
	    <<bzsetf<<irbz(irconback(1,i),3)<<"  to "
	    <<bzsetf<<irbz(irconback(2,i),1)<<" "
	    <<bzsetf<<irbz(irconback(2,i),2)<<" "
	    <<bzsetf<<irbz(irconback(2,i),3)<<endl;
      }
      fout<<Ncon<<" #segments for front bz"<<endl;
      for(i=1;i<Ncon+1;i++) {
	fout<<bzsetf<<bz(con(1,i),1)<<" "
	    <<bzsetf<<bz(con(1,i),2)<<" "
	    <<bzsetf<<bz(con(1,i),3)<<"  to "
	    <<bzsetf<<bz(con(2,i),1)<<" "
	    <<bzsetf<<bz(con(2,i),2)<<" "
	    <<bzsetf<<bz(con(2,i),3)<<endl;
      }
      fout<<Nconback<<" #segments for back bz (fractional coordinates)"<<endl;
      for(i=1;i<Nconback+1;i++) {
	fout<<bzsetf<<bz(conback(1,i),1)<<" "
	    <<bzsetf<<bz(conback(1,i),2)<<" "
	    <<bzsetf<<bz(conback(1,i),3)<<"  to "
	    <<bzsetf<<bz(conback(2,i),1)<<" "
	    <<bzsetf<<bz(conback(2,i),2)<<" "
	    <<bzsetf<<bz(conback(2,i),3)<<endl;
      }
      fout.close();
    
      fout.open("plotbz.sh");
      cerr<<"Generating gnuplot script plotbz.sh"<<endl;
    
      fout<<"#!/bin/ksh"<<endl
	  <<"# A script to plot brillouin zone and kpath."<<endl
	  <<"# It needs aflow with --bzplotdata option"<<endl
	  <<"# input files: bzplot.dat"<<endl
	  <<"# output files: bzplot.eps bzplot.png"<<endl
	  <<"#"<<endl<<"# written: 2010 wahyu@alumni.duke.edu"<<endl<<endl
	  <<"datfile=bzplot.dat"<<endl
	  <<"recipfile=recip.dat #reciprocal vectors"<<endl
	  <<"kptsfile=kpts.dat #kpts"<<endl
	  <<"kpathfile=kpath.dat #kpath"<<endl
	  <<"nonkpathfile=nonkpath.dat #irbz non kpath"<<endl
	  <<"fbzfile=fbz.dat #front bz"<<endl
	  <<"bbzfile=bbz.dat #back bz"<<endl<<endl
	  <<"#-----reading datfile-----"<<endl
	  <<"#---reciprocal lattice vectors b1 b2 b3---"<<endl
	  <<"exec 0< $datfile"<<endl
	  <<"$IFS read glatt"<<endl
	  <<"$IFS read gpath"<<endl
	  <<"$IFS read xrotview zrotview srest"<<endl
	  <<"read srest"<<endl
	  <<"$IFS read Lb1 Lb2 Lb3 srest"<<endl
	  <<"$IFS read Tb1 Tb2 Tb3 srest"<<endl
	  <<"read srest"<<endl
	  <<"$IFS read b11 b12 b13"<<endl
	  <<"$IFS read b21 b22 b23"<<endl
	  <<"$IFS read b31 b32 b33"<<endl
	  <<"((hb11=b11/2.0)); ((hb12=b12/2.0)); ((hb13=b13/2.0)); "<<endl
	  <<"((hb21=b21/2.0)); ((hb22=b22/2.0)); ((hb23=b23/2.0)); "<<endl
	  <<"((hb31=b31/2.0)); ((hb32=b32/2.0)); ((hb33=b33/2.0)); "<<endl
	  <<"#--length of b1 b2 b3 arrows--"<<endl
	  <<"((db11=b11*Lb1)); ((db12=b12*Lb1)); ((db13=b13*Lb1)); "<<endl
	  <<"((db21=b21*Lb2)); ((db22=b22*Lb2)); ((db23=b23*Lb2)); "<<endl
	  <<"((db31=b31*Lb3)); ((db32=b32*Lb3)); ((db33=b33*Lb3)); "<<endl
	  <<"echo $hb11 $hb12 $hb13 $db11 $db12 $db13 > $recipfile"<<endl
	  <<"echo $hb21 $hb22 $hb23 $db21 $db22 $db23 >> $recipfile"<<endl
	  <<"echo $hb31 $hb32 $hb33 $db31 $db32 $db33 >> $recipfile"<<endl
	  <<"((tb11=hb11+b11*Tb1)); ((tb12=hb12+b12*Tb1)); ((tb13=hb13+b13*Tb1));"<<endl
	  <<"((tb21=hb21+b21*Tb2)); ((tb22=hb22+b22*Tb2)); ((tb23=hb23+b23*Tb2));"<<endl
	  <<"((tb31=hb31+b31*Tb3)); ((tb32=hb32+b32*Tb3)); ((tb33=hb33+b33*Tb3));"<<endl
	  <<""<<endl
	  <<"#---kpts---"<<endl
	  <<"$IFS read Nkpts srest"<<endl
	  <<"rm -f $kptsfile"<<endl
	  <<"for ((i=0;i<Nkpts;i++))"<<endl
	  <<"do"<<endl
	  <<"  read kptsx[$i] kptsy[$i] kptsz[$i] klab[$i]"<<endl
	  <<"  echo ${kptsx[$i]} ${kptsy[$i]} ${kptsz[$i]} ${klab[$i]} >> $kptsfile"<<endl
	  <<"done"<<endl
	  <<""<<endl
	  <<"#---kpath---"<<endl
	  <<"read srest"<<endl
	  <<"read srest"<<endl
	  <<"$IFS read N srest"<<endl
	  <<"rm -f $kpathfile"<<endl
	  <<"for ((i=0;i<N;i++))"<<endl
	  <<"do"<<endl
	  <<"  read srest"<<endl
	  <<"  echo $srest >> $kpathfile"<<endl
	  <<"done"<<endl
	  <<"#---irbz non kpath---"<<endl
	  <<"$IFS read Nnonkpath srest"<<endl
	  <<"rm -f $nonkpathfile"<<endl
	  <<"for ((i=0;i<Nnonkpath;i++))"<<endl
	  <<"do"<<endl
	  <<"  read srest"<<endl
	  <<"  echo $srest >> $nonkpathfile"<<endl
	  <<"done"<<endl
	  <<"#---front bz---"<<endl
	  <<"$IFS read N srest"<<endl
	  <<"rm -f $fbzfile"<<endl
	  <<"for ((i=0;i<N;i++))"<<endl
	  <<"do"<<endl
	  <<"  read srest"<<endl
	  <<"  echo $srest >> $fbzfile"<<endl
	  <<"done"<<endl
	  <<"#---back bz---"<<endl
	  <<"$IFS read N srest"<<endl
	  <<"rm -f $bbzfile"<<endl
	  <<"for ((i=0;i<N;i++))"<<endl
	  <<"do"<<endl
	  <<"  read srest"<<endl
	  <<"  echo $srest >> $bbzfile"<<endl
	  <<"done"<<endl
	  <<""<<endl;
      /*
	fout<<"#getting the structure name and icsd number"<<endl
	<<"name=$(pwd)"<<endl
	<<"echo $name | sed \"s/\//\n/g\" | grep '_ICSD_' > wahyutmp"<<endl
	<<"exec 0< wahyutmp"<<endl
	<<"$IFS read name"<<endl
	<<"rm -f wahyutmp"<<endl
	<<""<<endl;
      */
      fout<<"echo \"set terminal postscript eps color enhanced \\\"" << DEFAULT_GNUPLOT_EPS_FONT << "\\\" 18 \"> plotbz.gnu"<<endl;
      //fout<<"echo \"set output \\\"bz_$name.eps\\\" \">> plotbz.gnu"<<endl
      fout<<"echo \"set output \\\"bzplot.eps\\\" \">> plotbz.gnu"<<endl
	  <<"echo \"set key off \">> plotbz.gnu"<<endl
	  <<"echo \"unset border \">> plotbz.gnu"<<endl
	  <<"echo \"unset xtics \">> plotbz.gnu"<<endl
	  <<"echo \"unset ytics \">> plotbz.gnu"<<endl
	  <<"echo \"unset ztics \">> plotbz.gnu"<<endl
	  <<"echo \"set view equal xyz \">> plotbz.gnu"<<endl
	  <<"echo \"set view $xrotview,$zrotview \">> plotbz.gnu"<<endl
	  <<"echo \"set title \\\"$glatt  path: $gpath\\\" \">> plotbz.gnu"<<endl
	  <<"for ((i=0;i<Nkpts;i++))"<<endl
	  <<"do"<<endl
	  <<"  echo \"set label \\\"${klab[$i]}\\\" at ${kptsx[i]},${kptsy[i]},${kptsz[i]}\" >> plotbz.gnu"<<endl
	  <<"done"<<endl
	  <<"echo \"set label \\\"b{/Symbol _1}\\\" at $tb11,$tb12,$tb13 font \\\"" << DEFAULT_GNUPLOT_EPS_FONT << ",14\\\" \">> plotbz.gnu"<<endl
	  <<"echo \"set label \\\"b{/Symbol _2}\\\" at $tb21,$tb22,$tb23 font \\\"" << DEFAULT_GNUPLOT_EPS_FONT << ",14\\\" \">> plotbz.gnu"<<endl
	  <<"echo \"set label \\\"b{/Symbol _3}\\\" at $tb31,$tb32,$tb33 font \\\"" << DEFAULT_GNUPLOT_EPS_FONT << ",14\\\" \">> plotbz.gnu"<<endl
	  <<"echo \"splot '$recipfile' using 1:2:3:4:5:6 with vectors head filled lt 1 lw 1 lc rgb \\\"#000000\\\", \\\\\" >> plotbz.gnu"<<endl
	  <<"echo \"'$fbzfile' using 1:2:3:4:5:6 with vector nohead lt 1 lw 1 lc rgb \\\"#000000\\\", \\\\\" >> plotbz.gnu"<<endl
	  <<"echo \"'$bbzfile' using 1:2:3:4:5:6 with vector nohead lt 0 lw 1.5 lc rgb \\\"#000000\\\", \\\\\" >> plotbz.gnu"<<endl
	  <<"echo \"'$kpathfile' using 1:2:3:4:5:6 with vector nohead lt 1 lw 1 lc rgb \\\"#FF0000\\\", \\\\\" >> plotbz.gnu"<<endl
	  <<"if((Nnonkpath>0))"<<endl
	  <<"then"<<endl
	  <<"  echo \"'$nonkpathfile' using 1:2:3:4:5:6 with vector nohead lt 2 lw 1 lc rgb \\\"#FF0000\\\", \\\\\" >> plotbz.gnu"<<endl
	  <<"fi"<<endl
	  <<"echo \"'$kptsfile' using 1:2:3 with points pt 7 ps 0.8 lc rgb \\\"#FF0000\\\" \" >> plotbz.gnu"<<endl
	  <<""<<endl
	  <<XHOST.command("gnuplot")<< " plotbz.gnu"<<endl;
      //fout<<XHOST.command("convert")<<" -density 200 bz_$name.eps bz_$name.png"<<endl;
      fout<<XHOST.command("convert")<<" -density 200 bzplot.eps bzplot.png"<<endl;
      fout<<"rm -f plotbz.gnu $recipfile $kptsfile $kpathfile $nonkpathfile $fbzfile $bbzfile"<<endl;
      fout.close();
      stringstream command;
      command.clear();command.str(std::string());
      command << "ksh plotbz.sh" << endl;
      aurostd::execute(command);
      command.clear();command.str(std::string());
    }
  }
}
// ***************************************************************************

// ***************************************************************************************************
namespace LATTICE {
  string KPOINTS_Directions(xstructure str_in,double grid,bool &foundBZ) {
    xstructure str_sp,str_sc;
    bool full_sym = false; //DX 8/29/17 - Speed increase
    //DX 8/29/17 [OBSOLETE] LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc);
    LATTICE::Standard_Lattice_StructureDefault(str_in,str_sp,str_sc,full_sym); //DX 8/29/17 - Speed increase
    string lattice_type=str_sp.bravais_lattice_variation_type;
    return LATTICE::KPOINTS_Directions(lattice_type,str_sp.lattice,grid,str_sp.iomode,foundBZ); 
  }
}

// ***************************************************************************************************
namespace LATTICE {
  bool string2kpoint(string _Kstr,xvector<double>& Kx,string& Ks) {
    string Kstr(_Kstr);aurostd::StringSubst(Kstr,"!"," ");
    vector<string> tokens;
    Kx[1]=0.0;Kx[2]=0.0;Kx[3]=0.0;Ks="";
    aurostd::string2tokens(Kstr,tokens," ");
    if(tokens.size()>=1) Kx[1]=aurostd::string2utype<double>(tokens.at(0));
    if(tokens.size()>=2) Kx[2]=aurostd::string2utype<double>(tokens.at(1));
    if(tokens.size()>=3) Kx[3]=aurostd::string2utype<double>(tokens.at(2));
    if(tokens.size()>=4) Ks=tokens.at(3);
    //  for(uint i=0;i<tokens.size();i++)    cerr << tokens.at(i) << endl;		 exit(0);
    return TRUE;
  }
  uint kpoint2stream(stringstream& oss, xvector<double> Kx,string Ks,uint &nkpoint) {
    oss.setf(std::ios::fixed);
    oss.precision(4);
    nkpoint++;
    oss <<"  " << (Kx[1]>=0.0?" ":"") << Kx[1] 
	<< " " << (Kx[2]>=0.0?" ":"") << Kx[2] 
	<< " " << (Kx[3]>=0.0?" ":"") << Kx[3] 
	<< "   " << double(1.0) << "   ! " << aurostd::PaddedPOST(Ks,20," ") << " // nkpoint=" << nkpoint << endl;
    return nkpoint;
  }
  uint kpoint2stream(stringstream& oss, bool isVASP, bool isQE, string A, string B, double grid,uint &nkpoint) {
    // prepare points
    xvector<double> Ax(3),Bx(3);
    string As="",Bs="";
    string2kpoint(A,Ax,As);
    string2kpoint(B,Bx,Bs);
    // prepare grid
    uint igrid=1;
    if(grid>0) igrid=uint(round(grid));
    if(grid<0) igrid=uint(double(round(max(1.0,-modulus(Bx-Ax)/grid))));
    if(grid==0.0) igrid=1;
    xvector<double> delta(3);delta=(Bx-Ax)/double(igrid);
    if(isVASP && grid>0) { // old VASP with with grid>0
      oss << A << endl << B << endl;
      nkpoint++;nkpoint++;
    }
    if(isQE || (isVASP && grid<0)) { // works with QE with grid>0 and <0
      LATTICE::kpoint2stream(oss,Ax,As,nkpoint); 
      for(uint i=1;i<(uint) igrid;i++)
	LATTICE::kpoint2stream(oss,xvector<double>(Ax+delta*double(i)),string(As+"-"+Bs+" ("+aurostd::utype2string(i)+"/"+aurostd::utype2string(igrid)+")"),nkpoint);
      LATTICE::kpoint2stream(oss,Bx,Bs,nkpoint);
    }
    return nkpoint;
  }
} // namespace LATTICE 

// ***************************************************************************************************
namespace LATTICE {
  string KPOINTS_Directions(string lattice_type,xmatrix<double> sp, double _grid,int iomode,bool &foundBZ) {
    // This function returns a string of kpath symmetry lines of the irreducible brillouin zone
    // in vasp's KPOINTS format.
    // lattice_type: bravais_lattice_variation_type
    // sp: standard primitive lattice
    
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool isVASP=(iomode==IOVASP_AUTO || iomode==IOVASP_POSCAR || iomode==IOVASP_ABCCAR || iomode==IOVASP_WYCKCAR);
    bool isQE=(iomode==IOQE_AUTO || iomode==IOQE_GEOM);
    if(!isVASP && !isQE) {
      cerr << "WARNING: LATTICE::KPOINTS_Directions unrecognized iomode=" << iomode << ", hence taking VASP" << endl;
      isVASP=TRUE;
    }
    double grid=_grid;
    //  if(grid>-1.0 && grid<0.0) grid=-1.0;
    if(grid>=0.0 && grid<1.0) grid=1.0;
    
    if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions iomode=" << iomode << endl;
    if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions isQE=" << isQE << endl;
    if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions isVASP=" << isVASP << endl;
    if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions grid=" << grid << endl;
    
    xvector<double> cdata(6),kdata(6);
    xmatrix<double> klattice(3,3),sc(3,3);
    klattice=ReciprocalLattice(sp);
    sc=LATTICE::sp2sc(sp,lattice_type,0);//from standard primitive to standard conventional
    cdata=Getabc_angles(sc,DEGREES);//data of conventional lattice_type
    //  rdata=LATTICE::Getabc_angles_Conventional(sp,lattice_type,DEGREES); // multiply with the lattice matrices..
    kdata=Getabc_angles(klattice,DEGREES);
    
    double a=cdata[1],b=cdata[2],c=cdata[3],alpha=cdata[4];//,beta=cdata[5],gamma=cdata[6];
    //  double kalpha=kdata[4],kbeta=kdata[5],kgamma=kdata[6];
    stringstream oss;
  
    foundBZ=FALSE;
    uint nkpoint=0;
    /*
      1. TRI order: kalpha,kbeta,kgamma  > 90 (kgamma<kalpha, kgamma<kbeta)
      or kalpha,kbeta,kgamma  < 90 (kgamma>kalpha, kgamma>kbeta)
      special case when kgamma=90
      "TRI1a" kalpha>90 kbeta>90 kgamma>90
      "TRI1b" kalpha<90 kbeta<90 kgamma<90
      "TRI2a" kalpha>90 kbeta>90 kgamma=90
      "TRI2b" kalpha<90 kbeta<90 kgamma=90
      2. "MCL" unique (order b<=c)
      3. MCLC (order alpha<90)
      "MCLC1"  kgamma>90
      "MCLC2"  kgamma=90
      "MCLC3"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 < 1
      "MCLC4"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 = 1
      "MCLC5"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 > 1
      4. "ORC" unique (order a<b<c)
      5. "ORCC" unique (order a<b)
      6. ORCF (order a<b<c)
      "ORCF1" "ORCF_invb2+invc2<inva2"  for 1/a^2 > 1/b^2 + 1/c^2
      "ORCF2" "ORCF_inva2<invb2+invc2"  for 1/a^2 < 1/b^2 + 1/c^2
      "ORCF3"                           for 1/a^2 = 1/b^2 + 1/c^2
      7. "ORCI" unique (order a<b<c)
      8. "TET" unique (order a a c)
      9. BCT (order a a c)
      "BCT1" "BCT_c<a" for c<a
      "BCT2" "BCT_c>a" for c>a
      10. "RHL1" alpha<90
      "RHL2" alpha>90
      11. "HEX" unique (order 60 90 90)
      12. "CUB" unique
      13. "FCC" unique (order 60 60 60)
      14. "BCC" unique
    */

    xvector<double> b1(3),b2(3),b3(3);
    b1(1)=klattice(1,1); b1(2)=klattice(1,2); b1(3)=klattice(1,3);
    b2(1)=klattice(2,1); b2(2)=klattice(2,2); b2(3)=klattice(2,3);
    b3(1)=klattice(3,1); b3(2)=klattice(3,2); b3(3)=klattice(3,3);
    
    //**************************** CUBIC (CUB) ******************************
    if(lattice_type=="CUB") {
      foundBZ=TRUE;
      string G,M,R,X;
      G="   0.000   0.000   0.000   ! \\Gamma";
      M="   0.500   0.500   0.000   ! M";
      R="   0.500   0.500   0.500   ! R";
      X="   0.000   0.500   0.000   ! X";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (simple cubic) G-X-M-G-R-X M-R"  << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,M,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,R,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,R,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,R,grid,nkpoint);
    }
    //****************************** FCC ************************************
    if(lattice_type=="FCC") {
      foundBZ=TRUE;
      string G,K,L,U,W,X;
      G="   0.000   0.000   0.000   ! \\Gamma ";
      K="   0.375   0.375   0.750   ! K ";
      L="   0.500   0.500   0.500   ! L ";
      U="   0.625   0.250   0.625   ! U ";
      W="   0.500   0.250   0.750   ! W ";
      X="   0.500   0.000   0.500   ! X ";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (face-centered cubic) G-X-W-K-G-L-U-W-L-K U-X" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,W,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,W,K,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,K,G,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,L,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,U,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,U,W,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,W,L,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,K,grid,nkpoint);            
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,U,X,grid,nkpoint);            

      //OLD: G-X-W-K-G-L-U-X  U-W-L-K
    }
    //************************ BCC ******************************************
    if(lattice_type=="BCC") {
      foundBZ=TRUE;
      string G,H,N,P;
      G="   0.000   0.000   0.000   ! \\Gamma ";
      H="   0.500  -0.500   0.500   ! H ";
      N="   0.000   0.000   0.500   ! N ";
      P="   0.250   0.250   0.250   ! P ";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (body-centered cubic) G-H-N-G-P-H P-N"  << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,H,grid,nkpoint);
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H,N,grid,nkpoint);
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,N,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,P,grid,nkpoint);
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,P,H,grid,nkpoint);
      if(isVASP && grid>0) oss << " "  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,P,N,grid,nkpoint);
    }  
    //************************* TETRAGONAL (TET) ****************************
    if(lattice_type=="TET") {
      foundBZ=TRUE;
      string G,A,M,R,X,Z;
      G="   0.000   0.000   0.000   ! \\Gamma";
      A="   0.500   0.500   0.500   ! A";
      M="   0.500   0.500   0.000   ! M";
      R="   0.000   0.500   0.500   ! R";
      X="   0.000   0.500   0.000   ! X";
      Z="   0.000   0.000   0.500   ! Z";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (tetragonal) G-X-M-G-Z-R-A-Z  X-R  M-A" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,M,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,G,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Z,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,R,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,R,A,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,A,Z,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,R,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,A,grid,nkpoint);            
    }
    // ************************* BCT  ***************************************
    if(lattice_type=="BCT1") {
      foundBZ=TRUE;
      float eta;
      string G,M,N,P,X;
      G="   0.000   0.000   0.000   ! \\Gamma";
      M="  -0.500   0.500   0.500   ! M";
      N="   0.000   0.500   0.000   ! N";
      P="   0.250   0.250   0.250   ! P";
      X="   0.000   0.000   0.500   ! X";
      eta=0.25*(1+c*c/a/a);
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      stringstream Z,Z1;
      Z << eta << "  " << eta << "  " << -eta << "  ! Z";
      Z1 << -eta << "  " << (1-eta) << "  " << eta << "  ! Z_1";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (body-centered tetragonal c < a) G-X-M-G-Z-P-N-Z1-M X-P" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,M,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Z.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z.str(),P,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,P,N,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,N,Z1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z1.str(),M,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,P,grid,nkpoint);
    
      //OLD: G-X-M-G-Z-P-N-Z1-P-X M-Z1 Z-N
    }
    if(lattice_type=="BCT2") {
      foundBZ=TRUE;
      string G,N,P,X,Z;
      G="   0.000   0.000   0.000   ! \\Gamma";
      N="   0.000   0.500   0.000   ! N";
      P="   0.250   0.250   0.250   ! P";
      X="   0.000   0.000   0.500   ! X";
      Z="   0.500   0.500  -0.500   ! Z";
      float eta,zeta;
      eta=0.25*(1+a*a/c/c);
      zeta=0.5*a*a/c/c;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;

      stringstream Sg1,Sg,Y1,Y;
      Sg << -eta << "  " << eta << "  " << eta << "  ! \\Sigma";
      Sg1 << eta << "  " << (1-eta) << "  " << -eta << "  ! \\Sigma_1";
      Y << -zeta << "  " << zeta << "   0.500   ! Y";
      Y1 << "   0.500   0.500  " << -zeta << "  ! Y_1";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (body-centered tetragonal a < c) G-X-Y-Sg-G-Z-Sg1-N-P-Y1-Z X-P" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,Y.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y.str(),Sg.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Sg.str(),G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,Sg1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Sg1.str(),N,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,N,P,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,P,Y1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y1.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,P,grid,nkpoint);
    
      //OLD: G-X-Y-Sg-G-Z-Y1-Sg1-N-P-Y X-P-Y1 Z-Sg1 Sg-N
    }
    //*********************** ORTHORHOMBIC (ORC) ****************************
    if(lattice_type=="ORC") {
      foundBZ=TRUE;
      string G,R,S,T,U,X,Y,Z;
      G="   0.000   0.000   0.000   ! \\Gamma";
      R="   0.500   0.500   0.500   ! R";
      S="   0.500   0.500   0.000   ! S";
      T="   0.000   0.500   0.500   ! T";
      U="   0.500   0.000   0.500   ! U";
      X="   0.500   0.000   0.000   ! X";
      Y="   0.000   0.500   0.000   ! Y";
      Z="   0.000   0.000   0.500   ! Z";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (orthorhombic) G-X-S-Y-G-Z-U-R-T-Z Y-T U-X S-R" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,S,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,S,Y,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y,G,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Z,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,U,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,U,R,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,R,T,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,T,Z,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y,T,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,U,X,grid,nkpoint);            
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,S,R,grid,nkpoint);
      //OLD: G-X-S-Y-G-Z-U-R-T-Z X-U Y-T S-R
    }  

    //********************** ORCF *******************************************
    if(lattice_type=="ORCF1" || lattice_type=="ORCF3") {
      foundBZ=TRUE;
      string G,Z,Y,L,T;
      G="   0.000   0.000   0.000   ! \\Gamma";
      L="   0.500   0.500   0.500   ! L";
      T="   1.000   0.500   0.500   ! T";
      Y="   0.500   0.000   0.500   ! Y";
      Z="   0.500   0.500   0.000   ! Z";
      stringstream A,A1,X,X1;
      float zeta,eta;
      zeta = 0.25*(1+a*a/b/b-a*a/c/c);
      eta = 0.25*(1+a*a/b/b+a*a/c/c);
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;

      A << "   0.500  " << 0.5+zeta << "  " <<  zeta << "   ! A";
      A1 << "   0.500  " <<  0.5-zeta << "  " <<  1.0-zeta << "   ! A_1";
      X << "   0.000  " <<  eta << "  " <<  eta << "   ! X";
      X1 << "   1.000  " <<  1.0-eta << "  " <<  1.0-eta << "   ! X_1";
      if(lattice_type=="ORCF1") {
	if(isQE) oss << "K_POINTS  crystal ! ";
	// [FIX]  if(isVASP) oss << "KPOINTS: ";
	oss << lattice_type << " (face-centered orthorhombic 1/a^2 > 1/b^2+1/c^2)) G-Y-T-Z-G-X-A1-Y T-X1 X-A-Z L-G" << endl;
      }
      else if(lattice_type=="ORCF3") {
	if(isQE) oss << "K_POINTS  crystal ! ";
	// [FIX]  if(isVASP) oss << "KPOINTS: ";
	oss << lattice_type << " (face-centered orthorhombic 1/a^2 = 1/b^2+1/c^2)) G-Y-T-Z-G-X-A1-Y X-A-Z L-G" << endl;
      }
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y,T,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,T,Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X.str(),A1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,A1.str(),Y,grid,nkpoint);
      if(lattice_type=="ORCF1") {
	if(isVASP && grid>0) oss << " " << endl;
	LATTICE::kpoint2stream(oss,isVASP,isQE,T,X1.str(),grid,nkpoint);
      }
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE, X.str(),A.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,A.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,G,grid,nkpoint);
    
      //OLD:
      //ORCF1: G-Y-T-Z-G-X-A1-X1-A-X Y-A1 T-X1 Z-A
      //ORCF3: G-Y-T-Z-G-X-A1-X1-A-X Y-A1 Z-A
    }
    if(lattice_type=="ORCF2") {
      foundBZ=TRUE;
      string G,Z,Y,X,L;
      G="   0.000   0.000   0.000   ! \\Gamma";
      L="   0.500   0.500   0.500   ! L";
      X="   0.000   0.500   0.500   ! X";
      Y="   0.500   0.000   0.500   ! Y";
      Z="   0.500   0.500   0.000   ! Z";
      stringstream H1,D,C,H,D1,C1;
      float eta,phi,delta;
      eta = 0.25*(1+a*a/b/b-a*a/c/c);
      phi = 0.25*(1+c*c/b/b-c*c/a/a);
      delta = 0.25*(1+b*b/a/a-b*b/c/c);
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: phi     = " << phi << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: delta   = " << delta << endl;

      C << "   0.500  " <<  0.5-eta  << "  " << 1.0-eta << "   ! C";
      C1 << "   0.500   " <<  0.5+eta << "  " <<  eta << "   ! C_1";
      D <<  0.5-delta << "   0.500  " << 1.0-delta << "   ! D";
      D1 <<  0.5+delta << "   0.500   " << delta << "   ! D_1";
      H <<  1.0-phi << "  " <<  0.5-phi << "   0.500   ! H";
      H1 <<  phi << "  " <<  phi+0.5 << "  " << "   0.5   ! H_1";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (face-centered orthorhombic 1/a^2 < 1/b^2+1/c^2) G-Y-C-D-X-G-Z-D1-H-C C1-Z X-H1 H-Y L-G" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y,C.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,C.str(),D.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,D.str(),X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,D1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,D1.str(),H.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H.str(),C.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,C1.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,H1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H.str(),Y,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,G,grid,nkpoint);

      //OLD: G-Y-C-D-X-G-Z-D1-C1-Z Y-H-C X-H1-D H1-C1 D1-H
    }
    //************************ ORCI *****************************************
    if(lattice_type=="ORCI") {
      foundBZ=TRUE;
      string G,Z,S,R,T,W;
      G="   0.000   0.000   0.000   ! \\Gamma";
      R="   0.000   0.500   0.000   ! R";
      S="   0.500   0.000   0.000   ! S";
      T="   0.000   0.000   0.500   ! T";
      W="   0.250   0.250   0.250   ! W";
      Z="   0.500   0.500  -0.500   ! Z";
      stringstream X1,X,Y,Y1,L2,L,L1;
      float zeta,eta,delta,mu;
      zeta = 0.25*(1+a*a/c/c);
      eta = 0.25*(1+b*b/c/c);
      delta = 0.25*(b*b-a*a)/c/c;
      mu = 0.25*(a*a+b*b)/c/c;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: delta   = " << delta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: mu      = " << mu << endl;

      L <<  -mu << "  " <<  mu << "  " <<  0.5-delta << "   ! L";
      L1 <<  mu << "  " <<  -mu << "  " <<  0.5+delta << "   ! L_1";    
      L2 <<  0.5-delta << "  " <<  0.5+delta << "  " <<  -mu << "   ! L_2";
      X <<  -zeta << "  " <<  zeta << "  " <<  zeta << "   ! X";
      X1 <<  zeta << "  " <<  (1.0-zeta) << "  " <<  -zeta << "   ! X_1";
      Y <<  eta << "  " <<  -eta << "  " <<  eta << "   ! Y";
      Y1 <<  1.0-eta << "  " <<  eta << "  " <<  -eta << "   ! Y_1";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (body-centered orthorhombc a < b < c) G-X-L-T-W-R-X1-Z-G-Y-S-W L1-Y Y1-Z" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X.str(),L.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L.str(),T,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,T,W,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,W,R,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,R,X1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X1.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y.str(),S,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,S,W,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L1.str(),Y.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y1.str(),Z,grid,nkpoint);

      //OLD: G-X-L-T-L1-Y-G-Z-X1-L2-Y1-S-W-R-X T-W R-X1 Y1-Z Y-S
    }
    //************************ ORCC *****************************************
    if(lattice_type=="ORCC") {
      foundBZ=TRUE;
      float zeta=(a*a+b*b)/(4.0*b*b);
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;

      string G,R,S,T,Y,Z;
      G="   0.000   0.000   0.000   ! \\Gamma";
      R="   0.000   0.500   0.500   ! R";
      S="   0.000   0.500   0.000   ! S";
      T="  -0.500   0.500   0.500   ! T";
      Y="  -0.500   0.500   0.000   ! Y";
      Z="   0.000   0.000   0.500   ! Z";
      stringstream A,X,X1,A1;
      A << zeta << "  " << zeta << "  0.500  ! A";
      A1 << (-zeta) << "  " << (1.0-zeta) << "  0.500  ! A_1";
      X << zeta << "  " << zeta << "  0.000  ! X";
      X1 << (-zeta) << "  " << (1.0-zeta) << "  0.000  ! X_1";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (C-centered orthorhombic a < b) G-X-S-R-A-Z-G-Y-X1-A1-T-Y Z-T" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X.str(),S,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,S,R,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,R,A.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,A.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y,X1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X1.str(),A1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,A1.str(),T,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,T,Y,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,T,grid,nkpoint);

      //OLD: G-X-S-X1-Y-G-Z-A-R-A1-T-Z X-A X1-A1 Y-T R-S
    }
    // ***************** HEXAGONAL (HEX) ************************************
    if(lattice_type=="HEX") {
      foundBZ=TRUE;
      string G,A,H,K,L,M;
      G="   0.000   0.000   0.000   ! \\Gamma";
      A="   0.000   0.000   0.500   ! A";
      H="   0.3333  0.3333  0.500   ! H";
      K="   0.3333  0.3333  0.000   ! K";
      L="   0.500   0.000   0.500   ! L";
      M="   0.500   0.000   0.000   ! M";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (hexagonal) G-M-K-G-A-L-H-A L-M K-H" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,M,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,K,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,K,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,A,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,A,L,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,H,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H,A,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,M,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,K,H,grid,nkpoint);

      // OLD: G-M-K-G-A-L-H-A  M-L  K-H
    }
    // ************************* RHOMBOHEDRAL (RHL) *************************
    if(lattice_type=="RHL1") {
      foundBZ=TRUE;
      string G,L,L1,Z,F;
      stringstream B,P,P1,P2,B1,X,Q;
      float alpharad,ap,h2,h,eta,nu;
      alpharad=alpha*PI/180.0;
      ap=2.0*a*sin(alpharad/2.0);
      h2=(a*a*cos(alpharad/2.0)*cos(alpharad/2.0))-(ap*ap/12.0);
      h=sqrt(h2);
      eta= 5.0/6.0 - ap*ap/h/h/18.0;
      nu=0.5*(1.5-eta);
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: nu      = " << nu << endl;
 
      G="   0.000   0.000   0.000    ! \\Gamma";
      F="   0.500   0.500   0.000    ! F";
      L="   0.500   0.000   0.000    ! L";
      L1="   0.000   0.000  -0.500    ! L_1";
      Z="   0.500   0.500   0.500    ! Z";
      B << eta << "  0.500  " << (1.0-eta) << "  ! B";
      B1 << "  0.500  " << (1.0-eta) << "  " << (eta-1.0) << "  ! B_1";
      P << eta << "  " << nu << "  " << nu << "  ! P";
      P1 << (1.0-nu) << "  " << (1.0-nu) << "  " << (1.0-eta) << "  ! P_1";
      P2 << nu << "  " << nu << "  " << (eta-1.0) << "  ! P_2";
      Q << (1.0-nu) << "  " << nu << "  0.000   ! Q";
      X << nu << "  0.000  " << -nu << "  ! X";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (rhombohedral alpha < 90) G-L-B1 B-Z-G-X Q-F-P1-Z L-P" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,L,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,B1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,B.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Q.str(),F,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,F,P1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,P1.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,P.str(),grid,nkpoint);
    
      // old: G-Z-P-L-G-L1-P2-F-P1-Z-B F-Q L-B1 X-G
    }
    if(lattice_type=="RHL2") {
      foundBZ=TRUE;
      string G,L,F,Z;
      stringstream Q,Q1,P1,P;
      float alpharad,ap,h2,h,ka,kalpha,kap,kh2,kh,eta,nu;
      alpharad=alpha*PI/180.0;
      ap=2.0*a*sin(alpharad/2.0);
      h2=(a*a*cos(alpharad/2.0)*cos(alpharad/2.0))-(ap*ap/12.0);
      h=sqrt(h2);
      ka=sqrt(4.0/3/ap/ap+1.0/9/h/h);
      kalpha=acos((1.0/9/h/h-2.0/3/ap/ap)/ka/ka);
      kap=2.0*ka*sin(kalpha/2.0);
      kh2=(ka*ka*cos(kalpha/2.0)*cos(kalpha/2.0))-(kap*kap/12.0);
      kh=sqrt(kh2);
      eta=0.5*ka*ka*h/kh;
      nu=(1.5-eta)*0.5;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: nu      = " << nu << endl;
    
      G="   0.000   0.000   0.000    ! \\Gamma";
      F="   0.500  -0.500   0.000    ! F";
      L="   0.500   0.000   0.000    ! L";
      Z="   0.500  -0.500   0.500    ! Z";
      Q << eta << "  " << eta << "  " << eta << "   ! Q";
      Q1 << (1-eta) << "  " << -eta << "  " << -eta << "   ! Q_1";
      P << (1.0-nu) << "  " << -nu << "  " << (1.0-nu) << "  ! P";
      P1 << nu << "  " << (nu-1.0) << "  " << (nu-1.0) << "  ! P_1";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (rhombohedral alpha > 90) G-P-Z-Q-G-F-P1-Q1-L-Z" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,P.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,P.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,Q.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Q.str(),G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,F,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,F,P1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,P1.str(),Q1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Q1.str(),L,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,Z,grid,nkpoint);
    
      //old: G-Q-Z-P-G-P1-Q1-L-Q P-Z-Q1 P-F-P1 F-G
    }
    //************************** MONOCLINIC (MCL) ***************************
    if(lattice_type=="MCL") {
      foundBZ=TRUE;
      string G,A,C,D,D1,E,X,Y,Y1,Z;
      G="   0.000   0.000   0.000   ! \\Gamma";
      A="   0.500   0.500   0.000   ! A";
      C="   0.000   0.500   0.500   ! C";
      D="   0.500   0.000   0.500   ! D";
      D1="   0.500   0.000  -0.500   ! D_1";
      E="   0.500   0.500   0.500   ! E";
      X="   0.000   0.500   0.000   ! X";
      Y="   0.000   0.000   0.500   ! Y";
      Y1="   0.000   0.000  -0.500   ! Y_1";
      Z="   0.500   0.000   0.000   ! Z";
      stringstream H2,H1,H,M2,M1,M;
      float alpharad,eta,nu;
      alpharad=alpha*PI/180.0;
      eta=0.5*(1-b*cos(alpharad)/c)/sin(alpharad)/sin(alpharad);
      nu=0.5-eta*c*cos(alpharad)/b;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: nu      = " << nu << endl;

      H << "   0.000  " << eta << "  " << (1-nu) << "  ! H";
      H1 << "   0.000  " << (1-eta) << "  " << nu << "  ! H_1";
      H2 << "   0.000  " << eta << "  " << -nu << "  ! H_2";
      M << "   0.500  " << eta << "  " << (1-nu) << "  ! M";
      M1 << "   0.500  " << (1-eta) << "  " << nu << "  ! M_1";
      M2 << "   0.500  " << eta << "  " << -nu << "  ! M_2";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (monoclinic) G-Y-H-C-E-M1-A-X-G-Z-D-M Z-A D-Y X-H1"  << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y,H.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H.str(),C,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,C,E,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,E,M1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M1.str(),A,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,A,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,D,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,D,M.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,A,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,D,Y,grid,nkpoint);	
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,H1.str(),grid,nkpoint);
	
      //old G-Y-H-C-E-M1-A-X-H1 M-D-Z Y-D
      //oldold G-Y-H-C-H1-X-H2 Y-D-Z-D1-M2-A-M1-E-M C-E X-A
    }
    // MCLC ---------------------------------------------------------------------
    if(lattice_type=="MCLC1" || lattice_type=="MCLC2") {
      foundBZ=TRUE;
      string G,N,N1,L,M,Z,Y,Y1;
      G="   0.000   0.000   0.000   ! \\Gamma";
      L="   0.500   0.500   0.500   ! L";
      M="   0.500   0.000   0.500   ! M";
      N="   0.500   0.000   0.000   ! N";
      N1="   0.000  -0.500   0.000   ! N_1";
      Y="   0.500   0.500   0.000   ! Y";
      Y1="  -0.500  -0.500   0.000   ! Y_1";
      Z="   0.000   0.000   0.500   ! Z";
      stringstream F,F1,F2,H,H1,I,I1,J,Q,X,X1,X2;
      float alpharad,zeta,eta,psi,phi,mu,delta;
      alpharad=alpha*PI/180.0;
      zeta = 0.5/sin(alpharad)/sin(alpharad) - 0.25*b/c/tan(alpharad)/sin(alpharad);
      eta = 0.5 + 2*c*cos(alpharad)*zeta/b;
      psi = 0.75 - 0.25*a*a/b/b/sin(alpharad)/sin(alpharad);
      phi = psi + 0.25*a*a/b/c/tan(alpharad)/sin(alpharad);
      mu = psi + 0.25*(2-b/c/cos(alpharad))/tan(alpharad)/tan(alpharad);
      delta = 1.0 - 0.5*b/tan(alpharad)*(1/c/sin(alpharad)-2/b/tan(alpharad)) - mu;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: psi     = " << psi << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: phi     = " << phi << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: mu      = " << mu << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: delta   = " << delta << endl;

      F << (1.0-zeta) << "  " << (1.0-zeta) << "  " << (1.0-eta) << "   ! F";
      F1 << zeta << "  " << zeta << "  " << eta << "   ! F_1";
      F2 << (-zeta) << "  " << (-zeta) << "  " << (1.0-eta) << "   ! F_2";
      H << (1.0-delta) << "  " << (1.0-mu) << "  " << (1.0-eta) << "   ! H";
      H1 << (-delta) << "  " << (-mu) << "  " << (1.0-eta) << "   ! H_1";
      I << phi << "  " << (1.0-phi) << "   0.500   ! I";
      I1 << (1.0-phi) << "  " << (phi-1.0) << "   0.500   ! I_1";
      J << (1.0-mu) << "  " << (-delta) << "  " << (1.0-eta) << "   ! J";
      Q << mu << "  " << delta << "  " << eta << "   ! Q";
      X << (1.0-psi) << "  " << (psi-1.0) << "   0.000   ! X";
      X1 << (psi) << "  " << (1-psi) << "   0.000   ! X_1";
      X2 << (psi-1.0) << "  " << (-psi) << "   0.000   ! X_2";
    
      if(lattice_type=="MCLC1") {
	if(isQE) oss << "K_POINTS  crystal ! ";
	// [FIX]  if(isVASP) oss << "KPOINTS: ";
	oss << lattice_type << " (C-centered monoclinic kgamma > 90) G-Y-F-L-I I1-Z-G-X X1-Y M-G-N Z-F1" << endl;
      }
      else if(lattice_type=="MCLC2") {
	if(isQE) oss << "K_POINTS  crystal ! ";
	// [FIX]  if(isVASP) oss << "KPOINTS: ";
	oss << lattice_type << " (C-centered monoclinic kgamma = 90) G-Y-F-L-I I1-Z-G-M N-G Z-F1" << endl;
      }
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y,F.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,F.str(),L,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,I.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,I1.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      if(lattice_type=="MCLC1") {
	LATTICE::kpoint2stream(oss,isVASP,isQE,G,X.str(),grid,nkpoint);
	if(isVASP && grid>0) oss << " " << endl;
	LATTICE::kpoint2stream(oss,isVASP,isQE,X1.str(),Y,grid,nkpoint);
	if(isVASP && grid>0) oss << " " << endl;
	LATTICE::kpoint2stream(oss,isVASP,isQE,M,G,grid,nkpoint);
	if(isVASP && grid>0) oss << " " << endl;
	LATTICE::kpoint2stream(oss,isVASP,isQE,G,N,grid,nkpoint);
      }
      else if(lattice_type=="MCLC2") {
	LATTICE::kpoint2stream(oss,isVASP,isQE,G,M,grid,nkpoint);
	if(isVASP && grid>0) oss << " " << endl;
	LATTICE::kpoint2stream(oss,isVASP,isQE,N,G,grid,nkpoint);
      }
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,F1.str(),grid,nkpoint);
    
      //old:
      //MCLC1: G-Y-F-L-I I1-Z-F1 Y-X1 X-G-N M-G
      //MCLC2: G-Y-F-L-I I1-Z-F1 N-G-M
      //oldold:
      //MCLC1: G-Y-F-L-F1-Z-F2 X1-Y I-L I1-Z X-G
      //MCLC2: G-Y-F-L-F1-Z-F2 I-L I1-Z X-G
    }
    if(lattice_type=="MCLC3" || lattice_type=="MCLC4") {
      foundBZ=TRUE;
      string G,N,N1,M,I,X,Z;
      G="   0.000   0.000   0.000   ! \\Gamma";
      I="   0.500  -0.500   0.500   ! I";
      M="   0.500   0.000   0.500   ! M";
      N="   0.500   0.000   0.000   ! N";
      N1="   0.000  -0.500   0.000   ! N_1";
      X="   0.500  -0.500   0.000   ! X";
      Z="   0.000   0.000   0.500   ! Z";
      stringstream Y1,F,F2,H1,H2,H,F1,Y,Y3,Y2;
      float alpharad,zeta,eta,psi,phi,mu,delta;
      alpharad=alpha*PI/180.0;
      zeta = 0.25/sin(alpharad)/sin(alpharad) + 0.25*b*b/a/a - 0.25*b/c/sin(alpharad)/tan(alpharad);
      eta = 0.5 + 2*c/b*zeta*cos(alpharad);
      mu = 0.25 + 0.25*b*b/a/a;
      delta = (2*mu-0.5)*c/b*cos(alpharad);
      phi = 0.75 - 0.25*b*b/a/a - 0.25*b/tan(alpharad)*(1/c/sin(alpharad)-1/b/tan(alpharad));
      psi = 0.5 + (2*phi-1)*c/b*cos(alpharad);
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: psi     = " << psi << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: phi     = " << phi << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: mu      = " << mu << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: delta   = " << delta << endl;

      F <<  (1.0-phi) << "  " <<  (1.0-phi) << "  " <<  (1.0-psi) << "   ! F";
      F1 <<  phi << "  " <<  (phi-1.0) << "  "  << psi << "   ! F_1";
      F2 <<  (1.0-phi) << "  " <<  (-phi) << "  " <<  (1.0-psi) << "   ! F_2";
      H <<  zeta << "  " <<  zeta << "  " <<  eta << "   ! H";
      H1 <<  (1.0-zeta) << "  " <<  (-zeta) << "  " <<  (1.0-eta) << "   ! H_1";
      H2 <<  (-zeta) << "  " <<  (-zeta) << "  " <<  (1.0-eta) << "   ! H_2";
      Y << mu << "  " <<  mu << "  " <<  delta << "   ! Y";
      Y1 <<  (1.0-mu) << "  " <<  (-mu) << "  " <<  (-delta) << "   ! Y_1";
      Y2 <<  (-mu) << "  " <<  (-mu) << "  " <<  (-delta) << "   ! Y_2";
      Y3 <<  mu << "  " <<  (mu-1.0) << "  " <<  delta << "   ! Y_3";
      if(lattice_type=="MCLC3") {
	if(isQE) oss << "K_POINTS  crystal ! ";
	// [FIX]  if(isVASP) oss << "KPOINTS: ";
	oss << lattice_type << " (kgamma < 90, bcos(alfa)/c+(bsin(alfa)/a)^2<1) G-Y-F-H-Z-I-X-G-Z M-G-N X-Y1-H1 I-F1" << endl;
      }
      else if(lattice_type=="MCLC4") {
	if(isQE) oss << "K_POINTS  crystal ! ";
	// [FIX]  if(isVASP) oss << "KPOINTS: ";
	oss << lattice_type << " (kgamma < 90, bcos(alfa)/c+(bsin(alfa)/a)^2=1) G-Y-F-H-Z-I-X-G-Z M-G-N X-Y1-H1" << endl;
      }
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y.str(),F.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,F.str(),H.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,I,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,I,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,N,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,Y1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y1.str(),H1.str(),grid,nkpoint);
      if(lattice_type=="MCLC3") {
	if(isVASP && grid>0) oss << " " << endl;
	LATTICE::kpoint2stream(oss,isVASP,isQE,I,F1.str(),grid,nkpoint);
      }
    
      //old:
      //MCLC3: G-Y-F-H-Z-I-F1 H1-Y1-X-G-N M-G
      //MCLC4: G-Y-F-H-Z-I H1-Y1-X-G-N M-G
      //oldold:
      //MCLC3: G-X-Y1-H1-F1-I-F2-Y3-X Z-I Y2-H2-Z-H-F-Y-G
      //MCLC4: G-X-Y1-H1-I-F2-Y3-X Z-I Y2-H2-Z-H-F-Y-G
    }
    if(lattice_type=="MCLC5") {
      foundBZ=TRUE;
      string G,N,N1,M,L,X,Z;
      G="   0.000   0.000   0.000   ! \\Gamma";
      L="   0.500   0.500   0.500   ! L";
      M="   0.500   0.000   0.500   ! M";
      N="   0.500   0.000   0.000   ! N";
      N1="   0.000  -0.500   0.000   ! N_1";
      X="   0.500  -0.500   0.000   ! X";
      Z="   0.000   0.000   0.500   ! Z";
      stringstream Y1,S,S1,H1,H2,H,Q,I,F1,F,F2,I1,Y,Y3,Y2;
      float alpharad,zeta,eta,psi,phi,mu,delta,nu,omega,rho;
      alpharad=alpha*PI/180.0;
      zeta = 0.25/sin(alpharad)/sin(alpharad) + 0.25*b*b/a/a - 0.25*b/c/sin(alpharad)/tan(alpharad);
      eta = 0.5 + 2.0*c/b*zeta*cos(alpharad);
      mu = 0.25 + 0.25*b*b/a/a - 0.25/tan(alpharad)/tan(alpharad) + 0.25*c*cos(alpharad)/b*(1.0/sin(alpharad)/sin(alpharad)-b*b/a/a);
      phi = 0.75 - 0.25*b*b/a/a - 0.25*b/tan(alpharad)*(1.0/c/sin(alpharad)-1.0/b/tan(alpharad));
      psi = 0.5 + (2.0*phi-1.0)*c/b*cos(alpharad);
      nu = 0.25 + 0.25*(b*b/a/a+b/c/sin(alpharad)/tan(alpharad)-3/tan(alpharad)/tan(alpharad)) + 0.5*c*cos(alpharad)/b*(1/sin(alpharad)/sin(alpharad)-b*b/a/a);
      omega = c*cos(alpharad)/b*(2.0*nu-0.5) + c*sin(alpharad)*tan(alpharad)/b*(2.0*nu-0.5*b*b/a/a-0.5);
      rho = 0.75 + 0.25*a*a/b/b/sin(alpharad)*(b/c/tan(alpharad)-1/sin(alpharad));
      delta = (2.0*mu-nu)*c*cos(alpharad)/b + 0.5*omega - 0.25;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: psi     = " << psi << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: phi     = " << phi << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: mu      = " << mu << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: delta   = " << delta << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: nu      = " << nu << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: rho     = " << rho << endl;
      if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: omega   = " << omega << endl;

      F <<  nu << "  " <<  nu << "  " <<  omega << "   ! F";
      F1 << (1.0-nu) << "  " << (1.0-nu) << "  " << (1.0-omega) << "   ! F_1";
      F2 <<  nu << "  " <<  (nu-1.0) << "  " <<  omega << "   ! F_2";
      H <<  zeta << "  " <<  zeta << "  " <<  eta << "   ! H";
      H1 <<  (1.0-zeta) << "  " <<  (-zeta) << "  " <<  (1.0-eta) << "   ! H_1";
      H2 <<  (-zeta) << "  " <<  (-zeta) << "  " <<  (1.0-eta) << "   ! H_2";
      I <<  rho  << "  " << (1.0-rho) << "   0.500   ! I";
      I1 << (1.0-rho) << "  " <<  (rho-1.0) << "   0.5   ! I_1";
      Q <<  phi << "  " <<  (phi-1.0) << "  "  << psi << "   ! Q";
      S <<  (1.0-phi) << "  " <<  (1.0-phi) << "  " <<  (1.0-psi) << "   ! S";
      S1 <<  (1.0-phi) << "  " <<  (-phi) << "  " <<  (1.0-psi) << "   ! S_1";
      Y << mu << "  " <<  mu << "  " <<  delta << "   ! Y";
      Y1 <<  (1.0-mu) << "  " <<  (-mu) << "  " <<  (-delta) << "   ! Y_1";
      Y2 <<  (-mu) << "  " <<  (-mu) << "  " <<  (-delta) << "   ! Y_2";
      Y3 <<  mu << "  " <<  (mu-1.0) << "  " <<  delta << "   ! Y_3";
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (kgamma < 90, bcos(alfa)/c+(bsin(alfa)/a)^2>1) G-Y-F-L-I I1-Z-G-X-Y1-H1 H-F1 F2-X M-G-N H-Z" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,Y.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y.str(),F.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,F.str(),L,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L,I.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,I1.str(),Z,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Z,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,X,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X,Y1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,Y1.str(),H1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H.str(),F1.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,F2.str(),X,grid,nkpoint);	
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,M,G,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G,N,grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,H.str(),Z,grid,nkpoint);
	
      //old: G-Y-F-L-I I1-Z-H-F1 H1-Y1-X-G-N M-G
      //oldold: G-X-Y1-H1-F2-Y3-X Z-I1 L-I Y2-H2-Z-H-F1-L-F-Y-G
    }
    //************** TRICLINIC ********************************************
    if( lattice_type=="TRI1A" || lattice_type=="TRI1a" || lattice_type=="TRI1B" || 
        lattice_type=="TRI1b" || lattice_type=="TRI2A" || lattice_type=="TRI2a" || 
        lattice_type=="TRI2B" || lattice_type=="TRI2b") {		//corey and david
      stringstream G,Y,A,D,H1,F2,B2,Z1,E,F,B,L,X,D1,H,O,K,M1,Y1,A1,F1,B1,Z,A2,D2,E1,M,R,N;
      double zeta,phi,psi,rho,tau,nu,theta,lambda,omega,eta,mu,delta;
      double aa,bb,cc;
    
      G << "   0.000   0.000   0.000   ! \\Gamma";
      if(lattice_type=="TRI1A" || lattice_type=="TRI1a" || lattice_type=="TRI2A" || lattice_type=="TRI2a") {		//corey and david
	foundBZ=TRUE;
	zeta=1.0+(b1(3)+b2(3))*(b1(3)+b2(3)+b3(3))/b1(1)/b1(1);
	zeta=zeta+ (b1(2)+b2(2))/b1(1)/b1(1)*(b1(2)-b2(3)*(b2(3)+b3(3))/b2(2));
	zeta=zeta*0.5;
	phi=0.5+0.5*b2(3)*(b2(3)+b3(3))/b2(2)/b2(2) - zeta*b1(2)/b2(2);
	psi=0.5 + zeta*b1(3)/b3(3)+phi*b2(3)/b3(3);
	rho=1.0+b1(2)*(b1(2)+b2(2))/b1(1)/b1(1);
	rho=rho+b1(2)/b1(1)/b1(1)*(b1(3)*(b1(3)+b3(3))/b1(2)+b2(3)*(b2(3)+b3(3))/b2(2));
	rho=rho*0.5;
	tau=0.5+0.5*b2(3)*(b2(3)+b3(3))/b2(2)/b2(2) +rho*b1(2)/b2(2);
	nu=0.5 + (rho*b1(3)-tau*b2(3))/b3(3);
	theta=1.0 -b1(2)*(b1(2)+b2(2))/b1(1)/b1(1);
	theta=theta+(b1(3)*(b1(3)+b3(3))*(b1(2)+b2(2)) - b1(2)*(b1(3)+b2(3))*(b1(3)+b2(3)+b3(3)))/b1(1)/b1(1)/b2(2);
	theta=theta*0.5;
	lambda=0.5+b1(2)/b2(2)*(1-theta)+0.5/b2(2)/b2(2)*((b1(3)+b2(3))*(b1(3)+b2(3)+b3(3))-b1(3)*(b1(3)+b3(3)));
	omega=0.5+(theta*b1(3)+lambda*b2(3))/b3(3);
	eta=0.5*(zeta-rho);
	mu=0.5+0.5*b2(3)*(b2(3)+b3(3))/b2(2)/b2(2)-eta*b1(2)/b2(2);
	delta=0.5+(eta*b1(3)+mu*b2(3))/b3(3);
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: phi     = " << phi << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: psi     = " << psi << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: rho     = " << rho << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: tau     = " << tau << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: nu      = " << nu << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: theta   = " << theta << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: lambda  = " << lambda << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: omega   = " << omega << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: mu      = " << mu << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: delta   = " << delta << endl;
      
	// 10 center of faces
	L << "   0.500   0.500   0.000   ! L";
	M << "   0.000   0.500   0.500   ! M";
	M1 << "   0.000  -0.500  -0.500   ! M_1";
	N << "   0.500   0.000   0.500   ! N";
	R << "   0.500   0.500   0.500   ! R";
	X << "   0.500   0.000   0.000   ! X";
	Y << "   0.000   0.500   0.000   ! Y";
	Y1 << "   0.000  -0.500   0.000   ! Y_1";
	Z << "   0.000   0.000   0.500   ! Z";
	Z1 << "   0.000   0.000  -0.500   ! Z_1";
	// 6 edges
	A << eta << "   " << mu << "   " << -delta << "   ! A";
	A1 << eta << "   " << mu-1.0 << "   " << -delta << "   ! A_1";
	A2 << eta << "   " << mu << "   " << 1.0-delta << "   ! A_2";
	B << -eta << "   " << -mu << "   " << delta-1.0 << "   ! B";
	B1 << -eta << "   " << -mu << "   " << delta << "   ! B_1";
	B2 << -eta << "   " << 1.0-mu << "   " << delta << "   ! B_2";
	// 12 corners
	D << zeta << "   " << phi << "   " << -psi << "   ! D";
	D1 << zeta << "   " << phi-1.0 << "   " << -psi << "   ! D_1";
	D2 << zeta << "   " << phi << "   " << 1.0-psi << "   ! D_2";
	E << theta << "   " << lambda << "   " << -omega << "   ! E";
	E1 << theta << "   " << lambda << "   " << 1.0-omega << "   ! E_1";
	F << rho << "   " << -tau << "   " << -nu << "   ! F";
	F1 << rho << "   " << -tau << "   " << 1.0-nu << "   ! F_1";
	F2 << rho << "   " << 1.0-tau << "   " << 1.0-nu << "   ! F_2";
	H << 1.0-theta << "   " << -lambda << "   " << omega << "   ! H";
	H1 << 1.0-theta << "   " << 1.0-lambda << "   " << omega << "   ! H_1";
	K << 1.0-zeta << "   " << 1.0-phi << "   " << psi << "   ! K";
	O << 1.0-rho << "   " << tau << "   " << nu << "   ! O";
      }
      if(lattice_type=="TRI1B" || lattice_type=="TRI1b" || lattice_type=="TRI2B" || lattice_type=="TRI2b") {		//corey and david
	foundBZ=TRUE;
	rho= 1.0 + b1(2)*(b1(2)+b2(2))/b1(1)/b1(1) - b1(3)*(b3(3)-b1(3))/b1(1)/b1(1);
	rho=rho - b2(3)*(b3(3)-b2(3))*b1(2)/b1(1)/b1(1)/b2(2);
	rho=0.5*rho;
	tau=0.5+rho*b1(2)/b2(2)-0.5*b2(3)*(b3(3)-b2(3))/b2(2)/b2(2);
	nu=0.5+(rho*b1(3)-tau*b2(3))/b3(3);
	zeta= 1.0 + b1(3)*(b3(3)-b1(3))/b1(1)/b1(1) - b2(3)*(b3(3)-b2(3))*b1(2)/b1(1)/b1(1)/b2(2);
	zeta=zeta + b1(2)*(b2(2)-b1(2))/b1(1)/b1(1);
	zeta=0.5*zeta;
	phi=0.5 + (zeta-1.0)*b1(2)/b2(2) - 0.5*b2(3)*(b3(3)-b2(3))/b2(2)/b2(2);
	psi=0.5 + ((zeta-1.0)*b1(3)-phi*b2(3))/b3(3);
	eta=0.5*(1.0+zeta-rho);
	mu=0.5*(tau-phi);
	delta=0.5*(1.0+psi-nu);
	aa=b1(3)/(b1(3)*b2(2)-b1(2)*b2(3));
	bb=(b3(3)-b1(3)-b2(3))/(b2(2)*(b3(3)-b1(3)-b2(3))+b2(3)*(b1(2)+b2(2)));
	cc=0.5*aa*(b1(3)-b2(3)+(b1(2)*b1(2)+b1(1)*b1(1))/b1(3)-b2(2)*b2(2)/b2(3));
	cc=cc+0.5*bb*(b3(3)-b1(3)+b2(2)*b2(2)/b2(3)+((b1(2)+b2(2))*(b1(2)+b2(2))+b1(1)*b1(1))/(b3(3)-b1(3)-b2(3)));
	theta=cc/b1(1)/(aa/b1(3)+bb/(b3(3)-b1(3)-b2(3)));
	theta=theta/b1(1);
	lambda=0.5*(b1(3)-b2(3)) - b1(1)*b1(1)*(theta-0.5)/b1(3) + 0.5*b1(2)*b1(2)/b1(3) - 0.5*b2(2)*b2(2)/b2(3);
	lambda=1.0 + theta*b1(2)/b2(2) + lambda*aa*b2(3)/b2(2);
	omega=0.5*(b1(3)-b2(3)) - b1(1)*b1(1)*(theta-0.5)/b1(3) + 0.5*b1(2)*b1(2)/b1(3) - 0.5*b2(2)*b2(2)/b2(3);
	omega= -0.5*b2(3)-b2(2)*(omega*aa*b2(3)+0.5*b2(2))/b2(3)+theta*b1(3)+(1.0-lambda)*b2(3);
	omega=omega/b3(3);
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: zeta    = " << zeta << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: phi     = " << phi << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: psi     = " << psi << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: rho     = " << rho << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: tau     = " << tau << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: nu      = " << nu << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: theta   = " << theta << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: lambda  = " << lambda << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: omega   = " << omega << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: eta     = " << eta << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: mu      = " << mu << endl;
	if(LDEBUG) cerr << "LATTICE::KPOINTS_Directions: delta   = " << delta << endl;

	// 10 center of faces
	L << "   0.500  -0.500   0.000   ! L";
	M << "   0.000   0.000   0.500   ! M";
	M1 << "   0.000   0.000  -0.500   ! M_1";
	N << "  -0.500  -0.500   0.500   ! N";
	R << "   0.000  -0.500   0.500   ! R";
	X << "   0.000  -0.500   0.000   ! X";
	Y << "   0.500   0.000   0.000   ! Y";
	Y1 << "  -0.500   0.000   0.000   ! Y_1";
	Z << "  -0.500   0.000   0.500   ! Z";
	Z1 << "   0.500   0.000  -0.500   ! Z_1";
	// 6 edges
	A << eta << "   " << mu << "   " << -delta << "   ! A";
	A1 << eta-1.0 << "   " << mu << "   " << -delta << "   ! A_1";
	A2 << eta-1.0 << "   " << mu << "   " << 1.0-delta << "   ! A_2";
	B << 1.0-eta << "   " << -mu << "   " << delta-1.0 << "   ! B";
	B1 << -eta << "   " << -mu << "   " << delta << "   ! B_1";
	B2 << 1.0-eta << "   " << -mu << "   " << delta << "   ! B_2";
	// 12 corners
	D << zeta << "   " << -phi << "   " << -psi << "   ! D";
	D1 << zeta-1.0 << "   " << -phi << "   " << -psi << "   ! D_1";
	D2 << zeta-1.0 << "   " << -phi << "   " << 1.0-psi << "   ! D_2";
	E << theta << "   " << -lambda << "   " << -omega << "   ! E";
	E1 << theta-1.0 << "   " << -lambda << "   " << 1.0-omega << "   ! E_1";
	F << rho << "   " << -tau << "   " << -nu << "   ! F";
	F1 << rho-1.0 << "   " << -tau << "   " << 1.0-nu << "   ! F_1";
	F2 << rho << "   " << -tau << "   " << 1.0-nu << "   ! F_2";
	H << -theta << "   " << lambda-1.0 << "   " << omega << "   ! H";
	H1 << 1.0-theta << "   " << lambda-1.0 << "   " << omega << "   ! H_1";
	K << 1.0-zeta << "   " << phi-1.0 << "   " << psi << "   ! K";
	O << -rho << "   " << tau-1.0 << "   " << nu << "   ! O";
      }
      if(isQE) oss << "K_POINTS  crystal ! ";
      // [FIX]  if(isVASP) oss << "KPOINTS: ";
      oss << lattice_type << " (triclinic) X-G-Y L-G-Z N-G-M R-G" << endl;
      if(isVASP && grid>0) oss << round(grid) << "   ! " << round(grid) << " grids "  << endl;
      if(isVASP && grid>0) oss << "Line-mode"  << endl;
      if(isVASP && grid<0) oss << "NKPOINTS  ! number of kpoints"  << endl;
      if(isVASP) oss << "reciprocal"  << endl;
      if(isQE) oss << "NKPOINTS  ! number of kpoints"  << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,X.str(),G.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G.str(),Y.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,L.str(),G.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G.str(),Z.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,N.str(),G.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,G.str(),M.str(),grid,nkpoint);
      if(isVASP && grid>0) oss << " " << endl;
      LATTICE::kpoint2stream(oss,isVASP,isQE,R.str(),G.str(),grid,nkpoint);
    }
    if(foundBZ==FALSE) {
      cerr << "WARNING: LATTICE::KPOINTS_Directions, lattice_type=" << lattice_type << " not found in aflow_lattice.cpp " << endl;
      cout << "WARNING: LATTICE::KPOINTS_Directions, lattice_type=" << lattice_type << " not found in aflow_lattice.cpp " << endl;
      oss  << "WARNING: LATTICE::KPOINTS_Directions, lattice_type=" << lattice_type << " not found in aflow_lattice.cpp " << endl;
    }
    string output(oss.str());
    if(grid<0) aurostd::StringSubst(output,"NKPOINTS",aurostd::utype2string(nkpoint));
    return output;
  }
} // namespace LATTICE

// ***************************************************************************


#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
