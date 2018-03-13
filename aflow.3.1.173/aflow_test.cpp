// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

#include "aflow.h"
//#define  __XOPTIMIZE
//#include "aflow_array.h"

// ***************************************************************************

// bool AFLOW_PTHREADS::FLAG;
// int  AFLOW_PTHREADS::MAX_PTHREADS;
// int  AFLOW_PTHREADS::RUNNING;
// pthread_t thread[MAX_ALLOCATABLE_PTHREADS];
// int iret[MAX_ALLOCATABLE_PTHREADS];
// bool thread_busy[MAX_ALLOCATABLE_PTHREADS];

void PERFORM_TESTJ(ostream& oss) {
  // load ICSD
  xmatrix<double> A(5,5),v(5,5);
  xvector<double> d(5);
  A(1,1)=1;A(1,2)=2;A(1,3)=3;A(1,4)=4;A(1,5)=5;
  A(2,1)=2;A(2,2)=2;A(2,3)=6;A(2,4)=7;A(2,5)=8;
  A(3,1)=3;A(3,2)=6;A(3,3)=3;A(3,4)=2;A(3,5)=1;
  A(4,1)=4;A(4,2)=7;A(4,3)=2;A(4,4)=4;A(4,5)=1;
  A(5,1)=5;A(5,2)=8;A(5,3)=1;A(5,4)=1;A(5,5)=5;
  oss << A << endl;
  /*
  int jacobi(const xmatrix<utype> &ain,xvector<utype> &d,xmatrix<utype> &v) {
    // Computes all eigenvalues and eigenvectors of a real symmetric xmatrix a[1..n][1..n].
    // On output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
    // v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
    // a. The function returns the number of Jacobi rotations that were required.
    */
  aurostd::jacobi(A,d,v);


  oss << aurostd::trasp(d) << endl;
 
  oss << v << endl;

}


void PERFORM_PRX(ostream& oss) {
  vector<string> vsystem;aurostd::string2tokens("AgPd,AgPt,AuPd,CdPd,CdPt,CdRh,CoPd,CoPt,CrIr,CrOs,CrPd,CrPt,CrRh,CuPd,CuPt,CuRh,FeIr,FePd,FePt,FeRh,HfIr,HfOs,HfPd,HfPt,HfRh,HfRu,HgPd,HgPt,HgRh,IrMn,IrMo,IrNb,IrNi,IrOs,IrRe,IrRh,IrRu,IrSc,IrTa,IrTc,IrTi,IrV,IrW,IrY,IrZn,IrZr,MnOs,MnPd,MnPt,MnRh,MnRu,MoOs,MoPd,MoPt,MoRh,MoRu,NbOs,NbPd,NbPt,NbRh,NbRu,NiPd,NiPt,OsRe,OsRh,OsRu,OsSc,OsTa,OsTc,OsTi,OsV,OsW,OsY,OsZr,PdPt,PdRe,PdSc,PdTa,PdTc,PdTi,PdV,PdW,PdY,PdZn,PdZr,PtRe,PtRh,PtRu,PtSc,PtTa,PtTc,PtTi,PtV,PtW,PtY,PtZn,PtZr,ReRh,ReRu,RhRu,RhSc,RhTa,RhTc,RhTi,RhV,RhW,RhY,RhZn,RhZr,RuSc,RuTa,RuTc,RuTi,RuV,RuW,RuY,RuZn,RuZr",vsystem,",");
  
  uint figN=5;
  for(uint i=0;i<vsystem.size();i+=2) {
    oss << "\\begin{figure*} ";// << endl;
    oss << "\\includegraphics[width=0.93\\linewidth,height=116mm]{shulls/fig_" << vsystem.at(i) << ".eps} ";// << endl;
    oss << "\\includegraphics[width=0.93\\linewidth,height=116mm]{shulls/fig_" << vsystem.at(i+1) << ".eps} ";// << endl;
    oss << "\\caption{\\small Convex hulls for the systems " << vsystem.at(i) << " and " << vsystem.at(i+1) << ".} ";// << endl;
    oss << "\\label{fig" << figN << "} ";// << endl;
    oss << "\\end{figure*}" << endl;
    figN++;
  }
  oss << endl;
  oss << "\\def\\allfigures{{\\cref{";
  for(uint i=5;i<=figN;i++) oss << "fig" << i << (i<figN?",":"");
  oss << "}}}" << endl;  
}

bool isPGM(string element) {
  if(element=="Os" || element=="Ru" || element=="Ir" || element=="Rh" || element=="Pd" || element=="Pt") return TRUE;
  return FALSE;
}


void PERFORM_TEST_ALLOYS(ostream& oss) {
  // load ICSD
  vector<string> vprotos;
  // aurostd::string2vectorstring(aflowlib::PrototypesIcsdReadme(),vprotos);oss << "vprotos.size()=" << vprotos.size() << endl;
  vector<string> vprotos2; for(uint i=0;i<vprotos.size();i++) if(aurostd::substring2bool(vprotos.at(i),"(2)")) vprotos2.push_back(vprotos.at(i));oss << "vprotos2.size()=" << vprotos2.size() << endl;
  vector<string> vprotos3; for(uint i=0;i<vprotos.size();i++) if(aurostd::substring2bool(vprotos.at(i),"(3)")) vprotos3.push_back(vprotos.at(i));oss << "vprotos3.size()=" << vprotos3.size() << endl;
  vector<string> vprotos4; for(uint i=0;i<vprotos.size();i++) if(aurostd::substring2bool(vprotos.at(i),"(4)")) vprotos4.push_back(vprotos.at(i));oss << "vprotos4.size()=" << vprotos4.size() << endl;
  vector<string> vprotos5; for(uint i=0;i<vprotos.size();i++) if(aurostd::substring2bool(vprotos.at(i),"(5)")) vprotos5.push_back(vprotos.at(i));oss << "vprotos5.size()=" << vprotos5.size() << endl;
  vector<string> vprotos6; for(uint i=0;i<vprotos.size();i++) if(aurostd::substring2bool(vprotos.at(i),"(6)")) vprotos6.push_back(vprotos.at(i));oss << "vprotos6.size()=" << vprotos6.size() << endl;
  vector<string> velement;
  // aurostd::string2tokens(string("Ag,Au,Cd,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // Cd is bad
  // aurostd::string2tokens(string("Ag,Au,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // Hg is bad
  // aurostd::string2tokens(string("Ag,Al,Au,Co,Cr,Cu,Fe,Hf,Ir,La,Mg,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); //-Cd-Hg +Al+Mg
  // aurostd::string2tokens(string("Ag,Al,Au,Cd,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mg,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // +Al+Mg
  //   aurostd::string2tokens(string("As,Bi,Ba,Be,Ag,Al,Au,Cd,Co,Cr,Cu,Fe,Ga,Ge,Hf,Hg,Ir,La,Li,Mg,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // +Al+Mg
  //  aurostd::string2tokens(string("Ag,Au,Cd,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,",");
  aurostd::string2tokens(string("Ag,Al,As,Au,B,Ba,Be,Bi,Br,Ca,Cd,Cl,Co,Cr,Cu,Fe,Ga,Ge,Hf,Hg,In,Ir,K,La,Li,Mg,Mn,Mo,Na,Nb,Ni,Os,P,Pb,Pd,Pt,Re,Rh,Ru,Sb,Sc,Se,Si,Sn,Sr,Ta,Tc,Te,Ti,Tl,V,W,Y,Zn,Zr"),velement,",");
  aurostd::sort(velement);

  
 //  for(uint i=0;i<velement.size();i++) 
//     for(uint j=i+1;j<velement.size();j++) 
//       for(uint k=j+1;k<velement.size();k++) 
// 	cout << "/home/auro/work/AFLOW3/aflow --terdata " << velement.at(i) << " " << velement.at(j) << " " << velement.at(k) << endl;
//   exit(0);

  vector<string> tokens,vspecies;
  bool found=FALSE;
  
  vector<string> vicsd2;
  for(uint i=0;i<velement.size();i++) {
    for(uint j=i+1;j<velement.size();j++) {  // if(isPGM(velement[i]) || isPGM(velement[j])) 
      {
	for(uint iproto=0;iproto<vprotos2.size();iproto++) {
	  if(aurostd::substring2bool(vprotos2[iproto],velement[i]) &&
	     aurostd::substring2bool(vprotos2[iproto],velement[j]))  {
	    aurostd::string2tokens(vprotos2[iproto],tokens);
	    XATOM_SplitAlloySpecies(string(tokens.at(0)),vspecies);
	    if(vspecies.at(0)==velement[i] && vspecies.at(1)==velement[j]) {
	      found=FALSE;
	      vector<string> tokensicsd;
	      for(uint n=0;n<vicsd2.size()&&!found;n++) {
		if(aurostd::substring2bool(vicsd2.at(n),velement[i]) && 
		   aurostd::substring2bool(vicsd2.at(n),velement[j])) {
		  aurostd::string2tokens(vicsd2.at(n),tokensicsd);
		  // cerr << tokens.at(0) << ":" << tokensicsd.at(0) << " - " << tokens.at(2) << ":" << tokensicsd.at(2) << endl;
		  found=found || (tokens.at(0)==tokensicsd.at(0) && tokens.at(2)==tokensicsd.at(2)); // alloy composition and pearson 
		}
	      }
	    }
	    if(!found) {
	      //   oss << vprotos2[iproto] << endl;
	      //      oss << "aflow --missing --aflow_proto=" << velement[i] << velement[j] << "/" << tokens.at(tokens.size()-2) << ".AB" << endl;
	      vicsd2.push_back(vprotos2[iproto]);
	    }
	  }
	}
      }
    }
  }
  oss << "vicsd2.size()=" << vicsd2.size() << endl;

  vector<string> vicsd3;
  for(uint i=0;i<velement.size();i++) {
    for(uint j=i+1;j<velement.size();j++) {
      for(uint k=j+1;k<velement.size();k++) {
	// if(isPGM(velement[i]) || isPGM(velement[j]) || isPGM(velement[k]))
	  {
	  for(uint iproto=0;iproto<vprotos3.size();iproto++) {
	    if(aurostd::substring2bool(vprotos3[iproto],velement[i]) &&
	       aurostd::substring2bool(vprotos3[iproto],velement[j]) &&
	       aurostd::substring2bool(vprotos3[iproto],velement[k]))  {
	      aurostd::string2tokens(vprotos3[iproto],tokens);
	      XATOM_SplitAlloySpecies(string(tokens.at(0)),vspecies);
	      if(vspecies.at(0)==velement[i] &&
		 vspecies.at(1)==velement[j] && 
		 vspecies.at(2)==velement[k]) {
		found=FALSE;
		vector<string> tokensicsd;
		for(uint n=0;n<vicsd3.size()&&!found;n++) {
		  if(aurostd::substring2bool(vicsd3.at(n),velement[i]) &&
		     aurostd::substring2bool(vicsd3.at(n),velement[j]) &&
		     aurostd::substring2bool(vicsd3.at(n),velement[k])) {
		    aurostd::string2tokens(vicsd3.at(n),tokensicsd);
		    found=found || (tokens.at(0)==tokensicsd.at(0) && tokens.at(2)==tokensicsd.at(2)); // alloy composition and pearson 
		  }
		}
	      }
	      if(!found) {
		//	oss << vprotos3[iproto] << endl;
		oss << "aflow --missing --aflow_proto=" << velement[i] << velement[j] << velement[k] << "/" << tokens.at(tokens.size()-2) << ".ABC" << endl;
 		vicsd3.push_back(vprotos3[iproto]);
	      }
	    }
	  }
	}
      }
    }
  }
  oss << "vicsd3.size()=" << vicsd3.size() << endl;
  exit(0);

  if(0) {
  vector<string> vicsd4;
  for(uint i=0;i<velement.size();i++) {
    for(uint j=i+1;j<velement.size();j++) {
      for(uint k=j+1;k<velement.size();k++) {
	for(uint l=k+1;l<velement.size();l++) {
	  if(isPGM(velement[i]) || isPGM(velement[j]) || isPGM(velement[k]) || isPGM(velement[l])) {
	    for(uint iproto=0;iproto<vprotos4.size();iproto++) {
	      if(aurostd::substring2bool(vprotos4[iproto],velement[i]) && 
		 aurostd::substring2bool(vprotos4[iproto],velement[j]) && 
		 aurostd::substring2bool(vprotos4[iproto],velement[k]) &&
		 aurostd::substring2bool(vprotos4[iproto],velement[l]))  {
		aurostd::string2tokens(vprotos4[iproto],tokens);
		XATOM_SplitAlloySpecies(string(tokens.at(0)),vspecies);
		if(vspecies.at(0)==velement[i] && 
		   vspecies.at(1)==velement[j] && 
		   vspecies.at(2)==velement[k] && 
		   vspecies.at(3)==velement[l]) {
		  found=FALSE;
		  vector<string> tokensicsd;
		  for(uint n=0;n<vicsd4.size()&&!found;n++) {
		    if(aurostd::substring2bool(vicsd4.at(n),velement[i]) && 
		       aurostd::substring2bool(vicsd4.at(n),velement[j]) && 
		       aurostd::substring2bool(vicsd4.at(n),velement[k]) && 
		       aurostd::substring2bool(vicsd4.at(n),velement[l])) {
		      aurostd::string2tokens(vicsd4.at(n),tokensicsd);
		  found=found || (tokens.at(0)==tokensicsd.at(0) && tokens.at(2)==tokensicsd.at(2)); // alloy composition and pearson 
		    }
		  }
		}
		if(!found) {
		  oss << vprotos4[iproto] << endl;
		  vicsd4.push_back(vprotos4[iproto]);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  oss << "vicsd4.size()=" << vicsd4.size() << endl;
  }
}

void __PERFORM_TEST(ostream& oss) {
  vector<string> velement;
  aurostd::string2tokens(string("Y,Sc,Zr,Hf,Ti,Tc,Re,Os,Ru,Co,Mg,Cd,Zn,Be,Tl"),velement,",");
  aurostd::sort(velement);
  for(uint i=0;i<velement.size();i++) 
    for(uint j=i+1;j<velement.size();j++)  
      oss << velement.at(i) << velement.at(j) << ",";
}

void _PERFORM_TEST(ostream& oss) {
  vector<string> vprotos;
  // FIX  aurostd::string2vectorstring(aflowlib::PrototypesIcsdReadme(),vprotos);
  vector<string> vprotos2;
  for(uint i=0;i<vprotos.size();i++)
    if(aurostd::substring2bool(vprotos.at(i),"(2)")) {
      aurostd::StringSubst(vprotos.at(i),"\t"," ");
      // cerr << "****" << vprotos.at(i) << "****" << endl;
      vprotos2.push_back(vprotos.at(i));
    }

  oss << vprotos2.size() << endl;

  vector<string> velement;
  aurostd::string2tokens(string("Ag,Al,As,Au,B,Ba,Be,Bi,Br,Ca,Cd,Cl,Co,Cr,Cu,Fe,Ga,Ge,Hf,Hg,In,Ir,K,La,Li,Mg,Mn,Mo,Na,Nb,Ni,Os,P,Pb,Pd,Pt,Re,Rh,Ru,Sb,Sc,Se,Si,Sn,Sr,Ta,Tc,Te,Ti,Tl,V,W,Y,Zn,Zr"),velement,",");
  aurostd::sort(velement);

  vector<string> vicsd;
  for(uint i=0;i<velement.size();i++) // if(velement.at(i)=="V")
    for(uint j=i+1;j<velement.size();j++)  
      for(uint l=0;l<vprotos2.size();l++) {
	if(aurostd::substring2bool(vprotos2.at(l),velement.at(i)))
	  if(aurostd::substring2bool(vprotos2.at(l),velement.at(j))) {
	      vector<string> tokens,vspecies;
	      //	      cerr << "[" << vprotos2.at(l) << "]" << endl;
	      aurostd::string2tokens(vprotos2.at(l),tokens);
	      // oss << tokens.size() << endl;
	      XATOM_SplitAlloySpecies(string(tokens.at(0)),vspecies);
	      // oss << vspecies.size() << endl;
	      if(vspecies.at(0)==velement.at(i))
		if(vspecies.at(1)==velement.at(j)) {
		  bool found=FALSE;
		  vector<string> tokensicsd;
		  for(uint n=0;n<vicsd.size()&&!found;n++) {
		    if(aurostd::substring2bool(vicsd.at(n),velement.at(i)))
		      if(aurostd::substring2bool(vicsd.at(n),velement.at(j))) {
			aurostd::string2tokens(vicsd.at(n),tokensicsd);
			// cerr << tokens.at(0) << ":" << tokensicsd.at(0) << " - " << tokens.at(2) << ":" << tokensicsd.at(2) << endl;
			if(tokens.at(0)==tokensicsd.at(0) && tokens.at(2)==tokensicsd.at(2)) found=TRUE; // alloy composition and pearson 
		      }
		  }
		  if(!found) {
		    //	    oss << velement.at(i) << velement.at(j) << " " << vprotos2.at(l) << endl;
		    oss << vprotos2.at(l) << endl;
		    vicsd.push_back(vprotos2.at(l));
		  }
		}
	  }
      }
  oss << vicsd.size() << endl;
  exit(0);
}

void PERFORM_TEST3(ostream& oss) {
  vector<string> vprotos;
  // FIX aurostd::string2vectorstring(aflowlib::PrototypesIcsdReadme(),vprotos);
  vector<string> vprotos3;
  for(uint i=0;i<vprotos.size();i++)
    if(aurostd::substring2bool(vprotos.at(i),"(3)"))
      vprotos3.push_back(vprotos.at(i));

  oss << vprotos3.size() << endl;

  vector<string> velement;
  aurostd::string2tokens(string("Ag,Al,As,Au,B,Ba,Be,Bi,Br,Ca,Cd,Cl,Co,Cr,Cu,Fe,Ga,Ge,Hf,Hg,In,Ir,K,La,Li,Mg,Mn,Mo,Na,Nb,Ni,Os,P,Pb,Pd,Pt,Re,Rh,Ru,Sb,Sc,Se,Si,Sn,Sr,Ta,Tc,Te,Ti,Tl,V,W,Y,Zn,Zr"),velement,",");
  aurostd::sort(velement);

  vector<string> vicsd;
  for(uint i=0;i<velement.size();i++) 
    for(uint j=i+1;j<velement.size();j++)  
      for(uint k=j+1;k<velement.size();k++) {
	bool found105=FALSE;
	// the 105 list from Jesus
	if(velement.at(i)=="Al" && velement.at(j)=="Au" && velement.at(k)=="Hf") found105=TRUE;
	if(velement.at(i)=="Al" && velement.at(j)=="Ge" && velement.at(k)=="Li") found105=TRUE;
	if(velement.at(i)=="Al" && velement.at(j)=="Li" && velement.at(k)=="Si") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Co" && velement.at(k)=="Hf") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Co" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Co" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Fe" && velement.at(k)=="Nb") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Fe" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Ir" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Ir" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Nb" && velement.at(k)=="Ru") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Ni" && velement.at(k)=="Sc") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Rh" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Rh" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="As" && velement.at(j)=="Ru" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Ba" && velement.at(j)=="Bi" && velement.at(k)=="K") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Co" && velement.at(k)=="Hf") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Co" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Co" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Hf" && velement.at(k)=="Rh") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Ir" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Ni" && velement.at(k)=="Sc") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Ni" && velement.at(k)=="Y") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Pd" && velement.at(k)=="Sc") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Rh" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="Bi" && velement.at(j)=="Rh" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="B" && velement.at(j)=="Li" && velement.at(k)=="Si") found105=TRUE;
	if(velement.at(i)=="Cd" && velement.at(j)=="Na" && velement.at(k)=="P") found105=TRUE;
	if(velement.at(i)=="Cl" && velement.at(j)=="La" && velement.at(k)=="Se") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Ge" && velement.at(k)=="Nb") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Ge" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Ge" && velement.at(k)=="V") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Hf" && velement.at(k)=="Sb") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Nb" && velement.at(k)=="Si") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Nb" && velement.at(k)=="Sn") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Sb" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Sb" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Si" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Sn" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Co" && velement.at(j)=="Sn" && velement.at(k)=="V") found105=TRUE;
	if(velement.at(i)=="Fe" && velement.at(j)=="Ge" && velement.at(k)=="W") found105=TRUE;
	if(velement.at(i)=="Fe" && velement.at(j)=="Nb" && velement.at(k)=="Sb") found105=TRUE;
	if(velement.at(i)=="Fe" && velement.at(j)=="Sb" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Fe" && velement.at(j)=="Sb" && velement.at(k)=="V") found105=TRUE;
	if(velement.at(i)=="Fe" && velement.at(j)=="Te" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="Ga" && velement.at(j)=="Nb" && velement.at(k)=="Ni") found105=TRUE;
	if(velement.at(i)=="Ga" && velement.at(j)=="Pt" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Hf" && velement.at(k)=="Ni") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Ir" && velement.at(k)=="Nb") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Ir" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Ir" && velement.at(k)=="V") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Ni" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Ni" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Pd" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Pt" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="Ge" && velement.at(j)=="Pt" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Hf" && velement.at(j)=="Ir" && velement.at(k)=="Sb") found105=TRUE;
	if(velement.at(i)=="Hf" && velement.at(j)=="Ni" && velement.at(k)=="Sn") found105=TRUE;
	if(velement.at(i)=="Hf" && velement.at(j)=="Pd" && velement.at(k)=="Sn") found105=TRUE;
	if(velement.at(i)=="Ir" && velement.at(j)=="Nb" && velement.at(k)=="Sn") found105=TRUE;
	if(velement.at(i)=="Ir" && velement.at(j)=="Sn" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="La" && velement.at(j)=="Pt" && velement.at(k)=="Sb") found105=TRUE;
	if(velement.at(i)=="La" && velement.at(j)=="Rh" && velement.at(k)=="Te") found105=TRUE;
	if(velement.at(i)=="Li" && velement.at(j)=="Sb" && velement.at(k)=="Zn") found105=TRUE;
	if(velement.at(i)=="Na" && velement.at(j)=="P" && velement.at(k)=="Sr") found105=TRUE;
	if(velement.at(i)=="Na" && velement.at(j)=="Sb" && velement.at(k)=="Sr") found105=TRUE;
	if(velement.at(i)=="Nb" && velement.at(j)=="Os" && velement.at(k)=="Sb") found105=TRUE;
	if(velement.at(i)=="Nb" && velement.at(j)=="Rh" && velement.at(k)=="Sn") found105=TRUE;
	if(velement.at(i)=="Nb" && velement.at(j)=="Ru" && velement.at(k)=="Sb") found105=TRUE;
	if(velement.at(i)=="Ni" && velement.at(j)=="Pb" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Ni" && velement.at(j)=="Sn" && velement.at(k)=="Ti") found105=TRUE;
	if(velement.at(i)=="Ni" && velement.at(j)=="Sn" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Os" && velement.at(j)=="Sb" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Pb" && velement.at(j)=="Pd" && velement.at(k)=="Zr") found105=TRUE;
	if(velement.at(i)=="Rh" && velement.at(j)=="Sn" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Ru" && velement.at(j)=="Sb" && velement.at(k)=="Ta") found105=TRUE;
	if(velement.at(i)=="Ru" && velement.at(j)=="Te" && velement.at(k)=="Zr") found105=TRUE;
	if(found105) {
	  for(uint l=0;l<vprotos3.size();l++) {
	    if(!aurostd::substring2bool(vprotos3.at(l),"As1Be1Li1_ICSD_100004") && 
	       !aurostd::substring2bool(vprotos3.at(l),"Be1Li1P1_ICSD_42037 ICSD_42037"))
	      {
	      if(aurostd::substring2bool(vprotos3.at(l),velement.at(i)))
		if(aurostd::substring2bool(vprotos3.at(l),velement.at(j)))
		  if(aurostd::substring2bool(vprotos3.at(l),velement.at(k))) {
		    vector<string> tokens,vspecies;
		    aurostd::string2tokens(vprotos3.at(l),tokens);
		    // oss << tokens.size() << endl;
		    XATOM_SplitAlloySpecies(tokens.at(0),vspecies);
		    // oss << vspecies.size() << endl;
		    if(vspecies.at(0)==velement.at(i))
		      if(vspecies.at(1)==velement.at(j))
			if(vspecies.at(2)==velement.at(k)) {
			  bool found=FALSE;
			  vector<string> tokensicsd;
			  for(uint n=0;n<vicsd.size()&&!found;n++) {
			    if(aurostd::substring2bool(vicsd.at(n),velement.at(i)))
			      if(aurostd::substring2bool(vicsd.at(n),velement.at(j)))
				if(aurostd::substring2bool(vicsd.at(n),velement.at(k))) {
				  aurostd::string2tokens(vicsd.at(n),tokensicsd);
				  if(tokens.at(0)==tokensicsd.at(0) && tokens.at(1)==tokensicsd.at(1) && tokens.at(3)==tokensicsd.at(3)) found=TRUE; // alloy composition and pearson 
				}
			  }
			  if(!found) {
			    oss << velement.at(i) << velement.at(j) << velement.at(k) << " " << vprotos3.at(l) << endl;
			    vicsd.push_back(vprotos3.at(l));
			  }
			}
		  }
	    }
	  }
	}
      }
  oss << vicsd.size() << endl;
  vector<string> tokens,tokens2;
  for(uint i=0;i<vicsd.size();i++) {
    aurostd::string2tokens(vicsd.at(i),tokens);
    oss << "/home/auro/work/AFLOW3/aflow --noldau --aflow_proto=" << tokens.at(tokens.size()-3) << endl;
  }
  for(uint i=0;i<vicsd.size();i++) {
    aurostd::string2tokens(vicsd.at(i),tokens);
    // oss << "mv ./ICSD/" << tokens.at(tokens.size()-3) << "/"+_AFLOWIN_;
    aurostd::string2tokens(string(tokens.at(tokens.size()-3)),tokens,"_");
    aurostd::StringSubst(tokens.at(0),"V1","V_sv");   
    aurostd::StringSubst(tokens.at(0),"Y1","Y_sv");   
    aurostd::StringSubst(tokens.at(0),"0","");aurostd::StringSubst(tokens.at(0),"1","");aurostd::StringSubst(tokens.at(0),"2","");
    aurostd::StringSubst(tokens.at(0),"3","");aurostd::StringSubst(tokens.at(0),"4","");aurostd::StringSubst(tokens.at(0),"5","");
    aurostd::StringSubst(tokens.at(0),"6","");aurostd::StringSubst(tokens.at(0),"7","");aurostd::StringSubst(tokens.at(0),"8","");
    aurostd::StringSubst(tokens.at(0),"9","");
    
    aurostd::StringSubst(tokens.at(0),"Ba","Ba_sv");aurostd::StringSubst(tokens.at(0),"Be","Be_sv");aurostd::StringSubst(tokens.at(0),"Bi","Bi_d");
    aurostd::StringSubst(tokens.at(0),"Ca","Ca_sv");aurostd::StringSubst(tokens.at(0),"Cr","Cr_pv");aurostd::StringSubst(tokens.at(0),"Cu","Cu_pv");
    aurostd::StringSubst(tokens.at(0),"Fe","Fe_pv");aurostd::StringSubst(tokens.at(0),"Ga","Ga_h");aurostd::StringSubst(tokens.at(0),"Ge","Ge_h");
    aurostd::StringSubst(tokens.at(0),"Hf","Hf_pv");aurostd::StringSubst(tokens.at(0),"In","In_d");aurostd::StringSubst(tokens.at(0),"Li","Li_sv");
    aurostd::StringSubst(tokens.at(0),"Mg","Mg_pv");aurostd::StringSubst(tokens.at(0),"Mn","Mn_pv");aurostd::StringSubst(tokens.at(0),"Mo","Mo_pv");
    aurostd::StringSubst(tokens.at(0),"Na","Na_sv");aurostd::StringSubst(tokens.at(0),"Nb","Nb_sv");aurostd::StringSubst(tokens.at(0),"Ni","Ni_pv");
    aurostd::StringSubst(tokens.at(0),"Os","Os_pv");aurostd::StringSubst(tokens.at(0),"Pb","Pb_d");aurostd::StringSubst(tokens.at(0),"Pd","Pd_pv");
    aurostd::StringSubst(tokens.at(0),"Re","Re_pv");aurostd::StringSubst(tokens.at(0),"Rh","Rh_pv");aurostd::StringSubst(tokens.at(0),"Ru","Ru_pv");
    aurostd::StringSubst(tokens.at(0),"Sc","Sc_sv");aurostd::StringSubst(tokens.at(0),"Sr","Sr_sv");aurostd::StringSubst(tokens.at(0),"Ta","Ta_pv");
    aurostd::StringSubst(tokens.at(0),"Tc","Tc_pv");aurostd::StringSubst(tokens.at(0),"Ti","Ti_sv");aurostd::StringSubst(tokens.at(0),"Tl","Tl_d");
    aurostd::StringSubst(tokens.at(0),"Zr","Zr_sv");
    
    oss << " ./LIB3/LIB/" << tokens.at(0) << "/" << tokens.at(1) << "_" << tokens.at(2) << ".ABC" << endl;
    //  oss << "mkdir  ./LIB3/" << " && " << "mkdir  ./LIB3/LIB/" << " && " << "mkdir  ./LIB3/LIB/" << tokens.at(0) << " && " "mkdir  ./LIB3/LIB/" << tokens.at(0) << " && " << "mkdir  ./LIB3/LIB/" << tokens.at(0) << "/" << tokens.at(1) << "_" << tokens.at(2) << ".ABC" << endl;
  }
  // oss << vprotos3.size() << endl;
  exit(0);
}

//extern vector<string> vAUID,vAURL;
vector<string> vAUID_new,vAURL_new;

void PERFORM_TEST1(ostream& oss) {
  aflowlib::auid2present(); // empty=- load
  
  oss << vAUID.size() << " " << vAURL.size() << endl; 
  
  for(uint i=0;i<vAUID.size();i++) { 
    uint j;
    if((j=aflowlib::auid2present(vAUID.at(i)))) {
      oss << "found duplicate i=" << i << " j=" << j << " " << vAUID.at(i) << " " << vAURL.at(i) << " " << vAURL_new.at(j) << endl;
    } else { 
      vAUID_new.push_back(vAUID.at(i)); // survives all
      vAURL_new.push_back(vAURL.at(i)); // survives all
      //      if(!aurostd::mod((int) vAUID_new.size(),20000))	oss << vAUID_new.size() << " " << vAURL_new.size() << endl; 
    }
  }
  oss << vAUID_new.size() << " " << vAURL_new.size() << endl; 
  
}

#define NTERM 5
#define SPREAD 0.1
#define NPT 100

//    void lfit(xvector<utype> x, xvector<utype> y, xvector<utype> sig, 
//                                  xvector<utype>& a, xvector<int> ia, 
//                                  xmatrix<utype>& covar, utype& chisq, 
//                                  void (*funcs)(utype, xvector<utype>));

void pinku_funcs(float x,aurostd::xvector<float>& afunc) {
        int i;
        afunc[1]=1.0;
        afunc[2]=x;
        for (i=3;i<=afunc.urows;i++) afunc[i]=sin(i*x);
}

int pinku_main(void) {
        int i,j;
        float chisq;

        xvector<int> ia(1,NTERM);
	xvector<float> a(1,NTERM);
        xvector<float> x(1,NPT);
        xvector<float> y(1,NPT);
        xvector<float> sig(1,NPT);
        xmatrix<float> covar(1,NTERM,1,NTERM);

        for (i=1;i<=NPT;i++) {
                x[i]=0.1*i;
                pinku_funcs(x[i],a);
                y[i]=0.0;
                for (j=1;j<=NTERM;j++) y[i] += j*a[j];
                y[i] += SPREAD*aurostd::ran0();
                sig[i]=SPREAD;
        }
        for (i=1;i<=NTERM;i++) ia[i]=1;
	aurostd::lfit(x,y,sig,a,ia,covar,chisq,pinku_funcs);
        printf("\n%11s %21s\n","parameter","uncertainty");
        for (i=1;i<=NTERM;i++)
                printf("  a[%1d] = %8.6f %12.6f\n",
                        i,a[i],sqrt(covar[i][i]));
        printf("chi-squared = %12f\n",chisq);
        printf("full covariance matrix\n");
        for (i=1;i<=NTERM;i++) {
                for (j=1;j<=NTERM;j++) printf("%12f",covar[i][j]);
                printf("\n");
        }
        printf("\npress RETURN to continue...\n");
        (void) getchar();
        /* Now check results of restricting fit parameters */
        for (i=2;i<=NTERM;i+=2) ia[i]=0;
        lfit(x,y,sig,a,ia,covar,chisq,pinku_funcs);
        printf("\n%11s %21s\n","parameter","uncertainty");
        for (i=1;i<=NTERM;i++)
                printf("  a[%1d] = %8.6f %12.6f\n",
                        i,a[i],sqrt(covar[i][i]));
        printf("chi-squared = %12f\n",chisq);
        printf("full covariance matrix\n");
        for (i=1;i<=NTERM;i++) {
                for (j=1;j<=NTERM;j++) printf("%12f",covar[i][j]);
                printf("\n");
        }
        printf("\n");
       return 0;
}



// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************

//aurostd::StringSubst(tokens.at(0),"","B_h");
//aurostd::StringSubst(tokens.at(0),"","K_sv");
//aurostd::StringSubst(tokens.at(0),"","V_sv");
//aurostd::StringSubst(tokens.at(0),"","W_pv");
//aurostd::StringSubst(tokens.at(0),"","Y_sv");



