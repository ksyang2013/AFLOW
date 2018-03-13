// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// Dane Morgan

#ifndef _AFLOW_PFLOW_FUNCS_CPP_
#define _AFLOW_PFLOW_FUNCS_CPP_

#include "aflow_pflow.h"

// ***************************************************************************
// XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY XRAY
// ***************************************************************************
// Function DebyeWallerFactor
// ***************************************************************************
// All data must be sent in SI units.
// This formula is from B.E. Warren, X-ray Diffraction, eq. 5.9.
// It is based on the Debye approximation and only holds for T>>TDebye.
// Note that this function returns the sqrt of what we usually call the
// DW factor.  This is because this term multiplies the scattering
// factors f, not f^2.  DW^2 would be the appropriate term to modulate
// the intensity, but DW is appropriate to modulate the scattering factors
// (see Warren,eq.3.24).
double DebyeWallerFactor(const double& theta, const double& lambda,
			 const double& temp, const double& debye_temp,
			 const double& mass) {
  double st=sin(theta);
  double h=PLANKSCONSTANT_h;
  double twoB=h*h*temp*12.0/(mass*KBOLTZ*debye_temp*debye_temp);
  double twoM=twoB*st*st/(lambda*lambda);
  double DWfact=exp(-1.0*twoM/2.0); // Use e^-M, not e^-2M.
  return DWfact;
}

// ***************************************************************************
// balanceEquation()
// ***************************************************************************
xvector<double> balanceChemicalEquation(const vector<xvector<double> >& _lhs,const vector<xvector<double> >& _rhs,
    bool reduce,bool normalize,double tol){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="balanceChemicalEquation():";
  xvector<double> dummy;
  if(_lhs.size()==0 || _rhs.size()==0){return dummy;}
  int dim=_lhs[0].rows;
  for(uint i=1;i<_lhs.size();i++){if(_lhs[i].rows!=dim){return dummy;}}
  for(uint i=0;i<_rhs.size();i++){if(_rhs[i].rows!=dim){return dummy;}}

  if(LDEBUG){
    cerr << soliloquy << " lhs=";
    for(uint i=0;i<_lhs.size();i++){
      cerr << _lhs[i] << (i!=_lhs.size()-1?", ":"");
    }
    cerr << endl;
    cerr << soliloquy << " rhs=";
    for(uint i=0;i<_rhs.size();i++){
      cerr << _rhs[i] << (i!=_rhs.size()-1?", ":"");
    }
    cerr << endl;
  }

  vector<xvector<double> > lhs;
  vector<xvector<double> > rhs;
  if(reduce){
    for(uint i=0;i<_lhs.size();i++){lhs.push_back(aurostd::reduceByGCD(_lhs[i],tol));}
    for(uint i=0;i<_rhs.size();i++){rhs.push_back(aurostd::reduceByGCD(_rhs[i],tol));}
    if(LDEBUG){
      cerr << soliloquy << " lhs REDUCED=";
      for(uint i=0;i<lhs.size();i++){
        cerr << lhs[i] << (i!=lhs.size()-1?", ":"");
      }
      cerr << endl;
      cerr << soliloquy << " rhs REDUCED=";
      for(uint i=0;i<rhs.size();i++){
        cerr << rhs[i] << (i!=rhs.size()-1?", ":"");
      }
      cerr << endl;
    }
  }else{lhs=_lhs;rhs=_rhs;}
  
  // needs to organized in the following way
  // [[Mn=2,Cu=1,Fe=0],   //left_hand_side
  // [Mn=0,Cu=5,Fe=3],    //left_hand_side
  // [Mn=-1,Cu=-1,Fe=0],  //right_hand_side
  // [Mn=-2,Cu=-1,Fe=0]]  //right_hand_side
  // i.e. compounds in rows, elements in cols
  // left_hand_side of reaction is POSITIVE
  // right_hand_side of reaction is NEGATIVE

  xmatrix<double> composition_matrix(lhs.size()+rhs.size(),dim);
  int counter;
  for(uint i=0;i<lhs.size();i++){counter=1;for(int j=lhs[i].lrows;j<=lhs[i].urows;j++){composition_matrix(i+1,counter++)=lhs[i][j];}}             //lhs
  for(uint i=0;i<rhs.size();i++){counter=1;for(int j=rhs[i].lrows;j<=rhs[i].urows;j++){composition_matrix(i+1+lhs.size(),counter++)=-rhs[i][j];}} //rhs

  if(LDEBUG){
    cerr << soliloquy << " composition matrix:" << endl;
    cerr << composition_matrix << endl;
  }

  return balanceChemicalEquation(composition_matrix,normalize,tol);
}

//normalize === set first coefficient to 1
xvector<double> balanceChemicalEquation(const xmatrix<double>& _composition_matrix,bool normalize,double tol){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="balanceChemicalEquation():";
  if(_composition_matrix.rows<=_composition_matrix.cols){
    cerr << soliloquy << " composition matrix (m<=n) will NOT yield a viable null space for this analysis" << endl;
    exit(1);
  }
  xmatrix<double> composition_matrix=_composition_matrix;
  xmatrix<double> Q=aurostd::generalHouseHolderQRDecomposition(composition_matrix);
  if(LDEBUG){
    cerr << soliloquy << " Q:" << endl;
    cerr << Q << endl;
  }
  xvector<double> coef(composition_matrix.rows);
  for(int i=1;i<Q.rows+1;i++){coef[i]=Q(i,Q.ucols);}
  if(LDEBUG){cerr << soliloquy << " PRE-coefficients =" << coef << endl;}
  if(normalize){
    if(abs(coef[1])<tol){return coef;} //avoid divide by 0
    coef/=coef[1];
  }
  if(LDEBUG){cerr << soliloquy << " POST-coefficients=" << coef << endl;}
  if(LDEBUG){cerr << soliloquy << " checking: all scalar products should yield 0" << endl;}
  double sum;
  for(int i=1;i<_composition_matrix.cols+1;i++){
    if(LDEBUG){cerr << "component-" << i << ": sum=";}
    sum=0.0;
    for(int j=1;j<_composition_matrix.rows+1;j++){sum+=coef(j)*_composition_matrix(j,i);}
    if(LDEBUG){cerr << sum << endl;}
    if(abs(sum)>_ZERO_TOL_){cerr << soliloquy << " ERRROR! chemical equation was not balanced (run with --debug to see)" << endl;exit(1);}
  }
  return coef;
}

// ***************************************************************************
// ParseChemFormula
// ***************************************************************************
void ParseChemFormula(string& ChemFormula,vector<string>& ChemNameRef,vector<float>& ChemConcRef) {
  /*
    ChemFormula=MgB2.3 -> ChemNameRef=["Mg","B"] and ChemConcRef=[1,2.3]

    input:
    ChemFormula

    output:
    ChemNameRef, ChemConcRef

    wahyu@alumni.duke.edu
  */
  uint nchar;
  //  int i,j,k,L,itmp;
  float AtomConc;
  string AtomSymbol;

  ChemNameRef.clear(); ChemConcRef.clear();
  while(ChemFormula.size()>0) {
    nchar=3;
    if(ChemFormula.size()>nchar-1) {
      //check 3-character atom symbol
      ParseChemFormulaIndividual(3,ChemFormula,AtomSymbol,AtomConc);
      if(AtomSymbol!="") {
        ChemNameRef.push_back(AtomSymbol);
        ChemConcRef.push_back(AtomConc);
      }//if valid 3-character atomic symbol found
      else{
        //check 2-character atom symbol
        ParseChemFormulaIndividual(2,ChemFormula,AtomSymbol,AtomConc);
        if(AtomSymbol!="") {
          ChemNameRef.push_back(AtomSymbol);
          ChemConcRef.push_back(AtomConc);
        }//if valid 2-character atomic symbol found
        else{
          //check 1-character atom symbol
          ParseChemFormulaIndividual(1,ChemFormula,AtomSymbol,AtomConc);
          if(AtomSymbol!="") {
            ChemNameRef.push_back(AtomSymbol);
            ChemConcRef.push_back(AtomConc);
          }//if valid 1-character atomic symbol found
        }//else nchar 1
      }//else nchar 2
    }//if ChemFormula.size()>2
    else{
      nchar=2;
      if(ChemFormula.size()>nchar-1) {
        //check 2-character atom symbol
        ParseChemFormulaIndividual(2,ChemFormula,AtomSymbol,AtomConc);
        if(AtomSymbol!="") {
          ChemNameRef.push_back(AtomSymbol);
          ChemConcRef.push_back(AtomConc);
        }//if valid 2-character atomic symbol found
        else{
          //check 1-character atom symbol
          ParseChemFormulaIndividual(1,ChemFormula,AtomSymbol,AtomConc);
          if(AtomSymbol!="") {
            ChemNameRef.push_back(AtomSymbol);
            ChemConcRef.push_back(AtomConc);
          }//if valid 1-character atomic symbol found
        }//else nchar 1
      }//if ChemFormula.size()>1
      else{
        nchar=1;
        ParseChemFormulaIndividual(1,ChemFormula,AtomSymbol,AtomConc);
        if(AtomSymbol!="") {
          ChemNameRef.push_back(AtomSymbol);
          ChemConcRef.push_back(AtomConc);
        }//if valid 1-character atomic symbol found    
      }//if ChemFormula.size()==1
    }
  }//while not empty

}
// ***************************************************************************
// ParseChemFormulaIndividual
// ***************************************************************************
void ParseChemFormulaIndividual(uint nchar, string& ChemFormula, string& AtomSymbol, float& AtomConc) {
  /*
    get 1 valid pair of atom symbol and its concentration from ChemFormula.
    start from the beginning of ChemFormula. If found, return the AtomSymbol and AtomConc
    and truncate it from the ChemFormula.

    input:
    nchar     : the number of characters in the atom symbol

    wahyu@alumni.duke.edu
   */

  uint iret;
  int i;
  string sN,sC;
  
  AtomSymbol="";
  AtomConc=0.0;
  sN=ChemFormula.substr(0,nchar);
  iret=GetAtomNumber(sN);
  if(iret>0) {//valid atomic symbol
    AtomSymbol=sN;
    if(ChemFormula.size()==nchar) {//last atom with no concentration
      ChemFormula="";
      AtomConc=1.0;
    }
    else{
      ChemFormula=ChemFormula.substr(nchar,ChemFormula.size());
      if((ChemFormula[0]>'0' && ChemFormula[0]<'9') || ChemFormula[0]=='.') {
        //parse the concentration
        for(i=0;i<(int)ChemFormula.size();i++) {
          if(!((ChemFormula[i]>'0' && ChemFormula[i]<'9') || ChemFormula[i]=='.')) break;
        }
        sC=ChemFormula.substr(0,i);
        if(sC[0]=='.') sC="0"+sC;//append 0 in case the first character is '.'
        if(sC[sC.size()-1]=='.') sC=sC+"0";
        AtomConc=aurostd::string2utype<float>(sC);
        if(i==(int)ChemFormula.size()) ChemFormula="";
        else ChemFormula=ChemFormula.substr(i,ChemFormula.size());
      }
      else  AtomConc=1.0;//set concentration to 1.0 if atom symbol is valid but concentration is not given
    }
  }//if valid atomic symbol found  
}

// ***************************************************************************
// GetXRAY
// ***************************************************************************
// This function gets XRAY following the convasp framework.
// Dane Morgan
namespace pflow {
  void GetXray(const xstructure& str, vector<double>& dist,vector<double>& sf,
	       const double& lambda, vector<double>& scatt_fact,
	       vector<double>& mass, vector<double>& twoB_vec) {
    // Get data from str.
    // Set scale to 1 so you don't need to rescale coordinates.
    xstructure sstr=str;
    sstr=ReScale(sstr,1.0);
  
    xmatrix<double> rlat(3,3);rlat=sstr.klattice;
    int num_atoms=sstr.atoms.size();
    //  matrix<double> fpos=GetFpos(sstr);

    // Set parameters for Debye Waller factors.
    double temp=300;
    double debye_temp=300;
    double dw,theta;

    // Get scattering factors and masses and Debye-Waller 2B values.
    double h=PLANKSCONSTANT_h;
    for(int i=0;i<num_atoms;i++) {
      _atom atom=sstr.atoms.at(i); atom.CleanName();
      scatt_fact[i]=GetXrayScattFactor(atom.cleanname,lambda);
      mass[i]=GetAtomMass(atom.cleanname);
      twoB_vec[i]=h*h*temp*12.0/(mass[i]*KBOLTZ*debye_temp*debye_temp);
    }
    //  Get max h,k,l value that gives any scattering (theta<=90).
    int kmx0=(int) (4.0*PI/(modulus(rlat(1))*lambda)+1);
    int kmx1=(int) (4.0*PI/(modulus(rlat(2))*lambda)+1);
    int kmx2=(int) (4.0*PI/(modulus(rlat(3))*lambda)+1);
    int kmx=max(kmx0,kmx1);
    kmx=max(kmx2,kmx);
    int len=2*kmx+1;
    int tlen=len*len*len;  
    double sfr,sfi;
    dist = vector<double> (tlen,0.0);
    sf = vector<double> (tlen,0.0);
  
    for(int i0=-kmx;i0<=kmx;i0++) {
      for(int i1=-kmx;i1<=kmx;i1++) {
	for(int i2=-kmx;i2<=kmx;i2++) {
	  int ii0=i0+kmx;
	  int ii1=i1+kmx;
	  int ii2=i2+kmx;
	  int id=ii2+ii1*len+ii0*len*len;
	  sfr=0.0;
	  sfi=0.0;
	  xvector<double> rv(3);
	  for(int ic=1;ic<=3;ic++)
	    rv(ic)=i0*rlat(1,ic)+i1*rlat(2,ic)+i2*rlat(3,ic);
	  double rvnorm=modulus(rv);
	  if(rvnorm>0.0) dist[id]=TWOPI/rvnorm;
	  for(int iat=0;iat<num_atoms;iat++) {
	    // Get Debye Waller factor.
	    dw=1;
	    if(dist[id]>0) {
	      double term=lambda/(2.0*dist[id]);
	      if(term<=1) {
		theta=std::asin(term);
		theta=theta*360.0/TWOPI; // rad->degrees
		dw=DebyeWallerFactor(theta,lambda*1.0E-10,temp,debye_temp,mass[iat]);
	      }
	    }
	    //  double gdotr=i0*fpos[iat][0]+i1*fpos[iat][1]+i2*fpos[iat][2];
	    double gdotr=i0*sstr.atoms.at(iat).fpos(1)+i1*sstr.atoms.at(iat).fpos(2)+i2*sstr.atoms.at(iat).fpos(3);
	    gdotr=TWOPI*gdotr;
	    sfr=sfr+dw*scatt_fact[iat]*std::cos(gdotr);
	    sfi=sfi-dw*scatt_fact[iat]*std::sin(gdotr);
	  }
	  sf[id]=sfr*sfr+sfi*sfi;
	} // i2
      } // i1
    } // i0
  }
}

// ***************************************************************************
// RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF RDF
// ***************************************************************************
// Function GetRDF
// ***************************************************************************
// This function gets the radial distribution functions (RDF).
// The RDF for an atom is just a binned histogram of the
// number of atoms (possibly of a given type) a distance r away.  
// The RDF for a type is just the average RDF's for all the atoms
// of that type.
// The rdf is stored as follows: I=(0,natoms-1),J,K=(0,ntypes)
// row: J+(ntypes+1)*I = atom I / Type J RDF.
//    (if row: J=ntypes then the rdf is the atom I / All types RDF.)
// I don't think the following ever got coded.
// lines (ntypes+1)*natoms+K+(ntypes+1)*J = type J / type K average RDF.
//    (if row: K=ntypes then the rdf is the type I / All types RDF.)
// Dane Morgan, Modified by Stefano Curtarolo
namespace pflow {
  void GetRDF(xstructure str, const double& rmax,
	      const int& nbins, matrix<double>& rdf_all) {
    int natoms=str.atoms.size();
    std::deque<int> num_each_type=str.num_each_type;
    int ntyp=num_each_type.size();
    rdf_all=matrix<double> ((ntyp+1)*natoms,nbins);
    deque<deque<_atom> > nmat;
    // [OBSOLETE]    pflow::GetStrNeighData(str,rmax,nmat);
   str.GetStrNeighData(rmax,nmat);   // once GetRD goes in xstructure I can remove the copy
    // for(int i=0;i<nmat[0].size();i++) cout << AtomDist(nmat[0][0],nmat[0][i]) << " "; cout << endl; exit(0);
    for(int I1=0;I1<(int)nmat.size();I1++) { // Each atom for which we find RDF.
      int I2=1;
      double dist=0;
      int s2=(int)nmat.at(I1).size();
      _atom a1=nmat.at(I1).at(0);
      // Check each atom neighbor for I1 in proper range.
      dist=AtomDist(a1,nmat.at(I1).at(I2));
      while (dist<rmax && I2<s2) {
	int ib=int((dist/rmax)*nbins); // first bin is [0,rmax/nbins).
	int J2=nmat.at(I1).at(I2).type;     // CONVASP_MODE
	//      cout << I1 << " " << I2 << " " <<  J2 << " " << s2 << " " << dist << " " << ib << endl;
	rdf_all[(ntyp+1)*I1+J2][ib]++; // Does all binning for atom/type RDFs.
	I2++;
	if(I2<s2) dist=AtomDist(a1,nmat[I1][I2]);
      }  
      // Get sum over all types
      for(int ib=0;ib<nbins;ib++) {
	for(int it=0;it<ntyp;it++) {
	  rdf_all[(ntyp+1)*I1+ntyp][ib]+=rdf_all[(ntyp+1)*I1+it][ib];
	}
      }
    } // I1 loop
  }
}
// ***************************************************************************
// Function GetRDFShells
// ***************************************************************************
//   This function gets the neighbor shells from the
//   radial distribution functions.  Shells are found
//   by the following method.  The code calculates
//   the change in rdf from drdf=rdf[i+1]-rdf[i].  We
//   step through drdf starting adding up rdf into shell1
//   while drdf>=0 and then while drdf<=0.  When drdf gets
//   >=0 again we start a new shell until drdf passes
//   through <=0 and is >=0 again, at which point we again
//   start a new shell.  In other words, we increment the
//   shell we are adding to when drdf goes from >=0 to <0 to
//  >=0 again.
namespace pflow {
  void GetRDFShells(const xstructure& str,const double& rmax,const int& nbins,
		    const int& smooth_width,const matrix<double>& rdf,
		    matrix<double>& rdfsh,matrix<double>& rdfsh_loc) {
    // double TOL=1e-5; // DANE not used
    // int natom=(int)str.atoms.size(); // DANE not used
    if(smooth_width) {;} // phony just to keep smooth_width busy
    if(str.atoms.size()) {;} // phony just to keep str busy

    int _rdi=(int)rdf.size();
    rdfsh = matrix<double>  (_rdi,0); pflow::VVset(rdfsh,0);
    rdfsh_loc =  matrix<double> (_rdi,0);  pflow::VVset(rdfsh_loc,0);
    double dr=(rmax/(double)nbins);
    for(int i=0;i<_rdi;i++) {

      /*
      // get smoothed rdf_sm
      vector<double> rdf_sm(nbins,0.0);
      int cnt=0;
      for(int ib=0;ib<nbins;ib++) {
      for(int ism=-(smooth_width-1)/2;ism<=(smooth_width-1)/2;ism++) {
      int id=ib+ism;
      if(id>0 && id<nbins) {
      cnt++;
      rdf_sm[ib]+=rdf[i][ib+ism];
      }
      }
      if(cnt>0) rdf_sm[ib]=rdf_sm[ib]/cnt;
      cnt=0;
      }
      rdf_sm=SmoothFunc(rdf[i],smooth_width);
      */

      // get drdf,ddrdf of rdf
      vector<double> drdf(nbins,0.0);
      vector<double> ddrdf(nbins,0.0);
      for(int ib=1;ib<nbins-1;ib++) {
	drdf[ib]=(rdf[i][ib+1]-rdf[i][ib-1]);
      }
      for(int ib=1;ib<nbins-1;ib++) {
	//      ddrdf[ib]=(rdf[i][ib+1]+rdf[i][ib-1]-2*rdf[i][ib]);
	ddrdf[ib]=(drdf[ib+1]-drdf[ib-1]);
      }

      // get all zeros of drdf with ddrf!=0
      // -999 means nothing,-1 means ddrdf<0, 0 means ddrdf=0, +1 means ddrdf>0.
      vector<int> drdf_zeros(nbins-1,-999);
      for(int ib=1;ib<nbins-1;ib++) {
	// If drdf==0
	if(drdf[ib]==0 && ddrdf[ib]<0) drdf_zeros[ib]=-1;
	if(drdf[ib]==0 && ddrdf[ib]==0) drdf_zeros[ib]=0;
	if(drdf[ib]==0 && ddrdf[ib]==1) drdf_zeros[ib]=1;
	// If drdf== is passing through 0 from +->- or -->+
	if(drdf[ib]!=0) {
	  if(drdf[ib-1]>0 && drdf[ib+1]<0) drdf_zeros[ib]=-1;
	  if(drdf[ib-1]<0 && drdf[ib+1]>0) drdf_zeros[ib]=+1;
	  // If drdf== is becoming or changinf from 0 through
	  // 0->+/- or +/-->0.
	  if(drdf[ib-1]==0 && ddrdf[ib]<0) drdf_zeros[ib]=-1;
	  if(drdf[ib-1]==0 && ddrdf[ib]==0) drdf_zeros[ib]=0;
	  if(drdf[ib-1]==0 && ddrdf[ib]>0) drdf_zeros[ib]=1;
	  if(drdf[ib+1]==0 && ddrdf[ib]<0) drdf_zeros[ib]=-1;
	  if(drdf[ib+1]==0 && ddrdf[ib]==0) drdf_zeros[ib]=0;
	  if(drdf[ib+1]==0 && ddrdf[ib]>0) drdf_zeros[ib]=1;
	}

	// tpx
	//      cout << "drdf_zeros " << ib << " " << rdf[i][ib] << " " << drdf[ib] << " " << ddrdf[ib] << " " << drdf_zeros[ib] << endl;
      }

      // get the actual shell atoms counts and avg radius.
      double shell=0;
      double rsh_avg=0;
      double state_dn=0;
      for(int ib=0;ib<nbins-2;ib++) {
	//tpx
	//           cout << i << " " << ib << " "
	//	         << rdf[i][ib] << " " << drdf[ib] << " " << shell << " " << endl;
	shell=shell+rdf[i][ib];
	double rad=(dr*(double)ib)+dr/2.0;
	rsh_avg=rsh_avg+rdf[i][ib]*rad;
	if(drdf_zeros[ib]==-1) { // Set state_dn
	  state_dn=1;
	}
	if(state_dn && drdf_zeros[ib]==1) { // New shell
	  if(shell>0) rsh_avg=rsh_avg/shell;
	  rdfsh[i].push_back(shell);
	  rdfsh_loc[i].push_back(rsh_avg);
	  shell=0;
	  rsh_avg=0;
	  state_dn=0;
	}
      } // for ib
    } // for i
  }
}

// ***************************************************************************
// Function RdfSh_RMS
// ***************************************************************************
// This function compares the rdf shells of two
// atoms for each type and returns the rms.  
namespace pflow {
  double RdfSh_RMS(const int iaA, const int iaB, const int nsh_max,const int nt,
		   const matrix<double>& rdfsh_all_A
		   ,const matrix<double>& rdfsh_all_B) {
    double rms=0;
    int cnt=0;
    for(int it=0;it<nt;it++) {
      int idA=(nt+1)*iaA+it;
      int idB=(nt+1)*iaB+it;
      int tempsh=min((int) rdfsh_all_A[idA].size(),(int) rdfsh_all_B[idB].size());
      int nsh = min(tempsh,nsh_max);
      for(int ish=0;ish<nsh;ish++) {
	cnt++;
	rms+=(rdfsh_all_A[idA][ish]-rdfsh_all_B[idB][ish])
	  *(rdfsh_all_A[idA][ish]-rdfsh_all_B[idB][ish]);
	// tpx
	//      cout << "iaA,iaB,it,ish,cnt,rms " <<iaA<<" "<<iaB<<" "<<it<<" "<<ish<<" "<< cnt << " " << rms << endl;
      }
    }
    if(cnt>0) rms=sqrt(rms/cnt);
    return rms;
  }
}

// ***************************************************************************
// Function CmpRDFShell
// ***************************************************************************
// This function compares the rdf shells of two structures.
// Two atoms similarity are based on the rms error between
// their first nsh nn shells (for every type).  Atoms must
// be the same type to even be compared.
// The results are given in:
// best_match: For each atom in str_A gives the best matching atom
//             in str_B of the smae type.
// rms_mat:  Gives the rms error between every pair of atoms.
// For the best matches I start with the first atom of str_A and
// compare to all of B.  Each successive match for an atom of str_A
// is performed excluding the previously matched atoms of str_B.
// str_A and str_B must have the same number of each type of atom.
namespace pflow {
  void CmpRDFShells(const xstructure& str_A, const xstructure& str_B,
		    const matrix<double>& rdfsh_all_A,
		    const matrix<double>& rdfsh_all_B,
		    const int nsh, vector<int>& best_match,
		    matrix<double>& rms_mat) {
    double TOL=1e-15;
    std::deque<int> netype_A=str_A.num_each_type;
    std::deque<int> netype_B=str_B.num_each_type;
    // Exit if A and B have different numbers of any types of atoms.
    if(!pflow::VVequal(netype_A,netype_B) || netype_A.size()!=netype_B.size()) {
      cerr << "ERROR: in CmpRDFShell" << endl;
      cerr << "ERROR: structures A and B do not have the same "<< endl;
      cerr << "ERROR: number of each type of atom. "<< endl;
      cerr << "ERROR: Exiting" << endl;
      exit(1);
    }
    int nat=str_A.atoms.size();
    int nt=netype_A.size();

    // Get rms_mat
    rms_mat = matrix<double> (nat,nat); pflow::VVset(rms_mat,-1.0);

    for(int iaA=0;iaA<nat;iaA++) {
      for(int iaB=0;iaB<nat;iaB++) {
	rms_mat[iaA][iaB]=RdfSh_RMS(iaA,iaB,nsh,nt,rdfsh_all_A,rdfsh_all_B);
      }
    }
    // Get best_match matrix.
    // For each atom in A check the rms error for all atoms of that
    // type in B.  Then take the best one and then remove it from
    // the list for further comparisons.
    best_match=vector<int> (nat,-1);
    vector<int> cand_Batoms(nat,0);
    for(int i=0;i<nat;i++) cand_Batoms[i]=i;
    for(int iaA=0;iaA<(int)rms_mat.size();iaA++) {
      double rms=-10;
      double trms=-10;
      int best_at_id=-1;
      int typeA=str_A.atoms.at(iaA).type;
      // Check each B atom and if it is the right type and lowest rms
      // then track its id.
      for(int ipB=0;ipB<(int)cand_Batoms.size();ipB++) {
	int iaB=cand_Batoms[ipB];
	int typeB=str_B.atoms.at(iaB).type;
	if(typeA==typeB) {
	  trms=rms_mat[iaA][iaB];
	  if(rms<-1) rms=trms+1.0;
	  if(trms<rms-TOL) {
	    best_at_id=iaB;
	    rms=trms;
	  }
	}
      }
      best_match[iaA]=best_at_id;
      // Remove best_at_id from future comparisons.
      cand_Batoms.erase(find(cand_Batoms.begin(),cand_Batoms.end(),best_at_id));
    }
  }
}

// ***************************************************************************
// Function GetSmoothRDF
// ***************************************************************************
// This function gets a smoothed RDF.
namespace pflow {
  matrix<double> GetSmoothRDF(const matrix<double>& rdf,
			      const double& sigma) {
    matrix<double> rdf_sm=rdf;
    for(int i=0;i<(int)rdf.size();i++) {
      rdf_sm[i]=pflow::SmoothFunc(rdf[i],sigma);
    }
    return rdf_sm;
  }
}
// ***************************************************************************
// COMPARE STRUCTURES COMPARE STRUCTURES COMPARE STRUCTURES COMPARE STRUCTURES
// ***************************************************************************
// ***************************************************************************
// CmpStrDist
// ***************************************************************************
// Compares the distances within rcut between str1 and str2.
// Assigns matrix of distances.  Each row is for a different
// pair type. Pair types are identified by pair i,j in row
// j+ntypes*i.  1<= i,j <=ntypes. I always take j>=i.
// Original by Dane Morgan, modified by STefano Curtarolo for type shift !
namespace pflow {
  void CmpStrDist(xstructure str1, xstructure str2,const double& cutoff,
		  matrix<double>& dist1, matrix<double>& dist2,
		  matrix<double>& dist_diff,matrix<double>& dist_diff_n) {

    deque<deque<_atom> > neigh_mat1;
    deque<deque<_atom> > neigh_mat2;
    // [OBSOLETE] pflow::GetStrNeighData(str1,cutoff,neigh_mat1);
    // [OBSOLETE] pflow::GetStrNeighData(str2,cutoff,neigh_mat2);
    str1.GetStrNeighData(cutoff,neigh_mat1);
    str2.GetStrNeighData(cutoff,neigh_mat2);

    int ntypes1=(int)str1.num_each_type.size();
    int ntypes2=(int)str2.num_each_type.size();
    int nprs1=(ntypes1*(ntypes1+1))/2;
    int nprs2=(ntypes2*(ntypes2+1))/2;
    vector<double> tmpvec;
    dist1=matrix<double> (nprs1,tmpvec);
    dist2=matrix<double> (nprs2,tmpvec);

    // Collect distances.
    // dist1
    for(int ia=0;ia<(int)neigh_mat1.size();ia++) {
      _atom a = neigh_mat1[ia][0];
      for(int in=1;in<(int)neigh_mat1[ia].size();in++) {
	_atom an = neigh_mat1[ia][in];
	int ta=a.type;              // CONVASP_MODE
	int tna=an.type;            // CONVASP_MODE
	int i=std::min(ta,tna);
	int j=std::max(ta,tna);
	int id=j-i+ntypes1*i-max(0,i*(i-1)/2);
	dist1[id].push_back(AtomDist(a,an));
      } // in
    } // ia
    // dist2
    for(int ia=0;ia<(int)neigh_mat2.size();ia++) {
      _atom a = neigh_mat2[ia][0];
      for(int in=1;in<(int)neigh_mat2[ia].size();in++) {
	_atom an = neigh_mat2[ia][in];
	int ta=a.type;            // CONVASP_MODE
	int tna=an.type;          // CONVASP_MODE
	int i=std::min(ta,tna);
	int j=std::max(ta,tna);
	int id=j-i+ntypes2*i-max(0,i*(i-1)/2);
	dist2[id].push_back(AtomDist(a,an));
      } // in
    } // ia

    // Sort dist vectors.
    for(int ip=0;ip<nprs1;ip++) {
      sort(dist1[ip].begin(),dist1[ip].end());  
    }
    for(int ip=0;ip<nprs2;ip++) {
      sort(dist2[ip].begin(),dist2[ip].end());  
    }

    // tpx
    //  pflow::Vout(dist1[0],cout);
    //  pflow::Vout(dist2[0],cout);
    // Assign distance difference matrix
    int ntypes_min=std::min(ntypes1,ntypes2);
    int nprs_min=ntypes_min*(ntypes_min+1)/2;
    dist_diff = matrix<double> (nprs_min,tmpvec);
    dist_diff_n = matrix<double> (nprs_min,tmpvec);
    for(int i=0;i<ntypes_min;i++) {
      for(int j=i;j<ntypes_min;j++) {
	int id1=j-i+ntypes1*i-max(0,i*(i-1)/2);
	int id2=j-i+ntypes2*i-max(0,i*(i-1)/2);
	int idmin=j-i+ntypes_min*i-max(0,i*(i-1)/2);
	int num_dist=std::min((int)dist1[id1].size(),(int)dist2[id2].size());
	for(int ip=0;ip<num_dist;ip++) {
	  double d1=dist1[id1][ip];
	  double d2=dist2[id2][ip];
	  dist_diff[idmin].push_back(d2-d1);
	  dist_diff_n[idmin].push_back(2*(d2-d1)/(d1+d2));
	}// ip
      }// j
    }// i
  } // end routine
}

// ****************************************************************************************************
// PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA P
// ****************************************************************************************************
namespace pflow {
  // Constructors
  pdosdata::pdosdata() {
  }
  void pdosdata::PrintParams(ostream& outf, const std::vector<string>& Ltotnames) {
    outf << "EMIN = " << emin << endl;
    outf << "EMAX = " << emax << endl;
    outf << "EFERMI = " << efermi << endl;
    outf << "NBINS = " << nbins << endl;
    outf << "SPIN = " << spin << endl;
    outf << "NLM = " << nlm << endl;
    outf << "NATOMS = " << natoms << endl;
    outf << "SMOOTH_SIGMA = " << smooth_sigma << endl;
    outf << "PRINT_PARAMS = " << print_params << endl;
    for(int i=0;i<(int)pdos_at.size();i++) {
      outf << "# CASE = " << i+1 << endl;
      outf << "  ATOMS = ";
      Vout(pdos_at[i],outf);
      if(pdos_k.size()>0) {
	outf << "  KPOINTS = ";
	Vout(pdos_k[i],outf);
      }
      if(pdos_b.size()>0) {
	outf << "  BANDS = ";
	Vout(pdos_b[i],outf);
      }
      outf << "  LMVALUES = ";
      Vout(pdos_lm[i],outf);
      outf << "  # ( LMVALUES = ";
      for(int j=0;j<(int)pdos_lm[i].size();j++) {
	outf << " " << Ltotnames[pdos_lm[i][j]-1];
      }
      outf << " )" << endl;
    }
  }

  void pdosdata::PrintPDOS(ostream& outf, const int& sp) {
    outf.precision(5);
    outf.setf(std::ios::fixed,std::ios::floatfield);
    outf.setf(std::ios::left,std::ios::adjustfield);
    
    if(sp==0) outf << "# Energy         Up            Cumulative_Up" << endl;
    if(sp==1) outf << "# Energy         Up             Dn             Up-Dn        Cumulative_Up   Cumulative_Dn   Cumulative_Up-Dn" << endl;
    for(int ib=0;ib<(int)pdos.size();ib++) {
      for(int i=0;i<(int)pdos[ib].size();i++) {
	outf << setw(10) << pdos[ib][i] << "     ";
      }
      outf << endl;
    }
  }
}

// ****************************************************************************************************
// RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY
// ****************************************************************************************************
namespace pflow {
  void rtparams::free() {
  }
  
  void rtparams::copy(const rtparams& b) {
    resx=b.resx;
    resy=b.resy;
    viewdir=b.viewdir;
    viewdir_s=b.viewdir_s;
    updir=b.updir;
    updir_s=b.updir_s;
    zoom=b.zoom;
    aspectratio=b.aspectratio;
    antialiasing=b.antialiasing;
    raydepth=b.raydepth;
    center=b.center;
    center_guide=b.center_guide;
    center_s=b.center_s;
    background=b.background;
    lightcenter=b.lightcenter;
    lightrad=b.lightrad;
    lightcolor=b.lightcolor;
    sphtex_tex=b.sphtex_tex;
    sphtex_tex_def=sphtex_tex_def;
    sphtex_color=b.sphtex_color;
    sphtex_color_def=b.sphtex_color_def;
    sphtex_phong=b.sphtex_phong;
    sphtex_phong_def=b.sphtex_phong_def;
    sphtex_names=b.sphtex_names;
    sph_rad=b.sph_rad;
    sph_rad_def=b.sph_rad_def;
    shading=b.shading;
    outfile=b.outfile;
    sc=b.sc;
    sc_s=b.sc_s;
    calc_type=b.calc_type;
    input_files=b.input_files;
    first_set=b.first_set;
    insert_file=b.insert_file;
    rotation=b.rotation;
    // Plane variables
    plane=b.plane;
    plane_s=b.plane_s;
    plane_center=b.plane_center;
    plane_center_def=b.plane_center_def;
    plane_center_s=b.plane_center_s;
    plane_normal=b.plane_normal;
    plane_normal_def=b.plane_normal_def;
    plane_normal_s=b.plane_normal_s;
    plane_color=b.plane_color;
    plane_color_def=b.plane_color_def;
    plane_color_s=b.plane_color_s;
    planetex_tex=b.planetex_tex;
    planetex_tex_def=b.planetex_tex_def;
    planetex_tex_s=b.planetex_tex_s;
  }

  // Constructors
  rtparams::rtparams() {
    calc_type = 0;
    resx=600;
    resy=600;
    viewdir = vector<double> (3,0.0);
    viewdir[1]=1;
    viewdir_s = 0;
    updir = vector<double> (3,0.0);
    updir[2]=1;
    updir_s=0;
    zoom=1;
    aspectratio=1;
    antialiasing=0;
    raydepth=12;
    center = vector<double> (3,0.0);
    center_guide = vector<double> (6,0.0);
    center_s = 0;
    background = vector<double> (3,1.0);
    lightcenter = matrix<double> (1,center);
    lightrad = vector<double> (1,0.001);
    vector<double> color (3,0.0);
    lightcolor = matrix<double> (1,color);
    sphtex_tex_def = vector<double> (4);
    sphtex_tex_def[0] = 0.3; // ambient
    sphtex_tex_def[1] = 0.6; // diffuse
    sphtex_tex_def[2] = 0.2; // specular
    sphtex_tex_def[3] = 1.0; // opacity
    sphtex_color_def = vector<double> (3,0.5);
    sphtex_phong_def = vector<double> (2,0.0);
    sphtex_phong_def[1] = 10000;
    sph_rad_def = 1;
    shading = "fullshade";
    outfile = "POSCAR_RT";
    sc=matrix<double> (3,3); pflow::VVset(sc,0.0);
    sc[0][0]=1;
    sc[1][1]=1;
    sc[2][2]=1;
    sc_s=0;
    first_set=1;
    insert_file="NO_INSERT_FILE";
    rotation = vector<double> (6,0.0);
    struct_origin = vector<double> (3,0.0);
    //  int struct_origin_s=0;
    // Plane variables
    plane=0;
    plane_s=0;
    plane_center = vector<double> (3);
    plane_normal = vector<double> (3);
    plane_color = vector<double> (3);
    planetex_tex = vector<double> (4);
    plane_center_def = vector<double> (3);
    plane_normal_def = vector<double> (3);
    plane_color_def = vector<double> (3);
    planetex_tex_def = vector<double> (4);
    plane_center_s=0;
    plane_normal_s=0;
    plane_color_s=0;
    planetex_tex_s=0;
  }
  
  rtparams::rtparams(const rtparams& b) {copy(b);}

  const rtparams& rtparams::operator=(const rtparams& b) {
    if(this != &b) {
      free();
      copy(b);
    }
    return *this;
  }

  void rtparams::Print() const {
    cout << "CENTER = " << center_guide[0] << " " << center_guide[1] << " " << center_guide[2] << " " << center_guide[3] << " " << center_guide[4] << " " << center_guide[5] << endl;
    cout << "VIEWDIR = " << viewdir[0] << " " << viewdir[1] << " " << viewdir[2] << endl;
    cout << "UPDIR = " << updir[0] << " " << updir[1] << " " << updir[2] << endl;
    cout << "STRUCTURE_ORIGIN = " << struct_origin[0] << " " << struct_origin[1] << " " << struct_origin[2] << endl;
    cout << "ROTATION = " << rotation[0] << " " << rotation[1]
	 << " " << rotation[2] << " " << rotation[3]
	 << " " << rotation[4] << " " << rotation[5] << endl;
  }
}

// ***************************************************************************
// GetDatFromOUTCAR
// ***************************************************************************
// This gets the lattice vectors and num_each_type from an
// OUTCAR file.
// Dane Morgan style
namespace pflow {
  void GetDatFromOutcar(vector<matrix<double> >& lat_vec,
			deque<int>& num_each_type, ifstream& outcar_inf) {
    int first_lat=1;
    vector<string> a(3,"Z");
    lat_vec.clear();
    string sdum;
    while (outcar_inf >> a[2]) {
      string key= (a[0]+" "+a[1]+" "+a[2]);
      if(key=="reciprocal lattice vectors") {
	matrix<double> lat(3,3);
	for(int ic=0;ic<3;ic++) {
	  outcar_inf >> lat[ic][0] >> lat[ic][1] >> lat[ic][2];
	  outcar_inf >>sdum>>sdum>>sdum;
	}
	if(!first_lat) { // skip first lat.
	  lat_vec.push_back(lat);
	}
	first_lat=0;
      }
      if(key=="ions per type") {
	string s;
	int id=0;
	getline(outcar_inf,s);
	aurostd::GetNextVal(s,id).c_str(); // remove an "=" sign.
	while (id<(int)s.size()) {
	  int n=atoi(aurostd::GetNextVal(s,id).c_str());
	  num_each_type.push_back(n);
	}
      }
      a[0]=a[1];
      a[1]=a[2];
    }
  }

}

// ***************************************************************************
// GetDatFromXDATCAR
// ***************************************************************************
// This gets the fpos from XDATCAR.
namespace pflow {
  void GetDatFromXdatcar(vector<matrix<double> >& fpos_vec,
			 ifstream& xdatcar_inf)
  {
    fpos_vec.clear();
    string s;
    matrix<double> dpmat(0);
    vector<double> dpvec(3);
    // Read in header lines
    for(int il=0;il<6;il++) {
      getline(xdatcar_inf,s);
    }
    // Read in all sets of fpos.
    int keep_reading=1;
    while (keep_reading) {
      getline(xdatcar_inf,s);
      int id;
      string key;
      while (s.size()>1 && key!="Konfig=") {
	id=0;
	for(int ic=0;ic<3;ic++) {
	  dpvec[ic]=atof(aurostd::GetNextVal(s,id).c_str());
	}
	dpmat.push_back(dpvec);
	if(!getline(xdatcar_inf,s)) keep_reading=0;
	id=0;
	key=aurostd::GetNextVal(s,id);
      }
      fpos_vec.push_back(dpmat);
      dpmat=matrix<double> (0);
    }
  }
}

// ***************************************************************************
// GetStrVecFromOUTCAR_XDATCAR
// ***************************************************************************
namespace pflow {
  vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(ifstream& outcar_inf, ifstream& xdatcar_inf) {
    vector<matrix<double> > fpos_vec;
    vector<matrix<double> > lat_vec;
    deque<int> num_each_type;
    pflow::GetDatFromOutcar(lat_vec,num_each_type,outcar_inf);
    pflow::GetDatFromXdatcar(fpos_vec,xdatcar_inf);
    if(lat_vec.size()!=fpos_vec.size()) {
      cout << endl;
      cerr << "WARNING: RayTraceFuncs.cc/GetStrVecFromOUTCAR_XDATCAR" << endl;
      cerr << "WARNING: Number of lattice vector and positions of images are not equal." << endl;
      cerr << "WARNING: Number of lattices: " << lat_vec.size()<< endl;
      cerr << "WARNING: Number of positions: " << fpos_vec.size()<< endl;
      cerr << "WARNING: Resizing number of lattices to match number of positions and using last lattice to fill out vector if needed. " << endl;
      cout << endl;
    }
    if(lat_vec.size()>fpos_vec.size()) {lat_vec.resize(fpos_vec.size());}
    if(lat_vec.size()<fpos_vec.size()) {
      matrix<double> tlat=lat_vec[lat_vec.size()-1];
      for(int i=(int)lat_vec.size();i<(int)fpos_vec.size();i++) {
	lat_vec.push_back(tlat);
      }
    }
    vector<xstructure> vstr;
    xstructure str;
    int nstr=(int)lat_vec.size();
    for(int is=0;is<nstr;is++) {
      int nat=fpos_vec[0].size();
      vector<string> names(nat,"H");
      vector<int> names_were_given(nat,FALSE);
      str=pflow::SetLat(str,lat_vec[is]);
      str=pflow::SetNumEachType(str,num_each_type);
      str=pflow::AddAllAtomPos(str,fpos_vec[is],0);
      str=pflow::SetAllAtomNames(str,names);
      str=pflow::SetNamesWereGiven(str,names_were_given);
      vstr.push_back(str);
    }
    return vstr;
  }
}

// ***************************************************************************
// PrintStrVec
// ***************************************************************************
namespace pflow {
  void PrintStrVec(const vector<xstructure>& vstr, ostream& outf) {
    if(vstr.size()==0) return;
    for(int ist=0;ist<(int)vstr.size()-1;ist++) {
      outf << vstr[ist];
      outf << endl;
    }
    outf << vstr[vstr.size()-1];
  }
}

// ***************************************************************************
// ReadInRTParams
// ***************************************************************************
// Reads the ray tracing parameters and stores them in the rtparam object.
// Input data is stored in rtinfile using tokens of the form
// TOKEN=value.  The input can have arbitrary spaces and blank lines.
// The input can also have arbitrary stuff after value as long as it
// is separated by a space.  The value can sometimes be more than one
// field.  Comments are any sequence of text on a single line that
// follows a "#".
namespace pflow {
  void ReadInRTParams(ifstream& rtinfile, pflow::rtparams& rtp) {

    // Read in all the tokens.
    string s,s_ns;
    vector<string> token_vec;
    vector<string> val_vec;
    while (!rtinfile.eof()) {
      getline(rtinfile,s);
      s_ns=aurostd::RemoveSpaces(s); // Get string with no spaces.
      if(s_ns.size()>0) { // Make sure line was not blank
	if(s_ns[0]!='#') { // Exclude comment lines
	  string token,sval;
	  int id=s.find('=');
	  if(id>=(int)s.length()) {
	    cout << "ERROR: The following token is incorrectly formatted: " << s << endl;
	    cout << "ERROR: Exiting" << endl;
	    exit(1);
	  }
	  token=s.substr(0,id);
	  token=aurostd::RemoveSpaces(token);
	  int i_f=std::min(s.length(),s.find('#')); // End of string
	  sval=s.substr(id+1,i_f-id-1);
	  token_vec.push_back(token);
	  val_vec.push_back(sval);
	} //if s_ns[0]!=#
      } // if s_ns.size>0

    } // while !rtinfile.eof

    // Read in values for variables in the rtparams object.
    int first_light = 1;
    vector<double> tmp;
    for(int i=0;i<(int)token_vec.size();i++) {
      int id=0;
      string tok=token_vec[i];
      string sval;
      double val;
      int found_token=0;
      if(tok=="RESX") {
	id=0;
	val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	rtp.resx=val;
	found_token=1;
      }// RESX
      if(tok=="RESY") {
	id=0;
	val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	rtp.resy=val;
	found_token=1;
      }// RESY
      if(tok=="VIEWDIR") {
	// Note that there are 3 values here.
	int id=0;
	for(int ic=0;ic<3;ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.viewdir[ic]=val;
	}
	found_token=1;
	rtp.viewdir_s=1;
      }// VIEWDIR
      if(tok=="UPDIR") {
	// Note that there are 3 values here.
	int id=0;
	for(int ic=0;ic<3;ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.updir[ic]=val;
	}
	found_token=1;
	rtp.updir_s=1;
      }// UPDIR
      if(tok=="ZOOM") {
	id=0;
	val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	rtp.zoom=val;
	found_token=1;
      }// ZOOM
      if(tok=="ASPECTRATIO") {
	id=0;
	val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	rtp.aspectratio=val;
	found_token=1;
      }// ASPECTRATIO
      if(tok=="ANTIALIASING") {
	id=0;
	val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	rtp.antialiasing=val;
	found_token=1;
      }// ANTIALIASING
      if(tok=="RAYDEPTH") {
	id=0;
	val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	rtp.raydepth=val;
	found_token=1;
      }// RAYDEPTH
      if(tok=="CENTER") {
	// Note that there are 6 values here.
	int id=0;
	for(int ic=0;ic<6;ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.center_guide[ic]=val;
	}
	rtp.center_s=1;
	found_token=1;
      }// CENTER
      if(tok=="LIGHT") {
	/* This is a little confusing.  There can be multiple
	   LIGHT tokens.  Each time one is found a new LIGHT is
	   added to the rtparams object.  Note that there are
	   7 values here.  3 for center, 1 for rad, 3 for color.
	   These values must be pushed back onto the correct vectors
	   and matrices.  However, for the first LIGHT token, one must
	   overwrite the default rather than add a new light.  
	*/
	vector<double> lcen(3);
	double lrad;
	vector<double> lcolor(3);
	int id=0;
	for(int ic=0;ic<3;ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  lcen[ic]=val;
	}
	val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	lrad=val;
	for(int ic=0;ic<3;ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  lcolor[ic]=val;
	}
	if(first_light) {
	  first_light=0;
	  rtp.lightcenter[0]=lcen;
	  rtp.lightrad[0]=lrad;
	  rtp.lightcolor[0]=lcolor;
	}
	else{
	  rtp.lightcenter.push_back(lcen);
	  rtp.lightrad.push_back(lrad);
	  rtp.lightcolor.push_back(lcolor);
	}
	found_token=1;
      }// LIGHT
      if(tok=="BACKGROUND") {
	// Note that there are 3 values here.
	int id=0;
	for(int ic=0;ic<3;ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.background[ic]=val;
	}
	found_token=1;
      }// BACKGROUND
      if(tok=="ATOMCOLOR") {
	// Note that there are 4 values here.
	int id=0;
	tmp = vector<double> (4);
	for(int ic=0;ic<(int)tmp.size();ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  tmp[ic]=val;
	}
	rtp.sphtex_color.push_back(tmp);
	found_token=1;
      }// ATOMCOLOR
      if(tok=="ATOMTEXTURE") {
	// Note that there are 5 values here.
	int id=0;
	tmp = vector<double> (5);
	for(int ic=0;ic<(int)tmp.size();ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  tmp[ic]=val;
	}
	rtp.sphtex_tex.push_back(tmp);
	found_token=1;
      }// ATOMTEXTURE
      if(tok=="ATOMRAD") {
	// Note that there are 2 values here.
	int id=0;
	tmp = vector<double> (2);
	for(int ic=0;ic<(int)tmp.size();ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.sph_rad.push_back(val);
	}
	found_token=1;
      }//ATOMRAD
      if(tok=="SHADING") {
	int id=0;
	rtp.shading = aurostd::GetNextVal(val_vec[i],id);
	found_token=1;
      }// SHADING
      if(tok=="OUTFILE") {
	int id=0;
	rtp.outfile = aurostd::GetNextVal(val_vec[i],id);
	found_token=1;
      }// OUTFILE
      if(tok=="SUPERCELL") {
	// Note that there are 9 values here.
	int id=0;
	for(int ic=0;ic<3;ic++) {
	  for(int jc=0;jc<3;jc++) {
	    val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	    rtp.sc[ic][jc]=val;
	  }
	}
	found_token=1;
	rtp.sc_s=1;
      }//SUPERCELL
      if(tok=="CALCTYPE") {
	int id=0;
	rtp.calc_type=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }//CALCTYPE
      if(tok=="INFILE") {
	int id=0;
	while (id<(int)val_vec[i].size()) {
	  string s;
	  s=aurostd::GetNextVal(val_vec[i],id);
	  if(s.size()>0) rtp.input_files.push_back(s);
	}
	found_token=1;
      }//INSERT_FILE
      if(tok=="INSERT_FILE") {
	int id=0;
	string s;
	s=aurostd::GetNextVal(val_vec[i],id);
	rtp.insert_file=s;
	found_token=1;
      }//INSERT_FILE
      if(tok=="ROTATION") {
	int id=0;
	for(int ic=0;ic<6;ic++) {
	  string s=aurostd::GetNextVal(val_vec[i],id);
	  rtp.rotation[ic] = atof(s.c_str());
	  //rtp.rotation[ic] = atof(aurostd::GetNextVal(val_vec[ic],id).c_str());
	}
	found_token=1;
      }//ROTATION
      if(tok=="STRUCTURE_ORIGIN") {
	int id=0;
	for(int ic=0;ic<3;ic++) {
	  string s=aurostd::GetNextVal(val_vec[i],id);
	  rtp.struct_origin[ic] = atof(s.c_str());
	}
	rtp.struct_origin_s=1;
	found_token=1;
      }//STRUCTURE_ORIGIN
      if(tok=="PLANE") {
	int id=0;
	rtp.plane=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }//PLANE
      if(tok=="PLANECENTER") {
	// Note that there are 3 values here.
	int id=0;
	for(int ic=0;ic<(int)rtp.plane_center.size();ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.plane_center[ic]=val;
	}
	found_token=1;
	rtp.updir_s=1;
      }// PLANECENTER
      if(tok=="PLANENORMAL") {
	// Note that there are 3 values here.
	int id=0;
	for(int ic=0;ic<(int)rtp.plane_normal.size();ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.plane_normal[ic]=val;
	}
	found_token=1;
      }// PLANENORMAL
      if(tok=="PLANETEXTURE") {
	// Note that there are 4 values here.
	int id=0;
	for(int ic=0;ic<(int)rtp.planetex_tex.size();ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.planetex_tex[ic]=val;
	}
	found_token=1;
      }// PLANETEXTURE
      if(tok=="PLANECOLOR") {
	// Note that there are 3 values here.
	int id=0;
	for(int ic=0;ic<(int)rtp.plane_color.size();ic++) {
	  val=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	  rtp.plane_color[ic]=val;
	}
	found_token=1;
      }// PLANECOLOR

      if(!found_token) {
	cerr << "ERROR: You have input a token " << tok <<endl;
	cerr << "ERROR: This token is not recognized. Exiting! " << endl;
	exit(1);
      }

    } // for i
  }// end routine
}

// ***************************************************************************
// ReadInStrVec
// ***************************************************************************
namespace pflow {
  void ReadInStrVec(vector<xstructure>& vstr, ifstream& strlist_inf)  {
    vstr.clear();
    stringstream sss;
    vector<string> vline,vtmp;
    aurostd::stream2vectorstring(strlist_inf,vline);
    uint iline=0;
    while (iline<vline.size()) {
      if(aurostd::RemoveWhiteSpaces(vline.at(iline))=="") {
	xstructure str(sss);
	vstr.push_back(str);
	sss.clear();sss.str("");
      } else {
	sss << vline.at(iline) << endl;
      }
      iline++;
    } 
    // [OBSOLETE]    //  int cnt=0;
    // [OBSOLETE]       xstructure str;
    // [OBSOLETE]    // Read in first structure.
    // [OBSOLETE]    strlist_inf >> str;
    // [OBSOLETE]    //  cnt++;
    // [OBSOLETE]    //  cout << "Read in structure: " << cnt << endl;
    // [OBSOLETE]    vstr.push_back(str);
    // [OBSOLETE]    // Read in all remaining structures.
    // [OBSOLETE]    while (getline(strlist_inf,dum)) {
    // [OBSOLETE]      xstructure str;
    // [OBSOLETE]      strlist_inf >> str;
    // [OBSOLETE]      //    cnt++;
    // [OBSOLETE]      //    cout << "Read in structure: " << cnt << endl;
    // [OBSOLETE]      vstr.push_back(str);
    // [OBSOLETE]    }
  }
}

// ***************************************************************************
// JoinStrVec
// ***************************************************************************
// Note that this can be made more memory efficient
// quite easily if that is needed.
namespace pflow {
  void JoinStrVec(vector<xstructure> vstr1,
		  vector<xstructure> vstr2,
		  vector<xstructure>& vstrtot) {
    // cout << "In JOIN" << endl;
    // cout << "strlist1 size " << vstr1.size() << endl;
    // cout << "strlist2 size " << vstr2.size() << endl;
    int size1=vstr1.size();
    int size2=vstr2.size();
    if(size1<size2) {// Pad vstr1
      xstructure str=vstr1[size1-1];
      for(int is=size1;is<size2;is++) {
	vstr1.push_back(str);
      }
    }
    if(size2<size1) {// Pad vstr2
      xstructure str=vstr2[size2-1];
      for(int is=size2;is<size1;is++) {
	vstr2.push_back(str);
      }
    }
    for(int is=0;is<size1;is++) {
      xstructure str1=vstr1[is];
      xstructure str2=vstr2[is];
      matrix<double> lat1=pflow::GetLat(str1);
      xmatrix<double> xlat1(3,3);xlat1=str1.lattice;
      matrix<double> cpos2=pflow::GetCpos(str2);
      deque<int> num_each_type_1=pflow::GetNumEachType(str1);
      int num_types_1=num_each_type_1.size();
      deque<int> num_each_type_2=pflow::GetNumEachType(str2);
      int num_types_2=num_each_type_2.size();
      // Loop over str2 atoms, assign them a type and c/d positions in an atom.
      int cnt=0;
      for(int it=0;it<num_types_2;it++) {
	for(int ia=0;ia<num_each_type_2.at(it);ia++) {
	  _atom a;
	  a=pflow::SetCpos(a,cpos2[cnt]);
	  a.fpos=C2F(xlat1,a.cpos);
	  //	a=pflow::SetFpos(a,C2F(lat1,cpos2[cnt]));
	  a=pflow::SetType(a,num_types_1+it);
	  str1.AddAtom(a);
	  cnt++;
	}
      }
      //      cout << str1.atoms.size() << endl;
      vstrtot.push_back(str1);
    }//is
  }//end routine
}

// ***************************************************************************
// SetStrFromRTParams
// ***************************************************************************
// This assumes you have read in the rtparams with
// SetRTParams.  It uses both a structure and the rtparams
// and alters the structure based on rtparams.  Here the
// structure is made into a supercell if needed and the
// struc_origin is set.
namespace pflow {
  void SetStrFromRTParams(xstructure& str, pflow::rtparams& rtp) {
    // Make structral adjustments
    if(rtp.sc_s) {
      xmatrix<double> _supercell(3,3);
      _supercell=pflow::matrix2xmatrix(rtp.sc);
      // Make a supercell if it was set by user.  Note that this moves all atoms into the unit cell.
      str=GetSuperCell(str,_supercell);
      // Reset rtparams now that you have a supercell.
      pflow::SetRTParams(str,rtp);
    }
    str=pflow::SetOrigin(str,rtp.struct_origin);
  }
}
// ***************************************************************************
// SuperCellStrVec
// ***************************************************************************
// Gets supercells for every strucutre in the structure vector.
namespace pflow {
  void SuperCellStrVec(vector<xstructure>& vstr, const matrix<double>& sc) {
    if(vstr.size()==0) return;
    for(int is=0;is<(int)vstr.size();is++) {
      xmatrix<double> _supercell(3,3);
      _supercell=pflow::matrix2xmatrix(sc);
      vstr[is]=GetSuperCell(vstr[is],_supercell); // Makes a supercell.
    }
  }
}
// ***************************************************************************
// UpDateRTParams
// *************************************************
// This simply makes changes to rtparams as the program steps
// through different frames.  This is use for values that must
// evolve during the movie (e.g., the center).
namespace pflow {
  void UpDateRTParams(pflow::rtparams& rtp, const int& istr, int nstr) {
    // Update center
    // tpx tpx ??? (cut time to 1/2 movie for moving center by /2).
    //nstr=(int) (nstr/1.5);
    double r=0;
    if(nstr>1) {
      r = (double)istr/(double)(nstr-1);
    }
    if(r>1) r=1;
    if(rtp.center_s) {
      for(int ic=0;ic<3;ic++) {
	rtp.center[ic]=rtp.center_guide[2*ic]+r*(rtp.center_guide[2*ic+1]-rtp.center_guide[2*ic]);
      }
    }
    //tpx
    cout << "CENTER " ;
    pflow::Vout(rtp.center,cout);
  }
}

// ***************************************************************************
// SetRTParams
// ***************************************************************************
// This assumes you have read in the rtparams with
// SetRTParams.  It uses both a structure and the rtparams
// and alters rtparams based on what has been read in and the
// structure.  The rtparam defaults are set based on the
// ntypes and other structural characteristics.
namespace pflow {
  void SetRTParams(xstructure& str, pflow::rtparams& rtp) {
    // Now set all defaults for rtp that
    // use structural information.
    str=ReScale(str,1.0);
    matrix<double> lat=pflow::GetLat(str);
    matrix<double> cpos=pflow::GetCpos(str);
    int nat=cpos.size();

    // For all items that must be set for each type we do that here.
    // These include sphtex_* and sph_rad.
    // Note that we protect against reading in data for atom types that do not exist.
    // Also note that this makes of unique formats for the rtparam
    // data that only exist the when they are read in.  Therefore,
    // we must only do this the first time this set function is
    // called.  This is controlled by the first_set parameter.
    if(rtp.first_set) {
      int ntypes = pflow::GetNumEachType(str).size();
      matrix<double> tmp;
      int size;
      // Set sphtex_tex
      size=rtp.sphtex_tex_def.size();
      tmp=matrix<double> (ntypes,rtp.sphtex_tex_def);
      for(int i=0;i<(int)rtp.sphtex_tex.size();i++) {
	int id=(int)rtp.sphtex_tex[i][0]-1;
	for(int j=1;j<size+1;j++) {
	  if(id>=0 && id<ntypes) tmp[id][j-1]=rtp.sphtex_tex[i][j];
	}
      }
      rtp.sphtex_tex=tmp;
      // Set sphtex_color
      size=rtp.sphtex_color_def.size();
      tmp=matrix<double> (ntypes,rtp.sphtex_color_def);
      for(int i=0;i<(int)rtp.sphtex_color.size();i++) {
	int id=(int)rtp.sphtex_color[i][0]-1;
	for(int j=1;j<size+1;j++) {
	  if(id>=0 && id<ntypes) tmp[id][j-1]=rtp.sphtex_color[i][j];
	}
      }
      rtp.sphtex_color=tmp;
      // Set sphtex_phong
      size=rtp.sphtex_phong_def.size();
      tmp=matrix<double> (ntypes,rtp.sphtex_phong_def);
      for(int i=0;i<(int)rtp.sphtex_phong.size();i++) {
	int id=(int)rtp.sphtex_phong[i][0]-1;
	for(int j=1;j<size+1;j++) {
	  if(id>=0 && id<ntypes) tmp[id][j-1]=rtp.sphtex_phong[i][j];
	}
      }
      rtp.sphtex_phong=tmp;
      // Set sph_rad
      size=1;
      vector<double> vtmp(ntypes,rtp.sph_rad_def);
      for(int i=0;i<(int)rtp.sph_rad.size()/2;i++) {
	int id=(int)rtp.sph_rad[2*i]-1;
	if(id>=0 && id<ntypes) vtmp[id]=rtp.sph_rad[2*i+1];
      }
      rtp.sph_rad=vtmp;
      // Set sphtex_names
      rtp.sphtex_names = vector<string> (ntypes,"txt_atom_type_");
      for(int it=0;it<ntypes;it++) {
	ostringstream os;
	os << it+1 << ends;
	string s(os.str());
	rtp.sphtex_names[it]=rtp.sphtex_names[it]+s;
      }
      rtp.first_set=0;
    }// if first_set    

    //  Set view direction along the 0 axis.
    //  if(!rtp.viewdir_s) {rtp.viewdir=lat[0];}

    //  Set up direction along the 1 axis.
    //  if(!rtp.updir_s) {rtp.updir=lat[0];}

    /*
      Set up struct_origin, which is the point around which rotation
      will take place.  Here we default to the first moment of the atom
      positions.
    */
    if(! rtp.struct_origin_s) {// If struct_orig was not set by user
      // rtp.struct_origin=SVprod(0.5,VVsum(lat[0],VVsum(lat[1],lat[2])));
      vector<double> mom1=xvector2vector(GetMom1(str));
      rtp.struct_origin=mom1;
    }

    /*
      Set up center (where camera is located).
      Define the center to be displaced from the centroid along the view
      direction. Define center to be located a distance d=-2.0*R/tan(45)=-2.0*R
      from the centroid, where R is the the largest projected distance
      of an atom from the centroid (first moment) in the image plane (the
      plane defined to be perpendicular to the view direction and containing the
      centroid).  So define
      mom1=first moment
      vd=unit vector along view direction
      up=unit vector along up direction
      upp=unit vector consisting of projection of up vector into the image plane.
      vdp=unit vector in image plane perpendicular to vd and updir.

      Actually, instead of mom1 being the centroid just set it to the
      struct_origin.  All else is the same.
    */
    if(! rtp.center_s) {// If center was not set by user
      vector<double> mom1=rtp.struct_origin;
      vector<double> vd=pflow::SVprod(1.0/pflow::norm(rtp.viewdir),rtp.viewdir);
      vector<double> up=pflow::SVprod(1.0/pflow::norm(rtp.updir),rtp.updir);
      vector<double> upp;
      upp=pflow::VVsum(up,pflow::SVprod(-pflow::VVprod(up,vd),vd));
      upp=pflow::SVprod(1.0/pflow::norm(upp),upp);
      vector<double> vdp=pflow::VVcross(vd,upp);
      double proj=0;
      //      double max_iat=0;
      for(int iat=0;iat<nat;iat++) {
	double p1=pflow::VVprod(pflow::VVdiff(cpos[iat],mom1),vdp);
	double p2=pflow::VVprod(pflow::VVdiff(cpos[iat],mom1),upp);
	double nproj=sqrt(p1*p1+p2*p2);
	if(nproj>proj) {
	  proj=nproj;
	  //	  max_iat=iat;
	}
      }
      double displ=-2.0*proj;
      rtp.center=pflow::VVsum(mom1,pflow::SVprod(displ,vd));
    }// if ! center_s
    else{ // if center is set by user
      for(int ic=0;ic<3;ic++) {
	rtp.center[ic]=rtp.center_guide[2*ic];
      }
    }
  }
}

// ***************************************************************************
// GetRTDatFile
// ***************************************************************************
namespace pflow {
  void GetRTDatFile(xstructure str, const pflow::rtparams& rtp,
		    ostringstream& rtdat_file) {
    rtdat_file << "BEGIN_SCENE" << endl;
    rtdat_file << "  RESOLUTION " << rtp.resx << " " << rtp.resy << endl;
    /*
      string tga="tga";
      string outfile;
      outfile=rtp.outfile+"_"+tga;
      rtdat_file << "  OUTFILE " << outfile << endl;
    */
    rtdat_file << endl;
    rtdat_file << "CAMERA" << endl;
    rtdat_file << "  ZOOM " << rtp.zoom << endl;
    rtdat_file << "  ASPECTRATIO " << rtp.aspectratio <<  endl;
    rtdat_file << "  ANTIALIASING " << rtp.antialiasing << endl;
    rtdat_file << "  RAYDEPTH " << rtp.raydepth << endl;
    rtdat_file << "  CENTER " << rtp.center[0] << " " << rtp.center[1] << " " << rtp.center[2] << endl;
    rtdat_file << "  VIEWDIR " << rtp.viewdir[0] << " " << rtp.viewdir[1] << " " << rtp.viewdir[2] << endl;
    rtdat_file << "  UPDIR " << rtp.updir[0] << " " << rtp.updir[1] << " " << rtp.updir[2] << endl;
    rtdat_file <<  endl;
    rtdat_file << "END_CAMERA" << endl;
    rtdat_file <<  endl;
    rtdat_file << "BACKGROUND " << rtp.background[0] << " " << rtp.background[1] << " " << rtp.background[2] << endl;
    rtdat_file <<  endl;
    for(int il=0;il<(int)rtp.lightcenter.size();il++) {
      rtdat_file << "LIGHT ";
      rtdat_file << "CENTER " << rtp.lightcenter[il][0] << " " << rtp.lightcenter[il][1] << " " << rtp.lightcenter[il][2];
      rtdat_file << " RAD " << rtp.lightrad[il];
      rtdat_file << " COLOR " << rtp.lightcolor[il][0] << " " << rtp.lightcolor[il][1] << " " << rtp.lightcolor[il][2];
      rtdat_file <<  endl << endl;
    }
    for(int iat=0;iat<(int) rtp.sphtex_names.size();iat++) {
      rtdat_file << "TEXDEF " << rtp.sphtex_names[iat] << " ";
      rtdat_file << "AMBIENT " << rtp.sphtex_tex[iat][0] << " ";
      rtdat_file << "DIFFUSE " << rtp.sphtex_tex[iat][1] << " ";
      rtdat_file << "SPECULAR " << rtp.sphtex_tex[iat][2] << " ";
      rtdat_file << "OPACITY " << rtp.sphtex_tex[iat][3] << " ";
      rtdat_file << endl;
      rtdat_file << "PHONG PLASTIC " << rtp.sphtex_phong[iat][0] << " " << "PHONG_SIZE " << rtp.sphtex_phong[iat][1] << endl;
      rtdat_file << "  COLOR " << rtp.sphtex_color[iat][0] << " " << rtp.sphtex_color[iat][1] << " " << rtp.sphtex_color[iat][2] << endl;
      rtdat_file << "  TEXFUNC 0" << endl;
      rtdat_file << endl;
    }
    // Do atoms
    str=ReScale(str,1.0);
    matrix<double> cpos=pflow::GetCpos(str);
    deque<int> num_each_type=pflow::GetNumEachType(str);
    //int nat=(int)cpos.size();
    int ntype=(int)num_each_type.size();
    int cnt=0;
    for(int it=0;it<ntype;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
	rtdat_file << "SPHERE CENTER " << cpos[cnt][0] << " " << cpos[cnt][1] << " " << cpos[cnt][2];
	rtdat_file << " RAD " << rtp.sph_rad[it];
	rtdat_file << " " << rtp.sphtex_names[it] << endl;
	rtdat_file << endl;
	cnt++;
      }
    }
    // Do plane
    if(rtp.plane==1) {
      rtdat_file << "PLANE" << endl;
      rtdat_file << "  CENTER " << rtp.plane_center[0] << " " << rtp.plane_center[1] << " " << rtp.plane_center[2] << endl;
      rtdat_file << "  NORMAL " << rtp.plane_center[0] << " " << rtp.plane_center[1] << " " << rtp.plane_center[2] << endl;
      rtdat_file << "  TEXTURE " << endl;
      rtdat_file << "    AMBIENT " << rtp.planetex_tex[0] << " ";
      rtdat_file << "DIFFUSE " << rtp.planetex_tex[1] << " ";
      rtdat_file << "SPECULAR " << rtp.planetex_tex[2] << " ";
      rtdat_file << "OPACITY " << rtp.planetex_tex[3] << " ";
      rtdat_file << "  COLOR " << rtp.plane_color[0] << " " << rtp.plane_color[1] << " " << rtp.plane_color[2] << endl;
      rtdat_file << "  TEXFUNC 0" << endl;
      cout << endl;
    }

    rtdat_file << "END_SCENE" << endl;
    rtdat_file << ends;
  }
}

// ***************************************************************************
// PrintRDatFile
// ***************************************************************************
namespace pflow {
  string PrintRTDatFile(ostringstream& rtdat_file, const pflow::rtparams& rt_params) {  
    string filename = rt_params.outfile+".dat";
    ostringstream outss;
    outss << filename << ends;
    ofstream outf(outss.str().c_str());
    outf << rtdat_file.str();
    return filename;
  }
}

// ***************************************************************************
// CreateRTtgaFile
// ***************************************************************************
namespace pflow {
  string CreateRTtgaFile(const string& datfile, const pflow::rtparams& rt_params) {  
    // Run tachyon on datfile.
    ostringstream tachyon_stream;
    tachyon_stream << "tachyon " << "-" << rt_params.shading << " " << datfile << ends;
    // system(tachyon_stream.str().c_str());
    aurostd::execute(tachyon_stream);
    // Move outfile.tga to user defined name.
    ostringstream mv_stream;
    string tgafile = rt_params.outfile+".tga";
    mv_stream << "mv outfile.tga " << tgafile << ends;
    aurostd::execute(mv_stream);
    //  system(mv_stream.str().c_str());
    return tgafile;
  }
}

// ***************************************************************************
// CreateRTjpgFile
// ***************************************************************************
namespace pflow {
  string CreateRTjpgFile(const string& tgafile, const pflow::rtparams& rt_params) {  
    // Run convert.
    ostringstream convert_stream;
    string jpgfile = rt_params.outfile+".jpg";
    convert_stream << XHOST.command("convert") << " " << tgafile << " " << jpgfile << ends;
    aurostd::execute(convert_stream);
    // system(convert_stream.str().c_str());
    return jpgfile;
  }
}

// ***************************************************************************
// GetRTencFile
// ***************************************************************************
namespace pflow {
  void GetRTencFile(const pflow::rtparams& rtp, const int nim,
		    ostringstream& os) {
    os << "PATTERN         IBBPBBPBBPBBPBB" << endl;
    os << "OUTPUT          "<< rtp.outfile << ".mpg" << endl;
    os << "GOP_SIZE        16" << endl;
    os << "SLICES_PER_FRAME        1" << endl;
    os << "BASE_FILE_FORMAT        PPM" << endl;
    os << endl;
    os << "INPUT_CONVERT   djpeg *" << endl;
    os << "INPUT_DIR       ." << endl;
    os << "INPUT" << endl;
    os << rtp.outfile << "_*.jpg       [" << aurostd::PaddedNumString(0,5) << "-" << nim-1 << "+1]" << endl;
    os << "END_INPUT" << endl;
    os << endl;
    os << "IQSCALE         8" << endl;
    os << "PQSCALE         10" << endl;
    os << "BQSCALE         25" << endl;
    os << "PIXEL           HALF" << endl;
    os << "RANGE           10" << endl;
    os << "PSEARCH_ALG     LOGARITHMIC" << endl;
    os << "BSEARCH_ALG     CROSS2" << endl;
    os << "REFERENCE_FRAME DECODED" << endl;
    os << ends;
  }
}

// ***************************************************************************
// PrintRTencFile
// ***************************************************************************
namespace pflow {
  string PrintRTencFile(const pflow::rtparams& rt_params, ostringstream& rtenc_file) {
    string filename = rt_params.outfile+".enc";
    ostringstream outss;
    outss << filename << ends;
    ofstream outf(outss.str().c_str());
    outf << rtenc_file.str();
    return filename;
  }
}

// ***************************************************************************
// CreateRTmgpFile
// ***************************************************************************
namespace pflow {
  string CreateRTmpgFile(const pflow::rtparams& rt_params, const string& encfile) {
    // Run mpeg_encoder on enc_file.
    ostringstream encoder_stream;
    encoder_stream << "mpeg_encode " << " " << encfile << ends;
    /*  char* encoder_char = encoder_stream.str().c_str();
	system(encoder_char);*/
    system(encoder_stream.str().c_str());
    // Get outfile.mpg name
    string file;
    file = rt_params.outfile+".mpg";
    return file;
  }
}

// ***************************************************************************
// RayTraceManager
// ***************************************************************************
namespace pflow {
  void RayTraceManager(vector<string> argv) {
    ifstream rtinfile(argv.at(2).c_str()); // File where RT params are input.
    aurostd::InFileExistCheck("RayTraceFuncs.cc/RayTraceManager",argv.at(2),rtinfile,cerr);
    pflow::rtparams rtp; // Object that stores RT params.
    ReadInRTParams(rtinfile,rtp); // Sets RT params object from input.
    switch (rtp.calc_type) {
    case 0:{ // A single file of structures.
      // Read in structure list.
      if(rtp.input_files.size()<1) {
	cerr << "ERROR: RayTraceManager" << endl;
	cerr << "ERROR: You must specify a structure list file to open with the token INFILE." << endl;
	cerr << "ERROR: For example:  INFILE = STRLIST ." << endl;
	cerr << "ERROR: Exiting." << endl;
	exit(1);
      }
      ifstream strlist_inf(rtp.input_files[0].c_str());
      aurostd::InFileExistCheck("RayTraceFuncs.cc/RayTraceManager",rtp.input_files[0].c_str(),strlist_inf,cerr);
      vector<xstructure> vstr;
      cout << endl;
      cout << "Reading in structure list." << endl;
      cout << endl;
      pflow::ReadInStrVec(vstr,strlist_inf);
      int nstr=vstr.size();
      cout << "We are working with this many structures: " << nstr << endl;
      // Set rtparams.
      pflow::SetRTParams(vstr[0],rtp);
      // Set structures from rtparams.
      for(int is=0;is<nstr;is++) {
	pflow::SetStrFromRTParams(vstr[is],rtp);
      }
      cout << endl;
      cout << "Here are the values of the ray trace parameters." << endl;
      rtp.Print();
      cout << endl;
      // Rotate vstr
      RotateStrVec(vstr,rtp.rotation);
      // Create dat,tga,jpg files.
      string datfile,tgafile,jpgfile;
      for(int is=0;is<nstr;is++) {
	cout << endl;
	cout << "Creating dat/tga/jpg file num " << is+1 << " out of a total of " << nstr << endl;
	cout << endl;
	// Update rtparams
	pflow::UpDateRTParams(rtp,is,nstr);
	// Get dat file stringstream.
	ostringstream rtdat_file;
	xstructure str=vstr[is];
	pflow::GetRTDatFile(str,rtp,rtdat_file);       // Puts formatted .dat file into a stringstream.
	datfile=pflow::PrintRTDatFile(rtdat_file,rtp); // Prints outfile.dat file and returns name.
	tgafile=pflow::CreateRTtgaFile(datfile,rtp);   // Creates outfile.tga file and returns name.
	jpgfile=pflow::CreateRTjpgFile(tgafile,rtp);   // Creates outfile.jpg file and returns name.
	// Copy files to proper names (outfile.num.suffix).
	ostringstream mvcmd_dat;
	ostringstream mvcmd_tga;
	ostringstream mvcmd_jpg;
	string name;
	string num = aurostd::PaddedNumString(is,5);
	string datname = (rtp.outfile+"_"+num+".dat");
	mvcmd_dat << "/bin/mv " << datfile << " " << datname << ends;
	system(mvcmd_dat.str().c_str());
	string tganame = (rtp.outfile+"_"+num+".tga");
	mvcmd_tga << "/bin/mv " << tgafile << " " << tganame << ends;
	system(mvcmd_tga.str().c_str());
	string jpgname = (rtp.outfile+"_"+num+".jpg");
	mvcmd_jpg << "/bin/mv " << jpgfile << " " << jpgname << ends;
	system(mvcmd_jpg.str().c_str());
	if(nstr>1) { // Remove dat,tga files if there are a lot of them.
	  ostringstream rmcmd_dat;
	  ostringstream rmcmd_tga;
	  rmcmd_dat << "/bin/rm " << datname << ends;
	  system(rmcmd_dat.str().c_str());
	  //	rmcmd_tga << "/bin/rm " << tganame << ends;
	  //	system(rmcmd_tga.str());
	}
      }
      // Get and print enc file for mpeg encoder.
      ostringstream rtenc_file;
      GetRTencFile(rtp,nstr,rtenc_file);
      string encfile;
      encfile=PrintRTencFile(rtp,rtenc_file);
      // Create mpg file with mpeg encoder.
      cout << endl;
      cout << "Creating mpg file " << rtp.outfile << ".mpg" << endl;
      cout << endl;
      string mpgfile;
      mpgfile=CreateRTmpgFile(rtp,encfile);
      break;
    }
    default:{
      cout << endl;
      cerr << "ERROR: RTManager" <<endl;
      cerr << "ERROR: You have set CALCTYPE= " << rtp.calc_type << endl;
      cerr << "ERROR: This does not correspond to any function. "<< endl;
      cerr << "ERROR: exiting" << endl;
      cout << endl;
      exit(1);
    }
    }//switch
  }
}

// ***************************************************************************
// GetRotationMatrix
// ***************************************************************************
// This gets a rotation matrix from 3 angles assumed
// to represent a rotation around x, then y, then z.
// Angles are assumed to be in radians.
namespace pflow {
  matrix<double> GetRotationMatrix(const vector<double>& angles) {
    // Sin and cos.
    vector<double> sn(3,0.0);
    vector<double> cs(3,0.0);
    for(int ic=0;ic<3;ic++) {
      sn[ic]=sin(angles[ic]);
      cs[ic]=cos(angles[ic]);
    }
    // Set rotation matrix (do x, then y, then z rotation).
    matrix<double> xm(3,3);pflow::VVset(xm,0.0);
    matrix<double> ym(3,3);pflow::VVset(ym,0.0);
    matrix<double> zm(3,3);pflow::VVset(zm,0.0);

    xm[0][0]=1;
    xm[1][1]=cs[0];
    xm[1][2]=-sn[0];
    xm[2][1]=sn[0];
    xm[2][2]=cs[0];

    ym[0][0]=cs[1];
    ym[0][2]=sn[1];
    ym[1][1]=1;
    ym[2][0]=-sn[1];
    ym[2][2]=cs[1];

    zm[0][0]=cs[2];
    zm[0][1]=-sn[2];
    zm[1][0]=sn[2];
    zm[1][1]=cs[2];
    zm[2][2]=1;

    matrix<double> rm;
    rm=pflow::MMmult(zm,pflow::MMmult(ym,xm));
    return rm;
  }
}

// ***************************************************************************
// RotateStrVec
// ***************************************************************************
// This rotates each structure in the vstr.
// The rotation goes from an initial to a final set of angles
// in steps of the primitive rotation.
// The primitive rotation is the rotation around x, y, then z
// by the amount of change given in rot divided by the number
// of structures - 1.  All frames are initially rotated
// according to the initial rotation angles.  The rotation is
// then done as one primitive rotation per frame.
namespace pflow {
  void RotateStrVec(vector<xstructure>& vstr, const vector<double>& rot) {
    // Get initial rotation matrix.
    vector<double> angles(3);
    for(int ic=0;ic<3;ic++) {
      angles[ic]=rot[2*ic];
      angles[ic]=angles[ic]*TWOPI/360.0;
    }
    matrix<double> irm=GetRotationMatrix(angles);
  
    // get primitive rotation matrix.
    int s=vstr.size()-1;
    if(s<1) s=1;
    for(int ic=0;ic<3;ic++) {
      angles[ic]=(rot[2*ic+1]-rot[2*ic])/(s);
      angles[ic]=angles[ic]*TWOPI/360.0;
    }
    matrix<double> prm=GetRotationMatrix(angles);
    matrix<double> rm=irm;
    xmatrix<double> xprm(3,3); xprm=pflow::matrix2xmatrix(prm);
    xmatrix<double> xrm(3,3);  xrm=pflow::matrix2xmatrix(rm);
    for(int is=0;is<(int)vstr.size();is++) {
      //    xrm=pflow::matrix2xmatrix(rm);
      //   vstr[is]=Rotate(vstr[is],xrm);
      vstr[is]=Rotate(vstr[is],pflow::matrix2xmatrix(rm));
      rm=pflow::MMmult(prm,rm);
    }
  }
}

// ***************************************************************************
// PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PRO
// ***************************************************************************
namespace pflow {
  // Constructors
  projdata::projdata()
  {
    sp=0;
    rspin=2;
    nl_max=4; // for s,p,d,f orbitals
    nlm_max=16; // for s,p,d,f orbitals
    nlmtot_max=20; // for s,p,d,f orbitals + p,d,f,all totals
    nl=nl_max; // Default to max values
    nlm=nlm_max; // Default to max values
    nlmtot=nlmtot_max; // Default to max values
    LMnames = vector<string> (nlm);
    LMnames[0]="S     ";
    LMnames[1]="Py    ";
    LMnames[2]="Pz    ";
    LMnames[3]="Px    ";
    LMnames[4]="Dxy   ";
    LMnames[5]="Dyz   ";
    LMnames[6]="Dz2   ";
    LMnames[7]="Dxz   ";
    LMnames[8]="Dx2-y2";
    LMnames[9]="F1    ";
    LMnames[10]="F2    ";
    LMnames[11]="F3    ";
    LMnames[12]="F4    ";
    LMnames[13]="F5    ";
    LMnames[14]="F6    ";
    LMnames[15]="F7    ";
    Lnames = vector<string> (nl);
    Lnames[0]="S     ";
    Lnames[1]="P     ";
    Lnames[2]="D     ";
    Lnames[3]="F     ";
    LLMnames = vector<string> (nlm+4);
    LLMnames[0]="S     ";
    LLMnames[1]="Py    ";
    LLMnames[2]="Pz    ";
    LLMnames[3]="Px    ";
    LLMnames[4]="Ptot  ";
    LLMnames[5]="Dxy   ";
    LLMnames[6]="Dyz   ";
    LLMnames[7]="Dz2   ";
    LLMnames[8]="Dxz   ";
    LLMnames[9]="Dx2-y2";
    LLMnames[10]="Dtot  ";
    LLMnames[11]="F1    ";
    LLMnames[12]="F2    ";
    LLMnames[13]="F3    ";
    LLMnames[14]="F4    ";
    LLMnames[15]="F5    ";
    LLMnames[16]="F6    ";
    LLMnames[17]="F7    ";
    LLMnames[18]="Ftot  ";
    LLMnames[19]="Tot   ";
  }
  
  void projdata::Print(ostream& outf) {
    outf.precision(3);
    outf.setf(std::ios::fixed,std::ios::floatfield);
    outf.setf(std::ios::left,std::ios::adjustfield);
  
    int ioncnt;
    int w1=5;
    int w2=6;
    string key;

    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "WARNING:  There is some confusing business about these projections." << endl;
    outf << "WARNING:  To get results that match the vasp summations you may need" << endl;
    outf << "WARNING:  to alter the vasp source code.  See the convasp help by" << endl;
    outf << "WARNING:  typing convasp -h .  Read the section on -pocc." << endl;
    outf << endl;
    outf << "The basic data" << endl;
    outf << "nlm " << nlm << endl;
    outf << "nkpts " << nkpts << endl;
    outf << "nbands " << nbands << endl;
    outf << "nions " << nions << endl;
    outf << "ntypes " << ntypes << endl;
    outf << "type  " << "num of that type " << endl;
    for(int it=0;it<ntypes;it++) {
      outf << it << "         " << num_each_type.at(it) << endl;
    }
    outf << "Spin polarization (0=no,1=yes) " << sp << endl;
    outf << "rspin (spin degeneracy) " << rspin << endl;
    outf << "Fermi Weights (nbands X nkpts)" << endl;
    outf << "Band  kpt1  kpt2 ..." << endl;
    outf << "Up spin" << endl;
    for(int ib=0;ib<(int)wfermi_u.size();ib++) {
      outf << " " << ib;
      for(int ik=0;ik<(int)wfermi_u[ib].size();ik++) {
	outf << "   " << wfermi_u[ib][ik];
      }
      outf << endl;
    }
    outf << "Down spin" << endl;
    for(int ib=0;ib<(int)wfermi_d.size();ib++) {
      outf << " " << ib;
      for(int ik=0;ik<(int)wfermi_d[ib].size();ik++) {
	outf << "   " << wfermi_d[ib][ik];
      }
      outf << endl;
    }

    outf << "kpt values and weights" << endl;
    for(int ik=0;ik<(int)wkpt.size();ik++) {
      outf << "  " << ik+1 << "  "
	   << kpts[ik][0] << "  "
	   << kpts[ik][1] << "  "
	   << kpts[ik][2] << "      "
	   << wkpt[ik] << endl;
    }

    outf << "Band energies" << endl;
    outf << "   " << "kpt   " << "band  " << "energy (up/dn) " << endl;
    for(int ik=0;ik<(int)ener_k_b_u.size();ik++) {
      for(int ib=0;ib<(int)ener_k_b_u[ik].size();ib++) {
	outf << "   " << setw(5) << ik+1 << " " << setw(5) << ib+1 << " " << setprecision(5) << ener_k_b_u[ik][ib];
	if(sp) outf << "  " << setprecision(5) << ener_k_b_d[ik][ib];
	outf << endl;
      }
    }

    outf.precision(3);

    outf << "Projections" << endl;
    outf << "   " << "type  " << "ion   " << "kpt   " << "band  " << "lm    " << "projection (up/dn)" << endl;
    ioncnt=-1;
    for(int it=0;it<(int)pdat_u.size();it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
	ioncnt++;
	for(int ik=0;ik<(int)pdat_u[it].size();ik++) {
	  for(int ib=0;ib<(int)pdat_u[it][ik][iat].size();ib++) {
	    for(int ilm=0;ilm<(int)pdat_u[it][ik][iat][ib].size();ilm++) {
	      outf << "   ";
	      outf << setw(w2) << it+1;
	      outf << setw(w2) << ioncnt+1;
	      outf << setw(w2) << ik+1;
	      outf << setw(w2) << ib+1;
	      outf << "   " << LMnames[ilm];
	      outf << "  " << setw(w2) << pdat_u[it][ik][iat][ib][ilm];
	      if(sp) outf << "  " << setw(w2) << pdat_d[it][ik][iat][ib][ilm];
	      outf << endl;
	    }
	  }
	}
      }
    }

    // Occupations ION:KPT:BND:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:KPT:BND:LM" << endl;
    outf << endl;
    key="ION:KPT:BND:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  " << "kpt   " << "band  ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
	ioncnt++;
	for(int ik=0;ik<nkpts;ik++) {
	  for(int ib=0;ib<nbands;ib++) {
	    outf << "   " << setw(w1) << ioncnt+1;
	    outf << " " << setw(w1) << it+1;
	    outf << " " << setw(w1) << ik+1;
	    outf << " " << setw(w1) << ib+1;
	    outf << " " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][0];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][1];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][2];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][3];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][1];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][4];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][5];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][6];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][7];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][8];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][2];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][9];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][10];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][11];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][12];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][13];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][14];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][15];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][3];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][4];
	    outf << " " << key;
	    outf << endl;
	  }
	}
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  for(int ik=0;ik<nkpts;ik++) {
	    for(int ib=0;ib<nbands;ib++) {
	      outf << "   " << setw(w1) << ioncnt+1;
	      outf << " " << setw(w1) << it+1;
	      outf << " " << setw(w1) << ik+1;
	      outf << " " << setw(w1) << ib+1;
	      outf << " " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][0];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][1];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][2];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][3];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][1];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][4];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][5];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][6];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][7];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][8];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][2];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][9];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][10];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][11];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][12];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][13];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][14];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][15];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][3];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][4];
	      outf << " " << key;
	      outf << endl;
	    }
	  }
	}
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  for(int ik=0;ik<nkpts;ik++) {
	    for(int ib=0;ib<nbands;ib++) {
	      outf << "   " << setw(w1) << ioncnt+1;
	      outf << " " << setw(w1) << it+1;
	      outf << " " << setw(w1) << ik+1;
	      outf << " " << setw(w1) << ib+1;
	      outf << " " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][0]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][0];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][1]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][1];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][2]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][2];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][3]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][3];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][1]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][1];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][4]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][4];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][5]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][5];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][6]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][6];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][7]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][7];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][8]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][8];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][2]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][2];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][9]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][9];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][10]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][10];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][11]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][11];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][12]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][12];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][13]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][13];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][14]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][14];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][15]-occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][15];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][3]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][3];
	      outf << "  " << setw(w2) << occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][4]-occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][4];
	      outf << " " << key;
	      outf << endl;
	    }
	  }
	}
      }
    }

    // Occupations ION:BND:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:BND:LM" << endl;
    outf << endl;
    key="ION:BND:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  " << "band  ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
	ioncnt++;
	for(int ib=0;ib<nbands;ib++) {
	  outf << "   " << setw(w1) << ioncnt+1;
	  outf << " " << setw(w1) << it+1;
	  outf << " " << setw(w1) << ib+1;
	  outf << " " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][0];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][1];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][2];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][3];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][1];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][4];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][5];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][6];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][7];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][8];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][2];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][9];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][10];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][11];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][12];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][13];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][14];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][15];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][3];
	  outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][4];
	  outf << " " << key;
	  outf << endl;
	}
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  for(int ib=0;ib<nbands;ib++) {
	    outf << "   " << setw(w1) << ioncnt+1;
	    outf << " " << setw(w1) << it+1;
	    outf << " " << setw(w1) << ib+1;
	    outf << " " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][0];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][1];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][2];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][3];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][1];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][4];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][5];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][6];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][7];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][8];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][2];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][9];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][10];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][11];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][12];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][13];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][14];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_d[ioncnt][ib][15];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][3];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_d[ioncnt][ib][4];
	    outf << " " << key;
	    outf << endl;
	  }
	}
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "band  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  for(int ib=0;ib<nbands;ib++) {
	    outf << "   " << setw(w1) << ioncnt+1;
	    outf << " " << setw(w1) << it+1;
	    outf << " " << setw(w1) << ib+1;
	    outf << " " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][0]-occ_vs_ion_bnd_lm_d[ioncnt][ib][0];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][1]-occ_vs_ion_bnd_lm_d[ioncnt][ib][1];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][2]-occ_vs_ion_bnd_lm_d[ioncnt][ib][2];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][3]-occ_vs_ion_bnd_lm_d[ioncnt][ib][3];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][1]-occ_vs_ion_bnd_l_d[ioncnt][ib][1];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][4]-occ_vs_ion_bnd_lm_d[ioncnt][ib][4];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][5]-occ_vs_ion_bnd_lm_d[ioncnt][ib][5];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][6]-occ_vs_ion_bnd_lm_d[ioncnt][ib][6];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][7]-occ_vs_ion_bnd_lm_d[ioncnt][ib][7];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][8]-occ_vs_ion_bnd_lm_d[ioncnt][ib][8];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][2]-occ_vs_ion_bnd_l_d[ioncnt][ib][2];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][9]-occ_vs_ion_bnd_lm_d[ioncnt][ib][9];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][10]-occ_vs_ion_bnd_lm_d[ioncnt][ib][10];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][11]-occ_vs_ion_bnd_lm_d[ioncnt][ib][11];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][12]-occ_vs_ion_bnd_lm_d[ioncnt][ib][12];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][13]-occ_vs_ion_bnd_lm_d[ioncnt][ib][13];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][14]-occ_vs_ion_bnd_lm_d[ioncnt][ib][14];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_lm_u[ioncnt][ib][15]-occ_vs_ion_bnd_lm_d[ioncnt][ib][15];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][3]-occ_vs_ion_bnd_l_d[ioncnt][ib][3];
	    outf << "  " << setw(w2) << occ_vs_ion_bnd_l_u[ioncnt][ib][4]-occ_vs_ion_bnd_l_d[ioncnt][ib][4];
	    outf << " " << key;
	    outf << endl;
	  }
	}
      }
    }// if sp

    // Occupations ION:KPT:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:KPT:LM" << endl;
    outf << endl;
    key="ION:KPT:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  " << "kpt   ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
	ioncnt++;
	for(int ik=0;ik<nkpts;ik++) {
	  outf << "   " << setw(w1) << ioncnt+1;
	  outf << " " << setw(w1) << it+1;
	  outf << " " << setw(w1) << ik+1;
	  outf << " " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][0];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][1];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][2];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][3];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][1];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][4];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][5];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][6];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][7];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][8];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][2];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][9];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][10];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][11];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][12];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][13];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][14];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][15];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][3];
	  outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][4];
	  outf << " " << key;
	  outf << endl;
	}
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  for(int ik=0;ik<nkpts;ik++) {
	    outf << "   " << setw(w1) << ioncnt+1;
	    outf << " " << setw(w1) << it+1;
	    outf << " " << setw(w1) << ik+1;
	    outf << " " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][0];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][1];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][2];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][3];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][1];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][4];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][5];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][6];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][7];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][8];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][2];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][9];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][10];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][11];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][12];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][13];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][14];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_d[ioncnt][ik][15];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][3];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_d[ioncnt][ik][4];
	    outf << " " << key;
	    outf << endl;
	  }
	}
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  " << "kpt   ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  for(int ik=0;ik<nkpts;ik++) {
	    outf << "   " << setw(5) << ioncnt+1;
	    outf << " " << setw(5) << it+1;
	    outf << " " << setw(5) << ik+1;
	    outf << " " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][0]-occ_vs_ion_kpt_lm_d[ioncnt][ik][0];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][1]-occ_vs_ion_kpt_lm_d[ioncnt][ik][1];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][2]-occ_vs_ion_kpt_lm_d[ioncnt][ik][2];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][3]-occ_vs_ion_kpt_lm_d[ioncnt][ik][3];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][1]-occ_vs_ion_kpt_l_d[ioncnt][ik][1];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][4]-occ_vs_ion_kpt_lm_d[ioncnt][ik][4];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][5]-occ_vs_ion_kpt_lm_d[ioncnt][ik][5];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][6]-occ_vs_ion_kpt_lm_d[ioncnt][ik][6];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][7]-occ_vs_ion_kpt_lm_d[ioncnt][ik][7];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][8]-occ_vs_ion_kpt_lm_d[ioncnt][ik][8];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][2]-occ_vs_ion_kpt_l_d[ioncnt][ik][2];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][9]-occ_vs_ion_kpt_lm_d[ioncnt][ik][9];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][10]-occ_vs_ion_kpt_lm_d[ioncnt][ik][10];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][11]-occ_vs_ion_kpt_lm_d[ioncnt][ik][11];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][12]-occ_vs_ion_kpt_lm_d[ioncnt][ik][12];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][13]-occ_vs_ion_kpt_lm_d[ioncnt][ik][13];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][14]-occ_vs_ion_kpt_lm_d[ioncnt][ik][14];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_lm_u[ioncnt][ik][15]-occ_vs_ion_kpt_lm_d[ioncnt][ik][15];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][3]-occ_vs_ion_kpt_l_d[ioncnt][ik][3];
	    outf << "  " << setw(w2) << occ_vs_ion_kpt_l_u[ioncnt][ik][4]-occ_vs_ion_kpt_l_d[ioncnt][ik][4];
	    outf << " " << key;
	    outf << endl;
	  }
	}
      }
    }

    // Occupations ION:LM
    outf << endl;
    outf << "**************************************************" << endl;
    outf << endl;
    outf << "Occupations vs. ION:LM" << endl;
    outf << endl;
    key="ION:LM";
    outf << "Up spin"<< endl;
    outf << "   " << "ion   " << "type  ";
    for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
      outf << LLMnames[ilm] << "  ";
    }
    outf << endl;
    ioncnt=-1;
    for(int it=0;it<ntypes;it++) {
      for(int iat=0;iat<num_each_type.at(it);iat++) {
	ioncnt++;
	outf << "   " << setw(w1) << ioncnt+1;
	outf << " " << setw(w1) << it+1;
	outf << " " << setw(w2) << occ_vs_ion_lm_u[ioncnt][0];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][1];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][2];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][3];
	outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][1];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][4];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][5];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][6];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][7];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][8];
	outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][2];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][9];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][10];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][11];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][12];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][13];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][14];
	outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][15];
	outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][3];
	outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][4];
	outf << " " << key;
	outf << endl;
      }
    }
    if(sp) {
      outf << "Down spin"<< endl;
      outf << "   " << "ion   " << "type  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  outf << "   " << setw(w1) << ioncnt+1;
	  outf << " " << setw(w1) << it+1;
	  outf << " " << setw(w2) << occ_vs_ion_lm_d[ioncnt][0];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][1];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][2];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][3];
	  outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][1];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][4];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][5];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][6];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][7];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][8];
	  outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][2];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][9];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][10];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][11];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][12];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][13];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][14];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_d[ioncnt][15];
	  outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][3];
	  outf << "  " << setw(w2) << occ_vs_ion_l_d[ioncnt][4];
	  outf << " " << key;
	  outf << endl;
	}
      }
      outf << "Up-Down spin"<< endl;
      outf << "   " << "ion   " << "type  ";
      for(int ilm=0;ilm<(int)LLMnames.size();ilm++) {
	outf << LLMnames[ilm] << "  ";
      }
      outf << endl;
      ioncnt=-1;
      for(int it=0;it<ntypes;it++) {
	for(int iat=0;iat<num_each_type.at(it);iat++) {
	  ioncnt++;
	  outf << "   " << setw(5) << ioncnt+1;
	  outf << " " << setw(5) << it+1;
	  outf << " " << setw(w2) << occ_vs_ion_lm_u[ioncnt][0]-occ_vs_ion_lm_d[ioncnt][0];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][1]-occ_vs_ion_lm_d[ioncnt][1];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][2]-occ_vs_ion_lm_d[ioncnt][2];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][3]-occ_vs_ion_lm_d[ioncnt][3];
	  outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][1]-occ_vs_ion_l_d[ioncnt][1];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][4]-occ_vs_ion_lm_d[ioncnt][4];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][5]-occ_vs_ion_lm_d[ioncnt][5];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][6]-occ_vs_ion_lm_d[ioncnt][6];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][7]-occ_vs_ion_lm_d[ioncnt][7];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][8]-occ_vs_ion_lm_d[ioncnt][8];
	  outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][2]-occ_vs_ion_l_d[ioncnt][2];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][9]-occ_vs_ion_lm_d[ioncnt][9];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][10]-occ_vs_ion_lm_d[ioncnt][10];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][11]-occ_vs_ion_lm_d[ioncnt][11];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][12]-occ_vs_ion_lm_d[ioncnt][12];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][13]-occ_vs_ion_lm_d[ioncnt][13];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][14]-occ_vs_ion_lm_d[ioncnt][14];
	  outf << "  " << setw(w2) << occ_vs_ion_lm_u[ioncnt][15]-occ_vs_ion_lm_d[ioncnt][15];
	  outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][3]-occ_vs_ion_l_d[ioncnt][3];
	  outf << "  " << setw(w2) << occ_vs_ion_l_u[ioncnt][4]-occ_vs_ion_l_d[ioncnt][4];
	  outf << " " << key;
	  outf << endl;
	}
      }
    }// sp
  }
}

// ***************************************************************************
// PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJF
// ***************************************************************************

// ***************************************************************************
// ProcessProjection
// ***************************************************************************
// The projections might be a complex amplitude, a true
// real probability, or a probability with a complex
// phase.  Depending on which, you need to process
// the projection differently.  Hence this routine, so
// all processing can be changed at once.
namespace pflow {
  std::complex<double> ProcessProjection(const std::complex<double>& proj) {
    std::complex<double> p;
    p=proj;
    // For proj=probability with phase.
    // This allows for only positive probabilities (does not exactly match vasp output).
    // return sqrt(p*conj(p));
    // For proj=real probability with complex part always zero.
    // This allows for negative probabilities (does exactly match vasp output).
    return p;
  }
}
// ***************************************************************************
// ReadInProj
// ***************************************************************************
// Need to read everything in twice if it is spin polarized.
// rspin=1 for spin polarized, rspin=2 for non-spin polarized.
namespace pflow {
  void ReadInProj(projdata& pd) {
    int have_aug;
    string s;
    char c;
    ifstream infile(pd.PROOUTinfile.c_str());
    aurostd::InFileExistCheck("ReadInProj",pd.PROOUTinfile.c_str(),infile,cerr);
    string sdum;
    // Initial data
    infile >> sdum; // Title
    infile >> sdum >> sdum >> sdum >> pd.nkpts;
    infile >> sdum >> sdum >> sdum >> pd.nbands;
    infile >> sdum >> sdum >> sdum >> pd.nions;
    infile >> pd.ntypes >> pd.ntypes;
    pd.num_each_type=vector<int> (pd.ntypes);
    for(int it=0;it<pd.ntypes;it++) {
      infile >> pd.num_each_type.at(it);
    }
    // For getting kpoint energies.
    int skip=2;
    if(pd.nions==1) skip=1;

    // Fermi weights (kpts (inner loop) and bands (outer loop)).
    pd.wfermi_u = matrix<double> (pd.nbands,pd.nkpts,0.0);
    for(int ib=0;ib<pd.nbands;ib++) {
      for(int ik=0;ik<pd.nkpts;ik++) {
	infile >> pd.wfermi_u[ib][ik];
	// tpx
	// cout<<"ib,ik,wfermi_u "<<ib<<" "<<ik<<" "<<pd.wfermi_u[ib][ik]<<endl;
      }
    }

    // Check to see if there are spd or spdf electrons (nl=3 or nl=4).
    if(aurostd::FindIfStringInStream("s      p      d      f",infile)) {
      pd.nl=4;
    }
    else{
      pd.nl=3;
    }
    if(pd.nl==3) {
      pd.nlm=9; // for s,p,d orbitals
      pd.nlmtot=12; // for s,p,d orbitals + p,d,all totals
    }

    // Projections
    matrix<std::complex<double> > m(pd.nbands,pd.nlm,0.0);
    matrix<matrix<std::complex<double> > > mm(pd.nkpts,pd.nions,m);
    pd.pdat_u=vector<matrix<matrix<std::complex<double> > > > (pd.ntypes,mm);
    vector<double> rval(pd.nlm);
    vector<double> ival(pd.nlm);
    for(int it=0;it<pd.ntypes;it++) {
      for(int ik=0;ik<pd.nkpts;ik++) {
	for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
	  for(int ib=0;ib<pd.nbands;ib++) {
	    for(int ilm=0;ilm<pd.nlm;ilm++) {
	      infile >> rval[ilm] >> ival[ilm];
	      std::complex<double> ctmp (rval[ilm],ival[ilm]);
	      pd.pdat_u[it][ik][iat][ib][ilm]=ctmp;
	      // tpx
	      //  cout << "PROJ ib,ilm " << ib << " " << ilm << " " <<  rval[ilm] << " " <<  ival[ilm] << " " << pd.pdat_u[it][ik][iat][ib][ilm] << endl;
	    }
	  }
	}
      }
    }

    // Augmented Projections
    getline(infile,s); // Gets final carriage return from previous line.
    getline(infile,s); // Gets "augmentation part"
    c = infile.peek();
    if(c=='#') {
      have_aug=0; // There are no augmented projections.
    }
    else{
      have_aug=1; // There are augmented projections.
    }
    if(have_aug) {
      for(int ik=0;ik<pd.nkpts;ik++) {
	for(int ib=0;ib<pd.nbands;ib++) {
	  for(int it=0;it<pd.ntypes;it++) {
	    for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
	      for(int ilm=0;ilm<pd.nlm;ilm++) {
		infile >> rval[ilm];
		//tpx
		// cout << "Rval " << rval[ilm] << " " << ilm  << endl;
	      }
	      for(int ilm=0;ilm<pd.nlm;ilm++) {
		pd.pdat_u[it][ik][iat][ib][ilm]=
		  pd.pdat_u[it][ik][iat][ib][ilm]+rval[ilm];
	      }
	    }
	  }
	}
      }
      getline(infile,s); // Gets final carriage return last line of augmentation charges.
    }//if have_aug

    // Kpt weights and values and energy vs. kpt,bnd
    pd.wkpt = vector<double> (pd.nkpts,0.0);
    pd.kpts = matrix<double> (pd.nkpts,3,0.0);
    pd.ener_k_b_u = matrix<double> (pd.nkpts,pd.nbands,0.0);
    getline(infile,s); // Gets line "# of k-points: ..."
    for(int ik=0;ik<pd.nkpts;ik++) {
      getline(infile,s); // Gets blank line.
      infile >>sdum>>sdum>>sdum>>pd.kpts[ik][0]>>pd.kpts[ik][1]>>pd.kpts[ik][2]>>sdum>>sdum>>pd.wkpt[ik];
      getline(infile,s); // Gets final carriage return from previous line.
      getline(infile,s); // Gets line.
      for(int ib=0;ib<pd.nbands;ib++) {
	infile >>sdum>>sdum>>sdum>>sdum>>pd.ener_k_b_u[ik][ib]>>sdum>>sdum>>sdum;
	for(int ii=0;ii<pd.nions+skip;ii++) {
	  // Gets lines (nbands*(nions+skip)).  Note that
	  // the numer of blank lines varies in different vasp versions
	  // so this loop does not count blank lines.
	  getline(infile,s);
	  if(aurostd::RemoveSpaces(s).size()==0) {
	    ii=ii-1; // Don't count blank lines.
	  }
	}
      }
    }

    // Do again if spin polarized.  Look for another PROOUT title line.
    // At this point lines could be blank, end of file, or PROOUT.  We want
    // to keep reading if they are blank, stop if they are EOF, and do
    // spin-polarized if they are PROOUT.
    infile >> s; // Skip all whitespace and get PROOUT title if it exists.
    if(aurostd::RemoveSpaces(s)=="PROOUT") {
      // Set spin polarization dependent variables.
      pd.sp=1;
      pd.rspin=1;
    
      getline(infile,s); // Gets endline after PROOUT.
      getline(infile,s); // get line.
      getline(infile,s); // get line.
    
      // Fermi weights (kpts (inner loop) and bands (outer loop)).
      pd.wfermi_d = matrix<double> (pd.nbands,pd.nkpts,0.0);
      for(int ib=0;ib<pd.nbands;ib++) {
	for(int ik=0;ik<pd.nkpts;ik++) {
	  infile >> pd.wfermi_d[ib][ik];
	}
      }

      // Projections
      matrix<std::complex<double> > m(pd.nbands,pd.nlm,0.0);
      matrix<matrix<std::complex<double> > > mm(pd.nkpts,pd.nions,m);
      pd.pdat_d=vector<matrix<matrix<std::complex<double> > > > (pd.ntypes,mm);
      vector<double> rval(pd.nlm);
      vector<double> ival(pd.nlm);
      for(int it=0;it<pd.ntypes;it++) {
	for(int ik=0;ik<pd.nkpts;ik++) {
	  for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
	    for(int ib=0;ib<pd.nbands;ib++) {
	      for(int ilm=0;ilm<pd.nlm;ilm++) {
		infile >> rval[ilm] >> ival[ilm];
	      }
	      for(int ilm=0;ilm<pd.nlm;ilm++) {
		std::complex<double> ctmp (rval[ilm],ival[ilm]);
		pd.pdat_d[it][ik][iat][ib][ilm]=ctmp;
	      }
	    }
	  }
	}
      }

      // Augmented Projections
      string s;
      getline(infile,s); // Gets final carriage return from previous line.
      getline(infile,s); // Gets "augmentation part"
      c = infile.peek();
      if(c=='#') {
	have_aug=0; // There are no augmented projections.
      }
      else{
	have_aug=1; // There are augmented projections.
      }
      if(have_aug) {
	for(int ik=0;ik<pd.nkpts;ik++) {
	  for(int ib=0;ib<pd.nbands;ib++) {
	    for(int it=0;it<pd.ntypes;it++) {
	      for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
		for(int ilm=0;ilm<pd.nlm;ilm++) {
		  infile >> rval[ilm];
		}
		for(int ilm=0;ilm<pd.nlm;ilm++) {
		  pd.pdat_d[it][ik][iat][ib][ilm]=
		    pd.pdat_d[it][ik][iat][ib][ilm]+rval[ilm];
		}
	      }
	    }
	  }
	}
	getline(infile,s); // Gets final carriage return last line of augmentation charges.
      }//if have_aug

      // Kpt weights and values and energy vs. kpt,bnd
      pd.ener_k_b_d = matrix<double> (pd.nkpts,pd.nbands,0.0);
      getline(infile,s); // Gets line "# of k-points: ..."
      for(int ik=0;ik<pd.nkpts;ik++) {
	getline(infile,s); // Gets blank line.
	infile >>sdum>>sdum>>sdum>>pd.kpts[ik][0]>>pd.kpts[ik][1]>>pd.kpts[ik][2]>>sdum>>sdum>>pd.wkpt[ik];
	getline(infile,s); // Gets final carriage return from previous line.
	getline(infile,s); // Gets line.
	for(int ib=0;ib<pd.nbands;ib++) {
	  infile >>sdum>>sdum>>sdum>>sdum>>pd.ener_k_b_d[ik][ib]>>sdum>>sdum>>sdum;
	  for(int ii=0;ii<pd.nions+skip;ii++) {
	    // Gets lines (nbands*(nions+skip)).  Note that
	    // the numer of blank lines varies in different vasp versions
	    // so this loop does not count blank lines.
	    getline(infile,s);
	    if(aurostd::RemoveSpaces(s).size()==0) {
	      ii=ii-1; // Don't count blank lines.
	    }
	  }
	}
      }

    }// if spin polarized

  }// end
}

// ***************************************************************************
// CalcNeatProj
// ***************************************************************************
// Calculates all the basic occupancies.  If only_occ=true
// then it does this only for occupied states.  If only_occ=false
// then it does this for all states, which does not really
// give occupations, but sets things up for PDOS calculations.
namespace pflow {
  void CalcNeatProj(projdata& pd, int only_occ) {
    int ioncnt;

    // Get occ_vs_ion_kpt_bnd_lm_(u,d)
    // Get occ_vs_ion_kpt_lm_(u,d)
    // Get occ_vs_ion_bnd_lm_(u,d)
    // Get occ_vs_ion_lm_(u,d)
    matrix<double> m1 (pd.nbands,pd.nlm_max,0.0);
    matrix<double> m2 (pd.nkpts,pd.nlm_max,0.0);
    matrix<double> m3 (pd.nbands,pd.nlm_max,0.0);
    pd.occ_vs_ion_kpt_bnd_lm_u = matrix<matrix<double> > (pd.nions,pd.nkpts,m1);
    pd.occ_vs_ion_kpt_lm_u = vector<matrix<double> > (pd.nions,m2);
    pd.occ_vs_ion_bnd_lm_u = vector<matrix<double> > (pd.nions,m3);
    pd.occ_vs_ion_lm_u = matrix<double> (pd.nions,pd.nlm_max,0.0);  
    if(pd.sp) {
      pd.occ_vs_ion_kpt_bnd_lm_d = matrix<matrix<double> > (pd.nions,pd.nkpts,m1);
      pd.occ_vs_ion_kpt_lm_d = vector<matrix<double> > (pd.nions,m2);
      pd.occ_vs_ion_bnd_lm_d = vector<matrix<double> > (pd.nions,m3);
      pd.occ_vs_ion_lm_d = matrix<double> (pd.nions,pd.nlm_max,0.0);  
    }
    ioncnt=-1;
    for(int it=0;it<pd.ntypes;it++) {
      for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
	ioncnt++;
	for(int ik=0;ik<pd.nkpts;ik++) {
	  for(int ib=0;ib<pd.nbands;ib++) {
	    for(int ilm=0;ilm<pd.nlm;ilm++) {
	      double temp;
	      temp=ProcessProjection(pd.pdat_u[it][ik][iat][ib][ilm]).real();
	      if(only_occ) {
		pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm]=pd.rspin*pd.wfermi_u[ib][ik]*(temp);
	      }
	      else{
		pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm]=pd.rspin*(temp);
	      }
	      pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][ilm]=
		pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][ilm]+
		pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm];
	      pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][ilm]=
		pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][ilm]+
		pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm];
	      pd.occ_vs_ion_lm_u[ioncnt][ilm]=
		pd.occ_vs_ion_lm_u[ioncnt][ilm]+
		pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][ilm];
	      if(pd.sp) {
		temp=ProcessProjection(pd.pdat_d[it][ik][iat][ib][ilm]).real();
		if(only_occ) {
		  pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm]=pd.rspin*pd.wfermi_d[ib][ik]*(temp);
		}
		else{
		  pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm]=pd.rspin*(temp);
		}
		pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][ilm]=
		  pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][ilm]+
		  pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm];
		pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][ilm]=
		  pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][ilm]+
		  pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm];
		pd.occ_vs_ion_lm_d[ioncnt][ilm]=
		  pd.occ_vs_ion_lm_d[ioncnt][ilm]+
		  pd.wkpt[ik]*pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][ilm];
	      }
	    }
	  }
	}
      }
    }

    // Get sums over m.
    // Get occ_vs_ion_l_(u,d)
    // Get occ_vs_ion_kpt_l_(u,d)
    // Get occ_vs_ion_bnd_l_(u,d)
    // Get occ_vs_ion_kpt_bnd_l_(u,d)
    // This assumes the usual spd orbitals.
    matrix<double> n1 (pd.nbands,pd.nl_max+1,0.0);
    matrix<double> n2 (pd.nkpts,pd.nl_max+1,0.0);
    matrix<double> n3 (pd.nbands,pd.nl_max+1,0.0);
    pd.occ_vs_ion_kpt_bnd_l_u = matrix<matrix<double> > (pd.nions,pd.nkpts,n1);
    pd.occ_vs_ion_kpt_l_u = vector<matrix<double> > (pd.nions,n2);
    pd.occ_vs_ion_bnd_l_u = vector<matrix<double> > (pd.nions,n3);
    pd.occ_vs_ion_l_u = matrix<double> (pd.nions,pd.nl_max+1,0.0);  
    if(pd.sp) {
      pd.occ_vs_ion_kpt_bnd_l_d = matrix<matrix<double> > (pd.nions,pd.nkpts,n1);
      pd.occ_vs_ion_kpt_l_d = vector<matrix<double> > (pd.nions,n2);
      pd.occ_vs_ion_bnd_l_d = vector<matrix<double> > (pd.nions,n3);
      pd.occ_vs_ion_l_d = matrix<double> (pd.nions,pd.nl_max+1,0.0);  
    }

    ioncnt=-1;
    for(int it=0;it<pd.ntypes;it++) {
      for(int iat=0;iat<pd.num_each_type.at(it);iat++) {
	ioncnt++;
	for(int ik=0;ik<pd.nkpts;ik++) {
	  for(int ib=0;ib<pd.nbands;ib++) {

	    // S
	    pd.occ_vs_ion_l_u[ioncnt][0]+=pd.occ_vs_ion_lm_u[ioncnt][0];
	    pd.occ_vs_ion_kpt_l_u[ioncnt][ik][0]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][0];
	    pd.occ_vs_ion_bnd_l_u[ioncnt][ib][0]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][0];
	    pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][0]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][0];
	    if(pd.sp) {
	      pd.occ_vs_ion_l_d[ioncnt][0]+=pd.occ_vs_ion_lm_d[ioncnt][0];
	      pd.occ_vs_ion_kpt_l_d[ioncnt][ik][0]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][0];
	      pd.occ_vs_ion_bnd_l_d[ioncnt][ib][0]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][0];
	      pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][0]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][0];
	    }

	    // P
	    for(int ip=0;ip<3;ip++) {
	      pd.occ_vs_ion_l_u[ioncnt][1]+=pd.occ_vs_ion_lm_u[ioncnt][1+ip];
	      pd.occ_vs_ion_kpt_l_u[ioncnt][ik][1]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][1+ip];
	      pd.occ_vs_ion_bnd_l_u[ioncnt][ib][1]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][1+ip];
	      pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][1]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][1+ip];
	      if(pd.sp) {
		pd.occ_vs_ion_l_d[ioncnt][1]+=pd.occ_vs_ion_lm_d[ioncnt][1+ip];
		pd.occ_vs_ion_kpt_l_d[ioncnt][ik][1]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][1+ip];
		pd.occ_vs_ion_bnd_l_d[ioncnt][ib][1]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][1+ip];
		pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][1]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][1+ip];
	      }
	    }

	    // D
	    for(int id=0;id<5;id++) {
	      pd.occ_vs_ion_l_u[ioncnt][2]+=pd.occ_vs_ion_lm_u[ioncnt][4+id];
	      pd.occ_vs_ion_kpt_l_u[ioncnt][ik][2]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][4+id];
	      pd.occ_vs_ion_bnd_l_u[ioncnt][ib][2]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][4+id];
	      pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][2]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][4+id];
	      if(pd.sp) {
		pd.occ_vs_ion_l_d[ioncnt][2]+=pd.occ_vs_ion_lm_d[ioncnt][4+id];
		pd.occ_vs_ion_kpt_l_d[ioncnt][ik][2]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][4+id];
		pd.occ_vs_ion_bnd_l_d[ioncnt][ib][2]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][4+id];
		pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][2]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][4+id];
	      }
	    }

	    // F
	    for(int id=0;id<7;id++) {
	      pd.occ_vs_ion_l_u[ioncnt][3]+=pd.occ_vs_ion_lm_u[ioncnt][9+id];
	      pd.occ_vs_ion_kpt_l_u[ioncnt][ik][3]+=pd.occ_vs_ion_kpt_lm_u[ioncnt][ik][9+id];
	      pd.occ_vs_ion_bnd_l_u[ioncnt][ib][3]+=pd.occ_vs_ion_bnd_lm_u[ioncnt][ib][9+id];
	      pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][3]+=pd.occ_vs_ion_kpt_bnd_lm_u[ioncnt][ik][ib][9+id];
	      if(pd.sp) {
		pd.occ_vs_ion_l_d[ioncnt][3]+=pd.occ_vs_ion_lm_d[ioncnt][9+id];
		pd.occ_vs_ion_kpt_l_d[ioncnt][ik][3]+=pd.occ_vs_ion_kpt_lm_d[ioncnt][ik][9+id];
		pd.occ_vs_ion_bnd_l_d[ioncnt][ib][3]+=pd.occ_vs_ion_bnd_lm_d[ioncnt][ib][9+id];
		pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][3]+=pd.occ_vs_ion_kpt_bnd_lm_d[ioncnt][ik][ib][9+id];
	      }
	    }

	  }//ib
	}//ik

	// Normalize sums over m for multiple counting.
	for(int il=0;il<pd.nl;il++) {
	  pd.occ_vs_ion_l_u[ioncnt][il]/=(pd.nbands*pd.nkpts);
	  if(pd.sp)  pd.occ_vs_ion_l_d[ioncnt][il]/=(pd.nbands*pd.nkpts);
	  for(int ik=0;ik<pd.nkpts;ik++) {
	    pd.occ_vs_ion_kpt_l_u[ioncnt][ik][il]/=(pd.nbands);

	    if(pd.sp) pd.occ_vs_ion_kpt_l_d[ioncnt][ik][il]/=(pd.nbands);
	  }
	  for(int ib=0;ib<pd.nbands;ib++) {
	    pd.occ_vs_ion_bnd_l_u[ioncnt][ib][il]/=(pd.nkpts);
	    if(pd.sp) pd.occ_vs_ion_bnd_l_d[ioncnt][ib][il]/=(pd.nkpts);
	  }
	}

	// Calculate total sums over L.

	int iid=pd.nl_max;
	for(int id=0;id<pd.nl_max;id++) {
	  pd.occ_vs_ion_l_u[ioncnt][iid]+=pd.occ_vs_ion_l_u[ioncnt][id];	
	  if(pd.sp) pd.occ_vs_ion_l_d[ioncnt][iid]+=pd.occ_vs_ion_l_d[ioncnt][id];
	  for(int ik=0;ik<pd.nkpts;ik++) {
	    pd.occ_vs_ion_kpt_l_u[ioncnt][ik][iid]+=pd.occ_vs_ion_kpt_l_u[ioncnt][ik][id];
	    if(pd.sp) pd.occ_vs_ion_kpt_l_d[ioncnt][ik][iid]+=pd.occ_vs_ion_kpt_l_d[ioncnt][ik][id];
	  }
	  for(int ib=0;ib<pd.nbands;ib++) {
	    pd.occ_vs_ion_bnd_l_u[ioncnt][ib][iid]+=pd.occ_vs_ion_bnd_l_u[ioncnt][ib][id];
	    if(pd.sp) pd.occ_vs_ion_bnd_l_d[ioncnt][ib][iid]+=pd.occ_vs_ion_bnd_l_d[ioncnt][ib][id];
	  }
	  for(int ik=0;ik<pd.nkpts;ik++) {
	    for(int ib=0;ib<pd.nbands;ib++) {
	      pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][iid]+=pd.occ_vs_ion_kpt_bnd_l_u[ioncnt][ik][ib][id];
	      if(pd.sp) pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][iid]+=pd.occ_vs_ion_kpt_bnd_l_d[ioncnt][ik][ib][id];
	    }
	  }
	}

      }// at
    }// type

    // Create single occ variable for all the l,m and sums over m together (for lmtot).
    // This way you can just loop over lm from 0 to lmtot-1.
    m1 = matrix<double> (pd.nbands,pd.nlmtot,0.0);
    pd.occ_vs_ion_kpt_bnd_lmtot_u = matrix<matrix<double> > (pd.nions,pd.nkpts,m1);
    if(pd.sp) pd.occ_vs_ion_kpt_bnd_lmtot_d = matrix<matrix<double> > (pd.nions,pd.nkpts,m1);
    for(int ii=0;ii<pd.nions;ii++) {
      for(int ik=0;ik<pd.nkpts;ik++) {
	for(int ib=0;ib<pd.nbands;ib++) {
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][0]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][0];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][1]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][1];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][2]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][2];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][3]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][3];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][4]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][1];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][5]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][4];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][6]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][5];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][7]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][6];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][8]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][7];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][9]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][8];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][10]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][2];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][11]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][9];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][12]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][10];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][13]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][11];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][14]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][12];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][15]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][13];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][16]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][14];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][17]=pd.occ_vs_ion_kpt_bnd_lm_u[ii][ik][ib][15];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][18]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][3];
	  pd.occ_vs_ion_kpt_bnd_lmtot_u[ii][ik][ib][19]=pd.occ_vs_ion_kpt_bnd_l_u[ii][ik][ib][4];
	  if(pd.sp) {
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][0]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][0];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][1]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][1];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][2]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][2];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][3]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][3];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][4]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][1];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][5]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][4];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][6]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][5];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][7]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][6];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][8]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][7];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][9]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][8];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][10]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][2];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][11]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][9];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][12]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][10];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][13]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][11];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][14]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][12];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][15]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][13];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][16]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][14];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][17]=pd.occ_vs_ion_kpt_bnd_lm_d[ii][ik][ib][15];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][18]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][3];
	    pd.occ_vs_ion_kpt_bnd_lmtot_d[ii][ik][ib][19]=pd.occ_vs_ion_kpt_bnd_l_d[ii][ik][ib][4];
	  }
	}
      }
    }
  }
}

// ***************************************************************************
// ReadInPDOSData
// ***************************************************************************
// Reads in parameters that determine what partial DOS
//   will be calculated.
namespace pflow {
  void ReadInPDOSData(const pflow::projdata& prd, pflow::pdosdata& pdd) {

    ifstream infile(pdd.PDOSinfile.c_str());
    aurostd::InFileExistCheck("ReadInPDOSData",pdd.PDOSinfile.c_str(),infile,cerr);

    // Defaults (emin,emax,nbins,smooth_sigma,print_params)
    double emn=prd.ener_k_b_u[0][0];
    double emx=emn;
    double e;
    for(int ik=0;ik<prd.nkpts;ik++) {
      for(int ib=0;ib<prd.nbands;ib++) {
	e = prd.ener_k_b_u[ik][ib];
	if(emn>e) emn=e;
	if(emx<e) emx=e;
	if(prd.sp) e = prd.ener_k_b_d[ik][ib];
	if(emn>e) emn=e;
	if(emx<e) emx=e;
      }
    }
    pdd.emin=emn-0.5;
    pdd.emax=emx+0.5;
    pdd.nbins=300;
    pdd.smooth_sigma=(pdd.emax-pdd.emin)/pdd.nbins;
    pdd.print_params=0;
	  
    // Read in all the tokens.

    string s,s_ns;
    vector<string> token_vec;
    vector<string> val_vec;
    while (!infile.eof()) {
      getline(infile,s);
      s_ns=aurostd::RemoveSpaces(s); // Get string with no spaces.
      if(s_ns.size()>0) { // Make sure line was not blank
	if(s_ns[0]!='#') { // Exclude comment lines
	  string token;
	  string sval;
	  int id=s.find('=');
	  if(id>=(int)s.length()) {
	    cout << "ERROR: The following token is incorrectly formatted: " << s << endl;
	    cout << "ERROR: Exiting" << endl;
	    exit(1);
	  }
	  token=s.substr(0,id);
	  token=aurostd::RemoveSpaces(token);
	  int i_f=std::min((int)s.length(),(int)s.find('#')); // End of string
	  sval=s.substr(id+1,i_f-id-1);
	  token_vec.push_back(token);
	  val_vec.push_back(sval);
	} //if s_ns[0]!=#
      } // if s_ns.size>0

    } // while !infile.eof

    // Read in values for variables in the pdosdata object.
    vector<int> tmp;
    int atom_cnt=0;
    for(int i=0;i<(int)token_vec.size();i++) {
      int id=0;
      string tok=token_vec[i];
      string sval;
      // double dval;
      int ival;
      int found_token=0;
      if(tok=="ATOMS") {
	id=0;
	atom_cnt++;
	tmp.clear();
	// set up atoms matrix
	while (id<(int)val_vec[i].size()) {
	  string s;
	  s=aurostd::GetNextVal(val_vec[i],id);
	  if(s.size()>0) {
	    ival=atoi(s.c_str());
	    tmp.push_back(ival);
	  }
	}
	pdd.pdos_at.push_back(tmp);
	// add vector to pdos_k,pdos_b,pdos_lm
	tmp.clear();
	pdd.pdos_k.push_back(tmp);
	pdd.pdos_b.push_back(tmp);
	pdd.pdos_lm.push_back(tmp);
	found_token=1;
      }// ATOMS
      if(tok=="KPOINTS") {
	if((int)pdd.pdos_k.size()!=atom_cnt || atom_cnt==0) AtomCntError(tok,(int)pdd.pdos_k.size(),atom_cnt);
	id=0;
	tmp.clear();
	// set up atoms matrix
	while (id<(int)val_vec[i].size()) {
	  string s;
	  s=aurostd::GetNextVal(val_vec[i],id);
	  if(s.size()>0) {
	    ival=atoi(s.c_str());
	    pdd.pdos_k[pdd.pdos_k.size()-1].push_back(ival);
	  }
	}
	found_token=1;
      }// KPOINTS
      if(tok=="BANDS") {
	id=0;
	tmp.clear();
	// set up atoms matrix
	while (id<(int)val_vec[i].size()) {
	  string s;
	  s=aurostd::GetNextVal(val_vec[i],id);
	  if(s.size()>0) {
	    ival=atoi(s.c_str());
	    pdd.pdos_b[pdd.pdos_b.size()-1].push_back(ival);
	  }
	}
	found_token=1;
      }// LMVALUES
      if(tok=="LMVALUES") {
	if((int)pdd.pdos_lm.size()!=atom_cnt || atom_cnt==0) AtomCntError(tok,(int)pdd.pdos_lm.size(),atom_cnt);
	id=0;
	tmp.clear();
	// set up atoms matrix
	while (id<(int)val_vec[i].size()) {
	  string s;
	  s=aurostd::GetNextVal(val_vec[i],id);
	  if(s.size()>0) {
	    ival=atoi(s.c_str());
	    pdd.pdos_lm[pdd.pdos_lm.size()-1].push_back(ival);
	  }
	}
	found_token=1;
      }// LMVALUES
      if(tok=="EMIN") {
	id=0;
	pdd.emin=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// EMIN
      if(tok=="EMAX") {
	id=0;
	pdd.emax=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// EMAX
      if(tok=="NBINS") {
	id=0;
	pdd.nbins=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// NBIN
      if(tok=="SMOOTH_SIGMA") {
	id=0;
	pdd.smooth_sigma=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// SMOOTH_SIGMA
      if(tok=="PRINT_PARAMS") {
	id=0;
	pdd.print_params=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// PRINT_PARAMS

      if(!found_token) {
	cerr << "ERROR: You have input a token " << tok <<endl;
	cerr << "ERROR: This token is not recognized. Exiting! " << endl;
	exit(1);
      }

    } // for i

    // Defaults (all atoms, all kpts, all bands, lm = nlmtot = all l+m)
    vector<int> vv(4);
    vv[0]=prd.nions;
    vv[1]=prd.nkpts;
    vv[2]=prd.nbands;
    vv[3]=prd.nlmtot;
    int mx = *max_element(vv.begin(),vv.end());
    vector<int> cnt(mx);
    for(int i=1;i<=mx;i++) {
      cnt[i-1]=i;
    }
    //  int mx = max(prd.nions,prd.nkpts,prd.nbands,prd.nlmtot);
    vector<int> atv(prd.nions);
    for(int i=0;i<(int)pdd.pdos_at.size();i++) {
      if(pdd.pdos_at[i].size()==0) pdd.pdos_at[i] = vector<int> (cnt.begin(),cnt.begin()+prd.nions);
      if(pdd.pdos_k[i].size()==0) pdd.pdos_k[i] = vector<int> (cnt.begin(),cnt.begin()+prd.nkpts);
      if(pdd.pdos_b[i].size()==0) pdd.pdos_b[i] = vector<int> (cnt.begin(),cnt.begin()+prd.nbands);
      if(pdd.pdos_lm[i].size()==0) pdd.pdos_lm[i] = vector<int> (1,prd.nlmtot);
    }

    // Check that all data is within bounds
    for(int i=0;i<(int)pdd.pdos_at.size();i++) { // cases
      for(int ia=0;ia<(int)pdd.pdos_at[i].size();ia++) {
	int atp=pdd.pdos_at[i][ia];
	if(atp<1 || atp>prd.nions) {
	  cerr << "ERROR: ProjFuncs/ReadPDOSdat" << endl;
	  cerr << "ERROR: Error in case: " << i+1 << endl;
	  cerr << "ERROR: You have entered too low or high an atom number in entry: " << ia+1 << endl;
	  cerr << "ERROR: Bad value is: " << atp << endl;
	  cerr << "ERROR: Exiting." << endl;
	  exit(1);
	}
      }
      for(int ik=0;ik<(int)pdd.pdos_k[i].size();ik++) {
	int kp=pdd.pdos_k[i][ik];
	if(kp<1 || kp>prd.nkpts) {
	  cerr << "ERROR: ProjFuncs/ReadPDOSdat" << endl;
	  cerr << "ERROR: Error in case: " << i+1 << endl;
	  cerr << "ERROR: You have entered too low or high a kpt number in entry: " << ik+1 << endl;
	  cerr << "ERROR: Bad value is: " << kp << endl;
	  cerr << "ERROR: Exiting." << endl;
	  exit(1);
	}
      }
      for(int ib=0;ib<(int)pdd.pdos_b[i].size();ib++) {
	int bp=pdd.pdos_b[i][ib];
	if(bp<1 || bp>prd.nbands) {
	  cerr << "ERROR: ProjFuncs/ReadPDOSdat" << endl;
	  cerr << "ERROR: Error in case: " << i+1 << endl;
	  cerr << "ERROR: You have entered too low or high a bnd number in entry: " << ib+1 << endl;
	  cerr << "ERROR: Bad value is: " << bp << endl;
	  cerr << "ERROR: Exiting." << endl;
	  exit(1);
	}
      }
      for(int ilm=0;ilm<(int)pdd.pdos_lm[i].size();ilm++) {
	int lmp=pdd.pdos_lm[i][ilm];
	if(lmp<1 || lmp>prd.nlmtot) {
	  cerr << "ERROR: ProjFuncs/ReadPDOSdat" << endl;
	  cerr << "ERROR: Error in case: " << i+1 << endl;
	  cerr << "ERROR: You have entered too low or high an lm number in entry: " << ilm+1 << endl;
	  cerr << "ERROR: Bad value is: " << lmp << endl;
	  cerr << "ERROR: Exiting." << endl;
	  exit(1);
	}
      }
    }
  }// end routine
}

// ***************************************************************************
// CalcPDOS
// ***************************************************************************
// Calculates the total PDOS based on the settings in pdd.
namespace pflow {
  void CalcPDOS(const projdata& prd, pdosdata& pdd) {
    double de=(pdd.emax-pdd.emin)/pdd.nbins;
    double e;
    int ibin;
    int sp_size=1;
    if(prd.sp) sp_size=3;
    pdd.pdos = matrix<double> (pdd.nbins,2*sp_size+1,0.0);
    for(int ib=0;ib<pdd.nbins;ib++) {
      pdd.pdos[ib][0]=pdd.emin+(double)ib*de;
    }
    for(int i=0;i<(int)pdd.pdos_at.size();i++) { // cases
      for(int ia=0;ia<(int)pdd.pdos_at[i].size();ia++) {
	int at=pdd.pdos_at[i][ia]-1;
	for(int ik=0;ik<(int)pdd.pdos_k[i].size();ik++) {
	  int kp=pdd.pdos_k[i][ik]-1;
	  for(int ib=0;ib<(int)pdd.pdos_b[i].size();ib++) {
	    int bp=pdd.pdos_b[i][ib]-1;
	    for(int ilm=0;ilm<(int)pdd.pdos_lm[i].size();ilm++) {
	      int lmp=pdd.pdos_lm[i][ilm]-1;
	      e = prd.ener_k_b_u[kp][bp];
	      ibin=int((e-pdd.emin)/de);
	      if(ibin>=0 && ibin<pdd.nbins) {
		pdd.pdos[ibin][1]+=prd.occ_vs_ion_kpt_bnd_lmtot_u[at][kp][bp][lmp]*prd.wkpt[kp];
		if(prd.sp) {
		  e = prd.ener_k_b_d[kp][bp];
		  ibin=int((e-pdd.emin)/de);
		  pdd.pdos[ibin][2]+=prd.occ_vs_ion_kpt_bnd_lmtot_d[at][kp][bp][lmp]*prd.wkpt[kp];
		}
	      }// ibin in range
	    }//ilm
	  }//ik
	}//ib
      }//iat
    }//cases
    // Now do Up-Dn if spin polarized.
    if(prd.sp) {
      for(int ibin=0;ibin<pdd.nbins;ibin++) {
	pdd.pdos[ibin][3]=pdd.pdos[ibin][1]-pdd.pdos[ibin][2];
      }
    }

    // Normalize DOS to a density
    for(int ib=0;ib<pdd.nbins;ib++) {
      for(int i=1;i<sp_size+1;i++) {
	pdd.pdos[ib][i]/=de;
      }
    }

    // Smooth DOS
    SmoothPDOS(prd,pdd);

    // Get integrated DOS
    for(int i=1;i<sp_size+1;i++) {
      pdd.pdos[0][i+sp_size]=pdd.pdos[0][i]*de;
    }
    for(int ib=1;ib<pdd.nbins;ib++) {
      for(int i=1;i<sp_size+1;i++) {
	pdd.pdos[ib][i+sp_size]=pdd.pdos[ib-1][i+sp_size]+pdd.pdos[ib][i]*de;
      }
    }
  }//end routine
}

// ***************************************************************************
// SmoothPDOS
// ***************************************************************************
// Smooths the total PDOS based on Gaussian smoothing.
namespace pflow {
  void SmoothPDOS(const projdata& prd, pdosdata& pdd) {
    double de=(pdd.emax-pdd.emin)/pdd.nbins;
    //   double e; // DANE not used
    //  int ibin; // DANE not used
    int sp_size=1;
    if(prd.sp) sp_size=3;

    // Get smoothed pdos.
    for(int is=1;is<=sp_size;is++) {
      vector<double> tdos(pdd.nbins);
      for(int ib=0;ib<pdd.nbins;ib++) {
	tdos[ib]=pdd.pdos[ib][is];
      }
      double sig=pdd.smooth_sigma/de;
      tdos=SmoothFunc(tdos,sig);
      for(int ib=0;ib<pdd.nbins;ib++) {
	pdd.pdos[ib][is]=tdos[ib];
      }
    }

    /* Old smooth DOS calc - does not seem to work right!
       int range=(int)(5.0*pdd.smooth_sigma/de); // Uses gaussian out to 5 sigma.
       matrix<double> wt(pdd.nbins,2*range+1,0.0);

       // Get Gaussian weights
       vector<double> norm(pdd.nbins,0.0);
       for(int ib=0;ib<pdd.nbins;ib++) {
       for(int i=-range;i<=range;i++) {
       double x=(double)(i)*de+de/2.0;
       if((ib+i)>=0 && (ib+i)<pdd.nbins) {
       wt[ib][i+range]=Normal(x,de/2,pdd.smooth_sigma);
       norm[ib]=norm[ib]+wt[ib][i+range];
       }
       }
       }
       // Normalize to one
       for(int ib=0;ib<pdd.nbins;ib++) {
       for(int i=-range;i<=range;i++) {
       wt[ib][i+range]=wt[ib][i+range]/norm[ib];
       }
       }

       // Average in weighted nearby bins.
       for(int is=1;is<=sp_size;is++) {
       for(int ib=0;ib<pdd.nbins;ib++) {
       double tdos=0;
       for(int i=-range;i<=range;i++) {
       if((ib+i)>0 && (ib+i)<pdd.nbins) {
       tdos=tdos+wt[ib][i+range]*pdd.pdos[ib+i][is];
       }
       }
       pdd.pdos[ib][is]=tdos;
       }
       }
    */
  }
}

// ***************************************************************************
// AtomCntError
// ***************************************************************************
// Smooths the total PDOS based on Gaussian smoothing.
namespace pflow {
  void AtomCntError(const string& tok, const int tokcnt, const int atom_cnt) {
    cerr << "ERROR: in AtomCntError " <<endl;
    cerr << "ERROR: The token " << tok << " has been defined (possibly by default) this many times: " << tokcnt << endl;
    cerr << "ERROR: You have used the token ATOMS this many times: " << atom_cnt << endl;
    cerr << "ERROR: These must be equal and >0 - the token ATOMS must precede each case and must occurr at least once. "<< endl;
    cerr << "ERROR: This error is probably due to your having forgot to use the ATOMS token. Please type ATOMS=*** as the first line for each case."<< endl;
    cerr << "ERROR: Exiting!" << endl;
    exit(1);
  }
}

// ***************************************************************************
// SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUN
// ***************************************************************************

// ***************************************************************************
// ReadSumDOSParams
// ***************************************************************************
// Reads in parameters from PDOSParams_infile
// that decide what summations to perform over pdos.
// This routine makes use of the pdosdata object which
// contains some parameters that are not used (parameters
// for kpoints, bands, etc.)

namespace pflow {
  void ReadSumDOSParams(ifstream& infile, pflow::pdosdata& pdd) {
    // Defaults (emin,emax,nbins,smooth_sigma,print_params,efermi)
    pdd.emin=0.0;
    pdd.emax=0.0;
    pdd.nbins=301;
    pdd.smooth_sigma=(pdd.emax-pdd.emin)/pdd.nbins;
    pdd.print_params=0;
    pdd.efermi=-999;
    pdd.spin=1; // 1: non spin polarized, 2 spin polarized.          
    pdd.nlm=9; // all the possible l and m:  9 for s,p,d,  16 for s,p,d,f

    // Read in all the tokens.
    string s,s_ns;
    vector<string> token_vec;
    vector<string> val_vec;
    while (!infile.eof()) {
      getline(infile,s);
      s_ns=aurostd::RemoveSpaces(s); // Get string with no spaces.
      if(s_ns.size()>0) { // Make sure line was not blank
	if(s_ns[0]!='#') { // Exclude comment lines
	  string token;
	  string sval;
	  int id=s.find('=');
	  if(id>=(int)s.length()) {
	    cout << "ERROR: The following token is incorrectly formatted: " << s << endl;
	    cout << "ERROR: Exiting" << endl;
	    exit(1);
	  }
	  token=s.substr(0,id);
	  token=aurostd::RemoveSpaces(token);
	  int i_f=std::min((int)s.length(),(int)s.find('#')); // End of string
	  sval=s.substr(id+1,i_f-id-1);
	  token_vec.push_back(token);
	  val_vec.push_back(sval);
	} //if s_ns[0]!=#
      } // if s_ns.size>0
    } // while !infile.eof

    // Read in values for variables in the pdosdata object.
    vector<int> tmp;
    int atom_cnt=0;
    for(int i=0;i<(int)token_vec.size();i++) {
      int id=0;
      string tok=token_vec[i];
      string sval;
      // double dval; // DANE not used
      int ival;
      int found_token=0;
      if(tok=="ATOMS") {
	id=0;
	atom_cnt++;
	tmp.clear();
	// set up atoms matrix
	while (id<(int)val_vec[i].size()) {
	  string s;
	  s=aurostd::GetNextVal(val_vec[i],id);
	  if(s.size()>0) {
	    ival=atoi(s.c_str());
	    tmp.push_back(ival);
	  }
	}
	pdd.pdos_at.push_back(tmp);
	// add vector to pdos_lm
	tmp.clear();
	pdd.pdos_lm.push_back(tmp);
	found_token=1;
      }// ATOMS
      if(tok=="LMVALUES") {
	if((int)pdd.pdos_lm.size()!=atom_cnt || atom_cnt==0) AtomCntError(tok,(int)pdd.pdos_lm.size(),atom_cnt);
	id=0;
	tmp.clear();
	// set up atoms matrix
	while (id<(int)val_vec[i].size()) {
	  string s;
	  s=aurostd::GetNextVal(val_vec[i],id);
	  if(s.size()>0) {
	    ival=atoi(s.c_str());
	    pdd.pdos_lm[pdd.pdos_lm.size()-1].push_back(ival);
	  }
	}
	found_token=1;
      }// LMVALUES
      if(tok=="SMOOTH_SIGMA") {
	id=0;
	pdd.smooth_sigma=atof(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// SMOOTH_SIGMA
      if(tok=="PRINT_PARAMS") {
	id=0;
	pdd.print_params=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// PRINT_PARAMS
      if(tok=="EFERMI") {
	id=0;
	pdd.efermi=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// EFERMI
      if(tok=="SPIN") {
	id=0;
	pdd.spin=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// SPIN
      if(tok=="NLM") {
	id=0;
	pdd.nlm=atoi(aurostd::GetNextVal(val_vec[i],id).c_str());
	found_token=1;
      }// NLM
      if(!found_token) {
	cerr << "ERROR: You have input a token " << tok <<endl;
	cerr << "ERROR: This token is not recognized. Exiting! " << endl;
	exit(1);
      }
    } // for i over tokens

    // Check that all data is within bounds
    for(int i=0;i<(int)pdd.pdos_at.size();i++) { // cases
      for(int ilm=0;ilm<(int)pdd.pdos_lm[i].size();ilm++) {
	int lmp=pdd.pdos_lm[i][ilm];
	if(lmp<1 || lmp>pdd.nlm) {
	  cerr << "ERROR: SumPDOSFuncs/ReadSumDOSParams" << endl;
	  cerr << "ERROR: Error in case: " << i+1 << endl;
	  cerr << "ERROR: You have entered too low or high an lm number in entry: " << ilm+1 << endl;
	  cerr << "ERROR: Bad value is: " << lmp << endl;
	  cerr << "ERROR: Exiting." << endl;
	  exit(1);
	}
      }
    }
  }// end routine
}

// ***************************************************************************
// ReadSumDOSParams
// ***************************************************************************
// Reads in pdos data from PDOSinfile.
namespace pflow {
  void ReadInPDOSData(pflow::matrix<matrix<double> >& allpdos, pflow::pdosdata& pdd,
		      ifstream& infile) {
    string sdum;
    // Get natoms
    infile >> pdd.natoms;
    getline(infile,sdum);
    // Skip lines
    for(int i=0;i<4;i++) {
      getline(infile,sdum);
    }
    // Get number of bins and E Fermi
    double tmpe;
    infile >> sdum >> sdum >> pdd.nbins >> tmpe;
    getline(infile,sdum);
    if(abs(pdd.efermi+999)<1e-6) { // Reset efermi if it has special -999 value
      pdd.efermi=tmpe;
    }
    // Skip total DOS
    for(int i=0;i<pdd.nbins;i++) {
      getline(infile,sdum);
    }
    // Initialize size of allpdos[atoms][spin][lm][bin]
    matrix<double> tmp(pdd.nlm+1,pdd.nbins);
    allpdos = matrix<matrix<double> > (pdd.natoms,pdd.spin,tmp);
    // tpx
    //cout << spin << " " << pdd.nlm << " " << pdd.natoms << " " << pdd.nbins << endl;
    // Loop over each atom
    for(int ia=0;ia<pdd.natoms;ia++) {
      // Throw away first line since it is not DOS info
      getline(infile,sdum);
      // This is needed since getline just gets the end of the previous line when
      // the previous line has been read in by >>.  Only for atom 0 is the
      // previous line read in by getline, so only one more getline is needed.
      if(ia>0) {
	getline(infile,sdum);
      }
      // Loop over each bin
      for(int ib=0;ib<pdd.nbins;ib++) {
	// Read in energy into first spin and nlm array elements.
	infile >> allpdos[ia][0][0][ib];
	// Loop over each lm
	for(int ilm=1;ilm<=pdd.nlm;ilm++) {
	  // Loop over each spin
	  for(int is=0;is<pdd.spin;is++) {
	    infile >> allpdos[ia][is][ilm][ib];
	    // tpx
	    /*
	      cout << "ia ib is ilm allpdos " << ia << " "
	      << ib << " "
	      << is << " "
	      << ilm << " "
	      << allpdos[ia][is][ilm][ib] << endl;
	    */
	  } // spin
	} // lm
      } // bin
    } // at
  } // end routine
}

// ***************************************************************************
// SumPDOS
// ***************************************************************************
// Sums the PDOS accoring to the parameters specified.
namespace pflow {
  void SumPDOS(const pflow::matrix<matrix<double> >& allpdos, pflow::pdosdata& pdd) {
    if(pdd.spin==1) {
      pdd.pdos = matrix<double> (pdd.nbins,3);
    }
    else{
      pdd.pdos = matrix<double> (pdd.nbins,7);
    }
    // Set energies
    for(int ib=0;ib<pdd.nbins;ib++) {
      pdd.pdos[ib][0]=allpdos[0][0][0][ib]-pdd.efermi;
    }
    // Loop over flagged atoms and lm and sum
    for(int ic=0;ic<(int)pdd.pdos_at.size();ic++) {
      for(int ida=0;ida<(int)pdd.pdos_at[ic].size();ida++) {
	int ia=pdd.pdos_at[ic][ida]-1;
	// Error check
	if(ia<0 || ia>=pdd.natoms) {
	  cerr << "ERROR: SumPDOSFuncs/SumPDOS" << endl;
	  cerr << "ERROR: Error in case: " << ic+1 << endl;
	  cerr << "ERROR: You have entered too low or high an atom number in entry: " << ida+1 << endl;
	  cerr << "ERROR: Min/Max allowed values are: " << "1 / " << pdd.natoms << endl;
	  cerr << "ERROR: Bad value is: " << ia+1 << endl;
	  cerr << "ERROR: Exiting." << endl;
	  exit(1);
	}
	// Loop over each lm
	for(int idlm=0;idlm<(int)pdd.pdos_lm[ic].size();idlm++) {
	  int ilm=pdd.pdos_lm[ic][idlm];
	  // Error check
	  if(ilm<1 || ilm>pdd.nlm) {
	    cerr << "ERROR: SumPDOSFuncs/SumPDOS" << endl;
	    cerr << "ERROR: Error in case: " << ic+1 << endl;
	    cerr << "ERROR: You have entered too low or high an lm number in entry: " << idlm+1 << endl;
	    cerr << "ERROR: Min/Max allowed values are: " << "1 / " << pdd.nlm << endl;
	    cerr << "ERROR: Bad value is: " << ilm << endl;
	    cerr << "ERROR: Exiting." << endl;
	    exit(1);
	  }
	  // Loop over each bin
	  for(int ib=0;ib<pdd.nbins;ib++) {
	    // Do spin sum
	    //  double val; // DANE not used
	    // Note: allpdos[atoms][spin][lm][bin]
	    switch (pdd.spin) {
	    case 1:{ // non spin polarized
	      pdd.pdos[ib][1]=pdd.pdos[ib][1]+allpdos[ia][0][ilm][ib];
	      break;
	    }
	    case 2:{ // spin polarized
	      pdd.pdos[ib][1]=pdd.pdos[ib][1]+allpdos[ia][0][ilm][ib];
	      pdd.pdos[ib][2]=pdd.pdos[ib][2]+allpdos[ia][1][ilm][ib];
	      break;
	    }
	    default:
	      cerr << "ERROR: SumPDOSFuncs/SumPDOS" << endl;
	      cerr << "ERROR: Error in case: " << ic+1 << endl;
	      cerr << "ERROR: You have entered too low or high an spin number" << endl;
	      cerr << "ERROR: Min/Max allowed values are: " << "1 / 2" << endl;
	      cerr << "ERROR: Bad value is: " << pdd.spin << endl;
	      cerr << "ERROR: Exiting." << endl;
	      exit(1);
	    } // switch spin
	  } // bin
	} // lm
      } // at
    } // case
      // Get cumulative and difference values
    double dE=pdd.pdos[1][0]-pdd.pdos[0][0];
    if(pdd.spin==1) {
      pdd.pdos[0][2]=0;
      for(int ib=1;ib<pdd.nbins;ib++) {
	pdd.pdos[ib][2]=pdd.pdos[ib-1][2]+pdd.pdos[ib][1]*dE;
      }
    }
    else{
      pdd.pdos[0][3]=pdd.pdos[0][1]-pdd.pdos[0][2];
      pdd.pdos[0][4]=0;
      pdd.pdos[0][5]=0;
      pdd.pdos[0][6]=0;
      for(int ib=1;ib<pdd.nbins;ib++) {
	pdd.pdos[ib][3]=pdd.pdos[ib][1]-pdd.pdos[ib][2];
	pdd.pdos[ib][4]=pdd.pdos[ib-1][4]+pdd.pdos[ib][1]*dE;
	pdd.pdos[ib][5]=pdd.pdos[ib-1][5]+pdd.pdos[ib][2]*dE;
	pdd.pdos[ib][6]=pdd.pdos[ib-1][6]+pdd.pdos[ib][3]*dE;
      }
    }
  }
}

// ***************************************************************************
// RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBF
// ***************************************************************************

// ***************************************************************************
// Function TotalAtomDist
// ***************************************************************************
// NOTE: THIS FUNCTION IS OUTDATED BY RBPOSCAR DISP FUNC BELOW.
// This function gets the total distances between
// two structures by taking the sqrt of the sum of
// the squared distances between each atom.  
// If path_flag = N/n then the nearest
// images are used to measure distance.
namespace pflow {
  double TotalAtomDist(xstructure str, xstructure str00, const string& path_flag) {
    str=ReScale(str,1);
    str00=ReScale(str00,1);
    matrix<double> cpos=pflow::GetCpos(str);
    matrix<double> cpos00=pflow::GetCpos(str00);
    matrix<double> lat=pflow::GetLat(str);
    int nat=std::min(cpos.size(),cpos00.size());
    double dtot=0;
    for(int iat=0;iat<nat;iat++) {
      vector<double> dp=pflow::VVdiff(cpos[iat],cpos00[iat]);
      // If path_flag is n/N then the path is taken between
      // nearest images.  Otherwise path is between the atoms given.
      if(path_flag=="n" || path_flag=="N") {
	vector<double> ddp=pflow::vecC2F(lat,dp);
	for(int ic=0;ic<3;ic++) {
	  ddp[ic]=ddp[ic]-Nint(ddp[ic]);
	}
	dp=pflow::vecF2C(lat,ddp);
      }
      dtot=dtot+pflow::norm(dp)*pflow::norm(dp);
    }
    return sqrt(dtot);
  }
}

// ***************************************************************************
// Function GetRBDir
// ***************************************************************************
// This function gets the directories of the rubber
// band run.  The directories are assumed to
// go from 00 to 0(nim+1) (or 00 to (nim+1) if nim>8).
namespace pflow {
  vector<string> GetRBDir(const int& nim) {
    vector<string> rbdir(nim+2);
    for(int im=0;im<nim+2;im++) {
      ostringstream tmp;
      if(im<10) {
	tmp << "0" << im << ends;
      }
      else{
	tmp << im << ends;
      }
      rbdir[im]=tmp.str();
    }
    return rbdir;
  }
}

// ***************************************************************************
// Function GetRBEner
// ***************************************************************************
// This function gets the energy from the directories
// of a rubber band run. The energies are looked for
// in the OSZICAR file in the line before the last line.
namespace pflow {
  vector<double> GetRBEner(const int& nim) {
    vector<double> ener((nim+2),0.0);
    std::string tmpfile=aurostd::TmpFileCreate();
    vector<string> rbdir=GetRBDir(nim);
    for(int im=0;im<nim+2;im++) {
      ostringstream cmd;
      cmd << "cd " << rbdir[im] << ";";
      cmd << "tail -2 OSZICAR | head -1 | awk '{if($1==\"CG\") {print $4;}} {if($1==\"DAV:\") {print $3;}} {if($1==\"RMM:\") {print $3;}}' >> " << tmpfile << ";" 	<< "cd ..;" << ends;
      system(cmd.str().c_str());
      ifstream ifaustream(tmpfile.c_str());
      ifaustream >> ener[im];
      aurostd::RemoveFile(tmpfile);
    }
    return ener;
  }
}

// ***************************************************************************
// Function GetRBStruct
// ***************************************************************************
// This function gets the structures from a rubber
// band run.  The directories are assumed to go from 00 to 0(nim+1)
// (or 00 to (nim+1) if nim>8).  
namespace pflow {
  vector<xstructure> GetRBStruct(const int& nim) {
    vector<xstructure> vstr(nim+2);
    vector<string> rbdir=GetRBDir(nim);
    for(int im=0;im<nim+2;im++) {
      ostringstream file;
      if(im==0 || im==(nim+1)) { // First and last files are POSCAR
	file << rbdir[im] << "/POSCAR" << ends;
      }
      else{
	file << rbdir[im] << "/CONTCAR" << ends;
      }
      ifstream ifaustream(file.str().c_str());
      ifaustream >> vstr[im];
      ifaustream.close();
      ifaustream.clear();
    }
    return vstr;
  }
}

// ***************************************************************************
// Function GetRBDistCum
// ***************************************************************************
// This function gets the distances between images
// from the directories of a rubber band run. The initial distance
// for image i is taken by calculating the sum of the Euclidean
// distances between all atoms from image i to image i-1.
// The distances are then summed, so that d[i] is the cumulative
// distance to structure i.
namespace pflow {
  vector<double> GetRBDistCum(const vector<xstructure>& vstr, const string& path_flag) {
    int nim = vstr.size()-2;
    vector<double> d((nim+2),0.0);
    xstructure str=vstr[0];
    xstructure last_str=vstr[0];
    xstructure diffstr;
    matrix<double> cm;

    for(int im=0;im<nim+2;im++) {
      str=vstr[im];
      double dist=0.0;
      RBPoscarDisp(last_str,str,diffstr,dist,cm,path_flag);
      d[im]=dist;
      last_str=str;
    }

    // Turn d into an cumulative distance.
    for(int i=1;i<(int) d.size();i++) {
      d[i]=d[i]+d[i-1];
    }
    return d;
  }
}

// ***************************************************************************
// Function GetRBDistFromStrI
// ***************************************************************************
// This function gets the distances between images
// from the directories of a rubber band run and a given
// structure, I (dist[j] = strI - str[j]).
namespace pflow {
  vector<double> GetRBDistFromStrI(const vector<xstructure>& vstr,
				   const xstructure& strI,
				   const string& path_flag) {
    int nim = vstr.size()-2;
    vector<double> d((nim+2),0.0);
    xstructure diffstr;
    matrix<double> cm;
    for(int im=0;im<nim+2;im++) {
      RBPoscarDisp(vstr[im],strI,diffstr,d[im],cm,path_flag);
    }
    return d;
  }
}

// ***************************************************************************
//  Function RBPoscarDisp
// ***************************************************************************
// This function gets the displacement between two POSCAR files (str2-str1).  
namespace pflow {
  void RBPoscarDisp(const xstructure& str1in, const xstructure& str2in,
		    xstructure& diffstr, double& totdist, matrix<double>& cm,
		    const string& path_flag) {
    diffstr=str1in;
    xstructure str1=str1in;
    xstructure str2=str2in;
    str1=ReScale(str1,1);
    str2=ReScale(str2,1);
    matrix<double> cpos1=pflow::GetCpos(str1);
    matrix<double> cpos2=pflow::GetCpos(str2);
    matrix<double> lat1=pflow::GetLat(str1);
    matrix<double> lat2=pflow::GetLat(str1);
    matrix<double> latdiff(3,3);
    cm = matrix<double> (2);
    int nat=min(cpos1.size(),cpos2.size());
    matrix<double> cposdiff(nat,3);
    totdist=0;
    for(int ic=0;ic<3;ic++) {
      latdiff[ic]=pflow::VVdiff(lat2[ic],lat1[ic]);
    }
    for(int iat=0;iat<nat;iat++) {
      vector<double> dp=pflow::VVdiff(cpos2[iat],cpos1[iat]);
      // If path_flag is n/N then the path is taken between
      // nearest images.  Otherwise path is between the atoms given.
      if(path_flag=="n" || path_flag=="N") {
	vector<double> ddp=pflow::vecC2F(lat1,dp);
	for(int ic=0;ic<3;ic++) {
	  ddp[ic]=ddp[ic]-Nint(ddp[ic]);
	}
	dp=pflow::vecF2C(lat1,ddp);
      }
      totdist=totdist+pflow::norm(dp)*pflow::norm(dp);
      cposdiff[iat]=dp;
      // tpx I think
      //    Vout(cposdiff[iat],cout);
      //    cout << "  " << iat << endl;
    }
    totdist=sqrt(totdist);
    diffstr=pflow::SetAllAtomPos(diffstr,cposdiff,1);
    //  diffstr.SetLat(latdiff); // Just keep lat=lat1.
    cm[0]=xvector2vector(GetMom1(str1));
    cm[1]=xvector2vector(GetMom1(str2));
  }
}

// ***************************************************************************
// CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUN
// ***************************************************************************

// **************************************************
// Function CompAperpB
// **************************************************
// This function returns the component of A perpendicular to B.

namespace pflow {
  vector<double> CompAperpB(const vector<double>& a, const vector<double>& b) {
    double dp;
    dp=pflow::VVprod(a,b);
    double nb=pflow::norm(b);
    return pflow::VVdiff(a,pflow::SVprod(dp/(nb*nb),b));
  }
}

namespace pflow {
  vector<double> GetDispToAtomDir(const xstructure& a, const vector<double>& in_fpos, const int at_num) {
    // returns in_fpos-(fpos of atom at_num) in direct coordinates.  in_fpos should be in direct coords.
    vector<double> disp(3,0.0);
    for(int ic=0;ic<3;ic++) {
      //   disp[ic]=in_fpos[ic]-fpos[at_num][ic]; // CONVASP
      disp[ic]=in_fpos[ic]-a.atoms.at(at_num).fpos(ic+1); // AFLOW STRUCTURE
      //  if(disp[ic]>0.5) disp[ic]=disp[ic]-1;
      //  if(disp[ic]<-0.5) disp[ic]=disp[ic]+1;
    }
    return disp;
  }
  double GetDistToAtomDir(const xstructure& a, const vector<double>& in_fpos, const int at_num) {
    // returns norm(in_fpos-(pos of atom at_num)).  in_fpos should be in direct coords.
    vector<double> disp(3);
    disp=GetDispToAtomDir(a,in_fpos,at_num);
    disp=pflow::vecF2C(xmatrix2matrix(a.lattice),disp);
    return a.scale*pflow::norm(disp);
  }
  vector<double> GetDispToAtomImageDir(const xstructure& a,const vector<double>& in_fpos, const int at_num) {
    // returns {in_fpos}-{fpos of atom at_num} in direct coordinates, where the difference is taken
    // between the nearest images (braces denote the fact that each set point is really a set of points
    // in different image cells.  in_fpos should be in direct coords.
    vector<double> disp(3,0.0);
    for(int ic=0;ic<3;ic++) {
      //      disp[ic]=in_fpos[ic]-fpos[at_num][ic]; // CONVASP
      disp[ic]=in_fpos[ic]-a.atoms.at(at_num).fpos(ic+1); // AFLOW
      disp[ic]=disp[ic]-(int)disp[ic];
      if(disp[ic]>0.5) disp[ic]=disp[ic]-1;
      if(disp[ic]<-0.5) disp[ic]=disp[ic]+1;
    }
    return disp;
  }
  double GetDistToAtomImageDir(const xstructure& a,const vector<double>& in_fpos, const int at_num) {
    // returns {in_fpos}-{fpos of atom at_num} in direct coordinates, where the difference is taken
    // between the nearest images (braces denote the fact that each set point is really a set of points
    // in different image cells.  in_fpos should be in direct coords.
    vector<double> disp(3);
    disp=pflow::GetDispToAtomImageDir(a,in_fpos,at_num);
    disp=pflow::vecF2C(xmatrix2matrix(a.lattice),disp);
    return a.scale*norm(disp);
  }
}

namespace pflow {
  void GetChgInt(vector<matrix<double> >& rad_chg_int, matrix<double>& vor_chg_int) {
    // Read in CHGCAR
    // Read in intial POSCAR format at beginning of file
    xstructure str;
    cin >> str;
    // Read in the grid for total charge
    vector<int> ngrid(3);
    int npts;
    cin >> ngrid[0] >> ngrid[1] >> ngrid[2];
    npts=ngrid[0]*ngrid[1]*ngrid[2];
    // Read in the total charge
    vector<double> chg_tot(npts,0.0);
    for(int i=0;i<npts;i++) {
      cin >> chg_tot[i];
    }
    // Read in dummy fields
    // If first field is augmentation then read 4*natoms+natoms-1 more fields. For vasp 4.4.5.
    // If first field is not augmentation then read natoms-1 more fields. For vasp 4.4.1.
    string ds;
    int natoms=pflow::GetNumAtoms(str);
    cin >> ds;
    string keyword = "augmentation";
    if(ds==keyword) {
      for(int i=0;i<(5*natoms-1);i++) {
	cin >> ds;
      }
    }
    else{
      for(int i=0;i<(natoms-1);i++) {
	cin >> ds;
      }
    }

    // Read in the grid for diff charge
    cin >> ngrid[0] >> ngrid[1] >> ngrid[2];
    npts=ngrid[0]*ngrid[1]*ngrid[2];
    // Read in the diff charge
    vector<double> chg_diff(npts,0.0);
    for(int i=0;i<npts;i++) {
      cin >> chg_diff[i];
    }
    // Read in 4*natoms=4*natoms dummy fields.  This is
    // not necessary and the lines only exist in vasp 4.4.5.
    /*
      for(int i=0;i<4*natoms;i++) {
      cin >> ds;
      }
    */

    // Loop over grid and bin charge densities.

    // Radial charge integration.
    // Integral over a sphere will be done by setting up NRBIN bins, each
    // of witdh DR, our to RMAX.  Then each charge at each point will be put is the
    // appropriate bin for each atom.  Note that the points and all possible relevant
    // images must be looped over.
    double DR=0.05; // width of radial bins.
    double RMAX=3; // Max distance to bin to.
    int NRBIN=(int)(RMAX/DR); // First bin is 0->DR, last bin is RMAX-DR->RMAX
    matrix<double> dum_mat(NRBIN,5,0.0);
    rad_chg_int = vector<matrix<double> > (natoms,dum_mat);
    vector<int> ig(3);
    // Determine the number of images along each lattice vector which must be considered.
    // You will center the charge density on each atom so the max number of images along a
    // given lattice param should be no more than however many lattice params are needed
    // to make sure that you are more than RMAX from an atom along that lattice param direction.
    vector<int> imax(3);
    double scale=pflow::GetScale(str);
    matrix<double> lat=pflow::GetLat(str);
    for(int ic=0;ic<3;ic++) {
      imax[ic]=(int)(RMAX/(scale*pflow::norm(lat[ic]))+1)*ngrid[ic];
    }
    //    double sumchg=0;  // DANE not used
    // Atom loop
    matrix<double> at_fpos=pflow::GetFpos(str);
    for(int iat=0;iat<natoms;iat++) {
      // Initialize radial positions in rad_chg_int    
      for(int ib=0;ib<NRBIN;ib++) {
	rad_chg_int[iat][ib][0]=(ib+1)*DR;
      }
      // Find location of atom on grid.
      vector<int> at_grid_pos(3);
      for(int ic=0;ic<3;ic++) {
	at_grid_pos[ic]=Nint(at_fpos[iat][ic]*ngrid[ic]);
      }
      // Grid loop: This loops over grid points shifted to center around grid point approximation to
      // position of iat.  The ig values therefore need to be modified to be true grid points
      // (true_ig) by undoing the shift (adding back the at_grid_pos).
      for(ig[0]=-imax[0];ig[0]<imax[0];ig[0]++) {
	for(ig[1]=-imax[1];ig[1]<imax[1];ig[1]++) {
	  for(ig[2]=-imax[2];ig[2]<imax[2];ig[2]++) {
	    // Get true grid points
	    vector<int> true_ig(3);
	    vector<double> chg_fpos(3);
	    for(int ic=0;ic<3;ic++) {
	      true_ig[ic]=ig[ic]+at_grid_pos[ic];
	      // Get true position of this grid point in direct coords.
	      chg_fpos[ic]=((float)true_ig[ic]/(float)ngrid[ic]);
	      // Shift true grid point back into the 000 cell.
	      true_ig[ic]=true_ig[ic]-(int)((float)true_ig[ic]/(float)ngrid[ic])*ngrid[ic];
	      if(true_ig[ic]<0) true_ig[ic]=true_ig[ic]+ngrid[ic];
	    } // ic
	    // Get charge values for each grid position
	    int id=true_ig[0]+true_ig[1]*ngrid[0]+true_ig[2]*ngrid[0]*ngrid[1];
	    double chgtot=chg_tot[id];
	    double chgdiff=chg_diff[id];
	    double chgup=0.5*(chg_tot[id]+chg_diff[id]);
	    double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
	    vector<double> fpos(3);
	    // Get distance from iat to true grid point position.
	    double dist=pflow::GetDistToAtomDir(str,chg_fpos,iat);
	    // Get bin for this atom and increment radial charge density.
	    int ibin=(int)(dist/DR);	  
	    if(ibin<NRBIN) {
	      rad_chg_int[iat][ibin][1]=rad_chg_int[iat][ibin][1]+chgtot;
	      rad_chg_int[iat][ibin][2]=rad_chg_int[iat][ibin][2]+chgdiff;
	      rad_chg_int[iat][ibin][3]=rad_chg_int[iat][ibin][3]+chgup;
	      rad_chg_int[iat][ibin][4]=rad_chg_int[iat][ibin][4]+chgdn;
	    } // if ibin<NRBIN
	  } // ig[2]
	} // ig[1]
      } // ig[0]
    } // iat

    // Voronoi charge integration
    vor_chg_int = matrix<double> (natoms,4,0.0);
    for(ig[0]=0;ig[0]<ngrid[0];ig[0]++) {
      for(ig[1]=0;ig[1]<ngrid[1];ig[1]++) {
	for(ig[2]=0;ig[2]<ngrid[2];ig[2]++) {
	  // Get charge values for this grid position
	  int id=ig[0]+ig[1]*ngrid[0]+ig[2]*ngrid[0]*ngrid[1];
	  double chgtot=chg_tot[id];
	  double chgdiff=chg_diff[id];
	  double chgup=0.5*(chg_tot[id]+chg_diff[id]);
	  double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
	  vector<double> fpos(3);
	  // Set up grid positions in direct coordinates.
	  for(int ic=0;ic<3;ic++) {
	    fpos[ic]=((float)ig[ic])/(float)ngrid[ic];
	  } // ic
	  // Loop over atoms
	  double dist;
	  double min_dist=pflow::GetDistToAtomImageDir(str,fpos,0);
	  int min_dist_id=0;
	  for(int iat=0;iat<natoms;iat++) {
	    // Get "image" distance to each atom (distance to closest image)
	    dist=pflow::GetDistToAtomImageDir(str,fpos,iat);
	    // Keep track of which atom is the minimum "image" distance away
	    if(dist<min_dist) {
	      min_dist=dist;
	      min_dist_id=iat;
	    } // if dist<min_dist
	  } // iat
	  vor_chg_int[min_dist_id][0]=vor_chg_int[min_dist_id][0]+chgtot;
	  vor_chg_int[min_dist_id][1]=vor_chg_int[min_dist_id][1]+chgdiff;
	  vor_chg_int[min_dist_id][2]=vor_chg_int[min_dist_id][2]+chgup;
	  vor_chg_int[min_dist_id][3]=vor_chg_int[min_dist_id][3]+chgdn;
	} // ig[2]
      } // ig[1]
    } // ig[0]

    // Normalize and turn rad_chg_int from chg in each radial shell to cumulative integral.
    for(int iat=0;iat<natoms;iat++) {
      for(int i=0;i<4;i++) {
	vor_chg_int[iat][i]=vor_chg_int[iat][i]/npts;
      }
      for(int i=1;i<5;i++) {
	for(int ib=0;ib<NRBIN;ib++) {
	  rad_chg_int[iat][ib][i]=rad_chg_int[iat][ib][i]/npts;
	  if(ib>0) rad_chg_int[iat][ib][i]=(rad_chg_int[iat][ib][i]+rad_chg_int[iat][ib-1][i]);
	} // ib
      } // for loop over tot,diff,up,dn
    } // iat
  }
}

namespace pflow {
  void pd_params::Print(ostream& outf) const{
    outf << type << " # Type " << endl;
    outf << scale << " # scale " << endl;
    for(int ic=0;ic<3;ic++) {
      for(int jc=0;jc<3;jc++) {
	outf << pts[ic][jc] << " ";
      }
      outf << " # pts " << endl;
    }
    for(int ic=0;ic<3;ic++) {
      for(int jc=0;jc<3;jc++) {
	outf << dpts[ic][jc] << " ";
      }
      outf << " # dpts " << endl;
    }
    outf << orig_loc << " # origin_loc " << endl;
    outf << ortho << " # ortho " << endl;
  }
}

// ***************************************************************************
// pflow::ReadCHGCAR
// ***************************************************************************
namespace pflow {
  bool ReadCHGCAR(xstructure& str,
      stringstream& chgcar_header,
      vector<int>& ngrid,
      vector<int>& format_dim,    
      vector<double>& chg_tot,
      vector<double>& chg_diff,
      stringstream& chgcar_ss,
      ostream& oss) {
    // This funtion reads in the CHGCAR or AECCAR file from a file.
    // Operates under assumption that it was formatted by v46s (updated by Corey)
    // v46s is standard for AFLOWLIB
    // Format:
    //      Starts with POSCAR
    //      newline
    //      NX NY NZ
    //      NX*NY*NZ entries with a newline after each one (grouped by 5 usually, but 
    //          this code can read any grouping) (Tot dens = pUP+pDn)
    //      PAW augmentation occupancies, one for each atom
    //      If spin-polarized:
    //          0's line
    //          NX NY NZ
    //          NX*NY*NZ entries with a newline after each one (grouped by 5) 
    //              (Diff dens = pUP-pDN)
    //          PAW augmentation occupancies, one for each atom
    // AECCAR's will stop after total density
    //
    // PAW augmentation occupancies are not used at all in AFLOW, and all other ReadCHGCAR 
    // scripts simply ignore it
    // there also seems to be some major discrepancies between how different versions of 
    // VASP format the CHGCAR
    // OBSOLETE:
    // here's text for previous versions (OBSOLETE):
    //      If first field is augmentation then read 4*natoms+natoms-1 more fields. For vasp 4.4.5.
    //      If first field is not augmentation then read natoms-1 more fields. For vasp 4.4.1.
    // previous version text (OBSOLETE):
    //      Here is the format of a CHGCAR file (for v445 with
    //      spin polarization from DEC linux).
    //      Note that for Intel linux vasp445 the newlines after
    //      each chg density line seem to missing, however, the
    //      code still reads this in correctly.
    //      newline
    //      Starts with a POSCAR file
    //      newline
    //      NX NY NZ
    //      NX*NY*NZ entries with a newline after each one (grouped by 5) (Tot dens = pUP+pDn)
    //      Natom lines, each with 4 fields, with a newline after each one
    //      Natom entries (grouped by 5).
    //      NX NY NZ
    //      NX*NY*NZ entries with a newline after each one (grouped by 5) (Diff dens = pUP-pDN)
    //      Natom lines, each with 4 fields, with a newline after each one
    //      Natom entries (grouped by 5). // I don't think these are actually in the CHGCAR.

    string soliloquy="pflow::ReadCHGCAR():  ";    // so you know who's talking

    // DEBUG
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << soliloquy << "BEGIN" << endl;

    string content,line, keyword="augmentation";
    vector<string> vcontent, sum_tokens, line_tokens;
    vector<int> ngrid_compare;
    content=chgcar_ss.str();
    aurostd::string2vectorstring(content,vcontent,true,true); //consecutive=true, trim_edges=true, CO 170613, we have issues with "\n \n" vs. "\n\n", we want both to be separated
    stringstream poscar;
    chgcar_header.str("");
    int grid_point=0;
    uint natoms=0, linecount=0, npts=1, npts_compare=1, numcolumns=0, num_skip=0;
    bool got_chgtot=false;      //so we know that we're now reading chgdiff
    if(LDEBUG) oss << soliloquy << "CHECKING FORMAT OF SYSTEM NAME AND SCALING FACTOR" << endl;
    //do some quick checks to make sure it's format correctly, otherwise kill the script
    //first line is comment line, skip check
    //next line contains scaling factor
    aurostd::string2tokens(vcontent.at(1),line_tokens," ");
    if(line_tokens.size()!=1) {
      oss << endl;
      oss << soliloquy << "ERROR: Unrecognized CHGCAR header format (scaling factor)." << endl;
      oss << endl;
      return FALSE;
    }
    line_tokens.clear();
    if(LDEBUG) oss << soliloquy << "FORMAT FOR SCALING FACTOR LOOKS OK" << endl;

    //next three lines should contain lattice vectors
    if(LDEBUG) oss << soliloquy << "CHECKING FORMAT FOR LATTICE VECTORS" << endl;
    for(uint i=2;i<5;i++) {
      aurostd::string2tokens(vcontent.at(i),line_tokens," ");
      if(line_tokens.size()!=3) {
        oss << endl;
        oss << soliloquy << "ERROR: Unrecognized CHGCAR header format (lattice vectors)." << endl;
        oss << endl;
        return FALSE;
      }
    }
    line_tokens.clear();
    if(LDEBUG) oss << soliloquy << "FORMAT FOR LATTICE VECTORS LOOKS OK" << endl;

    //if we get here, assume CHGCAR is formatted correctly for now
    //read in header, scaling factor, lattice vectors, number of atoms, coordinate type
    if(LDEBUG) oss << soliloquy << "READING POSCAR" << endl;
    for(uint i=0;i<7;i++) {
      poscar << vcontent.at(i) << endl;
      chgcar_header << vcontent.at(i) << endl;
      linecount=i;
    }
    //read in number of atoms
    aurostd::string2tokens(vcontent.at(5),sum_tokens," ");
    for(uint i=0;i<sum_tokens.size();i++) {
      natoms+=aurostd::string2utype<uint>(sum_tokens.at(i));
    }
    linecount++;   //start on coordinate line
    //read in coordinates
    for(uint i=0;i<natoms;i++) {
      poscar << vcontent.at(linecount) << endl;
      chgcar_header << vcontent.at(linecount) << endl;
      linecount++;
    }
    //define xstruture
    if(LDEBUG) oss << soliloquy << "DEFINE XSTRUCTURE" << endl;
    str=poscar;
    if(LDEBUG) oss << str << endl;
    //skip empty line
    linecount++;
    //now we get chg_tot and chg_diff (if it exists)
    if(LDEBUG) oss << soliloquy << "READ CHG VALUES" << endl;
    while(vcontent.size()>linecount) {
      //see if it's augmentation occupancies that we're reading, and skip appropriate number of lines
      if(vcontent.at(linecount).compare(0,12,keyword)==0) {
        if(LDEBUG) oss << soliloquy << "GET FORMAT FOR AUGMENTATION OCCUPANCIES" << endl;
        for(uint i=0;i<natoms;i++) {
          //check formatting of line first, 'augmentation occupancies natom npts'
          aurostd::string2tokens(vcontent.at(linecount),line_tokens," ");
          if(line_tokens.size()!=4) {
            oss << endl;
            oss << soliloquy << "ERROR: Unrecognized CHGCAR format (augmentation occupancies)." << endl;
            oss << endl;
            return FALSE;
          }
          npts=aurostd::string2utype<int>(line_tokens.at(3));
          line_tokens.clear();
          aurostd::string2tokens(vcontent.at(linecount+1),line_tokens," ");
          numcolumns=line_tokens.size();
          // this npts (augmentation occupancies) will be different than chgtot, grab it for comparsion
          format_dim.push_back(npts);
          format_dim.push_back(numcolumns);
          num_skip=(uint)(std::ceil((double)npts/(double)numcolumns));    //skip right number of liens
          linecount+=num_skip+1;  //+1 to get to next augmentation line
          if(LDEBUG) oss << soliloquy << "SUCCESSFULLY GATHERED FORMAT FOR AUGMENTATION OCCUPANCIES" << endl;
        }
      } else {    //must be chg_tot or chg_diff
        npts=1; ngrid.clear();      //clear for next round
        line_tokens.clear();
        //if gathering chg_diff, we have to ignore extra line of 0's
        if(got_chgtot) {
          linecount++;
        }
        //get npts
        sum_tokens.clear();
        aurostd::string2tokens(vcontent.at(linecount),sum_tokens, " ");
        //check formatting, this line should contain npts
        if(sum_tokens.size()!=3) {
          oss << endl;
          oss << soliloquy << "ERROR: Unrecognized CHGCAR format (number of grid points)." << endl;
          oss << endl;
          return FALSE;
        }
        for(uint i=0;i<sum_tokens.size();i++) {
          grid_point=aurostd::string2utype<int>(sum_tokens.at(i));   
          npts=npts*grid_point;
          ngrid.push_back(grid_point);
        }
        //initialize once
        if(!got_chgtot) {
          if(LDEBUG) oss << soliloquy << "READING CHG_DIFF" << endl;
          //now with npts, initialize chg_tot and chg_diff
          chg_tot=vector<double>(npts,0.0);
          ngrid_compare=ngrid;
          npts_compare=npts;  //save for comparison later with chg_diff
        } else {
          if(LDEBUG) oss << soliloquy << "READING CHG_TOT" << endl;
          chg_diff=vector<double>(npts,1.0);
          if(npts!=npts_compare) {    //chg_tot and chg_diff have different number of points
            oss << endl;
            oss << soliloquy << "ERROR: Number of grid points for CHG_tot and CHG_diff do not match. " << endl;
            oss << soliloquy << "ERROR: npts_total: " << npts << endl;
            oss << soliloquy << "ERROR: npts_diff: " << npts_compare << endl;
            oss << soliloquy << "ERROR: This will give nonsense!! " << endl;
            oss << endl;
            return FALSE;
          }
          if(ngrid!=ngrid_compare) {  //chg_tot and chg_diff have different grid
            oss << endl;
            oss << soliloquy << "ERROR: Grids for CHG_tot and CHG_diff do not match. " << endl;
            oss << soliloquy << "ERROR: ngrid_total: " << ngrid.at(0) << " " << ngrid.at(1) << " " << ngrid.at(2) << endl;
            oss << soliloquy << "ERROR: ngrid_diff: " << ngrid_compare.at(0) << " " << ngrid_compare.at(1) << " " << ngrid_compare.at(2) << endl;
            oss << soliloquy << "ERROR: This will give nonsense!! " << endl;
            oss << endl;
            return FALSE;
          }
        }
        //assume VASP CHGCAR format (number of columns) is constant, grab number
        linecount++;
        aurostd::string2tokens(vcontent.at(linecount),line_tokens," ");
        numcolumns=line_tokens.size();
        //already grab npts for comparsion, we just need numcolumns
        format_dim.push_back(numcolumns);

        for(uint i=0;i<(uint) npts;i+=numcolumns) {
          line_tokens.clear();
          aurostd::string2tokens(vcontent.at(linecount),line_tokens," ");
          numcolumns=line_tokens.size();  //there may be lines at end with only 1 or 2 entries
          for(uint j=0;j<numcolumns;j++) {
            if(!got_chgtot) {
              chg_tot.at(i+j)=aurostd::string2utype<double>(line_tokens.at(j));
            } else {
              chg_diff.at(i+j)=aurostd::string2utype<double>(line_tokens.at(j));
            }
          }
          linecount++;
        }
        if(!got_chgtot) {
          got_chgtot=true;
        }
        line_tokens.clear();
        if(LDEBUG) oss << soliloquy << "SUCCESSFULLY READ CHG VALUES" << endl;
      }
    }
    if(LDEBUG) oss << soliloquy << "DONE" << endl;
    return TRUE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ReadCHGCAR
// ***************************************************************************
namespace pflow {
  bool ReadChg(xstructure& str,
      vector<int>& ngrid,
      vector<double>& chg_tot,
      vector<double>& chg_diff,
      istream& chgfile) {
    //OLD READCHG FORMAT, just keeping it to not screw up other functions that depend on it
    stringstream chgcar_header,chgcar_ss;
    vector<int> format_dim;
    chgcar_ss << chgfile.rdbuf();
    ostringstream oss;
    return ReadCHGCAR(str,chgcar_header,ngrid,format_dim,chg_tot,chg_diff,chgcar_ss,oss);
  }
} // namespace pflow

// ***************************************************************************
// GetChgInt
// ***************************************************************************
//  This funtion reads in the CHGCAR file from standard input
//  and returns a matrix of charge integrals in the Voronoi volume around
//  each atom and a vector of matrices which gives, for each atom,
//  the integrated charges in a sphere vs. the radius of the sphere.
namespace pflow {
  void GetChgInt(vector<matrix<double> >& rad_chg_int,
		 matrix<double>& vor_chg_int,
		 xstructure& str,
		 vector<int>& ngrid,
		 vector<double>& chg_tot,
		 vector<double>& chg_diff) {
    // Loop over grid and bin charge densities.

    // Radial charge integration.
    // Integral over a sphere will be done by setting up NRBIN bins, each
    // of witdh DR, our to RMAX.  Then each charge at each point will be put is the
    // appropriate bin for each atom.  Note that the points and all possible relevant
    // images must be looped over.
    double DR=0.05; // width of radial bins.
    double RMAX=3; // Max distance to bin to.
    int NRBIN=(int)(RMAX/DR); // First bin is 0->DR, last bin is RMAX-DR->RMAX
    matrix<double> dum_mat(NRBIN,5,0.0);
    int natoms=pflow::GetNumAtoms(str);
    int npts=ngrid[0]*ngrid[1]*ngrid[2];
    rad_chg_int = vector<matrix<double> > (natoms,dum_mat);
    vector<int> ig(3);
    // Determine the number of images along each lattice vector which must be considered.
    // You will center the charge density on each atom so the max number of images along a
    // given lattice param should be no more than however many lattice params are needed
    // to make sure that you are more than RMAX from an atom along that lattice param direction.
    vector<int> imax(3);
    double scale=pflow::GetScale(str);
    matrix<double> lat=pflow::GetLat(str);
    for(int ic=0;ic<3;ic++) {
      imax[ic]=(int)(RMAX/(scale*norm(lat[ic]))+1)*ngrid[ic];
    }
    //   double sumchg=0;   // DANE not used
    // Atom loop
    matrix<double> at_fpos=pflow::GetFpos(str);
    for(int iat=0;iat<natoms;iat++) {
      // Initialize radial positions in rad_chg_int    
      for(int ib=0;ib<NRBIN;ib++) {
	rad_chg_int[iat][ib][0]=(ib+1)*DR;
      }
      // Find location of atom on grid.
      vector<int> at_grid_pos(3);
      for(int ic=0;ic<3;ic++) {
	at_grid_pos[ic]=Nint(at_fpos[iat][ic]*ngrid[ic]);
      }
      // Grid loop: This loops over grid points shifted to center around grid point approximation to
      // position of iat.  The ig values therefore need to be modified to be true grid points
      // (true_ig) by undoing the shift (adding back the at_grid_pos).
      for(ig[0]=-imax[0];ig[0]<imax[0];ig[0]++) {
	for(ig[1]=-imax[1];ig[1]<imax[1];ig[1]++) {
	  for(ig[2]=-imax[2];ig[2]<imax[2];ig[2]++) {
	    // Get true grid points
	    vector<int> true_ig(3);
	    vector<double> chg_fpos(3);
	    for(int ic=0;ic<3;ic++) {
	      true_ig[ic]=ig[ic]+at_grid_pos[ic];
	      // Get true position of this grid point in direct coords.
	      chg_fpos[ic]=((float)true_ig[ic]/(float)ngrid[ic]);
	      // Shift true grid point back into the 000 cell.
	      true_ig[ic]=true_ig[ic]-(int)((float)true_ig[ic]/(float)ngrid[ic])*ngrid[ic];
	      if(true_ig[ic]<0) true_ig[ic]=true_ig[ic]+ngrid[ic];
	    } // ic
	      // Get charge values for each grid position
	    int id=true_ig[0]+true_ig[1]*ngrid[0]+true_ig[2]*ngrid[0]*ngrid[1];
	    double chgtot=chg_tot[id];
	    double chgdiff=chg_diff[id];
	    double chgup=0.5*(chg_tot[id]+chg_diff[id]);
	    double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
	    vector<double> fpos(3);
	    // Get distance from iat to true grid point position.
	    double dist=pflow::GetDistToAtomDir(str,chg_fpos,iat);
	    // Get bin for this atom and increment radial charge density.
	    int ibin=(int)(dist/DR);	  
	    if(ibin<NRBIN) {
	      rad_chg_int[iat][ibin][1]=rad_chg_int[iat][ibin][1]+chgtot;
	      rad_chg_int[iat][ibin][2]=rad_chg_int[iat][ibin][2]+chgdiff;
	      rad_chg_int[iat][ibin][3]=rad_chg_int[iat][ibin][3]+chgup;
	      rad_chg_int[iat][ibin][4]=rad_chg_int[iat][ibin][4]+chgdn;
	    } // if ibin<NRBIN
	  } // ig[2]
	} // ig[1]
      } // ig[0]
    } // iat

      // Voronoi charge integration
    vor_chg_int = matrix<double> (natoms,4,0.0);
    for(ig[0]=0;ig[0]<ngrid[0];ig[0]++) {
      for(ig[1]=0;ig[1]<ngrid[1];ig[1]++) {
	for(ig[2]=0;ig[2]<ngrid[2];ig[2]++) {
	  // Get charge values for this grid position
	  int id=ig[0]+ig[1]*ngrid[0]+ig[2]*ngrid[0]*ngrid[1];
	  double chgtot=chg_tot[id];
	  double chgdiff=chg_diff[id];
	  double chgup=0.5*(chg_tot[id]+chg_diff[id]);
	  double chgdn=0.5*(chg_tot[id]-chg_diff[id]);
	  vector<double> fpos(3);
	  // Set up grid positions in direct coordinates.
	  for(int ic=0;ic<3;ic++) {
	    fpos[ic]=((float)ig[ic])/(float)ngrid[ic];
	  } // ic
	    // Loop over atoms
	  double dist;
	  double min_dist=pflow::GetDistToAtomImageDir(str,fpos,0);
	  int min_dist_id=0;
	  for(int iat=0;iat<natoms;iat++) {
	    // Get "image" distance to each atom (distance to closest image)
	    dist=pflow::GetDistToAtomImageDir(str,fpos,iat);
	    // Keep track of which atom is the minimum "image" distance away
	    if(dist<min_dist) {
	      min_dist=dist;
	      min_dist_id=iat;
	    } // if dist<min_dist
	  } // iat
	  vor_chg_int[min_dist_id][0]=vor_chg_int[min_dist_id][0]+chgtot;
	  vor_chg_int[min_dist_id][1]=vor_chg_int[min_dist_id][1]+chgdiff;
	  vor_chg_int[min_dist_id][2]=vor_chg_int[min_dist_id][2]+chgup;
	  vor_chg_int[min_dist_id][3]=vor_chg_int[min_dist_id][3]+chgdn;
	} // ig[2]
      } // ig[1]
    } // ig[0]

      // Normalize and turn rad_chg_int from chg in each radial shell to cumulative integral.
    for(int iat=0;iat<natoms;iat++) {
      for(int i=0;i<4;i++) {
	vor_chg_int[iat][i]=vor_chg_int[iat][i]/npts;
      }
      for(int i=1;i<5;i++) {
	for(int ib=0;ib<NRBIN;ib++) {
	  rad_chg_int[iat][ib][i]=rad_chg_int[iat][ib][i]/npts;
	  if(ib>0) rad_chg_int[iat][ib][i]=(rad_chg_int[iat][ib][i]+rad_chg_int[iat][ib-1][i]);
	} // ib
      } // for loop over tot,diff,up,dn
    } // iat
  }
}

// ***************************************************************************
// ReadPlaneDensParams
// ***************************************************************************
// Reads in parameters.
// Input file format:

// D # Coordinates for following points (Direct/Cartesian)
// scale # Scale factor - edges of plane get mult. by this (but not origin).
// x y z # origin point
// x y z # "X" axis
// x y z # "Y" axis
// Nx Ny # Number of X and Y grid points
// Middle # Location for origin (Middle/Corner).
// Ortho # Whether to use Y orthogonal to X (Ortho/Strict).
namespace pflow {

  void ReadPlaneDensParams(const xstructure& str, pd_params& pdp, istream& infile) {
    // Note that this converts everything to direct coordinates.
    string s;
    double TOL=1e-7;
    pdp.pts = matrix<double> (3,3);
    pdp.dpts = matrix<double> (3,3);
    infile >> pdp.type;
    getline(infile,s);
    infile >> pdp.scale;
    getline(infile,s);
    for(int ic=0;ic<3;ic++) {
      for(int jc=0;jc<3;jc++) {
	infile >> pdp.pts[ic][jc];
      }
      getline(infile,s);
      // transform points to cart coordinates if necesary.
      if(pdp.type[0]=='D' || pdp.type[0]=='d') pdp.pts[ic]=pflow::vecF2C(pflow::GetLat(str),pdp.pts[ic]);
    }
    infile >> pdp.Nx >> pdp.Ny;
    getline(infile,s);
    infile >> pdp.orig_loc;
    getline(infile,s);
    infile >> pdp.ortho;
    getline(infile,s);
    pdp.dpts[0]=pdp.pts[0];
    pdp.dpts[1]=VVdiff(pdp.pts[1],pdp.pts[0]);
    pdp.dpts[2]=VVdiff(pdp.pts[2],pdp.pts[0]);
    // Orthogonalize if necessary
    // Finds vector of same length as dpts[2] but perp. to dpts[1].
    if(pdp.ortho[0]=='O' || pdp.ortho[0]=='o') {
      double norm_orig=norm(pdp.dpts[2]);
      pdp.dpts[2]=pflow::CompAperpB(pdp.dpts[2],pdp.dpts[1]);
      double norm_new=norm(pdp.dpts[2]);
      if(norm_new<TOL) {
	cout << "ERROR: ChgFuncs/ReadPlaneDensParams" << endl;
	cout << "Your two axis are parallel and no plane can be determined." << endl;
	cout << "Exiting." << endl;
	exit(1);
      }
      pdp.dpts[2]=SVprod(norm_orig/norm_new,pdp.dpts[2]);
    }
    // Shift origin to middle if necessary.
    // I always place dpts[0] at the corner.  So to put origin
    // at the middle shift dpts[0] by -0.5*(dpts[1]+dpts[2]).
    if(pdp.orig_loc[0]=='M' || pdp.orig_loc[0]=='m') {
      vector<double> shift = pflow::SVprod(0.5,VVsum(pdp.dpts[1],pdp.dpts[2]));
      pdp.dpts[0]=pflow::VVdiff(pdp.dpts[0],shift);
    }
    // Convert everything to direct coordinates
    for(int ic=0;ic<3;ic++) {
      pdp.pts[ic]=pflow::vecC2F(pflow::GetLat(str),pdp.pts[ic]);
      pdp.dpts[ic]=pflow::vecC2F(pflow::GetLat(str),pdp.dpts[ic]);
      // rescale coordinates
      if(ic>0) pdp.dpts[ic]=pflow::SVprod(pdp.scale,pdp.dpts[ic]);
    }
  }
}

// ***************************************************************************
// GetPlaneDens
// ***************************************************************************
// Gets density in plane.  Work always within direct coordinates.
// The original densities are stored in chg_tot, chg_diff such that
// if the direct coordinates of a grid point are i/Nx,j/Ny,k/Nz then
// the index for the charge at that point is id=i+j*Nx+k*Nx*Ny.  The
// charge in the plane will be stored the same way.  Let Npx,Npy =
// plane grid.  Then if the coordinates of a grid point in the plane are
// i/Npx, j/Npy then the index for the charge at that point in id=i+j*Npx.
// You need to be able to access charge densities and structural
// parameters for this routine.
namespace pflow {
  void GetPlaneDens(const pd_params& pdp,
		    vector<double>& dens2d_tot,
		    vector<double>& dens2d_diff,
		    const xstructure& str,
		    const vector<int>& ngrid,
		    const vector<double>& chg_tot,
		    const vector<double>& chg_diff) {
    if(str.atoms.size()) {;} // phony just to keep str busy
    int Nx=pdp.Nx;
    int Ny=pdp.Ny;
    int npts=Nx*Ny;
    dens2d_tot=vector<double> (npts,0.0);
    dens2d_diff=vector<double> (npts,0.0);
    // Loop over each point in the 2d plane
    for(int i=0;i<Nx;i++) {
      for(int j=0;j<Ny;j++) {
	// Get id2d for this point.
	int id2d =i+Nx*j;
	// Get direct coordinates for this point.
	vector<double> pt(3,0.0);
	vector<double> Xpt(3,0.0);
	vector<double> Ypt(3,0.0);
	Xpt=SVprod(float(i)/float(Nx-1),pdp.dpts[1]);
	Ypt=SVprod(float(j)/float(Ny-1),pdp.dpts[2]);
	pt=VVsum(Xpt,Ypt);
	pt=VVsum(pdp.dpts[0],pt);
	// Get charge at this point.
	// Simply use value of charge at nearest point
	// on the 3d charge grid.
	int ii=Nint(ngrid[0]*pt[0])%ngrid[0];
	int jj=Nint(ngrid[1]*pt[1])%ngrid[1];
	int kk=Nint(ngrid[2]*pt[2])%ngrid[2];
	if(ii<0) ii=ii+ngrid[0];
	if(jj<0) jj=jj+ngrid[1];
	if(kk<0) kk=kk+ngrid[2];
	int id3d = ii+ngrid[0]*jj+ngrid[0]*ngrid[1]*kk;
	dens2d_tot[id2d]=chg_tot[id3d];
	dens2d_diff[id2d]=chg_diff[id3d];
      }
    }  
  }  
}

// ***************************************************************************
// PrintPlaneDens
// ***************************************************************************
namespace pflow {
  void PrintPlaneDens(const pd_params& pdp,
		      const vector<double>& dens2d_tot,
		      const vector<double>& dens2d_diff,
		      const xstructure& str) {  
    if(str.atoms.size()) {;} // phony just to keep str busy
    
    ofstream outf_tot("dens2d.tot.out");
    ofstream outf_up("dens2d.up.out");
    ofstream outf_dn("dens2d.dn.out");
    ofstream outf_diff("dens2d.diff.out");

    int Nx=pdp.Nx;
    int Ny=pdp.Ny;
    for(int i=0;i<Nx;i++) {
      for(int j=0;j<Ny;j++) {
	int id=i+Nx*j;
	outf_tot << dens2d_tot[id] << " ";
	outf_diff << dens2d_diff[id] << " ";
	outf_up << (dens2d_tot[id]+dens2d_diff[id])/2.0 << " ";
	outf_dn << (dens2d_tot[id]-dens2d_diff[id])/2.0 << " ";
      }
      outf_tot << endl;
      outf_diff << endl;
      outf_up << endl;
      outf_dn << endl;
    }
  }
}

// ***************************************************************************
// EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWA
// ***************************************************************************
double CONV_FACT=(1.0E+10*E_ELECTRON/(4.0*PI*EPS_VACUUM));

// ***************************************************************************
// Ewald
// ***************************************************************************
// Calculates Ewald Sum.  Based on routine from Eric Wu.

// ******************************************************************
// **                   Subroutine ewald                           **
// **                                                              **
// **   Subroutine ewald calculates the ewald sum in eV            **
// **   Also, this calculates the ewald gradients in ev   since    **
// **         it is just a small calculation to do- ewald scales   **
// **         as N^2, extra force calculation  adds as  N          **
// **                                                              **
// **   NB if you want to convert forces to ev/A you need the      **
// **   inverse of the lattice vectors matrix (A) then do a        **
// **    Fx=ainvx*fx+ainvy*fy+ainvz*fz                             **
// **                 where Fx is in ev/A fx is reduced            **
// **   to get the components along each cooredinate (I think)     **
// **                                                              **
// **                                                              **
// **   References:                                                **
// **     www.dl.ac.uk/TCS/Software/DL_POLY/USRMAN/node62.html     **
// **     www.ee.duke.edu/~ayt/ewaldpaper/node8.html               **
// **     Moldy (molecular dynamics) manual                        **
// **   http://www.earth.ox.ac.uk/~keith/moldy-manual/node11.html  **
// **     Allen+tildesly p159 (for cubic)                          **
// **     Kittel Appendix B                                        **
// **                                                              **
// **   Formula:                                                   **
// **                                                              **
// **   For the energy:                                            **
// **      Etotal = Erecip + Ereal + Epoint                        **
// **                                                              **
// **      Epoint = [1/(4*pi*eo)]                                  **
// **                   * sqrt(eta)/sqrt(pi) * sum(i) qi^2         **
// **              +[1/(4*pi*e0)] * [pi/(2*vol*eta)]               **
// **                   * [sum(i) qi]^2   (for charged cell)       **
// **      Ereal = [1/(4*pi*e0)]* 0.5 *                            **
// **                      sum(ij) (qi*qj*erfc(sqrt(eta)*rij)/rij  **      
// **      Erecip = [1/(4*pi*e0)] *0.5* [4*pi/vol] *               **
// **               sum(G) 1/G^2 * exp[-G^2/(4*eta)] *             **
// **               |sum(i) (qi*exp(-iGri)|^2                      **
// **                                                              **
// **                                                              **
// **  For the forces:                                             **
// **      procedure: take the neg. derivatives.                   **
// **                 The point is const,                          **
// **                 and disappears.  So                          **
// **                                                              **
// **      Fj=  Fjrecip + Fjreal                                   **
// **                                                              **
// **      Fjxrecip = +[1/(4*pi*e0)] * [4*pi/vol] *                **
// **                 qj*  sum(G) Gx * G*1/G^2 *                   **
// **                  exp[-G^2/(4*eta)]*                          **
// **                             [-sin(Grj)*sum(i) qicos(Gri)]    **
// **                            +[cos(Grj)]*sum(i) qisin(Gri)]    **
// **       to get this, recall exp(-iGr)=cos(Gr)+isin(Gr)         **
// **       and put in Erecip.  Then take the derivatives of       **
// **       cos^2(Gr) and sin^2(Gr)                                **
// **                                                              **
// **      Fjxreal =   -[1/(4*pi*e0)]  * qj *                      **
// **                  Sum(i) qi/rij^3*                            **
// **                         [2*sqrt(eta)*rij/sqrt(pi)            **
// **                                  *exp(-eta*rij^2)            **
// **                          +erfc(sqrt(eta)*rij)]*rxij          **
// **                recall                                        **
// **        erf(x)=2/sqrt(pi)integral(0,x)exp(-t^2)dt             **
// **                                                              **
// **                                                              **
// **   input variables                                            **
// **                                                              **
// **     lat(3,3) real space lattice vectors input (ax ay az      **
// **                                              bx by bz        **
// **                                              cx cy cz)       **
// **                                               but stored     **
// **                                       ax bx cx...         .  **
// **                                       as per fort convention **
// **     rlat(3,3) recip space lattice vectors (see above)        **
// **     natoms  number of atoms                                  **
// **     ntype   # of types of atoms                              **
// **     vol     cell volume (A^3)                                **
// **     atfpos(natoms,3)   position of atoms in fractional       **
// **     atchg(natoms)      nominal charge of atoms by TYPE       **
// **     attyp(natoms)        type of atom (type of at 1, etc.    **
// **                                                              **
// **   output variables                                           **
// **                                                              **
// **     eewalde   ewald energy in   eV                           **
// **     ewaldf(3,natoms)   negative derivative of ewald          **
// **                        energy                                **
// **                          aka ewald forces ev/A               **
// ******************************************************************

namespace pflow {
  void Ewald(const xstructure& in_str,double& epoint,double& ereal,
	     double& erecip,double& eewald,double& eta,const double& SUMTOL) {
    xstructure str=in_str;
    str=PutInCell(str);
    matrix<double> lat = pflow::GetScaledLat(str);
    matrix<double> rlat = pflow::RecipLat(lat);
    int natoms=pflow::GetNumAtoms(str);
    double vol=pflow::GetVol(lat);
    matrix<double> atfpos=pflow::GetFpos(str);
    vector<int> attyp=pflow::GetTypes(str);
    vector<string> names=pflow::GetNames(str);
    // Charges should be stored in names - turn into vector of doubles.
    vector<double> atchg(natoms);
    double totchg=0;
    for(int i=0;i<(int) names.size();i++) {
      atchg[i]=atof(names[i].c_str());
      totchg=totchg+atchg[i];
    }

    // Constants
    double TOL=1e-6;

    // Check chg neutral
    if(abs(totchg)>TOL) {
      cerr << "WARNING: Cell is not neutral" << endl;
      cerr << "WARNING: Total charge = " << totchg << endl;
      cerr << "WARNING: Convasp will use jellium background to force neutrality."<< endl;
    }

    // Get eta if it was not input (it is <0 if it needs to be found).
    if(eta<=0) {eta=GetEta(natoms,vol);}
    // Constants
    double rteta=sqrt(eta);
    // Point energy
    epoint=GetPointEner(rteta,atchg,vol);
    // Recipricol energy
    erecip=GetRecipEner(eta,atchg,vol,rlat,atfpos,SUMTOL);
    // Real space energy
    ereal=GetRealEner(eta,atchg,vol,lat,atfpos,SUMTOL);
    // Sum all terms
    eewald=(ereal+erecip-epoint);
  }
}

// ***************************************************************************
// GetEta
// ***************************************************************************
namespace pflow {
  double GetEta(const int& natoms,const double& vol) {
    double term=std::pow((double) (5.5*natoms/vol/vol),(double) 1.0/6.0);
    double eta=sqrt(PI)/1.2*term;
    eta=eta*eta;
    return eta;
  }
}

// ***************************************************************************
// GetPointEner
// ***************************************************************************
// Gets point energy.
namespace pflow {
  double GetPointEner(const double& rteta, const vector<double>& atchg,
		      const double& vol) {
    double ept=0;
    double chg_cell_corr=0;
    for(int ia=0;ia<(int) atchg.size();ia++) {
      ept=ept+atchg[ia]*atchg[ia];
      chg_cell_corr=chg_cell_corr+atchg[ia];
    }
    ept=ept*rteta/RTPI;
    chg_cell_corr=chg_cell_corr*PI/(2.0*vol*rteta*rteta);
    ept=ept+chg_cell_corr;
    ept=CONV_FACT*ept;
    return ept;
  }
}

// ***************************************************************************
// GetRecipEner
// ***************************************************************************
// Gets recipricol energy.
namespace pflow {
  double GetRecipEner(const double& eta,const vector<double>& atchg,
		      const double& vol,const matrix<double>& rlat,
		      const matrix<double>& atfpos,const double& SUMTOL) {
    double TOL=1e-16;
    double log_eps=-30; // This assumes sf<
    double arg=0;
    double erecip=0;
    int gcnt=1;
    int gcont=1; // whether or not to add more shells.
    while(gcont) {
      double maxterm=0;
      vector<int> ig(3);
      for(ig[0]=-gcnt;ig[0]<=gcnt;ig[0]++) {
	for(ig[1]=-gcnt;ig[1]<=gcnt;ig[1]++) {
	  for(ig[2]=-gcnt;ig[2]<=gcnt;ig[2]++) {
	    if( (abs(ig[0])==gcnt) || (abs(ig[1])==gcnt) ||
		 (abs(ig[2])==gcnt) || gcnt==1) { // Just does new shells
	      vector<double> gvec=pflow::VMmult(ig,rlat);
	      double gsq=pflow::VVdot(gvec,gvec);
	      if(gsq>TOL) { // Avoid gsq=0.
		// get structure factor
		double sfr=0;
		double sfi=0;
		for(int ia=0;ia<(int) atfpos.size();ia++) {
		  double exparg=TWOPI*pflow::VVprod(atfpos[ia],ig);
		  sfr=sfr+atchg[ia]*cos(exparg);
		  sfi=sfi+atchg[ia]*sin(exparg);
		}
		if(fabs(sfr)<1e-16) {sfr=0.0;}
		if(fabs(sfi)<1e-16) {sfi=0.0;}
		double sf=sfr*sfr+sfi*sfi;
		// Get expval.  In order not to get floating point
		// errors due to manipulating small numbers we must
		// not evaluate the exp when it produces irrelevantly
		// small arguments.  This happens when
		// expval*sf << eps (where safe eps=1e-30), or equivalently,
		// when ln(sf)-ln(g^2)-g^2/(4eta) < ln(eps)
		arg=gsq/(4.0*eta);
		double expval;
		double term1=0;
		if(sf>0) {
		  term1=log(sf)-log(gsq)-gsq/(4*eta);
		}
		if(term1<log_eps) { // Avoids some floating point exceptions.
		  expval=0;
		}
		else{
		  expval=exp(-arg)/gsq;
		}
		double term=expval*sf;
		//	      if(term<1e-16) {term=0;}
		erecip=erecip+term;
		// Set maxterm to max of abs(term) and expval.
		// The expval term is needed since in some shells
		// term may be very small due to coincidently small sf,
		// even before convergence.
		if(fabs(term)>maxterm) {maxterm=fabs(term);}
		if(expval>maxterm) {maxterm=expval;}
		// ADD FORCES ???
	      } // If gsq>TOL
	    } // If to do only new shells (ig==gcnt)
	  } // ig0
	} // ig1
      } // ig2
      gcnt++; // Add new shell
      if(maxterm<SUMTOL) {gcont=0;} // Do not add another shell
    } // while gcont
    erecip=erecip*0.5*4.0*PI/vol; // CGS units, convert to eV at end.
    erecip=CONV_FACT*erecip;
    // FORCES STUFF HERE ???				  
    return erecip;				  
  } // End GetRecipEner
}

// ***************************************************************************
// GetRealEner
// ***************************************************************************
// Gets real space energy.
namespace pflow {
  double GetRealEner(const double& eta,const vector<double>& atchg,
		     const double& vol,const matrix<double>& lat,
		     const matrix<double>& atfpos,const double& SUMTOL) {
    if(vol) {;} // phony just to keep vol busy

    double ereal=0;
    double TOL=1e-5;
    double log_eps=-30;
    double erfcarg=0;
    double rteta=sqrt(eta);
    int nat=atfpos.size();
    int rcnt=1;
    int rcont=1;
    while(rcont) {
      double maxterm=0;
      vector<int> ir(3);
      for(ir[0]=-rcnt;ir[0]<=rcnt;ir[0]++) {
	for(ir[1]=-rcnt;ir[1]<=rcnt;ir[1]++) {
	  for(ir[2]=-rcnt;ir[2]<=rcnt;ir[2]++) {
	    if( (abs(ir[0])==rcnt) || (abs(ir[1])==rcnt) ||
		 (abs(ir[2])==rcnt) || rcnt==1) { // Just does new shells
	      for(int ia=0;ia<nat;ia++) { // All atoms in unit cell
		for(int ja=0;ja<nat;ja++) { // All neighbors in cell given by ir
		  // Get displacement between atoms ia(cell 0) and ja(cell ir).
		  vector<double> disp=pflow::VVdiff(atfpos[ja],atfpos[ia]);
		  disp=pflow::VVsum(disp,ir);
		  disp=pflow::vecF2C(lat,disp);
		  double dist=pflow::norm(disp);
		  if(dist>TOL) { // Avoid atom dist to itself.
		    erfcarg=rteta*dist;
		    // Get erfcval.  In order not to get floating point
		    // errors due to manipulating small numbers we must
		    // not evaluate the erfc when it produces irrelevantly
		    // small arguments.  This happens when
		    // erfc(x)*abs(q) << eps (where safe eps=1e-30), or equivalently,
		    // when ln(abs(q))-x^2-ln(x)-0.5*ln(pi) < ln(eps) (where I have
		    // used the asymptotic form erfc(x)~e^-x^2/(x*rt(pi))
		    // for large x.
		    double erfcval;
		    double term1=atchg[ia]*atchg[ja]/dist;
		    double term2=0;
		    if(abs(term1)>0) {
		      term2=log(abs(term1))-erfcarg*erfcarg-log(erfcarg)-0.5*log(PI);
		    }
		    if(term2<log_eps) { // Avoids some floating point exceptions.
		      erfcval=0;
		    }
		    else{
		      erfcval=erfc(erfcarg);
		    }
		    // tpx
		    // cout << "term1 erfcval " << term1 << " " << erfcval << endl;
		    double term=term1*erfcval;
		    ereal=ereal+term;
		    if(fabs(term)>maxterm) {maxterm=fabs(term);}
		    // FORCES ???
		  } // rdist>TOL
		} // ja
	      } // ia
	    } // if to do only new shells (ir==rcnt)
	  } // rg0
	} // rg1
      } // rg2
      rcnt++; // Add new shell
      if(maxterm<SUMTOL) {rcont=0;} // Do not add another shell
    } // while rcont
    ereal=ereal*0.5; // CGS units, convert to eV at end.
    ereal=CONV_FACT*ereal;
    return ereal;				  
  }  
}

// ***************************************************************************
// GetSreenedESEner
// ***************************************************************************
// Gets electrostatic energy with screening (does sum in real space).
namespace pflow {
  double ScreenedESEner(const xstructure& in_str,const double& Ks,const double& SUMTOL) {
    xstructure str=in_str;
    str=PutInCell(str);
    matrix<double> lat = pflow::GetScaledLat(str);
    int natoms=pflow::GetNumAtoms(str);
    // double vol=pflow::GetVol(lat);  // DANE not used
    matrix<double> atfpos=pflow::GetFpos(str);
    vector<int> attyp=pflow::GetTypes(str);
    vector<string> names=pflow::GetNames(str);
    // Charges should be stored in names - turn into vector of doubles.
    vector<double> atchg(natoms);
    double totchg=0;
    for(int i=0;i<(int)names.size();i++) {
      atchg[i]=atof(names[i].c_str());
      totchg=totchg+atchg[i];
    }
    double ereal=0;
    double TOL=1e-5;
    int nat=atfpos.size();
    int rcnt=1;
    int rcont=1;
    while(rcont) {
      double maxterm=0;
      vector<int> ir(3);
      for(ir[0]=-rcnt;ir[0]<=rcnt;ir[0]++) {
	for(ir[1]=-rcnt;ir[1]<=rcnt;ir[1]++) {
	  for(ir[2]=-rcnt;ir[2]<=rcnt;ir[2]++) {
	    if( (abs(ir[0])==rcnt) || (abs(ir[1])==rcnt) ||
		 (abs(ir[2])==rcnt) || rcnt==1) { // Just does new shells
	      for(int ia=0;ia<nat;ia++) { // All atoms in unit cell
		for(int ja=0;ja<nat;ja++) { // All neighbors in cell given by ir
		  // Get displacement between atoms ia(cell 0) and ja(cell ir).
		  vector<double> disp=pflow::VVdiff(atfpos[ja],atfpos[ia]);
		  disp=pflow::VVsum(disp,ir);
		  disp=pflow::vecF2C(lat,disp);
		  double dist=pflow::norm(disp);
		  if(dist>TOL) { // Avoid atom dist to itself.
		    double term=exp(-Ks*dist)*atchg[ia]*atchg[ja]/dist;
		    ereal=ereal+term;
		    if(fabs(term)>maxterm) {maxterm=fabs(term);}
		    // FORCES ???
		  } // rdist>TOL
		} // ja
	      } // ia
	    } // if to do only new shells (ir==rcnt)
	  } // rg0
	} // rg1
      } // rg2
      rcnt++; // Add new shell
      if(maxterm<SUMTOL) {rcont=0;} // Do not add another shell
    } // while rcont
    ereal=ereal*0.5; // CGS units, convert to eV at end.
    ereal=CONV_FACT*ereal;
    return ereal;				  
  }  
}

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
#endif
