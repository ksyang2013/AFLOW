// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
//
// Acknowledgements to Dane Morgan and Anton Van der Ven for discussions and
// suggestions.

#include "aflow.h"

#define cdebug cerr
#define _EPS_ 0.001
#define _EPS_roundoff_ 1.0e-8

//  int myints[] = {16,2,77,29};
//  vector<int> fifth (myints, myints + sizeof(myints) / sizeof(int));

template<class utype>
vector<utype> svector2vector(utype svec[]) {
  vector<utype> vec(svec, svec + sizeof(svec) / sizeof(utype));
  return vec;
}

// --------------------------------------------------------------------------
// ----------------------------------------------------- NEIGHBOURS OPERATION
bool StepNeighboursPerform(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
  xstructure str;
  str=BringInCell(a);
  ostringstream aus;
  bool Krun=TRUE;
  bool LVERBOSE=TRUE;
  // DO THE NEIGHBOURS CALCULATION
  if(aurostd::substring2bool(AflowIn,"[AFLOW_NEIGHBOURS]CALC",TRUE)) {
    kflags.KBIN_NEIGHBOURS_CALCULATION=TRUE;
    kflags.KBIN_NEIGHBOURS_WRITE=aurostd::substring2bool(AflowIn,"[AFLOW_NEIGHBOURS]WRITE",TRUE);
    kflags.KBIN_NEIGHBOURS_RADIUS=aurostd::substring2utype<double>(AflowIn,"[AFLOW_NEIGHBOURS]]RADIUS=",TRUE);
    kflags.KBIN_NEIGHBOURS_DRADIUS=aurostd::substring2utype<double>(AflowIn,"[AFLOW_NEIGHBOURS]DRADIUS=",TRUE);
  }
  if(kflags.KBIN_NEIGHBOURS_CALCULATION==FALSE) return TRUE;
  
  // CALC = TRUE then prepare and act
  if(kflags.KBIN_NEIGHBOURS_RADIUS>0.0) {
    if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS found RADIUS="
		     << kflags.KBIN_NEIGHBOURS_RADIUS<<" " << Message(aflags,"user,host,time") << endl;
    if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  } else {
    kflags.KBIN_NEIGHBOURS_RADIUS=KBIN_NEIGHBOURS_RADIUS_DEFAULT;
    if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS default RADIUS="
		     << kflags.KBIN_NEIGHBOURS_RADIUS<<" " << Message(aflags,"user,host,time") << endl;
    if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  }
  if(kflags.KBIN_NEIGHBOURS_DRADIUS>0.0) {
    if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS found DRADIUS="
		     << kflags.KBIN_NEIGHBOURS_DRADIUS<<" " << Message(aflags,"user,host,time") << endl;
    if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
    kflags.KBIN_NEIGHBOURS_DRADIUS=KBIN_NEIGHBOURS_DRADIUS_DEFAULT;
    if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS default DRADIUS="
		     << kflags.KBIN_NEIGHBOURS_DRADIUS<<" " << Message(aflags,"user,host,time") << endl;
    if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  }
  
  str.neighbours_radius=kflags.KBIN_NEIGHBOURS_RADIUS;
  str.neighbours_dradius=kflags.KBIN_NEIGHBOURS_DRADIUS;
  
  xmatrix<double> lattice(3,3); lattice=str.lattice;
  xvector<double> a1(3),a2(3),a3(3);
  a1=str.lattice(1);a2=str.lattice(2);a3=str.lattice(3);
  // bool atom_found;
  
  /*
  str.nbins=(int) ceil(str.neighbours_radius/str.neighbours_dradius)+1;
  uint numatoms=str.atoms.size();
  int ibin;
  str.ndims=LatticeDimensionSphere(lattice,str.neighbours_radius);
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS inside sphere, dimensions = ["
		   << str.ndims[1] << "," << str.ndims[2] << "," << str.ndims[3] << "] "
		   << Message(aflags,"user,host,time") << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  
  _atom atom;
  uint iat,iat1,iat2;

  str.natoms=deque<deque<_atom> >(numatoms);
  str.rshell=vector<vector<double> >(numatoms);
  str.nshell=vector<vector<int> >(numatoms);
  for(iat=0;iat<numatoms;iat++) {
    //  str.natoms.push_back(deque<_atom>(10));
    str.rshell.at(iat)=vector<double>(str.nbins,0.0);
    str.nshell.at(iat)=vector<int>(str.nbins,0);
  }
  // --------------------------- preparing list
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS preparing list " << Message(aflags,"user,host,time") << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  for(iat1=0;iat1<numatoms;iat1++) {
    for(int i=-str.ndims[1];i<=str.ndims[1];i++)
      for(int j=-str.ndims[2];j<=str.ndims[2];j++)
	for(int k=-str.ndims[3];k<=str.ndims[3];k++)
	  for(iat2=0;iat2<numatoms;iat2++) {
	    atom=str.atoms.at(iat2);
	    // atom.cpos=atom.cpos+i*a1+j*a2+k*a3;
	    // atom.fpos(1)+=i;atom.fpos(2)+=j;atom.fpos(3)+=k;
	    atom.ijk(1)=i;atom.ijk(2)=j;atom.ijk(3)=k;
	    // atom.reference=str.scale*AtomDist(str.atoms.at(iat1),atom);
	    atom.reference=AtomDist(str,str.atoms.at(iat1),atom);  // without the ijk
	    if(atom.reference<=str.neighbours_radius) {
	      ibin=(int) ((double) atom.reference/str.neighbours_dradius);  // get index
	      // cerr << " atom.reference="  << atom.reference << " str.neighbours_dradius="  << str.neighbours_dradius << " ibin=" << ibin << " str.nbins=" << str.nbins << endl;
	      str.nshell.at(iat1).at(ibin)+=1;                        // add to the bin
	      str.rshell.at(iat1).at(ibin)+=atom.reference;           // make the shell (must be divided)
	      atom.print_cartesian=FALSE;
	      atom.verbose=FALSE;
	      str.natoms.at(iat1).push_back(atom);
	      // cout << atom << endl;
	    }
	  }
  }
  //  --------------------------- sort the natoms by distance
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS sorting by distances " << Message(aflags,"user,host,time") << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  for(iat1=0;iat1<numatoms;iat1++)                                                      // sort the natoms by distance
    sort(str.natoms.at(iat1).begin(),str.natoms.at(iat1).end(),_atom_reference_cmp());  // sort the natoms by distance
  
  //  --------------------------- normalize shell
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS normalize shell " << Message(aflags,"user,host,time") << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS number: ";
  for(iat1=0;iat1<numatoms;iat1++) {                                                    // normalize shell
    if(LVERBOSE) aus << str.natoms.at(iat1).size() << " ";                              // normalize shell
    for(ibin=0;ibin<str.nbins;ibin++)                                                   // normalize shell
      if(str.nshell.at(iat1).at(ibin)>0)                                                // normalize shell
	str.rshell.at(iat1).at(ibin)/=str.nshell.at(iat1).at(ibin);                     // normalize shell
  }
  if(LVERBOSE) aus << Message(aflags,"user,host,time") << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  
  //  --------------------------- create unique list
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS create unique list of atoms " << Message(aflags,"user,host,time") << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  if(LVERBOSE) cout << "00000  MESSAGE NEIGHBOURS working size(before)=";
  
  for(iat1=0;iat1<numatoms;iat1++)                                                      // create unique list
    for(iat2=0;iat2<str.natoms.at(iat1).size();iat2++)                                  // create unique list
      str.ashell.push_back(str.natoms.at(iat1).at(iat2));                               // create unique list
  
  atom=str.atoms.at(0); // just a reference.
  atom.cpos.clear();    // clear the atom, just a reference.
  atom.fpos.clear();    // clear the atom, just a reference.
  atom.ijk.clear();     // clear the atom, just a reference.
  atom.cpos=str.origin;
  for(iat=0;iat<str.ashell.size();iat++)
    str.ashell.at(iat).reference=AtomDist(str,str.ashell.at(iat),atom);
  if(LVERBOSE) {cout << str.ashell.size() << " "; cout.flush();} // << endl;
  sort(str.ashell.begin(),str.ashell.end(),_atom_reference_cmp());  // sort the natoms by distance
  
  int _max_dist=0;
  int _mod_dist=(int) str.natoms.at(0).size()/10; // str.ashell.size()/1000;
  for(iat1=0;iat1<numatoms;iat1++)                                                     // normalize shell
    for(ibin=0;ibin<str.nbins;ibin++)                                                  // normalize shell
      if(str.nshell.at(iat1).at(ibin)>=_max_dist)                                      // normalize shell
	_max_dist=str.nshell.at(iat1).at(ibin);                                        // normalize shell
  _max_dist+=10+numatoms;                                                              // seems to work
  
  //  for(i=0;i<3;i++)
  { // must be repeated  
    for(iat1=0;iat1<str.ashell.size();iat1++) {                                        // create unique list
      if(!mod((int) iat1,_mod_dist)) if(LVERBOSE) {cout << ".";cout.flush();}          // create unique list
      for(iat2=iat1+1;iat2<iat1+_max_dist+2 && iat2<str.ashell.size();iat2++)          // create unique list
	if(iat1!=iat2)                                                                 // create unique list
	  if(SameAtom(str,str.ashell.at(iat1),str.ashell.at(iat2)))                    // create unique list
	    str.ashell.erase(str.ashell.begin()+iat2--);                               // create unique list
    }
    sort(str.ashell.begin(),str.ashell.end(),_atom_reference_cmp());                   // create unique list
  }
  if(LVERBOSE) aus << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  sort(str.ashell.begin(),str.ashell.end(),_atom_reference_cmp());                     // sort natoms by distance
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS working size(after)=";
  if(LVERBOSE) aus << str.ashell.size() << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  if(LVERBOSE) cout << "00000  MESSAGE NEIGHBOURS working size(double checking)";
  
  for(iat1=0;iat1<str.ashell.size();iat1++) {                                          // be double sure
    if(!mod((int)iat1,(int)str.ashell.size()/25)) if(LVERBOSE) {cout<<".";cout.flush();}// be double sure
    for(iat2=0;iat2<str.ashell.size();iat2++)                                          // be double sure
      if(iat2!=iat1)                                                                   // be double sure
	if(SameAtom(str,str.ashell.at(iat1),str.ashell.at(iat2)))                      // be double sure
	  str.ashell.erase(str.ashell.begin()+iat2--);                                 // be double sure
  }                                                                                    // be double sure
  
  if(LVERBOSE) aus << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  sort(str.ashell.begin(),str.ashell.end(),_atom_reference_cmp());  // sort the natoms by distance
  if(LVERBOSE) aus << "00000  MESSAGE NEIGHBOURS working size(after)=";
  if(LVERBOSE) aus << str.ashell.size() << endl;
  if(LVERBOSE) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  
  if(str.ashell.size()>KBIN_NEIGHBOURS_MAX_NUMBER) {
    aus << "EEEEE  ERROR: StepNeighboursPerform cannot deal with more than " << KBIN_NEIGHBOURS_MAX_NUMBER << " atoms." << endl;
    aus << "EEEEE  ERROR: You need to reduce the radius of the neighbours calculation (RADIUS=XXX) or" << endl;
    aus << "EEEEE  ERROR: modify KBIN_NEIGHBOURS_MAX_NUMBER in aflow.h and recompile." << endl;
    aus << "EEEEE  ERROR: WARNING: by increasing KBIN_NEIGHBOURS_MAX_NUMBER you can use a lot of memory" << endl;
    aus << "EEEEE  ERROR: and disk space if you decide to save the space group file (WRITE)." << endl;
    aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
    return FALSE;
  }
  
  if(!str.pgroup_calculated) {
    aus << "EEEEE  ERROR: Symmetry must be calculated to perform neightbour symmetry list." << endl;
    aus << "EEEEE  ERROR: add \"[AFLOW_NEIGHBOURS]CALCULATION\" to your " << _AFLOWIN_ << " and restart." << endl;
    aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
    return FALSE;
  }
  
  {
    vector<vector<int> > ashell_pgroup(str.ashell.size(),str.pgroup.size());
    for(iat=0;iat<30;iat++) {                                 // for every atom
      //    for(iat=0;iat<str.ashell.size();iat++) {                                 // for every atom
      atom = str.ashell.at(iat);
      cerr << "X " << atom << endl;
      //      for(iop=0;iop<str.pgroup.size();iop++) {                               // for every operation
      for(uint iop=0;iop<1;iop++) {                               // for every operation
	//	atom=SYM_Apply(str.ashell.at(iat),str.pgroup.at(iop),str);
	cerr << iop << " " << atom << endl;
      }
      cerr << endl;
    }
  }
  
  exit(0);
  
  //    for(ii=0;ii<str.ashell.size();ii++) cout << str.ashell.at(ii) << " " << str.ashell.at(ii).reference << endl;

  //      for(iat1=0;iat1<str.natoms.size();iat1++)
  //        for(iat2=0;iat2<str.natoms.at(iat1).size();iat2++)
  // 	 cout << str.natoms.at(iat1).at(iat2) << endl;

  //     for(iat1=0;iat1<numatoms;iat1++) {
  //       for(ibin=0;ibin<str.nbins;ibin++)
  //     	cout << str.nshell.at(iat1).at(ibin) << " ";
  //       cout << endl;
  //     }

  for(iat1=0;iat1<numatoms;iat1++) {
    cout << str.natoms.at(iat1).size() << endl;
    //  cout << str.atoms.at(iat1).type << endl;
  }
  cout << setprecision(10);
  for(iat1=0;iat1<numatoms;iat1++) {
    //       // cout << str.natoms.at(iat1).size() << " " ;for(iat2=0;iat2<10;iat2++) {cout << str.natoms.at(iat1).at(iat2).reference << " ";} cout << endl;
    // for(iat2=0;iat2<10;iat2++) {cout << str.natoms.at(iat1).at(iat2).number+1 << "/" << str.natoms.at(iat1).at(iat2).reference << " ";} cout << endl;
  }

  str.write_DEBUG_flag=TRUE;
  cout << a << endl;

  exit(0);
  str.neighbours_calculated=TRUE;
  // Krun=(Krun && SYM::CalculateSpaceGroup(FileMESSAGE,a,aflags,kflags.KBIN_NEIGHBOURS_WRITE));
  //  a=str;




  */

  return Krun;
}

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
