// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *           Aflow COREY OSES - Duke University 2013-2017                  *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XCOMBOS_CPP_
#define _AUROSTD_XCOMBOS_CPP_

#include "aurostd_xcombos.h"

namespace aurostd {
//--------------------------------------------------------------------------------
// class xcombos
//--------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// constructor
//------------------------------------------------------------------------------
xcombos::xcombos() {free();}
xcombos::xcombos(const std::vector<int>& vec,bool sort) {reset(vec,sort);}
xcombos::xcombos(int choice_count,int choose_count) {reset(choice_count,choose_count);}
xcombos::xcombos(const xcombos& b) {copy(b);} // copy PUBLIC
xcombos::~xcombos() {free();}

void xcombos::free() {
  m_initialized=FALSE;
  m_input.clear();
  n_choices=0;
  m_choose=0;
  m_permutation=FALSE;
  m_sort=FALSE;
  m_started=FALSE;
  m_exhausted=FALSE;
  m_current.clear();
  m_p.clear();
  m_x=0;
  m_y=0;
}

void xcombos::copy(const xcombos& b) {  //copy PRIVATE
  m_initialized=b.m_initialized;
  m_input=b.m_input;
  n_choices=b.n_choices;
  m_choose=b.m_choose;
  m_permutation=b.m_permutation;
  m_sort=b.m_sort;
  m_started=b.m_started;
  m_exhausted=b.m_exhausted;
  m_current=b.m_current;
  m_p=b.m_p;
  m_x=b.m_x;
  m_y=b.m_y;
}

const xcombos& xcombos::operator=(const xcombos& other) {
  if(this!=&other) {free();copy(other);}
  return *this;
}

xcombos& xcombos::operator++() {  //remember, this is PREFIX (faster than POSTFIX)
  if(!m_initialized) {
    cerr << "xcombos::operator++(): ERROR - cannot increment uninitializd xcombos class" << endl;
    exit(1);
  }
  if(m_permutation) {incrementPermutation();}
  else {incrementCombinations();}
  return *this;
}

void xcombos::reset() {  //clear PUBLIC
  if(m_permutation) {reset(m_input,m_sort);}
  else {reset(n_choices,m_choose);}
}

void xcombos::reset(std::vector<int> vec,bool sort) { //do NOT make input vec a const &, this will screw up reset()
  free();
  m_input=vec;
  m_permutation=TRUE;
  m_sort=sort;
  m_started=FALSE;
  if(!m_input.size()){m_exhausted=TRUE; return;}  //safe
  initialize();
}

void xcombos::reset(int choice_count,int choose_count) {
  free();
  n_choices=choice_count;
  m_choose=choose_count;
  m_permutation=FALSE;
  m_started=FALSE;
  if(n_choices<=0 || m_choose<=0 || m_choose>n_choices) {m_exhausted=TRUE; return;}  //safe
  initialize();
}

void xcombos::initialize() {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(m_permutation) {
    m_current=m_input;
    if(m_sort) {std::sort(m_current.begin(),m_current.end());}
  }else {
    if(LDEBUG) {cerr << "combinations: " << n_choices << " choose " << m_choose << endl;}
    m_current.resize(n_choices);
    for(uint i=0;i<m_current.size();i++) {
      if((int)i<m_choose) {m_current[i]=1;}
      else {m_current[i]=0;}
    }
		initializeCombinationsP();
  }
  m_initialized=TRUE;
}

std::vector<int> xcombos::getCombo() const {return m_current;}
int xcombos::getN() const {return n_choices;}
int xcombos::getM() const {return m_choose;}

std::vector<uint> xcombos::getIndices() const {
  std::vector<uint> vout;
  for(uint i=0;i<m_current.size();i++) {
    if(m_current[i]==1) {vout.push_back(i);}
  }
  return vout;
}

bool xcombos::increment() {
  if(!m_initialized) {return FALSE;}     //safety
  if(m_exhausted==TRUE) {return FALSE;}  //safety
  if(!m_started) {m_started=TRUE;return TRUE;} //allows us to use increment() for first grab too!
  ++(*this);
  if(m_exhausted==TRUE) {return FALSE;}
  return TRUE;
}

void xcombos::incrementPermutation() {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  if(m_exhausted) {return;}
  //Shen, MK. BIT (1962) 2(228). doi:10.1007/BF01940170
  //note this will generate next permutation in lexicographical order (left to right)
  int _i=-1;
  int _j=-1;
  for(int i=1;i<(int)m_current.size();i++) {if(m_current[i-1]<m_current[i]&&(i>_i)){_i=i;}}
  if(_i==-1) {m_exhausted=TRUE; return;} //stop condition
  for(int j=0;j<(int)m_current.size();j++) {if(m_current[_i-1]<m_current[j]&&(j>_j)){_j=j;}}
  if(LDEBUG) {cerr << "xcombos::incrementPermutation(): i=" << _i << "  j=" << _j << " " << endl;}
  std::swap(m_current[_i-1],m_current[_j]);
  for(int i=0;i<((int)m_current.size()-_i+1)/2;i++) {std::swap(m_current[_i+i],m_current[m_current.size()-i-1]);}
}

void xcombos::incrementCombinations() {
  if(m_exhausted) {return;}
	setCombinationsIncrementParameters();
  m_current[m_x]=1;
  m_current[m_y]=0;
  //for(uint i=0;i<m_current.size();i++) {cerr << m_current[i];}
  //cerr << endl;
}

void xcombos::initializeCombinationsP() {  //combinations only
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  int i=0;
  m_p.resize(n_choices+2);
  m_p[i++]=-2;
  while(m_choose-i+1>0) {m_p[i]=m_choose-i+1;i++;}
  while(i<n_choices+1) {m_p[i++]=0;}
  m_p[n_choices+1]=n_choices+1;
  if(m_choose==0) {m_p[n_choices]=1;}
  if(LDEBUG) {
    cerr << "xcombos::initializeCombinationsP(): p(0): ";
    for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
    cerr << endl;
  }
}

void xcombos::setCombinationsIncrementParameters() { //combinations only
  if(!m_initialized) {return;}
  if(m_exhausted) {return;}
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //Chase, PJ. ACM (1970) 13(6). doi:10.1145/362384.362502
  //implementation modified from twiddle.c by Matthew Belmonte
  //http://www.netlib.no/netlib/toms/382
  //note this will generate next permutation in lexicographical order (left to right)
  //originally, Belmonte programmed right to left
	int i, j, k;
  i=j=k=-1;
	j=n_choices;
  if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): j(0)=" << j << endl;}
	while(m_p[j] <= 0) {j--;}
  if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): j(1)=" << j << endl;}
	if(m_p[j+1]==0) {
		for(i=j+1; i != n_choices; i++) {m_p[i]=-1;}
		m_p[j]=0;
		m_x=n_choices-1;
		m_p[n_choices]=1;
		m_y=j-1;
    if(LDEBUG) {
      cerr << "xcombos::setCombinationsIncrementParameters(): p(0): ";
      for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
      cerr << endl;
      cerr << "xcombos::setCombinationsIncrementParameters(): x=" << m_x << "(1)" << endl;
      cerr << "xcombos::setCombinationsIncrementParameters(): y=" << m_y << "(0)" << endl;
    }
	}else {
		if(j<n_choices) {m_p[j+1]=0;}
		j--;
		while(m_p[j]>0) {j--;}
    if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): j(2)=" << j << endl;}
		k=j+1;
		i=j;
    if(LDEBUG) {
      cerr << "xcombos::setCombinationsIncrementParameters(): i(0)=" << i << endl;
      cerr << "xcombos::setCombinationsIncrementParameters(): k(0)=" << k << endl;
    }
		while(m_p[i]==0) {m_p[i--]=-1;}
    if(LDEBUG) {cerr << "xcombos::setCombinationsIncrementParameters(): i(1)=" << i << endl;}
		if(m_p[i]==-1) {
			m_p[i]=m_p[k];
			m_x=i-1;
			m_y=k-1;
			m_p[k]=-1;
      if(LDEBUG) {
        cerr << "xcombos::setCombinationsIncrementParameters(): p(1): ";
        for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
        cerr << endl;
        cerr << "x=" << m_x << "(1)" << endl;
        cerr << "y=" << m_y << "(0)" << endl;
      }
		}else {
      if(LDEBUG) {
        cerr << "xcombos::setCombinationsIncrementParameters(): i=" << i << endl;
        cerr << "xcombos::setCombinationsIncrementParameters(): p(2): ";
        for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
        cerr << endl;
      }
			if(i==0) {m_exhausted=TRUE; return;}
			//if(i==m_p[n_choices+1]) {m_exhausted=TRUE; cerr << "HERE" << endl;return;}
			else {
				m_p[j]=m_p[i];
				m_p[i]=0;
				m_x=j-1;
				m_y=i-1;
        if(LDEBUG) {
          cerr << "xcombos::setCombinationsIncrementParameters(): p(3): ";
          for(uint l=0;l<m_p.size();l++) {cerr << m_p[l] << " ";} 
          cerr << endl;
          cerr << "x=" << m_x << "(1)" << endl;
          cerr << "y=" << m_y << "(0)" << endl;
          //exit(0);
        }
			}
		}
	}
  return;
}

} // namespace aurostd

#endif  // _AUROSTD_XCOMBOS_CPP_
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *           Aflow COREY OSES - Duke University 2013-2017                  *
// *                                                                         *
// ***************************************************************************
