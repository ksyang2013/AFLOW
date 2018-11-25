// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 2013-2014
// added template<class utype> bool xoption::args2addattachedscheme SC 2017
// streamline schemes SC 2017

#ifndef _AUROSTD_XOPTION_CPP_
#define _AUROSTD_XOPTION_CPP_
//#include "aflow.h"
using std::ostream;

namespace aurostd {
  // ***************************************************************************

  // constructure
  xoption::xoption() {
    keyword=""; //DX 20180824 - missing from constructor
    isentry=FALSE; 
    content_string="";
    content_double=0.0;
    content_int=0;
    content_uint=0;
    option=FALSE;
    option_default=FALSE;
    xscheme="";
    vxscheme.clear();
    vxsghost.clear();
    preserved=FALSE;
    LDEBUG=FALSE;
  }

  // destructor
  xoption::~xoption() {
    free();
  }

  // free
  void xoption::free() {
  }

  // copy fuction
  void xoption::copy(const xoption& b) {
    keyword=b.keyword; //DX 20180824 - missing from copy constructor
    isentry=b.isentry;
    content_string=b.content_string;
    content_double=b.content_double;
    content_int=b.content_int;
    content_uint=b.content_uint;
    option=b.option;
    option_default=b.option_default;
    xscheme=b.xscheme;
    vxscheme.clear();for(uint i=0;i<b.vxscheme.size();i++) vxscheme.push_back(b.vxscheme.at(i));
    vxsghost.clear();for(uint i=0;i<b.vxsghost.size();i++) vxsghost.push_back(b.vxsghost.at(i));
    preserved=b.preserved;
    LDEBUG=b.LDEBUG;
  }

  // copy conctructor
  xoption::xoption(const xoption& b) {
    //  free();
    // *this=b;
    copy(b);
  }

  // copy operator b=a
  const xoption& xoption::operator=(const xoption& b) {  // operator=
    if(this!=&b) {
      free();
      copy(b);}
    return *this;
  }
  
  std::ostream& operator<<(std::ostream& oss,const xoption& a) {
    for(uint i=0;i<a.vxscheme.size();i++)
      oss << a.vxscheme.at(i) << (i<a.vxscheme.size()-1?",":"");
    return oss;
  }
  
  void xoption::clear() {
    xoption aflow_option_temp;
    copy(aflow_option_temp);
  }

  // **************************************************************************
  void xoption::options2entry(string options_FILE,string input_keyword,int _option_DEFAULT,string xscheme_DEFAULT) {
    clear();
    LDEBUG=FALSE;
    bool option_DEFAULT=FALSE; 
    if(_option_DEFAULT==0) option_DEFAULT=FALSE; // it is a int.. it might be -1
    if(_option_DEFAULT==1) option_DEFAULT=TRUE; // it is a int.. it might be -1
    isentry=option_DEFAULT;option=option_DEFAULT;content_string=xscheme_DEFAULT;xscheme=xscheme_DEFAULT;preserved=FALSE;   // DEFAULT
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: BEGIN " << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: input_keyword=\"" << input_keyword << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: option_DEFAULT=" << (option_DEFAULT?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: xscheme_DEFAULT=\"" << xscheme_DEFAULT << "\"" << endl;
    // start the scan
    //string keyword; //CO 180404 - now a member of the object
    vector<string> vkeyword;
    // tokenize the option
    aurostd::string2tokens(input_keyword,vkeyword,"|"); 
    if(LDEBUG) for(uint i=0;i<vkeyword.size();i++) cerr << "\"" << vkeyword.at(i) << "\"" << endl;
    // loop through the scan
    if(vkeyword.size()>0) {
      // some default
      keyword=vkeyword.at(0);
      for(uint i=0;i<vkeyword.size();i++) 
	if(aurostd::substring2bool(options_FILE,vkeyword.at(i),TRUE))
	  keyword=vkeyword.at(i);
      // found one keyword
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: keyword=\"" << keyword << "\"" << endl;
      // LOOK FOR EXIST/!EXIST ENTRY
      if(_option_DEFAULT==aurostd_xoptionONOFF) {
	isentry=aurostd::substring2bool(options_FILE,keyword,TRUE);
	if(isentry) {option=TRUE;content_string="ON";}
	if(!isentry) {option=FALSE;content_string="OFF";}
      } // aurostd_xoptionONOFF exit/~exit
      // LOOK FOR ON/OFF MODE WITH strings/schemes.
      if(_option_DEFAULT==0 || _option_DEFAULT==1) {
	if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: LOOK FOR ON/OFF MODE WITH strings/schemes" << endl;
	// start the scan
	isentry=aurostd::substring2bool(options_FILE,keyword,TRUE);
	if(isentry && xscheme_DEFAULT.empty()) {
	  content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,FALSE));
	  string saus=content_string;content_string="";
	  if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: found saus=" << saus << endl;
	  vector<string> tokens;aurostd::string2tokens(saus,tokens,",");
	  for(uint i=0;i<tokens.size();i++) { //      c<< tokens.at(i) << endl;
	    if(tokens.at(i)=="ON" || tokens.at(i)[0]=='T' || tokens.at(i)[0]=='t' || tokens.at(i)[0]=='1' || tokens.at(i)[0]=='Y' || tokens.at(i)[0]=='y') {
	      option=TRUE;content_string=saus;} // modify option and value
	    if(tokens.at(i)=="OFF"|| tokens.at(i)[0]=='F' || tokens.at(i)[0]=='f' || tokens.at(i)[0]=='0' || tokens.at(i)[0]=='N' || tokens.at(i)[0]=='n') {
	      option=FALSE;content_string=saus;} // modify option and value
	    // if(tokens.at(i)=="REMOVE_RELAX_1") {option=TRUE;content_string=saus;} // modify option and value // compatibility with SPIN
	    // if(tokens.at(i)=="REMOVE_RELAX_2") {option=TRUE;content_string=saus;} // modify option and value // compatibility with SPIN
	    if(tokens.at(i)=="REMOVE_RELAX_1") {content_string=saus;} // modify option and value // compatibility with SPIN but dont touch ON/OFF
	    if(tokens.at(i)=="REMOVE_RELAX_2") {content_string=saus;} // modify option and value // compatibility with SPIN but dont touch ON/OFF
	  }
	}
	// SCHEME MODE
	if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: xscheme_DEFAULT=\"" << xscheme_DEFAULT << "\"" << endl;
	if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: xscheme_DEFAULT.empty()=" << xscheme_DEFAULT.empty() << endl;
	if(isentry && !xscheme_DEFAULT.empty()) {
	  if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: SCHEME MODE" << endl;
	  content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,FALSE));
	  option=isentry;
	}
	if(isentry && (xscheme_DEFAULT.empty() && content_string.empty())) {
	  if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: SCHEME MODE EMPTY DEFAULT STILL EMPTY CONTENT" << endl;
	  content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,FALSE));
	  option=isentry;
	}
	if(!isentry && option_DEFAULT) {
	  option=TRUE;content_string="ON";
	}
      } // 0/1 on/off mode
      // LOOK FOR EXIST/!EXIST ENTRY
      if(_option_DEFAULT==aurostd_xoptionMULTI) {
	vector<string> voptions_FILE,vcontent;
	aurostd::string2vectorstring(options_FILE,voptions_FILE);
	isentry=FALSE;content_string="";
	for(uint i=0;i<voptions_FILE.size();i++) {
	  if(aurostd::substring2bool(voptions_FILE.at(i),keyword,TRUE)) {
	    vector<string> vstrcheck;
	    string strcheck=aurostd::toupper(aurostd::RemoveWhiteSpaces(aurostd::substring2string(voptions_FILE.at(i),keyword,FALSE)));
	    aurostd::StringSubst(strcheck,";",",");
	    aurostd::string2tokens(strcheck,vstrcheck,","); 
	    for(uint j=0;j<vstrcheck.size();j++) {
	      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: BEFORE keyword=" << keyword << "   " << "vstrcheck.at(j)=" << vstrcheck.at(j) << endl;
	      if(aurostd::substring2bool(keyword,"KPOINTS")) {
		if(vstrcheck.at(j)=="A") vstrcheck.at(j)="AUTO";
		if(vstrcheck.at(j)=="G") vstrcheck.at(j)="GAMMA";
		if(vstrcheck.at(j)=="M") vstrcheck.at(j)="MONKHORST_PACK";
	      }	 
	      if(aurostd::substring2bool(keyword,"IGNORE_AFIX")) {
		;
	      }; // dummy load
	      if(aurostd::substring2bool(keyword,"CONVERT_UNIT_CELL")) {
		if(vstrcheck.at(j)=="SPRIM") vstrcheck.at(j)="STANDARD_PRIMITIVE";
		if(vstrcheck.at(j)=="STD_PRIM") vstrcheck.at(j)="STANDARD_PRIMITIVE";
		if(vstrcheck.at(j)=="STANDARD_PRIMITIVE") vstrcheck.at(j)="STANDARD_PRIMITIVE";
		if(vstrcheck.at(j)=="SCONV") vstrcheck.at(j)="STANDARD_CONVENTIONAL";
		if(vstrcheck.at(j)=="STD_CONV") vstrcheck.at(j)="STANDARD_CONVENTIONAL";
		if(vstrcheck.at(j)=="STANDARD_CONVENTIONAL") vstrcheck.at(j)="STANDARD_CONVENTIONAL";
		if(vstrcheck.at(j)=="NIGGLI") vstrcheck.at(j)="NIGGLI";
		if(vstrcheck.at(j)=="MINK") vstrcheck.at(j)="MINKOWSKI";
		if(vstrcheck.at(j)=="MINKOWSKI") vstrcheck.at(j)="MINKOWSKI";
		if(vstrcheck.at(j)=="INCELL") vstrcheck.at(j)="INCELL";
		if(vstrcheck.at(j)=="COMPACT") vstrcheck.at(j)="COMPACT";
		if(vstrcheck.at(j)=="INCOMPACT") vstrcheck.at(j)="COMPACT";
		if(vstrcheck.at(j)=="INWIGNERSEITZ") vstrcheck.at(j)="WIGNERSEITZ";
		if(vstrcheck.at(j)=="WS") vstrcheck.at(j)="WIGNERSEITZ";
		if(vstrcheck.at(j)=="WIGNERSEITZ") vstrcheck.at(j)="WIGNERSEITZ";
		if(vstrcheck.at(j)=="C") vstrcheck.at(j)="CARTESIAN";
		if(vstrcheck.at(j)=="CART") vstrcheck.at(j)="CARTESIAN";
		if(vstrcheck.at(j)=="CARTESIAN") vstrcheck.at(j)="CARTESIAN";
		if(vstrcheck.at(j)=="F") vstrcheck.at(j)="FRACTIONAL";
		if(vstrcheck.at(j)=="FRAC") vstrcheck.at(j)="FRACTIONAL";
		if(vstrcheck.at(j)=="FRACTIONAL") vstrcheck.at(j)="FRACTIONAL";
		if(vstrcheck.at(j)=="D") vstrcheck.at(j)="DIRECT";
		if(vstrcheck.at(j)=="DIR") vstrcheck.at(j)="DIRECT";
		if(vstrcheck.at(j)=="DIRECT") vstrcheck.at(j)="DIRECT";
		if(vstrcheck.at(j)=="PRE") vstrcheck.at(j)="PRESERVE";
		if(vstrcheck.at(j)=="PRES") vstrcheck.at(j)="PRESERVE";
		if(vstrcheck.at(j)=="PRESERVE") vstrcheck.at(j)="PRESERVE";
	      }	
	      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: AFTER keyword=" << keyword << "   " << "vstrcheck.at(j)=" << vstrcheck.at(j) << endl;
	      vcontent.push_back(vstrcheck.at(j));
	    }
	  }
	}
	for(uint i=0;i<vcontent.size();i++)
	  content_string+=vcontent.at(i)+(i<vcontent.size()-1?",":"");
	aurostd::StringSubst(content_string,"=","_");aurostd::StringSubst(content_string,";",",");
	if(vcontent.size()) isentry=TRUE;
      } // aurostd_xoptionMULTI list
    }
    content_double=aurostd::string2utype<double>(content_string);
    content_int=aurostd::string2utype<int>(content_string);
    content_uint=aurostd::string2utype<uint>(content_string);
    xscheme=content_string;
    aurostd::string2tokens(xscheme,vxscheme,","); 
    if(LDEBUG) if(_option_DEFAULT==aurostd_xoptionMULTI) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::options2entry: vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;

    preserved=FALSE;
    for(uint i=0;i<vxscheme.size()&&!preserved;i++) preserved=(vxscheme.at(i)=="PRESERVED");
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: isentry=" << (isentry?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: content_string=\"" << content_string << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: content_double=\"" << content_double << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: content_int=\"" << content_int << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: content_uint=\"" << content_uint << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: option=" << (option?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: preserved=" << (preserved?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: xscheme=" << xscheme << endl << endl;
    if(isentry && content_string.empty()) {
      // check for errors
      cerr << "ERROR - aurostd::xoption::options2entry: content_string=" <<  content_string << endl;
      cerr << "ERROR - aurostd::xoption::options2entry: content_double=" <<  content_double << endl;
      cerr << "ERROR - aurostd::xoption::options2entry: content_int=" <<  content_int << endl;
      cerr << "ERROR - aurostd::xoption::options2entry: content_uint=" <<  content_uint << endl;
      cerr << "ERROR - aurostd::xoption::options2entry: keyword=" << keyword  << endl;
      cerr << "ERROR - aurostd::xoption::options2entry: isentry=" << isentry  << endl;
      // check for errors
      cout << "ERROR - aurostd::xoption::options2entry: content_string=" <<  content_string << endl;
      cout << "ERROR - aurostd::xoption::options2entry: content_double=" <<  content_double << endl;
      cout << "ERROR - aurostd::xoption::options2entry: content_int=" <<  content_int << endl;
      cout << "ERROR - aurostd::xoption::options2entry: content_uint=" <<  content_uint << endl;
      cout << "ERROR - aurostd::xoption::options2entry: keyword=" << keyword  << endl;
      cout << "ERROR - aurostd::xoption::options2entry: isentry=" << isentry  << endl;
      exit(0);
    }
    if(LDEBUG) cerr << "DEBUG - aurostd::xoption::options2entry: END" << endl;
    //  exit(0);
    // return isentry;
  }

  void xoption::scheme2scheme(char c,string s) {
    for(uint i=0;i<vxscheme.size();i++) {
      if(vxscheme.at(i).at(0)==c ||
	 vxscheme.at(i).at(0)==aurostd::tolower(c) ||
	 vxscheme.at(i).at(0)==aurostd::toupper(c)) {
	xscheme=s;
      }
    }
  }

  void xoption::scheme2scheme(string s1,string s2) {
    for(uint i=0;i<vxscheme.size();i++) {
      if(vxscheme.at(i)==s1 ||
	 vxscheme.at(i)==aurostd::tolower(s1) ||
	 vxscheme.at(i)==aurostd::toupper(s1)) {
	xscheme=s2;
	//  for(uint i=0;i<vxscheme.size();i++) if(vxscheme.at(i)==s1) scheme=s2;
      }
    }
  }

  bool xoption::isscheme(string check) const { //CO 180101
    string a,b;
    for(uint i=0;i<vxscheme.size();i++) {
      a=aurostd::toupper(vxscheme.at(i));b=aurostd::toupper(check);                                  // shortcuts
      aurostd::StringSubst(a,"GAMMA","G");aurostd::StringSubst(b,"GAMMA","G");                      // shortcuts
      aurostd::StringSubst(a,"MONKHORST_PACK","M");aurostd::StringSubst(b,"MONKHORST_PACK","M");    // shortcuts
      aurostd::StringSubst(a,"MP","M");aurostd::StringSubst(b,"MP","M");                            // shortcuts
      aurostd::StringSubst(a,"AUTO","A");aurostd::StringSubst(b,"AUTO","A");                        // shortcuts
      if(a==b) return TRUE;
    }
    return FALSE;
  }

  bool xoption::refresh() {
    content_string="";
    for(uint i=0;i<vxscheme.size();i++)
      content_string+=vxscheme.at(i)+(i<vxscheme.size()-1?",":"");
    aurostd::StringSubst(content_string,"=","_");aurostd::StringSubst(content_string,";",",");
    content_double=aurostd::string2utype<double>(content_string);
    content_int=aurostd::string2utype<int>(content_string);
    content_uint=aurostd::string2utype<uint>(content_string);
    xscheme=content_string;
    return TRUE;
  }

  // [OBSOLETE] uint xoption::addscheme(string _xscheme)   { return opscheme(_xscheme,TRUE); }
  // [OBSOLETE] uint xoption::purgescheme(string _xscheme) { return opscheme(_xscheme,FALSE); }
  uint xoption::push(string _xscheme)        { return opscheme(_xscheme,TRUE); }
  uint xoption::pop(string _xscheme)         { return opscheme(_xscheme,FALSE); }

  uint xoption::opscheme(string _xscheme,bool operation) {
    if(operation==TRUE) {
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::opscheme: ADD=" << aurostd::toupper(_xscheme) << endl;
      if(LDEBUG) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: ADD_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
      vxscheme.push_back(aurostd::toupper(_xscheme));
      if(LDEBUG) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: ADD_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
    } else {
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::opscheme: PURGE=" << aurostd::toupper(_xscheme) << endl;
      if(LDEBUG) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: PURGE_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
      vector<string> _vxscheme(vxscheme);
      vxscheme.clear();
      for(uint i=0;i<_vxscheme.size();i++) {
	if(aurostd::toupper(_vxscheme.at(i))!=aurostd::toupper(_xscheme)) vxscheme.push_back(_vxscheme.at(i));
      }
      if(LDEBUG) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: PURGE_AFTER vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
    }
    refresh();
    return vxscheme.size();
  }
 
  bool xoption::flag(string _xscheme,bool operation) {
    // [OBSOLETE]    if(operation) push(_xscheme);
    // [OBSOLETE]    if(!operation) pop(_xscheme);
    if(operation) opscheme(_xscheme,TRUE);  // push
    if(!operation) opscheme(_xscheme,FALSE);  // pop
    return operation;
  }
  
  bool xoption::flag(string xscheme) const { return isscheme(xscheme);  } // same as ischeme
  
  bool xoption::flag(void) const {  // same as ischeme
    if(vxscheme.size()>0) return TRUE;
    return FALSE;
  }
  
  string xoption::getattachedscheme(string xscheme) const {
    if(vxsghost.size()==0) return "";
    for(uint i=0;i<vxsghost.size()-1;i+=2) {
      if(LDEBUG) cerr << i << " --- [" << aurostd::toupper(xscheme) << "] --- [" << aurostd::toupper(vxsghost.at(i)) << "] --- [" << aurostd::toupper(vxsghost.at(i+1)) << "]" << endl;
      if(aurostd::toupper(xscheme)==aurostd::toupper(vxsghost.at(i))) 
	return vxsghost.at(i+1);
    }
    return "";
  }
  template<class utype> utype xoption::getattachedutype(string xscheme) {
    return aurostd::string2utype<utype>(getattachedscheme(xscheme));
  }

  uint xoption::opattachedscheme(string _xscheme,string attached,bool operation) {
    if(operation==TRUE) {
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::opattachedscheme: ADD=" << aurostd::toupper(_xscheme) << endl;
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::opattachedscheme: GHOST=" << attached << endl;
      vxsghost.push_back(aurostd::toupper(_xscheme));
      vxsghost.push_back(attached);
    }  else {
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::opattachedscheme: PURGE=" << aurostd::toupper(_xscheme) << endl;
      vector<string> _vxsghost(vxsghost);
      vxsghost.clear();
      for(uint i=0;i<_vxsghost.size();i+=2) {
	if(aurostd::toupper(_vxsghost.at(i))!=aurostd::toupper(_xscheme)) {
	  vxsghost.push_back(_vxsghost.at(i));vxsghost.push_back(_vxsghost.at(i+1));
	}
      }
      if(LDEBUG) for(uint i=0;i<vxsghost.size();i++) cerr << "PURGEATTACHED_AFTER vxsghost.at(" << i << ")=" << vxsghost.at(i) << endl;
    } 
    refresh();
    return vxsghost.size();
  }  
    
  uint xoption::addattachedscheme(string _xscheme,string attached,bool operation) {
    if(operation) return opattachedscheme(_xscheme,attached,TRUE);
    return vxsghost.size();
   }
  
  // [OBSOLETE] uint xoption::purgeattachedscheme(string _xscheme) {
  // [OBSOLETE]   return opattachedscheme(_xscheme,"",FALSE);
  // [OBSOLETE] }
    
 uint xoption::push_attached(string _xscheme,string attached) {
    return opattachedscheme(_xscheme,attached,TRUE);
  }
  
  uint xoption::pop_attached(string _xscheme) {
    return opattachedscheme(_xscheme,"",FALSE);
  }
  

 
  bool xoption::args2addattachedscheme(vector<string> &argv,const string _xscheme,const string& _s_search,string string_default) {
    vector<string> cmds;
    return args2addattachedscheme(argv,cmds,_xscheme,_s_search,string_default);
  }
  
  bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string xscheme,const string& _s_search,string string_default) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string s_search(_s_search);
    if(aurostd::args2attachedflag(argv,cmds,s_search)) {
      flag(xscheme,TRUE);
      addattachedscheme(xscheme,aurostd::args2attachedstring(argv,s_search,string_default),TRUE);
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " attached=" << aurostd::args2attachedstring(argv,s_search,string_default) << endl;
      return TRUE;
    } 
    aurostd::StringSubst(s_search,"=","");
    if(aurostd::args2flag(argv,cmds,s_search)) {
      //    cerr << aurostd::args2string(argv,s_search,string_default) << endl;
      flag(xscheme,TRUE);
      // [OBSOLETE]      addattachedscheme(xscheme,string_default,TRUE);
      // [OBSOLETE] if(LDEBUG) cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " taking=" << string_default << endl;
      addattachedscheme(xscheme,aurostd::args2string(argv,s_search,string_default),TRUE);
      if(LDEBUG) cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " taking=" << aurostd::args2string(argv,s_search,string_default) << endl;
      return TRUE;
    }
    return FALSE;
  }

  bool xoption::args2addattachedscheme(vector<string> &argv,const string xscheme,const string& _s_search,char const* string_default) {
    return args2addattachedscheme(argv,xscheme,_s_search,string(string_default));
  }
  
  bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string xscheme,const string& _s_search,char const* string_default) {
    return args2addattachedscheme(argv,cmds,xscheme,_s_search,string(string_default));
  }
  
   template<class utype> bool xoption::args2addattachedscheme(vector<string> &argv,const string xscheme,const string& _s_search,utype utype_default) {
    return args2addattachedscheme(argv,xscheme,_s_search,aurostd::utype2string(utype_default));
  }
  
  template<class utype> bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string xscheme,const string& _s_search,utype utype_default) {
    return args2addattachedscheme(argv,cmds,xscheme,_s_search,aurostd::utype2string(utype_default));
  }
    
}

#endif  // _AUROSTD_XOPTION_CPP_


// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************
