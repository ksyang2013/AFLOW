// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 2013-2014
// added template<class utype> bool xoption::args2addattachedscheme SC 2017

#ifndef _AUROSTD_XOPTION_CPP_
#define _AUROSTD_XOPTION_CPP_
//#include "aflow.h"
using std::ostream;

namespace aurostd {
  // ***************************************************************************

  // constructure
  xoption::xoption() {
    isentry=FALSE; 
    content_string="";
    content_double=0.0;
    content_int=0;
    content_uint=0;
    option=FALSE;
    option_default=FALSE;
    scheme="";
    vscheme.clear();
    vsghost.clear();
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
    isentry=b.isentry;
    content_string=b.content_string;
    content_double=b.content_double;
    content_int=b.content_int;
    content_uint=b.content_uint;
    option=b.option;
    option_default=b.option_default;
    scheme=b.scheme;
    vscheme.clear();for(uint i=0;i<b.vscheme.size();i++) vscheme.push_back(b.vscheme.at(i));
    vsghost.clear();for(uint i=0;i<b.vsghost.size();i++) vsghost.push_back(b.vsghost.at(i));
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
    for(uint i=0;i<a.vscheme.size();i++)
      oss << a.vscheme.at(i) << (i<a.vscheme.size()-1?",":"");
    return oss;
  }
  
  void xoption::clear() {
    xoption aflow_option_temp;
    copy(aflow_option_temp);
  }

  // **************************************************************************
  void xoption::options2entry(string options_FILE,string _keyword,int _option_DEFAULT,string scheme_DEFAULT) {
    clear();
    LDEBUG=FALSE;
    bool option_DEFAULT=FALSE; 
    if(_option_DEFAULT==0) option_DEFAULT=FALSE; // it is a int.. it might be -1
    if(_option_DEFAULT==1) option_DEFAULT=TRUE; // it is a int.. it might be -1
    isentry=option_DEFAULT;option=option_DEFAULT;content_string=scheme_DEFAULT;scheme=scheme_DEFAULT;preserved=FALSE;   // DEFAULT
    if(LDEBUG) cerr << "DEBUG: KBIN_options2entry" << endl;
    if(LDEBUG) cerr << "DEBUG: _keyword=\"" << _keyword << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG: option_DEFAULT=" << (option_DEFAULT?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG: scheme_DEFAULT=\"" << scheme_DEFAULT << "\"" << endl;
    // start the scan
    string keyword;
    vector<string> vkeyword;
    // tokenize the option
    aurostd::string2tokens(_keyword,vkeyword,"|"); 
    if(LDEBUG) for(uint i=0;i<vkeyword.size();i++) cerr << "\"" << vkeyword.at(i) << "\"" << endl;
    // loop through the scan
    if(vkeyword.size()>0) {
      // some default
      keyword=vkeyword.at(0);
      for(uint i=0;i<vkeyword.size();i++) 
	if(aurostd::substring2bool(options_FILE,vkeyword.at(i),TRUE))
	  keyword=vkeyword.at(i);
      // found one keyword
      if(LDEBUG) cerr << "DEBUG: keyword=\"" << keyword << "\"" << endl;
      // LOOK FOR EXIST/!EXIST ENTRY
      if(_option_DEFAULT==aurostd_xoptionONOFF) {
	isentry=aurostd::substring2bool(options_FILE,keyword,TRUE);
	if(isentry) {option=TRUE;content_string="ON";}
	if(!isentry) {option=FALSE;content_string="OFF";}
      } // aurostd_xoptionONOFF exit/~exit
      // LOOK FOR ON/OFF MODE WITH strings/schemes.
      if(_option_DEFAULT==0 || _option_DEFAULT==1) {
	if(LDEBUG) cerr << "DEBUG: LOOK FOR ON/OFF MODE WITH strings/schemes" << endl;
	// start the scan
	isentry=aurostd::substring2bool(options_FILE,keyword,TRUE);
	if(isentry && scheme_DEFAULT.empty()) {
	  content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,FALSE));
	  string saus=content_string;content_string="";
	  if(LDEBUG) cerr << "DEBUG: found saus=" << saus << endl;
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
	if(LDEBUG) cerr << "DEBUG: scheme_DEFAULT=\"" << scheme_DEFAULT << "\"" << endl;
	if(LDEBUG) cerr << "DEBUG: scheme_DEFAULT.empty()=" << scheme_DEFAULT.empty() << endl;
	if(isentry && !scheme_DEFAULT.empty()) {
	  if(LDEBUG) cerr << "DEBUG: SCHEME MODE" << endl;
	  content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,FALSE));
	  option=isentry;
	}
	if(isentry && (scheme_DEFAULT.empty() && content_string.empty())) {
	  if(LDEBUG) cerr << "DEBUG: SCHEME MODE EMPTY DEFAULT STILL EMPTY CONTENT" << endl;
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
	      if(LDEBUG) cerr << "DEBUG: BEFORE keyword=" << keyword << "   " << "vstrcheck.at(j)=" << vstrcheck.at(j) << endl;
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
	      if(LDEBUG) cerr << "DEBUG: AFTER keyword=" << keyword << "   " << "vstrcheck.at(j)=" << vstrcheck.at(j) << endl;
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
    scheme=content_string;
    aurostd::string2tokens(scheme,vscheme,","); 
    if(LDEBUG) if(_option_DEFAULT==aurostd_xoptionMULTI) for(uint i=0;i<vscheme.size();i++) cerr << "DEBUG: vscheme.at(" << i << ")=" << vscheme.at(i) << endl;

    preserved=FALSE;
    for(uint i=0;i<vscheme.size()&&!preserved;i++) preserved=(vscheme.at(i)=="PRESERVED");
    if(LDEBUG) cerr << "DEBUG: isentry=" << (isentry?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG: content_string=\"" << content_string << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG: content_double=\"" << content_double << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG: content_int=\"" << content_int << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG: content_uint=\"" << content_uint << "\"" << endl;
    if(LDEBUG) cerr << "DEBUG: option=" << (option?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG: preserved=" << (preserved?"TRUE":"FALSE") << endl;
    if(LDEBUG) cerr << "DEBUG: scheme=" << scheme << endl << endl;
    if(isentry && content_string.empty()) {   // check for errors
      cerr << "ERROR: KBIN_options2entry:  content_string=" <<  content_string << endl;
      cerr << "ERROR:                      content_double=" <<  content_double << endl;
      cerr << "ERROR:                      content_int=" <<  content_int << endl;
      cerr << "ERROR:                      content_uint=" <<  content_uint << endl;
      cerr << "ERROR:                      keyword=" << keyword  << endl;
      cerr << "ERROR:                      isentry=" << isentry  << endl;
      cout << "ERROR: KBIN_options2entry:  content_string=" <<  content_string << endl;
      cout << "ERROR:                      content_double=" <<  content_double << endl;
      cout << "ERROR:                      content_int=" <<  content_int << endl;
      cout << "ERROR:                      content_uint=" <<  content_uint << endl;
      cout << "ERROR:                      keyword=" << keyword  << endl;
      cout << "ERROR:                      isentry=" << isentry  << endl;
      exit(0);
    }
    //  exit(0);
    // return isentry;
  }

  void xoption::string2scheme(char c,string s) {
    for(uint i=0;i<vscheme.size();i++) if(vscheme.at(i).at(0)==c || vscheme.at(i).at(0)==aurostd::tolower(c) || vscheme.at(i).at(0)==aurostd::toupper(c)) scheme=s;
  }

  void xoption::string2scheme(string s1,string s2) {
    for(uint i=0;i<vscheme.size();i++) if(vscheme.at(i)==s1 || vscheme.at(i)==aurostd::tolower(s1) || vscheme.at(i)==aurostd::toupper(s1)) scheme=s2;
    //  for(uint i=0;i<vscheme.size();i++) if(vscheme.at(i)==s1) scheme=s2;
  }

  bool xoption::isscheme(string check) const { //CO 180101
    string a,b;
    for(uint i=0;i<vscheme.size();i++) {
      a=aurostd::toupper(vscheme.at(i));b=aurostd::toupper(check);                                  // shortcuts
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
    for(uint i=0;i<vscheme.size();i++)
      content_string+=vscheme.at(i)+(i<vscheme.size()-1?",":"");
    aurostd::StringSubst(content_string,"=","_");aurostd::StringSubst(content_string,";",",");
    content_double=aurostd::string2utype<double>(content_string);
    content_int=aurostd::string2utype<int>(content_string);
    content_uint=aurostd::string2utype<uint>(content_string);
    scheme=content_string;
    return TRUE;
  }

  void xoption::addscheme(string check) {
    if(LDEBUG) cerr << "ADD=" << aurostd::toupper(check) << endl;
    if(LDEBUG) for(uint i=0;i<vscheme.size();i++) cerr << "ADD_BEFORE vscheme.at(" << i << ")=" << vscheme.at(i) << endl;
    vscheme.push_back(aurostd::toupper(check));
    if(LDEBUG) for(uint i=0;i<vscheme.size();i++) cerr << "ADD_BEFORE vscheme.at(" << i << ")=" << vscheme.at(i) << endl;
    refresh();
  }

  uint xoption::purgescheme(string check) {
    if(LDEBUG) cerr << "PURGE=" << aurostd::toupper(check) << endl;
    if(LDEBUG) for(uint i=0;i<vscheme.size();i++) cerr << "PURGE_BEFORE vscheme.at(" << i << ")=" << vscheme.at(i) << endl;
    vector<string> _vscheme(vscheme);
    vscheme.clear();
    for(uint i=0;i<_vscheme.size();i++) {
      if(aurostd::toupper(_vscheme.at(i))!=aurostd::toupper(check)) vscheme.push_back(_vscheme.at(i));
    }
    if(LDEBUG) for(uint i=0;i<vscheme.size();i++) cerr << "PURGE_AFTER vscheme.at(" << i << ")=" << vscheme.at(i) << endl;
    refresh();
    return vscheme.size();
  }

  bool xoption::flag(string scheme,bool flag) {
    if(flag) addscheme(scheme);
    if(!flag) purgescheme(scheme);
    return flag;
  }
  bool xoption::addattachedscheme(string scheme,string attached,bool _flag_) {
    if(_flag_) {
      // if(!flag(scheme)) flag(scheme,TRUE);
      vsghost.push_back(aurostd::toupper(scheme));
      vsghost.push_back(attached);}
    return _flag_;
  }

  string xoption::getattachedscheme(string scheme) const {
    if(vsghost.size()==0) return "";
    for(uint i=0;i<vsghost.size()-1;i+=2) {
      if(LDEBUG) cerr << i << " --- [" << aurostd::toupper(scheme) << "] --- [" << aurostd::toupper(vsghost.at(i)) << "] --- [" << aurostd::toupper(vsghost.at(i+1)) << "]" << endl;
      if(aurostd::toupper(scheme)==aurostd::toupper(vsghost.at(i))) 
	return vsghost.at(i+1);
    }
    return "";
  }
  template<class utype> utype xoption::getattachedutype(string scheme) {
    return aurostd::string2utype<utype>(getattachedscheme(scheme));
  }
  
  uint xoption::purgeattachedscheme(string check) {
    if(LDEBUG) cerr << "PURGEATTACHED=" << aurostd::toupper(check) << endl;
    if(LDEBUG) for(uint i=0;i<vsghost.size();i++) cerr << "PURGEATTACHED_BEFORE vsghost.at(" << i << ")=" << vsghost.at(i) << endl;
    vector<string> _vsghost(vsghost);
    vsghost.clear();
    for(uint i=0;i<_vsghost.size();i+=2) {
      if(aurostd::toupper(_vsghost.at(i))!=aurostd::toupper(check)) {
	vsghost.push_back(_vsghost.at(i));vsghost.push_back(_vsghost.at(i+1));
      }
    }
    if(LDEBUG) for(uint i=0;i<vsghost.size();i++) cerr << "PURGEATTACHED_AFTER vsghost.at(" << i << ")=" << vsghost.at(i) << endl;
    refresh();
    return vsghost.size();
  }
  
  bool xoption::args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,string string_default) {
    vector<string> cmds;
    return args2addattachedscheme(argv,cmds,scheme,_s_search,string_default);
  }
  
  bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,string string_default) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string s_search(_s_search);
    if(aurostd::args2attachedflag(argv,cmds,s_search)) {
      flag(scheme,TRUE);
      addattachedscheme(scheme,aurostd::args2attachedstring(argv,s_search,string_default),TRUE);
      if(LDEBUG) cerr << "xoption::args2addscheme: scheme=" << scheme << " s_search=" << s_search << " attached=" << aurostd::args2attachedstring(argv,s_search,string_default) << endl;
      return TRUE;
    } 
    aurostd::StringSubst(s_search,"=","");
    if(aurostd::args2flag(argv,cmds,s_search)) {
      //    cerr << aurostd::args2string(argv,s_search,string_default) << endl;
      flag(scheme,TRUE);
      // [OBSOLETE]      addattachedscheme(scheme,string_default,TRUE);
      // [OBSOLETE] if(LDEBUG) cerr << "xoption::args2addscheme: scheme=" << scheme << " s_search=" << s_search << " taking=" << string_default << endl;
      addattachedscheme(scheme,aurostd::args2string(argv,s_search,string_default),TRUE);
      if(LDEBUG) cerr << "xoption::args2addscheme: scheme=" << scheme << " s_search=" << s_search << " taking=" << aurostd::args2string(argv,s_search,string_default) << endl;
      return TRUE;
    }
    return FALSE;
  }

  bool xoption::args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,char const* string_default) {
    return args2addattachedscheme(argv,scheme,_s_search,string(string_default));
  }
  
  bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,char const* string_default) {
    return args2addattachedscheme(argv,cmds,scheme,_s_search,string(string_default));
  }
  
   template<class utype> bool xoption::args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,utype utype_default) {
    return args2addattachedscheme(argv,scheme,_s_search,aurostd::utype2string(utype_default));
  }
  
  template<class utype> bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,utype utype_default) {
    return args2addattachedscheme(argv,cmds,scheme,_s_search,aurostd::utype2string(utype_default));
  }
  
  
  bool xoption::flag(string scheme) const {  // same as ischeme
    return isscheme(scheme); 
  }

  bool xoption::flag(void) const {  // same as ischeme
    if(vscheme.size()>0) return TRUE;
    return FALSE;
  }


  
}

#endif  // _AUROSTD_XOPTION_CPP_


// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2015              *
// *                                                                        *
// **************************************************************************
