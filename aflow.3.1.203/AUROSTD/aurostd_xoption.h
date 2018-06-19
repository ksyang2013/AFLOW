// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 2013-2014
// added template<class utype> bool xoption::args2addattachedscheme SC 2017

#ifndef _AUROSTD_XOPTION_H_
#define _AUROSTD_XOPTION_H_

// --------------------------------------------------------------------------
// general flag for xoption to take/manipulate options
#define aurostd_xoptionONOFF int(-1)
#define aurostd_xoptionMULTI int(-2)

using std::ostream;
using std::vector;
using std::string;

namespace aurostd {
  // namespace aurostd
  class xoption {
  public:
    // trivial constructurs/destuctors/operators
    xoption();                                         // default, just allocate
    ~xoption();                                        // kill everything
    xoption(const xoption& b);                         // constructor copy
    const xoption& operator=(const xoption &b);        // copy
    friend std::ostream& operator<<(std::ostream&,const xoption&);       // ostream
    void clear(void);                                  // clear
    // CONTENT
    string keyword;            // the keyword found (we can provide a bunch with |) //CO 180404
    bool isentry;              // is the entry available
    string content_string;     // the content
    double content_double;     // the content
    int content_int;           // the content
    uint content_uint;         // the content
    bool option;               // the output
    bool option_default;       // the default 
    string xscheme;            // the content
    vector<string> vxscheme;   // tokenized "," content
    vector<string> vxsghost;   // tokenized "," content
    bool preserved;            // the output
    bool LDEBUG;               // to print or not to print
    // LOAD BOOLS FUNCTIONS
    void options2entry(string,string,int=aurostd_xoptionONOFF,string="");
    void scheme2scheme(char,string);
    void scheme2scheme(string,string);
    bool isscheme(string) const; // check if available //CO 180101
    // [OBSOLETE] uint addscheme(string);      // add scheme then returns vscheme.size()
    // [OBSOLETE] uint purgescheme(string);    // remove scheme then returns vscheme.size()
    uint opscheme(string,bool);  // add/remove scheme then returns vscheme.size()
    uint push(string);           // add scheme then returns vscheme.size()
    uint pop(string);            // remove scheme then returns vscheme.size()
    // for plain flags
    bool flag(string,bool);      // if bool=TRUE/FALSE => add/remove "string" 
    bool flag(string) const;     // interrogate=TRUE/FALSE, same as ischeme // CO 180101
    bool flag(void) const;       // return if there is any scheme inside // CO 180101
    // attached stuff..
    uint opattachedscheme(string,string,bool);                       // add/remove attached_scheme if flag=TRUE, then returns vghost.size()
    uint addattachedscheme(string scheme,string attached,bool flag); // add attached_scheme if flag=TRUE, then returns vghost.size()
    // [OBSOLETE] uint purgeattachedscheme(string check);            // remove attached_scheme, then returns vghost.size() - same as pop_attached
    uint push_attached(string scheme,string attached);               // add attached_scheme, then returns vghost.size() - like addattachedscheme with flag=TRUE
    uint pop_attached(string check);                                 // remove attached_scheme, then returns vghost.size()
    string getattachedscheme(string scheme) const; // CO 180101
    template<class utype> utype getattachedutype(string scheme);
    bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,string string_default); 
    bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,string string_default);
    bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,char const* string_default); 
    bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,char const* string_default);
    template<class utype> bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,utype utype_default); 
    template<class utype> bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,utype utype_default);
    bool refresh(void);
  private:                                              //
    void free();                                        // free space
    void copy(const xoption& b);                        //
  };
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XOPTION_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

