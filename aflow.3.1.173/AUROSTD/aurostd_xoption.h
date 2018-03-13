// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
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
    bool isentry;              // is the entry available
    string content_string;     // the content
    double content_double;     // the content
    int content_int;           // the content
    uint content_uint;         // the content
    bool option;               // the output
    bool option_default;       // the default 
    string scheme;             // the content
    vector<string> vscheme;    // tokenized "," content
    vector<string> vsghost;    // tokenized "," content
    bool preserved;            // the output
    bool LDEBUG;               // to print or not to print
    // LOAD BOOLS 
    // FUNCTIONS
    void options2entry(string,string,int=aurostd_xoptionONOFF,string="");
    void string2scheme(char,string);
    void string2scheme(string,string);
    bool isscheme(string) const; // check if available //CO 180101
    // [OBSOLETE] bool isscheme(string,string); // check one or more
    void addscheme(string); // add if appropriate
    // [OBSOLETE] void addscheme(string,string); // add one or more
    uint purgescheme(string); // remove if appropriate    returns number of vscheme
    // [OBSOLETE] void purgescheme(string,string); // purge one or more
    // for plain flags
    bool flag(string,bool); // if bool=TRUE/FALSE => add/remove "string" 
    bool flag(string) const; // interrogate=TRUE/FALSE, same as ischeme //CO 180101
    bool flag(void) const; // return if there is any scheme inside //CO 180101
    bool addattachedscheme(string scheme,string attached,bool flag);
    string getattachedscheme(string scheme) const; //CO 180101
    template<class utype> utype getattachedutype(string scheme);
    uint purgeattachedscheme(string check);  //   returns number of vghost
    bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,string string_default); 
    bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,string string_default);
    bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,char const* string_default); 
    bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,char const* string_default);
    template<class utype> bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,utype utype_default); 
    template<class utype> bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,utype utype_default);
    bool refresh(void);
  private:                                             //
    void free();                                        // free space
    void copy(const xoption& b);                        //
  };
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XOPTION_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                                                                         *
// ***************************************************************************

