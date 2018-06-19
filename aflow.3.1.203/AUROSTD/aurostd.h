// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_MAIN_H_
#define _AUROSTD_MAIN_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <fstream>
#include <grp.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <pthread.h>
#include <pwd.h>
#include <queue>
#include <sstream>
#include <stdarg.h>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <sys/sysctl.h>
// #include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <sys/wait.h>
#include <time.h>
#include <typeinfo>
#include <unistd.h>
#include <vector>
#include <list> //CO 170806 - need for POCC
#include <netdb.h>  //CO 180321 - frisco needs for AFLUX

#ifdef _USE_AFLOW_H_
//#include "aflow.h"
#endif

using std::abs;
using std::cerr;
using std::cin;
using std::cout;
using std::deque;
using std::endl;
using std::ends;
using std::ifstream;
using std::ios_base;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;
using std::vector;

#ifndef SWAP
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#endif

#ifndef AUROSTD_NAN
#define AUROSTD_NAN 1E9
#endif

#ifndef AUROSTD_DEFAULT_PRECISION
#define AUROSTD_DEFAULT_PRECISION 20
#endif

//CO 171215 - more global control
//this is for PRINTING, not for algorithmic zeroing
#ifndef AUROSTD_ROUNDOFF_TOL
#define AUROSTD_ROUNDOFF_TOL 1e-6
#endif

//CO 171215 - more global control
//this is SMALLER than _ZERO_TOL_ in aflow.h 
//(aurostd default is less conservative to match read/write roundoff error)
#ifndef AUROSTD_IDENTITY_TOL
#define AUROSTD_IDENTITY_TOL 1e-6
#endif

//CO 171002 - USEFUL!
#ifndef AUROSTD_MAX_UINT
#define AUROSTD_MAX_UINT std::numeric_limits<uint>::max()
#endif

//CO 171002 - USEFUL!
#ifndef AUROSTD_MAX_DOUBLE
#define AUROSTD_MAX_DOUBLE std::numeric_limits<double>::max()
#endif

//CO 180101 - stream2stream modes
#ifndef DEFAULT_STREAM
#define DEFAULT_STREAM 'D'
#endif
#ifndef FIXED_STREAM
#define FIXED_STREAM 'F'
#endif
#ifndef SCIENTIFIC_STREAM
#define SCIENTIFIC_STREAM 'S'
#endif


#define AUROSTD_ZIP_BIN "xz"
#define AUROSTD_ZIP_EXT ".xz"

#define __GNU_CPP
#ifndef __xprototype
#ifdef GNU
#define __xprototype __attribute__((const))
#else
#define __xprototype
#endif
#endif

#include "aurostd_xscalar.h"
#include "aurostd_xcomplex.h"
#include "aurostd_xvector.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xtensor.h"
#include "aurostd_xrandom.h"
#include "aurostd_xoption.h"
#include "aurostd_argv.h"
#include "aurostd_xcombos.h"

using aurostd::min;
using aurostd::max;
using aurostd::mod;
using aurostd::_isodd;
using aurostd::_iseven;
using aurostd::_isfloat;
using aurostd::_iscomplex;
using aurostd::sign;
using aurostd::nint;
using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::clear;
using aurostd::xvector;
using aurostd::xtensor3;
using aurostd::xtensor4;
using aurostd::xtensor5;
using aurostd::xtensor6;
using aurostd::xtensor7;
using aurostd::xtensor8;
using aurostd::xoption;

using aurostd::xvector;
using aurostd::xmatrix;
using aurostd::xcombos;

#ifndef uint
typedef unsigned uint;
#endif

// ----------------------------------------------------------------------------

#ifndef __STRICT_ANSI__
#define __STRICT_ANSI__
#endif

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#ifndef OFF
#define OFF FALSE
#endif
#ifndef ON
#define ON TRUE
#endif


#define _AUROSTD_XLIBS_ERROR_    string("ERROR - AUROSTD_XLIBS++ [***]:")
#define _AUROSTD_XLIBS_WARNING_  string("WARNING - AUROSTD_XLIBS++ [***]:")

#define EMPTY_WORDING string("blablabla")

//extern bool QUIET,DEBUG;
//extern class _XHOST XHOST;
#include "../aflow.h"
#include "../SQLITE/sqlite3.h"

template<class utype> std::ostream& operator<<(std::ostream&,const std::vector<utype>&);// __xprototype;
template<class utype> std::ostream& operator<<(std::ostream&,const std::deque<utype>&);// __xprototype;

// ----------------------------------------------------------------------------
// TIME stuff
namespace aurostd {
  int get_day(void);
  int get_month(void);
  int get_year(void);
  long int get_date(void);
  int get_hour(void);
  int get_min(void);
  int get_sec(void);
  long double get_seconds(void);
  long double get_seconds(long double reference_seconds);
  long double get_delta_seconds(long double& seconds_begin);
  long double get_useconds(void);
  long double get_useconds(long double reference_useconds);
  long double get_delta_useconds(long double& useconds_begin);
  string get_time(void);
  string get_datetime(void);
  string get_datetime_formatted(bool separate_date_time=true);
  bool beep(uint=2000,uint=100); // standard values
}
// ----------------------------------------------------------------------------
namespace aurostd {
  void sizes(void) __xprototype;
  // aflow_aurostd.cpp
  //int Nint(double x);
  //int Sign(const double& x);
  //int SignNoZero(const double& x);
  template<class utype> void aswap(utype &a,utype &b);                     // some dumb algebra
  template<class utype> utype max(const std::vector<utype> vec);           // some dumb algebra
  template<class utype> utype max(const std::deque<utype> vec);            // some dumb algebra
  template<class utype> utype max(const std::vector<vector<utype> >mat);   // some dumb algebra
  template<class utype> utype min(const std::vector<utype> vec);           // some dumb algebra
  template<class utype> utype min(const std::deque<utype> vec);            // some dumb algebra
  template<class utype> utype min(const std::vector<vector<utype> >mat);   // some dumb algebra
  template<class utype> utype sum(const std::vector<utype> vec);           // some dumb algebra
  template<class utype> utype sum(const std::deque<utype> vec);            // some dumb algebra
  template<class utype> utype sum(const std::vector<vector<utype> >mat);   // some dumb algebra
  template<class utype> utype mean(const std::vector<utype> vec);          // some dumb algebra
  template<class utype> utype mean(const std::deque<utype> vec);           // some dumb algebra
  template<class utype> utype mean(const std::vector<vector<utype> >mat);  // some dumb algebra
  template<class utype> vector<utype> reset(std::vector<utype>& v);        // some dumb algebra
  template<class utype> deque<utype> reset(std::deque<utype>& v);          // some dumb algebra
  template<class utype> vector<vector<utype> > reset(std::vector<std::vector<utype> > m);  // some dumb algebra
  template<class utype> vector<utype> clear(std::vector<utype>& v);        // some dumb algebra
  template<class utype> deque<utype> clear(std::deque<utype>& v);          // some dumb algebra
  template<class utype> vector<vector<utype> > clear(std::vector<std::vector<utype> > m);  // some dumb algebra
  template<class utype> void random_shuffle(std::vector<utype>& v);        // some dumb algebra
  template<class utype> void random_shuffle(std::deque<utype>& v);         // some dumb algebra
  template<class utype> bool identical(std::vector<utype> v1,std::vector<utype> v2,utype epsilon);
  template<class utype> bool identical(std::deque<utype> v1,std::deque<utype> v2,utype epsilon);
  bool identical(std::vector<int> v1,std::vector<int> v2,int epsilon);
  bool identical(std::deque<int> v1,std::deque<int> v2,int epsilon);
  template<class utype> bool identical(const std::vector<vector<utype> >& m1,const std::vector<vector<utype> >& m2,utype epsilon);
  string toupper(const string& in)  __xprototype;
  string tolower(const string& in)  __xprototype;
  char toupper(const char& in)  __xprototype;
  char tolower(const char& in)  __xprototype;
  int GetNumFields(const string& s);
  string GetNextVal(const string& s,int& id);
  string PaddedNumString(const int num,const int ndigits);
  template<class utype> string PaddedPRE(utype,int,string=" ");
  string PaddedPRE(string,int,string=" ");
  template<class utype> string PaddedPOST(utype,int,string=" ");
  string PaddedPOST(string,int,string=" ");
  template<class utype> string PaddedCENTER(utype,int,string=" ");
  string PaddedCENTER(string,int,string=" ");
  //about cleaning up strings
  string CleanStringASCII(const string& s) __xprototype;
  string CGI_StringClean(const string& stringIN) __xprototype;
  string RemoveWhiteSpaces(const string& s) __xprototype;
  string RemoveWhiteSpaces(const string& s, const char toogle) __xprototype;
  string RemoveWhiteSpacesFromTheBack(const string& s) __xprototype;
  string RemoveWhiteSpacesFromTheFront(const string& s) __xprototype;
  string RemoveWhiteSpacesFromTheFrontAndBack(const string& s) __xprototype;
  string RemoveSpaces(const string& s) __xprototype;
  string RemoveSpaces(const string& s, const char toogle) __xprototype;
  string RemoveSpacesFromTheBack(const string& s) __xprototype;
  string RemoveTabs(const string& s) __xprototype;
  string RemoveTabs(const string& s, const char toogle) __xprototype;
  string RemoveTabsFromTheBack(const string& s) __xprototype;
  string RemoveComments(const string& s) __xprototype;
  string RemoveCharacter(const string& s, const char character) __xprototype;
  string RemoveNumbers(const string& s) __xprototype;
  string RemoveRounding(const string& s) __xprototype;
  //string RemoveCharacter(const string& s, const char character, const char toogle) __xprototype;
  string RemoveSubStringFirst(const string& str_orig, const string& str_rm) __xprototype;
  string RemoveSubString(const string& str_orig, const string& str_rm) __xprototype;
  // about directories and file existing or not
  bool DirectoryMake(string Directory);
  bool SSH_DirectoryMake(string user, string machine,string Directory);
  bool DirectoryChmod(string chmod_string,string Directory);
  bool DirectoryLS(string Directory,vector<string> &vfiles);
  bool DirectoryLocked(string directory,string="LOCK");
  bool DirectorySkipped(string directory);
  bool DirectoryWritable(string directory);
  bool DirectoryUnwritable(string directory);
  string TmpFileCreate(string prefix);
  string TmpFileCreate(void);
  string TmpDirectoryCreate(string prefix);
  string TmpDirectoryCreate(void);
  string CleanFileName(string fileIN);
  string ProperFileName(string fileIN);
  bool CopyFile(string file_from,string file_to);
  //CO - START
  bool MatchCompressed(const string& CompressedFileName,const string& FileNameOUT);
  // [OBSOLETE]  bool DecompressFile(const string& CompressedFileName);
  bool efile2tempfile(string _FileNameIN, string& FileNameOUT); //CO 180220
  bool IsCompressed(string FileNameIn,string& FileNameOut);
  bool IsCompressed(string FileNameIn);
  string GetCompressionExtension(const string& CompressedFileName);
  //CO - END
  bool UncompressFile(const string& FileName,const string& command);  bool UncompressFile(const string& FileName); // with guess  
  bool CompressFile(const string& FileName,const string& command=AUROSTD_ZIP_BIN); // with default
  bool ZIP2ZIP(string dir,string from,string to,bool=TRUE);
  bool BZ2XZ(string dir,bool=TRUE); bool GZ2XZ(string dir,bool=TRUE);
  // [OBSOLETE]  bool BunzipFile(const string& FileName); bool BzipFile(const string& FileName);
  // [OBSOLETE]  bool GunzipFile(const string& FileName); bool GzipFile(const string& FileName); 
  // [OBSOLETE]  bool XunzipFile(const string& FileName); bool XzipFile(const string& FileName); 
  // [OBSOLETE]  bool UnzipFile(const string& FileName);  bool ZipFile(const string& FileName);
  bool FileExist(const string& FileName);  bool FileExist(const string& FileName,string &FileNameOut);
  bool EFileExist(const string& FileName); bool EFileExist(const string& FileName,string &FileNameOut);
  int  FileSize(const string& FileName);
  bool FileEmpty(const string& FileName);
  bool FileNotEmpty(const string& FileName);
  string FileToString(const string& FileName);
  void InFileExistCheck(const string& routine,const string& filename,ifstream& file_to_check,std::ostream& outf);
  bool IsCommandAvailable(const string& command);
  bool IsCommandAvailable(const string& command, string& position);
  bool IsCommandAvailableModify(string& command);
  bool CommandRequired(const string& command);
  bool CommandRequired(const string& command, string& position);
  bool IsExecutableAvailable(const string& executable);
  bool IsExecutableAvailable(const string& executable, string& position);
  bool ExecutableRequired(const string& executable);
  bool ExecutableRequired(const string& executable, string& position);
  // about cleaning, higiene is important
  void StringstreamClean(ostringstream& stream);
  void StringstreamClean(stringstream& stream);
  int FindIfStringInStream(const string& key,std::istream& instream);
  // remove comments
  
  // about printing
  void PrintMessageStream(ofstream& FileERROR,ostringstream& stream,const bool& quiet);
  void PrintMessageStream(std::ostream& FileERROR,ostringstream& stream,const bool& quiet);
  void PrintMessageStream(ofstream& FileERROR,ostringstream& stream,const bool& quiet,const bool& osswrite,std::ostream& oss);
  void PrintMessageStream(ostringstream& stream,const bool& quiet);
  void PrintErrorStream(ofstream& FileERROR,ostringstream& stream,const bool& quiet);
  void PrintErrorStream(std::ostream& FileERROR,ostringstream& stream,const bool& quiet);
  void PrintErrorStream(ofstream& FileERROR,ostringstream& stream,const bool& quiet,const bool& osswrite,std::ostream& oss);
  void PrintErrorStream(ostringstream& stream,const bool& quiet);
  void PrintWarningStream(ofstream& FileERROR,ostringstream& stream,const bool& quiet);
  void PrintWarningStream(std::ostream& FileERROR,ostringstream& stream,const bool& quiet);
  void PrintWarningStream(ofstream& FileERROR,ostringstream& stream,const bool& quiet,const bool& osswrite,std::ostream& oss);
  void PrintWarningStream(ostringstream& stream,const bool& quiet);

  void PrintMessageStream(ofstream& FileERROR,stringstream& stream,const bool& quiet);
  void PrintMessageStream(std::ostream& FileERROR,stringstream& stream,const bool& quiet);
  void PrintMessageStream(ofstream& FileERROR,stringstream& stream,const bool& quiet,const bool& osswrite,std::ostream& oss);
  void PrintMessageStream(stringstream& stream,const bool& quiet);
  void PrintErrorStream(ofstream& FileERROR,stringstream& stream,const bool& quiet);
  void PrintErrorStream(std::ostream& FileERROR,stringstream& stream,const bool& quiet);
  void PrintErrorStream(ofstream& FileERROR,stringstream& stream,const bool& quiet,const bool& osswrite,std::ostream& oss);
  void PrintErrorStream(stringstream& stream,const bool& quiet);
  void PrintWarningStream(ofstream& FileERROR,stringstream& stream,const bool& quiet);
  void PrintWarningStream(std::ostream& FileERROR,stringstream& stream,const bool& quiet);
  void PrintWarningStream(ofstream& FileERROR,stringstream& stream,const bool& quiet,const bool& osswrite,std::ostream& oss);
  void PrintWarningStream(stringstream& stream,const bool& quiet);

  // about executing
  bool execute(ostringstream& command);
  bool execute(stringstream& command);
  bool execute(string command);
  bool execute(vector<string> vcommand);
  bool execute(deque<string> dcommand);
#ifdef _stringcharstar_
  bool execute(char* command);
#endif
  // Execute and report
  string execute2string(ostringstream& command);
  string execute2string(stringstream& command);
  string execute2string(string command);
  vector<string> execute2string(vector<string> vcommand);
  deque<string> execute2string(deque<string> dcommand);
#ifdef _stringcharstar_
  string execute2string(char* command);
#endif
  template<class utype> utype execute2utype(ostringstream& command);
  template<class utype> utype execute2utype(stringstream& command);
  template<class utype> utype execute2utype(string command);
  template<class utype> vector<utype> execute2utype(vector<utype> vcommand);
  template<class utype> deque<utype> execute2utype(deque<utype> dcommand);
#ifdef _stringcharstar_
  template<class utype> utype execute2utype(char* command);
#endif
  // about sleeping
  unsigned int Sleep(unsigned int seconds);
  // about extracting from to files
  bool ExtractToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,string Keyword);
  bool ExtractToFileEXPLICIT(string StringIN,string FileNameOUTPUT,string Keyword);
  bool ExtractToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractToFileEXPLICIT(string StringIN,string FileNameOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,string Keyword);
  bool ExtractToStringEXPLICIT(string StringIN,string& StringOUTPUT,string Keyword);
  bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractToStringEXPLICIT(string StringIN,string& StringOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword);
  bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(stringstream StringStreamIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword);
  // take the last
  bool ExtractLastToStringstreamEXPLICIT(ifstream &FileIN,stringstream& StringstreamOUTPUT,string Keyword);
  bool ExtractLastToStringstreamEXPLICIT(ifstream &FileIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractLastToStringstreamEXPLICIT(stringstream StringStreamIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop);
  bool ExtractLastToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop);
  // take just after
  bool ExtractJustAfterToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,string Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,string Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(stringstream StringStreamIN,stringstream& StringstreamOUTPUT,string Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(string StringIN,string& StringOUTPUT,string Keyword_start);
  // about taking in istreams and stringstream and strings
  uint stream2vectorstring(std::istream& istreamIN,vector<string> &vstringout);
  uint stream2vectorstring(std::ifstream& ifstreamIN,vector<string> &vstringout);
  uint stream2vectorstring(stringstream& stringstreamIN,vector<string> &vstringout);
  uint string2vectorstring(const string& stringIN,vector<string> &vstringout,bool consecutive=false,bool trim_edges=false); //CO 170613, defaults to usual string2tokens() behavior
  vector<string> stream2vectorstring(std::istream& istreamIN);
  vector<string> stream2vectorstring(std::ifstream& ifstreamIN);
  vector<string> stream2vectorstring(stringstream& stringstreamIN);
  vector<string> string2vectorstring(const string& stringIN,bool consecutive=false,bool trim_edges=false);  //CO 170613, defaults to usual string2tokens() behavior
  string liststring2string(string="",string="",string="",string="",string="",string="",string="",string="",
			   string="",string="",string="",string="",string="",string="",string="",string="",
			   string="",string="",string="",string="",string="",string="",string="",string="",
			   string="",string="",string="",string="",string="",string="",string="",string="",
			   string="",string="",string="",string="",string="",string="",string="",string="",
			   string="",string="",string="",string="",string="",string="",string="",string="");
  uint stream2dequestring(std::istream& istreamIN,deque<string> &vstringout);
  uint stream2dequestring(std::ifstream& ifstreamIN,deque<string> &vstringout);
  uint stream2dequestring(stringstream& stringstreamIN,deque<string> &vstringout);
  uint string2dequestring(const string& stringIN,deque<string> &vstringout);
  deque<string> stream2dequestring(std::istream& istreamIN);
  deque<string> stream2dequestring(std::ifstream& ifstreamIN);
  deque<string> stream2dequestring(stringstream& stringstreamIN);
  deque<string> string2dequestring(const string& stringIN);
  // about writing files
  bool string2file(const string& StringOUTPUT,const string& FileNameOUTPUT,string="");
  bool string2compressfile(const string& commamnd,const string& StringOUTPUT,const string& FileNameOUTPUT,string="");
  bool string2gzfile(const string& StringOUTPUT,const string& FileNameOUTPUT,string="");
  bool string2bz2file(const string& StringOUTPUT,const string& FileNameOUTPUT,string="");
  bool string2xzfile(const string& StringOUTPUT,const string& FileNameOUTPUT,string="");
  bool stringstream2file(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,string="");
  bool stringstream2compressfile(const string& commamnd,const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,string="");
  bool stringstream2gzfile(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,string="");
  bool stringstream2bz2file(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,string="");
  bool stringstream2xzfile(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,string="");
  // file2string
  uint file2string(string FileNameIN,string& StringIN);
  uint bz2file2string(string FileNameIN,string& StringIN);
  uint gzfile2string(string FileNameIN,string& StringIN);
  uint xzfile2string(string FileNameIN,string& StringIN);
  uint zipfile2string(string _FileNameIN,string& StringIN); //CO
  uint efile2string(string FileNameIN,string& StringIN);
  // file2vectorstring  
  uint file2vectorstring(string FileNameIN,vector<string>& vlines);
  uint bz2file2vectorstring(string FileNameIN,vector<string>& vlines);
  uint gzfile2vectorstring(string FileNameIN,vector<string>& vlines);
  uint xzfile2vectorstring(string FileNameIN,vector<string>& vlines);
  uint efile2vectorstring(string FileNameIN,vector<string>& vlines);
  // file2dequestring  
  uint file2dequestring(string FileNameIN,deque<string>& vlines);
  uint bz2file2dequestring(string FileNameIN,deque<string>& vlines);
  uint gzfile2dequestring(string FileNameIN,deque<string>& vlines);
  uint xzfile2dequestring(string FileNameIN,deque<string>& vlines);
  uint efile2dequestring(string FileNameIN,deque<string>& vlines);
  // file2stringstream
  bool file2stringstream(string FileNameIN,stringstream& StringstreamIN);
  bool bz2file2stringstream(string FileNameIN,stringstream& StringstreamIN);
  bool gzfile2stringstream(string FileNameIN,stringstream& StringstreamIN);
  bool xzfile2stringstream(string FileNameIN,stringstream& StringstreamIN);
  bool zipfile2stringstream(string _FileNameIN,stringstream& StringstreamIN); //CO
  bool efile2stringstream(string FileNameIN,stringstream& StringstreamIN);
  // return directly the string
  string file2string(string FileNameIN);
  string bz2file2string(string FileNameIN);
  string gzfile2string(string FileNameIN);
  string xzfile2string(string FileNameIN);
  string efile2string(string FileNameIN);
  // reading url to string/stringstream/tokens/vector/deque
  bool url2file(string url,string& fileIN,bool=FALSE)  __xprototype;   // bool = verbose
  bool url2string(string url,string& stringIN,bool=FALSE)  __xprototype;   // bool = verbose
  bool url2stringstream(string url,stringstream& stringstreamIN,bool=FALSE)  __xprototype;  // bool = verbose
  bool url2vectorstring(string url,vector<string>& vlines,bool=FALSE)  __xprototype;  // bool = verbose
  bool url2dequestring(string url,deque<string>& vlines,bool=FALSE)  __xprototype;  // bool = verbose
  template<typename utype> uint url2tokens(string url,vector<utype>& tokens,const string& delimiters = " ")  __xprototype;
  template<typename utype> uint url2tokens(string url,deque<utype>& tokens,const string& delimiters = " ")  __xprototype;
  string url2string(string url)  __xprototype;
  // about istream/ostream
  string ostream2string(std::ostream& oss);
  uint stream2string(std::istream& istreamIN,string &vstringout);
  uint stream2string(std::ifstream& ifstreamIN,string &vstringout);
  uint stream2string(stringstream& stringstreamIN,string &vstringout);
  // about environments
  string getenv2string(const string& str);
  int getenv2int(const string& str);
  uint getenv2uint(const string& str);
  double getenv2double(const string& str);
  // about operating on files
  bool ChmodFile(string chmod_string,string FileNameOUTPUT);
  //CO - START
  bool MoveFile(const string& _file,const string& destination);
  bool MoveFiles(const vector<string>& files,const string& destination);
  bool RenameFile(const string& _file,const string& destination); //CO 171025
  bool IsDirectory(string path);
  bool IsFile(string path);
  //CO - END
  bool RemoveFile(string FileNameOUTPUT);
  bool RemoveFiles(const vector<string>& files); //CO
  bool RemoveDirectory(const string& dir); //CO
  // about getting info from strings
  uint string2tokens(const string& str,vector<string>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype;  //CO 170613, defaults to usual string2tokens() behavior
  uint string2tokens(const string& str,deque<string>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype; //CO 170613, defaults to usual string2tokens() behavior
  template<class utype> uint string2tokens(const string& str,std::vector<utype>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype;  //CO 170613, defaults to usual string2tokens() behavior
  template<class utype> uint string2tokens(const string& str,std::deque<utype>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype; //CO 170613, defaults to usual string2tokens() behavior
  uint string2tokensAdd(const string& str,vector<string>& tokens,const string& delimiters = " ") __xprototype;
  uint string2tokensAdd(const string& str,deque<string>& tokens,const string& delimiters = " ") __xprototype;
  template<class utype> uint string2tokensAdd(const string& str,std::vector<utype>& tokens,const string& delimiters = " ") __xprototype;
  template<class utype> uint string2tokensAdd(const string& str,std::deque<utype>& tokens,const string& delimiters = " ") __xprototype;

  template<typename typeTo, typename typeFrom> typeTo StringStreamConvert(typeFrom from) __xprototype;
  //  template<typename typeFrom> string StringStreamConvert(typeFrom from) __xprototype;
  template<typename typeFrom> string StringConvert(typeFrom from) __xprototype;
  template<typename typeTo, typename typeFrom> typeTo NumberStreamConvert(typeFrom from) __xprototype;

  // [OBSOLETE]  double string2double(const string& from) __xprototype;
  vector<double> vectorstring2vectordouble(vector<string> from) __xprototype;
  // [OBSOLETE]  long double string2longdouble(const string& from) __xprototype;
  // [OBSOLETE]  int string2int(const string& from) __xprototype;
  string string2string(const string& from) __xprototype;
  template<typename utype> utype string2utype(const string& from) __xprototype;
  vector<int> vectorstring2vectorint(vector<string> from) __xprototype;
  // [OBSOLETE] uint string2uint(const string& from) __xprototype;
  vector<uint> vectorstring2vectoruint(vector<string> from) __xprototype;
  // [OBSOLETE] long int string2longint(const string& from) __xprototype;
  // [OBSOLETE] float string2float(const string& from) __xprototype;

  vector<float> vectorstring2vectorfloat(vector<string> from) __xprototype;
  string vectorstring2string(const vector<string>& vstrings) __xprototype;
  string vectorstring2string(const deque<string>& vstrings) __xprototype;

  // [OBSOLETE] string double2string(double from) __xprototype;
  // [OBSOLETE] string double2string(double from,int precision) __xprototype;
  // [OBSOLETE] string longdouble2string(long double from) __xprototype;
  // [OBSOLETE] string int2string(int from) __xprototype;
  // [OBSOLETE] string uint2string(uint from) __xprototype;
  // [OBSOLETE] string longint2string(long int from) __xprototype;
  // [OBSOLETE] string float2string(float from) __xprototype;
  template<typename utype> string utype2string(const utype& from) __xprototype;
  template<typename utype> string utype2string(const utype& from,int precision) __xprototype;
  string utype2string(double from,bool roff);
  string utype2string(double from,int precision,bool roff);
  string utype2string(double from,bool roff,double tol);
  string utype2string(double from,int precision,bool roff,double tol);
  string utype2string(double from,bool roff,char FORMAT);
  string utype2string(double from,int precision,bool roff,char FORMAT);
  string utype2string(double from,bool roff,double tol,char FORMAT);
  string utype2string(double from,int precision,bool roff,double tol,char FORMAT);
  // [OBSOLETE]  template<class u1> string utype2string(u1) __xprototype;
  // [OBSOLETE]  template<class u1, class u2> string utype2string(u1,u2) __xprototype;
  // [OBSOLETE]  template<class u1, class u2, class u3> string utype2string(u1,u2,u3) __xprototype;
  // [OBSOLETE]  template<class u1, class u2, class u3, class u4> string utype2string(u1,u2,u3,u4) __xprototype;
  
  template<class utype> deque<utype> utypes2deque(utype u1) __xprototype;
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2) __xprototype;
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3) __xprototype;
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3,utype u4) __xprototype;

  void StringCommasColumsVectorInt(string vstring,vector<int> &vint) __xprototype;
  void StringCommasColumsVectorUnsignedInt(string vstring,vector<uint> &vuint) __xprototype;
  void StringCommasColumsVectorFloat(string vstring,vector<float> &vfloat) __xprototype;
  void StringCommasColumsVectorDouble(string vstring,vector<double> &vdouble) __xprototype;
  int GetNLinesString(const string& strstream) __xprototype;
  int GetNLinesString(const stringstream& strstream) __xprototype;
  int GetNLinesFile(const string& file_name) __xprototype;
  string GetLineString(const string& strstream, const int& line) __xprototype;
  string GetLineString(const stringstream& strstream, const int& line) __xprototype;
  stringstream GetLineStringstream(const stringstream& strstream, const int& line) __xprototype;
  // substitute strings in strings and stringstreams
  bool StringsAlphabetic(string A,string B);
  string StringSubst(string &strstring, const string &strfind, const string &strreplace);
//  string StringSubst(string &strstring, const string &strfind0, const string &strfind1, const string &strfind2, const string &strfind3, const string &strreplace);
  string StringSubst(string &strstring, const char &charfind, const char &charreplace);
  void StringStreamSubst(string &strstring, const string &strfind, const string &strreplace);
  // about present substrings
  bool substring2bool(const string& strstream, const string& strsub1, bool CLEAN);
  bool substring2bool(const vector<string>& vstrstream, const string& strsub1, bool CLEAN);
  bool substring2bool(const deque<string>& vstrstream, const string& strsub1, bool CLEAN);
  bool substring2bool(const stringstream& strstream, const string& strsub1, bool CLEAN);
  bool substring_present_file(const string& FileName, const string& strsub1, bool CLEAN);
  bool substring_present_file_FAST(const string& FileName, const string& strsub1, bool CLEAN);
  bool substring2bool(const string& strstream, const string& strsub1);
  bool substring2bool(const vector<string>& vstrstream, const string& strsub1);
  bool substring2bool(const deque<string>& vstrstream, const string& strsub1);
  bool substring2bool(const stringstream& strstream, const string& strsub1);
  bool withinList(const vector<string>& list,const string& str);
  bool substring_present_file(const string& FileName, const string& strsub1) ;
  bool substring_present_file_FAST(const string& FileName, const string& strsub1) ;
  // about present substrings and taking off the value
  string substring2string(const string& strstream, const string& strsub1) __xprototype;
  string substring2string(const string& strstream, const string& strsub1, bool CLEAN) __xprototype;
  string substring2string(const string& strstream, const string& strsub1, const string& strsub2) __xprototype;
  string substring2string(const string& strstream, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;

  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1) __xprototype;
  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1, bool CLEAN) __xprototype;
  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1, const string& strsub2) __xprototype;
  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  // [OBSOLETE] int substring2integer(const string& strstream, const string& strsub1) __xprototype;
  // [OBSOLETE] int substring2integer(const string& strstream, const string& strsub1, bool CLEAN) __xprototype;
  // [OBSOLETE] int substring2integer(const string& strstream, const string& strsub1, const string& strsub2) __xprototype;
  // [OBSOLETE] int substring2integer(const string& strstream, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  // [OBSOLETE] int substring2int(const string& strstream, const string& strsub1) __xprototype;
  // [OBSOLETE] int substring2int(const string& strstream, const string& strsub1, bool CLEAN) __xprototype;
  // [OBSOLETE] int substring2int(const string& strstream, const string& strsub1, const string& strsub2) __xprototype;
  // [OBSOLETE] int substring2int(const string& strstream, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  // [OBSOLETE] double substring2double(const string& strstream, const string& strsub1) __xprototype;
  // [OBSOLETE] double substring2double(const string& strstream, const string& strsub1, bool CLEAN) __xprototype;
  // [OBSOLETE] double substring2double(const string& strstream, const string& strsub1, const string& strsub2) __xprototype;
  // [OBSOLETE] double substring2double(const string& strstream, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  // about present substrings and taking off the values in vectors

  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1) __xprototype;
  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, bool CLEAN) __xprototype;
  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2) __xprototype;
  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  template<typename utype> uint substring2utypes(const string& strstream, vector<string> &vstringout, const string& strsub1) __xprototype;
  template<typename utype> uint substring2utypes(const string& strstream, vector<string> &vstringout, const string& strsub1, bool CLEAN) __xprototype;
  template<typename utype> uint substring2utypes(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  template<typename utype> uint substring2utypes(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2) __xprototype;
  // [OBSOLETE] uint substring2integers(const string& strstream, vector<int> &vintout, const string& strsub1) __xprototype;
  // [OBSOLETE] uint substring2integers(const string& strstream, vector<int> &vintout, const string& strsub1, bool CLEAN) __xprototype;
  // [OBSOLETE] uint substring2integers(const string& strstream, vector<int> &vintout, const string& strsub1, const string& strsub2) __xprototype;
  // [OBSOLETE] uint substring2integers(const string& strstream, vector<int> &vintout, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  // [OBSOLETE] uint substring2ints(const string& strstream, vector<int> &vintout, const string& strsub1) __xprototype;
  // [OBSOLETE] uint substring2ints(const string& strstream, vector<int> &vintout, const string& strsub1, bool CLEAN) __xprototype;
  // [OBSOLETE] uint substring2ints(const string& strstream, vector<int> &vintout, const string& strsub1, const string& strsub2) __xprototype;
  // [OBSOLETE] uint substring2ints(const string& strstream, vector<int> &vintout, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;
  // [OBSOLETE] uint substring2doubles(const string& strstream, vector<double> &vdoubleout, const string& strsub1) __xprototype;
  // [OBSOLETE] uint substring2doubles(const string& strstream, vector<double> &vdoubleout, const string& strsub1, bool CLEAN) __xprototype;
  // [OBSOLETE] uint substring2doubles(const string& strstream, vector<double> &vdoubleout, const string& strsub1, const string& strsub2) __xprototype;
  // [OBSOLETE] uint substring2doubles(const string& strstream, vector<double> &vdoubleout, const string& strsub1, const string& strsub2, bool CLEAN) __xprototype;

}

// ***************************************************************************
namespace aurostd { // aurostd_crc64.cpp
  uint64_t crc64(uint64_t crc, const unsigned char *s, uint64_t l);
  uint64_t crc64(uint64_t crc, const string s);
  string crc2string(uint64_t crc);
  int crc64_main(void);
}

// ***************************************************************************

namespace aurostd {
  string html2latex(const string& str) __xprototype;
  string html2txt(const string& str) __xprototype;
  string string2latex(const string& str) __xprototype;
  string latex2html(const string& str) __xprototype;
  string latex2txt(const string& str) __xprototype;
}


// ***************************************************************************
// SORT WORLD
// double
// ----------------------------------------------------------------------------
// sort for vectors

namespace aurostd {
  template<class utype1> void sort(vector<utype1>& arr);
  template<class utype1> void sort_remove_duplicates(vector<utype1>& arr);
  template<class utype1,class utype2> void sort(vector<utype1>& arr, vector<utype2>& brr);
  template<class utype1,class utype2,class utype3> void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr);
  template<class utype1,class utype2,class utype3,class utype4> void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr, vector<utype4>& drr);
}

namespace aurostd {   // DOUBLE
  class _sort_double_value0 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]<v2[0]); } };
  class _isort_double_value0 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]>v2[0]); } };
  class _sort_double_value1 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[1]<v2[1]); } };
  class _isort_double_value1 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[1]>v2[1]); } };
  class _sort_double_value2 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[2]<v2[2]); } };
  class _isort_double_value2 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[2]>v2[2]); } };
  class _sort_double_value3 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[3]<v2[3]); } };
  class _isort_double_value3 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[3]>v2[3]); } };
  class _sort_double_value4 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[4]<v2[4]); } };
  class _isort_double_value4 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[4]>v2[4]); } };
  class _sort_double_value5 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[5]<v2[5]); } };
  class _isort_double_value5 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[5]>v2[5]); } };
  class _sort_double_value6 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[6]<v2[6]); } };
  class _isort_double_value6 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[6]>v2[6]); } };
  class _sort_double_value7 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[7]<v2[7]); } };
  class _isort_double_value7 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[7]>v2[7]); } };
  class _sort_double_value8 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[8]<v2[8]); } };
  class _isort_double_value8 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[8]>v2[8]); } };
  class _sort_double_value9 {                    // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[9]<v2[9]); } };
  class _isort_double_value9 {                   // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[9]>v2[9]); } };
  class _sort_double_value01 {                // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]<v2[0]+v2[1]); } };
  class _isort_double_value01 {               // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]>v2[0]+v2[1]); } };
  class _sort_double_value012 {                // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]<v2[0]+v2[1]+v2[2]); } };
  class _isort_double_value012 {               // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]>v2[0]+v2[1]+v2[2]); } };
  class _sort_double_value0123 {                // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]<v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _isort_double_value0123 {               // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]>v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _sort_double_value01234 {                // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]<v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
  class _isort_double_value01234 {               // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]>v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
}
// int
namespace aurostd {   // INT
  class _sort_int_value0 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]<v2[0]); } };
  class _isort_int_value0 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]>v2[0]); } };
  class _sort_int_value1 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[1]<v2[1]); } };
  class _isort_int_value1 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[1]>v2[1]); } };
  class _sort_int_value2 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[2]<v2[2]); } };
  class _isort_int_value2 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[2]>v2[2]); } };
  class _sort_int_value3 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[3]<v2[3]); } };
  class _isort_int_value3 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[3]>v2[3]); } };
  class _sort_int_value4 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[4]<v2[4]); } };
  class _isort_int_value4 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[4]>v2[4]); } };
  class _sort_int_value5 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[5]<v2[5]); } };
  class _isort_int_value5 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[5]>v2[5]); } };
  class _sort_int_value6 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[6]<v2[6]); } };
  class _isort_int_value6 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[6]>v2[6]); } };
  class _sort_int_value7 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[7]<v2[7]); } };
  class _isort_int_value7 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[7]>v2[7]); } };
  class _sort_int_value8 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[8]<v2[8]); } };
  class _isort_int_value8 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[8]>v2[8]); } };
  class _sort_int_value9 {                    // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[9]<v2[9]); } };
  class _isort_int_value9 {                   // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[9]>v2[9]); } };
  class _sort_int_value01 {                // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]<v2[0]+v2[1]); } };
  class _isort_int_value01 {               // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]>v2[0]+v2[1]); } };
  class _sort_int_value012 {                // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]<v2[0]+v2[1]+v2[2]); } };
  class _isort_int_value012 {               // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]>v2[0]+v2[1]+v2[2]); } };
  class _sort_int_value0123 {                // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]<v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _isort_int_value0123 {               // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]>v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _sort_int_value01234 {                // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]<v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
  class _isort_int_value01234 {               // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]>v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
  // STRING
  void sort(vector<string>& arg);
  void sort(deque<string>& arg);
  void sort_remove_duplicates(vector<string>& arg);
  void sort_remove_duplicates(deque<string>& arg);
  class _sort_string_ {                    // sorting through reference
  public:
    bool operator()(const string& str1, const string& str2) const { return (bool) (str1<str2);}
  };
  void rsort(vector<string>& arg);
  void rsort(deque<string>& arg);
  void rsort_remove_duplicates(vector<string>& arg);
  void rsort_remove_duplicates(deque<string>& arg);

  
  // _STRING_INT_
  void sort(vector<string>& varg1,vector<int>& varg2);
  void sort(deque<string>& varg1,deque<int>& varg2);
  struct _string_int_ {
    string arg1;
    int arg2;
  };
  class _sort_string_int_ {                    // sorting through reference
  public:
    bool operator()(const _string_int_& x1, const _string_int_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_DOUBLE_
  void sort(vector<string>& varg1,vector<double>& varg2);
  void sort(deque<string>& varg1,deque<double>& varg2);
  struct _string_double_ {
    string arg1;
    double arg2;
  };
  class _sort_string_double_ {                    // sorting through reference
  public:
    bool operator()(const _string_double_& x1, const _string_double_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_
  void sort(vector<string>& varg1,vector<string>& varg2);
  void sort(deque<string>& varg1,deque<string>& varg2);
  struct _string_string_ {
    string arg1;
    string arg2;
  };
  class _sort_string_string_ {                    // sorting through reference
  public:
    bool operator()(const _string_string_& x1, const _string_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _DOUBLE_INT_
  void sort(vector<double>& varg1,vector<int>& varg2);
  void sort(deque<double>& varg1,deque<int>& varg2);
  struct _double_int_ {
    double arg1;
    int arg2;
  };
  class _sort_double_int_ {                    // sorting through reference
  public:
    bool operator()(const _double_int_& x1, const _double_int_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _DOUBLE_DOUBLE_
  void sort(vector<double>& varg1,vector<double>& varg2);
  void sort(deque<double>& varg1,deque<double>& varg2);
  struct _double_double_ {
    double arg1;
    double arg2;
  };
  class _sort_double_double_ {                    // sorting through reference
  public:
    bool operator()(const _double_double_& x1, const _double_double_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _DOUBLE_STRING_
  void sort(vector<double>& varg1,vector<string>& varg2);
  void sort(deque<double>& varg1,deque<string>& varg2);
  struct _double_string_ {
    double arg1;
    string arg2;
  };
  class _sort_double_string_ {                    // sorting through reference
  public:
    bool operator()(const _double_string_& x1, const _double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_INT_STRING
  void sort(vector<string>& varg1,vector<int>& varg2,vector<string>& varg3);
  void sort(deque<string>& varg1,deque<int>& varg2,deque<string>& varg3);
  struct _string_int_string_ {
    string arg1;
    int arg2;
    string arg3;
  };
  class _sort_string_int_string_ {             // sorting through reference
  public:
    bool operator()(const _string_int_string_& x1, const _string_int_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_DOUBLE_STRING
  void sort(vector<string>& varg1,vector<double>& varg2,vector<string>& varg3);
  void sort(deque<string>& varg1,deque<double>& varg2,deque<string>& varg3);
  struct _string_double_string_ {
    string arg1;
    double arg2;
    string arg3;
  };
  class _sort_string_double_string_ {             // sorting through reference
  public:
    bool operator()(const _string_double_string_& x1, const _string_double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_STRING
  void sort(vector<string>& varg1,vector<string>& varg2,vector<string>& varg3);
  void sort(deque<string>& varg1,deque<string>& varg2,deque<string>& varg3);
  struct _string_string_string_ {
    string arg1;
    string arg2;
    string arg3;
  };
  class _sort_string_string_string_ {             // sorting through reference
  public:
    bool operator()(const _string_string_string_& x1, const _string_string_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_DOUBLE_STRING
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<string>& varg4);
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<string>& varg4);
  struct _string_string_double_string_ {
    string arg1;
    string arg2;
    double arg3;
    string arg4;
  };
  class _sort_string_string_double_string_ {             // sorting through reference
  public:
    bool operator()(const _string_string_double_string_& x1, const _string_string_double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_DOUBLE_DOUBLE_STRING
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<double>& varg4,vector<string>& varg5);
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<double>& varg4,deque<string>& varg5);
  struct _string_string_double_double_string_ {
    string arg1;
    string arg2;
    double arg3;
    double arg4;
    string arg5;
  };
  class _sort_string_string_double_double_string_ {             // sorting through reference
  public:
    bool operator()(const _string_string_double_double_string_& x1, const _string_string_double_double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
}

// ***************************************************************************
// some statistical stuff
template<class utype> utype combinations(utype n,utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination
template<class utype> utype Cnk(utype n,utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination

// ***************************************************************************
namespace aurostd {
  /* template <typename utype> utype sum(vector<utype>& a) {
     utype result = 0;
     for (unsigned int i=0; i<a.size();i++) result += a.at(i);
     return result;
     } */
  vector<vector<double> > ShiftFirstColumn(const vector<vector<double> >& vva, const double& value);
  vector<vector<double> > ShrinkValuesExceptFirstColumn(const vector<vector<double> >& vva, const double& Fi);
  vector<vector<double> > NormalizeAndSum3DVector(const vector<vector<vector<double> > >& vvva, const vector<double>& vFi);
  vector<vector<double> > Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva);
  vector<vector<double> > Sum2DVectorExceptFirstColumn(const vector<vector<double> >& vva, const vector<vector<double> >& vvb);
  string vector2string(const vector<vector<double> >& vva);
  vector<vector<double> > ReduceVector(const vector<vector<double> >& vva, const int& n);
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n);
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n, const double& Emin, const double& Emax);
  double CalculateIntegrate(const vector<vector<double> >& vva);
  double CalculateIntegrate(const vector<vector<double> >& vva, const double& Emin, const double& Emax);
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva);
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva, const double& min, const double& max);
  double FindMaxInTDOS(const vector<vector<double> >& vva, const double& min, const double& max);
}

// ***************************************************************************


// ***************************************************************************
bool initialize_templates_never_call_this_procedure(bool flag);

////////////////////////////////////////////////////////////////////////////////
namespace aurostd {
  // joins int/string type of objects together by a delimiter
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter,
			const char& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter);
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter,
			const char& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter);
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter);
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter);
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter);
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter);
}
////////////////////////////////////////////////////////////////////////////////

// DX 1/18/18 - Add xcomplex to json
namespace aurostd {
 template<typename utype> string _xcomplex2json(xcomplex<utype>& number) __xprototype;
 string xcomplex2json(xcomplex<double>& number);
}

// DX 8/3/17 - Add Matrix print out
namespace aurostd {
  // [OBSOLETE] string xmatDouble2String(const xmatrix<double>& xmat_in, bool roff=false);
  string xmatDouble2String(const xmatrix<double>& xmat_in, int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
}

//CO 171215 - more json functionality
namespace aurostd {
  string wrapString(const string& input,const string& wrapper);
  string wrapString(const string& input,const string& wrapper_start,const string& wrapper_end);
}

namespace aurostd {
  // [OBSOLETE] vector<string> vecDouble2vecString(const vector<double>& vin, bool roff=false);
  vector<string> vecDouble2vecString(const vector<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
  // [OBSOLETE] vector<string> xvecDouble2vecString(const xvector<double>& vin, bool roff=false);
  vector<string> xvecDouble2vecString(const xvector<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
  // [OBSOLETE] deque<string> deqDouble2deqString(const deque<double>& vin, bool roff=false);
  deque<string> deqDouble2deqString(const deque<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
}

namespace aurostd {
  vector<string> wrapVecEntries(const vector<string>& vin,string wrap);
  vector<string> wrapVecEntries(const vector<string>& vin,string wrap_start,string wrap_end);
  deque<string> wrapDeqEntries(const deque<string>& vin,string wrap);
  deque<string> wrapDeqEntries(const deque<string>& vin,string wrap_start,string wrap_end);
}

//base64 stuff
//CO - START
namespace aurostd {
  //http://www.adp-gmbh.ch/cpp/common/base64.html for base64 encoding/decoding
  //http://stackoverflow.com/questions/535444/custom-manipulator-for-c-iostream for fancy stream manipulators
  //static inline bool isBase64(unsigned char c); //determines if char is base64
  inline bool isBase64(unsigned char c); //determines if char is base64
  std::string base64Encoder(unsigned char const* bytes_to_encode, unsigned int in_len); //encodes bytes to base64
  std::string base64Decoder(std::string const& encoded_string); //decodes base64 to bytes
  bool bin2base64(const std::string& b_file, std::string& b64String); //converts binary file to base64 string
  bool base642bin(const std::string& b64String, const std::string& b_file); //converts base64 string to binary file

  //structure that allows "cout << b64_encoder << ifstream"
  struct b64_encoder_proxy {
    explicit b64_encoder_proxy(std::ostream & os):os(os){}

    template<typename Rhs>
    friend std::ostream & operator<<(b64_encoder_proxy const& b, 
                                     Rhs const& rhs) {
      return b.os << rhs;
    }

    friend std::ostream & operator<<(b64_encoder_proxy const& b, 
                                     ifstream& file) {
      // Stop eating new lines in binary mode!!!
      file.unsetf(std::ios::skipws);

      // get its size:
      std::streampos fileSize;

      file.seekg(0, std::ios::end);
      fileSize = file.tellg();
      file.seekg(0, std::ios::beg);

      // read the data:
      vector<unsigned char> vec((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));

      return b.os << base64Encoder(reinterpret_cast<const unsigned char*>(&vec[0]), vec.size());
    }
  
  private:
    std::ostream & os;
  };

  //structure that allows "ofstream << b64_decoder << string" 
  //WARNING, SPITS BINARY CHARACTERS OUT IF COUT
  struct b64_decoder_proxy {
    explicit b64_decoder_proxy(std::ostream & os):os(os){}

    template<typename Rhs>
    friend std::ostream & operator<<(b64_decoder_proxy const& b, 
                                     Rhs const& rhs) {
      return b.os << rhs;
    }

    friend std::ostream & operator<<(b64_decoder_proxy const& b, 
                                     std::string const& rhs) {
      return b.os << base64Decoder(rhs);
    }

  private:
    std::ostream & os;
  };
  
  //structure that allows "cout << b64_encoder << ifstream"
  struct b64_encoder_creator { };
  extern b64_encoder_creator b64_encoder;
  b64_encoder_proxy operator<<(std::ostream & os, b64_encoder_creator);
  
  //structure that allows "ofstream << b64_decoder << string" 
  //WARNING, SPITS BINARY CHARACTERS OUT IF COUT
  struct b64_decoder_creator { };
  extern b64_decoder_creator b64_decoder;
  b64_decoder_proxy operator<<(std::ostream & os, b64_decoder_creator);
}

//binary to base64 conversion
const std::string base64_chars = 
  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  "abcdefghijklmnopqrstuvwxyz"
  "0123456789+/";
//CO - END

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
