// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AUROSTD_MAIN_CPP_
#define _AUROSTD_MAIN_CPP_
//#include "aflow.h"
#include "aurostd.h"

#define _CIN_LINE_BUFFER_LENGTH_     16384
#ifndef  CHMOD_BIN
#define  CHMOD_BIN  string("chmod")
#endif

using std::vector;   // for pennsy
using std::deque;   // for pennsy
using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;
using std::string;
using std::cerr;

using aurostd::utype2string;
using aurostd::xvector;
using aurostd::xmatrix;
using aurostd::sign;
using aurostd::ran0;

#define COMMENT_NEGLECT_1 string("#")
//#define COMMENT_NEGLECT_2 string("// ")
#define COMMENT_NEGLECT_2 string("//")
#define COMMENT_NEGLECT_3 string("!")

//CO 171215 - moved to xscalar
// ***************************************************************************
// ROUNDOFF for scalars
//namespace aurostd { //DX add roundoff for scalar values
//  template<class utype>
//  utype roundoff(const utype& x, utype _tol_){
//    return ((abs(x)<(utype) _tol_) ? (utype) 0.0 : x);
//  }
//  double _aurostd_initialize_roundoff(const double& x,double y) {return roundoff(x,y);}
//  float _aurostd_initialize_roundoff(const float& x,float y) {return roundoff(x,y);}
//  int _aurostd_initialize_roundoff(const int& x,int y) {return roundoff(x,y);}
//}

// ***************************************************************************
// TIME evolution stuff
namespace aurostd {
  int get_day(void) { time_t t=time(0);struct tm *now=localtime(&t);return now->tm_mday;}
  int get_month(void) { time_t t=time(0);struct tm *now=localtime(&t);return now->tm_mon+1;}
  int get_year(void) { time_t t=time(0);struct tm *now=localtime(&t);return now->tm_year+1900;}
  long int get_date(void) { return aurostd::get_year()*10000+aurostd::get_month()*100+aurostd::get_day();}
  int get_hour(void) { time_t t=time(0);struct tm *now=localtime(&t);return now->tm_hour;}
  int get_min(void) { time_t t=time(0);struct tm *now=localtime(&t);return now->tm_min;}
  int get_sec(void) { time_t t=time(0);struct tm *now=localtime(&t);return now->tm_sec;}
  long double get_seconds(void) {timeval tim;gettimeofday(&tim,NULL);return tim.tv_sec+tim.tv_usec/1e6;}
  long double get_seconds(long double reference_seconds) { return get_seconds()-reference_seconds;}
  long double get_delta_seconds(long double& seconds_begin) {long double out=get_seconds()-seconds_begin;seconds_begin=get_seconds();return out;}
  long double get_useconds(void) {timeval tim;gettimeofday(&tim,NULL);return tim.tv_usec;}
  long double get_useconds(long double reference_useconds) { return get_useconds()-reference_useconds;}
  long double get_delta_useconds(long double& useconds_begin) {long double out=get_useconds()-useconds_begin;useconds_begin=get_useconds();return out;}
  string get_time(void) {int h=get_hour(),m=get_min(),s=get_sec();return (h<10?"0":"")+aurostd::utype2string(h)+":"+(m<10?"0":"")+aurostd::utype2string(m)+":"+(s<10?"0":"")+aurostd::utype2string(s);}
  string get_datetime(void) { return utype2string(get_date())+"_"+get_time();}
  string get_datetime_formatted(bool separate_date_time){  //CO 171215
    stringstream misc_ss;
    int y=aurostd::get_year(),b=aurostd::get_month(),d=aurostd::get_day(),h=get_hour(),m=get_min(),s=get_sec();
    misc_ss << y << "/" << (b<10?"0":"") << b << "/" << (d<10?"0":"") << d << (separate_date_time?" ":"-") << (h<10?"0":"") << h << ":" << (m<10?"0":"") << m << ":" << (s<10?"0":"") << s;
    return misc_ss.str();
  }

  bool beep(uint freq,uint duration) {
    return aurostd::execute("beep -f "+aurostd::utype2string<uint>(freq)+" -l "+aurostd::utype2string<uint>(duration));
  }
}

// ***************************************************************************
// FILES creation/destruction
namespace aurostd {
  string TmpFileCreate(string identifier) {
    string str=XHOST.Tmpfs+"/_aflow_"+identifier+"."+XHOST.User+".pid"+XHOST.ostrPID.str()+".a"+AFLOW_VERSION+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) get_useconds()))+".tmp";
    // cerr << str << endl;
    return str;
  }
  string TmpFileCreate(void) {
    return TmpFileCreate("");}
  string TmpDirectoryCreate(string identifier) {
    string dir=XHOST.Tmpfs+"/_aflow_"+identifier+"_"+XHOST.User+"_pid"+XHOST.ostrPID.str()+"_a"+AFLOW_VERSION+"_rnd"+aurostd::utype2string(uint((double) std::floor((double) 100000*aurostd::ran0())))+"_u"+aurostd::utype2string(uint((double) get_useconds()))+"_tmp";
    DirectoryMake(dir);
    return dir;}
  string TmpDirectoryCreate(void) {
    return TmpDirectoryCreate("");}
}

// ***************************************************************************
// Function extra operator << for vector
template<class utype>                            // operator <<  vector<>
std::ostream& operator<< (std::ostream& buf,const std::vector<utype>& x) {
  for(uint i=0;i<x.size();i++) {
    buf << x.at(i) << " ";
  }
  return buf;
}
// ***************************************************************************
// Function extra operator << for deque
template<class utype>                            // operator <<  deque<>
std::ostream& operator<< (std::ostream& buf,const std::deque<utype>& x) {
  for(uint i=0;i<x.size();i++) {
    buf << x.at(i) << " ";
  }
  return buf;
}
// ***************************************************************************
ostream& operator<< (ostream& b,const vector<uint>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const deque<uint>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const vector<char>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const deque<char>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const vector<int>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const deque<int>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const vector<long>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const deque<long>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const vector<double>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const deque<double>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const vector<long double>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const deque<long double>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const vector<string>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}
ostream& operator<< (ostream& b,const deque<string>& x) {for(uint i=0;i<x.size();i++) b << x.at(i) << " "; return b;}

namespace aurostd {
  // ***************************************************************************
  // Function aswap
  // ***************************************************************************
  // namespace aurostd
  template<class utype> void aswap(utype &a,utype &b) {utype temp=a;a=b;b=temp;}
  void _aurostd_initialize_aswap(bool& x,bool& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(char& x,char& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(int& x,int& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(uint& x,uint& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(float& x,float& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(double& x,double& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(string& x,string& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(long int& x,long int& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(long long int& x,long long int& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(long double& x,long double& y) {aswap(x,y);}
#ifdef _AUROSTD_XCOMPLEX_
  //   void _aurostd_initialize_aswap(xcomplex<float>& x,xcomplex<float>& y) {aswap(x,y);}
  //   void _aurostd_initialize_aswap(xcomplex<double>& x,xcomplex<double>& y) {aswap(x,y);}
  //   void _aurostd_initialize_aswap(xcomplex<long double>& x,xcomplex<long double>& y) {aswap(x,y);}
#endif

  // ***************************************************************************
  // Function max of a vector/deque
  // ***************************************************************************
  // Stefano
  template<class utype> utype max(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec.at(i)>=out) out=vec.at(i);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_max(const vector<bool> vec) { return max(vec);}
  char _aurostd_initialize_max(const vector<char> vec) { return max(vec);}
  string _aurostd_initialize_max(const vector<string> vec) { return max(vec);}
  int _aurostd_initialize_max(const vector<int> vec) { return max(vec);}
  long _aurostd_initialize_max(const vector<long> vec) { return max(vec);}
  uint _aurostd_initialize_max(const vector<uint> vec) { return max(vec);}
  float _aurostd_initialize_max(const vector<float> vec) { return max(vec);}
  double _aurostd_initialize_max(const vector<double> vec) { return max(vec);}
  long double _aurostd_initialize_max(const vector<long double> vec) { return max(vec);}

  template<class utype> utype max(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec.at(i)>=out) out=vec.at(i);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_max(const deque<bool> vec) { return max(vec);}
  char _aurostd_initialize_max(const deque<char> vec) { return max(vec);}
  string _aurostd_initialize_max(const deque<string> vec) { return max(vec);}
  int _aurostd_initialize_max(const deque<int> vec) { return max(vec);}
  long _aurostd_initialize_max(const deque<long> vec) { return max(vec);}
  uint _aurostd_initialize_max(const deque<uint> vec) { return max(vec);}
  float _aurostd_initialize_max(const deque<float> vec) { return max(vec);}
  double _aurostd_initialize_max(const deque<double> vec) { return max(vec);}
  long double _aurostd_initialize_max(const deque<long double> vec) { return max(vec);}

  // ***************************************************************************
  // Function max of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype max(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=mat.at(0).at(0);
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat.at(i).size();j++)
	if(mat.at(i).at(j)>=out) out=mat.at(i).at(j);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_max(const vector<vector<bool> > mat) { return max(mat);}
  char _aurostd_initialize_max(const vector<vector<char> > mat) { return max(mat);}
  string _aurostd_initialize_max(const vector<vector<string> > mat) { return max(mat);}
  int _aurostd_initialize_max(const vector<vector<int> > mat) { return max(mat);}
  long _aurostd_initialize_max(const vector<vector<long> > mat) { return max(mat);}
  uint _aurostd_initialize_max(const vector<vector<uint> > mat) { return max(mat);}
  float _aurostd_initialize_max(const vector<vector<float> > mat) { return max(mat);}
  double _aurostd_initialize_max(const vector<vector<double> > mat) { return max(mat);}
  long double _aurostd_initialize_max(const vector<vector<long double> > mat) { return max(mat);}

  // ***************************************************************************
  // Function min of a vector/deque
  // ***************************************************************************
  // Stefano
  template<class utype> utype min(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec.at(i)<=out) out=vec.at(i);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_min(const vector<bool> vec) { return min(vec);}
  char _aurostd_initialize_min(const vector<char> vec) { return min(vec);}
  string _aurostd_initialize_min(const vector<string> vec) { return min(vec);}
  int _aurostd_initialize_min(const vector<int> vec) { return min(vec);}
  long _aurostd_initialize_min(const vector<long> vec) { return min(vec);}
  uint _aurostd_initialize_min(const vector<uint> vec) { return min(vec);}
  float _aurostd_initialize_min(const vector<float> vec) { return min(vec);}
  double _aurostd_initialize_min(const vector<double> vec) { return min(vec);}
  long double _aurostd_initialize_min(const vector<long double> vec) { return min(vec);}

  // Stefano
  template<class utype> utype min(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec.at(i)<=out) out=vec.at(i);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_min(const deque<bool> vec) { return min(vec);}
  char _aurostd_initialize_min(const deque<char> vec) { return min(vec);}
  string _aurostd_initialize_min(const deque<string> vec) { return min(vec);}
  int _aurostd_initialize_min(const deque<int> vec) { return min(vec);}
  long _aurostd_initialize_min(const deque<long> vec) { return min(vec);}
  uint _aurostd_initialize_min(const deque<uint> vec) { return min(vec);}
  float _aurostd_initialize_min(const deque<float> vec) { return min(vec);}
  double _aurostd_initialize_min(const deque<double> vec) { return min(vec);}
  long double _aurostd_initialize_min(const deque<long double> vec) { return min(vec);}

  // ***************************************************************************
  // Function min of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype min(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=mat.at(0).at(0);
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat.at(i).size();j++)
	if(mat.at(i).at(j)<=out) out=mat.at(i).at(j);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_min(const vector<vector<bool> > mat) { return min(mat);}
  char _aurostd_initialize_min(const vector<vector<char> > mat) { return min(mat);}
  string _aurostd_initialize_min(const vector<vector<string> > mat) { return min(mat);}
  int _aurostd_initialize_min(const vector<vector<int> > mat) { return min(mat);}
  long _aurostd_initialize_min(const vector<vector<long> > mat) { return min(mat);}
  uint _aurostd_initialize_min(const vector<vector<uint> > mat) { return min(mat);}
  float _aurostd_initialize_min(const vector<vector<float> > mat) { return min(mat);}
  double _aurostd_initialize_min(const vector<vector<double> > mat) { return min(mat);}
  long double _aurostd_initialize_min(const vector<vector<long double> > mat) { return min(mat);}

  // ***************************************************************************
  // Function sum of a vector/deque
  // ***************************************************************************
  template<class utype> utype sum(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec.at(i);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_sum(const vector<bool> vec) { return sum(vec);}
  char _aurostd_initialize_sum(const vector<char> vec) { return sum(vec);}
  string _aurostd_initialize_sum(const vector<string> vec) { return sum(vec);}
  int _aurostd_initialize_sum(const vector<int> vec) { return sum(vec);}
  long _aurostd_initialize_sum(const vector<long> vec) { return sum(vec);}
  uint _aurostd_initialize_sum(const vector<uint> vec) { return sum(vec);}
  float _aurostd_initialize_sum(const vector<float> vec) { return sum(vec);}
  double _aurostd_initialize_sum(const vector<double> vec) { return sum(vec);}
  long double _aurostd_initialize_sum(const vector<long double> vec) { return sum(vec);}

  template<class utype> utype sum(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec.at(i);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_sum(const deque<bool> vec) { return sum(vec);}
  char _aurostd_initialize_sum(const deque<char> vec) { return sum(vec);}
  string _aurostd_initialize_sum(const deque<string> vec) { return sum(vec);}
  int _aurostd_initialize_sum(const deque<int> vec) { return sum(vec);}
  long _aurostd_initialize_sum(const deque<long> vec) { return sum(vec);}
  uint _aurostd_initialize_sum(const deque<uint> vec) { return sum(vec);}
  float _aurostd_initialize_sum(const deque<float> vec) { return sum(vec);}
  double _aurostd_initialize_sum(const deque<double> vec) { return sum(vec);}
  long double _aurostd_initialize_sum(const deque<long double> vec) { return sum(vec);}

  // ***************************************************************************
  // Function sum of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype sum(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat.at(i).size();j++)
	out+=mat.at(i).at(j);
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_sum(const vector<vector<bool> > mat) { return sum(mat);}
  char _aurostd_initialize_sum(const vector<vector<char> > mat) { return sum(mat);}
  string _aurostd_initialize_sum(const vector<vector<string> > mat) { return sum(mat);}
  int _aurostd_initialize_sum(const vector<vector<int> > mat) { return sum(mat);}
  long _aurostd_initialize_sum(const vector<vector<long> > mat) { return sum(mat);}
  uint _aurostd_initialize_sum(const vector<vector<uint> > mat) { return sum(mat);}
  float _aurostd_initialize_sum(const vector<vector<float> > mat) { return sum(mat);}
  double _aurostd_initialize_sum(const vector<vector<double> > mat) { return sum(mat);}
  long double _aurostd_initialize_sum(const vector<vector<long double> > mat) { return sum(mat);}

  // ***************************************************************************
  // Function mean of a vector/deque
  // ***************************************************************************
  template<class utype> utype mean(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec.at(i);
    return (utype) out/((utype) vec.size());
  }
  // overload to force compiling
  bool _aurostd_initialize_mean(const vector<bool> vec) { return mean(vec);}
  char _aurostd_initialize_mean(const vector<char> vec) { return mean(vec);}
  int _aurostd_initialize_mean(const vector<int> vec) { return mean(vec);}
  long _aurostd_initialize_mean(const vector<long> vec) { return mean(vec);}
  uint _aurostd_initialize_mean(const vector<uint> vec) { return mean(vec);}
  float _aurostd_initialize_mean(const vector<float> vec) { return mean(vec);}
  double _aurostd_initialize_mean(const vector<double> vec) { return mean(vec);}
  long double _aurostd_initialize_mean(const vector<long double> vec) { return mean(vec);}

  template<class utype> utype mean(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec.at(i);
    return (utype) out/((utype) vec.size());
  }
  // overload to force compiling
  bool _aurostd_initialize_mean(const deque<bool> vec) { return mean(vec);}
  char _aurostd_initialize_mean(const deque<char> vec) { return mean(vec);}
  int _aurostd_initialize_mean(const deque<int> vec) { return mean(vec);}
  long _aurostd_initialize_mean(const deque<long> vec) { return mean(vec);}
  uint _aurostd_initialize_mean(const deque<uint> vec) { return mean(vec);}
  float _aurostd_initialize_mean(const deque<float> vec) { return mean(vec);}
  double _aurostd_initialize_mean(const deque<double> vec) { return mean(vec);}
  long double _aurostd_initialize_mean(const deque<long double> vec) { return mean(vec);}

  // ***************************************************************************
  // Function mean of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype mean(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat.at(i).size();j++)
	out+=mat.at(i).at(j);
    return (utype) out/((utype) mat.size()*mat.at(0).size());
  }
  // overload to force compiling
  bool _aurostd_initialize_mean(const vector<vector<bool> > mat) { return mean(mat);}
  char _aurostd_initialize_mean(const vector<vector<char> > mat) { return mean(mat);}
  int _aurostd_initialize_mean(const vector<vector<int> > mat) { return mean(mat);}
  long _aurostd_initialize_mean(const vector<vector<long> > mat) { return mean(mat);}
  uint _aurostd_initialize_mean(const vector<vector<uint> > mat) { return mean(mat);}
  float _aurostd_initialize_mean(const vector<vector<float> > mat) { return mean(mat);}
  double _aurostd_initialize_mean(const vector<vector<double> > mat) { return mean(mat);}
  long double _aurostd_initialize_mean(const vector<vector<long double> > mat) { return mean(mat);}

  // ***************************************************************************
  // Function reset of a vector/deque
  // ***************************************************************************
  template<class utype> vector<utype> reset(vector<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec.at(i)=(utype) 0;
    return vec;
  }
  // overload to force compiling
  vector<bool>  _aurostd_initialize_reset(vector<bool>& vec) { return reset(vec);}
  vector<char>  _aurostd_initialize_reset(vector<char>& vec) { return reset(vec);}
  vector<string>  _aurostd_initialize_reset(vector<string>& vec) { return reset(vec);}
  vector<int>   _aurostd_initialize_reset(vector<int>& vec) { return reset(vec);}
  vector<long>  _aurostd_initialize_reset(vector<long>& vec) { return reset(vec);}
  vector<uint>  _aurostd_initialize_reset(vector<uint>& vec) { return reset(vec);}
  vector<float> _aurostd_initialize_reset(vector<float>& vec) { return reset(vec);}
  vector<double>  _aurostd_initialize_reset(vector<double>& vec) { return reset(vec);}
  vector<long double>  _aurostd_initialize_reset(vector<long double>& vec) { return reset(vec);}

  template<class utype> deque<utype> reset(deque<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec.at(i)=(utype) 0;
    return vec;
  }
  // overload to force compiling
  deque<bool>  _aurostd_initialize_reset(deque<bool>& vec) { return reset(vec);}
  deque<char>  _aurostd_initialize_reset(deque<char>& vec) { return reset(vec);}
  deque<string>  _aurostd_initialize_reset(deque<string>& vec) { return reset(vec);}
  deque<int>   _aurostd_initialize_reset(deque<int>& vec) { return reset(vec);}
  deque<long>  _aurostd_initialize_reset(deque<long>& vec) { return reset(vec);}
  deque<uint>  _aurostd_initialize_reset(deque<uint>& vec) { return reset(vec);}
  deque<float> _aurostd_initialize_reset(deque<float>& vec) { return reset(vec);}
  deque<double>  _aurostd_initialize_reset(deque<double>& vec) { return reset(vec);}
  deque<long double>  _aurostd_initialize_reset(deque<long double>& vec) { return reset(vec);}

  // ***************************************************************************
  // Function reset of a vector<vector<>>
  // ***************************************************************************
  template<class utype> vector<vector<utype> > reset(vector<vector<utype> > mat) {
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat.at(i).size();j++)
	mat.at(i).at(j)=(utype) 0;
    return mat;
  }
  // overload to force compiling
  vector<vector<bool> > _aurostd_initialize_reset(vector<vector<bool> > mat) { return reset(mat);}
  vector<vector<char> > _aurostd_initialize_reset(vector<vector<char> > mat) { return reset(mat);}
  vector<vector<string> > _aurostd_initialize_reset(vector<vector<string> > mat) { return reset(mat);}
  vector<vector<int> > _aurostd_initialize_reset(vector<vector<int> > mat) { return reset(mat);}
  vector<vector<long> > _aurostd_initialize_reset(vector<vector<long> > mat) { return reset(mat);}
  vector<vector<uint> > _aurostd_initialize_reset(vector<vector<uint> > mat) { return reset(mat);}
  vector<vector<float> > _aurostd_initialize_reset(vector<vector<float> > mat) { return reset(mat);}
  vector<vector<double> > _aurostd_initialize_reset(vector<vector<double> > mat) { return reset(mat);}
  vector<vector<long double> > _aurostd_initialize_reset(vector<vector<long double> > mat) { return reset(mat);}

  // ***************************************************************************
  // Function clear of a vector/deque
  // ***************************************************************************
  template<class utype> vector<utype> clear(vector<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec.at(i)=(utype) 0;
    return vec;
  }
  // overload to force compiling
  vector<bool>  _aurostd_initialize_clear(vector<bool>& vec) { return clear(vec);}
  vector<char>  _aurostd_initialize_clear(vector<char>& vec) { return clear(vec);}
  vector<string>  _aurostd_initialize_clear(vector<string>& vec) { return clear(vec);}
  vector<int>   _aurostd_initialize_clear(vector<int>& vec) { return clear(vec);}
  vector<long>  _aurostd_initialize_clear(vector<long>& vec) { return clear(vec);}
  vector<uint>  _aurostd_initialize_clear(vector<uint>& vec) { return clear(vec);}
  vector<float> _aurostd_initialize_clear(vector<float>& vec) { return clear(vec);}
  vector<double>  _aurostd_initialize_clear(vector<double>& vec) { return clear(vec);}
  vector<long double>  _aurostd_initialize_clear(vector<long double>& vec) { return clear(vec);}

  template<class utype> deque<utype> clear(deque<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec.at(i)=(utype) 0;
    return vec;
  }
  // overload to force compiling
  deque<bool>  _aurostd_initialize_clear(deque<bool>& vec) { return clear(vec);}
  deque<char>  _aurostd_initialize_clear(deque<char>& vec) { return clear(vec);}
  deque<string>  _aurostd_initialize_clear(deque<string>& vec) { return clear(vec);}
  deque<int>   _aurostd_initialize_clear(deque<int>& vec) { return clear(vec);}
  deque<long>  _aurostd_initialize_clear(deque<long>& vec) { return clear(vec);}
  deque<uint>  _aurostd_initialize_clear(deque<uint>& vec) { return clear(vec);}
  deque<float> _aurostd_initialize_clear(deque<float>& vec) { return clear(vec);}
  deque<double>  _aurostd_initialize_clear(deque<double>& vec) { return clear(vec);}
  deque<long double>  _aurostd_initialize_clear(deque<long double>& vec) { return clear(vec);}

  // ***************************************************************************
  // Function clear of a vector<vector<>>
  // ***************************************************************************
  template<class utype> vector<vector<utype> > clear(vector<vector<utype> > mat) {
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat.at(i).size();j++)
	mat.at(i).at(j)=(utype) 0;
    return mat;
  }
  // overload to force compiling
  vector<vector<bool> > _aurostd_initialize_clear(vector<vector<bool> > mat) { return clear(mat);}
  vector<vector<char> > _aurostd_initialize_clear(vector<vector<char> > mat) { return clear(mat);}
  vector<vector<string> > _aurostd_initialize_clear(vector<vector<string> > mat) { return clear(mat);}
  vector<vector<int> > _aurostd_initialize_clear(vector<vector<int> > mat) { return clear(mat);}
  vector<vector<long> > _aurostd_initialize_clear(vector<vector<long> > mat) { return clear(mat);}
  vector<vector<uint> > _aurostd_initialize_clear(vector<vector<uint> > mat) { return clear(mat);}
  vector<vector<float> > _aurostd_initialize_clear(vector<vector<float> > mat) { return clear(mat);}
  vector<vector<double> > _aurostd_initialize_clear(vector<vector<double> > mat) { return clear(mat);}
  vector<vector<long double> > _aurostd_initialize_clear(vector<vector<long double> > mat) { return clear(mat);}

  // ***************************************************************************
  // Function random_shuffle of a vector/deque
  // ***************************************************************************
  template<class utype> void random_shuffle(vector<utype>& vec) {
    std::random_shuffle(vec.begin(),vec.end());
  }
  // overload to force compiling
  void _aurostd_initialize_random_shuffle(vector<bool>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<char>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<string>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<int>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<long>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<uint>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<float>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<double>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<long double>& vec) {random_shuffle(vec);}

  template<class utype> void random_shuffle(deque<utype>& vec) {
    std::random_shuffle(vec.begin(),vec.end());
  }
  // overload to force compiling
  void _aurostd_initialize_random_shuffle(deque<bool>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<char>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<string>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<int>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<long>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<uint>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<float>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<double>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<long double>& vec) {random_shuffle(vec);}

  // ***************************************************************************
  // Function isequal of vector vector
  // ***************************************************************************
  template<class utype> bool identical(vector<utype> vec1,vector<utype> vec2,utype epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }

  template<class utype> bool identical(deque<utype> vec1,deque<utype> vec2,utype epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }

  bool identical(vector<int> vec1,vector<int> vec2,int epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }

  bool identical(deque<int> vec1,deque<int> vec2,int epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }

  //   // overload to force compiling
  //  bool _aurostd_initialize_isequal(vector<bool> v1,vector<bool> v2,bool epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<char> v1,vector<char> v2,char epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<int> v1,vector<int> v2,int epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<long> v1,vector<long> v2,long epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<uint> v1,vector<uint> v2,uint epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<float> v1,vector<float> v2,float epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<double> v1,vector<double> v2,double epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<long double> v1,vector<long double> v2,long double epsilon) { return isequal(v1,v2,epsilon);}

  // ***************************************************************************
  // Function isequal of vector<vector<>> vector<vector<>>
  // ***************************************************************************
  template<class utype> bool identical(const vector<vector<utype> >& mat1,const vector<vector<utype> >& mat2,utype epsilon) {
    if(mat1.size()!=mat2.size()) return FALSE;
    for(uint i=0;i<mat1.size();i++) {
      if(mat1[i].size()!=mat2[i].size()) return FALSE;
      for(uint j=0;j<mat1[i].size();j++)
	if(aurostd::abs(mat1[i][j]-mat2[i][j])>epsilon) return FALSE;
    }
    return TRUE;
  }

  //   // overload to force compiling
  //   bool _aurostd_initialize_isequal(const vector<vector<bool> >& m1,const vector<vector<bool> >& m2,bool epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<char> >& m1,const vector<vector<char> >& m2,char epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<int> >& m1,const vector<vector<int> >& m2,int epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<long> >& m1,const vector<vector<long> >& m2,long epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<uint> >& m1,const vector<vector<uint> >& m2,uint epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<float> >& m1,const vector<vector<float> >& m2,float epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<double> >& m1,const vector<vector<double> >& m2,double epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<long double> >& m1,const vector<vector<long double> >& m2,long double epsilon) { return isequal(m1,m2,epsilon);}

  // ***************************************************************************
  // Function isequal of vector<vector<vector<>>> vector<vector<vector<>>>
  // ***************************************************************************
  template<class utype> bool identical(const vector<vector<vector<utype> > >& t1,const vector<vector<vector<utype> > >& t2,utype epsilon) {
    if(t1.size()!=t2.size()) return FALSE;
    for(uint i=0;i<t1.size();i++) {
      if(t1[i].size()!=t2[i].size()) return FALSE;
      for(uint j=0;j<t1[i].size();j++) {
	if(t1[i][j].size()!=t2[i][j].size()) return FALSE;
	for(uint k=0;k<t1[i][i].size();k++) 	
	  if(aurostd::abs(t1[i][j][k]-t2[i][j][k])>epsilon) return FALSE;
      }
    }
    return TRUE;
  }
  //   // overload to force compiling
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<bool> > >& t1,const vector<vector<vector<bool> > >& t2,bool epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<char> > >& t1,const vector<vector<vector<char> > >& t2,char epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<int> > >& t1,const vector<vector<vector<int> > >& t2,int epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<long> > >& t1,const vector<vector<vector<long> > >& t2,long epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<uint> > >& t1,const vector<vector<vector<uint> > >& t2,uint epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<float> > >& t1,const vector<vector<vector<float> > >& t2,float epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<double> > >& t1,const vector<vector<vector<double> > >& t2,double epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<long double> > >& t1,const vector<vector<vector<long double> > >& t2,long double epsilon) { return isequal(t1,t2,epsilon);}

  // ***************************************************************************
  // Function toupper/tolower
  // ***************************************************************************
  string toupper(const string& in) {
    string out(in);
    for(uint i=0;i<out.length();i++)
      out.at(i)=std::toupper(out.at(i));
    return out;
  }
  
  string tolower(const string& in) {
    string out(in);
    for(uint i=0;i<out.length();i++)
      out.at(i)=std::tolower(out.at(i));
    return out;
  }

  char toupper(const char& in) {
    return std::toupper(in);
  }
  
  char tolower(const char& in) {
    return std::tolower(in);
  }
  

  // ***************************************************************************
  // Function GetNumFields
  // ***************************************************************************
  // Dane Morgan
  int GetNumFields(const string& s) {
    int nf=0;
    int in_a_field=0;
    for(uint i=0;i<s.size();i++) {
      if(!(in_a_field) && s[i]!=' ') {
	in_a_field=1;
	nf++;
      }
      if(in_a_field && s[i]==' ') {
	in_a_field=0;
      }
    }
    return nf;
  }

  // ***************************************************************************
  // Function GetNextVal
  // ***************************************************************************
  // Dane Morgan - Stefano Curtarolo
  string GetNextVal(const string& s, int& id) {
    string ss;
    int i=id;
    while (s[i]==' ') {i++;} // ignore leading spaces.
    while (i<(int) s.size() && s[i]!=' ') { // pull out all text until next space.
      ss+=s[i]; i++;
    }
    id=i;
    return ss;
  }

  // ***************************************************************************
  // Function PaddedNumString
  // ***************************************************************************
  // Dane Morgan - Stefano Curtarolo
  string PaddedNumString(const int num, const int ndigits) {
    ostringstream oss;
    oss << std::setw(ndigits) << std::setfill('0') << num;// << ends;
    return oss.str();
  }

  // ***************************************************************************
  // Function PaddedPRE
  // ***************************************************************************
  // Add PRE characters to pad
  string PaddedPRE(string input,int depth,string ch) {
    stringstream aus("");
    aus << input;
    string strout="";
    for(int i=0;i<depth-(int) aus.str().length();i++) { strout+=ch; }
    strout+=aus.str();
    return strout;
  } 
  template<class utype> string PaddedPRE(utype input,int depth,string ch) {
    stringstream sss;
    sss << input;
    return PaddedPRE(sss.str(),depth,ch);
  }

  // ***************************************************************************
  // Function PaddedPOST
  // ***************************************************************************
  // Add POST characters to pad
  string PaddedPOST(string input,int depth,string ch) {
    stringstream aus("");
    aus << input;
    string strout=aus.str();
    for(int i=0;i<depth-(int) aus.str().length();i++) { strout+=ch; }
    return strout;
  }
  template<class utype> string PaddedPOST(utype input,int depth,string ch) {
    stringstream sss;
    sss << input;
    return PaddedPOST(sss.str(),depth,ch);
  }

  // ***************************************************************************
  // Function PaddedCENTER
  // ***************************************************************************
  // Add PRE AND POST characters to pad so that string is in the center
  string PaddedCENTER(string input,int depth,string ch) {
    stringstream aus("");
    int pre=(depth-(int) input.length())/2;
    int post=depth-pre-(int) input.length();
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: input.length()=" << input.length() << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: pre=" << pre << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: post=" << post << endl;   
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: depth=" << depth << endl;   
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: pre+post+input.length()=" << pre+post+input.length() << endl;   
    for(int i=1;i<pre;i++) { aus << ch; }  
    aus << input;
    for(int i=1;i<post;i++) { aus << ch; }  
    return aus.str();
  }
  template<class utype> string PaddedCENTER(utype input,int depth,string ch) {
    stringstream sss;
    sss << input;
    return PaddedCENTER(sss.str(),depth,ch);
  }
// ***************************************************************************
  // Function CleanStringASCII
  // ***************************************************************************
  // Clean a string from ASCII junk
  // Stefano Curtarolo
  string CleanStringASCII(const string& s) {
    string ss="";
    for(uint i=0;i<s.length();i++) {
      if(s[i]>='A' && s[i]<='Z') ss+=s[i];  // LETTERS
      if(s[i]>='a' && s[i]<='z') ss+=s[i];  // letters
      if(s[i]>='0' && s[i]<='9') ss+=s[i];  // numbers
      if(s[i]=='.' || s[i]=='+' || s[i]=='-' || s[i]=='*' || s[i]=='/') ss+=s[i];  // operations
      if(s[i]=='_' || s[i]=='#' || s[i]=='&' || s[i]==':' || s[i]==',' || s[i]=='@' || s[i]=='$') ss+=s[i];  // underscore
      if(s[i]=='=' || s[i]=='|' || s[i]=='\'' || s[i]==' ') ss+=s[i];  // underscore
    }
    return ss;
  }

  // ***************************************************************************
  // Function CGI_StringClean
  // ***************************************************************************
  // Clean a string from CGI junk
  string CGI_StringClean(const string& stringIN) {
    string stringOUT=stringIN;
    aurostd::StringSubst(stringOUT,"%0D%0A","\n");aurostd::StringSubst(stringOUT,"%0d%0a","\n");   // newlines
    aurostd::StringSubst(stringOUT,"+"," ");    // spaces
    aurostd::StringSubst(stringOUT,"%28","(");aurostd::StringSubst(stringOUT,"%29",")");   // ()
    aurostd::StringSubst(stringOUT,"%5B","[");aurostd::StringSubst(stringOUT,"%5D","]");   // []
    aurostd::StringSubst(stringOUT,"%7B","{");aurostd::StringSubst(stringOUT,"%7D","}");   // {}
    aurostd::StringSubst(stringOUT,"%2B","+");aurostd::StringSubst(stringOUT,"%2F","/");   //  operations
    aurostd::StringSubst(stringOUT,"%23","#");aurostd::StringSubst(stringOUT,"%21","!");
    aurostd::StringSubst(stringOUT,"%3F","?");aurostd::StringSubst(stringOUT,"%2C",",");
    aurostd::StringSubst(stringOUT,"%3A",":");aurostd::StringSubst(stringOUT,"%3B",";");
    aurostd::StringSubst(stringOUT,"%27","'");aurostd::StringSubst(stringOUT,"%22","\"");
    aurostd::StringSubst(stringOUT,"%60","`");aurostd::StringSubst(stringOUT,"%40","@");
    aurostd::StringSubst(stringOUT,"%24","$");aurostd::StringSubst(stringOUT,"%25","%");
    aurostd::StringSubst(stringOUT,"%5E","^");aurostd::StringSubst(stringOUT,"%26","&");
    aurostd::StringSubst(stringOUT,"%3D","=");aurostd::StringSubst(stringOUT,"%7E","~");
    aurostd::StringSubst(stringOUT,"%5C","\\");aurostd::StringSubst(stringOUT,"%7C","|");
    aurostd::StringSubst(stringOUT,"%3C","<");aurostd::StringSubst(stringOUT,"%3E",">");  // <>
    aurostd::StringSubst(stringOUT,"\n\n","\n");
    return stringOUT;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpaces
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpaces(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) if(s[i]!=' ' && s[i]!='\t') ss+=s[i];
    return ss;
  }
  string RemoveWhiteSpaces(const string& s, const char toogle) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    bool copy=TRUE;
    for (uint i=0;i<s.size();i++) {
      if(s[i]==toogle) copy=!copy;
      if(copy) if(s[i]!=' ' && s[i]!='\t') ss+=s[i];
      if(!copy) ss+=s[i];
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpacesFromTheBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[ss.size()-1]==' ' || ss[ss.size()-1]=='\t') {
      ss.erase(ss.size()-1,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheFront
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Oses
  string RemoveWhiteSpacesFromTheFront(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[0]==' ' || ss[ss.size()-1]=='\t') {
      ss.erase(0,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }
  
  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheFrontAndBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Oses
  string RemoveWhiteSpacesFromTheFrontAndBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    ss=RemoveWhiteSpacesFromTheBack(ss);
    ss=RemoveWhiteSpacesFromTheFront(ss);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSpaces
  // ***************************************************************************
  // Removes all spaces from a string. Morgan / Curtarolo
  string RemoveSpaces(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) if(s[i]!=' ') ss+=s[i];
    return ss;
  }
  string RemoveSpaces(const string& s, const char toogle) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    bool copy=TRUE;
    for (uint i=0;i<s.size();i++) {
      if(s[i]==toogle) copy=!copy;
      if(copy) if(s[i]!=' ') ss+=s[i];
      if(!copy) ss+=s[i];
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSpacesFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveSpacesFromTheBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[ss.size()-1]==' ') {
      ss.erase(ss.size()-1,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTabs
  // ***************************************************************************
  // Removes all tabs from a string. Stefano Curtarolo
  string RemoveTabs(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) if(s[i]!='\t') ss+=s[i];
    return ss;
  }
  string RemoveTabs(const string& s, const char toogle) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    bool copy=TRUE;
    for (uint i=0;i<s.size();i++) {
      if(s[i]==toogle) copy=!copy;
      if(copy) if(s[i]!='\t') ss+=s[i];
      if(!copy) ss+=s[i];
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTabsFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string.
  // Dane Morgan / Stefano Curtarolo
  string RemoveTabsFromTheBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[ss.size()-1]=='\t') {
      ss.erase(ss.size()-1,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveComments
  // ***************************************************************************
  // Removes all comments from a string.
  // Stefano Curtarolo
  //   string RemoveComments(const string& s) {
  //     if(s.size()==0) return s;  // nothing to do
  //     string ss;
  //     bool copy=TRUE;
  //     for (uint i=0;i<s.size();i++) {
  //       if(s[i]=='#')  copy=FALSE;
  //       if(s[i]=='\n') copy=TRUE;
  //       if(copy) ss+=s[i];
  //     }
  //     return ss;
  //   }
  
  string RemoveComments(const string &strin) {
    string strout=strin;
    vector<string> vstrout;aurostd::string2vectorstring(strout,vstrout);
    strout.clear();
    for(uint i=0;i<vstrout.size();i++) {
      if(vstrout.at(i).find(COMMENT_NEGLECT_1)!=string::npos) 
	vstrout.at(i)=vstrout.at(i).substr(0,vstrout.at(i).find(COMMENT_NEGLECT_1));  
      if(vstrout.at(i).find(COMMENT_NEGLECT_2)!=string::npos && vstrout.at(i).find(":"+COMMENT_NEGLECT_2)==string::npos)  // look for // but dont touch ://
	vstrout.at(i)=vstrout.at(i).substr(0,vstrout.at(i).find(COMMENT_NEGLECT_2));
      if(vstrout.at(i).find(COMMENT_NEGLECT_3)!=string::npos)
	vstrout.at(i)=vstrout.at(i).substr(0,vstrout.at(i).find(COMMENT_NEGLECT_3));
      //  if(!vstrout.at(i).empty()) cout << vstrout.at(i) << endl;
      if(!vstrout.at(i).empty()) 
	strout+=vstrout.at(i)+"\n";
    }  
    return strout;
  }
  
  // ***************************************************************************
  // Function RemoveCharacter
  // ***************************************************************************
  // Removes charecters from string
  // Stefano Curtarolo
  string RemoveCharacter(const string& s, const char character) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) {
      if(s[i]!=character) ss+=s[i];
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveNumbers
  // ***************************************************************************
  // Removes numbers from string
  // Stefano Curtarolo
  string RemoveNumbers(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    ss=RemoveCharacter(ss,'0');ss=RemoveCharacter(ss,'1');ss=RemoveCharacter(ss,'2');
    ss=RemoveCharacter(ss,'3');ss=RemoveCharacter(ss,'4');ss=RemoveCharacter(ss,'5');
    ss=RemoveCharacter(ss,'6');ss=RemoveCharacter(ss,'7');ss=RemoveCharacter(ss,'8');
    ss=RemoveCharacter(ss,'9');ss=RemoveCharacter(ss,'.');
    return ss;
  }

  // ***************************************************************************
  // Function RemoveRounding
  // ***************************************************************************
  // Removes rounding from string
  // Stefano Curtarolo
  string RemoveRounding(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    ss=RemoveSubString(ss,"(0)");ss=RemoveSubString(ss,"(1)");ss=RemoveSubString(ss,"(2)");
    ss=RemoveSubString(ss,"(3)");ss=RemoveSubString(ss,"(4)");ss=RemoveSubString(ss,"(5)");
    ss=RemoveSubString(ss,"(6)");ss=RemoveSubString(ss,"(7)");ss=RemoveSubString(ss,"(8)");
    ss=RemoveSubString(ss,"(9)");
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSubStringFirst
  // ***************************************************************************
  // Removes the first substring from string
  // Stefano Curtarolo
  string RemoveSubStringFirst(const string& str_orig, const string& str_rm) {
    string t=str_orig;
    std::string::size_type i = t.find(str_rm);
    if(i != std::string::npos)
      t.erase(i, str_rm.length( ));
    return t;
  }

  // ***************************************************************************
  // Function RemoveSubString
  // ***************************************************************************
  // Removes all instances of substring from string
  // Stefano Curtarolo
  string RemoveSubString(const string& str_orig, const string& str_rm) {
    string t=str_orig;
    string::size_type i;
    while(t.find(str_rm)!=string::npos) {
      i = t.find(str_rm);
      if(i != std::string::npos) t.erase(i, str_rm.length( ));
    };
    return t;
  }

  
  // ***************************************************************************
  // Function DirectoryMake
  // ***************************************************************************
  // Stefano Curtarolo
  // Make a directory by splitting each part.. returns false if something wrong
  // happens.
  bool DirectoryMake(string _Directory) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string Directory(CleanFileName(_Directory));
    std::vector<string> tokens;
    string dir_tmp,command;
    aurostd::string2tokens(Directory,tokens,"/");
    if(Directory.at(0)=='/') tokens[0]="/"+tokens[0]; // fix the root thing
    dir_tmp="";
    for(uint i=0;i<tokens.size();i++) {
      dir_tmp+=tokens.at(i)+"/";
      if(!aurostd::FileExist(dir_tmp)) {
	command="mkdir "+dir_tmp;
       	if(LDEBUG) cerr << "DirectoryMake creating directory=" <<  command << endl;
	aurostd::execute(command);
	if(!aurostd::FileExist(dir_tmp)) {
	  // if(!QUIET) cout << "EEEEE   can not make directory: " <<  command << endl;
	  return FALSE; // found some error in making directory
	}
      }
    }
    return TRUE;
  }

  // ***************************************************************************
  // Function SSH_DirectoryMake
  // ***************************************************************************
  // Stefano Curtarolo
  // Make a directory by splitting each part.. cant check much remotely
  // it starts from the assumption that you DO HAVE access to the remove account
  // which does not pretend a password
  bool SSH_DirectoryMake(string user, string machine,string _Directory) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string Directory(CleanFileName(_Directory));
    std::vector<string> tokens;
    string dir_tmp,command;
    aurostd::string2tokens(Directory,tokens,"/");
    if(Directory.at(0)=='/') tokens[0]="/"+tokens[0]; // fix the root thing
    dir_tmp="";
    for(uint i=0;i<tokens.size();i++) {
      dir_tmp+=tokens.at(i)+"/";
      command="ssh "+user+"@"+machine+" mkdir "+dir_tmp;
      if(LDEBUG) cerr << "SSH_DirectoryMake creating directory=" <<  command << endl;
      aurostd::execute(command);
    }
    return TRUE;
  }

  // ***************************************************************************
  // Function DirectoryChmod
  // ***************************************************************************
  // Stefano Curtarolo
  // return FALSE if something got messed up
  bool DirectoryChmod(string chmod_string,string _Directory) {
    string Directory(CleanFileName(_Directory));
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << CHMOD_BIN << " " << chmod_string << " " << Directory << endl;
    aurostd::execute(aus);
    return TRUE;
  }

  // ***************************************************************************
  // Function DirectoryLS
  // ***************************************************************************
  // Stefano Curtarolo
  // Returns the content of a directory without "." and ".."
  bool DirectoryLS(string _Directory,vector<string> &vfiles) {
    vfiles.clear();
    string Directory(CleanFileName(_Directory));
    DIR *dp;
    struct dirent *ep;
    string file;
    dp=opendir(Directory.c_str());
    if(dp!=NULL) {
      while((ep=readdir(dp))) {
	file=(ep->d_name);
	if(file!="." && file!="..")
	  vfiles.push_back(file);
      }
      (void) closedir (dp);
    } else {
      cerr << "ERROR: aurostd::DirectoryLS Couldn't open the directory: " << Directory << endl;
      return FALSE;
    }
    return TRUE;
  }

  // ***************************************************************************
  // Function DirectoryLocked
  // ***************************************************************************
  bool DirectoryLocked(string directory,string LOCK) {
    if(FileExist(directory+"/"+LOCK)) return TRUE;
    if(FileExist(directory+"/"+LOCK+".gz")) return TRUE;
    if(FileExist(directory+"/"+LOCK+".bz2")) return TRUE;
    return FALSE;
  }

  // ***************************************************************************
  // Function DirectorySkipped
  // ***************************************************************************
  bool DirectorySkipped(string directory) {
    if(FileExist(directory+"/SKIP")) return TRUE;
    if(FileExist(directory+"/SKIP.gz")) return TRUE;
    if(FileExist(directory+"/SKIP.bz2")) return TRUE;
    return FALSE;
  }

  // ***************************************************************************
  // Function DirectoryWritable and DirectoryUnwritable
  // ***************************************************************************
  bool DirectoryWritable(string _Directory) {
    string Directory(CleanFileName(_Directory));
    string filename=string(Directory+"/aflow.writable."+aurostd::utype2string(uint((double) std::floor(100000*aurostd::ran0())))+".test");
    string2file("DirectoryWritable",filename);
    if(!FileExist(filename)) return FALSE;
    string command=string("rm -f "+filename);
    execute(command);
    return TRUE;
  }
  bool DirectoryUnwritable(string Directory) {
    return !DirectoryWritable(Directory);
  }
  
  // ***************************************************************************
  // Function CleanFileName
  // ***************************************************************************
  // Stefano Curtarolo
  // cleans file names from obvious things
  string CleanFileName(string fileIN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string fileOUT=fileIN;
    if(LDEBUG) cerr << "aurostd::CleanFileName: " << fileOUT << endl;
// [OBSOLETE] interferes with ~/.aflow.rc   if(aurostd::substring2bool(fileOUT,"~/")) aurostd::StringSubst(fileOUT,"~/","/home/"+XHOST.User+"/");
    aurostd::StringSubst(fileOUT,"//","/");
    aurostd::StringSubst(fileOUT,"/./","/");
    if(LDEBUG) cerr << "aurostd::CleanFileName: " << fileOUT << endl;
    return fileOUT;
  }
  
  // ***************************************************************************
  // Function ProperFileName
  // ***************************************************************************
  // Stefano Curtarolo
  // fix file names from obvious things
  string ProperFileName(string fileIN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string fileOUT=fileIN;
    if(LDEBUG) cerr << "aurostd::ProperFileName: " << fileOUT << endl;
    aurostd::StringSubst(fileOUT,"//",".");
    aurostd::StringSubst(fileOUT,"/",".");
    if(LDEBUG) cerr << "aurostd::ProperFileName: " << fileOUT << endl;
    return fileOUT;
  }
  
  // ***************************************************************************
  // Function CopyFile
  // ***************************************************************************
  // Stefano Curtarolo
  // copy the file but does not check if the directory can be made or not...
  bool CopyFile(string from,string to) {
    stringstream command;
    command << "cp -f " << CleanFileName(from) << " " << CleanFileName(to) << " " << endl;
    aurostd::execute(command);
    return TRUE;
  }
  
  //CO - START
  //***************************************************************************//
  // aurostd::MatchCompressed
  //***************************************************************************//
  bool MatchCompressed(const string& CompressedFileName,const string& FileNameOUT) {
    //Corey Oses
    //Will try to mimic compression of a given file, useful for overwriting
    if(!IsCompressed(CompressedFileName)&&!IsCompressed(FileNameOUT)) {return TRUE;}
    if(IsCompressed(FileNameOUT)) {
      if(GetCompressionExtension(CompressedFileName)==GetCompressionExtension(FileNameOUT)) {return TRUE;}
      return FALSE;}
    if(substring2bool(CompressedFileName,".bz2")) {BzipFile(FileNameOUT);return TRUE;}
    if(substring2bool(CompressedFileName,".gz")) {GzipFile(FileNameOUT);return TRUE;}
    if(substring2bool(CompressedFileName,".zip")) {ZipFile(FileNameOUT);return TRUE;}
    return FALSE;
  }

  //***************************************************************************//
  // aurostd::DecompressFile
  //***************************************************************************//
  bool DecompressFile(const string& CompressedFileName){
    //Corey Oses
    //try to decompress the file the ways we know how
    //try not to use unless necessary, prefer efile2stringstream, etc.
    //where you might use it: open binary files
    if(!IsCompressed(CompressedFileName)){return TRUE;}
    if(substring2bool(CompressedFileName,".bz2")) {BzipFile(CompressedFileName);return TRUE;}
    if(substring2bool(CompressedFileName,".gz")) {GzipFile(CompressedFileName);return TRUE;}
    if(substring2bool(CompressedFileName,".zip")) {ZipFile(CompressedFileName);return TRUE;}
    return FALSE;
  }

  //***************************************************************************//
  // aurostd::efile2tempfile
  //***************************************************************************//
  bool efile2tempfile(string _FileNameIN, string& FileNameOUT) {
    //Corey Oses
    //Decompresses file to temp file, and return its path
    //SAFE:  Will return FileNameIn if compression not needed
    string FileNameIN = aurostd::CleanFileName(_FileNameIN);
    if((!aurostd::FileExist(FileNameIN)) || aurostd::FileEmpty(FileNameIN)) {
      cerr << endl;
      cerr << "ERROR - efile2tempfile:  file=" << FileNameIN << " not present/is empty !" << endl;
      cerr << endl;
      return FALSE;
    }
    //check if compressed, and decompress
    if(!aurostd::IsCompressed(FileNameIN)) {
      FileNameOUT = FileNameIN;
      return TRUE;
    }
    stringstream FileNameIN_ss;
    aurostd::efile2stringstream(FileNameIN, FileNameIN_ss);
    FileNameOUT = aurostd::TmpFileCreate();
    aurostd::stringstream2file(FileNameIN_ss, FileNameOUT);
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::IsCompressed
  //***************************************************************************//
  bool IsCompressed(string FileNameIN,string& FileNameOUT) {
    //Corey Oses
    //Given a File, it will return name of decompressed variant
    //NB:  Does not actually decompress
    //FileNameOUT="";
    if(substring2bool(FileNameIN,".bz2")) {FileNameOUT=StringSubst(FileNameIN,".bz2","");return TRUE;}
    if(substring2bool(FileNameIN,".tar.gz")) {FileNameOUT=StringSubst(FileNameIN,".tar.gz","");return TRUE;}
    if(substring2bool(FileNameIN,".gz")) {FileNameOUT=StringSubst(FileNameIN,".gz","");return TRUE;}
    if(substring2bool(FileNameIN,".zip")) {FileNameOUT=StringSubst(FileNameIN,".zip","");return TRUE;}
    FileNameOUT=FileNameIN;return FALSE;
  }

  //***************************************************************************//
  // aurostd::IsCompressed
  //***************************************************************************//
  bool IsCompressed(string FileNameIN) {
    string FileNameOUT;
    return IsCompressed(FileNameIN,FileNameOUT);
  }

  //***************************************************************************//
  // aurostd::GetCompressionExtension
  //***************************************************************************//
  string GetCompressionExtension(const string& CompressedFileName) {
    //Corey Oses
    //Will determine zipped extension
    string extension="";
    if(!IsCompressed(CompressedFileName)) {return extension;}
    if(substring2bool(CompressedFileName,".bz2")) {return ".bz2";}
    if(substring2bool(CompressedFileName,".gz")) {return ".gz";}
    if(substring2bool(CompressedFileName,".zip")) {return ".zip";}
    return extension;
  }
  //CO - END

  // ***************************************************************************
  // Function Bunzip BzipFile
  // ***************************************************************************
  // Stefano Curtarolo
  // Bzip the file (does not check)
  bool BunzipFile(string _file) {
    string file(CleanFileName(_file));
    if(!aurostd::IsCommandAvailable("bzip2")) {
      cerr << "ERROR - BunzipFile:  command \"bzip2\" is necessary !" << endl;
      return FALSE;
    }   
    // string file(_file); //  cerr << file << endl;
    if(aurostd::substring2bool(file,".bz2")) aurostd::execute("bzip2 -dqf "+file);
    return TRUE;
  }
  bool BzipFile(string _file) {
    string file(CleanFileName(_file));
    if(!aurostd::IsCommandAvailable("bzip2")) {
      cerr << "ERROR - BzipFile:  command \"bzip2\" is necessary !" << endl;
      return FALSE;
    }   
    if(FileExist(file+".bz2")) {cerr << file << ".bz2" << endl;}
    if(FileExist(file+".bz2")) {aurostd::execute("rm -f "+file+".bz2");}
    if(!aurostd::substring2bool(file,".bz2")) aurostd::execute("bzip2 -9qf "+file);
    return TRUE;
  }

  //CO - START
  //***************************************************************************//
  // aurostd::GzipFile
  //***************************************************************************//
  bool GzipFile(const string& _file) {
    string file(CleanFileName(_file));
    if(!aurostd::IsCommandAvailable("gzip")) {
      cerr << "ERROR - GzipFile:  command \"gzip\" is necessary !" << endl;
      return FALSE;
    }
    if(FileExist(file+".gz")) {cerr << file << ".gz" << endl;}
    if(FileExist(file+".gz")) {aurostd::execute("rm -f "+file+".gz");}
    if(!aurostd::substring2bool(file,".gz")) aurostd::execute("gzip -9f "+file);
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::ZipFile
  //***************************************************************************//
  bool ZipFile(const string& _file) {
    string file(CleanFileName(_file));
    if(!aurostd::IsCommandAvailable("zip")) {
      cerr << "ERROR - ZipFile:  command \"zip\" is necessary !" << endl;
      return FALSE;
    }
    if(FileExist(file+".zip")) {cerr << file << ".zip" << endl;}
    if(FileExist(file+".zip")) {aurostd::execute("rm -f "+file+".zip");}
    if(!aurostd::substring2bool(file,".zip")) aurostd::execute("zip -9qm "+file+".zip "+file);
    return TRUE;
  }
  //CO - END

  // ***************************************************************************
  // Function FileExist
  // ***************************************************************************
  // Stefano Curtarolo
  // return a simple bool, nothing else, from a string (which is supposed to be
  // a file name)
  bool FileExist(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    //   cerr << FileName << endl;
    bool exist=FALSE;
    ifstream FileStream;
    FileStream.open(FileName.c_str(),std::ios::in);
    FileStream.clear();
    FileStream.close();
    if(FileStream.good()) {exist=TRUE;} else {exist=FALSE;}
    return exist;
  }

  // ***************************************************************************
  // Function EFileExist
  // ***************************************************************************
  // Stefano Curtarolo
  bool EFileExist(const string& FileName, string& FileNameOut) {
    //FileNameOut=FileName; //CO
    if(FileExist(FileName)) {FileNameOut=FileName;return TRUE;}
    return CompressedFileExist(FileName,FileNameOut); //CO
  }

  // ***************************************************************************
  // Function EFileExist
  // ***************************************************************************
  // Stefano Curtarolo
  bool EFileExist(const string& FileName) {
    string FileNameOut;
    return EFileExist(FileName,FileNameOut);
  }

  //CO - START
  // ***************************************************************************
  // Function CompressedFileExist
  // ***************************************************************************
  // Corey Oses
  // tells you if the compressed variant exists
  bool CompressedFileExist(const string& FileName, string& FileNameOut){
    if(FileExist(FileName+".bz2")) {FileNameOut=FileName+".bz2";return TRUE;}
    if(FileExist(FileName+".gz")) {FileNameOut=FileName+".gz";return TRUE;}
    if(FileExist(FileName+".zip")) {FileNameOut=FileName+".zip";return TRUE;}
    FileNameOut=FileName;return FALSE;
  }

  // ***************************************************************************
  // Function CompressedFileExist
  // ***************************************************************************
  // Corey Oses
  // tells you if the compressed variant exists
  bool CompressedFileExist(const string& FileName){
    string FileNameOut;
    return CompressedFileExist(FileName,FileNameOut);
  }
  //CO - END

  // ***************************************************************************
  // Function FileSize
  // ***************************************************************************
  // Stefano Curtarolo - jan 08
  // returns in bytes the size of a file
  int FileSize(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    ifstream FileStream;
    int sizeout;
    FileStream.open(FileName.c_str(),std::ios::in);
    if(!FileStream.good()) {
      sizeout=0;
    } else {
      string FileString; FileString="";char c; while (FileStream.get(c)) FileString+=c;
      sizeout=FileString.length();
    }
    FileStream.close();
    return sizeout;
  }

  // ***************************************************************************
  // Function FileEmpty && FileNotEmpty
  // ***************************************************************************
  // Stefano Curtarolo - jan 08
  // returns in bytes the size of a file
  bool FileEmpty(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    if(FileExist(FileName)==FALSE) return TRUE;  // does not exist hence empty
    // it exists
    if(1) {
      int i=0;
      ifstream FileStream;
      FileStream.open(FileName.c_str(),std::ios::in);
      char c; 
      while (FileStream.get(c)&&i<256) {i++;};
      // count no more that 16... it is not worth to count more
      FileStream.close();
      if(i>0) return FALSE;
      else return TRUE;
    }
    if(0) {
      if(FileSize(FileName)<=1) return TRUE;
      else return FALSE;
    }
    return FALSE;
  }
  bool FileNotEmpty(const string& FileName) {
    return !FileEmpty(FileName);
  }

  // ***************************************************************************
  // Function FileToString
  // ***************************************************************************
  // Loat the content of a file into a string
  string FileToString(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    ifstream FileStream;
    stringstream strstreamout;
    string strline;
    strstreamout.str(std::string());
    //  cerr << FileName.c_str() << endl;  // DEBUG
    FileStream.open(FileName.c_str(),std::ios::in);
    if(FileStream.good()) { // !=NULL) {
      while(getline(FileStream,strline)) {
	strstreamout << strline << endl;
      }
    }
    FileStream.clear();FileStream.close();
    return strstreamout.str();
  }

  // ***************************************************************************
  // Function InFileExistCheck
  // ***************************************************************************
  // Dane Morgan
  void InFileExistCheck(const string& routine, const string& _FileName,
			ifstream& file_to_check, std::ostream& outf) {
    string FileName(CleanFileName(_FileName));
    if(!file_to_check) {
      outf << "Aflow VERSION "<< string(AFLOW_VERSION) << endl;
      outf << " ERROR: in " << routine << endl;
      outf << " ERROR: Cannot open file: \"" << FileName << "\"" << endl;
      outf << " ERROR: Exiting" << endl;
      exit(1);
    }
  }

  // ***************************************************************************
  // Function IsCommandAvailable
  // ***************************************************************************
  // tells you if the command is available
  bool IsCommandAvailable(const string& command, string& position) {
    // position=aurostd::execute2string("which "+command+" 2>&1 2> /dev/null");
    position=aurostd::execute2string("bash -c \"which "+command+" 2>&1 2> /dev/null\"");
    // cerr << position.length() << endl;
    aurostd::StringSubst(position,"\n","");
    if(position.length()>0) return TRUE;
    if(aurostd::FileExist("./"+command)) {position="./"+command;return TRUE;}
    if(aurostd::FileExist("/bin/"+command)) {position="/bin/"+command;return TRUE;}
    if(aurostd::FileExist("/sbin/"+command)) {position="/sbin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/bin/"+command)) {position="/usr/bin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/sbin/"+command)) {position="/usr/sbin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/local/bin/"+command)) {position="/usr/local/bin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/local/sbin/"+command)) {position="/usr/local/sbin/"+command;return TRUE;}  // go around path
    position="";
    return FALSE;
  }

  bool IsCommandAvailable(const string& command) {
    string position;
    return aurostd::IsCommandAvailable(command,position);
  }

  bool IsCommandAvailableModify(string& location) {
    if(!aurostd::IsCommandAvailable(location)) return FALSE;
    else return aurostd::IsCommandAvailable(location,location);
  }

  // ***************************************************************************
  // Function CommandRequired
  // ***************************************************************************
  // tells you if the command is available
  bool CommandRequired(const string& command, string& position) {
    position=aurostd::execute2string("which "+command);
    aurostd::StringSubst(position,"\n","");
    if(position.length()>0) return TRUE;
    cerr << "ERROR: CommandRequired \"" << command << "\" is not available" << endl;
    exit(0);
    return FALSE;
  }

  bool CommandRequired(const string& command) {
    string position;
    return CommandRequired(command,position);
  }

  // ***************************************************************************
  // Function IsExecutableAvailable
  // ***************************************************************************
  // tells you if the executable is available
  bool IsExecutableAvailable(const string& executable, string& position) {
    return IsCommandAvailable(executable,position);
  }

  bool IsExecutableAvailable(const string& executable) {
    string position;
    return IsCommandAvailable(executable,position);
  }

  // ***************************************************************************
  // Function ExecutableRequired
  // ***************************************************************************
  // tells you if the executable is available
  bool ExecutableRequired(const string& executable, string& position) {
    return CommandRequired(executable,position);
  }

  bool ExecutableRequired(const string& executable) {
    string position;
    return CommandRequired(executable,position);
  }

  // ***************************************************************************
  // DeleteOstringStreams
  // ***************************************************************************
  void StringstreamClean(ostringstream &aus) {
    aus.seekp(0,ios_base::beg);       // RESET
    for(int i=0;i<BUFFER_MAXLEN;i++)  // RESET
      aus<<(char)0;                   // RESET
    aus.seekp(0,ios_base::beg);       // RESET
    // aus.str(std::string());
  }
  void StringstreamClean(stringstream &aus) {
    aus.seekp(0,ios_base::beg);       // RESET
    for(int i=0;i<BUFFER_MAXLEN;i++)  // RESET
      aus<<(char)0;                   // RESET
    aus.seekp(0,ios_base::beg);       // RESET
    // aus.str(std::string());
  }


  // ***************************************************************************
  // Function FindIfStringInStream
  // ***************************************************************************
  //  This function returns true if string is in stream
  //  (on one line), or otherwise false. The search starts
  //  at the present file pointer location.  Note that this
  //  does alter the input string, resetting the file pointer
  //  to the input value at the end.
  // Dane Morgan style

  int FindIfStringInStream(const string& key, std::istream& instream) {
    // Get file pointer location at entry
    int loc=instream.tellg();
    int found_match=0;
    string s;
    getline(instream,s);
    int cont=0;
    if(getline(instream,s)) cont=1;
    while (cont) {
      int id=s.find(key);
      if(id!=(int) s.npos) { // Found key
	cont=0;
	found_match=1;
      }
      if(!getline(instream,s)) cont=0;
    }
    // Clear any fail bits associated with searching.
    instream.clear();
    // Set file pointer location to entry value
    instream.seekg(loc);
    return found_match;
  }

  // ***************************************************************************
  // Print Messages Errors and Warnings on and off streams.
  // ***************************************************************************
#define ErrorBarString "EEEEE  ---------------------------------------------------------------------------------------------------------------------------- "

  // with ostringstream
  void PrintMessageStream(ofstream &FileERROR,ostringstream &stream,const bool &quiet) {
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(!quiet) {cout << stream.str().c_str();cout.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintMessageStream(std::ostream &FileERROR,ostringstream &stream,const bool &quiet) {
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(!quiet) {cout << stream.str().c_str();cout.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintMessageStream(ofstream &FileERROR,ostringstream &stream,const bool &quiet,const bool& osswrite,std::ostream& oss) {
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(osswrite) {if(!quiet) {oss << stream.str().c_str();oss.flush();}}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintMessageStream(ostringstream &stream,const bool &quiet) {
    if(!quiet) {cout << stream.str().c_str();cout.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(ofstream &FileERROR,ostringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileERROR.flush();
    cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(std::ostream &FileERROR,ostringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileERROR.flush();
    cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(ofstream &FileERROR,ostringstream &stream,const bool &quiet,const bool& osswrite,std::ostream& oss) {
    if(quiet) {;} // phony just to keep quiet busy
    if(osswrite) {;} // phony just to keep quiet busy
    FileERROR << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileERROR.flush();
    oss << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(ostringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(ofstream &FileERROR,ostringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << stream.str().c_str(); FileERROR.flush();
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(std::ostream &FileERROR,ostringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << stream.str().c_str(); FileERROR.flush();
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(ofstream &FileERROR,ostringstream &stream,const bool &quiet,const bool& osswrite,std::ostream& oss) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(osswrite) {oss << stream.str().c_str();oss.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(ostringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  // with stringstream
  void PrintMessageStream(ofstream &FileERROR,stringstream &stream,const bool &quiet) {
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(!quiet) {cout << stream.str().c_str();cout.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintMessageStream(std::ostream &FileERROR,stringstream &stream,const bool &quiet) {
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(!quiet) {cout << stream.str().c_str();cout.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintMessageStream(ofstream &FileERROR,stringstream &stream,const bool &quiet,const bool& osswrite,std::ostream& oss) {
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(osswrite) {if(!quiet) {oss << stream.str().c_str();oss.flush();}}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintMessageStream(stringstream &stream,const bool &quiet) {
    if(!quiet) {cout << stream.str().c_str();cout.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(ofstream &FileERROR,stringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileERROR.flush();
    cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(std::ostream &FileERROR,stringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileERROR.flush();
    cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(ofstream &FileERROR,stringstream &stream,const bool &quiet,const bool& osswrite,std::ostream& oss) {
    if(quiet) {;} // phony just to keep quiet busy
    if(osswrite) {;} // phony just to keep quiet busy
    FileERROR << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileERROR.flush();
    oss << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintErrorStream(stringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(ofstream &FileERROR,stringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << stream.str().c_str(); FileERROR.flush();
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(std::ostream &FileERROR,stringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << stream.str().c_str(); FileERROR.flush();
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(ofstream &FileERROR,stringstream &stream,const bool &quiet,const bool& osswrite,std::ostream& oss) {
    if(quiet) {;} // phony just to keep quiet busy
    FileERROR << stream.str().c_str(); FileERROR.flush();
    if(osswrite) {oss << stream.str().c_str();oss.flush();}
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  void PrintWarningStream(stringstream &stream,const bool &quiet) {
    if(quiet) {;} // phony just to keep quiet busy
    cout << stream.str().c_str();cout.flush();
    // cerr << stream.str().c_str(); cerr.flush();
    aurostd::StringstreamClean(stream);
  }

  // ***************************************************************************
  // Execute Streams/Strings/C_strings
  // ***************************************************************************
  bool execute(ostringstream &command) {
    // cerr << "COMMAND " <<  command.str().c_str() << endl;
    system(command.str().c_str());
    aurostd::StringstreamClean(command);
    return TRUE;
  }

  bool execute(stringstream &command) {
    // cerr << "COMMAND " <<  command.str().c_str() << endl;
    system(command.str().c_str());
    aurostd::StringstreamClean(command);
    return TRUE;
  }

  bool execute(string command) {
    // cerr << "COMMAND " <<  command.c_str() << endl;
    system(command.c_str());
    //   command="";
    return TRUE;
  }

#ifdef _stringcharstar_
  bool execute(char* command) {
    // cerr << "COMMAND " <<  command << endl;
    system(command);
    return TRUE;
  }
#endif

  // ***************************************************************************
  // Execute vectors/deque of Strings
  // ***************************************************************************
  bool execute(deque<string> vcommand) {
    for(uint i=0;i<vcommand.size();i++)
      execute(vcommand.at(i));
    return TRUE;
  }
  bool execute(vector<string> vcommand) {
    for(uint i=0;i<vcommand.size();i++)
      execute(vcommand.at(i));
    return TRUE;
  }

  // ***************************************************************************
  // Execute & Report Streams/Strings/C_strings
  // ***************************************************************************
  string execute2string(string command) {
    // bool INIT_VERBOSE=TRUE;
    // cerr << "COMMAND " <<  command << endl;
    stringstream strstream,cmdstream;
    string file=aurostd::TmpFileCreate("execute_report");
    cmdstream << command << " > " << file <<  endl;
    system(cmdstream.str().c_str());
    // command="";
    strstream << aurostd::file2string(file);
    cmdstream.clear();cmdstream.str(std::string());
#ifndef _AFLOW_TEMP_PRESERVE_
    aurostd::RemoveFile(file);
#endif
    string strout=strstream.str();
    if(strout.length()>0)
      if(strout.at(strout.length()-1)=='\n')
	strout.erase(strout.length()-1);
    return strout;
  }

  string execute2string(ostringstream &command) {
    string command_str=command.str();
    aurostd::StringstreamClean(command);
    return execute2string(command_str);
  }

  string execute2string(stringstream &command) {
    string command_str=command.str();
    aurostd::StringstreamClean(command);
    return execute2string(command_str);
  }

  vector<string> execute2string(vector<string> vcommand) {
    vector<string> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back(execute2string(vcommand.at(i)));
    return out;
  }

  deque<string> execute2string(deque<string> vcommand) {
    deque<string> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back(execute2string(vcommand.at(i)));
    return out;
  }

#ifdef _stringcharstar_
  string execute2string(char* command) {
    string command_str=string(command);
    return execute2string(command_str);
  }
#endif

  // ***************************************************************************
  // Execute & Report Int Streams/Strings/C_strings
  // ***************************************************************************
  template<class utype> utype execute2utype(ostringstream &command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template<class utype> utype execute2utype(stringstream &command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template<class utype> utype execute2utype(string command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template<class utype> vector<utype> execute2utype(vector<utype> vcommand) {
    vector<utype> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back((utype) execute2utype<utype>(vcommand.at(i)));
    return out;
  }

  template<class utype> deque<utype> execute2utype(deque<utype> vcommand) {
    deque<utype> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back((utype) execute2utype<utype>(vcommand.at(i)));
    return out;
  }

#ifdef _stringcharstar_
  template<class utype> utype execute2utype(char* command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }
#endif

  // ***************************************************************************
  // Sleep
  // ***************************************************************************
  unsigned int Sleep(unsigned int seconds) {
    //  ostringstream aus;
    // aus << "sleep " << (int) seconds << " " << endl;
    // aurostd::execute(aus);
    return sleep(seconds);
  }

  // *******************************************************************************************
  // *******************************************************************************************
  bool ExtractToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,string Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ofstream FileOUTPUT;
    FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
    string strline,subS2;
    subS2=Keyword; // STRING TO SEARCH
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline))
      if(aurostd::substring2bool(strline,subS2)) {
	FileOUTPUT << strline.substr(strline.find(subS2)+subS2.length()) << endl;
	status=TRUE;
      }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
    return status;  // return FALSE if the keyword was never found
  }

  bool ExtractToFileEXPLICIT(string StringIN,string FileNameOUTPUT,string Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ofstream FileOUTPUT;
    FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
    string strline,subS2;
    subS2=Keyword; // STRING TO SEARCH
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      strline=tokens.at(i);
      if(aurostd::substring2bool(strline,subS2)) {
	FileOUTPUT << strline.substr(strline.find(subS2)+subS2.length()) << endl;
	status=TRUE;
      }
    }
    FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
    return status;  // return FALSE if the keyword was never found
  }

  // *******************************************************************************************
  bool ExtractToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,string Keyword_start,string Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ofstream FileOUTPUT;
    FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
    string strline,subS2;
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline)) {
      if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
      if(status) FileOUTPUT << strline << endl;
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
    return status;  // return FALSE if something got messed up
  }

  bool ExtractToFileEXPLICIT(string StringIN,string FileNameOUTPUT,string Keyword_start,string Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ofstream FileOUTPUT;
    FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
    string strline,subS2;
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      strline=tokens.at(i);
      if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
      if(status) FileOUTPUT << strline << endl;
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
    return status;  // return FALSE if something got messed up
  }

  // *******************************************************************************************
  // *******************************************************************************************
  bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,string Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string strline,subS2;
    subS2=Keyword; // STRING TO SEARCH
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline))
      if(aurostd::substring2bool(strline,subS2)) {
	StringOUTPUT=StringOUTPUT+strline.substr(strline.find(subS2)+subS2.length())+"\n";
	status=TRUE;
      }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status;  // return FALSE if the keyword was never found
  }

  bool ExtractToStringEXPLICIT(string StringIN,string& StringOUTPUT,string Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string strline,subS2;
    subS2=Keyword; // STRING TO SEARCH
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      strline=tokens.at(i);
      if(aurostd::substring2bool(strline,subS2)) {
	StringOUTPUT=StringOUTPUT+strline.substr(strline.find(subS2)+subS2.length())+"\n";
	status=TRUE;
      }
    }
    return status;  // return FALSE if the keyword was never found
  }

  // *******************************************************************************************
  bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,string Keyword_start,string Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string strline,subS2;
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline)) {
      if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
      if(status) StringOUTPUT=StringOUTPUT+strline+"\n";
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status;  // return FALSE if something got messed up
  }

  bool ExtractToStringEXPLICIT(string StringIN,string& StringOUTPUT,string Keyword_start,string Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string strline,subS2;
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      strline=tokens.at(i);
      if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
      if(status) StringOUTPUT=StringOUTPUT+strline+"\n";
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    return status;  // return FALSE if something got messed up
  }

  // *******************************************************************************************
  // *******************************************************************************************
  bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    string strline,subS2;
    subS2=Keyword; // STRING TO SEARCH
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline))
      if(aurostd::substring2bool(strline,subS2)) {
	StringstreamOUTPUT << strline.substr(strline.find(subS2)+subS2.length()) << endl;
	status=TRUE;
      }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status;  // return FALSE if the keyword was never found
  }

  // *******************************************************************************************
  bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop) {    // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    string strline;
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline)) {
      if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
      if(status) StringstreamOUTPUT << strline << endl;
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status;  // return FALSE if something got messed up
  }

  bool ExtractToStringstreamEXPLICIT(stringstream StringStreamIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop) { // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    string StringIN=StringStreamIN.str();
    return ExtractToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start,Keyword_stop);
  }

  bool ExtractToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop) {  // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      if(aurostd::substring2bool(tokens.at(i),Keyword_stop))  status=FALSE;
      if(status) StringstreamOUTPUT << tokens.at(i) << endl;
      if(aurostd::substring2bool(tokens.at(i),Keyword_start)) status=TRUE;
    }
    return status;  // return FALSE if something got messed up
  }

  bool ExtractToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword) {  // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      if(aurostd::substring2bool(tokens.at(i),Keyword)) StringstreamOUTPUT << tokens.at(i).substr(tokens.at(i).find(Keyword)+Keyword.length()) << endl;
      if(aurostd::substring2bool(tokens.at(i),Keyword)) status=TRUE;
    }
    return status;  // return FALSE if something got messed up
  }

  // *******************************************************************************************
  // *******************************************************************************************
  bool ExtractLastToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword) {   // AFLOW_FUNCTION_IMPLEMENTATION
    return ExtractToStringstreamEXPLICIT(FileIN,StringstreamOUTPUT,Keyword);
  }

  // *******************************************************************************************
  bool ExtractLastToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop) { // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "LDEBUG: ExtractLastToStringstreamEXPLICIT" << endl;
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    vector<string> tokens;
    aurostd::stream2vectorstring(FileIN,tokens);
    int istart=-1,istop=-1;
    for(int i=(int) tokens.size()-1;i>0;i--) {
      if(aurostd::substring2bool(tokens.at(i),Keyword_stop)  && istop<1)  istop=i-1;
      if(aurostd::substring2bool(tokens.at(i),Keyword_start) && istart<1) istart=i+1;
    }
    if(LDEBUG) cerr << "LDEBUG: " << istart << " " << istop << endl;
    if(istart>0 && istop>0) {
      for(int i=istart;i<=istop;i++) StringstreamOUTPUT << tokens.at(i) << endl;
      if(LDEBUG) cerr << "LDEBUG: " << StringstreamOUTPUT.str() << endl;
      return TRUE;
    }
    return FALSE;
  }


  bool ExtractLastToStringstreamEXPLICIT(stringstream StringStreamIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop) { // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    string StringIN=StringStreamIN.str();
    return ExtractLastToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start,Keyword_stop);
  }

  bool ExtractLastToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword_start,string Keyword_stop) {  // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "LDEBUG: ExtractLastToStringstreamEXPLICIT" << endl;
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    vector<string> tokens;
    aurostd::string2vectorstring(StringIN,tokens);
    int istart=-1,istop=-1;
    for(int i=(int) tokens.size()-1;i>0;i--) {
      if(aurostd::substring2bool(tokens.at(i),Keyword_stop)  && istop<1)  istop=i-1;
      if(aurostd::substring2bool(tokens.at(i),Keyword_start) && istart<1) istart=i+1;
    }
    if(LDEBUG) cerr << "LDEBUG: " << istart << " " << istop << endl;
    if(istart>0 && istop>0) {
      for(int i=istart;i<=istop;i++) StringstreamOUTPUT << tokens.at(i) << endl;
      if(LDEBUG) cerr << "LDEBUG: " << StringstreamOUTPUT.str() << endl;
      return TRUE;
    }
    return FALSE;
  }

  // *******************************************************************************************
  // *******************************************************************************************

  bool ExtractJustAfterToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,string Keyword_start) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ofstream FileOUTPUT;
    FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
    string strline,subS2;
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline)) {
      if(status) FileOUTPUT << strline << endl;
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
    return status;  // return FALSE if something got messed up
  }

  bool ExtractJustAfterToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,string Keyword_start) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string strline,subS2;
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline)) {
      if(status) StringOUTPUT=StringOUTPUT+strline+"\n";
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status;  // return FALSE if something got messed up
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,string Keyword_start) {    // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    string strline;
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline)) {
      if(status) StringstreamOUTPUT << strline << endl;
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status;  // return FALSE if something got messed up
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(stringstream StringStreamIN,stringstream& StringstreamOUTPUT,string Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    string StringIN=StringStreamIN.str();
    return ExtractJustAfterToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start);
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(string StringIN,stringstream& StringstreamOUTPUT,string Keyword_start) {  // AFLOW_FUNCTION_IMPLEMENTATION
    StringstreamOUTPUT.clear();StringstreamOUTPUT.str(std::string());
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      if(status) StringstreamOUTPUT << tokens.at(i) << endl;
      if(aurostd::substring2bool(tokens.at(i),Keyword_start)) status=TRUE;
    }
    return status;  // return FALSE if something got messed up
  }

  bool ExtractJustAfterToStringEXPLICIT(string StringIN,string& StringOUTPUT,string Keyword_start) {  // AFLOW_FUNCTION_IMPLEMENTATION
    stringstream StringstreamOUTPUT;
    bool out=ExtractJustAfterToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start);
    StringOUTPUT=StringstreamOUTPUT.str();
    return out;  // return FALSE if something got messed up
  }

  // ***************************************************************************
  // Function stream2vectorstring return UINT
  // ***************************************************************************
  // take istream into a vector strings - Stefano Curtarolo
  uint stream2vectorstring(std::istream& istreamIN,vector<string> &vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  // take ifstream into a vector strings - Stefano Curtarolo
  uint stream2vectorstring(std::ifstream& ifstreamIN,vector<string> &vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  } 
  uint stream2vectorstring(std::stringstream& stringstreamIN,vector<string> &vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  uint string2vectorstring(const string& stringIN,vector<string> &vstringout,bool consecutive,bool trim_edges) {  //CO 170613
    //CO mods 170613
    //we are adding functionality here, because string2tokens will treat "\n\n" same as "\n", but not "\n \n"
    //consecutive will do the following:  "sssss" -> <"s","s",...>
    //trim_edges will remove delimiters from beginning and end, similar to consecutive=false behavior
    //return aurostd::string2tokens(stringIN,vstringout,"\n",true);
    uint count=aurostd::string2tokens(stringIN,vstringout,"\n",consecutive);
    if(trim_edges){
      //start with front
      while(vstringout.size()){
        if(!aurostd::RemoveWhiteSpaces(vstringout.front()).empty()){
          break;
        }
        vstringout.erase(vstringout.begin());
        count--;
      }
      //now back
      while(vstringout.size()){
        if(!aurostd::RemoveWhiteSpaces(vstringout.back()).empty()){
          break;
        }
        vstringout.pop_back();
        count--;
      }
    }
    return count;
  }

  // ***************************************************************************
  // Function string2vectorstring return VECTOR
  // ***************************************************************************
  // take sitring into a vector strings - Stefano Curtarolo
  vector<string> stream2vectorstring(std::istream& istreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(istreamIN,vstringout);
    return vstringout;
  }
  vector<string> stream2vectorstring(std::ifstream& iftreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(iftreamIN,vstringout);
    return vstringout;
  }
  vector<string> stream2vectorstring(std::stringstream& stringstreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(stringstreamIN,vstringout);
    return vstringout;
  }
  vector<string> string2vectorstring(const string& stringIN,bool consecutive,bool trim_edges) { //CO 170613
    vector<string> vstringout;
    aurostd::string2vectorstring(stringIN,vstringout,consecutive,trim_edges); //CO 170613
    return vstringout;
  }

  // ***************************************************************************
  // Function liststring2string return string
  // ***************************************************************************
  string liststring2string(string s00,string s01,string s02,string s03,string s04,string s05,string s06,string s07,
			   string s08,string s09,string s0A,string s0B,string s0C,string s0D,string s0E,string s0F,
			   string s10,string s11,string s12,string s13,string s14,string s15,string s16,string s17,
			   string s18,string s19,string s1A,string s1B,string s1C,string s1D,string s1E,string s1F,
			   string s20,string s21,string s22,string s23,string s24,string s25,string s26,string s27,
			   string s28,string s29,string s2A,string s2B,string s2C,string s2D,string s2E,string s2F) {
    string out="";
    if(s00!="") out+=s00+"\n";
    if(s01!="") out+=s01+"\n";
    if(s02!="") out+=s02+"\n";
    if(s03!="") out+=s03+"\n";
    if(s04!="") out+=s04+"\n";
    if(s05!="") out+=s05+"\n";
    if(s06!="") out+=s06+"\n";
    if(s07!="") out+=s07+"\n";
    if(s08!="") out+=s08+"\n";
    if(s09!="") out+=s09+"\n";
    if(s0A!="") out+=s0A+"\n";
    if(s0B!="") out+=s0B+"\n";
    if(s0C!="") out+=s0C+"\n";
    if(s0D!="") out+=s0D+"\n";
    if(s0E!="") out+=s0E+"\n";
    if(s0F!="") out+=s0F+"\n";
    if(s10!="") out+=s10+"\n";
    if(s11!="") out+=s11+"\n";
    if(s12!="") out+=s12+"\n";
    if(s13!="") out+=s13+"\n";
    if(s14!="") out+=s14+"\n";
    if(s15!="") out+=s15+"\n";
    if(s16!="") out+=s16+"\n";
    if(s17!="") out+=s17+"\n";
    if(s18!="") out+=s18+"\n";
    if(s19!="") out+=s19+"\n";
    if(s1A!="") out+=s1A+"\n";
    if(s1B!="") out+=s1B+"\n";
    if(s1C!="") out+=s1C+"\n";
    if(s1D!="") out+=s1D+"\n";
    if(s1E!="") out+=s1E+"\n";
    if(s1F!="") out+=s1F+"\n";
    if(s20!="") out+=s20+"\n";
    if(s21!="") out+=s21+"\n";
    if(s22!="") out+=s22+"\n";
    if(s23!="") out+=s23+"\n";
    if(s24!="") out+=s24+"\n";
    if(s25!="") out+=s25+"\n";
    if(s26!="") out+=s26+"\n";
    if(s27!="") out+=s27+"\n";
    if(s28!="") out+=s28+"\n";
    if(s29!="") out+=s29+"\n";
    if(s2A!="") out+=s2A+"\n";
    if(s2B!="") out+=s2B+"\n";
    if(s2C!="") out+=s2C+"\n";
    if(s2D!="") out+=s2D+"\n";
    if(s2E!="") out+=s2E+"\n";
    if(s2F!="") out+=s2F+"\n";
    return out;
  }
  
  // ***************************************************************************
  // Function stream2dequestring return UINT
  // ***************************************************************************
  // take istream into a deque strings - Stefano Curtarolo
  uint stream2dequestring(std::istream& istreamIN,deque<string> &vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  // take ifstream into a deque strings - Stefano Curtarolo
  uint stream2dequestring(std::ifstream& ifstreamIN,deque<string> &vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  uint stream2dequestring(std::stringstream& stringstreamIN,deque<string> &vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  uint string2dequestring(const string& stringIN,deque<string> &vstringout) {
    return aurostd::string2tokens(stringIN,vstringout,"\n");
  }

  // ***************************************************************************
  // Function string2dequestring return DEQUE
  // ***************************************************************************
  // take sitring into a deque strings - Stefano Curtarolo
  deque<string> stream2dequestring(std::istream& istreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(istreamIN,vstringout);
    return vstringout;
  }
  deque<string> stream2dequestring(std::ifstream& iftreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(iftreamIN,vstringout);
    return vstringout;
  }
  deque<string> stream2dequestring(std::stringstream& stringstreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(stringstreamIN,vstringout);
    return vstringout;
  }
  deque<string> string2dequestring(const string& stringIN) {
    deque<string> vstringout;
    aurostd::string2dequestring(stringIN,vstringout);
    return vstringout;
  }

  // ***************************************************************************
  // Function string2file string2gzfile string2bz2file
  // ***************************************************************************
  // write string to file - Stefano Curtarolo
  bool string2file(const string& StringOUTPUT,string FileNameOUTPUT,string mode) {
    if(mode=="POST" || mode=="APPEND") {
      stringstream FileINPUT;
      aurostd::file2stringstream(FileNameOUTPUT,FileINPUT);
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      FileOUTPUT << FileINPUT.str();
      FileOUTPUT << StringOUTPUT;    
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return TRUE;  // return FALSE if something got messed up
    }
    if(mode=="PRE") {
      stringstream FileINPUT;
      aurostd::file2stringstream(FileNameOUTPUT,FileINPUT);
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      FileOUTPUT << StringOUTPUT;    
      FileOUTPUT << FileINPUT.str();
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return TRUE;  // return FALSE if something got messed up
    }
    if(mode=="WRITE" || mode=="") {
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      FileOUTPUT << StringOUTPUT;    
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return TRUE;  // return FALSE if something got messed up
    }   
    return FALSE;
  }

  bool string2gzfile(const string& StringOUTPUT,string FileNameOUTPUT,string mode) {
    bool out=string2file(StringOUTPUT,FileNameOUTPUT,mode);
    aurostd::execute("gzip -9fq "+FileNameOUTPUT);
    return out;
  }

  bool string2bz2file(const string& StringOUTPUT,string FileNameOUTPUT,string mode) {
    bool out=string2file(StringOUTPUT,FileNameOUTPUT,mode);
    aurostd::execute("bzip2 -9fq "+FileNameOUTPUT);
    return out;
  }

  // ***************************************************************************
  // Function stringstream2file stringstream2gzfile stringstream2bz2file
  // ***************************************************************************
  // write string to file - Stefano Curtarolo
  bool stringstream2file(const stringstream& StringstreamOUTPUT,string FileNameOUTPUT,string mode) {
    if(mode=="POST" || mode=="APPEND") {
      stringstream FileINPUT;
      aurostd::file2stringstream(FileNameOUTPUT,FileINPUT);
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      FileOUTPUT << FileINPUT.str();
      FileOUTPUT << StringstreamOUTPUT.str();
      // FileOUTPUT << StringstreamOUTPUT.rdbuf();
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return TRUE;  // return FALSE if something got messed up
    }
    if(mode=="PRE") {
      stringstream FileINPUT;
      aurostd::file2stringstream(FileNameOUTPUT,FileINPUT);
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      FileOUTPUT << StringstreamOUTPUT.str();
      // FileOUTPUT << StringstreamOUTPUT.rdbuf();
      FileOUTPUT << FileINPUT.str();
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return TRUE;  // return FALSE if something got messed up
    }
    if(mode=="WRITE" || mode=="") {
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      FileOUTPUT << StringstreamOUTPUT.str();
      // FileOUTPUT << StringstreamOUTPUT.rdbuf();
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return TRUE;  // return FALSE if something got messed up
    }   
    return FALSE;
  }

  bool stringstream2gzfile(const stringstream& StringstreamOUTPUT,string FileNameOUTPUT,string mode) {
    bool out=stringstream2file(StringstreamOUTPUT,FileNameOUTPUT,mode);
    aurostd::execute("gzip -9fq "+FileNameOUTPUT);
    return out;
  }

  bool stringstream2bz2file(const stringstream& StringstreamOUTPUT,string FileNameOUTPUT,string mode) {
    bool out=stringstream2file(StringstreamOUTPUT,FileNameOUTPUT,mode);
    aurostd::execute("bzip2 -9fq "+FileNameOUTPUT);
    return out;
  }
  
  // ***************************************************************************
  // Function ostream2string  istream2string   
  // ***************************************************************************
  // convert ostream/istream to string - Stefano Curtarolo
  std::string ostream2string(std::ostream& oss) {
    std::stringstream soss;
    soss << oss.rdbuf();
    return soss.str();
  }

  // ***************************************************************************
  // Function stream2string
  // ***************************************************************************
  // take istream into a  strings - Stefano Curtarolo
  uint stream2string(std::istream& istreamIN,string &vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout+=string(tmp)+"/n";
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length();  // return FALSE if something got messed up
  }

  // take ifstream into a  strings - Stefano Curtarolo
  uint stream2string(std::ifstream& ifstreamIN,string &vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout+=string(tmp)+"/n";
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length();  // return FALSE if something got messed up
  }

  uint stream2string(std::stringstream& stringstreamIN,string &vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout+=string(tmp)+"/n";
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length();  // return FALSE if something got messed up
  }

  // ***************************************************************************
  // Function getenv2string getenv2int getenv2uint getenv2double
  // ***************************************************************************
  // convert environments to string;
  string getenv2string(const string& str) {
    if(getenv(str.c_str())==NULL) return string("");
    return string(getenv(str.c_str()));
  }
  int getenv2int(const string& str) {
    if(getenv(str.c_str())==NULL) return int(0);
    return aurostd::string2utype<int>(getenv(str.c_str()));
  }
  uint getenv2uint(const string& str) {
    if(getenv(str.c_str())==NULL) return uint(0);
    return aurostd::string2utype<uint>(getenv(str.c_str()));
  }
  double getenv2double(const string& str) {
    if(getenv(str.c_str())==NULL) return double(0);
    return aurostd::string2utype<double>(getenv(str.c_str()));
  }

  // ***************************************************************************
  // Function file2string bz2file2string gzfile2string zipfile2string efile2string
  // ***************************************************************************
  // write file to string - Stefano Curtarolo
  uint file2string(string _FileNameIN,string& StringIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - file2string:  file=" << FileNameIN << " not present !" << endl;
      return 0;}   
    ifstream FileIN;
    FileIN.open(FileNameIN.c_str(),std::ios::in);
    char c; while (FileIN.get(c)) StringIN+=c;
    FileIN.clear();FileIN.close();
    return StringIN.length();  // return 0 if something got messed up
  }
  uint bz2file2string(string _FileNameIN,string& StringIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    // cerr << "bz2file2string; BEGIN" << endl;
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - bz2file2string:  file=" << FileNameIN << " not present !" << endl;
      return 0;}   
    if(!aurostd::IsCommandAvailable("bzcat")) {
      // cerr << "ERROR - bz2file2string:  command \"bzcat\" is necessary !" << endl;
      return 0;}   
    StringIN=aurostd::execute2string("bzcat "+FileNameIN);
    return StringIN.length();  // return 0 if something got messed up
  }
  uint gzfile2string(string _FileNameIN,string& StringIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    // cerr << "gzfile2string; BEGIN" << endl;
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - gzfile2string:  file=" << FileNameIN << " not present !" << endl;
      return 0;}   
    if(!aurostd::IsCommandAvailable("zcat")) {
      // cerr << "ERROR - gzfile2string:  command \"zcat\" is necessary !" << endl;
      return 0;}   
    StringIN=aurostd::execute2string("zcat "+FileNameIN);
    return StringIN.length();  // return 0 if something got messed up
  }
  //CO - START
  uint zipfile2string(string _FileNameIN,string& StringIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      return 0;}
    if(!aurostd::IsCommandAvailable("unzip")) {
      return 0;}
    StringIN=aurostd::execute2string("unzip -p "+FileNameIN);
    return StringIN.length();  // return 0 if something got messed up
  }
  //CO - END
  uint efile2string(string _FileNameIN,string& StringIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN),FileNameOUT;  //corey
    // cerr << "efile2string; BEGIN FileNameIN=[" << FileNameIN << "]" << endl;
    if(!EFileExist(FileNameIN,FileNameOUT)) {
      // cerr << "ERROR - efile2string:  file=" << FileNameIN << " not present !" << endl;
      return 0;}   
    // cerr << aurostd::substring2bool(FileNameIN,".bz2") << endl;
    if(aurostd::substring2bool(FileNameOUT,".bz2")) {
      // cerr << "efile2string: found .bz2, using bz2file2string" << endl;
      return aurostd::bz2file2string(FileNameOUT,StringIN);
    }
    if(aurostd::substring2bool(FileNameOUT,".gz"))  {
      // cerr << "efile2string: found .gz, using gzfile2string" << endl;
      return aurostd::gzfile2string(FileNameOUT,StringIN);
    }
    // cerr << "efile2string: file2string" << endl;
    return aurostd::file2string(FileNameOUT,StringIN);
  }
  
  // ***************************************************************************
  // Function file2vectorstring bz2file2vectorstring gzfile2vectorstring efile2vectorstring 
  // ***************************************************************************
  // write file to vector string - Stefano Curtarolo
  uint file2vectorstring(string FileNameIN,vector<string>& vlines) {
    return aurostd::string2vectorstring(file2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  uint bz2file2vectorstring(string FileNameIN,vector<string>& vlines) {
    return aurostd::string2vectorstring(bz2file2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  uint gzfile2vectorstring(string FileNameIN,vector<string>& vlines) {
    return aurostd::string2vectorstring(gzfile2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  uint efile2vectorstring(string FileNameIN,vector<string>& vlines) {
    return aurostd::string2vectorstring(efile2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  // ***************************************************************************
  // Function file2dequestring bz2file2dequestring gzfile2dequestring efile2dequestring 
  // ***************************************************************************
  // write file to deque string - Stefano Curtarolo
  uint file2dequestring(string FileNameIN,deque<string>& vlines) {
    return aurostd::string2dequestring(file2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  uint bz2file2dequestring(string FileNameIN,deque<string>& vlines) {
    return aurostd::string2dequestring(bz2file2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  uint gzfile2dequestring(string FileNameIN,deque<string>& vlines) {
    return aurostd::string2dequestring(gzfile2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  uint efile2dequestring(string FileNameIN,deque<string>& vlines) {
    return aurostd::string2dequestring(efile2string(aurostd::CleanFileName(FileNameIN)),vlines);
  }

  // ***************************************************************************
  // Function file2stringstream bz2file2stringstream gzfile2stringstream zipfile2stringstream efile2stringstream
  // ***************************************************************************
  // write file to stringstream - Stefano Curtarolo
  bool file2stringstream(string _FileNameIN,stringstream& StringstreamIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - file2stringstream:  file=" << FileNameIN << " not present !" << endl;
      return FALSE;}   
    ifstream FileIN;
    StringstreamIN.clear();StringstreamIN.str(std::string());
    FileIN.open(FileNameIN.c_str(),std::ios::in);
    char c; while (FileIN.get(c)) StringstreamIN.put(c);
    FileIN.clear();FileIN.close();
    return TRUE;  // return FALSE if something got messed up
  }
  bool bz2file2stringstream(string _FileNameIN,stringstream& StringstreamIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - bz2file2stringstream:  file=" << FileNameIN << " not present !" << endl;
      return FALSE;}   
    if(!aurostd::IsCommandAvailable("bzcat")) {
      cerr << "ERROR - bz2file2stringstream:  command \"bzcat\" is necessary !" << endl;
      return FALSE;}   
    StringstreamIN.clear();StringstreamIN.str(std::string());
    StringstreamIN << execute2string("bzcat "+FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  bool gzfile2stringstream(string _FileNameIN,stringstream& StringstreamIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - gzfile2stringstream:  file=" << FileNameIN << " not present !" << endl;
      return FALSE;}   
    if(!aurostd::IsCommandAvailable("zcat")) {
      cerr << "ERROR - gzfile2stringstream:  command \"zcat\" is necessary !" << endl;
      return FALSE;}   
    StringstreamIN.clear();StringstreamIN.str(std::string());
    StringstreamIN << execute2string("zcat "+FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  //CO - START
  //zipfile 2 a string
  bool zipfile2stringstream(string _FileNameIN,stringstream& StringstreamIN) {
    string FileNameIN=CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - zipfile2stringstream:  file=" << FileNameIN << " not present !" << endl;
      return FALSE;}
    if(!aurostd::IsCommandAvailable("zcat")) {
      cerr << "ERROR - zipfile2stringstream:  command \"unzip\" is necessary !" << endl;
      return FALSE;}
    StringstreamIN.clear();StringstreamIN.str(std::string());
    StringstreamIN << execute2string("unzip -p "+FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  //CO - END
  bool efile2stringstream(string _FileNameIN,stringstream& StringstreamIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN),FileNameOUT;  //corey
    if(!EFileExist(FileNameIN,FileNameOUT)) {
      cerr << "ERROR - efile2stringstream:  file=" << FileNameIN << " not present !" << endl;
      return FALSE;}   
    if(aurostd::substring2bool(FileNameOUT,".bz2")) return aurostd::bz2file2stringstream(FileNameOUT,StringstreamIN);
    if(aurostd::substring2bool(FileNameOUT,".gz"))  return aurostd::gzfile2stringstream(FileNameOUT,StringstreamIN);
    return aurostd::file2stringstream(FileNameOUT,StringstreamIN);
  }
  
  // ***************************************************************************
  // Function file2string
  // ***************************************************************************
  // write file to string - Stefano Curtarolo
  string file2string(string FileNameIN) {
    ifstream FileIN;
    string StringIN;
    FileIN.open(FileNameIN.c_str(),std::ios::in);
    char c; while (FileIN.get(c)) StringIN+=c;
    FileIN.clear();FileIN.close();
    return StringIN;
  }
  string bz2file2string(string FileNameIN) {
    string StringIN;
    bz2file2string(FileNameIN,StringIN);
    return StringIN;
  }
  string gzfile2string(string FileNameIN) {
    string StringIN;
    gzfile2string(FileNameIN,StringIN);
    return StringIN;
  }
  string efile2string(string FileNameIN) {
    string StringIN;
    efile2string(FileNameIN,StringIN);
    return StringIN;
  }

  // ***************************************************************************
  // Function url2file
  // ***************************************************************************
  // wget URL to string - Stefano Curtarolo
  bool url2file(string url,string& fileIN,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(!aurostd::IsCommandAvailable("wget")) {
      cerr << "ERROR - aurostd::url2file(): command \"wget\" is necessary !" << endl;
      return FALSE;}	
    string _url=url;
    aurostd::StringSubst(_url,"http://","");
    aurostd::StringSubst(_url,"//","/");
    if(LDEBUG) cerr << "aurostd::url2file(): Loading url=" << _url << endl;
    if(verbose) cout << "aurostd::url2file(): Loading url=" << _url << endl;
#ifndef _MACOSX_
    aurostd::execute("wget --quiet --no-cache -O "+fileIN+" http://"+_url);
#else
    aurostd::execute("wget --quiet -O "+fileIN+" http://"+_url);
#endif    
    if(aurostd::FileEmpty(fileIN)) {
      aurostd::StringSubst(_url,":AFLOW","/AFLOW");
#ifndef _MACOSX_
      aurostd::execute("wget --quiet --no-cache -O "+fileIN+" http://"+_url);
#else
      aurostd::execute("wget --quiet -O "+fileIN+" http://"+_url); // _MACOSX_
#endif    
      if(aurostd::FileEmpty(fileIN)) {
	cerr << "ERROR - aurostd::url2file(): URL not found http://" << _url << endl;
	return FALSE;
      }
    }
    return TRUE;
  }
  
  // ***************************************************************************
  // Function url2string
  // ***************************************************************************
  // wget URL to string - Stefano Curtarolo
  bool url2string(string url,string& stringIN,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(!aurostd::IsCommandAvailable("wget")) {
      cerr << "ERROR - aurostd::url2string(): command \"wget\" is necessary !" << endl;
      return FALSE;}	
    string _url=url;
    aurostd::StringSubst(_url,"http://","");
    aurostd::StringSubst(_url,"//","/");
    if(LDEBUG) cerr << "aurostd::url2string(): Loading url=" << _url << endl;
    if(verbose) cout << "aurostd::url2string(): Loading url=" << _url << endl;
#ifndef _MACOSX_
    stringIN=aurostd::execute2string("wget --quiet --no-cache -O /dev/stdout http://"+_url);
#else
    stringIN=aurostd::execute2string("wget --quiet -O /dev/stdout http://"+_url); // _MACOSX_
#endif    
    if(stringIN=="") {
      aurostd::StringSubst(_url,":AFLOW","/AFLOW");
#ifndef _MACOSX_
      stringIN=aurostd::execute2string("wget --quiet --no-cache -O /dev/stdout http://"+_url);
#else
      stringIN=aurostd::execute2string("wget --quiet -O /dev/stdout http://"+_url); // _MACOSX_
#endif    
      if(stringIN=="") {
	cerr << "ERROR - aurostd::url2string(): URL not found http://" << _url << endl;
	return FALSE;
      }
    }
    //    aurostd::StringSubst(stringIN,"h1h","h1,"); // old Frisco php error patch
    return TRUE;
  }
  
  // ***************************************************************************
  // Function url2stringstream
  // ***************************************************************************
  // wget URL to stringstream - Stefano Curtarolo
  bool url2stringstream(string url,stringstream& stringstreamIN,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aurostd::url2stringstream() Loading url=" << url << endl;
    if(verbose) cout << "aurostd::url2stringstream() Loading url=" << url << endl;
    string stringIN;
    bool out=url2string(url,stringIN,verbose);
    stringstreamIN.clear(); stringstreamIN << stringIN;
    return out;
  }

  // ***************************************************************************
  // Function url2vectorstring
  // ***************************************************************************
  // wget URL to vectorstring - Stefano Curtarolo
  bool url2vectorstring(string url,vector<string>& vlines,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aurostd::url2vectorstring() Loading url=" << url << endl;
    if(verbose) cout << "aurostd::url2vectorstring() Loading url=" << url << endl;
    string stringIN;
    bool out=url2string(url,stringIN,verbose);
    aurostd::string2tokens(stringIN,vlines);
    return out;
  }
  
  // ***************************************************************************
  // Function url2dequestring
  // ***************************************************************************
  // wget URL to dequestring - Stefano Curtarolo
  bool url2dequestring(string url,deque<string>& vlines,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aurostd::url2dequestring() Loading url=" << url << endl;
    if(verbose) cout << "aurostd::url2dequestring() Loading url=" << url << endl;
    string stringIN;
    bool out=url2string(url,stringIN,verbose);
    aurostd::string2tokens(stringIN,vlines);
    return out;
  }
  
  // ***************************************************************************
  // Function url2tokens
  // ***************************************************************************
  // wget URL to vector of tokens - Stefano Curtarolo
  template<typename utype> uint url2tokens(string url,vector<utype>& tokens,const string& delimiters) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "aurostd::url2tokens<utype>(vector) Loading url=" << url << endl;
    if(!aurostd::IsCommandAvailable("wget")) {
      cerr << "ERROR - aurostd::url2tokens<utype>(vector): command \"wget\" is necessary !" << endl;
      return 0;}	
    tokens.clear(); 
    string content;
    aurostd::url2string(url,content);
    if(content=="") {cerr << "ERROR - aurostd::url2tokens<utype>(vector): URL empty http://" << url << endl;return 0;}
    vector<string> stokens;
    aurostd::string2tokens(content,stokens,delimiters);
    for(uint i=0;i<stokens.size();i++)
      if(stokens.at(i)!="") 
	tokens.push_back(aurostd::string2utype<utype>(stokens.at(i)));
    if(LDEBUG) cerr << "aurostd::url2tokens() [5] tokens.size()=" << tokens.size() << endl;
    return tokens.size();
  }

  // ***************************************************************************
  // Function url2tokens
  // ***************************************************************************
  // wget URL to deque of tokens - Stefano Curtarolo
  template<typename utype> uint url2tokens(string url,deque<utype>& tokens,const string& delimiters) {
    vector<utype> vtokens;
    aurostd::url2tokens(url,vtokens,delimiters);
    for(uint i=0;i<vtokens.size();i++) tokens.push_back(vtokens.at(i));
    return tokens.size();
  }

  // ***************************************************************************
  // Function url2string
  // ***************************************************************************
  // wget URL to stringstream - Stefano Curtarolo
  string url2string(string url) {
    string stringIN;
    url2string(url,stringIN);
    return stringIN;
  }

  // ***************************************************************************
  // Function ChmodFile
  // ***************************************************************************
  // change mod of a file - Stefano Curtarolo
  bool ChmodFile(string chmod_string,string file) {
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << CHMOD_BIN << " " << chmod_string << " " << file << endl;
    aurostd::execute(aus);
    return TRUE;  // return FALSE if something got messed up
  }

  //CO - START
  //***************************************************************************//
  // aurostd::MoveFile
  //***************************************************************************//
  // Corey Oses
  // Move file to destination
  bool MoveFile(const string& _file,const string& destination) {
    string file=CleanFileName(_file);
    if(!aurostd::IsCommandAvailable("mv")) {
      cerr << "ERROR - MoveFile:  command \"mv\" is necessary !" << endl;
      return FALSE;
    }
    if(!aurostd::FileExist(file)) {
      cerr << "ERROR - MoveFile:  " << file << " cannot be found !" << endl;
      return FALSE;
    }
    if(!aurostd::IsDirectory(destination)){ //CO 180220 //FileExist(destination)) {
      cerr << "ERROR - MoveFile:  " << destination << " cannot be found !" << endl;
      return FALSE;
    }
    aurostd::execute("mv "+file+" "+destination+"/");
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::MoveFiles
  //***************************************************************************//
  // Corey Oses
  // Move files to destination, do one at a time to check if files exist
  bool MoveFiles(const vector<string>& files,const string& destination) {
    for(uint i=0;i<files.size();i++) {
      if(!MoveFile(files.at(i),destination)) {return FALSE;}
    }
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::RenameFile
  //***************************************************************************//
  // Corey Oses
  bool RenameFile(const string& _file,const string& destination) {
    string file=CleanFileName(_file);
    if(!aurostd::IsCommandAvailable("mv")) {
      cerr << "ERROR - RenameFile:  command \"mv\" is necessary !" << endl;
      return FALSE;
    }
    if(!aurostd::FileExist(file)) {
      cerr << "ERROR - RenameFile:  " << file << " cannot be found !" << endl;
      return FALSE;
    }
    aurostd::execute("mv "+file+" "+destination);
    return TRUE;
  }

  // ***************************************************************************
  // Function IsDirectory
  // ***************************************************************************
  // true if path is directory
  // http://stackoverflow.com/questions/146924/how-can-i-tell-if-a-given-path-is-a-directory-or-a-file-c-c
  bool IsDirectory(string path){
    struct stat s;
    if( stat(path.c_str(),&s) != 0 ){return FALSE;} //error
    if( s.st_mode & S_IFDIR ){return TRUE;}
    return FALSE;
  }

  // ***************************************************************************
  // Function IsFile
  // ***************************************************************************
  // true if path is file
  // http://stackoverflow.com/questions/146924/how-can-i-tell-if-a-given-path-is-a-directory-or-a-file-c-c
  bool IsFile(string path){
    struct stat s;
    if( stat(path.c_str(),&s) != 0 ){return FALSE;} //error
    if( s.st_mode & S_IFREG ){return TRUE;}
    return FALSE;
  }
  //CO - END

  // ***************************************************************************
  // Function RemoveFile
  // ***************************************************************************
  // change remove a file - Stefano Curtarolo
  bool RemoveFile(string file) {
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    remove(file.c_str());
    //    aus << "rm -f \"" << file << "\"" << endl; //   cerr << aus.str() << endl;
    //    aurostd::execute(aus);
    return TRUE;  // return FALSE if something got messed up
  }

  //CO - START
  //***************************************************************************//
  // aurostd::RemoveFiles
  //***************************************************************************//
  bool RemoveFiles(const vector<string>& files) {
    for(uint i=0;i<files.size();i++) {
      if(!RemoveFile(files.at(i))) {return FALSE;}
    }
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::RemoveDirectory
  //***************************************************************************//
  bool RemoveDirectory(const string& dir) {
    if(!IsCommandAvailable("rm")) {
      cerr << "ERROR - rm:  command \"rm\" is necessary !" << endl;
      return FALSE;
    }
    aurostd::execute("rm -r "+dir);
    return TRUE;
  }
  //CO - END

  // ***************************************************************************
  // Function string2tokens string2tokens<utype>
  // ***************************************************************************
  // Finds string2tokens to split strings in tokens
  // Stefano Curtarolo
  // void string2tokens(const string& str,vector<string>& tokens,const string& delimiters=" ") {
  uint string2tokens(const string& str,std::vector<string>& tokens,const string& delimiters,bool consecutive) { //CO 170613
    //CO mods 170613
    //we are adding functionality here, because string2tokens will treat "\n\n" same as "\n", but not "\n \n"
    //consecutive will do the following:  "sssss" -> <"s","s",...>
    //trim_edges will remove delimiters from beginning and end, similar to consecutive=false behavior
    //return aurostd::string2tokens(stringIN,vstringout,"\n",true);
    tokens.clear(); // clear in the case there was something already in!
    string::size_type lastPos=(consecutive ? 0 : str.find_first_not_of(delimiters,0));   // Skip delimiters at beginning.
    string::size_type pos=str.find_first_of(delimiters,lastPos);     // Find first "non-delimiter".
    while (string::npos!=pos || string::npos!=lastPos) {
      tokens.push_back(str.substr(lastPos,pos-lastPos));             // Found a token, add it to the vector.
      if(consecutive){lastPos=(string::npos!=pos ? pos+1 : string::npos);}
      else{lastPos=str.find_first_not_of(delimiters,pos);}  // Skip delimiters.  Note the "not_of"
      pos=str.find_first_of(delimiters,lastPos);                     // Find next "non-delimiter"
    }
    return tokens.size();
  }
  
  // void string2tokens(const string& str,deque<string>& tokens,const string& delimiters=" ") {
  uint string2tokens(const string& str,std::deque<string>& tokens,const string& delimiters,bool consecutive) { //CO 170613
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters,consecutive);  //CO 170613
    tokens.clear();
    for(i=0;i<vtokens.size();i++) tokens.push_back(vtokens.at(i));
    return tokens.size();
  }
  template<class utype> uint string2tokens(const string& str,std::vector<utype>& tokens,const string& delimiters,bool consecutive) {  //CO 170613
    vector<string> stokens;
    uint out=aurostd::string2tokens(str,stokens, delimiters,consecutive); //CO 170613
    tokens.clear();
    for(uint i=0;i<stokens.size();i++)
      tokens.push_back(aurostd::string2utype<utype>(stokens.at(i)));
    return out;
  }
  template<class utype> uint string2tokens(const string& str,std::deque<utype>& tokens,const string& delimiters,bool consecutive) { //CO 170613
    deque<string> stokens;
    uint out=aurostd::string2tokens(str,stokens, delimiters,consecutive); //CO 170613
    tokens.clear();
    for(uint i=0;i<stokens.size();i++)
      tokens.push_back(aurostd::string2utype<utype>(stokens.at(i)));
    return out;
  }
  
  // ***************************************************************************
  // Function string2tokensAdd string2tokensAdd<utype>
  // ***************************************************************************

  uint string2tokensAdd(const string& str,std::vector<string>& tokens,const string& delimiters) {
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(vtokens.at(i));
    return tokens.size();
  }
  uint string2tokensAdd(const string& str,std::deque<string>& tokens,const string& delimiters) {
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(vtokens.at(i));
    return tokens.size();
  }
  template<class utype> uint string2tokensAdd(const string& str,std::vector<utype>& tokens,const string& delimiters) {
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(aurostd::string2utype<utype>(vtokens.at(i)));
    return tokens.size();
  }
  template<class utype> uint string2tokensAdd(const string& str,std::deque<utype>& tokens,const string& delimiters) {
    deque<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(aurostd::string2utype<utype>(vtokens.at(i)));
    return tokens.size();
  }
  
  // ***************************************************************************
  // Function StringStreamConvert
  // ***************************************************************************
  // convert whatever into a string !
  template<typename typeTo, typename typeFrom> typeTo StringStreamConvert(typeFrom from) {
    std::stringstream temp;
    temp << from;
    typeTo to=typeTo();
    temp >> to;
    return to;
  }

  template<typename typeFrom> std::string StringConvert(typeFrom from) {
    return StringStreamConvert<std::string>((typeFrom) from);
  }
  // to initialize the templates....
  void _StringConvert(void) {
    cerr << StringStreamConvert<std::string>((int) 1);
    cerr << StringStreamConvert<std::string>((float) 1);
    cerr << StringStreamConvert<std::string>((double) 1);
    cerr << StringConvert<int>((int) 1);
    cerr << StringConvert<float>((float) 1);
    cerr << StringConvert<double>((double) 1);
  }

  // ***************************************************************************
  // Function stream2stream
  // ***************************************************************************
  // convert whatever into a string !
  template<typename typeTo, typename typeFrom> typeTo stream2stream(typeFrom from,int precision,char FORMAT) {
    std::stringstream temp;
    if(FORMAT==DEFAULT_STREAM){;} //default
    if(FORMAT==FIXED_STREAM){temp << std::fixed;}
    if(FORMAT==SCIENTIFIC_STREAM){temp << std::scientific;}
    temp.precision(precision);
    temp << from;
    typeTo to=typeTo();
    temp >> to;
    return to;
  }
  template<typename typeTo, typename typeFrom> typeTo stream2stream(typeFrom from,int precision) {
    return (typeTo) stream2stream<typeTo>(from,precision,DEFAULT_STREAM);
  }
  template<typename typeTo, typename typeFrom> typeTo stream2stream(typeFrom from) {
    return (typeTo) stream2stream<typeTo>(from,AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);
  }
  template<typename utype> utype string2utype(const string& from) {
    return (utype) stream2stream<utype>(from,AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);
  }

  vector<double> vectorstring2vectordouble(vector<string> from) {
    vector<double> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<double>(from.at(i)));return vout;
  }

  string string2string(const string& from) {
    return from;
  }

  vector<int> vectorstring2vectorint(vector<string> from) {
    vector<int> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<int>(from.at(i)));return vout;
  }

  vector<uint> vectorstring2vectoruint(vector<string> from) {
    vector<uint> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<uint>(from.at(i)));return vout;
  }
  
  vector<float> vectorstring2vectorfloat(vector<string> from) {
    vector<float> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<float>(from.at(i)));return vout;
  }

  string vectorstring2string(const vector<string>& vstrings) {
    string out="";
    for(uint istr=0;istr<vstrings.size();istr++) out+=vstrings.at(istr);
    return out;
  }
  string vectorstring2string(const deque<string>& vstrings) {
    string out="";
    for(uint istr=0;istr<vstrings.size();istr++) out+=vstrings.at(istr);
    return out;
  }
  
  // ***************************************************************************
  // Function utype2string
  // ***************************************************************************
  
  template<typename utype> string utype2string(const utype& from) {
    return (string) stream2stream<string>(from);
  }
  template<typename utype> string utype2string(const utype& from,int precision) {
    return (string) stream2stream<string>(from,precision);
  }
  template<typename utype> string utype2string(const utype& from,int precision,char FORMAT) { //see DEFAULT_STREAM, FIXED_STREAM, SCIENTIFIC_STREAM
    return (string) stream2stream<string>(from,precision,FORMAT);
  }

  //cannot template this the same as others, char's don't make sense with roff
  string utype2string(double from,bool roff) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,DEFAULT_STREAM);}
  string utype2string(double from,int precision,bool roff) {return utype2string(from,precision,roff,AUROSTD_ROUNDOFF_TOL,DEFAULT_STREAM);}
  string utype2string(double from,bool roff,double tol) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,tol,DEFAULT_STREAM);}
  string utype2string(double from,int precision,bool roff,double tol) {return utype2string(from,precision,roff,tol,DEFAULT_STREAM);}
  string utype2string(double from,bool roff,char FORMAT) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,FORMAT);}
  string utype2string(double from,int precision,bool roff,char FORMAT) {return utype2string(from,precision,roff,AUROSTD_ROUNDOFF_TOL,FORMAT);}
  string utype2string(double from,bool roff,double tol,char FORMAT) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,tol,FORMAT);}
  string utype2string(double from,int precision,bool roff,double tol,char FORMAT) {
      double tmp=from;
      if(roff){tmp=roundoff(from,tol);}
      return (string) stream2stream<string>(tmp,precision,FORMAT);
  }

  // ***************************************************************************
  // Function utypes2deque
  // ***************************************************************************
  template<class utype> deque<utype> utypes2deque(utype u1) {
    deque<utype> out;
    out.push_back(u1);
    return out;}
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2) {
    deque<utype> out;
    out.push_back(u1);out.push_back(u2);
    return out;}
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3) {
    deque<utype> out;
    out.push_back(u1);out.push_back(u2);out.push_back(u3);
    return out;}
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3,utype u4) {
    deque<utype> out;
    out.push_back(u1);out.push_back(u2);out.push_back(u3);out.push_back(u4);
    return out;}

  // ***************************************************************************
  // Function StringCommasColumsVectorInt
  // ***************************************************************************
  void StringCommasColumsVectorInt(string vstring,vector<int> &vint) {
    vector<string> tokens_commas,tokens_colums;
    vint.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas.at(i),":")) {
	string2tokens(tokens_commas.at(i),tokens_colums,":");
	for(int j=aurostd::string2utype<int>(tokens_colums.at(0));
	    j<=aurostd::string2utype<int>(tokens_colums.at(tokens_colums.size()-1));
	    j++)
	  vint.push_back(j);
      } else {
	vint.push_back(aurostd::string2utype<int>(tokens_commas.at(i)));
      }	
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorUnsignedInt
  // ***************************************************************************
  void StringCommasColumsVectorUnsignedInt(string vstring,vector<uint> &vuint) {
    vector<string> tokens_commas,tokens_colums;
    vuint.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas.at(i),":")) {
	string2tokens(tokens_commas.at(i),tokens_colums,":");
	for(uint j=aurostd::string2utype<uint>(tokens_colums.at(0));
	    j<=(uint) aurostd::string2utype<uint>(tokens_colums.at(tokens_colums.size()-1));
	    j++)
	  vuint.push_back(j);
      } else {
	vuint.push_back(aurostd::string2utype<uint>(tokens_commas.at(i)));
      }	
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorFloat
  // ***************************************************************************
  void StringCommasColumsVectorFloat(string vstring,vector<float> &vfloat) {
    vector<string> tokens_commas,tokens_colums;
    vfloat.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas.at(i),":")) {
	string2tokens(tokens_commas.at(i),tokens_colums,":");
	for(float j=aurostd::string2utype<float>(tokens_colums.at(0));
	    j<=aurostd::string2utype<float>(tokens_colums.at(tokens_colums.size()-1));
	    j++)
	  vfloat.push_back(j);
      } else {
	vfloat.push_back(aurostd::string2utype<float>(tokens_commas.at(i)));
      }	
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorDouble
  // ***************************************************************************
  void StringCommasColumsVectorDouble(string vstring,vector<double> &vdouble) {
    vector<string> tokens_commas,tokens_colums;
    vdouble.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas.at(i),":")) {
	string2tokens(tokens_commas.at(i),tokens_colums,":");
	for(double j=aurostd::string2utype<double>(tokens_colums.at(0));
	    j<=aurostd::string2utype<double>(tokens_colums.at(tokens_colums.size()-1));
	    j++)
	  vdouble.push_back(j);
      } else {
	vdouble.push_back(aurostd::string2utype<double>(tokens_commas.at(i)));
      }	
    }
  }

  // ***************************************************************************
  // Function StringsAlphabetic
  // ***************************************************************************
  // says if two strings are in alphabetical order
  bool StringsAlphabetic(string A,string B) {
    if(A<B) return TRUE; // cerr << "A<B" << "  " << A << "<" << B << endl;
    if(A>B) return FALSE; //cerr << "A>B" << "  " << A << ">" << B << endl;
    if(A==B) return TRUE; //cerr << "A==B" << "  " << A << "==" << B << endl;
    return TRUE;
  }

  // ***************************************************************************
  // Function StringSubst
  // ***************************************************************************
  // Substitute strings here and there
  // Stefano Curtarolo
  string StringSubst(string &strstring, const string &strfind, const string &strreplace) {
    if(strfind.empty()) return strstring;
    string::size_type pos=0;
    while((pos=strstring.find(strfind, pos))!=string::npos) {
      strstring.erase(pos, strfind.length());
      strstring.insert(pos, strreplace);
      pos+=strreplace.length();
    }
    return strstring;
  }

  string StringSubst(string &strstring, const char &charfind, const char &charreplace) {
    string stroutput;
    for (uint i=0;i<strstring.size();i++)
      if(strstring[i]==charfind)
	stroutput+=charreplace;
      else
	stroutput+=strstring[i];
    strstring=stroutput;
    return strstring;
  }
  
  void StringStreamSubst(stringstream &strstringsteram, const string &strfind, const string &strreplace) {
    string strstring=strstringsteram.str();
    StringSubst(strstring,strfind,strreplace);
    strstringsteram.clear();
    strstringsteram.str(std::string());
    strstringsteram << strstring;
  }

  string StringSubst(string &strstring, const string &strfind0, const string &strfind1, const string &strfind2, const string &strfind3, const string &strreplace) {
    StringSubst(strstring,strfind0,strreplace);
    StringSubst(strstring,strfind1,strreplace);
    StringSubst(strstring,strfind2,strreplace);
    StringSubst(strstring,strfind3,strreplace);
    return strstring;
  }

  
  // ***************************************************************************
  // Function SubStrings
  // ***************************************************************************
  // Finds strings here and there.
  // Stefano Curtarolo

  int GetNLinesString(const string& str) {
    // SLOW
    //     string _str(str);
    //     int N=1;
    //     while(_str.find("\n")!=string::npos) {
    //       N++;_str=_str.substr(_str.find("\n")+1);
    //     }
    //     return N;
    stringstream strstream(str);
    return GetNLinesString(strstream);
  }
  int GetNLinesString(const stringstream& strstream) {
    // VERY SLOW    return aurostd::GetNLinesString(strstream.str());
    // FAST
    stringstream strstream_new(strstream.str());  //copy stringstream
    string line;
    int count_line = 0;
    while(getline(strstream_new, line)) {count_line++;}
    return count_line;
  }
 
  int GetNLinesFile(const string& file_name) {
    stringstream streamFILE;
    aurostd::file2stringstream(file_name, streamFILE);
    return GetNLinesString(streamFILE);
  }

  string GetLineString(const string& strstream, const int& line) {
    string _strstream(strstream),_strline;
    //  if(line>aurostd::GetNLinesString(_strstream)) return (string) "";   // TOO SLOW IF THE STRING IS LONG !
    for(int i=0;i<line;i++) {
      _strline=_strstream.substr(0,_strstream.find("\n"));
      _strstream=_strstream.substr(_strstream.find("\n")+1);
    }
    return _strline;
  }
  string GetLineString(const stringstream& strstream, const int& line) {
    return aurostd::GetLineString(strstream.str(),line);
  }

  // ***************************************************************************
  // Function SubStringsPresent
  // ***************************************************************************
  bool substring2bool(const string& strstream, const string& strsub1, bool CLEAN) {
    bool LDEBUG=FALSE;//TRUE;
    if(LDEBUG) cerr << "bool substring2bool(const string& strstream, const string& strsub1, bool CLEAN)" << endl;
    if(LDEBUG) cerr << "DEBUG substring2bool: (BEGIN) [" << strsub1 << "] " << CLEAN << endl;
    string _strstream(strstream),_strline,_strsub1(strsub1);
    string strout="";
    string::size_type idxS1;
    if(CLEAN==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    if(LDEBUG) cerr << "DEBUG substring2bool: 1 [" << strstream << "], [" << _strsub1 << "]" << endl; 
    if(_strstream.find(_strsub1)==string::npos) return (bool) FALSE;
    if(LDEBUG) cerr << "DEBUG substring2bool: 2 [" << strstream << "], [" << _strsub1 << "]" << endl; 
    vector<string> tokens;
    aurostd::string2tokens(_strstream,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      _strline=tokens.at(i);
      if(_strline.find("#")!=string::npos  
	 && _strline.find("#1")==string::npos && _strline.find("#2")==string::npos && _strline.find("#3")==string::npos // avoid space groups
	 && _strline.find("#4")==string::npos && _strline.find("#5")==string::npos && _strline.find("#6")==string::npos // avoid space groups
	 && _strline.find("#7")==string::npos && _strline.find("#8")==string::npos && _strline.find("#9")==string::npos // avoid space groups
	 ) _strline=_strline.substr(0,_strline.find("#")); // avoid space groups
      if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));
      if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));
      idxS1=_strline.find(_strsub1);
      if(idxS1!=string::npos) {
	strout=_strline.substr(_strline.find(_strsub1)+_strsub1.length());
	strout=aurostd::RemoveWhiteSpacesFromTheBack(strout);
	if(LDEBUG) cerr << "DEBUG substring2bool: (END) " << strsub1 << " " << CLEAN << endl;
 	return (bool) TRUE;
      }
    }
    return (bool) FALSE;
  }
  
  bool substring2bool(const vector<string>& vstrstream, const string& strsub1, bool CLEAN) {
    //  if(LDEBUG) cerr << "bool substring2bool(const vector<string>& vstrstream, const string& strsub1, bool CLEAN)" << endl;
    for(uint i=0;i<vstrstream.size();i++)
      if(aurostd::substring2bool(vstrstream.at(i),strsub1,CLEAN)) return TRUE;
    return FALSE;
  }
  bool substring2bool(const deque<string>& vstrstream, const string& strsub1, bool CLEAN) {
    //  if(LDEBUG) cerr << "bool substring2bool(const deque<string>& vstrstream, const string& strsub1, bool CLEAN)" << endl;
    for(uint i=0;i<vstrstream.size();i++)
      if(aurostd::substring2bool(vstrstream.at(i),strsub1,CLEAN)) return TRUE;
    return FALSE;
  }

  bool substring2bool(const stringstream& strstream, const string& strsub1, bool CLEAN) {
    //  if(LDEBUG) cerr << "bool substring2bool(const stringstream& strstream, const string& strsub1, bool CLEAN)" << endl;
    return aurostd::substring2bool(strstream.str(),strsub1,CLEAN);
  }

  bool substring2bool(const string& strstream, const string& strsub1) {
    //  if(LDEBUG)   cerr << "bool substring2bool(const string& strstream, const string& strsub1)" << endl;
    return (bool) aurostd::substring2bool(strstream,strsub1,FALSE);
  }

  bool substring2bool(const vector<string>& vstrstream, const string& strsub1) {
    //  if(LDEBUG) cerr << "bool substring2bool(const vector<string>& vstrstream, const string& strsub1)" << endl;
    for(uint i=0;i<vstrstream.size();i++)
      if(aurostd::substring2bool(vstrstream.at(i),strsub1)) return TRUE;
    return FALSE;
  }
  bool substring2bool(const deque<string>& vstrstream, const string& strsub1) {
    //  if(LDEBUG) cerr << "bool substring2bool(const deque<string>& vstrstream, const string& strsub1)" << endl;
    for(uint i=0;i<vstrstream.size();i++)
      if(aurostd::substring2bool(vstrstream.at(i),strsub1)) return TRUE;
    return FALSE;
  }

  bool substring2bool(const stringstream& strstream, const string& strsub1) {
    //  if(LDEBUG) cerr << "bool substring2bool(const stringstream& strstream, const string& strsub1)" << endl;
    return (bool) aurostd::substring2bool(strstream.str(),strsub1,FALSE);
  }
  bool substring_present_file(const string& FileName, const string& strsub1) {
    //  if(LDEBUG) cerr << "bool substring_present_file(const string& FileName, const string& strsub1)" << endl;
    return (bool) aurostd::substring2bool(FileName,strsub1,FALSE);
  }
  bool substring_present_file_FAST(const string& FileName, const string& strsub1) {
    //  if(LDEBUG) cerr << "bool substring_present_file_FAST(const string& FileName, const string& strsub1)" << endl;
    // be careful, this does not filter-out # comments
    return (bool) substring_present_file_FAST(FileName,strsub1,FALSE);
  }

  bool withinList(const vector<string>& list,const string& str) {
    for(uint i=0;i<list.size();i++){if(list[i]==str){return true;}}
    return false;
  }

  // ***************************************************************************
  bool substring_present_file(const string& FileName, const string& strsub1, bool CLEAN) {
    string StringFile;
    ifstream FileFile;
    FileFile.open(FileName.c_str(),std::ios::in);
    FileFile.clear();FileFile.seekg(0);
    if(!FileFile) {
      cerr << "ERROR  FileName=" << FileName << "   not found" << endl;
      return FALSE;
      //    exit(0);
    }
    StringFile="";char c; while (FileFile.get(c)) StringFile+=c;
    FileFile.close();
    return aurostd::substring2bool(StringFile,strsub1,CLEAN);
  }

  bool substring_present_file_FAST(const string& FileName, const string& strsub1, bool CLEAN) {
    // be careful, this does not filter-out # comments
    string temp_file=aurostd::TmpFileCreate("substring");
    ostringstream aus;
    ifstream FileFile;
    int found=0;
    aurostd::StringstreamClean(aus);
    aus << "rm -f " << temp_file  << endl;
    if(CLEAN==FALSE) aus << "cat " << FileName << " | grep -c \"" << strsub1 << "\" > " << temp_file  << endl;
    if(CLEAN==TRUE ) aus << "cat " << FileName << " | sed \"s/ //g\" | sed \"s/\t//g\"  | grep -c \"" << aurostd::RemoveWhiteSpaces(strsub1) << "\" > " << temp_file  << endl;
    aus << "echo >> " << temp_file << endl; // to give EOL
    aurostd::execute(aus);
    if(!aurostd::FileExist(FileName)) {
      cerr << "ERROR: substring_present_file_FAST: file input not found =" << FileName << endl;
      exit(0);
    }
    if(!aurostd::FileExist(temp_file)) {
      cerr << "ERROR: substring_present_file_FAST: file output not found =" << temp_file << endl;
      exit(0);
    }
    FileFile.open(temp_file.c_str(),std::ios::in);
    FileFile.clear();FileFile.seekg(0);
    FileFile >> found;
    FileFile.close();
    aurostd::StringstreamClean(aus);
#ifndef _AFLOW_TEMP_PRESERVE_
    aurostd::RemoveFile(temp_file);
#endif
    if(found>0) return TRUE;
    return FALSE;
  }

  // ***************************************************************************
  // Function SubStringsPresent and EXTRACT
  // ***************************************************************************
  string substring2string(const string& strstream, const string& strsub1, bool CLEAN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "DEBUG substring2string3: (BEGIN) " << strsub1 << " " << CLEAN << endl;
    string _strstream(strstream),_strline,_strsub1(strsub1);
    string strout="";
    string::size_type idxS1;
    if(CLEAN==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    if(_strstream.find(_strsub1)==string::npos) return (string) strout;
    //  transform(_strstream.begin(),_strstream.end(),_strstream.begin(),toupper); // pout everything UPPER
    vector<string> tokens;
    aurostd::string2tokens(_strstream,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      _strline=tokens.at(i);
      if(_strline.find(COMMENT_NEGLECT_1)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_1));
      if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));
      if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));
      idxS1=_strline.find(_strsub1);
      if(idxS1!=string::npos) {
	strout=_strline.substr(_strline.find(_strsub1)+_strsub1.length());
	strout=aurostd::RemoveWhiteSpacesFromTheBack(strout);
	if(LDEBUG) cerr << "DEBUG substring2string3: (END) " << strsub1 << " " << CLEAN << endl;
	return (string) strout;
      }
    }
    return (string) strout;
  }

  string substring2string(const string& strstream, const string& strsub1, const string& strsub2, bool CLEAN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "DEBUG substring2string5: (BEGIN) " << strsub1 << " " << CLEAN << endl;
    string _strstream(strstream),_strline,_strsub1(strsub1),_strsub2(strsub2);
    string strout="";
    string::size_type idxS1,idxS2;
    if(CLEAN==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    if(_strstream.find(_strsub1)==string::npos) return (string) strout;
    if(_strstream.find(_strsub2)==string::npos) return (string) strout;
    //  transform(_strstream.begin(),_strstream.end(),_strstream.begin(),toupper); // pout everything UPPER
    vector<string> tokens;
    aurostd::string2tokens(_strstream,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      _strline=tokens.at(i);
      if(_strline.find(COMMENT_NEGLECT_1)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_1));
      if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));
      if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));
      idxS1=_strline.find(_strsub1);
      idxS2=_strline.find(_strsub2);
      if(idxS1!=string::npos && idxS2!=string::npos && idxS1<=idxS2) {
	strout=_strline.substr(std::max(_strline.find(_strsub2)+_strsub2.length(),_strline.find(_strsub1)+_strsub1.length()));
	strout=aurostd::RemoveWhiteSpacesFromTheBack(strout);
	if(LDEBUG) cerr << "DEBUG substring2string5: (END) " << strsub1 << " " << CLEAN << endl;
	return (string) strout;
      }
    }
    return (string) strout;
  }

  string substring2string(const string& strstream, const string& strsub1) {
    return (string) substring2string(strstream,strsub1,FALSE);
  }
  string substring2string(const string& strstream, const string& strsub1, const string& strsub2) {
    return (string) substring2string(strstream,strsub1,strsub2,FALSE);
  }

  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1, bool CLEAN) {
    return string2utype<utype>(substring2string(strstream,strsub1,CLEAN));
  }
  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1) {
    return string2utype<utype>(substring2string(strstream,strsub1));
  }
  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1, const string& strsub2, bool CLEAN) {
    return string2utype<utype>(substring2string(strstream,strsub1,strsub2,CLEAN));
  }
  template<typename utype> utype substring2utype(const string& strstream, const string& strsub1, const string& strsub2) {
    return string2utype<utype>(substring2string(strstream,strsub1,strsub2));
  }
      
  // ***************************************************************************
  // Function SubStringsPresentExtractString and other
  // ***************************************************************************
  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, bool CLEAN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "DEBUG substring2strings3: (BEGIN) " << strsub1 << " " << CLEAN << endl;
    string _strstream(strstream),_strline,_strsub1(strsub1);
    string::size_type idxS1;
    vstringout.clear(); // clear so it is empty
    if(CLEAN==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    if(_strstream.find(_strsub1)==string::npos) return 0; // there is not
    //  transform(_strstream.begin(),_strstream.end(),_strstream.begin(),toupper); // pout everything UPPER
    vector<string> tokens;
    aurostd::string2tokens(_strstream,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      _strline=tokens.at(i);
      if(_strline.find(COMMENT_NEGLECT_1)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_1));
      if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));
      if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));
      idxS1=_strline.find(_strsub1);
      if(idxS1!=string::npos) {
	vstringout.push_back(aurostd::RemoveWhiteSpacesFromTheBack(_strline.substr(_strline.find(_strsub1)+_strsub1.length())));
      }
    }
    if(LDEBUG) cerr << "DEBUG substring2strings3: (END) " << strsub1 << " " << vstringout.size() << " " << CLEAN << endl;
    return vstringout.size();
  }

  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2, bool CLEAN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "DEBUG substring2strings5: (BEGIN) " << strsub1 << " " << strsub2 << " " << CLEAN << endl;
    string _strstream(strstream),_strline,_strsub1(strsub1),_strsub2(strsub2);
    string::size_type idxS1,idxS2;
    vstringout.clear(); // clear so it is empty
    if(CLEAN==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    if(_strstream.find(_strsub1)==string::npos) return 0; // there is not
    if(_strstream.find(_strsub2)==string::npos) return 0; // there is not
    //  transform(_strstream.begin(),_strstream.end(),_strstream.begin(),toupper); // pout everything UPPER
    vector<string> tokens;
    aurostd::string2tokens(_strstream,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      _strline=tokens.at(i);
      if(_strline.find(COMMENT_NEGLECT_1)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_1));
      if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));
      if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));
      idxS1=_strline.find(_strsub1);
      idxS2=_strline.find(_strsub2);
      if(idxS1!=string::npos && idxS2!=string::npos) {
	vstringout.push_back(aurostd::RemoveWhiteSpacesFromTheBack(_strline.substr(std::max(_strline.find(_strsub2)+_strsub2.length(),_strline.find(_strsub1)+_strsub1.length()))));
      }
    }
    if(LDEBUG) cerr << "DEBUG substring2string5: (END) " << strsub1 << " " << strsub2 << " " << vstringout.size() << " " << CLEAN << endl;
    return vstringout.size();
  }
  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1) {
    return substring2strings(strstream,vstringout,strsub1,FALSE);
  }
  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2) {
    return substring2strings(strstream,vstringout,strsub1,strsub2,FALSE);
  }

  template<typename utype> uint substring2utypes(const string& strstream, vector<int> &vintout, const string& strsub1, bool CLEAN) {
    string _strstream(strstream),_strsub1(strsub1);
    vintout.clear();
    vector<string> vstringout;
    uint i=aurostd::substring2strings(_strstream,vstringout,_strsub1,CLEAN);
    for(i=0;i<vstringout.size();i++) vintout.push_back(string2utype<utype>(vstringout.at(i)));
    return vintout.size();
  }
  template<typename utype> uint substring2utypes(const string& strstream, vector<int> &vintout, const string& strsub1) {
    return substring2utypes<utype>(strstream,vintout,strsub1,FALSE);
  }
  template<typename utype> uint substring2utypes(const string& strstream, vector<int> &vintout, const string& strsub1, const string& strsub2, bool CLEAN) {
    string _strstream(strstream),_strsub1(strsub1),_strsub2(strsub2);
    vintout.clear();
    vector<string> vstringout;
    uint i=aurostd::substring2strings(_strstream,vstringout,_strsub1,_strsub2,CLEAN);
    for(i=0;i<vstringout.size();i++) vintout.push_back(string2utype<utype>(vstringout.at(i)));
    return vintout.size();
  }
  template<typename utype> uint substring2utypes(const string& strstream, vector<int> &vintout, const string& strsub1, const string& strsub2) {
    return substring2utypes<utype>(strstream,vintout,strsub1,strsub2,FALSE);
  }
}

// ***************************************************************************
// FUNCTION HTML LATEX TXT

namespace aurostd {

  // http://www.w3schools.com/tags/ref_entities.asp
  // http://en.wikibooks.org/wiki/LaTeX/Accents

  string html2latex(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"_","\\_");
    aurostd::StringSubst(out,"<sub>","$_{");aurostd::StringSubst(out,"</sub>","}$");
    aurostd::StringSubst(out,"<i>","{\\it ");aurostd::StringSubst(out,"</i>","}");
    aurostd::StringSubst(out,"<b>","{\\bf "); aurostd::StringSubst(out,"</b>","}");
    aurostd::StringSubst(out,"<blink>","{\\bf "); aurostd::StringSubst(out,"</blink>","}");
    aurostd::StringSubst(out,"MgB2","MgB$_2$");
    aurostd::StringSubst(out,"Schuttler","Sch\\\"uttler");
    aurostd::StringSubst(out,"Csányi","Cs\\'anyi");aurostd::StringSubst(out,"Csanyi","Cs\\'anyi");
    aurostd::StringSubst(out,"Pólya","P\\'{o}lya");
    if(!aurostd::substring2bool(out,"Rosenbrock")) aurostd::StringSubst(out,"Rosen","Ros\\'en");
    // string bar="";//;bar.at(0)=92;

    // http://en.wikibooks.org/wiki/LaTeX/Accents
    // umlaut
    aurostd::StringSubst(out,"&auml;","\\\"{a}");aurostd::StringSubst(out,"&Auml;","\\\"{A}");
    aurostd::StringSubst(out,"&euml;","\\\"{e}");aurostd::StringSubst(out,"&Euml;","\\\"{E}");
    aurostd::StringSubst(out,"&iuml;","\\\"{i}");aurostd::StringSubst(out,"&Iuml;","\\\"{I}");
    aurostd::StringSubst(out,"&ouml;","\\\"{o}");aurostd::StringSubst(out,"&Ouml;","\\\"{O}");
    aurostd::StringSubst(out,"&uuml;","\\\"{u}");aurostd::StringSubst(out,"&Uuml;","\\\"{U}");
    // grave accent
    aurostd::StringSubst(out,"&agrave;","\\`{a}");aurostd::StringSubst(out,"&Agrave;","\\`{A}");
    aurostd::StringSubst(out,"&egrave;","\\`{e}");aurostd::StringSubst(out,"&Egrave;","\\`{E}");
    aurostd::StringSubst(out,"&igrave;","\\`{i}");aurostd::StringSubst(out,"&Igrave;","\\`{I}");
    aurostd::StringSubst(out,"&ograve;","\\`{o}");aurostd::StringSubst(out,"&Ograve;","\\`{O}");
    aurostd::StringSubst(out,"&ugrave;","\\`{u}");aurostd::StringSubst(out,"&Ugrave;","\\`{U}");
    // acute accent
    aurostd::StringSubst(out,"&aacute;","\\'{a}");aurostd::StringSubst(out,"&Aacute;","\\'{A}");
    aurostd::StringSubst(out,"&eacute;","\\'{e}");aurostd::StringSubst(out,"&Eacute;","\\'{E}");
    aurostd::StringSubst(out,"&iacute;","\\'{i}");aurostd::StringSubst(out,"&Iacute;","\\'{I}");
    aurostd::StringSubst(out,"&oacute;","\\'{o}");aurostd::StringSubst(out,"&Oacute;","\\'{O}");
    aurostd::StringSubst(out,"&uacute;","\\'{u}");aurostd::StringSubst(out,"&Uacute;","\\'{U}");
    // tilde
    aurostd::StringSubst(out,"&atilde;","\\~{a}");aurostd::StringSubst(out,"&Atilde;","\\~{A}");
    aurostd::StringSubst(out,"&etilde;","\\~{e}");aurostd::StringSubst(out,"&Etilde;","\\~{E}");
    aurostd::StringSubst(out,"&itilde;","\\~{i}");aurostd::StringSubst(out,"&Itilde;","\\~{I}");
    aurostd::StringSubst(out,"&otilde;","\\~{o}");aurostd::StringSubst(out,"&Otilde;","\\~{O}");
    aurostd::StringSubst(out,"&utilde;","\\~{u}");aurostd::StringSubst(out,"&Utilde;","\\~{U}");
    // circ
    aurostd::StringSubst(out,"&acirc;","\\^{a}");aurostd::StringSubst(out,"&Acirc;","\\^{A}");
    aurostd::StringSubst(out,"&ecirc;","\\^{e}");aurostd::StringSubst(out,"&Ecirc;","\\^{E}");
    aurostd::StringSubst(out,"&icirc;","\\^{i}");aurostd::StringSubst(out,"&Icirc;","\\^{I}");
    aurostd::StringSubst(out,"&ocirc;","\\^{o}");aurostd::StringSubst(out,"&Ocirc;","\\^{O}");
    aurostd::StringSubst(out,"&ucirc;","\\^{u}");aurostd::StringSubst(out,"&Ucirc;","\\^{U}");
    // ring
    aurostd::StringSubst(out,"&aring;","\\r{a}");aurostd::StringSubst(out,"&Aring;","\\r{A}");
    aurostd::StringSubst(out,"&ering;","\\r{e}");aurostd::StringSubst(out,"&Ering;","\\r{E}");
    aurostd::StringSubst(out,"&iring;","\\r{i}");aurostd::StringSubst(out,"&Iring;","\\r{I}");
    aurostd::StringSubst(out,"&oring;","\\r{o}");aurostd::StringSubst(out,"&Oring;","\\r{O}");
    aurostd::StringSubst(out,"&uring;","\\r{u}");aurostd::StringSubst(out,"&Uring;","\\r{U}");
    // cedil
    aurostd::StringSubst(out,"&acedil;","\\c{a}");aurostd::StringSubst(out,"&Acedil;","\\c{A}");
    aurostd::StringSubst(out,"&ecedil;","\\c{e}");aurostd::StringSubst(out,"&Ecedil;","\\c{E}");
    aurostd::StringSubst(out,"&icedil;","\\c{i}");aurostd::StringSubst(out,"&Icedil;","\\c{I}");
    aurostd::StringSubst(out,"&ocedil;","\\c{o}");aurostd::StringSubst(out,"&Ocedil;","\\c{O}");
    aurostd::StringSubst(out,"&ucedil;","\\c{u}");aurostd::StringSubst(out,"&Ucedil;","\\c{U}");
    // math
    aurostd::StringSubst(out,"&Alpha;","$\\Alpha$");aurostd::StringSubst(out,"&alpha;","$\\alpha$");
    aurostd::StringSubst(out,"&Beta;","$\\Βeta$");aurostd::StringSubst(out,"&beta;","$\\beta$");
    aurostd::StringSubst(out,"&Gamma;","$\\Gamma$");aurostd::StringSubst(out,"&gamma;","$\\gamma$");
    aurostd::StringSubst(out,"&Delta;","$\\Delta$");aurostd::StringSubst(out,"&delta;","$\\delta$");
    aurostd::StringSubst(out,"&Epsilon;","$\\Εpsilon$");aurostd::StringSubst(out,"&epsilon;","$\\epsilon$");
    aurostd::StringSubst(out,"&Zeta;","$\\Ζeta$");aurostd::StringSubst(out,"&zeta;","$\\zeta$");
    aurostd::StringSubst(out,"&Eta;","$\\Eta$");aurostd::StringSubst(out,"&eta;","$\\eta$");
    aurostd::StringSubst(out,"&Theta;","$\\Theta$");aurostd::StringSubst(out,"&theta;","$\\theta$");
    aurostd::StringSubst(out,"&Iota;","$\\Ιiota$");aurostd::StringSubst(out,"&iota;","$\\iota$");
    aurostd::StringSubst(out,"&Kappa;","$\\Kappa$");aurostd::StringSubst(out,"&kappa;","$\\kappa$");
    aurostd::StringSubst(out,"&Lambda;","$\\Lambda$");aurostd::StringSubst(out,"&lambda;","$\\lambda$");
    aurostd::StringSubst(out,"&Mu;","$\\Mu$");aurostd::StringSubst(out,"&mu;","$\\mu$");
    aurostd::StringSubst(out,"&Nu;","$\\Νu$");aurostd::StringSubst(out,"&nu;","$\\nu$");
    aurostd::StringSubst(out,"&Xi;","$\\Xi$");aurostd::StringSubst(out,"&xi;","$\\xi$");
    aurostd::StringSubst(out,"&Omicron;","$\\Omicron$");aurostd::StringSubst(out,"&omicron;","$\\omicron$");
    aurostd::StringSubst(out,"&Pi;","$\\Pi$");aurostd::StringSubst(out,"&pi;","$\\pi$");
    aurostd::StringSubst(out,"&Rho;","$\\Rho$");aurostd::StringSubst(out,"&rho;","$\\rho$");
    aurostd::StringSubst(out,"&Sigma;","$\\Sigma$");aurostd::StringSubst(out,"&sigma;","$\\sigma$");
    aurostd::StringSubst(out,"&Tau;","$\\Tau$");aurostd::StringSubst(out,"&tau;","$\\tau$");
    aurostd::StringSubst(out,"&Upsilon;","$\\Upsilon$");aurostd::StringSubst(out,"&upsilon;","$\\upsilon$");
    aurostd::StringSubst(out,"&Phi;","$\\Phi$");aurostd::StringSubst(out,"&phi;","$\\phi$");
    aurostd::StringSubst(out,"&Chi;","$\\Chi$");aurostd::StringSubst(out,"&chi;","$\\chi$");
    aurostd::StringSubst(out,"&Psi;","$\\Psi$");aurostd::StringSubst(out,"&psi;","$\\psi$");
    aurostd::StringSubst(out,"&Omega;","$\\Omega$");aurostd::StringSubst(out,"&omega;","$\\omega$");
    aurostd::StringSubst(out,"&thetasym","$\\thetasym$");
    // FINAL
    aurostd::StringSubst(out,"&","\\&");

    return out;
  }

  string html2txt(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"<sub>","");aurostd::StringSubst(out,"</sub>","");
    aurostd::StringSubst(out,"<i>","");aurostd::StringSubst(out,"</i>","");
    aurostd::StringSubst(out,"<b>",""); aurostd::StringSubst(out,"</b>","");
    aurostd::StringSubst(out,"MgB2","MgB2");
    aurostd::StringSubst(out,"&","&");
    aurostd::StringSubst(out,"_","");aurostd::StringSubst(out,"\\","");
    return out;
  }


  // ***************************************************************************
  // Function aurostd::string2latex
  // ***************************************************************************
  string string2latex(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"_pv","_{pv}");aurostd::StringSubst(out,"_sv","_{sv}");aurostd::StringSubst(out,"_h","_{h}");
    aurostd::StringSubst(out,"_d","_{d}");aurostd::StringSubst(out,"_s","_{s}");
    aurostd::StringSubst(out,"_1","_{1}");aurostd::StringSubst(out,"_2","_{2}");aurostd::StringSubst(out,"_3","_{3}");
    return out;
  }
  
  // ***************************************************************************
  // Function aurostd::latex2html
  // ***************************************************************************
  string latex2html(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"\\alpha","&alpha;");aurostd::StringSubst(out,"\\Alpha","&Alpha;");
    aurostd::StringSubst(out,"\\beta","&beta;");aurostd::StringSubst(out,"\\Beta","&Beta;");
    aurostd::StringSubst(out,"\\epsilon","&epsilon;");aurostd::StringSubst(out,"\\Epsilon","&Epsilon;");
    aurostd::StringSubst(out,"\\eta","&eta;");aurostd::StringSubst(out,"\\Eta","&Eta;");
    aurostd::StringSubst(out,"\\gamma","&gamma;");aurostd::StringSubst(out,"\\Gamma","&Gamma;");
    aurostd::StringSubst(out,"\\delta","&delta;");aurostd::StringSubst(out,"\\Delta","&Delta;");
    aurostd::StringSubst(out,"\\omega","&omega;");aurostd::StringSubst(out,"\\Omega","&Omega;");
    aurostd::StringSubst(out,"\\sigma","&sigma;");aurostd::StringSubst(out,"\\Sigma","&Sigma;");
    aurostd::StringSubst(out,"_{a}","<sub>a</sub>");aurostd::StringSubst(out,"_a","<sub>a</sub>"); 
    aurostd::StringSubst(out,"_{b}","<sub>b</sub>");aurostd::StringSubst(out,"_b","<sub>b</sub>");
    aurostd::StringSubst(out,"_{c}","<sub>d</sub>");aurostd::StringSubst(out,"_c","<sub>d</sub>");
    aurostd::StringSubst(out,"_{d}","<sub>d</sub>");aurostd::StringSubst(out,"_d","<sub>d</sub>");
    aurostd::StringSubst(out,"_{h}","<sub>h</sub>");aurostd::StringSubst(out,"_h","<sub>h</sub>");
    aurostd::StringSubst(out,"_{s}","<sub>s</sub>");aurostd::StringSubst(out,"_s","<sub>s</sub>");
    aurostd::StringSubst(out,"_{v}","<sub>v</sub>");aurostd::StringSubst(out,"_v","<sub>v</sub>");
    aurostd::StringSubst(out,"_{AB}","<sub>AB</sub>");
    aurostd::StringSubst(out,"_{AB2}","<sub>AB2</sub>");
    aurostd::StringSubst(out,"_{A2B2}","<sub>A2B2</sub>");
    aurostd::StringSubst(out,"_{AB3}","<sub>AB3</sub>"); 
    for(uint i=0;i<100;i++) aurostd::StringSubst(out,"_{"+aurostd::utype2string(i)+"}","<sub>"+aurostd::utype2string(i)+"</sub>");
    for(uint i=0;i<10;i++) aurostd::StringSubst(out,"_"+aurostd::utype2string(i)+"","<sub>"+aurostd::utype2string(i)+"</sub>"); // patch
    for(uint i1=0;i1<=3;i1++)
      for(uint i2=0;i2<=3;i2++)
	for(uint i3=0;i3<=3;i3++)
	  aurostd::StringSubst(out,
			       "^{["+aurostd::utype2string(i1)+aurostd::utype2string(i2)+aurostd::utype2string(i3)+"]}",
			       "<sup>"+aurostd::utype2string(i1)+aurostd::utype2string(i2)+aurostd::utype2string(i3)+"</sup>");
    string s="AB";
    stringstream ss;
    for(uint i1=0;i1<=1;i1++)
      for(uint i2=0;i2<=1;i2++)
	for(uint i3=0;i3<=1;i3++)
	  for(uint i4=0;i4<=1;i4++)
	    for(uint i5=0;i5<=1;i5++) {
	      ss.clear();ss.str("");
	      ss << s.at(i1) << s.at(i2) << s.at(i3) << s.at(i4) << s.at(i5); 
	      aurostd::StringSubst(out,"_{"+ss.str()+"}","<sub>"+ss.str()+"</sub>");
	    }
    //    return out;
    //  }  
    //  string latex2html(const string& str) {
    // string out=str;
    aurostd::StringSubst(out,"\\&","&");
    aurostd::StringSubst(out,"MgB$_2$","MgB<sub>2</sub>");
    //  aurostd::StringSubst(out,"<sub>","$_{");aurostd::StringSubst(out,"</sub>","}$");
    //  aurostd::StringSubst(out,"<i>","{\\it ");aurostd::StringSubst(out,"</i>","}");
    // aurostd::StringSubst(out,"<b>","{\\bf "); aurostd::StringSubst(out,"</b>","}");
    // aurostd::StringSubst(out,"&","\\&");
    //  aurostd::StringSubst(out,"Schuttler","Sch\\\"uttler");
    //  aurostd::StringSubst(out,"Csányi","Cs\\'anyi");aurostd::StringSubst(out,"Csanyi","Cs\\'anyi");
    // aurostd::StringSubst(out,"Rosen","Ros\\'en");
    // umlaut
    aurostd::StringSubst(out,"\\:a","&auml;");aurostd::StringSubst(out,"\\:A","&Auml;");
    aurostd::StringSubst(out,"\\:e","&euml;");aurostd::StringSubst(out,"\\:E","&Euml;");
    aurostd::StringSubst(out,"\\:i","&iuml;");aurostd::StringSubst(out,"\\:I","&Iuml;");
    aurostd::StringSubst(out,"\\:o","&ouml;");aurostd::StringSubst(out,"\\:O","&Ouml;");
    aurostd::StringSubst(out,"\\:u","&uuml;");aurostd::StringSubst(out,"\\:U","&Uuml;");
    // grave accent
    aurostd::StringSubst(out,"\\`a","&agrave;");aurostd::StringSubst(out,"\\`A","&Agrave;");
    aurostd::StringSubst(out,"\\`e","&egrave;");aurostd::StringSubst(out,"\\`E","&Egrave;");
    aurostd::StringSubst(out,"\\`i","&igrave;");aurostd::StringSubst(out,"\\`I","&Igrave;");
    aurostd::StringSubst(out,"\\`o","&ograve;");aurostd::StringSubst(out,"\\`O","&Ograve;");
    aurostd::StringSubst(out,"\\`u","&ugrave;");aurostd::StringSubst(out,"\\`U","&Ugrave;");
    // acute accent
    aurostd::StringSubst(out,"\\'a","&aacute;");aurostd::StringSubst(out,"\\'A","&Aacute;");
    aurostd::StringSubst(out,"\\'e","&eacute;");aurostd::StringSubst(out,"\\'E","&Eacute;");
    aurostd::StringSubst(out,"\\'i","&iacute;");aurostd::StringSubst(out,"\\'I","&Iacute;");
    aurostd::StringSubst(out,"\\'o","&oacute;");aurostd::StringSubst(out,"\\'O","&Oacute;");
    aurostd::StringSubst(out,"\\'u","&uacute;");aurostd::StringSubst(out,"\\'U","&Uacute;");
    // tilde
    aurostd::StringSubst(out,"\\~a","&atilde;");aurostd::StringSubst(out,"\\~A","&Atilde;");
    aurostd::StringSubst(out,"\\~e","&etilde;");aurostd::StringSubst(out,"\\~E","&Etilde;");
    aurostd::StringSubst(out,"\\~i","&itilde;");aurostd::StringSubst(out,"\\~I","&Itilde;");
    aurostd::StringSubst(out,"\\~o","&otilde;");aurostd::StringSubst(out,"\\~O","&Otilde;");
    aurostd::StringSubst(out,"\\~u","&utilde;");aurostd::StringSubst(out,"\\~U","&Utilde;");
    return out;
  }

  string latex2txt(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"\\&","&");
    aurostd::StringSubst(out,"MgB$_2$","MgB2");
    aurostd::StringSubst(out,"<sub>","");aurostd::StringSubst(out,"</sub>","");
    aurostd::StringSubst(out,"<i>","");aurostd::StringSubst(out,"</i>","");
    aurostd::StringSubst(out,"<b>","");aurostd::StringSubst(out,"</b>","");
    return out;
  }
}

// ***************************************************************************
// SORT WORLD
// ----------------------------------------------------------------------------
// sort for vector (starting from xvector)

namespace aurostd {
  template<class utype1> // function quicksort
  void sort(vector<utype1>& arr) {
    xvector<utype1> xarr(arr.size());
    for(uint i=0;i<arr.size();i++) xarr[i+1]=arr.at(i);
    //  aurostd::sort(xarr.rows,xarr);
    aurostd::sort(xarr);
    arr.clear();
    for(int i=0;i<xarr.rows;i++) {
      arr.push_back(xarr[i+1]);
    }
  }

  template<class utype1,class utype2> // function quicksort
  void sort(vector<utype1>& arr, vector<utype2>& brr) {
    xvector<utype1> xarr(arr.size());
    xvector<utype2> xbrr(brr.size());
    for(uint i=0;i<arr.size();i++) xarr[i+1]=arr.at(i);
    for(uint i=0;i<brr.size();i++) xbrr[i+1]=brr.at(i);
    aurostd::sort2(xarr.rows,xarr,xbrr);
    // aurostd::sort2(xarr,xbrr);
    arr.clear();brr.clear();
    for(int i=0;i<xarr.rows;i++) {
      arr.push_back(xarr[i+1]);
      brr.push_back(xbrr[i+1]);
    }
  }

  template<class utype1,class utype2,class utype3> // function quicksort
  void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr) {
    xvector<utype1> xarr(arr.size());
    xvector<utype2> xbrr(brr.size());
    xvector<utype3> xcrr(crr.size());
    for(uint i=0;i<arr.size();i++) xarr[i+1]=arr.at(i);
    for(uint i=0;i<brr.size();i++) xbrr[i+1]=brr.at(i);
    for(uint i=0;i<crr.size();i++) xcrr[i+1]=crr.at(i);
    aurostd::sort3(xarr.rows,xarr,xbrr,xcrr);
    // aurostd::sort3(xarr,xbrr,xcrr);
    arr.clear();brr.clear();crr.clear();
    for(int i=0;i<xarr.rows;i++) {
      arr.push_back(xarr[i+1]);
      brr.push_back(xbrr[i+1]);
      crr.push_back(xcrr[i+1]);
    }
  }

  template<class utype1,class utype2,class utype3,class utype4> // function quicksort
  void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr, vector<utype4>& drr) {
    xvector<utype1> xarr(arr.size());
    xvector<utype2> xbrr(brr.size());
    xvector<utype3> xcrr(crr.size());
    xvector<utype4> xdrr(drr.size());
    for(uint i=0;i<arr.size();i++) xarr[i+1]=arr.at(i);
    for(uint i=0;i<brr.size();i++) xbrr[i+1]=brr.at(i);
    for(uint i=0;i<crr.size();i++) xcrr[i+1]=crr.at(i);
    for(uint i=0;i<drr.size();i++) xdrr[i+1]=drr.at(i);
    //    aurostd::sort4(xarr.rows,xarr,xbrr,xcrr,xdrr);
    aurostd::sort4(xarr,xbrr,xcrr,xdrr);
    arr.clear();brr.clear();crr.clear();drr.clear();
    for(int i=0;i<xarr.rows;i++) {
      arr.push_back(xarr[i+1]);
      brr.push_back(xbrr[i+1]);
      crr.push_back(xcrr[i+1]);
      drr.push_back(xdrr[i+1]);
    }
  }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// sort for vector of strings
namespace aurostd {
  void sort(vector<string>& arg) {
    sort(arg.begin(),arg.end(),aurostd::_sort_string_());
  }
  void sort(deque<string>& arg) {
    std::sort(arg.begin(),arg.end(),aurostd::_sort_string_());
  }
  void rsort(vector<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
  }
  void rsort(deque<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
  }
}

// sort_remove_duplicates for vector of strings
namespace aurostd {
  void sort_remove_duplicates(vector<string>& arg) {
    sort(arg.begin(),arg.end(),aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
  void sort_remove_duplicates(deque<string>& arg) {
    std::sort(arg.begin(),arg.end(),aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
  void rsort_remove_duplicates(vector<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
  void rsort_remove_duplicates(deque<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
}


// ----------------------------------------------------------------------------
// sort for vector/deque of string_int
namespace aurostd {
  void sort(vector<string>& varg1,vector<int>& varg2) {
    vector<aurostd::_string_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_int_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
  void sort(deque<string>& varg1,deque<int>& varg2) {
    deque<aurostd::_string_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_int_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_double
namespace aurostd {
  void sort(vector<string>& varg1,vector<double>& varg2) {
    vector<aurostd::_string_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_double_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
  void sort(deque<string>& varg1,deque<double>& varg2) {
    deque<aurostd::_string_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_double_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2) {
    vector<aurostd::_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2) {
    deque<aurostd::_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
}


// ----------------------------------------------------------------------------
// sort for vector/deque of double_int
// HERE THEY ARE

namespace aurostd {
  void sort(vector<double>& varg1,vector<int>& varg2) {
    vector<aurostd::_double_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_double_int_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
  void sort(deque<double>& varg1,deque<int>& varg2) {
    deque<aurostd::_double_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_double_int_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
    //   }
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of double_double
namespace aurostd {
  void sort(vector<double>& varg1,vector<double>& varg2) {
    vector<aurostd::_double_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_double_double_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
  void sort(deque<double>& varg1,deque<double>& varg2) {
    deque<aurostd::_double_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_double_double_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of double_string
namespace aurostd {
  void sort(vector<double>& varg1,vector<string>& varg2) {
    vector<aurostd::_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
  void sort(deque<double>& varg1,deque<string>& varg2) {
    deque<aurostd::_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);}
    sort(vv.begin(),vv.end(),_sort_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_int_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<int>& varg2,vector<string>& varg3) {
    vector<aurostd::_string_int_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_int_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;}
  }
  void sort(deque<string>& varg1,deque<int>& varg2,deque<string>& varg3) {
    deque<aurostd::_string_int_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_int_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_double_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<double>& varg2,vector<string>& varg3) {
    vector<aurostd::_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;}
  }
  void sort(deque<string>& varg1,deque<double>& varg2,deque<string>& varg3) {
    deque<aurostd::_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2,vector<string>& varg3) {
    vector<aurostd::_string_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2,deque<string>& varg3) {
    deque<aurostd::_string_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_double_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<string>& varg4) {
    vector<aurostd::_string_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);vv.at(i).arg4=varg4.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;varg4.at(i)=vv.at(i).arg4;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<string>& varg4) {
    deque<aurostd::_string_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);vv.at(i).arg4=varg4.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;varg4.at(i)=vv.at(i).arg4;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_double_double_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<double>& varg4,vector<string>& varg5) {
    vector<aurostd::_string_string_double_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);vv.at(i).arg4=varg4.at(i);vv.at(i).arg5=varg5.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_double_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;varg4.at(i)=vv.at(i).arg4;varg5.at(i)=vv.at(i).arg5;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<double>& varg4,deque<string>& varg5) {
    deque<aurostd::_string_string_double_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv.at(i).arg1=varg1.at(i);vv.at(i).arg2=varg2.at(i);vv.at(i).arg3=varg3.at(i);vv.at(i).arg4=varg4.at(i);vv.at(i).arg5=varg5.at(i);}
    sort(vv.begin(),vv.end(),_sort_string_string_double_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1.at(i)=vv.at(i).arg1;varg2.at(i)=vv.at(i).arg2;varg3.at(i)=vv.at(i).arg3;varg4.at(i)=vv.at(i).arg4;varg5.at(i)=vv.at(i).arg5;}
  }
}

// ***************************************************************************
// Function some statistical stuff
// combinations
// ***************************************************************************
template<class utype> utype combinations(utype n,utype k) { // http://en.wikipedia.org/wiki/Combination // C^n_k=n!/k!(n-k)!   hard to calculate
  double cnk=1.0;
  for(utype i=0;i<=k-1;i++) cnk=cnk*(n-i)/(k-i);
  return (utype) cnk;
}

template<class utype> utype Cnk(utype n,utype k) { return combinations(n,k);}  // http://en.wikipedia.org/wiki/Combination


// ***************************************************************************
// aurostd::ShiftFirstColumn(const vector<vector<double> >& a, const double& value)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > ShiftFirstColumn(const vector<vector<double> >& vva, const double& value) {
    //change value in the first column (usually menas energy in DOS)
    vector<vector<double> > vvb=vva;
    for (uint i=0; i<vvb.size(); i++) {
      vvb.at(i).at(0)=vvb.at(i).at(0) - value;
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::ShrinkValuesExceptFirstColumn(const vector<vector<double> >& vva, const double& value)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > ShrinkValuesExceptFirstColumn(const vector<vector<double> >& vva, const double& Fi) {
    //shrink Fis (usually means DOS Fis in DOS); Fi means probability
    vector<vector<double> > vvb=vva;
    for (uint i=0; i<vvb.size(); i++) {
      for (uint j=1; j<vvb.at(i).size();j++) {
	vvb.at(i).at(j)*=Fi;
      }
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// vector<vector<double> > aurostd::NormalizeAndSum3DVector(const vector<vector<vector<double> > >& vvva, const vector<vector<double> >& vFi)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > NormalizeAndSum3DVector(const vector<vector<vector<double> > >& vvva, const vector<double>& vFi) {
    //normalize DOS and sum
    if(vvva.size()!=vFi.size()) {cerr << "Vector sizes are not equal! Aborting!" << endl; exit(2);}
    vector<vector<double> > vvb, vv_tmp, vv_tmp_shrinked;
    vector<vector<vector<double> > > vvvc;
    double Fi;
    for (uint i=0; i<vvva.size();i++) {
      vv_tmp=vvva.at(i);
      Fi=vFi.at(i);
      vv_tmp_shrinked=aurostd::ShrinkValuesExceptFirstColumn(vv_tmp, Fi);
      vvvc.push_back(vv_tmp_shrinked);
    }
    vvb=aurostd::Sum3DVectorAndReduce2D(vvvc);
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva) {
    //The first column will not change! (For example, PDOS into TOTALPDOS)
    vector<vector<double> > vvtmp, vv_sum; 
    vv_sum=vvva.at(0);
    for (uint i=1; i<vvva.size();i++) {
      vvtmp=vvva.at(i);
      vv_sum=aurostd::Sum2DVectorExceptFirstColumn(vv_sum, vvtmp);
    }
    return vv_sum;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > Sum2DVectorExceptFirstColumn(const vector<vector<double> >& vva, const vector<vector<double> >& vvb) {
    if((vva.size()!=vvb.size()) && (vva.at(0).size() != vvb.at(0).size())) {cerr << "Error, two vectors! Aborting! " << endl; exit(2);}

    vector<vector<double> > vv_sum; vv_sum.resize(vva.size());
    for (uint i=0; i<vva.size(); i++) {
      int N=vva.at(i).size();
      vv_sum.at(i).resize(N);
    }

    for (uint i=0; i<vva.size();i++) {
      vv_sum[i][0]=vva.at(i).at(0);
      for (uint j=1; j<vva.at(i).size();j++) {
	vv_sum[i][j]=vva[i][j] + vvb[i][j];
      }
    }
    return vv_sum;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::ReduceVector(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > ReduceVector(const vector<vector<double> >& vva, const int& n) {
    //Pick up the first (begin from 0) and the nth column of 2D vector
    vector<vector<double> > vvb; vvb.clear();
    vector<double> vtmp;
    for (uint i=0; i<vva.size();i++) {
      vtmp.clear();
      vtmp.push_back(vva.at(i).at(0));
      vtmp.push_back(vva.at(i).at(n));
      vvb.push_back(vtmp);
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n) {
    //Calculate integration of vva, the 0st column is x0, x1..., the n column is y1, y2 ...
    //begin from 0
    vector<vector<double> > vvb=aurostd::ReduceVector(vva, n);
    return aurostd::CalculateIntegrate(vvb);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n, const double& Emin, const double& Emax) {
    //Calculate integration of vva, the 0st column is x0, x1..., the n column is y1, y2 ...
    //begin from 0
    vector<vector<double> > vvb=aurostd::ReduceVector(vva, n);
    return aurostd::CalculateIntegrate(vvb, Emin, Emax);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva) {
    double Emin=-100; double Emax=0.0; //default setting
    return aurostd::CalculateIntegrate(vva, Emin, Emax);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva, const double& Emin, const double& Emax) {
    //Integral function
    //format of vva: x0, y0; x1, y1; x2, y2
    double integral_result=0.0;
    double area_tmp =0.0;
    double xbeg, xend, ybeg, yend;
    for (uint i=0; i<vva.size()-1;i++) {
      xbeg=vva.at(i).at(0); xend=vva.at(i+1).at(0);
      ybeg=vva.at(i).at(1); yend=vva.at(i+1).at(1);
      if(xbeg >= Emin && xend <= Emax) {
	area_tmp=0.5*(ybeg + yend)*(xend - xbeg);
	integral_result += area_tmp;
      }
    }
    return integral_result;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::vector2string(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  string vector2string(const vector<vector<double> >& vva) {
    stringstream ss_vva; ss_vva.str(std::string());
    ss_vva << std::scientific;
    for (uint i=0; i<vva.size();i++) {
      for (uint j=0; j<vva.at(i).size();j++) {
	ss_vva << vva.at(i).at(j) << "   ";
      }
      ss_vva << endl;
    }
    return ss_vva.str();
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva) {
    double min=-10;  //default
    double max=10;
    return aurostd::FindMaxIn2DvectorExcept1stColumn(vva, min, max);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxIn2DvectorExcept1stColumn(const vector<vector<double>& vva, const double& min, const double& max)
// ***************************************************************************
namespace aurostd  {
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva, const double& min, const double& max) {
    double max_value=0.0;
    for (uint i=0; i<vva.size();i++) {
      double E_tmp=vva.at(i).at(0);
      if(E_tmp >= min && E_tmp <= max) {
	for (uint j=1; j<vva.at(i).size();j++) {
	  double db_tmp=vva.at(i).at(j);
	  if(abs(db_tmp) > max_value) max_value=abs(db_tmp);
	}
      }
    }
    return max_value;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxInTDOS(const vector<vector<double> >& vva, const double& min, const double& max)
// ***************************************************************************
namespace aurostd  {
  double FindMaxInTDOS(const vector<vector<double> >& vva, const double& min, const double& max) {
    double max_value=0.0;
    for (uint i=0; i<vva.size();i++) {
      double E_tmp=vva.at(i).at(0);
      if(E_tmp >= min && E_tmp <= max) {
	int column_max;
	if(vva.at(0).size()==3) column_max=2; //get rid of the sum of TDOS
	if(vva.at(0).size()==5) column_max=3;
	for (int j=1; j<column_max;j++) {
	  double db_tmp=vva.at(i).at(j);
	  if(abs(db_tmp) > max_value) max_value=db_tmp;
	}
      }
    }
    return max_value;
  }
} // namespace aurostd


namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(vector<uint>& uientries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters int/uint type of objects together by a delimiter
  // no point for double objects, faster to just do it on the spot with
  // setprecision,fixed, etc.
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();

    if (ientries.rows > 2) {
      for (int i =ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) {  //CO 180216 - added -1
          output << lDelim;
        } else if (i !=ientries.urows) {
          output << delim;
        }
      }
    } else {
      for (int i = ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) {  //CO 180216 - added -1
          output << mDelim;
        } else if (i != ientries.urows) {
          output << delim;
        }
      }
    }
    return output.str();
  }
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();

    if (ientries.size() > 2) {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries.at(i);
        if (i == ientries.size() - 2) {
          output << lDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries.at(i);
        if (i == ientries.size() - 2) {
          output << mDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter) {
    return joinWDelimiter(uientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    if (uientries.size() > 2) {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries.at(i);
        if (i == uientries.size() - 2) {
          output << lDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries.at(i);
        if (i == uientries.size() - 2) {
          output << mDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(vector<string>& _sentries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters string type of objects together by a delimiter
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries,
			const stringstream& delimiter) {
    return joinWDelimiter(_sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter) {
    stringstream output;
    vector<string> sentries;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    // go through once to eliminate empty strings
    for (uint i = 0; i < _sentries.size(); i++) {
      if (_sentries.at(i).length()) {
        sentries.push_back(_sentries.at(i));
      }
    }
    if (sentries.size() > 2) {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries.at(i);
        if (i == sentries.size() - 2) {
          output << lDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries.at(i);
        if (i == sentries.size() - 2) {
          output << mDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
} // namespace

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(deque<uint>& uientries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters int/uint type of objects together by a delimiter
  // no point for double objects, faster to just do it on the spot with
  // setprecision,fixed, etc.
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    if (ientries.size() > 2) {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries.at(i);
        if (i == ientries.size() - 2) {
          output << lDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries.at(i);
        if (i == ientries.size() - 2) {
          output << mDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter) {
    return joinWDelimiter(uientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    if (uientries.size() > 2) {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries.at(i);
        if (i == uientries.size() - 2) {
          output << lDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries.at(i);
        if (i == uientries.size() - 2) {
          output << mDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
}

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(deque<string>& _sentries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters string type of objects together by a delimiter
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter,
			const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter,
			const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter) {
    return joinWDelimiter(_sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter,
			const stringstream& m_delimiter,
			const stringstream& l_delimiter) {
    stringstream output;
    vector<string> sentries;  // no point working with deque
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    // go through once to eliminate empty strings
    for (uint i = 0; i < _sentries.size(); i++) {
      //DX - Should not be an "!"; we want it to push back if it has a length [OBSOLETE] if (!_sentries.at(i).length()) {
      if (_sentries.at(i).length()) {
        sentries.push_back(_sentries.at(i));
      }
    }
    if (sentries.size() > 2) {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries.at(i);
        if (i == sentries.size() - 2) {
          output << lDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries.at(i);
        if (i == sentries.size() - 2) {
          output << mDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
}

namespace aurostd {
  string wrapString(const string& input,const string& wrapper){return wrapString(input,wrapper,wrapper);}
  string wrapString(const string& input,const string& wrapper_start,const string& wrapper_end){
    if(input.empty()){return input;}
    return wrapper_start+input+wrapper_end;
  }
}

//DX 1/18/18 - START: XCOMPLEX TO JSON
namespace aurostd {
  //***************************************************************************//
  // aurostd::xcomplex2json
  //***************************************************************************//
  template<typename utype> string _xcomplex2json(xcomplex<utype>& number){
    string eendl="";
    bool roff=true; //round off
    stringstream sss;
    stringstream sscontent_json;
    vector<string> vcontent_json;
    // real
    sscontent_json << "\"real\":\"" << aurostd::utype2string(number.re,5,roff) << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    // imaginary
    sscontent_json << "\"imag\":\"" << aurostd::utype2string(number.im,5,roff) << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << eendl;
    return sss.str();    
  }
}

//Need to initalize 
namespace aurostd {
  string xcomplex2json(xcomplex<double>& number){ return _xcomplex2json(number); }
}

//DX 1/18/18 - END: XCOMPLEX TO JSON

//DX 8/3/17 - START: Matrix to JSON
namespace aurostd {
  //***************************************************************************//
  // aurostd::xmatDouble2String(xmatrix<double>& xmat_in)
  //***************************************************************************//
  // converts xmatrix<double> to json string
  // [OBSOLTE] string xmatDouble2String(const xmatrix<double>& xmat_in, bool roff){
  // [OBSOLTE]   stringstream output;
  // [OBSOLTE]   vector<string> rows;
  // [OBSOLTE]   for(uint i=1;i<(uint)xmat_in.rows+1;i++){
  // [OBSOLTE]     stringstream row;
  // [OBSOLTE]     xvector<double> xvec = xmat_in(i); //DX 8/22/17 - added roundoff
  // [OBSOLTE]     if(roff){ xvec = roundoff(xvec,1e-8);} //DX 8/22/17 - added roundoff
  // [OBSOLTE]     row << "[" << joinWDelimiter(xvecDouble2vecString(xvec),",") << "]";
  // [OBSOLTE]     rows.push_back(row.str());
  // [OBSOLTE]   }
  // [OBSOLTE]   output << joinWDelimiter(rows,",");
  // [OBSOLTE]   return output.str();
  // [OBSOLTE] }
  string xmatDouble2String(const xmatrix<double>& xmat_in, int precision, bool roff, double tol, char FORMAT){
    stringstream output;
    vector<string> rows;
    for(int i=xmat_in.urows;i<=xmat_in.urows;i++){
      stringstream row;
      xvector<double> xvec = xmat_in(i); //DX 8/22/17 - added roundoff
      //if(roff){ xvec = roundoff(xvec,tol);} //DX 8/22/17 - added roundoff
      row << "[" << joinWDelimiter(xvecDouble2vecString(xvec,precision,roff,tol,FORMAT),",") << "]";
      rows.push_back(row.str());
    }
    output << joinWDelimiter(rows,",");
    return output.str();
  }
}
//DX 8/3/17 - START: Matrix to END

namespace aurostd {
  //***************************************************************************//
  // aurostd::vecDouble2vecString(vector<double>& vin,int precision)
  //***************************************************************************//
  // converts vector<double> to vector<string> with precision
  // also works for xvectors and deques
  // [OBSOLTE] vector<string> vecDouble2vecString(const vector<double>& vin, bool roff) {
  // [OBSOLTE]   vector<string> vout;
  // [OBSOLTE]   for(uint i=0;i<vin.size();i++){
  // [OBSOLTE]     double tmp = vin.at(i); //DX 8/22/17 - add roundoff
  // [OBSOLTE]     if(roff){ tmp=scalar_roundoff(tmp,1e-8); } //DX 8/22/17 - add roundoff
  // [OBSOLTE]     vout.push_back(aurostd::utype2string(tmp)); //DX 8/22/17 - add roundoff
  // [OBSOLTE]   }
  // [OBSOLTE]   return vout;
  // [OBSOLTE] }
  vector<string> vecDouble2vecString(const vector<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    vector<string> vout;
    for(uint i=0;i<vin.size();i++){
      //double tmp = vin.at(i); //DX 8/22/17 - add roundoff
      //if(roff){ tmp=roundoff(tmp,tol); } //DX 8/22/17 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i],precision,roff,tol,FORMAT)); //DX 8/22/17 - add roundoff
    }
    return vout;
  }
  // [OBSOLTE] vector<string> xvecDouble2vecString(const xvector<double>& vin, bool roff) {
  // [OBSOLTE]   vector<string> vout;
  // [OBSOLTE]   for(uint i=1;i<(uint)vin.rows+1;i++){
  // [OBSOLTE]     double tmp = vin(i); //DX 8/22/17 - add roundoff
  // [OBSOLTE]    if(roff){ tmp=scalar_roundoff(tmp,1e-8); } //DX 8/22/17 - add roundoff
  // [OBSOLTE]     vout.push_back(aurostd::utype2string(tmp)); //DX 8/22/17 - add roundoff
  // [OBSOLTE]   }
  // [OBSOLTE]   return vout;
  // [OBSOLTE] }
  vector<string> xvecDouble2vecString(const xvector<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    vector<string> vout;
    for(int i=vin.lrows;i<=vin.urows;i++){
      //double tmp = vin(i); //DX 8/22/17 - add roundoff
      //if(roff){ tmp=roundoff(tmp,tol); } //DX 8/22/17 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i],precision,roff,tol,FORMAT)); //DX 8/22/17 - add roundoff
    }
    return vout;
  }
  // [OBSOLTE] deque<string> deqDouble2deqString(const deque<double>& vin, bool roff) {
  // [OBSOLTE]   deque<string> vout;
  // [OBSOLTE]   for(uint i=0;i<vin.size();i++){
  // [OBSOLTE]     double tmp = vin.at(i); //DX 8/22/17 - add roundoff
  // [OBSOLTE]     if(roff){ tmp=scalar_roundoff(tmp,1e-8); } //DX 8/22/17 - add roundoff
  // [OBSOLTE]     vout.push_back(aurostd::utype2string(tmp)); //DX 8/22/17 - add roundoff
  // [OBSOLTE]   }
  // [OBSOLTE]   return vout;
  // [OBSOLTE] }
  deque<string> deqDouble2deqString(const deque<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    deque<string> vout;
    for(uint i=0;i<vin.size();i++){
      //double tmp = vin.at(i); //DX 8/22/17 - add roundoff
      //if(roff){ tmp=roundoff(tmp,tol); } //DX 8/22/17 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i],precision,roff,tol,FORMAT)); //DX 8/22/17 - add roundoff
    }
    return vout;
  }
}

namespace aurostd {
  //***************************************************************************//
  // aurostd::wrapVecEntries(vector<string>& vin,string wrap)
  //***************************************************************************//
  // individually wraps entries of vector with specified string
  // converts <a,b,c> to <'a','b','c'>
  // also works for deques
  vector<string> wrapVecEntries(const vector<string>& vin,string wrap){
    return wrapVecEntries(vin,wrap,wrap);
  }
  vector<string> wrapVecEntries(const vector<string>& vin,string wrap_start,string wrap_end){
    vector<string> vout;
    for(uint i=0;i<vin.size();i++){
      if(vin.at(i).length()){
        vout.push_back(wrap_start+vin.at(i)+wrap_end);
      }
    }
    return vout;
  }
  deque<string> wrapDeqEntries(const deque<string>& vin,string wrap){
    return wrapDeqEntries(vin,wrap,wrap);
  }
  deque<string> wrapDeqEntries(const deque<string>& vin,string wrap_start,string wrap_end){
    deque<string> vout;
    for(uint i=0;i<vin.size();i++){
      if(vin.at(i).length()){
        vout.push_back(wrap_start+vin.at(i)+wrap_end);
      }
    }
    return vout;
  }
}

//base64 stuff
//CO - START
namespace aurostd {
  // ***************************************************************************
  // aurostd::isBase64(unsigned char c)
  // ***************************************************************************
  // determines if char is base64
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  //static inline bool isBase64(unsigned char c) {
  inline bool isBase64(unsigned char c) {
    return (isalnum(c) || (c == '+') || (c == '/'));
  }

  // ***************************************************************************
  // aurostd::base64Encoder(unsigned char const* bytes_to_encode, unsigned int in_len)
  // ***************************************************************************
  // encodes bytes to base64
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  std::string base64Encoder(unsigned char const* bytes_to_encode, unsigned int in_len) {
    std::string ret;
    int i = 0;
    int j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
      char_array_3[i++] = *(bytes_to_encode++);
      if (i == 3) {
        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for(i = 0; (i <4) ; i++) {
          ret += base64_chars[char_array_4[i]];
        }
        i = 0;
      }
    }

    if (i) {
      for(j = i; j < 3; j++) {
        char_array_3[j] = '\0';
      }

      char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
      char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
      char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
      char_array_4[3] = char_array_3[2] & 0x3f;

      for (j = 0; (j < i + 1); j++) {
        ret += base64_chars[char_array_4[j]];
      }

      while((i++ < 3)) {
        ret += '=';
      }
    }

    return ret;
  }

  // ***************************************************************************
  // aurostd::base64Decoder(std::string const& encoded_string)
  // ***************************************************************************
  // decodes base64 to bytes
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  std::string base64Decoder(std::string const& encoded_string) {
    int in_len = encoded_string.size();
    int i = 0;
    int j = 0;
    int in_ = 0;
    unsigned char char_array_4[4], char_array_3[3];
    std::string ret;

    while (in_len-- && ( encoded_string[in_] != '=') && isBase64(encoded_string[in_])) {
      char_array_4[i++] = encoded_string[in_]; in_++;
      if (i ==4) {
        for (i = 0; i <4; i++) {
          char_array_4[i] = base64_chars.find(char_array_4[i]);
        }

        char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (i = 0; (i < 3); i++) {
          ret += char_array_3[i];
        }
        i = 0;
      }
    }

    if (i) {
      for (j = i; j <4; j++) {
        char_array_4[j] = 0;
      }

      for (j = 0; j <4; j++) {
        char_array_4[j] = base64_chars.find(char_array_4[j]);
      }

      char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
      char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
      char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

      for (j = 0; (j < i - 1); j++) {
        ret += char_array_3[j];
      }
    }

    return ret;
  }

  // ***************************************************************************
  // aurostd::bin2base64(std::string& b_file, std::string& b64String)
  // ***************************************************************************
  // converts binary file to base64 string
  bool bin2base64(std::string& b_file, std::string& b64String) {
    stringstream output;
    if (!aurostd::FileExist(b_file)) {
      cerr << "ERROR - Binary file " << b_file << " does not exist!";
      return FALSE;
    }
    ifstream file(b_file.c_str(), std::ios::in | std::ios::binary );
    output << b64_encoder << file;
    b64String=output.str();
    return TRUE;
  }
  
  // ***************************************************************************
  // aurostd::base642bin(const std::string& b64String, std::string& b_file)
  // ***************************************************************************
  // converts base64 string to binary file
  bool base642bin(const std::string& b64String, std::string& b_file) {
    ofstream output;
    output.open(b_file.c_str(),std::ios::out | std::ios::binary);
    output << b64_decoder << b64String;
    output.flush();output.clear();output.close();
    return TRUE;
  }
  
  b64_encoder_proxy operator<<(std::ostream & os, b64_encoder_creator) {
    return b64_encoder_proxy(os);
  }

  b64_decoder_proxy operator<<(std::ostream & os, b64_decoder_creator) {
    return b64_decoder_proxy(os);
  }

}  // namespace aurostd
//CO - END

#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                                                                         *
// ***************************************************************************

