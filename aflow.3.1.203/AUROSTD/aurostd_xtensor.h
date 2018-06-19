// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************


// ---------------------------------------------------------------------------
// ------------------ implementation for template<class utype> xtensorX<utype>

#ifndef _AUROSTD_XTENSOR_H_
#define _AUROSTD_XTENSOR_H_

#define _AUROSTD_XTENSOR_DEFAULT_SIZE_ 3
#define _AUROSTD_XTENSOR_TOLERANCE_IDENTITY_ 1.0e-6
#define _AUROSTD_XTENSOR_TOLERANCE_ROUNDOFF_ 1.0e-6

#define _AUROSTD_XTENSOR3_INDEX_   3
#define _AUROSTD_XTENSOR4_INDEX_   4
#define _AUROSTD_XTENSOR5_INDEX_   5
#define _AUROSTD_XTENSOR6_INDEX_   6
#define _AUROSTD_XTENSOR7_INDEX_   7
#define _AUROSTD_XTENSOR8_INDEX_   8

// ---------------------------------------------------------------------------

namespace aurostd {
  // namespace aurostd
  int eijk(int,int,int) __xprototype;                       // eijk Ricci Tensor (function)
  int eijk(xvector<int>) __xprototype;                       // eijk Ricci Tensor (function)
  int estarijk(int,int,int) __xprototype;                       // eijk Ricci Tensor (function)
  int estarijk(xvector<int>) __xprototype;                       // eijk Ricci Tensor (function)
}

// ------------------------------------------------------------ class xtensor3
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    class xtensor3  {
  public:
    xtensor3 (int=3,int=3,int=3,                            // default constructor  
	      int=1,int=1,int=1);                           // default constructor  
    xtensor3 (const xtensor3<utype>&);                      // copy constructor  
    xtensor3<utype>& operator=(const xtensor3<utype>&);     // assignment
    ~xtensor3(void);                                        // default destructor  
    utype** operator[](int) const;                          // indicize i,j,k ONE
    utype& operator()(int,int,int) const;                   // indicize i,j,k ONE
    utype& operator()(int,int,int,bool) const;              //indicize boundary c
    // math operators
    xtensor3<utype>& operator +=(const xtensor3<utype>&);   //
    xtensor3<utype>& operator -=(const xtensor3<utype>&);   //
    // ostream operator
    //  friend ostream& operator<<<utype>(ostream&,const xtensor3<utype>&);
    int  index[_AUROSTD_XTENSOR3_INDEX_+1];                         // index  (from 1 to D)
    int lindex[_AUROSTD_XTENSOR3_INDEX_+1];                         // lindex (from 1 to D)
    int uindex[_AUROSTD_XTENSOR3_INDEX_+1];                         // uindex (from 1 to D)
    bool iscubic;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    utype*** corpus;
    bool isfloat,iscomplex;
    char size;
    long int tsize;
  };
}

// ------------------------------------------------------------ class xtensor4
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    class xtensor4  {
  public:
    xtensor4 (int=3,int=3,int=3,int=3,                      // default constructor  
	      int=1,int=1,int=1,int=1);                     // default constructor  
    xtensor4 (const xtensor4<utype>&);                      // copy constructor  
    xtensor4<utype>& operator=(const xtensor4<utype>&);     // assignment
    ~xtensor4(void);                                        // default destructor  
    utype*** operator[](int) const;                         // indicize i,j,k,l ONE
    utype& operator()(int,int,int,int) const;               // indicize i,j,k,l ONE
    utype& operator()(int,int,int,int,bool) const;          //indicize boundary c
    // math operators
    xtensor4<utype>& operator +=(const xtensor4<utype>&);   //
    xtensor4<utype>& operator -=(const xtensor4<utype>&);   //
    // ostream operator
    //  friend ostream& operator<<<utype>(ostream&,const xtensor4<utype>&);
    int  index[_AUROSTD_XTENSOR4_INDEX_+1];                         // index  (from 1 to D)
    int lindex[_AUROSTD_XTENSOR4_INDEX_+1];                         // lindex (from 1 to D)
    int uindex[_AUROSTD_XTENSOR4_INDEX_+1];                         // uindex (from 1 to D)
    bool iscubic;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    utype**** corpus;
    bool isfloat,iscomplex;
    char size;
    long int tsize;
  };
}

// ------------------------------------------------------------ class xtensor5
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    class xtensor5  {
  public:
    xtensor5 (int=3,int=3,int=3,int=3,int=3,                // default constructor  
	      int=1,int=1,int=1,int=1,int=1);               // default constructor  
    xtensor5 (const xtensor5<utype>&);                      // copy constructor  
    xtensor5<utype>& operator=(const xtensor5<utype>&);     // assignment
    ~xtensor5(void);                                        // default destructor  
    utype**** operator[](int) const;                        // indicize i,j,k,l,m ONE
    utype& operator()(int,int,int,int,int) const;           // indicize i,j,k,l,m ONE
    utype& operator()(int,int,int,int,int,bool) const;       //indicize boundary c
    // math operators
    xtensor5<utype>& operator +=(const xtensor5<utype>&);   //
    xtensor5<utype>& operator -=(const xtensor5<utype>&);   //
    // ostream operator
    //  friend ostream& operator<<<utype>(ostream&,const xtensor5<utype>&);
    int  index[_AUROSTD_XTENSOR5_INDEX_+1];                         // index  (from 1 to D)
    int lindex[_AUROSTD_XTENSOR5_INDEX_+1];                         // lindex (from 1 to D)
    int uindex[_AUROSTD_XTENSOR5_INDEX_+1];                         // uindex (from 1 to D)
    bool iscubic;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    utype***** corpus;
    bool isfloat,iscomplex;
    char size;
    long int tsize;
  };
}

// ------------------------------------------------------------ class xtensor6
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    class xtensor6  {
  public:
    xtensor6 (int=3,int=3,int=3,int=3,int=3,int=3,          // default constructor  
	      int=1,int=1,int=1,int=1,int=1,int=1);         // default constructor  
    xtensor6 (const xtensor6<utype>&);                      // copy constructor  
    xtensor6<utype>& operator=(const xtensor6<utype>&);     // assignment
    ~xtensor6(void);                                        // default destructor  
    utype***** operator[](int) const;                       // indicize i,j,k,l,m,n ONE
    utype& operator()(int,int,int,int,int,int) const;       // indicize i,j,k,l,m,n ONE
    utype& operator()(int,int,int,int,int,int,bool) const;  //indicize boundary c
    // math operators
    xtensor6<utype>& operator +=(const xtensor6<utype>&);   //
    xtensor6<utype>& operator -=(const xtensor6<utype>&);   //
    // ostream operator
    //  friend ostream& operator<<<utype>(ostream&,const xtensor6<utype>&);
    int  index[_AUROSTD_XTENSOR6_INDEX_+1];                         // index  (from 1 to D)
    int lindex[_AUROSTD_XTENSOR6_INDEX_+1];                         // lindex (from 1 to D)
    int uindex[_AUROSTD_XTENSOR6_INDEX_+1];                         // uindex (from 1 to D)
    bool iscubic;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    utype****** corpus;
    bool isfloat,iscomplex;
    char size;
    long int tsize;
  };
}

// ------------------------------------------------------------ class xtensor7
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    class xtensor7  {
  public:
    xtensor7 (int=3,int=3,int=3,int=3,int=3,int=3,int=3,    // default constructor  
	      int=1,int=1,int=1,int=1,int=1,int=1,int=1);   // default constructor  
    xtensor7 (const xtensor7<utype>&);                      // copy constructor  
    xtensor7<utype>& operator=(const xtensor7<utype>&);     // assignment
    ~xtensor7(void);                                        // default destructor  
    utype****** operator[](int) const;                      // indicize i,j,k,l,m,n ONE
    utype& operator()(int,int,int,int,int,int,int) const;   // indicize i,j,k,l,m,n ONE
    utype& operator()(int,int,int,int,int,int,int,bool) const;  //indicize boundary c
    // math operators
    xtensor7<utype>& operator +=(const xtensor7<utype>&);   //
    xtensor7<utype>& operator -=(const xtensor7<utype>&);   //
    // ostream operator
    //  friend ostream& operator<<<utype>(ostream&,const xtensor7<utype>&);
    int  index[_AUROSTD_XTENSOR7_INDEX_+1];                         // index  (from 1 to D)
    int lindex[_AUROSTD_XTENSOR7_INDEX_+1];                         // lindex (from 1 to D)
    int uindex[_AUROSTD_XTENSOR7_INDEX_+1];                         // uindex (from 1 to D)
    bool iscubic;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    utype******* corpus;
    bool isfloat,iscomplex;
    char size;
    long int tsize;
  };
}

// ------------------------------------------------------------ class xtensor8
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    class xtensor8  {
  public:
    xtensor8 (int=3,int=3,int=3,int=3,int=3,int=3,int=3,int=3, //default constructor  
	      int=1,int=1,int=1,int=1,int=1,int=1,int=1,int=1);//default constructor  
    xtensor8 (const xtensor8<utype>&);                      // copy constructor  
    xtensor8<utype>& operator=(const xtensor8<utype>&);     // assignment
    ~xtensor8(void);                                        // default destructor  
    utype******* operator[](int) const;                     // indicize i,j,k,l,m,n ONE
    utype& operator()(int,int,int,int,int,int,int,int) const; // indicize i,j,k,l,m,n ONE
    utype& operator()(int,int,int,int,int,int,int,int,bool) const; //indicize boundary c
    // math operators
    xtensor8<utype>& operator +=(const xtensor8<utype>&);   //
    xtensor8<utype>& operator -=(const xtensor8<utype>&);   //
    // ostream operator
    //  friend ostream& operator<<<utype>(ostream&,const xtensor8<utype>&);
    int  index[_AUROSTD_XTENSOR8_INDEX_+1];                         // index  (from 1 to D)
    int lindex[_AUROSTD_XTENSOR8_INDEX_+1];                         // lindex (from 1 to D)
    int uindex[_AUROSTD_XTENSOR8_INDEX_+1];                         // uindex (from 1 to D)
    bool iscubic;
    // operations
    void set(const utype&);
    void reset(void);
    void clear(void);
  private:
    utype******** corpus;
    bool isfloat,iscomplex;
    char size;
    long int tsize;
  };
}

// -------------------------------------------------------------------- debug
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    void xtensor3debug(const xtensor3<utype>& t,const string& str) __xprototype;

  template<class utype>
    void xtensor4debug(const xtensor4<utype>& t,const string& str) __xprototype;

  template<class utype>
    void xtensor5debug(const xtensor5<utype>& t,const string& str) __xprototype;

  template<class utype>
    void xtensor6debug(const xtensor6<utype>& t,const string& str) __xprototype;

  template<class utype>
    void xtensor7debug(const xtensor7<utype>& t,const string& str) __xprototype;

  template<class utype>
    void xtensor8debug(const xtensor8<utype>& t,const string& str) __xprototype;
}

// -------------------------------------------------------------------- allocate
namespace aurostd {  
  // namespace aurostd
  template<class utype>
    void allocate_xtensor3corpus(utype*** &corpus,int* lindex,int* uindex,int* index);
  
  template<class utype>
    void allocate_xtensor4corpus(utype**** &corpus,int* lindex,int* uindex,int* index);
  
  template<class utype>
    void allocate_xtensor5corpus(utype***** &corpus,int* lindex,int* uindex,int* index);
  
  template<class utype>
    void allocate_xtensor6corpus(utype****** &corpus,int* lindex,int* uindex,int* index);
  
  template<class utype>
    void allocate_xtensor7corpus(utype******* &corpus,int* lindex,int* uindex,int* index);
  
  template<class utype>
    void allocate_xtensor8corpus(utype******** &corpus,int* lindex,int* uindex,int* index);
}
  
// --------------------------------------------------------------------------
// ---------------------- primitives for template<class utype> xtensorX<utype>

// ------------------------------------------------- unary and binary operators
namespace aurostd {  
  // namespace aurostd

  template<class utype> xtensor3<utype>
    operator+(const xtensor3<utype>&,const xtensor3<utype>&) __xprototype;

  template<class utype> xtensor4<utype>
    operator+(const xtensor4<utype>&,const xtensor4<utype>&) __xprototype;

  template<class utype> xtensor5<utype>
    operator+(const xtensor5<utype>&,const xtensor5<utype>&) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator+(const xtensor6<utype>&,const xtensor6<utype>&) __xprototype;
  
  template<class utype> xtensor7<utype>
    operator+(const xtensor7<utype>&,const xtensor7<utype>&) __xprototype;
  
  template<class utype> xtensor8<utype>
    operator+(const xtensor8<utype>&,const xtensor8<utype>&) __xprototype;
  


  template<class utype> xtensor3<utype>
    operator-(const xtensor3<utype>&,const xtensor3<utype>&) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator-(const xtensor4<utype>&,const xtensor4<utype>&) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator-(const xtensor5<utype>&,const xtensor5<utype>&) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator-(const xtensor6<utype>&,const xtensor6<utype>&) __xprototype;
  
  template<class utype> xtensor7<utype>
    operator-(const xtensor7<utype>&,const xtensor7<utype>&) __xprototype;
  
  template<class utype> xtensor8<utype>
    operator-(const xtensor8<utype>&,const xtensor8<utype>&) __xprototype;
  


  template<class utype> xtensor3<utype>
    operator/(const xtensor3<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator/(const xtensor4<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator/(const xtensor5<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator/(const xtensor6<utype>&,const utype) __xprototype;

  template<class utype> xtensor7<utype>
    operator/(const xtensor7<utype>&,const utype) __xprototype;

  template<class utype> xtensor8<utype>
    operator/(const xtensor8<utype>&,const utype) __xprototype;

  

  template<class utype> xtensor3<utype>
    operator*(const xtensor3<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator*(const xtensor4<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator*(const xtensor5<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator*(const xtensor6<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor7<utype>
    operator*(const xtensor7<utype>&,const utype) __xprototype;
  
  template<class utype> xtensor8<utype>
    operator*(const xtensor8<utype>&,const utype) __xprototype;
  


  template<class utype> xtensor3<utype>
    operator*(const utype,const xtensor3<utype>&) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator*(const utype,const xtensor4<utype>&) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator*(const utype,const xtensor5<utype>&) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator*(const utype,const xtensor6<utype>&) __xprototype;
  
  template<class utype> xtensor7<utype>
    operator*(const utype,const xtensor7<utype>&) __xprototype;
  
  template<class utype> xtensor8<utype>
    operator*(const utype,const xtensor8<utype>&) __xprototype;
  


  template<class utype>
    ostream& operator<<(ostream&,const xtensor3<utype>&) __xprototype;
  
  template<class utype>
    ostream& operator<<(ostream&,const xtensor4<utype>&) __xprototype;
  
  template<class utype>
    ostream& operator<<(ostream&,const xtensor5<utype>&) __xprototype;
  
  template<class utype>
    ostream& operator<<(ostream&,const xtensor6<utype>&) __xprototype;

  template<class utype>
    ostream& operator<<(ostream&,const xtensor7<utype>&) __xprototype;

  template<class utype>
    ostream& operator<<(ostream&,const xtensor8<utype>&) __xprototype;



  template<class utype> xtensor3<utype>
    operator+(const xtensor3<utype>&) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator+(const xtensor4<utype>&) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator+(const xtensor5<utype>&) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator+(const xtensor6<utype>&) __xprototype;

  template<class utype> xtensor7<utype>
    operator+(const xtensor7<utype>&) __xprototype;

  template<class utype> xtensor8<utype>
    operator+(const xtensor8<utype>&) __xprototype;


  
  template<class utype> xtensor3<utype>
    operator-(const xtensor3<utype>&) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator-(const xtensor4<utype>&) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator-(const xtensor5<utype>&) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator-(const xtensor6<utype>&) __xprototype;
  
  template<class utype> xtensor7<utype>
    operator-(const xtensor7<utype>&) __xprototype;
  
  template<class utype> xtensor8<utype>
    operator-(const xtensor8<utype>&) __xprototype;

}
  
// ---------------------------------------------------- operators and functions
namespace aurostd {  
  // namespace aurostd

  template<class utype> xtensor3<utype>
    operator+(const xtensor3<utype>& a,const xtensor3<utype>& b) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator+(const xtensor4<utype>& a,const xtensor4<utype>& b) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator+(const xtensor5<utype>& a,const xtensor5<utype>& b) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator+(const xtensor6<utype>& a,const xtensor6<utype>& b) __xprototype;

  template<class utype> xtensor7<utype>
    operator+(const xtensor7<utype>& a,const xtensor7<utype>& b) __xprototype;

  template<class utype> xtensor8<utype>
    operator+(const xtensor8<utype>& a,const xtensor8<utype>& b) __xprototype;


  
  template<class utype> xtensor3<utype>
    operator-(const xtensor3<utype>& a,const xtensor3<utype>& b) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator-(const xtensor4<utype>& a,const xtensor4<utype>& b) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator-(const xtensor5<utype>& a,const xtensor5<utype>& b) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator-(const xtensor6<utype>& a,const xtensor6<utype>& b) __xprototype;
  
  template<class utype> xtensor7<utype>
    operator-(const xtensor7<utype>& a,const xtensor7<utype>& b) __xprototype;
  
  template<class utype> xtensor8<utype>
    operator-(const xtensor8<utype>& a,const xtensor8<utype>& b) __xprototype;
  


  template<class utype> xtensor3<utype>
    operator*(const utype s,const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator*(const utype s,const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator*(const utype s,const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator*(const utype s,const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    operator*(const utype s,const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    operator*(const utype s,const xtensor8<utype>& a) __xprototype;


  
  template<class utype> xtensor3<utype>
    operator*(const xtensor3<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator*(const xtensor4<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator*(const xtensor5<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator*(const xtensor6<utype>& a,const utype s) __xprototype;

  template<class utype> xtensor7<utype>
    operator*(const xtensor7<utype>& a,const utype s) __xprototype;

  template<class utype> xtensor8<utype>
    operator*(const xtensor8<utype>& a,const utype s) __xprototype;


  
  template<class utype> xtensor3<utype>
    operator/(const xtensor3<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor4<utype>
    operator/(const xtensor4<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor5<utype>
    operator/(const xtensor5<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor6<utype>
    operator/(const xtensor6<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor7<utype>
    operator/(const xtensor7<utype>& a,const utype s) __xprototype;
  
  template<class utype> xtensor8<utype>
    operator/(const xtensor8<utype>& a,const utype s) __xprototype;
  


  template<class utype> xtensor3<utype>
    sign(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    sign(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    sign(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    sign(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    sign(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    sign(const xtensor8<utype>& a) __xprototype;



  template<class utype> xtensor3<utype>
    nint(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    nint(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    nint(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    nint(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    nint(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    nint(const xtensor8<utype>& a) __xprototype;



  template<class utype> xtensor3<utype>
    abs(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    abs(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    abs(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    abs(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    abs(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    abs(const xtensor8<utype>& a) __xprototype;


  template<class utype> xtensor3<utype>
    floor(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    floor(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    floor(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    floor(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    floor(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    floor(const xtensor8<utype>& a) __xprototype;


  template<class utype> xtensor3<utype>
    ceil(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    ceil(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    ceil(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    ceil(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    ceil(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    ceil(const xtensor8<utype>& a) __xprototype;


  template<class utype> xtensor3<utype>
    round(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    round(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    round(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    round(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    round(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    round(const xtensor8<utype>& a) __xprototype;



  template<class utype> xtensor3<utype>
    trunc(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> xtensor4<utype>
    trunc(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> xtensor5<utype>
    trunc(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> xtensor6<utype>
    trunc(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    trunc(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    trunc(const xtensor8<utype>& a) __xprototype;

  
  
  template<class utype> void
    reset(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> void
    reset(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> void
    reset(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> void
    reset(const xtensor6<utype>& a) __xprototype;

  template<class utype> void
    reset(const xtensor7<utype>& a) __xprototype;

  template<class utype> void
    reset(const xtensor8<utype>& a) __xprototype;


  
 template<class utype> void
    clear(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> void
    clear(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> void
    clear(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> void
    clear(const xtensor6<utype>& a) __xprototype;

  template<class utype> void
    clear(const xtensor7<utype>& a) __xprototype;

  template<class utype> void
    clear(const xtensor8<utype>& a) __xprototype;


  
  template<class utype> void
    set(const xtensor3<utype>& a,const utype& s) __xprototype;
  
  template<class utype> void
    set(const xtensor4<utype>& a,const utype& s) __xprototype;
  
  template<class utype> void
    set(const xtensor5<utype>& a,const utype& s) __xprototype;
  
  template<class utype> void
    set(const xtensor6<utype>& a,const utype& s) __xprototype;

  template<class utype> void
    set(const xtensor7<utype>& a,const utype& s) __xprototype;

  template<class utype> void
    set(const xtensor8<utype>& a,const utype& s) __xprototype;



  template<class utype> utype
    sum(const xtensor3<utype>& a) __xprototype;

  template<class utype> utype
    sum(const xtensor4<utype>& a) __xprototype;

  template<class utype> utype
    sum(const xtensor5<utype>& a) __xprototype;

  template<class utype> utype
    sum(const xtensor6<utype>& a) __xprototype;

  template<class utype> utype
    sum(const xtensor7<utype>& a) __xprototype;

  template<class utype> utype
    sum(const xtensor8<utype>& a) __xprototype;



  template<class utype> utype
    min(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> utype
    min(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> utype
    min(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> utype
    min(const xtensor6<utype>& a) __xprototype;
  
  template<class utype> utype
    min(const xtensor7<utype>& a) __xprototype;
  
  template<class utype> utype
    min(const xtensor8<utype>& a) __xprototype;
  


  template<class utype> utype
    min(const xtensor3<utype>& a,int& i,int& j,int& k) __xprototype;
  
  template<class utype> utype
    min(const xtensor4<utype>& a,int& i,int& j,int& k,int& l) __xprototype;
  
  template<class utype> utype
    min(const xtensor5<utype>& a,int& i,int& j,int& k,int& l,int& m) __xprototype;
  
  template<class utype> utype
    min(const xtensor6<utype>& a,int& i,int& j,int& k,int& l,int& m,int& n) __xprototype;

  template<class utype> utype
    min(const xtensor7<utype>& a,int& i,int& j,int& k,int& l,int& m,int& n,int& o) __xprototype;

  template<class utype> utype
    min(const xtensor8<utype>& a,int& i,int& j,int& k,int& l,int& m,int& n,int& o,int& p) __xprototype;



  template<class utype> utype
    max(const xtensor3<utype>& a) __xprototype;
  
  template<class utype> utype
    max(const xtensor4<utype>& a) __xprototype;
  
  template<class utype> utype
    max(const xtensor5<utype>& a) __xprototype;
  
  template<class utype> utype
    max(const xtensor6<utype>& a) __xprototype;
  
  template<class utype> utype
    max(const xtensor7<utype>& a) __xprototype;
  
  template<class utype> utype
    max(const xtensor8<utype>& a) __xprototype;
  


  template<class utype> utype
    max(const xtensor3<utype>& a,int& i,int& j,int& k) __xprototype;
  
  template<class utype> utype
    max(const xtensor4<utype>& a,int& i,int& j,int& k,int& l) __xprototype;
  
  template<class utype> utype
    max(const xtensor5<utype>& a,int& i,int& j,int& k,int& l,int& m) __xprototype;
  
  template<class utype> utype
    max(const xtensor6<utype>& a,int& i,int& j,int& k,int& l,int& m,int& n) __xprototype;

  template<class utype> utype
    max(const xtensor7<utype>& a,int& i,int& j,int& k,int& l,int& m,int& n,int& o) __xprototype;

  template<class utype> utype
    max(const xtensor8<utype>& a,int& i,int& j,int& k,int& l,int& m,int& n,int& o,int& p) __xprototype;



  template<class utype> utype
    trace(const xtensor3<utype>& a) __xprototype;

  template<class utype> utype
    trace(const xtensor4<utype>& a) __xprototype;

  template<class utype> utype
    trace(const xtensor5<utype>& a) __xprototype;

  template<class utype> utype
    trace(const xtensor6<utype>& a) __xprototype;

  template<class utype> utype
    trace(const xtensor7<utype>& a) __xprototype;

  template<class utype> utype
    trace(const xtensor8<utype>& a) __xprototype;



  template<class utype> xtensor3<utype>
    identity(const xtensor3<utype>& a) __xprototype;

  template<class utype> xtensor4<utype>
    identity(const xtensor4<utype>& a) __xprototype;

  template<class utype> xtensor5<utype>
    identity(const xtensor5<utype>& a) __xprototype;

  template<class utype> xtensor6<utype>
    identity(const xtensor6<utype>& a) __xprototype;

  template<class utype> xtensor7<utype>
    identity(const xtensor7<utype>& a) __xprototype;

  template<class utype> xtensor8<utype>
    identity(const xtensor8<utype>& a) __xprototype;

}

// **************************************************************************

#endif

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************
