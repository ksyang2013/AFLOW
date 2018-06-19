// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

// ****************************************************************************
//#warning "xtensor template TFUNC"

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function TFUNC xtensor3<>
  xtensor3<utype> TFUNC(const xtensor3<utype>& a) {
    aurostd::xtensor3debug(a,"function TFUNC");
    xtensor3<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          c[i][j][k]=(utype) TFUNC(a[i][j][k]);
    return (xtensor3<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function TFUNC xtensor4<>
  xtensor4<utype> TFUNC(const xtensor4<utype>& a) {
    aurostd::xtensor4debug(a,"function TFUNC");
    xtensor4<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            c[i][j][k][l]=(utype) TFUNC(a[i][j][k][l]);
    return (xtensor4<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function TFUNC xtensor5<>
  xtensor5<utype> TFUNC(const xtensor5<utype>& a) {
    aurostd::xtensor5debug(a,"function TFUNC");
    xtensor5<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              c[i][j][k][l][m]=(utype) TFUNC(a[i][j][k][l][m]);
    return (xtensor5<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function TFUNC xtensor6<>
  xtensor6<utype> TFUNC(const xtensor6<utype>& a) {
    aurostd::xtensor6debug(a,"function TFUNC");
    xtensor6<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                c[i][j][k][l][m][n]=(utype) TFUNC(a[i][j][k][l][m][n]);
    return (xtensor6<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function TFUNC xtensor7<>
  xtensor7<utype> TFUNC(const xtensor7<utype>& a) {
    aurostd::xtensor7debug(a,"function TFUNC");
    xtensor7<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  c[i][j][k][l][m][n][o]=(utype) TFUNC(a[i][j][k][l][m][n][o]);
    return (xtensor7<utype>) c;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>                            // function TFUNC xtensor8<>
  xtensor8<utype> TFUNC(const xtensor8<utype>& a) {
    aurostd::xtensor8debug(a,"function TFUNC");
    xtensor8<utype> c(a);
    for(int i=a.lindex[1];i<=a.uindex[1];i++)
      for(int j=a.lindex[2];j<=a.uindex[2];j++)
        for(int k=a.lindex[3];k<=a.uindex[3];k++)
          for(int l=a.lindex[4];l<=a.uindex[4];l++)
            for(int m=a.lindex[5];m<=a.uindex[5];m++)
              for(int n=a.lindex[6];n<=a.uindex[6];n++)
                for(int o=a.lindex[7];o<=a.uindex[7];o++)
                  for(int p=a.lindex[8];p<=a.uindex[8];p++)
                    c[i][j][k][l][m][n][o][p]=(utype) TFUNC(a[i][j][k][l][m][n][o][p]);
    return (xtensor8<utype>) c;
  }
}


// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

