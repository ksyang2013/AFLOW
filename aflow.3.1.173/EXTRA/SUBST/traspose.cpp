/******************************************************************
* (C) 1999 - Stefano Curtarolo                                    *
* send bugs to auro@mit.edu                                       *
******************************************************************/

#include<iostream.h>
#include"xlibs++/xmath.h"
#include<string.h>

#define FLAG char
#define REL "beta 0.99.2.pl1"

int main(int argc, char **argv){

  /* declarations */
  char ctemp[100];
  char cc[1],c;

  FILE *in_file_pointer;
  unsigned int in_file_lenght,n;
  unsigned int i,j,imax,jmax,eol,eof;
  double x;

  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */
  if (argc<=1) {
    cerr << "traspose version " << REL << endl;
    cerr << "Usage: traspose inputfile " << endl;
    cerr << "       traspose matrix-file.  " << endl;
    cerr << endl;
    cerr << "Send bugs to stefano@mit.edu " << endl;
    exit(1);
  }
  //  cerr << argv[1] << endl;
  
  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */
  
  /* ------------------------------------------- initialization */
  in_file_pointer=fopen(argv[1],"r");
  if (in_file_pointer==NULL) {
    cerr << "traspose error:" << argv[1] << " not found" << endl;
    exit(0);
  }
  imax=0;jmax=0;
  eol=0;eof=0;
  for(;!eof;) {
    fscanf(in_file_pointer,"%c",cc);c=*cc;
    if(c==10) { 
      imax++;fscanf(in_file_pointer,"%c",cc);c=*cc;
      if(c==10) eof=1;
    } 
    if(c==0) eof=1; 
  }
  //  cout << imax << endl; 
   if (fclose(in_file_pointer)) {
    cerr << "traspose error: error closing " << argv[1]  << endl;
    exit(0);
  }

  /* ---------------------------------------------------------- */

  in_file_pointer=fopen(argv[1],"r");
  if (in_file_pointer==NULL) {
    cerr << "traspose error:" << argv[1] << " not found" << endl;
    exit(0);
  }
  for(jmax=0;fscanf(in_file_pointer,"%s",ctemp)!=EOF;jmax++);
  jmax/=imax;
  // cout << jmax << endl;
  if (fclose(in_file_pointer)) {
    cerr << "traspose error: error closing " << argv[1]  << endl;
    exit(0);
  }

  /* ---------------------------------------------------------- */

  in_file_pointer=fopen(argv[1],"r");
  if (in_file_pointer==NULL) {
    cerr << "traspose error:" << argv[1] << " not found" << endl;
    exit(0);
  }

  matrix<double> mat(imax,jmax);
  for (i=1;i<=imax;i++)
    for (j=1;j<=jmax;j++) {
      fscanf(in_file_pointer,"%s",ctemp);
      mat(i,j)=atof(ctemp);
    }
  
  cout << trasp(mat) << endl;

  /* ----------------------------------------- close everything */
  if (fclose(in_file_pointer)) {
    cerr << "traspose error: error closing " <<argv[1]  << endl;
    exit(0);
  }
  exit(1);
}

  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */


