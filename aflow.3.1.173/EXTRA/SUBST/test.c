/******************************************************************
* (C) 1999 - Stefano Curtarolo                                    *
* send bugs to auro@mit.edu                                       *
******************************************************************/

#include<stdio.h>
#include<string.h>
#define FLAG char
#define REL "beta 0.99.2.pl1"

char char_good(char charin){
  if(charin>=48 && charin<=57) return (char) 1; 
  if(charin>=40 && charin<=41) return (char) 1; 
  if(charin==43 || charin==45) return (char) 1; 
  if(charin==20 || charin==95) return (char) 1; 
  return (char) 0; 
};

char * extract_number(char * string_in) {
  static char string_out[1024];
  int i,j;
  for (i=0,j=0;i<=strlen(string_in);i++) 
    if(char_good(string_in[i])) {
      string_out[j++]=string_in[i];
      if(i<=strlen(string_in) && !char_good(string_in[i+1]))
	i=strlen(string_in)+10;
    }
  string_out[j]=0;
  return string_out;
}

int main(int argc, char **argv){

  /* declarations */
  char *in_file,*out_file;

  FILE *in_file_pointer,*out_file_pointer;
  char *in_string,*out_string;
  unsigned int in_string_lenght,out_string_lenght,nfiles,n;
  FLAG verbose,found,backup;
  char c;
  unsigned int i=0,ii;
  
  /*from inputs*/
// for(c=20;c<=155;c++) printf("%c %i \n",c,c);
// exit(1);

  if(argc==1) {
    printf("Extract photos for AUROTRAINS \n");
    printf("need at least one file name \n");
    exit(1);
  } else {
    for (i=0;i<argc;i++) 
      printf("%s",extract_number(argv[i]));
  }
 // printf("%s \n",argv[0]);
 // printf("%s \n",argv[1]);

}


  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */


