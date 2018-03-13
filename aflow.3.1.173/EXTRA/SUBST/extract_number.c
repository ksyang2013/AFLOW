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

int main(int argc, char **argv){

  /* declarations */
  char *in_file,*out_file,ctemp[1024];

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
    printf("ciao");
    exit(1);
  } else {
    for (i=0;i<=strlen(argv[1]);i++) {
      if(char_good(argv[1][i])) {
	printf("%c",argv[1][i]);
	if(i<=strlen(argv[1]) && !char_good(argv[1][i+1]))
	  i=strlen(argv[1])+10;
      }
    }
  }
 // printf("%s \n",argv[0]);
 // printf("%s \n",argv[1]);

}


  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */


