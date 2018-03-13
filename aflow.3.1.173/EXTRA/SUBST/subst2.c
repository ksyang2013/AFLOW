/******************************************************************
* (C) 2000 - Stefano Curtarolo                                    *
* send bugs to auro@mit.edu                                       *
******************************************************************/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#define FLAG char
#define REL "beta 0.99.2.pl1"

int main(int argc, char **argv){

  /* declarations */
  char *in_file,*out_file,ctemp[1024];

  FILE *in_file_pointer,*out_file_pointer;
  char *in_string,*out_string;
  unsigned int in_string_lenght,out_string_lenght,nfiles,n;
  FLAG verbose,found,backup;
  char c;
  unsigned int i=0,ii,j;
  
  /*from inputs*/
  backup=0;
  verbose=0;
//  for(i=70;i<=85;i++) printf("%i %c \n",i,i); exit(1);
  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */
//  printf("number: %i \n",argc);
  found=0;
  if (argc<=3) {
    printf("subst2 version %s\n",REL);
    printf("Usage: subst2 string1 string2 inputfile1 [inputfile2 ...]\n");
    printf("       replace str1 with str2 deleting len(str2) chars. \n");
    printf("\n");
    printf("-h, --help            this help\n");
    printf("-v, --verbose         verbose\n");
    printf("-V, --version         version\n");
    printf("\n");
    printf("Send bugs to stefano@mit.edu\n");
    exit(0);
    }
  if (argc>=4) {
    nfiles=argc-3;
    in_string=(char*) malloc(strlen(argv[1])+1);
    if (in_string==NULL) {printf("Subst memory error");exit(1);}
    strcpy(in_string,(char*) argv[1]);
//    printf("%s ",in_string);

    out_string=(char*) malloc(strlen(argv[2])+1);
    if (out_string==NULL) {
       printf("Subst memory error");
       exit(1);
    }
    strcpy(out_string,(char*) argv[2]);
//    printf("%s ",out_string);
  }

  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */

  for(n=0;n<nfiles;n++) {
    in_file=(char*) malloc(strlen(argv[n+3])+1);
    out_file=(char*) malloc(strlen(argv[n+3])+2);
    if (in_file==NULL || out_file==NULL) {
       printf("Subst memory error");
       exit(1);
    }
    strcpy(in_file,(char*) argv[n+3]);
    strcpy(out_file,(char*) argv[n+3]);
    out_file[strlen(argv[n+3])]=126;
    out_file[strlen(argv[n+3])+1]=0;

    if (verbose){ 
      printf("Subst message: substituting \"%s\" with \"%s\" in file %s\n",
             in_string,out_string,in_file);
    }

  /* ------------------------------------------- initialization */
    in_file_pointer=fopen(in_file,"r");
    out_file_pointer=fopen(out_file,"w");
    if (in_file_pointer==NULL) {
      printf("Subst error: %s not found\n",in_file);
      exit(1);
    }
    if (out_file_pointer==NULL) {
      printf("Subst error: error opening %s \n",out_file);
      exit(1);
    }
  /* --------------------------------------------------- backup */
 
    c=(char) fgetc(in_file_pointer);
    while (c!=EOF){ 
       fprintf(out_file_pointer,"%c",c);
       c=(char) fgetc(in_file_pointer);
    } 
    if (fclose(in_file_pointer)) {
      printf("Subst error: error closing %s \n",in_file);
      exit(1);
    }
    if (fclose(out_file_pointer)) {
      printf("Subst error: error closing %s \n",out_file);
      exit(1);
    }
    
    in_file_pointer=fopen(out_file,"r");
    out_file_pointer=fopen(in_file,"w");
    if (in_file_pointer==NULL) {
       printf("Subst error: %s not found\n",in_file);
       exit(1);
    }
    if (out_file_pointer==NULL) {
       printf("Subst error: error opening %s \n",out_file);
       exit(1);
    }

  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */

    found=1;
    in_string_lenght=strlen(in_string);
    out_string_lenght=strlen(out_string);

    c=(char) fgetc(in_file_pointer);
    while (c!=EOF) {
      found=1;
      i=0;
      while (found && c!=EOF){
        ctemp[i]=c;
        if(c!=in_string[i])
	   found=0;
        else 
	   i=i+1;     
        if (!found) 
	   for(ii=0;ii<=i;ii++) {
	      fprintf(out_file_pointer,"%c",ctemp[ii]);
              if (verbose) printf("%c",ctemp[ii]);
	}
        if (found && i==in_string_lenght){
	   fprintf(out_file_pointer,"%s",out_string);
	   if (verbose) printf("%s",out_string);
	   found=0;
	   for (j=1;j<strlen(out_string);j++)
	     c=(char) fgetc(in_file_pointer);
        }
        c=(char) fgetc(in_file_pointer);
      }
    }

  /* ----------------------------------------- close everything */
    if (fclose(in_file_pointer)) {
      printf("Subst error: error closing %s \n",in_file);
      exit(1);
    }
    if (fclose(out_file_pointer)) {
      printf("Subst error: error closing %s \n",out_file);
      exit(1);
    }
    free(in_file);
    free(out_file);
  }
  free(in_string);
  free(out_string);
  exit(0);
}

  /* ---------------------------------------------------------- */
  /* ---------------------------------------------------------- */


