#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

int main(int argc,char **argv)
{

	if ( argc != 3 ){
		printf("%s input_file_name output_file_name\n", argv[0]);
		exit(-1);
	}
  /* Open the input file given from CLI for input */
  FILE * Fin= fopen(argv[1], "r");
  /* Open the output file given from CLI for output */
  FILE * Fout= fopen(argv[2], "w");
   
	if ( !Fin ){
		printf("File %s could not be opened!\n", argv[1]);
		exit(-1);
	}
	
	if ( !Fout ){
		printf("File %s could not be opened!\n", argv[2]);
		exit(-1);
	}
	
  int i;
  int Line;
  int GC=0;
  int AT=0;
  int max_lines = 0;
  char ch;
  /* Count the number of lines in the file */
  while ((ch = fgetc(Fin)) != EOF)
  {
    if (ch == '\n')
		max_lines++;
  }

  fseek(Fin, 0, SEEK_SET);
  
   /* malloc a 2d array with max_lines (lines of the file) as the rows*/
  char ** buffer;
  buffer=(char**)malloc(sizeof(char*)*max_lines );


  size_t len = 0;

  /* Read line-by-line the file and store each in the array named buffer */
  for(Line=0;Line<max_lines;Line++)
	getline(&buffer[Line], &len, Fin);



  /* The number of characters in the second line of the record should be equal to the number
  of characters to the last line of  the record.*/
  for(i = 1; i < max_lines ; i+=4){
	if (strlen(buffer[i+2])!=strlen(buffer[i]))
	{
		printf("ERROR Lines %d and %d have different length\n",i, i+2);
		exit(-1);
	}
  }


  /* Count the number of GC and non-GC nucleotides
  per read*/
  for( i = 1; i < max_lines; i+=4){
	GC = 0;
	AT = 0;
	for (int k=0;k< strlen(buffer[i])-1;k++)
	{
		switch (buffer[i][k])
		{
		  case 'G':
		  case 'C':GC++;
		  break;
		  case 'A':
		  case 'T': AT++;
		  break;
		}
	  }
	  fprintf(Fout,"%f\t%s", 1.0*GC/(GC+AT),buffer[i]);
	}


  /*close the files opened*/
  fclose(Fin);
  fclose(Fout);

  /* free the allocated memory */
  for (i=0;i<max_lines;i++)
  free(buffer[i]);
  free(buffer);

  exit(0);

}
