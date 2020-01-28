#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

typedef struct {
	double GC;
	char * sequence;
}pair;

int main(int argc,char **argv)
{

	if ( argc != 4 ){
		printf("%s\n", "./fastq_parallel input_file output_file num_threads");
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
	
	int threads = atoi(argv[3]);
	if( threads <= 0 ){
		printf("Invalid number of threads!\n");
		exit(-1);
	}
	

  int Line;
  int max_lines = 0;
  char ch;
  size_t index = 0;
  
  /*Count the lines in the file*/
  while ((ch = fgetc(Fin)) != EOF)
  {
    if (ch == '\n')
		max_lines++;
  }
  
  /*Restart file pointer*/
  fseek(Fin, 0, SEEK_SET);
  
  /* malloc a 2d array with max_lines (lines of the file) as the rows*/
  char ** buffer;
  buffer=(char**)malloc(sizeof(char*)*max_lines );

  pair *pairs;
  /*Each record consists of 4 lines, so the total number of records will be max_lines / 4*/
  size_t pairs_size = max_lines / 4;
  /*So, in order to group each line in one struct variable, an array of pairs will be created*/
  pairs = (pair *) malloc(sizeof(pair)*pairs_size);

  size_t len = 0;

  /* Read line-by-line the file and store each in the array named buffer */
  for(Line=0;Line<max_lines;Line++)
	getline(&buffer[Line], &len, Fin);

   int i;
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

  int k = 0;
  int GC=0;
  int AT=0;
  #pragma omp parallel for num_threads(threads) default(none) private(i, AT, GC) shared(pairs, index, max_lines, buffer) schedule(dynamic)
  for( i = 1; i < max_lines; i+=4){
	AT = 0;
	GC = 0;
	#pragma omp parallel for num_threads(2) default(none) private(k) shared(buffer, i) reduction(+:AT,GC)
	for (int k=0;k< strlen(buffer[i])-1;k++ )
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
	  #pragma omp critical
	  {
		pairs[index].GC =  1.0*GC/(GC+AT);
		pairs[index].sequence = buffer[i];
		index++;
	  }
	}



  /*output the nucleotide sequence and the GC content fastq to the output file*/
  for( i = 0; i < pairs_size; i++)
	fprintf(Fout,"%f\t%s", pairs[i].GC,pairs[i].sequence);

  /*close the files opened*/
  fclose(Fin);
  fclose(Fout);

  /* free the allocated memory */
  for (i=0;i<max_lines;i++)
  free(buffer[i]);
  free(buffer);

  exit(0);

}

