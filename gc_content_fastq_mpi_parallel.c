#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


typedef struct data_pair{
	char   *dna_sequence;
	char   *phred_sequence;
	int dna_line;
	int phred_line;
} data_pair;

typedef data_pair* Data_pair;
/*parameters:
	Fin: a file descriptor.
	lines: the memory address where the value of the number of lines of the file will be stored.
	max_column_length: the memory address where the value of number of characters contained in the longest line in file will be stored.
	A function to count the lines of the file and the number of characters contained in the line with the maximum width.*/
void countLinesAndMaxColumnLength(FILE *Fin, int *lines, int *max_column_length);
/*parameters:
	rows: number of rows
	cols: number of columns
A function that allocates memory for rows*cols contiguous characters*/
char **alloc_2d_char(int rows, int cols);

int main(int argc, char **argv){
    
    if ( argc != 3 ){
		printf("%s input_file output_file\n", argv[0]);
		MPI_Finalize();
		exit(-1);
	}

    FILE * Fin= fopen(argv[1], "r");
	if ( !Fin ){
		printf("File %s could not be opened!\n", argv[1]);
		exit(-1);
	}
	
	/*Number of lines read from the file*/
	int lines;
	/*Number of characters contained in the line with the maximum width in file.*/
	int max_column_length;
    countLinesAndMaxColumnLength(Fin, &lines, &max_column_length);
	
	/*Restart file pointer*/
	fseek(Fin, 0, SEEK_SET);
	/*Allocate memory for 2d array of characters*/
    char **buffer=(char**)malloc(sizeof(char*)*lines);
	for(int i = 0; i < lines ; i ++){
		buffer[i] = (char *)malloc(sizeof(char)*max_column_length);
	}
	
	size_t len = 0;
	/*Read line-by-line the file*/
	for(int Line=0;Line < lines;Line++)
	    getline(&buffer[Line], &len, Fin); 
	
	/*Calculate the number of records( one record has four lines)*/
	int chunk_size = lines  / 4;
	/*Allocate memory for each record*/
	Data_pair data_pairs = (Data_pair)malloc(sizeof(data_pair)*chunk_size);
	for(int i = 0 ; i < chunk_size ; i++){
		data_pairs[i].dna_sequence = (char *) malloc(sizeof(char)*max_column_length);
		data_pairs[i].phred_sequence = (char *) malloc(sizeof(char)*max_column_length);
	}
	/*Initialize each record with the values stored in the buffer*/
	for(int i = 0, j = 1 ; i < chunk_size ; i++, j+=4){
		strcpy(data_pairs[i].dna_sequence, buffer[j]);
		strcpy(data_pairs[i].phred_sequence, buffer[j+2]);
		data_pairs[i].dna_line = j+1;
		data_pairs[i].phred_line = j+2+1;
	}

	/*Free buffer*/
	for(int i = 0 ; i < lines ; i++)
		free(buffer[i]);
	free(buffer);


	int rank;
	int size;
	

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	int local_n = chunk_size / size;
	int START = rank * local_n;
	int END = START + local_n;
	if(rank == size - 1)
		END = chunk_size;
	
	/* The number of characters in the second line of the record should be equal to the number
	of characters to the last line of  the record.*/
	for(int i = START; i < END ; i++){
		if (strlen(data_pairs[i].dna_sequence) !=strlen (data_pairs[i].phred_sequence))
		{
			printf("ERROR Lines %d and %d have different length\n", data_pairs[i].dna_line, data_pairs[i].phred_line);
			exit(-1);
		}
	}


	int GC = 0;
	int AT = 0;
	int offset = END - START;
	int message_length = max_column_length+15;
	char ** send_buffer = alloc_2d_char(offset, message_length);
	int j = 0;
	/* Count the number of GC and non-GC nucleotides
	per read*/
	for(int i = START; i < END ; i++){
		GC = 0;
		AT = 0;
		for (int k=0;k< strlen(data_pairs[i].dna_sequence);k++)
		{

			switch (data_pairs[i].dna_sequence[k])
			{
			case 'G':
			case 'C':GC++;
			break;
			case 'A':
			case 'T': AT++;
			break;
			}
		}
		/*Convert value of gc (double) to string.*/
		char gcd[9];
	  	sprintf(gcd,"%.6lf", 1.0*GC/(GC+AT));
		gcd[8] = '\0';
		/*Create message format and store it to send_buffer*/
		sprintf(send_buffer[j], "%s\t%s", gcd, data_pairs[i].dna_sequence);
		j++;
	}
	
	if( rank > 0){
		/*Send the number of messages stored in send_buffer to process with rank 0*/
		MPI_Send(&offset, 1, MPI_INT, 0, 2 , MPI_COMM_WORLD);
		/*Send the whole array to process with rank 0*/
		MPI_Send(&(send_buffer[0][0]), offset*message_length, MPI_CHAR, 0 , 1, MPI_COMM_WORLD);
		/*Release the allocated memory*/
		free(send_buffer[0]);
		free(send_buffer);
		free(data_pairs);
	}

	if(rank == 0){
		FILE * Fout= fopen(argv[2], "w");
		if(!Fout){
			printf("File %s could not be opened!\n", argv[2]);
			exit(-1);
		}
		/*Write to file the content of send buffer that process with rank 0 has stored*/
		for(int i = 0 ; i < offset ; i++){
			fprintf(Fout ,"%s", send_buffer[i]);
		}
		free(send_buffer[0]);
		free(send_buffer);

		MPI_Status status;
		MPI_Status status2;
		int count;
		char ** tmp_dna_sequences;

		int *offsets = (int *)malloc(sizeof(int)*size);
		/*Get the number of messages that each process stores to send_buffer*/
		for(int i = 1, offset_local = 0 ; i < size ; i++){
			MPI_Recv(&offset_local, 1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
			offsets[status.MPI_SOURCE] = offset_local; 
		}
		/*Get the 2d array from each process and write the content to the file*/
		for(int i = 1 ; i < size ; i++){
			tmp_dna_sequences = alloc_2d_char(offsets[i], message_length);
			MPI_Recv(&(tmp_dna_sequences[0][0]), offsets[i]*message_length, MPI_CHAR, i, 1, MPI_COMM_WORLD, &status2);
			for(int p = 0 ; p < offsets[i]; p++)
				fprintf(Fout,"%s", tmp_dna_sequences[p]);
			free(tmp_dna_sequences[0]);
			free(tmp_dna_sequences);
		}
		free(offsets);

		fclose(Fin);
  		fclose(Fout);
	}

	MPI_Finalize();
    
    return 0;
}

void countLinesAndMaxColumnLength(FILE *Fin, int *lines, int *max_column_length){
	*lines = 0;
	int column_length = 0;
	*max_column_length = 0;
	char ch;
	while ((ch = fgetc(Fin)) != EOF)
	{
		column_length++;
	    if (ch == '\n'){
			(*lines)++;
			if( *max_column_length < column_length)
				*max_column_length = column_length;
			column_length = 0;
		}
	}
}

char **alloc_2d_char(int rows, int cols) {
    char *data = (char *)malloc(rows*cols*sizeof(char));
    char **array= (char **)malloc(rows*sizeof(char*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}