#include<iostream>
#include<stdio.h>
#include<ctype.h>
#include<stdlib.h>
#include<getopt.h>
#include "include/simulation.h"

using namespace std;

double t=1., j=1., H=.0;
int n=20, d=50, print_int=1,  opt;
const char* filename=NULL;

void parse_input_args(int argc, char **argv){
	// cout<<"parsing input args";
	while ((opt = getopt(argc,argv,"t:j:H:n:d:f:p:")) != EOF)
	{
		switch(opt)
		{
			case 't':
				t = strtod(optarg, NULL);
				break;
			case 'j':
				j = strtod(optarg, NULL);
				break;
			case 'd':
				d = stoi (optarg);
				break;
			case 'H':
				H = strtod(optarg, NULL);
				break;
			case 'f':
				filename = optarg;
				break;
			case 'n':
				n = stoi(optarg);
				break;
			case 'p':
				print_int = stoi(optarg);
				break;
			case '?':
				fprintf(stderr, "usuage is \n\
	-t : for setting temperature, \n\
	-j : for setting interaction energy, \n\
	-H : for setting external field,\n\
	-d : for setting number of advances in simulation,\n\
	-f : for setting output filename. otherwise will print to stdout\n\
	-p : for setting the frequency of printing\n\
	-n : side of square lattice");
				break;
			default:
				cout<<endl;return;
		}
	}
}


int main(int argc, char **argv)
{
	FILE* fp = stdout;
	parse_input_args(argc, argv);
	if (filename != NULL){
		fp = fopen(filename, "w");
	}
	simulation s = simulation (n, t, j, H);
	s.print_interval=print_int;
	s.advance(d, fp);
	fclose(fp);
	return 0;
}
