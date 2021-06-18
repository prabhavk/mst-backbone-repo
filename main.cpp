#include <stdio.h>
#include <string>
#include "mstBackbone.h"
#include <iostream>
#include <sys/stat.h>
using namespace std;

int main(int argc, char **argv)
{	
	if (argc < 2) {
		cout << "Input file is not specified" << endl;		
		cout << "Default syntax is: mst-backbone <tab> sequencesFileName" << endl;
	} else {	
		string sequenceFileName = argv[1];		
		struct stat buffer;
		if (stat (sequenceFileName.c_str(), &buffer) != 0){
			cout << "Please check if the input filename is correct" << endl;
		} else {			
			int numberOfLargeEdgesInSubtree;
			if (argc > 2) {
				numberOfLargeEdgesInSubtree = stoi(argv[2]);
			} else {
				numberOfLargeEdgesInSubtree = 10;
			}	

			if (numberOfLargeEdgesInSubtree < 10){
				numberOfLargeEdgesInSubtree = 10;
			}			
			// Fasta format (done)
			// Phylip format
			// Ambiguous characters
			// Gaps
			// logDet distance = LogDet(s1,s2)
			// e = smallest distance that is greater than zero 
			// n = total number of vertices
			// perturbed distance d = LogDet(s1,s2) + (e/(2*n-3))*(min(s1.id+1,s2.id+1)) + (e/n^2)*(max(s1.id,s2.id))
			MSTBackbone MSTBackboneObj(sequenceFileName, numberOfLargeEdgesInSubtree);	
		}
	}
	return 0;
} 