#include <stdio.h>
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include "modelSelector.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <chrono>
#include "utilities_MS.h"

using namespace std;
int main(int argc, char **argv) {
	
	
	if (argc < 2) {
		cout << "Input files are not specified" << endl;
		cout << "Default syntax is: markovModelSelector <tab> ";
		cout << "undirectedEdgeListFileName <tab> ";
		cout << "sequenceFileName <tab> ";
		cout << "u_name <tab> ";
		cout << "v_name" << endl;
	} else {		
		string undirectedEdgeListFileName = argv[1];
		string sequenceFileName = argv[2];
		string u_name = argv[3];
		string v_name = argv[4];
		struct stat buffer;
		if (stat (sequenceFileName.c_str(), &buffer) != 0) {
			cout << "Please check if the input filenames are correct" << endl;			
		} else {						
//			cout << "printing numeric limits " << endl;
//			cout << "min value for double is " << numeric_limits<double>::min()<< endl;
//			cout << "log min value for double is " << log(numeric_limits<double>::min())<< endl;
			ModelSelector * MS = new ModelSelector(sequenceFileName, undirectedEdgeListFileName, u_name, v_name);
			MS->PerformModelSelection();
			MS->WriteTotalComputeTime();
			delete MS;
		}
	}
	return 0;
}
