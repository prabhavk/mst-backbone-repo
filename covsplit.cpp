#include <stdio.h>
#include <string>
#include "mstBackbone.h"
#include <iostream>
#include <experimental/filesystem>
// #include <boost/filesystem.hpp>
#include <sys/stat.h>
using namespace std;
namespace fs = std::experimental::filesystem;

int main(int argc, char **argv)
{
    string path_to_alignment_file;
	string patch_name;	
	fs::path alignment_file_path_obj;
	fs::path prefix_path_obj;
    string prefix_for_output_files;
	string path_to_directory;
    int size_of_subtree;
 	struct stat buffer;
	// bool alignment_file_not_set = 0;
	bool flag_size_of_subtree = 0;
	bool flag_prefix = 0;
	// bool localPhyloOnly = 0;
	bool flag_modelSelection = 0;
	// bool useChowLiu = 1;
	string modelForRooting = "UNREST";
	string root_supertree = "no";
	string arg_localPhyloOnly;
	string arg_modelSelection;
	string arg_chowliu;
	string string_verbose;
	string distance_measure_for_NJ = "logDet";
	bool flag_verbose = 0;
	MSTBackbone * MSTBackboneObj;
    if (argc < 2) {        
        cerr << "Usage: " << argv[0] << " --seq alignment.fas --constraint_size size_of_subtree --distance_measure_for_NJ LogDet --out prefix_for_output --root_supertree no" << endl;
		cerr << endl;
        return (-1);
    } else {        
        // parse arguments            
        for (int i = 1; i < argc ; i++) {
        // path to multiple sequence alignment file
			if (strcmp(argv[i], "--patch") == 0) {		
                if (i < argc -1) {					
                    patch_name = argv[++i];
					cout << "Applying patch " << patch_name << endl;
                }        
            } else if (strcmp(argv[i], "--seq") == 0) {				
                if (i < argc -1) {					
                    path_to_alignment_file = argv[++i];
					cout << "Alignment file name is " << path_to_alignment_file << endl;
					if (stat (path_to_alignment_file.c_str(), &buffer) != 0) { // check if input file exists
						cout << "Please check if the input filename is correct" << endl;
						exit (-1);
					}
					alignment_file_path_obj = path_to_alignment_file;
                }
            }
        }
    }
}