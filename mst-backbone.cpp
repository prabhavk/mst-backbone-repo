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
		// duplicate_sequences
		// highly similar sequnces (Hamming distance of 1; try upto 5)
		// incremental construction		
		// genotype-recombination network
		// Select subset of vertices based on phenotypes
		// Amino acid/Drug resistance/Phenotype network
		// Network layout
		// Dynamic network layout?		
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
	string distance_measure_for_NJ = "Hamming";
	bool flag_verbose = 0;
	MSTBackbone * MSTBackboneObj;
    if (argc < 2) {        
        cerr << "Usage: " << argv[0] << " --seq alignment.fas --constraint_size size_of_subtree --distance_measure_for_NJ logDet --out prefix_for_output --root_supertree no" << endl;
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
        // prefix_to_output_files (e.g. prefix is "covid-19_size_42" and the alignment is in directory "/foo/bar/data/")		
            } else if (strcmp(argv[i], "--distance_measure_for_NJ") == 0) {
                if (i < argc -1) {
					distance_measure_for_NJ = argv[++i];
					if (distance_measure_for_NJ == "logDet" || distance_measure_for_NJ == "Jukes-Cantor" ||distance_measure_for_NJ == "Hamming") {
						continue;
					} else {
						cout << "Enter one of the following distance measures: logDet, Jukes-Cantor, Hamming" << endl;
						exit (-1);
					}
                }
        // prefix_to_output_files (e.g. prefix is "covid-19_size_42" and the alignment is in directory "/foo/bar/data/")		
            } else if (strcmp(argv[i], "--out") == 0) {
                if (i < argc -1) {
					flag_prefix = 1;
                    prefix_for_output_files = argv[++i];
					// prefix_for_output_files = alignment_file_path_obj.parent_path().string() + prefix_for_output_files;
					prefix_path_obj =  alignment_file_path_obj.parent_path();
					prefix_path_obj /= prefix_for_output_files;
					cout << "prefix for output files is " << prefix_for_output_files << endl;					
                }
        // size of subtree Vs
            } else if (strcmp(argv[i], "--verbose") == 0) {
                if (i < argc -1) {
			string_verbose = argv[++i];
			if (string_verbose == "True") {
				flag_verbose = 1;
			}			
                }        
            } else if (strcmp(argv[i], "--constraint_size") == 0) {
                if (i < argc -1) {
					flag_size_of_subtree = 1;
                    size_of_subtree = stoi(argv[++i]);
                }
			} else if (strcmp(argv[i], "--help") == 0) {
				cout << "Example for mst-constrained tree construction: " << argv[0] << " --seq alignment.fas --constraint_size size --distance_measure_for_NJ logDet --out prefix_for_output --root_supertree no" << endl;
			} else if (strcmp(argv[i], "--root_super_tree") == 0) {
				if (i < argc -1) {					
					root_supertree = argv[++i];
					if (root_supertree == "yes" || root_supertree == "no") {
						continue;
					} else {
						cout << "Enter one of the following options for root_supertree: no, yes" << endl;
						exit (-1);
					}					
                }
			}
        }

		if (!flag_size_of_subtree) {
			size_of_subtree = 10;
		}

		if (!flag_prefix) {
			prefix_path_obj =  alignment_file_path_obj.parent_path();
			prefix_path_obj /= "mstbackbone_output";
			// prefix_for_output_files = alignment_file_path_obj.parent_path().string() + "_mstbackbone";
		}		
		MSTBackboneObj = new MSTBackbone(path_to_alignment_file, size_of_subtree, prefix_path_obj.string(),patch_name,distance_measure_for_NJ,flag_verbose,root_supertree);		
		delete MSTBackboneObj;
		// MSTBackbone MSTBackboneObj(path_to_alignment_file, size_of_subtree, prefix_path_obj.string(),localPhyloOnly,modelSelection,modelForRooting,useChowLiu);
    }

	return 0;
} 
