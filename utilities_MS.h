#ifndef UTILITIES_H
#define UTILITIES_H
#include <vector>
#include <tuple>
#include <map>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

//#include "asa047.hpp"
using namespace Eigen;
using namespace std;

Matrix4f ConstructRateMatrix(vector <float> parameters);
Matrix4f ConstructRateMatrix(vector <float> parameters){
	float p1 = parameters[0];
	float p2 = parameters[1];
	float p3 = parameters[2];
	float a = parameters[3];
	float b = parameters[4];
	float c = parameters[5];
	float d = parameters[6];
	float e = parameters[7];
	float f = parameters[8];
	float g = parameters[9];
	float h = parameters[10];
	float p4 = 1-(p1+p2+p3);
	float j = -1*(-p1*(a+b+c)+p2*d+p3*g)/p4;	
	float k = -1*(-p2*(d+e+f)+p1*a+p3*h)/p4;	
	float i = (1-(p1*(a+c)+p2*(d+f)+2*p3*(g+h)+p4*(j+k)))/(2*p3);
	float l = -1*(-p3*(g+h+i)+p1*b+p2*e)/p4;	
	Matrix4f Q;
	Q	<< -(a+b+c) , a , b , c,
		 d , -(d+e+f) , e , f,
		 g , h , -(g+h+i) , i,
		 j , k , l , -(j+k+l);
	return Q;
}

string EncodeAsDNA(vector<unsigned char> sequence);

string EncodeAsDNA(vector<unsigned char> sequence){
	string allDNA = "AGTC";
	string dnaSequence = "";
	for (unsigned char s : sequence){
		dnaSequence += allDNA[s];
	}
	return dnaSequence;
}
//
//tuple <vector<vector<unsigned char>>,vector<int>,vector<vector<int>>> MST_tree::GetCompressedSequencesSiteWeightsAndSiteRepeats(vector<int> vertexIdList){	
//	vector <vector<unsigned char>> compressedSequences;
//	vector <int> sitePatternWeights_ptr;
//	vector <vector <int>> sitePatternRepeats_ptr;
//	vector <vector<unsigned char>> distinctPatterns;
//	map <vector<unsigned char>,vector<int>> distinctPatternsToSitesWherePatternRepeats;
//	vector <MST_vertex*> vertexPtrList;
//	for (unsigned int i = 0; i < vertexIdList.size(); i++){		
//		MST_vertex* v_ptr = vertexMap[vertexIdList[i]];
//		vertexPtrList.push_back(v_ptr);
//		vector<unsigned char> compressedSequence;
//		compressedSequences.push_back(compressedSequence);
//	}
//	int numberOfSites = vertexPtrList[0]->sequence.size();
//	vector<unsigned char> sitePattern;
//	for(int site=0; site < numberOfSites; site++){
//		sitePattern.clear();
//		for (MST_vertex* v_ptr: vertexPtrList){
//			sitePattern.push_back(v_ptr->sequence[site]);}
//		if (find(distinctPatterns.begin(),distinctPatterns.end(),sitePattern)!=distinctPatterns.end()){
//			distinctPatternsToSitesWherePatternRepeats[sitePattern].push_back(site);
//			
//		} else {
//			distinctPatterns.push_back(sitePattern);	
//			vector<int> sitePatternRepeats;
//			sitePatternRepeats.push_back(site);
//			distinctPatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
//			for (unsigned int i=0; i < sitePattern.size();i++){				
//				compressedSequences[i].push_back(sitePattern[i]);
//			}
//		}
//	}
//	for (vector<unsigned char> sitePattern : distinctPatterns){
//		int sitePatternWeight = distinctPatternsToSitesWherePatternRepeats[sitePattern].size();
//		sitePatternWeights_ptr.push_back(sitePatternWeight);		
//		sitePatternRepeats_ptr.push_back(distinctPatternsToSitesWherePatternRepeats[sitePattern]);
//	}
//	return make_tuple(compressedSequences,sitePatternWeights_ptr,sitePatternRepeats_ptr);
//}

float ComputeNormalizedHammingDistance(vector<unsigned char> * recodedSeq1, vector<unsigned char> * recodedSeq2) {
	float hammingDistance = 0;
	int sequenceLength = recodedSeq1->size();
	for (int i=0;i<sequenceLength;i++){
		if ((*recodedSeq1)[i] != (*recodedSeq2)[i]){
			hammingDistance += 1;
		}		
	}	
	hammingDistance /= float(sequenceLength);	
	return (hammingDistance);
};

array <float, 4> GetBaseFreq(vector <unsigned char> * seq);

array <float, 4> GetBaseFreq(vector <unsigned char> * seq){
	array <float, 4> rootProbability;
	for (int i = 0; i < 4; i++){
		rootProbability[i] = 0;
	}
	int sequenceLength = seq->size();
	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
		rootProbability[(*seq)[sitePos]] += 1;
	}
	for (unsigned char dna = 0; dna < 4; dna ++){
		rootProbability[dna] /= float(sequenceLength);
	}
	return (rootProbability);
}

float ComputeDeterminantOfFourByFourMatrix(array<array<float,4>,4>* M);
float ComputeDeterminantOfFourByFourMatrix(array<array<float,4>,4>* M){
	float a = (*M)[0][0];
	float b = (*M)[0][1];
	float c = (*M)[0][2];
	float d = (*M)[0][3];
	float e = (*M)[1][0];
	float f = (*M)[1][1];
	float g = (*M)[1][2];
	float h = (*M)[1][3];
	float i = (*M)[2][0];
	float j = (*M)[2][1];
	float k = (*M)[2][2];
	float l = (*M)[2][3];
	float m = (*M)[3][0];
	float n = (*M)[3][1];
	float o = (*M)[3][2];
	float p = (*M)[3][3];
	float determinant;
	determinant = a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n - b*e*k*p + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m + c*e*j*p - c*e*l*n - c*f*i*p + c*f*l*m + c*h*i*n - c*h*j*m - d*e*j*o + d*e*k*n + d*f*i*o - d*f*k*m - d*g*i*n + d*g*j*m;
	return (determinant);
}

float ComputeLogDetDistance(vector<unsigned char> * seq1, vector <unsigned char> * seq2);

float ComputeLogDetDistance(vector<unsigned char> * seq1, vector <unsigned char> * seq2){
	array<array<float,4>,4> F;
	float logDet = 0.0;
	int seqLength = seq1->size();
	float d = 1/float(seqLength);
	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			F[i][j] = 0;
		}
	}
	for (int pos = 0 ; pos < seqLength; pos ++){
		F[(*seq1)[pos]][(*seq2)[pos]] += d;
	}
//	f1_a, f1_c, f1_g, f1_t = GetBaseFreq(seq1)
	array <float, 4> f1 = GetBaseFreq(seq1);
//   g1 = f1_a*f1_c*f1_g*f1_t
	float g1 = f1[0]*f1[1]*f1[2]*f1[3];
//  f2_a, f2_c, f2_g, f2_t = GetBaseFreq(seq2)
	array <float, 4> f2 = GetBaseFreq(seq2);    
//  g2 = f2_a*f2_c*f2_g*f2_t
	float g2 = f2[0]*f2[1]*f2[2]*f2[3]; 
//    ld = -0.25*(m.log(np.linalg.det(F))-0.5*(m.log(g1*g2)))
	float determinant = ComputeDeterminantOfFourByFourMatrix(&F);
	if (determinant > 0){
		logDet = -0.25*(log(determinant)-0.5*(log(g1*g2)));
	}	
	return (logDet);
}

vector<unsigned char> DecompressSequence(vector<unsigned char>* compressedSequence, vector<vector<int>>* sitePatternRepeats);

vector<unsigned char> DecompressSequence(vector<unsigned char>* compressedSequence, vector<vector<int>>* sitePatternRepeats){
	int totalSequenceLength = 0;
	for (vector<int> sitePatternRepeat: *sitePatternRepeats){
		totalSequenceLength += int(sitePatternRepeat.size());
	}
	vector <unsigned char> decompressedSequence;
	for (int v_ind = 0; v_ind < totalSequenceLength; v_ind++){
		decompressedSequence.push_back(char(0));
	}
	unsigned char dnaToAdd;
	for (int sitePatternIndex = 0; sitePatternIndex < int(compressedSequence->size()); sitePatternIndex++){
		dnaToAdd = (*compressedSequence)[sitePatternIndex];
		for (int pos: (*sitePatternRepeats)[sitePatternIndex]){
			decompressedSequence[pos] = dnaToAdd;
		}
	}
	return (decompressedSequence);	
}

void DecodeSequencesAndWriteToFasta(map<int,vector<unsigned char>> integerEncodedSequenceMap, string fileName){
	ofstream fastaFile;
	fastaFile.open(fileName);
	string characterSequence;
	for (pair<int,vector<unsigned char>> idAndSequence: integerEncodedSequenceMap){
		fastaFile << ">" << idAndSequence.first << endl;
		characterSequence = EncodeAsDNA(idAndSequence.second);
		fastaFile << characterSequence << endl;
	}
	fastaFile.close();
}

// void WriteSequencesToFasta(map<int,string> sequencesMap, string fileName);

// void WriteSequencesToFasta(map<int,string> sequencesMap, string fileName){
// 	ofstream fastaFile;	
// 	fastaFile.open(fileName);	
// 	for (pair<int,string> idAndSequence : sequencesMap){
// 		fastaFile << ">" << idAndSequence.first << "\n";
// 		fastaFile << idAndSequence.second << "\n";
// 	}
// 	fastaFile.close();
// };

// double ComputeLikelihoodOfObservingSequence(vector<unsigned char> * seq);

// double ComputeLikelihoodOfObservingSequence(vector<unsigned char> * seq){
// 	array <float, 4> rootProbability;
// 	double logLikelihoodForObservingSequence = 0;
// 	for (int i = 0; i < 4; i++){
// 		rootProbability[i] = 0;
// 	}
// 	int sequenceLength = seq->size();
// 	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
// 		rootProbability[(*seq)[sitePos]] += 1;
// 	}
// 	for (unsigned char dna = 0; dna < 4; dna ++){
// 		rootProbability[dna] /= float(sequenceLength);
// 	}
// 	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
// 		logLikelihoodForObservingSequence += log(rootProbability[(*seq)[sitePos]]);
// 	}
// 	return (logLikelihoodForObservingSequence);
// };

// double ComputeEdgeLogLikelihood(vector<unsigned char> * seq_p, vector<unsigned char> * seq_c);

// double ComputeEdgeLogLikelihood(vector<unsigned char> * seq_p, vector<unsigned char> * seq_c){
// 	double edgeLogLikelihood = 0;
// 	array <array<float,4>,4> transitionMatrix;
// 	int sequenceLength = seq_p->size();
// 	unsigned char char_p; unsigned char char_c;
// 	float rowSum;
// 	for (char_p = 0; char_p < 4; char_p++){
// 		for (char_c = 0; char_c < 4; char_c++){
// 			transitionMatrix[char_p][char_c] = 0.0;
// 		}
// 	}
// 	for (int sitePos =0; sitePos < sequenceLength; sitePos++){
// 		char_p = (*seq_p)[sitePos];
// 		char_c = (*seq_c)[sitePos];
// 		transitionMatrix[char_p][char_c] += 1;
// 	}
// 	for (char_p = 0; char_p < 4; char_p++){
// 		rowSum = 0;
// 		for (char_c = 0; char_c < 4; char_c++){
// 			rowSum += transitionMatrix[char_p][char_c];
// 		}
// 		for (char_c = 0; char_c < 4; char_c++){
// 			transitionMatrix[char_p][char_c] /= float(rowSum);
// 		}
// 	}
// 	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
// 		char_p = (*seq_p)[sitePos];
// 		char_c = (*seq_c)[sitePos];
// 		edgeLogLikelihood += log(transitionMatrix[char_p][char_c]);
// 	}
// 	return edgeLogLikelihood;
// };

// tuple <vector<vector<unsigned char>>,vector<int>> GetCompressedSequencesAndSiteWeights(vector<vector<unsigned char>> fullSequences){		
// 	vector <vector <unsigned char>> compressedSequences;
// 	vector <int> sitePatternWeights;
// 	vector <vector <int> > sitePatternRepeats;
// 	vector <vector <unsigned char>> distinctPatterns;
// 	int sitePatternWeight;
// 	map <vector <unsigned char>,vector<int>> distinctPatternsToSitesWherePatternRepeats;
// 	for (unsigned int v_ind = 0; v_ind < fullSequences.size(); v_ind++){
// 		vector <unsigned char> compressedSequence;
// 		compressedSequences.push_back(compressedSequence);
// 	}
// 	int numberOfSites = fullSequences[0].size();
// 	vector<unsigned char> sitePattern;
// 	for (int site = 0; site < numberOfSites; site++){
// 		sitePattern.clear();
// 		for (unsigned int v_ind = 0; v_ind < fullSequences.size(); v_ind++){
// 			sitePattern.push_back(fullSequences[v_ind][site]);
// 		}
// 		if (find(distinctPatterns.begin(),distinctPatterns.end(),sitePattern)!=distinctPatterns.end()){
// 			distinctPatternsToSitesWherePatternRepeats[sitePattern].push_back(site);
			
// 		} else {
// 			distinctPatterns.push_back(sitePattern);
// 			vector <int> sitePatternRepeats;
// 			sitePatternRepeats.push_back(site);
// 			distinctPatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
// 			for (unsigned int i = 0; i < sitePattern.size(); i++){
// 				compressedSequences[i].push_back(sitePattern[i]);
// 			}
// 		}
// 	}
// 	for (vector<unsigned char> sitePattern : distinctPatterns){
// 		sitePatternWeight = distinctPatternsToSitesWherePatternRepeats[sitePattern].size();
// 		sitePatternWeights.push_back(sitePatternWeight);		
// 		sitePatternRepeats.push_back(distinctPatternsToSitesWherePatternRepeats[sitePattern]);
// 	}
// 	return make_tuple(compressedSequences,sitePatternWeights);	
// }

#endif