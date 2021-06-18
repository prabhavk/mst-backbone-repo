#ifndef EM_H
#define EM_H

#include <random>
#include <iostream>
#include <iomanip>
#include <chrono>

class EM_vertex{
public:
	string name;
	vector <EM_vertex*> neighbors;
	int timesVisited = 0;
	int sizeOfSubtree = 1;
	int degree = 0;
	int id;	
	void AddNeighbor(EM_vertex* v_ptr);
	
	EM_vertex(int idToAdd){
		id = idToAdd;
	}
	~EM_vertex(){		
	}
};

void EM_vertex::AddNeighbor(EM_vertex* v_ptr){
	degree += 1;
	neighbors.push_back(v_ptr);
};

class EM_tree{
private:
	default_random_engine generator;	
	double loglikelihood;
	double maximumLoglikelihood = 0.0;
	float rootProbability [4];	
	int numberOfSites;
	int numberOfSitePatterns;
	vector <vector<int>>* sitePatternRepeats;
	vector <int>* sitePatternWeights;
	vector <pair<int,int>>* edgesOrderedFromRootToLeaves;
	vector <array<double,4>*>* conditionalLikelihoodVector;
	vector <array<array<float,4>,4>*>* transitionMatrices;	
// Initialization
	void ComputeMaximumParsimonyEstimateOfAncestralSequences();
//	Expectation
	void ComputeConditionalLogLikelihoods();
	void ComputeMaximumAPosterioriEstimateOfAncestralSequences();
// Maximization
	void ComputeMaximumLikelihoodEstimateOfTransitionMatrices();
	void ComputeMaximumLikelihoodEstimateOfRootProbability();
//  In a tree with n leaves
//  sequences at indices 0 to n-1 corresponds to leaves 
//  sequences at indices n to 2n-3 corresponds to ancestral sequences updated with most recent EM run 	
//  Ancestral sequences for vertex rooting that maximizes likelihood.

public:
	int numberOfLeaves;
	int numberOfVerticesInMST;
	vector <vector <unsigned char>>* maximumLikelihoodAncestralSequences;
	vector <unsigned char>* rootSequence;
 	vector <vector<unsigned char>>* sequences;
	vector <EM_vertex*> vertices;
	vector <int> indsOfVerticesOfInterest;
	vector <pair<int,int>> edgesOfInterest;
	vector <int> indsOfVerticesToKeepInMST;
	int numberOfVerticesInSubtree;	
	void SetNumberOfVerticesInSubtree(int numberOfVerticesInSubtreeToAdd){
		numberOfVerticesInSubtree = numberOfVerticesInSubtreeToAdd;
	}
// Sequences are stored in the same order as the vertices in the unrooted tree
	void AddSequences(vector<vector<unsigned char>>* sequencesToAdd) {
		array<double,4>* conditionalLikelihood;
		sequences = sequencesToAdd;
		numberOfLeaves = sequences->size();
		for (int v_ind = 0; v_ind < 2*numberOfLeaves-2; v_ind++){
			conditionalLikelihood = new array<double, 4>;
			conditionalLikelihoodVector->push_back(conditionalLikelihood);
		}
		for (int edgeIndex = 0; edgeIndex < 2*numberOfLeaves-2; edgeIndex++){
			array<array<float,4>,4>* transitionMatrix = new array<array<float,4>,4>;
			transitionMatrices->push_back(transitionMatrix);
		}
		
	}
	
	void AddSitePatternWeights(vector <int>* sitePatternWeightsToAdd){
		sitePatternWeights = sitePatternWeightsToAdd;
		numberOfSitePatterns = sitePatternWeights->size();
		numberOfSites = 0;
		for (int sitePatternWeight: *sitePatternWeights){
			numberOfSites += sitePatternWeight;
		}
	}
	void AddSitePatternRepeats(vector<vector <int>>* sitePatternRepeatsToAdd) {sitePatternRepeats = sitePatternRepeatsToAdd;}
	void ComputeNJTree();
	void SelectIndsOfVerticesOfInterestAndEdgesOfInterest();	
	void PerformEM();	
// Iterates over edge e in t and
// (i) sets edges ordered from root to leaves
// (ii) performs EM for current edge set
// (iii) updates non-root ancestral sequences if likelihood score of current EM run is highest
	void SetEdgesOrderedFromRootToLeaves(pair<int,int> edge);
// Root ind is -1
// (i) performs EM for current edge set
// (ii) updates ancestral sequences if current EM run has highest likelihood score
	void PerformEMForCurrentEdgeSet();
	void PerformStructuralEM();
	void ComputeMutualInformation();
	// For each vertex pair compute mutual information 
	void ComputeChowLiuTree();
	// Compute maximum-weight spanning tree using mutual information
	void SelectVertexForRooting();
	// Select vertex for rooting
	void ContractHiddenVerticesWithInDegreeOne();
	// Contract hidden vertices with in-degree one
	void AddHiddenVertices();
	// Add hidden vertices such that there are n-1 hidden vertices	
	EM_tree(int numberOfVerticesInMSTToAdd){
		numberOfVerticesInMST = numberOfVerticesInMSTToAdd;
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		generator = default_random_engine(seed);
		edgesOrderedFromRootToLeaves = new vector<pair<int,int>>;
		rootSequence = new vector<unsigned char>;
		transitionMatrices = new vector <array<array<float,4>,4>*>;
		conditionalLikelihoodVector = new vector <array<double,4>*>;
		maximumLikelihoodAncestralSequences = new vector<vector<unsigned char>>;
	}
	~EM_tree(){
		for (EM_vertex* v : vertices){
			delete v;
		}
		sequences->clear();
		edgesOrderedFromRootToLeaves->clear();
		for (array<double,4>* conditionalLikelihood : *conditionalLikelihoodVector){
			delete conditionalLikelihood;
		}
		conditionalLikelihoodVector->clear();
		for (array<array<float,4>,4>* transitionMatrix : *transitionMatrices){
			delete transitionMatrix;
		}
		transitionMatrices->clear();
		for (array<double,4>* conditionalLikelihood : *conditionalLikelihoodVector){
			delete conditionalLikelihood;
		}
		
		maximumLikelihoodAncestralSequences->clear();
		delete maximumLikelihoodAncestralSequences;
		delete edgesOrderedFromRootToLeaves;
		delete rootSequence;
		delete transitionMatrices;
		delete conditionalLikelihoodVector;
	}
};

//void EM_tree::PerformStructuralEM(){
//	// Given a rooted tree and estimated of ancestral states
//	EM_tree:ComputeMutualInformation();
//	EM_tree:ComputeChowLiuTree();
//}
	

void EM_tree::ComputeNJTree(){
	map <pair<int,int>,float> distanceMap;
	map <int,float> R;
	vector <int> vertexIndsForIterating;
	unsigned int n = sequences->size();
	for (unsigned int i = 0; i < n; i++){
		R[i] = 0.0;
		vertexIndsForIterating.push_back(i);
		EM_vertex* v_ptr = new EM_vertex (i);
		vertices.push_back(v_ptr);
	}
	float distance;
	vector <unsigned char> seq_i; vector <unsigned char> seq_j;
	for (unsigned int i = 0; i < n; i++){		
		for (unsigned int j = i+1; j < n; j++){
			distance = 0.0;			
			seq_i = (*sequences)[i];
			seq_j = (*sequences)[j];
			for (int site=0; site < numberOfSitePatterns; site++){
				if (seq_i[site]!=seq_j[site]){
					distance += float((*sitePatternWeights)[site]);					
				}
			}			
			distanceMap[pair<int,int>(i,j)]=distance;
			R[i] += distance;
			R[j] += distance;
		}
	}
	for (unsigned int i = 0; i < n; i++){R[i] /= (n-2);}
	float neighborDist;
	float neighborDist_current;
	float newDist;
	int i; int j; int k;	
	int i_selected; int j_selected;
	int h = n;
	
	while (R.size()>3){
		neighborDist = 10000.0;
        for (unsigned int i_ind = 0; i_ind < n; i_ind++){
			for (unsigned int j_ind = i_ind+1; j_ind < n; j_ind++){
				i = vertexIndsForIterating[i_ind];
				j = vertexIndsForIterating[j_ind];				
				neighborDist_current = distanceMap[pair<int,int>(i,j)]-R[i]-R[j];
				if (neighborDist_current < neighborDist){
					neighborDist = neighborDist_current;
                    i_selected = i;
                    j_selected = j;
				}  
			}
		}
		vertexIndsForIterating.erase(remove(vertexIndsForIterating.begin(),vertexIndsForIterating.end(),i_selected),vertexIndsForIterating.end());
		vertexIndsForIterating.erase(remove(vertexIndsForIterating.begin(),vertexIndsForIterating.end(),j_selected),vertexIndsForIterating.end());		
		vertexIndsForIterating.push_back(h);
		EM_vertex* vertex_h_ptr = new EM_vertex (h);
		vertices.push_back(vertex_h_ptr);
		vertices[i_selected]->AddNeighbor(vertices[h]);
		vertices[j_selected]->AddNeighbor(vertices[h]);
		vertices[h]->AddNeighbor(vertices[i_selected]);
		vertices[h]->AddNeighbor(vertices[j_selected]);
        R.erase(i_selected);
        R.erase(j_selected);
        R[h]=0.0;
		n-=1;		
		for (unsigned int k_ind = 0; k_ind < n-1; k_ind++){
			k = vertexIndsForIterating[k_ind];
			if (k < i_selected){
				newDist = distanceMap[pair<int,int>(k,i_selected)] + distanceMap[pair<int,int>(k,j_selected)];
                newDist -= distanceMap[pair<int,int>(i_selected,j_selected)];
                newDist*=0.5;
                R[k] = float(R[k]*(n-1)-distanceMap[pair<int,int>(k,i_selected)]-distanceMap[pair<int,int>(k,j_selected)] + newDist)/float(n-2);
                distanceMap.erase(pair<int,int>(k,i_selected));
				distanceMap.erase(pair<int,int>(k,j_selected));               
			} else if (j_selected < k){
				newDist = distanceMap[pair<int,int>(i_selected,k)] + distanceMap[pair<int,int>(j_selected,k)];
                newDist -= distanceMap[pair<int,int>(i_selected,j_selected)];
                newDist*=0.5;
                R[k] = float(R[k]*(n-1)-distanceMap[pair<int,int>(i_selected,k)]-distanceMap[pair<int,int>(j_selected,k)] + newDist)/float(n-2);
                distanceMap.erase(pair<int,int>(i_selected,k));
                distanceMap.erase(pair<int,int>(j_selected,k));
			} else {
			    newDist = distanceMap[pair<int,int>(i_selected,k)] + distanceMap[pair<int,int>(k,j_selected)];
                newDist -= distanceMap[pair<int,int>(i_selected,j_selected)];
                newDist*=0.5;
                R[k] = float(R[k]*(n-1)-distanceMap[pair<int,int>(i_selected,k)]-distanceMap[pair<int,int>(k,j_selected)] + newDist)/float(n-2);
                distanceMap.erase(pair<int,int>(i_selected,k));
                distanceMap.erase(pair<int,int>(k,j_selected));
			}            
			distanceMap[pair<int,int>(k,h)] = newDist;
            R[h]+=newDist;
		}
        R[h]/=float(n-2);
		h+=1;
        distanceMap.erase(pair<int,int>(i_selected,j_selected));
	}
	EM_vertex* vertex_h_ptr = new EM_vertex (h);
	vertices.push_back(vertex_h_ptr);
	for (int v:vertexIndsForIterating){
		vertices[v]->AddNeighbor(vertices[h]);
		vertices[h]->AddNeighbor(vertices[v]);
	}	
};

void EM_tree::SelectIndsOfVerticesOfInterestAndEdgesOfInterest(){
	vector <EM_vertex*> verticesToVisit;
	vector <EM_vertex*> verticesVisited;
	int n_ind; EM_vertex* v;
	pair <int,int> edgeToAdd;
	for (int v_ind = 0; v_ind < numberOfVerticesInSubtree; v_ind++){
		verticesToVisit.push_back(vertices[v_ind]);
	}
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);		
		numberOfVerticesToVisit -= 1;
		for (EM_vertex* n: v->neighbors){
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()){
				n->timesVisited += 1;
				if ((n->degree - n->timesVisited) == 1){
					verticesToVisit.push_back(n);
					numberOfVerticesToVisit +=1;
				}
			} else {
				edgesOfInterest.push_back(pair<int,int>(v->id,n->id));
				n_ind = find(verticesVisited.begin(),verticesVisited.end(),n) - verticesVisited.begin();
				verticesVisited.erase(verticesVisited.begin()+n_ind);
			}
		}
	}	
	for (EM_vertex* v: verticesVisited){		
		if (v->id > numberOfLeaves -1){
			indsOfVerticesOfInterest.push_back(v->id);
		} else {
			indsOfVerticesToKeepInMST.push_back(v->id);
		}
	}
}

void EM_tree::SetEdgesOrderedFromRootToLeaves(pair<int,int> edge){
	int u_ind; int v_ind; int n_ind;
	EM_vertex* v;
	vector <EM_vertex*> verticesVisited;
	vector <EM_vertex*> verticesToVisit;
	verticesVisited.push_back(vertices[edge.first]);
	verticesVisited.push_back(vertices[edge.second]);	
	verticesToVisit.push_back(vertices[edge.first]);
	verticesToVisit.push_back(vertices[edge.second]);
	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
	u_ind = find(vertices.begin(), vertices.end(), vertices[edge.first]) - vertices.begin();
	v_ind = find(vertices.begin(), vertices.end(), vertices[edge.second]) - vertices.begin();	
	edgesOrderedFromRootToLeaves->clear();
	edgesOrderedFromRootToLeaves->push_back(pair<int,int>(-1,u_ind));
	edgesOrderedFromRootToLeaves->push_back(pair<int,int>(-1,v_ind));	
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		v_ind = find(vertices.begin(), vertices.end(), v) - vertices.begin();	
		for (EM_vertex* n: v->neighbors){
			if (find(verticesVisited.begin(), verticesVisited.end(), n) == verticesVisited.end()){				
				n_ind = find(vertices.begin(), vertices.end(), n) - vertices.begin();				
				edgesOrderedFromRootToLeaves->push_back(pair<int,int>(v_ind,n_ind));
				verticesVisited.push_back(vertices[n_ind]);
				verticesToVisit.push_back(vertices[n_ind]);
				numberOfVerticesToVisit += 1;
			}			
		}		
	}
}

void EM_tree::PerformEMForCurrentEdgeSet(){
	//	Initialize sequences using maximum parsimony
	ComputeMaximumParsimonyEstimateOfAncestralSequences();
	//	Iterate the following until loglikelihood value stops increasing
	int numberOfIterations = 0;
	double currentLogLikelihood = 0;
	this->loglikelihood = 0;
	bool continueEM = 1;
	while (continueEM and numberOfIterations < 10){
		numberOfIterations += 1;
		currentLogLikelihood = this->loglikelihood;
		ComputeMaximumLikelihoodEstimateOfRootProbability();
//		cout << "Computed maximum likelihood estimate of root probability";
		//	Estimate root probability
		ComputeMaximumLikelihoodEstimateOfTransitionMatrices();
//		cout << "Computed maximum likelihood estimate of transition matrices";
		//	Estimate transition matrices
		ComputeMaximumAPosterioriEstimateOfAncestralSequences();
//		cout << "Computed maximum a priori estimate of ancestral sequences";
		//	Compute expected states
		if (numberOfIterations < 2 or currentLogLikelihood < this->loglikelihood or abs(currentLogLikelihood - this->loglikelihood) > 0.001){
			continueEM = 1;
		} else {
			continueEM = 0;
		}
	}
	if (maximumLoglikelihood == 0 or maximumLoglikelihood < this->loglikelihood){
		maximumLoglikelihood = this->loglikelihood;
//		cout << maximumLoglikelihood << endl;
	// update ancestral sequences
		for (int v_ind = numberOfLeaves; v_ind < 2*numberOfLeaves-2; v_ind++){
			(*maximumLikelihoodAncestralSequences)[v_ind-numberOfLeaves] = (*sequences)[v_ind];
		}
	}
}


void EM_tree::PerformEM(){
	vector <int> aList;
	vector <int> bList;
	vector <pair<int,int>> edgesVisited;
	pair <int,int> edge;
	bool dontStop = 1;
	//	Initialize ancestral sequences
	for (int v_ind = numberOfLeaves; v_ind < 2*numberOfLeaves-2; v_ind++) {
		vector <unsigned char> sequenceToAdd;
		sequences->push_back(sequenceToAdd);		
	}
	// Initialize maximum likelihood ancestral sequences
	for (int v_ind = numberOfLeaves; v_ind < 2*numberOfLeaves-2; v_ind++) {
		vector <unsigned char> sequenceToAdd;
		maximumLikelihoodAncestralSequences->push_back(sequenceToAdd);		
	}
	for (EM_vertex* v: vertices){
		for (EM_vertex* n: v->neighbors){
			if (v->id < n->id){
				edge = pair<int,int>(v->id,n->id);
			} else {
				edge = pair<int,int>(n->id,v->id);
			}
			if (find(edgesVisited.begin(), edgesVisited.end(), edge) == edgesVisited.end() and dontStop){
				edgesVisited.push_back(edge);
				SetEdgesOrderedFromRootToLeaves(edge);
				PerformEMForCurrentEdgeSet();
			}
		}
	}
	edgesVisited.clear();
}

void EM_tree::ComputeMaximumParsimonyEstimateOfAncestralSequences(){
	int p_ind; int c_ind;	
	vector <unsigned char> V;
	unsigned char V_root;
	unsigned char stateAtParent;
	vector <vector<unsigned char>> VU;
	vector <unsigned char> VU_root;
	map <unsigned char, unsigned char> dnaCount;
	unsigned char maxCount = 0; unsigned char numberOfPossibleStates;
//	Initialize sequences for ancestors
	rootSequence->clear();
	for (int v_ind = numberOfLeaves; v_ind < 2*numberOfLeaves-2; v_ind++){
		(*sequences)[v_ind].clear();
	}
//	Compute V and VU for leaves
	for (int site = 0; site < numberOfSitePatterns; site++){
		V.clear();
		VU.clear();
		for (int v_ind = 0; v_ind < numberOfLeaves; v_ind++){
			unsigned char dna_base = ((*sequences)[v_ind])[site];
			V.push_back(dna_base);
			vector <unsigned char> vectorToAdd;
			vectorToAdd.push_back(dna_base);
			VU.push_back(vectorToAdd);
		}
//	Set values of V and VU for ancestral vertices
 		for (int v_ind = numberOfLeaves; v_ind < 2*numberOfLeaves-2; v_ind++){
			V.push_back(char(0));
			vector <unsigned char> vectorToAdd;			
			VU.push_back(vectorToAdd);
		}
// Given a root r and vertices u and v that are adjacent to r
// Iterate edges from leaves to u, v and compute VU for all ancestors
 		for (int edgeIndex = 2*numberOfLeaves-3; edgeIndex > 1; edgeIndex -= 2){
			maxCount = 0;
			dnaCount[0]=0; dnaCount[1]=0; dnaCount[2]=0; dnaCount[3]=0;
			tie (p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex];		
			for (unsigned char s: VU[c_ind]){dnaCount[s] +=1;}						
			tie (p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex-1];
			for (unsigned char s: VU[c_ind]){dnaCount[s] +=1;}			
			for (pair<unsigned char,unsigned char> dnaCountMember: dnaCount){
				if (maxCount < dnaCountMember.second){
					maxCount = dnaCountMember.second;}}
			VU[p_ind].clear();
			for (pair<unsigned char,unsigned char> dnaCountMember: dnaCount){
				if (dnaCountMember.second == maxCount){
					VU[p_ind].push_back(dnaCountMember.first);
				}				
			}			
		}		
		maxCount = 0;
		dnaCount[0]=0; dnaCount[1]=0; dnaCount[2]=0; dnaCount[3]=0;
		for (int edgeIndex = 1; edgeIndex > -1; edgeIndex--){
			tie (p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex];
			for (unsigned char s: VU[c_ind]){dnaCount[s] +=1;}}
		for (pair<unsigned char,unsigned char> dnaCountMember: dnaCount){
			if (maxCount < dnaCountMember.second){
				maxCount = dnaCountMember.second;}}
//	Compute VU for r
		VU_root.clear();
		for (pair<unsigned char,unsigned char> dnaCountMember: dnaCount){
			if (dnaCountMember.second == maxCount){
				VU_root.push_back(dnaCountMember.first);}}
//	Compute V for r
		if (VU_root.size()>1){
			numberOfPossibleStates = VU_root.size();
			uniform_int_distribution<int> distribution(0,numberOfPossibleStates-1);
			int pos = distribution(generator);
			V_root = VU_root[pos];			
		} else {
			V_root = VU_root[0];
		}
//	Select final values of V for all non-leaf vertices
		for (pair<int,int> edge : *edgesOrderedFromRootToLeaves){
			tie (p_ind, c_ind) = edge;
			// if V[p_ind] \in VU[c_ind] set V[c_ind] to V[p_ind]
			// else select V[c_ind] at random from VU[c_ind]
			if (p_ind == -1){
				stateAtParent = V_root;
			} else {
				stateAtParent = V[p_ind];
			}
			if(c_ind > numberOfLeaves-1){
				if (find(VU[c_ind].begin(),VU[c_ind].end(),stateAtParent)!=VU[c_ind].end()){
					V[c_ind] = stateAtParent;
				} else {
					numberOfPossibleStates = VU[c_ind].size();
					if (numberOfPossibleStates >1){
						uniform_int_distribution<int> distribution(0,numberOfPossibleStates-1);
						int pos = distribution(generator);
						V[c_ind] = VU[c_ind][pos];
					} else {						
						V[c_ind] = VU[c_ind][0];
					}					
				}
			}
		}
//	Push value of root ancestral sequence
		rootSequence->push_back(V_root);
//	Push values of all non-root ancestral sequences
		for (int v_ind = numberOfLeaves; v_ind < 2*numberOfLeaves-2; v_ind++){
			(*sequences)[v_ind].push_back(V[v_ind]);
		}	
	}
}

void EM_tree::ComputeMaximumLikelihoodEstimateOfRootProbability(){
	unsigned char dna_ind;
	map <unsigned char,int> dnaCount;
	dnaCount[0]=0; dnaCount[1]=0; dnaCount[2]=0; dnaCount[3]=0;
	for (int site =0; site < numberOfSitePatterns ; site++){
		dna_ind = (*rootSequence)[site];
		dnaCount[dna_ind] += (*sitePatternWeights)[site];
	}
	for (unsigned char dna_ind = 0; dna_ind < 4; dna_ind++){		
		rootProbability[dna_ind] = float(dnaCount[dna_ind])/float(numberOfSites);		
	}
}

void EM_tree::ComputeMaximumLikelihoodEstimateOfTransitionMatrices(){
	int p_ind; int c_ind; int site; float rowSum;	
	unsigned char dna_ind_p; unsigned char dna_ind_c;
	vector <unsigned char> seq_p; vector <unsigned char> seq_c; 	
	for (int edgeIndex = 0; edgeIndex < 2*numberOfLeaves-2; edgeIndex++){
		tie (p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex];		
		if (p_ind == -1){
			seq_p = *rootSequence;
		} else {
			seq_p = (*sequences)[p_ind];
		}
		seq_c = (*sequences)[c_ind];				
		for (dna_ind_p = 0; dna_ind_p < 4; dna_ind_p++){
			for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c++){				
				(*(*transitionMatrices)[edgeIndex])[dna_ind_p][dna_ind_c] = float(0);
			}
		}
		for (site = 0; site < numberOfSitePatterns;  site++){
			(*(*transitionMatrices)[edgeIndex])[seq_p[site]][seq_c[site]] += (*sitePatternWeights)[site];
		}
		for (dna_ind_p = 0; dna_ind_p < 4; dna_ind_p ++){
			rowSum = 0.0;
			for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
				rowSum += (*(*transitionMatrices)[edgeIndex])[dna_ind_p][dna_ind_c];
			}
			for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
				(*(*transitionMatrices)[edgeIndex])[dna_ind_p][dna_ind_c] /= rowSum;
			}
		}
	}
}

void EM_tree::ComputeMaximumAPosterioriEstimateOfAncestralSequences(){
	int p_ind; int c_ind; unsigned char dna_ind_p; unsigned char dna_ind_c;
	long double partialLikelihood; long double maxProbability; unsigned char stateWithMaxProbability;
	long double currentProbability; long double siteLikelihood;
	array <long double,4>* conditionalLikelihood_root = new array <long double,4>;
//		Iterate over site patterns
	loglikelihood = 0;
	for (int site = 0; site < numberOfSitePatterns; site++){
//	    Initialize conditional likelihoods
		for (int v_ind = 0; v_ind < numberOfLeaves; v_ind++){
			for (unsigned char dna_ind = 0; dna_ind < 4; dna_ind++){
				(*(*conditionalLikelihoodVector)[v_ind])[dna_ind] = 0.0;
			}
			(*(*conditionalLikelihoodVector)[v_ind])[(*sequences)[v_ind][site]] = 1.0;
		}
		for (int v_ind = numberOfLeaves; v_ind < 2*numberOfLeaves-2; v_ind++){
			for (unsigned char dna_ind = 0; dna_ind < 4; dna_ind++){
				(*(*conditionalLikelihoodVector)[v_ind])[dna_ind] = 1.0;
			}
		}
		for (unsigned char dna_ind = 0; dna_ind < 4; dna_ind++){
			(*conditionalLikelihood_root)[dna_ind] = 1.0;
		}
//    	Compute conditional likelihoods for all ancestors
		for (int edgeIndex = 2*numberOfLeaves-3; edgeIndex > 1; edgeIndex --){
			tie(p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex];
			for (unsigned char dna_ind_p = 0; dna_ind_p < 4; dna_ind_p++){
				partialLikelihood = 0.0;
				for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c++){
					partialLikelihood += (*(*transitionMatrices)[edgeIndex])[dna_ind_p][dna_ind_c] * (*(*conditionalLikelihoodVector)[c_ind])[dna_ind_c];
				}
				(*(*conditionalLikelihoodVector)[p_ind])[dna_ind_p] *= partialLikelihood;
			}
		}
		for (int edgeIndex = 1; edgeIndex > 0; edgeIndex --){
			tie(p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex];
			for (unsigned char dna_ind_p = 0; dna_ind_p < 4; dna_ind_p++){
				partialLikelihood = 0.0;
				for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c++){
					partialLikelihood += (*(*transitionMatrices)[edgeIndex])[dna_ind_p][dna_ind_c] * (*(*conditionalLikelihoodVector)[c_ind])[dna_ind_c];
				}
				(*conditionalLikelihood_root)[dna_ind_p] *= partialLikelihood;
			}
		}
		
//      Compute MAP estimate for root sequence
		maxProbability = -1; stateWithMaxProbability = 10;	
		for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c++) {
			currentProbability = rootProbability[dna_ind_c];
			currentProbability *= (*conditionalLikelihood_root)[dna_ind_c];
			if (currentProbability > maxProbability) {
				maxProbability = currentProbability;
				stateWithMaxProbability = dna_ind_c;
			}
		}
//		cout << "conditional likelihood is" << (*conditionalLikelihood_root)[0] << endl;
//		cout << "root probability is " << rootProbability[0] << endl;
		if (stateWithMaxProbability > 3) {
//			cout << maxProbability << "\tError in computing maximum a posterior estimate for ancestor vertex\n";
		} else {
			(*rootSequence)[site] = stateWithMaxProbability;
		}
//		Compute MAP estimate for each ancestral sequence
		for (int edgeIndex = 0; edgeIndex < 2; edgeIndex ++) {
			tie(p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex];
			if (c_ind > numberOfLeaves-1) {
				maxProbability = -1; stateWithMaxProbability = 10;
				dna_ind_p = (*rootSequence)[site];
				for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){ 
					currentProbability = (*(*transitionMatrices)[edgeIndex])[dna_ind_p][dna_ind_c];
					currentProbability *= (*(*conditionalLikelihoodVector)[c_ind])[dna_ind_c];
					if (currentProbability > maxProbability) {
						maxProbability = currentProbability;
						stateWithMaxProbability = dna_ind_c;
					}
				}
				if (stateWithMaxProbability > 3) {
//					cout << "Error in computing maximum a posterior estimate for ancestor vertex";
				} else {
					(*sequences)[c_ind][site] = stateWithMaxProbability;
				}
			}
		}
		for (int edgeIndex = 2; edgeIndex < 2*numberOfLeaves-2; edgeIndex ++) {
			tie(p_ind, c_ind) = (*edgesOrderedFromRootToLeaves)[edgeIndex];	
			if (c_ind > numberOfLeaves-1) {
				maxProbability = -1; stateWithMaxProbability = 10;
				dna_ind_p = (*sequences)[p_ind][site];
				for (dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++) {
					currentProbability = (*(*transitionMatrices)[edgeIndex])[dna_ind_p][dna_ind_c];
					currentProbability *= (*(*conditionalLikelihoodVector)[c_ind])[dna_ind_c];
					if (currentProbability > maxProbability) {
						maxProbability = currentProbability;
						stateWithMaxProbability = dna_ind_c;
					}
				}
				if (stateWithMaxProbability > 3) {
//					cout << "Error in computing maximum a posteriori estimate for ancestor vertex";
				} else {
					(*sequences)[c_ind][site] = stateWithMaxProbability;
				}
			}
		}
		siteLikelihood = 0.0;
		for (unsigned char dna_ind_p = 0; dna_ind_p < 4; dna_ind_p++){
			siteLikelihood += rootProbability[dna_ind_p] * (*conditionalLikelihood_root)[dna_ind_p];
		}
		try {
			loglikelihood += log(siteLikelihood) * (*sitePatternWeights)[site];
		} catch (const exception& e) {
			cout << "numerical computation error: error in computing log likelihood";
		}	
	}
	delete conditionalLikelihood_root;
}
#endif