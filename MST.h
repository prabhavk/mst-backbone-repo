#ifndef MST_H
#define MST_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <stdio.h>
#include "utilities.h"
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <chrono>
using namespace std;

enum dna {a,c,g,t,n};
enum aa {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
class mut {
public:
	int id;
	string dna_mut_id;
	string full_mut_id;
	vector <pair <int,dna>> dna_mut_list; // store pos and DNA substitution
	vector <pair <int,aa>> aa_mut_list; // store pos and AA substitution
	vector <pair<int,int>> gap_list; // count ambiguous and N as gap store start and end pos
};
class MST_vertex {
public:
	string name;
	int degree = 0;
	int numberOfLargeEdgesInSubtree = 0;
	int timesVisited = 0;
	int id;
	int idOfExternalVertex;
	int rank = 0;
	map <int,mut> mut_id_2_mut_obj;
	vector <unsigned char> sequence;
	vector <unsigned char> globallyCompressedSequence;
	vector <int> idsOfVerticesInSubtree;
	vector <MST_vertex *> neighbors;
	void AddNeighbor(MST_vertex * v_ptr);
	void Get_mut_id(vector <pair<int,dna>> mut_list);
	MST_vertex(int idToAdd, string nameToAdd, vector<unsigned char> sequenceToAdd) {
		id = idToAdd;
		sequence = sequenceToAdd;
		name = nameToAdd;
		idsOfVerticesInSubtree.push_back(idToAdd);
		idOfExternalVertex = -1;
	}
	~MST_vertex() {
		
	}
};

void MST_vertex::AddNeighbor(MST_vertex* v_ptr) {
	degree += 1;
	neighbors.push_back(v_ptr);
}

class MST_tree {
private:	
	int largestVertexIndex;
	int edgeWeightThreshold = 0;
	vector <MST_vertex*> verticesToVisit;
	bool ContainsVertex(int vertex_id);
	map <pair<int,int>,int> * allEdgeWeights;
	chrono::system_clock::time_point current_time;
	chrono::system_clock::time_point start_time;
	chrono::system_clock::time_point time_to_compute_MST;
	bool build_MST_incrementally = false;
public:
	int maxDegree;
	vector <int> siteWeights;
	string sequenceFileName;
	int v_ind;
	int numberOfLargeEdgesThreshold_input = 0;	
	int numberOfLargeEdgesThreshold = 0;
	int numberOfNonZeroWeightEdges = 0;
	vector <int> idsOfExternalVertices;
	map <string, unsigned char> mapDNAtoInteger;
	MST_vertex * subtree_v_ptr;
	map <int, MST_vertex *> * vertexMap;
	map <pair <int, int>, int> edgeWeightsMap;
	void AddEdgeWeight(int u_id, int v_id, int edgeWeight);
	void RemoveEdgeWeight(int u_id, int v_id);
	void AddEdgeWeightToDistanceGraph(int u_id, int v_id, int edgeWeight);
	void RemoveWeightedEdgesIncidentToVertexInDistanceGraph(int u_id);
	void SetCompressedSequencesAndSiteWeightsForInputSequences();
	void SetNumberOfLargeEdgesThreshold(int numberOfLargeEdges);
	void SetEdgeWeightThreshold(int edgeWeight){edgeWeightThreshold = edgeWeight;}
	void AddVertex(string name, vector <unsigned char> sequence);
	void AddVertexWithId(int id, string name, vector <unsigned char> sequence);
	void RemoveVertex(int vertex_id);
	void AddEdge(int u_id, int v_id, int edgeWeight);
	void RemoveEdge(int u_id, int v_id);
	void ResetVertexAttributesForSubtreeSearch();
	void UpdateMSTWithMultipleExternalVertices(vector <int> idsOfVerticesToKeep, vector <int> idsOfVerticesToRemove, vector <tuple<int,string,vector<unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd, vector <int> idsOfExternalVertices);
    void UpdateMSTWithMultipleExternalVertices_simple(vector <int> idsOfVerticesToKeep, vector <int> idsOfVerticesToRemove, vector <tuple<int,string,vector<unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd, vector <int> idsOfExternalVertices);	
	void UpdateMaxDegree();
	void UpdateMSTWithOneExternalVertex(vector <int> idsOfVerticesToRemove, string nameOfSequenceToAdd, vector <unsigned char> sequenceToAdd);
	bool ContainsEdge(int u_id, int v_id);
	int GetEdgeWeight(int u_id, int v_id);	
	int GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices);
	int GetNumberOfVertices();
	void ReadSequences(string sequenceFileNameToSet);
	void ComputeMST();
	void ComputeChowLiuTree();
	void ComputeMST_nonACGT();
	void ResetSubtreeSizeThreshold();
	void DoubleSubtreeSizeThreshold();
	int ComputeHammingDistance(vector <unsigned char> recodedSeq1, vector <unsigned char> recodedSeq2);			
	pair <vector <int>, vector <int>> GetIdsForSubtreeVerticesAndExternalVertices();	
	pair <bool, MST_vertex *> GetPtrToVertexSubtendingSubtree();
	pair <vector <int>,vector <int>> GetSubtreeVerticesAndExternalVertices();
	tuple <vector <string>, vector <vector <unsigned char>>, vector <int>, vector <vector<int>>> GetCompressedSequencesSiteWeightsAndSiteRepeats(vector <int> vertexIdList);
	vector <int> GetIdsOfClosestUnvisitedVertices(MST_vertex* u_ptr);
	void SetIdsOfExternalVertices();
	bool ShouldIComputeALocalPhylogeneticTree();
	void WriteToFile(string fileName);
	unsigned char ConvertDNAToChar(char dna);
	MST_tree() {		
		this->v_ind = 0;
		vector <unsigned char> emptySequence;
		this->allEdgeWeights = new map <pair<int,int>,int> ; 
		this->vertexMap = new map <int, MST_vertex *>;
		this->mapDNAtoInteger["A"] = 0;
		this->mapDNAtoInteger["C"] = 1;
		this->mapDNAtoInteger["G"] = 2;
		this->mapDNAtoInteger["T"] = 3;		
		this->mapDNAtoInteger["-"] = 4;
		this->mapDNAtoInteger["N"] = 4;
		this->mapDNAtoInteger["W"] = 4;
		this->mapDNAtoInteger["S"] = 4;
		this->mapDNAtoInteger["M"] = 4;
		this->mapDNAtoInteger["K"] = 4;
		this->mapDNAtoInteger["R"] = 4;
		this->mapDNAtoInteger["Y"] = 4;
		this->mapDNAtoInteger["B"] = 4;
		this->mapDNAtoInteger["D"] = 4;
		this->mapDNAtoInteger["H"] = 4;
		this->mapDNAtoInteger["V"] = 4;		
		//mapDNAtoInteger["N"] = 4;
		//mapDNAtoInteger["N"] = 4;
		//mapDNAtoInteger["N"] = 4;
		//mapDNAtoInteger["N"] = 4;
		//mapDNAtoInteger["N"] = 4;
		//mapDNAtoInteger["N"] = 4;		
	}
	~MST_tree() {		
		for (pair<int,MST_vertex*> VptrMap: *this->vertexMap){			
			delete VptrMap.second;
		}
		delete this->vertexMap;
		delete this->allEdgeWeights;
	}
};

void MST_tree::AddEdgeWeightToDistanceGraph(int u_id, int v_id, int edgeWeight) {
	if (u_id < v_id) {
		this->allEdgeWeights->insert(make_pair(make_pair(u_id,v_id),edgeWeight));
	} else {
		this->allEdgeWeights->insert(make_pair(make_pair(v_id,u_id),edgeWeight));
	}
}

void MST_tree::RemoveWeightedEdgesIncidentToVertexInDistanceGraph(int u_id) {
	
}

void MST_tree::SetIdsOfExternalVertices() {
	this->idsOfExternalVertices.clear();
	this->idsOfExternalVertices = this->GetIdsOfClosestUnvisitedVertices(this->subtree_v_ptr);
}

int MST_tree::GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices) {
	int edgeIndex;
	edgeIndex = numberOfVertices * (numberOfVertices-1)/2;
	edgeIndex -= (numberOfVertices - vertexIndex1) * (numberOfVertices-vertexIndex1-1)/2;
	edgeIndex += vertexIndex2 - vertexIndex1 - 1;
	return edgeIndex;
}

void MST_tree::SetNumberOfLargeEdgesThreshold(int numberOfLargeEdges_toSet) {
	this->numberOfLargeEdgesThreshold_input = numberOfLargeEdges_toSet;
	this->numberOfLargeEdgesThreshold = numberOfLargeEdges_toSet;
}

bool MST_tree::ShouldIComputeALocalPhylogeneticTree() {
	bool valueToReturn;
	bool verbose = 0;
	bool subtreeExtractionPossible;
	int numberOfNonZeroWeightEdgesInVmWithoutVs;
	tie (subtreeExtractionPossible, this->subtree_v_ptr) = this->GetPtrToVertexSubtendingSubtree();	
	if (subtreeExtractionPossible) {
		numberOfNonZeroWeightEdgesInVmWithoutVs = this->numberOfNonZeroWeightEdges - this->subtree_v_ptr->numberOfLargeEdgesInSubtree;
		if (numberOfNonZeroWeightEdgesInVmWithoutVs > this->numberOfLargeEdgesThreshold) {
			valueToReturn = 1;			
		} else {
			if (verbose) {
				cout << "Case 1: subtree extraction possible but number of external vertices is too small" << endl;
			}
			valueToReturn = 0;
			
		}
	} else {
		if (verbose) {
			cout << "Case 1: subtree extraction is not possible" << endl;
		}
		valueToReturn = 0;
	}
	return (valueToReturn);
}



void MST_tree::ResetSubtreeSizeThreshold() {
	this->numberOfLargeEdgesThreshold = this->numberOfLargeEdgesThreshold_input;
}

void MST_tree::DoubleSubtreeSizeThreshold() {
	this->numberOfLargeEdgesThreshold = this->numberOfLargeEdgesThreshold * 2;
}
int MST_tree::ComputeHammingDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2) {
	int hammingDistance = 0;
	for (unsigned int i=0;i<recodedSeq1.size();i++){
		if (recodedSeq1[i] != recodedSeq2[i]){
			hammingDistance+=1;
		}		
	}
	return (hammingDistance);
};

int MST_tree::GetNumberOfVertices() {
	return this->vertexMap->size();
};


bool MST_tree::ContainsVertex(int vertex_id) {
	return this->vertexMap->find(vertex_id)!=vertexMap->end();
}

void MST_tree::AddVertex(string name, vector <unsigned char> sequence) {
	MST_vertex * v = new MST_vertex(this->v_ind, name, sequence);
	this->vertexMap->insert(pair<int,MST_vertex*>(this->v_ind,v));	
	this->v_ind += 1;
}

void MST_tree::AddVertexWithId(int id, string name, vector <unsigned char> sequence) {
	MST_vertex * v = new MST_vertex(id, name, sequence);
	this->vertexMap->insert(pair<int,MST_vertex*>(id,v));
}

void MST_tree::RemoveVertex(int vertex_id) {
	MST_vertex* v = (*this->vertexMap)[vertex_id];
	for (MST_vertex* n: v->neighbors) {
		if (n->id < v->id) {
			this->edgeWeightsMap.erase(pair<int,int>(n->id,v->id));
		} else {
			this->edgeWeightsMap.erase(pair<int,int>(v->id,n->id));
		}
		n->neighbors.erase(remove(n->neighbors.begin(),n->neighbors.end(),v),n->neighbors.end());
		n->degree -= 1;
	}
	v->neighbors.clear();
	v->sequence.clear();
	v->idsOfVerticesInSubtree.clear();
	this->vertexMap->erase(vertex_id);
	delete v;
}

void MST_tree::AddEdge(int u_id, int v_id, int edgeWeight) {
	MST_vertex* u_ptr = (*this->vertexMap)[u_id];
	MST_vertex* v_ptr = (*this->vertexMap)[v_id];
	u_ptr->AddNeighbor(v_ptr);
	v_ptr->AddNeighbor(u_ptr);
	this->AddEdgeWeight(u_id,v_id,edgeWeight);
//	if (u_id < v_id) {
//		this->edgeWeightsMap[pair<int,int>(u_id,v_id)] = edgeWeight;
//	} else { 
//		this->edgeWeightsMap[pair<int,int>(v_id,u_id)] = edgeWeight;
//	}
	if (edgeWeight > 0) {
		this->numberOfNonZeroWeightEdges += 1;
	}
};

void MST_tree::RemoveEdge(int u_id, int v_id) {	
	if (u_id < v_id) {
		this->edgeWeightsMap.erase(pair<int,int>(u_id, v_id));
	} else {
		this->edgeWeightsMap.erase(pair<int,int>(v_id, u_id));
	}
	MST_vertex * u = (*this->vertexMap)[u_id];
	MST_vertex * v = (*this->vertexMap)[v_id];
	u->neighbors.erase(remove(u->neighbors.begin(),u->neighbors.end(),v),u->neighbors.end());
	u->degree -= 1;
	v->neighbors.erase(remove(v->neighbors.begin(),v->neighbors.end(),u),v->neighbors.end());
	v->degree -= 1;
}

bool MST_tree::ContainsEdge(int u_id, int v_id) {
	if (u_id < v_id) {
		return (this->edgeWeightsMap.find(pair<int,int>(u_id,v_id)) != this->edgeWeightsMap.end());
	} else {
		return (this->edgeWeightsMap.find(pair<int,int>(v_id,u_id)) != this->edgeWeightsMap.end());
	}	
}

int MST_tree::GetEdgeWeight(int u_id, int v_id) {
	if (u_id < v_id) {
		return this->edgeWeightsMap[pair<int,int>(u_id,v_id)];
	} else {
		return this->edgeWeightsMap[pair<int,int>(v_id,u_id)];
	}
}

void MST_tree::AddEdgeWeight(int u_id, int v_id, int edgeWeight) {
	pair<int,int> edge ;
	if (u_id < v_id){
		edge = make_pair(u_id,v_id);
	} else {
		edge = make_pair(v_id,u_id);
	}
	if (this->edgeWeightsMap.find(edge) != this->edgeWeightsMap.end()) {
		this->edgeWeightsMap[edge] = edgeWeight;
	} else {
		this->edgeWeightsMap.insert(make_pair(edge,edgeWeight));
	}	
}

void MST_tree::RemoveEdgeWeight(int u_id, int v_id) {
	pair <int, int> edge;
	if (u_id < v_id){
		edge = make_pair(u_id,v_id);
	} else {
		edge = make_pair(v_id,u_id);
	}
	this->edgeWeightsMap.erase(edge);
}

vector<int> MST_tree::GetIdsOfClosestUnvisitedVertices(MST_vertex* v_ptr) {
	int numberOfLargeEdgesEncountered = 0;
	vector <int> idsOfClosestUnvisitedVertices;
	vector <MST_vertex*> verticesInCurrentLevel;
	for (MST_vertex* n_ptr: v_ptr->neighbors) {		
		if (find(v_ptr->idsOfVerticesInSubtree.begin(),v_ptr->idsOfVerticesInSubtree.end(),n_ptr->id)==v_ptr->idsOfVerticesInSubtree.end()) {
			idsOfClosestUnvisitedVertices.push_back(n_ptr->id);
			if (this->GetEdgeWeight(v_ptr->id,n_ptr->id)  > edgeWeightThreshold) {
				numberOfLargeEdgesEncountered+=1;
			}
			if (numberOfLargeEdgesEncountered < numberOfLargeEdgesThreshold) {
				verticesInCurrentLevel.push_back(n_ptr);
			}
		}
	}
	vector <MST_vertex *> verticesInNextLevel;
	while (numberOfLargeEdgesEncountered < numberOfLargeEdgesThreshold and verticesInCurrentLevel.size() > 0) {
		for(MST_vertex * x_ptr:verticesInCurrentLevel) {
			for (MST_vertex * n_ptr : x_ptr->neighbors) {
				if (find(idsOfClosestUnvisitedVertices.begin(),idsOfClosestUnvisitedVertices.end(),n_ptr->id)==idsOfClosestUnvisitedVertices.end() and n_ptr->id!=v_ptr->id) {
					idsOfClosestUnvisitedVertices.push_back(n_ptr->id);
					if(this->GetEdgeWeight(x_ptr->id,n_ptr->id) > edgeWeightThreshold) {
						numberOfLargeEdgesEncountered+=1;
					}
					if (numberOfLargeEdgesEncountered < numberOfLargeEdgesThreshold) {
						verticesInNextLevel.push_back(n_ptr);	
					}
				}
			}
		}
		verticesInCurrentLevel = verticesInNextLevel;
		verticesInNextLevel.clear();
	}
	return idsOfClosestUnvisitedVertices;
}

pair <bool, MST_vertex *> MST_tree::GetPtrToVertexSubtendingSubtree() {
	this->ResetVertexAttributesForSubtreeSearch();	
	bool subTreeFound = 0;
	verticesToVisit.clear();
	for (pair <int, MST_vertex *> VptrMap: * this->vertexMap) {
		if (VptrMap.second->degree == 1) {
			verticesToVisit.push_back(VptrMap.second);
		}
	}
	vector <MST_vertex *> verticesVisited;
	int vertex_ind = verticesToVisit.size() -1;	
	while(vertex_ind > -1 and !subTreeFound) {
		this->subtree_v_ptr = verticesToVisit[vertex_ind];
		verticesToVisit.pop_back();
		vertex_ind -= 1;
		this->subtree_v_ptr->timesVisited += 1;
		for (MST_vertex * neighbor_ptr : this->subtree_v_ptr->neighbors) {
			if (neighbor_ptr->timesVisited < neighbor_ptr->degree) {
				neighbor_ptr->timesVisited += 1;
				for (int n_id : this->subtree_v_ptr->idsOfVerticesInSubtree) {
					neighbor_ptr->idsOfVerticesInSubtree.push_back(n_id);
				}
				neighbor_ptr->numberOfLargeEdgesInSubtree += this->subtree_v_ptr->numberOfLargeEdgesInSubtree;
				if (GetEdgeWeight(this->subtree_v_ptr->id,neighbor_ptr->id) > edgeWeightThreshold) {
					neighbor_ptr->numberOfLargeEdgesInSubtree+=1;
				}
				if (neighbor_ptr->degree - neighbor_ptr->timesVisited == 1) {
					if (neighbor_ptr->numberOfLargeEdgesInSubtree > numberOfLargeEdgesThreshold) {
						subTreeFound = 1;
						// set id to external vertex
						for (MST_vertex * v : neighbor_ptr->neighbors) {
							if (v->timesVisited < v->degree) {
								neighbor_ptr->idOfExternalVertex = v->id;
							}							
						}						
						return pair <bool, MST_vertex *> (subTreeFound,neighbor_ptr);						
					}
					verticesToVisit.push_back(neighbor_ptr);
					vertex_ind+=1;
				}
			}
		}
	}	
	return pair <bool, MST_vertex *> (subTreeFound,this->subtree_v_ptr);
}

void MST_tree::ResetVertexAttributesForSubtreeSearch() {
	for (pair <int, MST_vertex *> VIdAndPtr: *this->vertexMap) {
		VIdAndPtr.second->numberOfLargeEdgesInSubtree = 0;
		VIdAndPtr.second->idsOfVerticesInSubtree.clear();
		VIdAndPtr.second->idsOfVerticesInSubtree.push_back(VIdAndPtr.second->id);
		VIdAndPtr.second->timesVisited = 0;
	}
}


void MST_tree::UpdateMSTWithOneExternalVertex(vector<int> idsOfVerticesToRemove, string nameOfSequenceToAdd, vector <unsigned char> sequenceToAdd) {
	//	Remove vertices		
	for (int v_id: idsOfVerticesToRemove) {	
		this->RemoveVertex(v_id);
	}
	//	Remove neighbors and reset vertex attributes
	for (pair<int,MST_vertex*> VIdAndPtr: *this->vertexMap) {
		VIdAndPtr.second->numberOfLargeEdgesInSubtree = 0;
		VIdAndPtr.second->idsOfVerticesInSubtree.clear();
		VIdAndPtr.second->idsOfVerticesInSubtree.push_back(VIdAndPtr.second->id);
		VIdAndPtr.second->timesVisited = 0;
		VIdAndPtr.second->neighbors.clear();
		VIdAndPtr.second->degree = 0;
	}
	this->numberOfNonZeroWeightEdges = 0;
	// Add vertex
	int indexOfVertexToAdd = this->v_ind;
	this->AddVertex(nameOfSequenceToAdd, sequenceToAdd);
	MST_vertex * v_add; MST_vertex * v_inMST;
	v_add = (*this->vertexMap)[indexOfVertexToAdd];
	assert(v_add->name == nameOfSequenceToAdd);
	int edgeWeight;
	for (pair <int, MST_vertex *> idPtrPair : *this->vertexMap) {
		if (idPtrPair.first != indexOfVertexToAdd) {
			v_inMST = idPtrPair.second;
			edgeWeight = ComputeHammingDistance(v_add->sequence, v_inMST->sequence);
			if (v_inMST->id < v_add->id){					
				this->edgeWeightsMap[pair<int,int>(v_inMST->id,v_add->id)] = edgeWeight;					
			} else {					
				this->edgeWeightsMap[pair<int,int>(v_add->id,v_inMST->id)] = edgeWeight;					
			}
		}
	}	
	int numberOfVertices = int(this->vertexMap->size());	
	const int numberOfEdges = int(this->edgeWeightsMap.size());
	
	int * weights;
	weights = new int [numberOfEdges];
	
	typedef pair <int,int > E;
	E * edges;
	edges = new E [numberOfEdges];
	
	int edgeIndex = 0;
    for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
		edges[edgeIndex] = E(edgeAndWeight.first.first,edgeAndWeight.first.second);
		weights[edgeIndex] = edgeAndWeight.second;
		edgeIndex += 1;		
	}
	
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_distance_t, int>, boost::property < boost::edge_weight_t, int> > Graph;
	Graph g(edges, edges + numberOfEdges, weights, numberOfVertices);
	vector < boost::graph_traits < Graph >::vertex_descriptor >  p(num_vertices(g));
	prim_minimum_spanning_tree(g, &p[0]);
	delete[] edges;		
	delete[] weights;		
	vector <pair<int,int>> edgeWeightsToKeep;	
	vector <pair<int,int>> edgeWeightsToRemove;	
	for (size_t u = 0; u != p.size(); ++u) {
		if (p[u] != u) {		
			this->AddEdge(u,p[u],this->GetEdgeWeight(u,p[u]));
			if (u < p[u]) {
				edgeWeightsToKeep.push_back(pair<int,int>(u,p[u]));
			} else {
				edgeWeightsToKeep.push_back(pair<int,int>(p[u],u));
			}
		}
	}
	for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
		if (find(edgeWeightsToKeep.begin(),edgeWeightsToKeep.end(),edgeAndWeight.first) == edgeWeightsToKeep.end()){
			edgeWeightsToRemove.push_back(edgeAndWeight.first);
		}
	}
	for (pair<int,int> edge: edgeWeightsToRemove) {
		this->edgeWeightsMap.erase(edge);
	}	
}

void MST_tree::UpdateMaxDegree() {
	this->maxDegree = 0;
	for (pair <int,MST_vertex*> VIdAndPtr: *this->vertexMap) {
		if (this->maxDegree	< VIdAndPtr.second->degree) {
			this->maxDegree	= VIdAndPtr.second->degree;
		}
	}
}

void MST_tree::UpdateMSTWithMultipleExternalVertices_simple(vector<int> idsOfVerticesToKeep, vector<int> idsOfVerticesToRemove, vector<tuple<int,string,vector<unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd, vector<int> idsOfExternalVertices) {	
	// Remove weights of edges incident to vertices to removess
	// Add weights of edges 
	// Remove vertices
	
}

void MST_tree::UpdateMSTWithMultipleExternalVertices(vector<int> idsOfVerticesToKeep, vector<int> idsOfVerticesToRemove, vector<tuple<int,string,vector<unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd, vector<int> idsOfExternalVertices) {
	// Remove weights of edges incident to vertex
//	MST_vertex * v;
//	for (int v_id: idsOfVerticesToRemove) {
//		v = (*this->vertexMap)[v_id];
//		for (MST_vertex * n: v->neighbors) {
//			this->RemoveEdgeWeight(n->id, v_id);
//		}	
//	}	
	//	Remove vertices		
	for (int v_id: idsOfVerticesToRemove) {
		this->RemoveVertex(v_id);
	}
//	int numOfEdgesInCurrentMST = (int) this->edgeWeightsMap.size();
	//	Remove all edges in MST and reset all attributes for each vertex
	for (pair <int,MST_vertex*> VIdAndPtr: *this->vertexMap) {
		VIdAndPtr.second->numberOfLargeEdgesInSubtree = 0;
		VIdAndPtr.second->idsOfVerticesInSubtree.clear();
		VIdAndPtr.second->idsOfVerticesInSubtree.push_back(VIdAndPtr.second->id);
		VIdAndPtr.second->timesVisited = 0;
		VIdAndPtr.second->neighbors.clear();
		VIdAndPtr.second->degree = 0;
	}
	this->numberOfNonZeroWeightEdges = 0;
	int u_id; int v_id; int edgeWeight;
	vector <unsigned char> seq_u; vector <unsigned char> seq_v;
	string u_name; string v_name;
	
	int numberOfVerticesToKeep = int(idsOfVerticesToKeep.size());
	int numberOfVerticesToAdd = int(idAndNameAndSeqTupleForVerticesToAdd.size());
	int numberOfExternalVertices = int(idsOfExternalVertices.size());
	
	if (numberOfVerticesToAdd > 1) {
		for (int u_ind = 0; u_ind < numberOfVerticesToAdd -1; u_ind++) {
			tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];
			for (int v_ind = u_ind + 1; v_ind < numberOfVerticesToAdd; v_ind++) {
				tie (v_id, v_name, seq_v) = idAndNameAndSeqTupleForVerticesToAdd[v_ind];
				edgeWeight = ComputeHammingDistance(seq_u, seq_v);
				this->AddEdgeWeight(u_id,v_id,edgeWeight);
//				if (u_id < v_id){
//					this->edgeWeightsMap[pair<int,int>(u_id,v_id)] = edgeWeight;
//				} else {
//					this->edgeWeightsMap[pair<int,int>(v_id,u_id)] = edgeWeight;
//				}
			}
		}
	}
// Add newly introduced vertices
	for (int u_ind = 0; u_ind < numberOfVerticesToAdd; u_ind++) {
		tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];
		this->AddVertexWithId(u_id, u_name, seq_u);
	}
// Add edge weights for vertices to add to external vertices
	for (int u_ind = 0; u_ind < numberOfVerticesToAdd; u_ind++) {
		tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];
		for (int v_ind = 0; v_ind < numberOfExternalVertices; v_ind++) {
			v_id = idsOfExternalVertices[v_ind];
			seq_v = (*this->vertexMap)[v_id]->sequence;
			edgeWeight = ComputeHammingDistance(seq_u,seq_v);
			this->AddEdgeWeight(u_id,v_id,edgeWeight);
//			if (u_id < v_id){
//				this->edgeWeightsMap[pair<int,int>(u_id,v_id)] = edgeWeight;
//			} else {
//				this->edgeWeightsMap[pair<int,int>(v_id,u_id)] = edgeWeight;
//			}
		}
	}
	
	// Add edge weights for vertices to add to vertices to keep
	for (int u_ind = 0; u_ind < numberOfVerticesToAdd; u_ind++) {
		tie (u_id, u_name, seq_u) = idAndNameAndSeqTupleForVerticesToAdd[u_ind];		
		for (int v_ind = 0; v_ind < numberOfVerticesToKeep; v_ind++) {
			v_id = idsOfVerticesToKeep[v_ind];
			seq_v = (*this->vertexMap)[v_id]->sequence;
			edgeWeight = ComputeHammingDistance(seq_u,seq_v);
			this->AddEdgeWeight(u_id,v_id,edgeWeight);
//			if (u_id < v_id){
//				this->edgeWeightsMap[pair<int,int>(u_id,v_id)] = edgeWeight;
//			} else {
//				this->edgeWeightsMap[pair<int,int>(v_id,u_id)] = edgeWeight;
//			}			
		}		
	}

	// Add edge weights for vertices to keep to external vertices
	for (int u_ind = 0; u_ind < numberOfVerticesToKeep; u_ind++) {
			u_id = idsOfVerticesToKeep[u_ind];
			seq_u = (*this->vertexMap)[u_id]->sequence;
		for (int v_ind = 0; v_ind < numberOfExternalVertices; v_ind++) {
			v_id = idsOfExternalVertices[v_ind];
			seq_v = (*this->vertexMap)[v_id]->sequence;
			edgeWeight = ComputeHammingDistance(seq_u,seq_v);
			this->AddEdgeWeight(u_id,v_id,edgeWeight);
//			if (u_id < v_id){
//				this->edgeWeightsMap[pair<int,int>(u_id,v_id)] = edgeWeight;
//			} else {
//				this->edgeWeightsMap[pair<int,int>(v_id,u_id)] = edgeWeight;
//			}
		}
	}
	
// number of edges should be 
	
	vector <int> mstIds;
//	map<pair<int,int>> primId2MSTId;
	map<int,int> mstId2PrimId;
	int primId = 0;
	int mstId;
	for (pair <int,MST_vertex*> idPtrPair : *this->vertexMap) {
		mstId = idPtrPair.first;
//		if (find(mstIds.begin(),mstIds.end(),mstId) == mstIds.end()) {
		mstIds.push_back(mstId);
		mstId2PrimId.insert(make_pair(mstId,primId));
		primId += 1;
//		}
	}
	int numberOfVertices = int(this->vertexMap->size());	
	assert(primId == numberOfVertices);
	const int numberOfEdges = int(this->edgeWeightsMap.size());
//	cout << "Number of vertices in distance graph is " << numberOfVertices << endl;
//	cout << "number of edges in distance graph is " << numberOfEdges << endl;
	int * weights;
	weights = new int [numberOfEdges];
	
	typedef pair <int,int > E;
	E * edges;
	edges = new E [numberOfEdges];
		
	
	int edgeIndex = 0;	
    for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
//		edges[edgeIndex] = E(edgeAndWeight.first.first,edgeAndWeight.first.second);
		tie (u_id, v_id) = edgeAndWeight.first;
		edges[edgeIndex] = E(mstId2PrimId[u_id],mstId2PrimId[v_id]);
		weights[edgeIndex] = edgeAndWeight.second;
		edgeIndex += 1;		
	}
	
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_distance_t, int>, boost::property < boost::edge_weight_t, int> > Graph;
	Graph g(edges, edges + numberOfEdges, weights, numberOfVertices);
	vector < boost::graph_traits < Graph >::vertex_descriptor >  p(num_vertices(g));
	
//	int start_vertex_for_prim = 0;
//	for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
//		start_vertex_for_prim = edgeAndWeight.first.second;
//		break;
//	}	
//	cout << "start vertex for prim is \t" << start_vertex_for_prim << endl;
	prim_minimum_spanning_tree(g, &p[0]);
		
	vector <pair<int,int>> edgeWeightsToKeep;	
	vector <pair<int,int>> edgeWeightsToRemove;	
//	if (this->maxDegree == 0) {
//	cout << "edges in updated MST are" << endl;				
//	}
	for (size_t u = 0; u != p.size(); ++u) {
//		if (this->maxDegree == 0) {
//			cout << p[u] << "\t" << u << endl;				
//		}		
		if (p[u] != u){
			u_id = mstIds[p[u]];
			v_id = mstIds[u];
//			this->AddEdge(u,p[u],this->GetEdgeWeight(u,p[u]));
			this->AddEdge(u_id,v_id,this->GetEdgeWeight(u_id,v_id));
			if (u_id < v_id){
				edgeWeightsToKeep.push_back(pair<int,int>(u_id,v_id));
			} else {
				edgeWeightsToKeep.push_back(pair<int,int>(v_id,u_id));
			}
		}
	}
//	cout << "-----------------------------" << endl;
	this->UpdateMaxDegree();
	
	
//	cout << "maximum degree is " << this->maxDegree << endl;
//	cout << "largest element in weights is " << *max_element(weights,weights+this->edgeWeightsMap.size()) << endl;
	if (this->maxDegree == 0) {
		ofstream edgeListFile;
		edgeListFile.open(sequenceFileName+".debugEdgeList");
		for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
			edgeListFile << edgeAndWeight.first.first << "\t";
			edgeListFile << edgeAndWeight.first.second << "\t";
			edgeListFile << edgeAndWeight.second << endl;
		}
//	cout << "number of edges to keep is " << edgeWeightsToKeep.size() << endl;	
//	cout << "largest element in weights is " << *max_element(weights,weights+this->edgeWeightsMap.size()) << endl;
		edgeListFile.close();
	}
	
	delete[] edges;
	delete[] weights;
	for (pair<pair<int,int>,int> edgeAndWeight : this->edgeWeightsMap) {
		if (find(edgeWeightsToKeep.begin(),edgeWeightsToKeep.end(),edgeAndWeight.first) == edgeWeightsToKeep.end()){
			edgeWeightsToRemove.push_back(edgeAndWeight.first);
		}
	}
	for (pair<int,int> edge: edgeWeightsToRemove) {
		this->edgeWeightsMap.erase(edge);
	}
//	cout << "Number of vertices after MST update is " << this->vertexMap->size() << endl;
//	cout << "Number of edges after MST update is " << this->edgeWeightsMap.size() << endl;
//	this->UpdateMaxDegree();
}

tuple <vector<string>,vector<vector<unsigned char>>,vector<int>,vector<vector<int>>> MST_tree::GetCompressedSequencesSiteWeightsAndSiteRepeats(vector<int> vertexIdList){	
	vector <string> names;
	vector <vector<unsigned char>> compressedSequences;
	vector <int> sitePatternWeights_ptr;
	vector <vector <int>> sitePatternRepeats_ptr;
	vector <vector<unsigned char>> distinctPatterns;
	map <vector<unsigned char>,vector<int>> distinctPatternsToSitesWherePatternRepeats;
	vector <MST_vertex*> vertexPtrList;
	for (unsigned int i = 0; i < vertexIdList.size(); i++){		
		MST_vertex* v_ptr = (*this->vertexMap)[vertexIdList[i]];
		vertexPtrList.push_back(v_ptr);
		vector <unsigned char> compressedSequence;
		compressedSequences.push_back(compressedSequence);
		names.push_back(v_ptr->name);
	}
	int numberOfSites = vertexPtrList[0]->sequence.size();
	vector<unsigned char> sitePattern;
	for(int site=0; site < numberOfSites; site++){
		sitePattern.clear();
		for (MST_vertex* v_ptr: vertexPtrList){
			sitePattern.push_back(v_ptr->sequence[site]);}
		if (find(distinctPatterns.begin(),distinctPatterns.end(),sitePattern)!=distinctPatterns.end()){
			distinctPatternsToSitesWherePatternRepeats[sitePattern].push_back(site);
			
		} else {
			distinctPatterns.push_back(sitePattern);	
			vector<int> sitePatternRepeats;
			sitePatternRepeats.push_back(site);
			distinctPatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
			for (unsigned int i = 0; i < sitePattern.size(); i++){
				compressedSequences[i].push_back(sitePattern[i]);
			}
		}
	}
	for (vector<unsigned char> sitePattern : distinctPatterns){
		int sitePatternWeight = distinctPatternsToSitesWherePatternRepeats[sitePattern].size();
		sitePatternWeights_ptr.push_back(sitePatternWeight);		
		sitePatternRepeats_ptr.push_back(distinctPatternsToSitesWherePatternRepeats[sitePattern]);
	}
	return make_tuple(names,compressedSequences,sitePatternWeights_ptr,sitePatternRepeats_ptr);
}

void MST_tree::SetCompressedSequencesAndSiteWeightsForInputSequences() {				
	vector <vector<unsigned char>> distinctPatterns;
	map <vector<unsigned char>,vector<int>> distinctPatternsToSitesWherePatternRepeats;	
	int numberOfSites = (*this->vertexMap)[0]->sequence.size();
	int numberOfInputSequences = this->vertexMap->size();
	int sitePatternWeight; int v_id; int site;
	vector<unsigned char> sitePattern;
	for(site=0; site < numberOfSites; site++) {
		sitePattern.clear();
		for (v_id = 0; v_id < numberOfInputSequences; v_id ++) {
			sitePattern.push_back((*this->vertexMap)[v_id]->sequence[site]);
			}
		if (find(distinctPatterns.begin(),distinctPatterns.end(),sitePattern)!=distinctPatterns.end()) {
			distinctPatternsToSitesWherePatternRepeats[sitePattern].push_back(site);			
		} else {
			distinctPatterns.push_back(sitePattern);	
			vector<int> sitePatternRepeats;
			sitePatternRepeats.push_back(site);
			distinctPatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
			for (v_id = 0; v_id < numberOfInputSequences; v_id ++) {				
				(*this->vertexMap)[v_id]->globallyCompressedSequence.push_back(sitePattern[v_id]);
			}
		}
	}
	for (vector<unsigned char> sitePattern: distinctPatterns) {
		sitePatternWeight = distinctPatternsToSitesWherePatternRepeats[sitePattern].size();		
		this->siteWeights.push_back(sitePatternWeight);		
	}
}

void MST_tree::WriteToFile(string FileName) {
	ofstream mstFile;	
	mstFile.open(FileName);	
	MST_vertex * v;
	for (pair <int, MST_vertex *> vIdAndPtr: *this->vertexMap) {
		v = vIdAndPtr.second;
		for (MST_vertex * n: v->neighbors) {
			if (v->id < n->id) {
				mstFile << v->name << "\t" << n->name << "\t" << this->GetEdgeWeight(v->id, n->id) << endl; 
			}
		}
	}
	mstFile.close();
}

unsigned char MST_tree::ConvertDNAToChar(char dna) {
	string dna_upper = string(1,toupper(dna));
	unsigned char dna_char = 4;
	if (this->mapDNAtoInteger.find(dna_upper) != this->mapDNAtoInteger.end()) {
		dna_char = this->mapDNAtoInteger[dna_upper];
	} else {
		if (isspace(dna)) {
			cout << "DNA character is a whitespace" << endl;
		}
		cout << "DNA character " << dna_upper << " is not in dictionary keys" << endl;
	}	
	return (dna_char);
}


void MST_tree::ComputeMST_nonACGT() {



}

void MST_tree::ReadSequences(string sequenceFileNameToSet){
	this->sequenceFileName = sequenceFileNameToSet;
	vector <unsigned char> recodedSequence;
	recodedSequence.clear();
	unsigned int site = 0;
	unsigned char dna_char;
	int num_amb = 0;
	int num_non_amb = 0;
//	cout << "Sequence file name is " << this->sequenceFileName << endl;
	ifstream inputFile(this->sequenceFileName.c_str());
	string seqName;
	// vector mut_list;
	string seq = "";	
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {
//				sequenceNames.push_back(seqName);
				for (char const dna: seq) {
					if (!isspace(dna)) {
						dna_char = this->ConvertDNAToChar(dna);
						if (dna_char > 3) { // FIX_AMB
							num_amb += 1;
							dna_char = 3;
						} else {
							num_non_amb += 1;
						}
						recodedSequence.push_back(dna_char);					
						site += 1;							
						}						
					// dna_char = this->mapDNAtoInteger[string(1,toupper(dna))];										
				}
				this->AddVertex(seqName,recodedSequence);
				recodedSequence.clear();
			} 
			seqName = line.substr(1,line.length());
			seq = "";
			site = 0;			
		}
		else {
			seq += line ;
		}		
	}		
	for (char const dna: seq) {
		if (!isspace(dna)) {
			// dna_char = this->mapDNAtoInteger[string(1,toupper(dna))];
			dna_char = this->ConvertDNAToChar(dna);
			if (dna_char > 3) { // FIX_AMB
				num_amb += 1;
				dna_char = 3;
			} else {
				num_non_amb += 1;
			}
			recodedSequence.push_back(dna_char);
			site += 1;
		}
	}
	this->AddVertex(seqName,recodedSequence);
	recodedSequence.clear();
//	sequenceNames.push_back(seqName);
	inputFile.close();
	cout << "Number of ambiguous characters is " << float(num_amb) << "\tNumber of nonambiguous characters is " << float(num_non_amb) << endl;
	cout << "Fraction of ambiguous characters is " << float(num_amb)/float(num_amb + num_non_amb) << endl;
}

void MST_tree::ComputeMST() {	

// Construct MST incrementally
// Create distance matrix 

	int numberOfVertices = (this->v_ind);		
	const int numberOfEdges = numberOfVertices*(numberOfVertices-1)/2;		
	
	int * weights;
	weights = new int [numberOfEdges];
		
	int edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {			
			weights[edgeIndex] = ComputeHammingDistance((*this->vertexMap)[i]->sequence,(*this->vertexMap)[j]->sequence);
			edgeIndex += 1;
		}
	}
	typedef pair <int,int > E;

	E * edges;
	edges = new E [numberOfEdges];
	edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {
			edges[edgeIndex] = E(i,j);
			edgeIndex += 1;
		}
	}
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_distance_t, int>, boost::property < boost::edge_weight_t, int> > Graph;
	Graph g(edges, edges + numberOfEdges, weights, numberOfVertices);

	vector < boost::graph_traits < Graph >::vertex_descriptor >  p(num_vertices(g));
	prim_minimum_spanning_tree(g, &p[0]);
	delete[] edges;		
	int edgeCount = 0;
//	ofstream MSTFile;
//	MSTFile.open(sequenceFileName+".mst");
	for (size_t u = 0; u != p.size(); u++) {
		if (p[u] != u) {
			edgeCount += 1;
			if (u < p[u]) {
				edgeIndex = GetEdgeIndex(u,p[u],numberOfVertices);
			} else {
				edgeIndex = GetEdgeIndex(p[u],u,numberOfVertices);
			}
			this->AddEdge(u, p[u], weights[edgeIndex]);
//			MSTFile << (*this->vertexMap)[u]->name << "\t" << (*this->vertexMap)[p[u]]->name << "\t" << weights[edgeIndex] << endl;
		}
	}
	this->UpdateMaxDegree();
//	MSTFile.close();
	delete[] weights;
}

void MST_tree::ComputeChowLiuTree() {
// UNREST
// Fit Q
// Fit Q




}


#endif
