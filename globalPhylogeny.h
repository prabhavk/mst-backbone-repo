#ifndef phylogeny_H
#define phylogeny_H
#include <boost/algorithm/string.hpp>
#include <tuple>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include"utilities.h"

using namespace Eigen;

// vertex 
class vertex {
public:	
	int degree;
	int inDegree;
	int outDegree;
	int timesVisited;
	float vertexLogLikelihood;	
	float sumOfEdgeLogLikelihoods;
	bool observed;
	std::string name;
	std::string newickLabel;
	std::vector <unsigned char> fullSequence;
	std::vector <unsigned char> compressedSequence;
	vertex * parent;
	std::vector <vertex *> neighbors;	
	std::vector <vertex *> children;
	void AddNeighbor(vertex * v);
	void AddParent(vertex * v);
	void AddChild(vertex * v);
	vertex (std::string v_name, std::vector <unsigned char> sequenceToAdd) {
		this->name = v_name;
		this->fullSequence = sequenceToAdd;
		this->newickLabel = "";
		this->degree = 0;
		this->inDegree = 0;
		this->outDegree = 0;
		this->timesVisited = 0;
		this->vertexLogLikelihood = 0;
		this->sumOfEdgeLogLikelihoods = 0;
		this->observed = 0;
	}
	~vertex () {
		
	}
};

void vertex::AddNeighbor(vertex * v) {
	this->neighbors.push_back(v);
	this->degree += 1;
}

void vertex::AddParent(vertex * v) {
	this->parent = v;
	this->inDegree += 1;
}

void vertex::AddChild(vertex * v) {
	this->children.push_back(v);
	this->outDegree += 1;
}

// tree

class globalPhylogeny {
public:
	vertex * root;
	std::string sequenceFileName;
	std::vector <int> siteWeights;
	std::map <std::pair<vertex *, vertex *>, float> edgeLengths;
	std::map <std::pair<vertex *, vertex *>, float> edgeLogLikelihoods;
	float sequenceLength;
	void CompressSequencesAndSetSiteWeights();
	void SetExpectedCounts();
	std::map <std::pair<vertex *, vertex *>,Matrix4f> expectedCountsForVertexPair;
	std::map <vertex *, std::array <float, 4>> expectedCountsForVertex;
	void PerformModelSelection();
	void AddSequences(std::map<std::string,std::vector<unsigned char>> compressedSequencesList, std::vector <std::vector <int>> sitePatternRepeats);	
	void ParseEdgeListString(std::string weightedEdgeList);
	void AddWeightedEdges(std::vector<std::tuple<std::string,std::string,float>> weightedEdges);
	void SelectVertexForRooting();	
	void ComputeEdgeLoglikelihood(vertex * u, vertex * v);
	void ComputeEdgeLength(vertex * u, vertex * v);
	std::array <float, 4> GetProbabilityDistribution(vertex * u);
	Matrix4f GetTransitionProbability(vertex * u, vertex * v);
	void ComputeVertexLoglikelihood(vertex * u);
	void AddVertex(std::string v_name, std::vector<unsigned char> fullSequence);
	bool ContainsVertex(std::string v_name);
	void AddUndirectedEdge(vertex * u, vertex * v, float edgeLength);
	void AddDirectedEdge(vertex * u, vertex * v, float edgeLength);
	void RootTree();
	void SetEdgesForTreeTraversal();
	void ComputeNewickLabel();
	void WriteNewickLabel();
	void WriteEdges();
	std::vector <std::pair<vertex *, vertex *>> edgesForPreOrderTreeTraversal;
	std::vector <std::pair<vertex *, vertex *>> edgesForPostOrderTreeTraversal;
	std::vector <vertex *> preOrderVerticesWithoutLeaves;
	std::vector <vertex *> leaves;	
	std::vector <vertex *> nonLeafVertices;
	void ResetParentChildRelationships();
	void SetEdgesForPreOrderTraversal();
	void SetVerticesForPreOrderTraversalWithoutLeaves();
	void SetEdgesForPostOrderTraversal();
	void SetLeavesAndNonLeafVertices();	
	float GetEdgeLength(vertex * u, vertex * v);
	void ResetTimesVisited();
	std::map <std::string,vertex *> vertices;
	globalPhylogeny (std::string sequenceFileNameToSet) {
		this->sequenceFileName = sequenceFileNameToSet;
		
	}
	~globalPhylogeny () {
		for (std::pair<std::string,vertex *> nameAndPtr : this->vertices) {
			delete nameAndPtr.second;			
		}
		this->vertices.clear();
	}
};

float globalPhylogeny::GetEdgeLength(vertex * u, vertex * v) {
	float t;
	if (u < v) {
		t = this->edgeLengths[std::make_pair(u,v)];
	} else {
		t = this->edgeLengths[std::make_pair(v,u)];
	}
	return (t);
}

std::array <float, 4> globalPhylogeny::GetProbabilityDistribution(vertex * u) {
	std::array <float, 4> pi;
	for (int i = 0; i < 4; i ++) {
		pi[i] = 0.0;
	}
	int dna;
//	cout << "Sequence length is " << this->sequenceLength << endl;
//	cout << "Full sequence is " << EncodeAsDNA(u->fullSequence) << endl;
	for (int site = 0; site < (int) this->sequenceLength; site ++) {
		dna = u->fullSequence[site];
		pi[dna] += 1.0 ;
	}		
	for (int i = 0; i < 4; i ++) {
		pi[i] /= this->sequenceLength;
	}
	return (pi);
}

Matrix4f globalPhylogeny::GetTransitionProbability(vertex * u, vertex * v) {
	Matrix4f P = ArrayXXf::Zero(4,4);
	int dna_u; int dna_v;
	for (int site = 0; site < (int) this->sequenceLength; site ++) {
		dna_u = u->fullSequence[site];
		dna_v = v->fullSequence[site];
		P(dna_u,dna_v) += 1.0;
	}
	float rowSum;
	for (int i = 0; i < 4; i ++) {
		rowSum = 0;
		for (int j = 0; j < 4; j ++) {
			rowSum += P(i,j);
		}
		for (int j = 0; j < 4; j ++) {
			P(i,j) /= rowSum;
		}
	}
	return P;
}

void globalPhylogeny::ComputeEdgeLength(vertex * u, vertex * v) {
	float edgeLength = 0;
	int dna_u; int dna_v;
	for (int site = 0; site < (int) this->sequenceLength; site ++) {
		dna_u = u->fullSequence[site];
		dna_v = v->fullSequence[site];
		edgeLength += 1.0;		
	}
	std::cout << "Edge length before normalizing is " << edgeLength << std::endl;
	edgeLength /= this->sequenceLength;
	std::cout << "Edge length after normalizing is " <<  edgeLength << std::endl;
	this->edgeLengths.insert(std::make_pair(std::make_pair(u,v),edgeLength));
	return;	
}

void globalPhylogeny::ComputeEdgeLoglikelihood(vertex * u, vertex * v) {
	float edgeLoglikelihood = 0;	
	Matrix4f P = this->GetTransitionProbability(u,v);
	int dna_u; int dna_v;
	for (int site = 0; site < (int) this->sequenceLength; site ++) {
		dna_u = u->fullSequence[site];
		dna_v = v->fullSequence[site];
		edgeLoglikelihood += log(P(dna_u,dna_v));		
	}
	this->edgeLogLikelihoods.insert(std::make_pair(std::make_pair(u,v),edgeLoglikelihood));
	return;
}

void globalPhylogeny::ComputeVertexLoglikelihood(vertex * u) {
	std::array <float, 4> pi = this->GetProbabilityDistribution(u);
	u->vertexLogLikelihood = 0;
	int dna;
	for (int site = 0; site < (int) this->sequenceLength; site ++ ) {
		dna = u->fullSequence[site];
		u->vertexLogLikelihood += log(pi[dna]);
	}
	return;
}

void globalPhylogeny::AddVertex(std::string v_name, std::vector <unsigned char> fullSequence) {
	vertex * v = new vertex(v_name, fullSequence);	
	this->vertices[v_name] = v;
	if(boost::algorithm::starts_with(v_name, "h_")){
		v->observed = 0;		
	} else {
		v->observed = 1;
		v->newickLabel = v->name;
	}	
	return;
}

bool globalPhylogeny::ContainsVertex(std::string v_name) {
	bool returnValue;
	if (this->vertices.find(v_name) == this->vertices.end()) {
		returnValue = 0;
	} else {
		returnValue = 1;
	}
	return (returnValue);
}


void globalPhylogeny::AddSequences(std::map<std::string,std::vector<unsigned char>> sequencesList, std::vector <std::vector <int>> sitePatternRepeats) {
	std::string v_name;
	std::vector <unsigned char> compressedSequence;
	std::vector <unsigned char> fullSequence;	
	for (std::pair <std::string, std::vector<unsigned char>> nameAndSeq : sequencesList) {
		v_name = nameAndSeq.first;		
		if (!this->ContainsVertex(v_name)) {
			compressedSequence = nameAndSeq.second;
			fullSequence = DecompressSequence(&compressedSequence, &sitePatternRepeats);
			this->AddVertex(v_name, fullSequence);
			this->ComputeVertexLoglikelihood(this->vertices[v_name]);
		}
	}
}

void globalPhylogeny::AddWeightedEdges(vector<tuple<string,string,float>> weightedEdges) {
	string u_name; string v_name; float t;
	vertex * u; vertex * v;
//	cout << "Adding weighted edges" << endl;
	for (tuple <string, string, float> weightedEdge : weightedEdges) {		
		tie (u_name, v_name, t) = weightedEdge;
//		cout << "Attempting to add edge " << u_name << "\t" << v_name << endl;		
		u = this->vertices[u_name];
//		cout << "Found " << u->name << endl;		
		v = this->vertices[v_name];
//		cout << "Found " << v->name << endl;
		u->AddNeighbor(v);
		v->AddNeighbor(u);
//		cout << "Edge length for " << u->name << "\t" << v->name << "\t is \t"<< t << endl;
		if (u < v) {
			this->edgeLengths[pair<vertex *,vertex *>(u,v)] = t;
		} else {
			this->edgeLengths[pair<vertex *,vertex *>(v,u)] = t;
		}		
		this->ComputeEdgeLoglikelihood(u,v);
		this->ComputeEdgeLoglikelihood(v,u);
	}	
}

// Compress sequences and set site weights
void globalPhylogeny::CompressSequencesAndSetSiteWeights() {
	vector <vector<unsigned char>> distinctPatterns;
	map <vector<unsigned char>, vector<int>> distinctPatternsToSitesWherePatternRepeats;
	vector <vertex *> vertexList;
	vertex * v;
	for (pair <string, vertex *> namePtrPair : this->vertices) {
		v = namePtrPair.second;
		if (v->observed) {
			vertexList.push_back(v);
		}		
	}	
	int numberOfObservedSequences = vertexList.size();
//	cout << "number of observed sequences is " << numberOfObservedSequences << endl;
	int sitePatternWeight;
	vector <unsigned char> sitePattern;
	for (int site=0; site < (int) this->sequenceLength; site++) {
		sitePattern.clear();
		for (int v_id = 0; v_id < (int) vertexList.size(); v_id ++) {
			sitePattern.push_back(vertexList[v_id]->fullSequence[site]);
		}
//		if (site == 0) {
//			cout << EncodeAsDNA(sitePattern) << endl;
//		}		
		if (find(distinctPatterns.begin(),distinctPatterns.end(),sitePattern)!=distinctPatterns.end()){
			distinctPatternsToSitesWherePatternRepeats[sitePattern].push_back(site);			
		} else {
			distinctPatterns.push_back(sitePattern);	
			vector <int> sitePatternRepeats;
			sitePatternRepeats.push_back(site);
			distinctPatternsToSitesWherePatternRepeats[sitePattern] = sitePatternRepeats;						
			for (int v_id = 0; v_id < (int) vertexList.size(); v_id ++){				
				vertexList[v_id]->compressedSequence.push_back(sitePattern[v_id]);
			}
		}
	}
	for (vector<unsigned char> sitePattern: distinctPatterns){
		sitePatternWeight = distinctPatternsToSitesWherePatternRepeats[sitePattern].size();		
		this->siteWeights.push_back(sitePatternWeight);		
	}
//	cout << "Number of distinct site patterns is " << this->siteWeights.size() << endl;
}

void globalPhylogeny::RootTree() {
	this->SetLeavesAndNonLeafVertices();	
//	cout << "Number of leaf vertices is " << this->leaves.size() << endl;
//	cout << "Number of non-leaf vertices is " << this->nonLeafVertices.size() << endl;
	this->SelectVertexForRooting();
	this->WriteEdges();
	this->WriteNewickLabel();
	// If you're fitting w.r.t 
	// this->CompressSequencesAndSetSiteWeights();
}

void globalPhylogeny::SelectVertexForRooting() {
	// Select vertex that maximizes vertex logLikelihood + sum of edge logLikelihoods;
	float logLikelihoodForRooting = 0;
	float maxLogLikelihoodForRooting = -1 * pow(10,10);
	vertex * p; vertex * c;
	vertex * finalRoot;
	finalRoot = this->nonLeafVertices[0];		
	for (vertex * v: this->nonLeafVertices) {
		this->root = v;
		this->SetEdgesForTreeTraversal();		
		logLikelihoodForRooting = this->root->vertexLogLikelihood;		
		for (pair <vertex *, vertex *> edge : this->edgesForPreOrderTreeTraversal) {			
			logLikelihoodForRooting += this->edgeLogLikelihoods[edge];
		}
		if (logLikelihoodForRooting > maxLogLikelihoodForRooting) {
			maxLogLikelihoodForRooting = logLikelihoodForRooting;
			finalRoot = this->root;
		}		
	}	
	this->root = finalRoot;
	this->SetEdgesForTreeTraversal();
}

void globalPhylogeny::ResetParentChildRelationships() {
	vertex * v;
	for (pair <string, vertex *> namePtrPair : this->vertices) {
		v = namePtrPair.second;
		v->parent = v;
		v->children.clear();
		v->inDegree = 0;
		v->outDegree = 0;
	}	
}

void globalPhylogeny::SetEdgesForTreeTraversal() {	
	this->SetEdgesForPreOrderTraversal();
	this->SetVerticesForPreOrderTraversalWithoutLeaves();
	this->SetEdgesForPostOrderTraversal();
	this->ResetParentChildRelationships();
//	cout << "Number of edges in pre order tree traversal operation are " << endl;
//	cout << this->edgesForPreOrderTreeTraversal.size() << endl;
	vertex * p; vertex * c;
	for (pair <vertex *, vertex *> edge : this->edgesForPreOrderTreeTraversal) {
		tie (p,c) = edge;
		p->AddChild(c);
		c->AddParent(p);
	}
}

void globalPhylogeny::SetEdgesForPreOrderTraversal() {	
	this->edgesForPreOrderTreeTraversal.clear();
	vector <vertex *> verticesToVisit;
	vertex * p;	
	verticesToVisit.push_back(this->root);	
	int numberOfVerticesToVisit = verticesToVisit.size();
	int totalNumberOfedgesVisited = 0;
	while (numberOfVerticesToVisit > 0) {
		p = verticesToVisit[numberOfVerticesToVisit-1];
		p->timesVisited += 1;
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (vertex * c : p->neighbors) {
			totalNumberOfedgesVisited += 1;						
//			cout << "Name = " << c->name << "\t" << "Times visited = " << c->timesVisited << endl;
			if (c->timesVisited == 0) {				
				this->edgesForPreOrderTreeTraversal.push_back(pair <vertex*, vertex*>(p,c));
//				cout << p->name << "\t" << c->name << endl;
				verticesToVisit.push_back(c);
				numberOfVerticesToVisit += 1;
			}
		}
	}
//	cout << "Number of edges visited is " << totalNumberOfedgesVisited << endl;
	this->ResetTimesVisited();
}

void globalPhylogeny::SetEdgesForPostOrderTraversal() {
	this->edgesForPostOrderTreeTraversal.clear();
	int numberOfEdges = (int) this->edgesForPreOrderTreeTraversal.size();
	pair <vertex *, vertex *> edge;
	for (int edgeInd = numberOfEdges - 1; edgeInd > -1; edgeInd --) {
		edge = this->edgesForPreOrderTreeTraversal[edgeInd];
		this->edgesForPostOrderTreeTraversal.push_back(edge);
	}
}


void globalPhylogeny::SetVerticesForPreOrderTraversalWithoutLeaves() {
	this->preOrderVerticesWithoutLeaves.clear();
	for (pair <vertex*, vertex*> edge : this->edgesForPreOrderTreeTraversal) {
		if (find(this->preOrderVerticesWithoutLeaves.begin(),this->preOrderVerticesWithoutLeaves.end(),edge.first) == this->preOrderVerticesWithoutLeaves.end()){
			this->preOrderVerticesWithoutLeaves.push_back(edge.first);
		}
	}
}

void globalPhylogeny::SetLeavesAndNonLeafVertices() {
	if (this->leaves.size() == 0) {
		this->nonLeafVertices.clear();
		vertex * v;
		for (pair <string, vertex*> namePtrPair : this->vertices) {
			v = namePtrPair.second;
			if (v->observed) {
				this->leaves.push_back(v);
			} else {
				this->nonLeafVertices.push_back(v);
			}
		}
	}
//	cout << "Number of leaves is " << this->leaves.size() << endl; 
//	cout << "Number of non-leaf vertices is " << this->nonLeafVertices.size() << endl;
}

void globalPhylogeny::ResetTimesVisited() {
	vertex * v;
	for (pair <string, vertex*> namePtrPair: this->vertices){
		v = namePtrPair.second;
		v->timesVisited = 0;		
	}
}

void globalPhylogeny::ComputeNewickLabel() {		
	this->ResetTimesVisited();
	vector <vertex *> verticesToVisit = leaves;
	vertex * c;
	vertex * p;
	float edgeLength;
//	cout << "Number of edges for pre order tree traversal = " << this->edgesForPreOrderTreeTraversal.size() << endl;
//	cout << "Number of edges for post order tree traversal = " << this->edgesForPostOrderTreeTraversal.size() << endl;
	for (pair <vertex *, vertex *> edge : this->edgesForPostOrderTreeTraversal) {
		tie (p,c) = edge;
		p->timesVisited += 1;
		edgeLength = this->GetEdgeLength(p,c);
		if (p->timesVisited == p->outDegree) {
			p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength) + ")";
		} else if (p->timesVisited == 1) {
			p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
		} else {			
			p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength);
		}
	}
	this->root->newickLabel += ";";
	this->ResetTimesVisited();
}

void globalPhylogeny::WriteNewickLabel() {
	this->ComputeNewickLabel();
	ofstream newickFile;
	newickFile.open(this->sequenceFileName + ".newick");
	newickFile << this->root->newickLabel + "\n";
	newickFile.close();
}


void globalPhylogeny::WriteEdges() {
	ofstream edgeListFile;
	edgeListFile.open(this->sequenceFileName + ".edges");
	vertex * u; vertex * v; float t;
	for (pair < pair <vertex *, vertex *>, float > edgeAndLengthPair : this->edgeLengths ) {
		t = edgeAndLengthPair.second;
		tie (u, v) = edgeAndLengthPair.first;
		edgeListFile << u->name << "\t" << v->name << "\t" << t << endl;
	}	
//	vertex * u;
//	for (pair <string, vertex *> namePtrPair : this->vertices) {
//		u = namePtrPair.second;
//		for (vertex * v : u->neighbors) {
//			if (u->name < v->name) {
//				edgeListFile << u->name << "\t" << v->name << "\t" << this->edgeLengths[make_pair(u,v)] << endl;
//			}
//		}
//	}
	edgeListFile.close();
}



#endif