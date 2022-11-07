#ifndef rootedPhylogeny_H
#define rootedPhylogeny_H
#include <iostream>
#include "utilities.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/minima.hpp>
//#include <thread>
//#include <pthread.h> 
#include <random>
#include <map>
#include <tuple>
#include <chrono>
#include <iomanip>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "SEM.h"
#include "asa047.hpp"
using namespace Eigen;
//using namespace boost;

//class LogLikelihood : public Problem<double> {
//  public:
//    using typename cppoptlib::Problem<double>::Scalar;
//    using typename cppoptlib::Problem<double>::TVector;
//
//    double value(vector<float> parameters) {
//        
//        return   t1 * t1 + 100 * t2 * t2;
//    }
//    void gradient(const TVector &x, TVector &grad) {
//        grad[0]  = -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
//        grad[1]  = 200 * (x[1] - x[0] * x[0]);
//    }
//};

class rootedPhylogeny_vertex{	
public:
	int parent_id;	
	int numberOfDescendants = 0;
	int timesVisited = 0;
	int rateCategory = -1;
	int optimalRateCategoryForRootedTree = -1;
	int globallyOptimalRateCategoryForRootedTree = -1;
	float logScalingFactors;
	float logScalingFactors_firstDer;
	float logScalingFactors_secondDer;
	vector <int> children_id;
	vector <int> neighbors_id;
	int id;
	string name;	
	string newickLabel;
	vector <unsigned char> sequence;
	vector <unsigned char> compressedSequence;
	array <float, 4> baseFreq;
//	vector <array<int,4>> conditionalLikelihood_int;
	vector <array<double,4>> conditionalLikelihood;
	void AddParent(int p_id);
	void AddChild(int c_id);
	void AddNeighbor(int n_id);
	void AddNumberOfObservedSequences(int numberOfObservedSequencesToAdd);
	rootedPhylogeny_vertex(int idToAdd, string nameToAdd, vector <unsigned char> sequenceToAdd){
		this->id = idToAdd;
		this->parent_id = this->id;
		this->name = nameToAdd;
		this->sequence = sequenceToAdd;
		for (int dna = 0; dna < 4; dna++) {
			this->baseFreq[dna] = 0;
		}
		this->newickLabel = "";
	}
	~rootedPhylogeny_vertex(){
		
	}
};


void rootedPhylogeny_vertex::AddParent(int p_id){
	this->parent_id = p_id;
}

void rootedPhylogeny_vertex::AddChild(int c_id){
	this->children_id.push_back(c_id);
}

void rootedPhylogeny_vertex::AddNeighbor(int n_id) {
	this->neighbors_id.push_back(n_id);
}

class rootedPhylogeny_tree{
//	using namespace boost::multiprecision;
private:
	int numberOfObservedSequences;
	default_random_engine generator;		
	map <int, float> * scalingFactorForRateCategory;
	map <int, MatrixXf> * stationaryDistributionForCategory;
	map <int, int> * changePointForRateCategory;
	rootedPhylogeny_vertex * root;
	map <rootedPhylogeny_vertex*, float> * baseFreqChangeForEachDescendant;
	vector <rootedPhylogeny_vertex*> * changePoints;	
	int rateCategoryForOptimization;
	vector <float> unorderedWeightedSiteLogLikelihoods;
	vector <float> * thresholdList;
	double logLikelihoodConvergenceThreshold = 0.01;
	MatrixXd searchDirection;	
//	boost::multiprecision::cpp_dec_float_100 smallestMaxValueForConditionalLikelihoodOfAncestors;
public:
	int numberOfEdgesTried = 0;
	string sequenceFileName;
	ofstream logFile;
	vector <pair<int, int>> * edgesForRooting;
	int sequenceLength;
	map <string,unsigned char> mapDNAtoInteger;	
	map <int, MatrixXd> * parametersPerRateCategory;
	map <int, Matrix4f> * rateMatrixPerRateCategory;
	map <int, Matrix4f> * optimalRateMatrixPerRateCategory;
	map <int, Matrix4f> * optimalRateMatrixPerRateCategoryForRootedTree;
	map <int, Matrix4f> * globallyOptimalRateMatrixPerRateCategory;
	map <int, Matrix4f> * bestRateMatrixPerRateCategoryOverMultipleEMRuns;
	map <pair<int, int>, float> * optimalEdgeLengths;
	map <pair<int, int>, float> * optimalEdgeLengthsForRootedTree;
	map <pair<int, int>, float> * globallyOptimalEdgeLengths;
	vector <pair <int, int>> * globallyOptimalDirectedEdges;
	map <pair<int, int>, float> * bestEdgeLengthsOverMultipleEMRuns;
	float optimalThreshold;	
	float BICForCurrentThreshold = pow(10,10);	
	float minimum_BIC_for_rooted_tree = pow(10,10);
	float globally_minimum_BIC = pow(10,10);
	float optimal_BIC;
	double logLikelihood;
	array <double, 11> jacobian_dep;
	float scalingFactor;
	MatrixXd JacobianForRateMatrix;
	MatrixXd HessianInverserForRateMatrix;
	MatrixXd HessianDep;
	array <float, 4> rootProbability;
	vector <int> siteWeights;
	Matrix4f rateMatrix;	
	VectorXd initialEstimateForFreeParametersOfQ;
	pair <float, array<float,4>> GetLogAbsLargestElementAndScaledArray(array<float,4> originalArray);
	double GetAbsLargestElement(array<double,4> originalArray);
	void ScaleArrayByLargestElement(array<double,4> * arrayToBeScaled);
	pair <double, array<double,4>> GetLogAbsLargestElementAndScaledArray(array<double,4> originalArray);
	VectorXd initialEstimateForFreeParameters;
	VectorXd freeParametersExcBaseFreq;
	MatrixXd freeParametersIncBaseFreq;
	int numberOfRateCategories;	
	vector <rootedPhylogeny_vertex*> listOfChangePoints;
//	vector <int> * leaf_ids;
	map <pair<int, int>, float> * edgeLengthsMap;
	vector <pair<int, int>> * edgesForPostOrderTreeTraversal;
	vector <rootedPhylogeny_vertex*> * verticesForPreOrderTreeTraversal;
	map <int, rootedPhylogeny_vertex*>* vertexMap;
	vector <rootedPhylogeny_vertex *> * leaves;
//	vector <pair<int,int>> * directedEdgeList;
	float ComputeNTDiff(rootedPhylogeny_vertex* p, rootedPhylogeny_vertex* c);
	float logSmallestValueForDouble;
	float logLargestValueForDouble;
	void AddDirectedEdges(vector <pair<int, int>> * directedEdgeList_ptr);
	void AddEdge(int p_id, int c_id);
	void AddEdgeLength(int p_id, int c_id, float edgeLength);
	void SetEdgeLength(int p_id, int c_id, float edgeLength);
	void RemoveEdges();
	float GetEdgeLength(int u_id, int v_id);
	void ComputeAndSetEdgeLength(int u_id, int v_id);
	void ComputeEdgeLengths();
	void AddVertex(int idToAdd, string nameToAdd, vector <unsigned char> sequenceToAdd);
	void AddVertex(int idToAdd, string nameToAdd, vector <unsigned char> sequenceToAdd, vector <unsigned char> compressdSequenceToAdd);
	void ContractEdge(int idOfVertexToKeep, int idOfVertexToRemove);
	bool IsEdgeContractionFeasbile(int vertex_id_1, int vertex_id_2);
	bool DoesEdgeContractionReduceBIC(int idOfVertexToKeep, int idOfVertexToRemove);
	bool DoesEdgeContractionReduceAIC(int vertex_id_1, int vertex_id_2);
	void ContractOutdegreeZeroLeafIncidentZeroLengthEdges();
	void ContractZeroLengthEdges();
	void ContractEdgesBasedOnAIC();
	void SetEdgeLengthsToZeroSuchThatBICIsMinimized();
	void StoreDirectedEdges();
	void StoreOptimalRateCategoriesForRootedTree();
	void StoreGloballyOptimalRateCategories();
	void RestoreRateCategories();
	void StoreOptimalEdgeLengthsForRootedTree();
	void StoreOptimalRateMatricesForRootedTree();
	void RestoreRateMatrices();
	void StoreGloballyOptimalRateMatrices();
	void StoreGloballyOptimalEdgeLengths();
	void RestoreDirectedEdges();
	void RestoreEdgeLengths();
	float ComputeScalingFactor(Matrix4f Q);
	MatrixXf ComputeStationaryDistribution(Matrix4f Q);
	void AddNumberOfObservedSequences(int numberOfObservedSequencesToAdd);
	void WriteEdgeList(string sequenceFileName);
	void WriteAncestralSequences(string sequenceFileName);
	void ComputeBaseFreq();
	void SetThresholds();
	double local_min_brent();
	static double powell (double x[4]);
//	void testNelderMead_1();
	void testNelderMead_2();
	void WriteNewickFile(string sequenceFileName);
//	void SetLeafIds();
	void ComputeNumberOfDescendants();
	void ReadDirectedEdgeListForBifurcatingTree(string treeFileName);
	void ReadUndirectedEdgeList(string treeFileName);
	void ReadSequenceFile(string sequenceFileName);
	void SetRootProbabilityUsingAncestralSequence(); 
	void SetParameters(VectorXd x);
	void SetSequenceFileName(string sequenceFileName);
	void SetRateMatrixForRateCategory(Matrix4f Q, int rateCategory);
	void AddStationaryDistributionForCategory(MatrixXf stationaryDistribution ,int rateCategory);
	void AddRateCategoryForVertex(int v_id, int rateCategory);	
	void SetChangePointForRateCategory(int v_id, int rateCategory);
//	void ComputeLogLikelihoodUsingAncestralStates();
	MatrixXd GetFreeParametersDep(Matrix4f Q);
	Matrix4f GetRateMatrixForFreeParameters(MatrixXd B);
	void SetRateMatrixUsingFreeParametersDep(MatrixXd B);
	void InitializeModelParametersForNelderMead();
	void InitializeModelParametersForBFGS();
	void OptimizeModelParametersDep();
	void OptimizeModelParametersForRateCatDep(int rateCat);
	void ComputeMLEOfRootProbability();
	void OptimizeModelParametersForAGivenThresholdUsingNelderMead(float threshold);
	void EstimateAncestralSequencesByFittingTheGMMUsingMultiprecision();
	void EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision();
	void EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision();
	void ComputeMPEstimateOfAncestralSequences();
	void ComputeMAPEstimateOfAncestralSequencesUsingLongDouble();
	void ComputeInitialEstimateForRateMatrix();
	void ScaleEdgeLengths();
	int GetVertexId(string v_name);
	void SetMinLengthOfEdges();
	void ComputeMLEOfRateMatrices();
	void ComputeMLEOfRateMatricesForLeafLabeledTrees();
	void ComputeMLEOfEdgeLengths();
	void NelderMeadForPowell(int n, double start[], double xmin[], double *ynewlo,
		 double reqmin, double step[], int konvge, int kcount, 
		 int *icount, int *numres, int *ifault );
	void NelderMeadForOptimizingRateParametersForRateCat(int rateCat, int n, double start[], double xmin[], 
		 double *ynewlo, double reqmin, double step[], int konvge,
		 int kcount, int *icount, int *numres, int *ifault);
	void NelderMeadForOptimizingRateParametersForRateCatForFullyLabeledTree(int rateCat, int n, double start[], double xmin[], 
		 double *ynewlo, double reqmin, double step[], int konvge,
		 int kcount, int *icount, int *numres, int *ifault);
	void NelderMeadForOptimizingParametersOfRootProb(int n, double start[], double xmin[], 
		 double *ynewlo, double reqmin, double step[], int konvge,
		 int kcount, int *icount, int *numres, int *ifault);
	void ComputeJacobianForFullyLabeledTree();
	void ComputeHessianForRateMatrixParametersDep();
	void ComputeInitialEstimateForFreeParametersIncBaseFreq();
	void ComputeConditionalLikelihoods();
	MatrixXd GetJacobianForRateCategory(int rateCategoryForOptimization);
	void ComputeJacobianForRateMatrixForRateCategoryDep(int rateCat);
	MatrixXd GetJacobianOfLogLikelihoodForRootProbParameters();
	MatrixXd GetHessianOfLogLikelihoodForRootProbParameters();
	void ComputeLogLikelihoodUsingStoredConditionals();
	void ComputeLogLikelihoodUsingLongDouble();
	void ComputeLogLikelihood();
	void ComputeLogLikelihoodUsingMultiThreading(int numThreads);
	void ResetWeightedSiteLogLikelihoods();
	void ComputeWeightedSiteLogLikelihood(int site);
	void AddWeightedSiteLogLikelihoods();
	void ResetLogScalingFactors();
	double GetBICForGenerallyLabeledTreeUsingLongDouble();
	void ComputeLogLikelihoodForGenerallyLabeledTreeUsingLongDouble();
	void ComputeLogLikelihoodUsingDouble();
	double GetBIC();
	void PerformModelSelectionUsingNelderMead();
	void PerformModelSelectionUsingMultiThreadedNelderMead();
	void PerformModelSelectionUsingNelderMeadViaEM();
	void FitAHomogenousIrrerversibleMarkovModelViaNelderMead();
	void FitAHomogenousIrrerversibleMarkovModelViaBFGS();
	void FitAHomogenousIrrerversibleMarkovModelViaNelderMeadUsingEM();
	void PerformModelSelectionUsingBFGS();
	bool DoParametersEncodeAFeasibleRateMatrix(MatrixXd parameters);
	double GetNegLogLikelihoodForStepSize (double stepSize);	
	float GetNegLogLikelihoodForEdgeLength (rootedPhylogeny_vertex * c, float t);
	double GetNegLogLikelihoodForRateParametersForRateCatForNelderMead (double x[], int rateCat);
	double GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree (double x[], int rateCat);
	double GetNegLogLikelihoodForParametersOfRootProb (double x[]);
	float ComputeAbsDifferenceInBaseFreq(array <float, 4> baseFreq_1, array <float, 4> baseFreq_2);
	double SampleFunctionForMinimization(double x);
	double BrentLineSearchForSelectingStepSize(double a, double b, double t, double &x);
	void OptimizeEdgeLengthUsingBrentsLineSearch(rootedPhylogeny_vertex* child, float a, float b, float t, float &x);
	double BrentLineSearchForSelectingStepSizeForRateCat(double a, double b, double t, double &x, int rateCat);
	double GetOptimalStepSizeUsingLineSearch(MatrixXd B_current);
	double GetOptimalStepSizeForBFGS();
	double GetNormOfJacobian(MatrixXd JacobianForRateMatrix);
	void ComputeLogLikelihoodForFullyLabeledTree();
	void SetEdgesForRooting();
	void SetEdgesForTreeTraversalOperations();
	void RootTreeAtEdge(int u_id, int v_id);
	int GetNumberOfNonZeroLengthEdges();
	void ResetTimesVisited();
	void ClearDirectedEdges();
	void RemoveLengthsOfEdgesExitingRoot();
	void SetDefaultLengthsOfEdgesLeavingRoot();
	void RootUnrootedTreeAtEachEdgeAndSelectRootingAndModelThatMinimizesBIC();
	void SetOptimalRooting();
	void AddNeighborsForUnrootedTree();
	void SetParametersForRateMatrixForNelderMead(double x[], int rateCat);
	vector <double> GetParametersForRateCat(int rateCat);
	void SetEdgesForPostOrderTreeTraversal();
	void SetEdgeLengths();
	void SetLeaves();
	void SetVerticesForPreOrderTreeTraversal();
	void SetParametersForRateMatrixForBFGS(MatrixXd parameters, int rateCat);
	void SetRateCategories(float baseFreqThreshold);	
	void WriteRateCategoryPerVertexAndModelParameters(string sequenceFileName);
	void OptimizeModelParametersUsingNelderMead();
	void OptimizeModelParametersUsingMultiThreadedNelderMead();
	void OptimizeModelParametersUsingNelderMeadViaEM();
	void OptimizeModelParametersUsingBFGS();
	void OptimizeParametersForRateMatrixUsingNelderMead(int rateCat);
	void OptimizeParametersForRateMatrixUsingNelderMeadForFullyLabeledTree(int rateCat);
	void OptimizeParametersForRateMatrixUsingBFGS(int rateCat);
	void OptimizeEdgeLengthsForRateCat(int rateCat);
	void OptimizeEdgeLengthsForRateCatUsingNRWithScalers(int rateCat);

//	void OptimizeEdgeLengthUsingBrent(pair<int,int> edgeIds, int rateCat);
	void OptimizeEdgeLengthsForRateCatUsingBrent(int rateCat);
	void OptimizeEdgeLengthUsingBrent(rootedPhylogeny_vertex * c);
	void OptimizeEdgeLengthsForRateCatForFullyLabeledTree(int rateCat);
	void OptimizeModelParametersForRootProbUsingNelderMead();
	void OptimizeModelParametersForRootProbUsingNewtonRaphson();
	double ComputeEdgeLogLikelihood(rootedPhylogeny_vertex * p, rootedPhylogeny_vertex * c, Matrix4f P);
	void SetSiteWeights(vector <int> siteWeightsToSet);
	array <double, 4> GetLikelihoodArray(double elem_0, double elem_1, double elem_2, double elem_3);
	void InitializeConditionalLikelihoods();	
	void ResetConditionalLikelihoodsForAncestors();	
	Matrix4f ComputeFirstDerivativeOfRateMatrix(int rateCat, int par);
	Matrix4f ComputeFirstDerivativeOfMatrixExponential(float t, int rateCat, int par); 	
	Matrix4f ComputeFirstDerivativeOfRateMatrixDep(int par);
	Matrix4f ComputeDerivativeOfMatrixExponentialDep(float t, int par);
	Matrix4f ComputeSecondDerivativeOfRateMatrix(int par_1, int par_2);
	Matrix4f ComputeSecondDerivativeOfMatrixExponential(float t, int par_1, int par_2);	
	rootedPhylogeny_tree(){
		this->mapDNAtoInteger.insert(pair<string,unsigned char>("A",0));
		this->mapDNAtoInteger.insert(pair<string,unsigned char>("C",1));
		this->mapDNAtoInteger.insert(pair<string,unsigned char>("G",2));
		this->mapDNAtoInteger.insert(pair<string,unsigned char>("T",3));
		logSmallestValueForDouble = log(numeric_limits<double>::min());
		logLargestValueForDouble = log(numeric_limits<double>::max());
//		this->mapDNAtoInteger["C"] = 1;
//		this->mapDNAtoInteger["G"] = 2;
//		this->mapDNAtoInteger["T"] = 3;
		this->optimal_BIC = pow(10,10);
		this->optimalEdgeLengths = new map <pair<int,int>, float>;		
		this->globallyOptimalEdgeLengths = new map <pair<int,int>, float>;
		this->optimalEdgeLengthsForRootedTree = new map <pair<int,int>, float>;
		this->globallyOptimalDirectedEdges = new vector <pair<int,int>>;
		this->bestEdgeLengthsOverMultipleEMRuns = new map <pair<int,int>, float>;
		this->leaves = new vector <rootedPhylogeny_vertex *>;
		this->thresholdList = new vector <float>;
		this->initialEstimateForFreeParameters = VectorXd(11);
		this->freeParametersExcBaseFreq = VectorXd(11);
		this->rateMatrix = ArrayXXf::Zero(4,4);		
		this->JacobianForRateMatrix = ArrayXXd::Zero(11,1);
		this->HessianDep = ArrayXXd::Zero(11,11);
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		this->generator = default_random_engine(seed);	
		vector <unsigned char> root_sequence;
		this->verticesForPreOrderTreeTraversal = new vector <rootedPhylogeny_vertex*>;
		this->baseFreqChangeForEachDescendant = new map <rootedPhylogeny_vertex*, float>;
		this->root = new rootedPhylogeny_vertex(-1, "h_root", root_sequence);
		this->vertexMap = new map <int, rootedPhylogeny_vertex*>;
		this->changePoints= new vector <rootedPhylogeny_vertex*>;
		this->vertexMap->insert(pair<int, rootedPhylogeny_vertex*>(-1, this->root));
		this->edgesForPostOrderTreeTraversal = new vector<pair<int, int>>;
		this->edgeLengthsMap = new map <pair<int, int>, float>;
		this->parametersPerRateCategory = new map <int, MatrixXd>;
		this->rateMatrixPerRateCategory = new map <int, Matrix4f>;
		this->optimalRateMatrixPerRateCategory = new map <int, Matrix4f>;
		this->optimalRateMatrixPerRateCategoryForRootedTree = new map <int, Matrix4f>;		
		this->globallyOptimalRateMatrixPerRateCategory = new map <int, Matrix4f>;
		this->bestRateMatrixPerRateCategoryOverMultipleEMRuns = new map <int, Matrix4f>;
		this->stationaryDistributionForCategory = new map <int, MatrixXf>;
		this->changePointForRateCategory = new map <int, int>;
		this->scalingFactorForRateCategory = new map <int, float>;
		this->edgesForRooting = new vector <pair<int, int>>;
	}
	~rootedPhylogeny_tree(){
		for (pair<int,rootedPhylogeny_vertex*> idAndVertexPtrPair : (*this->vertexMap)){
			delete idAndVertexPtrPair.second;
		}
		delete this->vertexMap;
		delete this->thresholdList;
		delete this->leaves;
		delete this->optimalEdgeLengths;
		delete this->bestEdgeLengthsOverMultipleEMRuns;
		delete this->edgesForPostOrderTreeTraversal;
		delete this->edgeLengthsMap;
		delete this->parametersPerRateCategory;
		delete this->rateMatrixPerRateCategory;
		delete this->optimalRateMatrixPerRateCategory;
		delete this->optimalRateMatrixPerRateCategoryForRootedTree;
		delete this->optimalEdgeLengthsForRootedTree;
		delete this->globallyOptimalEdgeLengths;
		delete this->globallyOptimalDirectedEdges;
		delete this->globallyOptimalRateMatrixPerRateCategory;
		delete this->bestRateMatrixPerRateCategoryOverMultipleEMRuns;
		delete this->stationaryDistributionForCategory;
		delete this->verticesForPreOrderTreeTraversal;
		delete this->changePointForRateCategory;
		delete this->scalingFactorForRateCategory;	
		delete this->baseFreqChangeForEachDescendant;
		delete this->changePoints;
		delete this->edgesForRooting;
	}
};


void rootedPhylogeny_tree::ScaleArrayByLargestElement(array <double, 4> * arrayToBeScaled) {
	double largestElement = this->GetAbsLargestElement(*arrayToBeScaled);
	for (int i = 0; i < 4; i ++) {
		(*arrayToBeScaled)[i] /= largestElement;
	}
	return;
}

pair <float, array <float,4>> rootedPhylogeny_tree::GetLogAbsLargestElementAndScaledArray(array <float,4> originalArray) {
	float largestElement = 0;	
	for (int i = 0; i < 4; i ++) {
		if (largestElement < originalArray[i]){
			largestElement = originalArray[i];
		}
	}
	array <float,4> rescaledArray;
	if (largestElement != 0) {
		for (int i = 0; i < 4; i ++) {
		rescaledArray[i] = originalArray[i]/largestElement;
		}
	} else {
		cout << "Element with largest absolute value is zero" << endl;
		exit(-1);
	}		
	return make_pair(log(largestElement),rescaledArray);
}

pair <double, array <double,4>> rootedPhylogeny_tree::GetLogAbsLargestElementAndScaledArray(array <double,4> originalArray) {
	double largestElement = 0;	
	for (int i = 0; i < 4; i ++) {
		if (largestElement < abs(originalArray[i])){
			largestElement = abs(originalArray[i]);
		}
	}
	array <double,4> rescaledArray;
	if (largestElement != 0) {
		for (int i = 0; i < 4; i ++) {
		rescaledArray[i] = originalArray[i]/largestElement;
		}
	} else {
		cout << "Element with largest absolute value is zero" << endl;
		exit(-1);
	}		
	return make_pair(log(largestElement),rescaledArray);
}

double rootedPhylogeny_tree::GetAbsLargestElement(array <double, 4> originalArray) {
	double largestElement = 0;
	for (int i = 0; i < 4; i ++) {
		if (largestElement < abs(originalArray[i])) {
			largestElement = abs(originalArray[i]);
		}
	}
	return (largestElement);
}

void rootedPhylogeny_tree::SetSequenceFileName(string sequenceFileNameToSet) {
	this->sequenceFileName = sequenceFileNameToSet;
}

void rootedPhylogeny_tree::RootUnrootedTreeAtEachEdgeAndSelectRootingAndModelThatMinimizesBIC() {
	pair <int, int> currentEdge;
	this->logFile.open(sequenceFileName + ".modelSelection_log");
	cout << "Performing model selection by rooting tree at each edge " << endl;
	this->logFile << "Performing model selection by rooting tree at each edge " << endl;
	this->logFile << "u_name\tv_name\tBIC" << endl;	
	pair <string,string> optimalEdgeForRooting;
	string skipEdge;
	int numberOfEdgesTried = 0;
	for (pair <int, int> edge : * this->edgesForRooting) {
		numberOfEdgesTried += 1;
		cout << "Rooting along edge " << numberOfEdgesTried << " of " << this->edgesForRooting->size() << ":\t" << (*this->vertexMap)[edge.first]->name << "\t" << (*this->vertexMap)[edge.second]->name << endl;		
//		if ((*this->vertexMap)[edge.first]->name == "T6" or (*this->vertexMap)[edge.second]->name == "T6") {
//		cout << "Should edge be skipped? enter y or n" << endl;
//		cin >> skipEdge ;
		skipEdge = "n";
		if (skipEdge == "n") {
			this->RootTreeAtEdge(edge.first, edge.second);		
			this->PerformModelSelectionUsingNelderMead();
			cout << (*this->vertexMap)[edge.first]->name << "\t" << (*this->vertexMap)[edge.second]->name << "\t" << this->minimum_BIC_for_rooted_tree << endl;
			this->logFile << (*this->vertexMap)[edge.first]->name << "\t" << (*this->vertexMap)[edge.second]->name << "\t" << this->minimum_BIC_for_rooted_tree << endl;
			if (this->globally_minimum_BIC > this->minimum_BIC_for_rooted_tree or numberOfEdgesTried == 1) {
				this->globally_minimum_BIC = this->minimum_BIC_for_rooted_tree;
				this->StoreDirectedEdges();
				this->StoreGloballyOptimalEdgeLengths();
				this->StoreGloballyOptimalRateMatrices();
				this->StoreGloballyOptimalRateCategories();
				optimalEdgeForRooting.first = (*this->vertexMap)[edge.first]->name;
				optimalEdgeForRooting.second = (*this->vertexMap)[edge.second]->name;
	//			}			
			}
		}		
	}
	// Set topology and edge lengths for optimal rooted tree
	this->RestoreDirectedEdges();
	this->RestoreEdgeLengths();
	this->RestoreRateCategories();
	this->RestoreRateMatrices();	
	this->logFile << "Optimal edge for rooting is " << optimalEdgeForRooting.first << "\t" << optimalEdgeForRooting.second << endl;
	cout << "Optimal edge for rooting is " << optimalEdgeForRooting.first << "\t" << optimalEdgeForRooting.second << endl;
	this->logFile.close();
}

void rootedPhylogeny_tree::SetEdgesForRooting() {
	this->edgesForRooting->clear();
	int p_id; int c_id;
//	rootedPhylogeny_vertex * p;
//	rootedPhylogeny_vertex * c;
	for (pair<pair<int,int>,float> edgeIdAndLength : *this->edgeLengthsMap) {
		p_id = edgeIdAndLength.first.first;
		c_id = edgeIdAndLength.first.second;		
		this->edgesForRooting->push_back(make_pair(p_id,c_id));
	}
//	cout << "Number of edges for rooting is " << this->edgesForRooting->size() << endl;
}

void rootedPhylogeny_tree::StoreGloballyOptimalEdgeLengths() {
	this->globallyOptimalEdgeLengths->clear();
	for (pair<pair<int,int>,float> edgeIdAndLength : * this->optimalEdgeLengthsForRootedTree) {
		this->globallyOptimalEdgeLengths->insert(edgeIdAndLength);
	}
}

void rootedPhylogeny_tree::StoreOptimalEdgeLengthsForRootedTree() {
	this->optimalEdgeLengthsForRootedTree->clear();
	float t;
	int p_id; int c_id;
	rootedPhylogeny_vertex * c;	
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;
		c_id = c->id;
		p_id = c->parent_id;
		if (p_id != c_id) {
			t = this->GetEdgeLength(p_id,c_id);
			if (p_id < c_id) {
				this->optimalEdgeLengthsForRootedTree->insert(make_pair(make_pair(p_id,c_id),t));
			} else {
				this->optimalEdgeLengthsForRootedTree->insert(make_pair(make_pair(c_id,p_id),t));
			}
		}
	}
}

void rootedPhylogeny_tree::StoreOptimalRateMatricesForRootedTree() {
	this->optimalRateMatrixPerRateCategoryForRootedTree->clear();		
	for (pair <int, Matrix4f> rateCatRateMatrixPair : *this->rateMatrixPerRateCategory){
		this->optimalRateMatrixPerRateCategoryForRootedTree->insert(rateCatRateMatrixPair);
	}
}

void rootedPhylogeny_tree::StoreGloballyOptimalRateMatrices() {
	this->globallyOptimalRateMatrixPerRateCategory->clear();
//	cout << "Number of rate matrices stored is " << this->optimalRateMatrixPerRateCategoryForRootedTree->size() << endl;
	for (pair <int, Matrix4f> rateCatRateMatrixPair : *this->optimalRateMatrixPerRateCategoryForRootedTree){
		this->globallyOptimalRateMatrixPerRateCategory->insert(rateCatRateMatrixPair);
	}	
}

void rootedPhylogeny_tree::RestoreRateMatrices() {
	this->rateMatrixPerRateCategory->clear();
	for (pair <int, Matrix4f> rateCatRateMatrixPair : *this->globallyOptimalRateMatrixPerRateCategory){
		this->rateMatrixPerRateCategory->insert(rateCatRateMatrixPair);
	}
}

void rootedPhylogeny_tree::StoreDirectedEdges() {
	this->globallyOptimalDirectedEdges->clear();
	int p_id; int c_id;
	rootedPhylogeny_vertex * c;	
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;
		c_id = c->id;
		p_id = c->parent_id;
		if (p_id != c_id) {
			this->globallyOptimalDirectedEdges->push_back(pair<int,int>(p_id,c_id));
		}
	}
}

void rootedPhylogeny_tree::StoreOptimalRateCategoriesForRootedTree() {
	rootedPhylogeny_vertex * v;
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->optimalRateCategoryForRootedTree = v->rateCategory;
	}
}

void rootedPhylogeny_tree::StoreGloballyOptimalRateCategories() {
	rootedPhylogeny_vertex * v;
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->globallyOptimalRateCategoryForRootedTree = v->optimalRateCategoryForRootedTree;
	}
}
void rootedPhylogeny_tree::RestoreRateCategories() {
	rootedPhylogeny_vertex * v;
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->rateCategory = v->globallyOptimalRateCategoryForRootedTree;
	}
}

void rootedPhylogeny_tree::RestoreEdgeLengths() {
	this->edgeLengthsMap->clear();
	float t;
	int p_id; int c_id;
	rootedPhylogeny_vertex * c;	
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		c = idPtrPair.second;
		c_id = c->id;
		p_id = c->parent_id;
		if (p_id != c->id) {			
			if (p_id < c_id) {
				t = (*this->globallyOptimalEdgeLengths)[make_pair(p_id,c_id)];
			} else {
				t = (*this->globallyOptimalEdgeLengths)[make_pair(c_id,p_id)];
			}
			this->AddEdgeLength(p_id,c_id,t);
		}
	}
}

void rootedPhylogeny_tree::RestoreDirectedEdges(){
	this->ClearDirectedEdges();
	rootedPhylogeny_vertex * p; rootedPhylogeny_vertex * c;
	for (pair<int,int> parentChildId : *this->globallyOptimalDirectedEdges) {
		p = (*this->vertexMap)[parentChildId.first];
		c = (*this->vertexMap)[parentChildId.second];
		p->AddChild(c->id);
		c->AddParent(p->id);
	}
}

void rootedPhylogeny_tree::RemoveLengthsOfEdgesExitingRoot() {
	if (this->root->children_id.size() == 2) {
		int c_l_id = this->root->children_id[0];
		int c_r_id = this->root->children_id[1];	
		this->edgeLengthsMap->erase(make_pair(-1,c_l_id));
		this->edgeLengthsMap->erase(make_pair(-1,c_r_id));		
	}
}

void rootedPhylogeny_tree::SetDefaultLengthsOfEdgesLeavingRoot() {
	int c_l_id = this->root->children_id[0];
	int c_r_id = this->root->children_id[1];
	float default_edge_length = 0.001;
	this->edgeLengthsMap->insert(make_pair(make_pair(-1,c_l_id),default_edge_length));
	this->edgeLengthsMap->insert(make_pair(make_pair(-1,c_r_id),default_edge_length));
}

void rootedPhylogeny_tree::ClearDirectedEdges() {
	this->RemoveLengthsOfEdgesExitingRoot();	
	rootedPhylogeny_vertex * v;	
	for (pair <int, rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;
		v->parent_id = v->id;
		v->children_id.clear();		
//		v->inDegree = 0;
//		v->outDegree = 0;
		v->timesVisited = 0;
	}
}

void rootedPhylogeny_tree::ResetTimesVisited() {
	rootedPhylogeny_vertex * v;
	for (pair <int, rootedPhylogeny_vertex *> idPtrPair : * this->vertexMap) {
		v = idPtrPair.second;		
		v->timesVisited = 0;		
	}
}

void rootedPhylogeny_tree::SetEdgesForTreeTraversalOperations() {
	this->ResetTimesVisited();
	this->SetVerticesForPreOrderTreeTraversal();
	this->SetEdgesForPostOrderTreeTraversal();

}


void rootedPhylogeny_tree::RootTreeAtEdge(int u_id, int v_id){
//	int i = 0;	
	this->ClearDirectedEdges();
	this->root->AddChild(u_id);
	this->root->AddChild(v_id);
	rootedPhylogeny_vertex * u = (*this->vertexMap)[u_id];
	rootedPhylogeny_vertex * v = (*this->vertexMap)[v_id];
	u->AddParent(-1);
	v->AddParent(-1);
	rootedPhylogeny_vertex * n;
	u->timesVisited = 1; v->timesVisited = 1;	
	vector <rootedPhylogeny_vertex *> verticesToVisit;
	verticesToVisit.push_back(u); verticesToVisit.push_back(v);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0) {
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (int n_id : v->neighbors_id) {
			n = (*this->vertexMap)[n_id];
			if (n->timesVisited == 0) {
				v->AddChild(n->id);
				n->AddParent(v->id);
				n->timesVisited = 1;
				verticesToVisit.push_back(n);
				numberOfVerticesToVisit += 1;
			}
		}
	}
	this->SetDefaultLengthsOfEdgesLeavingRoot();
	this->SetEdgesForTreeTraversalOperations();
}

void rootedPhylogeny_tree::ComputeLogLikelihoodForGenerallyLabeledTreeUsingLongDouble() {
	
}

bool rootedPhylogeny_tree::DoParametersEncodeAFeasibleRateMatrix(MatrixXd parameters){
	// check if pi_1, pi_2 and pi_3 are feasible
//	cout << "Parameters being checked for feasibility are " << endl;
//	cout << parameters << endl;
	bool feasibility = 1;
	for (int i = 0; i < 3; i++) {
		if (parameters(i,0) <= 0.0 or parameters(i,0) >= 1.0) { 
			feasibility = 0;
		}
	}
	// check if pi_4 is feasible
	if (1 - (parameters(0,0) + parameters(1,0) + parameters(2,0)) <= 0.0) {
		feasibility = 0;
	} else if (1 - (parameters(0,0) + parameters(1,0) + parameters(2,0)) >= 1.0) {
		feasibility = 0;
	}
	// check if a through h are feasible
	for (int i = 3; i < 11; i++) {
		if (parameters(i,0) <= 0.0) { 
			feasibility = 0;
		}
	}
	// check if i, j, k and l are feasible
	Matrix4f Q = this->GetRateMatrixForFreeParameters(parameters);
//	cout << "Generated rate matrix is " << endl;
//	cout << Q << endl;
	for (int i = 0; i < 4; i ++){
		for (int j = 0; j < 4; j ++){
			if (i != j){
				if (Q(i,j) <= 0.0){
					feasibility = 0;
				}
			}
		}
	}
	return (feasibility);
}

int rootedPhylogeny_tree::GetNumberOfNonZeroLengthEdges(){
	int numOfNonZeroLengthEdges = 0;
	for (pair<pair<int,int>,float> edgeIdsAndLengthPair : *this->edgeLengthsMap){
		if (edgeIdsAndLengthPair.second > 0){
			numOfNonZeroLengthEdges += 1;
		}
	}
	return (numOfNonZeroLengthEdges);
}

void rootedPhylogeny_tree::SetRootProbabilityUsingAncestralSequence(){
	for (int dna = 0; dna < 4; dna ++){
		this->rootProbability[dna] = 0;
	}
	unsigned char dna;
	float sequenceLength = 0;
	for (int site = 0; site < (int) this->siteWeights.size(); site++){
		dna = this->root->compressedSequence[site];
		this->rootProbability[dna] += this->siteWeights[site];		
		sequenceLength += this->siteWeights[site];
	}
	for (int dna = 0; dna < 4; dna ++){
		this->rootProbability[dna] /= sequenceLength;
	}	
}

bool rootedPhylogeny_tree::DoesEdgeContractionReduceBIC(int idOfVertexToKeep, int idOfVertexToRemove) {	
//	vector <tuple<int,int,float>> edgesToRemove;
//	vector <tuple<int,int,float>> edgesToAdd;
//	rootedPhylogeny_vertex * k = (*this->vertexMap)[idOfVertexToKeep];
//	rootedPhylogeny_vertex * r = (*this->vertexMap)[idOfVertexToRemove];
//	float edgeLengthToAdd;
//	int ind;
	// Compute log likelihood for generally labeled tree
	// Remove edges
	// Add edges
	// Compute log likelihood for generally labeled tree
	// Check if BIC increases
//	if (r->parent_id == k->id) {
//		// Case 1: k is parent of r
//		// Remove r from list of children of k
//		ind = find(k->children_id.begin(),k->children_id.end(),r->id) - k->children_id.begin();
//		edgesToRemove.insert(make_tuple(k->id, r->id, this->GetEdgeLength(k->id,r->id)));
//		k->children_id.erase(k->children_id.begin()+ind);
//		for (int childOfR_id: r->children_id) {
//			// Set parent of childOfR to k
//			(*this->vertexMap)[childOfR_id]->parent_id = k->id;
//			// Add edge from k to childOfR
//			k->children_id.push_back(childOfR_id);
//			edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->id,childOfR_id)];
//			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->id,childOfR_id),edgeLengthToAdd));			
//			edgesToAdd.insert(make_tuple(k->id, childOfR_id, this->GetEdgeLength(r->id,childOfR_id)));
//			// Remove edge from r to child of r
////			this->edgeLengthsMap->erase(pair<int,int>(r->id,childOfR_id));
//			edgesToRemove.insert(make_tuple(r->id, childOfR_id, this->GetEdgeLength(r->id,childOfR_id)));
//		}		
//	} else {
//		// Case 2: k is child of r		
//		for (int childOfR_id: r->children_id) {
//			if (childOfR_id != k->id) {
//				// Set parent of childOfR to k
//				(*this->vertexMap)[childOfR_id]->parent_id = k->id;
//				// Add edge from k to childOfR
//				k->children_id.push_back(childOfR_id);
//				edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->id,childOfR_id)];
//				this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->id,childOfR_id),edgeLengthToAdd));
//				// Remove edge from r to child of r
////				this->edgeLengthsMap->erase(pair<int,int>(r->id,childOfR_id));
//				edgesToRemove.insert(make_tuple(r->id, childOfR_id, this->GetEdgeLength(r->id,childOfR_id)));
//			}
//		}
//		// Set parent of k to parent of r
//		k->parent_id = r->parent_id;
//		// If k is not the root then:
//		if (k->parent_id != -1) {
//			// Add edge from parent of k to k
//			edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->parent_id,r->id)];
//			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->parent_id,k->id),edgeLengthToAdd));
//			// Remove edge from parent of r to r
//			this->edgeLengthsMap->erase(pair<int,int>(r->parent_id,r->id));	
//			edgesToRemove.insert(make_tuple(r->parent_id,r->id, this->GetEdgeLength(r->parent_id,r->id)));
//		}		
//	}
//	// Remove all children of r
//	r->children_id.clear();
//	// Set parent of r to r
//	r->parent_id = r->id;
	bool F = 0;
	return (F);
}

vector <double> rootedPhylogeny_tree::GetParametersForRateCat(int rateCat) {
	vector <double> parametersToReturn;
	for (int i = 0; i < 11; i++){
		parametersToReturn.push_back(0);
	}
	Matrix4f Q;
	Q = (*this->rateMatrixPerRateCategory)[rateCat];
	// a
	parametersToReturn[0] = Q(0,1);
	// b
	parametersToReturn[1] = Q(0,2);
	// c
	parametersToReturn[2] = Q(0,3);
	// d
	parametersToReturn[3] = Q(1,0);
	// e
	parametersToReturn[4] = Q(1,2);
	// f
	parametersToReturn[5] = Q(1,3);
	// g
	parametersToReturn[6] = Q(2,0);
	// h
	parametersToReturn[7] = Q(2,1);
	// i
	parametersToReturn[8] = Q(2,3);
	// j
	parametersToReturn[9] = Q(3,0);
	// k
	parametersToReturn[10] = Q(3,1);	
	return (parametersToReturn);
}

void rootedPhylogeny_tree::SetLeaves(){
	this->leaves->clear();
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		if (idPtrPair.second->children_id.size() == 0){
			this->leaves->push_back(idPtrPair.second);
		}
	}
}

void rootedPhylogeny_tree::WriteRateCategoryPerVertexAndModelParameters(string fileNamePrefix){
	ofstream modelParametersFile;
	modelParametersFile.open(fileNamePrefix+".modelParameters");
	modelParametersFile << "Vertex name" << "\t" << "Rate category" << endl;
	rootedPhylogeny_vertex * v;
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : (*this->vertexMap)){
		v = idPtrPair.second;
		modelParametersFile << v->name << "\t" << v->rateCategory << endl;
	}
	modelParametersFile << "==============================================" << endl;
	modelParametersFile << "Root probability:" << endl;
	// A : 0, C : 1, G : 2, T : 3
	string dnaString = "ACGT";
	for (int i = 0; i < 4; i++){
		modelParametersFile << "P(" << dnaString[i] << ") = " ;
		modelParametersFile << this->rootProbability[i] << "\t";
	}
	modelParametersFile << endl;
	modelParametersFile << "==============================================" << endl;
	modelParametersFile << "Rate matrix parameters:" << endl;
	int rateCat; Matrix4f Q;
	for (pair<int,Matrix4f> catAndRateMatrixPair : *this->rateMatrixPerRateCategory){
		rateCat = catAndRateMatrixPair.first;
		Q = catAndRateMatrixPair.second;
		Q /= this->ComputeScalingFactor(Q);
		modelParametersFile << "Rate matrix for category " << rateCat << " is " << endl;
		modelParametersFile << Q << endl;
	}
	modelParametersFile << "==============================================" << endl;
	modelParametersFile.close();	
}

float rootedPhylogeny_tree::ComputeNTDiff(rootedPhylogeny_vertex* p, rootedPhylogeny_vertex* c){
	float ntDiff = 0;
	float seqLength = 0;
	unsigned char dna_p; unsigned char dna_c;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		dna_p = p->compressedSequence[site];
		dna_c = c->compressedSequence[site];
		if (dna_p != dna_c){
			ntDiff += this->siteWeights[site];
		}
		seqLength += this->siteWeights[site];
	}
	ntDiff /= seqLength;
	return (ntDiff);
}

void rootedPhylogeny_tree::SetEdgeLengths(){
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	float edgeLength;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			edgeLength = this->ComputeNTDiff(p,c);
			if (p->id < c->id){
				(*this->edgeLengthsMap)[pair<int,int>(p->id,c->id)] = edgeLength;
			} else {
				(*this->edgeLengthsMap)[pair<int,int>(c->id,p->id)] = edgeLength;
			}
		}	
	}
}

void rootedPhylogeny_tree::PerformModelSelectionUsingBFGS(){
	// Compute ancestral states by fitting the GMM
//	this->EstimateAncestralSequencesByFittingTheGMM();
//	this->ComputeMPEstimateOfAncestralSequences();
	this->SetLeaves();
	// Set edge lengths
	this->SetEdgeLengths();	
//	cout << this->GetEdgeLength(this->GetVertexId("h_root"),this->GetVertexId("h_11")) << endl;
	this->SetMinLengthOfEdges();
	this->SetVerticesForPreOrderTreeTraversal();
	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;
	double current_BIC = 0;
	this->optimal_BIC = 0;
	for (float threshold : (*this->thresholdList)){
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		cout << "Trying threshold " << noOfThresholdsTried << "\t";
		cout << "Number of rate categories is " << this->numberOfRateCategories << endl;
		this->InitializeModelParametersForBFGS();
		this->OptimizeModelParametersUsingBFGS();
		current_BIC = this->GetBIC();
		cout << "Current BIC is " << current_BIC << endl;
		if (this->optimal_BIC > current_BIC or noOfThresholdsTried == 1) {
			this->optimalRateMatrixPerRateCategory->clear();
			this->optimalEdgeLengths->clear();
			this->optimal_BIC = current_BIC;
			this->optimalThreshold = threshold;
			for (pair<int, Matrix4f> rateCatRateMatrixPair : *this->rateMatrixPerRateCategory){
				this->optimalRateMatrixPerRateCategory->insert(rateCatRateMatrixPair);
			}
			for (pair< pair <int,int>, float> edgeIdsAndLength : *this->edgeLengthsMap){
				this->optimalEdgeLengths->insert(edgeIdsAndLength);
			}
//			cout << "number of rate categories for optimal threshold is " << this->numberOfRateCategories << endl;
//			cout << "optimal BIC is " << this->optimal_BIC << endl;
		} else {
			this->SetRateCategories(this->optimalThreshold);
//			this->OptimizeModelParametersForAGivenThresholdUsingNelderMead(this->optimalThreshold);
			break;
		}
	}
}

void rootedPhylogeny_tree::PerformModelSelectionUsingMultiThreadedNelderMead(){
	//	this->EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision();
	//	this->EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision();	
	this->ComputeMPEstimateOfAncestralSequences();	
	// Set edge lengths
	this->SetEdgeLengths();	
	this->SetMinLengthOfEdges();
	//	this->SetVerticesForPreOrderTreeTraversal();
	//	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	int optimalNumberOfRateCategories = 0;
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;
	this->minimum_BIC_for_rooted_tree=pow(10,10);
	for (float threshold : (*this->thresholdList)) {
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		cout << "Trying threshold " << noOfThresholdsTried << "\t";
		cout << "Number of rate categories is " << this->numberOfRateCategories << endl;
		this->InitializeModelParametersForNelderMead();
		this->OptimizeModelParametersUsingNelderMead();
		this->BICForCurrentThreshold = this->GetBIC();
		cout << "Current BIC is " << this->BICForCurrentThreshold << endl;
		if (this->minimum_BIC_for_rooted_tree > this->BICForCurrentThreshold or noOfThresholdsTried == 1) {
			this->minimum_BIC_for_rooted_tree = this->BICForCurrentThreshold;
			optimalNumberOfRateCategories = this->numberOfRateCategories;
			this->StoreOptimalEdgeLengthsForRootedTree();
			this->StoreOptimalRateMatricesForRootedTree();
			this->StoreOptimalRateCategoriesForRootedTree();
			this->StoreGloballyOptimalEdgeLengths();
			this->StoreGloballyOptimalRateMatrices();
			this->StoreGloballyOptimalRateCategories();
		} else {
			if (optimalNumberOfRateCategories > 1) {
				cout << "Optimal number of rate categories is more than one" << endl;
			}
			break;	
		}
	}
	this->RestoreRateMatrices();
	this->RestoreRateCategories();
	this->RestoreEdgeLengths();
}

void rootedPhylogeny_tree::PerformModelSelectionUsingNelderMead(){
	//	this->EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision();
	//	this->EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision();	
	this->ComputeMPEstimateOfAncestralSequences();	
	// Set edge lengths
	this->SetEdgeLengths();	
	this->SetMinLengthOfEdges();
	//	this->SetVerticesForPreOrderTreeTraversal();
	//	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	int optimalNumberOfRateCategories = 0;
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;
	this->minimum_BIC_for_rooted_tree=pow(10,10);
	for (float threshold : (*this->thresholdList)) {
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		cout << "Trying threshold " << noOfThresholdsTried << "\t";
		cout << "Number of rate categories is " << this->numberOfRateCategories << endl;
		this->InitializeModelParametersForNelderMead();
		this->OptimizeModelParametersUsingNelderMead();
		this->BICForCurrentThreshold = this->GetBIC();
		cout << "Current BIC is " << this->BICForCurrentThreshold << endl;
		if (this->minimum_BIC_for_rooted_tree > this->BICForCurrentThreshold or noOfThresholdsTried == 1) {
			this->minimum_BIC_for_rooted_tree = this->BICForCurrentThreshold;
			optimalNumberOfRateCategories = this->numberOfRateCategories;
			this->StoreOptimalEdgeLengthsForRootedTree();
			this->StoreOptimalRateMatricesForRootedTree();
			this->StoreOptimalRateCategoriesForRootedTree();
			this->StoreGloballyOptimalEdgeLengths();
			this->StoreGloballyOptimalRateMatrices();
			this->StoreGloballyOptimalRateCategories();
		} else {
			if (optimalNumberOfRateCategories > 1) {
				cout << "Optimal number of rate categories is more than one" << endl;
			}
			break;	
		}
	}
	this->RestoreRateMatrices();
	this->RestoreRateCategories();
	this->RestoreEdgeLengths();
}

void rootedPhylogeny_tree::FitAHomogenousIrrerversibleMarkovModelViaNelderMead() {
//	this->EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision();
//	this->EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision();	
	this->ComputeMPEstimateOfAncestralSequences();	
	// Set edge lengths
	this->SetEdgeLengths();	
	this->SetMinLengthOfEdges();
//	this->SetVerticesForPreOrderTreeTraversal();
//	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;	
	for (float threshold : (*this->thresholdList)) {
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		cout << "Trying threshold " << noOfThresholdsTried << "\t";
		cout << "Number of rate categories is " << this->numberOfRateCategories << endl;
		this->InitializeModelParametersForNelderMead();
		this->OptimizeModelParametersUsingNelderMead();		
		this->BICForCurrentThreshold = this->GetBIC();
//		cout << "Current BIC is " << current_BIC << endl;
		if (this->minimum_BIC_for_rooted_tree > this->BICForCurrentThreshold or noOfThresholdsTried == 1) {
			this->minimum_BIC_for_rooted_tree = this->BICForCurrentThreshold;
			this->StoreOptimalEdgeLengthsForRootedTree();
			this->optimalRateMatrixPerRateCategory->clear();		
			this->optimalThreshold = threshold;	
			for (pair< int, Matrix4f> rateCatRateMatrixPair : *this->rateMatrixPerRateCategory){
				this->optimalRateMatrixPerRateCategory->insert(rateCatRateMatrixPair);
			}
			this->SetRateCategories(this->optimalThreshold);
		} else {
			break;	
		}		
	}
}

void rootedPhylogeny_tree::FitAHomogenousIrrerversibleMarkovModelViaNelderMeadUsingEM(){
//	this->EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision();	
//	this->EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision();	
	this->ComputeMPEstimateOfAncestralSequences();
	this->SetLeaves();
	// Set edge lengths
	this->SetEdgeLengths();	
	this->SetMinLengthOfEdges();
	this->SetVerticesForPreOrderTreeTraversal();
	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;
	double current_BIC = 0;
	this->optimal_BIC = 0;	
	for (float threshold : (*this->thresholdList)) {
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		cout << "Trying threshold " << noOfThresholdsTried << "\t";
		cout << "Number of rate categories is " << this->numberOfRateCategories << endl;
		this->InitializeModelParametersForNelderMead();
		this->OptimizeModelParametersUsingNelderMeadViaEM();
//		current_BIC = this->GetBIC();
//		cout << "Current BIC is " << current_BIC << endl;
//		if (this->optimal_BIC > current_BIC or noOfThresholdsTried == 1){			
		this->optimalRateMatrixPerRateCategory->clear();
		this->optimal_BIC = current_BIC;
		this->optimalThreshold = threshold;	
		for (pair< int, Matrix4f> rateCatRateMatrixPair : *this->rateMatrixPerRateCategory){
			this->optimalRateMatrixPerRateCategory->insert(rateCatRateMatrixPair);
		}
		this->SetRateCategories(this->optimalThreshold);
		break;
	}
}

void rootedPhylogeny_tree::FitAHomogenousIrrerversibleMarkovModelViaBFGS(){
	this->EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision();	
//	this->EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision();	
//	this->ComputeMPEstimateOfAncestralSequences();
	this->SetLeaves();
	// Set edge lengths
	this->SetEdgeLengths();	
	this->SetMinLengthOfEdges();
	this->SetVerticesForPreOrderTreeTraversal();
	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;
	double current_BIC = 0;
	this->optimal_BIC = 0;	
	for (float threshold : (*this->thresholdList)) {
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		cout << "Trying threshold " << noOfThresholdsTried << "\t";
		cout << "Number of rate categories is " << this->numberOfRateCategories << endl;
		cout << "Initializing model parameters for BFGS" << endl;
		this->InitializeModelParametersForBFGS();		
		cout << "Optimizing model parameters for BFGS" << endl;
		this->OptimizeModelParametersUsingBFGS();		
//		current_BIC = this->GetBIC();
//		cout << "Current BIC is " << current_BIC << endl;
//		if (this->optimal_BIC > current_BIC or noOfThresholdsTried == 1){			
		this->optimalRateMatrixPerRateCategory->clear();
		this->optimal_BIC = current_BIC;
		this->optimalThreshold = threshold;	
		for (pair< int, Matrix4f> rateCatRateMatrixPair : *this->rateMatrixPerRateCategory){
			this->optimalRateMatrixPerRateCategory->insert(rateCatRateMatrixPair);
		}
		this->SetRateCategories(this->optimalThreshold);
		break;
	}
}


void rootedPhylogeny_tree::PerformModelSelectionUsingNelderMeadViaEM(){
	// Compute ancestral states by fitting the GMM
//	this->EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision();
//	this->EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision();
	this->ComputeMPEstimateOfAncestralSequences();
	this->SetLeaves();
	// Set edge lengths
	this->SetEdgeLengths();
	this->SetMinLengthOfEdges();
	this->SetVerticesForPreOrderTreeTraversal();
	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();
	// Compute list of thresholds
	this->SetThresholds();
	// Iterate over thresholds in order of decreasing value
	int noOfThresholdsTried = 0;
	double current_BIC = 0;
	this->optimal_BIC = 0;
	bool initializeBICAndParameters = 0;
	int maxIter = 5;
	double bestBICForMultipleEMRuns = 0;
	for (float threshold : (*this->thresholdList)) {
		initializeBICAndParameters = 0;
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		cout << "Trying threshold " << noOfThresholdsTried << "\t";
		cout << "Number of rate categories is " << this->numberOfRateCategories << endl;		
		if (noOfThresholdsTried == 1){
			initializeBICAndParameters = 1;
		}
		// Try 10 times
		for (int iter = 0; iter < maxIter; iter++) {
			this->InitializeModelParametersForNelderMead();
			this->OptimizeModelParametersUsingNelderMeadViaEM();
			current_BIC = this->GetBIC();
			if (bestBICForMultipleEMRuns > current_BIC or iter == 0) {
				bestBICForMultipleEMRuns = current_BIC;
				this->bestRateMatrixPerRateCategoryOverMultipleEMRuns->clear();
				this->bestEdgeLengthsOverMultipleEMRuns->clear();
			}
			for (pair< int, Matrix4f> rateCatRateMatrixPair : *this->rateMatrixPerRateCategory) {
				this->bestRateMatrixPerRateCategoryOverMultipleEMRuns->insert(rateCatRateMatrixPair);
			}
			for (pair< pair <int,int>, float> edgeIdsAndLength : *this->edgeLengthsMap) {
				this->bestEdgeLengthsOverMultipleEMRuns->insert(edgeIdsAndLength);
			}
		}		
//		current_BIC = this->GetBIC();
		cout << "Minimal BIC for current threshold is " << bestBICForMultipleEMRuns << endl;
		if (this->optimal_BIC > bestBICForMultipleEMRuns or initializeBICAndParameters) {
			this->optimalRateMatrixPerRateCategory->clear();
			this->optimalEdgeLengths->clear();
			this->optimal_BIC = bestBICForMultipleEMRuns;
			this->optimalThreshold = threshold;
			for (pair< int, Matrix4f> rateCatRateMatrixPair : *this->bestRateMatrixPerRateCategoryOverMultipleEMRuns) {
				this->optimalRateMatrixPerRateCategory->insert(rateCatRateMatrixPair);
			}
			for (pair< pair <int,int>, float> edgeIdsAndLength : *this->bestEdgeLengthsOverMultipleEMRuns) {
				this->optimalEdgeLengths->insert(edgeIdsAndLength);
			}		
		} else {
			this->SetRateCategories(this->optimalThreshold);
			this->edgeLengthsMap->clear();
			for (pair< pair <int,int>, float> edgeIdsAndLength : *this->optimalEdgeLengths) {
				this->edgeLengthsMap->insert(edgeIdsAndLength);
			}
			break;
		}
	}
}


void rootedPhylogeny_tree::OptimizeModelParametersForAGivenThresholdUsingNelderMead(float threshold){
//	this->EstimateAncestralSequencesByFittingTheGMM();
	// Set edge lengths
//	this->SetEdgeLengths();
//	this->SetMinLengthOfEdges();
	// Compute base freq changes
//	this->ComputeBaseFreq();
	this->SetRateCategories(threshold);
	this->InitializeModelParametersForNelderMead();
	this->OptimizeModelParametersUsingNelderMead();	
}

double rootedPhylogeny_tree::GetBIC(){
	float sequenceLength = 0;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		sequenceLength += this->siteWeights[site];
	}
	double numberOfParameters;
//	cout << "current logLikelihood is " << this->logLikelihood << endl;
//	cout << "number of rate categories is " << this->numberOfRateCategories << endl;
	// If no of vertices in the same category as the root is one then 
	rootedPhylogeny_vertex * c_l;
	rootedPhylogeny_vertex * c_r;
	c_l = (*this->vertexMap)[this->root->children_id[0]];
	c_r = (*this->vertexMap)[this->root->children_id[1]];
	if (this->root->rateCategory != c_l->rateCategory and this->root->rateCategory != c_r->rateCategory){
		numberOfParameters = (11.0 * (this->numberOfRateCategories -1)) + 3.0;
	} else {
		numberOfParameters = 11.0 * (this->numberOfRateCategories);
	}
	double penaltyPerParameter = 1.0;
	double optimal_BIC = ( log(sequenceLength) * penaltyPerParameter * numberOfParameters ) - ( 2.0 * this->logLikelihood );	
	return optimal_BIC;
}

void rootedPhylogeny_tree::SetParametersForRateMatrixForNelderMead(double x[], int rateCat){
	for (int par = 0; par < 11; par ++) {
		(*this->parametersPerRateCategory)[rateCat](par,0) = x[par];
	}
	Matrix4f Q;
	// a
	Q(0,1) = x[0];
	// b
	Q(0,2) = x[1];
	// c
	Q(0,3) = x[2];
	// d
	Q(1,0) = x[3];
	// e
	Q(1,2) = x[4];
	// f
	Q(1,3) = x[5];
	// g
	Q(2,0) = x[6];
	// h
	Q(2,1) = x[7];
	// i
	Q(2,3) = x[8];
	// j
	Q(3,0) = x[9];
	// k
	Q(3,1) = x[10];
	// set l to 1
	Q(3,2) = 1.0;
	// D1 
	Q(0,0) = - (Q(0,1) + Q(0,2) + Q(0,3));
	// D2 
	Q(1,1) = - (Q(1,0) + Q(1,2) + Q(1,3));
	// D3 
	Q(2,2) = - (Q(2,0) + Q(2,1) + Q(2,3));
	// D4 
	Q(3,3) = - (Q(3,0) + Q(3,1) + Q(3,2));
	
	(*this->rateMatrixPerRateCategory)[rateCat] = Q;
	float scalingFactor = this->ComputeScalingFactor(Q);
	(*this->scalingFactorForRateCategory)[rateCat] = scalingFactor;
	if (this->root->rateCategory == rateCat){
		MatrixXf stationaryDistribution = this->ComputeStationaryDistribution(Q);
		for (int i = 0; i < 4; i++){
			if (stationaryDistribution(i,0) < 0){
				cout << "Stationary distribution has negative entry" << endl;
				cout << "Rate matrix is " << endl << Q << endl;
			}
			this->rootProbability[i] = stationaryDistribution(i,0);
		}
	}
	return;
}

double rootedPhylogeny_tree::GetNegLogLikelihoodForRateParametersForRateCatForNelderMead(double x[], int rateCat){	
	// Set parameters for rate cat
	// Check if constraints are being satisfied
	double valueToReturn;
	bool setParameters = 1;
	for (int par = 0; par < 11; par++){
		if (x[par] < pow(10,-4) or x[par] > pow(10,4)){
			setParameters = 0;
			valueToReturn = pow(10,20);
//			cout << "parameter out of bounds" << endl;			
		}
	}
	if (setParameters){
		this->SetParametersForRateMatrixForNelderMead(x, rateCat);	
//		this->ComputeLogLikelihoodUsingDouble();
		this->ComputeLogLikelihoodUsingLongDouble();
		valueToReturn = -1 * this->logLikelihood;		
	}
//	cout << "negative log likelihood is " << valueToReturn << endl;
	return (valueToReturn);
}

double rootedPhylogeny_tree::GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree(double x[], int rateCat){
	double valueToReturn;
	bool setParameters = 1;
	for (int par = 0; par < 11; par++){
		if (x[par] < pow(10,-4) or x[par] > pow(10,4)){
			setParameters = 0;
			valueToReturn = pow(10,20);		
		}
	}
	if (setParameters){
		this->SetParametersForRateMatrixForNelderMead(x, rateCat);	
		this->ComputeLogLikelihoodForFullyLabeledTree();		
		valueToReturn = -1 * this->logLikelihood;
	}
	return (valueToReturn);	
}

double rootedPhylogeny_tree::GetNegLogLikelihoodForParametersOfRootProb(double x[]){	
	// Set parameters for root prob
	// Check if constraints are being satisfied
	for (int par = 0; par < 3; par++){
		if (x[par] < 0.0 or x[par] > 1.0){
			return (pow(10,20));
		}
	}
	if ((x[0] + x[1] + x[3]) > 1.0){
		return (pow(10,20));
	}
	this->rootProbability[0] = x[0];
	this->rootProbability[1] = x[1];
	this->rootProbability[2] = x[2];
	this->rootProbability[3] = (1 - this->rootProbability[0] - this->rootProbability[1] - this->rootProbability[2]);
//	this->ComputeLogLikelihoodUsingDouble();
	this->ComputeLogLikelihoodUsingLongDouble();
	// return neg log likelihood
	return (-1 * this->logLikelihood);
}

void rootedPhylogeny_tree::OptimizeModelParametersUsingBFGS(){
	cout << "Optimizing model parameters" << endl;
	double currentLogLikelihood = 0;
	double updatedLogLikelihood = 0;
	bool convergenceNotReached = 1;	
	bool restart = 1;	
	MatrixXd B = ArrayXXd::Zero(11,1);
	MatrixXd B_current = ArrayXXd::Zero(11,1);
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++) {
		if (this->root->rateCategory == rateCat) {
			rootedPhylogeny_vertex * c_l;
			rootedPhylogeny_vertex * c_r;
			c_l = (*this->vertexMap)[this->root->children_id[0]];
			c_r = (*this->vertexMap)[this->root->children_id[1]];
			if (c_l->rateCategory != rateCat and c_r->rateCategory != rateCat) {
				// Construct a modified version of Nelder Mead for the following case
				// no. of vertices that are in the same rate cat as the root equals one		
				this->OptimizeModelParametersForRootProbUsingNelderMead();
			} else {
				convergenceNotReached = 1;
//				this->ComputeLogLikelihoodUsingDouble();				
				this->ComputeLogLikelihoodUsingLongDouble();				
				currentLogLikelihood = this->logLikelihood;				
				cout << "Current loglikelihood is " << this->logLikelihood << endl;
				while (convergenceNotReached) {
					restart = 1;
					// current parameters
					while (restart) {
						this->ComputeLogLikelihoodUsingLongDouble();
						currentLogLikelihood = this->logLikelihood;						
						B_current = (*this->parametersPerRateCategory)[rateCat];
						this->OptimizeParametersForRateMatrixUsingBFGS(rateCat);
						this->ComputeLogLikelihoodUsingLongDouble();				
						updatedLogLikelihood = this->logLikelihood;
						if (updatedLogLikelihood < currentLogLikelihood and abs(updatedLogLikelihood - currentLogLikelihood) > this->logLikelihoodConvergenceThreshold) {
							restart = 1;
							cout << "Restarting" << endl;
							cout << "Current loglikelihood is " << currentLogLikelihood << endl;
							cout << "Updated loglikelihood is " << updatedLogLikelihood << endl;
							this->SetParametersForRateMatrixForBFGS(B_current, rateCat);
//							this->SetParametersForRateMatrixForBFGS(x, rateCat);
						} else {
							restart = 0;
						}
					}
					this->ComputeLogLikelihoodUsingLongDouble();
					currentLogLikelihood = this->logLikelihood;
//					auto time_B = std::chrono::high_resolution_clock::now();
//					cout << "CPU time used for optimizing rate matrix is " << chrono::duration_cast<chrono::seconds>(time_B-time_A).count() << " second(s)\n";
					this->OptimizeEdgeLengthsForRateCat(rateCat);								
//					auto time_C = std::chrono::high_resolution_clock::now();
//					cout << "CPU time used for optimizing edge lengths is " << chrono::duration_cast<chrono::seconds>(time_C-time_B).count() << " second(s)\n";
					this->ComputeLogLikelihoodUsingLongDouble();
					updatedLogLikelihood = this->logLikelihood;
					if (abs(updatedLogLikelihood - currentLogLikelihood) < this->logLikelihoodConvergenceThreshold) {
						convergenceNotReached = 0;
					}					
					currentLogLikelihood = updatedLogLikelihood;
					cout << "Updated loglikelihood is " << this->logLikelihood << endl;
				}			
//				cout << "Edge length optimization complete" << endl;				
//				cout << "updated loglikelihood 2 is " << this->logLikelihood << endl;
			}
		} else {			
			convergenceNotReached = 1;
//			this->ComputeLogLikelihoodUsingDouble();
			this->ComputeLogLikelihoodUsingLongDouble();				
			currentLogLikelihood = this->logLikelihood;
			cout << "Current loglikelihood is " << this->logLikelihood << endl;
			while (convergenceNotReached) {
				restart = 1;
				while (restart) {
//					this->ComputeLogLikelihoodUsingDouble();
					this->ComputeLogLikelihoodUsingLongDouble();
					currentLogLikelihood = this->logLikelihood;
					B_current = (*this->parametersPerRateCategory)[rateCat];
//					currentParameters = this->GetParametersForRateCat(rateCat);
					this->OptimizeParametersForRateMatrixUsingBFGS(rateCat);
					this->ComputeLogLikelihoodUsingLongDouble();
//					this->ComputeLogLikelihoodUsingDouble();
					updatedLogLikelihood = this->logLikelihood;
					if (updatedLogLikelihood < currentLogLikelihood and abs(updatedLogLikelihood - currentLogLikelihood) > this->logLikelihoodConvergenceThreshold) {
						restart = 1;
						cout << "Restarting" << endl;
						cout << "Current loglikelihood is " << currentLogLikelihood << endl;
						cout << "Updated loglikelihood is " << updatedLogLikelihood << endl;
						// Reset parameters
						this->SetParametersForRateMatrixForBFGS(B_current, rateCat);
//						for (int par = 0; par < 11; par ++){
//							x[par] = currentParameters[par];
//						}						
//						this->SetParametersForRateMatrixForBFGS(x, rateCat);
					} else {
						restart = 0;
					}
				}
				this->ComputeLogLikelihoodUsingLongDouble();
//				this->ComputeLogLikelihoodUsingDouble();
				currentLogLikelihood = this->logLikelihood;
				this->OptimizeEdgeLengthsForRateCat(rateCat);								
				this->ComputeLogLikelihoodUsingLongDouble();
//				this->ComputeLogLikelihoodUsingDouble();
				updatedLogLikelihood = this->logLikelihood;
				if (abs(updatedLogLikelihood - currentLogLikelihood) < this->logLikelihoodConvergenceThreshold) {
						convergenceNotReached = 0;
				}
				currentLogLikelihood = updatedLogLikelihood;
				cout << "Updated loglikelihood is " << this->logLikelihood << endl;
			}
		}		
	}
	this->ComputeLogLikelihoodUsingLongDouble();
}


void rootedPhylogeny_tree::OptimizeModelParametersUsingNelderMeadViaEM() {
	double currentLogLikelihoodForEM = 0;
	double currentLogLikelihoodForQT = 0;
	double updatedLogLikelihoodForEM = 0;
	double updatedLogLikelihoodForQT = 0;
	bool convergenceForEMNotReached = 1;	
	bool convergenceForOptQandtNotReached = 1;	
	bool repeat = 1;
	int maxRepeats = 10;
	int numOfRepeats = 0;
	vector <double> currentParameters;
	double x[11];
	for (int i = 0; i < 11; i++) {
		x[i] = 0;
	}
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++) {
		if (this->root->rateCategory == rateCat) {
			rootedPhylogeny_vertex * c_l;
			rootedPhylogeny_vertex * c_r;
			c_l = (*this->vertexMap)[this->root->children_id[0]];
			c_r = (*this->vertexMap)[this->root->children_id[1]];
			if (c_l->rateCategory != rateCat and c_r->rateCategory != rateCat) {												
				this->SetRootProbabilityUsingAncestralSequence();				
			} else {				
				convergenceForEMNotReached = 1;
				this->ComputeLogLikelihoodUsingLongDouble();
				currentLogLikelihoodForEM = this->logLikelihood;
				cout << "Current loglikelihood for EM is ";
				cout << currentLogLikelihoodForEM << endl;								
//				Iterate E-M steps till convergence
				while (convergenceForEMNotReached){
//					Iterate opt Q and opt t till convergence
					convergenceForOptQandtNotReached = 1;
					this->ComputeLogLikelihoodForFullyLabeledTree();
					currentLogLikelihoodForQT = this->logLikelihood;
					cout << "Current loglikelihood for Q" << endl;
					cout << currentLogLikelihoodForQT << endl;
					while (convergenceForOptQandtNotReached) {
						repeat = 1;
						numOfRepeats = 1;
//						MLE estimate of rate matrix
						while (repeat) {							
							numOfRepeats += 1;
							currentParameters = this->GetParametersForRateCat(rateCat);
							this->OptimizeParametersForRateMatrixUsingNelderMeadForFullyLabeledTree(rateCat);
							this->ComputeLogLikelihoodForFullyLabeledTree();
							updatedLogLikelihoodForQT = this->logLikelihood;
							if (updatedLogLikelihoodForQT < currentLogLikelihoodForQT and abs(updatedLogLikelihoodForQT - currentLogLikelihoodForQT) > this->logLikelihoodConvergenceThreshold and numOfRepeats < maxRepeats) {
								repeat = 1;
								cout << "Restarting optimization of rate matrices" << endl;
								cout << "Current loglikelihood is " << currentLogLikelihoodForQT << endl;
								cout << "Updated loglikelihood is " << updatedLogLikelihoodForQT << endl;
								cout << "Difference in logLikelihood is ";
								cout << setprecision(10) << updatedLogLikelihoodForQT - currentLogLikelihoodForQT << endl;
								// reset parameters
								for (int par = 0; par < 11; par ++){
									x[par] = currentParameters[par];
								}
								this->SetParametersForRateMatrixForNelderMead(x, rateCat);
							} else {
								repeat = 0;
							}
						}
						cout << "Updated loglikelihood for Q" << endl;
						cout << updatedLogLikelihoodForQT << endl;
//							MLE estimate of edge lengths
						this->OptimizeEdgeLengthsForRateCatForFullyLabeledTree(rateCat);
						this->ComputeLogLikelihoodForFullyLabeledTree();
						updatedLogLikelihoodForQT = this->logLikelihood;
						cout << "Updated loglikelihood for t" << endl;
						cout << updatedLogLikelihoodForQT << endl;
						if (abs(currentLogLikelihoodForQT - updatedLogLikelihoodForQT) < this->logLikelihoodConvergenceThreshold){
							convergenceForOptQandtNotReached = 0;
						} else {
							currentLogLikelihoodForQT = updatedLogLikelihoodForQT;
						}
					}
//						MAP estimate of states
					this->ComputeMAPEstimateOfAncestralSequencesUsingLongDouble();
					updatedLogLikelihoodForEM = this->logLikelihood;
					cout << "Updated loglikelihood for EM is ";
					cout << updatedLogLikelihoodForEM << endl;
					if (abs(currentLogLikelihoodForEM - updatedLogLikelihoodForEM) < this->logLikelihoodConvergenceThreshold){
						convergenceForEMNotReached = 0;
					} else {
						currentLogLikelihoodForEM = updatedLogLikelihoodForEM;
					}
				}
				cout << "===================================================" << endl;
			}
		} else {
			convergenceForEMNotReached = 1;
			this->ComputeLogLikelihoodUsingLongDouble();
			currentLogLikelihoodForEM = this->logLikelihood;
			cout << "Current loglikelihood for EM is ";
			cout << currentLogLikelihoodForEM << endl;								
//			Iterate E-M steps till convergence
			while (convergenceForEMNotReached) {
//				Iterate opt Q and opt t till convergence
				convergenceForOptQandtNotReached = 1;
				this->ComputeLogLikelihoodForFullyLabeledTree();
				currentLogLikelihoodForQT = this->logLikelihood;
				cout << "Current loglikelihood for Q" << endl;
				cout << currentLogLikelihoodForQT << endl;
				while (convergenceForOptQandtNotReached) {
					repeat = 1;
					numOfRepeats = 1;
//					MLE estimate of rate matrix
					while (repeat) {							
						numOfRepeats += 1;
						currentParameters = this->GetParametersForRateCat(rateCat);
						this->OptimizeParametersForRateMatrixUsingNelderMeadForFullyLabeledTree(rateCat);
						this->ComputeLogLikelihoodForFullyLabeledTree();
						updatedLogLikelihoodForQT = this->logLikelihood;
						if (updatedLogLikelihoodForQT < currentLogLikelihoodForQT and abs(updatedLogLikelihoodForQT - currentLogLikelihoodForQT) > this->logLikelihoodConvergenceThreshold and numOfRepeats < maxRepeats) {
							repeat = 1;
							cout << "Restarting optimization of rate matrices" << endl;
							cout << "Current loglikelihood is " << currentLogLikelihoodForQT << endl;
							cout << "Updated loglikelihood is " << updatedLogLikelihoodForQT << endl;
							cout << "Difference in logLikelihood is ";
							cout << setprecision(10) << updatedLogLikelihoodForQT - currentLogLikelihoodForQT << endl;
							// reset parameters
							for (int par = 0; par < 11; par ++){
								x[par] = currentParameters[par];
							}
							this->SetParametersForRateMatrixForNelderMead(x, rateCat);
						} else {
							repeat = 0;
						}
					}
					cout << "Updated loglikelihood for Q" << endl;
					cout << updatedLogLikelihoodForQT << endl;
//						MLE estimate of edge lengths
					this->OptimizeEdgeLengthsForRateCatForFullyLabeledTree(rateCat);
					this->ComputeLogLikelihoodForFullyLabeledTree();
					updatedLogLikelihoodForQT = this->logLikelihood;
					cout << "Updated loglikelihood for t" << endl;
					cout << updatedLogLikelihoodForQT << endl;
					if (abs(currentLogLikelihoodForQT - updatedLogLikelihoodForQT) < this->logLikelihoodConvergenceThreshold){
						convergenceForOptQandtNotReached = 0;
					} else {
						currentLogLikelihoodForQT = updatedLogLikelihoodForQT;
					}
				}
//					MAP estimate of states
				this->ComputeMAPEstimateOfAncestralSequencesUsingLongDouble();
				updatedLogLikelihoodForEM = this->logLikelihood;
				cout << "Updated loglikelihood for EM is ";
				cout << updatedLogLikelihoodForEM << endl;
				if (abs(currentLogLikelihoodForEM - updatedLogLikelihoodForEM) < this->logLikelihoodConvergenceThreshold){
					convergenceForEMNotReached = 0;
				} else {
					currentLogLikelihoodForEM = updatedLogLikelihoodForEM;
				}
			}
			cout << "===================================================" << endl;
		}
	}
}

void rootedPhylogeny_tree::OptimizeModelParametersUsingNelderMead() {
	bool verbose = 1;
	int timesRestarted = 0;
	if (verbose) {
		cout << "Optimizing model parameters" << endl;
	}	
	double currentLogLikelihood = 0;
	double updatedLogLikelihood = 0;
	bool convergenceNotReached = 1;	
	bool restart = 1;
	vector <double> currentParameters;
	double x[11];
	for (int i = 0; i < 11; i++) {
		x[i] = 0;
	}
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++) {
		if (this->root->rateCategory == rateCat) {
			rootedPhylogeny_vertex * c_l;
			rootedPhylogeny_vertex * c_r;
			c_l = (*this->vertexMap)[this->root->children_id[0]];
			c_r = (*this->vertexMap)[this->root->children_id[1]];
			if (c_l->rateCategory != rateCat and c_r->rateCategory != rateCat) {
				// Construct a modified version of Nelder Mead for the following case
				// no. of vertices that are in the same rate cat as the root equals one						
				this->OptimizeModelParametersForRootProbUsingNelderMead();
			} else {
				auto time_A = std::chrono::high_resolution_clock::now();
				convergenceNotReached = 1;
				this->ComputeLogLikelihood();
				currentLogLikelihood = this->logLikelihood;
				if (verbose) {				
					cout << "Current loglikelihood is " << this->logLikelihood << endl;
				}
				while (convergenceNotReached) {
					restart = 1;
					timesRestarted = 0;
					// current parameters
					while (restart) {
						this->ComputeLogLikelihood();
						currentLogLikelihood = this->logLikelihood;						
						this->OptimizeParametersForRateMatrixUsingNelderMead(rateCat);
						this->ComputeLogLikelihood();
						updatedLogLikelihood = this->logLikelihood;
						if (timesRestarted < 5 and updatedLogLikelihood < currentLogLikelihood and abs(updatedLogLikelihood - currentLogLikelihood) > this->logLikelihoodConvergenceThreshold) {
							restart = 1;
							timesRestarted += 1;
							cout << "Restarting" << endl;
							cout << "Current loglikelihood is " << currentLogLikelihood << endl;
							cout << "Updated loglikelihood is " << updatedLogLikelihood << endl;
							// reset parameters
							for (int par = 0; par < 11; par ++){
								x[par] = currentParameters[par];
							}
							this->SetParametersForRateMatrixForNelderMead(x, rateCat);
						} else {							
							restart = 0;
							currentParameters = this->GetParametersForRateCat(rateCat);
						}
					}
					this->ComputeLogLikelihood();
					currentLogLikelihood = this->logLikelihood;
					auto time_B = std::chrono::high_resolution_clock::now();
					if (verbose) {					
						cout << "CPU time used for optimizing rate matrix is " << chrono::duration<double>(time_B-time_A).count() << " second(s)\n";	
						cout << "Optimizing edge lengths" << endl;
					}					
//					this->OptimizeEdgeLengthsForRateCat(rateCat);
//					this->OptimizeEdgeLengthsForRateCatUsingBrent(rateCat);
					this->OptimizeEdgeLengthsForRateCatUsingNRWithScalers(rateCat);
					auto time_C = std::chrono::high_resolution_clock::now();
					if (verbose) {					
						cout << "CPU time used for optimizing edge lengths is " << chrono::duration<double>(time_C-time_B).count() << " second(s)\n";			
					}
					this->ComputeLogLikelihood();
					updatedLogLikelihood = this->logLikelihood;
					if (updatedLogLikelihood - currentLogLikelihood < this->logLikelihoodConvergenceThreshold) {
//					if (abs(updatedLogLikelihood - currentLogLikelihood) < this->logLikelihoodConvergenceThreshold) {
						convergenceNotReached = 0;						
					}					
					currentLogLikelihood = updatedLogLikelihood;
					if (verbose){
						cout << "Updated loglikelihood is " << this->logLikelihood << endl;
					}					
				}			
//				cout << "Edge length optimization complete" << endl;				
//				cout << "updated loglikelihood 2 is " << this->logLikelihood << endl;
			}
		} else {
			convergenceNotReached = 1;
			this->ComputeLogLikelihood();				
			currentLogLikelihood = this->logLikelihood;
			if (verbose) {
				cout << "Current loglikelihood is " << this->logLikelihood << endl;
			}
			while (convergenceNotReached) {
				restart = 1;
				timesRestarted = 0;
				while (restart) {
					this->ComputeLogLikelihood();
					currentLogLikelihood = this->logLikelihood;
					currentParameters = this->GetParametersForRateCat(rateCat);
					this->OptimizeParametersForRateMatrixUsingNelderMead(rateCat);
					this->ComputeLogLikelihood();
					updatedLogLikelihood = this->logLikelihood;
					if (timesRestarted < 5 and updatedLogLikelihood < currentLogLikelihood and abs(updatedLogLikelihood - currentLogLikelihood) > this->logLikelihoodConvergenceThreshold) {
						restart = 1;
						timesRestarted += 1;
						cout << "Restarting" << endl;
						cout << "Current loglikelihood is " << currentLogLikelihood << endl;
						cout << "Updated loglikelihood is " << updatedLogLikelihood << endl;
						// reset parameters
						for (int par = 0; par < 11; par ++){
							x[par] = currentParameters[par];
						}
						this->SetParametersForRateMatrixForNelderMead(x, rateCat);
					} else {
						restart = 0;
					}
				}
				this->ComputeLogLikelihood();
//				this->ComputeLogLikelihoodUsingDouble();
				currentLogLikelihood = this->logLikelihood;
//				this->OptimizeEdgeLengthsForRateCat(rateCat);	
//				this->OptimizeEdgeLengthsForRateCatUsingBrent(rateCat);							
				this->OptimizeEdgeLengthsForRateCatUsingNRWithScalers(rateCat);
				this->ComputeLogLikelihood();
//				this->ComputeLogLikelihoodUsingDouble();
				updatedLogLikelihood = this->logLikelihood;
				if (verbose) {
					cout << "Updated loglikelihood is " << this->logLikelihood << endl;
				}
				if (updatedLogLikelihood - currentLogLikelihood < this->logLikelihoodConvergenceThreshold) {
//				if (abs(updatedLogLikelihood - currentLogLikelihood) < this->logLikelihoodConvergenceThreshold) {
						convergenceNotReached = 0;
				} else {
					currentLogLikelihood = updatedLogLikelihood;	
				}								
			}
		}		
	}
	this->ComputeLogLikelihood();
//	this->ComputeLogLikelihoodUsingDouble();
}

void rootedPhylogeny_tree::OptimizeModelParametersUsingMultiThreadedNelderMead() {
	bool verbose = 1;
	int timesRestarted = 0;
	if (verbose) {
		cout << "Optimizing model parameters" << endl;
	}	
	double currentLogLikelihood = 0;
	double updatedLogLikelihood = 0;
	bool convergenceNotReached = 1;	
	bool restart = 1;
	vector <double> currentParameters;
	double x[11];
	for (int i = 0; i < 11; i++) {
		x[i] = 0;
	}
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++) {
		if (this->root->rateCategory == rateCat) {
			rootedPhylogeny_vertex * c_l;
			rootedPhylogeny_vertex * c_r;
			c_l = (*this->vertexMap)[this->root->children_id[0]];
			c_r = (*this->vertexMap)[this->root->children_id[1]];
			if (c_l->rateCategory != rateCat and c_r->rateCategory != rateCat) {
				// Construct a modified version of Nelder Mead for the following case
				// no. of vertices that are in the same rate cat as the root equals one						
				this->OptimizeModelParametersForRootProbUsingNelderMead();
			} else {
				auto time_A = std::chrono::high_resolution_clock::now();
				convergenceNotReached = 1;
				this->ComputeLogLikelihood();
				currentLogLikelihood = this->logLikelihood;
				if (verbose) {				
					cout << "Current loglikelihood is " << this->logLikelihood << endl;
				}
				while (convergenceNotReached) {
					restart = 1;
					timesRestarted = 0;
					// current parameters
					while (restart) {
						this->ComputeLogLikelihood();
						currentLogLikelihood = this->logLikelihood;						
						this->OptimizeParametersForRateMatrixUsingNelderMead(rateCat);
						this->ComputeLogLikelihood();
						updatedLogLikelihood = this->logLikelihood;
						if (timesRestarted < 5 and updatedLogLikelihood < currentLogLikelihood and abs(updatedLogLikelihood - currentLogLikelihood) > this->logLikelihoodConvergenceThreshold) {
							restart = 1;
							timesRestarted += 1;
							cout << "Restarting" << endl;
							cout << "Current loglikelihood is " << currentLogLikelihood << endl;
							cout << "Updated loglikelihood is " << updatedLogLikelihood << endl;
							// reset parameters
							for (int par = 0; par < 11; par ++){
								x[par] = currentParameters[par];
							}
							this->SetParametersForRateMatrixForNelderMead(x, rateCat);
						} else {							
							restart = 0;
							currentParameters = this->GetParametersForRateCat(rateCat);
						}
					}
					this->ComputeLogLikelihood();
					currentLogLikelihood = this->logLikelihood;
					auto time_B = std::chrono::high_resolution_clock::now();
					if (verbose) {					
						cout << "CPU time used for optimizing rate matrix is " << chrono::duration<double>(time_B-time_A).count() << " second(s)\n";	
						cout << "Optimizing edge lengths" << endl;
					}					
//					this->OptimizeEdgeLengthsForRateCat(rateCat);
//					this->OptimizeEdgeLengthsForRateCatUsingBrent(rateCat);
					this->OptimizeEdgeLengthsForRateCatUsingNRWithScalers(rateCat);
					auto time_C = std::chrono::high_resolution_clock::now();
					if (verbose) {					
						cout << "CPU time used for optimizing edge lengths is " << chrono::duration<double>(time_C-time_B).count() << " second(s)\n";			
					}
					this->ComputeLogLikelihood();
					updatedLogLikelihood = this->logLikelihood;
					if (updatedLogLikelihood - currentLogLikelihood < this->logLikelihoodConvergenceThreshold) {
//					if (abs(updatedLogLikelihood - currentLogLikelihood) < this->logLikelihoodConvergenceThreshold) {
						convergenceNotReached = 0;						
					}					
					currentLogLikelihood = updatedLogLikelihood;
					if (verbose){
						cout << "Updated loglikelihood is " << this->logLikelihood << endl;
					}					
				}			
//				cout << "Edge length optimization complete" << endl;				
//				cout << "updated loglikelihood 2 is " << this->logLikelihood << endl;
			}
		} else {
			convergenceNotReached = 1;
			this->ComputeLogLikelihood();				
			currentLogLikelihood = this->logLikelihood;
			if (verbose) {
				cout << "Current loglikelihood is " << this->logLikelihood << endl;
			}
			while (convergenceNotReached) {
				restart = 1;
				timesRestarted = 0;
				while (restart) {
					this->ComputeLogLikelihood();
					currentLogLikelihood = this->logLikelihood;
					currentParameters = this->GetParametersForRateCat(rateCat);
					this->OptimizeParametersForRateMatrixUsingNelderMead(rateCat);
					this->ComputeLogLikelihood();
					updatedLogLikelihood = this->logLikelihood;
					if (timesRestarted < 5 and updatedLogLikelihood < currentLogLikelihood and abs(updatedLogLikelihood - currentLogLikelihood) > this->logLikelihoodConvergenceThreshold) {
						restart = 1;
						timesRestarted += 1;
						cout << "Restarting" << endl;
						cout << "Current loglikelihood is " << currentLogLikelihood << endl;
						cout << "Updated loglikelihood is " << updatedLogLikelihood << endl;
						// reset parameters
						for (int par = 0; par < 11; par ++){
							x[par] = currentParameters[par];
						}
						this->SetParametersForRateMatrixForNelderMead(x, rateCat);
					} else {
						restart = 0;
					}
				}
				this->ComputeLogLikelihood();
//				this->ComputeLogLikelihoodUsingDouble();
				currentLogLikelihood = this->logLikelihood;
//				this->OptimizeEdgeLengthsForRateCat(rateCat);	
//				this->OptimizeEdgeLengthsForRateCatUsingBrent(rateCat);							
				this->OptimizeEdgeLengthsForRateCatUsingNRWithScalers(rateCat);
				this->ComputeLogLikelihood();
//				this->ComputeLogLikelihoodUsingDouble();
				updatedLogLikelihood = this->logLikelihood;
				if (verbose) {
					cout << "Updated loglikelihood is " << this->logLikelihood << endl;
				}
				if (updatedLogLikelihood - currentLogLikelihood < this->logLikelihoodConvergenceThreshold) {
//				if (abs(updatedLogLikelihood - currentLogLikelihood) < this->logLikelihoodConvergenceThreshold) {
						convergenceNotReached = 0;
				} else {
					currentLogLikelihood = updatedLogLikelihood;	
				}								
			}
		}		
	}
	this->ComputeLogLikelihood();
//	this->ComputeLogLikelihoodUsingDouble();
}


void rootedPhylogeny_tree::OptimizeModelParametersForRootProbUsingNewtonRaphson(){
    // double tol = pow(10,-8);	
	// Compute Jacobian for root prob parameters
	// Compute Hessian for root prob parameters
	// Compute search direction
	// Select step size using Brent's line search // impose bounds using penalty
	// Iterate till convergence
}

void rootedPhylogeny_tree::OptimizeModelParametersForRootProbUsingNelderMead(){
	int i;
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 3;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];
	
	// Initial estimate of parameters
	for (i = 0; i < n; i ++){
		start[i] = this->rootProbability[i];
	}
		
	reqmin = 1.0E-08;
	double stepSize = pow(10,-2);
    for (i = 0; i < n; i++){
		step[i] = stepSize;
	}
	
	konvge = 10;
	kcount = 500;
	
	ynewlo = this->GetNegLogLikelihoodForParametersOfRootProb( start );
		
	this->NelderMeadForOptimizingParametersOfRootProb( n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );
}

void rootedPhylogeny_tree::OptimizeParametersForRateMatrixUsingBFGS(int rateCat){	
	this->rateCategoryForOptimization = rateCat;
	MatrixXd Jacobian_current;
	MatrixXd Jacobian_updated;
	MatrixXd InverseHessian_current = ArrayXXd::Zero(11,11);
	MatrixXd IdentityMatrix = ArrayXXd::Zero(11,11);
	MatrixXd InverseHessian_updated;
	MatrixXd B_current;
	MatrixXd B_updated;
	MatrixXd s;
	MatrixXd y;
	MatrixXd rho = ArrayXXd::Zero(11,11);;
	MatrixXd TempMatrix_l = ArrayXXd::Zero(11,11);
	MatrixXd TempMatrix_r = ArrayXXd::Zero(11,11);
	double currentLogLikelihood;
	double updatedLogLikelihood;
	B_current = (*this->parametersPerRateCategory)[rateCat];
	bool convergenceNotReached = 1;
	double stepSize;
//	double normOfJacobian;
	int iter = 0;
//	int maxIter = 100;
	double scalingFactorforInitialInverseHessian;
	for (int par_1 = 0; par_1 < 11; par_1 ++) {
		InverseHessian_current(par_1, par_1) = 1.0;
		IdentityMatrix(par_1, par_1) = 1.0;
	}
	// Compute gradient
	Jacobian_current = this->GetJacobianForRateCategory(rateCat);
	scalingFactorforInitialInverseHessian = abs(Jacobian_current(0,0)) + abs(Jacobian_current(1,0)) + abs(Jacobian_current(2,0));
	scalingFactorforInitialInverseHessian /= 3.0;	
	InverseHessian_current = InverseHessian_current / scalingFactorforInitialInverseHessian;
//	normOfJacobian = this->GetNormOfJacobian(Jacobian_current);	
	this->ComputeLogLikelihoodUsingLongDouble();
	currentLogLikelihood = this->logLikelihood;
	
	while (convergenceNotReached) {
		currentLogLikelihood = this->logLikelihood;
		iter += 1;
//		cout << "Computing search direction " << endl;
		this->searchDirection = -1 * (InverseHessian_current * Jacobian_current);	
//		cout << "Search direction is " << endl;
//		cout << this->searchDirection << endl;
		// Select step size
		// Impose bounds using a penalty
//		cout << "Computing step size " << endl;
		stepSize = this->GetOptimalStepSizeForBFGS();	
//		cout << "Step size is " << endl;
//		cout << stepSize << endl;
//		// Update parameters	
//		cout << "Updating parameters " << endl;
		B_updated = B_current + stepSize * this->searchDirection;
//		cout << "Updated parameters are " << endl;
//		cout << B_updated << endl;
		this->SetParametersForRateMatrixForBFGS(B_updated, rateCat);
		Jacobian_updated = this->GetJacobianForRateCategory(rateCat);
		y = Jacobian_updated - Jacobian_current;
		s = B_updated - B_current;
//		cout << "y" << endl;
//		cout << y << endl;
//		cout << "s" << endl;
//		cout << s << endl;
		rho = y.transpose()*s;
		rho = rho.inverse();
//		cout << "rho is " << endl;
//		cout << rho << endl;
//		cout << s * y.transpose() << endl;
		TempMatrix_l = IdentityMatrix - rho(0,0) * ( s * y.transpose() );
		TempMatrix_r = IdentityMatrix - rho(0,0) * ( y * s.transpose() );
		InverseHessian_updated = TempMatrix_l * InverseHessian_current;
		InverseHessian_updated = InverseHessian_updated * TempMatrix_r;
		InverseHessian_updated += rho(0,0) * ( s * s.transpose() );
		B_current = B_updated;
		InverseHessian_current = InverseHessian_updated;
		Jacobian_current = Jacobian_updated;
//		normOfJacobian = this->GetNormOfJacobian(Jacobian_current);
		this->ComputeLogLikelihoodUsingLongDouble();		
		updatedLogLikelihood = this->logLikelihood;
		if (abs(updatedLogLikelihood - currentLogLikelihood) < this->logLikelihoodConvergenceThreshold){
			convergenceNotReached = 0;
		}
		currentLogLikelihood = updatedLogLikelihood;	
	}
}

void rootedPhylogeny_tree::OptimizeParametersForRateMatrixUsingNelderMeadForFullyLabeledTree(int rateCat){
	int i;
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 11;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];
	
	// Initial estimate of parameters
	for (i = 0; i < n; i ++){
		start[i] = (*this->parametersPerRateCategory)[rateCat](i,0);
	}
		
	reqmin = 1.0E-08;
	double stepSize = pow(10,-2);
    for (i = 0; i < n; i++){
		step[i] = stepSize;
	}
	
	konvge = 10;
	kcount = 500;
	ynewlo = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree(start, rateCat);	
	
	this->NelderMeadForOptimizingRateParametersForRateCatForFullyLabeledTree( rateCat, n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );
}

void rootedPhylogeny_tree::OptimizeParametersForRateMatrixUsingNelderMead(int rateCat){
	int i;
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 11;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];
	
	// Initial estimate of parameters
	for (i = 0; i < n; i ++){
		start[i] = (*this->parametersPerRateCategory)[rateCat](i,0);
	}
		
	reqmin = 1.0E-08;
	double stepSize = pow(10,-2);
    for (i = 0; i < n; i++){
		step[i] = stepSize;
	}
	
	konvge = 10;
	kcount = 500;
	
	ynewlo = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( start , rateCat);
	
	this->NelderMeadForOptimizingRateParametersForRateCat( rateCat, n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );
	
	delete [] start;
	delete [] step;
	delete [] xmin;
}

void rootedPhylogeny_tree::NelderMeadForOptimizingParametersOfRootProb( int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }    
	y[n] = this->GetNegLogLikelihoodForParametersOfRootProb( start );
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = this->GetNegLogLikelihoodForParametersOfRootProb( start );
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = this->GetNegLogLikelihoodForParametersOfRootProb( pstar );
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
		y2star = this->GetNegLogLikelihoodForParametersOfRootProb( p2star );
//        y2star = this->powell ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForParametersOfRootProb( p2star );
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
			  y[j] = this->GetNegLogLikelihoodForParametersOfRootProb( xmin );
//              y[j] = this->powell ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForParametersOfRootProb( p2star );
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
	  z = this->GetNegLogLikelihoodForParametersOfRootProb( xmin );
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
	  z = this->GetNegLogLikelihoodForParametersOfRootProb( xmin );
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}

void rootedPhylogeny_tree::NelderMeadForOptimizingRateParametersForRateCat( int rateCat, int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    y[n] = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( start , rateCat);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( start , rateCat);
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( pstar , rateCat);;
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
		y2star = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( p2star , rateCat);
//        y2star = this->powell ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( p2star , rateCat);
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
			  y[j] = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( xmin , rateCat);
//              y[j] = this->powell ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( p2star , rateCat);
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
	  z = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( xmin , rateCat);
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
	  z = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMead( xmin , rateCat);
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}

void rootedPhylogeny_tree::NelderMeadForOptimizingRateParametersForRateCatForFullyLabeledTree( int rateCat, int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    y[n] = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( start , rateCat);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( start , rateCat);
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( pstar , rateCat);;
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
		y2star = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( p2star , rateCat);
//        y2star = this->powell ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( p2star , rateCat);
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
			  y[j] = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( xmin , rateCat);
//              y[j] = this->powell ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( p2star , rateCat);
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
	  z = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( xmin , rateCat);
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
	  z = this->GetNegLogLikelihoodForRateParametersForRateCatForNelderMeadForFullyLabeledTree( xmin , rateCat);
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}

void rootedPhylogeny_tree::NelderMeadForPowell( int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    y[n] = this->powell( start );
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = this->powell ( start );
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = this->powell ( pstar );
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        y2star = this->powell ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
              y[j] = this->powell ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}

//void rootedPhylogeny_tree::testNelderMead_1(){
//	int i;
//	int icount;
//	int ifault;
//	int kcount;
//	int konvge;
//	int n;
//	int numres;
//	double reqmin;
//	double *start;
//	double *step;
//	double *xmin;
//	double ynewlo;
//	
//	n = 4;
//
//	start = new double[n];
//	step = new double[n];
//	xmin = new double[n];
//
//	cout << "\n";
//	cout << "TEST02\n";
//	cout << "  Apply NELMIN to POWELL quartic function.\n";
//
//	start[0] =   3.0;
//	start[1] = - 1.0;
//	start[2] =   0.0;
//	start[3] =   1.0;
//
//	reqmin = 1.0E-08;
//
//	step[0] = 1.0;
//	step[1] = 1.0;
//	step[2] = 1.0;
//	step[3] = 1.0;
//
//	konvge = 10;
//	kcount = 500;
//
//	cout << "\n";
//	cout << "  Starting point X:\n";
//	cout << "\n";
//	
//	for ( i = 0; i < n; i++ )
//	{
//	cout << "  " << setw(14) << start[i] << "\n";
//	}
//	
//	ynewlo = powell ( start );
//
//	cout << "\n";
//	cout << "  F(X) = " << ynewlo << "\n";
//
//	nelmin ( powell, n, start, xmin, &ynewlo, reqmin, step,
//	konvge, kcount, &icount, &numres, &ifault );
//
//	cout << "\n";
//	cout << "  Return code IFAULT = " << ifault << "\n";
//	cout << "\n";
//	cout << "  Estimate of minimizing value X*:\n";
//	cout << "\n";
//	for ( i = 0; i < n; i++ )
//	{
//	cout << "  " << setw(14) << xmin[i] << "\n";
//	}
//
//	cout << "\n";
//	cout << "  F(X*) = " << ynewlo << "\n";
//
//	cout << "\n";
//	cout << "  Number of iterations = " << icount << "\n";
//	cout << "  Number of restarts =   " << numres << "\n";
//
//	delete [] start;
//	delete [] step;
//	delete [] xmin;
//	
//	return ;
//}

void rootedPhylogeny_tree::testNelderMead_2(){
	int i;
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 4;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];

	cout << "\n";
	cout << "TEST02\n";
	cout << "  Apply NELMIN to POWELL quartic function.\n";

	start[0] =   3.0;
	start[1] = - 1.0;
	start[2] =   0.0;
	start[3] =   1.0;

	reqmin = 1.0E-08;

	step[0] = 1.0;
	step[1] = 1.0;
	step[2] = 1.0;
	step[3] = 1.0;

	konvge = 10;
	kcount = 500;

	cout << "\n";
	cout << "  Starting point X:\n";
	cout << "\n";
	
	for ( i = 0; i < n; i++ )
	{
	cout << "  " << setw(14) << start[i] << "\n";
	}
	
	ynewlo = this->powell ( start );

	cout << "\n";
	cout << "  F(X) = " << ynewlo << "\n";
	
	this->NelderMeadForPowell( n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );

	cout << "\n";
	cout << "  Return code IFAULT = " << ifault << "\n";
	cout << "\n";
	cout << "  Estimate of minimizing value X*:\n";
	cout << "\n";
	
	for ( i = 0; i < n; i++ )
	{
	cout << "  " << setw(14) << xmin[i] << "\n";
	}

	cout << "\n";
	cout << "  F(X*) = " << ynewlo << "\n";

	cout << "\n";
	cout << "  Number of iterations = " << icount << "\n";
	cout << "  Number of restarts =   " << numres << "\n";

	delete [] start;
	delete [] step;
	delete [] xmin;
	
	return ;
}


double rootedPhylogeny_tree::powell(double x[4]){
	double fx;
	double fx1;
	double fx2;
	double fx3;
	double fx4;

	fx1 = x[0] + 10.0 * x[1];
	fx2 = x[2] - x[3];
	fx3 = x[1] - 2.0 * x[2];
	fx4 = x[0] - x[3];

	fx =            fx1 * fx1
	 +  5.0 * fx2 * fx2
	 +            fx3 * fx3 * fx3 * fx3
	 + 10.0 * fx4 * fx4 * fx4 * fx4;;

	return fx;
}

double rootedPhylogeny_tree::GetBICForGenerallyLabeledTreeUsingLongDouble(){	
	Matrix4f Q;
//	float scalingFactor;		
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f Q_scaled; Matrix4f P;
//	float t;	
	// Compute log likelihood
	
	// Compute BIC
//	this->SetMinLengthOfEdges();
//	long double maxValueForConditionalLikelihood;
//	long double partialLikelihood;
//	long double siteLikelihood;
//	long double logSiteLikelihood_mpf_float_1000;
	
	map <int, array <long double, 4>> conditionalLikelihoodMap;
	array <long double, 4> conditionalLikelihood;
//	long double maxConditionalLikelihood = 1;
//	double logSiteLikelihood_double;
//	bool showConditionalLikelihood = 0;	
	for (int i = 0; i < 4; i++){
		if (this->rootProbability[i] < 0){
			cout << "Root prob is negative" << endl;
			cout << "Rate matrix is " << endl;
			cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
			break;
		}
//		cout << "root prob for nt " << i << " is " << this->rootProbability[i] << endl;
	}
	this->logLikelihood = 0;
	
	for (unsigned int site = 0; site < this->siteWeights.size(); site++) {
		conditionalLikelihoodMap.clear();
		for (pair<int,int> edge_ids : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edge_ids.first];
			c = (*this->vertexMap)[edge_ids.second];
			if (c->id < this->numberOfObservedSequences) {				
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 0;
				}				
				conditionalLikelihood[c->compressedSequence[site]] = 1;
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for unobserved vertex			
			if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()) {
				
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 1;
				}
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));
			}
		}
	}
//		conditionalLikelihoodMap.clear();
//		for (pair<int,int> edge_ids : *this->edgesForPostOrderTreeTraversal) {
//			p = (*this->vertexMap)[edge_ids.first];
//			c = (*this->vertexMap)[edge_ids.second];
////			cout << "parent name is " << p->name << "\t" << "child name is " << c->name << endl;
//			// Initialize conditional likelihood for child if child is a leaf
//			if (c->id < this->numberOfObservedSequences) {				
//				for (int dna = 0; dna < 4; dna++) {
//					conditionalLikelihood[dna] = 0;
//				}				
//				conditionalLikelihood[c->compressedSequence[site]] = 1;
////				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
//				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c->id,conditionalLikelihood));
//			}
//			// Initialize conditional likelihood for parent			
//			if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()) {				
//				for (int dna = 0; dna < 4; dna++) {
//					conditionalLikelihood[dna] = 1;
////					if (p->name == "h_8"){
////						cout << "conditional likelihood for h_8 for dna " << dna << " is " << conditionalLikelihood[dna] << endl;
////					}
//				}
//				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));
////				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));
//			}
//						
//			partialLikelihood = 0;
//			t = this->GetEdgeLength(p->id, c->id);
////			cout << "edge length for " << p->id << " - " << c->id << " is " << t << endl;
//			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];			
//			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
////			cout << "scaling Factor is " << scalingFactor << endl;
//			Q_scaled = Q*(t/scalingFactor);
//			P = Q_scaled.exp();			
////			if (p->name == "h_8"){
////				cout << "Transition matrix is " << endl;
////				cout << P << endl;
////			}
////			maxValueForConditionalLikelihood = 0;
//			for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
//				partialLikelihood = 0;
//				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
//					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
//				}
//				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;				
////				if (conditionalLikelihoodMap[p->id][dna_p] > maxValueForConditionalLikelihood){
////					maxValueForConditionalLikelihood = conditionalLikelihoodMap[p->id][dna_p];				
////				}	
////				if (p->name == "h_1377"){
////					cout << "child name is " << c->name << endl;
////					cout << "partialLikelihood for dna " << static_cast <unsigned> (dna_p) << " is " << partialLikelihood << endl;
////				}
////				cout << "partial likelihood for dna " << static_cast<unsigned> (dna_p) << " is ";
////				cout << partialLikelihood << endl;
//			}
////			if (this->smallestMaxValueForConditionalLikelihoodOfAncestors > maxValueForConditionalLikelihood){
////				this->smallestMaxValueForConditionalLikelihoodOfAncestors = maxValueForConditionalLikelihood;
////			}			
//			// Erase conditional likelihood for child
////			if(conditionalLikelihoodMap.find(c->id) != conditionalLikelihoodMap.end()) {								
////				conditionalLikelihoodMap.erase(c->id);				
////			} else {
////				cout << "Child not present in conditional likelihood map" << endl;
////			}			
//		}		
//		for (rootedPhylogeny_vertex * v: *this->verticesForPreOrderTreeTraversal) {
//			showConditionalLikelihood = 0;
//			maxConditionalLikelihood = 0;
//			for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
//				if (maxConditionalLikelihood < conditionalLikelihoodMap[v->id][dna_c]) {
//					maxConditionalLikelihood = conditionalLikelihoodMap[v->id][dna_c];
//				}
//			}
////			if (maxConditionalLikelihood < pow(10,-10)){
////				showConditionalLikelihood = 1;
////			}
//			if (showConditionalLikelihood) {
//				cout << (*this->vertexMap)[v->parent_id]->name << " - " << v->name << endl;
//			}			
//		}
//		if (conditionalLikelihoodMap.find(-1) == conditionalLikelihoodMap.end()) {
//			cout << "Conditional likelihood for root is not computed" << endl;
//		}		
//		siteLikelihood = 0;
//		for (unsigned char dna = 0; dna < 4; dna++) {			
//			siteLikelihood += this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
//		}
//		if (site == 0 or site == 10 or site == 20 or site == 30 or site == 40 or site == 50) {
////			cout << "Likelihood: sitelikelihood for site " << site << " is " << siteLikelihood << endl;
////			for (int i = 0; i < 4; i++){
////				cout << "Root probability for " << i << " is ";
////				cout << this->rootProbability[i] << endl;
////				cout << "Conditional likelihood for root for " << i << " is ";				
////				cout << conditionalLikelihoodMap[-1][i] << endl;
////				cout << "Conditional likelihood for left child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[0]]->name << "of root for " << i << " is ";				
////				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[0]][i] << endl;
////				cout << "Edge length for left child is " << this->GetEdgeLength(-1,(*this->vertexMap)[-1]->children_id[0]) << endl;				
////				cout << "Conditional likelihood for right child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[1]]->name << " of root for " << i << " is ";				
////				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[1]][i] << endl;
////				cout << "Edge length for right child is " << this->GetEdgeLength(-1,(*this->vertexMap)[-1]->children_id[1]) << endl;
////				cout << "Conditional likelihood for vertex 17 for " << i << " is ";
////				cout << conditionalLikelihoodMap[17][i] << endl;				
////				cout << "Rate matrix is" << endl;
////				cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
////				cout << "Number of rate categories is" << endl;
////				cout << this->numberOfRateCategories << endl;
////			}
//		}
//		if (siteLikelihood == 0) {		
//			cout << "Problem with computing siteLikelihood for site " << site << endl;			
//			for (int i = 0; i < 4; i++){
////				cout << "Root probability for " << i << " is ";
////				cout << this->rootProbability[i] << endl;
////				cout << "Conditional likelihood for root for " << i << " is ";				
////				cout << conditionalLikelihoodMap[-1][i] << endl;
////				cout << "Conditional likelihood for left child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[0]]->name << "of root for " << i << " is ";				
////				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[0]][i] << endl;
////				cout << "Conditional likelihood for right child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[1]]->name << " of root for " << i << " is ";				
////				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[1]][i] << endl;
////				cout << "Conditional likelihood for vertex 17 for " << i << " is ";
////				cout << conditionalLikelihoodMap[17][i] << endl;				
////				cout << "Rate matrix is" << endl;
////				cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
////				cout << "Number of rate categories is" << endl;
////				cout << this->numberOfRateCategories << endl;
//			}
////			cout << "Edge length for vertex 17 is " << this->GetEdgeLength(17,(*this->vertexMap)[17]->parent_id) << endl;			
////			cout << "root probability for 0 is " << this->rootProbability[0] << endl;
//			break;
//		}
////		cout << "smallest maximum value of conditional likelihood is " << setprecision(100) << this->smallestMaxValueForConditionalLikelihoodOfAncestors << endl;
////		logSiteLikelihood_mpf_float_1000 = log(siteLikelihood);
////		logSiteLikelihood_double = logSiteLikelihood_mpf_float_1000.convert_to<double>();
//		this->logLikelihood += log(siteLikelihood) * this->siteWeights[site];
//		// Erase conditional likelihood for root
//		conditionalLikelihoodMap.erase(-1);
	double BIC_to_return = 0.0;
	return (BIC_to_return);
}

void rootedPhylogeny_tree::OptimizeEdgeLengthsForRateCatForFullyLabeledTree(int rateCat){
	rootedPhylogeny_vertex * p;	
	rootedPhylogeny_vertex * c;	
	// Iterate over edges for rate cat
	Matrix4f Q; Matrix4f Q_norm; Matrix4f Q_scaled; Matrix4f P;	
	Matrix4f P_first_der; Matrix4f P_second_der;
	double firstDerivative;
	double secondDerivative;
	double firstDerivative_perSite;
	double secondDerivative_perSite;
	float scalingFactor;
	unsigned char dna_p; unsigned char dna_c;
	float t;
	bool convergenceNotReached;
	float stepSize;
	float proposedChange_float;
	for (pair <int, rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap) {
		c = idPtrPair.second;
		if (c->id != -1){
			if (c->rateCategory == rateCat) {				
				p = (*this->vertexMap)[c->parent_id];
				scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
				Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
				Q_norm = Q/scalingFactor;
				// Iterate till convergence
				convergenceNotReached = 1;
				while (convergenceNotReached) {
					t = this->GetEdgeLength(p->id,c->id);
					Q_scaled = Q_norm * t;				
					P = Q_scaled.exp();
					P_first_der = Q_norm * P;
					P_second_der = Q_norm * P_first_der;
					// Compute first and second derivatives
					firstDerivative = 0;
					secondDerivative = 0;
					for (int site = 0; site < (int) this->siteWeights.size(); site++) {
						dna_p = p->compressedSequence[site];
						dna_c = c->compressedSequence[site];
						firstDerivative_perSite = P_first_der(dna_p,dna_c)/P(dna_p,dna_c);
						firstDerivative += firstDerivative_perSite * this->siteWeights[site];
						secondDerivative_perSite = P_second_der(dna_p,dna_c)/P(dna_p,dna_c) - pow(firstDerivative_perSite,2);
						secondDerivative += secondDerivative_perSite * this->siteWeights[site];
					}								
					stepSize = 1.0;
					proposedChange_float = firstDerivative/secondDerivative;					
					if (abs(proposedChange_float) < pow(10,-4) or t < pow(10,-5)) {
						convergenceNotReached = 0;		
					} else {
						if (t - proposedChange_float < 0.0) {
//							cout << "Case 1: proposed edge length is negative" << endl;
							while (t - stepSize * (proposedChange_float) < 0.0) {
								stepSize /= 2.0;
							}					
							t -= stepSize * (proposedChange_float);
						} else if (t - proposedChange_float > 1.0){
//							cout << "Case 2: proposed edge length is greater than one" << endl;														
//							cout << "vertex id is ";
//							cout << c->id << endl;
//							cout << "Edge length before updating is " << t << "\n";
//							cout << "Proposed change is " << -1 * proposedChange_float << "\n";
							t = pow(10,-5);							
						} else {							
							t -= proposedChange_float;
						}
//							cout << "Edge length after updating is " << t << "\n";
						this->SetEdgeLength(p->id,c->id,t);
					}
				}				
			}
		}
	}
}

void rootedPhylogeny_tree::OptimizeEdgeLengthUsingBrent(rootedPhylogeny_vertex * c) {
	rootedPhylogeny_vertex * p;
	p = (*this->vertexMap)[c->parent_id];
	float t = this->GetEdgeLength(p->id,c->id);
	float limit_lower = pow(10,-7);
	float limit_upper = 1.0;
	float thresh = pow(10,-4);
	this->OptimizeEdgeLengthUsingBrentsLineSearch(c,limit_lower,limit_upper,thresh,t);
	this->SetEdgeLength(p->id,c->id,t);
}

void rootedPhylogeny_tree::OptimizeEdgeLengthsForRateCatUsingBrent(int rateCat) {		
	vector <rootedPhylogeny_vertex*> verticesForRateCat;
	bool convergenceNotReached = 1;
	float currentLogLikelihood;
	float updatedLogLikelihood;
	for (rootedPhylogeny_vertex* b: *this->verticesForPreOrderTreeTraversal){
		if (b->id != -1 and b->rateCategory == rateCat){
			verticesForRateCat.push_back(b);
		}
	}
	this->ComputeLogLikelihood();
	currentLogLikelihood = this->logLikelihood;
	while (convergenceNotReached) {		
		for (rootedPhylogeny_vertex * child : verticesForRateCat) {		
			this->OptimizeEdgeLengthUsingBrent(child);
		}
		this->ComputeLogLikelihood();
		updatedLogLikelihood = this->logLikelihood;
		cout << "updatedLogLikelihood - currentLogLikelihood is " << updatedLogLikelihood - currentLogLikelihood << endl;
		if (updatedLogLikelihood - currentLogLikelihood < this->logLikelihoodConvergenceThreshold) {
			convergenceNotReached = 0;		
		} else {
			currentLogLikelihood = updatedLogLikelihood;
		}		
	}
}

void rootedPhylogeny_tree::OptimizeEdgeLengthsForRateCat(int rateCat){	
	map <int, array <long double, 4>> conditionalLikelihoodMap;
	array <long double, 4> conditionalLikelihood;
	map <int, array <long double, 4>> firstDerivativeOfConditionalLikelihoodMap;
	map <int, array <long double, 4>> secondDerivativeOfConditionalLikelihoodMap;
	array <long double, 4> firstDerivativeOfConditionalLikelihood;
	array <long double, 4> secondDerivativeOfConditionalLikelihood;
	long double firstDerivativeOfpartialLikelihood;
	long double secondDerivativeOfpartialLikelihood;
	long double partialLikelihood;
	long double siteLikelihood;
	long double firstDerivativeOfSiteLikelihood;
	long double secondDerivativeOfSiteLikelihood;
	long double firstDerivativeOfLogLikelihood;
	long double secondDerivativeOfLogLikelihood;
	float proposedChange_float;
	this->logLikelihood = 0;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * a;
//	rootedPhylogeny_vertex * b;
//	rootedPhylogeny_vertex * u;
	rootedPhylogeny_vertex * v;
	rootedPhylogeny_vertex * w;
//	rootedPhylogeny_vertex * parentOfA;
	float stepSize;
	Matrix4f Q; Matrix4f Q_norm; Matrix4f Q_scaled; float t; Matrix4f P;
	Matrix4f Q_ab; Matrix4f Q_ab_scaled; Matrix4f Q_ab_norm; float t_ab;
	Matrix4f Q_ac; Matrix4f Q_ac_scaled; Matrix4f Q_ac_norm; float t_ac;
	Matrix4f Q_uv; Matrix4f Q_uv_scaled; Matrix4f Q_uv_norm; float t_uv;
	Matrix4f Q_uw; Matrix4f Q_uw_scaled; Matrix4f Q_uw_norm; float t_uw;
	Matrix4f P_uv; Matrix4f P_uw;
	Matrix4f P_ab; Matrix4f P_ac; Matrix4f P_ab_der_1; Matrix4f P_ab_der_2;
//	float t_ab_new; float t_ab_old; 
	vector <rootedPhylogeny_vertex*> verticesForRateCat;	
	for (rootedPhylogeny_vertex* b: *this->verticesForPreOrderTreeTraversal){
		if (b->id != -1 and b->rateCategory == rateCat){
			verticesForRateCat.push_back(b);
		}
	}
	bool child_in_path;
	vector <rootedPhylogeny_vertex *> verticesInPathFromRootToANotIncA;
	// Iterate over vertices
	bool convergenceNotReached = 1;
	for (rootedPhylogeny_vertex * b : verticesForRateCat){		
		convergenceNotReached = 1;
		a = (*this->vertexMap)[b->parent_id];
//		cout << "optimizing length for edge " << a->name << " - " << b->name << endl;
//		cout << "current edge length is " << this->GetEdgeLength(a->id,b->id) << endl;
		// Length of edge (a, b) will be optimized
		verticesInPathFromRootToANotIncA.clear();
		while (a->id != -1 and a->parent_id != -1) {
			a = (*this->vertexMap)[a->parent_id];
			verticesInPathFromRootToANotIncA.push_back(a);
		}
		a = (*this->vertexMap)[b->parent_id];
		if (a->id != -1) {
			verticesInPathFromRootToANotIncA.push_back(this->root);
		}
		// Iterate till convergence
		while (convergenceNotReached){
			// Iterate over sites
			firstDerivativeOfLogLikelihood = 0;
			secondDerivativeOfLogLikelihood = 0;
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){									
				// Compute conditional likelihoods
				conditionalLikelihoodMap.clear();			
				for (pair <int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){
					p = (*this->vertexMap)[edgeIds.first];
					c = (*this->vertexMap)[edgeIds.second];
					Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
					t = this->GetEdgeLength(p->id,c->id);
					Q_norm = Q/(*this->scalingFactorForRateCategory)[c->rateCategory];
					Q_scaled = Q_norm * t;
					P = Q_scaled.exp();
					// Initialize conditional likelihood for leaves
					if (c->children_id.size()==0) {
						for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {
							conditionalLikelihood[dna_c] = 0.0;
						}
						conditionalLikelihood[c->compressedSequence[site]] = 1.0;
						conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c->id,conditionalLikelihood));
					}
					// Initialize conditional likelihood for ancestors
					if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
						conditionalLikelihood[dna_c] = 1.0;
						}				
						conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));					
					}		
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
						partialLikelihood = 0.0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++) {
							partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						}
						conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
					}
				}		
				// Compute first and second derivative of conditional likelihood
				firstDerivativeOfConditionalLikelihoodMap.clear();
				secondDerivativeOfConditionalLikelihoodMap.clear();
				// Reset a and c
				a = (*this->vertexMap)[b->parent_id];
				for (int child_id : a->children_id) {
					if (child_id != b->id){
						c = (*this->vertexMap)[child_id];
					}
				}
				// Initialize first and second derivative
				for (int dna = 0; dna < 4; dna ++) {
					firstDerivativeOfConditionalLikelihood[dna] = 1.0;
					secondDerivativeOfConditionalLikelihood[dna]  = 1.0;
				}
				
				for (int dna_p = 0; dna_p < 4; dna_p ++) {				
					// child b. length of edge (a, b) is being optimized.
					Q_ab = (*this->rateMatrixPerRateCategory)[b->rateCategory];
					Q_ab_norm = Q_ab/(*this->scalingFactorForRateCategory)[b->rateCategory];
					t_ab = this->GetEdgeLength(a->id, b->id);
					Q_ab_scaled = Q_ab_norm * t_ab;
					P_ab = Q_ab_scaled.exp();
					P_ab_der_1 = Q_ab_norm * P_ab;
					P_ab_der_2 = Q_ab_norm * P_ab_der_1;
//					cout << P_ab_der_2 << endl;
					firstDerivativeOfpartialLikelihood = 0;
					secondDerivativeOfpartialLikelihood = 0;
					for (int dna_c = 0; dna_c < 4; dna_c ++) {
						firstDerivativeOfpartialLikelihood += P_ab_der_1(dna_p,dna_c)*conditionalLikelihoodMap[b->id][dna_c];
						secondDerivativeOfpartialLikelihood += P_ab_der_2(dna_p,dna_c)*conditionalLikelihoodMap[b->id][dna_c];
					}
					firstDerivativeOfConditionalLikelihood[dna_p] *= firstDerivativeOfpartialLikelihood;				
					secondDerivativeOfConditionalLikelihood[dna_p] *= secondDerivativeOfpartialLikelihood;				
//					cout << "second der of partial likelihood is " << secondDerivativeOfpartialLikelihood << endl;
					// child c
					Q_ac = (*this->rateMatrixPerRateCategory)[c->rateCategory];
					Q_ac_norm = Q_ac/(*this->scalingFactorForRateCategory)[c->rateCategory];
					t_ac = this->GetEdgeLength(a->id, c->id);
					Q_ac_scaled = Q_ac_norm * t_ac;
					P_ac = Q_ac_scaled.exp();
					partialLikelihood = 0;
					for (int dna_c = 0; dna_c < 4; dna_c ++) {
						partialLikelihood += P_ac(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
					}
					firstDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
					secondDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;				
				}				
				firstDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<long double,4>>(a->id,firstDerivativeOfConditionalLikelihood));
				secondDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<long double,4>>(a->id,secondDerivativeOfConditionalLikelihood));
				// Compute the first derivative for ancestors of a
				for (rootedPhylogeny_vertex * u : verticesInPathFromRootToANotIncA){
					for (int dna = 0; dna < 4; dna ++) {
						firstDerivativeOfConditionalLikelihood[dna] = 1.0;
						secondDerivativeOfConditionalLikelihood[dna] = 1.0;
					}
					for (int dna_p = 0; dna_p < 4; dna_p++) {
						for (int child_id : u->children_id) {
							child_in_path = find(verticesInPathFromRootToANotIncA.begin(),verticesInPathFromRootToANotIncA.end(),(*this->vertexMap)[child_id]) != verticesInPathFromRootToANotIncA.end();
							if (child_id == a->id or child_in_path) {							
								v = (*this->vertexMap)[child_id];
								Q_uv = (*this->rateMatrixPerRateCategory)[v->rateCategory];
								Q_uv_norm = Q_uv/(*this->scalingFactorForRateCategory)[v->rateCategory];
								t_uv = this->GetEdgeLength(u->id, v->id);
								Q_uv_scaled = Q_uv_norm * t_uv;
								P_uv = Q_uv_scaled.exp();
								firstDerivativeOfpartialLikelihood = 0;
								secondDerivativeOfpartialLikelihood = 0;
								for (int dna_c = 0; dna_c < 4; dna_c ++) {
									firstDerivativeOfpartialLikelihood += P_uv(dna_p,dna_c)*firstDerivativeOfConditionalLikelihoodMap[v->id][dna_c];
									secondDerivativeOfpartialLikelihood += P_uv(dna_p,dna_c)*secondDerivativeOfConditionalLikelihoodMap[v->id][dna_c];
								}
							} else {							
								w = (*this->vertexMap)[child_id];
								Q_uw = (*this->rateMatrixPerRateCategory)[w->rateCategory];
								Q_uw_norm = Q_uw/(*this->scalingFactorForRateCategory)[w->rateCategory];
								t_uw = this->GetEdgeLength(u->id, w->id);
								Q_uw_scaled = Q_uw_norm * t_uw;
								P_uw = Q_uw_scaled.exp();
								partialLikelihood = 0;
								for (int dna_c = 0; dna_c < 4; dna_c ++) {
									partialLikelihood += P_uw(dna_p,dna_c)*conditionalLikelihoodMap[w->id][dna_c];
								}
							}
						}
						firstDerivativeOfConditionalLikelihood[dna_p] *= firstDerivativeOfpartialLikelihood;
						firstDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
						secondDerivativeOfConditionalLikelihood[dna_p] *= secondDerivativeOfpartialLikelihood;
						secondDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
					}
					firstDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<long double,4>>(u->id,firstDerivativeOfConditionalLikelihood));
					secondDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<long double,4>>(u->id,secondDerivativeOfConditionalLikelihood));
				}
				// compute site Likelihood
				siteLikelihood = 0;
				firstDerivativeOfSiteLikelihood = 0;
				secondDerivativeOfSiteLikelihood = 0;
				for (int dna = 0; dna < 4; dna ++){
					siteLikelihood += this->rootProbability[dna] * conditionalLikelihoodMap[-1][dna];
					firstDerivativeOfSiteLikelihood += this->rootProbability[dna] * firstDerivativeOfConditionalLikelihoodMap[-1][dna];
					secondDerivativeOfSiteLikelihood += this->rootProbability[dna] * secondDerivativeOfConditionalLikelihoodMap[-1][dna];
				}
//				cout << "Second der of site likelihood is " << secondDerivativeOfSiteLikelihood << endl;
				// compute first derivative of log likelihood
				firstDerivativeOfLogLikelihood += (firstDerivativeOfSiteLikelihood/siteLikelihood) * this->siteWeights[site];
				// compute second derivative of log likelihood
				secondDerivativeOfLogLikelihood += (secondDerivativeOfSiteLikelihood/siteLikelihood)  * this->siteWeights[site];
				secondDerivativeOfLogLikelihood -= (pow((firstDerivativeOfSiteLikelihood/siteLikelihood),2.0))  * this->siteWeights[site];			
			}
			t_ab = this->GetEdgeLength(a->id, b->id);
			stepSize = 1.0;
//			proposedChange_100 = firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood;
//			proposedChange_double = proposedChange_100.convert_to<double>();
			proposedChange_float = firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood;
//			cout << "Edge length before updating is " << t_ab << "\n";
//			cout << "first derivative is " << firstDerivativeOfLogLikelihood << endl;
//			cout << "second derivative is " << secondDerivativeOfLogLikelihood << endl;
			if (abs(proposedChange_float) < pow(10,-4) or t_ab < pow(10,-5)){
				convergenceNotReached = 0;
//				cout << "convergence reached" << endl;
			} else {
				if (t_ab - proposedChange_float < 0){
//					cout << "Case 1" << endl;
					while (t_ab - stepSize * (proposedChange_float) < 0){
						stepSize /= 2.0;
					}					
					t_ab -= stepSize * (proposedChange_float);
				} else {
//					cout << "Case 2" << endl;
					t_ab -= proposedChange_float;
				}
//				cout << "Edge length after updating is " << t_ab << "\n";
				this->SetEdgeLength(a->id,b->id,t_ab);
			}					
		}				
	}		
}


void rootedPhylogeny_tree::OptimizeEdgeLengthsForRateCatUsingNRWithScalers(int rateCat){	
	bool verbose = 0;
	map <rootedPhylogeny_vertex*, array <double, 4>> conditionalLikelihoodMap;
	array <double, 4> conditionalLikelihood;
	map <rootedPhylogeny_vertex*, array <double, 4>> firstDerivativeOfConditionalLikelihoodMap;
	map <rootedPhylogeny_vertex*, array <double, 4>> secondDerivativeOfConditionalLikelihoodMap;
	array <double, 4> firstDerivativeOfConditionalLikelihood;
	array <double, 4> secondDerivativeOfConditionalLikelihood;
	double firstDerivativeOfpartialLikelihood;
	double secondDerivativeOfpartialLikelihood;
	double partialLikelihood;
	double siteLikelihood;
	double firstDerivativeOfSiteLikelihood;
	double secondDerivativeOfSiteLikelihood;
	double firstDerivativeOfLogLikelihood;
	double secondDerivativeOfLogLikelihood;
	double firstDerivativeOfLogLikelihoodPerSite;
	double secondDerivativeOfLogLikelihoodPerSite;
	double proposedChange_float;
	double minDiffInLogScalingFactors_firstDer = -1 * pow(10,10);
	double minDiffInLogScalingFactors_secondDer = -1 * pow(10,10);
	double logLargestElement;
	float exponentForFirstDerivative;
	float exponentForSecondDerivative;
	this->logLikelihood = 0;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * a;
//	rootedPhylogeny_vertex * b;
//	rootedPhylogeny_vertex * u;
	rootedPhylogeny_vertex * v;
	rootedPhylogeny_vertex * w;
//	rootedPhylogeny_vertex * parentOfA;
	float stepSize;
	Matrix4f Q; Matrix4f Q_norm; Matrix4f Q_scaled; float t; Matrix4f P;
	Matrix4f Q_ab; Matrix4f Q_ab_scaled; Matrix4f Q_ab_norm; float t_ab;
	Matrix4f Q_ac; Matrix4f Q_ac_scaled; Matrix4f Q_ac_norm; float t_ac;
	Matrix4f Q_uv; Matrix4f Q_uv_scaled; Matrix4f Q_uv_norm; float t_uv;
	Matrix4f Q_uw; Matrix4f Q_uw_scaled; Matrix4f Q_uw_norm; float t_uw;
	Matrix4f P_uv; Matrix4f P_uw;
	Matrix4f P_ab; Matrix4f P_ac; Matrix4f P_ab_der_1; Matrix4f P_ab_der_2;
	vector <rootedPhylogeny_vertex*> verticesForRateCat;	
	for (rootedPhylogeny_vertex* b: *this->verticesForPreOrderTreeTraversal){
		if (b->id != -1 and b->rateCategory == rateCat){
			verticesForRateCat.push_back(b);
		}
	}
	bool child_in_path;
	vector <rootedPhylogeny_vertex *> verticesInPathFromRootToANotIncA;
	// Iterate over vertices
	bool convergenceNotReached = 1;
	for (rootedPhylogeny_vertex * b : verticesForRateCat){		
		convergenceNotReached = 1;
		a = (*this->vertexMap)[b->parent_id];
		if (verbose) {
			cout << "Optimizing length for edge " << a->name << " - " << b->name << endl;
			cout << "current edge length is " << this->GetEdgeLength(a->id,b->id) << endl;
		}
		// Length of edge (a, b) will be optimized
		verticesInPathFromRootToANotIncA.clear();
		while (a->id != -1 and a->parent_id != -1) {
			a = (*this->vertexMap)[a->parent_id];
			verticesInPathFromRootToANotIncA.push_back(a);
		}
		a = (*this->vertexMap)[b->parent_id];
		if (a->id != -1) {
			verticesInPathFromRootToANotIncA.push_back(this->root);
		}
		// Iterate till convergence
		while (convergenceNotReached){			
			// Iterate over sites
			firstDerivativeOfLogLikelihood = 0;
			secondDerivativeOfLogLikelihood = 0;
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){									
				this->ResetLogScalingFactors();			
				// Compute conditional likelihoods
				conditionalLikelihoodMap.clear();
				firstDerivativeOfConditionalLikelihoodMap.clear();
				secondDerivativeOfConditionalLikelihoodMap.clear();
				for (pair <int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){					
					p = (*this->vertexMap)[edgeIds.first];
					c = (*this->vertexMap)[edgeIds.second];
//					cout << "p: " << p->name << "\t" << isnan(p->logScalingFactors) << endl;
//					cout << "c: " << c->name << "\t" << isnan(c->logScalingFactors) << endl;
					Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
					t = this->GetEdgeLength(p->id,c->id);
					Q_norm = Q/(*this->scalingFactorForRateCategory)[c->rateCategory];
					Q_scaled = Q_norm * t;
					P = Q_scaled.exp();
					// Initialize conditional likelihood for leaves
					if (c->children_id.size()==0) {
						for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {
							conditionalLikelihood[dna_c] = 0.0;
						}
						conditionalLikelihood[c->compressedSequence[site]] = 1.0;
						conditionalLikelihoodMap.insert(pair<rootedPhylogeny_vertex*,array<double,4>>(c,conditionalLikelihood));
					}
					// Initialize conditional likelihood for ancestors
					if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
						conditionalLikelihood[dna_c] = 1.0;
						}				
						conditionalLikelihoodMap.insert(pair<rootedPhylogeny_vertex*,array<double,4>>(p,conditionalLikelihood));					
					}
					p->logScalingFactors += c->logScalingFactors;
//					cout << "p: " << p->name << "\t" << isnan(p->logScalingFactors) << endl;
//					cout << "c: " << c->name << "\t" << isnan(c->logScalingFactors) << endl;
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
						partialLikelihood = 0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++) {
							partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c][dna_c];
						}
						conditionalLikelihoodMap[p][dna_p] *= partialLikelihood;						
					}
//					if (this->GetAbsLargestElement(conditionalLikelihoodMap[p]) == 0) {
//						cout << "Largest element in conditional likelihood vector for " << p->name << "is zero " << endl;
//					}					
					tie (logLargestElement,conditionalLikelihoodMap[p]) = GetLogAbsLargestElementAndScaledArray(conditionalLikelihoodMap[p]);
//					cout << "Before adding log transformed largest element" << endl;
//					cout << "p: " << p->name << "\t" << isnan(p->logScalingFactors) << endl;
//					if (isnan(logLargestElement)) {
//						cout << "logLargestElement for " << p->name << " is nan" << endl;
//					}
					p->logScalingFactors += logLargestElement;
//					cout << "After adding log transformed largest element" << endl;
//					cout << "p: " << p->name << "\t" << isnan(p->logScalingFactors) << endl;
//					cout << "isnan(p->logScalingFactors)\t" << isnan(p->logScalingFactors) << endl;
				}
//				cout << "Computed conditional likelihoods" << endl;
				// Compute first and second derivative of conditional likelihood
				// Reset a and 
//				cout << "isnan(b->logScalingFactors)\t" << isnan(b->logScalingFactors) << endl;
				a = (*this->vertexMap)[b->parent_id];
				for (int child_id : a->children_id) {
					if (child_id != b->id){
						c = (*this->vertexMap)[child_id];
					}
				}
				// Initialize first and second derivative
				for (int dna = 0; dna < 4; dna ++) {
					firstDerivativeOfConditionalLikelihood[dna] = 1.0;
					secondDerivativeOfConditionalLikelihood[dna]  = 1.0;
				}
				
				for (int dna_p = 0; dna_p < 4; dna_p ++) {				
					// child b. length of edge (a, b) is being optimized.
					Q_ab = (*this->rateMatrixPerRateCategory)[b->rateCategory];
					Q_ab_norm = Q_ab/(*this->scalingFactorForRateCategory)[b->rateCategory];
					t_ab = this->GetEdgeLength(a->id, b->id);
					Q_ab_scaled = Q_ab_norm * t_ab;
					P_ab = Q_ab_scaled.exp();
					P_ab_der_1 = Q_ab_norm * P_ab;
					P_ab_der_2 = Q_ab_norm * P_ab_der_1;
//					cout << P_ab_der_2 << endl;
					firstDerivativeOfpartialLikelihood = 0;
					secondDerivativeOfpartialLikelihood = 0;
					for (int dna_c = 0; dna_c < 4; dna_c ++) {
						firstDerivativeOfpartialLikelihood += P_ab_der_1(dna_p,dna_c)*conditionalLikelihoodMap[b][dna_c];
						secondDerivativeOfpartialLikelihood += P_ab_der_2(dna_p,dna_c)*conditionalLikelihoodMap[b][dna_c];
					}
//					cout << "First derivative of partial likelihood is " << firstDerivativeOfpartialLikelihood << endl;
					firstDerivativeOfConditionalLikelihood[dna_p] *= firstDerivativeOfpartialLikelihood;				
					secondDerivativeOfConditionalLikelihood[dna_p] *= secondDerivativeOfpartialLikelihood;				
//					cout << "second der of partial likelihood is " << secondDerivativeOfpartialLikelihood << endl;
					// child c
					Q_ac = (*this->rateMatrixPerRateCategory)[c->rateCategory];
					Q_ac_norm = Q_ac/(*this->scalingFactorForRateCategory)[c->rateCategory];
					t_ac = this->GetEdgeLength(a->id, c->id);
					Q_ac_scaled = Q_ac_norm * t_ac;
					P_ac = Q_ac_scaled.exp();
					partialLikelihood = 0;
					for (int dna_c = 0; dna_c < 4; dna_c ++) {
						partialLikelihood += P_ac(dna_p,dna_c)*conditionalLikelihoodMap[c][dna_c];
					}
//					cout << "Partial likelihood is " << partialLikelihood << endl;
					firstDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
					secondDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;				
				}
//				if (this->GetAbsLargestElement(firstDerivativeOfConditionalLikelihood) == 0.0) {
//					cout << "Element with largest absolute value in first derivative of conditional likelihood vector for " << a->name << "is zero " << endl;
//				}
				tie(logLargestElement,firstDerivativeOfConditionalLikelihood) = GetLogAbsLargestElementAndScaledArray(firstDerivativeOfConditionalLikelihood);				
				firstDerivativeOfConditionalLikelihoodMap.insert(pair<rootedPhylogeny_vertex *,array<double,4>>(a,firstDerivativeOfConditionalLikelihood));
				a->logScalingFactors_firstDer += b->logScalingFactors + c->logScalingFactors;
				a->logScalingFactors_firstDer += logLargestElement;				
//				if (this->GetAbsLargestElement(secondDerivativeOfConditionalLikelihood) == 0) {
//					cout << "Element with largest absolute value in second derivative of conditional likelihood vector for " << a->name << "is zero " << endl;
//				}
				tie(logLargestElement,secondDerivativeOfConditionalLikelihood) = GetLogAbsLargestElementAndScaledArray(secondDerivativeOfConditionalLikelihood);				
				secondDerivativeOfConditionalLikelihoodMap.insert(pair<rootedPhylogeny_vertex *,array<double,4>>(a,secondDerivativeOfConditionalLikelihood));
				a->logScalingFactors_secondDer += b->logScalingFactors + c->logScalingFactors;
				a->logScalingFactors_secondDer += logLargestElement;
// 				Compute the first derivative for ancestors of a
//				cout << "Computing derivatives for vertices in path from root to A (not inc a)" << endl;
				for (rootedPhylogeny_vertex * u : verticesInPathFromRootToANotIncA){
//					cout << u->name << endl;
					for (int dna = 0; dna < 4; dna ++) {
						firstDerivativeOfConditionalLikelihood[dna] = 1.0;
						secondDerivativeOfConditionalLikelihood[dna] = 1.0;
					}
					for (int dna_p = 0; dna_p < 4; dna_p++) {
						for (int child_id : u->children_id) {
							child_in_path = find(verticesInPathFromRootToANotIncA.begin(),verticesInPathFromRootToANotIncA.end(),(*this->vertexMap)[child_id]) != verticesInPathFromRootToANotIncA.end();
							if (child_id == a->id or child_in_path) {							
								v = (*this->vertexMap)[child_id];
								Q_uv = (*this->rateMatrixPerRateCategory)[v->rateCategory];
								Q_uv_norm = Q_uv/(*this->scalingFactorForRateCategory)[v->rateCategory];
								t_uv = this->GetEdgeLength(u->id, v->id);
								Q_uv_scaled = Q_uv_norm * t_uv;
								P_uv = Q_uv_scaled.exp();
								firstDerivativeOfpartialLikelihood = 0;
								secondDerivativeOfpartialLikelihood = 0;
								for (int dna_c = 0; dna_c < 4; dna_c ++) {
									firstDerivativeOfpartialLikelihood += P_uv(dna_p,dna_c)*firstDerivativeOfConditionalLikelihoodMap[v][dna_c];
									secondDerivativeOfpartialLikelihood += P_uv(dna_p,dna_c)*secondDerivativeOfConditionalLikelihoodMap[v][dna_c];
								}
//								cout << "First derivative of partial likelihood is " << firstDerivativeOfpartialLikelihood << endl;
							} else {		
								w = (*this->vertexMap)[child_id];
								Q_uw = (*this->rateMatrixPerRateCategory)[w->rateCategory];
								Q_uw_norm = Q_uw/(*this->scalingFactorForRateCategory)[w->rateCategory];
								t_uw = this->GetEdgeLength(u->id, w->id);
								Q_uw_scaled = Q_uw_norm * t_uw;
								P_uw = Q_uw_scaled.exp();
								partialLikelihood = 0;
								for (int dna_c = 0; dna_c < 4; dna_c ++) {
									partialLikelihood += P_uw(dna_p,dna_c)*conditionalLikelihoodMap[w][dna_c];
								}
//								cout << "Partial likelihood is " << partialLikelihood << endl;
							}
						}
						firstDerivativeOfConditionalLikelihood[dna_p] *= firstDerivativeOfpartialLikelihood;
						firstDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
//						cout << "firstDerivativeOfConditionalLikelihood[dna_p] is " << firstDerivativeOfConditionalLikelihood[dna_p] << endl;
						secondDerivativeOfConditionalLikelihood[dna_p] *= secondDerivativeOfpartialLikelihood;
						secondDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
					}
//					if (this->GetAbsLargestElement(firstDerivativeOfConditionalLikelihood) == 0) {
//						cout << "Site is " << site << endl;
//						cout << "Element with largest absolute value in first derivative of conditional likelihood vector for " << u->name << " is zero " << endl;
//					}
					tie(logLargestElement,firstDerivativeOfConditionalLikelihood) = GetLogAbsLargestElementAndScaledArray(firstDerivativeOfConditionalLikelihood);
					firstDerivativeOfConditionalLikelihoodMap.insert(pair<rootedPhylogeny_vertex *,array<double,4>>(u,firstDerivativeOfConditionalLikelihood));
					u->logScalingFactors_firstDer += v->logScalingFactors_firstDer + w->logScalingFactors;
					u->logScalingFactors_firstDer += logLargestElement;
//					if (this->GetAbsLargestElement(secondDerivativeOfConditionalLikelihood) == 0) {
//						cout << "Site is " << site << endl;
//						cout << "Element with largest absolute value in second derivative of conditional likelihood vector for " << u->name << " is zero " << endl;
//					}
					tie(logLargestElement,secondDerivativeOfConditionalLikelihood) = GetLogAbsLargestElementAndScaledArray(secondDerivativeOfConditionalLikelihood);				
					secondDerivativeOfConditionalLikelihoodMap.insert(pair<rootedPhylogeny_vertex *,array<double,4>>(u,secondDerivativeOfConditionalLikelihood));
					u->logScalingFactors_secondDer += v->logScalingFactors_secondDer + w->logScalingFactors;
					u->logScalingFactors_secondDer += logLargestElement;
				}
				// compute site Likelihood
				siteLikelihood = 0;
				firstDerivativeOfSiteLikelihood = 0;
				secondDerivativeOfSiteLikelihood = 0;
				for (int dna = 0; dna < 4; dna ++){
					siteLikelihood += this->rootProbability[dna] * conditionalLikelihoodMap[this->root][dna];
					firstDerivativeOfSiteLikelihood += this->rootProbability[dna] * firstDerivativeOfConditionalLikelihoodMap[this->root][dna];
					secondDerivativeOfSiteLikelihood += this->rootProbability[dna] * secondDerivativeOfConditionalLikelihoodMap[this->root][dna];
				}
//				cout << "Second der of site likelihood is " << secondDerivativeOfSiteLikelihood << endl;
				// compute first derivative of log likelihood
//				assert(this->root->logScalingFactors_firstDer - this->root->logScalingFactors > logSmallestValueForDouble);
				exponentForFirstDerivative = this->root->logScalingFactors_firstDer - this->root->logScalingFactors;
				assert(exponentForFirstDerivative > logSmallestValueForDouble);
				assert(exponentForFirstDerivative < logLargestValueForDouble);

				firstDerivativeOfLogLikelihoodPerSite = (exp(this->root->logScalingFactors_firstDer - this->root->logScalingFactors) * (firstDerivativeOfSiteLikelihood/siteLikelihood));				
				assert (log(abs(firstDerivativeOfLogLikelihoodPerSite))+log(this->sequenceLength) < logLargestValueForDouble);				
				
				firstDerivativeOfLogLikelihood +=  firstDerivativeOfLogLikelihoodPerSite * this->siteWeights[site];

				exponentForSecondDerivative = this->root->logScalingFactors_secondDer - this->root->logScalingFactors;
				assert(exponentForSecondDerivative > logSmallestValueForDouble);
				assert(exponentForSecondDerivative < logLargestValueForDouble);
			
				secondDerivativeOfLogLikelihoodPerSite = (exp(this->root->logScalingFactors_secondDer - this->root->logScalingFactors) * (secondDerivativeOfSiteLikelihood/siteLikelihood));
				secondDerivativeOfLogLikelihoodPerSite -= pow(firstDerivativeOfLogLikelihoodPerSite,2.0);
				assert (log(abs(secondDerivativeOfLogLikelihoodPerSite))+log(this->sequenceLength) < logLargestValueForDouble);				

				secondDerivativeOfLogLikelihood += secondDerivativeOfLogLikelihoodPerSite * this->siteWeights[site];
			}
			t_ab = this->GetEdgeLength(a->id, b->id);
			stepSize = 1.0;
//			proposedChange_100 = firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood;
//			proposedChange_double = proposedChange_100.convert_to<double>();
			proposedChange_float = firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood;
			if (verbose) {
				cout << "Edge length before updating is " << t_ab << "\n";
			}
//			cout << "first derivative is " << firstDerivativeOfLogLikelihood << endl;
//			cout << "second derivative is " << secondDerivativeOfLogLikelihood << endl;
			if (abs(proposedChange_float) < pow(10,-4) or t_ab < pow(10,-5)){
				convergenceNotReached = 0;
//				cout << "convergence reached" << endl;
			} else {
				if (t_ab - proposedChange_float < 0){
//					cout << "Case 1" << endl;
					while (t_ab - stepSize * (proposedChange_float) < 0){
						stepSize /= 2.0;
					}					
					t_ab -= stepSize * (proposedChange_float);
				} else {
//					cout << "Case 2" << endl;					
					t_ab -= proposedChange_float;
				}
				if (verbose) {
					cout << "Edge length after updating is " << t_ab << "\n";
				}
				if (t_ab > 1.0) {
					t_ab = pow(10,-5);
				}
				this->SetEdgeLength(a->id,b->id,t_ab);
			}					
		}				
	}		
}
double rootedPhylogeny_tree::GetNormOfJacobian(MatrixXd JacobianForRateMatrix){
	double norm = 0;
	for (int par = 0; par < 11; par++){
		norm += pow(JacobianForRateMatrix(par,0),2);
	}
	return (sqrt(norm));
}

void rootedPhylogeny_tree::SetParametersForRateMatrixForBFGS(MatrixXd parameters, int rateCat){
	(*this->parametersPerRateCategory)[rateCat] = parameters;
	Matrix4f Q = this->GetRateMatrixForFreeParameters(parameters);
	(*this->rateMatrixPerRateCategory)[rateCat] = Q;
	(*this->scalingFactorForRateCategory)[rateCat] = 1.0;
}

double rootedPhylogeny_tree::GetOptimalStepSizeForBFGS() {
	double stepSize;
	double t = pow(10,-4);
	double lower_limit = 0.0;
	double upper_limit = 0.0;
	double inc = pow(10,-3);
	MatrixXd B_current = (*this->parametersPerRateCategory)[this->rateCategoryForOptimization];
	MatrixXd B_updated;	
	// Get lower limit on step size such that solution at lower limit is feasible
	B_updated = B_current + lower_limit * this->searchDirection;
	while (this->DoParametersEncodeAFeasibleRateMatrix(B_updated)) {
		lower_limit -= inc;
		B_updated = B_current + lower_limit * this->searchDirection;
	}
	lower_limit += inc;
	B_updated = B_current + lower_limit * this->searchDirection;
//	if (this->DoParametersEncodeAFeasibleRateMatrix(B_updated)){
//		cout << "feasible lower limit found " << endl;
//		cout << "feasible lower limit is " << lower_limit << endl;
//	}	
	// Get upper limit on step size such that solution at upper limit is feasible
	B_updated = B_current + upper_limit * this->searchDirection;
	while (this->DoParametersEncodeAFeasibleRateMatrix(B_updated)){
		upper_limit += inc;
		B_updated = B_current + upper_limit * this->searchDirection;
	}
	upper_limit -= inc;
//	B_updated = B_current + upper_limit * this->searchDirection;
//	if (this->DoParametersEncodeAFeasibleRateMatrix(B_updated)){
//		cout << "feasible upper limit found " << endl;
//		cout << "feasible upper limit is " << upper_limit << endl;
//	}
	double result;
	result = this->BrentLineSearchForSelectingStepSize(lower_limit, upper_limit, t, stepSize);
	return (stepSize);	
}

double rootedPhylogeny_tree::GetOptimalStepSizeUsingLineSearch(MatrixXd B_current){
	double t = pow(10,-4);
	double lower_limit = 0.0;
	double upper_limit = 100.0;
	// par 0 is pi_1, par 1 is pi_2, par 2 is pi_3
	for (int par = 0; par < 3; par ++){
		if (this->searchDirection(par,0) < 0){
			upper_limit = min(upper_limit,(-1.0 * B_current(par,0))/(this->searchDirection(par,0)));
		} else {
			upper_limit = min(upper_limit,(1.0 - B_current(par,0))/(this->searchDirection(par,0)));
		}
	}
	double directionForPi_4 = this->searchDirection(0,0) + this->searchDirection(1,0) + this->searchDirection(2,0);
	double sum_pi_1_pi_2_pi_3 = B_current(0,0) + B_current(1,0) + B_current(2,0);	
	// check for pi_4
	if (directionForPi_4 < 0){
		upper_limit = min(upper_limit, (-1.0 * sum_pi_1_pi_2_pi_3)/directionForPi_4);
	} else {
		upper_limit = min(upper_limit, (1.0 - sum_pi_1_pi_2_pi_3)/directionForPi_4);
	}
	for (int par = 3; par < 10; par++){
		if (this->searchDirection(par,0) < 0){
			upper_limit = min(upper_limit,(-1 * B_current(par,0))/this->searchDirection(par,0));
		}
	}
	cout << "upper limit is " << upper_limit << endl;
	double stepSize;
	double result;
	result = this->BrentLineSearchForSelectingStepSize(lower_limit, upper_limit, t, stepSize);
	if (result > pow(10,10)){
		cout << result << endl;
	}
	return (stepSize);	
	
}

void rootedPhylogeny_tree::OptimizeEdgeLengthUsingBrentsLineSearch(rootedPhylogeny_vertex* child, float a, float b, float t, float &x){
	float c;
	float d;
	float e;
	float eps;
	float fu;
	float fv;
	float fw;
	float fx;
	float m;
	float p;
	float q;
	float r;
	float sa;
	float sb;
	float t2;
	float tol;
	float u;
	float v;
	float w;
	
	//	const double value = 2.220446049250313E-016;
	//
	//  C is the square of the inverse of the golden ratio.
	//
	c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

	eps = sqrt ( 2.220446049250313E-016 );

	sa = a;
	sb = b;
	x = sa + c * ( b - a );
	w = x;
	v = w;
	e = 0.0;
	fx = this->GetNegLogLikelihoodForEdgeLength(child, x);
//	fx = this->GetNegLogLikelihoodForStepSize( x );
//	fx = this->SampleFunctionForMinimization( x );
	fw = fx;
	fv = fw;
	
	for ( ; ; )
  {
    m = 0.5 * ( sa + sb ) ;
    tol = eps * fabs ( x ) + t;
    t2 = 2.0 * tol;
//
//  Check the stopping criterion.
//
    if ( fabs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
    {
      break;
    }
//
//  Fit a parabola.
//
    r = 0.0;
    q = r;
    p = q;

    if ( tol < fabs ( e ) )
    {
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = - p;
      }
      q = fabs ( q );
      r = e;
      e = d;
    }

    if ( fabs ( p ) < fabs ( 0.5 * q * r ) &&
         q * ( sa - x ) < p &&
         p < q * ( sb - x ) )
    {
//
//  Take the parabolic interpolation step.
//
      d = p / q;
      u = x + d;
//
//  F must not be evaluated too close to A or B.
//
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
      {
        if ( x < m )
        {
          d = tol;
        }
        else
        {
          d = - tol;
        }
      }
    }
//
//  A golden-section step.
//
    else
    {
      if ( x < m )
      {
        e = sb - x;
      }
      else
      {
        e = sa - x;
      }
      d = c * e;
    }
//
//  F must not be evaluated too close to X.
//
    if ( tol <= fabs ( d ) )
    {
      u = x + d;
    }
    else if ( 0.0 < d )
    {
      u = x + tol;
    }
    else
    {
      u = x - tol;
    }
	fu = this->GetNegLogLikelihoodForEdgeLength(child, u);
//	fu = this->GetNegLogLikelihoodForStepSize( u );
//    fu = this->SampleFunctionForMinimization( u );
//
//  Update A, B, V, W, and X.
//
    if ( fu <= fx )
    {
      if ( u < x )
      {
        sb = x;
      }
      else
      {
        sa = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x )
      {
        sa = u;
      }
      else
      {
        sb = u;
      }

      if ( fu <= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }
  return ;
}


double rootedPhylogeny_tree::BrentLineSearchForSelectingStepSize(double a, double b, double t, double &x){
	double c;
	double d;
	double e;
	double eps;
	double fu;
	double fv;
	double fw;
	double fx;
	double m;
	double p;
	double q;
	double r;
	double sa;
	double sb;
	double t2;
	double tol;
	double u;
	double v;
	double w;
	
//	const double value = 2.220446049250313E-016;
	//
	//  C is the square of the inverse of the golden ratio.
	//
	c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

	eps = sqrt ( 2.220446049250313E-016 );

	sa = a;
	sb = b;
	x = sa + c * ( b - a );
	w = x;
	v = w;
	e = 0.0;
	fx = this->GetNegLogLikelihoodForStepSize( x );
//	fx = this->SampleFunctionForMinimization( x );
	fw = fx;
	fv = fw;
	
	for ( ; ; )
  {
    m = 0.5 * ( sa + sb ) ;
    tol = eps * fabs ( x ) + t;
    t2 = 2.0 * tol;
//
//  Check the stopping criterion.
//
    if ( fabs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
    {
      break;
    }
//
//  Fit a parabola.
//
    r = 0.0;
    q = r;
    p = q;

    if ( tol < fabs ( e ) )
    {
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = - p;
      }
      q = fabs ( q );
      r = e;
      e = d;
    }

    if ( fabs ( p ) < fabs ( 0.5 * q * r ) &&
         q * ( sa - x ) < p &&
         p < q * ( sb - x ) )
    {
//
//  Take the parabolic interpolation step.
//
      d = p / q;
      u = x + d;
//
//  F must not be evaluated too close to A or B.
//
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
      {
        if ( x < m )
        {
          d = tol;
        }
        else
        {
          d = - tol;
        }
      }
    }
//
//  A golden-section step.
//
    else
    {
      if ( x < m )
      {
        e = sb - x;
      }
      else
      {
        e = sa - x;
      }
      d = c * e;
    }
//
//  F must not be evaluated too close to X.
//
    if ( tol <= fabs ( d ) )
    {
      u = x + d;
    }
    else if ( 0.0 < d )
    {
      u = x + tol;
    }
    else
    {
      u = x - tol;
    }
	fu = this->GetNegLogLikelihoodForStepSize( u );
//    fu = this->SampleFunctionForMinimization( u );
//
//  Update A, B, V, W, and X.
//
    if ( fu <= fx )
    {
      if ( u < x )
      {
        sb = x;
      }
      else
      {
        sa = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x )
      {
        sa = u;
      }
      else
      {
        sb = u;
      }

      if ( fu <= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}

double rootedPhylogeny_tree::GetNegLogLikelihoodForStepSize(double stepSize){	
	int rateCat = this->rateCategoryForOptimization;
	MatrixXd B_current = (*this->parametersPerRateCategory)[rateCat];
	MatrixXd B_updated = B_current + stepSize*this->searchDirection;
	bool feasible = this->DoParametersEncodeAFeasibleRateMatrix(B_updated);
	double valueToReturn;
	if (feasible){
//		cout << "feasible solution found" << endl;
		this->SetParametersForRateMatrixForBFGS(B_updated, rateCat);
		this->ComputeLogLikelihoodUsingLongDouble();		
		valueToReturn = -1 * this->logLikelihood;
		this->SetParametersForRateMatrixForBFGS(B_current, rateCat);
	} else {
//		cout << "feasible solution not found" << endl;
		valueToReturn = pow(10,20);
	}	
	return (valueToReturn);
}


float rootedPhylogeny_tree::GetNegLogLikelihoodForEdgeLength(rootedPhylogeny_vertex * c, float t) {
	assert(c->id != c->parent_id);
	rootedPhylogeny_vertex * p = (*this->vertexMap)[c->parent_id];
	float t_current = this->GetEdgeLength(p->id,c->id);
	if (t < 0) {
		t *= -1;
	}
	if (t > 10) {
		t = pow(10,-4);
	}
	this->SetEdgeLength(p->id,c->id,t);
	this->ComputeLogLikelihood();
	float valueToReturn = -1 * this->logLikelihood;
	return (valueToReturn);
}

Matrix4f rootedPhylogeny_tree::GetRateMatrixForFreeParameters(MatrixXd B){
	Matrix4f Q;
	float pi_1; float pi_2; float pi_3; float pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;float i; float j; float k; float l;
	pi_1 = B(0,0); pi_2 = B(1,0); pi_3 = B(2,0);
	pi_4 = (1 - pi_1 - pi_2 - pi_3);
	// Construct Q
	a = B(3,0); b = B(4,0); c = B(5,0); d = B(6,0); e = B(7,0);
	f = B(8,0); g = B(9,0); h = B(10,0);

	i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3);
	j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4;
	k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4;
	l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4);
	
	Q << -(a+b+c), a, b, c,
		  d, -(d+e+f), e, f,
		  g, h, -(g+h+i), i,
		  j, k, l, -(j+k+l);
	
	return (Q);
}

void rootedPhylogeny_tree::SetVerticesForPreOrderTreeTraversal(){
	this->verticesForPreOrderTreeTraversal->clear();
	this->verticesForPreOrderTreeTraversal->push_back(this->root);
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	verticesToVisit.push_back(this->root);
	int numberOfVerticesToVisit = verticesToVisit.size();
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	while (numberOfVerticesToVisit > 0){
		p = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (int c_id : p->children_id){
			c = (*this->vertexMap)[c_id];
			this->verticesForPreOrderTreeTraversal->push_back(c);
			if (c->children_id.size() > 0){
				verticesToVisit.push_back(c);				
				numberOfVerticesToVisit += 1;
			}
		}
	}
}

void rootedPhylogeny_tree::SetRateCategories(float baseFreqThreshold){
	this->SetVerticesForPreOrderTreeTraversal();
	// Select change points;
	this->changePoints->clear();
	this->changePoints->push_back(this->root);
	for (pair<rootedPhylogeny_vertex*,float> descendantAndBaseFreqChange : *this->baseFreqChangeForEachDescendant){		
		if (descendantAndBaseFreqChange.second > baseFreqThreshold){
			this->changePoints->push_back(descendantAndBaseFreqChange.first);
		}
	}
	this->numberOfRateCategories = this->changePoints->size();
//	cout << "number of rate categories is " << this->numberOfRateCategories << endl;
	int rateCategory = -1;
	
	for (rootedPhylogeny_vertex * v : (*this->verticesForPreOrderTreeTraversal)){
		vector <rootedPhylogeny_vertex*>::iterator it = find(this->changePoints->begin(),this->changePoints->end(),v);
		if (it == this->changePoints->end()){
			rateCategory = (*this->vertexMap)[v->parent_id]->rateCategory;
			v->rateCategory = rateCategory;			
			
		} else {			
			rateCategory = distance(this->changePoints->begin(),it);
			v->rateCategory = rateCategory;
		}
		if (rateCategory == -1) {
			cout << "check assignment of rate categories" << endl;
		}
	}
}

void rootedPhylogeny_tree::OptimizeModelParametersDep(){	
	// Compute initial estimate of model parameters
	this->InitializeModelParametersForNelderMead();
	MatrixXd B;
	int rateCat;
	Matrix4f Q;
	for (pair<int,MatrixXd> rateCatParamPair : *this->parametersPerRateCategory){
		rateCat = rateCatParamPair.first;
		B = rateCatParamPair.second;		
		Q = this->GetRateMatrixForFreeParameters(B);
		this->rateMatrixPerRateCategory->insert(pair<int,Matrix4f>(rateCat,Q));
	}	
	this->SetMinLengthOfEdges(); // 10^-7;
	//	cout << "here 1" << endl;
	// Iterate over rate category
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++){
		this->OptimizeModelParametersForRateCatDep(rateCat);
	}
}

void rootedPhylogeny_tree::OptimizeModelParametersForRateCatDep(int rateCat){
	//	using namespace brent;
	if (1){
		double convergenceTolerance = pow(10,-4);
		MatrixXd InverseHessian_current = ArrayXXd::Zero(11,11);
		MatrixXd InverseHessian_updated;
		MatrixXd B_current; MatrixXd B_updated;
		MatrixXd Jacobian_current; MatrixXd Jacobian_updated;
		MatrixXd s;
		MatrixXd y;
		MatrixXd rho = ArrayXXd::Zero(11,11);;
		MatrixXd IdentityMatrix = ArrayXXd::Zero(11,11);
		MatrixXd TempMatrix_l = ArrayXXd::Zero(11,11);
		MatrixXd TempMatrix_r = ArrayXXd::Zero(11,11);
		int iter;
		double normOfJacobian;
		double stepSize;
		this->rateCategoryForOptimization = rateCat;
		
		for (int par_1 = 0; par_1 < 11; par_1 ++){		
			InverseHessian_current(par_1, par_1) = 1.0;
			IdentityMatrix(par_1, par_1) = 1.0;
		}
		
		B_current = (*this->parametersPerRateCategory)[rateCat];
		Jacobian_current = this->GetJacobianForRateCategory(rateCat);
		normOfJacobian = this->GetNormOfJacobian(Jacobian_current);
		iter = 0;
//		cout << "norm of Jacobian for iteration " << iter << " is " << normOfJacobian << endl;
		while (normOfJacobian > convergenceTolerance and iter < 100){
			cout << "B_current is \n" << B_current << endl;
//			cout << "norm of Jacobian for iteration " << iter << " is " << normOfJacobian << endl;
			iter += 1;
//			cout << "Jacobian current is " << endl;
//			cout << Jacobian_current << endl;
			this->searchDirection = -1 * (InverseHessian_current * Jacobian_current);
//			cout << "search direction" << endl;
//			cout << this->searchDirection << endl;
//			cout << "B_current is \n" << B_current << endl;
			stepSize = this->GetOptimalStepSizeUsingLineSearch(B_current);
			cout << "step size is " << stepSize << endl;			
			B_updated = B_current + stepSize * this->searchDirection;					
			this->SetParametersForRateMatrixForBFGS(B_updated,rateCat);
			Jacobian_updated = this->GetJacobianForRateCategory(rateCat);
			y = Jacobian_updated - Jacobian_current;
			s = B_updated - B_current;
//			cout << "y" << endl;
//			cout << y << endl;
//			cout << "s" << endl;
//			cout << s << endl;
			rho = y.transpose()*s;
			rho = rho.inverse();
//			cout << "rho is " << endl;
//			cout << rho << endl;
//			cout << s * y.transpose() << endl;
			TempMatrix_l = IdentityMatrix - rho(0,0) * ( s * y.transpose() );
			TempMatrix_r = IdentityMatrix - rho(0,0) * ( y * s.transpose() );
			InverseHessian_updated = TempMatrix_l * InverseHessian_current;
			InverseHessian_updated = InverseHessian_updated * TempMatrix_r;
			InverseHessian_updated += rho(0,0) * ( s * s.transpose() );
			B_current = B_updated;
			InverseHessian_current = InverseHessian_updated;
			normOfJacobian = this->GetNormOfJacobian(Jacobian_updated);
		}
		
	} else {
//		double t = pow(10,-4);
//		double lower_limit = -4.0;
//		double upper_limit = 4.0/3.0;
		double x = 0.1;
		double result;
//		result = this->BrentLineSearchForSelectingStepSize(lower_limit, upper_limit, t, x);
//		result = local_min(lower_limit, upper_limit, t, this->SampleFunctionForMinimization(x), x);
//		result = local_min (lower_limit, upper_limit, t, this->SampleFunctionForMinimization, x);
//		result = local_min (lower_limit, upper_limit, t, this->SampleFunctionForMinimization, x);
		cout << "optimal x is " << x << endl;
		cout << "value of function at minima is " << result << endl;
	}
//	double stepSize;
//	double negLogLikelihood;
//	negLogLikelihood = local_min (lower_limit, upper_limit, t, this->SampleFunctionForMinimization, stepSize);	
}

void rootedPhylogeny_tree::InitializeModelParametersForBFGS(){
	this->rateMatrixPerRateCategory->clear();
	this->scalingFactorForRateCategory->clear();
	this->parametersPerRateCategory->clear();
	bool resetAllElementsOrRateMatrix;
	// Initialize substitution matrices for each rate category
	map <int, Matrix4f> substitutionMatrixPerRateCategory;
	MatrixXf stationaryDistribution;  
	Matrix4f P;
	// Populate substitution matrices for each rate category
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	int rateCat; Matrix4f Q;
	unsigned char dna_p; unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			if (substitutionMatrixPerRateCategory.find(c->rateCategory) == substitutionMatrixPerRateCategory.end()){
				P = ArrayXXf::Zero(4,4);
				substitutionMatrixPerRateCategory.insert(pair<int,Matrix4f>(c->rateCategory,P));
			}
			p = (*this->vertexMap)[c->parent_id];
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){
				dna_p = p->compressedSequence[site];
				dna_c = c->compressedSequence[site];
				substitutionMatrixPerRateCategory[c->rateCategory](dna_p,dna_c) += this->siteWeights[site];
			}
		}
	}	
	float rowSum; float scalingFactor; MatrixXd parameters;
	for (pair<int,Matrix4f> rateCatSubsMatrixPair : substitutionMatrixPerRateCategory){
		resetAllElementsOrRateMatrix = 0;
		rateCat = rateCatSubsMatrixPair.first;
		P = rateCatSubsMatrixPair.second;
		for (int row = 0; row < 4; row++){
			rowSum = 0;
			for (int col = 0; col < 4; col++){
				rowSum += P(row,col);
			}
			for (int col = 0; col < 4; col++){
				P(row,col) /= rowSum;
			}
		}
				
		Q = P.log();			
		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col){
					if (Q(row,col) < 0){
						Q(row,col) *= -1;
					}
				}
			}
		}
		Q /= Q(3,2);
//		cout << Q << endl;
//		cout << "Is element Q(0,1) nan ? " << isnan(Q(0,1));		
		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col){
					if (isnan(Q(row,col)) or isinf(Q(row,col))){
						resetAllElementsOrRateMatrix = 1;
					}
				} else {
					if (abs(Q(row,col)) < pow(10,-6)) {
						resetAllElementsOrRateMatrix = 1;
					}
				}
			}
		}
		
		uniform_real_distribution<> distribution(0.0, 1.0);	
		if (resetAllElementsOrRateMatrix){
			for (int row = 0; row < 4; row ++){
				for (int col = 0; col < 4; col ++){
					if (row != col){						
						Q(row,col) = distribution(generator);
					}
				}
			}		
		}		
		Q(0,0) = - (Q(0,1) + Q(0,2) + Q(0,3));
		Q(1,1) = - (Q(1,0) + Q(1,2) + Q(1,3));
		Q(2,2) = - (Q(2,0) + Q(2,1) + Q(2,3));
		Q(3,3) = - (Q(3,0) + Q(3,1) + Q(3,2));
		scalingFactor = this->ComputeScalingFactor(Q);
		Q/=scalingFactor;
		this->rateMatrixPerRateCategory->insert(pair<int,Matrix4f>(rateCat,Q));
//		cout << "Initial rate matrix is " << endl;
//		cout << Q << endl;
//		cout << "Scaling factor is " << endl;
		scalingFactor = this->ComputeScalingFactor(Q);
//		cout << scalingFactor << endl;
		this->scalingFactorForRateCategory->insert(pair<int,float>(rateCat,scalingFactor));
//		cout << "Stationary distribution is " << endl;
		stationaryDistribution = this->ComputeStationaryDistribution(Q);
		cout << stationaryDistribution << endl;
		parameters = ArrayXXd::Zero(11,1);
		// pi_1
		parameters(0,0) = stationaryDistribution(0,0); 
		// pi_2
		parameters(1,0) = stationaryDistribution(1,0);
		// pi_3
		parameters(2,0) = stationaryDistribution(2,0);
		// a
		parameters(3,0) = Q(0,1); 
		// b
		parameters(4,0) = Q(0,2);
		// c
		parameters(5,0) = Q(0,3);
		// d
		parameters(6,0) = Q(1,0);
		// e
		parameters(7,0) = Q(1,2);
		// f
		parameters(8,0) = Q(1,3);
		// g
		parameters(9,0) = Q(2,0);
		// h
		parameters(10,0) = Q(2,1);		
		for (int par = 0; par < 11; par++) {
			if (parameters(par,0) < 0) {
				cout << "Initialized rate matrix has negative entries on non-diagonal elements" << endl;
				cout << Q << endl;
			}
		}
		this->parametersPerRateCategory->insert(pair<int,MatrixXd>(rateCat,parameters));
	}
	
	if (rateMatrixPerRateCategory->find(this->root->rateCategory) == rateMatrixPerRateCategory->end()){
//		cout << "rate cat for root is " << this->root->rateCategory << endl;
//		cout << "rate cat for left child of root is " << (*this->vertexMap)[this->root->children_id[0]]->rateCategory << endl;
//		cout << "rate cat for right child of root is " << (*this->vertexMap)[this->root->children_id[1]]->rateCategory << endl;
		// There is no rate matrix that is associated with the rate category root.cat
//		cout << " There is no rate matrix that is associated with the rate category root.cat" << endl;
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] = 0.0;
		}
		unsigned char dna_p;
		float sequenceLength = 0;
		for (unsigned int site = 0; site < this->siteWeights.size(); site++){
			dna_p = this->root->compressedSequence[site];
			this->rootProbability[dna_p] += this->siteWeights[site];
			sequenceLength += this->siteWeights[site];
		}
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] /= sequenceLength;
		}
		
	} else {
		// There is one rate matrix that is associated with the rate category root.cat
//		cout << "Rate matrix associated with the root is " << endl;
		Q = (*this->rateMatrixPerRateCategory)[this->root->rateCategory];		
//		cout << Q << endl;
//		cout << "Stationary distribution that is associated with the root is " << endl;
		MatrixXf stationaryDistribution = this->ComputeStationaryDistribution(Q);		
//		cout << stationaryDistribution << endl;
		this->rootProbability[0] = stationaryDistribution(0,0);
		this->rootProbability[1] = stationaryDistribution(1,0);
		this->rootProbability[2] = stationaryDistribution(2,0);
		this->rootProbability[3] = stationaryDistribution(3,0);
	}
}

void rootedPhylogeny_tree::InitializeModelParametersForNelderMead(){
	bool verbose = 1;
	this->rateMatrixPerRateCategory->clear();
	this->scalingFactorForRateCategory->clear();
	this->parametersPerRateCategory->clear();
	bool resetAllElementsOrRateMatrix = 0;
	// Initialize substitution matrices for each rate category
	map <int, Matrix4f> substitutionMatrixPerRateCategory;
	Matrix4f P;
	// Populate substitution matrices for each rate category
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	int rateCat; Matrix4f Q;
	unsigned char dna_p; unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			if (substitutionMatrixPerRateCategory.find(c->rateCategory) == substitutionMatrixPerRateCategory.end()){
				P = ArrayXXf::Zero(4,4);
				substitutionMatrixPerRateCategory.insert(pair<int,Matrix4f>(c->rateCategory,P));
			}
			p = (*this->vertexMap)[c->parent_id];			
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){
				dna_p = p->compressedSequence[site];
				dna_c = c->compressedSequence[site];
				substitutionMatrixPerRateCategory[c->rateCategory](dna_p,dna_c) += this->siteWeights[site];
			}			
		}
	}	
	float rowSum; float scalingFactor; MatrixXd parameters;
	for (pair<int,Matrix4f> rateCatSubsMatrixPair : substitutionMatrixPerRateCategory){
		resetAllElementsOrRateMatrix = 0;
		rateCat = rateCatSubsMatrixPair.first;
		P = rateCatSubsMatrixPair.second;
		for (int row = 0; row < 4; row++){
			rowSum = 0;
			for (int col = 0; col < 4; col++){
				rowSum += P(row,col);
			}
			for (int col = 0; col < 4; col++){
				P(row,col) /= rowSum;
			}
		}
		
		float maxNonDiagElementOfUnscaledQ = 0;
		Q = P.log();
		Q(0,0) = - (Q(0,1) + Q(0,2) + Q(0,3));
		Q(1,1) = - (Q(1,0) + Q(1,2) + Q(1,3));
		Q(2,2) = - (Q(2,0) + Q(2,1) + Q(2,3));
		Q(3,3) = - (Q(3,0) + Q(3,1) + Q(3,2));			
		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col){
					if (Q(row,col) < 0){
						Q(row,col) *= -1;
					}
					if (maxNonDiagElementOfUnscaledQ < Q(row,col) and !isnan(Q(row,col))) {
						maxNonDiagElementOfUnscaledQ = Q(row,col);
					}
				}
			}
		}
		uniform_real_distribution<> distribution_allZeroRow(0.0, maxNonDiagElementOfUnscaledQ);	
		// if each element of a row contains 0 then set each element in row to a uniform (0,maxNonDiagElementOfUnscaledQ)
		for (int row = 0; row < 4; row++) {
			if (abs(Q(row,row)) < pow(10,-8)) {
				for (int col = 0; col < 4; col ++) {
					Q(row,col) = distribution_allZeroRow(generator);
				}
			}
		}
		
		Q /= Q(3,2);
		 
//		cout << Q << endl;
//		cout << "Is element Q(0,1) nan ? " << isnan(Q(0,1));		
		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col){
					if (isnan(Q(row,col)) or isinf(Q(row,col)) or abs(Q(row,col)) > pow(10,2)){
						resetAllElementsOrRateMatrix = 1;
						cout << "Q contains a nan element" << endl;
					}
				} else {
					if (abs(Q(row,col)) < pow(10,-8)) {
						cout << "Q contains a zero on the diagonal" << endl;
						resetAllElementsOrRateMatrix = 1;
					}
				}
			}
		}
//		cout << "Rate matrix for category " << rateCat << "before resetting elements is " << endl;
//		cout << Q << endl;
		uniform_real_distribution<> distribution(0.0, 1.0);	
		if (resetAllElementsOrRateMatrix){
			for (int row = 0; row < 4; row ++){
				for (int col = 0; col < 4; col ++){
					if (row != col){						
						Q(row,col) = distribution(generator);
					}
				}
			}		
		}

		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col and Q(row,col) < pow(10,-8)){						
					Q(row,col) = distribution(generator);
				}
			}
		}		
		Q /= Q(3,2);				
		Q(0,0) = - (Q(0,1) + Q(0,2) + Q(0,3));
		Q(1,1) = - (Q(1,0) + Q(1,2) + Q(1,3));
		Q(2,2) = - (Q(2,0) + Q(2,1) + Q(2,3));
		Q(3,3) = - (Q(3,0) + Q(3,1) + Q(3,2));
		if (verbose) {
			cout << "Initial rate matrix for category " << rateCat << " is " << endl;
			cout << Q << endl;	
		}
		for (int row = 0; row < 4; row ++) {
			for (int col = 0; col < 4; col ++) {
				if (isnan(Q(row,col))) {
					cout << "Initialized rate matrix contains nan" << endl;
				}
			}
		}		
//		cout << "Is element Q(0,1) nan ? " << isnan(Q(0,1)) << endl;
		this->rateMatrixPerRateCategory->insert(pair<int,Matrix4f>(rateCat,Q));
		scalingFactor = this->ComputeScalingFactor(Q);		
		this->scalingFactorForRateCategory->insert(pair<int,float>(rateCat,scalingFactor));
		parameters = ArrayXXd::Zero(11,1);
		// a
		parameters(0,0) = Q(0,1); 
		// b
		parameters(1,0) = Q(0,2);
		// c
		parameters(2,0) = Q(0,3);
		// d
		parameters(3,0) = Q(1,0);
		// e
		parameters(4,0) = Q(1,2);
		// f
		parameters(5,0) = Q(1,3);
		// g
		parameters(6,0) = Q(2,0);
		// h
		parameters(7,0) = Q(2,1);
		// i
		parameters(8,0) = Q(2,3);
		// j
		parameters(9,0) = Q(3,0);	
		// k
		parameters(10,0) = Q(3,1);
		for (int par = 0; par < 11; par++){
			if (parameters(par,0) < 0){
				cout << "Initialized rate matrix has negative entries on non-diagonal elements" << endl;
				cout << Q << endl;
			}
		}
		this->parametersPerRateCategory->insert(pair<int,MatrixXd>(rateCat,parameters));
	}
	
	if (rateMatrixPerRateCategory->find(this->root->rateCategory) == rateMatrixPerRateCategory->end()){
//		cout << "rate cat for root is " << this->root->rateCategory << endl;
//		cout << "rate cat for left child of root is " << (*this->vertexMap)[this->root->children_id[0]]->rateCategory << endl;
//		cout << "rate cat for right child of root is " << (*this->vertexMap)[this->root->children_id[1]]->rateCategory << endl;
		// There is no rate matrix that is associated with the rate category root.cat
//		cout << " There is no rate matrix that is associated with the rate category root.cat" << endl;
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] = 0.0;
		}
		unsigned char dna_p;
		float sequenceLength = 0;
		for (unsigned int site = 0; site < this->siteWeights.size(); site++){
			dna_p = this->root->compressedSequence[site];
			this->rootProbability[dna_p] += this->siteWeights[site];
			sequenceLength += this->siteWeights[site];
		}
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] /= sequenceLength;
		}
		
	} else {
		// There is one rate matrix that is associated with the rate category root.cat
		Q = (*this->rateMatrixPerRateCategory)[this->root->rateCategory];
		MatrixXf stationaryDistribution = this->ComputeStationaryDistribution(Q);
		this->rootProbability[0] = stationaryDistribution(0,0);
		this->rootProbability[1] = stationaryDistribution(1,0);
		this->rootProbability[2] = stationaryDistribution(2,0);
		this->rootProbability[3] = stationaryDistribution(3,0);
	}
}

float rootedPhylogeny_tree::ComputeAbsDifferenceInBaseFreq(array<float, 4> baseFreq_1, array<float, 4> baseFreq_2){	
	// A is 0, C is 1, G is 2, T is 3
	float abs_nt_diff = 0;
	for (int dna = 0; dna < 4; dna ++){
		abs_nt_diff += abs(baseFreq_1[dna] - baseFreq_2[dna]);
	}	
	return (abs_nt_diff);	
}


void rootedPhylogeny_tree::SetThresholds(){
	this->baseFreqChangeForEachDescendant->clear();	
	this->thresholdList->clear();
	vector <pair<float,rootedPhylogeny_vertex*>> negBaseFreqDiffAndChangePointList;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	float negBaseFreqDiff;
	float baseFreqDiff;
	// goodness-of-clustering vs goodness-of-fit
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			negBaseFreqDiff = -1 * this->ComputeAbsDifferenceInBaseFreq(p->baseFreq, c->baseFreq);
			this->baseFreqChangeForEachDescendant->insert(pair<rootedPhylogeny_vertex*,float>(c,this->ComputeAbsDifferenceInBaseFreq(p->baseFreq, c->baseFreq)));
			negBaseFreqDiffAndChangePointList.push_back(pair<float,rootedPhylogeny_vertex*>(negBaseFreqDiff,c));
		}		
	}
	sort(negBaseFreqDiffAndChangePointList.begin(), negBaseFreqDiffAndChangePointList.end());
//	int numberOfRateCategories = 0;
	for (pair<float,rootedPhylogeny_vertex*> negBaseFreqDiffAndChangePoint : negBaseFreqDiffAndChangePointList){
		baseFreqDiff = -1 * negBaseFreqDiffAndChangePoint.first;
		if (find(this->thresholdList->begin(),this->thresholdList->end(),baseFreqDiff) == this->thresholdList->end()){
//			cout << baseFreqDiff << endl;
			this->thresholdList->push_back(baseFreqDiff);
		}
	}
	this->thresholdList->push_back(-1.0);
}

void rootedPhylogeny_tree::ComputeBaseFreq(){
	rootedPhylogeny_vertex * v;
	int dna;
	float sequenceLength = 0;
	for (unsigned int site = 0; site < this->siteWeights.size(); site ++){
		sequenceLength += this->siteWeights[site];
	}	
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		v = idPtrPair.second;
		for (unsigned int site = 0; site < this->siteWeights.size(); site ++){
			dna = v->compressedSequence[site];
			v->baseFreq[dna] += this->siteWeights[site];
		}
		for (int dna = 0; dna < 4; dna ++){
			v->baseFreq[dna] /= sequenceLength;
		}
	}	
}

double rootedPhylogeny_tree::SampleFunctionForMinimization(double x){	
	return ((x+3)*pow((x-1),2.0));
}

double rootedPhylogeny_tree::ComputeEdgeLogLikelihood(rootedPhylogeny_vertex * p, rootedPhylogeny_vertex * c, Matrix4f P){
	unsigned char dna_p; unsigned char dna_c;
	double edgeLogLik = 0;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		dna_p = p->compressedSequence[site];
		dna_c = c->compressedSequence[site];
		edgeLogLik += log(P(dna_p,dna_c))*this->siteWeights[site];
	}
	return (edgeLogLik);
}

void rootedPhylogeny_tree::AddEdgeLength(int p_id, int c_id, float edgeLength){
	if (p_id < c_id) {
		this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(p_id,c_id),edgeLength));
	} else {
		this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(c_id,p_id),edgeLength));
	}
	
}

void rootedPhylogeny_tree::SetEdgeLength(int p_id, int c_id, float edgeLength){
	if (p_id < c_id){
		(*this->edgeLengthsMap)[pair<int,int>(p_id,c_id)] = edgeLength;
	} else {
		(*this->edgeLengthsMap)[pair<int,int>(c_id,p_id)] = edgeLength;
	}
}

int rootedPhylogeny_tree::GetVertexId(string v_name){
	rootedPhylogeny_vertex * v;
	int idToReturn = -10;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		v = idPtrPair.second;
		if (v->name == v_name){
			idToReturn = v->id;						
		}
	}
	if (idToReturn == -10){
		cout << "check id to return" << endl;
	}
	return (idToReturn);
}

void rootedPhylogeny_tree::ReadSequenceFile(string sequenceFileName){
	map <string, vector <unsigned char>> fullSequencesMap;
	map <string, vector <unsigned char>> compressedSequencesMap;
	vector <vector <unsigned char>> compressedSequences;
	vector <unsigned char> recodedSequence;
//	unsigned char recodedDna;
	unsigned int site = 0;
	ifstream inputFile(sequenceFileName.c_str());
	string seqName;
	string seq = "";	
	int vertex_ind = 0;
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != ""){
//				sequenceNames.push_back(seqName);
				for (char const dna: seq){
					recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);
//					recodedDna = recodedSequence[site];
					site += 1;
					}
				fullSequencesMap.insert(pair<string,vector<unsigned char>>(seqName,recodedSequence));
//				cout << seqName << "\t";
//				cout << EncodeAsDNA(fullSequencesMap[seqName]) << endl;
//				MST_ptr->AddVertex(vertex_ind,seqName,recodedSequence);
				vertex_ind += 1;
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
	for (char const dna: seq){				
		recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);
//		recodedDna = recodedSequence[site];
		site += 1;
	}
	this->sequenceLength = (int) recodedSequence.size();
	this->logLikelihoodConvergenceThreshold = pow(10,-5) * this->sequenceLength;
//	int sequenceLength = recodedSequence.size();
//	MST_ptr->AddVertex(vertex_ind,seqName,recodedSequence);
	fullSequencesMap.insert(pair<string,vector<unsigned char>>(seqName,recodedSequence));
	recodedSequence.clear();
//	sequenceNames.push_back(seqName);
	inputFile.close();
//	cout << "number of sequences is " << fullSequencesMap.size() << endl;
	vector<vector<unsigned char>> fullSequences;
	vector <string> sequenceNames;
	for (pair<string,vector<unsigned char>> seqNameAndSeq : fullSequencesMap){
		fullSequences.push_back(seqNameAndSeq.second);
		sequenceNames.push_back(seqNameAndSeq.first);
	}	
	string v_name;
	vector <unsigned char> compressedSequence;
	int v_id;
	tie (compressedSequences, this->siteWeights) = GetCompressedSequencesAndSiteWeights(fullSequences);	
	for (unsigned int v_num = 0; v_num < fullSequences.size(); v_num ++){
		compressedSequence = compressedSequences[v_num];
		v_name = sequenceNames[v_num];		
		v_id = this->GetVertexId(v_name);
//		cout << "vertex name is " << v_name << "\t" << "vertex id is " << v_id << endl;
//		cout << "compressed sequence is " << EncodeAsDNA(compressedSequence) << endl;
		(*this->vertexMap)[v_id]->compressedSequence = compressedSequence;		
	}	
}


void rootedPhylogeny_tree::ReadUndirectedEdgeList(string treeFileName) {
	string u_name;
	string v_name;
	int u_id;
	int v_id;
	float edgeLength;
	vector <string> splitLine;
	vector <string> leafNames;
	vector <string> ancestorNames;
	vector <string> nonRootVertexNames;	
//	string rootName = "";
	vector <unsigned char> emptySequence;
	v_id = 0;
	ifstream edgeListFile(treeFileName.c_str());
	for (string line; getline(edgeListFile, line);){		
		boost::split(splitLine, line, [](char c){return c == '\t';});		
		u_name = splitLine[0];		
		v_name = splitLine[1];
		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),v_name) == nonRootVertexNames.end()){
			nonRootVertexNames.push_back(v_name);
		}		
		if (find(ancestorNames.begin(),ancestorNames.end(),u_name)==ancestorNames.end()){
			ancestorNames.push_back(u_name);
		}
		if (find(leafNames.begin(),leafNames.end(),v_name)==leafNames.end()){
			if(!boost::starts_with(v_name, "h_")){
				leafNames.push_back(v_name);
			}
		}
	}
	for (string name: leafNames){
		rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(v_id,name,emptySequence);
		this->vertexMap->insert(pair<int,rootedPhylogeny_vertex*>(v_id,v));
		v_id += 1;
	}
	// Remove root from ancestor names
//	for (string name: ancestorNames){
//		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),name)==nonRootVertexNames.end()){
//			rootName = name;
//		}		
//	}
	this->numberOfObservedSequences = leafNames.size();
//	cout << "number of observed sequences is " << this->numberOfObservedSequences << endl;
	// Change root name
//	(*this->vertexMap)[-1]->name = rootName;
//	ancestorNames.erase(remove(ancestorNames.begin(), ancestorNames.end(), rootName), ancestorNames.end());	
	for (string name: ancestorNames){
		rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(v_id,name,emptySequence);
		this->vertexMap->insert(pair<int,rootedPhylogeny_vertex*>(v_id,v));
		v_id += 1;
	}	
	edgeListFile.clear();
	edgeListFile.seekg(0, ios::beg);	
	for (string line; getline(edgeListFile, line);){
		boost::split(splitLine, line, [](char c){return c == '\t';});
		u_name = splitLine[0];
		v_name = splitLine[1];
		u_id = this->GetVertexId(u_name);
		v_id = this->GetVertexId(v_name);		
		(*this->vertexMap)[u_id]->AddNeighbor(v_id);
		(*this->vertexMap)[v_id]->AddNeighbor(u_id);
		edgeLength = stof(splitLine[2]);
		this->AddEdgeLength(u_id,v_id,edgeLength);
	}
	edgeListFile.close();
}

void rootedPhylogeny_tree::ReadDirectedEdgeListForBifurcatingTree(string treeFileName){
	string u_name;
	string v_name;
	int u_id;
	int v_id;
	float edgeLength;
	vector <string> splitLine;
	vector <string> leafNames;
	vector <string> ancestorNames;
	vector <string> nonRootVertexNames;	
	string rootName = "";
	vector <unsigned char> emptySequence;
	v_id = 0;
	ifstream edgeListFile(treeFileName.c_str());
	for (string line; getline(edgeListFile, line);){		
		boost::split(splitLine, line, [](char c){return c == '\t';});		
		u_name = splitLine[0];		
		v_name = splitLine[1];
		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),v_name) == nonRootVertexNames.end()){
			nonRootVertexNames.push_back(v_name);
		}		
		if (find(ancestorNames.begin(),ancestorNames.end(),u_name)==ancestorNames.end()){
			ancestorNames.push_back(u_name);
		}
		if (find(leafNames.begin(),leafNames.end(),v_name)==leafNames.end()){
			if(!boost::starts_with(v_name, "h_")){
				leafNames.push_back(v_name);
			}
		}
	}	
	for (string name: leafNames){
		rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(v_id,name,emptySequence);
		this->vertexMap->insert(pair<int,rootedPhylogeny_vertex*>(v_id,v));
		v_id += 1;
	}
	// Remove root from ancestor names
	for (string name: ancestorNames){
		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),name)==nonRootVertexNames.end()){
			rootName = name;
		}		
	}
	this->numberOfObservedSequences = leafNames.size();
//	cout << "number of observed sequences is " << this->numberOfObservedSequences << endl;
	// Change root name
	(*this->vertexMap)[-1]->name = rootName;
	ancestorNames.erase(remove(ancestorNames.begin(), ancestorNames.end(), rootName), ancestorNames.end());	
	for (string name: ancestorNames){		
		rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(v_id,name,emptySequence);
		this->vertexMap->insert(pair<int,rootedPhylogeny_vertex*>(v_id,v));
		v_id += 1;
	}	
	edgeListFile.clear();
	edgeListFile.seekg(0, ios::beg);	
	for (string line; getline(edgeListFile, line);){
		boost::split(splitLine, line, [](char c){return c == '\t';});
		u_name = splitLine[0];
		v_name = splitLine[1];
		u_id = this->GetVertexId(u_name);
		v_id = this->GetVertexId(v_name);
		(*this->vertexMap)[u_id]->AddChild(v_id);
		(*this->vertexMap)[v_id]->AddParent(u_id);
		edgeLength = stof(splitLine[2]);
		this->AddEdgeLength(u_id,v_id,edgeLength);
	}
	edgeListFile.close();
	this->SetLeaves();
}


void rootedPhylogeny_tree::ComputeLogLikelihoodForFullyLabeledTree(){
	this->logLikelihood = 0;
	Matrix4f P; Matrix4f Q; Matrix4f Q_scaled;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	float t;	
	unsigned char dna_p; unsigned char dna_c;
	unsigned char dna;
	float scalingFactor;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id, c->id);			
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];			
			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
			Q_scaled = Q*(t/scalingFactor);
			P = Q_scaled.exp();			
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){
				dna_p = p->compressedSequence[site];
				dna_c = c->compressedSequence[site];			
				this->logLikelihood += log(P(dna_p,dna_c))*this->siteWeights[site];				
			}
		}
	}
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		dna = (*this->vertexMap)[-1]->compressedSequence[site];
		this->logLikelihood += log(this->rootProbability[dna])*this->siteWeights[site];
	}
}

void rootedPhylogeny_tree::SetMinLengthOfEdges(){
	for (pair<pair<int,int>,float> edgeIdsAndLengthsPair: *this->edgeLengthsMap){		
		if (edgeIdsAndLengthsPair.second < pow(10,-7)){
			(*this->edgeLengthsMap)[edgeIdsAndLengthsPair.first]  = pow(10,-7);
		}
	}
}

void rootedPhylogeny_tree::EstimateAncestralSequencesByFittingTheGMMUsingLongDoublePrecision(){	
	this->ComputeMPEstimateOfAncestralSequences();
	map <int,Matrix4f> transitionMatrices;	
	map <int,array<long double,4>> conditionalLikelihoodMap;
	array <long double,4> conditionalLikelihood;	
	long double partialLikelihood;
	long double siteLikelihood;
	long double currentLogLikelihood = 0;
	long double previousLogLikelihood = 0;
	int iter = 0;
	int maxIter = 10;
	pair<int,int> edgeIds;
	unsigned char dna_p; unsigned char dna_c;
	char maxProbState;
	float rowSum;
	long double currentProb;
	long double maxProb;
	int numberOfVerticesToVisit;
	bool continueEM = 1;
	vector <rootedPhylogeny_vertex*> verticesToVisit;	
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f P;
	this->SetEdgesForPostOrderTreeTraversal();
	// Iterate till convergence of log likelihood
		while (continueEM and iter < maxIter) {
			iter += 1;			
//			cout << "root sequence is " << endl;
//			cout << EncodeAsDNA(this->root->compressedSequence) << endl;
			currentLogLikelihood = 0;
			// Estimate root probablity
			this->ComputeMLEOfRootProbability();		
			// Estimate transition matrices	
//			cout << "here 1" << endl;
			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
				c = idPtrPair.second;
				if (c->id != -1){
					p = (*this->vertexMap)[c->parent_id];
					P = ArrayXXf::Zero(4,4);
					for (unsigned int site = 0; site < siteWeights.size(); site++){
						dna_p = p->compressedSequence[site];
						dna_c = c->compressedSequence[site];
						P(dna_p,dna_c) += this->siteWeights[site];
					}
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
						rowSum = 0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							rowSum += P(dna_p,dna_c);
						}
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							 P(dna_p,dna_c) /= rowSum;
						}
					}
					transitionMatrices.insert(pair<int,Matrix4f>(c->id,P));
				}
			}
//			cout << "here 2" << endl;
			// Estimate ancestral sequences
			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
				c = idPtrPair.second;
				if (c->children_id.size() > 0){
					c->compressedSequence.clear();
				}
			}		
			// Iterate over sites		
			for (unsigned int site = 0 ; site < this->siteWeights.size(); site++){
				conditionalLikelihoodMap.clear();			
//				cout << "site is " << site << endl;
				for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){
					p = (*this->vertexMap)[edgeIds.first];
					c = (*this->vertexMap)[edgeIds.second];
					P = transitionMatrices[c->id];					
					// Initialize conditional likelihood for leaves
					if (c->children_id.size()==0){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
							conditionalLikelihood[dna_c] = 0;
						}
						conditionalLikelihood[c->compressedSequence[site]] = 1;
						conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c->id,conditionalLikelihood));
					}
					// Initialize conditional likelihood for ancestors
					if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
						conditionalLikelihood[dna_c] = 1;
						}				
						conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));					
					}		
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
						partialLikelihood = 0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						}
						conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
					}
				}
				maxProbState = -1;
				maxProb = 0;
				siteLikelihood = 0; 							
//				cout << prior << endl;
				for (int dna = 0; dna < 4; dna++) {
					currentProb = this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
					siteLikelihood += currentProb;
					if (currentProb > maxProb) {
						maxProb = currentProb;
						maxProbState = dna;		
					}
				}				
				currentLogLikelihood += log(siteLikelihood) * this->siteWeights[site];
				if (maxProbState == -1) {
					cout << "check state estimation" << endl;					
				}
//				for (int dna = 0; dna < 4; dna++){
//					cout << "conditional likelihood for dna " << dna << " is " << conditionalLikelihoodMap[-1][dna] << endl;
//				}
				(*this->vertexMap)[-1]->compressedSequence.push_back(maxProbState);
				verticesToVisit.clear();			
				for (int c_id: (*this->vertexMap)[-1]->children_id) {
					if ((*this->vertexMap)[c_id]->children_id.size()>0) {
						verticesToVisit.push_back((*this->vertexMap)[c_id]);
					}
				}
				numberOfVerticesToVisit = verticesToVisit.size();
				while (numberOfVerticesToVisit > 0) {
					c = verticesToVisit[numberOfVerticesToVisit-1];				
					verticesToVisit.pop_back();
					numberOfVerticesToVisit -= 1;
					p = (*this->vertexMap)[c->parent_id];
					P = transitionMatrices[c->id];
					dna_p = p->compressedSequence[site];
					maxProbState = -1;
					maxProb = 0;
					for (int dna_c = 0; dna_c < 4; dna_c++){
						currentProb = P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						if (currentProb > maxProb){
							maxProb = currentProb;
							maxProbState = dna_c;
						}					
					}
					if (maxProbState == -1){
						cout << "check state estimation" << endl;
					}
					c->compressedSequence.push_back(maxProbState);
					for (int c_id : c->children_id) {
						if ((*this->vertexMap)[c_id]->children_id.size()>0) {
							verticesToVisit.push_back((*this->vertexMap)[c_id]);
							numberOfVerticesToVisit += 1;
						}
					}					
				}
			}
			if (iter < 2 or currentLogLikelihood > previousLogLikelihood or abs(currentLogLikelihood-previousLogLikelihood) > 0.001){
				continueEM = 1;
				previousLogLikelihood = currentLogLikelihood;
			} else {
				continueEM = 0;
			}
		}
		this->logLikelihood = currentLogLikelihood;
}

void rootedPhylogeny_tree::EstimateAncestralSequencesByFittingTheGMMUsingDoublePrecision(){	
	this->ComputeMPEstimateOfAncestralSequences();
	map <int,Matrix4f> transitionMatrices;	
	map <int,array<double,4>> conditionalLikelihoodMap;
	array <double,4> conditionalLikelihood;	
	double partialLikelihood;
	double siteLikelihood;
	double currentLogLikelihood = 0;
	double previousLogLikelihood = 0;
	int iter = 0;
	int maxIter = 10;
	pair<int,int> edgeIds;
	unsigned char dna_p; unsigned char dna_c;
	char maxProbState;
	float rowSum;
	double currentProb;
	double maxProb;
	int numberOfVerticesToVisit;
	bool continueEM = 1;
	vector <rootedPhylogeny_vertex*> verticesToVisit;	
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f P;
	this->SetEdgesForPostOrderTreeTraversal();
	// Iterate till convergence of log likelihood
		while (continueEM and iter < maxIter) {
			iter += 1;			
//			cout << "root sequence is " << endl;
//			cout << EncodeAsDNA(this->root->compressedSequence) << endl;
			currentLogLikelihood = 0;
			// Estimate root probablity
			this->ComputeMLEOfRootProbability();		
			// Estimate transition matrices	
//			cout << "here 1" << endl;
			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
				c = idPtrPair.second;
				if (c->id != -1){
					p = (*this->vertexMap)[c->parent_id];
					P = ArrayXXf::Zero(4,4);
					for (unsigned int site = 0; site < siteWeights.size(); site++){
						dna_p = p->compressedSequence[site];
						dna_c = c->compressedSequence[site];
						P(dna_p,dna_c) += this->siteWeights[site];
					}
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
						rowSum = 0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							rowSum += P(dna_p,dna_c);
						}
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							 P(dna_p,dna_c) /= rowSum;
						}
					}
					transitionMatrices.insert(pair<int,Matrix4f>(c->id,P));
				}
			}
//			cout << "here 2" << endl;
			// Estimate ancestral sequences
			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
				c = idPtrPair.second;
				if (c->children_id.size() > 0){
					c->compressedSequence.clear();
				}
			}		
			// Iterate over sites		
			for (unsigned int site = 0 ; site < this->siteWeights.size(); site++){
				conditionalLikelihoodMap.clear();			
//				cout << "site is " << site << endl;
				for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){
					p = (*this->vertexMap)[edgeIds.first];
					c = (*this->vertexMap)[edgeIds.second];
					P = transitionMatrices[c->id];					
					// Initialize conditional likelihood for leaves
					if (c->children_id.size()==0){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
							conditionalLikelihood[dna_c] = 0;
						}
						conditionalLikelihood[c->compressedSequence[site]] = 1;
						conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
					}
					// Initialize conditional likelihood for ancestors
					if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
						conditionalLikelihood[dna_c] = 1;
						}				
						conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));					
					}		
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
						partialLikelihood = 0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						}
						conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
					}
				}
				maxProbState = -1;
				maxProb = 0;
				siteLikelihood = 0; 							
//				cout << prior << endl;
				for (int dna = 0; dna < 4; dna++) {
//					cout << "prior for dna:" << dna << " is " << prior[dna] << endl;
//					cout << "conditional for dna:" << dna << " is " << conditionalLikelihoodMap[-1][dna] << endl;
					currentProb = this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
					siteLikelihood += currentProb;
					if (currentProb > maxProb) {
						maxProb = currentProb;
						maxProbState = dna;		
					}
				}				
				currentLogLikelihood += log(siteLikelihood) * this->siteWeights[site];
				if (maxProbState == -1) {
					cout << "check state estimation" << endl;					
				}
//				for (int dna = 0; dna < 4; dna++){
//					cout << "conditional likelihood for dna " << dna << " is " << conditionalLikelihoodMap[-1][dna] << endl;
//				}
				(*this->vertexMap)[-1]->compressedSequence.push_back(maxProbState);
				verticesToVisit.clear();			
				for (int c_id: (*this->vertexMap)[-1]->children_id) {
					if ((*this->vertexMap)[c_id]->children_id.size()>0) {
						verticesToVisit.push_back((*this->vertexMap)[c_id]);
					}
				}
				numberOfVerticesToVisit = verticesToVisit.size();
				while (numberOfVerticesToVisit > 0) {
					c = verticesToVisit[numberOfVerticesToVisit-1];				
					verticesToVisit.pop_back();
					numberOfVerticesToVisit -= 1;
					p = (*this->vertexMap)[c->parent_id];
					P = transitionMatrices[c->id];
					dna_p = p->compressedSequence[site];
					maxProbState = -1;
					maxProb = 0;
					for (int dna_c = 0; dna_c < 4; dna_c++){
						currentProb = P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						if (currentProb > maxProb){
							maxProb = currentProb;
							maxProbState = dna_c;
						}					
					}
					if (maxProbState == -1){
						cout << "check state estimation" << endl;
					}
					c->compressedSequence.push_back(maxProbState);
					for (int c_id : c->children_id) {
						if ((*this->vertexMap)[c_id]->children_id.size()>0) {
							verticesToVisit.push_back((*this->vertexMap)[c_id]);
							numberOfVerticesToVisit += 1;
						}
					}					
				}
			}
//			cout << "iteration:" << iter << "\t";
//			cout << "previousLogLikelihood:" << previousLogLikelihood << endl;
//			cout << "currentLogLikelihood:" << currentLogLikelihood << endl;			
			if (iter < 2 or currentLogLikelihood > previousLogLikelihood or abs(currentLogLikelihood-previousLogLikelihood) > 0.001){
				continueEM = 1;
				previousLogLikelihood = currentLogLikelihood;
			} else {
				continueEM = 0;
			}
		}
		this->logLikelihood = currentLogLikelihood;
}


//void rootedPhylogeny_tree::EstimateAncestralSequencesByFittingTheGMMUsingMultiprecision(){
//	using namespace boost::multiprecision;
//	this->ComputeMPEstimateOfAncestralSequences();
//	map <int,Matrix4f> transitionMatrices;
//	array <float,4> prior;
//	map <int,array<cpp_dec_float_100,4>> conditionalLikelihoodMap;
//	array <cpp_dec_float_100,4> conditionalLikelihood;	
//	cpp_dec_float_100 partialLikelihood;
//	cpp_dec_float_100 siteLikelihood;
//	cpp_dec_float_100 currentLogLikelihood = 0;
//	cpp_dec_float_100 previousLogLikelihood = 0;
//	cpp_dec_float_100 maxProb;
//	cpp_dec_float_100 currentProb;
//	cpp_dec_float_100 maxValueForConditionalLikelihood;
////	map <int,array<double,4>> conditionalLikelihoodMap;
////	array <double,4> conditionalLikelihood;	
////	double partialLikelihood;
////	double siteLikelihood;
////	double currentLogLikelihood = 0;
////	double previousLogLikelihood = 0;
//	int iter = 0;
//	int maxIter = 10;
//	pair<int,int> edgeIds;
//	unsigned char dna_p; unsigned char dna_c;
//	char maxProbState;
//	float rowSum;
////	float currentProb;
////	float maxProb;
//	int numberOfVerticesToVisit;
//	bool continueEM = 1;
//	vector <rootedPhylogeny_vertex*> verticesToVisit;	
//	rootedPhylogeny_vertex * p;
//	rootedPhylogeny_vertex * c;
//	Matrix4f P;
//	this->SetEdgesForPostOrderTreeTraversal();
//	// Iterate till convergence of log likelihood
//		while (continueEM and iter < maxIter) {
//			iter += 1;			
////			cout << "root sequence is " << endl;
////			cout << EncodeAsDNA(this->root->compressedSequence) << endl;
//			currentLogLikelihood = 0;
//			// Estimate root probablity
//			this->ComputeMLEOfRootProbability();		
//			// Estimate transition matrices	
////			cout << "here 1" << endl;
//			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
//				c = idPtrPair.second;
//				if (c->id != -1){
//					p = (*this->vertexMap)[c->parent_id];
//					P = ArrayXXf::Zero(4,4);
//					for (unsigned int site = 0; site < siteWeights.size(); site++){
//						dna_p = p->compressedSequence[site];
//						dna_c = c->compressedSequence[site];
//						P(dna_p,dna_c) += this->siteWeights[site];
//					}
//					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
//						rowSum = 0;
//						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
//							rowSum += P(dna_p,dna_c);
//						}
//						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
//							 P(dna_p,dna_c) /= rowSum;
//						}
//					}
//					transitionMatrices.insert(pair<int,Matrix4f>(c->id,P));
//				}
//			}
////			cout << "here 2" << endl;
//			// Estimate ancestral sequences
//			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
//				c = idPtrPair.second;
//				if (c->children_id.size() > 0){
//					c->compressedSequence.clear();
//				}
//			}		
//			// Iterate over sites		
//			for (unsigned int site = 0 ; site < this->siteWeights.size(); site++){
//				conditionalLikelihoodMap.clear();			
////				cout << "site is " << site << endl;
//				for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){
//					p = (*this->vertexMap)[edgeIds.first];
//					c = (*this->vertexMap)[edgeIds.second];
//					P = transitionMatrices[c->id];					
//					maxValueForConditionalLikelihood = 0;
//					// Initialize conditional likelihood for leaves
//					if (c->children_id.size()==0){
//						for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
//							conditionalLikelihood[dna_c] = 0;
//						}
//						conditionalLikelihood[c->compressedSequence[site]] = 1;
//						conditionalLikelihoodMap.insert(pair<int,array<cpp_dec_float_100,4>>(c->id,conditionalLikelihood));
//					}
//					// Initialize conditional likelihood for ancestors
//					if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){
//						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
//						conditionalLikelihood[dna_c] = 1;
//						}				
//						conditionalLikelihoodMap.insert(pair<int,array<cpp_dec_float_100,4>>(p->id,conditionalLikelihood));					
//					}		
//					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
//						partialLikelihood = 0;
//						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
//							partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
//						}
//						conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
//						if (maxValueForConditionalLikelihood < conditionalLikelihoodMap[p->id][dna_p]){
//							maxValueForConditionalLikelihood = conditionalLikelihoodMap[p->id][dna_p];
//						}
//					}
//					if (maxValueForConditionalLikelihood < this->smallestMaxValueForConditionalLikelihoodOfAncestors){
//						this->smallestMaxValueForConditionalLikelihoodOfAncestors = maxValueForConditionalLikelihood;
//					}
//				}
//				maxProbState = -1;
//				maxProb = 0;
//				siteLikelihood = 0; 			
//				prior = this->rootProbability;
////				cout << prior << endl;
//				for (int dna = 0; dna < 4; dna++) {
////					cout << "prior for dna:" << dna << " is " << prior[dna] << endl;
////					cout << "conditional for dna:" << dna << " is " << conditionalLikelihoodMap[-1][dna] << endl;
//					currentProb = prior[dna]*conditionalLikelihoodMap[-1][dna];
//					siteLikelihood += currentProb;
//					if (currentProb > maxProb) {
//						maxProb = currentProb;
//						maxProbState = dna;		
//					} else {
////						cout << currentProb << endl;
//					}
//				}				
//				currentLogLikelihood += log(siteLikelihood) * this->siteWeights[site];
//				if (maxProbState == -1) {
//					cout << "check state estimation" << endl;
//				}
////				if (site < 5){
////					cout << "site is " << site << endl;
////					for (int dna = 0; dna < 4; dna++){
////						cout << "conditional likelihood for root for dna " << dna << " is " << conditionalLikelihoodMap[-1][dna] << endl;
////					}
////				}
//				(*this->vertexMap)[-1]->compressedSequence.push_back(maxProbState);
//				verticesToVisit.clear();			
//				for (int c_id: (*this->vertexMap)[-1]->children_id) {
//					if ((*this->vertexMap)[c_id]->children_id.size()>0) {
//						verticesToVisit.push_back((*this->vertexMap)[c_id]);
//					}
//				}
//				numberOfVerticesToVisit = verticesToVisit.size();
//				while (numberOfVerticesToVisit > 0) {
//					c = verticesToVisit[numberOfVerticesToVisit-1];				
//					verticesToVisit.pop_back();
//					numberOfVerticesToVisit -= 1;
//					p = (*this->vertexMap)[c->parent_id];
//					P = transitionMatrices[c->id];
//					dna_p = p->compressedSequence[site];
//					maxProbState = -1;
//					maxProb = 0;
//					for (int dna_c = 0; dna_c < 4; dna_c++){
//						currentProb = P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
//						if (currentProb > maxProb){
//							maxProb = currentProb;
//							maxProbState = dna_c;
//						}					
//					}
//					if (maxProbState == -1){
//						cout << "check state estimation" << endl;
//					}
//					c->compressedSequence.push_back(maxProbState);
//					for (int c_id : c->children_id) {
//						if ((*this->vertexMap)[c_id]->children_id.size()>0) {
//							verticesToVisit.push_back((*this->vertexMap)[c_id]);
//							numberOfVerticesToVisit += 1;
//						}
//					}					
//				}
//			}
////			cout << "iteration:" << iter << "\t";
////			cout << "previousLogLikelihood:" << previousLogLikelihood << endl;
////			cout << "currentLogLikelihood:" << currentLogLikelihood << endl;			
//			if (iter < 2 or currentLogLikelihood > previousLogLikelihood or abs(currentLogLikelihood-previousLogLikelihood) > 0.001){
//				continueEM = 1;
//				previousLogLikelihood = currentLogLikelihood;
//			} else {
//				continueEM = 0;
//			}
//		}
//		if (1){
//			cout << "smallest max value for conditional likelihood is " << this->smallestMaxValueForConditionalLikelihoodOfAncestors << endl;
//		}		
//		this->logLikelihood = currentLogLikelihood.convert_to<double>();
////		cout << "max log likelihood for fitting the GMM is " << currentLogLikelihood << endl;
//}

void rootedPhylogeny_tree::ComputeMLEOfEdgeLengths(){
	Matrix4f Q = this->rateMatrix;
	Matrix4f Q_scaled;
	Matrix4f P;
	Matrix4f P_der_1;
	Matrix4f P_der_2;
	float t_new;
	float t_old;
//	float length_reduction;
	float proposedChange;
	float s = 1;
	double firstDerivative;
	double secondDerivative;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
//	double edgeLogLik_old;
	double edgeLogLik_current;	
	double edgeLogLik_updated;
	unsigned char dna_p;
	unsigned char dna_c;	
	this->SetMinLengthOfEdges();
	bool continueOptimizingEdgeLength;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t_old = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t_old;
			P = Q_scaled.exp();
//			edgeLogLik_old = this->ComputeEdgeLogLikelihood(p,c,P);			
			continueOptimizingEdgeLength = 1;
			while (continueOptimizingEdgeLength){
				s = 1.0;
				Q_scaled = Q*t_old;
				P = Q_scaled.exp();
				P_der_1 = Q*P;
				P_der_2 = Q*Q;
				P_der_2 = P_der_2*P;
				firstDerivative = 0;
				secondDerivative = 0;
				edgeLogLik_current = this->ComputeEdgeLogLikelihood(p,c,P);
				for (unsigned int site = 0; site < siteWeights.size(); site++){
					dna_p = p->compressedSequence[site];
					dna_c = c->compressedSequence[site];
					firstDerivative += (P_der_1(dna_p,dna_c)/P(dna_p,dna_c))*(this->siteWeights[site]);
					secondDerivative += (P_der_2(dna_p,dna_c)/P(dna_p,dna_c) - pow(P_der_1(dna_p,dna_c)/P(dna_p,dna_c),2.0))*(this->siteWeights[site]);
				}				
				proposedChange = firstDerivative/secondDerivative; 				
				while (t_old - s*proposedChange < 0){
					s/=2.0;					
				}
				t_new = t_old - s*proposedChange;
				Q_scaled = Q*t_new;
				P = Q_scaled.exp();
				edgeLogLik_updated = this->ComputeEdgeLogLikelihood(p,c,P);
				while (edgeLogLik_updated < edgeLogLik_current and abs(edgeLogLik_updated - edgeLogLik_current) > 0.001) {
					s/=2.0;					
					t_new = t_old - s*proposedChange;
					Q_scaled = Q*t_new;
					P = Q_scaled.exp();
					edgeLogLik_updated = this->ComputeEdgeLogLikelihood(p,c,P); 
				}
				if (abs(edgeLogLik_current-edgeLogLik_updated) < 0.001){
					continueOptimizingEdgeLength = 0;
				}
				t_old = t_new;
			}
			if (p->id < c->id) {
				(*this->edgeLengthsMap)[pair<int,int>(p->id,c->id)] = t_old;
			} else {
				(*this->edgeLengthsMap)[pair<int,int>(c->id,p->id)] = t_old;
			}
		}
	}
}

void rootedPhylogeny_tree::ComputeMLEOfRateMatricesForLeafLabeledTrees(){	
	
	//	this->ComputeJacobianOfLogLikelihood();
//	int bits = std::numeric_limits<double>::digits;
//	using namespace brent;
//	double t = 0.001;
//	float scalingFactor;
//	double a = -4.0;
//	double b = 4.0/3.0;	
//	int status = 0;
//	double value;
//	double arg;
//	double result;
//	value = this->SampleFunctionForMinimization(arg);
//	result = local_min (a, b, t, this->SampleFunctionForMinimization, arg);
//	cout << "result is " << result << endl;
//	cout << "arg is " << arg << endl;
//	arg = local_min (a, b, status, this->SampleFunctionForMinimization());
//	cout << arg << endl;
//	MatrixXd searchDirection = ArrayXXd::Zero(11,1);
//	MatrixXd B_old = this->GetFreeParameters(this->rateMatrix);
//	for (int par_1 = 0; par_1 < 11; par_1++){
//		for (int par_2 = 0; par_2 < 11; par_2++){
//			this->Hessian(par_1,par_2) = 1;
//		}
//	}
//	searchDirection = -1*(this->Hessian.inverse()*this->Jacobian);
//	pair <double, double> r2 = brent_find_minima((*this->operator()), -4.0, 4.0/3.0, bits);
//	cout << "x at minima is " << r2.first << endl;
	// select scalingFactor using brent's algorithm
//	searchDirection = -1*this->Hessian.inverse()*this->Jacobian;
	
}


void rootedPhylogeny_tree::ScaleEdgeLengths(){ // Modify to include rate categories
	float scalingFactor = this->ComputeScalingFactor(this->rateMatrix);
	for (pair<pair<int,int>,float> edgeIdsAndLengthsPair: *this->edgeLengthsMap){
		(*this->edgeLengthsMap)[edgeIdsAndLengthsPair.first] = edgeIdsAndLengthsPair.second/scalingFactor;
	}
}


void rootedPhylogeny_tree::ComputeMAPEstimateOfAncestralSequencesUsingLongDouble(){
	map <int, array<long double,4>> conditionalLikelihoodMap;
	array <long double, 4> conditionalLikelihood;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	int maxProbState;
	long double maxProb;
	long double partialLikelihood;	
	long double siteLikelihood;
	long double currentProb;
	unsigned char dna_p;
	vector <rootedPhylogeny_vertex *> verticesToVisit;
	// Clear ancestral sequences
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->children_id.size() > 0){
			c->compressedSequence.clear();
		}
	}
	Matrix4f Q;
	Matrix4f Q_scaled;
	Matrix4f P;
	float t;
	float scalingFactor;
	int numberOfVerticesToVisit;
	double currentLogLikelihood = 0;
	// Iterate over sites		
	for (unsigned int site = 0 ; site < this->siteWeights.size(); site++) {
		conditionalLikelihoodMap.clear();			
		for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edgeIds.first];
			c = (*this->vertexMap)[edgeIds.second];
			t = this->GetEdgeLength(p->id,c->id);
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
			Q_scaled = Q*(t/scalingFactor);			
			P = Q_scaled.exp();			
			// Initialize conditional likelihood for leaves
			if (c->children_id.size()==0) {
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
					conditionalLikelihood[dna_c] = 0;
				}
				conditionalLikelihood[c->compressedSequence[site]] = 1;
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for ancestors
			if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){
				for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
				conditionalLikelihood[dna_c] = 1;
				}				
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));					
			}		
			for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				}
				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
			}
		}
		maxProbState = -1;
		maxProb = 0;
		siteLikelihood = 0; 			
		for (int dna = 0; dna < 4; dna++) {
			currentProb = this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
			siteLikelihood += currentProb;
			if (currentProb > maxProb){
				maxProb = currentProb;
				maxProbState = dna;		
			}
		}
		currentLogLikelihood += log(siteLikelihood) * this->siteWeights[site];
		if (maxProbState == -1){
			cout << "check state estimation" << endl;
		}			
		(*this->vertexMap)[-1]->compressedSequence.push_back(maxProbState);
		verticesToVisit.clear();			
		for (int c_id: (*this->vertexMap)[-1]->children_id) {
			if ((*this->vertexMap)[c_id]->children_id.size() > 0) {
				verticesToVisit.push_back((*this->vertexMap)[c_id]);
			}
		}
		numberOfVerticesToVisit = verticesToVisit.size();
		while (numberOfVerticesToVisit > 0) {
			c = verticesToVisit[numberOfVerticesToVisit-1];
			verticesToVisit.pop_back();
			numberOfVerticesToVisit -= 1;
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id,c->id);
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
			Q_scaled = Q*(t/scalingFactor);
			P = Q_scaled.exp();
			dna_p = p->compressedSequence[site];
			maxProbState = -1;
			maxProb = 0;
			for (int dna_c = 0; dna_c < 4; dna_c++){
				currentProb = P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				if (currentProb > maxProb){
					maxProb = currentProb;
					maxProbState = dna_c;
				}
			}
			if (maxProbState == -1){
				cout << "check state estimation" << endl;
			}
			c->compressedSequence.push_back(maxProbState);
			for (int c_id : c->children_id) {
				if ((*this->vertexMap)[c_id]->children_id.size()>0) {
					verticesToVisit.push_back((*this->vertexMap)[c_id]);
					numberOfVerticesToVisit += 1;
				}
			}
		}
	}
	this->logLikelihood = currentLogLikelihood;
}


Matrix4f rootedPhylogeny_tree::ComputeFirstDerivativeOfRateMatrix(int rateCat, int par){
	Matrix4f Q_der;	
	MatrixXd B = (*this->parametersPerRateCategory)[rateCat];
	
	float pi_1; float pi_2; float pi_3; float pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;
	
	pi_1 = B(0,0); pi_2 = B(1,0); pi_3 = B(2,0);
	pi_4 = (1 - pi_1 - pi_2 - pi_3);
	
	a = B(3,0); b = B(4,0); c = B(5,0); d = B(6,0); e = B(7,0);
	f = B(8,0); g = B(9,0); h = B(10,0);
		  
	float a_dot = 0; float b_dot = 0; float c_dot = 0; float d_dot = 0; 
	float e_dot = 0; float f_dot = 0; float g_dot = 0; float h_dot = 0; 
	float i_dot = 0; float j_dot = 0; float k_dot = 0; float l_dot = 0; 
	float D1_dot = 0; float D2_dot = 0; float D3_dot = 0; float D4_dot = 0; 	
	
	if (par == 0){
		// par is pi_1
		i_dot = -(a+b+2*c)/(2*pi_3);
		j_dot = (a+b+c)/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = -a/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/pow(pi_4,2);
		l_dot = -(a+2*c+3*b)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e) )/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (a+b+2*c)/(2*pi_3);
		D4_dot = (a+b)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 1){
		// par is pi_2
		i_dot = -(d+e+2*f)/(2*pi_3);
		j_dot = -d/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = (d+e+f)/(pi_4) + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = -(d+2*f+3*e)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (d+e+2*f)/(2*pi_3);
		D4_dot = (d+e)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 2){
		// par is pi_3
		i_dot = -(g+h)/(2*pi_3) + (pi_1*(a+b+2*c) + pi_2*(d+e+2*f) + pi_3*(g+h) -1)/(2*pow(pi_3,2));
		j_dot = -g/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/pow(pi_4,2);
		k_dot = -h/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = (g+h)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = -(g+h)/(2*pi_3) -(pi_1*(a+b+2*c) + pi_2*(d+e+2*f) - pi_3*(g+h) -1)/(2*pow(pi_3,2));
		D4_dot = (g+h)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 3){
		// par is a
		a_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = -pi_1/pi_4;
		l_dot = -pi_1/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);		
	} else if (par == 4){
		// par is b
		b_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = 0;
		l_dot = -(3*pi_1)/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);
	} else if (par == 5){
		// par is c
		c_dot = 1.0;
		i_dot = -pi_1/pi_3;
		j_dot = pi_1/pi_4;
		k_dot = 0.0;
		l_dot = -pi_1/pi_4;
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/pi_3;
		D4_dot = 0.0;
	} else if (par == 6){
		// par is d
		d_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = -pi_2/pi_4;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);
	} else if (par == 7){
		// par is e
		e_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = 0.0;
		k_dot = pi_2/pi_4;
		l_dot = -(3*pi_2)/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);	
	} else if (par == 8){
		// par is f
		f_dot = 1.0;
		i_dot = -pi_2/pi_3;
		j_dot = 0;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/pi_4;
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/pi_3;
		D4_dot = 0.0;
	} else if (par == 9){
		// par is g
		g_dot = 1.0;
		i_dot = -0.5;
		j_dot = -pi_3/pi_4;
		k_dot = 0.0;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	} else if (par == 10){
		// par is h
		h_dot = 1.0;
		i_dot = -0.5;
		j_dot = 0.0;
		k_dot = -pi_3/pi_4;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	}	
	Q_der << D1_dot, a_dot, b_dot, c_dot,
			 d_dot, D2_dot, e_dot, f_dot,
			 g_dot, h_dot, D3_dot, i_dot,
			 j_dot, k_dot, l_dot, D4_dot;	
	return (Q_der);
}

Matrix4f rootedPhylogeny_tree::ComputeFirstDerivativeOfRateMatrixDep(int par){
	Matrix4f Q_der = ArrayXXf::Zero(4,4);
	if (par == 0){
		// par is a
		Q_der(0,0) = -1;
		Q_der(0,1) = 1;
		
	} else if (par == 1){
		// par is b
		Q_der(0,0) = -1;
		Q_der(0,2) = 1;
		
	} else if (par == 2){
		// par is c
		Q_der(0,0) = -1;
		Q_der(0,3) = 1;
		
	} else if (par == 3){
		// par is d
		Q_der(1,1) = -1;
		Q_der(1,0) = 1;
		
	} else if (par == 4){
		// par is e
		Q_der(1,1) = -1;
		Q_der(1,2) = 1;
		
	} else if (par == 5){
		// par is f
		Q_der(1,1) = -1;
		Q_der(1,3) = 1;
		
	} else if (par == 6){
		// par is g
		Q_der(2,2) = -1;
		Q_der(2,0) = 1;
		
	} else if (par == 7){
		// par is h
		Q_der(2,2) = -1;
		Q_der(2,1) = 1;
		
	} else if (par == 8){
		// par is i
		Q_der(2,2) = -1;
		Q_der(2,3) = 1;
		
	} else if (par == 9){
		// par is j
		Q_der(3,3) = -1;
		Q_der(3,0) = 1;
		
	} else if (par == 10){
		// par is k
		Q_der(3,3) = -1;
		Q_der(3,1) = 1;		
	}
	return Q_der;
}

Matrix4f rootedPhylogeny_tree::ComputeFirstDerivativeOfMatrixExponential(float t, int rateCat, int par){
	Matrix4f Q; Matrix4f Q_der; MatrixXf Q_aug; MatrixXf Q_aug_scaled;
	MatrixXf P_aug; Matrix4f P_der;
	Q = (*this->rateMatrixPerRateCategory)[rateCat];
	Q_der = this->ComputeFirstDerivativeOfRateMatrix(rateCat, par);
	Q_aug = ArrayXXf::Zero(8,8);
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			Q_aug(row,col) = Q(row,col);
			Q_aug(row+4,col+4) = Q(row,col);
			Q_aug(row,col+4) = Q_der(row,col);
		}
	}
	Q_aug_scaled = Q_aug*t;
	P_aug = Q_aug.exp();
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			P_der(row,col) = P_aug(row,col+4);
		}
	}
	return P_der;
}

Matrix4f rootedPhylogeny_tree::ComputeSecondDerivativeOfMatrixExponential(float t, int par_1, int par_2){
	Matrix4f Q; Matrix4f Q_der_1; Matrix4f Q_der_2;
	Matrix4f Q_der_1_2; Matrix4f P_der;
	MatrixXf Q_aug; MatrixXf P_aug;
	Q_der_1_2 = ArrayXXf::Zero(4,4);
	Q = this->rateMatrix;
	Q_der_1 = this->ComputeFirstDerivativeOfRateMatrixDep(par_1);
	Q_der_2 = this->ComputeFirstDerivativeOfRateMatrixDep(par_2);
	Q_aug = ArrayXXf::Zero(16,16);
	for (int row = 0; row < 4; row ++){
		for (int col = 0; col < 4; col++){
			Q_aug(row,col) = Q(row,col);
			Q_aug(row+4,col+4) = Q(row,col);
			Q_aug(row+8,col+8) = Q(row,col);
			Q_aug(row+12,col+12) = Q(row,col);
			Q_aug(row,col+4) = Q_der_1(row,col);
			Q_aug(row+8,col+12) = Q_der_1(row,col);
			Q_aug(row,col+8) = Q_der_2(row,col);
			Q_aug(row+4,col+12) = Q_der_2(row,col);
			Q_aug(row,col+12) = Q_der_1_2(row,col);
		}
	}
	Q_aug = Q_aug*t;
	P_aug = Q_aug.exp();
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			P_der(row,col) = P_aug(row,col+12);
		}
	}
	return P_der;
}

MatrixXd rootedPhylogeny_tree::GetFreeParametersDep(Matrix4f Q){
	MatrixXd B;
	B = ArrayXXd::Zero(11,1);
	B(0,0) = Q(0,1);
	B(1,0) = Q(0,2);
	B(2,0) = Q(0,3);
	B(3,0) = Q(1,0);
	B(4,0) = Q(1,2);
	B(5,0) = Q(1,3);
	B(6,0) = Q(2,0);
	B(7,0) = Q(2,1);
	B(8,0) = Q(2,3);
	B(9,0) = Q(3,0);
	B(10,0) = Q(3,1);
	return B;
}

void rootedPhylogeny_tree::SetRateMatrixUsingFreeParametersDep(MatrixXd B){
	Matrix4f Q;
	Q(0,1) = B(0,0);
	Q(0,2) = B(1,0);
	Q(0,3) = B(2,0);
	Q(0,0) = -(Q(0,1) + Q(0,2) + Q(0,3));
	Q(1,0) = B(3,0);
	Q(1,2) = B(4,0);
	Q(1,3) = B(5,0);
	Q(1,1) = -(Q(1,0) + Q(1,2) + Q(1,3));
	Q(2,0) = B(6,0);
	Q(2,1) = B(7,0);
	Q(2,3) = B(8,0);
	Q(2,2) = -(Q(2,0) + Q(2,1) + Q(2,3));
	Q(3,0) = B(9,0);
	Q(3,1) = B(10,0);
	Q(3,2) = 1;
	Q(3,3) = -(Q(3,0) + Q(3,1) + 1);
	this->rateMatrix = Q;
}

void rootedPhylogeny_tree::ComputeMLEOfRateMatrices(){
	MatrixXd B_new;
	MatrixXd JacobianForRateMatrix; MatrixXd HessianDep; MatrixXd proposedChange;
	proposedChange = ArrayXXd::Zero(11,1);
	JacobianForRateMatrix = ArrayXXd::Zero(11,1);
	HessianDep = ArrayXXd::Zero(11,11);
	B_new = ArrayXXd::Zero(11,1);	
	double logLikelihood_updated;
	double logLikelihood_current;
	bool continueOptimization = 1;
	float s;	
	MatrixXd B_old = this->GetFreeParametersDep(this->rateMatrix);
	while (continueOptimization){
		this->ComputeLogLikelihoodForFullyLabeledTree();
		logLikelihood_current = this->logLikelihood;
		this->ComputeJacobianForFullyLabeledTree();
		this->ComputeHessianForRateMatrixParametersDep();
		proposedChange = this->HessianDep.inverse()*this->JacobianForRateMatrix;
		// Improve step size selection
		s = 1.0;
		for (int par = 0; par < 11; par++){
			while (B_old(par,0) - proposedChange(par,0)*s < 0){
				s /= 2.0;
			}
		}
		B_new = B_old - proposedChange*s;
		this->SetRateMatrixUsingFreeParametersDep(B_new);
		this->ComputeLogLikelihoodForFullyLabeledTree();
		logLikelihood_updated = this->logLikelihood;
		while (logLikelihood_updated < logLikelihood_current and abs(logLikelihood_updated - logLikelihood_current) > 0.001){
			s /= 2.0;
			B_new = B_old - proposedChange*s;
			this->SetRateMatrixUsingFreeParametersDep(B_new);
			this->ComputeLogLikelihoodForFullyLabeledTree();
			logLikelihood_updated = this->logLikelihood;
		}
		B_new = B_old - proposedChange*s;
		B_old = B_new;
		this->SetRateMatrixUsingFreeParametersDep(B_old);
		if (abs(logLikelihood_updated - logLikelihood_current) < 0.001) {
			continueOptimization = 0;
		}
	}
}

void rootedPhylogeny_tree::ComputeHessianForRateMatrixParametersDep(){
	Matrix4f P; Matrix4f P_der_1; Matrix4f P_der_2; Matrix4f P_der_1_2;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	int rateCat = 0;
	for (int par_1 = 0; par_1 < 11; par_1 ++){
		for (int par_2 = 0; par_2 < 11; par_2++){
			this->HessianDep(par_1,par_2) = 0;
		}
	}
	Matrix4f Q = this->rateMatrix;
	Matrix4f Q_scaled;
	float t;
	unsigned char dna_p;
	unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			for (int par_1 = 0; par_1 < 11; par_1++){
				P_der_1 = this->ComputeFirstDerivativeOfMatrixExponential(t, rateCat, par_1);
				for (int par_2 = 0; par_2 < 11; par_2++){
					P_der_2 = this->ComputeFirstDerivativeOfMatrixExponential(t, rateCat, par_2);
					P_der_1_2 = this->ComputeSecondDerivativeOfMatrixExponential(t, par_1, par_2);
					for (unsigned int site = 0; site < this->siteWeights.size(); site++){
						dna_p = p->compressedSequence[site];
						dna_c = c->compressedSequence[site];
						this->HessianDep(par_1,par_2) += ((P_der_1_2(dna_p,dna_c)/P(dna_p,dna_c)) - (P_der_1(dna_p,dna_c)*P_der_2(dna_p,dna_c))/pow(P(dna_p,dna_c),2))*this->siteWeights[site];
					}
				}
			}
		}
	}	
}

void rootedPhylogeny_tree::ComputeJacobianForFullyLabeledTree(){
	for (int par = 0; par < 11; par ++){
		this->JacobianForRateMatrix(par,0) = 0;
	}
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;	
	Matrix4f P; Matrix4f Q; Matrix4f Q_scaled; Matrix4f P_der;
	int rateCat = 0;
	Q = this->rateMatrix;	
	float t;
	unsigned char dna_p;
	unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			for (int par = 0; par < 11; par++){
				P_der = this->ComputeFirstDerivativeOfMatrixExponential(t, rateCat, par);					
				for (unsigned int site = 0; site < this->siteWeights.size(); site++){
					dna_p = p->compressedSequence[site];
					dna_c = c->compressedSequence[site];
					this->JacobianForRateMatrix(par,0) += (this->siteWeights[site]*P_der(dna_p,dna_c))/P(dna_p,dna_c);
				}
			}				
		}
	}
}

void rootedPhylogeny_tree::ComputeMLEOfRootProbability(){
	array <float,4> prob;
	for (int dna = 0; dna < 4; dna++){
		prob[dna] = 0;
	} 
	rootedPhylogeny_vertex * r = (*this->vertexMap)[-1];
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		prob[r->compressedSequence[site]] += this->siteWeights[site];		
	}
	float sum = 0;
	for (int dna = 0; dna < 4; dna++){
		sum += prob[dna];
	}
	for (int dna = 0; dna < 4; dna++){
		prob[dna] /= sum;
	}
	this->rootProbability = prob;
}



Matrix4f rootedPhylogeny_tree::ComputeDerivativeOfMatrixExponentialDep(float t, int par){
	VectorXd x(11);
	x = this->freeParametersExcBaseFreq;
	float pi_1; float pi_2; float pi_3; float pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;float i; float j; float k; float l;
	pi_1 = x[0]; pi_2 = x[1]; pi_3 = x[2];
	pi_4 = (1 - pi_1 - pi_2 - pi_3);
	// Construct Q
	a = x[3]; b = x[4]; c = x[5]; d = x[6]; e = x[7];
	f = x[8]; g = x[9]; h = x[10];

	i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3);
	j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4;
	k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4;
	l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4);
	Matrix4f Q;
	Q << -(a+b+c), a, b, c,
		  d, -(d+e+f), e, f,
		  g, h, -(g+h+i), i,
		  j, k, l, -(j+k+l);
		  
	float a_dot = 0; float b_dot = 0; float c_dot = 0; float d_dot = 0; 
	float e_dot = 0; float f_dot = 0; float g_dot = 0; float h_dot = 0; 
	float i_dot = 0; float j_dot = 0; float k_dot = 0; float l_dot = 0; 
	float D1_dot = 0; float D2_dot = 0; float D3_dot = 0; float D4_dot = 0; 	
	if (par == 0){
		// par is pi_1
		i_dot = -(a+b+2*c)/(2*pi_3);
		j_dot = (a+b+c)/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = -a/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/pow(pi_4,2);
		l_dot = -(a+2*c+3*b)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e) )/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (a+b+2*c)/(2*pi_3);
		D4_dot = (a+b)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 1){
		// par is pi_2
		i_dot = -(d+e+2*f)/(2*pi_3);
		j_dot = -d/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = (d+e+f)/(pi_4) + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = -(d+2*f+3*e)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (d+e+2*f)/(2*pi_3);
		D4_dot = (d+e)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 2){
		// par is pi_3
		i_dot = -(g+h)/(2*pi_3) + (pi_1*(a+b+2*c) + pi_2*(d+e+2*f) + pi_3*(g+h) -1)/(2*pow(pi_3,2));
		j_dot = -g/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/pow(pi_4,2);
		k_dot = -h/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = (g+h)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = -(g+h)/(2*pi_3) -(pi_1*(a+b+2*c) + pi_2*(d+e+2*f) - pi_3*(g+h) -1)/(2*pow(pi_3,2));
		D4_dot = (g+h)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 3){
		// par is a
		a_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = -pi_1/pi_4;
		l_dot = -pi_1/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);		
	} else if (par == 4){
		// par is b
		b_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = 0;
		l_dot = -(3*pi_1)/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);
	} else if (par == 5){
		// par is c
		c_dot = 1.0;
		i_dot = -pi_1/pi_3;
		j_dot = pi_1/pi_4;
		k_dot = 0.0;
		l_dot = -pi_1/pi_4;
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/pi_3;
		D4_dot = 0.0;
	} else if (par == 6){
		// par is d
		d_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = -pi_2/pi_4;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);
	} else if (par == 7){
		// par is e
		e_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = 0.0;
		k_dot = pi_2/pi_4;
		l_dot = -(3*pi_2)/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);	
	} else if (par == 8){
		// par is f
		f_dot = 1.0;
		i_dot = -pi_2/pi_3;
		j_dot = 0;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/pi_4;
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/pi_3;
		D4_dot = 0.0;
	} else if (par == 9){
		// par is g
		g_dot = 1.0;
		i_dot = -0.5;
		j_dot = -pi_3/pi_4;
		k_dot = 0.0;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	} else if (par == 10){
		// par is h
		h_dot = 1.0;
		i_dot = -0.5;
		j_dot = 0.0;
		k_dot = -pi_3/pi_4;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	}
	// Construct Q_der
	Matrix4f Q_der;
	Q_der << D1_dot, a_dot, b_dot, c_dot,
			 d_dot, D2_dot, e_dot, f_dot,
			 g_dot, h_dot, D3_dot, i_dot,
			 j_dot, k_dot, l_dot, D4_dot;
	// Construct Q_aug		 
	MatrixXf Q_aug;		
	Q_aug = ArrayXXf::Zero(8,8);			 
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			Q_aug(row,col) = Q(row,col);
			Q_aug(row+4,col+4) = Q(row,col);
			Q_aug(row,col+4) = Q_der(row,col);
		}
	}	
	MatrixXf Q_aug_scaled = Q_aug*t;
	MatrixXf P_aug = Q_aug_scaled.exp();
	Matrix4f P_der;
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			P_der(row,col) = P_aug(row,col+4);
		}
	}
	return P_der;
}

void rootedPhylogeny_tree::SetParameters(VectorXd x){
	for (int i = 0; i < 11; i++){
		this->freeParametersExcBaseFreq[i] = x[i];
	}	
	float pi_1; float pi_2; float pi_3; float pi_4;
	pi_1 = x[0]; pi_2 = x[1]; pi_3 = x[2]; pi_4 = 1-(pi_1+pi_2+pi_3);
	this->rootProbability[0] = pi_1;
	this->rootProbability[1] = pi_2;
	this->rootProbability[2] = pi_3;
	// check if x[0] + x[1] + x[2] < 1
	if (pi_4 < 0){
		cout << "check free parameters of stationary distribution" << endl;
	}
	this->rootProbability[3] = pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;float i; float j; float k; float l;
	a = x[3]; b = x[4]; c = x[5]; d = x[6]; e = x[7];
	f = x[8]; g = x[9]; h = x[10];
	i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3);
	j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4;
	k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4;
	l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4);
	this->rateMatrix << -(a+b+c), a, b, c,
		  d, -(d+e+f), e, f,
		  g, h, -(g+h+i), i,
		  j, k, l, -(j+k+l);
}

float rootedPhylogeny_tree::ComputeScalingFactor(Matrix4f Q){
	MatrixXf Q_aug = ArrayXXf::Zero(4,5);
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			Q_aug(row, col) = Q(row, col);
		}
	}
	for (int row = 0; row < 4; row++){
		Q_aug(row, 4) = 1;
	}	
	MatrixXf b = ArrayXXf::Zero(5,1);
	for (int row = 0; row < 4; row++){
		b(row,0) = 0;
	}
	b(4,0) = 1;	
	MatrixXf pi = ArrayXXf::Zero(1,4);
	pi = Q_aug.transpose().colPivHouseholderQr().solve(b).transpose();
	float scalingFactor = 0;
	for (int i = 0; i < 4; i++){
		scalingFactor -= pi(0,i) * Q(i,i);
	}
	return scalingFactor;
}

MatrixXf rootedPhylogeny_tree::ComputeStationaryDistribution(Matrix4f Q){
	MatrixXf Q_aug = ArrayXXf::Zero(4,5);
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			Q_aug(row, col) = Q(row, col);
		}
	}
	for (int row = 0; row < 4; row++){
		Q_aug(row, 4) = 1;
	}	
	MatrixXf b = ArrayXXf::Zero(5,1);
	for (int row = 0; row < 4; row++){
		b(row,0) = 0;
	}
	b(4,0) = 1;	
	MatrixXf pi = ArrayXXf::Zero(4,1);
	pi = Q_aug.transpose().colPivHouseholderQr().solve(b);
	return pi;	
}

void rootedPhylogeny_tree::ComputeMPEstimateOfAncestralSequences(){
	vector <int> preOrderVerticesWithoutLeaves_ids;
	vector <int> leafIds;
	preOrderVerticesWithoutLeaves_ids.push_back(-1);
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	rootedPhylogeny_vertex * v;
	verticesToVisit.push_back(this->root);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		v->compressedSequence.clear();
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (int c_id : v->children_id){
			if ((*this->vertexMap)[c_id]->children_id.size() > 0){
				preOrderVerticesWithoutLeaves_ids.push_back(c_id);
				verticesToVisit.push_back((*this->vertexMap)[c_id]);
				numberOfVerticesToVisit += 1;
			} else {
				leafIds.push_back(c_id);
			}
		}
	}
	int p_id; int numberOfPossibleStates; int pos;
	map <int,vector<unsigned char>> VU;
	map <int,unsigned char> V;
	for (unsigned int site = 0; site < siteWeights.size(); site++){
		VU.clear(); V.clear();
		// Set VU and V for leaves;
		for (int v_id : leafIds){
			V[v_id] = (*this->vertexMap)[v_id]->compressedSequence[site];
			VU[v_id].push_back((*this->vertexMap)[v_id]->compressedSequence[site]);
		}
		// Set VU for ancestors
		for (int v_ind = preOrderVerticesWithoutLeaves_ids.size()-1; v_ind > -1; v_ind--){
			p_id = preOrderVerticesWithoutLeaves_ids[v_ind];
			map <unsigned char, int> dnaCount;
			for (unsigned char dna = 0; dna < 4; dna++){
				dnaCount[dna] = 0;
			}
			for (int c_id : (*this->vertexMap)[p_id]->children_id){
				for (unsigned char dna: VU[c_id]){
					dnaCount[dna] += 1;
				}
			}
			int maxCount = 0;
			for (pair<unsigned char, int> dnaCountPair: dnaCount){
				if (dnaCountPair.second > maxCount){
					maxCount = dnaCountPair.second;
				}
			}			
			for (pair<unsigned char, int> dnaCountPair: dnaCount){
				if (dnaCountPair.second == maxCount){
					VU[p_id].push_back(dnaCountPair.first);					
				}
			}			
		}
		// Set V for ancestors
		for (int v_id : preOrderVerticesWithoutLeaves_ids){
			if (v_id == -1){
			// Set V for root
				if (VU[-1].size()==1){
					V[-1] = VU[-1][0];
				} else {
					numberOfPossibleStates = VU[-1].size();
					uniform_int_distribution<int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V[-1] = VU[-1][pos];
				}				
			} else {
				p_id = (*this->vertexMap)[v_id]->parent_id;
				if (find(VU[v_id].begin(),VU[v_id].end(),V[p_id])==VU[v_id].end()){
					numberOfPossibleStates = VU[v_id].size();
					uniform_int_distribution<int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V[v_id] = VU[v_id][pos];					
				} else {
					V[v_id] = V[p_id];
				}				
			}
			// push states to compressedSequence
			(*this->vertexMap)[v_id]->compressedSequence.push_back(V[v_id]);
		}		
	}	
}

array <double, 4> rootedPhylogeny_tree::GetLikelihoodArray(double elem_0, double elem_1, double elem_2, double elem_3){
	array <double, 4> likelihoodArray;		
	likelihoodArray[0] = elem_0;
	likelihoodArray[1] = elem_1;
	likelihoodArray[2] = elem_2;
	likelihoodArray[3] = elem_3;
	return likelihoodArray;
}

void rootedPhylogeny_tree::ComputeInitialEstimateForRateMatrix(){
	Matrix4f P = ArrayXXf::Zero(4, 4);	
	Matrix4f Q;	
	int dna_p; int dna_c;	
	int rateCategory = 1;	
	rootedPhylogeny_vertex * p = (*this->vertexMap)[-1];	    
	rootedPhylogeny_vertex * c;	
	// Select change-point for rateCategory
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	verticesToVisit.push_back(p);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		p = verticesToVisit[numberOfVerticesToVisit-1];		
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (int c_id: p->children_id){
			c = (*this->vertexMap)[c_id];
			if (c->rateCategory == rateCategory or 1){
				verticesToVisit.push_back(c);
				numberOfVerticesToVisit += 1;
				for (unsigned int site = 0; site < siteWeights.size(); site++){
					dna_p = p->compressedSequence[site];
					dna_c = c->compressedSequence[site];
					P(dna_p,dna_c) += this->siteWeights[site];
				}
			} else {
					cout << "rate category mismatch" << endl;
			}
		}
	}
	float rowSum;
	for (int row = 0; row < 4; row++){
		rowSum = 0;
		for (int col = 0; col < 4; col++){
			rowSum += P(row,col);
		}
		for (int col = 0; col < 4; col++){
			P(row,col)/=rowSum;
		}
	}	
	Q = P.log();	
	float l = Q(3,2);
	Q = Q/l;	
	this->rateMatrix = Q;
	
}

void rootedPhylogeny_tree::ComputeInitialEstimateForFreeParametersIncBaseFreq(){
//	int rateCategory = 1;
	Matrix4f P = ArrayXXf::Zero(4,4);
	Matrix4f Q;
//	int dna_p; int dna_c;
	this->ComputeMPEstimateOfAncestralSequences();
//	rootedPhylogeny_vertex * p;
//	rootedPhylogeny_vertex * c;
	float rowSum;
	for (int row = 0; row < 4; row++){
		rowSum = 0;
		for (int col = 0; col < 4; col++){
			rowSum += P(row,col);
		}
		for (int col = 0; col < 4; col++){
			P(row,col)/=rowSum;
		}
	}
	Q = P.log();
	MatrixXf pi;
	pi = this->ComputeStationaryDistribution(Q);
	float mu = 0;
	for (int i=0; i<4; i++){
		mu -= pi(i,0)*Q(i,i);
	}	
	Q = Q/mu;
	this->initialEstimateForFreeParameters << pi(0,0), pi(1,0), pi(2,0),
											   Q(0,1), Q(0,2), Q(0,3),
											   Q(1,0), Q(1,2), Q(1,3),
											   Q(2,0), Q(2,1);
}


void rootedPhylogeny_tree::InitializeConditionalLikelihoods(){
	unsigned char dna;
	rootedPhylogeny_vertex * v;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		this->root->conditionalLikelihood.push_back(this->GetLikelihoodArray(1,1,1,1));		
		for (int v_id = 0; v_id < this->numberOfObservedSequences; v_id++){
			v = (*this->vertexMap)[v_id];
			dna = v->compressedSequence[site];
			v->conditionalLikelihood.push_back(this->GetLikelihoodArray(0,0,0,0));			
			v->conditionalLikelihood[site][dna] = 1;
		}
		for (unsigned int v_id = this->numberOfObservedSequences; v_id < this->vertexMap->size()-1; v_id++){
			v = (*this->vertexMap)[v_id];
			v->conditionalLikelihood.push_back(this->GetLikelihoodArray(1,1,1,1));
		}
	}
}

void rootedPhylogeny_tree::ResetConditionalLikelihoodsForAncestors() {
	rootedPhylogeny_vertex * v;	
	for (unsigned int site = 0; site < this->siteWeights.size(); site++) {
		this->root->conditionalLikelihood[site] = this->GetLikelihoodArray(1,1,1,1);		
		for (unsigned int v_id = this->numberOfObservedSequences; v_id < this->vertexMap->size()-1; v_id++) {
			v = (*this->vertexMap)[v_id];
			v->conditionalLikelihood[site] = this->GetLikelihoodArray(1,1,1,1);
		}
	}
}

void rootedPhylogeny_tree::SetSiteWeights(vector <int> siteWeightsToSet) {
	this->siteWeights = siteWeightsToSet;
}

void rootedPhylogeny_tree::SetEdgesForPostOrderTreeTraversal() {
	this->edgesForPostOrderTreeTraversal->clear();
	vector < rootedPhylogeny_vertex *> verticesToVisit;
	for (pair <int, rootedPhylogeny_vertex*> idPtrPair: *this->vertexMap) {
		idPtrPair.second->timesVisited = 0;		
		if (idPtrPair.second->children_id.size()==0) {
			verticesToVisit.push_back(idPtrPair.second);
		}
	}	
	int numberOfVerticesToVisit = verticesToVisit.size();
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	while (numberOfVerticesToVisit > 0) {
		c = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		if (c->id != -1) {
			p = (*this->vertexMap)[c->parent_id];			
			this->edgesForPostOrderTreeTraversal->push_back(pair<int,int>(p->id,c->id));
			p->timesVisited += 1;
			if (p->timesVisited == (int) p->children_id.size()) {
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			}
		}
	}	
}

void rootedPhylogeny_tree::ComputeLogLikelihoodUsingStoredConditionals() {
	array <float,4> pi = this->rootProbability;
	this->logLikelihood = 0;
	double siteLikelihood;
	for (unsigned int site=0; site < this->siteWeights.size(); site++) {
		siteLikelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++) {
			siteLikelihood += pi[dna]*this->root->conditionalLikelihood[site][dna];
		}
		this->logLikelihood += log(siteLikelihood) * this->siteWeights[site];
	}
}

void rootedPhylogeny_tree::ComputeLogLikelihoodUsingDouble() {
	Matrix4f Q;	
	float scalingFactor;		
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f Q_scaled; Matrix4f P;
	float t;
//	using namespace boost::multiprecision;
//	this->SetMinLengthOfEdges();
//	double partialLikelihood;	
//	map <int, array<double,4>> conditionalLikelihoodMap;
//	array <double,4> conditionalLikelihood;
//	double maxConditionalLikelihood = 1;
//	double siteLikelihood;
//	double maxValueForConditionalLikelihood;
	double partialLikelihood;
	double siteLikelihood;
//	double logSiteLikelihood_mpf_float_1000;	
	map <int, array<double,4>> conditionalLikelihoodMap;
	array <double,4> conditionalLikelihood;
	double maxConditionalLikelihood = 1;
//	double logSiteLikelihood_double;
	bool showConditionalLikelihood = 0;	
	for (int i = 0; i < 4; i++){
		if (this->rootProbability[i] < 0){
			cout << "root prob is negative" << endl;
			cout << "Rate matrix is " << endl;
			cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
			break;
		}
//		cout << "root prob for nt " << i << " is " << this->rootProbability[i] << endl;
	}
	this->logLikelihood = 0;	
	for (unsigned int site = 0; site < this->siteWeights.size(); site++) {
		conditionalLikelihoodMap.clear();
		for (pair<int,int> edge_ids : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edge_ids.first];
			c = (*this->vertexMap)[edge_ids.second];
//			cout << "parent name is " << p->name << "\t" << "child name is " << c->name << endl;
			// Initialize conditional likelihood for child if child is a leaf
			if (c->id < this->numberOfObservedSequences) {				
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 0;
				}				
				conditionalLikelihood[c->compressedSequence[site]] = 1;
//				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for parent			
			if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()) {				
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 1;
//					if (p->name == "h_8"){
//						cout << "conditional likelihood for h_8 for dna " << dna << " is " << conditionalLikelihood[dna] << endl;
//					}
				}
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));
//				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));
			}
						
			partialLikelihood = 0;
			t = this->GetEdgeLength(p->id, c->id);
//			cout << "edge length for " << p->id << " - " << c->id << " is " << t << endl;
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];			
			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
//			cout << "scaling Factor is " << scalingFactor << endl;
			Q_scaled = Q*(t/scalingFactor);
			P = Q_scaled.exp();			
//			if (p->name == "h_8"){
//				cout << "Transition matrix is " << endl;
//				cout << P << endl;
//			}
//			maxValueForConditionalLikelihood = 0;
			for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				}
				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;				
//				if (conditionalLikelihoodMap[p->id][dna_p] > maxValueForConditionalLikelihood){
//					maxValueForConditionalLikelihood = conditionalLikelihoodMap[p->id][dna_p];				
//				}	
//				if (p->name == "h_1377"){
//					cout << "child name is " << c->name << endl;
//					cout << "partialLikelihood for dna " << static_cast <unsigned> (dna_p) << " is " << partialLikelihood << endl;
//				}
//				cout << "partial likelihood for dna " << static_cast<unsigned> (dna_p) << " is ";
//				cout << partialLikelihood << endl;
			}
//			if (this->smallestMaxValueForConditionalLikelihoodOfAncestors > maxValueForConditionalLikelihood){
//				this->smallestMaxValueForConditionalLikelihoodOfAncestors = maxValueForConditionalLikelihood;
//			}			
			// Erase conditional likelihood for child
//			if(conditionalLikelihoodMap.find(c->id) != conditionalLikelihoodMap.end()) {								
//				conditionalLikelihoodMap.erase(c->id);				
//			} else {
//				cout << "Child not present in conditional likelihood map" << endl;
//			}			
		}		
		for (rootedPhylogeny_vertex * v: *this->verticesForPreOrderTreeTraversal) {
			showConditionalLikelihood = 0;
			maxConditionalLikelihood = 0;
			for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
				if (maxConditionalLikelihood < conditionalLikelihoodMap[v->id][dna_c]) {
					maxConditionalLikelihood = conditionalLikelihoodMap[v->id][dna_c];
				}
			}
//			if (maxConditionalLikelihood < pow(10,-10)){
//				showConditionalLikelihood = 1;
//			}
			if (showConditionalLikelihood) {
				cout << (*this->vertexMap)[v->parent_id]->name << " - " << v->name << endl;
			}			
		}
		if (conditionalLikelihoodMap.find(-1) == conditionalLikelihoodMap.end()) {
			cout << "Conditional likelihood for root is not computed" << endl;
		}		
		siteLikelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++) {			
			siteLikelihood += this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
		}
		if (site == 0 or site == 10 or site == 20 or site == 30 or site == 40 or site == 50) {
//			cout << "Likelihood: sitelikelihood for site " << site << " is " << siteLikelihood << endl;
		}
		if (siteLikelihood == 0) {		
			cout << "Problem with computing siteLikelihood for site " << site << endl;			
			for (int i = 0; i < 4; i++){
//				cout << "Root probability for " << i << " is ";
//				cout << this->rootProbability[i] << endl;
//				cout << "Conditional likelihood for root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[-1][i] << endl;
//				cout << "Conditional likelihood for left child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[0]]->name << "of root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[0]][i] << endl;
//				cout << "Conditional likelihood for right child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[1]]->name << " of root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[1]][i] << endl;
//				cout << "Conditional likelihood for vertex 17 for " << i << " is ";
//				cout << conditionalLikelihoodMap[17][i] << endl;				
//				cout << "Rate matrix is" << endl;
//				cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
//				cout << "Number of rate categories is" << endl;
//				cout << this->numberOfRateCategories << endl;
			}
//			cout << "Edge length for vertex 17 is " << this->GetEdgeLength(17,(*this->vertexMap)[17]->parent_id) << endl;			
//			cout << "root probability for 0 is " << this->rootProbability[0] << endl;
			break;
		}
//		cout << "smallest maximum value of conditional likelihood is " << setprecision(100) << this->smallestMaxValueForConditionalLikelihoodOfAncestors << endl;
//		logSiteLikelihood_mpf_float_1000 = log(siteLikelihood);
//		logSiteLikelihood_double = logSiteLikelihood_mpf_float_1000.convert_to<double>();
		this->logLikelihood += log(siteLikelihood) * this->siteWeights[site];
		// Erase conditional likelihood for root
		conditionalLikelihoodMap.erase(-1);
	}
}

void rootedPhylogeny_tree::ResetLogScalingFactors() {
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : *this->vertexMap){
		idPtrPair.second->logScalingFactors = 0;
		idPtrPair.second->logScalingFactors_firstDer = 0;
		idPtrPair.second->logScalingFactors_secondDer = 0;
	}
}

void rootedPhylogeny_tree::ResetWeightedSiteLogLikelihoods() {
	this->unorderedWeightedSiteLogLikelihoods.clear();
}

void rootedPhylogeny_tree::AddWeightedSiteLogLikelihoods(){
	this->logLikelihood = 0;
	for (float weightedSiteLogLikelihood : this->unorderedWeightedSiteLogLikelihoods) {
		this->logLikelihood += weightedSiteLogLikelihood;
	}
}

void rootedPhylogeny_tree::ComputeWeightedSiteLogLikelihood(int site) {
	map <rootedPhylogeny_vertex*,array<float,4>> conditionalLikelihoodMap;
	array <float,4> conditionalLikelihood;
	float partialLikelihood;
	float siteLikelihood;
	float weightedSiteLikelihood;
	float largestConditionalLikelihood = 0;
	float currentProb;			
	vector <rootedPhylogeny_vertex *> verticesToVisit;	
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f P; Matrix4f Q; Matrix4f Q_scaled;
	float t;
	conditionalLikelihoodMap.clear();
	this->ResetLogScalingFactors();
	for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal) {
		p = (*this->vertexMap)[edgeIds.first];
		c = (*this->vertexMap)[edgeIds.second];
		t = this->GetEdgeLength(p->id, c->id);
		Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
		for (int row = 0; row < 4; row ++) {
			for (int col = 0; col < 4; col ++) {
				if (isnan(Q(row,col))) {
					cout << "Rate matrix contains nan" << endl;
					cout << "Rate category is " << c->rateCategory << endl;
				}
			}
		}			
		scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
//			cout << "scaling Factor is " << scalingFactor << endl;
		Q_scaled = Q*(t/scalingFactor);
		P = Q_scaled.exp();
		p->logScalingFactors += c->logScalingFactors;				
		// Initialize conditional likelihood for leaves
		if (c->children_id.size()==0) {
			for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
				conditionalLikelihood[dna_c] = 0;
			}
			conditionalLikelihood[c->compressedSequence[site]] = 1;
			conditionalLikelihoodMap.insert(pair <rootedPhylogeny_vertex *,array<float,4>>(c,conditionalLikelihood));
		}
		// Initialize conditional likelihood for ancestors			
		if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()){
			// Case 1: Ancestor is not an observed vertex
			if (p->id > this->numberOfObservedSequences -1) {
				for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
					conditionalLikelihood[dna_c] = 1;
				}
			} else {
			// Case 2: Ancestor is an observed vertex
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
					conditionalLikelihood[dna_c] = 0;
				}
				conditionalLikelihood[p->compressedSequence[site]] = 1;
			}								
			conditionalLikelihoodMap.insert(pair <rootedPhylogeny_vertex *,array<float,4>>(p,conditionalLikelihood));					
		}			
		if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()){
			for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
			conditionalLikelihood[dna_c] = 1;
			}				
			conditionalLikelihoodMap.insert(pair <rootedPhylogeny_vertex *,array<float,4>>(p,conditionalLikelihood));					
		}
		largestConditionalLikelihood = 0;		
		for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
			partialLikelihood = 0;
			for (unsigned char dna_c = 0; dna_c < 4; dna_c++) {
//					if (P(dna_p,dna_c)*conditionalLikelihoodMap[c][dna_c] == 0 and P(dna_p,dna_c) > 0 and conditionalLikelihoodMap[c][dna_c] > 0) {
//						cout << "Numerical underflow in computing partial likelihood" << endl;
//						cout << "P(y|x) is " << P(dna_p,dna_c) << endl;
//						cout << "L(y) is " << conditionalLikelihoodMap[c][dna_c] << endl;								
//						cout << "2^-256 is " << 1.0/pow(2,256) << endl;
//					}
				partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c][dna_c];
			}
			conditionalLikelihoodMap[p][dna_p] *= partialLikelihood;
			if (conditionalLikelihoodMap[p][dna_p] > largestConditionalLikelihood) {
				largestConditionalLikelihood = conditionalLikelihoodMap[p][dna_p];
			}
		}
		if (largestConditionalLikelihood != 0){
			for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
				conditionalLikelihoodMap[p][dna_p] /= largestConditionalLikelihood;
			}
			p->logScalingFactors += log(largestConditionalLikelihood);
		} else {
			cout << "Largest conditional likelihood value is zero" << endl;
			cout << "Transition matrix P is" << endl;
			cout << P << endl;
			cout << "Processing edge " << endl;
			cout << p->name << "\t" << c->name << endl;
			cout << "OutDegrees" << endl;
			cout << p->children_id.size() << "\t" << c->children_id.size() << endl;				
			exit(-1);						
		}
	}
	siteLikelihood = 0;						
	for (int dna = 0; dna < 4; dna++) {
		currentProb = this->rootProbability[dna] * conditionalLikelihoodMap[this->root][dna];
		siteLikelihood += currentProb;
	}
	weightedSiteLikelihood = (this->root->logScalingFactors + log(siteLikelihood)) * this->siteWeights[site];
	this->unorderedWeightedSiteLogLikelihoods.push_back(weightedSiteLikelihood); 
}

void rootedPhylogeny_tree::ComputeLogLikelihoodUsingMultiThreading(int numThreads){
	this->ResetWeightedSiteLogLikelihoods();
//	thread th(&rootedPhylogeny_tree::ComputeWeightedSiteLogLikelihood, this, 3);
//	th.join();
	cout << this->unorderedWeightedSiteLogLikelihoods.size();
}


void rootedPhylogeny_tree::ComputeLogLikelihood() {
	this->logLikelihood = 0;
	map <rootedPhylogeny_vertex*,array<float,4>> conditionalLikelihoodMap;
	array <float,4> conditionalLikelihood;
	float partialLikelihood;
	float siteLikelihood;	
	float largestConditionalLikelihood = 0;
	float currentProb;			
	vector <rootedPhylogeny_vertex *> verticesToVisit;	
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f P; Matrix4f Q; Matrix4f Q_scaled;
	float t;
	for (int site = 0; site < (int) this->siteWeights.size(); site++){
		conditionalLikelihoodMap.clear();
		this->ResetLogScalingFactors();
		for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edgeIds.first];
			c = (*this->vertexMap)[edgeIds.second];
			t = this->GetEdgeLength(p->id, c->id);
//			cout << "edge length for " << p->id << " - " << c->id << " is " << t << endl;
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
			for (int row = 0; row < 4; row ++) {
				for (int col = 0; col < 4; col ++) {
					if (isnan(Q(row,col))) {
						cout << "Rate matrix contains nan" << endl;
						cout << "Rate category is " << c->rateCategory << endl;
					}
				}
			}			
			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
//			cout << "scaling Factor is " << scalingFactor << endl;
			Q_scaled = Q*(t/scalingFactor);
			P = Q_scaled.exp();
			p->logScalingFactors += c->logScalingFactors;				
			// Initialize conditional likelihood for leaves
			if (c->children_id.size()==0) {
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
					conditionalLikelihood[dna_c] = 0;
				}
				conditionalLikelihood[c->compressedSequence[site]] = 1;
				conditionalLikelihoodMap.insert(pair <rootedPhylogeny_vertex *,array<float,4>>(c,conditionalLikelihood));
			}
			// Initialize conditional likelihood for ancestors			
			if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()){
				// Case 1: Ancestor is not an observed vertex
				if (p->id > this->numberOfObservedSequences -1) {
					for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
						conditionalLikelihood[dna_c] = 1;
					}
				} else {
				// Case 2: Ancestor is an observed vertex
					for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
						conditionalLikelihood[dna_c] = 0;
					}
					conditionalLikelihood[p->compressedSequence[site]] = 1;
				}								
				conditionalLikelihoodMap.insert(pair <rootedPhylogeny_vertex *,array<float,4>>(p,conditionalLikelihood));					
			}			
			if (conditionalLikelihoodMap.find(p) == conditionalLikelihoodMap.end()){
				for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
				conditionalLikelihood[dna_c] = 1;
				}				
				conditionalLikelihoodMap.insert(pair <rootedPhylogeny_vertex *,array<float,4>>(p,conditionalLikelihood));					
			}
			largestConditionalLikelihood = 0;		
			for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c++) {
//					if (P(dna_p,dna_c)*conditionalLikelihoodMap[c][dna_c] == 0 and P(dna_p,dna_c) > 0 and conditionalLikelihoodMap[c][dna_c] > 0) {
//						cout << "Numerical underflow in computing partial likelihood" << endl;
//						cout << "P(y|x) is " << P(dna_p,dna_c) << endl;
//						cout << "L(y) is " << conditionalLikelihoodMap[c][dna_c] << endl;								
//						cout << "2^-256 is " << 1.0/pow(2,256) << endl;
//					}
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c][dna_c];
				}
				conditionalLikelihoodMap[p][dna_p] *= partialLikelihood;
				if (conditionalLikelihoodMap[p][dna_p] > largestConditionalLikelihood) {
					largestConditionalLikelihood = conditionalLikelihoodMap[p][dna_p];
				}
			}
			if (largestConditionalLikelihood != 0){
				for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
					conditionalLikelihoodMap[p][dna_p] /= largestConditionalLikelihood;
				}
				p->logScalingFactors += log(largestConditionalLikelihood);
			} else {
				cout << "Largest conditional likelihood value is zero" << endl;
				cout << "Transition matrix P is" << endl;
				cout << P << endl;
				cout << "Processing edge " << endl;
				cout << p->name << "\t" << c->name << endl;
				cout << "OutDegrees" << endl;
				cout << p->children_id.size() << "\t" << c->children_id.size() << endl;				
				exit(-1);						
			}
		}
		siteLikelihood = 0;						
		for (int dna = 0; dna < 4; dna++) {
			currentProb = this->rootProbability[dna] * conditionalLikelihoodMap[this->root][dna];
			siteLikelihood += currentProb;
		}
//		if (site == 0) {
//			cout << "Root probability is" << endl;
//			for (int i = 0; i < 4; i++) {
//				cout << this->rootProbability[i] << "\t";
//			}
//			cout << endl;
//		}
//		if (site == 0) {
//			cout << "LogLikelihood for site 0 is " ;
//			cout << this->root->logScalingFactors + log(siteLikelihood) << endl;
//		}
		this->logLikelihood += (this->root->logScalingFactors + log(siteLikelihood)) * this->siteWeights[site];				
	}
}

void rootedPhylogeny_tree::ComputeLogLikelihoodUsingLongDouble() {
	Matrix4f Q;	
	float scalingFactor;		
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f Q_scaled; Matrix4f P;
	float t;
//	using namespace boost::multiprecision;
//	this->SetMinLengthOfEdges();
//	double partialLikelihood;	
//	map <int, array<double,4>> conditionalLikelihoodMap;
//	array <double,4> conditionalLikelihood;
//	double maxConditionalLikelihood = 1;
//	double siteLikelihood;
//	long double maxValueForConditionalLikelihood;
	long double partialLikelihood;
	long double siteLikelihood;
//	long double logSiteLikelihood_mpf_float_1000;	
	map <int, array<long double,4>> conditionalLikelihoodMap;
	array <long double,4> conditionalLikelihood;
	long double maxConditionalLikelihood = 1;
//	double logSiteLikelihood_double;
	bool showConditionalLikelihood = 0;	
	for (int i = 0; i < 4; i++) {
		if (this->rootProbability[i] < 0){
			cout << "Root prob is negative" << endl;
			cout << "Rate matrix is " << endl;
			cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
			break;
		}
//		cout << "root prob for nt " << i << " is " << this->rootProbability[i] << endl;
	}
	this->logLikelihood = 0;	
	for (unsigned int site = 0; site < this->siteWeights.size(); site++) {
		conditionalLikelihoodMap.clear();
		for (pair<int,int> edge_ids : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edge_ids.first];
			c = (*this->vertexMap)[edge_ids.second];
//			cout << "parent name is " << p->name << "\t" << "child name is " << c->name << endl;
			// Initialize conditional likelihood for child if child is a leaf
			if (c->id < this->numberOfObservedSequences) {				
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 0;
				}				
				conditionalLikelihood[c->compressedSequence[site]] = 1;
//				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for parent			
			if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()) {				
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 1;
//					if (p->name == "h_8"){
//						cout << "conditional likelihood for h_8 for dna " << dna << " is " << conditionalLikelihood[dna] << endl;
//					}
				}
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));
//				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));
			}
						
			partialLikelihood = 0;
			t = this->GetEdgeLength(p->id, c->id);
//			cout << "edge length for " << p->id << " - " << c->id << " is " << t << endl;
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];			
			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
//			cout << "scaling Factor is " << scalingFactor << endl;
			Q_scaled = Q*(t/scalingFactor);
			P = Q_scaled.exp();			
//			if (p->name == "h_8"){
//				cout << "Transition matrix is " << endl;
//				cout << P << endl;
//			}
//			maxValueForConditionalLikelihood = 0;
			for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				}
				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;				
//				if (conditionalLikelihoodMap[p->id][dna_p] > maxValueForConditionalLikelihood){
//					maxValueForConditionalLikelihood = conditionalLikelihoodMap[p->id][dna_p];				
//				}	
//				if (p->name == "h_1377"){
//					cout << "child name is " << c->name << endl;
//					cout << "partialLikelihood for dna " << static_cast <unsigned> (dna_p) << " is " << partialLikelihood << endl;
//				}
//				cout << "partial likelihood for dna " << static_cast<unsigned> (dna_p) << " is ";
//				cout << partialLikelihood << endl;
			}
//			if (this->smallestMaxValueForConditionalLikelihoodOfAncestors > maxValueForConditionalLikelihood){
//				this->smallestMaxValueForConditionalLikelihoodOfAncestors = maxValueForConditionalLikelihood;
//			}			
			// Erase conditional likelihood for child
//			if(conditionalLikelihoodMap.find(c->id) != conditionalLikelihoodMap.end()) {								
//				conditionalLikelihoodMap.erase(c->id);				
//			} else {
//				cout << "Child not present in conditional likelihood map" << endl;
//			}			
		}		
		for (rootedPhylogeny_vertex * v: *this->verticesForPreOrderTreeTraversal) {
			showConditionalLikelihood = 0;
			maxConditionalLikelihood = 0;
			for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
				if (maxConditionalLikelihood < conditionalLikelihoodMap[v->id][dna_c]) {
					maxConditionalLikelihood = conditionalLikelihoodMap[v->id][dna_c];
				}
			}
//			if (maxConditionalLikelihood < pow(10,-10)){
//				showConditionalLikelihood = 1;
//			}
			if (showConditionalLikelihood) {
				cout << (*this->vertexMap)[v->parent_id]->name << " - " << v->name << endl;
			}			
		}
		if (conditionalLikelihoodMap.find(-1) == conditionalLikelihoodMap.end()) {
			cout << "Conditional likelihood for root is not computed" << endl;
		}		
		siteLikelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++) {			
			siteLikelihood += this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
		}
		if (site == 0 or site == 10 or site == 20 or site == 30 or site == 40 or site == 50) {
//			cout << "Likelihood: sitelikelihood for site " << site << " is " << siteLikelihood << endl;
//			for (int i = 0; i < 4; i++){
//				cout << "Root probability for " << i << " is ";
//				cout << this->rootProbability[i] << endl;
//				cout << "Conditional likelihood for root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[-1][i] << endl;
//				cout << "Conditional likelihood for left child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[0]]->name << "of root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[0]][i] << endl;
//				cout << "Edge length for left child is " << this->GetEdgeLength(-1,(*this->vertexMap)[-1]->children_id[0]) << endl;				
//				cout << "Conditional likelihood for right child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[1]]->name << " of root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[1]][i] << endl;
//				cout << "Edge length for right child is " << this->GetEdgeLength(-1,(*this->vertexMap)[-1]->children_id[1]) << endl;
//				cout << "Conditional likelihood for vertex 17 for " << i << " is ";
//				cout << conditionalLikelihoodMap[17][i] << endl;				
//				cout << "Rate matrix is" << endl;
//				cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
//				cout << "Number of rate categories is" << endl;
//				cout << this->numberOfRateCategories << endl;
//			}
		}
		if (siteLikelihood == 0) {		
			cout << "Problem with computing siteLikelihood for site " << site << endl;			
			for (int i = 0; i < 4; i++){
//				cout << "Root probability for " << i << " is ";
//				cout << this->rootProbability[i] << endl;
//				cout << "Conditional likelihood for root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[-1][i] << endl;
//				cout << "Conditional likelihood for left child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[0]]->name << "of root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[0]][i] << endl;
//				cout << "Conditional likelihood for right child " << (*this->vertexMap)[(*this->vertexMap)[-1]->children_id[1]]->name << " of root for " << i << " is ";				
//				cout << conditionalLikelihoodMap[(*this->vertexMap)[-1]->children_id[1]][i] << endl;
//				cout << "Conditional likelihood for vertex 17 for " << i << " is ";
//				cout << conditionalLikelihoodMap[17][i] << endl;				
//				cout << "Rate matrix is" << endl;
//				cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
//				cout << "Number of rate categories is" << endl;
//				cout << this->numberOfRateCategories << endl;
			}
//			cout << "Edge length for vertex 17 is " << this->GetEdgeLength(17,(*this->vertexMap)[17]->parent_id) << endl;			
//			cout << "root probability for 0 is " << this->rootProbability[0] << endl;
			break;
		}
//		cout << "smallest maximum value of conditional likelihood is " << setprecision(100) << this->smallestMaxValueForConditionalLikelihoodOfAncestors << endl;
//		logSiteLikelihood_mpf_float_1000 = log(siteLikelihood);
//		logSiteLikelihood_double = logSiteLikelihood_mpf_float_1000.convert_to<double>();
		this->logLikelihood += log(siteLikelihood) * this->siteWeights[site];
		// Erase conditional likelihood for root
		conditionalLikelihoodMap.erase(-1);
	}
}

MatrixXd rootedPhylogeny_tree::GetJacobianForRateCategory(int rateCategoryForOptimization) {		
//	using namespace boost::multiprecision;	
	map <int, array <long double,4>> conditionalLikelihoodMap;	
	map <int, array <array <long double,4>,11> > derivativeOfConditionalLikelihoodMap;
	array <long double,4> conditionalLikelihood;
	array <float,4> pi_root;
	array <array <long double,4>,11> derivativeOfConditionalLikelihood;
	long double siteLikelihood;
	long double partialLikelihood;
	long double partialLikelihood_l;
	long double partialLikelihood_r;
	long double derivativeOfPartialLikelihood_l;
	long double derivativeOfPartialLikelihood_r;
	long double derivativeOfSiteLikelihood;
//	long double Jacobian_per_site_per_par;	
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * v;
	Matrix4f Q;	
	Matrix4f Q_scaled; Matrix4f P;
	Matrix4f Q_l; Matrix4f Q_r ; 
	Matrix4f Q_l_scaled; Matrix4f Q_r_scaled;
	Matrix4f P_l; Matrix4f P_r;
	Matrix4f P_l_der; Matrix4f P_r_der;
	float t; float t_l; float t_r;			
	int numberOfVerticesToVisit;
	rootedPhylogeny_vertex * c_l; rootedPhylogeny_vertex * c_r;
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	// Initialize Jacobian matrix
	MatrixXd Jacobian = MatrixXd::Zero(11,1);
	MatrixXd B_rate_cat_root;
	B_rate_cat_root = (*this->parametersPerRateCategory)[this->root->rateCategory];
	pi_root[0] = B_rate_cat_root(0,0);
	pi_root[1] = B_rate_cat_root(1,0);
	pi_root[2] = B_rate_cat_root(2,0);
	pi_root[3] = (1 - pi_root[0] - pi_root[1] - pi_root[2]);
	for (unsigned int site = 0; site < this->siteWeights.size(); site++) {
		// Compute conditional likelihoods
		for (pair<int,int> edge_ids : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edge_ids.first];
			c = (*this->vertexMap)[edge_ids.second];			
			// Set conditional likelihood for child if child is a leaf
			if (c->id < this->numberOfObservedSequences) {
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 0;
				}				
				conditionalLikelihood[c->compressedSequence[site]] = 1;
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for parent
			if(conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){				
				for (int dna = 0; dna < 4; dna++){
					conditionalLikelihood[dna] = 1;
				}
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));
			}
			partialLikelihood = 0;
			t = this->GetEdgeLength(p->id, c->id);
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				}
				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
			}
		}
		// Set derivative of conditional likelihood for leaves		
		for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap) {
			v = idPtrPair.second;
			if (v->children_id.size()==0) {
				for (int par = 0; par < 11; par ++) {
					for (int dna = 0; dna < 4; dna++) {		
						derivativeOfConditionalLikelihood[par][dna] = 0;
					}
				}
				derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<long double,4>,11>>(v->id,derivativeOfConditionalLikelihood));
			}			
		}
		verticesToVisit.clear();
		for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap) {
			v = idPtrPair.second;
			v->timesVisited = 0;
			if (v->children_id.size()==0) {
				verticesToVisit.push_back(v);
			}
		}
		numberOfVerticesToVisit = verticesToVisit.size();
		while (numberOfVerticesToVisit > 0) {
			c = verticesToVisit[numberOfVerticesToVisit-1];
			numberOfVerticesToVisit -= 1;
			verticesToVisit.pop_back();
			if (c->id != -1) {
				p = (*this->vertexMap)[c->parent_id];
				p->timesVisited += 1;
				if (p->timesVisited == 2) {
					c_l = (*this->vertexMap)[p->children_id[0]];
					c_r = (*this->vertexMap)[p->children_id[1]];
					t_l = this->GetEdgeLength(p->id,c_l->id);
					t_r = this->GetEdgeLength(p->id,c_r->id);
					Q_l = (*this->rateMatrixPerRateCategory)[c_l->rateCategory];
					Q_r = (*this->rateMatrixPerRateCategory)[c_r->rateCategory];
					Q_l_scaled = Q_l * t_l;
					Q_r_scaled = Q_r * t_r;
					P_l = Q_l_scaled.exp();
					P_r = Q_r_scaled.exp();					
					for (int par = 0; par < 11; par++) {
						for (int dna_p = 0; dna_p < 4; dna_p++) {
							derivativeOfConditionalLikelihood[par][dna_p] = 0;
						}
					}
					for (int dna_p = 0; dna_p < 4; dna_p ++) {
						partialLikelihood_l = 0;
						partialLikelihood_r = 0;
						for (int dna_c = 0; dna_c < 4; dna_c ++) {
							partialLikelihood_l += P_l(dna_p,dna_c)*conditionalLikelihoodMap[c_l->id][dna_c];
							partialLikelihood_r += P_r(dna_p,dna_c)*conditionalLikelihoodMap[c_r->id][dna_c];
						}
						for (int par = 0; par < 11; par++) {
							if (c_l->rateCategory == rateCategoryForOptimization) {
								P_l_der = this->ComputeFirstDerivativeOfMatrixExponential(t_l, c_l->rateCategory, par);
							} else {
								P_l_der = MatrixXf::Zero(4,4);
							}
							if (c_r->rateCategory == rateCategoryForOptimization) {
								P_r_der = this->ComputeFirstDerivativeOfMatrixExponential(t_r, c_r->rateCategory, par);
							} else {
								P_r_der = MatrixXf::Zero(4,4);
							}						
							derivativeOfPartialLikelihood_l = 0;
							derivativeOfPartialLikelihood_r = 0;
							for (int dna_c = 0; dna_c < 4; dna_c ++) {
								derivativeOfPartialLikelihood_l += P_l_der(dna_p,dna_c)*conditionalLikelihoodMap[c_l->id][dna_c];
								derivativeOfPartialLikelihood_l += P_l(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_l->id][par][dna_c];
								derivativeOfPartialLikelihood_r += P_r_der(dna_p,dna_c)*conditionalLikelihoodMap[c_r->id][dna_c];
								derivativeOfPartialLikelihood_r += P_r(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_r->id][par][dna_c];
							}
							derivativeOfConditionalLikelihood[par][dna_p] = derivativeOfPartialLikelihood_l*partialLikelihood_r;
							derivativeOfConditionalLikelihood[par][dna_p] += derivativeOfPartialLikelihood_r*partialLikelihood_l;
						}
					}
					derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<long double,4>,11>>(p->id,derivativeOfConditionalLikelihood));					
					verticesToVisit.push_back(p);
					numberOfVerticesToVisit += 1;
				}
			}
		}		
		if (conditionalLikelihoodMap.find(-1) == conditionalLikelihoodMap.end()) {
				cout << "conditional likelihood for root is not computed" << endl;
		}		
		siteLikelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++) {
			siteLikelihood += pi_root[dna]*conditionalLikelihoodMap[-1][dna];
		}		
		if (derivativeOfConditionalLikelihoodMap.find(-1) == derivativeOfConditionalLikelihoodMap.end()) {
				cout << "derivative of conditional likelihood for root is not computed" << endl;
		}		
		for	(int par = 0; par < 11; par++) {
			derivativeOfSiteLikelihood = 0;			
			if (this->root->rateCategory == rateCategoryForOptimization){
				if (par == 0){
					// par is pi_1
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][0];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				} else if (par == 1){
					// par is pi_2
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][1];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				} else if (par == 2){
					// par is pi_3
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][2];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				}			
			}			
			for (int dna = 0; dna < 4; dna++) {
				derivativeOfSiteLikelihood += this->rootProbability[dna]*derivativeOfConditionalLikelihoodMap[-1][par][dna];
			}			
			Jacobian(par,0) += (derivativeOfSiteLikelihood/siteLikelihood)*this->siteWeights[site];
		}
		conditionalLikelihoodMap.clear();
		derivativeOfConditionalLikelihoodMap.clear();
	}
	return (Jacobian);
}

void rootedPhylogeny_tree::ComputeJacobianForRateMatrixForRateCategoryDep(int rateCat) {
	Matrix4f Q;
	Q = this->rateMatrix;
//	for (int par = 0; par < 11; par++){
//		this->jacobian_dep[par] = 0;
//	}	
	map <int, array <long double,4>> conditionalLikelihoodMap;
	map <int, array <array <long double,4>,11>> derivativeOfConditionalLikelihoodMap;
	array <long double,4> conditionalLikelihood;
	array <array <long double,4>,11> derivativeOfConditionalLikelihood;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
//	rootedPhylogeny_vertex * v;	
//	array <float,4> pi = this->rootProbability;	
	vector <rootedPhylogeny_vertex *> verticesToVisit;
	for (pair <int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		idPtrPair.second->timesVisited = 0;
		if (idPtrPair.second->children_id.size() == 0){
			verticesToVisit.push_back(idPtrPair.second);
		}
	}
	vector <rootedPhylogeny_vertex *> postOrderVerticesToVisit;	
	int numberOfVerticesToVisit = verticesToVisit.size();	
	while (numberOfVerticesToVisit > 0){
		c = verticesToVisit[numberOfVerticesToVisit-1];	
		numberOfVerticesToVisit -= 1;
		verticesToVisit.pop_back();
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			p->timesVisited += 1;			
			if (p->timesVisited == 2){
				postOrderVerticesToVisit.push_back(p);
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			}
		}
	}	
	int c_l_id; int c_r_id;
//	rootedPhylogeny_vertex* c_l; rootedPhylogeny_vertex* c_r;	
	float t_l; float t_r;
	Matrix4f Q_l_scaled; Matrix4f Q_r_scaled;
	Matrix4f P_l; Matrix4f P_r;
	Matrix4f P_l_der; Matrix4f P_r_der;
	P_l_der = ArrayXXf::Zero(4,4);
	P_r_der = ArrayXXf::Zero(4,4);
	long double partialLikelihood_l; long double partialLikelihood_r;
	long double partialDerivativeOfLikelihood_l; long double partialDerivativeOfLikelihood_r; 
	long double siteLikelihood; long double derivativeOfSiteLikelihood;
//	long double derivativeOfRootProbability;
	for (int site = 0; site < this->numberOfObservedSequences; site++){
		conditionalLikelihoodMap.clear();
		derivativeOfConditionalLikelihoodMap.clear();
		for (rootedPhylogeny_vertex * p: postOrderVerticesToVisit){
			c_l_id = p->children_id[0];
			c_r_id = p->children_id[1];
//				c_l = (*this->vertexMap)[c_l_id];
//				c_r = (*this->vertexMap)[c_r_id];
			for (int c_id : p->children_id){
				if (c_id < this->numberOfObservedSequences){
					for (int dna = 0; dna < 4; dna++){
						conditionalLikelihood[dna] = 0;
					}
					conditionalLikelihood[(*this->vertexMap)[c_id]->compressedSequence[site]] = 1;
					conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c_id,conditionalLikelihood));
					for (int par = 0; par < 11; par ++){
						for (int dna = 0; dna < 4; dna++){					
							derivativeOfConditionalLikelihood[par][dna] = 0;
						}
					}
					derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<long double,4>,11>>(c_id,derivativeOfConditionalLikelihood));
				} else{
					// check if child is present in conditional and derivative of conditional maps
					
				}
			}			
			// Compute conditional likelihood for parent
			t_l = this->GetEdgeLength(p->id,c_l_id);
			t_r = this->GetEdgeLength(p->id,c_r_id);
			Q_l_scaled = Q*t_l;
			Q_r_scaled = Q*t_r;
			P_l = Q_l_scaled.exp();
			P_r = Q_r_scaled.exp();
			partialLikelihood_l = 0; partialLikelihood_r = 0;
			// Initialize conditional likelihood for parent
			for (int dna_p = 0; dna_p < 4; dna_p++){
				conditionalLikelihood[dna_p] = 0;
			}
			for (int dna_p = 0; dna_p < 4; dna_p++){
				for (int dna_c = 0; dna_c < 4; dna_c++){
					partialLikelihood_l += P_l(dna_p,dna_c)*conditionalLikelihoodMap[c_l_id][dna_c];
				}
				for (int dna_c = 0; dna_c < 4; dna_c++){
					partialLikelihood_r += P_r(dna_p,dna_c)*conditionalLikelihoodMap[c_r_id][dna_c];
				}
				conditionalLikelihood[dna_p] = partialLikelihood_l*partialLikelihood_r;
			}			
			conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));			
			// Initialize derivative of conditional likelihood for parent			
			for (int par = 0; par < 11; par ++){
				for (int dna = 0; dna < 4; dna++){		
					derivativeOfConditionalLikelihood[par][dna] = 0;
				}
			}
			derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<long double,4>,11>>(p->id,derivativeOfConditionalLikelihood));
			// Iterate over parameters
			for (int par = 0; par < 11; par++){
				P_l_der = this->ComputeDerivativeOfMatrixExponentialDep(t_l, par);
				P_r_der = this->ComputeDerivativeOfMatrixExponentialDep(t_r, par);
				for (int dna_p = 0; dna_p < 4; dna_p++){
					partialLikelihood_l = 0;
					for (int dna_c = 0; dna_c < 4; dna_c++){
						partialLikelihood_l += P_l(dna_p,dna_c)*conditionalLikelihoodMap[c_l_id][dna_c];
					}
					partialLikelihood_r = 0;
					for (int dna_c = 0; dna_c < 4; dna_c++){
						partialLikelihood_r += P_r(dna_p,dna_c)*conditionalLikelihoodMap[c_r_id][dna_c];
					}
					partialDerivativeOfLikelihood_l = 0; 
					for (int dna_c = 0; dna_c < 4; dna_c++){
						partialDerivativeOfLikelihood_l += P_l_der(dna_p,dna_c)*conditionalLikelihoodMap[c_l_id][dna_c];
						partialDerivativeOfLikelihood_l += P_l(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_l_id][par][dna_c];
					}
					partialDerivativeOfLikelihood_r = 0;
					for (int dna_c = 0; dna_c < 4; dna_c++){
						partialDerivativeOfLikelihood_r += P_r_der(dna_p,dna_c)*conditionalLikelihoodMap[c_r_id][dna_c];
						partialDerivativeOfLikelihood_r += P_r(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_r_id][par][dna_c];
					}
					derivativeOfConditionalLikelihoodMap[p->id][par][dna_p] = partialDerivativeOfLikelihood_l*partialLikelihood_r;
					derivativeOfConditionalLikelihoodMap[p->id][par][dna_p] += partialLikelihood_l*partialDerivativeOfLikelihood_r;					
				}					
			}
			// Remove conditional likelihood for children
			conditionalLikelihoodMap.erase(c_l_id);
			conditionalLikelihoodMap.erase(c_r_id);				
			derivativeOfConditionalLikelihoodMap.erase(c_l_id);
			derivativeOfConditionalLikelihoodMap.erase(c_r_id);
		}
		siteLikelihood = 0;
		for (int dna = 0; dna < 4; dna++){
			siteLikelihood += this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
		}
		if (site == 0 or site == 10 or site == 20 or site == 30 or site == 40 or site == 50){
			cout << "Jacobian: sitelikelihood for site " << site << " is " << siteLikelihood << endl;
		}
//			cout << "site likelihood for site " << site << "is " << siteLikelihood << endl;
		for	(int par = 0; par < 11; par++){
			derivativeOfSiteLikelihood = 0;			
			if (par == 0){
				// par is pi_1
				derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][0];
				derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
			} else if (par == 1){
				// par is pi_2
				derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][1];
				derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
			} else if (par == 2){
				// par is pi_3
				derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][2];
				derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
			}			
			for (int dna = 0; dna < 4; dna++){
				derivativeOfSiteLikelihood += this->rootProbability[dna]*derivativeOfConditionalLikelihoodMap[-1][par][dna];
			}
			this->JacobianForRateMatrix(par,0) += (derivativeOfSiteLikelihood/siteLikelihood)*this->siteWeights[site];
		}
		
	}
	
}
	
//	double logLikelihood_current = 0;
//	double logLikelihood_updated;
//	int iter = 0;
//	int maxNumberOfIterations = 20;
//	bool continueOptimization = 1;	
//	this->ComputeMPEstimateOfAncestralSequences();
//	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
//		if (idPtrPair.second->compressedSequence.size() != this->siteWeights.size()){
//			cout << "sequence length mismatch" << endl;
//		}
//	}
//	this->ComputeInitialEstimateForRateMatrix();	
//	this->ScaleEdgeLengths();
//	this->SetMinLengthOfEdges();
//	// Initial log likelihood
//	this->ComputeLogLikelihood();
//	//  cout << "initial loglikelihood is " << this->logLikelihood << endl;
//	while (continueOptimization){
//		iter += 1;
//		logLikelihood_current = this->logLikelihood;
//		this->ComputeMLEOfRootProbability();		
//		for (int i = 0; i < 5; i++){
//			// Optimize rate matrices
//			this->ComputeMLEOfRateMatrices();
//			// Optimize edge lengths
//			this->ComputeMLEOfEdgeLengths();
//			// Compute expected states	
//		}		
//		this->ComputeMAPEstimateOfAncestralSequences();
//		logLikelihood_updated = this->logLikelihood;
//		cout << "updated log likelihood is " << logLikelihood_updated << endl;
//	//	cout << "updated loglikelihood is " << this->logLikelihood << endl;
//		for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
//			if (idPtrPair.second->compressedSequence.size() != this->siteWeights.size()){
//				cout << "sequence length mismatch" << endl;
//			}
//		}
//		if (iter < 5 or abs(logLikelihood_updated-logLikelihood_current)>0.001){
//			continueOptimization = 1;
//		} else {
//			continueOptimization = 0;
//		}		
//	}	
	
//	this->ComputeLogLikelihoodForFullyLabeledTree();
//	cout << "current logLikelihood is " << endl;
//	cout << this->logLikelihood << endl;
	// ML estimate of rate matrices
//	this->ComputeMLEOfRateMatrices();
//	this->ComputeLogLikelihoodForFullyLabeledTree();
//	cout << "logLikelihood after optimizing rate matrices is " << endl;
//	cout << this->logLikelihood << endl;
	// ML estimate of edge lengths
//	this->ComputeMLEOfEdgeLengths();
//	this->ComputeLogLikelihoodForFullyLabeledTree();
//	cout << "logLikelihood after optimizing edge lengths is " << endl;
	// check convergence of log-likelihood
//	cout << this->logLikelihood << endl;
//	if (currentLogLikelihood < this->logLikelihood){
//		currentLogLikelihood = this->logLikelihood;
//	}

void rootedPhylogeny_tree::AddRateCategoryForVertex(int v_id, int rateCategory){
	(*this->vertexMap)[v_id]->rateCategory = rateCategory;	
}

void rootedPhylogeny_tree::SetChangePointForRateCategory(int v_id, int rateCategory){
	(*this->changePointForRateCategory)[rateCategory] = v_id;
}

void rootedPhylogeny_tree::SetRateMatrixForRateCategory(Matrix4f Q ,int rateCategory){
//	this->parametersPerRateCategory->insert(pair<int,Matrix4f>(rateCategory,Q));	
//	float scalingFactor = this->ComputeScalingFactor(Q);
//	this->scalingFactorForRateCategory->insert(pair<int,float>(rateCategory,scalingFactor));
}

void rootedPhylogeny_tree::AddStationaryDistributionForCategory(MatrixXf  stationaryDistribution ,int rateCategory){
	this->stationaryDistributionForCategory->insert(pair<int,MatrixXf >(rateCategory,stationaryDistribution));
}

void rootedPhylogeny_tree::WriteAncestralSequences(string sequenceFileName){
	ofstream ancestralSequencesFile;
	ancestralSequencesFile.open(sequenceFileName+".modelSelection_ancestralSequences");
	for (pair<int,rootedPhylogeny_vertex*> vIdAndPtr:(*this->vertexMap)){
		if (boost::algorithm::starts_with(vIdAndPtr.second->name , "h_")){			
			ancestralSequencesFile << ">" << vIdAndPtr.second->name << endl;
			ancestralSequencesFile << EncodeAsDNA(vIdAndPtr.second->sequence) << endl;
		} else {			
		}	
	}	
	ancestralSequencesFile.close();
}

bool rootedPhylogeny_tree::IsEdgeContractionFeasbile(int vertex_id_1, int vertex_id_2){	
	bool edgeContractionIsFeasible;
	if (vertex_id_1 > this->numberOfObservedSequences or vertex_id_2 > this->numberOfObservedSequences){
		edgeContractionIsFeasible = 1;
	} else {
		edgeContractionIsFeasible = 0;
	}
	return (edgeContractionIsFeasible);
}

void rootedPhylogeny_tree::ContractEdge(int idOfVertexToKeep, int idOfVertexToRemove){
	rootedPhylogeny_vertex * k = (*this->vertexMap)[idOfVertexToKeep];
	rootedPhylogeny_vertex * r = (*this->vertexMap)[idOfVertexToRemove];
	float edgeLengthToAdd;
	int ind;	
	if (r->parent_id == k->id){
		// Case 1: k is parent of r
		// Remove r from list of children of k
		ind = find(k->children_id.begin(),k->children_id.end(),r->id) - k->children_id.begin();
		k->children_id.erase(k->children_id.begin()+ind);
		for (int childOfR_id: r->children_id){
			// Set parent of childOfR to k
			(*this->vertexMap)[childOfR_id]->parent_id = k->id;			
			// Add edge from k to childOfR
			k->children_id.push_back(childOfR_id);
			edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->id,childOfR_id)];
			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->id,childOfR_id),edgeLengthToAdd));
			// Remove edge from r to child of r
			this->edgeLengthsMap->erase(pair<int,int>(r->id,childOfR_id));
		}		
	} else {
		// Case 2: k is child of r		
		for (int childOfR_id: r->children_id){
			if (childOfR_id != k->id){
				// Set parent of childOfR to k
				(*this->vertexMap)[childOfR_id]->parent_id = k->id;
				// Add edge from k to childOfR
				k->children_id.push_back(childOfR_id);
				edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->id,childOfR_id)];
				this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->id,childOfR_id),edgeLengthToAdd));
				// Remove edge from r to child of r
				this->edgeLengthsMap->erase(pair<int,int>(r->id,childOfR_id));
			}
		}
		// Set parent of k to parent of r
		k->parent_id = r->parent_id;
		// If k is not the root then:
		if (k->parent_id != -1){
			// Add edge from parent of k to k
			edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->parent_id,r->id)];
			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->parent_id,k->id),edgeLengthToAdd));
			// Remove edge from parent of r to r
			this->edgeLengthsMap->erase(pair<int,int>(r->parent_id,r->id));	
		}		
	}
	// Remove all children of r
	r->children_id.clear();
	// Set parent of r to r
	r->parent_id = r->id;
}

void rootedPhylogeny_tree::AddNumberOfObservedSequences(int numberOfObservedSequencesToAdd){
	this->numberOfObservedSequences = numberOfObservedSequencesToAdd;
}

void rootedPhylogeny_tree::AddVertex(int idToAdd, string nameToAdd, vector<unsigned char> sequenceToAdd){
	rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(idToAdd, nameToAdd, sequenceToAdd);
	(*this->vertexMap)[idToAdd] = v;
	
}

void rootedPhylogeny_tree::AddVertex(int idToAdd, string nameToAdd, vector<unsigned char> sequenceToAdd, vector<unsigned char> compressedSequenceToAdd){
	rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(idToAdd, nameToAdd, sequenceToAdd);
	v->compressedSequence = compressedSequenceToAdd;
	(*this->vertexMap)[idToAdd] = v;
}


void rootedPhylogeny_tree::RemoveEdges(){	
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		idPtrPair.second->parent_id = -1;
		idPtrPair.second->children_id.clear();
	}
}

void rootedPhylogeny_tree::AddDirectedEdges(vector <pair<int,int>> * directedEdgeList_ptr){
	// Remove edge lengths for root
	for (int child_id : this->root->children_id){
		if(this->edgeLengthsMap->find(pair<int,int>(this->root->id,child_id)) != this->edgeLengthsMap->end()){
			this->edgeLengthsMap->erase(pair<int,int>(this->root->id,child_id));
		}
	}
	this->RemoveEdges();	
	for (pair<int,int> pIdAndcIdPair : *directedEdgeList_ptr){		
		if (this->vertexMap->find(pIdAndcIdPair.second) == this->vertexMap->end()){
			cout << "vertex " << pIdAndcIdPair.second << " not in vertex map" << endl;
		}
		this->AddEdge(pIdAndcIdPair.first,pIdAndcIdPair.second);
		if (pIdAndcIdPair.first == -1){			
		}
	}
	int u_id; int v_id;
	// Set edge lengths for root by midpoint rooting	
	u_id = this->root->children_id[0];	
	v_id = this->root->children_id[1];	
	float edgeLength = this->GetEdgeLength(u_id,v_id);	
	this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(this->root->id,u_id),edgeLength/float(2)));
	this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(this->root->id,v_id),edgeLength/float(2)));
	this->SetEdgesForPostOrderTreeTraversal();
	
}

void rootedPhylogeny_tree::AddEdge(int p_id, int c_id){		
	(*this->vertexMap)[p_id]->AddChild(c_id);
	(*this->vertexMap)[c_id]->AddParent(p_id);	
}

void rootedPhylogeny_tree::ComputeAndSetEdgeLength(int u_id, int v_id){
	float edgeLength;	
	rootedPhylogeny_vertex * u = (*this->vertexMap)[u_id];
	rootedPhylogeny_vertex * v = (*this->vertexMap)[v_id];		
	edgeLength = ComputeNormalizedHammingDistance(&u->sequence,&v->sequence);	
	if (u_id < v_id){
		this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(u_id,v_id),edgeLength));
	} else {
		this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(v_id,u_id),edgeLength));
	}
}

float rootedPhylogeny_tree::GetEdgeLength(int u_id, int v_id){
	float edgeLength;
	if (u_id < v_id){
		edgeLength = (*this->edgeLengthsMap)[pair<int,int>(u_id,v_id)];
	} else {
		edgeLength = (*this->edgeLengthsMap)[pair<int,int>(v_id,u_id)];
	}
	return edgeLength;
}
void rootedPhylogeny_tree::ComputeEdgeLengths(){
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	float edgeLength;
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair: *this->vertexMap){
		c = (*this->vertexMap)[idPtrPair.first];
		if (c->parent_id != -1){
			p = (*this->vertexMap)[c->parent_id];
			edgeLength = ComputeNormalizedHammingDistance(&p->sequence,&c->sequence);
			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(p->id,c->id),edgeLength));
		}		
	}	
}

void rootedPhylogeny_tree::WriteEdgeList(string outputFilePrefix){
	ofstream edgeListFile;
	edgeListFile.open(outputFilePrefix + ".edgeList");
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	for (pair<int,rootedPhylogeny_vertex*> vIdVPtrPair : (*this->vertexMap)){
		c = vIdVPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			edgeListFile << p->name;
			edgeListFile << "\t" << c->name;
			edgeListFile << "\t" << this->GetEdgeLength(p->id,c->id) << endl;
			
		}
	}
	edgeListFile.close();
}

void rootedPhylogeny_tree::ComputeNumberOfDescendants(){
	vector <rootedPhylogeny_vertex *> verticesToVisit;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	for (pair<int,rootedPhylogeny_vertex *> vIdVPtrPair : (*this->vertexMap)){
		c = vIdVPtrPair.second;
		if (c->children_id.size() == 0){
			verticesToVisit.push_back(c);
			c->numberOfDescendants = 1;
		}
	}	
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		c = verticesToVisit[numberOfVerticesToVisit-1];		
		p = (*this->vertexMap)[c->parent_id];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;		
		p->timesVisited += 1;
		p->numberOfDescendants += c->numberOfDescendants;
		if (p->timesVisited == int(p->children_id.size()) and p->parent_id != -1){
			verticesToVisit.push_back(p);
			numberOfVerticesToVisit += 1;
		}
	}
}

void rootedPhylogeny_tree::WriteNewickFile(string sequenceFileName){
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;	
	float edgeLength;
	for (pair<int,rootedPhylogeny_vertex*> idAndVertex: *this->vertexMap){
		idAndVertex.second->timesVisited = 0;
		if (idAndVertex.second->children_id.size() == 0){
			idAndVertex.second->newickLabel = idAndVertex.second->name;
			verticesToVisit.push_back(idAndVertex.second);
		} else {
			idAndVertex.second->newickLabel = "";
		}
	}
	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		c = verticesToVisit[numberOfVerticesToVisit-1];	
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;		
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];			
			p->timesVisited += 1;
			edgeLength = GetEdgeLength(p->id, c->id);
			if (p->timesVisited == int(p->children_id.size())){
				p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength) + ")";
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			} else {
				p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
			}			
		}		
	}	
	ofstream newickFile;
	newickFile.open(sequenceFileName+".newick");
	newickFile << this->root->newickLabel << ";" << endl;
	newickFile.close();	
}

//void rootedPhylogeny_tree::WriteNewickFile(string sequenceFileName){
//	vector <rootedPhylogeny_vertex*> verticesToVisit;
//	rootedPhylogeny_vertex * c;
//	rootedPhylogeny_vertex * p;	
//	string newickLabelForTree;
//	float edgeLength;
//	for (pair<int,rootedPhylogeny_vertex*> idAndVertex: *this->vertexMap){
//		idAndVertex.second->timesVisited = 0;
//		if (idAndVertex.second->children_id.size() == 0){
//			idAndVertex.second->newickLabel = idAndVertex.second->name;
//			verticesToVisit.push_back(idAndVertex.second);
//		}
//	}
//	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
//	while (numberOfVerticesToVisit > 0){
//		c = verticesToVisit[numberOfVerticesToVisit-1];	
//		p = (*this->vertexMap)[c->parent_id];
//		verticesToVisit.pop_back();
//		numberOfVerticesToVisit -= 1;
//		if (p->parent_id != -1){
//			p->timesVisited += 1;
////			edgeLength = ComputeNormalizedHammingsDistance(&p->sequence, &c->sequence);
//			edgeLength = GetEdgeLength(p->id, c->id);
//			if (p->timesVisited == int(p->children_id.size())){
//				p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength) + ")";
//				verticesToVisit.push_back(p);
//				numberOfVerticesToVisit += 1;
//			} else {
//				p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
//			}
//		}
//	}	
//	if (p->children_id.size()==3){
//		p->newickLabel += "(";		
//		for (int i=0; i<2; i++){
//			c = (*this->vertexMap)[p->children_id[i]];
//			edgeLength = GetEdgeLength(p->id, c->id);
////			edgeLength = ComputeNormalizedHammingDistance(&p->sequence, &c->sequence);
//			p->newickLabel += c->newickLabel + ":" + to_string(edgeLength) + ",";			
//		}
//		c = (*this->vertexMap)[p->children_id[2]];
//		edgeLength = GetEdgeLength(p->id, c->id);
////		edgeLength = ComputeNormalizedHammingDistance(&p->sequence, &c->sequence);
//		p->newickLabel += c->newickLabel + ":" + to_string(edgeLength) + ");";		
//	} else {
//		c = (*this->vertexMap)[p->children_id[0]];
//		edgeLength = GetEdgeLength(p->id, c->id);
////		edgeLength = ComputeNormalizedHammingDistance(&p->sequence, &c->sequence);
//		p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
//		p->newickLabel += "," + p->name + ":0.00);";		
//	}
//	newickLabelForTree = p->newickLabel;
//	ofstream newickFile;
//	newickFile.open(sequenceFileName+".newick");
//	newickFile << newickLabelForTree << endl;
//	newickFile.close();
//}
void rootedPhylogeny_tree::ContractZeroLengthEdges(){
	// Traverse tree from root to leaves and for each visited vertex p.	
	// Check if edge should be contracted, i.e., at least one of p or c is not an observed vertex.
	// Contract edge (p,v) if length of edge for (p,v) is zero.
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	for (pair<int,rootedPhylogeny_vertex*> idVPtrPair : (*this->vertexMap)){
		if (idVPtrPair.second->parent_id == -1){
			p = idVPtrPair.second;
		}
	}
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	verticesToVisit.push_back(p);
	int numberOfVerticesToVisit = verticesToVisit.size();	
	bool edgeContractedInThisRound = 0;
	float edgeLength;
	while (numberOfVerticesToVisit > 0){
		edgeContractedInThisRound = 0;
		p = verticesToVisit[numberOfVerticesToVisit-1];		
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;		
		vector <int> idsOfChildrenToVisit;		
		idsOfChildrenToVisit = p->children_id;
		for (int child_id: idsOfChildrenToVisit){
			c = (*this->vertexMap)[child_id];
			edgeLength = (*this->edgeLengthsMap)[pair<int,int>(p->id,c->id)];
			if (edgeLength == 0){
				if (IsEdgeContractionFeasbile(p->id,c->id)){
					edgeContractedInThisRound = 1;
					if (p->id < this->numberOfObservedSequences){
						this->ContractEdge(p->id,c->id);
						verticesToVisit.push_back(p);
						numberOfVerticesToVisit += 1;
					} else {
						this->ContractEdge(c->id,p->id);
						verticesToVisit.push_back(c);
						numberOfVerticesToVisit += 1;
					}
				}
			}
			if (edgeContractedInThisRound){
				break;
			}
		}
		if (!edgeContractedInThisRound){
			for (int child_id: idsOfChildrenToVisit){
				c = (*this->vertexMap)[child_id];
				verticesToVisit.push_back(c);
				numberOfVerticesToVisit += 1;				
			}
		}
	}
}

void rootedPhylogeny_tree::SetEdgeLengthsToZeroSuchThatBICIsMinimized() {
	// Traverse tree from root to leaves and for each visited vertex p.	
	// Check if edge should be contracted, i.e., at least one of p or c is not an observed vertex.
	// Contract edge (p,v) if edge contraction reduces BIC
	
	// Check if setting edge length to zero would reduce BIC
	double BIC_opt;
	double BIC_current;
//	double logLik;
	int numberOfThresholdsTried = 0;
	double numberOfParameters;
	this->ComputeLogLikelihoodUsingLongDouble();
	numberOfParameters = this->GetNumberOfNonZeroLengthEdges();
	cout << "number of non zero length edges is " << numberOfParameters << endl;
	float sequenceLength = 0;
	for (int site = 0; site < (int) this->siteWeights.size(); site++){		
		sequenceLength += this->siteWeights[site];
	}
	BIC_opt = (log(sequenceLength) * numberOfParameters ) - ( 2.0 * this->logLikelihood );	
	vector <tuple<float,int,int>> sortedEdges;
	float t;
	int p_id; int c_id;
	for (pair<pair<int,int>,float> edgeIdAndLengthPair : *this->edgeLengthsMap){
		sortedEdges.push_back(make_tuple(edgeIdAndLengthPair.second,edgeIdAndLengthPair.first.first,edgeIdAndLengthPair.first.second));		
	}
	sort(sortedEdges.begin(), sortedEdges.end());
	vector <float> edgeLengthThresholdList;
	for (tuple<float,int,int> edge : sortedEdges){
		tie(t, p_id, c_id) = edge;
		if (t > 0.0 and find(edgeLengthThresholdList.begin(),edgeLengthThresholdList.end(),t) == edgeLengthThresholdList.end()){
			edgeLengthThresholdList.push_back(t);
		}
	}
	vector <tuple<float,int,int>> nonZeroLengthEdgesSelected;
	for (float threshold : edgeLengthThresholdList){
		break;
		numberOfThresholdsTried += 1;
		for (tuple<float,int,int> edge : sortedEdges){
			tie(t, p_id, c_id) = edge;
			if (t > 0 and t < threshold){
				if (p_id < c_id){
					(*this->edgeLengthsMap)[pair<int,int>(p_id, c_id)] = 0;
				} else {
					(*this->edgeLengthsMap)[pair<int,int>(c_id, p_id)] = 0;
				}				
				nonZeroLengthEdgesSelected.push_back(edge);
			}
		}
		this->ComputeLogLikelihoodUsingLongDouble();
		numberOfParameters = this->GetNumberOfNonZeroLengthEdges();
		BIC_current = (log(sequenceLength) * numberOfParameters ) - ( 2.0 * this->logLikelihood );	
		if (BIC_opt < BIC_current){
			for (tuple<float,int,int> edge : nonZeroLengthEdgesSelected){
				tie(t, p_id, c_id) = edge;
				if (p_id < c_id){
					(*this->edgeLengthsMap)[pair<int,int>(p_id, c_id)] = threshold;
				} else {
					(*this->edgeLengthsMap)[pair<int,int>(c_id, p_id)] = threshold;
				}
			}
			break;
			cout << " optimal threshold is " << threshold << endl;
		}		
	}
//	numberOfParameters = (11.0 * (this->numberOfRateCategories -1)) + 3.0;		
//	BIC = ( log(sequenceLength) * numberOfParameters ) - ( 2.0 * this->logLikelihood );	
}



#endif