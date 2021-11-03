void MSTBackbone::MSTBackboneWithFullSEMAndMultipleExternalVertices() {
	vector <string> names;
	vector <vector <unsigned char> > sequences;
	vector <int> sitePatternWeights;
	vector <vector <int> > sitePatternRepetitions;	
	vector <int> idsOfVerticesToRemove;
	vector <int> idsOfVerticesToKeep;
	vector <int> idsOfExternalVertices;
	vector <int> idsOfVerticesForSEM;
	vector <tuple <int, string, vector <unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd;	
	//----##############################################################---//
	//	1.	Initialize the global phylogenetic tree T as the empty graph   //
	//----##############################################################---//
//	cout << "Starting MST-backbone" << endl;
//	cout << "1.	Initialize the global phylogenetic tree T as the empty graph" << endl;
	int numberOfInputSequences = (int) this->M->vertexMap->size();	
	//this->T = new SEM(1);
	current_time = chrono::high_resolution_clock::now();
	timeTakenToComputeEdgeAndVertexLogLikelihoods = chrono::duration_cast<chrono::seconds>(current_time-current_time);
	
	// Initialize global phylogeny
	// idsOfVerticesForSEM.clear();
	// for (pair <int, MST_vertex *> vIdAndPtr : * this->M->vertexMap) {
	// 	idsOfVerticesForSEM.push_back(vIdAndPtr.first);
	// }
	// tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);	
	// this->T->sequenceFileName = this->sequenceFileName;
	// this->T->AddSequences(sequences);
	// this->T->AddNames(names);
	// this->T->AddSitePatternWeights(sitePatternWeights);
	// this->T->SetNumberOfInputSequences(numberOfInputSequences);	
	// this->T->numberOfObservedVertices = numberOfInputSequences;
	
	int largestIdOfVertexInMST = numberOfInputSequences;
	
	bool computeLocalPhylogeneticTree = 1;
	bool numberOfNonSingletonComponentsIsGreaterThanZero = 0;	
	
	while (computeLocalPhylogeneticTree) {
		cout << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
		this->mstBackboneLogFile << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
//		cout << "Max vertex degree in MST is " << this->M->maxDegree << endl;
		//----####################################################################---//
		//	2.	Compute the size of the smallest subtree ts = (Vs,Es) of M s.t.		 //
		//		|Vs| > s. Check if |Vm\Vs| > s.								   		 //
		// 		If yes then go to step 3 else go to step 9					   		 //
		// 		Bootstrapped alignments may contain zero-weight edges		   		 //
		//      If so then replace |Vs| with |{non-zero weighted edges in Es}| 		 //
		//      Additionally replace |Vm\Vs| with |{non-zero weighted edges in Es}|  //
		//----####################################################################---//				
		computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
//		cout << "2. Checking if local phylogenetic tree should be computed" << endl;
		if (computeLocalPhylogeneticTree) {
			//----####################################################################---//
			//	3.	Extract vertices inducing subtree (Vs), and external vertices (Ve)	 //
			//----####################################################################---//				
//			cout << "3. Extract vertices inducing subtree (Vs), and external vertices (Ve)" << endl;
			this->M->SetIdsOfExternalVertices();
			idsOfExternalVertices = this->M->idsOfExternalVertices;
			idsOfVerticesForSEM = this->M->subtree_v_ptr->idsOfVerticesInSubtree;
			this->numberOfVerticesInSubtree = this->M->subtree_v_ptr->idsOfVerticesInSubtree.size();
			for (int id: idsOfExternalVertices) {
				idsOfVerticesForSEM.push_back(id);
			}
			tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
			//----########################################################---//
			//	4.	Compute local phylogeny t over (Vs U Ve) via SEM      	 //
			//----########################################################---//
//			cout << "4.	Compute local phylogeny t over (Vs U Ve) via SEM" << endl;			
			this->t = new SEM(largestIdOfVertexInMST);
			this->t->AddSequences(sequences);
			this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);
			this->t->SetNumberOfInputSequences(numberOfInputSequences);
			this->t->AddRootVertex();
			this->t->AddNames(names);
			this->t->AddGlobalIds(idsOfVerticesForSEM);
			this->t->AddSitePatternWeights(sitePatternWeights);
			this->t->AddSitePatternRepeats(sitePatternRepetitions);			
			this->t->OptimizeTopologyAndParametersOfGMM();			
//			timeTakenToComputeUnrootedPhylogeny += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
			//----##################################################################---//	
			//  5.	Check if # of non-singleton components of forest f in t that       //
			//		is induced by Vs is greater than zero.							   //
			//		i.e., Does local phylogeny contain vertices/edges of interest?	   //
			//----##################################################################---//
			this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
			numberOfNonSingletonComponentsIsGreaterThanZero = this->t->IsNumberOfNonSingletonComponentsGreaterThanZero();			
//			cout << "5. Checking if there are any vertices of interest" << endl;
			if (!numberOfNonSingletonComponentsIsGreaterThanZero) {
				//----####################################################---//	
				//  6.	If no then double subtree size and go to step 2 	 //
				//		else reset subtree size and go to to step 7		     //
				//----####################################################---//		
				this->M->DoubleSubtreeSizeThreshold();
//				cout << "6. Doubling subtree size" << endl;
			} else {
				this->M->ResetSubtreeSizeThreshold();		
				//----################################---//
				//  7.	Add vertices/edges in f to T     //
				//----################################---//
//				cout << "7. Adding vertices/edges in f to T" << endl;
				this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
				//this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();				
				//this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);
				this->t->SetAncestralSequencesString();
				this->ancestralSequencesString += this->t->ancestralSequencesString;				
				t_start_time = chrono::high_resolution_clock::now();
				this->t->SetEdgeAndVertexLogLikelihoods();				
				// this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
				// this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
				t_end_time = chrono::high_resolution_clock::now();
				timeTakenToComputeEdgeAndVertexLogLikelihoods += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
				// Add vertex logLikelihoods
				// Add edge logLikelihoods
				largestIdOfVertexInMST = this->t->largestIdOfVertexInMST;
				//----##############################---//
				//  8.	Update M and go to step 1	   //
				//----##############################---//
//				cout << "8. Updating MST" << endl;
				this->t->SetInfoForVerticesToAddToMST();				
				this->M->UpdateMSTWithMultipleExternalVertices(t->idsOfVerticesToKeepInMST, t->idsOfVerticesToRemove, t->idAndNameAndSeqTuple, idsOfExternalVertices);								
				delete this->t;
			}			
			computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
		}		
//		cout << "CPU time used for computing local phylogeny is " << chrono::duration_cast<chrono::seconds>(t_end_time-t_start_time).count() << " second(s)\n";
//		this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration_cast<chrono::seconds>(t_end_time-t_start_time).count() << " second(s)\n";			
	}	
	//----########################################################---//
	//	9.	Compute phylogenetic tree t over vertices in M, and      //
	//		add vertices/edges in t to T							 //
	//----########################################################---//
	cout << "Computing phylogenetic tree over all vertices in MST" << endl;
//	this->M->UpdateMaxDegree();
//	cout << "Max vertex degree in MST is " << this->M->maxDegree << endl;
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> idPtrPair: * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(idPtrPair.first);
	}
//	cout << "Number of vertices in MST is " << idsOfVerticesForSEM.size() << endl;
//	cout << "Number of edges in MST is " << this->M->edgeWeightsMap.size() << endl;
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
	this->numberOfVerticesInSubtree = sequences.size();
	this->t = new SEM(largestIdOfVertexInMST);
	this->t->SetFlagForFinalIterationOfSEM();
	this->t->AddSequences(sequences);
	this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);
	this->t->SetNumberOfInputSequences(numberOfInputSequences);
	this->t->AddRootVertex();
	this->t->AddNames(names);
	this->t->AddGlobalIds(idsOfVerticesForSEM);
	this->t->AddSitePatternWeights(sitePatternWeights);
	this->t->AddSitePatternRepeats(sitePatternRepetitions);	
	this->t->OptimizeTopologyAndParametersOfGMM();			
//	timeTakenToComputeUnrootedPhylogeny += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
	this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
	this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
	this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();
	this->t->SetAncestralSequencesString();
	this->ancestralSequencesString += this->t->ancestralSequencesString;
	this->t->WriteAncestralSequences();	
	//this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);	
	t_start_time = chrono::high_resolution_clock::now();
	this->t->SetEdgeAndVertexLogLikelihoods();
	//this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
	//this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
	t_end_time = chrono::high_resolution_clock::now();
	timeTakenToComputeEdgeAndVertexLogLikelihoods += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
	delete this->t;
//	this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration_cast<chrono::seconds>(t_end_time-t_start_time).count() << " second(s)\n";		
	// assert that T is a tree
//	cout << "Number of vertices in T is " << this->T->vertexMap->size() << endl;
//	cout << "Number of edges in T is " << this->T->edgeLengths.size() << endl;
	//assert(this->T->vertexMap->size() == this->T->edgeLengths.size() + 1);
	//----##############---//		
	//	10.	Root T via EM  //
	//----##############---//
	current_time = chrono::high_resolution_clock::now();
	timeTakenToComputeGlobalUnrootedPhylogeneticTree = chrono::duration_cast<chrono::seconds>(current_time-start_time);
	timeTakenToComputeGlobalUnrootedPhylogeneticTree -= timeTakenToComputeEdgeAndVertexLogLikelihoods;
	cout << "CPU time used for computing global unrooted phylogenetic tree T is " << timeTakenToComputeGlobalUnrootedPhylogeneticTree.count() << " second(s)\n";
	this->mstBackboneLogFile << "CPU time used for computing global unrooted phylogenetic tree T is " << timeTakenToComputeGlobalUnrootedPhylogeneticTree.count() << " second(s)\n";		
	cout << "Fitting a general Markov (GM) model to T using reconstructed ancestral sequences" << endl;
	this->mstBackboneLogFile << "Fitting a general Markov (GM) model to T using reconstructed ancestral sequences" << endl;		
	//this->T->RootTreeBySumOfExpectedLogLikelihoods();	
	current_time = chrono::high_resolution_clock::now();
	timeTakenToRootViaEdgeLoglikelihoods = chrono::duration_cast<chrono::seconds>(current_time-start_time);
	timeTakenToRootViaEdgeLoglikelihoods -= timeTakenToComputeGlobalUnrootedPhylogeneticTree;
	cout << "CPU time used for fitting a GM model to fully labeled T is " << timeTakenToRootViaEdgeLoglikelihoods.count() << " second(s)\n";
	this->mstBackboneLogFile << "CPU time used for fitting a GM model to fully labeled T is " << timeTakenToRootViaEdgeLoglikelihoods.count() << " second(s)\n";
	//cout << "Log likelihood of fitting a GM model to fully labeled T is " << this->T->maxSumOfExpectedLogLikelihoods << endl;
	//this->mstBackboneLogFile << "Log likelihood of fitting a GM model to fully labeled T is " << this->T->maxSumOfExpectedLogLikelihoods << endl;
	cout << "Writing rooted tree in edge list format" << endl;
	//this->T->WriteRootedTreeAsEdgeList(sequenceFileName + ".edgeList_fullyLabeledRooting");
	//this->T->WriteRootedTreeInNewickFormat(sequenceFileName + ".newick_fullyLabeledRooting");
	//cout << "Rooting T by fitting a GM model using restricted SEM" << endl;
	//this->mstBackboneLogFile << "Rooting T by fitting a GM model using restricted SEM" << endl;
	//this->T->RootTreeByFittingAGMMViaEM();
	// current_time = chrono::high_resolution_clock::now();
	// timeTakenToRootViaRestrictedSEM = chrono::duration_cast<chrono::seconds>(current_time-start_time);	
	// timeTakenToRootViaRestrictedSEM -= timeTakenToComputeGlobalUnrootedPhylogeneticTree;
	// timeTakenToRootViaRestrictedSEM -= timeTakenToRootViaEdgeLoglikelihoods;
	// cout << "CPU time used for rooting leaf-labeled T via restricted SEM is " << timeTakenToRootViaRestrictedSEM.count() << " second(s)\n";
	// this->mstBackboneLogFile << "CPU time used for rooting leaf-labeled T via restricted SEM is " << timeTakenToRootViaRestrictedSEM.count() << " second(s)\n";	
	// cout << "Log likelihood of fitting a GM model to leaf-labeled T via restricted SEM is " << this->T->logLikelihood << endl;
	// this->mstBackboneLogFile << "Log likelihood of fitting a GM model to leaf-labeled T via restricted SEM is " << this->T->logLikelihood << endl;
	// cout << "Writing rooted tree in edge list format and newick format" << endl;
	// this->T->WriteRootedTreeAsEdgeList(sequenceFileName + ".edgeList_leafLabeledRooting");
	// this->T->WriteRootedTreeInNewickFormat(sequenceFileName + ".newick_leafLabeledRooting");
}