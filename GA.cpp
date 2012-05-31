
#include <pthread.h>
#include "GA.h"
#include "FkModel.h"
#include "Mat.h"
#include "Log.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

FILE* pLogFile;
char strLogSourceName[1024];
char strLogFileName[1024];
char strLog[1024];

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Ga ga;
int Ga::nNumberOfMachines = 0;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Ga::Ga()
{
	m_nCurIteration = 0;
	m_pTargetMeasurement1 = new double[Nw];
	m_pTargetMeasurement2 = new double[Nh];
	MinCost = NULL;
	Rank = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::InitGa()
{		
	for(int iPop = 0; iPop < Npop; ++iPop)
	{		
		Candidate* pCandidate = new Candidate(iPop);
		Population.push_back(pCandidate);
	}

	MinCost = new unsigned long int[MaxIterations];
	for(int iteration=0; iteration < MaxIterations; ++iteration)
	{
		MinCost[iteration] = 0;
	}
	
	int nRankSum = 0;
	for(int iRankSum = 1; iRankSum <= NsurvivingPopulation; ++iRankSum)
	{
		nRankSum += iRankSum;
	}

	Rank = new double[NsurvivingPopulation];
	for (int iRank = 1 ; iRank <= NsurvivingPopulation; ++iRank)
	{
		Rank[iRank-1] = double(NsurvivingPopulation - iRank + 1)/double(nRankSum);
		if(iRank != 1)
		{
			Rank[iRank-1] = Rank[iRank-1] + Rank[iRank-2];
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateTargetMeasurements()
{
	LOG("CreateTargetMeasurements");
	MPI_Status status;
	Candidate* pCandidate = new Candidate(-1);
	SaveMatToFile(pCandidate->m_pFibroblastMat, "TargetFibroblastMat.txt");
	
	// Execute the simulation for all protocols.
	CFkModel* pModel = new CFkModel();
	LOG("Executing 1st protocol");
	char modelOutput[1024] = {0};
	sprintf(modelOutput, "%s/Output", LOG_FOLDER);
	S1Protocol s1;			
	pModel->ExecuteModel(pCandidate->m_pFibroblastMat, pCandidate->m_pResult1, s1, modelOutput);
	LOG("Executing 2nd protocol");
	S2Protocol s2;
	pModel->ExecuteModel(pCandidate->m_pFibroblastMat, pCandidate->m_pResult2, s2);
	
	// Save the results to our log files.
	LOG("Saving results to file");
	SaveMatToFile(pCandidate->m_pResult1, "TargetFibroblastMatResults1.txt");
	SaveMatToFile(pCandidate->m_pResult2, "TargetFibroblastMatResults2.txt");
	
	// Save the measurements.
	for (int iW = 1; iW < Nw+1; ++iW)
	{
		m_pTargetMeasurement1[iW-1] = pCandidate->m_pResult1[Nh][iW];
	}
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		m_pTargetMeasurement2[iH-1] = pCandidate->m_pResult2[iH][Nw];
	}
	
	delete pCandidate;
	pCandidate = NULL;
	delete pModel;
	pModel = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateJobs()
{
	LOG("CreateJobs started");
	for (int nPopIndex = 0; nPopIndex < Npop; ++nPopIndex)
	{			
		// We don't need to recalculate the first candidate because
		// it hasn't been mutated.
		// bool bAddJob = (Population[nPopIndex]->m_nIndex != 0);
		
		// But we do need to calculate its cost if this is 
		// the first iteration ever.
		// bAddJob |= (m_nCurIteration == 0);
		
		// if(bAddJob)
		{	
			Job* pJob1 = new Job();
			pJob1->m_pCandidate = Population[nPopIndex];
			pJob1->m_nJobType = MPI_JOB_1_TAG;
			pJob1->m_pResultsMat = pJob1->m_pCandidate->m_pResult1;
			m_jobVector.AddJob(pJob1);
			
			Job* pJob2 = new Job();
			pJob2->m_pCandidate = Population[nPopIndex];
			pJob2->m_nJobType = MPI_JOB_2_TAG;
			pJob2->m_pResultsMat = pJob2->m_pCandidate->m_pResult2;
			m_jobVector.AddJob(pJob2);
		}
	}
	LOG1("CreateJobs finished, total number of jobs: %d", m_jobVector.GetSize());
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::ProcessJobs()
{
	LOG("------------- ProcessJobs started -------------");
	
	int nCurProcess = 1;
	int nNumberOfSlaveMachines = (nNumberOfMachines-1);
	MPI_Request* requestArray = new MPI_Request[nNumberOfSlaveMachines];
	MPI_Status* statusArray = new MPI_Status[nNumberOfSlaveMachines];

	while(true)
	{	
		Job* pJob = m_jobVector.GetJob();
		if(pJob == NULL)
		{
			LOG("ProcessJobs finished");
			break;
		}
		
		LOG2("ProcessJobs, processing job - Candidate: #%d, JobType: %d", pJob->m_pCandidate->m_nIndex, pJob->m_nJobType);
		LOG1("ProcessJobs, Sending job to process #%d", nCurProcess);
		MPI_Send(&(pJob->m_pCandidate->m_pFibroblastMat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, nCurProcess, pJob->m_nJobType, MPI::COMM_WORLD);
		LOG1("ProcessJobs, Async wait for job from process #%d", nCurProcess);
		MPI_Irecv(&(pJob->m_pResultsMat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, nCurProcess, MPI_RESULT_TAG, MPI::COMM_WORLD, &requestArray[nCurProcess-1]);
		
		if(nCurProcess == (nNumberOfSlaveMachines))
		{
			LOG("ProcessJobs, Wait for all pending jobs to finish");
			MPI_Waitall(nNumberOfSlaveMachines, requestArray, statusArray);
			nCurProcess = 1;
			LOG("ProcessJobs, Wait is over, continue processing jobs");
		}
		else
		{
			nCurProcess++;
		}
		
		delete pJob;
		pJob == NULL;
	}
	
	LOG("------------- ProcessJobs ended -------------");
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::ProcessResults()
{
	LOG("ProcessResults started");
	for (int nPopIndex = 0; nPopIndex < Npop; ++nPopIndex)
	{
		Candidate* pCandidate = Population[nPopIndex];
				
		// Calculate the cost for this candidate.
		pCandidate->m_cost = 0;
		
		// Calculate cost from result 1.
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			int nCandidateResult = (int)ceil(pCandidate->m_pResult1[Nh][iW]);
			int nTargetResult = (int)ceil(m_pTargetMeasurement1[iW]);
			pCandidate->m_cost += (unsigned long int)std::abs(nCandidateResult - nTargetResult);
		}

		// Calculate cost from result 2.
		for (int iH = 1; iH < Nh+1; ++iH)
		{
			int nCandidateResult = (int)ceil(pCandidate->m_pResult2[iH][Nw]);
			int nTargetResult = (int)ceil(m_pTargetMeasurement2[iH]);
			pCandidate->m_cost += (unsigned long int)std::abs(nCandidateResult - nTargetResult);
		}
		
		LOG3("ProcessResults, current cost for candidate #%d, %s: %d", pCandidate->m_nIndex, pCandidate->GetFullName(), pCandidate->m_cost);
	}
}	
	
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int Ga::GetMate()
{
	double parentRand = (double)rand()/(double)RAND_MAX;        
    for (int iRank = 0 ; iRank < (NsurvivingPopulation-1); ++iRank)
	{
        if(Rank[iRank] >= parentRand)
		{
            return iRank;
		}
	}
	return (NsurvivingPopulation-1);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateChild(Candidate* pParent1, Candidate* pParent2, Candidate* pChild)
{
	if((pChild == pParent1) || (pChild == pParent2))
	{
		throw;
	}
		
	LOG2("Parent1, Candidate #%d: %s", pParent1->m_nIndex, pParent1->GetFullName());	
	LOG2("Parent2, Candidate #%d: %s", pParent2->m_nIndex, pParent2->GetFullName());

	pChild->m_nCenterH = pParent1->m_nCenterH;
	pChild->m_nCenterW = pParent1->m_nCenterW;
	pChild->m_nHeight = pParent1->m_nHeight;
	pChild->m_nWidth = pParent1->m_nWidth;
	
	int nCrossPoint1 = rand()%4;
	switch(nCrossPoint1)
	{
		case 0:
			pChild->m_nCenterH = pParent2->m_nCenterH;	
			break;
		case 1:
			pChild->m_nCenterW = pParent2->m_nCenterW;
			break;
		case 2:
			pChild->m_nHeight = pParent2->m_nHeight;			
			break;
		case 3:
			pChild->m_nWidth = pParent2->m_nWidth;
			break;
	}

	int nCrossPoint2 = rand()%4;
	switch(nCrossPoint2)
	{
		case 0:
			pChild->m_nCenterH = pParent2->m_nCenterH;	
			break;
		case 1:
			pChild->m_nCenterW = pParent2->m_nCenterW;
			break;
		case 2:
			pChild->m_nHeight = pParent2->m_nHeight;			
			break;
		case 3:
			pChild->m_nWidth = pParent2->m_nWidth;
			break;
	}

	pChild->CreateFibroblasts();
	LOG2("Child,   Candidate #%d: %s", pChild->m_nIndex, pChild->GetFullName());	
	LOG("-");
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool Ga::FindSimilarCandidate(int nIndex)
{
	for(int iPop = 0; iPop < Npop; ++iPop)
	{
		if(iPop != nIndex)
		{
			if(Population[iPop]->m_nCenterH == Population[nIndex]->m_nCenterH &&
			   Population[iPop]->m_nCenterW == Population[nIndex]->m_nCenterW &&
			   Population[iPop]->m_nHeight == Population[nIndex]->m_nHeight &&
			   Population[iPop]->m_nWidth == Population[nIndex]->m_nWidth)
			{
				LOG2("Candiates are identical, #%d & #%d", nIndex, iPop);
				return true;
			}
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateMutations()
{
	LOG("Starting mutations...");
	LOG("--");	
	for(int iPop = 1; iPop < Npop; ++iPop)
	{
		LOG2("Candiate #%d before mutation: %s", iPop, Population[iPop]->GetFullName());
		Population[iPop]->Mutate();
		while(FindSimilarCandidate(iPop))
		{
			LOG2("Candiate #%d before unmutation: %s", iPop, Population[iPop]->GetFullName());
			Population[iPop]->UnMutate();
			LOG2("Candiate #%d after unmutation : %s", iPop, Population[iPop]->GetFullName());
			Population[iPop]->Mutate();
		}
		Population[iPop]->CreateFibroblasts();
		Population[iPop]->m_cost = 0;
		LOG2("Candiate #%d after mutation : %s", iPop, Population[iPop]->GetFullName());
		LOG("-");
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::RunGa()
{
	LOG("Init GA...");
	InitGa();

	// Create and measure our target mat.
	CreateTargetMeasurements();
	
	// Create min cost file
	char minCostFileName[1024] = {0};
	sprintf(minCostFileName, "%s/MinCost.txt", LOG_FOLDER);
	FILE* pMinCostFile = fopen(minCostFileName, "w");
	if(pMinCostFile != NULL)
	{
		fprintf(pMinCostFile, "Iteration | MinCost | Mat\n");
		fclose(pMinCostFile);
	}
	
	LOG("------------------");
	int nLastMinCostCounter = 0;
	while(m_nCurIteration <= MaxIterations)
	{
		LOG1("Starting iteration #%d", m_nCurIteration);
		LOG("------------------");
		clock_t iterationStartingTime = clock();
	
		// Determine the cost function for each member of the population.
		CreateJobs();
		ProcessJobs();
		ProcessResults();

		// Sort according to the cost and then mate according to the rank.
		LOG("Sorting the population...");
		std::sort(Population.begin(), Population.end(), CandidateCompare);
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			Population[iPop]->m_nIndex = iPop;
			LOG3("ProcessResults, current cost for candidate #%d, %s: %d", Population[iPop]->m_nIndex, Population[iPop]->GetFullName(), Population[iPop]->m_cost);
		}
		
		unsigned long int nTotalCost = 0;
		for(int iPop = 0; iPop < Npop; ++iPop)
		{	
			if(iPop == 0) // To save space of txt files...
			{
				char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
				sprintf(riseTimeCandidateFileName, "RiseTime_Itr%d_Indx%d_Prot1.txt", m_nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult1, riseTimeCandidateFileName);
				sprintf(riseTimeCandidateFileName, "RiseTime_Itr%d_Indx%d_Prot2.txt", m_nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult2, riseTimeCandidateFileName);
			}
			else
			{
				char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
				sprintf(riseTimeCandidateFileName, "Extra/RiseTime_Itr%d_Indx%d_Prot1.txt", m_nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult1, riseTimeCandidateFileName);
				sprintf(riseTimeCandidateFileName, "Extra/RiseTime_Itr%d_Indx%d_Prot2.txt", m_nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult2, riseTimeCandidateFileName);
			}
			nTotalCost += Population[iPop]->m_cost;
		}
		double dAvgCost = nTotalCost/Npop;

		// Create the next generation.		
		LOG("--");
		LOG("Creating the next generation...");
		LOG1("Need to create %d offspring in total", Nmates);
		LOG("--");
		int iOffspring = 0;
		while(iOffspring < Nmates)
		{	
			if(iOffspring == 0)
			{
				LOG("The first offspring will be a (mutated) clone of the best candidate");				
				CreateChild(Population[0], Population[0], Population[NsurvivingPopulation]);
				++iOffspring;
			}
			else
			{
				int nFirstParent = GetMate();
				int nSecondParent = GetMate();
				if(nFirstParent != nSecondParent)
				{
					LOG1("Creating offspring %d", iOffspring+1);
					Candidate* Parent1 = Population[nFirstParent];
					Candidate* Parent2 = Population[nSecondParent];
				
					// Choose a cross over point.
					// Add the new offspirng to the population.
					CreateChild(Parent1, Parent2, Population[NsurvivingPopulation + iOffspring]);
					++iOffspring;
				}
			}
		}
    		         
		CreateMutations();		
		
		MinCost[m_nCurIteration] = Population[0]->m_cost;
		LOG1("Avg Cost: %.3f", dAvgCost);
		LOG2("Min Cost: %d for %s", MinCost[m_nCurIteration], Population[0]->GetFullName());
		pMinCostFile = fopen(minCostFileName, "a");
		if(pMinCostFile != NULL)
		{
			fprintf(pMinCostFile, "%d | %d | %s\n", m_nCurIteration, MinCost[m_nCurIteration], Population[0]->GetFullName());
			fclose(pMinCostFile);
		}
		if(MinCost[m_nCurIteration] == 0)
		{
			LOG("Reached zero cost. Breaking the iterations !");
			break;
		}
		if(m_nCurIteration != 0)
		{
			if(MinCost[m_nCurIteration] == MinCost[m_nCurIteration-1])
			{
				nLastMinCostCounter++;
				/*
				if(nLastMinCostCounter >= 5)
				{
					LOG("We're starting to reach a dead end in the costs. Creating an additional clone of the best candidate !");					
					CreateChild(Population[0], Population[0], Population[Npop-1]);
					// Need to mutate and create mat...
				}
				*/
				if(nLastMinCostCounter == 10)
				{
					LOG("Reached a dead end in the costs. Breaking the iterations !");
					break;
				}				
			}
			else
			{
				nLastMinCostCounter = 0;
			}			
		}
		++m_nCurIteration;
		
		clock_t iterationEndingTime = clock();
		double iterationRunningTime = (iterationEndingTime - iterationStartingTime)/double(CLOCKS_PER_SEC);
		LOG1("Iteration start: %d", iterationStartingTime);
		LOG1("Iteration end: %d", iterationEndingTime);
		LOG1("Iteration duration: %.3f seconds", iterationRunningTime);
		LOG("------------------");
	}	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartMainProcess(char* sMachineName)
{
	LOG("Starting main process");
	
	// This is for the controlling master process.
	srand((unsigned int)time(NULL));
	
	// Start timing.
	clock_t mainStartingTime = clock();
		
	// Run main program.
	ga.RunGa();

	for(int iProcess = 1; iProcess < Ga::nNumberOfMachines; ++iProcess)
	{
		MPI_Send(NULL, 0, MPI_DOUBLE, iProcess, MPI_DIE_TAG, MPI::COMM_WORLD);
		LOG1("Sending quit flag to process %d", iProcess);
	}
	
	// End timing.
	clock_t mainEndingTime = clock();
	double mainRunningTime = (mainEndingTime - mainStartingTime)/double(CLOCKS_PER_SEC);
	LOG1("Main duration: %.3f seconds", mainRunningTime);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartSlaveProcess(int nProcess, char* sMachineName)
{
	LOG("Starting slave process");
		
	// This is for all the other processes which are not the master.		
	// Create all the vars we will use for data transfer.
	MPI_Status status; // can save resources by using the predefined constant MPI_STATUS_IGNORE as a special value for the status argument.	
	S1Protocol s1;
	S2Protocol s2;
	
	// Start the infinite loop (until the master tells us to quit).
	LOG("Starting loop");
	while(true)
	{			
		CFkModel* pModel = new CFkModel();
		double** fibroblast_mat = CreateMat();
		double** result_mat = CreateMat();
		
		LOG("Waiting to receive mat");		
		int nRet = MPI_Recv(&(fibroblast_mat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, MPI_MASTER, MPI_ANY_TAG, MPI::COMM_WORLD, &status);
		LOG1("Mat received, res: %d", nRet);
		
		// Start timing.
		clock_t startingTime = clock();
					
		if (status.MPI_TAG == MPI_JOB_1_TAG) 
		{
			LOG("Executing 1st protocol");
			pModel->ExecuteModel(fibroblast_mat, result_mat, s1);
		}
		else if (status.MPI_TAG == MPI_JOB_2_TAG) 
		{
			LOG("Executing 2nd protocol");
			pModel->ExecuteModel(fibroblast_mat, result_mat, s2);
		}
		else if (status.MPI_TAG == MPI_DIE_TAG) 
		{
			LOG("Got the die tag. Exiting...");
			break;
		}
		else
		{
			LOG("Got an invalid tag. Aborting !");
			break;
		}
		
		// End timing.
		clock_t endingTime = clock();
		double runningTime = (endingTime - startingTime)/double(CLOCKS_PER_SEC);
		
		LOG1("Finished executing protocol after %.3f seconds, sending results.", runningTime);
		MPI_Send(&(result_mat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, MPI_MASTER, MPI_RESULT_TAG, MPI::COMM_WORLD);
		LOG("Results were sent.");
		
		// Clear up the matrix we use for data transfer.
		delete pModel;
		pModel = NULL;
		DestroyMat(fibroblast_mat);
		DestroyMat(result_mat);
		
	} // End of loop.
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	int nProvided = -1;
	int nMpiIntRet = MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &nProvided);	
	srand((unsigned int)time(NULL));	
	
	// Getting general info about MPI.
	int nCpuNameLen = MPI_MAX_PROCESSOR_NAME;
	char sMachineName[MPI_MAX_PROCESSOR_NAME] = {0};
	MPI_Get_processor_name(sMachineName, &nCpuNameLen);
	int nCurProcess = MPI::COMM_WORLD.Get_rank(); // same as MPI_Comm_rank(MPI_COMM_WORLD, &tid)
	Ga::nNumberOfMachines = MPI::COMM_WORLD.Get_size(); // same as MPI_Comm_size(MPI_COMM_WORLD, &nthreads)
	
	// Create the log file name.
	sprintf(strLogFileName, "%s/Log_%i_on_%s.txt", LOG_FOLDER, nCurProcess, sMachineName);
	sprintf(strLogSourceName, "Process %i on %s | ", nCurProcess, sMachineName);
	
	// Clear the current log file.
	pLogFile = fopen(strLogFileName, "w");
	if(pLogFile != NULL) { fclose(pLogFile); }
	
	// Start the appropriate process: master/slave
	if(nCurProcess == MPI_MASTER)
	{		
		StartMainProcess(sMachineName);		
	}
	else
	{		
		StartSlaveProcess(nCurProcess, sMachineName);
	}
	
	LOG("Ending process");
	MPI_Finalize();
		
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
