
#include "GA.h"
#include "FkModel.h"
#include "SBModel.h"
#include "Mat.h"
#include "Log.h"

#include <iostream>
#include <fstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

FILE* pLogFile;
char strLogSourceName[1024];
char strLogFileName[1024];
char strLog[1024];

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

CGA::CGA()
{
	m_nNumberOfMachines = MPI::COMM_WORLD.Get_size(); // same as MPI_Comm_size(MPI_COMM_WORLD, &nthreads)

	m_pTargetMeasurement1 = new double[Nw_with_border];
	for(int iW = 0; iW < Nw_with_border; ++iW)
	{
		m_pTargetMeasurement1[iW] = 0.0;
	}
	m_pTargetMeasurement2 = new double[Nh_with_border];
	for(int iH = 0; iH < Nh_with_border; ++iH)
	{
		m_pTargetMeasurement2[iH] = 0.0;
	}
	
	m_pTargetFibroblastMat = CreateMat();
	
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

CGA::~CGA()
{
	delete m_pTargetMeasurement1;
	m_pTargetMeasurement1 = NULL;
	
	delete m_pTargetMeasurement2;
	m_pTargetMeasurement2 = NULL;
	
	DestroyMat(m_pTargetFibroblastMat);
	
	delete MinCost;
	MinCost = NULL;
	
	delete Rank;
	Rank = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool ReadIntoArray(ifstream& myfile, double** arr)
{
	int iH = 1;
    int iW = 1;
	char s[10] = {0};
	while(!myfile.eof())
    {
		myfile.getline(s, 10, ' ');
		arr[iH][iW] = atof(s);
		iW++;
		if(iW == Nw+1)
		{           
			iH++;
			if(iH == Nh+1)
			{
				return true;
			}
			iW = 1;
		}
    }
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// We read the rise time matrix of the forward model from the files and
// store it in our measurement vectors.
bool CGA::ReadTargetMeasurements()
{		
	LOG("ReadTargetMeasurements started");
    ifstream myfile;

	// Reading the target fibroblast mat.
	LOG("Opening file TargetFibroblastMat.txt");
	myfile.open("TargetFibroblastMat.txt");
	if(!myfile.is_open())
    {        
    	LOG("Failed to open file TargetFibroblastMat.txt");
        return false;
    }
    if(!ReadIntoArray(myfile, m_pTargetFibroblastMat))
	{
		LOG("Failed to read array out of file TargetFibroblastMat.txt");
        return false;
	}
	myfile.close();
	
	// Use this mat to read the rise time and then measure on the edges.
	double** arr = CreateMat();
	
	// Reading the first protocol measurements.
	LOG("Opening file TargetFibroblastMatResults1.txt");
	myfile.open("TargetFibroblastMatResults1.txt");
	if(!myfile.is_open())
    {        
    	LOG("Failed to open file TargetFibroblastMatResults1.txt");
		DestroyMat(arr);
        return false;
    }
    if(!ReadIntoArray(myfile, arr))
	{
		LOG("Failed to read array out of file TargetFibroblastMatResults1.txt");
		DestroyMat(arr);
        return false;
	}
	//SaveMatToFile(arr, "ReadTargetFibroblastMatResults1.txt");
	for (int iW = 1; iW <= Nw; ++iW)
	{		
		m_pTargetMeasurement1[iW] = arr[Nh-MeasurementMarginIndexes][iW] - arr[MeasurementMarginIndexes + 1][iW];
		//LOG3("Reading target measurements Prot1, (%d,%d): %.3f", Nh-MeasurementMarginIndexes, iW, m_pTargetMeasurement1[iW]);
	}
	myfile.close();

	// Reading the second protocol measurements.
	LOG("Opening file TargetFibroblastMatResults2.txt");
	myfile.open("TargetFibroblastMatResults2.txt");
    if(!myfile.is_open())
    {        
    	LOG("Failed to open file TargetFibroblastMatResults2.txt");
		DestroyMat(arr);
        return false;
    }
    if(!ReadIntoArray(myfile, arr))
	{
		LOG("Failed to read array out of file TargetFibroblastMatResults2.txt");
		DestroyMat(arr);
        return false;
	}
	//SaveMatToFile(arr, "ReadTargetFibroblastMatResults2.txt");
	for (int iH = 1; iH <= Nh; ++iH)
	{		
		m_pTargetMeasurement2[iH] = arr[iH][Nw-MeasurementMarginIndexes] - arr[iH][MeasurementMarginIndexes + 1];
		//LOG3("Reading target measurements Prot2, (%d,%d): %.3f", iH, Nw-MeasurementMarginIndexes, m_pTargetMeasurement2[iH]);
	}
	myfile.close();

	DestroyMat(arr);
	
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::CreateJobs()
{
	LOG("CreateJobs started");
	for (int nPopIndex = 0; nPopIndex < Npop; ++nPopIndex)
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
	LOG1("CreateJobs finished, total number of jobs: %d", m_jobVector.GetSize());
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::ProcessJobs()
{
	LOG("------------- ProcessJobs started -------------");
	
	int nCurProcess = 1;
	int nNumberOfSlaveMachines = (m_nNumberOfMachines-1);
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

void CGA::ProcessResults(Candidate* pCandidate)
{
	// Calculate the cost for this candidate.
	int nCost1 = 0;
	int nCost2 = 0;
	pCandidate->m_cost = 0;
	
	// Calculate cost from result 1.
	for (int iW = 1; iW < Nw+1; iW += SAMPLING_INTERVALS)
	{
		double dCandidateResult = pCandidate->m_pResult1[Nh-MeasurementMarginIndexes][iW] - pCandidate->m_pResult1[MeasurementMarginIndexes + 1][iW];
		int nCandidateResult = (int)ceil(dCandidateResult * 1000);
		int nTargetResult = (int)ceil(m_pTargetMeasurement1[iW] * 1000);
		int nCost = std::abs(nCandidateResult - nTargetResult);
		pCandidate->m_cost += (unsigned long int)nCost;
		nCost1 += (unsigned long int)nCost;
		LOG5("ProcessResults1 at %d, nCandidateResult: %d, nTargetResult: %d, diff: %d, nCost: %d", iW, nCandidateResult, nTargetResult, nCandidateResult - nTargetResult, nCost);
	}
	LOG1("ProcessResults1, total cost for protocol1: ", nCost1);
	
	// Calculate cost from result 2.
	for (int iH = 1; iH < Nh+1; iH += SAMPLING_INTERVALS)
	{
		double dCandidateResult = pCandidate->m_pResult2[iH][Nw-MeasurementMarginIndexes] - pCandidate->m_pResult2[iH][MeasurementMarginIndexes + 1];
		int nCandidateResult = (int)ceil(dCandidateResult * 1000);
		int nTargetResult = (int)ceil(m_pTargetMeasurement2[iH] * 1000);
		int nCost = std::abs(nCandidateResult - nTargetResult);
		pCandidate->m_cost += (unsigned long int)nCost;
		nCost2 += (unsigned long int)nCost;
		LOG5("ProcessResults2 at %d, nCandidateResult: %d, nTargetResult: %d, diff: %d, nCost: %d", iH, nCandidateResult, nTargetResult, nCandidateResult - nTargetResult, nCost);
	}
	LOG1("ProcessResults2, total cost for protocol2: ", nCost2);
	
	LOG3("- ProcessResults, current cost for candidate #%d, %s: %u", pCandidate->m_nIndex, pCandidate->GetFullName(), pCandidate->m_cost);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::ProcessResults()
{
	LOG("ProcessResults started");
	for (int nPopIndex = 0; nPopIndex < Npop; ++nPopIndex)
	{
		Candidate* pCandidate = Population[nPopIndex];
		ProcessResults(pCandidate);
	}
	LOG("ProcessResults ended");
}	
	
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int CGA::ChooseRandomParent()
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

void CGA::CreateRandomPopulation()
{
	for(int iPop = 0; iPop < Npop; ++iPop)
	{		
		Population.push_back(CreateRandomCandidate(iPop));
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate* CGA::CreateRandomCandidate(int nIndex)
{
	int nHStart = rand()%(Max_h_Fibroblast - Min_h_Fibroblast + 1) + Min_h_Fibroblast;
	int nWStart = rand()%(Max_w_Fibroblast - Min_w_Fibroblast + 1) + Min_w_Fibroblast;
	int nHEnd = rand()%(Max_h_Fibroblast - nHStart + 1) + nHStart;
	int nWEnd = rand()%(Max_w_Fibroblast - nWStart + 1) + nWStart;
	return new Candidate(nIndex, nHStart, nWStart, nHEnd, nWEnd);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int CGA::Mutate(int nInValue)
{
	return (nInValue + (rand()%9) - 4);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CGA::IsNewChild(Candidate* pChild)
{
	for(int iPastChild = 0; iPastChild < (int)PastCandidates.size(); iPastChild++)
	{
		if(strcmp(PastCandidates[iPastChild], pChild->GetFullName()) == 0)
		{
			return false;
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::ClearPastCandidates()
{
	for(int iPastChild = 0; iPastChild < (int)PastCandidates.size(); iPastChild++)
	{
		delete PastCandidates[iPastChild];		
	}
	PastCandidates.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::AddChildToPastCandidates(Candidate* pChild)
{
	char* pPastCandidate = new char[CANDIDATE_MAX_NAME_LENGTH];
	memcpy(pPastCandidate, pChild->GetFullName(), CANDIDATE_MAX_NAME_LENGTH);
	PastCandidates.push_back(pPastCandidate);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate* CGA::CreateChild(Candidate* pParent1, Candidate* pParent2, int nIndex)
{
	LOG2("Parent1, Candidate #%d: %s", pParent1->m_nIndex, pParent1->GetFullName());	
	LOG2("Parent2, Candidate #%d: %s", pParent2->m_nIndex, pParent2->GetFullName());	
	int nHStart = Mutate((rand()%2 == 1) ? pParent1->m_nHStart : pParent2->m_nHStart);
	int nWStart = Mutate((rand()%2 == 1) ? pParent1->m_nWStart : pParent2->m_nWStart);
	int nHEnd = Mutate((rand()%2 == 1) ? pParent1->m_nHEnd : pParent2->m_nHEnd);
	int nWEnd = Mutate((rand()%2 == 1) ? pParent1->m_nWEnd : pParent2->m_nWEnd);
	Candidate* pChild = new Candidate(nIndex, nHStart, nWStart, nHEnd, nWEnd);
	LOG2("Child,   Candidate #%d: %s", pChild->m_nIndex, pChild->GetFullName());
	
	// Check if we already had this child.
	if(!IsNewChild(pChild))
	{
		LOG("Already had this child, Don't use this child and try again");
		delete pChild;
		return NULL;
	}		
	LOG("-");
	AddChildToPastCandidates(pChild);
	return pChild;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int CGA::CalculateError(double** pBestMatch)
{
	int nError = 0;
	for(int iH = 1; iH < Nh+1; ++iH)
	{
		for(int iW = 1; iW < Nw+1; ++iW)
		{
			int nTarget = (int)ceil(m_pTargetFibroblastMat[iH][iW]);
			int nMatch = (int)ceil(pBestMatch[iH][iW]);
			if((nTarget == 0) != (nMatch == 0))
			{
				nError++;
			}
		}
	}
	return nError;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double CGA::CalculateTargetCoverage(double** pBestMatch)
{
	double dTargetSize = 0.0;
	double dCoverage = 0.0;
	for(int iH = 1; iH < Nh+1; ++iH)
	{
		for(int iW = 1; iW < Nw+1; ++iW)
		{
			int nTarget = (int)ceil(m_pTargetFibroblastMat[iH][iW]);
			int nMatch = (int)ceil(pBestMatch[iH][iW]);
			if(nTarget != 0)
			{
				//printf("Target=%d at iH=%d, iW=%d\n", nTarget, iH, iW);
				dTargetSize++;
				if(nMatch != 0)
				{
					//printf("nMatch=%d at iH=%d, iW=%d\n", nTarget, iH, iW);
					dCoverage++;
				}				
			}
		}
	}
	if(dCoverage == 0.0) return 0.0;
	return 100.0*(dCoverage/dTargetSize);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::RunGA()
{
	LOG("Init GA...");
	if(!ReadTargetMeasurements())
	{
		LOG("Failed to read the target measurements. Aborting.");
		return;
	}
	
	// Create the starting population of candidates.
	CreateRandomPopulation();
	
	for(int iPop = 0; iPop < Npop; ++iPop)
	{
		AddChildToPastCandidates(Population[iPop]);		
	}
	
	// Create min cost file.
	char minCostFileName[1024] = {0};
	sprintf(minCostFileName, "%s/MinCost.txt", LOG_FOLDER);
	FILE* pMinCostFile = fopen(minCostFileName, "w");
	if(pMinCostFile != NULL)
	{
		fprintf(pMinCostFile, "Iteration | MinCost | Mat | Error | Coverage\n");
		fclose(pMinCostFile);
	}
	
	// Start the loop.
	LOG("------------------");
	int nLastMinCostCounter = 0;
	int nCurIteration = 0;
	while(nCurIteration <= MaxIterations)
	{
		LOG1("Starting iteration #%d", nCurIteration);
		LOG("------------------");
		clock_t iterationStartingTime = clock();
	
		LOG("Candidates for this iteration:");
		for(int iPop = 0; iPop < Npop; ++iPop)
		{			
			LOG2("Candidate #%d, %s", Population[iPop]->m_nIndex, Population[iPop]->GetFullName());
		}
	
		// Process candidates - Execute the simulation with the backward model
		// for each of candidates.
		CreateJobs();
		ProcessJobs();
		
		// Go over the results of the simulations and determine 
		// the cost of each of candidates.
		ProcessResults();

		// Sort candidates according to the cost.
		LOG("--");
		LOG("Sorting the population...");
		std::sort(Population.begin(), Population.end(), CandidateCompare);
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			Population[iPop]->m_nIndex = iPop;
			LOG3("ProcessResults, current cost for candidate #%d, %s: %d", Population[iPop]->m_nIndex, Population[iPop]->GetFullName(), Population[iPop]->m_cost);
		}
		
		// Write output files (for presentation and debugging).
		unsigned long int nTotalCost = 0;
		for(int iPop = 0; iPop < Npop; ++iPop)
		{	
			if(iPop == 0) // To save space of txt files...
			{
				char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
				sprintf(riseTimeCandidateFileName, "RiseTime_Itr%d_Indx%d_Prot1.txt", nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult1, riseTimeCandidateFileName);
				sprintf(riseTimeCandidateFileName, "RiseTime_Itr%d_Indx%d_Prot2.txt", nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult2, riseTimeCandidateFileName);
			}
			else
			{
				char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
				sprintf(riseTimeCandidateFileName, "Extra/RiseTime_Itr%d_Indx%d_Prot1.txt", nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult1, riseTimeCandidateFileName);
				sprintf(riseTimeCandidateFileName, "Extra/RiseTime_Itr%d_Indx%d_Prot2.txt", nCurIteration, iPop);
				SaveMatToFile(Population[iPop]->m_pResult2, riseTimeCandidateFileName);
			}
			nTotalCost += Population[iPop]->m_cost;
		}
		double dAvgCost = nTotalCost/Npop;

		// Update the min cost file.
		MinCost[nCurIteration] = Population[0]->m_cost;
		LOG1("Avg Cost: %.3f", dAvgCost);
		LOG2("Min Cost: %d for %s", MinCost[nCurIteration], Population[0]->GetFullName());
		int nError = CalculateError(Population[0]->m_pFibroblastMat);
		LOG1("The (geomatric) error for the best match is: %d", nError);
		double dCoverage = CalculateTargetCoverage(Population[0]->m_pFibroblastMat);
		LOG1("The target coverage for the best match is: %.3f%%", dCoverage);
		pMinCostFile = fopen(minCostFileName, "a");
		if(pMinCostFile != NULL)
		{
			fprintf(pMinCostFile, "%d | %d | %s | %d | %.3f\n", nCurIteration, MinCost[nCurIteration], Population[0]->GetFullName(), nError, dCoverage);
			fclose(pMinCostFile);
		}
		
		// Check break conditions.
		if(MinCost[nCurIteration] == 0)
		{
			LOG("Reached zero cost. Breaking the iterations !");
			break;
		}
		if(nCurIteration != 0)
		{
			if(MinCost[nCurIteration] == MinCost[nCurIteration-1])
			{
				nLastMinCostCounter++;				
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
		++nCurIteration;
		
		clock_t iterationEndingTime = clock();
		double iterationRunningTime = (iterationEndingTime - iterationStartingTime)/double(CLOCKS_PER_SEC);
		LOG1("Iteration start: %d", iterationStartingTime);
		LOG1("Iteration end: %d", iterationEndingTime);
		LOG1("Iteration duration: %.3f seconds", iterationRunningTime);
				
		// Create the next generation.
		LOG("--");
		LOG("Creating the next generation...");
		LOG1("Need to create %d offsprings in total", Nmates);
		LOG("--");
		int iOffspring = 0;
		while(iOffspring < Nmates)
		{	
			int nFirstParent = 0;
			int nSecondParent = 0;
			if(iOffspring != 0)
			{
				while(nFirstParent == nSecondParent)
				{
					nFirstParent = ChooseRandomParent();
					nSecondParent = ChooseRandomParent();
				}				
			}
			
			int nNewCandidateIndex = NsurvivingPopulation + iOffspring;			
			delete Population[nNewCandidateIndex];
			LOG4("Creating offspring %d (candidate #%d) from parents %d and %d", iOffspring+1, nNewCandidateIndex, nFirstParent, nSecondParent);
			Population[nNewCandidateIndex] = CreateChild(Population[nFirstParent], Population[nSecondParent], nNewCandidateIndex);
			if(Population[nNewCandidateIndex] != NULL) // Only if we got a new child we can continue.
			{
				++iOffspring;
			}
		}
    		         
		//CreateMutations();
		LOG("------------------");
	}
	
	// Finished executing the GA.
	LOG2("The best match found is: %s, with cost: ", Population[0]->GetFullName(), MinCost[nCurIteration]);
	int nError = CalculateError(Population[0]->m_pFibroblastMat);
	LOG1("The error for this match is: %d", nError);
	
	ClearPastCandidates();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartMainProcess()
{
	LOG("Starting main process");
	
	// This is for the controlling master process.
	srand((unsigned int)time(NULL));
	
	// Start timing.
	clock_t mainStartingTime = clock();
		
	// Run main program.
	CGA ga;
	ga.RunGA();

	int nNumberOfMachines = MPI::COMM_WORLD.Get_size(); // same as MPI_Comm_size(MPI_COMM_WORLD, &nthreads)
	for(int iProcess = 1; iProcess < nNumberOfMachines; ++iProcess)
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

void StartSlaveProcess()
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
		//CFkModel* pModel = new CFkModel();
		CSBModel* pModel = new CSBModel();
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

#include "SBModelDefs.h"

void CGA::Test()
{
	//for(double dDiff = -0.0000020; dDiff <= 0.0000020; dDiff += 0.0000001)
	//for(double dj = -0.01; dj <= 0.01; dj += 0.0005)
	{

	LOG("----------------------- Test started -----------------------");
	
	int nIndex = 0;
	//Candidate* pCandidate = CreateRandomCandidate(nIndex);
	//Candidate* pCandidate = new Candidate(nIndex, 0, 0, 0, 0); // No mat.
	Candidate* pCandidate = new Candidate(nIndex, 50, 40, 64, 64);
	//Candidate* pCandidate = new Candidate(nIndex, 58,49,68,66);
	
	SaveMatToFile(pCandidate->m_pFibroblastMat, "TestMat.txt");
	LOG1("Testing cost with candidate: %s", pCandidate->GetFullName());
	S1Protocol s1;
	S2Protocol s2;
	
	CSBModel* pModel = new CSBModel(); // = new CFkModel();
	//pModel->SetDiffusion(Diffusion + dDiff);
	//LOG1("SetDiffusion=%.8f", Diffusion + dDiff);
	//pModel->SetJ(j_var + dj);
	//LOG1("SetJ=%.5f", j_var + dj);
	
	clock_t startingTime;
	clock_t endingTime;
	double runningTime = 0.0;
	
	LOG("Executing 1st protocol");
	startingTime = clock();
	pModel->ExecuteModel(pCandidate->m_pFibroblastMat, pCandidate->m_pResult1, s1, "/a/home/cc/students/enginer/shaishif/Logs/Output");	
	SaveMatToFile(pCandidate->m_pResult1, "TestMatResults1.txt");
	endingTime = clock();
	runningTime = (endingTime - startingTime)/double(CLOCKS_PER_SEC);	
	LOG1("Finished executing 1st protocol after %.3f seconds", runningTime);
	
	//if(false)
	{
		LOG("Executing 2nd protocol");
		startingTime = clock();
		pModel->ExecuteModel(pCandidate->m_pFibroblastMat, pCandidate->m_pResult2, s2);	
		SaveMatToFile(pCandidate->m_pResult2, "TestMatResults2.txt");
		endingTime = clock();
		runningTime = (endingTime - startingTime)/double(CLOCKS_PER_SEC);	
		LOG1("Finished executing 2nd protocol after %.3f seconds", runningTime);
	}	
	
	// Compare results to the target.
	//if(false)
	{
		if(ReadTargetMeasurements())
		{
			ProcessResults(pCandidate);				
			int nError = CalculateError(pCandidate->m_pFibroblastMat);
			LOG1("The (geomatric) error for this match is: %d", nError);
			
			double dCoverage = CalculateTargetCoverage(pCandidate->m_pFibroblastMat);
			LOG1("The target coverage for this match is: %.3f%%", dCoverage);	
		}
		else
		{
			LOG("Failed to read the target measurements. Aborting.");
		}
	}

	delete pModel;
	delete pCandidate;
	
	LOG("----------------------- Test ended -----------------------");
	
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CreateLogFiles()
{
	// Getting general info about MPI.
	int nCpuNameLen = MPI_MAX_PROCESSOR_NAME;
	char sMachineName[MPI_MAX_PROCESSOR_NAME] = {0};
	MPI_Get_processor_name(sMachineName, &nCpuNameLen);
	int nCurProcess = MPI::COMM_WORLD.Get_rank(); // same as MPI_Comm_rank(MPI_COMM_WORLD, &tid)
		
	// Create the log file name.
	sprintf(strLogFileName, "%s/Log_%i_on_%s.txt", LOG_FOLDER, nCurProcess, sMachineName);
	sprintf(strLogSourceName, "Process %i on %s | ", nCurProcess, sMachineName);
	
	// Clear the current log file.
	pLogFile = fopen(strLogFileName, "w");
	if(pLogFile != NULL) { fclose(pLogFile); }	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

//#define TESTING

int main(int argc, char *argv[])
{
	int nProvided = -1;
	int nMpiIntRet = MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &nProvided);	
	srand((unsigned int)time(NULL));	
		
	// Start the appropriate process: master/slave
	if(MPI::COMM_WORLD.Get_rank() == MPI_MASTER)
	{	
		CreateLogFiles();
	#ifdef TESTING
		CGA ga;
		ga.Test();
	}
	#else
		StartMainProcess();	
	}
	else
	{	
		CreateLogFiles();
		StartSlaveProcess();
	}
	LOG("Ending process");
	#endif // TESTING
		
	MPI_Finalize();
		
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
