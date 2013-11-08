
#include "GA.h"
//#include "FkModel.h"//
#include "SBModel.h"
#include "SBModelDefs.h"
#include "Mat.h"
#include "Log.h"

#include <iostream>
#include <fstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

CGA::CGA()
{
	bLogToFileOnly = false;
	m_nNumberOfMachines = MPI::COMM_WORLD.Get_size(); // same as MPI_Comm_size(MPI_COMM_WORLD, &nthreads)
	m_nNumberOfSlaveMachines = (m_nNumberOfMachines-1);

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
	delete [] m_pTargetMeasurement1;
	m_pTargetMeasurement1 = NULL;

	delete [] m_pTargetMeasurement2;
	m_pTargetMeasurement2 = NULL;

	DestroyMat(m_pTargetFibroblastMat);

	delete [] MinCost;
	MinCost = NULL;

	delete [] Rank;
	Rank = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#define MAX_FILE_ELEMENT_LENGTH 10

bool ReadIntoArray(ifstream& myfile, double** arr)
{
	int iH = 1;
    int iW = 1;
	char s[MAX_FILE_ELEMENT_LENGTH] = {0};
	while(!myfile.eof())
    {
		myfile.getline(s, MAX_FILE_ELEMENT_LENGTH, ' ');
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
	//SaveMatToFile(m_pTargetFibroblastMat, "ReadTargetFibroblastMat.txt");

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
		Job job1;
		job1.m_pCandidate = m_population[nPopIndex];
		job1.m_nJobType = MPI_JOB_1_TAG;
		job1.m_pResultsMat = job1.m_pCandidate->m_pResult1;
		m_jobVector.push_back(job1);

		Job job2;
		job2.m_pCandidate = m_population[nPopIndex];
		job2.m_nJobType = MPI_JOB_2_TAG;
		job2.m_pResultsMat = job2.m_pCandidate->m_pResult2;
		m_jobVector.push_back(job2);
	}
	LOG1("CreateJobs finished, total number of jobs: %d", m_jobVector.size());
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::ProcessJobs()
{
	LOG("------------- ProcessJobs started -------------");
	LOG2("Total of %d jobs to send to %d slave processes", m_jobVector.size(), m_nNumberOfSlaveMachines);

	int nCurProcess = 0;

	while(true)
	{
		if(m_jobVector.size() == 0)
		{
			break;
		}
		Job job = m_jobVector.back();
		m_jobVector.pop_back();

		//LOG2("ProcessJobs, processing job - Candidate: #%d, JobType: %d", pJob->m_pCandidate->m_nIndex, pJob->m_nJobType);
		//LOG1("ProcessJobs, Sending job to process #%d", nCurProcess);
		MPI_Send(&(job.m_pCandidate->m_pFibroblastMat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, nCurProcess+1, job.m_nJobType, MPI::COMM_WORLD);
		//LOG1("ProcessJobs, Async wait for job from process #%d", nCurProcess);
		MPI_Irecv(&(job.m_pResultsMat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, nCurProcess+1, MPI_RESULT_TAG, MPI::COMM_WORLD, &m_mpiRequestArray[nCurProcess]);

		nCurProcess++;
		if(nCurProcess == (m_nNumberOfSlaveMachines))
		{
			//LOG("ProcessJobs, Wait for all pending jobs to finish");
			MPI_Waitall(nCurProcess, m_mpiRequestArray, m_mpiStatusArray);
			nCurProcess = 0;
			//LOG("ProcessJobs, Wait is over, continue processing jobs");
		}
	}

	if(nCurProcess != 0)
	{
		MPI_Waitall(nCurProcess, m_mpiRequestArray, m_mpiStatusArray);
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
		//LOG5("ProcessResults1 at %d, nCandidateResult: %d, nTargetResult: %d, diff: %d, nCost: %d", iW, nCandidateResult, nTargetResult, nCandidateResult - nTargetResult, nCost);
	}
	//LOG1("ProcessResults1, total cost for protocol1: %d", nCost1);

	// Calculate cost from result 2.
	for (int iH = 1; iH < Nh+1; iH += SAMPLING_INTERVALS)
	{
		double dCandidateResult = pCandidate->m_pResult2[iH][Nw-MeasurementMarginIndexes] - pCandidate->m_pResult2[iH][MeasurementMarginIndexes + 1];
		int nCandidateResult = (int)ceil(dCandidateResult * 1000);
		int nTargetResult = (int)ceil(m_pTargetMeasurement2[iH] * 1000);
		int nCost = std::abs(nCandidateResult - nTargetResult);
		pCandidate->m_cost += (unsigned long int)nCost;
		nCost2 += (unsigned long int)nCost;
		//LOG5("ProcessResults2 at %d, nCandidateResult: %d, nTargetResult: %d, diff: %d, nCost: %d", iH, nCandidateResult, nTargetResult, nCandidateResult - nTargetResult, nCost);
	}
	//LOG1("ProcessResults2, total cost for protocol2: %d", nCost2);

	LOG3("- ProcessResults, current cost for candidate #%d, %s: %u", pCandidate->m_nIndex, pCandidate->GetFullName(), pCandidate->m_cost);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::ProcessResults()
{
	LOG("ProcessResults started");
	for (int nPopIndex = 0; nPopIndex < Npop; ++nPopIndex)
	{
		Candidate* pCandidate = m_population[nPopIndex];
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
		Candidate* pCandidate = CreateRandomCandidate(iPop);
		m_population.push_back(pCandidate);
		AddChildToPastCandidates(pCandidate);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double CGA::CalculatePatchDistance(const FibroblastPatch& newFibroblastPatch)
{
	// Calculate the center of the patch from (0,0):
	double dHDistance = (double)newFibroblastPatch.m_nHStart +
		(newFibroblastPatch.m_nHEnd - newFibroblastPatch.m_nHStart)/2;

	double dWDistance = (double)newFibroblastPatch.m_nWStart +
		(newFibroblastPatch.m_nWEnd - newFibroblastPatch.m_nWStart)/2;

	double dDistance = sqrt(pow(dHDistance,2) + pow(dWDistance,2));
	/*
	LOG5("CalculatePatchDistance (%d,%d)-(%d,%d) = %f",
		newFibroblastPatch.m_nHStart,
		newFibroblastPatch.m_nWStart,
		newFibroblastPatch.m_nHEnd,
		newFibroblastPatch.m_nWEnd,
		dDistance)
	*/
	return dDistance;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::AddSortedToPatchesVector(FibroblastPatchVector& vecFibroblastPatchVector, FibroblastPatch newFibroblastPatch)
{
	bool bAdded = false;
	FibroblastPatchVector vecNewFibroblastPatchVector;
	double dDistance = CalculatePatchDistance(newFibroblastPatch);

	for(int iPatch = 0; iPatch < vecFibroblastPatchVector.size(); ++iPatch)
	{
		const FibroblastPatch& curFibroblastPatch = vecFibroblastPatchVector[iPatch];
		double dCurDistance = CalculatePatchDistance(curFibroblastPatch);
		if(dDistance <= dCurDistance && !bAdded)
		{
			vecNewFibroblastPatchVector.push_back(newFibroblastPatch);
			bAdded = true;
		}
		vecNewFibroblastPatchVector.push_back(curFibroblastPatch);
	}
	if(!bAdded)
	{
		vecNewFibroblastPatchVector.push_back(newFibroblastPatch);
	}
	vecFibroblastPatchVector = vecNewFibroblastPatchVector;
}

Candidate* CGA::CreateRandomCandidate(int nIndex)
{
	FibroblastPatchVector vecFibroblastPatchVector;
	for(int iPatch = 0; iPatch < NUMBER_OF_FIBROBLAST_PATCHES; ++iPatch)
	{
		int nPartitionMinH = Min_h_Fibroblast;
		int nPartitionMinW = Min_w_Fibroblast;
		int nPartitionMaxH = Max_h_Fibroblast;
		int nPartitionMaxW = Max_w_Fibroblast;
		/*
		switch(iPatch)
		{
			case 0:
				nPartitionMaxH = nPartitionMinH + nHPartitionSize;
				nPartitionMaxW = nPartitionMinW + nWPartitionSize;
				break;
			case 1:
				nPartitionMaxH = nPartitionMinH + nHPartitionSize;
				nPartitionMinW = nPartitionMaxW - nWPartitionSize;
				break;
			case 2:
				nPartitionMinH = nPartitionMaxH - nHPartitionSize;
				nPartitionMaxW = nPartitionMinW + nWPartitionSize;
				break;
			case 3:
				nPartitionMinH = nPartitionMaxH - nHPartitionSize;
				nPartitionMinW = nPartitionMaxW - nWPartitionSize;
				break;
			default:
				break;
		}
		*/

		//LOG5("CreateRandomCandidate partition #%d sizes: (%d,%d) -> (%d,%d)", iPatch, nPartitionMinH, nPartitionMinW, nPartitionMaxH, nPartitionMaxW);
		int nHStart = rand()%(nPartitionMaxH - nPartitionMinH + 1) + nPartitionMinH;
		int nWStart = rand()%(nPartitionMaxW - nPartitionMinW + 1) + nPartitionMinW;
		int nHEnd = rand()%(nPartitionMaxH - nHStart + 1) + nHStart;
		int nWEnd = rand()%(nPartitionMaxW - nWStart + 1) + nWStart;
		AddSortedToPatchesVector(vecFibroblastPatchVector, FibroblastPatch(nHStart, nWStart, nHEnd, nWEnd));
	}
	return new Candidate(nIndex, vecFibroblastPatchVector); //nHStart, nWStart, nHEnd, nWEnd);
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

void CGA::ClearPopulation()
{
	for(int iPop = 0; iPop < Npop; ++iPop)
	{
		delete m_population[iPop];
	}
	m_population.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::ClearPastCandidates()
{
	for(int iPastChild = 0; iPastChild < (int)PastCandidates.size(); iPastChild++)
	{
		char* pCandidateName = PastCandidates[iPastChild];
		delete [] pCandidateName;
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

void CGA::CreateNextGeneration()
{
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
		delete m_population[nNewCandidateIndex];
		LOG4("Creating offspring %d (candidate #%d) from parents %d and %d", iOffspring+1, nNewCandidateIndex, nFirstParent, nSecondParent);
		m_population[nNewCandidateIndex] = CreateChild(m_population[nFirstParent], m_population[nSecondParent], nNewCandidateIndex);
		if(m_population[nNewCandidateIndex] != NULL) // Only if we got a new child we can continue.
		{
			++iOffspring;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate* CGA::CreateChild(Candidate* pParent1, Candidate* pParent2, int nIndex)
{
	LOG2("Parent1, Candidate #%d: %s", pParent1->m_nIndex, pParent1->GetFullName());
	LOG2("Parent2, Candidate #%d: %s", pParent2->m_nIndex, pParent2->GetFullName());
	FibroblastPatchVector vecFibroblastPatchVector;
	for(int iPatch = 0; iPatch < NUMBER_OF_FIBROBLAST_PATCHES; ++iPatch)
	{
		int nHStart = Mutate((rand()%2 == 1) ? pParent1->GetFibroblastPatch(iPatch).m_nHStart : pParent2->GetFibroblastPatch(iPatch).m_nHStart);
		int nWStart = Mutate((rand()%2 == 1) ? pParent1->GetFibroblastPatch(iPatch).m_nWStart : pParent2->GetFibroblastPatch(iPatch).m_nWStart);
		int nHEnd = Mutate((rand()%2 == 1) ? pParent1->GetFibroblastPatch(iPatch).m_nHEnd : pParent2->GetFibroblastPatch(iPatch).m_nHEnd);
		int nWEnd = Mutate((rand()%2 == 1) ? pParent1->GetFibroblastPatch(iPatch).m_nWEnd : pParent2->GetFibroblastPatch(iPatch).m_nWEnd);
		AddSortedToPatchesVector(vecFibroblastPatchVector, FibroblastPatch(nHStart, nWStart, nHEnd, nWEnd));
	}
	Candidate* pChild = new Candidate(nIndex, vecFibroblastPatchVector);
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
				//printf("Target=%d at iH=%d, iW=%d\n", nTarget, iH, iW);
				//printf("Match=%d at iH=%d, iW=%d\n", nMatch, iH, iW);
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
					//printf("nMatch=%d at iH=%d, iW=%d\n", nMatch, iH, iW);
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

void CGA::RunGA(int iAlgoIndex, double** pCombinedFibroblastMat)
{
	LOG("Init GA...");
	if(!ReadTargetMeasurements())
	{
		LOG("Failed to read the target measurements. Aborting.");
		return;
	}

	LOG("***************************************");
	LOG("***************************************");
	LOG1("Starting algorithm run %d", iAlgoIndex);
	LOG("***************************************");
	LOG("***************************************");

	// Create the starting population of candidates.
	CreateRandomPopulation();

	// Create min cost file.
	char minCostFileName[FILE_NAME_BUFFER_SIZE] = {0};
	sprintf(minCostFileName, "%s/MinCost_%d.txt", LOG_FOLDER, iAlgoIndex);
	FILE* pMinCostFile = fopen(minCostFileName, "w");
	if(pMinCostFile != NULL)
	{
		fprintf(pMinCostFile, "Iteration | MinCost | Mat | Error | Coverage\n");
		fclose(pMinCostFile);
	}

	char matlabFileName[FILE_NAME_BUFFER_SIZE] = {0};
	sprintf(matlabFileName, "%s/Matlab_%d.txt", LOG_FOLDER, iAlgoIndex);
	FILE* pMatlabFile = fopen(matlabFileName, "w");
	if(pMatlabFile != NULL)
	{
		fclose(pMatlabFile);
	}

	// Start the loop.
	LOG("------------------");
	int nLastMinCostCounter = 0;
	int nCurIteration = 0;
	while(nCurIteration <= MaxIterations)
	{
		LOG2("Starting iteration #%d of algo run #%d", nCurIteration, iAlgoIndex);
		LOG("------------------");
		clock_t iterationStartingTime = clock();

		LOG("Candidates for this iteration:");
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			LOG2("Candidate #%d, %s", m_population[iPop]->m_nIndex, m_population[iPop]->GetFullName());
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
		std::sort(m_population.begin(), m_population.end(), CandidateCompare);
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			m_population[iPop]->m_nIndex = iPop;
			LOG3("ProcessResults, current cost for candidate #%d, %s: %d", m_population[iPop]->m_nIndex, m_population[iPop]->GetFullName(), m_population[iPop]->m_cost);
		}

		// Write output files (for presentation and debugging).
		unsigned long int nTotalCost = 0;
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			if(iPop == 0) // To save space of txt files...
			{
				char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
				sprintf(riseTimeCandidateFileName, "BestRiseTimes/RiseTime_Itr%d_Indx%d_Prot1.txt", nCurIteration, iPop);
				SaveMatToFile(m_population[iPop]->m_pResult1, riseTimeCandidateFileName);
				sprintf(riseTimeCandidateFileName, "BestRiseTimes/RiseTime_Itr%d_Indx%d_Prot2.txt", nCurIteration, iPop);
				SaveMatToFile(m_population[iPop]->m_pResult2, riseTimeCandidateFileName);
			}
			else
			{
				char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
				sprintf(riseTimeCandidateFileName, "OtherRiseTimes/RiseTime_Itr%d_Indx%d_Prot1.txt", nCurIteration, iPop);
				SaveMatToFile(m_population[iPop]->m_pResult1, riseTimeCandidateFileName);
				sprintf(riseTimeCandidateFileName, "OtherRiseTimes/RiseTime_Itr%d_Indx%d_Prot2.txt", nCurIteration, iPop);
				SaveMatToFile(m_population[iPop]->m_pResult2, riseTimeCandidateFileName);
			}
			nTotalCost += m_population[iPop]->m_cost;
		}
		double dAvgCost = nTotalCost/Npop;

		// Update the min cost file.
		MinCost[nCurIteration] = m_population[0]->m_cost;
		LOG1("Avg Cost: %.3f", dAvgCost);
		LOG2("Min Cost: %d for %s", MinCost[nCurIteration], m_population[0]->GetFullName());
		int nError = CalculateError(m_population[0]->m_pFibroblastMat);
		LOG1("The (geomatric) error for the best match is: %d", nError);
		double dCoverage = CalculateTargetCoverage(m_population[0]->m_pFibroblastMat);
		LOG1("The target coverage for the best match is: %.3f%%", dCoverage);
		pMinCostFile = fopen(minCostFileName, "a");
		if(pMinCostFile != NULL)
		{
			fprintf(pMinCostFile, "%d | %d | %s | %d | %.3f\n", nCurIteration, MinCost[nCurIteration], m_population[0]->GetFullName(), nError, dCoverage);
			fclose(pMinCostFile);
		}

		pMatlabFile = fopen(matlabFileName, "a");
		if(pMatlabFile != NULL)
		{
			fprintf(pMatlabFile, "\t%% %d | %d | %s | %d | %.3f\n", nCurIteration, MinCost[nCurIteration], m_population[0]->GetFullName(), nError, dCoverage);
			for(int iPatch = 0; iPatch < NUMBER_OF_FIBROBLAST_PATCHES; ++iPatch)
			{
				fprintf(pMatlabFile, "\tif(iPatch == %d)\n", iPatch);
				fprintf(pMatlabFile, "\t\tfound_startIndexes = [%d,%d];\n", m_population[0]->GetFibroblastPatch(iPatch).m_nHStart, m_population[0]->GetFibroblastPatch(iPatch).m_nWStart);
				fprintf(pMatlabFile, "\t\tfound_endIndexes = [%d,%d];\n", m_population[0]->GetFibroblastPatch(iPatch).m_nHEnd, m_population[0]->GetFibroblastPatch(iPatch).m_nWEnd);
				fprintf(pMatlabFile, "\tend\n");
			}
			fprintf(pMatlabFile, "\n%%----------------------------\n\n");
			fclose(pMatlabFile);
		}

		// Check break conditions.
		if(MinCost[nCurIteration] == 0)
		{
			// This won't really happen...
			LOG("Reached zero cost. Breaking the iterations !");
			break;
		}
		if(nCurIteration != 0)
		{
			if(MinCost[nCurIteration] == MinCost[nCurIteration-1])
			{
				nLastMinCostCounter++;
				if(nLastMinCostCounter == MAX_REPEATING_COSTS_FOR_DEAD_END)
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

		// Create the next generation.
		CreateNextGeneration();

		clock_t iterationEndingTime = clock();
		double iterationRunningTime = (iterationEndingTime - iterationStartingTime)/double(CLOCKS_PER_SEC);
		LOG1("Iteration start: %d", iterationStartingTime);
		LOG1("Iteration end: %d", iterationEndingTime);
		LOG1("Iteration duration: %.3f seconds", iterationRunningTime);

		LOG("------------------");
	}

	// Finished executing the GA.
	LOG2("The best match found is: %s, with cost: %d", m_population[0]->GetFullName(), MinCost[nCurIteration]);
	int nError = CalculateError(m_population[0]->m_pFibroblastMat);
	LOG1("The error for this match is: %d", nError);

	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{
			pCombinedFibroblastMat[iH][iW] += m_population[0]->m_pFibroblastMat[iH][iW];
		}
	}

	ClearPopulation();
	ClearPastCandidates();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void CGA::Test()
{
	//for(double dDiff = -0.0000020; dDiff <= 0.0000020; dDiff += 0.0000001)
	//for(double dj = -0.01; dj <= 0.01; dj += 0.0005)
	{

	LOG("----------------------- Test started -----------------------");

	int nIndex = 0;
	int nHStart[NUMBER_OF_FIBROBLAST_PATCHES] = {0};
	int nWStart[NUMBER_OF_FIBROBLAST_PATCHES] = {0};
	int nHEnd[NUMBER_OF_FIBROBLAST_PATCHES] = {0};
	int nWEnd[NUMBER_OF_FIBROBLAST_PATCHES] = {0};

	nHStart[0] = 38;
	nWStart[0] = 38;
	nHEnd[0] = 52;
	nWEnd[0] = 52;

	nHStart[1] = 88;
	nWStart[1] = 88;
	nHEnd[1] = 102;
	nWEnd[1] = 102;

	//Candidate* pCandidate = CreateRandomCandidate(nIndex);
	FibroblastPatchVector vecFibroblastPatch;
	Candidate* pCandidate = new Candidate(nIndex, vecFibroblastPatch); //nHStart, nWStart, nHEnd, nWEnd);

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
	pModel->ExecuteModel(pCandidate->m_pFibroblastMat, pCandidate->m_pResult1, s1, "/a/home/cc/students/enginer/shaishif/Output/ModelLogs");
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
