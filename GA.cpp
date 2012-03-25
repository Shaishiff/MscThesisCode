
#include <pthread.h>
#include "GA.h"
#include "FkModel.h"

#define LOG_FOLDER "/a/home/cc/students/enginer/shaishif/Logs"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

double** CreateMat()
{
	double* data = (double *)malloc(Nh_with_border*Nw_with_border*sizeof(double));
    double** mat = (double **)malloc(Nh_with_border*sizeof(double*));
    for (int i = 0; i < Nh_with_border; i++)
	{
        mat[i] = &(data[Nw_with_border*i]);
	}
	
	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{
			mat[iH][iW] = 0.0;
		}
	}
	
	return mat;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void DestroyMat(double** mat)
{
	free(mat[0]);
	free(mat);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void PrintMat(double** mat)
{
	printf("Printing mat\n"); fflush(stdout);
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			printf("%d,%d: %f\n", iH, iW, mat[iH][iW]); 
			fflush(stdout);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool SaveMatToFile(double** mat, char* fileName)
{
	char fullFileName[1024] = {0};
	sprintf(fullFileName, "%s/%s", LOG_FOLDER, fileName);
	FILE* pFile = fopen(fullFileName, "w");
	if(pFile == NULL)
	{
		//printf("ERROR ------------ Failed to open file: %s ------------ ERROR\n", fullFileName);
		return false;
	}
	
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		for (int iW = 1; iW < Nw+1; ++iW)
		{
			fprintf(pFile, "%4.0f ", mat[iH][iW]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	//printf("Saved file: %s\n", fullFileName);
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

FILE* pLogFile;
char strLogSourceName[1024];
char strLogFileName[1024];
char strLog[1024];
#define LOG printf("%s%s", strLogSourceName, strLog); fflush(stdout); \
	pLogFile = fopen(strLogFileName, "a"); \
	if(pLogFile != NULL){ fprintf(pLogFile,"%s%s", strLogSourceName, strLog); fclose(pLogFile); }

Ga ga;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int Ga::nNumberOfMachines = 0;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void* ThreadFunction(void* arg)
{
	int* pThreadIndex = (int*)arg;
	ga.JobProcessingThreadFunc(*pThreadIndex);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Ga::Ga()
{
	m_nCurIteration = 0;
	m_pTargetMeasurement1 = new double[Nw];
	m_pTargetMeasurement2 = new double[Nh];
	MinCost = NULL;
	Rank = NULL;
	for(int iThread = 0; iThread < MAX_NUMBER_OF_THREADS; iThread++)
	{
		m_threadIndexArray[iThread] = iThread;
		m_threadFlagArray[iThread] = true;
		m_threadHandleArray[iThread] = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::InitGa()
{		
	for(int iPop = 0; iPop < Npop; ++iPop)
	{		
		Candidate* pCandidate = new Candidate(iPop);
		Population.push_back(pCandidate);

		char fileName[FILE_NAME_BUFFER_SIZE] = {0};
		sprintf(fileName, "OriginalCandidate_%d.txt", iPop);
		SaveMatToFile(pCandidate->m_pFibroblastMat, fileName);
	}

	MinCost = new double[MaxIterations];
	for(int iteration=0; iteration < MaxIterations; ++iteration)
	{
		MinCost[iteration] = 0.0;
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
	
	/*
	// Create the threads which will allocate the jobs to the
	// different machines. We start from 1 (and not 0) because
	// 0 is the index of the master machine (which runs the main program).
	sprintf(strLog, "Creating %d threads to communicate with the slave processes\n", nNumberOfMachines-1);
	LOG
    for(int iThread = 1; iThread < nNumberOfMachines; iThread++)
	{		
		int nRet = pthread_create(&m_threadHandleArray[iThread], NULL, ThreadFunction, (void*)(&m_threadIndexArray[iThread]));
	} 
	*/	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CreateTargetMeasurements()
{
	sprintf(strLog, "CreateTargetMeasurements\n");
	LOG	
	MPI_Status status;
	Candidate* pCandidate = new Candidate(-1);			
	
	sprintf(strLog, "Sending fibroblasts\n");
	LOG	
	MPI_Send(&(pCandidate->m_pFibroblastMat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, 1, MPI_JOB_1_TAG, MPI::COMM_WORLD);
	SaveMatToFile(pCandidate->m_pFibroblastMat, "TargetFibroblastMat.txt");
	
	// We wait for the result.
	sprintf(strLog, "Waiting to receive result mat 1 | Machine: %d\n", 1); 
	LOG
	MPI_Recv(&(pCandidate->m_pResult1[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, 1, MPI_RESULT_TAG, MPI::COMM_WORLD, &status);
	sprintf(strLog, "Received result mat 1 | Machine: %d\n", 1);
	LOG
	SaveMatToFile(pCandidate->m_pResult1, "TargetFibroblastMatResults1.txt");

	sprintf(strLog, "Sending fibroblasts");
	LOG	
	MPI_Send(&(pCandidate->m_pFibroblastMat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, 1, MPI_JOB_2_TAG, MPI::COMM_WORLD);
	
	sprintf(strLog, "Waiting to receive result mat 2 | Machine: %d\n", 1);
	LOG
	MPI_Recv(&(pCandidate->m_pResult2[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, 1, MPI_RESULT_TAG, MPI::COMM_WORLD, &status);
	sprintf(strLog, "Received result mat 2 | Machine: %d\n", 1);
	LOG
	SaveMatToFile(pCandidate->m_pResult2, "TargetFibroblastMatResults2.txt");
	
	for (int iW = 1; iW < Nw+1; ++iW)
	{
		m_pTargetMeasurement1[iW-1] = pCandidate->m_pResult1[Nh][iW];
	}
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		m_pTargetMeasurement2[iH-1] = pCandidate->m_pResult2[iH][Nw];
	}
	
	/*
	S1Protocol s1;
	S2Protocol s2;
	CFkModel* pModel = new CFkModel();	
	Candidate* pTargetCandidate = new Candidate(-1);		
	
	pModel->ExecuteModel(pTargetCandidate->GetFibroblastMat(), pTargetCandidate->m_pResult1, s1);	
	for (int iW = 1; iW < Nw+1; ++iW)
	{
		m_pTargetMeasurement1[iW-1] = pTargetCandidate->m_pResult1[Nh][iW];
	}
	
	pModel->ExecuteModel(pTargetCandidate->GetFibroblastMat(), pTargetCandidate->m_pResult2, s2);
	for (int iH = 1; iH < Nh+1; ++iH)
	{
		m_pTargetMeasurement2[iH-1] = pTargetCandidate->m_pResult2[iH][Nw];
	}
	*/
	/*
	SaveMatToFile(pCandidate->m_pResult1, "TargetRiseTime1.txt");
	SaveMatToFile(pCandidate->m_pResult2, "TargetRiseTime2.txt");
	sprintf(strLog, "Target: %s\n", pCandidate->GetFullName());
	LOG
	*/
	/*
	delete pModel;
	pModel = NULL;
	*/
	delete pCandidate;
	pCandidate = NULL;	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::ProcessJobs(int nStartPopIndex, int nEndPopIndex)
{
	sprintf(strLog, "ProcessJobs from %d to %d\n", nStartPopIndex, nEndPopIndex);
	LOG
	for (int nState = 0; nState < 2; ++nState)
	{
		for (int nPopIndex = nStartPopIndex; nPopIndex < nEndPopIndex; ++nPopIndex)
		{			
			// We don't need to recalculate the first candidate because
			// it hasn't been mutated.
			bool bAddJob = (Population[nPopIndex]->m_nIndex != 0);
			
			// But we do need to calculate is cost if this is 
			// the first iteration ever.
			bAddJob |= (m_nCurIteration == 0);
			
			if(bAddJob)
			{	
				Candidate* pCandidate = Population[nPopIndex];			
				int nCurMachine = nPopIndex%(nNumberOfMachines-1) + 1;
				sprintf(strLog, "Job for machine %d, candidate #%d: %s\n", nCurMachine, Population[nPopIndex]->m_nIndex, Population[nPopIndex]->GetFullName());
				LOG
								
				switch(nState)
				{
					case 0:
					{
						//double dFlag = MPI_FLAG_START_JOB;
						printf("ProcessJobs | Sending start job flag | Machine: %d\n", nCurMachine); fflush(stdout);
						//MPI_Send(&dFlag, 1, MPI_DOUBLE, nCurMachine, MPI_FLAG_MSG_TAG, MPI::COMM_WORLD);
						sleep(1);

						// Send the actual matrix to process.						
						double** mat = pCandidate->m_pFibroblastMat;
						printf("ProcessJobs | Processing candidate %d | Machine: %d\n", pCandidate->m_nIndex, nCurMachine); fflush(stdout);
						printf("ProcessJobs | Sending mat for processing | Machine: %d\n", nCurMachine); fflush(stdout);
						//MPI_Send(mat, (Nw+2)*(Nh+2), MPI_DOUBLE, nCurMachine, MPI_JOB_MSG_TAG, MPI::COMM_WORLD);
						sleep(1);
						break;
					}
					case 1:
					{
						MPI_Status status;
		
						// We wait for the result.
						printf("ProcessJobs | Waiting to receive result mat 1 | Machine: %d\n", nCurMachine); fflush(stdout);
						//MPI_Recv(pCandidate->m_pResult1, (Nw+2)*(Nh+2), MPI_DOUBLE, nCurMachine, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD, &status);
						printf("ProcessJobs | Received result mat 1 | Machine: %d\n", nCurMachine);  fflush(stdout);
						
						printf("ProcessJobs | Waiting to receive result mat 2 | Machine: %d\n", nCurMachine);  fflush(stdout);
						//MPI_Recv(pCandidate->m_pResult2, (Nw+2)*(Nh+2), MPI_DOUBLE, nCurMachine, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD, &status);
						printf("ProcessJobs | Received result mat 2 | Machine: %d\n", nCurMachine);  fflush(stdout);
						
						/*	
						char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
						sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_1.txt", 0, pCandidate->m_nIndex);
						CFkModel::SaveToOutputFile(pCandidate->m_pResult1, riseTimeCandidateFileName);
						sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_2.txt", 0, pCandidate->m_nIndex);
						CFkModel::SaveToOutputFile(pCandidate->m_pResult2, riseTimeCandidateFileName);
						
						// Calculate the cost for this candidate.
						pCandidate->m_cost = 0;
						
						// Calculate cost from result 1.
						for (int iW = 1; iW < Nw+1; ++iW)
						{
							pCandidate->m_cost += std::abs((long int)(pCandidate->m_pResult1[Nh][iW] - m_pTargetMeasurement1[iW]));
						}

						// Calculate cost from result 2.
						for (int iH = 1; iH < Nh+1; ++iH)
						{
							pCandidate->m_cost += std::abs((long int)(pCandidate->m_pResult2[iH][Nw] - m_pTargetMeasurement2[iH]));
						}
						*/
						//sprintf("ProcessJobs | calculated cost: %.2f, candidate: %d, Machine: %d\n", pCandidate->m_cost, pCandidate->m_nIndex, nCurMachine);
						//delete status;
					}
				}			
			}
		}
	}
}

void Ga::JobProcessingThreadFunc(int nThreadIndex)
{	
/*
	printf("JobProcessingThreadFunc | Started | Thread: %d\n", nThreadIndex);				
	while(ga.GetThreadFlag(nThreadIndex))
	{
		Job* pJob = ga.GetJob();
		if(pJob != NULL)
		{	
			printf("JobProcessingThreadFunc | Read job | Thread: %d\n", nThreadIndex);
						
			// Tell the corresponding machine that we want to start a new job.
			double dFlag = MPI_FLAG_START_JOB;
			printf("JobProcessingThreadFunc | Sending start job flag | Thread: %d\n", nThreadIndex);
			MPI_Send(&dFlag, 1, MPI_DOUBLE, nThreadIndex, MPI_FLAG_MSG_TAG, MPI::COMM_WORLD);
			sleep(1);
			
			// Send the actual matrix to process.
			Candidate* pCandidate = pJob->m_pCandidate;
			double** mat = pCandidate->m_pFibroblastMat;
			printf("JobProcessingThreadFunc | Processing candidate %d | Thread: %d\n", pCandidate->m_nIndex, nThreadIndex);
			
			printf("JobProcessingThreadFunc | Sending mat for processing | Thread: %d\n", nThreadIndex);
			MPI_Send(mat, (Nw+2)*(Nh+2), MPI_DOUBLE, nThreadIndex, MPI_JOB_MSG_TAG, MPI::COMM_WORLD);
			sleep(1);
		
			MPI_Status* status = new MPI_Status();
		
			// We wait for the result.
			printf("JobProcessingThreadFunc | Waiting to receive result mat 1 | Thread: %d\n", nThreadIndex);
			MPI_Recv(pCandidate->m_pResult1, (Nw+2)*(Nh+2), MPI_DOUBLE, nThreadIndex, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD, status);
			printf("JobProcessingThreadFunc | Received result mat 1 | Thread: %d\n", nThreadIndex);
			sleep(1);
			
			printf("JobProcessingThreadFunc | Waiting to receive result mat 2 | Thread: %d\n", nThreadIndex);
			MPI_Recv(pCandidate->m_pResult2, (Nw+2)*(Nh+2), MPI_DOUBLE, nThreadIndex, MPI_RESULT_MSG_TAG, MPI::COMM_WORLD, status);
			printf("JobProcessingThreadFunc | Received result mat 2 | Thread: %d\n", nThreadIndex);
			sleep(1);
					
			char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_1.txt", 0, pCandidate->m_nIndex);
			SaveMatToFile(pCandidate->m_pResult1, riseTimeCandidateFileName);
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_2.txt", 0, pCandidate->m_nIndex);
			SaveMatToFile(pCandidate->m_pResult2, riseTimeCandidateFileName);
			
			// Calculate the cost for this candidate.
			pCandidate->m_cost = 0;
			
			// Calculate cost from result 1.
			for (int iW = 1; iW < Nw+1; ++iW)
			{
				pCandidate->m_cost += std::abs((long int)(pCandidate->m_pResult1[Nh][iW] - m_pTargetMeasurement1[iW]));
			}

			// Calculate cost from result 2.
			for (int iH = 1; iH < Nh+1; ++iH)
			{
				pCandidate->m_cost += std::abs((long int)(pCandidate->m_pResult2[iH][Nw] - m_pTargetMeasurement2[iH]));
			}
			
			printf("JobProcessingThreadFunc | calculated cost: %.2f, candidate: %d, thread: %d\n", pCandidate->m_cost, pCandidate->m_nIndex, nThreadIndex);
			
			delete status;
			status = NULL;
			
			// Delete this job and move on to the next one.
			delete pJob;
			pJob = NULL;
		}
		usleep(500);
	}
	
	// No more jobs to process.
	printf("JobProcessingThreadFunc | No more jobs, exiting | Thread: %d\n", nThreadIndex);
	*/
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::CalculateCosts()
{
	sprintf(strLog, "CalculatingCosts started\n"); 
	LOG
		
	// Create the jobs.
	for (int nPopIndex = 0; nPopIndex < Npop; ++nPopIndex)
	{			
		// We don't need to recalculate the first candidate because
		// it hasn't been mutated.
		bool bAddJob = (Population[nPopIndex]->m_nIndex != 0);
		
		// But we do need to calculate is cost if this is 
		// the first iteration ever.
		bAddJob |= (m_nCurIteration == 0);
		
		if(bAddJob)
		{	
			sprintf(strLog, "Creating a job for candidate #%d: %s\n", Population[nPopIndex]->m_nIndex, Population[nPopIndex]->GetFullName());
			LOG
		
			Job* pJob = new Job();
			pJob->m_pCandidate = Population[nPopIndex];
			m_jobVector.AddJob(pJob);
		}
	}
	
	while(!m_jobVector.IsEmpty())
	{
		usleep(100);
	}
	
	sprintf(strLog, "CalculatingCosts finished\n"); 
	LOG
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
		
	sprintf(strLog,"Parent1, Candidate #%d: %s\n", pParent1->m_nIndex, pParent1->GetFullName());
	LOG
	sprintf(strLog,"Parent2, Candidate #%d: %s\n", pParent2->m_nIndex, pParent2->GetFullName());
	LOG

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
	sprintf(strLog,"Child,   Candidate #%d: %s\n", pChild->m_nIndex, pChild->GetFullName());	
	LOG
	sprintf(strLog,"-\n");
	LOG
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
				sprintf(strLog,"Candiates are identical, #%d & #%d\n", nIndex, iPop);
				LOG
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
	sprintf(strLog,"Starting mutations...\n");
	LOG
	sprintf(strLog,"--\n");
	LOG
	for(int iPop = 1; iPop < Npop; ++iPop)
	{
		sprintf(strLog,"Candiate #%d before mutation: %s\n", iPop, Population[iPop]->GetFullName());
		LOG
		Population[iPop]->Mutate();
		while(FindSimilarCandidate(iPop))
		{
			Population[iPop]->UnMutate();
			Population[iPop]->Mutate();
		}
		Population[iPop]->CreateFibroblasts();
		Population[iPop]->m_cost = 0.0;
		sprintf(strLog,"Candiate #%d after mutation : %s\n", iPop, Population[iPop]->GetFullName());
		LOG
		sprintf(strLog,"-\n");
		LOG
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Ga::RunGa()
{
	sprintf(strLog, "Init GA...\n");
	LOG
	InitGa();

	// Create and measure our target mat.
	CreateTargetMeasurements();
	
	sprintf(strLog, "------------------\n");	
	LOG
	while(m_nCurIteration <= MaxIterations)
	{
		sprintf(strLog, "Starting iteration #%d\n", m_nCurIteration);
		LOG
		sprintf(strLog, "------------------\n");
		LOG
		clock_t iterationStartingTime = clock();
	
		// Determine the cost function for each member of the population.
		//CalculateCosts();
		int nStartPopIndex = 0;
		while(nStartPopIndex < Npop)
		{
			int nEndPopIndex = nStartPopIndex + (nNumberOfMachines-1);
			ProcessJobs(nStartPopIndex, nEndPopIndex);
			nStartPopIndex = nEndPopIndex;
		}

		// Sort according to the cost and then mate according to the rank.
		sprintf(strLog, "Sorting the population...\n");
		LOG
		std::sort(Population.begin(), Population.end(), CandidateCompare);
		for(int iPop = 0; iPop < Npop; ++iPop)
		{
			Population[iPop]->m_nIndex = iPop;
		}

		char bestCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
		sprintf(bestCandidateFileName, "BestCandidate_%d.txt", m_nCurIteration);
		SaveMatToFile(Population[0]->m_pFibroblastMat, bestCandidateFileName);

		double dAvgCost = 0.0;
		for(int iPop = 0; iPop < Npop; ++iPop)
		{			
			char fileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(fileName, "Candidate_%d_%d.txt", m_nCurIteration, iPop);
			SaveMatToFile(Population[iPop]->m_pFibroblastMat, fileName);

			char riseTimeCandidateFileName[FILE_NAME_BUFFER_SIZE] = {0};
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_1.txt", m_nCurIteration, iPop);
			SaveMatToFile(Population[iPop]->m_pResult1, riseTimeCandidateFileName);
			sprintf(riseTimeCandidateFileName, "CandidateRiseTime_%d_%d_2.txt", m_nCurIteration, iPop);
			SaveMatToFile(Population[iPop]->m_pResult2, riseTimeCandidateFileName);

			dAvgCost += Population[iPop]->m_cost;
		}
		dAvgCost = dAvgCost/Npop;

		// Create the next generation.		
		sprintf(strLog,"--\n");
		LOG
		sprintf(strLog,"Creating the next generation...\n");
		LOG
		sprintf(strLog,"--\n");
		LOG
		int iOffspring = 0;
		while(iOffspring < Nmates)
		{	
			if(iOffspring == 0)
			{
				Candidate* Parent1 = Population[0];
				Candidate* Parent2 = Population[0];
				CreateChild(Parent1, Parent2, Population[NsurvivingPopulation + iOffspring]);
				++iOffspring;
			}
			else
			{
				int nFirstParent = GetMate();
				int nSecondParent = GetMate();
				if(nFirstParent != nSecondParent)
				{
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
		++m_nCurIteration;

		MinCost[m_nCurIteration] = Population[0]->m_cost;
		sprintf(strLog,"Avg Cost: %.3f\n", dAvgCost);
		LOG
		sprintf(strLog,"Min Cost: %.3f for %s\n", MinCost[m_nCurIteration], Population[0]->GetFullName());
		LOG
		if(MinCost[m_nCurIteration] <= 0.0)
		{
			break;
		}

		clock_t iterationEndingTime = clock();
		double iterationRunningTime = (iterationEndingTime - iterationStartingTime)/double(CLOCKS_PER_SEC);
		sprintf(strLog,"Iteration duration: %.3f seconds\n", iterationRunningTime);
		LOG
		sprintf(strLog,"------------------\n");
		LOG
	}
		
    // Terminate the threads.
	for(int iThread = 1; iThread < nNumberOfMachines; iThread++)
	{		
		m_threadFlagArray[iThread] = false;
	}		
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartMainProcess(char* sMachineName)
{
	sprintf(strLog, "Starting the main process on %s\n", sMachineName);
	LOG
	
	// This is for the controlling master process.
	srand((unsigned int)time(NULL));
	
	// Start timing.
	clock_t mainStartingTime = clock();
		
	// Run main program.
	ga.RunGa();

	for(int iProcess = 1; iProcess < Ga::nNumberOfMachines; ++iProcess)
	{
		MPI_Send(NULL, 0, MPI_DOUBLE, iProcess, MPI_DIE_TAG, MPI::COMM_WORLD);
		sprintf(strLog, "Sending quit flag to process %d\n", iProcess);
		LOG
	}
	
	// End timing.
	clock_t mainEndingTime = clock();
	double mainRunningTime = (mainEndingTime - mainStartingTime)/double(CLOCKS_PER_SEC);
	sprintf(strLog, "Main duration: %.3f seconds\n", mainRunningTime);
	LOG
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void StartSlaveProcess(int nProcess, char* sMachineName)
{
	sprintf(strLog, "Process %i on %s | Starting\n", nProcess, sMachineName);
	LOG
		
	// This is for all the other processes which are not the master.		
	// Create all the vars we will use for data transfer.
	MPI_Status status; // can save resources by using the predefined constant MPI_STATUS_IGNORE as a special value for the status argument.	
	S1Protocol s1;
	S2Protocol s2;
	
	// Start the infinite loop (until the master tells us to quit).
	sprintf(strLog, "Process %i on %s | Starting loop\n", nProcess, sMachineName);
	LOG
	while(true)
	{			
		CFkModel* pModel = new CFkModel();	
		double** fibroblast_mat = CreateMat();
		double** result_mat = CreateMat();
		
		sprintf(strLog, "Process %i on %s | Waiting to receive mat\n", nProcess, sMachineName);
		LOG
		int nRet = MPI_Recv(&(fibroblast_mat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, MPI_MASTER, MPI_ANY_TAG, MPI::COMM_WORLD, &status);
		sprintf(strLog, "Process %i on %s | Mat received, res: %d\n", nProcess, sMachineName, nRet);
		LOG
		
		if (status.MPI_TAG == MPI_JOB_1_TAG) 
		{
			sprintf(strLog, "Process %i on %s | Executing 1st protocol.\n", nProcess, sMachineName);
			LOG
			pModel->ExecuteModel(fibroblast_mat, result_mat, s1);
		}
		else if (status.MPI_TAG == MPI_JOB_2_TAG) 
		{
			sprintf(strLog, "Process %i on %s | Executing 2nd protocol.\n", nProcess, sMachineName);
			LOG
			pModel->ExecuteModel(fibroblast_mat, result_mat, s2);
		}
		else if (status.MPI_TAG == MPI_DIE_TAG) 
		{
			sprintf(strLog, "Process %i on %s | Got the die tag. Aborting.\n", nProcess, sMachineName);
			LOG
			break;
		}
		else
		{
			sprintf(strLog, "Process %i on %s | Got an invalid tag. Aborting.\n", nProcess, sMachineName);
			LOG
			break;
		}
		
		sprintf(strLog, "Process %i on %s | Finished executing protocol, sending results.\n", nProcess, sMachineName);
		LOG		
		MPI_Send(&(result_mat[0][0]), Nh_with_border*Nw_with_border, MPI_DOUBLE, MPI_MASTER, MPI_RESULT_TAG, MPI::COMM_WORLD);
		sprintf(strLog, "Process %i on %s | Results were sent.\n", nProcess, sMachineName);
		LOG	
		
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
	/*
	if(nProvided == MPI_THREAD_SINGLE)
	{
		printf("MPI_THREAD_SINGLE\n");
	}
	else if(nProvided == MPI_THREAD_FUNNELED)
	{
		printf("MPI_THREAD_FUNNELED\n");
	}
	else if(nProvided == MPI_THREAD_SERIALIZED)
	{
		printf("MPI_THREAD_SERIALIZED\n");
	}
	else if(nProvided == MPI_THREAD_MULTIPLE)
	{
		printf("MPI_THREAD_MULTIPLE\n");
	}
	*/
	
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
		
	sprintf(strLog, "Starting process = %i on %s, out of %i processes\n", nCurProcess, sMachineName, Ga::nNumberOfMachines);
	LOG
	
	sprintf(strLog, "MPI Init ret: %d, provided thread support: %d for process = %i on %s, out of %i processes\n", nMpiIntRet, nProvided, nCurProcess, sMachineName, Ga::nNumberOfMachines);
	LOG
	
	if(nCurProcess == MPI_MASTER)
	{		
		StartMainProcess(sMachineName);		
	}
	else
	{		
		StartSlaveProcess(nCurProcess, sMachineName);
	}
	
	printf("Ending process = %i on %s, out of %i processes\n", nCurProcess, sMachineName, Ga::nNumberOfMachines);
	MPI_Finalize();
		
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// TODO:
// Create threads the same number as the number of slave processes.
// Each thread will take a job from the queue and process it.
// Processing involves: Executing model, saving results and measurements and updating the cost
// of the specific jobs candidate.
// The thread will exit when the queue is empty.
// On exit the thread will signal that it's done and this is how the master process
// will know he can continue on.
// Don't forget that when master process exists he needs to signal all other slave
// processes to exit as well by sending them the appropriate flag.
