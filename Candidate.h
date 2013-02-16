
#ifndef __CANDIDATE__
#define __CANDIDATE__

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#define CANDIDATE_MAX_NAME_LENGTH 		1024
#define NUMBER_OF_FIBROBLAST_PATCHES	1

class Candidate
{
public:
	Candidate(int nIndex, 
		int nHStart[NUMBER_OF_FIBROBLAST_PATCHES], 
		int nWStart[NUMBER_OF_FIBROBLAST_PATCHES], 
		int nHEnd[NUMBER_OF_FIBROBLAST_PATCHES], 
		int nWEnd[NUMBER_OF_FIBROBLAST_PATCHES]);
	virtual ~Candidate();	
	char* GetFullName();

private:
	void CreateFibroblastBorders();
	void CreateFibroblastPatch();
	
public:	
	double** m_pFibroblastMat;
	double** m_pResult1;
	double** m_pResult2;		
	int m_nIndex;
	unsigned long int m_cost;
	int m_nHStart[NUMBER_OF_FIBROBLAST_PATCHES];
	int m_nWStart[NUMBER_OF_FIBROBLAST_PATCHES];
	int m_nHEnd[NUMBER_OF_FIBROBLAST_PATCHES];
	int m_nWEnd[NUMBER_OF_FIBROBLAST_PATCHES];
	
private:
	char m_cCandidateFullName[CANDIDATE_MAX_NAME_LENGTH];
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __CANDIDATE__

