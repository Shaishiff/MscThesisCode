
#ifndef __CANDIDATE__
#define __CANDIDATE__

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

struct FibroblastPatch
{
	int m_nHStart;
	int m_nWStart;
	int m_nHEnd;
	int m_nWEnd;

	FibroblastPatch()
	{
		m_nHStart = 0;
		m_nWStart = 0;
		m_nHEnd = 0;
		m_nWEnd = 0;
	}

	FibroblastPatch(int nHStart,
		int nWStart,
		int nHEnd,
		int nWEnd)
	{
		m_nHStart = nHStart;
		m_nWStart = nWStart;
		m_nHEnd = nHEnd;
		m_nWEnd = nWEnd;
	}
};

typedef std::vector<FibroblastPatch> FibroblastPatchVector;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#define CANDIDATE_MAX_NAME_LENGTH 		1024
#define NUMBER_OF_FIBROBLAST_PATCHES	4

class Candidate
{
public:
	Candidate(int nIndex,
		FibroblastPatchVector vecFibroblastPatch);
		//int nHStart[NUMBER_OF_FIBROBLAST_PATCHES], 
		//int nWStart[NUMBER_OF_FIBROBLAST_PATCHES], 
		//int nHEnd[NUMBER_OF_FIBROBLAST_PATCHES], 
		//int nWEnd[NUMBER_OF_FIBROBLAST_PATCHES]);
	virtual ~Candidate();	
	char* GetFullName();
	FibroblastPatch GetFibroblastPatch(int iPatch) const;

private:
	void CreateFibroblastBorders();
	void CreateFibroblastPatch();
	
public:	
	double** m_pFibroblastMat;
	double** m_pResult1;
	double** m_pResult2;		
	int m_nIndex;
	unsigned long int m_cost;

private:
	FibroblastPatchVector m_vecFibroblastPatch;

	//int m_nHStart[NUMBER_OF_FIBROBLAST_PATCHES];
	//int m_nWStart[NUMBER_OF_FIBROBLAST_PATCHES];
	//int m_nHEnd[NUMBER_OF_FIBROBLAST_PATCHES];
	//int m_nWEnd[NUMBER_OF_FIBROBLAST_PATCHES];
	
private:
	char m_cCandidateFullName[CANDIDATE_MAX_NAME_LENGTH];
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __CANDIDATE__

