
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

class Candidate
{
public:
	Candidate(int nIndex, FibroblastPatchVector vecFibroblastPatch);
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

private:
	#define CANDIDATE_MAX_NAME_LENGTH 		1024
	char m_cCandidateFullName[CANDIDATE_MAX_NAME_LENGTH];
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __CANDIDATE__

