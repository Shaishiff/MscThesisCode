
#ifndef __CANDIDATE__
#define __CANDIDATE__

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class Candidate
{
public:
	Candidate(int nIndex, int nHStart, int nWStart, int nHEnd, int nWEnd);
	virtual ~Candidate();	
	/*
	int GetHStart() const {return m_nHStart;}
	int GetWStart() const {return m_nWStart;}
	int GetHEnd() const {return m_nHEnd;}
	int GetWEnd() const {return m_nWEnd;}
	*/
	
	//void Mutate();
	//void UnMutate();	
	//double** GetFibroblastMat() {return m_pFibroblastMat;}
	
	char* GetFullName();

private:
	void ClearMat();
	void CreateFibroblastMat();
	//void Mutate(int nParam, int nVal);

public:	
	double** m_pFibroblastMat;
	double** m_pResult1;
	double** m_pResult2;		
	int m_nIndex;
	unsigned long int m_cost;
	int m_nHStart;
	int m_nWStart;
	int m_nHEnd;
	int m_nWEnd;
	
	// Mutation.
	/*
	int m_nParam1;
	int m_nParam2;
	int m_nVal1;
	int m_nVal2;
	*/
private:
	char m_cCandidateFullName[256];
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate);

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __CANDIDATE__

