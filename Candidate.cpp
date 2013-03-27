
#include "defs.h"
#include "Candidate.h"
#include "Mat.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate::Candidate(int nIndex, 
	int nHStart[NUMBER_OF_FIBROBLAST_PATCHES], 
	int nWStart[NUMBER_OF_FIBROBLAST_PATCHES], 
	int nHEnd[NUMBER_OF_FIBROBLAST_PATCHES], 
	int nWEnd[NUMBER_OF_FIBROBLAST_PATCHES])
{
	// Create the mats.
	m_pFibroblastMat = CreateMat();	
	m_pResult1 = CreateMat();
	m_pResult2 = CreateMat();
		
	// Init the vars.
	m_nIndex = nIndex;
	m_cost = 0;
	memset(m_cCandidateFullName, 0, sizeof(m_cCandidateFullName));
	int nTemp = 0;
	for(int i = 0; i < NUMBER_OF_FIBROBLAST_PATCHES; ++i)
	{
		m_nHStart[i] = nHStart[i];
		m_nWStart[i] = nWStart[i];
		m_nHEnd[i] = nHEnd[i];
		m_nWEnd[i] = nWEnd[i];
		/*if(m_nHStart[i] > m_nHEnd[i])
		{
			nTemp = m_nHStart[i];
			m_nHStart[i] = m_nHEnd[i];
			m_nHEnd[i] = nTemp;
		}
		if(m_nWStart[i] > m_nWEnd[i])
		{
			nTemp = m_nWStart[i];
			m_nWStart[i] = m_nWEnd[i];
			m_nWEnd[i] = nTemp;
		}*/
	}	
	
	// Now create the fibroblats.
	CreateFibroblastBorders();
	CreateFibroblastPatch();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate::~Candidate()
{
	DestroyMat(m_pFibroblastMat);
	DestroyMat(m_pResult1);
	DestroyMat(m_pResult2);
}
	
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::CreateFibroblastBorders()
{
	for (int iH = 0; iH < Nh_with_border; ++iH)
	{
		for (int iW = 0; iW < Nw_with_border; ++iW)
		{								
			if((iH == 0) || (iH == Nh+1) || (iW == 0) || (iW == Nw+1))
			{
				m_pFibroblastMat[iH][iW] = 1.0;
			}
			else
			{
				m_pFibroblastMat[iH][iW] = 0.0;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::CreateFibroblastPatch()
{
	for(int iPatch = 0; iPatch < NUMBER_OF_FIBROBLAST_PATCHES; ++iPatch)
	{
		for (int iH = m_nHStart[iPatch]; iH <= m_nHEnd[iPatch]; ++iH)
		{	
			for (int iW = m_nWStart[iPatch]; iW <= m_nWEnd[iPatch]; ++iW)
			{
				m_pFibroblastMat[iH][iW] = 1.0;
			}
		}	
	}	
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

char* Candidate::GetFullName()
{ 
	char temp[CANDIDATE_MAX_NAME_LENGTH] = {0};
	for(int iPatch = 0; iPatch < NUMBER_OF_FIBROBLAST_PATCHES; ++iPatch)
	{
		sprintf(temp,"%s", m_cCandidateFullName);
		if(iPatch == 0)
		{
			sprintf(m_cCandidateFullName,"(%d,%d) -> (%d,%d)", m_nHStart[iPatch], m_nWStart[iPatch], m_nHEnd[iPatch], m_nWEnd[iPatch]);
		}
		else
		{
			sprintf(m_cCandidateFullName,"%s + (%d,%d) -> (%d,%d)", temp, m_nHStart[iPatch], m_nWStart[iPatch], m_nHEnd[iPatch], m_nWEnd[iPatch]);
		}		
	}	
	return m_cCandidateFullName; 
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate)
{
	return (pFirstCandidate->m_cost < pSecondCandidate->m_cost);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

