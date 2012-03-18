
#include "defs.h"
#include "Candidate.h"
#include <algorithm>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate::Candidate(int nIndex)
{
	m_pResult1 = new double*[Nh+2];
	m_pResult2 = new double*[Nh+2];
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		m_pResult1[iH] = new double[Nw+2];
		m_pResult2[iH] = new double[Nw+2];
		for(int iW = 0; iW < Nw+2; ++iW)
		{											
			m_pResult1[iH][iW] = 0;
			m_pResult2[iH][iW] = 0;
		}
	}

	m_nCenterH = 0;
	m_nCenterW = 0;
	m_nHeight = 0;
	m_nWidth = 0;
	m_pModel = new CFkModel();
	m_pModel2 = new CFkModel();
	m_nIndex = nIndex;
	m_cost = 0.0;
	m_mat = new double*[Nh+2];
	ClearMat();
	if(nIndex == -1)
	{
		m_nCenterH = TargetCenterH;
		m_nCenterW = TargetCenterW;
		m_nHeight = TargetHeight;
		m_nWidth = TargetWidth;
		CreateFibroblastPatch();
	}
	else
	{
		AddFibroblastPatch(1);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

Candidate::~Candidate()
{
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		delete [] m_mat[iH];
		delete [] m_pResult1[iH];
		delete [] m_pResult2[iH];
	}
	delete [] m_mat;
	delete [] m_pResult1;
	delete [] m_pResult2;
	delete m_pModel;
	delete m_pModel2;
}
	
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::ClearMat()
{
	for(int iH = 0; iH < Nh+2; ++iH)
	{
		m_mat[iH] = new double[Nw+2];		
		for(int iW = 0; iW < Nw+2; ++iW)
		{								
			if((iH == 0) || (iH == Nh+1) || (iW == 0) || (iW == Nw+1))
			{
				m_mat[iH][iW] = 1;
			}
			else
			{
				m_mat[iH][iW] = 0;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::CreateFibroblastPatch()
{
	ClearMat();
	int nhStart = max((m_nCenterH),1);
	int nhEnd = min((m_nCenterH+m_nHeight),Nh);
	int nwStart = max((m_nCenterW),1);
	int nwEnd = min((m_nCenterW+m_nWidth),Nw);
	for (int iH=nhStart; iH < nhEnd; ++iH)
	{	
		for (int iW=nwStart; iW < nwEnd; ++iW)
		{
			m_mat[iH][iW] = 1;
		}
	}
	sprintf(m_str,"(%d,%d) %dX%d", m_nCenterH, m_nCenterW, m_nHeight, m_nWidth);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::Mutate()
{
	m_nParam1 = rand()%4;
	m_nVal1 = rand()%9 - 4;
	m_nParam2 = rand()%4;
	while(m_nParam2 == m_nParam1) { m_nParam2 = rand()%4; }
	m_nVal2 = rand()%9 - 4;

	switch(m_nParam1)
	{
		case 0:
			m_nCenterH += m_nVal1;
			break;
		case 1:
			m_nCenterW += m_nVal1;
			break;
		case 2:
			m_nHeight += m_nVal1;
			break;
		case 3:
			m_nWidth += m_nVal1;
			break;
	}		
	switch(m_nParam2)
	{
		case 0:
			m_nCenterH += m_nVal2;
			break;
		case 1:
			m_nCenterW += m_nVal2;
			break;
		case 2:
			m_nHeight += m_nVal2;
			break;
		case 3:
			m_nWidth += m_nVal2;
			break;
	}

	m_nCenterH = min(max(m_nCenterH, 1), Nh);
	m_nCenterW = min(max(m_nCenterW, 1), Nw);
	m_nHeight = min(max(m_nHeight, 1), Nh);
	m_nWidth = min(max(m_nWidth, 1), Nw);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::UnMutate()
{
	switch(m_nParam1)
	{
		case 0:
			m_nCenterH -= m_nVal1;
			break;
		case 1:
			m_nCenterW -= m_nVal1;
			break;
		case 2:
			m_nHeight -= m_nVal1;
			break;
		case 3:
			m_nWidth -= m_nVal1;
			break;
	}		
	switch(m_nParam2)
	{
		case 0:
			m_nCenterH -= m_nVal2;
			break;
		case 1:
			m_nCenterW -= m_nVal2;
			break;
		case 2:
			m_nHeight -= m_nVal2;
			break;
		case 3:
			m_nWidth -= m_nVal2;
			break;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Candidate::AddFibroblastPatch(int nFibroblast)
{
	m_nCenterH = rand()%int(Nh) + 1;
	m_nCenterW = rand()%int(Nw) + 1;
	m_nHeight = rand()%int(Nh) + 1;
	m_nWidth = rand()%int(Nw) + 1;
	CreateFibroblastPatch();
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

bool CandidateCompare(Candidate* pFirstCandidate, Candidate* pSecondCandidate)
{
	return (pFirstCandidate->m_cost < pSecondCandidate->m_cost);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

