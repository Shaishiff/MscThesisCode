
#ifndef __SAFE_VEC__
#define __SAFE_VEC__	

#include <pthread.h>
#include "defs.h"
#include "Candidate.h"

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

struct Job
{
	Candidate* m_pCandidate;
	int m_nJobType;
	double** m_pResultsMat;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

class CSafeJobVector
{
public:
	CSafeJobVector()
	{		
	}
	int GetSize()
	{
		return (int)m_vec.size();
	}
	bool IsEmpty()
	{
		pthread_mutex_lock(&m_mutex);
		bool bIsEmpty = m_vec.empty();		
		pthread_mutex_unlock(&m_mutex);
		return bIsEmpty;
	}	
	void AddJob(Job* pJob)
	{
		pthread_mutex_lock(&m_mutex);
		m_vec.push_back(pJob);
		pthread_mutex_unlock(&m_mutex);
	}
	Job* GetJob()
	{
		pthread_mutex_lock(&m_mutex);
		if(m_vec.empty())
		{
			pthread_mutex_unlock(&m_mutex);
			return NULL;
		}
		Job* pJob = m_vec.back();
		m_vec.pop_back();
		pthread_mutex_unlock(&m_mutex);
		return pJob;
	}

private:
	vector<Job*> m_vec;
	static pthread_mutex_t m_mutex;
};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

#endif // __SAFE_VEC__
