/******************************************************************
 * muensterTPCsim
 * 
 * Simulations of the Muenster TPC
 * 
 * @author Lutz Althüser
 * @date   2015-04-14
 *
 * @update 2015-11-02 - added comments
 *
 * @comment
 ******************************************************************/
#ifndef __muensterTPCPEVENTDATA_H__
#define __muensterTPCPEVENTDATA_H__

#include <string>
#include <vector>

using std::string;
using std::vector;

class muensterTPCEventData {
public:
	 muensterTPCEventData();
	~muensterTPCEventData();

public:
	void Clear();

public:
	int m_iEventId;								// the event ID
	int m_iNbTopPmtHits;					// number of top pmt hits
	int m_iNbBottomPmtHits;				// number of bottom pmt hits
	int m_iNbTopVetoPmtHits;			// number of top veto pmt hits
	int m_iNbBottomVetoPmtHits;		// number of bottom veto pmt hits
	vector<int> *m_pPmtHits;			// number of photon hits per pmt
	float m_fTotalEnergyDeposited;// total energy deposited in the ScintSD
	int m_iNbSteps;								// number of energy depositing steps
	vector<int> *m_pTrackId;			// id of the particle
	vector<int> *m_pParentId;			// id of the parent particle
	vector<int> *m_pParticleTypeCode;
	vector<string> *m_pParticleType;			// type of particle
	vector<string> *m_pParentType;				// type of particle
	vector<string> *m_pCreatorProcess;		// interaction
	vector<string> *m_pDepositingProcess;	// energy depositing process
	vector<float> *m_pX;					// position of the step
	vector<float> *m_pY;
	vector<float> *m_pZ;
	vector<float> *m_pEnergyDeposited; 			// energy deposited in the step
	vector<float> *m_pTrackInitEnergy;
	vector<float> *m_pKineticEnergy;	// particle kinetic energy after the step			
	vector<float> *m_pTime;						// time of the step
	vector<string> *m_pPrimaryParticleType;		// type of particle
	float m_fPrimaryEnergy;						// energy of the primary particle
	float m_fPrimaryX;								// position of the primary particle
	float m_fPrimaryY;
	float m_fPrimaryZ;	
};

#endif // __muensterTPCPEVENTDATA_H__

