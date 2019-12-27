/******************************************************************
 * muensterTPCsim
 * 
 * Simulations of the Muenster TPC
 * 
 * @author Lutz Althüser, based on muensterTPC (Levy) and Xenon100
 * @date   2015-04-14
 *
 * @comment 
 ******************************************************************/
#include "muensterTPCEventData.hh"

muensterTPCEventData::muensterTPCEventData()
{
	m_iEventId = 0;
	m_iNbTopPmtHits = 0;
	m_iNbBottomPmtHits = 0;
	m_pPmtHits = new vector<int>;

	m_fTotalEnergyDeposited = 0.;
	m_iNbSteps = 0;

	m_pTrackId = new vector<int>;
	m_pParentId = new vector<int>;
	m_pParticleTypeCode = new vector<int>;
	m_pParticleType = new vector<string>;
	m_pParentType = new vector<string>;
	m_pCreatorProcess = new vector<string>;
	m_pDepositingProcess = new vector<string>;
	m_pX = new vector<float>;
	m_pY = new vector<float>;
	m_pZ = new vector<float>;
	m_pEnergyDeposited = new vector<float>;
	m_pTrackInitEnergy = new vector<float>;
	m_pKineticEnergy = new vector<float>;
	m_pTime = new vector<float>;

	m_pPrimaryParticleType = new vector<string>;
	m_fPrimaryEnergy = 0.;
	m_fPrimaryX = 0.;
	m_fPrimaryY = 0.;
	m_fPrimaryZ = 0.;	
}

muensterTPCEventData::~muensterTPCEventData()
{
	delete m_pPmtHits;
	delete m_pTrackId;
	delete m_pParentId;
	delete m_pParticleTypeCode;
	delete m_pParticleType;
	delete m_pParentType;
	delete m_pCreatorProcess;
	delete m_pDepositingProcess;
	delete m_pX;
	delete m_pY;
	delete m_pZ;
	delete m_pEnergyDeposited;
	delete m_pKineticEnergy;
	delete m_pTime;
	delete m_pTrackInitEnergy;
	delete m_pPrimaryParticleType;
}

void
muensterTPCEventData::Clear()
{
	m_iEventId = 0;
	m_iNbTopPmtHits = 0;
	m_iNbBottomPmtHits = 0;

	m_pPmtHits->clear();

	m_fTotalEnergyDeposited = 0.0;
	m_iNbSteps = 0;

	m_pTrackId->clear();
	m_pParentId->clear();
	m_pParticleTypeCode->clear();
	m_pParticleType->clear();
	m_pParentType->clear();
	m_pCreatorProcess->clear();
	m_pDepositingProcess->clear();
	m_pX->clear();
	m_pY->clear();
	m_pZ->clear();
	m_pEnergyDeposited->clear();
	m_pTrackInitEnergy->clear();
	m_pKineticEnergy->clear();
	m_pTime->clear();

	m_pPrimaryParticleType->clear();
	m_fPrimaryEnergy = 0.;
	m_fPrimaryX = 0.;
	m_fPrimaryY = 0.;
	m_fPrimaryZ = 0.;	
}

