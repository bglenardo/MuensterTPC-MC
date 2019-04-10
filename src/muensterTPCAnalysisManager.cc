/******************************************************************
 * muensterTPCsim
 * 
 * Simulations of the Muenster TPC
 * 
 * @author 	Lutz Althüser, based on Xenon100 simulation
 * @date   	2015-04-14
 * 
 * @update 	2015-10-01
 *					added descriptions and bug fixes
 * @update	2015-07-07 
 *					new TTree write and AutoSave ability
 ******************************************************************/

// include GEANT4 classes
#include <G4SDManager.hh>
#include <G4Run.hh>
#include <G4Event.hh>
#include <G4HCofThisEvent.hh>
#include <G4SystemOfUnits.hh>
#include <G4Version.hh>

// include C++ classes
#include <numeric>

// include ROOT classes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TDirectory.h>

// include Muenster TPC classes
#include "muensterTPCPrimaryGeneratorAction.hh"
#include "muensterTPCAnalysisManager.hh"
#include "muensterTPCEventData.hh"
#include "muensterTPCLXeHit.hh"
#include "muensterTPCPmtHit.hh"
#include "muensterTPCDetectorConstruction.hh"

//******************************************************************/
// creation and initialization of the AnalysisManager
//******************************************************************/
muensterTPCAnalysisManager::muensterTPCAnalysisManager(muensterTPCPrimaryGeneratorAction *pPrimaryGeneratorAction) {
  G4int threadID = G4Threading::IsMultithreadedApplication() ? G4Threading::G4GetThreadId() : 0;
	
	// initialization of the HitsCollectionID variables 
	m_iLXeHitsCollectionID = -1;
	m_iPmtHitsCollectionID = -1;

	// default output file name (which should be redifined in the main class)
	m_hDataFilename = "events.root";

	// reference to the initialized PrimaryGeneratorAction
	m_pPrimaryGeneratorAction = pPrimaryGeneratorAction;
	// declaration of the EventData class
	m_pEventData = new muensterTPCEventData();
	writeEmptyEvents = kFALSE;
}


//******************************************************************/
//
//******************************************************************/
muensterTPCAnalysisManager::~muensterTPCAnalysisManager(){
}

//******************************************************************/
// define the BeginOfRun actions / prepare the output file for data output
//******************************************************************/
void muensterTPCAnalysisManager::BeginOfRun(const G4Run *pRun) {
		// do we write empty events or not?
		writeEmptyEvents = m_pPrimaryGeneratorAction->GetWriteEmpty();
  
		// create output file
		m_pTreeFile = new TFile(m_hDataFilename.c_str(), "RECREATE", "File containing event data for muensterTPCsim");
		TNamed *G4version = new TNamed("G4VERSION_TAG",G4VERSION_TAG);
		G4version->Write();
		TNamed *G4MCname = new TNamed("MC_TAG","muensterTPC");
		G4MCname->Write();
		TNamed *G4MCVersion = new TNamed("MCVERSION_TAG","X.Y.Z");
		G4MCname->Write();
		
		_events = m_pTreeFile->mkdir("events");
		_events->cd();
		
		// create ROOT Tree for the simulation data
		m_pTree = new TTree("events", "Tree containing event data for muensterTPCsim");

		// include missing ROOT classes
		gROOT->ProcessLine("#include <vector>");

		// initialize all tree branches for the different data types

		//******************************************************************/	
		// branches for each main event/particle beam:
		//******************************************************************/
		// eventid:	The eventid contains the individual ID of each event. Since 'event' is used in different situations,
		//					it means, every particle created by 'beamOn' starts a detector event. All following events will be
		//					be saved within the same eventid in specific branches (see below).
		//					Acces to the eventid in ROOT: int eventid;
		//																				T1->SetBranchAddress("eventid", &eventid);
		m_pTree->Branch("eventid", &m_pEventData->m_iEventId, "eventid/I");
		// ntpmthits:	total amount of top PMT hits for a specific eventid/particle beam
		//						Acces in ROOT: 	int ntpmthits;
		//														T1->SetBranchAddress("ntpmthits", &ntpmthits);
		m_pTree->Branch("ntpmthits", &m_pEventData->m_iNbTopPmtHits, "ntpmthits/I");
		// nbpmthits:	total amount of bottom PMT hits for a specific eventid/particle beam
		//						Acces in ROOT: 	int nbpmthits;
		//														T1->SetBranchAddress("nbpmthits", &nbpmthits);
		m_pTree->Branch("nbpmthits", &m_pEventData->m_iNbBottomPmtHits, "nbpmthits/I");

		//m_pTree->Branch("ntvetopmthits", &m_pEventData->m_iNbTopVetoPmtHits, "ntvetopmthits/I");
		//m_pTree->Branch("nbvetopmthits", &m_pEventData->m_iNbBottomVetoPmtHits, "nbvetopmthits/I");

		// pmthits:	total amount of PMT hits for a specific eventid and for each PMT
		//						Acces in ROOT: 	vector<int> *pmthits= new vector<int>;
		//														T1->SetBranchAddress("pmthits", &pmthits);
		//						Note: Do not access pmthits without calling a specifig vector element.
		m_pTree->Branch("pmthits", "vector<int>", &m_pEventData->m_pPmtHits);
		// etot:	Amount of energy, which is deopsited during this eventid/particle run.
		//				Acces in ROOT: 	float etot;
		//												T1->SetBranchAddress("etot", &etot);
		m_pTree->Branch("etot", &m_pEventData->m_fTotalEnergyDeposited, "etot/F");
		// nbpmthits:	total amount of bottom PMT hits for a specific eventid/particle beam
		//						Acces in ROOT: 	int nbpmthits;
		//														T1->SetBranchAddress("nbpmthits", &nbpmthits);
		m_pTree->Branch("nsteps", &m_pEventData->m_iNbSteps, "nsteps/I");
	
		//******************************************************************/	
		// branches for each event/particle which is created by the main event
		//******************************************************************/
		// trackid:	ID of the current event in the event track. All listed events are
		//					generated within the main eventid. (e.g. emitted gammas)
		//					Acces in ROOT: 	vector<int> *trackid= new vector<int>;
		//													T1->SetBranchAddress("trackid", &trackid);
		m_pTree->Branch("trackid", "vector<int>", &m_pEventData->m_pTrackId);
		// type:	type of the particles in the event track
		//				Acces in ROOT: 	vector<string> *type= new vector<string>;
		//												T1->SetBranchAddress("type", &type);
		m_pTree->Branch("type", "vector<string>", &m_pEventData->m_pParticleType);
		// parentid:	trackid of the parent track event
		//						Acces in ROOT: 	vector<int> *parentid= new vector<int>;
		//														T1->SetBranchAddress("parentid", &parentid);
		m_pTree->Branch("parentid", "vector<int>", &m_pEventData->m_pParentId);
		// parenttype:	parenttype of the parent track event
		//							Acces in ROOT: 	vector<string> *parenttype= new vector<string>;
		//															T1->SetBranchAddress("parenttype", &parenttype);
		m_pTree->Branch("parenttype", "vector<string>", &m_pEventData->m_pParentType);
		// creaproc:	name of the creation process of the track particle/trackid
		//						Acces in ROOT: 	vector<string> *creaproc= new vector<string>;
		//														T1->SetBranchAddress("creaproc", &creaproc);
		m_pTree->Branch("creaproc", "vector<string>", &m_pEventData->m_pCreatorProcess);
		// edproc:	name of the energy deposition process of the track particle/trackid
		//					Acces in ROOT: 	vector<string> *edproc= new vector<string>;
		//													T1->SetBranchAddress("edproc", &edproc);
		m_pTree->Branch("edproc", "vector<string>", &m_pEventData->m_pDepositingProcess);
		// Positions of the current particle/trackid
		// 		Acces in ROOT: 		vector<float> *xp= new vector<float>;
		//											T1->SetBranchAddress("xp", &xp);
		m_pTree->Branch("xp", "vector<float>", &m_pEventData->m_pX);
		// 		Acces in ROOT: 		vector<float> *yp= new vector<float>;
		//											T1->SetBranchAddress("yp", &yp);
		m_pTree->Branch("yp", "vector<float>", &m_pEventData->m_pY);
		// 		Acces in ROOT: 		vector<float> *zp= new vector<float>;
		//											T1->SetBranchAddress("zp", &zp);
		m_pTree->Branch("zp", "vector<float>", &m_pEventData->m_pZ);
		// ed:	energy deposition of the current particle/trackid
		// 			Acces in ROOT: 		vector<float> *ed= new vector<float>;
		//												T1->SetBranchAddress("ed", &ed);
		m_pTree->Branch("ed", "vector<float>", &m_pEventData->m_pEnergyDeposited);
		// time:	timestamp of the current particle/trackid
		// 				Acces in ROOT: 		vector<float> *time= new vector<float>;
		//													T1->SetBranchAddress("time", &time);
		m_pTree->Branch("time", "vector<float>", &m_pEventData->m_pTime);

		//******************************************************************/	
		// branches for each event/particle which contain information about the primary particle
		//******************************************************************/
		// type_pri:	type of the primary event/main event
		//						Acces in ROOT: 	vector<string> *type_pri= new vector<string>;
		//														T1->SetBranchAddress("type_pri", &type_pri);
		m_pTree->Branch("type_pri", "vector<string>", &m_pEventData->m_pPrimaryParticleType);
		// Energy and positions of the current particle/trackid
		// 		Acces in ROOT:	vector<float> *e_pri= new vector<float>;
		//										T1->SetBranchAddress("e_pri", &e_pri);
		m_pTree->Branch("e_pri", &m_pEventData->m_fPrimaryEnergy, "e_pri/F");
		// 		Acces in ROOT:	vector<float> *xp_pri= new vector<float>;
		//										T1->SetBranchAddress("xp_pri", &xp_pri);
		m_pTree->Branch("xp_pri", &m_pEventData->m_fPrimaryX, "xp_pri/F");
		// 		Acces in ROOT:	vector<float> *yp_pri= new vector<float>;
		//										T1->SetBranchAddress("yp_pri", &yp_pri);	
		m_pTree->Branch("yp_pri", &m_pEventData->m_fPrimaryY, "yp_pri/F");
		// 		Acces in ROOT:	vector<float> *zp_pri= new vector<float>;
		//										T1->SetBranchAddress("zp_pri", &zp_pri);
		m_pTree->Branch("zp_pri", &m_pEventData->m_fPrimaryZ, "zp_pri/F");

		//m_pTree->SetMaxTreeSize(10e9); /previous
		m_pTree->SetMaxTreeSize(1000*Long64_t(2000000000)); //2TB
		m_pTree->AutoSave();
	
		// Write the number of events in the output file
		m_iNbEventsToSimulate = pRun->GetNumberOfEventToBeProcessed();
		m_pNbEventsToSimulateParameter = new TParameter<int>("nbevents", m_iNbEventsToSimulate);
		m_pNbEventsToSimulateParameter->Write();
}

//******************************************************************/
// EndOfRun action/end of the simulation
//******************************************************************/
void muensterTPCAnalysisManager::EndOfRun(const G4Run *pRun) {
		// write and remove old revisions
		m_pTreeFile->Write(0,TObject::kOverwrite);
		//m_pTreeFile->Write();
		m_pTreeFile->Close();
}

//******************************************************************/
//	BeginOfEvent action - for each beamed particle
//******************************************************************/
void muensterTPCAnalysisManager::BeginOfEvent(const G4Event *pEvent) {
	// initialize the HitCollections if needed (-1 = not initialized)
	//G4cout << "Begin of Event" << G4endl;
	if(m_iLXeHitsCollectionID == -1)
	{
		G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
		m_iLXeHitsCollectionID = pSDManager->GetCollectionID("LXeHitsCollection");
	} 

//	if(m_iPmtHitsCollectionID == -1)
//	{
//		G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
//		m_iPmtHitsCollectionID = pSDManager->GetCollectionID("PmtHitsCollection");
//	}
}

//******************************************************************/
// EndOfEvent action - getting all event data
//******************************************************************/
void muensterTPCAnalysisManager::EndOfEvent(const G4Event *pEvent) {
	G4HCofThisEvent* pHCofThisEvent = pEvent->GetHCofThisEvent();
	muensterTPCLXeHitsCollection* pLXeHitsCollection = 0;
	muensterTPCPmtHitsCollection* pPmtHitsCollection = 0;

	G4int iNbLXeHits = 0, iNbPmtHits = 0;
	
	if(pHCofThisEvent)
	{
		if(m_iLXeHitsCollectionID != -1)
		{
			pLXeHitsCollection = (muensterTPCLXeHitsCollection *)(pHCofThisEvent->GetHC(m_iLXeHitsCollectionID));
			iNbLXeHits = (pLXeHitsCollection)?(pLXeHitsCollection->entries()):(0);
		}

		if(m_iPmtHitsCollectionID != -1)
		{
			pPmtHitsCollection = (muensterTPCPmtHitsCollection *)(pHCofThisEvent->GetHC(m_iPmtHitsCollectionID));
			iNbPmtHits = (pPmtHitsCollection)?(pPmtHitsCollection->entries()):(0);
		}
	}

	m_pEventData->m_iEventId = pEvent->GetEventID();

	m_pEventData->m_pPrimaryParticleType->push_back(m_pPrimaryGeneratorAction->GetParticleTypeOfPrimary());

	m_pEventData->m_fPrimaryEnergy = m_pPrimaryGeneratorAction->GetEnergyOfPrimary()/keV;
	m_pEventData->m_fPrimaryX = m_pPrimaryGeneratorAction->GetPositionOfPrimary().x()/mm;
	m_pEventData->m_fPrimaryY = m_pPrimaryGeneratorAction->GetPositionOfPrimary().y()/mm;
	m_pEventData->m_fPrimaryZ = m_pPrimaryGeneratorAction->GetPositionOfPrimary().z()/mm;

	G4int iNbSteps = 0;
	G4float fTotalEnergyDeposited = 0.;
	
	if(iNbLXeHits || iNbPmtHits)
	{
		// LXe hits
		for(G4int i=0; i<iNbLXeHits; i++)
		{
			muensterTPCLXeHit *pHit = (*pLXeHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pTrackId->push_back(pHit->GetTrackId());
				m_pEventData->m_pParentId->push_back(pHit->GetParentId());

				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());
				m_pEventData->m_pParentType->push_back(pHit->GetParentType());
				m_pEventData->m_pCreatorProcess->push_back(pHit->GetCreatorProcess());
				m_pEventData->m_pDepositingProcess->push_back(pHit->GetDepositingProcess());

				m_pEventData->m_pX->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pY->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZ->push_back(pHit->GetPosition().z()/mm);

				fTotalEnergyDeposited += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDeposited->push_back(pHit->GetEnergyDeposited()/keV);

				m_pEventData->m_pKineticEnergy->push_back(pHit->GetKineticEnergy()/keV);
				m_pEventData->m_pTime->push_back(pHit->GetTime()/second);

				iNbSteps++;
			}
		}

		m_pEventData->m_iNbSteps = iNbSteps;
		m_pEventData->m_fTotalEnergyDeposited = fTotalEnergyDeposited;

		G4int iNbTopPmts = (G4int) muensterTPCDetectorConstruction::GetGeometryParameter("NbTopPmts");
		G4int iNbBottomPmts = (G4int) muensterTPCDetectorConstruction::GetGeometryParameter("NbBottomPmts");
		G4int iNbTopVetoPmts = (G4int) muensterTPCDetectorConstruction::GetGeometryParameter("NbTopVetoPmts");
		G4int iNbBottomVetoPmts = (G4int) muensterTPCDetectorConstruction::GetGeometryParameter("NbBottomVetoPmts");

		m_pEventData->m_pPmtHits->resize(iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts+iNbBottomVetoPmts, 0);

		// Pmt hits
		for(G4int i=0; i<iNbPmtHits; i++)
			(*(m_pEventData->m_pPmtHits))[(*pPmtHitsCollection)[i]->GetPmtNb()]++;

		m_pEventData->m_iNbTopPmtHits =	accumulate(m_pEventData->m_pPmtHits->begin(), m_pEventData->m_pPmtHits->begin()+iNbTopPmts, 0);
		m_pEventData->m_iNbBottomPmtHits = accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts, m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts, 0);
		// m_pEventData->m_iNbTopVetoPmtHits = accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts, m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts, 0);
		// m_pEventData->m_iNbBottomVetoPmtHits =	accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts, m_pEventData->m_pPmtHits->end(), 0);

		//if((fTotalEnergyDeposited > 0. || iNbPmtHits > 0) && !FilterEvent(m_pEventData))
		
	    // save only energy depositing events
	    if(writeEmptyEvents) {
			m_pTree->Fill(); // write all events to the tree
	    } else {
		    if(fTotalEnergyDeposited > 0. || iNbPmtHits > 0) m_pTree->Fill(); // only events with some activity are written to the tree
	    }

		// auto save functionality to avoid data loss/ROOT can recover aborted simulations
		if ( pEvent->GetEventID() % 10000 == 0)
			m_pTree->AutoSave();

		m_pEventData->Clear();
	}
}

//******************************************************************/
//
//******************************************************************/
void muensterTPCAnalysisManager::Step(const G4Step *pStep) {
}

//******************************************************************/
// Event filer function which can be called above
//******************************************************************/
G4bool muensterTPCAnalysisManager::FilterEvent(muensterTPCEventData *pEventData) {
	G4double dEnergyDepositedSensitiveRegion = 0.;

	vector<float> *pX = pEventData->m_pX;
	vector<float> *pY = pEventData->m_pY;
	vector<float> *pZ = pEventData->m_pZ;
	vector<float> *pEnergyDeposited = pEventData->m_pEnergyDeposited;

	const G4double dDriftLength = muensterTPCDetectorConstruction::GetGeometryParameter("DriftLength");
	const G4double dRadius = muensterTPCDetectorConstruction::GetGeometryParameter("TeflonCylinderInnerRadius");

	for(G4int i=0; i<pEnergyDeposited->size(); i++)
	{
		if((*pZ)[i] < 0. && (*pZ)[i] > -dDriftLength && std::sqrt((*pX)[i]*(*pX)[i] + (*pY)[i]*(*pY)[i]) < dRadius)
			dEnergyDepositedSensitiveRegion += (*pEnergyDeposited)[i];
	}

//if(dEnergyDepositedSensitiveRegion > 0. && dEnergyDepositedSensitiveRegion < 100.)
	if(dEnergyDepositedSensitiveRegion > 0.)
		return false;
	else
		return true;
}

	
