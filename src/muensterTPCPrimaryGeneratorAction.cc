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
#include <globals.hh>
#include <G4RunManagerKernel.hh>
#include <G4Event.hh>
#include <Randomize.hh>

#include "muensterTPCParticleSource.hh"
#include "muensterTPCPrimaryGeneratorAction.hh"
#include "muensterTPCPrimaryGeneratorMessenger.hh"

muensterTPCPrimaryGeneratorAction::muensterTPCPrimaryGeneratorAction()
{
	m_pMessenger = new muensterTPCPrimaryGeneratorMessenger(this);
	m_pParticleSource = new muensterTPCParticleSource();

	m_hParticleTypeOfPrimary = "";
	m_dEnergyOfPrimary = 0.;
	m_hPositionOfPrimary = G4ThreeVector(0., 0., 0.);

	m_lSeeds[0] = -1;
	m_lSeeds[1] = -1;
}

muensterTPCPrimaryGeneratorAction::~muensterTPCPrimaryGeneratorAction()
{
	delete m_pParticleSource;
    
    //  write the source settings to file
    if (m_pParticleSource->GetEventInputFile() == "fromDecay0File") {
        TNamed *TypePar = new TNamed("SourceType", "Decay0File");
        TypePar->Write();
    } else {
        TNamed *TypePar =
        new TNamed("SourceType",
                   m_pParticleSource->GetParticleDefinition()->GetParticleName());
        TypePar->Write();
        
    }
}

void muensterTPCPrimaryGeneratorAction::GeneratePrimariesDecay0(G4Event *pEvent) {
    // generate a single particle
    m_lSeeds[0] = *(CLHEP::HepRandom::getTheSeeds());
    m_lSeeds[1] = *(CLHEP::HepRandom::getTheSeeds()+1);
    
    m_pParticleSource->GeneratePrimaryVertex(pEvent);
    /*
     // particle name of primary
     m_hParticleTypeOfPrimary = particleTable->FindParticle(particle->GetG4code())->GetParticleName();
     // kinetic energy of primary
     m_dEnergyOfPrimary       = particle->GetKineticEnergy();
     // position of primary
     m_hPositionOfPrimary     = vertex->GetPosition();
     // direction of primary
     m_hDirectionOfPrimary    = particle->GetMomentum();
     
     FillHistograms();
     
     numberOfAcceptedPrimaries += 1.0;
     */
}
    
void muensterTPCPrimaryGeneratorAction::GeneratePrimariesStandard(G4Event *pEvent)
{
	m_lSeeds[0] = *(CLHEP::HepRandom::getTheSeeds());
	m_lSeeds[1] = *(CLHEP::HepRandom::getTheSeeds()+1);

	G4StackManager *pStackManager = (G4RunManagerKernel::GetRunManagerKernel())->GetStackManager();

//    G4cout << "PrimaryGeneratorAction: track status: "
//        << pStackManager->GetNUrgentTrack() << " urgent, "
//        << pStackManager->GetNWaitingTrack() << " waiting, "
//        << pStackManager->GetNPostponedTrack() << " postponed"
//        << G4endl;

	if(!pStackManager->GetNPostponedTrack())
	{
		m_pParticleSource->GeneratePrimaryVertex(pEvent);
	}
	else
	{
		pStackManager->TransferStackedTracks(fPostpone, fUrgent);
		G4VTrajectory* pTrajectory;
		G4Track *pTrack = pStackManager->PopNextTrack(&pTrajectory);

		m_pParticleSource->GeneratePrimaryVertexFromTrack(pTrack, pEvent);

		delete pTrack;
	}
	G4PrimaryVertex *pVertex = pEvent->GetPrimaryVertex();
	G4PrimaryParticle *pPrimaryParticle = pVertex->GetPrimary();

	m_hParticleTypeOfPrimary = pPrimaryParticle->GetG4code()->GetParticleName();

	G4double dP = pPrimaryParticle->GetMomentum().mag();
	G4double dMass = pPrimaryParticle->GetMass();

	m_dEnergyOfPrimary = std::sqrt(dP*dP + dMass*dMass) - dMass;
	m_hPositionOfPrimary = pVertex->GetPosition();
}


void muensterTPCPrimaryGeneratorAction::GeneratePrimaries(G4Event *pEvent) {
    if (m_pParticleSource->GetEventInputFile() == "fromDecay0File") {
        GeneratePrimariesDecay0(pEvent);
    } else {
        GeneratePrimariesStandard(pEvent);
    }
}
