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
#include <G4ios.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4Track.hh>
#include <G4Event.hh>
#include <G4VProcess.hh>
#include <G4StackManager.hh>

#include "muensterTPCAnalysisManager.hh"

#include "muensterTPCStackingAction.hh"

muensterTPCStackingAction::muensterTPCStackingAction(muensterTPCAnalysisManager *pAnalysisManager)
{
	m_pAnalysisManager = pAnalysisManager;
}

muensterTPCStackingAction::~muensterTPCStackingAction()
{
}

G4ClassificationOfNewTrack
muensterTPCStackingAction::ClassifyNewTrack(const G4Track *pTrack)
{
	G4ClassificationOfNewTrack hTrackClassification = fUrgent;

	if(pTrack->GetDefinition()->GetParticleType() == "nucleus" && !pTrack->GetDefinition()->GetPDGStable())
	{
		if(pTrack->GetParentID() > 0 && pTrack->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay")
			//hTrackClassification = fPostpone;
			hTrackClassification = fKill;
	}

	return hTrackClassification;
}

void
muensterTPCStackingAction::NewStage()
{
}

void
muensterTPCStackingAction::PrepareNewEvent()
{ 
}

