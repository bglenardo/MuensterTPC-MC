/******************************************************************
 * muensterTPCsim
 * 
 * Simulations of the Muenster TPC
 * 
 * @author Lutz Althüser, based on muensterTPC (Levy) and Xenon100
 * @date   2015-04-14
 * @update 2016-02-16
 *
 * @comment - ready for G4 V10.02
 *					- added some festures of the new G4 versions (modular physics lists)
 ******************************************************************/
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleWithCuts.hh>
#include <G4ParticleTypes.hh>
#include <G4ParticleTable.hh>
#include <G4HadronCaptureProcess.hh>
#include <G4UserLimits.hh>
#include "G4UserSpecialCuts.hh"
#include <G4ios.hh>
#include <globals.hh>

#include "G4SystemOfUnits.hh"
#include "G4PhysicsListHelper.hh"

#include <iomanip>

#include <G4EmStandardPhysics.hh>
#include <G4EmLivermorePhysics.hh>
#include <G4EmPenelopePhysics.hh>

#include "muensterTPCPhysicsList.hh"
#include "muensterTPCPhysicsMessenger.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"

// Constructor /////////////////////////////////////////////////////////////
muensterTPCPhysicsList::muensterTPCPhysicsList():G4VUserPhysicsList() {

	defaultCutValue = 1.0 * mm;	//
	cutForGamma = defaultCutValue;
	cutForElectron = defaultCutValue;
	cutForPositron = defaultCutValue;

	VerboseLevel = 0;
	OpVerbLevel = 0;

	SetVerboseLevel(VerboseLevel);

  m_pMessenger = new muensterTPCPhysicsMessenger(this);
	
  particleList = new G4DecayPhysics("decays");
}

// Destructor //////////////////////////////////////////////////////////////
muensterTPCPhysicsList::~muensterTPCPhysicsList()
{  
  delete emPhysicsList;
  delete particleList;
  //	delete opPhysicsList;
  delete m_pMessenger;
  
  for(size_t i=0; i<hadronPhys.size(); i++) {
    delete hadronPhys[i];
  }
}

// Construct Particles /////////////////////////////////////////////////////
void muensterTPCPhysicsList::ConstructParticle() {
	// In this method, static member functions should be called
	// for all particles which you want to use.
	// This ensures that objects of these particle types will be
	// created in the program. 

  //if (m_hHadronicModel == "custom"){
	//G4cout << "----- ConstructMyBosons" << G4endl;
	ConstructMyBosons();
  //G4cout << "----- ConstructMyLeptons" << G4endl;
	ConstructMyLeptons();
  //G4cout << "----- ConstructMyHadrons" << G4endl;
	ConstructMyHadrons();
  //G4cout << "----- ConstructMyShortLiveds" << G4endl;
	ConstructMyShortLiveds();
  //} else if (m_hHadronicModel == "QGSP_BERT_HP") {
  particleList->ConstructParticle();
  //}
}

// construct Bosons://///////////////////////////////////////////////////
void muensterTPCPhysicsList::ConstructMyBosons() {
	// pseudo-particles
	G4Geantino::GeantinoDefinition();
	G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
	G4Gamma::GammaDefinition();

  //OpticalPhotons
	G4OpticalPhoton::OpticalPhotonDefinition();

}

// construct Leptons://///////////////////////////////////////////////////
void muensterTPCPhysicsList::ConstructMyLeptons() {
  // leptons
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	G4MuonPlus::MuonPlusDefinition();
	G4MuonMinus::MuonMinusDefinition();

	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

// construct Hadrons://///////////////////////////////////////////////////
void muensterTPCPhysicsList::ConstructMyHadrons()
{
 	//  mesons
	G4MesonConstructor mConstructor;
	mConstructor.ConstructParticle();

	// baryons
	G4BaryonConstructor bConstructor;
	bConstructor.ConstructParticle();

	// ions
	G4IonConstructor iConstructor;
	iConstructor.ConstructParticle();

}

// construct Shortliveds://///////////////////////////////////////////////////
// this can be added if necessary
//#include "G4ShortLivedConstructor.hh"

void muensterTPCPhysicsList::ConstructMyShortLiveds()
{  
	//G4ShortLivedConstructor slConstructor;
  	//slConstructor.ConstructParticle();
}

// Construct Processes //////////////////////////////////////////////////////
void muensterTPCPhysicsList::ConstructProcess()
{
  G4cout <<"MuensterTPCPhysicsList::MuensterTPCPhysicsList() EM physics: "<< m_hEMlowEnergyModel << G4endl;
  G4cout <<"MuensterTPCPhysicsList::MuensterTPCPhysicsList() Cerenkov: "<< m_bCerenkov << G4endl;

	AddTransportation();

	// EM physics
  if (m_hEMlowEnergyModel == "emstandard") {
    emPhysicsList = new G4EmStandardPhysics(); 
  } else if (m_hEMlowEnergyModel == "emlivermore"){
    emPhysicsList = new G4EmLivermorePhysics(); 
  } else if (m_hEMlowEnergyModel == "empenelope"){
    emPhysicsList = new G4EmPenelopePhysics(); 
  } else if (m_hEMlowEnergyModel == "old") {
    G4cout << "MuensterTPCPhysicsList::ConstructProcess() WARNING: Old version of low energy EM processes ... "<<G4endl;
  } else {
    G4cout <<"MuensterTPCPhysicsList::MuensterTPCPhysicsList() FATAL: Bad EM physics list chosen: "<<m_hEMlowEnergyModel<<G4endl;
    G4String msg = " Available choices are: <emstandard> <emlivermore (default)> <empenelope> <old>";
    G4Exception("MuensterTPCPhysicsList::ConstructProcess()","PhysicsList",FatalException,msg);
  }

	  // add the physics processes
  if (m_hEMlowEnergyModel != "old"){
    emPhysicsList->ConstructProcess();
  } else {
    ConstructEM(); // obsolete in GEANT4_VERSION > geant4.9.4.p02
  }

	  // construct optical physics...... is there a G4 standard for this one as well?
  ConstructOp();
  //opPhysicsList = new G4OpticalPhysics();
  //	opPhysicsList->ConstructProcess();

	// construct the Hadronic physics models
  hadronPhys.clear();
  if (m_hHadronicModel == "custom") {
    // custom hadronic physics list
    ConstructHad();
  } else if (m_hHadronicModel == "QGSP_BERT") {
    // implemented QGSP_BERT: is it done in the right way?
    // this follows the recipe from examples/extended/hadronic/Hadr01
    SetBuilderList1(false);
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BERT());
  } else if (m_hHadronicModel == "QGSP_BERT_HP") {
    // implemented QGSP_BERT_HP: is it done in the right way?
    // this follows the recipe from examples/extended/hadronic/Hadr01
    SetBuilderList1(true);
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BERT_HP());
  } else {
    G4String msg = "MuensterTPCPhysicsList::MuensterTPCPhysicsList() Available choices for Hadronic Physics are: <custom> <QGSP_BERT> <QGSP_BERT_HP>";
    G4Exception("MuensterTPCPhysicsList::ConstructProcess()","PhysicsList",FatalException,msg);
  }
  
  // construct processes
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i]->ConstructProcess();
  }
  
  // some other stuff
  if        (m_hHadronicModel == "custom"){
    ConstructGeneral();
  } else if (m_hHadronicModel == "QGSP_BERT_HP" ||
             m_hHadronicModel == "QGSP_BERT"   ) {
    particleList->ConstructProcess();
    ConstructGeneral();
  }
}

void muensterTPCPhysicsList::SetBuilderList1(G4bool flagHP)
{
  hadronPhys.push_back( new G4EmExtraPhysics());
  if(flagHP) {
    hadronPhys.push_back( new G4HadronElasticPhysicsHP() );
  } else {
    hadronPhys.push_back( new G4HadronElasticPhysics() );
  }
  hadronPhys.push_back( new G4StoppingPhysics());
  hadronPhys.push_back( new G4IonPhysics());
  hadronPhys.push_back( new G4NeutronTrackingCut());
}

// Transportation ///////////////////////////////////////////////////////////

//#include "XeMaxTimeCuts.hh"
//#include "XeMinEkineCuts.hh"

void muensterTPCPhysicsList::AddTransportation()
{
	//G4cout << "----- AddTransportation" << G4endl;
	G4VUserPhysicsList::AddTransportation();

//    GetParticleIterator()->reset();
//    while((*GetParticleIterator()) ())
//    {
//        G4ParticleDefinition *particle = GetParticleIterator()->value();
//        G4ProcessManager *pmanager = particle->GetProcessManager();
//        G4String particleName = particle->GetParticleName();

//        if(particleName == "neutron")
//            pmanager->AddDiscreteProcess(new XeMaxTimeCuts());
//        pmanager->AddDiscreteProcess(new XeMinEkineCuts());
//    }
}

// Electromagnetic Processes ////////////////////////////////////////////////
// all charged particles

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"


// e-
#include "G4eMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"


// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"


// alpha and GenericIon and deuterons, triton, He3:

//muon:
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCapture.hh"

//OTHERS:
#include "G4hIonisation.hh" 
#include "G4hMultipleScattering.hh"
#include "G4hBremsstrahlung.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"

//em process options to allow msc step-limitation to be switched off
#include "G4EmProcessOptions.hh"

void
muensterTPCPhysicsList::ConstructEM()
{
	//G4cout << "----- ConstructEM" << G4endl;
	// Is there a fluorcut missing?
	// at the end the fluorescence is switched to off (end of this chapter)
	//G4double fluorcut = 250 * eV;

  //set a finer grid of the physic tables in order to improve precision
  //former LowEnergy models have 200 bins up to 100 GeV
  G4EmProcessOptions opt;
  opt.SetMaxEnergy(100*GeV);
  opt.SetDEDXBinning(200);
  opt.SetLambdaBinning(200);

	GetParticleIterator()->reset();
	while((*GetParticleIterator()) ())
	{
		G4ParticleDefinition *particle = GetParticleIterator()->value();
		G4ProcessManager *pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		G4String particleType = particle->GetParticleType();
		G4double charge = particle->GetPDGCharge();

		//the Livermore models are not strictly necessary
		if(particleName == "gamma")
		{
			//gamma
			G4RayleighScattering* theRayleigh = new G4RayleighScattering();
			//theRayleigh->SetEmModel(new G4LivermoreRayleighModel());  
			pmanager->AddDiscreteProcess(theRayleigh);

			G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
			//thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
			//thePhotoElectricEffect->SetCutForLowEnSecPhotons(fluorcut);
			pmanager->AddDiscreteProcess(thePhotoElectricEffect);
	
			G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
			//theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
			pmanager->AddDiscreteProcess(theComptonScattering);
	
			G4GammaConversion* theGammaConversion = new G4GammaConversion();
			//theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
			pmanager->AddDiscreteProcess(theGammaConversion);

		}
		else if(particleName == "e-")
		{
			//electron
			// process ordering: AddProcess(name, at rest, along step, post step)
			// -1 = not implemented, then ordering
			// Multiple scattering
			G4eMultipleScattering* msc = new G4eMultipleScattering();
			msc->SetStepLimitType(fUseDistanceToBoundary);
			pmanager->AddProcess(msc,-1, 1, 1);

			// Ionisation
			G4eIonisation* eIonisation = new G4eIonisation();
			//eIonisation->SetCutForLowEnSecPhotons(fluorcut);
			//eIonisation->SetEmModel(new G4LivermoreIonisationModel());
			eIonisation->SetStepFunction(0.2, 100*um); //improved precision in tracking  
			pmanager->AddProcess(eIonisation,-1, 2, 2);
	
			// Bremsstrahlung
			G4eBremsstrahlung* eBremsstrahlung = new G4eBremsstrahlung();
			//eBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());
			//eBremsstrahlung->SetCutForLowEnSecPhotons(fluorcut);
			pmanager->AddProcess(eBremsstrahlung, -1,-1, 3);
		}
		else if(particleName == "e+")
		{
			//positron	
			G4eMultipleScattering* msc = new G4eMultipleScattering();
			msc->SetStepLimitType(fUseDistanceToBoundary);
			pmanager->AddProcess(msc,-1, 1, 1);
	
			// Ionisation
			G4eIonisation* eIonisation = new G4eIonisation();
			eIonisation->SetStepFunction(0.2, 100*um);     
			pmanager->AddProcess(eIonisation, -1, 2, 2);
			//Bremsstrahlung (use default, no low-energy available)
			pmanager->AddProcess(new G4eBremsstrahlung(), -1,-1, 3);

			//Annihilation
			pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);     
		}
		else if(particleName == "mu+" || particleName == "mu-")
		{
			//muon  
			pmanager->AddProcess(new G4eMultipleScattering,-1, 1, 1);
			pmanager->AddProcess(new G4MuIonisation(),-1, 2, 2);
			pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1, 3);
			pmanager->AddProcess(new G4MuPairProduction(),-1,-1, 4);
			if( particleName == "mu-" )
			  pmanager->AddProcess(new G4MuonMinusCapture(), 0,-1,-1);
		}
		// nucleus is not used yet		
		// (particleType == "nucleus" && charge != 0)
    		else if (particleName == "proton" || 
	     		 particleName == "pi+" || 
	     		 particleName == "pi-")
      		{
			//multiple scattering
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      
			//ionisation
			G4hIonisation* hIonisation = new G4hIonisation();
			hIonisation->SetStepFunction(0.2, 50*um);
			pmanager->AddProcess(hIonisation,-1, 2, 2);      
	
			//bremmstrahlung
			pmanager->AddProcess(new G4hBremsstrahlung,-1,-3, 3);
      		}
    		else if (particleName == "alpha"      ||
	    		 particleName == "deuteron"   ||
	     		 particleName == "triton"     ||
	      		 particleName == "He3")
      		{
			//multiple scattering
			pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
	
			//ionisation
			G4ionIonisation* ionIoni = new G4ionIonisation();
			ionIoni->SetStepFunction(0.1, 20*um);
			pmanager->AddProcess(ionIoni,-1, 2, 2);
      		}
    		else if (particleName == "GenericIon")
      		{
			// OBJECT may be dynamically created as either a GenericIon or nucleus
			// G4Nucleus exists and therefore has particle type nucleus
			// genericIon:
	
			//multiple scattering
			pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);

			//ionisation
			G4ionIonisation* ionIoni = new G4ionIonisation();
			ionIoni->SetEmModel(new G4IonParametrisedLossModel());
			ionIoni->SetStepFunction(0.1, 20*um);
			pmanager->AddProcess(ionIoni,-1, 2, 2);
      		}
    		else if ((!particle->IsShortLived()) &&
	     		 (charge != 0.0) && 
	     		 (particle->GetParticleName() != "chargedgeantino")) 
      		{
			//all others charged particles except geantino
        		G4hMultipleScattering* aMultipleScattering = new G4hMultipleScattering();
        		G4hIonisation* ahadronIon = new G4hIonisation();
	
			//multiple scattering
			pmanager->AddProcess(aMultipleScattering,-1,1,1);

			//ionisation
			pmanager->AddProcess(ahadronIon,-1,2,2);      
      		}
	}

  	// turn off msc step-limitation - especially as electron cut 1nm
  	opt.SetMscStepLimitation(fMinimal);

  	// switch on fluorescence, PIXE and Auger:
  	opt.SetFluo(false);
  	opt.SetPIXE(false);
  	opt.SetAuger(false);

}

// Optical Processes ////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Cerenkov.hh"

void
muensterTPCPhysicsList::ConstructOp()
{
	//G4cout << "----- ConstructOp" << G4endl;
	// default scintillation process
	G4Scintillation *theScintProcessDef =
		new G4Scintillation("Scintillation");
	// theScintProcessDef->DumpPhysicsTable();
	theScintProcessDef->SetTrackSecondariesFirst(true);
	theScintProcessDef->SetScintillationYieldFactor(1.0);	
	theScintProcessDef->SetScintillationExcitationRatio(0.0);	
	theScintProcessDef->SetVerboseLevel(OpVerbLevel);

	// scintillation process for alpha:
	G4Scintillation *theScintProcessAlpha =
		new G4Scintillation("Scintillation");
	// theScintProcessNuc->DumpPhysicsTable();
	theScintProcessAlpha->SetTrackSecondariesFirst(true);
	theScintProcessAlpha->SetScintillationYieldFactor(1.1);
	theScintProcessAlpha->SetScintillationExcitationRatio(1.0);
	theScintProcessAlpha->SetVerboseLevel(OpVerbLevel);

	// scintillation process for heavy nuclei
	G4Scintillation *theScintProcessNuc =
		new G4Scintillation("Scintillation");
	// theScintProcessNuc->DumpPhysicsTable();
	theScintProcessNuc->SetTrackSecondariesFirst(true);
	theScintProcessNuc->SetScintillationYieldFactor(0.2);
	theScintProcessNuc->SetScintillationExcitationRatio(1.0);
	theScintProcessNuc->SetVerboseLevel(OpVerbLevel);

  // add Cerenkov
  G4Cerenkov *fCerenkovProcess = new G4Cerenkov("Cerenkov");
  
  if (m_bCerenkov) {
    G4cout <<"muensterTPCPhysicsList::ConstructOp() Define Cerenkov .... "<<G4endl;
    // the maximum numer of photons per GEANT4 step.... 
    G4double fMaxNumPhotons = 100; // same as in G4OpticalPhysics.cc I think, but no guarantees
    // the maximum change in beta=v/c in percent
    G4double fMaxBetaChange = 10;  // same as in G4OpticalPhysics.cc
    // tracks secondaries before continuing with the original particle
    G4bool fTrackSecondariesFirst = true; // same as in G4OpticalPhysics.cc
    
    fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotons);
    fCerenkovProcess->SetMaxBetaChangePerStep(fMaxBetaChange);
    fCerenkovProcess->SetTrackSecondariesFirst(fTrackSecondariesFirst);
  } else {
    G4cout <<"muensterTPCPhysicsList::ConstructOp() Disable Cerenkov .... "<<G4endl;
  }

	// optical processes
	G4OpAbsorption *theAbsorptionProcess = new G4OpAbsorption();
	G4OpRayleigh *theRayleighScatteringProcess = new G4OpRayleigh();
	G4OpBoundaryProcess *theBoundaryProcess = new G4OpBoundaryProcess();

	//theAbsorptionProcess->DumpPhysicsTable();
	//theRayleighScatteringProcess->DumpPhysicsTable();
	theAbsorptionProcess->SetVerboseLevel(OpVerbLevel);
	theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
	theBoundaryProcess->SetVerboseLevel(OpVerbLevel);
 	
	//20th Jun 2012 P.Gumplinger (op-V09-05-04)
        //remove methods: SetModel/GetModel from G4OpBoundaryProcess class
	
	// -> remove:
	//G4OpticalSurfaceModel themodel = unified;
	//theBoundaryProcess->SetModel(themodel);

	GetParticleIterator()->reset();
	while((*GetParticleIterator()) ())
	{
		G4ParticleDefinition *particle = GetParticleIterator()->value();
		G4ProcessManager *pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();

		if(theScintProcessDef->IsApplicable(*particle))
		{
			//      if(particle->GetPDGMass() > 5.0*GeV) 
			if(particle->GetParticleName() == "GenericIon")
			{
				pmanager->AddProcess(theScintProcessNuc);	// AtRestDiscrete
				pmanager->SetProcessOrderingToLast(theScintProcessNuc,
					idxAtRest);
				pmanager->SetProcessOrderingToLast(theScintProcessNuc,
					idxPostStep);
			}
			else if(particle->GetParticleName() == "alpha")
			{
				pmanager->AddProcess(theScintProcessAlpha);
				pmanager->SetProcessOrderingToLast(theScintProcessAlpha,
					idxAtRest);
				pmanager->SetProcessOrderingToLast(theScintProcessAlpha,
					idxPostStep);
			}
			else
			{
				pmanager->AddProcess(theScintProcessDef);
				pmanager->SetProcessOrderingToLast(theScintProcessDef,
					idxAtRest);
				pmanager->SetProcessOrderingToLast(theScintProcessDef,
					idxPostStep);
			}
		}

		if(particleName == "opticalphoton")
		{
			// Step limitation seen as a process
		    pmanager->AddProcess(new G4UserSpecialCuts(),-1,-1,1);
			pmanager->AddDiscreteProcess(theAbsorptionProcess);
			pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
			pmanager->AddDiscreteProcess(theBoundaryProcess);
		}
		// ... and give those particles that need it a bit of Cerenkov.... and only if you want to
    if(fCerenkovProcess->IsApplicable(*particle) && m_bCerenkov){
      pmanager->AddProcess(fCerenkovProcess);
      pmanager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    }
	}
}

// Hadronic processes ////////////////////////////////////////////////////////
// for more information see: DMXPhysicsList.cc (Geant4.10 example source)

// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
//Version check 10.01, 10.01.p01, 10.01.p02, > 10.02
#ifdef G4VERSION_NUMBER <= 1020
	#include "G4GGNuclNuclCrossSection.hh"
#else
	#include "G4ComponentGGNuclNuclXsc.hh"
#endif

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"


void muensterTPCPhysicsList::ConstructHad()
{
	//G4cout << "----- ConstructHad" << G4endl;
	//Elastic models
	const G4double elastic_elimitPi = 1.0*GeV;

	G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
	G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
	elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 
  elastic_he->SetMinEnergy( elastic_elimitPi );

  
  // Inelastic scattering
 	const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  const G4double theHPMin =      0.0*GeV;
 	const G4double theHPMax =     20.0*MeV;

  	G4FTFModel * theStringModel = new G4FTFModel;
  	G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  	theStringModel->SetFragmentationModel( theStringDecay );
  	G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  	G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  	G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  	theFTFModel0->SetHighEnergyGenerator( theStringModel );
  	theFTFModel0->SetTransport( theCascade );
  	theFTFModel0->SetMinEnergy( theFTFMin0 );
  	theFTFModel0->SetMaxEnergy( theFTFMax );

  	G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  	theFTFModel1->SetHighEnergyGenerator( theStringModel );
  	theFTFModel1->SetTransport( theCascade );
  	theFTFModel1->SetMinEnergy( theFTFMin1 );
  	theFTFModel1->SetMaxEnergy( theFTFMax );

  	G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  	theBERTModel0->SetMinEnergy( theBERTMin0 );
  	theBERTModel0->SetMaxEnergy( theBERTMax );

  	G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  	theBERTModel1->SetMinEnergy( theBERTMin1 );
  	theBERTModel1->SetMaxEnergy( theBERTMax );

  	G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  	G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
			//Version check 10.01, 10.01.p01, 10.01.p02, > 10.02
	#ifdef G4VERSION_NUMBER <= 1020
		G4VCrossSectionDataSet * theGGNuclNuclData = G4CrossSectionDataSetRegistry::Instance()->
		  	GetCrossSectionDataSet(G4ComponentGGNuclNuclXsc::Default_Name());
	#else
		// name changed to G4ComponentGGNuclNuclXsc due to v10.02
		G4ComponentGGNuclNuclXsc* ionElasticXS = new G4ComponentGGNuclNuclXsc;
		 	G4VCrossSectionDataSet* theGGNuclNuclData = new G4CrossSectionElastic(ionElasticXS);
	#endif

  	GetParticleIterator()->reset();
  	while ((*GetParticleIterator())()) 
    	{
      		G4ParticleDefinition* particle = GetParticleIterator()->value();
      		G4ProcessManager* pmanager = particle->GetProcessManager();
      		G4String particleName = particle->GetParticleName();

      		if (particleName == "pi+") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          		theElasticProcess->RegisterMe( elastic_lhep1 );
          		theElasticProcess->RegisterMe( elastic_he );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
	  		//Inelastic scattering
	  		G4PionPlusInelasticProcess* theInelasticProcess = 
			  new G4PionPlusInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( thePiData );
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
		} 

      		else if (particleName == "pi-") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          		theElasticProcess->RegisterMe( elastic_lhep1 );
          		theElasticProcess->RegisterMe( elastic_he );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
	  		//Inelastic scattering
	  		G4PionMinusInelasticProcess* theInelasticProcess = 
			  new G4PionMinusInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( thePiData );
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );	  
	  		//Absorption
	  		pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
		}
      
      		else if (particleName == "kaon+") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
          		// Inelastic scattering	
	  		G4KaonPlusInelasticProcess* theInelasticProcess = 
			  new G4KaonPlusInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
			  GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
		}
      
      		else if (particleName == "kaon0S") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
          		// Inelastic scattering	 
	  		G4KaonZeroSInelasticProcess* theInelasticProcess = 
			  new G4KaonZeroSInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
			  GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );	  
			}

      		else if (particleName == "kaon0L") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
	  		// Inelastic scattering
	  		G4KaonZeroLInelasticProcess* theInelasticProcess = 
	    		  new G4KaonZeroLInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
			  GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 ); 
	  		pmanager->AddDiscreteProcess( theInelasticProcess );	  
		}

      		else if (particleName == "kaon-") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
          		// Inelastic scattering
	  		G4KaonMinusInelasticProcess* theInelasticProcess = 
	    		  new G4KaonMinusInelasticProcess("inelastic");	
          		theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
			  GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
	  		pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
		}

      		else if (particleName == "proton") 
		{
	 	 	// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
			  GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          		theElasticProcess->RegisterMe( elastic_chip );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
	  		// Inelastic scattering
	  		G4ProtonInelasticProcess* theInelasticProcess = 
	  		  new G4ProtonInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
		}
      		else if (particleName == "anti_proton") 
		{
	  		// Elastic scattering
          		const G4double elastic_elimitAntiNuc = 100.0*MeV;
          		G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          		elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          		G4CrossSectionElastic* elastic_anucxs = 
			  new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
          		G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
          		elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->AddDataSet( elastic_anucxs );
          		theElasticProcess->RegisterMe( elastic_lhep2 );
          		theElasticProcess->RegisterMe( elastic_anuc );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
	  		// Inelastic scattering
	  		G4AntiProtonInelasticProcess* theInelasticProcess = 
	  		  new G4AntiProtonInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( theAntiNucleonData );
	  		theInelasticProcess->RegisterMe( theFTFModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
	  		// Absorption
	  		pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
		}

      		else if (particleName == "neutron") 
		{
			// elastic scattering
			G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        		theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
			  GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        		G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
			elastic_neutronChipsModel->SetMinEnergy( 19.0*MeV );
        		theElasticProcess->RegisterMe( elastic_neutronChipsModel );
			G4NeutronHPElastic * theElasticNeutronHP = new G4NeutronHPElastic;
        		theElasticNeutronHP->SetMinEnergy( theHPMin );
        		theElasticNeutronHP->SetMaxEnergy( theHPMax );
			theElasticProcess->RegisterMe( theElasticNeutronHP );
			theElasticProcess->AddDataSet( new G4NeutronHPElasticData );
			pmanager->AddDiscreteProcess( theElasticProcess );
			// inelastic scattering		
			G4NeutronInelasticProcess* theInelasticProcess =
	 		  new G4NeutronInelasticProcess("inelastic");
			theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
			theInelasticProcess->RegisterMe( theFTFModel1 );
        		theInelasticProcess->RegisterMe( theBERTModel1 );
			G4NeutronHPInelastic * theNeutronInelasticHPModel = new G4NeutronHPInelastic;
        		theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        		theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
			theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
			theInelasticProcess->AddDataSet( new G4NeutronHPInelasticData );
			pmanager->AddDiscreteProcess(theInelasticProcess);
			// capture
			G4HadronCaptureProcess* theCaptureProcess =
			  new G4HadronCaptureProcess;
			G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
			theLENeutronCaptureModel->SetMinEnergy(theHPMin);
			theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
			theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
			theCaptureProcess->AddDataSet( new G4NeutronHPCaptureData);
			pmanager->AddDiscreteProcess(theCaptureProcess);
		}
      		else if (particleName == "anti_neutron") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
          		// Inelastic scattering (include annihilation on-fly)
	  		G4AntiNeutronInelasticProcess* theInelasticProcess = 
	  		  new G4AntiNeutronInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( theAntiNucleonData );
	  		theInelasticProcess->RegisterMe( theFTFModel0 );
	 	 	pmanager->AddDiscreteProcess( theInelasticProcess );	  
		}

     	 	else if (particleName == "deuteron") 
		{
	 	 	// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
          		// Inelastic scattering
	  		G4DeuteronInelasticProcess* theInelasticProcess = 
	  		  new G4DeuteronInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
		}
      
      		else if (particleName == "triton") 
		{
	  		// Elastic scattering
          		G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
          		// Inelastic scattering
	  		G4TritonInelasticProcess* theInelasticProcess = 
	  		  new G4TritonInelasticProcess("inelastic");
	  		theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
		}
      		else if (particleName == "alpha") 
		{
	  		// Elastic scattering
         	 	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          		theElasticProcess->RegisterMe( elastic_lhep0 );
	  		pmanager->AddDiscreteProcess( theElasticProcess );
          		// Inelastic scattering
	  		G4AlphaInelasticProcess* theInelasticProcess = 
	  	  new G4AlphaInelasticProcess("inelastic");	 
          		theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  		theInelasticProcess->RegisterMe( theFTFModel1 );
          		theInelasticProcess->RegisterMe( theBERTModel0 );
	  		pmanager->AddDiscreteProcess( theInelasticProcess );
		}
    	}
}

// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

void muensterTPCPhysicsList::ConstructGeneral()
{
	//G4cout << "----- ConstructGeneral" << G4endl;
	// Add Decay Process
	G4Decay *theDecayProcess = new G4Decay();

	GetParticleIterator()->reset();
	while((*GetParticleIterator()) ())
	{
		G4ParticleDefinition *particle = GetParticleIterator()->value();
		G4ProcessManager *pmanager = particle->GetProcessManager();

		if(theDecayProcess->IsApplicable(*particle)
			&& !particle->IsShortLived())
		{
			pmanager->AddProcess(theDecayProcess);
			// set ordering for PostStepDoIt and AtRestDoIt
			pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
			pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
		}
	}

	// Declare radioactive decay to the GenericIon in the IonTable.
	const G4IonTable *theIonTable =
		G4ParticleTable::GetParticleTable()->GetIonTable();
	G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();

	for(G4int i = 0; i < theIonTable->Entries(); i++)
	{
		G4String particleName =
			theIonTable->GetParticle(i)->GetParticleName();
		G4String particleType =
			theIonTable->GetParticle(i)->GetParticleType();

		if(particleName == "GenericIon")
		{
			G4ProcessManager *pmanager =
				theIonTable->GetParticle(i)->GetProcessManager();
			pmanager->SetVerboseLevel(VerboseLevel);
			pmanager->AddProcess(theRadioactiveDecay);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
			pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
		}
	}
}

// Cuts /////////////////////////////////////////////////////////////////////
void
muensterTPCPhysicsList::SetCuts()
{

	if(verboseLevel > 1)
		G4cout << "muensterTPCPhysicsList::SetCuts:";

	if(verboseLevel > 0)
	{
		G4cout << "muensterTPCPhysicsList::SetCuts:";
		G4cout << "CutLength : "
			<< G4BestUnit(defaultCutValue, "Length") << G4endl;
	}

	//special for low energy physics
	G4double lowlimit = 250 * eV;

	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,
		100. * GeV);

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma 
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");

	if(verboseLevel > 0) DumpCutValuesTable();
}

