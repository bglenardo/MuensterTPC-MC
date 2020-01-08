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
#ifndef __muensterTPCPPRIMARYGENERATORACTION_H__
#define __muensterTPCPPRIMARYGENERATORACTION_H__

// Additional Header Files
#include <globals.hh>

// Root Header Files
#include <TFile.h>
#include <TH1.h>
#include <TParameter.h>
#include <TRandom3.h>


// G4Header Files
#include <G4EmCalculator.hh>
#include <G4GeneralParticleSource.hh>
#include <G4Material.hh>
#include <G4Navigator.hh>
#include <G4ParticleGun.hh>
#include <G4ThreeVector.hh>
#include <G4VUserPrimaryGeneratorAction.hh>


#include "muensterTPCPrimaryGeneratorMessenger.hh"

class muensterTPCParticleSource;

class G4Event;

class muensterTPCPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
	muensterTPCPrimaryGeneratorAction();
	~muensterTPCPrimaryGeneratorAction();

public:
	const long *GetEventSeeds() { return m_lSeeds; }
	const G4String &GetParticleTypeOfPrimary() { return m_hParticleTypeOfPrimary; }
	G4double GetEnergyOfPrimary() { return m_dEnergyOfPrimary; }
	G4ThreeVector GetPositionOfPrimary() { return m_hPositionOfPrimary; }

    void GeneratePrimariesStandard(G4Event *pEvent);
    void GeneratePrimariesDecay0(G4Event *pEvent);
    void GeneratePrimariesXe131NeutrinoCapture(G4Event *pEvent);
    void GeneratePrimaries(G4Event *pEvent);

	void     SetWriteEmpty(G4bool doit){writeEmpty = doit;};
	G4bool   GetWriteEmpty(){return writeEmpty;};

private:
	muensterTPCPrimaryGeneratorMessenger *m_pMessenger;
	long m_lSeeds[2];
	G4bool	writeEmpty;
	G4String m_hParticleTypeOfPrimary;
	G4double m_dEnergyOfPrimary;
	G4ThreeVector m_hPositionOfPrimary;

  
    
	muensterTPCParticleSource *m_pParticleSource;
};

#endif // __muensterTPCPPRIMARYGENERATORACTION_H__

