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
#ifndef __muensterTPCPPARTICLESOURCE_H__
#define __muensterTPCPPARTICLESOURCE_H__

#include <G4VPrimaryGenerator.hh>
#include <G4Navigator.hh>
#include <G4ParticleMomentum.hh>
#include <G4ParticleDefinition.hh>
#include <G4Track.hh>

#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <vector>

//ROOT Header files
#include <TFile.h>
#include "TH1.h"

using std::set;
using std::vector;

#include "muensterTPCParticleSourceMessenger.hh"

class muensterTPCParticleSource: public G4VPrimaryGenerator {
public:
	muensterTPCParticleSource();
	~muensterTPCParticleSource();

public:
	void GeneratePrimaryVertex(G4Event *pEvent);
	void GeneratePrimaryVertexFromTrack(G4Track *pTrack, G4Event *pEvent);

	void SetPosDisType(G4String hSourcePosType) { m_hSourcePosType = hSourcePosType; }
	void SetPosDisShape(G4String hShape) { m_hShape = hShape; }
	void SetCenterCoords(G4ThreeVector hCenterCoords) { m_hCenterCoords = hCenterCoords; }
	void SetHalfZ(G4double dHalfz) { m_dHalfz = dHalfz; }
	void SetInnerHalfZ(G4double dInnerHalfz) { m_dInnerHalfz = dInnerHalfz; }
	void SetRadius(G4double dRadius) { m_dRadius = dRadius; }
	void SetInnerRadius(G4double dInnerRadius) { m_dInnerRadius = dInnerRadius; }

	void SetAngDistType(G4String hAngDistType) { m_hAngDistType = hAngDistType; }
	void SetParticleMomentumDirection(G4ParticleMomentum hMomentum) { m_hParticleMomentumDirection = hMomentum.unit(); }
    void SetParticleMomentum(std::vector<G4ParticleMomentum> *hMomentum) {
        m_hParticleMomentum = hMomentum;
    }
	void SetEnergyDisType(G4String hEnergyDisType) { m_hEnergyDisType = hEnergyDisType; }
	void SetEnergyFile(G4String hEnergyFile);
	void SetMonoEnergy(G4double dMonoEnergy) { m_dMonoEnergy = dMonoEnergy; }

  void SetNumberOfParticlesToBeGenerated(G4int iNumParticles) { m_iNumberOfParticlesToBeGenerated = iNumParticles; }

	void SetParticleDefinition(G4ParticleDefinition *pParticleDefinition);
	inline void SetParticleCharge(G4double dCharge) { m_dParticleCharge = dCharge; }

	void SetVerbosity(G4int iVerbosityLevel) { m_iVerbosityLevel = iVerbosityLevel; }

	const G4String &GetParticleType() { return m_pParticleDefinition->GetParticleName(); }
	const G4double GetParticleEnergy() { return m_dParticleEnergy; }
    inline G4ParticleDefinition *GetParticleDefinition() {
        return m_pParticleDefinition;
    };
    const G4ThreeVector &GetParticlePosition() { return m_hParticlePosition; }
    
    inline std::vector<G4ParticleMomentum>* GetParticleMomentum(){
        return m_hParticleMomentum;
    };
    
	G4bool ReadEnergySpectrum();
	void GeneratePointSource();
	void GeneratePointsInVolume();
	G4bool IsSourceConfined();
	void ConfineSourceToVolume(G4String);

	void GenerateIsotropicFlux();

	void GenerateMonoEnergetic();
	void GenerateEnergyFromSpectrum();
    
    void ReadEventFromDecay0File();
    int ConvertGeant3toGeant4ParticleCode(int G3ParticleCode);
    
    G4String GetEventInputFile() {
        return m_hEventInputFile;
    }
    void SetEventInputFile(G4String hEventInputFile) {
        m_hEventInputFile = hEventInputFile;
    }
    
    void SetInputFileName(G4String hInputFileName) {
        m_hInputFileName = hInputFileName;
    }

    void SetNeutrinoScatterElectronEnergy( G4double m_dNeutrinoScatterElectronEnergy_ ){
	m_dNeutrinoScatterElectronEnergy = m_dNeutrinoScatterElectronEnergy_;
    }
        
        
private:
	G4String m_hSourcePosType;
	G4String m_hShape;
	G4ThreeVector m_hCenterCoords;
	G4double m_dHalfz;
	G4double m_dInnerHalfz;
	G4double m_dRadius;
	G4double m_dInnerRadius;
	G4bool m_bConfine;
	set<G4String> m_hVolumeNames;
	G4String m_hAngDistType;
	G4double m_dMinTheta, m_dMaxTheta, m_dMinPhi, m_dMaxPhi;
	G4double m_dTheta, m_dPhi;
	G4String m_hEnergyDisType;
	G4String m_hEnergyFile;
	G4double m_dMonoEnergy;
    
    std::ifstream raw;
    G4String linebuffer;
    G4String m_hEventInputFile;
    G4String m_hInputFileName;
    G4bool first;
    G4bool stats;
    G4int EvID;
    G4int m_iNumberOfInteractionSites;
    
	G4int m_iNumberOfParticlesToBeGenerated;
	G4ParticleDefinition *m_pParticleDefinition;
	G4ParticleMomentum m_hParticleMomentumDirection;
	G4double m_dParticleEnergy;
	G4double m_dNeutrinoScatterElectronEnergy;
	G4double m_dParticleCharge;
	G4ThreeVector m_hParticlePosition;
	G4double m_dParticleTime;
	G4ThreeVector m_hParticlePolarization;
    
    vector<G4ThreeVector> *m_hInteractionVertexPosition;
    vector<G4ParticleMomentum>    *m_hParticleMomentum;
    vector<int> *m_hNumberOfProducedParticles;
    vector<int> *m_hParticleCode;
    vector<double> *m_hParticleTime;
    vector<double> *m_hParticleEnergy;

	G4int m_iVerbosityLevel;
 	TH1D m_hEnergySpectrum;

	muensterTPCParticleSourceMessenger *m_pMessenger;
	G4Navigator *m_pNavigator;
};

#endif // __muensterTPCPPARTICLESOURCE_H__

