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
#ifndef __muensterTPCPPARTICLESOURCEMESSENGER_H__
#define __muensterTPCPPARTICLESOURCEMESSENGER_H__

#include <G4UImessenger.hh>
#include <globals.hh>

class muensterTPCParticleSource;

class G4ParticleTable;
class G4IonTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class muensterTPCParticleSourceMessenger: public G4UImessenger {
public:
	muensterTPCParticleSourceMessenger(muensterTPCParticleSource *pParticleSource);
	~muensterTPCParticleSourceMessenger();
  
	void SetNewValue(G4UIcommand *pCommand, G4String hNewValues);
 
private:
	muensterTPCParticleSource *m_pParticleSource;
	G4ParticleTable *m_pParticleTable;
	G4IonTable *m_pIonTable;
    
  G4UIdirectory              *m_pDirectory;

  G4UIcmdWithAString         *m_pTypeCmd;
  G4UIcmdWithAString         *m_pShapeCmd;
  G4UIcmdWith3VectorAndUnit  *m_pCenterCmd;
  G4UIcmdWithADoubleAndUnit  *m_pHalfzCmd;
  G4UIcmdWithADoubleAndUnit  *m_pInnerHalfzCmd;
  G4UIcmdWithADoubleAndUnit  *m_pRadiusCmd;
  G4UIcmdWithADoubleAndUnit  *m_pInnerRadiusCmd;
  G4UIcmdWithAString         *m_pConfineCmd;         
  G4UIcmdWithAString         *m_pAngTypeCmd;
  G4UIcmdWithAString         *m_pEnergyTypeCmd;
  G4UIcmdWithAString         *m_pEnergyFileCmd;
  G4UIcmdWithAnInteger       *m_pVerbosityCmd;
  G4UIcommand                *m_pIonCmd;
  G4UIcommand                *m_pXe131NeutrinoCaptureCmd;
  G4UIcmdWithAString         *m_pParticleCmd;
  G4UIcmdWith3VectorAndUnit  *m_pPositionCmd;
  G4UIcmdWith3Vector         *m_pDirectionCmd;
  G4UIcmdWithADoubleAndUnit  *m_pEnergyCmd;
  G4UIcmdWithoutParameter    *m_pListCmd;
  G4UIcmdWithAnInteger       *m_pNumberOfParticlesToBeGeneratedCmd;
  G4UIcmdWithAString *m_pDecay0EventFromFileCmd;

  G4bool   m_bShootIon; 
  G4int    m_iAtomicNumber;
  G4int    m_iAtomicMass;
  G4int    m_iIonCharge;
  G4double m_dIonExciteEnergy; 
  G4double m_dNeutrinoScatterElectronEnergy;
};

#endif // __muensterTPCPPARTICLESOURCEMESSENGER_H__

