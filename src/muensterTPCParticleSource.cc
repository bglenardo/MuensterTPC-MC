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
#include <G4PrimaryParticle.hh>
#include <G4Event.hh>
#include <G4TransportationManager.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4IonTable.hh>
#include <G4Ions.hh>
#include <G4TrackingManager.hh>
#include <G4Track.hh>
#include <Randomize.hh>
#include "TH1.h" 
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <sstream>
#include <cmath>
#include <vector>

using std::stringstream;
using std::vector;

#include "muensterTPCParticleSource.hh"

muensterTPCParticleSource::muensterTPCParticleSource()
{
    first = true;
        stats = true;
    
    m_iNumberOfParticlesToBeGenerated = 1;
	m_pParticleDefinition = 0;
	G4ThreeVector hZero(0., 0., 0.);

	m_hParticleMomentumDirection = G4ParticleMomentum(1., 0., 0.);
    m_hParticleMomentum = new vector<G4ParticleMomentum>;
    m_dParticleEnergy = 1.0*MeV;
	m_hParticlePosition = hZero;
	m_dParticleTime = 0.0;
	m_hParticlePolarization = hZero;
	m_dParticleCharge = 0.0;

    EvID = 0;
    m_hInteractionVertexPosition = new vector<G4ThreeVector>;
    m_hNumberOfProducedParticles = new vector<int>;
    m_hParticleCode = new vector<int>;
    m_hParticleTime = new vector<double>;
    m_hParticleEnergy = new vector<double>;

	m_hSourcePosType = "Volume";
	m_hShape = "NULL";
	m_dHalfz = 0.;
	m_dInnerHalfz = 0.;
	m_dRadius = 0.;
        m_dInnerRadius = 0.;
	m_hCenterCoords = hZero;
	m_bConfine = false;
	m_hVolumeNames.clear();

	m_hAngDistType = "iso";
	m_dMinTheta = 0.;
	m_dMaxTheta = pi;
	m_dMinPhi = 0.;
	m_dMaxPhi = twopi;

	m_hEnergyDisType = "Mono";
	m_dMonoEnergy = 1*MeV;
	m_hEnergyFile = "";
	m_hEnergySpectrum = TH1D("EnergySpectrum", "", 1, 0.999, 1.001);
	m_hEnergySpectrum.SetBinContent(1, 1.);

	m_iVerbosityLevel = 0;

	m_pMessenger = new muensterTPCParticleSourceMessenger(this);
	m_pNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
}

muensterTPCParticleSource::~muensterTPCParticleSource()
{
	delete m_pMessenger;
    delete m_hInteractionVertexPosition;
    delete m_hNumberOfProducedParticles;
    delete m_hParticleCode;
    delete m_hParticleTime;
    delete m_hParticleEnergy;
    delete m_hParticleMomentum;

}

void
muensterTPCParticleSource::SetParticleDefinition(G4ParticleDefinition * aParticleDefinition)
{
	m_pParticleDefinition = aParticleDefinition;
	m_dParticleCharge = m_pParticleDefinition->GetPDGCharge();
}

void
muensterTPCParticleSource::SetEnergyFile(G4String hEnergyFile)
{
	m_hEnergyFile = hEnergyFile;

	ReadEnergySpectrum();
}

G4bool
muensterTPCParticleSource::ReadEnergySpectrum()
{
	// read the energy spectrum from the file
	std::ifstream hIn(m_hEnergyFile.c_str());

	if(hIn.fail())
	{
		G4cout << "Error: cannot open energy spectrum file " << m_hEnergyFile << "!" << G4endl;
		return false;
	}

	if(m_iVerbosityLevel >= 1)
		G4cout << "Source energy spectrum from file: " << m_hEnergyFile << G4endl;

	// read the header
	G4String hEnergyUnit;
	while(!hIn.eof())
	{
		G4String hHeader;
		hIn >> hHeader;

		if(hHeader == "unit:")
			hIn >> hEnergyUnit;
		else if(hHeader == "spectrum:")
			break;
		else
		{
			G4cout << "Error: unknown tag in spectrum file!" << G4endl;
			return false;
		}
	}

	G4double dFactor;
	if(hEnergyUnit == "eV")
		dFactor = MeV/eV;
	else if(hEnergyUnit == "keV")
		dFactor = MeV/keV;
	else if(hEnergyUnit == "MeV")
		dFactor = MeV/MeV;
	else if(hEnergyUnit == "GeV")
		dFactor = MeV/GeV;

	vector<G4double> hEnergyBins;
	vector<G4double> hProbabilities;

	while(!hIn.eof())
	{
		G4double dBinEnergy = 0., dProbability = 0.;

		hIn >> dBinEnergy >> dProbability;

		if(hIn.good())
		{
			if(m_iVerbosityLevel >= 2)
				G4cout << std::setprecision(3) << std::scientific << dBinEnergy << "  " << dProbability << G4endl;
			
			hEnergyBins.push_back(dBinEnergy*dFactor);
			hProbabilities.push_back(dProbability);
		}
	}

	G4int iNbBins = hEnergyBins.size();
	G4double dMin = hEnergyBins.front();
	G4double dMax = hEnergyBins.back();
	G4double dBinWidth = (dMax-dMin)/(iNbBins-1);

	m_hEnergySpectrum.Reset();
	m_hEnergySpectrum.SetBins(iNbBins, dMin-0.5*dBinWidth, dMax+0.5*dBinWidth);
	vector<G4double>::iterator pIt;
	for(pIt = hProbabilities.begin(); pIt != hProbabilities.end(); pIt++)
		m_hEnergySpectrum.SetBinContent((pIt-hProbabilities.begin())+1, *pIt);

	return true;
}

void
muensterTPCParticleSource::ConfineSourceToVolume(G4String hVolumeList)
{
	stringstream hStream;
	hStream.str(hVolumeList);
	G4String hVolumeName;

	// store all the volume names
	while(!hStream.eof())
	{
		hStream >> hVolumeName;
		m_hVolumeNames.insert(hVolumeName);
	}

	// checks if the selected volumes exist and store all volumes that match
	G4PhysicalVolumeStore *PVStore = G4PhysicalVolumeStore::GetInstance();
	G4bool bFoundAll = true;

	set<G4String> hActualVolumeNames;
	for(set<G4String>::iterator pIt = m_hVolumeNames.begin(); pIt != m_hVolumeNames.end(); pIt++)
	{
		G4String hRequiredVolumeName = *pIt;
		G4bool bMatch = false;

		if(bMatch = (hRequiredVolumeName.last('*') != std::string::npos))
			hRequiredVolumeName = hRequiredVolumeName.strip(G4String::trailing, '*');

		G4bool bFoundOne = false;
		for(G4int iIndex = 0; iIndex < (G4int) PVStore->size(); iIndex++)
		{
			G4String hName = (*PVStore)[iIndex]->GetName();

			if((bMatch && (hName.substr(0, hRequiredVolumeName.size())) == hRequiredVolumeName) || hName == hRequiredVolumeName)
			{
				hActualVolumeNames.insert(hName);
				bFoundOne = true;
			}
		}

		bFoundAll = bFoundAll && bFoundOne;
	}

	if(bFoundAll)
	{
		m_hVolumeNames = hActualVolumeNames;
		m_bConfine = true;

		if(m_iVerbosityLevel >= 1)
			G4cout << "Source confined to volumes: " << hVolumeList << G4endl;

		if(m_iVerbosityLevel >= 2)
		{
			G4cout << "Volume list: " << G4endl;

			for(set<G4String>::iterator pIt = m_hVolumeNames.begin(); pIt != m_hVolumeNames.end(); pIt++)
				G4cout << *pIt << G4endl;
		}
	}
	else if(m_hVolumeNames.empty())
		m_bConfine = false;
	else
	{
		G4cout << " **** Error: One or more volumes do not exist **** " << G4endl;
		G4cout << " Ignoring confine condition" << G4endl;
		m_hVolumeNames.clear();
		m_bConfine = false;
	}
}

void
muensterTPCParticleSource::GeneratePointSource()
{
	// Generates Points given the point source.
	if(m_hSourcePosType == "Point")
		m_hParticlePosition = m_hCenterCoords;
	else if(m_iVerbosityLevel >= 1)
		G4cout << "Error SourcePosType is not set to Point" << G4endl;
}

void
muensterTPCParticleSource::GeneratePointsInVolume()
{
	G4ThreeVector RandPos;
	G4double x = 0., y = 0., z = 0.;

	if(m_hSourcePosType != "Volume" && m_iVerbosityLevel >= 1)
		G4cout << "Error SourcePosType not Volume" << G4endl;

	if(m_hShape == "Sphere")
	{
		x = m_dRadius * 2.;
		y = m_dRadius * 2.;
		z = m_dRadius * 2.;
		while(((x * x) + (y * y) + (z * z)) > (m_dRadius * m_dRadius))
		{
			x = G4UniformRand();
			y = G4UniformRand();
			z = G4UniformRand();

			x = (x * 2. * m_dRadius) - m_dRadius;
			y = (y * 2. * m_dRadius) - m_dRadius;
			z = (z * 2. * m_dRadius) - m_dRadius;
		}
	}

	else if(m_hShape == "Cylinder")
	{
		x = m_dRadius * 2.;
		y = m_dRadius * 2.;
		while(((x * x) + (y * y)) > (m_dRadius * m_dRadius))
		{
			x = G4UniformRand();
			y = G4UniformRand();
			z = G4UniformRand();
			x = (x * 2. * m_dRadius) - m_dRadius;
			y = (y * 2. * m_dRadius) - m_dRadius;
			z = (z * 2. * m_dHalfz) - m_dHalfz;
		}
	}
        else if(m_hShape == "Shell")
        {
                x = m_dRadius*2.;
                y = m_dRadius*2.;
                z = G4UniformRand();
                z = (z*2.*m_dHalfz) - m_dHalfz;

                if( (z*z) > (m_dInnerHalfz*m_dInnerHalfz) )
                {
                   while( ( ((x*x) + (y*y)) > (m_dRadius*m_dRadius) ) )
                   {
                           x = G4UniformRand();
                           y = G4UniformRand();
                           x = (x*2.*m_dRadius) - m_dRadius;
                           y = (y*2.*m_dRadius) - m_dRadius;
                   }
                } 
                else
                {
                   while( ( ((x*x) + (y*y)) > (m_dRadius*m_dRadius) ) || 
                          ( ((x*x) + (y*y)) < (m_dInnerRadius*m_dInnerRadius) ) )
                   {
                           x = G4UniformRand();
                           y = G4UniformRand();
                           x = (x*2.*m_dRadius) - m_dRadius;
                           y = (y*2.*m_dRadius) - m_dRadius;
                   }
               } 
                
        }
	else
		G4cout << "Error: Volume Shape Does Not Exist" << G4endl;

	RandPos.setX(x);
	RandPos.setY(y);
	RandPos.setZ(z);
	m_hParticlePosition = m_hCenterCoords + RandPos;

}

G4bool
muensterTPCParticleSource::IsSourceConfined()
{
	// Method to check point is within the volume specified
	if(m_bConfine == false)
		G4cout << "Error: Confine is false" << G4endl;
	G4ThreeVector null(0., 0., 0.);
	G4ThreeVector *ptr;

	ptr = &null;

	// Check m_hParticlePosition is within a volume in our list
	G4VPhysicalVolume *theVolume;

	theVolume = m_pNavigator->LocateGlobalPointAndSetup(m_hParticlePosition, ptr, true);
	G4String theVolName = theVolume->GetName();

	set<G4String>::iterator pIt;
	if((pIt = m_hVolumeNames.find(theVolName)) != m_hVolumeNames.end())
	{
		if(m_iVerbosityLevel >= 1)
			G4cout << "Particle is in volume " << *pIt << G4endl;
		return (true);
	}
	else
		return (false);
}

void
muensterTPCParticleSource::GenerateIsotropicFlux()
{
	G4double rndm, rndm2;
	G4double px, py, pz;

	G4double sintheta, sinphi, costheta, cosphi;

	rndm = G4UniformRand();
	costheta = std::cos(m_dMinTheta) - rndm * (std::cos(m_dMinTheta) - std::cos(m_dMaxTheta));
	sintheta = std::sqrt(1. - costheta * costheta);

	rndm2 = G4UniformRand();
	m_dPhi = m_dMinPhi + (m_dMaxPhi - m_dMinPhi) * rndm2;
	sinphi = std::sin(m_dPhi);
	cosphi = std::cos(m_dPhi);

	px = -sintheta * cosphi;
	py = -sintheta * sinphi;
	pz = -costheta;

	G4double ResMag = std::sqrt((px * px) + (py * py) + (pz * pz));

	px = px / ResMag;
	py = py / ResMag;
	pz = pz / ResMag;

	m_hParticleMomentumDirection.setX(px);
	m_hParticleMomentumDirection.setY(py);
	m_hParticleMomentumDirection.setZ(pz);

	// m_hParticleMomentumDirection now holds unit momentum vector.
	if(m_iVerbosityLevel >= 2)
		G4cout << "Generating isotropic vector: " <<
			m_hParticleMomentumDirection << G4endl;
}

void
muensterTPCParticleSource::GenerateMonoEnergetic()
{
	m_dParticleEnergy = m_dMonoEnergy;
}

void
muensterTPCParticleSource::GenerateEnergyFromSpectrum()
{
	m_dParticleEnergy = m_hEnergySpectrum.GetRandom()*MeV;
}

void
muensterTPCParticleSource::GeneratePrimaryVertex(G4Event * evt)
{

    if (m_hEventInputFile == "fromDecay0File") { //reads events from DECAY0 file
        ReadEventFromDecay0File();

    } else if(m_pParticleDefinition == 0)
	{
		G4cout << "No particle has been defined!" << G4endl;
		return;
	}

	// Position
	G4bool srcconf = false;
	G4int LoopCount = 0;

	while(srcconf == false)
	{
		if(m_hSourcePosType == "Point")
			GeneratePointSource();
		else if(m_hSourcePosType == "Volume")
			GeneratePointsInVolume();
		else
		{
			G4cout << "Error: SourcePosType undefined" << G4endl;
			G4cout << "Generating point source" << G4endl;
			GeneratePointSource();
		}

		if(m_bConfine == true)
		{
			srcconf = IsSourceConfined();
			// if source in confined srcconf = true terminating the loop
			// if source isnt confined srcconf = false and loop continues
		}
		else if(m_bConfine == false)
			srcconf = true;		// terminate loop

		LoopCount++;
		if(LoopCount == 1000000)
		{
			G4cout << "*************************************" << G4endl;
			G4cout << "LoopCount = 1000000" << G4endl;
			G4cout << "Either the source distribution >> confinement" << G4endl;
			G4cout << "or any confining volume may not overlap with" << G4endl;
			G4cout << "the source distribution or any confining volumes" << G4endl;
			G4cout << "may not exist" << G4endl;
			G4cout << "If you have set confine then this will be ignored" << G4endl;
			G4cout << "for this event." << G4endl;
			G4cout << "*************************************" << G4endl;
			srcconf = true;		//Avoids an infinite loop
		}
	}

	// Angular stuff
	if(m_hAngDistType == "iso")
		GenerateIsotropicFlux();
	else if(m_hAngDistType == "direction")
		SetParticleMomentumDirection(m_hParticleMomentumDirection);
	else
		G4cout << "Error: AngDistType has unusual value" << G4endl;
	// Energy stuff
	if(m_hEnergyDisType == "Mono")
		GenerateMonoEnergetic();
	else if(m_hEnergyDisType == "Spectrum")
		GenerateEnergyFromSpectrum();
	else
		G4cout << "Error: EnergyDisType has unusual value" << G4endl;

	if ((m_iVerbosityLevel >= 1) && (m_hEventInputFile != "fromDecay0File"))
	{
		G4cout << "Particle name: " << m_pParticleDefinition->GetParticleName() << G4endl;
		G4cout << "       Energy: " << m_dParticleEnergy << G4endl;
		G4cout << "     Position: " << m_hParticlePosition << G4endl;
		G4cout << "    Direction: " << m_hParticleMomentumDirection << G4endl;
		G4cout << " NumberOfParticlesToBeGenerated: " << m_iNumberOfParticlesToBeGenerated << G4endl;
	}
    
    if(m_hEventInputFile == "fromDecay0File"){     // DECAY0 file
        if(m_iVerbosityLevel >= 1)
        G4cout << "fromDecay0File" << G4endl;
        
        evt->SetEventID(EvID);
        EvID += 1;
        
        //G4cout << "EvID is" << EvID << G4endl;
        
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle_def;
        G4PrimaryVertex *vertex =
        new G4PrimaryVertex(m_hParticlePosition, m_dParticleTime);
        
        for(int ip=0; ip<m_hNumberOfProducedParticles->back(); ip++){
            particle_def = particleTable->FindParticle((*m_hParticleCode)[ip]);
            // Set mass, direction and momentum (kinE)
            G4PrimaryParticle *particle =
            new G4PrimaryParticle(particle_def,
                                  (*m_hParticleMomentum)[ip].x(),
                                  (*m_hParticleMomentum)[ip].y(),
                                  (*m_hParticleMomentum)[ip].z());
            
            if(m_iVerbosityLevel > 1) {
                G4cout << "************* Particle definition (" << ip << "): " << G4endl;
                G4cout << "Particle name: " << particle_def->GetParticleName() << G4endl;
                G4cout << "Particle mass: " << particle->GetMass() << G4endl;
                G4cout << "     Position: " << vertex->GetPosition() << G4endl;
                G4cout << "     Momentum: " << particle->GetMomentum() << G4endl;
            }
            
            vertex->SetPrimary(particle);
        }
        evt->AddPrimaryVertex(vertex);

    } else if( m_hEventInputFile == "Xe131NeutrinoCapture" ) {

        evt->SetEventID(EvID);
        EvID += 1;

	if(m_iVerbosityLevel >= 1)
		G4cout << "Xe131 neutrino caputre generator" << G4endl;

	// First, we generate the Cs131 ion	

        G4PrimaryVertex *vertex = new G4PrimaryVertex(m_hParticlePosition, m_dParticleTime);

        if(m_iVerbosityLevel >= 2)
        G4cout << "Creating primaries and assigning to vertex" << G4endl;
        // create new primaries and set them to the vertex
        G4double ionmass = m_pParticleDefinition->GetPDGMass();
        G4double ionenergy = 0.000001*keV + ionmass;
        G4double ionpmom = std::sqrt(ionenergy * ionenergy - ionmass * ionmass);
        G4double ionpx = ionpmom * m_hParticleMomentumDirection.x();
        G4double ionpy = ionpmom * m_hParticleMomentumDirection.y();
        G4double ionpz = ionpmom * m_hParticleMomentumDirection.z();
        
        //G4cout << "Ionpmom: " << ionpmom << G4endl;

        G4PrimaryParticle *ion = new G4PrimaryParticle(m_pParticleDefinition, ionpx, ionpy, ionpz);
        ion->SetMass(ionmass);
        ion->SetCharge(m_dParticleCharge);
        ion->SetPolarization(m_hParticlePolarization.x(), m_hParticlePolarization.y(), m_hParticlePolarization.z());
	G4cout << "Ion energy: " << ionenergy << G4endl;
	G4cout << "Ion mass: " << ionmass << G4endl;
	G4cout << "Ion pmom: " << ionpmom << G4endl;

        vertex->SetPrimary(ion);

        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle_def;
        
        particle_def = particleTable->FindParticle(11);
        // Set mass, direction and momentum (kinE)
	G4double electronmass = particle_def->GetPDGMass();
        G4double electronenergy = m_dNeutrinoScatterElectronEnergy + electronmass;
	G4double electronpmom = std::sqrt(electronenergy * electronenergy - electronmass * electronmass);
        G4double electronpx = electronpmom * m_hParticleMomentumDirection.x();
        G4double electronpy = electronpmom * m_hParticleMomentumDirection.y();
        G4double electronpz = electronpmom * m_hParticleMomentumDirection.z();
        G4PrimaryParticle *primaryelectron =
        new G4PrimaryParticle(particle_def,
                                  electronpx,
                                  electronpy,
                                  electronpz);
        primaryelectron->SetPolarization(m_hParticlePolarization.x(), m_hParticlePolarization.y(), m_hParticlePolarization.z());
	G4cout << "Electron energy: " << electronenergy << G4endl;
	G4cout << "Electron mass: " << electronmass << G4endl;
	G4cout << "Electron pmom: " << electronpmom << G4endl;
            
            if(m_iVerbosityLevel > 1) {
                G4cout << "************* Particle definition: 11 " << G4endl;
                G4cout << "Particle name: " << particle_def->GetParticleName() << G4endl;
                G4cout << "Particle mass: " << primaryelectron->GetMass() << G4endl;
                G4cout << "     Position: " << vertex->GetPosition() << G4endl;
                G4cout << "     Momentum: " << primaryelectron->GetMomentum() << G4endl;
            }
        //    
        vertex->SetPrimary(primaryelectron);
      	G4cout << "Setting primaryelectron into primary vertex." << G4endl;
        evt->AddPrimaryVertex(vertex);
	G4cout << "Added primary vertex." << G4endl;
        if(m_iVerbosityLevel > 1)
            G4cout << " Primary Vetex generated " << G4endl; 
    
    } else {
        
        G4PrimaryVertex *vertex = new G4PrimaryVertex(m_hParticlePosition, m_dParticleTime);
        
        if(m_iVerbosityLevel >= 2)
        G4cout << "Creating primaries and assigning to vertex" << G4endl;
        // create new primaries and set them to the vertex
        G4double mass = m_pParticleDefinition->GetPDGMass();
        G4double energy = m_dParticleEnergy + mass;
        G4double pmom = std::sqrt(energy * energy - mass * mass);
        G4double px = pmom * m_hParticleMomentumDirection.x();
        G4double py = pmom * m_hParticleMomentumDirection.y();
        G4double pz = pmom * m_hParticleMomentumDirection.z();
        
        //G4cout << m_hParticlePosition << G4endl;

        for(G4int i = 0; i < m_iNumberOfParticlesToBeGenerated; i++)
        {
            G4PrimaryParticle *particle = new G4PrimaryParticle(m_pParticleDefinition, px, py, pz);
            particle->SetMass(mass);
            particle->SetCharge(m_dParticleCharge);
            particle->SetPolarization(m_hParticlePolarization.x(), m_hParticlePolarization.y(), m_hParticlePolarization.z());
            vertex->SetPrimary(particle);
        }
        evt->AddPrimaryVertex(vertex);
        if(m_iVerbosityLevel > 1)
            G4cout << " Primary Vetex generated " << G4endl;
    }
}


void
muensterTPCParticleSource::GeneratePrimaryVertexFromTrack(G4Track *pTrack, G4Event *pEvent)
{
	G4double dPX = pTrack->GetMomentum().x();
	G4double dPY = pTrack->GetMomentum().y();
	G4double dPZ = pTrack->GetMomentum().z();

	G4PrimaryVertex *pVertex = new G4PrimaryVertex(pTrack->GetPosition(), m_dParticleTime);

	G4PrimaryParticle *pPrimary = new G4PrimaryParticle(pTrack->GetDefinition(), dPX, dPY, dPZ);
	pPrimary->SetMass(pTrack->GetDefinition()->GetPDGMass());
	pPrimary->SetCharge(pTrack->GetDefinition()->GetPDGCharge());

	pVertex->SetPrimary(pPrimary);

	pEvent->AddPrimaryVertex(pVertex);
}


void muensterTPCParticleSource::ReadEventFromDecay0File() {
    // Read raw file and determine data structure
    const char* DELIMITER = " ";
    char* token[100] = {}; // initialize to 0
    char* buf = 0;
    
    if ((first == false) && (!raw.good())) {
        raw.close();
        if (m_iVerbosityLevel >= 1) {
            G4cout << "EndOfFile - not enough DECAY0 events!" << G4endl;
        }
        first = true;
    }
    
    // only the first time open file
    if(first){
        
        raw.open(m_hInputFileName);
        //G4cout << "Opened input file." << G4endl;
        // for each event    - event's number, time of event's start,
        //                     number of emitted particles;
        // for each particle - GEANT number of particle,
        //                     (x,y,z) components of momentum,
        //                     time shift from previous time
        // check if "DECAY0" is the first tag
        while ((getline(raw, linebuffer)) &&
               (linebuffer[strspn(linebuffer, " \t\v\r\n")] == '\0')) {
            // Skip blank lines
        }
        buf = strdup(linebuffer.c_str());
        token[0] = strtok(buf, DELIMITER); // first token
        if (strcmp (token[0],"DECAY0") != 0) {
            G4cout << "Error: DECAY0 file has unusual structure" << G4endl;
            return;
        }
        // skip header and read "first event" and "last event"
        while (strcmp (token[0],"First") != 0) {
            getline(raw, linebuffer);
            if (linebuffer[strspn(linebuffer, " \t\v\r\n")] != '\0') { // not blank
                buf = strdup(linebuffer.c_str());
                token[0] = strtok(buf, DELIMITER);
            }
        }
        getline(raw, linebuffer);
        buf = strdup(linebuffer.c_str());
        token[0] = strtok(buf, DELIMITER);
        token[1] = strtok(0, DELIMITER);
        first = false;
        
        if ((m_iVerbosityLevel >= 1) || (stats)) {
            G4cout << "****************** DECAY0 GENERATOR IN ACTION ***************"
            << G4endl;
            G4cout << "********************* Open file " << m_hInputFileName <<  G4endl;
            G4cout << "*************************************************************"
            << G4endl;
            G4cout << "First event:      " << token[0] << G4endl;
            G4cout << "Full # of events: " << token[1] << G4endl;
            G4cout << "*************************************************************"
            << G4endl;
            stats = false;
        }
        
        
        
        while ((getline(raw, linebuffer)) &&
               (linebuffer[strspn(linebuffer, " \t\v\r\n")] == '\0')) {
            // Skip blank lines
        }

    }
    if (m_iVerbosityLevel >= 1) {
        G4cout << "Read Event!" << G4endl;
    }
    m_iNumberOfInteractionSites = 1; // only one interaction site
    m_hNumberOfProducedParticles->clear();
    m_hParticleCode->clear();
    m_hParticleMomentum->clear();
    m_hParticleTime->clear();

    // read event
    buf = strdup(linebuffer.c_str());
    token[0] = strtok(buf, DELIMITER);
    token[1] = strtok(0, DELIMITER);
    token[2] = strtok(0, DELIMITER);
    if (m_iVerbosityLevel >= 1) {
        G4cout << "Event #: " << token[0] << G4endl;
    }
    m_hNumberOfProducedParticles->push_back(atoi(token[2]));

    for(int ip=0; ip<m_hNumberOfProducedParticles->back(); ip++){
        getline(raw, linebuffer);
        buf = strdup(linebuffer.c_str());
        token[0] = strtok(buf, DELIMITER);
        token[1] = strtok(0, DELIMITER);
        token[2] = strtok(0, DELIMITER);
        token[3] = strtok(0, DELIMITER);
        token[4] = strtok(0, DELIMITER);
        m_hParticleCode->push_back(ConvertGeant3toGeant4ParticleCode(atoi(token[0])));
        m_hParticleMomentum->push_back(G4ThreeVector(atof(token[1])*MeV,
                                                     atof(token[2])*MeV, atof(token[3])*MeV));
        m_hParticleTime->push_back(atof(token[4])*s);
    }

    while ((getline(raw, linebuffer)) &&
           (linebuffer[strspn(linebuffer, " \t\v\r\n")] == '\0')) {
        // Skip blank lines
    }
}

int muensterTPCParticleSource::ConvertGeant3toGeant4ParticleCode(int G3ParticleCode) {
    // create several ParticleGuns
    /***************
     DECAY0 particle code - in accordance with GEANT 3.21 manual of October, 1994:
     1 - gamma         2 - positron     3 - electron
     4 - neutrino      5 - muon+        6 - muon-
     7 - pion0         8 - pion+        9 - pion-
     10 - kaon0 long   11 - kaon+       12 - kaon-
     13 - neutron      14 - proton      15 - antiproton
     16 - kaon0 short  17 - eta         18 - lambda
     19 - sigma+       20 - sigma0      21 - sigma-
     22 - xi0          23 - xi-         24 - omega
     25 - antineutron  26 - antilambda  27 - antisigma-
     28 - antisigma0   29 - antisigma+  30 - antixi0
     31 - antixi+      32 - antiomega+  45 - deuteron
     46 - tritium      47 - alpha       48 - geantino
     49 - He3          50 - Cerenkov
     0-dummy
     **************/
    
    switch(G3ParticleCode) {
        case 1   : return 22;       // photon
        case 25  : return -2112;    // anti-neutron
        case 2   : return -11;      // e+
        case 26  : return -3122;    // anti-Lambda
        case 3   : return 11;       // e-
        case 27  : return -3222;    // Sigma-
        case 4   : return 12;       // e-neutrino (NB: flavour undefined by Geant)
        case 28  : return -3212;    // Sigma0
        case 5   : return -13;      // mu+
        case 29  : return -3112;    // Sigma+ (PB)*/
        case 6   : return 13;       // mu-
        case 30  : return -3322;    // Xi0
        case 7   : return 111;      // pi0
        case 31  : return -3312;    // Xi+
        case 8   : return 211;      // pi+
        case 32  : return -3334;    // Omega+ (PB)
        case 9   : return -211;     // pi-
        case 33  : return -15;      // tau+
        case 10  : return 130;      // K long
        case 34  : return 15;       // tau-
        case 11  : return 321;      // K+
        case 35  : return 411;      // D+
        case 12  : return -321;     // K-
        case 36  : return -411;     // D-
        case 13  : return 2112;     // n
        case 37  : return 421;      // D0
        case 14  : return 2212;     // p
        case 38  : return -421;     // D0
        case 15  : return -2212;    // anti-proton
        case 39  : return 431;      // Ds+
        case 16  : return 310;      // K short
        case 40  : return -431;     // anti Ds-
        case 17  : return 221;      // eta
        case 41  : return 4122;     // Lamba_c+
        case 18  : return 3122;     // Lambda
        case 42  : return 24;       // W+
        case 19  : return 3222;     // Sigma+
        case 43  : return -24;      // W-
        case 20  : return 3212;     // Sigma0
        case 44  : return 23;       // Z
        case 21  : return 3112;     // Sigma-
        case 45  : return 0;        // deuteron
        case 22  : return 3322;     // Xi0
        case 46  : return 0;        // triton
        case 23  : return 3312;     // Xi-
        case 47  : return 0;        // alpha
        case 24  : return 3334;     // Omega- (PB)
        case 48  : return 0;        // G nu ? PDG ID 0 is undefined
        default  : return 0;
    }
}

