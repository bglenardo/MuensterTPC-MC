// include GEANT4 classes
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Orb.hh>
#include <G4Polyhedra.hh>
#include <G4Trd.hh>
#include <G4Cons.hh>
#include <G4UnionSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4SDManager.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VPhysicalVolume.hh>
#include <G4GeometryManager.hh>
#include <G4UnitsTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4UserLimits.hh>
#include <G4RunManager.hh>

// include C++ classes
#include <globals.hh>
#include <vector>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using std::vector;
using std::stringstream;
using std::max;

// include Muenster TPC classes
#include "muensterTPCLXeSensitiveDetector.hh"
#include "muensterTPCPmtSensitiveDetector.hh"
#include "muensterTPCDetectorConstruction.hh"
#include "muensterTPCDetectorMessenger.hh"

map<G4String, G4double> muensterTPCDetectorConstruction::m_hGeometryParameters;

double OuterVesselAlphaChannel = 0.15;
double InnerVesselAlphaChannel = 0.15;
double PTFECylinderAlphaChannel = 0.25;
double CopperRingsAlphaChannel = 0.25;
double InnerFlangeAlphaChannel = 0.75;
double PMTHolderAlphaChannel = 0.75;
double DetectorMaterialAlphaChannel = 0.25;

//******************************************************************/
// Initialization
//******************************************************************/
muensterTPCDetectorConstruction::muensterTPCDetectorConstruction() {
  m_pRotationXPlus225 = new G4RotationMatrix();
  m_pRotationXPlus225->rotateX(22.5*deg);
        
  m_pRotationXMinus225 = new G4RotationMatrix();
  m_pRotationXMinus225->rotateX(-22.5*deg);
        
  m_pRotationXPlus45 = new G4RotationMatrix();
  m_pRotationXPlus45->rotateX(45.*deg);

  m_pRotationXMinus45 = new G4RotationMatrix();
  m_pRotationXMinus45->rotateX(-45.*deg);
        
  m_pRotationXPlus90 = new G4RotationMatrix();
  m_pRotationXPlus90->rotateX(90.*deg);
        
  m_pRotationXMinus90 = new G4RotationMatrix();
  m_pRotationXMinus90->rotateX(-90.*deg);
        
  m_pRotationX180 = new G4RotationMatrix();
  m_pRotationX180->rotateX(180.*deg);

  m_pRotationX0 = new G4RotationMatrix();
  m_pRotationX0->rotateX(0.*deg);
        
  m_pDetectorMessenger = new muensterTPCDetectorMessenger(this);
}

//******************************************************************/
// Destructor
//******************************************************************/
muensterTPCDetectorConstruction::~muensterTPCDetectorConstruction() {
  delete m_pDetectorMessenger;
}

//******************************************************************/
// Constructor
//******************************************************************/
G4VPhysicalVolume* muensterTPCDetectorConstruction::Construct() {

  DefineMaterials();
	LXeMaterial = G4Material::GetMaterial("LXe", false);

  DefineGeometryParameters();

  ConstructLaboratory();

  //ConstructOuterCryostat();

  //ConstructInnerCryostat();

  ConstructLiquidXenonCylinder();

  //ConstructXenon();

  //ConstructTPC();

  //ConstructFieldCage();

  //ConstructPmtArrays();
  
  //PrintPhysicalVolumes();

  //PrintGeometryInformation();
  
	return m_pLabPhysicalVolume;
}

//******************************************************************/
// GetGeometryParameter
//******************************************************************/
G4double muensterTPCDetectorConstruction::GetGeometryParameter(const char *szParameter) {
  return m_hGeometryParameters[szParameter];
}

//******************************************************************/
// DefineMaterials
//******************************************************************/
void muensterTPCDetectorConstruction::DefineMaterials() {
  //================================== elements ===================================
  G4Element *Xe = new G4Element("Xenon",     "Xe", 54., 131.293*g/mole);
  G4Element *H  = new G4Element("Hydrogen",  "H",  1.,  1.0079*g/mole);
  G4Element *C  = new G4Element("Carbon",    "C",  6.,  12.011*g/mole);
  G4Element *N  = new G4Element("Nitrogen",  "N",  7.,  14.007*g/mole);
  G4Element *O  = new G4Element("Oxygen",    "O",  8.,  15.999*g/mole);
  G4Element *F  = new G4Element("Fluorine",  "F",  9.,  18.998*g/mole);
  G4Element *Al = new G4Element("Aluminium", "Al", 13., 26.982*g/mole);
  G4Element *Si = new G4Element("Silicon",   "Si", 14., 28.086*g/mole);
  G4Element *Cr = new G4Element("Chromium",  "Cr", 24., 51.996*g/mole);
  G4Element *Mn = new G4Element("Manganese", "Mn", 25., 54.938*g/mole);
  G4Element *Fe = new G4Element("Iron",      "Fe", 26., 55.85*g/mole);
  G4Element *Ni = new G4Element("Nickel",    "Ni", 28., 58.693*g/mole);
  G4Element *Cu = new G4Element("Copper",    "Cu", 29., 63.546*g/mole);
  G4Element *Pb = new G4Element("Lead",      "Pb", 82., 207.2*g/mole);
  G4Element *Mo = new G4Element("Molybdenum","Mo", 42., 95.96*g/mole);
	
  //================================== materials ================================== 
  G4NistManager* pNistManager = G4NistManager::Instance();

  //------------------------------------- air -------------------------------------
  pNistManager->FindOrBuildMaterial("G4_AIR");

  //----------------------------------- vacuum ------------------------------------
  G4Material *Vacuum = new G4Material("Vacuum", 1.e-20*g/cm3, 2, kStateGas);
  Vacuum->AddElement(N, 0.755);
  Vacuum->AddElement(O, 0.245);
	
  //-------------------------------- liquid xenon ---------------------------------
  G4Material *LXe = new G4Material("LXe", 2.9172*g/cm3, 1, kStateLiquid, 168.15*kelvin, 1.5*atmosphere);
  LXe->AddElement(Xe, 1);

  const G4int iNbEntries = 3;

  G4double pdLXePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV}; //178nm 
  G4double pdLXeScintillation[]    = {0.1,     1.0,     0.1};
  G4double pdLXeRefractiveIndex[]  = {1.63,    1.61,    1.58};
  //measured at some point in the past:
  G4double pdLXeAbsorbtionLength[] = {100.*cm, 100.*cm, 100.*cm}; 
  G4double pdLXeScatteringLength[] = {30.*cm,  30.*cm,  30.*cm};
		
  G4MaterialPropertiesTable *pLXePropertiesTable = new G4MaterialPropertiesTable();
	
  pLXePropertiesTable->AddProperty("FASTCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
  pLXePropertiesTable->AddProperty("SLOWCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
  pLXePropertiesTable->AddProperty("RINDEX", pdLXePhotonMomentum, pdLXeRefractiveIndex, iNbEntries);
  pLXePropertiesTable->AddProperty("ABSLENGTH", pdLXePhotonMomentum, pdLXeAbsorbtionLength, iNbEntries);
  pLXePropertiesTable->AddProperty("RAYLEIGH", pdLXePhotonMomentum, pdLXeScatteringLength, iNbEntries);
	
  pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 0./(21.6*eV));
  pLXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 0);
  pLXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
  pLXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
  pLXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);//ratio btw fast time constant and slow time constant
	
  LXe->SetMaterialPropertiesTable(pLXePropertiesTable);

  //-------------------------------- gaseous xenon --------------------------------
  G4Material *GXe = new G4Material("GXe", 0.005887*g/cm3, 1, kStateGas, 173.15*kelvin, 1.5*atmosphere);
  GXe->AddElement(Xe, 1);

  G4double pdGXePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdGXeScintillation[]    = {0.1,     1.0,     0.1};
  G4double pdGXeRefractiveIndex[]  = {1.00,    1.00,    1.00};
  G4double pdGXeAbsorbtionLength[] = {100*m,   100*m,   100*m};
  G4double pdGXeScatteringLength[] = {100*m,   100*m,   100*m};

  G4MaterialPropertiesTable *pGXePropertiesTable = new G4MaterialPropertiesTable();

  pGXePropertiesTable->AddProperty("FASTCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
  pGXePropertiesTable->AddProperty("SLOWCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
  pGXePropertiesTable->AddProperty("RINDEX", pdGXePhotonMomentum, pdGXeRefractiveIndex, iNbEntries);
  pGXePropertiesTable->AddProperty("ABSLENGTH", pdGXePhotonMomentum, pdGXeAbsorbtionLength, iNbEntries);
  pGXePropertiesTable->AddProperty("RAYLEIGH", pdGXePhotonMomentum, pdGXeScatteringLength, iNbEntries);

  pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 0./(21.6*eV));
  pGXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 0);
  pGXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
  pGXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
  pGXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);

  GXe->SetMaterialPropertiesTable(pGXePropertiesTable);	

  //----------------------------------- quartz ------------------------------------
  // ref: http://www.sciner.com/Opticsland/FS.htm
  G4Material *Quartz = new G4Material("Quartz", 2.201*g/cm3, 2, kStateSolid, 168.15*kelvin, 1.5*atmosphere);
  Quartz->AddElement(Si, 1);
  Quartz->AddElement(O, 2);

  G4double pdQuartzPhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdQuartzRefractiveIndex[]  = {1.50,    1.56,    1.60};
  G4double pdQuartzAbsorbtionLength[] = {30*m,    30*m,    30*m};

  G4MaterialPropertiesTable *pQuartzPropertiesTable = new G4MaterialPropertiesTable();

  pQuartzPropertiesTable->AddProperty("RINDEX", pdQuartzPhotonMomentum, pdQuartzRefractiveIndex, iNbEntries);
  pQuartzPropertiesTable->AddProperty("ABSLENGTH", pdQuartzPhotonMomentum, pdQuartzAbsorbtionLength, iNbEntries);

  Quartz->SetMaterialPropertiesTable(pQuartzPropertiesTable);

  // //------------------------------- stainless steel -------------------------------
  G4Material *SS304LSteel = new G4Material("SS304LSteel", 8.00*g/cm3, 5, kStateSolid);
  SS304LSteel->AddElement(Fe, 0.65);
  SS304LSteel->AddElement(Cr, 0.20);
  SS304LSteel->AddElement(Ni, 0.12);
  SS304LSteel->AddElement(Mn, 0.02);
  SS304LSteel->AddElement(Si, 0.01);

  G4double pdSteelPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdSteelReflectivity[]   = {0.15,    0.2,     0.15};
  G4MaterialPropertiesTable *pSteelPropertiesTable = new G4MaterialPropertiesTable();

  pSteelPropertiesTable->AddProperty("REFLECTIVITY", pdSteelPhotonMomentum, pdSteelReflectivity, iNbEntries);
  SS304LSteel->SetMaterialPropertiesTable(pSteelPropertiesTable);

  //------------------------------- stainless steel -------------------------------
  G4Material *SS316LSteel = new G4Material("SS316LSteel", 8.00*g/cm3, 6, kStateSolid);
  SS316LSteel->AddElement(Fe, 0.682);
  SS316LSteel->AddElement(Cr, 0.172);
  SS316LSteel->AddElement(Ni, 0.109);
  SS316LSteel->AddElement(Mn, 0.016);
  SS316LSteel->AddElement(C, 0.0002);
  SS316LSteel->AddElement(Mo, 0.021);

	//As defined above:
  //G4double pdSteelPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
  //G4double pdSteelReflectivity[]   = {0.15,    0.2,     0.15};
  //G4MaterialPropertiesTable *pSteelPropertiesTable = new G4MaterialPropertiesTable();

  pSteelPropertiesTable->AddProperty("REFLECTIVITY", pdSteelPhotonMomentum, pdSteelReflectivity, iNbEntries);
  SS316LSteel->SetMaterialPropertiesTable(pSteelPropertiesTable);

  //---------------------------------- aluminium ----------------------------------
  pNistManager->FindOrBuildMaterial("G4_Al");

  //---------------------------- photocathode aluminium ---------------------------
  G4Material *PhotoCathodeAluminium = new G4Material("PhotoCathodeAluminium", 8.00*g/cm3, 1, kStateSolid);
  PhotoCathodeAluminium->AddElement(Al, 1);

  G4double pdPhotoCathodePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdPhotoCathodeRefractiveIndex[]  = {1.50,    1.56,    1.60};
  G4double pdPhotoCathodeAbsorbtionLength[] = {1.*nm,   1.*nm,   1.*nm};

  G4MaterialPropertiesTable *pPhotoCathodePropertiesTable = new G4MaterialPropertiesTable();

  pPhotoCathodePropertiesTable->AddProperty("RINDEX", pdPhotoCathodePhotonMomentum, pdPhotoCathodeRefractiveIndex, iNbEntries);
  pPhotoCathodePropertiesTable->AddProperty("ABSLENGTH", pdPhotoCathodePhotonMomentum, pdPhotoCathodeAbsorbtionLength, iNbEntries);

  PhotoCathodeAluminium->SetMaterialPropertiesTable(pPhotoCathodePropertiesTable);

  //------------------------------------- cirlex ---------------------------------
  G4Material *Cirlex = new G4Material("Cirlex", 1.43*g/cm3, 4, kStateSolid); //PMT base material
  Cirlex->AddElement(C, 22);
  Cirlex->AddElement(H, 10);
  Cirlex->AddElement(N, 2);
  Cirlex->AddElement(O, 5);

  G4double pdCirlexPhotonMomentum[]  = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdCirlexReflectivity[]   = {0.5,    0.5,     0.5};

  G4MaterialPropertiesTable *pCirlexPropertiesTable = new G4MaterialPropertiesTable();

  pCirlexPropertiesTable->AddProperty("REFLECTIVITY", pdCirlexPhotonMomentum, pdCirlexReflectivity, iNbEntries);
  Cirlex->SetMaterialPropertiesTable(pCirlexPropertiesTable);
  
	//----------------------------- grid mesh steel in LXe------------------------------
  	G4Material *GridMeshSS316LSteelLXe = new G4Material("GridMeshSS316LSteelLXe", 8.00*g/cm3, 6, kStateSolid);
  GridMeshSS316LSteelLXe->AddElement(Fe, 0.682);
  GridMeshSS316LSteelLXe->AddElement(Cr, 0.172);
  GridMeshSS316LSteelLXe->AddElement(Ni, 0.109);
  GridMeshSS316LSteelLXe->AddElement(Mn, 0.016);
  GridMeshSS316LSteelLXe->AddElement(C, 0.0002);
  GridMeshSS316LSteelLXe->AddElement(Mo, 0.021);

  G4double pdGridMeshPhotonMomentumLXe[] = {6.91*eV, 6.98*eV, 7.05*eV};
  
  G4double *pdGridMeshRefractiveIndexLXe = pdLXeRefractiveIndex;
  G4double dAbsortionLengthLXe = 1.35*mm; // exp(-GridMeshThickness/AbsortionLength) = transparency  lambda=-thickness/ln(transparency) transparency=1-(AreadGrid/AreaVoid)=0.895
  G4double pdGridMeshAbsortionLengthLXe[] = {dAbsortionLengthLXe, dAbsortionLengthLXe, dAbsortionLengthLXe};

 	G4MaterialPropertiesTable *pGridMeshLXePropertiesTable = new G4MaterialPropertiesTable();

 	pGridMeshLXePropertiesTable->AddProperty("RINDEX", pdGridMeshPhotonMomentumLXe, pdGridMeshRefractiveIndexLXe, iNbEntries);
  pGridMeshLXePropertiesTable->AddProperty("ABSLENGTH", pdGridMeshPhotonMomentumLXe, pdGridMeshAbsortionLengthLXe, iNbEntries);
  GridMeshSS316LSteelLXe->SetMaterialPropertiesTable(pGridMeshLXePropertiesTable);

  //----------------------------- grid mesh steel in GXe------------------------------
  G4Material *GridMeshSS316LSteelGXe = new G4Material("GridMeshSS316LSteelGXe", 8.00*g/cm3, 6, kStateSolid);
  GridMeshSS316LSteelGXe->AddElement(Fe, 0.682);
  GridMeshSS316LSteelGXe->AddElement(Cr, 0.172);
  GridMeshSS316LSteelGXe->AddElement(Ni, 0.109);
  GridMeshSS316LSteelGXe->AddElement(Mn, 0.016);
  GridMeshSS316LSteelGXe->AddElement(C, 0.0002);
  GridMeshSS316LSteelGXe->AddElement(Mo, 0.021);

  G4double pdGridMeshPhotonMomentumGXe[] = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double *pdGridMeshRefractiveIndexGXe = pdGXeRefractiveIndex;
  G4double dAbsortionLengthGXe = 1.35*mm; // exp(-GridMeshThickness/AbsortionLength) = transparency  lambda=-thickness/ln(transparency) transparency=1-(AreadGrid/AreaVoid)=0.895
  G4double pdGridMeshAbsortionLengthGXe[] = {dAbsortionLengthGXe, dAbsortionLengthGXe, dAbsortionLengthGXe};

 	G4MaterialPropertiesTable *pGridMeshGXePropertiesTable = new G4MaterialPropertiesTable();

 	pGridMeshGXePropertiesTable->AddProperty("RINDEX", pdGridMeshPhotonMomentumGXe, pdGridMeshRefractiveIndexGXe, iNbEntries);
  pGridMeshGXePropertiesTable->AddProperty("ABSLENGTH", pdGridMeshPhotonMomentumGXe, pdGridMeshAbsortionLengthGXe, iNbEntries);
  GridMeshSS316LSteelGXe->SetMaterialPropertiesTable(pGridMeshGXePropertiesTable);

  //------------------------------------ teflon -----------------------------------
  G4Material* Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
  Teflon->AddElement(C, 0.240183);
  Teflon->AddElement(F, 0.759817);

  G4double pdTeflonPhotonMomentum[]  = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdTeflonRefractiveIndex[] = {1.63,    1.61,    1.58};
  G4double pdTeflonReflectivity[]    = {0.95,    0.95,    0.95};
  G4double pdTeflonSpecularLobe[]    = {0.01,    0.01,    0.01};
  G4double pdTeflonSpecularSpike[]   = {0.01,    0.01,    0.01};
  G4double pdTeflonBackscatter[]     = {0.01,    0.01,    0.01};
  G4double pdTeflonEfficiency[]      = {1.0,     1.0,     1.0};
  G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

  pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, iNbEntries);
  pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, iNbEntries);
  pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularLobe, iNbEntries);
  pTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularSpike, iNbEntries);
  pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter, iNbEntries);
  pTeflonPropertiesTable->AddProperty("EFFICIENCY", pdTeflonPhotonMomentum, pdTeflonEfficiency, iNbEntries);
  Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

  //------------------------------------ GXe teflon -----------------------------------
  G4Material* GXeTeflon = new G4Material("GXeTeflon", 2.2*g/cm3, 2, kStateSolid);
  GXeTeflon->AddElement(C, 0.240183);
  GXeTeflon->AddElement(F, 0.759817);
  
  G4double pdGXeTeflonPhotonMomentum[iNbEntries]  = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdGXeTeflonRefractiveIndex[iNbEntries] = {1.63,    1.61,    1.58};
  G4double pdGXeTeflonReflectivity[iNbEntries]    = {0.95,    0.95,    0.95};
  G4double pdGXeTeflonSpecularLobe[iNbEntries]    = {0.01,    0.01,    0.01};
  G4double pdGXeTeflonSpecularSpike[iNbEntries]   = {0.01,    0.01,    0.01};
  G4double pdGXeTeflonBackscatter[iNbEntries]     = {0.01,    0.01,    0.01};
  G4double pdGXeTeflonEfficiency[iNbEntries]      = {1.0,     1.0,     1.0};
  
  G4MaterialPropertiesTable *pGXeTeflonPropertiesTable = new G4MaterialPropertiesTable();
  
  pGXeTeflonPropertiesTable->AddProperty("RINDEX", pdGXeTeflonPhotonMomentum, pdGXeTeflonRefractiveIndex, iNbEntries);
  pGXeTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdGXeTeflonPhotonMomentum, pdGXeTeflonReflectivity, iNbEntries);
  pGXeTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdGXeTeflonPhotonMomentum, pdGXeTeflonSpecularLobe, iNbEntries);
  pGXeTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", pdGXeTeflonPhotonMomentum, pdGXeTeflonSpecularSpike, iNbEntries);
  pGXeTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdGXeTeflonPhotonMomentum, pdGXeTeflonBackscatter, iNbEntries);
  pGXeTeflonPropertiesTable->AddProperty("EFFICIENCY", pdGXeTeflonPhotonMomentum, pdGXeTeflonEfficiency, iNbEntries);
  
  GXeTeflon->SetMaterialPropertiesTable(pGXeTeflonPropertiesTable);

  //------------------------------------- lead ------------------------------------
  G4Material *Lead = new G4Material("Lead", 11.34*g/cm3, 1);
  Lead->AddElement(Pb, 1); 

  //--------------------------------- polyethylene --------------------------------
  G4Material *Polyethylene = new G4Material("Polyethylene", 0.94*g/cm3, 2, kStateSolid);
  Polyethylene->AddElement(C, 2);
  Polyethylene->AddElement(H, 4);	

  //------------------------------------ copper -----------------------------------
  G4Material *Copper = new G4Material("Copper", 8.92*g/cm3, 1);
  Copper->AddElement(Cu, 1); 

  G4double pdCopperPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
  G4double pdCopperReflectivity[]   = {0.15,    0.2,     0.15};
  G4MaterialPropertiesTable *pCopperPropertiesTable = new G4MaterialPropertiesTable();

  pCopperPropertiesTable->AddProperty("REFLECTIVITY", pdCopperPhotonMomentum, pdCopperReflectivity, iNbEntries);
  Copper->SetMaterialPropertiesTable(pCopperPropertiesTable);
}

//******************************************************************/
// DefineGeometryParameters
//******************************************************************/
void muensterTPCDetectorConstruction::DefineGeometryParameters() {
  //TPC
  m_hGeometryParameters["ScreeningMeshToAnode"]   = 7.85*mm;
  m_hGeometryParameters["AnodeToLiquidLevel"]    = 3.15*mm;
  m_hGeometryParameters["LiquidLevelToGateMesh"]     = 2.0*mm;
  m_hGeometryParameters["LiquidLevelToScreenMesh"]     = 11.15*mm;

  m_hGeometryParameters["CathodeToScreenMesh"]     = 19.4*mm;
  m_hGeometryParameters["TopBottomPMTToScreeningMesh"]     = 3.0*mm;
  m_hGeometryParameters["BottomTopPMTToScreeningMesh"]     = 3.0*mm;

  m_hGeometryParameters["TeflonCentralCylinderInnerRadius"] = 40.0*mm;
  m_hGeometryParameters["TeflonCentralCylinderOuterRadius"] = 50.5*mm;
  m_hGeometryParameters["TeflonCentralCylinderHeight"] = 167.0*mm;
  m_hGeometryParameters["TeflonCentralCylinderThickness"] = 10.5*mm;
  
  m_hGeometryParameters["TeflonOuterCylinderInnerRadius"] = 50.5*mm;
  m_hGeometryParameters["TeflonOuterCylinderOuterRadius"] = 74.0*mm;
  m_hGeometryParameters["TeflonOuterCylinderHeight"] = 167.0*mm;
  m_hGeometryParameters["TeflonOuterCylinderThickness"] = 23.5*mm;
  
  //bottom part TPC
  m_hGeometryParameters["BottomPMTHolderTopPartThickness"] = 15.4*mm; // goes from middle of pmt to top of pmt NOT to the cathode	
  m_hGeometryParameters["BottomPMTHolderBottomPartThickness"] = 16.6*mm;	
  m_hGeometryParameters["BottomPMTBaseHolderThickness"] = 12.0*mm;	
  m_hGeometryParameters["BottomPMTHolderPlateRadius"]    = 74.0*mm;	
  
  m_hGeometryParameters["BottomPTFESlabThickness"] = 24.4*mm; //slab between top of bottom pmt and the cathode
  m_hGeometryParameters["BottomPTFESlabInnerRadius"] = 40.0*mm; 
  m_hGeometryParameters["BottomPTFESlabOuterRadius"] = 74.0*mm; 

  m_hGeometryParameters["BottomSteelRingThickness"] = 2.0*mm; //slab between top of bottom pmt and the cathode
  m_hGeometryParameters["BottomSteelRingInnerRadius"] = 47.5*mm; 
  m_hGeometryParameters["BottomSteelRingOuterRadius"] = 74.0*mm; 
  m_hGeometryParameters["BottomFlangeCopperRingToTopSteelRing"]=10.2*mm;

  //top part TPC
  m_hGeometryParameters["TopPMTHolderPlateThickness"] = 34.6*mm;
  m_hGeometryParameters["TopPMTHolderTopPartThickness"] = 16.6*mm; // goes from middle of pmt to top of pmt NOT to the cathode	
  m_hGeometryParameters["TopPMTHolderBottomPartThickness"] = 18.0*mm;
  m_hGeometryParameters["TopPMTHolderPlateRadius"]    = 74.0*mm;

  m_hGeometryParameters["TopPTFESlabThickness"] = 13.5*mm; //slab between top of bottom pmt and the cathode
  m_hGeometryParameters["TopPTFESlabInnerRadius"] = 40.0*mm; 
  m_hGeometryParameters["TopPTFESlabOuterRadius"] = 74.0*mm; 

  //copper rings
  m_hGeometryParameters["CopperRingsOuterRadius"] = 60.25*mm;
  m_hGeometryParameters["CopperRingsInnerRadius"] = 50.9*mm;
  m_hGeometryParameters["CopperRingsThickness"] = 2.0*mm;
  m_hGeometryParameters["CopperRingsNumber"] = 12;
  m_hGeometryParameters["CopperRingsSpacing"] = 11.0*mm;

  //meshes
  m_hGeometryParameters["GridMeshThickness"] = 0.15*mm;
  m_hGeometryParameters["GridMeshRadius"] = 50.*mm;
  m_hGeometryParameters["GridMeshNumber"] = 5;

  m_hGeometryParameters["GridMeshBorderInnerRadius"] = 50.*mm;
  m_hGeometryParameters["GridMeshBorderOuterRadius"] = 55*mm;

  m_hGeometryParameters["GateMeshRingThickness"] = 2.0*mm; //in reality this is the support

  m_hGeometryParameters["GridMeshSupportThickness"] = 2.0*mm;
  m_hGeometryParameters["GridMeshSupportInnerRadius"] = 50.9*mm;
  m_hGeometryParameters["GridMeshSupportOuterRadius"] = 60.25*mm;
  m_hGeometryParameters["GridMeshSupportNumber"] = 4;

  m_hGeometryParameters["CathodeGridMeshSupportThickness"] = 2.0*mm;
  m_hGeometryParameters["CathodeGridMeshSupportInnerRadius"] = 45.0*mm;
  m_hGeometryParameters["CathodeGridMeshSupportOuterRadius"] = 60.25*mm;

  //pmt
  m_hGeometryParameters["PmtWidth"]                 = 25.4*mm;
  m_hGeometryParameters["PmtSpacing"]               = 2.0*mm;
  m_hGeometryParameters["PmtWindowWidth"]           = 25.00*mm;
  m_hGeometryParameters["PmtWindowThickness"]       = 1.50*mm;
  m_hGeometryParameters["PmtCasingWidth"]           = 25.40*mm;
  m_hGeometryParameters["PmtCasingHeight"]          = 27.00*mm;
  m_hGeometryParameters["PmtCasingThickness"]       = 0.50*mm;
  m_hGeometryParameters["PmtPhotoCathodeWidth"]     = 22.00*mm;
  m_hGeometryParameters["PmtPhotoCathodeThickness"] = 0.50*mm;
  m_hGeometryParameters["PmtBaseThickness"]         = 1.50*mm;
  m_hGeometryParameters["PmtToPmtBase"]             = 3.00*mm;

  m_hGeometryParameters["NbTopPmtsFirstRow"]   = 2;
  m_hGeometryParameters["NbTopPmtsSecondRow"]  = 3;
  m_hGeometryParameters["NbTopPmtsThirdRow"]  = 2;
  m_hGeometryParameters["NbTopPmts"] = 7;

  m_hGeometryParameters["NbBottomPmtsFirstRow"]   = 2;
  m_hGeometryParameters["NbBottomPmtsSecondRow"]  = 3;
  m_hGeometryParameters["NbBottomPmtsThirdRow"]  = 2;
  m_hGeometryParameters["NbBottomPmts"] = 7;

  m_hGeometryParameters["SpaceBelowTopPMTHeight"] = 3.05*mm;

  // inner cryostat
  m_hGeometryParameters["InnerCryostatBlindFlangeThickness"]     = 22.2*mm;
  m_hGeometryParameters["InnerCryostatBoreFlangeThickness"]     = 22.2*mm;
  m_hGeometryParameters["InnerCryostatBoreFlangeOuterRadius"]     = 101.15*mm;
  m_hGeometryParameters["InnerCryostatBlindFlangeOuterRadius"]     = 101.15*mm;
  m_hGeometryParameters["BoreFlangeNumber"] =2 ;
  m_hGeometryParameters["BlindFlangeNumber"] =2 ;

  m_hGeometryParameters["InnerCryostatVesselOuterRadius"]     = 77.5*mm;
  m_hGeometryParameters["InnerCryostatVesselInnerRadius"]     = 74.0*mm;
  m_hGeometryParameters["InnerCryostatVesselHeight"]     = 334.0*mm;
  m_hGeometryParameters["InnerCryostatVesselThickness"]     = 3.5*mm;

  m_hGeometryParameters["CryostatToTopPMTCasingSpacing"] = 18.45*mm;

  m_hGeometryParameters["CryostatBottomDistanceBetweenLXeAndPTFE"]=10.2*mm; //this small distance at the bottom of the cryostat where you have LXe and screws etc etc

  m_hGeometryParameters["TopBoreFlangeInnerCryostatToLiquidLevel"]     = 86.4*mm;

  m_hGeometryParameters["InnerCryostatCopperRingsOuterRadius"]     = 85.5*mm;
  m_hGeometryParameters["InnerCryostatCopperRingsInnerRadius"]     = 76.5*mm;
  m_hGeometryParameters["InnerCryostatCopperRingsThickness"]     = 2.0*mm;
  m_hGeometryParameters["InnerCryostatCopperRingsNumber"] =2 ;

  // outer cryostat
  m_hGeometryParameters["OuterCryostatBlindFlangeThickness"]     = 16.0*mm;
  m_hGeometryParameters["OuterCryostatBoreFlangeThickness"]     = 16.0*mm;
  m_hGeometryParameters["OuterCryostatBoreFlangeOuterRadius"]     = 167.5*mm;
  m_hGeometryParameters["OuterCryostatBoreFlangeInnerRadius"]     = 130.5*mm;
  m_hGeometryParameters["OuterCryostatBlindFlangeOuterRadius"]     = 167.5*mm;
  m_hGeometryParameters["OuterCryostatBoreFlangeNumber"] =2 ;
  m_hGeometryParameters["OuterCryostatBlindFlangeNumber"] =2 ;

	m_hGeometryParameters["OuterCryostatVesselOuterRadius"]     = 136.5*mm;
  m_hGeometryParameters["OuterCryostatVesselInnerRadius"]     = 133.5*mm;
  m_hGeometryParameters["OuterCryostatVesselHeight"]     = 1512.0*mm; //inc bore flanges
  m_hGeometryParameters["OuterCryostatVesselThickness"]     = 3.0*mm;

  m_hGeometryParameters["OuterCryostatToTurboPumpTubeInnerRadius"]=50.0*mm;
  m_hGeometryParameters["OuterCryostatToTurboPumpTubeOuterRadius"]=52.0*mm;
  //m_hGeometryParameters["OuterCryostatToTurboPumpTubeHeight"]=
  m_hGeometryParameters["OuterCryostatToTurboPumpTubeFlangeInnerRadius"]=50.0*mm;
  m_hGeometryParameters["OuterCryostatToTurboPumpTubeFlangeOuterRadius"]=75.75*mm;
  m_hGeometryParameters["OuterCryostatToTurboPumpTubeFlangeThickness"]=20.0*mm;

  m_hGeometryParameters["OuterCryostatBottomToTurboPumpTube"]=242.0*mm;
	m_hGeometryParameters["OuterCryostatTopToTopFlangeTPCCryostat"]=795.0*mm;
 	m_hGeometryParameters["OuterCryostatBottomToBottomFlangeTPCCryostat"]=264.6*mm;

  // derived quantities
  //LXeHeight is defined as the cryostat so that the LXe can be used as mother volume
 	m_hGeometryParameters["LXeHeight"] = GetGeometryParameter("InnerCryostatVesselHeight");

  m_hGeometryParameters["GXeHeight"] = GetGeometryParameter("InnerCryostatBoreFlangeThickness")+GetGeometryParameter("CryostatToTopPMTCasingSpacing")+GetGeometryParameter("TopPMTHolderPlateThickness")+GetGeometryParameter("ScreeningMeshToAnode")+ GetGeometryParameter("AnodeToLiquidLevel")+0.01;

  m_hGeometryParameters["LXeLevelZeroOffset"] = 0.5*GetGeometryParameter("LXeHeight")-GetGeometryParameter("InnerCryostatBoreFlangeThickness")-GetGeometryParameter("CryostatToTopPMTCasingSpacing")-GetGeometryParameter("TopPMTHolderPlateThickness")-GetGeometryParameter("ScreeningMeshToAnode")-GetGeometryParameter("AnodeToLiquidLevel");//LXe is usually mother volume. Defining liquid level as a zero offset, a reference point for the calculation of all the other offsets when the LXe is the mother volume

  const G4double dLXeLevelZeroOffset =GetGeometryParameter("LXeLevelZeroOffset");
  //G4cout<<"dLXeLevelZeroOffset: " << dLXeLevelZeroOffset << G4endl;

   m_hGeometryParameters["InnerCryostatVesselOffsetZ"] = -(0.5*GetGeometryParameter("OuterCryostatVesselHeight")-GetGeometryParameter("OuterCryostatBottomToBottomFlangeTPCCryostat")-GetGeometryParameter("InnerCryostatBlindFlangeThickness")-0.5*GetGeometryParameter("InnerCryostatVesselHeight")); //wrt center of outer cryostat

  m_hGeometryParameters["LXeLevelToTopTPCCryostat"] = GetGeometryParameter("AnodeToLiquidLevel")+ GetGeometryParameter("ScreeningMeshToAnode")+ GetGeometryParameter("TopPMTHolderPlateThickness")+ GetGeometryParameter("CryostatToTopPMTCasingSpacing")+ GetGeometryParameter("InnerCryostatBoreFlangeThickness");

  // verifications
  // assert(GetGeometryParameter("TeflonPanelHeight") == GetGeometryParameter("DriftLength")-GetGeometryParameter("GridRingThickness"));
}

//******************************************************************/
// ConstructLaboratory
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructLaboratory() {
  //---------------------------------- laboratory ---------------------------------
  const G4double dLabHalfX = 1.5*m;
  const G4double dLabHalfY = 1.5*m;
  const G4double dLabHalfZ = 1.5*m;

  G4Box *pLabBox = new G4Box("LabBox", dLabHalfX, dLabHalfY, dLabHalfZ);
	
  G4Material *Air = G4Material::GetMaterial("G4_AIR");

  m_pLabLogicalVolume = new G4LogicalVolume(pLabBox, Air, "LabLogicalVolume", 0, 0, 0);

  //m_pLabLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

  G4Colour LabColor(1.0,0.0,0.0,1.0); //red
  G4VisAttributes *pLabVisAtt = new G4VisAttributes(LabColor);
  //G4VisAttributes *pLabVisAtt = new G4VisAttributes();
  pLabVisAtt->SetVisibility(false);
  m_pLabLogicalVolume->SetVisAttributes(pLabVisAtt);

  m_pLabPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(), m_pLabLogicalVolume, "Lab", 0, false, 0);
	
  //m_pMotherLogicalVolume = m_pLabLogicalVolume;
}

//******************************************************************/
// ConstructOuterCryostat
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructOuterCryostat() {
   G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");
   G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
   G4Material *Teflon = G4Material::GetMaterial("Teflon");
   G4Material *SS316LSteel = G4Material::GetMaterial("SS316LSteel");
   G4Material *Copper = G4Material::GetMaterial("Copper");
   
   const G4double dLXeHeight = GetGeometryParameter("LXeHeight");
   const G4double dLXeHalfZ = 0.5*dLXeHeight;
   const G4double dOuterCryostatBlindFlangeThickness = GetGeometryParameter("OuterCryostatBlindFlangeThickness");
   const G4double dOuterCryostatBoreFlangeThickness = GetGeometryParameter("OuterCryostatBoreFlangeThickness");
   const G4double dOuterCryostatBoreFlangeOuterRadius = GetGeometryParameter("OuterCryostatBoreFlangeOuterRadius");
   const G4double dOuterCryostatBoreFlangeInnerRadius = GetGeometryParameter("OuterCryostatBoreFlangeInnerRadius");
   const G4double dOuterCryostatBlindFlangeOuterRadius = GetGeometryParameter("OuterCryostatBlindFlangeOuterRadius");
   const G4int dOuterCryostatBoreFlangeNumber = GetGeometryParameter("OuterCryostatBoreFlangeNumber");
   const G4int dOuterCryostatBlindFlangeNumber = GetGeometryParameter("OuterCryostatBlindFlangeNumber");
   const G4double dOuterCryostatVesselOuterRadius = GetGeometryParameter("OuterCryostatVesselOuterRadius");
   const G4double dOuterCryostatVesselInnerRadius = GetGeometryParameter("OuterCryostatVesselInnerRadius");
   const G4double dOuterCryostatVesselHeight = GetGeometryParameter("OuterCryostatVesselHeight");
   const G4double dOuterCryostatVesselThickness = GetGeometryParameter("OuterCryostatVesselThickness");
   const G4double dOuterCryostatToTurboPumpTubeInnerRadius = GetGeometryParameter("OuterCryostatToTurboPumpTubeInnerRadius");
   const G4double dOuterCryostatToTurboPumpTubeOuterRadius = GetGeometryParameter("OuterCryostatToTurboPumpTubeOuterRadius");
   const G4double dOuterCryostatToTurboPumpTubeFlangeInnerRadius = GetGeometryParameter("OuterCryostatToTurboPumpTubeFlangeInnerRadius");
   const G4double dOuterCryostatToTurboPumpTubeFlangeOuterRadius = GetGeometryParameter("OuterCryostatToTurboPumpTubeFlangeOuterRadius");
   const G4double dOuterCryostatToTurboPumpTubeFlangeThickness= GetGeometryParameter("OuterCryostatToTurboPumpTubeFlangeThickness");
   
   const G4double dOuterCryostatBottomToTurboPumpTube = GetGeometryParameter("OuterCryostatBottomToTurboPumpTube");
   const G4double dOuterCryostatTopToTopFlangeTPCCryostat = GetGeometryParameter("OuterCryostatTopToTopFlangeTPCCryostat");
   const G4double dOuterCryostatBottomToBottomFlangeTPCCryostat = GetGeometryParameter("OuterCryostatBottomToBottomFlangeTPCCryostat");

   const G4double dInnerCryostatBlindFlangeThickness = GetGeometryParameter("InnerCryostatBlindFlangeThickness");

   const G4double dLXeLevelZeroOffset = GetGeometryParameter("LXeLevelZeroOffset");
   const G4double dScreeningMeshToAnode= GetGeometryParameter("ScreeningMeshToAnode");
   const G4double dAnodeToLiquidLevel= GetGeometryParameter("AnodeToLiquidLevel");
   const G4double dTopPMTHolderPlateThickness= GetGeometryParameter("TopPMTHolderPlateThickness");
   const G4double dCryostatToTopPMTCasingSpacing= GetGeometryParameter("CryostatToTopPMTCasingSpacing");
 
   const G4double dLXeLevelToTopTPCCryostat= GetGeometryParameter("LXeLevelToTopTPCCryostat");
   
   //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< outer cryostat vessel>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
   const G4double dOuterCryostatVesselHalfZ=0.5*dOuterCryostatVesselHeight;
   
   const G4double dOuterCryostatVesselOffsetZ=(dOuterCryostatVesselHalfZ-dOuterCryostatBottomToBottomFlangeTPCCryostat-dInnerCryostatBlindFlangeThickness-dLXeHeight+dLXeLevelToTopTPCCryostat);

   //G4cout << "OuterCryostatVesselOffsetZ: " << dOuterCryostatVesselOffsetZ << G4endl;
   
   G4Tubs *pOuterCryostatVesselTubs = new G4Tubs("OuterCryostatVesselTubs", dOuterCryostatVesselInnerRadius,dOuterCryostatVesselOuterRadius,dOuterCryostatVesselHalfZ , 0.*deg, 360.*deg);
   
   m_pOuterCryostatVesselLogicalVolume = new G4LogicalVolume(pOuterCryostatVesselTubs, SS304LSteel, "OuterCryostatVesselLogicalVolume", 0, 0, 0);
   
   m_pOuterCryostatVesselPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dOuterCryostatVesselOffsetZ),
							   m_pOuterCryostatVesselLogicalVolume, "OuterCryostatVessel", m_pLabLogicalVolume, false, 0);

	// =========================== bore flanges  ====================================
  const G4double dOuterCryostatBoreFlangeHalfZ = dOuterCryostatBoreFlangeThickness*0.5;
  
  const G4double dOuterCryostatBottomBoreFlangeOffsetZ = (dOuterCryostatVesselOffsetZ-dOuterCryostatVesselHalfZ+dOuterCryostatBoreFlangeHalfZ);
   //const G4double dOuterCryostatBottomBoreFlangeOffsetZ = -(dLXeHalfZ+dInnerCryostatBlindFlangeThickness+dOuterCryostatBottomToBottomFlangeTPCCryostat+dOuterCryostatBoreFlangeHalfZ);// -(dLXeHalfZ+dTopBoreFlangeOuterCryostatToLiquidLevel-dOuterCryostatBoreFlangeHalfZ);

  G4Tubs *pOuterCryostatBoreFlangeTubs = new G4Tubs("OuterCryostatBoreFlangeTubs", dOuterCryostatVesselInnerRadius,dOuterCryostatBoreFlangeOuterRadius,dOuterCryostatBoreFlangeHalfZ , 0.*deg, 360.*deg);

  m_pOuterCryostatBoreFlangeLogicalVolume = new G4LogicalVolume(pOuterCryostatBoreFlangeTubs, SS316LSteel, "OuterCryostatBoreFlangeLogicalVolume", 0, 0, 0);
       
  G4double pdOuterCryostatBoreFlangeZ[dOuterCryostatBoreFlangeNumber];
	
  for(int i = 0; i < dOuterCryostatBoreFlangeNumber; i++)
    {
	  
      pdOuterCryostatBoreFlangeZ[i]=dOuterCryostatBottomBoreFlangeOffsetZ+i*(dOuterCryostatVesselHeight-2*dOuterCryostatBoreFlangeHalfZ);//334-22.2;

      m_hOuterCryostatBoreFlangePhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdOuterCryostatBoreFlangeZ[i]),
			   m_pOuterCryostatBoreFlangeLogicalVolume,"OuterCryostatBoreFlange",m_pLabLogicalVolume, false, i));

      //G4cout<<"pdOuterCryostatBoreFlangeZ[i]: "<<pdOuterCryostatBoreFlangeZ[i]<<G4endl;
                                                                                    
    }

  // ======================== blind flanges ==========================================
  const G4double dOuterCryostatBlindFlangeHalfZ = dOuterCryostatBlindFlangeThickness*0.5;
  
	const G4double dOuterCryostatBottomBlindFlangeOffsetZ = dOuterCryostatBottomBoreFlangeOffsetZ-dOuterCryostatBoreFlangeHalfZ-dOuterCryostatBlindFlangeHalfZ;//dOuterCryostatVesselOffsetZ-dOuterCryostatVesselHalfZ-dOuterCryostatBlindFlangeHalfZ;

//const G4double dOuterCryostatBottomBlindFlangeOffsetZ = -(dLXeHalfZ+dInnerCryostatBlindFlangeThickness+dOuterCryostatBottomToBottomFlangeTPCCryostat+dOuterCryostatBoreFlangeThickness+dOuterCryostatBlindFlangeHalfZ);

  G4Tubs *pOuterCryostatBlindFlangeTubs = new G4Tubs("OuterCryostatBlindFlangeTubs", 0,dOuterCryostatBlindFlangeOuterRadius,dOuterCryostatBlindFlangeHalfZ , 0.*deg, 360.*deg);

  m_pOuterCryostatBlindFlangeLogicalVolume = new G4LogicalVolume(pOuterCryostatBlindFlangeTubs, SS316LSteel, "OuterCryostatBlindFlangeLogicalVolume", 0, 0, 0);
       
  G4double pdOuterCryostatBlindFlangeZ[dOuterCryostatBlindFlangeNumber];
	
  for(int i = 0; i < dOuterCryostatBlindFlangeNumber; i++)
    {
	  
      pdOuterCryostatBlindFlangeZ[i]=dOuterCryostatBottomBlindFlangeOffsetZ+i*(dOuterCryostatVesselHeight+2*dOuterCryostatBlindFlangeHalfZ);

      m_hOuterCryostatBlindFlangePhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdOuterCryostatBlindFlangeZ[i]),
			   m_pOuterCryostatBlindFlangeLogicalVolume,"OuterCryostatBlindFlange",m_pLabLogicalVolume, false, i));

      //G4cout<<"pdOuterCryostatBlindFlangeZ[i]: "<<pdOuterCryostatBlindFlangeZ[i]<<G4endl;
                                                                                    
    }

	// ======================== vacuum between cryostats==========================================
  G4Tubs *pOuterCryostatVacuumTubs = new G4Tubs("OuterCryostatVacuumTubs",0,dOuterCryostatVesselInnerRadius,dOuterCryostatVesselHalfZ , 0.*deg, 360.*deg);
   
  m_pOuterCryostatVacuumLogicalVolume = new G4LogicalVolume(pOuterCryostatVacuumTubs, Vacuum, "OuterCryostatVacuumLogicalVolume", 0, 0, 0);
   
  m_pOuterCryostatVacuumPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dOuterCryostatVesselOffsetZ), m_pOuterCryostatVacuumLogicalVolume, "OuterCryostatVacuum", m_pLabLogicalVolume, false, 0);

  //================================== attributes =================================
  G4Colour hStainlessSteelColor(0.500, 0.500, 0.500, 0.1);
  G4Colour hVacuumColor(0.600, 0.200, 0.250, 0.1);
  G4Colour hRedColor(1.,0., 0., 0.5);
  
  G4Colour hOuterCryostatVesselColor(0.,0.392,0,OuterVesselAlphaChannel);  //0 100 0 dark green
  G4VisAttributes *pOuterCryostatVesselColorVisAtt = new G4VisAttributes(hOuterCryostatVesselColor);
  pOuterCryostatVesselColorVisAtt->SetVisibility(true);
  m_pOuterCryostatVesselLogicalVolume->SetVisAttributes(pOuterCryostatVesselColorVisAtt);

  G4Colour hOuterCryostatBoreFlangeColor(1.0,0.549, 0);  //255 140 0 orange
  G4VisAttributes *pOuterCryostatBoreFlangeColorVisAtt = new G4VisAttributes(hOuterCryostatBoreFlangeColor);
  pOuterCryostatBoreFlangeColorVisAtt->SetVisibility(true);
  m_pOuterCryostatBoreFlangeLogicalVolume->SetVisAttributes(pOuterCryostatBoreFlangeColorVisAtt);

  G4Colour hOuterCryostatBlindFlangeColor(0.545,0.353, 0);  //139 90 0 dark orange
  G4VisAttributes *pOuterCryostatBlindFlangeColorVisAtt = new G4VisAttributes(hOuterCryostatBlindFlangeColor);
  pOuterCryostatBlindFlangeColorVisAtt->SetVisibility(true);
  m_pOuterCryostatBlindFlangeLogicalVolume->SetVisAttributes(pOuterCryostatBlindFlangeColorVisAtt);

  G4Colour hOuterCryostatVacuumColor(1,1,1,OuterVesselAlphaChannel);  //1 1 1  white
  G4VisAttributes *pOuterCryostatVacuumColorVisAtt = new G4VisAttributes(hOuterCryostatVacuumColor);
  pOuterCryostatVacuumColorVisAtt->SetVisibility(true);
  m_pOuterCryostatVacuumLogicalVolume->SetVisAttributes(pOuterCryostatVacuumColorVisAtt);
 }

//******************************************************************/
// ConstructInnerCryostat
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructInnerCryostat() {
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< cryostat >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  const G4double dLXeHeight = GetGeometryParameter("LXeHeight");
  const G4double dLXeHalfZ = 0.5*dLXeHeight;
  //const G4double d = GetGeometryParameter("");

  const G4double dInnerCryostatVesselOuterRadius = GetGeometryParameter("InnerCryostatVesselOuterRadius");
  const G4double dInnerCryostatVesselInnerRadius = GetGeometryParameter("InnerCryostatVesselInnerRadius");
  const G4double dInnerCryostatVesselHeight = GetGeometryParameter("InnerCryostatVesselHeight");
  const G4double dOuterCryostatVesselHeight = GetGeometryParameter("OuterCryostatVesselHeight");
  //const G4double dInnerCryostatVesselThickness = GetGeometryParameter("InnerCryostatVesselThickness");

  //const G4double dCryostatToTopPMTCasingSpacing = GetGeometryParameter("CryostatToTopPMTCasingSpacing");
  //const G4double dCryostatBottomDistanceBetweenLXeAndPTFE = GetGeometryParameter("CryostatBottomDistanceBetweenLXeAndPTFE");

  const G4double dInnerCryostatBlindFlangeThickness = GetGeometryParameter("InnerCryostatBlindFlangeThickness");
  const G4double dInnerCryostatBlindFlangeOuterRadius= GetGeometryParameter("InnerCryostatBlindFlangeOuterRadius");

  const G4double dInnerCryostatBoreFlangeThickness = GetGeometryParameter("InnerCryostatBoreFlangeThickness");
  const G4double dInnerCryostatBoreFlangeOuterRadius= GetGeometryParameter("InnerCryostatBoreFlangeOuterRadius");

  //const G4double dTopBoreFlangeInnerCryostatToLiquidLevel= GetGeometryParameter("dTopBoreFlangeInnerCryostatToLiquidLevel");

  const G4double dInnerCryostatVesselHalfZ = dInnerCryostatVesselHeight*0.5;
  const G4double dOuterCryostatVesselHalfZ = dOuterCryostatVesselHeight*0.5;
  const G4double dTotalTPCHalfZ = dOuterCryostatVesselHeight*0.5;
  const G4double dOuterCryostatBottomToBottomFlangeTPCCryostat= GetGeometryParameter("OuterCryostatBottomToBottomFlangeTPCCryostat");

  G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");
  G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
  G4Material *Teflon = G4Material::GetMaterial("Teflon");
  G4Material *SS316LSteel = G4Material::GetMaterial("SS316LSteel");
  G4Material *Copper = G4Material::GetMaterial("Copper");
  G4Material *GXeTeflon = G4Material::GetMaterial("GXeTeflon");

  //=============================== cryostat vessel ===============================
  const G4double dInnerCryostatVesselOffsetZ = GetGeometryParameter("InnerCryostatVesselOffsetZ"); 
 
  //G4cout << "InnerCryostatVesselOffsetZ: " << dInnerCryostatVesselOffsetZ << G4endl;

  G4Tubs *pInnerCryostatVesselTubs = new G4Tubs("InnerCryostatVesselTubs", dInnerCryostatVesselInnerRadius,dInnerCryostatVesselOuterRadius,dInnerCryostatVesselHalfZ , 0.*deg, 360.*deg);

  m_pInnerCryostatVesselLogicalVolume = new G4LogicalVolume(pInnerCryostatVesselTubs, SS316LSteel, "InnerCryostatVesselLogicalVolume", 0, 0, 0);

  // m_pInnerCryostatVesselPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dInnerCryostatVesselOffsetZ), m_pInnerCryostatVesselLogicalVolume, "InnerCryostat", m_pLXeLogicalVolume, false, 0);
  //m_pInnerCryostatVesselPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dInnerCryostatVesselOffsetZ), m_pInnerCryostatVesselLogicalVolume, "InnerCryostat", m_pLabLogicalVolume, false, 0);
  m_pInnerCryostatVesselPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dInnerCryostatVesselOffsetZ), m_pInnerCryostatVesselLogicalVolume, "InnerCryostat", m_pOuterCryostatVacuumLogicalVolume, false, 0);

  // =========================== bore flanges  ====================================
  const G4double dInnerCryostatBoreFlangeHalfZ = dInnerCryostatBoreFlangeThickness*0.5;
  const G4int dBoreFlangeNumber = GetGeometryParameter("BoreFlangeNumber");
	
  const G4double dInnerCryostatBottomBoreFlangeOffsetZ = (dInnerCryostatVesselOffsetZ-dLXeHalfZ+dInnerCryostatBoreFlangeHalfZ);//-(dLXeHalfZ-dInnerCryostatBoreFlangeHalfZ);// -(dLXeHalfZ+dTopBoreFlangeInnerCryostatToLiquidLevel-dInnerCryostatBoreFlangeHalfZ);

  G4Tubs *pInnerCryostatBoreFlangeTubs = new G4Tubs("InnerCryostatBoreFlangeTubs", dInnerCryostatVesselOuterRadius,dInnerCryostatBoreFlangeOuterRadius,dInnerCryostatBoreFlangeHalfZ , 0.*deg, 360.*deg);

  m_pInnerCryostatBoreFlangeLogicalVolume = new G4LogicalVolume(pInnerCryostatBoreFlangeTubs, SS316LSteel, "InnerCryostatBoreFlangeLogicalVolume", 0, 0, 0);
       
  G4double pdBoreFlangeZ[dBoreFlangeNumber];
	
  for(int i = 0; i < dBoreFlangeNumber; i++)
    {
	  
      pdBoreFlangeZ[i]=dInnerCryostatBottomBoreFlangeOffsetZ+i*(dInnerCryostatVesselHeight-dInnerCryostatBoreFlangeThickness);//334-22.2;

      // m_hInnerCryostatBoreFlangePhysicalVolumes.push_back
      // 	(new G4PVPlacement(0,G4ThreeVector(0,0,pdBoreFlangeZ[i]),
      // 			   m_pInnerCryostatBoreFlangeLogicalVolume,"InnerCryostatBoreFlange",m_pLXeLogicalVolume, false, i));
      // m_hInnerCryostatBoreFlangePhysicalVolumes.push_back
      // 	(new G4PVPlacement(0,G4ThreeVector(0,0,pdBoreFlangeZ[i]),
      // 			   m_pInnerCryostatBoreFlangeLogicalVolume,"InnerCryostatBoreFlange",m_pLabLogicalVolume, false, i));

 m_hInnerCryostatBoreFlangePhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdBoreFlangeZ[i]),
			   m_pInnerCryostatBoreFlangeLogicalVolume,"InnerCryostatBoreFlange", m_pOuterCryostatVacuumLogicalVolume, false, i));

      //G4cout<<"pdBoreFlangeZ[i]: "<<pdBoreFlangeZ[i]<<G4endl;
                                                                                    
    }

  // ======================== blind flanges ==========================================
  const G4double dInnerCryostatBlindFlangeHalfZ = dInnerCryostatBlindFlangeThickness*0.5;
  const G4int dBlindFlangeNumber = GetGeometryParameter("BlindFlangeNumber");
	
  const G4double dInnerCryostatBottomBlindFlangeOffsetZ = (dInnerCryostatVesselOffsetZ-dLXeHalfZ-dInnerCryostatBoreFlangeThickness+dInnerCryostatBlindFlangeHalfZ);//(dLXeHalfZ+dInnerCryostatBlindFlangeHalfZ);// -(dLXeHalfZ+dTopBlindFlangeInnerCryostatToLiquidLevel-dInnerCryostatBlindFlangeHalfZ);

  G4Tubs *pInnerCryostatBlindFlangeTubs = new G4Tubs("InnerCryostatBlindFlangeTubs", 0,dInnerCryostatBlindFlangeOuterRadius,dInnerCryostatBlindFlangeHalfZ , 0.*deg, 360.*deg);

  m_pInnerCryostatBlindFlangeLogicalVolume = new G4LogicalVolume(pInnerCryostatBlindFlangeTubs, SS316LSteel, "InnerCryostatBlindFlangeLogicalVolume", 0, 0, 0);
       
  G4double pdBlindFlangeZ[dBlindFlangeNumber];
	
  for(int i = 0; i < dBlindFlangeNumber; i++)
    {
	  
      pdBlindFlangeZ[i]=dInnerCryostatBottomBlindFlangeOffsetZ+i*(dInnerCryostatVesselHeight+dInnerCryostatBlindFlangeThickness);//334-22.2;

      // m_hInnerCryostatBlindFlangePhysicalVolumes.push_back
      // 	(new G4PVPlacement(0,G4ThreeVector(0,0,pdBlindFlangeZ[i]),
      // 			   m_pInnerCryostatBlindFlangeLogicalVolume,"InnerCryostatBlindFlange",m_pLXeLogicalVolume, false, i));

      // m_hInnerCryostatBlindFlangePhysicalVolumes.push_back
      // 	(new G4PVPlacement(0,G4ThreeVector(0,0,pdBlindFlangeZ[i]),
      // 				   m_pInnerCryostatBlindFlangeLogicalVolume,"InnerCryostatBlindFlange",m_pLabLogicalVolume, false, i));
      m_hInnerCryostatBlindFlangePhysicalVolumes.push_back
       	(new G4PVPlacement(0,G4ThreeVector(0,0,pdBlindFlangeZ[i]),
       				   m_pInnerCryostatBlindFlangeLogicalVolume,"InnerCryostatBlindFlange",m_pOuterCryostatVacuumLogicalVolume, false, i));
      //G4cout<<"pdBlindFlangeZ[i]: "<<pdBlindFlangeZ[i]<<G4endl;                                                                    
    }

  // ========================== copper rings =============================
  //the 2 copper rings between the blind and bore flanges

  const G4double dInnerCryostatCopperRingsOuterRadius = GetGeometryParameter("InnerCryostatCopperRingsOuterRadius");
  const G4double dInnerCryostatCopperRingsInnerRadius = GetGeometryParameter("InnerCryostatCopperRingsInnerRadius");
  const G4double dInnerCryostatCopperRingsThickness = GetGeometryParameter("InnerCryostatCopperRingsThickness");
  const G4int dInnerCryostatCopperRingsNumber= GetGeometryParameter("InnerCryostatCopperRingsNumber");

  const G4double dInnerCryostatCopperRingsHalfZ = 0.5*GetGeometryParameter("InnerCryostatCopperRingsThickness");

  const G4double dInnerCryostatBottomCopperRingsOffsetZ =(dInnerCryostatVesselOffsetZ-dLXeHalfZ+dInnerCryostatCopperRingsHalfZ);// -(dLXeHalfZ+dInnerCryostatCopperRingsHalfZ);// -(dLXeHalfZ+dTopCopperRingsInnerCryostatToLiquidLevel-dInnerCryostatCopperRingsHalfZ);
  
  G4Tubs *pInnerCryostatCopperRingsTubs = new G4Tubs("InnerCryostatCopperRingsTubs", dInnerCryostatCopperRingsInnerRadius,dInnerCryostatCopperRingsOuterRadius,dInnerCryostatCopperRingsHalfZ , 0.*deg, 360.*deg);
  
     m_pInnerCryostatCopperRingsLogicalVolume = new G4LogicalVolume(pInnerCryostatCopperRingsTubs, Copper, "InnerCryostatCopperRingsLogicalVolume", 0, 0, 0);
       
     G4double pdInnerCryostatCopperRingsZ[dInnerCryostatCopperRingsNumber];
	
  for(int i = 0; i < dInnerCryostatCopperRingsNumber; i++)
    {
	  
      pdInnerCryostatCopperRingsZ[i]=dInnerCryostatBottomCopperRingsOffsetZ+i*(dInnerCryostatVesselHeight-dInnerCryostatCopperRingsThickness);

      // m_hInnerCryostatCopperRingsPhysicalVolumes.push_back
      // 	(new G4PVPlacement(0,G4ThreeVector(0,0,pdInnerCryostatCopperRingsZ[i]),
      // 			   m_pInnerCryostatCopperRingsLogicalVolume,"InnerCryostatCopperRings",m_pLabLogicalVolume, false, i));
  m_hInnerCryostatCopperRingsPhysicalVolumes.push_back
      	(new G4PVPlacement(0,G4ThreeVector(0,0,pdInnerCryostatCopperRingsZ[i]),
      			   m_pInnerCryostatCopperRingsLogicalVolume,"InnerCryostatCopperRings",m_pOuterCryostatVacuumLogicalVolume, false, i));
      //G4cout<<"pdInnerCryostatCopperRingsZ[i]: "<<pdInnerCryostatCopperRingsZ[i]<<G4endl;
                                                                                 
    }

  //================================== attributes =================================
  G4Colour hStainlessSteelColor(0.500, 0.500, 0.500, 0.1);
  G4Colour hVacuumColor(0.600, 0.200, 0.250, 0.1);
  G4Colour hRedColor(1.,0., 0., 0.5);
  G4VisAttributes *pRedVisAtt = new G4VisAttributes(hRedColor);
  pRedVisAtt->SetVisibility(true);

  G4Colour hInnerCryostatVesselColor(0.012,0.705,0.8,InnerVesselAlphaChannel);  // 3 180 204 true iris blu
  G4VisAttributes *pInnerCryostatVesselColorVisAtt = new G4VisAttributes(hInnerCryostatVesselColor);
  pInnerCryostatVesselColorVisAtt->SetVisibility(true);
  m_pInnerCryostatVesselLogicalVolume->SetVisAttributes(pInnerCryostatVesselColorVisAtt);

  G4Colour hInnerCryostatBoreFlangeColor(1.0,0.549, 0,InnerFlangeAlphaChannel);  //255 140 0 orange
  G4VisAttributes *pInnerCryostatBoreFlangeColorVisAtt = new G4VisAttributes(hInnerCryostatBoreFlangeColor);
  pInnerCryostatBoreFlangeColorVisAtt->SetVisibility(true);
  m_pInnerCryostatBoreFlangeLogicalVolume->SetVisAttributes(pInnerCryostatBoreFlangeColorVisAtt);

  G4Colour hInnerCryostatBlindFlangeColor(0.545,0.353, 0,InnerFlangeAlphaChannel);  //139 90 0 dark orange
  G4VisAttributes *pInnerCryostatBlindFlangeColorVisAtt = new G4VisAttributes(hInnerCryostatBlindFlangeColor);
  pInnerCryostatBlindFlangeColorVisAtt->SetVisibility(true);
  m_pInnerCryostatBlindFlangeLogicalVolume->SetVisAttributes(pInnerCryostatBlindFlangeColorVisAtt);

  G4Colour hInnerCryostatCopperRingsColor(0.635,0.804, 0.353,CopperRingsAlphaChannel);  //162 205 90  dark olive green
  G4VisAttributes *pInnerCryostatCopperRingsColorVisAtt = new G4VisAttributes(hInnerCryostatCopperRingsColor);
  pInnerCryostatCopperRingsColorVisAtt->SetVisibility(true);
  m_pInnerCryostatCopperRingsLogicalVolume->SetVisAttributes(pInnerCryostatCopperRingsColorVisAtt);
}

//******************************************************************/
// ConstructXenon
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructXenon() {
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< xenon >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  const G4double dLXeRadius = GetGeometryParameter("InnerCryostatVesselOuterRadius");

  //G4cout << "LXe radius: " << dLXeRadius << G4endl;

  const G4double dLXeHeight = GetGeometryParameter("LXeHeight");
  const G4double dGXeHeight = GetGeometryParameter("GXeHeight");

  //G4cout << "dLXeHeight: " << dLXeHeight << G4endl;
  //G4cout << "dGXeHeight: " << dGXeHeight << G4endl;

  // const G4double dZeroOffset = (m_pMotherLogicalVolume==m_pLXeLogicalVolume)?(GetGeometryParameter("GateMeshZeroOffset")):(0.);
  // G4cout<<"dZeroOffset: " << dZeroOffset << G4endl;

  const G4double dLiquidLevelToGateMesh = GetGeometryParameter("LiquidLevelToGateMesh");
  const G4double dInnerCryostatVesselOffsetZ = GetGeometryParameter("InnerCryostatVesselOffsetZ"); 

  //G4Material *LXe = G4Material::GetMaterial("LXe");
	G4Material *LXe = LXeMaterial;
  G4Material *GXe = G4Material::GetMaterial("GXe");

  //================================ liquid xenon =================================
  //---------------------- tube the size of the cryostat -----------------------
  const G4double dLXeHalfZ = 0.5*dLXeHeight;
  //G4cout<<"dLXeHalfZ: " << dLXeHalfZ << G4endl;

    const G4double dLXeOffsetZ =dInnerCryostatVesselOffsetZ;
  //G4cout<<"LXeOffsetZ: "<<dLXeOffsetZ<<G4endl;

  G4Tubs *pLXeTubs = new G4Tubs("LXeTubs", 0.*cm, dLXeRadius, dLXeHalfZ, 0.*deg, 360.*deg);

  m_pLXeLogicalVolume = new G4LogicalVolume(pLXeTubs, LXe, "LXeVolume", 0, 0, 0);

  // m_pLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLXeOffsetZ),
  // 					   m_pLXeLogicalVolume, "LXe", m_pMotherLogicalVolume, false, 0);
  m_pLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLXeOffsetZ),
					   m_pLXeLogicalVolume, LXeMaterial->GetName(), m_pOuterCryostatVacuumLogicalVolume, false, 0);

	//TackLimit for opticalphotons (see physics list)
  G4double maxTrack = 100000;
  fStepLimit = new G4UserLimits(DBL_MAX,maxTrack);
  m_pLXeLogicalVolume->SetUserLimits(fStepLimit);
	
  //================================ gaseous xenon ================================
  //defined as a G4Tubs inside the LXe G4Tubs
		
  const G4double dGXeRadius = dLXeRadius;//GetGeometryParameter("InnerCryostatVesselInnerRadius");
  const G4double dGXeHalfZ = 0.5*dGXeHeight;
  const G4double dGXeOffsetZ = dLXeHalfZ-dGXeHalfZ;

  //G4cout<<"dGXeHalfZ: " << dGXeHalfZ << G4endl;
  //G4cout<<"GXeOffsetZ: "<< dGXeOffsetZ <<G4endl;

  G4Tubs *pGXeTubs = new G4Tubs("GXeTube", 0.*cm, dGXeRadius, dGXeHalfZ, 0.*deg, 360.*deg);

  m_pGXeLogicalVolume = new G4LogicalVolume(pGXeTubs, GXe, "GXeLogicalVolume", 0, 0, 0);

  m_pGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dGXeOffsetZ),
					   m_pGXeLogicalVolume, "GXe", m_pLXeLogicalVolume, false, 0);


  //------------------------------ xenon sensitivity ------------------------------
  G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

  muensterTPCLXeSensitiveDetector *pLXeSD = new muensterTPCLXeSensitiveDetector("muensterTPC/LXeSD");
  pSDManager->AddNewDetector(pLXeSD);
  m_pLXeLogicalVolume->SetSensitiveDetector(pLXeSD);
  m_pGXeLogicalVolume->SetSensitiveDetector(pLXeSD);

  //================================== attributes =================================
  G4Colour hLXeColor(0.0,0.0,1.0,DetectorMaterialAlphaChannel); //blue
  G4Colour hGXeColor(0.0,1.0,1.0,DetectorMaterialAlphaChannel); //cyan

  G4VisAttributes *pLXeVisAtt = new G4VisAttributes(hLXeColor);
  pLXeVisAtt->SetVisibility(true);
  m_pLXeLogicalVolume->SetVisAttributes(pLXeVisAtt);
	
  G4VisAttributes *pGXeVisAtt = new G4VisAttributes(hGXeColor);
  pGXeVisAtt->SetVisibility(true);
  m_pGXeLogicalVolume->SetVisAttributes(pGXeVisAtt);
}

//******************************************************************/
// ConstructLiquidXenonCylinder
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructLiquidXenonCylinder() {
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< xenon >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//nEXO-like detector
   
  const G4double dLXeRadius = 60.*cm;

  G4cout << "LXe radius: " << dLXeRadius << G4endl;

  const G4double dLXeHeight = 120.*cm;

    
//DARWIN-like detector
  /* 
  const G4double dLXeRadius = 125.*cm;

  G4cout << "LXe radius: " << dLXeRadius << G4endl;

  const G4double dLXeHeight = 250.*cm;
  */
    
  G4Material *LXe = LXeMaterial;
  G4Material *GXe = G4Material::GetMaterial("GXe");
  G4Material *Copper = G4Material::GetMaterial("Copper");

  //================================ liquid xenon =================================
  //---------------------- tube the size of the cryostat -----------------------
  const G4double dLXeHalfZ = 0.5*dLXeHeight;
  //G4cout<<"dLXeHalfZ: " << dLXeHalfZ << G4endl;

    const G4double dLXeOffsetZ = 0.*cm;
  //G4cout<<"LXeOffsetZ: "<<dLXeOffsetZ<<G4endl;

  G4Tubs *pLXeTubs = new G4Tubs("LXeTubs", 0.*cm, dLXeRadius, dLXeHalfZ, 0.*deg, 360.*deg);

  m_pLXeLogicalVolume = new G4LogicalVolume(pLXeTubs, LXe, "LXeVolume", 0, 0, 0);

  // m_pLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLXeOffsetZ),
  // 					   m_pLXeLogicalVolume, "LXe", m_pMotherLogicalVolume, false, 0);
  m_pLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLXeOffsetZ),
					   m_pLXeLogicalVolume, LXeMaterial->GetName(), m_pLabLogicalVolume, false, 0);

	//TackLimit for opticalphotons (see physics list)
  G4double maxTrack = 100000;
  fStepLimit = new G4UserLimits(DBL_MAX,maxTrack);
  m_pLXeLogicalVolume->SetUserLimits(fStepLimit);
	
  //------------------------------ xenon sensitivity ------------------------------
  G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

  muensterTPCLXeSensitiveDetector *pLXeSD = new muensterTPCLXeSensitiveDetector("muensterTPC/LXeSD");
  pSDManager->AddNewDetector(pLXeSD);
  m_pLXeLogicalVolume->SetSensitiveDetector(pLXeSD);

  //================================== attributes =================================
  G4Colour hLXeColor(0.0,0.0,1.0,DetectorMaterialAlphaChannel); //blue

  G4VisAttributes *pLXeVisAtt = new G4VisAttributes(hLXeColor);
  pLXeVisAtt->SetVisibility(true);
  m_pLXeLogicalVolume->SetVisAttributes(pLXeVisAtt);

  G4cout << "Constructed cylinder of LXe..." << G4endl;


  //=================== Outer shell around LXe volume ============================
  G4Tubs * outerShellTubs = new G4Tubs("OuterShellTubs",dLXeRadius,dLXeRadius + 1.0*cm, 
                                       dLXeHalfZ + 1.0*cm, 0.*deg, 360.*deg);
  m_pOuterShellLogicalVolume = new G4LogicalVolume(outerShellTubs, Copper, "OuterShellLog",0,0,0);

  m_pOuterShellPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dLXeOffsetZ),
					   m_pOuterShellLogicalVolume, "CopperShell", m_pLabLogicalVolume, false, 0);


	
}

//******************************************************************/
// ConstructTPC
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructTPC() {
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TPC >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  const G4double dLXeLevelZeroOffset =GetGeometryParameter("LXeLevelZeroOffset");

  const G4double dLXeHeight = GetGeometryParameter("LXeHeight");
  const G4double dLXeHalfZ = 0.5*dLXeHeight;
  const G4double dLiquidLevelToGateMesh = GetGeometryParameter("LiquidLevelToGateMesh");

  const G4double dGXeHeight = GetGeometryParameter("GXeHeight");
  const G4double dGXeHalfZ = 0.5*dGXeHeight;

  const G4double dTeflonCentralCylinderInnerRadius= GetGeometryParameter("TeflonCentralCylinderInnerRadius");
  const G4double dTeflonCentralCylinderOuterRadius= GetGeometryParameter("TeflonCentralCylinderOuterRadius");
  const G4double dTeflonCentralCylinderHeight =  GetGeometryParameter("TeflonCentralCylinderHeight");
  const G4double dTeflonCentralCylinderHalfZ =  0.5*dTeflonCentralCylinderHeight;
  
  //G4cout << "TeflonCentralCylinderHeight: " << dTeflonCentralCylinderHeight << G4endl;
  //G4cout << "TeflonCentralCylinderInnerRadius: " << dTeflonCentralCylinderInnerRadius << G4endl;
  //G4cout << "TeflonCentralCylinderOuterRadius: " << dTeflonCentralCylinderOuterRadius << G4endl;

  const G4double dTeflonOuterCylinderInnerRadius= GetGeometryParameter("TeflonOuterCylinderInnerRadius");
  const G4double dTeflonOuterCylinderOuterRadius= GetGeometryParameter("TeflonOuterCylinderOuterRadius");
  const G4double dTeflonOuterCylinderHeight =  GetGeometryParameter("TeflonOuterCylinderHeight");
  const G4double dTeflonOuterCylinderHalfZ =  0.5*dTeflonOuterCylinderHeight;
  
  //G4cout << "TeflonOuterCylinderHeight: " << dTeflonOuterCylinderHeight << G4endl;
  //G4cout << "TeflonOuterCylinderInnerRadius: " << dTeflonOuterCylinderInnerRadius << G4endl;
  //G4cout << "TeflonOuterCylinderOuterRadius: " << dTeflonOuterCylinderOuterRadius << G4endl;

  const G4double dCryostatBottomDistanceBetweenLXeAndPTFE= GetGeometryParameter("CryostatBottomDistanceBetweenLXeAndPTFE");

  G4Material *SS316LSteel = G4Material::GetMaterial("SS316LSteel");
  G4Material *PTFE = G4Material::GetMaterial("Teflon");
  G4Material *GXeTeflon = G4Material::GetMaterial("GXeTeflon");

  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< PTFE Inner cylinder >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  const G4double dPTFEInnerCylinderOffsetZ = dLXeLevelZeroOffset-0.5*dTeflonCentralCylinderHeight-dLiquidLevelToGateMesh;
  //G4cout << "PTFEInnerCylinderOffsetZ: " << dPTFEInnerCylinderOffsetZ << G4endl;

  G4Tubs *pPTFEInnerCylinderTubs = new G4Tubs("PTFEInnerCylinderTubs", dTeflonCentralCylinderInnerRadius, dTeflonCentralCylinderOuterRadius, dTeflonCentralCylinderHalfZ, 0.*deg, 360.*deg);

  m_pPTFEInnerCylinderLogicalVolume = new G4LogicalVolume(pPTFEInnerCylinderTubs, PTFE, "PTFEInnerCylinderLogicalVolume", 0, 0, 0);

  m_pPTFEInnerCylinderPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPTFEInnerCylinderOffsetZ),
							 m_pPTFEInnerCylinderLogicalVolume, "PTFEInnerCylinder", m_pLXeLogicalVolume, false, 0);


  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< PTFE Outer cylinder >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // In reality these are slabs put on top of one another but it can be approximated very well to a cylinder
  const G4double dPTFEOuterCylinderOffsetZ = dLXeLevelZeroOffset-0.5*dTeflonOuterCylinderHeight-dLiquidLevelToGateMesh;
  //G4cout << "PTFEOuterCylinderOffsetZ: " << dPTFEOuterCylinderOffsetZ << G4endl;

  G4Tubs *pPTFEOuterCylinderTubs = new G4Tubs("PTFEOuterCylinderTubs", dTeflonOuterCylinderInnerRadius, dTeflonOuterCylinderOuterRadius, dTeflonOuterCylinderHalfZ, 0.*deg, 360.*deg);

  m_pPTFEOuterCylinderLogicalVolume = new G4LogicalVolume(pPTFEOuterCylinderTubs, PTFE, "PTFEOuterCylinderLogicalVolume", 0, 0, 0);

  m_pPTFEOuterCylinderPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPTFEOuterCylinderOffsetZ),
  							 m_pPTFEOuterCylinderLogicalVolume, "PTFEOuterCylinder", m_pLXeLogicalVolume, false, 0);

  
  //==========================  PTFE bottom PMT holder ==========================
  //includes in reality 3 slabs. top part of pmt holder, bottom part of pmt holder, pmt base holder
  //const G4double = GetGeometryParameter("");
  const G4double dBottomPMTHolderTopPartThickness = GetGeometryParameter("BottomPMTHolderTopPartThickness");
  const G4double dBottomPMTHolderBottomPartThickness= GetGeometryParameter("BottomPMTHolderBottomPartThickness");
  const G4double dBottomPMTBaseHolderThickness= GetGeometryParameter("BottomPMTBaseHolderThickness");
  
  const G4double dPTFEBottomPMTHolderHalfZ = 0.5*(dBottomPMTHolderTopPartThickness+dBottomPMTHolderBottomPartThickness+dBottomPMTBaseHolderThickness);
  //G4cout << "PTFEBottomPMTHolderHalfZ: " << dPTFEBottomPMTHolderHalfZ << G4endl;

  const G4double dPTFEBottomPMTHolderOffsetZ = -(dLXeHalfZ-dPTFEBottomPMTHolderHalfZ-dCryostatBottomDistanceBetweenLXeAndPTFE);
  //G4cout << "PTFEBottomPMTHolderOffsetZ: " << dPTFEBottomPMTHolderOffsetZ << G4endl;

  G4Tubs *pPTFEBottomPMTHolderTubs = new G4Tubs("PTFEBottomPMTHolderTubs", 0, dTeflonOuterCylinderOuterRadius, dPTFEBottomPMTHolderHalfZ, 0.*deg, 360.*deg);

  m_pPTFEBottomPMTHolderLogicalVolume = new G4LogicalVolume(pPTFEBottomPMTHolderTubs, PTFE, "PTFEBottomPMTHolderLogicalVolume", 0, 0, 0);

  m_pPTFEBottomPMTHolderPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPTFEBottomPMTHolderOffsetZ),
							   m_pPTFEBottomPMTHolderLogicalVolume, "PTFEBottomPMTHolder", m_pLXeLogicalVolume, false, 0);

  
  //=========================== PTFE slab between bottom pmt and cathode ==========================
  const G4double dBottomPTFESlabThickness= GetGeometryParameter("BottomPTFESlabThickness");
  const G4double dBottomPTFESlabInnerRadius = GetGeometryParameter("BottomPTFESlabInnerRadius");
  const G4double dBottomPTFESlabOuterRadius = GetGeometryParameter("BottomPTFESlabOuterRadius");

  const G4double dBottomPTFESlabHalfZ = dBottomPTFESlabThickness*0.5;
  //G4cout << "BottomPTFESlabHalfZ: " << dBottomPTFESlabHalfZ << G4endl;

  const G4double dBottomPTFESlabOffsetZ= -(dLXeHalfZ-2*dPTFEBottomPMTHolderHalfZ-dBottomPTFESlabHalfZ-dCryostatBottomDistanceBetweenLXeAndPTFE);
  //G4cout << "BottomPTFESlabOffsetZ: " <<  dBottomPTFESlabOffsetZ << G4endl;

  G4Tubs *pBottomPTFESlabTubs = new G4Tubs("BottomPTFESlabTubs", dBottomPTFESlabInnerRadius,dBottomPTFESlabOuterRadius, dBottomPTFESlabHalfZ, 0.*deg, 360.*deg);

  m_pBottomPTFESlabLogicalVolume = new G4LogicalVolume(pBottomPTFESlabTubs, PTFE, "BottomPTFESlabLogicalVolume", 0, 0, 0);

  m_pBottomPTFESlabPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBottomPTFESlabOffsetZ),
						      m_pBottomPTFESlabLogicalVolume, "BottomPTFESlab", m_pLXeLogicalVolume, false, 0);
  

  //==================PTFE slab between gate mesh and top pmt ========================
  // in reality two small slabs
  const G4double dTopPTFESlabThickness= GetGeometryParameter("TopPTFESlabThickness");
  const G4double dTopPTFESlabInnerRadius = GetGeometryParameter("TopPTFESlabInnerRadius");
  const G4double dTopPTFESlabOuterRadius = GetGeometryParameter("TopPTFESlabOuterRadius");

  const G4double dTopPTFESlabHalfZ = dTopPTFESlabThickness*0.5;
  //G4cout << "TopPTFESlabHalfZ: " << dTopPTFESlabHalfZ << G4endl;

  //const G4double dTopPTFESlabOffsetZ= dLXeLevelZeroOffset+dTopPTFESlabHalfZ-dLiquidLevelToGateMesh;//wrt LXe
  const G4double dTopPTFESlabOffsetZ= -dGXeHalfZ+dTopPTFESlabHalfZ-dLiquidLevelToGateMesh;//wrt GXe
  //G4cout << "TopPTFESlabOffsetZ: " <<  dTopPTFESlabOffsetZ << G4endl;

  G4Tubs *pTopPTFESlabTubs = new G4Tubs("TopPTFESlabTubs", dTopPTFESlabInnerRadius,dTopPTFESlabOuterRadius, dTopPTFESlabHalfZ, 0.*deg, 360.*deg);

  m_pTopPTFESlabLogicalVolume = new G4LogicalVolume(pTopPTFESlabTubs, GXeTeflon, "TopPTFESlabLogicalVolume", 0, 0, 0);

  //m_pTopPTFESlabPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTopPTFESlabOffsetZ),
  // m_pTopPTFESlabLogicalVolume, "TopPTFESlab", m_pLXeLogicalVolume, false, 0);

  m_pTopPTFESlabPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dTopPTFESlabOffsetZ), m_pTopPTFESlabLogicalVolume, "TopPTFESlab", m_pGXeLogicalVolume, false, 0);


  //================ PTFE Top PMT Holder====================
  // includes two slabs to hold pmt
  const G4double dTopPMTHolderTopPartThickness = GetGeometryParameter("TopPMTHolderTopPartThickness");
  const G4double dTopPMTHolderBottomPartThickness= GetGeometryParameter("TopPMTHolderBottomPartThickness");
  const G4double  dTopPMTHolderPlateRadius = GetGeometryParameter("TopPMTHolderPlateRadius");
  const G4double  dLiquidLevelToScreenMesh= GetGeometryParameter("LiquidLevelToScreenMesh");

  const G4double dPTFETopPMTHolderHalfZ = 0.5*(dTopPMTHolderTopPartThickness+dTopPMTHolderBottomPartThickness);
  //G4cout << "PTFETopPMTHolderHalfZ: " << dPTFETopPMTHolderHalfZ << G4endl;

  // const G4double dPTFETopPMTHolderOffsetZ = dLXeLevelZeroOffset-dLiquidLevelToGateMesh+dTopPTFESlabThickness+dPTFETopPMTHolderHalfZ;//(dLXeHalfZ+dPTFETopPMTHolderHalfZ+dLiquidLevelToScreenMesh);
  const G4double dPTFETopPMTHolderOffsetZ = -dGXeHalfZ-dLiquidLevelToGateMesh+dTopPTFESlabThickness+dPTFETopPMTHolderHalfZ; //wrt GXe
  //G4cout << "PTFETopPMTHolderOffsetZ: " << dPTFETopPMTHolderOffsetZ << G4endl;

  G4Tubs *pPTFETopPMTHolderTubs = new G4Tubs("PTFETopPMTHolderTubs", 0, dTopPMTHolderPlateRadius, dPTFETopPMTHolderHalfZ, 0.*deg, 360.*deg);

  m_pPTFETopPMTHolderLogicalVolume = new G4LogicalVolume(pPTFETopPMTHolderTubs, GXeTeflon, "PTFETopPMTHolderLogicalVolume", 0, 0, 0);

  // m_pPTFETopPMTHolderPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPTFETopPMTHolderOffsetZ),
  // 							m_pPTFETopPMTHolderLogicalVolume, "PTFETopPMTHolder", m_pLXeLogicalVolume, false, 0);
  m_pPTFETopPMTHolderPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dPTFETopPMTHolderOffsetZ),
   							m_pPTFETopPMTHolderLogicalVolume, "PTFETopPMTHolder", m_pGXeLogicalVolume, false, 0);


 //================ Steel ring at the bottom on which the PTFE slabs rest ===================
  const G4double dBottomSteelRingThickness=GetGeometryParameter("BottomSteelRingThickness");
	const G4double dBottomSteelRingInnerRadius=GetGeometryParameter("BottomSteelRingInnerRadius");
  const G4double dBottomSteelRingOuterRadius=GetGeometryParameter("BottomSteelRingOuterRadius");
  const G4double dBottomFlangeCopperRingToTopSteelRing=GetGeometryParameter("BottomFlangeCopperRingToTopSteelRing");

  const G4double dBottomSteelRingHalfZ = 0.5*dBottomSteelRingThickness;
  //G4cout << "BottomSteelRingHalfZ: " << dBottomSteelRingHalfZ << G4endl;

  const G4double dBottomSteelRingOffsetZ = -(dLXeHalfZ-dBottomFlangeCopperRingToTopSteelRing+dBottomSteelRingHalfZ);
  //G4cout << "BottomSteelRingOffsetZ: " << dBottomSteelRingOffsetZ << G4endl;

  G4Tubs *pBottomSteelRingTubs = new G4Tubs("BottomSteelRingTubs", dBottomSteelRingInnerRadius,dBottomSteelRingOuterRadius, dBottomSteelRingHalfZ, 0.*deg, 360.*deg);

  m_pBottomSteelRingLogicalVolume = new G4LogicalVolume(pBottomSteelRingTubs, SS316LSteel, "BottomSteelRingLogicalVolume", 0, 0, 0);

  m_pBottomSteelRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dBottomSteelRingOffsetZ), m_pBottomSteelRingLogicalVolume, "SS316LSteelBottomSteelRing", m_pLXeLogicalVolume, false, 0);

  //=============================== optical surfaces ==============================
  G4double dSigmaAlpha = 0.1;
	G4OpticalSurface *pTeflonOpticalSurface = new G4OpticalSurface("TeflonOpticalSurface",
		unified, ground, dielectric_metal, dSigmaAlpha);
	
	G4OpticalSurface *pSS316LSteelOpticalSurface = new G4OpticalSurface("SS316LSteelOpticalSurface",
		unified, polished, dielectric_metal, 0.);
		
	G4OpticalSurface *pGXeTeflonOpticalSurface = new G4OpticalSurface("GXeTeflonOpticalSurface", 
		unified, groundbackpainted, dielectric_dielectric, dSigmaAlpha);
		
	pTeflonOpticalSurface->SetMaterialPropertiesTable(PTFE->GetMaterialPropertiesTable());
	pSS316LSteelOpticalSurface->SetMaterialPropertiesTable(SS316LSteel->GetMaterialPropertiesTable());
	pGXeTeflonOpticalSurface->SetMaterialPropertiesTable(GXeTeflon->GetMaterialPropertiesTable());
	
	new G4LogicalBorderSurface("LXeTeflonCylinderLogicalBorderSurface",
		m_pLXePhysicalVolume, m_pPTFEInnerCylinderPhysicalVolume, pTeflonOpticalSurface);

	new G4LogicalBorderSurface("LXeTeflonSlabLogicalBorderSurface",
		m_pLXePhysicalVolume, m_pBottomPTFESlabPhysicalVolume, pTeflonOpticalSurface);

	new G4LogicalBorderSurface("LXeTeflonPMTHolderLogicalBorderSurface",
		m_pLXePhysicalVolume, m_pPTFEBottomPMTHolderPhysicalVolume, pTeflonOpticalSurface);

	new G4LogicalBorderSurface("GXeTeflonPMTHolderLogicalBorderSurface",
		m_pGXePhysicalVolume, m_pPTFETopPMTHolderPhysicalVolume, pGXeTeflonOpticalSurface);

	new G4LogicalBorderSurface("GXeTeflonSlabLogicalBorderSurface",
		m_pGXePhysicalVolume, m_pTopPTFESlabPhysicalVolume, pGXeTeflonOpticalSurface);
	
  //================================== attributes =================================
  G4Colour hCopperColor(0.835, 0.424, 0.059, CopperRingsAlphaChannel);
  G4Colour hTeflonColor(1., 0., 1.,PTFECylinderAlphaChannel); //magenta
  G4Colour hOutTeflon(0.545, 0.227, 0.384,PTFECylinderAlphaChannel); //139 58 98 hotpink
  G4Colour hRedColor(1.,0., 0., PMTHolderAlphaChannel); //
  G4Colour hYellowColor(1.,1., 0.,0.25); //yello
  G4Colour hSteelColor(0.87,0.72, 0.53,0.25); //222 184 135 burlywood
  // G4VisAttributes *pRedVisAtt = new G4VisAttributes(hRedColor);
  // pRedVisAtt->SetVisibility(true);

  G4VisAttributes *pPTFEInnerCylinderVisAtt = new G4VisAttributes(hTeflonColor);
  pPTFEInnerCylinderVisAtt->SetVisibility(true);
  m_pPTFEInnerCylinderLogicalVolume->SetVisAttributes(pPTFEInnerCylinderVisAtt);

  G4VisAttributes *pPTFEOuterCylinderVisAtt = new G4VisAttributes(hTeflonColor);
  pPTFEOuterCylinderVisAtt->SetVisibility(true);
  m_pPTFEOuterCylinderLogicalVolume->SetVisAttributes(pPTFEOuterCylinderVisAtt);

  G4VisAttributes *pPTFEBottomPMTHolderVisAtt = new G4VisAttributes(hRedColor);
  pPTFEBottomPMTHolderVisAtt->SetVisibility(true);
  m_pPTFEBottomPMTHolderLogicalVolume->SetVisAttributes(pPTFEBottomPMTHolderVisAtt);

  G4VisAttributes *pBottomPTFESlabVisAtt = new G4VisAttributes(hYellowColor);
  pBottomPTFESlabVisAtt->SetVisibility(true);
  m_pBottomPTFESlabLogicalVolume->SetVisAttributes(pBottomPTFESlabVisAtt);

  G4VisAttributes *pPTFETopPMTHolderVisAtt = new G4VisAttributes(hRedColor);
  pPTFETopPMTHolderVisAtt->SetVisibility(true);
  m_pPTFETopPMTHolderLogicalVolume->SetVisAttributes(pPTFETopPMTHolderVisAtt);
  
  G4VisAttributes *pTopPTFESlabVisAtt = new G4VisAttributes(hYellowColor);
  pTopPTFESlabVisAtt->SetVisibility(true);
  m_pTopPTFESlabLogicalVolume->SetVisAttributes(pTopPTFESlabVisAtt);

  G4VisAttributes *pBottomSteelRingVisAtt = new G4VisAttributes(hSteelColor);
  pBottomSteelRingVisAtt->SetVisibility(true);
  m_pBottomSteelRingLogicalVolume->SetVisAttributes(pBottomSteelRingVisAtt);
}

//******************************************************************/
// ConstructFieldCage
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructFieldCage() {
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< field cage >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  const G4double dLXeHeight = GetGeometryParameter("LXeHeight");
  const G4double dLXeHalfZ = 0.5*dLXeHeight;
  const G4double dLiquidLevelToGateMesh = GetGeometryParameter("LiquidLevelToGateMesh");
 
  G4Material *SS316LSteel = G4Material::GetMaterial("SS316LSteel");
  G4Material *Teflon = G4Material::GetMaterial("Teflon");
  G4Material *Copper = G4Material::GetMaterial("Copper");
  G4Material *GridMeshSS316LSteelLXe = G4Material::GetMaterial("GridMeshSS316LSteelLXe");
  G4Material *GridMeshSS316LSteelGXe = G4Material::GetMaterial("GridMeshSS316LSteelGXe");

  const G4double dLXeLevelZeroOffset = GetGeometryParameter("LXeLevelZeroOffset");

  //================================= CopperRings ===============================
  const G4double dCopperRingsOuterRadius= GetGeometryParameter("CopperRingsOuterRadius");
  const G4double dCopperRingsInnerRadius= GetGeometryParameter("CopperRingsInnerRadius");
  const G4double dCopperRingsThickness= GetGeometryParameter("CopperRingsThickness");
  const G4int dCopperRingsNumber= GetGeometryParameter("CopperRingsNumber");
  const G4double dCopperRingsSpacing= GetGeometryParameter("CopperRingsSpacing");
  const G4double dGateMeshRingThickness= GetGeometryParameter("GateMeshRingThickness");

  const G4double dFirstCopperRingOffsetZ = dLXeLevelZeroOffset-dLiquidLevelToGateMesh-dGateMeshRingThickness-dCopperRingsSpacing-0.5*dCopperRingsThickness;//(dLXeHalfZ-dLiquidLevelToGateMesh-dGateMeshRingThickness-dCopperRingsSpacing-0.5*dCopperRingsThickness);

  G4Tubs *pCopperRingsTubs = new G4Tubs("CopperRingsTube", dCopperRingsInnerRadius,dCopperRingsOuterRadius, 0.5*dCopperRingsThickness, 0.*deg, 360.*deg);

  m_pCopperRingsLogicalVolume = new G4LogicalVolume(pCopperRingsTubs,Copper,"CopperRingsLogicalVolume");

  G4double pdCopperRingsZ[dCopperRingsNumber];
	
  for(int i = 0; i < dCopperRingsNumber; i++)
    {
      pdCopperRingsZ[i] = dFirstCopperRingOffsetZ-i*(dCopperRingsSpacing+dCopperRingsThickness);

      m_hCopperRingsPhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdCopperRingsZ[i]),
			   m_pCopperRingsLogicalVolume,"CopperRings",m_pLXeLogicalVolume, false, i));

      //G4cout<<"pdCopperRingsZ[i]: "<<pdCopperRingsZ[i]<<G4endl;                                                                               
    }

  const G4double dGridMeshThickness= GetGeometryParameter("GridMeshThickness");
  const G4double dGridMeshRadius= GetGeometryParameter("GridMeshRadius");
  const G4int dGridMeshNumber= GetGeometryParameter("GridMeshNumber");

  const G4double dGridMeshBorderInnerRadius= GetGeometryParameter("GridMeshBorderInnerRadius");
  const G4double dGridMeshBorderOuterRadius= GetGeometryParameter("GridMeshBorderOuterRadius");

  const G4double dGridMeshSupportThickness= GetGeometryParameter("GridMeshSupportThickness");
  const G4double dGridMeshSupportInnerRadius= GetGeometryParameter("GridMeshSupportInnerRadius");
  const G4double dGridMeshSupportOuterRadius= GetGeometryParameter("GridMeshSupportOuterRadius");
  const G4int dGridMeshSupportNumber= GetGeometryParameter("GridMeshSupportNumber");
		
  const G4double dCathodeGridMeshSupportThickness= GetGeometryParameter("CathodeGridMeshSupportThickness");
  const G4double dCathodeGridMeshSupportInnerRadius = GetGeometryParameter("CathodeGridMeshSupportInnerRadius");
  const G4double dCathodeGridMeshSupportOuterRadius = GetGeometryParameter("CathodeGridMeshSupportOuterRadius");
  //const G4double d = GetGeometryParameter("");

  const G4double dAnodeToLiquidLevel = GetGeometryParameter("AnodeToLiquidLevel");
  const G4double dScreeningMeshToAnode = GetGeometryParameter("ScreeningMeshToAnode");
  const G4double dCathodeToScreenMesh = GetGeometryParameter("CathodeToScreenMesh");

  const G4double dCryostatBottomDistanceBetweenLXeAndPTFE = GetGeometryParameter("CryostatBottomDistanceBetweenLXeAndPTFE");
  const G4double dBottomPMTHolderBottomPartThickness = GetGeometryParameter("BottomPMTHolderBottomPartThickness");
  const G4double dBottomPMTHolderTopPartThickness = GetGeometryParameter("BottomPMTHolderTopPartThickness");
  const G4double dBottomPMTBaseHolderThickness = GetGeometryParameter("BottomPMTBaseHolderThickness");
  const G4double dTopBottomPMTToScreeningMesh = GetGeometryParameter("TopBottomPMTToScreeningMesh");

  //================================= Meshes Support ===============================	
  //---------------------------screen meshes + gate mesh + anode -------------------------------
  const G4double dFirstGridMeshSupportOffsetZ = dLXeLevelZeroOffset+dAnodeToLiquidLevel+dScreeningMeshToAnode-0.5*dGridMeshSupportThickness;//(dLXeHalfZ+dAnodeToLiquidLevel+dScreeningMeshToAnode-0.5*dGridMeshSupportThickness);

  G4Tubs *pGridMeshSupportTubs = new G4Tubs("GridMeshSupportTube", dGridMeshSupportInnerRadius,dGridMeshSupportOuterRadius, 0.5*dGridMeshSupportThickness, 0.*deg, 360.*deg);

  m_pGridMeshSupportLogicalVolume = new G4LogicalVolume(pGridMeshSupportTubs,SS316LSteel,"GridMeshSupportLogicalVolume");

  G4double pdGridMeshSupportZ[dGridMeshSupportNumber];

  pdGridMeshSupportZ[0] = dFirstGridMeshSupportOffsetZ;
  pdGridMeshSupportZ[1] = pdGridMeshSupportZ[0]- dScreeningMeshToAnode + 2*0.5*dGridMeshSupportThickness;
  pdGridMeshSupportZ[2] = pdGridMeshSupportZ[1]-0.5*dGridMeshSupportThickness-dAnodeToLiquidLevel - dLiquidLevelToGateMesh-0.5*dGridMeshSupportThickness;
  pdGridMeshSupportZ[3] = pdGridMeshSupportZ[2]+0.5*dGridMeshSupportThickness-(dCopperRingsNumber)*(dCopperRingsSpacing+dCopperRingsThickness)-dCopperRingsSpacing-dCathodeGridMeshSupportThickness-dCathodeToScreenMesh-0.5*dGridMeshSupportThickness;//201.55;

  for(int i = 0; i < dGridMeshSupportNumber; i++)
    {
	
      m_hGridMeshSupportPhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdGridMeshSupportZ[i]),
			   m_pGridMeshSupportLogicalVolume,"GridMeshSupport",m_pLXeLogicalVolume, false, i));

      //G4cout<<"pdGridMeshSupportZ[i]: "<<pdGridMeshSupportZ[i]<<G4endl;                                                                             
    }

  //--------------------------------cathode support------------------------------
  const G4double dCathodeGridMeshSupportOffsetZ = -(dLXeHalfZ-dCryostatBottomDistanceBetweenLXeAndPTFE-dBottomPMTBaseHolderThickness-dBottomPMTHolderBottomPartThickness-dBottomPMTHolderTopPartThickness-dTopBottomPMTToScreeningMesh-dCathodeToScreenMesh-0.5*dCathodeGridMeshSupportThickness);

  G4Tubs *pCathodeGridMeshSupportTubs = new G4Tubs("CathodeGridMeshSupportTube", dCathodeGridMeshSupportInnerRadius,dCathodeGridMeshSupportOuterRadius, 0.5*dCathodeGridMeshSupportThickness, 0.*deg, 360.*deg);

  m_pCathodeGridMeshSupportLogicalVolume = new G4LogicalVolume(pCathodeGridMeshSupportTubs,SS316LSteel,"CathodeGridMeshSupportLogicalVolume");

  m_pCathodeGridMeshSupportPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dCathodeGridMeshSupportOffsetZ),
							      m_pCathodeGridMeshSupportLogicalVolume, "CathodeGridMeshSupport", m_pLXeLogicalVolume, false, 0);

  //G4cout << "Cathode Offset: " <<  dCathodeGridMeshSupportOffsetZ <<G4endl;
  
  //================================= Meshes Borders ===============================	

  const G4double dFirstGridMeshOffsetZ = dLXeLevelZeroOffset+dAnodeToLiquidLevel+dScreeningMeshToAnode-0.5*dGridMeshThickness;//(dLXeHalfZ+dAnodeToLiquidLevel+dScreeningMeshToAnode-0.5*dGridMeshThickness); //top screening mesh

  G4Tubs *pGridMeshBorderTubs = new G4Tubs("GridMeshBorderTube", dGridMeshBorderInnerRadius,dGridMeshBorderOuterRadius, 0.5*dGridMeshThickness, 0.*deg, 360.*deg);

  m_pGridMeshBorderLogicalVolume = new G4LogicalVolume(pGridMeshBorderTubs,SS316LSteel,"GridMeshBorderLogicalVolume");

  G4double pdGridMeshBorderZ[dGridMeshNumber]; 

  pdGridMeshBorderZ[0] = dFirstGridMeshOffsetZ;
  pdGridMeshBorderZ[1] = pdGridMeshBorderZ[0]- dScreeningMeshToAnode+0.5*dGridMeshThickness;
  pdGridMeshBorderZ[2] = pdGridMeshBorderZ[1]-0.5*dGridMeshThickness-dAnodeToLiquidLevel - dLiquidLevelToGateMesh-0.5*dGridMeshThickness;
  pdGridMeshBorderZ[3] = pdGridMeshBorderZ[2]+0.5*dGridMeshThickness-(dCopperRingsNumber)*(dCopperRingsSpacing+dCopperRingsThickness)-dCopperRingsSpacing-dCathodeGridMeshSupportThickness-0.5*dGridMeshThickness;//-dCathodeToScreenMesh-0.5*dGridMeshSupportThickness;
  pdGridMeshBorderZ[4] = pdGridMeshBorderZ[3]-0.5*dGridMeshThickness-dCathodeToScreenMesh-0.5*dGridMeshThickness;//201.55;

	for(int i = 0; i < dGridMeshNumber; i++)
    {
	
      m_hGridMeshBorderPhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdGridMeshBorderZ[i]),
			   m_pGridMeshBorderLogicalVolume,"GridMeshBorder",m_pLXeLogicalVolume, false, i));

      //G4cout<<"pdGridMeshBorderZ[i]: "<<pdGridMeshBorderZ[i]<<G4endl;                                                                            
    }

  //================================= Meshes ===============================	
  //the actual grid - same for all
	const G4double dGXeHalfZ = 0.5*GetGeometryParameter("GXeHeight");
  //	const G4double dFirstGridMeshOffsetZ = (dLXeHalfZ+dAnodeToLiquidLevel+dScreeningMeshToAnode-0.5*dGridMeshThickness);

  G4Tubs *pGridMeshTubs = new G4Tubs("GridMeshTube", 0,dGridMeshRadius, 0.5*dGridMeshThickness, 0.*deg, 360.*deg);

   m_pGridMeshGXeLogicalVolume = new G4LogicalVolume(pGridMeshTubs,GridMeshSS316LSteelGXe,"GXeGridMeshLogicalVolume");
   m_pGridMeshLXeLogicalVolume = new G4LogicalVolume(pGridMeshTubs,GridMeshSS316LSteelLXe,"LXeGridMeshLogicalVolume");

  G4double pdGridMeshZ[dGridMeshNumber]; 

  pdGridMeshZ[0] = dFirstGridMeshOffsetZ; //top screening mesh wrt LXe
  pdGridMeshZ[1] = pdGridMeshZ[0]- dScreeningMeshToAnode+0.5*dGridMeshThickness;
  pdGridMeshZ[2] = pdGridMeshZ[1]-0.5*dGridMeshThickness-dAnodeToLiquidLevel - dLiquidLevelToGateMesh-0.5*dGridMeshThickness;
  pdGridMeshZ[3] = pdGridMeshZ[2]+0.5*dGridMeshThickness-(dCopperRingsNumber)*(dCopperRingsSpacing+dCopperRingsThickness)-dCopperRingsSpacing-dCathodeGridMeshSupportThickness-0.5*dGridMeshThickness;
  pdGridMeshZ[4] = pdGridMeshZ[3]-0.5*dGridMeshThickness-dCathodeToScreenMesh-0.5*dGridMeshThickness;//201.55;

  pdGridMeshZ[0] =(-dGXeHalfZ+dAnodeToLiquidLevel+dScreeningMeshToAnode-0.5*dGridMeshThickness);//top screening mesh wrt GXe
pdGridMeshZ[1] = pdGridMeshZ[0]- dScreeningMeshToAnode+0.5*dGridMeshThickness;
	
  for(int i = 0; i < dGridMeshNumber; i++)
    {

      if (i<2)
		m_hGridMeshPhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdGridMeshZ[i]),
			   m_pGridMeshGXeLogicalVolume,"GXeGridMesh",m_pGXeLogicalVolume, false, i));
      else
  	m_hGridMeshPhysicalVolumes.push_back
	(new G4PVPlacement(0,G4ThreeVector(0,0,pdGridMeshZ[i]),
			   m_pGridMeshLXeLogicalVolume,"LXeGridMesh",m_pLXeLogicalVolume, false, i));

      //G4cout<<"pdGridMeshZ[i]: "<<pdGridMeshZ[i]<<G4endl;                                                                
    }

  //=============================== optical surfaces ==============================
  // no need for meshes as same index of refraction as LXe

  //================================== attributes =================================
  G4Colour hGridMeshColor(0.5, 0.5, 0.5,0.5); //gray
  G4Colour hGridRingColor(0.04, 0.72, 0.0, 1.);
  G4Colour hGridSupportColor(0.04, 0.72, 0.0, 1.);//green
  G4Colour hTeflonColor(1., 1., 1., PTFECylinderAlphaChannel);
  G4Colour hCopperColor(0.835, 0.424, 0.059, CopperRingsAlphaChannel);
  G4Colour hRedColor(1.,0., 0., 0.5);
  // G4VisAttributes *pRedVisAtt = new G4VisAttributes(hRedColor);
  // pRedVisAtt->SetVisibility(true);

  G4VisAttributes *pCopperRingsVisAtt = new G4VisAttributes(hCopperColor);
  pCopperRingsVisAtt->SetVisibility(true);
  m_pCopperRingsLogicalVolume->SetVisAttributes(pCopperRingsVisAtt);

  G4VisAttributes *pGridMeshSupportVisAtt = new G4VisAttributes(hGridSupportColor);
  pGridMeshSupportVisAtt->SetVisibility(true);
  m_pGridMeshSupportLogicalVolume->SetVisAttributes(pGridMeshSupportVisAtt);

  G4VisAttributes *pCathodeGridMeshSupportVisAtt = new G4VisAttributes(hGridSupportColor);
  pCathodeGridMeshSupportVisAtt->SetVisibility(true);
  m_pCathodeGridMeshSupportLogicalVolume->SetVisAttributes(pCathodeGridMeshSupportVisAtt);

  G4VisAttributes *pGridMeshVisAtt = new G4VisAttributes(hGridMeshColor);
  pGridMeshVisAtt->SetVisibility(true);
  m_pGridMeshLXeLogicalVolume->SetVisAttributes(pGridMeshVisAtt);
  m_pGridMeshGXeLogicalVolume->SetVisAttributes(pGridMeshVisAtt);

  G4VisAttributes *pGridMeshBorderVisAtt = new G4VisAttributes(hGridMeshColor);
  pGridMeshBorderVisAtt->SetVisibility(true);
  m_pGridMeshBorderLogicalVolume->SetVisAttributes(pGridMeshBorderVisAtt);
}

//******************************************************************/
// ConstructPmtArrays
//******************************************************************/
void muensterTPCDetectorConstruction::ConstructPmtArrays() {
  //G4cout << "----- ConstructPMTs" << G4endl;
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< pmt arrays >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  const G4double dPmtWindowWidth = GetGeometryParameter("PmtWindowWidth");
  const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");

  const G4double dPmtCasingWidth = GetGeometryParameter("PmtCasingWidth");
  const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
  const G4double dPmtCasingThickness = GetGeometryParameter("PmtCasingThickness");

  const G4double dPmtPhotoCathodeWidth = GetGeometryParameter("PmtPhotoCathodeWidth");
  const G4double dPmtPhotoCathodeThickness = GetGeometryParameter("PmtPhotoCathodeThickness");

  const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");

  const G4double dSpaceBelowTopPMTHeight= GetGeometryParameter("SpaceBelowTopPMTHeight");

  G4Material *Quartz = G4Material::GetMaterial("Quartz"); //for the window
  G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel"); //for the frames
  G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
  G4Material *Aluminium = G4Material::GetMaterial("PhotoCathodeAluminium"); //for the photocathode
  G4Material *Cirlex = G4Material::GetMaterial("Cirlex"); //for the base
  G4Material *GXe = G4Material::GetMaterial("GXe");
  G4Material *GXeTeflon = G4Material::GetMaterial("GXeTeflon");

  //==================================== pmts =====================================
  //G4cout << "----- PMTs" << G4endl;
  //--------------------------------- pmt window ----------------------------------
  const G4double dPmtWindowHalfX = 0.5*dPmtWindowWidth;
  const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;

  G4Box *pPmtWindowBox = new G4Box("PmtWindowBox", dPmtWindowHalfX, dPmtWindowHalfX, dPmtWindowHalfZ);

  m_pPmtWindowLogicalVolume = new G4LogicalVolume(pPmtWindowBox, Quartz, "PmtWindowLogicalVolume", 0, 0, 0);

  //--------------------------------- pmt casing ----------------------------------
  const G4double dPmtCasingHalfX = 0.5*dPmtCasingWidth;
  const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;

  G4Box *pPmtCasingBox = new G4Box("PmtCasingBox", dPmtCasingHalfX, dPmtCasingHalfX, dPmtCasingHalfZ);

  m_pPmtCasingLogicalVolume = new G4LogicalVolume(pPmtCasingBox, SS304LSteel, "PmtCasingLogicalVolume", 0, 0, 0);

  //-------------------------------- pmt interior ---------------------------------
  const G4double dPmtInteriorHalfX = dPmtCasingHalfX-dPmtCasingThickness;
  const G4double dPmtInteriorHalfZ = dPmtCasingHalfZ-dPmtCasingThickness;

  G4Box *pPmtInteriorBox = new G4Box("PmtInteriorBox", dPmtInteriorHalfX, dPmtInteriorHalfX, dPmtInteriorHalfZ);

  m_pPmtInteriorLogicalVolume = new G4LogicalVolume(pPmtInteriorBox, Vacuum, "PmtInteriorLogicalVolume", 0, 0, 0);

  m_pPmtInteriorPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
						   m_pPmtInteriorLogicalVolume, "PmtInterior", m_pPmtCasingLogicalVolume, false, 0);

  //------------------------------ pmt photocathode -------------------------------
  const G4double dPmtPhotoCathodeHalfX = 0.5*dPmtPhotoCathodeWidth;
  const G4double dPmtPhotoCathodeHalfZ = 0.5*dPmtPhotoCathodeThickness;
  const G4double dPmtPhotoCathodeOffsetZ = dPmtCasingHalfZ-dPmtPhotoCathodeHalfZ;

  G4Box *pPmtPhotoCathodeBox = new G4Box("PmtPhotoCathodeBox", dPmtPhotoCathodeHalfX, dPmtPhotoCathodeHalfX, dPmtPhotoCathodeHalfZ);

  m_pPmtPhotoCathodeLogicalVolume = new G4LogicalVolume(pPmtPhotoCathodeBox, Aluminium, "PmtPhotoCathodeLogicalVolume", 0, 0, 0);

  m_pPmtPhotoCathodePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., -dPmtPhotoCathodeOffsetZ),
						       m_pPmtPhotoCathodeLogicalVolume, "PmtPhotoCathode", m_pPmtCasingLogicalVolume, false, 0);

  //---------------------------------- pmt base -----------------------------------
  const G4double dPmtBaseHalfX = dPmtCasingHalfX;
  const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;

  G4Box *pPmtBaseBox = new G4Box("PmtBaseBox", dPmtBaseHalfX, dPmtBaseHalfX, dPmtBaseHalfZ);

  m_pPmtBaseLogicalVolume = new G4LogicalVolume(pPmtBaseBox, Cirlex, "PmtBaseLogicalVolume", 0, 0, 0);

  //------------------------------------space below top Pmts -----------------------------
  const G4double dSpaceBelowTopPMTHalfZ = 0.5*dSpaceBelowTopPMTHeight;

  G4Box *pSpaceBelowTopPMTBox = new G4Box("SpaceBelowTopPMTBox", dPmtCasingHalfX, dPmtCasingHalfX, dSpaceBelowTopPMTHalfZ);

  m_pSpaceBelowTopPMTLogicalVolume = new G4LogicalVolume(pSpaceBelowTopPMTBox, GXe, "SpaceBelowTopPMTLogicalVolume", 0, 0, 0);

  //================================== top array ==================================
  G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
	
  stringstream hVolumeName;
  for(G4int iPmtNb=0; iPmtNb<iNbTopPmts; iPmtNb++)
    {
      //G4cout << "----- TopPMT " << iPmtNb << " of " << iNbTopPmts-1 << G4endl;
      hVolumeName.str(""); hVolumeName << "PmtWindowNo" << iPmtNb;
      m_hPmtWindowPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_WINDOW),
							      GetPmtPosition(iPmtNb, PMT_WINDOW), 
							      m_pPmtWindowLogicalVolume,
							      hVolumeName.str(), 
							      m_pPTFETopPMTHolderLogicalVolume,
							      false, iPmtNb));

      hVolumeName.str(""); hVolumeName << "PmtCasingNo" << iPmtNb;

      m_hPmtCasingPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_CASING),
							      GetPmtPosition(iPmtNb, PMT_CASING), 
							      m_pPmtCasingLogicalVolume,
							      hVolumeName.str(), 
							      m_pPTFETopPMTHolderLogicalVolume,
							      false, iPmtNb));

      hVolumeName.str(""); hVolumeName << "SpaceBelowTopPMTNo" << iPmtNb;

      m_hSpaceBelowTopPMTPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_SPACE),
							      GetPmtPosition(iPmtNb, PMT_SPACE), 
							      m_pSpaceBelowTopPMTLogicalVolume,
							      hVolumeName.str(), 
							      m_pPTFETopPMTHolderLogicalVolume, 
							      false, iPmtNb));

      hVolumeName.str(""); hVolumeName << "PmtBaseNo" << iPmtNb;

      m_hPmtBasePhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_BASE),
							    GetPmtPosition(iPmtNb, PMT_BASE), 
							    m_pPmtBaseLogicalVolume, hVolumeName.str(),
							    m_pPTFETopPMTHolderLogicalVolume, 
							    false, iPmtNb));
    }
	
  //================================ bottom array =================================
  G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");

  for(G4int iPmtNb=iNbTopPmts; iPmtNb<iNbTopPmts+iNbBottomPmts; iPmtNb++)
    {
      //G4cout << "----- BottomPMT " << iPmtNb << " of " << iNbTopPmts+iNbBottomPmts-1 << G4endl;
      hVolumeName.str(""); hVolumeName << "PmtWindowNo" << iPmtNb;
      m_hPmtWindowPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_WINDOW),
							      GetPmtPosition(iPmtNb, PMT_WINDOW), 
							      m_pPmtWindowLogicalVolume,
							      hVolumeName.str(), 
							      m_pPTFEBottomPMTHolderLogicalVolume, 
							      false, iPmtNb));

      hVolumeName.str(""); hVolumeName << "PmtCasingNo" << iPmtNb;

      m_hPmtCasingPhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_CASING),
							      GetPmtPosition(iPmtNb, PMT_CASING), 
							      m_pPmtCasingLogicalVolume,
							      hVolumeName.str(),
							      m_pPTFEBottomPMTHolderLogicalVolume, 
							      false, iPmtNb));

      hVolumeName.str(""); hVolumeName << "PmtBaseNo" << iPmtNb;

      m_hPmtBasePhysicalVolumes.push_back(new G4PVPlacement(GetPmtRotation(iPmtNb, PMT_BASE),
      							    GetPmtPosition(iPmtNb, PMT_BASE), 
							    m_pPmtBaseLogicalVolume, hVolumeName.str(),
      							    m_pPTFEBottomPMTHolderLogicalVolume, 
							    false, iPmtNb));

      // m_pLXeLogicalVolume,
      //G4cout << hVolumeName.str() << G4endl;
    }

  //------------------------------- pmt sensitivity -------------------------------
  //G4cout << "----- PMT sensitivity " << G4endl;
  G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

  muensterTPCPmtSensitiveDetector *pPmtSD = new muensterTPCPmtSensitiveDetector("muensterTPC/PmtSD");
  pSDManager->AddNewDetector(pPmtSD);
  m_pPmtPhotoCathodeLogicalVolume->SetSensitiveDetector(pPmtSD);

  //================================== optical surface =================================	
  //G4cout << "----- optical surface " << G4endl;
	G4OpticalSurface *pSS304LSteelOpticalSurface = new G4OpticalSurface("SS304LSteelOpticalSurface",
		unified, polished, dielectric_metal);
        pSS304LSteelOpticalSurface->SetMaterialPropertiesTable(SS304LSteel->GetMaterialPropertiesTable());
	
	stringstream hStream;
	for(G4int iPmtNb=iNbTopPmts; iPmtNb<iNbTopPmts+iNbBottomPmts; iPmtNb++)
	{
		hStream.str("PmtCasingLogicalBorderSurfacePmt");
		hStream << iPmtNb;
		new G4LogicalBorderSurface(hStream.str(),
		   m_pLXePhysicalVolume, m_hPmtCasingPhysicalVolumes[iPmtNb], pSS304LSteelOpticalSurface);
	}

	stringstream hStream2;
	for(G4int iPmtNb=0; iPmtNb<iNbTopPmts; iPmtNb++)
	{
		hStream2.str("PmtCasingGXeLogicalBorderSurfacePmt");
		hStream2 << iPmtNb;
		new G4LogicalBorderSurface(hStream2.str(),
		   m_pGXePhysicalVolume, m_hPmtCasingPhysicalVolumes[iPmtNb], pSS304LSteelOpticalSurface);
	}

  //---------------------------------- attributes ---------------------------------
  //G4cout << "----- attributes " << G4endl;
  G4Colour hPmtWindowColor(0.4,0.804, 0.666,0.75);  //102 205 170 aquamarine
  G4VisAttributes *pPmtWindowVisAtt = new G4VisAttributes(hPmtWindowColor);
  pPmtWindowVisAtt->SetVisibility(true);
  m_pPmtWindowLogicalVolume->SetVisAttributes(pPmtWindowVisAtt);

  G4Colour hPmtCasingColor(0.623, 0.713, 0.804,0.75); //159 182 205 slate gray
  G4VisAttributes *pPmtCasingVisAtt = new G4VisAttributes(hPmtCasingColor);
  pPmtCasingVisAtt->SetVisibility(true);
  m_pPmtCasingLogicalVolume->SetVisAttributes(pPmtCasingVisAtt);

  G4VisAttributes *pPmtInteriorVisAtt = new G4VisAttributes();
  pPmtInteriorVisAtt->SetVisibility(true);
  m_pPmtInteriorLogicalVolume->SetVisAttributes(pPmtInteriorVisAtt);	

  G4Colour hPmtPhotoCathodeColor(1., 0.725, 0.059,1.0); //255 185 15 gold
  G4VisAttributes *pPmtPhotoCathodeVisAtt = new G4VisAttributes(hPmtPhotoCathodeColor);
  pPmtPhotoCathodeVisAtt->SetVisibility(true);
  m_pPmtPhotoCathodeLogicalVolume->SetVisAttributes(pPmtPhotoCathodeVisAtt);

  G4Colour hPmtBaseColor(1., 0.204, 0.702,0.75); //255;52;179 pink
  G4VisAttributes *pPmtBaseVisAtt = new G4VisAttributes(hPmtBaseColor);
  pPmtBaseVisAtt->SetVisibility(true);
  m_pPmtBaseLogicalVolume->SetVisAttributes(pPmtBaseVisAtt);

  G4Colour hPmtSpaceColor(0, 1, 1);
  G4VisAttributes *pPmtSpaceVisAtt = new G4VisAttributes(hPmtSpaceColor);
  pPmtSpaceVisAtt->SetVisibility(true);
  m_pSpaceBelowTopPMTLogicalVolume->SetVisAttributes(pPmtSpaceVisAtt);
  //G4cout << "----- ConstructPMTs -- END." << G4endl;
}

//******************************************************************/
// GetPmtPosition
//******************************************************************/
G4ThreeVector muensterTPCDetectorConstruction::GetPmtPosition(G4int iPmtNb, PmtPart ePmtPart) {
//iPmtNb is the pmt number so ranging from 0 to 6 in both top and bottom array, ePmtPart is the different cases found in GetPmtPositionTopArray i.e PMT_WINDOW, PMT_CASING etc etc...
  const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
  const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");

  G4ThreeVector hPos;

  if(iPmtNb < iNbTopPmts)
    hPos = GetPmtPositionTopArray(iPmtNb, ePmtPart);
  else if(iPmtNb < iNbTopPmts+iNbBottomPmts)
    hPos = GetPmtPositionBottomArray(iPmtNb, ePmtPart);
	
  return hPos;
}

//******************************************************************/
// GetPmtPositionTopArray
//******************************************************************/
G4ThreeVector
muensterTPCDetectorConstruction::GetPmtPositionTopArray(G4int iPmtNb, PmtPart ePmtPart)
{

  const G4double dLXeLevelZeroOffset = GetGeometryParameter("LXeLevelZeroOffset");

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< top array >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  const G4double dLXeHeight = GetGeometryParameter("LXeHeight");

  const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");
  const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
  const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");

  const G4double dPmtToPmtBase= GetGeometryParameter("PmtToPmtBase");
  const G4double dBottomTopPMTToScreeningMesh= GetGeometryParameter("BottomTopPMTToScreeningMesh");
 
  const G4double dSpaceBelowTopPMTHeight= GetGeometryParameter("SpaceBelowTopPMTHeight");

  const G4double dAnodeToLiquidLevel = GetGeometryParameter("AnodeToLiquidLevel");
  const G4double dScreeningMeshToAnode = GetGeometryParameter("ScreeningMeshToAnode");
  //	const G4double= GetGeometryParameter("");
  //	const G4double= GetGeometryParameter("");

  vector<G4int> hTopPmtsPerRow;
  hTopPmtsPerRow.push_back((G4int) GetGeometryParameter("NbTopPmtsFirstRow"));
  hTopPmtsPerRow.push_back((G4int) GetGeometryParameter("NbTopPmtsSecondRow"));
  hTopPmtsPerRow.push_back((G4int) GetGeometryParameter("NbTopPmtsThirdRow"));

  //================================ top array ================================-
  const G4double dLXeHalfZ = 0.5*dLXeHeight;

  const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;
  const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;
  const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;
  const G4double dPmtSpaceHalfZ = 0.5*dSpaceBelowTopPMTHeight;

  const G4double dTopPmtBaseOffsetZ = 0.5*dPmtCasingHeight+dPmtToPmtBase+dPmtBaseHalfZ;//dLXeLevelZeroOffset+dAnodeToLiquidLevel+dScreeningMeshToAnode+dBottomTopPMTToScreeningMesh+dPmtWindowThickness+dPmtCasingHeight+dPmtToPmtBase+dPmtBaseHalfZ;
  const G4double dTopPmtCasingOffsetZ =(dTopPmtBaseOffsetZ-dPmtBaseHalfZ-dPmtToPmtBase-dPmtCasingHalfZ);
  const G4double dTopPmtWindowOffsetZ =dTopPmtCasingOffsetZ-dPmtCasingHalfZ-dPmtWindowHalfZ;
  const G4double dTopPmtSpaceOffsetZ=dTopPmtCasingOffsetZ-dPmtCasingHalfZ-dPmtWindowHalfZ-dPmtSpaceHalfZ;

  //G4cout << "dTopPmtBaseOffsetZ: " << dTopPmtBaseOffsetZ<< G4endl;
  //G4cout << "dTopPmtCasingOffsetZ: " << dTopPmtCasingOffsetZ<< G4endl;
  //G4cout << "dTopPmtWindowOffsetZ: " << dTopPmtWindowOffsetZ<< G4endl;

  G4ThreeVector hPos = ComputeXYPmtPositionForGridPattern(iPmtNb, hTopPmtsPerRow);

  switch(ePmtPart)
    {
    case PMT_WINDOW:
      hPos.setZ(dTopPmtWindowOffsetZ);
      break;

    case PMT_CASING:
      hPos.setZ(dTopPmtCasingOffsetZ);
      break;

    case PMT_BASE:
      hPos.setZ(dTopPmtBaseOffsetZ);
      break;

    case PMT_SPACE:
      hPos.setZ(dTopPmtSpaceOffsetZ);
      break;

    }

  return hPos;
}

//******************************************************************/
// GetPmtPositionBottomArray
//******************************************************************/
G4ThreeVector muensterTPCDetectorConstruction::GetPmtPositionBottomArray(G4int iPmtNb, PmtPart ePmtPart) {
  //const G4double dLXeLevelZeroOffset =GetGeometryParameter("LXeLevelZeroOffset");

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bottom array >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  const G4double dLXeHeight = GetGeometryParameter("LXeHeight");

  const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");

  const G4double dPmtWindowThickness = GetGeometryParameter("PmtWindowThickness");
  const G4double dPmtCasingHeight = GetGeometryParameter("PmtCasingHeight");
  const G4double dPmtBaseThickness = GetGeometryParameter("PmtBaseThickness");
  const G4double dPmtPhotoCathodeThickness= GetGeometryParameter("PmtPhotoCathodeThickness");
  
  const G4double dPmtToPmtBase= GetGeometryParameter("PmtToPmtBase");
  
  const G4double dBottomPMTBaseHolderThickness= GetGeometryParameter("BottomPMTBaseHolderThickness");
  const G4double dBottomPMTHolderTopPartThickness= GetGeometryParameter("BottomPMTHolderTopPartThickness");
  const G4double dBottomPMTHolderBottomPartThickness= GetGeometryParameter("BottomPMTHolderBottomPartThickness");
  const G4double dCryostatBottomDistanceBetweenLXeAndPTFE= GetGeometryParameter("CryostatBottomDistanceBetweenLXeAndPTFE");

  const G4double dPTFEBottomPMTHolderHalfZ = 0.5*(dBottomPMTHolderTopPartThickness+dBottomPMTHolderBottomPartThickness+dBottomPMTBaseHolderThickness);

  vector<G4int> hBottomPmtsPerRow;
  hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsFirstRow"));
  hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsSecondRow"));
  hBottomPmtsPerRow.push_back((G4int) GetGeometryParameter("NbBottomPmtsThirdRow"));

  //================================ bottom array ================================-
  const G4double dLXeHalfZ = 0.5*dLXeHeight;

  const G4double dPmtWindowHalfZ = 0.5*dPmtWindowThickness;
  const G4double dPmtCasingHalfZ = 0.5*dPmtCasingHeight;
  const G4double dPmtBaseHalfZ = 0.5*dPmtBaseThickness;

  //const G4double dBottomPmtBaseOffsetZ =dPTFEBottomPMTHolderHalfZ-dPmtCasingHalfZ-dPmtWindowThickness-dPmtPhotoCathodeThickness-dPmtBaseThickness;
  const G4double dBottomPmtBaseOffsetZ =(dPTFEBottomPMTHolderHalfZ-dPmtCasingHeight-dPmtWindowThickness-dPmtToPmtBase-dPmtBaseHalfZ);

  const G4double dBottomPmtCasingOffsetZ = (dBottomPmtBaseOffsetZ+dPmtToPmtBase+dPmtBaseHalfZ+dPmtCasingHalfZ);
  const G4double dBottomPmtWindowOffsetZ = dBottomPmtCasingOffsetZ+dPmtCasingHalfZ+dPmtWindowHalfZ;

  //G4cout << "dBottomPmtBaseOffsetZ: " << dBottomPmtBaseOffsetZ<< G4endl;
  //G4cout << "dBottomPmtCasingOffsetZ: " << dBottomPmtCasingOffsetZ<< G4endl;
  //G4cout << "dBottomPmtWindowOffsetZ: " << dBottomPmtWindowOffsetZ<< G4endl;

  G4ThreeVector hPos = ComputeXYPmtPositionForGridPattern(iPmtNb-iNbTopPmts, hBottomPmtsPerRow);

  switch(ePmtPart)
    {
    case PMT_WINDOW:
      hPos.setZ(dBottomPmtWindowOffsetZ);
      break;

    case PMT_CASING:
      hPos.setZ(dBottomPmtCasingOffsetZ);
      break;

    case PMT_BASE:
      hPos.setZ(dBottomPmtBaseOffsetZ);
      break;
    }

  return hPos;
}

//******************************************************************/
// GetPmtRotation
//******************************************************************/
G4RotationMatrix* muensterTPCDetectorConstruction::GetPmtRotation(G4int iPmtNb, PmtPart ePmtPart) {
  const G4int iNbTopPmts = (G4int) GetGeometryParameter("NbTopPmts");
  const G4int iNbBottomPmts = (G4int) GetGeometryParameter("NbBottomPmts");


  G4RotationMatrix *pRotationMatrix;

  if(iPmtNb < iNbTopPmts )
     pRotationMatrix = m_pRotationX0;
  else if( iPmtNb < iNbTopPmts+iNbBottomPmts)
    pRotationMatrix = m_pRotationX180;

  return pRotationMatrix;
}

//******************************************************************/
// ComputeXYPmtPositionForGridPattern
//******************************************************************/
G4ThreeVector muensterTPCDetectorConstruction::ComputeXYPmtPositionForGridPattern(G4int iPmtNb, const vector<G4int> &hPmtsPerRow) {
  G4double dPmtWidth = muensterTPCDetectorConstruction::GetGeometryParameter("PmtWidth");
  G4double dPmtSpacing = muensterTPCDetectorConstruction::GetGeometryParameter("PmtSpacing");

  G4int iRow = 0, iColumn;
  G4int iTotal = hPmtsPerRow[0];

  while(iPmtNb > iTotal-1)
    iTotal += hPmtsPerRow[++iRow];
	
  iColumn = hPmtsPerRow[iRow]-(iTotal-iPmtNb);

  G4double dX = -(hPmtsPerRow[iRow]-1)*(dPmtWidth+dPmtSpacing)/2. + iColumn*(dPmtWidth+dPmtSpacing);
  G4double dY = (hPmtsPerRow.size()-1)*(dPmtWidth+dPmtSpacing)/2. - iRow*(dPmtWidth+dPmtSpacing);

  return G4ThreeVector(dX, dY, 0.);
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Material Reflectivity <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//******************************************************************/
// SetTeflonReflectivity
//******************************************************************/
void muensterTPCDetectorConstruction::SetTeflonReflectivity(G4double dReflectivity) {
  G4Material *pTeflonMaterial = G4Material::GetMaterial(G4String("Teflon"));

  if(pTeflonMaterial)
    {
      G4cout << "----> Setting Teflon reflectivity to " << dReflectivity << G4endl;

      G4MaterialPropertiesTable *pTeflonPropertiesTable = pTeflonMaterial->GetMaterialPropertiesTable();
		
      G4double teflon_PP[] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
      G4double teflon_REFL[] = {dReflectivity, dReflectivity, dReflectivity};
      pTeflonPropertiesTable->RemoveProperty("REFLECTIVITY");
      pTeflonPropertiesTable->AddProperty("REFLECTIVITY", teflon_PP, teflon_REFL, 3);
    }
  else
    {
      G4cout << "!!!!> Teflon material not found!" << G4endl;
      exit(-1);
    }
}

//******************************************************************/
// SetSS304LSteelReflectivity
//******************************************************************/
void muensterTPCDetectorConstruction::SetSS304LSteelReflectivity(G4double dSteelReflectivity) {
  G4Material *pSteelMaterial = G4Material::GetMaterial(G4String("SS304LSteel"));

  if(pSteelMaterial)
    {
      G4cout << "----> Setting SS304LSteel reflectivity to " << dSteelReflectivity << G4endl;

      G4MaterialPropertiesTable *pSteelPropertiesTable = pSteelMaterial->GetMaterialPropertiesTable();
		
      G4double Steel_PP[] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
      G4double Steel_REFL[] = {dSteelReflectivity, dSteelReflectivity, dSteelReflectivity};
      pSteelPropertiesTable->RemoveProperty("REFLECTIVITY");
      pSteelPropertiesTable->AddProperty("REFLECTIVITY", Steel_PP, Steel_REFL, 3);
    }
  else
    {
      G4cout << "!!!!> SS304LSteel material not found!" << G4endl;
      exit(-1);
    }
}


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Material Dependent Physical Processes <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//******************************************************************/
// SetLXeScintillation
//******************************************************************/
void muensterTPCDetectorConstruction::SetLXeScintillation(G4bool bScintillation) {
  G4cout << "----> Setting LXe(GXe) scintillation to " << bScintillation << G4endl;
			
	G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));
	if(pLXeMaterial) 
	{	
		G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
		if(bScintillation)
		{
			//pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 10./(1.0*keV));
  		pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1./(21.6*eV));
			//pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 100000./(1.0*keV));//create 1e5 photons per keV
		}    
	}
  else
  {
  	G4cout << "ls!> LXe materials not found!" << G4endl;
    exit(-1);
  }
	
  G4Material *pGXeMaterial = G4Material::GetMaterial(G4String("GXe"));
  if(pGXeMaterial)
  {	
  	G4MaterialPropertiesTable *pGXePropertiesTable = pGXeMaterial->GetMaterialPropertiesTable();
    if(bScintillation)
		{
			//pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1000/(1.0*keV));
			pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1./(21.6*eV));
		}  
	}
  else
  {
  	G4cout << "ls!> GXe materials not found!" << G4endl;
    exit(-1);
  }
}

//******************************************************************/
// SetLXeLevel
//******************************************************************/
void muensterTPCDetectorConstruction::SetLXeLevel(G4double dlevel) {
	G4cout << "----> Setting LXe(GXe) level to " << dlevel/cm << " cm" << G4endl;
	G4cout << "----> default: 0 mm | no LXe: X mm | no GXe: -Y mm" << G4endl;
	
 	// Open geometry for the physical volume to be modified ...
 	//G4GeometryManager::OpenGeometry(m_pGXePhysicalVolume);
 
 	// Modify dimension of the solid ...
 	//m_pGXePhysicalVolume->GetLogicalVolume()->GetSolid()->SetXHalfLength(dlevel);
 
 	// Close geometry for the portion modified ...
 	//G4GeometryManager::CloseGeometry(m_pGXePhysicalVolume);
	
	//m_pGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dGXeOffsetZ),
	//				   m_pGXeLogicalVolume, "GXe", m_pLXeLogicalVolume, false, 0);

}

//******************************************************************/
// SetLXeMeshMaterial
//******************************************************************/
void muensterTPCDetectorConstruction::SetLXeMeshMaterial(const G4String& name) {
	// search the material by its name
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if(mat && mat != LXeMaterial) {
		m_pGridMeshLXeLogicalVolume->SetMaterial(mat);
    UpdateGeometry();
		G4cout << "----> New LXeMesh material " << mat->GetName() << G4endl;
  }

  if(!mat) {
    G4cout << "--> warning from DetectorConstruction::SetLXeMeshMaterial : "
           << name << " not found" << G4endl;  
  } 
}

//******************************************************************/
// SetGXeMeshMaterial
//******************************************************************/
void muensterTPCDetectorConstruction::SetGXeMeshMaterial(const G4String& name) {
	// search the material by its name
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if(mat && mat != LXeMaterial) {
		m_pGridMeshGXeLogicalVolume->SetMaterial(mat);
    UpdateGeometry();
		G4cout << "----> New GXeMesh material " << mat->GetName() << G4endl;
  }

  if(!mat) {
    G4cout << "--> warning from DetectorConstruction::SetGXeMeshMaterial : "
           << name << " not found" << G4endl;  
  } 
}

//******************************************************************/
// SetMaterial
//******************************************************************/
void muensterTPCDetectorConstruction::SetMaterial(const G4String& name) {
  // search the material by its name
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if(mat && mat != LXeMaterial) {
    LXeMaterial = mat;
		m_pLXeLogicalVolume->SetMaterial(LXeMaterial);
    m_pLXePhysicalVolume->SetName(LXeMaterial->GetName());
		//m_hGridMeshPhysicalVolumes.at(2)->SetLogicalVolume(m_pGridMeshGXeLogicalVolume);
		//m_hGridMeshPhysicalVolumes.at(2)->SetName(m_hGridMeshPhysicalVolumes.at(0)->GetName());
		//m_hGridMeshPhysicalVolumes.at(3)->SetLogicalVolume(m_pGridMeshGXeLogicalVolume);
		//m_hGridMeshPhysicalVolumes.at(3)->SetName(m_hGridMeshPhysicalVolumes.at(0)->GetName());
		//m_hGridMeshPhysicalVolumes.at(4)->SetLogicalVolume(m_pGridMeshGXeLogicalVolume);
		//m_hGridMeshPhysicalVolumes.at(4)->SetName(m_hGridMeshPhysicalVolumes.at(0)->GetName());
    UpdateGeometry();
		G4cout << "----> New material " << LXeMaterial->GetName() << G4endl;
  }

  if(!mat) {
    G4cout << "--> warning from DetectorConstruction::SetMaterial : "
           << name << " not found" << G4endl;  
  } 
}

//******************************************************************/
// SetLXeAbsorbtionLength
//******************************************************************/
void muensterTPCDetectorConstruction::SetLXeAbsorbtionLength(G4double dAbsorbtionLength) {
  G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));

  if(pLXeMaterial)
    {
      G4cout << "----> Setting LXe absorbtion length to " << dAbsorbtionLength/cm << " cm" << G4endl;

      G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
			
      G4double LXe_PP[] = {6.91*eV, 6.98*eV, 7.05*eV};
      G4double LXe_ABSL[] = {dAbsorbtionLength, dAbsorbtionLength, dAbsorbtionLength};
      pLXePropertiesTable->RemoveProperty("ABSLENGTH");
      pLXePropertiesTable->AddProperty("ABSLENGTH", LXe_PP, LXe_ABSL, 3);
    }
  else
    {
      G4cout << "ls!> LXe materials not found!" << G4endl;
      exit(-1);
    }
}

//******************************************************************/
// SetLXeRayScatterLength
//******************************************************************/
void muensterTPCDetectorConstruction::SetLXeRayScatterLength(G4double dRayScatterLength) {
  G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));
  
  if(pLXeMaterial)
    {
      
      G4cout << "----> Setting LXe scattering length to " << dRayScatterLength/cm << " cm" << G4endl;
      
      G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
              
      G4double LXe_PP[] = {6.91*eV, 6.98*eV, 7.05*eV};
      G4double LXe_SCAT[] = {dRayScatterLength, dRayScatterLength, dRayScatterLength};
      pLXePropertiesTable->RemoveProperty("RAYLEIGH");
      pLXePropertiesTable->AddProperty("RAYLEIGH", LXe_PP, LXe_SCAT, 3);
 
    }
  else
    {
      G4cout << "ls!> LXe materials not found!" << G4endl;
      exit(-1);
    }

}

void muensterTPCDetectorConstruction::SetGXeMeshTransparency(G4double dTransparency)
{
  G4Material *pMeshMaterial = G4Material::GetMaterial(G4String("GridMeshSS316LSteelGXe"));
  
  if(pMeshMaterial)
    {
      G4cout << "----> Setting grid transparency to " << dTransparency*100 << " %" << G4endl;
      G4double dAbsorptionLength =  ((G4double)GetGeometryParameter("GridMeshThickness"))/(-log(dTransparency)); 
      G4MaterialPropertiesTable *pGateMeshPropertiesTable = pMeshMaterial->GetMaterialPropertiesTable();
      const G4int iNbEntries = 3;
      G4double pdMeshPhotonMomentum[iNbEntries] = {6.91*eV, 6.98*eV, 7.05*eV};
      G4double pdMeshAbsorptionLength[iNbEntries] = {dAbsorptionLength, dAbsorptionLength, dAbsorptionLength};
      pGateMeshPropertiesTable->RemoveProperty("ABSLENGTH");
      pGateMeshPropertiesTable->AddProperty("ABSLENGTH", pdMeshPhotonMomentum, pdMeshAbsorptionLength, iNbEntries);
   }
  else
    {
      G4cout << "ls!> GridMeshSS316LSteelGXe not found!" << G4endl;
      exit(-1);
    }
}

void muensterTPCDetectorConstruction::SetLXeMeshTransparency(G4double dTransparency)
{
  G4Material *pMeshMaterial = G4Material::GetMaterial(G4String("GridMeshSS316LSteelLXe"));
  
  if(pMeshMaterial)
    {
      G4cout << "----> Setting grid transparency to " << dTransparency*100 << " %" << G4endl;
      G4double dAbsorptionLength =  ((G4double)GetGeometryParameter("GridMeshThickness"))/(-log(dTransparency)); 
      G4MaterialPropertiesTable *pGateMeshPropertiesTable = pMeshMaterial->GetMaterialPropertiesTable();
      const G4int iNbEntries = 3;
      G4double pdMeshPhotonMomentum[iNbEntries] = {6.91*eV, 6.98*eV, 7.05*eV};
      G4double pdMeshAbsorptionLength[iNbEntries] = {dAbsorptionLength, dAbsorptionLength, dAbsorptionLength};
      pGateMeshPropertiesTable->RemoveProperty("ABSLENGTH");
      pGateMeshPropertiesTable->AddProperty("ABSLENGTH", pdMeshPhotonMomentum, pdMeshAbsorptionLength, iNbEntries);
   }
  else
    {
      G4cout << "ls!> GridMeshSS316LSteelLXe not found!" << G4endl;
      exit(-1);
    }
}

void muensterTPCDetectorConstruction::SetGXeAbsorbtionLength(G4double dAbsorbtionLength)
{
  G4Material *pGXeMaterial = G4Material::GetMaterial(G4String("GXe"));
  
  if(pGXeMaterial)
    {
      G4cout << "----> Setting GXe absorbtion length to " << dAbsorbtionLength/m << " m" << G4endl;
    
      G4MaterialPropertiesTable *pGXePropertiesTable = pGXeMaterial->GetMaterialPropertiesTable();
    
      const G4int iNbEntries = 3;
    
      G4double GXe_PP[iNbEntries] = {6.91*eV, 6.98*eV, 7.05*eV};
      G4double GXe_ABSL[iNbEntries] = {dAbsorbtionLength, dAbsorbtionLength, dAbsorbtionLength};
      pGXePropertiesTable->RemoveProperty("ABSLENGTH");
      pGXePropertiesTable->AddProperty("ABSLENGTH", GXe_PP, GXe_ABSL, iNbEntries);
    }
  else
    {
      G4cout << "ls!> GXe materials not found!" << G4endl;
      exit(-1);
    }
}

void muensterTPCDetectorConstruction::SetLXeRefractionIndex(G4double dRefractionIndex)
{
  G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));
  
  if(pLXeMaterial)
    {
      G4cout << "----> Setting LXe refraction index to " << dRefractionIndex << G4endl;
    
      G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
	
      const G4int iNbEntries = 3;
    
      G4double LXe_PP[iNbEntries] = {6.91*eV, 6.98*eV, 7.05*eV};
      G4double LXe_RI[iNbEntries] = {dRefractionIndex, dRefractionIndex, dRefractionIndex};
      pLXePropertiesTable->RemoveProperty("RINDEX");
	  pLXePropertiesTable->AddProperty("RINDEX", LXe_PP, LXe_RI, iNbEntries);
    }
  else
    {
      G4cout << "ls!> LXe materials not found!" << G4endl;
      exit(-1);
    }
}

void muensterTPCDetectorConstruction::SetGXeTeflonReflectivity(G4double dGXeReflectivity)
{
  G4Material *pGXeTeflonMaterial = G4Material::GetMaterial(G4String("GXeTeflon"));
  
  if(pGXeTeflonMaterial)
    {
      G4cout << "----> Setting GXe Teflon reflectivity to " << dGXeReflectivity << G4endl;
    
      G4MaterialPropertiesTable *pGXeTeflonPropertiesTable = pGXeTeflonMaterial->GetMaterialPropertiesTable();
    
      const G4int iNbEntries = 3;
    
      G4double teflon_PP[iNbEntries] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
      G4double teflon_REFL[iNbEntries] = {dGXeReflectivity, dGXeReflectivity, dGXeReflectivity};
      pGXeTeflonPropertiesTable->RemoveProperty("REFLECTIVITY");
      pGXeTeflonPropertiesTable->AddProperty("REFLECTIVITY", teflon_PP, teflon_REFL, iNbEntries);
      //pTeflonPropertiesTable->DumpTable();
   }
  else
    {
      G4cout << "ls!> Teflon material not found!" << G4endl;
      exit(-1);
    }
}

//******************************************************************/
// UpdateGeometry
//******************************************************************/
void muensterTPCDetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


