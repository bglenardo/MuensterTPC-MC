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
#ifndef __muensterTPCPDETECTORCONSTRUCTION_H__
#define __muensterTPCPDETECTORCONSTRUCTION_H__

#include <globals.hh>

#include <vector>
#include <map>

using std::vector;
using std::map;

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class muensterTPCDetectorMessenger;

#include <G4VUserDetectorConstruction.hh>

class muensterTPCDetectorConstruction: public G4VUserDetectorConstruction {
public:
	muensterTPCDetectorConstruction();
	~muensterTPCDetectorConstruction();

	G4VPhysicalVolume* Construct();

	void SetTeflonReflectivity(G4double dReflectivity);
	void SetGXeTeflonReflectivity(G4double dGXeReflectivity);
	void SetSS304LSteelReflectivity(G4double dReflectivity);
	void SetLXeScintillation(G4bool dScintillation);
	void SetLXeLevel(G4double dlevel);
	void SetMaterial (const G4String&);
	void SetLXeMeshMaterial (const G4String&);
	void SetGXeMeshMaterial (const G4String&);
	void SetLXeAbsorbtionLength(G4double dAbsorbtionLength);
	void SetGXeAbsorbtionLength(G4double dAbsorbtionLength);
	void SetLXeRayScatterLength(G4double dRayScatterLength);
	void SetLXeRefractionIndex(G4double dRefractionIndex);
	
	void SetLXeMeshTransparency(G4double dTransparency); 
	void SetGXeMeshTransparency(G4double dTransparency); 

	static G4double GetGeometryParameter(const char *szParameter);

public:
	G4Material*	GetMaterial()	{return LXeMaterial;};

private:
	void DefineMaterials();
	void DefineGeometryParameters();
	void UpdateGeometry();

	G4Material* LXeMaterial;
	G4UserLimits* fStepLimit;       // pointer to user step limits

	void ConstructLaboratory();
	void ConstructShield();
	void ConstructXenon();
	void ConstructLiquidXenonCylinder();
	void ConstructBell();
	void ConstructFieldCage();
	void ConstructTPC();
	void ConstructPmtSupports();
	void ConstructPmtArrays();
	void ConstructCryostat();
	void ConstructInnerCryostat();
	void ConstructOuterCryostat();

	void PrintGeometryInformation();
	void PrintPhysicalVolumes();

  typedef enum {PMT_CASING, PMT_WINDOW, PMT_BASE, PMT_SPACE} PmtPart;

	G4ThreeVector GetPmtPosition(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionTopArray(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionBottomArray(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionTopVetoArray(G4int iPmtNb, PmtPart ePmtPart);
	G4ThreeVector GetPmtPositionBottomVetoArray(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotation(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotationTopArray(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotationTopVetoArray(G4int iPmtNb, PmtPart ePmtPart);
	G4RotationMatrix *GetPmtRotationBottomVetoArray(G4int iPmtNb, PmtPart ePmtPart);

	vector<G4int> ComputePmtPartition(G4double dRadius, G4double dHeightFraction);
	G4ThreeVector ComputeXYPmtPositionForGridPattern(G4int iPmtNb, const vector<G4int> &hPmtsPerRow);
	G4ThreeVector ComputeXYPmtPositionForRingPattern(G4int iPmtNb, G4int iNbPmts, G4double dRingRadius);

private:
	// rotation matrices
	G4RotationMatrix *m_pRotationXPlus225;
	G4RotationMatrix *m_pRotationXMinus225;
	G4RotationMatrix *m_pRotationXPlus45;
	G4RotationMatrix *m_pRotationXMinus45;
	G4RotationMatrix *m_pRotationXPlus90;
	G4RotationMatrix *m_pRotationXMinus90;
	G4RotationMatrix *m_pRotationX180;
	G4RotationMatrix *m_pRotationX0;

	G4RotationMatrix *m_pRotationZTeflonLedBlock;

	// logical volumes
	G4LogicalVolume *m_pMotherLogicalVolume;
	G4LogicalVolume *m_pLabLogicalVolume;
	G4LogicalVolume *m_pShieldPolishLeadLogicalVolume;
	G4LogicalVolume *m_pShieldFrenchLeadLogicalVolume;
	G4LogicalVolume *m_pShieldPolyethyleneLogicalVolume;
	G4LogicalVolume *m_pShieldCopperLogicalVolume;
	G4LogicalVolume *m_pShieldCavityLogicalVolume;
	G4LogicalVolume *m_pShieldSteelBarLogicalVolume;

  //mystuff
	G4LogicalVolume *m_pLXeLogicalVolume;
	G4LogicalVolume *m_pGXeLogicalVolume;
	G4LogicalVolume *m_pVetoGXeLogicalVolume;
	G4LogicalVolume *m_pLidInTubeGXeLogicalVolume;
	G4LogicalVolume *m_pInnerCryostatVesselLogicalVolume;
 	G4LogicalVolume *m_pInnerCryostatBoreFlangeLogicalVolume;
	G4LogicalVolume *m_pInnerCryostatBlindFlangeLogicalVolume;
	G4LogicalVolume *m_pInnerCryostatCopperRingsLogicalVolume;
	G4LogicalVolume *m_pOuterCryostatVesselLogicalVolume;
	G4LogicalVolume *m_pOuterCryostatVacuumLogicalVolume;

	G4LogicalVolume *m_pPTFEInnerCylinderLogicalVolume;
	G4LogicalVolume *m_pPTFEOuterCylinderLogicalVolume;
	G4LogicalVolume *m_pPTFEBottomPMTHolderLogicalVolume;
	G4LogicalVolume *m_pPTFETopPMTHolderLogicalVolume;
	G4LogicalVolume *m_pBottomPTFESlabLogicalVolume;
	G4LogicalVolume *m_pTopPTFESlabLogicalVolume;
	G4LogicalVolume *m_pCopperRingsLogicalVolume;	
	G4LogicalVolume *m_pCathodeGridMeshSupportLogicalVolume;
	G4LogicalVolume *m_pGridMeshSupportLogicalVolume;
	G4LogicalVolume *m_pGridMeshLXeLogicalVolume;
	G4LogicalVolume *m_pGridMeshGXeLogicalVolume;
	G4LogicalVolume *m_pGridMeshBorderLogicalVolume;
	G4LogicalVolume *m_pBottomSteelRingLogicalVolume;
	G4LogicalVolume *m_pOuterCryostatBoreFlangeLogicalVolume;
	G4LogicalVolume *m_pOuterCryostatBlindFlangeLogicalVolume;
	G4LogicalVolume *m_pSpaceBelowTopPMTLogicalVolume;

	G4LogicalVolume *m_pBellLogicalVolume;	
	G4LogicalVolume *m_pTopPmtTeflonHolderLogicalVolume;
	G4LogicalVolume *m_pBellPmtSupportRingLogicalVolume;
	G4LogicalVolume *m_pBellSupportCylinderInLXeLogicalVolume;
	G4LogicalVolume *m_pBellSupportCylinderInVetoGXe1LogicalVolume;
	G4LogicalVolume *m_pBellSupportCylinderInVetoGXe2LogicalVolume;
	G4LogicalVolume *m_pTopVetoAngleLogicalVolume;
	G4LogicalVolume *m_pUpperSideVetoAngleLogicalVolume;

	G4LogicalVolume *m_pTopGridMeshLogicalVolume;
	G4LogicalVolume *m_pTopGridRingLogicalVolume;
	G4LogicalVolume *m_pAnodeGridMeshLogicalVolume;
	G4LogicalVolume *m_pAnodeGridRingLogicalVolume;
	G4LogicalVolume *m_pBottomGridMeshLogicalVolume;
	G4LogicalVolume *m_pBottomGridRingLogicalVolume;
	G4LogicalVolume *m_pCathodeGridMeshLogicalVolume;
	G4LogicalVolume *m_pCathodeGridRingLogicalVolume;
	G4LogicalVolume *m_pScreenGridRingLogicalVolume;
	G4LogicalVolume *m_pScreenGridMeshLogicalVolume;
	G4LogicalVolume *m_pFieldShapingWireRingsLogicalVolume;

	G4LogicalVolume *m_pGridHolderInLXeLogicalVolume;
	G4LogicalVolume *m_pGridHolderInGXeLogicalVolume;

	G4LogicalVolume *m_pCopperScreenLogicalVolume;
	G4LogicalVolume *m_pTopPlateLogicalVolume;	
	G4LogicalVolume *m_pTeflonPanelLogicalVolume;
	G4LogicalVolume *m_pTeflonRodsLogicalVolume;
	G4LogicalVolume *m_pBottomPlateLogicalVolume;	
	G4LogicalVolume *m_pUpperBasePlateLogicalVolume;	
	G4LogicalVolume *m_pLowerTeflonPanelLogicalVolume;
	G4LogicalVolume *m_pLowerBasePlateLogicalVolume;	
	G4LogicalVolume *m_pBottomPmtPlateLogicalVolume;	
	G4LogicalVolume *m_pBottomVetoAngleLogicalVolume;
	G4LogicalVolume *m_pLowerSideVetoAngleLogicalVolume;
	G4LogicalVolume *m_pTeflonSideVetoLiningLogicalVolume;
	G4LogicalVolume *m_pTeflonLedBlockLogicalVolume;
	G4LogicalVolume *m_pTeflonBottomReflectorTopLogicalVolume;
	G4LogicalVolume *m_pTeflonBottomReflectorBottomLogicalVolume;
	G4LogicalVolume *m_pBottomReflectorHolderLogicalVolume;

	G4LogicalVolume *m_pPmtWindowLogicalVolume;
	G4LogicalVolume *m_pPmtCasingLogicalVolume;
	G4LogicalVolume *m_pPmtInteriorLogicalVolume;
	G4LogicalVolume *m_pPmtPhotoCathodeLogicalVolume;
	G4LogicalVolume *m_pPmtBaseLogicalVolume;

	G4LogicalVolume *m_pCryostatVesselFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselInnerWallLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselInnerWallDomeLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselOuterWallLogicalVolume;
	G4LogicalVolume *m_pCryostatVesselOuterWallDomeLogicalVolume;
	G4LogicalVolume *m_pVacuumCryostatVesselLogicalVolume;

	G4LogicalVolume *m_pCryostatLidFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInnerWallLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInnerWallDomeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOuterWallLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOuterWallDomeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeInsideVetoGXeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeInsideVetoLXeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeElbowLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeTopLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeTopFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeTopToShieldLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeJacketLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeJacketCoverLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeOuterJacketLogicalVolume;
	G4LogicalVolume *m_pCryostatLidInTubeOuterJacketCoverLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeElbowLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeTopLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeTopFlangeLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeTopToShieldLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeJacketLogicalVolume;
	G4LogicalVolume *m_pCryostatLidOutTubeJacketCoverLogicalVolume;

	G4LogicalVolume *m_pCryostatLidLogicalVolume;
	G4LogicalVolume *m_pVacuumCryostatLidLogicalVolume;

        G4LogicalVolume *m_pOuterShellLogicalVolume;


  // physical volumes
	G4VPhysicalVolume *m_pLabPhysicalVolume;

	G4VPhysicalVolume *m_pShieldPolishLeadPhysicalVolume;
	G4VPhysicalVolume *m_pShieldFrenchLeadPhysicalVolume;
	G4VPhysicalVolume *m_pShieldPolyethylenePhysicalVolume;
	G4VPhysicalVolume *m_pShieldCopperPhysicalVolume;
	G4VPhysicalVolume *m_pShieldCavityPhysicalVolume;
	G4VPhysicalVolume *m_pShieldSteelBarPhysicalVolume;

  //my stuff
	G4VPhysicalVolume *m_pLXePhysicalVolume;
        G4VPhysicalVolume *m_pOuterShellPhysicalVolume;
	G4VPhysicalVolume *m_pGXePhysicalVolume;
	G4VPhysicalVolume *m_pVetoGXePhysicalVolume;
	G4VPhysicalVolume *m_pLidInTubeGXePhysicalVolume;

	G4VPhysicalVolume *m_pPTFEInnerCylinderPhysicalVolume;
	G4VPhysicalVolume *m_pPTFEOuterCylinderPhysicalVolume;
	G4VPhysicalVolume *m_pPTFEBottomPMTHolderPhysicalVolume;
	G4VPhysicalVolume *m_pPTFETopPMTHolderPhysicalVolume;
	G4VPhysicalVolume *m_pBottomPTFESlabPhysicalVolume;
	G4VPhysicalVolume *m_pTopPTFESlabPhysicalVolume;
	G4VPhysicalVolume *m_pCopperRingsPhysicalVolume;
	G4VPhysicalVolume *m_pInnerCryostatVesselPhysicalVolume;
	G4VPhysicalVolume *m_pBottomSteelRingPhysicalVolume;
	G4VPhysicalVolume *m_pOuterCryostatVesselPhysicalVolume;
	G4VPhysicalVolume *m_pOuterCryostatVacuumPhysicalVolume;

	G4VPhysicalVolume *m_pCathodeGridMeshSupportPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hCopperRingsPhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hGridMeshSupportPhysicalVolumes;	
	vector<G4VPhysicalVolume *> m_hGridMeshPhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hGridMeshBorderPhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hInnerCryostatBoreFlangePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hInnerCryostatBlindFlangePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hInnerCryostatCopperRingsPhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hOuterCryostatBoreFlangePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hOuterCryostatBlindFlangePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hSpaceBelowTopPMTPhysicalVolumes;

	G4VPhysicalVolume *m_pBellPhysicalVolume;
	G4VPhysicalVolume *m_pTopPmtTeflonHolderPhysicalVolume;
	G4VPhysicalVolume *m_pBellPmtSupportRingPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hTopVetoAnglePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hUpperSideVetoAnglePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hBellSupportCylinderInLXePhysicalVolumes;
	G4VPhysicalVolume *m_pBellSupportCylinderInVetoGXe1PhysicalVolume;
	G4VPhysicalVolume *m_pBellSupportCylinderInVetoGXe2PhysicalVolume;

	G4VPhysicalVolume *m_pTopGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pTopGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pAnodeGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pAnodeGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pBottomGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pBottomGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pCathodeGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pCathodeGridRingPhysicalVolume;	
	G4VPhysicalVolume *m_pScreenGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pScreenGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pGridHolderInLXePhysicalVolume;
	G4VPhysicalVolume *m_pGridHolderInGXePhysicalVolume;

	vector<G4VPhysicalVolume *> m_hFieldShapingWireRingsPhysicalVolumes;

	G4VPhysicalVolume *m_pCopperScreenPhysicalVolume;
	G4VPhysicalVolume *m_pTopPlatePhysicalVolume;
	G4VPhysicalVolume *m_pTeflonPanelPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hTeflonRodsPhysicalVolumes;
	G4VPhysicalVolume *m_pBottomPlatePhysicalVolume;
	G4VPhysicalVolume *m_pUpperBasePlatePhysicalVolume;
	G4VPhysicalVolume *m_pLowerTeflonPanelPhysicalVolume;
	G4VPhysicalVolume *m_pLowerBasePlatePhysicalVolume;
	G4VPhysicalVolume *m_pBottomPmtPlatePhysicalVolume;
	vector<G4VPhysicalVolume *> m_hBottomVetoAnglePhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hLowerSideVetoAnglePhysicalVolumes;
	G4VPhysicalVolume *m_pTeflonSideVetoLiningPhysicalVolume;

	vector<G4VPhysicalVolume *> m_hTeflonLedBlockPhysicalVolumes;

	vector<G4VPhysicalVolume *> m_hPmtWindowPhysicalVolumes;
	vector<G4VPhysicalVolume *> m_hPmtCasingPhysicalVolumes;
	G4VPhysicalVolume *m_pPmtInteriorPhysicalVolume;
	G4VPhysicalVolume *m_pPmtPhotoCathodePhysicalVolume;
	vector<G4VPhysicalVolume *> m_hPmtBasePhysicalVolumes;

	G4VPhysicalVolume *m_pCryostatVesselFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselInnerWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselInnerWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselOuterWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVesselOuterWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pVacuumCryostatVesselPhysicalVolume;

	G4VPhysicalVolume *m_pCryostatLidFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInnerWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInnerWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOuterWallPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOuterWallDomePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeInsideVetoGXePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeInsideVetoLXePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeElbowPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeTopPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeTopFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeTopToShieldPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeJacketPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeJacketCoverPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeOuterJacketPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidInTubeOuterJacketCoverPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeElbowPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeTopPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeTopFlangePhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeTopToShieldPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeJacketPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatLidOutTubeJacketCoverPhysicalVolume;

	G4VPhysicalVolume *m_pCryostatLidPhysicalVolume;
	G4VPhysicalVolume *m_pVacuumCryostatLidPhysicalVolume;

	vector<G4VPhysicalVolume *> m_hBottomReflectorHolderPhysicalVolumes;
	G4VPhysicalVolume *m_pTeflonBottomReflectorTopPhysicalVolume;
	G4VPhysicalVolume *m_pTeflonBottomReflectorBottomPhysicalVolume;
	
	static map<G4String, G4double> m_hGeometryParameters;
	
	muensterTPCDetectorMessenger *m_pDetectorMessenger;
};

#endif // __muensterTPCPDETECTORCONSTRUCTION_H__

