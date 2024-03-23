//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file FlashDetectorConstruction.cc
/// \brief Implementation of the FlashDetectorConstruction class


#include "FlashDetectorConstruction.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "G4SDManager.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"

#include "Applicator.hh"


#include "G4MaterialPropertiesTable.hh"

#include "G4PSEnergyDeposit.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VisAttributes.hh"
#include "FlashDetectorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::FlashDetectorConstruction()
    : G4VUserDetectorConstruction(), physicalTreatmentRoom(0),logicTreatmentRoom(0), Collimator(0), fPhantom(0), fPhantomLogicalVolume(0), fPhant_phys(0), fCheckOverlaps(true) {

    DefineMaterials();
    fDetectorMessenger = new FlashDetectorMessenger(this);

    SetPhantomSize(100. *cm, 100. *cm, 100. *cm);
    SetAirGap(0*cm); // Set the air gap between the water phantom and the end of the applicator
    SetDetectorThickness(10*um);
    SetDetector_subThickness(370*um);
    SetDetectorWidth(5*mm);
    SetAirGap_water_detector(0*cm); // Set the air gap between the end of the water phantom and the entrance of the detector
}


FlashDetectorConstruction::~FlashDetectorConstruction() {

    delete fDetectorMessenger;
}


void FlashDetectorConstruction::DefineMaterials() {
    nist = G4NistManager::Instance();
    //write here a function to define custom materials
    G4bool isotopes = false;
    Si = nist->FindOrBuildElement("Si", isotopes);
    C = nist->FindOrBuildElement("C", isotopes);

  }


G4VPhysicalVolume *
FlashDetectorConstruction::ConstructPhantom(G4double CollPos) {
    //This function creates a cubic phantom with the point Collpos on the surface of the cube.

    fPhantomMaterial = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");//(EJ200


    std::vector<G4double> energy     = {0.49195596066714914 * eV, 0.4890084978884109 * eV, 0.48623128113710545 * eV, 0.483807713078551 * eV, 0.4820654122719176 * eV, 0.48088894414499217 * eV, 0.48003524673921893 * eV, 0.4785898125980352 * eV, 0.4792718279591098 * eV, 0.47780040336958346 * eV, 0.477036340757137 * eV, 0.4761936051791291 * eV, 0.4767787354207054 * eV, 0.47555089522562244 * eV, 0.47498126256188955 * eV, 0.4744727467360512 * eV, 0.47322038252337034 * eV, 0.4738008478505348 * eV, 0.47171360936683193 * eV, 0.4699727964917785 * eV, 0.4681554086365279 * eV, 0.4657899251706182 * eV, 0.46325050725853345 * eV, 0.4620432658977899 * eV, 0.4610700170427911 * eV, 0.4599122249027489 * eV, 0.45823526156414973 * eV, 0.45618471582746395 * eV, 0.45367321338727334 * eV, 0.4512818724373308 * eV, 0.4489156091281936 * eV, 0.4468918977104331 * eV, 0.44425675389476715 * eV, 0.4419634013542076 * eV, 0.43846715233480865 * eV, 0.4372439045475469 * eV, 0.43552515056067403 * eV, 0.433719960440586 * eV, 0.4318397097023332 * eV, 0.42985959485490316 * eV, 0.4278524548553194 * eV, 0.4253937528097755 * eV, 0.42329055480951455 * eV, 0.42120805142155404 * eV, 0.4191459387013957 * eV, 0.41710391862764096 * eV, 0.4150816989584075 * eV, 0.41350050550466466 * eV, 0.411634575621106 * eV, 0.40982667685077445 * eV, 0.4077140169332079 * eV, 0.40578161048269684 * eV, 0.4037580786132559 * eV, 0.4020030039006597 * eV, 0.40009275659344784 * eV, 0.3982317536941048 * eV, 0.39638798328337305 * eV, 0.3945612071105923 * eV, 0.3928126229219143 * eV, 0.39095770623551396 * eV, 0.38918052649443385 * eV, 0.387419430722183 * eV, 0.3856742015558006 * eV, 0.38394462553142794 * eV, 0.3822304929972698 * eV, 0.38053159802887787 * eV, 0.37907647753930257 * eV, 0.44012186115214175 * eV};

    std::vector<G4double> rindex     = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58};
        
    std::vector<G4double> absorption = {250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm};
        
    std::vector<G4double> scint_spectrum = {0.05823104, 0.08459519, 0.12368447, 0.16171644, 0.2004555 , 0.23829181, 0.28591094, 0.39202739, 0.33083689, 0.45293744, 0.50300763, 0.59858826, 0.54285467, 0.6623958 , 0.7095871 , 0.74558667, 0.82240405, 0.78620626, 0.88281501, 0.93160535, 0.96721042, 0.99279729, 0.98753113, 0.94173171, 0.89653785, 0.85309237, 0.79523671, 0.75072293, 0.7165109 , 0.70802769, 0.70453461, 0.70386926, 0.67891865, 0.64947693, 0.57685102, 0.53492592, 0.4923089 , 0.45173506, 0.41767365, 0.38724143, 0.35779866, 0.32811311, 0.30682193, 0.28819214, 0.27072672, 0.25475833, 0.2386236 , 0.22649603, 0.21389315, 0.20329625, 0.19105111, 0.18188671, 0.16754422, 0.15650855, 0.14564101, 0.1358271 , 0.12767657, 0.11902703, 0.1121864 , 0.10522102, 0.09856753, 0.09357741, 0.0877556 , 0.08259914, 0.07677733, 0.07278523, 0.06712976, 0.62144578};


        G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
        MPT->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV);
        MPT->AddProperty("RINDEX", energy, rindex);
        MPT->AddProperty("ABSLENGTH", energy, absorption);
        MPT-> AddProperty("SCINTILLATIONCOMPONENT1", energy, scint_spectrum);
        MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
        MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20.*ns);
        MPT->AddConstProperty("SCINTILLATIONYIELD1", 1.0);

        fPhantomMaterial->SetMaterialPropertiesTable(MPT);
        G4cout << "Phantom G4MaterialPropertiesTable:" << G4endl;
        MPT->DumpTable();




    fPosition_coefficient = CollPos;
    
    fPhantom_coordinateX = (fPosition_coefficient * mm + fPhantomSizeX / 2);

    fPhantomPosition =  G4ThreeVector(fPhantom_coordinateX, 0. * mm, 0. * mm); //phantom is constructed with the entrance surface attached to the applicator 
    

    // Definition of the solid volume of the Phantom
    fPhantom = new G4Box("Phantom", fPhantomSizeX / 2, fPhantomSizeY / 2, fPhantomSizeZ / 2);

    // Definition of the logical volume of the Phantom
    fPhantomLogicalVolume =
        new G4LogicalVolume(fPhantom, fPhantomMaterial, "phantomLog", 0, 0, 0);

    // Definition of the physical volume of the Phantom
    fPhant_phys =
        new G4PVPlacement(0, fPhantomPosition, "phantomPhys", fPhantomLogicalVolume,
                            physicalTreatmentRoom, false, 0);
    //define the region to set cuts in FlashPhysicsList.cc and step limit
    G4Region *PhantomRegion = new G4Region("Phantom_reg");
    fPhantomLogicalVolume->SetRegion(PhantomRegion);
    PhantomRegion->AddRootLogicalVolume(fPhantomLogicalVolume);

    // Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
    red->SetVisibility(true);

    blue = new G4VisAttributes(G4Colour(0 / 255., 0. / 255., 255. / 255.));
    blue->SetVisibility(true);

    fPhantomLogicalVolume->SetVisAttributes(red);
    //set step limit in phantom
    G4double maxStep = 0.1 * mm;
    fStepLimit = new G4UserLimits(maxStep);
    fPhantomLogicalVolume->SetUserLimits(fStepLimit);

    return fPhant_phys;
}


G4VPhysicalVolume *
FlashDetectorConstruction::ConstructDetector(){
    //Detector
    G4double fDensity_SiC=3.22*g/cm3;

    SiC=new G4Material("SiC", fDensity_SiC,2);
    SiC->AddElement(Si,1);
    SiC->AddElement(C,1);
    
    fDetectorMaterial=SiC;


    fDetectorPosition=fPhantom_coordinateX+fAirGap+fPhantomSizeX/2+fDet_thickness/2+fAirGap_phantom_det;
    
    fDet_box = new G4Box("Detector",fDet_thickness/2,fDet_width/2,fDet_width/2);
    
    // Definition of the logical volume of the Detector
    fDetLogicalVolume =
        new G4LogicalVolume(fDet_box, fDetectorMaterial, "DetectorLog", 0, 0, 0);
    fDet_phys = new G4PVPlacement(0,G4ThreeVector(fDetectorPosition, 0. * mm, 0. * mm), "DetPhys",fDetLogicalVolume,physicalTreatmentRoom,false, 0, fCheckOverlaps);

    
    fDet_sub = new G4Box("Det_sub",fDet_sub_thickness/2,fDet_width/2,fDet_width/2);
    
    // Definition of the logical volume of the Detector 
    fDet_sub_LogicalVolume =
        new G4LogicalVolume(fDet_sub, fDetectorMaterial, "Det_sub_Log", 0, 0, 0);
    fDet_sub_phys = new G4PVPlacement(0,G4ThreeVector(fDetectorPosition+fDet_thickness+fDet_sub_thickness/2, 0. * mm, 0. * mm), "Det_sub_Phys",fDet_sub_LogicalVolume,physicalTreatmentRoom,false, 0, fCheckOverlaps);


    return fDet_phys;

}


G4VPhysicalVolume *FlashDetectorConstruction::Construct() {
    // -----------------------------
    // Treatment room - World volume
    //------------------------------
    // Treatment room sizes 
    const G4double worldX = 400.0 * cm;
    const G4double worldY = 400.0 * cm;
    const G4double worldZ = 400.0 * cm;
    G4bool isotopes = false;
    // Filled with air
    airNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);

    G4Box *treatmentRoom = new G4Box("TreatmentRoom", worldX, worldY, worldZ);
    logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, airNist,
                                            "logicTreatmentRoom", 0, 0, 0);
    physicalTreatmentRoom =
        new G4PVPlacement(0, G4ThreeVector(), "physicalTreatmentRoom",
                            logicTreatmentRoom, 0, false, 0);

    // The treatment room is invisible in the Visualisation
    logicTreatmentRoom->SetVisAttributes(G4VisAttributes::GetInvisible());

    // -----------------------------
    // Applicator + phantom + Default dimensions
    //------------------------------
    Collimator = new Applicator(physicalTreatmentRoom);
    
    fPhantom_physical =
            ConstructPhantom(Collimator->fFinalApplicatorXPositionFlash +
    Collimator->fHightFinalApplicatorFlash+fAirGap);
    //ConstructDetector();
    
     return physicalTreatmentRoom;
}


void FlashDetectorConstruction::ConstructSDandField() {
//modify this function if you want to insert a sensitive detector
}


/////MESSANGER ///

G4bool FlashDetectorConstruction::SetPhantomMaterial(G4String material)
{

    if (G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	fPhantomMaterial  = pMat;

	if (fPhantomLogicalVolume) 
	{
	    
	    fPhantomLogicalVolume ->  SetMaterial(pMat);

	    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    G4cout << "The material of Phantom/Detector has been changed to " << material << G4endl;
	}
    }
    else
    {
	G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
	G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl; 
	return false;
    }

    return true;
}


void FlashDetectorConstruction::SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) fPhantomSizeX = sizeX;
    if (sizeY > 0.) fPhantomSizeY = sizeY;
    if (sizeZ > 0.) fPhantomSizeZ = sizeZ;
}


void FlashDetectorConstruction::SetAirGap(G4double displ)
{
   fAirGap=displ;
}


G4bool FlashDetectorConstruction::SetDetectorMaterial(G4String material)
{
    if (G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	fDetectorMaterial  = pMat;

	if (fDetLogicalVolume) 
	{
	    
	    fDetLogicalVolume ->  SetMaterial(pMat);

	    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    G4cout << "The material of Phantom/Detector has been changed to " << material << G4endl;
	}
    }
    else
    {
	G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
	G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl; 
	return false;
    }

    return true;
}




void FlashDetectorConstruction::SetDetectorThickness(G4double thickness)
{
   fDet_thickness=thickness;
}


void FlashDetectorConstruction::SetDetectorWidth(G4double width)
{
   fDet_width=width;
}


void FlashDetectorConstruction::SetDetector_subThickness(G4double thickness_sub)
{
    fDet_sub_thickness= thickness_sub;
}


void FlashDetectorConstruction::SetAirGap_water_detector(G4double spost)
{
  fAirGap_phantom_det=spost;
}


void FlashDetectorConstruction::SetDetectorPosition(G4double position)
{
   fDetectorPosition=position;
}
