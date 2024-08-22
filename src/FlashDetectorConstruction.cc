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
#include "G4SubtractionSolid.hh"
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

    SetPhantomSize(10. *cm, 10. *cm, 10. *cm);
    SetPinholeDistance(20. *cm);
    SetDetectorDistance(20.*cm);

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
    
    std::vector<G4double> energy     = {00.5 * eV, 0.35 *eV};
    std::vector<G4double> rindex     = {1.58, 1.58};
    std::vector<G4double> absorption = {250.*cm, 250.*cm};
    std::vector<G4double> scint_spectrum = {0.5, 0.5};

        G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
        MPT->AddConstProperty("SCINTILLATIONYIELD", 1000./MeV);
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
FlashDetectorConstruction::ConstructPinhole() {
    //This function creates a ....

    //Material properties: 
    PinholeMaterial = nist->FindOrBuildMaterial("G4_AIR");

    // Outer and inner cilinder dimensions 
    G4double outerRadius = 0.5*m;
    G4double height = 0.5*mm;
    G4double startAngle = 0.*deg;
    G4double spanningAngle = 360.*deg;
    G4double innerRadius = 1.*mm;  
    G4Tubs* outerCylinder = new G4Tubs("OuterCylinder", 0, outerRadius, height/2, startAngle, spanningAngle);
    G4Tubs* innerCylinder = new G4Tubs("InnerCylinder", 0, innerRadius, height/2, startAngle, spanningAngle);
    G4SubtractionSolid* PinholeCilinder = new G4SubtractionSolid("Pinhole", outerCylinder, innerCylinder);
    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateY(90.*deg); // Ruota di 90 gradi attorno all'asse Y

    // Definition of the logical volume of the Pinhole
    PinholeLogicalVolume = new G4LogicalVolume(PinholeCilinder, PinholeMaterial, "pinholeLog", 0, 0, 0);
    PinholePosition =  G4ThreeVector(fPhantomSizeX + PinholeDistance, 0. * mm, 0. * mm); 
    // Definition of the physical volume of the Pinhole
    Pihole_phys = new G4PVPlacement(rotationMatrix, PinholePosition, "pinholePhys", PinholeLogicalVolume, physicalTreatmentRoom, false, 0);
    
    // Visualisation attributes of the phantom
    gray = new G4VisAttributes(G4Colour(211 / 255., 211 / 255., 211 / 255.));
    gray->SetVisibility(true);
    PinholeLogicalVolume->SetVisAttributes(gray);

    return Pihole_phys;

}


/*G4VPhysicalVolume *
FlashDetectorConstruction::ConstructDetector() {
    //This function creates a ....

    //Material properties: 
   // DetectorMaterial = nist->FindOrBuildMaterial("G4_AIR");


    // Definition of the logical volume of the Pinhole
    DetectorLogicalVolume = new G4LogicalVolume(DetectorSheet, PinholeMaterial, "pinholeLog", 0, 0, 0);
    DetectorPosition =  G4ThreeVector(fPhantomSizeX + PinholeDistance, 0. * mm, 0. * mm); 
    // Definition of the physical volume of the Pinhole
    Pihole_phys = new G4PVPlacement(0, DetectorPosition, "detectorPhys", DetectorLogicalVolume, physicalTreatmentRoom, false, 0);
    
    // Visualisation attributes of the phantom
    gray = new G4VisAttributes(G4Colour(211 / 255., 211 / 255., 211 / 255.));
    gray->SetVisibility(true);
    DetectorLogicalVolume->SetVisAttributes(gray);

    return Detector_phys;

}*/



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
    std::vector<G4double> energy     = {00.5 * eV, 0.35 *eV};
    std::vector<G4double> rindex     = {1., 1.};
    G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
    MPT->AddProperty("RINDEX", energy, rindex);
    airNist->SetMaterialPropertiesTable(MPT);
    G4cout << "Air G4MaterialPropertiesTable:" << G4endl;
    MPT->DumpTable();

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



    // -----------------------------
    // Pinhole camera
    //------------------------------

    Pihole_physical = ConstructPinhole();

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


void FlashDetectorConstruction::SetPinholeDistance(G4double distanceX)
{
    if (distanceX > 0.) PinholeDistance = distanceX;

}

void FlashDetectorConstruction::SetDetectorDistance(G4double distanceX)
{
    if (distanceX > 0.) DetectorDistance = distanceX;

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
