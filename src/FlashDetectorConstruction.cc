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
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
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

    SetAirGap(0*cm); // Set the air gap between the water phantom and the end of the applicator
    SetPhantomSize(10. *cm, 10. *cm, 10. *cm);
    SetPinholeDistance(25. *cm);
    SetDetectorDistance(10.*cm); // Set the air gap between the water phantom and the end of the detector

    SetDetectorThickness(1*mm);
    SetDetectorWidth(15*cm);
        
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


G4VPhysicalVolume *FlashDetectorConstruction::ConstructPhantom(G4double CollPos) {
    //This function creates a cubic phantom with the point Collpos on the surface of the cube.

    fPhantomMaterial = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");//(EJ200
    
    std::vector<G4double> energy     = {2.48 * eV, 3.1 * eV}; //sto generando solo tra i 400 e i 500 nm. lambda [nm] = 1240/E[eV]
    std::vector<G4double> rindex     = {1.58, 1.58};
    std::vector<G4double> absorption = {380.*cm, 380.*cm};
    std::vector<G4double> scint_spectrum = {0.5, 0.5};

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
    fPhantomLogicalVolume = new G4LogicalVolume(fPhantom, fPhantomMaterial, "phantomLog", 0, 0, 0);

    // Definition of the physical volume of the Phantom
    fPhant_phys = new G4PVPlacement(0, fPhantomPosition, "phantomPhys", fPhantomLogicalVolume, physicalTreatmentRoom, false, 0);
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


    // Creazione e configurazione della superficie ottica
    G4OpticalSurface* opticalSurface = new G4OpticalSurface("OpticalSurface");
    opticalSurface->SetType(dielectric_dielectric);  // Tipo di superficie
    opticalSurface->SetFinish(ground);  // Finitura opaca
    opticalSurface->SetModel(unified);  // Modello di riflessione e trasmissione
    opticalSurface->SetPolish(0.0);  // Coefficiente di riflessione a zero

    // Associa la superficie ottica ai confini
    new G4LogicalBorderSurface("BorderSurface", fPhant_phys, physicalTreatmentRoom, opticalSurface);
    return fPhant_phys;
}


G4VPhysicalVolume *FlashDetectorConstruction::ConstructPinhole() {
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

    G4ThreeVector PinholePosition = G4ThreeVector(fPhantom_coordinateX + fPhantomSizeX/2 + fDet_thickness/2 + PinholeDistance, 0. * mm, 0. * mm); //sicura di fDet_thickness/2?!

    // Definition of the physical volume of the Pinhole
    Pihole_phys = new G4PVPlacement(rotationMatrix, PinholePosition, "pinholePhys", PinholeLogicalVolume, physicalTreatmentRoom, false, 0);
    
    // Visualisation attributes of the pinhole
    gray = new G4VisAttributes(G4Colour(211 / 255., 211 / 255., 211 / 255.));
    gray->SetVisibility(true);
    PinholeLogicalVolume->SetVisAttributes(gray);

    return Pihole_phys;

}



G4VPhysicalVolume *FlashDetectorConstruction::ConstructDetector(){
    //Detector

    //poi andrà cambiato    DetectorMaterial = nist->FindOrBuildMaterial("G4_AIR");
    G4double fDensity_SiC=3.22*g/cm3;
    SiC=new G4Material("SiC", fDensity_SiC,2);
    SiC->AddElement(Si,1);
    SiC->AddElement(C,1);
    DetectorMaterial=SiC;
    
    Det_box = new G4Box("Detector", fDet_thickness/2, fDet_width/2,fDet_width/2);
    fDetectorPosition_t = fPhantomSizeX + DetectorDistance + fDet_thickness/2;
    fDetectorPosition_l = fPhantomSizeX*0.5 + DetectorDistance + fDet_thickness/2;

    G4int offset = fPhantomSizeX * 0.5 ; 

    G4RotationMatrix* rotationMatrix_z = new G4RotationMatrix();
    rotationMatrix_z->rotateZ(90.*deg); // Ruota di 90 gradi attorno all'asse Z
    G4RotationMatrix* rotationMatrix_y = new G4RotationMatrix();
    rotationMatrix_y->rotateY(90.*deg); // Ruota di 90 gradi attorno all'asse Y
    

    // Definition of the logical volume of the Detector
    fDetLogicalVolume = new G4LogicalVolume(Det_box, DetectorMaterial, "DetectorLog", 0, 0, 0);
    fDet_phys1 = new G4PVPlacement(0,G4ThreeVector(fDetectorPosition_t, 0., 0.), "DetPhys",fDetLogicalVolume,physicalTreatmentRoom,false, 0, fCheckOverlaps);
    fDet_phys2 = new G4PVPlacement(rotationMatrix_z,G4ThreeVector(offset, fDetectorPosition_l, 0.), "DetPhys",fDetLogicalVolume,physicalTreatmentRoom,false, 0, fCheckOverlaps);
    fDet_phys3 = new G4PVPlacement(rotationMatrix_z,G4ThreeVector(offset, -fDetectorPosition_l, 0.), "DetPhys",fDetLogicalVolume,physicalTreatmentRoom,false, 0, fCheckOverlaps);
    fDet_phys4 = new G4PVPlacement(rotationMatrix_y, G4ThreeVector(offset, 0., fDetectorPosition_l), "DetPhys", fDetLogicalVolume, physicalTreatmentRoom, false, 0, fCheckOverlaps);
    fDet_phys5 = new G4PVPlacement(rotationMatrix_y, G4ThreeVector(offset, 0., -fDetectorPosition_l), "DetPhys", fDetLogicalVolume, physicalTreatmentRoom, false, 0, fCheckOverlaps);
    

    // Visualisation attributes of the detector
    gray = new G4VisAttributes(G4Colour(211 / 255., 211 / 255., 211 / 255.));
    gray->SetVisibility(true);
    fDetLogicalVolume->SetVisAttributes(gray);

    return fDet_phys1, fDet_phys2, fDet_phys3, fDet_phys4, fDet_phys5;

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
    std::vector<G4double> energy     = {2.48 * eV, 3.1 * eV};

    //std::vector<G4double> energy     = {00.5 * eV, 0.35 *eV};
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
    fPhantom_physical = ConstructPhantom(Collimator->fFinalApplicatorXPositionFlash +
    Collimator->fHightFinalApplicatorFlash+fAirGap);

    // -----------------------------
    // Pinhole camera
    //------------------------------
    //Pihole_physical = ConstructPinhole();

    // -----------------------------
    // Detector pannel
    //------------------------------
    ConstructDetector();
    
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
	DetectorMaterial  = pMat;

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




