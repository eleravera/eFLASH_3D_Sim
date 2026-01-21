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
#include "G4Trd.hh"

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
#include "G4PhysicsOrderedFreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::FlashDetectorConstruction()
    : G4VUserDetectorConstruction(),logicTreatmentRoom(0), physicalTreatmentRoom(0), Collimator(0), fPhantom(0), fPhantomLogicalVolume(0), fPhant_phys(0), fCheckOverlaps(true) {
    /*
    This constructor initializes the FlashDetectorConstruction class, which is derived from `G4VUserDetectorConstruction`.
   It sets up the materials, the treatment room, and other parameters required for constructing the detector, phantom, 
   and collimator. It also sets the necessary physical and logical volumes, and prepares the detector messenger for 
   communication with the user interface.

   Initialization:
   1. `DefineMaterials()` is called to define the materials used for the phantom, collimator, and detector.
   2. A new instance of `FlashDetectorMessenger` is created (?)
   3. Physical parameters like air gap, phantom size, pinhole distance, and detector distance are set through setter methods.
   4. The detector's thickness and width are also initialized.

   Parameters:
   - AirGap: Air gap between the phantom and applicator or detector.
   - Phantom Size: Dimensions of the phantom (10 cm x 10 cm x 10 cm).
   - PinholeDistance: Distance for pinhole collimator setup (5 cm).
   - DetectorDistance: Distance between the phantom and the detector (5 cm).
   - DetectorThickness: Thickness of the detector (1 mm).
   - DetectorWidth: Width of the detector (15 cm).

   Notes:
   - The constructor ensures that the necessary parameters are initialized before the actual construction of the detector 
     and phantom. It also prepares the room and interaction setup (e.g., overlaps checking).
     */

    DefineMaterials();
    fDetectorMessenger = new FlashDetectorMessenger(this);

    SetAirGap(0.011*cm); // Set the air gap between the water phantom and the end of the applicator
    SetPhantomSize(10. *cm, 10. *cm, 10. *cm);
    SetPinholeDistance(20 *cm); // Set the air gap between the water phantom and the pinhole
    SetDetectorDistance(25*cm); // Set the air gap between the water phantom and the detector

}


FlashDetectorConstruction::~FlashDetectorConstruction() {
    delete fDetectorMessenger;
}


void FlashDetectorConstruction::DefineMaterials() {


    std::vector<G4double> energy_eV =  {2.48 * eV, 2.48554 * eV, 2.49111 * eV, 2.49671 * eV, 2.50232 * eV,
                                        2.50797 * eV, 2.51364 * eV, 2.51933 * eV, 2.52505 * eV, 2.5308 * eV,
                                        2.53657 * eV, 2.54237 * eV, 2.5482 * eV, 2.55405 * eV, 2.55993 * eV,
                                        2.56584 * eV, 2.57177 * eV, 2.57774 * eV, 2.58372 * eV, 2.58974 * eV,
                                        2.59579 * eV, 2.60186 * eV, 2.60796 * eV, 2.6141 * eV, 2.62026 * eV,
                                        2.62644 * eV, 2.63266 * eV, 2.63891 * eV, 2.64519 * eV, 2.6515 * eV,
                                        2.65783 * eV, 2.6642 * eV, 2.6706 * eV, 2.67703 * eV, 2.68349 * eV,
                                        2.68998 * eV, 2.6965 * eV, 2.70306 * eV, 2.70965 * eV, 2.71627 * eV,
                                        2.72292 * eV, 2.7296 * eV, 2.73632 * eV, 2.74307 * eV, 2.74985 * eV,
                                        2.75667 * eV, 2.76352 * eV, 2.77041 * eV, 2.77733 * eV, 2.78428 * eV,
                                        2.79127 * eV, 2.79829 * eV, 2.80535 * eV, 2.81245 * eV, 2.81958 * eV,
                                        2.82675 * eV, 2.83395 * eV, 2.84119 * eV, 2.84847 * eV, 2.85579 * eV,
                                        2.86314 * eV, 2.87053 * eV, 2.87796 * eV, 2.88543 * eV, 2.89294 * eV,
                                        2.90048 * eV, 2.90807 * eV, 2.91569 * eV, 2.92336 * eV, 2.93106 * eV,
                                        2.93881 * eV, 2.9466 * eV, 2.95443 * eV, 2.9623 * eV, 2.97021 * eV,
                                        2.97817 * eV, 2.98616 * eV, 2.9942 * eV, 3.00229 * eV, 3.01042 * eV,
                                        3.01859 * eV, 3.02681 * eV, 3.03507 * eV, 3.04337 * eV, 3.05173 * eV,
                                        3.06013 * eV, 3.06857 * eV, 3.07706 * eV, 3.0856 * eV, 3.09419 * eV,
                                        3.10282 * eV, 3.1115 * eV, 3.12023 * eV, 3.12901 * eV, 3.13784 * eV,
                                        3.14672 * eV, 3.15565 * eV, 3.16464 * eV, 3.17367 * eV, 3.18275 * eV};


    // Filled with air
    std::vector<G4double> rindex_air(energy_eV.size(), 1.0);
    airNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", false);
    G4MaterialPropertiesTable* MPT_Air = new G4MaterialPropertiesTable();
    MPT_Air->AddProperty("RINDEX", energy_eV, rindex_air);
    airNist->SetMaterialPropertiesTable(MPT_Air);

    //Detector Material:
    nist = G4NistManager::Instance();
    G4bool isotopes = false;
    Si = nist->FindOrBuildElement("Si", isotopes);
    C = nist->FindOrBuildElement("C", isotopes);
    G4double fDensity_SiC=3.22*g/cm3;
    SiC=new G4Material("SiC", fDensity_SiC,2);
    SiC->AddElement(Si,1);
    SiC->AddElement(C,1);
    DetectorMaterial=SiC;


    //Pinhole Material: 
    PinholeMaterial = DetectorMaterial; 
    std::vector<G4double> rindex_pinhole(energy_eV.size(), 0.);
    G4MaterialPropertiesTable* MPT_Pinhole = new G4MaterialPropertiesTable();
    MPT_Pinhole->AddProperty("RINDEX", energy_eV, rindex_pinhole);


    //Phantom Material 
    fPhantomMaterial = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");//(EJ200
    std::vector<G4double> rindex_phantom(energy_eV.size(), 1.58);
    std::vector<G4double> absorption_phantom(energy_eV.size(), 380*cm);
    std::vector<G4double> scint_spectrum = {
    1.29325823e-03, 1.37725624e-03, 1.49392523e-03, 1.62872116e-03,
    1.76709994e-03, 1.89678442e-03, 2.02925932e-03, 2.18866545e-03,
    2.37045550e-03, 2.54485733e-03, 2.70392592e-03, 2.86984544e-03,
    3.05539900e-03, 3.25037808e-03, 3.44686216e-03, 3.65245055e-03,
    3.87781381e-03, 4.13403911e-03, 4.42723566e-03, 4.72360702e-03,
    4.97654538e-03, 5.21865096e-03, 5.53591475e-03, 5.94006294e-03,
    6.36366509e-03, 6.78060178e-03, 7.23600177e-03, 7.75049234e-03,
    8.27073322e-03, 8.74432801e-03, 9.16317289e-03, 9.52955226e-03,
    9.84583558e-03, 1.01158770e-02, 1.03576573e-02, 1.05959763e-02,
    1.08362121e-02, 1.10680638e-02, 1.13086734e-02, 1.16099146e-02,
    1.19864184e-02, 1.23784360e-02, 1.27445080e-02, 1.31200345e-02,
    1.35390363e-02, 1.39661756e-02, 1.43594450e-02, 1.47721692e-02,
    1.52923095e-02, 1.59151659e-02, 1.65747371e-02, 1.72259721e-02,
    1.78464899e-02, 1.84457880e-02, 1.90794663e-02, 1.97759354e-02,
    2.04575732e-02, 2.10574741e-02, 2.16260419e-02, 2.22218972e-02,
    2.28063408e-02, 2.33154118e-02, 2.37603399e-02, 2.41874676e-02,
    2.45679923e-02, 2.48154196e-02, 2.49202857e-02, 2.49653928e-02,
    2.49451327e-02, 2.46872565e-02, 2.40419991e-02, 2.29673901e-02,
    2.14553942e-02, 1.96407502e-02, 1.77243405e-02, 1.58008403e-02,
    1.39000570e-02, 1.20654388e-02, 1.03490333e-02, 8.79389959e-03,
    7.43346122e-03, 6.23781585e-03, 5.11630945e-03, 4.01780240e-03,
    2.98042520e-03, 2.05731965e-03, 1.31345370e-03, 8.00507676e-04,
    4.88459520e-04, 3.20924352e-04, 2.42354569e-04, 1.99351353e-04,
    1.63441527e-04, 1.25298758e-04, 8.73247099e-05, 6.44140957e-05,
    6.23114217e-05, 7.03130337e-05, 7.61461225e-05, 6.75378790e-05 };

    G4MaterialPropertiesTable* MPT_Phantom = new G4MaterialPropertiesTable();
    MPT_Phantom->AddProperty("RINDEX", energy_eV, rindex_phantom);
    MPT_Phantom->AddProperty("ABSLENGTH", energy_eV, absorption_phantom);
    MPT_Phantom->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV);
    MPT_Phantom-> AddProperty("SCINTILLATIONCOMPONENT1", energy_eV, scint_spectrum);
    MPT_Phantom->AddConstProperty("RESOLUTIONSCALE", 1.0);
    MPT_Phantom->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1*ns);
    MPT_Phantom->AddConstProperty("SCINTILLATIONRISETIME1", 0.9*ns);   
    
    std::vector<G4double> rayleigh_length_phantom(energy_eV.size(), 1000*mm);
    MPT_Phantom->AddProperty(
        "RAYLEIGH",
        energy_eV.data(),
        rayleigh_length_phantom.data(),
        energy_eV.size()
    ); 

    fPhantomMaterial->SetMaterialPropertiesTable(MPT_Phantom);

    G4PhysicsOrderedFreeVector* spectrum = MPT_Phantom->GetProperty("SCINTILLATIONCOMPONENT1");
    std::cout << "Scint. spectrum: "<< std::endl;
    for (size_t i = 0; i < spectrum->GetVectorLength(); ++i) {
        G4double energy = spectrum->Energy(i);
        G4double value = (*spectrum)[i];
        G4cout << energy / eV << " eV, " << value << G4endl;
        }

    G4cout << "----- Material properties table printed by DetectorConstruction: -----" << G4endl;
    G4cout << "Phantom G4MaterialPropertiesTable:" << G4endl;
    MPT_Phantom->DumpTable();
        G4cout << "Pinhole G4MaterialPropertiesTable:" << G4endl;
    MPT_Pinhole->DumpTable();
        G4cout << "Air G4MaterialPropertiesTable:" << G4endl;
    MPT_Air->DumpTable();
   }



void FlashDetectorConstruction::DefineSurfaces(){

    //Surface: Phantom-world
    std::vector<G4double> reflectivity_phantom = {1., 1.};
    std::vector<G4double> energy     = {2.48 * eV,  3.18275 * eV}; //lamda in range 400-500 nm; lambda [nm] = 1240/E[eV]

    PhantomOpticalSurface = new G4OpticalSurface("PhantomOpticalSurface");
    PhantomOpticalSurface->SetType(dielectric_dielectric);  
    PhantomOpticalSurface->SetModel(unified); 
    PhantomOpticalSurface->SetFinish(polished);  
    G4MaterialPropertiesTable* WrappingProperty = new G4MaterialPropertiesTable();
    WrappingProperty->AddProperty("REFLECTIVITY", energy, reflectivity_phantom); 
    PhantomOpticalSurface->SetMaterialPropertiesTable(WrappingProperty);
    PhantomSurface = new G4LogicalBorderSurface("PhantomOpticalSurface", fPhantom_physical, physicalTreatmentRoom, PhantomOpticalSurface);

}




G4VPhysicalVolume *FlashDetectorConstruction::ConstructPhantom(G4double CollPos) {
    /*
   This function creates a cubic phantom with material properties defined in the `DefineMaterials` method. 
   The phantom's dimensions, position, and other properties are set based on member variables. 
   The phantom is placed in the treatment room volume and is associated with a dedicated region 
   to control simulation settings like step limits and physics cuts.

   Parameters:
   - CollPos: A coefficient that determines the phantom's position along the X-axis, relative to the applicator.

   Key Steps:
   1. Calculate the phantom's position (`fPhantomPosition`) based on `CollPos`, 
      including a shift to account for the air gap.
   2. Construct a G4Box object for the phantom geometry.
   3. Create a logical volume (`fPhantomLogicalVolume`) using the predefined phantom material.
   4. Place the logical volume within the treatment room as a physical volume.
   5. Define a dedicated region (`PhantomRegion`) for the phantom to apply specific physics cuts.
   6. Set a step limit for the phantom to control the maximum step size during tracking.
   7. Apply visualization attributes to the phantom for graphical rendering.

   Returns:
   - A pointer to the physical volume (`fPhant_phys`) representing the phantom in the simulation.
    */

    fPosition_coefficient = CollPos;
    fPhantom_coordinateX = (fPosition_coefficient * mm + fPhantomSizeX / 2);
    fPhantomPosition =  G4ThreeVector(fPhantom_coordinateX, 0. * mm, 0. * mm); //phantom is constructed with the entrance surface attached to the applicator 
    
    fPhantom = new G4Box("Phantom", fPhantomSizeX / 2, fPhantomSizeY / 2, fPhantomSizeZ / 2);
    fPhantomLogicalVolume = new G4LogicalVolume(fPhantom, fPhantomMaterial, "phantomLog", 0, 0, 0);
    fPhant_phys = new G4PVPlacement(0, fPhantomPosition, "phantomPhys", fPhantomLogicalVolume, physicalTreatmentRoom, false, 0);
    
    //define the region to set cuts in FlashPhysicsList.cc and step limit
    G4Region *PhantomRegion = new G4Region("Phantom_reg");
    fPhantomLogicalVolume->SetRegion(PhantomRegion);
    PhantomRegion->AddRootLogicalVolume(fPhantomLogicalVolume);

    //set step limit in phantom
    G4double maxStep = 0.1 * mm;
    fStepLimit = new G4UserLimits(maxStep);
    fPhantomLogicalVolume->SetUserLimits(fStepLimit);

    // Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
    red->SetVisibility(true);
    blue = new G4VisAttributes(G4Colour(0 / 255., 0. / 255., 255. / 255.));
    blue->SetVisibility(true);
    fPhantomLogicalVolume->SetVisAttributes(red);

    return fPhant_phys;

}

std::vector<G4VPhysicalVolume*> FlashDetectorConstruction::ConstructPinhole(G4double CollRadius) {
    /* 
    This function constructs a pinhole collimator composed of multiple sheets, each featuring a central hole. 
   The sheets are modeled as truncated pyramids (`G4Trd`) with a cylindrical hole (`G4Tubs`) subtracted 
   to create the pinhole structure. The sheets are positioned and rotated to form a cubical arrangement around the phantom.

   Parameters:
   - collimatorThickness: Thickness of each collimator sheet.
   - squareSize: Size of the square cross-section of the sheet, influenced by `PinholeDistance`.
   - innerRadius: Radius of the central hole.

   Key Steps:
   1. Define the geometry of a single pinhole sheet by subtracting a cylinder from a truncated pyramid.
   2. Create a logical volume (`PinholeLogicalVolume`) using the defined geometry and material properties.
   3. Calculate the positions and orientations for each sheet.
   4. Place the sheets in the treatment room, forming a cubical arrangement around the phantom, using `G4PVPlacement` 
      and appropriate rotation matrices.
   5. Assign visualization attributes for graphical representation.
   6. Return a vector of physical volumes representing the constructed pinhole collimator sheets.

   Returns:
   - A `std::vector<G4VPhysicalVolume*>` containing the physical volumes of all pinhole sheets.

   Notes:
   - Ensures proper alignment and non-overlapping placement of sheets.
   - Uses rotation matrices to orient the sheets correctly in the 3D space.
    */
    G4double pinholeThickness = 0.11*mm;
    G4double pinholeSquareSize = fPhantomSizeX + PinholeDistance * 2;
    G4double innerRadius = 0.250 * mm;

    // Geometry
    G4Trd* squareSolid = new G4Trd("BlackSheet", (pinholeSquareSize + pinholeThickness) / 2, 
                                   (pinholeSquareSize - pinholeThickness) / 2, 
                                   (pinholeSquareSize + pinholeThickness) / 2, 
                                   (pinholeSquareSize - pinholeThickness) / 2, 
                                   pinholeThickness / 2);
    G4Tubs* innerCylinder = new G4Tubs("InnerCylinder", 0, innerRadius, pinholeThickness / 2, 0. * deg, 360. * deg);
    G4SubtractionSolid* PinholeCilinder = new G4SubtractionSolid("Pinhole", squareSolid, innerCylinder);
    PinholeLogicalVolume = new G4LogicalVolume(PinholeCilinder, PinholeMaterial, "pinholeLog", 0, 0, 0);


    G4Tubs* innerCylinder_back = new G4Tubs("InnerCylinder", 0, CollRadius , pinholeThickness / 2, 0. * deg, 360. * deg);
    G4SubtractionSolid* PinholeCilinder_back = new G4SubtractionSolid("PinholeBack", squareSolid, innerCylinder_back);
    PinholeLogicalVolume_back = new G4LogicalVolume(PinholeCilinder_back, PinholeMaterial, "pinholeLogBack", 0, 0, 0);

    G4double PinholePosition_l = fPhantomSizeX * 0.5 + PinholeDistance;
    G4double PinholePosition_t = fPhantomSizeX + PinholeDistance;

    // Pinhole 1
    G4RotationMatrix* rotationMatrix_y = new G4RotationMatrix();
    rotationMatrix_y->rotateY(90. * deg); // Rotazione attorno all'asse Y
    G4VPhysicalVolume* Pinhole_phys1 = new G4PVPlacement(rotationMatrix_y, 
                                                         G4ThreeVector(PinholePosition_t + fAirGap, 0., 0.), 
                                                         "pinholePhys", 
                                                         PinholeLogicalVolume, 
                                                         physicalTreatmentRoom, false, 0, fCheckOverlaps);
    std::cout<< "Pinhole on the front: "<< std::endl; 
    std::cout<< G4ThreeVector(PinholePosition_t + fAirGap, 0., 0.)<< std::endl<< std::endl; 

    // Pinhole 2
    G4RotationMatrix* rotationMatrix_x1 = new G4RotationMatrix();
    rotationMatrix_x1->rotateX(-90. * deg); // Rotazione attorno all'asse X
    G4VPhysicalVolume* Pinhole_phys2 = new G4PVPlacement(rotationMatrix_x1, 
                                                         G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, PinholePosition_l, 0.), 
                                                         "pinholePhys", 
                                                         PinholeLogicalVolume, 
                                                         physicalTreatmentRoom, false, 0, fCheckOverlaps);
    std::cout<< "Pinhole on the side: "<< std::endl; 
    std::cout<< G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, PinholePosition_l, 0.)<< std::endl<< std::endl; 


    // Pinhole 3
    G4RotationMatrix* rotationMatrix_x2 = new G4RotationMatrix();
    rotationMatrix_x2->rotateX(90. * deg); // Rotazione attorno all'asse X
    G4VPhysicalVolume* Pinhole_phys3 = new G4PVPlacement(rotationMatrix_x2, 
                                                         G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, -PinholePosition_l, 0.), 
                                                         "pinholePhys", 
                                                         PinholeLogicalVolume, 
                                                         physicalTreatmentRoom, false, 0, fCheckOverlaps);

    // Pinhole 4
    G4RotationMatrix* reverse = new G4RotationMatrix();
    reverse->rotateX(180. * deg); // Rotazione inversa attorno all'asse X
    G4VPhysicalVolume* Pinhole_phys4 = new G4PVPlacement(reverse, 
                                                         G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, 0., PinholePosition_l), 
                                                         "pinholePhys", 
                                                         PinholeLogicalVolume, 
                                                         physicalTreatmentRoom, false, 0, fCheckOverlaps);

    // Pinhole 5
    G4VPhysicalVolume* Pinhole_phys5 = new G4PVPlacement(0, 
                                                         G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, 0., -PinholePosition_l), 
                                                         "pinholePhys", 
                                                         PinholeLogicalVolume, 
                                                         physicalTreatmentRoom, false, 0, fCheckOverlaps);


    // Pinhole 6
    G4RotationMatrix* rotationMatrix_y2 = new G4RotationMatrix();
    rotationMatrix_y2->rotateY(-90. * deg); // Rotazione attorno all'asse Y
    G4VPhysicalVolume* Pinhole_phys6 = new G4PVPlacement(rotationMatrix_y2, 
                                                         G4ThreeVector(-PinholeDistance + fAirGap, 0., 0.), 
                                                         "pinholePhys", 
                                                         PinholeLogicalVolume_back, 
                                                         physicalTreatmentRoom, false, 0, fCheckOverlaps);

    // Visualization
    gray = new G4VisAttributes(G4VisAttributes::GetInvisible());
    gray->SetVisibility(false);
    PinholeLogicalVolume->SetVisAttributes(G4VisAttributes::GetInvisible());
    PinholeLogicalVolume_back->SetVisAttributes(G4VisAttributes::GetInvisible());


    std::vector<G4VPhysicalVolume*> fPinhole_Phys = {Pinhole_phys1, Pinhole_phys2, Pinhole_phys3, Pinhole_phys4, Pinhole_phys5,Pinhole_phys6};
    return fPinhole_Phys;
}



std::vector<G4VPhysicalVolume*> FlashDetectorConstruction::ConstructDetector(){
    /*
    This function constructs a detector with a box shape (`G4Box`) and places it in five different positions around 
   the phantom to form a cubical arrangement. The detector is constructed using the material properties defined 
   in the `DefineMaterials()` function. The detector is rotated and positioned using `G4RotationMatrix` to ensure 
   proper orientation in the treatment room.

   Parameters:
   - Det_thickness: Thickness of the detector.
   - Det_width: Width of the detector.
   - DetectorMaterial: Material used to construct the detector.
   - DetectorDistance: Distance from the phantom to the detector.
   - AirGap: Gap between the phantom and the detector.
   - fPhantomSizeX: Size of the phantom in the X direction.

   Key Steps:
   1. Define the geometry of the detector using `G4Box` with specified thickness and width.
   2. Calculate the positions of the detector, taking into account the phantom size and specified distances.
   3. Create rotation matrices to rotate the detector into the correct orientation.
   4. Create logical volumes for the detector and place them in the treatment room using `G4PVPlacement`.
   5. Assign visualization attributes to the detector for graphical representation.

   Returns:
   - A `std::vector<G4VPhysicalVolume*>` containing the physical volumes of the detector placed in the treatment room.

   Notes:
   - The detector is placed in five different positions: along the X, Y, and Z axes with various rotations to form 
     a cubical arrangement around the phantom.
   - The air gap and the distance between the phantom and the detector are considered when calculating the positions.
   */
    
    G4double detectorSquareSize = fPhantomSizeX + DetectorDistance * 2;
    G4double detectorThickness = 0.1 * cm;

    // Definizione della geometria del detector
    G4Trd* squareSolid = new G4Trd("Detector", 
                                   (detectorSquareSize + detectorThickness) / 2, 
                                   (detectorSquareSize - detectorThickness) / 2, 
                                   (detectorSquareSize + detectorThickness) / 2, 
                                   (detectorSquareSize - detectorThickness) / 2, 
                                   detectorThickness / 2);

    fDetectorPosition_l = fPhantomSizeX * 0.5 + DetectorDistance;
    fDetectorPosition_t = fPhantomSizeX + DetectorDistance;

    fDetLogicalVolume = new G4LogicalVolume(squareSolid, DetectorMaterial, "DetectorLog", 0, 0, 0);

    // Detector 1
    G4RotationMatrix* rotationMatrix_y = new G4RotationMatrix();
    rotationMatrix_y->rotateY(90. * deg); // Rotazione attorno all'asse Y
    fDet_phys1 = new G4PVPlacement(rotationMatrix_y, G4ThreeVector(fDetectorPosition_t + fAirGap, 0., 0.), 
                                                       "DetPhys", 
                                                       fDetLogicalVolume, 
                                                       physicalTreatmentRoom, false, 0, fCheckOverlaps);

    std::cout<< "Detector on the front: "<< std::endl; 
    std::cout<< G4ThreeVector(fDetectorPosition_t + fAirGap, 0., 0.) << std::endl<< std::endl; 

    // Detector 2
    G4RotationMatrix* rotationMatrix_x1 = new G4RotationMatrix();
    rotationMatrix_x1->rotateX(-90. * deg);
    fDet_phys2 = new G4PVPlacement(rotationMatrix_x1, G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, fDetectorPosition_l, 0.), 
                                                       "DetPhys", 
                                                       fDetLogicalVolume, 
                                                       physicalTreatmentRoom, false, 0, fCheckOverlaps);

    std::cout<< "Detector on the side: "<< std::endl; 
    std::cout<< G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, fDetectorPosition_l, 0.) << std::endl<< std::endl; 


    // Detector 3
    G4RotationMatrix* rotationMatrix_x2 = new G4RotationMatrix();
    rotationMatrix_x2->rotateX(90. * deg); // Rotazione attorno all'asse X
    fDet_phys3 = new G4PVPlacement(rotationMatrix_x2, G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, -fDetectorPosition_l, 0.), 
                                                       "DetPhys", 
                                                       fDetLogicalVolume, 
                                                       physicalTreatmentRoom, false, 0, fCheckOverlaps);

    // Detector 4
    G4RotationMatrix* reverse = new G4RotationMatrix();
    reverse->rotateX(180. * deg); // Rotazione inversa attorno all'asse X
    fDet_phys4 = new G4PVPlacement(reverse, G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, 0., fDetectorPosition_l), 
                                                       "DetPhys", 
                                                       fDetLogicalVolume, 
                                                       physicalTreatmentRoom, false, 0, fCheckOverlaps);

    // Detector 5
    fDet_phys5 = new G4PVPlacement(0, G4ThreeVector(fPhantomSizeX * 0.5 + fAirGap, 0., -fDetectorPosition_l), 
                                                       "DetPhys", 
                                                       fDetLogicalVolume, 
                                                       physicalTreatmentRoom, false, 0, fCheckOverlaps);

    // Attributi di visualizzazione
    gray = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    gray->SetVisibility(true);
    fDetLogicalVolume->SetVisAttributes(gray);

    // Ritorno dei volumi fisici in ordine
    std::vector<G4VPhysicalVolume*> fDet_Phys = {fDet_phys1, fDet_phys2, fDet_phys3, fDet_phys4, fDet_phys5};
    return fDet_Phys;
}



G4VPhysicalVolume *FlashDetectorConstruction::Construct() {
    // -----------------------------
    // Treatment room - World volume
    //------------------------------
    const G4double worldX = 400.0 * cm;
    const G4double worldY = 400.0 * cm;
    const G4double worldZ = 400.0 * cm;

    G4Box *treatmentRoom = new G4Box("TreatmentRoom", worldX, worldY, worldZ);
    logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, airNist, "logicTreatmentRoom", 0, 0, 0);
    physicalTreatmentRoom = new G4PVPlacement(0, G4ThreeVector(), "physicalTreatmentRoom", logicTreatmentRoom, 0, false, 0);

    // The treatment room is invisible in the Visualisation
    logicTreatmentRoom->SetVisAttributes(G4VisAttributes::GetInvisible());

    // -----------------------------
    // Applicator + phantom + Default dimensions
    //------------------------------
    Collimator = new Applicator(physicalTreatmentRoom);
    fPhantom_physical = ConstructPhantom(Collimator->fFinalApplicatorXPositionFlash +Collimator->fHightFinalApplicatorFlash+fAirGap);


    // -----------------------------
    // Pinhole camera
    //------------------------------  
    Pihole_physical = ConstructPinhole( Collimator->fOuterRadiusFirstApplicatorFlash);

    // -----------------------------
    // Detector pannel
    //------------------------------
    Detector_physical = ConstructDetector();
    
    DefineSurfaces(); 
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


