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
/// \file FlashDetectorConstruction.hh

#ifndef FlashDetectorConstruction_h
#define FlashDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "tls.hh"
#include "G4ThreeVector.hh"
#include "G4UserLimits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class Applicator;
class G4VSensitiveDetector;
class G4NistManager;
class G4Tubs;
class G4Box;
class G4Element;
class G4VisAttributes;


class FlashDetectorMessenger;


class FlashDetectorConstruction : public G4VUserDetectorConstruction {
public:
  //Constructor and destructor for FlashDetectorConstruction class
  FlashDetectorConstruction();
  virtual ~FlashDetectorConstruction();

  //Logical and Physical volume of the treatment room
  G4LogicalVolume *logicTreatmentRoom;
  G4VPhysicalVolume *physicalTreatmentRoom;

  //The functions to create Phantom, Pinhole and Photodetector - and so set the carachteristic of the detector
  G4VPhysicalVolume *ConstructPhantom(G4double CollPos);
  std::vector<G4VPhysicalVolume*> ConstructPinhole(G4double CollRadius);
  std::vector<G4VPhysicalVolume*> ConstructDetector();
  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField(); // ? 

  G4bool  SetPhantomMaterial(G4String material);
  G4bool  SetDetectorMaterial(G4String material);
  void SetAirGap(G4double position);
  void SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetPinholeDistance(G4double distanceX);
  void SetDetectorDistance(G4double distanceX);
  void DefineMaterials();
  void DefineSurfaces(); 

  //Visualization attributes
  G4VisAttributes *skyBlue;
  G4VisAttributes *red;
  G4VisAttributes *blue;
  G4VisAttributes *green;
  G4VisAttributes *gray;

  
private:

  FlashDetectorMessenger* fDetectorMessenger;

  G4Material *airNist;
  G4Material *fPhantomMaterial;
  G4Material *PinholeMaterial;
  Applicator *Collimator;

  G4OpticalSurface *PhantomOpticalSurface;
  G4LogicalBorderSurface *PhantomSurface; 

  //Attributes set by some public functions, such as SetAirGap, SetPhantomSize.. 
  G4double fAirGap;
  G4double fPhantomSizeX, fPhantomSizeY, fPhantomSizeZ; 
  G4double fPhantom_coordinateX, fPosition_coefficient;
  G4ThreeVector fPhantomPosition;
  G4double PinholeDistance;
  G4double DetectorDistance;

  G4double fDet_thickness,fDet_width, fDet_sub_thickness; 
  G4double fDetectorPosition_t, fDetectorPosition_l, fAirGap_phantom_det;
  //Materials 
  G4Element *Si;
  G4Element *C;
  G4Material *SiC;
  G4Material *DetectorMaterial;

  G4Box *Det_box;
  G4LogicalVolume *fDetLogicalVolume;
  G4VPhysicalVolume *fDet_phys1, *fDet_phys2, *fDet_phys3, *fDet_phys4, *fDet_phys5;
  std::vector<G4VPhysicalVolume*> Detector_physical; 

  //Phantom
  G4Box *fPhantom;
  G4LogicalVolume *fPhantomLogicalVolume;
  G4VPhysicalVolume *fPhant_phys;
  G4VPhysicalVolume *fPhantom_physical;


  //Pinhole 
  G4LogicalVolume *PinholeLogicalVolume;
  G4LogicalVolume *PinholeLogicalVolume_back;
  G4VPhysicalVolume *Pihole_phys1, *Pihole_phys2, *Pihole_phys3, *Pihole_phys4, *Pihole_phys5;
  std::vector<G4VPhysicalVolume*> Pihole_physical;


  G4UserLimits *fStepLimit;
  G4bool fCheckOverlaps;
 
  G4NistManager *nist;

};

#endif
