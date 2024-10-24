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
#include "Applicator.hh"

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
  G4VPhysicalVolume *physicalTreatmentRoom;
  G4LogicalVolume *logicTreatmentRoom;
  G4VPhysicalVolume *ConstructPhantom(G4double CollPos, G4VPhysicalVolume *fPhysiFirstApplicatorFlash);
  G4VPhysicalVolume *ConstructPinhole();

  G4VPhysicalVolume *ConstructDetector();

  FlashDetectorConstruction();
  virtual ~FlashDetectorConstruction();

  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();

  G4bool  SetPhantomMaterial(G4String material);
  G4bool  SetDetectorMaterial(G4String material);
  void SetAirGap(G4double position);
  void SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetPinholeDistance(G4double distanceX);
  void SetDetectorDistance(G4double distanceX);
  void SetDetectorThickness(G4double thickness);
  void SetDetectorWidth(G4double width);

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

  G4double fAirGap;
  G4double fPhantomSizeX, fPhantomSizeY, fPhantomSizeZ, fPhantom_coordinateX,fPosition_coefficient;
  G4ThreeVector fPhantomPosition;
  G4double PinholeDistance;
  G4double DetectorDistance;

  G4double fDet_thickness,fDet_width,fDet_sub_thickness,fDetectorPosition_t,fDetectorPosition_l,fAirGap_phantom_det;
  G4Element *Si;
  G4Element *C;
  G4Material *SiC;
  G4Material *DetectorMaterial;

  G4Box *fPhantom;
  //G4SubtractionSolid *PinholeCilinder;

  G4Box *Det_box;
  G4LogicalVolume *fDetLogicalVolume;
  G4VPhysicalVolume *fDet_phys1, *fDet_phys2, *fDet_phys3;//, *fDet_phys4, *fDet_phys5;

  //G4Box *fDet_sub;
  //G4LogicalVolume *fDet_sub_LogicalVolume;
  //G4VPhysicalVolume *fDet_sub_phys;

  
  void DefineMaterials();


  G4LogicalVolume *fPhantomLogicalVolume;
  G4VPhysicalVolume *fPhant_phys;
  G4VPhysicalVolume *fPhantom_physical;
    
  G4LogicalVolume *PinholeLogicalVolume;
  G4VPhysicalVolume *Pihole_phys;
  G4VPhysicalVolume *Pihole_physical;


  G4UserLimits *fStepLimit;
  G4bool fCheckOverlaps;
 
  G4NistManager *nist;

};

#endif
