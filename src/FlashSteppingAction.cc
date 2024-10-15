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
/// \file FlashSteppingAction.cc
/// \brief Implementation of the FlashSteppingAction class

#include "FlashSteppingAction.hh"
#include "FlashDetectorConstruction.hh"
#include "FlashEventAction.hh"
#include "FlashRunAction.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include <sstream>
#include <string>
#include "G4OpticalPhoton.hh"
#include <common.hh>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

double FlashSteppingAction::killedPhotonCount = 0;
double FlashSteppingAction::killedPhotonCount_applicator = 0;
double FlashSteppingAction::killedPhotonCount_phantom = 0;


FlashSteppingAction::FlashSteppingAction(FlashEventAction *)
    : G4UserSteppingAction() {}

FlashSteppingAction::~FlashSteppingAction() {}

void FlashSteppingAction::UserSteppingAction(const G4Step *aStep)
{

  G4int eventid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  G4StepPoint *postStep = aStep->GetPostStepPoint();
  G4StepPoint *preStep = aStep->GetPreStepPoint();


  if (postStep->GetStepStatus() == fGeomBoundary)
    {
      G4String volumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
      G4String prevolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
      G4Track* track = aStep->GetTrack();
      const G4ParticleDefinition* particleDef = track->GetDefinition();
      
      
      if(particleDef == G4OpticalPhoton::OpticalPhotonDefinition()) { 
        
          //Kill photons that exit non parallel from the phantom 
          if (prevolumeName == "phantomLog" && volumeName ==  "logicTreatmentRoom") {
            G4ThreeVector photonDirection = track->GetMomentumDirection();  
            // Define a small tolerance value as cosThetaMax
            const G4double cosThetaMax = std::cos(0.1 * CLHEP::pi / 180.0);  // conv degree in radians

            // Calcola il coseno dell'angolo rispetto agli assi
            G4double cosThetaX = photonDirection.x();  // Prodotto scalare con (1,0,0)
            G4double cosThetaY = photonDirection.y();  // Prodotto scalare con (0,1,0)
            G4double cosThetaZ = photonDirection.z();  // Prodotto scalare con (0,0,1)

            // Se il fotone è entro 1 grado rispetto a uno degli assi principali killalo
            if (!(cosThetaX > cosThetaMax || cosThetaY > cosThetaMax || cosThetaZ > cosThetaMax)) {
                track->SetTrackStatus(fStopAndKill); 
                killedPhotonCount++;

            } else {
                // std::cout << "A photon passed the selection" << std::endl;
            }

            }

          if (prevolumeName == "FirstApplicatorFlash" && volumeName ==  "phantomLog") {
                track->SetTrackStatus(fStopAndKill); 
                killedPhotonCount_phantom++;
                //std::cout<<"Photons phantom Log-> FirstApplicatorFlash has been killed " << std::endl; 

          }

          if (prevolumeName == "phantomLog" && volumeName ==  "FirstApplicatorFlash") {
                track->SetTrackStatus(fStopAndKill); 
                killedPhotonCount_applicator++;
                //std::cout<<"Photons phantom Log-> FirstApplicatorFlash has been killed " << std::endl; 

          }



          //Save photons 
          if (prevolumeName == "logicTreatmentRoom" && volumeName ==  "DetectorLog") {
            G4double pos_x = aStep->GetTrack()->GetPosition().x();
            G4double pos_y = aStep->GetTrack()->GetPosition().y();
            G4double pos_z = aStep->GetTrack()->GetPosition().z();

            // append to detection_vector the current info
            detection photon_maps =  detection(pos_x/mm, pos_y/mm, pos_z/mm);
            detection_vector1.push_back(photon_maps);
            //std::cout << "A photon has been saved on file" << std::endl;  
            //photon_maps.print();
            } 
        
       }
    }


}