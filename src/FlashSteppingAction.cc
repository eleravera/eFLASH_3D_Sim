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
#include "G4OpBoundaryProcess.hh"
#include <common.hh>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "FlashRunAction.hh"
#include "G4AnalysisManager.hh"
#include <unordered_set>

G4int FlashSteppingAction::TransmissionCount = 0 ;
G4int FlashSteppingAction::FresnelRefractionCount = 0 ; 
G4int FlashSteppingAction::FresnelReflectionCount = 0 ;
G4int FlashSteppingAction::TotalInternalReflectionCount = 0 ;
G4int FlashSteppingAction::LambertianReflectionCount = 0 ; 
G4int FlashSteppingAction::LobeReflectionCount = 0 ;
G4int FlashSteppingAction::SpikeReflectionCount = 0 ; 
G4int FlashSteppingAction::BackScatteringCount = 0 ; 
G4int FlashSteppingAction::AbsorptionCount = 0 ; 
G4int FlashSteppingAction::PhotonTotalInternalReflectionCount = 0 ;
G4int FlashSteppingAction::PhotonRefractionCount = 0 ;
G4int FlashSteppingAction::PhotonReflectionCount = 0 ;
G4int FlashSteppingAction::PhotonsOutOfWorld = 0 ;



FlashSteppingAction::FlashSteppingAction(FlashEventAction *)
    : G4UserSteppingAction() {}

FlashSteppingAction::~FlashSteppingAction() {}

void FlashSteppingAction::UserSteppingAction(const G4Step *aStep)
{
    //G4AnalysisManager* analysisMan = G4AnalysisManager::Instance(); //DA IMPLEMENTARE PER VEDERE ANGOLI E RISOLUZIONE. 

    G4Track* track = aStep->GetTrack();
    G4StepPoint *postStep = aStep->GetPostStepPoint();
    G4StepPoint *preStep = aStep->GetPreStepPoint();
    static G4ParticleDefinition* opticalphoton = G4OpticalPhoton::OpticalPhotonDefinition();

    G4OpBoundaryProcessStatus theStatus = Undefined;
    static G4ThreadLocal G4OpBoundaryProcess* boundary = NULL; 
    //find the boundary process only once
    if(!boundary){
      G4ProcessManager* pm = aStep->GetTrack()->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
      G4ProcessVector* pv = pm->GetProcessList();
      G4int i;
      for( i=0;i<nprocesses;i++){
        if((*pv)[i]->GetProcessName()=="OpBoundary"){
          boundary = (G4OpBoundaryProcess*)(*pv)[i];
          break;
          }
        }
      }


    if(track->GetDefinition() == opticalphoton) { 
      if(!(aStep->GetPostStepPoint()->GetPhysicalVolume())){//out of world
        PhotonsOutOfWorld++; 
        return;}
    
      if(postStep->GetStepStatus() == fGeomBoundary) { // if at boundary
        theStatus = boundary->GetStatus();

        switch(theStatus){
            case Absorption: 
              AbsorptionCount++;
              break;
            case FresnelRefraction:
              FresnelRefractionCount++;
              PhotonRefractionCount++;
              break;
            case FresnelReflection:
              FresnelReflectionCount++;
              PhotonReflectionCount++;
              break;
            case TotalInternalReflection: 
              TotalInternalReflectionCount++;
              PhotonTotalInternalReflectionCount++; 
              if (PhotonTotalInternalReflectionCount > 10) { // Kill photons with internal reflection very high
                track->SetTrackStatus(fStopAndKill);
              } 
              break;
            case LambertianReflection:
              LambertianReflectionCount++;
              break;
            case LobeReflection:
              LobeReflectionCount++;
              break;
            case SpikeReflection:
              SpikeReflectionCount++;
              break;
            case BackScattering:
              BackScatteringCount++;
              break;
            default: 
              break;
            }
        }

    G4String volumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    G4String prevolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

    /*if (prevolumeName == "logicTreatmentRoom" && volumeName ==  "pinholeLog"){

      G4double pos_x = aStep->GetTrack()->GetPosition().x();
      G4double pos_y = aStep->GetTrack()->GetPosition().y();
      G4double pos_z = aStep->GetTrack()->GetPosition().z();
      
      detection photon_maps =  detection(pos_x/mm, pos_y/mm, pos_z/mm);
      detection_vector.push_back(photon_maps);
      std::cout << "PINHOLE: " << std::endl; 
      photon_maps.print();
      std::cout << std::endl << std::endl; 

    }*/

    if (prevolumeName == "logicTreatmentRoom" && volumeName ==  "DetectorLog"){

      G4double pos_x = aStep->GetTrack()->GetPosition().x();
      G4double pos_y = aStep->GetTrack()->GetPosition().y();
      G4double pos_z = aStep->GetTrack()->GetPosition().z();
      
      detection photon_maps =  detection(pos_x/mm, pos_y/mm, pos_z/mm);
      detection_vector.push_back(photon_maps);
      /*std::cout << "DETECTOR: " << std::endl; 
      photon_maps.print();
      std::cout << std::endl << std::endl; */
    }
    
    }   /* end of if optical photon */


    


}







/*//Save photons 
          G4String volumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
          G4String prevolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

          if (prevolumeName == "phantomLog" && volumeName ==  "logicTreatmentRoom") {
            G4ThreeVector photonDirection = track->GetMomentumDirection();  
            // Define a small tolerance value as cosThetaMax
            const G4double cosThetaMax = std::cos(0.5 * CLHEP::pi / 180.0);  // conv degree in radians

            // Calcola il coseno dell'angolo rispetto agli assi
            G4double cosThetaX = photonDirection.x();  // Prodotto scalare con (1,0,0)
            G4double cosThetaY = photonDirection.y();  // Prodotto scalare con (0,1,0)
            G4double cosThetaZ = photonDirection.z();  // Prodotto scalare con (0,0,1)

            // Se il fotone è entro 1 grado rispetto a uno degli assi principali killalo
            if (!(cosThetaX > cosThetaMax || cosThetaY > cosThetaMax || cosThetaZ > cosThetaMax)) {
                track->SetTrackStatus(fStopAndKill); 

            } 
            }

      if (prevolumeName == "logicTreatmentRoom" && volumeName ==  "DetectorLog") {
        G4double pos_x = aStep->GetTrack()->GetPosition().x();
        G4double pos_y = aStep->GetTrack()->GetPosition().y();
        G4double pos_z = aStep->GetTrack()->GetPosition().z();

        // append to detection_vector the current info
        detection photon_maps =  detection(pos_x/mm, pos_y/mm, pos_z/mm);
        detection_vector1.push_back(photon_maps);
        //std::cout << "A photon has been saved on file" << std::endl;  
        //photon_maps.print();
        }    */