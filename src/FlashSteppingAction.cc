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

G4int FlashSteppingAction::TransmissionCount = 0 ; 
G4int FlashSteppingAction::FresnelRefractionCount = 0 ; 
G4int FlashSteppingAction::FresnelReflectionCount = 0 ;
G4int FlashSteppingAction::TotalInternalReflectionCount = 0 ;
G4int FlashSteppingAction::LambertianReflectionCount = 0 ; 
G4int FlashSteppingAction::LobeReflectionCount = 0 ;
G4int FlashSteppingAction::SpikeReflectionCount = 0 ; 
G4int FlashSteppingAction::BackScatteringCount = 0 ; 
G4int FlashSteppingAction::AbsorptionCount = 0 ; 

int i = 0 ;

FlashSteppingAction::FlashSteppingAction(FlashEventAction *)
    : G4UserSteppingAction() {}

FlashSteppingAction::~FlashSteppingAction() {}

void FlashSteppingAction::UserSteppingAction(const G4Step *aStep)
{
    G4AnalysisManager* analysisMan = G4AnalysisManager::Instance();


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


    /*if (track->GetDefinition() == opticalphoton && postStep->GetStepStatus() == fGeomBoundary) {

        G4String volumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
        G4String prevolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
      }*/



    if(track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) { 
    
      if(!(aStep->GetPostStepPoint()->GetPhysicalVolume())){//out of world
        return;}
    
      G4String thePrePV  = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
      G4String thePostPV = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
   
    if(postStep->GetStepStatus() == fGeomBoundary){ // if at boundary
        
      theStatus = boundary->GetStatus();

      std::cout << "\n pre-volume and post-volume: " << thePrePV << "   " << thePostPV << std::endl; 
      
      /* Count detected photons */
      switch(theStatus){
      case Absorption: // absorption at every boundary
        AbsorptionCount++;
        std::cout << "AbsorptionCount ++: " << AbsorptionCount << std::endl; 
        break;

      case FresnelRefraction:
        FresnelRefractionCount++;
        std::cout << "FresnelRefractionCount ++: " << FresnelRefractionCount << std::endl; 

            break;
      
      case FresnelReflection:
        FresnelReflectionCount++;
        std::cout << "FresnelReflectionCount ++: " << FresnelReflectionCount << std::endl; 
            break;
      
      case TotalInternalReflection:
        TotalInternalReflectionCount++;
        std::cout << "TotalInternalReflectionCount ++: " << TotalInternalReflectionCount << std::endl; 
        track->SetTrackStatus(fStopAndKill);
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
  
    } /* end of switch */
   }  /* end of if at boundary */
  }   /* end of if optical photon */

    

}


        /*if (prevolumeName == "logicTreatmentRoom" && volumeName ==  "DetectorLog") {
            G4double pos_x = aStep->GetTrack()->GetPosition().x();
            G4double pos_y = aStep->GetTrack()->GetPosition().y();
            G4double pos_z = aStep->GetTrack()->GetPosition().z();

            // append to detection_vector the current info
            //detection photon_maps =  detection(pos_x/mm, pos_y/mm, pos_z/mm);
            //detection_vector1.push_back(photon_maps);
            std::cout << "A photon has been saved on file: positions are : " << pos_x << " , " << pos_y << " , " << pos_z  << std::endl;  
            i++;
            std::cout<<"i: " << i<< std::endl; 
            //photon_maps.print();
            } */  