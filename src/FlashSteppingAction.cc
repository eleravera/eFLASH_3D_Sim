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

G4int FlashSteppingAction::TransmissionCount = 0 ; //ma sto facendo bene a fare questa cosa o è una porcheria???
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

static std::unordered_set<G4int> trappedPhoton;


FlashSteppingAction::FlashSteppingAction(FlashEventAction *)
    : G4UserSteppingAction() {}

FlashSteppingAction::~FlashSteppingAction() {}

void FlashSteppingAction::UserSteppingAction(const G4Step *aStep)
{
    //G4AnalysisManager* analysisMan = G4AnalysisManager::Instance();

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
      //reset this variable for every photon
    
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
              //std::cout << "REFRACTION: " << track->GetTrackID()  << std::endl; 
              break;
            case FresnelReflection:
              FresnelReflectionCount++;
              PhotonReflectionCount++;
              //std::cout << "REFLECTION: " << track->GetTrackID()  << std::endl; 
              break;
            case TotalInternalReflection: 
              TotalInternalReflectionCount++;
              PhotonTotalInternalReflectionCount++; 
              trappedPhoton.insert(track->GetTrackID()); // Add photon to trapped set if it undergoes total internal reflection
              //std::cout << "TOTAL INTERNAL REFLECTION: " << track->GetTrackID()  << std::endl; 
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
        
        if (track->GetTrackStatus() == fStopAndKill) { // Check if the photon is absorbed within the volume and had previously undergone TIR
            G4String thePrePV = preStep->GetPhysicalVolume()->GetName();

            photonProcess::AbsorptionLocation loc;
            G4ThreeVector genPosition = track->GetVertexPosition();
            G4ThreeVector momentumDir = track->GetVertexMomentumDirection();
            G4double theta = momentumDir.angle(G4ThreeVector(0, 1, 0)); // angle with the y-axis
            G4double phi = std::atan2(momentumDir.x(), momentumDir.z()); // angle with the zx-plane. It should be zero if on z-axis

            if (trappedPhoton.find(track->GetTrackID()) != trappedPhoton.end()) {// If the track ID exists in the set, the photon was previously "trapped"
                
                if (thePrePV == "phantomPhys"){
                  //std::cout << "Photon absorbed in volume: " << thePrePV << " with a track id = " << track->GetTrackID() << std::endl;
                  loc = photonProcess::PHANTOM;
                }

                else if (thePrePV == "physicalTreatmentRoom") {
                  //std::cout << "Photon absorbed in volume: " << thePrePV << " with a track id = " << track->GetTrackID() << std::endl;
                  loc = photonProcess::TREATMENT_ROOM;
                }

                else {
                  //std::cout << "Photon absorbed in ANOTHER volume: " << thePrePV << " with a track id = " << track->GetTrackID() << std::endl;
                  loc = photonProcess::OTHER;
                }
                
                photonProcess photon = photonProcess(track->GetTrackID(), genPosition.x()/mm, genPosition.y()/mm, genPosition.z()/mm, theta, phi, PhotonTotalInternalReflectionCount, PhotonReflectionCount, PhotonRefractionCount, loc) ; 
                //photon.print(); 
                photonProcess_vector.push_back(photon);

                trappedPhoton.erase(track->GetTrackID()); // Remove from the set as it’s now absorbed
                PhotonTotalInternalReflectionCount = 0; 
                PhotonRefractionCount = 0; 
                PhotonReflectionCount = 0; 
            }
        } /* end of if killed */
    }   /* end of if optical photon */

}

