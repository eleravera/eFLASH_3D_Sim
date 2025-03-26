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

G4int FlashSteppingAction::FresnelRefractionCount = 0 ; 
G4int FlashSteppingAction::FresnelReflectionCount = 0 ;
G4int FlashSteppingAction::TotalInternalReflectionCount = 0 ;
G4int FlashSteppingAction::AbsorptionCount = 0 ; 
G4int FlashSteppingAction::PhotonTotalInternalReflectionCount = 0 ;
G4int FlashSteppingAction::PhotonRefractionCount = 0 ;
G4int FlashSteppingAction::PhotonReflectionCount = 0 ;
long long FlashSteppingAction::TotalPhotonGeneratedCount = 0; 
long long FlashSteppingAction::PhotonExitingPhantomCount = 0; 



FlashSteppingAction::FlashSteppingAction(FlashEventAction *)
    : G4UserSteppingAction() {}

FlashSteppingAction::~FlashSteppingAction() {}



void FlashSteppingAction::HandleBoundaryProcesses(const G4Step* aStep, G4StepPoint* preStep, G4StepPoint* postStep) {
    G4Track* track = aStep->GetTrack();
    G4OpBoundaryProcessStatus theStatus = Undefined;
    static G4ThreadLocal G4OpBoundaryProcess* boundary = NULL;

    if (!boundary) {
        G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
        G4int nprocesses = pm->GetProcessListLength();
        G4ProcessVector* pv = pm->GetProcessList();
        for (G4int i = 0; i < nprocesses; i++) {
            if ((*pv)[i]->GetProcessName() == "OpBoundary") {
                boundary = (G4OpBoundaryProcess*)(*pv)[i];
                break;
            }
        }
    }
    
    if (postStep->GetStepStatus() == fGeomBoundary) {
        switch (boundary->GetStatus()) {
            case Absorption: AbsorptionCount++; break;
            case FresnelRefraction: FresnelRefractionCount++; PhotonRefractionCount++; break;
            case FresnelReflection: FresnelReflectionCount++; PhotonReflectionCount++; break;
            case TotalInternalReflection:
                TotalInternalReflectionCount++; PhotonTotalInternalReflectionCount++;
                if (PhotonTotalInternalReflectionCount > 10) {
                    track->SetTrackStatus(fStopAndKill);
                }
                break;
            default: break;
        }
    }
}


void FlashSteppingAction::CheckPhotonExit(const G4Step* aStep, G4StepPoint* preStep, G4StepPoint* postStep) {
    // Controllo se il passo è valido
    if (!aStep) {
        G4cerr << "Errore: Step non valido!" << G4endl;
        return;
    }

    // Recupero del fotone
    G4Track* track = aStep->GetTrack();
    if (!track) {
        G4cerr << "Errore: Track non valido!" << G4endl;
        return;
    }

    // Verifica dello stato del fotone (non deve essere fermato)
    if (track->GetTrackStatus() == fStopAndKill) {
        G4cerr << "Errore: Il fotone è fermato." << G4endl;
        return;
    }



    // Get the names of the volumes
    G4String volumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    G4String prevolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

    // Check if the photon is exiting the 'phantomLog' to the 'logicTreatmentRoom'
    if (prevolumeName == "phantomLog" && volumeName == "logicTreatmentRoom") {
        // Get the position where the photon exits the 'phantomLog' (before it enters 'logicTreatmentRoom')
        G4ThreeVector exitPosition = preStep->GetPosition();
        
        // Get the momentum of the photon before it exits (inside 'phantomLog')
        G4ThreeVector preExitMomentum = preStep->GetMomentum();
        
        // Get the momentum of the photon after it exits ('logicTreatmentRoom')
        G4ThreeVector postExitMomentum = aStep->GetTrack()->GetMomentum();

        // Log or save the positions and momentum
        G4cout << "Photon exited the phantom at position: "
               << exitPosition.x() / mm << " mm, "
               << exitPosition.y() / mm << " mm, "
               << exitPosition.z() / mm << " mm" << G4endl;

        G4cout << "Photon momentum before exiting (inside phantom): "
               << preExitMomentum.x() / GeV << " GeV/c, "
               << preExitMomentum.y() / GeV << " GeV/c, "
               << preExitMomentum.z() / GeV << " GeV/c" << G4endl;

        G4cout << "Photon momentum after exiting (inside treatment room): "
               << postExitMomentum.x() / GeV << " GeV/c, "
               << postExitMomentum.y() / GeV << " GeV/c, "
               << postExitMomentum.z() / GeV << " GeV/c" << G4endl;

        // Also print the information to standard output for convenience
        std::cout << "Exit Position: " << exitPosition.x() / mm << " " << exitPosition.y() / mm << " " << exitPosition.z() / mm << " mm\n";
        std::cout << "Momentum before exiting: " << preExitMomentum.x() / GeV << " " << preExitMomentum.y() / GeV << " " << preExitMomentum.z() / GeV << " GeV/c\n";
        std::cout << "Momentum after exiting: " << postExitMomentum.x() / GeV << " " << postExitMomentum.y() / GeV << " " << postExitMomentum.z() / GeV << " GeV/c\n"<< G4endl << G4endl;
        
        // Increment the count of photons exiting the phantom
        PhotonExitingPhantomCount++;
    }
}


void FlashSteppingAction::HandlePhotonDetection(const G4Step* aStep, G4StepPoint* preStep, G4StepPoint* postStep) {
    G4String preVolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    G4String postVolumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    if (preVolumeName == "logicTreatmentRoom" && postVolumeName == "DetectorLog") {
        G4ThreeVector position = aStep->GetTrack()->GetPosition();
        detection photon_maps(position.x() / mm, position.y() / mm, position.z() / mm);
        detection_vector.push_back(photon_maps);
        //photon_maps.print();
    }
}

void FlashSteppingAction::UserSteppingAction(const G4Step *aStep) {
    G4Track* track = aStep->GetTrack();
    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) { 
        TotalPhotonGeneratedCount++;
        G4StepPoint* preStep = aStep->GetPreStepPoint();
        G4StepPoint* postStep = aStep->GetPostStepPoint();
        //HandleBoundaryProcesses(aStep, preStep, postStep);

        //CheckPhotonExit(aStep, preStep, postStep); //- > per studiare Snell 
        //HandlePhotonDetection(aStep, preStep, postStep); //-> per salvare i dati e le mappe. 
    }
}


/*//Save photons with telecentric
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
            } */