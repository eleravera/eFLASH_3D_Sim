#include "FlashTrackingAction.hh"
#include "FlashSteppingAction.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include <common.hh>


FlashTrackingAction::FlashTrackingAction(FlashSteppingAction* steppingAction)
    : G4UserTrackingAction(), fSteppingAction(steppingAction) {}

FlashTrackingAction::~FlashTrackingAction() {}


void FlashTrackingAction::PreUserTrackingAction(const G4Track*){

}

void FlashTrackingAction::PostUserTrackingAction(const G4Track* aTrack){

  /*const G4ParticleDefinition* particle = aTrack->GetParticleDefinition();
  
  if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
    const G4Step* aStep = aTrack->GetStep();
    G4StepPoint *postStep = aStep->GetPostStepPoint();
    G4StepPoint *preStep = aStep->GetPreStepPoint();
    photonProcess::AbsorptionLocation loc;
    G4ThreeVector genPosition = aTrack->GetVertexPosition();
    G4ThreeVector momentumDir = aTrack->GetVertexMomentumDirection();
    G4double theta = momentumDir.angle(G4ThreeVector(0, 1, 0)); // angle with the y-axis
    G4double phi = std::atan2(momentumDir.x(), momentumDir.z());

    if (postStep->GetPhysicalVolume()) {
      G4String PostvolumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
      G4String PreVolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

      if (PostvolumeName == "logicTreatmentRoom"){
        loc = photonProcess::TREATMENT_ROOM;
      }
      else if(PostvolumeName == "phantomLog"){
        loc = photonProcess::PHANTOM;
      }
      else if(PostvolumeName == "pinholeLog"){
        loc = photonProcess::PINHOLE;
      }
      else if(PostvolumeName == "DetectorLog"){
        loc == photonProcess::DETECTOR;
      }
      else {
        loc = photonProcess::OTHER;
      }

    } 
    else if(!(aStep->GetPostStepPoint()->GetPhysicalVolume())) {
      loc = photonProcess::OUT_OF_WORLD;
    }

    else {
      G4cout << "WARNING" << G4endl;
    }
    photonProcess photon = photonProcess(aTrack->GetTrackID(), genPosition.x()/mm, genPosition.y()/mm, genPosition.z()/mm, theta, phi, fSteppingAction->PhotonTotalInternalReflectionCount, fSteppingAction->PhotonReflectionCount, fSteppingAction->PhotonRefractionCount, loc) ; 
    //photon.print(); 
    photonProcess_vector.push_back(photon);
    fSteppingAction->PhotonTotalInternalReflectionCount = 0; 
    fSteppingAction->PhotonRefractionCount = 0; 
    fSteppingAction->PhotonReflectionCount = 0; 
    fSteppingAction->PhotonsOutOfWorld = 0 ;
  
  }*/
}
