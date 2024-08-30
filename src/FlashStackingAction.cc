#include "FlashStackingAction.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"

FlashStackingAction::FlashStackingAction()
    : G4UserStackingAction() {}

FlashStackingAction::~FlashStackingAction() {}


G4ClassificationOfNewTrack FlashStackingAction::ClassifyNewTrack(const G4Track *aTrack)
  {

    /*const G4ParticleDefinition* particleDef = aTrack->GetDefinition();

    G4ThreeVector photonDirection = aTrack->GetMomentumDirection();  

    // Define a small tolerance value
    const G4double tolerance = 1e-3;

    // Check if the photon is parallel to the x-axis
    if (!(std::abs(photonDirection.x()) > 1.0 - tolerance &&
        std::abs(photonDirection.y()) < tolerance &&
        std::abs(photonDirection.z()) < tolerance))
        {
          std::cout << "A photon passed the selection" << std::endl;
          // Kill the photon
          return fKill;
        }*/

  return fUrgent;
  }


  void FlashStackingAction::NewStage(){}

  void FlashStackingAction::PrepareNewEvent(){}