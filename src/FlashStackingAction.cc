#include "FlashStackingAction.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"

FlashStackingAction::FlashStackingAction()
    : G4UserStackingAction() {}

FlashStackingAction::~FlashStackingAction() {}


G4ClassificationOfNewTrack FlashStackingAction::ClassifyNewTrack(const G4Track *aTrack)
  {

    //const G4Track* track = aStep->GetTrack();
    const G4ParticleDefinition* particleDef = aTrack->GetDefinition();

    //std::cout << "particle = " << particleDef->GetParticleName() << " with GetParentID = " << aTrack->GetParentID() << " and  track->GetTrackID() = " << aTrack->GetTrackID() << std::endl;    

  return fUrgent;
  }


  void FlashStackingAction::NewStage(){}

  void FlashStackingAction::PrepareNewEvent(){}