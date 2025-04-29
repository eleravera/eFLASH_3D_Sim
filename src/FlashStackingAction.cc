#include "FlashStackingAction.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include <common.hh>

long long FlashStackingAction::TotalPhotonGeneratedCount = 0;
float FlashStackingAction::PhotonEnergy;

FlashStackingAction::FlashStackingAction()
    : G4UserStackingAction() {}

FlashStackingAction::~FlashStackingAction() {}


G4ClassificationOfNewTrack FlashStackingAction::ClassifyNewTrack(const G4Track *aTrack)
  {
    const G4ParticleDefinition* particleDef = aTrack->GetDefinition();

    if (particleDef == G4OpticalPhoton::OpticalPhotonDefinition()) {
        TotalPhotonGeneratedCount++;
        PhotonEnergy = aTrack->GetKineticEnergy() / eV;
        //G4cout << "[Optical photon] Energy: " << PhotonEnergy << " eV" << G4endl;
    }    

  return fUrgent;
  }


  void FlashStackingAction::NewStage(){}

  void FlashStackingAction::PrepareNewEvent(){}