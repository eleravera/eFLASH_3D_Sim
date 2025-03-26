#include "FlashTrackingAction.hh"
#include "FlashSteppingAction.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include <common.hh>




FlashTrackingAction::FlashTrackingAction(FlashSteppingAction* steppingAction)
    : G4UserTrackingAction(), 
      fSteppingAction(steppingAction),
      count_phantom(0), 
      count_treatmentRoom(0),
      count_pinhole(0),
      count_detector(0),
      count_wrap(0),
      count_other(0),
      count_outOfWorld(0) {}


FlashTrackingAction::~FlashTrackingAction() {}


void FlashTrackingAction::PreUserTrackingAction(const G4Track*){

}


void FlashTrackingAction::ProcessPhotonData(const G4Track* aTrack, const G4Step* aStep) {
 /*This function is called once at the end of the tracking of a photon, when the track is completed (i.e., the photon is absorbed, exits the world, or is otherwise terminated). It first checks if the particle is an optical photon. 
  If the particle is an optical photon, the function retrieves information about the photon's current step (pre-step and post-step).
  It then calculates the photon's generation position and momentum direction, from which the angles theta and phi are derived.
  
  The function checks if the photon has entered a new physical volume, indicating it may have interacted with an interface (e.g., between materials or volumes).
  Based on the name of the post-step volume, it determines the photon’s absorption location (such as "TREATMENT_ROOM", "PHANTOM", "PINHOLE", "DETECTOR", or "OTHER"). If the photon is no longer inside the simulation world, its location is set to "OUT_OF_WORLD".
  
  A new photonProcess object is then created to store the photon’s data, including its track ID, position, direction, reflection counts, and absorption location.
  After recording the data, reflection, refraction, and total internal reflection counts are reset to 0.
  Finally, the photonProcess object is added to the photonProcess_vector for further processing, and the photon counts are reset for the next photon.*/


    const G4StepPoint* postStep = aStep->GetPostStepPoint();
    const G4StepPoint* preStep = aStep->GetPreStepPoint();

    photonProcess::AbsorptionLocation loc;
    G4ThreeVector genPosition = aTrack->GetVertexPosition();
    G4ThreeVector momentumDir = aTrack->GetVertexMomentumDirection();
    G4double theta = momentumDir.angle(G4ThreeVector(0, 1, 0)); // angle with the y-axis
    G4double phi = std::atan2(momentumDir.x(), momentumDir.z());

    // Determine the absorption location based on the post-step volume
    if (postStep->GetPhysicalVolume()) {
        G4String PostvolumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
        G4String PreVolumeName = preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

        if (PostvolumeName == "logicTreatmentRoom") {
            loc = photonProcess::TREATMENT_ROOM;
        } else if (PostvolumeName == "phantomLog") {
            loc = photonProcess::PHANTOM;
        } else if (PostvolumeName == "pinholeLog") {
            loc = photonProcess::PINHOLE;
        } else if (PostvolumeName == "DetectorLog") {
            loc = photonProcess::DETECTOR;
        } else {
            loc = photonProcess::OTHER;
        }
    } else if (!(aStep->GetPostStepPoint()->GetPhysicalVolume())) {
        loc = photonProcess::OUT_OF_WORLD;
    } else {
        G4cout << "WARNING: No valid post-step volume found!" << G4endl;
        loc = photonProcess::OTHER;
    }

    // Create and store the photon process data
    photonProcess photon = photonProcess(
        aTrack->GetTrackID(),
        genPosition.x() / mm,
        genPosition.y() / mm,
        genPosition.z() / mm,
        theta,
        phi,
        fSteppingAction->PhotonTotalInternalReflectionCount,
        fSteppingAction->PhotonReflectionCount,
        fSteppingAction->PhotonRefractionCount,
        loc
    );
    
    photonProcess_vector.push_back(photon);

    // Reset photon counts for the next photon
    fSteppingAction->PhotonTotalInternalReflectionCount = 0;
    fSteppingAction->PhotonRefractionCount = 0;
    fSteppingAction->PhotonReflectionCount = 0;


    //check per sapere dove viene assorbito il fotone
    if (fSteppingAction->PhotonTotalInternalReflectionCount==0 &&  fSteppingAction->PhotonReflectionCount ==0   &&  fSteppingAction->PhotonRefractionCount==0){
        G4ThreeVector finalPosition = postStep->GetPosition();
        G4double x_final_mm = finalPosition.x() / mm;
        G4double y_final_mm = finalPosition.y() / mm;
        G4double z_final_mm = finalPosition.z() / mm;

        std::cout << "All counts = 0; final loc: " << loc << std::endl
          << "Gen position: " << genPosition.x() / mm << " " 
          << genPosition.y() / mm << " " 
          << genPosition.z() / mm << std::endl
          << "Final position: " << x_final_mm << " " 
          << y_final_mm << " " 
          << z_final_mm << std::endl 
          << std::endl;

    }
}


void FlashTrackingAction::DetermineAbsorptionLocation(const G4Step* aStep, photonProcess::AbsorptionLocation& loc) {
    const G4StepPoint* postStep = aStep->GetPostStepPoint();
    if (postStep->GetPhysicalVolume()) {
        G4String PostvolumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

        G4cout << "Photon is in volume (dentro la funzione): " << PostvolumeName << G4endl;


        if (PostvolumeName == "phantomLog") {
            loc = photonProcess::PHANTOM;
            count_phantom++;  // Increment counter for phantom
        } else if (PostvolumeName == "logicTreatmentRoom") {
            loc = photonProcess::TREATMENT_ROOM;
            count_treatmentRoom++;  // Increment counter for treatment room
        } else if (PostvolumeName == "pinholeLog") {
            loc = photonProcess::PINHOLE;
            count_pinhole++;  // Increment counter for pinhole
        } else if (PostvolumeName == "DetectorLog") {
            loc = photonProcess::DETECTOR;
            count_detector++;  // Increment counter for detector
        } else if (PostvolumeName == "WrapLog") {
            loc = photonProcess::WRAP;
            count_wrap++;  // Increment counter for wrap
        } else {
            loc = photonProcess::OTHER;
            count_other++;  // Increment counter for other
        }
    } else {
        loc = photonProcess::OUT_OF_WORLD;
        count_outOfWorld++;  // Increment counter for photons out of world
    }
}




void FlashTrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
    /*const G4ParticleDefinition* particle = aTrack->GetParticleDefinition();

    if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
        const G4Step* aStep = aTrack->GetStep();
        const G4StepPoint* postStep = aStep->GetPostStepPoint();
        photonProcess::AbsorptionLocation loc;
        // Call the refactored function to process photon data
        ProcessPhotonData(aTrack, aStep);
        // Optionally, you can include any additional processing logic here if needed
        DetermineAbsorptionLocation(aStep, loc);
    }

  G4cout << "Photons absorbed in Phantom: " << count_phantom << G4endl;
  G4cout << "Photons absorbed in Treatment Room: " << count_treatmentRoom << G4endl;
  G4cout << "Photons absorbed in Pinhole: " << count_pinhole << G4endl;
  G4cout << "Photons absorbed in Detector: " << count_detector << G4endl;
  G4cout << "Photons absorbed in Wrap: " << count_wrap << G4endl;
  G4cout << "Photons absorbed in Other location: " << count_other << G4endl;
  G4cout << "Photons out of world: " << count_outOfWorld << G4endl<<G4endl;
  
  */

}