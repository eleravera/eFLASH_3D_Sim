#include "FlashTrackingAction.hh"
#include "FlashSteppingAction.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include <common.hh>

FlashTrackingAction::FlashTrackingAction(FlashSteppingAction* steppingAction, FlashStackingAction* stackingAction)
    : G4UserTrackingAction(), 
      fSteppingAction(steppingAction),
      fStackingAction(stackingAction),
      count_phantom(0), 
      count_treatmentRoom(0),
      count_pinhole(0),
      count_detector(0),
      count_other(0),
      count_outOfWorld(0),
      fRayleighCount(0) {}

FlashTrackingAction::~FlashTrackingAction() {}

void FlashTrackingAction::PreUserTrackingAction(const G4Track* aTrack){
    if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        fRayleighCount = 0;
        FlashSteppingAction::FirstRayleighTheta = -1.0;
        FlashSteppingAction::HasRayleigh = false;
    }
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
        fStackingAction->PhotonEnergy,
        fSteppingAction->PhotonTotalInternalReflectionCount,
        fSteppingAction->PhotonReflectionCount,
        fSteppingAction->PhotonRefractionCount,
        fRayleighCount,
        FlashSteppingAction::FirstRayleighTheta,    
        FlashSteppingAction::HasRayleigh,           
        loc,
        fSteppingAction->PhantomExitingPosition.x() / mm,
        fSteppingAction->PhantomExitingPosition.y() / mm,
        fSteppingAction->PhantomExitingPosition.z() / mm,
        fSteppingAction->MomentumDirectionInside.x(),
        fSteppingAction->MomentumDirectionInside.y(),
        fSteppingAction->MomentumDirectionInside.z(),
        fSteppingAction->MomentumDirectionOutside.x(),
        fSteppingAction->MomentumDirectionOutside.y(),
        fSteppingAction->MomentumDirectionOutside.z()
    );
    photonProcess_vector.push_back(photon); //Save the class info into the vector 
    //photon.print();




    G4ThreeVector tmpExitPos = fSteppingAction->PhantomExitingPosition;
    G4ThreeVector tmpMomIn = fSteppingAction->MomentumDirectionInside;
    G4ThreeVector tmpMomOut = fSteppingAction->MomentumDirectionOutside;

    // Reset photon counts for the next photon
    fSteppingAction->PhotonTotalInternalReflectionCount = 0;
    fSteppingAction->PhotonRefractionCount = 0;
    fSteppingAction->PhotonReflectionCount = 0;

    fSteppingAction->PhantomExitingPosition.setX(-1000.);
    fSteppingAction->PhantomExitingPosition.setY(-1000.);
    fSteppingAction->PhantomExitingPosition.setZ(-1000.);

    fSteppingAction->MomentumDirectionInside.setX(-1000.);
    fSteppingAction->MomentumDirectionInside.setY(-1000.);
    fSteppingAction->MomentumDirectionInside.setZ(-1000.);

    fSteppingAction->MomentumDirectionOutside.setX(-1000.);
    fSteppingAction->MomentumDirectionOutside.setY(-1000.);
    fSteppingAction->MomentumDirectionOutside.setZ(-1000.);


}



void FlashTrackingAction::DetermineAbsorptionLocation(const G4Step* aStep, photonProcess::AbsorptionLocation& loc) {
    const G4StepPoint* postStep = aStep->GetPostStepPoint();
    if (postStep->GetPhysicalVolume()) {
        G4String PostvolumeName = postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

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

    const G4ParticleDefinition* particle = aTrack->GetParticleDefinition();

    if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
        const G4Step* aStep = aTrack->GetStep();
        const G4StepPoint* postStep = aStep->GetPostStepPoint();
        photonProcess::AbsorptionLocation loc;
        
        //ProcessPhotonData(aTrack, aStep); // Call the function to process photon data
        //DetermineAbsorptionLocation(aStep, loc); 
    
    }

}

