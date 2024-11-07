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
//  flash.cc
// Authors: J.Pensavalle, G.Milluzzo, F. Romano
// jake.pensavalle@pi.infn.it, giuliana.milluzzo@ct.infn.it, francesco.romano@ct.infn.it



#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "Applicator.hh"
#include "FlashActionInitialization.hh"
#include "FlashDetectorConstruction.hh"
#include "FlashPhysicsList.hh"
#include "G4ScoringManager.hh"

#include "Randomize.hh"
#include <common.hh>

// concurrent vector to write output in multithread mode without conflicts 
tbb::concurrent_vector<detection> detection_vector1;
tbb::concurrent_vector<detection> detection_vector2;

int main(int argc, char **argv) {

  auto *runManager=G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 1;
  runManager->SetNumberOfThreads(nThreads);
 
  G4Random::setTheSeed(45691);

  runManager->SetUserInitialization(new FlashDetectorConstruction);

  runManager->SetUserInitialization(new FlashPhysicsList);

  runManager->SetUserInitialization(new FlashActionInitialization);

  G4VisManager *visManager = new G4VisExecutive;

  visManager->Initialize();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  G4ScoringManager::GetScoringManager();
  

  // clears output vectors before run
  detection_vector1.clear();    
  detection_vector2.clear();

  G4UIExecutive *ui = 0;
    if (argc == 1) {
      ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute init_vis.mac");
      ui->SessionStart();
      delete ui;
    }
    else
      {
      G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);

    }

  // Write results to output
    
    std::ofstream file_out2("./photon_dist/photon_position_generation.raw");
    for (uint32_t i=0; i<detection_vector1.size(); i++) {
      file_out2.write(reinterpret_cast<char*>(&detection_vector1[i]), sizeof(detection));
    
      //std::cout<< reinterpret_cast<char*>(&detection_vector1[i]) << std::endl; 
    }

    file_out2.close();


    std::ofstream file_out3("./photon_dist/photon_angle_generation.raw");
    for (uint32_t i=0; i<detection_vector2.size(); i++) {
      file_out3.write(reinterpret_cast<char*>(&detection_vector2[i]), sizeof(detection));
    
      //std::cout<< reinterpret_cast<char*>(&detection_vector1[i]) << std::endl; 
    }

    file_out3.close();

  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
