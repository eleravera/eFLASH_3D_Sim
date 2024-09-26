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
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include <G4Timer.hh>

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
// currently multithread not active
tbb::concurrent_vector<detection> detection_vector1;
tbb::concurrent_vector<detection> detection_vector2;

int main(int argc, char **argv) {

  

  //  G4Random::setTheEngine(new CLHEP::MTwistEngine);

  //auto *runManager=G4RunManagerFactory::CreateRunManager();
  auto *runManager = new G4MTRunManager();
  G4int nThreads = 6;
  runManager->SetNumberOfThreads(nThreads);
 
  G4Random::setTheSeed(45698);

  runManager->SetUserInitialization(new FlashDetectorConstruction);

  runManager->SetUserInitialization(new FlashPhysicsList);

  runManager->SetUserInitialization(new FlashActionInitialization);



  G4VisManager *visManager = new G4VisExecutive;

  visManager->Initialize();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  G4ScoringManager::GetScoringManager();
  

  // clears output vectors before run
  detection_vector1.clear();    


  G4Timer timer;
  timer.Start();
  G4UIExecutive *ui = 0;
    if (argc == 1) {
      ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute init_vis.mac");
      ui->SessionStart(); //per lanciare simulazione senza UI basta disalibilitare questa parte. 
      delete ui;
    }
    else
      {
      G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);

    }

  //runManager->BeamOn(100);
  timer.Stop();
  
  // Write results to output
  std::ofstream file_out2("./photon_dist/telecentric/photon_maps_01.raw");
  for (uint32_t i=0; i<detection_vector1.size(); i++) {
      file_out2.write(reinterpret_cast<char*>(&detection_vector1[i]), sizeof(detection));
    }
  
  std::cout<< "detection_vector1.size(): " << detection_vector1.size() << std::endl;
  file_out2.close();


  std::cout << "Number of threds: " << runManager->GetNumberOfThreads() << std::endl;
  std::cout << "Elapsed time: " << timer.GetRealElapsed() << " seconds" << std::endl;



  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
