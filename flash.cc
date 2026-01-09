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
// flash.cc
// Authors: Eleonora Ravera, eleonora.ravera@phd.unipi.it
//Based on a GEANT4 Example


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
#include "G4Timer.hh"

#include "Randomize.hh"
#include <common.hh>

// concurrent vector to write output in multithread mode without conflicts 
tbb::concurrent_vector<photonProcess> photonProcess_vector;
tbb::concurrent_vector<detection> detection_vector;

int main(int argc, char **argv) {

  if (argc < 4) {
        G4cerr << "Usage: " << argv[0] << " <macro_file> <seed> <output_file>" << G4endl;
        return 1;
    }

  G4String macroFile = argv[1];
  int seed = std::stoi(argv[2]);
  G4String outputFileName = argv[3];

  G4Random::setTheSeed(seed);
  auto *runManager = new G4MTRunManager();
  G4int nThreads = 235;
  runManager->SetNumberOfThreads(nThreads);
 
  runManager->SetUserInitialization(new FlashDetectorConstruction);
  runManager->SetUserInitialization(new FlashPhysicsList);
  runManager->SetUserInitialization(new FlashActionInitialization);

  G4VisManager *visManager = new G4VisExecutive;
  visManager->Initialize();

  // clears output vectors before run
  detection_vector.clear();  
  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  G4Timer timer;
  timer.Start();  
  UImanager->ApplyCommand("/control/execute " + macroFile);
  timer.Stop();

  std::ofstream file_out2(outputFileName.c_str());
  if (!file_out2.is_open()) {
    G4cerr << "Error: unable to open file " << outputFileName << G4endl;
  }
  else{
  // Write results to output
  for (uint32_t i=0; i<detection_vector.size(); i++) {
      file_out2.write(reinterpret_cast<char*>(&detection_vector[i]), sizeof(detection));
    }
  }
  
  //print some interesting information
  std::cout<< "Photons passing the selection: " << detection_vector.size() << std::endl;
  file_out2.close();

  std::cout << "        ------      " << std::endl;
  std::cout << "Number of threds: " << runManager->GetNumberOfThreads() << std::endl;
  std::cout << "Elapsed time: " << timer.GetRealElapsed() << " seconds" << std::endl;
  std::cout<<"outputFileName: " << outputFileName <<std::endl;
  std::cout<<"Seed: " << seed << std::endl;

  delete visManager;
  delete runManager;
  return 0;
}





//QUESTO MAIL È PER ME PER RUNNARE TEST SUL MIO PC
/*int main(int argc, char **argv) {

  auto *runManager=G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 1;
  runManager->SetNumberOfThreads(nThreads);
 
  G4Random::setTheSeed(6009);

  runManager->SetUserInitialization(new FlashDetectorConstruction);

  runManager->SetUserInitialization(new FlashPhysicsList);

  runManager->SetUserInitialization(new FlashActionInitialization);

  G4VisManager *visManager = new G4VisExecutive;

  visManager->Initialize();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  G4ScoringManager::GetScoringManager();
  

  G4Timer timer;
  timer.Start();  

  // clears output vectors before run
  photonProcess_vector.clear();    
  detection_vector.clear();    


  G4UIExecutive *ui = 0;
    if (argc == 1) {
      ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute init_vis.mac");
      ui->SessionStart(); // If you want to start the interactive session
      delete ui;
    }
    else
      {
      G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);

    }
    
    timer.Stop();
    */


  // Write results to output
    /*std::ofstream file_out1("./optical_properties/n_158_Rayleight/test.raw");
    for (const photonProcess& p : photonProcess_vector) {

        // PART 1 – campi base del fotone
        file_out1.write(reinterpret_cast<const char*>(&p.event_id), sizeof(uint32_t));

        file_out1.write(reinterpret_cast<const char*>(&p.x), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.y), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.z), sizeof(float));

        file_out1.write(reinterpret_cast<const char*>(&p.theta), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.phi), sizeof(float));

        file_out1.write(reinterpret_cast<const char*>(&p.energy), sizeof(float));

        // PART 2 – contatori
        file_out1.write(reinterpret_cast<const char*>(&p.tirCount), sizeof(uint32_t));
        file_out1.write(reinterpret_cast<const char*>(&p.reflectionCount), sizeof(uint32_t));
        file_out1.write(reinterpret_cast<const char*>(&p.refractionCount), sizeof(uint32_t));
        file_out1.write(reinterpret_cast<const char*>(&p.rayleighCount), sizeof(uint32_t));

        // PART 3 – Rayleigh extra info
        file_out1.write(reinterpret_cast<const char*>(&p.firstRayleighTheta), sizeof(float));

        uint32_t hasR = p.hasRayleigh ? 1 : 0;
        file_out1.write(reinterpret_cast<const char*>(&hasR), sizeof(uint32_t));

        // PART 4 – absorption location (enum → uint32)
        uint32_t v = static_cast<uint32_t>(p.volume);
        file_out1.write(reinterpret_cast<const char*>(&v), sizeof(uint32_t));

        // PART 5 – exit position
        file_out1.write(reinterpret_cast<const char*>(&p.x_exit), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.y_exit), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.z_exit), sizeof(float));

        // PART 6 – momentum inside
        file_out1.write(reinterpret_cast<const char*>(&p.px_inside), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.py_inside), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.pz_inside), sizeof(float));

        // PART 7 – momentum outside
        file_out1.write(reinterpret_cast<const char*>(&p.px_outside), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.py_outside), sizeof(float));
        file_out1.write(reinterpret_cast<const char*>(&p.pz_outside), sizeof(float));
    }
    file_out1.close();*/


    // Write results to output - SEMPRE PER RUNNURE SUL MIO PC
    /*std::ofstream file_out2("./photon_dist/test.raw");
    for (uint32_t i=0; i<detection_vector.size(); i++) {
      file_out2.write(reinterpret_cast<char*>(&detection_vector[i]), sizeof(detection));
    }
    file_out2.close();


  std::cout << "Elapsed time: " << timer.GetRealElapsed() << " seconds" << std::endl;

  delete visManager;
  delete runManager;
  return 0;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
