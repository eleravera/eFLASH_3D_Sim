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
tbb::concurrent_vector<detection> detection_vector_0;
tbb::concurrent_vector<detection> detection_vector_1;
tbb::concurrent_vector<detection> detection_vector_2;
tbb::concurrent_vector<detection> detection_vector_3;
tbb::concurrent_vector<detection> detection_vector_ge4;

int main(int argc, char **argv) {

  if (argc < 4) {
    G4cerr << "Usage: " << argv[0] << " <macro_file> <seed> <output_file>" << G4endl;
    return 1;
  }

  G4String macroFile = argv[1];
  int seed = std::stoi(argv[2]);
  G4String outputFileName = argv[3];

  G4Random::setTheSeed(seed);

  auto* runManager = new G4MTRunManager();
  G4int nThreads = 235;
  runManager->SetNumberOfThreads(nThreads);

  runManager->SetUserInitialization(new FlashDetectorConstruction);
  runManager->SetUserInitialization(new FlashPhysicsList);
  runManager->SetUserInitialization(new FlashActionInitialization);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Clear output vectors before run
  photonProcess_vector.clear();
  detection_vector.clear();
  detection_vector_0.clear();
  detection_vector_1.clear();
  detection_vector_2.clear();
  detection_vector_3.clear();
  detection_vector_ge4.clear();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4Timer timer;
  timer.Start();
  UImanager->ApplyCommand("/control/execute " + macroFile);
  timer.Stop();

  auto writeDetectionFile = [](const G4String& filename,
                               const tbb::concurrent_vector<detection>& vec)
  {
    std::ofstream fout(filename.c_str(), std::ios::binary);
    if (!fout.is_open()) {
      G4cerr << "Error: unable to open file " << filename << G4endl;
      return;
    }

    for (uint32_t i = 0; i < vec.size(); i++) {
      fout.write(reinterpret_cast<const char*>(&vec[i]), sizeof(detection));
    }
    fout.close();
  };

  // Remove .raw extension if present
  G4String outputStem = outputFileName;
  if (outputStem.size() >= 4 &&
      outputStem.substr(outputStem.size() - 4, 4) == ".raw") {
    outputStem = outputStem.substr(0, outputStem.size() - 4);
  }

  // Write all hits + separated files by number of reflections
  writeDetectionFile(outputStem + ".raw",      detection_vector);
  writeDetectionFile(outputStem + "_0.raw",    detection_vector_0);
  writeDetectionFile(outputStem + "_1.raw",    detection_vector_1);
  writeDetectionFile(outputStem + "_2.raw",    detection_vector_2);
  writeDetectionFile(outputStem + "_3.raw",    detection_vector_3);
  writeDetectionFile(outputStem + "_ge4.raw",  detection_vector_ge4);

  // Print summary
  std::cout << "Photons passing the selection: " << detection_vector.size() << std::endl;
  std::cout << "Hits with 0 reflections:  " << detection_vector_0.size() << std::endl;
  std::cout << "Hits with 1 reflection:   " << detection_vector_1.size() << std::endl;
  std::cout << "Hits with 2 reflections:  " << detection_vector_2.size() << std::endl;
  std::cout << "Hits with 3 reflections:  " << detection_vector_3.size() << std::endl;
  std::cout << "Hits with >=4 reflections:" << detection_vector_ge4.size() << std::endl;

  std::cout << "        ------      " << std::endl;
  std::cout << "Number of threads: " << runManager->GetNumberOfThreads() << std::endl;
  std::cout << "Elapsed time: " << timer.GetRealElapsed() << " seconds" << std::endl;
  std::cout << "outputFileName: " << outputFileName << std::endl;
  std::cout << "Seed: " << seed << std::endl;

  delete visManager;
  delete runManager;
  return 0;
}





//QUESTO MAIL È PER ME PER RUNNARE TEST SUL MIO PC
/*int main(int argc, char **argv) {

  auto *runManager=G4RunManagerFactory::CreateRunManager();
  G4int nThreads = 10;
  runManager->SetNumberOfThreads(nThreads);
 
  G4Random::setTheSeed(6003);

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
  detection_vector_0.clear();
  detection_vector_1.clear();
  detection_vector_2.clear();
  detection_vector_3.clear();
  detection_vector_ge4.clear();


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
    


  // Write results to output
  std::ofstream file_out1("./optical_properties/n_158_Rayleight/1m/seed6009_100evt.raw");
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
    file_out1.close();

    G4cout << "WRITE DEBUG: all = " << detection_vector.size() << G4endl;
    G4cout << "WRITE DEBUG: 0   = " << detection_vector_0.size() << G4endl;
    G4cout << "WRITE DEBUG: 1   = " << detection_vector_1.size() << G4endl;
    G4cout << "WRITE DEBUG: 2   = " << detection_vector_2.size() << G4endl;
    G4cout << "WRITE DEBUG: 3   = " << detection_vector_3.size() << G4endl;
    G4cout << "WRITE DEBUG: ge4 = " << detection_vector_ge4.size() << G4endl;

    // Write results to output
    std::ofstream file_out2("./photon_dist/pinhole/test_new_maps.raw");
    for (uint32_t i=0; i<detection_vector.size(); i++) {
      file_out2.write(reinterpret_cast<char*>(&detection_vector[i]), sizeof(detection));
    }
    file_out2.close();


    std::ofstream f0("./photon_dist/pinhole/test_new_maps_0.raw", std::ios::binary);
    for (uint32_t i = 0; i < detection_vector_0.size(); i++) {
        f0.write(reinterpret_cast<char*>(&detection_vector_0[i]), sizeof(detection));
    }
    f0.close();

  std::ofstream f1("./photon_dist/pinhole/test_new_maps_1.raw", std::ios::binary);
  for (uint32_t i = 0; i < detection_vector_1.size(); i++) {
      f1.write(reinterpret_cast<char*>(&detection_vector_1[i]), sizeof(detection));
  }
  f1.close();

  std::ofstream f2("./photon_dist/pinhole/test_new_maps_2.raw", std::ios::binary);
  for (uint32_t i = 0; i < detection_vector_2.size(); i++) {
      f2.write(reinterpret_cast<char*>(&detection_vector_2[i]), sizeof(detection));
  }
  f2.close();

  std::ofstream f3("./photon_dist/pinhole/test_new_maps_3.raw", std::ios::binary);
  for (uint32_t i = 0; i < detection_vector_3.size(); i++) {
      f3.write(reinterpret_cast<char*>(&detection_vector_3[i]), sizeof(detection));
  }
  f3.close();

  std::ofstream f4("./photon_dist/pinhole/test_new_maps_ge4.raw", std::ios::binary);
  for (uint32_t i = 0; i < detection_vector_ge4.size(); i++) {
      f4.write(reinterpret_cast<char*>(&detection_vector_ge4[i]), sizeof(detection));
  }
  f4.close();


  std::cout << "Elapsed time: " << timer.GetRealElapsed() << " seconds" << std::endl;

  delete visManager;
  delete runManager;
  return 0;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
