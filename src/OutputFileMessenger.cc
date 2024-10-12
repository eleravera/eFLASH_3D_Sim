#include "OutputFileMessenger.hh"
#include "G4UImanager.hh"
//#include "G4Random.hh"

OutputFileMessenger::OutputFileMessenger() {
    // Comando per settare il file di output
    outputFileCmd = new G4UIcmdWithAString("/myApp/outputFile", this);
    outputFileCmd->SetGuidance("Set output file name.");
    outputFileCmd->SetParameterName("outputFile", false);  // "false" indica che è obbligatorio

    // Comando per settare il seed
    seedCmd = new G4UIcmdWithAnInteger("/myApp/setSeed", this);
    seedCmd->SetGuidance("Set random seed.");
    seedCmd->SetParameterName("seed", false);  // "false" indica che è obbligatorio
}

OutputFileMessenger::~OutputFileMessenger() {
    delete outputFileCmd;
    delete seedCmd;
}


void OutputFileMessenger::CheckParameters() const {
    if (outputFileName.empty()) {
        G4cerr << "Error: Output file name must be set using /myApp/outputFile" << G4endl;
        exit(1);  // Esce se non è stato specificato il file di output
    }

    if (seed == -1) {  // supponendo che -1 indichi un seed non inizializzato
        G4cerr << "Error: Random seed must be set using /myApp/setSeed" << G4endl;
        exit(1);  // Esce se non è stato specificato il seed
    }
}

void OutputFileMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == outputFileCmd) {
        outputFileName = newValue;
        G4cout << "Output file set to: " << outputFileName << G4endl;
    } else if (command == seedCmd) {
        seed = G4UIcommand::ConvertToInt(newValue);
        G4cout << "Seed set to: " << seed << G4endl;
    }
}

G4String OutputFileMessenger::GetOutputFileName() const {
    return outputFileName;
}
