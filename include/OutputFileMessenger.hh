#ifndef OUTPUTFILEMESSENGER_HH
#define OUTPUTFILEMESSENGER_HH

#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "Randomize.hh"

class OutputFileMessenger : public G4UImessenger {
public:
    OutputFileMessenger();
    virtual ~OutputFileMessenger();

    virtual void SetNewValue(G4UIcommand* command, G4String newValue);

    G4String GetOutputFileName() const;

    // Funzione per controllare i parametri obbligatori
    void CheckParameters() const;

private:
    G4UIcmdWithAString* outputFileCmd;
    G4UIcmdWithAnInteger* seedCmd;
    G4String outputFileName;
    G4int seed;
};

#endif
