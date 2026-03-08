#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <tbb/tbb.h>
#include "G4VUserTrackInformation.hh"


/* Detection of photons detected by the imaging system.
   The class has 3 attributes and a print member. 
*/
/*class detection {
    float x;
    float y;
    float z;

public:
    detection(float x_p, float y_p, float z_p) : x(x_p), y(y_p), z(z_p) {}

    void print() const {
        std::cout << "Detection coordinates:\n"
                  << "  x: " << std::fixed << std::setprecision(2) << x << "\n"
                  << "  y: " << std::fixed << std::setprecision(2) << y << "\n"
                  << "  z: " << std::fixed << std::setprecision(2) << z << "\n";
    }
};*/


class detection {
public:
    detection(double x_, double y_, double z_, int nReflections_)
        : x(x_), y(y_), z(z_), nReflections(nReflections_) {}

    double x;
    double y;
    double z;
    int nReflections;
};

/* photonProcess --.
    ....
*/
class photonProcess {
public:
    enum AbsorptionLocation {
        PHANTOM = 0,
        TREATMENT_ROOM = 1,
        PINHOLE = 2,
        DETECTOR = 3,
        OUT_OF_WORLD = 4,
        OTHER = 5
    };

    uint32_t event_id;
    float x;                       //Coordinate della scintillazione 
    float y;
    float z;
    float theta;                    //Angoli di emissione al punto di generazione
    float phi;
    float energy;
    uint32_t tirCount;            // Contatore per TIR = total internal reflection
    uint32_t reflectionCount;      // Contatore per riflessione 
    uint32_t refractionCount;      // Contatore per rifrazione
    uint32_t rayleighCount;   // Numero di scattering Rayleigh
    float firstRayleighTheta;
    bool hasRayleigh;
    AbsorptionLocation volume;
    float x_exit;                       //Coordinate all'interfaccia scintillatore/treatment room 
    float y_exit;
    float z_exit;
    float px_inside;                       //Coordinate del momento del fotone all'arrivo all'interfaccia
    float py_inside;
    float pz_inside;    
    float px_outside;                       //Coordinate del momento del fotone all'uscita dell'interfaccia
    float py_outside;
    float pz_outside;


    photonProcess(uint32_t e, float x_p, float y_p, float z_p, float theta_p, float phi_p, float en,
                uint32_t tir_count, uint32_t reflection_count, uint32_t refraction_count, 
                uint32_t ray_count, float first_rayleigh_theta, bool has_rayleigh,   // <-- AGGIUNTI QUI
                AbsorptionLocation v,
                float x_i, float y_i, float z_i,
                float px_i, float py_i, float pz_i, 
                float px_o, float py_o, float pz_o)
        : event_id(e),
        x(x_p), y(y_p), z(z_p), 
        theta(theta_p), phi(phi_p), energy(en), 
        tirCount(tir_count), reflectionCount(reflection_count),
        refractionCount(refraction_count), rayleighCount(ray_count),
        firstRayleighTheta(first_rayleigh_theta),
        hasRayleigh(has_rayleigh ? 1 : 0),        
        volume(v),
        x_exit(x_i), y_exit(y_i), z_exit(z_i), 
        px_inside(px_i), py_inside(py_i), pz_inside(pz_i),
        px_outside(px_o), py_outside(py_o), pz_outside(pz_o)
    {}


    void print() {
        std::cout << "Event ID: " << event_id << "\n"
                  << "Position of generation: (" << x << ", " << y << ", " << z << ") mm \n"
                  << "Angles emission: Theta = " << theta << ", Phi = " << phi << "\n"
                  << "Photon's energy = " << energy << " eV \n"
                  << "Total Internal Reflection Count: " << tirCount << "\n"
                  << "Reflection Count: " << reflectionCount << "\n"
                  << "Refraction Count: " << refractionCount << "\n"
                  << "Absorption Location: ";
                    switch (volume) {
                        case AbsorptionLocation::PHANTOM:        std::cout << "PHANTOM"; break;
                        case AbsorptionLocation::TREATMENT_ROOM: std::cout << "TREATMENT_ROOM"; break;
                        case AbsorptionLocation::PINHOLE:        std::cout << "PINHOLE"; break;
                        case AbsorptionLocation::DETECTOR:       std::cout << "DETECTOR"; break;
                        case AbsorptionLocation::OTHER:          std::cout << "OTHER"; break;
                        case AbsorptionLocation::OUT_OF_WORLD:   std::cout << "OUT_OF_WORLD"; break;
                    } 
        std::cout << "\n";
        std::cout << "Exit Position: (" << x_exit << ", " << y_exit << ", " << z_exit << ") mm \n"
              << "Momentum Inside:  (" << px_inside << ", " << py_inside << ", " << pz_inside << ")\n"
              << "Momentum Outside: (" << px_outside << ", " << py_outside << ", " << pz_outside << ")\n"
              << "Rayleigh count: " << rayleighCount << "\n";

        std::cout << std::endl << std::endl;
    }
};


extern tbb::concurrent_vector<detection> detection_vector;
extern tbb::concurrent_vector<photonProcess> photonProcess_vector;



class FlashPhotonTrackInfo : public G4VUserTrackInformation {
public:
    FlashPhotonTrackInfo() = default;
    virtual ~FlashPhotonTrackInfo() = default;

    int nInternalReflections = 0;

    bool exitedPhantomByRefraction = false;
};

#endif
