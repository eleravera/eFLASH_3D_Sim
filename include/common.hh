#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <tbb/tbb.h>

/* Detection of photons detected by the imaging system.
   The class has 3 attributes and a print member. 
*/
class detection {
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
    uint32_t tirCount;            // Contatore per TIR = total internal reflection
    uint32_t reflectionCount;      // Contatore per riflessione 
    uint32_t refractionCount;      // Contatore per rifrazione
    AbsorptionLocation volume;

    photonProcess(uint32_t e, float x_p, float y_p, float z_p, float theta_p, float phi_p, 
                  uint32_t tir_count, uint32_t reflection_count, uint32_t refraction_count, 
                  AbsorptionLocation v)
        : event_id(e), x(x_p), y(y_p), z(z_p), theta(theta_p), phi(phi_p), 
          tirCount(tir_count), reflectionCount(reflection_count), refractionCount(refraction_count), 
          volume(v) {}

    void print() {
        std::cout << "Event ID: " << event_id << "\n"
                  << "Position: (" << x << ", " << y << ", " << z << ")\n"
                  << "Angles: Theta = " << theta << ", Phi = " << phi << "\n"
                  << "Total Internal Reflection Count: " << tirCount << "\n"
                  << "Reflection Count: " << reflectionCount << "\n"
                  << "Refraction Count: " << refractionCount << "\n"
                  << "Absorption Location: ";

        // Stampa il valore corrispondente di absorption location
        switch(volume) {
            case PHANTOM:
                std::cout << "PHANTOM";
                break;
            case TREATMENT_ROOM:
                std::cout << "TREATMENT_ROOM";
                break;
            case PINHOLE:
                std::cout << "PINHOLE";
                break; 
            case DETECTOR:
                std::cout << "DETECTOR";
                break;
            case OUT_OF_WORLD:
                std::cout << "OUT_OF_WORLD";
                break;
            case OTHER:
                std::cout << "OTHER";
                break;
        }
        std::cout << std::endl << std::endl;
    }
};





extern tbb::concurrent_vector<detection> detection_vector;
extern tbb::concurrent_vector<photonProcess> photonProcess_vector;

#endif
