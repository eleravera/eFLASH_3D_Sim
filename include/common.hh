#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <tbb/tbb.h>

/* Detection of photons detected by the imaging system.
   The class has 3 attributes and a print member. 
*/
class detection{
    float x;
    float y;
    float z; 
    
    public:
        detection(float x_p, float y_p, float z_p){
            x = x_p;
            y = y_p;
            z = z_p;
        }

    void print() {
        std::cout<< x << " " << y << " " << z << std::endl; 
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
        OTHER = 2
    };

    uint32_t event_id;
    float x;
    float y;
    float z;
    float theta;
    float phi;
    uint32_t tirCount;            // TIR = total internal reflection
    uint32_t reflectionCount;      // Contatore per riflessione
    uint32_t refractionCount;      // Contatore per rifrazione
    AbsorptionLocation volume;

    // Costruttore aggiornato per includere reflectionCount e refractionCount
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
            case OTHER:
                std::cout << "OTHER";
                break;
        }
        std::cout << std::endl << std::endl;
    }
};





extern tbb::concurrent_vector<detection> detection_vector1;
extern tbb::concurrent_vector<detection> detection_vector2;
extern tbb::concurrent_vector<photonProcess> photonProcess_vector;

#endif
