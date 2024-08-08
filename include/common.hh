#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <tbb/tbb.h>


class detection{
    //uint32_t event_id;
    //float event_id;
    float x;
    float y;
    float z; 
    //float cos_x; 
    //float cos_y;
    //float cos_z; 
    
public:
    //detection(float e, float x_p, float y_p, float z_p)
    detection(float x_p, float y_p, float z_p)

    {
        //event_id=e;
        x = x_p;
        y = y_p;
        z = z_p;
        //cosx = cosx_p; 
        //cosy = cosy_p; 
        //cosz = cosz_p; 

    }

void print() {
    std::cout<< x << " " << y << " " << z << std::endl; 
}
};

extern tbb::concurrent_vector<detection> detection_vector1;
extern tbb::concurrent_vector<detection> detection_vector2;

#endif
