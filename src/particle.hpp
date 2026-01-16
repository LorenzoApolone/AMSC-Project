#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <limits> 

struct Particle {
    // Position
    std::vector<double> position;
    
    // Velocity 
    std::vector<double> velocity;

    // Best position achieved
    std::vector<double> best_position;
    
    // Current f(x)
    double current_value; 

    // Best f(x)
    double best_value;    

    // Constructor
    Particle(int dimension) {
        position.resize(dimension);
        velocity.resize(dimension);
        best_position.resize(dimension);
        
        best_value = std::numeric_limits<double>::infinity(); 
    }
};

#endif