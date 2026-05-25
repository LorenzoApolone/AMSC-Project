#ifndef PARTICLE_HPP
#define PARTICLE_HPP
/**
 * @file particle.hpp
 * @brief Particle state used by the standard PSO implementations.
 */

#include <vector>
#include <limits> 

/**
 * @struct Particle
 * @brief Stores the mutable state of one PSO particle.
 *
 * A particle keeps its current position and velocity, together with the best
 * position it has found so far and the corresponding fitness value. The same
 * structure is used by the serial and MPI versions of the standard PSO solver.
 */
struct Particle {
    ///< Current position in the search space.
    std::vector<double> position;
    
    ///< Current velocity vector.
    std::vector<double> velocity;

    ///< Best position found by this particle.
    std::vector<double> best_position;
    
    ///< Fitness value at the current position.
    double current_value; 

    ///< Best fitness value found by this particle.
    double best_value;    

    /**
     * @brief Construct a particle with vectors sized to the problem dimension.
     *
     * @param dimension Number of coordinates in the search space.
     */
    Particle(int dimension) {
        position.resize(dimension);
        velocity.resize(dimension);
        best_position.resize(dimension);
        
        best_value = std::numeric_limits<double>::infinity(); 
    }
};

#endif
