#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>

//path e fitness
struct Individual {
    std::vector<int> path;
    double fitness;
    
    Individual(const std::vector<int>& path = {}) : path(path), fitness(0.0) {}
};

#endif
