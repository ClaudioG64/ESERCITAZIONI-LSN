#ifndef TSPGA_H
#define TSPGA_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib> // Per std::abs
#include "City.h"
#include "Individual.h"
#include "random.h"

class TSPGA {
private:
    int numCities;
    int populationSize;
    bool citiesOnCircle;
    bool useL1;
    
    std::vector<City> cities;
    std::vector<Individual> population;
    Random rand;
    
    double pMutationPair;
    double pMutationShift;
    double pMutationSwap;
    double pMutationInversion;
    double pCrossover;
    
public:
    TSPGA(int numCities, int populationSize, bool citiesOnCircle = true, bool useL1 = false);
    
    void initialize();
    void run(int generations);
    
    void generateCities();
    void initializePopulation();
    bool checkIndividual(const Individual& ind) const;
    double calculateFitness(const Individual& ind) const;
    
    void mutatePair(Individual& ind);
    void mutateShift(Individual& ind);
    void mutateSwap(Individual& ind);
    void mutateInversion(Individual& ind);
    Individual crossover(const Individual& parent1, const Individual& parent2);
    
    Individual selectParent();
    void evolve();
    
    void printBest() const;
    Individual getBestIndividual() const;
    
    void saveCities(const std::string& filename) const;
    void saveBestPath(const std::string& filename) const;
    void saveFitnessHistory(const std::string& filename) const;
    
    int randomInt(int min, int max);
    double randomDouble(double min, double max);
};

#endif