#ifndef TSPGA_H
#define TSPGA_H

#include <vector>
#include <iostream>
#include <fstream>
#include "City.h"
#include "Individual.h"
#include "random.h"

class TSPGA {
private:
    int numCities;
    int populationSize;
    bool useL1;
    
    std::vector<City> cities;
    std::vector<Individual> population;
    Random rand;
    
    double pMutationPair;
    double pMutationShift;
    double pMutationSwap;
    double pMutationInversion;
    double pCrossover;
    
    // Aggiungi questa dichiarazione
    void initializePopulation();
    
public:
    TSPGA(int populationSize, bool useL1, int rank);
    
    void setCities(const std::vector<City>& cities_list);
    void initialize();
    void nextGeneration();
    void run(int generations);
    
    bool checkIndividual(const Individual& ind) const;
    double calculateFitness(const Individual& ind) const;
    
    void mutatePair(Individual& ind);
    void mutateShift(Individual& ind);
    void mutateSwap(Individual& ind);
    void mutateInversion(Individual& ind);
    Individual crossover(const Individual& parent1, const Individual& parent2);
    
    Individual selectParent();
    void evolve();
    
    void replaceWorstIndividual(const Individual& ind);
    
    Individual getBestIndividual() const;
    void printBest() const;
    
    void saveCities(const std::string& filename) const;
    void saveBestPath(const std::string& filename) const;
    void saveFitnessHistory(const std::string& filename) const;
    
    int randomInt(int min, int max);
    double randomDouble(double min, double max);
};

#endif