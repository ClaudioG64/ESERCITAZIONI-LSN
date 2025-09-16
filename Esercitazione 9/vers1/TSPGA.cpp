#include "TSPGA.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <ctime>
#include <fstream>
#include <iostream>
#include <functional>
#include <cstdlib>

TSPGA::TSPGA(int numCities, int populationSize, bool citiesOnCircle, bool useL1)
    : numCities(numCities), populationSize(populationSize),
      citiesOnCircle(citiesOnCircle), useL1(useL1),
      pMutationPair(0.1), pMutationShift(0.08),
      pMutationSwap(0.07), pMutationInversion(0.09),
      pCrossover(0.6) {
    
    int seed[4];
    int p1 = 0, p2 = 0;
    std::time_t now = std::time(0);
    seed[0] = static_cast<int>(now % 10000);
    seed[1] = static_cast<int>((now / 10000) % 10000);
    seed[2] = static_cast<int>((now / 1000000) % 10000);
    seed[3] = static_cast<int>((now / 100000000) % 10000);
    p1 = static_cast<int>(now % 1000);
    p2 = static_cast<int>((now / 1000) % 1000);
    rand.SetRandom(seed, p1, p2);
}

void TSPGA::initialize() {
    generateCities();
    initializePopulation();
    for (const auto& ind : population) {
        if (!checkIndividual(ind)) {
            std::cerr << "Error: Invalid individual in population!" << std::endl;
            exit(1);
        }
    }
    std::cout << "✓ All individuals are valid!" << std::endl;
}

int TSPGA::randomInt(int min, int max) {
    return static_cast<int>(rand.Rannyu(min, max + 1));
}

double TSPGA::randomDouble(double min, double max) {
    return rand.Rannyu(min, max);
}

void TSPGA::generateCities() {
    cities.clear();
    if (citiesOnCircle) {
        for (int i = 0; i < numCities; ++i) {
            double angle = 2.0 * M_PI * randomDouble(0.0, 1.0);
            cities.emplace_back(std::cos(angle), std::sin(angle));
        }
    } else {
        for (int i = 0; i < numCities; ++i) {
            cities.emplace_back(randomDouble(0.0, 1.0), randomDouble(0.0, 1.0));
        }
    }
}

void TSPGA::initializePopulation() {
    population.clear();
    for (int i = 0; i < populationSize; ++i) {
        std::vector<int> path(numCities);
        path[0] = 1;
        std::vector<int> otherCities;
        for (int j = 2; j <= numCities; ++j) {
            otherCities.push_back(j);
        }
        for (int k = otherCities.size() - 1; k > 0; --k) {
            int j = randomInt(0, k);
            std::swap(otherCities[k], otherCities[j]);
        }
        for (int j = 1; j < numCities; ++j) {
            path[j] = otherCities[j - 1];
        }
        population.emplace_back(path);
    }
}

bool TSPGA::checkIndividual(const Individual& ind) const {
    if (ind.path.size() != numCities) return false;
    if (ind.path[0] != 1) return false;
    std::vector<bool> visited(numCities + 1, false);
    for (int city : ind.path) {
        if (city < 1 || city > numCities || visited[city]) return false;
        visited[city] = true;
    }
    return true;
}

double TSPGA::calculateFitness(const Individual& ind) const {
    double totalDistance = 0.0;
    for (int i = 0; i < numCities; ++i) {
        int city1 = ind.path[i] - 1;
        int city2 = ind.path[(i + 1) % numCities] - 1;
        double dx = cities[city1].x - cities[city2].x;
        double dy = cities[city1].y - cities[city2].y;
        if (useL1) {
            totalDistance += std::abs(dx) + std::abs(dy);
        } else {
            totalDistance += dx * dx + dy * dy;
        }
    }
    return totalDistance;
}

void TSPGA::mutatePair(Individual& ind) {
    if (randomDouble(0.0, 1.0) < pMutationPair) {
        int idx1 = randomInt(1, numCities - 1);
        int idx2 = randomInt(1, numCities - 1);
        
        if (idx1 != idx2) {
            std::swap(ind.path[idx1], ind.path[idx2]);
        }
    }
}

void TSPGA::mutateShift(Individual& ind) {
    if (randomDouble(0.0, 1.0) < pMutationShift) {
        int start = randomInt(1, numCities - 2);
        int length = randomInt(1, numCities - start - 1);
        int shift = randomInt(1, numCities - start - length - 1); // -1 per evitare overflow
        
        std::vector<int> block(ind.path.begin() + start, ind.path.begin() + start + length);
        ind.path.erase(ind.path.begin() + start, ind.path.begin() + start + length);
        ind.path.insert(ind.path.begin() + start + shift, block.begin(), block.end());
    }
}

void TSPGA::mutateSwap(Individual& ind) {
    if (randomDouble(0.0, 1.0) < pMutationSwap) {
        int start1 = randomInt(1, numCities - 2);
        int length1 = randomInt(1, std::min(3, numCities - start1 - 1));
        
        int start2 = randomInt(1, numCities - 2);
        int max_attempts = 100;
        int attempts = 0;
        
        // Assicurati che i blocchi non si sovrappongano
        while (abs(start2 - start1) < length1 && attempts < max_attempts) {
            start2 = randomInt(1, numCities - 2);
            attempts++;
        }
        
        if (attempts >= max_attempts) return; // Non trovare un blocco valido
        
        int length2 = randomInt(1, std::min(3, numCities - start2 - 1));
        
        // Estrai i blocchi
        std::vector<int> block1(ind.path.begin() + start1, ind.path.begin() + start1 + length1);
        std::vector<int> block2(ind.path.begin() + start2, ind.path.begin() + start2 + length2);
        
        // Scambia i blocchi
        for (int i = 0; i < length1; i++) {
            ind.path[start1 + i] = block2[i];
        }
        for (int i = 0; i < length2; i++) {
            ind.path[start2 + i] = block1[i];
        }
    }
}

void TSPGA::mutateInversion(Individual& ind) {
    if (randomDouble(0.0, 1.0) < pMutationInversion) {
        int start = randomInt(1, numCities - 2);
        int end = randomInt(start + 1, numCities - 1);
        
        std::reverse(ind.path.begin() + start, ind.path.begin() + end + 1);
    }
}

Individual TSPGA::crossover(const Individual& parent1, const Individual& parent2) {
    int cutPoint = randomInt(1, numCities - 2);
    std::vector<int> childPath(parent1.path.begin(), parent1.path.begin() + cutPoint);
    
    // Aggiungi le città mancanti nell'ordine del secondo genitore
    for (int city : parent2.path) {
        if (std::find(childPath.begin(), childPath.end(), city) == childPath.end()) {
            childPath.push_back(city);
        }
    }
    
    return Individual(childPath);
}

Individual TSPGA::selectParent() {
    double r = randomDouble(0.0, 1.0);
    int idx = static_cast<int>(populationSize * std::pow(r, 2));
    return population[std::min(idx, populationSize - 1)];
}

void TSPGA::evolve() {
    std::vector<Individual> newPopulation;
    newPopulation.push_back(population[0]);
    while (newPopulation.size() < populationSize) {
        if (randomDouble(0.0, 1.0) < pCrossover) {
            Individual parent1 = selectParent();
            Individual parent2 = selectParent();
            Individual child = crossover(parent1, parent2);
            mutatePair(child);
            mutateShift(child);
            mutateSwap(child);
            mutateInversion(child);
            if (checkIndividual(child)) {
                newPopulation.push_back(child);
            }
        } else {
            Individual individual = selectParent();
            mutatePair(individual);
            mutateShift(individual);
            mutateSwap(individual);
            mutateInversion(individual);
            if (checkIndividual(individual)) {
                newPopulation.push_back(individual);
            }
        }
    }
    population = newPopulation;
}

void TSPGA::run(int generations) {
    for (auto& ind : population) {
        ind.fitness = calculateFitness(ind);
    }
    std::sort(population.begin(), population.end(), 
              [](const Individual& a, const Individual& b) {
                  return a.fitness < b.fitness;
              });
    std::ofstream fitnessFile("fitness_history.dat");
    for (int gen = 0; gen < generations; ++gen) {
        evolve();
        for (auto& ind : population) {
            ind.fitness = calculateFitness(ind);
        }
        std::sort(population.begin(), population.end(), 
                  [](const Individual& a, const Individual& b) {
                      return a.fitness < b.fitness;
                  });
        double avgFitness = 0.0;
        for (int i = 0; i < populationSize / 2; ++i) {
            avgFitness += population[i].fitness;
        }
        avgFitness /= (populationSize / 2);
        fitnessFile << gen << " " << population[0].fitness << " " << avgFitness << std::endl;
        if (gen % 100 == 0) {
            std::cout << "Generation " << gen << ": Best = " << population[0].fitness 
                      << ", Avg = " << avgFitness << std::endl;
        }
    }
    fitnessFile.close();
}

void TSPGA::saveCities(const std::string& filename) const {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& city : cities) {
            file << city.x << " " << city.y << std::endl;
        }
        file.close();
        std::cout << "Cities saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void TSPGA::saveBestPath(const std::string& filename) const {
    Individual best = getBestIndividual();
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int cityIndex : best.path) {
            file << cities[cityIndex-1].x << " " << cities[cityIndex-1].y << std::endl;
        }
        file << cities[best.path[0]-1].x << " " << cities[best.path[0]-1].y << std::endl;
        file.close();
        std::cout << "Best path saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void TSPGA::saveFitnessHistory(const std::string& filename) const {
    std::ifstream inFile("fitness_history.dat");
    std::ofstream outFile(filename);
    if (inFile.is_open() && outFile.is_open()) {
        outFile << inFile.rdbuf();
        inFile.close();
        outFile.close();
        std::cout << "Fitness history saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to copy fitness history file" << std::endl;
    }
}

Individual TSPGA::getBestIndividual() const {
    return population[0];
}

void TSPGA::printBest() const {
    Individual best = getBestIndividual();
    std::cout << "Best path: ";
    for (int city : best.path) {
        std::cout << city << " ";
    }
    std::cout << "\nFitness: " << best.fitness << std::endl;
}