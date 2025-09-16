#include <iostream>
#include "TSPGA.h"


int main() {
    std::cout << "=== TSP Genetic Algorithm with 34 Cities ===" << std::endl;
    
    // Esperimento 1: 34 città sulla circonferenza
    std::cout << "\n1. 34 cities on circumference:" << std::endl;
    TSPGA tsp1(34, 100, true, false);
    tsp1.initialize();
    tsp1.saveCities("cities_circle.dat");  // Salva le città
    tsp1.run(1000);
    tsp1.saveBestPath("best_path_circle.dat");  // Salva il percorso migliore
    tsp1.saveFitnessHistory("fitness_circle.dat");  // Salva la storia della fitness
    
    // Esperimento 2: 34 città dentro un quadrato
    std::cout << "\n2. 34 cities inside square:" << std::endl;
    TSPGA tsp2(34, 100, false, false);
    tsp2.initialize();
    tsp2.saveCities("cities_square.dat");
    tsp2.run(1000);
    tsp2.saveBestPath("best_path_square.dat");
    tsp2.saveFitnessHistory("fitness_square.dat");
    
    std::cout << "\nExperiments completed! Check fitness_history.dat and best_path.dat" << std::endl;
    return 0;
}