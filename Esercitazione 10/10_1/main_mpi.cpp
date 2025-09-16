#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "TSPGA.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int N_migr = 10; // Migrazione ogni 10 generazioni
    const int generations = 1000;
    const int populationSize = 100;
    const bool useL1 = false; // Usa L2 distance
    const int numCities = 34; // Per esercizio 10.1
    const bool citiesOnCircle = false;  //False => città in un quadrato, true => città sulla circonferenza

    // Inizializza il random seed per ogni processo in modo diverso
    srand(time(0) + world_rank);

    // Genera le città casualmente (solo rank 0)
    std::vector<City> cities;
    if (world_rank == 0) {
        std::cout << "Generating " << numCities << " cities on " 
                  << (citiesOnCircle ? "circumference" : "square") << std::endl;
        
        if (citiesOnCircle) {
            // Città sulla circonferenza
            for (int i = 0; i < numCities; ++i) {
                double angle = 2.0 * M_PI * rand() / RAND_MAX;
                cities.emplace_back(std::cos(angle), std::sin(angle));
            }
        } else {
            // Città dentro un quadrato [0,1]x[0,1]
            for (int i = 0; i < numCities; ++i) {
                double x = static_cast<double>(rand()) / RAND_MAX;
                double y = static_cast<double>(rand()) / RAND_MAX;
                cities.emplace_back(x, y);
            }
        }
    }

    // Broadcast di numCities a tutti i processi
    int n_cities = numCities;
    MPI_Bcast(&n_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Preparazione degli array per le coordinate
    double* x_coords = new double[n_cities];
    double* y_coords = new double[n_cities];
    
    if (world_rank == 0) {
        for (int i = 0; i < n_cities; i++) {
            x_coords[i] = cities[i].x;
            y_coords[i] = cities[i].y;
        }
    }

    // Broadcast delle coordinate
    MPI_Bcast(x_coords, n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y_coords, n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Ricostruzione del vettore cities negli altri processi
    if (world_rank != 0) {
        for (int i = 0; i < n_cities; i++) {
            cities.emplace_back(x_coords[i], y_coords[i]);
        }
    }
    delete[] x_coords;
    delete[] y_coords;

    // Inizializzazione del TSPGA
    TSPGA tsp(populationSize, useL1, world_rank);
    tsp.setCities(cities);
    tsp.initialize();

    // Loop principale delle generazioni
    for (int gen = 0; gen < generations; gen++) {
        tsp.nextGeneration();

        if (gen % N_migr == 0 && world_size > 1) {
            // Migrazione: scambio dei migliori individui
            Individual best = tsp.getBestIndividual();

            // Assicurati che il best individual sia valido
            if (!tsp.checkIndividual(best)) {
                if (world_rank == 0) {
                    std::cerr << "Error: Best individual is invalid!" << std::endl;
                }
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            int* best_path = new int[n_cities];
            for (int i = 0; i < n_cities; i++) {
                best_path[i] = best.path[i];
            }

            int* all_bests = new int[world_size * n_cities];
            MPI_Allgather(best_path, n_cities, MPI_INT, all_bests, n_cities, MPI_INT, MPI_COMM_WORLD);

            // Scelta casuale di un individuo da un altro processo
            int index;
            do {
                index = rand() % world_size;
            } while (index == world_rank);

            Individual received;
            received.path.resize(n_cities);
            for (int i = 0; i < n_cities; i++) {
                received.path[i] = all_bests[index * n_cities + i];
            }

            // Controlla se l'individuo ricevuto è valido prima di sostituire
            if (tsp.checkIndividual(received)) {
                tsp.replaceWorstIndividual(received);
            } else {
                if (world_rank == 0) {
                    std::cerr << "Warning: Received invalid individual during migration from process " << index << std::endl;
                }
            }

            delete[] best_path;
            delete[] all_bests;
        }

        if (world_rank == 0 && gen % 100 == 0) {
            std::cout << "Generation " << gen << " completed" << std::endl;
        }
    }

    // Raccolta del miglior individuo globale
    Individual local_best = tsp.getBestIndividual();
    double local_best_fitness = local_best.fitness;
    int* local_best_path = new int[n_cities];
    for (int i = 0; i < n_cities; i++) {
        local_best_path[i] = local_best.path[i];
    }

    double* all_fitnesses = nullptr;
    int* all_paths = nullptr;
    if (world_rank == 0) {
        all_fitnesses = new double[world_size];
        all_paths = new int[world_size * n_cities];
    }

    MPI_Gather(&local_best_fitness, 1, MPI_DOUBLE, all_fitnesses, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(local_best_path, n_cities, MPI_INT, all_paths, n_cities, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        int best_index = 0;
        for (int i = 1; i < world_size; i++) {
            if (all_fitnesses[i] < all_fitnesses[best_index]) {
                best_index = i;
            }
        }
        
        std::cout << "Global best fitness: " << all_fitnesses[best_index] << std::endl;
        
        // Salva il miglior percorso globale
        Individual global_best;
        global_best.path.resize(n_cities);
        for (int i = 0; i < n_cities; i++) {
            global_best.path[i] = all_paths[best_index * n_cities + i];
        }
        
        // Usiamo una istanza temporanea di TSPGA per salvare
        TSPGA temp_tsp(populationSize, useL1, 0);
        temp_tsp.setCities(cities);
        // Creiamo una popolazione temporanea con solo il miglior individuo
        temp_tsp.initialize(); // Inizializza una popolazione
        temp_tsp.replaceWorstIndividual(global_best); // Sostituisce il peggiore con il globale
        temp_tsp.saveBestPath("global_best_path.dat");
        temp_tsp.saveCities("cities.dat");
        
        delete[] all_fitnesses;
        delete[] all_paths;
    }

    delete[] local_best_path;
    MPI_Finalize();
    return 0;
}