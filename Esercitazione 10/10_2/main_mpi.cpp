#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <limits>
#include "TSPGA.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int N_migr = 10; // Default: migrazione ogni 10 generazioni
    if (argc > 1) {
        N_migr = atoi(argv[1]);
    }

    const int generations = 5000;
    const int populationSize = 200;
    const bool useL1 = false;

    srand(time(0) + world_rank);

    // Lettura delle città dal file (solo rank 0)
    std::vector<City> cities;
    int n_cities = 0;
    if (world_rank == 0) {
        std::ifstream file("cap_prov_ita.dat");
        if (file.is_open()) {
            double x, y;
            while (file >> x >> y) {
                cities.emplace_back(x, y);
            }
            file.close();
            n_cities = cities.size();
            std::cout << "Loaded " << n_cities << " cities from file." << std::endl;

            // Normalizza le coordinate
            double min_x = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::lowest();
            double min_y = std::numeric_limits<double>::max();
            double max_y = std::numeric_limits<double>::lowest();
            
            for (const auto& city : cities) {
                min_x = std::min(min_x, city.x);
                max_x = std::max(max_x, city.x);
                min_y = std::min(min_y, city.y);
                max_y = std::max(max_y, city.y);
            }
            
            for (auto& city : cities) {
                city.x = (city.x - min_x) / (max_x - min_x);
                city.y = (city.y - min_y) / (max_y - min_y);
            }
        } else {
            std::cerr << "Unable to open file cap_prov_ita.dat" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&n_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double* x_coords = new double[n_cities];
    double* y_coords = new double[n_cities];
    
    if (world_rank == 0) {
        for (int i = 0; i < n_cities; i++) {
            x_coords[i] = cities[i].x;
            y_coords[i] = cities[i].y;
        }
    }

    MPI_Bcast(x_coords, n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y_coords, n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

    // File per salvare la storia della fitness (solo rank 0)
    std::ofstream fitnessFile;
    if (world_rank == 0) {
        fitnessFile.open("fitness_history_ita.dat");
    }

    // Loop principale
    for (int gen = 0; gen < generations; gen++) {
        tsp.nextGeneration();

        if (gen % N_migr == 0 && world_size > 1) {
            // Migrazione
            Individual best = tsp.getBestIndividual();
            int* best_path = new int[n_cities];
            for (int i = 0; i < n_cities; i++) {
                best_path[i] = best.path[i];
            }

            int* all_bests = new int[world_size * n_cities];
            MPI_Allgather(best_path, n_cities, MPI_INT, all_bests, n_cities, MPI_INT, MPI_COMM_WORLD);

            int index;
            do {
                index = rand() % world_size;
            } while (index == world_rank);

            Individual received;
            received.path.resize(n_cities);
            for (int i = 0; i < n_cities; i++) {
                received.path[i] = all_bests[index * n_cities + i];
            }

            if (tsp.checkIndividual(received)) {
                tsp.replaceWorstIndividual(received);
            }

            delete[] best_path;
            delete[] all_bests;
        }

        // Raccolta statistiche fitness (solo rank 0)
        if (world_rank == 0) {
            double best_fitness = tsp.getBestIndividual().fitness;
            double avg_fitness = tsp.getAverageFitness();

            fitnessFile << gen << " " << best_fitness << " " << avg_fitness << std::endl;
        }

        if (world_rank == 0 && gen % 100 == 0) {
            std::cout << "Generation " << gen << " completed" << std::endl;
        }
    }

    if (world_rank == 0) {
        fitnessFile.close();
    }

    // Raccolta risultato finale
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
        
        Individual global_best;
        global_best.path.resize(n_cities);
        for (int i = 0; i < n_cities; i++) {
            global_best.path[i] = all_paths[best_index * n_cities + i];
        }
        
        TSPGA temp_tsp(populationSize, useL1, 0);
        temp_tsp.setCities(cities);
        temp_tsp.initialize();
        temp_tsp.replaceWorstIndividual(global_best);
        temp_tsp.saveBestPath("global_best_path_ita.dat");
        temp_tsp.saveCities("cities_ita.dat");
        
        delete[] all_fitnesses;
        delete[] all_paths;
    }

    delete[] local_best_path;
    MPI_Finalize();
    return 0;
}