#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include "random.h"

using namespace std;

const int Nsteps = 100;   // Number of steps
const int M = 10000;      // Number of random walks (trials)
const int Nblocks = 100;  // Number of blocks for statistical uncertainty
const double a = 1.0;     // Step size

// Function to simulate the lattice random walk
void RandomWalkLattice(Random &rnd, vector<double> &r_lattice) {
    int x = 0, y = 0, z = 0;  // Initial position (origin)
    for (int i = 0; i < Nsteps; ++i) {
        int direction = rnd.Rannyu(0, 3);  // Choose a direction (0=x, 1=y, 2=z)
        int step = rnd.Rannyu(0, 2) * 2 - 1;  // Choose forward (+1) or backward (-1)
        
        if (direction == 0) x += step;  // x direction
        else if (direction == 1) y += step;  // y direction
        else z += step;  // z direction
        
        r_lattice[i] += sqrt(x*x + y*y + z*z);  // Update the distance
    }
}

// Function to simulate the continuum random walk
void RandomWalkContinuum(Random &rnd, vector<double> &r_continuum) {
    double theta, phi;
    double x = 0, y = 0, z = 0;  // Initial position (origin)
    
    for (int i = 0; i < Nsteps; ++i) {
        theta = rnd.Rannyu(0, M_PI);  // Uniformly distributed theta
        phi = rnd.Rannyu(0, 2*M_PI);  // Uniformly distributed phi
        
        // Calculate the step components in the x, y, and z directions
        double step_x = a * sin(theta) * cos(phi);
        double step_y = a * sin(theta) * sin(phi);
        double step_z = a * cos(theta);
        
        // Update position
        x += step_x;
        y += step_y;
        z += step_z;
        
        r_continuum[i] += sqrt(x*x + y*y + z*z);  // Update the distance
    }
}

// Function to calculate the statistical uncertainty
void ComputeUncertainty(vector<vector<double>> &r_lattice_block, vector<vector<double>> &r_continuum_block, vector<double> &avg_r_lattice, vector<double> &avg_r_continuum) {
    for (int i = 0; i < Nsteps; ++i) {
        // Calculate average and uncertainty for lattice walk
        double mean_lattice = 0.0;
        double mean_continuum = 0.0;
        for (int j = 0; j < Nblocks; ++j) {
            mean_lattice += r_lattice_block[j][i];
            mean_continuum += r_continuum_block[j][i];
        }
        mean_lattice /= Nblocks;
        mean_continuum /= Nblocks;
        
        avg_r_lattice[i] = mean_lattice;
        avg_r_continuum[i] = mean_continuum;
    }
}

// Main function
int main() {
    Random rnd;
    
    // Set the random seed
    int seed[4] = {1, 2, 3, 4};
    rnd.SetRandom(seed, 0, 0);

    // Vectors to store distances for each random walk
    vector<double> r_lattice(Nsteps, 0.0);
    vector<double> r_continuum(Nsteps, 0.0);
    
    // Vectors to store distances in blocks for uncertainty computation
    vector<vector<double>> r_lattice_block(Nblocks, vector<double>(Nsteps, 0.0));
    vector<vector<double>> r_continuum_block(Nblocks, vector<double>(Nsteps, 0.0));

    // Perform M trials and accumulate the distances in blocks
    for (int j = 0; j < M; ++j) {
        RandomWalkLattice(rnd, r_lattice);   // Perform lattice walk
        RandomWalkContinuum(rnd, r_continuum); // Perform continuum walk
        
        // Store results in blocks for statistical uncertainty calculation
        for (int i = 0; i < Nsteps; ++i) {
            r_lattice_block[j % Nblocks][i] = r_lattice[i];
            r_continuum_block[j % Nblocks][i] = r_continuum[i];
        }
    }

    // Average the distances and compute the uncertainty
    vector<double> avg_r_lattice(Nsteps, 0.0);
    vector<double> avg_r_continuum(Nsteps, 0.0);
    ComputeUncertainty(r_lattice_block, r_continuum_block, avg_r_lattice, avg_r_continuum);

    // Write the results to a file for plotting
    ofstream outFile("random_walk_results.dat");
    for (int i = 0; i < Nsteps; ++i) {
        outFile << i << " " << avg_r_lattice[i] << " " << avg_r_continuum[i] << endl;
    }
    outFile.close();

    return 0;
}
