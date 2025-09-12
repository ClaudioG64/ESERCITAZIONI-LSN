#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

// Funzione da integrare
double f(double x) {
    return (M_PI / 2) * cos(M_PI * x / 2);
}

// PDF Importance Sampling: g(x) = 2(1-x)
double g(double x) {
    return 2 * (1 - x);
}

// Inversa della cumulativa di g(x)
double inverse_g(Random &rnd) {
    double r = rnd.Rannyu();
    return 1 - sqrt(1 - r);
}

// Calcolo della media progressiva e dell'errore
void progressiveMean(const vector<double> &data, vector<double> &mean, vector<double> &err) {
    int N = data.size();
    mean.resize(N, 0.0);
    err.resize(N, 0.0);
    
    vector<double> sum_prog(N, 0.0), sum2_prog(N, 0.0);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            sum_prog[i] += data[j];
            sum2_prog[i] += data[j] * data[j];
        }
        sum_prog[i] /= (i + 1);
        sum2_prog[i] /= (i + 1);
        mean[i] = sum_prog[i];
        if (i > 0) {
            err[i] = sqrt((sum2_prog[i] - sum_prog[i] * sum_prog[i]) / i);
        }
    }
}

int main() {

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while (!input.eof()){
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    int M = 100000; // Numero totale di punti
    int N = 100;    // Numero di blocchi
    int L = M / N;  // Punti per blocco
    
    vector<double> results_uniform, results_importance;
    
    for (int i = 0; i < N; i++) {
        double sum_uniform = 0.0;
        double sum_importance = 0.0;
        
        for (int j = 0; j < L; j++) {
            double x_uniform = rnd.Rannyu();
            sum_uniform += f(x_uniform);
            
            double x_importance = inverse_g(rnd);
            sum_importance += f(x_importance) / g(x_importance);
        }
        
        results_uniform.push_back(sum_uniform / L);
        results_importance.push_back(sum_importance / L);
    }
    
    vector<double> mean_uniform, err_uniform;
    vector<double> mean_importance, err_importance;
    
    progressiveMean(results_uniform, mean_uniform, err_uniform);
    progressiveMean(results_importance, mean_importance, err_importance);
    
    ofstream output("output.dat");
    if (!output) {
        cerr << "PROBLEM: Unable to open output.dat" << endl;
        return 1;
    }
    
    // Scrive intestazione nel file
    output << "# N  mean_uniform  error_uniform  mean_importance  error_importance" << endl;
    
    for (int i = 0; i < N; i++) {
        output << i + 1 << " "
               << mean_uniform[i] << " " << err_uniform[i] << " "
               << mean_importance[i] << " " << err_importance[i] << endl;
    }
    
    output.close();
    cout << "Dati salvati in output.dat" << endl;

    rnd.SaveSeed();
    
    return 0;
}
