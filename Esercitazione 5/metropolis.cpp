#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

// Impostazioni della simulazione
const int M = 1000000; // Numero di step Metropolis
const int N = 100;    // Numero di blocchi
const double delta = 1.2; // Ampiezza dello step uniforme

// Struttura per la posizione
struct Vec3 {
    double x, y, z;
};

// Funzione d'onda del ground state dell'idrogeno (1s)
double Psi_1s(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);
    return exp(-r) / sqrt(M_PI);
}

// Funzione per il sampling di Metropolis con passo uniforme
Vec3 MetropolisStep(Vec3 r, Random &rnd) {
    Vec3 r_new;
    r_new.x = r.x + rnd.Rannyu(-delta, delta);  //cambia Rannyu -> Gauss(0, sigma)
    r_new.y = r.y + rnd.Rannyu(-delta, delta);
    r_new.z = r.z + rnd.Rannyu(-delta, delta);

    double p_old = pow(Psi_1s(r.x, r.y, r.z), 2);
    double p_new = pow(Psi_1s(r_new.x, r_new.y, r_new.z), 2);

    if (p_new / p_old >= rnd.Rannyu()) {
        return r_new; // Accetta la nuova posizione
    } else {
        return r; // Rimani nella vecchia posizione
    }
}

// Funzione per calcolare l'errore statistico
double error(double av, double av2, int n) {
    if (n == 0) return 0;
    return sqrt((av2 - av * av) / n);
}

// Funzione principale
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

    // Output
    ofstream output("r_mean.dat");
    ofstream positions("positions.dat");  // Aggiunto per il salvataggio delle posizioni

    // Posizione iniziale
    Vec3 r = {1.0, 1.0, 1.0};

    double sum_r = 0.0, sum_r2 = 0.0;
    int accepted = 0;
    
    int L = M / N;  // Numero di step per blocco

    // Algoritmo di Metropolis
    for (int i = 0; i < N; i++) {
        double block_sum = 0;
        double block_sum2 = 0;

        for (int j = 0; j < L; j++) {
            Vec3 r_new = MetropolisStep(r, rnd);
            if (r_new.x != r.x || r_new.y != r.y || r_new.z != r.z) {
                accepted++;
            }
            r = r_new;

            // Salvataggio della posizione accettata nel file positions.dat
            positions << r.x << " " << r.y << " " << r.z << endl;

            double r_mod = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
            block_sum += r_mod;
            block_sum2 += r_mod * r_mod;
        }

        // Media per il blocco
        double mean = block_sum / L;
        double mean2 = block_sum2 / L;
        
        sum_r += mean;
        sum_r2 += mean2;

        double prog_mean = sum_r / (i + 1);
        double prog_mean2 = sum_r2 / (i + 1);
        double err = error(prog_mean, prog_mean2, i);

        output << i + 1 << " " << prog_mean << " " << err << endl;
    }

    cout << "Accettanza: " << (double)accepted / M << endl;
    output.close();
    positions.close();  // Chiudi il file delle posizioni
    rnd.SaveSeed();
    return 0;
}
