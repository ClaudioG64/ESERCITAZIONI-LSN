#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"

using namespace std;

// Funzione d'onda psi_210 (in coordinate cartesiane)
double psi_210(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);
    if (r == 0) return 0;  // Evita problemi numerici
    double theta = acos(z / r);
    double a0 = 1.0;  // unità di Bohr (LJ units)

    double psi = (1.0 / (8.0 * sqrt(2.0 * M_PI))) * (r / a0) * exp(-r / (2.0 * a0)) * cos(theta);
    return psi * psi;  // Densità di probabilità (modulo quadro)
}

int main() {
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // Parametri Metropolis
    int M = 100000;  // Numero di step totali
    int N = 100;  // Numero di blocchi
    int L = M / N;  // Step per blocco
    double delta = 1.2;  // Passo massimo per la transizione uniforme

    // Posizione iniziale arbitraria
    double x = 2.0, y = 2.0, z = 2.0;
    double x_new, y_new, z_new;

    ofstream output("r_mean_psi210.dat");
    ofstream positions("positions_psi210.dat");  // 🔴 NUOVO FILE PER LE POSIZIONI

    double sum_r = 0.0, sum_r2 = 0.0;
    for (int i = 0; i < N; i++) {
        double sum_block = 0.0;

        for (int j = 0; j < L; j++) {
            x_new = x + rnd.Rannyu(-delta, delta);
            y_new = y + rnd.Rannyu(-delta, delta);
            z_new = z + rnd.Rannyu(-delta, delta);

            double p_old = psi_210(x, y, z);
            double p_new = psi_210(x_new, y_new, z_new);
            double alpha = min(1.0, p_new / p_old);

            if (rnd.Rannyu() < alpha) {  // Accetta la nuova posizione
                x = x_new;
                y = y_new;
                z = z_new;
                positions << x << " " << y << " " << z << endl;  // 🔴 Salva la posizione accettata
            }

            sum_block += sqrt(x*x + y*y + z*z);
        }

        double mean_r = sum_block / L;
        sum_r += mean_r;
        sum_r2 += mean_r * mean_r;

        double mean_prog = sum_r / (i + 1);
        double mean2_prog = sum_r2 / (i + 1);
        double sigma = sqrt((mean2_prog - mean_prog * mean_prog) / i);

        output << i + 1 << " " << mean_prog << " " << sigma << endl;
    }

    output.close();
    positions.close();  // 🔴 Chiude il file delle posizioni
    rnd.SaveSeed();
    return 0;
}
