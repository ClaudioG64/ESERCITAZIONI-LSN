#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

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
            if(property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // Parametri
    int N = 100;      // numero di sottointervalli
    int M = 10000;    // numeri estratti per ogni test
    int ntests = 100; // numero di test indipendenti

    double expected = (double)M / N;

    ofstream output("output1_3.dat");

    for(int t=0; t<ntests; t++){
        vector<int> counts(N, 0);

        // Estrazioni
        for(int i=0; i<M; i++){
            double r = rnd.Rannyu();
            int bin = int(r * N); // assegno il numero casuale a un intervallo
            counts[bin]++;
        }

        // Calcolo chi-quadro
        double chi2 = 0.0;
        for(int i=0; i<N; i++){
            chi2 += pow(counts[i] - expected, 2) / expected;
        }

        output << t+1 << " " << chi2 << endl;
    }

    output.close();
    rnd.SaveSeed();

    return 0;
}
