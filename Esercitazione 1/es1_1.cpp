#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

double error(double av, double av2, int n) {
    if (n == 0) return 0;
    return sqrt((av2 - av*av)/n);
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
            if(property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // Parametri
    int M = 100000;   // Numero totale di numeri casuali
    int N = 100;      // Numero di blocchi
    int L = M/N;      // Numeri per blocco

    // Uso vector al posto degli array statici
    vector<double> ave(N), av2(N);

    for(int i=0; i<N; i++){
        double sum = 0.0;
        for(int j=0; j<L; j++){
            sum += rnd.Rannyu();
        }
        ave[i] = sum/L;      // media del blocco
        av2[i] = ave[i]*ave[i];
    }

    // Medie progressive
    ofstream output("output1_1.dat");
    double sum_prog = 0.0, su2_prog = 0.0;
    for(int i=0; i<N; i++){
        sum_prog += ave[i];
        su2_prog += av2[i];
        double mean = sum_prog/(i+1);
        double mean2 = su2_prog/(i+1);
        double err = error(mean, mean2, i);
        output << i+1 << " " << mean << " " << err << endl;
    }
    output.close();

    rnd.SaveSeed();
    return 0;
}
