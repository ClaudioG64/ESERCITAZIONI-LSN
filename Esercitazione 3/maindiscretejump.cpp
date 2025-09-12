//questo codice per rispondeere alla richiesta 1 (time step unico discreto)

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

// Funzione per calcolare il payoff di un'opzione call
double call_payoff(double S, double K) {
    return max(0.0, S - K);
}

// Funzione per calcolare il payoff di un'opzione put
double put_payoff(double S, double K) {
    return max(0.0, K - S);
}

// Funzione per simulare il prezzo finale dell'asset
double simulate_asset_price(double S0, double mu, double sigma, double T, Random &rnd) {
    double gauss_bm = rnd.Gauss(0, 1);
    return S0 * exp((mu - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * gauss_bm);
}

int main() {
    // Parametri dell'opzione e del mercato
    double S0 = 100.0;   // Prezzo iniziale dell'asset
    double K = 100.0;    // Prezzo di esercizio
    double T = 1.0;      // Tempo alla scadenza in anni
    double r = 0.1;     // Tasso di interesse privo di rischio
    double sigma = 0.25;  // Volatilità dell'asset

    int M = 100000;      // Numero di simulazioni
    int N = 100;         // Numero di blocchi
    int L = M / N;       // Numero di simulazioni per blocco

    // Inizializzazione del generatore di numeri casuali
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
        return 1;
    }
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
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        return 1;
    }

    // Vettori per memorizzare i risultati
    vector<double> call_prices(N, 0.0);
    vector<double> put_prices(N, 0.0);

    // Simulazione Monte Carlo
    for (int i = 0; i < N; ++i) {
        double call_sum = 0.0;
        double put_sum = 0.0;
        for (int j = 0; j < L; ++j) {
            double S_T = simulate_asset_price(S0, r, sigma, T, rnd);
            call_sum += exp(-r * T) * call_payoff(S_T, K);
            put_sum += exp(-r * T) * put_payoff(S_T, K);
        }
        call_prices[i] = call_sum / L;
        put_prices[i] = put_sum / L;
    }

    // Calcolo delle medie progressive e degli errori
    vector<double> call_avg(N, 0.0), call_err(N, 0.0);
    vector<double> put_avg(N, 0.0), put_err(N, 0.0);

    for (int i = 0; i < N; ++i) {
        double call_sum_prog = 0.0, call_sum2_prog = 0.0;
        double put_sum_prog = 0.0, put_sum2_prog = 0.0;
        for (int j = 0; j <= i; ++j) {
            call_sum_prog += call_prices[j];
            call_sum2_prog += call_prices[j] * call_prices[j];
            put_sum_prog += put_prices[j];
            put_sum2_prog += put_prices[j] * put_prices[j];
        }
        call_avg[i] = call_sum_prog / (i + 1);
        put_avg[i] = put_sum_prog / (i + 1);
        if (i > 0) {
            call_err[i] = sqrt((call_sum2_prog / (i + 1) - call_avg[i] * call_avg[i]) / i);
            put_err[i] = sqrt((put_sum2_prog / (i + 1) - put_avg[i] * put_avg[i]) / i);
        }
    }

    // Salvataggio dei risultati su file
    ofstream output("option_prices.dat");
    if (output.is_open()) {
        output << "# Block  Call_Price  Call_Error  Put_Price  Put_Error" << endl;
        for (int i = 0; i < N; ++i) {
            output << i + 1 << " "
                   << call_avg[i] << " " << call_err[i] << " "
                   << put_avg[i] << " " << put_err[i] << endl;
        }
        output.close();
    } else {
        cerr << "PROBLEM: Unable to open option_prices.dat" << endl;
        return 1;
    }

    rnd.SaveSeed();
    return 0;
}
