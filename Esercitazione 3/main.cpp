#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

// payoff
inline double call_payoff(double S, double K) { return (S > K) ? (S - K) : 0.0; }
inline double put_payoff(double S, double K)  { return (S < K) ? (K - S) : 0.0; }

// errore statistico (data blocking)
inline double error(double mean, double mean2, int n) {
    if (n == 0) return 0.0;
    return sqrt((mean2 - mean*mean)/n);
}

// simulazione di un path GBM discretizzato con  'steps'  sottointervalli
vector<double> simulate_asset_path(double S0, double mu, double sigma, double T, int steps, Random &rnd) {
    vector<double> path(steps + 1);
    path[0] = S0;
    double dt = T / steps;
    for (int i = 1; i <= steps; ++i) {
        double z = rnd.Gauss(0.0, 1.0);
        path[i] = path[i-1] * exp((mu - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*z);
    }
    return path;
}

int main() {
    // Parametri opzione/mercato
    const double S0 = 100.0;
    const double K  = 100.0;
    const double T  = 1.0;
    const double r  = 0.1;
    const double sigma = 0.25;

    // Monte Carlo + Data Blocking
    const int M = 100000;  // lanci totali
    const int N = 100;     // blocchi
    const int L = M / N;   // lanci per blocco
    const int steps = 100; // discretizzazione richiesta

    // RNG (come da file del prof)
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) { Primes >> p1 >> p2; }
    else { cerr << "PROBLEM: Unable to open Primes\n"; return 1; }
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while (!input.eof()){
            input >> property;
            if (property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else { cerr << "PROBLEM: Unable to open seed.in\n"; return 1; }

    // --- vettori per i blocchi (DIRECT) ---
    vector<double> call_direct_block(N, 0.0), put_direct_block(N, 0.0);

    // --- vettori per i blocchi (DISCRETIZED PATH) ---
    vector<double> call_disc_block(N, 0.0),  put_disc_block(N, 0.0);

    // --------- MONTE CARLO ----------
    for (int i = 0; i < N; ++i) {
        double sum_call_direct = 0.0, sum_put_direct = 0.0;
        double sum_call_disc   = 0.0, sum_put_disc   = 0.0;

        for (int j = 0; j < L; ++j) {
            // 1) DIRECT SAMPLING OF S(T)
            // GBM: S_T = S0 * exp((r - 0.5*sigma^2) T + sigma sqrt(T) Z)
            double z = rnd.Gauss(0.0, 1.0);
            double S_T_direct = S0 * exp((r - 0.5*sigma*sigma)*T + sigma*sqrt(T)*z);

            sum_call_direct += exp(-r*T) * call_payoff(S_T_direct, K);
            sum_put_direct  += exp(-r*T) * put_payoff(S_T_direct, K);

            // 2) DISCRETIZED GBM PATH (100 step)
            auto path = simulate_asset_path(S0, r, sigma, T, steps, rnd);
            double S_T_disc = path.back();

            sum_call_disc += exp(-r*T) * call_payoff(S_T_disc, K);
            sum_put_disc  += exp(-r*T) * put_payoff(S_T_disc, K);
        }

        call_direct_block[i] = sum_call_direct / L;
        put_direct_block[i]  = sum_put_direct  / L;

        call_disc_block[i]   = sum_call_disc   / L;
        put_disc_block[i]    = sum_put_disc    / L;
    }

    // --------- DATA BLOCKING: medie progressive + errori ---------
    ofstream out_direct("option_prices_direct.dat");
    ofstream out_disc("option_prices_discretized.dat");

    if (!out_direct || !out_disc) {
        cerr << "PROBLEM: Unable to open output files\n";
        return 1;
    }

    out_direct << "# Block  Call_Price  Call_Error  Put_Price  Put_Error\n";
    out_disc   << "# Block  Call_Price  Call_Error  Put_Price  Put_Error\n";

    double sum_c_d=0, sum2_c_d=0, sum_p_d=0, sum2_p_d=0; // direct
    double sum_c_s=0, sum2_c_s=0, sum_p_s=0, sum2_p_s=0; // discretized

    for (int i = 0; i < N; ++i) {
        // DIRECT
        sum_c_d  += call_direct_block[i];
        sum2_c_d += call_direct_block[i]*call_direct_block[i];
        sum_p_d  += put_direct_block[i];
        sum2_p_d += put_direct_block[i]*put_direct_block[i];

        double mean_c_d  = sum_c_d/(i+1);
        double mean2_c_d = sum2_c_d/(i+1);
        double err_c_d   = error(mean_c_d, mean2_c_d, i);

        double mean_p_d  = sum_p_d/(i+1);
        double mean2_p_d = sum2_p_d/(i+1);
        double err_p_d   = error(mean_p_d, mean2_p_d, i);

        out_direct << (i+1) << " "
                   << mean_c_d << " " << err_c_d << " "
                   << mean_p_d << " " << err_p_d << "\n";

        // DISCRETIZED
        sum_c_s  += call_disc_block[i];
        sum2_c_s += call_disc_block[i]*call_disc_block[i];
        sum_p_s  += put_disc_block[i];
        sum2_p_s += put_disc_block[i]*put_disc_block[i];

        double mean_c_s  = sum_c_s/(i+1);
        double mean2_c_s = sum2_c_s/(i+1);
        double err_c_s   = error(mean_c_s, mean2_c_s, i);

        double mean_p_s  = sum_p_s/(i+1);
        double mean2_p_s = sum2_p_s/(i+1);
        double err_p_s   = error(mean_p_s, mean2_p_s, i);

        out_disc << (i+1) << " "
                 << mean_c_s << " " << err_c_s << " "
                 << mean_p_s << " " << err_p_s << "\n";
    }

    out_direct.close();
    out_disc.close();

    rnd.SaveSeed();
    return 0;
}
