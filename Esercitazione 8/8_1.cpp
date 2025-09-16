#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <limits>
#include "random.h"

using namespace std;

// Potential
double V(double x) { return pow(x,4) - 2.5*pow(x,2); }

// log g(x;a,sigma) = - (x-a)^2 / (2 sigma^2)
inline double log_g(double x, double a, double sigma) {
    double dx = x - a;
    return - (dx*dx) / (2.0 * sigma * sigma);
}

// log psi using log-sum-exp for numerical stability
double log_psi(double x, double sigma, double mu) {
    double l1 = log_g(x, +mu, sigma);
    double l2 = log_g(x, -mu, sigma);
    double m = std::max(l1, l2);
    double w = exp(l1 - m) + exp(l2 - m);
    return m + log(w);
}

// acceptance ratio computed safely from logpsi difference
inline double acceptance_ratio_from_logs(double logpsi_new, double logpsi_old) {
    double log_alpha = 2.0 * (logpsi_new - logpsi_old); // because psi^2
    if (log_alpha >= 0.0) return 1.0;
    return exp(log_alpha);
}

// stable local energy: compute psi''/psi without huge exponentials
double local_energy_stable(double x, double mu, double sigma) {
    // compute log-weights and scaled weights
    double l1 = log_g(x, mu, sigma);   // for g(x,a=mu)
    double l2 = log_g(x, -mu, sigma);  // for g(x,a=-mu)
    double m = std::max(l1, l2);
    double w1 = exp(l1 - m);
    double w2 = exp(l2 - m);
    double denom = w1 + w2;

    // factors for g''/g : (x-a)^2/sigma^4 - 1/sigma^2
    double xp = x - mu;
    double xm = x + mu; // since a=-mu -> x - (-mu) = x + mu

    double factor1 = (xp*xp) / (pow(sigma,4)) - 1.0/(sigma*sigma);
    double factor2 = (xm*xm) / (pow(sigma,4)) - 1.0/(sigma*sigma);

    // psi''/psi = (factor1*g1 + factor2*g2) / (g1 + g2)
    // but scaled: g1 ~ w1*exp(m), g2 ~ w2*exp(m) -> exp(m) cancels
    double psi_dd_over_psi = (factor1 * w1 + factor2 * w2) / denom;

    double kinetic = -0.5 * psi_dd_over_psi;
    return kinetic + V(x);
}

// block-averaging error from block means
void block_stats_from_samples(const vector<double> &samples, int block_size, double &mean, double &err, int &nblocks) {
    int N = samples.size();
    if (N <= 0) { mean = 0.0; err = 0.0; nblocks = 0; return; }
    nblocks = N / block_size;
    if (nblocks <= 0) { // not enough samples for one block
        // fallback: simple mean and naive error
        double s=0, s2=0;
        for (double v : samples) { s+=v; s2+=v*v; }
        mean = s/N;
        double var = s2/N - mean*mean;
        if (var < 0) var = 0;
        err = sqrt(var / N);
        return;
    }
    vector<double> blocks(nblocks, 0.0);
    for (int b=0; b<nblocks; ++b) {
        double s=0;
        for (int i=0; i<block_size; ++i) s += samples[b*block_size + i];
        blocks[b] = s / block_size;
    }
    // mean over blocks
    double sb=0;
    for (double v: blocks) sb += v;
    mean = sb / nblocks;
    // sample standard deviation of block means
    double sd2 = 0.0;
    for (double v: blocks) sd2 += (v - mean)*(v - mean);
    if (nblocks > 1) {
        double sd = sqrt(sd2 / (nblocks - 1));
        err = sd / sqrt(double(nblocks));
    } else {
        err = 0.0;
    }
}

int main() {
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (!Primes.is_open()) { cerr<<"Unable to open Primes"<<endl; return 1; }
    Primes >> p1 >> p2;
    Primes.close();
    ifstream input("seed.in");
    if (!input.is_open()) { cerr<<"Unable to open seed.in"<<endl; return 1; }
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    input.close();
    rnd.SetRandom(seed, p1, p2);

    ofstream out("energies_scan.dat");
    out << "# mu sigma energy error acceptance" << endl;

    // grid params (tweak if you want)
    double mu_min = 0.5, mu_max = 1.0, dmu = 0.02;
    double sigma_min = 0.4, sigma_max = 0.8, dsigma = 0.02;

    int M_total = 20000;    // total Metropolis steps per parameter set
    int burn_in = 2000;
    double delta = 1.5;     // step size for uniform proposal; tune to ~0.5 acceptance
    int block_size = 100;   // block size for blocking analysis

    double bestE = numeric_limits<double>::infinity();
    double bestErr = 0;
    double bestMu = 0, bestSigma = 0;

    for (double mu = mu_min; mu <= mu_max + 1e-9; mu += dmu) {
        for (double sigma = sigma_min; sigma <= sigma_max + 1e-9; sigma += dsigma) {

            double x = 0.0;
            double logpsi_old = log_psi(x, sigma, mu);
            int accepted = 0;

            vector<double> samples;
            samples.reserve(M_total - burn_in);

            for (int step=0; step < M_total; ++step) {
                // propose
                double x_new = x + rnd.Rannyu(-delta, delta);
                double logpsi_new = log_psi(x_new, sigma, mu);

                double alpha = acceptance_ratio_from_logs(logpsi_new, logpsi_old);
                if (alpha >= 1.0 || rnd.Rannyu() < alpha) {
                    x = x_new;
                    logpsi_old = logpsi_new;
                    ++accepted;
                }
                if (step >= burn_in) {
                    double e = local_energy_stable(x, mu, sigma);
                    samples.push_back(e);
                }
            }

            int N = samples.size();
            double meanE, errE;
            int nblocks;
            block_stats_from_samples(samples, block_size, meanE, errE, nblocks);

            double acc_rate = double(accepted) / double(M_total);
            out << mu << " " << sigma << " " << setprecision(10) << meanE << " " << errE << " " << acc_rate << endl;
            cout << "mu="<<mu<<" sigma="<<sigma<<" E="<<meanE<<" +/- "<<errE<<" acc="<<acc_rate<< " blocks="<<nblocks<<endl;

            if (meanE < bestE) {
                bestE = meanE;
                bestErr = errE;
                bestMu = mu;
                bestSigma = sigma;
            }
        }
    }

    cout << "\n=== BEST PARAMETERS FOUND ===\n";
    cout << "mu = " << bestMu << "  sigma = " << bestSigma << endl;
    cout << "E_min = " << bestE << " +/- " << bestErr << endl;

    out << "# === BEST ===" << endl;
    out << "# mu_best sigma_best E_min E_err" << endl;
    out << bestMu << " " << bestSigma << " " << bestE << " " << bestErr << endl;

    out.close();
    rnd.SaveSeed();
    return 0;
}
