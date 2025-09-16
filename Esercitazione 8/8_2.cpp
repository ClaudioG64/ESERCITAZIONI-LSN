#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

/* ---------------------------
   FUNCTIONS FROM YOUR FILE
   --------------------------- */

double V(double x) {
    return pow(x,4) - 2.5*pow(x,2);
}

double g(double x, double a, double sigma) {
    return exp(-pow(x-a,2) / (2.0*sigma*sigma));
}

double psi(double x, double sigma, double mu) {
    return g(x,+mu,sigma) + g(x,-mu,sigma);
}

double psi2(double x, double sigma, double mu) {
    double p = psi(x,sigma,mu);
    return p*p;
}

// seconda derivata di g/g
double g_dd_over_g(double x, double a, double sigma) {
    double z = (x-a)/sigma;
    return (z*z)/(sigma*sigma) - 1.0/(sigma*sigma);
}

// stable local energy implementation
double local_energy_stable(double x, double mu, double sigma) {
    // compute log weights to avoid overflow
    double l1 = - ( (x - mu)*(x - mu) ) / (2.0 * sigma * sigma);
    double l2 = - ( (x + mu)*(x + mu) ) / (2.0 * sigma * sigma);
    double m = std::max(l1, l2);
    double w1 = exp(l1 - m); // scaled g1
    double w2 = exp(l2 - m); // scaled g2
    double denom = w1 + w2;  // proportional to psi (scaled)

    // factors: g''/g = (x-a)^2/sigma^4 - 1/sigma^2
    double factor1 = ( (x - mu)*(x - mu) ) / (pow(sigma,4)) - 1.0/(sigma*sigma);
    double factor2 = ( (x + mu)*(x + mu) ) / (pow(sigma,4)) - 1.0/(sigma*sigma);

    double psi_dd_over_psi = ( factor1 * w1 + factor2 * w2 ) / denom;
    double kinetic = -0.5 * psi_dd_over_psi;
    return kinetic + V(x);
}

// --- wrapper compatibile con il resto del codice ---
// mantiene il nome `local_energy` usato altrove ma usa l'implementazione stabile
double local_energy(double x, double mu, double sigma) {
    return local_energy_stable(x, mu, sigma);
}


/* ---------------------------
   HELPER: errore statistico
   --------------------------- */
double block_error(const vector<double>& v){
    int n = v.size();
    if(n <= 1) return 0.0;
    double sum=0, sum2=0;
    for(double x : v){ sum += x; sum2 += x*x; }
    double mean = sum/double(n);
    double var = (sum2/double(n) - mean*mean);
    if (var < 0) var = 0;
    return sqrt(var / double(n));
}

/* ---------------------------
   METROPOLIS SAMPLER + BLOCKING
   -> sample from |psi|^2 (Metropolis on x),
      compute local_energy(x) at each sampled x after burn-in
   Returns: mean energy, standard error, and the list of samples if requested
   --------------------------- */
struct EnergyResult {
    double mean;
    double error;
    double acceptance_rate;
    vector<double> samples; // energies (per measurement) or x samples if requested
};

EnergyResult estimate_energy(Random &rnd,
                             double mu, double sigma,
                             int M_total, int burn_in,
                             double delta_move,
                             int n_blocks,
                             bool save_x_samples=false,
                             vector<double> *out_x_samples = nullptr)
{
    // set up
    int M_eff = M_total - burn_in;
    if(M_eff <= 0) {
        cerr << "[estimate_energy] M_total <= burn_in !\n";
        exit(EXIT_FAILURE);
    }
    if(n_blocks <= 0) n_blocks = 1;
    int block_size = M_eff / n_blocks;
    if(block_size <= 0) block_size = M_eff; // fallback

    vector<double> block_means;
    block_means.reserve(n_blocks);

    vector<double> energies_per_measure; energies_per_measure.reserve(M_eff);
    vector<double> x_samples; if(save_x_samples && out_x_samples) x_samples.reserve(M_eff);

    // Metropolis on x
    double x = 0.0; // initial
    int accepted = 0;
    int attempted = 0;

    // thermalization
    for(int i=0; i<burn_in; ++i){
        double xnew = x + rnd.Rannyu(-delta_move, delta_move);
        double a = psi2(xnew, sigma, mu) / psi2(x, sigma, mu);
        if(a >= 1.0 || rnd.Rannyu() < a) { x = xnew; }
    }

    // sampling
    for(int i=0; i<M_eff; ++i){
        double xnew = x + rnd.Rannyu(-delta_move, delta_move);
        double a = psi2(xnew, sigma, mu) / psi2(x, sigma, mu);
        attempted++;
        if(a >= 1.0 || rnd.Rannyu() < a) { x = xnew; accepted++; }

        double e_local = local_energy(x, mu, sigma);
        energies_per_measure.push_back(e_local);
        if(save_x_samples && out_x_samples) x_samples.push_back(x);
    }

    double acc_rate = attempted>0 ? double(accepted)/double(attempted) : 0.0;

    // blocking: compute block means
    int idx = 0;
    for(int b=0; b<n_blocks; ++b){
        double sum=0;
        int count = 0;
        for(int k=0; k<block_size && idx < (int)energies_per_measure.size(); ++k, ++idx){
            sum += energies_per_measure[idx];
            count++;
        }
        if(count>0) block_means.push_back( sum/double(count) );
    }
    // if leftover measures (due to integer division), fold them into last block
    if(idx < (int)energies_per_measure.size()){
        double sum = 0;
        int count = 0;
        for(; idx < (int)energies_per_measure.size(); ++idx, ++count) sum += energies_per_measure[idx];
        if(count>0) {
            if(block_means.empty()) block_means.push_back(sum/double(count));
            else {
                // merge: average with last block weighted
                double last_mean = block_means.back();
                int last_count = block_size;
                double combined = ( last_mean*last_count + sum ) / double(last_count + count);
                block_means.back() = combined;
            }
        }
    }

    // final statistics
    double mean = 0;
    for(double b : block_means) mean += b;
    mean /= double(block_means.size());
    double err = block_error(block_means);

    EnergyResult res;
    res.mean = mean;
    res.error = err;
    res.acceptance_rate = acc_rate;
    if(save_x_samples && out_x_samples) {
        *out_x_samples = x_samples;
    }
    return res;
}

/* ---------------------------
   MAIN: Simulated Annealing driver
   --------------------------- */
int main(){
    // RNG init
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if(!Primes.is_open()){ cerr<<"Unable to open Primes\n"; return 1;}
    Primes >> p1 >> p2;
    Primes.close();
    ifstream input("seed.in");
    if(!input.is_open()){ cerr<<"Unable to open seed.in\n"; return 1;}
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    input.close();
    rnd.SetRandom(seed, p1, p2);

    // ----- SA parameters (tweak these) -----
    int SA_steps = 200;                 // number of SA steps (parameter updates)
    double T0 = 1.0;                    // initial SA temperature
    double alpha = 0.95;                // geometric cooling (T <- alpha * T)
    double sigma_step = 0.05;           // proposal width for sigma (gaussian)
    double mu_step = 0.05;              // proposal width for mu (gaussian)

    // ----- MC sampling parameters for energy estimation (inner loop) -----
    int M_total = 50000;      // total MC steps to estimate <H> per trial (including burn)
    int burn_in = 5000;       // burn-in for each evaluation
    double delta_move = 2.5;  // Metropolis step size on x (tune to get decent acceptance)
    int n_blocks = 50;        // blocks used for error estimate per evaluation

    /*
    // initial variational params
    double mu = 0.78;    // initial guess (from previous 8.1)
    double sigma = 0.66; // initial guess
    */

    // initial variational params
    double mu = 0.5;    // initial guess (from previous 8.1)
    double sigma = 0.1; // initial guess

    // Prepare outputs
    ofstream sa_out("sa_history.dat");
    sa_out << "# step  T_sa   H_avg   H_err   mu   sigma   accept_flag   acc_rate_MC" << endl;

    // compute initial energy
    EnergyResult cur = estimate_energy(rnd, mu, sigma, M_total, burn_in, delta_move, n_blocks);
    double cur_H = cur.mean;

    cout << fixed << setprecision(8);
    cout << "[SA] start: mu="<<mu<<" sigma="<<sigma<<" <H>="<<cur_H<<" err="<<cur.error
         <<" acc_MC="<<cur.acceptance_rate << endl;
    sa_out << 0 << " " << T0 << " " << cur.mean << " " << cur.error << " " << mu << " " << sigma << " " << 1 << " " << cur.acceptance_rate << endl;

    // SA loop
    double Tsa = T0;
    for(int step=1; step<=SA_steps; ++step){
        // propose new parameters (Gaussian proposals)
        double mu_prop = mu + rnd.Gauss(0.0, mu_step);
        double sigma_prop = sigma + rnd.Gauss(0.0, sigma_step);
        if(sigma_prop <= 0.01) sigma_prop = 0.01; // constraint

        // estimate energy at proposed params
        EnergyResult prop = estimate_energy(rnd, mu_prop, sigma_prop, M_total, burn_in, delta_move, n_blocks);
        double prop_H = prop.mean;

        // SA acceptance: Metropolis on parameters with temperature Tsa
        double dE = prop_H - cur_H;
        bool accept = false;
        if(dE < 0) accept = true;
        else {
            double pacc = exp(-dE / Tsa);
            if(rnd.Rannyu() < pacc) accept = true;
        }

        if(accept){
            mu = mu_prop;
            sigma = sigma_prop;
            cur_H = prop_H;
            cur = prop;
        }

        // save history
        sa_out << step << " " << Tsa << " " << cur.mean << " " << cur.error << " "
               << mu << " " << sigma << " " << accept << " " << prop.acceptance_rate << endl;

        // print some progress
        if(step % 10 == 0 || step==1 || step==SA_steps){
            cout << "[SA] step="<<step<<" Tsa="<<Tsa<<" <H>="<<cur.mean<<" +- "<<cur.error
                 <<" mu="<<mu<<" sigma="<<sigma<<" accept="<<(accept? "Y":"N")<<endl;
        }

        // cool down
        Tsa *= alpha;
    }

    sa_out.close();

    // Final optimal parameters
    double mu_opt = mu;
    double sigma_opt = sigma;
    cout << "\n[SA] finished. Optimal: mu="<<mu_opt<<" sigma="<<sigma_opt<<" <H>="<<cur.mean<<" +- "<<cur.error<<endl;

    // ------------------------------
    // Produce final detailed data for the optimal parameters
    // 1) compute <H> and errors vs number of blocks (for plotting)
    // We'll evaluate energies with a large M and then recompute blocking for several L.
    // ------------------------------
    int M_big = 200000; // total samples (including burn)
    int burn_big = 20000;
    double delta_big = delta_move;
    // produce energies and x samples
    vector<double> x_samples;
    EnergyResult final_res = estimate_energy(rnd, mu_opt, sigma_opt, M_big, burn_big, delta_big, 1, true, &x_samples);

    // write samples to file for histogramming
    ofstream samp_out("samples_opt.dat");
    for(double xs : x_samples) samp_out << xs << "\n";
    samp_out.close();

    // blocking analysis: vary L from small to large
    // We'll compute mean and error for different number of blocks Nblocks = M_eff / L
    int M_eff_big = (M_big - burn_big);
    ofstream eb_out("energy_blocks_opt.dat");
    eb_out << "# L   Nblocks   H_mean   H_error" << endl;
    for(int L=10; L<=5000; L = (L<100? L+10 : (L<500? L+50 : L+500))){
        if(L > M_eff_big) break;
        int Nblocks = M_eff_big / L;
        if(Nblocks < 2) continue;
        // compute block means
        vector<double> block_means;
        block_means.reserve(Nblocks);
        int idx = 0;
        for(int b=0; b<Nblocks; ++b){
            double sum=0;
            for(int k=0; k<L; ++k, ++idx) sum += local_energy(x_samples[idx], mu_opt, sigma_opt);
            block_means.push_back(sum/double(L));
        }
        // compute stats
        double mean=0; for(double v: block_means) mean+=v; mean/=double(block_means.size());
        double err = block_error(block_means);
        eb_out << L << " " << Nblocks << " " << mean << " " << err << "\n";
    }
    eb_out.close();

    // produce a coarse histogram file from samples_opt (optional)
    int nbin = 100;
    double xmin = -5.0, xmax = 5.0;
    vector<int> hist(nbin,0);
    for(double xs : x_samples){
        if(xs < xmin || xs >= xmax) continue;
        int b = int( (xs - xmin) / (xmax - xmin) * nbin );
        if(b>=0 && b<nbin) hist[b]++;
    }
    ofstream hist_out("hist_opt.dat");
    hist_out << "# x_mid  count" << endl;
    for(int b=0;b<nbin;++b){
        double xmid = xmin + (b+0.5)*(xmax-xmin)/nbin;
        hist_out << xmid << " " << hist[b] << "\n";
    }
    hist_out.close();

    // Save RNG seed and finish
    rnd.SaveSeed();
    cout << "[DONE] Files written: sa_history.dat samples_opt.dat energy_blocks_opt.dat hist_opt.dat\n";
    return 0;
}
