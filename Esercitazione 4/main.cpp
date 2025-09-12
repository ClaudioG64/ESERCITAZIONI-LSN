#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;


double error(double av, double av2, int n) {
    if (n == 0) return 0;
    return sqrt((av2 - av * av) / n);
}
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


   //ESERCIZIO 1
   int throwsperblock = 1000;  // 10^3 lanci per blocco

   ofstream output("output.dat");  // File di output

   for (int i = 0; i < 50; i++) {
      int nblocks = 200 * (i + 1);  // Numero di blocchi
      int M = nblocks * throwsperblock;  // Numero totale di lanci

      double sum_r = 0, sum_r2 = 0;  // Per media e varianza

      for (int b = 0; b < nblocks; b++) {
         double r_i = 0;  // Media del blocco

         for (int j = 0; j < throwsperblock; j++) {
            r_i += rnd.Rannyu();
         }
         r_i /= throwsperblock;  // Media di un blocco

         sum_r += r_i;
         sum_r2 += r_i * r_i;
      }

      double mean_r = sum_r / nblocks;
      double mean_r2 = sum_r2 / nblocks;
      double sigma = error(mean_r, mean_r2, nblocks);

      cout << "Iterazione: " << i + 1 
         << " | Blocchi: " << nblocks 
         << " | M: " << M
         << " | Media: " << mean_r 
         << " | Errore: " << sigma << endl;

      output << nblocks << " " << mean_r << " " << sigma << endl;
    }

   output.close();
   rnd.SaveSeed();
    
    return 0;
}
