/*
 * Metropolis.cpp
 *
 * Run Monte Carlo simulation using the metropolis algorithm.
 *
 * Write out the energy, the magnetization, the specific heat, and the
 * susceptibility.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_log.h>
#include "lattice.h"


int main(int argc, char** const argv) {
  if (argc != 6) {
    fprintf(stderr, "Usage %s xDim yDim equil steps temp\n", argv[0]);
    fflush(stderr);
    exit(1);
  }

  unsigned int xDim   = atoi(argv[1]);  // Lattice size in the X dimension.
  unsigned int yDim   = atoi(argv[2]);  // Lattice size in the Y dimension.
  unsigned int equil  = atoi(argv[3]);  // Equlibration steps.
  unsigned int steps  = atoi(argv[4]);  // Number of iterations.
  unsigned int seed   = atoi(argv[5]);  // 100x Temperature.

  float temp = (float) seed / 100;
  unsigned int latticeSize = xDim * yDim;

  double aveE           = 0.0;
  double aveM           = 0.0;
  double aveAbsM        = 0.0;
  double sqrE           = 0.0;
  double sqrM           = 0.0;
  double specificHeat   = 0.0;
  double susceptibility = 0.0;

  double energyData[steps];
  double magnetizationData[steps];
  double autocorrelation[steps];
  double scaleFactor          = 0.0;
  double autocorrelationTime  = 0.0;
  double stdE                 = 0.0;
  double stdM                 = 0.0;

  Lattice* lattice = new Lattice(xDim, yDim, seed);

  // Initialize and equilibrate the lattice.
  // We perform a metropolis 'step' at random trying to get every member of the
  // lattice an 'equil' number of times.
  unsigned int randomSite;
  for (unsigned int i = 0; i < equil * latticeSize; i++) {
    randomSite = (int) floor(latticeSize * gsl_rng_uniform(lattice->generator));
    lattice->metropolis(randomSite);
  }

  // Make measurements every N 'latticeSize' cycles (this is the number of laps).
  // This should be based on the autocorrelation time.
  unsigned int laps = 5;
  unsigned int counter  = 0;
  for(unsigned int i = 0; i < steps * latticeSize * laps; i++) {
    randomSite = (int) floor(latticeSize * gsl_rng_uniform(lattice->generator));
    lattice->metropolis(randomSite);

    // Save data every N laps.
    if (i % (latticeSize * laps) == 0) {
      energyData[counter] = lattice->calculateTotalEnergy();
      aveE               += energyData[counter];

      magnetizationData[counter]    = lattice->calculateMagnetization();
      aveM                         += magnetizationData[counter];
      magnetizationData[counter]    = fabs(magnetizationData[counter]);
      aveAbsM                      += magnetizationData[counter];

      sqrE += energyData[counter] * energyData[counter];
      sqrM += magnetizationData[counter] * magnetizationData[counter];

      counter++;
    }
  }

  // Average all data.
  aveE     /= steps;
  aveM     /= steps;
  aveAbsM  /= steps;
  sqrE /= steps;
  sqrM /= steps;

  // Do some bootstrapping.
  specificHeat   = lattice->calculateSpecificHeat(aveE, sqrE);
  susceptibility = lattice->calculateSusceptibility(aveAbsM, sqrM);

  lattice->PrintLattice();
  delete lattice;

  // Generate Chi[0] for scaling.
  for (unsigned int i = 0; i < steps; i++)
    scaleFactor += (magnetizationData[i] * magnetizationData[i]);
  scaleFactor /= steps;
  scaleFactor -= (aveAbsM * aveAbsM);

  autocorrelation[0] = 1.0;

  // Calculate autocorrelation function using the magnetization.
  for (unsigned int t = 1; t < steps; t++) {
    autocorrelation[t] = 0.0;

    for (unsigned int i = 0; i < steps - t; i++)
      autocorrelation[t] += (magnetizationData[i] * magnetizationData[i+t]);
    autocorrelation[t] /= (steps - t);

    autocorrelation[t] -= (aveAbsM * aveAbsM);
    autocorrelation[t] /= scaleFactor;
  }

  // Generate autocorrelation time.
  double tempd;
  unsigned int i = 1;
  while (autocorrelation[i] > 0
      && autocorrelation[i] < autocorrelation[i-1]
      && i < steps) {
    
    tempd = -gsl_sf_log(autocorrelation[i]);
    tempd = i / tempd;
    autocorrelationTime += tempd;
  
    i++;
  }

  // Calculate standard deviations.
  if (i == 1) {  // Assume the autocorrelation time is 0.
    autocorrelationTime = 0.0;
    stdE = 0.0;
    stdM = 0.0;
  } else {
    autocorrelationTime /= (i - 1);

    stdE  = 2 * autocorrelationTime / steps;
    stdE += sqrE - aveE * aveE;
    stdE  = sqrt(stdE);

    stdM  = 2 * autocorrelationTime / steps;
    stdM += sqrM - aveAbsM * aveAbsM;
    stdM  = sqrt(stdM);
  }

  printf("%f\t%lf\t", temp, autocorrelationTime);
  printf("%lf\t%lf\t", aveE, stdE);
  printf("%lf\t%lf\t", aveAbsM, stdM);
  printf("%lf\t%lf\t", specificHeat, susceptibility);
  printf("%lf\t%lf\n", aveM, scaleFactor);

  return 0;
}
