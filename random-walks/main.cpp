/*
 * Random walk in 3D
 *
 * This is a small experiment to test to empirically observe the behaviour of a
 * random walk.
 * Specifically, we want to observe the following condition:
 *
 *    < x^2 > = N * l^2
 *
 *   Where x is the distance, N is the number of steps, and l is the step
 *   length.
 *
 *  More specifically, the above expression tells us that the average of the
 *  square distance grows linearly with the number of steps take.
 *  
 *  We want to observe this by conductin a number of experiments in which we
 *  take N number of steps.
 *
 * This program will output a CSV file of where each line will be the number of
 * steps taken and the average square distance from the origin (0,0).
 */
#include <iostream>
#include <gsl/gsl_rng.h>

int main() {
  FILE* output;
  output = fopen("walk.csv", "w");

  int sampleSize = 100000;
  int minLength = 1;
  int maxLength = 150;
  int stepLength = 1;

  double dist = 0;
  double distTotal = 0;
  double distAve = 0;

  // Mersenne Twister Generator.
  gsl_rng* generator = gsl_rng_alloc(gsl_rng_mt19937);

  int direction = 0;  // Used to determine the direction of a step.
  int x = 0;
  int y = 0;
  int z = 0;

  // Here we want to see the relationship between the number of steps taken and
  // the average distance traveled.
  // We will do a 'maxLength' number of random walks from 'minLength'-sized
  // random walks, to 'maxLength'-sized random walks.
  for (int totalSteps = minLength; totalSteps <= maxLength; totalSteps++) {
    distTotal = 0;

    // For each 'totalSteps'-sized random walk we will perform 'sampleSize'
    // number of them and take the average square distance.
    for (int test = 0; test < sampleSize; test++) {
      x = 0;
      y = 0;
      z = 0;

      // Here we walk the walk.
      // We take 'totalSteps' in any of the 3 dimensions (6 directions).
      // All steps ae of length 'stepLength'.
      for (int step = 0; step < totalSteps; step++) {
        direction = gsl_rng_get(generator) % 6;
        if      (direction == 0) x += stepLength;
        else if (direction == 1) y += stepLength;
        else if (direction == 2) z += stepLength;
        else if (direction == 3) x -= stepLength;
        else if (direction == 4) y -= stepLength;
        else if (direction == 5) z -= stepLength;
      }

      // Compute distance square for each walk.
      dist = x*x + y*y + z*z;
      distTotal += dist;
    }

    // Compute average distance taken in all walks.
    distAve = distTotal / sampleSize;
    std::cout << totalSteps << " " << distAve << " " << distAve/totalSteps << "\n";

    fprintf(output, "%i,\t%lf\n", totalSteps, distAve);
  }

  gsl_rng_free(generator);
  fclose(output);
  return 0;
}
