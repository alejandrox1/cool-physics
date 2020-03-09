/*
 * lattice.cpp
 */
#include "lattice.h"


Lattice::Lattice(unsigned int x, unsigned int y, unsigned int seed) {
  generator = gsl_rng_alloc(gsl_rng_mt19937);  // Mersenne Twister.
  gsl_rng_set(generator, seed);

  temp = (float)seed / 100;
  beta = 1.0 / temp;
  xDim = x;
  yDim = y;
  latticeSize = x * y;

  // Instantiate lattice.
  for (unsigned int i = 0; i < latticeSize; i++) {
    if (gsl_rng_uniform(generator) < 0.5)
      lattice.push_back(-1);
    else
      lattice.push_back(1);
  }

  // Calculate exponential factors.
  exponentials[0] = gsl_sf_exp(-beta * 4);
  exponentials[1] = gsl_sf_exp(-beta * 8);

  randomU = gsl_rng_uniform(generator);
  probability = 1 - gsl_sf_exp(-2 * beta);

  totalEnergy = calculateTotalEnergy();
}

Lattice::Lattice() {
  Lattice(32, 32, 227);
}

Lattice::~Lattice() {}

/*
 * Print methods.
 */
void Lattice::PrintLattice() {
  for (unsigned int i = 0; i < latticeSize; i++) {
    if (i % xDim == 0)
      printf("\n");

    if (lattice[i] == -1)
      printf(" ");
    else if (lattice[i] == 1)
      printf("x");
    else
      printf("ERROR");
  }
  printf("\n");
  fflush(stdout);
}

void Lattice::PrintCluster() {
  for (unsigned int i = 0; i < latticeSize; i++) {
    if (i % xDim == 0)
      printf("\n");

    // If key not found in map iterator to end is returned.
    if (cluster.find(i) != cluster.end())
      printf("x");
    else
      printf(" ");
  }
  printf("\n");
  fflush(stdout);
}

/*
 * Periodic Boundary Conditions. 
 */
// getHalfNeighbors with helical periodic boundary conditions.
void Lattice::getHalfNeighbors(unsigned int site) {
  if (site < latticeSize - xDim) {
    nextX = site + 1;
    nextY = site + xDim;
  } else if (site < latticeSize - 1) {
    nextX = site + 1;
    nextY = site + xDim - latticeSize;
  } else {  // site == latticeSize - 1;
    nextX = 0;
    nextY = xDim;
  }
}

void Lattice::getNeighbors(unsigned int site) {
  getHalfNeighbors(site);

  if (site > xDim) {
    prevX = site - 1;
    prevY = site - xDim;
  } else if (site > 0) {
    prevX = site - 1;
    prevY = site + latticeSize - xDim;
  } else {  // site = 0;
    prevX = latticeSize - 1;
    prevY = latticeSize - xDim;
  }
}

/*
 * Physical Quantitites.
 */
int Lattice::calculateHalfEnergy(unsigned int site) {
  getHalfNeighbors(site);

  return -lattice[site] * (lattice[nextX] + lattice[nextY]);
}

int Lattice::calculateEnergy(unsigned int site) {
  getNeighbors(site);

  return -lattice[site] * (lattice[prevX] + lattice[nextX] + lattice[prevY] + lattice[nextY]);
}

double Lattice::calculateTotalEnergy() {
  double totalEnergy = 0.0;
  for (unsigned int i = 0; i < latticeSize; i++)
    totalEnergy += calculateHalfEnergy(i);
  return totalEnergy / latticeSize;
}

double Lattice::calculateMagnetization() {
  double magnetization = 0.0;
  for (unsigned int i = 0; i < latticeSize; i++)
    magnetization += lattice[i];
  return magnetization / latticeSize;
}

double Lattice::calculateSpecificHeat(double averageE, double squaredE) {
  double diff = squaredE - (averageE * averageE);
  return beta * beta * diff * latticeSize;
}

double Lattice::calculateSusceptibility(double averageM, double squaredM) {
  double diff = squaredM - (averageM * averageM);
  return beta * diff * latticeSize;
}

/*
 * Simulation Algorithms.
 */
bool Lattice::metropolis(unsigned int site) {
  // This is half the change.
  int difference = -calculateEnergy(site);
  //printf("%d %d %0.8lf\n", difference, difference / 2 - 1, probability);

  if (difference <= 0.0) {  // Lower energy encountered, so flip.
    lattice[site] *= -1;
    return true;
  } else {  // Probabilistic acceptance.
    randomU = gsl_rng_uniform(generator);
    probability = exponentials[difference / 2 - 1];

    if (randomU < probability) {
      lattice[site] *= -1;
      return true;
    } else {
      return false;
    }
  }
}

void Lattice::growCluster(unsigned int site, int spin) {
  getNeighbors(site);
  unsigned int curPrevX = prevX;
  unsigned int curNextX = nextX;
  unsigned int curPrevY = prevY;
  unsigned int curNextY = nextY;

  if (lattice[curPrevX] == spin && (cluster.find(curPrevX) == cluster.end())) {
    randomU = gsl_rng_uniform(generator);
    if (randomU < probability) {
      cluster[curPrevX] = curPrevX;
      growCluster(curPrevX, spin);
    }
  }

  if (lattice[curNextX] == spin && (cluster.find(curNextX) == cluster.end())) {
    randomU = gsl_rng_uniform(generator);
    if (randomU < probability) {
      cluster[curNextX] = curNextX;
      growCluster(curNextX, spin);
    }
  }

  if (lattice[curPrevY] == spin && (cluster.find(curPrevY) == cluster.end())) {
    randomU = gsl_rng_uniform(generator);
    if (randomU < probability) {
      cluster[curPrevY] = curPrevY;
      growCluster(curPrevY, spin);
    }
  }

  if (lattice[curNextY] == spin && (cluster.find(curNextY) == cluster.end())) {
    randomU = gsl_rng_uniform(generator);
    if (randomU < probability) {
      cluster[curNextY] = curNextY;
      growCluster(curNextY, spin);
    }
  }
}

void Lattice::flipCluster() {
  std::unordered_map<unsigned int, unsigned int>::iterator got;
  for (got = cluster.begin(); got != cluster.end(); ++got)
    got->second *= -1;
}

// flipComplement flips the spins of all sites that are not part of the
// cluster.
void Lattice::flipComplement() {
  for (unsigned int i = 0; i < latticeSize; i++)
    if (cluster.find(i) == cluster.end())
      lattice[i] *= -1;
}

// wolff returns the size of the cluster.
unsigned int Lattice::wolff(unsigned int site) {
  cluster[site] = site;
  growCluster(site, lattice[site]);

  if (cluster.size() > (latticeSize / 2))
    flipComplement();
  else
    flipCluster();

  unsigned int size = cluster.size();
  return size;
}
