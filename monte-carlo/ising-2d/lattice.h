/*
 * lattice.h
 *
 */
#ifndef _LATTICE_H
#define _LATTICE_H

#include <unordered_map>     // Clusters.
#include <vector>            // Lattice - a vector of vectors.
#include <gsl/gsl_rng.h>     // Random number generators.
#include <gsl/gsl_sf_exp.h>  // Exponential functions.


class Lattice {
  public:
    Lattice(unsigned int x, unsigned int y, unsigned int seed);
    Lattice();
    ~Lattice();

    /* 
     * Print methods. 
     */
    void PrintLattice();
    void PrintCluster();

    /* 
     * Periodic boundary conditions.
     */
    void getNeighbors(unsigned int site);
    void getHalfNeighbors(unsigned int site);

    /* 
     * Calculate quantities.
     */
    int calculateEnergy(unsigned int site);
    int calculateHalfEnergy(unsigned int site);
    double calculateTotalEnergy();
    double calculateMagnetization();
    double calculateSpecificHeat(double averageE, double squaredE);
    double calculateSusceptibility(double average, double squared);

    /*
     * Simulation algorithms.
     */
    // metropolis returns true if spin flipped, falseotherwise.
    bool metropolis(unsigned int site);
    unsigned int wolff(unsigned int site);
    void growCluster(unsigned int site, int spin);
    void flipCluster();
    void flipComplement();

    /*
     * Data.
     */
    gsl_rng* generator;
    // The cluster map is only meant to keep track of the sites that have
    // already been visited. Thus, only the key matters.
    std::unordered_map<unsigned int, unsigned int> cluster;

    std::vector<int> lattice;  // Magnetic dipole values {+1, -1}.
    unsigned int latticeSize;  // Number lattice sites.
    unsigned int xDim;
    unsigned int yDim;

    // Neighbord lattice sites.
    unsigned int prevX;
    unsigned int nextX;
    unsigned int prevY;
    unsigned int nextY;

    float temp;          // kT in natural units (k == 1).
    float beta;          // Inverse temperature.
    double randomU;      // Random number in the range [0,1).
    double probability;  // Probability of flipping a dipole.

    double exponentials[2];
    double totalEnergy;
};

#endif // _LATTICE_H
