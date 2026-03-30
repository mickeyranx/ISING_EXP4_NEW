#pragma once
#include <vector>
#include <random>
#include <iostream>
using namespace std;
class Simulation {
public:
    
    static mt19937 gen;
    static uniform_real_distribution<double> distr1;
    static uniform_int_distribution<int> distr2;
    static random_device rand_device;

    static unsigned int try_count;
    static unsigned int acceptance_count;
    

    static vector<int> initializeLatticeCold(int L);

    static vector<int> initializeLatticeHot(int L);

    //return a small list of 4 indices of Positions of the neighbors (top, right, bottom, left)
    static vector<int> getNeighborPos(int i, int L);

    static double changeInEnergy(vector<int>& config_1, int L, int pos, int s, double h);

    static void sweepMetropolis(vector<int>& config_1,int L ,double beta, double h);

    static void sweepMetropolisMultihit(vector<int>& config_1, int L, double beta, double h, int tries);

    static void sweepHeatbath(vector<int>& config_1, int L, double beta, double h);

    //this is method is for Metropolis multihit, it performs <draw_interval> sweeps
    static void draw(vector<int>& config, int L, double beta, double h, int draw_interval, int multi_hit);

    //this method is for Heathbath, it performs <draw_interval> sweeps
    static void draw(vector<int>& config, int L, double beta, double h, int draw_interval);

    static double averageEnergy(vector<int>& config, int K, int L, double h);


    static double averageMagnetisation(int M, vector<int>& config);

    

    static void printConfig(vector<int>& config);

};

