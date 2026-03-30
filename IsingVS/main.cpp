#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <set>
#include <tuple>
#include <iomanip>
#include "Simulation.h"
#include <random>
#include <string>
#include <fstream>
#include <numeric>
#include <format> //works only with C++20
using namespace std;



// RNG with the box m�ller method
static double boxMueller(double mu, double sigma) {
    //Zufallszahlen aus GLeichverteilung
    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_real_distribution<double> unif_distr(0.0, 1.0);
    double m = unif_distr(gen);
    double n = unif_distr(gen);
    double r = sigma * sqrt(-2.0 * log(m));

    return r * cos(2.0 * M_PI * n) + mu;

}

// approximation of a certain gaussian Integral with MC Integration (A1 b)
static double approxGaussianIntegral(int sample_size) {
    double upper = 1.5;
    double lower = -1.5;
    double sigma = 1;
    double sum = 0;
    set<double> x_values;
    //calculate with box-m�ller
    for (int i = 0; i < sample_size; i++)
    {
        double x = boxMueller(0, sigma);
        while (x < lower || x > upper) { //only allow samples within boundaries
            x = boxMueller(0, sigma);
        }
        x_values.insert(x);
    }
    //cout << x_values.size() << endl;
    for (double x : x_values)
    {
        double n = exp(-pow(x, 2) / 2);
        sum += n;
    }

    //approximate with integration and summs
    return (upper - lower) * sum/sample_size;
    
}
//approximation of PI witch MC (A1 a)
static double approximatonOfPi(int sample_size) {
    random_device rand_dev; 
    mt19937 gen(rand_dev());
    uniform_real_distribution<double> unif_distr(-1.0, 1.0);
    int square = 0;
    int circle = 0;
    for (int i = 0; i < sample_size; i++)
    {
        double x = unif_distr(gen);
        double y = unif_distr(gen);
        if (pow(x,2) + pow(y,2)  <= 1) {
            circle++;
            
        }
        square++;

    }

    return 4 * circle / double(square);

}

//backtrack algorithm to invoke all possible configurations of a LxL lattice 
static void backtrack(int L, int counter, vector<int> config, vector<vector<int>>& list_of_configs) {
    if (counter >= pow(L, 2)) { //if the counter reaches the maximum amount of spins possible
        list_of_configs.push_back(config);
        return;
    }
    config.push_back(1);
    counter++;
    backtrack(L, counter, config, list_of_configs); //explore possibilites for 1 at this position until limit is reached
    
    
    counter--;
    config.pop_back(); //remove the previouls inserted 1 at this position
    config.push_back(-1);
    counter++;
    backtrack(L, counter, config, list_of_configs);//explore possibilites for -1 at this position until limit is reached
    return;

}

//calculates the mean energy pp, mean absolute magnetism pp 
//and mean magnetism pp explicitly with the partition function over all possible configurations
static void explicitIsing(int L) {
    //-----------------------------------
    //               setup
    //-----------------------------------
    //list of all possible configurations
    vector<vector<int>> list_of_configs = {};
    //create all configs with backtrack-algorithm
    backtrack(L, 0, {}, list_of_configs);
    cout << list_of_configs.size() << endl;
    //initialize list of betas
    vector<double> betas = {0.1, 0.2, 0.3, 0.35, 0.4, 0.415, 0.44, 0.4406868, 0.445, 0.45, 0.465, 0.48, 0.5, 0.55, 0.6, 0.7, 0.9};
    double K = (double) L * L;
    //-----------------------------------
    //      calculate observables
    //-----------------------------------
    ofstream File("explicit_vals_L=" + to_string(L) + ".txt");
    File << "beta" << "\t" << "<e>" << "\t" << "<m>" << "\t" << "<|m|>" << "\n";
    File << fixed << setprecision(7);
    for (double beta : betas) {
        //1.partition-function
        double Z = 0;
        vector<double> energies = {};
        for (vector<int> config : list_of_configs) {
            double H = Simulation::averageEnergy(config, L*L,L , 0 )*L*L;
            energies.push_back(H);
            Z += exp(-beta * H);
        }

        //2. mean energy pp, mean magnetism pp, mean absolute magnetism pp
        double mean_energy = 0;
        double mean_mag = 0;
        double mean_abs_mag = 0;
        int i = 0;
        for (vector<int> config : list_of_configs) {
            double H_i = energies[i];
            mean_energy += 1.0 / Z * exp(-beta * H_i) * H_i / ((double)L * L);
            double M_i = Simulation::averageMagnetisation(K, config);
            mean_mag += 1.0 / (K*Z) *  M_i * exp(-beta * H_i);
            mean_abs_mag += 1.0 / (K*Z) * abs(M_i) * exp(-beta * H_i);
            i++;
        }

        File << beta << "\t" << mean_energy << "\t" << mean_mag << "\t" << mean_abs_mag << "\n";


    }
    
    File.close();





}


//Simulation for exercises 3 and 4
static vector<double> startSimulationMetropolis(int L, double beta, double h,int therm_steps, int N, int draw_interval, bool hot_start, int multihit) {
    //------------------------------
    //            setup
    //------------------------------
    vector<int> config = {};
    if (hot_start) {
        config = Simulation::initializeLatticeHot(L);

    }
    else {
        config = Simulation::initializeLatticeCold(L);
    }
    
   
    for (int i = 0; i < therm_steps; i++) //thermalisation
    {
        Simulation::sweepMetropolisMultihit(config, L, beta, h, multihit);
    }
    cout << "therm done" << endl;
    //Observables
    vector<double> energies1(N, 0);
    vector<double> energies2(N, 0);
    vector<double> absmag1(N, 0);
    vector<double> absmag2(N, 0);
    vector<double> mag1(N, 0);
    vector<double> mag2(N, 0);
    //------------------------------
    //            sweeps
    //------------------------------
    int K = L * L;
    int counter = 0;
    for (int i = 0; i < N; i++)
    {
        Simulation::draw(config, L ,beta , h , draw_interval, multihit);
        double E = Simulation::averageEnergy(config, K, L, h);
        energies1[i] = E;
        energies2[i] = E * E;
        double M = Simulation::averageMagnetisation(L * L, config);
        double absmag_i = abs(M) / K;
        absmag1[i] = absmag_i;
        absmag2[i] = absmag_i * absmag_i;
        mag1[i] = M / K;
        mag2[i] = M / K * M / K;
        
        cout << counter << endl;

        counter++;
    }
   
    //------------------------------
    //    calculate observables
    //------------------------------
    vector<double> vals = {};

    double e_mean = accumulate(energies1.begin(), energies1.end(), 0.0) / N; //accumulate(start, end) is a sum over all members
    double e2_mean = accumulate(energies2.begin(), energies2.end(), 0.0) / N;
    double m_mean = accumulate(mag1.begin(), mag1.end(), 0.0) / N;
    double m2_mean = accumulate(mag2.begin(), mag2.end(), 0.0) / N;
    double absm_mean = accumulate(absmag1.begin(), absmag1.end(), 0.0) / N;
    double absm2_mean = accumulate(absmag2.begin(), absmag2.end(), 0.0) / N;
    vals.push_back(e_mean);
    vals.push_back(sqrt((e2_mean - e_mean * e_mean) / (N - 1)));
    vals.push_back(m_mean);
    vals.push_back(sqrt((m2_mean - m_mean * m_mean) / (N - 1)));
    vals.push_back(absm_mean);
    vals.push_back(sqrt((absm2_mean - absm_mean * absm_mean) / (N - 1)));
    vals.push_back(beta * beta * (e2_mean - e_mean * e_mean)); //specific heat
    vals.push_back(e2_mean);
    vals.push_back(m2_mean);
    return vals;

   
}


//Simulation for exercises 3 and 4
static vector<double> startSimulationHeatbath(int L, double beta, double h, int therm_steps, int N, int draw_interval, bool hot_start) {
    //------------------------------
    //            setup
    //------------------------------
    vector<int> config = {};
    if (hot_start) {
        config = Simulation::initializeLatticeHot(L);

    }
    else {
        config = Simulation::initializeLatticeCold(L);
    }


    for (int i = 0; i < therm_steps; i++)
    {
        Simulation::sweepHeatbath(config, L, beta, h);
    }

    //Observables
    vector<double> energies1(N, 0);
    vector<double> energies2(N, 0);
    vector<double> absmag1(N, 0);
    vector<double> absmag2(N, 0);
    vector<double> mag1(N, 0);
    vector<double> mag2(N, 0);
    //------------------------------
    //            sweeps
    //------------------------------
    int K = L * L;
    for (int i = 0; i < N; i++)
    {
        Simulation::draw(config, L, beta, h, draw_interval);
        double E = Simulation::averageEnergy(config, K, L, h);
        energies1[i] = E;
        energies2[i] = E * E;
        double M = Simulation::averageMagnetisation(L * L, config);
        double absmag_i = abs(M) / K;
        absmag1[i] = absmag_i;
        absmag2[i] = absmag_i * absmag_i;
        mag1[i] = M / K;
        mag2[i] = M / K * M / K;


    }
   

    //------------------------------
    //    calculate observables
    //------------------------------
    vector<double> vals = {};

    double e_mean = accumulate(energies1.begin(), energies1.end(), 0.0) / N; //accumulate(start, end) is a sum over all members
    double e2_mean = accumulate(energies2.begin(), energies2.end(), 0.0) / N;
    double m_mean = accumulate(mag1.begin(), mag1.end(), 0.0) / N;
    double m2_mean = accumulate(mag2.begin(), mag2.end(), 0.0) / N;
    double absm_mean = accumulate(absmag1.begin(), absmag1.end(), 0.0) / N;
    double absm2_mean = accumulate(absmag2.begin(), absmag2.end(), 0.0) / N;
    vals.push_back(e_mean);
    vals.push_back( sqrt( (e2_mean - e_mean * e_mean) / (N - 1) ));
    vals.push_back(m_mean);
    vals.push_back(sqrt( (m2_mean - m_mean * m_mean) / (N - 1) ));
    vals.push_back(absm_mean);
    vals.push_back(sqrt( (absm2_mean - absm_mean * absm_mean) / (N - 1) ));
    vals.push_back(beta * beta * (e2_mean - e_mean * e_mean)); //specific heat
    vals.push_back(e2_mean);
    vals.push_back(m2_mean);
    return vals;

}

//testing verions of startSimulation(..) for optimal parameters
static void parameterTesting(string output_filename, int L, double beta, double h, int N, int draw_interval, bool hot_start, bool Metropolis, int multi_hit, int block_size, double epsilon, int therm_time ) {
    //------------------------------
    //            setup
    //------------------------------
    vector<int> config = {};
    if (hot_start) { //initialise 
        config = Simulation::initializeLatticeHot(L);

    }
    else {
        config = Simulation::initializeLatticeCold(L);
    }


    for (int i = 0; i < therm_time; i++)//thermalisation steps
    {
        cout << "i" << endl;
        Simulation::draw(config, L, beta, h, draw_interval, multi_hit);
    }


    ofstream File(output_filename);
    File << fixed << setprecision(5);
    File << "beta" << "\t" << "e" << "\t" << "m" << "\t" << "|m|" << "\n";
    //------------------------------
    //            sweeps
    //------------------------------
    int K = L * L; //no. lattice points

    if (Metropolis) {
        //initialise block means for energy and absolute value of magnetisation
        double prev_e_mean = 0;
        double prev_abs_m_mean = 0;
        clock_t start_2 = clock();

        for (int i = 0; i < N; i++)// N x block_size is the total amount of sweeps
        {
            
            
            vector<double> e_block(block_size, 0); //initialise block vectors
            vector <double> abs_m_block(block_size, 0);
            

            for (int j = 0; j < block_size; j++) { //block

                Simulation::draw(config, L, beta, h, draw_interval, multi_hit);
                double e = Simulation::averageEnergy(config, K, L, h); 
               
                double M = Simulation::averageMagnetisation(L * L, config);
                double absmag_i = abs(M) / K;
                

                double m = M / K;
                
                e_block[j] = e;
                abs_m_block[j] = absmag_i;

                //writing the important data for testing thermalisation to a file
                File << beta << "\t" << e << "\t" << m << "\t" << absmag_i << "\n";


            }

            double next_e_mean = accumulate(e_block.begin(), e_block.end(), 0.0) / block_size;
            double next_abs_m_mean = accumulate(abs_m_block.begin(), abs_m_block.end(), 0.0) / block_size;
            
           
            if (i > 0) {
                cout << "comparing " << next_abs_m_mean << " with " << prev_abs_m_mean << " : " << abs(next_abs_m_mean - prev_abs_m_mean) << endl;
                clock_t end_2 = clock();
                double elapsed = double(end_2 - start_2) / CLOCKS_PER_SEC;
                cout << elapsed << endl;
                if (abs(next_abs_m_mean - prev_abs_m_mean) < epsilon && abs(next_e_mean - prev_e_mean) < epsilon) {
                    cout << "convergence reached" << to_string(i * block_size) << endl;
                   
                }

            }
            
            prev_abs_m_mean = next_abs_m_mean;
            prev_e_mean = next_e_mean;

        }
       
    }
    else
    {
        //initialise block means for energy and absolute value of magnetisation
        double prev_e_mean = 0;
        double prev_abs_m_mean = 0;
        clock_t start_2 = clock();

        for (int i = 0; i < N; i++)// N x block_size is the total amount of sweeps
        {


            vector<double> e_block(block_size, 0); //initialise block vectors
            vector <double> abs_m_block(block_size, 0);


            for (int j = 0; j < block_size; j++) { //block

                Simulation::draw(config, L, beta, h, draw_interval, multi_hit);
                double e = Simulation::averageEnergy(config, K, L, h);

                double M = Simulation::averageMagnetisation(L * L, config);
                double absmag_i = abs(M) / K;


                double m = M / K;

                e_block[j] = e;
                abs_m_block[j] = absmag_i;

                //writing the important data for testing thermalisation to a file
                File << beta << "\t" << e << "\t" << m << "\t" << absmag_i << "\n";


            }

            double next_e_mean = accumulate(e_block.begin(), e_block.end(), 0.0) / block_size;
            double next_abs_m_mean = accumulate(abs_m_block.begin(), abs_m_block.end(), 0.0) / block_size;


            if (i > 0) {
                cout << "comparing " << next_abs_m_mean << " with " << prev_abs_m_mean << " : " << abs(next_abs_m_mean - prev_abs_m_mean) << endl;
                clock_t end_2 = clock();
                double elapsed = double(end_2 - start_2) / CLOCKS_PER_SEC;
                cout << elapsed << endl;
                if (abs(next_abs_m_mean - prev_abs_m_mean) < epsilon && abs(next_e_mean - prev_e_mean) < epsilon) {
                    cout << "convergence reached" << to_string(i * block_size) << endl;

                }

            }

            prev_abs_m_mean = next_abs_m_mean;
            prev_e_mean = next_e_mean;

        }
    }

    File.close();
    return;

}





int main()
{
    clock_t start = clock();
    //--------------------------------------------
    //                 exercise 1
    //--------------------------------------------
    /*
    printf("approximation of PI: %.3f", approximatonOfPi(10000));
    printf("\napproximation of gaussian integral: %.4f", approxGaussianIntegral(50000));
    */

    //--------------------------------------------
    //                 exercise 2
    //--------------------------------------------
    /*
    explicitIsing(2);
    explicitIsing(3);
    explicitIsing(4);
    */
    //--------------------------------------------
    //          exercises 3 & 4
    //--------------------------------------------
    
    /* 
    to alter algorithm of choice 
    use startSimulationMetropolis() and or startSimulationHeathbath
    */

    /*
    //temperature
    //double beta = 0.55;
    //external magnetic field
    double h = 0.0;
    //number of thermalize sweeps
    int therm_steps = 0;
    //sweeps between drawing
    int draw_interval = 1;
    //lattice size
    int L = 128; //actual size is LxL
    //number of draws
    int N = 200; //actual number of sweeps is draw_interval * N
    vector<double> betas = {0.44}; //beta with coldstart
    //vector<double> betas = { 0.1, 0.2, 0.3,0.35 }; //bet with hotstart
    //vector<double> betas = {0.5,0.65, 0.8, 0.95}; //beta near critial coldstart
    

    string output_path = "output_files/testing_therm_params/testing_metro_L=" + std::to_string(L) + "_";
    int mh = 2; //multihit parameter
    bool hot = false;
    for (double beta : betas) {

        parameterTesting(output_path + std::format("beta={:.3f}_MH={}_start={}.txt", beta, mh, hot), L, beta, h, N, draw_interval, hot, true, mh, 100,0.001, therm_steps);
        


    }
    cout << (double) Simulation::acceptance_count/ (double) Simulation::try_count << endl;
    
    */
    //temperature
    //double beta = 0.55;
    //external magnetic field
    double h = 0.0;
    //number of thermalize sweeps
    int therm_steps = 20000;
    //sweeps between drawing
    int draw_interval = 200;
    //lattice size
    int L = 128; //actual size is LxL
    //number of draws
    int N = 100; //actual number of sweeps is draw_interval * N
    //vector<double> betas = {0.48, 0.46, 0.45, 0.445,0.44, 0.42, 0.40};
    //vector<int> seeds = { 9381, 76372, 91023, 81237, 381237, 81238, 12388, 81238, 2314, 63123, 381, 238, 23, 7653 , 9139, 3241};
    vector<double> betas = { 0.44 };
    vector<double> seeds = {642891, 3169};
    int multi_hit = 1;
    bool hot = true;
    //vector<double> betas = {0.4406868};
    int counter = 0;
    for (double beta : betas)
    {
        Simulation::gen = mt19937(seeds[counter]);
        string path = "output_files/metro_output/";
        string filename = to_string(L) + "_Metropolis_beta=" + to_string(beta) + "_start=" + to_string(hot) + ".txt";
        vector<double> vals = startSimulationMetropolis(L, beta, h, therm_steps, N, draw_interval, hot, multi_hit);
        ofstream File(path + filename);
        File << fixed << setprecision(12);
        File << "beta" << "\t" << "<e>" << "\t" << "de" << "\t" << "<m>" << "\t" << "dm" << "\t";
        File << "<|m|>" << "\t" << "d|m|" << "\t" << "c_v/(L^2)" << "\t" << "<e^2>" << "\t" << "<m^2>" << "\t";
        File << "t_therm" << "\t" << "draw_interval" << "\t" << "draws" << "\t" << "seed" << "\n";
        File << beta << "\t" << vals[0] << "\t" << vals[1] << "\t" << vals[2] << "\t" << vals[3] << "\t" << vals[4] << "\t";
        File << vals[5] << "\t" << vals[6] << "\t" << vals[7] << "\t" << vals[8] << "\t";
        File << therm_steps << "\t" << draw_interval << "\t" << N << "\t" << seeds[counter] << "\n";
        File.close();

        
        counter++;
        
    }
    hot = false;
    for (double beta : betas)
    {
        Simulation::gen = mt19937(seeds[counter]);
        string path = "output_files/metro_output/";
        string filename = to_string(L) + "_Metropolis_beta=" + to_string(beta) + "_start=" + to_string(hot) + ".txt";
        vector<double> vals = startSimulationMetropolis(L, beta, h, therm_steps, N, draw_interval, hot, multi_hit);
        ofstream File(path + filename);
        File << fixed << setprecision(12);
        File << "beta" << "\t" << "<e>" << "\t" << "de" << "\t" << "<m>" << "\t" << "dm" << "\t";
        File << "<|m|>" << "\t" << "d|m|" << "\t" << "c_v/(L^2)" << "\t" << "<e^2>" << "\t" << "<m^2>" << "\t";
        File << "t_therm" << "\t" << "draw_interval" << "\t" << "draws" << "\t" << "seed" << "\n";
        File << beta << "\t" << vals[0] << "\t" << vals[1] << "\t" << vals[2] << "\t" << vals[3] << "\t" << vals[4] << "\t";
        File << vals[5] << "\t" << vals[6] << "\t" << vals[7] << "\t" << vals[8] << "\t";
        File << therm_steps << "\t" << draw_interval << "\t" << N << "\t" << seeds[counter] << "\n";
        File.close();
        
        counter++;

    }
   
    
    

   
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    printf("execution time: %.3f sec", elapsed);
    
}

