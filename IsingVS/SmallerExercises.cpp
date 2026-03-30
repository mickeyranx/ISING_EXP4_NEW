#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <set>
using namespace std;


class SmallerExercises {
public:
    static vector<int> Prng(int seed, int iterations);
    static double approxPi(int N);

};

vector<int> SmallerExercises::Prng(int seed, int iterations) {
    int seed_length = to_string(seed).length();

    vector<int> output;

    //string number = to_string(seed);
    int squared = pow(seed, 2);
    string number = to_string(squared);

    for (int i = 0; i < iterations; i++)
    {



        int current_number_length = number.length();
        //add a zero to the number if its a three digit number
        if (current_number_length < 2 * seed_length) {
            int diff = 2 * seed_length - current_number_length;
            for (int j = 0; j < diff; j++)
            {
                number = "1" + number;
            }

        }

        //remove first and last character
        string middle = number.substr(seed_length / 2, seed_length);

        //stoi converts string to int
        int new_number = pow(stoi(middle), 2);
        number = to_string(new_number);
        output.push_back(stoi(middle));
    }

    //transformiere Zahl zwischen 0 und 1

    return output;
}


