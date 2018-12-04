#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
// #include "mex.h"

// #include "gnuplot-iostream.h"

using namespace std;

int main()
{
    clock_t start;
    double duration;

    start = clock();
    //Setting up random number generation
//     random_device rd;
//     mt19937 gen(rd());
//     uniform_real_distribution<> dis(0.0, 1.0);

    int N=5000000;
    int I=50;
    int S = N-I;
    double b=0.2;
    double g=0.15;
    double t = 0;
    double tmax = 700;

    double rand1 = 0;
    double rand2 = 0;
    double dt = 0;
    double P = 0;
    double R1 = 0;
    double R2 = 0;
    double Rtotal = 0;

    double maxRand = (double) RAND_MAX;

    vector<int> Sarray;
    vector<int> Iarray;
    vector<double> tarray;

    Sarray.push_back(S);
    Iarray.push_back(I);
    tarray.push_back(0);

    while( t<tmax ) {
//         cout << t << endl;
        R1 = b*S*I/N;
        R2 = g*I;
        Rtotal = R1 + R2;
        rand1 = rand()/ maxRand;
        rand2 = rand()/ maxRand;
//         rand1 = dis(gen);
//         rand2 = dis(gen);
        dt = -log (rand1)/Rtotal;
        // double dt = 0.0001;
        // double rand2 = 0.3;
        P = rand2*Rtotal;
        if ( P <= R1 ) {
            S = S - 1;
            I = I + 1;
        } else {
            S = S + 1;
            I = I - 1;
        }
        Sarray.push_back(S);
        Iarray.push_back(I);
        t = t + dt;
        tarray.push_back(t);
    }
    cout << I << endl;
    duration = (clock() - start) / (double) CLOCKS_PER_SEC;
    // Gnuplot gp;
    // gp << "plot" << gp.file1d(tarray) << gp.file1d(Iarray) << endl;
    cout<< "printf: " << duration << '\n';
    return 0;
}
