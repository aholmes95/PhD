#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include <omp.h>
#include <fstream>
#include <string>
#include <iomanip>
#include "boost/random.hpp"
// #include "mex.h"

// #include "gnuplot-iostream.h"

using namespace std;

double mean(vector<double> arr)
{
    double length = arr.size();
    double sum = 0;
    for (int i =0; i < length; i++)
    {
        sum = sum + arr[i];
    }
    double av = sum/length;
    return av;
}

double halfStatMean(vector<int> &arr)
{
    double length = arr.size();
    double newLength = round(length/2);
    vector<double> arr2(newLength);
    for (int i=0; i < newLength; i++)
    {
        arr2[i] = arr[newLength-i];
    }
    double newAv = mean(arr2);
    return newAv;
}

int findGreaterP(double P, vector<double> &arr)
{
    if (arr[0] >= P) {
        return 0;
    }
    else {
        for (int i=1; i < arr.size(); i++)
        {
            if (arr[i] >= P && arr[i-1] < P) {
                return i;
            }
        }
    }
}

void printVector(vector<double> &myArr)
{
    cout << "| " << myArr[0];
    for (int i = 1; i < myArr.size() - 1; i++)
    {
        cout << " | " <<myArr[i];
    }
    cout << " | " << myArr.back() << " |" << endl;
}

void cumSum(vector<double> &arr, int N, vector<double> &arr2)
{
    arr2[0] = arr[0];
    for (int i = 1; i < N; i++)
    {
        arr2[i] = arr2[i-1] + arr[i];
    }
}

int sum(vector<int> arr)
{
    double sum = 0;
    for (int i = 0; i < arr.size(); i++)
    {
        sum = sum + arr[i];
    }
    return sum;
}


double SIATreatment(double alpha, double beta, double gamma, int numHouseholds, int householdSize, int index)
{
    cout << setprecision(16);
    // cin.get();
    clock_t start;
    double duration;

    start = clock();

    //Setting up random number generation
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    double tmax = 700;
    double eps = 0;

    int initInfNum = numHouseholds/10;
    int numHouses = numHouseholds;

    int n=0;
    double rand1 = 0;
    double rand2 = 0;
    double dt = 0;
    double P = 0;
    double remainP = 0;
    double R1 = 0;
    double R2 = 0;
    double Rtotal = 0;
    double oldSum = 0;
    double newSum = 0;
    int eventHouseNum = 0;
    int eventNum = 0;
    double t1 = 0;
    double diff = 0;
    int flag = 0;

    vector<int> Sarray;
    vector<int> Iarray;
    vector<double> tarray;

    vector<int> N(numHouses);
    vector<int> S(numHouses);
    vector<int> I(numHouses);

    for (int i = 0; i < numHouses; i++)
    {
        N[i] = householdSize;
        S[i] = householdSize;
        I[i] = 0;
    }

    int numPeop = sum(N);
    int numS = numPeop - initInfNum;
    int numI = initInfNum;
    for (int i = 0; i < initInfNum; i++)
    {
        S[i] = S[i] - 1;
        I[i] = I[i] + 1;
    }

    double rateMat[numHouses][2];

    Sarray.push_back(numPeop-initInfNum);
    Iarray.push_back(initInfNum);
    tarray.push_back(0);

    vector<double> sumMat(numHouses);
    vector<double> cumTotMat(numHouses);
    vector<double> withinCumTotMat(2);
    vector<double> tempStorage(2);

    while( t1<tmax ) {
        // cout << t1 << endl;
        eps = alpha*sum(I)/sum(N);
        for (int i = 0; i < numHouses; i++)
        {
            if (N[i] == 1)
            {
                rateMat[i][0] = (eps + beta*I[i])*S[i];
            }
            else
            {
                rateMat[i][0] = (eps+beta*I[i]/(N[i]-1))*S[i];
            }
            rateMat[i][1] = gamma*I[i];
        }

        double tempSum2 = 0;
        for (int i = 0; i < numHouses; i++)
        {
            tempSum2 = 0;
            for (int j = 0; j < 2; j++)
            {
                tempSum2 = tempSum2 + rateMat[i][j];
            }
            sumMat[i] = tempSum2;
        }
        cumSum(sumMat, numHouses, cumTotMat);
        Rtotal = cumTotMat.back();
        rand1 = dis(gen);
        rand2 = dis(gen);
        dt = -log(rand1)/Rtotal;

        P = rand2*Rtotal;
        eventHouseNum = findGreaterP(P, cumTotMat);

        if (eventHouseNum == 0) {
            remainP = P;
        } else {
            remainP = P - cumTotMat[eventHouseNum-1];
        }

        for (int i = 0; i<2; i++)
        {
            tempStorage[i] = rateMat[eventHouseNum][i];
        }

        cumSum(tempStorage, 2, withinCumTotMat);
        eventNum = findGreaterP(remainP, withinCumTotMat);


        if (eventNum == 0) {
            S[eventHouseNum] = S[eventHouseNum] - 1;
            I[eventHouseNum] = I[eventHouseNum] + 1;
        } else if (eventNum == 1) {
            S[eventHouseNum] = S[eventHouseNum] + 1;
            I[eventHouseNum] = I[eventHouseNum] - 1;
        }

        numI = sum(I);
        numS = sum(S);

        if (numI == 0)
        {
	    cout << "Disease died out" << endl;
            SIATreatment(alpha, beta, gamma, numHouseholds, householdSize, index);
            return 0;
        }
        n = n + 1;
        t1 = t1+dt;
        Iarray.push_back(numI);
        Sarray.push_back(numS);
        tarray.push_back(t1);
    }
    ofstream outFile("IFiles2/" + to_string(index) + "myFileI.csv");
    ofstream outFile2("tFiles2/" + to_string(index) + "myFilet.csv");
    cout << "Should be saving" << endl;
    for (const auto &e : Iarray) outFile << setprecision(16) << e << "\n";
    for (const auto &e : tarray) outFile2 << setprecision(16) << e << "\n";
    return 0;

}

int main()
{

    int numRuns = 1;
    int numHouseholds = 5000;
    int householdSize = 2;
    double alpha = 0.1;
    double beta = 0.9;
    double gamma = 0.328304;
    int temp = 0;

    #pragma omp parallel for
    for (int i = 0; i < numRuns; i++)
    {
        cout << i << endl;
        temp = SIATreatment(alpha, beta, gamma, numHouseholds, householdSize, i);
    }
    return 0;
}
