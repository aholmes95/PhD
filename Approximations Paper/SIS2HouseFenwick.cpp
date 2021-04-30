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

double getSum(vector<double> BITree, int index)
{
    double outSum = 0;
    index = index + 1;
    while( index > 0 )
    {
        outSum = outSum + BITree[index];
        index -= (index & (-index));
    }
    return outSum;
}

void updateBIT(vector<double> &BITree, int i, double v)
{
    i = i + 1;
    int n = BITree.size() - 1;
    while( i <= n )
    {
        BITree[i] = BITree[i] + v;
        i += i & (-i);
    }
}

vector<double> construct(vector<double> arr, int n)
{
    vector<double> BITree(n+1);
    for( int i = 0; i < n; i++ )
    {
        updateBIT(BITree, i, arr[i]);
    }
    return BITree;
}

int set_bit(int value, int bit)
{
    return value | (1<<bit);
}

int clear_bit(int value, int bit)
{
    return value & ~(1<<bit);
}

int search(vector<double> BITree, double x, int n)
{
    int i = 0;
    int numBits = ceil(log2(n));
    for( int b = numBits; b > 0; b--)
    {
        i = set_bit(i, b-1);
        if( i < n )
        {
            if( BITree[i] < x )
            {
                x -= BITree[i];
            } else
            {
                i = clear_bit(i, b-1);
            }
        } else
        {
            i = clear_bit(i, b-1);
        }
    }
    return i;
}


double SIATreatment(double eps, double beta, double gamma, int numHouseholds, int householdSize, int index)
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

    double tmax = 500;

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

    for (int i = 0; i < numHouses; i++)
    {
        rateMat[i][0] = (eps+beta*I[i])*S[i];
        rateMat[i][1] = gamma*I[i];
    }

    Sarray.push_back(numPeop-initInfNum);
    Iarray.push_back(initInfNum);
    tarray.push_back(0);

    vector<double> sumMat(numHouses);

    double tempSum2 = 0;
    for (int i = 0; i < numHouses; i++) //Calculate vector containing sum of rates for each house
    {
        tempSum2 = 0;
        for (int j = 0; j < 2; j++)
        {
            tempSum2 = tempSum2 + rateMat[i][j];
        }
        sumMat[i] = tempSum2;
    }

    vector<double> BITree = construct(sumMat, numHouses);

    while( t1<tmax ) {
        // cout << t1 << endl;
        Rtotal =getSum(BITree, numHouses-1);
        rand1 = dis(gen);
        rand2 = dis(gen);
        dt = -log(rand1)/Rtotal;
        P = rand2*Rtotal;
        eventHouseNum = search(BITree, P, numHouses);

        if (eventHouseNum == 0)
        {
            remainP = P;
        }
        else
        {
            remainP = P - getSum(BITree, eventHouseNum-1);
        }

        if (remainP <= rateMat[eventHouseNum][0])
        {
            eventNum = 0;
        }
        else
        {
            eventNum = 1;
        }


        if (eventNum == 0) {
            S[eventHouseNum] = S[eventHouseNum] - 1;
            I[eventHouseNum] = I[eventHouseNum] + 1;
            numI = numI+1;
        } else if (eventNum == 1) {
            S[eventHouseNum] = S[eventHouseNum] + 1;
            I[eventHouseNum] = I[eventHouseNum] - 1;
            numI = numI-1;
        }

        if (numI == 0)
        {
            cin.get();
        }

        double newRate0 = (eps+beta*I[eventHouseNum])*S[eventHouseNum];
        double newRate1 = gamma*I[eventHouseNum];
        double totalChange = newRate0 + newRate1 - rateMat[eventHouseNum][0] - rateMat[eventHouseNum][1];
        rateMat[eventHouseNum][0] = newRate0;
        rateMat[eventHouseNum][1] = newRate1;
        updateBIT(BITree, eventHouseNum, totalChange);
        n = n + 1;
        t1 = t1+dt;
        Iarray.push_back(numI);
        tarray.push_back(t1);
    }
    ofstream outFile("IFiles/" + to_string(index) + "myFileI.csv");
    ofstream outFile2("tFiles/" + to_string(index) + "myFilet.csv");
    cout << "Should be saving" << endl;
    for (const auto &e : Iarray) outFile << setprecision(16) << e << "\n";
    for (const auto &e : tarray) outFile2 << setprecision(16) << e << "\n";
    return 0;

}

int main()
{

    int numRuns = 500;
    int numHouseholds = 1500;
    int householdSize = 2;
    double eps = 0.02;
    double beta = 0.06;
    double gamma = 0.065;
    int temp = 0;

    #pragma omp parallel for
    for (int i = 0; i < numRuns; i++)
    {
        cout << i << endl;
        temp = SIATreatment(eps, beta, gamma, numHouseholds, householdSize, i);
    }
    return 0;
}
