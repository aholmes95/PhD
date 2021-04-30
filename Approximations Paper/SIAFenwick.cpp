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
#include <math.h>
#include <stdio.h>
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

double SIAFenwick(int index)
{
    // cin.get();
    clock_t start;
    double duration;

    start = clock();

    //Setting up random number generation
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    double beta=0.0516;
    double delta1=0.0513;
    double delta2 = 0.0513;
    double rho = 0.0165;
    double lambda = 0.185;

    // double beta=0.15;
    // double delta1=0.1;
    // double delta2 = 0;
    // double rho = 0;
    // double lambda = 0;

    double tmax = 700;
    double eps = 0.004;
    int initInfNum = 160;
    int numHouses = 10;
    int householdSize = 600;

    int n=0;
    double rand1 = 0;
    double rand2 = 0;
    double dt = 0;
    double P = 0;
    double remainP = 0;
    double R1 = 0;
    double R2 = 0;
    double Rtotal = 0;
    int numPeop = numHouses*householdSize;
    int eventHouseNum = 0;
    int eventNum = 0;
    int numS = numPeop - initInfNum;
    int numI = initInfNum;
    int numA = 0;
    double t1 = 0;
    double diff = 0;

    vector<int> Sarray;
    vector<int> Iarray;
    vector<double> tarray;

    vector<int> N(numHouses);
    vector<int> S(numHouses);
    vector<int> I(numHouses);
    vector<int> A(numHouses);

    for (int i = 0; i < numHouses; i++)
    {
        N[i] = householdSize;
        S[i] = householdSize;
        I[i] = 0;
        A[i] = 0;
    }

    for (int i = 0; i < numHouses; i++)
    {
        S[i] = S[i] - 16;
        I[i] = I[i] + 16;
    }

    double rateMat[numHouses][5];

    Sarray.push_back(numPeop-initInfNum);
    Iarray.push_back(initInfNum);
    tarray.push_back(0);

    vector<double> sumMat(numHouses);
    vector<double> cumTotMat(numHouses);
    vector<double> withinCumTotMat(5);
    vector<double> tempStorage(5);

    for (int i = 0; i < numHouses; i++)
    {
        rateMat[i][0] = (eps+beta*I[i]/(N[i]-1))*S[i];
        rateMat[i][1] = delta1*I[i];
        rateMat[i][2] = delta2*A[i];
        rateMat[i][3] = rho*A[i];
        rateMat[i][4] = lambda*I[i];
    }



    double tempSum2 = 0;
    for (int i = 0; i < numHouses; i++) //Calculate vector containing sum of rates for each house
    {
        tempSum2 = 0;
        for (int j = 0; j < 5; j++)
        {
            tempSum2 = tempSum2 + rateMat[i][j];
        }
        sumMat[i] = tempSum2;
    }

    vector<double> BITree = construct(sumMat, numHouses);

    while( t1<tmax ) {

        Rtotal = getSum(BITree, numHouses-1);
        if(Rtotal < 0)
        {
            cin.get();
        }
        rand1 = dis(gen);
        rand2 = dis(gen);
        dt = -log(rand1)/Rtotal; //time until next event
        P = rand2*Rtotal;
        eventHouseNum = search(BITree, P, numHouses); //Determine which house the event happens in
        if (eventHouseNum == 0) {
            remainP = P;
        } else {
            remainP = P - getSum(BITree, eventHouseNum-1);
        }

        if (remainP <= rateMat[eventHouseNum][0])
        {
            eventNum = 0;
        } else if (remainP <= rateMat[eventHouseNum][0] + rateMat[eventHouseNum][1])
        {
            eventNum = 1;
        } else if (remainP <= rateMat[eventHouseNum][0] + rateMat[eventHouseNum][1] + rateMat[eventHouseNum][2])
        {
            eventNum = 2;
        } else if (remainP <= rateMat[eventHouseNum][0] + rateMat[eventHouseNum][1] + rateMat[eventHouseNum][2] + rateMat[eventHouseNum][3])
        {
            eventNum= 3;
        } else
        {
            eventNum = 4;
        }

        // Update the state of the house depending on which event it was that happened.
        // Also update counter containing the total number of I and A to stop us having to re-sum it every time.
        if (eventNum == 0) {
            S[eventHouseNum] = S[eventHouseNum] - 1;
            I[eventHouseNum] = I[eventHouseNum] + 1;
            numI = numI + 1;
            numS = numS - 1;
        } else if (eventNum == 1) {
            S[eventHouseNum] = S[eventHouseNum] + 1;
            I[eventHouseNum] = I[eventHouseNum] - 1;
            numI = numI - 1;
            numS = numS + 1;
        } else if (eventNum == 2) {
            S[eventHouseNum] = S[eventHouseNum] + 1;
            A[eventHouseNum] = A[eventHouseNum] - 1;
            numA = numA - 1;
            numS = numS + 1;
        } else if (eventNum == 3) {
            I[eventHouseNum] = I[eventHouseNum] + 1;
            A[eventHouseNum] = A[eventHouseNum] - 1;
            numI = numI + 1;
            numA = numA - 1;
        } else if (eventNum == 4) {
            A[eventHouseNum] = A[eventHouseNum] + 1;
            I[eventHouseNum] = I[eventHouseNum] - 1;
            numI = numI - 1;
            numA = numA + 1;
        }
        // Check for eradication of the disease
        // if (numA + numI == 0)
        // {
        //     cout << "There are no infected (I or A) individuals left. The disease has been eradicated." << endl;
        //     // return 0;
        // }

        double newRate0 = (eps+beta*I[eventHouseNum]/(N[eventHouseNum]-1))*S[eventHouseNum];
        double newRate1 = delta1*I[eventHouseNum];
        double newRate2 = delta2*A[eventHouseNum];
        double newRate3 = rho*A[eventHouseNum];
        double newRate4 = lambda*I[eventHouseNum];
        double totalChange = newRate0 + newRate1 + newRate2 + newRate3 + newRate4 - rateMat[eventHouseNum][0] - rateMat[eventHouseNum][1]- rateMat[eventHouseNum][2] - rateMat[eventHouseNum][3] - rateMat[eventHouseNum][4];
        rateMat[eventHouseNum][0] = newRate0;
        rateMat[eventHouseNum][1] = newRate1;
        rateMat[eventHouseNum][2] = newRate2;
        rateMat[eventHouseNum][3] = newRate3;
        rateMat[eventHouseNum][4] = newRate4;
        updateBIT(BITree, eventHouseNum, totalChange);
        n = n + 1;
        t1 = t1+dt;
        Iarray.push_back(numI);
        tarray.push_back(t1);
    }

    double infAv = halfStatMean(Iarray);

    // Save files containing simulation results to plot in matlab/ python
    ofstream outFile("IFilesFenwickSize2/" + to_string(index) + "myFileIFenwick.csv");
    ofstream outFile2("tFilesFenwickSize2/" + to_string(index) + "myFiletFenwick.csv");
    for (const auto &e : Iarray) outFile << e << "\n";
    for (const auto &e : tarray) outFile2 << e << "\n";

    return infAv; //Return average infectious steady state

}

int main()
{
    int numRuns = 500;
    vector<double> myVec(numRuns);

    // Run in parallel
    #pragma omp parallel for
    for (int i = 0; i < numRuns; i++)
    {
        cout << "i = " << i << endl;
        myVec[i] = SIAFenwick(i);
        // cout << myVec[i] << endl;
    }
    // double relError = myMean/expState - 1;
    // cout << "The average endemic state for the standard model is: " << myMean << endl;
    // cout << "The relative error for this is: " << relError*100 << endl;
    return 0;
}
