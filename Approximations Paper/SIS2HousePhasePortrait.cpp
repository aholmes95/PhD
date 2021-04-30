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
    return arr.size();
}

void printVector(vector<double> &myArr)
{
    cout << "| " << myArr[0];
    for (int i = 1; i < myArr.size() - 1; i++)
    {
        cout << " | " <<myArr[i];
    }
    cout << " | " << myArr.back() << " |" << endl;
    return;
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


double SIATreatment(double alpha, double beta, double gamma, int numHouseholds, int householdSize, int index, int init1, int init2, int init3)
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

    double tmax = 10000;
    double eps = 0;

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

    vector<int> state1Arr;
    vector<int> state2Arr;
    vector<int> state3Arr;

    int state1 = 0;
    int state2 = 0;
    int state3 = 0;

    for (int i = 0; i < numHouses; i++)
    {
        N[i] = householdSize;
        S[i] = householdSize;
        I[i] = 0;
    }

    int numPeop = sum(N);
    int numS = 0;
    int numI = 0;
    for (int i = 0; i < init1; i++)
    {
        S[i] = 2;
        I[i] = 0;
    }

    for (int i = init1+1; i < init1+init2; i++)
    {
        S[i] = 1;
        I[i] = 1;
    }

    for (int i = init1+init2+1; i < init1+init2+init3; i++)
    {
        S[i] = 0;
        I[i] = 2;
    }
    state1 = init1;
    state2 = init2;
    state3 = init3;

    double rateMat[numHouses][2];

    Sarray.push_back(2*state1+state2);
    Iarray.push_back(state2+2*state3);
    tarray.push_back(0);
    state1Arr.push_back(state1);
    state2Arr.push_back(state2);
    state3Arr.push_back(state3);

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
            if (S[eventHouseNum] == 1)
            {
                state1 = state1 - 1;
                state2 = state2 + 1;
            } else if (S[eventHouseNum] == 0)
            {
                state2 = state2 - 1;
                state3 = state3 + 1;
            } else
            {
                cout << "Something has gone wrong around line 240ish" << endl;
                cin.get();
            }
        } else if (eventNum == 1) {
            S[eventHouseNum] = S[eventHouseNum] + 1;
            I[eventHouseNum] = I[eventHouseNum] - 1;
            if (S[eventHouseNum] == 2)
            {
                state1 = state1 + 1;
                state2 = state2 - 1;
            } else if (S[eventHouseNum] == 1)
            {
                state2 = state2 + 1;
                state3 = state3 - 1;
            } else
            {
                cout << "Something has gone wrong around line 250ish" << endl;
                cin.get();
            }
        }

        numI = sum(I);
        numS = sum(S);

        if (numI == 0)
        {
	    cout << "Disease died out" << endl;
            SIATreatment(alpha, beta, gamma, numHouseholds, householdSize, index, init1, init2, init3);
            return 0;
        }
        n = n + 1;
        t1 = t1+dt;
        Iarray.push_back(numI);
        Sarray.push_back(numS);
        state1Arr.push_back(state1);
        state2Arr.push_back(state2);
        state3Arr.push_back(state3);
        tarray.push_back(t1);
    }
    ofstream outFile("IFiles2/" + to_string(index) + "myFileI.csv");
    ofstream outFile2("tFiles2/" + to_string(index) + "myFilet.csv");
    ofstream outFile3("State Files/" + to_string(index) + "state1.csv");
    ofstream outFile4("State Files/" + to_string(index) + "state2.csv");
    ofstream outFile5("State Files/" + to_string(index) + "state3.csv");
    cout << "Should be saving" << endl;
    for (const auto &e : Iarray) outFile << setprecision(16) << e << "\n";
    for (const auto &e : tarray) outFile2 << setprecision(16) << e << "\n";
    for (const auto &e : state1Arr) outFile3 << setprecision(16) << e << "\n";
    for (const auto &e : state2Arr) outFile4 << setprecision(16) << e << "\n";
    for (const auto &e : state3Arr) outFile5 << setprecision(16) << e << "\n";

    return 0;

}

int main(int argc, char* argv[])
{
    double alpha = atof(argv[1]);
    double beta = atof(argv[2]);
    double gamma = atof(argv[3]);
    int N = atof(argv[4]);
    vector<int> initState1 = {100,200,100,200,300,50,30,150};
    vector<int> initState2 = {100,100,50,10,100,200,100,100};
    vector<int> initState3 = {N-200,N-300,N-150,N-210,N-400,N-250,N-130,N-250};
//     int numRuns = initState1.size();
    int numRuns = 1;
    int numHouseholds = N;
    int householdSize = 2;
    int temp = 0;
    int init1 = 0;
    int init2 = 0;
    int init3 = 0;
    double progress = 0;
    int counter = 0;
//     #pragma omp parallel for
    for (int i = 0; i < numRuns; i++)
    {
        progress = counter/numRuns;
//         cout << "progress = " << progress << "%" << endl;
        init1 = initState1[i];
        init2 = initState2[i];
        init3 = initState3[i];
        temp = SIATreatment(alpha, beta, gamma, numHouseholds, householdSize, i, init1, init2, init3);
        counter = counter+1;
    }
    return 0;
}
