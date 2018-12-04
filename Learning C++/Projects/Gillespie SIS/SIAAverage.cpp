#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
// #include "mex.h"

// #include "gnuplot-iostream.h"

using namespace std;

void createHouseMat(int arr[1500][10])
{
    for( int i=0; i<1500; i++ )
     {
        arr[i][0] = 4;
        arr[i][1] = 4;
        arr[i][2] = 0;
        arr[i][3] = 0;
        arr[i][4] = 0;
        arr[i][5] = 0;
        arr[i][6] = 0;
        arr[i][7] = 0;
        arr[i][8] = 0;
        arr[i][9] = 0;
    }
}

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
    int n = 0;
    for (int i=0; i < arr.size(); i++)
    {
        if (arr[i] >= P) {
            return n;
        } else {
            n = n + 1;
        }
    }
    return n;
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
    // vector<double> arr2(N);
    // vector<double> cumVec(1500);
    // double a = arr[0];
    arr2[0] = arr[0];
    for (int i = 1; i < N; i++)
    {
        arr2[i] = arr2[i-1] + arr[i];
    }
}

void interpolate(vector<double> &templateVec, vector<int> &actualVec, vector<double> &timeVec, double tmax)
{
    templateVec[0] = actualVec[0];
    double timeStep = tmax/templateVec.size();
    int n = 0;
    for (int i=1; i<templateVec.size(); i++) {
        double curVal = actualVec[i];
        cout << "timeVec[i] " << timeVec[i] << endl;
        cout << "timeStep*i " << timeStep*i << endl; 
        if (timeVec[i] > timeStep*i)
        {
            templateVec[i] = actualVec[n];
        } else if (timeVec[i] == timeStep*i)
        {
            templateVec[i] = actualVec[i+1];
        } else if (timeVec[i] < timeStep*i)
        {
            int n2 = findGreaterP(timeStep*i, timeVec);
            templateVec[i] = actualVec[n2-1];
        }
    }
}




int main()
{
    // vector<double> myArr = {1,1,1,2,2,3,3,3,5,8,10};
    // printVector(myArr);
    // vector<double> arr2 = cumSum(myArr, myArr.size());
    // printVector(arr2);
    // int testn = findGreaterP(14, arr2);
    // cout << testn << endl;

    // cin.get();
    // clock_t start;
    // double duration;
    vector<double> vec1 = {0, 0, 0, 0, 0};
    vector<int> vec2 = {1,2,5,10,9,8};
    vector<double> vecTime = {0, 0.2, 0.3, 0.8, 1.2, 2};
    double maxtime = 2;
    interpolate(vec1,vec2, vecTime, maxtime);
    printVector(vec1);
    cin.get();

    // start = clock();
    //Setting up random number generation
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    double beta=0.0516;
    double delta1=0.0513;
    double delta2 = 0.0513;
    double rho = 0.0165;
    double lambda = 0.185;
    double tmax = 1000000;
    double eps = 0.004;
    double alpha = 0.164865777630342;
    int initInfNum = 160;
    int numHouses = 1500;
    int numTimes = 1000;

    vector<double> sumInf;

    int n=0;
    // double maxRand = RAND_MAX;

    for (int numSims = 0; numSims < numTimes; numSims++)
    {

        double rand1 = 0;
        double rand2 = 0;
        double dt = 0;
        double P = 0;
        double remainP = 0;
        double R1 = 0;
        double R2 = 0;
        double Rtotal = 0;
        int numPeop = numHouses*4;
        double oldSum = 0;
        double newSum = 0;
        int eventHouseNum = 0;
        int eventNum = 0;
        int numI = initInfNum;
        double t1 = 0;
        double diff = 0;


        // double maxRand = (double) RAND_MAX;

        vector<int> Sarray;
        vector<int> Iarray;
        vector<double> tarray;

        int houseMat[numHouses][10];

        createHouseMat(houseMat);

        int initInf[initInfNum];

        for (int i=0; i<initInfNum; i++) initInf[i] = i;

        for (int i=0; i<initInfNum; i++)
        {
            houseMat[initInf[i]][1] = houseMat[initInf[i]][1] - 1;
            houseMat[initInf[i]][2] = houseMat[initInf[i]][2] + 1;
        }

        vector<int> N(numHouses);
        vector<int> S(numHouses);
        vector<int> I(numHouses);
        vector<int> A(numHouses);

        double rateMat[numHouses][5];
        for( int i=0; i<numHouses; i++ )
        {
            rateMat[i][0] = houseMat[i][4];
            rateMat[i][1] = houseMat[i][5];
            rateMat[i][2] = houseMat[i][6];
            rateMat[i][3] = houseMat[i][7];
            rateMat[i][4] = houseMat[i][8];
        }

        for (int i=0; i <numHouses; i++)
        {
            N[i] = houseMat[i][0];
            S[i] = houseMat[i][1];
            I[i] = houseMat[i][2];
            A[i] = houseMat[i][3];
        }
        Sarray.push_back(numHouses-initInfNum);
        Iarray.push_back(initInfNum);
        tarray.push_back(0);

        eps = alpha*numI/numPeop;

        for (int i=0; i<numHouses; i++)
        {
            rateMat[i][0] = (eps + beta*I[i]/(N[i]-1))*S[i];
            rateMat[i][1] = delta1*I[i];
            rateMat[i][2] = delta2*A[i];
            rateMat[i][3] = rho*A[i];
            rateMat[i][4] = lambda*I[i];
            // cout << rateMat[i][0] << endl;
            // cout << i << endl;
            // cin.get();
        }

        vector<double> newHouseRate(5);
        vector<double> sumMat(numHouses);
        vector<double> cumTotMat(numHouses);
        vector<double> withinCumTotMat(5);
        vector<double> tempStorage(5);

        for (int i = 0; i < numHouses; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                // cout << sumMat[i] << endl;
                // cin.get();
                sumMat[i] = sumMat[i] + rateMat[i][j];
                // cout << sumMat[i] << endl;
            }
        }
        // printVector(sumMat);
        // cin.get();

        // partial_sum(sumMat[0], sumMat[1499], cumTotMat);
        // cumTotMat = cumSum1500(sumMat);
        cumSum(sumMat, numHouses, cumTotMat);
        // cout << "cumTotMat[234]: " << cumTotMat[234] << endl;

        while( t1<tmax ) {
            // cout << t1 << endl;
            if( n > 0 ) {
                // // cout << eventNum << endl;
                if (n > 1) {
                    eps = alpha*numI/numPeop;
                    if (eventNum == 0) {
                        newHouseRate[0] = (eps + beta*I[eventHouseNum]/(N[eventHouseNum]-1))*S[eventHouseNum];
                        newHouseRate[1] = delta1*I[eventHouseNum];
                        newHouseRate[4] = lambda*I[eventHouseNum];
                    } else if (eventNum == 1) {
                        newHouseRate[0] = (eps + beta*I[eventHouseNum]/(N[eventHouseNum]-1))*S[eventHouseNum];
                        newHouseRate[1] = delta1*I[eventHouseNum];
                        newHouseRate[4] = lambda*I[eventHouseNum];
                    } else if (eventNum == 2) {
                        newHouseRate[0] = (eps + beta*I[eventHouseNum]/(N[eventHouseNum]-1))*S[eventHouseNum];
                        newHouseRate[2] = delta2*A[eventHouseNum];
                        newHouseRate[3] = rho*A[eventHouseNum];
                    } else {
                        newHouseRate[0] = (eps + beta*I[eventHouseNum]/(N[eventHouseNum]-1))*S[eventHouseNum];
                        newHouseRate[1] = delta1*I[eventHouseNum];
                        newHouseRate[2] = delta2*A[eventHouseNum];
                        newHouseRate[3] = rho*A[eventHouseNum];
                        newHouseRate[4] = lambda*I[eventHouseNum];
                    }
                }
                    newHouseRate[0] = (eps + beta*I[eventHouseNum]/(N[eventHouseNum]-1))*S[eventHouseNum];
                    newHouseRate[1] = delta1*I[eventHouseNum];
                    newHouseRate[2] = delta2*A[eventHouseNum];
                    newHouseRate[3] = rho*A[eventHouseNum];
                    newHouseRate[4] = lambda*I[eventHouseNum];
                // for (int i = 0; i < 1500; i++)
                // {
                //     cout << "N" << N[i] << endl;
                //     cout << "S" << S[i] << endl;
                //     cout << "I" << I[i] << endl;
                //     cout << "A" << A[i] << endl;
                //     cin.get();
                // }
                newSum = 0;
                oldSum = 0;
                for (int i=0; i<5; i++)
                {
                    newSum = newSum + newHouseRate[i];
                    oldSum = oldSum + rateMat[eventHouseNum][i];
                }
                diff = newSum - oldSum;
                sumMat[eventHouseNum] = sumMat[eventHouseNum] + diff;
                rateMat[eventHouseNum][0] = newHouseRate[0];
                rateMat[eventHouseNum][1] = newHouseRate[1];
                rateMat[eventHouseNum][2] = newHouseRate[2];
                rateMat[eventHouseNum][3] = newHouseRate[3];
                rateMat[eventHouseNum][4] = newHouseRate[4];
                // cumTotMat = cumSum1500(sumMat);
                // cumSum(sumMat, 5000, cumTotMat);

                for (int i = eventHouseNum; i < numHouses; i++)
                {
                    cumTotMat[i] = cumTotMat[i] + diff;
                }
                // cout << "Hello" << endl;
                // cout << cumTotMat[234] << endl;
            }
            // cout << "Hello 1!" << endl;
            // cout << "Rtotal: " << Rtotal << endl;
            Rtotal = cumTotMat.back();
            // cout << "Rtotal: " << Rtotal << endl;
            // cout << Rtotal << endl;
            // cin.get();
            rand1 = dis(gen);
            rand2 = dis(gen);
            // cout << rand() << endl;
            // cout << maxRand << endl;
            // cout << rand()/maxRand << " Is this zero? " << endl;
            // rand1 = rand()/maxRand;
            // rand2 = rand()/maxRand;
            // cout << rand1 << "   " << rand2 << endl;
            // cout << log(2.718281828) << endl;
            // cout << log(10) << endl;
            // cin.get();
            dt = -log(rand1)/Rtotal;
            // cout << "dt: " << dt << endl;
            // cout << dt << endl;
            // cin.get();
            P = rand2*Rtotal;
            eventHouseNum = findGreaterP(P, cumTotMat);
            // cout << "eventhouseNum: " << eventHouseNum << endl;
            // cout << P << endl;
            // cin.get();
            // cout << eventHouseNum << endl;
            if (eventHouseNum == 0) {
                remainP = P;
            } else {
                remainP = P - cumTotMat[eventHouseNum-1];
            }
            for (int i = 0; i<5; i++)
            {
                tempStorage[i] = rateMat[eventHouseNum][i];
            }
            // withinCumTotMat = cumSum5(tempStorage);
            cumSum(tempStorage, 5, withinCumTotMat);
            eventNum = findGreaterP(remainP, withinCumTotMat);

            // for (int i=0; i < withinCumTotMat.size(); i++)
            // {
            //     // cout << "House Number: " << eventHouseNum << endl;
            //     // cout << "ELEMENT i OF WITHINCUMTOT VECTOR: " << withinCumTotMat[i] << endl;
            //     // cout << "i: " << i << endl;
            //     // cin.get();
            // }
            // cout << "eventNum: " << eventNum << endl;
            if (eventNum == 0) {
                S[eventHouseNum] = S[eventHouseNum] - 1;
                I[eventHouseNum] = I[eventHouseNum] + 1;
                numI = numI + 1;
            } else if (eventNum == 1) {
                S[eventHouseNum] = S[eventHouseNum] + 1;
                I[eventHouseNum] = I[eventHouseNum] - 1;
                numI = numI - 1;
            } else if (eventNum == 2) {
                S[eventHouseNum] = S[eventHouseNum] + 1;
                A[eventHouseNum] = A[eventHouseNum] - 1;
            } else if (eventNum == 3) {
                I[eventHouseNum] = I[eventHouseNum] + 1;
                A[eventHouseNum] = A[eventHouseNum] - 1;
                numI = numI + 1;
            } else if (eventNum == 4) {
                A[eventHouseNum] = A[eventHouseNum] + 1;
                I[eventHouseNum] = I[eventHouseNum] - 1;
                numI = numI - 1;
                // cout << "Event 5! " << endl;
            }

            n = n + 1;
            t1 = t1+dt;
            // cout << rand1 << endl;
            // cout << rand2 << endl;
            Iarray.push_back(numI);
            tarray.push_back(t1);
            // cout << "Number of infecteds at time step n: " << numI << endl;
            // cout << dt << endl;
        }
        // cout << I << endl;
        // duration = (clock() - start) / (double) CLOCKS_PER_SEC;
        // Gnuplot gp;
        // gp << "plot" << gp.file1d(tarray) << gp.file1d(Iarray) << endl;
        // cout<< "printf: " << duration << '\n';
        // for (int i=0; i < Iarray.size(); i++)
        // {
        //     cout << Iarray[i] << endl;
        //     cin.get();
        // }
        // double infAv = halfStatMean(Iarray);
        // cout << "The average number of infecteds at steady state is: " << infAv << endl;
        // cout <<  "Total Number Infected at steady state: " << numI << endl;
        // cout<< "printf: " << duration << '\n';
    }
    return 0;
}
