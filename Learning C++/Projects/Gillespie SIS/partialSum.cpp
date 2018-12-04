#include <cstdlib>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>

using namespace std;

vector<double> cumSum5(vector<double> arr)
{
    double curSum = 0;
    vector<double> arr2(5);
    // vector<double> cumVec(5);
    arr2[0] = arr[0];
    for (int i = 1; i < 5; i++)
    {
        arr2[i] = arr2[i-1] + arr[i];
    }
    return arr2;
}

int main() {
    vector<double> myVector = {1,2,3,4,5};
    vector<double> myCumVector(5);
    // partial_sum(myVector.begin(), myVector.end(), myCumVector);
    vector<double> myArrayNew(5);
    myArrayNew = cumSum5(myVector);
    cout << "Hello" << endl;
    for (int i=0; i<5; i++)
    {
        cout << myArrayNew[i] << endl;
    }
    return 0;
}
