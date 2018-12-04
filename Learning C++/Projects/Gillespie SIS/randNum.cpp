// #include <cstdlib>
// #include <cmath>
#include <iostream>
#include <random>

// using namespace std;

int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int n = 0; n < 10; n++) {
        std::cout << dis(gen) << ' ';
    }
    // double randNum = rand()/ (double) RAND_MAX;
    // double randNum2 = rand()/ (double) RAND_MAX;
    // cout << randNum << endl;
    // cout << randNum2 << endl;
    // cout << RAND_MAX << endl;

    return 0;
}
