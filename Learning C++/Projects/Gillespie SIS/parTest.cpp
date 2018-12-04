#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include <omp.h>

using namespace std;
int main()
{
    clock_t start;
    double duration;
    vector<double> vect(1000000000);
    #pragma omp parallel for
    for (int i=0; i<1000000000; i++)
         // int N = 100000;
        {
            double v1 = rand() % 100;
            double v2 = rand() % 100;
            double v3 = rand() % 100;
            double v4 = rand() % 100;
            double v5 = rand() % 100;
            double v6 = rand() % 100;
            double v7 = rand() % 100;
            double v8 = rand() % 100;
            double v9 = rand() % 100;
            double v10 = rand() % 100;
            double v11 = rand() % 100;
            double v12 = rand() % 100;
            double v13 = rand() % 100;
            double v14 = rand() % 100;
            double v15 = rand() % 100;
            double v16 = rand() % 100;
            double vv1 = ((v1*v2+1)/(v3*v4+1));
            double vv2 = ((v5*v6+1)/(v7*v8+1));
            double vv3 = ((v9*v10+1)/(v11*v12+1));
            double vv4 = ((v13*v14+1)/(v15*v16+1));
            vect[i] = ((vv1*vv2+1)/(vv3*vv4+1));
            // cout << vect[i] << endl;
        }

         // printf("https://helloacm.com\n");
         // int id = omp_get_thread_num();
         // int data = id;
         // int total = omp_get_num_threads();
         // printf("Greetings from process %d out of %d with Data %d\n", id, total, data);
     // printf("parallel for ends.\n");
    // cout << vect[1000] << endl;
    duration = (clock() - start) / (double) CLOCKS_PER_SEC;
    cout<< "printf: " << duration << '\n';
    return 0;
}
