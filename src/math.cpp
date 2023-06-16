#include <iostream>
#include <string>
#include <cmath>

using namespace std;

namespace math {
    int factorial(int x) {
        int ans = 1;
        if(x == 0) return 1;
        if(x > 0) {
            for(int i = 1; i <= x; i++) {
                ans = ans * i;
            }
            return ans;
        } else {
            for(int i = 1; i <= x; i++) {
                ans = ans * i;
            }
            return ans*-1;
        }
    }

    double newton(double x_n, auto f, auto f1) {
        double x_n1 = x_n;
        double x_absChange = 9999;
        double accuracy = pow(1, -10);

        while(x_absChange > accuracy) {
            x_n = x_n1;

            x_n1 = x_n - (f(x_n)/f1(x_n));

            x_absChange = abs(x_n1 - x_n);
        }

        return x_n1;
    }
}