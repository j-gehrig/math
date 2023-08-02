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

    double sin_x_reduce(double x) {
        if(x <= (M_PI / 2)) {
        } else if(x <= M_PI) {
            x = M_PI - x;
        } else if(x <= (1.5*M_PI)) {
            x = M_PI - x;
        } else if(x <= (2*M_PI)) {
            x = (-2*M_PI) + x;
        }
        return x;
    }
    double sin_taylor(double x) {
        double n = 0;
        double e = 1e-17;
        double sa;
        double sb = 9999;
        double s = 0;
        while(sb > e) {
            sa = s;

            s += pow(-1, n) * (pow(x, 2*n + 1) / factorial(2*n + 1));
            sb = s - sa;
            if(sb < 0.0) sb = -1.0*sb;

            n++;
        }
        return s;
    }
    double cos(double x) {
        double u = floor(x / (M_PI * 2));
        double ans;

        if(u < -1) {
            u = u * -1;
            x = (u - floor(u)) * (M_PI * 2) * -1;
        } else if(u > 1) {
            x = (u - floor(u)) * (M_PI * 2);
        }
        if(abs(x) <= M_PI/2 && abs(x) >= M_PI/4) {
            ans = sin(M_PI/2 - x);
        } else if(x >= 0) {
            ans = sin(M_PI/2 + x);
        } else {
            ans = sin(M_PI/2 + abs(x)) * -1;
        }
        return ans;
    }
    double cos_taylor(double x) {
        double n = 0;
        double e = 1e-17;
        double sa;
        double sb = 9999;
        double s = 0;
        while(sb > e) {
            sa = s;

            s += pow(-1, n) * (pow(x, 2*n) / factorial(2*n));
            sb = s - sa;
            if(sb < 0.0) sb = -1.0*sb;

            n++;
        }
        return s;
    }
    double sin(double x) {
        double u = floor(x / (M_PI * 2));
        double ans;

        if(u < -1) {
            u = u * -1;
            x = (u - floor(u)) * (M_PI * 2) * -1;
        } else if(u > 1) {
            x = (u - floor(u)) * (M_PI * 2);
        }
        if(x < 0) {
            x = sin_x_reduce(x * -1);
            if(x >= M_PI/4 && x < M_PI/2) {
                ans = cos_taylor(M_PI/2 - x) * -1;
            } else {
                ans = sin_taylor(x) * -1;
            }
        } else {
            x = sin_x_reduce(x);
            if(x >= M_PI/4 && x < M_PI/2) {
                ans = cos_taylor(M_PI/2 - x);
            } else {
                ans = sin_taylor(x);
            }
        }

        return ans;
    }
    double tan(double x) {
        return sin(x)/cos(x);
    }


    template<typename T>
    double newton(double x_n, double (*f)(T, T, T), double (*f1)(T, T, T), double param = 0, double param1 = 0) {
        double x_n1 = x_n;
        double x_absChange = 9999;
        double accuracy = pow(1, -10);

        while(x_absChange > accuracy) {
            x_n = x_n1;

            x_n1 = x_n - (f(x_n, param, param1)/f1(x_n, param, param1));

            x_absChange = abs(x_n1 - x_n);
        }

        return x_n1;
    }

    double root_f(double x, double n, double a) {
        return (pow(x, n) - a);
    }
    double root_f1(double x, double n, double param1) {
        return (n * pow(x, n-1));
    }
    double root(double n, double a) {
        if(a == 0) return 0;
        if(n == 1) return a;
        if(n == 0) return 1;

        double x_0 = n+a;
        return newton(x_0, root_f, root_f1, n, a);
    }

    double sqroot_Heron(double a, double e) {
        double x, xa, xb;
        x = 1; xb = 9999.0;
        while(xb > e) {
            xa = x;

            x = 0.5 * (x + a/x);

            xb = x - xa;
            if(xb < 0) { xb = -1*xb; }
        }
        return x;
    }

    double arcsin_f(double x, double a, double param1) {
        return (sin(x) - a);
    }
    double arcsin_f1(double x, double param, double param1) {
        return cos(x);
    }
    double arcsin(double x) {
        if(x > 1) return 0;

        x = newton(0, arcsin_f, arcsin_f1, x);
        return x;
    }
    double arccos(double x) {
        x = arcsin(x) - M_PI/2;
        return x;
    }

    double arctan_f(double x, double a, double param1) {
        return (tan(x) - a);
    }
    double arctan_f1(double x, double param, double param1) {
        return (1/pow(cos(x),2));
    }
    double arctan(double x) {
        if(x > 1) return 0;

        x = newton(0, arctan_f, arctan_f1, x);
        return x;
    }

    double exp_taylor(double x) {
        double e = 1e-15, n = 1.0, s = 1.0, sa, sb = 9999.0;
        while(sb > e) {
            sa = s;

            s = s + pow(x, n) / factorial(n);

            sb = s - sa;
            if(sb < 0.0) sb = -1.0*sb;

            n++;
        }
        return s;
    }
    double exp(double x) {
        double x_dez, x_betrag, x_ger;
        if(x < 0) {
            x_betrag = x * -1;
            x_dez = -1 * (x_betrag - floor(x_betrag));
            x_ger = -1 * floor(x_betrag);
        } else {
            x_dez = x - floor(x);
            x_ger = floor(x);
        }
        return exp_taylor(x_ger) * exp_taylor(x_dez);
    }

    double ln_f(double x, double a, double param1) {
        return (exp(x) - a);
    }
    double ln_f1(double x, double param, double param1) {
        return exp(x);
    }
    double ln(double x) {
        double x_0 = -1 + sqroot_Heron(2*x-1, 1e-5);
        return newton(x_0, ln_f, ln_f1, x);
    }
    double log(double b, double x) {
        return ln(x)/ln(b);
    }

    double pow(double a, double n) {
        double ans = 1;

        float n2;
        if(n == 0) {
            return 1;
        } else if(n == 1) {
            return a;
        } else if(a == 0) {
            return 0;
        } else if(a == 1) {
            return 1;
        } else if(modf(n, &n2) == 0) { // if exponent has floating point
            if(a > 0) {
                ans = exp(ln(a)*n);
            } else {
                ans = -1*exp(ln(-1*a)*n);
            }
        } else if(n > 0) {
            for(int i = 0; i < n; i++) {
                ans = ans * a;
            }
        } else {
            for(int i = 0; i > n; i--) {
                ans = ans / a;
            }
        }

        return ans;
    }
}
