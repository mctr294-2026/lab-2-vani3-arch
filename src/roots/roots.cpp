#include <iostream> 
#include "roots.hpp"
#include <functional>
#include <cmath>
using namespace std;

double tol = 1e-9;

/* Tries to find a zero crossing in f() in the interval [a,b] with the bisection method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
 */
bool bisection(function<double(double)> f,
               double a, double b,
               double *root) {
    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) < tol) {
        *root = a;
        return true;
    }

    if (std::abs(fb) < tol) {
        *root = b;
        return true;
    }

    // Root must be bracketed
    if (fa * fb > 0)
        return false;

    double c = (a + b) / 2.0;

    // Bisection converges because interval width halves every step
    for(int i = 0; abs(f(c))> tol; i++) {
        double fc = f(c);

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }

        c = (a + b) / 2.0;
    }
            

    *root = c;
    return true;
}
/* Tries to find a zero crossing in f() in the interval [a,b] with the
 * false positive / regula falsi method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
 */
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root){
    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) < tol) {
        *root = a;
        return true;
    }

    if (std::abs(fb) < tol) {
        *root = b;
        return true;
    }

    // Root must be bracketed
    if (fa * fb > 0)
        return false;

    double c, fc;

    while (true) {
        // False position formula
        c = (a * fb - b * fa) / (fb - fa);
        fc = f(c);

        if (std::abs(fc) < tol) {
            *root = c;
            return true;
        }

        // Keep the sign change
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
}

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the netwon-raphson method, given a function that computes the
 * derivative g() and a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or the derivative is zero.
 */
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    double x = c;
    const int max_iter = 1000;

    if (std::abs(f(x)) < tol) {
        *root = x;
        return true;
    }

    for (int i = 0; i < max_iter; ++i) //Loop up to 1000 times to try to get closer to the root
    {
        double gx = g(x); //derivative at current guess

        if (std::abs(gx) < 1e-14)
            return false;

        double x_new = x - f(x) / gx;

        if (x_new < a || x_new > b)
            return false;

        if (std::abs(x_new - x) < tol)
        {
            *root = x_new;
            return true;
        }

        x = x_new;
    }

    return false;
}

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the secant method, given a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or derivative is zero.
 */
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root){
    const double tol = 1e-10;  

    double x0 = c;
    double x1 = b;

    double f0 = f(x0);
    double f1 = f(x1);

    if (std::abs(f0) < tol) {
        *root = x0;
        return true;
    }

    if (std::abs(f1) < tol) {
        *root = x1;
        return true;
    }

    while (true) {
        if (std::abs(f1 - f0) < 1e-14)
            return false;

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        double f2 = f(x2);

        if (std::abs(f2) < tol) {
            *root = x2;
            return true;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }
}