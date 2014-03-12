#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// Step size
#define H 0.001
// Start time
#define T_start 0.0
// Finish time
#define T_finish 4.0
// Eye movement amount
#define Delta_g 2.0
// Time constant
#define T1 0.15
#define T2 0.012
#define TN 25
// Initial eye movement amount
#define G_init 0.0
// Velocity of eye movement
#define V_init 20.0
// Firing rate of neural integrator
#define N_init 0.5
// Firing rate of burst cells
#define B_init 5.0
// Firing rate of resettable neural integrator
#define S_init 1.0

// Runge-Kutta function
void d1(int i, double **solution,
        double (*saccade)(double, double, double, double, double)) {
  solution[1][i] = saccade(solution[0][0], solution[0][1], solution[0][2],
                           solution[0][3], solution[0][4]);
}

void d2(int i, double **solution,
        double (*saccade)(double, double, double, double, double)) {
  solution[2][i] = saccade(solution[0][0] + (H / 2.0) * solution[1][0],
                           solution[0][1] + (H / 2.0) * solution[1][1],
                           solution[0][2] + (H / 2.0) * solution[1][2],
                           solution[0][3] + (H / 2.0) * solution[1][3],
                           solution[0][4] + (H / 2.0) * solution[1][4]);
}

void d3(int i, double **solution,
        double (*saccade)(double, double, double, double, double)) {
  solution[3][i] = saccade(solution[0][0] + (H / 2.0) * solution[2][0],
                           solution[0][1] + (H / 2.0) * solution[2][1],
                           solution[0][2] + (H / 2.0) * solution[2][2],
                           solution[0][3] + (H / 2.0) * solution[2][3],
                           solution[0][4] + (H / 2.0) * solution[2][4]);
}

void d4(int i, double **solution,
        double (*saccade)(double, double, double, double, double)) {
  solution[4][i] = saccade(
      solution[0][0] + H * solution[3][0], solution[0][1] + H * solution[3][1],
      solution[0][2] + H * solution[3][2], solution[0][3] + H * solution[3][3],
      solution[0][4] + H * solution[3][4]);
}

// Weighting
void weighting(int i, double **solution) {
  solution[0][i] =
      solution[0][i] + H * (solution[1][i] + 2.0 * solution[2][i] +
                            2.0 * solution[3][i] + solution[4][i]) /
                           6.0;
}

// Approximation function for saccade response of brainstem cells
double cells(double x) {
  if (x < 0.0) {
    x = -800 * (1.0 - exp(-fabs(x) / 6.0));
  } else {
    x = 800 * (1.0 - exp(-fabs(x) / 6.0));
  }

  return (x);
}

// g, v, n, b, s differential equation
double fg(double g, double v, double n, double b, double s) { return (v); }

double fv(double g, double v, double n, double b, double s) {
  return (-v * (1.0 / T1 + 1.0 / T2) + (-g + n + (T1 + T2) * b) / (T1 * T2));
}

double fn(double g, double v, double n, double b, double s) {
  return (-n / TN + b);
}

double fb(double g, double v, double n, double b, double s) {
  return ((-b + cells(Delta_g - s)) * s);
}

double fs(double g, double v, double n, double b, double s) { return (b); }

int main() {
  int i, j;
  double (*saccade[5])(double, double, double, double, double) = {fg, fv, fn,
                                                                  fb, fs};
  void (*rungekutta[4])(int, double * *solution,
                        double (*saccade)(double, double, double, double,
                                          double)) = {d1, d2, d3, d4};

  double **solution = new double *[5];
  for (i = 0; i < 5; ++i) {
    solution[i] = new double[5];
  }

  double t = T_start;
  solution[0][0] = G_init;
  solution[0][1] = V_init;
  solution[0][2] = N_init;
  solution[0][3] = B_init;
  solution[0][4] = S_init;

  ofstream out("saccade.dat");

  for (t = 0.0; t <= T_finish; t += H) {
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 5; ++j) {
        // Runge-Kutta method
        rungekutta[i](j, solution, saccade[j]);
      }
    }
    for (j = 0; j < 5; ++j) {
      // Weighting
      weighting(j, solution);
    }
    out << t << " " << solution[0][0] << endl;
  }

  out.close();

  for (i = 0; i < 5; ++i) {
    delete[] solution[i];
  }
  delete[] solution;

  return 0;
}
