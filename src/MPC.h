#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

const size_t N  = 15;
const double dt = 0.05;
const int latency_dt = 2;

// Cost scaling factors
const double scaling_acceleration = 11.0, scaling_delta_d = 400.0, scaling_acceleration_d  = 2.0, scaling_delta  = 1.0;

struct Prediction {
    vector<double>  X, Y, psi, v, cte, e_psi, delta, a;
};

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  Prediction Pred(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  double last_delta {0}, last_acc {0.1};

};

#endif /* MPC_H */
