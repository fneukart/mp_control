#include "MPC.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Set the timestep length and duration
// size_t N = 0;
// double dt = 0;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

const double cte_reference = 0;
const double e_psi_reference = 0;
const double v_reference  = 75;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


size_t  x_init     = 0;
size_t  y_init     = x_init + N;
size_t  psi_init   = y_init + N;
size_t  v_start     = psi_init + N;
size_t  cte_start   = v_start + N;
size_t  e_psi_init  = cte_start + N;
size_t  delta_init = e_psi_init + N;
size_t  a_init     = delta_init + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  vector<double> actuation_history;

  FG_eval(Eigen::VectorXd coeffs, vector<double> history ) { 
    this->coeffs = coeffs; 
    this->actuation_history = history;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    fg[0] = 0;
    for (int t = 0; t < N; t++) {
      //cte
      fg[0] +=  1.0 * CppAD::pow(vars[cte_start + t] - cte_reference, 2);
      //error in orientation
      fg[0] +=  1.0 * CppAD::pow(vars[e_psi_init + t] - e_psi_reference, 2);
      //error in velocity
      fg[0] +=  1.0 * CppAD::pow(vars[v_start + t] - v_reference, 2);
    }

    for (int t = 0; t < N - 1; t++) {
      //steering angle
      fg[0] +=  scaling_delta * CppAD::pow(vars[delta_init + t], 2);
      //acc
      fg[0] +=  scaling_acceleration   * CppAD::pow(vars[a_init + t], 2);
    }

    for (int t = 0; t < N - 2; t++) {
      fg[0] +=  scaling_delta_d * CppAD::pow(vars[delta_init + t + 1] - vars[delta_init + t], 2);
      fg[0] +=  scaling_acceleration_d   * CppAD::pow(vars[a_init + t + 1]  - vars[a_init + t], 2);
    }


    fg[1 + x_init]     = vars[x_init];
    fg[1 + y_init]     = vars[y_init];
    fg[1 + psi_init]   = vars[psi_init];
    fg[1 + v_start]     = vars[v_start];
    fg[1 + cte_start]   = vars[cte_start];
    fg[1 + e_psi_init]  = vars[e_psi_init];


    for (int t = 1; t < N ; t++) {
      AD<double>  x1    = vars[x_init + t];
      AD<double>  y1    = vars[y_init + t];
      AD<double>  psi1  = vars[psi_init + t];
      AD<double>  v1    = vars[v_start + t];
      AD<double>  cte1  = vars[cte_start + t];
      AD<double>  e_psi1 = vars[e_psi_init + t];

      AD<double>  x0    = vars[x_init + t - 1];
      AD<double>  y0    = vars[y_init + t - 1];
      AD<double>  psi0  = vars[psi_init + t - 1];
      AD<double>  v0    = vars[v_start + t - 1];
      AD<double>  cte0  = vars[cte_start + t - 1];
      AD<double>  e_psi0 = vars[e_psi_init + t - 1];

      AD<double>  delta0  = vars[delta_init + t - 1];
      AD<double>  a0      = vars[a_init + t - 1];

      AD<double>  f0      = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3]*x0 * x0 * x0;
      AD<double>  psides0 = CppAD::atan(coeffs[1] + 2*coeffs[2]*x0  + 3*coeffs[3]*x0*x0);

      fg[1 + x_init + t]   = x1    - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_init + t]   = y1    - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_init + t] = psi1  - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t]   = v1    - (v0 + a0 * dt);
      fg[1 + cte_start + t] = cte1  - ((f0 - y0) + (v0 * CppAD::sin(e_psi0) * dt));
      fg[1 + e_psi_init + t]= e_psi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }

  }
};
//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

Prediction MPC::Pred(Eigen::VectorXd x0, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  // N timesteps is N-1 actuations
  size_t n_vars = N * 6 + (N - 1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  double x    = x0[0];
  double y    = x0[1];
  double psi  = x0[2];
  double v    = x0[3];
  double cte  = x0[4];
  double e_psi = x0[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }

  vars[x_init]   = x;
  vars[y_init]   = y;
  vars[psi_init] = psi;
  vars[v_start]   = v;
  vars[cte_start] = cte;
  vars[e_psi_init]= e_psi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.

  for (int i = 0; i < delta_init; i++) {
    vars_lowerbound[i]  = -1.0e19;
    vars_upperbound[i]  = 1.0e19;
  }

  for (int i = delta_init; i < a_init; i++) {
    vars_lowerbound[i]  = -deg2rad(25);
    vars_upperbound[i]  = deg2rad(25);
  }

  for (int i = delta_init; i < delta_init + latency_dt; i++) {
    vars_lowerbound[i] = last_delta;
    vars_upperbound[i] = last_delta;
  }

  for (int i = a_init; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  for (int i = a_init; i < a_init + latency_dt; i++) {
    vars_lowerbound[i] = last_acc;
    vars_upperbound[i] = last_acc;
  }

  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_init]     = x;
  constraints_lowerbound[y_init]     = y;
  constraints_lowerbound[psi_init]   = psi;
  constraints_lowerbound[v_start]     = v;
  constraints_lowerbound[cte_start]   = cte;
  constraints_lowerbound[e_psi_init]  = e_psi;


  constraints_upperbound[x_init]     = x;
  constraints_upperbound[y_init]     = y;
  constraints_upperbound[psi_init]   = psi;
  constraints_upperbound[v_start]     = v;
  constraints_upperbound[cte_start]   = cte;
  constraints_upperbound[e_psi_init]  = e_psi;

  vector<double> actuation_history = { last_delta, last_acc};
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, actuation_history);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // 7 - place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // 8 - solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  Prediction pred;
  for (auto k = 0; k < N - 1; k++) {

	//relevant for the given simulator
    pred.X.push_back(solution.x[x_init + k]);
    pred.Y.push_back(solution.x[y_init + k]);
    pred.delta.push_back(solution.x[delta_init + k]);
    pred.a.push_back(solution.x[a_init + k]);

    //irrelevant for the given simulator
    pred.psi.push_back(solution.x[psi_init + k]);
    pred.v.push_back(solution.x[v_start + k]);
    pred.cte.push_back(solution.x[cte_start + k]);
    pred.e_psi.push_back(solution.x[e_psi_init + k]);
  }

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost     " << cost << std::endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  //
  return pred;
}
