
#include <cmath>
#include <cstdlib>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>

#include <deal.II/base/numbers.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/grid/grid_in.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace dealii;

template <int dim> class LaplaceSolver {
public:
  LaplaceSolver(Triangulation<dim> &, int);
  void run(const std::map<types::boundary_id, const Function<dim> *> &, int,
           const Function<dim> &);
  void estimate_error(const Function<dim> &, int);

private:
  void setup_system();
  void
  assemble_system(const std::map<types::boundary_id, const Function<dim> *> &);
  void solve();
  void output_results(int, int);

  int problem;

  Triangulation<dim> triangulation;
  const FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
};

template <int dim>
LaplaceSolver<dim>::LaplaceSolver(Triangulation<dim> &triangulation_,
                                  int problem_)
    : problem(problem_), fe(1), dof_handler(triangulation_){};

template <int dim> void LaplaceSolver<dim>::setup_system() {

  dof_handler.distribute_dofs(fe);

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
};

template <int dim>
void LaplaceSolver<dim>::assemble_system(
    const std::map<types::boundary_id, const Function<dim> *>
        &boundary_function) {

  const QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    fe_values.reinit(cell);

    cell_matrix = 0.;

    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) +=
              (fe_values.shape_grad(i, q_index) *
               fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index));
      }
    }

    cell->get_dof_indices(local_dof_indices);
    for (const unsigned int i : fe_values.dof_indices()) {
      for (const unsigned int j : fe_values.dof_indices())
        system_matrix.add(local_dof_indices[i], local_dof_indices[j],
                          cell_matrix(i, j));
    }
  }

  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler, boundary_function,
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution,
                                     system_rhs);
}

template <int dim> void LaplaceSolver<dim>::solve() {

  SolverControl solver_control(5000, 1e-7 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}

template <int dim>
void LaplaceSolver<dim>::output_results(int problem, int cycle) {

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  {
    std::ofstream output("solution-" + std::to_string(problem) + "-cycle-" +
                         std::to_string(cycle) + ".vtu");
    data_out.write_vtu(output);
  }
}

template <int dim>
void LaplaceSolver<dim>::run(
    const std::map<types::boundary_id, const Function<dim> *>
        &boundary_function,
    int cycle, const Function<dim> &exact_solution) {

  setup_system();
  assemble_system(boundary_function);
  solve();
  output_results(problem, cycle);
  estimate_error(exact_solution, cycle);
}

template <int dim>
void LaplaceSolver<dim>::estimate_error(const Function<dim> &exact_solution,
                                        int cycle) {

  Vector<float> difference_per_cell(triangulation.n_active_cells());
  VectorTools::integrate_difference(
      dof_handler, solution, exact_solution, difference_per_cell,
      QGauss<dim>(fe.degree + 2), VectorTools::H1_norm);

  double H1_error = difference_per_cell.l2_norm();
  std::cout << "H1 error = " << H1_error << std::endl;

  VectorTools::integrate_difference(
      dof_handler,
      solution,       // numerical solution u_h
      exact_solution, // Function<dim>
      difference_per_cell, QGauss<dim>(fe.degree + 2), VectorTools::L2_norm);

  double L2_error = difference_per_cell.l2_norm();
  std::cout << "L2 error = " << L2_error << std::endl;

  std::ofstream file("problem-" + std::to_string(problem) + ".txt",
                     std::ios::app); // open in append mode

  if (file.is_open()) {
    file << cycle << "," << L2_error << "," << H1_error << "\n";
    file.close();
  } else {
    std::cerr << "Unable to open file\n";
  }
}

/*
 *  Problem 1
 *  Infinitely long cuboidal conductor (z-invariant) with
 *  rectangular cross-section 0<x<w, 0<y<h.
 *
 *  Solve Laplace equation:
 *      ∂²u/∂x² + ∂²u/∂y² = 0
 *
 *  Boundary conditions:
 *      u(0,y) = 0
 *      u(w,y) = 0
 *      u(x,0) = 0
 *      u(x,h) = sin(x)
 *
 *  Find u(x,y) using separation of variables.
 *
 *
 *  Theoretical solution: u(x,y) = (sinh(y)/sinh(h)) sin(x) if w = pi
 */

int height = 2 * numbers::PI;
int width = numbers::PI;

class Boundary1 : public Function<2> {
public:
  virtual double value(const Point<2> &p, const unsigned int) const override {

    if (abs(p[1] - height) < 1e-6) {
      return std::sin(p[0]);
    }
    return 0;
  }
};

class ExactSolution1 : public Function<2> {
public:
  virtual double value(const Point<2> &p, const unsigned int) const override {

    return sin(p[0]) * std::sinh(p[1]) / std::sinh(height);
  }
  virtual Tensor<1, 2>
  gradient(const Point<2> &p, const unsigned int component = 0) const override {
    Tensor<1, 2> grad;

    grad[0] = std::sinh(p[1]) / std::sinh(height) * std::cos(p[0]); // du/dx
    grad[1] = std::cosh(p[1]) / std::sinh(height) * std::sin(p[0]); // du/dy

    return grad;
  }
};

void make_grid1(Triangulation<2> &triangulation) {

  GridGenerator::hyper_rectangle(triangulation, Point<2>(0, 0),
                                 Point<2>(numbers::PI, height));
  triangulation.refine_global(1);
}

void make_grid2(Triangulation<3> &triangulation) {
  GridGenerator::hyper_cube<3>(triangulation);
  triangulation.refine_global();
}

class ExactSolution2 : public Function<3> {
public:
  virtual double value(const Point<3> &p, const unsigned int) const override {

    return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]) *
           std::sinh(numbers::PI * p[2] * numbers::SQRT2);
  }
  virtual Tensor<1, 3>
  gradient(const Point<3> &p, const unsigned int component = 0) const override {

    Tensor<1, 3> grad;

    grad[0] = numbers::PI * std::cos(numbers::PI * p[0]) *
              std::sin(numbers::PI * p[1]) *
              std::sinh(std::sqrt(2.0) * numbers::PI * p[2]);

    grad[1] = numbers::PI * std::sin(numbers::PI * p[0]) *
              std::cos(numbers::PI * p[1]) *
              std::sinh(std::sqrt(2.0) * numbers::PI * p[2]);

    grad[2] = std::sqrt(2.0) * numbers::PI * std::sin(numbers::PI * p[0]) *
              std::sin(numbers::PI * p[1]) *
              std::cosh(std::sqrt(2.0) * numbers::PI * p[2]);

    return grad;
  }
};

int main() {

  if (false) {

    Triangulation<2> triangulation;
    make_grid1(triangulation);
    LaplaceSolver<2> problem1(triangulation, 1);

    std::map<types::boundary_id, const Function<2> *> boundary_map;
    Boundary1 boundary;
    boundary_map[0] = &boundary;

    ExactSolution1 exact_solution;

    for (int i = 1; i < 15; i++) {
      problem1.run(boundary_map, i, exact_solution);
      triangulation.refine_global();
    }
  }

  {
    Triangulation<3> triangulation;
    make_grid2(triangulation);
    LaplaceSolver<3> problem2(triangulation, 2);
    std::map<types::boundary_id, const Function<3> *> boundary_map;

    ExactSolution2 exact_solution;

    boundary_map[0] = &exact_solution;

    for (int i = 1; i < 9; i++) {
      problem2.run(boundary_map, i, exact_solution);
      triangulation.refine_global();
    }
  }
}
