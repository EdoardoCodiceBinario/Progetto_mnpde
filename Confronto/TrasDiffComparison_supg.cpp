/*
Questo file risolve il problema 
-a*laplac(u)+b*grad(u))+c*u=f(x,y) su Omega
dove b è il termine di trasporto, c quello di reazione e a il termine di diffusione.
In particolare in questo file andiamo a introdurre due parametri al fine di confrontare il
fenomeno di dominazione per diffusione e per trasporto.
Noteremo che in questo caso, cioè dove la stabilizzazione è stata fatta mediante il metodo
Streamline Upwind Petrov Galerkin, il fenomeno di dominazione mediante trasporto è stabilizzato
in maniera migliore di quanto aspettato rispetto alla stima data dal Lemma di Cea e 
dal Lemma Bramble-Hilbert.
*/




#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_convergence_table.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/base/numbers.h>
#include <chrono>



using namespace dealii;

template <int dim>
struct TrasdiffParameters
{
  TrasdiffParameters()
  {
    prm.enter_subsection("Problem parameters");
    {
      prm.add_parameter("Finite element degree", fe_degree);
      prm.add_parameter("Initial refinement", initial_refinement);
      prm.add_parameter("Number of cycles", n_cycles);
      prm.add_parameter("Exact solution expression", exact_solution_expression);
      prm.add_parameter("Right hand side expression", rhs_expression);
      prm.add_parameter("Coefficiente diffusione", diffusione_coefficiente);
      prm.add_parameter("Coefficiente trasporto", trasporto_coefficiente);

    }
    prm.leave_subsection();

    prm.enter_subsection("Convergence table");
    convergence_table.add_parameters(prm);
    prm.leave_subsection();

    try
      {
        prm.parse_input("parametri_comp_" + std::to_string(dim) + "d.prm");
      }
    catch (std::exception &exc)
      {
        prm.print_parameters("parametri_comp_" + std::to_string(dim) + "d.prm");
        prm.parse_input("parametri_comp" + std::to_string(dim) + "d.prm");
      }
    std::map<std::string, double> constants;
    constants["pi"] = numbers::PI;
    constants["a"] = diffusione_coefficiente;
    constants["b"] = trasporto_coefficiente;
    exact_solution.initialize(FunctionParser<dim>::default_variable_names(),
                              {exact_solution_expression},
                              constants);
    rhs_function.initialize(FunctionParser<dim>::default_variable_names(),
                            {rhs_expression},
                            constants);                                              
  }
  unsigned int fe_degree                 = 1;
  unsigned int initial_refinement        = 3;
  unsigned int n_cycles                  = 1;
  std::string  exact_solution_expression = "exp(-(x+y)/a)";
  std::string  rhs_expression            = "-(2+2*b)*exp(-(x+y)/a)/a";

  float diffusione_coefficiente = 1.0; 
  float trasporto_coefficiente=0.1; 
  FunctionParser<dim> exact_solution;
  FunctionParser<dim> rhs_function;

  mutable ParsedConvergenceTable convergence_table;

  ParameterHandler prm;
};

template <int dim>
Tensor<1, dim> b_coefficient(const Point<dim> &p)
{

  Tensor<1, dim> b;
  (void) p;  
  b[0]=1;
  b[1]=1;
  return b;
}

template <int dim>
double c_coefficient(const Point<dim> &p)
{
  (void) p;
  return 1.0;
}



float funzione_peclet(float peclet)
{
    if (peclet < 1e-10) {
        return 0.0;
    }
    return (1.0 / std::tanh(peclet)) - 1.0 / peclet;
}

template <int dim>
class Trasdiff
{
public:
  Trasdiff(const TrasdiffParameters<dim> &parameters);
  void
  run();

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results(const unsigned int cycle) const;

  const TrasdiffParameters<dim> &par;

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

  AffineConstraints<double> constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
};



template <int dim>
Trasdiff<dim>::Trasdiff(const TrasdiffParameters<dim> &par)
  : par(par)
  , fe(par.fe_degree)
  , dof_handler(triangulation)
{}



template <int dim>
void
Trasdiff<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(par.initial_refinement);

  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}


template <int dim>
void
Trasdiff<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  constraints.clear();
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           par.exact_solution,
                                           constraints);


  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void
Trasdiff<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                          update_hessians |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  float tau=1;  
  float h=1;
  float peclet;
  float bnorm=1;

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;
      h=cell->diameter();
      

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        for (const unsigned int i : fe_values.dof_indices())
          {
            const auto &x_q = fe_values.quadrature_point(q_index);
            
            bnorm=par.trasporto_coefficiente*b_coefficient(x_q).norm();
            peclet=(h*bnorm)/(2*par.diffusione_coefficiente);      
            tau=(h/(2*bnorm))*funzione_peclet(peclet);
                
            for (const unsigned int j : fe_values.dof_indices())
              {

                // termine di diffusione 
               cell_matrix(i, j) +=
               par.diffusione_coefficiente*
               (fe_values.shape_grad(j, q_index) * // grad phi_i(x_q)
                fe_values.shape_grad(i, q_index) * // grad_phi_j(x_q)
                fe_values.JxW(q_index));           // dx
                         
              
               // termine di trasporto
              cell_matrix(i, j) +=
                (par.trasporto_coefficiente*b_coefficient(x_q)*
                fe_values.shape_grad(j, q_index) *   // grad phi_i(x_q)
                 fe_values.shape_value(i, q_index) * // phi_j(x_q)
                 fe_values.JxW(q_index));            // dx

              // termine reazione  
                cell_matrix(i, j) +=
                c_coefficient(x_q) *                // c(x_q)
               (fe_values.shape_value(j, q_index) * // grad phi_i(x_q)
                fe_values.shape_value(i, q_index) * // grad_phi_j(x_q)
                fe_values.JxW(q_index));       

                

               //INIZIO PARTE STABILIZZAZIONE 
               

                cell_matrix(i,j)+=
                tau*                                                                                // tau
                (par.trasporto_coefficiente*b_coefficient(x_q) * fe_values.shape_grad(i,q_index))*  // b* gradi
                (par.trasporto_coefficiente*b_coefficient(x_q) * fe_values.shape_grad(j,q_index))*  // b* gradj
                fe_values.JxW(q_index);                                                             //dx

                
                cell_matrix(i,j)+=
                tau*                                                                                //tau 
                (par.trasporto_coefficiente*b_coefficient(x_q) * fe_values.shape_grad(i,q_index))*  // b*gradi   
                c_coefficient(x_q)*fe_values.shape_value(j,q_index)*                                // c*phi_j
                fe_values.JxW(q_index);                                                             // dx
                

              // FINE PARTE STABILIZZAZIONE
              
              }   

            cell_rhs(i) += (fe_values.shape_value(i, q_index)  * // phi_i(x_q)
                            par.rhs_function.value(x_q) *       // f(x_q)
                            fe_values.JxW(q_index));            // dx

            // INIZIO PARTE STABILIZZAZIONE 
                         
              
            cell_rhs(i) +=tau*                                                                                 // tau
                          par.rhs_function.value(x_q)*                                                         // f(x_q)
                          par.trasporto_coefficiente*
                          (b_coefficient(x_q) * fe_values.shape_grad(i, q_index))*                              // b*grad_i
                          fe_values.JxW(q_index);                                                              // dx     
                          
            // FINE PARTE STABILIZZAZIONE        
            
          }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}





template <int dim>
void Trasdiff<dim>::solve()
{
  PreconditionIdentity preconditioner;
  SolverControl solver_control(100000000, 1e-12);
  SolverGMRES<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs,preconditioner);
  constraints.distribute(solution);  
  std::cout << solver_control.last_step()
            << " GMRES iterations needed to obtain convergence." << std::endl;
}
                                      



template <int dim>
void
Trasdiff<dim>::output_results(const unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  auto fname =
    "solution-comparison-supg-" + std::to_string(dim) + "d_" + std::to_string(cycle) + ".vtu";

  std::ofstream output(fname);
  data_out.write_vtu(output);

  static std::vector<std::pair<double, std::string>> times_and_names;
  times_and_names.push_back({cycle, fname});

  std::ofstream pvd_output("solution-comparison-supg-" + std::to_string(dim) + "d.pvd");

  DataOutBase::write_pvd_record(pvd_output,

 times_and_names);
}



template <int dim>
void
Trasdiff<dim>::run()
{
  Point<dim> p(-1,-1);
  Point<dim> center(0,0);
  std::cout << "Solving problem in " << dim << " space dimensions."
            << std::endl;
  for (unsigned int cycle = 0; cycle < par.n_cycles; ++cycle)
    {
      std::cout << "ciclo: " << cycle << std::endl;


      if (cycle == 0)

        make_grid();
      else 
        {
    triangulation.refine_global(1);
        }
      setup_system();
      assemble_system();
      solve();
      output_results(cycle);
      par.convergence_table.error_from_exact(dof_handler,
        solution,
        par.exact_solution);                         
                                             
                                                                                      
                                               
    }
    par.convergence_table.output_table(std::cout);                         
                                             
                                                                                      
                            
}


int
main()
{
  {
    TrasdiffParameters<2> par;
    Trasdiff<2>           laplace_problem_2d(par);
    laplace_problem_2d.run();
  }

  return 0;
}