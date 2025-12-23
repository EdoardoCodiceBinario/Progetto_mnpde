/*
Problema di trasporto–diffusione–reazione con campo di trasporto intenso
e diffusione ridotta. 
Non sappiamo la soluzione esatta.
In questo caso il metodo di stabilizzazione è applicato ed
è GLS (Galerkin Least Squares).
Notiamo che abbiamo accesso a più valori di a
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

// librerie per es 2 
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/base/symmetric_tensor.h>



using namespace dealii;

template <int dim>
struct DiffTrasParameters
{
  DiffTrasParameters()
  {
    prm.enter_subsection("Problem parameters");
    {
      prm.add_parameter("Finite element degree", fe_degree);
      prm.add_parameter("Initial refinement", initial_refinement);
      prm.add_parameter("Number of cycles", n_cycles);
      prm.add_parameter("Right hand side expression", rhs_expression);
      prm.add_parameter("Coefficiente diffusione", diffusione_coefficiente);
    }
    prm.leave_subsection();
    try
      {
        prm.parse_input("parametri_boh_" + std::to_string(dim) + "d.prm");
      }
    catch (std::exception &exc)
      {
        prm.print_parameters("parametri_boh_" + std::to_string(dim) + "d.prm");
        prm.parse_input("parametri_boh_" + std::to_string(dim) + "d.prm");
      }
      
    
    std::cout << "Valore di a (diffusione_coefficiente): " << diffusione_coefficiente << std::endl;
    
    std::cout << "Right hand side expression: " << rhs_expression << std::endl;
    std::cout << "CASO GLS" << std::endl;
    
    std::map<std::string, double> constants;
    constants["pi"] = numbers::PI;
    constants["a"] = diffusione_coefficiente;
    
    rhs_function.initialize(FunctionParser<dim>::default_variable_names(),
                            {rhs_expression},
                            constants);                                            
  }
  
  unsigned int fe_degree                 = 2;
  unsigned int initial_refinement        = 2;
  unsigned int n_cycles                  = 6;
  std::string  rhs_expression            = "-3*exp(-(x+y)/g)/g";
  double  diffusione_coefficiente = 0.001; 
  
  FunctionParser<dim> rhs_function;

  ParameterHandler prm;
};


template <int dim>
double c_coefficient(const Point<dim> &p)
{
    double cx = std::exp(-50 * ((p[0]-0.5)*(p[0]-0.5) + (p[1]-0.5)*(p[1]-0.5)));
    double cy = 0.5 * std::sin(5 * numbers::PI * p[0]) * std::cos(5 * numbers::PI * p[1]);
    return 50.0*cx + 30.0*cy;
}


double funzione_peclet(double peclet)
{
    if (std::abs(peclet) < 1e-10) {
        return 0.0;
    }
    return 1.0 / std::tanh(peclet) - 1.0 / peclet;
}
template <int dim>
Tensor<1, dim> b_coefficient(const Point<dim> &p)
{
    Tensor<1, dim> b;
    const double x0 = 0.5;
    const double y0 = 0.5;
    const double strength = 30.0;
    const double k = 5.0 * numbers::PI;

    // componente vortice
    double vort_x =  strength * (p[1]-y0);
    double vort_y = -strength * (p[0]-x0);

    // componente ondulata
    double onda_x = 50 * std::sin(k * p[1]);
    double onda_y = 50 * std::cos(k * p[0]);

    b[0] = vort_x + onda_x;
    b[1] = vort_y + onda_y;
    return b;
}




template <int dim>
class DiffTras
{
public:
  DiffTras(const DiffTrasParameters<dim> &parameters);
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results(const unsigned int cycle) const;

  const DiffTrasParameters<dim> &par;

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
DiffTras<dim>::DiffTras(const DiffTrasParameters<dim> &par)
  : par(par)
  , fe(par.fe_degree)
  , dof_handler(triangulation)
{}



template <int dim>
void DiffTras<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(par.initial_refinement);

  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}


template <int dim>
void DiffTras<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  constraints.clear();
  
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);


  VectorTools::interpolate_boundary_values(dof_handler,
    0,         
    Functions::ZeroFunction<dim>(),
    constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
  sparsity_pattern.copy_from(dsp);


  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}



template <int dim>
void DiffTras<dim>::assemble_system()
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

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  float tau=1;  
  float h=1;
  float peclet;
  float bnorm=1;  
  Tensor<1,dim> b;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      h=cell->diameter();
 

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        for (const unsigned int i : fe_values.dof_indices())
          {

            float lap_i=trace(fe_values.shape_hessian(i, q_index));

            const auto &x_q = fe_values.quadrature_point(q_index);

            b = b_coefficient(x_q);
            bnorm = std::sqrt(b*b);
            peclet=(h*bnorm)/(2*par.diffusione_coefficiente); 
           tau = (h / (2*bnorm)) * funzione_peclet(peclet);
            for (const unsigned int j : fe_values.dof_indices())
              {
                // termine di diffusione 
                cell_matrix(i, j) +=
                  par.diffusione_coefficiente *
                  (fe_values.shape_grad(j, q_index) * 
                   fe_values.shape_grad(i, q_index) * 
                   fe_values.JxW(q_index));
                         
                // termine di trasporto
                cell_matrix(i, j) +=
                  (b *
                   fe_values.shape_grad(j, q_index) *
                   fe_values.shape_value(i, q_index) * 
                   fe_values.JxW(q_index));

                // termine reazione  
                cell_matrix(i, j) +=
                  c_coefficient(x_q) *
                  (fe_values.shape_value(j, q_index) * 
                   fe_values.shape_value(i, q_index) * 
                   fe_values.JxW(q_index)); 
                   
                    float lap_j=trace(fe_values.shape_hessian(j, q_index));

               //INIZIO PARTE STABILIZZAZIONE 
               
               cell_matrix(i,j) += 
               par.diffusione_coefficiente* par.diffusione_coefficiente* // a^2
               tau*                                                     //tau
               lap_i*                       //lap_phi_i(x_q)
               lap_j*                       //lap_phi_j(x_q)   
               fe_values.JxW(q_index);      //dx

               cell_matrix(i,j) += 
               -par.diffusione_coefficiente*                             //-a
               tau*                                                     //tau
               lap_i*                                                   //lap_phi_i(x_q)  
               (b*fe_values.shape_grad(j, q_index))*   // b*grad_phi_j(x_q)
               fe_values.JxW(q_index); 


               cell_matrix(i,j) += 
               -par.diffusione_coefficiente*                             //-a
               tau*                                                     //tau
               lap_j*                                                   //lap_phi_j(x_q)  
               (b*fe_values.shape_grad(i, q_index))*   // b*grad_phi_i(x_q) 
               fe_values.JxW(q_index); 


                cell_matrix(i,j)+=
                tau*
                (b * fe_values.shape_grad(i,q_index))*
                (b * fe_values.shape_grad(j,q_index))*
                fe_values.JxW(q_index); 
                    
              // FINE PARTE STABILIZZAZIONE
              
              // FINE PARTE STABILIZZAZIONE
              }   

              
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * 
                            par.rhs_function.value(x_q) *
                            fe_values.JxW(q_index));

             // INIZIO PARTE STABILIZZAZIONE 
                   
             cell_rhs(i) += -par.diffusione_coefficiente*
             tau*
             par.rhs_function.value(x_q)*
             lap_i*
             fe_values.JxW(q_index);       
              
cell_rhs(i) +=tau*
             par.rhs_function.value(x_q)*
             (b* fe_values.shape_grad(i, q_index))*
             fe_values.JxW(q_index);      
                         
            // FINE PARTE STABILIZZAZIONE                        
          }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}



template <int dim>
void DiffTras<dim>::solve()
{
  PreconditionIdentity preconditioner;
  SolverControl solver_control(10000000, 1e-12);
  SolverGMRES<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, preconditioner);
  constraints.distribute(solution);  
  std::cout << solver_control.last_step()
            << " GMRES iterations needed to obtain convergence." << std::endl;
}
                                      

template <int dim>
void DiffTras<dim>::output_results(const unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  auto fname =
    "solution-boh-GLS-" + std::to_string(dim) + "d_" + std::to_string(cycle) + ".vtu";

  std::ofstream output(fname);
  data_out.write_vtu(output);

  static std::vector<std::pair<double, std::string>> times_and_names;
  times_and_names.push_back({static_cast<double>(cycle), fname});

  std::ofstream pvd_output("solution-boh-gls" + std::to_string(dim) + "d.pvd");
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}



template <int dim>
void DiffTras<dim>::run()
{
  std::cout << "Solving problem in " << dim << " space dimensions."
            << std::endl;
  for (unsigned int cycle = 0; cycle < par.n_cycles; ++cycle)
    {
      std::cout << "ciclo: " << cycle << std::endl;

      if (cycle == 0)
        make_grid();
      else
        {
          Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
          triangulation.refine_global(1);
        }
      setup_system();
      assemble_system();
      solve();
      output_results(cycle);
    }
}


int main()
{
  {
    DiffTrasParameters<2> par;
    DiffTras<2>           laplace_problem_2d(par);
    laplace_problem_2d.run();
  }

  return 0;
}