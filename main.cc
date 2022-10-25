//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Basic framework for modeling action potential propagation in neurons using the single-cable model.
//Created May 2021
//authors: Rahul Gulati (2021)
//

//deal.II headers
#include "include/headers.h"
//input parameter headers
#include "parameters.h"
//physics headers
#include "include/mechanics.h"

//Namespace
namespace elasticity1
{
  using namespace dealii;

    //Initial conditions
  template <int dim>
  class InitalConditions: public Function<dim>{
  public:
    InitalConditions (): Function<dim>(DIMS){}
    void vector_value (const Point<dim> &p , Vector<double>   &values) const{
      Assert (values.size() == DIMS, ExcDimensionMismatch (values.size(),DIMS));
      values(0) = V_rest; 
      values(1) = 0.0;
    }
    };

  
  template <int dim>
  class elasticity{
  public:
    elasticity ();
    ~elasticity ();
    void run ();

  private:
    void applyBoundaryConditions(const unsigned int increment);
    void setup_system ();
    void assemble_system (const unsigned int increment, const unsigned int iteration);
    void solveIteration ();
    void solve (const unsigned int increment);
    void L2_projection();
    void L2_solveIteration();
    void refine_grid ();
    void output_results (const unsigned int increment);
    Triangulation<dim>                        triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    AffineConstraints<double>                      constraints, constraints2, constraints3, constraints_L2;

    SparsityPattern                         sparsity_pattern, sparsity_pattern_L2;
    SparseMatrix<double>                    system_matrix, mass_matrix_L2;
    Vector<double>                          dU, Un, U;
    Vector<double>                          system_rhs, system_rhs_current, system_rhs_n, system_rhs_m, system_rhs_h, U_L2, U_current, U_n_iter, U_m_iter, U_h_iter;
  
    //solution variables
    unsigned int currentIncrement, currentIteration;
    double totalTime, currentTime, dt;
    std::vector<std::string> nodal_solution_names; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
     std::vector<std::string> nodal_solution_names_L2; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation_L2;
    //history variables
    std::map<typename Triangulation<dim>::active_cell_iterator, std::vector<history<dim>*> >  history_variables; 
  };

  template <int dim>
  elasticity<dim>::elasticity ():
    fe(FE_Q<dim>(2),DIMS),
    dof_handler (triangulation){
    currentIncrement=0; currentTime=0;
    
    //nodal Solution names
    for (unsigned int i=0; i<dim; ++i){
      nodal_solution_names.push_back("u"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
       }
    //current
    char buffer[100];
    for(char i=0;i<1;i++){
      sprintf(buffer, "Current");
      nodal_solution_names.push_back(buffer);nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    }

      nodal_solution_names_L2.push_back("n_value"); nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);
  }
  
  template <int dim>
  elasticity<dim>::~elasticity (){
    dof_handler.clear ();
  }

  //Apply boundary conditions
  template <int dim>
  void elasticity<dim>::applyBoundaryConditions(const unsigned int increment){
    constraints.clear (); constraints2.clear (); constraints3.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints2);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints3);

    //L2 boundary conditions
    constraints_L2.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, constraints_L2);
    constraints_L2.close ();
    
    //Setup boundary conditions
      std::vector<bool> uBCX0 (DIMS, false); uBCX0[0]=true; 
      VectorTools::interpolate_boundary_values (dof_handler, 0, ConstantFunction<dim>(VBoundary, DIMS), constraints, uBCX0);
      VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(DIMS), constraints3, uBCX0);
      //VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints);
      //VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints2);
    constraints.close ();
    constraints2.close ();
    constraints3.close ();
  }
  
  //Setup
  template <int dim>
  void elasticity<dim>::setup_system (){
    //TimerOutput::Scope t(computing_timer, "setup");
    dof_handler.distribute_dofs (fe);
    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions(0);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    
    dU.reinit(dof_handler.n_dofs());
    U.reinit(dof_handler.n_dofs());
    Un.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    //L2 setup
    system_rhs_current.reinit (dof_handler.n_dofs());
    system_rhs_n.reinit (dof_handler.n_dofs());
    system_rhs_m.reinit (dof_handler.n_dofs());
    system_rhs_h.reinit (dof_handler.n_dofs());
    U_L2.reinit (dof_handler.n_dofs());
    U_current.reinit (dof_handler.n_dofs());
    U_n_iter.reinit (dof_handler.n_dofs());
    U_m_iter.reinit (dof_handler.n_dofs());
    U_h_iter.reinit (dof_handler.n_dofs());
    DynamicSparsityPattern dsp_L2 (dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp_L2); //check
    sparsity_pattern_L2.copy_from(dsp_L2);
    mass_matrix_L2.reinit(sparsity_pattern_L2);
    
    //setup history variables
    const QGauss<dim> quadrature_formula(3);
    FEValues<dim> fe_values (fe, quadrature_formula, update_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc= dof_handler.end();
    for (; cell!=endc; ++cell){
      for (unsigned int q=0; q<fe_values.n_quadrature_points; q++){
	history_variables[cell].push_back(new history<dim>);
      }
    }
      
}

  //Assembly
  template <int dim>
  void elasticity<dim>::assemble_system (const unsigned int increment, const unsigned int iteration){
    //TimerOutput::Scope t(computing_timer, "assembly");
    system_rhs=0.0; system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(3);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
     for (; cell!=endc; ++cell){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  if (std::abs(U(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  else{ULocal[i]=U(local_dof_indices[i]);}
	  ULocal[i].diff (i, dofs_per_cell);
	  ULocalConv[i]= Un(local_dof_indices[i]);
	}
	//get defomration map
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points);

	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell);
	Table<2, double> KMatrix(dofs_per_cell, dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i) {R[i]=0.0;}
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell; ++j) {
	  KMatrix[i][j]=0.0; }
	}
	//populate residual vector 
	residualForMechanics(fe_values, 0, ULocal, ULocalConv, R, KMatrix, defMap, cell, increment, history_variables[cell], iteration);
	//evaluate Residual(R) and Jacobian(R')
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell; ++j){
	    // R' by AD
	    //local_matrix(i,j)= R[i].fastAccessDx(j);
	    local_matrix(i,j)= KMatrix[i][j];
	  }
	  //R
	  local_rhs(i) = -R[i].val();
	}
	//
	if ((currentIteration==0) && ((currentIncrement==1) || (currentIncrement==1000)) ){
	  constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	}
	  else if ((currentIncrement==1)  || (currentIncrement==1000) ){
	   constraints3.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	}

	else{
	  constraints2.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	}
     }
    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }

    template <int dim>
  void elasticity<dim>::L2_projection (){
      //TimerOutput::Scope t(computing_timer, "L2_projection");
    system_rhs_current=0.0; system_rhs_n=0.0; system_rhs_m=0.0; system_rhs_h=0.0;
    mass_matrix_L2=0.0;
    const QGauss<dim>  quadrature_formula(3);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   |
                             update_quadrature_points |
                             update_JxW_values);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs_current (dofs_per_cell), local_rhs_n (dofs_per_cell), local_rhs_m (dofs_per_cell), local_rhs_h (dofs_per_cell); 
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs_current = 0;  local_rhs_n = 0;  local_rhs_m = 0;  local_rhs_h = 0; 
	cell->get_dof_indices (local_dof_indices);
	
	for (unsigned int q=0; q<n_q_points; q++){
	  for(unsigned int i=0;i<dofs_per_cell;i++){
	    const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first;
	    if (ci==0){
	      local_rhs_current(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->current*fe_values.JxW(q);
	      local_rhs_n(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->n_iter*fe_values.JxW(q);
	      local_rhs_m(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->m_iter*fe_values.JxW(q);
	      local_rhs_h(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->h_iter*fe_values.JxW(q);
	    }
	    for(unsigned int j=0;j<dofs_per_cell;j++){
	      const unsigned int cj = fe_values.get_fe().system_to_component_index(j).first;
	      if (ci==cj){
		local_matrix(i,j) += fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q);
	      }
	    }
	  }
	}
	constraints_L2.distribute_local_to_global (local_matrix, local_rhs_current, local_dof_indices, mass_matrix_L2, system_rhs_current);
	constraints_L2.distribute_local_to_global (local_rhs_n, local_dof_indices, system_rhs_n);
	constraints_L2.distribute_local_to_global (local_rhs_m, local_dof_indices, system_rhs_m);
	constraints_L2.distribute_local_to_global (local_rhs_h, local_dof_indices, system_rhs_h);
      }
    mass_matrix_L2.compress (VectorOperation::add);
    system_rhs_current.compress (VectorOperation::add);
    system_rhs_n.compress (VectorOperation::add);
    system_rhs_m.compress (VectorOperation::add);
    system_rhs_h.compress (VectorOperation::add);
    }


  

  //Solve
  template <int dim>
  void elasticity<dim>::solveIteration(){
    /*
    //TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    /*    
    //Iterative solvers from Petsc and Trilinos
    SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
    LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
    LA::SolverGMRES solver(solver_control);
#endif
    LA::MPI::PreconditionAMG preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA
    //data.symmetric_operator = true;
#else
    // Trilinos defaults are good 
#endif
    preconditioner.initialize(system_matrix, data);
    solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
    
    //Direct solver MUMPS
    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);
    //
    */
    // Works great:
    SolverControl         solver_control(1000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control);
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
    cg.solve(system_matrix, dU, system_rhs, preconditioner);
    
    /*    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(dU, system_rhs);
    */
    if ( (currentIteration==0) && ((currentIncrement==1) || (currentIncrement==1000)) ){
      constraints.distribute (dU);
    }
    else if ((currentIncrement==1) || (currentIncrement==1000)){
      constraints3.distribute (dU);
    }
    else{
      constraints2.distribute (dU);
      }
    //locally_relevant_solution = completely_distributed_solution;
    //dU = completely_distributed_solution; 
  }

   //Solve
  template <int dim>
  void elasticity<dim>::L2_solveIteration(){
    // Works great:
    SolverControl         solver_control_L2(1000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control_L2);
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(mass_matrix_L2, 1.2);
    //current
    cg.solve(mass_matrix_L2, U_current, system_rhs_current, preconditioner);
    constraints_L2.distribute (U_current);
    //n
    cg.solve(mass_matrix_L2, U_n_iter, system_rhs_n, preconditioner);
    constraints_L2.distribute (U_n_iter);
    //m
    cg.solve(mass_matrix_L2, U_m_iter, system_rhs_m, preconditioner);
    constraints_L2.distribute (U_m_iter);
    //h
    cg.solve(mass_matrix_L2, U_h_iter, system_rhs_h, preconditioner);
    constraints_L2.distribute (U_h_iter);
  }


  
  //Solve
  template <int dim>
  void elasticity<dim>::solve(const unsigned int increment){
    double res=1, tol=1.0e-10, abs_tol=1.0e-12, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];
    while (true){
      if (currentIteration>=3){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); std::cout<<buffer; break; exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); std::cout<<buffer; break; exit (1);}
      assemble_system(increment, currentIteration);
      current_norm=system_rhs.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); std::cout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); std::cout<<buffer; break;}
      solveIteration();
      U+=dU; //UGhost=U; 
      ++currentIteration;
    }
    Un=U; //UnGhost=Un;
  }

  //Output
  template <int dim>
  void elasticity<dim>::output_results (const unsigned int cycle) {
    //TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (Un, "Volt_mag"); //nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
    //data_out.add_data_vector (Un, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
    data_out.add_data_vector (U_current, "Current_mag"); //nodal_solution_names_L2, DataOut<dim>::type_dof_data, nodal_data_component_interpretation_L2);
    data_out.add_data_vector (U_n_iter, "n_mag");
    data_out.add_data_vector (U_m_iter, "m_mag");
    data_out.add_data_vector (U_h_iter, "h_mag");
    
    data_out.build_patches ();
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
  }

  //Solve problem
  template <int dim>
  void elasticity<dim>::run (){
    //setup problem geometry and mesh
    GridGenerator::hyper_cube (triangulation, 0, problemWidth, true);
    triangulation.refine_global (refinementFactor);
    
    setup_system ();
    std::cout << "   Number of active cells:       "
	  << triangulation.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;
    //setup initial conditions
    VectorTools::interpolate(dof_handler, InitalConditions<dim>(), U); Un=U;
    output_results (0);
    //Time stepping
    currentIncrement=0;
    for (currentTime=0; currentTime<Total_time; currentTime+=DT){
      currentIncrement++;
      applyBoundaryConditions(currentIncrement);
      solve(currentIncrement);
      L2_projection();
      L2_solveIteration();
      output_results(currentIncrement);
      std::cout << std::endl;
    }
    //computing_timer.print_summary ();
  }
}


int main(int argc, char *argv[]){
  try
    {
      using namespace dealii;
      using namespace elasticity1;
      //Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      elasticity<dimension> problem;
      problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
