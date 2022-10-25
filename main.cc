//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Basic framework for modeling spatio-temporal action potential propagation in neurons using electro-diffusive unmyelinated PNP model.
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
      
      if ( (  (p[1] <= (1.6+memThick)) && (p[1] >= (1.6-memThick))  ) || ( (p[1] <= (0.5+memThick)) && (p[1] >= (0.5-memThick))          )   ){
        values(0) = 0.00;
        values(1) = 0.00;
        values(2) = 0.00;
        values(3) = 0.0;
      }

      else if ( (p[1]<0.5) || (p[1]>1.6) ){
        values(0) = initialExtraConcNa;
        values(1) = initialExtraConcK;
        values(2) = initialExtraConcA;
        values(3) = 0.00;
      }
      else {
        values(0) = initialIntraConcNa;
        values(1) = initialIntraConcK;
        values(2) = initialIntraConcA;
        values(3) = V_rest;
      }
      
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
    //void setup_system_history_variables ();
    void assemble_system (const unsigned int increment, const unsigned int iteration);
    void solveIteration ();
    void solve (const unsigned int increment);
    void L2_projection();
    void L2_solveIteration();
    void refine_grid ();
    void output_results (const unsigned int increment);
    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    //Triangulation<dim>                        triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;
    AffineConstraints<double>               constraints, constraints2, constraints3, constraints_L2;
    LA::MPI::SparseMatrix                     system_matrix, mass_matrix_L2;
    LA::MPI::Vector                           locally_relevant_solution, U, Un, UGhost, UnGhost, dU, UGhost_L2;
    LA::MPI::Vector                           system_rhs, system_rhs_L2, locally_relevant_solution_L2, system_rhs_current, system_rhs_n, system_rhs_m, system_rhs_h, U_L2, U_current, U_n_iter, U_m_iter, U_h_iter;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
  
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
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
    fe(FE_Q<dim>(1),DIMS),
    dof_handler (triangulation),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times){
    currentIncrement=0; currentTime=0;
    
    //nodal Solution names
      nodal_solution_names.push_back("c_Na"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
      nodal_solution_names.push_back("c_K"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
      nodal_solution_names.push_back("c_A"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
      nodal_solution_names.push_back("volt"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

      nodal_solution_names_L2.push_back("n_value");nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar); 
      nodal_solution_names_L2.push_back("m_value");nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);
      nodal_solution_names_L2.push_back("h_value");nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);
      nodal_solution_names_L2.push_back("h2_value");nodal_data_component_interpretation_L2.push_back(DataComponentInterpretation::component_is_scalar);
        
  }
  
  template <int dim>
  elasticity<dim>::~elasticity (){
    dof_handler.clear ();
  }

  //Apply boundary conditions
  template <int dim>
  void elasticity<dim>::applyBoundaryConditions(const unsigned int increment){
    constraints.clear (); constraints2.clear (); constraints3.clear ();
    constraints.reinit (locally_relevant_dofs);
    constraints2.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints2);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints3);

    //L2 boundary conditions
    constraints_L2.clear ();
    constraints_L2.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints_L2);
    constraints_L2.close ();
    
    //Setup boundary conditions
      std::vector<bool> uBCX0 (DIMS, false); uBCX0[3]=true; 
      //VectorTools::interpolate_boundary_values (dof_handler, 2, ConstantFunction<dim>(VBoundary, DIMS), constraints, uBCX0);
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints, uBCX0);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints, uBCX0);
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints2, uBCX0);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints2, uBCX0);
      //VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints);
      //VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(dim), constraints2);

      std::vector<bool> uBCX01 (DIMS, false); uBCX01[0]=true;
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints, uBCX01); 
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints2, uBCX01);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints, uBCX01);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints2, uBCX01);

      std::vector<bool> uBCX02 (DIMS, false); uBCX02[1]=true;
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints, uBCX02);
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints2, uBCX02);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints, uBCX02);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints2, uBCX02);

      std::vector<bool> uBCX03 (DIMS, false); uBCX03[2]=true;
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints, uBCX03);
      VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints2, uBCX03);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints, uBCX03);
      VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS), constraints2, uBCX03);
    constraints.close ();
    constraints2.close ();
    constraints3.close ();
  }
  
  //Setup
  template <int dim>
  void elasticity<dim>::setup_system (){
    TimerOutput::Scope t(computing_timer, "setup");
    dof_handler.distribute_dofs (fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    
    locally_relevant_solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    //Non-ghost vectors
    system_rhs.reinit (locally_owned_dofs, mpi_communicator);
    system_rhs_L2.reinit (locally_owned_dofs, mpi_communicator);
    U.reinit (locally_owned_dofs, mpi_communicator);
    Un.reinit (locally_owned_dofs, mpi_communicator);
    dU.reinit (locally_owned_dofs, mpi_communicator);
    //Ghost vectors
    UGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    
    
    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions(0);

    system_rhs_current.reinit (locally_owned_dofs, mpi_communicator);
    system_rhs_n.reinit (locally_owned_dofs, mpi_communicator);
    system_rhs_m.reinit (locally_owned_dofs, mpi_communicator);
    system_rhs_h.reinit (locally_owned_dofs, mpi_communicator);
    U_L2.reinit (locally_owned_dofs, mpi_communicator);
    U_current.reinit (locally_owned_dofs, mpi_communicator);
    U_n_iter.reinit (locally_owned_dofs, mpi_communicator);
    U_m_iter.reinit (locally_owned_dofs, mpi_communicator);
    U_h_iter.reinit (locally_owned_dofs, mpi_communicator);
    UGhost_L2.reinit (locally_owned_dofs, mpi_communicator);
    locally_relevant_solution_L2.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    DynamicSparsityPattern dsp_L2 (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp_L2, constraints_L2, false);
    SparsityTools::distribute_sparsity_pattern (dsp_L2, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    mass_matrix_L2.reinit (locally_owned_dofs, locally_owned_dofs, dsp_L2, mpi_communicator);

    
    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints2, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

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
    TimerOutput::Scope t(computing_timer, "assembly");
    system_rhs=0.0; system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(3);
    const QGauss<dim-1>	face_quadrature_formula (2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
     for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  if (std::abs(UGhost(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  else{ULocal[i]=UGhost(local_dof_indices[i]);}
	  ULocal[i].diff (i, dofs_per_cell);
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);
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
	residualForMechanics(fe_values, fe_face_values, 0, ULocal, ULocalConv, R, KMatrix, defMap, cell, increment, history_variables[cell], iteration);
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell; ++j){
	    // R' by AD
	    local_matrix(i,j)= R[i].fastAccessDx(j);
	  }
	  //R
	  local_rhs(i) = -R[i].val();
	}
	if ((currentIteration==0) && (currentIncrement==1) ){
	  constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
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
      TimerOutput::Scope t(computing_timer, "L2_projection");
    system_rhs_L2 =0.0;
    mass_matrix_L2=0.0;
    const QGauss<dim>  quadrature_formula(3);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   |
                             update_quadrature_points |
                             update_JxW_values);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix_L2 (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs_L2 (dofs_per_cell); 
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix_L2 = 0; local_rhs_L2 = 0; 
	cell->get_dof_indices (local_dof_indices);
	
	for (unsigned int q=0; q<n_q_points; q++){
	  for(unsigned int i=0;i<dofs_per_cell;i++){
	    const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first;
	    if (ci==0){
	      //local_rhs_current(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->current*fe_values.JxW(q);
	      local_rhs_L2(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->n_iter*fe_values.JxW(q);
	      //local_rhs_m(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->m_iter*fe_values.JxW(q);
	      //local_rhs_h(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->h_iter*fe_values.JxW(q);
	    }
	    else if (ci==1){
	      local_rhs_L2(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->m_iter*fe_values.JxW(q);
	    }
	     else if (ci==2){
	      local_rhs_L2(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->h_iter*fe_values.JxW(q);
	      }
	    else if (ci==3){
              local_rhs_L2(i) += fe_values.shape_value(i,q)*history_variables[cell][q]->h_iter*fe_values.JxW(q);
              }
	    for(unsigned int j=0;j<dofs_per_cell;j++){
	      const unsigned int cj = fe_values.get_fe().system_to_component_index(j).first;
	      if (ci==cj){
		local_matrix_L2(i,j) += fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q);
	      }
	    }
	  }
	}
	constraints_L2.distribute_local_to_global (local_matrix_L2, local_rhs_L2, local_dof_indices, mass_matrix_L2, system_rhs_L2);
	//constraints_L2.distribute_local_to_global (local_rhs_n, local_dof_indices, system_rhs_n);
	//constraints_L2.distribute_local_to_global (local_rhs_m, local_dof_indices, system_rhs_m);
	//constraints_L2.distribute_local_to_global (local_rhs_h, local_dof_indices, system_rhs_h);
      }
    mass_matrix_L2.compress (VectorOperation::add);
    system_rhs_L2.compress (VectorOperation::add);
    //system_rhs_n.compress (VectorOperation::add);
    //system_rhs_m.compress (VectorOperation::add);
    //system_rhs_h.compress (VectorOperation::add);
    }


  

  //Solve
  template <int dim>
  void elasticity<dim>::solveIteration(){
    
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
    */
    Timer time;
    time.start();
    //Direct solver MUMPS
    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    solver.set_symmetric_mode(false);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);
    //
    time.stop();
    pcout << "   Solved in " << time.last_wall_time()<<" seconds "<<std::endl;
    /*
    // Works great:
    SolverControl         solver_control(10000, 1e-8);
    SolverCG<Vector<double>> cg(solver_control);
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
    cg.solve(system_matrix, dU, system_rhs, preconditioner);
    */
    /*    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(dU, system_rhs);
    */
    if ( (currentIteration==0) && ((currentIncrement==1)) ){
      constraints.distribute (completely_distributed_solution);
      }
      //else if ((currentIncrement==1) || (currentIncrement==1000)){
      //constraints3.distribute (dU);
      //}
      else{
      constraints2.distribute (completely_distributed_solution);
      }
      
    locally_relevant_solution = completely_distributed_solution;
    dU = completely_distributed_solution; 
  }

   //Solve
  template <int dim>
  void elasticity<dim>::L2_solveIteration(){
    // Works great:
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    SolverControl cn1;
    PETScWrappers::SparseDirectMUMPS solver(cn1, mpi_communicator);
    solver.set_symmetric_mode(true);
    solver.solve(mass_matrix_L2, completely_distributed_solution, system_rhs_L2);
    constraints_L2.distribute (completely_distributed_solution);
    locally_relevant_solution_L2 = completely_distributed_solution;
    U_L2 = locally_relevant_solution_L2;
    UGhost_L2=U_L2;
    
    /*
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
    */
    }


  
  //Solve
  template <int dim>
  void elasticity<dim>::solve(const unsigned int increment){
    double res=1, tol=1.0e-8, abs_tol=1.0e-9, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];
    while (true){
      if (currentIteration>=4){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break; exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system(increment, currentIteration);
      current_norm=system_rhs.l2_norm();
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration();
      U+=dU; UGhost=U; 
      ++currentIteration;
    }
    Un=U; UnGhost=Un;
  }

  //Adaptive grid refinement
  template <int dim>
  void elasticity<dim>::refine_grid (){
    TimerOutput::Scope t(computing_timer, "adaptiveRefinement");
    
    bool checkForFurtherRefinement=true;
    while (checkForFurtherRefinement){
      bool isMeshRefined=false;
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();
      for (;cell!=endc; ++cell){
	if (cell->is_locally_owned()){
	  //limit the maximal and minimal refinement depth of the mesh
	  unsigned int current_level = t_cell->level();
	  
	  // Mark qPoins where refinement is to be done using bool.
	  bool mark_refine = false;
	  if (((cell->face(3)->center()[1]>1.6) && (1.6> cell->face(2)->center()[1]) )  ||  (cell->face(3)->center()[1]>=(0.5) ) && ( (0.5)>= cell->face(2)->center()[1])    ){
	      mark_refine=true; //set refine
	      
	  }

	  if ( (mark_refine && (current_level < 14))){
	    cell->set_refine_flag(); isMeshRefined=true; //refine
	  }
	}
	++t_cell;
      }

      //check for blocking in MPI                                                                                                                                                                                                                                               
      double checkSum=0.0;
      if (isMeshRefined){checkSum=1.0;}
      checkSum= Utilities::MPI::sum(checkSum, mpi_communicator);
      if (checkSum>0.0){
	
	parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector> soltrans(dof_handler);
	// prepare the triangulation,
	triangulation.prepare_coarsening_and_refinement();
	
	// prepare the SolutionTransfer object for coarsening and refinement
	// and give the solution vector that we intend to interpolate later,
	//Define a vector of vectors to store prev step  and prev to prev step ghosted variables
	std::vector<const LA::MPI::Vector*> InputGhosted(2);
	InputGhosted[0]=&UnGhost;	     InputGhosted[1]=&UGhost_L2;
	soltrans.prepare_for_coarsening_and_refinement(InputGhosted);  
	triangulation.execute_coarsening_and_refinement();
	setup_system(); 
	std::vector< LA::MPI::Vector*> tmp(2);
	tmp[0]=(&Un);
	tmp[1]=(&U_L2);
	
	soltrans.interpolate(tmp);
	
	UnGhost=(*tmp[0]);	     
	UGhost_L2=(*tmp[1]);
	UnGhost.update_ghost_values();
	UGhost_L2.update_ghost_values();
      }
      //set flag for another check of refinement
      checkForFurtherRefinement=false;
    }
  }
  
  
  //Output
  template <int dim>
  void elasticity<dim>::output_results (const unsigned int cycle) {
    TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (UnGhost, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
    //data_out.add_data_vector (Un, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
    data_out.add_data_vector (UGhost_L2, nodal_solution_names_L2, DataOut<dim>::type_dof_data, nodal_data_component_interpretation_L2);
    /*data_out.add_data_vector (U_n_iter, "n_mag");
    data_out.add_data_vector (U_m_iter, "m_mag");
    data_out.add_data_vector (U_h_iter, "h_mag");
    */
    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
      data_out.add_data_vector (subdomain, "subdomain"); 
      //std::cout<<"Stage 8"<<std::endl;
    data_out.build_patches ();
    //std::cout<<"Stage 9"<<std::endl;
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2) +
				                                "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      std::vector<std::string> filenames;
      for (unsigned int i=0;
	   i<Utilities::MPI::n_mpi_processes(mpi_communicator);
	   ++i)
	filenames.push_back ("solution-" +
			     Utilities::int_to_string (cycle, 2) +
			     "." +
			     Utilities::int_to_string (i, 4) +
			     ".vtu");
      
      std::ofstream master_output (("solution-" +
				    Utilities::int_to_string (cycle, 2) +
				    ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
      }
  }

  //Solve problem
  template <int dim>
  void elasticity<dim>::run (){
    //setup problem geometry and mesh
      std::vector<unsigned int> n_subdivisions;
      n_subdivisions.push_back(64*5);//128*12 //1500, 100
      n_subdivisions.push_back(32); //128
      GridGenerator::subdivided_hyper_rectangle(triangulation, n_subdivisions,
						(dim==3 ? Point<dim>(0.0, 0.0, 0.0) : Point<dim>(0.0, 0.0)),
						(dim==3 ? Point<dim>(5.0, 1.0, 1.0) : Point<dim>(10000.0, 2.1)),
						true);
  

    setup_system();
    refine_grid();
    refine_grid();
    refine_grid();
    pcout << "   Number of active cells:       "
	  << triangulation.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;
    //setup initial conditions
    VectorTools::interpolate(dof_handler, InitalConditions<dim>(), U); Un=U;
    //sync ghost vectors to non-ghost vectors
    UGhost=U;  UnGhost=Un;
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
      pcout << std::endl;
    }
    //computing_timer.print_summary ();
  }
}


int main(int argc, char *argv[]){
  try
    {
      using namespace dealii;
      using namespace elasticity1;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
