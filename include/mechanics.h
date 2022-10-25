//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2021
//authors: Rahul Gulati (2021)
//
#ifndef MECHANICS_H_
#define MECHANICS_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "deformationMap.h"

//history variables
template <int dim>
class history{
public:
  history() : n_iter( (0.1/(std::exp(1.0) -1.0))/( (0.1/(std::exp(1.0) -1.0)) +0.125) ), h_iter( 0.07/(0.07 + (1.0/(1.0+std::exp(3.0))) )  ), m_iter( (2.5/(std::exp(2.5)-1.0))/((2.5/(std::exp(2.5)-1.0))+4.0)   ), n_inc( (0.1/(std::exp(1.0) -1.0))/( (0.1/(std::exp(1.0) -1.0)) +0.125)), h_inc( 0.07/(0.07 + (1.0/(1.0+std::exp(3.0))) )  ), m_inc(  (2.5/(std::exp(2.5)-1.0))/((2.5/(std::exp(2.5)-1.0))+4.0)  ), current(0.0) {};
  double n_iter,n_inc, h_iter,h_inc, m_iter, m_inc, current;
};

template<int dim>
void evaluateGatingConstants(double& V, Table<1,double>& alpha, Table<1,double>& beta){
  alpha[0]= 0.07*std::exp((V_rest-V)/20.0); //h, m, n
  alpha[1]= (2.5-0.1*(V-V_rest))/(std::exp(2.5-0.1*(V-V_rest)) - 1.0);
  alpha[2]= (0.1-0.01*(V-V_rest))/(std::exp(1.0 -0.1*(V-V_rest)) - 1.0);
  beta[0]=  1.0/(1.0+std::exp(3.0-0.1*(V-V_rest)));
  beta[1]=  4.0*std::exp((V_rest-V)/18.0);
  beta[2]=  0.125*std::exp((V_rest-V)/80.0);
}

//Mechanics implementation
template <class T, int dim>
void evaluateStress(FEValues<dim>& fe_values,unsigned int DOF, Table<1, T>& ULocal, Table<1, double>& ULocalConv, Table<1, T>& Volt, Table<1, double>& Term2, Table<2, T>& dVolt,Table<1, double>& ConstantTerm, deformationMap<T, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator& cell, const unsigned int currentIncrement, std::vector<history<dim>*>& history_variables, const unsigned int currentIteration, Table<1, Sacado::Fad::DFad<double>>& phi){

  //number of quadrature points
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    if ( currentIteration == 0 ){
      history_variables[q]->n_inc = history_variables[q]->n_iter;
      history_variables[q]->h_inc = history_variables[q]->h_iter; 
      history_variables[q]->m_inc = history_variables[q]->m_iter;   
    }
        

    Table<2, T> gradV(n_q_points, dim);
    evaluateScalarFunctionGradient<T, dim>(fe_values, DOF, ULocal, gradV, defMap, false);
    Table<1, T> V_T(n_q_points);
    evaluateScalarFunction<T, dim>(fe_values, DOF, ULocal, V_T);
    Table<1, double> V_ConvT(n_q_points);
    evaluateScalarFunction<double, dim>(fe_values, DOF, ULocalConv, V_ConvT);

    //Gating Constants
    double V_value= V_ConvT[q]; //V_T[q].val();
    Table<1, double> alpha (3), beta(3);
    evaluateGatingConstants<dim>(V_value, alpha, beta);
	history_variables[q]->n_iter =  history_variables[q]->n_inc + DT*(alpha[2]*(1.0-history_variables[q]->n_inc)-beta[2]*history_variables[q]->n_inc); 
	history_variables[q]->m_iter =  history_variables[q]->m_inc + DT*(alpha[1]*(1.0-history_variables[q]->m_inc)-beta[1]*history_variables[q]->m_inc); 
	history_variables[q]->h_iter =  history_variables[q]->h_inc + DT*(alpha[0]*(1.0-history_variables[q]->h_inc)-beta[0]*history_variables[q]->h_inc); 
	

    Table<1, double> V_conv(n_q_points);
    for (unsigned int i=0; i<dim; ++i){
	dVolt[q][i] = gradV[q][i];
    }
	Volt[q] = V_T[q];
	V_conv[q]= V_ConvT[q];

    //Residual_Terms
    Point<dim> point1;
    point1=fe_values.quadrature_point(q);
    double g_Na1=g_Na, g_K1=g_K, g_L1=g_L,C_m1=C_m;
    if (Myelin){
      g_Na1=0.0; g_K1=0.0; g_L1=1.0e-6/33.0; C_m1=MC_m;
    if ((std::abs(std::remainder(point1[0],2.0)-0.0)< 0.1) && (point1[0]>1.5)) {
    g_Na1=g_Na; g_K1=g_K; g_L1=g_L; C_m1=C_m;
    }
    }
      
    ConstantTerm[q]= (C_m1/DT) + g_K1*(std::pow(history_variables[q]->n_iter,4)) + g_Na1*(std::pow(history_variables[q]->m_iter,3))*history_variables[q]->h_iter + g_L1;

  Term2[q] = V_conv[q]*(C_m1/DT) + g_K1*(std::pow(history_variables[q]->n_iter,4.0))*E_K + g_Na1*(std::pow(history_variables[q]->m_iter,3))*history_variables[q]->h_iter*E_Na + g_L1*E_L;
    
    //Current for L2 projection
    history_variables[q]->current = C_m1*(Volt[q].val()-V_conv[q])*(1.0/DT) + g_K1*(Volt[q].val()-E_K)*(std::pow(history_variables[q]->n_iter,4.0))+g_Na1*(Volt[q].val() - E_Na)*(std::pow(history_variables[q]->m_iter,3))*history_variables[q]->h_iter+ g_L1*(Volt[q].val()-E_L);
    for (unsigned int index1=0; index1<1; ++index1){
      phi[q]=0.0;
      for (unsigned int k=0; k<fe_values.dofs_per_cell; ++k){
	if (fe_values.get_fe().system_to_component_index(k).first==1){
	  phi[q]+=ULocal[k]*fe_values.shape_value(k, q); //U
	}
      }
    }


    
  }
}

//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double> >& ULocal, Table<1, double>& ULocalConv, Table<1, Sacado::Fad::DFad<double> >& R, Table<2, double>& KMatrix, deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator& cell, const unsigned int currentIncrement, std::vector<history<dim>*>& history_variables, const unsigned int currentIteration){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  
  //temporary arrays
  Table<1, Sacado::Fad::DFad<double> > Volt (n_q_points);
  Table<1, double> Term2 (n_q_points);
  Table<2, Sacado::Fad::DFad<double>> dVolt (n_q_points, dim);
  Table<1, double> diffusivity (dim);
  diffusivity[0]= 1.0*a_RADIUS*(1.0/(2.0*R_SP_RESIST));
  if (dim>1){
    diffusivity[1]= diffusivity[0];
  }
  if (dim>2){
    diffusivity[2]= diffusivity[0];
  }
  Table<1, double> Constant_Term (n_q_points);
  Table<1, Sacado::Fad::DFad<double>> phi(n_q_points);
  evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, DOF, ULocal, ULocalConv, Volt, Term2, dVolt, Constant_Term, defMap, cell, currentIncrement, history_variables, currentIteration, phi);
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      //if (ck>=0 && ck<dim){
      if (ck==0){
    
	for (unsigned int d = 0; d < dim; d++){
	  R[i] += diffusivity[d]*dVolt[q][d]*fe_values.shape_grad(i,q)[d]*fe_values.JxW(q);
	}
	R[i] += 1.0*(Volt[q]*Constant_Term[q] - Term2[q])*fe_values.shape_value(i,q)*fe_values.JxW(q);
      }
      else if (ck==1){
	R[i] += fe_values.shape_value(i,q)*phi[q]*fe_values.JxW(q);
	for (unsigned int d = 0; d < dim; d++){
	  R[i] += (a_RADIUS/(2.0*R_SP_RESIST))*dVolt[q][d]*fe_values.shape_value(i,q)*fe_values.JxW(q);
	}
      }
    }
  }
  
  //KMatrix
  for(unsigned int i=0; i<dofs_per_cell; ++i){
    for(unsigned int j=0; j<dofs_per_cell; ++j){
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
      const unsigned int cj = fe_values.get_fe().system_to_component_index(j).first - DOF;
      for (unsigned int q=0; q<n_q_points; ++q){
	if (ck ==0 && cj==ck){ 

	  for(unsigned int d=0; d < dim; d++){
	    KMatrix[i][j] += diffusivity[d]*fe_values.shape_grad(i,q)[d]*fe_values.shape_grad(j,q)[d]*fe_values.JxW(q);
	  }
	  KMatrix[i][j] += 1.0*Constant_Term[q]*fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q);
	}
	else if (ck==1 && cj==ck){
	  KMatrix[i][j] += fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q);
	}
      }
    }
  }
}

#endif /* MECHANICS_H_ */
