//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2011
//authors: Rahul, Rudraa (2011, 2018)
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
  alpha[0]= 0.07*std::exp((-V)/20.0); //h, m, n 
  if (std::abs(V-25.0)<1.0e-5){
    alpha[1]=0.0;
  }
  else{
    alpha[1]= (2.5-0.1*(V))/(std::exp(2.5-0.1*(V)) - 1.0);
  }
  if (std::abs(V-10.0)<1.0e-5){
    alpha[2]=0.0;
  }
  else{
    alpha[2]= (0.1-0.01*(V))/(std::exp(1.0 -0.1*(V)) - 1.0);
  } 
  beta[0]=  1.0/(1.0+std::exp(3.0-0.1*(V)));
  beta[1]=  4.0*std::exp((-V)/18.0);
  beta[2]=  0.125*std::exp((-V)/80.0);
}

//Mechanics implementation
template <class T, int dim>
void evaluateStress(FEValues<dim>& fe_values,unsigned int DOF, Table<1, T>& ULocal, Table<1, double>& ULocalConv, Table<2, Sacado::Fad::DFad<double>>& Term11, Table<1, Sacado::Fad::DFad<double>>& Term12, Table<1, Sacado::Fad::DFad<double>>& Term13, double& Term14, Table<2, Sacado::Fad::DFad<double>>& Term21,Table<1, Sacado::Fad::DFad<double>>& Term22,  double& Term23, Table<1, Sacado::Fad::DFad<double>>& Term24, deformationMap<T, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator& cell, const unsigned int currentIncrement, std::vector<history<dim>*>& history_variables, const unsigned int currentIteration){

  //number of quadrature points
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    if ( currentIteration == 0 ){
      history_variables[q]->n_inc = history_variables[q]->n_iter;
      history_variables[q]->h_inc = history_variables[q]->h_iter; 
      history_variables[q]->m_inc = history_variables[q]->m_iter;   
    }
        
    //Get V_A, V_AConv
    Table<1, Sacado::Fad::DFad<double>> V_A (n_q_points);
    Table<1, double> V_AConv(n_q_points);
    Table<2, Sacado::Fad::DFad<double>> gradV_A(n_q_points, dim);
    for(unsigned int k=0; k<dim;++k){
      gradV_A[q][k]=0.0;
    }
    V_A[q]=0.0; V_AConv[q]=0.0;
    for (unsigned int k=0; k<fe_values.dofs_per_cell; ++k){
      if (fe_values.get_fe().system_to_component_index(k).first==0){
	V_A[q]+=ULocal[k]*fe_values.shape_value(k, q); //U
	V_AConv[q] += ULocalConv[k]*fe_values.shape_value(k, q); //U
	for (unsigned int i=0; i<dim; ++i){
	  gradV_A[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	}
      }
    }

    
    //Get V_M, V_MConv
    Table<1, Sacado::Fad::DFad<double>> V_M (n_q_points);
    Table<1, double> V_MConv(n_q_points);
    Table<2, Sacado::Fad::DFad<double>> gradV_M(n_q_points, dim);
    for(unsigned int k=0; k<dim;++k){
      gradV_M[q][k]=0.0;
    }
    V_M[q]=0.0; V_MConv[q]=0.0;
    for (unsigned int k=0; k<fe_values.dofs_per_cell; ++k){
      if (fe_values.get_fe().system_to_component_index(k).first==1){
	V_M[q]+=ULocal[k]*fe_values.shape_value(k, q); //U
	V_MConv[q] += ULocalConv[k]*fe_values.shape_value(k, q); //U
	for (unsigned int i=0; i<dim; ++i){
	  gradV_M[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	}
      }
    }
    
    //Gating Constants
    double V_value;
    if (dim==1){
      V_value= V_AConv[q]-V_rest; 
    }
    else{
      V_value= V_AConv[q];
    }
    Table<1, double> alpha (3), beta(3);
    evaluateGatingConstants<dim>(V_value, alpha, beta);
    history_variables[q]->n_iter =  history_variables[q]->n_inc + DT*(alpha[2]*(1.0-history_variables[q]->n_inc)-beta[2]*history_variables[q]->n_inc); 
    history_variables[q]->m_iter =  history_variables[q]->m_inc + DT*(alpha[1]*(1.0-history_variables[q]->m_inc)-beta[1]*history_variables[q]->m_inc); 
    history_variables[q]->h_iter =  history_variables[q]->h_inc + DT*(alpha[0]*(1.0-history_variables[q]->h_inc)-beta[0]*history_variables[q]->h_inc); 
    
    //Residual_Terms
    Point<dim> point1;
    point1=fe_values.quadrature_point(q);
    double g_Na1=g_Na, g_K1=g_K, g_L1=g_L,C_m1=C_m, C_my1= C_my, E_L1=E_L,R_i1=R_i,R_my1=R_my, R_m1=R_m, Inj=0.0, rPA1=rPA;
    if (point1[0]<0.1) {
      Inj=0.0;
    }
    
    //R_my1=0.00000001;
    if (Myelin){ 
      if ((std::abs(std::remainder(point1[0]*10000.0,70.0)-0.0)< 2.3)) {
	R_my1=0.00000001; rPA1=rPN;
      }
    }
    
    Table<2, Sacado::Fad::DFad<double>> gradV_Total(n_q_points, dim);
    for(unsigned int k=0; k<dim; ++k){
      gradV_Total[q][k] = gradV_A[q][k] + gradV_M[q][k];
    }
    Term14=(RADIUS/(2.0*R_i1));
    for(unsigned int k=0; k<dim; ++k){
      Term11[q][k]= Term14*(gradV_Total[q][k]);
    }
    Term12[q]= (C_m1/DT)*(V_A[q]-V_AConv[q]) +  g_K1*(std::pow(history_variables[q]->n_iter,4))*(V_A[q]-E_K) + g_Na1*(std::pow(history_variables[q]->m_iter,3))*history_variables[q]->h_iter*(V_A[q]-E_Na) + g_L1*(V_A[q]-E_L1) -Inj; // + (V_A[q]); 
    Term13[q]=(C_m1/DT)+  g_K1*(std::pow(history_variables[q]->n_iter,4)) + g_Na1*(std::pow(history_variables[q]->m_iter,3))*history_variables[q]->h_iter + g_L1;// +1.0;


    
    Term24[q] = (RADIUS*RADIUS/(2.0*R_i1*periRadius)) + (1.0/(2.0*3.14*periRadius*rPA1));
    for(unsigned int k=0; k<dim; ++k){
      Term21[q][k] = Term24[q]*gradV_Total[q][k] ;
    }
    
    Term22[q] = V_M[q]*(1.0/R_my1) + (C_my1/DT)*(V_M[q]-V_MConv[q]);
    Term23= (C_my1/DT) + (1.0/R_my1);
  
    
    //Current for L2 projection
    history_variables[q]->current = 0.0; //C_m1*(Volt[q].val()-V_conv[q])*(1.0/DT) + g_K1*(Volt[q].val()-E_K)*(std::pow(history_variables[q]->n_iter,4.0))+g_Na1*(Volt[q].val() - E_Na)*(std::pow(history_variables[q]->m_iter,3))*history_variables[q]->h_iter+ g_L1*(Volt[q].val()-E_L1);
    //std::cout<<"V_A="<<V_A[q].val()<<" ,V_AConv="<<V_AConv[q]<<" ,gradV_A="<<gradV_A[q][0].val()<<" ,V_M="<<V_M[q].val()<<" ,V_MConv="<<V_MConv[q]<<" ,gradV_M="<<gradV_M[q][0].val()<<std::endl;
  }
}

//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double>& ULocalConv, Table<1, Sacado::Fad::DFad<double> >& R, Table<2, double>& KMatrix, deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator& cell, const unsigned int currentIncrement, std::vector<history<dim>*>& history_variables, const unsigned int currentIteration){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  
  //temporary arrays
  //Table<1, Sacado::Fad::DFad<double> > Volt (n_q_points);
  Table<2, Sacado::Fad::DFad<double>> Term11 (n_q_points, dim);
  Table<1, Sacado::Fad::DFad<double>> Term12 (n_q_points);
  Table<1, Sacado::Fad::DFad<double>> Term13 (n_q_points);
  double Term14;
  Table<1, Sacado::Fad::DFad<double>> Term22 (n_q_points);
  Table<2, Sacado::Fad::DFad<double>> Term21 (n_q_points, dim);
  double Term23;
  Table<1, Sacado::Fad::DFad<double>> Term24(n_q_points);
  
  //evaluate stress
  evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, DOF, ULocal, ULocalConv, Term11, Term12, Term13, Term14, Term21, Term22, Term23, Term24, defMap, cell, currentIncrement, history_variables, currentIteration);

  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
	for (unsigned int d = 0; d < dim; d++){
	  R[i] += Term11[q][d]*fe_values.shape_grad(i,q)[d]*fe_values.JxW(q);
	}
	R[i] += Term12[q]*fe_values.shape_value(i,q)*fe_values.JxW(q);
      }
      else if (ck==1){
	R[i] += fe_values.shape_value(i,q)*Term22[q]*fe_values.JxW(q);
	for (unsigned int d = 0; d < dim; d++){
	  R[i] += Term21[q][d]*fe_values.shape_grad(i,q)[d]*fe_values.JxW(q);
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
	  KMatrix[i][j] += Term13[q].val()*fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q);
	  for(unsigned int d=0; d<dim; ++d){
	    KMatrix[i][j] += Term14*fe_values.shape_grad(i,q)[d]*fe_values.shape_grad(j,q)[d]*fe_values.JxW(q);
	  }
	}
	  else if (ck ==0 && cj==1){
	    for(unsigned int d=0; d<dim; ++d){
	      KMatrix[i][j] += Term14*fe_values.shape_grad(i,q)[d]*fe_values.shape_grad(j,q)[d]*fe_values.JxW(q);
	    } 
	  }
	else if (ck==1 && cj==ck){
	  KMatrix[i][j] += Term23*fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q);
	   for(unsigned int d=0; d<dim; ++d){
	     KMatrix[i][j] += Term24[q].val()*fe_values.shape_grad(i,q)[d]*fe_values.shape_grad(j,q)[d]*fe_values.JxW(q);
	  }
	}
	 else if (ck ==1 && cj==0){
	    for(unsigned int d=0; d<dim; ++d){
	      KMatrix[i][j] += Term24[q].val()*fe_values.shape_grad(i,q)[d]*fe_values.shape_grad(j,q)[d]*fe_values.JxW(q);
	    } 
	  }
      }
    }
  }
}

#endif /* MECHANICS_H_ */
