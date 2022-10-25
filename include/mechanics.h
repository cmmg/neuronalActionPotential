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
void evaluateGatingConstants(double V, Table<1,double>& alpha, Table<1,double>& beta){
  alpha[0]= 0.07*std::exp((V_rest-V)/20.0); //h, m, n                                                                                                                                      
  alpha[1]= (2.5-0.1*(V-V_rest))/(std::exp(2.5-0.1*(V-V_rest)) - 1.0);
  alpha[2]= (0.1-0.01*(V-V_rest))/(std::exp(1.0 -0.1*(V-V_rest)) - 1.0);
  beta[0]=  1.0/(1.0+std::exp(3.0-0.1*(V-V_rest)));
  beta[1]=  4.0*std::exp((V_rest-V)/18.0);
  beta[2]=  0.125*std::exp((V_rest-V)/80.0);
}

//Mechanics implementation
template <class T, int dim>
void evaluateStress(FEValues<dim>& fe_values, FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<1, double>& ULocalConv, Table<1, Sacado::Fad::DFad<double>>& Term11, Table<2, Sacado::Fad::DFad<double>>& Term12, Table<1, Sacado::Fad::DFad<double>>& Term13, Table<2, Sacado::Fad::DFad<double>>& Term14, Table<1, Sacado::Fad::DFad<double>>& Term21, Table<2, Sacado::Fad::DFad<double>>& Term22, Table<1, Sacado::Fad::DFad<double>>&  Term23, Table<2, Sacado::Fad::DFad<double>>& Term24, Table<1, Sacado::Fad::DFad<double>>& Term31, Table<2, Sacado::Fad::DFad<double>>& Term32, Table<1, Sacado::Fad::DFad<double>>& Term33, Table<2, Sacado::Fad::DFad<double>>& Term34, Table<1, Sacado::Fad::DFad<double>>& Term41, Table<2, Sacado::Fad::DFad<double>>& Term42, Table<1, Sacado::Fad::DFad<double>>& Term50, Table<1, Sacado::Fad::DFad<double>>& Term51, Table<1, Sacado::Fad::DFad<double>>& Term52, deformationMap<T, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator& cell, const unsigned int currentIncrement, std::vector<history<dim>*>& history_variables, const unsigned int currentIteration){
  //number of quadrature points
  unsigned int n_q_points= fe_values.n_quadrature_points;
   Table<1, Sacado::Fad::DFad<double>> concNa (n_q_points), concK (n_q_points), concA (n_q_points), volt (n_q_points);
   Table<1, double> concConvNa (n_q_points), concConvK (n_q_points), concConvA (n_q_points), voltConv (n_q_points);
   Table<2, Sacado::Fad::DFad<double>> gradConcNa (n_q_points, dim), gradConcK (n_q_points, dim), gradConcA (n_q_points, dim), gradVolt (n_q_points, dim);
  //loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    if ( currentIteration == 0 ){
      history_variables[q]->n_inc = history_variables[q]->n_iter;
      history_variables[q]->h_inc = history_variables[q]->h_iter;
      history_variables[q]->m_inc = history_variables[q]->m_iter;   
    }
   
    for(unsigned int k=0; k<dim;++k){
      gradConcNa[q][k]=0.0;
      gradConcK[q][k]=0.0;
      gradConcA[q][k]=0.0;
      gradVolt[q][k]=0.0;
    }
    concNa[q]=0.0; concConvNa[q]=0.0;
    concK[q]=0.0; concConvK[q]=0.0;
    concA[q]=0.0; concConvA[q]=0.0;
    volt[q] = 0.0; voltConv[q] = 0.0;
    for (unsigned int k=0; k<fe_values.dofs_per_cell; ++k){
      if (fe_values.get_fe().system_to_component_index(k).first==0){
	concNa[q]+=ULocal[k]*fe_values.shape_value(k, q); //U
	concConvNa[q] += ULocalConv[k]*fe_values.shape_value(k, q); //U
	for (unsigned int i=0; i<dim; ++i){
	  gradConcNa[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	}
      }
      if (fe_values.get_fe().system_to_component_index(k).first==1){
        concK[q]+=ULocal[k]*fe_values.shape_value(k, q); //U 
        concConvK[q] += ULocalConv[k]*fe_values.shape_value(k, q); //U 
        for (unsigned int i=0; i<dim; ++i){
          gradConcK[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU  
        }
      }
      if (fe_values.get_fe().system_to_component_index(k).first==2){
        concA[q]+=ULocal[k]*fe_values.shape_value(k, q); //U 
        concConvA[q] += ULocalConv[k]*fe_values.shape_value(k, q); //U 
        for (unsigned int i=0; i<dim; ++i){
          gradConcA[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU 
        }
      }
      if (fe_values.get_fe().system_to_component_index(k).first==3){
        volt[q]+=ULocal[k]*fe_values.shape_value(k, q); //U  
        voltConv[q] += ULocalConv[k]*fe_values.shape_value(k, q); //U  
        for (unsigned int i=0; i<dim; ++i){
          gradVolt[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU 
        }
      }
    }

    Point<dim> point1;
    point1=fe_values.quadrature_point(q);
    unsigned int faceID, faceID1=2;
    if ((cell->face(3)->center()[1]>(0.5) ) && ((0.5)> cell->face(2)->center()[1])){
      faceID = 3;
    }
    else{
      faceID = 2;
    }
    fe_face_values.reinit (cell, faceID);
    unsigned int n_face_q_points= fe_face_values.n_quadrature_points;
    Table<1, Sacado::Fad::DFad<double> > voltFace(n_face_q_points), voltFaceConv(n_face_q_points);
    unsigned int face_dofs_per_cell= fe_face_values.dofs_per_cell;
    for (unsigned int q_face=0; q_face<n_face_q_points; ++q_face){
      voltFace[q_face]=0.0; voltFaceConv[q_face]=0.0;
      for (unsigned int k=0; k<face_dofs_per_cell; ++k){
        unsigned int ck = fe_face_values.get_fe().system_to_component_index(k).first;
        if (ck==3){
          voltFace[q_face]+=ULocal[k]*fe_face_values.shape_value(k, q_face);
	  voltFaceConv[q_face]+=ULocalConv[k]*fe_face_values.shape_value(k, q_face);
        }
      }
    }
    //get cNa & cK in & out for the Lower membrane Interface                                                                                                                                                    
    if ((cell->face(3)->center()[1]>(0.5) ) && ((0.5)> cell->face(2)->center()[1])){
      faceID = 2; faceID1 = 3;
    }
    else if  ((cell->face(3)->center()[1]>(1.6) ) && ((1.6)> cell->face(2)->center()[1])){
      faceID = 3; faceID1 =2;
    }
    fe_face_values.reinit (cell, faceID);
    n_face_q_points= fe_face_values.n_quadrature_points;
    Table<1, Sacado::Fad::DFad<double> > cNaOut(n_face_q_points), cKOut(n_face_q_points);
    face_dofs_per_cell= fe_face_values.dofs_per_cell;
    for (unsigned int q_face=0; q_face<n_face_q_points; ++q_face){
      cNaOut[q_face]=0.0; cKOut[q_face]=0.0;
      for (unsigned int k=0; k<face_dofs_per_cell; ++k){
        unsigned int ck = fe_face_values.get_fe().system_to_component_index(k).first;
        if (ck==0){
          cNaOut[q_face] += ULocal[k]*fe_face_values.shape_value(k, q_face);
          //voltFaceConv[q_face]+=ULocalConv[k]*fe_face_values.shape_value(k, q_face);                                                                                                                          
        }
        if (ck==1){
          cKOut[q_face] += ULocal[k]*fe_face_values.shape_value(k, q_face);
        }

      }
    }
    fe_face_values.reinit (cell, faceID1);
    n_face_q_points= fe_face_values.n_quadrature_points;
    Table<1, Sacado::Fad::DFad<double> > cNaIn(n_face_q_points), cKIn(n_face_q_points);
    face_dofs_per_cell= fe_face_values.dofs_per_cell;
    for (unsigned int q_face=0; q_face<n_face_q_points; ++q_face){
      cNaIn[q_face]=0.0; cKIn[q_face]=0.0;
      for (unsigned int k=0; k<face_dofs_per_cell; ++k){
        unsigned int ck = fe_face_values.get_fe().system_to_component_index(k).first;
        if (ck==0){
          cNaIn[q_face] += ULocal[k]*fe_face_values.shape_value(k, q_face);
          //voltFaceConv[q_face]+=ULocalConv[k]*fe_face_values.shape_value(k, q_face);                                                                                                                          
        }
        if (ck==1){
          cKIn[q_face] += ULocal[k]*fe_face_values.shape_value(k, q_face);
        }

      }
    }
    Table<1, double> alpha (3), beta(3), inf(3), tau(3);
    evaluateGatingConstants<dim>(voltFace[0].val(), alpha, beta);
    for (unsigned int k=0; k<3; ++k){
      tau[k]= 1.0 /(alpha[k]+beta[k]);
      inf[k]= alpha[k]/(alpha[k]+beta[k]);
    }

    history_variables[q]->n_iter = inf[2] -(inf[2]-history_variables[q]->n_inc)*std::exp(-DT/tau[2]);
    history_variables[q]->m_iter = inf[1] -(inf[1]-history_variables[q]->m_inc)*std::exp(-DT/tau[1]);                                                                                                                           
    history_variables[q]->h_iter = inf[0] -(inf[0]-history_variables[q]->h_inc)*std::exp(-DT/tau[0]);
    
    double memDist, capacitance, myelinCap;
    if ((cell->face(3)->center()[1]>=(1.6)) && ( (1.6)>= cell->face(2)->center()[1]) && (cell->face(2)->center()[0]>= (6.2) ) && (cell->face(2)->center()[0]<=(6.3) ) ){
      memDist = cell->face(3)->center()[1] - cell->face(2)->center()[1];
      capacitance = 8.88541*memRelativePerm/(10000.0*memDist);
      myelinCap = 8.88541*memRelativePerm/(10000.0*(dimH-0.016));
      std::cout<<"Point="<<point1[0]<<","<<point1[1]<<" , voltFace="<<voltFace[0].val()<<" ,faceid="<<faceID<<" ,n="<<history_variables[q]->n_iter<<", m="<<history_variables[q]->m_iter<<", h="<<history_variables[q]->h_iter<<" ,memDist="<<memDist<<" ,capacitance="<<capacitance<<" ,myelinCapacitance="<<myelinCap<<std::endl;
    }
    if ((cell->face(3)->center()[1]>=(0.5)) && ( (0.5)>= cell->face(2)->center()[1]) && (cell->face(2)->center()[0]>= (6.2) ) && (cell->face(2)->center()[0]<=(6.3) ) ){
      memDist = cell->face(3)->center()[1] - cell->face(2)->center()[1];
      capacitance = 8.88541*memRelativePerm/(10000.0*memDist);
      std::cout<<"Point="<<point1[0]<<","<<point1[1]<<" , voltFace="<<voltFace[0].val()<<" ,faceid="<<faceID<<" ,n="<<history_variables[q]->n_iter<<", m="<<history_variables[q]->m_iter<<", h="<<history_variables[q]->h_iter<<" ,memDist="<<memDist<<" ,capacitance="<<capacitance<<std::endl;
      }
    
    //Residual_Terms
    double gNa1=gNa, gK1=gK;

    double diffusivityNa = diffusivityConstNa, diffusivityK = diffusivityConstK, diffusivityA = diffusivityConstA, permittivity = permittivityW;
   
    double start = 5.0, dimb= dimB;

    if( ((point1[1]>1.6 && point1[1] < (1.6+dimH) ) || (point1[1]>(0.5-dimH) && point1[1] < 0.5)) && point1[0]<start ){
      diffusivityNa = 0.0; diffusivityK = 0.0; diffusivityA = 0.0; permittivity = permittivityM;
    } 


    while (start<neuronLength){
      if( ((point1[1]>1.6 && point1[1] < (1.6+dimH) ) || (point1[1]>(0.5-dimH) && point1[1] < 0.5)) && point1[0]>(start+dimb) && point1[0]< (start+dimL) ){
        diffusivityNa = 0.0; diffusivityK = 0.0; diffusivityA = 0.0; permittivity = permittivityM;
      }
      start +=dimL;
    }

    if( ((point1[1]>1.6 && point1[1] < (1.6+dimH) ) || (point1[1]>(0.5-dimH) && point1[1] < 0.5))  && point1[0]>(start+dimb) ){
      diffusivityNa = 0.0; diffusivityK = 0.0; diffusivityA = 0.0; permittivity = permittivityM;
    } 

    if ((cell->face(3)->center()[1]>=(1.6)) && ( (1.6)>= cell->face(2)->center()[1])){
       diffusivityNa = 0.0; diffusivityK = 0.0; diffusivityA = 0.0; permittivity = permittivityM;
    }

    if ((cell->face(3)->center()[1]>=(0.5)) && ( (0.5)>= cell->face(2)->center()[1])){
    diffusivityNa = 0.0; diffusivityK = 0.0; diffusivityA = 0.0; permittivity = permittivityM;
    }

    //DC network
    double start1, dimb1= dimB;
    start1= 5.0+dimB+5.0;

    while(start1<neuronLength){
      if( (  ( cell->face(3)->center()[1]>=(1.6) && cell->face(3)->center()[1]<=(1.617)  ) || ( cell->face(2)->center()[1]<=(0.5) && cell->face(2)->center()[1]>=(0.483)  )) && point1[0]>start1 && point1[0]< (start1+dimb1) ){
	diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }

      if( (  ( cell->face(3)->center()[1]>=(1.609) && cell->face(2)->center()[1]<=(1.609)  ) || (cell->face(3)->center()[1]>=(0.491) && cell->face(2)->center()[1]<=(0.491)) )  && point1[0]> start1 && point1[0]< (start1+dimL-3*dimb1) ){
	diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }

      if( (  ( cell->face(3)->center()[1]>=(1.611) && cell->face(2)->center()[1]<=(1.611)  ) || (cell->face(3)->center()[1]>=(0.489) && cell->face(2)->center()[1]<=(0.489)) )  && point1[0]> start1 && point1[0]< (start1+dimL-3*dimb1) ){
        diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }

      if( (  ( cell->face(3)->center()[1]>=(1.607) && cell->face(2)->center()[1]<=(1.607)  ) || (cell->face(3)->center()[1]>=(0.493) && cell->face(2)->center()[1]<=(0.493)) )  && point1[0]> start1 && point1[0]< (start1+dimL-3*dimb1) ){
        diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }
      if( (  ( cell->face(3)->center()[1]>=(1.605) && cell->face(2)->center()[1]<=(1.605)  ) || (cell->face(3)->center()[1]>=(0.495) && cell->face(2)->center()[1]<=(0.495)) )  && point1[0]> start1 && point1[0]< (start1+dimL-3*dimb1) ){
        diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }
      if( (  ( cell->face(3)->center()[1]>=(1.613) && cell->face(2)->center()[1]<=(1.613)  ) || (cell->face(3)->center()[1]>=(0.487) && cell->face(2)->center()[1]<=(0.487)) )  && point1[0]> start1 && point1[0]< (start1+dimL-3*dimb1) ){
        diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }
      if( (  ( cell->face(3)->center()[1]>=(1.615) && cell->face(2)->center()[1]<=(1.615)  ) || (cell->face(3)->center()[1]>=(0.485) && cell->face(2)->center()[1]<=(0.485)) )  && point1[0]> start1 && point1[0]< (start1+dimL-3*dimb1) ){
        diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }


      if( (  ( cell->face(3)->center()[1]>=(1.6) && cell->face(3)->center()[1]<=(1.617)  ) || ( cell->face(2)->center()[1]<=(0.5) && cell->face(2)->center()[1]>=(0.483)  )  ) && point1[0]> (start1+dimL-4.0*dimb1) && point1[0]< (start1+dimL-3.0*dimb1) ){
        diffusivityNa = diffusivityConstNa; diffusivityK = diffusivityConstK; diffusivityA = diffusivityConstA; permittivity = permittivityW;
      }

      start1 += dimL;
    }

    
    Term11[q] = (concNa[q]-concConvNa[q])/deltaT;
    Term21[q] = (concK[q]-concConvK[q])/deltaT;
    Term31[q] = (concA[q]-concConvA[q])/deltaT;

    for(unsigned int k=0; k<dim; ++k){
      Term12[q][k] = 1.0*diffusivityNa*gradConcNa[q][k];
      Term22[q][k] = 1.0*diffusivityK*gradConcK[q][k];
      Term32[q][k] = 1.0*diffusivityA*gradConcA[q][k];
    }
    
    Term13[q] = diffusivityNa;
    Term23[q] = diffusivityK;
    Term33[q] = diffusivityA;
    double alphaNa, alphaK, alphaA;
    double pumpCurrent, maxPumpCurrent=0.1, k1=3.0, k2=9.0;
    if (includePump==1){
      pumpCurrent = maxPumpCurrent/(std::pow(1.0+k1/cKOut[0].val(),2)*std::pow(1.0 + k2/cNaIn[0].val(),3));
    }
    else{
      pumpCurrent = 0.0;
    }
    alphaNa = gasConstantR*temperature*1000.0/(faradayConstant);
    alphaK = alphaNa;
    alphaA = -1.0*alphaNa;
    for(unsigned int k=0; k<dim; ++k){
      Term14[q][k] = diffusivityNa*concNa[q]*gradVolt[q][k]/alphaNa;
      Term24[q][k] = diffusivityK*concK[q]*gradVolt[q][k]/alphaK;
      Term34[q][k] = diffusivityA*concA[q]*gradVolt[q][k]/alphaA;
    }

    Term41[q] = -1.0*faradayConstant*(1.0/permittivity)*(concNa[q]+concK[q]-concA[q]);
    for (unsigned int k=0; k<dim; ++k){
      Term42[q][k] = gradVolt[q][k];
    }
    for (unsigned int q_face=0; q_face<n_face_q_points; ++q_face){
      Term50[q_face] = 225.0*200.0/(faradayConstant);
      Term51[q_face] =( (gNa1*(std::pow(history_variables[q]->m_iter,3))*(std::pow(history_variables[q]->h_iter,1))+gNaLeak*includeLeak)*(voltFace[q_face].val()-E_Na)- 3.0*pumpCurrent)/faradayConstant;
      Term52[q_face] = ((gK1*(std::pow(history_variables[q]->n_iter,4))+gKLeak*includeLeak)+2.0*pumpCurrent)*(voltFace[q_face].val()-E_K)/faradayConstant;
    }
  }
}

//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values, FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double>& ULocalConv, Table<1, Sacado::Fad::DFad<double> >& R, Table<2, double>& KMatrix, deformationMap<Sacado::Fad::DFad<double>, dim>& defMap, typename DoFHandler<dim>::active_cell_iterator& cell, const unsigned int currentIncrement, std::vector<history<dim>*>& history_variables, const unsigned int currentIteration){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  unsigned int n_face_q_points=fe_face_values.n_quadrature_points;
  //temporary arrays
  Table<1, Sacado::Fad::DFad<double>> Term11 (n_q_points), Term21 (n_q_points), Term31 (n_q_points), Term41 (n_q_points);
  Table<2, Sacado::Fad::DFad<double>> Term12 (n_q_points, dim), Term14 (n_q_points, dim), Term22 (n_q_points, dim), Term24 (n_q_points, dim), Term32 (n_q_points, dim), Term34 (n_q_points, dim), Term42 (n_q_points, dim);
  Table<1, Sacado::Fad::DFad<double>> Term13 (n_q_points), Term23 (n_q_points), Term33 (n_q_points), Term50 (n_face_q_points), Term51 (n_face_q_points), Term52 (n_face_q_points);
  //evaluate stress
  evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, fe_face_values, DOF, ULocal, ULocalConv, Term11, Term12, Term13, Term14, Term21, Term22, Term23, Term24, Term31, Term32, Term33, Term34, Term41, Term42, Term50, Term51, Term52, defMap, cell, currentIncrement, history_variables, currentIteration);
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck==0){
	for (unsigned int d = 0; d < dim; d++){
	  R[i] += (Term12[q][d] + Term14[q][d] ) * fe_values.shape_grad(i,q)[d]*fe_values.JxW(q);
	}
	R[i] += Term11[q]*fe_values.shape_value(i,q)*fe_values.JxW(q);
      }
      else if (ck==1){
	R[i] += fe_values.shape_value(i,q)*Term21[q]*fe_values.JxW(q);
	for (unsigned int d = 0; d < dim; d++){
	  R[i] += (Term22[q][d] + Term24[q][d])*fe_values.shape_grad(i,q)[d]*fe_values.JxW(q);
	}
      }
      else if (ck==2){
        R[i] += fe_values.shape_value(i,q)*Term31[q]*fe_values.JxW(q);
        for (unsigned int d = 0; d < dim; d++){
          R[i] += (Term32[q][d] + Term34[q][d])*fe_values.shape_grad(i,q)[d]*fe_values.JxW(q);
        }
      }
      else if (ck==3){
        R[i] += fe_values.shape_value(i,q)*Term41[q]*fe_values.JxW(q);
        for (unsigned int d = 0; d < dim; d++){
          R[i] += Term42[q][d]*fe_values.shape_grad(i,q)[d]*fe_values.JxW(q);
        }
      }
      
    }
  }

  if (currentIncrement>15){
    unsigned int faceID;
    //Upper Membrane flux
    int check2 = 0;
    double start = 5.0;
    while (start<neuronLength){
    
      if ( (cell->face(2)->center()[0]> start) && (cell->face(2)->center()[0]<(start+dimB)) ){
	check2=1;
      }
      start += dimL;
    }
    
    if ((cell->face(3)->center()[1]>1.6) && (1.6> cell->face(2)->center()[1]) && check2==1 ){
      faceID=2;
      fe_face_values.reinit (cell, faceID);
      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) { //fe face values ? no as rest equal 0                                                                                      
	//evaluate U                                                                                                                                                                      
	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	//const unsigned int dofs = fe_values.dofs_per_cell;                                                                                              
	if (ck==0){                                                                                                                                
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    if (currentIncrement>15 && currentIncrement<125){
	      if ((cell->face(3)->center()[1]>1.6) && (1.6> cell->face(2)->center()[1]) && (cell->face(2)->center()[0]>5.0) && (cell->face(2)->center()[0]<(5.0+dimB)) ){
		R[i] += -1.0*Term50[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral
	      }
	    }
	    R[i] += 1.0*Term51[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral
	  }
	}
	if (ck==1){
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    R[i] += 1.0*Term52[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral                                                                                                                              
	  }
	}
      }

      faceID=3;
      fe_face_values.reinit (cell, faceID);
      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) { //fe face values ? no as rest equal 0                                                                                                 
	//evaluate U
	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	if (ck==0){                                                                                                                                                                                                                      
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    R[i] += -1.0*Term51[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral, more -ve charge outside due to Na influx inward.           
	  }
	}
	if (ck==1){
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    R[i] += -1.0*Term52[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral, K conc increases outside.
	  }
	}
      }
    }
    
    
    //Lower Membrane flux
    if ((cell->face(3)->center()[1]>=(0.5) ) && ( (0.5)>= cell->face(2)->center()[1])  && check2 ==1   ){
      faceID=2;
      fe_face_values.reinit (cell, faceID);
      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) { //fe face values ? no as rest equal 0                                                                                                                                       
                                                                                                                                                                
	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	                                                                                                                                                                 
	if (ck==0){
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    R[i] += -1.0*Term51[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral                                                                                                                              
	  }
	}
	if (ck==1){
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    R[i] += -1.0*Term52[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral
	  }
	}
      }
      
      faceID=3;
      fe_face_values.reinit (cell, faceID);
      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) {
                                                                                                                                                        
	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	if (ck==0){
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    R[i] += 1.0*Term51[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral, makes +ve intra, Na flows in.
	  }
	}
	if (ck==1){
	  for (unsigned int q=0; q<n_face_q_points; ++q){
	    R[i] += 1.0*Term52[q]*fe_face_values.shape_value(i, q)*fe_face_values.JxW(q); //surface integral, makes intra -ve, K flows out.
	  }
	}
      }
    }
    
  }
  

}

#endif /* MECHANICS_H_ */
