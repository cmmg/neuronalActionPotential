//problem geometry, mesh control
#define DIMS 2 //Dimension of the parameter to be solved(voltage)
#define dimension 1 //3D problem
#define problemWidth 1.0  // in cm
#define refinementFactor 8


//Hodgkin-Huxley Parameters
#define C_m 1.45  //uF/cm^2
#define MC_m 3.2e-6
#define g_Na 120.0 //mS/cm^2 
#define g_K 36.0  
#define g_L 0.3 
#define V_rest -70.0  //mV 
#define E_Na (115.0 + V_rest)
#define E_K (-12.0 + V_rest)
#define E_L (10.6 + V_rest)
#define a_RADIUS 2.0*0.000647668
#define R_SP_RESIST 1.0 
#define VBoundary 200.0 //mV
#define Myelin 0


//time step controls
#define DT 0.01 // unit: ms
#define Total_time 20.0 //ms

//output controls
#define outputFileName "solution"
