//problem geometry, mesh control
#define DIMS 2 //Dimension of the parameter to be solved(voltage)
#define dimension 1 //3D problem
#define problemWidth 0.1 //
#define refinementFactor 9

//Hodgkin-Huxley Parameters
#define C_m 1.23  // uF/cm^2
#define C_my 0.113
#define R_m 24.8    // ohm cm
#define R_my 63.7
#define R_i 0.712
#define g_L (1/R_m) //1/ohm cm
#define g_Na 120.0
#define g_K 36.0
#define V_rest -65.0
#define E_Na (115.0 + V_rest)
#define E_K (-12.0 + V_rest)
#define E_L 0.0
#define RADIUS 0.55/1000.0
#define periRadius (0.55 + 0.0123)/1000.0   // um
#define VBoundary 100.0 
#define Myelin 1
#define rPA 96.3*100.0
#define rPN 321.0*100.0


//time step controls
#define DT 0.01 //ms
#define Total_time 10000 //ms

//output controls
#define outputFileName "solution"
