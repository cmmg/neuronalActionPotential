//problem geometry, mesh control
#define DIMS 2                              //Dimension of the parameter to be solved(voltage)
#define dimension 1                         //1D problem
#define problemWidth 0.1                    //length of neuronal axon in cm
#define refinementFactor 9                  //mesh refinement factor

//Hodgkin-Huxley Parameters
#define C_m 1.23                            //membrane capacitance
#define C_my 0.113                          //myelin capacitance
#define R_m 24.8                            //leak resistance
#define R_my 63.7                           //myelin resistance
#define R_i 0.712                           //axonal resistance
#define g_L (1/R_m)                         //conductivity of leak channels in mS/cm^2 
#define g_Na 120.0                          //peak conductivity of Na ion channels in mS/cm^2
#define g_K 36.0                            //peak conductivity of K ion channels in mS/cm^2
#define V_rest -65.0                        //membrane resting potential in mV
#define E_Na (115.0 + V_rest)               //Nernst potential of Na ions
#define E_K (-12.0 + V_rest)                //Nernst potential of K ions
#define E_L 0.0                             //potential of leak channels
#define RADIUS 0.55/1000.0                  //axonal radius
#define periRadius (0.55 + 0.0123)/1000.0   //peri-axonal radius
#define VBoundary 100.0                     //initial dirichlet boundary condition on the left boundary in mV
#define Myelin 1                            //Myelin present or absent
#define rPA 96.3*100.0                      //peri-axonal resistance
#define rPN 321.0*100.0                     //paranodal resistance


//time step controls
#define DT 0.01 //ms                        //timestep in milli-seconds
#define Total_time 10000 //ms               //total time in milli-seconds

//output controls
#define outputFileName "solution"
