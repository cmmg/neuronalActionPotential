//problem geometry, mesh control
#define DIMS 2                              //Dimension of the parameter to be solved(voltage)
#define dimension 1                         //1D problem
#define problemWidth 1.0                    //length of neuronal axon in cm
#define refinementFactor 8                  //mesh refinement factor


//Hodgkin-Huxley Parameters
#define C_m 1.45                            //membrane capacitance
#define MC_m 3.2e-6                         //myelin capacitance  
#define g_Na 120.0                          //peak conductivity of Na ion channels in mS/cm^2  
#define g_K 36.0                            //peak conductivity of K ion channels in mS/cm^2
#define g_L 0.3                             //conductivity of leak channels in mS/cm^2
#define V_rest -70.0                        //membrane resting potential in mV 
#define E_Na (115.0 + V_rest)               //Nernst potential of Na ions
#define E_K (-12.0 + V_rest)                //Nernst potential of K ions
#define E_L (10.6 + V_rest)                 //potential of leak channels
#define a_RADIUS 2.0*0.000647668            //normalized axonal radius keeping R_SP_RESIST = 1.0
#define R_SP_RESIST 1.0                     //Resistivity
#define VBoundary 200.0                     //initial dirichlet boundary condition on the left boundary in mV            
#define Myelin 0                            //Myelin present or absent


//time step controls
#define DT 0.01                             //timestep in milli-seconds
#define Total_time 20.0                     //total time in milli-seconds

//output controls
#define outputFileName "solution"
