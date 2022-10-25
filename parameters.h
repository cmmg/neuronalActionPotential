//problem geometry, mesh control
#define DIMS 4                                                //Dimension of the parameters to be solved(c_Na, c_K, c_A, voltage)
#define dimension 2                                           //2D problem
#define dimL 70.0                                             //Internodal distance in micro-meter
#define dimB 5.0                                              //span of node of Ranvier in micro-meter
#define dimH 0.03583                                          //used to define the thickness of myelin sheath in the internodal region  
#define neuronLength 1500.0                                   //total length of the neuron in micro-meter

//Simple diffusion
#define deltaT 0.01                                           //timestep in milli-seconds
#define diffusivityConstNa 1.33                               //Diffusion coefficient of Na ions
#define diffusivityConstK 1.96                                //Diffusion coefficient of K ions
#define diffusivityConstA 2.00                                //Diffusion coefficient of Cl ions
#define memRelativePerm 13.39*0.25                            //relative permittivity of membrane
#define permittivityM memRelativePerm*8.88541/1000.0          //permittivity of membrane. 
#define permittivityW 80.0*8.88541/1000.0                     //permittivity of the intra/extra- cellular regions
#define faradayConstant 96485.0                               //Faraday constant
#define gasConstantR 8.31454                                  //gas constant R
#define temperature 279.450                                   //Temperature in Kelvin
#define gNa 120.0*10.0                                        //peak conductivity of Na ion channels
#define gK 36.0*10.0;                                         //peak conductivity of K ion channels
#define gNaLeak 0.065*10.0                                    //leak conductance of Na ions
#define gKLeak 0.435*10.0                                     //leak conductance of K ions
#define gA 0.0                                                //leak conductance of Cl ions  
#define V_rest -70.0                                          //membrane resting potential
#define E_Na (115.0 + V_rest)                                 //Nernst potential of Na ions
#define E_K (-12.0 + V_rest)                                  //Nernst potential of K ions
#define initialIntraConcNa 10.0                               //initial intra-cellular sodium ion concentration (mM)
#define initialIntraConcK 140.0                               //initial intra-cellular potassium ion concentration (mM)
#define initialIntraConcA 150.0107                            //initial intra-cellular chloride ion concentration (mM)  
#define initialExtraConcNa 155.0                              //initial extra-cellular sodium ion concentration (mM)
#define initialExtraConcK 3.5                                 //initial extra-cellular potassium ion concentration (mM)
#define initialExtraConcA 158.5                               //initial extra-cellular chloride ion concentration (mM)

#define includePump 0                                         //0= no pump term, 1= include ion pump
#define includeLeak 1.0                                       //0 = no leak current, 1 = include leak current

//time step controls
#define DT deltaT //ms                                        //same as deltaT i.e. timestep in milli-seconds
#define Total_time 8 //ms                                     //total simulation time in milli-seconds

//output controls
#define outputFileName "solution"
