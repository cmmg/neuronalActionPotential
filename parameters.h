//problem geometry, mesh control
#define DIMS 4 //Dimension of the parameter to be solved(c_Na, c_K, c_A, voltage)
#define dimension 2 //2D problem
#define problemWidth 1.0 //
#define refinementFactor 7
#define dimL 70.0
#define dimB 5.0
#define dimH 0.01792
#define neuronLength 1500.0

//Simple diffusion
#define deltaT 0.01
#define diffusivityConstNa 1.33
#define diffusivityConstK 1.96
#define diffusivityConstA 2.00
#define memRelativePerm 13.39*0.25
#define permittivityM memRelativePerm*8.88541/1000.0
#define permittivityW 80.0*8.88541/1000.0
#define faradayConstant 96485.0
#define gasConstantR 8.31454
#define temperature 279.450
#define gNa 120.0*10.0
#define gK 36.0*10.0;
#define gNaLeak 0.065*10.0
#define gKLeak 0.435*10.0
#define gA 0.0
#define V_rest -70.0
#define E_Na (115.0 + V_rest)
#define E_K (-12.0 + V_rest)
#define initialIntraConcNa 10.0
#define initialIntraConcK 140.0
#define initialIntraConcA 150.0128
#define initialExtraConcNa 155.0
#define initialExtraConcK 3.5
#define initialExtraConcA 158.5

#define includePump 0
#define includeLeak 1.0

//time step controls
#define DT deltaT //ms
#define Total_time 15 //ms


//output controls
#define outputFileName "solution"
