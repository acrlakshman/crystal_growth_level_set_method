//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

const int m=51, n=51, l=51;

RANGE<TV> domain(-1.5,1.5,-1.5,1.5,-1.5,1.5);
//RANGE<TV> domain(-2.,2.,-2.,2.,-2.,2.);
//RANGE<TV> domain(0.,1.,.0,1.,.0,1.);
T_GRID grid(m,n,l,domain);
T_ARRAYS phi(grid.Domain_Indices());
T_ARRAYS temperature(grid.Domain_Indices());
TV_ARRAYS grad_temperature(grid.Domain_Indices());
T_ARRAYS T_interface(temperature);
TV_ARRAYS root(grid.Domain_Indices());
TV_ARRAYS psi(grid.Domain_Indices());
T_ARRAYS Curvature(grid.Domain_Indices());
TV_ARRAYS velocity(grid.Domain_Indices());
T_ARRAYS Iphi(phi); // Interpolated phi
TV_ARRAYS Ipsi(psi); // Interpolated psi
T_ARRAYS Rphi(phi); // Reinitialized phi
TV_ARRAYS GradUDotPsi(grid.Domain_Indices());
T_ARRAYS temperature_P(grid.Domain_Indices());
T_ARRAYS temperature_N(grid.Domain_Indices());
TV_ARRAYS grad_temperature_P(grid.Domain_Indices());
TV_ARRAYS grad_temperature_N(grid.Domain_Indices());
LEVELSET_3D<T_GRID> levelset(grid,phi,0);
DUALCONTOUR_3D<T> dualcontour(levelset); // object to construct triangulated mesh using dualcontour algorithm
TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create(); // variable to store triangulated surface
T dt;
T time=(T).0;

int order=2; // Gradient order numerically
int Reinit_Steps=3; // Reinitialization steps // Currently this uses fast marching method upto distance 5*dx
int Extrapolate_Steps=5; // Extrapolation steps
T_ARRAYS velocity_x(grid.Domain_Indices());
T_ARRAYS velocity_y(grid.Domain_Indices());
T_ARRAYS velocity_z(grid.Domain_Indices());
TV_ARRAYS velocity_x_dxdydz(grid.Domain_Indices());
TV_ARRAYS velocity_y_dxdydz(grid.Domain_Indices());
TV_ARRAYS velocity_z_dxdydz(grid.Domain_Indices());
T_ARRAYS velocity_n(grid.Domain_Indices());

const T diffusion_coefficient=(T)1.;
const T epsilon_c=(T).002;
const T epsilon_v=(T).002;
T_ARRAYS temperature_new(grid.Domain_Indices());

/**********************************************************/
// Initialization
/**********************************************************/
for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
  const T_INDEX& index=iterator.Index();

  const TV X=grid.X(index);
  root(index)=X;

  // Levelset field
  // ## Spheres
  //phi(index)=(X-TV::Axis_Vector(1)*.0-TV::Axis_Vector(2)*.0-TV::Axis_Vector(3)*.0).Magnitude()-(T).2;
  //phi(index)=(X-TV::Axis_Vector(1)*.35-TV::Axis_Vector(2)*.35-TV::Axis_Vector(3)*.35).Magnitude()-(T).15;
  //phi(index)=(X-TV::Axis_Vector(1)*2.-TV::Axis_Vector(2)*2.75-TV::Axis_Vector(3)*2.).Magnitude()-(T).5;

  // ## Ellipsoid
  //phi(index)=sqrt(sqr(X.x-(T).5)/sqr(.2)+sqr(X.y-(T).5)/sqr(.1)+sqr(X.z-(T).5)/sqr(.15))-(T)1.;

  // ## Sphere with varying gradient [Anumolu and Trujillo, APS DFD 2012]
  //phi(index)=(sqr(X.x-(T)1.)+sqr(X.y-(T)1.)+sqr(X.z-(T)1.)+(T).1)*(sqrt(sqr(X.x)+sqr(X.y)+sqr(X.z))-1);

  // 6 spheres for Stephan's problem
  // In literature each sphere is of radius .1
  T phi_1_6[]={.0,.0,.0,.0,.0,.0};
  phi_1_6[0]=(X-TV::Axis_Vector(1)*.05).Magnitude()-(T).2;
  phi_1_6[1]=(X+TV::Axis_Vector(1)*.05).Magnitude()-(T).2;
  phi_1_6[2]=(X-TV::Axis_Vector(2)*.05).Magnitude()-(T).2;
  phi_1_6[3]=(X+TV::Axis_Vector(2)*.05).Magnitude()-(T).2;
  phi_1_6[4]=(X-TV::Axis_Vector(3)*.05).Magnitude()-(T).2;
  phi_1_6[5]=(X+TV::Axis_Vector(3)*.05).Magnitude()-(T).2;

  // ## min(6 spheres from above)
  //phi(index)=min(phi_1_6[0],min(phi_1_6[1],min(phi_1_6[2],min(phi_1_6[3],min(phi_1_6[4],phi_1_6[5])))));
  // ## min(2 spheres from above)
  phi(index)=min(phi_1_6[0],phi_1_6[1]);

  // Gradient of levelset field
  // Analytical equation for sphere (NOT USED, since gradient is computed numerically later)
  //  psi(index)=(X-TV::Axis_Vector(1)*.5-TV::Axis_Vector(2)*.5-TV::Axis_Vector(3)*.5)/(phi(index)+(T).25);

  // Temperature field
  //temperature(index)=(T)500;
  //temperature(index)=1000*sin(atan(1)*4*X(1));
  
  if(phi(index)<=0)
    temperature(index)=(T)0.;
  else
    temperature(index)=(T)-.5;
  
  temperature_new(index)=temperature(index);
  temperature_P(index)=temperature(index);
  temperature_N(index)=temperature(index);
  T_interface(index)=temperature(index);

  // Velocity field
  velocity(index)=TV();
  velocity_x(index)=velocity(index).x;
  velocity_y(index)=velocity(index).y;
  velocity_z(index)=velocity(index).z;

  velocity_x_dxdydz(index)=TV();
  velocity_y_dxdydz(index)=TV();
  velocity_z_dxdydz(index)=TV();

  velocity_n(index)=T();

  // Initializing GradUDotPsi
  GradUDotPsi(index)=TV();
}

Rphi=phi;
// Gradient of levelset field (Numerically)
cout<<"Numerically computing gradient"<<endl;
Compute_Gradient(grid,phi,psi,order);
Compute_Curvature(grid,psi,Curvature);

// Gradient of temperature fields (Numerically)
Compute_Gradient(grid,temperature,grad_temperature,order);
Compute_Gradient(grid,temperature_P,grad_temperature_P,order);
Compute_Gradient(grid,temperature_N,grad_temperature_N,order);

// Creating object for HermiteInterpolation
HermiteInterpolation<T> interpolant(grid,root,phi,psi);
Reinitialize<T> signed_distance(grid,phi,psi,Reinit_Steps);

// Reinitialization at initial time
//signed_distance.Reinitialize_Phi();
//Rphi=signed_distance.Rphi_F();
levelset.Get_Signed_Distance_Using_FMM(Rphi,.0,5*grid.dX.x);

surface=dualcontour.Create_Triangulated_Surface_From_Levelset(levelset);
//

Extrapolate<T> extrapolate_temperature(grid,temperature,phi,psi,Extrapolate_Steps);
// objects for temperature_P and temperature_N: useful to compute its gradient at interfacial locations

Reinitialize<T> Get_Interface_Location_temperature_P(grid,phi,psi,Reinit_Steps);
Reinitialize<T> Get_Interface_Location_temperature_N(grid,phi,psi,Reinit_Steps);
TV_ARRAYS& root_temperature_P=Get_Interface_Location_temperature_P.root;
HermiteInterpolation<T> interpolant_temperature_P_phi(grid,root_temperature_P,phi,psi);
TV_ARRAYS& root_temperature_N=Get_Interface_Location_temperature_N.root;
HermiteInterpolation<T> interpolant_temperature_N_phi(grid,root_temperature_N,phi,psi);

// Interpolant to interpolate temperature_P and temperature_N
HermiteInterpolation<T> interpolant_temperature_P(grid,root_temperature_P,temperature_P,grad_temperature_P);
HermiteInterpolation<T> interpolant_temperature_N(grid,root_temperature_N,temperature_N,grad_temperature_N);
cout<<"Initialized all fields"<<endl;
