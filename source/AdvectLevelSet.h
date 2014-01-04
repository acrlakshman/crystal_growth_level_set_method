//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

/**********************************************************/
// Using semi-Lagrangian approach
/**********************************************************/

// Reinitializing level set and its gradient fields
Rphi=phi;
//Timing
//clock_t start, end;
//double cpuTime;

//start = clock();
//signed_distance.Reinitialize_Phi(); // PDE approach
//end=clock();
//cpuTime=(end-start)*1000.0/(CLOCKS_PER_SEC);
//cout<<"cputime = "<<cpuTime<<endl;

//Rphi=signed_distance.Rphi_F(); // Get signed distance obtained using PDE approach
//psi=signed_distance.Rpsi_F(); // Get gradient of signed distance function
levelset.Get_Signed_Distance_Using_FMM(Rphi,.0,5*grid.dX.x); // Using FMM
phi=Rphi;
Compute_Gradient(grid,phi,psi,order);

// Obtain root of the characteristic
for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
  const T_INDEX& index=iterator.Index();
  const TV X=grid.X(index);

  T TPeriod=2.5; // For Enright test case
  //velocity(index)=TV();//TV::All_Ones_Vector();

  // horizontal field
  //velocity(index).x=(T)1.;

  // circular field
  //velocity(index).x=-atan(1)*(T)4.*(X(2)-(T).5)/(T)3.14;
  //velocity(index).y=atan(1)*(T)4.*(X(1)-(T).5)/(T)3.14;

  // single vortex field
  /*
  velocity(index).x=(T)2.*sqr(sin(atan(1)*(T)4.*X(1)))*sin((T)2.*atan(1)*(T)4.*X(2))*sin((T)2.*atan(1)*(T)4.*X(3))*cos(atan(1)*(T)4.*time/TPeriod);
  velocity(index).y=-sqr(sin(atan(1)*(T)4.*X(2)))*sin((T)2.*atan(1)*(T)4.*X(1))*sin((T)2.*atan(1)*(T)4.*X(3))*cos(atan(1)*(T)4.*time/TPeriod);
  velocity(index).z=-sqr(sin(atan(1)*(T)4.*X(3)))*sin((T)2.*atan(1)*(T)4.*X(2))*sin((T)2.*atan(1)*(T)4.*X(1))*cos(atan(1)*(T)4.*time/TPeriod);
  */
  velocity_x(index)=velocity(index).x;
  velocity_y(index)=velocity(index).y;
  velocity_z(index)=velocity(index).z;

  root(index)=X-velocity(index)*dt;
 }
Compute_Gradient(grid,velocity_x,velocity_x_dxdydz,order);
Compute_Gradient(grid,velocity_y,velocity_y_dxdydz,order);
Compute_Gradient(grid,velocity_z,velocity_z_dxdydz,order);
Compute_InnerProduct(grid,'T','V',velocity_x_dxdydz,velocity_y_dxdydz,velocity_z_dxdydz,psi,GradUDotPsi);

// Interpolate phi and psi to root
interpolant.Interpolate_Phi();

Iphi=interpolant.Iphi_F();
Ipsi=interpolant.Ipsi_F();

// Update levelset and its gradient
for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
  const T_INDEX& index=iterator.Index();

  // d(phi)/dt=0, hence update is not needed
  phi(index)=Iphi(index);

  // d(psi)/dt=-\grad(velocity) \cdot \psi
  psi(index)=Ipsi(index)-GradUDotPsi(index);
 }
//Compute_Gradient(grid,phi,psi,order); // Use this, if not willing to use gradient obtained using gradient augmented level set method
Compute_Curvature(grid,psi,Curvature);

