//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "Extrapolate.h"

using namespace std;
using namespace PhysBAM;

/*****************************************************/
// Constructor
/*****************************************************/
template<class T>
Extrapolate<T>::Extrapolate(const T_GRID& grid_, T_ARRAYS& u_, T_ARRAYS& phi_, TV_ARRAYS& psi_, int Extrapolate_Steps_)
  :grid(grid_),u(u_),phi(phi_),psi(psi_),Extrapolate_Steps(Extrapolate_Steps_)
{
  x_cells=grid.Numbers_Of_Nodes().x;
  y_cells=grid.Numbers_Of_Nodes().y;
  z_cells=grid.Numbers_Of_Nodes().z;
  u_temp=u,u_n=u,grad_u_n=psi;
  Iu_n_P=u,Eu_n_P=u,Egrad_u_n_P=psi;
  Iu_P=u,Eu_P=u,Egrad_u_P=psi,grad_u=psi;
  Iu_n_N=u,Eu_n_N=u,Egrad_u_n_N=psi;
  Iu_N=u,Eu_N=u,Egrad_u_N=psi;
  root=psi;
  velocity=psi;

  dt=0.5*min(grid.dX.x,grid.dX.y,grid.dX.z);
  cout<<"Extrapolation time step:="<<dt<<endl;

  Normal=psi;
  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    Normal(index)=psi(index)/psi(index).Magnitude();
    root(index)=grid.X(index);
  }
}

/*****************************************************/
// Destructor
/*****************************************************/
template<class T>
Extrapolate<T>::~Extrapolate()
{}

/*****************************************************/
// Function: Extrapolate_Phi
/*****************************************************/
template<class T>
void Extrapolate<T>::Extrapolate_Phi()
{
  /*** ### Constant extrapolation ### ***/
  // START COMMENT
  // Computing gradient of 'u' using central difference scheme
  u_temp=u;
  Compute_Gradient(grad_u,2);
  Eu_P=u,Egrad_u_P=grad_u; // Extrapolating to positive \phi
  Eu_N=u,Egrad_u_N=grad_u; // Extrapolating to negative \phi

  // ### Solving for u_P ### //
  cout<<"Solving for u_P"<<endl;
  // Creating object for HermiteInterpolation
  HermiteInterpolation<T> interpolant_u_P(grid,root,Eu_P,Egrad_u_P);

  int Extrapolate_Counter=0;
  while(Extrapolate_Counter<Extrapolate_Steps){
    // Find root of the characteristic
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      velocity(index)=(T)Heaviside(phi(index))*Normal(index);
      root(index)=X-velocity(index)*dt;
    }

    interpolant_u_P.Interpolate_Phi();
    Iu_P=interpolant_u_P.Iphi_F();

    // Solving for u using du/dt=0
    Eu_P=Iu_P;

    u_temp=Eu_P; // This is used to compute any griadient
    Compute_Gradient(Egrad_u_P,2);
    Extrapolate_Counter++;
  } // End u_P extrapolation

  // ### Solving for u_N ### //
  cout<<"Solving for u_N"<<endl;
  // Creating object for HermiteInterpolation
  HermiteInterpolation<T> interpolant_u_N(grid,root,Eu_N,Egrad_u_N);

  Extrapolate_Counter=0;
  while(Extrapolate_Counter<Extrapolate_Steps){
    // Find root of the characteristic
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      velocity(index)=(T)Heaviside(-phi(index))*Normal(index);
      root(index)=X-velocity(index)*dt;
    }

    interpolant_u_N.Interpolate_Phi();
    Iu_N=interpolant_u_N.Iphi_F();

    // Solving for u using du/dt=0
    Eu_N=Iu_N;

    u_temp=Eu_N;
    Compute_Gradient(Egrad_u_N,2);
    Extrapolate_Counter++;
  } // End u_N extrapolation
  // END COMMENT
  /*** ### END Constant extrapolation ### ***/

  /*** ### Linear extrapolation ### ***/
  /* // START COMMENT
  // Computing gradient of 'u' using central difference scheme
  u_temp=u;
  Compute_Gradient(grad_u,2);
  Eu_P=u,Egrad_u_P=grad_u;
  Eu_N=u,Egrad_u_N=grad_u;

  // Computing u_n and its gradient using central difference scheme
  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    u_n(index)=(Normal(index).x*grad_u(index).x)+(Normal(index).y*grad_u(index).y);
  }
  u_temp=u_n;
  Compute_Gradient(grad_u_n,2);
  Eu_n_P=u_n,Egrad_u_n_P=grad_u_n;

  // ### Solving for u_n_P ### //
  cout<<"Solving for u_n_P"<<endl;
  // Creating object for HermiteInterpolation
  HermiteInterpolation<T> interpolant_u_n_P(grid,root,Eu_n_P,Egrad_u_n_P);

  int Extrapolate_Counter=0;
  while(Extrapolate_Counter<2*Extrapolate_Steps){
    // Find root of the characteristic
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      velocity(index)=(T)Heaviside(phi(index)+grid.dX.x)*Normal(index);
      root(index)=X-velocity(index)*dt;
    }

    interpolant_u_n_P.Interpolate_Phi();
    Iu_n_P=interpolant_u_n_P.Iphi_F();

    // Solving for u_n using du_n/dt=0
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      Eu_n_P(index)=Iu_n_P(index);
    }
    u_temp=Eu_n_P;
    Compute_Gradient(Egrad_u_n_P,2);
    Extrapolate_Counter++;
  } // End u_n_P extrapolation

  // ### Solving for u_P ### //
  cout<<"Solving for u_P"<<endl;
  // Creating object for HermiteInterpolation
  HermiteInterpolation<T> interpolant_u_P(grid,root,Eu_P,Egrad_u_P);

  Extrapolate_Counter=0;
  while(Extrapolate_Counter<Extrapolate_Steps){
    // Find root of the characteristic
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      velocity(index)=(T)Heaviside(phi(index))*Normal(index);
      root(index)=X-velocity(index)*dt;
    }

    interpolant_u_P.Interpolate_Phi();
    Iu_P=interpolant_u_P.Iphi_F();

    // Solving for u using du/dt=H(\phi)*u_n
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      Eu_P(index)=Iu_P(index)+(T)Heaviside(phi(index))*Eu_n_P(index)*dt;
    }
    u_temp=Eu_P;
    Compute_Gradient(Egrad_u_P,2);
    Extrapolate_Counter++;
  } // End u_P extrapolation

  // ### Computing u_n_N ### //
  cout<<"Solving for u_n_N"<<endl;
  Eu_n_N=u_n,Egrad_u_n_N=grad_u_n;

  // ### Solving for u_n_N ### //
  // Creating object for HermiteInterpolation
  HermiteInterpolation<T> interpolant_u_n_N(grid,root,Eu_n_N,Egrad_u_n_N);

  Extrapolate_Counter=0;
  while(Extrapolate_Counter<2*Extrapolate_Steps){
    // Find root of the characteristic
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      velocity(index)=(T)Heaviside(-phi(index)+grid.dX.x)*Normal(index);
      root(index)=X-velocity(index)*dt;
    }

    interpolant_u_n_N.Interpolate_Phi();
    Iu_n_N=interpolant_u_n_N.Iphi_F();

    // Solving for u_n using du_n/dt=0
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      Eu_n_N(index)=Iu_n_N(index);
    }
    u_temp=Eu_n_N;
    Compute_Gradient(Egrad_u_n_N,2);
    Extrapolate_Counter++;
  } // End u_n_N extrapolation

  // ### Solving for u_N ### //
  cout<<"Solving for u_N"<<endl;
  // Creating object for HermiteInterpolation
  HermiteInterpolation<T> interpolant_u_N(grid,root,Eu_N,Egrad_u_N);

  Extrapolate_Counter=0;
  while(Extrapolate_Counter<Extrapolate_Steps){
    // Find root of the characteristic
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      velocity(index)=(T)Heaviside(-phi(index))*Normal(index);
      root(index)=X-velocity(index)*dt;
    }

    interpolant_u_N.Interpolate_Phi();
    Iu_N=interpolant_u_N.Iphi_F();

    // Solving for u using du/dt=H(\phi)*u_n
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      Eu_N(index)=Iu_N(index)+(T)Heaviside(-phi(index))*Eu_n_N(index)*dt;
    }
    u_temp=Eu_N;
    Compute_Gradient(Egrad_u_N,2);
    Extrapolate_Counter++;
  } // End u_N extrapolation
  */ // END COMMENT
}

/*****************************************************/
// Function: Gradient
/*****************************************************/
// derivative: 1-'x', 2-'y', 3-'z'
template<class T>
void Extrapolate<T>::Compute_Gradient(TV_ARRAYS &psi,int order)
{
  //cout<<"Computing Gradient"<<endl;
  typedef VECTOR<int,d> T_INDEX;
  typedef ARRAY<T,T_INDEX> T_ARRAYS;

  T_ARRAYS psi_x,psi_y,psi_z;
  psi_x=u_temp;
  psi_y=u_temp;
  psi_z=u_temp;

  // Computing to 2nd order accuracy
  Calculate_Gradient(psi_x,1,2); // psi_x
  Calculate_Gradient(psi_y,2,2); // psi_y
  Calculate_Gradient(psi_z,3,2); // psi_z
  if(order==3){
    Calculate_Gradient(psi_x,1,order); // psi_x
    Calculate_Gradient(psi_y,2,order); // psi_y
    Calculate_Gradient(psi_z,3,order); // psi_z
  }

 for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
   const T_INDEX& index=iterator.Index();
   psi(index)(1)=psi_x(index);
   psi(index)(2)=psi_y(index);
   psi(index)(3)=psi_z(index);
 }
}

/************* Calculating gradient ***************/
// derivative: 1-'x', 2-'y', 3-'z'
template<class T>
void Extrapolate<T>::Calculate_Gradient(T_ARRAYS &psi_p,int component,int order)
{
  typedef VECTOR<int,d> T_INDEX;

  int x_cells=grid.Numbers_Of_Nodes().x;
  int y_cells=grid.Numbers_Of_Nodes().y;
  int z_cells=grid.Numbers_Of_Nodes().z;

  switch (component)
    {
    case 1:
      if(order==2){
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  psi_p(index)=T();
	  if(index(1)==1){
	    psi_p(index)+=(u_temp(index+T_INDEX::Axis_Vector(1))-u_temp(index))*(grid.One_Over_DX()(1));
	  }
	  else if(index(1)==x_cells){
	    psi_p(index)+=(u_temp(index)-u_temp(index-T_INDEX::Axis_Vector(1)))*(grid.One_Over_DX()(1));
	  }
	  else{
	    psi_p(index)+=(u_temp(index+T_INDEX::Axis_Vector(1))-u_temp(index-T_INDEX::Axis_Vector(1)))*(T)0.5*(grid.One_Over_DX()(1));
	  }
	}
      }
      else{
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  if(index(1)==2){
	    psi_p(index)=T();
	    psi_p(index)+=((-u_temp(index+T_INDEX::Axis_Vector(1)*2))+((T)6.*u_temp(index+T_INDEX::Axis_Vector(1)))-((T)3.*u_temp(index))-((T)2.*u_temp(index-T_INDEX::Axis_Vector(1))))*((T)(1./6.))*(grid.One_Over_DX()(1));
	  }
	  else if(index(1)>2 && index(1)<=(x_cells-1)){
	    psi_p(index)=T();
	    psi_p(index)+=(((T)2.*u_temp(index+T_INDEX::Axis_Vector(1)))+((T)3.*u_temp(index))-((T)6.*u_temp(index-T_INDEX::Axis_Vector(1)))+(u_temp(index-T_INDEX::Axis_Vector(1)*2)))*((T)(1./6.))*(grid.One_Over_DX()(1));
	  }
	  else{}
	}
      }
      break;

    case 2:
      if(order==2){
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  psi_p(index)=T();
	  if(index(2)==1){
	    psi_p(index)+=(u_temp(index+T_INDEX::Axis_Vector(2))-u_temp(index))*(grid.One_Over_DX()(2));
	  }
	  else if(index(2)==y_cells){
	    psi_p(index)+=(u_temp(index)-u_temp(index-T_INDEX::Axis_Vector(2)))*(grid.One_Over_DX()(2));
	  }
	  else{
	    psi_p(index)+=(u_temp(index+T_INDEX::Axis_Vector(2))-u_temp(index-T_INDEX::Axis_Vector(2)))*(T)0.5*(grid.One_Over_DX()(2));
	  }
	}
      }
      else{
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  if(index(2)==2){
	    psi_p(index)=T();
	    psi_p(index)+=((-u_temp(index+T_INDEX::Axis_Vector(2)*2))+((T)6.*u_temp(index+T_INDEX::Axis_Vector(2)))-((T)3.*u_temp(index))-((T)2.*u_temp(index-T_INDEX::Axis_Vector(2))))*((T)(1./6.))*(grid.One_Over_DX()(2));
	  }
	  else if(index(2)>2 && index(2)<=(y_cells-1)){
	    psi_p(index)=T();
	    psi_p(index)+=(((T)2.*u_temp(index+T_INDEX::Axis_Vector(2)))+((T)3.*u_temp(index))-((T)6.*u_temp(index-T_INDEX::Axis_Vector(2)))+(u_temp(index-T_INDEX::Axis_Vector(2)*2)))*((T)(1./6.))*(grid.One_Over_DX()(2));
	  }
	  else{}
	}
      }
      break;

    case 3:
      if(order==2){
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  psi_p(index)=T();
	  if(index(3)==1){
	    psi_p(index)+=(u_temp(index+T_INDEX::Axis_Vector(3))-u_temp(index))*(grid.One_Over_DX()(3));
	  }
	  else if(index(3)==z_cells){
	    psi_p(index)+=(u(index)-u(index-T_INDEX::Axis_Vector(3)))*(grid.One_Over_DX()(3));
	  }
	  else{
	    psi_p(index)+=(u_temp(index+T_INDEX::Axis_Vector(3))-u_temp(index-T_INDEX::Axis_Vector(3)))*(T)0.5*(grid.One_Over_DX()(3));
	  }
	}
      }
      else{
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  if(index(3)==2){
	    psi_p(index)=T();
	    psi_p(index)+=((-u_temp(index+T_INDEX::Axis_Vector(3)*2))+((T)6.*u_temp(index+T_INDEX::Axis_Vector(3)))-((T)3.*u_temp(index))-((T)2.*u_temp(index-T_INDEX::Axis_Vector(3))))*((T)(1./6.))*(grid.One_Over_DX()(3));
	  }
	  else if(index(3)>2 && index(3)<=(z_cells-1)){
	    psi_p(index)=T();
	    psi_p(index)+=(((T)2.*u_temp(index+T_INDEX::Axis_Vector(3)))+((T)3.*u_temp(index))-((T)6.*u_temp(index-T_INDEX::Axis_Vector(3)))+(u_temp(index-T_INDEX::Axis_Vector(3)*2)))*((T)(1./6.))*(grid.One_Over_DX()(3));
	  }
	  else{}
	}
      }
      break;
    }
}

template class Extrapolate<float>;
