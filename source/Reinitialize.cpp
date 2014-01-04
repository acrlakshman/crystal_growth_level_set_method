//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "Reinitialize.h"

using namespace std;
using namespace PhysBAM;

/*****************************************************/
// Constructor
/*****************************************************/
template<class T>
Reinitialize<T>::Reinitialize(const T_GRID& grid_, T_ARRAYS& phi_, TV_ARRAYS& psi_, int Reinit_Steps_)
  :grid(grid_),phi(phi_),psi(psi_),Reinit_Steps(Reinit_Steps_)
{
  x_cells=grid.Numbers_Of_Nodes().x;
  y_cells=grid.Numbers_Of_Nodes().y;
  z_cells=grid.Numbers_Of_Nodes().z;
  Iphi=phi,Rphi=phi,Rpsi=psi;
  velocity=psi;
  root=psi;

  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    const TV X=grid.X(index);
    root(index)=X;
  }

  Interfacial_Rectangles=phi;
  Convergence_Flag=Interfacial_Rectangles;

  dt=0.9*min(grid.dX.x,grid.dX.y,grid.dX.z);
  cout<<"Reinitialization time step:="<<dt<<endl;

  N_Rectangles=0;
  criterion=0.001*grid.dX.x*grid.dX.y*grid.dX.z; // convergence criterion for modified Newtons method
  criterion=5.0e-7; // hard coded this
}

/*****************************************************/
// Destructor
/*****************************************************/
template<class T>
Reinitialize<T>::~Reinitialize()
{}

/*****************************************************/
// Function: Reinitialize_Phi
/*****************************************************/
template<class T>
void Reinitialize<T>::Reinitialize_Phi()
{
  cout<<"Reinitializing"<<endl;
  // Find root of the characteristic
  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    const TV X=grid.X(index);
    velocity(index)=sign(phi(index))*psi(index)/psi(index).Magnitude();
    root(index)=X-velocity(index)*dt;
  }

  Rphi=phi,Rpsi=psi;
  // Creating object for HermiteInterpolation
  HermiteInterpolation<T> interpolant(grid,root,Rphi,Rpsi);

  // interpolate phi, psi to root locations
  interpolant.Interpolate_Phi();
  Iphi=interpolant.Iphi_F();

  cout<<"In Newtons method"<<endl;
  Locate_Interface_Modified_Newtons_Method(interpolant);
  cout<<"Exited Newtons method"<<endl;

  // Iterating for few times now
  // Find root of the characteristic
  int Reinit_Counter=0;
  while(Reinit_Counter<Reinit_Steps){
    cout<<"Reinitializing iteration: "<<Reinit_Counter+1<<endl;
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      const TV X=grid.X(index);
      velocity(index)=sign(phi(index))*Rpsi(index)/Rpsi(index).Magnitude();
      root(index)=X-velocity(index)*dt;
    }

    // interpolate phi, psi to root locations
    interpolant.Interpolate_Phi();
    Iphi=interpolant.Iphi_F();
    TV_ARRAYS Rpsi_temp=interpolant.Ipsi_F();

    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      if(Interfacial_Rectangles(index)!=(T)1.){
	Rphi(index)=Iphi(index)+sign(phi(index))*dt;
	Rpsi(index)=Rpsi_temp(index);
      }
    }
    Reinit_Counter++;
  }
  cout<<"Reinitialized"<<endl;
}

// Locate_Interface_Modified_Newtons_Method
template<class T>
void Reinitialize<T>::Locate_Interface_Modified_Newtons_Method(HermiteInterpolation<T>& interpolant)
{
  // Locate interfacial rectangles
  Locate_Interfacial_Rectangles(phi,Iphi);

  // Computing root for tagged nodes
  // NOTE: boundary nodes excluded to facilitate search algorithm
  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices().Thickened(-1));iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    if(Interfacial_Rectangles(index)==(T)1.){
      //cout<<"Interfacial_Rectangles("<<index<<")"<<Interfacial_Rectangles(index)<<endl;
      // Modified Newtons method 3D
      Modified_Newtons_Method(interpolant,index);
      if(Convergence_Flag(index)==0.);
	//cout<<"Not converged"<<endl;
    }
  }
}

// Modified_Newtons_Method (Chopp SIAM 2001)
template<class T>
void Reinitialize<T>::Modified_Newtons_Method(HermiteInterpolation<T>& interpolant,T_INDEX index)
{
  typedef double D;
  int iter=1,MaxIter=100;
  T error=1.;
  TV X=grid.X(index),X_k=grid.X(index),X_khalf=grid.X(index),psi_k=TV();
  T phi_k=.0,psi_xk=.0,psi_yk=.0,psi_zk=.0;
  T delta_1x=.0,delta_1y=.0,delta_1z=.0,delta_2x=.0,delta_2y=.0,delta_2z=.0,Delta2DotProduct=.0;
  T_INDEX reference_index=index;

  reference_index=grid.Cell(root(index),0);

  root(index)=X_k;

  while(error>criterion && iter<MaxIter){

    phi_k=interpolant.Calculate_Interpolant(index,reference_index,0,(T).0);
    psi_xk=interpolant.Calculate_Interpolant(index,reference_index,1,grid.dX.x);
    psi_yk=interpolant.Calculate_Interpolant(index,reference_index,2,grid.dX.y);
    psi_zk=interpolant.Calculate_Interpolant(index,reference_index,3,grid.dX.z);

    psi_k(1)=psi_xk,psi_k(2)=psi_yk,psi_k(3)=psi_zk;

    // Computing delta1
    delta_1x=-phi_k*psi_xk/(sqr(psi_xk)+sqr(psi_yk)+sqr(psi_zk));
    delta_1y=-phi_k*psi_yk/(sqr(psi_xk)+sqr(psi_yk)+sqr(psi_zk));
    delta_1z=-phi_k*psi_zk/(sqr(psi_xk)+sqr(psi_yk)+sqr(psi_zk));
    // Computing X_khalf
    X_khalf(1)=X_k(1)+delta_1x;
    X_khalf(2)=X_k(2)+delta_1y;
    X_khalf(3)=X_k(3)+delta_1z;
    Delta2DotProduct=Dot_Product(X-X_k,psi_k);

    // Computing delta2
    delta_2x=(X(1)-X_k(1))-(Delta2DotProduct*psi_xk/(sqr(psi_xk)+sqr(psi_yk)+sqr(psi_zk)));
    delta_2y=(X(2)-X_k(2))-(Delta2DotProduct*psi_yk/(sqr(psi_xk)+sqr(psi_yk)+sqr(psi_zk)));
    delta_2z=(X(3)-X_k(3))-(Delta2DotProduct*psi_zk/(sqr(psi_xk)+sqr(psi_yk)+sqr(psi_zk)));
    // Updating X_k
    X_k(1)=X_khalf(1)+delta_2x;
    X_k(2)=X_khalf(2)+delta_2y;
    X_k(3)=X_khalf(3)+delta_2z;
    iter++;

    // Updating root(index) to X_k for interpolation
    root(index)=X_k;
    reference_index=grid.Cell(root(index),0);

    // Calculate error
    error=sqrt(sqr(delta_1x)+sqr(delta_1y)+sqr(delta_1z)+sqr(delta_2x)+sqr(delta_2y)+sqr(delta_2z));

  }

  Convergence_Flag(index)=0;
  // Update convergence flag
  if(fabs(error)<=criterion){
    Convergence_Flag(index)=1;
  }
  else if(fabs(error)>criterion && iter==MaxIter){
    Convergence_Flag(index)=0;
    cout<<"index="<<index<<endl;
    cout<<"Not Converged by Modified Newton's method"<<endl;
  }
  else{}

  //#########IF NOT CONVERGED ABOVE, WE FIND BY MOVING ALONG INSTANTANEOUS NORMAL
  // Losasso et al. (2006) Computers and Fluids
   // START COMMENT
  if(!Convergence_Flag(index)){
    T alpha=1;
    int quitloop=0;
    T phi_Goal=.0;
    error=1,iter=1;
    X_k=grid.X(index);
    root(index)=X_k;
    phi_k=phi(index),psi_xk=psi(index).x,psi_yk=psi(index).y,psi_zk=psi(index).z;
    psi_k(1)=psi_xk,psi_k(2)=psi_yk,psi_k(3)=psi_zk;

    while(quitloop!=1 && iter<MaxIter){
      root(index)=X_k-fabs(phi_k)*sign(phi_k)*psi_k/psi_k.Magnitude();
      reference_index=grid.Cell(root(index),0);
      if (reference_index(1)<1 || reference_index(1)>x_cells || reference_index(2)<1 || reference_index(2)>y_cells || reference_index(3)<1 || reference_index(3)>z_cells)
	break; // break the loop and continue to line search algorithm

      phi_k=interpolant.Calculate_Interpolant(index,reference_index,0,(T).0);
      psi_xk=interpolant.Calculate_Interpolant(index,reference_index,1,grid.dX.x);
      psi_yk=interpolant.Calculate_Interpolant(index,reference_index,2,grid.dX.y);
      psi_zk=interpolant.Calculate_Interpolant(index,reference_index,3,grid.dX.z);
      error=fabs(phi_k-phi_Goal);
      X_k=root(index);
      if(error<criterion)
	quitloop=1;

      iter++;
    }

    // Update convergence flag
    if(fabs(error)<=criterion){
      Convergence_Flag(index)=1;
      cout<<"Converged by moving along instantaneous normal"<<endl;
    }
    else if(fabs(error)>criterion && iter==MaxIter){
      Convergence_Flag(index)=0;
      cout<<"Not converged by moving along instantaneous normal"<<endl;
    }
    else{}
  }
  //#########END
   // END COMMENT

  //#########IF NOT CONVERGED ABOVE, WE FIND BY USING LINE SEARCH ALGORITHM
  if(!Convergence_Flag(index)){
    T alpha=1,alphaV=1;
    int quitloop=0;
    T phi_Goal=.0,phi_S=phi(index);
    error=1,iter=1;
    X_k=grid.X(index);
    root(index)=X_k;
    phi_k=phi(index),psi_xk=psi(index).x,psi_yk=psi(index).y,psi_zk=psi(index).z;
    psi_k(1)=psi_xk,psi_k(2)=psi_yk,psi_k(3)=psi_zk;

    while(quitloop!=1 && iter<MaxIter){
      root(index)=X_k+alpha*(phi_Goal-phi_S)*psi_k/psi_k.Magnitude();

      reference_index=grid.Cell(root(index),0);

      phi_k=interpolant.Calculate_Interpolant(index,reference_index,0,(T).0);
      psi_xk=interpolant.Calculate_Interpolant(index,reference_index,1,grid.dX.x);
      psi_yk=interpolant.Calculate_Interpolant(index,reference_index,2,grid.dX.y);
      psi_zk=interpolant.Calculate_Interpolant(index,reference_index,3,grid.dX.z);

      if(phi(index)>0){
	error=fabs(phi_k-phi_Goal);
	if(error<=criterion)
	  quitloop=1;
	else if(phi_k<0)
	  alpha/=(T)2.;
	else if(phi_k>0){
	  X_k=root(index);
	  phi_S=phi_k;
	}
	else{}
      }
      else if(phi(index)<0){
	error=fabs(phi_k-phi_Goal);
	if(error<=criterion)
	  quitloop=1;
	else if(phi_k>0)
	  alpha/=(T)2.;
	else if(phi_k<0){
	  X_k=root(index);
	  phi_S=phi_k;
	}
	else{}
      }

      iter++;
    }

    // Update convergence flag
    if(fabs(error)<=criterion){
      Convergence_Flag(index)=1;
      cout<<"Converged by using line search algorithm"<<endl;
    }
    else if(fabs(error)>criterion && iter==MaxIter){
      Convergence_Flag(index)=0;
      cout<<"Not converged by using line search algorithm"<<endl;
    }
    else{}
  }
  //#########END

  // Updating psi at interfacial rectangles (nodes)
  reference_index=grid.Cell(root(index),0);

  phi_k=interpolant.Calculate_Interpolant(index,reference_index,0,(T).0);
  psi_xk=interpolant.Calculate_Interpolant(index,reference_index,1,grid.dX.x);
  psi_yk=interpolant.Calculate_Interpolant(index,reference_index,2,grid.dX.y);
  psi_zk=interpolant.Calculate_Interpolant(index,reference_index,3,grid.dX.z);

  Rphi(index)=sign(phi(index))*(X-root(index)).Magnitude();
  Rpsi(index).x=psi_xk;
  Rpsi(index).y=psi_yk;
  Rpsi(index).z=psi_zk;
  Rpsi(index)/=Rpsi(index).Magnitude();
}

// Locate_Interfacial_Rectangles
template<class T>
void Reinitialize<T>::Locate_Interfacial_Rectangles(const T_ARRAYS& phi,T_ARRAYS& Iphi)
{
  N_Rectangles=0;
  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    Interfacial_Rectangles(index)=0.;
    if(phi(index)*Iphi(index)<0){
      N_Rectangles++;
      Interfacial_Rectangles(index)=1.;
    }
  }
}

template class Reinitialize<float>;
