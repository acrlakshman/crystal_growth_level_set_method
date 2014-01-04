//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "HermiteInterpolation.h"

using namespace std;
using namespace PhysBAM;

/*****************************************************/
// Constructor
/*****************************************************/
template<class T>
HermiteInterpolation<T>::HermiteInterpolation(const T_GRID& grid_, TV_ARRAYS& root_, T_ARRAYS& phi_, TV_ARRAYS& psi_)
  :grid(grid_),root(root_),phi(phi_),psi(psi_)
{
  x_cells=grid.Numbers_Of_Nodes().x;
  y_cells=grid.Numbers_Of_Nodes().y;
  z_cells=grid.Numbers_Of_Nodes().z;
  Iphi=phi;
  Ipsi=psi;
  psi_xy=phi;
  psi_xz=phi;
  psi_yz=phi;
  psi_xyz=phi;
}

/*****************************************************/
// Destructor
/*****************************************************/
template<class T>
HermiteInterpolation<T>::~HermiteInterpolation()
{}

/*****************************************************/
// Function: Interpolate_Phi
/*****************************************************/
template<class T>
void HermiteInterpolation<T>::Interpolate_Phi()
{
  // Compute cross derivatives
  Compute_Cross_Derivative(1); // psi_xy
  Compute_Cross_Derivative(2); // psi_xz
  Compute_Cross_Derivative(3); // psi_yz
  Compute_Cross_Derivative(4); // psi_xyz

  // Iterate through all grid nodes (except the border cells)
  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices().Thickened(-1));iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();

    reference_index=grid.Cell(root(index),0);

    Calculate_Iphi(index,reference_index);
    Calculate_Ipsi(index,reference_index);;
  }
}

// Calculate Iphi at root(index) location
template<class T>
void HermiteInterpolation<T>::Calculate_Iphi(T_INDEX index,T_INDEX reference_index)
{
  int derivative=0;
  T hC=.0;
  Iphi(index)=(T).0;
  Iphi(index)=Calculate_Interpolant(index,reference_index,derivative,hC);
}

// Calculate Ipsi at root(index) location
template<class T>
void HermiteInterpolation<T>::Calculate_Ipsi(T_INDEX index,T_INDEX reference_index)
{
  int derivative=1; //- for 'x' derivative
  T hx=grid.dX.x, hy=grid.dX.y, hz=grid.dX.z;
  T hC=hx; //- for 'x' derivative
  Ipsi(index).x=(T).0;
  Ipsi(index).x=Calculate_Interpolant(index,reference_index,derivative,hC);

  derivative=2; //- for 'y' derivative
  hC=hy; //- for 'y' derivative
  Ipsi(index).y=(T).0;
  Ipsi(index).y=Calculate_Interpolant(index,reference_index,derivative,hC);

  derivative=3; //- for 'z' derivative
  hC=hz; //- for 'z' derivative
  Ipsi(index).z=(T).0;
  Ipsi(index).z=Calculate_Interpolant(index,reference_index,derivative,hC);
}

// Calculate Iphi at root(index) location
template<class T>
T HermiteInterpolation<T>::Calculate_Interpolant(T_INDEX index,T_INDEX reference_index,int derivative,T hC)
{
  T_INDEX index_x=reference_index+T_INDEX::Axis_Vector(1);
  T_INDEX index_y=reference_index+T_INDEX::Axis_Vector(2);
  T_INDEX index_z=reference_index+T_INDEX::Axis_Vector(3);
  T_INDEX index_xy=reference_index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2);
  T_INDEX index_xz=reference_index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(3);
  T_INDEX index_yz=reference_index+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3);
  T_INDEX index_xyz=reference_index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3);
  T hx=grid.dX.x, hy=grid.dX.y, hz=grid.dX.z;
  T hxy=hx*hy, hxz=hx*hz, hyz=hy*hz, hxyz=hx*hy*hz;
  T etax=(root(index)(1)-grid.X(reference_index)(1))/hx;
  T etay=(root(index)(2)-grid.X(reference_index)(2))/hy;
  T etaz=(root(index)(3)-grid.X(reference_index)(3))/hz;

  T interpolated_value=(T).0;
  counterg=0;

  interpolated_value=(phi(reference_index)*BasisPoly(0,0,0,0,0,0,derivative,hC,etax,etay,etaz)) +
    (phi(index_z)*BasisPoly(0,0,1,0,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(reference_index).z*hz*BasisPoly(0,0,0,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_z).z*hz*BasisPoly(0,0,1,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(reference_index).y*hy*BasisPoly(0,0,0,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_z).y*hy*BasisPoly(0,0,1,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_yz(reference_index)*hyz*BasisPoly(0,0,0,0,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_yz(index_z)*hyz*BasisPoly(0,0,1,0,1,1,derivative,hC,etax,etay,etaz)) +
    (phi(index_y)*BasisPoly(0,1,0,0,0,0,derivative,hC,etax,etay,etaz)) +
    (phi(index_yz)*BasisPoly(0,1,1,0,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_y).z*hz*BasisPoly(0,1,0,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_yz).z*hz*BasisPoly(0,1,1,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_y).y*hy*BasisPoly(0,1,0,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_yz).y*hy*BasisPoly(0,1,1,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_yz(index_y)*hyz*BasisPoly(0,1,0,0,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_yz(index_yz)*hyz*BasisPoly(0,1,1,0,1,1,derivative,hC,etax,etay,etaz)) +
    (psi(reference_index).x*hx*BasisPoly(0,0,0,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_z).x*hx*BasisPoly(0,0,1,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi_xz(reference_index)*hxz*BasisPoly(0,0,0,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xz(index_z)*hxz*BasisPoly(0,0,1,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xy(reference_index)*hxy*BasisPoly(0,0,0,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xy(index_z)*hxy*BasisPoly(0,0,1,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(reference_index)*hxyz*BasisPoly(0,0,0,1,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(index_z)*hxyz*BasisPoly(0,0,1,1,1,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_y).x*hx*BasisPoly(0,1,0,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_yz).x*hx*BasisPoly(0,1,1,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi_xz(index_y)*hxz*BasisPoly(0,1,0,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xz(index_yz)*hxz*BasisPoly(0,1,1,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xy(index_y)*hxy*BasisPoly(0,1,0,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xy(index_yz)*hxy*BasisPoly(0,1,1,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(index_y)*hxyz*BasisPoly(0,1,0,1,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(index_yz)*hxyz*BasisPoly(0,1,1,1,1,1,derivative,hC,etax,etay,etaz)) +
    (phi(index_x)*BasisPoly(1,0,0,0,0,0,derivative,hC,etax,etay,etaz)) +
    (phi(index_xz)*BasisPoly(1,0,1,0,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_x).z*hz*BasisPoly(1,0,0,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_xz).z*hz*BasisPoly(1,0,1,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_x).y*hy*BasisPoly(1,0,0,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_xz).y*hy*BasisPoly(1,0,1,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_yz(index_x)*hyz*BasisPoly(1,0,0,0,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_yz(index_xz)*hyz*BasisPoly(1,0,1,0,1,1,derivative,hC,etax,etay,etaz)) +
    (phi(index_xy)*BasisPoly(1,1,0,0,0,0,derivative,hC,etax,etay,etaz)) +
    (phi(index_xyz)*BasisPoly(1,1,1,0,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_xy).z*hz*BasisPoly(1,1,0,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_xyz).z*hz*BasisPoly(1,1,1,0,0,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_xy).y*hy*BasisPoly(1,1,0,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_xyz).y*hy*BasisPoly(1,1,1,0,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_yz(index_xy)*hyz*BasisPoly(1,1,0,0,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_yz(index_xyz)*hyz*BasisPoly(1,1,1,0,1,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_x).x*hx*BasisPoly(1,0,0,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_xz).x*hx*BasisPoly(1,0,1,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi_xz(index_x)*hxz*BasisPoly(1,0,0,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xz(index_xz)*hxz*BasisPoly(1,0,1,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xy(index_x)*hxy*BasisPoly(1,0,0,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xy(index_xz)*hxy*BasisPoly(1,0,1,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(index_x)*hxyz*BasisPoly(1,0,0,1,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(index_xz)*hxyz*BasisPoly(1,0,1,1,1,1,derivative,hC,etax,etay,etaz)) +
    (psi(index_xy).x*hx*BasisPoly(1,1,0,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi(index_xyz).x*hx*BasisPoly(1,1,1,1,0,0,derivative,hC,etax,etay,etaz)) +
    (psi_xz(index_xy)*hxz*BasisPoly(1,1,0,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xz(index_xyz)*hxz*BasisPoly(1,1,1,1,0,1,derivative,hC,etax,etay,etaz)) +
    (psi_xy(index_xy)*hxy*BasisPoly(1,1,0,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xy(index_xyz)*hxy*BasisPoly(1,1,1,1,1,0,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(index_xy)*hxyz*BasisPoly(1,1,0,1,1,1,derivative,hC,etax,etay,etaz)) +
    (psi_xyz(index_xyz)*hxyz*BasisPoly(1,1,1,1,1,1,derivative,hC,etax,etay,etaz));

  return interpolated_value;
}

// component=1(psi_xy); 2(psi_xz); 3(psi_yz); 4(psi_xyz)
template<class T>
void HermiteInterpolation<T>::Compute_Cross_Derivative(int component)
{
  switch (component)
    {
    case 1: //psi_xy
      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	psi_xy(index)=T();
	if(index(1)==1){
	  psi_xy(index)+=(psi(index+T_INDEX::Axis_Vector(1)).y-psi(index).y)*(grid.One_Over_DX()(1));
	}
	else if(index(1)==x_cells){
	  psi_xy(index)+=(psi(index).y-psi(index-T_INDEX::Axis_Vector(1)).y)*(grid.One_Over_DX()(1));
	}
	else{
	  psi_xy(index)+=(psi(index+T_INDEX::Axis_Vector(1)).y-psi(index-T_INDEX::Axis_Vector(1)).y)*(T)0.5*(grid.One_Over_DX()(1));
	}
      }
      break;

    case 2: //psi_xz
      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	psi_xz(index)=T();
	if(index(3)==1){
	  psi_xz(index)+=(psi(index+T_INDEX::Axis_Vector(3)).x-psi(index).x)*(grid.One_Over_DX()(3));
	}
	else if(index(3)==z_cells){
	  psi_xz(index)+=(psi(index).x-psi(index-T_INDEX::Axis_Vector(3)).x)*(grid.One_Over_DX()(3));
	}
	else{
	  psi_xz(index)+=(psi(index+T_INDEX::Axis_Vector(3)).x-psi(index-T_INDEX::Axis_Vector(3)).x)*(T)0.5*(grid.One_Over_DX()(3));
	}
      }
      break;

    case 3: //psi_yz
      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	psi_yz(index)=T();
	if(index(2)==1){
	  psi_yz(index)+=(psi(index+T_INDEX::Axis_Vector(2)).z-psi(index).z)*(grid.One_Over_DX()(2));
	}
	else if(index(2)==y_cells){
	  psi_yz(index)+=(psi(index).z-psi(index-T_INDEX::Axis_Vector(2)).z)*(grid.One_Over_DX()(2));
	}
	else{
	  psi_yz(index)+=(psi(index+T_INDEX::Axis_Vector(2)).z-psi(index-T_INDEX::Axis_Vector(2)).z)*(T)0.5*(grid.One_Over_DX()(2));
	}
      }
      break;

    case 4: //psi_xyz
      T One_Over_dxdy=grid.One_Over_DX()(1)*grid.One_Over_DX()(2);
      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	psi_xyz(index)=T();
	if(index(2)==1){
	  if(index(1)==1){
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(1)).z+psi(index).z)*One_Over_dxdy;
	  }
	  else if(index(1)==x_cells){
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(2)-T_INDEX::Axis_Vector(1)).z-psi(index).z+psi(index-T_INDEX::Axis_Vector(1)).z)*One_Over_dxdy;
	  }
	  else{
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(2)-T_INDEX::Axis_Vector(1)).z-psi(index+T_INDEX::Axis_Vector(1)).z+psi(index-T_INDEX::Axis_Vector(1)).z)*(T)0.5*One_Over_dxdy;
	  }
	}
	else if(index(2)==y_cells){
	  if(index(1)==1){
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(1)).z-psi(index).z-psi(index+T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z+psi(index-T_INDEX::Axis_Vector(2)).z)*One_Over_dxdy;
	  }
	  else if(index(1)==x_cells){
	    psi_xyz(index)+=(psi(index).z-psi(index-T_INDEX::Axis_Vector(1)).z-psi(index-T_INDEX::Axis_Vector(2)).z+psi(index-T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z)*One_Over_dxdy;
	  }
	  else{
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(1)).z-psi(index-T_INDEX::Axis_Vector(1)).z-psi(index+T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z+psi(index-T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z)*(T)0.5*One_Over_dxdy;
	  }
	}
	else{
	  if(index(1)==1){
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z+psi(index-T_INDEX::Axis_Vector(2)).z)*(T).5*One_Over_dxdy;
	  }
	  else if(index(1)==x_cells){
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(2)-T_INDEX::Axis_Vector(1)).z-psi(index-T_INDEX::Axis_Vector(2)).z+psi(index-T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z)*(T).5*One_Over_dxdy;
	  }
	  else{
	    psi_xyz(index)+=(psi(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)).z-psi(index-T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)).z-psi(index+T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z+psi(index-T_INDEX::Axis_Vector(1)-T_INDEX::Axis_Vector(2)).z)*(T)0.25*One_Over_dxdy;
	  }
	}
      }
      break;
    }
}

// basis polynomial: B^v_d.
// derivative:- 'none':0, 'x':1, 'y':2, 'z':3
template<class T>
T HermiteInterpolation<T>::BasisPoly(int v1,int v2,int v3,int d1,int d2,int d3,int derivative,T hC,T etax,T etay,T etaz)
{
  counterg++;

  T poly=.0;

  if (derivative==0){ // derivative='none'
    if (v1==0){
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    }
    else{ // if (v1==1)
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    } // end if (v1==0)-else
  }
  else if (derivative==1){ // derivative='x'
     if (v1==0){
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)1./hC)*CubicPoly(etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)1./hC)*CubicPoly(etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    }
    else{ // if (v1==1)
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*CubicPoly((T)1.-etax,3)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly(etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly((T)1.-etay,1)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly((T)1.-etay,1)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*CubicPoly((T)1.-etay,1)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etax,4)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    } // end if (v1==0)-else
  }
  else if (derivative==2){ // derivative='y'
     if (v1==0){
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)1./hC)*CubicPoly(etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)1./hC)*CubicPoly(etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    }
    else{ // if (v1==1)
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)1./hC)*CubicPoly(etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)1./hC)*CubicPoly(etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*CubicPoly((T)1.-etay,3)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,1);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly((T)1.-etaz,1);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*CubicPoly(etaz,2);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etay,4)*((T)-1.)*CubicPoly((T)1.-etaz,2);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    } // end if (v1==0)-else
  }
  else{ // derivative='z'
     if (v1==0){
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly(etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,1)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly(etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,2)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly(etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    }
    else{ // if (v1==1)
      if (d1==0){
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly(etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=CubicPoly((T)1.-etax,1)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      }
      else{ // if (d1==1)
	if (v2==0){
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly(etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	}
	else{ // if (v2==1)
	  if (d2==0){
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*CubicPoly((T)1.-etay,1)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  }
	  else{ // if (d2==1)
	    if (d3==0){
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,3);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*CubicPoly((T)1.-etaz,3);
	      } // end if (v3==0)-else
	    }
	    else{ // if (d3==1)
	      if (v3==0){
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)1./hC)*CubicPoly(etaz,4);
	      }
	      else{
		poly=((T)-1.)*CubicPoly((T)1.-etax,2)*((T)-1.)*CubicPoly((T)1.-etay,2)*((T)-1./hC)*((T)-1.)*CubicPoly((T)1.-etaz,4);
	      } // end if (v3==0)-else
	    } // end if (d3==0)-else
	  } // end if (d2==0)-else
	} // end if (v2==0)-else
      } // end if (d1==0)-else
    } // end if (v1==0)-else
  }

  return poly;
}

// Cubic polynomials. 'f':1, 'g':2, 'fd':3, 'gd':4
template<class T>
T HermiteInterpolation<T>::CubicPoly(T eta,int func)
{
  if (func==1)
    return ( (T)1.-((T)3.*eta*eta)+((T)2.*eta*eta*eta) );
  else if (func==2)
    return ( eta*((T)1.-eta)*((T)1.-eta) );
  else if (func==3)
    return ( -(T)6.*eta + (T)6.*eta*eta );
  else
    return ( ((T)1.-eta)*((T)1.-eta) - ((T)2.*eta*((T)1.-eta)) );
}

template class HermiteInterpolation<float>;
