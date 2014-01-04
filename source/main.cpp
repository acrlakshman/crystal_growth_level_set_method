//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Title: Weak Formulation of Gradient Augmented Level set method to Stephan type problems
// Course: CS838 - Fall 2012
// Instructor: Prof. Eftychios Sifakis ( pages.cs.wisc.edu/~sifakis/ )
// University of Wisconsin Madison
//#####################################################################

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include "RANGE_ITERATOR.h"
#include "HermiteInterpolation.h"
#include "Reinitialize.h"
#include "Extrapolate.h"

#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include "OBJ_FILE.h"

#include <ctime>
//#include <algorithm>

using namespace std;
using namespace PhysBAM;

/*****************************************************/
// Functions: min and max
/*****************************************************/
inline float min(float a, float b) { return a>=b ? b : a; }
inline float max(float a, float b) { return a>=b ? a : b; }

/*****************************************************/
// Function: Curvature
/*****************************************************/

template<class T,int d>
void Compute_Curvature(const GRID<VECTOR<T,d> >& grid,const ARRAY<VECTOR<T,d>,VECTOR<int,d> >&psi,ARRAY<T,VECTOR<int,d> >& curvature)
{
  //cout<<"Computing Gradient"<<endl;
  typedef VECTOR<int,d> T_INDEX;
  typedef ARRAY<T,T_INDEX> T_ARRAYS;
  typedef VECTOR<T,d> TV;
  typedef ARRAY<TV,T_INDEX> TV_ARRAYS;

  T_ARRAYS normalX(grid.Domain_Indices());
  T_ARRAYS normalY(grid.Domain_Indices());
  T_ARRAYS normalZ(grid.Domain_Indices());
  T_ARRAYS normalX_x(grid.Domain_Indices());
  T_ARRAYS normalY_y(grid.Domain_Indices());
  T_ARRAYS normalZ_z(grid.Domain_Indices());

  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    normalX(index)=psi(index).x/psi(index).Magnitude();
    normalY(index)=psi(index).y/psi(index).Magnitude();
    normalZ(index)=psi(index).z/psi(index).Magnitude();
  }
  normalX_x=normalX; // initialization only
  normalY_y=normalY; // initialization only
  normalZ_z=normalZ; // initialization only

  Calculate_Gradient(grid,normalX,normalX_x,1,2);
  Calculate_Gradient(grid,normalY,normalY_y,2,2);
  Calculate_Gradient(grid,normalZ,normalZ_z,3,2);

  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
    const T_INDEX& index=iterator.Index();
    curvature(index)=T();
    curvature(index)=normalX_x(index)+normalY_y(index)+normalZ_z(index);
  }
}

/*****************************************************/
// Function: Inner product (of tensor and vector)
/*****************************************************/
template<class T,int d>
void Compute_InnerProduct(const GRID<VECTOR<T,d> >& grid,char Type1,char Type2,const ARRAY<VECTOR<T,d>,VECTOR<int,d> >& velocity_x_dxdydz,const ARRAY<VECTOR<T,d>,VECTOR<int,d> >& velocity_y_dxdydz,const ARRAY<VECTOR<T,d>,VECTOR<int,d> >& velocity_z_dxdydz,const ARRAY<VECTOR<T,d>,VECTOR<int,d> >& psi,ARRAY<VECTOR<T,d>,VECTOR<int,d> >& GradUDotPsi)
{
  typedef VECTOR<int,d> T_INDEX;
  if(Type1=='T' && Type2=='V'){ // Dot product of tensor and vector
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      GradUDotPsi(index).x=velocity_x_dxdydz(index).x*psi(index).x+velocity_y_dxdydz(index).x*psi(index).y+velocity_z_dxdydz(index).x*psi(index).z;
      GradUDotPsi(index).y=velocity_x_dxdydz(index).y*psi(index).x+velocity_y_dxdydz(index).y*psi(index).y+velocity_z_dxdydz(index).y*psi(index).z;
      GradUDotPsi(index).z=velocity_x_dxdydz(index).z*psi(index).x+velocity_y_dxdydz(index).z*psi(index).y+velocity_z_dxdydz(index).z*psi(index).z;
    }
  }
}

/*****************************************************/
// Function: Gradient
/*****************************************************/
// derivative: 1-'x', 2-'y', 3-'z'
template<class T,int d>
void Compute_Gradient(const GRID<VECTOR<T,d> >& grid,const ARRAY<T,VECTOR<int,d> >& phi,ARRAY<VECTOR<T,d>,VECTOR<int,d> >&psi,int order)
{
  //cout<<"Computing Gradient"<<endl;
  typedef VECTOR<int,d> T_INDEX;
  typedef ARRAY<T,T_INDEX> T_ARRAYS;

  T_ARRAYS psi_x,psi_y,psi_z;
  psi_x=phi;
  psi_y=phi;
  psi_z=phi;

  // Computing to 2nd order accuracy
  Calculate_Gradient(grid,phi,psi_x,1,2); // psi_x
  Calculate_Gradient(grid,phi,psi_y,2,2); // psi_y
  Calculate_Gradient(grid,phi,psi_z,3,2); // psi_z
  if(order==3){
    Calculate_Gradient(grid,phi,psi_x,1,order); // psi_x
    Calculate_Gradient(grid,phi,psi_y,2,order); // psi_y
    Calculate_Gradient(grid,phi,psi_z,3,order); // psi_z
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
template<class T,int d>
void Calculate_Gradient(const GRID<VECTOR<T,d> >& grid,const ARRAY<T,VECTOR<int,d> >& phi,ARRAY<T,VECTOR<int,d> >&psi_p,int component,int order)
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
	    psi_p(index)+=(phi(index+T_INDEX::Axis_Vector(1))-phi(index))*(grid.One_Over_DX()(1));
	  }
	  else if(index(1)==x_cells){
	    psi_p(index)+=(phi(index)-phi(index-T_INDEX::Axis_Vector(1)))*(grid.One_Over_DX()(1));
	  }
	  else{
	    psi_p(index)+=(phi(index+T_INDEX::Axis_Vector(1))-phi(index-T_INDEX::Axis_Vector(1)))*(T)0.5*(grid.One_Over_DX()(1));
	  }
	}
      }
      else{
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  if(index(1)==2){
	    psi_p(index)=T();
	    psi_p(index)+=((-phi(index+T_INDEX::Axis_Vector(1)*2))+((T)6.*phi(index+T_INDEX::Axis_Vector(1)))-((T)3.*phi(index))-((T)2.*phi(index-T_INDEX::Axis_Vector(1))))*((T)(1./6.))*(grid.One_Over_DX()(1));
	  }
	  else if(index(1)>2 && index(1)<=(x_cells-1)){
	    psi_p(index)=T();
	    psi_p(index)+=(((T)2.*phi(index+T_INDEX::Axis_Vector(1)))+((T)3.*phi(index))-((T)6.*phi(index-T_INDEX::Axis_Vector(1)))+(phi(index-T_INDEX::Axis_Vector(1)*2)))*((T)(1./6.))*(grid.One_Over_DX()(1));
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
	    psi_p(index)+=(phi(index+T_INDEX::Axis_Vector(2))-phi(index))*(grid.One_Over_DX()(2));
	  }
	  else if(index(2)==y_cells){
	    psi_p(index)+=(phi(index)-phi(index-T_INDEX::Axis_Vector(2)))*(grid.One_Over_DX()(2));
	  }
	  else{
	    psi_p(index)+=(phi(index+T_INDEX::Axis_Vector(2))-phi(index-T_INDEX::Axis_Vector(2)))*(T)0.5*(grid.One_Over_DX()(2));
	  }
	}
      }
      else{
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  if(index(2)==2){
	    psi_p(index)=T();
	    psi_p(index)+=((-phi(index+T_INDEX::Axis_Vector(2)*2))+((T)6.*phi(index+T_INDEX::Axis_Vector(2)))-((T)3.*phi(index))-((T)2.*phi(index-T_INDEX::Axis_Vector(2))))*((T)(1./6.))*(grid.One_Over_DX()(2));
	  }
	  else if(index(2)>2 && index(2)<=(y_cells-1)){
	    psi_p(index)=T();
	    psi_p(index)+=(((T)2.*phi(index+T_INDEX::Axis_Vector(2)))+((T)3.*phi(index))-((T)6.*phi(index-T_INDEX::Axis_Vector(2)))+(phi(index-T_INDEX::Axis_Vector(2)*2)))*((T)(1./6.))*(grid.One_Over_DX()(2));
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
	    psi_p(index)+=(phi(index+T_INDEX::Axis_Vector(3))-phi(index))*(grid.One_Over_DX()(3));
	  }
	  else if(index(3)==z_cells){
	    psi_p(index)+=(phi(index)-phi(index-T_INDEX::Axis_Vector(3)))*(grid.One_Over_DX()(3));
	  }
	  else{
	    psi_p(index)+=(phi(index+T_INDEX::Axis_Vector(3))-phi(index-T_INDEX::Axis_Vector(3)))*(T)0.5*(grid.One_Over_DX()(3));
	  }
	}
      }
      else{
	for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	  const T_INDEX& index=iterator.Index();
	  if(index(3)==2){
	    psi_p(index)=T();
	    psi_p(index)+=((-phi(index+T_INDEX::Axis_Vector(3)*2))+((T)6.*phi(index+T_INDEX::Axis_Vector(3)))-((T)3.*phi(index))-((T)2.*phi(index-T_INDEX::Axis_Vector(3))))*((T)(1./6.))*(grid.One_Over_DX()(3));
	  }
	  else if(index(3)>2 && index(3)<=(z_cells-1)){
	    psi_p(index)=T();
	    psi_p(index)+=(((T)2.*phi(index+T_INDEX::Axis_Vector(3)))+((T)3.*phi(index))-((T)6.*phi(index-T_INDEX::Axis_Vector(3)))+(phi(index-T_INDEX::Axis_Vector(3)*2)))*((T)(1./6.))*(grid.One_Over_DX()(3));
	  }
	  else{}
	}
      }
      break;
    }
}

/*****************************************************/
// Function: Apply zero-gradient boundary condition for a scalar (phi)
/*****************************************************/
// INCOMPLETE (Not needed if reinitialization is used)
template<class T,int d>
void Zero_Gradient_Boundary_Condition(const GRID<VECTOR<T,d> >& grid,ARRAY<T,VECTOR<int,d> >& x)
{
  typedef VECTOR<int,d> T_INDEX;

  PHYSBAM_ASSERT(grid.Domain_Indices()==x.Domain_Indices());

  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(!grid.Domain_Indices().Thickened(-1).Lazy_Inside(index))
            x(index)=(T).0;
    }
}

/*****************************************************/
// Function: Apply boundary condition for temperature
/*****************************************************/

template<class T,int d>
void Apply_Boundary_Condition(const GRID<VECTOR<T,d> >& grid,ARRAY<T,VECTOR<int,d> >& x)
{
  typedef VECTOR<int,d> T_INDEX;

  PHYSBAM_ASSERT(grid.Domain_Indices()==x.Domain_Indices());

  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(!grid.Domain_Indices().Thickened(-1).Lazy_Inside(index))
            x(index)=(T)-.5;
    }
}

/**********************************************************/
// Function: Apply zero boundary condition for temperature
/**********************************************************/

template<class T,int d>
void Zero_Boundary_Condition(const GRID<VECTOR<T,d> >& grid,ARRAY<T,VECTOR<int,d> >& x)
{
  typedef VECTOR<int,d> T_INDEX;

  PHYSBAM_ASSERT(grid.Domain_Indices()==x.Domain_Indices());

  for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(!grid.Domain_Indices().Thickened(-1).Lazy_Inside(index))
            x(index)=(T).0;
    }
}

/**********************************************************/
// Function Compute_Laplacian
/**********************************************************/
/*
template<class T,int d>
void Compute_Laplacian(const GRID<VECTOR<T,d> >& grid,const ARRAY<T,VECTOR<int,d> >& x,const ARRAY<T,VECTOR<int,d> >& phi,const ARRAY<T,VECTOR<int,d> >& T_interface,ARRAY<T,VECTOR<int,d> >& Lx)
//void Compute_Laplacian(const GRID<VECTOR<T,d> >& grid,const ARRAY<T,VECTOR<int,d> >& x,ARRAY<T,VECTOR<int,d> >& Lx)
{
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(grid.Domain_Indices()==x.Domain_Indices());
    PHYSBAM_ASSERT(grid.Domain_Indices()==Lx.Domain_Indices());

    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices().Thickened(-1));iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        Lx(index)=T();
        for(int v=1;v<=d;v++)
            Lx(index)+=(x(index-T_INDEX::Axis_Vector(v))-(T)2.*x(index)+x(index+T_INDEX::Axis_Vector(v)))*sqr(grid.One_Over_DX()(v));}
}
*/

template<class T,int d>
void Compute_Laplacian(const GRID<VECTOR<T,d> >& grid,const ARRAY<T,VECTOR<int,d> >& x,const ARRAY<T,VECTOR<int,d> >& phi,const ARRAY<T,VECTOR<int,d> >& T_interface,ARRAY<T,VECTOR<int,d> >& Lx)
{
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(grid.Domain_Indices()==x.Domain_Indices());
    PHYSBAM_ASSERT(grid.Domain_Indices()==Lx.Domain_Indices());

    // ###Using linear extrapolation
    /* // START COMMENT
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices().Thickened(-1));iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        Lx(index)=T();
        for(int v=1;v<=d;v++){
	  if(phi(index)*phi(index+T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index+T_INDEX::Axis_Vector(v)))*grid.X(index)(v))+(fabs(phi(index))*grid.X(index+T_INDEX::Axis_Vector(v))(v)))/(fabs(phi(index+T_INDEX::Axis_Vector(v)))+fabs(phi(index)));

	    T x_frac=(x_I-grid.X(index)(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index+T_INDEX::Axis_Vector(v))-T_interface(index))*x_I*grid.One_Over_DX()(v) + (T_interface(index)*grid.X(index+T_INDEX::Axis_Vector(v))(v) - T_interface(index+T_INDEX::Axis_Vector(v))*grid.X(index)(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-0.5;

	    Lx(index)+=(x(index-T_INDEX::Axis_Vector(v))-(x(index)*((T)1.+((T)1./x_frac))))*sqr(grid.One_Over_DX()(v));
	  }
	  else if(phi(index)*phi(index-T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index))*grid.X(index-T_INDEX::Axis_Vector(v))(v))+(fabs(phi(index-T_INDEX::Axis_Vector(v)))*grid.X(index)(v)))/(fabs(phi(index))+fabs(phi(index-T_INDEX::Axis_Vector(v))));

	    T x_frac=(x_I-grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index)-T_interface(index-T_INDEX::Axis_Vector(v)))*x_I*grid.One_Over_DX()(v) + (T_interface(index-T_INDEX::Axis_Vector(v))*grid.X(index)(v) - T_interface(index)*grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-.5;

	    Lx(index)+=(x(index+T_INDEX::Axis_Vector(v))-(x(index)*(((T)2.-x_frac)/((T)1.-x_frac))))*sqr(grid.One_Over_DX()(v));
	  }
	  else{
	    Lx(index)+=(x(index-T_INDEX::Axis_Vector(v))-(T)2.*x(index)+x(index+T_INDEX::Axis_Vector(v)))*sqr(grid.One_Over_DX()(v));
	  }
	}
    }
    */ // END COMMENT

    // ### Using constant extrapolation
     // START COMMENT
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices().Thickened(-1));iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        Lx(index)=T();
        for(int v=1;v<=d;v++){
	  if(phi(index)*phi(index+T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index+T_INDEX::Axis_Vector(v)))*grid.X(index)(v))+(fabs(phi(index))*grid.X(index+T_INDEX::Axis_Vector(v))(v)))/(fabs(phi(index+T_INDEX::Axis_Vector(v)))+fabs(phi(index)));

	    T x_frac=(x_I-grid.X(index)(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index+T_INDEX::Axis_Vector(v))-T_interface(index))*x_I*grid.One_Over_DX()(v) + (T_interface(index)*grid.X(index+T_INDEX::Axis_Vector(v))(v) - T_interface(index+T_INDEX::Axis_Vector(v))*grid.X(index)(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-0.5;

	    Lx(index)+=(x(index-T_INDEX::Axis_Vector(v))-(T)2.*x(index))*sqr(grid.One_Over_DX()(v));
	  }
	  else if(phi(index)*phi(index-T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index))*grid.X(index-T_INDEX::Axis_Vector(v))(v))+(fabs(phi(index-T_INDEX::Axis_Vector(v)))*grid.X(index)(v)))/(fabs(phi(index))+fabs(phi(index-T_INDEX::Axis_Vector(v))));

	    T x_frac=(x_I-grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index)-T_interface(index-T_INDEX::Axis_Vector(v)))*x_I*grid.One_Over_DX()(v) + (T_interface(index-T_INDEX::Axis_Vector(v))*grid.X(index)(v) - T_interface(index)*grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-.5;

	    Lx(index)+=(x(index+T_INDEX::Axis_Vector(v))-(T)2.*x(index))*sqr(grid.One_Over_DX()(v));
	  }
	  else{
	    Lx(index)+=(x(index-T_INDEX::Axis_Vector(v))-(T)2.*x(index)+x(index+T_INDEX::Axis_Vector(v)))*sqr(grid.One_Over_DX()(v));
	  }
	}
    }
    // END COMMENT

    // ### Default Laplacian
    /* // START COMMENT
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices().Thickened(-1));iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        Lx(index)=T();
        for(int v=1;v<=d;v++)
	  Lx(index)+=(x(index-T_INDEX::Axis_Vector(v))-(T)2.*x(index)+x(index+T_INDEX::Axis_Vector(v)))*sqr(grid.One_Over_DX()(v));
    }
    */ // END COMMENT
}

/**********************************************************/
// Function: Compute_rhs (for backward implicit Euler)
/**********************************************************/

template<class T,int d>
void Compute_rhs(const GRID<VECTOR<T,d> >& grid,const ARRAY<T,VECTOR<int,d> >& x,const ARRAY<T,VECTOR<int,d> >& phi,const ARRAY<T,VECTOR<int,d> >& T_interface,ARRAY<T,VECTOR<int,d> >& Lx)
{
    typedef VECTOR<int,d> T_INDEX;

    PHYSBAM_ASSERT(grid.Domain_Indices()==x.Domain_Indices());
    PHYSBAM_ASSERT(grid.Domain_Indices()==Lx.Domain_Indices());

    
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices().Thickened(-1));iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        Lx(index)=T();

	// ### Using constant extrapolation
	// START COMMENT
        for(int v=1;v<=d;v++){
	  if(phi(index)*phi(index+T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index+T_INDEX::Axis_Vector(v)))*grid.X(index)(v))+(fabs(phi(index))*grid.X(index+T_INDEX::Axis_Vector(v))(v)))/(fabs(phi(index+T_INDEX::Axis_Vector(v)))+fabs(phi(index)));

	    T x_frac=(x_I-grid.X(index)(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index+T_INDEX::Axis_Vector(v))-T_interface(index))*x_I*grid.One_Over_DX()(v) + (T_interface(index)*grid.X(index+T_INDEX::Axis_Vector(v))(v) - T_interface(index+T_INDEX::Axis_Vector(v))*grid.X(index)(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-0.5;

	    Lx(index)+=T_int*sqr(grid.One_Over_DX()(v));
	  }
	  else if(phi(index)*phi(index-T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index))*grid.X(index-T_INDEX::Axis_Vector(v))(v))+(fabs(phi(index-T_INDEX::Axis_Vector(v)))*grid.X(index)(v)))/(fabs(phi(index))+fabs(phi(index-T_INDEX::Axis_Vector(v))));

	    T x_frac=(x_I-grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index)-T_interface(index-T_INDEX::Axis_Vector(v)))*x_I*grid.One_Over_DX()(v) + (T_interface(index-T_INDEX::Axis_Vector(v))*grid.X(index)(v) - T_interface(index)*grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-.5;

	    Lx(index)+=T_int*sqr(grid.One_Over_DX()(v));
	  }
	  else{
	  }
	}
	// END COMMENT

	// ### Using linear extrapolation
	/* // START COMMENT
        for(int v=1;v<=d;v++){
	  if(phi(index)*phi(index+T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index+T_INDEX::Axis_Vector(v)))*grid.X(index)(v))+(fabs(phi(index))*grid.X(index+T_INDEX::Axis_Vector(v))(v)))/(fabs(phi(index+T_INDEX::Axis_Vector(v)))+fabs(phi(index)));

	    T x_frac=(x_I-grid.X(index)(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index+T_INDEX::Axis_Vector(v))-T_interface(index))*x_I*grid.One_Over_DX()(v) + (T_interface(index)*grid.X(index+T_INDEX::Axis_Vector(v))(v) - T_interface(index+T_INDEX::Axis_Vector(v))*grid.X(index)(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-0.5;

	    Lx(index)+=(T_int/x_frac)*sqr(grid.One_Over_DX()(v));
	  }
	  else if(phi(index)*phi(index-T_INDEX::Axis_Vector(v))<0){
	    T x_I=((fabs(phi(index))*grid.X(index-T_INDEX::Axis_Vector(v))(v))+(fabs(phi(index-T_INDEX::Axis_Vector(v)))*grid.X(index)(v)))/(fabs(phi(index))+fabs(phi(index-T_INDEX::Axis_Vector(v))));

	    T x_frac=(x_I-grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    T T_int=(T_interface(index)-T_interface(index-T_INDEX::Axis_Vector(v)))*x_I*grid.One_Over_DX()(v) + (T_interface(index-T_INDEX::Axis_Vector(v))*grid.X(index)(v) - T_interface(index)*grid.X(index-T_INDEX::Axis_Vector(v))(v))*grid.One_Over_DX()(v);
	    //T_int=(T)-.5;

	    Lx(index)+=(T_int/((T)1. - x_frac))*sqr(grid.One_Over_DX()(v));
	  }
	  else{
	  }
	}
	*/ // END COMMENT
    }
}

/**********************************************************/
// Class CG_VECTOR
/**********************************************************/

template<class T>
class CG_VECTOR
    :public KRYLOV_VECTOR_BASE<T>
{
  typedef VECTOR<int,3> T_INDEX;
  typedef ARRAY<T,T_INDEX> T_ARRAYS;

  typedef KRYLOV_VECTOR_BASE<T> BASE;

  T_ARRAYS& array;

public:
  CG_VECTOR(T_ARRAYS& array_input) : array(array_input) {}

  static T_ARRAYS& Array(BASE& base_array)
  {return ((CG_VECTOR&)(base_array)).array;}

  static const T_ARRAYS& Array(const BASE& base_array)
  {return ((const CG_VECTOR&)(base_array)).array;}

  BASE& operator+=(const BASE& bv)
  {array+=Array(bv);return *this;}

  BASE& operator-=(const BASE& bv)
  {array-=Array(bv);return *this;}

  BASE& operator*=(const T a)
  {array*=a;return *this;}

  void Copy(const T c,const BASE& bv)
  {ARRAY<T,VECTOR<int,3> >::Copy(c,Array(bv),array);}

  void Copy(const T c1,const BASE& bv1,const BASE& bv2)
  {ARRAY<T,VECTOR<int,3> >::Copy(c1,Array(bv1),Array(bv2),array);}
  
  int Raw_Size() const
  {
    PHYSBAM_NOT_IMPLEMENTED();
  }
    
  T& Raw_Get(int i)
  {
    PHYSBAM_NOT_IMPLEMENTED();
  }
  
};

/**********************************************************/
// Class CG_SYSTEM
/**********************************************************/

template<class T>
class CG_SYSTEM
  :public KRYLOV_SYSTEM_BASE<T>
{
  static const int d=3;
  typedef VECTOR<int,d> T_INDEX;
  typedef ARRAY<T,T_INDEX> T_ARRAYS;
  typedef VECTOR<T,d> TV;
  typedef GRID<TV> T_GRID;
  typedef ARRAY<TV,T_INDEX> TV_ARRAYS;

  typedef KRYLOV_SYSTEM_BASE<T> BASE;
  typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;

  const T_GRID grid;

  const T time;
  const T dt;
  const T diffusion_coefficient;
  const T_ARRAYS phi,T_interface;

public:
  CG_SYSTEM(const T_GRID& grid_input,const T time_input,const T dt_input,const T diffusion_coefficient_input,const T_ARRAYS phi_input,const T_ARRAYS T_interface_input)
    :BASE(false,false),grid(grid_input),time(time_input),dt(dt_input),diffusion_coefficient(diffusion_coefficient_input),phi(phi_input),T_interface(T_interface_input) {}

  void Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const
  {
    const T_ARRAYS& v_array=CG_VECTOR<T>::Array(v);
    T_ARRAYS& result_array=CG_VECTOR<T>::Array(result);
    Compute_Laplacian(grid,v_array,phi,T_interface,result_array);
    //Compute_Laplacian(grid,v_array,result_array);
    result_array *= -(T)1.*dt*diffusion_coefficient; // Use .5 for CN method
    result_array += v_array;
  }

  double Inner_Product(const VECTOR_BASE& vx,const VECTOR_BASE& vy) const
  {
    const T_ARRAYS& vx_array=CG_VECTOR<T>::Array(vx);
    const T_ARRAYS& vy_array=CG_VECTOR<T>::Array(vy);

    double result=0.;
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      result+=(vx_array(index)*vy_array(index));
    }
    return result;
  }

  T Convergence_Norm(const VECTOR_BASE& vx) const
  {
    const T_ARRAYS& vx_array=CG_VECTOR<T>::Array(vx);
    T result=0.;
    for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
      const T_INDEX& index=iterator.Index();
      result=std::max(result,fabs(vx_array(index)));
    }
    return result;
  }

  void Project(VECTOR_BASE& vx) const
  {
    T_ARRAYS& vx_array=CG_VECTOR<T>::Array(vx);
    Zero_Boundary_Condition(grid,vx_array);
  }

  void Project_Nullspace(VECTOR_BASE& x) const
  {
    //PHYSBAM_NOT_IMPLEMENTED();
  }

  void Set_Boundary_Conditions(VECTOR_BASE& v) const
  {
    T_ARRAYS& v_array=CG_VECTOR<T>::Array(v);
    Apply_Boundary_Condition(grid,v_array);
  }

};

/**********************************************************/
// Function main
/**********************************************************/

int main(int argc,char* argv[])
{
    STREAM_TYPE stream_type((float)0);

    static const int d=3;
    typedef float T;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef GRID<TV> T_GRID;
    typedef ARRAY<T,T_INDEX> T_ARRAYS;
    typedef ARRAY<TV,T_INDEX> TV_ARRAYS;

    #include "DefineFields.h"

    // Create common directory to write grid
    FILE_UTILITIES::Create_Directory("output/common");
    FILE_UTILITIES::Write_To_File(stream_type,"output/common/grid",grid);
    // Create directory to write temperature
    FILE_UTILITIES::Create_Directory("output/0");
    FILE_UTILITIES::Write_To_File(stream_type,"output/0/temperature",temperature);
    FILE_UTILITIES::Write_To_File(stream_type,"output/0/levelset",grid,Rphi);

    FILE_UTILITIES::Create_Directory("objfiles");
    Write_Obj_File(*surface,"objfiles/levelset_0.obj");

    // Declaring required variables
    T_ARRAYS rhs(grid.Domain_Indices());
    T_ARRAYS temp_q(grid.Domain_Indices());
    T_ARRAYS temp_s(grid.Domain_Indices());
    T_ARRAYS temp_r(grid.Domain_Indices());
    T_ARRAYS temp_k(grid.Domain_Indices());
    T_ARRAYS temp_z(grid.Domain_Indices());

    #include "SetDeltaT.h"

    for(int frame=1;frame<300;frame++){
      // Update time
      time += dt;
      cout<<"Time = "<<time<<endl;

      /***********************************************/
            // START COMMENT
      // Extrapolate temperature field and compute its gradient at interface locations for nodes neighboring interface
      extrapolate_temperature.Extrapolate_Phi();
      temperature_P=extrapolate_temperature.Eu_P_F();
      temperature_N=extrapolate_temperature.Eu_N_F();
      grad_temperature_P=extrapolate_temperature.Egrad_u_P_F();
      grad_temperature_N=extrapolate_temperature.Egrad_u_N_F();

      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	Get_Interface_Location_temperature_P.Interfacial_Rectangles(index)=0;
	//Get_Interface_Location_temperature_N.Interfacial_Rectangles(index)=0;
	if(index(1)!=1 && index(2)!=1 && index(3)!=1 && index(1)!=m && index(2)!=n && index(3)!=l){
	  if(phi(index)*phi(index+T_INDEX::Axis_Vector(1))<0 || phi(index)*phi(index+T_INDEX::Axis_Vector(2))<0 || phi(index)*phi(index+T_INDEX::Axis_Vector(3))<0 || phi(index)*phi(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2))<0 || phi(index)*phi(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(3))<0 || phi(index)*phi(index+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3))<0 || phi(index)*phi(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3))<0){
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index)=1;
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1))=1;
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(2))=1;
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(3))=1;
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2))=1;
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(3))=1;
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3))=1;
	    Get_Interface_Location_temperature_P.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3))=1;

	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index)=1;
	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1))=1;
	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(2))=1;
	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(3))=1;
	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2))=1;
	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(3))=1;
	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3))=1;
	    //Get_Interface_Location_temperature_N.Interfacial_Rectangles(index+T_INDEX::Axis_Vector(1)+T_INDEX::Axis_Vector(2)+T_INDEX::Axis_Vector(3))=1;
	  }
	}
      }
      Get_Interface_Location_temperature_N.Interfacial_Rectangles=Get_Interface_Location_temperature_P.Interfacial_Rectangles;

      // Compute interfacial locations for interfacial nodes
      /* // START COMMENT FOR NEWTON'S METHOD WITH EXTRAPOLATION
      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	if(Get_Interface_Location_temperature_P.Interfacial_Rectangles(index)==1){
	  Get_Interface_Location_temperature_P.Modified_Newtons_Method(interpolant_temperature_P_phi,index);
	  //Get_Interface_Location_temperature_N.Modified_Newtons_Method(interpolant_temperature_N_phi,index);
	}
      }
      Get_Interface_Location_temperature_N.root=Get_Interface_Location_temperature_P.root;

      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	if(Get_Interface_Location_temperature_P.Interfacial_Rectangles(index)==1){
	  const T_INDEX& reference_index=grid.Cell(Get_Interface_Location_temperature_P.root(index),0);

	  interpolant_temperature_P.Iphi(index)=interpolant_temperature_P.Calculate_Interpolant(index,reference_index,0,(T).0);
	  interpolant_temperature_P.Ipsi(index).x=interpolant_temperature_P.Calculate_Interpolant(index,reference_index,1,grid.dX.x);
	  interpolant_temperature_P.Ipsi(index).y=interpolant_temperature_P.Calculate_Interpolant(index,reference_index,2,grid.dX.y);
	  interpolant_temperature_P.Ipsi(index).z=interpolant_temperature_P.Calculate_Interpolant(index,reference_index,3,grid.dX.z);

	  interpolant_temperature_N.Iphi(index)=interpolant_temperature_N.Calculate_Interpolant(index,reference_index,0,(T).0);
	  interpolant_temperature_N.Ipsi(index).x=interpolant_temperature_N.Calculate_Interpolant(index,reference_index,1,grid.dX.x);
	  interpolant_temperature_N.Ipsi(index).y=interpolant_temperature_N.Calculate_Interpolant(index,reference_index,2,grid.dX.y);
	  interpolant_temperature_N.Ipsi(index).z=interpolant_temperature_N.Calculate_Interpolant(index,reference_index,3,grid.dX.z);
	}
      }
      */ // END COMMENT FOR NEWTON'S METHOD WITH EXTRAPOLATION
       // END COMMENT

      /********* Compute T_interface **********/
       // START COMMENT
      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	T_interface(index)=-epsilon_c*Curvature(index);

	// Computing velocity_n and velocity of the interface by obtaining velocities at interface
	/* // START COMMENT
	velocity_n(index)=(interpolant_temperature_P.Ipsi(index).x-interpolant_temperature_N.Ipsi(index).x)*psi(index).x/psi(index).Magnitude()+
	  (interpolant_temperature_P.Ipsi(index).y-interpolant_temperature_N.Ipsi(index).y)*psi(index).y/psi(index).Magnitude()+
	  (interpolant_temperature_P.Ipsi(index).z-interpolant_temperature_N.Ipsi(index).z)*psi(index).z/psi(index).Magnitude();
	velocity_n(index)*=(T)-1.*epsilon_v;

	T_interface(index)+=velocity_n(index);

	// Computing velocity of the interface
	velocity(index)=(interpolant_temperature_P.Ipsi(index
)-interpolant_temperature_N.Ipsi(index));
	*/ // END COMMENT

	// Computing velocity_n and velocity using grad_temperature_P and grad_temperature_N
	// START COMMENT
	velocity_n(index)=(grad_temperature_P(index).x-grad_temperature_N(index).x)*psi(index).x/psi(index).Magnitude()+
	  (grad_temperature_P(index).y-grad_temperature_N(index).y)*psi(index).y/psi(index).Magnitude()+
	  (grad_temperature_P(index).z-grad_temperature_N(index).z)*psi(index).z/psi(index).Magnitude();
	velocity_n(index)*=(T)-1.*epsilon_v;

	T_interface(index)+=velocity_n(index);

	// Computing velocity of the interface
	velocity(index)=grad_temperature_P(index)-grad_temperature_N(index);
	// END COMMENT
      }
       //END COMMENT
      /**********************************************/

      // Adjust dt
      // START COMMENT
      cout<<"Before adjusting, dt="<<dt<<endl;
      dt=(T).5*min(grid.dX.x,grid.dX.y,grid.dX.z);
      T dt_temp=.0;
      for(RANGE_ITERATOR<d> iterator(grid.Domain_Indices());iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	if(Get_Interface_Location_temperature_P.Interfacial_Rectangles(index)==1){
	  dt_temp=(T).5/((velocity(index).x/grid.dX.x)+(velocity(index).y/grid.dX.y)+(velocity(index).z/grid.dX.z));
	  dt=min(dt,fabs(dt_temp));
	}
      }
      cout<<"After adjusting dt="<<dt<<endl;
      // END COMMENT

      /* Solve for diffusion equation (Stefan's problem) */

      // Set boundary conditions at time = t_n
      Apply_Boundary_Condition(grid,temperature);

      // Compute right-hand-side
      //Compute_Laplacian(grid,temperature,phi,T_interface,rhs);
      //Compute_Laplacian(grid,temperature,rhs);
      Compute_rhs(grid,temperature,phi,T_interface,rhs);
      rhs *= (T)1.*dt*diffusion_coefficient; // Use .5 for CN method
      rhs += temperature;

      // Define CG-related vectors
      CG_VECTOR<T> x(temperature_new);
      CG_VECTOR<T> b(rhs);
      CG_VECTOR<T> q(temp_q);
      CG_VECTOR<T> s(temp_s);
      CG_VECTOR<T> r(temp_r);
      CG_VECTOR<T> k(temp_k);
      CG_VECTOR<T> z(temp_z);

      // Scale laplacian to generate right-hand-side

      // Generate CG-formatted system object
      CG_SYSTEM<T> cg_system(grid,time,dt,diffusion_coefficient,phi,T_interface);
      //CG_SYSTEM<T> cg_system(grid,time,dt,diffusion_coefficient);

      // Generate Conjugate Gradients solver object
      CONJUGATE_GRADIENT<T> cg;
      cg.print_residuals=false;
      cg.print_diagnostics=true;

      //Solve linear system using CG
      // CHECK THE CONVERGENCE TOLERANCE
      cg.Solve(cg_system,x,b,q,s,r,k,z,1e-6,0,1000);

      temperature=temperature_new;

      // Overwrite boundary conditions for time = t_n
      Apply_Boundary_Condition(grid,temperature);

      // Advect levelset function
      #include "AdvectLevelSet.h"
      surface=dualcontour.Create_Triangulated_Surface_From_Levelset(levelset);
      Write_Obj_File(*surface,"objfiles/levelset_"+STRING_UTILITIES::Value_To_String(frame)+".obj");

      FILE_UTILITIES::Create_Directory("output/"+STRING_UTILITIES::Value_To_String(frame));
      // Writing temperature to directory
      FILE_UTILITIES::Write_To_File(stream_type,"output/"+STRING_UTILITIES::Value_To_String(frame)+"/temperature",temperature);
      // Writing levelset to directory
      FILE_UTILITIES::Write_To_File(stream_type,"output/"+STRING_UTILITIES::Value_To_String(frame)+"/levelset",grid,Rphi);
    }

    LOG::Initialize_Logging();
    LOG::Finish_Logging();
}
