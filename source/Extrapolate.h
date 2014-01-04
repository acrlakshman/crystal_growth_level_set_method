//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

/*********************************************/
// Extrapolates u from phi<0 to phi>0 and also
// from phi>0 to phi<0...

// LINEAR EXTRAPOLATION and CONSTANT EXTRAPOLATION
// Based on T. Aslam (2003) JCP BUT USING SEMI-LAGRAGIAN APPROACH
/*********************************************/
#pragma once

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include "RANGE_ITERATOR.h"
#include "PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h"
#include "HermiteInterpolation.h"

namespace PhysBAM{

  template<class T>
class Extrapolate
  {
  public:
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef GRID<TV> T_GRID;
    typedef ARRAY<T,T_INDEX> T_ARRAYS;
    typedef ARRAY<TV,T_INDEX> TV_ARRAYS;

    const T_GRID& grid;
    const T_ARRAYS& u;
    const T_ARRAYS& phi;
    const TV_ARRAYS& psi;
    TV_ARRAYS grad_u,Normal;
    TV_ARRAYS velocity,root;
    int Extrapolate_Steps;

    T dt;
    T_ARRAYS u_temp,u_n,Iu_n_P,Eu_n_P,Iu_P,Eu_P;
    T_ARRAYS Iu_n_N,Eu_n_N,Iu_N,Eu_N;
    TV_ARRAYS grad_u_n,Egrad_u_n_P,Egrad_u_P;
    TV_ARRAYS Egrad_u_n_N,Egrad_u_N;

    int x_cells,y_cells,z_cells;

    Extrapolate(const T_GRID& grid,T_ARRAYS& u_, T_ARRAYS& phi_, TV_ARRAYS& psi_, int Extrapolate_Steps_);
    ~Extrapolate();

    inline int Heaviside(T heav){ return ((heav>(T)0.) ? 1 : 0); }

    void Extrapolate_Phi();
    void Compute_Gradient(TV_ARRAYS &psi,int order);
    void Calculate_Gradient(T_ARRAYS &psi_p,int component,int order);

    // Return variables
    T_ARRAYS Iu_P_F(){ return Iu_P; }
    T_ARRAYS Eu_P_F(){ return Eu_P; }
    TV_ARRAYS Egrad_u_P_F(){ return Egrad_u_P; }
    T_ARRAYS Iu_N_F(){ return Iu_N; }
    T_ARRAYS Eu_N_F(){ return Eu_N; }
    TV_ARRAYS Egrad_u_N_F(){ return Egrad_u_N; }
  };
}
