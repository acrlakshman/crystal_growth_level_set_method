//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#pragma once

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include "RANGE_ITERATOR.h"

namespace PhysBAM{

  template<class T>
class HermiteInterpolation
  {
  public:
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef GRID<TV> T_GRID;
    typedef ARRAY<T,T_INDEX> T_ARRAYS;
    typedef ARRAY<TV,T_INDEX> TV_ARRAYS;

    const T_GRID& grid;
    TV_ARRAYS& root;
    T_ARRAYS& phi;
    TV_ARRAYS& psi;
    T_ARRAYS Iphi;
    TV_ARRAYS Ipsi;
    T_INDEX reference_index;

    int counterg;
    int x_cells; // number of cells in x-direction
    int y_cells; // number of cells in y-direction
    int z_cells; // number of cells in z-direction

    // Cross derivatives to form interpolating polynomial
    T_ARRAYS psi_xy;
    T_ARRAYS psi_xz;
    T_ARRAYS psi_yz;
    T_ARRAYS psi_xyz;

    HermiteInterpolation(const T_GRID& grid, TV_ARRAYS& root, T_ARRAYS& phi_, TV_ARRAYS& psi_);
    ~HermiteInterpolation();

    void Interpolate_Phi();
    void Compute_Cross_Derivative(int component); // component=1(psi_xy); 2(psi_xz); 3(psi_yz); 4(psi_xyz)
    void Calculate_Iphi(T_INDEX index,T_INDEX reference_index);
    void Calculate_Ipsi(T_INDEX index,T_INDEX reference_index);
    T Calculate_Interpolant(T_INDEX index,T_INDEX reference_index,int derivative,T hC);
    T BasisPoly(int v1,int v2,int v3,int d1,int d2,int d3,int derivative,T hC,T etax,T etay,T etaz);
    T CubicPoly(T eta,int func);

    T_ARRAYS Iphi_F(){ return Iphi; }
    TV_ARRAYS Ipsi_F(){ return Ipsi; }
    T_ARRAYS psi_xy_F(){ return psi_xy; }
    T_ARRAYS psi_xz_F(){ return psi_xz; }
    T_ARRAYS psi_yz_F(){ return psi_yz; }
    T_ARRAYS psi_xyz_F(){ return psi_xyz; }
  };
}
