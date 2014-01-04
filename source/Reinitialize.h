//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#pragma once

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include "RANGE_ITERATOR.h"
#include "HermiteInterpolation.h"

namespace PhysBAM{

  template<class T>
class Reinitialize
  {
  public:
    static const int d=3;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef GRID<TV> T_GRID;
    typedef ARRAY<T,T_INDEX> T_ARRAYS;
    typedef ARRAY<TV,T_INDEX> TV_ARRAYS;

    const T_GRID& grid;
    TV_ARRAYS root;
    const T_ARRAYS& phi;
    const TV_ARRAYS& psi;
    TV_ARRAYS velocity;
    int Reinit_Steps;

    T dt;
    TV_ARRAYS root_M;
    T_ARRAYS Iphi,Rphi;
    TV_ARRAYS Rpsi;
    T_ARRAYS Interfacial_Rectangles,Convergence_Flag; // both are defined as floats

    int N_Rectangles; // Number of nodes that are tagged
    int x_cells,y_cells,z_cells;
    T criterion;

    // Cross derivatives to form interpolating polynomial
    T_ARRAYS psi_xy;
    T_ARRAYS psi_xz;
    T_ARRAYS psi_yz;
    T_ARRAYS psi_xyz;

    Reinitialize(const T_GRID& grid, T_ARRAYS& phi_, TV_ARRAYS& psi_, int Reinit_Steps_);
    ~Reinitialize();

    void Reinitialize_Phi();
    void Locate_Interface_Modified_Newtons_Method(HermiteInterpolation<T>& interpolant);
    void Locate_Interfacial_Rectangles(const T_ARRAYS& phi,T_ARRAYS& Iphi);
    void Modified_Newtons_Method(HermiteInterpolation<T>& interpolant,T_INDEX index);

    // Return variables
    T_ARRAYS Iphi_F(){ return Iphi; }
    T_ARRAYS Rphi_F(){ return Rphi; }
    TV_ARRAYS Rpsi_F(){ return Rpsi; }
  };
}
