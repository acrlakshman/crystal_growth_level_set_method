//#####################################################################
// Copyright 2012, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

// set delta T

dt=min(grid.dX.x,grid.dX.y,grid.dX.z);
dt=(T).5*dt*dt;
cout<<"dt="<<dt<<endl;
