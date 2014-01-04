//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OBJ_FILE
//#####################################################################
#ifndef __OBJ_FILE__
#define __OBJ_FILE__

#include <fstream>

namespace PhysBAM{

template<class TV,class T_MESH> void
Write_Obj_File(const MESH_OBJECT<TV,T_MESH>& mesh_object,const std::string& filename)
{
    std::ofstream output(filename.c_str());
    if(!output.is_open())
        PHYSBAM_FATAL_ERROR("Could not open file "+filename+" for writing");
    static const int d1=TV::dimension;
    static const int d2=T_MESH::dimension+1;

    for(int p=1;p<=mesh_object.particles.array_collection->Size();p++){
        output<<"v";
        for(int i=1;i<=d1;i++)
            output<<" "<<mesh_object.particles.X(p)(i);
	if(d1==2) output<<" 0.0";
        output<<std::endl;}

    for(int e=1;e<=mesh_object.mesh.elements.m;e++){
        output<<"f";
        for(int i=1;i<=d2;i++)
            output<<" "<<mesh_object.mesh.elements(e)(i);
        output<<std::endl;}

    output.close();
    
}

}
#endif
