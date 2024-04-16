//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

#include "geometries/geometry.h"

#include "voxel_mesher_coloring.h"

namespace Kratos {

class ColorCellFacesBetweenColors: public VoxelMesherColoring {

public:
    ColorCellFacesBetweenColors(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
        VoxelMesherColoring(rModeler, ColoringParameters)
    {}

    void Apply() const override
    {
        Parameters parameters = GetParameters();
        const int outside_color = parameters["outside_color"].GetInt();
        const int interface_color = parameters["color"].GetInt();
        const int cell_color = parameters["cell_color"].GetInt();

        GetMeshColors().CalculateElementalFaceColorsBetweenColors(interface_color, outside_color, cell_color);
    }

};

}