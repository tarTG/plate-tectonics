/******************************************************************************
 *  plate-tectonics, a plate tectonics simulation library
 *  Copyright (C) 2012-2013 Lauri Viitanen
 *  Copyright (C) 2014-2015 Federico Tomassetti, Bret Curtis
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, see http://www.gnu.org/licenses/
 *****************************************************************************/

#include "mass.hpp"

#include <algorithm>

// ----------------------------------------------
// MassBuilder
// ----------------------------------------------

MassBuilder::MassBuilder(const HeightMap& map)
    : mass(0), center(0.f,0.f)
{
    uint32_t index = 0;
    for(const auto& data : map.getData())
    {
        addPoint(map.getDimension().coordOF(index),data);
        ++index;
    }

}

MassBuilder::MassBuilder()
    : mass(0), center(0.f,0.f)
{
}

void MassBuilder::addPoint(const Platec::Vector2D<uint32_t>& point,
                            const float crust)
{
    ASSERT(crust >= 0.0f, "Crust should be not negative");
    mass += crust;
    // Update the center coordinates weighted by mass.
    center.shift(Platec::Vector2D<float_t>(point.x(),point.y()) * crust);
}

Mass MassBuilder::build()
{
    if (mass <= 0) {
        return Mass(0, 0, 0);
    } else {
        ASSERT(mass > 0, "Mass was zero!");
        float inv_mass = 1 / mass;
        return Mass(mass, center.x() * inv_mass, center.y() * inv_mass);
    }
}

// ----------------------------------------------
// Mass
// ----------------------------------------------

Mass::Mass(float mass_, float cx_, float cy_)
    : mass(mass_), cx(cx_), cy(cy_), _totalX(0), _totalY(0)
{

}

void Mass::incMass(float delta)
{
    mass += delta;
    if (mass < 0.0f) {
        if (mass > -0.01f) {
            mass = 0.0f;
        } else {
            ASSERT(0, "A negative mass is not allowed");
        }
    }
}

float Mass::getMass() const
{
    return mass;
}

float Mass::getCx() const
{
    return cx;
}

float Mass::getCy() const
{
    return cy;
}

bool Mass::null() const
{
    return mass <= 0;
}

