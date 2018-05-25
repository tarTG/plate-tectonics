/******************************************************************************
 *  plate-tectonics, a plate tectonics simulation library
 *  Copyright (C) 2010 Craig McQueen (http://craig.mcqueen.id.au)
 *  Copyright (C) 2014-2015 Federico Tomassetti, Bret Curtis
 *
 *  This is code from the Simple Pseudo-random Number Generators
 *  Available on GitHub https://github.com/cmcqueen/simplerandom
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

#define NOMINMAX


#include "simplerandom.hpp"
#include <stddef.h>
#include "utils.hpp"



SimpleRandom::SimpleRandom(uint32_t seed)
{
    rng = std::mt19937(seed);
    uintdist =  std::uniform_int_distribution<uint32_t>(std::numeric_limits<uint32_t>::min(),std::numeric_limits<uint32_t>::max() );
    intdist = std::uniform_int_distribution<int32_t>(std::numeric_limits<int32_t>::min(),std::numeric_limits<int32_t>::max());
    doubledist = std::uniform_real_distribution<float>(0.0f, 1.0f);
    floatdist = std::uniform_real_distribution<float>(-0.5f, 0.5f);    
}



uint32_t SimpleRandom::next()
{
    return uintdist(rng);
}

float SimpleRandom::next_float()
{
    return doubledist(rng);
}

float SimpleRandom::next_float_signed()
{
    return floatdist(rng);
}

int32_t SimpleRandom::next_signed()
{
    return intdist(rng);
}