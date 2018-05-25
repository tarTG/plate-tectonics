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

#ifndef SIMPLE_RANDOM_HPP
#define SIMPLE_RANDOM_HPP

#include <random>



class SimpleRandom {
public:
    SimpleRandom(uint32_t seed);
    uint32_t next();
    int32_t next_signed();
    // Return a random value in [0.0, 1.0]
    float next_float();
    // Return a random value in [-0.5f, 0.5f]
    float next_float_signed();
private:
    std::mt19937 rng;
    std::uniform_int_distribution<uint32_t> uintdist;
    std::uniform_int_distribution<int32_t> intdist;
    std::uniform_real_distribution<float> doubledist;
    std::uniform_real_distribution<float> floatdist;
};

#endif
