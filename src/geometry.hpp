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

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#define NOMINMAX

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <type_traits>
#include "utils.hpp"
#include "Vector2D.h"


class WorldDimension;

/// Dimension of a Rectangle.
class Dimension {
    
protected:
    Platec::vec2ui dim;
public:

    /// Initialize the dimension with the given values
    Dimension(const uint32_t width,const uint32_t height);

    uint32_t getWidth() const {
        return dim.x();
    }
    uint32_t getHeight() const {
        return dim.y();
    }
    uint32_t getArea() const {
        return getWidth() * getHeight();
    }
    uint32_t indexOf(const uint32_t x, const uint32_t y) const;
    uint32_t indexOf(const Platec::vec2ui& point) const;
    
    uint32_t yFromIndex(const uint32_t index) const;
    uint32_t xFromIndex(const uint32_t index) const;    
    
    const Platec::vec2ui coordOF(const uint32_t index) const;

    template <class T>
    bool contains(const Platec::Vector2D<T>& p) const {
        //using std::floor here to avoid floating point inaccuarcy
       return (std::floor(p.x()) >= 0 && std::floor(p.x()) < getWidth()  
               && std::floor(p.y()) >= 0.0f && std::floor(p.y()) < getHeight());
    }
    void grow(Platec::vec2ui growSize);

};

class WorldDimension : public Dimension {
public:
    WorldDimension(const uint32_t width,const uint32_t height);
    uint32_t getMax() const;
    uint32_t xMod(const uint32_t x) const;
    Platec::vec2ui xMod(const Platec::vec2ui& point) const;
    uint32_t yMod(const uint32_t y) const;
    Platec::vec2ui yMod(const Platec::vec2ui& point) const;
    Platec::vec2ui pointMod(const Platec::vec2ui& point) const;
    Platec::vec2ui  normalize(const Platec::vec2ui& point) const;
    uint32_t lineIndex(const uint32_t y) const;
    uint32_t normalizedIndexOf(const uint32_t x, const uint32_t y) const;
    uint32_t normalizedIndexOf(const Platec::vec2ui& point) const; 
    uint32_t xCap(const uint32_t x) const;
    Platec::vec2ui xCap(const Platec::vec2ui& point) const;
    uint32_t yCap(const uint32_t y) const;
    Platec::vec2ui yCap(const Platec::vec2ui& point) const;
    
    template <class T>
    Platec::Vector2D<T> wrap(const Platec::Vector2D<T>& point) const
    {
        T xval = point.x(), yval = point.y();
        if(std::floor(xval) < 0)
        {
            xval += static_cast<T>(getWidth());
        }
        else if (std::floor(xval) > getWidth())
        {
            xval -= static_cast<T>(getWidth());
        }
        
        if(std::floor(yval) < 0)
        {
            yval += static_cast<T>(getHeight());
        }
        else if (std::floor(yval) > getHeight())
        {
            yval -= static_cast<T>(getHeight());
        }
        return Platec::Vector2D<T> (xval,yval);
    }

};

#endif
