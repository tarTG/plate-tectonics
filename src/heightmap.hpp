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

#ifndef HEIGHTMAP_HPP
#define HEIGHTMAP_HPP


#include <stdint.h>
#include <vector>

#include "dimension.h"


template <typename Value>
class Matrix
{
private:

    std::vector<Value> data;
    Dimension dimension;   
    
public:
    
    Matrix(const uint32_t width,const uint32_t height)
        : data(width*height,0.0), dimension(width, height)  {
    }

    Matrix(const std::vector<Value>& dataVal,
            const uint32_t width,const uint32_t height)
        : data(dataVal), dimension(width, height) {
    }
    
    Matrix(const Dimension& dim) : data(dim.getArea(),0.0), dimension(dim) {
    }
    
    Matrix(const std::vector<Value>& dataVal,const Dimension& dim)
        : data(dataVal), dimension(dim) {
    }    

    void set_all(const Value value) {
   //     std::vector<Value>(data.size(),value).swap(data);
        for(auto& d : data)
        {
            d = value;
        }
    }
    

    void set(const Platec::vec2ui& point, const Value value) {
        data.at(dimension.indexOf(point)) = value;
    }

    const Value get(const Platec::vec2ui& point) const {
        return data.at(dimension.indexOf(point));
    }
    
    const Value get(const uint32_t index) const {
        return data.at(index);
    }
    
    Value& operator[](uint32_t index) {
        return data.at(index);
    }
    const Value& operator[](uint32_t index) const {
        return data.at(index);;
    }        

    const Value& operator[](const Platec::vec2ui& point) const {
        return data.at(dimension.indexOf(point));;
    }
    
    Value& operator[](const Platec::vec2ui& point) {
        return data.at(dimension.indexOf(point));
    }

    const Value* raw_data() const {
        // cant use .data() because msvc
        return &data[0];
    }
    
    Value* raw_data() {
        // cant use .data() because msvc
        return &data[0];
    }
        
    uint32_t width() const {
        return dimension.getWidth();
    }
    
    uint32_t height() const {
       return dimension.getHeight();
    }
    uint32_t area() const {
        return dimension.getArea(); 
    }

    const Dimension& getDimension() const {
        return dimension;
    }

    const std::vector<Value>& getData() const {
        return data;
    }
    std::vector<Value>& getData() {
        return data;
    }
};

typedef Matrix<float_t> HeightMap;
typedef Matrix<uint32_t> AgeMap;
typedef Matrix<uint32_t> IndexMap;

#endif
