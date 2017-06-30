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

#define _USE_MATH_DEFINES // Winblow$.


#include <cfloat>    // FT_EPSILON
#ifdef __MINGW32__ // this is to avoid a problem with the hypot function which is messed up by Python...
#undef __STRICT_ANSI__
#endif
#include <cmath>     // sin, cos
#include <cstdlib>   // rand
#include <vector>
#include <stdexcept> // std::invalid_argument
#include <assert.h>
#include <unistd.h>
#include "plate.hpp"
#include "heightmap.hpp"
#include "utils.hpp"
#include "plate_functions.hpp"

using namespace std;

plate::plate(const long seed,const HeightMap&  m, 
            const Dimension& plateDimension,
            const Platec::vec2f& topLeftCorner,
         const uint32_t plate_age,const Dimension& worldDimension) :
    randsource(seed),
    mass(MassBuilder(m).build()),
    map(m),
    age_map(std::vector<uint32_t>(plateDimension.getArea(),0),plateDimension),
    worldDimension(worldDimension),
    movement(randsource) {
    const uint32_t plate_area = plateDimension.getArea();

    bounds = std::make_shared<Bounds>(worldDimension,
                    topLeftCorner, plateDimension);

    
    auto mapItr = map.getData().begin();
    
    std::replace_if(age_map.getData().begin(),age_map.getData().end(),
                    [&](const auto& val)
                    {
                        auto ret = *mapItr > 0;
                        ++mapItr;
                        return ret;
                    }, plate_age);


    segments = std::make_shared<Segments>(plate_area);
    mySegmentCreator = std::make_shared<MySegmentCreator>(bounds, segments, map);
    segments->setSegmentCreator(mySegmentCreator);
    segments->setBounds(bounds);
}


uint32_t plate::addCollision(const Platec::vec2ui& point)
{
    auto& seg = getContinentAt(point);
    seg.incCollCount();
    return seg.getArea();
}

void plate::addCrustByCollision(const Platec::vec2ui& point,const float_t z,
                            const uint32_t time,const ContinentId activeContinent)
{
    // Add crust. Extend plate if necessary.
    setCrust(point, getCrust(point) + z, time);

    auto index = bounds->getValidMapIndex(point);
    segments->setId(index.first, activeContinent);

    ISegmentData& data = segments->getSegmentData(activeContinent);
    data.incArea();
    data.enlarge_to_contain(index.second);
}

void plate::addCrustBySubduction(const Platec::vec2ui& originPoint,const float_t sediment,const uint32_t time,
                              const Platec::vec2f& dir)
{
    if(sediment <= 0) //should be checked earliers
    {
        return;
    }
    // TODO: Create an array of coordinate changes that would create
    //       a circle around current point. Array is static and it is
    //       initialized at the first call to this function.
    //       After all points of the circle are checked around subduction
    //       point the one with most land mass around it will be chosen as
    //       "most in land" point and subducting crust is added there.
    //       However to achieve a little more "natural" look normal
    //       distributed randomness is added around the "center" point.
    //       Benefits:
    //           NEVER adds crust outside plate.
    //           ALWAYS goes inland as much as possible
    //       Drawbacks:
    //           Additional logic required
    //           Might place crust on other continent on same plate!


    // Take vector difference only between plates that move more or less
    // to same direction. This makes subduction direction behave better.
    //
    // Use of "this" pointer is not necessary, but it make code clearer.
    // Cursed be those who use "m_" prefix in member names! >(
    auto dotDir = dir;
    if(movement.dot(dir) > 0)
    {
        dotDir = dotDir - movement.velocityVector();
    }
    
    //What the hell is he doin here. Why have we to calculete 10 * x +3 ????
    auto offset =std::pow((float_t)randsource.next_double(),3);
     offset = std::copysign(offset,  2 * (int)(randsource.next() % 2) - 1); 
    
    auto offset2 =std::copysign(std::pow((float_t)randsource.next_double(),3),
                                                    randsource.next_signed()); 
    
    dotDir = Platec::vec2f(10 * dotDir.x() + 3 * offset,10 * dotDir.y() + 3 * offset2);

    const auto p =  Platec::vec2ui(dotDir.x(),dotDir.y()) + bounds->getValidMapIndex(originPoint).second; 
    if (bounds->isInLimits(p) )
    {
        
        uint32_t tmpindex = bounds->index(p);
        auto mapItr = map.getData().begin() + tmpindex;
        auto ageMapItr = age_map.getData().begin() + tmpindex;

        if (*mapItr > 0 )
        {

            *ageMapItr = ((*mapItr) * (*ageMapItr) + sediment * time) / ((*mapItr) + sediment);

            (*mapItr) += sediment;
            mass.incMass(sediment);
        }
    }
}

float plate::aggregateCrust(plate* p, uint32_t wx, uint32_t wy)
{

    const auto index = bounds->getValidMapIndex(Platec::vec2ui(wx,wy));

    const ContinentId seg_id = segments->id(index.first);

    // This check forces the caller to do things in proper order!
    //
    // Usually continents collide at several locations simultaneously.
    // Thus if this segment that is being merged now is removed from
    // segmentation bookkeeping, then the next point of collision that is
    // processed during the same iteration step would cause the test
    // below to be true and system would experience a premature abort.
    //
    // Therefore, segmentation bookkeeping is left intact. It doesn't
    // cause significant problems because all crust is cleared and empty
    // points are not processed at all. (Test on (seg_id >= seg_data.size()) removed)

    // One continent may have many points of collision. If one of them
    // causes continent to aggregate then all successive collisions and
    // attempts of aggregation would necessarily change nothing at all,
    // because the continent was removed from this plate earlier!
    if (segments->getSegmentData(seg_id).isEmpty()) {
        return 0;   // Do not process empty continents.
    }

    ContinentId activeContinent = p->selectCollisionSegment(Platec::vec2ui(wx, wy));

    // Wrap coordinates around world edges to safeguard subtractions.
    wx += worldDimension.getWidth();
    wy += worldDimension.getHeight();

    // Aggregating segment [%u, %u]x[%u, %u] vs. [%u, %u]@[%u, %u]\n",
    //      seg_data[seg_id].x0, seg_data[seg_id].y0,
    //      seg_data[seg_id].x1, seg_data[seg_id].y1,
    //      _dimension.getWidth(), _dimension.getHeight(), lx, ly);

    float old_mass = mass.getMass();

    // Add all of the collided continent's crust to destination plate.
    for (uint32_t y = segments->getSegmentData(seg_id).getTop(); y <= segments->getSegmentData(seg_id).getBottom(); ++y)
    {
        for (uint32_t x = segments->getSegmentData(seg_id).getLeft(); x <= segments->getSegmentData(seg_id).getRight(); ++x)
        {
            const uint32_t i = y * bounds->width() + x;
            if ((segments->id(i) == seg_id) && (map[i] > 0))
            {
                p->addCrustByCollision(Platec::vec2ui(wx + x - index.second.x(), wy + y - index.second.y()),
                                       map[i], age_map[i], activeContinent);

                mass.incMass(-1.0f * map[i]);
                map[i] = 0.0f;
            }
        }
    }

    segments->getSegmentData(seg_id).markNonExistent(); // Mark segment as non-existent
    return old_mass - mass.getMass();
}

void plate::applyFriction(const float_t deformed_mass)
{
    // Remove the energy that deformation consumed from plate's kinetic
    // energy: F - dF = ma - dF => a = dF/m.
    if (!mass.isNull())
    {
        movement.applyFriction(deformed_mass, mass.getMass());
    }
}

void plate::collide(plate& p,const float_t coll_mass)
{
    if (!mass.isNull() && coll_mass > 0) {
        movement.collide(mass, p, coll_mass);
    }
}

const surroundingPoints plate::calculateCrust(const Platec::vec2ui& position, const float_t height)
{
    surroundingPoints ret;
    // Build masks for accessible directions (4-way).
    // Allow wrapping around map edges if plate has world wide dimensions.
    const uint32_t world_height = worldDimension.getHeight();
    bool height_bit = bounds->height() == world_height;
    
    if(position.x() != 0)
    {
        ret.westIndex = bounds->index(position- Platec::vec2ui(1,0));
        if(map[ret.westIndex] < height)
        {
          ret.westCrust = map[ret.westIndex];
        }
    }
    else
    {
        ret.westIndex = bounds->index(Platec::vec2ui(0,position.y()));
    }
    
  
    ret.eastIndex = bounds->index(bounds->getDimension().xMod(position +Platec::vec2ui(1,0))) ;
    if(map[ret.eastIndex] < height)
    {
      ret.eastCrust = map[ret.eastIndex];
    }

//    if(position.y() > 0| height_bit)
//    {
////        ret.northIndex = bounds->index(position - Platec::vec2ui(0,1));
////        if(map[ret.northIndex] < height)
////        {
////          ret.northCrust = map[ret.northIndex];
////        }
//
//        // Calculate the x and y offset of neighbour directions.
//        // If neighbour is out of plate edges, set it to zero. This protects
//        // map memory reads from segment faulting.
//        const uint32_t y_mod = bounds->getDimension().yMod(position.y());
//
//        uint32_t y_mod_minus_1 = y_mod == 0 ? 0 : y_mod - 1;
//        
//        ret.northIndex = bounds->index(Platec::vec2ui(position.x(),y_mod_minus_1));
//
//        ret.northCrust = map[ret.northIndex] * ((map[ret.northIndex] < height));
//    }
//    else
//    {
//        ret.northCrust = bounds->index(Platec::vec2ui(position.x(),0));
//    }
    
    
   ret.southIndex = bounds->index(bounds->getDimension().yMod(position +Platec::vec2ui(0,1))) ;
    if(map[ret.southIndex] < height)
    {
      ret.southCrust = map[ret.southIndex];
    }    
    
    uint32_t n_mask = -((position.y() > 0)   | height_bit       );

    // Calculate the x and y offset of neighbour directions.
    // If neighbour is out of plate edges, set it to zero. This protects
    // map memory reads from segment faulting.
    const uint32_t y_mod = bounds->getDimension().yMod(position.y());
            
    uint32_t y_mod_minus_1 = y_mod == 0 ? 0 : y_mod - 1;

     uint32_t       n = n_mask==-1 ? y_mod_minus_1 : 0;
    ret.northIndex = bounds->index(Platec::vec2ui(position.x(),n));

    ret.northCrust = map[ret.northIndex] * (n_mask & (map[ret.northIndex] < height));

    return ret;
}

void plate::findRiverSources(const float_t lower_bound, vector<uint32_t>& sources)
{
    const uint32_t bounds_height = bounds->height();
    const uint32_t bounds_width = bounds->width();
    

    // Find all tops.
    for (uint32_t y = 0; y < bounds_height; ++y) {
      //  const uint32_t y_width = y * bounds_width;
        for (uint32_t x = 0; x < bounds_width; ++x) {
            const uint32_t index = bounds->index(Platec::vec2ui(x, y));

            if (map[index] < lower_bound) {
                continue;
            }


            // This location is either at the edge of the plate or it is not the
            // tallest of its neightbours. Don't start a river from here.
            if (calculateCrust(Platec::vec2ui(x, y),map[index]).onIsLower()) {
                continue;
            }

            sources.push_back(index);
        }
    }
}

void plate::flowRivers(float lower_bound, vector<uint32_t>* sources, HeightMap& tmp)
{
    const uint32_t bounds_area = bounds->area();
    vector<uint32_t> sinks_data;

    static vector<bool> s_flowDone;
    if (s_flowDone.size() < bounds_area) {
        s_flowDone.resize(bounds_area);
    }
    fill(s_flowDone.begin(), s_flowDone.begin() + bounds_area, false);

    // From each top, start flowing water along the steepest slope.
    while (!sources->empty()) {
        while (!sources->empty()) {
            const uint32_t index = sources->back();
            const uint32_t y = index / bounds->width();
            const uint32_t x = index - y * bounds->width();

            sources->pop_back();

            if (map[index] < lower_bound) {
                continue;
            }

            surroundingPoints neighbors = calculateCrust(Platec::vec2ui(x, y), map[index]);

            // If this is the lowest part of its neighbourhood, stop.
            if (neighbors.oneIsHigher()) {
                continue;
            }

            neighbors.westCrust += (neighbors.westCrust == 0) * map[index];
            neighbors.eastCrust += (neighbors.eastCrust == 0) * map[index];
            neighbors.northCrust += (neighbors.northCrust == 0) * map[index];
            neighbors.southCrust += (neighbors.southCrust == 0) * map[index];

            // Find lowest neighbour.
            float lowest_crust = neighbors.westCrust;
            uint32_t dest = index - 1;

            if (neighbors.eastCrust < lowest_crust) {
                lowest_crust = neighbors.eastCrust;
                dest = index + 1;
            }

            if (neighbors.northCrust < lowest_crust) {
                lowest_crust = neighbors.northCrust;
                dest = index - bounds->width();
            }

            if (neighbors.southCrust < lowest_crust) {
                lowest_crust = neighbors.southCrust;
                dest = index + bounds->width();
            }

            // if it's not handled yet, add it as new sink.
            if (dest < bounds->area() && !s_flowDone[dest]) {
                sinks_data.push_back(dest);
                s_flowDone[dest] = true;
            }

            // Erode this location with the water flow.
            tmp[index] -= (tmp[index] - lower_bound) * 0.2;
        }


        sources->swap(sinks_data);
        sinks_data.clear();
    }
}

void plate::erode(float lower_bound)
{
    vector<uint32_t> sources_data;
    vector<uint32_t>* sources = &sources_data;

    HeightMap tmpHm(map);
    findRiverSources(lower_bound, sources_data);
    flowRivers(lower_bound, sources, tmpHm);

    // Add random noise (10 %) to heightmap.
    for (uint32_t i = 0; i < bounds->area(); ++i) {
        float alpha = 0.2 * (float)randsource.next_double();
        tmpHm[i] += 0.1 * tmpHm[i] - alpha * tmpHm[i];
    }

    map = tmpHm;
    tmpHm.set_all(0.0f);
    MassBuilder massBuilder;

    for (uint32_t y = 0; y < bounds->height(); ++y)
    {
        for (uint32_t x = 0; x < bounds->width(); ++x)
        {
            const uint32_t index = y * bounds->width() + x;
            massBuilder.addPoint(Platec::Vector2D<uint32_t>(x, y), map[index]);
            tmpHm[index] += map[index]; // Careful not to overwrite earlier amounts.

            if (map[index] < lower_bound)
                continue;

            surroundingPoints neighbors = calculateCrust(Platec::vec2ui(x, y), map[index]);
            
            // This location has no neighbours (ARTIFACT!) or it is the lowest
            // part of its area. In either case the work here is done.
            if (neighbors.oneIsHigher())
                continue;

            // The steeper the slope, the more water flows along it.
            // The more downhill (sources), the more water flows to here.
            // 1+1+10 = 12, avg = 4, stdev = sqrt((3*3+3*3+6*6)/3) = 4.2, var = 18,
            //  1*1+1*1+10*10 = 102, 102/4.2=24
            // 1+4+7 = 12, avg = 4, stdev = sqrt((3*3+0*0+3*3)/3) = 2.4, var = 6,
            //  1*1+4*4+7*7 = 66, 66/2.4 = 27
            // 4+4+4 = 12, avg = 4, stdev = sqrt((0*0+0*0+0*0)/3) = 0, var = 0,
            //  4*4+4*4+4*4 = 48, 48/0 = inf -> 48
            // If there's a source slope of height X then it will always cause
            // water erosion of amount Y. Then again from one spot only so much
            // water can flow.
            // Thus, the calculated non-linear flow value for this location is
            // multiplied by the "water erosion" constant.
            // The result is max(result, 1.0). New height of this location could
            // be e.g. h_lowest + (1 - 1 / result) * (h_0 - h_lowest).

            // Calculate the difference in height between this point and its
            // nbours that are lower than this point.
            float w_diff = map[index] - neighbors.westCrust;
            float e_diff = map[index] - neighbors.eastCrust;
            float n_diff = map[index] - neighbors.northCrust;
            float s_diff = map[index] - neighbors.southCrust;

            float min_diff = w_diff;
            min_diff -= (min_diff - e_diff) * (e_diff < min_diff);
            min_diff -= (min_diff - n_diff) * (n_diff < min_diff);
            min_diff -= (min_diff - s_diff) * (s_diff < min_diff);

            // Calculate the sum of difference between lower neighbours and
            // the TALLEST lower neighbour.
            float diff_sum = (w_diff - min_diff) * (neighbors.westCrust > 0) +
                             (e_diff - min_diff) * (neighbors.eastCrust > 0) +
                             (n_diff - min_diff) * (neighbors.northCrust > 0) +
                             (s_diff - min_diff) * (neighbors.southCrust > 0);

            // Erosion difference sum is negative!
            ASSERT(diff_sum >= 0, "Difference sum must be positive");

            if (diff_sum < min_diff)
            {
                // There's NOT enough room in neighbours to contain all the
                // crust from this peak so that it would be as tall as its
                // tallest lower neighbour. Thus first step is make ALL
                // lower neighbours and this point equally tall.
                tmpHm[neighbors.westIndex ] += (w_diff - min_diff) * (neighbors.westCrust > 0);
                tmpHm[neighbors.eastIndex ] += (e_diff - min_diff) * (neighbors.eastCrust > 0);
                tmpHm[neighbors.northIndex] += (n_diff - min_diff) * (neighbors.northCrust > 0);
                tmpHm[neighbors.southIndex] += (s_diff - min_diff) * (neighbors.southCrust > 0);
                tmpHm[index] -= min_diff;

                min_diff -= diff_sum;

                // Spread the remaining crust equally among all lower nbours.
                min_diff /= 1 + (neighbors.westCrust > 0) + (neighbors.eastCrust > 0) +
                            (neighbors.northCrust > 0) + (neighbors.southCrust > 0);

                tmpHm[neighbors.westIndex ] += min_diff * (neighbors.westCrust > 0);
                tmpHm[neighbors.eastIndex ] += min_diff * (neighbors.eastCrust > 0);
                tmpHm[neighbors.northIndex] += min_diff * (neighbors.northCrust > 0);
                tmpHm[neighbors.southIndex] += min_diff * (neighbors.southCrust > 0);
                tmpHm[index] += min_diff;
            }
            else
            {
                float unit = min_diff / diff_sum;

                // Remove all crust from this location making it as tall as
                // its tallest lower neighbour.
                tmpHm[index] -= min_diff;

                // Spread all removed crust among all other lower neighbours.
                tmpHm[neighbors.westIndex ] += unit * (w_diff - min_diff) * (neighbors.westCrust > 0);
                tmpHm[neighbors.eastIndex ] += unit * (e_diff - min_diff) * (neighbors.eastCrust > 0);
                tmpHm[neighbors.northIndex] += unit * (n_diff - min_diff) * (neighbors.northCrust > 0);
                tmpHm[neighbors.southIndex] += unit * (s_diff - min_diff) * (neighbors.southCrust > 0);
            }
        }
    }
    map = tmpHm;
    mass = massBuilder.build();
}

const std::pair<uint32_t,float_t> plate::getCollisionInfo(const Platec::vec2ui& point) const
{
    const ISegmentData& seg = getContinentAt(point);

    if(seg.getArea() == 0)
    {
        return  std::make_pair(seg.collCount(),seg.collCount());
    }
    return std::make_pair(seg.collCount(),seg.collCount() /(float_t)(seg.getArea()));


}

uint32_t plate::getContinentArea(const Platec::vec2ui& point) const
{
    const auto index = bounds->getValidMapIndex(point);
    ASSERT(segments->id(index.first) < segments->size(), "Segment index invalid");
    return segments->getSegmentData(segments->id(index.first)).getArea();
}

float_t plate::getCrust(const Platec::vec2ui& point) const
{
    const auto index = bounds->getMapIndex(point);
    return index.first != BAD_INDEX ? map[index.first] : 0;
}

uint32_t plate::getCrustTimestamp(const Platec::vec2ui& point) const
{
    const auto index = bounds->getMapIndex(point);
    return index.first != BAD_INDEX ? age_map[index.first] : 0;
}

void plate::getMap(const float** c, const uint32_t** t) const
{
    if (c) {
        *c = map.raw_data();
    }
    if (t) {
        *t = age_map.raw_data();
    }
}

void plate::move(const Dimension& worldDimension)
{
    movement.move(worldDimension);
    bounds->shift(movement.velocityVector());
}

void plate::resetSegments()
{
    ASSERT(bounds->area() == segments->getArea(), "Segments doesn't have the expected area");
    segments->reset();
}

void plate::setCrust(const Platec::vec2ui& point, float_t z, uint32_t t)
{
    
    // Do not accept negative values.
    z = std::max(0.f,z);


    auto index = bounds->getMapIndex(point);

    if (index.first == BAD_INDEX)
    {
        // Extending plate for nothing!
        ASSERT(z > 0, "Height value must be non-zero");

        const uint32_t ilft = bounds->left();
        const uint32_t itop = bounds->top();
        const uint32_t irgt = bounds->rightAsUintNonInclusive();
        const uint32_t ibtm = bounds->bottomAsUintNonInclusive();

        auto normPoint = worldDimension.normalize(point);
        

        
        // Calculate distance of new point from plate edges.
        const uint32_t _lft = ilft - normPoint.x();
        const uint32_t _rgt = (worldDimension.getWidth() & -(normPoint.x() < ilft)) + normPoint.x() - irgt;
        const uint32_t _top = itop - normPoint.y();
        const uint32_t _btm = (worldDimension.getHeight() & -(normPoint.y() < itop)) + normPoint.y() - ibtm;

        // Set larger of horizontal/vertical distance to zero.
        // A valid distance is NEVER larger than world's side's length!
        uint32_t d_lft = _lft & -(_lft <  _rgt) & -(_lft < worldDimension.getWidth());
        uint32_t d_rgt = _rgt & -(_rgt <= _lft) & -(_rgt < worldDimension.getWidth());
        uint32_t d_top = _top & -(_top <  _btm) & -(_top < worldDimension.getHeight());
        uint32_t d_btm = _btm & -(_btm <= _top) & -(_btm < worldDimension.getHeight());

        // Scale all changes to multiple of 8.
        d_lft = ((d_lft > 0) + (d_lft >> 3)) << 3;
        d_rgt = ((d_rgt > 0) + (d_rgt >> 3)) << 3;
        d_top = ((d_top > 0) + (d_top >> 3)) << 3;
        d_btm = ((d_btm > 0) + (d_btm >> 3)) << 3;

        // Make sure plate doesn't grow bigger than the system it's in!
        if (bounds->width() + d_lft + d_rgt > worldDimension.getWidth())
        {
            d_lft = 0;
            d_rgt = worldDimension.getWidth() - bounds->width();
        }

        if (bounds->height() + d_top + d_btm > worldDimension.getHeight())
        {
            d_top = 0;
            d_btm = worldDimension.getHeight() - bounds->height();
        }

        // Index out of bounds, but nowhere to grow!
        ASSERT(d_lft + d_rgt + d_top + d_btm != 0, "Invalid plate growth deltas");

        const uint32_t old_width  = bounds->width();
        const uint32_t old_height = bounds->height();

        bounds->shift(Platec::vec2f(-1.0*d_lft, -1.0*d_top));
        bounds->grow(Platec::vec2ui(d_lft + d_rgt, d_top + d_btm));

        HeightMap tmph = HeightMap(bounds->width(), bounds->height());
        AgeMap    tmpa = AgeMap(bounds->width(), bounds->height());
        uint32_t* tmps = new uint32_t[bounds->area()];
        tmph.set_all(0);
        tmpa.set_all(0);
        memset(tmps, 255, bounds->area()*sizeof(uint32_t));

        // copy old plate into new.
        for (uint32_t j = 0; j < old_height; ++j)
        {
            const uint32_t dest_i = (d_top + j) * bounds->width() + d_lft;
            const uint32_t src_i = j * old_width;
            memcpy(&tmph[dest_i], &map[src_i], old_width *
                   sizeof(float));
            memcpy(&tmpa[dest_i], &age_map[src_i], old_width *
                   sizeof(uint32_t));
            memcpy(&tmps[dest_i], &segments->id(src_i), old_width *
                   sizeof(uint32_t));
        }

        map     = tmph;
        age_map = tmpa;
        segments->reassign(bounds->area(),std::vector<uint32_t>(tmps, tmps +bounds->area()));

        // Shift all segment data to match new coordinates.
        segments->shift(Platec::vec2ui(d_lft, d_top));


        index = bounds->getValidMapIndex(normPoint);

        assert(index.first < bounds->area());
    }

    // Update crust's age.
    // If old crust exists, new age is mean of original and supplied ages.
    // If no new crust is added, original time remains intact.
    const uint32_t old_crust = -(map[index.first] > 0);
    const uint32_t new_crust = -(z > 0);
    t = (t & ~old_crust) | ((uint32_t)((map[index.first] * age_map[index.first] + z * t) /
                                       (map[index.first] + z)) & old_crust);
    age_map[index.first] = (t & new_crust) | (age_map[index.first] & ~new_crust);

    mass.incMass(-1.0f * map[index.first]);
    mass.incMass(z);      // Update mass counter.
    map[index.first] = z;     // Set new crust height to desired location.
}

ContinentId plate::selectCollisionSegment(const Platec::vec2ui& point) const
{
    return segments->id(bounds->getValidMapIndex(point).first);
}

///////////////////////////////////////////////////////////////////////////////
/// Private methods ///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

uint32_t plate::createSegment(const Platec::vec2ui& point)
{
    return mySegmentCreator->createSegment(point,worldDimension);
}

ISegmentData& plate::getContinentAt(const Platec::vec2ui& point)
{
    return segments->getSegmentData(segments->getContinentAt(point,worldDimension));
}

const ISegmentData& plate::getContinentAt(const Platec::vec2ui& point) const
{
    return segments->getSegmentData(segments->getContinentAt(point,worldDimension));
}
