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

#include <memory>

#include "segment_creator.hpp"
#include "movement.hpp"
#include "segments.hpp"
#include "bounds.hpp"

MySegmentCreator::MySegmentCreator(Bounds& bounds, ISegments* segments, HeightMap& map_)
: bounds(bounds), segments(segments), map(map_) {

}


uint32_t MySegmentCreator::calcDirection(const Platec::vec2ui& point, const uint32_t origin_index, const uint32_t ID) const
{
    if(bounds.isInLimits(point))
    {
        return ID;
    }

    if (hasLowerID(getLeftIndex(origin_index),ID)) {
        return segments->id(getLeftIndex(origin_index));
    } else if (hasLowerID(getRightIndex(origin_index),ID)) {
        return segments->id(getRightIndex(origin_index));
    } else if (hasLowerID(getTopIndex(origin_index),ID)) {
       return segments->id(getTopIndex(origin_index));
    } else if (hasLowerID(getBottomIndex(origin_index),ID)) {
      return segments->id(getBottomIndex(origin_index));
    }

    return ID;
}

Span MySegmentCreator::scanSpans(std::vector<uint32_t>& spans_todo, std::vector<uint32_t>& spans_done) const
{
    Span span;
    Span tmpSpan;
    do // Find an unscanned span on this line.
    {
        span.end = spans_todo.back();
        spans_todo.pop_back();

        span.start = spans_todo.back();
        spans_todo.pop_back();

        // Reduce any done spans from this span.
        for (uint32_t j = 0; j < spans_done.size();
                j += 2)
        {
            // Saved coordinates are AT the point
            // that was included last to the span.
            // That's why equalities matter.
            tmpSpan = Span(spans_done[j],spans_done[j+1]);

            if (span.start >= tmpSpan.start &&
                    span.start <= tmpSpan.end)
                span.start = getRightIndex(tmpSpan.end);

            if (span.end >= tmpSpan.start &&
                    span.end <= tmpSpan.end)
                span.end = getLeftIndex(tmpSpan.start);
        }


        // Required to fix the underflow of end - 1.
        if(span.end >= bounds.width())
        {
            span.start = std::numeric_limits<uint32_t>::max();
            --span.end;
        }
    } while (span.notValid() && spans_todo.empty());
    return span;
}

ContinentId MySegmentCreator::createSegment(const Platec::vec2ui& point, 
                                    const Dimension& worldDimension) const
{
    const uint32_t bounds_width = bounds.width();
    const uint32_t bounds_height = bounds.height();
    const uint32_t origin_index = bounds.index(point);
    const uint32_t ID = segments->size();

    if (segments->id(origin_index) < ID) {
        return segments->id(origin_index);
    }

    uint32_t nbour_id = calcDirection(point, origin_index, ID);

    if (nbour_id < ID)
    {
        segments->setId(origin_index, nbour_id);
        (*segments)[nbour_id].incArea();

        (*segments)[nbour_id].enlarge_to_contain(point);

        return nbour_id;
    }

    uint32_t lines_processed;
    SegmentData pData = SegmentData(point,point, 0);
    static std::vector<uint32_t>* spans_todo = NULL;
    static std::vector<uint32_t>* spans_done = NULL;
    static uint32_t spans_size = 0;
    // MK: This code was originally allocating the 2D arrays per function call.
    // This was eating up a tremendous amount of cpu.
    // They are now static and they grow as needed, which turns out to be seldom.
    if (spans_size < bounds_height) {
        delete[] spans_todo;
        delete[] spans_done;
        spans_todo = new std::vector<uint32_t>[bounds_height];
        spans_done = new std::vector<uint32_t>[bounds_height];
        spans_size = bounds_height;
    }
    segments->setId(origin_index, ID);
    spans_todo[point.y()].push_back(point.x());
    spans_todo[point.y()].push_back(point.x());

    do
    {
        lines_processed = 0;
        for (uint32_t line = 0; line < bounds_height; ++line)
        {
            
            if (spans_todo[line].empty())
                continue;

            Span span = scanSpans(spans_todo[line], spans_done[line]);

            if (span.notValid()) // Nothing to do here anymore...
                continue;
                        
            // Calculate line indices. Allow wrapping around map edges.
            const uint32_t row_above = ((line - 1) & -(line > 0)) |
                                       ((bounds_height - 1) & -(line == 0));
            const uint32_t row_below = (line + 1) & -(line < bounds_height - 1);
            const uint32_t line_here = line * bounds_width;
            const uint32_t line_above = row_above * bounds_width;
            const uint32_t line_below = row_below * bounds_width;

            // Extend the beginning of line.
            while (span.start != 0 && segments->id(line_here+span.start-1) > ID &&
                    map[line_here+span.start-1] >= CONT_BASE)
            {
                --span.start;
                segments->setId(line_here + span.start, ID);

                // Count volume of pixel...
            }

            // Extend the end of line.
            while (span.end < bounds_width - 1 &&
                    segments->id(line_here + span.end + 1) > ID &&
                    map[line_here + span.end + 1] >= CONT_BASE)
            {
                ++span.end;
                segments->setId(line_here + span.end, ID);

                // Count volume of pixel...
            }

            // Check if should wrap around left edge.
            if (bounds_width == worldDimension.getWidth() && span.start == 0 &&
                    segments->id(line_here+bounds_width-1) > ID &&
                    map[line_here+bounds_width-1] >= CONT_BASE)
            {
                segments->setId(line_here + bounds_width - 1, ID);
                spans_todo[line].push_back(bounds_width - 1);
                spans_todo[line].push_back(bounds_width - 1);

                // Count volume of pixel...
            }

            // Check if should wrap around right edge.
            if (bounds_width == worldDimension.getWidth() && span.end == bounds_width - 1 &&
                    segments->id(line_here+0) > ID &&
                    map[line_here+0] >= CONT_BASE)
            {
                segments->setId(line_here + 0, ID);
                spans_todo[line].push_back(0);
                spans_todo[line].push_back(0);

                // Count volume of pixel...
            }

            pData.incArea(1 + span.end - span.start); // Update segment area counter.

            // Record any changes in extreme dimensions.
            if (line < pData.getTop()) pData.setTop(line);
            if (line > pData.getBottom()) pData.setBottom(line);
            if (span.start < pData.getLeft()) pData.setLeft(span.start);
            if (span.end > pData.getRight()) pData.setRight(span.end);

            if (line != 0 || bounds_height == worldDimension.getHeight()) {
                for (uint32_t j = span.start; j <= span.end; ++j)
                    if (segments->id(line_above + j) > ID &&
                            map[line_above + j] >= CONT_BASE)
                    {
                        uint32_t a = j;
                        segments->setId(line_above + a, ID);

                        // Count volume of pixel...

                        while (++j < bounds_width &&
                                segments->id(line_above + j) > ID &&
                                map[line_above + j] >= CONT_BASE)
                        {
                            segments->setId(line_above + j, ID);

                            // Count volume of pixel...
                        }

                        uint32_t b = --j; // Last point is invalid.

                        spans_todo[row_above].push_back(a);
                        spans_todo[row_above].push_back(b);
                        ++j; // Skip the last scanned point.
                    }
            }

            if (line < bounds_height - 1 || bounds_height == worldDimension.getHeight()) {
                for (uint32_t j = span.start; j <= span.end; ++j)
                    if (segments->id(line_below + j) > ID &&
                            map[line_below + j] >= CONT_BASE)
                    {
                        uint32_t a = j;
                        segments->setId(line_below + a, ID);

                        // Count volume of pixel...

                        while (++j < bounds_width &&
                                segments->id(line_below + j) > ID &&
                                map[line_below + j] >= CONT_BASE)
                        {
                            segments->setId(line_below + j, ID);

                            // Count volume of pixel...
                        }

                        uint32_t b = --j; // Last point is invalid.

                        spans_todo[row_below].push_back(a);
                        spans_todo[row_below].push_back(b);
                        ++j; // Skip the last scanned point.
                    }
            }

            spans_done[line].push_back(span.start);
            spans_done[line].push_back(span.end);
            ++lines_processed;
        }
    } while (lines_processed > 0);

    for (uint32_t line = 0; line < bounds_height; line++) {
        spans_todo[line].clear();
        spans_done[line].clear();
    }
    segments->add(pData);

    return ID;
}

const uint32_t MySegmentCreator::getBottomIndex(const int32_t originIndex) const {
    return originIndex + bounds.width();
}

const uint32_t MySegmentCreator::getLeftIndex(const int32_t originIndex) const {
    return originIndex-1;
}

const uint32_t MySegmentCreator::getRightIndex(const int32_t originIndex) const {
    return originIndex+1;
}

const uint32_t MySegmentCreator::getTopIndex(const int32_t originIndex) const {
    return originIndex - bounds.width();
}

const bool MySegmentCreator::hasLowerID(const uint32_t index, const ContinentId ID) const {
    //check if the value of the index is higher than CONT_BASE and
    //if ID is lower than the given ID
    return map.get(index) >= CONT_BASE && segments->id(index) < ID;
}
