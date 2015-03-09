#ifndef SEGMENTS_HPP
#define SEGMENTS_HPP

#include <vector>
#include <cmath>     // sin, cos
#include "simplerandom.hpp"
#include "heightmap.hpp"
#include "rectangle.hpp"
#include "segment_data.hpp"
#include "utils.hpp"
#include "bounds.hpp"
#include "movement.hpp"
#include "mass.hpp"
#include "segment_creator.hpp"

typedef uint32_t ContinentId;

class ISegments
{
public:
    virtual uint32_t area() = 0;
    virtual void reset() = 0;
    virtual void reassign(uint32_t newarea, uint32_t* tmps) = 0;
    virtual void shift(uint32_t d_lft, uint32_t d_top) = 0;
    virtual uint32_t size() const = 0;
    virtual const SegmentData& operator[](uint32_t index) const = 0;
    virtual SegmentData& operator[](uint32_t index) = 0;
    virtual void add(const SegmentData& data) = 0;
    virtual const ContinentId& id(uint32_t index) const = 0;
    virtual ContinentId& id(uint32_t index) = 0;
    virtual void setId(uint32_t index, ContinentId id) const = 0;
    virtual ContinentId getContinentAt(int x, int y) const = 0;
};

class Segments : public ISegments
{
public:
    Segments(uint32_t plate_area);
    ~Segments();
    void setSegmentCreator(ISegmentCreator* segmentCreator)
    {
        _segmentCreator = segmentCreator;
    }
    void setBounds(Bounds* bounds)
    {
        _bounds = bounds;
    }
    uint32_t area();
    void reset();
    void reassign(uint32_t newarea, uint32_t* tmps);
    void shift(uint32_t d_lft, uint32_t d_top);
    uint32_t size() const;
    const SegmentData& operator[](uint32_t index) const;
    SegmentData& operator[](uint32_t index);
    void add(const SegmentData& data);
    const ContinentId& id(uint32_t index) const;
    ContinentId& id(uint32_t index);
    void setId(uint32_t index, ContinentId id) const;
    ContinentId getContinentAt(int x, int y) const;
private:
    std::vector<SegmentData> seg_data; ///< Details of each crust segment.
    ContinentId* segment;              ///< Segment ID of each piece of continental crust.
    int _area; /// Should be the same as the bounds area of the plate
    ISegmentCreator* _segmentCreator;
    Bounds* _bounds;
};

#endif
