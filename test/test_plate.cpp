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

#include "plate.hpp"
#include "gtest/gtest.h"
#include "noise.hpp"
#include "simplexnoise.hpp"

void initializeHeightmapWithNoise(long seed, float *heightmap, const Dimension& wd)
{
    createNoise(heightmap, wd, SimpleRandom(seed), true);
    for (uint32_t i=0; i<wd.getArea(); i++) {
        if (heightmap[i]<0.0f) {
            heightmap[i] *= -1.0f;
        }
    }
}

TEST(SimpleRandom, NextRepeatability)
{
    SimpleRandom sr1(1);
    EXPECT_EQ(81414, sr1.next());
    EXPECT_EQ(1328228615, sr1.next());
    EXPECT_EQ(3215746516, sr1.next());

    SimpleRandom sr999(999);
    EXPECT_EQ(69012276, sr999.next());
    EXPECT_EQ(3490172125, sr999.next());
    EXPECT_EQ(3364058674, sr999.next());
}

TEST(Noise, SimplexRawNoiseRepeatability)
{
    EXPECT_FLOAT_EQ(-0.12851511f, raw_noise_4d(0.3f, 0.78f, 1.677f, 0.99f));
    EXPECT_FLOAT_EQ(-0.83697641f, raw_noise_4d(-0.3f, 0.78f, 1.677f, 0.99f));
    EXPECT_FLOAT_EQ(-0.5346415f, raw_noise_4d(7339.3f, 0.78f, 1.677f, 0.99f));
    EXPECT_FLOAT_EQ(0.089452535f, raw_noise_4d(0.3f, 70.78f, 1.677f, 0.0009f));
    EXPECT_FLOAT_EQ(-0.063593678f, raw_noise_4d(0.3f, 500.78f, 1.677f, 500.99f));
}

TEST(Noise, SimplexNoiseRepeatability)
{
    Dimension wd(233, 111);
    float *heightmap = new float[wd.getArea()];
    initializeHeightmapWithNoise(123, heightmap, wd);

    EXPECT_FLOAT_EQ(0.50098729f, heightmap[0]);
    EXPECT_FLOAT_EQ(0.39222634f, heightmap[1000]);
    EXPECT_FLOAT_EQ(0.51659518f, heightmap[2000]);
    EXPECT_FLOAT_EQ(0.5479334f, heightmap[5000]);
    EXPECT_FLOAT_EQ(0.59222502f, heightmap[8000]);
    EXPECT_FLOAT_EQ(0.36362505f, heightmap[11000]);
    EXPECT_FLOAT_EQ(0.57599854f, heightmap[13000]);

    delete[] heightmap;
}

TEST(CreatePlate, SquareDoesNotExplode)
{
    float *heightmap = new float[40000]; // 200 x 200
    initializeHeightmapWithNoise(678, heightmap, Dimension(200, 200));
    
    HeightMap m = HeightMap(std::vector<float>(heightmap,heightmap+40000),200,200);
    plate p = plate(123, m, Dimension(100, 3),Platec::vec2f( 50, 23), 18, Dimension(200, 200));
}

TEST(CreatePlate, NotSquareDoesNotExplode)
{
    float *heightmap = new float[80000]; // 200 x 400
    initializeHeightmapWithNoise(678, heightmap, Dimension(200, 400));
    
    HeightMap m = HeightMap(std::vector<float>(heightmap,heightmap+80000),200,400);
    plate p = plate(123, m, Dimension(100, 3),Platec::vec2f( 50, 23), 18, Dimension(200, 400));
}

// TODO test also when plate is large as world
TEST(Plate, calculateCrust)
{
    float *heightmap = new float[256 * 128];
    initializeHeightmapWithNoise(678, heightmap, Dimension(256, 128));
    
    HeightMap m = HeightMap(std::vector<float>(heightmap,heightmap+(256 * 128)),256, 128);
    plate p = plate(123, m, Dimension(100, 3),Platec::vec2f( 50, 23), 18, Dimension(256, 128));
    uint32_t x, y, index;
    float w_crust, e_crust, n_crust, s_crust;
    uint32_t w, e, n, s;

    // top left corner
    x = 0;
    y = 0;
    index = 1;
    p.calculateCrust(x,y,index,w_crust,e_crust,n_crust,s_crust,w,e,n,s);
    EXPECT_EQ(0, w);
    EXPECT_EQ(1, e);
    EXPECT_EQ(0, n);
    EXPECT_EQ(100, s);

    // bottom right corner
    x = 99;
    y = 2;
    index = 1;
    p.calculateCrust(x,y,index,w_crust,e_crust,n_crust,s_crust,w,e,n,s);
    EXPECT_EQ(298, w);
    EXPECT_EQ(200, e);
    EXPECT_EQ(199, n);
    EXPECT_EQ(99, s);

    // point in the middle
    x = 50;
    y = 1;
    index = 1;
    p.calculateCrust(x,y,index,w_crust,e_crust,n_crust,s_crust,w,e,n,s);
    EXPECT_EQ(149, w);
    EXPECT_EQ(151, e);
    EXPECT_EQ(50, n);
    EXPECT_EQ(250, s);
}

class MockSegmentData : public ISegmentData
{
public:
    MockSegmentData(uint32_t collCount, uint32_t area)
        : _collCount(collCount), _area(area)
    {

    }
    virtual void incCollCount() {
        _collCount++;
    }
    virtual void incArea() {
        _area++;
    }
    virtual void enlarge_to_contain(const Platec::vec2ui& point) {
        _enlargePoint = point;
    }
    virtual void markNonExistent() {
        throw std::runtime_error("Not implemented");
    }
    virtual void shift(const Platec::vec2ui& shiftDir) {
        throw std::runtime_error("Not implemented");
    }

    virtual uint32_t getLeft() const {
        throw std::runtime_error("Not implemented");
    }
    virtual uint32_t getRight() const {
        throw std::runtime_error("Not implemented");
    }
    virtual uint32_t getTop() const {
        throw std::runtime_error("Not implemented");
    }
    virtual uint32_t getBottom() const {
        throw std::runtime_error("Not implemented");
    }
    virtual bool isEmpty() const {
        throw std::runtime_error("Not implemented");
    }
    virtual uint32_t getArea() const {
        return _area;
    }
    virtual uint32_t collCount() const {
        return _collCount;
    }
    Platec::vec2ui enlargedPoint() {
        return _enlargePoint;
    }
private:
    uint32_t _collCount;
    uint32_t _area;
    Platec::vec2ui _enlargePoint = Platec::vec2ui(0,0);
};

class MockSegments : public ISegments
{
public:
    MockSegments(Platec::vec2ui p, ContinentId id, ISegmentData* data)
        : _p(p), _id(id), _data(data)
    {

    }
    virtual const uint32_t getArea() const {
        throw std::runtime_error("Not implemented");
    }
    virtual void reset() {
        throw std::runtime_error("Not implemented");
    }
    virtual void reassign(const uint32_t newarea,const std::vector<uint32_t>& tmps) {
        throw std::runtime_error("Not implemented");
    }
    virtual void shift(const Platec::vec2ui& dir) {
        throw std::runtime_error("Not implemented");
    }
    virtual const uint32_t size() const {
        throw std::runtime_error("(MockSegments::size) Not implemented");
    }
    virtual const ISegmentData& operator[](uint32_t index) const {
        if (index == _id) {
            return *_data;
        } else {
            throw std::runtime_error("(MockSegments::operator[]) Unexpected call");
        }
    }
    virtual ISegmentData& operator[](uint32_t index) {
        if (index == _id) {
            return *_data;
        } else {
            throw std::runtime_error("(MockSegments::operator[]) Unexpected call");
        }
    }
    virtual void add(const SegmentData& data) {
        throw std::runtime_error("Not implemented");
    }
    virtual const ContinentId& id(uint32_t index) const {
        throw std::runtime_error("(MockSegments::id) Not implemented");
    }
    virtual ContinentId& id(uint32_t index) {
        throw std::runtime_error("(MockSegments::id) Not implemented");
    }
    virtual void setId(uint32_t index, ContinentId id) {
        throw std::runtime_error("Not implemented");
    }
    virtual ContinentId getContinentAt(const Platec::vec2ui& point,
                                    const Dimension& worldDimension) const {
        if (_p == point) {
            return _id;
        } else {
            throw std::runtime_error("(MockSegments::getContinentAt) Unexpected call");
        }
    }
private:
    Platec::vec2ui _p;
    ContinentId _id;
    ISegmentData* _data;
};

TEST(Plate, addCollision)
{
    float *heightmap = new float[256 * 128];
    initializeHeightmapWithNoise(678, heightmap, Dimension(256, 128));
    
    HeightMap m = HeightMap(std::vector<float>(heightmap,heightmap+(256 * 128)),256, 128);
    plate p = plate(123, m, Dimension(100, 3),Platec::vec2f( 50, 23), 18, Dimension(256, 128));

    MockSegmentData* mSeg = new MockSegmentData(7, 789);
    std::shared_ptr<MockSegments> mSegments = std::make_shared<MockSegments>(Platec::vec2ui(123, 78), 99, mSeg);
    p.injectSegments(mSegments);

    uint32_t area = p.addCollision(Platec::vec2ui(123, 78));

    EXPECT_EQ(789, area);
    EXPECT_EQ(8, mSeg->collCount());
}

class MockSegments2 : public ISegments
{
public:
    MockSegments2(Platec::vec2ui p, ContinentId id, ISegmentData* data, uint32_t index)
        : _p(p), _id(id), _data(data), _index(index)
    {

    }
    virtual const uint32_t getArea() const{
        throw std::runtime_error("(MockSegments2::area) Not implemented");
    }
    virtual void reset() {
        throw std::runtime_error("(MockSegments2::reset) Not implemented");
    }
    virtual void reassign(const uint32_t newarea,const std::vector<uint32_t>& tmps) {
        throw std::runtime_error("(MockSegments2::reassign) Not implemented");
    }
    virtual void shift(const Platec::vec2ui& dir) {
        throw std::runtime_error("(MockSegments2::shift) Not implemented");
    }
    virtual const uint32_t size() const {
        throw std::runtime_error("(MockSegments2::size) Not implemented");
    }
    virtual const ISegmentData& operator[](uint32_t index) const {
        if (index == _id) {
            return *_data;
        } else {
            throw std::runtime_error("(MockSegments2::operator[]) Unexpected call");
        }
    }
    virtual ISegmentData& operator[](uint32_t id) {
        if (id == _id) {
            return *_data;
        } else {
            throw std::runtime_error(std::string("(MockSegments2::operator[]) Unexpected call with id "
                                       + Platec::to_string(id)));
        }
    }
    virtual void add(const SegmentData& data) {
        throw std::runtime_error("(MockSegments2::add) Not implemented");
    }
    virtual const ContinentId& id(uint32_t index) const {
        if (_index == index) return _id;
        throw std::runtime_error(
            std::string("(MockSegments2::id) Unexpected value ")
            + Platec::to_string(index)
            + " expected was "
            + Platec::to_string(_index));
    }
    virtual ContinentId& id(uint32_t index) {
        if (_index == index) return _id;
        throw std::runtime_error(
            std::string("(MockSegments2::id) Unexpected value ")
            + Platec::to_string(index)
            + " expected was "
            + Platec::to_string(_index));
    }
    virtual void setId(uint32_t index, ContinentId id) {
        if (_index == index) {
            _id = id;
        } else {
            throw std::runtime_error(std::string("(MockSegments2::setId) Unexpected call with index ")
                                + Platec::to_string(index));
        }
    }
    virtual ContinentId getContinentAt(const Platec::vec2ui& point,
                                    const Dimension& worldDimension) const {
        if (_p == point) {
            return _id;
        } else {
            throw std::runtime_error("(MockSegments2::getContinentAt) Unexpected call");
        }
    }
private:
    Platec::vec2ui _p;
    ContinentId _id;
    ISegmentData* _data;
    uint32_t _index;
};

TEST(Plate, addCrustByCollision)
{
    const uint32_t worldWidth = 256;
    const uint32_t worldHeight = 128;

    float *heightmap = new float[worldWidth * worldHeight];
    Dimension wd(worldWidth, worldHeight);
    initializeHeightmapWithNoise(1, heightmap, wd);

    // Suppose the plate start at 170, 70 and ends at 250, 125
    HeightMap m = HeightMap(std::vector<float>(heightmap,heightmap+(256 * 128)),256, 128);
    plate p = plate(123, m, Dimension(80, 55),Platec::vec2f( 170, 70), 18, Dimension(256, 128));
    // the point of collision in world coordinates
    const uint32_t worldPointX = 240;
    const uint32_t worldPointY = 120;
    // the point of collision in plate coordinates
    const uint32_t platePointX = 70;
    const uint32_t platePointY = 50;
    uint32_t indexInPlate = platePointY * 80 + platePointX;

    MockSegmentData* mSeg = new MockSegmentData(7, 789);
    std::shared_ptr<MockSegments2> mSegments = std::make_shared<MockSegments2>(Platec::vec2ui(worldPointX, worldPointY), 99, mSeg, indexInPlate);
    p.injectSegments(mSegments);

    uint32_t timestampIn_240_120before = p.getCrustTimestamp(worldPointX, worldPointY);
    float crustIn_240_120before = p.getCrust(worldPointX, worldPointY);

    // Assumptions:
    // the point is in the plate bounds

    p.addCrustByCollision(
        worldPointX, worldPointY, // Point of impact
        0.8f, // Amount of crust
        123, // Current age
        99); // Active continent

    // Age of the point should be updated
    uint32_t timestampIn_240_120after = p.getCrustTimestamp(worldPointX, worldPointY);
    ASSERT_EQ(true, timestampIn_240_120after > timestampIn_240_120before);
    ASSERT_EQ(true, timestampIn_240_120after < 123 );

    // Crust should be increased
    float crustIn_240_120after = p.getCrust(worldPointX, worldPointY);
    EXPECT_FLOAT_EQ(crustIn_240_120before + 0.8f, crustIn_240_120after);

    // The activeContinent should now owns the point
    EXPECT_EQ(99, mSegments->getContinentAt(Platec::vec2ui(worldPointX, worldPointY),
                                                wd));

    // The activeContinent should contains the point
 //   ASSERT_EQ(false, NULL==mSeg->enlargedPoint());
    EXPECT_EQ(70, mSeg->enlargedPoint().x());
    EXPECT_EQ(50, mSeg->enlargedPoint().y());

    // The activeContinent are should be increased
    EXPECT_EQ(790, mSeg->getArea());
}

TEST(Plate, addCrustBySubduction)
{
    const uint32_t worldWidth = 256;
    const uint32_t worldHeight = 128;

    float *heightmap = new float[worldWidth * worldHeight];
    initializeHeightmapWithNoise(1, heightmap, Dimension(worldWidth, worldHeight));

    HeightMap m = HeightMap(std::vector<float>(heightmap,heightmap+(256 * 128)),256, 128);
    plate p = plate(123, m, Dimension(80, 55),Platec::vec2f( 170, 70), 18, Dimension(256, 128));
    // Suppose the plate start at 170, 70 and ends at 250, 125
    // the point of collision in world coordinates
    const uint32_t worldPointX = 240;
    const uint32_t worldPointY = 120;
    // the point of collision in plate coordinates
    const uint32_t platePointX = 70;
    const uint32_t platePointY = 50;
    uint32_t indexInPlate = platePointY * 80 + platePointX;

    MockSegmentData* mSeg = new MockSegmentData(7, 789);
    std::shared_ptr<MockSegments2>  mSegments = std::make_shared<MockSegments2>(Platec::vec2ui(worldPointX, worldPointY), 99, mSeg, indexInPlate);
    p.injectSegments(mSegments);

    uint32_t timestampIn_240_120before = p.getCrustTimestamp(worldPointX, worldPointY);
    float crustIn_240_120before = p.getCrust(worldPointX, worldPointY);

    float massBefore = p.getMass();

    // Assumptions:
    // the point is in the plate bounds
    float dx = 0.0f;
    float dy = 0.0f;
    p.addCrustBySubduction(
        worldPointX, worldPointY, // Point of impact
        0.8f, // Amount of crust
        123, // Current age
        dx, dy); // Direction of the subducting plate

    // Crust should be increased
    float crustIn_240_120after = p.getCrust(worldPointX, worldPointY);
    EXPECT_FLOAT_EQ(crustIn_240_120before + 0.8f, crustIn_240_120after);

    // The mass should be increased
    float massAfter = p.getMass();
    EXPECT_EQ(massBefore + 0.8f, massAfter);

    // Age of the point should be updated
    uint32_t timestampIn_240_120after = p.getCrustTimestamp(worldPointX, worldPointY);
    ASSERT_EQ(true, timestampIn_240_120after > timestampIn_240_120before);
    ASSERT_EQ(true, timestampIn_240_120after < 123 );
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}