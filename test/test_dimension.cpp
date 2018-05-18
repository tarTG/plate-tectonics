/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "gtest/gtest.h"
#include "dimension.h"

 const Dimension testDim(20,10);

TEST(Dimension,getWidth)
{
    EXPECT_EQ(20u,testDim.getWidth());
}
TEST(Dimension,getHeight)
{
    EXPECT_EQ(10u,testDim.getHeight());
}
TEST(Dimension,getArea)
{
    EXPECT_EQ(200u,testDim.getArea());
}
TEST(Dimension,indexOf)
{
    EXPECT_EQ(42u,testDim.indexOf(Platec::Vector2D<uint32_t>(2,2)));
}
    
TEST(Dimension,yFromIndex)
{
    EXPECT_EQ(1u,testDim.yFromIndex(22));   
}
TEST(Dimension,xFromIndex)
{
     EXPECT_EQ(2u,testDim.xFromIndex(22));
}
    
TEST(Dimension,coordOF)
{
    auto coord = testDim.coordOF(22);
    EXPECT_EQ(2u,coord.x());
     EXPECT_EQ(1u,coord.y());
}

TEST(Dimension,contains)
{
    auto pInside{Platec::Vector2D<uint32_t> (5,5)};
    auto pXoutside{Platec::Vector2D<uint32_t> (25,5)};
    auto pYoutside{Platec::Vector2D<uint32_t> (5,25)};   
    EXPECT_EQ(true, testDim.contains(pInside));
    EXPECT_EQ(false, testDim.contains(pXoutside));
    EXPECT_EQ(false, testDim.contains(pYoutside));
    
    //float tests
    auto pFlotaCornerTop{Platec::Vector2D<float> (19.f,9.f)};
    auto pFlotaCornerBot{Platec::Vector2D<float> (0.f,0.f)};    
    EXPECT_EQ(true, testDim.contains(pFlotaCornerTop));
    EXPECT_EQ(true, testDim.contains(pFlotaCornerBot));
}

TEST(Dimension,grow)
{
  Dimension growDim(20,10);
  growDim.grow(Platec::Vector2D<uint32_t> (5,5));
  
    EXPECT_EQ(25u, growDim.getWidth());
    EXPECT_EQ(15u, growDim.getHeight());  
  
   
}

TEST(Dimension,getMax)
{
    EXPECT_EQ(20u, testDim.getMax());    
}

TEST(Dimension,xModPoint)
{
    auto pXoutside{Platec::Vector2D<uint32_t> (25,5)};  
    auto p = testDim.xMod(pXoutside);
    EXPECT_EQ(5u, p.x());      
}
TEST(Dimension,xMod)
{
    EXPECT_EQ(5u, testDim.xMod(25));    
}
TEST(Dimension,yModPoint)
{
    auto pYoutside{Platec::Vector2D<uint32_t> (5,15)};  
    auto p = testDim.yMod(pYoutside);
    EXPECT_EQ(5u, p.y());         
}
TEST(Dimension,yMod)
{
    EXPECT_EQ(5u, testDim.yMod(15));     
}

TEST(Dimension,normalize) 
{
    auto pYoutside{Platec::Vector2D<uint32_t> (25,15)};
    auto p = testDim.normalize(pYoutside);
    EXPECT_EQ(5u, p.x());      
     EXPECT_EQ(5u, p.y());        
}

TEST(Dimension,lineIndex) 
{
    EXPECT_EQ(100u, testDim.lineIndex(5)); 
}

TEST(Dimension,xCap)
{
    EXPECT_EQ(19u, testDim.xCap(25));     
}


TEST(Dimension,yCap) 
{
     EXPECT_EQ(9u, testDim.yCap(25));     
}

TEST(Dimension,xCapPos)
{
    EXPECT_EQ(19u, testDim.xCap(Platec::Vector2D<uint32_t> (25,15)).x());     
}


TEST(Dimension,yCapPos) 
{
     EXPECT_EQ(9u, testDim.yCap(Platec::Vector2D<uint32_t> (25,25)).y());     
}

TEST(Dimension,wrap) 
{
    auto pInside{Platec::Vector2D<uint32_t> (5,5)};
    auto pXoutsideTop{Platec::Vector2D<uint32_t> (25,5)};
    auto pYoutsideTop{Platec::Vector2D<uint32_t> (5,15)};
    auto pXoutsideBot{Platec::Vector2D<int32_t> (-5,5)};
    auto pYoutsideBot{Platec::Vector2D<int32_t> (5,-5)};    
    
    EXPECT_EQ(5u,testDim.wrap(pInside).x());
    EXPECT_EQ(5u,testDim.wrap(pInside).y()); 

    EXPECT_EQ(5u,testDim.wrap(pXoutsideTop).x());
    EXPECT_EQ(5u,testDim.wrap(pYoutsideTop).y());

    EXPECT_EQ(15,testDim.wrap(pXoutsideBot).x());
    EXPECT_EQ(5,testDim.wrap(pYoutsideBot).y());    
}
    
TEST(Dimension, getDimensions) 
{
     EXPECT_EQ(20u,testDim.getDimensions().x());
    EXPECT_EQ(10u,testDim.getDimensions().y());   
}