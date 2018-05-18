/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include "gtest/gtest.h"
#include "vector2D.h"

using namespace Platec;

const  Vector2D<int> intVec(10,20); 
const  Vector2D<float> floatVec(10.f,20.f); 
  
TEST(Vector2D, x)
{
    EXPECT_EQ(10,intVec.x());
    EXPECT_EQ(10.f,floatVec.x());
}
TEST(Vector2D,y)
{
    EXPECT_EQ(20,intVec.y());
    EXPECT_EQ(20.f,floatVec.y());    
}
TEST(Vector2D,length) 
{
    EXPECT_EQ(22.360679775f,intVec.length());
    EXPECT_EQ(22.360679775f,floatVec.length()); 
}
    
TEST(Vector2D,opMinus)
{
    auto resultInt = intVec - Vector2D<int>(10,20);
    EXPECT_EQ(0,resultInt.x());
    EXPECT_EQ(0,resultInt.y());  
    auto resultFloat = floatVec - Vector2D<float>(10.f,20.f);
    EXPECT_EQ(0.f,resultFloat.x());
    EXPECT_EQ(0.f,resultFloat.y());      
}
    
TEST(Vector2D,dotProduct)  
{
    EXPECT_EQ(500,intVec.dotProduct(Vector2D<int>(10,20)));
    EXPECT_EQ(500.f,floatVec.dotProduct(Vector2D<float>(10.f,20.f))); 
}
    
TEST(Vector2D,opEquals)
 {
    ASSERT_EQ(true,intVec == Vector2D<int>(10,20) );
    ASSERT_EQ(false,intVec == Vector2D<int>(9,20) );   
    ASSERT_EQ(true,floatVec == Vector2D<float>(10,20) );
    ASSERT_EQ(false,floatVec == Vector2D<float>(9,20) );      
 }
    
TEST(Vector2D,opPuls) 
{
    auto resultInt = intVec + Vector2D<int>(10,20);
    EXPECT_EQ(20,resultInt.x());
    EXPECT_EQ(40,resultInt.y());  
    auto resultFloat = floatVec + Vector2D<float>(10.f,20.f);
    EXPECT_EQ(20.f,resultFloat.x());
    EXPECT_EQ(40.f,resultFloat.y());       
}
    
TEST(Vector2D,opMult) 
{
    auto resultInt = intVec * 10;
    EXPECT_EQ(100,resultInt.x());
    EXPECT_EQ(200,resultInt.y());  
    auto resultFloat = floatVec * 10;
    EXPECT_EQ(100.f,resultFloat.x());
    EXPECT_EQ(200.f,resultFloat.y());           
}
    
TEST(Vector2D,opDiv ) 
{
    auto resultInt = intVec / 10;
    EXPECT_EQ(1,resultInt.x());
    EXPECT_EQ(2,resultInt.y());  
    auto resultFloat = floatVec / 10;
    EXPECT_EQ(1.f,resultFloat.x());
    EXPECT_EQ(2.f,resultFloat.y());            
}    
    
TEST(Vector2D,shift) 
{
    Vector2D<int> intVecShift(10,20); 
    intVecShift.shift(Vector2D<int>(10,20));
    EXPECT_EQ(20,intVecShift.x());
    EXPECT_EQ(40,intVecShift.y());  
    
    
    Vector2D<float> floatVecShift(10.f,20.f); 
    floatVecShift.shift(Vector2D<float>(10.f,20.f));
    EXPECT_EQ(20.f,floatVecShift.x());
    EXPECT_EQ(40.f,floatVecShift.y());      
    
} 
    
TEST(Vector2D,getTopPosition) 
{
    auto resultInt = intVec.getTopPosition();
    EXPECT_EQ(10,resultInt.x());
    EXPECT_EQ(19,resultInt.y());  
    auto resultFloat = floatVec.getTopPosition();
    EXPECT_EQ(10.f,resultFloat.x());
    EXPECT_EQ(19.f,resultFloat.y());         
}
    
TEST(Vector2D,getBottomPosition) 
{
    auto resultInt = intVec.getBottomPosition();
    EXPECT_EQ(10,resultInt.x());
    EXPECT_EQ(21,resultInt.y());  
    auto resultFloat = floatVec.getBottomPosition();
    EXPECT_EQ(10.f,resultFloat.x());
    EXPECT_EQ(21.f,resultFloat.y());         
}
    
TEST(Vector2D,getLeftPosition) 
{
    auto resultInt = intVec.getLeftPosition();
    EXPECT_EQ(9,resultInt.x());
    EXPECT_EQ(20,resultInt.y());  
    auto resultFloat = floatVec.getLeftPosition();
    EXPECT_EQ(9.f,resultFloat.x());
    EXPECT_EQ(20.f,resultFloat.y());         
}
     
TEST(Vector2D,getRightPosition) 
{
    auto resultInt = intVec.getRightPosition();
    EXPECT_EQ(11,resultInt.x());
    EXPECT_EQ(20,resultInt.y());  
    auto resultFloat = floatVec.getRightPosition();
    EXPECT_EQ(11.f,resultFloat.x());
    EXPECT_EQ(20.f,resultFloat.y());     
}    