/**
*@file Geometry.h
*@brief 一些常用几何算法 http://dev.gameres.com/Program/Abstract/Geometry.htm
*@author DionysosLai，email: 906391500@qq.com
*@version 1.0
*@data 2015-1-6 17:41
*/

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "cocos2d.h"

class Geometry
{
public:
	Geometry();
	virtual ~Geometry();

private:

public:
///@brief 叉积
///@param[in] (x1, y1)---点1， (x2, y2)---点2
///@return 叉积值
	static double vectorProduct(const double& x1, const double& y1, const double& x2, const double& y2);  // 行列式 

///@brief 判断折线的拐向，
///以线段P1为参考，；
///@param[in] p1--线段1，p2--线段2
///@pre p0为公共端点，即判断p2p0 在 p1p0的方向
///@return 1--左边 0 -- 共线 -1---右边
	static int polyLineDerection(const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p2);

///@brief 判断线段与直线是否相交
///@param[in] p0,p1--线段两个端点， q0,q1--直线两个点值
///@return true---在线段上， false---不在线段上
	static bool straightLineIsIntersect(const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& q0, const cocos2d::CCPoint& q1);

///@brief 判断点是否在矩形中
///@param[in|out] p0--点 r1--矩形中心 width--矩形宽 heigth---矩形高
///@return true--在矩形内（包括点在矩形边上）， false---在矩形外
	static bool pointIsInRect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& r1, const float& width, const float& heigth );

///@brief 判断点是否在矩形内
///@param[in] p0--点，p1,p2--矩形的一条对角线端点
///@return true--在矩形内（包括点在矩形边上）， false---在矩形外
	static bool pointIsInRect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2 );

///@brief 判断点是否在线段上
///@param[in] p0--点，p1,p2--线段两个端点
///设点为Q，线段为P1P2 ，判断点Q在该线段上的依据是：( Q - P1 ) × ( P2 - P1 ) = 0 
///且 Q 在以 P1，P2为对角顶点的矩形内。前者保证Q点在直线P1P2上，后者是保证Q点不在
///线段P1P2的延长线或反向延长线上
///@return true---在线段上， false---不在线段上
	static bool pointIsAtSegment(const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2);

///@brief 点到线段最近一个点
///@param[in] p0--要判断点， s0，s1--线段两个端点
///@return crosspoint---最近的点
	static cocos2d::CCPoint nearestPointToSegmentLine( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& s0, const cocos2d::CCPoint& s1 );

///@brief 判断两个矩形是否相交
///@param[in] aa,bb--矩形1一个对角线端点， cc,dd--矩形2一个对角线端点
///@return true---相交， false---不相交
	static bool isRectsInterSect(const cocos2d::CCPoint& aa, const cocos2d::CCPoint& bb, const cocos2d::CCPoint& cc, const cocos2d::CCPoint& dd);

///@brief 判断点是否在多边形（包括点在边上）
///@param[in] p0--要判断点， poly--多边形点集合， numberOfPoints--多边形点数量
///@return true---点在多边形内， false---点不在多边形内
	static bool pointIsInPolygon(const cocos2d::CCPoint& p0, std::vector<cocos2d::CCPoint>* point);

///@brief 判断点是否在园内
///@param[in] p0--点， r0---圆心 radius---半径
///@return true---点在园内， false--点在多边形内
	static bool pointInInCircle(const cocos2d::CCPoint& p0, const cocos2d::CCPoint& r0, const float& radius);

///@brief 判断线段与线段是否相交
///@param[in] aa,bb--线段1两个端点， cc,dd--线段1两个端点
///@return true---在线段上， false---不在线段上
	static bool segmentLineIsIntersect(const cocos2d::CCPoint& aa, const cocos2d::CCPoint& bb, const cocos2d::CCPoint& cc, const cocos2d::CCPoint& dd);

///@brief 获取两点间的角度值
///@param[in] posBegin---起点  posEnd---终点
///@return 角度值
	static float getAngle( const cocos2d::CCPoint& posBegin, const cocos2d::CCPoint& posEnd );

///@brief 判断线段圆是否相交
///@param[in] r1---圆心 radius--半径 （p1, p2）---线段两个端点值
///@return true---在圆心内， false---不在圆心内
	static bool isCircleLineCollision(const cocos2d::CCPoint& r1, const float& radius, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2);

///@brief 判断两个圆是否碰撞
///@param[in] r1--圆1圆心 radius1---圆1半径 r2--圆2圆心 radius2---圆2半径
///@return ture---相交 false---不相交
	static bool isCircleCollision(const cocos2d::CCPoint& r1, const float& radius1, const cocos2d::CCPoint& r2, const float& radius2);

///@brief 判断两个矩形是否相交
///@param[in] a0, aWidth, bHeight---矩形1中点和宽度、高度 
///@return true---相交  false---不相交
	static bool isRectsCollision(const cocos2d::CCPoint& a0, const float& aWidth, const float& aHeight, 
		const cocos2d::CCPoint& b0, const float& bWidth, const float& bHeight);

///@brief 判断线段是否在多边形内
///@param[in] a0 a1 线段两个端点 point----多边形各个端点值
///@return true--在多边形内 false---不在多边形内
	static bool isSegmentLineInPoly(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, std::vector<cocos2d::CCPoint>* point);

///@brief 计算两条共线线段的交点
///@param[in] a0a1---线段1 b0b1---线段2
///@param[out] commomPoint----交点
///@pre 两条线段必须是共线---不是点
///@return 0----两条线段不共线 1----两条线段共线且无交点 2----两条线段共线且无数个交点  3----两条线段共线且只有一个交点
	static int commonPointSegments(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, 
		const cocos2d::CCPoint& b0, const cocos2d::CCPoint& b1, cocos2d::CCPoint& commomPoint);

///@brief 计算两条线段的交点
///@param[in] a0a1---线段1 b0b1---线段2
///@param[out] commomPoint----交点
///@pre 两条线段，不是点
///@return 0----两条线段没有交点 1----两条线段共线有无数个交点 2----两条线段共线且只有一个交点 3----两条线段不共线，且只有一个交点
	static int pointOfSegments(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, 
		const cocos2d::CCPoint& b0, const cocos2d::CCPoint& b1, cocos2d::CCPoint& commomPoint);

///@brief 计算线段与圆的交点
///@param[in] a0a1---线段1 r0---圆心， radius---半径 
///@param[out] commomPoint1----交点1 commomPoint2----交点2
///@pre 两条线段，不是点
///@return 0----没有交点 1----1个交点 2----2个交点
	static int pointOfSegmentCircle(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, const cocos2d::CCPoint& r0, 
		const float& radius, cocos2d::CCPoint& commomPoint1, cocos2d::CCPoint& commomPoint2);

///@brief 计算点集的凸包---使用葛立恒（Graham）扫描法 
///@param[in] pointIn---点集  pointOut---凸包点集
///@pre 至少2点 http://zh.wikipedia.org/wiki/%E5%87%B8%E5%8C%85#.E8.91.9B.E7.AB.8B.E6.81.92.EF.BC.88Graham.EF.BC.89.E6.89.AB.E6.8F.8F.E6.B3.95
/// http://blog.csdn.net/suwei19870312/article/details/5422818 http://dev.gameres.com/Program/Abstract/Geometry.htm
///@return 
	static void tubaoCalcute(const std::vector<cocos2d::CCPoint>* pointIn, std::vector<cocos2d::CCPoint>* pointOut);
};
#endif	///<(__GEOMETRY_H__)