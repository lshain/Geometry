/**
*@file Geometry.cpp
*@brief 一些常用几何算法 http://dev.gameres.com/Program/Abstract/Geometry.htm
*@author DionysosLai，email: 906391500@qq.com
*@version 1.0
*@data 2015-1-6 17:41
*/

#include "Geometry.h"
#include <limits>

USING_NS_CC;

Geometry::Geometry()
{

}

Geometry::~Geometry()
{

}

double Geometry::vectorProduct( const double& x1, const double& y1, const double& x2, const double& y2 )
{
	return (x1*y2-x2*y1);	///< 叉积 （与数学叉积相反）
}

int Geometry::polyLineDerection( const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p2 )
{
	float vectorProductResult = 0.0f;
	int derection = 0;

	vectorProductResult = (float)vectorProduct(p1.x-p0.x, p1.y-p0.y, p2.x-p0.x, p2.y-p0.y);

	if (vectorProductResult < 1 && vectorProductResult > -1)
	{
		derection = 0;	///< 共线
	}
	else if (vectorProductResult >= 1)
	{
		derection = -1;	///< p2p0在p1p0的右手边
	}
	else if(vectorProductResult <= -1)
	{
		derection = 1;	///< p2p0在p1p0的左手边
	}
	
	return derection;
}

bool Geometry::straightLineIsIntersect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& q0, const cocos2d::CCPoint& q1 )
{
	/// 先判断q0与q1能否组成一条直线
	if (!q0.equals(q1))
	{
		/// 当q0 == p0 q1 == p1时，结果为0。 只需要判断线段是否跨立直线即可。
		if (0 <= vectorProduct(p0.x-q0.x, p0.y-q0.y, q1.x-q0.x, q1.y-q0.y) * 
			vectorProduct(q1.x-q0.x, q1.y-q0.y, p1.x-q0.x, p1.y-q0.y))
		{
			/* CCLOG("相交");*/
			return true;
		}

		/* CCLOG("不相交");*/
		return false;
	}

	/*	CCLOG("点q0点q1不能构成一条直线");*/
	return false;
}

bool Geometry::pointIsInRect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& r1, const float& width, const float& heigth )
{
	float xLeft = r1.x - width/2.f;
	float xRigth = r1.x + width/2.f;
	float yUp = r1.y + heigth/2.f;
	float yBottom = r1.y - heigth/2.f;

	if (p0.x >= xLeft && p0.x<= xRigth && p0.y >= yBottom && p0.y <= yUp)
	{
		return true;
	}

	return false;
}

bool Geometry::pointIsInRect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2 )
{
	float xMax = 0, xMin, yMax = 0, yMin = 0;
	xMax = p1.x > p2.x ? p1.x : p2.x;
	xMin = p1.x > p2.x ? p2.x : p1.x;
	yMax = p1.y > p2.y ? p1.y : p2.y;
	yMin = p1.y > p2.y ? p2.y : p1.y;

	if ( p0.x >= xMin && p0.x <= xMax && p0.y >= yMin && p0.y <= yMax)
	{
		/* CCLOG("Point is at the recangle.");*/
		return true;
	}
	/* CCLOG("Point isn't at the recangle.");*/
	return false;
}


bool Geometry::pointIsAtSegment( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2 )
{
	/// 先判断是否点在以p1 p2为对角线的矩形内
	if (pointIsInRect(p0, p1, p2))
	{
		/// 判断p1p0, p2p0是否共线
		if (0 == polyLineDerection(p0, p1, p2))
		{
			return true;
		}

		/* CCLOG("Point isn't at the line.");*/
		return false;
	}
	//CCLOG("Point isn't at the line.");
	return false;
}

cocos2d::CCPoint Geometry::nearestPointToSegmentLine( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& s0, const cocos2d::CCPoint& s1 )
{
	/// 判断线段是否是一个点
	if (s0.equals(s1))
	{
		return s0;
	}

	/// 初始垂足为原点；
	CCPoint crossPoint = CCPointZero;
	do 
	{
		/// 判断线段是否平行于x轴
		if (s0.y == s1.y)
		{
			crossPoint = ccp(p0.x, s0.y);
			break;
		}
		/// 判断线段是否平行于y轴
		if (s0.x == s1.x)
		{
			crossPoint = ccp(s0.x, p0.y);
			break;
		}
		/// 如果线段不是特殊情况，则只能采用直线方程方式联立求解
		float k = (s1.y - s0.y)/(s1.x - s0.x);		///< 求得斜率
		/// 线段直线方程：	y = k* ( x - s0.x) + s0.y
		/// 垂线方程为：	y = (-1/k) * (x - p0.x) + p0.y 。
		/// 联立两直线方程解得
		float x = ( k*k * s0.x + k * (p0.y - s0.y ) + p0.x ) / ( k*k + 1);
		float y = k * ( x - s0.x) + s0.y;
		crossPoint = ccp(x, y);
		break;
	} while (0);

	/// 判断垂直是否在线段上
	if (pointIsAtSegment(crossPoint, s0, s1))
	{
		return crossPoint;
	}
	else
	{
		/// 如果不在则计算两端点到垂足的距离，选择距离垂足较近的端点返回。
		float distance1 = ccpDistance(crossPoint, s0);
		float distance2 = ccpDistance(crossPoint, s1);
		if (distance1 < distance2)
		{
			return s0;
		}
		else
		{
			return s1;
		}
	}
}

bool Geometry::isRectsInterSect( const cocos2d::CCPoint& aa, const cocos2d::CCPoint& bb, const cocos2d::CCPoint& cc, const cocos2d::CCPoint& dd )
{
	CCPoint	centre1 = ccp((aa.x + bb.x)/2.0f, (aa.y + bb.y)/2.0f);	///< 获得中点值
	CCPoint	centre2 = ccp((cc.x + dd.x)/2.0f, (cc.y + dd.y)/2.0f);
	float	lengthX	= abs(centre1.x - centre2.x);	///< 获得两个矩形中心的距离 
	float	lengthY	= abs(centre1.y - centre2.y); 
	float	lengthRect1X	= abs(aa.x - bb.x);		///< 获得两个矩形长和宽
	float	lengthRect1Y	= abs(aa.y - bb.y);
	float	lengthRect2X	= abs(cc.x - dd.x);
	float	lengthRect2Y	= abs(cc.y - dd.y);

	/// 这里减去1是调整误差用的。
	return  (lengthX < (lengthRect1X + lengthRect2X)/2.0f-1 && lengthY < (lengthRect1Y + lengthRect2Y)/2.0f-1) ? true : false;
}

bool Geometry::pointIsInPolygon( const cocos2d::CCPoint& p0, std::vector<cocos2d::CCPoint>* point )
{
	unsigned int count  = 0;		///< 用来标记射线L与多边形的交点数；
	cocos2d::CCSize	winsize = CCDirector::sharedDirector()->getWinSize();
	/// 已点p0向左向右做一条射线L；
	CCPoint leftPoint = ccp(-100.f, p0.y);
	CCPoint rightPoint = p0;


	/// 判断每条边
	unsigned int numberOfPoints = point->size();
	for (unsigned int i = 0; i < numberOfPoints; i++)
	{
		/// 先判断点p0是否在边s上；
		CCPoint s0 = point->at(i);
		CCPoint s1 = point->at((i+1)%(numberOfPoints));
		if (pointIsAtSegment(p0, s0, s1))
		{
			/*			CCLOG("Point is at the %dth line", i);*/

			return true;
		}

		/// 判断边s是否是平行线；
		if (s0.y != s1.y)
		{		
			do 
			{
				/// 判断边s的是否有端点在L上 同时 再判断该点是否是边s纵坐标较大的一个点
				if (pointIsAtSegment(s0, leftPoint, rightPoint))
				{
					if (s0.y > s1.y)
					{
						count += 1;
					}
					break;
				}	
				if (pointIsAtSegment(s1, leftPoint, rightPoint))
				{
					if (s0.y < s1.y)
					{
						count += 1;
					}

					break;
				}	

				/// 如果边s没有端点在L上，则判断s与L是否相交
				if (segmentLineIsIntersect(leftPoint, rightPoint, s0, s1))
				{
					count += 1;
				}	
			} while (0);
		}
	}
	if (count%2 == 1)
	{
		//		CCLOG("true");
		return true;
	}
	else
	{
		//		CCLOG("false");
		return false;
	}

	// 	if (1 == count%2)
	// 	{
	// /*		CCLOG("Point is in polygon!");*/
	// 		return true;
	// 	}
	// 	else
	// 	{
	// /*		CCLOG("Point is in polygon!");*/
	// 		return false;
	// 	}
}

bool Geometry::pointInInCircle( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& r0, const float& radius )
{
	if (ccpDistance(p0, r0) < radius)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Geometry::segmentLineIsIntersect( const cocos2d::CCPoint& aa, const cocos2d::CCPoint& bb, const cocos2d::CCPoint& cc, const cocos2d::CCPoint& dd )
{
	/// 以两条线段形成的矩形不重合，说明两条线段必然不相交
	if (!isRectsInterSect(aa, bb, cc, dd))
	{
		return false;
	}

	/// 必须二者互相跨立 注意"="的情况。
	if (0 < vectorProduct(aa.x-cc.x, aa.y-cc.y, dd.x-cc.x, dd.y-cc.y) * 
		vectorProduct(dd.x-cc.x, dd.y-cc.y, bb.x-cc.x, bb.y-cc.y) &&
		0 < vectorProduct(cc.x-aa.x, cc.y-aa.y, bb.x-aa.x, bb.y-aa.y) * 
		vectorProduct(bb.x-aa.x, bb.y-aa.y, dd.x-aa.x, dd.y-aa.y))
	{
		return true;
	}
	return false; 
}

float Geometry::getAngle( const cocos2d::CCPoint& posBegin, const cocos2d::CCPoint& posEnd )
{
	CCPoint deltaPoint = posEnd - posBegin;

	float angleRadians = atanf(deltaPoint.y / deltaPoint.x);
	float angleDegrees = CC_RADIANS_TO_DEGREES(angleRadians);
	float cocosAngle = angleDegrees;
	if (deltaPoint.x >= 0 && deltaPoint.y >= 0)	///< 第一象限
	{
		cocosAngle = cocosAngle;
	}
	if (deltaPoint.x < 0 && deltaPoint.y >= 0)	///< 第二象限
	{
		cocosAngle = 180.f + cocosAngle;
	}
	if (deltaPoint.x < 0 && deltaPoint.y < 0)	///< 第三象限
	{
		cocosAngle = 180.f + cocosAngle;
	}
	if (deltaPoint.x >= 0 && deltaPoint.y < 0)	///< 第四象限
	{
		cocosAngle = 360.f + cocosAngle;
	}	

	//	CCLOG("%f", cocosAngle);
	return cocosAngle;
}

bool Geometry::isCircleLineCollision( const cocos2d::CCPoint& r1, const float& radius, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2 )
{
	/// 判断线段是否是一个点
	float length = 0.f;
	if (p1.equals(p2))
	{
		length = ccpDistance(r1, p1);
	}

	/// 初始垂足为0；
	CCPoint crossPoint = CCPointZero;
	do 
	{
		/// 判断线段是否平行于x轴
		if (p1.y == p2.y)
		{
			crossPoint = ccp(r1.x, p1.y);
			break;
		}
		/// 判断线段是否平行于y轴
		if (p1.x == p2.x)
		{
			crossPoint = ccp(p1.x, r1.y);
			break;
		}
		/// 如果线段不是特殊情况，则只能采用直线方程方式联立求解
		float k = (p2.y - p1.y)/(p2.x - p1.x);		///< 求得斜率
		/// 线段直线方程：	y = k* ( x - s0.x) + s0.y
		/// 垂线方程为：	y = (-1/k) * (x - p0.x) + p0.y 。
		/// 联立两直线方程解得
		float x = ( k*k * p1.x + k * (r1.y - p1.y ) + r1.x ) / ( k*k + 1);
		float y = k * ( x - p1.x) + p1.y;
		crossPoint = ccp(x, y);

		/// 判断垂直是否在线段上

		if (pointIsAtSegment(crossPoint, p1, p2))
		{
			/*		return crossPoint;*/
		}
		else
		{
			/// 如果不在则计算两端点到垂足的距离，选择距离垂足较近的端点返回。
			float distance1 = ccpDistance(crossPoint, p1);
			float distance2 = ccpDistance(crossPoint, p2);
			if (distance1 < distance2)
			{
				crossPoint = p1;
			}
			else
			{
				crossPoint= p2;
			}
		}

		length = ccpDistance(r1, crossPoint);

		break;
	} while (0);

	if (length < radius)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Geometry::isCircleCollision( const cocos2d::CCPoint& r1, const float& radius1, const cocos2d::CCPoint& r2, const float& radius2 )
{
	float circleDistance = ccpDistance(r1, r2);
	if (circleDistance < (radius1+radius2))
	{
		return true;
	}
	return false;
}

bool Geometry::isRectsCollision( const cocos2d::CCPoint& a0, const float& aWidth, const float& aHeight, const cocos2d::CCPoint& b0, const float& bWidth, const float& bHeight )
{
	float	lengthX	= abs(a0.x - b0.x);	///< 获得两个矩形中心的距离 
	float	lengthY	= abs(a0.y - b0.y); 

	return  (lengthX < (aWidth + bWidth)/2.0f && lengthY < (aHeight + bHeight)/2.0f) ? true : false;
}

bool Geometry::isSegmentLineInPoly( const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, std::vector<cocos2d::CCPoint>* point )
{
		/*
	if 线端PQ的端点不都在多边形内 
		then return false;
	点集pointSet初始化为空;
	for 多边形的每条边s
		do if 线段的某个端点在s上
			then 将该端点加入pointSet;
		else if s的某个端点在线段PQ上
			then 将该端点加入pointSet;
		else if s和线段PQ相交 // 这时候已经可以肯定是内交了
			then return false;
		将pointSet中的点按照X-Y坐标排序;
		for pointSet中每两个相邻点 pointSet[i] , pointSet[ i+1]
		do if pointSet[i] , pointSet[ i+1] 的中点不在多边形中
			then return false;
		return true;
		*/
	if (!pointIsInPolygon(a0, point) || !pointIsInPolygon(a1, point))		///< 首先判断线段a0a1
	{
		return false;
	}
	
	/// 判断每条边
	unsigned int numberOfPoints = point->size();
	std::vector<CCPoint> pointSet;
	for (unsigned int i = 0; i < numberOfPoints-1; i++)
	{
		/// 先判断点p0是否在边s上；
		CCPoint s0 = point->at(i);
		CCPoint s1 = point->at((i+1)%(numberOfPoints));
		do 
		{
			/// 判断线段a0a1端点是否在s0s1上
			if (pointIsAtSegment(a0, s0, s1))
			{
				pointSet.push_back(a0);
				break;
			}
			if (pointIsAtSegment(a1, s0, s1))
			{
				pointSet.push_back(a1);
				break;
			}
			/// 判断线段s0s1是否在a0a1上
			if (pointIsAtSegment(s0, a0, a1))
			{
				pointSet.push_back(s0);
				break;
			}
			if (pointIsAtSegment(s1, a0, a1))
			{
				pointSet.push_back(s1);
				break;
			}
			/// 判断a0a1是否与s0s1相交，相交---则必然内交，必然判定线段a0a1不在多边形内
			if (segmentLineIsIntersect(a0, a1, s0, s1))
			{
				return false;
			}
		} while (0);
	}

	unsigned int numOfPointSet = pointSet.size();
	if (0 == numOfPointSet )
	{
		return true;
	}

	/// 排序poisntSet, 有小到大，优先级为x
	CCAssert(0 == numOfPointSet%2, "The pointSet's points num should be 2 times!");
	for (unsigned int i = 0; i < numOfPointSet-1; ++i)
	{
		/// 冒泡排序
		for (unsigned int j = 0; j < numOfPointSet; ++j)
		{
			CCPoint p0 = pointSet.at(i);
			CCPoint p1 = pointSet.at(j);
			if (p0.x > p1.x)
			{
				std::swap(pointSet.at(i), pointSet.at(j));
			}
			else if (p0.x == p1.x)
			{
				if (p0.y > p1.y)
				{
					std::swap(pointSet.at(i), pointSet.at(j));
				}
			}

		}
	}

	

	/// 判断相邻交点集的中点是否在多边形内
	for (unsigned int i = 0; i < numOfPointSet-1; ++i)
	{
		CCPoint p0 = pointSet.at(i);
		CCPoint p1 = pointSet.at((i+1)%numOfPointSet);
		CCPoint midP = ccpMidpoint(p0, p1);
		if (!pointIsInPolygon(midP, point))
		{
			return false;
		}
	}
	return true;
}

int Geometry::commonPointSegments( const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, const cocos2d::CCPoint& b0, const cocos2d::CCPoint& b1, cocos2d::CCPoint& commomPoint )
{
	CCAssert(a1.equals(a0) || b0.equals(b1), "a0 should not be equal to a1, this is same to b0 and b1!");

	CCPoint p0 = a0, p1 = a1;
	CCPoint q0 = b0, q1 = b1;
	commomPoint = CCPointZero;
	int commomType = 0;
	/// 先判断2条线段是否共线
	if (0 != polyLineDerection(p0, q0, q1) || 0 != polyLineDerection(p0, q0, q1))
	{
		commomType = 0;
	}
	else
	{
		/// 设置2条线段比较长的为p0p1
		float l0 = ccpLength(ccpSub(p1, p0));
		float l1 = ccpLength(ccpSub(q0, q1));
		if (l0 < l1)
		{
			/// 线段p0p1和线段q0q1调换
			CCPoint point;
			point = p0;
			p0 = q0;
			q0 = point;
			point = p1;
			p1 = q1;
			q1 = p1;
		}

		/// 判断线段q0q1两点是否在线段p0p1上
		int m = 0, n = 0;
		/// 判断q0在线段p0p1上情况
		if (pointIsAtSegment(q0, p0, p1))
		{
			if (q0.equals(p0))
			{
				m = 0;
			}
			else if (q0.equals(p1))
			{
				m = 0;
			}
			m = 1;
		}
		else
		{
			m = 2;
		}
		/// 判断q1在线段p0p1上情况
		if (pointIsAtSegment(q1, p0, p1))
		{
			if (q1.equals(p0))
			{
				n = 0;
			}
			else if (q1.equals(p1))
			{
				n = 0;
			}
			n = 1;
		}
		else
		{
			n = 2;
		}

		switch (m)
		{
		case 0:
			{
				switch (n)
				{
				case 0:
				case 1:
					commomType = 2;
					break;
				case 2:
					{
						commomType = 3;
						commomPoint = q0;
					}
					break;
				default:
					break;
				}
			}
			break;
		case 1:
			{
				commomType = 2;
			}
			break;
		case 2:
			{
				switch (n)
				{
				case 0:
					{
						commomType = 3;
						commomPoint = q1;
					}
					break;
				case 1:
					{
						commomType = 2;
					}
					break;
				case 2:
					{
						commomType = 1;
					}
					break;
				default:
					break;
				}
			}
			break;
		default:
			break;
		}
	}
	return commomType;
}

int Geometry::pointOfSegments( const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, const cocos2d::CCPoint& b0, const cocos2d::CCPoint& b1, cocos2d::CCPoint& commomPoint )
{
	CCAssert(a1.equals(a0) || b0.equals(b1), "a0 should not be equal to a1, this is same to b0 and b1!");

	CCPoint p0 = a0, p1 = a1;
	CCPoint q0 = b0, q1 = b1;
	int commomType = 0;
	commomPoint = CCPointZero;
	/// 首先判断2条线段是否相交，不相交，自然没有交点
	if (segmentLineIsIntersect(p0, p1, q0, q1))
	{
		commomType = 0;
	}
	else
	{
		/// 两条线段相交，将两条线段当做直线处理
		//////////////////////////////////////////////////////////////////////////
		/// 情况1 线段p0p1平行于y轴
		if (p0.x == p1.x)
		{
			/// 情况1.1 线段q0q1平行于y轴
			if (q0.x == q1.x)
			{
				/// 判断p0p1与q0q1是否可共线
				if (p0.x == q0.x)
				{
					commomType = commonPointSegments(p0, p1, q0, q1, commomPoint);
				}
				else
				{
					commomType = 0;
				}
			}
			/// 情况1.2 若q0q1不平行于Y轴，则交点横坐标为p0的横坐标，代入到q0q1的直线方程中可以计算出交点纵坐标；
			else
			{
				commomType = 3;
				commomPoint.x = p0.x;
				commomPoint.y = (q1.y-q0.y)/(q1.x-q0.x)*(p0.x-q0.x) + q0.y;
			}
		}
		/// 情况2 p0和p1横坐标不同，但是q0和q1横坐标相同，即q0q1平行于Y轴，则交点横坐标为q0的横坐标，代入到p0p1的直线方程中可以计算出交点纵坐标；
		else if (p0.x != p1.x && q0.x == q1.x)
		{
			commomType = 3;
			commomPoint.x = q0.x;
			commomPoint.y = (p1.y-p0.y)/(p1.x-p0.x)*(q0.x-p0.x) + p0.y;
		}
		/// 情况3 如果p0和p1纵坐标相同，即p0p1平行于X轴
		else if (p0.y == p1.y)
		{
			///< 情况3.1  若q0q1也平行于X轴
			if (q0.y == q1.y)
			{
				if (p0.y == q0.y)
				{
					commomType = commonPointSegments(p0, p1, q0, q1, commomPoint);
				}
				else
				{
					commomType = 0;
				}
			}
			///< 情况3.2 若q0q1不平行于X轴，则交点纵坐标为p0的纵坐标，代入到q0q1的直线方程中可以计算出交点横坐标；
			else
			{
				commomType = 3;
				commomPoint.y = p0.y;
				commomPoint.x = (p0.y+q0.y)*(q1.x-q0.x)/(q1.y-q0.y) + q0.x;
			}
		}
		/// 情况4 如果p0和p1纵坐标不同，但是q0和q1纵坐标相同，即q0q1平行于X轴，则交点纵坐标为q0的纵坐标，代入到p0p1的直线方程中可以计算出交点横坐标；
		else if (p0.y != p0.y && q0.y == q1.y)
		{
			commomType = 3;
			commomPoint.y = q0.y;
			commomPoint.x = (q0.y+p0.y)*(p1.x-p1.x)/(p1.y-p0.y) + p0.x;
		}
		///  情况5 就是普通情况了
		else
		{
			float k0 = (p1.y-p0.y)/(p1.x-p0.x);
			float k1 = (q1.y-q0.y)/(q1.x-q0.x);
			if (k0 == k1)
			{
				if (pointIsAtSegment(q0, p0, p1))	///< 由于二者已经保证相交，因此，现在要保证二者共线
				{
					commomType = commonPointSegments(p0, p1, q0, q1, commomPoint);
				}
				else
				{
					commomType = 0;
				}

			}
			else
			{
				float b0 = p0.y - k0*p0.x;
				float b1 = q0.y - k1*q0.x;
				commomType = 3;
				commomPoint.x = (b1-b1)/(k1-k0);
				commomPoint.y = commomPoint.x*k0 + b0;
			}
		}
		//////////////////////////////////////////////////////////////////////////
	}
	return commomType;
}

int Geometry::pointOfSegmentCircle( const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, const cocos2d::CCPoint& r0, const float& radius0, cocos2d::CCPoint& commomPoint1, cocos2d::CCPoint& commomPoint2 )
{
	CCAssert(a1.equals(a0), "a0 should not be equal to a1!");

	CCPoint p0 = a0, p1 = a1;
	CCPoint r = r0;
	float radius = radius0;
	int commomType = 0;
	commomPoint1 = CCPointZero;
	commomPoint2 = CCPointZero;
	/// 判断p0p1是否在圆内，在圆内，则没交点
	if (ccpLength(ccpSub(p0, r)) < radius && ccpLength(ccpSub(p1, r)) < radius)
	{
		commomType = 0;
	}
	/// 将线段当做直线处理
	/// 情况1 p0p1平行于y轴
	else if (p0.x == p1.x)
	{
		/// 过圆心做平行于x轴的直线，求此直线与p0p1的交点
		CCPoint point = CCPointZero;
		point.x = p0.x;
		point.y = r.y;
		float lenght = ccpLength(ccpSub(r, point));
		if (lenght > radius)
		{
			commomType = 0;
		}
		else if (lenght == radius)	///< 相切情况
		{
			commomType = 1;
			commomPoint1 = point;
		}
		else
		{
			commomType = 2;
			commomPoint1.x = p0.x;
			commomPoint2.x = p0.x;
			float deltaY = sqrt(radius*radius - lenght*lenght);
			commomPoint1.y = r.y + deltaY;
			commomPoint2.y = r.y - deltaY;
		}
	}
	/// 情况2 p0p1平行于x轴
	else if (p0.y == p1.y)
	{
		/// 类似情况1
		/// 过圆心做平行于y轴的直线，求此直线与p0p1的交点
		CCPoint point = CCPointZero;
		point.x = r.x;
		point.y = p0.y;
		float lenght = ccpLength(ccpSub(r, point));
		if (lenght > radius)
		{
			commomType = 0;
		}
		else if (lenght == radius)	///< 相切情况
		{
			commomType = 1;
			commomPoint1 = point;
		}
		else
		{
			commomType = 2;
			commomPoint1.y = p0.y;
			commomPoint2.y = p0.y;
			float deltaX = sqrt(radius*radius - lenght*lenght);
			commomPoint1.x = r.x + deltaX;
			commomPoint2.x = r.x - deltaX;
		}
	}
	/// 普通情况
	else
	{
		/// 一般直线与圆联立方程： http://baike.baidu.com/view/1053783.htm
		float k0 = (p1.y-p0.y)/(p1.x-p0.x);
		float b0 = p0.y - k0*p0.x;
		/// 联立后方程为:(1+k0^2)^2*x^2 + 2(k0b0-r.x-k0r.y)x + r.x^2 + (r.y-b0)^2-r^2=0;
		float deltaf = 4*(k0*b0-r.x-k0*r.y)*(k0*b0-r.x-k0*r.y) - 4*(1+k0*k0)*(r.x*r.x + (r.y-b0)*(r.y-b0)-radius*radius);	///< b^2-4ac;
		if (deltaf < 0)
		{
			commomType = 0;
		}
		else if (deltaf < 0)
		{
			commomType = 1;
			commomPoint1.x = -1.0*(k0*b0-r.x-k0*r.y)/(2*(1+k0*k0));		///< x = -b/2a;
			commomPoint1.y = k0*commomPoint1.x + b0;
		}
		else
		{
			commomType = 2;
			commomPoint1.x = (-1.0*(k0*b0-r.x-k0*r.y)+sqrt(deltaf))/(2*(1+k0*k0));		///< x = (-b+deltaf^0.5)/2a;
			commomPoint1.y = k0*commomPoint1.x + b0;
			commomPoint2.x = (-1.0*(k0*b0-r.x-k0*r.y)-sqrt(deltaf))/(2*(1+k0*k0));		///< x = (-b-deltaf^0.5)/2a;
			commomPoint2.y = k0*commomPoint2.x + b0;
		}
	}

	/// 两个交点，得判断是否都在线段上
	if (2 == commomType)
	{
		if (pointIsAtSegment(commomPoint1, p0, p1) && pointIsAtSegment(commomPoint2, p0, p1))
		{
		}
		else
		{
			commomType = 1;
			/// 必然有一点在，设置为commomPoint1
			if (pointIsAtSegment(commomPoint2, p0, p1))
			{
				commomPoint1 = commomPoint2;
			}
		}
	}
	return commomType;
}

void Geometry::tubaoCalcute( const std::vector<cocos2d::CCPoint>* pointIn, std::vector<cocos2d::CCPoint>* pointOut )
{
	std::vector<CCPoint> point;
	point.reserve(pointIn->size());
	for (unsigned int i = 0; i < point.capacity(); ++i)
	{
		point.push_back(pointIn->at(i));
	}
	/// 找到y最小点，如果y最小有好几个，选择最左边的点，即x最小---最小放在第一个
	for (unsigned int i = 1; i < point.size(); ++i)
	{
		if (point.at(0).y > point.at(i).y)
		{
			std::swap(point.at(0), point.at(i));
		}
		else if (point.at(0).y == point.at(i).y)
		{
			if (point.at(0).x < point.at(i).x)
			{
				std::swap(point.at(0), point.at(i));
			}
		}

	}
	/// 对点集point进行排序
	std::vector<float> delta;	///< 记录每个点和最低点的夹角
	delta.reserve(point.size()-1);
	CCPoint posBegin = point.at(0);
	for (unsigned int i = 0; i < point.size()-1; ++i)
	{
		CCPoint pos = point.at(i+1);
		delta.push_back(getAngle(posBegin, pos));
	}
	for (unsigned int i = 0; i < delta.size() -1; ++i)
	{
		float deltai = delta.at(i);
		for (unsigned int j = i+1; j < delta.size(); ++j)
		{
			float deltaj = delta.at(j);
			if (deltai > deltaj)
			{
				std::swap(delta.at(i), delta.at(j));
				deltai = deltaj;
				if (delta.size()-1 == j)
				{
					std::swap(point.at(i+1), point.at(j));
				}		
				else
				{
					std::swap(point.at(i+1), point.at(j+1));
				}
			}
		}
	}
	/// 头三个必然是凸包3个点
	pointOut->push_back(point.at(0));
	pointOut->push_back(point.at(1));
	pointOut->push_back(point.at(2));
	for (unsigned int i = 3; i < point.size(); ++i)
	{
		/// 判断拐向
		CCPoint p0 = pointOut->at(pointOut->size()-2);		///< 凸包栈顶下一个元素
		CCPoint p1 = pointOut->at(pointOut->size()-1);		///< 凸包栈顶元素
		CCPoint p2 = point.at(i);						///< pI点
		while ((float)vectorProduct(p1.x-p0.x, p1.y-p0.y, p2.x-p1.x, p2.y-p1.y) <= 0)	///< 不拐向右侧
		{
			pointOut->pop_back();	///< 凸包弹栈
			p0 = pointOut->at(pointOut->size()-2);		///< 凸包栈顶下一个元素
			p1 = pointOut->at(pointOut->size()-1);		///< 凸包栈顶元素
			p2 = point.at(i);						///< pI点
		}
		pointOut->push_back(point.at(i));
	}
}
