/**
*@file Geometry.cpp
*@brief һЩ���ü����㷨 http://dev.gameres.com/Program/Abstract/Geometry.htm
*@author DionysosLai��email: 906391500@qq.com
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
	return (x1*y2-x2*y1);	///< ��� ������ѧ����෴��
}

int Geometry::polyLineDerection( const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p2 )
{
	float vectorProductResult = 0.0f;
	int derection = 0;

	vectorProductResult = (float)vectorProduct(p1.x-p0.x, p1.y-p0.y, p2.x-p0.x, p2.y-p0.y);

	if (vectorProductResult < 1 && vectorProductResult > -1)
	{
		derection = 0;	///< ����
	}
	else if (vectorProductResult >= 1)
	{
		derection = -1;	///< p2p0��p1p0�����ֱ�
	}
	else if(vectorProductResult <= -1)
	{
		derection = 1;	///< p2p0��p1p0�����ֱ�
	}
	
	return derection;
}

bool Geometry::straightLineIsIntersect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& q0, const cocos2d::CCPoint& q1 )
{
	/// ���ж�q0��q1�ܷ����һ��ֱ��
	if (!q0.equals(q1))
	{
		/// ��q0 == p0 q1 == p1ʱ�����Ϊ0�� ֻ��Ҫ�ж��߶��Ƿ����ֱ�߼��ɡ�
		if (0 <= vectorProduct(p0.x-q0.x, p0.y-q0.y, q1.x-q0.x, q1.y-q0.y) * 
			vectorProduct(q1.x-q0.x, q1.y-q0.y, p1.x-q0.x, p1.y-q0.y))
		{
			/* CCLOG("�ཻ");*/
			return true;
		}

		/* CCLOG("���ཻ");*/
		return false;
	}

	/*	CCLOG("��q0��q1���ܹ���һ��ֱ��");*/
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
	/// ���ж��Ƿ������p1 p2Ϊ�Խ��ߵľ�����
	if (pointIsInRect(p0, p1, p2))
	{
		/// �ж�p1p0, p2p0�Ƿ���
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
	/// �ж��߶��Ƿ���һ����
	if (s0.equals(s1))
	{
		return s0;
	}

	/// ��ʼ����Ϊԭ�㣻
	CCPoint crossPoint = CCPointZero;
	do 
	{
		/// �ж��߶��Ƿ�ƽ����x��
		if (s0.y == s1.y)
		{
			crossPoint = ccp(p0.x, s0.y);
			break;
		}
		/// �ж��߶��Ƿ�ƽ����y��
		if (s0.x == s1.x)
		{
			crossPoint = ccp(s0.x, p0.y);
			break;
		}
		/// ����߶β��������������ֻ�ܲ���ֱ�߷��̷�ʽ�������
		float k = (s1.y - s0.y)/(s1.x - s0.x);		///< ���б��
		/// �߶�ֱ�߷��̣�	y = k* ( x - s0.x) + s0.y
		/// ���߷���Ϊ��	y = (-1/k) * (x - p0.x) + p0.y ��
		/// ������ֱ�߷��̽��
		float x = ( k*k * s0.x + k * (p0.y - s0.y ) + p0.x ) / ( k*k + 1);
		float y = k * ( x - s0.x) + s0.y;
		crossPoint = ccp(x, y);
		break;
	} while (0);

	/// �жϴ�ֱ�Ƿ����߶���
	if (pointIsAtSegment(crossPoint, s0, s1))
	{
		return crossPoint;
	}
	else
	{
		/// ���������������˵㵽����ľ��룬ѡ����봹��Ͻ��Ķ˵㷵�ء�
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
	CCPoint	centre1 = ccp((aa.x + bb.x)/2.0f, (aa.y + bb.y)/2.0f);	///< ����е�ֵ
	CCPoint	centre2 = ccp((cc.x + dd.x)/2.0f, (cc.y + dd.y)/2.0f);
	float	lengthX	= abs(centre1.x - centre2.x);	///< ��������������ĵľ��� 
	float	lengthY	= abs(centre1.y - centre2.y); 
	float	lengthRect1X	= abs(aa.x - bb.x);		///< ����������γ��Ϳ�
	float	lengthRect1Y	= abs(aa.y - bb.y);
	float	lengthRect2X	= abs(cc.x - dd.x);
	float	lengthRect2Y	= abs(cc.y - dd.y);

	/// �����ȥ1�ǵ�������õġ�
	return  (lengthX < (lengthRect1X + lengthRect2X)/2.0f-1 && lengthY < (lengthRect1Y + lengthRect2Y)/2.0f-1) ? true : false;
}

bool Geometry::pointIsInPolygon( const cocos2d::CCPoint& p0, std::vector<cocos2d::CCPoint>* point )
{
	unsigned int count  = 0;		///< �����������L�����εĽ�������
	cocos2d::CCSize	winsize = CCDirector::sharedDirector()->getWinSize();
	/// �ѵ�p0����������һ������L��
	CCPoint leftPoint = ccp(-100.f, p0.y);
	CCPoint rightPoint = p0;


	/// �ж�ÿ����
	unsigned int numberOfPoints = point->size();
	for (unsigned int i = 0; i < numberOfPoints; i++)
	{
		/// ���жϵ�p0�Ƿ��ڱ�s�ϣ�
		CCPoint s0 = point->at(i);
		CCPoint s1 = point->at((i+1)%(numberOfPoints));
		if (pointIsAtSegment(p0, s0, s1))
		{
			/*			CCLOG("Point is at the %dth line", i);*/

			return true;
		}

		/// �жϱ�s�Ƿ���ƽ���ߣ�
		if (s0.y != s1.y)
		{		
			do 
			{
				/// �жϱ�s���Ƿ��ж˵���L�� ͬʱ ���жϸõ��Ƿ��Ǳ�s������ϴ��һ����
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

				/// �����sû�ж˵���L�ϣ����ж�s��L�Ƿ��ཻ
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
	/// �������߶��γɵľ��β��غϣ�˵�������߶α�Ȼ���ཻ
	if (!isRectsInterSect(aa, bb, cc, dd))
	{
		return false;
	}

	/// ������߻������ ע��"="�������
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
	if (deltaPoint.x >= 0 && deltaPoint.y >= 0)	///< ��һ����
	{
		cocosAngle = cocosAngle;
	}
	if (deltaPoint.x < 0 && deltaPoint.y >= 0)	///< �ڶ�����
	{
		cocosAngle = 180.f + cocosAngle;
	}
	if (deltaPoint.x < 0 && deltaPoint.y < 0)	///< ��������
	{
		cocosAngle = 180.f + cocosAngle;
	}
	if (deltaPoint.x >= 0 && deltaPoint.y < 0)	///< ��������
	{
		cocosAngle = 360.f + cocosAngle;
	}	

	//	CCLOG("%f", cocosAngle);
	return cocosAngle;
}

bool Geometry::isCircleLineCollision( const cocos2d::CCPoint& r1, const float& radius, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2 )
{
	/// �ж��߶��Ƿ���һ����
	float length = 0.f;
	if (p1.equals(p2))
	{
		length = ccpDistance(r1, p1);
	}

	/// ��ʼ����Ϊ0��
	CCPoint crossPoint = CCPointZero;
	do 
	{
		/// �ж��߶��Ƿ�ƽ����x��
		if (p1.y == p2.y)
		{
			crossPoint = ccp(r1.x, p1.y);
			break;
		}
		/// �ж��߶��Ƿ�ƽ����y��
		if (p1.x == p2.x)
		{
			crossPoint = ccp(p1.x, r1.y);
			break;
		}
		/// ����߶β��������������ֻ�ܲ���ֱ�߷��̷�ʽ�������
		float k = (p2.y - p1.y)/(p2.x - p1.x);		///< ���б��
		/// �߶�ֱ�߷��̣�	y = k* ( x - s0.x) + s0.y
		/// ���߷���Ϊ��	y = (-1/k) * (x - p0.x) + p0.y ��
		/// ������ֱ�߷��̽��
		float x = ( k*k * p1.x + k * (r1.y - p1.y ) + r1.x ) / ( k*k + 1);
		float y = k * ( x - p1.x) + p1.y;
		crossPoint = ccp(x, y);

		/// �жϴ�ֱ�Ƿ����߶���

		if (pointIsAtSegment(crossPoint, p1, p2))
		{
			/*		return crossPoint;*/
		}
		else
		{
			/// ���������������˵㵽����ľ��룬ѡ����봹��Ͻ��Ķ˵㷵�ء�
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
	float	lengthX	= abs(a0.x - b0.x);	///< ��������������ĵľ��� 
	float	lengthY	= abs(a0.y - b0.y); 

	return  (lengthX < (aWidth + bWidth)/2.0f && lengthY < (aHeight + bHeight)/2.0f) ? true : false;
}

bool Geometry::isSegmentLineInPoly( const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, std::vector<cocos2d::CCPoint>* point )
{
		/*
	if �߶�PQ�Ķ˵㲻���ڶ������ 
		then return false;
	�㼯pointSet��ʼ��Ϊ��;
	for ����ε�ÿ����s
		do if �߶ε�ĳ���˵���s��
			then ���ö˵����pointSet;
		else if s��ĳ���˵����߶�PQ��
			then ���ö˵����pointSet;
		else if s���߶�PQ�ཻ // ��ʱ���Ѿ����Կ϶����ڽ���
			then return false;
		��pointSet�еĵ㰴��X-Y��������;
		for pointSet��ÿ�������ڵ� pointSet[i] , pointSet[ i+1]
		do if pointSet[i] , pointSet[ i+1] ���е㲻�ڶ������
			then return false;
		return true;
		*/
	if (!pointIsInPolygon(a0, point) || !pointIsInPolygon(a1, point))		///< �����ж��߶�a0a1
	{
		return false;
	}
	
	/// �ж�ÿ����
	unsigned int numberOfPoints = point->size();
	std::vector<CCPoint> pointSet;
	for (unsigned int i = 0; i < numberOfPoints-1; i++)
	{
		/// ���жϵ�p0�Ƿ��ڱ�s�ϣ�
		CCPoint s0 = point->at(i);
		CCPoint s1 = point->at((i+1)%(numberOfPoints));
		do 
		{
			/// �ж��߶�a0a1�˵��Ƿ���s0s1��
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
			/// �ж��߶�s0s1�Ƿ���a0a1��
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
			/// �ж�a0a1�Ƿ���s0s1�ཻ���ཻ---���Ȼ�ڽ�����Ȼ�ж��߶�a0a1���ڶ������
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

	/// ����poisntSet, ��С�������ȼ�Ϊx
	CCAssert(0 == numOfPointSet%2, "The pointSet's points num should be 2 times!");
	for (unsigned int i = 0; i < numOfPointSet-1; ++i)
	{
		/// ð������
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

	

	/// �ж����ڽ��㼯���е��Ƿ��ڶ������
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
	/// ���ж�2���߶��Ƿ���
	if (0 != polyLineDerection(p0, q0, q1) || 0 != polyLineDerection(p0, q0, q1))
	{
		commomType = 0;
	}
	else
	{
		/// ����2���߶αȽϳ���Ϊp0p1
		float l0 = ccpLength(ccpSub(p1, p0));
		float l1 = ccpLength(ccpSub(q0, q1));
		if (l0 < l1)
		{
			/// �߶�p0p1���߶�q0q1����
			CCPoint point;
			point = p0;
			p0 = q0;
			q0 = point;
			point = p1;
			p1 = q1;
			q1 = p1;
		}

		/// �ж��߶�q0q1�����Ƿ����߶�p0p1��
		int m = 0, n = 0;
		/// �ж�q0���߶�p0p1�����
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
		/// �ж�q1���߶�p0p1�����
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
	/// �����ж�2���߶��Ƿ��ཻ�����ཻ����Ȼû�н���
	if (segmentLineIsIntersect(p0, p1, q0, q1))
	{
		commomType = 0;
	}
	else
	{
		/// �����߶��ཻ���������߶ε���ֱ�ߴ���
		//////////////////////////////////////////////////////////////////////////
		/// ���1 �߶�p0p1ƽ����y��
		if (p0.x == p1.x)
		{
			/// ���1.1 �߶�q0q1ƽ����y��
			if (q0.x == q1.x)
			{
				/// �ж�p0p1��q0q1�Ƿ�ɹ���
				if (p0.x == q0.x)
				{
					commomType = commonPointSegments(p0, p1, q0, q1, commomPoint);
				}
				else
				{
					commomType = 0;
				}
			}
			/// ���1.2 ��q0q1��ƽ����Y�ᣬ�򽻵������Ϊp0�ĺ����꣬���뵽q0q1��ֱ�߷����п��Լ�������������ꣻ
			else
			{
				commomType = 3;
				commomPoint.x = p0.x;
				commomPoint.y = (q1.y-q0.y)/(q1.x-q0.x)*(p0.x-q0.x) + q0.y;
			}
		}
		/// ���2 p0��p1�����겻ͬ������q0��q1��������ͬ����q0q1ƽ����Y�ᣬ�򽻵������Ϊq0�ĺ����꣬���뵽p0p1��ֱ�߷����п��Լ�������������ꣻ
		else if (p0.x != p1.x && q0.x == q1.x)
		{
			commomType = 3;
			commomPoint.x = q0.x;
			commomPoint.y = (p1.y-p0.y)/(p1.x-p0.x)*(q0.x-p0.x) + p0.y;
		}
		/// ���3 ���p0��p1��������ͬ����p0p1ƽ����X��
		else if (p0.y == p1.y)
		{
			///< ���3.1  ��q0q1Ҳƽ����X��
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
			///< ���3.2 ��q0q1��ƽ����X�ᣬ�򽻵�������Ϊp0�������꣬���뵽q0q1��ֱ�߷����п��Լ������������ꣻ
			else
			{
				commomType = 3;
				commomPoint.y = p0.y;
				commomPoint.x = (p0.y+q0.y)*(q1.x-q0.x)/(q1.y-q0.y) + q0.x;
			}
		}
		/// ���4 ���p0��p1�����겻ͬ������q0��q1��������ͬ����q0q1ƽ����X�ᣬ�򽻵�������Ϊq0�������꣬���뵽p0p1��ֱ�߷����п��Լ������������ꣻ
		else if (p0.y != p0.y && q0.y == q1.y)
		{
			commomType = 3;
			commomPoint.y = q0.y;
			commomPoint.x = (q0.y+p0.y)*(p1.x-p1.x)/(p1.y-p0.y) + p0.x;
		}
		///  ���5 ������ͨ�����
		else
		{
			float k0 = (p1.y-p0.y)/(p1.x-p0.x);
			float k1 = (q1.y-q0.y)/(q1.x-q0.x);
			if (k0 == k1)
			{
				if (pointIsAtSegment(q0, p0, p1))	///< ���ڶ����Ѿ���֤�ཻ����ˣ�����Ҫ��֤���߹���
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
	/// �ж�p0p1�Ƿ���Բ�ڣ���Բ�ڣ���û����
	if (ccpLength(ccpSub(p0, r)) < radius && ccpLength(ccpSub(p1, r)) < radius)
	{
		commomType = 0;
	}
	/// ���߶ε���ֱ�ߴ���
	/// ���1 p0p1ƽ����y��
	else if (p0.x == p1.x)
	{
		/// ��Բ����ƽ����x���ֱ�ߣ����ֱ����p0p1�Ľ���
		CCPoint point = CCPointZero;
		point.x = p0.x;
		point.y = r.y;
		float lenght = ccpLength(ccpSub(r, point));
		if (lenght > radius)
		{
			commomType = 0;
		}
		else if (lenght == radius)	///< �������
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
	/// ���2 p0p1ƽ����x��
	else if (p0.y == p1.y)
	{
		/// �������1
		/// ��Բ����ƽ����y���ֱ�ߣ����ֱ����p0p1�Ľ���
		CCPoint point = CCPointZero;
		point.x = r.x;
		point.y = p0.y;
		float lenght = ccpLength(ccpSub(r, point));
		if (lenght > radius)
		{
			commomType = 0;
		}
		else if (lenght == radius)	///< �������
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
	/// ��ͨ���
	else
	{
		/// һ��ֱ����Բ�������̣� http://baike.baidu.com/view/1053783.htm
		float k0 = (p1.y-p0.y)/(p1.x-p0.x);
		float b0 = p0.y - k0*p0.x;
		/// �����󷽳�Ϊ:(1+k0^2)^2*x^2 + 2(k0b0-r.x-k0r.y)x + r.x^2 + (r.y-b0)^2-r^2=0;
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

	/// �������㣬���ж��Ƿ����߶���
	if (2 == commomType)
	{
		if (pointIsAtSegment(commomPoint1, p0, p1) && pointIsAtSegment(commomPoint2, p0, p1))
		{
		}
		else
		{
			commomType = 1;
			/// ��Ȼ��һ���ڣ�����ΪcommomPoint1
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
	/// �ҵ�y��С�㣬���y��С�кü�����ѡ������ߵĵ㣬��x��С---��С���ڵ�һ��
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
	/// �Ե㼯point��������
	std::vector<float> delta;	///< ��¼ÿ�������͵�ļн�
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
	/// ͷ������Ȼ��͹��3����
	pointOut->push_back(point.at(0));
	pointOut->push_back(point.at(1));
	pointOut->push_back(point.at(2));
	for (unsigned int i = 3; i < point.size(); ++i)
	{
		/// �жϹ���
		CCPoint p0 = pointOut->at(pointOut->size()-2);		///< ͹��ջ����һ��Ԫ��
		CCPoint p1 = pointOut->at(pointOut->size()-1);		///< ͹��ջ��Ԫ��
		CCPoint p2 = point.at(i);						///< pI��
		while ((float)vectorProduct(p1.x-p0.x, p1.y-p0.y, p2.x-p1.x, p2.y-p1.y) <= 0)	///< �������Ҳ�
		{
			pointOut->pop_back();	///< ͹����ջ
			p0 = pointOut->at(pointOut->size()-2);		///< ͹��ջ����һ��Ԫ��
			p1 = pointOut->at(pointOut->size()-1);		///< ͹��ջ��Ԫ��
			p2 = point.at(i);						///< pI��
		}
		pointOut->push_back(point.at(i));
	}
}
