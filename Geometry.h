/**
*@file Geometry.h
*@brief һЩ���ü����㷨 http://dev.gameres.com/Program/Abstract/Geometry.htm
*@author DionysosLai��email: 906391500@qq.com
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
///@brief ���
///@param[in] (x1, y1)---��1�� (x2, y2)---��2
///@return ���ֵ
	static double vectorProduct(const double& x1, const double& y1, const double& x2, const double& y2);  // ����ʽ 

///@brief �ж����ߵĹ���
///���߶�P1Ϊ�ο�����
///@param[in] p1--�߶�1��p2--�߶�2
///@pre p0Ϊ�����˵㣬���ж�p2p0 �� p1p0�ķ���
///@return 1--��� 0 -- ���� -1---�ұ�
	static int polyLineDerection(const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p2);

///@brief �ж��߶���ֱ���Ƿ��ཻ
///@param[in] p0,p1--�߶������˵㣬 q0,q1--ֱ��������ֵ
///@return true---���߶��ϣ� false---�����߶���
	static bool straightLineIsIntersect(const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& q0, const cocos2d::CCPoint& q1);

///@brief �жϵ��Ƿ��ھ�����
///@param[in|out] p0--�� r1--�������� width--���ο� heigth---���θ�
///@return true--�ھ����ڣ��������ھ��α��ϣ��� false---�ھ�����
	static bool pointIsInRect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& r1, const float& width, const float& heigth );

///@brief �жϵ��Ƿ��ھ�����
///@param[in] p0--�㣬p1,p2--���ε�һ���Խ��߶˵�
///@return true--�ھ����ڣ��������ھ��α��ϣ��� false---�ھ�����
	static bool pointIsInRect( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2 );

///@brief �жϵ��Ƿ����߶���
///@param[in] p0--�㣬p1,p2--�߶������˵�
///���ΪQ���߶�ΪP1P2 ���жϵ�Q�ڸ��߶��ϵ������ǣ�( Q - P1 ) �� ( P2 - P1 ) = 0 
///�� Q ���� P1��P2Ϊ�ԽǶ���ľ����ڡ�ǰ�߱�֤Q����ֱ��P1P2�ϣ������Ǳ�֤Q�㲻��
///�߶�P1P2���ӳ��߻����ӳ�����
///@return true---���߶��ϣ� false---�����߶���
	static bool pointIsAtSegment(const cocos2d::CCPoint& p0, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2);

///@brief �㵽�߶����һ����
///@param[in] p0--Ҫ�жϵ㣬 s0��s1--�߶������˵�
///@return crosspoint---����ĵ�
	static cocos2d::CCPoint nearestPointToSegmentLine( const cocos2d::CCPoint& p0, const cocos2d::CCPoint& s0, const cocos2d::CCPoint& s1 );

///@brief �ж����������Ƿ��ཻ
///@param[in] aa,bb--����1һ���Խ��߶˵㣬 cc,dd--����2һ���Խ��߶˵�
///@return true---�ཻ�� false---���ཻ
	static bool isRectsInterSect(const cocos2d::CCPoint& aa, const cocos2d::CCPoint& bb, const cocos2d::CCPoint& cc, const cocos2d::CCPoint& dd);

///@brief �жϵ��Ƿ��ڶ���Σ��������ڱ��ϣ�
///@param[in] p0--Ҫ�жϵ㣬 poly--����ε㼯�ϣ� numberOfPoints--����ε�����
///@return true---���ڶ�����ڣ� false---�㲻�ڶ������
	static bool pointIsInPolygon(const cocos2d::CCPoint& p0, std::vector<cocos2d::CCPoint>* point);

///@brief �жϵ��Ƿ���԰��
///@param[in] p0--�㣬 r0---Բ�� radius---�뾶
///@return true---����԰�ڣ� false--���ڶ������
	static bool pointInInCircle(const cocos2d::CCPoint& p0, const cocos2d::CCPoint& r0, const float& radius);

///@brief �ж��߶����߶��Ƿ��ཻ
///@param[in] aa,bb--�߶�1�����˵㣬 cc,dd--�߶�1�����˵�
///@return true---���߶��ϣ� false---�����߶���
	static bool segmentLineIsIntersect(const cocos2d::CCPoint& aa, const cocos2d::CCPoint& bb, const cocos2d::CCPoint& cc, const cocos2d::CCPoint& dd);

///@brief ��ȡ�����ĽǶ�ֵ
///@param[in] posBegin---���  posEnd---�յ�
///@return �Ƕ�ֵ
	static float getAngle( const cocos2d::CCPoint& posBegin, const cocos2d::CCPoint& posEnd );

///@brief �ж��߶�Բ�Ƿ��ཻ
///@param[in] r1---Բ�� radius--�뾶 ��p1, p2��---�߶������˵�ֵ
///@return true---��Բ���ڣ� false---����Բ����
	static bool isCircleLineCollision(const cocos2d::CCPoint& r1, const float& radius, const cocos2d::CCPoint& p1, const cocos2d::CCPoint& p2);

///@brief �ж�����Բ�Ƿ���ײ
///@param[in] r1--Բ1Բ�� radius1---Բ1�뾶 r2--Բ2Բ�� radius2---Բ2�뾶
///@return ture---�ཻ false---���ཻ
	static bool isCircleCollision(const cocos2d::CCPoint& r1, const float& radius1, const cocos2d::CCPoint& r2, const float& radius2);

///@brief �ж����������Ƿ��ཻ
///@param[in] a0, aWidth, bHeight---����1�е�Ϳ�ȡ��߶� 
///@return true---�ཻ  false---���ཻ
	static bool isRectsCollision(const cocos2d::CCPoint& a0, const float& aWidth, const float& aHeight, 
		const cocos2d::CCPoint& b0, const float& bWidth, const float& bHeight);

///@brief �ж��߶��Ƿ��ڶ������
///@param[in] a0 a1 �߶������˵� point----����θ����˵�ֵ
///@return true--�ڶ������ false---���ڶ������
	static bool isSegmentLineInPoly(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, std::vector<cocos2d::CCPoint>* point);

///@brief �������������߶εĽ���
///@param[in] a0a1---�߶�1 b0b1---�߶�2
///@param[out] commomPoint----����
///@pre �����߶α����ǹ���---���ǵ�
///@return 0----�����߶β����� 1----�����߶ι������޽��� 2----�����߶ι���������������  3----�����߶ι�����ֻ��һ������
	static int commonPointSegments(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, 
		const cocos2d::CCPoint& b0, const cocos2d::CCPoint& b1, cocos2d::CCPoint& commomPoint);

///@brief ���������߶εĽ���
///@param[in] a0a1---�߶�1 b0b1---�߶�2
///@param[out] commomPoint----����
///@pre �����߶Σ����ǵ�
///@return 0----�����߶�û�н��� 1----�����߶ι��������������� 2----�����߶ι�����ֻ��һ������ 3----�����߶β����ߣ���ֻ��һ������
	static int pointOfSegments(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, 
		const cocos2d::CCPoint& b0, const cocos2d::CCPoint& b1, cocos2d::CCPoint& commomPoint);

///@brief �����߶���Բ�Ľ���
///@param[in] a0a1---�߶�1 r0---Բ�ģ� radius---�뾶 
///@param[out] commomPoint1----����1 commomPoint2----����2
///@pre �����߶Σ����ǵ�
///@return 0----û�н��� 1----1������ 2----2������
	static int pointOfSegmentCircle(const cocos2d::CCPoint& a0, const cocos2d::CCPoint& a1, const cocos2d::CCPoint& r0, 
		const float& radius, cocos2d::CCPoint& commomPoint1, cocos2d::CCPoint& commomPoint2);

///@brief ����㼯��͹��---ʹ�ø����㣨Graham��ɨ�跨 
///@param[in] pointIn---�㼯  pointOut---͹���㼯
///@pre ����2�� http://zh.wikipedia.org/wiki/%E5%87%B8%E5%8C%85#.E8.91.9B.E7.AB.8B.E6.81.92.EF.BC.88Graham.EF.BC.89.E6.89.AB.E6.8F.8F.E6.B3.95
/// http://blog.csdn.net/suwei19870312/article/details/5422818 http://dev.gameres.com/Program/Abstract/Geometry.htm
///@return 
	static void tubaoCalcute(const std::vector<cocos2d::CCPoint>* pointIn, std::vector<cocos2d::CCPoint>* pointOut);
};
#endif	///<(__GEOMETRY_H__)