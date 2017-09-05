#pragma once
#include "geomodelline.h"
#include "GeoModelArray.h"
#include "GeoBase.h"

class GeoCalibration :
	public GeoModelLine
{
public:
	GeoCalibration(void);
	~GeoCalibration(void);
	GeoBase m_base;
	// 外方位元素检校
	void ExtOrientCali(GeoOrbit *orbit, GeoAttitude *att, GeoTime *time, GeoCamera *cam, StrModelParamInput input,vector<StrGCP>ZY3_GCP);
	//计算前后两帧偏置补偿矩阵
	bool calcOffsetMatrix(GeoModelArray* pModel, StrGCP* pGCP, int numGCP, OffsetAngle &angle);
	//根据修正后Ru计算角速度
	void CalcRealOmega(GeoModelArray *pModel, OffsetAngle Ru, structEOP* Eop, Gyro & omega);
	//计算残差
	void calcRMS(GeoModelArray* pModel, string workpath, StrGCP *pGCP, int numGCP);
	//根据真实控制点计算残差
	void calcGCPerr(GeoModelArray* pModel, string strImg, string out, vector<RMS>&acc, bool isPlus1);
};

