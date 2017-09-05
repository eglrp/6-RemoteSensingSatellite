#include "GeoCalibration.h"


GeoCalibration::GeoCalibration(void)
{
}


GeoCalibration::~GeoCalibration(void)
{
}

/////////////////////////////////////////////
//功能：检校线性模型
//输入：姿态、轨道、行时、相机参数、模型参数、控制点vector
//输出：更新Ru参数
//注意：目前未知参数为三个角度
//日期：2016.12.21
////////////////////////////////////////////
void GeoCalibration::ExtOrientCali(GeoOrbit *orbit, GeoAttitude *att, GeoTime *time, GeoCamera *cam, StrModelParamInput input,vector<StrGCP>ZY3_GCP)
{
	m_orbit = orbit;
	m_att = att;
	m_time = time;
	m_cam = cam;
	m_input = input;
	// 获取模型影像长度和宽度
	m_xnum = m_time->get_num();
	m_ynum = m_cam->get_ynum();
	// 获取模型所定位的基准椭球参数
	m_datum = m_orbit->get_datum();
	// 计算姿轨多项式模型以及整数位的姿轨以及常见的安装
	m_att->get_ROff(m_body2att.R);//相对于对Ru赋初值单位阵
	ReCalPosAndAtt();

	//下面开始定标过程
	double xpixel, ypixel, X, Y, Z, Xs, Ys, Zs, t, Xt[3],Xr[3];
	double Xbar, Ybar, Zbar;
	int i,j, num = ZY3_GCP.size();
	//构建Ru并设初值为0；
	double phi=0,omg=0,kap=0,sinphi,sinomg,sinkap,cosphi,cosomg,coskap;	
	//设置迭代判断初值
	double iter_phi = PI, iter_omg = PI, iter_kap = PI, iter_count = 0;
	//有关法化计算的数组
	double ATA[9], ATL[3], L;

	while (ZY3_GCP.size())
	{
		memset(ATA, 0, sizeof(double)*9);	memset(ATL, 0, sizeof(double)*3);
		sinphi = sin(phi); sinomg = sin(omg); sinkap = sin(kap);
		cosphi = cos(phi); cosomg = cos(omg); coskap = cos(kap);		
		m_base.Eulor2Matrix(phi,omg,kap,213,Ru);	
		for (i=0; i<num; i++)
		{
			double Orb[3], Att[9], Inner[3];
			xpixel = ZY3_GCP[i].y;
			ypixel = ZY3_GCP[i].x;
			//根据像素位置得到轨道姿态和内方位元素
			GetOrbitAttitudeInnerBaseXY(xpixel, ypixel, Orb, Att, Inner);
			m_base.Geograph2Rect(m_datum,ZY3_GCP[i].lat/180*PI,ZY3_GCP[i].lon/180*PI,ZY3_GCP[i].h,X,Y,Z);
			//根据地面点位置和卫星位置计算摄影光线指向
			Xt[0] = X - Orb[0];	Xt[1] = Y - Orb[1];	Xt[2] = Z - Orb[2];
			m_base.invers_matrix(Att,3);
			m_base.Multi(Att,Xt,Xr,3,3,1);
			
			//组合为X_ Y_ Z_
			Xbar = Ru[0]*Xr[0] + Ru[3]*Xr[1] + Ru[6]*Xr[2];
			Ybar = Ru[1]*Xr[0] + Ru[4]*Xr[1] + Ru[7]*Xr[2];
			Zbar = Ru[2]*Xr[0] + Ru[5]*Xr[1] + Ru[8]*Xr[2];
	
			//旋转矩阵对各项的偏微分
			//对phi角的偏微分
			double partial_a1phi = -sinphi*coskap + cosphi*sinomg*sinkap;
			double partial_a2phi = 0;
			double partial_a3phi = -cosphi*coskap - sinphi*sinomg*sinkap;
			double partial_b1phi = sinphi*sinkap + cosphi*sinomg*coskap;
			double partial_b2phi = 0;
			double partial_b3phi = cosphi*sinkap - sinphi*sinomg*coskap;
			double partial_c1phi = cosphi*cosomg;
			double partial_c2phi = 0;
			double partial_c3phi = -sinphi*cosomg;
			//对omega角的偏微分
			double partial_a1omg = sinphi*cosomg*sinkap;
			double partial_a2omg = -sinomg*sinkap;
			double partial_a3omg = cosphi*cosomg*sinkap;
			double partial_b1omg = sinphi*cosomg*coskap;
			double partial_b2omg = -sinomg*coskap;
			double partial_b3omg = cosphi*cosomg*coskap;
			double partial_c1omg = -sinphi*sinomg;
			double partial_c2omg = -cosomg;
			double partial_c3omg = -cosphi*sinomg;
			//对kappa角的偏微分
			double partial_a1kap = -cosphi*sinkap + sinphi*sinomg*coskap;
			double partial_a2kap = cosomg*coskap;
			double partial_a3kap = sinphi*sinkap + cosphi*sinomg*coskap;
			double partial_b1kap = -cosphi*coskap - sinphi*sinomg*sinkap;
			double partial_b2kap = -cosomg*sinkap;
			double partial_b3kap = sinphi*coskap - cosphi*sinomg*sinkap;
			double partial_c1kap = 0;
			double partial_c2kap = 0;
			double partial_c3kap = 0;
			//求Xbar,Ybar,Zbar三者的偏导数
			double partial_Xbar_phi = Xr[0]*partial_a1phi + Xr[1]*partial_b1phi + Xr[2]*partial_c1phi;
			double partial_Xbar_omg = Xr[0]*partial_a1omg + Xr[1]*partial_b1omg + Xr[2]*partial_c1omg;
			double partial_Xbar_kap = Xr[0]*partial_a1kap + Xr[1]*partial_b1kap + Xr[2]*partial_c1kap;
			double partial_Ybar_phi = Xr[0]*partial_a2phi + Xr[1]*partial_b2phi + Xr[2]*partial_c2phi;
			double partial_Ybar_omg = Xr[0]*partial_a2omg + Xr[1]*partial_b2omg + Xr[2]*partial_c2omg;
			double partial_Ybar_kap = Xr[0]*partial_a2kap + Xr[1]*partial_b2kap + Xr[2]*partial_c2kap;
			double partial_Zbar_phi = Xr[0]*partial_a3phi + Xr[1]*partial_b3phi + Xr[2]*partial_c3phi;
			double partial_Zbar_omg = Xr[0]*partial_a3omg + Xr[1]*partial_b3omg + Xr[2]*partial_c3omg;
			double partial_Zbar_kap = Xr[0]*partial_a3kap + Xr[1]*partial_b3kap + Xr[2]*partial_c3kap;
			//方程Vx法化
			double Ax[3],Ay[3];
			Ax[0] = 1/Zbar*partial_Xbar_phi - Xbar/Zbar/Zbar*partial_Zbar_phi;
			Ax[1] = 1/Zbar*partial_Xbar_omg - Xbar/Zbar/Zbar*partial_Zbar_omg;
			Ax[2] = 1/Zbar*partial_Xbar_kap - Xbar/Zbar/Zbar*partial_Zbar_kap;
			L = Inner[0] - Xbar/Zbar;
			m_base.pNormal(Ax, 3, L, ATA, ATL, 1.0);
			//方程Vy法化
			Ay[0] = 1/Zbar*partial_Ybar_phi - Ybar/Zbar/Zbar*partial_Zbar_omg;
			Ay[1] = 1/Zbar*partial_Ybar_omg - Ybar/Zbar/Zbar*partial_Zbar_omg;
			Ay[2] = 1/Zbar*partial_Ybar_kap - Ybar/Zbar/Zbar*partial_Zbar_kap;			
			L = Inner[1] - Ybar/Zbar;
			m_base.pNormal(Ay, 3, L, ATA, ATL, 1.0);
		}
		//迭代求解
		m_base.solve33(ATA,ATL);
		if (abs(ATL[0]>iter_phi&&ATL[1]>iter_omg&&ATL[3]>iter_kap||iter_count>20))
			break;
		iter_count++;
		phi += ATL[0]; omg += ATL[1]; 	kap += ATL[2];
		m_base.Eulor2Matrix(phi,omg,kap,213,Ru);
		iter_phi = abs(ATL[0]), iter_omg = abs(ATL[1]), iter_kap = abs(ATL[2]);
	}
	// 开始计算精度
	double mean = 0;
	m_base.Eulor2Matrix(phi, omg, kap, 213, Ru);
	double Ru0[9];
	m_base.Multi(m_body2att.R, Ru, Ru0, 3, 3, 3);
	memcpy(m_body2att.R, Ru0, sizeof(double)*9);
	ComputerGlobalParam();
	double dxy, dx, dy;
	for(int i=0; i<num; i++)
	{
		FromLatLon2XY(ZY3_GCP[i].lat/180*PI,ZY3_GCP[i].lon/180*PI, ZY3_GCP[i].h, dx, dy);
		dx = ZY3_GCP[i].y - dx;
		dy = ZY3_GCP[i].x - dy;
		mean += (dx*dx + dy*dy);
	}
	mean = sqrt(mean/num);
	printf("偏置补偿后精度：%lf像素\n", mean);	
}

/////////////////////////////////////////////
//功能：检校线性模型2
//输入：面阵几何模型，控制点
//输出：更新Ru参数
//注意：
//作者：JYH
//日期：2017.08.18
////////////////////////////////////////////
bool GeoCalibration::calcOffsetMatrix(GeoModelArray* pModel, StrGCP* pGCP, int numGCP, OffsetAngle &angle)
{
	//int UnknowNum = 3;
	double A_T_P_A[9], A_T_P_L[3];
	int i;

	angle.RuKappa = angle.RuOmega = angle.RuPhi = 0;

	double pParams[3];
	double RCam2wgs84[9], Rwgs842cam[9], Ru[9], XYZ[3], Rinstall[9];
	int inter = 0;

	pModel->GetCam2WGS84(RCam2wgs84);
	memcpy(Rwgs842cam, RCam2wgs84, sizeof(double) * 9);
	m_base.invers_matrix(Rwgs842cam, 3);
	pModel->GetCamPos(XYZ);
	pModel->GetCamInstall(Rinstall);

	do
	{
		memset(A_T_P_A, 0, sizeof(double) * 9);
		memset(A_T_P_L, 0, sizeof(double) * 3);
		//memset(pX,0,sizeof(double)*UnknowNum);	

		double x, y, X, Y, Z, XX, YY, ZZ;
		m_base.rot(angle.RuPhi, angle.RuOmega, angle.RuKappa, Ru);
		for (int i = 0; i < numGCP; i++)
		{
			memset(pParams, 0, sizeof(double) * 3);
			x = pGCP[i].x;
			y = pGCP[i].y;
			StrDATUM WGS84;
			m_base.Geograph2Rect(WGS84, pGCP[i].lat, pGCP[i].lon, pGCP[i].h, X, Y, Z);

			double POS[3] = { X - XYZ[0],Y - XYZ[1],Z - XYZ[2] };
			double POSCam[3];
			m_base.Multi(Rwgs842cam, POS, POSCam, 3, 3, 1);
			XX = Ru[0] * POSCam[0] + Ru[3] * POSCam[1] + Ru[6] * POSCam[2];
			YY = Ru[1] * POSCam[0] + Ru[4] * POSCam[1] + Ru[7] * POSCam[2];
			ZZ = Ru[2] * POSCam[0] + Ru[5] * POSCam[1] + Ru[8] * POSCam[2];
			pParams[0] = ((Ru[5] * YY - Ru[4] * ZZ)*ZZ - (Ru[4] * XX - Ru[3] * YY)*XX) / (ZZ*ZZ);
			pParams[1] = (sin(angle.RuKappa)*ZZ*ZZ + (sin(angle.RuKappa)*XX +
				cos(angle.RuKappa)*YY)*XX) / (ZZ*ZZ);
			pParams[2] = (ZZ*YY) / (ZZ*ZZ);

			double L, phiX, phiY;

			pModel->GetInner(x, y, phiX, phiY);

			double tv1[3] = { phiX,phiY,1 };
			double tv2[3];
			m_base.Multi(Rinstall, tv1, tv2, 3, 3, 1);

			phiX = tv2[0] / tv2[2];
			phiY = tv2[1] / tv2[2];

			L = phiX - XX / ZZ;
			m_base.pNormal(pParams, 3, L, A_T_P_A, A_T_P_L, 1.);

			
			memset(pParams, 0, sizeof(double) * 3);
			pParams[0] = ((Ru[3] * ZZ - Ru[5] * XX)*ZZ - (Ru[4] * XX - Ru[3] * YY)*YY) / (ZZ*ZZ);
			pParams[1] = (cos(angle.RuKappa)*ZZ*ZZ + (sin(angle.RuKappa)*XX +
				cos(angle.RuKappa)*YY)*YY) / (ZZ*ZZ);
			pParams[2] = -(ZZ*XX) / (ZZ*ZZ);
			L = phiY - YY / ZZ;
			m_base.pNormal(pParams, 3, L, A_T_P_A, A_T_P_L, 1.);
		}

		m_base.Gauss(A_T_P_A, A_T_P_L, 3);
		angle.RuPhi += A_T_P_L[0];
		angle.RuOmega += A_T_P_L[1];
		angle.RuKappa += A_T_P_L[2];
		//GaussExt(A_T_P_A,A_T_P_L,pX,UnknowNum);
		//for(int j = 0;j<UnknowNum;j++)
		//	pUnknow[j] += pX[j];
		inter++;
		if (inter > 10)
		{
			printf("=>out of  convergence for offset matrix!\n");
			break;
		}

	} while (!m_base.isExitlter(A_T_P_L, 1e-9, 3));

	return true;
}

/////////////////////////////////////////////
//功能：计算两景影像间的真实角速度
//输入：面阵几何模型，控制点
//输出：精度文件
//注意：
//作者：JYH
//日期：2017.08.30
////////////////////////////////////////////
void GeoCalibration::CalcRealOmega(GeoModelArray *pModel,OffsetAngle Ru, structEOP* Eop,  Gyro &omega)
{
	double Rbody2wgs84L[9], Rbody2wgs84R[9], Rbody2wgs84Rt[9], dRu[9], Res[9];

	Attitude attL = pModel[0].GetQuatCam2wgs84();
	m_base.Quat2Matrix(attL.Q1, attL.Q2, attL.Q3, attL.Q0, Rbody2wgs84L);
	m_base.invers_matrix(Rbody2wgs84L, 3);

	Attitude attR = pModel[1].GetQuatCam2wgs84();
	m_base.Quat2Matrix(attR.Q1, attR.Q2, attR.Q3, attR.Q0, Rbody2wgs84Rt);
	m_base.rot(Ru.RuPhi, Ru.RuOmega, Ru.RuKappa, dRu);
	m_base.Multi(Rbody2wgs84Rt, dRu, Rbody2wgs84R, 3, 3, 3);
	//m_base.invers_matrix(Rbody2wgs84R, 3);//Cbj

	double dt = attR.UTC - attL.UTC;
	double om1[3];
	m_base.Multi(Rbody2wgs84R, Rbody2wgs84L, Res, 3, 3, 3);//先转左边，再转右边
	omega.UT = attL.UTC;
	omega.wx= (Res[5] - Res[7]) / 2 / dt;
	omega.wy = (Res[6] - Res[2]) / 2 / dt;
	omega.wz = (Res[1] - Res[3]) / 2 / dt;

	string out = "D:\\2_ImageData\\0.1\\Point1\\LAC\\2.影像文件\\omega.txt";
	FILE *fp = fopen(out.c_str(), "a");
	fprintf(fp, "%.9f\t%.9f\t%.9f\n", omega.wx, omega.wy, omega.wz);
	fclose(fp);
}

/////////////////////////////////////////////
//功能：残差计算
//输入：面阵几何模型，控制点
//输出：精度文件
//注意：
//作者：JYH
//日期：2017.08.22
////////////////////////////////////////////
void GeoCalibration::calcRMS(GeoModelArray* pModel,string out,StrGCP *pGCP,int numGCP)
{
	FILE* fp = fopen(out.c_str(), "w");
	fprintf(fp, "%d\n", numGCP);
	double maxx, maxy, minx, miny, rmsx, rmsy, plane;
	rmsx = rmsy = 0.;
	maxx = maxy = 0;
	minx = miny = 999999999.;
	for (int i = 0; i < numGCP; i++)
	{
		double x, y, ex, ey;
		pModel->FromLatLon2XY(pGCP[i].lat, pGCP[i].lon, pGCP[i].h, x, y);
		ex = x - pGCP[i].x;
		ey = y - pGCP[i].y;

		minx = min(minx, fabs(ex));
		miny = min(miny, fabs(ey));

		maxx = max(maxx, fabs(ex));
		maxy = max(maxy, fabs(ey));

		rmsx += ex*ex / numGCP;
		rmsy += ey*ey / numGCP;

		fprintf(fp, "%04d\t%lf\t%lf\t0.\t0.\t%lf\t%lf\n", i, pGCP[i].x, pGCP[i].y, ex, ey);
	}
	plane = sqrt(rmsx + rmsy);
	rmsx = sqrt(rmsx);
	rmsy = sqrt(rmsy);
	fprintf(fp, "x: %lf\t%lf\t%lf\n", minx, maxx, rmsx);
	fprintf(fp, "y: %lf\t%lf\t%lf\n", miny, maxy, rmsy);
	fprintf(fp, "plane: %lf", plane);
	fclose(fp);
}

/////////////////////////////////////////////
//功能：真实控制点残差计算
//输入：面阵几何模型，控制点
//输出：精度比较文件
//注意：
//作者：GZC
//日期：2017.08.23
////////////////////////////////////////////
void GeoCalibration::calcGCPerr(GeoModelArray* pModel, string strImg,string out,vector<RMS>&acc,bool isPlus1)
{
	string tmp1 = strImg.substr(strImg.rfind('_') - 3, 3);
	string tmp2 = strImg.substr(0, strImg.rfind('_') - 4);
	int tmp3 = atoi(tmp1.c_str());// +1;
	if (isPlus1==1)
	{		tmp3 = atoi(tmp1.c_str()) +1;	}//正序递推
	else
	{		tmp3 = atoi(tmp1.c_str());	}//倒序递推
	
	char strgcp[512];
	sprintf(strgcp, "%s_%03d%s",tmp2.c_str(), tmp3, ".ctl");
	FILE* fp = fopen(strgcp, "r");
	string strRes = strImg.substr(0, strImg.rfind('_')) + out + "_rms.txt";
	FILE* fpres = fopen(strRes.c_str(), "w");
	int numGCP;
	fscanf(fp, "%d\n", &numGCP);
	double maxx, maxy, minx, miny, rmsx, rmsy, plane;
	rmsx = rmsy = 0.;
	maxx = maxy = 0;
	minx = miny = 999999999.;
	for (int i = 0; i < numGCP; i++)
	{
		double lx, ly, rx, ry, ex, ey, lat, lon, h;
		fscanf(fp, "%*d\t%lf\t%lf\t%lf\t%lf\t%lf\n",  &ly, &lx,&lat, &lon, &h);
		pModel->FromLatLon2XY(lat/180.*PI, lon / 180.*PI, h, rx, ry);
		ex = rx - lx;
		ey = ry - ly;

		minx = min(minx, fabs(ex));
		miny = min(miny, fabs(ey));

		maxx = max(maxx, fabs(ex));
		maxy = max(maxy, fabs(ey));

		rmsx += ex*ex / numGCP;
		rmsy += ey*ey / numGCP;

		fprintf(fpres, "%04d\t%lf\t%lf\t0.\t0.\t%lf\t%lf\n", i, lx, ly, ex, ey);
	}
	plane = sqrt(rmsx + rmsy);
	rmsx = sqrt(rmsx);
	rmsy = sqrt(rmsy);
	fprintf(fpres, "x: %lf\t%lf\t%lf\n", minx, maxx, rmsx);
	fprintf(fpres, "y: %lf\t%lf\t%lf\n", miny, maxy, rmsy);
	fprintf(fpres, "plane: %lf", plane);
	RMS tmpRms;
	tmpRms.rmsall = plane; tmpRms.rmsx = rmsx; tmpRms.rmsy = rmsy;
	acc.push_back(tmpRms);
	fclose(fp);
	fclose(fpres);
}