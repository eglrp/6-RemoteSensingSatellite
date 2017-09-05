// GeometryProcess.cpp : Defines the entry point for the console application.
//
#include <time.h>
#include <vector>
#include<iostream>
#include <iomanip>
#include "WorkFlow_ZY3.h"
#include <windows.h>
#pragma comment(lib, "winmm.lib")
using namespace std;

int main(int argc, char* argv[])
{	
	string argv2 = "C:\\Users\\wcsgz\\Documents\\2-CProject\\6-����ģ��\\ExtDlls\\EOP00.txt";
	//////////////////////////////////////////////////////////////////////////
	//���ܣ�����ΪС���������������֤
	//���ڣ�2017.08.14
	/////////////////////////////////////////////////////////////////////////
	WorkFlow_ZY3 *pflow = new WorkFlow_ZY3();
	pflow->getEOP(argv2);	
	//pflow->LittleArrayCamera(argv[1]);//����ģ���Լ�����rpc	
	//pflow->Image_registration_rpc(argv[1]);//����siftGPU��ƥ����Ƶ�
	//pflow->CalcRealMatchPoint(argv[1]);//���ݲ�������ģ������һ���������Ƶ�,Ŀǰ�������޲���ϼ�ϵͳ��
	//pflow->CalcRealAttitude(argv[1]);//���ݵ�һ֡��ƥ��������̬
	//pflow->CalcOmegaKalman(argv[1]);//����ƥ�����ÿ������˲�������̬
	pflow->CalcRealAttitude_sparse(argv[1]);


	//pflow->NADCamera(argv[1]);
	//pflow->ModelVerify();


	//////////////////////////////////////////////////////////////////////////
	//���ܣ�����Ϊ����02�Ƕ���Ͷ�λģ��
	//���ڣ�2017.04.25
	/////////////////////////////////////////////////////////////////////////
	//WorkFlow_ZY3 *pflow = new WorkFlow_ZY3();
	////�ⷽλԪ�ؼ�У
	////pflow->CalibrationModel2(argv[1], argv2);
	////����ģ�͹���
	//pflow->GenerateRigorousModel(argv[1],argv2);
	//pflow->AccuracyVerify(argv[1]);
	//double lat, lon, h = 0, x=1, y=1;
	//while (x!=0&&y!=0)
	//{
	//	cout << "������Sample��Line(��������0�������˳�)" << endl;
	//	cin >> x >> y;
	//	pflow->pModel->FromXY2LatLon(y, x, h, lat, lon);
	//	lat = lat * 180 / PI;
	//	lon = lon * 180 / PI;
	//	cout << "��γ�ȷֱ���" << endl;
	//	cout << setiosflags(ios::fixed);//������仰�����Ƶľ���С���ľ�����
	//	cout << setprecision(10) << lat << "," << setprecision(10) << lon << endl;
	//}
	

	//////////////////////////////////////////////////////////////////////////
	//���ܣ�����Ϊ����01�Ƕ���Ͷ�λģ��
	//���ڣ�2016.12.21
	/////////////////////////////////////////////////////////////////////////
	//WorkFlow_ZY3 *pflow = new WorkFlow_ZY3();
	//�ⷽλԪ�ؼ�У
	//pflow->CalibrationModel(argv[1],argv[2]);
	//����ģ�͹���
	//pflow->GenerateRigorousModel(argv[1],argv[2]);
	//double lat, lon, h, x, y, lat1, lon1, lat2, lon2, lat3, lon3;
	//h = 0;
	//FromXY2LatLon(x,y,h,lat,lon),x��line��y��sample
	//pflow->pModel->FromXY2LatLon( 7183, 3256,h, lat, lon);
	//lat1 = lat*180/PI;
	//lon1 = lon*180/PI;
	//pflow->pModel->FromXY2LatLon(14375, 1793,  h, lat, lon);
	//lat2 = lat*180/PI;
	//lon2 = lon*180/PI;
	//pflow->pModel->FromXY2LatLon(15279, 16913, h, lat, lon);
	//lat3 = lat*180/PI;
	//lon3 = lon*180/PI;

	//����381��
	//h = 326.5;
	//pflow->pModel->FromXY2LatLon(7579, 5267, h, lat, lon);
	//lat1 = lat*180/PI;
	//lon1 = lon*180/PI;
	//pflow->pModel->FromXY2LatLon(1571, 11433,  h, lat, lon);
	//lat2 = lat*180/PI;
	//lon2 = lon*180/PI;
	//pflow->pModel->FromXY2LatLon(2741, 22080, h, lat, lon);
	//lat3 = lat*180/PI;
	//lon3 = lon*180/PI;
	//pflow->pModel->FromXY2LatLon(5383, 14532, h, lat, lon);
	//lat1 = lat*180/PI;
	//lon1 = lon*180/PI;	
	//h = 413.3;
	//pflow->pModel->FromXY2LatLon(1683, 5562,  h, lat, lon);
	//lat2 = lat*180/PI;
	//lon2 = lon*180/PI;
	////pflow->pModel->FromLatLon2XY(lat, lon, h, x, y);
	if(pflow!=NULL)		delete []pflow;		pflow = NULL;
	PlaySound(TEXT("C:\\WINDOWS\\Media\\Alarm01.wav"), NULL, NULL);
	return 0;
}
