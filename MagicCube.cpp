//双缓冲无屏闪稳定动画
//多线程实现动画和输入的分离
//没有利用任何3D引擎，纯利用线性代数知识暴力计算实现3D效果。
//1920*1080 缩放125%
#include<iostream>
#include<algorithm>
#include<stack>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <windows.h>
#include<string>
#include <thread>
using namespace std;

//Definition
#define findNextPara(p) while(strIn[p] && !(strIn[p++] == '-'))
#define findNumber(x) while(strIn[p] <= '9' && strIn[p] >= '0'){ x *= 10;x += strIn[p] - '0';p++; }
typedef HWND(WINAPI *PROCGETCONSOLEWINDOW)();
HWND cmd; HDC hdc; HDC mydc; HDC tmpdc; HBITMAP hBitmap1; HBITMAP hBitmap2;
HPEN hPen; HBRUSH hBrush[7]; HPEN blackPen; HBRUSH blackBrush;
const double pi = acos(-1);
const double perPi = pi / 180;
const double esp = 5.0;
const int baseX = 270;
const int baseY = 270;
int isTraval; int fixConstant = 3; int finishFlag;
const int surfaceIndex[6][4] =
{
	{0,1,2,3},
	{4,5,6,7},
	{0,4,7,3},
	{1,5,6,2},
	{0,1,5,4},
	{3,2,6,7}
};//F,B,L,R,U,D
const int basePoint1[8][3] =
{
	{ 1,-1, 1},
	{ 1, 1, 1},
	{ 1, 1,-1},
	{ 1,-1,-1},
	{-1,-1, 1},
	{-1, 1, 1},
	{-1, 1,-1},
	{-1,-1,-1}
};
const int basePoint2[8][3] =
{
	{ 1, 0, 1 },
	{ 1, 1, 1 },
	{ 1, 1, 0 },
	{ 1, 0, 0 },
	{ 0, 0, 1 },
	{ 0, 1, 1 },
	{ 0, 1, 0 },
	{ 0, 0, 0 }
};
int basecolor[8] =
{
	0x000000FF,
	0x00EE861C,
	0x0000EE00,
	0x0000EEEE,
	0x0000A5FF,
	0x00D670DA,
	0x00555555,
	0x00AAAAAA
};//red,blue,green,yellow,orange,purple,innercolor,hpen
struct xpoint
{
	double x, y, z;
	xpoint operator+(xpoint a) { return { a.x + x,a.y + y,a.z + z }; }
	xpoint operator-(xpoint a) { return { x - a.x,y - a.y,z - a.z }; }
	xpoint operator*(double a) { return { x*a,y*a,z*a }; }
	xpoint operator/(double a) { return { x/a,y/a,z/a }; }
	bool operator<(xpoint a) { return x < a.x; }
};
struct Rotation//alpha天顶角, beta方位角
{
	double alpha, beta;
};
//Initialization
void initializePenAndBrush()
{
	for (int i = 0; i<7; i++) hBrush[i] = CreateSolidBrush(basecolor[i]);
	hPen = CreatePen(PS_SOLID, 4, basecolor[7]);
}
void initializeConsole()
{
	//设置窗口最大化
	HWND hwnd = GetForegroundWindow();
	int cx = GetSystemMetrics(SM_CXSCREEN);
	int cy = GetSystemMetrics(SM_CYSCREEN);
	LONG WinStyle = GetWindowLong(hwnd, GWL_STYLE);
	SetWindowLong(hwnd, GWL_STYLE, (WinStyle | WS_POPUP | WS_MAXIMIZE));
	SetWindowPos(hwnd, HWND_TOP, 0, 0, cx, cy, 0);
}
void initializeGDI()
{
	//get handle of consolewindow
	HMODULE hKernel32 = GetModuleHandleA("kernel32");

	PROCGETCONSOLEWINDOW GetConsoleWindow = (PROCGETCONSOLEWINDOW)GetProcAddress(hKernel32, "GetConsoleWindow");
	cmd = GetConsoleWindow();
	hdc = GetDC(cmd);
	mydc = CreateCompatibleDC(hdc);
	tmpdc = CreateCompatibleDC(hdc);
	hBitmap1 = CreateCompatibleBitmap(hdc, baseX * 2, baseY * 2);
	hBitmap2 = CreateCompatibleBitmap(hdc, baseX * 2, baseY * 2);
	SelectObject(mydc, hBitmap1);
	SelectObject(tmpdc, hBitmap2);
	//initialize pen and brush

	initializePenAndBrush();
	blackPen = CreatePen(PS_SOLID, 1, 0x000C0C0C);
	blackBrush = CreateSolidBrush(0x000C0C0C);
}
class faceOfCube
{
public:
	xpoint originalPoint[4];
	xpoint transferPoint[4];//投影绝对坐标
	xpoint originalCenterPoint;
	xpoint transferCenterPoint;//投影中心点
	int cmpIndex;//旋转中心cube的index
	int cubeIndex;//面归属的方块
	int colorIndex;
	void show()
	{
		POINT tempFace[4];
		for (int i = 0; i < 4; i++) tempFace[i] = { baseX + (LONG)transferPoint[i].z, baseY - (LONG)transferPoint[i].y };
		SelectObject(mydc, hPen);
		SelectObject(mydc, hBrush[colorIndex]);
		Polygon(mydc, tempFace, 4);
	}//注意这里画在虚拟位图上,显示需要使用bitblt
		
};
class cube
{
private:
	xpoint originalPoint[8];
	xpoint transferPoint[8];//投影绝对坐标
	xpoint originalCenterPoint;
	xpoint transferCenterPoint;//投影中心点
	int facesIndex[6];
	faceOfCube *faces;
	//检验基坐标变换
	//{ 1,-1,3 }, { 3,1,0 }, { 1,1,-1 }, { 2, 0, 3 }->{ 2, -1, 3 }
	//{ 1,-1,0 }, { 2,1,3 }, { 3,1,2 }, { 5, 0, 7 }->{ 2, 3, -1 }
	//{ 1,-1,0 }, { 2,1,3 }, { 3,1,2 }, { -9, -8, -13 }->{ 3, -3, -2 }
	xpoint BaseCoordinateTransformation(xpoint u, xpoint v, xpoint w, xpoint m)
	{
		double a = u.x, b = u.y, c = u.z;
		double d = v.x, e = v.y, f = v.z;
		double p = w.x, q = w.y, r = w.z;
		double x = m.x, y = m.y, z = m.z;
		return{ (f*q*x - e * r*x - f * p*y + d * r*y + e * p*z - d * q*z) /
			(c*e*p - b * f*p - c * d*q + a * f*q + b * d*r - a * e*r),
			(c*q*x - b * r*x - c * p*y + a * r*y + b * p*z - a * q*z) /
			(b*f*p - c * e*p + c * d*q - a * f*q - b * d*r + a * e*r),
			(c*e*x - b * f*x - c * d*y + a * f*y + b * d*z - a * e*z) /
			(c*e*p - b * f*p - c * d*q + a * f*q + b * d*r - a * e*r) };
	}
	xpoint chaMultiplication(xpoint u, xpoint v)
	{
		double l = u.x, m = u.y, n = u.z;
		double o = v.x, p = v.y, q = v.z;
		return { m * q - n * p,n * o - l * q,l * p - m * o };
	}
	void unitize(xpoint &vector)  //单位化
	{
		double a = vector.x, b = vector.y, c = vector.z;
		double mo = sqrt(a * a + b * b + c * c);
		vector = { a / mo,b / mo,c / mo };
	}
	void print(xpoint v) { cout << v.x << " " << v.y << " " << v.z << endl; }
	void CoordinateTransformation(xpoint viewVector)  //坐标变换
	{
		double a = viewVector.x, b = viewVector.y, c = viewVector.z;
		xpoint u = viewVector; xpoint v;
		if (c == 0) v = { 0,0,1 };
		if (c > 0) v = { -1 * a,-1 * b,(a * a + b * b) / c };
		else v = { a,b,-(a * a + b * b) / c };
		xpoint w = chaMultiplication(v, u);
		unitize(u); unitize(v); unitize(w);
		for (int i = 0; i<8; i++) transferPoint[i] = BaseCoordinateTransformation(u, v, w, originalPoint[i]);
		transferCenterPoint = BaseCoordinateTransformation(u, v, w, originalCenterPoint);
	}
	void setTransferFace()
	{
		for (int i = 0; i < 6; i++)
		{
			xpoint sum = { 0,0 };
			for (int j = 0; j < 4; j++)
			{
				faces[facesIndex[i]].transferPoint[j] = transferPoint[surfaceIndex[i][j]];
				sum = sum + transferPoint[surfaceIndex[i][j]];
			}
			faces[facesIndex[i]].transferCenterPoint = sum / 4;
		}
	}
public:
	void setOriginalCenterPoint(xpoint originalCenterPoint) { this->originalCenterPoint = originalCenterPoint; }
	void setOriginalPoint(xpoint originalPoint, int index) { this->transferPoint[index] = this->originalPoint[index] = originalPoint; }
	xpoint getOriginalCenterPoint() { return originalCenterPoint; }
	xpoint getTransferCenterPoint() { return transferCenterPoint; }
	xpoint getOriginalPoint(int index) { return originalPoint[index]; }
	void initialize(int index,faceOfCube *faces)
	{
		this->faces = faces;
		for (int i = 0; i < 6; i++)
		{
			//设置index,颜色,坐标,中心点
			facesIndex[i] = index * 6 + i;
			faces[facesIndex[i]].colorIndex = i;
			faces[facesIndex[i]].cubeIndex = index;
			faces[facesIndex[i]].cmpIndex = 0;
			xpoint sum = { 0,0 };
			for (int j = 0; j < 4; j++)
			{
				faces[facesIndex[i]].originalPoint[j] = originalPoint[surfaceIndex[i][j]];
				sum = sum + originalPoint[surfaceIndex[i][j]];
			}
			faces[facesIndex[i]].originalCenterPoint = sum / 4;
		}
		setTransferFace();
	}
	void shiftView(double alpha, double beta)  //传入旋转偏移量
	{
		double r = 1000;
		double x = r * sin(alpha)*cos(beta), y = r * sin(alpha)*sin(beta), z = r * cos(alpha);
		CoordinateTransformation({ x,y,z });
		setTransferFace();
	}
};
class magicCube
{
private:
	Rotation rotation;
	POINT basePoint;
	double baseR;
	double baseSize;
	stack<pair<char, int> > steps;//side,stateFlag
	cube *cubes;
	faceOfCube *faces;
	bool cmpIndex(int a, int b)
	{
		return cubes[faces[a].cmpIndex].getTransferCenterPoint()<cubes[faces[b].cmpIndex].getTransferCenterPoint() 
			|| faces[a].cmpIndex==faces[b].cmpIndex
			&& faces[a].transferCenterPoint.x < faces[b].transferCenterPoint.x;
	}
	xpoint calculateShiftUnitVecter(char coordinate,double degree)
	{
		if (coordinate == 'x')return{ 0,-1 * sin(degree),cos(degree) - 1 };
		if (coordinate == 'y')return{ sin(degree),0,cos(degree) - 1 };
		if (coordinate == 'z')return{ -1 * sin(degree),cos(degree) - 1,0 };
	}
	double calculateRelativeDegree(char coordinate,xpoint calcuPoint)
	{
		int x = calcuPoint.x, y = calcuPoint.y, z = calcuPoint.z;
		if (coordinate == 'x')return { atan2(z,y) - pi / 2 };
		if (coordinate == 'y')return { atan2(x,z) };
		if (coordinate == 'z')return { atan2(y,x) - pi / 2 };
	}
	double calculateMo(char coordinate, xpoint calcuPoint)
	{
		int x = calcuPoint.x, y = calcuPoint.y, z = calcuPoint.z;
		if (coordinate == 'x')return sqrt(y * y + z * z);
		if (coordinate == 'y')return sqrt(x * x + z * z);
		if (coordinate == 'z')return sqrt(x * x + y * y);
	}
	void turnDegree(char side, double degree)  //弧度,规定x-y(y=0),y-z(z=0),z-x(z=0)为弧度正方向,degree为相对转动,比如每次1°这样子转哈
	{
		xpoint symbolPoint;
		char symbolCoordinate;
		const double blankFlag = 0xABCDEF;
		double b = blankFlag;  //减少代码量且好看
		//小写表示中间旋转,如u=MU
		if (side == 'R') symbolPoint = { b,1,b };
		if (side == 'l') symbolPoint = { b,0,b };
		if (side == 'L') symbolPoint = { b,-1,b };
		if (side == 'U') symbolPoint = { b,b,1 };
		if (side == 'u') symbolPoint = { b,b,0 };
		if (side == 'D') symbolPoint = { b,b,-1 };
		if (side == 'F') symbolPoint = { 1,b,b };
		if (side == 'f') symbolPoint = { 0,b,b };
		if (side == 'B') symbolPoint = { -1,b,b };

		int index[9]; int p = 0;
		if (symbolPoint.x != b)
		{
			for (int i = 0; i < 27; i++)
				if (abs(cubes[i].getOriginalCenterPoint().x - symbolPoint.x * baseR) < esp) index[p++] = i;
			int tmpIndex[3];
			for (int i = 0; i < 27; i++)
			{
				if (cubes[i].getOriginalCenterPoint().y == 0 && cubes[i].getOriginalCenterPoint().z == 0)
				{
					if (cubes[i].getOriginalCenterPoint().x < 0)tmpIndex[0] = i;
					if (cubes[i].getOriginalCenterPoint().x == 0)tmpIndex[1] = i;
					if (cubes[i].getOriginalCenterPoint().x > 0)tmpIndex[2] = i;
				}
			}
			for (int i = 0; i < 162; i++)
			{
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().x < 0)faces[i].cmpIndex = tmpIndex[0];
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().x == 0)faces[i].cmpIndex = tmpIndex[1];
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().x > 0)faces[i].cmpIndex = tmpIndex[2];
			}
			symbolCoordinate = 'x';
		}
		if (symbolPoint.y != b)
		{
			for (int i = 0; i < 27; i++)
				if (abs(cubes[i].getOriginalCenterPoint().y - symbolPoint.y * baseR) < esp) index[p++] = i;
			int tmpIndex[3];
			for (int i = 0; i < 27; i++)
			{
				if (cubes[i].getOriginalCenterPoint().x == 0 && cubes[i].getOriginalCenterPoint().z == 0)
				{
					if (cubes[i].getOriginalCenterPoint().y < 0)tmpIndex[0] = i;
					if (cubes[i].getOriginalCenterPoint().y == 0)tmpIndex[1] = i;
					if (cubes[i].getOriginalCenterPoint().y > 0)tmpIndex[2] = i;
				}
			}
			for (int i = 0; i < 162; i++)
			{
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().y < 0)faces[i].cmpIndex = tmpIndex[0];
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().y == 0)faces[i].cmpIndex = tmpIndex[1];
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().y > 0)faces[i].cmpIndex = tmpIndex[2];
			}
			symbolCoordinate = 'y';
		}
		if (symbolPoint.z != b)
		{
			for (int i = 0; i < 27; i++)
				if (abs(cubes[i].getOriginalCenterPoint().z - symbolPoint.z * baseR) < esp) index[p++] = i;
			int tmpIndex[3];
			for (int i = 0; i < 27; i++)
			{
				if (cubes[i].getOriginalCenterPoint().x == 0 && cubes[i].getOriginalCenterPoint().y == 0)
				{
					if (cubes[i].getOriginalCenterPoint().z < 0)tmpIndex[0] = i;
					if (cubes[i].getOriginalCenterPoint().z == 0)tmpIndex[1] = i;
					if (cubes[i].getOriginalCenterPoint().z > 0)tmpIndex[2] = i;
				}
			}
			for (int i = 0; i < 162; i++)
			{
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().z < 0)faces[i].cmpIndex = tmpIndex[0];
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().z == 0)faces[i].cmpIndex = tmpIndex[1];
				if (cubes[faces[i].cubeIndex].getOriginalCenterPoint().z > 0)faces[i].cmpIndex = tmpIndex[2];
			}
			symbolCoordinate = 'z';
		}

		//计算每一个点的偏移
		for (int i = 0; i<9; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				xpoint nowXpoint = cubes[index[i]].getOriginalPoint(j);
				double oldDegree = calculateRelativeDegree(symbolCoordinate, nowXpoint);
				double nowDegree = oldDegree + degree;
				double mo = calculateMo(symbolCoordinate, nowXpoint);
				xpoint nowShiftVector = calculateShiftUnitVecter(symbolCoordinate, nowDegree) * mo;
				xpoint oldShiftVector = calculateShiftUnitVecter(symbolCoordinate, oldDegree) * mo;
				xpoint ShiftVector = nowShiftVector - oldShiftVector;
				cubes[index[i]].setOriginalPoint(nowXpoint + ShiftVector, j);
			}
			//中心点偏移
			xpoint nowXpoint = cubes[index[i]].getOriginalCenterPoint();
			double oldDegree = calculateRelativeDegree(symbolCoordinate, nowXpoint);
			double nowDegree = oldDegree + degree;
			double mo = calculateMo(symbolCoordinate, nowXpoint);
			xpoint nowShiftVector = calculateShiftUnitVecter(symbolCoordinate, nowDegree) * mo;
			xpoint oldShiftVector = calculateShiftUnitVecter(symbolCoordinate, oldDegree) * mo;
			xpoint ShiftVector = nowShiftVector - oldShiftVector;
			cubes[index[i]].setOriginalCenterPoint(nowXpoint + ShiftVector);
		}
	}
public:
	void printOriginalCenterPoint()
	{ 
		for (int i = 0; i < 27; i++)
			cout << i << ": " << cubes[i].getOriginalCenterPoint().x << " " << cubes[i].getOriginalCenterPoint().y << " " << cubes[i].getOriginalCenterPoint().z << endl;
		cout << endl;
	}
	void setbasePoint(POINT basePoint) { this->basePoint = basePoint; }
	void setRotation(Rotation rotation) { this->rotation = rotation; }
	void setBaseR(double baseR) { this->baseR = baseR; }
	void setBaseSize(double baseSize) { this->baseSize = baseSize; }
	Rotation getRotation() { return rotation; }
	POINT getBasePoint() { return basePoint; }
	void clear()
	{
		SelectObject(mydc, blackPen);
		SelectObject(mydc, blackBrush);
		Rectangle(mydc, 0, 0, baseX * 2, baseY * 2);
	}
	void initialize()
	{
		baseSize = 0.45;
		baseR =100.0;
		cubes = new cube[27];
		faces = new faceOfCube[162];
		for (int x = -1; x <=1; x++)
			for (int y = -1; y <= 1; y++)
				for (int z = -1; z <= 1; z++)
				{
					//cube的初始化
					int index = (x + 1) * 9 + (y + 1) * 3 + (z + 1);
					cubes[index].setOriginalCenterPoint({ x * baseR,y * baseR,z * baseR });
					for (int i = 0; i < 8; i++)
						cubes[index].setOriginalPoint({ (x + baseSize * basePoint1[i][0])*baseR,(y + baseSize * basePoint1[i][1])*baseR ,(z + baseSize * basePoint1[i][2])*baseR }, i);
					cubes[index].initialize(index,faces);
				}
		//设置魔方内部颜色
		for (int i = 0; i < 162; i++)
		{
			double x = faces[i].originalCenterPoint.x, y = faces[i].originalCenterPoint.y, z = faces[i].originalCenterPoint.z;
			if (x<(1 + baseSize)*baseR - esp && x>-((1 + baseSize)*baseR - esp)
				&& y<(1 + baseSize)*baseR - esp && y>-((1 + baseSize)*baseR - esp)
				&& z<(1 + baseSize)*baseR - esp && z>-((1 + baseSize)*baseR - esp))
				faces[i].colorIndex = 6;
		}
	}
	void shiftView(double alpha, double beta) { for (int i = 0; i < 27; i++) cubes[i].shiftView(alpha, beta); }
	void calibration()//校准原坐标
	{
		double calibrationOrdinaryConnstant[6] = { -(1 + baseSize)*baseR,-(1 - baseSize)*baseR,-baseSize * baseR,baseSize*baseR,(1 - baseSize)*baseR,(1 + baseSize)*baseR };
		double calibrationCenterConnstant[3] = { -baseR,0,baseR };
		for (int i = 0; i < 27; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				xpoint tmp = cubes[i].getOriginalPoint(j);
				for (int k = 0; k < 6; k++)
				{
					if (abs(tmp.x - calibrationOrdinaryConnstant[k]) < esp)tmp.x = calibrationOrdinaryConnstant[k];
					if (abs(tmp.y - calibrationOrdinaryConnstant[k]) < esp)tmp.y = calibrationOrdinaryConnstant[k];
					if (abs(tmp.z - calibrationOrdinaryConnstant[k]) < esp)tmp.z = calibrationOrdinaryConnstant[k];
					cubes[i].setOriginalPoint(tmp, j);
				}
			}
			xpoint tmp = cubes[i].getOriginalCenterPoint();
			for (int k = 0; k < 3; k++)
			{
				if (abs(tmp.x - calibrationCenterConnstant[k]) < esp)tmp.x = calibrationCenterConnstant[k];
				if (abs(tmp.y - calibrationCenterConnstant[k]) < esp)tmp.y = calibrationCenterConnstant[k];
				if (abs(tmp.z - calibrationCenterConnstant[k]) < esp)tmp.z = calibrationCenterConnstant[k];
				cubes[i].setOriginalCenterPoint(tmp);
			}
		}
	};
	void randomRotation(int step)
	{
		const int sleepTime = 20;
		const double turnSpeedIndex = 13.0;
		double rotationDegree = 0;
		if (step >= 20 && step < 50)rotationDegree = pi * 2;
		else if(step >= 50)rotationDegree = pi * 4;
		Rotation oldRotation = rotation;
		double perDegree1 = perPi * turnSpeedIndex;
		double perDegree2 = rotationDegree / (step * (90 / (int)(turnSpeedIndex)+1));
		srand(time(0));
		int oldNumber = -1;
		while (step--)
		{
			int randNumber = rand() % 18;
			int stateFlag = randNumber / 9;
			randNumber /= 2;
			if (randNumber == oldNumber)
			{
				step++;
				continue;
			}
			else
				oldNumber = randNumber;
			char side;
			if (randNumber == 0)side = 'R';
			if (randNumber == 1)side = 'l';
			if (randNumber == 2)side = 'L';
			if (randNumber == 3)side = 'U';
			if (randNumber == 4)side = 'u';
			if (randNumber == 5)side = 'D';
			if (randNumber == 6)side = 'F';
			if (randNumber == 7)side = 'f';
			if (randNumber == 8)side = 'B';

			const int normalFlag = stateFlag ? stateFlag : -1;
			steps.push({ side,normalFlag });
			double cnt = 0;
			while (cnt + turnSpeedIndex < 90)
			{
				turnDegree(side, perDegree1*normalFlag);
				oldRotation.beta += perDegree2;
				shiftView(rotation.alpha, oldRotation.beta);
				show();
				Sleep(sleepTime);
				cnt += turnSpeedIndex;
			}
			turnDegree(side, (90 - cnt)*perPi*normalFlag);
			oldRotation.beta += perDegree2;
			calibration();
			shiftView(rotation.alpha, oldRotation.beta);
			show();
		}

	}
	void turn(char side,int stateFlag)//-1逆1正-2两逆转2两顺转
	{
		steps.push({ side,stateFlag });
		const int sleepTime = 20;
		const double speedIndex = 7.0;
		double perDegree = perPi * speedIndex;
		double cnt = 0;
		while (cnt + speedIndex < 90)
		{
			turnDegree(side, perDegree*stateFlag);
			shiftView(rotation.alpha, rotation.beta);
			show();
			Sleep(sleepTime);
			cnt += speedIndex;
		}
		turnDegree(side, (90-cnt)*perPi*stateFlag);
		calibration();
		shiftView(rotation.alpha, rotation.beta);
		show();
	};
	void back(int step)
	{
		while (!steps.empty() && step--)
		{
			turn(steps.top().first, -steps.top().second);
			steps.pop(); steps.pop();//turn会重新加入。。。。。。
			Sleep(15);
		}
	}
	void verticalAndHorizontalRotation(double alpha, double beta, double useTime) //beta,相对角度
	{
		const int sleepTime = 10;
		Rotation oldRotation = rotation;
		rotation.alpha = alpha;
		rotation.beta += beta;
		double perDegree1 = (rotation.alpha - oldRotation.alpha) / useTime * sleepTime;
		double perDegree2 = beta / useTime * sleepTime;
		while (oldRotation.beta + perDegree2 < rotation.beta)
		{
			oldRotation.alpha += perDegree1;
			oldRotation.beta += perDegree2;
			shiftView(oldRotation.alpha, oldRotation.beta);
			show();
			Sleep(sleepTime);
		}
		shiftView(rotation.alpha, rotation.beta);
	}
	void verticalRotation(double alpha, double useTime)
	{
		const int sleepTime = 10;
		Rotation oldRotation = rotation;
		rotation.alpha = alpha;
		double perDegree = (rotation.alpha - oldRotation.alpha) / useTime * sleepTime;
		while (abs(oldRotation.alpha + perDegree - rotation.alpha)>abs(perDegree))
		{
			oldRotation.alpha += perDegree;
			shiftView(oldRotation.alpha, oldRotation.beta);
			show();
			Sleep(sleepTime);
		}
		shiftView(rotation.alpha, rotation.beta);
	}
	void showOneByOne(int sleepTime)
	{
		clear();
		shiftView(rotation.alpha, rotation.beta);
		int index[27] = { 6,7,8,15,16,17,24,25,26,3,4,5,12,13,14,21,22,23,0,1,2,9,10,11,18,19,20};
		int tmpFaceIndex[6];
		for (int k = 0; k < 27; k++)
		{
			int p = 0;
			for (int j = 0; j < 162; j++)if (faces[j].cubeIndex == index[k])tmpFaceIndex[p++] = j;
			for (int i = 5; i >= 0; i--)
				for (int j = 0; j < i; j++)
				{
					if (!(cmpIndex(tmpFaceIndex[j], tmpFaceIndex[j + 1])))
					{
						int tmp = tmpFaceIndex[j + 1];
						tmpFaceIndex[j + 1] = tmpFaceIndex[j];
						tmpFaceIndex[j] = tmp;
					}
				}
			if (k != 26)
			{
				for (int i = 0; i < 6; i++) faces[tmpFaceIndex[i]].show();
				draw();
				Sleep(sleepTime);
			}
		}
		BitBlt(tmpdc, 0, 0, baseX * 2, baseY * 2, mydc, 0, 0, SRCCOPY);//暂存在tmpdc
		//模拟抛物线方程
		xpoint facePoint[6][4];
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 4; j++)
				facePoint[i][j] = faces[tmpFaceIndex[i]].transferPoint[j];
		const double a =4000.0;
		const double minh = 4.0;
		const double maxh = 180.0;
		double k = 0.7;
		double deltak = 0.08;
		double h = maxh;
		double maxt = sqrt(2 * h / a);
		double t = 0;
		while (h > minh)
		{
			while (t < maxt)
			{
				BitBlt(mydc, 0, 0, baseX * 2, baseY * 2, tmpdc, 0, 0, SRCCOPY);
				double y = a * t * t / 2;
				for (int i = 0; i < 6; i++)
					for (int j = 0; j < 4; j++)
						faces[tmpFaceIndex[i]].transferPoint[j] = { facePoint[i][j].x,facePoint[i][j].y + h - y ,facePoint[i][j].z };
				for (int i = 0; i < 6; i++) faces[tmpFaceIndex[i]].show();
				draw();
				Sleep(sleepTime/4);
				t += 0.010;
			}
			h *= k;
			k -= deltak;
			maxt = sqrt(2 * h / a);
			t = -maxt;
		}
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 4; j++)
				faces[tmpFaceIndex[i]].transferPoint[j] = { facePoint[i][j].x,facePoint[i][j].y,facePoint[i][j].z };
		for (int i = 0; i < 6; i++) faces[tmpFaceIndex[i]].show();
		draw();
	}
	void show()
	{
		//每次绘图前清空dc
		clear();
		int index[162];
		for (int i = 0; i < 162; i++)index[i] = i;
		for(int i=162 - 1;i >= 0;i--)
			for (int j = 0; j < i; j++)
			{
				if (!(cmpIndex(index[j], index[j + 1])))
				{
					int tmp = index[j+1];
					index[j+1] = index[j];
					index[j] = tmp;
				}
			}
		for (int i = 0; i < 162; i++) faces[index[i]].show();
		draw();
	}
	void draw()		//从虚拟dc画到显示屏
	{
		BitBlt(hdc, basePoint.x - baseX, basePoint.y - baseY, baseX * 2, baseY * 2, mydc, 0, 0, SRCCOPY);
	}
};
class interaction
{
private:
	string strIn;
	int len;
	int unHideFlag;
	int preIndex = 7;
	magicCube myMagicCube;
public:
	interaction(magicCube myMagicCube) { this->myMagicCube = myMagicCube; }
	void operator()(int index)
	{
		if (index == 1)//start
		{
			//400,700->1440,500,40*26,40*(-5)
			//720°,40*(18*pi/180)
			//baseR:100->60,40*1
			int cnt = 0;
			const double useTime = 4.0;//其实是2.0，因为原定sleep(20)太卡改成了10
			const int maxCnt = 40;
			const int originalX = 400;
			const int originalY = 700;
			const double originalBaseR = 100.0;
			const double tmpBaseR = 60.0;
			for (double i = 1; i <= maxCnt * useTime; i++)
			{
				myMagicCube.setbasePoint({ originalX + (LONG)(i * 26 / useTime), originalY - (LONG)(i * 5 / useTime) });
				myMagicCube.shiftView(myMagicCube.getRotation().alpha, myMagicCube.getRotation().beta + i * 18 / useTime * pi / 180);
				myMagicCube.setBaseR(originalBaseR - i / useTime);
				myMagicCube.calibration();//使baseR的改变生效
				myMagicCube.show();
				Sleep(10);
			}
			for (double i = 1; i <= maxCnt * 4.0; i++)
			{
				double tmpAlpha = myMagicCube.getRotation().alpha + i * 2.25 * pi / 180.0;
				while (tmpAlpha > pi)tmpAlpha -= pi * 2;
				if (tmpAlpha > 0)
				{
					myMagicCube.shiftView(tmpAlpha, myMagicCube.getRotation().beta);
					hBrush[4] = CreateSolidBrush(basecolor[4]);
					hBrush[5] = CreateSolidBrush(basecolor[5]);
				}
				else
				{
					tmpAlpha += pi;
					myMagicCube.shiftView(tmpAlpha, myMagicCube.getRotation().beta + pi);
					hBrush[4] = CreateSolidBrush(basecolor[5]);
					hBrush[5] = CreateSolidBrush(basecolor[4]);
				}
				myMagicCube.setBaseR(tmpBaseR + i / 4.0);
				myMagicCube.calibration();
				myMagicCube.show();
				Sleep(10);
			}
			isTraval = 1;
			thread t(*this, 2);
			t.detach();
		}
		if (index == 2)//traval
		{
			finishFlag = 0;
			int alpha = myMagicCube.getRotation().alpha;
			int beta = myMagicCube.getRotation().alpha;
			const double speedConstant1 = 0.5;
			const double speedConstant2 = 0.5;
			int cnt = 0;
			while(isTraval)
			{
				double tmpAlpha = alpha + cnt * speedConstant1 * pi / 180;
				tmpAlpha -= (int)(tmpAlpha / (2 * pi)) * 2 * pi;//转换定义域到0-2pi
				tmpAlpha = (pi - abs(tmpAlpha - pi));//最终天顶角转换函数
				myMagicCube.shiftView(tmpAlpha, myMagicCube.getRotation().beta + cnt * speedConstant2 * pi / 180);
				myMagicCube.calibration();
				myMagicCube.show();
				cnt++;
				Sleep(10);
			}
			const int sleepTime = 10;
			for (int i = 0; i < baseY * 2; i+=10)
			{
				SelectObject(mydc, blackPen);
				SelectObject(mydc, blackBrush);
				Rectangle(mydc, 0, i, baseX * 2, i + 10);
				Sleep(sleepTime);
				myMagicCube.draw();
			}
			finishFlag = 1;
		}
		if (index >= 3)//固定像素
		{
			myMagicCube.shiftView(myMagicCube.getRotation().alpha, myMagicCube.getRotation().beta);
			myMagicCube.show();
			while (index==fixConstant) { myMagicCube.draw(); Sleep(20); }
		}
	}
	void clear() { strIn.clear(); len = 0; }
	void stopStaticDraw()
	{ 
		myMagicCube.clear();
		myMagicCube.draw();
		fixConstant++;
	}
	void staticDraw()
	{
		thread t(*this, fixConstant);
		t.detach();
	}
	void print(string strP[], int cnt) { for (int i = 0; i < cnt; i++)cout << "\t" << strP[i] << endl; cout << endl; }
	void strInPrint()
	{
		system("cls");
		cout << strIn << endl;
	}
	void errorPrint()
	{
		strInPrint();
		cout << "黑框框魔方：" << "无效输入" << endl;
	}
	void unShowPrint()
	{
		strInPrint();
		cout << "黑框框魔方：" << "魔方未显示,请使用appear命令" << endl;
	}
	void afterTravel()
	{
		while (!finishFlag)Sleep(1);
		unHideFlag = 0;
	}
	void readIn()
	{
		clear();
		while (!(len = strIn.size()))//非空输入
			getline(cin, strIn);
		if (isTraval)
		{
			isTraval = 0;
			afterTravel();
		}
		printJudge();
	}
	void printJudge()
	{
		const int cnt = 12;
		string str[cnt];
		str[0] = "appear";
		str[1] = "back";
		str[2] = "help";
		str[3] = "hide";
		str[4] = "play";
		str[5] = "setc";
		str[6] = "setv";
		str[7] = "show";
		str[8] = "sn";
		str[9] = "traval";
		str[10] = "upset";
		str[11] = "welcome";


		int len = strIn.size();
		if (strIn[0] == '/')
		{
			int commandIndex = -1;
			int shortIndex = -1;//-1表示无对应，0表示刚好，1表示有后续
			for (int i = 0; i < cnt; i++)
			{
				int pdlen = str[i].size();
				int flg = 1;
				for (int j = 0; j < pdlen; j++)
				{
					if (strIn[j + 1] != str[i][j])
					{
						flg = 0;
						break;
					}
				}
				if (flg)
				{
					commandIndex = i;
					if (strIn[pdlen + 1] == 0) shortIndex = 0;
					if (strIn[pdlen + 1] == ' ') shortIndex = 1;
					break;
				}
			}
			if (shortIndex == -1)errorPrint();
			else if (shortIndex == 0)//命令介绍
			{
				if (commandIndex == 0)appear();
				if (commandIndex == 1)back();
				if (commandIndex == 2)help();
				if (commandIndex == 3)hide();
				if (commandIndex == 4)play();
				if (commandIndex == 5)setc();
				if (commandIndex == 6)setv();
				if (commandIndex == 7)show();
				if (commandIndex == 8)sn();
				if (commandIndex == 9)traval();
				if (commandIndex == 10)upset();
				if (commandIndex == 11)welcome();
				preIndex = commandIndex;
			}
			else if (shortIndex == 1)//命令执行
			{
				stopStaticDraw();
				strInPrint();
				int p = 0;
				if (commandIndex == 0)
				{
					findNextPara(p);
					if (strIn[p] == '0' && strIn[p+1] == 0)
					{ 
						staticDraw();
						unHideFlag = 1;
						cout << "黑框框魔方：" << "显示成功" << endl;
						return;
					}
					else if (strIn[p] == '1' && strIn[p + 1] == ' ' || strIn[p + 1] == '-')
					{
						stopStaticDraw();
						int x = 0, y = 0;
						findNextPara(p);
						findNumber(x);
						findNextPara(p);
						findNumber(y);
						if (!x || !y) errorPrint();
						else
						{
							myMagicCube.setbasePoint({ x,y });
							staticDraw();
							unHideFlag = 1;
							cout << "黑框框魔方：" << "显示成功" << endl;
						}
						return;
					}
					errorPrint();
					staticDraw();
					return;
				}
				if (unHideFlag)
				{
					if (commandIndex == 1)
					{
						int n = 0;
						findNextPara(p);
						findNumber(n);
						if (n == 0) errorPrint();
						else
						{
							myMagicCube.back(n);
							cout << "还原成功" << endl;

						}
						staticDraw();
						return;
					}
					if (commandIndex == 5)
					{
						int x = 0, symbol = 0;
						findNextPara(p);
						if (strIn[p] >= '0' && strIn[p] <= '7' && strIn[p + 1] == ' ' || strIn[p + 1] == '-')
						{
							symbol = strIn[p] - '0';
							findNextPara(p);
							for (int i = 0; i < 6; i++)
							{
								int num = -1;
								if (strIn[p] <= '9' && strIn[p] >= '0')
									num = strIn[p] - '0';
								if (strIn[p] <= 'f' && strIn[p] >= 'a')
									num = 10 + strIn[p] - 'a';
								if (strIn[p] <= 'F' && strIn[p] >= 'A')
									num = 10 + strIn[p] - 'A';
								if (num == -1)
								{
									cout << "黑框框魔方：" << "颜色格式输入错误" << endl;
									staticDraw();
									return;
								}
								else
								{
									x *= 16;
									x += num;
								}
								p++;
							}
							basecolor[symbol] = x;
							initializePenAndBrush();
							myMagicCube.show();
							cout << "黑框框魔方：" << "设置成功" << endl;
							staticDraw();
							return;
						}
						errorPrint();
						staticDraw();
						return;
					}
					if (commandIndex == 6)
					{
						int aDegree = 0, bDegree = 0;
						findNextPara(p);
						findNumber(aDegree);
						findNextPara(p);
						findNumber(bDegree);
						if (!aDegree) errorPrint();
						else
						{
							if (aDegree >= 180) errorPrint();
							if (bDegree>360) errorPrint();
							if (aDegree < 180 && bDegree <= 360)
							{
								double alpha = aDegree * pi / 180;
								double beta = bDegree * pi / 180;

								if(beta==0) myMagicCube.verticalRotation(alpha, 800);
								else myMagicCube.verticalAndHorizontalRotation(alpha, beta, 800);
								double nowBeta = myMagicCube.getRotation().beta - (int)(myMagicCube.getRotation().beta / (pi * 2)) * pi * 2;
								cout << "黑框框魔方：" << "设置成功" << endl;
							}
						}
						staticDraw();
						return;
					}
					if (commandIndex == 10)
					{
						findNextPara(p);
						int n = 0;
						if (strIn[p] != 0 && strIn[p] >= '0'  && strIn[p] <= '9')
						{
							int tmpnum = strIn[p] - '0';
							if (strIn[p+1] != 0 )
							{
								if (strIn[p + 1] >= '0'  && strIn[p + 1] <= '9')
									n = tmpnum * 10 + strIn[p + 1] - '0';
								else errorPrint();
							}
							else n = tmpnum;
						}
						else errorPrint();
						myMagicCube.randomRotation(n);
						staticDraw();
						return;
					}
				}
				else if(commandIndex==1|| commandIndex == 5|| commandIndex == 6|| commandIndex == 10) unShowPrint();
			}
		}
		else//旋转魔方
		{
			if (unHideFlag)
			{
				stopStaticDraw();
				int p = 0;
				int flagP = -1;//防止进入死循环
				while (p < len)
				{
					if (flagP == p)
					{
						errorPrint();
						break;
					}
					else flagP = p;
					if (strIn[p] == ' ')
					{
						p++;
						continue;
					}
					if (strIn[p] == 'R')
					{
						int state = -1;
						if (strIn[p + 1] == '\'')state = 1;
						if (strIn[p + 1] == '2')state = -2;
						if (strIn[p + 1] == ' ')p--;
						myMagicCube.turn('R', state);
						p += 2;
						continue;
					}
					if (strIn[p] == 'L')
					{
						int state = 1;
						if (strIn[p + 1] == '\'')state = -1;
						if (strIn[p + 1] == '2')state = 2;
						if (strIn[p + 1] == ' ')p--;
						myMagicCube.turn('L', state);
						p += 2;
						continue;
					}
					if (strIn[p] == 'U')
					{
						int state = -1;
						if (strIn[p + 1] == '\'')state = 1;
						if (strIn[p + 1] == '2')state = -2;
						if (strIn[p + 1] == ' ')p--;
						myMagicCube.turn('U', state);
						p += 2;
						continue;
					}
					if (strIn[p] == 'D')
					{
						int state = 1;
						if (strIn[p + 1] == '\'')state = -1;
						if (strIn[p + 1] == '2')state = 2;
						if (strIn[p + 1] == ' ')p--;
						myMagicCube.turn('D', state);
						p += 2;
						continue;
					}
					if (strIn[p] == 'F')
					{
						int state = -1;
						if (strIn[p + 1] == '\'')state = 1;
						if (strIn[p + 1] == '2')state = -2;
						if (strIn[p + 1] == ' ')p--;
						myMagicCube.turn('F', state);
						p += 2;
						continue;
					}
					if (strIn[p] == 'B')
					{
						int state = 1;
						if (strIn[p + 1] == '\'')state = -1;
						if (strIn[p + 1] == '2')state = 2;
						if (strIn[p + 1] == ' ')p--;
						myMagicCube.turn('B', state);
						p += 2;
						continue;
					}
					if (strIn[p] == 'M')
					{
						if (strIn[p + 1] == 'L')
						{
							int state = 1;
							if (strIn[p + 2] == '\'')state = -1;
							if (strIn[p + 2] == '2')state = 2;
							if (strIn[p + 2] == ' ')p--;
							myMagicCube.turn('l', state);
							p += 3;
							continue;
						}
						if (strIn[p + 1] == 'U')
						{
							int state = -1;
							if (strIn[p + 2] == '\'')state = 1;
							if (strIn[p + 2] == '2')state = -2;
							if (strIn[p + 2] == ' ')p--;
							myMagicCube.turn('u', state);
							p += 3;
							continue;
						}
						if (strIn[p + 1] == 'F')
						{
							int state = -1;
							if (strIn[p + 2] == '\'')state = 1;
							if (strIn[p + 2] == '2')state = -2;
							if (strIn[p + 2] == ' ')p--;
							myMagicCube.turn('f', state);
							p += 3;
							continue;
						}
					}
				}
				staticDraw();
			}
			else unShowPrint();
		}
	}
	void appear()
	{
		system("cls");
		cout << "/appear" << endl;
		const int cnt = 7;
		string str[cnt];
		str[0] = "-flag，是否使用新坐标。0：使用之前隐藏时的坐标，1：使用新坐标";
		str[1] = "-x，屏幕的x坐标。原x：1440";
		str[2] = "-y，屏幕的y坐标。原y：500";
		str[3] = '\n';
		str[4] = "输入格式：/appear -flag 【-x】 【-y】";
		str[5] = "样例输入：/appear -0";
		str[6] = "样例说明：，如需快速演示输入样例，请输入/show";
		print(str, cnt);
	}
	void back()
	{

		system("cls");
		cout << "/back" << endl;
		const int cnt = 5;
		string str[cnt];
		str[0] = "-n，（n > 0）倒退步数，倒退至无法倒退则会自动停下";
		str[1] = '\n';
		str[2] = "输入格式：/back -n";
		str[3] = "样例输入：/back -2147483647";
		str[4] = "样例说明：如需快速演示输入样例，请输入/show";
		print(str, cnt);
	}
	void help()
	{
		system("cls");
		cout << "/help" << endl;
		const int cnt = 13;
		string str[cnt];
		str[0] = "/appear：在指定位置显示魔方";
		str[1] = "/back：倒退指定步数";
		str[2] = "/help：查看命令";
		str[3] = "/hide：隐藏魔方";
		str[4] = "/play：游戏玩法介绍 && 开始游戏";
		str[5] = "/setc：自定义魔方颜色";
		str[6] = "/setv：自定义魔方视角";
		str[7] = "/show：快速演示每个样例输入";
		str[8] = "/sn：辛马斯特标记介绍";
		str[9] = "/traval：缓动视角观察魔方";
		str[10] = "/upset：随机打乱魔方";
		str[11] = "/welcome：输出欢迎文字";
		str[12] = "注：分辨率1920 * 1080 缩放125%方可正常显示,若魔方显示在屏幕外请使用/appear";
		print(str, cnt);
	}
	void hide()
	{
		system("cls");
		cout << "/hide" << endl;
		if (unHideFlag)
		{
			stopStaticDraw();
			unHideFlag = 0;
			cout << "黑框框魔方：" << "隐藏成功" << endl;
		}
		else cout << "黑框框魔方：" << "我还没appear呢！！！" << endl;
	}
	void play()
	{
		isTraval = 0;
		Sleep(20);
		system("cls");
		cout << "/play" << endl;
		const int cnt = 6;
		string str[cnt];
		str[0] = "本游戏采用辛马斯特标记（Singmaster notation）表示旋转。";
		str[1] = '\n'+"视角为左视，即左面与前面显示在前"+'\n';
		str[2] = "输入格式：直接输入标记即可旋转魔方。";
		str[3] = "样例输入：MU MF ML MU' MF2 R L U D R' L' U2 D2";
		str[4] = "样例说明：M表示中层，必须大写（按下CUPS LOCK即可），角标 '表示逆时针旋转，2表示两圈";
		str[5] = "如需更详细了解辛马斯特标记，请输入/sn,如需快速演示输入样例，请输入/show";
		print(str, cnt);
		staticDraw();
		unHideFlag = 1;
	}
	void setc()
	{
		system("cls");
		cout << "/setc" << endl;
		const int cnt = 6;
		string str[cnt];
		str[0] = "-symbol：前后左右上下的symbol分别为012345，内表面代号6，边缘代号7。";
		str[1] = "-color：BBGGRR, BB：8bit表示蓝色；GG：8bit表示绿色；RR：8bit表示红色。";
		str[2] = '\n';
		str[3] = "输入格式：/setc -symbol -color";
		str[4] = "样例输入：/setc -0 -0055FF";
		str[5] = "样例说明：不区分大小写,如需快速演示输入样例，请输入/show";
		print(str, cnt);
	}
	void setv()
	{
		system("cls");
		cout << "/setv" << endl;
		const int cnt = 6;
		string str[cnt];
		str[0] = "-alpha：距z轴正方向的绝对角度。（0<alpha<180）。";
		str[1] = "-beta：俯视时顺时针旋转的相对角度。（0<=beta<=360）";
		str[2] = '\n';
		str[3] = "输入格式：/setv -alpha -beta";
		str[4] = "样例输入：/setv -140 -360";
		str[5] = "样例说明：如需快速演示输入样例，请输入/show";
		print(str, cnt);
	}
	void show()
	{
		clear();
		if (preIndex == 0) strIn = "/appear -0";
		if (preIndex == 1) strIn = "/back -2147483647";
		if (preIndex == 4) strIn = "MU MF ML MU' MF2 R L U D R' L' U2 D2";
		if (preIndex == 5) strIn = "/setc -0 -0055FF";
		if (preIndex == 6) strIn = "/setv -140 -360";
		if (preIndex == 10) strIn = "/upset -48";
		len = strIn.size();
		printJudge();
	}
	void sn()
	{
		system("cls");
		cout << "/sn" << endl;
		const int cnt = 13;
		string str[cnt];
		str[0] = "辛马斯特标记乃是通用标准，通常被俗称为“魔方公式符号”，于1978年12月由辛马斯特发明。";
		str[1] = "魔方各层以英文首字母指代。";
		str[2] = "R（Right）、L（Left）、U（Up）、D（Down）、F（Front）、B（Back）";
		str[3] = "分别指代右、左、顶（上）、底（下）、正（前）、背（后）层。";
		str[4] = "顺时针旋转90°：无后缀；逆时针旋转90°：后缀【'】；旋转180°：后缀【2】。";
		str[5] = '\n';
		str[6] = "完整的辛马斯特标记可以理解为【以面向指代层的视角，按方向进行旋转】。";
		str[7] = "例如：R，以面向右面视角，将右面顺时针旋转90°。从正面视角来看，即右面“向上”转90°。";
		str[8] = "又例如：D，以面向底面视角，将底面顺时针旋转90°。从正面视角来看，即右面“向右”转90°。";
		str[9] = "又例如：B'，以面向背面视角，将背面逆时针旋转90°。从正面视角来看，即背面“向右”转90°。";
		str[10] = '\n';
		str[11] = "除此之外，若要记录更加详细的魔方转动，还会用到：M（Middle）与U、F、L合用，指代各中层；";
		str[12] = "例如：MU，以顶面视角，将中间层顺时针旋转90°。从正面视角来看，即上数第二层“向左”转90°。";
		print(str, cnt);
	}
	void traval()
	{
		if (unHideFlag)
		{
			if (!isTraval)
			{
				stopStaticDraw();
				isTraval = 1;
				thread t(*this, 2);
				t.detach();
			}
		}
		else unShowPrint();
	}
	void upset()
	{
		system("cls");
		cout << "/upset" << endl;
		const int cnt = 5;
		string str[cnt];
		str[0] = "-n：随机旋转魔方n次，（n <= 99），若n>=20则可体验视角旋转效果。";
		str[1] = '\n';
		str[2] = "输入格式：/upset -n";
		str[3] = "样例输入：/upset -48";
		str[4] = "样例说明：如需快速演示输入样例，请输入/show";
		print(str, cnt);
	}
	void welcome()
	{
		system("cls");
		cout << "/welcome" << endl;
		const int cnt = 4;
		string str[cnt];
		str[0] = "欢迎来到控制台魔方世界!!!";
		str[1] = "如果您是一个魔方爱好者，同时又是一个程序员，那么此款游戏将绝对会符合您的胃口。";
		str[2] = "在这里，您可以通过命令输入的形式，自定义自己的魔方，并操控这个属于自己的魔方。";
		str[3] = "现在请您输入/help来查看有哪些命令吧!";
		print(str, cnt);
	}
	void welcomePrintOneByOne(int sleepTime)
	{
		const int cnt = 4;
		string str[cnt];
		str[0] = "Welcome to the world of Console-Magic-Cube!!!";
		str[1] = "如果您是一个魔方爱好者，同时又是一个程序员，那么此款游戏将绝对会符合您的胃口。";
		str[2] = "在这里，您可以通过命令输入的形式，自定义自己的魔方，并操控这个属于自己的魔方。";
		str[3] = "现在请您输入/help来查看有哪些命令吧!";
		for (int i = 0; i < cnt; i++)
		{
			int tmplen = str[i].size();
			cout << '\t' << '\t' << '\t' << '\t' << '\t' << '\t' << '\t' << '\t' << '\t';
			for (int j = 0; j < tmplen; j++) { cout << str[i][j]; if (j != tmplen - 1)Sleep(sleepTime); }
			cout << endl;
		}
		clear();
	}
	void start()
	{
		myMagicCube.showOneByOne(40);
		thread t(*this,1);
		t.detach();
		welcomePrintOneByOne(30);
		myMagicCube.setbasePoint({ 1440,500 });
		cout << '\t' << '\t' << '\t' << '\t' << '\t' << '\t' << '\t' << '\t' << '\t';
		system("pause");
		system("cls");
	}
};

int main()
{
	initializeGDI();
	initializeConsole();

	Sleep(1000);
	magicCube myMagicCube;
	myMagicCube.initialize();
	myMagicCube.setbasePoint({ 400,700 });
	myMagicCube.setRotation({ 1,-pi/4 });
	
	interaction myInteraction(myMagicCube);
	myInteraction.start();
	while (1)
	{
		myInteraction.readIn();
	}

	system("pause");
	return 0;
}