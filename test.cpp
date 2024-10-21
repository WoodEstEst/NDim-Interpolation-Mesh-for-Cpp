#include "test.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include "SimpleTensor.h"
#include "onedim_InterpMesh.h"
#include "ndim_InterpMesh.h"

std::vector<double> random_x(double min, double max, int div) {
	std::mt19937 gen;
	std::uniform_real_distribution<> dis(min, max);
	std::vector<double> x;
	x.push_back(min);
	x.push_back(max);
	for (int i = 0; i < div - 2; ++i) {
		x.push_back(dis(gen));
	}
	std::sort(x.begin(), x.end());
	////显示x
	//for (int i = 0; i < x.size(); ++i) {
	//	std::cout << x[i] << '\n';
	//}
	return x;
}

void test_onedim()
{
	//生成数据点
	std::vector<double> x = random_x(0, 2, 10);
	std::vector<std::array<double,1>> y;
	for (int i = 0; i < x.size(); ++i) { y.push_back({ std::sin(M_PI*x[i]) }); }
	for (int i = 0; i < x.size(); ++i) { std::cout << x[i] << '\n'; }
	std::cout << '\n';
	for (int i = 0; i < x.size(); ++i) { std::cout << y[i][0] << '\n'; }
	std::cout << '\n';
	//生成插值网格
	oneDim_InterpMesh<double,1> mesh(x, y, InterpolationType::Hermite);
	//采样0:0.01:2
	std::vector<double> x_sp;
	std::vector<std::array<double, 1>> y_sp;
	for (double xi = 0.0; xi <= 2.0; xi += 0.01) {
		x_sp.push_back(xi);
		y_sp.push_back(mesh.getPointVal(xi));
	}
	for (int i = 0; i < x_sp.size(); ++i) { std::cout << x_sp[i] << '\n'; }
	std::cout << '\n';
	for (int i = 0; i < x_sp.size(); ++i) { std::cout << y_sp[i][0] << '\n'; }
	std::cout << '\n';

}

void test_ndim()
{
	//生成数据点
	std::vector<double> x1 = random_x(0, 2, 10);
	std::vector<double> x2 = random_x(0, 1, 10);
	//产生展开的二维数据
	std::ofstream outp("points.csv");
	std::vector<std::array<double, 1>> y;
	for (int i1 = 0; i1 < x1.size(); ++i1) {
		for (int i2 = 0; i2 < x2.size(); ++i2) {
			y.push_back({ x2[i2] * x2[i2] * sin(M_PI * x1[i1]) });
			outp << x1[i1] << "," << x2[i2] << "," << y.back()[0] << '\n';
		}
	}
	outp.close();
	std::vector<std::vector<double>> x_ndim = { x1, x2 };
	//生成插值网格
	NDim_InterpMesh<double,2,1> mesh(x_ndim, y, InterpolationType::Akima);
	//采样0:0.01:2,0:0.1:1,结果输出到csv文件
	std::ofstream outYsp("meshgridY.csv");
	std::vector<double> x1_sp, x2_sp;
	std::vector<std::vector<std::array<double, 1>>> meshgrid_Y_sp;
	for (double xi1 = 0.0; xi1 <= 2.0; xi1 += 0.02) { x1_sp.push_back(xi1); }
	for (double xi2 = 0.0; xi2 <= 1.0; xi2 += 0.01) { x2_sp.push_back(xi2); }
	meshgrid_Y_sp.resize(x1_sp.size(), std::vector<std::array<double, 1>>(x2_sp.size()));
	for (int i1 = 0; i1 < x1_sp.size(); ++i1) {
		for (int i2 = 0; i2 < x2_sp.size(); ++i2) {
			meshgrid_Y_sp[i1][i2] = mesh.getPointVal({ x1_sp[i1],x2_sp[i2] });
			outYsp << meshgrid_Y_sp[i1][i2][0];
			if (i2 < x2_sp.size() - 1) { outYsp << ","; }
		}
		if (i1 < x1_sp.size() - 1) { outYsp << '\n'; }
	}
	outYsp.close();
}
