#pragma once
/*
说明：
基于GitHub lecopivo的N-dimensional interpolation library项目思路改造的n维插值功能代码，
插值方式有线性和makima插值，补充了边界条件。

参考项目：https://github.com/lecopivo/Interpolation
参考：
Hermite三次插值与修正Makima：
	https://blog.csdn.net/dsn0606/article/details/106389133
	https://ww2.mathworks.cn/help/matlab/ref/makima.html
	文档论文"akima1970.pdf"
三次样条与边界条件：
	https://www.cnblogs.com/flysun027/p/10371726.html
	文档"spline_note.pdf"
*/

#include <cmath>
#include <type_traits>
#include <vector>
#include <array>
#include <functional>

//插值类型枚举
enum class InterpolationType
{
	Linear,
	Hermite,
	Akima
};

//一维插值器
/*
函数说明：
val_map		(i)			-> sample_val
x_map		(i)			-> sample_x
val_mesh	(xq, loc)	-> interp_val_array
*/

//线性插值
auto LinearInterpolation = [](auto val_map, std::vector<double>& x_map) {
	//插值器生成：生成一个输入为 xq和i，输出为double vq的函数
	//i_loc为查询点的区间编号，即查询点前一个节点编号，在外部判断
	return [=](double xq, size_t i_loc) mutable {
		//vq_arr是val_map返回值
		using ValArrayType = std::remove_reference_t<decltype(val_map(0))>;
		ValArrayType vq_arr;
		ValArrayType v_arr0 = val_map(i_loc);
		ValArrayType v_arr1 = val_map(i_loc + 1);
		double x0 = x_map.at(i_loc);
		double x1 = x_map.at(i_loc + 1);
		double wx = (xq - x0) / (x1 - x0);
		for (size_t val_i = 0; val_i < v_arr0.size(); val_i++) {
			vq_arr.at(val_i) = wx * v_arr1.at(val_i) + (1 - wx) * v_arr0.at(val_i);
		}
		return vq_arr;
		};
	};

//Hermite三次插值 一维至少需要3个点
auto HermiteInterpolation = [](auto val_map, std::vector<double>& x_map) {
	return [=](double xq, size_t i_loc) mutable {
		//vq_arr是val_map返回值		
		using ValArrayType = std::remove_reference_t<decltype(val_map(0))>;
		ValArrayType vq_arr;
		std::array<ValArrayType, 4> v_arr;
		std::array<ValArrayType, 3> delta_arr; //区间斜率
		ValArrayType slope_arr1, slope_arr2; //点斜率
		std::array<double, 4> x;
		size_t i_ini = 0;
		size_t i_end = x_map.size() - 2;
		size_t boundary_margin_left = i_loc - i_ini;
		size_t boundary_margin_right = i_end - i_loc;
		//区间边界情况判断和提取值
		v_arr[1] = val_map(i_loc);
		v_arr[2] = val_map(i_loc + 1);
		x[1] = x_map.at(i_loc);
		x[2] = x_map.at(i_loc + 1);
		double hx = x[2] - x[1];
		double dx = xq - x[1];
		double wx = dx / hx;
		if (boundary_margin_right >= 1) {
			v_arr[3] = val_map(i_loc + 2);
			x[3] = x_map.at(i_loc + 2);
		}
		if (boundary_margin_left >= 1) {
			v_arr[0] = val_map(i_loc - 1);
			x[0] = x_map.at(i_loc - 1);
		}
		for (size_t val_i = 0; val_i < vq_arr.size(); val_i++) {
			//需要所在区间及往外一个区间的斜率, 边界情况考虑特殊方法插值
			delta_arr[1][val_i] = (v_arr[2][val_i] - v_arr[1][val_i]) / hx;
			if (boundary_margin_right >= 1) {
				delta_arr[2][val_i] = (v_arr[3][val_i] - v_arr[2][val_i]) / (x[3] - x[2]);
			}
			if (boundary_margin_left >= 1) {
				delta_arr[0][val_i] = (v_arr[1][val_i] - v_arr[0][val_i]) / (x[1] - x[0]);
			}
			//left+right一定>=1
			if (boundary_margin_right < 1) {
				delta_arr[2][val_i] = 2 * delta_arr[1][val_i] - delta_arr[0][val_i];
			}
			if (boundary_margin_left < 1) {
				delta_arr[0][val_i] = 2 * delta_arr[1][val_i] - delta_arr[2][val_i];
			}
			//估计点斜率
			slope_arr1[val_i] = (v_arr[2][val_i] - v_arr[1][val_i]) / (x[2] - x[1]);
			slope_arr2[val_i] = (v_arr[3][val_i] - v_arr[2][val_i]) / (x[3] - x[2]);
			//Hermite三次插值
			vq_arr[val_i] = pow(1 - wx, 2) * ((1 + 2 * wx) * v_arr[1][val_i] + dx * slope_arr1[val_i])
				+ pow(wx, 2) * ((1 + 2 * (1 - wx)) * v_arr[2][val_i] + (dx - hx) * slope_arr2[val_i]);
		}
		return vq_arr;
		};
	};

//修正akima插值 一维至少需要3个点
auto AkimaInterpolation = [](auto val_map, std::vector<double>& x_map) {
	return [=](double xq, size_t i_loc) mutable {
		//vq_arr是val_map返回值		
		using ValArrayType = std::remove_reference_t<decltype(val_map(0))>;
		ValArrayType vq_arr;
		std::array<ValArrayType, 6> v_arr;
		std::array<ValArrayType, 5> delta_arr; //区间斜率
		std::array<ValArrayType, 4> w_arr;
		std::array<ValArrayType, 4> normw_arr;
		ValArrayType slope_arr2, slope_arr3;//点斜率
		std::array<double, 6> x;
		size_t i_ini = 0;
		size_t i_end = x_map.size() - 2;
		size_t boundary_margin_left = i_loc - i_ini;  //插入区间到边界的距离
		size_t boundary_margin_right = i_end - i_loc; //插入区间到边界的距离
		//区间边界情况判断和提取值
		v_arr[2] = val_map(i_loc);
		v_arr[3] = val_map(i_loc + 1);
		x[2] = x_map.at(i_loc);
		x[3] = x_map.at(i_loc + 1);
		double hx = x[3] - x[2];
		double dx = xq - x[2];
		double wx = dx / hx;
		if (boundary_margin_right >= 1) {
			v_arr[4] = val_map(i_loc + 2);
			x[4] = x_map.at(i_loc + 2);
		}
		if (boundary_margin_left >= 1) {
			v_arr[1] = val_map(i_loc - 1);
			x[1] = x_map.at(i_loc - 1);
		}
		if (boundary_margin_right >= 2) {
			v_arr[5] = val_map(i_loc + 3);
			x[5] = x_map.at(i_loc + 3);
		}
		if (boundary_margin_left >= 2) {
			v_arr[0] = val_map(i_loc - 2);
			x[0] = x_map.at(i_loc - 2);
		}
		for (size_t val_i = 0; val_i < vq_arr.size(); val_i++) {
			//需要所在区间及往外两个区间的斜率, 边界情况考虑特殊方法插值
			delta_arr[2][val_i] = (v_arr[3][val_i] - v_arr[2][val_i]) / hx;
			if (boundary_margin_right >= 1)
				delta_arr[3][val_i] = (v_arr[4][val_i] - v_arr[3][val_i]) / (x[4] - x[3]);
			if (boundary_margin_left >= 1)
				delta_arr[1][val_i] = (v_arr[2][val_i] - v_arr[1][val_i]) / (x[2] - x[1]);
			//left+right一定>=1
			if (boundary_margin_right >= 2)
				delta_arr[4][val_i] = (v_arr[5][val_i] - v_arr[4][val_i]) / (x[5] - x[4]);
			else if (boundary_margin_right >= 1)
				delta_arr[4][val_i] = 2 * delta_arr[3][val_i] - delta_arr[2][val_i];
			else {
				delta_arr[3][val_i] = 2 * delta_arr[2][val_i] - delta_arr[1][val_i];
				delta_arr[4][val_i] = 3 * delta_arr[2][val_i] - 2 * delta_arr[1][val_i];
			}
			if (boundary_margin_left >= 2)
				delta_arr[0][val_i] = (v_arr[1][val_i] - v_arr[0][val_i]) / (x[1] - x[0]);
			else if (boundary_margin_left >= 1)
				delta_arr[0][val_i] = 2 * delta_arr[1][val_i] - delta_arr[2][val_i];
			else {
				delta_arr[1][val_i] = 2 * delta_arr[2][val_i] - delta_arr[3][val_i];
				delta_arr[0][val_i] = 3 * delta_arr[2][val_i] - 2 * delta_arr[3][val_i];
			}
		
			//估计点斜率(修正法的权重)
			w_arr[0][val_i] = abs(delta_arr[3][val_i] - delta_arr[2][val_i]) + abs(delta_arr[3][val_i] + delta_arr[2][val_i]) / 2;
			w_arr[1][val_i] = abs(delta_arr[1][val_i] - delta_arr[0][val_i]) + abs(delta_arr[1][val_i] + delta_arr[0][val_i]) / 2;
			w_arr[2][val_i] = abs(delta_arr[4][val_i] - delta_arr[3][val_i]) + abs(delta_arr[4][val_i] + delta_arr[3][val_i]) / 2;
			w_arr[3][val_i] = abs(delta_arr[2][val_i] - delta_arr[1][val_i]) + abs(delta_arr[2][val_i] + delta_arr[1][val_i]) / 2;
			//归一化权重(同时防止除数为0)
			if (w_arr[0][val_i] + w_arr[1][val_i] <= 1E-60) { w_arr[0][val_i] = 1.0; w_arr[1][val_i] = 1.0; }
			if (w_arr[2][val_i] + w_arr[3][val_i] <= 1E-60) { w_arr[2][val_i] = 1.0; w_arr[3][val_i] = 1.0; }
			normw_arr[0][val_i] = w_arr[0][val_i] / (w_arr[0][val_i] + w_arr[1][val_i]);
			normw_arr[1][val_i] = w_arr[1][val_i] / (w_arr[0][val_i] + w_arr[1][val_i]);
			normw_arr[2][val_i] = w_arr[2][val_i] / (w_arr[2][val_i] + w_arr[3][val_i]);
			normw_arr[3][val_i] = w_arr[3][val_i] / (w_arr[2][val_i] + w_arr[3][val_i]);
			//估计斜率
			slope_arr2[val_i] = (normw_arr[0][val_i] * delta_arr[1][val_i] + normw_arr[1][val_i] * delta_arr[2][val_i]);
			slope_arr3[val_i] = (normw_arr[2][val_i] * delta_arr[2][val_i] + normw_arr[3][val_i] * delta_arr[3][val_i]);
			//Hermite三次插值
			vq_arr[val_i] = pow(1 - wx, 2) * ((1 + 2 * wx) * v_arr[2][val_i] + dx * slope_arr2[val_i])
				+ pow(wx, 2) * ((1 + 2 * (1 - wx)) * v_arr[3][val_i] + (dx - hx) * slope_arr3[val_i]);
		}
		return vq_arr;
		};
	};


//n维线性插值器生成
/*
函数说明：
val_mapInitial	({i_dim1,i_dim2...})					-> sample_val_array
val_mapFinal	({xq_dim1,xq_dim1}, {i_dim1,i_dim2...})	-> interp_val_array
x_map			(i,dim)									-> ith sample_x at dim
*/


//对val_map第I维插值: 
//val_mapI({xq_dim1,...,xq_dimI-1}, {i_dimI,...,i_dimN})	
//转化为
//val_mapI_next({xq_dim1,...,xq_dimI}, {i_dimI+1,...,i_dimN})
template <size_t I, size_t DIM>
auto Interpolate_IthArg = [](auto val_mapI, std::vector<std::vector<double>>& x_map, auto interpolation) {
	//生成第I维为xq替换i的总映射
	auto val_mapI_next = [=,&x_map](
		std::array<double, DIM> xqs,
		std::array<size_t, DIM> iqs, 
		std::array<size_t, DIM> i_locs) mutable {
		//生成val_mapI在I维方向一维离散点映射
		auto dimI_val_map = [val_mapI, xqs, iqs, i_locs](size_t i) mutable {
			iqs.at(I) = i;
			return val_mapI(xqs, iqs, i_locs);
		};
		double xqI = xqs.at(I);
		size_t i_loc = i_locs.at(I);
		//生成一个xqI到vq的映射
		auto interp_Idim = interpolation(dimI_val_map, x_map.at(I));
		auto vq = interp_Idim(xqI, i_loc);
		return vq;
	};
	return val_mapI_next;
};


//插值器生成递归, 将val_map逐维度插值
template <size_t I, size_t DIM, class ValMapI, class Interpolation>
auto InterpolationDimWise_iter(ValMapI val_mapI, std::vector<std::vector<double>>& x_map, Interpolation interpolation) {
	if constexpr (I == DIM) {
		return val_mapI;
	}
	else {
		//生成第I维为xq替换i的映射
		auto next_val_mapI = Interpolate_IthArg<I, DIM>(val_mapI, x_map, interpolation);
		//递归调用插值下一维
		return InterpolationDimWise_iter<I + 1, DIM>(next_val_mapI, x_map, interpolation);
	}
};

//插值器生成函数
template <size_t DIM, class ValMap, class Interpolation>
auto InterpolationDimWise(ValMap val_map, std::vector<std::vector<double>>& x_map, Interpolation interpolation) {
	//构造下标、插值点坐标混合参数映射迭代函数参数
	auto val_mapInitial = [val_map](
		std::array<double, DIM> xqs,
		std::array<size_t, DIM> iqs,
		std::array<size_t, DIM> i_locs) mutable {
			return val_map(iqs);
		};
	//调用迭代插值
	auto val_mapFinal = InterpolationDimWise_iter<0, DIM>(val_mapInitial, x_map, interpolation);
	return val_mapFinal;
};

//剔除iqs下标参数，返回插值点坐标插值函数
template<size_t DIM, class ValMapI>
auto InterpolationRemoveIqs(ValMapI val_mapFinal) {
	auto val_interpMesh = [val_mapFinal](
		std::array<double, DIM> xqs,
		std::array<size_t, DIM> i_locs) mutable {
			std::array<size_t, DIM> iqs{ 0 };
			return val_mapFinal(xqs, iqs, i_locs);
		};
	return val_interpMesh;
}
