#pragma once
#include <functional>
#include <windows.h>
#include "Interpolator.h"
/*
说明：该类是对单维插值的封装，
*/

//一维插值网格
template<class DataType, size_t VALNUM>
class oneDim_InterpMesh
{
public:
	//默认构造的插值类型为Akima
	oneDim_InterpMesh(
		const std::vector<double>& x_vec,
		const std::vector<std::array<DataType, VALNUM>>& val_vec,
		InterpolationType interp_type = InterpolationType::Akima
	) :
		m_x(x_vec), m_val(val_vec)
	{
		//尺寸检查
		if (val_vec.size() != x_vec.size()) { exit(114514); }
		m_size = x_vec.size();
		//检查x序列的排序是否递增
		for (size_t i = 1; i < m_size; i++) {
			if (m_x[i] <= m_x[i - 1]) { exit(114514); }
		}
		//形成映射
		val_map = [this](const size_t iq) {
			return this->m_val[iq];
			};
		switch (interp_type) {
		case InterpolationType::Linear:
		{
			val_interp = LinearInterpolation(val_map, m_x);
			break;
		}
		case InterpolationType::Hermite:
		{
			val_interp = HermiteInterpolation(val_map, m_x);
			break;
		}
		case InterpolationType::Akima:
		{
			val_interp = AkimaInterpolation(val_map, m_x);
			break;
		}
		default:
			exit(114514);
			break;
		}
	};

protected:	
	typedef std::array<DataType, VALNUM>	ValArrayType;
	std::vector<ValArrayType>	m_val;		//数据数组	
	std::vector<double>			m_x;		//刻度数组	
	size_t						m_size;		//刻度个数
	size_t						x_loc_mark;	//动态插值坐标储存
	std::function<ValArrayType(size_t)> val_map;
	std::function<ValArrayType(double, size_t)> val_interp;
public:
	//获取插值点的值
	ValArrayType getPointVal(const double xq)
	{
		size_t x_loc = x_loc_mark;
		bool out_of_range = false;
		//动态插值位置判断
		if (!(m_x[x_loc] <= xq && m_x[x_loc + 1] > xq))
		{
			//遍历判断区间
			if (m_x[0] > xq || m_x[m_size - 1] < xq) {
				out_of_range = true;
			}
			else {
				for (size_t i = 1; i < m_size; i++) {
					if (m_x[i] > xq) {
						x_loc = i - 1;
						break;
					}
				}
			}
		}
		//出界返回报错
		if (out_of_range) {
			std::cout << "插值点越界" << std::endl;
			exit(114514);
		}
		x_loc_mark = x_loc;
		return val_interp(xq, x_loc);
	}
};

