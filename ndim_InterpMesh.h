#pragma once
#include <functional>
#include <windows.h>
#include "SimpleTensor.h"
#include "Interpolator.h"
/*
说明：该类是对多维插值的封装，
*/

//n维插值网格
template<class DataType, size_t DIM, size_t VALNUM>
class NDim_InterpMesh
{
public:
	//默认构造的插值类型为Akima
	NDim_InterpMesh(
		const std::vector<std::vector<double>>& x_ndim,
		const std::vector<std::array<DataType, VALNUM>>& val_flattensor,
		InterpolationType interp_type = InterpolationType::Akima
		) :
		m_x(x_ndim), m_size({ 0 }), x_locs_mark({ 0 })
	{
		//数据尺寸检查
		if (m_x.size() != DIM) {
			exit(114514);
		}
		for (size_t dim = 0; dim < DIM; dim++) {
			m_size.at(dim) = m_x.at(dim).size();
			if (interp_type == InterpolationType::Linear && m_size.at(dim) < 2) {
				exit(114514);
			}
			else if (m_size.at(dim) < 3) {
				exit(114514);
			}
		}
		//坐标升序检查
		for (auto& x : m_x) {
			for (size_t i = 1; i < x.size(); i++) {
				if (x.at(i) < x.at(i - 1)) exit(114514);
			}
		}
		//数据放入张量
		m_val_tensor = new SimpleTensor<DataType, DIM, VALNUM>(val_flattensor, m_size);
		//形成映射
		val_map = [this](const std::array<size_t, DIM>& iqs) {
			return this->m_val_tensor->at_idxs(iqs);
			};
		switch (interp_type) {
		case InterpolationType::Linear:
		{
			auto val_mapFinal = InterpolationDimWise<DIM>(val_map, m_x, LinearInterpolation);
			val_interpMesh = InterpolationRemoveIqs<DIM>(val_mapFinal);
			break;
		}

		case InterpolationType::Hermite:
		{
			auto val_mapFinal = InterpolationDimWise<DIM>(val_map, m_x, HermiteInterpolation);
			val_interpMesh = InterpolationRemoveIqs<DIM>(val_mapFinal);
			break;
		}

		case InterpolationType::Akima:
		{
			auto val_mapFinal = InterpolationDimWise<DIM>(val_map, m_x, AkimaInterpolation);
			val_interpMesh = InterpolationRemoveIqs<DIM>(val_mapFinal);
			break;
		}
		default:
			exit(114514);
			break;
		}
	};
	~NDim_InterpMesh() { delete m_val_tensor; }

protected:
	typedef std::array<DataType, VALNUM>	ValArrayType;
	SimpleTensor<DataType, DIM, VALNUM>*	m_val_tensor;	//数据张量指针
	std::vector<std::vector<double>>		m_x;			//网格每维坐标刻度
	std::array<size_t, DIM>					m_size;			//网格每维坐标数
	std::array<size_t, DIM>					x_locs_mark;	//动态插值坐标储存
	std::function<ValArrayType(const std::array<size_t, DIM>&)> val_map;
	std::function<ValArrayType(std::array<double, DIM>, std::array<size_t, DIM>)> val_interpMesh;

public:
	//获取插值点的值
	ValArrayType getPointVal(const std::vector<double> xqs_vec)
	{
		std::array<size_t, DIM> x_locs = x_locs_mark;
		std::array<bool, DIM>	need_new_loc;
		std::array<double, DIM> xqs = { 0.0 };
		std::copy(xqs_vec.begin(), xqs_vec.end(), xqs.begin());
		//动态插值位置判断
		for (size_t dim = 0; dim < DIM; dim++) {
			need_new_loc[dim] = true;
			if (m_x.at(dim).at(x_locs_mark[dim]) <= xqs[dim] && m_x.at(dim).at(x_locs_mark[dim] + 1) >= xqs[dim])
				need_new_loc[dim] = false;
		}
		//判断区间（每个维度坐标点不多，用遍历即可）
		bool out_of_range = false;
		for (size_t dim = 0; dim < DIM; dim++) {
			if (need_new_loc[dim]) {
				if (m_x.at(dim).at(0) > xqs[dim] || m_x.at(dim).at(m_size.at(dim) - 1) < xqs[dim]) {
					out_of_range = true;
					break;
				}
				for (size_t i = 1; i < m_size.at(dim); i++) {
					if (m_x.at(dim).at(i) >= xqs[dim]) {
						x_locs[dim] = i - 1;
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
		std::copy(x_locs.begin(), x_locs.end(), x_locs_mark.begin());
		return val_interpMesh(xqs, x_locs);
	}
};