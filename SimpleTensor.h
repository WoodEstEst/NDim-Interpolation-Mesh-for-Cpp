#pragma once
#include <vector>
#include <array>
#include <initializer_list>
#include <stdexcept>

//简单多值张量类型
template<class DataType, size_t DIM, size_t VALNUM>
class SimpleTensor
{
public:
	SimpleTensor(
		const std::vector<std::array<DataType, VALNUM>>& data,
		const std::array<size_t, DIM>& size) :
		m_size(size), m_data(data)
	{
		size_t diffCount = 1;
		for (size_t d = 0; d < DIM; d++) {
			size_t dim = DIM - 1 - d;
			m_diffs_peridx.at(dim) = diffCount;
			diffCount *= m_size.at(dim);
		}
		if (diffCount != m_data.size()) { 
			exit(114514);
			//throw std::runtime_error("Data size does not match tensor size.");
		}
	}
protected:
	typedef std::array<DataType, VALNUM> ValArrayType;
	std::vector<std::array<DataType, VALNUM>> m_data;
	std::array<size_t, DIM> m_size;
	std::array<size_t, DIM> m_diffs_peridx;
public:
	//方法：将多维下标数组转为展开一维数组下标
	size_t get_flat_idx(const std::array<size_t, DIM>& idxs) {
		size_t flat_idx = 0;
		for (size_t dim = 0; dim < DIM; dim++) flat_idx += idxs[dim] * m_diffs_peridx[dim];
		return flat_idx;
	}
	//方法：从多维下标获得值引用
	ValArrayType& at_idxs(const std::array<size_t, DIM>& idxs) {
		size_t flat_idx = get_flat_idx(idxs);
		return m_data.at(flat_idx);
	}
	//方法：从多维下标获得值引用(初始化列表参数)
	ValArrayType& at_idxs(const std::initializer_list<size_t>& idxs_list) {
		if (idxs_list.size() != DIM) {
			exit(114514);
		}
		std::array<size_t, DIM> idxs;
		std::copy(idxs_list.begin(), idxs_list.end(), idxs.begin());
		size_t flat_idx = get_flat_idx(idxs);
		return m_data.at(flat_idx);
	}
	//方法：从一维下标获得值引用
	ValArrayType&
		at_flat_idx(size_t flat_idx) {
		return m_data.at(flat_idx);
	}
	//获取数值数列
	std::vector<std::array<DataType, VALNUM>> getData() {
		return m_data;
	}
};