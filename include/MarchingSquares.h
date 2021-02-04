#ifndef _MARCHING_SQUARES_2D__
#define _MARCHING_SQUARES_2D__

#include "Grid2.h"
#include "OwnCustomAttribute.h"
#include <Eigen/Core>
#include <unordered_map>

using namespace std;


class MarchingSquares {

public:

	MarchingSquares(const Grid2& levelSet) : m_levelSet(levelSet) { };

	/* Extracts a polygon, starting with the cell indicated by initialCell */
	vector<vector<Eigen::Vector3d>> extract(float isoValue);

private:

	const Grid2 &m_levelSet;
	static const std::array<uint8_t, 16> edgeCrossingTable;
	static const std::array<std::array<int8_t, 5>, 16> edgeTable;

	void add_vertex(int i, int j, uint8_t edgeType, float isoValue, std::vector<std::tuple<Eigen::Vector3d, int, bool>>& vertices, unordered_map<int, int>& edgeToVertexMap);

};


#endif