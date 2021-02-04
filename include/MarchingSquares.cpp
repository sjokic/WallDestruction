#include "MarchingSquares.h"
#include <array>

void MarchingSquares::add_vertex(int x, int y, uint8_t edgeType, float isoValue, std::vector<std::tuple<Eigen::Vector3d, int, bool>>& vertices, std::unordered_map<int, int>& edgeToVertexMap) {
	const Array2d& levelSetValues = m_levelSet.x();
	int x1 = x;
	int y1 = y;
	int x2 = x + (edgeType == 0);
	int y2 = y + (edgeType == 1);
	double s1 = std::abs(levelSetValues(x1, y1) - isoValue);
	double s2 = std::abs(levelSetValues(x2, y2) - isoValue);
	double t = s1 / (s1 + s2);
	t = std::max<double>(t, double(1e-3));
	t = std::min<double>(t, double(1.0 - 1e-3));
	int resx = m_levelSet.getResolutionX();
	edgeToVertexMap[2 * (x1 + y1 * resx) + edgeType] = vertices.size();
	vertices.push_back(std::make_tuple(m_levelSet.getGridSpacing() * Eigen::Vector3d(x1 + (edgeType == 0 ? t : 0.0f), y1 + (edgeType == 1 ? t : 0.0f), 0), -1, false));
}

vector<vector<Eigen::Vector3d>> MarchingSquares::extract(float isoValue) {
	int resX = m_levelSet.getResolutionX();
	int resY = m_levelSet.getResolutionY();

	std::array<int, 4> edgeOffsets = { 0, 1, 3, 2 * resX };

	std::vector<std::tuple<Eigen::Vector3d, int, bool>> vertices;
	std::vector<std::pair<int, int>> edges;
	std::unordered_map<int, int> edgeToVertexMap;

	const Array2d& levelSetValues = m_levelSet.x();

	for (int x = 0; x < resX - 1; x++) {
		if ((levelSetValues(x, 0) >= isoValue) ^ (levelSetValues(x + 1, 0) >= isoValue)) {
			add_vertex(x, 0, 0, isoValue, vertices, edgeToVertexMap);
		}
	}

	for (int y = 0; y < resY - 1; y++) {
		if ((levelSetValues(0, y) >= isoValue) ^ (levelSetValues(0, y + 1) >= isoValue)) {
			add_vertex(0, y, 1, isoValue, vertices, edgeToVertexMap);
		}
	}

	for (int y = 0; y < resY - 1; y++) {
		for (int x = 0; x < resX - 1; x++) {
			uint8_t cubetype = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (levelSetValues(x + j, y + i) >= isoValue) {
						cubetype |= (0x1 << (i * 2 + j));
					}
				}
			}

			if (cubetype == 0x0 || cubetype == 0xf) {
				continue;
			}

			for (int i = 0; i < 2; i++) {
				if (edgeCrossingTable[cubetype] & (0x1 << i)) {
					std::array<int, 2> idx = { x + 1, y + 1 };
					idx[i]--;
					add_vertex(idx[0], idx[1], i, isoValue, vertices, edgeToVertexMap);
				}
			}

			int idx = 2 * (x + y * resX);
			for (int i = 0; edgeTable[cubetype][i] != -1; i += 2) {
				int v1 = edgeToVertexMap[idx + edgeOffsets[edgeTable[cubetype][i]]];
				int v2 = edgeToVertexMap[idx + edgeOffsets[edgeTable[cubetype][i + 1]]];
				std::get<1>(vertices[v1]) = v2;
				std::get<2>(vertices[v2]) = true;
			}
		}
	}

	vector<vector<Eigen::Vector3d>> extractedLines;
	std::vector<bool> visited(vertices.size(), false);
	for (auto& vertex : vertices) {
		if (std::get<2>(vertex)) {
			continue;
		}
		vector<Eigen::Vector3d> points;
		points.push_back(std::get<0>(vertex));
		int nextVertex = std::get<1>(vertex);
		while (nextVertex != -1) {
			visited[nextVertex] = true;
			auto nextVertexPoint = std::get<0>(vertices[nextVertex]);
			// two vertex positions can be equal when a grid vertex has
			// the value zero
			if (points.back() != nextVertexPoint) {
				points.push_back(nextVertexPoint);
			}
			nextVertex = std::get<1>(vertices[nextVertex]);
		}
		extractedLines.push_back(points);
	}

	for (auto& vertex : vertices) {
		if (!std::get<2>(vertex)) {
			continue;
		}
		if (visited[std::get<1>(vertex)]) {
			continue;
		}
		vector<Eigen::Vector3d> points;
		points.push_back(std::get<0>(vertex));
		int nextVertex = std::get<1>(vertex);
		while (!visited[nextVertex]) {
			visited[nextVertex] = true;
			auto nextVertexPoint = std::get<0>(vertices[nextVertex]);
			// two vertex positions can be equal when a grid vertex has
			// the value zero
			if (points.back() != nextVertexPoint) {
				points.push_back(nextVertexPoint);
			}
			nextVertex = std::get<1>(vertices[nextVertex]);
		}
		extractedLines.push_back(points);
	}

	return extractedLines;
}

const std::array<uint8_t, 16> MarchingSquares::edgeCrossingTable = { 0x0, 0x0, 0x2, 0x2, 0x1, 0x1, 0x3, 0x3, 0x3, 0x3, 0x1, 0x1, 0x2, 0x2, 0x0, 0x0 };

const std::array<std::array<int8_t, 5>, 16> MarchingSquares::edgeTable = { {
	{ { -1, -1, -1, -1, -1 } },
	{ { 1, 0, -1, -1, -1 } },
	{ { 0, 2, -1, -1, -1 } },
	{ { 1, 2, -1, -1, -1 } },
	{ { 3, 1, -1, -1, -1 } },
	{ { 3, 0, -1, -1, -1 } },
	{ { 0, 2, 3, 1, -1 } },
	{ { 3, 2, -1, -1, -1 } },
	{ { 2, 3, -1, -1, -1 } },
	{ { 1, 3, 2, 0, -1 } },
	{ { 0, 3, -1, -1, -1 } },
	{ { 1, 3, -1, -1, -1 } },
	{ { 2, 1, -1, -1, -1 } },
	{ { 2, 0, -1, -1, -1 } },
	{ { 0, 1, -1, -1, -1 } },
	{ { -1, -1, -1, -1, -1 } }
} };

