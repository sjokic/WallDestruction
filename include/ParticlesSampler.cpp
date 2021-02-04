#include "ParticlesSampler.h"

float safeRandom() {
	return (rand() / (float)RAND_MAX) - 1e-7;
}

#pragma region Functionalities
//template <class Eigen::Vector3d, template <class> class VoxelType>
//void ParticlesSamplerBase<Eigen::Vector3d, VoxelType>::interpolateVelocities(shared_ptr<Interpolation::Interpolant<Eigen::Vector3d, VoxelType, Eigen::Vector3d>> pInterpolant) {
//	for (int i = 0; i < m_pParticlesData->getVelocities().size(); i++) {
//		m_pParticlesData->getVelocities()[i] = pInterpolant->interpolate(m_pParticlesData->getPositions()[i]);
//	}
//}
//
//template <class Eigen::Vector3d, template <class> class VoxelType>
//bool ParticlesSamplerBase<Eigen::Vector3d, VoxelType>::interpolateVelocityAttributes(const string& attributeName, shared_ptr<Interpolation::Interpolant<Eigen::Vector3d, VoxelType, Eigen::Vector3d>> pInterpolant) {
//	if (!m_pParticlesData->hasVectorBasedAttribute(attributeName))
//		return false;
//
//	vector<Eigen::Vector3d>& vectorAttributes = m_pParticlesData->getVectorBasedAttribute(attributeName);
//	for (int i = 0; i < vectorAttributes.size(); i++) {
//		vectorAttributes[i] = pInterpolant->interpolate(m_pParticlesData->getPositions()[i]);
//	}
//
//	return true;
//}


//shared_ptr<ParticlesData<Eigen::Vector3d>> createSampledParticlesUniform(shared_ptr<QuadGridMesh<Eigen::Vector3d>> pQuadGridMesh, uint particlesPerCell) {
//	int totalNumberOfCells = (pQuadGridMesh->getGridDimensions().x - 2) * (pQuadGridMesh->getGridDimensions().y - 2);
//	shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData(new ParticlesData<Eigen::Vector3d>(totalNumberOfCells * particlesPerCell));
//
//	float dx = pQuadGridMesh->getGridSpacing();
//
//	vector<Eigen::Vector3d>& particlesPositions = pParticlesData->getPositions();
//	vector<Eigen::Vector3d>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//	for (int i = 1; i < pQuadGridMesh->getGridDimensions().x - 1; i++) {
//		for (int j = 1; j < pQuadGridMesh->getGridDimensions().y - 1; j++) {
//			for (int l = 0; l < particlesPerCell; l++) {
//				Vector2 particlePosition = Vector2(i + safeRandom(), j + safeRandom());
//				particlesPositions.push_back(particlePosition * dx);
//				particlesVelocities.push_back(Vector2());
//				particlesResampled.push_back(false);
//			}
//		}
//	}
//
//	return pParticlesData;
//}
//
#pragma endregion

#pragma region InitializationFunctions
void ParticlesSampler::initializeUniformSampling(bool jitter) {
	uint totalNumberOfCells = m_grid.getResolutionX()* m_grid.getResolutionY();
	float dx = m_grid.getGridSpacing();
	m_pParticlesData = make_shared<ParticlesData>(totalNumberOfCells * m_particlesPerCell);
	
	vector<Eigen::Vector3d>& particlesPositions = m_pParticlesData->getPositions();
	vector<Eigen::Vector3d>& particlesVelocities = m_pParticlesData->getVelocities();
	vector<bool>& particlesResampled = m_pParticlesData->getResampledParticles();

	for (uint i = 0; i < m_grid.getResolutionX(); i++) {
		for (uint j = 0; j < m_grid.getResolutionY(); j++) {
			Eigen::Vector3d particlePosition = Eigen::Vector3d(i + 0.5, j + 0.5, 0);
			if (jitter) {
				particlePosition.x() = particlePosition.x() + safeRandom() * 0.5;
				particlePosition.y() = particlePosition.y() + safeRandom() * 0.5;
			}
			particlesPositions.push_back(particlePosition * dx);
			particlesVelocities.push_back(Eigen::Vector3d::Zero());
			particlesResampled.push_back(false);
		}
	}
}
void ParticlesSampler::initializePoissonSampling() {
	uint totalNumberOfCells = m_grid.getResolutionX() * m_grid.getResolutionY();
	float dx = m_grid.getGridSpacing();
	m_pParticlesData = make_shared<ParticlesData>(totalNumberOfCells * m_particlesPerCell);

	vector<Eigen::Vector3d>& particlesPositions = m_pParticlesData->getPositions();
	vector<Eigen::Vector3d>& particlesVelocities = m_pParticlesData->getVelocities();
	vector<bool>& particlesResampled = m_pParticlesData->getResampledParticles();

	for (uint i = 0; i < m_grid.getResolutionX(); i++) {
		for (uint j = 0; j < m_grid.getResolutionY(); j++) {
			Eigen::Vector3d particlePosition = Eigen::Vector3d(i + 0.5, j + 0.5, 0);
			
			particlesPositions.push_back(particlePosition * dx);
			particlesVelocities.push_back(Eigen::Vector3d::Zero());
			particlesResampled.push_back(false);
		}
	}
}
#pragma endregion

#pragma region ApplySource
void ParticlesSampler::applyParticles(bool jitter) {
	uint totalNumberOfCells = m_gridResolution.x() * m_gridResolution.y() * m_gridResolution.z();
	m_pParticlesData = make_shared<ParticlesData>(totalNumberOfCells);

	vector<Eigen::Vector3d>& particlesPositions = m_pParticlesData->getPositions();
	vector<Eigen::Vector3d>& particlesVelocities = m_pParticlesData->getVelocities();
	vector<bool>& particlesResampled = m_pParticlesData->getResampledParticles();

	for (uint i = m_initialIndex.x(); i < m_gridResolution.x(); i++) {
		for (uint j = m_initialIndex.y(); j < m_gridResolution.y(); j++) {
			for (uint k = m_initialIndex.z(); k < m_gridResolution.z(); k++) {
				Eigen::Vector3d particlePosition = Eigen::Vector3d(i + 0.5, j + 0.5, k + 0.5);
				if (jitter) {
					particlePosition.x() += safeRandom() * 0.5;
					particlePosition.y() += safeRandom() * 0.5;
					particlePosition.z() += safeRandom() * 0.5;
				}

				particlesPositions.push_back(particlePosition * m_cellSize);
				particlesVelocities.push_back(Eigen::Vector3d::Zero());
				particlesResampled.push_back(false);
			}
		}
	}
}
#pragma endregion

//template <class Eigen::Vector3d>
//shared_ptr<ParticlesData<Eigen::Vector3d>> createSampledParticlesUniform(shared_ptr<HexaGridMesh<Eigen::Vector3d>> pHexaGridMesh, uint particlesPerCell) {
//	int totalNumberOfCells = (pHexaGridMesh->getGridDimensions().x - 2) * (pHexaGridMesh->getGridDimensions().y - 2) * (pHexaGridMesh->getGridDimensions().z - 2);
//	shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData(new ParticlesData<Eigen::Vector3d>(totalNumberOfCells * particlesPerCell));
//
//	float dx = pHexaGridMesh->getGridSpacing();
//
//	vector<Eigen::Vector3d>& particlesPositions = pParticlesData->getPositions();
//	vector<Eigen::Vector3d>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//	for (int i = 1; i < pHexaGridMesh->getGridDimensions().x - 1; i++) {
//		for (int j = 1; j < pHexaGridMesh->getGridDimensions().y - 1; j++) {
//			for (int k = 1; k < pHexaGridMesh->getGridDimensions().z - 1; k++) {
//				for (int l = 0; l < particlesPerCell; l++) {
//					Vector3 particlePosition = Vector3(i + safeRandom(), j + safeRandom(), k + safeRandom());
//					particlesPositions.push_back(particlePosition * dx);
//					particlesVelocities.push_back(Vector3(0, 0, 0));
//
//					particlesResampled.push_back(false);
//				}
//			}
//		}
//	}
//}
//
//template <class Eigen::Vector3d>
//shared_ptr<ParticlesData<Eigen::Vector3d>> createSampledParticlesLocalUniform(shared_ptr<QuadGridMesh<Eigen::Vector3d>> pQuadGridMesh, uint particlesPerCell) {
//	int totalNumberOfCells = (pQuadGridMesh->getGridDimensions().x - 2) * (pQuadGridMesh->getGridDimensions().y - 2);
//	shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData(new ParticlesData<Eigen::Vector3d>(totalNumberOfCells * particlesPerCell));
//
//	float dx = pQuadGridMesh->getGridSpacing();
//	//This will subdivide the cells according with number of particles per cell to provide better sampling
//	int refactorScale = floor(sqrt(particlesPerCell));
//
//	vector<Eigen::Vector3d>& particlesPositions = pParticlesData->getPositions();
//	vector<Eigen::Vector3d>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//	dimensions_t newGridDimensions(pQuadGridMesh->getGridDimensions().x * refactorScale, pQuadGridMesh->getGridDimensions().y * refactorScale);
//	for (int i = refactorScale; i < newGridDimensions.x - refactorScale; i++) {
//		for (int j = refactorScale; j < newGridDimensions.y - refactorScale; j++) {
//			Eigen::Vector3d particlePosition = Eigen::Vector3d(i + safeRandom(), j + safeRandom());
//			particlePosition /= refactorScale;
//			particlesPositions.push_back(particlePosition * dx);
//			particlesVelocities.push_back(Eigen::Vector3d());
//
//			particlesResampled.push_back(false);
//		}
//	}
//
//	return pParticlesData;
//}
//
//template <class Eigen::Vector3d>
//shared_ptr<ParticlesData<Eigen::Vector3d>> createSampledParticlesLocalUniform(shared_ptr<HexaGridMesh<Eigen::Vector3d>> pHexaGridMesh, uint particlesPerCell) {
//	int totalNumberOfCells = (pHexaGridMesh->getGridDimensions().x - 2) * (pHexaGridMesh->getGridDimensions().y - 2) * (pHexaGridMesh->getGridDimensions().z - 2);
//	shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData(new ParticlesData<Eigen::Vector3d>(totalNumberOfCells * particlesPerCell));
//
//	float dx = pHexaGridMesh->getGridSpacing();
//	//This will subdivide the cells according with number of particles per cell to provide better sampling
//	int refactorScale = floor(cbrt(particlesPerCell));
//
//	vector<Eigen::Vector3d>& particlesPositions = pParticlesData->getPositions();
//	vector<Eigen::Vector3d>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//	dimensions_t newGridDimensions(pHexaGridMesh->getGridDimensions().x * refactorScale,
//		pHexaGridMesh->getGridDimensions().y * refactorScale,
//		pHexaGridMesh->getGridDimensions().z * refactorScale);
//
//	for (int i = refactorScale; i < newGridDimensions.x - refactorScale; i++) {
//		for (int j = refactorScale; j < newGridDimensions.y - refactorScale; j++) {
//			for (int k = refactorScale; k < newGridDimensions.z - refactorScale; k++) {
//				Vector3 particlePosition = Vector3(i + safeRandom(), j + safeRandom(), k + safeRandom());
//				particlePosition /= refactorScale;
//				particlesPositions.push_back(particlePosition * dx);
//				particlesVelocities.push_back(Vector3(0, 0, 0));
//
//				particlesResampled.push_back(false);
//			}
//		}
//	}
//
//	return pParticlesData;
//}
//
//
//template <class Eigen::Vector3d>
//shared_ptr<ParticlesData<Eigen::Vector3d>> createSampledParticlesPoisson(shared_ptr<QuadGridMesh<Eigen::Vector3d>> pQuadGridMesh, uint particlesPerCell) {
//	int totalNumberOfCells = (pQuadGridMesh->getGridDimensions().x - 2) * (pQuadGridMesh->getGridDimensions().y - 2);
//	shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData(new ParticlesData<Eigen::Vector3d>(totalNumberOfCells * particlesPerCell));
//
//	float dx = pQuadGridMesh->getGridSpacing();
//	//This will subdivide the cells according with number of particles per cell to provide better sampling
//	int refactorScale = floor(sqrt(particlesPerCell));
//
//	vector<Eigen::Vector3d>& particlesPositions = pParticlesData->getPositions();
//	vector<Eigen::Vector3d>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//	float radius = (dx / refactorScale);
//	Eigen::Vector3d initialBoundary = pQuadGridMesh->getBounds().first + Eigen::Vector3d(1, 1) * dx;
//	Eigen::Vector3d finalBoundary = pQuadGridMesh->getBounds().second - Eigen::Vector3d(1, 1) * dx;
//	dimensions_t newGridDimensions((finalBoundary.x - initialBoundary.x) / radius, (finalBoundary.y - initialBoundary.y) / radius);
//
//	Array2D<int> particleIndices(newGridDimensions);
//	particleIndices.assign(-1);
//	for (int i = refactorScale; i < newGridDimensions.x - refactorScale; i++) {
//		for (int j = refactorScale; j < newGridDimensions.y - refactorScale; j++) {
//			bool closeToOtherParticles = true;
//			Vector2 particlePos;
//			if (closeToOtherParticles) {
//				int k = 0;
//				do
//				{
//					particlePos = Eigen::Vector3d(i + (rand() / (float)RAND_MAX), j + (rand() / (float)RAND_MAX)) * (dx / refactorScale);
//
//					if (i == refactorScale || i == newGridDimensions.x - refactorScale - 1) {
//						closeToOtherParticles = false; continue; k++;
//					}
//					else if (j == refactorScale || j == newGridDimensions.y - refactorScale - 1) {
//						closeToOtherParticles = false; continue; k++;
//					}
//
//					if (particleIndices(i - 1, j) != -1) {
//						Eigen::Vector3d leftParticlePosition = particlesPositions[particleIndices(i - 1, j)];
//						float neighborParticleDistance = (particlePos - leftParticlePosition).length();
//						if (neighborParticleDistance < radius) {
//							k++;
//							continue;
//						}
//					}
//					if (particleIndices(i, j - 1) != -1) {
//						Eigen::Vector3d bottomParticlePosition = particlesPositions[particleIndices(i, j - 1)];
//						float neighborParticleDistance = (particlePos - bottomParticlePosition).length();
//						if (neighborParticleDistance < radius) {
//							k++;
//							continue;
//						}
//					}
//					if (particleIndices(i - 1, j - 1) != -1) {
//						Eigen::Vector3d bottomLeftParticlePosition = particlesPositions[particleIndices(i - 1, j - 1)];
//						float neighborParticleDistance = (particlePos - bottomLeftParticlePosition).length();
//						if (neighborParticleDistance < radius) {
//							k++;
//							continue;
//						}
//					}
//
//					closeToOtherParticles = false;
//					k++;
//				} while (closeToOtherParticles && k < 30);
//			}
//
//			if (!closeToOtherParticles) {
//				particleIndices(i, j) = particlesPositions.size();
//
//				particlesPositions.push_back(particlePos);
//				particlesResampled.push_back(false);
//				particlesVelocities.push_back(Vector2(0, 0));
//			}
//
//		}
//	}
//
//	return pParticlesData;
//}
//
//template <class Eigen::Vector3d>
//shared_ptr<ParticlesData<Eigen::Vector3d>> createSampledParticlesPoisson(shared_ptr<HexaGridMesh<Eigen::Vector3d>> pHexaGridMesh, uint particlesPerCell) {
//	int totalNumberOfCells = (pHexaGridMesh->getGridDimensions().x - 2) * (pHexaGridMesh->getGridDimensions().y - 2) * (pHexaGridMesh->getGridDimensions().z - 2);
//	shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData(new ParticlesData<Eigen::Vector3d>(totalNumberOfCells * particlesPerCell));
//
//	float dx = pHexaGridMesh->getGridSpacing();
//	//This will subdivide the cells according with number of particles per cell to provide better sampling
//	int refactorScale = floor(cbrt(particlesPerCell));
//
//	vector<Eigen::Vector3d>& particlesPositions = pParticlesData->getPositions();
//	vector<Eigen::Vector3d>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//	float radius = (dx / refactorScale);
//	Eigen::Vector3d initialBoundary = pHexaGridMesh->getBounds().first + Eigen::Vector3d(1, 1, 1) * dx;
//	Eigen::Vector3d finalBoundary = pHexaGridMesh->getBounds().second - Eigen::Vector3d(1, 1, 1) * dx;
//	dimensions_t newGridDimensions((finalBoundary.x - initialBoundary.x) / radius,
//		(finalBoundary.y - initialBoundary.y) / radius,
//		(finalBoundary.z - initialBoundary.z) / radius);
//
//	Array3D<int> particleIndices(newGridDimensions);
//	particleIndices.assign(-1);
//	for (int i = refactorScale; i < newGridDimensions.x - refactorScale; i++) {
//		for (int j = refactorScale; j < newGridDimensions.y - refactorScale; j++) {
//			for (int k = refactorScale; k < newGridDimensions.z - refactorScale; k++) {
//				bool closeToOtherParticles = true;
//				Eigen::Vector3d particlePos;
//				if (closeToOtherParticles) {
//					int l = 0;
//					do
//					{
//						particlePos = Eigen::Vector3d(i + (rand() / (float)RAND_MAX), j + (rand() / (float)RAND_MAX), k + (rand() / (float)RAND_MAX)) * (dx / refactorScale);
//
//						if (i == refactorScale || i == newGridDimensions.x - refactorScale - 1) {
//							closeToOtherParticles = false; continue; k++;
//						}
//						else if (j == refactorScale || j == newGridDimensions.y - refactorScale - 1) {
//							closeToOtherParticles = false; continue; k++;
//						}
//
//						if (particleIndices(i - 1, j, k) != -1) {
//							Eigen::Vector3d leftParticlePosition = particlesPositions[particleIndices(i - 1, j, k)];
//							float neighborParticleDistance = (particlePos - leftParticlePosition).length();
//							if (neighborParticleDistance < radius) {
//								l++;
//								continue;
//							}
//						}
//						if (particleIndices(i, j - 1, k) != -1) {
//							Eigen::Vector3d bottomParticlePosition = particlesPositions[particleIndices(i, j - 1, k)];
//							float neighborParticleDistance = (particlePos - bottomParticlePosition).length();
//							if (neighborParticleDistance < radius) {
//								l++;
//								continue;
//							}
//						}
//						if (particleIndices(i - 1, j - 1, k) != -1) {
//							Eigen::Vector3d bottomLeftParticlePosition = particlesPositions[particleIndices(i - 1, j - 1, k)];
//							float neighborParticleDistance = (particlePos - bottomLeftParticlePosition).length();
//							if (neighborParticleDistance < radius) {
//								l++;
//								continue;
//							}
//						}
//
//						//k + 1 index
//						if (particleIndices(i - 1, j, k + 1) != -1) {
//							Eigen::Vector3d leftParticlePosition = particlesPositions[particleIndices(i - 1, j, k + 1)];
//							float neighborParticleDistance = (particlePos - leftParticlePosition).length();
//							if (neighborParticleDistance < radius) {
//								l++;
//								continue;
//							}
//						}
//						if (particleIndices(i, j - 1, k + 1) != -1) {
//							Eigen::Vector3d bottomParticlePosition = particlesPositions[particleIndices(i, j - 1, k + 1)];
//							float neighborParticleDistance = (particlePos - bottomParticlePosition).length();
//							if (neighborParticleDistance < radius) {
//								l++;
//								continue;
//							}
//						}
//						if (particleIndices(i - 1, j - 1, k + 1) != -1) {
//							Eigen::Vector3d bottomLeftParticlePosition = particlesPositions[particleIndices(i - 1, j - 1, k + 1)];
//							float neighborParticleDistance = (particlePos - bottomLeftParticlePosition).length();
//							if (neighborParticleDistance < radius) {
//								l++;
//								continue;
//							}
//						}
//
//
//
//						closeToOtherParticles = false;
//						l++;
//					} while (closeToOtherParticles && l < 30);
//				}
//
//				if (!closeToOtherParticles) {
//					particleIndices(i, j, k) = particlesPositions.size();
//
//					particlesPositions.push_back(particlePos);
//					particlesResampled.push_back(false);
//					particlesVelocities.push_back(Eigen::Vector3d(0, 0, 0));
//				}
//			}
//		}
//	}
//
//	return pParticlesData;
//}
//
//template <class Eigen::Vector3d>
//void resampleParticles(shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData, shared_ptr<QuadGridMesh<Eigen::Vector3d>> pQuadGridMesh, uint particlesPerCell) {
//	using cellParticleCount_t = typename ParticlesSamplerBase<Eigen::Vector3d>::cellParticleCount_t;
//	using CountCompareNode = typename ParticlesSamplerBase<Eigen::Vector3d>::CountCompareNode;
//	vector<Vector2>& particlesPositions = pParticlesData->getPositions();
//	vector<Vector2>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//	float dx = pQuadGridMesh->getGridSpacing();
//	Array2D<int> m_particlesCount(pQuadGridMesh->getGridSpacing());
//
//	m_particlesCount.assign(0);
//	for (int i = 0; i < particlesPositions.size(); i++) {
//		particlesResampled[i] = false;
//		boundaryResample(i, pParticlesData);
//
//		int indexI = particlesPositions[i].x / dx, indexJ = particlesPositions[i].y / dx;
//		m_particlesCount(indexI, indexJ) += 1;
//	}
//
//	priority_queue<cellParticleCount_t, vector<cellParticleCount_t>, CountCompareNode> pq;
//	for (int i = 1; i < m_particlesCount.getDimensions().x - 1; i++) {
//		for (int j = 1; j < m_particlesCount.getDimensions().y - 1; j++) {
//			pq.push(cellParticleCount_t(dimensions_t(i, j), m_particlesCount(i, j)));
//		}
//	}
//
//	for (int i = 0; i < particlesPositions.size(); i++) {
//		int maxParticles = ceil(particlesPerCell * 1.05f);
//		Vector2 gridSpacePosition = particlesPositions[i] / dx;
//		if (m_particlesCount(gridSpacePosition.x, gridSpacePosition.y) > maxParticles) {
//			cellParticleCount_t top = pq.top();
//			m_particlesCount(gridSpacePosition.x, gridSpacePosition.y) -= 1;
//			pq.pop();
//			particlesPositions[i] = Vector2(top.cellIndex.x + safeRandom(), top.cellIndex.y + safeRandom()) * dx;
//			particlesResampled[i] = true;
//			top.count++;
//			pq.push(top);
//		}
//	}
//}
//
//template <class Eigen::Vector3d>
//void resampleParticles(shared_ptr<ParticlesData<Eigen::Vector3d>> pParticlesData, shared_ptr<HexaGridMesh<Eigen::Vector3d>> pHexaGridMesh, uint particlesPerCell) {
//	using cellParticleCount_t = typename ParticlesSamplerBase<Eigen::Vector3d>::cellParticleCount_t;
//	using CountCompareNode = typename ParticlesSamplerBase<Eigen::Vector3d>::CountCompareNode;
//
//	vector<Eigen::Vector3d>& particlesPositions = pParticlesData->getPositions();
//	vector<Eigen::Vector3d>& particlesVelocities = pParticlesData->getVelocities();
//	vector<bool>& particlesResampled = pParticlesData->getResampledParticles();
//
//
//	float dx = pHexaGridMesh->getGridSpacing();
//	Array3D<int> m_particlesCount(pHexaGridMesh->getGridSpacing());
//
//	m_particlesCount.assign(0);
//	for (int i = 0; i < particlesPositions.size(); i++) {
//		particlesResampled[i] = false;
//		boundaryResample(i, pParticlesData);
//
//		int indexI = particlesPositions[i].x / dx, indexJ = particlesPositions[i].y / dx, indexK = particlesPositions[i].z / dx;
//		m_particlesCount(indexI, indexJ, indexK) += 1;
//	}
//
//	priority_queue<cellParticleCount_t, vector<cellParticleCount_t>, CountCompareNode> pq;
//	for (int i = 1; i < m_particlesCount.getDimensions().x - 1; i++) {
//		for (int j = 1; j < m_particlesCount.getDimensions().y - 1; j++) {
//			for (int k = 1; k < m_particlesCount.getDimensions().z - 1; k++) {
//				pq.push(cellParticleCount_t(dimensions_t(i, j, k), m_particlesCount(i, j, k)));
//			}
//		}
//	}
//
//	for (int i = 0; i < particlesPositions.size(); i++) {
//		int maxParticles = ceil(particlesPerCell * 1.05f);
//		Eigen::Vector3d gridSpacePosition = particlesPositions[i] / dx;
//		int currParCount = m_particlesCount(gridSpacePosition.x, gridSpacePosition.y, gridSpacePosition.z);
//		if (m_particlesCount(gridSpacePosition.x, gridSpacePosition.y, gridSpacePosition.z) > maxParticles) {
//			cellParticleCount_t top = pq.top();
//			m_particlesCount(gridSpacePosition.x, gridSpacePosition.y, gridSpacePosition.z) -= 1;
//			pq.pop();
//			particlesPositions[i] = Eigen::Vector3d(top.cellIndex.x + safeRandom(),
//				top.cellIndex.y + safeRandom(),
//				top.cellIndex.z + safeRandom()) * dx;
//			particlesResampled[i] = true;
//			top.count++;
//			pq.push(top);
//		}
//	}
//}



