#ifndef PARTICLES_SAMPLER_H_
#define PARTICLES_SAMPLER_H_
#pragma once

#include "ParticlesData.h"
#include "Grid2.h"

class ParticlesSampler {

	public:
	typedef enum samplingMethod_t {
		uniform,
		uniformJittered,
		poisson
	} samplingMethod_t;

	typedef enum boundaryResampleStrat_t {
		clamp, //clamps particles to be inside the gridMesh
		periodic, //Re-inserts the particle on the other side of the domain
		fullResample //Fully resamples the particle according with the resampling strategy
	} boundaryResampleStrat_t;

	struct cellParticleCount_t {
		dimensions_t cellIndex;
		int count;

		cellParticleCount_t(const dimensions_t& gCellIndex, int gCount) {
			cellIndex = gCellIndex;
			count = gCount;
		}
	};

	struct CountCompareNode : public std::binary_function<cellParticleCount_t, cellParticleCount_t, bool> {
		bool operator()(const cellParticleCount_t lhs, const cellParticleCount_t rhs) const
		{
			return lhs.count > rhs.count;
		}
	};

	#pragma region Constructors
	ParticlesSampler(const Grid2 &grid, samplingMethod_t samplingMethod, uint particlesPerCell, shared_ptr<ParticlesData> pParticlesData = nullptr)
		: m_grid(grid), m_pParticlesData(pParticlesData), m_particlesPerCell(particlesPerCell), m_samplingMethod(samplingMethod) {
		initializeParticles();
	};

	ParticlesSampler(const Eigen::Vector3i& gridResolution, double cellSize, const Eigen::Vector3i & initialIndex = Eigen::Vector3i(0, 0, 0))
		: m_gridResolution(gridResolution), m_initialIndex(initialIndex), m_cellSize(cellSize), m_internalGrid(make_shared<Grid2>()), m_grid(*m_internalGrid), 
		m_particlesPerCell(1), m_samplingMethod(uniform) {
		applyParticles();
	};
	#pragma endregion

	#pragma region AccessFunctions
	inline shared_ptr<ParticlesData> getParticlesData() {
		return m_pParticlesData;
	}

	inline int getParticlesPerCell() const {
		return m_particlesPerCell;
	}
	#pragma endregion

	#pragma region Functionalities
	/** Particle datas individual velocity update*/
	//void interpolateVelocities(shared_ptr<Interpolation::Interpolant<Eigen::Vector3d, VoxelType, Eigen::Vector3d>> pInterpolant);

	/** Particle datas individual vector-based attribute update. Returns false if the attribute is not found inside particles data.*/
	//bool interpolateVelocityAttributes(const string& attributeName, shared_ptr<Interpolation::Interpolant<Eigen::Vector3d, VoxelType, Eigen::Vector3d>> pInterpolant);
	#pragma endregion

	#pragma region Initialization
	inline void initializeParticles() {
		switch (m_samplingMethod) {
			case ParticlesSampler::uniform:
				initializeUniformSampling();
			break;
			case ParticlesSampler::uniformJittered:
				initializeUniformSampling(true);
			break;
			case ParticlesSampler::poisson:
				initializePoissonSampling();
			break;
			default:
				initializeUniformSampling();
			break;
		}
	}
	void initializeUniformSampling(bool jitter = false);
	void initializePoissonSampling();
	#pragma endregion

#pragma region ApplySource
	void applyParticles(bool jitter = false);
#pragma endregion

	protected:

	#pragma region ClassMembers

	// Grid created internally if none is passed to constructor, so that m_grid can reference it
	shared_ptr<Grid2> m_internalGrid = nullptr; 

	const Grid2 &m_grid;
	Eigen::Vector3i m_gridResolution;
	Eigen::Vector3i m_initialIndex;
	double m_cellSize;

	uint m_particlesPerCell;
	shared_ptr<ParticlesData> m_pParticlesData;

	samplingMethod_t m_samplingMethod;

	//Separate boundary conditions by location for fast checking
	boundaryResampleStrat_t m_boundaryResample[6];
	#pragma endregion
};




#endif
