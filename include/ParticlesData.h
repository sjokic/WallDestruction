
#ifndef PARTICLES_DATA_H_
#define PARTICLES_DATA_H_
#pragma once

#include "OwnCustomAttribute.h"
#include <Eigen/Core>

using namespace std;


class ParticlesData : public OwnCustomAttribute<vector<Eigen::Vector3d>>, public OwnCustomAttribute <vector<float>>,
	public OwnCustomAttribute<vector<int>> {

public:
	//Will only reserve sizes inside vectors, explicit position/velocity initialization is on user-side
	ParticlesData(int initialNumParticles) {
		m_numInitialParticles = initialNumParticles;
		m_positions.reserve(initialNumParticles);
		m_velocities.reserve(initialNumParticles);
	}

	#pragma region Acceess Functions
	/** Simple "getters and setters"*/
	int getNumParticles() const {
		return m_positions.size();
	}

	const vector<Eigen::Vector3d>& getPositions() const {
		return m_positions;
	}

	vector<Eigen::Vector3d>& getPositions() {
		return m_positions;
	}

	const vector<Eigen::Vector3d>& getVelocities() const {
		return m_velocities;
	}

	vector<Eigen::Vector3d>& getVelocities() {
		return m_velocities;
	}

	const vector<bool>& getResampledParticles() const {
		return m_resampledParticles;
	}

	vector<bool>& getResampledParticles() {
		return m_resampledParticles;
	}

	/** Adds a particle: initializes the position and adds placeholders for all custom allocated variables */
	void addParticle(const Eigen::Vector3d& position) {
		m_positions.push_back(position);
		m_velocities.push_back(Eigen::Vector3d());
		//This is a particle that is initally sampled
		m_resampledParticles.push_back(true);


		auto vectorBasedAttributes = OwnCustomAttribute<vector<Eigen::Vector3d>>::get(this)->getAttributesMap();
		for (auto iter = vectorBasedAttributes.begin(); iter != vectorBasedAttributes.end(); iter++) {
			iter->second.push_back(Eigen::Vector3d());
		}

		auto scalarBasedAttributes = OwnCustomAttribute<vector<float>>::get(this)->getAttributesMap();
		for (auto iter = scalarBasedAttributes.begin(); iter != scalarBasedAttributes.end(); iter++) {
			iter->second.push_back(0.0f);
		}

		auto intBasedAttributes = OwnCustomAttribute<vector<int>>::get(this)->getAttributesMap();
		for (auto iter = intBasedAttributes.begin(); iter != intBasedAttributes.end(); iter++) {
			iter->second.push_back(0);
		}
	}

	/** Use this function for resampling, since it resets all interior particle fields */
	void resampleParticle(uint particleIndex, const Eigen::Vector3d& position) {
		m_positions[particleIndex] = position;
		Eigen::Vector3d zeroVelocity;
		//Zero out velocity
		m_velocities[particleIndex] = zeroVelocity;

		auto vectorBasedAttributes = OwnCustomAttribute<vector<Eigen::Vector3d>>::get(this)->getAttributesMap();
		for (auto iter = vectorBasedAttributes.begin(); iter != vectorBasedAttributes.end(); iter++) {
			iter->second[particleIndex] = zeroVelocity;
		}

		auto scalarBasedAttributes = OwnCustomAttribute<vector<float>>::get(this)->getAttributesMap();
		for (auto iter = scalarBasedAttributes.begin(); iter != scalarBasedAttributes.end(); iter++) {
			iter->second[particleIndex] = 0.f;
		}

		auto intBasedAttributes = OwnCustomAttribute<vector<int>>::get(this)->getAttributesMap();
		for (auto iter = intBasedAttributes.begin(); iter != intBasedAttributes.end(); iter++) {
			iter->second[particleIndex] = 0;
		}

		m_resampledParticles[particleIndex] = true;
	}

	#pragma endregion
	protected:
		int m_numInitialParticles;
		/* Standard particles properties*/
		vector<Eigen::Vector3d> m_positions;
		vector<Eigen::Vector3d> m_velocities;
		/** Is important to tag if the particles were resampled, for both rendering & simulation purposes */
		vector<bool> m_resampledParticles;
};

#endif
