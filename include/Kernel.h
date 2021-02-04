#ifndef KERNEL_H
#define KERNEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>

class Kernel { // 2D
public:
	static float getRadius() { return m_h; }
	static void setRadius(float h) {
		m_h = h;
		m_h2 = h*h;

		m_poly6 = 315.0 / (64.0 * M_PI * pow(m_h, 9));
		m_spiky = -45.0 / (M_PI * pow(m_h, 6));
		m_visco = 45.0 / (M_PI * pow(m_h, 6));
		m_cubic = 40.0 / (7.0 * M_PI * m_h2);
		m_cubic_g = 240.0 / (7.0 * M_PI * m_h2);
	}

	static float W_poly6(const Eigen::Vector3d &r) {
		float r_ = r.norm();
		float res = 0;
		if (r_ <= m_h) res = m_poly6 * pow(m_h2 - r_*r_, 3);
		return res;
	}

	static Eigen::Vector3d gradW_spiky(const Eigen::Vector3d& r) {
		float r_ = r.norm();		
		Eigen::Vector3d res = Eigen::Vector3d::Zero();
		if (r_ > 0.001 && r_ <= m_h) res = m_spiky * pow(m_h - r_, 2) / r_ * r;
		return res;
	}

	static float lapW_visco(const Eigen::Vector3d& r) {
		float r_ = r.norm();
		float res = 0;
		if (r_ <= m_h) res = m_visco * (m_h - r_);
		return res;
	}

	static float W_cubic(const Eigen::Vector3d& r) {
		float r_ = r.norm();
		float q = r_ / m_h;
		float res = 0;
		if (q <= 1.0) {
			if (q <= 0.5) {
				float q2 = q * q;
				float q3 = q2 * q;
				res = m_cubic * (6.0*(q3 - q2) + 1.0);
			} 
			else {
				res = m_cubic * 2.0 * pow(1.0-q, 3);
			}
		}
		return res;
	}

	static Eigen::Vector3d gradW_cubic(const Eigen::Vector3d& r) {
		float r_ = r.norm();
		float q = r_ / m_h;
		Eigen::Vector3d res = Eigen::Vector3d::Zero();
		if (r_ > 1e-6 && q <= 1) {
			Eigen::Vector3d gradq = r / (r_*m_h);
			if (q <= 0.5) {
				res = m_cubic_g * q * (3.0* q - 2.0) * gradq;
			}
			else {
				float factor = 1.0 - q;
				res = m_cubic_g * (-factor * factor) * gradq;
			}
		}
		return res;
	}

protected:
	static float m_h;
	static float m_h2;
	static float m_poly6;
	static float m_spiky;
	static float m_visco;
	static float m_cubic;
	static float m_cubic_g;
};

#endif KERNEL_H