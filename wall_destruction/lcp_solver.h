#include <vector>
#include <limits>

#include <set>
#include <utility>
#include <vector>
#include "AABB.h"
#include "RigidObject.h"

// code is from Baraff: https://www.cs.cmu.edu/~baraff/papers/sig94.pdf
// and https://github.com/thesunking/mesh-test


class LCPSolver
{
private:
	int n; //total number of contact points

	std::vector<unsigned int> set;
	// set[i] == 1 -> set C, contact point
	// set[i] == 2 -> set NC, non contact point

	Eigen::MatrixXd matA; //matrix relating forces between contact points

	Eigen::VectorXd f; //forces at each contact point
	Eigen::VectorXd df; //direction of change for f
	Eigen::VectorXd a; //accelerations at each contact point
	Eigen::VectorXd da; //change in accelerations for each contact point
	Eigen::VectorXd b; //force-independent variables


	double s; //step-size
	unsigned int d; //index of contact with negative acceleration ("driving index")
	unsigned int j; //index of step-size limiting contact


	void calcForces(); //main loop for
	bool driveToZero(); //increases current contact's force until its acceleration is zero
	void fDirection(); //calculates df
	void maxStep(); //calculates maximum step size


public:
	LCPSolver(int number_of_contacts, Eigen::MatrixXd matA, std::vector<double> b);
	std::vector<double> getForces();
};