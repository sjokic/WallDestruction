#include "lcp_solver.h"


LCPSolver::LCPSolver(int number_of_contacts,
                     Eigen::MatrixXd arg_matA,
                     std::vector<double> arg_b)
{
	n = number_of_contacts;

	matA = Eigen::MatrixXd::Zero(n, n);
	b = Eigen::VectorXd::Zero(n);
	f = Eigen::VectorXd::Zero(n);
	df = Eigen::VectorXd::Zero(n);
	a = Eigen::VectorXd::Zero(n);
	da = Eigen::VectorXd::Zero(n);

	//prep the vectors
	for (unsigned int i = 0; i != n; ++i)
	{
		set.push_back(0);
	}

	matA = arg_matA;

	for (unsigned int i = 0; i != n; ++i)
	{
		b(i) = arg_b[i];
	}

}


void LCPSolver::calcForces()
{
	a = b;
	d = 0;

	for (unsigned int i = 0; i != n; ++i)
	{
		d = i;
		bool returnFlag = true;
		while (returnFlag)
		{
			returnFlag = driveToZero();
		}
	}
}

bool LCPSolver::driveToZero()
{
	while (true) {
		fDirection();
		da = matA * df;

		maxStep();
		f = f + s * df;
		a = a + s * da;

		if (set[j] == 1) {
			set[j] = 2;
			return true;
		}
		else if (set[j] == 2) {
			set[j] = 1;
			return true;
		}
		else {
			set[j] = 1;
			return false; //continue
		}
	}
}

void LCPSolver::fDirection()
{
	for (unsigned int i = 0; i != n; ++i)
	{
		df[i] = 0.0;
	}

	df[d] = 1.0;

	std::vector<unsigned int> indicesInC;
	for (unsigned int i = 0; i != n; ++i)
	{
		if (set[i] == 1)
		{
			indicesInC.push_back(i);
		}
	}

	unsigned int indSize = indicesInC.size();

	if (indSize == 0)
	{
		return;
	}

	//A11 = Acc
	Eigen::MatrixXd matAcc = Eigen::MatrixXd::Zero(indSize, indSize);
	for (unsigned int i = 0; i != indSize; ++i)
	{
		for (unsigned int j = 0; j != indSize; ++j)
		{
			matAcc(i, j) = matA(indicesInC[i], indicesInC[j]);
		}
	}

	//v1 = Acd
	Eigen::VectorXd vecAcd = Eigen::VectorXd::Zero(indSize);
	for (unsigned int i = 0; i != indSize; ++i)
	{
		vecAcd(i) = -matA(indicesInC[i], d);
	}

	Eigen::MatrixXd x = matAcc.ldlt().solve(vecAcd);

	for (unsigned int i = 0; i != indSize; ++i)
	{
		df(indicesInC[i]) = abs(x(i)) <= std::numeric_limits<float>::epsilon() ? 0 : x(i);
		
	}
}

void LCPSolver::maxStep()
{
	s = std::numeric_limits<double>::max();
	j = -1;

	if (da[d] > 0)
	{
		j = d;
		s = -a[d] / da[d];
	}

	for (unsigned int i = 0; i != n; ++i)
	{
		if ((set[i] == 1) && (df[i] < 0))
		{
			double temp_s = -f[i] / df[i];
			if (temp_s < s)
			{
				s = temp_s;
				j = i;
			}
		}

		if ((set[i] == 2) && (da[i] < 0))
		{
			double temp_s = -a[i] / da[i];
			if (temp_s < s)
			{
				s = temp_s;
				j = i;
			}
		}
	}
}

std::vector<double> LCPSolver::getForces()
{
    // return 0 for (nearly) zero vector, otherwise solution can't be found
	bool zeros = b.isZero();
	if (zeros) {
		return std::vector<double>(b.size());
	}

	calcForces();

	std::vector<double> temp_f;
	for (unsigned int i = 0; i != n; ++i)
	{
		temp_f.push_back(f[i]);
	}
	return temp_f;
}

