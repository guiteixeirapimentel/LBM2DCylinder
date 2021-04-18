#pragma once
#include <vector>
#include <cassert>

class ScalarField2D
{
public:
	ScalarField2D(int NX, int NY, double val = 0.0)
		:
		cNX(NX),
		cNY(NY)
	{
		cData.resize(cNX * cNY, val);
	}

	double& operator()(int i, int j)
	{
		const int index = i + (j*cNX);

		assert(i >= 0 && i < cNX && j >= 0 && j < cNY);

		return reinterpret_cast<double&>(cData[index]);
	}

	double operator()(int i, int j) const
	{
		const int index = i + (j*cNX);

		assert(i >= 0 && i < cNX && j >= 0 && j < cNY);

		return cData[index];
	}

public:
	const int cNX;
	const int cNY;

	std::vector<double> cData;
};