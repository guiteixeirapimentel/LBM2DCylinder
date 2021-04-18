#pragma once
#include <vector>
#include <cassert>

class FieldLBMD2Q9
{
public:
	FieldLBMD2Q9(int NX, int NY, double val = 0.0)
		:
		cNX(NX),
		cNY(NY)
	{
		cData.resize(NX * NY * 9, val);
	}

	~FieldLBMD2Q9()
	{}

	double& operator()(int i, int j, int b)
	{
		const int index = (i * 9) + (j * 9 * cNX) + b;
		assert(i < cNX && i >= 0 && j < cNY && j >= 0 && b < 9 && b >= 0 );
		return reinterpret_cast<double&>(cData[index]);
	}

	double operator()(int i, int j, int b) const
	{
		const int index = (i * 9) + (j * 9 * cNX) + b;
		assert(i < cNX && i >= 0 && j < cNY && j >= 0 && b < 9 && b >= 0);
		return cData[index];
	}
public:
	const int cNX;
	const int cNY;
	std::vector<double> cData;

};