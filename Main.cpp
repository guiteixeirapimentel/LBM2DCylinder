#include <iostream>
#include <fstream>

#include <string>

#include "FieldLBMD2Q9.h"
#include "ScalarField2D.h"

double fieq(int i, double rho, double ux, double uy);

//void setBC(FieldLBMD2Q9& )

void getMacroscopicMoments(const FieldLBMD2Q9& fi, ScalarField2D& rhoOut, ScalarField2D& uxOut, ScalarField2D& uyOut);
void getAfterColision(const ScalarField2D& rho, const ScalarField2D& ux, const ScalarField2D& uy, const FieldLBMD2Q9& fiold, FieldLBMD2Q9& fistarOut);
void stream(const FieldLBMD2Q9& fistar, FieldLBMD2Q9& finewOut);

void saveFieldToFile(const std::string& filename, const ScalarField2D& field);

void setInitialCond(FieldLBMD2Q9& fi);

constexpr double wi[] = {
	4.0 / 9.0,
	1.0 / 9.0,
	1.0 / 9.0,
	1.0 / 9.0,
	1.0 / 9.0,
	1.0 / 36.0,
	1.0 / 36.0,
	1.0 / 36.0,
	1.0 / 36.0
};

constexpr int cix[] = {
	0,
	1,
	0,
	-1,
	0,
	1,
	-1,
	-1,
	1
};

constexpr int ciy[] = {
	0,
	0,
	1,
	0,
	-1,
	1,
	1,
	-1,
	-1
};

constexpr double Lx = 0.05;
constexpr double Ly = 0.05;

constexpr int NX = 101;
constexpr int NY = 101;

constexpr double dx = Lx / (NX - 1);
constexpr double dy = Ly / (NY - 1);

constexpr double dt = dx*0.5;
constexpr double maxT = 5.0;

constexpr int maxStepsTime = ((maxT + dt) / dt);

constexpr double nu = 0.005;

constexpr double uxx = 0.5;

constexpr double Re = 1.0*uxx*Lx / nu;


constexpr double courantN = uxx * dt / dx;


constexpr double cs = 0.57735026919*dx/dt;
constexpr double cs2 = 1.0 / 3.0 * (dx*dx)/(dt*dt);
constexpr double cs4 = 1.0 / 9.0* (dx*dx) / (dt*dt)* (dx*dx) / (dt*dt);


constexpr double ratioCsUxx = uxx / cs;

constexpr double tau = nu / (cs2) + dt / 2;

FieldLBMD2Q9 fiold(NX, NY, 0.0);
FieldLBMD2Q9 finew(NX, NY, 0.0);

ScalarField2D rho(NX, NY, 0.0);
ScalarField2D ux(NX, NY, 0.0);
ScalarField2D uy(NX, NY, 0.0);

int main()
{
	double t = 0.0;
	FieldLBMD2Q9* pNewField = &finew;
	FieldLBMD2Q9* pOldField = &fiold;

	setInitialCond(*pOldField);
	setInitialCond(*pNewField);

	for (int nTime = 0; nTime < maxStepsTime; nTime++)
	{
		getMacroscopicMoments(*pOldField, rho, ux, uy);

		if (nTime % 1 == 0)
		{
			saveFieldToFile("ux/field_ux_" + std::to_string(t) + ".csv", ux);
			saveFieldToFile("uy/field_uy_" + std::to_string(t) + ".csv", uy);
			saveFieldToFile("rho/field_rho_" + std::to_string(t) + ".csv", rho);

			std::cout << t << std::endl;
		}
				

		getAfterColision(rho, ux, uy, *pOldField, *pNewField);

		stream(*pOldField, *pNewField);

		auto tmp = pOldField;
		pOldField = pNewField;
		pNewField = tmp;

		t += dt;
	}

	return 0;
}

double fieq(int i, double rho, double ux, double uy)
{
	assert(i >= 0 && i < 9);

	const double udotci = ux * cix[i] + uy * ciy[i];
	const double udotu = ux * ux + uy * uy;

	return wi[i] * rho * (1.0 + udotci / cs2 + udotci * udotci / (2.0 * cs4) - udotu / (2.0 * cs2));
}

void stream(const FieldLBMD2Q9& fistar, FieldLBMD2Q9& finewOut) {
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			for (int b = 0; b < 9; b++) {

				int nexti = i + cix[b];
				int nextj = j + ciy[b];

				if (nexti < 0) {
					nexti = NX - 1;
				}
				else if (nexti >= NX) {
					nexti = 0;
				}

				if (nextj < 0) {
					nextj = NY - 1;
				}
				else if (nextj >= NY) {
					nextj = 0;
				}
				
				finewOut(nexti, nextj, b) = fistar(i, j, b);

			}
		}
	}
}

void getMacroscopicMoments(const FieldLBMD2Q9& fi, ScalarField2D& rhoOut, ScalarField2D& uxOut, ScalarField2D& uyOut)
{
	for (int j = 0; j < NY; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			const double rho = fi(i, j, 0) + fi(i, j, 1) + fi(i, j, 2) + fi(i, j, 3) + fi(i, j, 4) + fi(i, j, 5) + fi(i, j, 6)
				+ fi(i, j, 7) + fi(i, j, 8);

			double ux = dx/dt * ((fi(i, j, 1) + fi(i, j, 5) + fi(i, j, 8)) - (fi(i, j, 3) + fi(i, j, 6) + fi(i, j, 7))) / rho;
			double uy = dx / dt * ((fi(i, j, 2) + fi(i, j, 5) + fi(i, j, 6)) - (fi(i, j, 4) + fi(i, j, 7) + fi(i, j, 8))) / rho;

			if (rho == 0.0) {
				ux = 0.0;
				uy = 0.0;
			}

			rhoOut(i, j) = rho;

			uxOut(i, j) = ux;
			uyOut(i, j) = uy;
		}
	}
}

void getAfterColision(const ScalarField2D& rho, const ScalarField2D& ux, const ScalarField2D& uy, const FieldLBMD2Q9& fiold, FieldLBMD2Q9& fistarOut)
{
	for (int j = 0; j < NY; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			for (int b = 0; b < 9; b++)
			{
				fistarOut(i, j, b) = fiold(i, j, b)*(1.0 - dt/tau) +  fieq(b, rho(i, j), ux(i, j), uy(i, j))*dt/tau;
			}
		}
	}
}

void saveFieldToFile(const std::string& filename, const ScalarField2D& field)
{
	std::ofstream arq(filename);

	for (int j = 0; j < NY; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			arq << field(i, j) << "\t";
		}
		arq << "\n";
	}

	arq.close();
}

void setInitialCond(FieldLBMD2Q9& fi)
{
	constexpr double rho0 = 1.225;
	constexpr double ux0 = 0.0;
	constexpr double uy0 = 0.0;

	// all field
	for (int j = 0; j < NY; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			for (int b = 0; b < 9; b++)
			{
				fi(i, j, b) = fieq(b, rho0, ux0, uy0);
			}			
		}
	}

	// center flow

	constexpr int top = NY / 4;
	constexpr int bot = NY / 2;

	constexpr int left = NX / 4;
	constexpr int right = NX / 2;

	for (int j = top; j < bot; j++)
	{
		for (int i = left; i < right; i++)
		{
			for (int b = 0; b < 9; b++)
			{
				fi(i, j, b) = fieq(b, rho0, uxx, uy0);
			}
		}
	}
}