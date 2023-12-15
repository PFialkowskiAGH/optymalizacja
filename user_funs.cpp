#include"user_funs.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
matrix* pom;

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	double& pomx = x(0, 0);
	double exponent = 0.1 * pomx - 2.0 * M_PI;
	exponent *= exponent;
	y = -cos(0.1 * pomx) * exp(-exponent) + 0.002 * pow(0.1 * pomx, 2.0);
	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(3, new double[3]{ 5, 1, 10 });
	matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	y = abs(max - 50);
	pom = Y;
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	double a = 0.98, b = 0.63, g = 9.81, PA = 1, TA = 90, PB = 1, DB = 0.00365665, Fin = 0.01, Tin = 10;
	matrix dY(3, 1);
	double FAout = Y(0) > 0 ? a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA) : 0;
	double FBout = Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;
	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA - Y(2));
	return dY;
}

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m * pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
	return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	double& x1 = x(0);
	double& x2 = x(1);
	y = pow(x1,2) + pow(x2,2) -cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;
	return y;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(2, 1), Y_ref(2, new double[2]{ 3.14,0 });
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Y_ref, x);
	int n = get_len(Y[0]);
	y = 0;
	for (int i = 0; i < n; ++i)
		y = y + 10 * pow(Y_ref(0) - Y[1](i, 0), 2) +
		pow(Y_ref(1) - Y[1](i, 1), 2) +
		pow(x(0) * (Y_ref(0) - Y[1](i, 0)) + x(1) * (Y_ref(1) - Y[1](i, 1)), 2);
	y = y * 0.1;
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	double mr = 1, mc = 10, l = 0.5, b = 0.5;
	double I = mr * l * l / 3 + mc * l * l;
	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;
	return dY;
}