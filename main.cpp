/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
*********************************************/

#include "opt_alg.h"

void lab1();
void lab2();
//void lab3();
//void lab4();
//void lab5();
//void lab6();

int main()
{
	try
	{
		//lab1();
		lab2();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl;
	}
	system("pause");
	return 0;
}

void lab1()
{
	double x0 = m2d(200 * rand_mat() - 100), d = 1, alpha = 2.8, epsilon = 1e-5, gamma = 1e-200;
	solution opt;
	int Nmax = 1000;

	////Funkcja testowa
	//double* ab;
	//cout << x0 << endl << endl;
	//ab = expansion(ff1T, x0, d, alpha, Nmax);
	//cout << ab[0] << '\t' << ab[1] << endl << endl;
	////double a = -100;
	////double b = 100;
	//solution::clear_calls();
	//opt = fib(ff1T, ab[0], ab[1], epsilon);
	////opt = fib(ff1T, a, b, epsilon);
	//cout << opt << endl << endl;
	//solution::clear_calls();
	//opt = lag(ff1T, ab[0], ab[1], epsilon, gamma, Nmax);
	////opt = lag(ff1T, a, b, epsilon, gamma, Nmax);
	//cout << opt << endl;
	//solution::clear_calls();

	//Zbiorniki
	epsilon = 1e-10;
	opt = fib(ff1R, 1e-4, 1e-2, epsilon);
	int n = get_len(pom[0]);
	for (int i = 0; i < n; ++i)
		cout << pom[0](i, 0) << "," << pom[1](i, 0) << "," << pom[1](i, 1) << "," << pom[1](i, 2) << endl;
	cout << opt << endl << endl;
	solution::clear_calls();
	opt = lag(ff1R, 1e-4, 1e-2, epsilon, gamma, Nmax);
	n = get_len(pom[0]);
	for (int i = 0; i < n; ++i)
		cout << pom[0](i, 0) << "," << pom[1](i, 0) << "," << pom[1](i, 1) << "," << pom[1](i, 2) << endl;
	cout << opt << endl << endl;
	solution::clear_calls();
}

void lab2()
{
	//Funkcja testowa
	double s = 0.7, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-3;
	int Nmax = 1000;
	solution opt;
	matrix x0, s0;

	s0 = matrix(2, 1, s);
	x0 = 2 * rand_mat(2, 1) - 1;
	//x0(1) = 0;
	//x0(0) = 0;
	cout << x0 << endl << endl;
	opt = HJ(ff2T, x0, s, alphaHJ, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();
	opt = Rosen(ff2T, x0, s0, alphaR, beta, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Ramie robota
	/*s = 2;

	x0 = 10 * rand_mat(2, 1);
	cout << x0 << endl << endl;
	opt = HJ(ff2R, x0, s, alphaHJ, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();
	s0 = matrix(2, 1, s);
	opt = Rosen(ff2R, x0, s0, alphaR, beta, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();*/
}