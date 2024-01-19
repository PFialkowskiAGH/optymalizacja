#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		int i = 0;
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X0.y == X1.y)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X0.y < X1.y)
		{
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y)
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x - d);
				return p;
			}
		}
		solution X2;
		while (true)
		{
			++i;
			X2.x = X0.x + pow(alpha, i) * d;
			X2.fit_fun(ff, ud1, ud2);
			if (X2.y > X1.y || solution::f_calls > Nmax)
				break;
			X0 = X1;
			X1 = X2;
		}
		if (d > 0)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X2.x);
		}
		else {
			p[0] = m2d(X2.x);
			p[1] = m2d(X0.x);
		}
		cout << "f_calls: " << solution::f_calls << endl;
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.ud = b - a;
		int n = (int)(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
		int* F = new int[n] {1, 1};
		for (int i = 2; i < n; ++i)
			F[i] = F[i - 2] + F[i - 1];
		solution A(a), B(b), C, D;
		C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x; //mocne za³o¿enie
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);
		for (int i = 0; i <= n - 3; ++i)
		{
			//cout << "a: " << A.x << endl;
			//cout << "b: " << B.x << endl;
			if (C.y < D.y)
				B.x(0, 0) = D.x(0, 0);
			else
				A.x(0, 0) = C.x(0, 0);
			C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
			D.x = A.x + B.x - C.x;
			C.fit_fun(ff, ud1, ud2);
			D.fit_fun(ff, ud1, ud2);

			Xopt.ud.add_row((B.x - A.x)());
		}
		Xopt = C;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.ud = b - a;
		solution A(a), B(b), C, D, D_old(a);
		C.x = (b + a) / 2;
		double l, m;
		while (true)
		{
			//cout << "a: " << A.x << endl;
			//cout << "b: " << B.x << endl;
			A.fit_fun(ff, ud1, ud2);
			B.fit_fun(ff, ud1, ud2);
			C.fit_fun(ff, ud1, ud2);
			l = m2d(A.y * (pow(B.x, 2) - pow(C.x, 2)) + B.y * (pow(C.x, 2) - pow(A.x, 2)) + C.y * (pow(A.x, 2) - pow(B.x, 2)));
			m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x));
			if (m <= 0)
			{
				Xopt = D_old;
				Xopt.flag = 2;
				return Xopt;
			}
			D.x = 0.5 * l / m;
			D.fit_fun(ff, ud1, ud2);
			if (A.x < D.x && D.x < C.x)
			{
				if (D.y < C.y)
				{
					B.x = C.x;
					C.x = D.x;
				}
				else {
					A.x = D.x;
				}
			}
			else if (C.x < D.x && D.x < B.x)
			{
				if (D.y < C.y)
				{
					A.x = C.x;
					C.x = D.x;
				}
				else {
					B.x = D.x;
				}
			}
			else
			{
				Xopt = D_old;
				Xopt.flag = 2;
				return Xopt;
			}

			Xopt.ud.add_row((B.x - A.x)());

			if (B.x - A.x < epsilon || abs(m2d(D.x) - m2d(D_old.x)) < gamma)
			{
				Xopt = D;
				Xopt.flag = 0;
				return Xopt;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = D;
				Xopt.flag = 1;
				break;
			}
			D_old = D;
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.ud = trans(x0);
		solution XB(x0), XB_old, X;
		XB.fit_fun(ff, ud1, ud2);
		//int i = 1;
		while (true)
		{
			/*cout << "it: " << i << endl;
			i++;
			cout << "XB: " << XB.x << endl;*/
			X = HJ_trial(ff, XB, s, ud1, ud2);
			if (X.y < XB.y)
			{
				while (true)
				{
					XB_old = XB;
					XB = X;

					Xopt.ud.add_row(trans(XB.x));

					X.x = 2 * XB.x - XB_old.x;
					X.fit_fun(ff, ud1, ud2);
					X = HJ_trial(ff, X, s, ud1, ud2);
					if (X.y >= XB.y)
						break;
					if (solution::f_calls > Nmax)
					{
						Xopt = XB;
						Xopt.flag = 0;
						return Xopt;
					}
				}
			}
			else
				s *= alpha;
			if (s < epsilon)
			{
				Xopt = XB;
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = XB;
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_dim(XB);
		matrix D = ident_mat(n);
		solution X;
		for (int i = 0; i < n; ++i)
		{
			X.x = XB.x + s * D[i];
			X.fit_fun(ff, ud1, ud2);
			if (X.y < XB.y)
				XB = X;
			else
			{
				X.x = XB.x - s * D[i];
				X.fit_fun(ff, ud1, ud2);
				if (X.y < XB.y)
					XB = X;
			}
		}
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.ud = trans(x0);
		solution XB(x0), X;
		int n = get_dim(XB);
		matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
		XB.fit_fun(ff, ud1, ud2);
		//int i = 1;
		while (true)
		{
			/*cout << "it: " << i << endl;
			i++;
			cout << "XB: " << XB.x << endl;*/
			for (int i = 0; i < n; ++i)
			{
				X.x = XB.x + s(i) * D[i];
				X.fit_fun(ff, ud1, ud2);
				if (XB.y > X.y)
				{
					XB = X;
					l(i) += s(i);
					s(i) *= alpha;
				}
				else
				{
					s(i) *= -beta;
					++p(i);
				}
			}

			Xopt.ud.add_row(trans(XB.x));

			bool change = true;
			for (int i = 0; i < n; ++i)
				if (p(i) == 0 || l(i) == 0)
				{
					change = false;
					break;
				}
			if (change)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = l(i);
				Q = D * Q;
				v = Q[0] / norm(Q[0]);
				D.set_col(v, 0);
				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; ++j)
						temp = temp + trans(Q[i]) * D[j] * D[j];
					v = (Q[i] - temp) / norm(Q[i] - temp);
					D.set_col(v, i);
				}
				s = s0;
				l = matrix(n, 1);
				p = matrix(n, 1);
			}
			double max_s = abs(s(0));
			for (int i = 1; i < n; ++i)
				if (max_s < abs(s(i)))
					max_s = abs(s(i));
			if (max_s < epsilon)
			{
				Xopt = XB;
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = XB;
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
		solution X(x0), X1;
		while (true)
		{
			X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
			if (norm(X1.x - X.x) < epsilon)
			{
				Xopt = X1;
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = X1;
				Xopt.flag = 0;
				break;
			}
			c = dc * c;
			X = X1;
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//pc=p__, S=p, D=e
		solution Xopt;
		int n = get_len(x0);
		matrix D = ident_mat(n);
		int N = n + 1;
		solution* S = new solution[N];
		S[0].x = x0;
		S[0].fit_fun(ff, ud1, ud2);
		for (int i = 1; i < N; ++i)
		{
			S[i].x = S[0].x + s * D(i-1);
			S[i].fit_fun(ff, ud1, ud2);
		}
		solution PO, PE, PZ;
		matrix pc;
		int i_min, i_max;
		while (true)
		{
			i_min = i_max = 0;
			for (int i = 1; i < N; ++i)
			{
				if (S[i_min].y > S[i].y)
					i_min = i;
				if (S[i_max].y < S[i].y)
					i_max = i;
			}
			pc = matrix(n, 1);
			for (int i = 0; i < N; ++i)
				if (i != i_max)
					pc = pc + S[i].x;

			pc = pc / (N - 1.0);
			PO.x = pc + alpha * (pc - S[i_max].x);

			PO.fit_fun(ff, ud1, ud2);

			if (S[i_min].y <= PO.y && PO.y < S[i_max].y)
				S[i_max] = PO;
			else if (PO.y < S[i_min].y)
			{
				PE.x = pc + gamma * (PO.x - pc);
				PE.fit_fun(ff, ud1, ud2);

				if (PE.y < PO.y)
					S[i_max] = PE;
				else
					S[i_max] = PO;
			}
			else
			{
				PZ.x = pc + beta * (S[i_max].x - pc);
				PZ.fit_fun(ff, ud1, ud2);
				if (PZ.y < S[i_max].y)
					S[i_max] = PZ;
				else
				{
					for (int i = 0; i < N; ++i)
						if (i != i_min)
						{
							S[i].x = delta * (S[i].x + S[i_min].x);
							S[i].fit_fun(ff, ud1, ud2);
						}
				}
			}

			double max_s = norm(S[0].x - S[i_min].x);
			for (int i = 1; i < N; ++i)
				if (max_s < norm(S[i].x - S[i_min].x))
					max_s = norm(S[i].x - S[i_min].x);
			if (max_s < epsilon)
			{
				Xopt = S[i_min];
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt = S[i_min];
				Xopt.flag = 0;
				break;
			}
		}
		delete[]S;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

