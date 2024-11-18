// AlgoPointTask.cpp : ;Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include<vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iomanip>


using namespace std;
int sing(double a) {
	if (a > 0) return 1;
	else if (a < 0) return -1;
	else return 0;
}
vector <vector<float>>Transpos(vector<vector<float>>& A) {
	vector <vector<float>>At(A[0].size(), vector<float>(A.size()));
	for (int i = 0; i < At.size(); i++)
	{
		for (int j = 0; j < At[0].size(); j++)
		{
			At[i][j] = A[j][i];
		}
	}
	return At;
}
vector < float > MxV(vector<vector<float>>& A, vector <float>& x) {
	vector<float> result(A.size());

	int temp = 0;

	for (int i = 0; i < A.size(); i++) {
		temp = 0;
		for (int j = 0; j < A[0].size(); j++)
		{
			temp += x[j] * A[i][j];
		}
		result[i] = temp;
	}
	return result;


}
vector<vector<float>>MxM(vector<vector<float>>& A, vector<vector<float>>& B)
{
	vector<vector<float>> T(A.size(), vector<float>(B[0].size()));
	for (int i = 0; i < A.size(); i++) {

		for (int j = 0; j < B[0].size(); j++)
		{
			T[i][j] = 0;
			for (int k = 0; k < B.size(); k++) T[i][j] += A[i][k] * B[k][j];

		}
	}
	return T;
}
float norma(vector<float>& r) {
	float result = 0;
	for (int i = 0; i < r.size(); i++) {
		result += r[i] * r[i];
	}
	return sqrt(result);
}
float Skalar(vector<float>& a, vector<float>& b)
{
	float res = 0;
	for (int i = 0; i < a.size(); i++) {
		res += a[i] * b[i];
	}
	return res;
}
vector<float> R(vector<vector<float>>& A, vector<float>& x, vector<float>& b, vector<float>& r) {
	vector<float>Ax = MxV(A, x);

	for (int i = 0; i < A.size(); i++)
		r[i] = b[i] - Ax[i];

	return r;
}
vector<vector<float>> Diag(vector<vector<float>>& D, vector<float>& x)
{

	for (int j = 0; j < x.size(); j++) {
		D[j][j] = x[j] * x[j];
	}
	return D;
}
vector<float > SolveSlau(vector<vector<float>>& T, vector<float>& q) {
	vector<float> d(T[0].size());
	vector<float> y(T[0].size());
	vector< vector<float>>S(T[0].size(), vector<float>(T[0].size()));

	vector<float> u(T[0].size());
	for (int k = 0; k < T.size(); k++) {
		double temp = 0.0;
		for (int i = 0; i < k - 1; i++)temp += d[i] * pow(abs(S[i][k]), 2);
		d[k] = sing(T[k][k] - temp);
		S[k][k] = sqrt(abs(T[k][k] - temp));
		double temp2 = 0.0;
		for (int j = k; j < T.size(); j++) {
			for (int i = 0; i < k - 1; i++)temp2 += d[i] * S[i][k] * S[i][j];
			S[k][j] = (T[k][j] - temp2) / (S[k][k] * d[k]);

		}
	}
	y[0] = q[0] / (T[0][0] * d[0]);
	for (int i = 0; i < T.size(); i++) {
		double temp = 0.0;
		for (int k = 0; k < i - 1; k++)temp += d[k] * y[k] * S[k][i];
		y[i] = (q[i] - temp) / (S[i][i] * d[i]);

	}
	u[T[0].size() - 1] = y[T[0].size() - 1] / S[T[0].size() - 1][T[0].size() - 1];
	for (int i = T[0].size() - 1; i >= 0; i--) {
		double temp = 0.0;
		for (int k = i + 1; k < T[0].size(); k++) temp += u[k] * S[i][k];
		u[i] = (y[i] - temp) / S[i][i];

	}
	return u;
}
vector<float> U(vector<vector<float>>& A, vector<vector<float>>& D, vector<float>& c, vector<float>& r) {
	vector<vector <float>> At;
	vector<vector<float>> T(A.size(), vector<float>(A[0].size()));
	vector<vector <float>> AD;
	vector<float>u;
	vector<float>ADc;
	vector<float>q(A[0].size());
	AD = MxM(A, D);
	ADc = MxV(AD, c);
	for (int i = 0; i < A.size(); i++) {
		q[i] = r[i] + ADc[i];
	}
	At = Transpos(A);
	T = MxM(AD, At);


	u = SolveSlau(T, q);
	return u;
}
vector<float> G(vector<vector<float>>& D, vector<vector<float>>& A, vector<float>& u, vector<float>& c) {
	vector<vector<float>>At = Transpos(A);
	vector<float>ATu = MxV(At, u);
	vector<float> g(A[0].size());
	vector<float> s(u.size());
	for (int i = 0; i < 6; i++) {
		g[i] = c[i] - ATu[i];
	}
	return g;
}
vector<float>S(vector<vector<float>>& D, vector<vector<float>>& A, vector<float>& u, vector<float>& c) {


	vector<float> g(A[0].size());
	vector<float> s(u.size());

	g = G(D, A, u, c);
	s = MxV(D, g);
	for (int i = 0; i < g.size(); i++) {
		s[i] = -s[i];
	}

	return s;
}
float MinBy2Num(float a, float b) {
	if (a < b)return a;
	else return b;
}
float findLamda(vector<float>& s, vector<float>& x, vector<float>& r) {

	float l_ = 1; //lamda с чертой
	float L = 0; //laMDA
	float min; //lamda штрих
	float gamma = 1 / 3.0;
	vector<float>MinSet;

	for (int i = 0; i < x.size(); i++) {

		if (s[i] >= 0) continue;
		MinSet.push_back(-x[i] / s[i]);


	}


	for (int i = 0; i < r.size(); i++) {
		if (!MinSet.empty()) { //у казызвает на условие что  не все !( s[i]>0) 
			min = *min_element(MinSet.begin(), MinSet.end());
			l_ = min * gamma;
			MinSet.erase(MinSet.begin(), MinSet.end());
		}
		else if (r[i] > 1e-6 && !MinSet.empty()) {
			L = MinBy2Num(1, l_);
		}
		else if (r[i] > 1e-6 && MinSet.empty())
		{
			L = 1;

		}
		else if (r[i] < 1e-6)
		{
			L = l_;

		}







	}
	return L;
}

void Solution(vector<vector<float>>& A, vector<float>& c, vector<float>& b, vector<float>& x) {
	float eps = 1e-6;
	int n = 9;
	int m = 6;


	vector<float> r(m);
	vector<float> E(n);
	vector<float> s(n);

	vector<vector<float>>D(n, vector<float>(n));


	vector<float> u;
	vector<float> g;

	float Lamda;
	r = R(A, x, b, r);
	D = Diag(D, x);

	u = U(A, D, c, r);

	s = S(D, A, u, c);

	while (norma(r) > 1e-6)
	{

		r = R(A, x, b, r);
		D = Diag(D, x);

		u = U(A, D, c, r);

		s = S(D, A, u, c);
		Lamda = findLamda(s, x, r);
		for (int i = 0; i < x.size(); i++) {
			x[i] = x[i] + Lamda * s[i];

		}
		for (int i = 0; i < r.size(); i++) {

			r[i] = (1 - Lamda) * r[i];
		}

	}



	std::cout << "[";
	for (int i = 0; i < x.size(); i++) {
		std::cout << x[i] << " ";

	}
	std::cout << "]";
	std::cout << endl;

	std::cout << "cTx = " << setprecision(9) << Skalar(c, x) << endl;
}
int main()
{


	int n = 9;
	int m = 6;
	vector<float> x(9, 1);
	vector < vector < float >> A = { {1, 1, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 1, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 1, 1},
		{1, 0, 0, 1, 0, 0, 1, 0, 0},
		{0, 1, 0, 0, 1, 0, 0, 1, 0},
		{0, 0, 1, 0, 0, 1, 0, 0, 1} };


	vector<float> c = { -5, 0, -3, -3, -7, -9, -3, -5, -2 };
	vector<float> b = { 1, 1, 1, 1, 1, 1 };

	Solution(A, c, b, x);
}

