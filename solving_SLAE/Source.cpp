#include <iostream>
#include <vector>
#include<time.h>
using namespace std;

void LU(vector <vector <double>> A, vector <vector <double>> &L,
	vector <vector <double>> &U, int n,vector<int> &transp) // transp - vector перестановок
{
	U = A;
	vector<double> tmp(n);
	transp.clear();
	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++) {
			if (U[i][i] != 0)L[j][i] = U[j][i] / U[i][i];//Если U[i][i]=0,то U[j]- в конец.
			else {
				transp.push_back(i);
				tmp = U[i];
				for (int idx = i; idx < n - 1; idx++)  //перестановка строки в конец
				{
					U[idx] = U[idx + 1];
				}
				U[n - 1] = tmp;
			}
		}
	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			for (int j = i; j < n; j++) {
				L[j][i] = U[j][i] / U[i][i];
			}
		for (int i = k; i < n; i++)
			for (int j = k - 1; j < n; j++)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}

}

void mult(vector <vector <double>> A, vector <vector <double>> B,
	vector <vector <double>> &R, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
}

void show(vector <vector <double>> A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << A[i][j] << "\t";
		}
		cout << endl;
	}
}
void slae_solution(vector <vector <double>> L,
	vector <vector <double>> U, vector<double>b, vector<double>&x, vector<int> transp, int n) {
	double sum = 0;
	vector <double> y(n);
	//согласование столбца b с матрицами L и U:
	double tmp;
	if (!transp.empty()) {
		for (int j = 0; j < transp.size(); j++) {
			tmp = b[transp[j]];
			for (int idx = transp[j]; idx < n - 1; idx++)  //перестановка строки в конец
			{
				b[idx] = b[idx + 1];
			}
			b[n - 1] = tmp;
		}
	}

	//1 step  Ly=b,  where Ux=y :
	y[0] = b[0];

	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			y[i] -= L[i][j] * y[j];
		}
		y[i] += b[i];
	}
	//2 step  Ux=y:
	for (int i = 0; i < n-1; i++)x[i] = 0;
	x[n - 1] = y[n - 1] / U[n - 1][n - 1];

	for (int i = n-2; i >= 0; i--)
	{
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= U[i][j] * x[j];
		}
		x[i] = (x[i] + y[i]) / U[i][i];
	}
}

void reverse_matrix(vector <vector <double>> L,
	vector <vector <double>> U, vector<double>b, vector < vector <double>>&A, vector<int> transp, int n) {
	
	vector<double> e(n);
	vector<double> x(n);
	vector<int> transportation;
	for (int i = 0; i < n; i++)
	{
		e[i] = 1;
		slae_solution(L, U, e, x, transportation, n);
		A[i] = x;
		e[i] = 0;
	}
	vector<double> tmp(n);
	int i;
	
	//перестановка в зависимости от вектора перестановок:
 	while (!transp.empty()) {
		i = transp.back();
		transp.pop_back();
		tmp = A[n - 1];
		for (int idx = n - 1; idx > i; idx--)  //перестановка строки в конец
		{
			A[idx] = A[idx - 1];
		}
		A[i] = tmp;
	}

	double temp_num;
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j <i; j++)
		{
			temp_num = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp_num;
		}
	}

}

double cond_number(vector <vector <double>> A,
	vector <vector <double>> B,int n) {
	double max_sum_A = -255;//!
	double curr_sum_A= 0;
	double max_sum_B = -255;//!
	double curr_sum_B = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			curr_sum_A += A[i][j];
			curr_sum_B += B[i][j];
		}
		if (curr_sum_A > max_sum_A)max_sum_A = curr_sum_A;
		if (curr_sum_B > max_sum_B)max_sum_B = curr_sum_B;
		curr_sum_A = 0;
		curr_sum_B = 0;
	}
	return max_sum_A*max_sum_B;
}
int main()
{
	const double n = 4;
	double det = 1.0;
	vector <vector <double>> A(n), L(n), U(n), R(n), L_(n), U_(n),A_1(n), TMP(n);
	vector <double> b(n);
	vector <double> y(n),x(n);
	srand(time(NULL));
	for (int i = 0; i < n; i++)
	{
		b[i] = rand() % 20;
		for (int j = 0; j < n; j++)
		{
			A[i].push_back(rand() % 20);			
			L[i].push_back(0);
			U[i].push_back(0);
			R[i].push_back(0);
		}
	}
	vector<int> transporation;
	cout << "Fisrt matrix" << endl;
	show(A, n);
	LU(A, L, U, n, transporation);
	for (int i = 0; i < n; i++)det *= U[i][i];
	cout << "U matrix" << endl;
	show(U, n);
	cout << "L matrix" << endl;
	show(L, n);
	mult(L, U, R, n);
	cout << "L*U matrix" << endl;
	show(R, n);
	cout << "detA= " <<det<<endl;
	cout << "b=  "<<endl;
	if (!transporation.empty())
		cout << "swap";
	for (int i = 0; i < n; i++) cout << b[i] << " " << endl;
	slae_solution(L, U, b, x, transporation, n);
	cout << "x=  " << endl;
	for (int i = 0; i < n; i++) cout << x[i] << " " << endl;

	reverse_matrix(L, U, b, A_1, transporation, n);
	
	cout << "A^-1=  " << endl;
	show(A_1, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			R[i][j] = 0;
		}
	}
	cout << "A * A^-1=  " << endl;
	mult(A, A_1, R,n);
	cout << '\n' << endl;
	show(R,n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			R[i][j] = 0;
		}
	}
	cout << "A^-1*A=  " << endl;
	mult(A_1, A, R, n);
	cout << '\n' << endl;
	show(R, n);

	cout << "cond_number=  " << endl;
	det=cond_number(A, A_1, n);
	cout << det << endl;
	system("pause");
	return 0;
}