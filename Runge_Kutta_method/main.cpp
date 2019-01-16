#include <iostream>
#include <cmath>
#include<fstream>
using namespace std;

double c1, c2, c3, t0, t1, t2;
//c1: Ex0/m, c2:Ey0/m, c3: qB0/m

double func1(double t, double x, double v, double y, double w)
{
	return v;
}

double func2(double t, double x, double v, double y, double w)
{
	return c1 * exp(-pow((t / t1), 2)) + c3 * w*tanh(t / t0);
}

double func3(double t, double x, double v, double y, double w)
{
	return w;
}

double func4(double t, double x, double v, double y, double w)
{
	return c2 * exp(-pow((t / t2), 2)) - c3 * v*tanh(t / t0);
}

int main()
{
	double k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4;
	double h;
	int n = 200;
	double *x, *v, *y, *w, *t;

	t = new double[n + 1];
	x = new double[n + 1];
	v = new double[n + 1];
	y = new double[n + 1];
	w = new double[n + 1];

	//常参数赋值
	c1 = 7.9, c2 = 7.7, c3 = 0.95, t0 = 1.0, t1 = 1.0, t2 = 1.0;

	//初值
	t[0] = 0.0;
	t[n] = 10.0;
	x[0] = 0;
	v[0] = 0;
	y[0] = 0;
	w[0] = 0;
	h = (t[n] - t[0]) / n;

	for (int i = 1; i <= n; i++)
	{
		t[i - 1] = t[0] + (i - 1)*h;

		k1 = func1(t[i - 1], x[i - 1], v[i - 1], y[i - 1], w[i - 1]);
		l1 = func2(t[i - 1], x[i - 1], v[i - 1], y[i - 1], w[i - 1]);
		m1 = func3(t[i - 1], x[i - 1], v[i - 1], y[i - 1], w[i - 1]);
		n1 = func4(t[i - 1], x[i - 1], v[i - 1], y[i - 1], w[i - 1]);

		k2 = func1(t[i - 1] + h / 2.0, x[i - 1] + h * k1 / 2.0, v[i - 1] + h * l1 / 2.0, y[i - 1] + h * m1 / 2.0, w[i - 1] + h * n1 / 2.0);
		l2 = func2(t[i - 1] + h / 2.0, x[i - 1] + h * k1 / 2.0, v[i - 1] + h * l1 / 2.0, y[i - 1] + h * m1 / 2.0, w[i - 1] + h * n1 / 2.0);
		m2 = func3(t[i - 1] + h / 2.0, x[i - 1] + h * k1 / 2.0, v[i - 1] + h * l1 / 2.0, y[i - 1] + h * m1 / 2.0, w[i - 1] + h * n1 / 2.0);
		n2 = func4(t[i - 1] + h / 2.0, x[i - 1] + h * k1 / 2.0, v[i - 1] + h * l1 / 2.0, y[i - 1] + h * m1 / 2.0, w[i - 1] + h * n1 / 2.0);

		k3 = func1(t[i - 1] + h / 2.0, x[i - 1] + h * k2 / 2.0, v[i - 1] + h * l2 / 2.0, y[i - 1] + h * m2 / 2.0, w[i - 1] + h * n2 / 2.0);
		l3 = func2(t[i - 1] + h / 2.0, x[i - 1] + h * k2 / 2.0, v[i - 1] + h * l2 / 2.0, y[i - 1] + h * m2 / 2.0, w[i - 1] + h * n2 / 2.0);
		m3 = func3(t[i - 1] + h / 2.0, x[i - 1] + h * k2 / 2.0, v[i - 1] + h * l2 / 2.0, y[i - 1] + h * m2 / 2.0, w[i - 1] + h * n2 / 2.0);
		n3 = func4(t[i - 1] + h / 2.0, x[i - 1] + h * k2 / 2.0, v[i - 1] + h * l2 / 2.0, y[i - 1] + h * m2 / 2.0, w[i - 1] + h * n2 / 2.0);

		k4 = func1(t[i - 1] + h, x[i - 1] + h * k3, v[i - 1] + h * l3, y[i - 1] + h * m3, w[i - 1] + h * n3);
		l4 = func2(t[i - 1] + h, x[i - 1] + h * k3, v[i - 1] + h * l3, y[i - 1] + h * m3, w[i - 1] + h * n3);
		m4 = func3(t[i - 1] + h, x[i - 1] + h * k3, v[i - 1] + h * l3, y[i - 1] + h * m3, w[i - 1] + h * n3);
		n4 = func4(t[i - 1] + h, x[i - 1] + h * k3, v[i - 1] + h * l3, y[i - 1] + h * m3, w[i - 1] + h * n3);

		x[i] = x[i - 1] + h / 6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		v[i] = v[i - 1] + h / 6.0*(l1 + 2.0*l2 + 2.0*l3 + l4);
		y[i] = y[i - 1] + h / 6.0*(m1 + 2.0*m2 + 2.0*m3 + m4);
		w[i] = w[i - 1] + h / 6.0*(n1 + 2.0*n2 + 2.0*n3 + n4);
	}

	ofstream outfile("data.csv", ios::ate);
	for (int i = 0; i <= n; i++)
		outfile << t[i] << ',' << x[i] << ',' << v[i] << ',' << y[i] << ',' << w[i] << endl;
	outfile.close();
}