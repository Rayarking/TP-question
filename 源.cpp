
#include   <fstream>   
#include   <iomanip>   
#include   <iostream>   
#include   <cmath>   
using   namespace   std;   
#define   a(j)   (* (C+(M-1)*N+j))         //   销量数组
#define   b(i)   (* (C+i*N+N-1))                 //   产量数组   
#define   c(i,j)   (* (C+i*N+j))                 //   运价数组   
#define   x(i,j)   (* (X+i*(N-1)+j))         //   运量数组   
const   double   BIG_NUM = 1.0E15;                             //   任意大数   
															   //   (   <   BIG_NUM:   基本解,   >=   BIG_NUM   :   运量为   0   )   
#define   s(i,j)   (*(S+i*(N-1)+j))   //   检验数数组Sij   */   
#define   u(i)   (*(U+i))                                         //   位势数组   Ui   
#define   v(i)   (*(V+i))                                         //   位势数组   Vi   
#define   cpi(k)   ((CP+k)->i)                                 //   闭回路点   i   标   
#define   cpj(k)   ((CP+k)->j)                                 //   闭回路点   j   标   
#define   cpf(k)   ((CP+k)->f)                                 //   闭回路点   f   标   

															   /*
															   f   =   0:   j++;
															   f   =   1:   i--;
															   f   =   2:   j--;
															   f   =   3:   i++;
															   */

void   TP(int M, int N, double *C, double *X);

int   main()
{
	int   M, N, i, j;

	double*   C;   //   存储运价,   产量及销量   
	double*   X;   //   存储运量分配方案   
	double   z;

	ifstream   infile;
	char   fn[80];

	double   sum;
	cout.setf(ios_base::left, ios_base::adjustfield);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout.precision(3);

	cout << "请输入数据文件名:   ";
	cin >> fn;

	infile.open(fn);
	if (!infile)
	{
		cout << "文件打开失败!\n";
		system("pause");
		exit(0);
	}

	infile >> M >> N;

	M++;
	N++;

	X = new double[sizeof(double)*(M - 1)*(N - 1)];
	C = new double[sizeof(double)*M*N];

	//   把运价,   供应量和需求量的数据读入到数组   c(i,j)   
	for (i = 0; i<M; ++i)
		for (j = 0; j<N; ++j)
		{
			infile >> z;
			c(i, j) = z;
		}
	int s1 = 0, s2 = 0;


	infile.close();

	cout << "\n=============   数据文件   ================\n";
	for (i = 0; i<M; ++i)
	{
		for (j = 0; j<N; ++j)
			cout << setw(10) << c(i, j);

		cout << endl;
	}
	for (i = 0; i < N-1; ++i)
	{
		s1 += a(i);
	}
	for (j = 0; j < M-1; j++)
	{
		s2 += b(j);
	}
	if (s1< s2) {
		cout << "无解\n";
		return 0;
	}
	system("pause");

	TP(M, N, C, X);

	//   输出产销分配方案   
	cout << "\n=============   最优解   ===================\n";
	sum = 0;
	for (i = 0; i<M - 1; ++i)
	{
		for (j = 0; j<N - 1; ++j)
			if (x(i, j) >= BIG_NUM)
				cout << setw(10) << "******";
			else
			{
				cout << setw(10) << x(i, j);
				sum += (x(i, j)*c(i, j));
			}
		cout << endl;
	}

	//cout<<"\n\n\tThe min Cost is: %-10.4f\n", sum);   
	cout << "\n\n\t最高产量:" << setw(10) << sum << endl; //我们现在是在求max，max=-min   


	free(X);
	free(C);

	system("pause");
	return   0;
}

//   记录闭回路点结构   
struct   PATH
{
	int i, j, f;
};

void   TP(int M, int N, double* C, double* X)
{
	double *U, *V, *S;
	int  MN1, m, n;
	struct   PATH*   CP;
	int   k, i, j, l, k1, l1, ip;
	double  Cmin, sum;
	int I0, J0, Imin, Jmin;
	int fi, fj, fc, f;

	MN1 = (M - 1) + (N - 1) - 1;
	m = M - 1;
	n = N - 1;

	S = new double[sizeof(double)*(M - 1)*(N - 1)];
	U = new double[sizeof(double)*M];
	V = new double[sizeof(double)*N];
	CP = new PATH[sizeof(struct PATH)*(MN1 + 1)];

	//   解初始化   Xij   =   BIG_NUM   
	for (i = 0; i<m; ++i)
		for (j = 0; j<n; ++j)
			x(i, j) = BIG_NUM;

	//   最小元素法求初始可行解   
	for (k = 0; k < MN1; ++k)
	{
		Cmin = BIG_NUM;

		for (i = 0; i < m; ++i)
		{
			fi = 0;
			for (l = 0; l < k; ++l)
				//   去除已经用过的行   
				if (i == cpi(l))
				{
					fi = 1;
					break;
				}

			if (fi == 1)
				continue;
			for (j = 0; j < n; ++j)
			{
				fj = 0;
				for (l = 0; l < k; ++l)
					//   去除已经用过的列   
					if (j == cpj(l))
					{
						fj = 1;
						break;
					}

				if (fj == 1)
					continue;
				if (Cmin   >   c(i, j))
				{
					Cmin = c(i, j);
					I0 = i;
					J0 = j;
				}
			}         //   end   for   j   
		}                 //   end   for   i   

						  //   得到了未划去的最小运价所在格的坐标(I0,J0)和最小运价Cmin   
		if (k   >   0)
			if (Cmin == BIG_NUM && cpi(k - 1) == 0)
			{
				for (l1 = 0; l1<m; l1++)
					if (x(l1, cpj(k - 1)) == BIG_NUM)
						x(l1, cpj(k - 1)) = 0;
			}
			else if (Cmin == BIG_NUM && cpi(k - 1) != 0)

				for (l1 = 0; l1<n; l1++)
					if (x(cpi(k - 1), l1) == BIG_NUM)
						x(cpi(k - 1), l1) = 0;
		if (b(I0)   <   a(J0))
		{
			cpi(k) = I0;
			cpj(k) = -1;
			x(I0, J0) = b(I0);
			a(J0) -= b(I0);
			b(I0) = 0;
		}
		else
		{
			cpi(k) = -1;
			cpj(k) = J0;
			x(I0, J0) = a(J0);
			b(I0) -= a(J0);
			a(J0) = 0;
		}
	}   //   end   for   k   用最小元素法求得了初使可行解   
		//   输出初始可行解   
	cout << "\n=============   初始解   ===================\n";
	sum = 0;

	for (i = 0; i < M - 1; i++)
	{
		for (j = 0; j < N - 1; j++)
			if (x(i, j) >= BIG_NUM)
				cout << setw(10) << "******";
			else
			{
				cout << setw(10) << x(i, j);
				sum += (x(i, j)   *   c(i, j));
			}
		cout << endl;
	}

	cout << "\n\n\t初始产量:" << setw(10) << sum << endl;//我们现在是在求max，max=-min   
	system("pause");
	while (true)
	{
		//   位势置初值   Ui,   Vi   =   BIG_NUM   
		for (i = 0; i < m; i++)
			u(i) = BIG_NUM;
		for (j = 0; j < n; j++)
			v(j) = BIG_NUM;
		//   求位势   
		l = 0;
		u(0) = 0;
		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++)
			{
				if (u(i) >= BIG_NUM && v(j) >= BIG_NUM && x(i, j)<BIG_NUM)
				{         //   记录未求过位势的位置   
					cpi(l) = i;
					cpj(l) = j;
					l++;
				}
				else if (x(i, j)<BIG_NUM && u(i)<BIG_NUM)
					v(j) = c(i, j) - u(i);
				else if (x(i, j)<BIG_NUM && v(j)<BIG_NUM)
					u(i) = c(i, j) - v(j);
			}
		//   按记录位置求其余位势   
		if (l>0)
			while (true)
			{
				ip = 0;
				for (k = 0; k<l; k++)
				{
					i = cpi(k);
					j = cpj(k);
					if (u(i) >= BIG_NUM && v(j) >= BIG_NUM && x(i, j)<BIG_NUM)
					{//记录未求过位势的位置   
						cpi(ip) = i;
						cpj(ip) = j;
						ip++;
					}
					else if (x(i, j)<BIG_NUM && u(i)<BIG_NUM)
						v(j) = c(i, j) - u(i);
					else if (x(i, j)<BIG_NUM && v(j)<BIG_NUM)
						u(i) = c(i, j) - v(j);
				}//end for k   
				if (ip == 0)
					break;
				l = ip;
			}   //   end   for   while   
				//   求检验数   
		for (i = 0; i < m; i++)
			for (j = 0; j<n; j++)
			{
				s(i, j) = BIG_NUM;
				if (x(i, j) >= BIG_NUM)
					s(i, j) = c(i, j) - u(i) - v(j);
			}
		//   求最小检验数   
		Cmin = BIG_NUM;
		for (i = 0; i<m; i++)
			for (j = 0; j<n; j++)
				if (Cmin>s(i, j))
				{
					Cmin = s(i, j);
					I0 = i;
					J0 = j;
				}
		if (Cmin >= 0)
			return;   //   已经得到最优解,返回主程序   
					  //   此时找到了入基变量   X(   I0,   J0   )   
					  //   求闭回路   
		for (k = 0; k < MN1; k++)
			cpf(k) = -1;
		cpi(0) = I0;
		cpj(0) = J0;
		k = 0;
		while (true)
		{
			f = cpf(k);

			//   设置闭回路搜索方向   
			while (true)
			{
				i = cpi(k);
				j = cpj(k);
				fc = 0;
				f++;
				cpf(k) = f;
				if (f >= 4)
					break;
				//   避免反向搜索   
				if (k>0)
				{
					if (f == 0 && cpf(k - 1) == 2)
						continue;
					else   if (f == 1 && cpf(k - 1) == 3)
						continue;
					else   if (f == 2 && cpf(k - 1) == 0)
						continue;
					else   if (f == 3 && cpf(k - 1) == 1)
						continue;
				}
				if (f == 0)
				{  //   沿   j+   方向搜索   
					while (true)
					{
						j++;
						if (j >= n)
						{
							fc = 2;   break;
						}
						if (i == I0 && j == J0)
						{
							fc = 1;   break;
						}
						if (s(i, j) >= BIG_NUM)
						{
							fc = 3;   break;
						}
					}
				}   //   end   for   j+   
				else   if (f == 1)
				{         //   沿   i-   方向搜索   
					while (true)
					{
						i--;
						if (i   <   0)
						{
							fc = 2;   break;
						}
						if (i == I0 && j == J0)
						{
							fc = 1;   break;
						}
						if (s(i, j) >= BIG_NUM)
						{
							fc = 3;   break;
						}
					}
				}         //   end   for   i-   
				else   if (f == 2)
				{         //   沿   j-   方向搜索   
					while (true)
					{
						j--;
						if (j   <   0)
						{
							fc = 2;   break;
						}
						if (i == I0 && j == J0)
						{
							fc = 1;   break;
						}
						if (s(i, j) >= BIG_NUM)
						{
							fc = 3;   break;
						}
					}
				}   //   end   for   j-   
				else   if (f == 3)
				{         //   沿   i+   方向搜索   
					while (true)
					{
						i++;
						if (i >= m)
						{
							fc = 2;   break;
						}
						if (i == I0 && j == J0)
						{
							fc = 1;   break;
						}
						if (s(i, j) >= BIG_NUM)
						{
							fc = 3;   break;
						}
					}
				}   //   end   for   i+   
				if (fc == 1 || fc == 3)
					break;
			}         //   end   for   while   flag   2   
			if (fc == 0)
				k--;         //   沿此方向搜索失败,退回到前一点   
			else   if (fc == 1)
				break;   //   搜索完成   
			else   if (fc == 3)
			{         //   沿此方向搜索成功,前进一点   
				k++;
				cpi(k) = i;
				cpj(k) = j;
				cpf(k) = -1;
			}
		}         //   end   while   
				  //   去除闭回路中的非转折点   
		l = 0;
		while (l   <   k - 1)
		{
			i = l + 1;
			while (i <= k)
			{
				if (cpf(l) == cpf(i))
					i++;
				else
					break;
			}
			if (i   >(l + 1))
			{
				j = l + 1;
				k1 = k - (i - j);
				//   如果某些点前进方向相同,则去除中间各点 
				while (i <= k)
				{
					cpi(j) = cpi(i);
					cpj(j) = cpj(i);
					cpf(j) = cpf(i);
					i++;
					j++;
				}
				l += 2;
				k = k1;
			}
			else
				l++;
		}         //   end   for   while   l   <   k   -   1   
				  //   查找闭回路上基本解的最小值   
		Cmin = x(cpi(1), cpj(1));
		Imin = cpi(1);
		Jmin = cpj(1);
		for (i = 3; i <= k; i += 2)
			if (Cmin   >   x(cpi(i), cpj(i)))
			{
				Cmin = x(cpi(i), cpj(i));
				Imin = cpi(i);
				Jmin = cpj(i);
			}
		//   换入基变量   
		x(I0, J0) = Cmin;
		for (i = 1; i <= k; i += 2)
		{
			x(cpi(i), cpj(i)) -= Cmin;
			cout << "hao" << x(cpi(i), cpj(i));
			if ((i + 1) <= k)
				x(cpi(i + 1), cpj(i + 1)) += Cmin;
		}
		x(Imin, Jmin) = BIG_NUM;
	}
	delete[]   CP;
	delete[]   V;
	delete[]   U;
	delete[]   S;
 }