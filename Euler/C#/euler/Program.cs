using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace euler
{
    class Program
    {
        static double h, x0, y0, X;
        static int size, n, vidu;
        static double[] CX;//Ngiem chinh xac
        static double[] A;
        static double[] B;
        static double[] C;
        static double[] Y; //Ma tran cac gia tri ban dau

        static void VD2()
        {
            vidu = 2;
            //h = 1;
            // h = 0.1;
            // h = 0.01;
            h = 0.001;
            X = 1;
            x0 = 0;
            y0 = 1;
            size = (int)((X - x0) / h + 1);
            A = new double[size];
            B = new double[size];
            C = new double[size];

            //Tính nghiệm chính xác   
            CX = new double[size];
            int k = 0;
            for (double i = 0; i <= X; i = i + h)
            {
                CX[k] = 2*Math.Pow(Math.E,i)-i-1;
                k++;
            }

            forward();
            backward();
            hinhthang();
        }
        static double f(double x, double y)//VD1 - VD2: giải phương trình vi phân cấp 1
        {
            // if (vidu == 1)
            return x + y; //e^(x) = y;
            // return y - 0.5 * Math.Pow(Math.E, 0.5 * x) * Math.Sin(5 * x) + 5 * Math.Pow(Math.E, 0.5 * x) * Math.Cos(5 * x);
        }
        static void VD3()
        {
            vidu = 3;
            h = 0.1;
            X = 1;
            x0 = 0;
            n = 2;
            size = (int)((X - x0) / h + 1);
            A = new double[size];
            Y = new double[n];
            //dieu kien ban dau: y(0) = 2; y'(0) = 1
            Y[0] = 2;
            Y[1] = 1;

            //Tính nghiệm chính xác
            CX = new double[size];
            int k = 0;
            for (double i = 0; i <= X; i = i + h)
            {
                CX[k] = Math.Pow(Math.E, -i) * (2 * Math.Cos(i) + 3 * Math.Sin(i));
                k++;
            }
            N_Euler_hinhthang();
            N_Euler_hien();
        }

        static void VD4()
        {
            vidu = 4;
            h = 0.1;
            X = 1;
            x0 = 0;
            n = 3;
            size = (int)((X - x0) / h + 1);
            A = new double[size];
            Y = new double[n];
            //dieu kien ban dau: y(0) = 1; y'(0) = 0          
            Y[0] = 1;
            Y[1] = 1;
            Y[2] = 1;

            //Tính nghiệm chính xác
            CX = new double[size];
            int k = 0;
            for (double i = 0; i <= X; i = i + h)
            {
                CX[k] = Math.Pow(Math.E, i);
                k++;
            }
            N_Euler_hinhthang();
            N_Euler_hien();
        }


        static double ptvp(double x, double[] f)//VD3 - VD4: giải phương trình vi phân cấp k
        {
            if (vidu == 3)
                return -2 * f[0] - 2 * f[1];
            return f[0] - f[1] + f[2];
        }

        static void forward()
        {
            double y = y0, x = x0;
            for (int i = 0; i < size; i++)
            {
                A[i] = y;
                y = y + h * f(x, y);
                x = x + h;
            }
        }
        static void backward()
        {
            double y = y0, x = x0, I;
            for (int i = 0; i < size; i++)
            {
                B[i] = y;
                I = h * f(x + h, y + h * f(x, y));
                x = x + h;
                y = y + I;
            }
        }
        static void hinhthang()
        {
            double y = y0, x = x0, I;
            for (int i = 0; i < size; i++)
            {
                C[i] = y;
                I = (h / 2) * (f(x, y) + f(x + h, y + h * f(x, y)));
                y = y + I;
                x = x + h;
            }
        }

        static void N_Euler_hien() // GIải ptvp cấp k bằng pp Euler hiện
        {
            double x = x0;
            double[] Z = Y;
            for (int i = 0; i < size; i++)
            {
                A[i] = Z[0];
                for (int j = 0; j < n - 1; j++)
                {
                    Z[j] = Z[j] + h * Z[j + 1];
                }
                Z[n - 1] = Z[n - 1] + h * ptvp(x, Z);
                x = x + h;
            }
        }

        static void N_Euler_hinhthang() // Giải ptvp cấp k bằng pp hình thang
        {
            B = new double[size];
            double x = x0;
            double[] Z = Y;
            double[] Z1 = new double[n];
            double[] Z2 = new double[n];
            for (int i = 0; i < size; i++)
            {
                B[i] = Z[0];
                for (int j = 0; j <= n - 2; j++)
                {
                    Z2[j] = Z[j] + h * Z[j + 1];
                }
                Z2[n - 1] = Z[n - 1] + h * ptvp(x, Z);
                for (int j = 0; j <= n - 2; j++)
                {
                    Z1[j] = Z[j] + 0.5 * h * (Z[j + 1] + Z2[j + 1]);
                }
                Z1[n - 1] = Z[n - 1] + 0.5 * h * (ptvp(x, Z) + ptvp(x + h, Z2));
                x = x + h;
                Z = Z1;
            }
        }
        static void Main(string[] args)
        {
            //VD1();
            VD2();
            Console.WriteLine("x           Chính xác     \t   Euler hien   \t   Euler an    \t   Hinh thang");
            for (int i = 0; i < CX.Length; i++)
            {
                Console.WriteLine("{0:F5}     {1:f5}     \t    {2:F5}      \t   {3:f5}   \t   {4:f5}", x0 + i * h, CX[i], A[i], B[i], C[i]);
            }

            //VD3();
            //VD4();
            //Console.WriteLine("x           \t  Chính xác  \t   Euler hien  \t  Hinh thang");
            //for (int i = 0; i < CX.Length; i++)
            //{
            //   Console.WriteLine("x: {0:F5}   \t   {1:F5}   \t   {2:F5}   \t   {3:F5}", x0 + i * h, CX[i], A[i], B[i]);
            // }
            Console.Read();
        }
    }
}
