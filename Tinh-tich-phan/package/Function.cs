using System;
using System.Collections.Generic;
using PolynomialFunction;
namespace Tinh_tich_phan
{
    public class Function{
        private const double eta = 1e-7;

        // private static Polynomial poly = null;
        private List<List<double>> extremes_of_fn = new List<List<double>>{
            new List<double>{}, // ! Extremes of f0
            new List<double>{}, // ! Extremes of f1
            new List<double>{}, // ! Extremes of f2
            new List<double>{}, // ! Extremes of f3
            new List<double>{}, // ! Extremes of f4
        };
        private static Func<double,double> f_value = f; 
        private static Func<double,double> f1_value = f1; 
        private Random rand = new Random();
        public static double f(double x){
            return 1.0/(x*x+1);
        }
        public double value(double x) => f(x);

        public Function(){

        }
        public void SwitchFunc(int n){
            switch (n)
            {
                case 0:
                    f_value = f;
                    f1_value = f1;
                break;
                case 1:
                    f_value = f1;
                    f1_value = f2;
                break;
                case 2:
                    f_value = f2;
                    f1_value = f3;
                break;
                case 3:
                    f_value = f3;
                    f1_value = f4;
                break;
                case 4:
                    f_value = f4;
                    f1_value = f5;
                break;
                default:
                    Console.WriteLine("ERROR!!");
                break;
            }
        }
        public static double f1(double x) => (f(x + eta) - f(x - eta))/(2*eta);
        public static double f2(double x) => (f1(x + eta) - f1(x - eta))/(2*eta);
        public static double f3(double x) => (f2(x + eta) - f2(x - eta))/(2*eta);
        public static double f4(double x) => (f3(x + eta) - f3(x - eta))/(2*eta);
        public static double f5(double x) => (f4(x + eta) - f4(x - eta))/(2*eta);

        public bool Sign_of(double value) => value>0;
        public double Random_num(double a, double b) => a + rand.NextDouble()*(b-a);
        public double FindMax_Fn(int n, double a, double b){       // !!!! Not short + complete
            
            SwitchFunc(n);
            int countExtremes = extremes_of_fn[n].Count; 
            // Console.WriteLine(countExtremes);
            double max = 0;
            double x = 0;
            double fa=0;
            double fb =0;
            if(countExtremes==0){
                extremes_of_fn[n] = ListExtremes_x_byGD(n,a,b);
                countExtremes = extremes_of_fn[n].Count; 
            }
            /* Sau khi thu được cực trị */
            if(countExtremes!=0){

                x = extremes_of_fn[n][0];
                max = Math.Abs(f_value(x));
                // max = fn(n,x);
                for(int i = 1; i < countExtremes; i++){
                    double temp = Math.Abs(f_value(extremes_of_fn[n][i]));
                    // Console.WriteLine(temp); 
                    // double temp = fn(n,extremes_of_fn[n][i]); 
                    if(max < temp){
                        max = temp;
                    }
                }
                return max;
            }else{
                // ! Nếu là hàm hằng hoặc là không có
                fa = Math.Abs(f_value(a));
                fb = Math.Abs(f_value(b));
                return (fa>fb)? fa:fb;
            }
        }

        List<double> ListExtremes_x_byGD(int n, double a, double b){
            List<double> extremes = new List<double>{};
            int MAX_LOOP = 10000;
            int loop=0;
            int counting = 0;
            /* for Gradient  */
            int sign = 0;
            double alp = 0;
            double x = 0;
            double f1x = 0;

            /* Kiểm tra hàm hằng */
            for(int i = 1; i<=100; i++){
                if(Math.Abs(f3(Random_num(a,b)))<1e-10){
                    counting++;
                }else{
                    break;
                }
                if(i==100) return null;
            }

            /* Gradient Descend */
            x = a;
            f1x = f1_value(x);
            // f1x = fn(n+1,x);
            sign = (f1x<0)? -1:1;
            alp = -((int)Math.Log10(Math.Abs(f1x))+1);
            alp = Math.Exp(Math.Log(10)*alp);
            while(x <= b && ++loop<MAX_LOOP){
                if(Math.Abs(f1x)<1e-10 || f1x==0){
                    extremes.Add(x);
                    x += alp;
                    f1x = f1_value(x);
                    // f1x = fn(n+1,x);
                    sign = (f1x<0)? -1:1;
                    continue;
                }
                alp = Change_alpha(n,alp,sign,x,f1x);
                x += sign * alp;
                f1x = f1_value(x);
                // f1x = fn(n+1,x);
                if(f1x * sign<0) sign*=(-1);
            }
            return extremes;
        }

        public double Change_alpha(int n, double alpha, double sign, double x, double f1x){
            double new_x = x + sign * alpha * f1x;
            double f1_new = f1_value(new_x);
            // double f1_new = fn(n+1,new_x);
            if(f1x * f1_new > 0 && alpha<1){//
                alpha*=1.5;
            }else if(f1x*f1_new < 0){
                alpha/=1.5;
            }
            return alpha;
        }
        // public Polynomial ThisPoly{
        //     get{return this.poly;}
        //     set{this.poly = value;}
        // }
        
    }
}