using System;
using System.Collections.Generic;
using PolynomialFunction;
namespace Tinh_tich_phan
{
    public abstract class CalcIntegral{
        protected const double eta = 1e-7;
        protected double a;
        protected double b;
        protected int n;
        protected double step;
        public double PRINT_TO_SCREEN = 1e-15;
        protected double epsilon = 1e-15;
        protected double Mn;
        protected List<Point> inputXY;

        // protected Polynomial newton;
        protected Function fx;

        public double f(double x) => ThisFx.value(x);
        public abstract double Calc_Integral();
        public abstract int ConvertEpsilon_toN(double eps);
        public abstract double ConvertN_toEpsilon(int n);
        public List<Point> DivideAB_byN(int n){
            List<Point> result = new List<Point>{};
            double x_0 = ThisA;
            double x_n = ThisB;
            double x_i = 0;
            double y_i = 0;

            if(n>0){
                ThisStep = (ThisB-ThisA)/n;
                result.Add(new Point(x_0,f(x_0)));
                for(int i = 1; i < n; i++){
                    x_i = x_0 + ThisStep*i;
                    y_i = f(x_i);
                    Point newPoint = new Point(x_i,y_i);
                    result.Add(newPoint);
                }
                result.Add(new Point(x_n,f(x_n)));
                return result;
            }else{
                return null;
            }
        }
        /**
         *  ! -------------------------------------- GETTER_SETTER ----------------------------------------------------  
        **/
        public double ThisA{
            get{return this.a;}
            set{this.a = value;}
        }
        public double ThisB{
            get{return this.b;}
            set{this.b = value;}
        }
        public double ThisStep{
            get{return this.step;}
            set{this.step = value;}
        }
        public double ThisEpsilon{
            get{return this.epsilon;}
            set{this.epsilon = value;}
        }
        public double ThisMn{
            get{return this.Mn;}
            set{this.Mn = value;}
        }
        public int ThisN{
            get{return this.n;}
            set{this.n = value;}
        }
        public List<Point> ThisInputXY{
            get{return this.inputXY;}
            set{this.inputXY = value;}
        }
        public Function ThisFx{
            get{return this.fx;}
            set{this.fx = value;}
        }
        // public Polynomial ThisLagrange{
        //     get{return this.lagrange;}
        //     set{this.lagrange = value;}
        // }
    }
}