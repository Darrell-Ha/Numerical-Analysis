using System.Collections.Generic;
using System;
namespace PolynomialFunction
{
    public class Point{
        private double x;
        private double y;

        public Point(){
            
        }
        public Point(double x, double y){
            this.x = x;
            this.y = y;
        }

        public Point(double x){ //// by function
            this.x = x;
            this.y = Math.Exp(-x)+2/(3*x*x+1);
        }

        public double ThisX{
            get {return this.x;}
            set {this.x = value;}
        }

        public double ThisY{
            get {return this.y;}
            set {this.y = value;}
        }
    }
    public class Polynomial{
        private List<double> coeffs;

        public Polynomial(List<double> coeffs){
            this.coeffs = coeffs;
        }
        public Polynomial(){
            this.coeffs = new List<double>{0};
        }
        public override string ToString(){
            string result="";
            int deg = coeffs.Count-1;
            double temp=0;
            
            for (int i=0; i<=deg; i++){
                temp = coeffs[i];
                if(i==0 && temp<0){}

                if(temp==0) continue;
                if(temp<0||i==0){
                    result+=$"{temp}";
                }else{
                    result+=$"+ {temp}";
                }
                if(deg-i!=0){
                    result+=$"   x^{deg-i}\n";
                }            
            }
            return result;
        }

        public List<double> ThisCoeffs{
            get{return this.coeffs;}
            set{this.coeffs=value;}
        }

        public static Polynomial operator +(Polynomial poly) => poly;
        public static Polynomial operator +(Polynomial poly, double num){
            Polynomial result = poly;
            int index_Deg_0 = result.ThisCoeffs.Count-1;
            result.ThisCoeffs[index_Deg_0] += num;
            return result;
        }

        public static Polynomial operator +(Polynomial poly1, Polynomial poly2){

            List<double> coef1 = poly1.ThisCoeffs;
            List<double> coef2 = poly2.ThisCoeffs;
            List<double> result = new List<double>{};
            int countAdd = 0;
            int deg1= coef1.Count-1;
            int deg2= coef2.Count-1;
            int degMax = (deg1>deg2)? deg1: deg2;      /* not completed */
            int degMin = (deg1<deg2)? deg1: deg2;
            countAdd = degMax-degMin;
            for(int k = 0; k<countAdd; k++){
                if(deg1>deg2){
                    coef2.Insert(0,0);
                }else{
                    coef1.Insert(0,0);
                }
            }
            for(int i=0; i<=degMax; i++){
                double plus = coef1[i]+coef2[i]; 
                result.Add(plus);
            }
            Polynomial resultPoly = new Polynomial(result);
            return resultPoly;
        }
        public static Polynomial operator *(double num, Polynomial poly){
            List<double> coef = new List<double>(poly.ThisCoeffs);
            int deg = coef.Count-1;
            for(int i = 0; i <= deg; i++){
                coef[i] *= num;
            }
            Polynomial result = new Polynomial(coef);
            return result;
        }
        public static Polynomial operator /(Polynomial poly, double num){
            if(num == 0){
                throw new DivideByZeroException();
            }else{
                return (1/num)*poly;
            } 
        }
        public static Polynomial operator -(Polynomial poly) => (-1)*poly;
        public static Polynomial operator -(Polynomial poly1, Polynomial poly2) => poly1 + (-poly2);

        public static Polynomial operator *(Polynomial poly1, Polynomial poly2){
            Polynomial temp1 = poly1;
            Polynomial temp2 = poly2;
            Polynomial result = new Polynomial();
            // int deg1 = temp1.ThisCoeffs.Count-1;
            int deg2 = temp2.ThisCoeffs.Count-1;

            for(int i = 0; i <= deg2; i++){
                Polynomial temp = temp2.ThisCoeffs[i]*temp1;
                for(int times = 1; times <=(deg2-i); times++){
                    temp.ThisCoeffs.Add(0);
                }
                result+=temp;
            }
            return result;
            
        }
        /** 
         * ? Horner: chia (division)
        ** function solving p(x)=(x-c)q(x)+p(c)
        *!  input:
        **      + coeficient of q(x)
        **      + tuple (x-c) = (1,-c) to calculate p(c)
        *!  return tuple (q(x), p(c))
        **/
        public static (Polynomial Qx,double Pc) operator/(Polynomial poly, (int x, double _c) x_c){
            
            List<double> tempPoly = new List<double>{};
            tempPoly = poly.ThisCoeffs;
            int deg = tempPoly.Count-1;
            double pc = 0;
            List<double> newCoeff =  new List<double>{};

            newCoeff.Add(tempPoly[0]);
            for(int i = 1; i<=deg; i++){
                newCoeff.Add(0);
                newCoeff[i]=tempPoly[i]+newCoeff[i-1]*(-x_c._c);
            }
            pc = newCoeff[deg];
            newCoeff.RemoveAt(deg);
            Polynomial qx = new Polynomial(newCoeff);
            return (qx, pc);
        }

        /** 
         * ? Horner: nhân (multiple)
        ** function solving p(x)=(x-c)q(x)
        *!  input:
        **      + coeficient of q(x)
        **      + tuple (x-c) = (1,-c)
        *!  return q(x)*(x-c) =p(x)
        **/
        public static Polynomial operator*(Polynomial poly, (int x, double _c) x_c){
            List<double> result = new List<double>{};
            result = poly.ThisCoeffs;
            int index_Deg_0 = 0;
            result.Add(0);
            index_Deg_0 = result.Count-1;
            for(int i = index_Deg_0; i>0; i--){
                result[i]=result[i]-result[i-1]*(-x_c._c);
            }
            Polynomial px = new Polynomial(result);
            return px;
        }

        public double f_At(double x){
            return (this/(1,-x)).Pc;
        }

        /**
         *  ? Tính đạo hàm cấp "degree" của f tại x
         *  !<Áp dụng Hocner>
        **/
        public double f_dLevel_at(int degree, double x){
            Polynomial temp = new Polynomial(this.coeffs);
            int deg = temp.ThisCoeffs.Count-1;
            if(degree>deg){
                return 0;
            }else{
                for(int times = 1; times <= degree; times++){
                    temp = (temp/(1,-x)).Qx;
                }
                return temp.f_At(x);
            }
        }
        /**
         *  ? Tính tích phân của đa thức  
        **/
        public double Integral(double a, double b){
            int deg = ThisCoeffs.Count-1;
            double result = 0;
            if(ThisCoeffs.Count!=0){
                for(int i = 0; i <= deg; i++){
                    int newExp = deg - i + 1;
                    result += ThisCoeffs[i]*(Math.Pow(b,newExp)-Math.Pow(a,newExp))/(newExp);
                }
            }else{
                throw new Exception("\nERROR!!: this poly empty => Can't calculate integral");
            }
            return result;
        }

    }
}