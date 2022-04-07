using System;
using System.Collections.Generic;
using PolynomialFunction;

namespace Tinh_tich_phan
{
    public class NewtonCotez : CalcIntegral{

        private List<double> cotezCoeff;

        public NewtonCotez(Function fx, double a, double b, int n){
            this.fx = fx;
            this.a = a;
            this.b = b;
            this.n = n;
            this.step = (b-a)/n;
            this.inputXY = DivideAB_byN(ThisN);
            this.cotezCoeff = GetCotezCoeff(ThisN);
        }
        
        /**
         * ? Newton - Cotez: Xây dựng đa thức nội suy trên toàn bộ miền, sau đó tính tích phân 
         * 
         *  * Công thức tính: I = h(A_0.y_0 + A_1.y_1 + ... + A_n.y_n)
         *  * Biến đổi: với x = x_0 + th
         *  *  A_i = int^n_0 \prod_{k!=i,k=0,1,..} (t-k)/(i-k)  (i = 0,1,2,...,n | k!=i)
         *          ! Tính chất A_i = A_{n-i} => chỉ cần tính một nửa rồi lấy đối xứng
         *  ?#################################################################################
         *  ?#                                                                               #
         *  ?! ### NEEDTODO ### =>>> Tính đc A_i ------>  A_i*y_i = I_i ------> I = \sum I_i #
         *  ?#                                                                               #
         *  ?#################################################################################
         * ? Cụ thể:
         *  TODO-1: Xác định y_i, E_i = \prod_{k!=i,k=0,1,..}(i-k), Q_i = \prod_{k!=i,k=0,1,..}(t-k)
         *      * Với Omega(t) = (t-0)(t-1)(t-2)...(t-n)           (lưu vào list hệ số coef của poly omega)
         *          ! Q_i(t) = Omega/(t-i)                         (lưu vào list hệ số coef của poly Q_i) 
         *          ! E_i = Q_i(i)
         *      * y_i = f(x_i) 
         *  TODO-2: Tính A_i = 1/E_i * (tích phân của Q_i từ 0 đến n =) int^n_0 Omega(t)/(t-i)
         *      * Hàm tính tích phân của đa thức Q_i = Q_i.Integral(0,n) (tùy có thể trình bày thêm)
         *              ! A_i = Q_i.ThisCoeffs[u]*(n^{deg-u+1})/(deg-u+1) (u = 0,1...,degQ_i)
         *              ! A_i = A_i/E_i         (Hệ số A_i lưu vào list)
         *  TODO-3: (tính được I_i = y_i * A_i)    
         *  TODO-4: I_n = h(I_0+I_1+...+I_n)     
         *    
        **/

        public List<double> GetCotezCoeff(int n){
            List<double> A = new List<double>{};
            Polynomial Omega = new Polynomial()+1;
            
            int mid = n/2;

            for(int i = 0; i <= n; i++){
                Omega = Omega*(1,-i);       
            }

            for(int i = 0; i <= mid; i++){
                Polynomial Q_i = (Omega/(1,-i)).Qx;
                double e_i = Q_i.f_At(i);
                double A_i = (1/e_i) *Q_i.Integral(0,n);

                A.Insert(i,A_i);
                if(i == mid && n%2 == 0) continue;
                A.Insert(A.Count-1-i,A_i);
            }
            return A;
        }

        public override double Calc_Integral(){
            double result = 0;

            for(int i = 0; i <= n; i++){
                result += ThisInputXY[i].ThisY*ThisCotezCoeff[i];
            }
            result*=ThisStep;
            return result;
        }
        public void PrintCotezCoeff(){
            List<double> list = ThisCotezCoeff;
            int countA = list.Count;
            double sum =0;
            for(int i = 0; i < countA; i++){
                Console.WriteLine($"A_({i}) = {list[i]}");
                sum += list[i];
            }
            Console.WriteLine($"\nsum = {sum}");
        }
        public override int ConvertEpsilon_toN(double eps){
            return 0;
        }
        public override double ConvertN_toEpsilon(int n){
            return 0;
        }

        public List<double> ThisCotezCoeff{
            get{return this.cotezCoeff;}
            set{this.cotezCoeff = value;}
        }
    }
}