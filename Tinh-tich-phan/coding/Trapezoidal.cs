
using System;
using System.Collections.Generic;
using PolynomialFunction;

namespace Tinh_tich_phan
{
    public class Trapezoidal : CalcIntegral{

        public Trapezoidal(Function fx, double a, double b, int n){     // ? MAIN 1: INPUT f,a,b,n
            this.fx      = fx;
            this.a       = a;
            this.b       = b;
            this.n       = n;
            this.Mn      = fx.FindMax_Fn(2,a,b);
            this.epsilon = ConvertN_toEpsilon(n);
            this.inputXY = DivideAB_byN(n);
            this.PRINT_TO_SCREEN = ThisEpsilon;
        }
        public Trapezoidal(Function fx, double a, double b, double eps){    // ? MAIN 2: INPUT f,a,b,eps
            this.fx         = fx;
            this.a          = a;
            this.b          = b;
            this.epsilon    = eps;
            this.Mn         = ThisFx.FindMax_Fn(2,ThisA,ThisB);     //! Tìm M2
            this.n          = ConvertEpsilon_toN(eps);
            this.inputXY    = DivideAB_byN(ThisN);
            this.PRINT_TO_SCREEN = ThisEpsilon;
        }
        public Trapezoidal(Function fx, double a, double b){
            this.fx         = fx;
            this.a          = a;
            this.b          = b;
            this.n          = (int)(1000000*(b-a));                 ///! HUGEEE
            // this.Mn         = fx.FindMax_Fn(2,a,b);             /// ! just test, not compulsory
            this.epsilon    = ConvertN_toEpsilon(ThisN);
            this.inputXY    = DivideAB_byN(ThisN);
        }
        // public Trapezoidal(List<Point> inputXY, double a, double b){    // ! Not completed
        //     this.inputXY = inputXY;
        //     this.a       = a;
        //     this.b       = b;
        //     this.n       = ThisInputXY.Count;

        //     this.epsilon = ConvertN_toEpsilon(n);
        // }
        /**************************************************************************************************** 
         * ? Ý tưởng: xấp xỉ đa thức cấp 1 đi qua 2 điểm (x_k,x_{k+1})
         *     * I_k = (y_k+y_{k+1})/2 (diện tích hình thang vuông có 2 đáy có độ dài y_k và y_{k+1})
         * 
         * 
         * ? Công thức tích phân hình thang
         *      ! I_n = h(0.5(y_0+y_n)+y_1+y_2+...+y_{n-1}) 
         *   
         * ? Cụ thể:
         *  *INPUT  : fx,(a,b), epsilon
         *  *OUTPUT : I_n
         * TODO
         *      ! Tính giá trị M2 = max|f"(x)| trên (a,b)   
         *          ? (Thông qua Gradient Descent)
         *      * Tính giá trị n để chia [a,b] thành n khoảng
         *          ? |I-I_n| <= M_2/2 * (b-a)^3/n^2  < epsilon
         *              ! => n = [sqrt(M2 * (b-a)^3/(12eps))] + 1 
         *      * Tính h = (b-a)/n
         *      * Sau khi chia n khoảng có độ dài h
         *          ? [a,b] chứa các điểm (x_0 = a < x_1 < ... < x_{n-1} < x_n = b)
         *              ! Tính tọa độ (x_i,y_i) = (x_i, f(x_i))
         *      * Có được n+1 điểm (x_i,y_i) (i=0,1,...n)
         * 
         *      ! return result = h(0.5(y_0+y_n)+y_1+y_2+...+y_{n-1})
        *******************************************************************************************************/
        public override double Calc_Integral(){
            
            double result = 0;
            if(ThisInputXY!=null){

                for(int i = 1; i < ThisN; i++){
                    result +=  ThisInputXY[i].ThisY;
                }
                result += 0.5*(ThisInputXY[0].ThisY + ThisInputXY[ThisN].ThisY);
                result *= ThisStep;
                return result;
            }else{
                return -1;
            }
        }
        /**
         * ? Công thức sai số của tích phân hình thang
         *  ! |I-I_n| <= M_2/12 * (b-a)^3/n^2 
        **/
        public override int ConvertEpsilon_toN(double eps){
            int n = 0;
            n = (int)Math.Sqrt(ThisMn*Math.Pow(ThisB-ThisA,3)/(12*ThisEpsilon)) + 1;
            if(Math.Abs(n)>=1e8){
                throw new Exception($"\nThisMn = {ThisMn}\nn = {n} too bigg => infty");
            }
            return n;
        }
        public override double ConvertN_toEpsilon(int n){
            double eps = 0;
            eps = ThisMn*Math.Pow(ThisB-ThisA,3)/(12*n*n);
            if(Math.Abs(n)>=1e8){
                throw new Exception($"\nn = {n} too bigg => infty");
            }
            return eps;
        }
    }
}