using System;
using System.Collections.Generic;
using PolynomialFunction;
namespace Tinh_tich_phan
{
    public class Simpson : CalcIntegral{

        public Simpson(Function fx, double a, double b, int n){ // ? MAIN 1: INPUT f,a,b,n
            this.fx      = fx;
            this.a       = a;
            this.b       = b;
            this.n       = n;
            this.Mn      = ThisFx.FindMax_Fn(4,a,b);
            this.epsilon = ConvertN_toEpsilon(ThisN);
            this.inputXY = DivideAB_byN(ThisN);
            this.PRINT_TO_SCREEN = ThisEpsilon;
        }
        public Simpson(Function fx, double a, double b, double eps){    // ? MAIN 2: INPUT f,a,b,eps
            this.fx      = fx;
            this.a       = a;
            this.b       = b;
            this.epsilon = eps;
            this.Mn      = ThisFx.FindMax_Fn(4,ThisA,ThisB) ;
            this.n       = ConvertEpsilon_toN(ThisEpsilon);
            this.inputXY = DivideAB_byN(ThisN);
            this.PRINT_TO_SCREEN = ThisEpsilon;
        }
        public Simpson(Function fx, double a, double b){
            this.fx         = fx;
            this.a          = a;
            this.b          = b;
            this.n          = (int)(100000*(b-a));
            this.n          = (ThisN%2==0)? ThisN+2:ThisN+1;
            this.Mn         = ThisFx.FindMax_Fn(4,a,b);
            this.epsilon    = ConvertN_toEpsilon(ThisN);
            this.inputXY    = DivideAB_byN(n);
        }
        public Simpson(List<Point> list, double step){
            ThisInputXY = list;
            ThisStep = step;
            this.n = list.Count-1;
        }
        /**************************************************************************************************** 
         * ? Ý tưởng: Giống hình thang nhưng
         *            !xấp xỉ đa thức cấp 2 đi qua 3 điểm (x_{2k},x_{2k+1},x_{2k+2})
         *     * I_k = h/3 * (y_{2k+2}+4y_{2k+1}+y_{2k}) 
         * 
         * 
         * ? Công thức tích phân Simpson
         *      ! I_n = h/3 * (y_0+y_{2k} + 2(y_2+y_4+...+y_{2n-2}) + 4(y_1+y_3+....+y_{2n-1}))
         *   
         * ? Cụ thể:
         *  *INPUT  : fx, (a,b), epsilon
         *  *OUTPUT : I_{2n}
         * TODO >>>>
         *      * Tính giá trị M4 = max|f^(4)(x)| trên (a,b)   
         *          ? (Thông qua Gradient Descent)
         *      * Tính giá trị 2n để chia [a,b] thành 2n khoảng
         *          ? |I-I_{2n}| <= M_4/90 * (b-a)^5/n^4  < epsilon
         *              ! => 2n = [sqrt[4](M4 * (b-a)^5/(90eps))] + 1:2 để thành chẵn 
         *      * Tính h = (b-a)/(2n)
         *      * Sau khi chia n khoảng có độ dài h
         *          ? [a,b] chứa các điểm (x_0 = a < x_1 < ... < x_{2n-1} < x_{2n} = b)
         *              ! Tính tọa độ (x_i,y_i) = (x_i, f(x_i))
         *      * Có được 2n+1 điểm (x_i,y_i) (i=0,1,...2n)
         * 
         *      ! return result = h/3 * (y_0+y_{2k} + 2(y_2+y_4+...+y_{2n-2}) + 4(y_1+y_3+....+y_{2n-1}))
        *******************************************************************************************************/
        public override double Calc_Integral(){
            
            double result = 0;
            if(ThisInputXY!=null){
                result += ThisInputXY[0].ThisY + ThisInputXY[ThisN].ThisY;
                for(int i = 1; i < ThisN; i++){
                    if(i%2==0){
                        result += 2*ThisInputXY[i].ThisY;
                    }else{
                        result += 4*ThisInputXY[i].ThisY;
                    }
                }
                result *= ThisStep/3;
                return result;
            }else{
                return -1;
            }
        }
        /**
         * ? Công thức sai số của tích phân Simpson
         *  ! |I-I_{2n}| <= M_4/90 * (b-a)^5/n^4  < epsilon 
        **/
        public override int ConvertEpsilon_toN(double eps){
            int n = 0;
            n = (int)Math.Exp(0.25*Math.Log(ThisMn*Math.Pow(ThisB-ThisA,5)/(90*ThisEpsilon)));
            n = (n%2==0)? n+2 : n+1;
            if(Math.Abs(n)>=1e8){
                throw new Exception($"\nThisMn = {ThisMn}\nn = {n} too bigg => infty");
            }
            return n;
        }
        public override double ConvertN_toEpsilon(int n){
            double eps = 0;
            eps = ThisMn*Math.Pow(ThisB-ThisA,5)/(90*Math.Pow(n,4));
            if(Math.Abs(n)>=1e8){
                throw new Exception($"\nn = {n} too bigg => infty");
            }
            return eps;
        }
    }


}