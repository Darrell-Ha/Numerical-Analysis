using System;
using System.Collections.Generic;
using PolynomialFunction;

namespace Noi_suy_newton
{
    public class Newton{
        private List<Point> inputXY;
        private Polynomial P_n;
        private Polynomial omega;
        private List<List<double>> diffs;
        private int deg; 
        private int status;
        private double factorial = 1;
        private double step;

        public Newton(List<Point> inputXY, int status = 0){
            this.status = status;
            this.deg = inputXY.Count-1;
            this.inputXY = inputXY;
            this.diffs = new List<List<double>>{};
            this.P_n = new Polynomial();        // * P_n = 0
            this.omega = new Polynomial(new List<double>{1});  // * poly omega = 1;
            this.step = inputXY[ThisDeg].ThisX - inputXY[ThisDeg-1].ThisX; //! Caution
        }

        public Newton(int status){
            this.status = status;
            this.diffs = new List<List<double>>{};
            this.P_n = new Polynomial();
        }

        public Polynomial Interpolation(){
            
            for(int i = 0; i <= ThisDeg; i++){
                
                if(status==0 && !Is_DuplicateXY(ThisInputXY[i])){
                    Newton_Base(ThisInputXY[i],i);
                }else if(status == 1 && !Is_DuplicateXY(ThisInputXY[i]) && Is_StableStep(i-1,ThisInputXY[i])){
                    Newton_Forward(ThisInputXY[i],i);
                }else if(status == -1 && !Is_DuplicateXY(ThisInputXY[ThisDeg-i]) && Is_StableStep(ThisDeg-i+1,ThisInputXY[ThisDeg-i])){
                    Newton_Backward(ThisInputXY[ThisDeg-i],i);
                }else{
                    break;
                }
            
            }
            return ThisP_n;
        }

        /**
         * ! Nguyên tắc áp dụng
         *  * Kiểm tra các điều kiện (lặp điểm x, thỏa mãn mốc cách đều,....)
         *  * Thêm điểm -> Update bảng sai phân/tỷ sai phân
         *  * Cập nhật Đa thức nội suy 
         *      * P_n = P_n + hệ số * Omega;
         *          ? Mốc bất kì:
         *                      * x = x 
         *                      ! P_n = y_0 + y[x_0,x_1](x-x_)+...+y[x_0,...,x_n](x-x_0)(x-x_1)...(x-x_{n-1})
         *          ? Mốc cách đều: 
         *              ? Tiến: 
         *                      * x = x_0 + th
         *                      ! P_n = y_0 + (\Delta y_0)/1! (t)+ (\Delta^2 y_0)/2! t(t-1)....(\Delta^n y_0)/n! t(t-1)...(t-n+1)
         *              ? Lùi: 
         *                      * x = x_n + th
         *                      ! P_n = y_n + (\Nelta y_n)/1! (t)+ (\Nabla^2 y_n)/2! t(t-1)....(\Nabla^n y_n)/n! t(t+1)...(t+n-1)
         *         
        **/
        public void Newton_Base(Point newPoint, int count){
            UpdateDiff(newPoint);
            if(count == 0){
                ThisP_n += ThisDiff[0][0];
            }else{
                ThisOmega *= (1,-ThisInputXY[count-1].ThisX);
                ThisP_n += ThisDiff[count][0]*ThisOmega;
            }
            // Console.WriteLine(ThisP_n.ToString()+"\n\n\n\n");
        }

        public void Newton_Forward(Point newPoint, int count){
            UpdateDiff(newPoint);
            if(count == 0){
                ThisP_n += ThisDiff[0][0];
            }else{
                ThisFactorial *= count;
                ThisOmega *= (1,-(count-1));
                ThisP_n += (ThisDiff[count][0]/ThisFactorial)*ThisOmega;
            }

        }

        public void Newton_Backward(Point newPoint, int count){
            UpdateDiff(newPoint);
            if(count == 0){
                ThisP_n += ThisDiff[0][0];
            }else{
                int index = ThisDiff[count].Count-1;
                ThisFactorial *= count;
                ThisOmega *= (1,(count-1));
                ThisP_n += (ThisDiff[count][index]/ThisFactorial)*ThisOmega;
            }
        }

        public void UpdateDiff(Point newPoint){

            int lastIndex = 0;
            double tempRate = 0;
            List<double> newColDiff = new List<double>{};
            ThisDiff.Add(newColDiff);
            if(ThisStatus!=-1){
                ThisDiff[0].Add(newPoint.ThisY);
            }else{
                ThisDiff[0].Insert(0,newPoint.ThisY);
            }
            int countPoint_now = ThisDiff[0].Count;
            // Console.WriteLine($"{countPoint_now}");
            for(int k = 1; k < countPoint_now; k++){
                lastIndex = (status!=-1)? ThisDiff[k-1].Count-1 : 1;
                tempRate = ThisDiff[k-1][lastIndex]-ThisDiff[k-1][lastIndex-1];
                if(ThisStatus == 0){
                    tempRate /= ThisInputXY[countPoint_now-1].ThisX - ThisInputXY[countPoint_now-1-k].ThisX;
                }
                if(ThisStatus!=-1){
                    ThisDiff[k].Add(tempRate);
                }else{
                    ThisDiff[k].Insert(0,tempRate);
                }
                
            }

        }
        public bool Is_StableStep(int index_lastNewPoint, Point newPoint){
            bool result = true;
            if(index_lastNewPoint>=0 && index_lastNewPoint <= ThisDeg){
                result = (ThisStep-Math.Abs(ThisInputXY[index_lastNewPoint].ThisX-newPoint.ThisX) < 1e-10);
            }
            if(!result){
                Console.WriteLine("ERROR!!: Unstable");
            }
            return result;
        }
        public bool Is_DuplicateXY(Point newPoint){
            bool result = false;
            if(ThisDiff.Count!=0){
                for(int i = 0; i<ThisDiff[0].Count; i++){
                    if((ThisStatus!=-1 && ThisInputXY[i].ThisX==newPoint.ThisX) || (ThisStatus==-1 && ThisInputXY[ThisDeg-i].ThisX==newPoint.ThisX)){
                        result = true;
                        Console.WriteLine($"ERROR!!: Duplicated at {i+1}: {newPoint.ThisX}");
                        
                        break;
                    }
                }
            }
            return result;
        }

        public double Exchange_x_To_t(double x){
            double t = 0;
            if(ThisStatus==0){
                t = x;
            }else if(ThisStatus>0){
                t = (x-ThisInputXY[0].ThisX)/ThisStep;
            }else{
                t = (x-ThisInputXY[ThisDeg].ThisX)/ThisStep;
            }
            return t;
        }

        public List<Point> ThisInputXY{
            get{return this.inputXY;}
            set{this.inputXY=value; this.deg=value.Count-1;}
        }

        public int ThisDeg{
            get{return this.deg;}
            set{this.deg = value;}
        }

        public Polynomial ThisP_n{
            get{return this.P_n;}
            set{this.P_n = value;}
        }
        public Polynomial ThisOmega{
            get{return this.omega;}
            set{this.omega = value;}
        }
        public int ThisStatus{
            get{return this.status;}
            set{this.status = value;}
        }
        public double ThisStep{
            get{return this.step;}
            set{this.step = value;}
        }
        public double ThisFactorial{
            get{return this.factorial;}
            set{this.factorial = value;}
        }
        public List<List<double>> ThisDiff{
            get{return this.diffs;}
            set{this.diffs=value;}
        }
    }
}