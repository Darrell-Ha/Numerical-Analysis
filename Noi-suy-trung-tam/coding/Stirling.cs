using System.Collections.Generic;
using System;
using PolynomialFunction;
namespace Noi_suy_trung_tam
{

    public class Stirling{
       
        private List<Point> inputXY;
        private Polynomial P_n;
        private Polynomial omega;
        private List<List<double>> diffs;
        private int deg; 
        private int status;
        private int index_zero = -1;
        private double factorial = 1;
        private double step;

        public Stirling(List<Point> inputXY){
            this.deg        = inputXY.Count-1;
            this.index_zero = deg/2;
            this.inputXY    = inputXY;
            this.diffs      = new List<List<double>>{new List<double>{ThisInputXY[ThisIndexZero].ThisY}};
            this.P_n        = new Polynomial(new List<double>{ThisInputXY[ThisIndexZero].ThisY});        // * P_n = y_0
            this.omega      = new Polynomial(new List<double>{1});                                       // * poly omega = 1;
            this.step       = inputXY[ThisDeg].ThisX - inputXY[ThisDeg-1].ThisX;                         // ! Caution
        }


        /**
         * ? Nguyên tắc áp dụng: (Gauss1 + Gauss2)/2 (dùng sai phân tiến là chính)
         *  * Gauss1 : y_0 => y_1 => y_-1 => y_2 =>....
         *          ! bắt đầu bằng Status = 1 
         *  * Gauss2 : y_0 => y_-1 => y_1 => y_-2=>....
         *          ! bắt đầu bằng Status = -1 
         * 
         * /------------------------------------------------------------/
         * ? Cụ thể:
         * * Check tính hợp lệ:
         *          ? Số điểm có là lẻ (2n+1)?
         *          ? Không bị lặp điểm? 
         *          ? Mốc cách đều ?
         * * Khởi tạo sẵn
         *      ! index_zero = [(2n+1)/2]
         *      ! factorial = 1
         *      ! P_n = y_[index_zero]
         *      ! omega = 1
         * * Duyệt input x_k và x_{-k}: k từ 1 -> n 
         *      * Thêm 2 điểm mới và Cập nhật bảng sai phân:
         *              ? điểm (x_k,y_k) với Status = 1 (tiến)
         *                  ! Thêm giá trị sai phân mới vào cuối cột sai phân từ cấp 0 (giá trị y_k) đến cấp 2k-1
         *              ? điểm (x_{-k},y_{-k}) với Status = -1 (lùi)
         *                  ! Thêm giá trị sai phân mới vào đầu cột sai phân từ cấp 0 (giá trị y_{-k}) đến cấp 2k
         *      * Sau khi cập nhật
         *              ? <Lấy giá trị của ô sai phân đầu tiên của cột cấp 2k>
         *                      ! coef_even = ThisDiff[2k][0]
         *              ? <Lấy trung bình của 2 sai phân đầu tiên của cột cấp 2k-1>
         *                      ! coef_odd = (ThisDiff[2k-1][0] + ThisDiff[2k-1][1])/2 
         * ???? <GHI CHÚ> Số cột = số điểm nạp vào ??? 
         *      * Tính: omega *= (t^2 - (k-1)^2) 
         * 
         *      * Cập nhật giá trị factorial*=(2k-1)
         *      * Tính: P_odd = (coef_odd/factorial) * omega/t
         * 
         *      * Cập nhật giá trị factorial*=2k
         *      * Tính: P_even = (coef_even/factorial) * omega
         * 
         *      * Tính: P_n += P_odd + P_even
         *  /////////////////////////////////////////////////////
         * * Kết thúc duyệt: Thu được đa thức P_n: return P_n;
         * 
        **/
        public virtual Polynomial Interpolation(){
            
            double coef_odd = 0;
            double coef_even = 0;
            Polynomial P_odd = new Polynomial();
            Polynomial P_even = new Polynomial();
            
            if((ThisDeg+1)%2!=0 && !Is_DuplicateXY() && Is_StableStep()){       // ? Số điểm có là số lẻ??

                for(int k = 1; k <= ThisDeg/2; k++){
                    
                    // ? Thêm điểm (x_k,y_k)
                    ThisStatus = 1;
                    UpdateDiff(ThisInputXY[ThisIndexZero+k]);
                    
                    // ? Thêm điểm (x_{-k},y_{-k})
                    ThisStatus = -1;
                    UpdateDiff(ThisInputXY[ThisIndexZero-k]);

                    // ? Lấy hệ số cho đa thức bậc lẻ và chẵn
                    coef_odd = (ThisDiff[2*k-1][0] + ThisDiff[2*k-1][1])/2;
                    coef_even = ThisDiff[2*k][0];

                    // ! Tính Omega*=(t^2-(k-1)^2)
                    ThisOmega *= (1,-(k-1)); 
                    ThisOmega *= (1,(k-1)); 

                    // ? Tính đa thức bậc lẻ
                    ThisFactorial *= (2*k-1);
                    P_odd = (coef_odd/factorial)*((ThisOmega/(1,0)).Qx); // ! lưu ý: Qx = omega/t
                    
                    // ? Tính đa thức bậc chẵn
                    ThisFactorial *= 2*k;
                    P_even = (coef_even/factorial)*ThisOmega;

                    // ! Cập nhật P_n
                    ThisP_n += P_odd + P_even;
                }

            }else{
                Console.WriteLine("ERROR!!: Need Odd Points!!");
            }
            return ThisP_n;
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
                // if(ThisStatus == 0){
                //     tempRate /= ThisInputXY[countPoint_now-1].ThisX - ThisInputXY[countPoint_now-1-k].ThisX;
                // }
                if(ThisStatus!=-1){
                    ThisDiff[k].Add(tempRate);
                }else{
                    ThisDiff[k].Insert(0,tempRate);
                }
                
            }

        }
        public bool Is_StableStep(){
            bool result = true;
            double tempStep = 0;
            for(int i = 0; i < ThisDeg; i++){
                tempStep = ThisInputXY[i+1].ThisX - ThisInputXY[i].ThisX;
                if(Math.Abs(tempStep-ThisStep)>=1e-15){
                    result = false;
                    Console.WriteLine("ERROR!!: Unstable");
                    break;
                }
            }
            return result;
        }
        public bool Is_DuplicateXY(){
            bool result = false;
            if(ThisInputXY.Count!=0){
                for(int i = 0; i < ThisDeg; i++){
                    double x_i = ThisInputXY[i].ThisX; 
                    for(int k = i+1; k <= ThisDeg; k++){
                        if(x_i == ThisInputXY[k].ThisX){
                            result = true;
                            Console.WriteLine($"ERROR!!: Duplicated at {i+1} and {k+1}: {x_i}");
                            break;
                        }
                    }
                }
            }
            return result;
        }
        public double Exchange_x_To_t(double x){
            double t = 0;
            if(ThisIndexZero!=-1){
                t = (x-ThisInputXY[ThisIndexZero].ThisX)/ThisStep;
            }else{
                Console.WriteLine("ERROR!!!");
            }
            return t;
        }
        public List<Point> ThisInputXY{
            get{return this.inputXY;}
            set{this.inputXY = value; this.deg = value.Count-1;}
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
        public int ThisIndexZero{
            get{return this.index_zero;}
            set{this.index_zero = value;}
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