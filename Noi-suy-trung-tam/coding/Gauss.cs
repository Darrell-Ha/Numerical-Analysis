
using System;
using System.Collections.Generic;
using PolynomialFunction;
namespace Noi_suy_trung_tam
{
    public class Gauss {

        private List<Point> inputXY;
        private Polynomial P_n;
        private Polynomial omega;
        private List<List<double>> diffs;
        private int deg; 
        private int status;
        private int index_zero=-1;
        private double factorial = 1;
        private double step;

        public Gauss(List<Point> inputXY, int status = 1){
            this.status = status;
            this.deg = inputXY.Count-1;
            this.index_zero = deg/2;
            this.inputXY = inputXY;
            this.diffs = new List<List<double>>{new List<double>{ThisInputXY[ThisIndexZero].ThisY}};
            this.P_n = new Polynomial(new List<double>{ThisInputXY[ThisIndexZero].ThisY});        // * P_n = y_0
            this.omega = new Polynomial(new List<double>{1});  // * poly omega = 1;
            this.step = inputXY[ThisDeg].ThisX - inputXY[ThisDeg-1].ThisX; //! Caution
        }


        /**
         * ? Nguyên tắc áp dụng: Newton tiến (dùng sai phân tiến là chính)
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
         *          ? Cách đều mốc?
         * * Khởi tạo sẵn
         *      ! index_zero = [(2n+1)/2]
         *      ! factorial = 1
         *      ! P_n = y_[index_zero]
         *      ! omega = 1
         *      ! index_curr = index_zero (vị trí xuất phát)
         * * Duyệt input x_k: k từ 1 -> 2n 
         *      * Cập nhật giá trị factorial*=k
         *      * Cập nhật vị trí điểm cần nạp hiện tại 
         *              ! index_curr += ThisStatus * k
         *      * Thêm điểm mới và cập nhật bảng sai phân:
         *              ? Nếu đang ở trạng thái Tiến: (Status = 1)
         *                  ! Thêm giá trị sai phân mới vào cuối cột sai phân từ cấp 0 (giá trị y_k) đến cấp k
         *              ? Nếu đang ở trạng thái Lùi:  (Status = -1)
         *                  ! Thêm giá trị sai phân mới vào đầu cột sai phân từ cấp 0 (giá trị y_k) đến cấp k
         *      * Sau khi cập nhật
         *              ! Lấy giá trị của ô sai phân đầu tiên của cột mới nhất = ThisDiff[k][0]
         *              ? <GHI CHÚ> Số cột = số điểm nạp vào 
         *      * Tính: omega *= (t - lastDistance) (với lastDistance = 0,+-1,+-2,...,+-n)
         *      * Tính: P_n += (ThisDiff[k][0]/factorial)*omega
         *      *
         *      * lưu lastDistance = index_curr - index_zero
         *      * Chuyển trạng thái status *= -1
         *  
         * * Kết thúc duyệt: Thu được đa thức P_n: return P_n;
         * 
        **/
        public virtual Polynomial Interpolation(){
            
            int index_curr = ThisIndexZero;
            int lastDistance = 0;
            if((ThisDeg+1)%2!=0 && !Is_DuplicateXY() && Is_StableStep()){       // ? Số điểm có là số lẻ??

                for(int k = 1; k <= ThisDeg; k++){

                    ThisFactorial *= k;
                    index_curr += ThisStatus*k;
                    UpdateDiff(ThisInputXY[index_curr]);

                    ThisOmega *= (1,-lastDistance);
                    ThisP_n += (ThisDiff[k][0]/ThisFactorial)*ThisOmega;

                    // ! Chuyển trạng thái tiến -> lùi hoặc ngược lại
                    ThisStatus*=-1;
                    lastDistance = index_curr - ThisIndexZero;
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