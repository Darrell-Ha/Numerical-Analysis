
using System;
using System.Collections.Generic;
using Noi_suy_newton;
using Noi_suy_lagrange;
using PolynomialFunction;

namespace Noi_suy_nguoc
{
    public class InverseInterpolation{

        private const int countNumPoint = 4;
        private double EPSILON = 1e-7;
        private double _y;
        private List<Point> inputXY;
        private List<(int left, int right)> spaces;              // ! các khoảng song ánh và chứa nghiệm
        private List<(bool flag,(int left, int right))> spacesInterpol;     // ! các khoảng song ánh tương đương để nội suy nghiệm
        public InverseInterpolation(List<Point> inputXY, double _y){
            this.inputXY = inputXY;
            this._y = _y;
            this.spaces = new List<(int left, int right)>{};
            this.spacesInterpol = new List<(bool flag,(int left, int right))>{};
        }

        public List<double> Result_Interpolation(){
            List<double> result = new List<double>{};
            int countBij_spaces = 0;

            /**
             * TODO-1: Tìm khoảng song ánh và chứa nghiệm 
            **/
            ThisSpaces = Find_Bijection_has_x0(ThisY);
            countBij_spaces = ThisSpaces.Count;
            /**
             * TODO-2: Với các khoảng song ánh tương ứng, rút ra số mốc phù hợp
            **/
            if(countBij_spaces > 0){
                for(int i = 0; i < countBij_spaces; i++){
                    ThisSpacesInterpol.Add(Find_InputSpace(ThisSpaces[i],ThisY));
                }
            /**
             * TODO-3: Nội suy với từng khoảng song ánh tìm được 
            **/
                for(int i = 0; i < countBij_spaces; i++){
                    double x = Find_x_by_Interpolation(ThisSpacesInterpol[i]);
                    result.Add(x);
                }
            }else{
                throw new Exception("\n\tNot contain Bij-space. Program Stopped\n");
            }
            return result;
        }

        /**
         * TODO-1
        **/
        public List<(int left, int right)> Find_Bijection_has_x0(double _y){
            List<(int left, int right)> list = new List<(int left, int right)>{};
            (int left, int right) bij_space = (0,0);
            int signCurr = 0;
            int left = 0;
            int right = 1;
            int countPoint = ThisInputXY.Count;

            double divided_Yi1_Yi =ThisInputXY[1].ThisY - ThisInputXY[0].ThisY;
            divided_Yi1_Yi /=ThisInputXY[1].ThisX - ThisInputXY[0].ThisX;
            int sign = SignOf(divided_Yi1_Yi); 

            for(int i = 1; i < countPoint-1; i++){
                divided_Yi1_Yi =ThisInputXY[i+1].ThisY - ThisInputXY[i].ThisY;
                divided_Yi1_Yi /=ThisInputXY[i+1].ThisX - ThisInputXY[i].ThisX;
                signCurr = SignOf(divided_Yi1_Yi);
                if(signCurr!= sign){
                    right = i;
                    bij_space = (left,right);
                    if((sign!=0) && IsIn(bij_space,_y)){
                        list.Add(bij_space);
                    }
                    left = i;
                    sign = signCurr;
                }
                if(i==countPoint-2){
                    right = i+1;
                    bij_space = (left,right);
                    list.Add(bij_space);
                }
            }
            return list;
        }

        /**
         * TODO-2
        **/
        public (bool flag,(int left, int right)) Find_InputSpace((int left, int right) space,double _y){
            bool flag = true;
            int left = space.left;
            int right = space.right;
            int indexMin = left;
            int count = right - left + 1;
            double minDist = Math.Abs(_y - ThisInputXY[left].ThisY);
            double tempMin = 0;

            if(count > countNumPoint){
                for(int i = left + 1; i <= right; i++){
                    tempMin = Math.Abs(_y - ThisInputXY[i].ThisY);
                    if(tempMin < minDist){
                        minDist = tempMin;
                        indexMin = i;
                    }
                }
                if(right - indexMin + 1 <= countNumPoint){
                    left = indexMin;
                }else if(indexMin - left + 1 <= countNumPoint){
                    right = indexMin;
                    flag = false;
                }else{
                    left = indexMin; 
                    right = indexMin + countNumPoint - 1;
                }
            }

            return (flag,(left,right));
        }

        public double Find_x_by_Interpolation((bool flag,(int left, int right) space) spaceInput){
            double x = 0;
            (int left, int right) space = (spaceInput.space.left,spaceInput.space.right);
            if(Is_StableStep(space)){
                if(spaceInput.flag){
                    x = Method_Stable_Forward(space,ThisY);
                }else{
                    x = Method_Stable_Backward(space,ThisY);
                }
            }else{
                x = Method_NotStable(space,ThisY);
            }
            return x;
        }

        public double Method_NotStable((int left, int right) space, double _y){
            double x = 0;
            Polynomial poly = new Polynomial();
            List<Point> inputYX = ConvertToList(0,space);
            Newton process = new Newton(inputYX,0);
            // Lagrange process = new Lagrange(inputYX);
            /**
             * !! x = \phi(y) 
            **/
            poly = process.Interpolation();
            x = poly.f_At(_y);              
            return x;
        }

        public double Method_Stable_Forward((int left, int right) space, double _y){
            List<Point> inputXY = ConvertToList(1,space); 
            Newton process = new Newton(inputXY,1);
            Polynomial eta = process.Interpolation();
            Polynomial R = new Polynomial(new List<double>{process.ThisDiff[1][0],process.ThisInputXY[0].ThisY});
            eta = eta - R;
            double h = process.ThisStep;
            double deltaY0 = process.ThisDiff[1][0];
            double t_0 = _y - process.ThisInputXY[0].ThisY;          //! t_0 = (*y - y_0)
            double t = 0;       
            double x_0 = process.ThisInputXY[0].ThisX;
            double x_bef = x_0;
            double x_curr = 0;
            double eps = 99;

            while(eps >= ThisEpsilon){ 
                t  = (t_0 - eta.f_At(t))/deltaY0;
                x_curr = x_0 + h*t;
                eps = Math.Abs((x_curr-x_bef)/x_curr);
                x_bef = x_curr;  
            }
            
            return x_curr;
        }
        public double Method_Stable_Backward((int left, int right) space, double _y){
            List<Point> inputXY = ConvertToList(1,space); 
            Newton process = new Newton(inputXY,-1);
            Polynomial eta = process.Interpolation();
            int indexNabla = process.ThisDiff[1].Count-1;
            int countPoint = inputXY.Count;
            double nablaY0 = process.ThisDiff[1][indexNabla];
            Polynomial R = new Polynomial(new List<double>{nablaY0,process.ThisInputXY[countPoint-1].ThisY});
            eta = eta - R;
            double h = process.ThisStep;
            double t_0 = _y - process.ThisInputXY[countPoint-1].ThisY;          //! t_0 = (*y - y_0)
            double t = 0;       
            double x_0 = process.ThisInputXY[countPoint-1].ThisX;
            double x_bef = x_0;
            double x_curr = 0;
            double eps = 99;

            while(eps >= ThisEpsilon){ 
                t  = (t_0 - eta.f_At(t))/nablaY0;
                x_curr = x_0 + h*t;
                eps = Math.Abs((x_curr-x_bef)/x_curr);
                x_bef = x_curr;  
            }
            
            return x_curr;
        }
        private int SignOf(double num){
            if(num > 0) return 1;
            else if(num < 0) return -1;
            else return 0;
        }

        private bool IsIn((int left, int right) space, double _y){
            bool res1 = (ThisInputXY[space.left].ThisY < _y) && (_y < ThisInputXY[space.right].ThisY);
            bool res2 = (ThisInputXY[space.left].ThisY > _y) && (_y > ThisInputXY[space.right].ThisY);
            return res1 || res2;
        }


        public bool Is_StableStep((int left, int right) space){
            bool result = true;
            int left = space.left;
            int right = space.right;
            double step = ThisInputXY[right].ThisX - ThisInputXY[right-1].ThisX; 
            double tempStep = 0;
            for(int i = left; i < right; i++){
                tempStep = ThisInputXY[i+1].ThisX - ThisInputXY[i].ThisX;
                if(Math.Abs(tempStep-step)>=1e-15){
                    result = false;
                    Console.WriteLine("ERROR!!: Unstable");
                    break;
                }
            }
            return result;
        }

        public List<Point> ConvertToList(int flag, (int left, int right) space){
            List<Point> result = new List<Point>{};
            int a = space.left;
            int b = space.right;
            for(int i = a; i <= b; i++){
                if(flag == 1){                      ///! Cách đều
                    result.Add(ThisInputXY[i]);
                }else{
                    Point point = new Point(ThisInputXY[i].ThisY, ThisInputXY[i].ThisX);
                    result.Add(point);
                }
            }
            return result;
        }
        public List<Point> ThisInputXY{
            get{return this.inputXY;}
            set{this.inputXY = value;}
        }

        public double ThisY{
            get{return this._y;}
            set{this._y = value;}
        }
        public double ThisEpsilon{
            get{return this.EPSILON;}
            set{this.EPSILON = value;}
        }

        public List<(int left, int right)> ThisSpaces{
            get{return this.spaces;}
            set{this.spaces = value;}
        }
        public List<(bool flag,(int left, int right))> ThisSpacesInterpol{
            get{return this.spacesInterpol;}
            set{this.spacesInterpol = value;}
        }

    }
}