using System;
using System.Collections.Generic;
using PolynomialFunction;

namespace Chebysev{
    
    public class NormalizeChebysev{
        private double a;
        private const double Pi=3.141592653589793;
        private double b;

        public NormalizeChebysev(double a, double b){
            this.a = a;
            this.b = b;
        }
        /**
         *  ? Lấy các điểm nghiệm của đa thức Chebysev làm mốc 
        **/
        public List<Point> OptimizePoint(int numPoint){
            List<Point> list = new List<Point>{};
            // double x=0;
            double t=0;
            for (int i=0; i<numPoint; i++){
                // x = Math.Cos(Pi/(2*numPoint)+i*Pi/numPoint);
                // x = (2*x-b-a)/(b-a);
                t = Math.Cos(Pi/(2*numPoint)+i*Pi/numPoint);
                t = (t*(b-a)+(a+b))/2;
                Point newPoint = new Point(t);
                list.Add(newPoint);
            }
            return list;
        }
        /**
         *  ? Random ra một số lượng điểm cách đều 
        **/
        public List<Point> RandomPoint(int numPoint){
            List<Point> result = new List<Point>{};
            double step = (b-a)/(numPoint-1);
            double x = a;
            for(int i = 1; i<=numPoint; i++){
                Point newPoint = new Point(x);
                result.Add(newPoint);
                x += step;
            }
            return result;
        }

    }
}