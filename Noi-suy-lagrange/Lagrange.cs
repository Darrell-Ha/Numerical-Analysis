using System;
using System.Collections.Generic;
using PolynomialFunction;

namespace Noi_suy_lagrange
{
    public class Lagrange{
        private List<Point> inputXY;
        private int deg; 
        private Polynomial P_n;

        public Lagrange (List<Point> inputXY){
            this.inputXY = inputXY;
            this.deg = inputXY.Count-1;
            this.P_n = new Polynomial();
            // ! initialized P_0(x) = 0;
        }

        public Polynomial Interpolation(){

            // ! Lagrange use formula 2

            // * Step 1: Calculate Omega
            // Polynomial Omega = new Polynomial(new List<double>{1,-ThisInputXY[0].ThisX});
            Polynomial Omega = new Polynomial(new List<double>{1});
            for(int i = 0; i <= ThisDeg; i++){
                Omega *= (1,-ThisInputXY[i].ThisX);
            }
            /**
             * * Step 2: Calculate coefficient of Polynomial
             *      * Q(x) = Omega/(x-x_k)  
             *      * D_k = Q(x_k)
             *      * y_k*L_k = (y_k/D_k)*Q(x)
            **/
            for(int k = 0; k <= ThisDeg; k++){
                double x_k = ThisInputXY[k].ThisX;
                double y_k = ThisInputXY[k].ThisY;
                
                Polynomial Qx = (Omega/(1,-x_k)).Qx;
                double D_k = Qx.f_At(x_k);
                Polynomial P_x_k = (y_k/D_k)*Qx;

                // * update P_n
                ThisP_n += P_x_k;  
            }

            return ThisP_n;

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
        public List<double> GetInputX(){
            List<double> inputX = new List<double>{};
            for(int i=0; i<=ThisDeg; i++){
                
                inputX.Add(ThisInputXY[i].ThisX);
            }
            return inputX;
        }
    }
}