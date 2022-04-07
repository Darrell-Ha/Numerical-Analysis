using System;
using System.IO;
using System.Collections.Generic;
using PolynomialFunction;

namespace Tinh_tich_phan
{
    class Program
    {
        static readonly string inputFile = @"input.txt";/* Mốc tùy ý */
        static readonly string outputFile = @"output.txt";/* Mốc cách đều, tiến lùi tùy */
        public static List<Point> ReadFromFile(){
            List<Point> input = new List<Point>{};
            
            if(File.Exists(inputFile)){
                string line="";
                using(StreamReader file = new StreamReader(inputFile)){
                    while((line=file.ReadLine())!=null){ //
                        string [] xy = line.Split("\t");
                        double x = Double.Parse(xy[0]);
                        double y = Double.Parse(xy[1]);
                        Point newPoint = new Point(x,y);
                        input.Add(newPoint);
                    }
                }
            }
            return input;
        }
        static void Main(string[] args)
        {
            List<Point> inputList =  ReadFromFile();

            Function fx = new Function();
            // Trapezoidal process = new Trapezoidal(fx,0,1,0.5*1e-6);
            // Simpson process = new Simpson(fx,0,1,0.5*1e-6);
            Simpson process = new Simpson(inputList,0.05);
            // NewtonCotez process = new NewtonCotez(fx,0,1,10);
            
            /**
             *  !! <------------------------- OUTPUT ----------------------------------> 
            **/
            double result = process.Calc_Integral();
            Console.WriteLine(result); 
            // Console.WriteLine(Math.Round(result,-(int)Math.Log10(process.PRINT_TO_SCREEN),MidpointRounding.ToEven));
            // Console.WriteLine($"n = {process.ThisN}");
            // Console.WriteLine($"Max Mn = {process.ThisMn}");
            // Console.WriteLine($"error = {process.ThisEpsilon}"); 

            /**
             * ?? In bảng hệ số Newton Cotez
            **/
            // process.PrintCotezCoeff();

        }
    }
}
