using System;
using System.IO;
using System.Collections.Generic;
using PolynomialFunction;
using Chebysev;

namespace Noi_suy_newton
{   

    class Program
    {
        static readonly string inputFile = @"inputXY.txt";/* Mốc tùy ý */
        static readonly string outputFile = @"output.txt";/* Mốc cách đều, tiến lùi tùy */
        public static List<Point> ReadFromFile(){
            List<Point> input = new List<Point>{};
            
            if(File.Exists(inputFile)){
                string line="";
                using(StreamReader file = new StreamReader(inputFile)){
                    while((line=file.ReadLine())!=null){ //
                        string [] xy = line.Split(" ");
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
            /* use random Point with Chebysev function */
            // NormalizeChebysev chebysev = new NormalizeChebysev(2,4);
            // List<Point> inputList = chebysev.OptimizePoint(10);

            /* use input from file */   
            List<Point> inputList = ReadFromFile();


            /* Newton: 
                mốc tùy ý = 0      Forward = 1    backward = -1 
            */
            Newton process = new Newton(inputList,0);
            
            /* ---------------------------------- Running -------------------------------------- */

            Polynomial result = process.Interpolation();
            /* Print Diff table */
            // foreach(List<double> list in process.ThisDiff){
            //     foreach (double num in list){
            //         Console.Write($"{num}   ");
            //     }
            //     Console.WriteLine();
            // }
            /* ---------------------------------------------------------------------------------- */
            /* Print result with errors*/
            using(StreamWriter fileWrite = new StreamWriter(outputFile)){
                
                fileWrite.WriteLine($"P_n(x) = \n {result.ToString()}\n\n");
                /* Error: Sai số */
                for(int k = 0; k<=process.ThisDeg; k++){
                    double x_k = process.ThisInputXY[k].ThisX;
                    double y_k = process.ThisInputXY[k].ThisY;
                    /* Applied Newton forward/backward */
                    double t_k = process.Exchange_x_To_t(x_k);
                    double Px_k = result.f_At(t_k);
                    /*  */
                    fileWrite.WriteLine($"*At x_{t_k} = {x_k}");
                    fileWrite.WriteLine($"\tP(x_{t_k}) - y_{t_k} = {Px_k - y_k}\n");
                }


            }
            /** 
             * ? calculate P_n(x)
            **/ 
            // double x = 0;
            // double px = 0;
            // char key = 'y';
            // while(key=='y'){
            //     Console.WriteLine($"Calculate P_n(x), x = ?");
            //     x = double.Parse(Console.ReadLine());
            //     x = process.Exchange_x_To_t(x);
            //     px = result.f_At(x);
            //     Console.WriteLine($"f({x}) ~~ {px}");
            //     Console.WriteLine("Continue? press 'y' || press 'n' if not");
            //     key = Char.Parse(Console.ReadLine());
            // }
        }
    }
}
