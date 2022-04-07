using System;
using System.IO;
using System.Collections.Generic;
using PolynomialFunction;
using Chebysev;

namespace Noi_suy_lagrange
{
    class Program
    {
        static readonly string inputFile = @"input.txt";
        static readonly string outputFile = @"output.txt";
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
            /* use Chebysev */
            // NormalizeChebysev norm = new NormalizeChebysev(2,3); /* input (a,b) */
            // List<Point> inputNormalize = norm.OptimizePoint(4);  /* input number of points */
            // Lagrange process = new Lagrange(inputNormalize);

            /* use file input */
            Lagrange process = new Lagrange(ReadFromFile());    /* input from file - not use Chebysev*/
            
            
            Polynomial result = process.Interpolation();

            /* Print result with errors*/
            using(StreamWriter fileWrite = new StreamWriter(outputFile)){
                
                fileWrite.WriteLine($"P_n(x) = \n {result.ToString()}\n\n");
                /* Error: Sai số */
                for(int k = 0; k<=process.ThisDeg; k++){
                    double x_k = process.ThisInputXY[k].ThisX;
                    double y_k = process.ThisInputXY[k].ThisY;
                    double Px_k = result.f_At(x_k);
                    fileWrite.Write($"*At x_{k} = {x_k}       ");
                    fileWrite.WriteLine($"Px_{k} - y_{k} = {Px_k - y_k}\n");
                }

                // * calculate P_n(x) */ 
                // double x = 0;
                // double px = 0;
                // char key = 'y';
                // while(key=='y'){
                //     Console.WriteLine($"Calculate P_n(x), x = ?");
                //     x = double.Parse(Console.ReadLine());
                //     px = result.f_At(x);
                //     Console.WriteLine($"f({x}) ~~ {px}");
                //     fileWrite.WriteLine($"f({x}) ~~ {px}");
                //     Console.WriteLine("Continue? press 'y' || press 'n' if not");
                //     key = Char.Parse(Console.ReadLine());
                // }

            }
        }
    }
}
