using System;
using System.IO;
using System.Collections.Generic;
using PolynomialFunction;
using Chebysev;

namespace Noi_suy_nguoc
{   

    class Program1
    {
        static readonly string inputFile = @"inputXY.txt";/* Mốc tùy ý */
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
            // ! <--------------------INPUT----------------------------->
            /**
             * 
             * ? 1. use random Point with Chebysev function 
             * 
            **/
            // NormalizeChebysev chebysev = new NormalizeChebysev(2,4);
            // List<Point> inputList = chebysev.RandomPoint(12);
            /**
             * 
             *  ? 2. use input from file
             * 
            **/   
            List<Point> inputList = ReadFromFile();
            
            double y = 0.92109;
            /**
             * 
             * 
            */
            // ! <---------------------------------------------------------------------------->
            ///////////////////////////////////////////////////////////////////////

            // ! <-------------------------------- METHODS -------------------------------------->
            InverseInterpolation process = new InverseInterpolation(inputList,y);

            
            // ? ---------------------------------- RUNNING -------------------------------------- */

            List<double> result = process.Result_Interpolation();

            // ! <---------------------------------------------------------------------------------->

            /**
             * 
             *  ? Print result
             * 
            **/
            int countBij_spaces = process.ThisSpaces.Count;
            List<(int left, int right)> spaces = process.ThisSpaces;
            List<(bool flag ,(int left, int right) space)> spacesInter = process.ThisSpacesInterpol;
            using(StreamWriter fileWrite = new StreamWriter(outputFile)){

                fileWrite.WriteLine($"-------------------------------------------------------------");
                fileWrite.WriteLine($"y = {y}");
                for(int i = 0; i < countBij_spaces; i++){
                    int left = spacesInter[i].space.left;
                    int right = spacesInter[i].space.right;
                    fileWrite.WriteLine($"Khoang {i+1}: ({inputList[spaces[i].left].ThisX}; {inputList[spaces[i].right].ThisX})");
                    fileWrite.Write($"Lay: [");
                    for(int j = left; j<=right;j++){
                        fileWrite.Write($"{inputList[j].ThisX} ");
                    }
                    fileWrite.WriteLine($"]\n\tx_{i+1} = {result[i]}\n\n");
                }
            }
        }
    }
}
