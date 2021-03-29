using System;
using System.Collections.Generic;

namespace CMLA_1_Matrix_Gauss
{

    public class Matrix
    {
        private double[,] data;

        public int n { get; }
        public int m { get; }

        public Matrix(int m, int n)
        {
            this.m = m;
            this.n = n;
            this.data = new double[m, n];
        }

        public double this[int i, int j]
        {
            get
            {
                return this.data[i, j];
            }
            set
            {
                this.data[i, j] = value;
            }
        }


        public void ProcessFunctionOutput(Action<int, int> func)
        {
            for (var i = 0; i < this.m; i++)
            {
                for (var j = 0; j < this.n; j++)
                {
                    func(i, j);
                }
                Console.WriteLine("");
            }
        }

        public void ProcessFunctionOverData(Action<int, int> func)
        {
            for (var i = 0; i < this.m; i++)
            {
                for (var j = 0; j < this.n; j++)
                {
                    func(i, j);
                }

            }
        }

        public void SwapRow(int k, int l)
        {
            for (var i = 0; i < this.n; i++)
            {
                double t = data[k, i];
                data[k, i] = data[l, i];
                data[l, i] = t;
            }
        }

        public int MaxElem(int index)
        {
            int FinalIndex = index;
            double MaxElem = data[index, index];
            for (var i = index + 1; i < this.n; i++)
            {
                if (Math.Abs(data[i, index]) > Math.Abs(MaxElem))
                {
                    MaxElem = data[i, index];
                    FinalIndex = i;
                }

            }
            return FinalIndex;
        }

        public void SingleMatrix()
        {
            ProcessFunctionOverData((i, j) =>
            {

                if (i == j)
                    data[i, j] = 1;
                else
                    data[i, j] = 0;

            });
        }

        public void InputElem()
        {
            ProcessFunctionOverData((i, j) =>
            {
                Console.Write($"data[{i + 1}][{j + 1}] = ");
                data[i, j] = Double.Parse(Console.ReadLine());
            });
        }

        public void OutputMatrix()
        {
            ProcessFunctionOutput((i, j) => {
                Console.Write(data[i, j] + " ");
            });
        }

        public void Round()
        {
            ProcessFunctionOutput((i, j) => {
                data[i, j] = Math.Round(data[i, j], 3);
                Console.Write(data[i, j] + " ");
            });
        }

        public void OutPutMatrixToInt()
        {
            ProcessFunctionOutput((i, j) => {
                Console.Write(Math.Round(Math.Abs(data[i, j])) + " ");
            });
  
        }

        public void CopyMatrix(Matrix matrix)
        {
            ProcessFunctionOverData((i, j) =>
            {
                data[i, j] = matrix[i, j];
            });
        }

        public static Matrix operator *(Matrix matrix, Matrix matrix2)
        {
            if (matrix.n != matrix2.m)
            {
                throw new ArgumentException("matrixes can not be multiplied");
            }
            var result = new Matrix(matrix.m, matrix2.n);
            result.ProcessFunctionOverData((i, j) => {
                for (var k = 0; k < matrix.n; k++)
                {
                    result[i, j] += matrix[i, k] * matrix2[k, j];
                }
            });
            return result;
        }

       

        private Matrix CreateMatrixWithoutRow(int row)
        {
            if (row < 0 || row >= this.m)
            {
                throw new ArgumentException("invalid row index");
            }
            var result = new Matrix(this.m - 1, this.n);
            result.ProcessFunctionOverData((i, j) => result[i, j] = i < row ? this[i, j] : this[i + 1, j]);
            return result;
        }

        private Matrix CreateMatrixWithoutColumn(int column)
        {
            if (column < 0 || column >= this.n)
            {
                throw new ArgumentException("invalid column index");
            }
            var result = new Matrix(this.n, this.m - 1);
            result.ProcessFunctionOverData((i, j) => result[i, j] = j < column ? this[i, j] : this[i, j + 1]);
            return result;
        }


        public double CalculateDeterminant()
        {
            
            if (this.n == 2)
            {
                return this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0];
            }
            double result = 0;
            for (var j = 0; j < this.n; j++)
            {
                result += (j % 2 == 1 ? 1 : -1) * this[1, j] * this.CreateMatrixWithoutColumn(j).CreateMatrixWithoutRow(1).CalculateDeterminant();
            }
            return result;
        }

        public bool IsSymmetric()
        {
            bool IsTrue = true;
            ProcessFunctionOverData((i, j) => {
                if (!(data[i, j] == data[j, i]))
                    IsTrue = false;
            });
            return IsTrue;
        }

        public Matrix Transpose()
        {
            Matrix result = new Matrix(n,n);
            ProcessFunctionOverData((i, j) => {
                result[j, i] = data[i, j];
            });
            return result; 
        }

    }


    class Program
    {
        static Matrix GaussMethod(Matrix A, Matrix E, int n)
        {
            for (var k = 0; k < n - 1; k++)
            {
                int MaxIndex = A.MaxElem(k);
                A.SwapRow(k, MaxIndex);
                E.SwapRow(k, MaxIndex);
    
                for (var i = k + 1; i < n; i++)
                {
                    double m = (-1) * (A[i, k] / A[k, k]);
                    for (var j = 0; j < n; j++)
                    {
                        A[i, j] = A[i, j] + m * A[k, j];
                        E[i, j] = E[i, j] + m * E[k, j];
                    }
                }
                //Console.WriteLine();
                //A.OutputMatrix(); 
            }
    
            Matrix X = new Matrix(n, n);
            for (var i = n - 1; i >= 0; i--)
            {
                for (var j = 0; j < n; j++)
                {
                    double s = 0;
                    for (var k = i + 1; k < n; k++)
                    {
                        s += A[i, k] * X[k, j];
                    }
                    X[i, j] = (E[i, j] - s) / A[i, i];
    
                }
            }
            return X;
        }



        static Matrix UU(Matrix A, Matrix U, int n, Matrix b)
        {
            double[] y = new double[n];
            Matrix x = new Matrix(n, 1);
            for (int i = 0; i < n; i++)
            {
                double SumK = 0;
                for (int k = 0; k <= i-1; k++)
                {
                   SumK +=(U[k, i] * U[k, i]);
                }
                if(SumK > A[i,i])
                    throw new Exception("Division by 0");

                U[i, i] = Math.Sqrt(A[i, i] - SumK);

                if (U[i, i] == 0)
                    throw new Exception("Diagonal element is 0");
                for (int j = i+1; j < n; j++)
                {
                    SumK = 0;
                    for (int k = 0; k <= i-1; k++)
                    {
                       SumK+=(U[k, i] * U[k, j]);
                    }
                    U[i, j] = (A[i, j] - SumK) / U[i, i];
                }

            }
            for (int i = 0; i < n; i++)
            {
                double SumK = 0;
                for (int k = 0; k <= i-1; k++)
                {
                    SumK += (U[k, i] * y[k]);
                }
                y[i] = (1 / U[i, i]) * (b[i,0] - SumK);

            }
            for (int i = n-1; i >= 0; i--)
            {
                double SumK = 0;
                for (int k = i+1; k < n; k++)
                {
                    SumK += (U[i, k] * x[k,0]);
                }
                x[i,0] = (1 / U[i, i]) * (y[i] - SumK);
            }

            return x;
            
        }

        static void Main(string[] args)
        {
            //Console.WriteLine("Input n: ");
            //int n = Int32.Parse(Console.ReadLine());
            //Matrix A = new Matrix(n, n);
            //Matrix AOriginal = new Matrix(n, n);
            //Matrix SM = new Matrix(n, n);
            //Matrix X = new Matrix(n, n);
            //A.InputElem();
            //Console.WriteLine("\nA^ :\n");
            //A.OutputMatrix();
            //if(A.CalculateDeterminant() == 0){
            //    Console.WriteLine("\ndet(A) = 0");
            //}
            //else
            //{
            //    AOriginal.CopyMatrix(A);
            //    SM.SingleMatrix();
            //    X = GaussMethod(A, SM, n);
            //    Console.WriteLine("\nA^(-1) :\n");
            //    X.Round();

            //    Console.WriteLine("\nA * A^(-1) = E\n");
            //    Matrix tmp = new Matrix(n, n);
            //    tmp = AOriginal * X;
            //    tmp.OutPutMatrixToInt();


            //}

            Console.WriteLine("Input n: ");
            int n = Int32.Parse(Console.ReadLine());
            Matrix A = new Matrix(n, n);
            Console.WriteLine("\nEnter A Elem:");
            A.InputElem();
            if (A.IsSymmetric())
            {
                Matrix U = new Matrix(n, n);
                Matrix b = new Matrix(n, 1);
                Console.WriteLine("\nEnter b Elem:");
                b.InputElem();
                Matrix x = new Matrix(n, 1);

                x = UU(A, U, 3, b);

                Console.WriteLine("\nx :");
                x.OutputMatrix();

                Console.WriteLine("\nU :");
                U.OutputMatrix();

                double detA = 0;
                for (int i = 0; i < n; i++)
                {
                    detA *= U[i, i];
                }
                Console.WriteLine($"\ndet(A) = det(U)^2 = {Math.Pow(detA,2.0F)}");

                Console.WriteLine("\nU^t * U :");
                Matrix tmpU = U.Transpose() * U;
                tmpU.OutputMatrix();

                Console.WriteLine("\nA :");
                A.OutputMatrix();

                Console.WriteLine("\nA * x :");
                Matrix tmp = A * x;
                tmp.OutputMatrix();

                Console.WriteLine("\nb :");
                b.OutputMatrix();

                Console.ReadKey();
            }
            else
            {
                throw new Exception("Matrix not symmetric");
            }


            Console.ReadKey();


        }
    }
    
}
