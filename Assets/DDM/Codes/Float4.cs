using MathNet.Numerics.LinearAlgebra;
using System;
using UnityEngine;

namespace DDM_Impl
{
    public struct Float4
    {
        public float x;
        public float y;
        public float z;
        public float w;
        public float this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return x;
                    case 1: return y;
                    case 2: return z;
                    case 3: return w;
                    default:
                        throw new IndexOutOfRangeException();
                }
            }
            set
            {
                switch (index)
                {
                    case 0: x = value; break;
                    case 1: y = value; break;
                    case 2: z = value; break;
                    case 3: w = value; break;
                    default:
                        break;
                }
            }
        }
    }
    public struct Float3
    {
        public float x;
        public float y;
        public float z;
        public float this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return x;
                    case 1: return y;
                    case 2: return z;
                    default:
                        throw new IndexOutOfRangeException();
                }
            }
            set
            {
                switch (index)
                {
                    case 0: x = value; break;
                    case 1: y = value; break;
                    case 2: z = value; break;
                    default:
                        break;
                }
            }
        }
        public Matrix<float> ToVertical()
        {
            Matrix<float> matrix = Matrix<float>.Build.Dense(3, 1);
            matrix[0, 0] = x;
            matrix[1, 0] = y;
            matrix[2, 0] = z;
            return matrix;
        }
        public Matrix<float> ToHorizontal()
        {
            Matrix<float> matrix = Matrix<float>.Build.Dense(1,3);
            matrix[0, 0] = x;
            matrix[1, 1] = y;
            matrix[0, 2] = z;
            return matrix;
        }
    }
}