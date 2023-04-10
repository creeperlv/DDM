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
    public struct Float3x3
    {
        public float x0y0;
        public float x0y1;
        public float x0y2;
        public float x1y0;
        public float x1y1;
        public float x1y2;
        public float x2y0;
        public float x2y1;
        public float x2y2;
        public float this[int x, int y]
        {
            get => this[x * 3 + y];
            set => this[x * 3 + y] = value;
        }
        public float this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0:
                        return x0y0;
                    case 1: return x0y1;
                    case 2: return x0y2;
                    case 3: return x1y0;
                    case 4: return x1y1;
                    case 5: return x1y2;
                    case 6: return x2y0;
                    case 7: return x2y1;
                    case 8: return x2y2;
                    default:
                        break;
                }
                throw new IndexOutOfRangeException("Invalid matrix index!");
            }
            set
            {
                switch (index)
                {
                    case 0: x0y0 = value; return;
                    case 1: x0y1 = value; return;
                    case 2: x0y2 = value; return;
                    case 3: x1y0 = value; return;
                    case 4: x1y1 = value; return;
                    case 5: x1y2 = value; return;
                    case 6: x2y0 = value; return;
                    case 7: x2y1 = value; return;
                    case 8: x2y2 = value; return;
                    default:
                        break;
                }
                throw new IndexOutOfRangeException("Invalid matrix index!");
            }
        }
        public Matrix<float> ToMatrix()
        {
            var mat = Matrix<float>.Build.Dense(3, 3);
            for (int x = 0; x < 3; x++)
            {
                for (int y = 0; y < 3; y++)
                {
                    mat[x, y] = this[x, y];
                }
            }
            return mat;
        }
        public static Float3x3 FromMatrix(Matrix<float> mat)
        {
            Float3x3 float3X3 = new Float3x3();

            for (int x = 0; x < 3; x++)
            {
                for (int y = 0; y < 3; y++)
                {
                    mat[x, y] = float3X3[x, y];
                }
            }
            return float3X3;
        }

        public override string ToString()
        {
            return $"{x0y0},{x0y1},{x0y2}\n{x1y0},{x1y1},{x1y2}\n{x2y0},{x2y1},{x2y2}";
        }
    }
}