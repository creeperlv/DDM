using MathNet.Numerics.LinearAlgebra;
using System;
using UnityEngine;

namespace DDM_Impl
{
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
                    float3X3[x, y] = mat[x, y];
                }
            }
            return float3X3;
        }

        public override string ToString()
        {
            return $"{x0y0},{x0y1},{x0y2}\n{x1y0},{x1y1},{x1y2}\n{x2y0},{x2y1},{x2y2}";
        }
    }
    public struct RowFirstFloat4x4
    {
        public float m00;
        public float m01;
        public float m02;
        public float m03;

        public float m10;
        public float m11;
        public float m12;
        public float m13;

        public float m20;
        public float m21;
        public float m22;
        public float m23;

        public float m30;
        public float m31;
        public float m32;
        public float m33;

        public RowFirstFloat4x4(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30, float m31, float m32,float m33) : this()
        {
            this.m00 = m00;
            this.m01 = m01;
            this.m02 = m02;
            this.m03 = m03;
            this.m10 = m10;
            this.m11 = m11;
            this.m12 = m12;
            this.m13 = m13;
            this.m20 = m20;
            this.m21 = m21;
            this.m22 = m22;
            this.m23 = m23;
            this.m30 = m30;
            this.m31 = m31;
            this.m32 = m32;
            this.m33 = m33;
        }

        public static RowFirstFloat4x4 FromMatrix(Matrix4x4 mat) {
            return new RowFirstFloat4x4(mat.m00, mat.m01, mat.m02, mat.m03,
                mat.m10, mat.m11, mat.m12, mat.m13,
                mat.m20, mat.m21, mat.m22, mat.m23,
                mat.m30, mat.m31, mat.m32, mat.m33
                );
        }
    }
}