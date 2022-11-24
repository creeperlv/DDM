using MathNet.Numerics.LinearAlgebra;
using System.Runtime.CompilerServices;
using UnityEngine;

namespace DDM_Impl
{
    public static class MatrixUtils
    {
        public static float[,] Diagonal(int[,] matrix, int n)
        {
            float[,] data = new float[n, n];
            for (int i = 0; i < n; i++)
            {
                data[i, i] = matrix[i, i];
            }
            return data;
        }
        [MethodImpl (MethodImplOptions.AggressiveInlining)] 
        public static Vector3 ToVector3(Matrix<float> mat)
        {
            return new Vector3(mat[0, 0], mat[1, 0], mat[2,0]);
        }
        public static Matrix4x4 RootSpaceMatrix(this Transform transform,Transform Root) {
            return Matrix4x4.TRS(transform.position - Root.position,Quaternion.FromToRotation(transform.forward,Root.forward)/* Quaternion.Euler(transform.eulerAngles-Root.eulerAngles)*/, 
                //transform.localScale
                Vector3.one
                );
        }
        public static Matrix4x4 LocalToMatrix(this Transform t)
        {
            return Matrix4x4.TRS(t.localPosition, t.localRotation, t.localScale);
        }
        public static Matrix4x4 GlobalToMatrix(this Transform t)
        {
            return Matrix4x4.TRS(t.position, t.rotation, t.lossyScale);
        }
        public static Matrix<float> ToMatrix (this Matrix4x4 m)
        {
            Matrix<float> matrix = Matrix<float>.Build.Dense(4, 4);
            for (int x = 0; x < 4; ++x)
            {
                for (int y = 0; y < 4; ++y)
                {
                    matrix[x, y] = m[x,y];
                }
            }
            return matrix;
            //return Matrix<float>.Build.DenseOfArray(new float[,] { { m[0, 0], m[0, 1], m[0, 2], m[0, 3] },
            //    { m[1, 0], m[1, 1], m[1, 2], m[1, 3] } ,
            //    { m[2, 0], m[2, 1], m[2, 2], m[2, 3] },
            //    { m[3, 0], m[3, 1], m[3, 2], m[3, 3] }});
            
        }
    }
}