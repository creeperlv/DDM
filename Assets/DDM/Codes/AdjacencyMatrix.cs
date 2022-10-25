using MathNet.Numerics.LinearAlgebra;
using UnityEngine;

namespace DDM_Impl
{
    public class DegMatrix
    {
        public int Size;
        public Matrix<float> Matrix;
        public static DegMatrix FromMesh(Mesh mesh)
        {
            DegMatrix degMatrix = new DegMatrix();
            degMatrix.Size = mesh.vertices.Length;
            var triangles = mesh.triangles;
            degMatrix.Matrix = Matrix<float>.Build.Sparse(degMatrix.Size, degMatrix.Size);
            for (int i = 0; i < triangles.Length - 2; i += 3)
            {
                var i_0 = triangles[i];
                var i_1 = triangles[i + 1];
                var i_2 = triangles[i + 2];
                {
                    degMatrix.Matrix[i_0, i_0] += 2;
                    degMatrix.Matrix[i_1, i_1] += 2;
                    degMatrix.Matrix[i_2, i_2] += 2;
                }
            }
            return degMatrix;
        }
    }
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
        public static Matrix<float> ToMatrix (this Matrix4x4 m)
        {
            return Matrix<float>.Build.DenseOfArray(new float[,] { { m[0, 0], m[0, 1], m[0, 2], m[0, 3] },
                { m[1, 0], m[1, 1], m[1, 2], m[1, 3] } ,
                { m[2, 0], m[2, 1], m[2, 2], m[2, 3] },
                { m[3, 0], m[3, 1], m[3, 2], m[3, 3] }});
            
        }
    }
    [System.Serializable]
    public class AdjacencyMatrix
    {
        public int n;
        public Matrix<float> data;
        public static AdjacencyMatrix FromMesh(Mesh mesh)
        {
            var vertices = mesh.vertices;
            var triangles = mesh.triangles;
            Matrix<float> data = Matrix<float>.Build.Sparse(vertices.Length, vertices.Length);
            for (int i = 0; i < triangles.Length - 2; i += 3)
            {
                var i_0 = triangles[i];
                var i_1 = triangles[i + 1];
                var i_2 = triangles[i + 2];
                {
                    data[i_0, i_1] = 1;
                    data[i_0, i_2] = 1;
                    data[i_1, i_2] = 1;
                }
            }
            AdjacencyMatrix adjacencyMatrix = new AdjacencyMatrix { n = vertices.Length, data = data };
            return adjacencyMatrix;
        }
        public static implicit operator float[,](AdjacencyMatrix matrix)
        {
            return matrix.data.ToArray();
        }
        public static implicit operator Matrix<float>(AdjacencyMatrix matrix)
        {
            return matrix.data;
        }
    }
}