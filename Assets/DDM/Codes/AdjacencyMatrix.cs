using MathNet.Numerics.LinearAlgebra;
using UnityEngine;
using UnityEngine.Profiling;

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
                    data[i_0, i_1] =1;
                    data[i_0, i_2] =1;
                    data[i_1, i_2] =1;
                }
            }
            AdjacencyMatrix adjacencyMatrix = new AdjacencyMatrix { n = vertices.Length, data = data };
            return adjacencyMatrix;
        }
        public static int[,] BuildAdjacencyMatrix(Vector3[] v, int[] t, int maxNeighbors)
        {
            var adj = new int[v.Length, maxNeighbors];
            for (int i = 0; i < adj.GetLength(0); ++i)
                for (int j = 0; j < adj.GetLength(1); ++j)
                    adj[i, j] = -1;

                for (int tri = 0; tri < t.Length; tri = tri + 3)
                {
                    AddEdgeToAdjacencyMatrixDirect(ref adj, t[tri], t[tri + 1]);
                    AddEdgeToAdjacencyMatrixDirect(ref adj, t[tri], t[tri + 2]);
                    AddEdgeToAdjacencyMatrixDirect(ref adj, t[tri + 1], t[tri + 2]);
                }
            
         

            return adj;
        }
        private static void AddEdgeToAdjacencyMatrixDirect(ref int[,] adjacencyMatrix, int v0, int v1)
        {
            AddVertexToAdjacencyMatrix(ref adjacencyMatrix, v0, v1);
            AddVertexToAdjacencyMatrix(ref adjacencyMatrix, v1, v0);
        }

        private static void AddVertexToAdjacencyMatrix(ref int[,] adjacencyMatrix, int from, int to)
        {
            var maxNeighbors = adjacencyMatrix.GetLength(1);
            for (int i = 0; i < maxNeighbors; i++)
            {
                if (adjacencyMatrix[from, i] == to)
                    break;

                if (adjacencyMatrix[from, i] == -1)
                {
                    adjacencyMatrix[from, i] = to;
                    break;
                }
            }
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