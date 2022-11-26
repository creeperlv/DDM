using MathNet.Numerics.LinearAlgebra;
using System.Runtime.CompilerServices;
using UnityEngine;
using UnityEngine.Profiling;

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
        public static int[,] BuildAdjacencyMatrix(Vector3[] v, int[] t, int MAX_WIDTH, float MIN_DIST = 0.00001f)
        {
            var adj = new int[v.Length, MAX_WIDTH];
            for (int i = 0; i < adj.GetLength(0); ++i)
                for (int j = 0; j < adj.GetLength(1); ++j)
                    adj[i, j] = -1;

            if (MIN_DIST == 0.0f)
            {
                for (int tri = 0; tri < t.Length; tri = tri + 3)
                {
                    AddEdgeToAdj_Directly(ref adj, t[tri], t[tri + 1]);
                    AddEdgeToAdj_Directly(ref adj, t[tri], t[tri + 2]);
                    AddEdgeToAdj_Directly(ref adj, t[tri + 1], t[tri + 2]);
                }
            }
            else
            {
                int[] mapToUnique = MapVertToPos(v, MIN_DIST);

                for (int tri = 0; tri < t.Length; tri = tri + 3)
                {
                    AddEdgeToAdj(ref adj, mapToUnique, t[tri], t[tri + 1]);
                    AddEdgeToAdj(ref adj, mapToUnique, t[tri], t[tri + 2]);
                    AddEdgeToAdj(ref adj, mapToUnique, t[tri + 1], t[tri + 2]);
                }

                BroadcastAdjacencyFromUniqueToAllVertices(ref adj, mapToUnique);
            }

            return adj;
        }

        private static void BroadcastAdjacencyFromUniqueToAllVertices(ref int[,] adjacencyMatrix, int[] mapToUnique)
        {

            var maxNeighbors = adjacencyMatrix.GetLength(1);

            for (int i = 0; i < mapToUnique.Length; ++i)
            {
                var u = mapToUnique[i];
                if (u == i)
                    continue;

                for (int j = 0; j < maxNeighbors && adjacencyMatrix[u, j] != -1; ++j)
                    adjacencyMatrix[i, j] = adjacencyMatrix[u, j];
            }
        }
        public static int[] MapVertToPos(Vector3[] v, float minSqrDistance = 0.00001f)
        {

            var mapToUnique = new int[v.Length];
            for (int i = 0; i < mapToUnique.Length; ++i)
                mapToUnique[i] = -1;

            for (int i = 0; i < v.Length; i++)
                for (int j = i; j < v.Length; j++)
                    if (mapToUnique[j] == -1)
                    {
                        var u = mapToUnique[i];
                        if (u == -1)
                            u = i;

                        var dx = v[u].x - v[j].x;
                        var dy = v[u].y - v[j].y;
                        var dz = v[u].z - v[j].z;
                        if (dx * dx + dy * dy + dz * dz <= minSqrDistance)
                        {
                            if (mapToUnique[i] == -1)
                                mapToUnique[i] = u; 
                            mapToUnique[j] = u;
                        }
                    }


            return mapToUnique;
        }

        private static void AddVert_Adj(ref int[,] adjacencyMatrix, int from, int to)
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

        private static void AddEdgeToAdj_Directly(ref int[,] adjacencyMatrix, int v0, int v1)
        {
            AddVert_Adj(ref adjacencyMatrix, v0, v1);
            AddVert_Adj(ref adjacencyMatrix, v1, v0);
        }

        private static void AddEdgeToAdj(ref int[,] adjacencyMatrix, int[] mapToUnique, int v0, int v1)
        {
            var u0 = mapToUnique[v0];
            var u1 = mapToUnique[v1];

            AddEdgeToAdj_Directly(ref adjacencyMatrix, u0, u1);
        }
        public static Matrix<float> BuildLaplacianMatrixFromAdjacentMatrix(
            int VertexLength, int[,] AdjacencyMatrix, bool normalize = true)
        {
            Matrix<float> lapl = Matrix<float>.Build.Sparse(VertexLength, VertexLength);
            int maxNeighbors = AdjacencyMatrix.GetLength(1);

            for (int vi = 0; vi < VertexLength; vi++)
            {
                int viDeg = 0;
                for (int j = 0; j < maxNeighbors; j++)
                {
                    int vj = AdjacencyMatrix[vi, j];
                    if (vj < 0)
                    {
                        break;
                    }
                    ++viDeg;
                    lapl[vi, vj] = -1;
                }

                if (!normalize)
                {
                    lapl[vi, vi] = viDeg;
                }
                else
                {
                    for (int j = 0; j < maxNeighbors; j++)
                    {
                        int vj = AdjacencyMatrix[vi, j];
                        if (vj < 0)
                        {
                            break;
                        }
                        lapl[vi, vj]= lapl[vi, vj] / viDeg;
                    }
                    lapl[vi, vi]= 1.0f;
                }

            }

            return lapl;
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3 ToVector3(Matrix<float> mat)
        {
            return new Vector3(mat[0, 0], mat[1, 0], mat[2, 0]);
        }
        public static Matrix4x4 RootSpaceMatrix(this Transform transform, Transform Root)
        {
            return Matrix4x4.TRS(transform.position - Root.position, Quaternion.FromToRotation(transform.forward, Root.forward)/* Quaternion.Euler(transform.eulerAngles-Root.eulerAngles)*/,
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
        public static void ToMatrix(this Matrix4x4 m, Matrix<float> target)
        {
            for (int x = 0; x < 4; x++)
            {
                for (int y = 0; y < 4; y++)
                {
                    target[x, y] = m[x, y];
                }
            }
        }
        public static Matrix<float> ToMatrix(this Matrix4x4 m)
        {
            Matrix<float> matrix = Matrix<float>.Build.Dense(4, 4);
            for (int x = 0; x < 4; ++x)
            {
                for (int y = 0; y < 4; ++y)
                {
                    matrix[x, y] = m[x, y];
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