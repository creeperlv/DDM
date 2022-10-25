using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

namespace DDM_Impl
{
    public class DDM_Surface : MonoBehaviour
    {
        public ComputeShader DDM_Stage_0;
        public Mesh Mesh;
        public SkinnedGeometry DDMMesh;
        public ComputeShader DDM_Stage_1;
        public SkinnedMeshRenderer SkinRenderer;
        public Image DebugImage;
        public float Lambda;
        public int iterations=1;
        // Start is called before the first frame update
        void Start()
        {
            init();
        }
        void init()
        {
            if (!SystemInfo.supportsComputeShaders)
            {
                Debug.Log("Computer Shaders is not supported.");
                Debug.Log("DDM will not work.");
                return;
            }
            Mesh = SkinRenderer.sharedMesh;

            DDMMesh = new SkinnedGeometry(Mesh.vertexCount);
            CalcuatePrecomputes();
        }
        Matrix<float> normalizedLaplacian;
        public void CalcuatePrecomputes()
        {
            CalcuateLaplace();
            CalcBMatrix(Lambda);
            BuildSkinnedGeometry();
        }
        public void BuildSkinnedGeometry()
        {
            /* 
             * Equation.1 
             */
            var bones = SkinRenderer.bones;
            Matrix<float>[] M = new Matrix<float>[bones.Length];
            var bp = Mesh.bindposes;
            for (int i = 0; i < bones.Length; i++)
            {
                M[i] = bones[i].worldToLocalMatrix.ToMatrix();

            }
            var bw = Mesh.boneWeights;
            var vert = Mesh.vertices;
            Vector3[] vs = new Vector3[vert.Length];
            Matrix<float>[] vs_mat = new Matrix<float>[vert.Length];
            Matrix<float>[] Smooth_vs_mat = new Matrix<float>[vert.Length];
            Matrix<float>[] us = new Matrix<float>[vert.Length];
            Matrix<float>[] Smooth_us = new Matrix<float>[vert.Length];
            for (int i = 0; i < bw.Length; i++)
            {
                var bw_i = bw[i];
                var vert_i = vert[i];

                Matrix<float> u = Matrix<float>.Build.DenseOfArray(new float[,] {
                    { vert_i.x },
                    { vert_i.y },
                    { vert_i.z },
                    { 1 }
                });
                us[i] = u;
                Matrix<float> v;
                {
                    v = M[bw_i.boneIndex0].Multiply(bw_i.weight0).Multiply(u);
                    v = v.Add(M[bw_i.boneIndex1].Multiply(bw_i.weight1).Multiply(u));
                    v = v.Add(M[bw_i.boneIndex2].Multiply(bw_i.weight2).Multiply(u));
                    v = v.Add(M[bw_i.boneIndex3].Multiply(bw_i.weight3).Multiply(u));
                }
                vs_mat[i] = v;
                vs[i]=new Vector3(v[0, 0], v[1, 0], v[2,0]);
            }

            for (int i = 0; i < vert.Length; i++)
            {
                Smooth_vs_mat[i] = IterativeCalcB(vs_mat[i], iterations);
                Smooth_us[i] = IterativeCalcB(us[i], iterations);
            }
        }
        Matrix<float> A;
        Matrix<float> B;
        public void CalcAMatrix(float Lambda)
        {
            A=Matrix<float>.Build.DiagonalIdentity(normalizedLaplacian.RowCount, normalizedLaplacian.ColumnCount).Subtract(normalizedLaplacian.Multiply(Lambda));
        }
        public void CalcBMatrix(float Lambda)
        {
            B = Matrix<float>.Build.DiagonalIdentity(normalizedLaplacian.RowCount, normalizedLaplacian.ColumnCount).Add(normalizedLaplacian.Multiply(Lambda));
        }
        public Matrix<float> IterativeCalcA(Matrix<float> m, int it)
        {
            if (it == 0) return m;
            it -= 1;
            return IterativeCalcA(m, it)*A;
        }
        public Matrix<float> IterativeCalcB(Matrix<float> m, int iteration)
        {
            if (iteration == 0) return m;
            iteration -= 1;
            return B.Solve(IterativeCalcB(m, iteration));
        }
        public void CalcuateLaplace()
        {

            AdjacencyMatrix _AdjacencyMatrix = AdjacencyMatrix.FromMesh(Mesh);
            var deg = DegMatrix.FromMesh(Mesh);
            Matrix<float> DegMat = deg.Matrix;
            DegMat = DegMat.Multiply(3);
            Matrix<float> LaplacianMatrix = DegMat.Subtract(_AdjacencyMatrix);
            var DL = LaplacianMatrix.Diagonal();
            var DL_m = Matrix<float>.Build.Diagonal(DL.ToArray());
            var DL_m_I = DL_m.Inverse();
            normalizedLaplacian = LaplacianMatrix.Multiply(DL_m_I);

            var t = DebugUtilities.FromMatrix(normalizedLaplacian.ToArray(), _AdjacencyMatrix.n, _AdjacencyMatrix.n, 3, 0.2f);
            var s = Sprite.Create(t, new Rect(0, 0, _AdjacencyMatrix.n, _AdjacencyMatrix.n), new Vector2(.5f, .5f));
            DebugImage.sprite = s;
        }
        // Update is called once per frame
        void Update()
        {

        }
        void UpdateCS()
        {

        }
    }
    public struct SkinnedGeometry
    {
        public SkinnedGeometry(int Count)
        {
            Vertices = new Vector3[Count];
            Normals = new Vector3[Count];
            DeltaNor = new Vector3[Count];
            DeltaVer = new Vector3[Count];
        }
        public Vector3[] Vertices;
        public Vector3[] Normals;
        public Vector3[] DeltaVer;
        public Vector3[] DeltaNor;
    }
}