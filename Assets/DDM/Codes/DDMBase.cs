using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
namespace DDM_Impl
{
    public class DDMBase : MonoBehaviour
    {
        public SkinnedMeshRenderer TargetSMR;
        public float Lambda;
        public int iterations = 1;
        Mesh Mesh;
        Transform[] Bones;
        void Start()
        {
            Mesh= TargetSMR.sharedMesh;
            Bones = TargetSMR.bones;
            Precompute();
        }

        Matrix<float> B;
        void Precompute()
        {
            CalcuateNormalizedLaplace();
            CalcBMatrix(Lambda);
        }
        // Update is called once per frame
        void Update()
        {

        }
        void OnFrame()
        {

        }
        Matrix<float> u_;
        Matrix<float> v_;
        Matrix<float> normalizedLaplacian;
        Matrix<float> _Smooth_vs_mat;
        Matrix<float> _Smooth_us;

        public void CalcBMatrix(float Lambda)
        {
            B = Matrix<float>.Build.DiagonalIdentity(normalizedLaplacian.RowCount, normalizedLaplacian.ColumnCount).Add(normalizedLaplacian.Multiply(Lambda));
        }
        public void BuildSkinnedGeometry()
        {
            /* 
             * Equation.1 
             */
            var bones = TargetSMR.bones;
            Matrix<float>[] M = new Matrix<float>[bones.Length];
            var bp = Mesh.bindposes;
            for (int i = 0; i < bones.Length; i++)
            {
                M[i] = bones[i].worldToLocalMatrix.ToMatrix();

            }
            var bw = Mesh.boneWeights;
            var vert = Mesh.vertices;
            u_ = Matrix<float>.Build.Sparse(vert.Length, 4);
            v_ = Matrix<float>.Build.Sparse(vert.Length, 4);
            Vector3[] vs = new Vector3[vert.Length];
            Matrix<float>[] vs_mat = new Matrix<float>[vert.Length];
            _Smooth_vs_mat = Matrix<float>.Build.Sparse(vert.Length, 4);
            Matrix<float>[] us = new Matrix<float>[vert.Length];
            _Smooth_us = Matrix<float>.Build.Sparse(vert.Length, 4);
            for (int i = 0; i < bw.Length; i++)
            {
                var bw_i = bw[i];
                var vert_i = vert[i];

                Matrix<float> u = Matrix<float>.Build.DenseOfColumnArrays(new float[] { vert_i.x, vert_i.y, vert_i.z, 1 });
                us[i] = u;
                Matrix<float> v;
                {
                    v = M[bw_i.boneIndex0].Multiply(bw_i.weight0).Multiply(u);
                    v = v.Add(M[bw_i.boneIndex1].Multiply(bw_i.weight1).Multiply(u));
                    v = v.Add(M[bw_i.boneIndex2].Multiply(bw_i.weight2).Multiply(u));
                    v = v.Add(M[bw_i.boneIndex3].Multiply(bw_i.weight3).Multiply(u));
                }
                vs_mat[i] = v;
                u_.SetRow(i, u.Column(0));
                v_.SetRow(i, v.Column(0));
                vs[i] = new Vector3(v[0, 0], v[1, 0], v[2, 0]);
            }
            _Smooth_vs_mat = IterativeCalcB(v_, iterations);
            _Smooth_us = IterativeCalcB(u_, iterations);
        }
        public Matrix<float> IterativeCalcB(Matrix<float> m, int iteration)
        {
            if (iteration == 0) return m;
            iteration -= 1;
            return B.Solve(IterativeCalcB(m, iteration));
        }
        public void CalcuateNormalizedLaplace()
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

        }
    }

}