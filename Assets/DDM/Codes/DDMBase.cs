using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
namespace DDM_Impl
{
    public class DDMBase : MonoBehaviour
    {
        public SkinnedMeshRenderer TargetSMR;
        public float Lambda;
        public int iterations = 1;
        public Transform Root;
        Mesh CurrentMesh;
        Vector3[] Vertices;
        Vector3[] Normals;
        Matrix4x4[] BindPoses;
        Vector3[] AlteredVertices;
        Vector3[] AlteredNormals;
        Transform[] Bones;
        void Start()
        {
            CurrentMesh = Mesh.Instantiate(TargetSMR.sharedMesh);
            var MF = this.GetComponent<MeshFilter>();
            MF.mesh = CurrentMesh;
            BindPoses = TargetSMR.sharedMesh.bindposes;
            //var obj = new GameObject("DDM Render");
            //obj.transform.parent = TargetSMR.transform.parent;
            //obj.transform.position= new Vector3(0,0,0);
            //obj.transform.rotation= Quaternion.identity;
            //obj.transform.localScale=TargetSMR.transform.localScale;
            //var MF=obj.AddComponent<MeshFilter>();
            //MF.mesh = CurrentMesh;
            CurrentMesh.MarkDynamic();
            //var MR=obj.AddComponent<MeshRenderer>();
            //MR.materials=TargetSMR.materials;
            //TargetSMR.sharedMesh.MarkDynamic();
            Vertices = CurrentMesh.vertices;
            Normals = CurrentMesh.normals;
            AlteredVertices = CurrentMesh.vertices;
            AlteredNormals = CurrentMesh.normals;
            Bones = TargetSMR.bones;
            boneM = new Matrix4x4[Bones.Length];
            boneM_ = new Matrix<float>[Bones.Length];
            Precompute();
        }

        Matrix<float> B;
        Matrix<float> B_p;
        Matrix<float>[,] Psis;
        void Precompute()
        {
            CalcuateNormalizedLaplace();

            CalcBMatrix(Lambda);
            BuildUs();
            CalcPsis();
        }
        void CalcPsis()
        {
            //Debug.Log($"u_:Row->{u_.RowCount},Column->{u_.ColumnCount}");
            Psis = new Matrix<float>[CurrentMesh.boneWeights.Length, 4];
            for (int i = 0; i < CurrentMesh.vertices.Length; i++)
            {
                // i : i;
                Matrix<float> Psi_i_0 = Matrix<float>.Build.Sparse(4, 4, 0);
                Matrix<float> Psi_i_1 = Matrix<float>.Build.Sparse(4, 4, 0);
                Matrix<float> Psi_i_2 = Matrix<float>.Build.Sparse(4, 4, 0);
                Matrix<float> Psi_i_3 = Matrix<float>.Build.Sparse(4, 4, 0);
                for (int k = 0; k < CurrentMesh.vertices.Length; k++)
                {
                    var wei = CurrentMesh.boneWeights[k];
                    //var v=Vertices[k];
                    //var u_k = Matrix<float>.Build.DenseOfColumns(new float[1][]{ new float[4] { v.x, v.y, v.z, 1 } });
                    var u_k = u_.Row(k).ToColumnMatrix();
                    var u_k_calced = u_k * u_k.Transpose();
                    if (B[k, i] == 0f)
                    {
                        Psi_i_0 += (wei.weight0 * u_k_calced);// [0, 0];
                        Psi_i_1 += (wei.weight1 * u_k_calced);// [0, 0];
                        Psi_i_2 += (wei.weight2 * u_k_calced);// [0, 0];
                        Psi_i_3 += (wei.weight3 * u_k_calced);// [0, 0];

                    }
                    else
                    {
                        Psi_i_0 += IterativeSolveAx_B(wei.weight0 * u_k_calced, iterations, B[k, i]);
                        Psi_i_1 += IterativeSolveAx_B(wei.weight1 * u_k_calced, iterations, B[k, i]);
                        Psi_i_2 += IterativeSolveAx_B(wei.weight2 * u_k_calced, iterations, B[k, i]);
                        Psi_i_3 += IterativeSolveAx_B(wei.weight3 * u_k_calced, iterations, B[k, i]);

                    }
                }
                Psis[i, 0] = Psi_i_0;
                Psis[i, 1] = Psi_i_1;
                Psis[i, 2] = Psi_i_2;
                Psis[i, 3] = Psi_i_3;
            }
            Debug.Log("Psi 0,0:");
            Debug.Log(Psis[0, 0]);
        }
        // Update is called once per frame
        void Update()
        {
            OnFrame();
        }
        Matrix4x4[] boneM;
        Matrix<float>[] boneM_;
        void OnFrame()
        {
            for (int i = 0; i < Bones.Length; i++)
            {
                boneM[i] = Bones[i].RootSpaceMatrix (Root)* BindPoses[i];
                boneM_[i] = boneM[i].ToMatrix();
            }
            for (int i = 0; i < Vertices.Length; i++)
            {
                var wei = CurrentMesh.boneWeights[i];

                Matrix<float> Omega = (boneM_[wei.boneIndex0]) * Psis[i, 0] +
                    (boneM_[wei.boneIndex1]) * Psis[i, 1] +
                    (boneM_[wei.boneIndex2]) * Psis[i, 2] +
                    (boneM_[wei.boneIndex3]) * Psis[i, 3];
                Matrix<float> Q = Omega.SubMatrix(0, 3, 0, 3);
                Matrix<float> pt = Omega.SubMatrix(3, 1, 0, 3);
                Matrix<float> p = pt.Transpose();
                //Matrix<float> pt=p.Transpose();
                Matrix<float> q = Omega.SubMatrix(0, 3, 3, 1);
                //Debug.Log($"P:Row->{pt.RowCount},Column->{pt.ColumnCount}");
                //Debug.Log($"Pt:Row->{pt.RowCount},Column->{pt.ColumnCount}");
                //Debug.Log($"q:Row->{q.RowCount},Column->{q.ColumnCount}");
                var USV = Q - q * pt;
                //Debug.Log($"USV:Row->{USV.RowCount},Column->{USV.ColumnCount}");
                try
                {
                    var SVD = USV.Svd();
                    var R = SVD.U * SVD.VT;
                    var T = q - R * p;
                    //Debug.Log($"q:{q}");
                    //Debug.Log($"p:{p}");
                    //Debug.Log(T);
                    //Matrix4x4 matrix4X4 = Matrix4x4.Translate(MatrixUtils.ToVector3(T));
                    //Matrix4x4 _r = new Matrix4x4(
                    //    new Vector4 { x = R[0, 0], y = R[1, 0], z = R[2, 0], w = T[0, 0] },
                    //    new Vector4 { x = R[0, 1], y = R[1, 1], z = R[2, 1], w = T[1, 0] },
                    //    new Vector4 { x = R[0, 2], y = R[1, 2], z = R[2, 2], w = T[2, 0] },
                    //    new Vector4 { x = 0, y = 0, z = 0, w = 1 }
                    //);
                    Matrix4x4 trans = Matrix4x4.Scale(Vector3.one);
                    for (int x = 0; x < 3; ++x)
                    {
                        for (int y = 0; y < 3; ++y)
                        {
                            trans[x, y] = R[x, y];
                        }
                    }
                    trans[0, 3] = T[0,0];
                    trans[1, 3] = T[1,0];
                    trans[2, 3] = T[2, 0];
                    trans[3, 3] = 1;
                    //new Vector4 { x = R[0, 0], y = R[0, 1], z = R[0, 2], w = T[0, 0] },
                    //    new Vector4 { x = R[1, 1], y = R[1, 1], z = R[1, 2], w = T[1, 0] },
                    //    new Vector4 { x = R[2, 2], y = R[2, 1], z = R[2, 2], w = T[2, 0] },
                    
                    AlteredVertices[i] =trans.MultiplyPoint(Vertices[i]);
                    AlteredNormals[i] = trans.MultiplyVector(Normals[i]);
                }
                catch (System.Exception)
                {
                    //Debug.Log($"Failed on SVD:{USV}");
                }

            }
            CurrentMesh.vertices = AlteredVertices;
            CurrentMesh.normals = AlteredNormals;
            //Graphics.DrawMesh(CurrentMesh, Matrix4x4.identity, TargetSMR.sharedMaterial, 0);
        }
        Matrix<float> u_;
        Matrix<float> v_;
        Matrix<float> normalizedLaplacian;
        Matrix<float> _Smooth_vs_mat;
        Matrix<float> _Smooth_us;

        public void CalcBMatrix(float Lambda)
        {
            B = Matrix<float>.Build.DiagonalIdentity(LaplacianMatrix.RowCount, LaplacianMatrix.ColumnCount)+(LaplacianMatrix*Lambda);
            //B_p = Matrix<float>.Build.Sparse(B.RowCount, B.ColumnCount);
            //B.Power(iterations,B_p);
        }
        public void BuildUs()
        {
            /* 
             * Equation.1 
             */
            //Matrix<float>[] M = new Matrix<float>[Bones.Length];
            //var bp = CurrentMesh.bindposes;
            //for (int i = 0; i < Bones.Length; i++)
            //{
            //    M[i] = Bones[i].worldToLocalMatrix.ToMatrix();

            //}
            //var bw = CurrentMesh.boneWeights;
            var vert = CurrentMesh.vertices;
            u_ = Matrix<float>.Build.Sparse(vert.Length, 4);
            //v_ = Matrix<float>.Build.Sparse(vert.Length, 4);
            //Vector3[] vs = new Vector3[vert.Length];
            //Matrix<float>[] vs_mat = new Matrix<float>[vert.Length];
            //_Smooth_vs_mat = Matrix<float>.Build.Sparse(vert.Length, 4);
            //Matrix<float>[] us = new Matrix<float>[vert.Length];
            //_Smooth_us = Matrix<float>.Build.Sparse(vert.Length, 4);
            for (int i = 0; i < vert.Length; i++)
            {
                //var bw_i = bw[i];
                var vert_i = vert[i];

                Matrix<float> u = Matrix<float>.Build.DenseOfColumnArrays(new float[] { vert_i.x, vert_i.y, vert_i.z, 1 });
                //us[i] = u;
                //Matrix<float> v;
                //{
                //    v = M[bw_i.boneIndex0].Multiply(bw_i.weight0).Multiply(u);
                //    v = v.Add(M[bw_i.boneIndex1].Multiply(bw_i.weight1).Multiply(u));
                //    v = v.Add(M[bw_i.boneIndex2].Multiply(bw_i.weight2).Multiply(u));
                //    v = v.Add(M[bw_i.boneIndex3].Multiply(bw_i.weight3).Multiply(u));
                //}
                //vs_mat[i] = v;
                u_.SetRow(i, u.Column(0));
                //v_.SetRow(i, v.Column(0));
                //vs[i] = new Vector3(v[0, 0], v[1, 0], v[2, 0]);
            }
            //_Smooth_vs_mat = IterativeCalcB(v_, iterations);
            //_Smooth_us = IterativeCalcB(u_, iterations);
        }
        public Matrix<float> IterativeSolveAx_B(Matrix<float> m, int iteration, float B)
        {
            if (iteration == 1) return m;
            iteration -= 1;
            return IterativeSolveAx_B(m, iteration, B) / B;
        }
        public Matrix<float> IterativeCalcB(Matrix<float> m, int iteration)
        {
            if (iteration == 0) return m;
            iteration -= 1;
            return B.Solve(IterativeCalcB(m, iteration));
        }
        Matrix<float> LaplacianMatrix;
        public void CalcuateNormalizedLaplace()
        {
            //BuildD();
            //LaplacianMatrix=Matrix<float>.Build.Sparse(Vertices.Length, Vertices.Length);
            //for (int i = 0; i < CurrentMesh.triangles.Length; i += 3)
            //{
            //    var a = CurrentMesh.triangles[i];
            //    var b = CurrentMesh.triangles[i+1];
            //    var c = CurrentMesh.triangles[i+2];
            //    LaplacianMatrix[a, b] = -1;
            //    LaplacianMatrix[a, c] = -1;
            //    LaplacianMatrix[b, c] = -1;
            //}
            //for (int i = 0; i < D.Length; i++)
            //{
            //    LaplacianMatrix[i, i] = D[i];

            //}
            AdjacencyMatrix _AdjacencyMatrix = AdjacencyMatrix.FromMesh(CurrentMesh);
            var deg = DegMatrix.FromMesh(CurrentMesh);
            Matrix<float> DegMat = deg.Matrix;
            DegMat = DegMat.Multiply(3);
            LaplacianMatrix = DegMat.Subtract(_AdjacencyMatrix);
            var DL = LaplacianMatrix.Diagonal();
            var DL_m = Matrix<float>.Build.Diagonal(DL.ToArray());
            var DL_m_I = DL_m.Inverse();
            normalizedLaplacian = LaplacianMatrix.Multiply(DL_m_I);

        }
        int[] D;
        void BuildD()
        {
            D = new int[Vertices.Length];
            for (int i = 0; i < CurrentMesh.triangles.Length - 2; i += 3)
            {
                D[CurrentMesh.triangles[i]] += 2;
                D[CurrentMesh.triangles[i + 1]] += 2;
                D[CurrentMesh.triangles[i + 2]] += 2;
            }
        }
    }

}