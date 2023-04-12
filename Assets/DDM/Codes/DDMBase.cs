using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Single;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Unity.VisualScripting;
using UnityEngine;
namespace DDM_Impl
{
    public class DDMBase : MonoBehaviour
    {
        public SkinnedMeshRenderer TargetSMR;
        public ComputeShader _CShader;
        ComputeBuffer OutVert;
        ComputeBuffer OutNor;
        ComputeBuffer CB_Rs;
        ComputeBuffer CB_p_s;
        ComputeBuffer CB_q_s;
        ComputeBuffer CB_USVs;
        Float3x3[] USVs;
        ComputeBuffer BoneMatrixs;
        ComputeBuffer BoneBindings;
        ComputeBuffer CBPsis;
        public bool UseComputeShader;
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
        Float4[] BoneBinding;
        Float3[] q_s;
        Float3[] p_s;
        int __KERNEL_FIRST_PASS;
        int __KERNEL_SECOND_PASS;
        bool FinalCS;
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
            for (int i = 0; i < Bones.Length; i++)
            {
                boneM_[i] = Matrix<float>.Build.Dense(4, 4);
            }
            Precompute();
            FinalCS = UseComputeShader && Constants.SupportComputerShader;
            Debug.Log($"Final Use Compute Shader:{FinalCS}");
            if (FinalCS)
                SetupComputerShader();
        }
        private void OnDestroy()
        {
            try
            {
                CBPsis.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                CB_p_s.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                CB_q_s.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                CB_Rs.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                CB_USVs.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                OutNor.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                OutVert.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                vert.Dispose();
            }
            catch (System.Exception)
            {
            }
            try
            {
                nors.Dispose();
            }
            catch (System.Exception)
            {
            }
        }
        ComputeBuffer vert;
        ComputeBuffer nors;
        void SetupComputerShader()
        {
            __KERNEL_FIRST_PASS = _CShader.FindKernel("FirstPass");
            __KERNEL_SECOND_PASS = _CShader.FindKernel("SecondPass");
            vert = new ComputeBuffer(Vertices.Length, sizeof(float) * 3);
            nors = new ComputeBuffer(Vertices.Length, sizeof(float) * 3);
            OutVert = new ComputeBuffer(Vertices.Length, sizeof(float) * 3);
            OutNor = new ComputeBuffer(Vertices.Length, sizeof(float) * 3);
            vert.SetData(Vertices);
            nors.SetData(Normals);
            CBPsis = new ComputeBuffer(CurrentMesh.boneWeights.Length * 4, sizeof(float) * 4 * 4);
            CSPsis = new Matrix4x4[CurrentMesh.boneWeights.Length * 4];
            Rs = new Float3x3[CurrentMesh.vertices.Length];
            CB_q_s = new ComputeBuffer(CurrentMesh.vertices.Length, sizeof(float) * 3);
            CB_p_s = new ComputeBuffer(CurrentMesh.vertices.Length, sizeof(float) * 3);
            CB_Rs = new ComputeBuffer(CurrentMesh.vertices.Length, sizeof(float) * 3 * 3);
            CB_USVs = new ComputeBuffer(CurrentMesh.vertices.Length, sizeof(float) * 3 * 3);
            BoneMatrixs = new ComputeBuffer(CurrentMesh.vertices.Length, sizeof(float) * 4 * 4);
            BoneBinding = new Float4[CurrentMesh.boneWeights.Length];
            BoneBindings = new ComputeBuffer(BoneBinding.Length, sizeof(float) * 4);
            q_s = new Float3[CurrentMesh.vertices.Length];
            p_s = new Float3[CurrentMesh.vertices.Length];

            for (int i = 0; i < CurrentMesh.boneWeights.Length; i++)
            {
                BoneBinding[i] = new Float4
                {
                    x = CurrentMesh.boneWeights[i].weight0,
                    y = CurrentMesh.boneWeights[i].weight1,
                    z = CurrentMesh.boneWeights[i].weight2,
                    w = CurrentMesh.boneWeights[i].weight3
                };

            }
            BoneBindings.SetData(BoneBinding);
            USVs = new Float3x3[CurrentMesh.vertices.Length];
            for (int i = 0; i < CurrentMesh.boneWeights.Length; i++)
            {
                for (int x = 0; x < 4; x++)
                {
                    CSPsis[i * 4 + x] = Psis[i, x].ToMatrix();
                }
            }
            //CBPsis.Set()
            CB_USVs.SetData(USVs);
            CBPsis.SetData(CSPsis);
            //_CShader.SetBuffer(__KERNEL_FIRST_PASS, "Vertices", vert);
            _CShader.SetBuffer(__KERNEL_FIRST_PASS, "Psis", CBPsis);
            _CShader.SetBuffer(__KERNEL_FIRST_PASS, "BoneBinding", BoneBindings);
            _CShader.SetBuffer(__KERNEL_FIRST_PASS, "USVs", CB_USVs);
            _CShader.SetBuffer(__KERNEL_FIRST_PASS, "BoneM", BoneMatrixs);
            _CShader.SetBuffer(__KERNEL_FIRST_PASS, "q_s", CB_q_s);
            _CShader.SetBuffer(__KERNEL_FIRST_PASS, "p_s", CB_p_s);

            _CShader.SetBuffer(__KERNEL_SECOND_PASS, "Rs", CB_Rs);
            _CShader.SetBuffer(__KERNEL_SECOND_PASS, "q_s", CB_q_s);
            _CShader.SetBuffer(__KERNEL_SECOND_PASS, "p_s", CB_p_s);
            _CShader.SetBuffer(__KERNEL_SECOND_PASS, "O_Vert", OutVert);
            _CShader.SetBuffer(__KERNEL_SECOND_PASS, "O_Nor", OutNor);
            _CShader.SetBuffer(__KERNEL_SECOND_PASS, "Vertices", vert);
            _CShader.SetBuffer(__KERNEL_SECOND_PASS, "Normals", nors);
        }
        Matrix<float> B;
        Matrix<float> A;
        Matrix<float> A_p;
        Matrix<float> B_p;
        Matrix<float>[,] Psis;
        Matrix4x4[] CSPsis;
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
            Matrix<float>[,] _Psis = new Matrix<float>[CurrentMesh.boneWeights.Length, 4];

            Matrix<float> Psi_i_0;
            Matrix<float> Psi_i_1;
            Matrix<float> Psi_i_2;
            Matrix<float> Psi_i_3;
            for (int i = 0; i < CurrentMesh.vertices.Length; i++)
            {
                {
                    var wei = CurrentMesh.boneWeights[i];
                    var u_col = u_.Row(i);
                    var u_k_calced = u_col.OuterProduct(u_col);
                    {
                        Psis[i, 0] = (wei.weight0 * u_k_calced);// [0, 0];
                        Psis[i, 1] = (wei.weight1 * u_k_calced);// [0, 0];
                        Psis[i, 2] = (wei.weight2 * u_k_calced);// [0, 0];
                        Psis[i, 3] = (wei.weight3 * u_k_calced);// [0, 0];
                    }
                }
                //Psis[i, 0] = Psi_i_0.Duplicate();
                //Psis[i, 1] = Psi_i_1.Duplicate();
                //Psis[i, 2] = Psi_i_2.Duplicate();
                //Psis[i, 3] = Psi_i_3.Duplicate();
            }
            for (int it_count = 0; it_count < iterations; it_count++)
            {
                var tmp = _Psis.Duplicate();
                _Psis = Psis.Duplicate();
                Psis = tmp;
                //Psis = tmp;
                for (int i = 0; i < CurrentMesh.vertices.Length; i++)
                {
                    // i : i;
                    {
                        Psi_i_0 = Matrix<float>.Build.Dense(4, 4, 0);
                        Psi_i_1 = Matrix<float>.Build.Dense(4, 4, 0);
                        Psi_i_2 = Matrix<float>.Build.Dense(4, 4, 0);
                        Psi_i_3 = Matrix<float>.Build.Dense(4, 4, 0);
                    }
                    for (int k = 0; k < CurrentMesh.vertices.Length; k++)
                    {
                        var AorB = B_p[i, k];
                        if (AorB != 0f)
                        {

                            Psi_i_0 += (AorB * _Psis[k, 0]);
                            Psi_i_1 += (AorB * _Psis[k, 1]);
                            Psi_i_2 += (AorB * _Psis[k, 2]);
                            Psi_i_3 += (AorB * _Psis[k, 3]);
                        }
                    }
                    Psis[i, 0] = Psi_i_0.Duplicate();
                    Psis[i, 1] = Psi_i_1.Duplicate();
                    Psis[i, 2] = Psi_i_2.Duplicate();
                    Psis[i, 3] = Psi_i_3.Duplicate();
                }
            }
        }
        // Update is called once per frame
        void Update()
        {
            OnFrame();
        }
        Matrix4x4[] boneM;
        Float3x3[] Rs;
        Matrix<float>[] boneM_;
        Matrix<float> Q = Matrix<float>.Build.Dense(3, 3);
        Matrix<float> pt = Matrix<float>.Build.Dense(1, 3);
        Matrix<float> p = Matrix<float>.Build.Dense(3, 1);
        Matrix<float> q = Matrix<float>.Build.Dense(3, 1);
        //[MethodImpl(MethodImplOptions.AggressiveInlining)]
        void OnComputeShader()
        {
            for (int i = 0; i < Bones.Length; i++)
            {
                boneM[i] = Bones[i].GlobalToMatrix() * BindPoses[i];
                boneM[i].ToMatrix(boneM_[i]);
            }
            //First Pass

            BoneMatrixs.SetData(boneM);

            _CShader.Dispatch(__KERNEL_FIRST_PASS, Vertices.Length, 1, 1);
            //Managed Pass
            CB_USVs.GetData(USVs);

            CB_p_s.GetData(p_s);
            CB_q_s.GetData(q_s);
            int ErrorCount = 0;
            for (int i = 0; i < USVs.Length; i++)
            {
                try
                {
                    var USV = USVs[i].ToMatrix();
                    var svd = USV.Svd();
                    //Debug.Log(USV);
                    //Rs[i] = Float3x3.FromMatrix(svd.U * svd.VT);
                    var q = q_s[i].ToVertical();
                    var p = q_s[i].ToVertical();
                    var R = svd.U * svd.VT;
                    var T = q - R * p;
                    Matrix4x4 trans = Matrix4x4.zero;
                    for (int x = 0; x < 3; x++)
                    {
                        for (int y = 0; y < 3; y++)
                        {
                            trans[x, y] = R[x, y];
                        }
                    }
                    trans[0, 3] = T[0, 0];
                    trans[1, 3] = T[1, 0];
                    trans[2, 3] = T[2, 0];
                    trans[3, 3] = 1;

                    AlteredVertices[i] = trans.MultiplyPoint3x4(Vertices[i]);
                    AlteredNormals[i] = trans.MultiplyVector(Normals[i]);
                }
                catch (System.Exception)
                {
                    ErrorCount++;
                }
            }
            //Debug.Log($"Error Count : {ErrorCount} out of {USVs.Length} USVs.");
            //CB_Rs.SetData(Rs);
            //CB_p_s.SetData(p_s);
            //CB_q_s.SetData(q_s);
            ////
            //_CShader.Dispatch(__KERNEL_SECOND_PASS, Vertices.Length, 1, 1);
            //OutVert.GetData(AlteredVertices);
            //OutNor.GetData(AlteredNormals);

            CurrentMesh.vertices = AlteredVertices;
            CurrentMesh.normals = AlteredNormals;
            Graphics.DrawMesh(CurrentMesh, Matrix4x4.identity, TargetSMR.sharedMaterial, 0);
        }
        void OnFrame()
        {
            if (FinalCS)
            {
                OnComputeShader();
                return;
            }
            for (int i = 0; i < Bones.Length; i++)
            {
                boneM[i] = Bones[i].GlobalToMatrix() * BindPoses[i];
                boneM[i].ToMatrix(boneM_[i]);
            }
            Matrix<float> Omega;
            for (int i = 0; i < Vertices.Length; i++)
            {
                var wei = CurrentMesh.boneWeights[i];
                Omega = (boneM_[wei.boneIndex0]) * Psis[i, 0] +
                   (boneM_[wei.boneIndex1]) * Psis[i, 1] +
                   (boneM_[wei.boneIndex2]) * Psis[i, 2] +
                   (boneM_[wei.boneIndex3]) * Psis[i, 3];

                for (int x = 0; x < 3; x++)
                {
                    for (int y = 0; y < 3; y++)
                    {
                        Q[x, y] = Omega[x, y];
                    }
                }
                q[0, 0] = Omega[0, 3];
                q[1, 0] = Omega[1, 3];
                q[2, 0] = Omega[2, 3];
                pt[0, 0] = Omega[3, 0];
                pt[0, 1] = Omega[3, 1];
                pt[0, 2] = Omega[3, 2];
                pt.Transpose(p);
                var USV = Q - q * pt;
                try
                {
                    MathNet.Numerics.LinearAlgebra.Factorization.Svd<float> SVD = USV.Svd();
                    var R = SVD.U * SVD.VT;
                    var T = q - R * p;
                    Matrix4x4 trans = Matrix4x4.zero;
                    for (int x = 0; x < 3; x++)
                    {
                        for (int y = 0; y < 3; y++)
                        {
                            trans[x, y] = R[x, y];
                        }
                    }
                    trans[0, 3] = T[0, 0];
                    trans[1, 3] = T[1, 0];
                    trans[2, 3] = T[2, 0];
                    trans[3, 3] = 1;

                    AlteredVertices[i] = trans.MultiplyPoint3x4(Vertices[i]);
                    AlteredNormals[i] = trans.MultiplyVector(Normals[i]);
                }
                catch (System.Exception)
                {
                    //Debug.Log($"Failed on SVD:{USV}");
                }

            }
            //if (GC_Count == 10)
            //{
            //    System.GC.Collect(0, System.GCCollectionMode.Optimized, false);
            //    GC_Count = 0;
            //}
            //GC_Count++;
            CurrentMesh.vertices = AlteredVertices;
            CurrentMesh.normals = AlteredNormals;
            Graphics.DrawMesh(CurrentMesh, Matrix4x4.identity, TargetSMR.sharedMaterial, 0);
        }
        int GC_Count = 0;
        Matrix<float> u_;
        Matrix<float>[] u__;
        Matrix<float> v_;
        Matrix<float> normalizedLaplacian;
        Matrix<float> _Smooth_vs_mat;
        Matrix<float> _Smooth_us;

        public void CalcAMatrix(float Lambda)
        {
            A = Matrix<float>.Build.DiagonalIdentity(LaplacianMatrix.RowCount, LaplacianMatrix.ColumnCount) - (LaplacianMatrix * Lambda);
            A_p = Matrix<float>.Build.Sparse(A.RowCount, A.ColumnCount);
            A.Power(iterations, A_p);
        }
        public void CalcBMatrix(float Lambda)
        {
            var w = LaplacianMatrix.RowCount;
            var h = LaplacianMatrix.ColumnCount;
            //var DL = Matrix<float>.Build.DiagonalIdentity(w, h) * Lambda;
            B = (Matrix<float>.Build.DiagonalIdentity(w, h) - (Lambda * LaplacianMatrix));//.Transpose(); //* DL);
            B_p = B;
            return;
            //B_p = Matrix<float>.Build.DiagonalIdentity(w, h);// + (Lambda * LaplacianMatrix); //* DL);
            //for (int i = 0; i < iterations; i++)
            //{
            //    B_p = iB.Solve(B_p);
            //}

            //B_p = Matrix<float>.Build.Sparse(w,h);
            //B.Power(iterations, B_p);
            //B_p = B_p.Inverse();
            //B_p = Matrix<float>.Build.Dense(B.RowCount, B.ColumnCount, 1);
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
            u__ = new Matrix<float>[vert.Length];
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
                u__[i] = u;
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
            //return m * B;
            //return m * (Mathf.Pow(B, iteration));
            if (iteration == 1) return m;
            iteration -= 1;
            return m + IterativeSolveAx_B(m, iteration, B) / B;
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
            var adj = MatrixUtils.BuildAdjacencyMatrix(Vertices, CurrentMesh.triangles, 32);
            normalizedLaplacian = LaplacianMatrix = MatrixUtils.BuildLaplacianMatrixFromAdjacentMatrix(Vertices.Length, adj, true);
            return;
            AdjacencyMatrix _AdjacencyMatrix = AdjacencyMatrix.FromMesh(CurrentMesh);
            var deg = DegMatrix.FromMesh(CurrentMesh);
            Matrix<float> DegMat = deg.Matrix;
            //DegMat = DegMat.Multiply(3);
            LaplacianMatrix = DegMat - (_AdjacencyMatrix);

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
    public static class Constants
    {
        public static bool SupportComputerShader = SystemInfo.supportsComputeShaders;
    }
}