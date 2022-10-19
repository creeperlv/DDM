using System.Collections;
using System.Collections.Generic;
using UnityEngine;
namespace DDM_Impl
{
    public class DDM_Surface : MonoBehaviour
    {
        public ComputeShader DDM_Stage_0;
        public Mesh Mesh;
        public DDMMesh DDMMesh;
        public SkinnedMeshRenderer SkinRenderer;
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
            DDMMesh = new DDMMesh(Mesh.vertexCount);
        }
        // Update is called once per frame
        void Update()
        {

        }
    }

    public struct DDMMesh
    {
        public DDMMesh(int Count)
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