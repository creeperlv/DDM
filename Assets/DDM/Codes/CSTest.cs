using System.Collections;
using System.Collections.Generic;
using UnityEngine;
namespace DDM_Impl
{
    public class CSTest : MonoBehaviour
    {
        int f;
        int s;
        public ComputeShader Shader;
        ComputeBuffer A;
        ComputeBuffer B;
        ComputeBuffer C;
        Float3x3[] data;
        Float3x3[] data2;
        Float3x3[] data3;
        private void OnDestroy()
        {
            A.Dispose(); B.Dispose(); C.Dispose();
        }
        void Start()
        {
            f = Shader.FindKernel("First");
            s = Shader.FindKernel("Second");
            int LENGTH = 16;
            data = new Float3x3[LENGTH];
            data2 = new Float3x3[LENGTH];
            data3 = new Float3x3[LENGTH];
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = new Float3x3();
                for (int x = 0; x < 3; x++)
                {
                    for (int y = 0; y < 3; y++)
                    {
                        data[i][x, y] = i * 100 + x * 10 + y;
                    }
                }
            }
            A = new ComputeBuffer(data.Length, sizeof(float) * 3 * 3);
            B = new ComputeBuffer(data.Length, sizeof(float) * 3 * 3);
            C = new ComputeBuffer(data.Length, sizeof(float) * 3 * 3);
            A.SetData(data);
            B.SetData(data2);
            C.SetData(data3);
            Shader.SetBuffer(f, "Rs", A);
            Shader.SetBuffer(f, "Out_Rs", B);
            Shader.SetBuffer(s, "Out_Rs", B);
            Shader.SetBuffer(s, "Out_Rs_2", C);
        }

        // Update is called once per frame
        void Update()
        {
            Shader.Dispatch(f, data.Length, 1, 1);
            B.GetData(data2);
            foreach (var item in data2)
            {
                Debug.Log("FIRST:" + item);
            }
            Shader.Dispatch(s, data.Length, 1, 1);
            C.GetData(data3);
            foreach (var item in data3)
            {
                Debug.Log("SECOND:" + item);
            }
        }
    }

}