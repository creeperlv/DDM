using DDM_Impl;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Threading.Tasks;
using Unity.Jobs.LowLevel.Unsafe;
using UnityEditor;
using UnityEditor.SceneManagement;
using UnityEngine;
namespace DDM_Impl.Editors
{
    [CustomEditor(typeof(DDMBase))]
    public class DDMPrecomputeBaker : Editor
    {
        bool IsRunning = false;
        DDMBase Target_Base;
        void Bake(string Path)
        {
            IsRunning = true;
            Debug.Log("Start Baking");
            Target_Base.Precompute();
            if (File.Exists(Path))
            {
                File.Delete(Path);
            }
            FileInfo fileInfo = new FileInfo(Path);
            Debug.Log("Write to:" + fileInfo.FullName);
            using (var fs = fileInfo.Create())
            {
                {
                    var W = Target_Base.Psis.GetLength(0);
                    for (int i = 0; i < W; i++)
                    {
                        for (int ii = 0; ii < 4; ii++)
                        {
                            for (int x = 0; x < 4; x++)
                            {
                                for (int y = 0; y < 4; y++)
                                {
                                    fs.Write(BitConverter.GetBytes(Target_Base.Psis[i, ii][x, y]));
                                }
                            }
                        }
                    }
                }
            }
            Debug.Log("Baking ended.");

            IsRunning = false;
        }
        public override void OnInspectorGUI()
        {

            //base.OnInspectorGUI();
            DrawDefaultInspector();
            Target_Base = (DDMBase)target;
            var path = GUILayout.TextField(EditorSceneManager.GetActiveScene().path.Replace(".unity", "") + "/" + target.GetInstanceID() + ".txt");
            if (IsRunning == false)
            {
                if (GUILayout.Button(new GUIContent { text = "Bake" }))
                {
                    Target_Base.PreBake();

                    Task.Run(() =>
                    {
                        Bake(path);
                    }
                    );
                }
            }
            else
            {
                var rect0 = EditorGUILayout.GetControlRect(false, EditorGUIUtility.singleLineHeight);
                var rect1 = EditorGUILayout.GetControlRect(false, EditorGUIUtility.singleLineHeight);

                if (Target_Base.PROG__STAGE != 0)
                {
                    EditorGUI.ProgressBar(rect0, Mathf.Min(Mathf.Max(0, ((float)Target_Base.PROG__STAGE / 4f)), 1), $"Stage at:{Target_Base.PROG__STAGE}");
                }
                else
                    EditorGUI.ProgressBar(rect0, 0, $"Stage at:{Target_Base.PROG__STAGE}");
                if (Target_Base.PROG__MAX != 0)
                {
                    float prog= (float)Target_Base.PROG__VALUE / ((float)Target_Base.PROG__MAX);
                    EditorGUI.ProgressBar(rect1,Mathf.Min(1,Mathf.Max(0, prog)), "Current Stage:"+prog*100+"%");
                }
                else
                {
                    EditorGUI.ProgressBar(rect1, 0.5f, "Current Stage (Undetermined)");
                }
            }
        }
    }

}
