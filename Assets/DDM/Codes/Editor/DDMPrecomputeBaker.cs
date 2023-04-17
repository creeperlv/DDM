using DDM_Impl;
using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
namespace DDM_Impl.Editors
{
    [CustomEditor(typeof(DDMBase))]
    public class DDMPrecomputeBaker : Editor
    {
        public override void OnInspectorGUI()
        {
            //base.OnInspectorGUI();
            DrawDefaultInspector();


            DDMBase Target_Base = (DDMBase)target;
            if (GUILayout.Button(new GUIContent { text = "Bake" }))
            {
                Target_Base.Precompute();
            }
        }
    }

}
