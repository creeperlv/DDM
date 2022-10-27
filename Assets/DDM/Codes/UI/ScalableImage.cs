using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ScalableImage : MonoBehaviour
{
    public RectTransform BindedTransform;
    public float Scale { get => BindedTransform.localScale.x; set => BindedTransform.localScale = new Vector3(value, value, value); }
    // Start is called before the first frame update
    void Start()
    {
        
    }
    public void SetScale(float s) => Scale = s;
    // Update is called once per frame
    void Update()
    {
        
    }
}
