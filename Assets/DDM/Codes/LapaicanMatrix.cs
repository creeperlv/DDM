using UnityEngine;

namespace DDM_Impl
{
    public class LapaicanMatrix
    {
        public int n;
        public int[] data;
        public LapaicanMatrix(int n)
        {
            this.n = n;
            data=new int[n*n];
        }
    }
    public class DebugUtilities
    {
        public static Texture2D FromMatrix(float[,] matrix,int Width,int Height,int offset=0,float intensity=1f)
        {
            var texture=new Texture2D(Width,Height, TextureFormat.RGBA32,true,true,false);
            texture.minimumMipmapLevel = 15;
            texture.filterMode = FilterMode.Point;
            for (int i_x = 0; i_x < Width; i_x++)
            {
                for (int i_y = 0; i_y < Height; i_y++)
                {
                    float i = matrix[i_x, i_y];
                    float c = (i+offset) * intensity;
                    texture.SetPixel(i_x, Height-i_y, new Color(c, c, c, (i == 0?0:1)));
                }
            }
            texture.Apply();
            return texture;
        }
    }
}