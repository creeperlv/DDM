# DDM on Unity3D

This project is a shcoolword. Aims on implement DDM in pure C# and Compute Shader, also provide a prototype bake tool.

This project is very experimental, should never directly be used for production!

# Further Works

Modify the way calculations work to avoid GC. e.g.: use Matrix<T>.Multiply(Matrix<T> other, Matrix<T> result) instead of overloaded * operator in Matrix<T>.

Properly implement a SVD method in Compute Shader.

Suggest to load Baked data on the game/scene load stage.


# Credit

[Math.NET](https://numerics.mathdotnet.com/)

[Uniry3D](https://unity.com/)

[DDM](https://www.ea.com/seed/news/siggraph2019-direct-delta-mush)