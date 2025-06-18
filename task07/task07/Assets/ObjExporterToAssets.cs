using UnityEngine;
using System.IO;
using System.Text;

public class ObjExporterToAssets : MonoBehaviour
{

    static public string MeshToObj(Mesh mesh, Transform transform)
    {
        StringBuilder sb = new StringBuilder();

        sb.AppendLine($"# Exported from Unity: {mesh.name}");

        // 頂点
        foreach (Vector3 v in mesh.vertices)
        {
            Vector3 worldV = transform.TransformPoint(v); // ワールド座標に変換
            sb.AppendLine($"v {worldV.x} {worldV.y} {worldV.z}");
        }

        // UV
        foreach (Vector2 uv in mesh.uv)
        {
            sb.AppendLine($"vt {uv.x} {uv.y}");
        }

        // 法線
        foreach (Vector3 n in mesh.normals)
        {
            Vector3 worldN = transform.TransformDirection(n);
            sb.AppendLine($"vn {worldN.x} {worldN.y} {worldN.z}");
        }

        sb.AppendLine("g " + mesh.name);

        // 三角形インデックス（1ベース）
        int[] triangles = mesh.triangles;
        for (int i = 0; i < triangles.Length; i += 3)
        {
            int a = triangles[i] + 1;
            int b = triangles[i + 1] + 1;
            int c = triangles[i + 2] + 1;
            sb.AppendLine($"f {a}/{a}/{a} {b}/{b}/{b} {c}/{c}/{c}");
        }

        return sb.ToString();
    }
}