using System.Collections.Generic;
using UnityEngine;

public class VertexTrajectory : MonoBehaviour
{
    int iVertex = 40;
    Vector3 posIni;

    LineRenderer line;
    readonly List<Vector3> points = new();

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        line = GetComponent<LineRenderer>();
        line.startWidth = 0.1f; // Set the start width of the line
        line.endWidth   = 0.1f; // Set the end width of the line
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] vertices = mesh.vertices; // Get the vertices of the mesh
        this.posIni = vertices[iVertex];
    }

    // Update is called once per frame
    void Update()
    {
        var pos = this.transform.TransformPoint(this.posIni);
        points.Add(pos);
        line.positionCount = points.Count;
        line.SetPositions(points.ToArray());    
    }
}
