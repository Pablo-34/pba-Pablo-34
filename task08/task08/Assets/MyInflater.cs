using UnityEngine;
using Unity.Mathematics;
using System.IO;

public class MyInflater : MonoBehaviour
{
    public float stiffness = 1.0f;
    public float volume_trg = 10.0f;
    float volume_ini = 0.0f;
    Mesh mesh;
    Vector3[] vtx2xyz_ini;
    int[] tri2vtx;
    BlockSparseMatrix bsm;
    int[] line2vtx;
    float lambda = 0.0f;

    void Start()
    {
        mesh = GetComponent<MeshFilter>().mesh;
        vtx2xyz_ini = mesh.vertices;
        tri2vtx = mesh.triangles;        
        var vtx2vtx = TopologyOfUniformMesh.Vtx2Vtx(tri2vtx, 3, vtx2xyz_ini.Length);
        Debug.Log($"num_vtx: {vtx2xyz_ini.Length}");
        // adding a vertex that connects to all vertices in the mesh. 
        // This vertex is used to specify the volume constraint using the Lagrange multiplier method.
        var vtx2vtx_extended = TopologyOfUniformMesh.AddVertexConnectingToAllVertexForVtx2Vtx(vtx2vtx.Item1, vtx2vtx.Item2);
        bsm = new BlockSparseMatrix();
        bsm.Initialize(vtx2vtx_extended.Item1, vtx2vtx_extended.Item2);
        // compute the list of edges in the mesh
        // each edge is represented by two vertex indices
        line2vtx = TopologyOfUniformMesh.Line2Vtx(tri2vtx, 3, vtx2xyz_ini.Length);

        // compute the initial volume of the mesh
        volume_ini = 0.0f;
        for (int i_tri = 0; i_tri < tri2vtx.Length / 3; ++i_tri)
        {
            int[] node2vtx = {
                tri2vtx[i_tri*3+0],
                tri2vtx[i_tri*3+1],
                tri2vtx[i_tri*3+2] };
            float3[] node2xyz = new float3[3] {
                vtx2xyz_ini[node2vtx[0]],
                vtx2xyz_ini[node2vtx[1]],
                vtx2xyz_ini[node2vtx[2]] };
            var w_dw = volume_tri_origin_graident(node2xyz);
            volume_ini += w_dw.Item1;
        }
        Debug.Log($"volume_ini {volume_ini}");        
        volume_trg = volume_ini;
    }

    void Update()
    {
        int i_frame = Time.frameCount;
        volume_trg = volume_ini * (1.0f + 0.3f * Mathf.Sin((float)i_frame * 0.05f)); // target volume
        int num_vtx = vtx2xyz_ini.Length;
        Vector3[] vtx2xyz = mesh.vertices; // the current vertex positions
        float3[] rhs_vector = new float3[num_vtx + 1];
        // initialize the right-hand side vector as zero
        for (int i_vtx = 0; i_vtx < num_vtx + 1; ++i_vtx)
        {
            rhs_vector[i_vtx] = float3.zero;
        }
        bsm.SetZero(); // initialize the block sparse matrix to zero
        for (int i_tri = 0; i_tri < tri2vtx.Length / 3; ++i_tri)
        {
            int[] node2vtx = {
                tri2vtx[i_tri*3+0],
                tri2vtx[i_tri*3+1],
                tri2vtx[i_tri*3+2] };
            float3[] node2xyz = new float3[3] {
                vtx2xyz[node2vtx[0]],
                vtx2xyz[node2vtx[1]],
                vtx2xyz[node2vtx[2]] };
            var w_dw = volume_tri_origin_graident(node2xyz);
            float volume = w_dw.Item1;
            float3[] grad_volume = w_dw.Item2;
            // -----------------------------
            // write some code below to set values in the linear system to set constraint to specify volume
            // use the lagrange multiplier method to set the volume constraint

            // RHS contributions from volume constraint for each vertex of the triangle
            rhs_vector[node2vtx[0]] += lambda * grad_volume[0];
            rhs_vector[node2vtx[1]] += lambda * grad_volume[1];
            rhs_vector[node2vtx[2]] += lambda * grad_volume[2];

            // Subtract current volume from the lambda row's RHS.
            rhs_vector[num_vtx].x -= volume;

            // Hessian contributions from volume constraint:
            // These blocks connect vertex degrees of freedom to the lambda degree of freedom.
            bsm.AddBlockAt(node2vtx[0], num_vtx, OuterProduct(grad_volume[0], new float3(1, 0, 0)));
            bsm.AddBlockAt(node2vtx[1], num_vtx, OuterProduct(grad_volume[1], new float3(1, 0, 0)));
            bsm.AddBlockAt(node2vtx[2], num_vtx, OuterProduct(grad_volume[2], new float3(1, 0, 0)));

            // Block (num_vtx, vtx_idx): Represents d^2 L / (d_lambda dx_k) = (dC / dx_k)^T.
            bsm.AddBlockAt(num_vtx, node2vtx[0], OuterProduct(new float3(1, 0, 0), grad_volume[0]));
            bsm.AddBlockAt(num_vtx, node2vtx[1], OuterProduct(new float3(1, 0, 0), grad_volume[1]));
            bsm.AddBlockAt(num_vtx, node2vtx[2], OuterProduct(new float3(1, 0, 0), grad_volume[2]));

            // end of the edit
            // ---------------------
        }
        rhs_vector[num_vtx].x += volume_trg; // this is an important line of code!
        // setting mass-spring system as a regularizer
        float energy = 0.0f;        
        for (int i_line = 0; i_line < line2vtx.Length / 2; ++i_line)
        {
            int i_vtx0 = line2vtx[i_line * 2 + 0];
            int i_vtx1 = line2vtx[i_line * 2 + 1];
            float length_ini = math.distance(vtx2xyz_ini[i_vtx0], vtx2xyz_ini[i_vtx1]);
            float3[] node2xyz = new float3[2] {
                vtx2xyz[i_vtx0],
                vtx2xyz[i_vtx1] };
            var w_dw_ddw = spring_energy_gradient_hessian(
                node2xyz, length_ini, stiffness);
            float w = w_dw_ddw.Item1; // energy
            float3[] dw = w_dw_ddw.Item2; // gradient
            float3x3[,] ddw = w_dw_ddw.Item3; // hessian
            energy += w;
            rhs_vector[i_vtx0] += dw[0];
            rhs_vector[i_vtx1] += dw[1];
            bsm.AddBlockAt(i_vtx0, i_vtx0, ddw[0, 0]);
            bsm.AddBlockAt(i_vtx0, i_vtx1, ddw[0, 1]);
            bsm.AddBlockAt(i_vtx1, i_vtx0, ddw[1, 0]);
            bsm.AddBlockAt(i_vtx1, i_vtx1, ddw[1, 1]);
        }
        // damp the system by adding a diagonal 
        for (int i_vtx = 0; i_vtx < num_vtx + 1; ++i_vtx)
        {
            bsm.AddBlockAt(i_vtx, i_vtx, float3x3.identity * 1.0f);
        }
        Debug.Log($"#frame: {i_frame},  energy: {energy}, lambda: {lambda}");
        float3[] delta = bsm.ConjugateGradientSolver(rhs_vector, 50, 1e-7f);
        for (int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
        {
            vtx2xyz[i_vtx] -= (Vector3)delta[i_vtx];
        }
        lambda -= delta[num_vtx].x;
        mesh.vertices = vtx2xyz;
        mesh.RecalculateNormals();
        // output to an OBJ file
        if (Time.frameCount == 160)
        {
            string objText = ObjExporterToAssets.MeshToObj(mesh, transform);
            string path = Application.dataPath + "/task08.obj";
            File.WriteAllText(path, objText);
            Debug.Log("Obj file written to " + path);
        }        
//        Debug.Assert(false);
    }

    static float3x3 OuterProduct(float3 a, float3 b) {
        return new float3x3(
            a[0] * b[0], a[0] * b[1], a[0] * b[2],
            a[1] * b[0], a[1] * b[1], a[1] * b[2],
            a[2] * b[0], a[2] * b[1], a[2] * b[2]);            
    }

    //  compute a spring's elastic potential energy, energy's gradient, and energy's hessian
    static (float, float3[], float3x3[,]) spring_energy_gradient_hessian(float3[] node2xyz, float length_ini, float stiffness)
    {
        float length = math.distance(node2xyz[0], node2xyz[1]); // distance between p0 and p1
        float C = length - length_ini; // the length differences.
        float3 u01 = math.normalize(node2xyz[1] - node2xyz[0]);
        float energy = 0.5f * stiffness * C * C; // Hooke's law. energy is square of length difference W=1/2*k*C*C   
        float3[] dC = new float3[] { -u01, u01 };
        float3[] gradients = new float3[2] {
            stiffness * C * dC[0],
            stiffness * C * dC[1]};
        float3x3 l = OuterProduct(u01, u01);
        float3x3 n = stiffness * l;
        float3x3 m = (stiffness * C / length) * (float3x3.identity - l);
        float3x3 o = n + m;
        float3x3[,] hessian = new float3x3[2, 2] {
            {o, -o},
            {-o, o}
        };
        return (energy, gradients, hessian);
    }

    // compute volume of a tetrahedron connecting vertices of the triangle and the origin
    // and the gradient of the volume with respect to the vertex positions.
    static (float, float3[]) volume_tri_origin_graident(float3[] node2xyz)
    {
        float volume = 1.0f / 6.0f * math.dot(math.cross(node2xyz[1], node2xyz[2]), node2xyz[0]);   
        float3[] gradient = new float3[3];
        // ---------        
        // write some code to compute the gradient of the volume with respect to the vertex positions
        gradient[0] = 1.0f / 6.0f * math.cross(node2xyz[1], node2xyz[2]);
        gradient[1] = 1.0f / 6.0f * math.cross(node2xyz[2], node2xyz[0]);
        gradient[2] = 1.0f / 6.0f * math.cross(node2xyz[0], node2xyz[1]);

        // end of editing
        // ---------------
        return (volume, gradient);
    }    
}
