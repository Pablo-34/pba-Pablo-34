using UnityEngine;
using Unity.Mathematics;
using System.IO;
using System.Text;

public class MyCloth : MonoBehaviour
{
    public float stiffness = 100.0f;
    public float mass_point = 1.0f;
    public Vector3 gravity = new Vector3(0.0f, -0.05f, 0.0f);
    public float timeStep = 0.10f; // Time step for simulation
    int[] line2vtx; // the vertex indicies of the edges
    Vector3[] vtx2xyz_ini; // initial vertex positions
    Vector3[] vtx2velo; // the current velocity of the vertices
    BlockSparseMatrix bsm;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        vtx2xyz_ini = mesh.vertices;        
        vtx2velo = new Vector3[vtx2xyz_ini.Length];        
        for(int i_vtx=0;i_vtx<vtx2xyz_ini.Length;++i_vtx){
            vtx2velo[i_vtx] = Vector3.zero;
        }
        //
        var tri2vtx = mesh.triangles;    
        line2vtx = TopologyOfUniformMesh.Line2Vtx(tri2vtx, 3, vtx2xyz_ini.Length); 
        //        
        {
            // initialize the block sparse matrix pattern
            var vtx2vtx = TopologyOfUniformMesh.Vtx2Vtx(tri2vtx, 3, vtx2xyz_ini.Length);             
            bsm = new BlockSparseMatrix();
            bsm.Initialize(vtx2vtx.Item1, vtx2vtx.Item2);
        }
    }

    // Update is called once per frame
    void Update()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] vtx2xyz = mesh.vertices; // the current vertex positions
        float energy = 0.0f;
        float3[] vtx2grad = new float3 [vtx2xyz.Length];        
        for (int i_vtx = 0; i_vtx < vtx2xyz.Length; ++i_vtx) {
            vtx2grad[i_vtx] = float3.zero;
        }
        bsm.SetZero();        
        // Update positions and velocities using semi-implicit Euler integration
        for (int i_vtx = 0; i_vtx < vtx2xyz.Length; ++i_vtx)
        {
            vtx2xyz[i_vtx] += (Vector3)(float3)vtx2velo[i_vtx] * timeStep;
        } 
        for(int i_line=0;i_line<line2vtx.Length/2;++i_line){            
            int[] node2vtx = new int [] {line2vtx[i_line*2+0], line2vtx[i_line*2+1]}; // another vertex index of this line element 
            float length_ini = Vector3.Distance(vtx2xyz_ini[node2vtx[0]], vtx2xyz_ini[node2vtx[1]]); // initial length of this line
            float3[] node2xyz = new[] { (float3)vtx2xyz[node2vtx[0]], (float3)vtx2xyz[node2vtx[1]]};
            float w = energy_spring(node2xyz, length_ini, stiffness);
            energy += w;
            var grad_w = gradient_spring(node2xyz, length_ini, stiffness);
            var hessian_w = hessian_spring(node2xyz, length_ini, stiffness);
            for(int i_node=0;i_node<2;++i_node){
                vtx2grad[node2vtx[i_node]] += grad_w[i_node]; 
                for(int j_node=0;j_node<2;++j_node){           
                   bsm.AddBlockAt(node2vtx[i_node], node2vtx[j_node], hessian_w[i_node,j_node]);
                }                
            }
        }
        for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
            float3x3 mass_matrix = float3x3.identity * (mass_point / (timeStep * timeStep));
            bsm.AddBlockAt(i_vtx, i_vtx, mass_matrix);
        }
        for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
            vtx2grad[i_vtx] -= (float3)gravity * mass_point;
            energy -= mass_point * Vector3.Dot(gravity, vtx2xyz[i_vtx]);
        }
        {
            // set fixed boundary condition
            SphereCollider sphereCollider = GetComponent<SphereCollider>();
            Vector3 center = sphereCollider.center;
            float radius =  sphereCollider.radius;
            for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
                Vector3 pos = vtx2xyz_ini[i_vtx];
                if( Vector3.Distance(pos, center) < radius ){ // inside sphere
                    vtx2grad[i_vtx] = float3.zero;
                    vtx2velo[i_vtx] = Vector3.zero;
                    bsm.AddBlockAt(i_vtx, i_vtx, 1.0f);
                    bsm.SetFixed(i_vtx);
                }
            }
        }
        // solve linear system
        float3[] delta_vtx2xyz = bsm.ConjugateGradientSolver(vtx2grad, 10, 1.0e-6f);
        for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
            vtx2velo[i_vtx] -= (Vector3)delta_vtx2xyz[i_vtx]/timeStep;
            vtx2xyz[i_vtx] -= (Vector3)(float3)delta_vtx2xyz[i_vtx];
        } 
        Debug.Log($"energy {energy}");
        mesh.vertices = vtx2xyz;
        mesh.RecalculateNormals();
        if (Time.frameCount == 500) {
            string objText = ObjExporterToAssets.MeshToObj(mesh, transform);
            string path = Application.dataPath + "/MyCloth.obj";
            File.WriteAllText(path, objText);
            Debug.Log("Obj file written to " + path);
        }
    }

    static float energy_spring(float3[] node2xyz, float length_ini, float stiffness) {
        float length = math.distance(node2xyz[0], node2xyz[1]); // distance between p0 and p1
        float C = length - length_ini; // the length differences.
        return 0.5f * stiffness * C * C; // Hooke's law. energy is square of length difference W=1/2*k*C*C
    }

    // compute the gradient of a spring elastic energy w.r.t their end positions
    static float3[] gradient_spring(float3[] node2xyz, float length_ini, float stiffness) {
        float length = math.distance(node2xyz[0], node2xyz[1]); // distance between p0 and p1
        float C = length - length_ini; // the length differences.
        float3 u01 = math.normalize(node2xyz[1] - node2xyz[0]);
        float3[] dC = new float3[] {-u01, u01};
        return new float3[2] {
            stiffness * C * dC[0], 
            stiffness * C * dC[1]};
    }

    static float3x3 OuterProduct(float3 a, float3 b) {
        return new float3x3(
            a[0] * b[0], a[0] * b[1], a[0] * b[2],
            a[1] * b[0], a[1] * b[1], a[1] * b[2],
            a[2] * b[0], a[2] * b[1], a[2] * b[2]);            
    }

    // compute the hessian of a spring elastic energy w.r.t their end positions
    static float3x3[,] hessian_spring(float3[] node2xyz, float length_ini, float stiffness) {
        float length = math.distance(node2xyz[0], node2xyz[1]); // distance between p0 and p1
        float C = length - length_ini; // the length differences.
        float3 u01 = math.normalize(node2xyz[1] - node2xyz[0]);
        float3x3 l = OuterProduct(u01, u01);
        float3x3 n = stiffness * l;
        // ----------------------
        // write some code below to modify the hessian

        // Add the second-order term for accurate Hessian computation
        float3x3 identity = float3x3.identity;
        float3x3 second_order = stiffness * (C / length) * (identity - l);

        // Complete Hessian: first-order + second-order terms
        float3x3 o = n + second_order;

        return new float3x3[2,2] {
            {o, -o}, 
            {-o, o}                   
        };
    }
}
