using UnityEngine;


public class MyCloth : MonoBehaviour
{
    float stiffness = 60.0f;
    float mass = 1.0f;
    float learning_rate = 1.0e-3f;
    Vector3 gravity = new Vector3(0.0f, -0.1f, 0.0f);
    int[] line2vtx; // the vertex indicies of the edges
    Vector3[] vtx2xyz_ini; // initial vertex positions
    Mesh mesh; // mesh filter component

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        mesh = GetComponent<MeshFilter>().mesh; // get component
        Vector3[] vtx2xyz = mesh.vertices; // get vertex positions
        var tri2vtx = mesh.triangles;  // get triangle index
        line2vtx = TopologyOfTriMesh.Line2Vtx(tri2vtx, vtx2xyz.Length);  // get index of vertex for edges in the mesh
        vtx2xyz_ini = mesh.vertices; // initial vertex positions
    }

    // Update is called once per frame
    void Update()
    {
        Vector3[] vtx2xyz = mesh.vertices; // the current vertex positions
        Debug.Log($"Energy={total_energy(vtx2xyz)}");
        for(int i=0;i<100;++i){
          step_gradient_descent(vtx2xyz); // update vertex positions by gradient descent
        }
        mesh.vertices = vtx2xyz; // update vertex positions
        mesh.RecalculateNormals(); // update normals
    }

    float total_energy(Vector3[] vtx2xyz) {
        float eng = 0.0F;
        for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
            eng += Vector3.Dot(vtx2xyz[i_vtx], -gravity) * mass;
        }
        for(int i_line=0;i_line<line2vtx.Length/2;++i_line){            
            int i0_vtx = line2vtx[i_line*2+0]; // one of the vertex index of this line element
            int i1_vtx = line2vtx[i_line*2+1]; // another vertex index of this line element 
            float length_ini = Vector3.Distance(vtx2xyz_ini[i0_vtx], vtx2xyz_ini[i1_vtx]); // initial length of this line
            Vector3[] node2xyz = new[] { vtx2xyz[i0_vtx], vtx2xyz[i1_vtx]};
            eng += energy_spring(node2xyz, length_ini, stiffness);
        }
        return eng;
    }

    void step_gradient_descent(Vector3[] vtx2xyz) {
        Vector3[] vtx2grad = new Vector3 [vtx2xyz.Length]; // get space for gradient
        for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
            vtx2grad[i_vtx] = - gravity * mass;
        }
        for(int i_line=0;i_line<line2vtx.Length/2;++i_line){            
            int i0_vtx = line2vtx[i_line*2+0]; // one of the vertex index of this line element
            int i1_vtx = line2vtx[i_line*2+1]; // another vertex index of this line element 
            float length_ini = Vector3.Distance(vtx2xyz_ini[i0_vtx], vtx2xyz_ini[i1_vtx]); // initial length of this line
            Vector3[] node2xyz = new[] { vtx2xyz[i0_vtx], vtx2xyz[i1_vtx]};
            // float w = energy_spring(node2xyz, length_ini, stiffness);
            var grad_w = gradient_spring(node2xyz, length_ini, stiffness);
            vtx2grad[i0_vtx] += grad_w[0];
            vtx2grad[i1_vtx] += grad_w[1];                        
        }
        {
            // set fixed boundary condition 
            SphereCollider sphereCollider = GetComponent<SphereCollider>();
            Vector3 center = sphereCollider.center;
            float radius =  sphereCollider.radius;
            for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
                Vector3 pos = vtx2xyz[i_vtx];
                if( Vector3.Distance(pos, center) < radius ){
                    vtx2grad[i_vtx] = Vector3.zero;
                }
            }
        }
        // gradient descent 
        for(int i_vtx=0;i_vtx<vtx2xyz.Length;++i_vtx){
            vtx2xyz[i_vtx] -= learning_rate * vtx2grad[i_vtx];
        }
    }
    
    // this function returns the spring elastic energy
    float energy_spring(Vector3[] node2xyz, float length_ini, float stiffness) {
        float length = Vector3.Distance(node2xyz[0], node2xyz[1]); // distance between p0 and p1
        float C = length - length_ini; // the length differences.
        return 0.5f * stiffness * C * C; // Hooke's law. energy is square of length difference W=1/2*k*C*C
    }

    // this function returns the gradient of a spring elastic energy w.r.t. the spring's end position
    Vector3[] gradient_spring(Vector3[] node2xyz, float length_ini, float stiffness) {
        float length = Vector3.Distance(node2xyz[0], node2xyz[1]); // distance between p0 and p1
        // ---------- write some code below

        // Avoid division by zero when nodes are at the same position
        if (length < 1e-8f)
        {
            return new Vector3[2] { Vector3.zero, Vector3.zero };
        }

        // Unit vector from p0 to p1
        Vector3 direction = (node2xyz[1] - node2xyz[0]) / length;

        // Spring force magnitude: k * (current_length - rest_length)
        float force_magnitude = stiffness * (length - length_ini);

        // Force vector
        Vector3 force = force_magnitude * direction;

        // Gradient is the negative of force (since we want to minimize energy)
        // For p0: force points away from p1 when spring is stretched
        // For p1: force points away from p0 when spring is stretched
        Vector3 grad_p0 = -force;  // Force on p0
        Vector3 grad_p1 = force;   // Force on p1 (opposite direction)

        return new Vector3[2] { grad_p0, grad_p1 };

        // return new Vector3[2] { Vector3.zero, Vector3.zero }; // comment out this line
    }
}
