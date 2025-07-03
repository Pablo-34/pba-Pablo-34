using UnityEngine;

public class MyHangerGD : MonoBehaviour
{
    public static Vector3 gravity = new Vector3(0, -1.0f, 0);
    public static float mass = 1.0f;

    Vector3 pin = new Vector3(-1f, 1f, 0); // pin position in local coordinates

    Mesh mesh;
    Vector3 cog;

    GameObject cylinder;
    GameObject sphere;


    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {

        Application.targetFrameRate = 60; // set target frame rate
        mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] vtx2xyz = mesh.vertices;
        int[] tri = mesh.triangles;
        //
        float volume = 0;
        cog = Vector3.zero; // the center of the gravity
        // -------------------------------------------------------
        // edit the code below to compute the center of the gravity (cog)
        for (int i_tri = 0; i_tri < tri.Length / 3; i_tri++)
        {
            int i0 = tri[i_tri * 3 + 0];
            int i1 = tri[i_tri * 3 + 1];
            int i2 = tri[i_tri * 3 + 2];
            Vector3 v0 = vtx2xyz[i0];
            Vector3 v1 = vtx2xyz[i1];
            Vector3 v2 = vtx2xyz[i2];
            float volume_tet = Vector3.Dot(Vector3.Cross(v0, v1), v2) / 6.0f;
            volume += volume_tet;

            // Calculate centroid of the tetrahedron
            Vector3 centroid_tet = (v0 + v1 + v2) / 4.0f;

            cog += centroid_tet * volume_tet;
        }
        // normalize the center of gravity by the total volume
        cog /= volume;
        // end of edit
        // ------------------------
        Debug.Log("Volume: " + volume + " CoG" + cog);

        cylinder = GameObject.CreatePrimitive(PrimitiveType.Cylinder);
        SetCylinderBetween(cylinder, pin, cog, 0.05f); // add a cylinder between pin and CoG
        sphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        sphere.transform.position = pin; // add a sphere at the pin position
        sphere.transform.localScale = new Vector3(0.2f, 0.2f, 0.2f); // scale the sphere
    }

    // Update is called once per frame
    void Update()
    {
        int i_frame = Time.frameCount;
        Vector3 pin_goal = pin + new Vector3(1f - 1f * Mathf.Cos(i_frame * 0.05f), 0f, 0f);
        Debug.Log("Frame: " + i_frame + " PinDelta: " + pin_goal);

        //
        float penalty = 100.0f;
        float lr = 0.001f;
        for (int iter = 0; iter < 100; ++iter)
        {
            Vector3 transl = this.transform.position; // transformation 
            Quaternion rot = this.transform.rotation; // rotation
            Vector3 pin_def = transl + rot * pin;
            float w = 0.5f * penalty * (pin_def - pin_goal).sqrMagnitude - mass * Vector3.Dot(transl + rot * cog, gravity);
            Debug.Log("Frame" + i_frame + "Energy: " + w);
            // ----------------------------------
            // edit below to compute the gradient
            Vector3 dwdt = Vector3.zero; // differentiation of the energy w.r.t. translation
            Vector3 dwdo = Vector3.zero; // differentiation of the rotation w.r.t. rotation

            // Gradient with respect to translation
            Vector3 pin_error = pin_def - pin_goal;
            dwdt = penalty * pin_error - mass * gravity;

            // Gradient with respect to rotation (using quaternion differentiation)
            // For rotation gradient, we need to consider how rotation affects both pin_def and cog_def
            Vector3 pin_local = pin;
            Vector3 cog_local = cog;

            // Derivative of pin_def with respect to rotation
            Vector3 dpin_drot = Vector3.Cross(pin_local, Vector3.zero);

            // More accurate approach using quaternion derivatives
            Vector3 rotated_pin = rot * pin_local;
            Vector3 rotated_cog = rot * cog_local;

            // Gradient components from pin constraint
            Vector3 dw_drot_pin = penalty * Vector3.Cross(rotated_pin, pin_error);

            // Gradient components from gravity potential
            Vector3 dw_drot_gravity = -mass * Vector3.Cross(rotated_cog, gravity);

            // Total rotation gradient
            dwdo = dw_drot_pin + dw_drot_gravity;
            // end of edit
            // -----------------------------------
            this.transform.position -= lr * dwdt;
            this.transform.rotation = Quaternion.AngleAxis(-dwdo.magnitude * lr * 180.0f / Mathf.PI, Vector3.Normalize(dwdo)) * rot;
        }

        {
            Vector3 transl = this.transform.position; // transformation 
            Quaternion rot = this.transform.rotation; // rotation
            Vector3 pin_def = transl + rot * pin;
            Vector3 cog_def = transl + rot * cog;
            SetCylinderBetween(cylinder, pin_def, cog_def, 0.05f); // add a cylinder between pin and CoG
        }
        sphere.transform.position = pin_goal;

    }

    public void SetCylinderBetween(
        GameObject cyl,
        Vector3 p0,
        Vector3 p1,
        float radius = 0.05f)
    {
        cyl.transform.position = 0.5f * (p0 + p1);
        Vector3 dir = p1 - p0;
        cyl.transform.rotation = Quaternion.FromToRotation(Vector3.up, dir);
        float length = dir.magnitude;
        cyl.transform.localScale = new Vector3(radius * 2, length * 0.5f, radius * 2);
    }



}



