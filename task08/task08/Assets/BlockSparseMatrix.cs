using Unity.Mathematics;
using System.Diagnostics;


public class BlockSparseMatrix
{
    public int[] row2idx;
    public int[] idx2col;
    public float3x3[] idx2val;
    float3[] p;
    float3[] Ap;

    public void Initialize(int[] _row2col, int[] _idx2col) {
        row2idx = _row2col;
        idx2col = _idx2col;
        idx2val = new float3x3 [idx2col.Length];
        int n = row2idx.Length - 1;
        p = new float3[n];
        Ap = new float3[n];
    }

    public void SetZero() {
        for(int idx=0;idx<idx2val.Length;++idx){
            idx2val[idx] = float3x3.zero;
        }
    }

    public void AddBlockAt(int i_row, int i_col, float3x3 val) {
        for(int idx=row2idx[i_row];idx<row2idx[i_row+1];++idx){
            if( idx2col[idx] == i_col ){
                idx2val[idx] += val;
                return;
            }
        }
        UnityEngine.Debug.Assert(false);
    }

    public void SetFixed(int i_vtx) {
        for(int j_vtx=0; j_vtx<row2idx.Length-1; ++j_vtx){
            if( j_vtx == i_vtx ){
                for(int idx=row2idx[j_vtx];idx<row2idx[j_vtx+1];++idx){
                    if( idx2col[idx] == i_vtx ){                
                        idx2val[idx] += float3x3.identity;
                    }　else{
                        idx2val[idx] = float3x3.zero;
                        // Debug.Log($"{j_vtx} {idx2col[idx]}, {idx2val[idx]}");
                    }
                }
            }
            else {
                for(int idx=row2idx[j_vtx];idx<row2idx[j_vtx+1];++idx){
                    if( idx2col[idx] == i_vtx ){                
                        idx2val[idx] = float3x3.zero;
                    }
                }
            }
        }
    }

    public void Multiply(float3[] x, float3[] result) {    
        // Multiply the sparse matrix with vector x and store the result in 'result'
        for (int i_row = 0; i_row < row2idx.Length - 1; i_row++) {
            result[i_row] = float3.zero;
            for (int idx = row2idx[i_row]; idx < row2idx[i_row + 1]; idx++) {
                int j_col = idx2col[idx];
                result[i_row] += math.mul(idx2val[idx], x[j_col]);
            }
        }
    }    

    public float3[] ConjugateGradientSolver(float3[] r, int maxIterations, float tolerance) {
        int n = r.Length;
        float3[] x = new float3[n];
        for(int i=0;i<n;++i){
            x[i] = float3.zero;
        }

        for (int i = 0; i < n; i++) {
            p[i] = r[i];
        }

        float rsOld = 0;
        for (int i = 0; i < n; i++) {
//            UnityEngine.Debug.Log($"{i}, {r[i]}");
            rsOld += math.dot(r[i], r[i]);
        }

        for (int iter = 0; iter < maxIterations; iter++) {            
            Multiply(p, Ap);

            float pAp = 0;
            for (int i = 0; i < n; i++) {
                // UnityEngine.Debug.Log($"{i}, {p[i]}, {Ap[i]}");
                pAp += math.dot(p[i], Ap[i]);
            }
            // UnityEngine.Debug.Assert(false);
            float alpha = rsOld / pAp;

            for (int i = 0; i < n; i++) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }

            float rsNew = 0;
            for (int i = 0; i < n; i++) {
                rsNew += math.dot(r[i], r[i]);
            }

            // UnityEngine.Debug.Log($"{iter}, {rsNew} {rsOld} {alpha}, {pAp}");
            if (math.sqrt(rsNew) < tolerance) {
                break;
            }

            float beta = rsNew / rsOld;
            for (int i = 0; i < n; i++) {
                p[i] = r[i] + beta * p[i];
            }

            rsOld = rsNew;
        }
        return x;
    }
}
