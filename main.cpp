#include <bits/stdc++.h>
#define loop(i, start, end) for(int i = start; i <= end; ++i)
#define rloop(i,end, start) for(int i = end; i >= start ; --i)
using namespace std;
using db = double;
template<typename T>
using V = vector<T>;
const db EPS = 1e-9;
V<V<db>> transpose(V<V<db>>& A){
    int n = A.size(), m = A[0].size();
    V<V<db>> At(m, V<db>(n, 0.0));
    loop(i, 0, n - 1) loop(j, 0, m - 1) At[j][i] = A[i][j];
    return At;
}
V<V<db>> matmul(V<V<db>>& A, V<V<db>>& B){
    int n = A.size(), m = A[0].size(), r = B[0].size();
    V<V<db>> C(n, V<db>(r, 0.0));
    loop(i, 0, n - 1){
        loop(j, 0, r - 1){
            loop(k, 0, m - 1){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}
void subtract_rows(V<V<db>>& A, int i, int j, db p, db q){
    int n = A.size(), m = A[0].size();
    loop(col, 0, m - 1){
        A[i][col] = p * A[i][col] - q * A[j][col];
    }
}
void RREF(V<V<db>>& A){//pivots are not necessarily 1
    int n = A.size(), m = A[0].size();
    int curr_row = 0, pivot_col = 0;
    while(curr_row < n && pivot_col < m){
        int pivot_row = curr_row;
        loop(i, curr_row, n - 1){
            if(abs(A[i][pivot_col]) > abs(A[curr_row][pivot_col])){
                pivot_row = i;
                swap(A[curr_row], A[pivot_row]);
            }
        }
        pivot_row = curr_row;
        if(abs(A[pivot_row][pivot_col]) > EPS /*means non zero*/){
            db x = A[pivot_row][pivot_col];
            loop(j, 0, m - 1) A[pivot_row][j] /= x;
            loop(r, 0, n - 1){
                if(r == pivot_row) continue;
                db m = A[r][pivot_col] / A[pivot_row][pivot_col];
                subtract_rows(A, r, pivot_row, 1, m);
            }
            curr_row++;
        }
        pivot_col++;
    }
}
int matrix_rank(V<V<db>> B){
    int n = B.size(), m = B[0].size();
    V<V<db>> A(n, V<db>(m, 0.0));
    loop(i, 0, n - 1) loop(j, 0, m - 1) A[i][j] = B[i][j];
    RREF(A);
    int rk = 0;
    loop(i, 0, n - 1){
        loop(j, 0, m - 1){
            if(A[i][j] != 0){
                rk++;
                break;
            }
        }
    }
    return rk;
}
V<V<db>> system_solve(V<V<db>>& A, V<V<db>>& b){
    int n = A.size(), m = A[0].size();
    V<V<db>> x(m, V<db>(1, 0.0));
    auto At = transpose(A);
    auto M = matmul(At, A);
    auto c = matmul(At, b);
    loop(i, 0, m - 1){
        M[i].emplace_back(c[i][0]);
    }
    RREF(M);
    rloop(i, m - 1, 0){
        x[i][0] = M[i][m];
    }
    return x;
}
V<V<db>> work(V<V<db>>& T, V<V<db>>& α){
    int n = α.size();
    V<V<db>> B(n, V<db>(1, 0.0)); //B is a basis for Z(α, T)
    loop(i, 0, n - 1){
        B[i][0] = α[i][0];
    }
    //linearly independent if rank < number of columns
    int cols = 1;
    if(matrix_rank(B) < cols){
        return {{0}};//Z(α, T) = 0
    }
    while(matrix_rank(B) == cols){
        V<V<db>> prev(n, V<db>(1, 0.0));
        loop(i, 0, n - 1){
            prev[i][0] = B[i][cols - 1];
        }
        auto new_col = matmul(T, prev);
        loop(i, 0, n - 1){
            B[i].emplace_back(new_col[i][0]);
        }
        cols++;
    }
    cols--;
    for(auto& rows: B) rows.pop_back();
    V<V<db>> prev(n, V<db>(1, 0.0));
    loop(i, 0, n - 1){
        prev[i][0] = B[i][cols - 1];
    }
    auto Tkα = matmul(T, prev);
    auto coeff = system_solve(B, Tkα);
    reverse(coeff.begin(), coeff.end());
    return coeff;
}
V<V<db>> companion_matrix(V<V<db>>& T, V<V<db>>& α){
    V<V<db>> p = work(T, α);
    int k = p.size();
    V<V<db>> companion(k, V<db>(k, 0.0));
    loop(j, 0, k - 2){
        companion[j + 1][j] = 1;
    }
    loop(i, 0, k - 1){
        companion[i][k - 1] = -p[k - 1 - i][0];
    }
    return companion;
}
int main(){
    V<V<db>> T = {{0, 0, 0}, {1, 0, -1}, {0, 1, 2}};
    V<V<db>> α = {{1}, {0}, {0}};

    V<V<db>> p = work(T, α);
    reverse(p.begin(), p.end());
    int k = p.size() - 1;
    cout << "x" << k + 1 << ": " << "1\n";
    for(auto& coeff: p){
        cout << "x" << k << ": " << coeff[0] << '\n';
    }
    k = p.size();
    auto cm = companion_matrix(T, α);
    cm[0][k - 1] = 1;
    loop(i, 1, k - 1){
        loop(j, 0, k - 1) cout << cm[i][j] << ' ';
        cout << '\n';
    }
    return 0;
}