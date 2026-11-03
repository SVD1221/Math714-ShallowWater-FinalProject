#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

void writeTo(ofstream *file, vector<float> *vec)
{
    for (int i = 0; i < (*vec).size(); i++){
        *file << (*vec)[i] << " ";
    }
    *file << "\n";
}

// g++ lax_friedrich.cpp -o lax_friedrich

// Shallow Water equations
// h_t + (hu)_x = 0
// (hu)_t + (huu + 0.5ghh) = 0
// ICs: h(x, 0) = h_0 & u(x, 0) = u_0   
// BCs: Periodic => h(0, t) = h(L, t) and h_x(0, t) = h_x(L, t)
// New BCs under d = uh & p = huu+0.5ghh:
//          => 
//          => 
// Leapfrog in time
// Center Difference in space

const int N = 20;
const int interations = 1000;
const float g = 9.8;

float speed(vector<float> q)
{
    return abs(q[1]/q[0]) + sqrt(g*q[0]);
}

void solve(vector<float> h0, vector<float> u0, float delta_t, float delta_x, int iterations)
{
    vector<float> mass_flux(N);
    vector<float> momentum_flux(N);

    for (int j = 0; j < N; j++){
        mass_flux[j] = h0[j]*u0[j];
        momentum_flux[j] = h0[j]*u0[j]*u0[j] + 0.5*g*h0[j]*h0[j];
    }


    float nu = delta_t/delta_x;

    vector<vector<float>> curr_q = {h0, mass_flux};
    vector<vector<float>> next_q(2, vector<float>(N, 0));
    vector<vector<float>> curr_F = {mass_flux, momentum_flux};
    vector<vector<float>> next_F(2, vector<float>(N, 0));

    for (int n = 0; n < iterations; n++){
        ofstream file;
        file.open("output/lax_friedrich" + std::to_string(n) + ".txt");

        vector<float> h(N, 0);
        vector<float> u(N, 0);
        for (int j = 0; j < N; j++){
            h[j] = curr_q[0][j];
            u[j] = curr_q[1][j]/curr_q[0][j];
        }

        writeTo(&file, &h);
        writeTo(&file, &u);
        file.close();

        vector<vector<float>> F_pos_half(2, vector<float>(N, 0));
        vector<vector<float>> F_neg_half(2, vector<float>(N, 0));

        for (int j = 0; j < N; j++){

            float s_pos_half = max(speed({curr_q[0][(j+1)%N], curr_q[1][(j+1)%N]}), speed({curr_q[0][j], curr_q[1][j]}));
            float u_pos_diff = curr_q[1][(j+1)%N]/curr_q[0][(j+1)%N] - curr_q[1][j]/curr_q[0][j];

            F_pos_half[0][j] = 0.5*(curr_F[0][(j+1)%N] - curr_F[0][j]) - 0.5*s_pos_half*u_pos_diff;
            F_pos_half[1][j] = 0.5*(curr_F[1][(j+1)%N] - curr_F[1][j]) - 0.5*s_pos_half*u_pos_diff;

            float s_neg_half = max(speed({curr_q[0][j], curr_q[1][j]}), speed({curr_q[0][(j-1+N)%N], curr_q[1][(j-1+N)%N]}));
            float u_neg_diff = curr_q[1][j]/curr_q[0][j] - curr_q[1][(j-1+N)%N]/curr_q[0][(j-1+N)%N];
            
            F_neg_half[0][j] = 0.5*(curr_F[0][j] - curr_F[0][(j-1+N)%N]) - 0.5*s_neg_half*u_neg_diff;
            F_neg_half[1][j] = 0.5*(curr_F[1][j] - curr_F[1][(j-1+N)%N]) - 0.5*s_neg_half*u_neg_diff;
        }

        for (int j = 0; j < N; j++){
            next_q[0][j] = curr_q[0][j] - nu*(F_pos_half[0][j] - F_neg_half[0][j]);
            next_q[1][j] = curr_q[1][j] - nu*(F_pos_half[1][j] - F_neg_half[1][j]);
        }

        for (int j = 0; j < N; j++){
            next_F[0][j] = next_q[1][j];
            next_F[1][j] = next_q[1][j]*next_q[1][j]/next_q[0][j] + 0.5*g*next_q[0][j]*next_q[0][j];
        }

        for (int j = 0; j < N; j++){
            curr_q[0][j] = next_q[0][j];
            curr_q[1][j] = next_q[1][j];

            curr_F[0][j] = next_F[0][j];
            curr_F[1][j] = next_F[1][j];
        }
    }

    ofstream file;
    file.open("output/lax_friedrich" + std::to_string(iterations) + ".txt");

    vector<float> h(N, 0);
    vector<float> u(N, 0);
    for (int j = 0; j < N; j++){
        h[j] = curr_q[0][j];
        u[j] = curr_q[1][j]/curr_q[0][j];
    }

    writeTo(&file, &h);
    writeTo(&file, &u);
    file.close();
}


int main()
{
    vector<float> h0(N);
    vector<float> u0(N);

    float delta_t = 0.001;
    float delta_x = 1.0/(N-1);

    for (int j = 0; j < N; j++){
        h0[j] = j*(N-1-j)/((N-1)*(N-1)/2/2.) + 1;
        u0[j] = 0;
    }

    solve(h0, u0, delta_t, delta_x, interations);

    return 0;
}