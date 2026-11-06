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
// (hu)_t + (huu + 0.5ghh)_x = 0
// ICs: h(x, 0) = h_0 & u(x, 0) = u_0   
// BCs: periodic

const float g = 9.8;

float speed(vector<float> q)
{
    return abs(q[1]/q[0]) + sqrt(g*q[0]);
}

void solve(vector<float> h0, vector<float> u0, float delta_x, float delta_t, int num_point, int iterations)
{
    const float nu = delta_t/delta_x;

    vector<float> mass_flux(num_point);
    vector<float> momentum_flux(num_point);

    for (int j = 0; j < num_point; j++){
        mass_flux[j] = h0[j]*u0[j];
        momentum_flux[j] = h0[j]*u0[j]*u0[j] + 0.5*g*h0[j]*h0[j];
    }

    vector<vector<float>> curr_q = {h0, mass_flux};
    vector<vector<float>> next_q(2, vector<float>(num_point, 0));
    vector<vector<float>> curr_F = {mass_flux, momentum_flux};
    vector<vector<float>> next_F(2, vector<float>(num_point, 0));

    for (int n = 0; n < iterations; n++){
        ofstream file;
        file.open("output/lax_friedrich" + std::to_string(n) + ".txt");

        vector<float> h(num_point, 0);
        vector<float> u(num_point, 0);
        for (int j = 0; j < num_point; j++){
            h[j] = curr_q[0][j];
            u[j] = curr_q[1][j]/curr_q[0][j];
        }

        writeTo(&file, &h);
        writeTo(&file, &u);
        file.close();

        //IMEX
        // implicitly solve: h_t + (hu)_x = 0
        // explicitly solve: (hu)_t + (huu + 0.5ghh)_x = 0


        //lax-friedrich
        vector<vector<float>> F_pos_half(2, vector<float>(num_point, 0));
        vector<vector<float>> F_neg_half(2, vector<float>(num_point, 0));

        for (int j = 0; j < num_point; j++) {

            // if (j == 0 || j == num_point-1) {
            //     float s_pos_half = max(speed({curr_q[0][(j+1)%num_point], curr_q[1][(j+1)%num_point]}), speed({curr_q[0][j], 0}));
            //     float u_pos_diff = curr_q[1][(j+1)%num_point]/curr_q[0][(j+1)%num_point];

            //     F_pos_half[0][j] = 0.5*(curr_F[0][(j+1)%num_point] - curr_F[0][j]) - 0.5*s_pos_half*u_pos_diff;
            //     F_pos_half[1][j] = 0.5*(curr_F[1][(j+1)%num_point] - curr_F[1][j]) - 0.5*s_pos_half*u_pos_diff;

            //     float s_neg_half = max(speed({curr_q[0][j], 0}), speed({curr_q[0][(j-1+num_point)%num_point], curr_q[1][(j-1+num_point)%num_point]}));
            //     float u_neg_diff = -curr_q[1][(j-1+num_point)%num_point]/curr_q[0][(j-1+num_point)%num_point];
                
            //     F_neg_half[0][j] = 0.5*(curr_F[0][j] - curr_F[0][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
            //     F_neg_half[1][j] = 0.5*(curr_F[1][j] - curr_F[1][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
            // }
            // else{
            //     float s_pos_half = max(speed({curr_q[0][(j+1)%num_point], curr_q[1][(j+1)%num_point]}), speed({curr_q[0][j], curr_q[1][j]}));
            //     float u_pos_diff = curr_q[1][(j+1)%num_point]/curr_q[0][(j+1)%num_point] - curr_q[1][j]/curr_q[0][j];

            //     F_pos_half[0][j] = 0.5*(curr_F[0][(j+1)%num_point] - curr_F[0][j]) - 0.5*s_pos_half*u_pos_diff;
            //     F_pos_half[1][j] = 0.5*(curr_F[1][(j+1)%num_point] - curr_F[1][j]) - 0.5*s_pos_half*u_pos_diff;

            //     float s_neg_half = max(speed({curr_q[0][j], curr_q[1][j]}), speed({curr_q[0][(j-1+num_point)%num_point], curr_q[1][(j-1+num_point)%num_point]}));
            //     float u_neg_diff = curr_q[1][j]/curr_q[0][j] - curr_q[1][(j-1+num_point)%num_point]/curr_q[0][(j-1+num_point)%num_point];
                
            //     F_neg_half[0][j] = 0.5*(curr_F[0][j] - curr_F[0][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
            //     F_neg_half[1][j] = 0.5*(curr_F[1][j] - curr_F[1][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
            // }

            float s_pos_half = max(speed({curr_q[0][(j+1)%num_point], curr_q[1][(j+1)%num_point]}), speed({curr_q[0][j], curr_q[1][j]}));
            float u_pos_diff = curr_q[1][(j+1)%num_point]/curr_q[0][(j+1)%num_point] - curr_q[1][j]/curr_q[0][j];

            F_pos_half[0][j] = 0.5*(curr_F[0][(j+1)%num_point] - curr_F[0][j]) - 0.5*s_pos_half*u_pos_diff;
            F_pos_half[1][j] = 0.5*(curr_F[1][(j+1)%num_point] - curr_F[1][j]) - 0.5*s_pos_half*u_pos_diff;

            float s_neg_half = max(speed({curr_q[0][j], curr_q[1][j]}), speed({curr_q[0][(j-1+num_point)%num_point], curr_q[1][(j-1+num_point)%num_point]}));
            float u_neg_diff = curr_q[1][j]/curr_q[0][j] - curr_q[1][(j-1+num_point)%num_point]/curr_q[0][(j-1+num_point)%num_point];
            
            F_neg_half[0][j] = 0.5*(curr_F[0][j] - curr_F[0][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
            F_neg_half[1][j] = 0.5*(curr_F[1][j] - curr_F[1][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
        }

        for (int j = 0; j < num_point; j++) {
            // if (j == 0 || j == num_point-1) {
            //     next_q[0][j] = curr_q[0][j] - nu*(F_pos_half[0][j] - F_neg_half[0][j]);
            //     next_q[1][j] = -nu*(F_pos_half[1][j] - F_neg_half[1][j]);
            // }
            // else{
            //     next_q[0][j] = curr_q[0][j] - nu*(F_pos_half[0][j] - F_neg_half[0][j]);
            //     next_q[1][j] = curr_q[1][j] - nu*(F_pos_half[1][j] - F_neg_half[1][j]);
            // }

            next_q[0][j] = curr_q[0][j] - nu*(F_pos_half[0][j] - F_neg_half[0][j]);
            next_q[1][j] = curr_q[1][j] - nu*(F_pos_half[1][j] - F_neg_half[1][j]);
        }

        for (int j = 0; j < num_point; j++) {
            // if (j == 0 || j == num_point-1) {
            //     next_F[0][j] = 0;
            //     next_F[1][j] = 0.5*g*next_q[0][j]*next_q[0][j];
            // }
            // else{
            //     next_F[0][j] = next_q[1][j];
            //     next_F[1][j] = next_q[1][j]*next_q[1][j]/next_q[0][j] + 0.5*g*next_q[0][j]*next_q[0][j];
            // }

            next_F[0][j] = next_q[1][j];
            next_F[1][j] = next_q[1][j]*next_q[1][j]/next_q[0][j] + 0.5*g*next_q[0][j]*next_q[0][j];
        }

        // Update step
        for (int j = 0; j < num_point; j++) {
            curr_q[0][j] = next_q[0][j];
            curr_q[1][j] = next_q[1][j];

            curr_F[0][j] = next_F[0][j];
            curr_F[1][j] = next_F[1][j];
        }
    }

    ofstream file;
    file.open("output/lax_friedrich" + std::to_string(iterations) + ".txt");

    vector<float> h(num_point, 0);
    vector<float> u(num_point, 0);
    for (int j = 0; j < num_point; j++){
        h[j] = curr_q[0][j];
        u[j] = curr_q[1][j]/curr_q[0][j];
    }

    writeTo(&file, &h);
    writeTo(&file, &u);
    file.close();
}


int main()
{
    const int N = 50;
    const int T = 5000;
    const float delta_x = 1./(N-1);
    const float delta_t = 0.00075;

    vector<float> h0(N);
    vector<float> u0(N);

    for (int j = 0; j < N; j++){
        float x = j*delta_x;
        // h0[j] = 4.0*x*(1-x) + 1;            // symmetric IC
        h0[j] = 6.75*x*x*(1-x)/2 + 1;       // asymmetric IC

        u0[j] = 0;
    }

    solve(h0, u0, delta_x, delta_t, N, T);

    return 0;
}