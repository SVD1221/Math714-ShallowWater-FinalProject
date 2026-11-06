#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

// g++ bathymetry.cpp -o bathymetry

void writeTo(ofstream *file, vector<float> *vec)
{
    for (int i = 0; i < (*vec).size(); i++){
        *file << (*vec)[i] << " ";
    }
    *file << "\n";
}

// Shallow Water equations with Bathymetry
// w_t + (wu)_x = 0
// (wu)_t + (wuu + 0.5gww)_x = 0
//  w(x,t) = h(x,t) + b(x), h(x,t) > 0
// ICs: w(x, 0) = w_0 & u(x, 0) = u_0   
// BCs: periodic

const float g = 9.8;
const float pi = 3.1415;

float speed(vector<float> q, float b)
{
    // wu/(h+b) + sqrt(g(h+b))
    return abs(q[1]/(q[0]+b)) + sqrt(g*(q[0]+b));
}

void solve(vector<float> b, vector<float> h0, vector<float> u0, float delta_x, float delta_t, int num_point, int iterations)
{
    //  w(x,t) = h(x,t) + b(x)

    const float nu = delta_t/delta_x;

    vector<float> mass_flux(num_point);
    vector<float> momentum_flux(num_point);
    vector<float> w(num_point);

    for (int j = 0; j < num_point; j++){
        w[j] = h0[j] + b[j];                                      //good
        mass_flux[j] = w[j]*u0[j];                              //good
        momentum_flux[j] = w[j]*u0[j]*u0[j] + 0.5*g*w[j]*w[j];  //good
    }

    vector<vector<float>> curr_q = {h0, mass_flux};             //good
    vector<vector<float>> next_q(2, vector<float>(num_point, 0));
    vector<vector<float>> curr_F = {mass_flux, momentum_flux};
    vector<vector<float>> next_F(2, vector<float>(num_point, 0));

    for (int n = 0; n < iterations; n++){
        ofstream file;
        file.open("output/bathymetry/lax_friedrich" + std::to_string(n) + ".txt");

        vector<float> h(num_point, 0);
        vector<float> u(num_point, 0);
        for (int j = 0; j < num_point; j++){
            h[j] = curr_q[0][j];
            u[j] = curr_q[1][j]/(curr_q[0][j]+b[j]);
        }

        writeTo(&file, &h);
        writeTo(&file, &u);
        file.close();

        //lax-friedrich
        vector<vector<float>> F_pos_half(2, vector<float>(num_point, 0));
        vector<vector<float>> F_neg_half(2, vector<float>(num_point, 0));

        for (int j = 0; j < num_point; j++) {
            float s_pos_half = max(speed({curr_q[0][(j+1)%num_point], curr_q[1][(j+1)%num_point]}, b[(j+1)%num_point]), 
                                    speed({curr_q[0][j], curr_q[1][j]}, b[j]));
            float u_pos_diff = curr_q[1][(j+1)%num_point]/w[(j+1)%num_point] - curr_q[1][j]/w[j];

            F_pos_half[0][j] = 0.5*(curr_F[0][(j+1)%num_point] - curr_F[0][j]) - 0.5*s_pos_half*u_pos_diff;
            F_pos_half[1][j] = 0.5*(curr_F[1][(j+1)%num_point] - curr_F[1][j]) - 0.5*s_pos_half*u_pos_diff;

            float s_neg_half = max(speed({curr_q[0][j], curr_q[1][j]}, b[j]), 
                                    speed({curr_q[0][(j-1+num_point)%num_point], curr_q[1][(j-1+num_point)%num_point]}, b[(j-1+num_point)%num_point]));
            float u_neg_diff = curr_q[1][j]/w[j] - curr_q[1][(j-1+num_point)%num_point]/w[(j-1+num_point)%num_point];
            
            F_neg_half[0][j] = 0.5*(curr_F[0][j] - curr_F[0][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
            F_neg_half[1][j] = 0.5*(curr_F[1][j] - curr_F[1][(j-1+num_point)%num_point]) - 0.5*s_neg_half*u_neg_diff;
        }

        for (int j = 0; j < num_point; j++) {
            next_q[0][j] = curr_q[0][j] - nu*(F_pos_half[0][j] - F_neg_half[0][j]);     // updates h
            next_q[1][j] = curr_q[1][j] - nu*(F_pos_half[1][j] - F_neg_half[1][j]);     // updates wu
        }

        for (int j = 0; j < num_point; j++) {
            next_F[0][j] = next_q[1][j];
            next_F[1][j] = next_q[1][j]*next_q[1][j]/w[j] + 0.5*g*(next_q[0][j]+b[j])*w[j];
        }

        // Update step
        for (int j = 0; j < num_point; j++) {
            curr_q[0][j] = next_q[0][j];
            curr_q[1][j] = next_q[1][j];

            curr_F[0][j] = next_F[0][j];
            curr_F[1][j] = next_F[1][j];

            w[j] = next_q[0][j] + b[j];
        }
    }

    ofstream file;
    file.open("output/bathymetry/lax_friedrich" + std::to_string(iterations) + ".txt");

    vector<float> u(num_point, 0);
    for (int j = 0; j < num_point; j++){
        w[j] = curr_q[0][j] + b[j];                 //calc w = h+b
        u[j] = curr_q[1][j]/w[j];
    }

    writeTo(&file, &w);
    writeTo(&file, &u);
    file.close();
}


int main()
{
    const int N = 10;
    const int T = 1000;
    const float delta_x = 1./(N-1);
    const float delta_t = 0.001;
    const int G = 1;

    vector<float> w0(N);
    vector<float> u0(N);
    vector<float> b(N);
    vector<float> h0(N);

    for (int j = 0; j < N; j++){
        float x = j*delta_x;
        // w0[j] = 4.0*x*(1-x) + 1;            // symmetric IC
        w0[j] = 6.75*x*x*(1-x)/2 + G;       // asymmetric IC
        b[j] = 0.25*sin(pi*x);
        h0[j] = w0[j] - b[j];

        u0[j] = 0;
    }

    solve(b, h0, u0, delta_x, delta_t, N, T);

    return 0;
}
