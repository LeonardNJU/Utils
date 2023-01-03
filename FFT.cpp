#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

void real2complex(vector<double> &real, vector<complex<double>> &complex) {
    complex.resize(real.size());
    for (int i = 0; i < real.size(); i++) {
        complex[i] = ::complex < double > (real[i], 0);
    }
}

void complex2real(vector<complex<double>> &complex, vector<double> &real) {
    real.resize(complex.size());
    for (int i = 0; i < complex.size(); i++) {
        real[i] = complex[i].real();
    }
}

void bit_reverse_copy(vector<complex<double>> &a, vector<complex<double>> &b) {
    int n = a.size();
    for (int i = 0; i < n; i++) {
        int j = 0;
        for (int k = 0; k < log2(n); k++) {
            j = (j << 1) | ((i >> k) & 1);
        }
        b[j] = a[i];
    }
}

//a is input
//fft stands for fft/ifft, true stands for fft, false stands for ifft
// A is an answer
//length should be handled by user
void FFT(vector<complex<double>> &a, vector<complex<double>> &A, bool fft = true) {
    int n = a.size();

    bit_reverse_copy(a, A);

    int m = 1;
    for (int s = 1; s <= log2(n); ++s) {
        m <<= 1;
        complex<double> w_m(cos(2 * M_PI / m), sin(2 * M_PI / m));
        if (!fft)w_m = complex<double>(cos(2 * M_PI / m), -sin(2 * M_PI / m));
        for (int k = 0; k < n; k += m) {
            complex<double> w(1);
            for (int j = 0; j < m / 2; ++j) {
                complex<double> t = w * A[k + j + m / 2];
                complex<double> u = A[k + j];
                A[k + j] = u + t;
                A[k + j + m / 2] = u - t;
                w *= w_m;
            }
        }
    }
    if (!fft) {
        for (int i = 0; i < n; ++i) {
            A[i] /= n;
        }
    }
}

vector<double> convolution(vector<double> &da, vector<double> &db) {
    vector<complex<double>> a, b;
    real2complex(da, a);
    real2complex(db, b);
    vector<complex<double>> A, B;
    vector<complex<double>> temp;

    int n = a.size() + b.size();
    int N;
    for (N = 1; N < n; N <<= 1);
    for (int i = a.size(); i < N; i++) {
        a.push_back(0);
    }
    for (int i = b.size(); i < N; i++) {
        b.push_back(0);
    }
    A.resize(N);
    B.resize(N);
    temp.resize(N);

    FFT(a, A);
    FFT(b, B);
    for (int i = 0; i < A.size(); ++i) {
        A[i] *= B[i];
    }
    FFT(A, temp, false);
    vector<double> ans;
    complex2real(temp, ans);
    return ans;
}

int test() {
    vector<double> a;
    vector<double> b;
    cout<<"input two polynomials(from const to high,2x+1 you type 1 2):"<<endl;
    cout<<"a:";
    double temp;
    while (cin >> temp) {
        if (!temp)break;
        a.push_back(temp);
    }
    cout<<"b:";
    while (cin>>temp){
        if (!temp)break;
        b.push_back(temp);
    }
    vector<double> ans = convolution(a, b);
    for (int i = ans.size() - 1; i >= 0; --i) {
        if (ans[i] >1e-2 && i != 0) {
            cout << ans[i] << "x^" << i << " +";
        } else if (ans[i] >1e-2 && i == 0) {
            cout << ans[i] << endl;
        }
    }
    return 0;
}