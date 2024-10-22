#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

using namespace std;

const int MAX_ITER = 1000;

struct InteriorPointResult {
    vector<double> x;
    double optimalValue;
};

vector<double> matVecMult(const vector<vector<double>>& A, const vector<double>& x, int mmV = 1) {
    int n = A.size();
    int m = x.size();
    vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[i] += mmV * A[i][j] * x[j];
        }
    }

    return result;
}

vector<vector<double>> matMatMult(const vector<vector<double>>& A, const vector<vector<double>>& D) {
    int m = A.size();
    int n = A[0].size();
    int p = D[0].size();

    vector<vector<double>> result(m, vector<double>(p, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * D[k][j];
            }
        }
    }

    return result;
}

vector<vector<double>> transpose(const vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    vector<vector<double>> transposed(cols, vector<double>(rows, 0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            transposed[j][i] = matrix[i][j];
        }
    }
    return transposed;
}

void swap(double& a, double& b) {
    double temp = a;
    a = b;
    b = temp;
}

void addMultipleOfRow(vector<vector<double>>& matrix, int sourceRow, double scaleFactor, int targetRow) {
    for (size_t i = 0; i < matrix[0].size(); ++i) {
        matrix[targetRow][i] += scaleFactor * matrix[sourceRow][i];
    }
}

vector<vector<double>> inverseMatrix(const vector<vector<double>>& matrix) {
    int n = matrix.size();

    vector<vector<double>> augmentedMatrix(n, vector<double>(2 * n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = matrix[i][j];
            if (i == j) {
                augmentedMatrix[i][j + n] = 1.0;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        double scale = augmentedMatrix[i][i];
        for (int j = 0; j < 2 * n; j++) {
            augmentedMatrix[i][j] /= scale;
        }
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = -augmentedMatrix[k][i];
                addMultipleOfRow(augmentedMatrix, i, factor, k);
            }
        }
    }

    vector<vector<double>> inverse(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmentedMatrix[i][j + n];
        }
    }

    return inverse;
}

vector<vector<double>> matrixDifference(const vector<vector<double>>& matrix1, const vector<vector<double>>& matrix2) {
    int rows = matrix1.size();
    int cols = matrix1[0].size();

    vector<vector<double>> result(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }

    return result;
}

vector<vector<double>> projectionMatrix(const vector<vector<double>>& A) {
    vector<vector<double>> A_transpose = transpose(A);
    vector<vector<double>> A_times_Atranspose = matMatMult(A, A_transpose);
    int n = A.size();

    vector<vector<double>> A_times_Atranspose_inverse = inverseMatrix(A_times_Atranspose);
    vector<vector<double>> preproj = matMatMult(matMatMult(A_transpose, A_times_Atranspose_inverse), A);
    vector<vector<double>> I(preproj.size(), vector<double>(preproj.size(), 0.0));
    for (int i = 0; i < preproj.size(); ++i) {
        I[i][i] = 1;
    }
    vector<vector<double>> projection = matrixDifference(I, preproj);
    return projection;
}

InteriorPointResult interiorPointMethod(const vector<double>& C, const vector<vector<double>>& A, 
const vector<double>& b, const vector<double>& initPoint, double alpha, double EPSILON, char maxMin) {
    int n = C.size();
    vector<double> x = initPoint;
    int iteration = 0;

    while (iteration < MAX_ITER) {
        vector<vector<double>> D(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            D[i][i] = x[i];
        }
        vector<double> c_tilde = (maxMin == 'm' ? matVecMult(D, C) : matVecMult(D, C, -1));
        vector<vector<double>> a_tilde = matMatMult(A, D); 

        vector<vector<double>> P = projectionMatrix(a_tilde);
        vector<double> c_p = matVecMult(P, c_tilde);
        int idx = -1;
        double maxAbsVal = 0.0;
        for (int i = 0; i < n; i++) {
            if (c_p[i] < 0 && abs(c_p[i]) > maxAbsVal) {
                maxAbsVal = abs(c_p[i]);
                idx = i;
            }
        }
        if (idx == -1) {
            break;
        }
        for (int i = 0; i < n; i++) {
            c_p[i] = 1 + (alpha * c_p[i])/maxAbsVal;
        }
        vector<double> x_old = x;
        x = matVecMult(D, c_p);

        double maxDiffNew = 0.0;
        for (int i = 0; i < n; ++i) {
            maxDiffNew = max(maxDiffNew, abs(x[i] - x_old[i]));
        }
        if (maxDiffNew < EPSILON) {
            break;
        }
        iteration++;
    }

    double optimalValue = 0.0;
    for (int i = 0; i < n; ++i) {
        optimalValue += C[i] * x[i];
    }
    if (maxMin == 'n') {
        for (double xi : x) {
            xi *= -1;
        }
    }
    return {x, optimalValue};
}

bool isFeasible(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b, double EPSILON) {
    vector<double> Ax = matVecMult(A, x);
    for (int i = 0; i < Ax.size(); i++) {
        if (Ax[i] > b[i]) {
            return false;
        }
    }
    return true;
}

int main() {

    /*int xSize = 2;
    vector<double> C = {1, 1, 0, 0};
    vector<vector<double>> A = {{2, 4, 1, 0}, {1, 3, 0, -1}};
    vector<double> b = {16, 9};
    vector<double> x_initial = {0.5, 3.5, 1, 2};
    char maxMin = 'm';
    double eps = 0.000001;*/

    /*int xSize = 2;
    vector<double> C = {3, 2, 0, 0};
    vector<vector<double>> A = {{1, 1, 1, 0}, {2, 1, 0, 1}};
    vector<double> b = {4, 5};
    vector<double> x_initial = {1, 1, 2, 2};
    char maxMin = 'm';
    double eps = 0.000001;*/

    /*int xSize = 3;
    vector<double> C = {2, 5, 7};
    vector<vector<double>> A = {{1, 2, 3}};
    vector<double> b = {6};
    vector<double> x_initial = {1, 1, 1};
    char maxMin = 'm';
    double eps = 0.000001;*/

    /*int xSize = 3;
    vector<double> C = {-2, 2, -6, 0, 0, 0};
    vector<vector<double>> A = {{2, 1, -2, 1, 0, 0}, {1, 2, 4, 0, 1, 0}, {1, -1, 2, 0, 0, 1}};
    vector<double> b = {24, 23, 10};
    vector<double> x_initial = {1, 1, 1, 1, 16, 8};
    char maxMin = 'n';
    double eps = 0.00001;*/

    /*int xSize = 3;
    vector<double> C = {9, 10, 16, 0, 0, 0};
    vector<vector<double>> A = {{18, 15, 12, 1, 0, 0}, {6, 4, 8, 0, 1, 0}, {5, 3, 3, 0, 0, 1}};
    vector<double> b = {360, 192, 180};
    vector<double> x_initial = {1, 1, 1, 315, 174, 169};
    char maxMin = 'm'; //m if maximaize, n if minimize
    double eps = 0.0001;*/
    double alpha1 = 0.5;
    double alpha2 = 0.9;

    if (!isFeasible(A, x_initial, b, eps)) {
        cout << "Initial point is not feasible!" << endl;
        return -1;
    }

    InteriorPointResult result1 = interiorPointMethod(C, A, b, x_initial, alpha1, eps, maxMin);
    if (result1.optimalValue == -numeric_limits<double>::infinity()) {
        cout << "The problem does not have a solution!" << endl;
    } else {
        cout << "Solution with alpha = 0.5:" << endl;
        for (int i = 0; i < xSize; i++) {
            cout << result1.x[i] << " ";
        }
        cout << endl << "Optimal value: " << result1.optimalValue << endl;
    }

    InteriorPointResult result2 = interiorPointMethod(C, A, b, x_initial, alpha2, eps, maxMin);
    if (result2.optimalValue == -numeric_limits<double>::infinity()) {
        cout << "The problem does not have a solution!" << endl;
    } else {
        cout << "Solution with alpha = 0.9:" << endl;
        for (int i = 0; i < xSize; i++) {
            cout << result2.x[i] << " ";
        }
        cout << endl << "Optimal value: " << result2.optimalValue << endl;
    }

    return 0;
}