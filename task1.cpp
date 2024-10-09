#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

class Simplex {
public:
    Simplex(const vector<vector<double> >& A, const vector<double>& b, const vector<double>& C, bool maximize, double epsilon);
    void printProblem();
    bool solve();
    void printSolution();

private:
    vector<vector<double> > tableau;
    vector<int> basicVariables;
    int numRows, numCols;
    bool maximize;
    double epsilon;
    int precision;

    bool checkOptimality();
    void pivot(int pivotRow, int pivotCol);
    int selectPivotColumn();
    int selectPivotRow(int pivotCol);
    
    bool checkInfeasibility();
    bool checkUnboundedness();
    bool isDegenerate();
};

Simplex::Simplex(const vector<vector<double> >& A, const vector<double>& b, const vector<double>& C, bool maximize, double epsilon) {
    numRows = A.size();
    numCols = A[0].size();
    this->maximize = maximize;
    this->epsilon = epsilon;

    precision = -log10(epsilon);
    tableau = vector<vector<double> >(numRows + 1, vector<double>(numCols + numRows + 1));

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            tableau[i][j] = A[i][j];
        }
        tableau[i][numCols + i] = 1;
        tableau[i][numCols + numRows] = b[i];
        basicVariables.push_back(numCols + i);
    }

    for (int j = 0; j < numCols; ++j) {
        tableau[numRows][j] = maximize ? -C[j] : C[j];
    }
}

bool Simplex::checkInfeasibility() {
    for (int i = 0; i < numRows; ++i) {
        if (tableau[i][numCols + numRows] < -epsilon) {
            bool allNonPositive = true;
            for (int j = 0; j < numCols; ++j) {
                if (tableau[i][j] > epsilon) {
                    allNonPositive = false;
                    break;
                }
            }
            if (allNonPositive) {
                cout << "The problem is infeasible!" << endl;
                return true;
            }
        }
    }
    return false;
}

bool Simplex::checkOptimality() {
    for (int j = 0; j < numCols + numRows; ++j) {
        if (tableau[numRows][j] < -epsilon) {
            return false;
        }
    }
    return true;
}

bool Simplex::isDegenerate() {
    for (int i = 0; i < numRows; ++i) {
        if (fabs(tableau[i][numCols + numRows]) < epsilon) {
            cout << "Degeneracy detected!" << endl;
            return true;
        }
    }
    return false;
}

void Simplex::pivot(int pivotRow, int pivotCol) {
    double pivotElement = tableau[pivotRow][pivotCol];

    for (int j = 0; j <= numCols + numRows; ++j) {
        tableau[pivotRow][j] /= pivotElement;
    }

    for (int i = 0; i <= numRows; ++i) {
        if (i != pivotRow) {
            double ratio = tableau[i][pivotCol];
            for (int j = 0; j <= numCols + numRows; ++j) {
                tableau[i][j] -= ratio * tableau[pivotRow][j];
            }
        }
    }

    basicVariables[pivotRow] = pivotCol;
}

int Simplex::selectPivotColumn() {
    int pivotCol = 0;
    for (int j = 1; j < numCols + numRows; ++j) {
        if (tableau[numRows][j] < tableau[numRows][pivotCol]) {
            pivotCol = j;
        }
    }
    return pivotCol;
}

int Simplex::selectPivotRow(int pivotCol) {
    int pivotRow = -1;
    double minRatio = numeric_limits<double>::max();

    for (int i = 0; i < numRows; ++i) {
        if (tableau[i][pivotCol] > epsilon) {
            double ratio = tableau[i][numCols + numRows] / tableau[i][pivotCol];
            if (ratio < minRatio) {
                minRatio = ratio;
                pivotRow = i;
            }
        }
    }

    return pivotRow;
}

void Simplex::printProblem() {
    cout << "\nProblem:\n";
    cout << (maximize ? "max" : "min") << " z = ";
    for (int k = 0; k < numCols; k++) {
        cout << (maximize ? -1 * tableau.at(tableau.size() - 1).at(k) : tableau.at(tableau.size() - 1).at(k)) << "*x_" << k + 1;
        (k != numCols - 1 ? cout << " + " : cout << "");
    }
    cout << "\nsubject to the constraints:\n";
    for (int j = 0; j < numRows; j++) {
        for (int i = 0; i < numCols; i++) {
            cout << tableau.at(j).at(i) << "*x_" << i + 1;
            (i != numCols - 1 ? cout << " + " : cout << " <= " << tableau.at(j).at(tableau.at(j).size() - 1));
        }
        cout << "\n";
    }
    cout << "\n";
}

bool Simplex::solve() {
    if (checkInfeasibility()) {
        return false;
    }
    
    while (!checkOptimality()) {
        if (isDegenerate()) {
            cout << "Warning: Degenerate solution!" << endl;
        }
        
        int pivotCol = selectPivotColumn();
        int pivotRow = selectPivotRow(pivotCol);
        if (pivotRow == -1) {
            cout << "The method is not applicable (unbounded solution)!" << endl;
            return false;
        }

        pivot(pivotRow, pivotCol);
    }

    return true;
}

void Simplex::printSolution() {
    cout << fixed << setprecision(precision);
    cout << "Solution: " << endl;
    
    vector<double> solution(numCols, 0.0);
    
    for (int i = 0; i < numRows; ++i) {
        if (basicVariables[i] < numCols) {
            solution[basicVariables[i]] = tableau[i][numCols + numRows];
        }
    }

    for (int j = 0; j < numCols; ++j) {
        cout << "x_" << j + 1 << " = " << solution[j] << endl;
    }

    cout << (maximize ? "Maximum" : "Minimum") << " value: " << (maximize ? 1 : -1) * tableau[numRows][numCols + numRows] << endl;
}

int main() {
    int m, n;
    cout << "Enter the number of constraints (m): ";
    cin >> m;
    cout << "Enter the number of variables (n): ";
    cin >> n;

    vector<vector<double> > A(m, vector<double>(n));
    vector<double> b(m);
    vector<double> C(n);

    cout << "Enter the coefficients of the constraints (A matrix):\n";
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the right-hand side values (b vector):\n";
    for (int i = 0; i < m; ++i) {
        cin >> b[i];
    }

    cout << "Enter the coefficients of the objective function (C vector):\n";
    for (int j = 0; j < n; ++j) {
        cin >> C[j];
    }

    bool maximize;
    char optType;
    cout << "Do you want to maximize or minimize the objective function? (m for maximize, n for minimize): ";
    cin >> optType;
    maximize = (optType == 'm');

    double epsilon;
    cout << "Enter approximation tolerance (epsilon): ";
    cin >> epsilon;

    Simplex simplex(A, b, C, maximize, epsilon);
    simplex.printProblem();
    if (simplex.solve()) {
        simplex.printSolution();
    }

    return 0;
}