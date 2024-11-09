#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <set>
using namespace std;

const int INF = 1000;

int Russell(vector<vector<int> >& grid, vector<int>& supply, vector<int>& demand, int& ans, int n, int m) {
    vector<vector<int> > delta(n, vector<int>(m));
    vector<int> rowMax(n,0);
    vector<int> colMax(m,0);
    set<int> rowsMarked;
    set<int> colsMarked;
    vector<vector<int> > allocation(n, vector<int>(m));
    while (rowsMarked.size() != n && colsMarked.size() != m) { // 
        for (int i = 0; i < n; ++i) {
            if (rowsMarked.find(i) != rowsMarked.end()) continue;
            int max = 0;
            for (int j = 0; j < m; ++j) {
                if (colsMarked.find(j) != colsMarked.end()) continue;
                if (grid[i][j] > max) max = grid[i][j];
            }
            rowMax[i] = max;
            //cout << rowMax[i] << " "; // TO BE DELETED
        }
        //cout << endl; // TO BE DELETED
        
        for (int j = 0; j < m; ++j) {
            if (colsMarked.find(j) != colsMarked.end()) continue;
            int max = 0;
            for (int i = 0; i < n; ++i) {
                if (rowsMarked.find(i) != rowsMarked.end()) continue;
                if (grid[i][j] > max) max = grid[i][j];
            }
            colMax[j] = max;
            //cout << colMax[j] << " "; // TO BE DELETED
        }
        //cout << endl; // TO BE DELETED
        
        int mostNegative = 0;
        pair<int, int> coordinate;
        for (int i = 0; i < n; ++i) {
            if (rowsMarked.find(i) != rowsMarked.end()) continue;
            for (int j = 0; j < m; ++j) {
                if (colsMarked.find(j) != colsMarked.end()) continue;
                delta[i][j] = grid[i][j] - (rowMax[i] + colMax[j]);
                if (delta[i][j] <= mostNegative) {
                    mostNegative = delta[i][j];
                    coordinate = {i, j};
                }
                //cout << delta[i][j] << " "; // TO BE DELETED
            }
            //cout << endl; // TO BE DELETED
        }
        //cout << mostNegative << endl; // TO BE DELETED
        allocation[coordinate.first][coordinate.second] = min(supply[coordinate.first],demand[coordinate.second]);
        //cout << coordinate.first + 1 << " " << coordinate.second + 1 << " " << allocation[coordinate.first][coordinate.second] << endl; // TO BE DELETED
        if (supply[coordinate.first] < demand[coordinate.second]) {
            rowsMarked.insert(coordinate.first);
            demand[coordinate.second] -= supply[coordinate.first];
            supply[coordinate.first] = 0;
            //cout << "Row " << coordinate.first + 1 << " was deleted" << endl; // TO BE DELETED
        } else {
            colsMarked.insert(coordinate.second);
            supply[coordinate.first] -= demand[coordinate.second];
            demand[coordinate.second] = 0;
            //cout << "Column " << coordinate.second + 1 << " was deleted" << endl; // TO BE DELETED
        }
    }
    
    cout << "Vectors of initial basic feasible solution (Russell's approximation method):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "[";
        for (int j = 0; j < m; j++) {
            if (j == m - 1) {
                cout << allocation[i][j];
            } else {
                cout << allocation[i][j] << " ";
            }
            ans += allocation[i][j]*grid[i][j];
        }
        if (i == n - 1) {
            cout << "]" << endl;
        } else {
            cout << "], ";
        }
    }
    cout << "The basic feasible solution is " << ans << endl;
    return ans;
}

void printTable(vector<vector<int> > grid, int n, int m, vector<int> supply, vector<int> demand) {
    cout << endl;
    cout << "Input parameter table:" << endl;
    int width = 10;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << setw(width) << grid[i][j];
        }
        cout << setw(width) << supply[i] << endl; 
    }
    cout << string((m + 1) * width, '-') << endl;
    for (int i = 0; i < m; i++) {
        cout << setw(width) << demand[i];
    }
    cout << endl << endl;
}

int main() {
    vector<vector<int> > grid(3, vector<int>(4));;
    vector<int> supply(3);
    vector<int> demand(4);
    int n = 3;
    int m = 4;
    int supSum = 0;
    int demSum = 0;
    bool checkNegative = false;
    cout << "Enter the vector of coefficients of supply:" << endl;
    for (int i = 0; i < n; i++) {
        cin >> supply[i];
        supSum += supply[i];
    }
    cout << "Enter the matrix of coefficients of costs (" << n << " rows and " << m << " columns):" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> grid[i][j];
            if (grid[i][j] < 0) {
                checkNegative = true;
            }
        }
    }
    cout << "Enter the vector of coefficients of demand:" << endl;
    for (int i = 0; i < m; i++) {
        cin >> demand[i];
        demSum += demand[i];
    }
    int ans = 0;

    if (checkNegative) {
        cout << "The method is not applicable!";
        return 0;
    }

    if (supSum != demSum) {
        cout << "The problem is not balanced!";
        return 0;
    }

    printTable(grid, n, m, supply, demand);

    Russell(grid, supply, demand, ans, n, m);
    
    
    return 0;
}