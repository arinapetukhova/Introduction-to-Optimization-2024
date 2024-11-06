//we consulted: https://www.geeksforgeeks.org/transportation-problem-set-4-vogels-approximation-method/

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

const int INF = 1000;

pair<vector<int>, vector<int> > findDiff(vector<vector<int> >& grid) {
    vector<int> rowDiff;
    vector<int> colDiff;

    for (int i = 0; i < grid.size(); i++) {
        vector<int> arr = grid[i];
        sort(arr.begin(), arr.end());
        rowDiff.push_back(arr[1] - arr[0]);
    }

    for (int col = 0; col < grid[0].size(); col++) {
        vector<int> arr;
        for (int i = 0; i < grid.size(); i++) {
            arr.push_back(grid[i][col]);
        }
        sort(arr.begin(), arr.end());
        colDiff.push_back(arr[1] - arr[0]);
    }

    return make_pair(rowDiff, colDiff);
}

int Vogel(vector<vector<int> >& grid, vector<int>& supply, vector<int>& demand, int& ans, int n, int m) {
    vector<vector<int> > x0(n, vector<int>(m));
    while (*max_element(supply.begin(), supply.end()) != 0 || *max_element(demand.begin(), demand.end()) != 0) {
        pair<vector<int>, vector<int> > diffs = findDiff(grid);
        vector<int> row = diffs.first;
        vector<int> col = diffs.second;

        int maxi1 = *max_element(row.begin(), row.end());
        int maxi2 = *max_element(col.begin(), col.end());

        if (maxi1 >= maxi2) {
            for (int ind = 0; ind < row.size(); ind++) {
                if (row[ind] == maxi1) {
                    int mini1 = *min_element(grid[ind].begin(), grid[ind].end());
                    for (int ind2 = 0; ind2 < grid[ind].size(); ind2++) {
                        if (grid[ind][ind2] == mini1) {
                            int mini2 = min(supply[ind], demand[ind2]);
                            ans += mini2 * mini1;
                            x0[ind][ind2] = mini2;
                            supply[ind] -= mini2;
                            demand[ind2] -= mini2;

                            if (demand[ind2] == 0) {
                                for (int r = 0; r < n; r++) {
                                    grid[r][ind2] = INF;
                                }
                            } else {
                                fill(grid[ind].begin(), grid[ind].end(), INF);
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        } else {
            for (int ind = 0; ind < col.size(); ind++) {
                if (col[ind] == maxi2) {
                    int mini1 = INF;
                    for (int j = 0; j < n; j++) {
                        mini1 = min(mini1, grid[j][ind]);
                    }

                    for (int ind2 = 0; ind2 < n; ind2++) {
                        if (grid[ind2][ind] == mini1) {
                            int mini2 = min(supply[ind2], demand[ind]);
                            ans += mini2 * mini1;
                            x0[ind2][ind] = mini2;
                            supply[ind2] -= mini2;
                            demand[ind] -= mini2;

                            if (demand[ind] == 0) {
                                for (int r = 0; r < n; r++) {
                                    grid[r][ind] = INF;
                                }
                            } else {
                                fill(grid[ind2].begin(), grid[ind2].end(), INF);
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    cout << "Vectors of initial basic feasible solution (Vogel's approximation method):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "[";
        for (int j = 0; j < m; j++) {
            if (j == m - 1) {
                cout << x0[i][j];
            } else {
                cout << x0[i][j] << " ";
            }
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

    Vogel(grid, supply, demand, ans, n, m);
    
    
    return 0;
}