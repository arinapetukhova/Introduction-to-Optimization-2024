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
    while (rowsMarked.size() != n && colsMarked.size() != m) {
        for (int i = 0; i < n; ++i) {
            if (rowsMarked.find(i) != rowsMarked.end()) continue;
            int max = 0;
            for (int j = 0; j < m; ++j) {
                if (colsMarked.find(j) != colsMarked.end()) continue;
                if (grid[i][j] > max) max = grid[i][j];
            }
            rowMax[i] = max;
        }        
        for (int j = 0; j < m; ++j) {
            if (colsMarked.find(j) != colsMarked.end()) continue;
            int max = 0;
            for (int i = 0; i < n; ++i) {
                if (rowsMarked.find(i) != rowsMarked.end()) continue;
                if (grid[i][j] > max) max = grid[i][j];
            }
            colMax[j] = max;
        }      
        int mostNegative = 0;
        pair<int, int> coordinate;
        for (int i = 0; i < n; ++i) {
            if (rowsMarked.find(i) != rowsMarked.end()) continue;
            for (int j = 0; j < m; ++j) {
                if (colsMarked.find(j) != colsMarked.end()) continue;
                delta[i][j] = grid[i][j] - (rowMax[i] + colMax[j]);
                if (delta[i][j] <= mostNegative) {
                    mostNegative = delta[i][j];
                    coordinate.first = i;
                    coordinate.second = j;
                }
            }
        }
        allocation[coordinate.first][coordinate.second] = min(supply[coordinate.first],demand[coordinate.second]);
        if (supply[coordinate.first] < demand[coordinate.second]) {
            rowsMarked.insert(coordinate.first);
            demand[coordinate.second] -= supply[coordinate.first];
            supply[coordinate.first] = 0;
        } else {
            colsMarked.insert(coordinate.second);
            supply[coordinate.first] -= demand[coordinate.second];
            demand[coordinate.second] = 0;
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
    return ans;
}

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
    return ans;
}

vector<vector<int> > northWestCorner(const vector<int>& supply, const vector<int>& demand) {
    vector<int> S = supply;
    vector<int> D = demand;
    vector<vector<int> > x0(S.size(), vector<int>(D.size(), 0));
    
    int i = 0, j = 0;
    while (i < S.size() && j < D.size()) {
        int allocation = min(S[i], D[j]);
        x0[i][j] = allocation;
        S[i] -= allocation;
        D[j] -= allocation;

        if (S[i] == 0) i++;
        if (D[j] == 0) j++;
    }

    return x0;
}

void northWestAlgorithm(const vector<int>& supply, const vector<vector<int> >& cost, const vector<int>& demand) {
    vector<vector<int> > nwSolution = northWestCorner(supply, demand);
    cout << "\nVectors of initial basic feasible solution (North-West Corner method):\n";
    int n = nwSolution.size();
    int m = nwSolution[0].size();
    for (int i = 0; i < n; i++) {
        cout << "[";
        for (int j = 0; j < m; j++) {
            if (j == m - 1) {
                cout << nwSolution[i][j];
            } else {
                cout << nwSolution[i][j] << " ";
            }
        }
        if (i == n - 1) {
            cout << "]" << endl;
        } else {
            cout << "], ";
        }
    }
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
    vector<int> supply2(3);
    vector<int> demand2(4);
    vector<vector<int> > grid2(3, vector<int>(4));;
    int n = 3;
    int m = 4;
    int ans = 0;
    int ans2 = 0;
    int supSum = 0;
    int demSum = 0;
    bool checkNegative = false;
    cout << "Enter the vector of coefficients of supply:" << endl;
    for (int i = 0; i < n; i++) {
        cin >> supply[i];
        supply2[i] = supply[i];
        supSum += supply[i];
    }
    cout << "Enter the matrix of coefficients of costs (" << n << " rows and " << m << " columns):" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> grid[i][j];
            grid2[i][j] = grid[i][j];
            if (grid[i][j] < 0) {
                checkNegative = true;
            }
        }
    }
    cout << "Enter the vector of coefficients of demand:" << endl;
    for (int i = 0; i < m; i++) {
        cin >> demand[i];
        demand2[i] = demand[i];
        demSum += demand[i];
    }

    if (checkNegative) {
        cout << "The method is not applicable!";
        return 0;
    }

    if (supSum != demSum) {
        cout << "The problem is not balanced!";
        return 0;
    }

    printTable(grid, n, m, supply, demand);

    northWestAlgorithm(supply, grid, demand);
    Vogel(grid, supply, demand, ans, n, m);
    Russell(grid2, supply2, demand2, ans2, n, m);    
    
    return 0;
}