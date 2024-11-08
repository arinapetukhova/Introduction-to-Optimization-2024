#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void printTable(const vector<int>& supply, const vector<int>& demand, const vector<vector<int>>& cost) {
    cout << "\nSupply/Demand Table:\n";
    for (int d : demand) cout << d << " " << setw(4);
    cout << "Demand/Supply";
    cout << "\n";

    for (size_t i = 0; i < supply.size(); i++) {
        for (size_t j = 0; j < demand.size(); j++) {
            cout << cost[i][j] << " " << setw(4);
        }
        cout << setw(4) << " | "<< supply[i];
        cout << "\n";
    }
}

vector<vector<int>> northWestCorner(const vector<int>& supply, const vector<int>& demand) {
    vector<int> S = supply;
    vector<int> D = demand;
    vector<vector<int>> x0(S.size(), vector<int>(D.size(), 0));
    
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

void northWestAlgorithm(const vector<int>& supply, const vector<vector<int>>& cost, const vector<int>& demand) {
    printTable(supply, demand, cost);

    vector<vector<int>> nwSolution = northWestCorner(supply, demand);
    cout << "\nNorth-West Corner Solution:\n";
    for (const auto& row : nwSolution) {
        for (int val : row) cout << setw(4) << val << " ";
        cout << "\n";
    }
}


int main() {
    const int supplySize = 3;
    const int demandSize = 4;

    vector<int> supply(supplySize);
    cout << "Enter the supply values vector: ";
    for (int i = 0; i < supplySize; i++) {
        cin >> supply[i];
    }

    vector<int> demand(demandSize);
    cout << "Enter the demand values vector: ";
    for (int i = 0; i < demandSize; i++) {
        cin >> demand[i];
    }

    vector<vector<int>> cost(supplySize, vector<int>(demandSize));
    cout << "Enter the cost values matrix:\n";
    for (int i = 0; i < supplySize; i++) {
        for (int j = 0; j < demandSize; j++) {
            cin >> cost[i][j];
        }
    }

    northWestAlgorithm(supply, cost, demand);

    return 0;
}
