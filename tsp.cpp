#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <random>
#include <stack>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace std::chrono;

const int INF = 1000000000;
struct Edge {
    int w = INF, to = -1;
};

int dist(vector<double> &p1, vector<double> &p2) {
    return round(sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2)));
}
vector<vector<int>> mst(int n, vector<vector<double>> &points, vector<vector<int>> adjW) {
    int total_weight = 0;
    vector<bool> selected(n, false);
    vector<Edge> min_e(n);
    vector<vector<int>> adj(n);
    min_e[0].w = 0;

    for (int i=0; i<n; ++i) {
        int v = -1;
        for (int j = 0; j < n; ++j) {
            if (!selected[j] && (v == -1 || min_e[j].w < min_e[v].w))
                v = j;
        }

        if (min_e[v].w == INF) {
            cout << "No MST!" << endl;
            exit(0);
        }

        selected[v] = true;
        total_weight += min_e[v].w;
        if (min_e[v].to != -1) {
            adj[v].push_back(min_e[v].to);
            adj[min_e[v].to].push_back(v);
        }

        for (int to = 0; to < n; ++to) {
            if (adjW[v][to] < min_e[to].w)
                min_e[to] = {adjW[v][to], v};
        }
    }
    
    return adj;
}
vector<int> cutShort(int n, vector<int> tour) {
    vector<int> shortcut;
    bool visited[n] = {false};
    for(int i = 0; i < tour.size(); i++) {
        if(!visited[tour[i]]) {
            shortcut.push_back(tour[i]);
            visited[tour[i]] = true;
        }
    }
    return shortcut;
}

vector<int> greedyTour(int n , vector<vector<double>> &points) {
    vector<int> tour(n);
    bool used[n] = {false};
    int ran = rand() % n;
    used[ran] = true;
    int best;
    for(int i = 1; i < n; i++) {
        best = -1;
        for(int j = 0; j < n; j++) {
            if(!used[j] && (best == -1 || (dist(points[tour[i-1]], points[j]) < dist(points[tour[i-1]], points[best])))) {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
    return tour;
}

vector<int> eulerTour(vector<vector<int>> &adj){
    
    vector<int> path;
    int index = 0;
    stack<int> st;
    st.push(index);
    while(!st.empty()){
        int current = st.top();
        if(adj[current].size() == 0) {
            path.push_back(current);
            st.pop();
        } else {
            st.push(adj[current][0]);
            int toRemove = adj[current][0];
            vector<int>::iterator it = find(adj[toRemove].begin(), adj[toRemove].end(), current);
	        int index = distance(adj[toRemove].begin(), it);
            adj[toRemove].erase(adj[toRemove].begin()+index);
            adj[current].erase(adj[current].begin());
        }
    }
    return path; 
}

void minMatching(vector<vector<int>> &adj, vector<vector<int>> &adjW) {

    vector<int> oddVert;
    vector<bool> oddBool;
    
    for (int i = 0; i < adj.size(); i++) {
        if (adj[i].size() % 2) {
            oddVert.push_back(i);
            oddBool.push_back(false);
        }
    }
   
   
    int min = INT32_MAX, tmp, minIdx;
    for(int i = 0; i < oddVert.size(); i++) {
        if (oddBool[i])
            continue;
        min = INT32_MAX;
        for(int j = 0; j < oddVert.size(); j++) { 
            tmp = adjW[oddVert[i]][oddVert[j]];
            if(tmp < min && !oddBool[i] && !oddBool[j]) {
                min = tmp;
                minIdx = j;
            }
        }
        adj[oddVert[i]].push_back(oddVert[minIdx]);
        adj[oddVert[minIdx]].push_back(oddVert[i]);
        oddBool[i] = true;
        oddBool[minIdx] = true;
    }
}

int totalDist(int n, vector<int> &tour, vector<vector<int>> &adjW) {
    int sum = 0;
    for(int i = 0; i < n - 1; i++) {
        sum += adjW[tour[i]][tour[i+1]];
    }
    return sum;
}
vector<int> swap(int s, int f, vector<int> &tour) {
    vector<int> res(tour.size());
    /*cout << "s: " << s << ", f: " << f << endl;
    for(int i = 0; i < tour.size(); i++) {
        cout << tour[i] << ' ';
    }
    cout << endl;*/
    res = tour;
    int k = f;
    for(s;s <= f; s++) {
        res[s] = tour[k];
        k--;
    }

    /*for(int i = 0; i < res.size(); i++) {
        cout << res[i] << ' ';
    }
    cout << endl;*/
    return res;
}
high_resolution_clock::time_point start;
bool deadline() {
    auto stop1 = high_resolution_clock::now(); 
    auto duration1 = duration_cast<microseconds>(stop1 - start); 
    if(duration1.count() > 1970000){
        return true;
    }
    return false;
}
vector<int> twoOpt(int n, vector<vector<int>> &adjW, vector<int> &tour) {
    int best = totalDist(n, tour, adjW);
    int tmp;
    bool loop = true;
    vector<int> currentTour;
    start = high_resolution_clock::now();
    while(loop) {
        if(deadline()){
            loop = false;
            break;
        }
        for(int i = 0; i < n-1 && loop; i++) {
            if(deadline()){
                loop = false;
                break;
            }
            for(int j = i+1; j < n && loop; j++) {
 
                if(deadline())  {
                    loop = false;
                    break;
                }
                currentTour = swap(i, j, tour);
                tmp = totalDist(n, currentTour, adjW);
                if(tmp < best){
                    //cout << "found better path" << endl;
                    best = tmp;
                    //cout << "new best " << tmp << endl;
                    tour = currentTour;
                }

            }
        }
    }
    return tour;
}


int main(int argc, char **argv) {
    int n;
    cin >> n;
    vector<vector<double>> points(n, vector<double>(2));
    for(int i = 0; i < n; i++) {
        cin >> points[i][0];
        cin >> points[i][1];
    }

    vector<vector<int>> adjW(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            if(i == j) {
                adjW[i][j] = INF;
            } else {
                adjW[i][j] = dist(points[i], points[j]);
                adjW[j][i] = adjW[i][j];
            }
        }
    }
    //cout << "MST:" << endl;
    //vector<vector<int>> adj = mst(n, points, adjW);
    /*for (int i = 0; i < adj.size(); i++) {
        for (int j = 0; j < adj[i].size(); j++) {
            cout << adj[i][j] << ' ';
        }
        cout << endl;
    }*/
    //cout << "MIN MATCHING:" << endl;
    //minMatching(adj, adjW);
    /*for (int i = 0; i < adj.size(); i++) {
        for (int j = 0; j < adj[i].size(); j++) {
            cout << adj[i][j] << ' ';
        }
        cout << endl;
    }*/
    //vector<int> tour = eulerTour(adj);
    /*for (int i = 0; i < tour.size(); i++) {
        cout << tour[i] << " ";
    }*/
    //vector<int> res = cutShort(n, tour);
    vector<int> res = greedyTour(n, points);
    
    /*for (int i = 0; i < res.size(); i++) {
        cout << res[i] << endl;
    }*/
    
    vector<int> tour = twoOpt(n, adjW, res);

    for (int i = 0; i < n; i++) {
        cout << tour[i] << endl;
    }
    //cout << endl;

    return 0;
}
