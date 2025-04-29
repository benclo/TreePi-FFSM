// graph class header file
#ifndef GRAPHFFSM_H
#define GRAPHFFSM_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <unordered_map>
// #include "p1_helper.h"

using namespace std;

struct GraphNode {
  string label;
  unordered_map<int, int> adj_list;   // First int for adjacent GraphNode index, second int for edge weight
  // GraphNode(string l) : label(l), adj_list() {} 
};

class GraphFFSM {
  private:
    string name;
    vector<GraphNode> nodes;
    int node_count;
    int edge_count;
    string canonical_label;
  public:
    // Constructor
    GraphFFSM();

    // Desctructor
    ~GraphFFSM();
    
    // Setter and Mutator
    void setName(string n);
    bool addNode(int idx, string l);
    bool addEdge(int n1, int n2, int edge_weight);
    void setMatrix(vector<vector<string>>& m);

    // Getter
    string getName();
    GraphNode getNode(int);
    string getNodeLabel(int node_idx) {return nodes[node_idx].label;}
    string getEdgeLabel(int row, int col);
    int getNodeCount() {return node_count;}
    int getEdgeCount() {return edge_count;}

    // Display functions (for Debugging purpose)
    void displayGraph();

    // Debug purpose
    string genGraphCode();

    // string canonicalLabel();
    // string backTrackCanonicalLabel(vector<vector<string>>& m, vector<vector<int>>& partition, vector<int>& mapping, vector<vector<string>>& m_mapped, int p_cur);
    
    friend bool isIsomorphism(GraphFFSM& g, GraphFFSM& s);
    friend bool isSameGraph(GraphFFSM& g, GraphFFSM& s);
    // Decide if current graph is subgraph of another graph
    friend bool isSubgraph(GraphFFSM& g, GraphFFSM& s);
    friend bool backTrackSearch(GraphFFSM& g, GraphFFSM& s, int idx, vector<int> mapped);

};

#endif