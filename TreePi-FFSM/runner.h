// Class definition of runner class
#ifndef RUNNER_H
#define RUNNER_H

#include <vector>
#include "graphFFSM.h"
#include "adj_matrix.h"
#include "Graph.h"

class Runner {
  private:
    vector<GraphFFSM*> graph_db;
    int minsup;
    int graph_count;
    vector<AdjMatrix> freq_nodes;
    vector<AdjMatrix> freq_edges;
    vector<AdjMatrix> output;
  public:
    // Constructor
    Runner();
    Runner(vector<GraphFFSM*>& db, int ms);

    void displayOutput(unordered_map<Graph, int, GraphHasher> &freqTrees) const;
    // Setter
    void setGraphDB(vector<GraphFFSM*>& db) {graph_db = db; graph_count = db.size();}
    void setMinsup(int ms) {minsup = ms;}

    void run();

    vector<AdjMatrix> generateFrequentNodes();
    vector<AdjMatrix> generateFrequentEdges();
    void ffsmExplore(vector<AdjMatrix>& sub_cam, vector<AdjMatrix>& freq_cam);
}; 


#endif
