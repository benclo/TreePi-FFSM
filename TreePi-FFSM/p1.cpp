#include <iostream>
#include <cstdlib>
#include <fstream>
#include <unordered_map>
#include "graphFFSM.h"
#include "p1_helper.h"
#include "runner.h"
#include <random>
#include <tuple>
#include <ctime>
#include <set>
#include <algorithm>
#include "FFSM.h"

using namespace std;

/* processInput function reads minsup value and graphs data from graph.txt
 * Input: vector of graph pointers, where to store graphs from input
 * Output: minsup value
 */
bool processInput(vector<GraphFFSM*>& graph_list, int& ms) {
    ifstream inputFile("graph.txt");  // Open the file
    if (!inputFile) {
        cout << "Error: Unable to open graph.txt!" << endl;
        return false;
    }

    int minsup;
    string holder;
    vector<string> read_input;
    GraphFFSM* new_graph = nullptr;

    // Read minsup value
    do {
        if (getline(inputFile, holder)) {
            read_input = splitString(holder);
        }
        else {
            return false;
        }
    } while (read_input.empty());

    if (read_input.size() == 1 && isMinsup(read_input)) {
        minsup = stoi(read_input[0], NULL, 10);
    }
    else {
        cout << "Minsup input is incorrect!\n";
        return false;
    }
    ms = minsup;

    int data_mode = 0; // Track input type: 1 = graph name, 2 = node, 3 = edge

    while (getline(inputFile, holder)) {
        read_input = splitString(holder);
        if (read_input.empty()) continue;  // Skip empty lines

        if (read_input.size() == 1) {
            new_graph = new GraphFFSM;
            new_graph->setName(read_input[0]);
            graph_list.push_back(new_graph);
            data_mode = 1;
        }
        else if (read_input.size() == 2) {
            if (data_mode != 1 && data_mode != 2) {
                cout << "Data_mode is incorrect, cannot add new node!\n";
                return false;
            }
            if (!isNode(read_input)) {
                cout << "Node input is incorrect!\n";
                return false;
            }
            new_graph->addNode(stoi(read_input[0], NULL, 10), read_input[1]);
            data_mode = 2;
        }
        else if (read_input.size() == 3) {
            if (data_mode != 2 && data_mode != 3) {
                cout << "Data_mode is incorrect, cannot add new edge!\n";
                return false;
            }
            if (!isEdge(read_input)) {
                cout << "Edge input is incorrect!\n";
                return false;
            }
            new_graph->addEdge(stoi(read_input[0], NULL, 10), stoi(read_input[1], NULL, 10), stoi(read_input[2], NULL, 10));
            data_mode = 3;
        }
    }

    // Validate read graphs
    for (GraphFFSM* g : graph_list) {
        int node_num = g->getNodeCount();
        for (int j = 0; j < node_num; j++) {
            if (g->getNodeLabel(j) == "") {
                cout << "Some nodes are missing from input file. Cannot proceed!\n";
                return false;
            }
        }
    }

    inputFile.close(); // Close the file
    return true;
}
// Function to generate connected graphs with numeric labels and write to a file
void generateChemicalGraphsFFSM(int minSupport = 1, int numGraphs = 250, int numVertices = 13, int numEdges = 15) {
    mt19937 rng(static_cast<unsigned>(time(0)));  // Random number generator
    uniform_int_distribution<int> vertexDist(0, numVertices - 1);  // Vertex selection
    uniform_int_distribution<int> nodeLabelDist(0, 3);  // Numeric node labels {0, 1, 2, 3}
    uniform_int_distribution<int> edgeWeightDist(1, 3); // Edge weights {1, 2, 3}

    // Open the output file
    ofstream outFile("graph.txt");

    // Write the minimum support at the beginning
    outFile << minSupport << "\n";

    // Generate the graphs
    for (int g = 0; g < numGraphs; ++g) {
        // Write the graph identifier
        outFile << "G" << g + 1 << "\n";

        // Node labels: Randomly assign numeric labels {0, 1, 2, 3} to each node
        vector<int> nodeLabels(numVertices);
        for (int i = 0; i < numVertices; ++i) {
            nodeLabels[i] = nodeLabelDist(rng);  // Labels are randomly chosen from {0, 1, 2, 3}
        }

        // Set to store unique edges
        set<pair<int, int>> edgeSet;
        vector<tuple<int, int, int>> edges;

        // Ensure connectivity by creating a spanning tree first
        vector<int> perm(numVertices);
        for (int i = 0; i < numVertices; ++i) perm[i] = i;
        shuffle(perm.begin(), perm.end(), rng);

        // Add edges to form a spanning tree
        for (int i = 1; i < numVertices; ++i) {
            int v1 = perm[i - 1];
            int v2 = perm[i];
            int weight = edgeWeightDist(rng);  // Random edge weight
            edges.emplace_back(v1, v2, weight);
            edgeSet.insert({ min(v1, v2), max(v1, v2) });
        }

        // Add remaining random edges ensuring no duplicates and no self-loops
        while (edges.size() < numEdges) {
            int v1 = vertexDist(rng);
            int v2 = vertexDist(rng);
            while (v1 == v2 || edgeSet.count({ min(v1, v2), max(v1, v2) })) {
                v1 = vertexDist(rng);
                v2 = vertexDist(rng);
            }
            int weight = edgeWeightDist(rng);  // Random edge weight
            edges.emplace_back(v1, v2, weight);
            edgeSet.insert({ min(v1, v2), max(v1, v2) });
        }

        // Write nodes with their labels
        for (int i = 0; i < numVertices; ++i) {
            outFile << i << " " << nodeLabels[i] << "\n";  // Node ID and label
        }

        // Write edges with weights
        for (const auto& edge : edges) {
            outFile << get<0>(edge) << " " << get<1>(edge) << " " << get<2>(edge) << "\n";
        }

        if (g < numGraphs - 1) outFile << "\n";  // Ensure no trailing newline after the last graph
    }

    // Close the output file
    outFile.close();
    cout << "Generated " << numGraphs << " connected graphs in the required format.\n";
}

unordered_map<Graph, int, GraphHasher> FFSM() {

  unordered_map<Graph, int, GraphHasher> freqTrees;
  generateChemicalGraphsFFSM();
  vector<GraphFFSM*> graph_ptrs;
  int minsup;
  bool read_successful = processInput(graph_ptrs, minsup);
  int n = graph_ptrs.size();

  if (read_successful == false) {
    for (int i = 0; i < n; i++) {
      delete graph_ptrs[i];
    }
    return freqTrees;
  }
  /*
  for (int i = 0; i < n; i++) {
    graph_ptrs[i]->displayGraph();
  }
  */
  Runner test_runner(graph_ptrs, minsup);
  cout << "Here1";
  test_runner.run();
  cout << "Here2";
  test_runner.displayOutput(freqTrees);

  for (int i = 0; i < n; i++) {
    delete graph_ptrs[i];
  }

  return freqTrees;
}
