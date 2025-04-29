#ifndef FFSM_RUNNER_H
#define FFSM_RUNNER_H

#include <vector>
#include <string>
#include "Graph.h"

using namespace std;

// Forward declaration
class GraphFFSM;

/**
 * Reads minimum support and graphs from graph.txt.
 * @param graph_list Vector where parsed Graph pointers will be stored.
 * @param ms Reference to store the minsup value.
 * @return True if parsing is successful, false otherwise.
 */
bool processInput(std::vector<GraphFFSM*>& graph_list, int& ms);

/**
 * Generates synthetic connected chemical-style graphs and writes to graph.txt.
 * @param minSupport Minimum support value to write as first line in file.
 * @param numGraphs Number of graphs to generate.
 * @param numVertices Number of vertices per graph.
 * @param numEdges Number of edges per graph.
 */
void generateChemicalGraphsFFSM(int minSupport, int numGraphs, int numVertices, int numEdges);

/**
 * Entry point for running the FFSM graph mining process.
 * Generates graphs, parses input, runs mining algorithm.
 * @return Exit status (0 for success).
 */
unordered_map<Graph, int, GraphHasher> FFSM();

#endif // FFSM_RUNNER_H
