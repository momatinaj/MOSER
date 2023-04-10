#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/Graph.h"
#include "Escape/GetAllCounts.h"
#include <iostream>
#include <set>
#include <list>
#include <map>
#include <vector>
#include <random>
#include <chrono>

using namespace Escape;

struct Edge
{
    int src;
    int dst;
    int idx;
};
void update_counts_using_delta(double *nonInd, double *delta, int size)
{
    for (int i = 0; i < size; i++)
        nonInd[i] += delta[i];
}
int get_node_degree(CGraph &cg, VertexIdx node)
{
    int degree = cg.offsets[node + 1] - cg.offsets[node];
    if (cg.nbors[cg.offsets[node + 1] - 1] < 0) // check for invalid edge!
        degree--;
    return degree;
}

bool check_switch(CGraph &g, Edge ab, Edge cd)
{
    int a = ab.src, b = ab.dst, c = cd.src, d = cd.dst;
    if (a == c || a == d || b == c || b == d)
        return false;

    // check ad or cb exists in graph or not
    for (int j = g.offsets[a]; j < g.offsets[a + 1]; ++j)
        if (g.nbors[j] == d)
            return false;
    for (int j = g.offsets[c]; j < g.offsets[c + 1]; ++j)
        if (g.nbors[j] == b)
            return false;
    return true;
}

Edge random_edge_picker(CGraph &cg)
{
    // picks a random edge from the graph (only gets the edge on the even indexes)
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<int> distr(0, (cg.nVertices - 1));

    // get the first edge
    int src = distr(generator);
    int src_deg = get_node_degree(cg, src);
    while (src_deg <= 0)
    { // make sure the node has at least one valid neighbour!
        src = distr(generator);
        src_deg = get_node_degree(cg, src);
    }

    std::uniform_int_distribution<int> distr2(cg.offsets[src], cg.offsets[src + 1] - 1);
    int dst = cg.nbors[distr2(generator)];
    while (dst < 0) // make sure we don't choose an invalid neighbour
        dst = cg.nbors[distr2(generator)];

    Edge edge;
    edge.src = src;
    edge.dst = dst;

    return edge;
}

Edge find_switch_candidate(CGraph &cg, Edge ab)
{
    Edge cd;
    int RETRY_COUNT = 10;
    while (RETRY_COUNT > 0)
    {
        cd = random_edge_picker(cg);
        if (check_switch(cg, ab, cd))
            break;
        else
            RETRY_COUNT--;
    }
    if (RETRY_COUNT)
        return cd;
    else
    {
        std::cout << "Could not find a switch candidate!\n";
        cd.src = -1;
        return cd;
    }
}

int count_triangles_around_the_edge(CGraph &cg, Edge edge)
{
    int currVertex = edge.src, triangles = 0;
    for (int i = cg.offsets[currVertex]; i < cg.offsets[currVertex + 1]; ++i)
    {
        int neighbour = cg.nbors[i];
        if (neighbour == edge.dst || neighbour < 0)
            continue;
        for (int j = cg.offsets[edge.dst]; j < cg.offsets[edge.dst + 1]; j++)
        {
            if (neighbour == cg.nbors[j] && cg.nbors[j] != edge.src && cg.nbors[j] >= 0)
            {
                triangles++;
                break;
            }
        }
    }
    return triangles;
}

void update_3node_deletion(CGraph &cg, Edge edge, double (&delta)[4])
{
    int n = cg.nVertices;
    int deg_src = get_node_degree(cg, edge.src);
    int deg_dst = get_node_degree(cg, edge.dst);

    // gets the non-induced counts and update them
    // 0. (n*(n-1)*(n-2))/6;  // number of independent sets (doesn't change)
    delta[0] = 0;
    // 1. m*(n-2);    // number of plain edges
    delta[1] = -(n - 2);
    // 2. w = w+(deg*(deg-1))/2; // updating total wedge count
    // for the two ends of the edge degree reduced by one
    delta[2] = -(deg_src + deg_dst - 2);
    // 3. trianlge counts (should be zero)
    delta[3] = -count_triangles_around_the_edge(cg, edge);
}
void update_3node_addition(CGraph &cg, Edge edge, double (&delta)[4])
{
    int n = cg.nVertices;
    int deg_src = get_node_degree(cg, edge.src);
    int deg_dst = get_node_degree(cg, edge.dst);
    // gets the non-induced counts and update them
    // 0. (n*(n-1)*(n-2))/6;  // number of independent sets (doesn't change)
    delta[0] = 0;
    // 1. m*(n-2);    // number of plain edges
    delta[1] = (n - 2);
    // 2. w = w+(deg*(deg-1))/2; // updating total wedge count
    // for the two ends of the edge degree reduced by one
    delta[2] = deg_src + deg_dst;
    // 3. trianlge counts (should be zero)
    delta[3] = count_triangles_around_the_edge(cg, edge);
}
void simulate_deletion_on_CG(CGraph &cg, Edge edge)
{
    // 1. find the position of the src in dst nbors
    // 2. swap the src with the last nbor of dst
    // 3. make the last neighbour an invalid neighbour
    // 4. do this on the other side (dst, src)
    for (int i = cg.offsets[edge.src]; i < cg.offsets[edge.src + 1]; i++)
    {
        // skip unwanted neighbors
        if (cg.nbors[i] != edge.dst)
            continue;

        // swap the deleted nbor with the last element
        int last_nbor_idx = cg.offsets[edge.src + 1] - 1;
        cg.nbors[i] = cg.nbors[last_nbor_idx];
        // make the last idx an invalid node
        cg.nbors[last_nbor_idx] = -23;
        break;
    }
    for (int i = cg.offsets[edge.dst]; i < cg.offsets[edge.dst + 1]; i++)
    {
        // skip unwanted neighbors
        if (cg.nbors[i] != edge.src)
            continue;

        // swap the deleted nbor with the last element
        int last_nbor_idx = cg.offsets[edge.dst + 1] - 1;
        cg.nbors[i] = cg.nbors[last_nbor_idx];
        // make the last idx an invalid node
        cg.nbors[last_nbor_idx] = -23;
        break;
    }
    cg.nEdges -= 2;
}

void simulate_addition_on_CG(CGraph &cg, Edge edge)
{
    int last_idx = cg.offsets[edge.src + 1] - 1;
    if (cg.nbors[last_idx] != -23)
        std::cout << "SOMETHING WENT REALLY WRONG!\n";
    else
        cg.nbors[last_idx] = edge.dst;

    last_idx = cg.offsets[edge.dst + 1] - 1;
    if (cg.nbors[last_idx] != -23)
        std::cout << "SOMETHING WENT REALLY WRONG!\n";
    else
        cg.nbors[last_idx] = edge.src;

    cg.nEdges += 2;
}

bool one_full_switch_tracking(CGraph &cg, double (&nonInd)[4])
{
    auto t_switch_begin = std::chrono::high_resolution_clock::now();
    Edge ab = random_edge_picker(cg);
    Edge cd = find_switch_candidate(cg, ab);
    auto t_switch_end = std::chrono::high_resolution_clock::now();
    auto t_switch = std::chrono::duration_cast<std::chrono::nanoseconds>(t_switch_end - t_switch_begin);

    if (cd.src == -1)
        return false;

    double delta3[4];

    update_3node_deletion(cg, ab, delta3);
    update_counts_using_delta(nonInd, delta3, 4);

    simulate_deletion_on_CG(cg, ab);

    update_3node_deletion(cg, cd, delta3);
    update_counts_using_delta(nonInd, delta3, 4);

    simulate_deletion_on_CG(cg, cd);

    // SWITCH EDGES//
    Edge ad = {ab.src, cd.dst};
    Edge cb = {cd.src, ab.dst};

    // ADDITION//
    update_3node_addition(cg, ad, delta3);
    update_counts_using_delta(nonInd, delta3, 4);

    simulate_addition_on_CG(cg, ad);

    update_3node_addition(cg, cb, delta3);
    update_counts_using_delta(nonInd, delta3, 4);

    simulate_addition_on_CG(cg, cb);

    return true;
}

int main(int argc, char *argv[])
{
    auto t_profile_begin = std::chrono::high_resolution_clock::now();
    TriangleInfo trinfo;
    Graph g;
    if (loadGraph(argv[1], g, 1, IOFormat::escape))
        exit(1);

    int NUMBER_OF_STEPS;

    if (argc < 3)
    {
        std::cout << "You can specify the number of steps!\n";
        std::cout << "The default value is 10K. \n";
        NUMBER_OF_STEPS = 10000;
    }
    else
    {
        NUMBER_OF_STEPS = atoi(argv[2]);
    }

    if (loadGraph(argv[1], g, 1, IOFormat::escape))
        exit(1);

    printf("Loaded graph\n");
    CGraph cg = makeCSR(g);
    cg.sortById();
    printf("Converted to CSR\n");

    printf("Relabeling graph\n");
    CGraph cg_relabel = cg.renameByDegreeOrder();
    cg_relabel.sortById();
    printf("Creating DAG\n");
    CDAG dag = degreeOrdered(&cg_relabel);

    (dag.outlist).sortById();
    (dag.inlist).sortById();

    auto t_graph_load_end = std::chrono::high_resolution_clock::now();
    auto t_graph_load = std::chrono::duration_cast<std::chrono::nanoseconds>(t_graph_load_end - t_profile_begin);
    printf("Graph loaded in: %.3f seconds.\n", t_graph_load.count() * 1e-9);

    double nonInd[4];

    auto t_3count_begin = std::chrono::high_resolution_clock::now();
    // printf("Counting 3-vertex\n");
    trinfo = getAllThree(&cg_relabel, &dag, nonInd, true);
    auto t_3count_end = std::chrono::high_resolution_clock::now();
    auto t_3count = std::chrono::duration_cast<std::chrono::nanoseconds>(t_3count_end - t_3count_begin);
    printf("3 Nodes Counted in: %.3f seconds.\n", t_3count.count() * 1e-9);

    double full_res[NUMBER_OF_STEPS + 1][4];
    for (int j = 0; j < 4; j++)
        full_res[0][j] = nonInd[j];

    auto t_dynamic_begin = std::chrono::high_resolution_clock::now();
    for (int i = 0; i <= NUMBER_OF_STEPS; i++)
    {
        one_full_switch_tracking(cg, nonInd);
        for (int j = 0; j < 4; j++)
            full_res[i][j] = nonInd[j];
    }

    auto t_dynamic_end = std::chrono::high_resolution_clock::now();
    auto t_dynamic = std::chrono::duration_cast<std::chrono::nanoseconds>(t_dynamic_end - t_dynamic_begin);
    auto t_full = std::chrono::duration_cast<std::chrono::nanoseconds>(t_dynamic_end - t_profile_begin);
    printf("Track and Count took: %.3f seconds.\n", t_dynamic.count() * 1e-9);
    printf("Full Algorithm took: %.3f seconds.\n", t_full.count() * 1e-9);

    FILE *f = fopen("out.txt", "w");
    if (!f)
    {
        printf("could not write to output to out.txt\n");
        return 0;
    }
    fprintf(f, "%lld\n", cg_relabel.nVertices);
    fprintf(f, "%lld\n", cg_relabel.nEdges);
    fprintf(f, "%lld\n", NUMBER_OF_STEPS);

    for (int i = 0; i < NUMBER_OF_STEPS; i++)
    {
        for (int j = 0; j < 4; j++)
            fprintf(f, "%f ", full_res[i][j]);
        fprintf(f, "\n");
    }
    fclose(f);
}
