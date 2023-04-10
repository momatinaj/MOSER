#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/FourVertex.h"
#include "Escape/Conversion.h"
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

struct FourInfo
{
    int fourCycleCount;
    int fourCliqueCount;
};

void print_edge(Edge edge)
{
    std::cout << edge.src << " " << edge.dst << "\n";
}

void printGraph(Graph &g)
{
    std::cout << "printing G(" << g.nVertices << ", " << g.nEdges << ")\n";
    for (int i = 0; i < g.nEdges; i++)
    {
        if (i % 2 == 0)
            std::cout << g.srcs[i] << "  " << g.dsts[i] << "\n";
    }
    std::cout << "--------------------\n";
}

void printCGraph(CGraph &g)
{
    std::cout << "printing CGraph(" << g.nVertices << ", " << g.nEdges << ")\n";
    for (int i = 0; i < g.nVertices; ++i)
    {
        std::cout << i << " :\n";
        for (int j = g.offsets[i]; j < g.offsets[i + 1]; ++j)
            std::cout << g.nbors[j] << " ";
        std::cout << "\n";
    }
    std::cout << "--------------------\n";
}

void print_count_results(double *result, int size)
{
    std::cout << "printing count results :\n";
    for (int i = 0; i < size; i++)
    {
        std::cout << result[i] << " ";
    }
    std::cout << "\n";
}

bool isEdge(CGraph &cg, VertexIdx src, VertexIdx dst)
{
    if (dst < 0 || src < 0)
        return false;

    for (int i = cg.offsets[src]; i < cg.offsets[src + 1]; i++)
        if (dst == cg.nbors[i])
            return true;

    return false;
}

int get_node_degree(CGraph &cg, VertexIdx node)
{
    int degree = cg.offsets[node + 1] - cg.offsets[node];
    if (cg.nbors[cg.offsets[node + 1] - 1] < 0) // check for invalid edge!
        degree--;
    return degree;
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

int count_3paths_around_the_edge(CGraph &cg, Edge edge)
{

    int currVertex = edge.src;
    int deg_dst = get_node_degree(cg, edge.dst);

    int d = 0;
    for (int i = cg.offsets[currVertex]; i < cg.offsets[currVertex + 1]; ++i)
    {
        int neighbour = cg.nbors[i];
        if (neighbour == edge.dst || neighbour < 0)
            continue;
        int deg_n = get_node_degree(cg, neighbour);

        d += -deg_n + 1;
    }

    currVertex = edge.dst;
    for (int i = cg.offsets[currVertex]; i < cg.offsets[currVertex + 1]; ++i)
    {
        int neighbour = cg.nbors[i];
        if (neighbour == edge.src || neighbour < 0)
            continue;

        int deg_n = get_node_degree(cg, neighbour);
        d += -deg_n + 1;
    }
    return d;
}

int degree_on_one_hop_neighbourhood(CGraph &cg, Edge edge)
{
    int currVertex = edge.src;
    int degrees = 0;
    for (int i = cg.offsets[currVertex]; i < cg.offsets[currVertex + 1]; ++i)
    {
        int neighbour = cg.nbors[i];
        if (neighbour == edge.dst || neighbour < 0)
            continue;
        for (int j = cg.offsets[edge.dst]; j < cg.offsets[edge.dst + 1]; j++)
        {
            if (neighbour == cg.nbors[j] && cg.nbors[j] != edge.src && cg.nbors[j] >= 0)
            {
                int deg_n = get_node_degree(cg, neighbour);
                degrees += deg_n - 2;
                break;
            }
        }
    }
    return degrees;
}

int edge_triangle_count_delta(CGraph &cg, Edge edge, bool add_mode = false)
{
    int currVertex = edge.src;
    int etd = 0; // edge_triangles_delta = (1 - triangles count around the edge)
    int triangles = 0;
    for (int i = cg.offsets[currVertex]; i < cg.offsets[currVertex + 1]; ++i)
    {
        int neighbour = cg.nbors[i];
        if (neighbour == edge.dst || neighbour < 0)
            continue;
        for (int j = cg.offsets[edge.dst]; j < cg.offsets[edge.dst + 1]; j++)
        {
            if (neighbour == cg.nbors[j] && cg.nbors[j] != edge.src && cg.nbors[j] >= 0)
            {
                Edge candidate_1 = {currVertex, neighbour};
                triangles = count_triangles_around_the_edge(cg, candidate_1);
                if (add_mode)
                    etd += triangles;
                else
                    etd += 1 - triangles;

                Edge candidate_2 = {edge.dst, neighbour};
                triangles = count_triangles_around_the_edge(cg, candidate_2);
                if (add_mode)
                    etd += triangles;
                else
                    etd += 1 - triangles;

                triangles = 0;
                break;
            }
        }
    }
    return etd;
}

int triangle_around_node(CGraph &cg, int node, int skip_node = -1)
{
    int triangles = 0;
    int skipped = 0;
    for (int i = cg.offsets[node]; i < cg.offsets[node + 1]; ++i)
    {

        int neighbour = cg.nbors[i];
        if (neighbour < 0)
            continue;
        for (int j = cg.offsets[neighbour]; j < cg.offsets[neighbour + 1]; ++j)
        {
            int endpoint = cg.nbors[j]; // now we have three nodes (node, neighbour, endpoint)
            if (endpoint == node || endpoint < 0)
                continue;
            if (isEdge(cg, node, endpoint))
                triangles++;
            else if (endpoint == skip_node)
                skipped++;
        }
    }
    return (triangles / 2) + skipped;
}

FourInfo fourInfo_delta_calc(CGraph &cg, Edge edge)
{
    int currVertex = edge.src;
    int cycles = 0;
    int cliques = 0;
    for (int i = cg.offsets[currVertex]; i < cg.offsets[currVertex + 1]; ++i)
    {
        int neighbour = cg.nbors[i];
        if (neighbour == edge.dst || neighbour < 0)
            continue;
        for (int j = cg.offsets[edge.dst]; j < cg.offsets[edge.dst + 1]; j++)
        {
            if (cg.nbors[j] >= 0 && neighbour != cg.nbors[j] && cg.nbors[j] != edge.src && isEdge(cg, neighbour, cg.nbors[j]))
            {
                cycles++;
                if (isEdge(cg, neighbour, edge.dst) && isEdge(cg, edge.src, cg.nbors[j]))
                {
                    cliques++;
                }
            }
        }
    }
    FourInfo result = {cycles, cliques / 2};
    return result;
}

void test_diamond_count_logic(CGraph &cg)
{
    std::cout << "---------- TESTING DIAMOND LOGIC ----------\n";
    int diamond_count = 0;
    for (int node = 0; node < cg.nVertices; node++)
    {
        for (int i = cg.offsets[node]; i < cg.offsets[node + 1]; ++i)
        {
            int neighbour = cg.nbors[i];
            Edge edge = {node, neighbour};
            int tri = count_triangles_around_the_edge(cg, edge);
            diamond_count += (tri * (tri - 1)) / 2;
            std::cout << ">>> on edge: <" << edge.src << ", " << edge.dst << ">\n";
            std::cout << ">> triangles: " << tri << "  DC: " << diamond_count << "\n";
        }
    }
    std::cout << "Diamond Count = " << diamond_count / 2 << "\n";
    std::cout << "---------- ---------- ---------- ----------\n";
}

int test_3path_logic(CGraph &cg, Edge edge)
{
    std::cout << "---------- TESTING 3PATH LOGIC ----------\n";
    int result = 0;
    for (int node = 0; node < cg.nVertices; node++)
    {
        int node_degree = cg.offsets[node + 1] - cg.offsets[node];
        if (node == edge.src || node == edge.dst)
        {
            node_degree--;
        }

        for (int j = cg.offsets[node]; j < cg.offsets[node + 1]; j++)
        {
            int neighbour = cg.nbors[j];
            int nei_degree = cg.offsets[neighbour + 1] - cg.offsets[neighbour];
            if (neighbour == edge.src || neighbour == edge.dst)
                nei_degree--;
            if ((neighbour == edge.src && node == edge.dst) || (neighbour == edge.dst && node == edge.src))
                continue;

            result += (node_degree - 1) * (nei_degree - 1);

            std::cout << "NODE FOUND (" << node << ", " << neighbour << ") => RESULT:" << result << "\n";
        }
    }
    return (result / 2);
    std::cout << "---------- ---------- ---------- ----------\n";
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

void update_4node_deletion(CGraph &cg, Edge edge, double (&nonInd3)[4], double (&delta)[11], double w1, double t1)
{
    // w1 is wedge count before deletion w2 is wedge count after deletion
    // t1 is triangle count before deletion t2 is triangle count after deletion
    int n = cg.nVertices;
    int m = cg.nEdges / 2;
    int deg_src = get_node_degree(cg, edge.src);
    int deg_dst = get_node_degree(cg, edge.dst);
    double w2 = nonInd3[2];
    double t2 = nonInd3[3];

    // nonInd[0] = (n*(n-1)*(n-2)*(n-3))/24;  // number of independent sets
    delta[0] = 0;
    // nonInd[1] = m*((n-2)*(n-3)/2);    // number of only edges
    delta[1] = -((n - 2) * (n - 3) / 2);
    // nonInd[2] = (m*(m-1)/2) - w; // number of matchings
    delta[2] = w1 + 1 - m - w2; // check this! answer is -14
    // nonInd[3] = w*(n-3);        // number of only wedges
    delta[3] = (w2 - w1) * (n - 3);
    // nonInd[4] = t*(n-3);        // number of only triangles
    delta[4] = (t2 - t1) * (n - 3);

    double deg_src_delta = 0.5 * (-(deg_src * deg_src) + 3 * deg_src - 2);
    double deg_dst_delta = 0.5 * (-(deg_dst * deg_dst) + 3 * deg_dst - 2);
    delta[5] = deg_dst_delta + deg_src_delta;

    delta[6] = count_3paths_around_the_edge(cg, edge) - ((deg_src - 1) * (deg_dst - 1)) - 3 * (t2 - t1);

    // ret.tailedtris += (degi - 2) * info.perVertex[i];
    int sigma_neighbourhood_degrees = -degree_on_one_hop_neighbourhood(cg, edge);
    // triangle count around the edge = t2-t1    // degree changes
    int src_change = (deg_src - 2) * (t2 - t1); // triangle_around_node_after_deletion(edge.src, edge, cg)
    int dst_change = (deg_dst - 2) * (t2 - t1);
    int tri_difference = triangle_around_node(cg, edge.src) + triangle_around_node(cg, edge.dst) - 2 * (t1 - t2); // triangle_around_node_after_deletion(edge.dst, edge, cg)
    delta[7] = sigma_neighbourhood_degrees + src_change + dst_change - tri_difference;

    // cycles
    FourInfo res4 = fourInfo_delta_calc(cg, edge);
    delta[8] = -res4.fourCycleCount;

    int t3 = t1 - t2;
    delta[9] = -((t3 * (t3 - 1)) / 2) + edge_triangle_count_delta(cg, edge);

    // cliques
    delta[10] = -res4.fourCliqueCount;
}

void update_4node_addition(CGraph &cg, Edge edge, double (&nonInd3)[4], double (&delta)[11], double w1, double t1)
{
    // w1 is wedge count before addition w2 is wedge count after addition
    // t1 is triangle count before addition t2 is triangle count after addition
    int n = cg.nVertices;
    int m = cg.nEdges / 2;
    int deg_src = get_node_degree(cg, edge.src);
    int deg_dst = get_node_degree(cg, edge.dst);
    double w2 = nonInd3[2];
    double t2 = nonInd3[3];

    // nonInd[0] = (n*(n-1)*(n-2)*(n-3))/24;  // number of independent sets
    delta[0] = 0;
    // nonInd[1] = m*((n-2)*(n-3)/2);    // number of only edges
    delta[1] = ((n - 2) * (n - 3) / 2);
    // nonInd[2] = (m*(m-1)/2) - w; // number of matchings
    delta[2] = m - (w2 - w1); // check this! answer is -14
    // nonInd[3] = w*(n-3);        // number of only wedges
    delta[3] = (w2 - w1) * (n - 3);
    // nonInd[4] = t*(n-3);        // number of only triangles
    delta[4] = (t2 - t1) * (n - 3);

    double deg_src_delta = 0.5 * ((deg_src * deg_src) - deg_src);
    double deg_dst_delta = 0.5 * ((deg_dst * deg_dst) - deg_dst);
    delta[5] = deg_dst_delta + deg_src_delta;

    delta[6] = -count_3paths_around_the_edge(cg, edge) + ((deg_src * deg_dst)) - 3 * (t2 - t1);

    // ret.tailedtris += (degi - 2) * info.perVertex[i];
    int sigma_neighbourhood_degrees = degree_on_one_hop_neighbourhood(cg, edge);
    // triangle count around the edge = t2-t1    // degree changes
    int src_change = (deg_src - 2) * (t2 - t1); // triangle_around_node_after_deletion(edge.src, edge, cg)
    int dst_change = (deg_dst - 2) * (t2 - t1);
    int tri_difference = triangle_around_node(cg, edge.src, edge.dst) + triangle_around_node(cg, edge.dst, edge.src);
    delta[7] = sigma_neighbourhood_degrees + src_change + dst_change + tri_difference;

    // cycles
    FourInfo res4 = fourInfo_delta_calc(cg, edge);
    delta[8] = res4.fourCycleCount;

    int t3 = t2 - t1;
    delta[9] = ((t3 * (t3 - 1)) / 2) + edge_triangle_count_delta(cg, edge, true);

    // cliques
    delta[10] = res4.fourCliqueCount;
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

void update_counts_using_delta(double *nonInd, double *delta, int size)
{
    for (int i = 0; i < size; i++)
        nonInd[i] += delta[i];
}

bool one_full_switch_tracking(CGraph &cg, double (&nonInd_three)[4], double (&nonInd_four)[11])
{
    auto t_switch_begin = std::chrono::high_resolution_clock::now();
    Edge ab = random_edge_picker(cg);
    Edge cd = find_switch_candidate(cg, ab);
    auto t_switch_end = std::chrono::high_resolution_clock::now();
    auto t_switch = std::chrono::duration_cast<std::chrono::nanoseconds>(t_switch_end - t_switch_begin);

    if (cd.src == -1)
        return false;

    double delta3[4], delta4[11];
    double w1 = nonInd_three[2], t1 = nonInd_three[3];

    update_3node_deletion(cg, ab, delta3);
    update_counts_using_delta(nonInd_three, delta3, 4);
    // print_count_results(nonInd_three, 4);

    update_4node_deletion(cg, ab, nonInd_three, delta4, w1, t1);
    update_counts_using_delta(nonInd_four, delta4, 11);
    // print_count_results(nonInd_four, 11);

    simulate_deletion_on_CG(cg, ab);
    // printCGraph(cg);

    w1 = nonInd_three[2];
    t1 = nonInd_three[3];
    update_3node_deletion(cg, cd, delta3);
    update_counts_using_delta(nonInd_three, delta3, 4);
    // print_count_results(nonInd_three, 4);

    update_4node_deletion(cg, cd, nonInd_three, delta4, w1, t1);
    update_counts_using_delta(nonInd_four, delta4, 11);
    // print_count_results(nonInd_four, 11);

    simulate_deletion_on_CG(cg, cd);
    // printCGraph(cg);

    // SWITCH EDGES//
    Edge ad = {ab.src, cd.dst};
    // print_edge(ad);
    Edge cb = {cd.src, ab.dst};
    // print_edge(cb);
    // ADDITION//

    w1 = nonInd_three[2];
    t1 = nonInd_three[3];
    update_3node_addition(cg, ad, delta3);
    update_counts_using_delta(nonInd_three, delta3, 4);
    // print_count_results(nonInd_three, 4);

    update_4node_addition(cg, ad, nonInd_three, delta4, w1, t1);
    update_counts_using_delta(nonInd_four, delta4, 11);
    // print_count_results(nonInd_four, 11);

    simulate_addition_on_CG(cg, ad);
    // printCGraph(cg);

    w1 = nonInd_three[2];
    t1 = nonInd_three[3];
    update_3node_addition(cg, cb, delta3);
    update_counts_using_delta(nonInd_three, delta3, 4);
    // print_count_results(nonInd_three, 4);

    update_4node_addition(cg, cb, nonInd_three, delta4, w1, t1);
    update_counts_using_delta(nonInd_four, delta4, 11);
    // print_count_results(nonInd_four, 11);

    simulate_addition_on_CG(cg, cb);
    // printCGraph(cg);

    return true;
}

void add_to_output_file(int remaining_steps, double full_res_three[][4], double full_res_four[][11])
{
    FILE *f = fopen("out.txt", "a");

    for (int i = 0; i < remaining_steps; i++)
    {
        for (int j = 0; j < 4; j++)
            fprintf(f, "%f ", full_res_three[i][j]);
        fprintf(f, "\n");
    }
    fprintf(f, "-------------------------- \n");
    for (int i = 0; i < remaining_steps; i++)
    {
        for (int j = 0; j < 11; j++)
            fprintf(f, "%f ", full_res_four[i][j]);
        fprintf(f, "\n");
    }
    fprintf(f, "-------------------------- \n");
    fclose(f);
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

    CGraph cg = makeCSR(g);
    cg.sortById();

    FILE *f = fopen("out.txt", "w");
    if (!f)
    {
        printf("could not write to output to out.txt\n");
        exit(1);
    }
    fprintf(f, "%s\n", argv[1]);
    fprintf(f, "%lld\n", NUMBER_OF_STEPS);
    fprintf(f, "%lld\n", cg.nVertices);
    fprintf(f, "%lld\n", cg.nEdges);
    fclose(f);

    CGraph cg_relabel = cg.renameByDegreeOrder();
    cg_relabel.sortById();

    CDAG dag = degreeOrdered(&cg_relabel);
    (dag.outlist).sortById();
    (dag.inlist).sortById();

    auto t_graph_load_end = std::chrono::high_resolution_clock::now();
    auto t_graph_load = std::chrono::duration_cast<std::chrono::nanoseconds>(t_graph_load_end - t_profile_begin);
    printf("Graph loaded in: %.3f seconds.\n", t_graph_load.count() * 1e-9);

    double nonInd_three[4], nonInd_four[11];

    // add code to bypass seg fault for too much mem alloc
    int memsize = std::min(NUMBER_OF_STEPS, 50000);
    int overflow = std::max(0, NUMBER_OF_STEPS - 50000);

    double full_res_four[memsize][11];
    double full_res_three[memsize][4];

    auto t_3count_begin = std::chrono::high_resolution_clock::now();
    trinfo = getAllThree(&cg_relabel, &dag, nonInd_three, true);
    auto t_3count_end = std::chrono::high_resolution_clock::now();
    auto t_3count = std::chrono::duration_cast<std::chrono::nanoseconds>(t_3count_end - t_3count_begin);
    printf("3 Nodes Counted in: %.3f seconds.\n", t_3count.count() * 1e-9);

    auto t_4count_begin = std::chrono::high_resolution_clock::now();
    getAllFour(&cg_relabel, &dag, nonInd_four, trinfo);
    auto t_4count_end = std::chrono::high_resolution_clock::now();
    auto t_4count = std::chrono::duration_cast<std::chrono::nanoseconds>(t_4count_end - t_4count_begin);
    printf("4 Nodes Counted in: %.3f seconds.\n", t_4count.count() * 1e-9);

    for (int j = 0; j < 4; j++)
        full_res_three[0][j] = nonInd_three[j];
    for (int j = 0; j < 11; j++)
        full_res_four[0][j] = nonInd_four[j];

    auto t_dynamic_begin = std::chrono::high_resolution_clock::now();
    int i = 1;
    while (true)
    {
        if (memsize == i)
        {
            add_to_output_file(memsize, full_res_three, full_res_four);
            i = 0;
            if (overflow == 0)
                break;
            else if (overflow - 50000 > 0)
            {
                memsize = 50000;
                overflow -= 50000;
            }
            else
            {
                memsize = overflow;
                overflow = 0;
            }
        }
        one_full_switch_tracking(cg, nonInd_three, nonInd_four);
        for (int j = 0; j < 4; j++)
            full_res_three[i][j] = nonInd_three[j];
        for (int j = 0; j < 11; j++)
            full_res_four[i][j] = nonInd_four[j];
        i++;
    }
    auto t_dynamic_end = std::chrono::high_resolution_clock::now();
    auto t_dynamic = std::chrono::duration_cast<std::chrono::nanoseconds>(t_dynamic_end - t_dynamic_begin);
    auto t_full = std::chrono::duration_cast<std::chrono::nanoseconds>(t_dynamic_end - t_profile_begin);
    printf("Track and Count took: %.3f seconds.\n", t_dynamic.count() * 1e-9);
    printf("Full Algorithm took: %.3f seconds.\n", t_full.count() * 1e-9);
}
