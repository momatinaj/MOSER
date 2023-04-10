#ifndef ESCAPE_GETALLCOUNTS_H_
#define ESCAPE_GETALLCOUNTS_H_

#include <algorithm>

#include "Escape/FourVertex.h"
#include "Escape/Utils.h"
#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/Graph.h"
#include "Escape/Conversion.h"
#include "Escape/FiveTrees.h"
#include "Escape/FiveFromCycleClique.h"
#include "Escape/FiveFromTriangles.h"
#include "Escape/AlmostFiveClique.h"
#include "Escape/FiveCycle.h"
#include "Escape/WedgeCollisions.h"
#include "Escape/TriangleProgram.h"

using namespace Escape;

// overloaded functions for 3-vertex and 4-vertex patterns
TriangleInfo getAllThree(CGraph *cg, CDAG *dag, double (&nonInd)[4], bool flag)
{
    double n, m, w;
    TriangleInfo info;

    n = cg->nVertices;
    m = 0;
    w = 0;

    for (VertexIdx i = 0; i < n; i++)
    {
        VertexIdx deg = cg->offsets[i + 1] - cg->offsets[i]; // degree of i
        m = m + deg;
        w = w + (deg * (deg - 1)) / 2; // updating total wedge count
    }
    m = m / 2;

    nonInd[0] = (n * (n - 1) * (n - 2)) / 6; // number of independent sets
    nonInd[1] = m * (n - 2);                 // number of plain edges
    nonInd[2] = w;                           // number of plain wedges

    info = betterWedgeEnumerator(&(dag->outlist));
    nonInd[3] = info.total;
    return info;
}

void getAllFour(CGraph *cg, CDAG *dag, double (&nonInd)[11], TriangleInfo &tri_info)
{
    double n, m, w, t;
    // TriangleInfo tri_info;

    n = cg->nVertices;

    m = 0;
    w = 0;
    // tri_info = betterWedgeEnumerator(&(dag->outlist));
    t = tri_info.total;

    for (VertexIdx i = 0; i < n; i++)
    {
        VertexIdx deg = cg->offsets[i + 1] - cg->offsets[i]; // degree of i
        m = m + deg;
        w = w + (deg * (deg - 1)) / 2; // updating total wedge count
    }
    m = m / 2;

    nonInd[0] = (n * (n - 1) * (n - 2) * (n - 3)) / 24; // number of independent sets
    nonInd[1] = m * ((n - 2) * (n - 3) / 2);            // number of only edges
    nonInd[2] = (m * (m - 1) / 2) - w;                  // number of matchings
    nonInd[3] = w * (n - 3);                            // number of only wedges
    nonInd[4] = t * (n - 3);                            // number of only triangles

    // printf("Getting easy four vertex patterns\n");
    SomeFourPatterns four_info = easyFourCounter(cg, &(dag->outlist));

    // printf("Getting four cycles\n");
    EdgeIdx fourcycles = fourCycleCounter(&(dag->outlist), &(dag->inlist));

    // printf("Getting four cliques\n");
    EdgeIdx fourcliques = fourCliqueCounter(&(dag->outlist));

    nonInd[5] = four_info.threestars;
    nonInd[6] = four_info.threepaths;
    nonInd[7] = four_info.tailedtris;
    nonInd[8] = fourcycles;
    nonInd[9] = four_info.chordalcycles;
    nonInd[10] = fourcliques;
}

// This function generates all non-induced counts for 3-vertex patterns.
// It is a wrapper function that calls the main algorithmic parts, and finally
// calls conversion functions to get induced counts.
//
// Input: pointer to CGraph, corresponding DAG, empty array nonInd with 4 entries
// No output: nonInd will have non-induced counts

void getAllThree(CGraph *cg, CDAG *dag, double (&nonInd)[4])
{
    double n, m, w;
    TriangleInfo info;

    n = cg->nVertices;
    m = 0;
    w = 0;

    for (VertexIdx i = 0; i < n; i++)
    {
        VertexIdx deg = cg->offsets[i + 1] - cg->offsets[i]; // degree of i
        m = m + deg;
        w = w + (deg * (deg - 1)) / 2; // updating total wedge count
    }
    m = m / 2;

    nonInd[0] = (n * (n - 1) * (n - 2)) / 6; // number of independent sets
    nonInd[1] = m * (n - 2);                 // number of plain edges
    nonInd[2] = w;                           // number of plain wedges

    info = betterWedgeEnumerator(&(dag->outlist));
    nonInd[3] = info.total;
}

// This function generates all non-induced counts for up to 4-vertex patterns.
// It is a wrapper function that calls the main algorithmic parts, and finally
// calls conversion functions to get induced counts.
//
// Input: pointer to CGraph, corresponding DAG, empty array nonInd with 11 entries
// No output: nonInd will have non-induced counts

void getAllFour(CGraph *cg, CDAG *dag, double (&nonInd)[11])
{
    double n, m, w, t;
    TriangleInfo tri_info;

    n = cg->nVertices;

    m = 0;
    w = 0;
    tri_info = betterWedgeEnumerator(&(dag->outlist));
    t = tri_info.total;

    for (VertexIdx i = 0; i < n; i++)
    {
        VertexIdx deg = cg->offsets[i + 1] - cg->offsets[i]; // degree of i
        m = m + deg;
        w = w + (deg * (deg - 1)) / 2; // updating total wedge count
    }
    m = m / 2;

    nonInd[0] = (n * (n - 1) * (n - 2) * (n - 3)) / 24; // number of independent sets
    nonInd[1] = m * ((n - 2) * (n - 3) / 2);            // number of only edges
    nonInd[2] = (m * (m - 1) / 2) - w;                  // number of matchings
    nonInd[3] = w * (n - 3);                            // number of only wedges
    nonInd[4] = t * (n - 3);                            // number of only triangles

    printf("Getting easy four vertex patterns\n");
    SomeFourPatterns four_info = easyFourCounter(cg, &(dag->outlist));

    printf("Getting four cycles\n");
    EdgeIdx fourcycles = fourCycleCounter(&(dag->outlist), &(dag->inlist));

    printf("Getting four cliques\n");
    EdgeIdx fourcliques = fourCliqueCounter(&(dag->outlist));

    nonInd[5] = four_info.threestars;
    nonInd[6] = four_info.threepaths;
    nonInd[7] = four_info.tailedtris;
    nonInd[8] = fourcycles;
    nonInd[9] = four_info.chordalcycles;
    nonInd[10] = fourcliques;
}

// This function generates all non-induced counts for 5-vertex patterns.
// It is a wrapper function that calls the main algorithmic parts, and finally
// calls conversion functions to get induced counts.
//
// Input: pointer to CGraph, corresponding DAG, array of non-induced 4-vertex counts, empty array nonInd with 34 entries
// No output: nonInd will have non-induced counts

void getAllFive(CGraph *cg, CDAG *dag, double (&nonIndFour)[11], double (&nonIndFive)[34])
{
    double n, m, w, t;
    NoninducedFourCounts nonIndFourStruct;
    SomeFourPatterns four_info;

    nonIndFourStruct.threestars = nonIndFour[5];
    nonIndFourStruct.threepaths = nonIndFour[6];
    nonIndFourStruct.tailedtris = nonIndFour[7];
    nonIndFourStruct.fourcycles = nonIndFour[8];
    nonIndFourStruct.chordalcycles = nonIndFour[9];
    nonIndFourStruct.fourcliques = nonIndFour[10];

    four_info.threestars = nonIndFourStruct.threestars;
    four_info.threepaths = nonIndFourStruct.threepaths;
    four_info.tailedtris = nonIndFourStruct.tailedtris;
    four_info.chordalcycles = nonIndFourStruct.chordalcycles;

    printf("Getting all triangles\n");
    TriangleInfo tri_info = betterWedgeEnumerator(&(dag->outlist));
    TriangleList allTris = storeAllTriangles(cg, tri_info.total);

    printf("Also getting reverse triangle info\n");
    TriangleInfo in_tri_info = moveOutToIn(&(dag->outlist), &(dag->inlist), &tri_info);

    n = cg->nVertices;

    m = 0;
    w = 0;

    for (VertexIdx i = 0; i < n; i++)
    {
        VertexIdx deg = cg->offsets[i + 1] - cg->offsets[i]; // degree of i
        m = m + deg;
        w = w + (deg * (deg - 1)) / 2; // updating total wedge count
    }
    m = m / 2;

    t = tri_info.total;

    nonIndFive[0] = (n * (n - 1) * (n - 2) * (n - 3) * (n - 4)) / 120;
    nonIndFive[1] = m * ((n - 2) * (n - 3) * (n - 4)) / 6;
    nonIndFive[2] = ((m * (m - 1)) / 2 - w) * (n - 4);
    nonIndFive[3] = w * ((n - 3) * (n - 4)) / 2;
    nonIndFive[4] = t * ((n - 3) * (n - 4)) / 2;

    for (int i = 5; i < 11; i++)
        nonIndFive[i] = nonIndFour[i] * (n - 4);

    nonIndFive[11] = w * (m - 2) - 3 * t - 3 * nonIndFourStruct.threestars - 2 * nonIndFourStruct.threepaths;
    nonIndFive[12] = t * (m - 3) - nonIndFourStruct.tailedtris;

    printf("Counting trees\n");
    FiveTrees tree_counts = fiveTreeCounter(cg, nonIndFourStruct, tri_info.total);

    printf("Counting triangle based patterns\n");
    FiveFromTriangles tri_based_counts = fiveFromTriCounter(cg, &(dag->outlist), &tri_info, four_info);

    Count hourglass = count5_Hourglass<false>(&(dag->outlist), &tri_info, 0);
    Count stingray = count5_Stingray<false>(&(dag->outlist), &(dag->inlist), &tri_info, 0);
    Count three_tri_col = count5_StellateTrident<false>(&(dag->outlist), &tri_info, 0);
    Count tri_strip = count5_TriangleStrip<true>(&(dag->outlist), nonIndFourStruct.fourcliques, &tri_info);
    Count cobra = count5_Cobra<true>(&(dag->outlist), &(dag->inlist), nonIndFourStruct.fourcliques);

    printf("Counting 4-cycle and 4-clique based patterns\n");
    CycleBased cycle_related = fourCycleBasedCounter(&(dag->outlist), &(dag->inlist), &tri_info, &in_tri_info, nonIndFourStruct.chordalcycles);
    CliqueBased clique_related = fourCliqueBasedCounter(cg, &(dag->outlist), &tri_info);

    printf("Counting five cycles\n");
    EdgeIdx five_cycle = fiveCycleCounter(&(dag->outlist), &(dag->inlist));

    printf("Counting collision patterns\n");
    CollisionPatterns collision_vals = fromTriangleList(cg, &allTris);

    printf("Counting almost cliques\n");
    EdgeIdx almost_clique = almostFiveClique(cg);

    nonIndFive[13] = tree_counts.fourstars;
    nonIndFive[14] = tree_counts.prongs;
    nonIndFive[15] = tree_counts.fourpaths;
    nonIndFive[16] = tri_based_counts.forktailedtris;
    nonIndFive[17] = tri_based_counts.longtailedtris;
    nonIndFive[18] = tri_based_counts.doubletailedtris;
    nonIndFive[19] = cycle_related.tailedfourcycles;
    nonIndFive[20] = five_cycle;
    nonIndFive[21] = hourglass;
    nonIndFive[22] = cobra;
    nonIndFive[23] = stingray;
    nonIndFive[24] = cycle_related.hattedfourcycles;
    nonIndFive[25] = collision_vals.threeWedgeCol;
    nonIndFive[26] = three_tri_col;
    nonIndFive[27] = clique_related.tailedfourcliques;
    nonIndFive[28] = tri_strip;
    nonIndFive[29] = collision_vals.chordalWedgeCol;
    nonIndFive[30] = collision_vals.wheel;
    nonIndFive[31] = clique_related.hattedfourcliques;
    nonIndFive[32] = almost_clique;
    nonIndFive[33] = clique_related.fivecliques;
}
#endif
