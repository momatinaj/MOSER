#ifndef ESCAPE_FIVETREES_H_
#define ESCAPE_FIVETREES_H_

#include <algorithm>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/FourVertex.h"
#include "Escape/Utils.h"

using namespace Escape;

// The following structure stores all the trees counts

struct FiveTrees
{
    //trees
    EdgeIdx fourstars;
    EdgeIdx prongs;
    EdgeIdx fourpaths;
};



FiveTrees fiveTreeCounter(CGraph *g, NoninducedFourCounts fourcounts, EdgeIdx triangle_count)
{
    FiveTrees ret;

    ret.fourstars = 0;
    ret.prongs = 0;
    ret.fourpaths = 0;

    // total number of fourstars = \sum_i {d_i \choose 4}
    for (VertexIdx i = 0; i < g->nVertices; i++)
    {
        VertexIdx degi;
        degi = g->offsets[i+1] - g->offsets[i]; //degree of i
        ret.fourstars += (degi*(degi-1)*(degi-2)*(degi-3))/24; // update total number fourstars
    }

    // total number of prongs = \sum_{e=(i,j)} ((d_i-1) \choose 2)*(d_j-1) - 2*#tailed-triangles
    for (VertexIdx i = 0; i < g->nVertices; i++)
    {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i]; // degree of i
        for (EdgeIdx posj = g->offsets[i]; posj < g->offsets[i+1]; posj++) //looping over neighbors of i
        {
            VertexIdx j = g->nbors[posj]; // neighbor j
            VertexIdx degj = g->offsets[j+1] - g->offsets[j]; // degree of j
//             printf("%ld %ld, %ld %ld\n",i,j,degi,degj);
            ret.prongs += ((degi-1)*(degj-1)*(degj-2))/2; // update prong count
        }
    }

    ret.prongs -= 2*fourcounts.tailedtris; // correction for prongs

    // we also set up the wedges array in this loop
    for (VertexIdx i = 0; i < g->nVertices; i++)
    {
        EdgeIdx wedges = 0; // initializing wedge count to 0
        EdgeIdx wedges_sqr = 0; // initializing wedge_sqr to 0
        for (EdgeIdx posj = g->offsets[i]; posj < g->offsets[i+1]; posj++) //looping over neighbors of i
        {
            VertexIdx j = g->nbors[posj]; // neighbor j
            VertexIdx degj = g->offsets[j+1] - g->offsets[j]; // degree of j

            wedges += degj-1; // update number of wedges ending at i
            wedges_sqr += (degj-1)*(degj-1); // update wedge_sqr ending at i
        }

//         printf("%ld %ld: %ld %ld\n",g->nVertices,i,wedges,wedges_sqr);

//         if (i == g->nVertices-1)
//         {
//             printf("ends %ld %ld\n",g->offsets[i],g->offsets[i+1]);
//             for (EdgeIdx posj = g->offsets[i]; posj < g->offsets[i+1]; posj++) //looping over neighbors of i
//             {
//                 VertexIdx j = g->nbors[posj]; // neighbor j
//                 VertexIdx degj = g->offsets[j+1] - g->offsets[j]; // degree of j
// 
//                 printf("pos = %ld, j = %ld, degj = %ld\n",posj,j,degj);    
//                 wedges += degj-1; // update number of wedges ending at i
//                 wedges_sqr += (degj-1)*(degj-1); // update wedge_sqr ending at i
//             }
//         }
// 

        ret.fourpaths += ((wedges*wedges) - wedges_sqr)/2;
    }

    ret.fourpaths -= 4*fourcounts.fourcycles;
    ret.fourpaths -= 2*fourcounts.tailedtris;
    ret.fourpaths -= 3*triangle_count;
    
    return ret;
}

#endif



