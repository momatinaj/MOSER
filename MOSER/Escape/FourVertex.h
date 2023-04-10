#ifndef ESCAPE_FOURVERTEX_H_
#define ESCAPE_FOURVERTEX_H_

#include <algorithm>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/Utils.h"

using namespace Escape;

// The following structures are simply lists to store the counts
// of numerous 4-vertex patterns. Nothing fancy here.
struct InducedFourCounts
{
    EdgeIdx threestars;
    EdgeIdx threepaths;
    EdgeIdx tailedtris;
    EdgeIdx fourcycles;
    EdgeIdx chordalcycles;
    EdgeIdx fourcliques;
};

struct NoninducedFourCounts
{
    EdgeIdx threestars;
    EdgeIdx threepaths;
    EdgeIdx tailedtris;
    EdgeIdx fourcycles;
    EdgeIdx chordalcycles;
    EdgeIdx fourcliques;
};

struct SomeFourPatterns
{
    EdgeIdx threestars;
    EdgeIdx threepaths;
    EdgeIdx tailedtris;
    EdgeIdx chordalcycles;
};


// Some of the four vertex patterns can be counted by simple combinatorial arguments.
// Let d_v be the degree of vertex v, t_v be the number of triangles incident. 
// For edge e, define t_e to be the number of triangles incident to e.
// The following formulae are for non-induced counts.
// 
// #3-stars = \sum_v {d_v \choose 3}
// #3-paths = \sum_{e = (u,v)} (d_u - 1)(d_v - 1) - 3*#triangles
// #tailed-triangles = \sum_v (d_v-2)t_v
// #chordal-cycles = \sum_e {t_e \choose 2}
//
// 

SomeFourPatterns easyFourCounter(CGraph *g, CGraph *gout)
{
    SomeFourPatterns ret;
    ret.threestars = 0;
    ret.threepaths = 0;
    ret.chordalcycles = 0;
    ret.tailedtris = 0;
    EdgeIdx degi=0, degj=0;

    TriangleInfo info;
    info = betterWedgeEnumerator(gout);
    
    for (VertexIdx i=0; i < g->nVertices; i++)
    {
        degi = g->offsets[i+1]-g->offsets[i];
        ret.threestars += degi*(degi-1)*(degi-2)/6;
        ret.tailedtris += (degi-2)*info.perVertex[i];
    }


    for(VertexIdx i=0; i < gout->nVertices; ++i)
        for(EdgeIdx posj = gout->offsets[i]; posj < gout->offsets[i+1]; ++posj)
        {
            VertexIdx j = gout->nbors[posj];
            degi = g->offsets[i+1] - g->offsets[i];
            degj = g->offsets[j+1] - g->offsets[j];
            ret.threepaths += (degi-1)*(degj-1);

            ret.chordalcycles += info.perEdge[posj]*(info.perEdge[posj]-1)/2;
        }
        
    ret.threepaths -= 3*info.total;
    
    return ret;
}


// Old version uses old wedge enumerator
// 

SomeFourPatterns oldeasyFourCounter(CGraph *g, CGraph *gout)
{
    SomeFourPatterns ret;
    ret.threestars = 0;
    ret.threepaths = 0;
    ret.chordalcycles = 0;
    ret.tailedtris = 0;
    EdgeIdx degi=0, degj=0;

    TriangleInfo info;
    info = wedgeEnumerator(gout);
    
    for (VertexIdx i=0; i < g->nVertices; i++)
    {
        degi = g->offsets[i+1]-g->offsets[i];
        ret.threestars += degi*(degi-1)*(degi-2)/6;
        ret.tailedtris += (degi-2)*info.perVertex[i];
    }


    for(VertexIdx i=0; i < gout->nVertices; ++i)
        for(EdgeIdx posj = gout->offsets[i]; posj < gout->offsets[i+1]; ++posj)
        {
            VertexIdx j = gout->nbors[posj];
            degi = g->offsets[i+1] - g->offsets[i];
            degj = g->offsets[j+1] - g->offsets[j];
            ret.threepaths += (degi-1)*(degj-1);

            ret.chordalcycles += info.perEdge[posj]*(info.perEdge[posj]-1)/2;
        }
        
    ret.threepaths -= 3*info.total;
    
    return ret;
}

EdgeIdx fourCycleCounter(CGraph *gout, CGraph *gin)
{
   EdgeIdx ret = 0;

   VertexIdx i,j,k;

   VertexIdx *wedge_count = new VertexIdx[gout->nVertices+1]; //stores number of wedges ending at a vertex
   for (i=0; i < gout->nVertices; i++) //initialize all wedge_count values to 0
       wedge_count[i] = 0;

   for (i=0; i < gin->nVertices; ++i) // loop over vertices
   {
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gout->offsets[j]; next < gout -> offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = gout->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k>=i)   // break ties to prevent overcount
                   continue;
               wedge_count[k]++; // increment number of wedges ending at k
           }
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives inout wedge
           {
               k = gin->nbors[next]; // i <- j <- k is inout wedge centered at j
               wedge_count[k]++; // increment number of wedges ending at k
           }
       }

       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gout->offsets[j]; next < gout -> offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = gout->nbors[next];  // i <- j -> k is outout wedge centered at j
//                printf("outout %ld %ld, count %ld\n",i,k,wedge_count[k]);
               ret += (wedge_count[k]*(wedge_count[k]-1))/2; // every pair of wedges ending at k yields a four-cycle
               wedge_count[k] = 0; //reset value of wedge_count
           }
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives inout wedge
           {
               k = gin->nbors[next]; // i <- j <- k is inout wedge centered at j
//                printf("inin %ld %ld, count %ld\n",i,k,wedge_count[k]);
               ret += (wedge_count[k]*(wedge_count[k]-1))/2; // every pair of wedges ending at k yields a four-cycle
               wedge_count[k] = 0; //reset value of wedge_count
           }
       }
   }
   return ret;
}

EdgeIdx oldfourCycleCounter(CGraph *gout, CGraph *gin)
{
   EdgeIdx ret = 0;

   EdgeIdx type1 = 0; // initializing counts
   EdgeIdx type2 = 0;
   EdgeIdx type3 = 0;
   EdgeIdx *outout = new EdgeIdx[gout->nEdges+1]; // stores the ends of outout wedges from a vertex
   EdgeIdx *inout = new EdgeIdx[gout->nEdges+1];  // stores the ends of inout wedges from a vertex
   EdgeIdx cur_outout;
   EdgeIdx cur_inout;
   VertexIdx i,j,k;
   EdgeIdx it, it2;
   VertexIdx *outout_runs = new VertexIdx[gout->nVertices+1]; // stores "runs" of outout wedges from a vertex with the same other end
   VertexIdx *outout_ends = new VertexIdx[gout->nVertices+1]; // stores "runs" of inout wedges from a vertex with the same other end

   for (i=0; i < gin->nVertices; ++i) // loop over vertices
   {
       cur_outout = 0;
       cur_inout = 0;
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gout->offsets[j]; next < gout -> offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = gout->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k==i)   // not a wedge, so ignore
                   continue;
               outout[cur_outout] = k; // store k as an end of outout wedge
               ++cur_outout; // increment position in outout[..] array
           }
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives inout wedge
           {
               k = gin->nbors[next]; // i <- j <- k is inout wedge centered at j
               inout[cur_inout] = k; // store k as end of inout wedge
               ++cur_inout; // increment position in inout[..] array
           }
       }

       // at this point, outout contains a list of other end of all outout wedges with one end at i.
       // By sorting outout, we can count the number of outout wedges (with one end at i) that end at some k.
       // (analogously for inout)

       std::sort(outout,outout+cur_outout);
       std::sort(inout,inout+cur_inout);

       // having sorted outout and inout, we determine how many times some vertex k appears in each of
       // these arrays. Such a contiguous sequence of k's in (say) outout is called a "run"

       outout[cur_outout] = -1;   // placing a dummy end at outout and inout. -1 is not a vertex, so this leads to a benign termination of loops
       inout[cur_inout] = -1;

       // we will determine the number of times prev appears in outout. observe that any time a vertex appears twice in outout, we get a type1 4-cycle.
       // so prev initializes to outout[0], and we start the run length at 0.

       EdgeIdx prev = outout[0]; 
       EdgeIdx run = 0;
       EdgeIdx counter = 0;

       for (it = 1; it <= cur_outout; ++it) // iterate over outout array
       {
           ++run; // increment run length. at this point, we have already seen this run length (even before looking at outout[it])
           if (outout[it] == prev)  // current element outout[it] is also prev, we need to increse the run. so we continue to next element.
               type1 += run;        // we update type1 wedges. The wedge represented by outout[it] forms a type1 4-cycle with every other wedge in this run so far
           else     // outout[it] != prev, so the run of prev has ended
           {
               outout_ends[counter] = prev;  // we store this information in outout_ends and outout_runs
               outout_runs[counter] = run;
               ++counter;   // just updating counter for these array
               run = 0;     // reinitialize run and prev appropriately
               prev = outout[it];
           }
       }

       // the next part is a little trickier. we want to count both type2 and type3 4-cycles. 
       // type2 4-cycles are formed by collisions between two inout wedges, so counting is identical to type1, simply using the array inout.

       // type3 4-cycles are formed by a collision of an outout wedge and inwedge. So we have to compare run lengths of outout and inout.
       // we stored the run-lengths in outout_ends (storing the vertices) and outout_runs (storing the actual lengths). 

       // initialize as before
       run = 0;
       prev = inout[0];
       it2 = 0; // this is to loops over outout_ends
       for (it = 1; it <= cur_inout; ++it) // iterate over inout array
       {
           ++run; // increment run length. at this point, we have already seen this run length (even before looking at inout[it])
           if (inout[it] == prev) // current element is also prev, so run should increments. So continue.
               type2 += run;     // update type2 wedges, since current wedge forms type2 4-cycle with everything in run so far
           else   // inout[it] != prev, so run has ended
           {
               while(outout_ends[it2] < prev && it2 < counter)  // loop over outout_ends to find if prev exists as end of outout wedge
                   ++it2;
               
               if (outout_ends[it2] == prev && it2 < counter)   // we found prev in outout_ends, so we update type3 counts
                   type3 += run*outout_runs[it2];   // this is product of run length for prev in outout and inout
               run = 0; // run of prev as ended, so reinitialize
               prev = inout[it];
           }
       }
   }
   type1 = type1/2; // every type1 4-cycle was counted twice
   ret = type1 + type2 + type3;
//    printf("type1 = %ld, type2 = %ld, type3 = %ld\n",type1,type2,type3);
   return ret;
}

// Four-clique counter
// Input: The out-DAG of a graph sorted by ID (It is *critical* that DAG is sorted by ID, not be degree.)
// Output: The number of 4-cliques

EdgeIdx fourCliqueCounter(CGraph *gout)
{
    EdgeIdx ret = 0; // return value
    VertexIdx *triends = new VertexIdx[gout->nVertices+1]; // array to store triangle ends

    for (VertexIdx i=0; i < gout->nVertices; ++i) // loop over vertices
        for (VertexIdx posj = gout->offsets[i]; posj < gout->offsets[i+1]; ++posj) // loop over out-neighbors of i
        {
            VertexIdx j = gout->nbors[posj]; // j is current out-neighbor
            VertexIdx count = 0;
            for (VertexIdx posk = posj+1; posk < gout->offsets[i+1]; ++posk) // loop over another out-neighbor of i, that is "ahead" of j in list of out-neighbors
            {
                VertexIdx k = gout->nbors[posk]; // k is next out-neighbor
//                 printf("looking at tri %ld %ld %ld\n",i,j,k);
                if (gout->isEdgeBinary(j,k)) // check if edge (j,k) is present
                {
                    triends[count] = k;  // so (i,j,k) form a triangle. we store the fact that k forms a triangle with edge (i,j) in digraph gout
                    ++count;
                }
            }

            for (VertexIdx posk = 0; posk < count; ++posk) // loop over all pairs of triangles formed by (i,j)
            {
                VertexIdx k = triends[posk]; // k is vertex as index posk in triends
                VertexIdx degk = gout->offsets[k+1] - gout->offsets[k]; // gettting degree of k in gout
                VertexIdx remaining = count-posk; // number of vertices that k needs to be checked with

//                 printf("looking at %ld %ld %ld\n",i,j,k);
                if (degk >= remaining)
                {   
                    // We will search all other vertices in triends in k's adj list
                    for (VertexIdx posell = posk+1; posell < count; ++posell)
                    {
                        VertexIdx ell = triends[posell]; 
//                         printf("Going to triends %ld %ld %ld %ld\n",i,j,k,ell);
                        if (gout->isEdgeBinary(k,ell)) // (k,ell) is an end, thus (i,j,k,ell) form a 4-clique
                        {
                            ++ret;
                        }
                    }
                }
                else
                {
                    // We will search all vertices in k's adj list in the remaining portion of triends
                    for (EdgeIdx posell = gout->offsets[k]; posell < gout->offsets[k+1]; posell++)
                    {
                        VertexIdx ell = gout->nbors[posell];
                        if (binarySearch(triends+posk+1,count-posk-1,ell) != -1)
                        {
                            ++ret;
                        }
                    }
                }
            }
        }
    return ret;
}

// Old four-clique counter uses linear search
// Also there is no optimization in finding pairs of triends that form edges (which finally forms the 4-clique)
EdgeIdx oldfourCliqueCounter(CGraph *gout)
{
    EdgeIdx ret = 0; // return value
    VertexIdx *triends = new VertexIdx[gout->nVertices+1]; // array to store triangle ends

    for (VertexIdx i=0; i < gout->nVertices; ++i) // loop over vertices
        for (VertexIdx posj = gout->offsets[i]; posj < gout->offsets[i+1]; ++posj) // loop over out-neighbors of i
        {
            VertexIdx j = gout->nbors[posj]; // j is current out-neighbor
            VertexIdx count = 0;
            for (VertexIdx posk = posj+1; posk < gout->offsets[i+1]; ++posk) // loop over another out-neighbor of i, that is "ahead" of j in list of out-neighbors
            {
                VertexIdx k = gout->nbors[posk]; // k is next out-neighbor
                if (gout->isEdge(j,k) != -1) // check if edge (j,k) is present
                {
                    triends[count] = k;  // so (i,j,k) form a triangle. we store the fact that k forms a triangle with edge (i,j) in digraph gout
                    ++count;
                }
            }

            for (VertexIdx posk = 0; posk < count; ++posk) // loop over all pairs of triangles formed by (i,j)
                for (VertexIdx posell = posk+1; posell < count; ++posell)
                {
                    VertexIdx k = triends[posk]; // k and ell are the vertices 
                    VertexIdx ell = triends[posell];
                    if (gout->isEdge(k,ell) != -1) // (k,ell) is an end, thus (i,j,k,ell) form a 4-clique
                    {
                        ++ret;
                        //if (ret%10000000 == 0)
                        //    printf("Found another 10M\n");
                    }
                }
        }
    return ret;
}


#endif



