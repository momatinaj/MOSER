#ifndef ESCAPE_FIVEFROMCYCLECLIQUE_H_
#define ESCAPE_FIVEFROMCYCLECLIQUE_H_

#include <algorithm>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"


using namespace Escape;

// The following structures are simply lists to store the counts
// of numerous 5-vertex patterns. We break it up into 4-cycle based and 4-clique based patterns
struct CycleBased
{
    EdgeIdx tailedfourcycles;
    EdgeIdx hattedfourcycles;
};

struct CliqueBased
{
    EdgeIdx tailedfourcliques;
    EdgeIdx hattedfourcliques;
    EdgeIdx fivecliques;
};




CycleBased initfourCycleBasedCounter(CGraph *gout, CGraph *gin, TriangleInfo *outinfo, TriangleInfo *ininfo)
{
   CycleBased ret;

   // BUG FIX: Sept 20, 2019
   // Bug found for user Beatrice on bitbucket. The following arrays were initialized to size gout->nEdges.
   //    This is a problem if the number of edges is less than the number of vertices. That never happened in the test data sets earlier,
   //    so we never detected this bug. 
   
   EdgeIdx *outout_count = new VertexIdx[gout->nVertices+1]; // stores the counts of outout wedges from a vertex
   EdgeIdx *inout_count = new VertexIdx[gout->nVertices+1];  // stores the counts of inout wedges from a vertex
   VertexIdx i,j,k;
   VertexIdx degi, degj, degk; 
   EdgeIdx total_tri;
   //EdgeIdx to_add;
   //EdgeIdx to_add_hatted;

   //int always_print = 0;
   int DEBUG = 0;

   //EdgeIdx numType1FourCycle, numType2FourCycle, numType3FourCycle;
   double type1tailed = 0, type2tailed = 0;
   EdgeIdx type1hatted = 0, type2hatted = 0, type3hatted = 0;
   EdgeIdx type3tailed = 0;


   // initializing counts to 0
   for (i=0; i < gin->nVertices; ++i)
   {
       outout_count[i] = 0;
       inout_count[i] = 0;
   }



   for (i=0; i < gin->nVertices; ++i) // loop over vertices
   {
       degi = (gout->offsets[i+1] - gout->offsets[i]) + (gin->offsets[i+1] - gin->offsets[i]);
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gout->offsets[j]; next < gout -> offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = gout->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k==i)   // not a wedge, so ignore
                   continue;
               ++outout_count[k];
           }
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives inout wedge
           {
               k = gin->nbors[next]; // i <- j <- k is inout wedge centered at j
               ++inout_count[k];
           }
       }

       // we loop over all outout and inout wedges again. By keeping track of the edges and triangles
       // incident to each wedge, and using outout_count and inout_count, we can get the final pattern count

       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos)
       {
           j = gin->nbors[pos];
           degj = (gout->offsets[j+1] - gout->offsets[j]) + (gin->offsets[j+1] - gin->offsets[j]);
           for (EdgeIdx next = gout->offsets[j]; next < gout->offsets[j+1]; ++next)
           {
               k = gout->nbors[next]; // i <- j -> k is outout wedge
               if (k==i)
                   continue;
               degk = (gout->offsets[k+1] - gout->offsets[k]) + (gin->offsets[k+1]-gin->offsets[k]);

               // to get type1 4-cycles, we look at outout_count[k]
               // there are outout_count[k]-1 type1 4-cycles with wedge (i,j,k), and we multiply this by the total number of triangles on (i,j,k)
               total_tri = ininfo->perEdge[pos] + outinfo->perEdge[next];
               type1hatted += (outout_count[k]-1)*total_tri;

               // to get type3 4-cycles, we look at inout_out[k]. That is also the number of such 4-cycles with (i,j,k)
               type3hatted += inout_count[k]*total_tri;

               if (k < i)
                    type1tailed += ((double)(degk-2)/2 + (double)(degi-2)/2 + (double)(degj-2))*(outout_count[k]-1);

               type3tailed += (degk-2 + degi-2 + degj-2)*inout_count[k];

if(DEBUG)
{
               printf("Consider %lld <- %lld -> %lld\n",i,j,k);
               printf("Info: outout_count[%lld] = %lld, inout_count[%lld] = %lld, total_tri = %lld\n",k,outout_count[k],k,inout_count[k],total_tri);
               printf("Degs: degi = %lld, degj = %lld, degk = %lld\n",degi,degj,degk);
               printf("Adding %lld to type1hatted, %lld to type3hatted\n",(outout_count[k]-1)*total_tri,inout_count[k]*total_tri);
               printf("Adding %lf to type1tailed, %lld to type3tailed\n\n",((double)(degk-2)/2 + (double)(degi-2)/2 + (double)(degj-2))*(outout_count[k]-1),(degk-2 + degi-2 + degj-2)*inout_count[k]);
}
            }
/*               to_add = (outout_count[k] - 1)*(degj - 2)/2 + inout_count[k]*(degj-2); // (i,j,k) form a type1 4-cycle with every other outout wedge. An edge incident to j gives the tail. This pattern is touched twice, so we divide by 2
               tailedfourcycles += to_add;

               total_tri = ininfo->perEdge[pos] + outinfo->perEdge[next];
               to_add_hatted = (outout_count[k] - 1)*total_tri + inout_count[k]*total_tri;
               hattedfourcycles += to_add_hatted;
               
               if (always_print || to_add > 0)
               {
                    printf("Incrementing tailed %lld in 1st\n",to_add);
                    printf("i = %lld, j = %lld, k = %lld\n",i,j,k);
                    printf("outout_count[k] = %lld, inout_count[k] = %lld, degj = %lld\n\n",outout_count[k], inout_count[k], degj);
               }
               if (always_print || to_add_hatted > 0)
               {
                    printf("Incrementing hatted %lld in 1st\n",to_add_hatted);
                    printf("i = %lld, j = %lld, k = %lld\n",i,j,k);
                    printf("outout_count[k] = %lld, inout_count[k] = %lld, total_tri = %lld\n\n",outout_count[k], inout_count[k], total_tri);
               }*/

           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next)
           {
               k = gin->nbors[next]; // i <- j <- k is inout wedge
               degk = (gout->offsets[k+1] - gout->offsets[k]) + (gin->offsets[k+1]-gin->offsets[k]);
 
               total_tri = ininfo->perEdge[pos] + ininfo->perEdge[next];

               type2hatted += (inout_count[k]-1)*total_tri;
               
               type3hatted += outout_count[k]*total_tri;
               
               type2tailed += ((double)(degk-2)/2 + (double)(degi-2)/2 + degj-2)*(inout_count[k]-1);

               type3tailed += (degj-2)*outout_count[k];
if(DEBUG)
{
               printf("Consider %lld <- %lld <- %lld\n",i,j,k);
               printf("Info: outout_count[%lld] = %lld, inout_count[%lld] = %lld, total_tri = %lld\n",k,outout_count[k],k,inout_count[k],total_tri);
               printf("Degs: degi = %lld, degj = %lld, degk = %lld\n",degi,degj,degk);
               printf("Adding %lld to type2hatted, %lld to type3hatted\n",(inout_count[k]-1)*total_tri,outout_count[k]*total_tri);
               printf("Adding %lf to type2tailed, %lld to type3tailed\n\n",((double)(degk-2)/2 + (double)(degi-2)/2 + degj-2)*(inout_count[k]-1),(degj-2)*outout_count[k]);
}
           } 
       }

       // clearing outout_count and inout_count
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gout->offsets[j]; next < gout -> offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = gout->nbors[next];  // i <- j -> k is outout wedge centered at j
               outout_count[k] = 0;
           }
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives inout wedge
           {
               k = gin->nbors[next];
               inout_count[k] = 0;
           }
       }
   }

//    printf("type1hatted = %lld, type2hatted = %lld, type3hatted = %lld\n",(long)type1hatted, (long)type2hatted, (long)type3hatted);
//    printf("type1tailed = %lld, type2tailed = %lld, type3tailed = %lld\n",(long)type1tailed, (long)type2tailed, (long)type3tailed);

   ret.tailedfourcycles = (long)(type1tailed + type2tailed + type3tailed + 0.5);
   ret.hattedfourcycles = (long)(type1hatted/2 + type2hatted + type3hatted + 0.5);
   return ret;
}

CycleBased fourCycleBasedCounter(CGraph *gout, CGraph *gin, TriangleInfo *outinfo, TriangleInfo *ininfo, EdgeIdx chordalcycles)
{
    CycleBased init = initfourCycleBasedCounter(gout, gin, outinfo, ininfo);
    CycleBased ret;
    
    ret.tailedfourcycles = init.tailedfourcycles - 2*chordalcycles;
    ret.hattedfourcycles = init.hattedfourcycles - 4*chordalcycles;
    return ret;

}

CliqueBased fourCliqueBasedCounter(CGraph *g, CGraph *gout, TriangleInfo *info)
{
    CliqueBased ret; // return value

    ret.tailedfourcliques = 0;
    ret.hattedfourcliques = 0;
    ret.fivecliques = 0;

    bool to_print = false;
   
    int DEBUG = 0; 
    VertexIdx *dummy = new VertexIdx[5];
    VertexIdx *triends = new VertexIdx[gout->nVertices+1]; // array to store triangle ends
    VertexIdx *k4ends = new VertexIdx[gout->nVertices+1]; // array to store 4-clique ends
    EdgeIdx *ind_from_i = new VertexIdx[gout->nVertices+1]; // array to store indices into nbors
    EdgeIdx *ind_from_j = new VertexIdx[gout->nVertices+1]; // array to store indices into nbors

    for (VertexIdx i=0; i < gout->nVertices; ++i) // loop over vertices
        for (VertexIdx posj = gout->offsets[i]; posj < gout->offsets[i+1]; ++posj) // loop over out-neighbors of i
        {
            VertexIdx j = gout->nbors[posj]; // j is current out-neighbor
            VertexIdx count = 0;
            for (VertexIdx posk = posj+1; posk < gout->offsets[i+1]; ++posk) // loop over another out-neighbor of i, that is "ahead" of j in list of out-neighbors
            {
                VertexIdx k = gout->nbors[posk]; // k is next out-neighbor
                EdgeIdx edge_ind_jk = gout->getEdgeBinary(j,k);// check if edge (j,k) is present
                if (edge_ind_jk != -1) // edge (j,k) is present
                {
                    triends[count] = k;  // so (i,j,k) form a triangle. we store the fact that k forms a triangle with edge (i,j) in digraph gout
                    ind_from_i[count] = posk;
                    ind_from_j[count] = edge_ind_jk;
                    ++count;
                }
            }

            for (VertexIdx posk = 0; posk < count; ++posk) // loop over all of triangles with edge (i,j)
            {
                VertexIdx count2 = 0;
                VertexIdx k = triends[posk]; // (i,j,k) is triangle
                for (VertexIdx posell = posk+1; posell < count; ++posell) // loop over another triangle with edge (i,j)
                {
                    VertexIdx ell = triends[posell]; // (i,j,ell) is next triangle
                    EdgeIdx edge_ind_kell = gout->getEdgeBinary(k,ell);
                    if (edge_ind_kell != -1) // (k,ell) is an end, thus (i,j,k,ell) form a 4-clique
                    {
                        k4ends[count2] = ell; // ell forms 4-clique with (i,j,k), so this is stored in k4ends
                        ++count2;

                        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
                        VertexIdx degj = g->offsets[j+1] - g->offsets[j];
                        VertexIdx degk = g->offsets[k+1] - g->offsets[k];
                        VertexIdx degell = g->offsets[ell+1] - g->offsets[ell];

                        VertexIdx total_edge = degi + degj + degk + degell - 12;

                        VertexIdx tri_ij = info->perEdge[posj];
                        VertexIdx tri_ik = info->perEdge[ind_from_i[posk]];
                        VertexIdx tri_iell = info->perEdge[ind_from_i[posell]];
                        VertexIdx tri_jk = info->perEdge[ind_from_j[posk]];
                        VertexIdx tri_jell = info->perEdge[ind_from_j[posell]];
                        VertexIdx tri_kell = info->perEdge[edge_ind_kell];

                        VertexIdx total_tri = tri_ij + tri_ik + tri_iell + tri_jk + tri_jell + tri_kell - 12;

                        ret.tailedfourcliques += total_edge;
                        ret.hattedfourcliques += total_tri;
                    }
                }

                for (VertexIdx posell = 0; posell < count2; ++posell) // loop over all pairs of 4-cliques with triangle (i,j,k)
                    for (VertexIdx posoh = posell+1; posoh < count2; ++posoh)
                    {
                        VertexIdx ell = k4ends[posell]; // ell and oh are the vertices
                        VertexIdx oh = k4ends[posoh];
//                         if (DEBUG)
//                             printf("%lld %lld %lld: checking %lld %lld\n",i,j,k,ell,oh);
                        if (gout->isEdge(ell,oh) != -1) // (ell,oh) is edge, so (i,j,k,ell,oh) form 5-clique
                        {
                            ret.fivecliques++;
                            if (DEBUG)
                            {
                                dummy[0] = i;
                                dummy[1] = j;
                                dummy[2] = k;
                                dummy[3] = ell;
                                dummy[4] = oh;
                                std::sort(dummy,dummy+5);
                                printf("%lld %lld %lld %lld %lld\n",dummy[0],dummy[1],dummy[2],dummy[3],dummy[4]);
                            }
                            if(ret.fivecliques%100000000 == 0 && to_print)
                                printf("Found another 100M 5-cliques\n");
                        }
                    }
            }
        }
    return ret;
}

#endif



