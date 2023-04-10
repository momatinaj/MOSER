#ifndef ESCAPE_FIVECYCLE_H_
#define ESCAPE_FIVECYCLE_H_

#include <algorithm>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"


using namespace Escape;

EdgeIdx fiveCycleCounter(CGraph *gout, CGraph *gin)
{
   EdgeIdx ret=0;

   VertexIdx i,j,k,ell;

   VertexIdx *wedge_count = new VertexIdx[gin->nVertices+1];

   int DEBUG = 0;

   for (i=0; i < gin->nVertices; i++)
       wedge_count[i] = 0;


   for (i=0; i < gin->nVertices; ++i) // loop over vertices
   {
       if (DEBUG)
           printf("----Handling %lld\n",i);

       // loop over inout wedges ending at i
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives an inout wedge
           {
               k = gin->nbors[next];  // i <- j <- k is wedge
               wedge_count[k]++;
           }

           for (EdgeIdx next = gout->offsets[j]; next < gout->offsets[j+1]; ++next)
           {
               k = gout->nbors[next]; // i <- j -> k is wedge
               if (k==i)
                   continue;
               wedge_count[k]++;
           }
       }

       // loop over inout wedges starting at i
       for (EdgeIdx pos = gout->offsets[i]; pos < gout->offsets[i+1]; ++pos)
       {
           j = gout->nbors[pos];
           for (EdgeIdx next = gout->offsets[j]; next < gout->offsets[j+1]; ++next)
           {
               k = gout->nbors[next]; // i -> j -> k is wedge
               wedge_count[k]++;
           }
       }


       // generate three paths ending at i
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos)
       {
           j = gin->nbors[pos];
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next)
           {
               k = gin->nbors[next]; // i <- j <- k is wedge
               for (EdgeIdx next2 = gout->offsets[k]; next2 < gout->offsets[k+1]; ++next2)
               {
                   ell = gout->nbors[next2]; // i <- j <- k -> ell is three path
                   if (ell == j || ell == i)
                       continue;
                   ret += wedge_count[ell];
                   if (gout->getEdgeBinary(i,k) != -1 || gout->getEdgeBinary(k,i) != -1)
                       ret--;
                   if (gout->getEdgeBinary(j,ell) != -1 || gout->getEdgeBinary(ell,j) != -1)
                       ret--;
                   if (DEBUG)
                   {
                       printf("Three path: %lld <- %lld <- %lld -> %lld\n",i,j,k,ell);
                       printf("Wedge count: %lld\n",wedge_count[ell]);
                       if (gout->getEdgeBinary(i,k) != -1 || gout->getEdgeBinary(k,i) != -1)
                           printf("Subtracted for i,k\n");
                       if (gout->getEdgeBinary(j,ell) != -1 || gout->getEdgeBinary(ell,j) != -1)
                           printf("Subtracted for j,ell\n");
                   }
               }
           }
       }

       // clear wedge_count
       // loop over inout wedges ending at i
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives an inout wedge
           {
               k = gin->nbors[next];  // i <- j <- k is wedge
               wedge_count[k] = 0;
           }

           for (EdgeIdx next = gout->offsets[j]; next < gout->offsets[j+1]; ++next)
           {
               k = gout->nbors[next];
               wedge_count[k] = 0;
           }
       }

       // loop over inout wedges starting at i
       for (EdgeIdx pos = gout->offsets[i]; pos < gout->offsets[i+1]; ++pos)
       {
           j = gout->nbors[pos];
           for (EdgeIdx next = gout->offsets[j]; next < gout->offsets[j+1]; ++next)
           {
               k = gout->nbors[next]; // i -> j -> k is wedge
               wedge_count[k] = 0;
           }
       }
   }

   return ret;
}

/*
EdgeIdx fiveCycleCounter(CGraph *gout, CGraph *gin, EdgeIdx tailedtris)
{
    EdgeIdx initial = rawFiveCycleCounter(gout,gin);
    return initial - tailedtris;
}
*/


#endif



