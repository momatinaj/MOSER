#ifndef ESCAPE_WEDGECOLLISIONS_H_
#define ESCAPE_WEDGECOLLISIONS_H_

#include <algorithm>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/Utils.h"

using namespace Escape;



struct CollisionPatterns
{
    EdgeIdx threeWedgeCol;
    EdgeIdx chordalWedgeCol;
    EdgeIdx wheel;
};

EdgeIdx threeWedgeCollision(CGraph *g)
{
   EdgeIdx ret=0;

   VertexIdx i,j,k;

   VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];

   int DEBUG = 0;

   for (i=0; i < g->nVertices; i++)
       wedge_count[i] = 0;


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
       if (DEBUG)
           printf("----Handling %lld\n",i);

       // loop over wedges to populate the ends of the wedges
       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               wedge_count[k]++;
           }
       }

       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               ret += wedge_count[k]*(wedge_count[k]-1)*(wedge_count[k]-2)/6;
               wedge_count[k] = 0;
           }
       }
   }
   return ret;
}

void generateSquare(CGraph *g, int thresh)
{
   VertexIdx i,j,k;

   VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];

   EdgeIdx nnz = 0;

   int DEBUG = 0;

   for (i=0; i < g->nVertices; i++)
       wedge_count[i] = 0;


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
       if (DEBUG)
           printf("----Handling %lld\n",i);

       // loop over wedges to populate the ends of the wedges
       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               wedge_count[k]++;
           }
       }

       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               if (wedge_count[k] >= thresh)
                   nnz++;
               wedge_count[k] = 0;
           }
       }
   }

   VertexIdx *sqr_nbors = new VertexIdx[nnz+1];
   VertexIdx *sqr_offsets = new VertexIdx[g->nVertices+1];
   VertexIdx *sqr_vals = new VertexIdx[nnz+1];

   EdgeIdx current = 0;

   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
       sqr_offsets[i] = current;
       // loop over wedges to populate the ends of the wedges
       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               wedge_count[k]++;
           }
       }

       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               if (wedge_count[k] >= thresh)
               {
                   sqr_nbors[current] = k;
                   sqr_vals[current] = wedge_count[k];
                   current++;
               }
               wedge_count[k] = 0;
           }
       }
   }
}

EdgeIdx squareDetails(CGraph *g)
{
   EdgeIdx ret=0;

   VertexIdx i,j,k;

   VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];

   int to_collect = 50;
   EdgeIdx nnz[to_collect+1], partial[to_collect+1];

   int DEBUG = 0;

   for (i=0; i < g->nVertices; i++)
       wedge_count[i] = 0;

   for (i=0; i <= to_collect; i++)
       nnz[i]=0;


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
       if (DEBUG)
           printf("----Handling %lld\n",i);

       // loop over wedges to populate the ends of the wedges
       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               wedge_count[k]++;
           }
       }

       for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = g->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
               if (k <= i) // k is higher in the order, so ignore wedge
                   continue;
               if (wedge_count[k] < to_collect)
                   nnz[wedge_count[k]]++;
               else
                   nnz[to_collect]++;
               ret += wedge_count[k]*(wedge_count[k]-1)*(wedge_count[k]-2)/6;
               wedge_count[k] = 0;
           }
       }
   }

   partial[to_collect] = nnz[to_collect];
   printf("partial[%d] = %lld\n",to_collect,partial[to_collect]);
   for (i=to_collect-1; i > 0; i--)
   {
       partial[i] = partial[i+1]+nnz[i];
       printf("partial[%lld] = %lld\n",i,partial[i]);
   }
   return ret;
}



EdgeIdx threeWedgeFromFourCycle(CGraph *gout, CGraph *gin)
{
   EdgeIdx ret=0;

   VertexIdx i,j,k;

   Pair *wedge_ends = new Pair[2*gout->nEdges];
   Pair *diagonals = new Pair[10*gout->nEdges];

   VertexIdx len;

   int DEBUG = 0;


   for (i=0; i < gin->nVertices; ++i) // loop over vertices
   {
       VertexIdx degi = (gout->offsets[i+1] - gout->offsets[i]) + (gin->offsets[i+1] - gin->offsets[i]);
       len = 0;

       if (DEBUG)
           printf("----Handling %lld\n",i);

       // loop over wedges to populate the ends of the wedges
       for (EdgeIdx pos = gin->offsets[i]; pos < gin->offsets[i+1]; ++pos) // loop over in-neighbors of i
       {
           j = gin->nbors[pos]; // j is current in-neighbor
           for (EdgeIdx next = gout->offsets[j]; next < gout -> offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
           {
               k = gout->nbors[next];  // i <- j -> k is outout wedge centered at j
               VertexIdx degk = (gout->offsets[k+1] - gout->offsets[k]) + (gin->offsets[k+1] - gin->offsets[k]);
               if (degk > degi || (degk == degi && k >= i)) // k is higher in the order, so ignore wedge
                   continue;
               wedge_ends[len].first = j;
               wedge_ends[len].second = k;
               ++len;
           }
           for (EdgeIdx next = gin->offsets[j]; next < gin->offsets[j+1]; ++next) // loop over in-neighbors of j, note this gives inout wedge
           {
               k = gin->nbors[next]; // i <- j <- k is inout wedge centered at j
               wedge_ends[len].first = j;
               wedge_ends[len].second = k;
               ++len;
           }
       }

       // sorting wedge_ends by second (that is k) and then by first (that is j)
       std::sort(wedge_ends,wedge_ends+len,pairCompareSecond);
       // now, all wedges ending at k appear contiguously

       //putting dummy at the end
       wedge_ends[len].first = -1;
       wedge_ends[len].second = -1;

       if (DEBUG)
       {
           printf("\nInitial: len = %lld\n",len);
           for (VertexIdx ind = 0; ind < len; ind++)
           {
               printf("%lld %lld\n",wedge_ends[ind].first,wedge_ends[ind].second);
           }
       }

       // we now compute run lengths by the second element in wedge_ends
       // each run provides exactly the length choose 3 new patterns

       //initialize
       VertexIdx run = 1;
       VertexIdx prev = 0;

       for (VertexIdx ind = 1; ind < len; ind++) //looping over wedge_ends
       {
           if(wedge_ends[prev].second == wedge_ends[ind].second) // if previous wedge_end and current (ind) wedge_end have same value of k, run continues
               run++;
           else  // a contiguous sequence of same k values has ended
           {
               if (DEBUG)
                   printf("Run for %lld is %lld\n",prev,run);
               ret += run*(run-1)*(run-2)/6;   // updating pattern count
               prev = ind;    // reinitialize of next run
               run = 1;
           }
       }
       ret += run*(run-1)*(run-2)/6;    // for the last run

       // now we generate four-cycles and populate the diagonals
       VertexIdx current = 0;
       for (VertexIdx ind1 = 0; ind1 < len; ind1++)
           for (VertexIdx ind2 = ind1+1; ind2 < len; ind2++)  // double loop over all wedges. note that ind2 > ind1
           {
               if (wedge_ends[ind1].second != wedge_ends[ind2].second) // if the wedges at ind1 and ind2 have different second values, they cannot generate the same four-cycle. since wedge_ends is sorted by second value, we have exhausted possible four-cycles with wedge_ends[ind1]. so break the loop
                   break;

               // we have found a 4-cycle, which we store in diagonals. we first sort the diagonal to maintain consistency
               if (wedge_ends[ind1].first < wedge_ends[ind2].first) 
               {
                   diagonals[current].first = wedge_ends[ind1].first;       //diagonals stores the opposite ends of the 4-cycles, with maintaining first < second
                   diagonals[current].second = wedge_ends[ind2].first;
               }
               else
               {
                   diagonals[current].second = wedge_ends[ind1].first;
                   diagonals[current].first = wedge_ends[ind2].first;
               }
               current++;
           }

       VertexIdx total_cycles = current;
       std::sort(diagonals,diagonals+current,pairCompareSecond);    // we sort the diagonals, to get all same values together

       if (DEBUG)
       {
           printf("\nDiagonals: total = %lld\n",total_cycles);
           for (VertexIdx ind = 0; ind < total_cycles; ind++)
           {
               printf("%lld %lld\n",diagonals[ind].first,diagonals[ind].second);
           }
       }

       // for each run of the same diagonal, we update the pattern count by run length choose 2

       // initialize run
       run = 1;
       prev = 0;

       for (VertexIdx ind = 1; ind < total_cycles; ind++)
       {
           if((diagonals[ind].first == diagonals[prev].first) && (diagonals[ind].second == diagonals[prev].second)) // so the 4-cycle at ind and at prev are the same. increment run length
               run++;
           else     // run has ended, so update pattern count and reinitialize run
           {
               if (DEBUG)
                   printf("Run for %lld is %lld\n",prev,run);
               ret += run*(run-1)/2;
               prev = ind;
               run = 1;
           }
       }
       ret += run*(run-1)/2;
               
   }

   return ret;
}

CollisionPatterns newWedgeBasedCollision(CGraph *g, CGraph *gout)
{
   CollisionPatterns ret;

   ret.threeWedgeCol = 0;
   ret.chordalWedgeCol = 0;
   int size = 20;

   VertexIdx i,j,k;

   VertexIdx **wedge = new VertexIdx*[g->nVertices+1];
   VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];

   int DEBUG = 0;

   for (i=0; i < g->nVertices; i++)
   {
       wedge_count[i] = 0;
       wedge[i] = new VertexIdx[size];
   }


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
//         VertexIdx degi = g->offsets[i+1] - g->offsets[i];
        if (DEBUG)
            printf("----Handling %lld\n",i);
 

        // loop over wedges to populate the ends of the wedges
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
//                 VertexIdx degk = g->offsets[k+1] - g->offsets[k];
//                 if (degk < degi || ((degk == degi) && k <= i)) // k is lower in degree ordering, so ignore wedge
//                     continue;
                if (wedge_count[k] < size)
                {
                    wedge[k][wedge_count[k]] = j;
                    wedge_count[k]++;
                }
            }
        } 

        // process each four-cycle to get chordal-4-cycles, and update pattern count 
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j

                for (VertexIdx ind1 = 0; ind1 < wedge_count[k]; ind1++)
                {
                    for (VertexIdx ind2 = ind1+1; ind2 < wedge_count[k]; ind2++)
                    {
                        if (g->getEdgeBinary(wedge[k][ind1],wedge[k][ind2]) != -1)
                            ret.chordalWedgeCol++;
                    }
                }

                //clearing memory and reinitializing values
                wedge_count[k] = 0;
            }
        }
   }
   return ret;
}

CollisionPatterns wedgeBasedCollision(CGraph *g, CGraph *gout)
{
   CollisionPatterns ret;

   ret.threeWedgeCol = 0;
   ret.chordalWedgeCol = 0;

   VertexIdx i,j,k;

   VertexIdx **wedge = new VertexIdx*[g->nVertices+1];
   VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];
   VertexIdx *total = new VertexIdx[g->nVertices+1];

   bool *mem_end_seen = new bool[g->nVertices+1];

   int DEBUG = 0;

   for (i=0; i < g->nVertices; i++)
   {
       wedge_count[i] = 0;
       mem_end_seen[i] = false;
       total[i] = 0;
   }


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
        if (DEBUG)
            printf("----Handling %lld\n",i);
 

        // loop over wedges to populate the ends of the wedges
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
//                 VertexIdx degk = g->offsets[k+1] - g->offsets[k];
//                 if (degk < degi || ((degk == degi) && k <= i)) // k is lower in degree ordering, so ignore wedge
//                     continue;
                wedge_count[k]++;
            }
        } 

        // allocating memory for wedge[k], for all k
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
                if (!mem_end_seen[k])
                {
                    wedge[k] = new VertexIdx[wedge_count[k]+1];
                    mem_end_seen[k] = true;
                }
            }
        }
        
        // populate wedge[k], for all k
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
                VertexIdx degk = g->offsets[k+1] - g->offsets[k];
                if (degk < degi || ((degk == degi) && k <= i)) // k is lower in degree ordering, so ignore wedge
                    continue;
                wedge[k][total[k]] = j;
                total[k]++;
            }
        }
         
        // process each four-cycle to get chordal-4-cycles, and update pattern count 
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
                if (wedge_count[k] == 0)
                    continue;
                

//                 if (wedge_count[k] < 10)
//                 {
//                 for (VertexIdx ind1 = 0; ind1 < wedge_count[k]; ind1++)
//                 {
//                     for (VertexIdx ind2 = ind1+1; ind2 < wedge_count[k]; ind2++)
//                     {
//                         ret.chordalWedgeCol++;
//                     }
//                 }
//                 }
// 
                //clearing memory and reinitializing values
                delete[] wedge[k];
                wedge_count[k] = 0;
                mem_end_seen[k] = false;
                total[k] = 0;
            }
        }
   }
   return ret;
}

CollisionPatterns oldWedgeBasedCollision(CGraph *g, CGraph *gout)
{
   CollisionPatterns ret;

   ret.threeWedgeCol = 0;
   ret.chordalWedgeCol = 0;

   VertexIdx i,j,k;

   Pair *wedge_ends = new Pair[10*gout->nEdges];

   VertexIdx *candidates = new VertexIdx[g->nVertices+1]; 

   EdgeIdx cand_count,count; 

   int DEBUG = 0;

//    for (i=0; i < g->nVertices; i++)
//        wedge_count[i] = 0;


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
        if (DEBUG)
            printf("----Handling %lld\n",i);
 
        count = 0;
        // loop over wedges to populate the ends of the wedges
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
                VertexIdx degk = g->offsets[k+1] - g->offsets[k];
                if (degk < degi || ((degk == degi) && k <= i)) // k is lower in degree ordering, so ignore wedge
                    continue;
                wedge_ends[count].first = j;
                wedge_ends[count].second = k;
                count++;
            }
        } 

        std::sort(wedge_ends,wedge_ends+count,pairCompareSecond);
       
        if (DEBUG)
        {
            printf("Wedge_ends: ");
            for (VertexIdx ind = 0; ind < count; ind++)
            {
                printf("(%lld, %lld) ",wedge_ends[ind].first,wedge_ends[ind].second);
            }
            printf("\n-----------\n");
        }
        
        VertexIdx ptr = 0; 
        while(ptr < count)
        {
            VertexIdx current_k = wedge_ends[ptr].second;
            cand_count = 0;
            while(ptr < count && (wedge_ends[ptr].second == current_k))
            {
                candidates[cand_count] = wedge_ends[ptr].first;
                cand_count++;
                ptr++;
            }
            if (DEBUG)
                printf("Looking at %lld %lld, candidates = %lld, at ptr = %lld\n",i,current_k,cand_count,ptr);

            if (cand_count < 3)
                continue;

            ret.threeWedgeCol += (cand_count)*(cand_count-1)*(cand_count-2)/6;

            for (EdgeIdx posj1 = 0; posj1 < cand_count; posj1++)
            {
                VertexIdx j1 = candidates[posj1];
                VertexIdx outdegj1 = gout->offsets[j1+1] - gout->offsets[j1];
                if (cand_count < outdegj1)
                {
                    for (EdgeIdx posj2 = 0; posj2 < cand_count; posj2++)
                    {
                        VertexIdx j2 = candidates[posj2];
                        if (g->getEdgeBinary(j1,j2) != -1)
                        {
                            ret.chordalWedgeCol += cand_count - 2;
                        }
                    }
                }
                else
                {
                    for (EdgeIdx posj2 = gout->offsets[j1]; posj2 < gout->offsets[j1+1]; posj2++)
                    {
                        VertexIdx j2 = gout->nbors[posj2];
                        if (binarySearch(candidates,cand_count,j2) != -1)
                        {
                            ret.chordalWedgeCol += cand_count - 2;
                        }
                    }
                }
            }

        }
   }
   return ret;
}

CollisionPatterns chordalCycleWedgeCollision(CGraph *g, CGraph *gout)
{
   CollisionPatterns ret;

   ret.chordalWedgeCol = 0;
   ret.threeWedgeCol = 0;

   VertexIdx i,j,k;

   VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];
   VertexIdx *triangles = new VertexIdx[g->nVertices+1];
   
   VertexIdx count; 

   int DEBUG = 0;

   for (i=0; i < g->nVertices; i++)
       wedge_count[i] = 0;


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
        if (DEBUG)
            printf("----Handling %lld\n",i);
 
        // loop over wedges to populate the ends of the wedges
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
                VertexIdx degk = g->offsets[k+1] - g->offsets[k];
                if ((degk < degi) || (degk==degi && k <= i)) // k is higher in the order, so ignore wedge
                    continue;
                wedge_count[k]++;
            }
        } 
         

        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos)
        {
            count = 0;
            j = g->nbors[pos];
            for (EdgeIdx next = gout->offsets[j]; next < gout->offsets[j+1]; ++next)
            {
                k = gout->nbors[next];
                if (g->getEdgeBinary(k,i) != -1)
                {
                    triangles[count] = k;
                    count++;
                }
            }

            for (VertexIdx ptr = 0; ptr < count; ++ptr)
            {
                k = triangles[ptr];

                if (DEBUG)
                    printf("Tri %lld %lld %lld\n",i,j,k);

                for (EdgeIdx ind = g->offsets[j]; ind < g->offsets[j+1]; ++ind)
                {
                    VertexIdx ell = g->nbors[ind];
                    VertexIdx degell = g->offsets[ell+1] - g->offsets[ell];
                    if ((degell < degi) || (degell==degi && ell <= i))
                        continue;
                    if (DEBUG)
                        printf("Checking %lld\n",ell);
                    if (g->isEdgeBinary(ell,k))
                    {
                        if (DEBUG)
                            printf("%lld %lld %lld %lld: count = %lld\n",i,j,k,ell,wedge_count[ell]);
                        ret.chordalWedgeCol += wedge_count[ell]-2;
                    }
                }
            }

        }

        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
                ret.threeWedgeCol += (wedge_count[k])*(wedge_count[k]-1)*(wedge_count[k]-2)/6;
                wedge_count[k] = 0;
            }
        }
   }
   return ret;
}


CollisionPatterns fromTriangleList(CGraph *g, TriangleList *allTris)
{
   CollisionPatterns ret;

   ret.chordalWedgeCol = 0;
   ret.threeWedgeCol = 0;
   ret.wheel = 0;

   VertexIdx i,j,k,ell;

   VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];
  
   VertexIdx *diamond_count = new VertexIdx[g->nVertices+1];
   
   EdgeIdx low, mid, high;

   int DEBUG = 0;

   for (i=0; i < g->nVertices; i++)
   {
       wedge_count[i] = 0;
       diamond_count[i] = 0;
   }


   for (i=0; i < g->nVertices; ++i) // loop over vertices
   {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
 
        // loop over wedges to populate the ends of the wedges
        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i -- j -- k is wedge centered at j
                if (k <= i)
                    continue;
                wedge_count[k]++;
            }
        } 
         

        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos)
        {
            j = g->nbors[pos];
            VertexIdx degj = g->offsets[j+1] - g->offsets[j];
            EdgeIdx edge_pos = 0;
            if ((degj < degi) || (degj==degi && j <= i))
                edge_pos = g->getEdgeBinary(j,i);
            else
                edge_pos = pos;


            for (EdgeIdx tri_pos = allTris->trioffsets[edge_pos]; tri_pos < allTris->trioffsets[edge_pos+1]; tri_pos++)
            {
                k = allTris->triangles[tri_pos];
                VertexIdx degk = g->offsets[k+1] - g->offsets[k];
                EdgeIdx new_pos = 0;
                if ((degk < degj) || (degk==degj && k <= j))
                    new_pos = g->getEdgeBinary(k,j);
                else
                    new_pos = g->getEdgeBinary(j,k);
                if (DEBUG)
                    printf("---Handling %lld %lld %lld\n",i,j,k);


                low = allTris->trioffsets[new_pos];
                high = allTris->trioffsets[new_pos+1]-1;


                while(low <= high)
                {
                     mid = (low+high)/2;
                     if (allTris->triangles[mid] == i)
                         break;
                     if (allTris->triangles[mid] > i)
                         high = mid-1;
                     else
                         low = mid+1;
                }
                if (low > high)
                {
                    printf("Error in binary search of allTris %lld %lld %lld\n",i,j,k);
                    printf("%lld %lld: ",i,j);
                    printf("%lld\n",allTris->trioffsets[edge_pos]);
                    for (EdgeIdx tri_pos2 = allTris->trioffsets[edge_pos]; tri_pos2 < allTris->trioffsets[edge_pos+1]; tri_pos2++)
                        printf("%lld ",allTris->triangles[tri_pos2]);
                    printf("\n");
                    printf("%lld %lld: ",j,k);
                    printf("%lld\n",allTris->trioffsets[new_pos]);
                    for (EdgeIdx tri_pos3 = allTris->trioffsets[new_pos]; tri_pos3 < allTris->trioffsets[new_pos+1]; tri_pos3++)
                        printf("%lld ",allTris->triangles[tri_pos3]);
                    printf("\n");
                    printf("%lld %lld: ",k,j);
                    EdgeIdx new_pos4 = g->getEdgeBinary(k,j);
                    printf("%lld\n",allTris->trioffsets[new_pos4]);
                    for (EdgeIdx tri_pos4 = allTris->trioffsets[new_pos4]; tri_pos4 < allTris->trioffsets[new_pos4+1]; tri_pos4++)
                        printf("%lld ",allTris->triangles[tri_pos4]);
                    printf("\n");
                    for (EdgeIdx ptr = g->offsets[j]; ptr < g->offsets[j+1]; ptr++)
                        printf("%lld ",g->nbors[ptr]);
                    printf("\n");
                    exit(EXIT_FAILURE);
                }
                
                for (EdgeIdx next_tri_pos = mid; next_tri_pos < allTris->trioffsets[new_pos+1]; next_tri_pos++)
                {
                    ell = allTris->triangles[next_tri_pos];
                    if (wedge_count[ell] > 2)
                    {
                        if (DEBUG)
                            printf("Consider %lld %lld %lld %lld: adding %lld\n",i,j,k,ell,wedge_count[ell]-2);

                        ret.chordalWedgeCol += wedge_count[ell]-2;
                    }
                    if (i < k && i < ell)
                    {
                        diamond_count[ell]++;
                    }
                }
            }
            
            for (EdgeIdx tri_pos = allTris->trioffsets[edge_pos]; tri_pos < allTris->trioffsets[edge_pos+1]; tri_pos++)
            {
                k = allTris->triangles[tri_pos];
                if (k <= i)
                    continue;
                VertexIdx degk = g->offsets[k+1] - g->offsets[k];
                VertexIdx new_pos = 0;
                if ((degk < degj) || (degk==degj && k <= j))
                    new_pos = g->getEdgeBinary(k,j);
                else
                    new_pos = g->getEdgeBinary(j,k);
                if (DEBUG)
                    printf("---Handling %lld %lld %lld\n",i,j,k);


                low = allTris->trioffsets[new_pos];
                high = allTris->trioffsets[new_pos+1]-1;

                while(low <= high)
                {
                     mid = (low+high)/2;
                     if (allTris->triangles[mid] == i)
                         break;
                     if (allTris->triangles[mid] > i)
                         high = mid-1;
                     else
                         low = mid+1;
                }
               
                for (EdgeIdx next_tri_pos = mid; next_tri_pos < allTris->trioffsets[new_pos+1]; next_tri_pos++)
                {
                    ell = allTris->triangles[next_tri_pos];
//                     printf("ell %lld, count %lld, total %lld\n",ell,diamond_count[ell],ret.wheel);
                    ret.wheel += diamond_count[ell]*(diamond_count[ell]-1)/2;
                    diamond_count[ell] = 0;
                }
            }
        }

        for (EdgeIdx pos = g->offsets[i]; pos < g->offsets[i+1]; ++pos) // loop over in-neighbors of i
        {
            j = g->nbors[pos]; // j is current in-neighbor
            for (EdgeIdx next = g->offsets[j]; next < g->offsets[j+1]; ++next) // loop over out-neighbors of j, note this gives an outout wedge
            {
                k = g->nbors[next];  // i <- j -> k is outout wedge centered at j
//                 if (k <= i) // k is higher in the order, so ignore wedge
//                     continue;
                if (DEBUG)
                    if (wedge_count[k] != 0)
                        printf("%lld %lld %lld: %lld\n",i,j,k,wedge_count[k]);
                ret.threeWedgeCol += (wedge_count[k])*(wedge_count[k]-1)*(wedge_count[k]-2)/6;
                wedge_count[k] = 0;
            }
        }
   }
   ret.chordalWedgeCol = ret.chordalWedgeCol/2;
   return ret;
}




#endif


