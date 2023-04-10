#ifndef ESCAPE_ALMOSTFIVECLIQUE_H_
#define ESCAPE_ALMOSTFIVECLIQUE_H_

#include <algorithm>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/Utils.h"

using namespace Escape;

// Requires gfull to be sorted by vertex ID
EdgeIdx almostFiveClique(CGraph *gsorted)
{
    EdgeIdx ret = 0;
    int DEBUG = 0;

    EdgeIdx *triangles = new VertexIdx[gsorted->nEdges+1];
    VertexIdx count = 0;

    for (VertexIdx i=0; i < gsorted->nVertices; i++)
    {
        VertexIdx degi = gsorted->offsets[i+1] - gsorted->offsets[i];
        for (EdgeIdx pos = gsorted->offsets[i]; pos < gsorted->offsets[i+1]; pos++)
        {
            VertexIdx j = gsorted->nbors[pos];
            VertexIdx degj = gsorted->offsets[j+1] - gsorted->offsets[j];
            
            if (degj < degi || (degj == degi && j < i))
                continue;
            count = 0;

            if (DEBUG)
                printf("Processing %lld %lld\n",i,j);

            for (EdgeIdx ptr = gsorted->offsets[i]; ptr < gsorted->offsets[i+1]; ptr++)
            {
                VertexIdx nbr = gsorted->nbors[ptr];
//                 VertexIdx degnbr = gsorted->offsets[nbr+1] - gsorted->offsets[nbr];
                if (DEBUG)
                    printf("Looking for %lld %lld\n",nbr,j);
                if (gsorted->getEdgeBinary(nbr,j) != -1)
                {
                    triangles[count] = nbr;
                    if (DEBUG)
                        printf("Found triangle %lld %lld %lld\n",i,j,nbr);
                    count++;
                }
            }

            std::sort(triangles,triangles+count);


            for (VertexIdx ptr = 0; ptr < count; ptr++)
            {
                VertexIdx nbr = triangles[ptr];
                VertexIdx degnbr = gsorted->offsets[nbr+1] - gsorted->offsets[nbr];

                if (degj < degnbr || (degj == degnbr && j < nbr))
                    continue;
                if (degnbr < degi || (degnbr == degi && nbr < i))
                    continue;

                if (DEBUG)
                    printf("Looking at %lld %lld %lld\n",i,j,nbr);

                VertexIdx fourclique = 0;
                
                if (count < degnbr)
                {
                    if (DEBUG)
                        printf("triangles < degnbr\n");
                    for (VertexIdx ptr2 = 0; ptr2 < count; ptr2++)
                    {
                        if (DEBUG)
                            printf("Looking for edge between %lld %lld\n",triangles[ptr2],nbr);
                        if (gsorted->getEdgeBinary(triangles[ptr2],nbr) != -1)
                        {
                            if (DEBUG)
                            {
                                printf("Four-clique %lld %lld %lld %lld\n",i,j,nbr,triangles[ptr2]);
                            }
                            fourclique++;
                        }
                    }
                }
                else
                {
                    if (DEBUG)
                        printf("triangles >= degnbr\n");
                    for (VertexIdx next = gsorted->offsets[nbr]; next < gsorted->offsets[nbr+1]; next++)
                    {
                        if (binarySearch(triangles,count,gsorted->nbors[next]) != -1)
                        {
                            fourclique++;
                            if (DEBUG)
                                printf("Four-clique %lld %lld %lld %lld\n",i,j,gsorted->nbors[next],triangles[binarySearch(triangles,count,gsorted->nbors[next])]);
                        }
                    }
                }
                if (DEBUG)
                    printf("Fourclique = %lld\n",fourclique);
                ret += fourclique*(fourclique-1)/2;
            }
        }
    }

    if (DEBUG)
        printf("Total found %lld\n",ret);
//     return ret - 10*fivecliques;
    return ret;
}

// Requires gfull to be sorted by vertex ID
EdgeIdx oldAlmostFiveClique(CGraph *gsorted)
{
    EdgeIdx ret = 0;
    int DEBUG = 0;

    EdgeIdx *triangles = new VertexIdx[gsorted->nEdges+1];
    VertexIdx count = 0;

    for (VertexIdx i=0; i < gsorted->nVertices; i++)
    {
        for (EdgeIdx pos = gsorted->offsets[i]; pos < gsorted->offsets[i+1]; pos++)
        {
            VertexIdx j = gsorted->nbors[pos];
            if (j < i)
                continue;
            count = 0;

            if (DEBUG)
                printf("Processing %lld %lld\n",i,j);

            EdgeIdx jpos = gsorted->offsets[j];
            for (EdgeIdx ipos = gsorted->offsets[i]; ipos < gsorted->offsets[i+1]; ipos++)
            {
                VertexIdx inbr = gsorted->nbors[ipos];
                while(jpos < gsorted->offsets[j+1] && (gsorted->nbors[jpos] < inbr))
                    jpos++;
                if(gsorted->nbors[jpos] == inbr)
                {
                    triangles[count] = inbr;
                    count++;
                }
            }

            std::sort(triangles,triangles+count);

            for (VertexIdx ptr = 0; ptr < count; ptr++)
            {
                VertexIdx k = triangles[ptr];
                if (i > k || k > j)
                    continue;
                VertexIdx fourclique = 0;

                for (VertexIdx ptr2 = 0; ptr2 < count; ptr2++)
                {
                    if(gsorted->isEdge(k,triangles[ptr2]) != -1)
                    {
                        if (DEBUG)
                            printf("4-clique: %lld %lld %lld %lld\n",i,j,k,triangles[ptr2]);

                        fourclique++;
                    }
                }
                if (DEBUG)
                    printf("Total for %lld %lld %lld is %lld\n",i,j,k,fourclique);
                ret += fourclique*(fourclique-1)/2;
            }
        }
    }

    if (DEBUG)
        printf("Total found %lld\n",ret);
//     return ret - 10*fivecliques;
    return ret;
}


EdgeIdx testingAlmostFiveClique(CGraph *gsorted)
{
    EdgeIdx ret = 0;
    int DEBUG = 0;

    EdgeIdx *triangles = new VertexIdx[gsorted->nEdges+1];
    VertexIdx count = 0;
    VertexIdx smaller, bigger;

    for (VertexIdx i=0; i < gsorted->nVertices; i++)
    {
        for (EdgeIdx pos = gsorted->offsets[i]; pos < gsorted->offsets[i+1]; pos++)
        {
            VertexIdx j = gsorted->nbors[pos];
            if (j < i)
                continue;
            count = 0;

            if (DEBUG)
                printf("Processing %lld %lld\n",i,j);

            VertexIdx degi = gsorted->offsets[i+1] - gsorted->offsets[i];
            VertexIdx degj = gsorted->offsets[j+1] - gsorted->offsets[j];

            if (degi <= degj)
            {
                smaller = i;
                bigger = j;
            }
            else
            {
                smaller = j;
                bigger = i;
            }

            for (EdgeIdx ptr = gsorted->offsets[smaller]; ptr < gsorted->offsets[smaller+1]; ptr++)
            {
                VertexIdx nbr = gsorted->nbors[ptr];
                if (nbr < i || nbr > j)
                    continue;
                if (gsorted->getEdgeBinary(nbr,bigger) != -1)
                {
                    triangles[count] = nbr;
                    count++;
                }
            }

            std::sort(triangles,triangles+count);

//             for (VertexIdx ptr = 0; ptr < count; ptr++)
//             {
//                 VertexIdx k = triangles[ptr];
//                 if (i > k || k > j)
//                     continue;
//                 VertexIdx fourclique = 0;
// 
//                 for (VertexIdx ptr2 = 0; ptr2 < count; ptr2++)
//                 {
//                     if(gsorted->isEdge(k,triangles[ptr2]) != -1)
//                     {
//                         if (DEBUG)
//                             printf("4-clique: %lld %lld %lld %lld\n",i,j,k,triangles[ptr2]);
// 
//                         fourclique++;
//                     }
//                 }
//                 if (DEBUG)
//                     printf("Total for %lld %lld %lld is %lld\n",i,j,k,fourclique);
//                 ret += fourclique*(fourclique-1)/2;
//             }
        }
    }

    if (DEBUG)
        printf("Total found %lld\n",ret);
//     return ret - 10*fivecliques;
    return ret;
}
#endif



