#ifndef ESCAPE_FIVE_FROM_TRIANGLES_H_
#define ESCAPE_FIVE_FROM_TRIANGLES_H_


#include "Escape/Utils.h"
#include "Escape/FourVertex.h"
#include "Escape/TriangleProgram.h"

namespace Escape
{

struct FiveFromTriangles
{
    EdgeIdx forktailedtris;
    EdgeIdx longtailedtris;
    EdgeIdx doubletailedtris;
};

FiveFromTriangles fiveFromTriCounter(CGraph *g, CGraph *gout, TriangleInfo *tri_info, SomeFourPatterns fourcounts)
{
    FiveFromTriangles ret;

    ret.forktailedtris = 0;
    ret.longtailedtris = 0;
    ret.doubletailedtris = 0;

    for (VertexIdx i = 0; i < g->nVertices; i++)
    {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];

        //forkedtailedtris = \sum_i tri[i]*(degi-2 \choose 2)
        ret.forktailedtris += tri_info->perVertex[i]*((degi-2)*(degi-3))/2;

        for (EdgeIdx posj = g->offsets[i]; posj < g->offsets[i+1]; posj++)
        {
            VertexIdx j = g->nbors[posj];
            VertexIdx degj = g->offsets[j+1] - g->offsets[j];

            //longtailedtris = \sum_{(i,j) in E} tri[j]*(degj-1) - 2*#chordal-cycles
//             printf("%ld %ld\n",tri_info->perVertex[i],degj-1);
            ret.longtailedtris += tri_info->perVertex[i]*(degj-1);
        }
    }
    ret.longtailedtris -= 4*fourcounts.chordalcycles;
    ret.longtailedtris -= 6*tri_info->total;
    ret.longtailedtris -= 2*fourcounts.tailedtris;


    for (VertexIdx i = 0; i < g->nVertices; i++)
    {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
        for (EdgeIdx posj = gout->offsets[i]; posj < gout->offsets[i+1]; posj++)
        {
            VertexIdx j = gout->nbors[posj];
            VertexIdx degj = g->offsets[j+1] - g->offsets[j];
//             printf("%ld %ld %ld\n",tri_info->perEdge[posj],degi-2,degj-2);
            ret.doubletailedtris += tri_info->perEdge[posj]*(degi-2)*(degj-2);
        }
    }

    ret.doubletailedtris -= 2*fourcounts.chordalcycles;
    return ret;
}

//Pattern descriptions: numbers are as per Sesh's notes
//
//pattern 9: "the hourglass": two triangles that intersect at one vertex
//and no edges in common.
//             o---o
//              \ /
//               o
//              / \                                                  c
//             o---o
//
//If v is the common vertex, the number of triangles intersecting at v is
//nTr(v) choose 2.  From this we have to subtract the number of cases where
//a pair of triangles share an edge.  Note that it is not possible for two
//distinct triangles to share more than 1 edge.
//
//We can optionally assign the count of this pattern to the central vertex.
//
//
//Template Params:
//  doPerVertex: set to true if you want to record per-vertex counts.  If false
//   the counts argument is ignored.
//
//Arguments:
//  g: a CGraph. (Could equally well work with a Graph)
//  tInfo: triangle info.  If not provided, will be computed here.
//  counts: an optional pointer to caller-allocated storage for a per-vertex
//         count (vertex corresponds to center in the diagram above). Template
//         parameter doPerVertex should be set.
//
//Return:
//  total count of pattern in the graph.
template <bool doPerVertex>
Count count5_Hourglass(const CGraph* g
  , const TriangleInfo* tInfo = 0
  , Count *counts = 0)
{
  if (!tInfo)
  {
    //compute it here, todo
    printf("unimplemented code path\n");
    exit(1);
  }

  //set counts to zero
  if (doPerVertex)
    std::fill(counts, counts + g->nVertices, 0);

  Count total = 0;
  for (VertexIdx u = 0; u < g->nVertices; ++u)
  {
    Count c = choose<2>(tInfo->perVertex[u]);
    total += c;
    if (doPerVertex)
      counts[u] += c;

    for (EdgeIdx e = g->offsets[u]; e < g->offsets[u + 1]; ++e)
    {
      auto v = g->nbors[e];
      Count c = choose<2>(tInfo->perEdge[e]);
      if (doPerVertex)
      {
        counts[u] -= c;
        counts[v] -= c;
      }
      total -= 2 * c;
    }
  }

  return total;
}


//Pattern 11, "the stingray"
//
//      o <------ u
//     /|\                                                            c
//    /.|.\                                                           c
//   o  |<-o----- e
//    \ | /
//     \|/
//      o <------ v
//      |
//      |
//      o
//
//
//Let e = (u, v) be the edge that forms the spine.  Then the number of such
//patterns is:
//    nTr(e) choose 2 * (max(d(v) - 3, 0) + max(d(u) - 3, 0))
//where d(u) is the total degree of u. This count can be associated with edge e.
//
//Arguments:
//  gOut: CSR graph
//  gIn:  CSC graph (only used for getting total degree).  We should precompute
//        and pass around total degree in the graph structure perhaps.
//  tInfo: triangle info.  If not provided, will be computed here.
//  perEdgeOut: an optional pointer to caller-allocated storage for a per-edge
//         count (edge e in the diagram above).
//
//Return:
//  total count of pattern in the graph.

template <bool doPerEdge>
Count count5_Stingray(const CGraph *gOut
  , const CGraph *gIn
  , const TriangleInfo *tInfo = 0
  , Count* counts = 0)
{
  if (!tInfo)
  {
    //compue here todo
    printf("unimplemented, please fix\n");
    exit(1);
  }

  if (doPerEdge)
    std::fill(counts, counts + gOut->nEdges, 0);

  Count total = 0;
  for (VertexIdx u = 0; u < gOut->nVertices; ++u)
  {
    auto du = gOut->degree(u) + gIn->degree(u);
    auto cu = czsub<Count>(du, 3);
    for (EdgeIdx e = gOut->offsets[u]; e < gOut->offsets[u + 1]; ++e)
    {
      auto v = gOut->nbors[e];
      auto dv = gOut->degree(v) + gIn->degree(v);
      auto cv = czsub<Count>(dv, 3);
      Count c = choose<2>(tInfo->perEdge[e]) * (cu + cv);
      total += c;
      if (doPerEdge)
        counts[e] = c;
    }
  }

  return total;
}



//Pattern 14: "the stellate trident" for want of a better name
//Dente == tooth, == triangle in analogy with saw-tooth functions
//Tri-dente because there are three
//Stellate because they are arranged in a star.
//Unfortunately does not look anything like an ordinary trident. :)
//
//  Anyone know the name of a children's toy that one holds in the hand and
//spins around? There are balls on the ends of three triangles which share
//the rotational axis and the momentum is transferred from ball to ball with
//a clacking sound as the child spins the toy with their hand.  Same topology
//anyhow.
//
//               .o-.
//              / |\ \                                                      c
//             o  | o o
//              \ |/ /
//               `o-'
//
//The count for this is nTr(e) choose 3 for the common edge e.
//
//Template Arguments
//  doPerEdge: whether a count per edge is desired (assigned to the common edge)
//
//Arguments:
//  g: CSR or CSC graph (COO would have been fine too)
//  tInfo: triangle info.  If not provided, will be computed here.
//  counts: optional caller-allocated storage for per-edge counts (if doPerEdge set)
//
//Return:
//  Count
template <bool doPerEdge>
Count count5_StellateTrident(const CGraph *g
  , const TriangleInfo *tInfo = 0
  , Count* counts = 0)
{
  if (!tInfo)
  {
    //todo
    printf("unimplemented\n");
    exit(1);
  }

  if (doPerEdge)
    std::fill(counts, counts + g->nEdges, 0);

  Count total = 0;
  for (VertexIdx u = 0; u < g->nVertices; ++u)
  {
    for (EdgeIdx e = g->offsets[u]; e < g->offsets[u + 1]; ++e)
    {
      Count c = choose<3>(tInfo->perEdge[e]);
      total += c;
      if (doPerEdge)
        counts[e] = c;
    }
  }
  return total;
}


//Pattern 16: "the triangle strip" (from opengl)
//
//     o---o
//    / \ / \                                c
//   o---o---o
//
//Method:
//
//Let tr(e) be the number of triangles that edge e participates in.
//For each triangle A, B, C in the graph:
//
//     A---C
//    / \ / \                                c
//   X---B---Y
//
//The number of strips as shown above is number of X's times the number of Y's.
//(tr(AB) - 1) * (tr(BC) - 1) is the number of such strips including the cases
//where X and Y are identical.  But if X and Y are identical, then the pattern
//is a 4-clique and not a triangle strip.  So the number of strips with distinct
//X and Y are:
//   (tr(AB) - 1) * (tr(BC) - 1) - k4(ABC)
//where k4 is the number of 4-cliques that triangle ABC participates in.
//
//Similarly, one can get the count for the strips of the form:
//       Y
//      / \                                                  c
//     A---C
//    / \ /
//   X---B
//   (tr(AB) - 1) * (tr(AC) - 1) - k4(ABC)
//
//and
//
//       Y
//      / \                                                  c
//     A---C
//      \ / \                                                c
//       B---X
//
//   (tr(AC) - 1) * (tr(BC) - 1) - k4(ABC)
//
//Now, if we sum these terms over all triangles t, we get the term sum_t 3 * k4(t).
//But if k4(G) is the total number of 4-cliques in G, then sum_t k4(t) = 4 * k4(G)
//because each clique contains 4 triangles.
//
//So the final formula is:
//
//count = sum_{triangles ABC}
//                 (tr(AB) - 1) * (tr(BC) - 1)
//               + (tr(BC) - 1) * (tr(AC) - 1)
//               + (tr(AC) - 1) * (tr(AB) - 1)
//       -12 * k4(G)
//
//
//Template arguments:
//  isDOG: if the input graph is already degree-ordered.
//  Note: presently only works if isDOG is true.
//
//Arguments:
//  gOut: graph in CSR representation
//  k4:   the number of 4-cliques in the graph.
//  tInfo: triangle information, optional.  Presently this is ignored if
//         the input graph is not already a DOG, since we don't have the
//         mapping between edges in the original graph and the DOG.
template <bool isDOG = false>
Count count5_TriangleStrip(const CGraph *gOut
  , Count k4
  , const TriangleInfo *tInfo_ = 0)
{
  static_assert(isDOG, "sorry isDOG must be true for now");
  CDAG dog;
  const CGraph *g;

  TriangleInfo tInfo;

  if (!isDOG)
  {
    //just make the DOG here, we can do better, but it's a complicated
    //mish-mash of template combinations, for later.
    dog = degreeOrdered((CGraph*) gOut);
    tInfo = newTriangleInfo(gOut);
    tInfo.total = countTriangles(&dog.outlist, true, true, &dog.inlist
      , tInfo.perVertex, tInfo.perEdge);
    g = &dog.outlist;
  }
  else
  {
    g = gOut;
    tInfo = *tInfo_;
  }

  //This is executed for each triangle.
  struct Functor
  {
    Count* perEdge;
    Count total;

    void operator ()(VertexIdx u, VertexIdx v, VertexIdx w, EdgeIdx vu, EdgeIdx vw, EdgeIdx uw)
    {
      //this may be slow - may want to compute these per edge outside the triangle loop instead.
      Count c0 = czsub<Count>(perEdge[vu], 1);
      Count c1 = czsub<Count>(perEdge[vw], 1);
      Count c2 = czsub<Count>(perEdge[uw], 1);
      total += c0 * c1 + c0 * c2 + c1 * c2;
    }
  };

  Functor f = {tInfo.perEdge, 0};
  triangleProgram<true, true, false>(g, f);

  //this just matches the conditions under which we did newTriangleInfo.  If you
  //change the handling of isDOG above, make this match!
  if (tInfo_ == 0 || !isDOG)
  {
    delTriangleInfo(tInfo);
    delCGraph(dog.outlist);
    delCGraph(dog.inlist);
  }

  return f.total - 12 * k4;
}


//Pattern 10: "the Cobra", not to be confused with the Stingray.  FYI some
//Cobras do have the horiz line between the fake 'eyes'.
//
//        o
//       / \                                                             c
//      o---o <--- e
//       \ /
//        o <---- u
//        |
//        o
//
//Method:
//  For each triangle that an edge e participates in, the number of pattern 10's
//with that triangle in the lower half of the figure above is:
//  (nTr(e) - 1) * (d(u) - 2)
//where d(u) is the total degree of u and d(u)-2 is clamped to zero.
//But this includes the case where the tail vertex is identical to the head
//vertex causing the resulting pattern to be a 4-clique and not a cobra. If edge
//e is part of a 4-clique, then it gets 2 fake cobras from that 4-clique.
//So the total number of these patterns for that edge e is:
//  [sum_tr(e) (nTr(e) - 1) * (d(u) - 2)] - 2 * k4(e)
//  = [(nTr(e) - 1) sum_tr(e) (d(u) - 2)] - 2 k4(e)
//where k4 is the number of 4-cliques in which e participates and d(u) is the
//total (in+out) degree of u.

//Currently we don't have k4 per edge, so we only return the sum of the
//above expression over all edges in the graph, but this would be easy to
//remedy if we track k4(e).
//
//Note that sum_e k4(e) = 6 * k4(G) since each 4-clique has 6 edges.
//
//Template parameters:
//  isDOG: set to true if the graph is already degree-ordered.
//
//Arguments:
//  gOut: CSR (out-edges) of the graph
//  gIn:  CSC (in-edges) of the graph, needed for total degree of u.
//  k4:  total number of 4-cliques in the graph
//
//Return:
//  Total count of pattern 10 in the graph

template <bool isDOG = false>
Count count5_Cobra(const CGraph* gOut, const CGraph* gIn, Count k4)
{
  static_assert(isDOG, "sorry isDOG must be true for now");

  Count *d2; //this is max(total_degree(u) - 2, 0) for each vertex u.

  d2 = new Count[gOut->nVertices];
  for (VertexIdx u = 0; u < gOut->nVertices; ++u)
    d2[u] = czsub<Count>(gOut->degree(u) + gIn->degree(u), 2);

  struct Functor
  {
    const Count *d2; //precalculated d2 values above.

    Count *tr; //the number of triangles per edge
    Count *du; //sum_u max(0, d(u) - 2) per edge


    Functor(EdgeIdx nEdges, Count* d2) : d2(d2)
    {
      tr = new Count[nEdges]();
      du = new Count[nEdges]();
    }

    ~Functor()
    {
      delete[] tr;
      delete[] du;
    }

    void operator ()(VertexIdx u, VertexIdx v, VertexIdx w, EdgeIdx vu, EdgeIdx vw, EdgeIdx uw)
    {
      tr[vu] += 1;
      du[vu] += d2[w];

      tr[vw] += 1;
      du[vw] += d2[u];

      tr[uw] += 1;
      du[uw] += d2[v];
    }
  };

  Functor f(gOut->nEdges, d2);
  triangleProgram<isDOG, true, false>(gOut, f, gIn);

  Count total = 0;
  for (EdgeIdx e = 0; e < gOut->nEdges; ++e)
    total += f.du[e] * czsub<Count>(f.tr[e], 1);
  total -= 12 * k4;

  delete[] d2;
  return total;
}










}

#endif
