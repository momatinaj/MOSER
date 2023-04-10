#ifndef ESCAPE_TRIANGLE_PROGRAM_H_
#define ESCAPE_TRIANGLE_PROGRAM_H_


#include "Graph.h"
#include "JointSort.h"
#include <algorithm>


namespace Escape
{

//Note: this is a slight generalization of wedgeEnumerator in Triadic.h
//
//Run some code for each instance of a triangle in a graph.  A triangle is
//defined as any set of three vertices such that all pairs are connected.  Note
//that this definition does not care about the directed-ness of the graph or
//the edges.
//
//The provided function is executed for each distinct triangle exactly once
//(but see parameters below).
//
//DAG construction: to avoid duplicate calls for different permutations of the
//same triangle, as well as to reduce the number of edge checks, we use an
//implicit DAG that is constructed from the input graph: for each edge of the
//input graph, (re)direct the edge so that it points from a vertex of lower
//total degree to a vertex of higher total degree. If tied in degree, the edge
//points from lower vertex id to higher vertex id.
//
//Note: instead of producing the entire degree-ordered graph, we keep a
//temporary patch of it in the neighborhood that we are currently exploring.
//
//Template parameters:
//
//isDOG: set to true if the input graph is already degree-ordered (because you
//       created the degree-ordered DAG explicitly using Digraph.h/degreeOrdered()
//
//doDegreeOrdering: if isDOG then this parameter is ignored.  Otherwise it can be
//                  used so enable/disable inline degree-ordering while enumerating
//                  triangles.  Note: if you disable degree-ordering and your
//                  input graph is directed with cycles, those cycles will _not_
//                  be enumerated as part of this procedure.  If the graph is
//                  undirected, then triangles will be enumerated an
//                  unspecified number of times (thrice in the current implementation).
//                  Bottom line: isDOG=false and doDegreeOrdering=false is a
//                  very unlikely corner case.
//
//isUndirected:     Set to true if the input graph is an undirected graph (with
//                  explicit reverse edges)
//
//f: a functor that has a method with the following signature:
//
//    operator ()(VertexIdx u, VertexIdx v, VertexIdx w, EdgeIdx vu, EdgeIdx vw, EdgeIdx uw)
//
//The edge index vu (and similarly vw and uw) refers to either an edge v->u
//or an edge u->v.  If the input graph is a DAG and you did not request degree
//ordering, then vu is precisely v -> u (and likewise for vw, uw).  If the
//input is an undirected graph, then you only get the edge index of one of each
//pair of directed dges for each side of the triangle - you can look up the
//others with g->getEdge.
//
//Todo: we need a better way to store and process undirected graphs.
//
//Inputs:
//  gOutList: CGraph with out-edges (CSR layout)
//  gInList:  CGraph with in-edges (CSC layout).  This is only needed if you
//    use doDegreeOrdering = true, or if the input graph is directed, otherwise
//    it can be null (or equal to gOutList)
//  
//
template <bool isDOG = false
  , bool doDegreeOrdering_ = true
  , bool isUndirected = false
  , class F>
void triangleProgram(const CGraph *gOutList, F& functor, const CGraph *gInList = 0)
{
  EdgeIdx *totalDegrees; //used to cache vertex degrees, only used if doDegreeOrdering
  VertexIdx *tmpNbors;   //temporary filtered and re-sorted list of neighbors for a vertex.
  EdgeIdx   *tmpOrigEdges; //track original edge indices in the reordered edge list.
  EdgeIdx    tmpnEdges;  //number of edges filtered by degree criterion

  //this just says that if isDOG, we will not do degree ordering regardless of
  //what was specified.
  constexpr bool doDegreeOrdering = !isDOG && doDegreeOrdering_;
  
  if (doDegreeOrdering)
  {
    //if CSC not provided, assume undirected graph with explicit reverse edges
    if (gInList == 0)
    {
      if (!isUndirected)
      {
        printf("you gave a directed graph, requested degree-ordering, but no inList\n");
        abort();
      }
      gInList = gOutList;
    }

    //store the in + out degree for each vertex and track the max out-degree at
    //the same time.
    EdgeIdx maxOutDegree = 0;
    totalDegrees = new EdgeIdx[gOutList->nVertices];
    for (VertexIdx v = 0; v < gOutList->nVertices; ++v)
    {
      auto td = gOutList->degree(v) + gInList->degree(v);
      maxOutDegree = std::max(maxOutDegree, td);
      totalDegrees[v] = td;
    }

    //make space for temporary nborlist
    tmpNbors = new VertexIdx[maxOutDegree];
    tmpOrigEdges = new EdgeIdx[maxOutDegree];
  }
  
  for (VertexIdx v = 0; v < gOutList->nVertices; ++v)
  {
    if (doDegreeOrdering)
    {

      //v1 is less than v2 if in the degree-ordered graph if
      //-  totalDegree(v1) < totalDegree(v2)
      //-  or totalDegrees are equal and index v1 < index v2
      auto dogLess = [&](VertexIdx v1, VertexIdx v2)
      {
        return totalDegrees[v1] < totalDegrees[v2]
          || (totalDegrees[v1] == totalDegrees[v2] && v1 < v2);
      };
      
      //populate the neighbors of v
      tmpnEdges = 0;
      for (auto e = gOutList->offsets[v]; e < gOutList->offsets[v + 1]; ++e)
      {
        auto nb = gOutList->nbors[e];
        if (dogLess(v, nb))
        {
          tmpNbors[tmpnEdges] = nb;
          tmpOrigEdges[tmpnEdges] = e;
          ++tmpnEdges;
        }
      }

      //sort the neighbors of v by same criterion, and keep track of the original
      //edge indices.
      auto it = JSIterator<VertexIdx, EdgeIdx> {tmpNbors, tmpOrigEdges};
      std::sort(it, it + tmpnEdges, dogLess);
    }
    else
    {
      tmpNbors   = gOutList->nbors + gOutList->offsets[v];
      tmpnEdges  = gOutList->degree(v);
    }

    //Now loop over the (filtered) out-edges of v
    for (EdgeIdx eu = 0; eu < tmpnEdges; ++eu)
    {
      VertexIdx u = tmpNbors[eu];

      for (EdgeIdx ew = eu + 1; ew < tmpnEdges; ++ew)
      {
        VertexIdx w = tmpNbors[ew];

        //check if u -> w is an edge.  For u -> w to be an edge in the
        //degree-ordered graph, either u -> w or u <- w must exist in
        //the original graph.
        //If the input graph is already degree-ordered, the check w -> u is unnecessary
        //if the input graph is undirected, the check w -> u is unnecessary
        auto uw = gOutList->getEdge(u, w);
        EdgeIdx vu, vw;

        if (doDegreeOrdering) //get edge index in original graph
        {
          vu = tmpOrigEdges[eu]; 
          vw = tmpOrigEdges[ew];
        }
        else //edge mapping is almost identity
        {
          vu = eu + gOutList->offsets[v];
          vw = ew + gOutList->offsets[v];
        }
        
        if (isDOG || isUndirected)
        {
          if (uw != invalidEdge)
          {
            functor(u, v, w, vu, vw, uw);
          }
        }
        else
        {
          auto wu = gOutList->getEdge(w, u);
          //flip orientation in functor call.
          if (wu != invalidEdge)
            functor(w, v, u, vw, vu, wu);
        }
      }
    }
  }

  if (doDegreeOrdering)
  {
    delete[] tmpNbors;
    delete[] tmpOrigEdges;
    delete[] totalDegrees;
  }
}

//get total, per-vertex and per-edge triangle counts.
//
//Arguments:
//  gOut: outbound edges (CSR)
//  degreeOrdered: whether graph is degree-ordered
//  directed: whether graph is directed (implied true if degreeOrdered)
//  gIn:  inbound edges (CSC), optional. If undirected, this argument is ignored.
//  perVertex: optional pointer to caller-allocated data for a per-vertex count
//  perEdge: optional pointer to caller-allocated data for a per-edge count
//
//Return:
//  total number of triangles in the graph, triangles being defined as in
//  TriangleProgram
int64_t countTriangles(const CGraph* gOut
  , bool degreeOrdered = false
  , bool directed = true
  , const CGraph *gIn=0
  , int64_t *perVertex = 0
  , int64_t *perEdge = 0);



}; //namespace

#endif
