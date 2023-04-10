#ifndef ESCAPE_EDGEHASH_H_
#define ESCAPE_EDGEHASH_H_

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"

#include <unordered_set>


namespace std
{
  template <> struct hash<std::pair<Escape::VertexIdx, Escape::VertexIdx>>
  {
    size_t operator () (const std::pair<Escape::VertexIdx, Escape::VertexIdx> &pair) const
    {
      //This is a bad hash combiner, just as an example.
      return std::hash<Escape::VertexIdx>()(pair.first)
        ^ std::hash<Escape::VertexIdx>()(pair.second);
    }
  };
}


namespace Escape
{

//This uses a hash to lookup whether an edge exists or not.
//
//Depending on how much precomputation can be done, using a minimal
//perfect hash, we can answer the question of whether an edge exists
//or not in O(1) time, using no additional memory beyond the original
//edge list.
//
//Currently just using the stdlib hash functions and is restricted
//to less than 2bn vertices.  Should be possible to get a good hashing
class EdgeHash
{
  using Edge = std::pair<VertexIdx, VertexIdx>;
  std::unordered_set<Edge> edgeSet;
  
  public:
    EdgeHash(const Graph& graph)
    {
      for (EdgeIdx i = 0; i < graph.nEdges; ++i)
        edgeSet.insert(std::make_pair(graph.srcs[i], graph.dsts[i]));
    }

    EdgeHash(const CGraph& graph)
    {
      //todo
    }
    
    ~EdgeHash() {} //nothing

    bool isEdge(VertexIdx a, VertexIdx b) const
    {
      return edgeSet.find(std::make_pair(a, b)) != edgeSet.end();
    }
};


}



#endif

