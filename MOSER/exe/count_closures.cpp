#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/Graph.h"
#include "Escape/GetAllCounts.h"


/* This code is executed as follows:

./count_closures <PATH FOR GRAPH>

It generates file out.txt in the following format. 
the first line is: n m
Every subsequent line has three numbers: i num_pairs num_closed

This means that there are num_pairs pairs of vertices with exactly i vertices in
common, and num_closed of them are closed (have an edge)
*/

using namespace Escape;

int main(int argc, char *argv[])
{
  Graph g;
  if (loadGraph(argv[1], g, 1, IOFormat::escape))
    exit(1);

  printf("Loaded graph\n");
  CGraph cg = makeCSR(g);  // generate CGraph
  cg.sortById();           // sort lists by Id
  printf("Converted to CSR\n");

  printf("Relabeling graph\n");
  CGraph cg_relabel = cg.renameByDegreeOrder();   // relabel by degree order
  cg_relabel.sortById();

  Count *common = new Count[cg_relabel.nVertices+1];  // initializing output arrays
  Count *closed = new Count[cg_relabel.nVertices+1];

  Count maximum = cClosure(&cg_relabel, common, closed);  // this gets all te closure information

  FILE* f = fopen("out.txt","w");
  if (!f)
  {
      printf("could not write to output to out.txt\n");
      return 0;
  }
  
  fprintf(f,"%lld %lld\n",cg_relabel.nVertices, cg_relabel.nEdges/2);  // print n m
  for(Count i = 1; i <= maximum; i++)   // print all the non-trivial closure information
  {
      if (common[i] != 0)
      fprintf(f,"%lld %lld %lld\n",i,common[i],closed[i]);
  }

  fclose(f);
}
