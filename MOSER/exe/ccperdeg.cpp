/*   ////////////////////////////////////////
The code outputs the clustering coefficient per degree. It is basically a wrapper,
that calls the function ccPerDeg in Triadic.h
USAGE:
        ./ccperdeg <INPUT FILE> <OUTPUT FILE>

   <INPUT FILE>: This is file with graph in Escape format.
   <OUTPUT FILE>: File where output is given.

   Each line of the output file will have 
        <degree> <average cc for degree> <number of vertices of degree>
   with a line for every degree for which the count is non-zero.

C. Seshadhri, May 2019
////////////////////////////////////////////////// */

#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/Graph.h"
#include "Escape/GetAllCounts.h"

using namespace Escape;

int main(int argc, char *argv[])
{
  Graph g; //variable that stores input graph
  if (loadGraph(argv[1], g, 1, IOFormat::escape)) //load graph from input file
    exit(1);

  printf("Loaded graph\n");
  CGraph cg = makeCSR(g); // convert graph into CSR
  cg.sortById(); // sort all the vertices by id, within each neighbor list
  printf("Converted to CSR\n");

  printf("Relabeling graph\n");
  CGraph cg_relabel = cg.renameByDegreeOrder();  // relabel graph, so that vertex id is actually the rank in degree list. Thus, 0 is min degree vertex, 1 is vertex with next degree, etc.
  cg_relabel.sortById();
  printf("Creating DAG\n");
  CDAG dag = degreeOrdered(&cg_relabel); // construct the degree ordered DAG

  (dag.outlist).sortById(); // sort both the outlists and inlists by ID (which is now rank in degree list), to simplify other parts of code
  (dag.inlist).sortById();

  float *ccdegarray; // array of floats for clustering coefficients per degree
  ccdegarray = new float[cg.nVertices+1]; // allocate array of floats, with length being number of vertices (trivial bound on the maximum degree)
  VertexIdx *degdistarray; // array for degree distribution
  degdistarray = new VertexIdx[cg.nVertices+1]; // allocate array of longs, with 

  VertexIdx maxdeg = ccPerDeg(&cg_relabel, &(dag.outlist), ccdegarray); // call ccPerDeg, which populates ccdegarray
  VertexIdx maxdeg2 = degDist(&cg_relabel, degdistarray); // get degree distribution
  
  if (maxdeg != maxdeg2) // some error, since each function reporting different maximum degree
  {
      printf("Error: ccPerDeg and degDist reporting different maximum degrees, %lld and %lld, respectively\n",maxdeg,maxdeg2); //print error message and abort
      return 0;
  }

  VertexIdx i; //loop variable

  FILE* f = fopen(argv[2],"w"); // write to output file
  if (!f) // some problem in opening file
  {
      printf("Could not write to output to %s\n",argv[2]); // print error message and abort
      return 0;
  }

  for(i=0; i<=maxdeg; ++i) // loop over all degrees
      if (degdistarray[i] != 0) // if there are vertices of degree i
          fprintf(f,"%lld %.4f %lld\n",i,ccdegarray[i],degdistarray[i]); // print out relevant info into file
  
  fclose(f);
}
