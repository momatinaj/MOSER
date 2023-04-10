/*   ////////////////////////////////////////
The code outputs the degree distributions of the degeneracy and degree ordered DAGS.
It is basically a wrapper, that calls the appropriate functions in Digraph.h
USAGE:
        ./dagdegdists<INPUT FILE> <OUTPUT FILE>

   <INPUT FILE>: This is file with graph in Escape format.
   <OUTPUT FILE>: File where output is given.

   The output file has two parts. The first part has the heading "Degree ordered".
   The following lines give the out and in-degree distributions of the degree ordered DAG.
   Each line of the output file will have 
        <degree> <number of vertices of this outdegree> <number of vertices of this indegree>
   with a line for every degree for which one of these counts is non-zero.

   The second part has the heading "Degeneracy ordered", and has the same information
   for the degeneracy ordered DAG.

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

  FILE* f = fopen(argv[2],"w"); // write to output file
  if (!f) // some problem in opening file
  {
      printf("Could not write to output to %s\n",argv[2]); // print error message and abort
      return 0;
  }

//   FILE* debug = fopen("test.txt","w"); // write some debug info to test file DEBUG ONLY
//   if (!debug) // some problem in opening file DEBUG ONLY
//   {
//       printf("Could not write to output to test.txt\n"); // print error message and abort
//       return 0;
//   }

  printf("Creating degree ordered DAG\n");
  CDAG dag = degreeOrdered(&cg); // construct the degree ordered DAG

  VertexIdx *outdegdistarray = new VertexIdx[cg.nVertices+1]; // allocate array of longs, with 
  VertexIdx *indegdistarray = new VertexIdx[cg.nVertices+1]; // allocate array of longs, with 

  VertexIdx maxoutdeg = degDist(&(dag.outlist), outdegdistarray); // get out-degree distribution 
  VertexIdx maxindeg = degDist(&(dag.inlist), indegdistarray); // get in-degree distribution 
  
  fprintf(f,"Degree ordered\n");
  for (VertexIdx i=0; i <= std::max(maxoutdeg, maxindeg); i++)
  {
      fprintf(f,"%lld %lld %lld\n",i,outdegdistarray[i],indegdistarray[i]);
  }

  printf("Creating degeneracy ordered DAG\n");
  dag = degenOrdered(&cg); // construct the degree ordered DAG

  maxoutdeg = degDist(&(dag.outlist), outdegdistarray); // get out-degree distribution 
  maxindeg = degDist(&(dag.inlist), indegdistarray); // get in-degree distribution 
 
  fprintf(f,"Degeneracy ordered\n");
  for (VertexIdx i=0; i <= std::max(maxoutdeg, maxindeg); i++)
  {
      fprintf(f,"%lld %lld %lld\n",i,outdegdistarray[i],indegdistarray[i]);
  }
//   (dag.outlist).print(debug); //print outlist of degeneracy dag into output file DEBUG ONLY
 
  fclose(f);
}
