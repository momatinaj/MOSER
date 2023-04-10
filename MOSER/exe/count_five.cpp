#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/FourVertex.h"
#include "Escape/Utils.h"
#include "Escape/GetAllCounts.h"

using namespace Escape;

int main(int argc, char *argv[])
{
  Graph g;
  printf("Loading graph\n");
  if (loadGraph(argv[1], g, 1, IOFormat::escape))
    exit(1);

  printf("Converting to CSR\n");
  CGraph cg = makeCSR(g);

  printf("Relabeling graph\n");
  CGraph cg_relabel = cg.renameByDegreeOrder();
  cg_relabel.sortById();
  printf("Creating DAG\n");
  CDAG dag = degreeOrdered(&cg_relabel);

  (dag.outlist).sortById();
  (dag.inlist).sortById();

  double nonInd_three[4], nonInd_four[11], nonInd_five[34];

  printf("Counting 3-vertex\n");
  getAllThree(&cg_relabel, &dag, nonInd_three);
  printf("Counting 4-vertex\n");
  getAllFour(&cg_relabel, &dag, nonInd_four);
  printf("Counting 5-vertex\n");
  getAllFive(&cg_relabel, &dag, nonInd_four, nonInd_five);

  FILE* f = fopen("out.txt","w");
  if (!f)
  {
      printf("could not write to output to out.txt\n");
      return 0;
  }
  fprintf(f,"%lld\n",cg_relabel.nVertices);
  fprintf(f,"%lld\n",cg_relabel.nEdges/2);
  for(int i = 0; i < 4; i++)
      fprintf(f,"%f\n",nonInd_three[i]);
  for(int i = 0; i < 11; i++)
      fprintf(f,"%f\n",nonInd_four[i]);
  for(int i = 0; i < 34; i++)
      fprintf(f,"%f\n",nonInd_five[i]);

  fclose(f);
}

