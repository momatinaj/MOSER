#ifndef ESCAPE_CONVERSION_H_
#define ESCAPE_CONVERSION_H_

#include <algorithm>

#include "Escape/FourVertex.h"
#include "Escape/Utils.h"

using namespace Escape;


// This matrix converts induced 3-vertex pattern counts to non-induced three pattern counts

const long ThreeIndToNonMatrix[4][4] = { 
    {1, 1, 1, 1},
    {0, 1, 2, 3},
    {0, 0, 1, 3},
    {0, 0, 0, 1}
};

// This matrix converts non-induced three pattern counts to induced three pattern counts

const long ThreeNonToIndMatrix[4][4] = { 
    {1, -1, -1, -1},
    {0, 1, -2, 3},
    {0, 0, 1, -3},
    {0, 0, 0, 1}
};

// This matrix converts induced 4-vertex pattern counts to non-induced three pattern counts

const long FourIndToNonMatrix[11][11] = { 
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    {0, 1, 2, 2, 3, 3, 3, 4, 4, 5, 6},
    {0, 0, 1, 0, 0, 0, 1, 1, 2, 2, 3},
    {0, 0, 0, 1, 3, 3, 2, 5, 4, 8, 12},
    {0, 0, 0, 0, 1, 0, 0, 1, 0, 2, 4},
    {0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 4},
    {0, 0, 0, 0, 0, 0, 1, 2, 4, 6, 12},
    {0, 0, 0, 0, 0, 0, 0, 1, 0, 4, 12},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 3},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 6},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}
};

// This matrix converts non-induced 4-vertex pattern counts to induced three pattern counts

const long FourNonToIndMatrix[11][11] = { 
    {  1,  -1,   1,   1,  -1,  -1,  -1,   1,   1,  -1,   1},
    {  0,   1,  -2,  -2,   3,   3,   3,  -4,  -4,   5,  -6},
    {  0,   0,   1,   0,   0,   0,  -1,   1,   2,  -2,   3},
    {  0,   0,   0,   1,  -3,  -3,  -2,   5,   4,  -8,  12},
    {  0,   0,   0,   0,   1,   0,   0,  -1,   0,   2,  -4},
    {  0,   0,   0,   0,   0,   1,   0,  -1,   0,   2,  -4},
    {  0,   0,   0,   0,   0,   0,   1,  -2,  -4,   6, -12},
    {  0,   0,   0,   0,   0,   0,   0,   1,   0,  -4,  12},
    {  0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   3},
    {  0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -6},
    {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1}
};





// This converts non-induced 3-vertex pattern counts to induced counts, by
// a simple linear transformation. The matrix is given above as ThreeNonToIndMatrix.
//
// Input: array of length 4 with non-induced counts
// Output: array of length 4 with induced counts

void getThreeInduced(EdgeIdx (&nonind)[4], EdgeIdx (&ind)[4])
{
    // Just multiplying matrix ThreeNonToIndMatrix with vector nonind
    for (int i=0; i<4; i++)
    {
        ind[i] = 0;
        for (int j=0; j<4; j++)
            ind[i] += ThreeNonToIndMatrix[i][j]*nonind[j];
    }    
}


// The converts non-induced counts to induced counts, by
// a simple linear transformation of the counts. The exact
// equations are given in Equation (1) of 
// "Path Sampling: A Fast and Provable Method for Estimating 4-Vertex Subgraph Counts" (http://arxiv.org/pdf/1411.4942v1.pdf)
//
// Input: NoninducedFourCounts struct
// Output: InducedFourCounts struct

InducedFourCounts convertNonToInd(NoninducedFourCounts nonind)
{
    InducedFourCounts ret;

    ret.fourcliques = nonind.fourcliques; 
    ret.chordalcycles = nonind.chordalcycles - 6*ret.fourcliques;
    ret.fourcycles = nonind.fourcycles - 3*ret.fourcliques - ret.chordalcycles;
    ret.tailedtris = nonind.tailedtris - 4*ret.chordalcycles - 12*ret.fourcliques;
    ret.threepaths = nonind.threepaths - 2*ret.tailedtris - 4*ret.fourcycles - 6*ret.chordalcycles - 12*ret.fourcliques;
    ret.threestars = nonind.threestars - ret.tailedtris - 2*ret.chordalcycles - 4*ret.fourcliques;

    return ret;
}

// This generates the closure matrix of 4-vertex counts
//
// Input: NoninducedFourCounts struct and InducedFourCounts struct, and 6X6 matrix
// Output: Nothing. The input 6X6 matrix will contain the closures, where i,j entry is the fraction of pattern i that induced pattern j

void closureMatrix(NoninducedFourCounts nonind, InducedFourCounts ind, double (&matrix)[6][6])
{
    EdgeIdx nonind_array[6], ind_array[6];

    nonind_array[0] = nonind.threestars;
    nonind_array[1] = nonind.threepaths;
    nonind_array[2] = nonind.tailedtris;
    nonind_array[3] = nonind.fourcycles;
    nonind_array[4] = nonind.chordalcycles;
    nonind_array[5] = nonind.fourcliques;

    ind_array[0] = ind.threestars;
    ind_array[1] = ind.threepaths;
    ind_array[2] = ind.tailedtris;
    ind_array[3] = ind.fourcycles;
    ind_array[4] = ind.chordalcycles;
    ind_array[5] = ind.fourcliques;

    int conversion[6][6];

    conversion[0][0] = 1;
    conversion[0][1] = 0;
    conversion[0][2] = 1;
    conversion[0][3] = 0;
    conversion[0][4] = 2;
    conversion[0][5] = 4;
    conversion[1][0] = 0;
    conversion[1][1] = 1;
    conversion[1][2] = 2;
    conversion[1][3] = 4;
    conversion[1][4] = 6;
    conversion[1][5] = 12;
    conversion[2][0] = 0;
    conversion[2][1] = 0;
    conversion[2][2] = 1;
    conversion[2][3] = 0;
    conversion[2][4] = 4;
    conversion[2][5] = 12;
    conversion[3][0] = 0;
    conversion[3][1] = 0;
    conversion[3][2] = 0;
    conversion[3][3] = 1;
    conversion[3][4] = 1;
    conversion[3][5] = 3;
    conversion[4][0] = 0;
    conversion[4][1] = 0;
    conversion[4][2] = 0;
    conversion[4][3] = 0;
    conversion[4][4] = 1;
    conversion[4][5] = 6;
    conversion[5][0] = 0;
    conversion[5][1] = 0;
    conversion[5][2] = 0;
    conversion[5][3] = 0;
    conversion[5][4] = 0;
    conversion[5][5] = 1;


    for (int i=0; i<6; i++)
        for (int j=0; j < 6; j++)
            matrix[i][j] = ((double)conversion[i][j])*((double)ind_array[j])/((double)nonind_array[i]);
}

// This creates the induced counts


void getFiveInduced(EdgeIdx (&nonind)[21], EdgeIdx (&ind)[21])
{
    int inverse[21][21];

    inverse[0][0] = 1;
    inverse[0][1] = 0;
    inverse[0][2] = 0;
    inverse[0][3] = -1;
    inverse[0][4] = 0;
    inverse[0][5] = 0;
    inverse[0][6] = 0;
    inverse[0][7] = 0;
    inverse[0][8] = 1;
    inverse[0][9] = 0;
    inverse[0][10] = 1;
    inverse[0][11] = 0;
    inverse[0][12] = 0;
    inverse[0][13] = -2;
    inverse[0][14] = -1;
    inverse[0][15] = -1;
    inverse[0][16] = 0;
    inverse[0][17] = 1;
    inverse[0][18] = 2;
    inverse[0][19] = -3;
    inverse[0][20] = 5;
    inverse[1][0] = 0;
    inverse[1][1] = 1;
    inverse[1][2] = 0;
    inverse[1][3] = -2;
    inverse[1][4] = -1;
    inverse[1][5] = -2;
    inverse[1][6] = -2;
    inverse[1][7] = 0;
    inverse[1][8] = 4;
    inverse[1][9] = 4;
    inverse[1][10] = 5;
    inverse[1][11] = 4;
    inverse[1][12] = 6;
    inverse[1][13] = -12;
    inverse[1][14] = -9;
    inverse[1][15] = -10;
    inverse[1][16] = -10;
    inverse[1][17] = 20;
    inverse[1][18] = 20;
    inverse[1][19] = -36;
    inverse[1][20] = 60;
    inverse[2][0] = 0;
    inverse[2][1] = 0;
    inverse[2][2] = 1;
    inverse[2][3] = 0;
    inverse[2][4] = -2;
    inverse[2][5] = -1;
    inverse[2][6] = -2;
    inverse[2][7] = -5;
    inverse[2][8] = 4;
    inverse[2][9] = 4;
    inverse[2][10] = 2;
    inverse[2][11] = 7;
    inverse[2][12] = 6;
    inverse[2][13] = -6;
    inverse[2][14] = -6;
    inverse[2][15] = -10;
    inverse[2][16] = -14;
    inverse[2][17] = 24;
    inverse[2][18] = 18;
    inverse[2][19] = -36;
    inverse[2][20] = 60;
    inverse[3][0] = 0;
    inverse[3][1] = 0;
    inverse[3][2] = 0;
    inverse[3][3] = 1;
    inverse[3][4] = 0;
    inverse[3][5] = 0;
    inverse[3][6] = 0;
    inverse[3][7] = 0;
    inverse[3][8] = -2;
    inverse[3][9] = 0;
    inverse[3][10] = -2;
    inverse[3][11] = 0;
    inverse[3][12] = 0;
    inverse[3][13] = 6;
    inverse[3][14] = 3;
    inverse[3][15] = 3;
    inverse[3][16] = 0;
    inverse[3][17] = -4;
    inverse[3][18] = -8;
    inverse[3][19] = 15;
    inverse[3][20] = -30;
    inverse[4][0] = 0;
    inverse[4][1] = 0;
    inverse[4][2] = 0;
    inverse[4][3] = 0;
    inverse[4][4] = 1;
    inverse[4][5] = 0;
    inverse[4][6] = 0;
    inverse[4][7] = 0;
    inverse[4][8] = -4;
    inverse[4][9] = -2;
    inverse[4][10] = 0;
    inverse[4][11] = -2;
    inverse[4][12] = 0;
    inverse[4][13] = 0;
    inverse[4][14] = 3;
    inverse[4][15] = 6;
    inverse[4][16] = 6;
    inverse[4][17] = -16;
    inverse[4][18] = -12;
    inverse[4][19] = 30;
    inverse[4][20] = -60;
    inverse[5][0] = 0;
    inverse[5][1] = 0;
    inverse[5][2] = 0;
    inverse[5][3] = 0;
    inverse[5][4] = 0;
    inverse[5][5] = 1;
    inverse[5][6] = 0;
    inverse[5][7] = 0;
    inverse[5][8] = 0;
    inverse[5][9] = -2;
    inverse[5][10] = -2;
    inverse[5][11] = -1;
    inverse[5][12] = 0;
    inverse[5][13] = 6;
    inverse[5][14] = 6;
    inverse[5][15] = 5;
    inverse[5][16] = 4;
    inverse[5][17] = -12;
    inverse[5][18] = -14;
    inverse[5][19] = 30;
    inverse[5][20] = -60;
    inverse[6][0] = 0;
    inverse[6][1] = 0;
    inverse[6][2] = 0;
    inverse[6][3] = 0;
    inverse[6][4] = 0;
    inverse[6][5] = 0;
    inverse[6][6] = 1;
    inverse[6][7] = 0;
    inverse[6][8] = 0;
    inverse[6][9] = -1;
    inverse[6][10] = -1;
    inverse[6][11] = -2;
    inverse[6][12] = -6;
    inverse[6][13] = 6;
    inverse[6][14] = 3;
    inverse[6][15] = 4;
    inverse[6][16] = 8;
    inverse[6][17] = -16;
    inverse[6][18] = -12;
    inverse[6][19] = 30;
    inverse[6][20] = -60;
    inverse[7][0] = 0;
    inverse[7][1] = 0;
    inverse[7][2] = 0;
    inverse[7][3] = 0;
    inverse[7][4] = 0;
    inverse[7][5] = 0;
    inverse[7][6] = 0;
    inverse[7][7] = 1;
    inverse[7][8] = 0;
    inverse[7][9] = 0;
    inverse[7][10] = 0;
    inverse[7][11] = -1;
    inverse[7][12] = 0;
    inverse[7][13] = 0;
    inverse[7][14] = 0;
    inverse[7][15] = 1;
    inverse[7][16] = 2;
    inverse[7][17] = -4;
    inverse[7][18] = -2;
    inverse[7][19] = 6;
    inverse[7][20] = -12;
    inverse[8][0] = 0;
    inverse[8][1] = 0;
    inverse[8][2] = 0;
    inverse[8][3] = 0;
    inverse[8][4] = 0;
    inverse[8][5] = 0;
    inverse[8][6] = 0;
    inverse[8][7] = 0;
    inverse[8][8] = 1;
    inverse[8][9] = 0;
    inverse[8][10] = 0;
    inverse[8][11] = 0;
    inverse[8][12] = 0;
    inverse[8][13] = 0;
    inverse[8][14] = 0;
    inverse[8][15] = -1;
    inverse[8][16] = 0;
    inverse[8][17] = 2;
    inverse[8][18] = 2;
    inverse[8][19] = -6;
    inverse[8][20] = 15;
    inverse[9][0] = 0;
    inverse[9][1] = 0;
    inverse[9][2] = 0;
    inverse[9][3] = 0;
    inverse[9][4] = 0;
    inverse[9][5] = 0;
    inverse[9][6] = 0;
    inverse[9][7] = 0;
    inverse[9][8] = 0;
    inverse[9][9] = 1;
    inverse[9][10] = 0;
    inverse[9][11] = 0;
    inverse[9][12] = 0;
    inverse[9][13] = 0;
    inverse[9][14] = -3;
    inverse[9][15] = -2;
    inverse[9][16] = -2;
    inverse[9][17] = 8;
    inverse[9][18] = 8;
    inverse[9][19] = -24;
    inverse[9][20] = 60;
    inverse[10][0] = 0;
    inverse[10][1] = 0;
    inverse[10][2] = 0;
    inverse[10][3] = 0;
    inverse[10][4] = 0;
    inverse[10][5] = 0;
    inverse[10][6] = 0;
    inverse[10][7] = 0;
    inverse[10][8] = 0;
    inverse[10][9] = 0;
    inverse[10][10] = 1;
    inverse[10][11] = 0;
    inverse[10][12] = 0;
    inverse[10][13] = -6;
    inverse[10][14] = -3;
    inverse[10][15] = -2;
    inverse[10][16] = 0;
    inverse[10][17] = 4;
    inverse[10][18] = 10;
    inverse[10][19] = -24;
    inverse[10][20] = 60;
    inverse[11][0] = 0;
    inverse[11][1] = 0;
    inverse[11][2] = 0;
    inverse[11][3] = 0;
    inverse[11][4] = 0;
    inverse[11][5] = 0;
    inverse[11][6] = 0;
    inverse[11][7] = 0;
    inverse[11][8] = 0;
    inverse[11][9] = 0;
    inverse[11][10] = 0;
    inverse[11][11] = 1;
    inverse[11][12] = 0;
    inverse[11][13] = 0;
    inverse[11][14] = 0;
    inverse[11][15] = -2;
    inverse[11][16] = -4;
    inverse[11][17] = 12;
    inverse[11][18] = 6;
    inverse[11][19] = -24;
    inverse[11][20] = 60;
    inverse[12][0] = 0;
    inverse[12][1] = 0;
    inverse[12][2] = 0;
    inverse[12][3] = 0;
    inverse[12][4] = 0;
    inverse[12][5] = 0;
    inverse[12][6] = 0;
    inverse[12][7] = 0;
    inverse[12][8] = 0;
    inverse[12][9] = 0;
    inverse[12][10] = 0;
    inverse[12][11] = 0;
    inverse[12][12] = 1;
    inverse[12][13] = -1;
    inverse[12][14] = 0;
    inverse[12][15] = 0;
    inverse[12][16] = -1;
    inverse[12][17] = 2;
    inverse[12][18] = 1;
    inverse[12][19] = -4;
    inverse[12][20] = 10;
    inverse[13][0] = 0;
    inverse[13][1] = 0;
    inverse[13][2] = 0;
    inverse[13][3] = 0;
    inverse[13][4] = 0;
    inverse[13][5] = 0;
    inverse[13][6] = 0;
    inverse[13][7] = 0;
    inverse[13][8] = 0;
    inverse[13][9] = 0;
    inverse[13][10] = 0;
    inverse[13][11] = 0;
    inverse[13][12] = 0;
    inverse[13][13] = 1;
    inverse[13][14] = 0;
    inverse[13][15] = 0;
    inverse[13][16] = 0;
    inverse[13][17] = 0;
    inverse[13][18] = -1;
    inverse[13][19] = 3;
    inverse[13][20] = -10;
    inverse[14][0] = 0;
    inverse[14][1] = 0;
    inverse[14][2] = 0;
    inverse[14][3] = 0;
    inverse[14][4] = 0;
    inverse[14][5] = 0;
    inverse[14][6] = 0;
    inverse[14][7] = 0;
    inverse[14][8] = 0;
    inverse[14][9] = 0;
    inverse[14][10] = 0;
    inverse[14][11] = 0;
    inverse[14][12] = 0;
    inverse[14][13] = 0;
    inverse[14][14] = 1;
    inverse[14][15] = 0;
    inverse[14][16] = 0;
    inverse[14][17] = 0;
    inverse[14][18] = -2;
    inverse[14][19] = 6;
    inverse[14][20] = -20;
    inverse[15][0] = 0;
    inverse[15][1] = 0;
    inverse[15][2] = 0;
    inverse[15][3] = 0;
    inverse[15][4] = 0;
    inverse[15][5] = 0;
    inverse[15][6] = 0;
    inverse[15][7] = 0;
    inverse[15][8] = 0;
    inverse[15][9] = 0;
    inverse[15][10] = 0;
    inverse[15][11] = 0;
    inverse[15][12] = 0;
    inverse[15][13] = 0;
    inverse[15][14] = 0;
    inverse[15][15] = 1;
    inverse[15][16] = 0;
    inverse[15][17] = -4;
    inverse[15][18] = -4;
    inverse[15][19] = 18;
    inverse[15][20] = -60;
    inverse[16][0] = 0;
    inverse[16][1] = 0;
    inverse[16][2] = 0;
    inverse[16][3] = 0;
    inverse[16][4] = 0;
    inverse[16][5] = 0;
    inverse[16][6] = 0;
    inverse[16][7] = 0;
    inverse[16][8] = 0;
    inverse[16][9] = 0;
    inverse[16][10] = 0;
    inverse[16][11] = 0;
    inverse[16][12] = 0;
    inverse[16][13] = 0;
    inverse[16][14] = 0;
    inverse[16][15] = 0;
    inverse[16][16] = 1;
    inverse[16][17] = -4;
    inverse[16][18] = -1;
    inverse[16][19] = 9;
    inverse[16][20] = -30;
    inverse[17][0] = 0;
    inverse[17][1] = 0;
    inverse[17][2] = 0;
    inverse[17][3] = 0;
    inverse[17][4] = 0;
    inverse[17][5] = 0;
    inverse[17][6] = 0;
    inverse[17][7] = 0;
    inverse[17][8] = 0;
    inverse[17][9] = 0;
    inverse[17][10] = 0;
    inverse[17][11] = 0;
    inverse[17][12] = 0;
    inverse[17][13] = 0;
    inverse[17][14] = 0;
    inverse[17][15] = 0;
    inverse[17][16] = 0;
    inverse[17][17] = 1;
    inverse[17][18] = 0;
    inverse[17][19] = -3;
    inverse[17][20] = 15;
    inverse[18][0] = 0;
    inverse[18][1] = 0;
    inverse[18][2] = 0;
    inverse[18][3] = 0;
    inverse[18][4] = 0;
    inverse[18][5] = 0;
    inverse[18][6] = 0;
    inverse[18][7] = 0;
    inverse[18][8] = 0;
    inverse[18][9] = 0;
    inverse[18][10] = 0;
    inverse[18][11] = 0;
    inverse[18][12] = 0;
    inverse[18][13] = 0;
    inverse[18][14] = 0;
    inverse[18][15] = 0;
    inverse[18][16] = 0;
    inverse[18][17] = 0;
    inverse[18][18] = 1;
    inverse[18][19] = -6;
    inverse[18][20] = 30;
    inverse[19][0] = 0;
    inverse[19][1] = 0;
    inverse[19][2] = 0;
    inverse[19][3] = 0;
    inverse[19][4] = 0;
    inverse[19][5] = 0;
    inverse[19][6] = 0;
    inverse[19][7] = 0;
    inverse[19][8] = 0;
    inverse[19][9] = 0;
    inverse[19][10] = 0;
    inverse[19][11] = 0;
    inverse[19][12] = 0;
    inverse[19][13] = 0;
    inverse[19][14] = 0;
    inverse[19][15] = 0;
    inverse[19][16] = 0;
    inverse[19][17] = 0;
    inverse[19][18] = 0;
    inverse[19][19] = 1;
    inverse[19][20] = -10;
    inverse[20][0] = 0;
    inverse[20][1] = 0;
    inverse[20][2] = 0;
    inverse[20][3] = 0;
    inverse[20][4] = 0;
    inverse[20][5] = 0;
    inverse[20][6] = 0;
    inverse[20][7] = 0;
    inverse[20][8] = 0;
    inverse[20][9] = 0;
    inverse[20][10] = 0;
    inverse[20][11] = 0;
    inverse[20][12] = 0;
    inverse[20][13] = 0;
    inverse[20][14] = 0;
    inverse[20][15] = 0;
    inverse[20][16] = 0;
    inverse[20][17] = 0;
    inverse[20][18] = 0;
    inverse[20][19] = 0;
    inverse[20][20] = 1;



    for (int i=0; i<21; i++)
    {
        ind[i] = 0;
        for (int j=0; j<21; j++)
            ind[i] += inverse[i][j]*nonind[j];
    }
}




// This generates the closure matrix of 5-vertex counts
//
// Input: Arrays of noninduced and induced counts, and 21X21 matrix
// Output: Nothing. The input 21X21 matrix will contain the closures, where i,j entry is the fraction of pattern i that induced pattern j

void closureMatrixFive(EdgeIdx (&nonind)[21], EdgeIdx (&ind)[21], double (&matrix)[21][21])
{
    int conversion[21][21];

    conversion[0][0] = 1;
    conversion[0][1] = 0;
    conversion[0][2] = 0;
    conversion[0][3] = 1;
    conversion[0][4] = 0;
    conversion[0][5] = 0;
    conversion[0][6] = 0;
    conversion[0][7] = 0;
    conversion[0][8] = 1;
    conversion[0][9] = 0;
    conversion[0][10] = 1;
    conversion[0][11] = 0;
    conversion[0][12] = 0;
    conversion[0][13] = 2;
    conversion[0][14] = 1;
    conversion[0][15] = 1;
    conversion[0][16] = 0;
    conversion[0][17] = 1;
    conversion[0][18] = 2;
    conversion[0][19] = 3;
    conversion[0][20] = 5;
    conversion[1][0] = 0;
    conversion[1][1] = 1;
    conversion[1][2] = 0;
    conversion[1][3] = 2;
    conversion[1][4] = 1;
    conversion[1][5] = 2;
    conversion[1][6] = 2;
    conversion[1][7] = 0;
    conversion[1][8] = 4;
    conversion[1][9] = 4;
    conversion[1][10] = 5;
    conversion[1][11] = 4;
    conversion[1][12] = 6;
    conversion[1][13] = 12;
    conversion[1][14] = 9;
    conversion[1][15] = 10;
    conversion[1][16] = 10;
    conversion[1][17] = 20;
    conversion[1][18] = 20;
    conversion[1][19] = 36;
    conversion[1][20] = 60;
    conversion[2][0] = 0;
    conversion[2][1] = 0;
    conversion[2][2] = 1;
    conversion[2][3] = 0;
    conversion[2][4] = 2;
    conversion[2][5] = 1;
    conversion[2][6] = 2;
    conversion[2][7] = 5;
    conversion[2][8] = 4;
    conversion[2][9] = 4;
    conversion[2][10] = 2;
    conversion[2][11] = 7;
    conversion[2][12] = 6;
    conversion[2][13] = 6;
    conversion[2][14] = 6;
    conversion[2][15] = 10;
    conversion[2][16] = 14;
    conversion[2][17] = 24;
    conversion[2][18] = 18;
    conversion[2][19] = 36;
    conversion[2][20] = 60;
    conversion[3][0] = 0;
    conversion[3][1] = 0;
    conversion[3][2] = 0;
    conversion[3][3] = 1;
    conversion[3][4] = 0;
    conversion[3][5] = 0;
    conversion[3][6] = 0;
    conversion[3][7] = 0;
    conversion[3][8] = 2;
    conversion[3][9] = 0;
    conversion[3][10] = 2;
    conversion[3][11] = 0;
    conversion[3][12] = 0;
    conversion[3][13] = 6;
    conversion[3][14] = 3;
    conversion[3][15] = 3;
    conversion[3][16] = 0;
    conversion[3][17] = 4;
    conversion[3][18] = 8;
    conversion[3][19] = 15;
    conversion[3][20] = 30;
    conversion[4][0] = 0;
    conversion[4][1] = 0;
    conversion[4][2] = 0;
    conversion[4][3] = 0;
    conversion[4][4] = 1;
    conversion[4][5] = 0;
    conversion[4][6] = 0;
    conversion[4][7] = 0;
    conversion[4][8] = 4;
    conversion[4][9] = 2;
    conversion[4][10] = 0;
    conversion[4][11] = 2;
    conversion[4][12] = 0;
    conversion[4][13] = 0;
    conversion[4][14] = 3;
    conversion[4][15] = 6;
    conversion[4][16] = 6;
    conversion[4][17] = 16;
    conversion[4][18] = 12;
    conversion[4][19] = 30;
    conversion[4][20] = 60;
    conversion[5][0] = 0;
    conversion[5][1] = 0;
    conversion[5][2] = 0;
    conversion[5][3] = 0;
    conversion[5][4] = 0;
    conversion[5][5] = 1;
    conversion[5][6] = 0;
    conversion[5][7] = 0;
    conversion[5][8] = 0;
    conversion[5][9] = 2;
    conversion[5][10] = 2;
    conversion[5][11] = 1;
    conversion[5][12] = 0;
    conversion[5][13] = 6;
    conversion[5][14] = 6;
    conversion[5][15] = 5;
    conversion[5][16] = 4;
    conversion[5][17] = 12;
    conversion[5][18] = 14;
    conversion[5][19] = 30;
    conversion[5][20] = 60;
    conversion[6][0] = 0;
    conversion[6][1] = 0;
    conversion[6][2] = 0;
    conversion[6][3] = 0;
    conversion[6][4] = 0;
    conversion[6][5] = 0;
    conversion[6][6] = 1;
    conversion[6][7] = 0;
    conversion[6][8] = 0;
    conversion[6][9] = 1;
    conversion[6][10] = 1;
    conversion[6][11] = 2;
    conversion[6][12] = 6;
    conversion[6][13] = 6;
    conversion[6][14] = 3;
    conversion[6][15] = 4;
    conversion[6][16] = 8;
    conversion[6][17] = 16;
    conversion[6][18] = 12;
    conversion[6][19] = 30;
    conversion[6][20] = 60;
    conversion[7][0] = 0;
    conversion[7][1] = 0;
    conversion[7][2] = 0;
    conversion[7][3] = 0;
    conversion[7][4] = 0;
    conversion[7][5] = 0;
    conversion[7][6] = 0;
    conversion[7][7] = 1;
    conversion[7][8] = 0;
    conversion[7][9] = 0;
    conversion[7][10] = 0;
    conversion[7][11] = 1;
    conversion[7][12] = 0;
    conversion[7][13] = 0;
    conversion[7][14] = 0;
    conversion[7][15] = 1;
    conversion[7][16] = 2;
    conversion[7][17] = 4;
    conversion[7][18] = 2;
    conversion[7][19] = 6;
    conversion[7][20] = 12;
    conversion[8][0] = 0;
    conversion[8][1] = 0;
    conversion[8][2] = 0;
    conversion[8][3] = 0;
    conversion[8][4] = 0;
    conversion[8][5] = 0;
    conversion[8][6] = 0;
    conversion[8][7] = 0;
    conversion[8][8] = 1;
    conversion[8][9] = 0;
    conversion[8][10] = 0;
    conversion[8][11] = 0;
    conversion[8][12] = 0;
    conversion[8][13] = 0;
    conversion[8][14] = 0;
    conversion[8][15] = 1;
    conversion[8][16] = 0;
    conversion[8][17] = 2;
    conversion[8][18] = 2;
    conversion[8][19] = 6;
    conversion[8][20] = 15;
    conversion[9][0] = 0;
    conversion[9][1] = 0;
    conversion[9][2] = 0;
    conversion[9][3] = 0;
    conversion[9][4] = 0;
    conversion[9][5] = 0;
    conversion[9][6] = 0;
    conversion[9][7] = 0;
    conversion[9][8] = 0;
    conversion[9][9] = 1;
    conversion[9][10] = 0;
    conversion[9][11] = 0;
    conversion[9][12] = 0;
    conversion[9][13] = 0;
    conversion[9][14] = 3;
    conversion[9][15] = 2;
    conversion[9][16] = 2;
    conversion[9][17] = 8;
    conversion[9][18] = 8;
    conversion[9][19] = 24;
    conversion[9][20] = 60;
    conversion[10][0] = 0;
    conversion[10][1] = 0;
    conversion[10][2] = 0;
    conversion[10][3] = 0;
    conversion[10][4] = 0;
    conversion[10][5] = 0;
    conversion[10][6] = 0;
    conversion[10][7] = 0;
    conversion[10][8] = 0;
    conversion[10][9] = 0;
    conversion[10][10] = 1;
    conversion[10][11] = 0;
    conversion[10][12] = 0;
    conversion[10][13] = 6;
    conversion[10][14] = 3;
    conversion[10][15] = 2;
    conversion[10][16] = 0;
    conversion[10][17] = 4;
    conversion[10][18] = 10;
    conversion[10][19] = 24;
    conversion[10][20] = 60;
    conversion[11][0] = 0;
    conversion[11][1] = 0;
    conversion[11][2] = 0;
    conversion[11][3] = 0;
    conversion[11][4] = 0;
    conversion[11][5] = 0;
    conversion[11][6] = 0;
    conversion[11][7] = 0;
    conversion[11][8] = 0;
    conversion[11][9] = 0;
    conversion[11][10] = 0;
    conversion[11][11] = 1;
    conversion[11][12] = 0;
    conversion[11][13] = 0;
    conversion[11][14] = 0;
    conversion[11][15] = 2;
    conversion[11][16] = 4;
    conversion[11][17] = 12;
    conversion[11][18] = 6;
    conversion[11][19] = 24;
    conversion[11][20] = 60;
    conversion[12][0] = 0;
    conversion[12][1] = 0;
    conversion[12][2] = 0;
    conversion[12][3] = 0;
    conversion[12][4] = 0;
    conversion[12][5] = 0;
    conversion[12][6] = 0;
    conversion[12][7] = 0;
    conversion[12][8] = 0;
    conversion[12][9] = 0;
    conversion[12][10] = 0;
    conversion[12][11] = 0;
    conversion[12][12] = 1;
    conversion[12][13] = 1;
    conversion[12][14] = 0;
    conversion[12][15] = 0;
    conversion[12][16] = 1;
    conversion[12][17] = 2;
    conversion[12][18] = 1;
    conversion[12][19] = 4;
    conversion[12][20] = 10;
    conversion[13][0] = 0;
    conversion[13][1] = 0;
    conversion[13][2] = 0;
    conversion[13][3] = 0;
    conversion[13][4] = 0;
    conversion[13][5] = 0;
    conversion[13][6] = 0;
    conversion[13][7] = 0;
    conversion[13][8] = 0;
    conversion[13][9] = 0;
    conversion[13][10] = 0;
    conversion[13][11] = 0;
    conversion[13][12] = 0;
    conversion[13][13] = 1;
    conversion[13][14] = 0;
    conversion[13][15] = 0;
    conversion[13][16] = 0;
    conversion[13][17] = 0;
    conversion[13][18] = 1;
    conversion[13][19] = 3;
    conversion[13][20] = 10;
    conversion[14][0] = 0;
    conversion[14][1] = 0;
    conversion[14][2] = 0;
    conversion[14][3] = 0;
    conversion[14][4] = 0;
    conversion[14][5] = 0;
    conversion[14][6] = 0;
    conversion[14][7] = 0;
    conversion[14][8] = 0;
    conversion[14][9] = 0;
    conversion[14][10] = 0;
    conversion[14][11] = 0;
    conversion[14][12] = 0;
    conversion[14][13] = 0;
    conversion[14][14] = 1;
    conversion[14][15] = 0;
    conversion[14][16] = 0;
    conversion[14][17] = 0;
    conversion[14][18] = 2;
    conversion[14][19] = 6;
    conversion[14][20] = 20;
    conversion[15][0] = 0;
    conversion[15][1] = 0;
    conversion[15][2] = 0;
    conversion[15][3] = 0;
    conversion[15][4] = 0;
    conversion[15][5] = 0;
    conversion[15][6] = 0;
    conversion[15][7] = 0;
    conversion[15][8] = 0;
    conversion[15][9] = 0;
    conversion[15][10] = 0;
    conversion[15][11] = 0;
    conversion[15][12] = 0;
    conversion[15][13] = 0;
    conversion[15][14] = 0;
    conversion[15][15] = 1;
    conversion[15][16] = 0;
    conversion[15][17] = 4;
    conversion[15][18] = 4;
    conversion[15][19] = 18;
    conversion[15][20] = 60;
    conversion[16][0] = 0;
    conversion[16][1] = 0;
    conversion[16][2] = 0;
    conversion[16][3] = 0;
    conversion[16][4] = 0;
    conversion[16][5] = 0;
    conversion[16][6] = 0;
    conversion[16][7] = 0;
    conversion[16][8] = 0;
    conversion[16][9] = 0;
    conversion[16][10] = 0;
    conversion[16][11] = 0;
    conversion[16][12] = 0;
    conversion[16][13] = 0;
    conversion[16][14] = 0;
    conversion[16][15] = 0;
    conversion[16][16] = 1;
    conversion[16][17] = 4;
    conversion[16][18] = 1;
    conversion[16][19] = 9;
    conversion[16][20] = 30;
    conversion[17][0] = 0;
    conversion[17][1] = 0;
    conversion[17][2] = 0;
    conversion[17][3] = 0;
    conversion[17][4] = 0;
    conversion[17][5] = 0;
    conversion[17][6] = 0;
    conversion[17][7] = 0;
    conversion[17][8] = 0;
    conversion[17][9] = 0;
    conversion[17][10] = 0;
    conversion[17][11] = 0;
    conversion[17][12] = 0;
    conversion[17][13] = 0;
    conversion[17][14] = 0;
    conversion[17][15] = 0;
    conversion[17][16] = 0;
    conversion[17][17] = 1;
    conversion[17][18] = 0;
    conversion[17][19] = 3;
    conversion[17][20] = 15;
    conversion[18][0] = 0;
    conversion[18][1] = 0;
    conversion[18][2] = 0;
    conversion[18][3] = 0;
    conversion[18][4] = 0;
    conversion[18][5] = 0;
    conversion[18][6] = 0;
    conversion[18][7] = 0;
    conversion[18][8] = 0;
    conversion[18][9] = 0;
    conversion[18][10] = 0;
    conversion[18][11] = 0;
    conversion[18][12] = 0;
    conversion[18][13] = 0;
    conversion[18][14] = 0;
    conversion[18][15] = 0;
    conversion[18][16] = 0;
    conversion[18][17] = 0;
    conversion[18][18] = 1;
    conversion[18][19] = 6;
    conversion[18][20] = 30;
    conversion[19][0] = 0;
    conversion[19][1] = 0;
    conversion[19][2] = 0;
    conversion[19][3] = 0;
    conversion[19][4] = 0;
    conversion[19][5] = 0;
    conversion[19][6] = 0;
    conversion[19][7] = 0;
    conversion[19][8] = 0;
    conversion[19][9] = 0;
    conversion[19][10] = 0;
    conversion[19][11] = 0;
    conversion[19][12] = 0;
    conversion[19][13] = 0;
    conversion[19][14] = 0;
    conversion[19][15] = 0;
    conversion[19][16] = 0;
    conversion[19][17] = 0;
    conversion[19][18] = 0;
    conversion[19][19] = 1;
    conversion[19][20] = 10;
    conversion[20][0] = 0;
    conversion[20][1] = 0;
    conversion[20][2] = 0;
    conversion[20][3] = 0;
    conversion[20][4] = 0;
    conversion[20][5] = 0;
    conversion[20][6] = 0;
    conversion[20][7] = 0;
    conversion[20][8] = 0;
    conversion[20][9] = 0;
    conversion[20][10] = 0;
    conversion[20][11] = 0;
    conversion[20][12] = 0;
    conversion[20][13] = 0;
    conversion[20][14] = 0;
    conversion[20][15] = 0;
    conversion[20][16] = 0;
    conversion[20][17] = 0;
    conversion[20][18] = 0;
    conversion[20][19] = 0;
    conversion[20][20] = 1;






    for (int i=0; i<21; i++)
        for (int j=0; j<21; j++)
            matrix[i][j] = ((double)conversion[i][j])*((double)ind[j])/((double)nonind[i]);

}
#endif



