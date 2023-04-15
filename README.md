# MOSER
 Scalable Network Motif Discovery using Serial Test

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Notes](#notes)


## Installation

1. Clone the repository to your local machine using: 
 ```Bash
 git clone https://github.com/momatinaj/MOSER.git
 ```

2. Navigate to the project directory using:
 ```Bash
 cd MOSER
 ```
3. Run the command `make`, to complile the project and create the exectuables.

## Usage

1. To run the ESCAPE subgraph counting algorithm, navigate to the wrappers folder:
```Bash
 cd wrappers/
 ```
and run the following command:
```Bash
 python3 subgraph_counts.py <INPUT GRAPH PATH> <SUBGRAPH SIZE> <OPTIONAL FLAGS>
 ```
For Example:
`python3 subgraph_counts.py ../graphs/ca-AstroPh.edges 4 -i`.

2. To run the MOSER+ or MOSER++ motif discovery algorithm, navigate to the wrappers folder and run the following command:
```Bash
 python3 moser+.py (or moser++.py) -g <ORIGINAL GRAPH PATH> -s <SUBGRAPH SIZE> -n <NUMBER OF STEPS>
 ```
For Example:
`python3 moser++.py -g ../graphs/ca-AstroPh.edges -s 4 -n 1000`.


## Notes

- SUBGRAPH SIZE = 3, 4, 5.

- OPTIONAL FLAGS: (-i)output counts as integers. Useful for small graphs, or for debugging.


