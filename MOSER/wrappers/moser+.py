"""

The file moser+.py is a wrapper for the MOSER+ program.
At every step, one single switch is performed on G and a bfs is performed over the affected area.
The change in the number of motifs is recorded.

"""
import random, argparse, numpy as np, networkx as nx, copy
from utils import run_command
from subgraph_counts import matrices, names
import moser_graph as graph

number_to_string = {3: "three", 4: "four", 5: "five"}


def parse_arguments():
    parser = argparse.ArgumentParser(description="MOSER+")

    parser.add_argument(
        "-g", "--graph", type=str, required=True, help="Path to the input graph"
    )
    parser.add_argument(
        "-s", "--motif-size", type=int, required=True, help="Motif size"
    )
    parser.add_argument(
        "-n",
        "--num-steps",
        type=int,
        default=10000,
        help="Number of steps (default: 10000)",
    )

    args = parser.parse_args()

    return args


def calc_p_value_upper(arr, n):
    return sum(row[n] > arr[0][n] for row in arr)


def calc_p_value(res, motif_size):
    uppers = [[] for _ in range(motif_size - 2)]
    for i in range(motif_size - 2):
        for j in range(len(res[i][0])):
            uppers[i].append(calc_p_value_upper(res[i], j))
    return uppers


def get_esc_count(graph_path, motif_size):
    command = (
        f"../exe/count_{number_to_string[motif_size]} {graph_path} {motif_size} -i"
    )
    output = run_command(command)
    motif_counts = [float(s.strip()) for s in open("out.txt").readlines()[2:]]
    return motif_counts


def single_step_track(g: graph.RandomGraph, e1, e2, motif_size):
    temp_file = "temp_esc"
    depth = motif_size - 2
    src1, dst1 = e1
    src2, dst2 = e2
    e3 = (src1, dst2)
    e4 = (src2, dst1)

    n1 = g.get_area_around_the_edge(e1, depth, True)
    n2 = g.get_area_around_the_edge(e2, depth, True)
    n3 = g.get_area_around_the_edge(e3, depth, True)
    n4 = g.get_area_around_the_edge(e4, depth, True)

    combined = list(set(n1 + n2 + n3 + n4))
    compressed_g = g.G.subgraph(list(combined))
    graph.save_graph(compressed_g, temp_file, form="ESCAPE")
    before = get_esc_count(temp_file, motif_size)

    after_removal = nx.Graph(compressed_g)
    after_removal.remove_edge(src1, dst1)
    after_removal.remove_edge(src2, dst2)
    after_removal.add_edge(src1, dst2)
    after_removal.add_edge(src2, dst1)

    graph.save_graph(after_removal, temp_file, form="ESCAPE")
    after = get_esc_count(temp_file, motif_size)

    delta = [a - b for a, b in zip(after, before)]

    return delta


def parse_non_induced_to_induced(non_induced, motif_size):
    return list(np.linalg.solve(matrices[motif_size], non_induced))


def result_parser(res, motif_size):
    if motif_size == 3:
        result = [[parse_non_induced_to_induced(r, 3) for r in res]]
    elif motif_size == 4:
        r_3 = [parse_non_induced_to_induced(r[:4], 3) for r in res]
        r_4 = [parse_non_induced_to_induced(r[4:], 4) for r in res]
        result = [r_3, r_4]
    elif motif_size == 5:
        r_3 = [parse_non_induced_to_induced(r[:4], 3) for r in res]
        r_4 = [parse_non_induced_to_induced(r[4:15], 4) for r in res]
        r_5 = [parse_non_induced_to_induced(r[15:], 5) for r in res]
        result = [r_3, r_4, r_5]
    return result


def full_trajecory(G, motif_size, num_steps, original_counts):
    result = [original_counts]
    for _ in range(num_steps - 1):
        # find a switch candidate
        e1, e2 = G.find_switch_candidates()
        if None in e1 or None in e2:
            continue
        # perform a single step track
        delta = single_step_track(G, e1, e2, motif_size)
        result.append([delta[i] + result[-1][i] for i in range(len(delta))])

        # perform the switch
        G.G.remove_edge(e1[0], e1[1])
        G.G.remove_edge(e2[0], e2[1])
        G.G.add_edge(e1[0], e2[1])
        G.G.add_edge(e2[0], e1[1])

    # replace every non-induced count with the induced count in result (handle different count cases)
    result = result_parser(result, motif_size)
    return result



def main():
    args = parse_arguments()
    G = graph.RandomGraph(args.graph, False)
    original_counts = get_esc_count(args.graph, args.motif_size)
    pivot = random.randint(1, args.num_steps)
    r1 = full_trajecory(copy.deepcopy(G), args.motif_size, pivot, original_counts)
    r2 = full_trajecory(
        copy.deepcopy(G), args.motif_size, args.num_steps - pivot, original_counts
    )
    uppers_1 = calc_p_value(r1, args.motif_size)
    uppers_2 = calc_p_value(r2, args.motif_size)
    for i in range(args.motif_size - 2):
        summation = [
            a + b + c
            for a, b, c in zip(
                uppers_1[i], uppers_2[i], [1 for _ in range(len(uppers_1[i]))]
            )
        ]
        result = [x / args.num_steps for x in summation]
        border = -2 if i == 0 else -6
        border = -21 if i == 2 else border
        print(f"patterns of size {i+3}: {names[i+3][border:]}")
        print(f"p-value for motifs of size {i+3}: {result[border:]}")

    # save the result
    with open(f"mp_out.txt", "w") as f:
        f.write(f"trajectory 1:\n")
        for j in range(args.motif_size - 2):
            for r in r1[j]:
                f.write(" ".join([str(i) for i in r]) + "\n")
        f.write(f"trajectory 2:\n")
        for j in range(args.motif_size - 2):
            for r in r2[j]:
                f.write(" ".join([str(i) for i in r]) + "\n")


if __name__ == "__main__":
    main()
