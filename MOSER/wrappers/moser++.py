from utils import run_command
from subgraph_counts import matrices, names
import random, argparse, numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(description="MOSER++")

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
    parser.add_argument(
        "-p", "--p-value", type=float, default=0.01, help="P-value (default: 0.01)"
    )

    args = parser.parse_args()

    return args


def calc_p_value_upper(arr, n):
    return sum(row[n] > arr[0][n] for row in arr)


def parse_non_induced_to_induced(non_induced, motif_size):
    return list(np.linalg.solve(matrices[motif_size], non_induced))


def parse_ATAC_output(motif_size):
    f = open("out.txt", "r")
    lines = f.readlines()
    f.close()

    num_steps = int(lines[1].strip())
    lines = lines[4:]
    res = []
    for i in range(motif_size - 2):
        res.append([])
        for j in range(num_steps):
            non_induced = [int(float(s)) for s in lines[j].strip().split(" ")]
            res[i].append(parse_non_induced_to_induced(non_induced, i + 3))
        lines = lines[num_steps + 1 :]
    return res


def calc_p_value(res, motif_size):
    uppers = [[] for _ in range(motif_size - 2)]
    for i in range(motif_size - 2):
        for j in range(len(res[i][0])):
            uppers[i].append(calc_p_value_upper(res[i], j))
    return uppers


def serial_test(args):
    pivot = random.randint(1, args.num_steps)

    cmd_1 = f"../exe/ATAC{args.motif_size} {args.graph} {pivot}"
    cmd_2 = f"../exe/ATAC{args.motif_size} {args.graph} {args.num_steps - pivot}"

    print(cmd_1)
    t1 = run_command(cmd_1)
    r1 = parse_ATAC_output(args.motif_size)
    uppers_1 = calc_p_value(r1, args.motif_size)
    # print(uppers_1)

    print(cmd_2)
    t2 = run_command(cmd_2)
    r2 = parse_ATAC_output(args.motif_size)
    uppers_2 = calc_p_value(r2, args.motif_size)
    # print(uppers_2)

    for i in range(args.motif_size - 2):
        summation = result = [
            a + b + c
            for a, b, c in zip(
                uppers_1[i], uppers_2[i], [1 for _ in range(len(uppers_1[i]))]
            )
        ]
        result = [x / args.num_steps for x in summation]
        border = -2 if i == 0 else -6
        print(f"patterns of size {i+3}: {names[i+3][border:]}")
        print(f"p-value for motifs of size {i+3}: {result[border:]}")


def main():
    args = parse_arguments()
    serial_test(args)


if __name__ == "__main__":
    main()
