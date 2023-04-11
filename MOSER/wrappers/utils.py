import subprocess
import psutil
import time


def run_command(command):
    result = {}
    start_time = time.time()
    experiment = subprocess.Popen(
        command,
        shell=True,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    pid = experiment.pid
    process = psutil.Process(pid)
    memory_usage = str(process.memory_info().rss / 1000) + " KB"

    out, err = experiment.communicate()
    if err != "":
        print(err)
    exec_time = time.time() - start_time
    result["output"] = str(out)
    result["time"] = str(exec_time)
    result["memory"] = memory_usage
    return result
