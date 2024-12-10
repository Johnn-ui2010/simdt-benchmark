import os
import sys
import subprocess
import time
if __name__ == '__main__':
    c_start = time.time()
    print("The current working directory is:", os.getcwd())

    #Compile the SGX part.
    subprocess.check_output(["make", "clean"], cwd="../")
    subprocess.check_output(["make", "SGX_PRERELEASE=1", "SGX_DEBUG=0", ], cwd="../", stderr=subprocess.DEVNULL)
    
    subprocess.check_output(["./app", "0" ], cwd="../")

    #Compile the native part.
    print("The current working directory is:", os.getcwd())

    subprocess.check_output(["make", "clean"], cwd="../")
    subprocess.check_output(["make", "native", ], cwd="../", stderr=subprocess.DEVNULL)
    
    subprocess.check_output(["./app", "1" ], cwd="../")
    c_end = time.time()
    print("total time: ", round(c_end-c_start, 3), " s")
