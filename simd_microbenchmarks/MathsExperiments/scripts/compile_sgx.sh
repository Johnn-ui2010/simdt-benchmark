#!/bin/bash
cd ..
make clean
make SGX_PRERELEASE=1 SGX_DEBUG=0
#make SGX_MODE=SIM Use it instead if there is no SGX hardware available.
./app 0