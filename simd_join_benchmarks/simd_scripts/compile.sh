#!/bin/bash

cd ..
make clean
make sgx
./app -a CHTfsimd
cd scripts