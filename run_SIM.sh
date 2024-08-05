#!/bin/bash

nohup ./zoupa < input_restart_testmuphi1 > output_muphi_test1 &

cd ./python 

python muphi_plot.py
