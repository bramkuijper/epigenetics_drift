# makefile to compile the indivual-based simulations

xepigenetics_drift : epigenetics_drift.cpp individual.cpp
	g++ -std=c++11 -Wall -O3 -o xepigenetics_drift epigenetics_drift.cpp individual.cpp


clean:
	rm -rf xepigenetics_drift

.PHONY: clean
