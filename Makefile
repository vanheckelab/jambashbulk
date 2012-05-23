all: debug opt

debug: jam2D_d

opt: jam2D

jam2D: jamBashbulk.cpp
	g++ -O3 -o jam2D jamBashbulk.cpp

jam2D_d: jamBashbulk.cpp
	g++ -g -o jam2D_d jamBashbulk.cpp	
