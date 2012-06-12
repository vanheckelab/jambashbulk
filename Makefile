opts = -Wall -Wextra
o3opts = $(opts) -O3
dopts = $(opts) -g

all: debug o3

debug: jam2D_d
o3: jam2D
win: jam2D.exe

jam2D: jamBashbulk.cpp
	g++ $(o3opts) -o jam2D jamBashbulk.cpp

jam2D_d: jamBashbulk.cpp
	g++ $(dopts) -o jam2D_d jamBashbulk.cpp	

jam2D.exe: jamBashbulk.cpp
	i586-mingw32msvc-g++ $(o3opts) -o jam2D.exe jamBashbulk.cpp

clean:
	rm -f jam2D_d jam2D
