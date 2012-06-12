src = src
bin = bin

opts = -Wall -Wextra
o3opts = $(opts) -O3
dopts = $(opts) -g

srcfiles = $(src)/jamBashbulk.cpp

binary = $(bin)/jam2D
binary_d = $(binary)_d

all: $(binary) $(binary_d)

clean: $(bin)
	rm -rf $(bin)

# normal binaries (-O3 and debug)

$(binary):  $(bin) $(srcfiles)
	g++ $(o3opts) -o $(bin)/jam2D $(srcfiles)

$(binary_d):  $(bin) $(srcfiles)
	g++ $(dopts) -o $(bin)/jam2D_d $(srcfiles)

$(bin):
	mkdir -p $(bin)


# specials

# win32

binary_exe = $(binary).exe
win32: $(binary_exe)

$(binary_exe): $(bin) $(srcfiles)
	i586-mingw32msvc-g++ $(o3opts) -o $(bin)/jam2D.exe $(srcfiles)

# dll to work from python
dll: $(bin) $(src)/j2d_dll.cpp $(srcfiles)
	g++ -shared -Wl,-soname,j2d_dll.so -o $(bin)/j2d_dll.so $(src)/j2d_dll.cpp