src = src
bin = bin

opts = -Wall -Wextra -Wconversion
o3opts = $(opts) -O3
dopts = $(opts) -g

srcfiles = $(src)/jamBashbulk.cpp
headers = $(src)/fheader.h
gitversion = $(src)/gitversion.h
allheaders = $(headers) $(gitversion)

binary = $(bin)/jam2D
binary_d = $(binary)_d

all: $(binary) $(binary_d)

clean: | $(bin)
	rm -rf $(bin)
	rm src/gitversion.h

# normal binaries (-O3 and debug)

$(binary): $(srcfiles) $(allheaders) | $(bin)
	g++ $(o3opts) -o $(bin)/jam2D $(srcfiles)

$(binary_d):  $(srcfiles) $(allheaders) | $(bin)
	g++ $(dopts) -o $(bin)/jam2D_d $(srcfiles)

$(bin):
	mkdir -p $(bin)

$(src)/gitversion.h: .git $(srcfiles) $(headers)
	echo "#define GIT_HEAD \"$(shell git rev-parse HEAD)\"" > $@
	echo "#define GIT_HEAD_DATE \"$(shell git log -1 --pretty=format:%ci HEAD $)\"" >> $@
	echo "#define GIT_CHANGED_FILES $(shell git status --porcelain src | grep '^ M' | wc -l)" >> $@

# specials

# win32

binary_exe = $(binary).exe
win32: $(binary_exe)

$(binary_exe): $(srcfiles) $(allheaders) | $(bin)
	i586-mingw32msvc-g++ $(o3opts) -o $(bin)/jam2D.exe $(srcfiles)

# dll to work from python
dll: $(src)/j2d_dll.cpp $(srcfiles) $(allheaders) | $(bin)
	g++ -shared -Wl,-soname,j2d_dll.so -o $(bin)/j2d_dll.so $(src)/j2d_dll.cpp
