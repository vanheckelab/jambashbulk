src = src
bin = bin

opts = 
warnopts = -Wall -Wextra -Wconversion -Wno-sign-conversion
o3opts = $(opts) -O2
dopts = $(warnopts) $(opts) -g

srcfiles = $(src)/jamBashbulk.cpp
headers = $(src)/fheader.h
gitversion = $(src)/gitversion.h
allheaders = $(headers) $(gitversion)

binary = $(bin)/jam2D
binary_d = $(binary)_d
dllname = jamBashbulk.so
dll = $(bin)/$(dllname)

all: $(binary) $(binary_d) $(dll)
dll: $(dll)

clean: | $(bin)
	rm -rf $(bin)
	rm src/gitversion.h

# normal binaries (-O3 and debug)

$(binary): $(srcfiles) $(allheaders) $(binary_d) | $(bin)
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

binary_win32 = $(binary).exe
dllname_win32 = jamBashbulk.dll
dll_win32 = $(bin)/$(dllname_win32)
win32: $(binary_win32) $(dll_win32)

$(binary_win32): $(srcfiles) $(allheaders) | $(bin)
	i586-mingw32msvc-g++ $(o3opts) -o $(bin)/jam2D.exe $(srcfiles)

$(dll_win32): $(srcfiles) $(allheaders) | $(bin)
	i586-mingw32msvc-g++ -shared -Wl,-soname,$(dllname_win32) $(o3pts) -o $(dll_win32) $(srcfiles)
	

# dll to work from python
$(dll): $(srcfiles) $(allheaders) | $(bin)
	g++ -fPIC -shared -Wl,-soname,$(dllname) $(o3pts) -o $(dll) $(srcfiles)
