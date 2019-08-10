PROJ_DIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
OCCA_DIR = /Users/wuxinhui/occa
#OCCA_DIR = /home/philipwxh/occa

ifndef OCCA_DIR
  include $(PROJ_DIR)/../../../scripts/Makefile
else
  include ${OCCA_DIR}/scripts/Makefile
endif

paths += -I./include 
flags += -DOCCA_GL_ENABLED=1 -g
#flags += -DOCCA_GL_ENABLED=1 -O3

ifeq ($(OS),OSX)
	linkerFlags += -framework OpenGL -framework GLUT
endif

ifeq ($(OS),LINUX)
	linkerFlags +=
# -lGLU -lglut
endif

compilerFlags += -std=c++11 -Wno-deprecated -Wno-unused-result


#---[ COMPILATION ]-------------------------------
headers = $(wildcard $(incPath)/*.hpp) $(wildcard $(incPath)/*.tpp)
sources = $(wildcard $(srcPath)/*.cpp)

objects  = $(subst $(srcPath)/,$(objPath)/,$(sources:.cpp=.o))

.PHONY: clean

all: main2d main2d_baseline runner2d

main2d_baseline: $(objects) $(headers) main2d_baseline.cpp
	$(compiler) $(compilerFlags) -o main2d_baseline $(flags) $(objects) main2d_baseline.cpp -L${OCCA_DIR}/lib $(paths) $(linkerFlags)

main2d: $(objects) $(headers) main2d.cpp
	$(compiler) $(compilerFlags) -o main2d $(flags) $(objects) main2d.cpp -L${OCCA_DIR}/lib $(paths) $(linkerFlags)

runner2d: runner2d.cpp
	g++ $(compilerFlags) runner2d.cpp -o runner2d;

$(objPath)/%.o:$(srcPath)/%.cpp $(wildcard $(subst $(srcPath)/,$(incPath)/,$(<:.cpp=.hpp))) $(wildcard $(subst $(srcPath)/,$(incPath)/,$(<:.cpp=.tpp)))
	$(compiler) $(compilerFlags) -o $@ $(flags) -c $(paths) $<

clean:
	rm -f $(objPath)/*.o;
	rm -f runner2d;
	rm -f runner3d;
	rm -f main;
	rm -f main2d;
	rm -f main2d_baseline;
	rm -f main3d;
	rm -f main3d_baseline;
	rm -rf main2d.dSYM/;
	rm -rf main2d_baseline.dSYM;
	rm -rf main3d.dSYM/;
	rm -rf main3d_baseline.dSYM;

#=================================================
