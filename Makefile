#ODEPath=/home/pzpzpzp1/Desktop/ode-0.13

all:
	g++ -DHAVE_CONFIG_H -I$(ODEPath)/ode/demo -I$(ODEPath)/ode/src  -I$(ODEPath)/include -I$(ODEPath)/include -DDRAWSTUFF_TEXTURE_PATH="\"$(ODEPath)/drawstuff/textures\"" -DdTRIMESH_ENABLED   -g -O2 -MT nest.o -MD -MP -MF ./nest.Tpo -c -o nest.o nest.cpp
	mv -f ./nest.Tpo ./nest.Po
	/bin/bash $(ODEPath)/libtool  --tag=CXX   --mode=link g++  -g -O2     -o nest nest.o $(ODEPath)/drawstuff/src/libdrawstuff.la $(ODEPath)/ode/src/libode.la -framework OpenGL -framework GLUT  -lm  -lpthread
#	$(ODEPath)/libccd/libtool: link: g++ -g -O2 -o nest nest.o  $(ODEPath)/drawstuff/src/.libs/libdrawstuff.a -lX11 $(ODEPath)/ode/src/.libs/libode.a -lGLU -lGL -lm -lpthread

start: ./nest -notex -noshadow

clean:
	rm nest