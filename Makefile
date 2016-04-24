all:
	g++  main.cpp  -o fluidx -framework GLUT -framework OpenGL -Wno-deprecated -O3

clean:
	rm -rf fluidx
