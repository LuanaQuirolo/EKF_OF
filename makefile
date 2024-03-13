GMM_FLAGS = -std=c++11 -Wall -Werror -lprofiler -Wconversion -pedantic -g -Wno-deprecated-declarations -o -fsingle-precision-constant -Wdouble-promotion -Wall -O3 -ffast-math -fsingle-precision-constant
VALGRIND_FLAGS = --leak-check=full --track-origins=yes --show-reachable=yes
NAME = main

default: exec

all: build valgrind exec

build: 
	gcc EKF.c matrix.c quaternions.c main.c -o  $(NAME) -fprofile-arcs -ftest-coverage -fsingle-precision-constant -Wdouble-promotion -Wall -O3 -ffast-math -fsingle-precision-constant

valgrind: build
	valgrind $(VALGRIND_FLAGS) ./$(NAME) 

exec: build
	./$(NAME) 

exec-all: build
#	./$(NAME) "quieto"
#	./$(NAME) "lineal"
#	./$(NAME) "rect"
#	./$(NAME) "rectangular"
#	./$(NAME) "rot_roll_30"
#	./$(NAME) "rot_pitch_30"
#	./$(NAME) "rot_yaw_30"
#	./$(NAME) "inf"
#	./$(NAME) "infinito"
#	./$(NAME) "altura"
#	./$(NAME) "u"
	./$(NAME) "circulo"