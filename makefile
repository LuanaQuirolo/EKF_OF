GMM_FLAGS = -std=c++11 -Wall -Werror -Wconversion -pedantic -g -Wno-deprecated-declarations -o
VALGRIND_FLAGS = --leak-check=full --track-origins=yes --show-reachable=yes
NAME = main

default: exec

all: build valgrind exec

build: 
	gcc *.c -o $(NAME)

valgrind: build
	valgrind $(VALGRIND_FLAGS) ./$(NAME) 

exec: build
	./$(NAME) 
