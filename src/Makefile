.PHONY: all clean test s21_math.a check gcov_report
CC=gcc
CFLAGS=-c
LDFLAGS= -Wall -Wextra -Werror
OS := $(shell uname -s)

default: all

all: s21_math.a test  

test: s21_math.a
	$(CC) $(LDFLAGS) test.c s21_math.a -o test  `pkg-config --cflags --libs check`
	./test

s21_math.a: s21_math.o
	ar rcs s21_math.a s21_math.o
	ranlib $@
	cp $@ lib$@
	ranlib lib$@

s21_math.o: s21_math.c
	$(CC) $(LDFLAGS) $(CFLAGS) -o s21_math.o s21_math.c

clean:
	rm -rf test *.a *.o *.out *.cfg fizz *.gc* *.info report

rebuild: clean all

check:
	cppcheck *.h *.c
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extension=c *.c *.h
	CK_FORK=no leaks --atExit -- ./test

gcov_report: 
	$(CC) --coverage $(LDFLAGS) test.c s21_math.c -o test `pkg-config --cflags --libs check`
	./test
	lcov -t "test" -o test.info -c -d ./
	genhtml -o report test.info
	open report/index.html
