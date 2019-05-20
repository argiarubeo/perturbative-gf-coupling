CC	= g++ -O2
LD  = -larmadillo

all: Nmag Nelec


Nmag:
	$(CC) ../src/Nmag.cpp -o Nmag $(LD)

Nelec:
	$(CC) ../src/Nelec.cpp -o Nelec $(LD)

clean:
	rm -f ./*.o
	rm -f ./Nmag
	rm -f ./Nelec
