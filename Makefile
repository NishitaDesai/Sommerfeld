

somm.x: SommCorr.cc
	g++ -O3 SommCorr.cc -o somm.x 	

toy.x: toy.cpp
	g++ -O3 toy.cpp -o toy.x 

.phony:
clean:
	rm somm.x
	rm *.txt
