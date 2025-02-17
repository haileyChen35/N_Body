# create executable file named program
program: simulation.cpp 
	g++ -std=c++11 -o program simulation.cpp

# plot multiple PDF
plot: simulation.cpp plot.py 
	python3 plot.py solar.tsv solar.pdf 10000

# combine multiple PDF into a GIF animation
animation: animation.py
	python3 animation.py

# clean
clean:
	rm *.pdf *.tsv *.gif program