all: toymc

clean:
	rm toymc

toymc: toyMC_backup.cc
	g++ -o run toyMC.cc `root-config --libs --cflags`

