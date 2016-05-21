#PLAT=mic
PLAT=cpu
APPNAME=SmithWatermanAffineGapMod
#APPNAME=SmithWatermanMod
#APPNAME=NeedlemanWunschAffineGapMod
#APPNAME=NeedlemanWunschMod

.PHONY: all
all:
ifeq ($(PLAT),cpu)
	icpc -O3 -std=c++11 -ip -xCORE-AVX2 -I../include -Iinclude -o ../build/Config_CPU.o -c ../src/Config_CPU.cpp
	icpc -O3 -std=c++11 -ip -xCORE-AVX2 -I../include -Iinclude -o ../build/$(APPNAME).o -c $(APPNAME).cpp
	icpc -O3 -std=c++11 -ip -xCORE-AVX2 -I../include -Iinclude -I../Utils -o $(APPNAME).out $(APPNAME)_main.cpp ../build/*.o
else
	icpc -O3 -std=c++11 -ip -mmic -I../include -Iinclude -o ../build/Config.o -c ../src/Config.cpp
	icpc -O3 -std=c++11 -ip -mmic -I../include -Iinclude -o ../build/$(APPNAME).o -c $(APPNAME).cpp
	icpc -O3 -std=c++11 -ip -mmic -I../include -Iinclude -I../Utils -o $(APPNAME).out $(APPNAME)_main.cpp ../build/*.o
endif

.PHONY: clean
clean:
	rm ../build/*.o *.out



