
ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = c++
endif

CXXFLAGS     = -std=c++11
CXXFLAGS     += `root-config --cflags`
CXXFLAGS     += `fastjet-config --cxxflags`
CXXFLAGS     += -I$(PYTHIA8)/include

LDFLAGS      = 
LDFLAGS     += `root-config --libs`
LDFLAGS     += `fastjet-config --libs --plugins`
LDFLAGS     += -L$(PYTHIA8)/lib -lpythia8

#all : PythiaHepMC

all: runBeamGasHepMC.exe

HEPMC_DIR = $(EICDIRECTORY)

#PythiaHepMC : PythiaHepMC.o
#	$(CXX) PythiaHepMC.o -o PythiaHepMC $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HE#PMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

#PythiaHepMC.o : PythiaHepMC.cxx CustomPythia8ToHepMC3.h
#	$(CXX) $(CXXFLAGS) -c PythiaHepMC.cxx -o PythiaHepMC.o -I$(HEPMC_DIR)/include -I.

#runPhoto.exe: PythiaPhotoProduction.o
#	$(CXX) PythiaPhotoProduction.o -o runPhoto.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

#PythiaPhotoProduction.o : PythiaPhotoProduction.cxx CustomPythia8ToHepMC3.h
#	$(CXX) $(CXXFLAGS) -c PythiaPhotoProduction.cxx -o PythiaPhotoProduction.o -I$(HEPMC_DIR)/include -I.

#runDIS.exe: PythiaDIS.o
#	$(CXX) PythiaDIS.o -o runDIS.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

#PythiaDIS.o : PythiaDIS.cxx CustomPythia8ToHepMC3.h
#	$(CXX) $(CXXFLAGS) -c PythiaDIS.cxx -o PythiaDIS.o -I$(HEPMC_DIR)/include -I.

runBeamGasHepMC.exe: BeamGasEvents.o eicBeamShape.o
	$(CXX) BeamGasEvents.o eicBeamShape.o -o runBeamGasHepMC.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

BeamGasEvents.o : BeamGasEvents.cxx eicBeamShape.o
	$(CXX) $(CXXFLAGS) -c BeamGasEvents.cxx -o BeamGasEvents.o -I$(HEPMC_DIR)/include -I.

eicBeamShape.o: eicBeamShape.cxx
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c eicBeamShape.cxx -o eicBeamShape.o

clean :
	rm -vf *.o *.exe *~
