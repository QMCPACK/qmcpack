#Utilities classes
$(ODIR)/OhmmsObject.o: $(OHMMSDIR)/src/Utilities/OhmmsObject.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/OhmmsObject.cpp -o $(ODIR)/OhmmsObject.o
$(ODIR)/RandomGenerator.o: $(OHMMSDIR)/src/Utilities/RandomGenerator.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/RandomGenerator.cpp -o $(ODIR)/RandomGenerator.o
$(ODIR)/SimpleParser.o: $(OHMMSDIR)/src/Utilities/SimpleParser.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/SimpleParser.cpp -o $(ODIR)/SimpleParser.o
$(ODIR)/RandRandom.o: $(OHMMSDIR)/src/Utilities/RandRandom.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/RandRandom.cpp -o $(ODIR)/RandRandom.o
$(ODIR)/OhmmsSpecies.o: $(OHMMSDIR)/src/Utilities/OhmmsSpecies.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/OhmmsSpecies.cpp -o $(ODIR)/OhmmsSpecies.o
$(ODIR)/SpeciesCollection.o: $(OHMMSDIR)/src/Utilities/SpeciesCollection.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/SpeciesCollection.cpp -o $(ODIR)/SpeciesCollection.o
$(ODIR)/MDParameters.o: $(OHMMSDIR)/src/Utilities/MDParameters.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/MDParameters.cpp -o $(ODIR)/MDParameters.o
$(ODIR)/OhmmsInform.o: $(OHMMSDIR)/src/Utilities/OhmmsInform.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/OhmmsInform.cpp -o $(ODIR)/OhmmsInform.o
$(ODIR)/OhmmsInfo.o: $(OHMMSDIR)/src/Utilities/OhmmsInfo.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/OhmmsInfo.cpp -o $(ODIR)/OhmmsInfo.o	
$(ODIR)/GlobalMDParameters.o: $(OHMMSDIR)/src/Utilities/GlobalMDParameters.cpp
	$(CXX) $(CXXFLAGSAPP) $(OHMMSINC) -c $(OHMMSDIR)/src/Utilities/GlobalMDParameters.cpp -o $(ODIR)/GlobalMDParameters.o
# DO NOT DELETE

OhmmsObject.o: $(OHMMSDIR)/src/Utilities/OhmmsObject.h

RandRandom.o: $(OHMMSDIR)/src/Utilities/RandRandom.h

RandomGenerator.o: $(OHMMSDIR)/src/Utilities/RandomGenerator.h
RandomGenerator.o: $(OHMMSDIR)/src/Utilities/RandRandom.h

OhmmsSpecies.o: $(OHMMSDIR)/src/Utilities/OhmmsSpecies.h

SimpleParser.o: $(OHMMSDIR)/src/Utilities/SimpleParser.h

SpeciesCollection.o: $(OHMMSDIR)/src/Utilities/SpeciesCollection.h
SpeciesCollection.o: $(OHMMSDIR)/src/Utilities/OhmmsSpecies.h

MDParameters.o: $(OHMMSDIR)/src/Utilities/MDParameters.h
MDParameters.o: $(OHMMSDIR)/src/Ohmms.h
MDParameters.o: $(OHMMSDIR)/src/OhmmsUnits.h
MDParameters.o: $(OHMMSDIR)/src/Utilities/OhmmsInfo.h
MDParameters.o: $(OHMMSDIR)/src/Utilities/OhmmsInform.h
MDParameters.o: $(OHMMSDIR)/src/Utilities/RandomGenerator.h
MDParameters.o: $(OHMMSDIR)/src/Utilities/RandRandom.h
MDParameters.o: $(OHMMSDIR)/src/Utilities/SimpleParser.h

OhmmsInform.o: $(OHMMSDIR)/src/Utilities/OhmmsInform.h

OhmmsInfo.o: $(OHMMSDIR)/src/Utilities/OhmmsInfo.h
OhmmsInfo.o: $(OHMMSDIR)/src/Utilities/OhmmsInform.h

GlobalMDParameters.o: $(OHMMSDIR)/src/Utilities/GlobalMDParameters.h
GlobalMDParameters.o: $(OHMMSDIR)/src/Utilities/MDParameters.h
GlobalMDParameters.o: $(OHMMSDIR)/src/Ohmms.h
GlobalMDParameters.o: $(OHMMSDIR)/src/Utilities/OhmmsInfo.h
GlobalMDParameters.o: $(OHMMSDIR)/src/Utilities/OhmmsInform.h
GlobalMDParameters.o: $(OHMMSDIR)/src/Utilities/RandomGenerator.h
GlobalMDParameters.o: $(OHMMSDIR)/src/Utilities/RandRandom.h
GlobalMDParameters.o: $(OHMMSDIR)/src/OhmmsUnits.h
