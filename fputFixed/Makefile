all: standards couplings state modeAmplitude tangent newton agp saveFile fputBreather fputAGPGradient fputAGPFlow fputAGPFlow2 fputAdiabatic fputAGP fputAGPNorm fputAGPNorm2A fputBreatherAmplitude fputAGPNormAvg fputEvolveRev fputAGPNormAvgSeed

standards: src/standards.cpp include/standards.hpp include/eigenClasses.hpp
	g++ -c src/standards.cpp -o lib/standards.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

couplings: src/couplings.cpp include/couplings.hpp include/constants.hpp include/eigenClasses.hpp include/standards.hpp
	g++ -c src/couplings.cpp -o lib/couplings.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

state: src/state.cpp include/state.hpp include/eigenClasses.hpp include/couplings.hpp include/standards.hpp include/modeAmplitude.hpp
	g++ -c src/state.cpp -o lib/state.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

modeAmplitude: src/modeAmplitude.cpp include/modeAmplitude.hpp include/eigenClasses.hpp include/standards.hpp include/couplings.hpp
	g++ -c src/modeAmplitude.cpp -o lib/modeAmplitude.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

tangent: src/tangent.cpp include/tangent.hpp include/eigenClasses.hpp include/standards.hpp include/couplings.hpp include/state.hpp
	g++ -c src/tangent.cpp -o lib/tangent.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

newton: src/newton.cpp include/newton.hpp include/eigenClasses.hpp include/standards.hpp include/couplings.hpp include/state.hpp include/tangent.hpp include/modeAmplitude.hpp
	g++ -c src/newton.cpp -o lib/newton.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

agp: src/agp.cpp include/agp.hpp include/eigenClasses.hpp include/standards.hpp include/state.hpp
	g++ -c src/agp.cpp -o lib/agp.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

saveFile: src/saveFile.cpp include/saveFile.hpp include/eigenClasses.hpp include/standards.hpp include/state.hpp include/agp.hpp
	g++ -c src/saveFile.cpp -o lib/saveFile.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include

fputBreather: src/fputBreather.cpp include/fputFixed.hpp
	g++ -c src/fputBreather.cpp -o lib/fputBreather.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputBreather.out lib/fputBreather.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAGPGradient: src/fputAGPGradient.cpp include/fputFixed.hpp
	g++ -c src/fputAGPGradient.cpp -o lib/fputAGPGradient.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGPGradient.out -O3 lib/fputAGPGradient.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAGPFlow: src/fputAGPFlow.cpp include/fputFixed.hpp
	g++ -c src/fputAGPFlow.cpp -o lib/fputAGPFlow.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGPFlow.out -O3 lib/fputAGPFlow.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAGPFlow2: src/fputAGPFlow2.cpp include/fputFixed.hpp
	g++ -c src/fputAGPFlow2.cpp -o lib/fputAGPFlow2.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGPFlow2.out -O3 lib/fputAGPFlow2.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAdiabatic: src/fputAdiabatic.cpp include/fputFixed.hpp
	g++ -c src/fputAdiabatic.cpp -o lib/fputAdiabatic.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAdiabatic.out -O3 lib/fputAdiabatic.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAGP: src/fputAGP.cpp include/fputFixed.hpp
	g++ -c src/fputAGP.cpp -o lib/fputAGP.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGP.out -O3 lib/fputAGP.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAGPNorm: src/fputAGPNorm.cpp include/fputFixed.hpp
	g++ -c src/fputAGPNorm.cpp -o lib/fputAGPNorm.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGPNorm.out -O3 lib/fputAGPNorm.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAGPNormAvg: src/fputAGPNormAvg.cpp include/fputFixed.hpp
	g++ -c src/fputAGPNormAvg.cpp -o lib/fputAGPNormAvg.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGPNormAvg.out -O3 lib/fputAGPNormAvg.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

fputAGPNorm2A: src/fputAGPNorm2A.cpp include/fputFixed.hpp
	g++ -c src/fputAGPNorm2A.cpp -o lib/fputAGPNorm2A.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGPNorm2A.out -O3 lib/fputAGPNorm2A.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o


fputBreatherAmplitude: src/fputBreatherAmplitude.cpp include/fputFixed.hpp
	g++ -c src/fputBreatherAmplitude.cpp -o lib/fputBreatherAmplitude.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputBreatherAmplitude.out -O3 lib/fputBreatherAmplitude.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o


fputEvolveRev: src/fputEvolveRev.cpp include/fputFixed.hpp
	g++ -c src/fputEvolveRev.cpp -o lib/fputEvolveRev.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputEvolveRev.out -O3 lib/fputEvolveRev.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o
 
fputAGPNormAvgSeed: src/fputAGPNormAvgSeed.cpp include/fputFixed.hpp
	g++ -c src/fputAGPNormAvgSeed.cpp -o lib/fputAGPNormAvgSeed.o -O3 -I/usr3/graduate/nachiket/Libraries/json/include -I/usr3/graduate/nachiket/Libraries/eigen -I./include
	g++ -o _fputAGPNormAvgSeed.out -O3 lib/fputAGPNormAvgSeed.o lib/saveFile.o lib/agp.o lib/newton.o lib/tangent.o lib/modeAmplitude.o lib/state.o lib/couplings.o lib/standards.o

