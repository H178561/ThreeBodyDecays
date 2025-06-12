# ThreeBodyDecays Project

Implementation of the [ThreeBodyDecays.jl](https://github.com/mmikhasenko/ThreeBodyDecays.jl) Julia class in C++ for usage in EvtGen.

For the usage in EvtGen with json files like in the [amplitude serialisation project](https://github.com/RUB-EP1/amplitude-serialization/tree/main) see [EvtThreeBodyDecay](https://github.com/H178561/EvtThreeBodyDecay)


## Build

To build the Model 

```
mkdir build
cd build
cmake ..
make
```

The file contains various G-Tests which can be tested with

```
./ThreeBodyDecaysTest
./ClebschGordanTest
./JsonModelTest
```

There is also an example of the usage in 

```
Examples/AmplitudeModelExample
```


## How to use in EvtGen

For the usage in EvtGen the Project has to be called in the EvtGen setup

in `CmakeList.txt`
