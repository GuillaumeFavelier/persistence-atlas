#include<LDistanceMatrix.h>

LDistanceMatrix::LDistanceMatrix():
  numberOfPoints_{0},
  numberOfScalarFields_{0},
  inputScalarFields_{nullptr},
  outputDistanceMatrix_{nullptr}
{}

LDistanceMatrix::~LDistanceMatrix(){
}

