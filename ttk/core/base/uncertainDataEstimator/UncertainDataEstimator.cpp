#include                  <UncertainDataEstimator.h>

UncertainDataEstimator::UncertainDataEstimator(){
  vertexNumber_ = 0;
  numberOfInputFields_ = 0;
  computeLowerBound_ = true;
  computeUpperBound_ = true;
  binCount_ = 0;
  inputData_ = nullptr;
  numberOfClusters_ = 0;
  clusteringData_ = nullptr;
  outputLowerBoundFields_ = nullptr;
  outputUpperBoundFields_ = nullptr;
  outputMeanFields_ = nullptr;
  outputProbabilityFields_ = nullptr;
}

UncertainDataEstimator::~UncertainDataEstimator(){
}

