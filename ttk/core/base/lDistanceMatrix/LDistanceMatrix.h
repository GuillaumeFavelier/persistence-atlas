/// \ingroup base
/// \class ttk::LDistanceMatrix 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %lDistanceMatrix processing package.
///
/// %LDistanceMatrix is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkLDistanceMatrix.cpp %for a usage example.

#pragma once

// base code includes
#include<LDistance.h>
#include<Wrapper.h>
#include<numeric>

namespace ttk{

  class LDistanceMatrix : public Debug{

    public:

      LDistanceMatrix();

      ~LDistanceMatrix();

      template <class dataType>
        int execute() const;

      inline int setDistanceType(const string distanceType){
        distanceType_=distanceType;
        return 0;
      }

      inline int setTransposeMatrix(const bool transposeMatrix){
        transposeMatrix_=transposeMatrix;
        return 0;
      }

      inline int setNumberOfPoints(const int numberOfPoints){
        numberOfPoints_=numberOfPoints;
        return 0;
      }

      inline int setNumberOfScalarFields(const int numberOfScalarFields){
        numberOfScalarFields_=numberOfScalarFields;
        return 0;
      }

      inline int setInputScalarFields(vector<void*>* inputScalarFields){
        inputScalarFields_=inputScalarFields;
        return 0;
      }

      inline int setOutputDistanceMatrix(void *data){
        outputDistanceMatrix_=data;
        return 0;
      }

    protected:

      string distanceType_;
      bool transposeMatrix_;
      int numberOfPoints_;
      int numberOfScalarFields_;
      vector<void*>* inputScalarFields_;
      void* outputDistanceMatrix_;
  };
}

template <class dataType>
int LDistanceMatrix::execute() const{
  Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPoints_<=0) return -1;
  if(!inputScalarFields_) return -1;
  if(!outputDistanceMatrix_) return -1;
#endif

  double* outputDistanceMatrix=static_cast<double*>(outputDistanceMatrix_);

  // compute matrix size
  int N=0;
  if(transposeMatrix_)
    N=numberOfPoints_;
  else
    N=numberOfScalarFields_;
#ifndef TTK_ENABLE_KAMIKAZE
  if(N<=1) return -1;
#endif

  // get scalar arrays
  vector<dataType*> inputScalarFields(N);
  for(int i=0; i<numberOfScalarFields_; ++i){
    inputScalarFields[i]=static_cast<dataType*>((*inputScalarFields_)[i]);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalarFields[i]) return -1;
#endif
  }

  const int n=atoi(distanceType_.c_str());

  const int matrixSize=N*N;

  std::fill(outputDistanceMatrix, outputDistanceMatrix+matrixSize, 0);

  const int K=N*(N-1)/2;
#ifdef TTK_ENABLE_OPENMP
# pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int k=0; k<K; ++k){
    // Linear index upper triangular matrix
    const int i = N - 2 - floor(sqrt(-8*k + 4*N*(N-1)-7)/2.0 - 0.5);
    const int j = k + i + 1 - N*(N-1)/2 + (N-i)*((N-i)-1)/2;

    LDistance dist;
    dist.setDebugLevel(fatalMsg);
    dist.setThreadNumber(1);

    if(transposeMatrix_){
      if(n>0)
        dist.computeLn<dataType>(inputScalarFields, i, j, nullptr, n, numberOfScalarFields_);
      else
        dist.computeLinf<dataType>(inputScalarFields, i, j, nullptr, numberOfScalarFields_);
    }
    else{
      dataType* scalarField1=inputScalarFields[i];
      dataType* scalarField2=inputScalarFields[j];

      if(n>0)
        dist.computeLn<dataType>(scalarField1, scalarField2, nullptr, n, numberOfPoints_);
      else
        dist.computeLinf<dataType>(scalarField1, scalarField2, nullptr, numberOfPoints_);
    }

    const dataType result=dist.getResult();
    outputDistanceMatrix[i*N+j]=result;
    outputDistanceMatrix[j*N+i]=result;
  }

  {
    stringstream msg;
    msg << "[LDistanceMatrix] Data-set (" << numberOfPoints_
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}
