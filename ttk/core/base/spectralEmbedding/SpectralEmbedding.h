/// \ingroup base
/// \class ttk::SpectralEmbedding 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %spectralEmbedding processing package.
///
/// %SpectralEmbedding is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkSpectralEmbedding.cpp %for a usage example.
#pragma once

// base code includes
#include<Wrapper.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include<Python.h>
#include<numpy/arrayobject.h>

namespace ttk{

  class SpectralEmbedding : public Debug{

    public:

      SpectralEmbedding();
      ~SpectralEmbedding();

      inline int setInputModulePath(const string& modulePath){
        modulePath_=modulePath;
        return 0;
      }

      inline int setInputModuleName(const string& moduleName){
        moduleName_=moduleName;
        return 0;
      }

      inline int setInputFunctionName(const string& functionName){
        functionName_=functionName;
        return 0;
      }

      inline int setInputMatrixDimension(const int matrixDimension){
        matrixDimension_=matrixDimension;
        return 0;
      }

      inline int setInputDistanceMatrix(void* data){
        distanceMatrix_=data;
        return 0;
      }

      inline int setInputEmbeddingDimension(int embeddingDimension){
        embeddingDimension_=embeddingDimension;
        return 0;
      }

      inline int setInputMinimumNumberOfComponents(int minimumNumberOfComponents){
        minimumNumberOfComponents_=minimumNumberOfComponents;
        return 0;
      }

      inline int setInputMaximumNumberOfComponents(int maximumNumberOfComponents){
        maximumNumberOfComponents_=maximumNumberOfComponents;
        return 0;
      }

      inline int setInputNumberOfNeighbors(int numberOfNeighbors){
        numberOfNeighbors_=numberOfNeighbors;
        return 0;
      }

      inline int setInputSigma(double sigma){
        sigma_=sigma;
        return 0;
      }

      inline int setOutputEmbedding(void* data){
        embedding_=data;
        return 0;
      }

      inline int setOutputEigenvalues(void* data){
        eigenvalues_=data;
        return 0;
      }

      inline int setOutputComponents(vector<void*>* data){
        components_=data;
        return 0;
      }

      int execute() const;

    protected:

      string modulePath_;
      string moduleName_;
      string functionName_;
      int matrixDimension_;
      int embeddingDimension_;
      int minimumNumberOfComponents_;
      int maximumNumberOfComponents_;
      int numberOfNeighbors_;
      double sigma_;
      void* distanceMatrix_;
      void* embedding_;
      void* eigenvalues_;
      vector<void*>* components_;
      char majorVersion_;
  };
}

