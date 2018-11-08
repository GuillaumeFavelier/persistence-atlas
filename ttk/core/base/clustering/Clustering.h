/// \ingroup base
/// \class ttk::Clustering 
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %clustering processing package.
///
/// %Clustering is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa vtkClustering.cpp %for a usage example.
#pragma once

// base code includes
#include<Wrapper.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include<Python.h>
#include<numpy/arrayobject.h>

namespace ttk{

  class Clustering : public Debug{

    public:

      Clustering();
      ~Clustering();

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

      inline int setInputNumberOfObjects(const int numberOfObjects){
        numberOfObjects_=numberOfObjects;
        return 0;
      }

      inline int setInputMaximumNumberOfClusters(const int maximumNumberOfClusters){
        maximumNumberOfClusters_=maximumNumberOfClusters;
        return 0;
      }

      inline int setInputMinimumNumberOfClusters(const int minimumNumberOfClusters){
        minimumNumberOfClusters_=minimumNumberOfClusters;
        return 0;
      }

      inline int setInputCurrentGap(const int currentGap){
        currentGap_=currentGap;
        return 0;
      }

      inline int setInputNumberOfComponents(const int numberOfComponents){
        numberOfComponents_=numberOfComponents;
        return 0;
      }

      inline int setInputMaps(vector<double>* data){
        maps_=data;
        return 0;
      }

      inline int setInputEigenvalues(void* data){
        eigenvalues_=data;
        return 0;
      }

      inline int setOutputAssignation(void* data){
        assignation_=data;
        return 0;
      }

      int execute() const;

    protected:

      string modulePath_;
      string moduleName_;
      string functionName_;
      int numberOfObjects_;
      int maximumNumberOfClusters_;
      int minimumNumberOfClusters_;
      int currentGap_;
      int numberOfComponents_;
      vector<double>* maps_;
      void* eigenvalues_;
      void* assignation_;
      char majorVersion_;
  };
}

