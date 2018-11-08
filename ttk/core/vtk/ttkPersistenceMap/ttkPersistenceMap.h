/// \ingroup vtk
/// \class ttkPersistenceMap
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date February 2018.
///
/// \brief TTK VTK-filter that wraps the persistenceMap processing package.
///
/// VTK wrapping code for the @PersistenceMap package.
/// 
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttk::PersistenceMap
#pragma once

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>

// ttk code includes
#include                  <PersistenceMap.h>
#include                  <ttkWrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkPersistenceMap
#else
class ttkPersistenceMap
#endif
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
    
    static ttkPersistenceMap* New();
    vtkTypeMacro(ttkPersistenceMap, vtkDataSetAlgorithm)
    
    // default ttk setters
    vtkSetMacro(debugLevel_, int);
    
    void SetThreadNumber(int threadNumber){\
      ThreadNumber = threadNumber;\
      SetThreads();\
    }\
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters
    
    vtkSetMacro(CriticalType, int);
    vtkGetMacro(CriticalType, int);
    
    vtkSetMacro(DefaultAmplitude, double);
    vtkGetMacro(DefaultAmplitude, double);
    
    vtkSetMacro(DefaultDeviation, double);
    vtkGetMacro(DefaultDeviation, double);
        
    vtkSetMacro(PersistenceScale, double);
    vtkGetMacro(PersistenceScale, double);
    
    vtkSetMacro(PersistenceSpread, double);
    vtkGetMacro(PersistenceSpread, double);
    
    vtkSetMacro(PersistenceThreshold, double);
    vtkGetMacro(PersistenceThreshold, double);
    
    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);

    vtkSetMacro(ReplaceScalarField, int);
    vtkGetMacro(ReplaceScalarField, int);

    vtkSetMacro(BatchMode, int);
    vtkGetMacro(BatchMode, int);

  protected:
   
    ttkPersistenceMap(){
      
        // init
      CriticalType = 2;
      PersistenceSpread = 0.1;
      PersistenceScale = 0.5;
      DefaultDeviation = 0;
      DefaultAmplitude = 0;
      PersistenceThreshold = 0.01;
      UseAllCores = false;
      ReplaceScalarField = 0;
      BatchMode = 0;
      ReplaceScalarField = 0;
      
      outputScalarField_ = NULL;
    }
    
    ~ttkPersistenceMap(){};
    
    TTK_SETUP();
    
    
  private:
    
    int                   CriticalType;
    double                PersistenceThreshold,
                          PersistenceScale, PersistenceSpread, 
                          DefaultAmplitude, DefaultDeviation;
    string                ScalarField;
    int                   ReplaceScalarField;
    int                   BatchMode;
    vtkSmartPointer<vtkDoubleArray> outputScalarField_;
    PersistenceMap        persistenceMap_;
    
};
