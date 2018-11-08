/// \ingroup vtk
/// \class ttkSpectralEmbedding
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the spectralEmbedding processing package.
///
/// VTK wrapping code for the @SpectralEmbedding package.
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
/// \sa ttk::SpectralEmbedding
#pragma once

// ttk code includes
#include<SpectralEmbedding.h>
#include<ttkWrapper.h>

// VTK includes
#include<vtkLongLongArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>
#include<vtkTable.h>
#include<vtkPoints.h>


class VTKFILTERSCORE_EXPORT ttkSpectralEmbedding 
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkSpectralEmbedding* New();
    vtkTypeMacro(ttkSpectralEmbedding, vtkDataSetAlgorithm)

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

    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);

    vtkSetMacro(ModulePath, string);
    vtkGetMacro(ModulePath, string);

    vtkSetMacro(ModuleName, string);
    vtkGetMacro(ModuleName, string);

    vtkSetMacro(FunctionName, string);
    vtkGetMacro(FunctionName, string);

    vtkSetMacro(EmbeddingDimension, int);
    vtkGetMacro(EmbeddingDimension, int);

    vtkSetMacro(MinimumNumberOfComponents, int);
    vtkGetMacro(MinimumNumberOfComponents, int);

    vtkSetMacro(MaximumNumberOfComponents, int);
    vtkGetMacro(MaximumNumberOfComponents, int);

    vtkSetMacro(NumberOfNeighbors, int);
    vtkGetMacro(NumberOfNeighbors, int);

    vtkSetMacro(Sigma, double);
    vtkGetMacro(Sigma, double);

    vtkSetMacro(MetaData, string);
    vtkGetMacro(MetaData, string);

    int FillInputPortInformation(int port, vtkInformation *info){

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
          break;
        default:
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info){

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
          break;
        default:
          break;
      }

      return 1;
    }

  protected:

    ttkSpectralEmbedding(){
      UseAllCores=false;
      ThreadNumber=1;
      EmbeddingDimension=2;
      MinimumNumberOfComponents=1;
      MaximumNumberOfComponents=2;
      NumberOfNeighbors=5;
      Sigma=1.0;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkSpectralEmbedding(){};

    TTK_SETUP();

  private:

    string ScalarField;
    string ModulePath;
    string ModuleName;
    string FunctionName;
    int EmbeddingDimension;
    int MinimumNumberOfComponents;
    int MaximumNumberOfComponents;
    int NumberOfNeighbors;
    double Sigma;
    string MetaData;
    SpectralEmbedding spectralEmbedding_;
};
