#include<ttkLDistanceMatrix.h>

vtkStandardNewMacro(ttkLDistanceMatrix)

  int ttkLDistanceMatrix::doIt(vector<vtkDataSet*>& inputs, vector<vtkDataSet*>& outputs){
    Memory m;

    vtkDataSet* input=inputs[0];
    vtkDataSet* output=outputs[0];

    vtkSmartPointer<vtkStringArray> metadata=vtkSmartPointer<vtkStringArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!metadata){
      cerr << "[ttkLDistanceMatrix] Error: vtkStringArray memory allocation problem." << endl;
      return -1;
    }
#endif

    // representative gives shared information
    vtkDataArray* representative=nullptr;
    int numberOfScalarFields=0;
    vector<void*> inputScalarFields;
    {
      vtkPointData* pointData=input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!pointData){
        cerr << "[ttkLDistanceMatrix] Error: input has no point data." << endl;
        return -1;
      }
#endif

      numberOfScalarFields=pointData->GetNumberOfArrays();
      metadata->SetNumberOfComponents(1);
      metadata->SetNumberOfTuples(numberOfScalarFields);
      metadata->SetName("ScalarField");

      if(numberOfScalarFields>1){
        inputScalarFields.resize(numberOfScalarFields, nullptr);
        for(int i=0; i<numberOfScalarFields; ++i){
          vtkDataArray* scalarField=pointData->GetArray(i);
#ifndef TTK_ENABLE_KAMIKAZE
          if(!scalarField){
            cerr << "[ttkLDistanceMatrix] Error: input scalar field " << i << " is not valid." << endl;
            return -1;
          }
#endif
          metadata->SetValue(i, scalarField->GetName());

          if(scalarField){
            inputScalarFields[i]=scalarField->GetVoidPointer(0);
            representative=scalarField;
          }
        }
      }
    }

    if(representative){
      const int numberOfPoints=representative->GetNumberOfTuples();
#ifndef TTK_ENABLE_KAMIKAZE
      if(numberOfPoints<=0){
        cerr << "[ttkLDistanceMatrix] Error: input has no point." << endl;
        return -1;
      }
#endif

      int matrixSize=0;
      if(TransposeMatrix)
        matrixSize=numberOfPoints*numberOfPoints;
      else
        matrixSize=numberOfScalarFields*numberOfScalarFields;

      vtkSmartPointer<vtkDoubleArray> distanceMatrix=vtkSmartPointer<vtkDoubleArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!distanceMatrix){
        cerr << "[ttkLDistanceMatrix] Error: vtkDoubleArray memory allocation problem." << endl;
        return -1;
      }
#endif
      distanceMatrix->SetNumberOfComponents(1);
      distanceMatrix->SetNumberOfTuples(matrixSize);
      distanceMatrix->SetName("DistanceMatrix");

      lDistanceMatrix_.setWrapper(this);
      lDistanceMatrix_.setDistanceType(DistanceType);
      lDistanceMatrix_.setTransposeMatrix(TransposeMatrix);
      lDistanceMatrix_.setNumberOfPoints(numberOfPoints);
      lDistanceMatrix_.setNumberOfScalarFields(numberOfScalarFields);
      lDistanceMatrix_.setInputScalarFields(&inputScalarFields);
      lDistanceMatrix_.setOutputDistanceMatrix(distanceMatrix->GetVoidPointer(0));

      switch(representative->GetDataType()){
        vtkTemplateMacro(({
              lDistanceMatrix_.execute<VTK_TT>();
              }));
      }

      vtkSmartPointer<vtkUnstructuredGrid> vtu=vtkSmartPointer<vtkUnstructuredGrid>::New();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!vtu){
        cerr << "[ttkLDistanceMatrix] Error: vtkUnstructuredGrid memory allocation problem." << endl;
        return -1;
      }
#endif

      {
        const int size=(TransposeMatrix)?numberOfPoints:numberOfScalarFields;
        vtkSmartPointer<vtkPoints> pts=vtkSmartPointer<vtkPoints>::New();
        pts->SetNumberOfPoints(size*size);

        const float delta=1.0f/size;
        vtkIdType ids[3];
        for(int i=0; i<size; ++i){
          const float x=i*delta;
          for(int j=0; j<size; ++j){
            const float y=j*delta;

            const int id=i*size+j;
            pts->SetPoint(id, x, y, 0.0f);

            if(i<size-1 and j<size-1){
              ids[0]=id;
              ids[1]=id+1;
              ids[2]=id+size+1;
              vtu->InsertNextCell(VTK_TRIANGLE,3,ids);

              ids[1]=id+size;
              vtu->InsertNextCell(VTK_TRIANGLE,3,ids);
            }
          }
        }
        vtu->SetPoints(pts);
      }

      vtkPointData* pointData=vtu->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!pointData){
        cerr << "[ttkLDistanceMatrix] Error: output vtkImageData has no point data." << endl;
        return -1;
      }
#endif
      pointData->AddArray(distanceMatrix);

      vtkFieldData* fieldData=vtu->GetFieldData();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!fieldData){
        cerr << "[ttkLDistanceMatrix] Error: output vtkImageData has no field data." << endl;
        return -1;
      }
#endif

      if(TransposeMatrix){
        vtkPointData* inputPointData=input->GetPointData();

        const int numberOfArrays=inputPointData->GetNumberOfArrays();
        for(int i=0; i<numberOfArrays; ++i){
          vtkDataArray* arr=inputPointData->GetArray(i);
          if(arr) fieldData->AddArray(arr);
        }
      }
      else
        fieldData->AddArray(metadata);

      output->ShallowCopy(vtu);
    }

    {
      stringstream msg;
      msg << "[ttkLDistanceMatrix] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
    }

    return 0;
  }
