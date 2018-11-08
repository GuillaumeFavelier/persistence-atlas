#include<SpectralEmbedding.h>
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)

SpectralEmbedding::SpectralEmbedding():
  matrixDimension_{0},
  embeddingDimension_{0},
  minimumNumberOfComponents_{0},
  maximumNumberOfComponents_{0},
  numberOfNeighbors_{5},
  distanceMatrix_{nullptr},
  embedding_{nullptr},
  eigenvalues_{nullptr},
  components_{nullptr}
{
  auto finalize_callback=[](){
    Py_Finalize();
  };

  if(!Py_IsInitialized()){
    Py_Initialize();
    atexit(finalize_callback);
  }

  const char* version=Py_GetVersion();
  if(version[0]>='3'){
    stringstream msg;
    msg << "[SpectralEmbedding] Initializing Python: " <<
      version[0] << version[1] << version[2] << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  else{
    cerr << "[SpectralEmbedding] Error: Python 3+ is required:\n" <<
      version << " is provided." << endl;
  }

  majorVersion_=version[0];
}

SpectralEmbedding::~SpectralEmbedding(){
}

int SpectralEmbedding::execute() const{
  Timer t;

  if(majorVersion_<'3') return -1;

#ifndef TTK_ENABLE_KAMIKAZE
  if(modulePath_.length()<=0) return -1;
  if(moduleName_.length()<=0) return -1;
  if(functionName_.length()<=0) return -1;
  if(!distanceMatrix_) return -1;
  if(!embedding_) return -1;
  if(!components_) return -1;
#endif

  // check embedding dimension
  int embeddingDimension=embeddingDimension_;
  embeddingDimension=std::max(2,embeddingDimension);
  embeddingDimension=std::min(3,embeddingDimension);
  if(embeddingDimension>=matrixDimension_-1)
    embeddingDimension=matrixDimension_-2;

  // check number of components
  int maximumNumberOfComponents=std::max(1,maximumNumberOfComponents_);
  maximumNumberOfComponents=std::min(matrixDimension_-2,maximumNumberOfComponents);

  int minimumNumberOfComponents=std::max(1,minimumNumberOfComponents_);
  minimumNumberOfComponents=std::min(maximumNumberOfComponents,minimumNumberOfComponents);

  // check number of neighbors
  int numberOfNeighbors=std::max(1,numberOfNeighbors_);

  // declared here to avoid crossing initialization with goto
  vector<PyObject*> gc;
  PyObject* pArray;
  PyObject* pPath;
  PyObject* pSys;
  PyObject* pName;
  PyObject* pModule;
  PyObject* pFunc;
  PyObject* pTrue;
  PyObject* pFalse;
  PyObject* pMinimumNumberOfComponents;
  PyObject* pMaximumNumberOfComponents;
  PyObject* pNumberOfNeighbors;
  PyObject* pSigma;
  PyObject* pJobs;
  // PyObject* pEmbedding;
  PyObject* pProjection;
  PyArrayObject* npArr;
  PyArrayObject* npEmbedding;
  PyArrayObject* npValues;
  PyArrayObject* npProjection;
  string modulePath;

  if(PyArray_API==NULL){
    import_array();
  }

  // convert the distance matrix into a NumPy array.
  const int numberOfDimensions=2;
  npy_intp dimensions[2]{matrixDimension_, matrixDimension_};

  pArray=PyArray_SimpleNewFromData(numberOfDimensions,
      dimensions, NPY_DOUBLE, distanceMatrix_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pArray){
    cerr << "[SpectralEmbedding] Python error: failed to convert the array." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pArray);

  npArr=reinterpret_cast<PyArrayObject*>(pArray);

  pSys=PyImport_ImportModule("sys");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pSys){
    cerr << "[SpectralEmbedding] Python error: failed to load the sys module." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pSys);

  pPath=PyObject_GetAttrString(pSys, "path");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pPath){
    cerr << "[SpectralEmbedding] Python error: failed to get the path variable." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pPath);

  if(modulePath_=="default")
    modulePath=VALUE(TTK_SCRIPTS_PATH);
  else
    modulePath=modulePath_;

  {
    stringstream msg;
    msg << "[SpectralEmbedding] Loading Python script from: "
      << modulePath << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  PyList_Append(pPath, PyUnicode_FromString(modulePath.data()));

  // set other parameters
  pFalse=PyLong_FromLong(0);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pFalse){
    cerr << "[SpectralEmbedding] Python error: cannot convert pFalse." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pFalse);

  pTrue=PyLong_FromLong(1);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pTrue){
    cerr << "[SpectralEmbedding] Python error: cannot convert pTrue." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pTrue);

  pMinimumNumberOfComponents=PyLong_FromLong(minimumNumberOfComponents);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pMinimumNumberOfComponents){
    cerr << "[SpectralEmbedding] Python error: cannot convert pMinimumNumberOfComponents." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pMinimumNumberOfComponents);

  pMaximumNumberOfComponents=PyLong_FromLong(maximumNumberOfComponents);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pMaximumNumberOfComponents){
    cerr << "[SpectralEmbedding] Python error: cannot convert pMaximumNumberOfComponents." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pMaximumNumberOfComponents);

  pNumberOfNeighbors=PyLong_FromLong(numberOfNeighbors);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pNumberOfNeighbors){
    cerr << "[SpectralEmbedding] Python error: cannot convert pNumberOfNeighbors." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pNumberOfNeighbors);

  pSigma=PyLong_FromLong(int(sigma_));
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pSigma){
    cerr << "[SpectralEmbedding] Python error: cannot convert pSigma." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pSigma);

  pJobs=PyLong_FromLong(threadNumber_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pJobs){
    cerr << "[SpectralEmbedding] Python error: cannot convert pJobs." << endl;
    goto collect_garbage;
  }
#endif

  // load module
  pName=PyUnicode_FromString(moduleName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pName){
    cerr << "[SpectralEmbedding] Python error: moduleName parsing failed." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pName);

  pModule=PyImport_Import(pName);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pModule){
    cerr << "[SpectralEmbedding] Python error: module import failed." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pModule);

  // configure embbeding function
  pFunc=PyObject_GetAttrString(pModule, functionName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pFunc){
    cerr << "[SpectralEmbedding] Python error: functionName parsing failed." << endl;
    goto collect_garbage;
  }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
  if(!PyCallable_Check(pFunc)){
    cerr << "[SpectralEmbedding] Python error: function call failed." << endl;
    goto collect_garbage;
  }
#endif

  /*
  pEmbedding=PyObject_CallFunctionObjArgs(pFunc, npArr, parameterF,
      parameterF, pNumberOfNeighbors, pSigma, pTrue, pJobs, NULL);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pEmbedding){
    cerr << "[SpectralEmbedding] Python error: function returned invalid object." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pEmbedding);

  npEmbedding=reinterpret_cast<PyArrayObject*>(PyList_GetItem(pEmbedding, 1));
#ifndef TTK_ENABLE_KAMIKAZE
  if(PyArray_NDIM(npEmbedding) != numberOfDimensions) {
    cerr << "[SpectralEmbedding] Python error: returned object is not an array or have wrong dimension." << endl;
    goto collect_garbage;
  }
#endif

  // convert embedding back to C++ array
  {
    double* embedding=static_cast<double*>(embedding_);

    const int size=PyArray_SIZE(npEmbedding);
    if(size > 0){
      if(PyArray_TYPE(npEmbedding) == NPY_DOUBLE){
        double* c_out=reinterpret_cast<double*>(PyArray_DATA(npEmbedding));
        const int numberOfTuples=PyArray_SHAPE(npEmbedding)[0];

        for(int i=0; i<numberOfTuples; i++){
          embedding[3*i]=c_out[i];
          embedding[3*i+1]=c_out[i+numberOfTuples];

          if(embeddingDimension==2)
            embedding[3*i+2]=0;
          else if(embeddingDimension==3)
            embedding[3*i+2]=c_out[i+numberOfTuples*2];
        }
      }
    }
  }
  */

  pProjection=PyObject_CallFunctionObjArgs(pFunc, npArr, pMinimumNumberOfComponents,
      pMaximumNumberOfComponents, pNumberOfNeighbors, pSigma, pJobs, NULL);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pProjection){
    cerr << "[SpectralEmbedding] Python error: function returned invalid object." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pProjection);

  npValues=reinterpret_cast<PyArrayObject*>(PyList_GetItem(pProjection, 0));
#ifndef TTK_ENABLE_KAMIKAZE
  if(PyArray_SIZE(npValues)<=0){
    cerr << "[SpectralEmbedding] Python error: returned object is not an array or have wrong dimension." << endl;
    goto collect_garbage;
  }
#endif

  // convert eigenvalues back to C++ array
  {
    if(PyArray_TYPE(npValues) == NPY_DOUBLE){
      double* c_out=reinterpret_cast<double*>(PyArray_DATA(npValues));
      double* eigenvalues=static_cast<double*>(eigenvalues_);

      std::fill(eigenvalues, eigenvalues+matrixDimension_, std::numeric_limits<double>::lowest());
      const int size=PyArray_SIZE(npValues);
      for(int i=0; i<size; ++i)
        eigenvalues[i]=c_out[i];
    }
  }

  npProjection=reinterpret_cast<PyArrayObject*>(PyList_GetItem(pProjection, 1));
#ifndef TTK_ENABLE_KAMIKAZE
  if(PyArray_NDIM(npProjection) != numberOfDimensions) {
    cerr << "[SpectralEmbedding] Python error: returned object is not an array or have wrong dimension." << endl;
    goto collect_garbage;
  }
#endif

  // convert projection back to C++ array
  {
    const int size=PyArray_SIZE(npProjection);
    if(size > 0){
      const int numberOfComponents=PyArray_SHAPE(npProjection)[1];

      vector<double*> components(numberOfComponents);
      for(int i=0; i<numberOfComponents; ++i)
        components[i]=static_cast<double*>((*components_)[i]);
      components_->resize(numberOfComponents);

      if(PyArray_TYPE(npProjection) == NPY_DOUBLE){
        double* c_out=reinterpret_cast<double*>(PyArray_DATA(npProjection));
        const int numberOfTuples=PyArray_SHAPE(npProjection)[0];

        for(int i=0; i<numberOfTuples; i++){
          for(int j=0; j<numberOfComponents; j++)
            components[j][i]=c_out[i+j*numberOfTuples];
        }
      }
    }
  }

  npEmbedding=reinterpret_cast<PyArrayObject*>(PyList_GetItem(pProjection, 2));
#ifndef TTK_ENABLE_KAMIKAZE
  if(PyArray_NDIM(npEmbedding) != numberOfDimensions) {
    cerr << "[SpectralEmbedding] Python error: returned object is not an array or have wrong dimension." << endl;
    goto collect_garbage;
  }
#endif

  // convert embedding back to C++ array
  {
    double* embedding=static_cast<double*>(embedding_);

    const int size=PyArray_SIZE(npEmbedding);
    if(size > 0){
      if(PyArray_TYPE(npEmbedding) == NPY_DOUBLE){
        double* c_out=reinterpret_cast<double*>(PyArray_DATA(npEmbedding));
        const int numberOfTuples=PyArray_SHAPE(npEmbedding)[0];

        for(int i=0; i<numberOfTuples; i++){
          embedding[3*i]=c_out[i];
          embedding[3*i+1]=c_out[i+numberOfTuples];

          if(embeddingDimension==2)
            embedding[3*i+2]=0;
          else if(embeddingDimension==3)
            embedding[3*i+2]=c_out[i+numberOfTuples*2];
        }
      }
    }
  }

#ifndef TTK_ENABLE_KAMIKAZE
collect_garbage:
#endif
  for(auto i : gc)
    Py_DECREF(i);

  {
    stringstream msg;
    msg << "[SpectralEmbedding] Data-set "
      << " processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

