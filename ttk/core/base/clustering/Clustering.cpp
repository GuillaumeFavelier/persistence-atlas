#include<Clustering.h>
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)

Clustering::Clustering():
  numberOfObjects_{0},
  maximumNumberOfClusters_{0},
  minimumNumberOfClusters_{0},
  currentGap_{0},
  numberOfComponents_{0},
  maps_{nullptr},
  eigenvalues_{nullptr},
  assignation_{nullptr}
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
    msg << "[Clustering] Initializing Python: " <<
      version[0] << version[1] << version[2] << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  else{
    cerr << "[Clustering] Error: Python 3+ is required:\n" <<
      version << " is provided." << endl;
  }

  majorVersion_=version[0];
}

Clustering::~Clustering(){
}

int Clustering::execute() const{
  Timer t;

  if(majorVersion_<'3') return -1;

#ifndef TTK_ENABLE_KAMIKAZE
  if(modulePath_.length()<=0) return -1;
  if(moduleName_.length()<=0) return -1;
  if(functionName_.length()<=0) return -1;
  if(!maps_) return -1;
  if(!assignation_) return -1;
#endif

  vector<double>& maps=*maps_;

  const int numberOfObjects=numberOfObjects_;

  // check number of clusters
  const int minimumNumberOfClusters=std::max(2, minimumNumberOfClusters_);
  const int maximumNumberOfClusters=std::max(minimumNumberOfClusters, maximumNumberOfClusters_);

  // check number of components
  int numberOfComponents=numberOfComponents_;
  if(numberOfComponents>=numberOfObjects_-1)
    numberOfComponents=numberOfObjects_-2;

  // declared here to avoid crossing initialization with goto
  vector<PyObject*> gc;
  vector<double> values;
  npy_intp mapDimensions[2];
  npy_intp valDimensions[2];
  PyObject* pArray;
  PyObject* pValues;
  PyObject* pPath;
  PyObject* pSys;
  PyObject* pName;
  PyObject* pModule;
  PyObject* pFunc;
  PyObject* pMaxClusters;
  PyObject* pMinClusters;
  PyObject* pCurrentGap;
  PyObject* pReturn;
  PyArrayObject* npArr;
  PyArrayObject* npReturn;
  string modulePath;

  if(PyArray_API==NULL){
    import_array();
  }

  // convert the maps into a NumPy array.
  mapDimensions[0]=numberOfObjects;
  mapDimensions[1]=numberOfComponents;
  pArray=PyArray_SimpleNewFromData(2,
      mapDimensions, NPY_DOUBLE, maps.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pArray){
    cerr << "[Clustering] Python error: failed to convert maps array." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pArray);

  npArr=reinterpret_cast<PyArrayObject*>(pArray);

  {
    double* eigenvalues=static_cast<double*>(eigenvalues_);
    for(int i=0; i<numberOfObjects_; ++i){
      if(eigenvalues[i]>std::numeric_limits<double>::lowest())
        values.push_back(eigenvalues[i]);
    }
  }

  valDimensions[0]=1;
  valDimensions[1]=values.size();
  pValues=PyArray_SimpleNewFromData(2,
      valDimensions, NPY_DOUBLE, values.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pValues){
    cerr << "[Clustering] Python error: failed to convert values array." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pValues);

  pSys=PyImport_ImportModule("sys");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pSys){
    cerr << "[Clustering] Python error: failed to load the sys module." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pSys);

  pPath=PyObject_GetAttrString(pSys, "path");
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pPath){
    cerr << "[Clustering] Python error: failed to get the path variable." << endl;
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
    msg << "[Clustering] Loading Python script from: "
      << modulePath << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
  PyList_Append(pPath, PyUnicode_FromString(modulePath.data()));

  // set other parameters
  pMaxClusters=PyLong_FromLong(maximumNumberOfClusters);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pMaxClusters){
    cerr << "[Clustering] Python error: cannot convert pMaxClusters." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pMaxClusters);

  pMinClusters=PyLong_FromLong(minimumNumberOfClusters);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pMinClusters){
    cerr << "[Clustering] Python error: cannot convert pMinClusters." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pMinClusters);

  pCurrentGap=PyLong_FromLong(currentGap_);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pCurrentGap){
    cerr << "[Clustering] Python error: cannot convert pCurrentGap." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pCurrentGap);


  // load module
  pName=PyUnicode_FromString(moduleName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pName){
    cerr << "[Clustering] Python error: moduleName parsing failed." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pName);

  pModule=PyImport_Import(pName);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pModule){
    cerr << "[Clustering] Python error: module import failed." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pModule);

  // configure function
  pFunc=PyObject_GetAttrString(pModule, functionName_.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pFunc){
    cerr << "[Clustering] Python error: functionName parsing failed." << endl;
    goto collect_garbage;
  }
  if(!PyCallable_Check(pFunc)){
    cerr << "[Clustering] Python error: function call failed." << endl;
    goto collect_garbage;
  }
#endif

  pReturn=PyObject_CallFunctionObjArgs(pFunc, npArr, pValues, pMinClusters, pMaxClusters, pCurrentGap, NULL);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!pReturn){
    cerr << "[Clustering] Python error: function returned invalid object." << endl;
    goto collect_garbage;
  }
#endif
  gc.push_back(pReturn);

  npReturn=reinterpret_cast<PyArrayObject*>(pReturn);
#ifndef TTK_ENABLE_KAMIKAZE
  if(PyArray_NDIM(npReturn) != 1) {
    cerr << "[Clustering] Python error: returned object is not an array or have wrong dimension." << endl;
    goto collect_garbage;
  }
#endif

  // convert back to C++ array
  {
    int* assignation=static_cast<int*>(assignation_);

    const int size=PyArray_SIZE(npReturn);
    if(size > 0){
      if(PyArray_TYPE(npReturn) == NPY_LONG){
        long int* c_out=reinterpret_cast<long int*>(PyArray_DATA(npReturn));

        for(int i=0; i<size; i++)
          assignation[i]=c_out[i];
      }
      else if(PyArray_TYPE(npReturn) == NPY_INT){
        int* c_out=reinterpret_cast<int*>(PyArray_DATA(npReturn));

        for(int i=0; i<size; i++)
          assignation[i]=c_out[i];
      }
      else if(PyArray_TYPE(npReturn) == NPY_LONGLONG){
        long long int* c_out=reinterpret_cast<long long int*>(PyArray_DATA(npReturn));

        for(int i=0; i<size; i++)
          assignation[i]=c_out[i];
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
    msg << "[Clustering] Data-set "
      << " processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

