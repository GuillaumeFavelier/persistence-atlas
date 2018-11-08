ttk_find_python()

ttk_add_vtk_library(ttkClustering
	SOURCES ttkClustering.cpp
	HEADERS ttkClustering.h
	LINK clustering ttkTriangulation)
