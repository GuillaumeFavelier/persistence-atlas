ttk_find_python()

ttk_add_vtk_library(ttkSpectralEmbedding
	SOURCES ttkSpectralEmbedding.cpp
	HEADERS ttkSpectralEmbedding.h
	LINK spectralEmbedding ttkTriangulation)
