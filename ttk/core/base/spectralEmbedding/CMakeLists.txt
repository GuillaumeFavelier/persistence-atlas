ttk_find_python()

ttk_add_base_library(spectralEmbedding
	SOURCES SpectralEmbedding.cpp
	HEADERS SpectralEmbedding.h
    LINK common)

install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/spectralEmbedding.py DESTINATION
  ${CMAKE_INSTALL_PREFIX}/share/scripts/ttk)

add_definitions(-DTTK_SCRIPTS_PATH=${CMAKE_INSTALL_PREFIX}/share/scripts/ttk)
