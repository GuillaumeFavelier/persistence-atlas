ttk_find_python()

ttk_add_base_library(clustering
	SOURCES Clustering.cpp
	HEADERS Clustering.h
    LINK common)

install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/clustering.py DESTINATION
  ${CMAKE_INSTALL_PREFIX}/share/scripts/ttk)

add_definitions(-DTTK_SCRIPTS_PATH=${CMAKE_INSTALL_PREFIX}/share/scripts/ttk)
