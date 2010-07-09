set (NAME grid)
set (SOURCES 
grid.cc
grid_averaging_metrics.cc
grid_create_elements.cc
grid_partition.cc
grid_readers.cc
grid_transform.cc
)

add_library(${NAME} STATIC ${SOURCES} )

install (FILES ${NAME}.h DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 
