if(NOT ParaView_SOURCE_DIR)
  cmake_minimum_required(VERSION 2.8)
  find_package(MPI REQUIRED)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

find_library(ZOLTAN_LIBRARY zoltan)
find_path(ZOLTAN_INCLUDE_DIR zoltan.h)

include_directories(${MPI_INCLUDE_PATH})
include_directories(${ZOLTAN_INCLUDE_DIR})
include_directories(${ParaView_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(servermanager_sources
  vtkParticlePartitionFilter.cxx
  vtkBoundsExtentTranslator.cxx)
set(servermanager_xml
  vtkParticlePartitionFilter.xml)
set(server_deps
  ${VTK_LIBRARIES}
  ${ZOLTAN_LIBRARY})

set(library_name ParticlePartitionPlugin)

add_paraview_plugin(${library_name} "1.0"
  SERVER_MANAGER_XML ${servermanager_xml}
  SERVER_MANAGER_SOURCES ${servermanager_sources}
  )

target_link_libraries(${library_name} ${server_deps})
