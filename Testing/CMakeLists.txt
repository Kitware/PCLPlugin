set(python_tests
  testPCLFilters.py
)

set(python_test args
)

set(vtkpython_exe ${ParaView_BINARY_DIR}/bin/vtkpython)

macro(add_python_test name)
  get_filename_component(base_name ${name} NAME_WE)
  set(test_name pcl_plugin_test_${base_name})
  add_test(${test_name} ${vtkpython_exe} ${CMAKE_CURRENT_SOURCE_DIR}/${name} ${python_test_args} ${${base_name}_extra_args})
endmacro()


foreach(name ${python_tests})
  add_python_test(${name})
endforeach()
