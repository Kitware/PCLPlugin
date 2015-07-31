# PointCloudLibrary Plugin for ParaView v1.1

# Authors:
  Pat Marion
  Roland Kwitt
  Brad Davis
  Casey Goodlett


# Installation and Build on Linux/MacOSX


## Introduction
This page details the minimum commands required to download and build the PCL Plugin for ParaView. We assume basic knowledge of cmake and the usual Unix build process.

We propose a three-stage build process to compile and run the PCL plugin for ParaView:
* Building a vanilla ParaView
* Building PCL (with ParaView's VTK)
* Building the PCL plugin

While other strategies are possible, we suggest to follow this process (for now) to ensure that both PCL and ParaView use the same VTK version.


## Pre-requirements

* A computer connected to the Internet
* Git
* CCMake
* gcc
* build-essential


## Get ParaView 4.2.0 (and VTK) source code

* git clone https://github.com/Kitware/ParaView.git
* cd ParaView
* git checkout master
* git checkout tags/v4.2.0
* git submodule init
* git submodule update


## Get PCL 1.7.2 source code

* git clone https://github.com/PointCloudLibrary/pcl.git
* cd pcl
* git checkout tags/pcl-1.7.2


## Get the ParaView PCL Plugin source code

* git clone https://github.com/Kitware/PCLPlugin.git
* cd PCLPlugin
* git checkout master


## Install required libraries (if not already installed)

You may need to install the following libraries:
* libqt4-dev (version >= 4.7, we tsedt it with 4.8.6)
* qt4-dev-tools (version >= 4.7, we tsedt it with 4.8.6)
* libflann-dev
* libeigen3-dev
* libboost-all-dev
* libxt-dev


## Build ParaView

* Move to the ParaView source code drectory
* mkdir Build
* You may want to add "Build" to .gitignore to avoid committing the whole build directory if you were to commit
* Move to the Build drectory
* ccmake ..
* Press Configure to build up the options
* Set BUILD_SHARED_LIBS to ON
* Set CMAKE_BUILD_TYPE to Release
* Set BUILD_TESTING to OFF
* Press Configure two times, and then Generate (pay attention to the error
  messages listed at the end of the Configure step, since some needed libraries
  may not be installed on your system)
* make
* In case an ugly "vtksys/ios/iostream iostream.h no such file" were to appear,
  try deleting the Build directory, and restart the procedure from the "ccmake .."
  step
* Since all ParaView 4.2 libraries have the "-pv4.2" suffix attache dto their
  names, some symbolic links have to be added:
  cd <FullPathToParaViewSource>/Build/lib/
  for i in *-pv*.so; do ln -s $i ${i%-pv*}.so; done

This will build ParaView as well as all the submodules (such as VTK). Optionally, you can specify -jN, e.g., make -j4 to speed up the compilation process by using more cores.

For more instructions refer to the [PCL build instructions](http://pointclouds.org/downloads/source.html),
MacOSX users may consider [this page](http://www.pcl-users.org/Problems-compiling-from-source-on-Mac-OS-X-Lion-td3862799.html) if they encounter compilation errors.


## Build PCL

Once we have built a vanilla ParaView, we can build PCL against ParaView's VTK version.

* Move to the PCL source code directory
* mkdir Build
* You may want to add "Build" to .gitignore to avoid committing the whole build directory if you were to commit
* cd Build
* ccmake ..
* Press Configure to build up the options
* Set VTK_DIR to <FullPathToParaViewSource>/Build/VTK
* Set CMAKE_CXX_FLAGS to -L<FullPathToParaViewSource>/Build/lib
* Press Configure two times, and then Generate (pay attention to the error
  messages listed at the end of the Configure step, since some needed libraries
  may not be installed on your system)
* make

For more instructions refer to the [ParaView build instructions](http://www.paraview.org/Wiki/ParaView:Build_And_Install)


## Build the ParaView PCL Plugin

* mkdir Build
* You may want to add "Build" to .gitignore to avoid committing the whole build directory if you were to commit
* cd Build
* ccmake ..
* Press Configure to build up the options
* Set BUILD_SHARED_LIBS to ON
* Set PCL_DIR to <FullPathToPCLSource>/Build
* Set ParaView_DIR to <FullPathToParaViewSource>/Build
* Configure and Generate (pay attention to the error messages listed at the
  end of the Configure step, since some needed libraries may not be installed
  on your system)
* make


## Importing the PCL plugin in ParaView

Now that we have built the plugin, we are ready to import it in ParaView.

Start the ParaView you have just built
* cd <FullPathToParaViewSource>/Build
* ./bin/paraview

Go to Tools/Manage Plugins. Press the Load New... button and navigate to ParaView's ParaView/Build/bin directory. Navigate to the <FullPathToPCLPluginSource>/Build/bin directory and select libvtkPCLFilters.so to load. You might also want to enable Auto Load so that the plugin gets automatically loaded. You should now be ready to use PCL algorithms in ParaView.



  