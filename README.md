## utcore
This is the utcore Ubitrack submodule.

## Description

The utcore contains fundamental datastructures and algorithms. Needed by all other modules.

## Dependencies

There are no dependencies for this module.

## For Users: Use this package

### Basic setup

    $ conan install ubitrack_core/1.3.0@ubitrack/stable

### Project setup

If you handle multiple dependencies in your project is better to add a *conanfile.txt*

    [requires]
    ubitrack_core/1.3.0@ubitrack/stable

    [generators]
    cmake
    txt

Complete the installation of requirements for your project running:

    $ mkdir build && cd build && conan install ..
    
Note: It is recommended that you run conan install from a build directory and not the root of the project directory.  This is because conan generates *conanbuildinfo* files specific to a single build configuration which by default comes from an autodetected default profile located in ~/.conan/profiles/default .  If you pass different build configuration options to conan install, it will generate different *conanbuildinfo* files.  Thus, they shoudl not be added to the root of the project, nor committed to git. 

## For Packagers: Publish this Package

The example below shows the commands used to publish to ulricheck conan repository. To publish to your own conan respository (for example, after forking this git repository), you will need to change the commands below accordingly. 

## Build  

This is a header only library, so nothing needs to be built.

## Package 

    $ conan create . ubitrack/stable
    
## Add Remote

    $ conan remote add camp "https://conan.campar.in.tum.de" True

## Upload

    $ conan upload -r camp ubitrack_core/1.3.0@ubitrack/stable
