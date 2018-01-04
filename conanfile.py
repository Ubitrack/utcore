import os
from conans import ConanFile, CMake
from conans.tools import download
from conans.tools import unzip


class UbitrackCoreConan(ConanFile):
    name = "ubitrack_core"
    version = "1.3.0"
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    requires = (
        "Boost.Bind/1.65.1@bincrafters/stable", 
        "Boost.Chrono/1.65.1@bincrafters/stable", 
        "Boost.Core/1.65.1@bincrafters/stable", 
        "Boost.Date_Time/1.65.1@bincrafters/stable",
        "Boost.Filesystem/1.65.1@bincrafters/stable", 
        "Boost.Iostreams/1.65.1@bincrafters/stable", 
        "Boost.Locale/1.65.1@bincrafters/stable", 
        "Boost.Math/1.65.1@bincrafters/stable", 
        "Boost.Mpl/1.65.1@bincrafters/stable", 
        "Boost.Numeric_Ublas/1.65.1@bincrafters/stable",
        "Boost.Program_Options/1.65.1@bincrafters/stable",
        "Boost.Random/1.65.1@bincrafters/stable", 
        "Boost.Regex/1.65.1@bincrafters/stable", 
        "Boost.Serialization/1.65.1@bincrafters/stable",
        "Boost.System/1.65.1@bincrafters/stable",
        "Boost.Test/1.65.1@bincrafters/stable",
        "Boost.Type_Traits/1.65.1@bincrafters/stable",
        "Boost.Utility/1.65.1@bincrafters/stable",

        "clapack/3.2.1@ulricheck/stable", 
        "msgpack/2.1.5@ulricheck/stable", 
        "ubitrack_boost_bindings/1.0@ulricheck/stable", 
        "ubitrack_tinyxml/2.5.3@ulricheck/stable", 
        "ubitrack_log4cpp/0.3.5@ulricheck/stable",
        )

    default_options = (
        "Boost.Chrono:shared=True", 
        "Boost.Date_Time:shared=True", 
        "Boost.Filesystem:shared=True", 
        "Boost.Iostreams:shared=True",  
        "Boost.Locale:shared=True", 
        "Boost.Program_Options:shared=True", 
        "Boost.Regex:shared=True", 
        "Boost.Serialization:shared=True", 
        "Boost.System:shared=True", 
        "opencv:shared=True", 
        "clapack:shared=True", 
        "msgpack:shared=True", 
        "ubitrack_log4cpp:shared=True",
        )

    # all sources are deployed with the package
    exports_sources = "cmake/*", "doc/*", "misc/*", "src/*", "tests/*", "CMakeLists.txt"

    def imports(self):
        self.copy(pattern="*.dll", dst="bin", src="bin") # From bin to bin
        self.copy(pattern="*.dylib*", dst="bin", src="lib") 
       
    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        self.copy("*.h", dst="include", src="src")
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.dylib*", dst="lib", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)

    def package_info(self):
        suffix = ""
        if self.settings.build_type == "Debug" and self.settings.os == "Windows":
            suffix = "d"
        self.cpp_info.libs.append("utcore%s" % (suffix))
