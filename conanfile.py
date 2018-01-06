import os
from conans import ConanFile, CMake
from conans.tools import download
from conans.tools import unzip



class UbitrackCoreConan(ConanFile):
    name = "ubitrack_core"
    version = "1.3.0"
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    options = {"shared": [True, False], "enable_tracing": [True, False]}
    requires = (
        "Boost/[>=1.64.0]@ulricheck/stable",

        "clapack/[>=3.2.1]@ulricheck/stable", 
        "msgpack/[>=2.1.5]@ulricheck/stable", 

        "ubitrack_boost_bindings/1.0@ulricheck/stable", 
        "ubitrack_tinyxml/2.5.3@ulricheck/stable", 
        "ubitrack_log4cpp/0.3.5@ulricheck/stable",
        )

    default_options = (
        "shared=True",
        "enable_tracing=True",
        )

    # all sources are deployed with the package
    exports_sources = "cmake/*", "doc/*", "misc/*", "src/*", "tests/*", "CMakeLists.txt"

    # UbitrackConfig.cmake should be reused by all depending packages
    exports = "cmake/UbitrackConfig.cmake"

    def configure(self):
        if self.options.shared:
            self.options['Boost'].shared = True 
            self.options['clapack'].shared = True 
            self.options['msgpack'].shared = True 
            self.options['ubitrack_log4cpp'].shared = True


    def imports(self):
        self.copy(pattern="*.dll", dst="bin", src="bin") # From bin to bin
        self.copy(pattern="*.dylib*", dst="lib", src="lib") 
       
    def build(self):
        cmake = CMake(self)
        cmake.definitions['BUILD_SHARED_LIBS'] = self.options.shared
        cmake.definitions['ENABLE_TRACING'] = self.options.enable_tracing
        cmake.configure()
        cmake.build()
        cmake.install()

    def package(self):
        # self.copy("*.h", dst="include", src="src")
        # self.copy("*.lib", dst="lib", keep_path=False)
        # self.copy("*.dll", dst="bin", keep_path=False)
        # self.copy("*.dylib*", dst="lib", keep_path=False)
        # self.copy("*.so", dst="lib", keep_path=False)
        # self.copy("*.a", dst="lib", keep_path=False)
        # self.copy("*", dst="bin", src="bin", keep_path=False)
        self.copy("UbitrackConfig.cmake", dst="cmake", src="cmake")

    def package_info(self):
        if self.options.enable_tracing:
            self.cpp_info.defines.append("ENABLE_EVENT_TRACING")
            if self.settings.os == "Windows":
                self.cpp_info.defines.append("HAVE_ETW")
            elif self.settings.os == "Linux":
                self.cpp_info.defines.append("HAVE_LTTNGUST")
            elif self.settings.os == "Macos":    
                self.cpp_info.defines.append("HAVE_DTRACE")


        suffix = ""
        if self.settings.build_type == "Debug" and self.settings.os == "Windows":
            suffix = "d"
        self.cpp_info.libs.append("utcore%s" % (suffix))
