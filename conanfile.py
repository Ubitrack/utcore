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
        "Boost.Config/1.65.1@bincrafters/stable",
        "Boost.Predef/1.65.1@bincrafters/stable",
        "Boost.Preprocessor/1.65.1@bincrafters/stable",
        "Boost.Assert/1.65.1@bincrafters/stable",
        "Boost.Io/1.65.1@bincrafters/stable",
        "Boost.Static_Assert/1.65.1@bincrafters/stable",
        "Boost.Winapi/1.65.1@bincrafters/stable",
        "Boost.Core/1.65.1@bincrafters/stable",
        "Boost.Throw_Exception/1.65.1@bincrafters/stable",
        "Boost.Array/1.65.1@bincrafters/stable",
        "Boost.Bind/1.65.1@bincrafters/stable",
        "Boost.Integer/1.65.1@bincrafters/stable",
        "Boost.Logic/1.65.1@bincrafters/stable",
        "Boost.Move/1.65.1@bincrafters/stable",
        "Boost.System/1.65.1@bincrafters/stable",
        "Boost.Type_Traits/1.65.1@bincrafters/stable",
        "Boost.Atomic/1.65.1@bincrafters/stable",
        "Boost.Smart_Ptr/1.65.1@bincrafters/stable",
        "Boost.Tuple/1.65.1@bincrafters/stable",
        "Boost.Exception/1.65.1@bincrafters/stable",
        "Boost.Level5Group/1.65.1@bincrafters/stable",
        "Boost.Concept_Check/1.65.1@bincrafters/stable",
        "Boost.Conversion/1.65.1@bincrafters/stable",
        "Boost.Detail/1.65.1@bincrafters/stable",
        "Boost.Function/1.65.1@bincrafters/stable",
        "Boost.Function_Types/1.65.1@bincrafters/stable",
        "Boost.Functional/1.65.1@bincrafters/stable",
        "Boost.Fusion/1.65.1@bincrafters/stable",
        "Boost.Iterator/1.65.1@bincrafters/stable",
        "Boost.Mpl/1.65.1@bincrafters/stable",
        "Boost.Optional/1.65.1@bincrafters/stable",
        "Boost.Type_Index/1.65.1@bincrafters/stable",
        "Boost.Typeof/1.65.1@bincrafters/stable",
        "Boost.Utility/1.65.1@bincrafters/stable",
        "Boost.Endian/1.65.1@bincrafters/stable",
        "Boost.Intrusive/1.65.1@bincrafters/stable",
        "Boost.Lambda/1.65.1@bincrafters/stable",
        "Boost.Numeric_Conversion/1.65.1@bincrafters/stable",
        "Boost.Numeric_Interval/1.65.1@bincrafters/stable",
        "Boost.Rational/1.65.1@bincrafters/stable",
        "Boost.Regex/1.65.1@bincrafters/stable",
        "Boost.Tokenizer/1.65.1@bincrafters/stable",
        "Boost.Tti/1.65.1@bincrafters/stable",
        "Boost.Container/1.65.1@bincrafters/stable",
        "Boost.Range/1.65.1@bincrafters/stable",
        "Boost.Ratio/1.65.1@bincrafters/stable",
        "Boost.Chrono/1.65.1@bincrafters/stable",
        "Boost.Filesystem/1.65.1@bincrafters/stable",
        "Boost.Foreach/1.65.1@bincrafters/stable",
        "Boost.Level8Group/1.65.1@bincrafters/stable",
        "Boost.Proto/1.65.1@bincrafters/stable",
        "Boost.Unordered/1.65.1@bincrafters/stable",
        "Boost.Algorithm/1.65.1@bincrafters/stable",
        "Boost.Lexical_Cast/1.65.1@bincrafters/stable",
        "Boost.Math/1.65.1@bincrafters/stable",
        "Boost.Phoenix/1.65.1@bincrafters/stable",
        "Boost.Timer/1.65.1@bincrafters/stable",
        "Boost.Random/1.65.1@bincrafters/stable",
        "Boost.Test/1.65.1@bincrafters/stable",
        "Boost.Variant/1.65.1@bincrafters/stable",
        "Boost.Iostreams/1.65.1@bincrafters/stable",
        "Boost.Level11Group/1.65.1@bincrafters/stable",
        "Boost.Serialization/1.65.1@bincrafters/stable",
        "Boost.Numeric_Ublas/1.65.1@bincrafters/stable",

        "clapack/3.2.1@ulricheck/stable", 
        "msgpack/2.1.5@ulricheck/stable", 
        "ubitrack_boost_bindings/1.0@ulricheck/stable", 
        "ubitrack_tinyxml/2.5.3@ulricheck/stable", 
        "ubitrack_log4cpp/0.3.5@ulricheck/stable",
        )

    default_options = (
        "Boost.Chrono:shared=True", 
        "Boost.Filesystem:shared=True", 
        "Boost.Iostreams:shared=True",  
        "Boost.Program_Options:shared=True", 
        "Boost.Regex:shared=True", 
        "Boost.System:shared=True", 
        "clapack:shared=True", 
        "msgpack:shared=True", 
        "ubitrack_log4cpp:shared=True",
        "shared=True",
        "enable_tracing=True",
        )

    # all sources are deployed with the package
    exports_sources = "cmake/*", "doc/*", "misc/*", "src/*", "tests/*", "CMakeLists.txt"

    # UbitrackConfig.cmake should be reused by all depending packages
    exports = "cmake/UbitrackConfig.cmake"

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
