from conans import ConanFile, CMake, tools
from conans.tools import os_info, SystemPackageTool

class UbitrackCoreConan(ConanFile):
    name = "ubitrack_core"
    version = "1.3.0"

    description = "Ubitrack Core Library"
    url = "https://github.com/Ubitrack/utcore.git"
    license = "GPL"

    short_paths = True
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    options = {
    "shared": [True, False], 
    "enable_tracing": [True, False], 
    "without_tests": [True, False]
    }
    
    requires = (
        "Boost/[>=1.59.0,<1.65.0]@camposs/stable",

        "clapack/[>=3.2.1]@camposs/stable", 
        "msgpack/[>=2.1.5]@camposs/stable", 

        "ubitrack_boost_bindings/1.0@ubitrack/stable", 
        "ubitrack_tinyxml/2.5.3@ubitrack/stable", 
        "ubitrack_log4cpp/0.3.5@ubitrack/stable",
        )

    default_options = (
        "shared=True",
        "enable_tracing=False",
        "without_tests=True"
        )

    # all sources are deployed with the package
    exports_sources = "cmake/*", "doc/*", "misc/*", "src/*", "tests/*", "CMakeLists.txt", "utcoreConfig.cmake"

    # UbitrackConfig.cmake should be reused by all depending packages
    exports = "cmake/UbitrackConfig.cmake"

    def configure(self):
        # Boost
        self.options["Boost"].without_atomic = True
        self.options["Boost"].without_container = True
        self.options["Boost"].without_context = True
        self.options["Boost"].without_coroutine = True
        self.options["Boost"].without_coroutine2 = True
        self.options["Boost"].without_exception = True
        # self.options["Boost"].without_fiber = True
        self.options["Boost"].without_graph = True
        self.options["Boost"].without_graph_parallel = True
        self.options["Boost"].without_locale = True
        self.options["Boost"].without_log = True
        # self.options["Boost"].without_metaparse = True
        self.options["Boost"].without_mpi = True
        self.options["Boost"].without_signals = True
        self.options["Boost"].without_timer = True
        # self.options["Boost"].without_type_erasure = True
        self.options["Boost"].without_wave = True

        if self.options.shared:
            self.options['Boost'].shared = True 
            self.options['clapack'].shared = True 
            self.options['msgpack'].shared = True 
            self.options['ubitrack_log4cpp'].shared = True
            self.options['zlib'].shared = True

    # @Todo Conan should never install system packages without asking !!!
    # def system_requirements(self):
    #     if self.options.enable_tracing:
    #         if os_info.is_linux:
    #             if os_info.with_apt:
    #                 installer = SystemPackageTool()
    #                 installer.update()
    #                 installer.install("python-software-properties")
    #                 try:
    #                     # XXX bit of a hack here .. it's using private members of package_tools ...
    #                     self.run("%s/usr/bin/apt-add-repository -y -u ppa:lttng/stable-2.10" % (installer._tool._sudo_str))
    #                 except:
    #                     self.output.warn("Could not add PPA for LTTNG 2.10!")

    #                 installer.update()
    #                 if self.settings.arch == "x86" and tools.detected_architecture() == "x86_64":
    #                     arch_suffix = ':i386'
    #                     installer.install("g++-multilib")
    #                 else:
    #                     arch_suffix = ''
    #                 installer.install("%s%s" % ("lttng-tools", arch_suffix))
    #                 # installer.install("lttng-modules-dkms")
    #                 installer.install("%s%s" % ("liblttng-ust-dev", arch_suffix))
    #             elif os_info.with_yum:
    #                 installer = SystemPackageTool()
    #                 if self.settings.arch == "x86" and tools.detected_architecture() == "x86_64":
    #                     arch_suffix = '.i686'
    #                     installer.install("glibc-devel.i686")
    #                 else:
    #                     arch_suffix = ''
    #                 installer.install("%s%s" % ("lttng-tools", arch_suffix))
    #                 installer.install("%s%s" % ("lttng-ust", arch_suffix))
    #             else:
    #                 self.output.warn("Could not determine package manager, skipping Linux system requirements installation.")

    def imports(self):
        self.copy(pattern="*.dll", dst="bin", src="bin") # From bin to bin
        self.copy(pattern="*.dylib*", dst="lib", src="lib") 
       
    def build(self):
        cmake = CMake(self)
        cmake.definitions['BUILD_SHARED_LIBS'] = self.options.shared
        cmake.definitions['ENABLE_TRACING'] = self.options.enable_tracing
        cmake.definitions['ENABLE_UNITTESTS'] = not self.options.without_tests
        cmake.definitions['CMAKE_POSITION_INDEPENDENT_CODE'] = True
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
        if self.settings.os == "Windows":
            suffix += self.version.replace(".", "")
            if self.settings.build_type == "Debug":
                suffix += "d"
        self.cpp_info.libs.append("utcore%s" % (suffix))
