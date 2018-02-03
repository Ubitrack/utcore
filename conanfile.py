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
    options = {"shared": [True, False], "enable_tracing": [True, False]}
    requires = (
        "Boost/[>=1.64.0]@camposs/stable",

        "clapack/[>=3.2.1]@camposs/stable", 
        "msgpack/[>=2.1.5]@camposs/stable", 

        "ubitrack_boost_bindings/1.0@camposs/stable", 
        "ubitrack_tinyxml/2.5.3@camposs/stable", 
        "ubitrack_log4cpp/0.3.5@camposs/stable",
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

    def system_requirements(self):
        if self.options.enable_tracing:
            if os_info.is_linux:
                if os_info.with_apt:
                    installer = SystemPackageTool()
                    installer.update()
                    installer.install("python-software-properties")
                    try:
                        # XXX bit of a hack here .. it's using private members of package_tools ...
                        self.run("%s/usr/bin/apt-add-repository -y -u ppa:lttng/stable" % (installer._tool._sudo_str))
                    except:
                        self.output.warn("Could not add PPA for LTTNG 2.10!")

                    installer.update()
                    if self.settings.arch == "x86" and tools.detected_architecture() == "x86_64":
                        arch_suffix = ':i386'
                        installer.install("g++-multilib")
                    else:
                        arch_suffix = ''
                    installer.install("%s%s" % ("lttng-tools", arch_suffix))
                    # installer.install("lttng-modules-dkms")
                    installer.install("%s%s" % ("liblttng-ust-dev", arch_suffix))
                elif os_info.with_yum:
                    installer = SystemPackageTool()
                    if self.settings.arch == "x86" and tools.detected_architecture() == "x86_64":
                        arch_suffix = '.i686'
                        installer.install("glibc-devel.i686")
                    else:
                        arch_suffix = ''
                    installer.install("%s%s" % ("lttng-tools", arch_suffix))
                    installer.install("%s%s" % ("lttng-ust", arch_suffix))
                else:
                    self.output.warn("Could not determine package manager, skipping Linux system requirements installation.")

    def imports(self):
        self.copy(pattern="*.dll", dst="bin", src="bin") # From bin to bin
        self.copy(pattern="*.dylib*", dst="lib", src="lib") 
       
    def build(self):
        cmake = CMake(self)
        cmake.definitions['BUILD_SHARED_LIBS'] = self.options.shared
        cmake.definitions['ENABLE_TRACING'] = self.options.enable_tracing
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
