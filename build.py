from conan.packager import ConanMultiPackager
import platform

if __name__ == "__main__":
    builder = ConanMultiPackager()
    builder.add_common_builds(shared_option_name="ubitrack_core:shared", pure_c=True)

    if platform.system() == "Windows":
        filtered_builds = []
        for settings, options, env_vars, build_requires, reference in builder.items:
            if settings["compiler"] != "Visual Studio" or options[name + ":shared"]:
                filtered_builds.append([settings, options, env_vars, build_requires, reference])
        builder.builds = filtered_builds

    builder.run()
