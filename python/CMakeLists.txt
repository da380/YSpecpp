# The Python module yspec._core, built by nanobind. In a wheel build,
# scikit-build-core finds nanobind and installs the module into the package;
# in a plain CMake build (-DYSPEC_PYTHON=ON), nanobind is found through the
# Python that CMake finds, which must have it installed.

find_package(Python 3.10 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

if(NOT nanobind_ROOT)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE nanobind_ROOT)
endif()
find_package(nanobind CONFIG REQUIRED)

# The stable ABI where Python has one for nanobind (3.12 on), so that one
# wheel serves every later Python too.
nanobind_add_module(_core STABLE_ABI NB_STATIC src/bindings.cpp)
target_link_libraries(_core PRIVATE YSpec::YSpec)
target_compile_definitions(_core PRIVATE YSPEC_VERSION="${PROJECT_VERSION}")
if(YSPEC_WITH_OPENMP)
  target_link_libraries(_core PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS _core LIBRARY DESTINATION yspec)

# The PREM decks of the examples, which yspec.Model.prem reads.
install(FILES ${PROJECT_SOURCE_DIR}/examples/prem.200
              ${PROJECT_SOURCE_DIR}/examples/prem.200.no
        DESTINATION yspec/models)
