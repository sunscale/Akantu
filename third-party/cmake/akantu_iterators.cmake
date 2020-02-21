set(AKANTU_ITERATORS_TARGETS_EXPORT ${AKANTU_TARGETS_EXPORT})

add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/akantu_iterators)

set(AKANTU_ITERATORS_VERSION "master" CACHE INTERNAL "")
set(AKANTU_ITERATORS_LIBRARIES "akantu_iterators" CACHE INTERNAL "")
package_add_to_export_list(akantu_iterators akantu_iterators)
