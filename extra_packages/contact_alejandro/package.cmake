set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/packages")

include(contact)
include(cpparray)
include(nlopt)
include(optimization)

#add_example(contact              "Examples on how to use contact within Akantu"        PACKAGE contact)
#add_example(optimization         "Optimization examples"                               PACKAGE optimization)
