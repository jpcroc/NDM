function(doxygen_doc)
    set(DOC_OUTPUTDIR ${CMAKE_CURRENT_BINARY_DIR})
    set(DOC_LANGUAGE English)
    set(QUIET_DOXYGEN YES)

    get_target_property(SOURCE_LIST ${TARGET_NAME} SOURCES)
    string(REPLACE ";" " " SOURCE_LIST_2 "${SOURCE_LIST}")

    find_package(Doxygen)

    if (DOXYGEN_FOUND)
        set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/cmake_files/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message(STATUS "Doxygen build started")

        add_custom_target(doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT} 2> ${DOC_OUTPUTDIR}/err.log
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Generating documentation with Doxygen"
            VERBATIM
        )
    else (DOXYGEN_FOUND)
        message(STATUS "Doxygen needs to be installed to generate the Doxygen documentation")
    endif (DOXYGEN_FOUND)
endfunction(doxygen_doc)
