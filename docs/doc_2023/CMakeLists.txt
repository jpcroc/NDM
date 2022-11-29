
log_info( "JUPYTER-BOOK compilation subdirectory docs/doc_2023" )

set(JUPY_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
# useless as build in source localize html and pdf in source.
set(JUPY_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/doc_2023_jupyter)

add_custom_target( doc_2023_html ALL
  COMMAND ${NDM_JUPYTER_BOOK} build ${JUPY_SOURCE}
  WORKING_DIRECTORY ${JUPY_SOURCE}
  COMMENT "Generating docs_2023 html documentation with jupyter-book"
)

add_custom_target( doc_2023_pdf ALL
  COMMAND ${NDM_JUPYTER_BOOK} build ${JUPY_SOURCE} --builder pdflatex
  WORKING_DIRECTORY ${JUPY_SOURCE}
  COMMENT "Generating docs_2023 pdf documentation with jupyter-book"
)

add_dependencies( doc_2023_html ndm_main.exe )
add_dependencies( doc_2023_pdf doc_2023_html )

install( FILES ${JUPY_SOURCE}/_build/latex/NDM.pdf DESTINATION ${JUPY_SOURCE}/.. )
install( FILES ${JUPY_SOURCE}/_build/latex/NDM.pdf DESTINATION ${CMAKE_INSTALL_PREFIX}/docs )

# printCmakeTrace()
