# look for mpirun location, compiler wrapper will be searched there first
macro(look_for_mpi project_name)
    # Find mpirun directory to look for wrapper
    set(${project_name}_MPIEXEC mpiexec mpirun)
    set(VAR_MPIEXEC "VAR_MPIEXEC-NOTFOUND" CACHE STRING "")
    foreach(mpi_prog IN ITEMS ${${project_name}_MPIEXEC})
        if (NOT VAR_MPIEXEC)
            if (DEFINED ENV{OMP_ROOT})
                find_program(VAR_MPIEXEC ${mpi_prog} PATHS $ENV{OMP_ROOT} NO_DEFAULT_PATH)
            elseif(DEFINED ENV{MPI_HOME})
	        find_program(VAR_MPIEXEC ${mpi_prog} PATHS $ENV{MPI_HOME} NO_DEFAULT_PATH)
		find_program(VAR_MPIEXEC ${mpi_prog} PATHS $ENV{MPI_HOME})
#                find_program(VAR_MPIEXEC ${mpi_prog} PATHS $ENV{MPI_HOME}/bin NO_DEFAULT_PATH)
            else()
                find_program(VAR_MPIEXEC ${mpi_prog})	# search anywhere
            endif()
        endif()
    endforeach()
    if (NOT VAR_MPIEXEC)
        message(STATUS "No mpirun found !")
    else()
        cmake_path(GET VAR_MPIEXEC PARENT_PATH mpi_path)
    endif()
endmacro(look_for_mpi)

# set CMAKE_<LANG>_COMPILER
macro(set_compiler project_name env_name lang)

    if (DEFINED ENV{OMPI_${env_name}})
        set(CMAKE_${lang}_COMPILER "$ENV{OMPI_${env_name}}" CACHE STRING "")
        message(STATUS "OMPI_${env_name} = $ENV{OMPI_${env_name}}")
    elseif(DEFINED ENV{MPI_${env_name}})
        set(CMAKE_${lang}_COMPILER "$ENV{MPI_${env_name}}" CACHE STRING "")
        message(STATUS "MPI_${env_name} = $ENV{MPI_${env_name}}")
    elseif(DEFINED ENV{${env_name}})
        set(CMAKE_${lang}_COMPILER "$ENV{${env_name}}" CACHE STRING "")
        message(STATUS "${env_name} = $ENV{${env_name}}")
    endif()

    # force <lang> compiler with CMAKE_<lang>_COMPILER
    if (DEFINED CMAKE_${lang}_COMPILER)
        foreach (config IN ITEMS ${compile_configuration_list})
            set(VAR_${config}_${lang} ${CMAKE_${lang}_COMPILER} CACHE STRING "" FORCE)
            #message(WARNING "CMAKE_${lang}_COMPILER = ${CMAKE_${lang}_COMPILER}")
        endforeach()
    endif()

    # force configuration with <project_name>_FORCE_CONFIG
    if (DEFINED ${project_name}_FORCE_CONFIG)
        if(NOT ${${project_name}_FORCE_CONFIG} IN_LIST compile_configuration_list)
            message(FATAL_ERROR "unknown configuration ${project_name}_FORCE_CONFIG: ${${project_name}_FORCE_CONFIG}\n
            available configurations are: ${compile_configuration_list}")
        else()
        foreach (config IN ITEMS ${compile_configuration_list})
            if (NOT ${config} STREQUAL ${${project_name}_FORCE_CONFIG})
                set(VAR_${config}_${lang} "NONE_${${project_name}_FORCE_CONFIG}_FORCED" CACHE STRING "" FORCE)
            elseif(DEFINED CMAKE_${lang}_COMPILER)
                cmake_path(GET CMAKE_${lang}_COMPILER FILENAME compiler_name)
                list(PREPEND ${config}_${lang}_COMPILERS ${compiler_name})
            endif()
        endforeach()
        endif()
    endif()
    
    look_for_mpi(${project_name})
    
    # For each possible configuration, look for a wrapper in the mpirun directory
    set(WRAPPER_FOUND False)
    foreach (config IN ITEMS ${compile_configuration_list})
        set(VAR_${config}_${lang} "VAR_${config}_${lang}-NOTFOUND" CACHE STRING "")
        foreach (wrapper IN ITEMS ${${config}_${lang}_WRAPPERS})
            if (NOT VAR_${config}_${lang})
                find_program(VAR_${config}_${lang} ${wrapper} PATHS ${mpi_path} NO_DEFAULT_PATH)
            endif()
        endforeach()
        if ((VAR_${config}_${lang}) AND (NOT VAR_${config}_${lang} STREQUAL "NONE_${${project_name}_FORCE_CONFIG}_FORCED")) # if a wrapper is found, get the underlying compiler        
	    set(WRAPPER_FOUND True)
	    execute_process(
                COMMAND ${VAR_${config}_${lang}} -show
                OUTPUT_VARIABLE SHOW_UNDERLYING ERROR_QUIET)
            if (SHOW_UNDERLYING)
                string(REPLACE " " ";" SHOW_UNDERLYING2 ${SHOW_UNDERLYING})
                list(GET SHOW_UNDERLYING2 0 ${config}_underlying_compiler)
            endif()
        endif()
    endforeach()

    # If no wrapper was found, look wider
    if (NOT ${WRAPPER_FOUND})
	    foreach (config IN ITEMS ${compile_configuration_list})
	        foreach (compiler IN ITEMS ${${config}_${lang}_COMPILERS})
		    if (NOT VAR_${config}_${lang})
		        find_program(VAR_${config}_${lang} ${compiler})
		    endif()
	        endforeach()
	        if (VAR_${config}_${lang}) # if a compiler is found, maybe it's a wrapper and we can get the underlying compiler
		    execute_process(
		        COMMAND ${VAR_${config}_${lang}} -show
		        OUTPUT_VARIABLE SHOW_UNDERLYING ERROR_QUIET)
		    if (SHOW_UNDERLYING)
		        string(REPLACE " " ";" SHOW_UNDERLYING2 ${SHOW_UNDERLYING})
		        list(GET SHOW_UNDERLYING2 0 ${config}_underlying_compiler)
		    endif()
	        endif()
	    endforeach()
    endif()
endmacro(set_compiler)

macro(set_config project_name lang)
    # Choose previously defined configuration
    set(config_set OFF)
    foreach (config IN ITEMS ${compile_configuration_list})
        if (NOT ${config_set} AND (((DEFINED ${project_name}_FORCE_CONFIG) AND (${config} STREQUAL ${${project_name}_FORCE_CONFIG})) OR (NOT (DEFINED ${project_name}_FORCE_CONFIG))))
            # if the underlying compiler is defined it means we have found a working wrapper for this config	
            if (DEFINED ${config}_underlying_compiler)
                if (${${config}_underlying_compiler} IN_LIST ${config}_${lang}_COMPILERS)
                    set(CMAKE_${lang}_COMPILER "${VAR_${config}_${lang}}" CACHE STRING "")
                    cmake_language(CALL CONFIG_${config})        #cmake 3.18
                    set(config_set ON)
                endif()
            elseif(DEFINED CMAKE_${lang}_COMPILER)
                cmake_path(GET CMAKE_${lang}_COMPILER FILENAME compiler_name)
                if (${compiler_name} IN_LIST ${config}_${lang}_COMPILERS)
                    set(CMAKE_${lang}_COMPILER "${VAR_${config}_${lang}}" CACHE STRING "")
                    cmake_language(CALL CONFIG_${config})        #cmake 3.18
                    set(config_set ON)
                endif()
            endif()
        endif()
    endforeach()
    
    # if still no config found, try serial
    if (NOT config_set)
        foreach (config IN ITEMS ${compile_configuration_list})
            if (VAR_${config}_${lang})
                set(CMAKE_${lang}_COMPILER "${VAR_${config}_${lang}}" CACHE STRING "")
		        cmake_language(CALL CONFIG_${config})        #cmake 3.18
                set(config_set ON)
            endif()
        endforeach()
    endif()

    # Display error/warning if no compiler was found
    if (NOT DEFINED CMAKE_${lang}_COMPILER)
        foreach (config IN ITEMS ${compile_configuration_list})
            if (VAR_${config}_${lang})
                message(STATUS "Compiler found for ${config}: ${VAR_${config}_${lang}}")
                if (DEFINED ${config}_underlying_compiler)
                    message(STATUS "Underlying compiler for ${VAR_${config}_${lang}}: ${${config}_underlying_compiler}")
                endif()
            endif()
        endforeach()
        message(FATAL_ERROR "No suitable ${lang} compiler found or provided, try to force a compiler with CMAKE_${lang}_COMPILER or a configuration with ${project_name}_FORCE_CONFIG")
    elseif(NOT ${config_set})
        message(WARNING "Unknown ${lang} compiler: ${CMAKE_${lang}_COMPILER}")
    endif()
endmacro(set_config)
