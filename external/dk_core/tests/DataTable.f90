program test_datatable
    use dk_datatable, only: DataTable
    use dk_datacolumn, only: DataColumn
    use ext_character, only: operator(//)
    use dk_test

    implicit none(external, type)

    type(DataTable), target :: table
    type(DataColumn), pointer :: col
    character(:), allocatable :: header, name, colname

    call check(table%num_columns()==0, "number of columns in an empty table", &
        "Found "//table%num_columns()//"instead of 0")
    call table%add_column("col_1")
    call check(table%num_columns()==1, "number of columns after adding a column", &
        "Found "//table%num_columns()//"instead of 1")
    call table%add_column("col_2")
    call check(table%num_columns()==2, "number of columns after adding two columns", &
        "Found "//table%num_columns()//"instead of 2")
    call table%get_header(header)

    call check(allocated(header), "header allocation status")
    if (allocated(header)) then
        call check(header == " col_1   col_2", "header value", &
        "Header: '"//header//"' instead of ' col_1   col_2'")
    end if

    name = table%column_name(1)
    call check(allocated(name), "column 1 name allocation status")
    if (allocated(name)) then
        call check(name == "col_1", "column 1 name value", &
        "Name: '"//name//"' instead of 'col_1'")
    end if

    name = table%column_name(2)
    call check(allocated(name), "column 2 name allocation status")
    if (allocated(name)) then
        call check(name == "col_2", "column 2 name value", &
        "Name: '"//name//"' instead of 'col_2'")
    end if

    call tests_finished
end program
