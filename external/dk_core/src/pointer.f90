module pointer
    implicit none(external, type)
    public

    interface
        function c_ptr_add(pointer, offset) bind(C, name="c_ptr_add")
            use iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr) :: c_ptr_add
            type(c_ptr), value :: pointer
            integer(c_size_t), value :: offset
        end function
    end interface
end module
