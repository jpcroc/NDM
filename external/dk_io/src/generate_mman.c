#include <sys/mman.h>
#include <errno.h>
#include <stdio.h>

int main() {
    printf("module mman\n");
    printf("    use iso_c_binding, only: c_int\n");
    printf("    implicit none(external, type)\n");
    printf("    integer(c_int), parameter :: PROT_EXEC = %i\n", PROT_EXEC);
    printf("    integer(c_int), parameter :: PROT_READ = %i\n", PROT_READ);
    printf("    integer(c_int), parameter :: PROT_WRITE = %i\n", PROT_WRITE);
    printf("    integer(c_int), parameter :: PROT_NONE = %i\n", PROT_NONE);
    printf("    integer(c_int), parameter :: MAP_SHARED = %i\n", MAP_SHARED);
//    printf("    integer(c_int), parameter :: MAP_SHARED_VALIDATE = %i\n", MAP_SHARED_VALIDATE);
    printf("    integer(c_int), parameter :: MAP_PRIVATE = %i\n", MAP_PRIVATE);
//    printf("    integer(c_int), parameter :: MAP_32BIT = %i\n", MAP_32BIT);
    printf("    integer(c_int), parameter :: MAP_ANON = %i\n", MAP_ANON);
    printf("    integer(c_int), parameter :: MAP_ANONYMOUS = %i\n", MAP_ANONYMOUS);
//    printf("    integer(c_int), parameter :: MAP_DENYWRITE = %i\n", MAP_DENYWRITE);
//    printf("    integer(c_int), parameter :: MAP_EXECUTABLE = %i\n", MAP_EXECUTABLE);
    printf("    integer(c_int), parameter :: MAP_FILE = %i\n", MAP_FILE);
    printf("    integer(c_int), parameter :: MAP_FIXED = %i\n", MAP_FIXED);
//    printf("    integer(c_int), parameter :: MAP_FIXED_NOREPLACE = %i\n", MAP_FIXED_NOREPLACE);
//    printf("    integer(c_int), parameter :: MAP_GROWSDOWN = %i\n", MAP_GROWSDOWN);
//    printf("    integer(c_int), parameter :: MAP_HUGETLB = %i\n", MAP_HUGETLB);
//    printf("    integer(c_int), parameter :: MAP_HUGE_2MB = %i\n", MAP_HUGE_2MB);
//    printf("    integer(c_int), parameter :: MAP_HUGE_1GB = %i\n", MAP_HUGE_1GB);
//    printf("    integer(c_int), parameter :: MAP_LOCKED = %i\n", MAP_LOCKED);
//    printf("    integer(c_int), parameter :: MAP_NONBLOCK = %i\n", MAP_NONBLOCK);
    printf("    integer(c_int), parameter :: MAP_NORESERVE = %i\n", MAP_NORESERVE);
//    printf("    integer(c_int), parameter :: MAP_POPULATE = %i\n", MAP_POPULATE);
//    printf("    integer(c_int), parameter :: MAP_SYNC = %i\n", MAP_SYNC);
//    printf("    integer(c_int), parameter :: MAP_UNINITIALIZED = %i\n", MAP_UNINITIALIZED);
    printf("    integer(c_int), parameter :: EACCES = %i\n", EACCES);
    printf("    integer(c_int), parameter :: EAGAIN = %i\n", EAGAIN);
    printf("    integer(c_int), parameter :: EBADF = %i\n", EBADF);
    printf("    integer(c_int), parameter :: EEXIST = %i\n", EEXIST);
    printf("    integer(c_int), parameter :: EINVAL = %i\n", EINVAL);
    printf("    integer(c_int), parameter :: ENFILE = %i\n", ENFILE);
    printf("    integer(c_int), parameter :: ENODEV = %i\n", ENODEV);
    printf("    integer(c_int), parameter :: ENOMEM = %i\n", ENOMEM);
    printf("    integer(c_int), parameter :: EOVERFLOW = %i\n", EOVERFLOW);
    printf("    integer(c_int), parameter :: EPERM = %i\n", EPERM);
    printf("    integer(c_int), parameter :: ETXTBSY = %i\n", ETXTBSY);
    printf("end module\n");
    return 0;
}
