#include <stdio.h>
#include <mpi.h>

main(argc, argv)
int argc;
char *argv[];
{
        char name[BUFSIZ];
        int length;

        MPI_Init(&argc, &argv);
        MPI_Get_processor_name(name, &length);
        printf("%s: hello world\n", name);
        MPI_Finalize();
}
