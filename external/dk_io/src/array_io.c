#include "grisu3.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double gay_strtod(const char* s00, char** se);
double fast_double_strtod(const char* p, char** end);
int s2d(const char* buffer, double* result);

// str must be null-terminated, and both str and work must be of the same size
// n is the size of x
// pos is the position of the character string that failed to be read as a number.
int dk_read_array(const char* str, double* x, size_t n, size_t* pos)
{
    if (str == NULL || x == NULL || n == 0) {
        *pos = 0;
        return 0;
    }

    size_t i, j;
    char* end;
    const char* start;

    start = str;

    for (i = 0; i < n; i++) {

        // Ignore spaces
        while (start[0] == 32) {
            start++;
        }

        // x[i] = strtod(start, &end);
        x[i] = gay_strtod(start, &end);

        // x[i] = fast_double_strtod(start, &end);

        // Return if the value could not be read
        if (start == end) {
            for (j = i; j < n; j++) {
                x[j] = NAN;
            }
            *pos = start - str + 1;
            return i;
        }

        start = end;

        // Ignore spaces
        while (start[0] == 32) {
            start++;
        }

        // Return if we found the end of line
        if (start[0] == 0 || start[0] == 10) {
            for (j = i + 1; j < n; j++) {
                x[j] = NAN;
            }
            *pos = start - str + 1;
            return i + 1;
        }
    }

    // Return the number of values that could be read
    *pos = start - str + 1;
    return i;
}

// char *dtoa_r(double dd, int mode, int ndigits, int *decpt, int *sign, char **rve, char *buf, size_t blen);
// int dtoa_grisu3(double v, char *dst);
int dtoa_grisu3(double v, char* dst);
char* gay_dtoa(double d, int mode, int ndigits, int* decpt, int* sign, char** rve);

int dk_write_array(char* str, double* x, size_t n, size_t buffsize)
{
    int i;
    int l;
    size_t s;

    // Fill the buffer with spaces
    memset(str, ' ', buffsize);

    l = 0;

    for (i = 0; i < n; i++) {

        s = dtoa_grisu3(x[i], str);

        str[s] = ' ';
        str = str + s + 1;
        l = l + s + 1;
    }

    // Return the number of elements subccessfully written
    return l;
}
