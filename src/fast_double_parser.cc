#include "fast_double_parser.h"
#include <cstdio>

extern "C" {
    double fast_double_strtod(const char *p, char **end);
}

// source file
double fast_double_strtod(const char *p, char **end) {
    double outDouble;
    // Use your existing implementation to parse the number.
    const char* endPtr = fast_double_parser::parse_number(p, &outDouble);
    // Set the output parameter `end` with the pointer returned by parse_number.
    *end = (char*)endPtr;  // Be careful with type casting. Make sure it's safe in your context.
    return outDouble;
}
