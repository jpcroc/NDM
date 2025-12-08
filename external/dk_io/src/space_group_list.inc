    type(SpaceGroup), dimension(*), parameter, private :: SPACE_GROUPS = [&
        SpaceGroup(1, "C1^1  ", "P 1             ", "P 1                            ", &
     "P 1                ", "P1        ", "     ", PRIMITIVE, 1), & !/*   1 */
        SpaceGroup(2, "Ci^1  ", "-P 1            ", "P -1                           ", &
     "P -1               ", "P-1       ", "     ", PRIMITIVE, 2), & !/*   2 */
        SpaceGroup(3, "C2^1  ", "P 2y            ", "P 2 = P 1 2 1                  ", &
     "P 1 2 1            ", "P2        ", "b    ", PRIMITIVE, 3), & !/*   3 */
        SpaceGroup(3, "C2^1  ", "P 2             ", "P 2 = P 1 1 2                  ", &
     "P 1 1 2            ", "P2        ", "c    ", PRIMITIVE, 3), & !/*   4 */
        SpaceGroup(3, "C2^1  ", "P 2x            ", "P 2 = P 2 1 1                  ", &
     "P 2 1 1            ", "P2        ", "a    ", PRIMITIVE, 3), & !/*   5 */
        SpaceGroup(4, "C2^2  ", "P 2yb           ", "P 2_1 = P 1 2_1 1              ", &
     "P 1 2_1 1          ", "P2_1      ", "b    ", PRIMITIVE, 3), & !/*   6 */
        SpaceGroup(4, "C2^2  ", "P 2c            ", "P 2_1 = P 1 1 2_1              ", &
     "P 1 1 2_1          ", "P2_1      ", "c    ", PRIMITIVE, 3), & !/*   7 */
        SpaceGroup(4, "C2^2  ", "P 2xa           ", "P 2_1 = P 2_1 1 1              ", &
     "P 2_1 1 1          ", "P2_1      ", "a    ", PRIMITIVE, 3), & !/*   8 */
        SpaceGroup(5, "C2^3  ", "C 2y            ", "C 2 = C 1 2 1                  ", &
     "C 1 2 1            ", "C2        ", "b1   ", C_FACE, 3), & !/*   9 */
        SpaceGroup(5, "C2^3  ", "A 2y            ", "C 2 = A 1 2 1                  ", &
     "A 1 2 1            ", "C2        ", "b2   ", A_FACE, 3), & !/*  10 */
        SpaceGroup(5, "C2^3  ", "I 2y            ", "C 2 = I 1 2 1                  ", &
     "I 1 2 1            ", "C2        ", "b3   ", BODY, 3), & !/*  11 */
        SpaceGroup(5, "C2^3  ", "A 2             ", "C 2 = A 1 1 2                  ", &
     "A 1 1 2            ", "C2        ", "c1   ", A_FACE, 3), & !/*  12 */
        SpaceGroup(5, "C2^3  ", "B 2             ", "C 2 = B 1 1 2 = B 2            ", &
     "B 1 1 2            ", "C2        ", "c2   ", PRIMITIVE, 3), & !/*  13 */
        SpaceGroup(5, "C2^3  ", "I 2             ", "C 2 = I 1 1 2                  ", &
     "I 1 1 2            ", "C2        ", "c3   ", BODY, 3), & !/*  14 */
        SpaceGroup(5, "C2^3  ", "B 2x            ", "C 2 = B 2 1 1                  ", &
     "B 2 1 1            ", "C2        ", "a1   ", PRIMITIVE, 3), & !/*  15 */
        SpaceGroup(5, "C2^3  ", "C 2x            ", "C 2 = C 2 1 1                  ", &
     "C 2 1 1            ", "C2        ", "a2   ", C_FACE, 3), & !/*  16 */
        SpaceGroup(5, "C2^3  ", "I 2x            ", "C 2 = I 2 1 1                  ", &
     "I 2 1 1            ", "C2        ", "a3   ", BODY, 3), & !/*  17 */
        SpaceGroup(6, "Cs^1  ", "P -2y           ", "P m = P 1 m 1                  ", &
     "P 1 m 1            ", "Pm        ", "b    ", PRIMITIVE, 4), & !/*  18 */
        SpaceGroup(6, "Cs^1  ", "P -2            ", "P m = P 1 1 m                  ", &
     "P 1 1 m            ", "Pm        ", "c    ", PRIMITIVE, 4), & !/*  19 */
        SpaceGroup(6, "Cs^1  ", "P -2x           ", "P m = P m 1 1                  ", &
     "P m 1 1            ", "Pm        ", "a    ", PRIMITIVE, 4), & !/*  20 */
        SpaceGroup(7, "Cs^2  ", "P -2yc          ", "P c = P 1 c 1                  ", &
     "P 1 c 1            ", "Pc        ", "b1   ", PRIMITIVE, 4), & !/*  21 */
        SpaceGroup(7, "Cs^2  ", "P -2yac         ", "P c = P 1 n 1                  ", &
     "P 1 n 1            ", "Pc        ", "b2   ", PRIMITIVE, 4), & !/*  22 */
        SpaceGroup(7, "Cs^2  ", "P -2ya          ", "P c = P 1 a 1                  ", &
     "P 1 a 1            ", "Pc        ", "b3   ", PRIMITIVE, 4), & !/*  23 */
        SpaceGroup(7, "Cs^2  ", "P -2a           ", "P c = P 1 1 a                  ", &
     "P 1 1 a            ", "Pc        ", "c1   ", PRIMITIVE, 4), & !/*  24 */
        SpaceGroup(7, "Cs^2  ", "P -2ab          ", "P c = P 1 1 n                  ", &
     "P 1 1 n            ", "Pc        ", "c2   ", PRIMITIVE, 4), & !/*  25 */
        SpaceGroup(7, "Cs^2  ", "P -2b           ", "P c = P 1 1 b = P b            ", &
     "P 1 1 b            ", "Pc        ", "c3   ", PRIMITIVE, 4), & !/*  26 */
        SpaceGroup(7, "Cs^2  ", "P -2xb          ", "P c = P b 1 1                  ", &
     "P b 1 1            ", "Pc        ", "a1   ", PRIMITIVE, 4), & !/*  27 */
        SpaceGroup(7, "Cs^2  ", "P -2xbc         ", "P c = P n 1 1                  ", &
     "P n 1 1            ", "Pc        ", "a2   ", PRIMITIVE, 4), & !/*  28 */
        SpaceGroup(7, "Cs^2  ", "P -2xc          ", "P c = P c 1 1                  ", &
     "P c 1 1            ", "Pc        ", "a3   ", PRIMITIVE, 4), & !/*  29 */
        SpaceGroup(8, "Cs^3  ", "C -2y           ", "C m = C 1 m 1                  ", &
     "C 1 m 1            ", "Cm        ", "b1   ", C_FACE, 4), & !/*  30 */
        SpaceGroup(8, "Cs^3  ", "A -2y           ", "C m = A 1 m 1                  ", &
     "A 1 m 1            ", "Cm        ", "b2   ", A_FACE, 4), & !/*  31 */
        SpaceGroup(8, "Cs^3  ", "I -2y           ", "C m = I 1 m 1                  ", &
     "I 1 m 1            ", "Cm        ", "b3   ", BODY, 4), & !/*  32 */
        SpaceGroup(8, "Cs^3  ", "A -2            ", "C m = A 1 1 m                  ", &
     "A 1 1 m            ", "Cm        ", "c1   ", A_FACE, 4), & !/*  33 */
        SpaceGroup(8, "Cs^3  ", "B -2            ", "C m = B 1 1 m = B m            ", &
     "B 1 1 m            ", "Cm        ", "c2   ", PRIMITIVE, 4), & !/*  34 */
        SpaceGroup(8, "Cs^3  ", "I -2            ", "C m = I 1 1 m                  ", &
     "I 1 1 m            ", "Cm        ", "c3   ", BODY, 4), & !/*  35 */
        SpaceGroup(8, "Cs^3  ", "B -2x           ", "C m = B m 1 1                  ", &
     "B m 1 1            ", "Cm        ", "a1   ", PRIMITIVE, 4), & !/*  36 */
        SpaceGroup(8, "Cs^3  ", "C -2x           ", "C m = C m 1 1                  ", &
     "C m 1 1            ", "Cm        ", "a2   ", C_FACE, 4), & !/*  37 */
        SpaceGroup(8, "Cs^3  ", "I -2x           ", "C m = I m 1 1                  ", &
     "I m 1 1            ", "Cm        ", "a3   ", BODY, 4), & !/*  38 */
        SpaceGroup(9, "Cs^4  ", "C -2yc          ", "C c = C 1 c 1                  ", &
     "C 1 c 1            ", "Cc        ", "b1   ", C_FACE, 4), & !/*  39 */
        SpaceGroup(9, "Cs^4  ", "A -2yab         ", "C c = A 1 n 1                  ", &
     "A 1 n 1            ", "Cc        ", "b2   ", A_FACE, 4), & !/*  40 */
        SpaceGroup(9, "Cs^4  ", "I -2ya          ", "C c = I 1 a 1                  ", &
     "I 1 a 1            ", "Cc        ", "b3   ", BODY, 4), & !/*  41 */
        SpaceGroup(9, "Cs^4  ", "A -2ya          ", "C c = A 1 a 1                  ", &
     "A 1 a 1            ", "Cc        ", "-b1  ", A_FACE, 4), & !/*  42 */
        SpaceGroup(9, "Cs^4  ", "C -2yac         ", "C c = C 1 n 1                  ", &
     "C 1 n 1            ", "Cc        ", "-b2  ", C_FACE, 4), & !/*  43 */
        SpaceGroup(9, "Cs^4  ", "I -2yc          ", "C c = I 1 c 1                  ", &
     "I 1 c 1            ", "Cc        ", "-b3  ", BODY, 4), & !/*  44 */
        SpaceGroup(9, "Cs^4  ", "A -2a           ", "C c = A 1 1 a                  ", &
     "A 1 1 a            ", "Cc        ", "c1   ", A_FACE, 4), & !/*  45 */
        SpaceGroup(9, "Cs^4  ", "B -2ab          ", "C c = B 1 1 n                  ", &
     "B 1 1 n            ", "Cc        ", "c2   ", PRIMITIVE, 4), & !/*  46 */
        SpaceGroup(9, "Cs^4  ", "I -2b           ", "C c = I 1 1 b                  ", &
     "I 1 1 b            ", "Cc        ", "c3   ", BODY, 4), & !/*  47 */
        SpaceGroup(9, "Cs^4  ", "B -2b           ", "C c = B 1 1 b = B b            ", &
     "B 1 1 b            ", "Cc        ", "-c1  ", PRIMITIVE, 4), & !/*  48 */
        SpaceGroup(9, "Cs^4  ", "A -2ab          ", "C c = A 1 1 n                  ", &
     "A 1 1 n            ", "Cc        ", "-c2  ", A_FACE, 4), & !/*  49 */
        SpaceGroup(9, "Cs^4  ", "I -2a           ", "C c = I 1 1 a                  ", &
     "I 1 1 a            ", "Cc        ", "-c3  ", BODY, 4), & !/*  50 */
        SpaceGroup(9, "Cs^4  ", "B -2xb          ", "C c = B b 1 1                  ", &
     "B b 1 1            ", "Cc        ", "a1   ", PRIMITIVE, 4), & !/*  51 */
        SpaceGroup(9, "Cs^4  ", "C -2xac         ", "C c = C n 1 1                  ", &
     "C n 1 1            ", "Cc        ", "a2   ", C_FACE, 4), & !/*  52 */
        SpaceGroup(9, "Cs^4  ", "I -2xc          ", "C c = I c 1 1                  ", &
     "I c 1 1            ", "Cc        ", "a3   ", BODY, 4), & !/*  53 */
        SpaceGroup(9, "Cs^4  ", "C -2xc          ", "C c = C c 1 1                  ", &
     "C c 1 1            ", "Cc        ", "-a1  ", C_FACE, 4), & !/*  54 */
        SpaceGroup(9, "Cs^4  ", "B -2xab         ", "C c = B n 1 1                  ", &
     "B n 1 1            ", "Cc        ", "-a2  ", PRIMITIVE, 4), & !/*  55 */
        SpaceGroup(9, "Cs^4  ", "I -2xb          ", "C c = I b 1 1                  ", &
     "I b 1 1            ", "Cc        ", "-a3  ", BODY, 4), & !/*  56 */
        SpaceGroup(10, "C2h^1 ", "-P 2y           ", "P 2/m = P 1 2/m 1              ", &
     "P 1 2/m 1          ", "P2/m      ", "b    ", PRIMITIVE, 5), & !/*  57 */
        SpaceGroup(10, "C2h^1 ", "-P 2            ", "P 2/m = P 1 1 2/m              ", &
     "P 1 1 2/m          ", "P2/m      ", "c    ", PRIMITIVE, 5), & !/*  58 */
        SpaceGroup(10, "C2h^1 ", "-P 2x           ", "P 2/m = P 2/m 1 1              ", &
     "P 2/m 1 1          ", "P2/m      ", "a    ", PRIMITIVE, 5), & !/*  59 */
        SpaceGroup(11, "C2h^2 ", "-P 2yb          ", "P 2_1/m = P 1 2_1/m 1          ", &
     "P 1 2_1/m 1        ", "P2_1/m    ", "b    ", PRIMITIVE, 5), & !/*  60 */
        SpaceGroup(11, "C2h^2 ", "-P 2c           ", "P 2_1/m = P 1 1 2_1/m          ", &
     "P 1 1 2_1/m        ", "P2_1/m    ", "c    ", PRIMITIVE, 5), & !/*  61 */
        SpaceGroup(11, "C2h^2 ", "-P 2xa          ", "P 2_1/m = P 2_1/m 1 1          ", &
     "P 2_1/m 1 1        ", "P2_1/m    ", "a    ", PRIMITIVE, 5), & !/*  62 */
        SpaceGroup(12, "C2h^3 ", "-C 2y           ", "C 2/m = C 1 2/m 1              ", &
     "C 1 2/m 1          ", "C2/m      ", "b1   ", C_FACE, 5), & !/*  63 */
        SpaceGroup(12, "C2h^3 ", "-A 2y           ", "C 2/m = A 1 2/m 1              ", &
     "A 1 2/m 1          ", "C2/m      ", "b2   ", A_FACE, 5), & !/*  64 */
        SpaceGroup(12, "C2h^3 ", "-I 2y           ", "C 2/m = I 1 2/m 1              ", &
     "I 1 2/m 1          ", "C2/m      ", "b3   ", BODY, 5), & !/*  65 */
        SpaceGroup(12, "C2h^3 ", "-A 2            ", "C 2/m = A 1 1 2/m              ", &
     "A 1 1 2/m          ", "C2/m      ", "c1   ", A_FACE, 5), & !/*  66 */
        SpaceGroup(12, "C2h^3 ", "-B 2            ", "C 2/m = B 1 1 2/m = B 2/m      ", &
     "B 1 1 2/m          ", "C2/m      ", "c2   ", PRIMITIVE, 5), & !/*  67 */
        SpaceGroup(12, "C2h^3 ", "-I 2            ", "C 2/m = I 1 1 2/m              ", &
     "I 1 1 2/m          ", "C2/m      ", "c3   ", BODY, 5), & !/*  68 */
        SpaceGroup(12, "C2h^3 ", "-B 2x           ", "C 2/m = B 2/m 1 1              ", &
     "B 2/m 1 1          ", "C2/m      ", "a1   ", PRIMITIVE, 5), & !/*  69 */
        SpaceGroup(12, "C2h^3 ", "-C 2x           ", "C 2/m = C 2/m 1 1              ", &
     "C 2/m 1 1          ", "C2/m      ", "a2   ", C_FACE, 5), & !/*  70 */
        SpaceGroup(12, "C2h^3 ", "-I 2x           ", "C 2/m = I 2/m 1 1              ", &
     "I 2/m 1 1          ", "C2/m      ", "a3   ", BODY, 5), & !/*  71 */
        SpaceGroup(13, "C2h^4 ", "-P 2yc          ", "P 2/c = P 1 2/c 1              ", &
     "P 1 2/c 1          ", "P2/c      ", "b1   ", PRIMITIVE, 5), & !/*  72 */
        SpaceGroup(13, "C2h^4 ", "-P 2yac         ", "P 2/c = P 1 2/n 1              ", &
     "P 1 2/n 1          ", "P2/c      ", "b2   ", PRIMITIVE, 5), & !/*  73 */
        SpaceGroup(13, "C2h^4 ", "-P 2ya          ", "P 2/c = P 1 2/a 1              ", &
     "P 1 2/a 1          ", "P2/c      ", "b3   ", PRIMITIVE, 5), & !/*  74 */
        SpaceGroup(13, "C2h^4 ", "-P 2a           ", "P 2/c = P 1 1 2/a              ", &
     "P 1 1 2/a          ", "P2/c      ", "c1   ", PRIMITIVE, 5), & !/*  75 */
        SpaceGroup(13, "C2h^4 ", "-P 2ab          ", "P 2/c = P 1 1 2/n              ", &
     "P 1 1 2/n          ", "P2/c      ", "c2   ", PRIMITIVE, 5), & !/*  76 */
        SpaceGroup(13, "C2h^4 ", "-P 2b           ", "P 2/c = P 1 1 2/b = P 2/b      ", &
     "P 1 1 2/b          ", "P2/c      ", "c3   ", PRIMITIVE, 5), & !/*  77 */
        SpaceGroup(13, "C2h^4 ", "-P 2xb          ", "P 2/c = P 2/b 1 1              ", &
     "P 2/b 1 1          ", "P2/c      ", "a1   ", PRIMITIVE, 5), & !/*  78 */
        SpaceGroup(13, "C2h^4 ", "-P 2xbc         ", "P 2/c = P 2/n 1 1              ", &
     "P 2/n 1 1          ", "P2/c      ", "a2   ", PRIMITIVE, 5), & !/*  79 */
        SpaceGroup(13, "C2h^4 ", "-P 2xc          ", "P 2/c = P 2/c 1 1              ", &
     "P 2/c 1 1          ", "P2/c      ", "a3   ", PRIMITIVE, 5), & !/*  80 */
        SpaceGroup(14, "C2h^5 ", "-P 2ybc         ", "P 2_1/c = P 1 2_1/c 1          ", &
     "P 1 2_1/c 1        ", "P2_1/c    ", "b1   ", PRIMITIVE, 5), & !/*  81 */
        SpaceGroup(14, "C2h^5 ", "-P 2yn          ", "P 2_1/c = P 1 2_1/n 1          ", &
     "P 1 2_1/n 1        ", "P2_1/c    ", "b2   ", PRIMITIVE, 5), & !/*  82 */
        SpaceGroup(14, "C2h^5 ", "-P 2yab         ", "P 2_1/c = P 1 2_1/a 1          ", &
     "P 1 2_1/a 1        ", "P2_1/c    ", "b3   ", PRIMITIVE, 5), & !/*  83 */
        SpaceGroup(14, "C2h^5 ", "-P 2ac          ", "P 2_1/c = P 1 1 2_1/a          ", &
     "P 1 1 2_1/a        ", "P2_1/c    ", "c1   ", PRIMITIVE, 5), & !/*  84 */
        SpaceGroup(14, "C2h^5 ", "-P 2n           ", "P 2_1/c = P 1 1 2_1/n          ", &
     "P 1 1 2_1/n        ", "P2_1/c    ", "c2   ", PRIMITIVE, 5), & !/*  85 */
        SpaceGroup(14, "C2h^5 ", "-P 2bc          ", "P 2_1/c = P 1 1 2_1/b = P 2_1/b", &
     "P 1 1 2_1/b        ", "P2_1/c    ", "c3   ", PRIMITIVE, 5), & !/*  86 */
        SpaceGroup(14, "C2h^5 ", "-P 2xab         ", "P 2_1/c = P 2_1/b 1 1          ", &
     "P 2_1/b 1 1        ", "P2_1/c    ", "a1   ", PRIMITIVE, 5), & !/*  87 */
        SpaceGroup(14, "C2h^5 ", "-P 2xn          ", "P 2_1/c = P 2_1/n 1 1          ", &
     "P 2_1/n 1 1        ", "P2_1/c    ", "a2   ", PRIMITIVE, 5), & !/*  88 */
        SpaceGroup(14, "C2h^5 ", "-P 2xac         ", "P 2_1/c = P 2_1/c 1 1          ", &
     "P 2_1/c 1 1        ", "P2_1/c    ", "a3   ", PRIMITIVE, 5), & !/*  89 */
        SpaceGroup(15, "C2h^6 ", "-C 2yc          ", "C 2/c = C 1 2/c 1              ", &
     "C 1 2/c 1          ", "C2/c      ", "b1   ", C_FACE, 5), & !/*  90 */
        SpaceGroup(15, "C2h^6 ", "-A 2yab         ", "C 2/c = A 1 2/n 1              ", &
     "A 1 2/n 1          ", "C2/c      ", "b2   ", A_FACE, 5), & !/*  91 */
        SpaceGroup(15, "C2h^6 ", "-I 2ya          ", "C 2/c = I 1 2/a 1              ", &
     "I 1 2/a 1          ", "C2/c      ", "b3   ", BODY, 5), & !/*  92 */
        SpaceGroup(15, "C2h^6 ", "-A 2ya          ", "C 2/c = A 1 2/a 1              ", &
     "A 1 2/a 1          ", "C2/c      ", "-b1  ", A_FACE, 5), & !/*  93 */
        SpaceGroup(15, "C2h^6 ", "-C 2yac         ", "C 2/c = C 1 2/n 1              ", &
     "C 1 2/n 1          ", "C2/c      ", "-b2  ", C_FACE, 5), & !/*  94 */
        SpaceGroup(15, "C2h^6 ", "-I 2yc          ", "C 2/c = I 1 2/c 1              ", &
     "I 1 2/c 1          ", "C2/c      ", "-b3  ", BODY, 5), & !/*  95 */
        SpaceGroup(15, "C2h^6 ", "-A 2a           ", "C 2/c = A 1 1 2/a              ", &
     "A 1 1 2/a          ", "C2/c      ", "c1   ", A_FACE, 5), & !/*  96 */
        SpaceGroup(15, "C2h^6 ", "-B 2ab          ", "C 2/c = B 1 1 2/n              ", &
     "B 1 1 2/n          ", "C2/c      ", "c2   ", PRIMITIVE, 5), & !/*  97 */
        SpaceGroup(15, "C2h^6 ", "-I 2b           ", "C 2/c = I 1 1 2/b              ", &
     "I 1 1 2/b          ", "C2/c      ", "c3   ", BODY, 5), & !/*  98 */
        SpaceGroup(15, "C2h^6 ", "-B 2b           ", "C 2/c = B 1 1 2/b = B 2/b      ", &
     "B 1 1 2/b          ", "C2/c      ", "-c1  ", PRIMITIVE, 5), & !/*  99 */
        SpaceGroup(15, "C2h^6 ", "-A 2ab          ", "C 2/c = A 1 1 2/n              ", &
     "A 1 1 2/n          ", "C2/c      ", "-c2  ", A_FACE, 5), & !/* 100 */
        SpaceGroup(15, "C2h^6 ", "-I 2a           ", "C 2/c = I 1 1 2/a              ", &
     "I 1 1 2/a          ", "C2/c      ", "-c3  ", BODY, 5), & !/* 101 */
        SpaceGroup(15, "C2h^6 ", "-B 2xb          ", "C 2/c = B 2/b 1 1              ", &
     "B 2/b 1 1          ", "C2/c      ", "a1   ", PRIMITIVE, 5), & !/* 102 */
        SpaceGroup(15, "C2h^6 ", "-C 2xac         ", "C 2/c = C 2/n 1 1              ", &
     "C 2/n 1 1          ", "C2/c      ", "a2   ", C_FACE, 5), & !/* 103 */
        SpaceGroup(15, "C2h^6 ", "-I 2xc          ", "C 2/c = I 2/c 1 1              ", &
     "I 2/c 1 1          ", "C2/c      ", "a3   ", BODY, 5), & !/* 104 */
        SpaceGroup(15, "C2h^6 ", "-C 2xc          ", "C 2/c = C 2/c 1 1              ", &
     "C 2/c 1 1          ", "C2/c      ", "-a1  ", C_FACE, 5), & !/* 105 */
        SpaceGroup(15, "C2h^6 ", "-B 2xab         ", "C 2/c = B 2/n 1 1              ", &
     "B 2/n 1 1          ", "C2/c      ", "-a2  ", PRIMITIVE, 5), & !/* 106 */
        SpaceGroup(15, "C2h^6 ", "-I 2xb          ", "C 2/c = I 2/b 1 1              ", &
     "I 2/b 1 1          ", "C2/c      ", "-a3  ", BODY, 5), & !/* 107 */
        SpaceGroup(16, "D2^1  ", "P 2 2           ", "P 2 2 2                        ", &
     "P 2 2 2            ", "P222      ", "     ", PRIMITIVE, 6), & !/* 108 */
        SpaceGroup(17, "D2^2  ", "P 2c 2          ", "P 2 2 2_1                      ", &
     "P 2 2 2_1          ", "P222_1    ", "     ", PRIMITIVE, 6), & !/* 109 */
        SpaceGroup(17, "D2^2  ", "P 2a 2a         ", "P 2_1 2 2                      ", &
     "P 2_1 2 2          ", "P2_122    ", "cab  ", PRIMITIVE, 6), & !/* 110 */
        SpaceGroup(17, "D2^2  ", "P 2 2b          ", "P 2 2_1 2                      ", &
     "P 2 2_1 2          ", "P22_12    ", "bca  ", PRIMITIVE, 6), & !/* 111 */
        SpaceGroup(18, "D2^3  ", "P 2 2ab         ", "P 2_1 2_1 2                    ", &
     "P 2_1 2_1 2        ", "P2_12_12  ", "     ", PRIMITIVE, 6), & !/* 112 */
        SpaceGroup(18, "D2^3  ", "P 2bc 2         ", "P 2 2_1 2_1                    ", &
     "P 2 2_1 2_1        ", "P22_12_1  ", "cab  ", PRIMITIVE, 6), & !/* 113 */
        SpaceGroup(18, "D2^3  ", "P 2ac 2ac       ", "P 2_1 2 2_1                    ", &
     "P 2_1 2 2_1        ", "P2_122_1  ", "bca  ", PRIMITIVE, 6), & !/* 114 */
        SpaceGroup(19, "D2^4  ", "P 2ac 2ab       ", "P 2_1 2_1 2_1                  ", &
     "P 2_1 2_1 2_1      ", "P2_12_12_1", "     ", PRIMITIVE, 6), & !/* 115 */
        SpaceGroup(20, "D2^5  ", "C 2c 2          ", "C 2 2 2_1                      ", &
     "C 2 2 2_1          ", "C222_1    ", "     ", C_FACE, 6), & !/* 116 */
        SpaceGroup(20, "D2^5  ", "A 2a 2a         ", "A 2_1 2 2                      ", &
     "A 2_1 2 2          ", "A2_122    ", "cab  ", A_FACE, 6), & !/* 117 */
        SpaceGroup(20, "D2^5  ", "B 2 2b          ", "B 2 2_1 2                      ", &
     "B 2 2_1 2          ", "B22_12    ", "bca  ", PRIMITIVE, 6), & !/* 118 */
        SpaceGroup(21, "D2^6  ", "C 2 2           ", "C 2 2 2                        ", &
     "C 2 2 2            ", "C222      ", "     ", C_FACE, 6), & !/* 119 */
        SpaceGroup(21, "D2^6  ", "A 2 2           ", "A 2 2 2                        ", &
     "A 2 2 2            ", "A222      ", "cab  ", A_FACE, 6), & !/* 120 */
        SpaceGroup(21, "D2^6  ", "B 2 2           ", "B 2 2 2                        ", &
     "B 2 2 2            ", "B222      ", "bca  ", PRIMITIVE, 6), & !/* 121 */
        SpaceGroup(22, "D2^7  ", "F 2 2           ", "F 2 2 2                        ", &
     "F 2 2 2            ", "F222      ", "     ", FACE, 6), & !/* 122 */
        SpaceGroup(23, "D2^8  ", "I 2 2           ", "I 2 2 2                        ", &
     "I 2 2 2            ", "I222      ", "     ", BODY, 6), & !/* 123 */
        SpaceGroup(24, "D2^9  ", "I 2b 2c         ", "I 2_1 2_1 2_1                  ", &
     "I 2_1 2_1 2_1      ", "I2_12_12_1", "     ", BODY, 6), & !/* 124 */
        SpaceGroup(25, "C2v^1 ", "P 2 -2          ", "P m m 2                        ", &
     "P m m 2            ", "Pmm2      ", "     ", PRIMITIVE, 7), & !/* 125 */
        SpaceGroup(25, "C2v^1 ", "P -2 2          ", "P 2 m m                        ", &
     "P 2 m m            ", "P2mm      ", "cab  ", PRIMITIVE, 7), & !/* 126 */
        SpaceGroup(25, "C2v^1 ", "P -2 -2         ", "P m 2 m                        ", &
     "P m 2 m            ", "Pm2m      ", "bca  ", PRIMITIVE, 7), & !/* 127 */
        SpaceGroup(26, "C2v^2 ", "P 2c -2         ", "P m c 2_1                      ", &
     "P m c 2_1          ", "Pmc2_1    ", "     ", PRIMITIVE, 7), & !/* 128 */
        SpaceGroup(26, "C2v^2 ", "P 2c -2c        ", "P c m 2_1                      ", &
     "P c m 2_1          ", "Pcm2_1    ", "ba-c ", PRIMITIVE, 7), & !/* 129 */
        SpaceGroup(26, "C2v^2 ", "P -2a 2a        ", "P 2_1 m a                      ", &
     "P 2_1 m a          ", "P2_1ma    ", "cab  ", PRIMITIVE, 7), & !/* 130 */
        SpaceGroup(26, "C2v^2 ", "P -2 2a         ", "P 2_1 a m                      ", &
     "P 2_1 a m          ", "P2_1am    ", "-cba ", PRIMITIVE, 7), & !/* 131 */
        SpaceGroup(26, "C2v^2 ", "P -2 -2b        ", "P b 2_1 m                      ", &
     "P b 2_1 m          ", "Pb2_1m    ", "bca  ", PRIMITIVE, 7), & !/* 132 */
        SpaceGroup(26, "C2v^2 ", "P -2b -2        ", "P m 2_1 b                      ", &
     "P m 2_1 b          ", "Pm2_1b    ", "a-cb ", PRIMITIVE, 7), & !/* 133 */
        SpaceGroup(27, "C2v^3 ", "P 2 -2c         ", "P c c 2                        ", &
     "P c c 2            ", "Pcc2      ", "     ", PRIMITIVE, 7), & !/* 134 */
        SpaceGroup(27, "C2v^3 ", "P -2a 2         ", "P 2 a a                        ", &
     "P 2 a a            ", "P2aa      ", "cab  ", PRIMITIVE, 7), & !/* 135 */
        SpaceGroup(27, "C2v^3 ", "P -2b -2b       ", "P b 2 b                        ", &
     "P b 2 b            ", "Pb2b      ", "bca  ", PRIMITIVE, 7), & !/* 136 */
        SpaceGroup(28, "C2v^4 ", "P 2 -2a         ", "P m a 2                        ", &
     "P m a 2            ", "Pma2      ", "     ", PRIMITIVE, 7), & !/* 137 */
        SpaceGroup(28, "C2v^4 ", "P 2 -2b         ", "P b m 2                        ", &
     "P b m 2            ", "Pbm2      ", "ba-c ", PRIMITIVE, 7), & !/* 138 */
        SpaceGroup(28, "C2v^4 ", "P -2b 2         ", "P 2 m b                        ", &
     "P 2 m b            ", "P2mb      ", "cab  ", PRIMITIVE, 7), & !/* 139 */
        SpaceGroup(28, "C2v^4 ", "P -2c 2         ", "P 2 c m                        ", &
     "P 2 c m            ", "P2cm      ", "-cba ", PRIMITIVE, 7), & !/* 140 */
        SpaceGroup(28, "C2v^4 ", "P -2c -2c       ", "P c 2 m                        ", &
     "P c 2 m            ", "Pc2m      ", "bca  ", PRIMITIVE, 7), & !/* 141 */
        SpaceGroup(28, "C2v^4 ", "P -2a -2a       ", "P m 2 a                        ", &
     "P m 2 a            ", "Pm2a      ", "a-cb ", PRIMITIVE, 7), & !/* 142 */
        SpaceGroup(29, "C2v^5 ", "P 2c -2ac       ", "P c a 2_1                      ", &
     "P c a 2_1          ", "Pca2_1    ", "     ", PRIMITIVE, 7), & !/* 143 */
        SpaceGroup(29, "C2v^5 ", "P 2c -2b        ", "P b c 2_1                      ", &
     "P b c 2_1          ", "Pbc2_1    ", "ba-c ", PRIMITIVE, 7), & !/* 144 */
        SpaceGroup(29, "C2v^5 ", "P -2b 2a        ", "P 2_1 a b                      ", &
     "P 2_1 a b          ", "P2_1ab    ", "cab  ", PRIMITIVE, 7), & !/* 145 */
        SpaceGroup(29, "C2v^5 ", "P -2ac 2a       ", "P 2_1 c a                      ", &
     "P 2_1 c a          ", "P2_1ca    ", "-cba ", PRIMITIVE, 7), & !/* 146 */
        SpaceGroup(29, "C2v^5 ", "P -2bc -2c      ", "P c 2_1 b                      ", &
     "P c 2_1 b          ", "Pc2_1b    ", "bca  ", PRIMITIVE, 7), & !/* 147 */
        SpaceGroup(29, "C2v^5 ", "P -2a -2ab      ", "P b 2_1 a                      ", &
     "P b 2_1 a          ", "Pb2_1a    ", "a-cb ", PRIMITIVE, 7), & !/* 148 */
        SpaceGroup(30, "C2v^6 ", "P 2 -2bc        ", "P n c 2                        ", &
     "P n c 2            ", "Pnc2      ", "     ", PRIMITIVE, 7), & !/* 149 */
        SpaceGroup(30, "C2v^6 ", "P 2 -2ac        ", "P c n 2                        ", &
     "P c n 2            ", "Pcn2      ", "ba-c ", PRIMITIVE, 7), & !/* 150 */
        SpaceGroup(30, "C2v^6 ", "P -2ac 2        ", "P 2 n a                        ", &
     "P 2 n a            ", "P2na      ", "cab  ", PRIMITIVE, 7), & !/* 151 */
        SpaceGroup(30, "C2v^6 ", "P -2ab 2        ", "P 2 a n                        ", &
     "P 2 a n            ", "P2an      ", "-cba ", PRIMITIVE, 7), & !/* 152 */
        SpaceGroup(30, "C2v^6 ", "P -2ab -2ab     ", "P b 2 n                        ", &
     "P b 2 n            ", "Pb2n      ", "bca  ", PRIMITIVE, 7), & !/* 153 */
        SpaceGroup(30, "C2v^6 ", "P -2bc -2bc     ", "P n 2 b                        ", &
     "P n 2 b            ", "Pn2b      ", "a-cb ", PRIMITIVE, 7), & !/* 154 */
        SpaceGroup(31, "C2v^7 ", "P 2ac -2        ", "P m n 2_1                      ", &
     "P m n 2_1          ", "Pmn2_1    ", "     ", PRIMITIVE, 7), & !/* 155 */
        SpaceGroup(31, "C2v^7 ", "P 2bc -2bc      ", "P n m 2_1                      ", &
     "P n m 2_1          ", "Pnm2_1    ", "ba-c ", PRIMITIVE, 7), & !/* 156 */
        SpaceGroup(31, "C2v^7 ", "P -2ab 2ab      ", "P 2_1 m n                      ", &
     "P 2_1 m n          ", "P2_1mn    ", "cab  ", PRIMITIVE, 7), & !/* 157 */
        SpaceGroup(31, "C2v^7 ", "P -2 2ac        ", "P 2_1 n m                      ", &
     "P 2_1 n m          ", "P2_1nm    ", "-cba ", PRIMITIVE, 7), & !/* 158 */
        SpaceGroup(31, "C2v^7 ", "P -2 -2bc       ", "P n 2_1 m                      ", &
     "P n 2_1 m          ", "Pn2_1m    ", "bca  ", PRIMITIVE, 7), & !/* 159 */
        SpaceGroup(31, "C2v^7 ", "P -2ab -2       ", "P m 2_1 n                      ", &
     "P m 2_1 n          ", "Pm2_1n    ", "a-cb ", PRIMITIVE, 7), & !/* 160 */
        SpaceGroup(32, "C2v^8 ", "P 2 -2ab        ", "P b a 2                        ", &
     "P b a 2            ", "Pba2      ", "     ", PRIMITIVE, 7), & !/* 161 */
        SpaceGroup(32, "C2v^8 ", "P -2bc 2        ", "P 2 c b                        ", &
     "P 2 c b            ", "P2cb      ", "cab  ", PRIMITIVE, 7), & !/* 162 */
        SpaceGroup(32, "C2v^8 ", "P -2ac -2ac     ", "P c 2 a                        ", &
     "P c 2 a            ", "Pc2a      ", "bca  ", PRIMITIVE, 7), & !/* 163 */
        SpaceGroup(33, "C2v^9 ", "P 2c -2n        ", "P n a 2_1                      ", &
     "P n a 2_1          ", "Pna2_1    ", "     ", PRIMITIVE, 7), & !/* 164 */
        SpaceGroup(33, "C2v^9 ", "P 2c -2ab       ", "P b n 2_1                      ", &
     "P b n 2_1          ", "Pbn2_1    ", "ba-c ", PRIMITIVE, 7), & !/* 165 */
        SpaceGroup(33, "C2v^9 ", "P -2bc 2a       ", "P 2_1 n b                      ", &
     "P 2_1 n b          ", "P2_1nb    ", "cab  ", PRIMITIVE, 7), & !/* 166 */
        SpaceGroup(33, "C2v^9 ", "P -2n 2a        ", "P 2_1 c n                      ", &
     "P 2_1 c n          ", "P2_1cn    ", "-cba ", PRIMITIVE, 7), & !/* 167 */
        SpaceGroup(33, "C2v^9 ", "P -2n -2ac      ", "P c 2_1 n                      ", &
     "P c 2_1 n          ", "Pc2_1n    ", "bca  ", PRIMITIVE, 7), & !/* 168 */
        SpaceGroup(33, "C2v^9 ", "P -2ac -2n      ", "P n 2_1 a                      ", &
     "P n 2_1 a          ", "Pn2_1a    ", "a-cb ", PRIMITIVE, 7), & !/* 169 */
        SpaceGroup(34, "C2v^10", "P 2 -2n         ", "P n n 2                        ", &
     "P n n 2            ", "Pnn2      ", "     ", PRIMITIVE, 7), & !/* 170 */
        SpaceGroup(34, "C2v^10", "P -2n 2         ", "P 2 n n                        ", &
     "P 2 n n            ", "P2nn      ", "cab  ", PRIMITIVE, 7), & !/* 171 */
        SpaceGroup(34, "C2v^10", "P -2n -2n       ", "P n 2 n                        ", &
     "P n 2 n            ", "Pn2n      ", "bca  ", PRIMITIVE, 7), & !/* 172 */
        SpaceGroup(35, "C2v^11", "C 2 -2          ", "C m m 2                        ", &
     "C m m 2            ", "Cmm2      ", "     ", C_FACE, 7), & !/* 173 */
        SpaceGroup(35, "C2v^11", "A -2 2          ", "A 2 m m                        ", &
     "A 2 m m            ", "A2mm      ", "cab  ", A_FACE, 7), & !/* 174 */
        SpaceGroup(35, "C2v^11", "B -2 -2         ", "B m 2 m                        ", &
     "B m 2 m            ", "Bm2m      ", "bca  ", PRIMITIVE, 7), & !/* 175 */
        SpaceGroup(36, "C2v^12", "C 2c -2         ", "C m c 2_1                      ", &
     "C m c 2_1          ", "Cmc2_1    ", "     ", C_FACE, 7), & !/* 176 */
        SpaceGroup(36, "C2v^12", "C 2c -2c        ", "C c m 2_1                      ", &
     "C c m 2_1          ", "Ccm2_1    ", "ba-c ", C_FACE, 7), & !/* 177 */
        SpaceGroup(36, "C2v^12", "A -2a 2a        ", "A 2_1 m a                      ", &
     "A 2_1 m a          ", "A2_1ma    ", "cab  ", A_FACE, 7), & !/* 178 */
        SpaceGroup(36, "C2v^12", "A -2 2a         ", "A 2_1 a m                      ", &
     "A 2_1 a m          ", "A2_1am    ", "-cba ", A_FACE, 7), & !/* 179 */
        SpaceGroup(36, "C2v^12", "B -2 -2b        ", "B b 2_1 m                      ", &
     "B b 2_1 m          ", "Bb2_1m    ", "bca  ", PRIMITIVE, 7), & !/* 180 */
        SpaceGroup(36, "C2v^12", "B -2b -2        ", "B m 2_1 b                      ", &
     "B m 2_1 b          ", "Bm2_1b    ", "a-cb ", PRIMITIVE, 7), & !/* 181 */
        SpaceGroup(37, "C2v^13", "C 2 -2c         ", "C c c 2                        ", &
     "C c c 2            ", "Ccc2      ", "     ", C_FACE, 7), & !/* 182 */
        SpaceGroup(37, "C2v^13", "A -2a 2         ", "A 2 a a                        ", &
     "A 2 a a            ", "A2aa      ", "cab  ", A_FACE, 7), & !/* 183 */
        SpaceGroup(37, "C2v^13", "B -2b -2b       ", "B b 2 b                        ", &
     "B b 2 b            ", "Bb2b      ", "bca  ", PRIMITIVE, 7), & !/* 184 */
        SpaceGroup(38, "C2v^14", "A 2 -2          ", "A m m 2                        ", &
     "A m m 2            ", "Amm2      ", "     ", A_FACE, 7), & !/* 185 */
        SpaceGroup(38, "C2v^14", "B 2 -2          ", "B m m 2                        ", &
     "B m m 2            ", "Bmm2      ", "ba-c ", PRIMITIVE, 7), & !/* 186 */
        SpaceGroup(38, "C2v^14", "B -2 2          ", "B 2 m m                        ", &
     "B 2 m m            ", "B2mm      ", "cab  ", PRIMITIVE, 7), & !/* 187 */
        SpaceGroup(38, "C2v^14", "C -2 2          ", "C 2 m m                        ", &
     "C 2 m m            ", "C2mm      ", "-cba ", C_FACE, 7), & !/* 188 */
        SpaceGroup(38, "C2v^14", "C -2 -2         ", "C m 2 m                        ", &
     "C m 2 m            ", "Cm2m      ", "bca  ", C_FACE, 7), & !/* 189 */
        SpaceGroup(38, "C2v^14", "A -2 -2         ", "A m 2 m                        ", &
     "A m 2 m            ", "Am2m      ", "a-cb ", A_FACE, 7), & !/* 190 */
        SpaceGroup(39, "C2v^15", "A 2 -2b         ", "A e m 2                        ", &
     "A e m 2            ", "Aem2      ", "     ", A_FACE, 7), & !/* 191 */
        SpaceGroup(39, "C2v^15", "B 2 -2a         ", "B m e 2                        ", &
     "B m e 2            ", "Bme2      ", "ba-c ", PRIMITIVE, 7), & !/* 192 */
        SpaceGroup(39, "C2v^15", "B -2a 2         ", "B 2 e m                        ", &
     "B 2 e m            ", "B2em      ", "cab  ", PRIMITIVE, 7), & !/* 193 */
        SpaceGroup(39, "C2v^15", "C -2a 2         ", "C 2 m e                        ", &
     "C 2 m e            ", "C2me      ", "-cba ", C_FACE, 7), & !/* 194 */
        SpaceGroup(39, "C2v^15", "C -2a -2a       ", "C m 2 e                        ", &
     "C m 2 e            ", "Cm2e      ", "bca  ", C_FACE, 7), & !/* 195 */
        SpaceGroup(39, "C2v^15", "A -2b -2b       ", "A e 2 m                        ", &
     "A e 2 m            ", "Ae2m      ", "a-cb ", A_FACE, 7), & !/* 196 */
        SpaceGroup(40, "C2v^16", "A 2 -2a         ", "A m a 2                        ", &
     "A m a 2            ", "Ama2      ", "     ", A_FACE, 7), & !/* 197 */
        SpaceGroup(40, "C2v^16", "B 2 -2b         ", "B b m 2                        ", &
     "B b m 2            ", "Bbm2      ", "ba-c ", PRIMITIVE, 7), & !/* 198 */
        SpaceGroup(40, "C2v^16", "B -2b 2         ", "B 2 m b                        ", &
     "B 2 m b            ", "B2mb      ", "cab  ", PRIMITIVE, 7), & !/* 199 */
        SpaceGroup(40, "C2v^16", "C -2c 2         ", "C 2 c m                        ", &
     "C 2 c m            ", "C2cm      ", "-cba ", C_FACE, 7), & !/* 200 */
        SpaceGroup(40, "C2v^16", "C -2c -2c       ", "C c 2 m                        ", &
     "C c 2 m            ", "Cc2m      ", "bca  ", C_FACE, 7), & !/* 201 */
        SpaceGroup(40, "C2v^16", "A -2a -2a       ", "A m 2 a                        ", &
     "A m 2 a            ", "Am2a      ", "a-cb ", A_FACE, 7), & !/* 202 */
        SpaceGroup(41, "C2v^17", "A 2 -2ab        ", "A e a 2                        ", &
     "A e a 2            ", "Aea2      ", "     ", A_FACE, 7), & !/* 203 */
        SpaceGroup(41, "C2v^17", "B 2 -2ab        ", "B b e 2                        ", &
     "B b e 2            ", "Bbe2      ", "ba-c ", PRIMITIVE, 7), & !/* 204 */
        SpaceGroup(41, "C2v^17", "B -2ab 2        ", "B 2 e b                        ", &
     "B 2 e b            ", "B2eb      ", "cab  ", PRIMITIVE, 7), & !/* 205 */
        SpaceGroup(41, "C2v^17", "C -2ac 2        ", "C 2 c e                        ", &
     "C 2 c e            ", "C2ce      ", "-cba ", C_FACE, 7), & !/* 206 */
        SpaceGroup(41, "C2v^17", "C -2ac -2ac     ", "C c 2 e                        ", &
     "C c 2 e            ", "Cc2e      ", "bca  ", C_FACE, 7), & !/* 207 */
        SpaceGroup(41, "C2v^17", "A -2ab -2ab     ", "A e 2 a                        ", &
     "A e 2 a            ", "Ae2a      ", "a-cb ", A_FACE, 7), & !/* 208 */
        SpaceGroup(42, "C2v^18", "F 2 -2          ", "F m m 2                        ", &
     "F m m 2            ", "Fmm2      ", "     ", FACE, 7), & !/* 209 */
        SpaceGroup(42, "C2v^18", "F -2 2          ", "F 2 m m                        ", &
     "F 2 m m            ", "F2mm      ", "cab  ", FACE, 7), & !/* 210 */
        SpaceGroup(42, "C2v^18", "F -2 -2         ", "F m 2 m                        ", &
     "F m 2 m            ", "Fm2m      ", "bca  ", FACE, 7), & !/* 211 */
        SpaceGroup(43, "C2v^19", "F 2 -2d         ", "F d d 2                        ", &
     "F d d 2            ", "Fdd2      ", "     ", FACE, 7), & !/* 212 */
        SpaceGroup(43, "C2v^19", "F -2d 2         ", "F 2 d d                        ", &
     "F 2 d d            ", "F2dd      ", "cab  ", FACE, 7), & !/* 213 */
        SpaceGroup(43, "C2v^19", "F -2d -2d       ", "F d 2 d                        ", &
     "F d 2 d            ", "Fd2d      ", "bca  ", FACE, 7), & !/* 214 */
        SpaceGroup(44, "C2v^20", "I 2 -2          ", "I m m 2                        ", &
     "I m m 2            ", "Imm2      ", "     ", BODY, 7), & !/* 215 */
        SpaceGroup(44, "C2v^20", "I -2 2          ", "I 2 m m                        ", &
     "I 2 m m            ", "I2mm      ", "cab  ", BODY, 7), & !/* 216 */
        SpaceGroup(44, "C2v^20", "I -2 -2         ", "I m 2 m                        ", &
     "I m 2 m            ", "Im2m      ", "bca  ", BODY, 7), & !/* 217 */
        SpaceGroup(45, "C2v^21", "I 2 -2c         ", "I b a 2                        ", &
     "I b a 2            ", "Iba2      ", "     ", BODY, 7), & !/* 218 */
        SpaceGroup(45, "C2v^21", "I -2a 2         ", "I 2 c b                        ", &
     "I 2 c b            ", "I2cb      ", "cab  ", BODY, 7), & !/* 219 */
        SpaceGroup(45, "C2v^21", "I -2b -2b       ", "I c 2 a                        ", &
     "I c 2 a            ", "Ic2a      ", "bca  ", BODY, 7), & !/* 220 */
        SpaceGroup(46, "C2v^22", "I 2 -2a         ", "I m a 2                        ", &
     "I m a 2            ", "Ima2      ", "     ", BODY, 7), & !/* 221 */
        SpaceGroup(46, "C2v^22", "I 2 -2b         ", "I b m 2                        ", &
     "I b m 2            ", "Ibm2      ", "ba-c ", BODY, 7), & !/* 222 */
        SpaceGroup(46, "C2v^22", "I -2b 2         ", "I 2 m b                        ", &
     "I 2 m b            ", "I2mb      ", "cab  ", BODY, 7), & !/* 223 */
        SpaceGroup(46, "C2v^22", "I -2c 2         ", "I 2 c m                        ", &
     "I 2 c m            ", "I2cm      ", "-cba ", BODY, 7), & !/* 224 */
        SpaceGroup(46, "C2v^22", "I -2c -2c       ", "I c 2 m                        ", &
     "I c 2 m            ", "Ic2m      ", "bca  ", BODY, 7), & !/* 225 */
        SpaceGroup(46, "C2v^22", "I -2a -2a       ", "I m 2 a                        ", &
     "I m 2 a            ", "Im2a      ", "a-cb ", BODY, 7), & !/* 226 */
        SpaceGroup(47, "D2h^1 ", "-P 2 2          ", "P m m m                        ", &
     "P 2/m 2/m 2/m      ", "Pmmm      ", "     ", PRIMITIVE, 8), & !/* 227 */
        SpaceGroup(48, "D2h^2 ", "P 2 2 -1n       ", "P n n n                        ", &
     "P 2/n 2/n 2/n      ", "Pnnn      ", "1    ", PRIMITIVE, 8), & !/* 228 */
        SpaceGroup(48, "D2h^2 ", "-P 2ab 2bc      ", "P n n n                        ", &
     "P 2/n 2/n 2/n      ", "Pnnn      ", "2    ", PRIMITIVE, 8), & !/* 229 */
        SpaceGroup(49, "D2h^3 ", "-P 2 2c         ", "P c c m                        ", &
     "P 2/c 2/c 2/m      ", "Pccm      ", "     ", PRIMITIVE, 8), & !/* 230 */
        SpaceGroup(49, "D2h^3 ", "-P 2a 2         ", "P m a a                        ", &
     "P 2/m 2/a 2/a      ", "Pmaa      ", "cab  ", PRIMITIVE, 8), & !/* 231 */
        SpaceGroup(49, "D2h^3 ", "-P 2b 2b        ", "P b m b                        ", &
     "P 2/b 2/m 2/b      ", "Pbmb      ", "bca  ", PRIMITIVE, 8), & !/* 232 */
        SpaceGroup(50, "D2h^4 ", "P 2 2 -1ab      ", "P b a n                        ", &
     "P 2/b 2/a 2/n      ", "Pban      ", "1    ", PRIMITIVE, 8), & !/* 233 */
        SpaceGroup(50, "D2h^4 ", "-P 2ab 2b       ", "P b a n                        ", &
     "P 2/b 2/a 2/n      ", "Pban      ", "2    ", PRIMITIVE, 8), & !/* 234 */
        SpaceGroup(50, "D2h^4 ", "P 2 2 -1bc      ", "P n c b                        ", &
     "P 2/n 2/c 2/b      ", "Pncb      ", "1cab ", PRIMITIVE, 8), & !/* 235 */
        SpaceGroup(50, "D2h^4 ", "-P 2b 2bc       ", "P n c b                        ", &
     "P 2/n 2/c 2/b      ", "Pncb      ", "2cab ", PRIMITIVE, 8), & !/* 236 */
        SpaceGroup(50, "D2h^4 ", "P 2 2 -1ac      ", "P c n a                        ", &
     "P 2/c 2/n 2/a      ", "Pcna      ", "1bca ", PRIMITIVE, 8), & !/* 237 */
        SpaceGroup(50, "D2h^4 ", "-P 2a 2c        ", "P c n a                        ", &
     "P 2/c 2/n 2/a      ", "Pcna      ", "2bca ", PRIMITIVE, 8), & !/* 238 */
        SpaceGroup(51, "D2h^5 ", "-P 2a 2a        ", "P m m a                        ", &
     "P 2_1/m 2/m 2/a    ", "Pmma      ", "     ", PRIMITIVE, 8), & !/* 239 */
        SpaceGroup(51, "D2h^5 ", "-P 2b 2         ", "P m m b                        ", &
     "P 2/m 2_1/m 2/b    ", "Pmmb      ", "ba-c ", PRIMITIVE, 8), & !/* 240 */
        SpaceGroup(51, "D2h^5 ", "-P 2 2b         ", "P b m m                        ", &
     "P 2/b 2_1/m 2/m    ", "Pbmm      ", "cab  ", PRIMITIVE, 8), & !/* 241 */
        SpaceGroup(51, "D2h^5 ", "-P 2c 2c        ", "P c m m                        ", &
     "P 2/c 2/m 2_1/m    ", "Pcmm      ", "-cba ", PRIMITIVE, 8), & !/* 242 */
        SpaceGroup(51, "D2h^5 ", "-P 2c 2         ", "P m c m                        ", &
     "P 2/m 2/c 2_1/m    ", "Pmcm      ", "bca  ", PRIMITIVE, 8), & !/* 243 */
        SpaceGroup(51, "D2h^5 ", "-P 2 2a         ", "P m a m                        ", &
     "P 2_1/m 2/a 2/m    ", "Pmam      ", "a-cb ", PRIMITIVE, 8), & !/* 244 */
        SpaceGroup(52, "D2h^6 ", "-P 2a 2bc       ", "P n n a                        ", &
     "P 2/n 2_1/n 2/a    ", "Pnna      ", "     ", PRIMITIVE, 8), & !/* 245 */
        SpaceGroup(52, "D2h^6 ", "-P 2b 2n        ", "P n n b                        ", &
     "P 2_1/n 2/n 2/b    ", "Pnnb      ", "ba-c ", PRIMITIVE, 8), & !/* 246 */
        SpaceGroup(52, "D2h^6 ", "-P 2n 2b        ", "P b n n                        ", &
     "P 2/b 2/n 2_1/n    ", "Pbnn      ", "cab  ", PRIMITIVE, 8), & !/* 247 */
        SpaceGroup(52, "D2h^6 ", "-P 2ab 2c       ", "P c n n                        ", &
     "P 2/c 2_1/n 2/n    ", "Pcnn      ", "-cba ", PRIMITIVE, 8), & !/* 248 */
        SpaceGroup(52, "D2h^6 ", "-P 2ab 2n       ", "P n c n                        ", &
     "P 2_1/n 2/c 2/n    ", "Pncn      ", "bca  ", PRIMITIVE, 8), & !/* 249 */
        SpaceGroup(52, "D2h^6 ", "-P 2n 2bc       ", "P n a n                        ", &
     "P 2/n 2/a 2_1/n    ", "Pnan      ", "a-cb ", PRIMITIVE, 8), & !/* 250 */
        SpaceGroup(53, "D2h^7 ", "-P 2ac 2        ", "P m n a                        ", &
     "P 2/m 2/n 2_1/a    ", "Pmna      ", "     ", PRIMITIVE, 8), & !/* 251 */
        SpaceGroup(53, "D2h^7 ", "-P 2bc 2bc      ", "P n m b                        ", &
     "P 2/n 2/m 2_1/b    ", "Pnmb      ", "ba-c ", PRIMITIVE, 8), & !/* 252 */
        SpaceGroup(53, "D2h^7 ", "-P 2ab 2ab      ", "P b m n                        ", &
     "P 2_1/b 2/m 2/n    ", "Pbmn      ", "cab  ", PRIMITIVE, 8), & !/* 253 */
        SpaceGroup(53, "D2h^7 ", "-P 2 2ac        ", "P c n m                        ", &
     "P 2_1/c 2/n 2/m    ", "Pcnm      ", "-cba ", PRIMITIVE, 8), & !/* 254 */
        SpaceGroup(53, "D2h^7 ", "-P 2 2bc        ", "P n c m                        ", &
     "P 2/n 2_1/c 2/m    ", "Pncm      ", "bca  ", PRIMITIVE, 8), & !/* 255 */
        SpaceGroup(53, "D2h^7 ", "-P 2ab 2        ", "P m a n                        ", &
     "P 2/m 2_1/a 2/n    ", "Pman      ", "a-cb ", PRIMITIVE, 8), & !/* 256 */
        SpaceGroup(54, "D2h^8 ", "-P 2a 2ac       ", "P c c a                        ", &
     "P 2_1/c 2/c 2/a    ", "Pcca      ", "     ", PRIMITIVE, 8), & !/* 257 */
        SpaceGroup(54, "D2h^8 ", "-P 2b 2c        ", "P c c b                        ", &
     "P 2/c 2_1/c 2/b    ", "Pccb      ", "ba-c ", PRIMITIVE, 8), & !/* 258 */
        SpaceGroup(54, "D2h^8 ", "-P 2a 2b        ", "P b a a                        ", &
     "P 2/b 2_1/a 2/a    ", "Pbaa      ", "cab  ", PRIMITIVE, 8), & !/* 259 */
        SpaceGroup(54, "D2h^8 ", "-P 2ac 2c       ", "P c a a                        ", &
     "P 2/c 2/a 2_1/a    ", "Pcaa      ", "-cba ", PRIMITIVE, 8), & !/* 260 */
        SpaceGroup(54, "D2h^8 ", "-P 2bc 2b       ", "P b c b                        ", &
     "P 2/b 2/c 2_1/b    ", "Pbcb      ", "bca  ", PRIMITIVE, 8), & !/* 261 */
        SpaceGroup(54, "D2h^8 ", "-P 2b 2ab       ", "P b a b                        ", &
     "P 2_1/b 2/a 2/b    ", "Pbab      ", "a-cb ", PRIMITIVE, 8), & !/* 262 */
        SpaceGroup(55, "D2h^9 ", "-P 2 2ab        ", "P b a m                        ", &
     "P 2_1/b 2_1/a 2/m  ", "Pbam      ", "     ", PRIMITIVE, 8), & !/* 263 */
        SpaceGroup(55, "D2h^9 ", "-P 2bc 2        ", "P m c b                        ", &
     "P 2/m 2_1/c 2_1/b  ", "Pmcb      ", "cab  ", PRIMITIVE, 8), & !/* 264 */
        SpaceGroup(55, "D2h^9 ", "-P 2ac 2ac      ", "P c m a                        ", &
     "P 2_1/c 2/m 2_1/a  ", "Pcma      ", "bca  ", PRIMITIVE, 8), & !/* 265 */
        SpaceGroup(56, "D2h^10", "-P 2ab 2ac      ", "P c c n                        ", &
     "P 2_1/c 2_1/c 2/n  ", "Pccn      ", "     ", PRIMITIVE, 8), & !/* 266 */
        SpaceGroup(56, "D2h^10", "-P 2ac 2bc      ", "P n a a                        ", &
     "P 2/n 2_1/a 2_1/a  ", "Pnaa      ", "cab  ", PRIMITIVE, 8), & !/* 267 */
        SpaceGroup(56, "D2h^10", "-P 2bc 2ab      ", "P b n b                        ", &
     "P 2_1/b 2/n 2_1/b  ", "Pbnb      ", "bca  ", PRIMITIVE, 8), & !/* 268 */
        SpaceGroup(57, "D2h^11", "-P 2c 2b        ", "P b c m                        ", &
     "P 2/b 2_1/c 2_1/m  ", "Pbcm      ", "     ", PRIMITIVE, 8), & !/* 269 */
        SpaceGroup(57, "D2h^11", "-P 2c 2ac       ", "P c a m                        ", &
     "P 2_1/c 2/a 2_1/m  ", "Pcam      ", "ba-c ", PRIMITIVE, 8), & !/* 270 */
        SpaceGroup(57, "D2h^11", "-P 2ac 2a       ", "P m c a                        ", &
     "P 2_1/m 2/c 2_1/a  ", "Pmca      ", "cab  ", PRIMITIVE, 8), & !/* 271 */
        SpaceGroup(57, "D2h^11", "-P 2b 2a        ", "P m a b                        ", &
     "P 2_1/m 2_1/a 2/b  ", "Pmab      ", "-cba ", PRIMITIVE, 8), & !/* 272 */
        SpaceGroup(57, "D2h^11", "-P 2a 2ab       ", "P b m a                        ", &
     "P 2_1/b 2_1/m 2/a  ", "Pbma      ", "bca  ", PRIMITIVE, 8), & !/* 273 */
        SpaceGroup(57, "D2h^11", "-P 2bc 2c       ", "P c m b                        ", &
     "P 2/c 2_1/m 2_1/b  ", "Pcmb      ", "a-cb ", PRIMITIVE, 8), & !/* 274 */
        SpaceGroup(58, "D2h^12", "-P 2 2n         ", "P n n m                        ", &
     "P 2_1/n 2_1/n 2/m  ", "Pnnm      ", "     ", PRIMITIVE, 8), & !/* 275 */
        SpaceGroup(58, "D2h^12", "-P 2n 2         ", "P m n n                        ", &
     "P 2/m 2_1/n 2_1/n  ", "Pmnn      ", "cab  ", PRIMITIVE, 8), & !/* 276 */
        SpaceGroup(58, "D2h^12", "-P 2n 2n        ", "P n m n                        ", &
     "P 2_1/n 2/m 2_1/n  ", "Pnmn      ", "bca  ", PRIMITIVE, 8), & !/* 277 */
        SpaceGroup(59, "D2h^13", "P 2 2ab -1ab    ", "P m m n                        ", &
     "P 2_1/m 2_1/m 2/n  ", "Pmmn      ", "1    ", PRIMITIVE, 8), & !/* 278 */
        SpaceGroup(59, "D2h^13", "-P 2ab 2a       ", "P m m n                        ", &
     "P 2_1/m 2_1/m 2/n  ", "Pmmn      ", "2    ", PRIMITIVE, 8), & !/* 279 */
        SpaceGroup(59, "D2h^13", "P 2bc 2 -1bc    ", "P n m m                        ", &
     "P 2/n 2_1/m 2_1/m  ", "Pnmm      ", "1cab ", PRIMITIVE, 8), & !/* 280 */
        SpaceGroup(59, "D2h^13", "-P 2c 2bc       ", "P n m m                        ", &
     "P 2/n 2_1/m 2_1/m  ", "Pnmm      ", "2cab ", PRIMITIVE, 8), & !/* 281 */
        SpaceGroup(59, "D2h^13", "P 2ac 2ac -1ac  ", "P m n m                        ", &
     "P 2_1/m 2/n 2_1/m  ", "Pmnm      ", "1bca ", PRIMITIVE, 8), & !/* 282 */
        SpaceGroup(59, "D2h^13", "-P 2c 2a        ", "P m n m                        ", &
     "P 2_1/m 2/n 2_1/m  ", "Pmnm      ", "2bca ", PRIMITIVE, 8), & !/* 283 */
        SpaceGroup(60, "D2h^14", "-P 2n 2ab       ", "P b c n                        ", &
     "P 2_1/b 2/c 2_1/n  ", "Pbcn      ", "     ", PRIMITIVE, 8), & !/* 284 */
        SpaceGroup(60, "D2h^14", "-P 2n 2c        ", "P c a n                        ", &
     "P 2/c 2_1/a 2_1/n  ", "Pcan      ", "ba-c ", PRIMITIVE, 8), & !/* 285 */
        SpaceGroup(60, "D2h^14", "-P 2a 2n        ", "P n c a                        ", &
     "P 2_1/n 2_1/c 2/a  ", "Pnca      ", "cab  ", PRIMITIVE, 8), & !/* 286 */
        SpaceGroup(60, "D2h^14", "-P 2bc 2n       ", "P n a b                        ", &
     "P 2_1/n 2/a 2_1/b  ", "Pnab      ", "-cba ", PRIMITIVE, 8), & !/* 287 */
        SpaceGroup(60, "D2h^14", "-P 2ac 2b       ", "P b n a                        ", &
     "P 2/b 2_1/n 2_1/a  ", "Pbna      ", "bca  ", PRIMITIVE, 8), & !/* 288 */
        SpaceGroup(60, "D2h^14", "-P 2b 2ac       ", "P c n b                        ", &
     "P 2_1/c 2_1/n 2/b  ", "Pcnb      ", "a-cb ", PRIMITIVE, 8), & !/* 289 */
        SpaceGroup(61, "D2h^15", "-P 2ac 2ab      ", "P b c a                        ", &
     "P 2_1/b 2_1/c 2_1/a", "Pbca      ", "     ", PRIMITIVE, 8), & !/* 290 */
        SpaceGroup(61, "D2h^15", "-P 2bc 2ac      ", "P c a b                        ", &
     "P 2_1/c 2_1/a 2_1/b", "Pcab      ", "ba-c ", PRIMITIVE, 8), & !/* 291 */
        SpaceGroup(62, "D2h^16", "-P 2ac 2n       ", "P n m a                        ", &
     "P 2_1/n 2_1/m 2_1/a", "Pnma      ", "     ", PRIMITIVE, 8), & !/* 292 */
        SpaceGroup(62, "D2h^16", "-P 2bc 2a       ", "P m n b                        ", &
     "P 2_1/m 2_1/n 2_1/b", "Pmnb      ", "ba-c ", PRIMITIVE, 8), & !/* 293 */
        SpaceGroup(62, "D2h^16", "-P 2c 2ab       ", "P b n m                        ", &
     "P 2_1/b 2_1/n 2_1/m", "Pbnm      ", "cab  ", PRIMITIVE, 8), & !/* 294 */
        SpaceGroup(62, "D2h^16", "-P 2n 2ac       ", "P c m n                        ", &
     "P 2_1/c 2_1/m 2_1/n", "Pcmn      ", "-cba ", PRIMITIVE, 8), & !/* 295 */
        SpaceGroup(62, "D2h^16", "-P 2n 2a        ", "P m c n                        ", &
     "P 2_1/m 2_1/c 2_1/n", "Pmcn      ", "bca  ", PRIMITIVE, 8), & !/* 296 */
        SpaceGroup(62, "D2h^16", "-P 2c 2n        ", "P n a m                        ", &
     "P 2_1/n 2_1/a 2_1/m", "Pnam      ", "a-cb ", PRIMITIVE, 8), & !/* 297 */
        SpaceGroup(63, "D2h^17", "-C 2c 2         ", "C m c m                        ", &
     "C 2/m 2/c 2_1/m    ", "Cmcm      ", "     ", C_FACE, 8), & !/* 298 */
        SpaceGroup(63, "D2h^17", "-C 2c 2c        ", "C c m m                        ", &
     "C 2/c 2/m 2_1/m    ", "Ccmm      ", "ba-c ", C_FACE, 8), & !/* 299 */
        SpaceGroup(63, "D2h^17", "-A 2a 2a        ", "A m m a                        ", &
     "A 2_1/m 2/m 2/a    ", "Amma      ", "cab  ", A_FACE, 8), & !/* 300 */
        SpaceGroup(63, "D2h^17", "-A 2 2a         ", "A m a m                        ", &
     "A 2_1/m 2/a 2/m    ", "Amam      ", "-cba ", A_FACE, 8), & !/* 301 */
        SpaceGroup(63, "D2h^17", "-B 2 2b         ", "B b m m                        ", &
     "B 2/b 2_1/m 2/m    ", "Bbmm      ", "bca  ", PRIMITIVE, 8), & !/* 302 */
        SpaceGroup(63, "D2h^17", "-B 2b 2         ", "B m m b                        ", &
     "B 2/m 2_1/m 2/b    ", "Bmmb      ", "a-cb ", PRIMITIVE, 8), & !/* 303 */
        SpaceGroup(64, "D2h^18", "-C 2ac 2        ", "C m c e                        ", &
     "C 2/m 2/c 2_1/e    ", "Cmce      ", "     ", C_FACE, 8), & !/* 304 */
        SpaceGroup(64, "D2h^18", "-C 2ac 2ac      ", "C c m e                        ", &
     "C 2/c 2/m 2_1/e    ", "Ccme      ", "ba-c ", C_FACE, 8), & !/* 305 */
        SpaceGroup(64, "D2h^18", "-A 2ab 2ab      ", "A e m a                        ", &
     "A 2_1/e 2/m 2/a    ", "Aema      ", "cab  ", A_FACE, 8), & !/* 306 */
        SpaceGroup(64, "D2h^18", "-A 2 2ab        ", "A e a m                        ", &
     "A 2_1/e 2/a 2/m    ", "Aeam      ", "-cba ", A_FACE, 8), & !/* 307 */
        SpaceGroup(64, "D2h^18", "-B 2 2ab        ", "B b e m                        ", &
     "B 2/b 2_1/e 2/m    ", "Bbem      ", "bca  ", PRIMITIVE, 8), & !/* 308 */
        SpaceGroup(64, "D2h^18", "-B 2ab 2        ", "B m e b                        ", &
     "B 2/m 2_1/e 2/b    ", "Bmeb      ", "a-cb ", PRIMITIVE, 8), & !/* 309 */
        SpaceGroup(65, "D2h^19", "-C 2 2          ", "C m m m                        ", &
     "C 2/m 2/m 2/m      ", "Cmmm      ", "     ", C_FACE, 8), & !/* 310 */
        SpaceGroup(65, "D2h^19", "-A 2 2          ", "A m m m                        ", &
     "A 2/m 2/m 2/m      ", "Ammm      ", "cab  ", A_FACE, 8), & !/* 311 */
        SpaceGroup(65, "D2h^19", "-B 2 2          ", "B m m m                        ", &
     "B 2/m 2/m 2/m      ", "Bmmm      ", "bca  ", PRIMITIVE, 8), & !/* 312 */
        SpaceGroup(66, "D2h^20", "-C 2 2c         ", "C c c m                        ", &
     "C 2/c 2/c 2/m      ", "Cccm      ", "     ", C_FACE, 8), & !/* 313 */
        SpaceGroup(66, "D2h^20", "-A 2a 2         ", "A m a a                        ", &
     "A 2/m 2/a 2/a      ", "Amaa      ", "cab  ", A_FACE, 8), & !/* 314 */
        SpaceGroup(66, "D2h^20", "-B 2b 2b        ", "B b m b                        ", &
     "B 2/b 2/m 2/b      ", "Bbmb      ", "bca  ", PRIMITIVE, 8), & !/* 315 */
        SpaceGroup(67, "D2h^21", "-C 2a 2         ", "C m m e                        ", &
     "C 2/m 2/m 2/e      ", "Cmme      ", "     ", C_FACE, 8), & !/* 316 */
        SpaceGroup(67, "D2h^21", "-C 2a 2a        ", "C m m e                        ", &
     "C 2/m 2/m 2/e      ", "Cmme      ", "ba-c ", C_FACE, 8), & !/* 317 */
        SpaceGroup(67, "D2h^21", "-A 2b 2b        ", "A e m m                        ", &
     "A 2/e 2/m 2/m      ", "Aemm      ", "cab  ", A_FACE, 8), & !/* 318 */
        SpaceGroup(67, "D2h^21", "-A 2 2b         ", "A e m m                        ", &
     "A 2/e 2/m 2/m      ", "Aemm      ", "-cba ", A_FACE, 8), & !/* 319 */
        SpaceGroup(67, "D2h^21", "-B 2 2a         ", "B m e m                        ", &
     "B 2/m 2/e 2/m      ", "Bmem      ", "bca  ", PRIMITIVE, 8), & !/* 320 */
        SpaceGroup(67, "D2h^21", "-B 2a 2         ", "B m e m                        ", &
     "B 2/m 2/e 2/m      ", "Bmem      ", "a-cb ", PRIMITIVE, 8), & !/* 321 */
        SpaceGroup(68, "D2h^22", "C 2 2 -1ac      ", "C c c e                        ", &
     "C 2/c 2/c 2/e      ", "Ccce      ", "1    ", C_FACE, 8), & !/* 322 */
        SpaceGroup(68, "D2h^22", "-C 2a 2ac       ", "C c c e                        ", &
     "C 2/c 2/c 2/e      ", "Ccce      ", "2    ", C_FACE, 8), & !/* 323 */
        SpaceGroup(68, "D2h^22", "C 2 2 -1ac      ", "C c c e                        ", &
     "C 2/c 2/c 2/e      ", "Ccce      ", "1ba-c", C_FACE, 8), & !/* 324 */
        SpaceGroup(68, "D2h^22", "-C 2a 2c        ", "C c c e                        ", &
     "C 2/c 2/c 2/e      ", "Ccce      ", "2ba-c", C_FACE, 8), & !/* 325 */
        SpaceGroup(68, "D2h^22", "A 2 2 -1ab      ", "A e a a                        ", &
     "A 2/e 2/a 2/a      ", "Aeaa      ", "1cab ", A_FACE, 8), & !/* 326 */
        SpaceGroup(68, "D2h^22", "-A 2a 2b        ", "A e a a                        ", &
     "A 2/e 2/a 2/a      ", "Aeaa      ", "2cab ", A_FACE, 8), & !/* 327 */
        SpaceGroup(68, "D2h^22", "A 2 2 -1ab      ", "A e a a                        ", &
     "A 2/e 2/a 2/a      ", "Aeaa      ", "1-cba", A_FACE, 8), & !/* 328 */
        SpaceGroup(68, "D2h^22", "-A 2ab 2b       ", "A e a a                        ", &
     "A 2/e 2/a 2/a      ", "Aeaa      ", "2-cba", A_FACE, 8), & !/* 329 */
        SpaceGroup(68, "D2h^22", "B 2 2 -1ab      ", "B b e b                        ", &
     "B 2/b 2/e 2/b      ", "Bbeb      ", "1bca ", PRIMITIVE, 8), & !/* 330 */
        SpaceGroup(68, "D2h^22", "-B 2ab 2b       ", "B b c b                        ", &
     "B 2/b 2/e 2/b      ", "Bbcb      ", "2bca ", PRIMITIVE, 8), & !/* 331 */
        SpaceGroup(68, "D2h^22", "B 2 2 -1ab      ", "B b e b                        ", &
     "B 2/b 2/e 2/b      ", "Bbeb      ", "1a-cb", PRIMITIVE, 8), & !/* 332 */
        SpaceGroup(68, "D2h^22", "-B 2b 2ab       ", "B b e b                        ", &
     "B 2/b 2/e 2/b      ", "Bbeb      ", "2a-cb", PRIMITIVE, 8), & !/* 333 */
        SpaceGroup(69, "D2h^23", "-F 2 2          ", "F m m m                        ", &
     "F 2/m 2/m 2/m      ", "Fmmm      ", "     ", FACE, 8), & !/* 334 */
        SpaceGroup(70, "D2h^24", "F 2 2 -1d       ", "F d d d                        ", &
     "F 2/d 2/d 2/d      ", "Fddd      ", "1    ", FACE, 8), & !/* 335 */
        SpaceGroup(70, "D2h^24", "-F 2uv 2vw      ", "F d d d                        ", &
     "F 2/d 2/d 2/d      ", "Fddd      ", "2    ", FACE, 8), & !/* 336 */
        SpaceGroup(71, "D2h^25", "-I 2 2          ", "I m m m                        ", &
     "I 2/m 2/m 2/m      ", "Immm      ", "     ", BODY, 8), & !/* 337 */
        SpaceGroup(72, "D2h^26", "-I 2 2c         ", "I b a m                        ", &
     "I 2/b 2/a 2/m      ", "Ibam      ", "     ", BODY, 8), & !/* 338 */
        SpaceGroup(72, "D2h^26", "-I 2a 2         ", "I m c b                        ", &
     "I 2/m 2/c 2/b      ", "Imcb      ", "cab  ", BODY, 8), & !/* 339 */
        SpaceGroup(72, "D2h^26", "-I 2b 2b        ", "I c m a                        ", &
     "I 2/c 2/m 2/a      ", "Icma      ", "bca  ", BODY, 8), & !/* 340 */
        SpaceGroup(73, "D2h^27", "-I 2b 2c        ", "I b c a                        ", &
     "I 2/b 2/c 2/a      ", "Ibca      ", "     ", BODY, 8), & !/* 341 */
        SpaceGroup(73, "D2h^27", "-I 2a 2b        ", "I c a b                        ", &
     "I 2/c 2/a 2/b      ", "Icab      ", "ba-c ", BODY, 8), & !/* 342 */
        SpaceGroup(74, "D2h^28", "-I 2b 2         ", "I m m a                        ", &
     "I 2/m 2/m 2/a      ", "Imma      ", "     ", BODY, 8), & !/* 343 */
        SpaceGroup(74, "D2h^28", "-I 2a 2a        ", "I m m b                        ", &
     "I 2/m 2/m 2/b      ", "Immb      ", "ba-c ", BODY, 8), & !/* 344 */
        SpaceGroup(74, "D2h^28", "-I 2c 2c        ", "I b m m                        ", &
     "I 2/b 2/m 2/m      ", "Ibmm      ", "cab  ", BODY, 8), & !/* 345 */
        SpaceGroup(74, "D2h^28", "-I 2 2b         ", "I c m m                        ", &
     "I 2/c 2/m 2/m      ", "Icmm      ", "-cba ", BODY, 8), & !/* 346 */
        SpaceGroup(74, "D2h^28", "-I 2 2a         ", "I m c m                        ", &
     "I 2/m 2/c 2/m      ", "Imcm      ", "bca  ", BODY, 8), & !/* 347 */
        SpaceGroup(74, "D2h^28", "-I 2c 2         ", "I m a m                        ", &
     "I 2/m 2/a 2/m      ", "Imam      ", "a-cb ", BODY, 8), & !/* 348 */
        SpaceGroup(75, "C4^1  ", "P 4             ", "P 4                            ", &
     "P 4                ", "P4        ", "     ", PRIMITIVE, 9), & !/* 349 */
        SpaceGroup(76, "C4^2  ", "P 4w            ", "P 4_1                          ", &
     "P 4_1              ", "P4_1      ", "     ", PRIMITIVE, 9), & !/* 350 */
        SpaceGroup(77, "C4^3  ", "P 4c            ", "P 4_2                          ", &
     "P 4_2              ", "P4_2      ", "     ", PRIMITIVE, 9), & !/* 351 */
        SpaceGroup(78, "C4^4  ", "P 4cw           ", "P 4_3                          ", &
     "P 4_3              ", "P4_3      ", "     ", PRIMITIVE, 9), & !/* 352 */
        SpaceGroup(79, "C4^5  ", "I 4             ", "I 4                            ", &
     "I 4                ", "I4        ", "     ", BODY, 9), & !/* 353 */
        SpaceGroup(80, "C4^6  ", "I 4bw           ", "I 4_1                          ", &
     "I 4_1              ", "I4_1      ", "     ", BODY, 9), & !/* 354 */
        SpaceGroup(81, "S4^1  ", "P -4            ", "P -4                           ", &
     "P -4               ", "P-4       ", "     ", PRIMITIVE, 10), & !/* 355 */
        SpaceGroup(82, "S4^2  ", "I -4            ", "I -4                           ", &
     "I -4               ", "I-4       ", "     ", BODY, 10), & !/* 356 */
        SpaceGroup(83, "C4h^1 ", "-P 4            ", "P 4/m                          ", &
     "P 4/m              ", "P4/m      ", "     ", PRIMITIVE, 11), & !/* 357 */
        SpaceGroup(84, "C4h^2 ", "-P 4c           ", "P 4_2/m                        ", &
     "P 4_2/m            ", "P4_2/m    ", "     ", PRIMITIVE, 11), & !/* 358 */
        SpaceGroup(85, "C4h^3 ", "P 4ab -1ab      ", "P 4/n                          ", &
     "P 4/n              ", "P4/n      ", "1    ", PRIMITIVE, 11), & !/* 359 */
        SpaceGroup(85, "C4h^3 ", "-P 4a           ", "P 4/n                          ", &
     "P 4/n              ", "P4/n      ", "2    ", PRIMITIVE, 11), & !/* 360 */
        SpaceGroup(86, "C4h^4 ", "P 4n -1n        ", "P 4_2/n                        ", &
     "P 4_2/n            ", "P4_2/n    ", "1    ", PRIMITIVE, 11), & !/* 361 */
        SpaceGroup(86, "C4h^4 ", "-P 4bc          ", "P 4_2/n                        ", &
     "P 4_2/n            ", "P4_2/n    ", "2    ", PRIMITIVE, 11), & !/* 362 */
        SpaceGroup(87, "C4h^5 ", "-I 4            ", "I 4/m                          ", &
     "I 4/m              ", "I4/m      ", "     ", BODY, 11), & !/* 363 */
        SpaceGroup(88, "C4h^6 ", "I 4bw -1bw      ", "I 4_1/a                        ", &
     "I 4_1/a            ", "I4_1/a    ", "1    ", BODY, 11), & !/* 364 */
        SpaceGroup(88, "C4h^6 ", "-I 4ad          ", "I 4_1/a                        ", &
     "I 4_1/a            ", "I4_1/a    ", "2    ", BODY, 11), & !/* 365 */
        SpaceGroup(89, "D4^1  ", "P 4 2           ", "P 4 2 2                        ", &
     "P 4 2 2            ", "P422      ", "     ", PRIMITIVE, 12), & !/* 366 */
        SpaceGroup(90, "D4^2  ", "P 4ab 2ab       ", "P 4 2_1 2                      ", &
     "P 4 2_1 2          ", "P42_12    ", "     ", PRIMITIVE, 12), & !/* 367 */
        SpaceGroup(91, "D4^3  ", "P 4w 2c         ", "P 4_1 2 2                      ", &
     "P 4_1 2 2          ", "P4_122    ", "     ", PRIMITIVE, 12), & !/* 368 */
        SpaceGroup(92, "D4^4  ", "P 4abw 2nw      ", "P 4_1 2_1 2                    ", &
     "P 4_1 2_1 2        ", "P4_12_12  ", "     ", PRIMITIVE, 12), & !/* 369 */
        SpaceGroup(93, "D4^5  ", "P 4c 2          ", "P 4_2 2 2                      ", &
     "P 4_2 2 2          ", "P4_222    ", "     ", PRIMITIVE, 12), & !/* 370 */
        SpaceGroup(94, "D4^6  ", "P 4n 2n         ", "P 4_2 2_1 2                    ", &
     "P 4_2 2_1 2        ", "P4_22_12  ", "     ", PRIMITIVE, 12), & !/* 371 */
        SpaceGroup(95, "D4^7  ", "P 4cw 2c        ", "P 4_3 2 2                      ", &
     "P 4_3 2 2          ", "P4_322    ", "     ", PRIMITIVE, 12), & !/* 372 */
        SpaceGroup(96, "D4^8  ", "P 4nw 2abw      ", "P 4_3 2_1 2                    ", &
     "P 4_3 2_1 2        ", "P4_32_12  ", "     ", PRIMITIVE, 12), & !/* 373 */
        SpaceGroup(97, "D4^9  ", "I 4 2           ", "I 4 2 2                        ", &
     "I 4 2 2            ", "I422      ", "     ", BODY, 12), & !/* 374 */
        SpaceGroup(98, "D4^10 ", "I 4bw 2bw       ", "I 4_1 2 2                      ", &
     "I 4_1 2 2          ", "I4_122    ", "     ", BODY, 12), & !/* 375 */
        SpaceGroup(99, "C4v^1 ", "P 4 -2          ", "P 4 m m                        ", &
     "P 4 m m            ", "P4mm      ", "     ", PRIMITIVE, 13), & !/* 376 */
        SpaceGroup(100, "C4v^2 ", "P 4 -2ab        ", "P 4 b m                        ", &
     "P 4 b m            ", "P4bm      ", "     ", PRIMITIVE, 13), & !/* 377 */
        SpaceGroup(101, "C4v^3 ", "P 4c -2c        ", "P 4_2 c m                      ", &
     "P 4_2 c m          ", "P4_2cm    ", "     ", PRIMITIVE, 13), & !/* 378 */
        SpaceGroup(102, "C4v^4 ", "P 4n -2n        ", "P 4_2 n m                      ", &
     "P 4_2 n m          ", "P4_2nm    ", "     ", PRIMITIVE, 13), & !/* 379 */
        SpaceGroup(103, "C4v^5 ", "P 4 -2c         ", "P 4 c c                        ", &
     "P 4 c c            ", "P4cc      ", "     ", PRIMITIVE, 13), & !/* 380 */
        SpaceGroup(104, "C4v^6 ", "P 4 -2n         ", "P 4 n c                        ", &
     "P 4 n c            ", "P4nc      ", "     ", PRIMITIVE, 13), & !/* 381 */
        SpaceGroup(105, "C4v^7 ", "P 4c -2         ", "P 4_2 m c                      ", &
     "P 4_2 m c          ", "P4_2mc    ", "     ", PRIMITIVE, 13), & !/* 382 */
        SpaceGroup(106, "C4v^8 ", "P 4c -2ab       ", "P 4_2 b c                      ", &
     "P 4_2 b c          ", "P4_2bc    ", "     ", PRIMITIVE, 13), & !/* 383 */
        SpaceGroup(107, "C4v^9 ", "I 4 -2          ", "I 4 m m                        ", &
     "I 4 m m            ", "I4mm      ", "     ", BODY, 13), & !/* 384 */
        SpaceGroup(108, "C4v^10", "I 4 -2c         ", "I 4 c m                        ", &
     "I 4 c m            ", "I4cm      ", "     ", BODY, 13), & !/* 385 */
        SpaceGroup(109, "C4v^11", "I 4bw -2        ", "I 4_1 m d                      ", &
     "I 4_1 m d          ", "I4_1md    ", "     ", BODY, 13), & !/* 386 */
        SpaceGroup(110, "C4v^12", "I 4bw -2c       ", "I 4_1 c d                      ", &
     "I 4_1 c d          ", "I4_1cd    ", "     ", BODY, 13), & !/* 387 */
        SpaceGroup(111, "D2d^1 ", "P -4 2          ", "P -4 2 m                       ", &
     "P -4 2 m           ", "P-42m     ", "     ", PRIMITIVE, 14), & !/* 388 */
        SpaceGroup(112, "D2d^2 ", "P -4 2c         ", "P -4 2 c                       ", &
     "P -4 2 c           ", "P-42c     ", "     ", PRIMITIVE, 14), & !/* 389 */
        SpaceGroup(113, "D2d^3 ", "P -4 2ab        ", "P -4 2_1 m                     ", &
     "P -4 2_1 m         ", "P-42_1m   ", "     ", PRIMITIVE, 14), & !/* 390 */
        SpaceGroup(114, "D2d^4 ", "P -4 2n         ", "P -4 2_1 c                     ", &
     "P -4 2_1 c         ", "P-42_1c   ", "     ", PRIMITIVE, 14), & !/* 391 */
        SpaceGroup(115, "D2d^5 ", "P -4 -2         ", "P -4 m 2                       ", &
     "P -4 m 2           ", "P-4m2     ", "     ", PRIMITIVE, 14), & !/* 392 */
        SpaceGroup(116, "D2d^6 ", "P -4 -2c        ", "P -4 c 2                       ", &
     "P -4 c 2           ", "P-4c2     ", "     ", PRIMITIVE, 14), & !/* 393 */
        SpaceGroup(117, "D2d^7 ", "P -4 -2ab       ", "P -4 b 2                       ", &
     "P -4 b 2           ", "P-4b2     ", "     ", PRIMITIVE, 14), & !/* 394 */
        SpaceGroup(118, "D2d^8 ", "P -4 -2n        ", "P -4 n 2                       ", &
     "P -4 n 2           ", "P-4n2     ", "     ", PRIMITIVE, 14), & !/* 395 */
        SpaceGroup(119, "D2d^9 ", "I -4 -2         ", "I -4 m 2                       ", &
     "I -4 m 2           ", "I-4m2     ", "     ", BODY, 14), & !/* 396 */
        SpaceGroup(120, "D2d^10", "I -4 -2c        ", "I -4 c 2                       ", &
     "I -4 c 2           ", "I-4c2     ", "     ", BODY, 14), & !/* 397 */
        SpaceGroup(121, "D2d^11", "I -4 2          ", "I -4 2 m                       ", &
     "I -4 2 m           ", "I-42m     ", "     ", BODY, 14), & !/* 398 */
        SpaceGroup(122, "D2d^12", "I -4 2bw        ", "I -4 2 d                       ", &
     "I -4 2 d           ", "I-42d     ", "     ", BODY, 14), & !/* 399 */
        SpaceGroup(123, "D4h^1 ", "-P 4 2          ", "P 4/m m m                      ", &
     "P 4/m 2/m 2/m      ", "P4/mmm    ", "     ", PRIMITIVE, 15), & !/* 400 */
        SpaceGroup(124, "D4h^2 ", "-P 4 2c         ", "P 4/m c c                      ", &
     "P 4/m 2/c 2/c      ", "P4/mcc    ", "     ", PRIMITIVE, 15), & !/* 401 */
        SpaceGroup(125, "D4h^3 ", "P 4 2 -1ab      ", "P 4/n b m                      ", &
     "P 4/n 2/b 2/m      ", "P4/nbm    ", "1    ", PRIMITIVE, 15), & !/* 402 */
        SpaceGroup(125, "D4h^3 ", "-P 4a 2b        ", "P 4/n b m                      ", &
     "P 4/n 2/b 2/m      ", "P4/nbm    ", "2    ", PRIMITIVE, 15), & !/* 403 */
        SpaceGroup(126, "D4h^4 ", "P 4 2 -1n       ", "P 4/n n c                      ", &
     "P 4/n 2/n 2/c      ", "P4/nnc    ", "1    ", PRIMITIVE, 15), & !/* 404 */
        SpaceGroup(126, "D4h^4 ", "-P 4a 2bc       ", "P 4/n n c                      ", &
     "P 4/n 2/n 2/c      ", "P4/nnc    ", "2    ", PRIMITIVE, 15), & !/* 405 */
        SpaceGroup(127, "D4h^5 ", "-P 4 2ab        ", "P 4/m b m                      ", &
     "P 4/m 2_1/b m      ", "P4/mbm    ", "     ", PRIMITIVE, 15), & !/* 406 */
        SpaceGroup(128, "D4h^6 ", "-P 4 2n         ", "P 4/m n c                      ", &
     "P 4/m 2_1/n c      ", "P4/mnc    ", "     ", PRIMITIVE, 15), & !/* 407 */
        SpaceGroup(129, "D4h^7 ", "P 4ab 2ab -1ab  ", "P 4/n m m                      ", &
     "P 4/n 2_1/m m      ", "P4/nmm    ", "1    ", PRIMITIVE, 15), & !/* 408 */
        SpaceGroup(129, "D4h^7 ", "-P 4a 2a        ", "P 4/n m m                      ", &
     "P 4/n 2_1/m m      ", "P4/nmm    ", "2    ", PRIMITIVE, 15), & !/* 409 */
        SpaceGroup(130, "D4h^8 ", "P 4ab 2n -1ab   ", "P 4/n c c                      ", &
     "P 4/n 2_1/c c      ", "P4/ncc    ", "1    ", PRIMITIVE, 15), & !/* 410 */
        SpaceGroup(130, "D4h^8 ", "-P 4a 2ac       ", "P 4/n c c                      ", &
     "P 4/n 2_1/c c      ", "P4/ncc    ", "2    ", PRIMITIVE, 15), & !/* 411 */
        SpaceGroup(131, "D4h^9 ", "-P 4c 2         ", "P 4_2/m m c                    ", &
     "P 4_2/m 2/m 2/c    ", "P4_2/mmc  ", "     ", PRIMITIVE, 15), & !/* 412 */
        SpaceGroup(132, "D4h^10", "-P 4c 2c        ", "P 4_2/m c m                    ", &
     "P 4_2/m 2/c 2/m    ", "P4_2/mcm  ", "     ", PRIMITIVE, 15), & !/* 413 */
        SpaceGroup(133, "D4h^11", "P 4n 2c -1n     ", "P 4_2/n b c                    ", &
     "P 4_2/n 2/b 2/c    ", "P4_2/nbc  ", "1    ", PRIMITIVE, 15), & !/* 414 */
        SpaceGroup(133, "D4h^11", "-P 4ac 2b       ", "P 4_2/n b c                    ", &
     "P 4_2/n 2/b 2/c    ", "P4_2/nbc  ", "2    ", PRIMITIVE, 15), & !/* 415 */
        SpaceGroup(134, "D4h^12", "P 4n 2 -1n      ", "P 4_2/n n m                    ", &
     "P 4_2/n 2/n 2/m    ", "P4_2/nnm  ", "1    ", PRIMITIVE, 15), & !/* 416 */
        SpaceGroup(134, "D4h^12", "-P 4ac 2bc      ", "P 4_2/n n m                    ", &
     "P 4_2/n 2/n 2/m    ", "P4_2/nnm  ", "2    ", PRIMITIVE, 15), & !/* 417 */
        SpaceGroup(135, "D4h^13", "-P 4c 2ab       ", "P 4_2/m b c                    ", &
     "P 4_2/m 2_1/b 2/c  ", "P4_2/mbc  ", "     ", PRIMITIVE, 15), & !/* 418 */
        SpaceGroup(136, "D4h^14", "-P 4n 2n        ", "P 4_2/m n m                    ", &
     "P 4_2/m 2_1/n 2/m  ", "P4_2/mnm  ", "     ", PRIMITIVE, 15), & !/* 419 */
        SpaceGroup(137, "D4h^15", "P 4n 2n -1n     ", "P 4_2/n m c                    ", &
     "P 4_2/n 2_1/m 2/c  ", "P4_2/nmc  ", "1    ", PRIMITIVE, 15), & !/* 420 */
        SpaceGroup(137, "D4h^15", "-P 4ac 2a       ", "P 4_2/n m c                    ", &
     "P 4_2/n 2_1/m 2/c  ", "P4_2/nmc  ", "2    ", PRIMITIVE, 15), & !/* 421 */
        SpaceGroup(138, "D4h^16", "P 4n 2ab -1n    ", "P 4_2/n c m                    ", &
     "P 4_2/n 2_1/c 2/m  ", "P4_2/ncm  ", "1    ", PRIMITIVE, 15), & !/* 422 */
        SpaceGroup(138, "D4h^16", "-P 4ac 2ac      ", "P 4_2/n c m                    ", &
     "P 4_2/n 2_1/c 2/m  ", "P4_2/ncm  ", "2    ", PRIMITIVE, 15), & !/* 423 */
        SpaceGroup(139, "D4h^17", "-I 4 2          ", "I 4/m m m                      ", &
     "I 4/m 2/m 2/m      ", "I4/mmm    ", "     ", BODY, 15), & !/* 424 */
        SpaceGroup(140, "D4h^18", "-I 4 2c         ", "I 4/m c m                      ", &
     "I 4/m 2/c 2/m      ", "I4/mcm    ", "     ", BODY, 15), & !/* 425 */
        SpaceGroup(141, "D4h^19", "I 4bw 2bw -1bw  ", "I 4_1/a m d                    ", &
     "I 4_1/a 2/m 2/d    ", "I4_1/amd  ", "1    ", BODY, 15), & !/* 426 */
        SpaceGroup(141, "D4h^19", "-I 4bd 2        ", "I 4_1/a m d                    ", &
     "I 4_1/a 2/m 2/d    ", "I4_1/amd  ", "2    ", BODY, 15), & !/* 427 */
        SpaceGroup(142, "D4h^20", "I 4bw 2aw -1bw  ", "I 4_1/a c d                    ", &
     "I 4_1/a 2/c 2/d    ", "I4_1/acd  ", "1    ", BODY, 15), & !/* 428 */
        SpaceGroup(142, "D4h^20", "-I 4bd 2c       ", "I 4_1/a c d                    ", &
     "I 4_1/a 2/c 2/d    ", "I4_1/acd  ", "2    ", BODY, 15), & !/* 429 */
        SpaceGroup(143, "C3^1  ", "P 3             ", "P 3                            ", &
     "P 3                ", "P3        ", "     ", PRIMITIVE, 16), & !/* 430 */
        SpaceGroup(144, "C3^2  ", "P 31            ", "P 3_1                          ", &
     "P 3_1              ", "P3_1      ", "     ", PRIMITIVE, 16), & !/* 431 */
        SpaceGroup(145, "C3^3  ", "P 32            ", "P 3_2                          ", &
     "P 3_2              ", "P3_2      ", "     ", PRIMITIVE, 16), & !/* 432 */
        SpaceGroup(146, "C3^4  ", "R 3             ", "R 3                            ", &
     "R 3                ", "R3        ", "H    ", R_CENTER, 16), & !/* 433 */
        SpaceGroup(146, "C3^4  ", "P 3*            ", "R 3                            ", &
     "R 3                ", "R3        ", "R    ", PRIMITIVE, 16), & !/* 434 */
        SpaceGroup(147, "C3i^1 ", "-P 3            ", "P -3                           ", &
     "P -3               ", "P-3       ", "     ", PRIMITIVE, 17), & !/* 435 */
        SpaceGroup(148, "C3i^2 ", "-R 3            ", "R -3                           ", &
     "R -3               ", "R-3       ", "H    ", R_CENTER, 17), & !/* 436 */
        SpaceGroup(148, "C3i^2 ", "-P 3*           ", "R -3                           ", &
     "R -3               ", "R-3       ", "R    ", PRIMITIVE, 17), & !/* 437 */
        SpaceGroup(149, "D3^1  ", "P 3 2           ", "P 3 1 2                        ", &
     "P 3 1 2            ", "P312      ", "     ", PRIMITIVE, 18), & !/* 438 */
        SpaceGroup(150, "D3^2  ", "P 3 2=          ", "P 3 2 1                        ", &
     "P 3 2 1            ", "P321      ", "     ", PRIMITIVE, 18), & !/* 439 */
        SpaceGroup(151, "D3^3  ", "P 31 2 (0 0 4)  ", "P 3_1 1 2                      ", &
     "P 3_1 1 2          ", "P3_112    ", "     ", PRIMITIVE, 18), & !/* 440 */
        SpaceGroup(152, "D3^4  ", "P 31 2=         ", "P 3_1 2 1                      ", &
     "P 3_1 2 1          ", "P3_121    ", "     ", PRIMITIVE, 18), & !/* 441 */
        SpaceGroup(153, "D3^5  ", "P 32 2 (0 0 2)  ", "P 3_2 1 2                      ", &
     "P 3_2 1 2          ", "P3_212    ", "     ", PRIMITIVE, 18), & !/* 442 */
        SpaceGroup(154, "D3^6  ", "P 32 2=         ", "P 3_2 2 1                      ", &
     "P 3_2 2 1          ", "P3_221    ", "     ", PRIMITIVE, 18), & !/* 443 */
        SpaceGroup(155, "D3^7  ", "R 3 2=          ", "R 3 2                          ", &
     "R 3 2              ", "R32       ", "H    ", R_CENTER, 18), & !/* 444 */
        SpaceGroup(155, "D3^7  ", "P 3* 2          ", "R 3 2                          ", &
     "R 3 2              ", "R32       ", "R    ", PRIMITIVE, 18), & !/* 445 */
        SpaceGroup(156, "C3v^1 ", "P 3 -2=         ", "P 3 m 1                        ", &
     "P 3 m 1            ", "P3m1      ", "     ", PRIMITIVE, 19), & !/* 446 */
        SpaceGroup(157, "C3v^2 ", "P 3 -2          ", "P 3 1 m                        ", &
     "P 3 1 m            ", "P31m      ", "     ", PRIMITIVE, 19), & !/* 447 */
        SpaceGroup(158, "C3v^3 ", "P 3 -2=c        ", "P 3 c 1                        ", &
     "P 3 c 1            ", "P3c1      ", "     ", PRIMITIVE, 19), & !/* 448 */
        SpaceGroup(159, "C3v^4 ", "P 3 -2c         ", "P 3 1 c                        ", &
     "P 3 1 c            ", "P31c      ", "     ", PRIMITIVE, 19), & !/* 449 */
        SpaceGroup(160, "C3v^5 ", "R 3 -2=         ", "R 3 m                          ", &
     "R 3 m              ", "R3m       ", "H    ", R_CENTER, 19), & !/* 450 */
        SpaceGroup(160, "C3v^5 ", "P 3* -2         ", "R 3 m                          ", &
     "R 3 m              ", "R3m       ", "R    ", PRIMITIVE, 19), & !/* 451 */
        SpaceGroup(161, "C3v^6 ", "R 3 -2=c        ", "R 3 c                          ", &
     "R 3 c              ", "R3c       ", "H    ", R_CENTER, 19), & !/* 452 */
        SpaceGroup(161, "C3v^6 ", "P 3* -2n        ", "R 3 c                          ", &
     "R 3 c              ", "R3c       ", "R    ", PRIMITIVE, 19), & !/* 453 */
        SpaceGroup(162, "D3d^1 ", "-P 3 2          ", "P -3 1 m                       ", &
     "P -3 1 2/m         ", "P-31m     ", "     ", PRIMITIVE, 20), & !/* 454 */
        SpaceGroup(163, "D3d^2 ", "-P 3 2c         ", "P -3 1 c                       ", &
     "P -3 1 2/c         ", "P-31c     ", "     ", PRIMITIVE, 20), & !/* 455 */
        SpaceGroup(164, "D3d^3 ", "-P 3 2=         ", "P -3 m 1                       ", &
     "P -3 2/m 1         ", "P-3m1     ", "     ", PRIMITIVE, 20), & !/* 456 */
        SpaceGroup(165, "D3d^4 ", "-P 3 2=c        ", "P -3 c 1                       ", &
     "P -3 2/c 1         ", "P-3c1     ", "     ", PRIMITIVE, 20), & !/* 457 */
        SpaceGroup(166, "D3d^5 ", "-R 3 2=         ", "R -3 m                         ", &
     "R -3 2/m           ", "R-3m      ", "H    ", R_CENTER, 20), & !/* 458 */
        SpaceGroup(166, "D3d^5 ", "-P 3* 2         ", "R -3 m                         ", &
     "R -3 2/m           ", "R-3m      ", "R    ", PRIMITIVE, 20), & !/* 459 */
        SpaceGroup(167, "D3d^6 ", "-R 3 2=c        ", "R -3 c                         ", &
     "R -3 2/c           ", "R-3c      ", "H    ", R_CENTER, 20), & !/* 460 */
        SpaceGroup(167, "D3d^6 ", "-P 3* 2n        ", "R -3 c                         ", &
     "R -3 2/c           ", "R-3c      ", "R    ", PRIMITIVE, 20), & !/* 461 */
        SpaceGroup(168, "C6^1  ", "P 6             ", "P 6                            ", &
     "P 6                ", "P6        ", "     ", PRIMITIVE, 21), & !/* 462 */
        SpaceGroup(169, "C6^2  ", "P 61            ", "P 6_1                          ", &
     "P 6_1              ", "P6_1      ", "     ", PRIMITIVE, 21), & !/* 463 */
        SpaceGroup(170, "C6^3  ", "P 65            ", "P 6_5                          ", &
     "P 6_5              ", "P6_5      ", "     ", PRIMITIVE, 21), & !/* 464 */
        SpaceGroup(171, "C6^4  ", "P 62            ", "P 6_2                          ", &
     "P 6_2              ", "P6_2      ", "     ", PRIMITIVE, 21), & !/* 465 */
        SpaceGroup(172, "C6^5  ", "P 64            ", "P 6_4                          ", &
     "P 6_4              ", "P6_4      ", "     ", PRIMITIVE, 21), & !/* 466 */
        SpaceGroup(173, "C6^6  ", "P 6c            ", "P 6_3                          ", &
     "P 6_3              ", "P6_3      ", "     ", PRIMITIVE, 21), & !/* 467 */
        SpaceGroup(174, "C3h^1 ", "P -6            ", "P -6                           ", &
     "P -6               ", "P-6       ", "     ", PRIMITIVE, 22), & !/* 468 */
        SpaceGroup(175, "C6h^1 ", "-P 6            ", "P 6/m                          ", &
     "P 6/m              ", "P6/m      ", "     ", PRIMITIVE, 23), & !/* 469 */
        SpaceGroup(176, "C6h^2 ", "-P 6c           ", "P 6_3/m                        ", &
     "P 6_3/m            ", "P6_3/m    ", "     ", PRIMITIVE, 23), & !/* 470 */
        SpaceGroup(177, "D6^1  ", "P 6 2           ", "P 6 2 2                        ", &
     "P 6 2 2            ", "P622      ", "     ", PRIMITIVE, 24), & !/* 471 */
        SpaceGroup(178, "D6^2  ", "P 61 2 (0 0 5)  ", "P 6_1 2 2                      ", &
     "P 6_1 2 2          ", "P6_122    ", "     ", PRIMITIVE, 24), & !/* 472 */
        SpaceGroup(179, "D6^3  ", "P 65 2 (0 0 1)  ", "P 6_5 2 2                      ", &
     "P 6_5 2 2          ", "P6_522    ", "     ", PRIMITIVE, 24), & !/* 473 */
        SpaceGroup(180, "D6^4  ", "P 62 2 (0 0 4)  ", "P 6_2 2 2                      ", &
     "P 6_2 2 2          ", "P6_222    ", "     ", PRIMITIVE, 24), & !/* 474 */
        SpaceGroup(181, "D6^5  ", "P 64 2 (0 0 2)  ", "P 6_4 2 2                      ", &
     "P 6_4 2 2          ", "P6_422    ", "     ", PRIMITIVE, 24), & !/* 475 */
        SpaceGroup(182, "D6^6  ", "P 6c 2c         ", "P 6_3 2 2                      ", &
     "P 6_3 2 2          ", "P6_322    ", "     ", PRIMITIVE, 24), & !/* 476 */
        SpaceGroup(183, "C6v^1 ", "P 6 -2          ", "P 6 m m                        ", &
     "P 6 m m            ", "P6mm      ", "     ", PRIMITIVE, 25), & !/* 477 */
        SpaceGroup(184, "C6v^2 ", "P 6 -2c         ", "P 6 c c                        ", &
     "P 6 c c            ", "P6cc      ", "     ", PRIMITIVE, 25), & !/* 478 */
        SpaceGroup(185, "C6v^3 ", "P 6c -2         ", "P 6_3 c m                      ", &
     "P 6_3 c m          ", "P6_3cm    ", "     ", PRIMITIVE, 25), & !/* 479 */
        SpaceGroup(186, "C6v^4 ", "P 6c -2c        ", "P 6_3 m c                      ", &
     "P 6_3 m c          ", "P6_3mc    ", "     ", PRIMITIVE, 25), & !/* 480 */
        SpaceGroup(187, "D3h^1 ", "P -6 2          ", "P -6 m 2                       ", &
     "P -6 m 2           ", "P-6m2     ", "     ", PRIMITIVE, 26), & !/* 481 */
        SpaceGroup(188, "D3h^2 ", "P -6c 2         ", "P -6 c 2                       ", &
     "P -6 c 2           ", "P-6c2     ", "     ", PRIMITIVE, 26), & !/* 482 */
        SpaceGroup(189, "D3h^3 ", "P -6 -2         ", "P -6 2 m                       ", &
     "P -6 2 m           ", "P-62m     ", "     ", PRIMITIVE, 26), & !/* 483 */
        SpaceGroup(190, "D3h^4 ", "P -6c -2c       ", "P -6 2 c                       ", &
     "P -6 2 c           ", "P-62c     ", "     ", PRIMITIVE, 26), & !/* 484 */
        SpaceGroup(191, "D6h^1 ", "-P 6 2          ", "P 6/m m m                      ", &
     "P 6/m 2/m 2/m      ", "P6/mmm    ", "     ", PRIMITIVE, 27), & !/* 485 */
        SpaceGroup(192, "D6h^2 ", "-P 6 2c         ", "P 6/m c c                      ", &
     "P 6/m 2/c 2/c      ", "P6/mcc    ", "     ", PRIMITIVE, 27), & !/* 486 */
        SpaceGroup(193, "D6h^3 ", "-P 6c 2         ", "P 6_3/m c m                    ", &
     "P 6_3/m 2/c 2/m    ", "P6_3/mcm  ", "     ", PRIMITIVE, 27), & !/* 487 */
        SpaceGroup(194, "D6h^4 ", "-P 6c 2c        ", "P 6_3/m m c                    ", &
     "P 6_3/m 2/m 2/c    ", "P6_3/mmc  ", "     ", PRIMITIVE, 27), & !/* 488 */
        SpaceGroup(195, "T^1   ", "P 2 2 3         ", "P 2 3                          ", &
     "P 2 3              ", "P23       ", "     ", PRIMITIVE, 28), & !/* 489 */
        SpaceGroup(196, "T^2   ", "F 2 2 3         ", "F 2 3                          ", &
     "F 2 3              ", "F23       ", "     ", FACE, 28), & !/* 490 */
        SpaceGroup(197, "T^3   ", "I 2 2 3         ", "I 2 3                          ", &
     "I 2 3              ", "I23       ", "     ", BODY, 28), & !/* 491 */
        SpaceGroup(198, "T^4   ", "P 2ac 2ab 3     ", "P 2_1 3                        ", &
     "P 2_1 3            ", "P2_13     ", "     ", PRIMITIVE, 28), & !/* 492 */
        SpaceGroup(199, "T^5   ", "I 2b 2c 3       ", "I 2_1 3                        ", &
     "I 2_1 3            ", "I2_13     ", "     ", BODY, 28), & !/* 493 */
        SpaceGroup(200, "Th^1  ", "-P 2 2 3        ", "P m -3                         ", &
     "P 2/m -3           ", "Pm-3      ", "     ", PRIMITIVE, 29), & !/* 494 */
        SpaceGroup(201, "Th^2  ", "P 2 2 3 -1n     ", "P n -3                         ", &
     "P 2/n -3           ", "Pn-3      ", "1    ", PRIMITIVE, 29), & !/* 495 */
        SpaceGroup(201, "Th^2  ", "-P 2ab 2bc 3    ", "P n -3                         ", &
     "P 2/n -3           ", "Pn-3      ", "2    ", PRIMITIVE, 29), & !/* 496 */
        SpaceGroup(202, "Th^3  ", "-F 2 2 3        ", "F m -3                         ", &
     "F 2/m -3           ", "Fm-3      ", "     ", FACE, 29), & !/* 497 */
        SpaceGroup(203, "Th^4  ", "F 2 2 3 -1d     ", "F d -3                         ", &
     "F 2/d -3           ", "Fd-3      ", "1    ", FACE, 29), & !/* 498 */
        SpaceGroup(203, "Th^4  ", "-F 2uv 2vw 3    ", "F d -3                         ", &
     "F 2/d -3           ", "Fd-3      ", "2    ", FACE, 29), & !/* 499 */
        SpaceGroup(204, "Th^5  ", "-I 2 2 3        ", "I m -3                         ", &
     "I 2/m -3           ", "Im-3      ", "     ", BODY, 29), & !/* 500 */
        SpaceGroup(205, "Th^6  ", "-P 2ac 2ab 3    ", "P a -3                         ", &
     "P 2_1/a -3         ", "Pa-3      ", "     ", PRIMITIVE, 29), & !/* 501 */
        SpaceGroup(206, "Th^7  ", "-I 2b 2c 3      ", "I a -3                         ", &
     "I 2_1/a -3         ", "Ia-3      ", "     ", BODY, 29), & !/* 502 */
        SpaceGroup(207, "O^1   ", "P 4 2 3         ", "P 4 3 2                        ", &
     "P 4 3 2            ", "P432      ", "     ", PRIMITIVE, 30), & !/* 503 */
        SpaceGroup(208, "O^2   ", "P 4n 2 3        ", "P 4_2 3 2                      ", &
     "P 4_2 3 2          ", "P4_232    ", "     ", PRIMITIVE, 30), & !/* 504 */
        SpaceGroup(209, "O^3   ", "F 4 2 3         ", "F 4 3 2                        ", &
     "F 4 3 2            ", "F432      ", "     ", FACE, 30), & !/* 505 */
        SpaceGroup(210, "O^4   ", "F 4d 2 3        ", "F 4_1 3 2                      ", &
     "F 4_1 3 2          ", "F4_132    ", "     ", FACE, 30), & !/* 506 */
        SpaceGroup(211, "O^5   ", "I 4 2 3         ", "I 4 3 2                        ", &
     "I 4 3 2            ", "I432      ", "     ", BODY, 30), & !/* 507 */
        SpaceGroup(212, "O^6   ", "P 4acd 2ab 3    ", "P 4_3 3 2                      ", &
     "P 4_3 3 2          ", "P4_332    ", "     ", PRIMITIVE, 30), & !/* 508 */
        SpaceGroup(213, "O^7   ", "P 4bd 2ab 3     ", "P 4_1 3 2                      ", &
     "P 4_1 3 2          ", "P4_132    ", "     ", PRIMITIVE, 30), & !/* 509 */
        SpaceGroup(214, "O^8   ", "I 4bd 2c 3      ", "I 4_1 3 2                      ", &
     "I 4_1 3 2          ", "I4_132    ", "     ", BODY, 30), & !/* 510 */
        SpaceGroup(215, "Td^1  ", "P -4 2 3        ", "P -4 3 m                       ", &
     "P -4 3 m           ", "P-43m     ", "     ", PRIMITIVE, 31), & !/* 511 */
        SpaceGroup(216, "Td^2  ", "F -4 2 3        ", "F -4 3 m                       ", &
     "F -4 3 m           ", "F-43m     ", "     ", FACE, 31), & !/* 512 */
        SpaceGroup(217, "Td^3  ", "I -4 2 3        ", "I -4 3 m                       ", &
     "I -4 3 m           ", "I-43m     ", "     ", BODY, 31), & !/* 513 */
        SpaceGroup(218, "Td^4  ", "P -4n 2 3       ", "P -4 3 n                       ", &
     "P -4 3 n           ", "P-43n     ", "     ", PRIMITIVE, 31), & !/* 514 */
        SpaceGroup(219, "Td^5  ", "F -4a 2 3       ", "F -4 3 c                       ", &
     "F -4 3 c           ", "F-43c     ", "     ", FACE, 31), & !/* 515 */
        SpaceGroup(220, "Td^6  ", "I -4bd 2c 3     ", "I -4 3 d                       ", &
     "I -4 3 d           ", "I-43d     ", "     ", BODY, 31), & !/* 516 */
        SpaceGroup(221, "Oh^1  ", "-P 4 2 3        ", "P m -3 m                       ", &
     "P 4/m -3 2/m       ", "Pm-3m     ", "     ", PRIMITIVE, 32), & !/* 517 */
        SpaceGroup(222, "Oh^2  ", "P 4 2 3 -1n     ", "P n -3 n                       ", &
     "P 4/n -3 2/n       ", "Pn-3n     ", "1    ", PRIMITIVE, 32), & !/* 518 */
        SpaceGroup(222, "Oh^2  ", "-P 4a 2bc 3     ", "P n -3 n                       ", &
     "P 4/n -3 2/n       ", "Pn-3n     ", "2    ", PRIMITIVE, 32), & !/* 519 */
        SpaceGroup(223, "Oh^3  ", "-P 4n 2 3       ", "P m -3 n                       ", &
     "P 4_2/m -3 2/n     ", "Pm-3n     ", "     ", PRIMITIVE, 32), & !/* 520 */
        SpaceGroup(224, "Oh^4  ", "P 4n 2 3 -1n    ", "P n -3 m                       ", &
     "P 4_2/n -3 2/m     ", "Pn-3m     ", "1    ", PRIMITIVE, 32), & !/* 521 */
        SpaceGroup(224, "Oh^4  ", "-P 4bc 2bc 3    ", "P n -3 m                       ", &
     "P 4_2/n -3 2/m     ", "Pn-3m     ", "2    ", PRIMITIVE, 32), & !/* 522 */
        SpaceGroup(225, "Oh^5  ", "-F 4 2 3        ", "F m -3 m                       ", &
     "F 4/m -3 2/m       ", "Fm-3m     ", "     ", FACE, 32), & !/* 523 */
        SpaceGroup(226, "Oh^6  ", "-F 4a 2 3       ", "F m -3 c                       ", &
     "F 4/m -3 2/c       ", "Fm-3c     ", "     ", FACE, 32), & !/* 524 */
        SpaceGroup(227, "Oh^7  ", "F 4d 2 3 -1d    ", "F d -3 m                       ", &
     "F 4_1/d -3 2/m     ", "Fd-3m     ", "1    ", FACE, 32), & !/* 525 */
        SpaceGroup(227, "Oh^7  ", "-F 4vw 2vw 3    ", "F d -3 m                       ", &
     "F 4_1/d -3 2/m     ", "Fd-3m     ", "2    ", FACE, 32), & !/* 526 */
        SpaceGroup(228, "Oh^8  ", "F 4d 2 3 -1ad   ", "F d -3 c                       ", &
     "F 4_1/d -3 2/c     ", "Fd-3c     ", "1    ", FACE, 32), & !/* 527 */
        SpaceGroup(228, "Oh^8  ", "-F 4ud 2vw 3    ", "F d -3 c                       ", &
     "F 4_1/d -3 2/c     ", "Fd-3c     ", "2    ", FACE, 32), & !/* 528 */
        SpaceGroup(229, "Oh^9  ", "-I 4 2 3        ", "I m -3 m                       ", &
     "I 4/m -3 2/m       ", "Im-3m     ", "     ", BODY, 32), & !/* 529 */
        SpaceGroup(230, "Oh^10 ", "-I 4bd 2c 3     ", "I a -3 d                       ", &
     "I 4_1/a -3 2/d     ", "Ia-3d     ", "     ", BODY, 32) & !/* 530 */
   ]
