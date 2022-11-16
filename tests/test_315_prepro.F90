
! https://www.nag.com/nagware/np/r61_doc/nag_f2008.html#AUTOTOC_3

program test

  use iso_fortran_env
  implicit none
  integer, parameter :: &
    double = REAL64, &
    simple = REAL32
  real(kind=double) :: a_vard
  real(kind=simple) :: a_vars
  character(len=:), allocatable :: &
    c_ver, c_opt, s_dat, s_time
  character(*), parameter :: &
    fma = '(a)', &
    fmesl = '(es23.15)', &  ! float long
    fmess = '(es12.5)', &   ! float short
    nl = achar(10), & ! NEW_LINE('A')
    tab = achar(9), &
    NULL = achar(0)   ! used as non visible character, non empty character(len=:), allocatable , and NULL as end of string C
  character(*), parameter :: banner = &
                             "              _ _           _       "//nl// &
                             "             (_) |         | |      "//nl// &
                             "    _ __ ___  _| | __ _  __| |_   _ "//nl// &
                             "   | '_ ` _ \| | |/ _` |/ _` | | | |"//nl// &
                             "   | | | | | | | | (_| | (_| | |_| |"//nl// &
                             "   |_| |_| |_|_|_|\__,_|\__,_|\__, |"//nl// &
                             "                               __/ |"//nl// &
                             "                              |___/ "//nl//nl// &
                             "   This code is under copyright & licence."//nl// &
                             "   The distribution of the package or parts of package is not allowed."//nl// &
                             "   For more details, please contact: mihai-cosmin.marinica@cea.fr"

  character(*), parameter :: banner_cosmin = &
                             " ,---.    ,---.-./`)   .---.       ____    ______        ____     __  "//nl// &
                             " |    \  /    \ .-.')  | ,_|     .'  __ `.|    _ `''.    \   \   /  / "//nl// &
                             " |  ,  \/  ,  / `-' \,-./  )    /   '  \  \ _ | ) _  \    \  _. /  '  "//nl// &
                             " |  |\_   /|  |`-'`'`\  '_ '`)  |___|  /  |( ''_'  ) |     _( )_ .'   "//nl// &
                             " |  _( )_/ |  |.---.  > (_)  )     _.-`   | . (_) `. | ___(_ o _)'    "//nl// &
                             " | (_ o _) |  ||   | (  .  .-'  .'   _    |(_    ._) '|   |(_,_)'     "//nl// &
                             " |  (_,_)  |  ||   |  `-'`-'|___|  _( )_  |  (_.\.' / |   `-'  /      "//nl// &
                             " |  |      |  ||   |   |        \ (_ o _) /       .'   \      /       "//nl// &
                             " '--'      '--''---'   `--------`'.(_,_).''-----'`      `-..-'        "//nl


  a_vard = 1/3d0
  a_vars = 1/3d0
  c_ver = compiler_version()
  c_opt = compiler_options()
  s_dat = __DATE__
  s_time = __TIME__


  write (6, fma) nl//banner//nl
  ! write (6, fma) nl//banner_cosmin//nl
  write (6, fma) nl//"compiled by "//c_ver
  write (6, fma) "compilation at "//__DATE__//" "//__TIME__//nl
  write (6, fma) nl//"compiler options was"//nl//c_opt
  write (6, "(a, es24.16)") nl//"a simple variable ", a_vars  ! to see limits 8 decimal digits
  write (6, "(a, es24.16)") nl//"a double variable ", a_vard  ! to see limits 16 decimal digits

end program test

