#!/usr/bin/env bash

# common milady bash utilities, as bash functions
# cd .../MILADY/scripts && source scripts/utils_milady.bash


function f_pause {
  read -p 'Continue? (ctrl-c to stop)' rep
}

function in_red () {
  # Red bold echo
  echo -e "\e[31m\e[1m"$@"\e[0m"
}

function in_red_pause () {
  # Red bold
  echo -e "\e[31m\e[1m"$@"\e[0m"
  f_pause
}

function in_green () {
  # Green bold echo
  echo -e "\e[32m\e[1m"$@"\e[0m"
}

function f_which {
  # which returning 'UNKNOWN' if not found
  # avoid stderr message if not found
  tmp=$(which ${1} 2> /dev/null) && echo ${tmp} || echo "UNKNOWN"
}


function f_nproc {
  # policy is milady compilation make -jx
  # x max ix 9 (cosmin choice 2022)
  # 4 ok for small desktop in 2022 if not nprocs present
  tmp=$(nproc 2> /dev/null) || tmp=4
  [[ "9" -gt ${tmp} ]] && echo ${tmp} || echo "9"
}

function f_cmake3 {
  tmp=$(\cmake --version)    # "\" avoid alias, caseof
  case "$tmp" in
    *'version 2'*)
      echo cmake3 ;;
    *'version 3'*)
      echo cmake ;;
    *)
      echo "UNKNOWN cmake" ;;
  esac
}
