# CD_MSG_SECTION(title)
#
# print a section title to standard output
#
AC_DEFUN([CD_MSG_SECTION], [
  dnl Make sure one argument is present
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  title="$1"
  title_length=${#title}
  N=$(( 72 - title_length ))
  spaces=`yes " " | head -n $N | tr -d '\n'`
  echo ""
  echo " ┌─────────────────────────────────────────────────────────────────────────┐"
  echo " │ ${title}${spaces}│"
  echo " └─────────────────────────────────────────────────────────────────────────┘"
  echo ""
]) # CD_MSG_SECTION
