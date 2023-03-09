if [ -z $1 ]; then
  echo "Please inform target test application (check apps at obj/test dir)"
  exit -1
fi
../obj/test/$1; exit $?