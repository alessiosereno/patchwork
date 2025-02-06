#!/bin/bash -
#===============================================================================
#
#          FILE: intstall.sh
#
#         USAGE: run "./install.sh [options]" from MOSE master directory
#
#   DESCRIPTION: A utility script that builds MOSE project
#===============================================================================

# DEBUGGING
set -e
set -C # noclobber

# INTERNAL VARIABLES AND INITIALIZATIONS
readonly PROJECT="patchwork"
readonly DIR=$(pwd)
readonly PROGRAM=`basename "$0"`

function usage () {
    echo "Install script of $PROJECT"
    echo "Usage:"
    echo
    echo "$PROGRAM --help|-?"
    echo "    Print this usage output and exit"
    echo
    echo "$PROGRAM --build  |-b"
    echo "    Build the whole project via CMake"
    echo
    echo "$PROGRAM --compile|-c <build>"
    echo "    Compile with build type (<build> = RELEASE,DEBUG,TESTING)"
    echo
    echo "$PROGRAM --setvars|-s "
    echo "    Set the project paths in the environment variables"
    echo
    echo "$PROGRAM --update |-u"
    echo "    Download the latest version of each submodule"
    echo
    echo "$PROGRAM --load |-l"
    echo "    Download the current version of each submodule"
    echo
}

function define_path () {
  rm -f .setvars.sh

  echo 'export PATCHWORKDIR='$DIR >> .setvars.sh
  echo 'function patchwork () { '$DIR'/bin/patchwork; }' >> .setvars.sh
  echo 'export -f patchwork' >> .setvars.sh

  grep -v "patchwork" $RCFILE > tmpfile && mv tmpfile $RCFILE
  echo 'source '$DIR'/.setvars.sh' >> $RCFILE
  source $RCFILE --force
}

function build_project () {
  
  rm -rf bin build && mkdir -p build
  if [[ $BUILD == standalone ]]; then
    echo 
    echo -e "\033[0;32m-- Stand-alone building \033[0m"
    echo
    git submodule update --init --recursive
    Master=None
  else
    echo
    echo -e "\033[0;32m-- Hydra building \033[0m"
    echo
    Master=hydra
  fi
  cd $DIR/build
  if command -v ifx &> /dev/null; then
      echo "Using Intel Fortran Compiler (ifx)"
      FC=ifx
  else
      echo "Intel Fortran Compiler (ifx) not found, falling back to GNU Fortran Compiler (gfortran)"
      FC=gfortran
  fi
  cmake .. -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_BUILD_TYPE=RELEASE -DMASTER=$Master
  make -j
}

function compile () {
  mkdir -p build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=$TYPE
  make
}

EXE=0
TYPE=0
SETVARS=0
UPDATE=0
LOAD=0
BUILD=0

# RETURN VALUES/EXIT STATUS CODES
readonly E_BAD_OPTION=254

# PROCESS COMMAND-LINE ARGUMENTS
if [ $# -eq 0 ]; then
  usage
  exit 0
fi

while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    # detect argument termination
    shift
    break
  fi
  case $1 in

    --build | -b )
      shift
      if (( $# > 0 )); then
        BUILD=$1
      else
        BUILD=standalone
      fi
      SETVARS=1
      ;;

    --compile | -c )
      shift
      if [ "$1"=="RELEASE" ] && [ "$1"=="DEBUG" ] && [ "$1"=="TESTING" ]; then
        TYPE="$1"
        EXE=1
      fi
      ;;

    --setvars | -s )
      shift
      SETVARS=1
      ;;

    --update | -u )
      shift
      UPDATE=1
      ;;

    --load | -l )
      shift
      LOAD=1
      ;;

    -? | --help )
      usage
      exit
      ;;

    -* )
      echo "Unrecognized option: $1" >&2
      usage
      exit $E_BAD_OPTION
      ;;

    * )
      break
      ;;
  esac
done

if [[ $SHELL == *"zsh"* ]]; then
  RCFILE=$HOME/.zshrc
elif [[ $SHELL == *"bash"* ]]; then
  RCFILE=$HOME/.bashrc
fi

if [ "$UPDATE" != "0" ]; then
  git submodule update --init --remote
elif [ "$LOAD" != "0" ]; then
  git submodule update --init
elif [[ "$BUILD" != "0" ]]; then
  build_project
elif [[ "$EXE" != "0" ]]; then
  compile
elif [ "$SETVARS" != "0" ]; then
  define_path $RCFILE
else
  usage
fi
