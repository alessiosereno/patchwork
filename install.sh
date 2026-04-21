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

function require_cmd () {
  # require_cmd <command> <hint>
  if ! command -v "$1" &> /dev/null; then
    echo -e "\033[0;31m-- ERROR: required command '$1' not found.\033[0m"
    echo "   $2"
    return 1
  fi
}

function check_prerequisites () {
  echo
  echo -e "\033[0;32m-- Checking prerequisites \033[0m"
  local missing=0
  require_cmd git    "Install git (e.g. 'brew install git' on macOS, 'apt install git' on Debian/Ubuntu)." || missing=1
  require_cmd cmake  "Install cmake (e.g. 'brew install cmake', 'apt install cmake')."                  || missing=1
  require_cmd make   "Install make (Xcode Command Line Tools on macOS, 'apt install build-essential')." || missing=1
  require_cmd python3 "Install Python 3 (e.g. 'brew install python', 'apt install python3')."           || missing=1

  if command -v python3 &> /dev/null; then
    if ! python3 -c "import venv" &> /dev/null; then
      echo -e "\033[0;31m-- ERROR: Python 'venv' module not available.\033[0m"
      echo "   On Debian/Ubuntu: 'apt install python3-venv'."
      missing=1
    fi
  fi

  if ! command -v ifx &> /dev/null && ! command -v gfortran &> /dev/null; then
    echo -e "\033[0;31m-- ERROR: no Fortran compiler found (neither 'ifx' nor 'gfortran').\033[0m"
    echo "   Install GNU Fortran: 'brew install gcc' on macOS, 'apt install gfortran' on Debian/Ubuntu."
    echo "   Or install the Intel oneAPI Fortran compiler (ifx)."
    missing=1
  fi

  if [ "$missing" -ne 0 ]; then
    echo
    echo -e "\033[0;31m-- Prerequisites check failed. Please install missing tools and retry.\033[0m"
    exit 1
  fi
  echo -e "\033[0;32m-- All prerequisites found \033[0m"
}

function define_path () {
  rm -f .setvars.sh

  echo 'export PATCHWORKDIR='$DIR >> .setvars.sh
  echo 'function patchwork () { '$DIR'/PATCHWORK.sh; }' >> .setvars.sh
  echo 'export -f patchwork' >> .setvars.sh

  grep -v "patchwork" $RCFILE > tmpfile && mv tmpfile $RCFILE
  echo 'source '$DIR'/.setvars.sh' >> $RCFILE
  source $RCFILE --force
}

function setup_python () {
  echo
  echo -e "\033[0;32m-- Setting up Python virtual environment \033[0m"
  echo
  python3 -m venv $DIR/.venv
  $DIR/.venv/bin/pip install --upgrade pip
  $DIR/.venv/bin/pip install -r $DIR/requirements.txt
}

function build_project () {

  check_prerequisites
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
  cd $DIR
  setup_python
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
