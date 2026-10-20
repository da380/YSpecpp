# Sourced by the scripts that run the examples, from their own directories.
#
# The programs are found in YSPEC_BIN if it is set, else in the build tree
# beside the examples (they are copied to build/examples), else on the PATH.
# The calculations use four threads unless OMP_NUM_THREADS says otherwise.

# In the build tree, where the sources are.
if [ -f ../source.sh ]; then
  . ../source.sh
fi

: "${OMP_NUM_THREADS:=4}"
export OMP_NUM_THREADS

program() {
  for directory in "${YSPEC_BIN:-}" ../../bin ../../build/bin; do
    if [ -n "$directory" ] && [ -x "$directory/$1" ]; then
      echo "$directory/$1"
      return
    fi
  done
  if command -v "$1" > /dev/null; then
    command -v "$1"
    return
  fi
  echo "cannot find $1: build YSpec, or set YSPEC_BIN to where it is" >&2
  exit 1
}

# Runs a program on an input file.
run() {
  name=$1
  shift
  executable=$(program "$name")
  echo "$name $*"
  "$executable" "$@"
}

# Plots with yspec-plot: from the PATH if it is there (in an activated
# environment, say), else through the poetry project in plotting/, found from
# YSPEC_PLOTTING, from the sources the build was made from, or beside these
# examples.
plot() {
  if command -v yspec-plot > /dev/null; then
    yspec-plot "$@"
    return
  fi
  project=${YSPEC_PLOTTING:-${YSPEC_SOURCE_DIR:-../..}/plotting}
  if command -v poetry > /dev/null &&
    poetry -P "$project" run yspec-plot --help > /dev/null 2>&1; then
    poetry -P "$project" run yspec-plot "$@"
  else
    echo "not plotting: set up the plotting environment with" \
      "'poetry -P $project install'" >&2
  fi
}
