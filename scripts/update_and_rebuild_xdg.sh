#!/usr/bin/env bash
#* This file is part of the MOOSE framework
#* https://mooseframework.inl.gov
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for i in "$@"
do
  shift
  if [[ "$i" == "-h" || "$i" == "--help" ]]; then
    help=1;
  fi

  if [ "$i" == "--fast" ]; then
    go_fast=1;
  fi

  if [ "$i" == "--skip-submodule-update" ]; then
    skip_sub_update=1
  else # Remove everything else before passing to cmake
    set -- "$@" "$i"
  fi
done

# Display help
if [ -n "$help" ]; then
  echo "Usage: $0 [-h | --help | --fast | --skip-submodule-update | <xdg cmake options> ]"
  echo
  echo "-h | --help              Display this message and list of available xdg options"
  echo "--fast                   Run XDG 'make' only, do NOT run CMake"
  echo "--skip-submodule-update  Do not update the XDG submodule, use the current version"
  echo
  echo "Influential variables"
  echo "EMBREE_DIR               Path to an Embree (v4) install; default: /opt/embree"
  echo "LIBMESH_DIR              Path to libmesh (for libmesh.pc); default: ../libmesh/installed"
  echo "METHOD                   libMesh build method used to pick libmesh.pc; default: \$METHOD or opt"
  echo "XDG_DIR                  XDG install prefix; default: ../framework/contrib/xdg/installed"
  echo "XDG_SRC_DIR              Path to XDG source; default: ../framework/contrib/xdg from submodule"
  exit 0
fi

if [[ -n "$go_fast" && $# != 1 ]]; then
  echo "Error: --fast can only be used by itself or with --skip-submodule-update."
  echo "Try again, removing either --fast or all other conflicting arguments!"
  exit 1;
fi

set -e

if [ -n "$XDG_SRC_DIR" ]; then
  skip_sub_update=1
else
  XDG_SRC_DIR=$(realpath "${SCRIPT_DIR}/../framework/contrib/xdg/.")
fi
XDG_BUILD_DIR_BASE="${XDG_SRC_DIR}/build"
if [ -n "$XDG_DIR" ]; then
  echo "INFO: XDG_DIR set - overriding default installed path"
  echo "INFO: No cleaning will be done in specified path"
else
  XDG_DIR="${XDG_SRC_DIR}/installed"
  rm -rf "$XDG_DIR"
fi

# Embree is a runtime dependency of the XDG ray tracer and is not bundled with
# XDG or MOOSE, so an existing install is required
: ${EMBREE_DIR:=/opt/embree}
if [ ! -d "$EMBREE_DIR" ]; then
  echo "Error: Could not find an Embree install at \$EMBREE_DIR=$EMBREE_DIR."
  echo "Install Embree v4 (e.g. 4.4.1) and pass its prefix via EMBREE_DIR."
  exit 1
fi

: ${LIBMESH_DIR:=$(realpath "${SCRIPT_DIR}/../libmesh/installed/.")}
if [ ! -d "$LIBMESH_DIR/lib/pkgconfig" ]; then
  echo "Error: Could not find libmesh pkg-config files at \$LIBMESH_DIR/lib/pkgconfig."
  echo "Build and install libmesh first, or set LIBMESH_DIR to its install prefix."
  exit 1
fi

: ${METHOD:=${METHOD:-opt}}

if [ -z "$skip_sub_update" ]; then
  git submodule update --init --checkout "${XDG_SRC_DIR}"
fi

# XDG's CMake only accepts Embree versions in [4.0.0, 4.1.0). Widen the range so
# newer Embree 4 releases (e.g. 4.4.1) are used. See the comment in
# CMakeLists.txt about Embree's lack of version-range support.
if ! grep -q 'find_package(embree 4.0.0\.\.\.<5\.0\.0' "$XDG_SRC_DIR/CMakeLists.txt"; then
  sed -i 's/find_package(embree 4.0.0\.\.\.<4\.1\.0/find_package(embree 4.0.0...<5.0.0/' "$XDG_SRC_DIR/CMakeLists.txt"
fi

# If we're not going fast, remove the build directory and reconfigure
if [ -z "$go_fast" ]; then
  rm -rf "$XDG_BUILD_DIR_BASE"
  mkdir -p "$XDG_BUILD_DIR_BASE"
  cd "$XDG_BUILD_DIR_BASE"

  PKG_CONFIG_PATH="$LIBMESH_DIR/lib/pkgconfig" \
  METHOD="$METHOD" \
  cmake "$XDG_SRC_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$XDG_DIR" \
    -DCMAKE_PREFIX_PATH="$EMBREE_DIR" \
    -DXDG_BUILD_TESTS=OFF \
    -DXDG_BUILD_TOOLS=OFF \
    -DXDG_ENABLE_EMBREE=ON \
    -DXDG_ENABLE_GPRT=OFF \
    -DXDG_ENABLE_LIBMESH=ON \
    -DXDG_ENABLE_MOAB=OFF \
    "$@"
fi

cd "$XDG_BUILD_DIR_BASE" ||
{ echo "Error: Need to run this script without --fast at least once."; exit 1; }

make -j ${MOOSE_JOBS:-4} install