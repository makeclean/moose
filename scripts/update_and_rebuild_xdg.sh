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
  echo "XDG_DEPS_DIR             Root directory for auto-built XDG dependencies; default: <moose>/framework/contrib/xdg-deps"
  echo "TBB_DIR                  Path to an oneTBB install; default: \$XDG_DEPS_DIR/oneTBB (built from source if missing)"
  echo "EMBREE_DIR               Path to an Embree (v4) install; default: \$XDG_DEPS_DIR/embree (built from source if missing)"
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

# Pinned versions for the auto-built dependencies. oneTBB v2023.0.0 is the
# release line Embree 4.4 was validated against; Embree v4.4.1 is the newest
# 4.x release and satisfies the widened [4.0.0, 5.0.0) range in XDG's
# CMakeLists.txt.
ONE_TBB_VERSION=v2023.0.0
EMBREE_VERSION=v4.4.1

# Build and install oneTBB (second-order XDG dependency: Embree's tasking
# library) into TBB_DIR; skipped if an install is already present there.
build_oneTBB() {
  if [ -d "$TBB_DIR/lib" ] || [ -d "$TBB_DIR/lib64" ]; then
    echo "INFO: Using existing oneTBB install at $TBB_DIR"
    return 0
  fi

  echo "INFO: Building oneTBB $ONE_TBB_VERSION into $TBB_DIR"
  local src="$XDG_DEPS_DIR/src/oneTBB"
  local bld="$XDG_DEPS_DIR/bld/oneTBB"
  if [ ! -d "$src" ]; then
    mkdir -p "$XDG_DEPS_DIR/src"
    git clone --depth 1 --branch "$ONE_TBB_VERSION" https://github.com/uxlfoundation/oneTBB "$src"
  fi
  rm -rf "$bld"
  mkdir -p "$bld"
  cd "$bld"
  # Force lib/ so the layout is identical on every platform
  # (GNUInstallDirs defaults to lib64 on RHEL-based systems)
  cmake "$src" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$TBB_DIR" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DTBB_TEST=OFF
  make -j ${MOOSE_JOBS:-4} install
  cd "$SCRIPT_DIR"
}

# Build and install Embree (first-order XDG dependency: the ray tracing
# backend) into EMBREE_DIR, linked against the oneTBB above; skipped if an
# install is already present there.
build_embree() {
  if [ -d "$EMBREE_DIR/lib" ] || [ -d "$EMBREE_DIR/lib64" ]; then
    echo "INFO: Using existing Embree install at $EMBREE_DIR"
    return 0
  fi

  echo "INFO: Building Embree $EMBREE_VERSION into $EMBREE_DIR"
  local src="$XDG_DEPS_DIR/src/embree"
  local bld="$XDG_DEPS_DIR/bld/embree"
  if [ ! -d "$src" ]; then
    mkdir -p "$XDG_DEPS_DIR/src"
    git clone --depth 1 --branch "$EMBREE_VERSION" https://github.com/embree/embree "$src"
  fi

  rm -rf "$bld"
  mkdir -p "$bld"
  cd "$bld"
  # EMBREE_TBB_ROOT steers Embree's TBB discovery at the oneTBB above (config
  # mode, with a module-mode fallback). Force lib/ so the layout is identical
  # on every platform (GNUInstallDirs defaults to lib64 on RHEL-based systems)
  cmake "$src" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$EMBREE_DIR" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DEMBREE_TBB_ROOT="$TBB_DIR" \
    -DEMBREE_ISPC_SUPPORT=OFF \
    -DEMBREE_TUTORIALS=OFF
  make -j ${MOOSE_JOBS:-4} install
  cd "$SCRIPT_DIR"

  # Embree sets the installed library's RUNPATH to $ORIGIN, expecting libtbb
  # to sit in the same install folder; symlink it in so libembree finds
  # oneTBB at runtime without LD_LIBRARY_PATH
  local embree_lib tbb_lib
  for d in "$EMBREE_DIR/lib" "$EMBREE_DIR/lib64"; do
    if [ -d "$d" ]; then
      embree_lib="$d"
      break
    fi
  done
  # oneTBB installs into lib64 on some platforms and lib on others
  for d in "$TBB_DIR/lib" "$TBB_DIR/lib64"; do
    if [ -d "$d" ]; then
      tbb_lib="$d"
      break
    fi
  done
  if [ -z "$tbb_lib" ] || [ -z "$embree_lib" ]; then
    echo "Error: Could not locate the oneTBB or Embree library directory."
    exit 1
  fi
  # libtbb.so.12 is the SONAME of oneTBB 2023.x
  ln -sf "$tbb_lib/libtbb.so.12" "$embree_lib/libtbb.so.12"
  ln -sf "$tbb_lib/libtbb.so" "$embree_lib/libtbb.so"
}

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

# XDG's dependency chain is XDG -> Embree -> oneTBB; neither Embree nor
# oneTBB is bundled with XDG or MOOSE, so each missing level is built from
# source and installed under XDG_DEPS_DIR unless an existing install is
# provided via TBB_DIR / EMBREE_DIR.
: ${XDG_DEPS_DIR:=$(realpath -m "${SCRIPT_DIR}/../framework/contrib/xdg-deps/.")}
: ${TBB_DIR:="${XDG_DEPS_DIR}/oneTBB"}
: ${EMBREE_DIR:="${XDG_DEPS_DIR}/embree"}

build_oneTBB
build_embree

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

  # MOOSE (moose.mk) looks for libxdg.so under $(XDG_DIR)/lib;
  # GNUInstallDirs defaults to lib64 on RHEL-based systems
  PKG_CONFIG_PATH="$LIBMESH_DIR/lib/pkgconfig" \
  METHOD="$METHOD" \
  cmake "$XDG_SRC_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$XDG_DIR" \
    -DCMAKE_INSTALL_LIBDIR=lib \
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