# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/LAGAN/CMakeLists.txt
#
# CMakeLists.txt file for LAGAN.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_LAGAN CXX)
message (STATUS "Configuring apps/LAGAN")

set (SEQAN_APP_VERSION "1.1.5")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
# Order of packages is important!

find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (LAGAN LAGAN.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (LAGAN ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (LAGAN)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install LAGAN in ${PREFIX}/bin directory
install (TARGETS LAGAN
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/LAGAN for SeqAn release builds.
#install (FILES LICENSE
#               README
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/small.fasta
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (LAGAN)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------
	
# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} LAGAN CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:LAGAN")
  set (CPACK_PACKAGE_NAME "LAGAN")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "LAGAN - Alignment Free Sequence Comparison")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (LAGAN "LAGAN")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:LAGAN")

