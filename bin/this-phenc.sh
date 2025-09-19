#!/bin/bash
# =============================================================================
# @file   this-phenc.sh
# @author Derek Anderson
# @date   09.19.2025
#
# Add PHENC include to your various include 
# paths (bash flavored).
#
# Usage: from the top-directory of the
# repo
#   source bin/this-phenc.sh
# =============================================================================

# function to add a path to another w/o duplication
add_this() {

  # function arguments
  local add="$1"
  local name="$2"

  # grab variable value
  eval "base=\${${name}:-}"

  # if directory not yet in path, add
  # add it to path
  if [[ $base != *"$add"* ]]; then
    eval "export $name=$base:$add" 
  fi

}  # end 'add_this()'

add_this $PWD/include PATH
add_this $PWD/include ROOT_INCLUDE_PATH

# end =========================================================================
