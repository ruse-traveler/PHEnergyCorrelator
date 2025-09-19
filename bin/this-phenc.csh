#!/bin/csh
# =============================================================================
# @file   this-phenc.csh
# @author Derek Anderson
# @date   09.19.2025
#
# Add PHENC include to your various include 
# paths (csh flavored).
#
# Usage: from the top-directory of the
# repo
#   source bin/this-phenc.csh
# =============================================================================

# FIXME this probably isn't gonna work
alias add_this '
  set add   = "$1"
  set name  = "$2"

  # Get the current value of the variable (empty if unset)
  if ( $?${name} ) then
    set base = "${${name}}"
  else
    set base = ""
  endif

  # Check whether $add is already present in $base.
  # We need to treat the value as a colon‑separated list.
  if ( "$base" !~ "*:${add}:*" && \
       "$base" !~ "*:${add}"   && \
       "$base" !~ "${add}:*"   && \
       "$base" != "${add}" ) then

    # Prepend the new directory (change to "${base}:${add}" to append)
    if ( "$base" == "" ) then
      set ${name} = "${add}"
    else
      set ${name} = "${add}:${base}"
    endif
     ${name}

  endif
'

add_this "$PWD/include" PATH
add_this "$PWD/include" ROOT_INCLUDE_PATH

# end =========================================================================
