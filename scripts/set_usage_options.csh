#!/bin/csh -v

#-----------------------------------------------------------------------
# script describes and defines usage options based on arguments given
# for zetac_run.csh
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# describe usage options
#-----------------------------------------------------------------------

  if ( $#argv == 0 ) then
    echo
    echo "zetac_run.csh [-c ] [-C ] [-D] [-G] [-mx=...] [-ny=...] [-M] [-o=...] [-q] [-s=...]"
    echo "              [-t=...] [-u] experiment_name"
    echo
    echo " -c                execute coldstart"
    echo " -C                generate checksums"
    echo " -D                turn off optimizations"
    echo " -G                re-generate grid_spec files for given .../<exp>"
    echo " -mx=              number of processors in x; overrides value in <exp/>resources file"
    echo " -ny=              number of processors in y; overrides value in <exp>/resources file"
    echo " -M                re-generate makefiles"
    echo " -o=               location of stdout used in batch submissions; default will be .../<exp>/jobs"
    echo " -q                execute in batch mode"
    echo " -s=               number of run segments; overrides value in .../<exp>/resources file"
    echo " -t=               time limit in hh:mm:ss; overrides value in .../<exp>/resources file"
    echo " experiment_name   name of experiment"
    exit
  endif

#-----------------------------------------------------------------------
# define usage options
#-----------------------------------------------------------------------
  
  set index = 1
  set limit = $#argv

#-----------------------------------------------------------------------
# define defaults
#-----------------------------------------------------------------------

  setenv progname               "zetac" 
  setenv check_tag              "false"
  setenv re_generate_grid       "false"
  setenv mx                         "0"
  setenv ny                         "0"
  setenv re_generate_make       "false"
  setenv environ          "INTERACTIVE"
  setenv segments                   "0"
  setenv resources_time             "0"
  
  while ( $index <= $limit )
  
    setenv argument "$argv[$index]"
    source $RUNDIR/scripts/set_arguments.csh 

#-----------------------------------------------------------------------
# define program name
#----------------------------------------------------------------------- 

    if ( "$argument" == "-c" ) then
      setenv progname "zetac_coldstart"
    endif

#-----------------------------------------------------------------------
# define checksums
#----------------------------------------------------------------------- 

    if ( "$argument" == "-C" ) then
      setenv check_tag .true.
    endif

#-----------------------------------------------------------------------
# define regenerate grid data
#----------------------------------------------------------------------- 

    if ( "$argument" == "-G" ) then
      setenv re_generate_grid "true"
    endif

#-----------------------------------------------------------------------
# define number of processors in x-direction
#----------------------------------------------------------------------- 

    if ( "$argument" == "-mx" ) then
      setenv mx "$variable"
    endif

#-----------------------------------------------------------------------
# define number of processors in y-direction
#----------------------------------------------------------------------- 

    if ( "$argument" == "-ny" ) then
      setenv ny "$variable"
    endif

#-----------------------------------------------------------------------
# define re-generate makefiles
#----------------------------------------------------------------------- 

    if ( "$argument" == "-M" ) then
      setenv re_generate_make "true"
    endif

#-----------------------------------------------------------------------
# define name of stdout in batch submissions file
#----------------------------------------------------------------------- 

    if ( "$argument" == "-o") then
      setenv STDOUT "$variable"
    endif

#-----------------------------------------------------------------------
# define interactive or batch mode
#----------------------------------------------------------------------- 

    if ( "$argument" == "-q" ) then
      setenv environ "BATCH"
    endif

#-----------------------------------------------------------------------
# define number of segments
#----------------------------------------------------------------------- 

    if ( "$argument" == "-s" ) then
      setenv segments "$variable"
    endif
  
#-----------------------------------------------------------------------
# define resources time-limit
#----------------------------------------------------------------------- 

    if ( "$argument" == "-t" ) then
      setenv resources_time "$variable"
    endif

    @ index++

  end

exit


