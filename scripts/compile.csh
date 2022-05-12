#!/bin/csh -f

  set expname = $argv[$#argv]

  setenv RUNDIR $cwd

  setenv EXEDIR $RUNDIR/exec
  setenv JOBDIR $RUNDIR/jobs/$expname
  setenv SRCDIR $RUNDIR/src
  setenv EXPDIR $RUNDIR/exp/$expname
  setenv ARCDIR /lustre/f2/dev/$USER/archive
  setenv TMPDIR /lustre/f2/scratch/$USER

#####################################################################
# module selection
#####################################################################

  source ./env.cshrc.intel
# set usage options based on arguments to zetac_run.csh
#####################################################################

  source $RUNDIR/scripts/set_usage_options.csh

#####################################################################
# check usage arguments
#####################################################################

  if ( $expname == '' ) goto USAGE

  if (! -d $EXPDIR ) goto NO_EXPERIMENT

  set cppDefs = ( -Duse_libMPI -Duse_netCDF -DLARGE_FILE )
  set srcList = ( shared/mpp/include_sma shared/mpp/include_nocommunication \
        shared/mpp/include_mpi shared/mpp/include_common shared/mpp/include )

#####################################################################
# check for existence of source directory
#####################################################################  

  if ( ! -d $SRCDIR ) goto NO_SOURCE

  cd $SRCDIR

#####################################################################  
# regenerate path_names files
#####################################################################  

  rm -f path_names*
  list_paths *

  sed -e /zetac_coldstart/d  path_names > tmp
  sed -e /zetac_write_time/d tmp        > path_names_zetac

  sed -e /coupler_main/d     path_names > tmp
  sed -e /zetac_write_time/d tmp        > path_names_zetac_coldstart

  sed -e /coupler_main/d     path_names > tmp
  sed -e /zetac_coldstart/d  tmp        > path_names_zetac_time

  rm -f tmp

  set template = $RUNDIR/etc/intel.mk

  alias make  "make -j 4"

#####################################################################
# when required create makefile for zetac
#####################################################################  

  cd $EXEDIR
    
  if ( $progname == "zetac" ) then
  
    if (! -e Makefile_zetac || $re_generate_make == "true" ) then
    
      echo ; echo " ...creating makefile for zetac "

      mkmf -a $SRCDIR -m Makefile_zetac -t $template \
              -p zetac.x -c "$cppDefs" $srcList \
              $SRCDIR/path_names_zetac
      setenv error $status
      if ( $error != 0 ) goto ZETAC_MAKE
      
      sed -e "s/atmos_zetac\/user/atmos_zetac\/.user/" \
                                       $SRCDIR/path_names_zetac > tmp
      cp tmp $SRCDIR/path_names_zetac
      
      sed -e "s/atmos_zetac\/user/atmos_zetac\/.user/" \
                                                 Makefile_zetac > tmp
      cp tmp Makefile_zetac
      
      rm tmp

    endif

#####################################################################
# compiling and linking zetac.x
#####################################################################  

    echo ; echo " ...compiling and linking zetac.x "

    make DEBUG=ON -f Makefile_zetac
#    make -f Makefile_zetac
    setenv error $status
    if ( $error != 0 ) goto ZETAC_EXEC

#####################################################################
# when required create makefile for zetac_coldstart
#####################################################################  

  else if ( $progname == "zetac_coldstart" ) then

    if (! -e Makefile_zetac_coldstart || $re_generate_make == "true" ) then

      echo ; echo " ...creating makefile for zetac_coldstart "

      mkmf -a $SRCDIR -m Makefile_zetac_coldstart -t $template \
              -p zetac_coldstart.x -c "$cppDefs" $srcList \
              $SRCDIR/path_names_zetac_coldstart
      setenv error $status
      if ( $error != 0 ) goto ZETAC_COLDSTART_MAKE

      sed -e "s/atmos_zetac\/user/atmos_zetac\/.user/" \
                             $SRCDIR/path_names_zetac_coldstart > tmp
      cp tmp $SRCDIR/path_names_zetac_coldstart

      sed -e "s/atmos_zetac\/user/atmos_zetac\/.user/" \
                                       Makefile_zetac_coldstart > tmp
      cp tmp Makefile_zetac_coldstart
      
      rm tmp

    endif

#####################################################################
# copy over user code from experiment directory to source directory
# unless otherwise specified
#####################################################################  

    if ( ! -e $EXPDIR/USER ) then
       setenv error 1
       goto NO_EXPERIMENT_USER
    else 
       if ( ! -e $SRCDIR/atmos_zetac/.user ) mkdir $SRCDIR/atmos_zetac/.user
       cp $EXPDIR/USER/*.f90 $SRCDIR/atmos_zetac/.user/
    endif

#####################################################################
# compiling and linking zetac_coldstart.x
#####################################################################  

    echo " ...compiling and linking zetac_coldstart.x "

    make -f Makefile_zetac_coldstart

    setenv error $status
    if ( $error != 0 ) goto ZETAC_COLDSTART_EXEC

  else

    exit

  endif

#####################################################################
# when required create makefile for write_time
#####################################################################

  if (! -e Makefile_write_time || $re_generate_make == "true" ) then

    echo ; echo " ...creating makefile for write_time "

    mkmf -a $SRCDIR -m Makefile_write_time -t $template \
           -p zetac_write_time.x -c "$cppDefs" $srcList \
          $SRCDIR/path_names_zetac_time

    setenv error $status
    if ( $error != 0 ) goto WRITE_TIME_MAKE

      sed -e "s/atmos_zetac\/user/atmos_zetac\/.user/" \
                                  $SRCDIR/path_names_zetac_time > tmp
      cp tmp $SRCDIR/path_names_zetac_time

      sed -e "s/atmos_zetac\/user/atmos_zetac\/.user/" \
                                            Makefile_write_time > tmp
      cp tmp Makefile_write_time
      
      rm tmp
 
   endif

#####################################################################
# compiling and linking zetac_write_time.x
#####################################################################  

  echo ; echo " ...compiling and linking zetac_write_time.x "

  make -f Makefile_write_time

  setenv error $status
  if ( $error != 0 ) goto WRITE_TIME_EXEC

#####################################################################
# set identity stamp
#####################################################################

  set identity = `date +%b%dh%Hm%Ms%S`
  set identity = $identity

#####################################################################
# location of stdout in batch submissions
#####################################################################

  setenv STDOUT $EXPDIR/JOBS/
  if ( ! -d $STDOUT ) mkdir -p $STDOUT

#####################################################################
# generate run script
#####################################################################

  echo " ...generating run script "

  setenv WRKDIR $TMPDIR/$identity

  if ( ! -d $WRKDIR ) mkdir -p $WRKDIR

  cd $WRKDIR

  set basname = `echo $RUNDIR | sed 's/\(\/\)/\\\//g'`
  set tmpname = `echo $TMPDIR | sed 's/\(\/\)/\\\//g'`
  set arcname = `echo $ARCDIR | sed 's/\(\/\)/\\\//g'`
  set outname = `echo $STDOUT | sed 's/\(\/\)/\\\//g'`

  set runscript = ${progname}.csh
  set template = ${progname}.template

echo "\
  s/BBBBBBB/$basname/g\
  s/CCCCCCC/$expname/g\
  s/DDDDDDD/$progname/g\
  s/GGGGGGG/$outname/g\
  s/HHHHHHH/$identity/g\
  s/IIIIIII/$num_proc_x/g\
  s/JJJJJJJ/$num_proc_y/g\
  s/KKKKKKK/$num_proc/g\
  s/LLLLLLL/$num_segs/g\
  s/TTTTTTT/$time_limit/g" \
  > sed_file

  cat $RUNDIR/bin/$template | sed -f sed_file > $RUNDIR/bin/$runscript

#####################################################################
# copy files to jobs directory
#####################################################################

  if (! -d $JOBDIR/$identity/bin )                   \
      mkdir -p $JOBDIR/$identity/bin

  if (! -d $JOBDIR/$identity/exp )                   \
      mkdir -p $JOBDIR/$identity/exp 

  if (! -d $JOBDIR/$identity/exec )                  \
      mkdir -p $JOBDIR/$identity/exec

  cp $RUNDIR/bin/$runscript $JOBDIR/$identity/bin
  cp $EXPDIR/*.nml          $JOBDIR/$identity/exp
  cp $EXPDIR/*_table        $JOBDIR/$identity/exp
  cp $EXPDIR/resources      $JOBDIR/$identity/exp

  gzip -cf $EXEDIR/${progname}.x >                   \
     $JOBDIR/$identity/exec/${progname}.x.gz

  gzip -cf $EXEDIR/zetac_write_time.x >              \
     $JOBDIR/$identity/exec/zetac_write_time.x.gz

#####################################################################
# batch run script submission
#####################################################################

  chmod +x $RUNDIR/bin/$runscript
  sbatch $RUNDIR/bin/$runscript
  setenv error $status

#####################################################################
# generate runscript for data long-term storage
#####################################################################

  echo " ...generating auxiliary run script "

  set runscript = ${progname}_store.csh
  set template = ${progname}_store.template

echo "\
  s/CCCCCC/$expname/g\
  s/GGGGGG/$outname/g\
  s/HHHHHH/$identity/g" \
  > sed_file

  cat $RUNDIR/bin/$template | sed -f sed_file > $JOBDIR/$identity/bin/$runscript
 
#####################################################################
# experiment completed successfully
#####################################################################

  if ( $error == 0 ) then
    echo ; echo " script exited normally "
  else
    echo ; echo " script exited abnormally "
  endif  

  exit

#####################################################################
# usage
#####################################################################

USAGE:
  echo
  exit

#####################################################################
# no valid experiment name
#####################################################################

NO_EXPERIMENT:
  echo ; echo " ERROR: must specify valid experiment name "
  exit

#####################################################################
# no valid experiment directory
#####################################################################

NO_EXPERIMENT_USER:
  echo ; echo " ERROR: must specify valid experiment user directory "
  exit

#####################################################################
# no valid source directory
#####################################################################

NO_SOURCE:
  echo ; echo " ERROR: no valid source directory present "
  exit

#####################################################################
# cannot create makefile for zetac
#####################################################################

ZETAC_MAKE:
  echo ; echo " ERROR: cannot create makefile for zetac "
  exit

#####################################################################
# cannot create makefile for zetac_colstart
#####################################################################

ZETAC_COLDSTART_MAKE:
  echo ; echo " ERROR: cannot create makefile for zetac_coldstart "
  exit

#####################################################################
# cannot create makefile for write_time
#####################################################################

WRITE_TIME_MAKE:
  echo ; echo " ERROR: cannot create makefile for write_time "
  exit

#####################################################################
# cannot create executable for zetac
#####################################################################

ZETAC_EXEC:
  echo ; echo "  ERROR: cannot create executable for zetac "
  exit

#####################################################################
# cannot create executable for zetac_coldstart
#####################################################################

ZETAC_COLDSTART_EXEC:
  echo ; echo " ERROR: cannot create executable for zetac_coldstart "
  exit

#####################################################################
# cannot create executable for write_time
#####################################################################

WRITE_TIME_EXEC:
  echo ; echo "  ERROR: cannot create executable for write_time "
  exit
