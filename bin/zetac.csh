#!/bin/csh -v

#####################################################################
#
# Purpose: script template to run zetac
#
# Authors: Steve.Garner@noaa.gov and Chris.Kerr@noaa.gov
#
# Version: 1.0
#
######################################################################

#scheduler-options
#SBATCH -D /ncrc/home2/Steve.Garner/SIGMA/exp/dp80L3/JOBS/
#SBATCH -o /ncrc/home2/Steve.Garner/SIGMA/exp/dp80L3/JOBS/%x.o%j
#SBATCH -J sigma
#SBATCH -A gfdl_w
#SBATCH -t 0:10:00
#SBATCH --cluster=c4
#SBATCH --partition=batch
#SBATCH --nodes=9
#SBATCH --mail-type=fail
#SBATCH --qos=urgent

  unset echo

#####################################################################
# retrieve runtime parameter list
#####################################################################

  set basname     = /ncrc/home2/Steve.Garner/SIGMA
  set expname     = dp80L3
  set jobname     = zetac

  set identity    = May04h15m17s52 
  set num_proc_x  = 16
  set num_proc_y  = 16
  set num_segs    = 1

  set time_limit  = 0:10:00

  set segment = 1

  @ num_proc  = $num_proc_x * $num_proc_y
  @ stack_size    = 100000 / $num_proc

#####################################################################
# set the directories
#####################################################################

  setenv RUNDIR $basname
  setenv ARCDIR /lustre/f2/dev/$USER
  setenv TMPDIR /lustre/f2/scratch/$USER
  setenv WRKDIR $TMPDIR/$identity
  setenv FMSDIR $ARCDIR/fms_data
  setenv ARCDIR $ARCDIR/zetac/$expname
  setenv EXPDIR $RUNDIR/exp/$expname
  setenv JOBDIR $ARCDIR/jobs/$identity
  setenv OUTDIR $EXPDIR/JOBS/$identity
  setenv HISDIR $ARCDIR/history
  setenv RESDIR $ARCDIR/restart

  if (! -d $WRKDIR)          mkdir -p $WRKDIR
  if (! -d $RESDIR)          mkdir -p $RESDIR
  if (! -d $HISDIR)          mkdir -p $HISDIR
  if (! -d $ARCDIR/nmlists)  mkdir -p $ARCDIR/nmlists
  if (! -d $ARCDIR/tables)   mkdir -p $ARCDIR/tables
  if (! -d $ARCDIR/log)      mkdir -p $ARCDIR/log
  if (! -d $OUTDIR)          mkdir -p $OUTDIR

  cd $WRKDIR

  if (! -d INPUT  ) mkdir INPUT 
  lfs setstripe INPUT   1048576 -1 12

  if (! -d OUTPUT ) mkdir OUTPUT
  lfs setstripe OUTPUT  1048576 -1 12

  if (! -d RESTART) mkdir RESTART
  lfs setstripe RESTART 1048576 -1 12

#####################################################################
# set alias
#####################################################################

  alias time_stamp ~$USER/ZETAC/scripts/time_stamp.csh -s -t:

#####################################################################
# load modules
#####################################################################

  source $MODULESHOME/init/tcsh
  module load fre/bronx-19

  setenv MPI_BUFFER_MAX 2000
  setenv MPICH_UNEX_BUFFER_SIZE 200000000
  setenv MPICH_MAX_SHORT_MSG_SIZE 64000
  setenv MPI_GROUP_MAX  100

#####################################################################
# generate time and grid parameters
#####################################################################

  cp $JOBDIR/exp/coupler.nml input.nml

  cat >> input.nml << ENDNL
&segments_nml
  num_segs = $num_segs
/
&zetac_layout_nml
    layout = 1, 1
/
ENDNL

  set executable = zetac_write_time.x
  gunzip -cf $JOBDIR/exec/${executable}.gz > $executable
  chmod +x $executable

  if (-e time_stamp.out) rm time_stamp.out
  if (-e horiz_grid.out) rm horiz_grid.out

  set outfile = $OUTDIR/zetac_time

  srun -n 1 ./$executable #> $outfile

  setenv error $status
  if ($error != 0) goto DONE

#####################################################################
# name the local input/output files
#####################################################################

  set log_file      = logfile.out
  set restart_file  = zetac.res.nc

#####################################################################
# assemble namelists
#####################################################################

  cp $JOBDIR/exp/coupler.nml input.nml

  set nml_files = `/bin/ls $JOBDIR/exp/*.nml`
  foreach file ($nml_files)
    if ( $file:t != coupler.nml ) then
      cat $file >> input.nml
    endif
  end

  cat >> input.nml << ENDNL  
&zetac_layout_nml
  layout = $num_proc_x, $num_proc_y
/
&fms_nml
  stack_size = $stack_size
  clock_grain = 'COMPONENT'
  print_memory_usage = .false.
/
&xgrid_nml
  make_exchange_reproduce = .false.
/
ENDNL

#####################################################################
# time stamps
#####################################################################

  mv time_stamp.out temp
  set dates = `head -2 temp`
  echo $dates[1-7] > time_stamp.out
  set date_name_beg = `time_stamp`
  echo $dates[8-14] > time_stamp.out
  set date_name_end = `time_stamp`
  mv temp time_stamp.out

#####################################################################
# collect input data
#####################################################################

  set topog_file = topography_2min.nc

  ln -s $FMSDIR/$topog_file INPUT/topography_data.nc

  set files = `ls -C1 $FMSDIR/data`

#  foreach file ($files)
#    cp $FMSDIR/data/$file INPUT
#  end

  echo
  echo "RUNNING EXPERIMENT: " $expname
  echo "CPU TIME LIMIT IS:  " $time_limit
  echo "PROCESSOR COUNT IS: " $num_proc
  echo "SEGMENT COUNT IS:   " $num_segs
  echo "START TIME IS:      " $date_name_beg
  echo "FINAL TIME IS:      " $date_name_end
  echo

#####################################################################
# name the input restart file
#####################################################################

  set restart_read = $date_name_beg".res.nc"

#####################################################################
# copy restart file to input directory
#####################################################################

  @ count = 2 + 2 * $segment
  head -$count time_stamp.out > time_stamps
  tail -2 time_stamps > temp
  head -1 temp > time_stamp.out
cat time_stamps
echo time_stamp.out
  set date_name = `time_stamp`
  set restart_read = ${date_name}.res.nc
  echo " input restart file is $restart_read "

  cd INPUT

  rm *.res *.res.nc
  if (-e $RESDIR/${restart_read}.tar ) then
    cp $RESDIR/${restart_read}.tar temp.tar
    tar -xf temp.tar ; rm -f temp.tar
    if (-e $restart_read) mv $restart_read $restart_file
  else if (-e $RESDIR/$restart_read) then
    cp $RESDIR/$restart_read $restart_file
  else
     goto NO_RESTART
  endif

  set resfiles = `/bin/ls zetac*.res.nc`
  if ($#resfiles == 0) goto NO_RESTART

#####################################################################
# copy tables to work directory
#####################################################################

  cd $WRKDIR
  
  cp $JOBDIR/exp/*table .

#####################################################################
# modify namelist to reflect time of current segment
#####################################################################

  if ($num_segs > "1") then
    $RUNDIR/scripts/modify_input.csh time_stamps input.nml
  endif

#  rm time_stamps

#####################################################################
# run ZETAC
#####################################################################

  echo ; echo " experiment is "$expname
  set executable = ${jobname}.x
  echo " Running "$executable
  date ; echo

  gunzip -cf $JOBDIR/exec/${executable}.gz > $executable
  chmod +x $executable

  limit stacksize unlimited

  set outfile = $OUTDIR/sigma_stdout
  if ( $num_segs > 1 ) set outfile = ${outfile}_${segment}

  srun -n $num_proc ./$executable > $outfile
  
  setenv error $status
  if ($error != 0) goto DONE

#####################################################################
# prefix ending time to archived files
#####################################################################

  tail -1 time_stamps > time_stamp.out
  set date_name = `time_stamp`

  set history_tar   = ${date_name}.nc.tar
  set restart_write = ${date_name}.res.nc
  set logfile_save  = ${date_name}.log
  set namelist_save = ${date_name}.nml
  set tables_save   = ${date_name}.tbl.tar

#####################################################################
# move restart files to archive
#####################################################################
set echo
  cd RESTART

  set ncfiles = `/bin/ls | grep nc.0000`
  alias nccombine  `which mppnccombine` -64

  unset echo
  foreach file ($ncfiles)
    set ncname = $file:r
    if (-e $ncname) rm $ncname
#    nccombine $ncname
    set infiles = ( `ls $ncname.????` )
#    combinenc $infiles $ncname && rm -f $infiles
    nccombine $ncname && rm -f $infiles
  end
set echo

  find -iname '*.res' > restart.list
  find -iname '*.res.nc' >> restart.list
  set resfiles = `wc -l restart.list | awk '{print $1}'`

  if ( $resfiles > 0 ) then

    set restart_write = ${restart_write}.tar
    tar -b 1000 -cf $RESDIR/$restart_write --files-from restart.list
    echo " restart file is $RESDIR/$restart_write "

  else

    set segment = $num_segs
    echo ; echo " no restart files were written "    
    goto NO_RESTART

  endif

#####################################################################
# move history files to archive
#####################################################################

  cd $WRKDIR

  set hisfiles = `/bin/ls *nc`

  if ($#hisfiles != 0) then
    tar -cf temp.tar $hisfiles
    mv temp.tar $HISDIR/$history_tar
    echo " history file is "$HISDIR/$history_tar
    goto RESUBMIT
  endif

  set hisfiles = `/bin/ls *nc.0000`
  if ($#hisfiles == 0) goto NO_HISTORY

#  unset echo
  foreach file ($hisfiles)
    set group = $file:r:r
    if ($group == "default") then
      set fname = ${date_name}.nc
    else
      set fname = ${date_name}.$group.nc
    endif
    set groupfiles = `/bin/ls ${group}.nc.????`
    foreach gfile ($groupfiles)
      mv $gfile ${fname}.$gfile:e
    end
  end

  set hisfiles = `/bin/ls *nc.????`
  set tarfile = history.tar

  set SEGDIR = $HISDIR/$identity
  if (! -d $SEGDIR) mkdir $SEGDIR
  if ( $num_segs > 1 ) then
    set SEGDIR = $SEGDIR/segment_$segment
    mkdir $SEGDIR
  endif

goto NOTAR

  tar -cf $tarfile $hisfiles
  setenv error $status
  if ($error != 0) goto DONE
  mv $tarfile $SEGDIR
  rm -f $hisfiles
  goto STORE

NOTAR:

  mv $hisfiles $SEGDIR

#####################################################################
# generate run script for long-term storage
#####################################################################

STORE:

  echo " ...generating run script for long-term storage"

  set runscript = ${jobname}_store.csh
  set template = $RUNDIR/bin/zetac_store.template

  set datname = $date_name
  set basname = `echo $basname | sed 's/\(\/\)/\\\//g'`

  echo "\
     s/BBBBBB/$basname/g\
     s/CCCCCC/$expname/g\
     s/HHHHHH/$identity/g\
     s/SSSSSS/$segment/g\
     s/ZZZZZZ/$datname/g" \
 > sed_file

  cat $template | sed -f sed_file > tmp
  cp tmp $JOBDIR/bin/$runscript
  rm sed_file tmp

###################################################################
# combine netcdf history files
#####################################################################

  cd $SEGDIR

  set runscript = combine.csh
  set template  = $RUNDIR/bin/combine.template

  echo "\
     s/BBBBBB/$basname/g\
     s/CCCCCC/$expname/g\
     s/HHHHHH/$identity/g\
     s/SSSSSS/$segment/g\
     s/NNNNNN/$num_segs/g" \
  > sed_file

  cat $template | sed -f sed_file > $runscript
  rm sed_file
  chmod +x $runscript

  echo " submitting nccombine for segment "$segment
  echo " history file will be "$ARCDIR/history/$history_tar:r

  sbatch $runscript > tmp

  set jobid = `tail -1 tmp`
  rm tmp
   echo "combine jobid is "$jobid
  checkjob -v $jobid

#####################################################################
# edit and resubmit runscript
#####################################################################

RESUBMIT:
  set runscript = $JOBDIR/bin/zetac.csh

  if ($segment < $num_segs) then
  
    @ next = $segment + 1

    sed -e "s/segment = $segment/segment = $next/g" $runscript > tmp

    cp tmp $runscript
    rm tmp

    echo " submitting segment $next of $num_segs "

    chmod +x $runscript
    sbatch $runscript > tmp

    set jobid = `tail -1 tmp`
    rm tmp
     echo "production jobid is "$jobid
    checkjob -v $jobid

  endif

#####################################################################
# save tables, logfile, namelist file to archive
#####################################################################

  if ( $segment == 1 ) then

    cd $JOBDIR/exp 

    set tbl_files = `/bin/ls *table`
    tar -cf temp.tar $tbl_files
    mv temp.tar $ARCDIR/tables/$tables_save

    cd $WRKDIR
    
    if (-e $log_file) then
      mv $log_file $ARCDIR/log/$logfile_save
    endif

    cp input.nml $ARCDIR/nmlists/$namelist_save

  endif

#####################################################################
# exit
#####################################################################

DONE:
  if ($error != 0) then
    echo ; echo " ERROR: abnormal job termination "
    echo ; echo " work directory is `pwd` "
  else
    if ($segment > $num_segs) then
      echo ; echo " cleaning up $JOBDIR "
      rm -rf $JOBDIR
    endif
  endif

  echo ; date

  exit

#####################################################################
# no grid file
#####################################################################

NO_GRID:
  echo ; echo " ERROR: grid generation failed "
  exit

#####################################################################
# NO OVERRIDE FILE
#####################################################################

NO_OVERRIDE:
  echo ; echo " ERROR: requested override file is out of range or not found "
  exit

#####################################################################
# NO HISTORY WRITTEN
#####################################################################

NO_HISTORY:
  echo ; echo " WARNING: no history files were written "
  exit

#####################################################################
# NO RESTART FILE
#####################################################################

NO_RESTART:
  echo ; echo " ERROR: no restart file could be found "
  setenv error 1
  exit

exit

