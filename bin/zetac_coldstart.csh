#!/bin/csh -v

#####################################################################
#
# Purpose: script template to run zetac_coldstart
#
# Authors: Steve.Garner@noaa.gov and Chris.Kerr@noaa.gov
#
######################################################################

#FRE scheduler-options
#SBATCH -D /ncrc/home2/Steve.Garner/SIGMA/exp/dp80L10/JOBS/
#SBATCH -o /ncrc/home2/Steve.Garner/SIGMA/exp/dp80L10/JOBS//%x.o%j
#SBATCH -J sigma_coldstart
#SBATCH -A gfdl_w
#SBATCH -t 0:15:00
#SBATCH --cluster=c4
#SBATCH --partition=batch
#SBATCH --nodes=9
#SBATCH --mail-type=fail
#SBATCH --qos=urgent

###########s##########################################################
# retrieve runtime parameter list
#####################################################################

  set basname     = /ncrc/home2/Steve.Garner/SIGMA
  set expname     = dp80L10
  set jobname     = zetac_coldstart

  set identity    = May03h02m22s03  
  set num_proc_x  = 16
  set num_proc_y  = 16
  set num_proc    = 9
  set num_segs    = 1

  set time_limit  = 0:15:00

  @ num_proc  = $num_proc_x * $num_proc_y 
  @ stack_size  = 100000 / $num_proc

  setenv MPICH_UNEX_BUFFER_SIZE 200000000
  setenv MPICH_MAX_SHORT_MSG_SIZE 64000
  
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
  
  if (! -d $WRKDIR)          mkdir -p $WRKDIR
  if (! -d $ARCDIR/restart)  mkdir -p $ARCDIR/restart
  if (! -d $ARCDIR/external) mkdir -p $ARCDIR/external
  if (! -d $ARCDIR/nmlists)  mkdir -p $ARCDIR/nmlists
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
  alias nccombine  `which mppnccombine` -64

#####################################################################
# load modules
#####################################################################

  module unload fre
  module load fre/bronx-19

  setenv MPI_BUFFER_MAX 2000
  setenv MPI_GROUP_MAX  100

#####################################################################
# gernerate time and grid parameters
#####################################################################

  cp $JOBDIR/exp/coupler.nml input.nml

  cat >> input.nml << ENDNL
&zetac_layout_nml
    layout = 1, 1
/
ENDNL

  set executable = zetac_write_time.x
  gunzip -cf $JOBDIR/exec/${executable}.gz > $executable
  chmod +x $executable

  if (-e time_stamp.out) rm time_stamp.out
  if (-e horiz_grid.out) rm horiz_grid.out

  srun -n 1 ./$executable

  setenv error $status
  if ($error != 0) goto DONE

#####################################################################
# name the local output files
#####################################################################

  set log_file      = logfile.000000.out
  set landmask_file = landmask.nc
  set sst_file      = sst.nc

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

echo $cwd
cat input.nml

#####################################################################
# collect input data
#####################################################################

  set topog_file = topography_2min.nc

  ln -s $FMSDIR/$topog_file INPUT/topography_data.nc

#####################################################################
# copy tables to work directory (diag table for base date)
#####################################################################

  cd $WRKDIR

  cp $JOBDIR/exp/diag_table .
  cp $JOBDIR/exp/field_table .

#####################################################################
# run zetac_coldstart
#####################################################################

  echo ; echo " experiment is $expname "
  set executable = ${jobname}.x
  echo " Running $executable ... "

  gunzip -cf $JOBDIR/exec/${executable}.gz > $executable
  chmod +x $executable

  srun -n $num_proc ./$executable > $OUTDIR/sigma_coldstart_stdout

  setenv error $status
  if ($error != 0) goto DONE

#####################################################################
# prefix starting time to archived files
#####################################################################

  set dates = `head -2 time_stamp.out`
  echo $dates[1-7] > time_stamp.out
  set date_name = `time_stamp`

  set restart_write = ${date_name}.res.nc
  set logfile_save  = ${date_name}.log
  set namelist_save = ${date_name}.nml

#####################################################################
# combine restart files in multi-fileset case
#####################################################################

  cd RESTART

  set multfiles = `/bin/ls *.res.nc.0000`

  if ($#multfiles >= 1) then

    foreach file ($multfiles)
      set ncname = $file:r
      rm -f $ncname
      nccombine $ncname
      setenv error $status
      if ($error != 0) then
        echo netcdf combine procedure failed for $ncname
      endif
    end    
    rm -f *.res.nc.*

  endif

#####################################################################
# move restart files to archive
#####################################################################

  set resfiles = `/bin/ls *.res.nc *.res`

  if ($#resfiles == 1) then

      if (-e $ARCDIR/restart/$restart_write) then
        mv $ARCDIR/restart/$restart_write $ARCDIR/restart/old_${restart_write}
        echo ; echo "renamed existing restart file"
      endif
    cp $resfiles $ARCDIR/restart/$restart_write

  else if ($#resfiles > 1) then

    set restart_write = $restart_write".tar"
      if (-e $ARCDIR/restart/$restart_write) then
        mv $ARCDIR/restart/$restart_write $ARCDIR/restart/old_${restart_write}
        echo "renamed existing restart tar file"
      endif
    tar -cf temp.tar $resfiles
    mv temp.tar $ARCDIR/restart/$restart_write

  else

    echo ; echo " no restart files were written "    
    goto CHECK
     
  endif

#####################################################################
# save logfile, namelist file to archive
#####################################################################

  cd $WRKDIR
  
  if (-e $log_file) then
     mv $log_file $ARCDIR/log/$logfile_save
  endif
  
  cp input.nml $ARCDIR/nmlists/$namelist_save
  setenv error $status

#####################################################################
# generate run script for long-term storage
#####################################################################

  echo " ...generating run script "

  set runscript = ${jobname}_store.csh

  set datname = $date_name
  set basname = `echo $basname | sed 's/\(\/\)/\\\//g'`

  echo "\
     s/BBBBBB/$basname/g\
     s/CCCCCC/$expname/g\
     s/HHHHHH/$identity/g\
     s/ZZZZZZ/$datname/g" \
 > sed_file
ls sed_file
  cat $JOBDIR/bin/$runscript | sed -f sed_file > tmp
  cp tmp $JOBDIR/bin/$runscript

  rm sed_file tmp

#####################################################################
# batch run script submission
#####################################################################

  chmod +x $JOBDIR/bin/$runscript

  sbatch $JOBDIR/bin/$runscript
  setenv error $status

#####################################################################
# exit
#####################################################################

DONE:
  if ($error != 0) then
     echo ; echo " script exited normally "
  else
    echo ; echo " script exited abnormally "    
  endif

  exit

