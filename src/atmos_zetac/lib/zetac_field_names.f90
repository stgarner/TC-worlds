module zetac_field_names_mod

implicit none

integer, public, parameter :: num_refs_fields=4
integer, public, parameter :: num_stat_fields=2
integer, public, parameter :: num_diag_fields=7
integer, public, parameter :: num_prog_fields=7
integer, public, parameter :: num_grid_fields=2
integer, public, parameter :: num_tend_fields=4
integer, public, parameter :: num_mean_fields=4

character*8, public, dimension(num_refs_fields) :: refs_names
character*8, public, dimension(num_stat_fields) :: stat_names
character*8, public, dimension(num_diag_fields) :: diag_names
character*8, public, dimension(num_prog_fields) :: prog_names
character*8, public, dimension(num_grid_fields) :: grid_names
character*8, public, dimension(num_tend_fields) :: tend_names
character*8, public, dimension(num_mean_fields) :: mean_names

integer, public, parameter :: ppref_=1, uuref_=2, taref_=3, qvref_=4

integer, public, parameter :: topog_=1, sst_=2

integer, public, parameter :: pp_=1, dp_=2, gz_=3, rh_=4, sd_=5, sm_=6, om_=7

integer, public, parameter :: ps_=1, uu_=2, vv_=3, oo_=4, ta_=5, qv_=6, qc_=7

integer, public, parameter :: zetam_=1, zetaw_=2

integer, public, parameter :: utend_=1, vtend_=2, ttend_=3, qtend_=4

integer, public, parameter :: umean_=1, vmean_=2, tmean_=3, qmean_=4

data refs_names /'ppref', 'uuref', 'taref', 'qvref'/

data stat_names /'topog', 'sst'/

data diag_names /'pp', 'dp', 'gz', 'rh', 'sd', 'sm', 'om'/

data prog_names /'ps', 'uu', 'vv', 'oo', 'ta', 'qv', 'qc'/

data grid_names /'zetam', 'zetaw'/

data tend_names /'utend', 'vtend', 'ttend', 'qtend'/

data mean_names /'umean', 'vmean', 'tmean', 'qmean'/

end module zetac_field_names_mod

