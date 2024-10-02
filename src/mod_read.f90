module mod_read

  use mod_para
  implicit none

  interface get_params
    module procedure get_params_int
    module procedure get_params_real
    module procedure get_params_real_int
    module procedure get_params_real8
    module procedure get_params_string
  end interface

contains

subroutine read_parameters()
  implicit none

  call get_params("problem",problem,"default")
  call get_params("mesh_dir",mesh_dir,"data")
  call get_params("data_dir",data_dir,"data")
  call get_params("simu_time_max",simu_time_max,1.0)
  call get_params("cfl_number",cfl_number,default_cfl_number)
  call get_params("timestep",timestep,-1.0)
  call get_params("export_grdsurf_velo",export_grdsurf_velo,0)
  call get_params("export_grdsurf_displ",export_grdsurf_displ,0)
  call get_params("export_grdsurf_strain",export_grdsurf_strain,0)
  call get_params("export_wave",export_wave,0)
  call get_params("export_wave_component",export_wave_component,1)
  call get_params("export_wave_timestep",export_wave_timestep,10)
  call get_params("export_media",export_media,0)
  call get_params("fault_snap_skip",fault_snap_skip,1)
  call get_params("wave_snap_skip",wave_snap_skip,1)
  call get_params("grdsurf_snap_skip",grdsurf_snap_skip,1)
  call get_params("flux_method",flux_method,2)
  call get_params("plasticity",plasticity,0)
  call get_params("use_damp",use_damp,0)
  call get_params("use_pml",use_pml,0)
  call get_params("thermalpressure",thermalpressure,0)
  call get_params("input_stress_type",input_stress_type,0)
  call get_params("friction_law",friction_law,0)
  call get_params("smooth_load",smooth_load,0)
  call get_params("smooth_load_time",smooth_load_time,0)
  call get_params("initial_condition_wave",initial_condition_wave,0)
  call get_params("src_loc_x",src_loc(1),0)
  call get_params("src_loc_y",src_loc(2),0)
  call get_params("src_loc_z",src_loc(3),0)
  call get_params("src_gaussian_width",src_gaussian_width,1e-16)
  call get_params("src_m0",src_m0,0)
  call get_params("src_mxx",src_mxx,1)
  call get_params("src_myy",src_myy,1)
  call get_params("src_mzz",src_mzz,1)
  call get_params("src_myz",src_myz,0)
  call get_params("src_mxz",src_mxz,0)
  call get_params("src_mxy",src_mxy,0)
  call get_params("ForcedRup",ForcedRup,0)
  call get_params("ForcedRup_x0",ForcedRup_x0,0)
  call get_params("ForcedRup_y0",ForcedRup_y0,0)
  call get_params("ForcedRup_Vr",ForcedRup_Vr,2400)
  call get_params("ForcedRup_rcrit",ForcedRup_rcrit,4e3)
  call get_params("ForcedRup_t0",ForcedRup_t0,0.5)
  call get_params("Tvisc",Tvisc,0.05)
  call get_params("cohesion",cohesion,0)
  call get_params("blkfric",blkfric,0.9)
  call get_params("SmaxAngle",SmaxAngle,45)
  !call get_params("coef_bxx",coef_bxx,0)
  !call get_params("coef_byy",coef_byy,0)
  !call get_params("coef_bxy",coef_bxy,0)
  !call get_params("fluidpres_profile_h1",fluidpres_profile_h1,0)
  !call get_params("fluidpres_profile_h2",fluidpres_profile_h2,0)
  !call get_params("fluidpres_profile_o1",fluidpres_profile_o1,1)
  !call get_params("fluidpres_profile_o2",fluidpres_profile_o2,0)
  call get_params("RS_V0",RS_V0,1e-6)
  call get_params("RS_f0",RS_f0,0.6)
  call get_params("RS_fw",RS_fw,0.2)

end subroutine

subroutine get_params_string(pnm,p,def_p)
  use, intrinsic :: iso_fortran_env, only:  output_unit
  use yaml, only: parse, error_length
  use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
  implicit none
  character(len=*),intent(in) :: pnm
  character(len=*),intent(in) :: def_p
  character(len=*),intent(out) :: p
  class(type_node), pointer :: root
  character(len=error_length) :: error
  type (type_error), pointer :: io_err

  root => parse("parameters.yaml", unit=100, error=error)
  select type (root)
  class is (type_dictionary)
    p = root%get_string(pnm,error=io_err)
    if (associated(io_err)) then
      p = def_p
      if(masternode) print*,trim(io_err%message)
    endif
    if(masternode) print*,pnm,' = ',trim(p)
  end select
  call root%finalize()
  deallocate(root)

end subroutine

subroutine get_params_int(pnm,p,def_p)
  use, intrinsic :: iso_fortran_env, only:  output_unit
  use yaml, only: parse, error_length
  use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
  implicit none
  character(len=*),intent(in) :: pnm
  integer,intent(out) :: p
  integer,intent(in) :: def_p
  class(type_node), pointer :: root
  character(len=error_length) :: error
  type (type_error), pointer :: io_err

  root => parse("parameters.yaml", unit=100, error=error)
  select type (root)
  class is (type_dictionary)
    p = root%get_integer(pnm,error=io_err)
    if (associated(io_err)) then
      p = def_p
      if(masternode) print*,trim(io_err%message)
    endif
    if(masternode) print*,pnm,' = ',p
  end select
  call root%finalize()
  deallocate(root)

end subroutine

subroutine get_params_real8(pnm,p,def_p)
  use, intrinsic :: iso_fortran_env, only:  output_unit
  use yaml, only: parse, error_length
  use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
  implicit none
  character(len=*),intent(in) :: pnm
  real(kind=rkind),intent(out) :: p
  real*8,intent(in) :: def_p
  class(type_node), pointer :: root
  character(len=error_length) :: error
  type (type_error), pointer :: io_err
  !integer :: p,def_p

  root => parse("parameters.yaml", unit=100, error=error)
  select type (root)
  class is (type_dictionary)
    p = root%get_real(pnm,error=io_err)
    if (associated(io_err)) then
      p = def_p
      if(masternode) print*,trim(io_err%message)
    endif
    if(masternode) print*,pnm,' = ',p
  end select
  call root%finalize()
  deallocate(root)

end subroutine

! for compatibility
subroutine get_params_real(pnm,p,def_p)
  implicit none
  character(len=*),intent(in) :: pnm
  real(kind=rkind),intent(out) :: p
  real*4,intent(in) :: def_p
  call get_params_real8(pnm,p,dble(def_p))
end subroutine

! for compatibility
subroutine get_params_real_int(pnm,p,def_p)
  implicit none
  character(len=*),intent(in) :: pnm
  real(kind=rkind),intent(out) :: p
  integer,intent(in) :: def_p
  call get_params_real8(pnm,p,dble(def_p))
end subroutine

end module
