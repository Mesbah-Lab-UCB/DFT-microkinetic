program test_2reac
  !
  ! declare variables and modules
  !
  use ZDPlasKin
  implicit none
  double precision, parameter :: gas_pressure     = 1.0d0, &  ! pressure, bar
                                 surface_site_density = 2.87d17

  double precision :: EN, NE, ER, gas_density, gas_temperature, N2_frac
  double precision :: time = 0.0d0, time_out = 0.0d0, time_end = 9.0d3, dtime, correction, h2init, n2init
  double precision :: charge_positive, charge_total
  integer :: i

  !-----------------------------------------------------------------------------------------
  ! Gas temperature read
  character(*), parameter     :: ele_in = 'Ele.dat'
  character(*), parameter     :: Tgas_in = 'Tgas.dat'
  character(*), parameter     :: other_para = 'other_para.dat'
  open(1,file=Tgas_in)
  read(1,*) ! exclude the header
  read(1,*) gas_temperature
  gas_density = 101325.0d0 * gas_pressure &
          / gas_temperature &
          / 1.38d-17  ! gas density, cm-3
  ! initial volume fractions of H2 and N2 (do not need to add to 1; relative values only required)
  !------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------
  open(81,file=ele_in)
  read(81,*) ! exclude the header
  read(81,*) EN, NE
  ER = (EN / ((gas_density)*1.d6)) * 1.d21
  ! the change of ER should not change that much after simulation (verified, around 1% difference)
  !------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------
  open(79,file=other_para)
  read(79,*) ! exclude the header
  read(79,*) N2_frac
  close(79)!

  ! initial volume fractions of H2 and N2 (do not need to add to 1; relative values only required)
  h2init = 1 - N2_frac
  n2init = N2_frac
  !
  ! print
  !
  write(*,'(/,A)') 'MODIFIED TWO-REACTION TEST CASE'
  write (*,*) 'Current Er:', ER
  write (*,*) 'Current ne:', NE
  !
  ! initialization of ZDPlasKin
  !
  call ZDPlasKin_init()
  !
  ! set the physical conditions of the calculation:
  !     the gas temperature and the reduced electric field
  !
  call ZDPlaskin_set_config(QTPLASKIN_SAVE=.true., ATOL=1.d-1)
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=ER )
  !
  ! set initial densities
  !
  call ZDPlasKin_set_density('N2',gas_density*n2init/(h2init + n2init))
  call ZDPlasKin_set_density('H2',gas_density*h2init/(h2init + n2init))
  call ZDPlasKin_set_density('Surf',surface_site_density)
  call ZDPlasKin_set_density('E',LDENS_CONST=.true.)
  call ZDPlasKin_set_density('E',NE)

  !
  ! save densities
  !
  !write(*,'(/,A,$)') 'TYPE 1 TO CONTINUE or 0 TO START FROM BEGINING: '
  !read(*,*) i
  i = 0
  if (i == 1) then
    open(1,file="densities_stat.bin",form='unformatted')
    read(1) density(:)
    close(1)
  endif
  !
  ! print column headers and initial values
  !
  write(*,'(8(A12))') ( trim(species_name(i)), i = 1, species_max )
  write(*,'(8(1pe12.4))') density(:)
  !
  ! time integration
  !
  !write(*,'(/,A,$)') 'ENTER NEW time_end = '
  !read(*,*) time_end
  do while(time .lt. time_end)
    if (time .lt. 100) then
      dtime = -1 ! ajustable timestep
    else if (time .gt. 1000) then
      dtime = 1.
      !write(*,'(1pe12.4)') time, dtime
    else
      dtime = 0.001
    end if


    call ZDPlasKin_timestep(time,dtime)

    time = time + dtime
    if (time .gt. 4000) then
      write(*,'(1pe12.4)') time
    end if

    if ( time > time_out ) then
      if (time .le. 4000) then
        write(*,'(1pe12.4)') time
      end if
      time_out = 2.0 * time
    endif
    !write(*,'(4(1pe12.4))') time, density(:) ! extremly slow
  enddo
  !
  ! for quasineutrality - charge normalization (by Sergey, 12/08/2015)
  !
  do while(time .lt. time_end)
    write(*,'(1pe12.4)') 'fuck?'
    if (time .lt. 1) then
      dtime = -1 ! YOU HAVE TO FIX IT AGAIN, adjustable dtime will work too slow
    else
      dtime = 0.5
    end if
    call ZDPlasKin_timestep(time, dtime)
    time = time + dtime
    call ZDPlasKin_get_density_total(ALL_ION_POSITIVE=charge_positive, ALL_CHARGE=charge_total)
    where ( species_charge(:) > 0 ) ! only positively-charged species
      density(:) = density(:) * ( 1.0d0 - charge_total / charge_positive ) ! quasi-neutrality
    end where
  enddo
  ! print column headers and final values
  !
  write(*,'(8(A12))') ( trim(species_name(i)), i = 1, species_max )
  write(*,'(8(1pe12.4))') density(:)
  !
  ! save densities
  !
  open(1,file="densities_stat.bin",form='unformatted')
  write(1) density(:)
  close(1)


  !
  ! end
  !
  !write(*,'(/,A,$)') 'PRESS ENTER TO EXIT ...'
  !read(*,*)
end program test_2reac