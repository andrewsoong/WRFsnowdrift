module snowsub

  !-----------------------------------------------------------------------
  ! Module to calculate the sublimation of snowdrift particles
  !-----------------------------------------------------------------------

  contains

  !Function to calculate the sublimate if both radiation and turbulence are on
  function calc_subl_both(pi, rmean, alpha_s, swdown, epsilon_s, &
    glw, kstb, t_phy, ls, rv, d, &
   psis, tracer, mmean, re, xkmv, xkhv, kt, esat, sigma, dt, ims, ime, jms, jme) result(snowsubl_limit_tracer)

     implicit none
     real, intent(in) :: pi,alpha_s, epsilon_s,kstb,ls, rv, d, tracer, mmean, re,dt, &
     esat ,kt,t_phy,rmean,xkmv,xkhv
     real,intent(inout) ::psis,sigma
     real,intent(in)::swdown ,glw
     integer, intent(in) :: ims, ime, jms, jme
     real :: sc, fv, qr_sw, qr_lwabs, qr_lwemiss, qr, n, z, dm, snowsubl_tmp, snowsubl_limit_tracer
     integer :: i, j, k

     !passing or a local variable ...check it  ...14.05.2025

     !shortwave absorption
     qr_sw = pi * rmean**2 * (1 - alpha_s) * swdown

     !longwave absorption
     qr_lwabs = pi * rmean**2 * epsilon_s * glw

     ! Longwave emission
     qr_lwemiss = 4.0 * pi * rmean**2 * epsilon_s * kstb * t_phy**4

     ! Net radiation, limit to positive values -> only apply, if radiation budget is positive
     qr = max(0.0, qr_sw + qr_lwabs- qr_lwemiss)

     sc = xkmv / xkhv !Sc=Schmidt number
     fv = 0.78 + 0.308 * sc**(1 / 3) * sqrt(re) !fv=Transfer coefficient 
     n = ((ls / (kt * t_phy)) * (ls / (rv * t_phy) - 1) + (rv * t_phy) / (d * esat)) / fv       !modication using eqn 21,schmid 2013 and eqn 11/12 sharma (2023)
     z = ((2.0 * pi * rmean * sigma) - ((qr / (kt * t_phy)) * ((ls / (rv * t_phy)) - 1.0)))
     dm = (z / n) * dt    ! mass loss of single particle (kg)

     !sublimation loss rate, Schmid 2021, equation (28)
     psis= dm / mmean

     !loss of sublimation per time step, Schmid 2021, equation (30)
     snowsubl_tmp = psis * tracer

     ! limit snowsubl with tracer and qvs as upper and lower boundary to snowsubl_limit, use that to update snowsubl, moist, tracer,
     ! temperature.., limit snowsubl to <0
     snowsubl_tmp = min(0.0, snowsubl_tmp) ! limit snowsubl to negative values, msaigger

     snowsubl_limit_tracer = max(-tracer, snowsubl_tmp) ! limit snowsubl to available snow particles, msaigger

  end function calc_subl_both
  function calc_subl_radiation(pi, rmean, alpha_s, swdown, epsilon_s, glw, kstb, t_phy, ls, rv, d, &
   kt, esat, sigma, nu, dt, tracer, mmean,psis) result(snowsubl_limit_tracer)
     implicit none
     ! Input Parameters
     real, intent(in) :: pi, rmean, alpha_s, swdown, epsilon_s, glw, kstb, t_phy
     real, intent(in) :: ls, rv, d, kt, esat, nu, dt, tracer, mmean
     real,intent(inout) ::sigma,psis
    
     ! Local Variables
     real :: qr_sw, qr_lwabs, qr_lwemiss, qr
     real :: n, z, snowsubl_tmp,dm
     real :: snowsubl_limit_tracer

     ! --- Shortwave absorption ---
     qr_sw = pi * rmean**2 * (1.0 - alpha_s) * swdown

     ! --- Longwave absorption ---
     qr_lwabs = pi * rmean**2 * epsilon_s * glw

     ! --- Longwave emission ---
     qr_lwemiss = 4.0 * pi * rmean**2 * epsilon_s * kstb * t_phy**4

     ! --- Net radiation, limit to positive values -> only apply, if radiation budget is positive ---
     qr = max(0.0, qr_sw + qr_lwabs - qr_lwemiss)

     ! --- Z term (including radiation effects) ---
     z = 2.0 * pi * rmean * sigma - (qr / (kt * nu * t_phy)) * (ls / (rv * t_phy) - 1.0)

     ! --- N term ---
     n = (ls / (kt * nu * t_phy)) * (ls / (rv * t_phy) - 1.0) + (rv * t_phy) / (nu * d * esat)

     ! --- Sublimation mass loss ---
     dm = (z / n) * dt

     !sublimation loss rate, Schmid 2021, equation (28)
     psis = dm / mmean

     !loss of sublimation per time step, Schmid 2021, equation (30)
     snowsubl_tmp = psis* tracer

     ! limit snowsubl with tracer and qvs as upper and lower boundary to snowsubl_limit, use that to update snowsubl, moist, tracer,
     ! temperature.., limit snowsubl to <0
     snowsubl_tmp = min(0.0, snowsubl_tmp) ! limit snowsubl to negative values, msaigger

     snowsubl_limit_tracer = max(-tracer, snowsubl_tmp) ! limit snowsubl to available snow particles, msaigger

  end function calc_subl_radiation

  function calc_subl_turbulence(pi, rmean, t_phy, ls, rv, d, kt, esat, sigma, re, &
   xkmv, xkhv, dt, tracer, mmean,psis) result(snowsubl_limit_tracer)
     implicit none
     ! Input Parameters
     real, intent(in) :: pi, rmean, t_phy, ls, rv, d, kt, esat, re, dt, tracer, mmean,xkmv,xkhv
     real, intent(inout) ::psis,sigma
   

     ! Local Variables
     real :: sc, fv, n, z, dm, snowsubl_tmp
     real :: snowsubl_limit_tracer

     ! --- Schmidt Number ---
     sc = xkmv/ xkhv ! Sc = Schmidt number

     ! --- Transfer Coefficient (fv) ---
     fv = 0.78 + 0.308 * sc**(1.0 / 3.0) * sqrt(re)

     ! --- N term ---
     n = ((ls / (kt * t_phy)) * ((ls / (rv * t_phy)) - 1.0) + (rv * t_phy) / (d * esat)) / fv

     ! --- Z term (Only Turbulence) ---
     z = 2.0 * pi * rmean * sigma

     ! --- Sublimation mass loss of single particle ---
     dm = (z / n) * dt

     !sublimation loss rate, Schmid 2021, equation (28)
     psis = dm / mmean

     !loss of sublimation per time step, Schmid 2021, equation (30)
     snowsubl_tmp = psis * tracer

     ! limit snowsubl with tracer and qvs as upper and lower boundary to snowsubl_limit, use that to update snowsubl, moist, tracer,
     ! temperature.., limit snowsubl to <0
     snowsubl_tmp = min(0.0, snowsubl_tmp) ! limit snowsubl to negative values, msaigger

     snowsubl_limit_tracer = max(-tracer, snowsubl_tmp) ! limit snowsubl to available snow particles, msaigger

  end function calc_subl_turbulence
  function calc_subl_no_rad_turb(pi, rmean, t_phy, ls, rv, d, kt, esat, sigma, re, dt, tracer, mmean,psis) result(snowsubl_limit_tracer)
     implicit none
     ! Input Parameters
     real, intent(in) :: pi, rmean, t_phy, ls, rv, d, kt, esat,re, dt, tracer, mmean
     real, intent(inout) ::sigma,psis    

     ! Local Variables
     real :: nu, n, z, dm, snowsubl_tmp
     real :: snowsubl_limit_tracer

     ! --- Nusselt Number Calculation ---
     if (re <= 10.0) then
        nu = 1.79 + 0.606 * sqrt(re)
     else
        nu = 1.88 + 0.580 * sqrt(re)
     endif

     ! --- Z term (No Radiation, No Turbulence) ---
     z = 2.0 * pi * rmean * sigma

     ! --- N term ---
     n = (ls / (kt * nu * t_phy)) * ((ls / (rv * t_phy)) - 1.0) + (rv * t_phy) / (nu * d * esat)

     ! --- Sublimation mass loss ---
     dm = (z / n) * dt

     !sublimation loss rate, Schmid 2021, equation (28)
     psis = dm / mmean

     !loss of sublimation per time step, Schmid 2021, equation (30)
     snowsubl_tmp = psis * tracer

     ! limit snowsubl with tracer and qvs as upper and lower boundary to snowsubl_limit, use that to update snowsubl, moist, tracer,
     ! temperature.., limit snowsubl to <0
     snowsubl_tmp = min(0.0, snowsubl_tmp) ! limit snowsubl to negative values, msaigger

     snowsubl_limit_tracer = max(-tracer, snowsubl_tmp) ! limit snowsubl to available snow particles, msaigger
     snowsubl_limit_tracer = min(0.0, snowsubl_limit_tracer)
     
  end function calc_subl_no_rad_turb





  subroutine snowsublimation(dt, config_flags, tracer, rho_phy,      &
  dz8w, t_phy, p_phy, rmean,                                         &
  alt, snowsubl, tempch, moistch, tempchacc, moistchacc,             &
  moist, swdown, psis, sigma, rhum,                &
  ustar, zmid, zlow, setvel, snowsublacc,glw, xkmv, xkhv,                              &
  ids, ide, jds, jde, kds, kde,                                      &
  ims, ime, jms, jme, kms, kme,                                      &
  its, ite, jts, jte, kts, kte)
     use module_configure         !/frame
     use module_state_description !/frame

     implicit none
     type(grid_config_rec_type), intent(in) :: config_flags
     integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
     ims, ime, jms, jme, kms, kme, &
     its, ite, jts, jte, kts, kte

     real, dimension(ims:ime, kms:kme, jms:jme, num_tracer), intent(inout) :: tracer
     real, dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(inout) :: moist
     real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: t_phy, snowsublacc,&
     tempchacc, moistchacc
     real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: tempch, moistch, psis, &
     sigma, rhum, snowsubl
     real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: rho_phy, alt, p_phy, zmid, &
     zlow, dz8w, setvel,rmean, xkmv, xkhv
     real, dimension(ims:ime, jms:jme), intent(in) :: ustar, swdown, glw
     real, intent(in) :: dt
     integer :: i, k, j

     !Local variables
     real :: rhoi, nue, kt, ls, rv, d, cp, pi,alpha_s, epsilon_s, kstb
     real :: mmean, re, nu, dm, qvs
     real :: esat, ew, rho_dry, snowsubl_limit, snowsubl_limit_tracer, snowsubl_tmp, &
     moist_tmp, t_tmp, delta_qv_limit_0, snowsubl_limit_moist_tmp, &
     delta_t_limit_tmp, snowsubl_limit_no_rad, &
     t_limit_tmp, esat_tmp, qvs_tmp, delta_qv_limit



     !---------------CONSTANTS------------------------------------------------
     rhoi = 900         !ice density [kg m-3]
     nue = 1.2e-5       !kinematic viscosity of air [m2 s-2]
     kt = 0.024         !thermal conductivity of air [J s-1 m-1 K-1]
     ls = 2.838e6       !latent heat of sublimation [J kg-1]
     rv = 461.5         !gas constant for water vapor [J kg-1 K-1]
     cp = 1005          !specific heat capacity of dry air [J kg-1 K-1]
     d = 2.25e-5        !molecular diffusivity of water vapor in air [m2 s-1]
     pi = 3.14159265359
     alpha_s = 0.9      !snow albedo
     epsilon_s = 0.99   ! snow emissivity
     kstb = 5.67e-8     ! Stefan-Boltzmann constant [W m^-2 K^-4]
     !-----------------------------------------------------------------------

     do j = jts, jte
        do i = its, ite
           do k = kts, kte !Height loop

              ! ----------------------------------------
              ! PREPARE AIR MOISTURE
              ! ----------------------------------------
              !if the snow concentration is very low, no sublimation of particles
              if (tracer(i, k, j, p_snowconc) .le. 1e-8) then
                 snowsubl(i, k, j) = 0.0
                 go to 105 !termination of snowdrift sublimation
              else
                 !calculate dry air mass density
                 rho_dry = 1.0 / alt(i, k, j)

                 !saturation water vapor pressure (Clausius-Clapeyron)
                 esat = 611.2 * exp((ls / rv) * ((1 / 273.15) - (1 / t_phy(i, k, j))))

                 !saturation mixing ratio
                 qvs = (0.622 * esat) / (p_phy(i, k, j) - esat)

                 !limit water vapor mixing ratio to posivive values, relative humidity to 100%
                 moist(i, k, j, p_qv) = max(0.0, moist(i, k, j, p_qv))
                 rhum(i, k, j) = min(100., (moist(i, k, j, p_qv) / qvs) * 100)

                 !calculate water vapor pressure
                 ew = (rhum(i, k, j) * esat) / 100

                 !water vapor deficit with respect to ice
                 sigma(i, k, j) = (ew / esat) - 1

                 if (sigma(i, k, j) .gt. -0.001) then !if relative humidity > 99.9 % no sublimation
                    snowsubl(i, k, j) = 0.0
                    go to 105 !termination of snowdrift sublimation
                 else

                    ! ---------------------------------------
                    ! SNOW PARTICLES
                    ! ---------------------------------------

                    !calculation of mean particle mass, Schmid 2021, equation (29)
                    mmean = (4 / 3) * pi * rhoi * rmean(i, k, j)**3

                    !reynolds particle number, Schmid 2021, equation (26)
                    re = 2 * rmean(i, k, j) * setvel(i, k, j) / nue


                    if (config_flags%subl_opt == 2) then

                       ! conservative approach: calcualte first moisture change and temperature change out of sublimation
                       ! -> restrict sublimation to saturation mixing ratio of new temperature
                       ! -> use that restricted sublimation for moisture change and temperature change

                       ! save current water vapor mixing ratio and temperature
                       moist_tmp = moist(i, k, j, p_qv)
                       t_tmp = t_phy(i, k, j)

                       ! limit sublimation in first step to current saturation mixing ratio
                       delta_qv_limit_0 = max(0.0, qvs - moist_tmp)



                       if (config_flags%subl_rad_opt == 1 .AND. config_flags%subl_turb_opt == 1) then
                          snowsubl_limit_tracer = calc_subl_both(pi, rmean(i,k,j), alpha_s, swdown(i,j), epsilon_s, &
                          glw(i,j), kstb, t_phy(i, k, j), ls, rv, d,psis(i, k, j), &
                          tracer(i, k, j,p_snowconc), mmean, re, xkmv(i,k,j), xkhv(i,k,j), kt, esat, sigma(i,k,j), dt, &
                          ims, ime, jms, jme)  !If both radiation and turbulence are on , the cooling effect is calculated using the volume sublimated in case of neither radiation nor turbulence

                          snowsubl_limit_no_rad = calc_subl_turbulence(pi, rmean(i,k,j), t_phy(i, k, j), &  
                          ls, rv, d, kt, esat, sigma(i,k,j), re, xkmv(i,k,j), xkhv(i,k,j), dt,tracer(i,k,j,p_snowconc),mmean,psis(i,k,j))  !If both radiation and turbulence are on , the cooling effect is 
                          
                          snowsubl_limit_moist_tmp = max(snowsubl_limit_no_rad, &
                          -delta_qv_limit_0 * rho_dry)
                          !only calculated due to the turbulence only(snowsubl_limit_no_rad  - this variable is taken as the volume of sublimated snow particles without the effect of the radiation)
                          ! calculate maximum temperature change
                          delta_t_limit_tmp = (ls * (-snowsubl_limit_moist_tmp)) / (cp * rho_dry)
                       else if (config_flags%subl_rad_opt == 1) then

                          snowsubl_limit_tracer = calc_subl_radiation(pi, rmean(i,k,j),&
                           alpha_s, swdown(i,j), epsilon_s, glw(i,j), kstb, t_phy(i, k, j), &
                           ls, rv, d, kt, esat, sigma(i,k,j),nu, dt,tracer(i,k,j,p_snowconc), mmean,psis(i,k,j))  !If only radiation is on , the cooling effect is calculated using the volume sublimated in case of neither radiation nor turbulence

                          snowsubl_limit_no_rad = calc_subl_no_rad_turb(pi, rmean(i,k,j), &
                          t_phy(i, k, j), ls, rv, d, kt, esat, sigma(i,k,j), re, dt,tracer(i,k,j,p_snowconc),mmean,psis(i,k,j))  !If only radiation is on , the cooling effect is calculated using the volume sublimated in case of neither radiation nor turbulence 
                          snowsubl_limit_moist_tmp = max(snowsubl_limit_no_rad, &
                          -delta_qv_limit_0 * rho_dry)
                          ! calculate maximum temperature change
                          delta_t_limit_tmp = (ls * (-snowsubl_limit_moist_tmp)) / (cp * rho_dry)

                       else if (config_flags%subl_turb_opt == 1) then
                          snowsubl_limit_tracer = calc_subl_turbulence(pi, rmean(i,k,j), t_phy(i, k, j), &  
                          ls, rv, d, kt, esat, sigma(i,k,j), re, xkmv(i,k,j), xkhv(i,k,j), dt,tracer(i,k,j,p_snowconc),mmean,psis(i,k,j)) 
                          
                          snowsubl_limit_moist_tmp = max(snowsubl_limit_tracer, &
                          -delta_qv_limit_0 * rho_dry)
                          ! calculate maximum temperature change
                          delta_t_limit_tmp = (ls * (-snowsubl_limit_moist_tmp)) / (cp * rho_dry)

                       else
                          snowsubl_limit_tracer = calc_subl_no_rad_turb(pi, rmean(i,k,j), &
                          t_phy(i, k, j), ls, rv, d, kt, esat, sigma(i,k,j), re, dt,tracer(i,k,j,p_snowconc),mmean,psis(i,k,j))
                          
                          snowsubl_limit_moist_tmp = max(snowsubl_limit_tracer, &
                          -delta_qv_limit_0 * rho_dry)
                          ! calculate maximum temperature change
                          delta_t_limit_tmp = (ls * (-snowsubl_limit_moist_tmp)) / (cp * rho_dry)
                       endif


                       ! calculate saturation mixing ratio for sublimation-cooled air
                       t_limit_tmp = t_tmp - delta_t_limit_tmp
                       esat_tmp = 611.2 * exp((ls / rv) * ((1 / 273.15) - (1 / t_limit_tmp)))  !saturated vapour pressure at the updated temperature
                       qvs_tmp = (0.622 * esat_tmp) / (p_phy(i, k, j) - esat_tmp)  !New saturation mixing ratio after temperature drop

                       ! limit sublimation to updated saturation mixing ratio
                       delta_qv_limit = qvs_tmp - moist_tmp
                       snowsubl_limit = max(snowsubl_limit_tracer, -delta_qv_limit * rho_dry)
                       snowsubl_limit = min(0.0, snowsubl_limit)

                       !! use limited sublimation to update mixing ratio and temperature
                       ! update moisture
                       moist(i, k, j, p_qv) = min(qvs_tmp, moist_tmp - snowsubl_limit / rho_dry)
                       !qvs_tmp is the new saturation mixing ratio after sublimation- that is qvs_tmp is the maximum possible vapour that can air hold in the new updated temperature.
                       !moist_tmp is the current water vapor mixing ratio before sublimation is applied.
                       !snowsubl_limit is the total mass of sublimated snow particles in kg m-3
                       !rho_dry is the dry air density in kg m-3
                       !snowsubl_limit/rho_dry gives us the sublimated snow converted to the mixing ratio
                       !snowsubl_limit is already negative ... therefore,
                       !moist_tmp - snowsubl_limit / rho_dry , we are adding the sublimated water vapour to the current moisture content


                       !write out rate of change of moisture due to snowdrift
                       moistch(i, k, j) = (moist(i, k, j, p_qv) - moist_tmp) / dt ! moisture change rate (kg kg-1 s-1)

                       ! accumulated moisture change
                       moistchacc(i, k, j) = moistchacc(i, k, j) + moistch(i, k, j) * dt ! acc moisture change (kg kg-1)
                       !moistch is the instantaneous rate of change , that is how much vapor mixing ratio is changing at that specific time step.
                       !moistchacc is the total amount of vapor added over the entire simulation time.

                       ! Update temperaure (Schmid 2021, equation (35) & (36))
                       t_phy(i, k, j) = t_tmp - (ls * (-snowsubl_limit) / (cp * rho_dry))

                       ! write out rate of temperature change
                       tempch(i, k, j) = (t_phy(i, k, j) - t_tmp) / dt ! temperature change rate (K s-1)

                       ! accumulate temperature change over time
                       tempchacc(i, k, j) = tempchacc(i, k, j) + tempch(i, k, j) * dt ! acc temperature change (K)

                    else ! no sublimation feedback
                       snowsubl_limit = snowsubl_limit_tracer
                    endif

                    ! write out sublimation rate
                    snowsubl(i, k, j) = snowsubl_limit / dt  ! sublimation rate (kg m-3 s-1)

                    ! accumulate snowsublimation
                    snowsublacc(i, k, j) = snowsublacc(i, k, j) + snowsubl_limit ! accumulated sublimation (kg m-3)

                    !update snow particle concentration
                    tracer(i, k, j, p_snowconc) = max(0.0, tracer(i, k, j, p_snowconc) + snowsubl_limit)

                 endif !end sigma > 0
              endif !end snowconc < d-8

              105 continue

           enddo !end k
        enddo !end i
     enddo !end j

  end subroutine snowsublimation
end module snowsub





