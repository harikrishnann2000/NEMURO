        subroutine relax_nutrient(loop, no3_relax_obs,
     &                           sioh4_relax_obs, dt)
C
C     Declarations
        implicit none
        real    no3_relax_obs(51),sioh4_relax_obs(51)
        integer loop, i , k, lpmb, dt
	integer mon, lpm
	integer dpm
	integer day2sec
	character*3 cpm
	dimension cpm(12), dpm(13)
	character*106 file_no3
	character*108 file_sioh4
C
C     Monthly climatology and time conversion
	data (dpm(i),i=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
	data (cpm(i),i=1,12)/'jan','feb','mar','apr','may','jun',
     &                       'jul','aug','sep','oct','nov','dec'/
	data day2sec/86400/
C
C     Shared month-tracking state
	common mon, lpm, lpmb
C
C     Select the climatological month for the current timestep
	if (loop .eq. 1) then
	mon = 1
	lpmb = 1
	lpm =  dpm(mon)*day2sec/dt
	endif
	if (loop .gt. lpm) then
	mon = mon + 1
	lpmb = lpm
	lpm = lpm + dpm(mon)*day2sec/dt
	endif
	if (mon .gt. 12) then
	mon = 1
	lpmb = lpm
	lpm = lpm + dpm(mon)*day2sec/dt
	endif
C
C     Construct monthly climatology filenames
        if ((loop.eq. 1)) then
	write (file_no3,'(a,a,a)')
     &"./INPUT/no3_clim_extract_",cpm(mon),"_argo.txt"
	write (file_sioh4,'(a,a,a)')
     &"./INPUT/sioh4_clim_extract_",cpm(mon),"_argo.txt"
	endif
        if ((loop.gt. lpmb).and.(loop.le. lpm)) then
	write (file_no3,'(a,a,a)')
     &"./INPUT/no3_clim_extract_",cpm(mon),"_argo.txt"
	write (file_sioh4,'(a,a,a)')
     &"./INPUT/sioh4_clim_extract_",cpm(mon),"_argo.txt"
        end if
C
C     Read the 51-level nitrate and silicate profiles
        open(50,file=file_no3)
        open(60,file=file_sioh4)
        do k=1,51
        read(50,*)no3_relax_obs(k)
        read(60,*)sioh4_relax_obs(k)
        end do
        close(50)
        close(60)
	return
        end
