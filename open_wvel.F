	include "netcdf.inc"
	parameter (lon=1, lat=1, idepth=51, itime=615)
	integer ncid, vdatid
	integer start(4)
	integer count(4)
	
	real dat(idepth,itime)
	double precision dat1(lon,lat,idepth,itime)
	data start/1,1,1,1/
	data count/lon,lat,idepth,itime/
c
c open netcdf file
	ncid = ncopn ('wvel_extracted_argo11.nc', ncnowrit, rcode)
c get variable id
	vdatid = ncvid (ncid,'W_VEL_FILL',rcode)
c
c get variable values
	call ncvgt (ncid, vdatid, start, count, dat1, rcode)
c
c close netcdf file
	 call ncclos(ncid, rocode)
c check data for land values
c
	do k=1,idepth
        do l=1,itime
	if (dat1(1,1,k,l) .lt. -999.999) then 
	 dat(k,l) = 0.0
	 else 
	 dat(k,l) = dat1(1,1,k,l) ! unit is m/s
	endif
	enddo
	enddo

c	do ii = 1,itime
c        write (*,*) (dat(ii,jj),jj=1,51)
c        enddo

c write into .bin format
	open
     &(unit=9,file='./INPUT/wvel_extracted_argo11.bin',
     &access='direct',form='unformatted',recl=51)
	irec = 1
	do ii = 1,itime
        write(9,rec=irec) ((dat(jj,ii))/10,jj=1,51)
        irec=irec+1
	enddo
	close (9)
	end
