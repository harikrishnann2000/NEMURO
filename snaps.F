	subroutine write_output (loop_ind)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       subroutine to write snapshots of model output in netcdf format.
c       output can be restircted to selected variables in Makefile by
c       enabling curresponding options in the option list.
c 
c       author: vinu k. v.            email: vinu@ees.hokudai.ac.jp
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "netcdf.inc"
#include "size.h"
#include "param.h"
	
	taun1 = 2
	
c	write(*,*)'Writing output for loop ind ',loop_ind
c	write(*,'(e34.28,f12.8)') rVolInt
	do nt=1,nn
	write(Tr(nt),'(a,i2.2)') 'Tr',nt
	enddo
	
	
	
	
c
c      storing starting dimension values
c
	start4(1) = 1
	start4(2) = 1
	start4(3) = 1
	start4(4) = loop_ind
	start3(1) = 1
	start3(2) = 1
	start3(3) = loop_ind
	startp(1) = 1
	startp(2) = loop_ind
	startinv(1) = 1
	startinv(2) = loop_ind

	rdepth(1) = 0.0
	do ii=1,km
	rdepth(ii)=zdepth(ii)*(-1.0)
	enddo
	

c define undef value
	undef = 0.0

c NOTE*** the array ave(km) is used to find inter-facial depth as
c anomaly from the average depth of the respective layers. i.e. the
c at each level the suface elevation is written as h-h(average.of.domain).
c NOTE***
	
c
c assign lon, lat, depth and time arrays
	rcount=0
	do ii=1,imt
	rcount=rcount+0.5
	rlon(ii)=  ii*dxd + reflon - dxd
c	rlon(icount)=phi(icount)/deg2rad
	enddo
	rcount=-75.75
	do ii=1,jmt
	rcount=rcount+0.5
	rlat(ii)= theta(ii)/deg2rad
c	rlat(ii)= ii*dyd + (reflat-jmt*dyd/2-0.5) + dyd
c	rlat(icount)=theta(icount)/deg2rad
	enddo
	rtime(1)=loop*1.0*dt/day2sec
	



	 ip_region_value = p_region_value
c open netcdf file
	if (loop_ind .eq. 1) then
       ncid = nccre ('./OUTPUTS/snapshots_exp3.nc',
     &ncclob,rcode)
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	write(var3d_file_name(loop3d),'(a,a,i2.2,a)')
     &'./OUTPUTS/snaps_',var3d_char(loop3d),ip_region_value,'.nc'
       ncid_var3d(loop3d) = nccre (var3d_file_name(loop3d), 
     &   ncclob,rcode)
	enddo
#endif
	endif
	if (loop_ind .ne. 1) then 
	ncid = ncopn ('./OUTPUTS/snapshots_exp3.nc',
     & ncwrite,rcode)
	call ncredf (ncid, rcode)
	vtimeid = ncvid (ncid,'time',rcode)
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	write(var3d_file_name(loop3d),'(a,a,i2.2,a)')
     &'./OUTPUTS/snaps_',var3d_char(loop3d),ip_region_value,'.nc'
	ncid_var3d(loop3d) = ncopn (var3d_file_name(loop3d), 
     &   ncwrite, rcode)
	call ncredf (ncid_var3d(loop3d), rcode)
	vtimeid_var3d(loop3d) = ncvid (ncid_var3d(loop3d),'time',rcode)
	enddo
#endif
#ifdef snap_vmix
	vtid = ncvid (ncid, 'vmix', rcode)
#endif
#ifdef snap_velocity
	vuid = ncvid (ncid, 'u', rcode)
	vvid = ncvid (ncid, 'v', rcode)
	vwid = ncvid (ncid, 'w', rcode)
#endif
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	vvar3did(loop3d) = ncvid (ncid_var3d(loop3d), 
     &                      var3d_char(loop3d), rcode)
	enddo
#endif
#ifdef snap_tracer
	do nt=3,nn
	vtempid(nt) = ncvid (ncid, Tr(nt), rcode)
	enddo
#endif
#ifdef snap_salt
	vsaltid = ncvid (ncid, 'salt1', rcode)
#endif
#ifdef snap_temp
	vdensid = ncvid (ncid, 'temp', rcode)
#endif
#ifdef snap_vorticity
	vvortid = ncvid (ncid, 'pvort', rcode)
#endif
#ifdef snap_tau
	vtauxid = ncvid (ncid, 'tx' ,rcode)
	vtauyid = ncvid (ncid, 'ty' ,rcode)
#endif
#ifdef snap_particle
	vpxid = ncvid (ncid, 'xp', rcode)
	vpyid = ncvid (ncid, 'yp', rcode)
#endif
#ifdef snap_inverse_model
	vthetaid = ncvid (ncid, 'Th',rcode)
#endif
	go to 101
	endif

c define dimension
	lonid = ncddef (ncid,'longitude',imt,rcode)
	latid = ncddef (ncid,'latitude',jmt,rcode)
	idepthid = ncddef (ncid,'depth',km,rcode)
	itimeid = ncddef (ncid,'time',ncunlim,rcode)
	pid = ncddef (ncid,'array',NPP,rcode)

#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	lonid_var3d(loop3d) = ncddef (ncid_var3d(loop3d),
     &'longitude',imt,rcode)
	latid_var3d(loop3d) = ncddef (ncid_var3d(loop3d),
     &'latitude',jmt,rcode)
	idepthid_var3d(loop3d) = ncddef (ncid_var3d(loop3d),
     &'depth',km,rcode)
	itimeid_var3d(loop3d) = ncddef (ncid_var3d(loop3d),
     &'time',ncunlim,rcode)
	enddo
#endif

c store dimension id into arrays
	vardim(1) = lonid
	vardim(2) = idepthid
	vardim(3) = latid
	vardim(4) = itimeid
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	vardim_var3d(1,loop3d) = lonid_var3d(loop3d)
	vardim_var3d(2,loop3d) = idepthid_var3d(loop3d)
	vardim_var3d(3,loop3d) = latid_var3d(loop3d)
	vardim_var3d(4,loop3d) = itimeid_var3d(loop3d)
	enddo
#endif
	
	varpdim (1) = pid
	varpdim (2) = itimeid


	vardim3(1) = lonid
	vardim3(2) = latid
	vardim3(3) = itimeid

	vardim2(1) = lonid
	vardim2(2) = latid

	vardim33(1) = lonid
	vardim33(2) = idepthid
	vardim33(3) = latid

c define variable id
	vlonid = ncvdef (ncid,'longitude',ncfloat,1,lonid,rcode)
	vlatid = ncvdef (ncid,'latitude',ncfloat,1,latid,rcode)
	vdepthid = ncvdef (ncid,'depth',ncfloat,1,idepthid,rcode)
	vtimeid = ncvdef (ncid,'time',ncfloat,1,itimeid,rcode)
#ifdef var3d_snapshots
	do loop3d = 1,nvar3d
	vlonid_var3d(loop3d) = ncvdef (ncid_var3d(loop3d),
     &'longitude',ncfloat,1,lonid_var3d(loop3d),rcode)
	vlatid_var3d(loop3d) = ncvdef (ncid_var3d(loop3d),
     &'latitude',ncfloat,1,latid_var3d(loop3d),rcode)
	vdepthid_var3d(loop3d) = ncvdef (ncid_var3d(loop3d),
     &'depth',ncfloat,1,idepthid_var3d(loop3d),rcode)
	vtimeid_var3d(loop3d) = ncvdef (ncid_var3d(loop3d),
     &'time',ncfloat,1,itimeid_var3d(loop3d),rcode)
	enddo
#endif
c define attribute
        call ncaptc(ncid,vlatid,'units',ncchar,13,
     $    'degrees_north',rcode)
	call ncaptc(ncid,vlonid,'units',ncchar,12,
     $    'degrees_east',rcode)
	 call ncaptc(ncid,vdepthid,'units',ncchar,6,
     $    'meters',rcode)
#ifdef var3d_snapshots
	 do loop3d = 1, nvar3d
        call ncaptc(ncid_var3d(loop3d),vlatid_var3d(loop3d),
     $   'units',ncchar,13, 'degrees_north',rcode)
	call ncaptc(ncid_var3d(loop3d),vlonid_var3d(loop3d),
     $   'units',ncchar,12, 'degrees_east',rcode)
	 call ncaptc(ncid_var3d(loop3d),vdepthid_var3d(loop3d),
     $   'units',ncchar,6, 'meters',rcode)
	 enddo
#endif
#ifdef snap_vmix
	vtid = ncvdef (ncid,'vmix',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vtid,'units',ncchar,2,'CM',rcode)
	call ncapt(ncid,vtid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_velocity
	vuid = ncvdef (ncid,'u',ncfloat,4,vardim,rcode)
	vvid = ncvdef (ncid,'v',ncfloat,4,vardim,rcode)
	vwid = ncvdef (ncid,'w',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vuid,'units',ncchar,3,'M/S',rcode)
	call ncaptc(ncid,vvid,'units',ncchar,3,'M/S',rcode)
	call ncaptc(ncid,vwid,'units',ncchar,3,'M/S',rcode)
	call ncapt(ncid,vuid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vvid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vwid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	vardim_var3d_pass (1) = vardim_var3d(1,loop3d)
	vardim_var3d_pass (2) = vardim_var3d(2,loop3d)
	vardim_var3d_pass (3) = vardim_var3d(3,loop3d)
	vardim_var3d_pass (4) = vardim_var3d(4,loop3d)
	vvar3did(loop3d) = ncvdef (ncid_var3d(loop3d), 
     & var3d_char(loop3d),ncfloat,4,
     & vardim_var3d_pass, rcode)
	call ncaptc(ncid_var3d(loop3d), vvar3did(loop3d),
     &   'units',ncchar,3,  var3d_units(loop3d), rcode )
	call ncapt(ncid_var3d(loop3d), vvar3did(loop3d),
     &   'missing_value',ncfloat ,1,undef,rcode)
	enddo
#endif
#ifdef snap_particle
	vpxid = ncvdef (ncid,'xp',ncfloat,2,varpdim,rcode)
	vpyid = ncvdef (ncid,'yp',ncfloat,2,varpdim,rcode)
	call ncaptc (ncid,vpxid,'units',ncchar,6,'number',rcode)
	call ncaptc (ncid,vpyid,'units',ncchar,6,'number',rcode)
#endif
#ifdef snap_inverse_model
	vthetaid = ncvdef (ncid,'Th',ncdouble,2,varinvdim,rcode)
	call ncaptc (ncid,vthetaid,'units',ncchar,4,'flux',rcode)
#endif
#ifdef snap_tau
	vtauxid = ncvdef (ncid,'tx',ncfloat,3,vardim3,rcode)
	vtauyid = ncvdef (ncid,'ty',ncfloat,3,vardim3,rcode)
	call ncaptc(ncid,vtauxid,'units',ncchar,9,'Dynes/CM2',rcode)
	call ncaptc(ncid,vtauyid,'units',ncchar,9,'Dynes/CM2',rcode)
	call ncapt(ncid,vtauxid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vtauyid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_tracer
	do nt=3,nn
	vtempid(nt) = ncvdef (ncid,Tr(nt),ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vtempid(nt),'units',ncchar,8,'mole/m^3',rcode)
	call ncapt(ncid,vtempid(nt),'missing_value',ncfloat
     $    ,1,undef,rcode)
	enddo
#endif
#ifdef snap_salt
	vsaltid = ncvdef (ncid,'salt1',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vsaltid,'units',ncchar,3,'psu',rcode)
	call ncapt(ncid,vsaltid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_temp
	vdensid = ncvdef (ncid,'temp',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vdensid,'units',ncchar,3,'psu',rcode)
	call ncapt(ncid,vdensid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_vorticity
	vvortid = ncvdef (ncid,'pvort',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vvortid,'units',ncchar,13,'pot.vorticity',rcode)
	call ncapt(ncid,vvortid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
	vmaskid = ncvdef (ncid,'mask',ncfloat,3,vardim33,rcode)
	call ncaptc (ncid,vmaskid,'units',ncchar,4,'mask',rcode)
	call ncapt(ncid,vmaskid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	
	
#ifdef monthly_climatology
        call ncaptc(ncid,vtimeid,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid,vtimeid,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
#endif
        call ncaptc(ncid,vtimeid,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid,vtimeid,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
        call ncaptc(ncid_var3d(loop3d),
     $     vtimeid_var3d(loop3d),'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid_var3d(loop3d),
     $     vtimeid_var3d(loop3d),'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
	enddo
#endif
#ifdef daily_wind
        call ncaptc(ncid,vtimeid,'units',ncchar,21,
     $     stamp_day_daily(iday_start_snap),rcode)
        call ncaptc(ncid,vtimeid,'time_origine',ncchar,20,
     $     stamp_time_daily(iday_start_snap),rcode)
#endif
	
c
c Global attributes.
c
c        call ncaptc(ncid,nc_global,'title',ncchar,8,
c     $    'OTTM.1.1',rcode)
c        call ncaptc(ncid,nc_global,'label',ncchar,36,
c     $    'Simulations using interannual inputs',rcode)
c        call ncaptc(ncid,nc_global,'history',ncchar,28,
c     $    'generated by module: snaps.F',rcode)
	
c leaving in define mode
	call ncendf (ncid, rcode)
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	call ncendf (ncid_var3d(loop3d), rcode)
	enddo
#endif
c puting values
	call ncvpt (ncid,vlonid,1,imt,rlon,rcode)
	call ncvpt (ncid,vlatid,1,jmt,rlat,rcode)
	call ncvpt (ncid,vdepthid,1,km,rdepth,rcode)
	call ncvpt (ncid,vmaskid,start3,count33,lmask,rcode)
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	call ncvpt (ncid_var3d(loop3d),vlonid_var3d(loop3d),
     &   1,imt,rlon,rcode)
	call ncvpt (ncid_var3d(loop3d),vlatid_var3d(loop3d),
     &   1,jmt,rlat,rcode)
	call ncvpt (ncid_var3d(loop3d),vdepthid_var3d(loop3d),
     &   1,km,rdepth,rcode)
	enddo
#endif
 101    continue
	 if ( loop_ind .gt. 1) call ncendf (ncid, rcode)
	call ncvpt (ncid,vtimeid,loop_ind,1,rtime(1),rcode)
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	 if ( loop_ind .gt. 1) call ncendf (ncid_var3d(loop3d), rcode)
	call ncvpt (ncid_var3d(loop3d) ,vtimeid_var3d(loop3d),
     &  loop_ind, 1, rtime(1), rcode)
	enddo
#endif
#ifdef snap_vmix
	do ii = 1,imt
	do kk = 1,km
	do jj = 1,jmt
	 var3d(ii,kk,jj) = vmix_coeff(ii,kk,jj)
	enddo
	enddo
	enddo
	call ncvpt (ncid,vtid,start4,count4,var3d,rcode)
#endif
#ifdef snap_velocity
	do ii = 1,imt
	do kk = 1,km
	do jj = 1,jmt
	 var3d(ii,kk,jj) = u(ii,kk,jj,taup) !*mask(ii,kk,jj)
	enddo
	enddo
	enddo
	call ncvpt (ncid,vuid,start4,count4,var3d,rcode)
	do ii = 1,imt
	do kk = 1,km
	do jj = 1,jmt
	 var3d(ii,kk,jj) = v(ii,kk,jj,taup)
	enddo
	enddo
	enddo
	call ncvpt (ncid,vvid,start4,count4,var3d,rcode)
	do ii = 1,imt
	do kk = 1,km
	do jj = 1,jmt
	 var3d(ii,kk,jj) = w(ii,kk,jj,taup)
	 var3d(1,kk,jj) = w(2,kk,jj,taup)
	 var3d(imt,kk,jj) = w(imt-1,kk,jj,taup)
	enddo
	enddo
	enddo
	call ncvpt (ncid,vwid,start4,count4,var3d,rcode)
#endif
#ifdef snap_tracer
	do nt=3,nn
	do ii = 1,imt
	do kk = 1,km
	do jj = 1,jmt
	 var3d(ii,kk,jj) = t(ii,kk,jj,nt,taun) !*lmask(ii,kk,jj)
	enddo
	enddo
	enddo
	call ncvpt (ncid,vtempid(nt),start4,count4,var3d,rcode)
	enddo
#endif
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	do ii = 1,imt
	do kk = 1,km
	do jj = 1,jmt
	 var3d(ii,kk,jj) = var3d_recv(ii,kk,jj,loop3d)*lmask(ii,kk,jj)
	enddo
	enddo
	enddo
	call ncvpt (ncid_var3d(loop3d), vvar3did(loop3d),
     &               start4,count4,var3d,rcode)
	enddo
#endif
#ifdef snap_salt
	do ii=1,imt
	do kk=1,km
	do jj=1,jmt
	 var3d(ii,kk,jj) = t(ii,kk,jj,2,taun)
	enddo
	enddo
	enddo
	
	call ncvpt (ncid,vsaltid,start4,count4,var3d,rcode)
#endif
#ifdef snap_temp
	do ii=1,imt
	do kk=1,km
	do jj=1,jmt
	 var3d(ii,kk,jj) = t(ii,kk,jj,1,taun)
	enddo
	enddo
	enddo
	call ncvpt (ncid,vdensid,start4,count4,var3d,rcode)
#endif
#ifdef snap_inverse_model
	call ncvpt (ncid,vthetaid,startinv,countinv,Thetainv,rcode)
#endif

c close netcdf file
	call ncclos (ncid, rcode)
#ifdef var3d_snapshots
	do loop3d = 1, nvar3d
	call ncclos (ncid_var3d(loop3d), rcode)
	enddo
#endif

	
	return
	
	end
