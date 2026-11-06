	subroutine implicit_vmix (KH, dz1, dt, km, t1, t2)

        real t1(km), dz1, t2(km), dt, KH(km)
	double precision t(km), dz(0:km+1), rA(km,km)
        double precision tM(km)
	double precision dtts, tmask(km)
	double precision rD1(km),  rD2(km)
	double precision rDDu(km), rDDl(km)
        double precision dzfact(km,km), dzfact_rev(km,km)
        double precision rheat1, rheat2
	real diffuse(km)


	
c
c initialize data
c

         dtts = dt
         do i=1,km
         dz(i) = dz1
         t(i) = t1(i)
         diffuse(i) = KH(i)
         enddo
         dz(0) = dz(1)
         dz(km+1) = dz(km)



        
	 do i=1,km
	 rDDu(i) =  1.0*dtts/((dz(i-0)+dz(i)))
	 rDDl(i) =  1.0*dtts/((dz(i)+dz(i+0)))
         if (i-1 .eq. 0) then
         rDDu(i) = 1.0*dtts /(dz(i-0)+dz(i))
         rDDl(i) = 1.0*dtts /(dz(i)+dz(i+0))
         endif
         rD1(i) = diffuse(i-1)/((dz(i-1)+dz(i))*0.5)
         rD2(i) = diffuse(i)/((dz(i)+dz(i+1))*0.5)
	 enddo
	
         rD1(1) = 0.0
         rD2(km) = 0.0
	
	
c
c form the matrix
c	
	 
	 do i=1,km
	 do j=1,km
         dzfact(i,j) = dz(i)/dz(j)
         dzfact_rev(i,j) = dz(j)/dz(i)
	  if ( j .eq. i-1 ) then 
	    rA(i,j) = -rDDu(i)*rD1(i)* dzfact(i,j)
	  elseif ( j .eq. i-0 ) then 
	    rA(i,j) = (1.0 + (rDDu(i)*rD1(i)+rDDl(i)*rD2(i)) )*
     &                                  dzfact(i,j)
	  elseif ( j .eq. i+1 ) then
	  rA(i,j) = -rDDl(i)*rD2(i) * dzfact(i,j)
	  else
	    rA(i,j) = 0.0
	  endif
	 enddo
	 enddo
	
	
	 iC = 1
	
	
	 call inverse (km, rA)
	
         do j=1,km
         do i=1,km
         rmid = dz(i)
         rA(j,i) = rA(j,i)*dzfact_rev(j,i)
         enddo
         enddo

	 
	 call multip (km, km, iC, rA, t, tM)  
	
         do i=1,km
         t2(i) = tM(i)
         enddo
	
	return
	end
	
	subroutine inverse (n, A)
        double precision A(n,n)
        double precision work(n)
        integer m, n, lda, info, lwork, ipiv(n)
        m = n
        lda = n
        lwork = n


        call dgetrf (m,n,A,lda,ipiv,info)

        call dgetri (n,A,lda,ipiv,work,lwork,info)

        return
        stop
        end
	
	subroutine multip (M, K, N, A, B, AB)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c a subroutine to multiply two matrix.
c
c        Usage: call multip (M, K, N, A, B, AB)
c
c        input: M = number of rows
c               K = number of coloumns
c               N = number of coloumns
c               A = array of dimension (M,K)
c               B = array of dimension (K,N)
c        output:AB = array of dimension (M,N)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	double precision A(M,K), B(K,N), AB(M,N)
        double precision rsum

        do i=1,M
        do ii=1,N
         rsum = 0.0
        do j=1,K
         rsum = rsum + A(i,j)*B(j,ii)
        enddo
         AB(i,ii) = rsum
        enddo
        enddo

        return
        stop
        end
