	subroutine porder(kn,dist, ntr,pos,nndist)
        integer kn, ntr
	double precision  dist(ntr), nndist(kn)
	integer pos(kn)
	   do 50 j = 1, ntr
	      if(j .le. kn) then
	         do 30 k = 1, j-1
	            if(dist(j) .lt. nndist(k)) then
	               do 20 k1 = j-1, k, -1
	                  nndist(k1+1) = nndist(k1)
20	                  pos(k1+1) = pos(k1)
	               nndist(k) = dist(j)
	               pos(k) = j
	               goto 50
	            endif
30	         continue
	         nndist(j) = dist(j)
	         pos(j) = j
	      else
	         if(dist(j) .ge. nndist(kn)) go to 50
	         do 40 k = 1, kn
	            if(dist(j) .lt. nndist(k)) then
	               do 35 k1 = kn-1, k, -1
	                  nndist(k1+1) = nndist(k1)
35	                  pos(k1+1) = pos(k1)
	               nndist(k) = dist(j)
	               pos(k) = j
	               goto 50
	            endif
40	         continue
  	      endif
50	   continue
	   return
	   end


C Output from Public domain Ratfor, version 1.0
      subroutine knnimp(x,ximp,p,n,imiss,irmiss,kn,workp,workn,iworkp,iw
     &orkn)
      integer p,n,i,j,k,kn,m
      integer imiss(p,n),irmiss(p)
      double precision x(p,n),ximp(p,n)
      integer iworkp(p),iworkn(n)
      double precision workp(p),workn(p)
      m=kn+1
      do23000 i=1,p
      if(irmiss(i) .eq. 0)then
      goto 23000
      endif
      do23004 j=1,n
      workn(j)=x(i,j)
      iworkn(j)=imiss(i,j)
23004 continue
23005 continue
      call misdis(workn,x,p,n,iworkn,imiss,workp,iworkp)
      call porder(m,workp,p,iworkp,workn)
      call misave(x,workn,p,n,iworkn,imiss,iworkp(2),kn)
      do23006 k=1,n
      if(iworkn(k).eq.0)then
      goto 23006
      endif
      ximp(i,k)=workn(k)
      if(iworkn(k).eq.2)then
      imiss(i,k)=2
      endif
23006 continue
23007 continue
23000 continue
23001 continue
      return
      end
      subroutine misdis(x0,x,p,n,imiss0,imiss,dis,iworkp)
      integer p,n
      integer imiss(p,n),iworkp(p),imiss0(n)
      double precision x0(n),x(p,n),dis(p),dismax
      dismax=1d10
      do23012 j=1,p 
      dis(j)=0d0
      iworkp(j)=0
23012 continue
23013 continue
      do23014 k=1,n
      if(imiss0(k).gt.0)then
      goto 23014
      endif
      do23018 j=1,p 
      if(imiss(j,k).gt.0)then
      goto 23018
      endif
      dis(j)= dis(j)+(x0(k)-x(j,k))**2
      iworkp(j)=iworkp(j)+1
23018 continue
23019 continue
23014 continue
23015 continue
      do23022 j=1,p
      if(iworkp(j).gt.0)then
      dis(j)=dis(j)/iworkp(j)
      else
      dis(j)=dismax
      endif
23022 continue
23023 continue
      return
      end
      subroutine misave(x,x0,p,n,imiss0,imiss,index,m)
      integer p,n,m,j,jj,k
      integer imiss0(n),imiss(p,n),index(m),ktot
      double precision x(p,n),x0(n)
      do23026 k=1,n
      x0(k)=0d0
      if(imiss0(k).eq.0)then
      goto 23026
      endif
      ktot=0
      do23030 j=1,m
      jj=index(j)
      if(imiss(jj,k).gt.0)then
      goto 23030
      endif
      x0(k)=x0(k)+x(jj,k)
      ktot=ktot+1
23030 continue
23031 continue
      if(ktot.gt.0)then
      x0(k)=x0(k)/ktot
      else
      imiss0(k)=2
      endif
23026 continue
23027 continue
      return
      end
C Output from Public domain Ratfor, version 1.0
      subroutine twomis(x,p,n,imiss,x0,imiss0,maxit,eps,istart,clust, ns
     &ize,dist,ratio,iter,iworkp,iworkn)
      integer p,n,imiss(p,n),imiss0(n,2),maxit,istart(2),clust(p,2),iwor
     &kp(p)
      double precision x(p,n),x0(n,2),eps,dist(p,2)
      integer nsize(2),iworkn(n)
      integer iter,imax
      double precision ratio,dold,dnew
      if(maxit.lt.1)then
      maxit=5
      endif
      do23002 i=1,2
      do23004 j=1,n
      x0(j,i)=x(istart(i),j)
      imiss0(j,i)=imiss(istart(i),j)
23004 continue
23005 continue
23002 continue
23003 continue
      iter =0
      ratio = 1d1
23006 if((iter.lt.maxit).and.(ratio.gt.eps))then
      iter=iter+1
      do23008 i=1,2
      call misdis(x0(1,i),x,p,n,imiss0(1,i),imiss,dist(1,i),iworkp)
      nsize(i)=0
23008 continue
23009 continue
      dnew=0d0
      do23010 j=1,p 
      if(dist(j,1).lt.dist(j,2))then
      imax=1
      else
      imax=2
      endif
      nsize(imax)=nsize(imax)+1
      clust(nsize(imax),imax)=j
      dnew=dnew+dist(j,imax)
23010 continue
23011 continue
      if(dnew.eq.0d0)then
      goto 23007
      endif
      if(iter.eq.1)then
      dold=dnew*1d1
      endif
      ratio=dabs(dnew-dold)/dold
      dold=dnew
      do23018 i=1,2
      do23020 j=1,n 
      iworkn(j)=1
23020 continue
23021 continue
      call misave(x,x0(1,i),p,n,iworkn,imiss,clust(1,i),nsize(i))
      do23022 j=1,n
      if(iworkn(j).eq.2)then
      imiss0(j,i)=1
      else
      imiss0(j,i)=0
      endif
23022 continue
23023 continue
23018 continue
23019 continue
      goto 23006
      endif
23007 continue
      return
      end
