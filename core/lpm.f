c----------------------------------------------------------------------
      subroutine lpm_stokes_particles_solver
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer             stage,nstage
      common /tstepstage/ stage,nstage

      nstage_part = 3
      if (abs(time_integ) .eq. 2) nstage_part = 1

      if(ipart_moving.eq.0) nstage_part = 1    ! Stationary particle, no need to iterate the particle position    
      if(ipart_moving.eq.2) nstage_part = 1    ! Forced moving

      if (istep.eq.0) then
         call lpm_init(0)
      else

         call set_tstep_coef_part(dt) ! in nek5000 with rk3
         do stage=1,nstage_part
            call lpm_usr_particles_solver
         enddo
         call compute_phig_qtl(dt,usrdiv) ! nek5000 (see Zwick 2018)
      endif

      if(mod(istep,iostep).eq.0.or. istep.eq.1) then
         call lpm_usr_particles_io
      endif

      return
      end

c----------------------------------------------------------------------
c     setup routines
c----------------------------------------------------------------------
      subroutine lpm_init(idum)
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LPM'

c     common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
c     integer gfirst, inoassignd, resetFindpts, pload

      icmtp = idum

      call set_part_pointers
c     call read_particle_input_par ! for lb code since no par file
      call set_bounds_box
      call set_part_params 
      call place_particles

      call update_particle_location   

      call set_check_spl_params 
      call move_particles_inproc

      if(ibm_debug_bin.eq.1 .and. n.ge.1) then
         iqt = 0
         ipt = 0
         do i = 1,n
            if(ipart(jrole,i).eq.1) iqt = iqt + 1 
            if(ipart(jrole,i).eq.2) ipt = ipt + 1 
         enddo
         if(mod(istep,20).eq.0)
     >    write(6,'(I4,A,3I6)') nid," contains ",n,ipt,iqt 
      endif
      
      if (two_way.gt.1) then
         if(ifibm .ge. 1) then
            call set_check_spl_params_ibm
            call compute_neighbor_el_proc
            call save_ibm_neighbor_bin
            if(nid.eq.0) print*, "===Save IBM Bin Structure Done==="
         endif
         ! 
         call set_check_spl_params
         call compute_neighbor_el_proc

         if (nid.eq.0) write(6,502)    "Queen  Bin size: ", 
     >  nlist_ibm,ndxgp_ibm,ndygp_ibm,ndzgp_ibm,
     >  rdxgp_ibm,rdygp_ibm,rdzgp_ibm
         if (nid.eq.0) write(6,502)    "Worker Bin size: ", 
     > nlist, ndxgp,ndygp,ndzgp, rdxgp,rdygp,rdzgp
 502  FORMAT(A20,4I6,3f12.5)     

!         if(ibm_debug_bin.eq.1)
        if(ibm_debug_bin.eq.1) call output_bin_structure
         !
         call create_extra_particles
         call send_ghost_particles
         call send_ghost_particles_ibm

         if(ibm_debug_worker.ge.1 .and. nfptsgp_ibm .gt.0) then
           write(6,'(A,I6,A,I6,A)')'Rank ', nid," receives ",
     $           nfptsgp_ibm," Ghost Queen"
           do i = 1, nfptsgp_ibm
              write(6,2025)"Ghost Queen in rank:", nid,
     $             ", from rank " , iptsgp_ibm(jgpps,i),
     $             ", jgp_back = "   , iptsgp_ibm(jgp_back,i),
     $             ", born Queen id ="   , iptsgp_ibm(jgp_queen,i),
     $             ", Pid1(born rank)=", iptsgp_ibm(jgp_pid1,i)
           enddo
         endif

         call spread_props_grid           
      endif

 2025 format(A,I6,A,I6,A,I6,A,I6,A,I6)

      call interp_props_part_location 

      if (time_integ .lt. 0) call pre_sim_collisions ! e.g., settling p
      if (time_integ .lt. 0) write(*,*) "call pre collision"

c     resetFindpts = 0
c     call computeRatio
c     call reinitialize
c     call printVerify

      return
      end

c----------------------------------------------------------------------
      subroutine place_particles
c
c     Place particles in this routine, also called for injection
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer icalld
      save    icalld
      data    icalld  /-1/

      real     unif_random,unif_random_norm,unif_random_cyl
     >        ,unif_random_sphere
      external unif_random,unif_random_norm,unif_random_cyl,
     >         unif_random_sphere

      integer        ndef
      common         ndef
      icalld = icalld + 1

      if(ipart_restartr .eq. 0) then

         if( IBM_Particle_shape .ge. 1 ) then   !!! IBM
            
c        set nwe based on ranks and IBM particle number

         n_IBMpart   = int(num_of_IBMpart/np)            ! num. part per proc
         nw_tmp      = iglsum(n_IBMpart,1)
         ndef        = num_of_IBMpart - nw_tmp         
         if( nid .lt. ndef  ) n_IBMpart = n_IBMpart + 1   ! add the remainder to proc before          
         if( n_IBMpart.eq.0 ) return

c        distribute IBM particle distribution & function
         if ( IBM_Particle_shape .eq. 1 ) then
            call cyl_marker_distribution
         else if  ( IBM_Particle_shape .eq. 2 ) then
            call sphere_marker_distribution
         else
            call lpm_user_marker_distribution
         endif 

         if(nid.eq.0) print*,"nid=",nid,",IBM_part=", n_IBMpart,
     $                                      ",nwe =", nwe         

         !  set initial ipart/rpart values Queen/worker markers
         n = 0
         do i_pt_part = 1,nwe
            n = n + 1            
            rpart(jrhop,n) = rho_p                               ! particle density 
            rpart(jtaup,n) = rpart(jdp,n)**2*rho_p/18.0d+0/mu_0  ! particle time scale

            rpart(jspl,n)  = rspl ! super particle loading
            rpart(jrpe,n)  = rpart(jspl,n)**(1./3.)*rpart(jdp,n)/2.
            rpart(jvol,n)  = rpart(jspl,n)*rpart(jvol,n)
         
            rpart(jtemp, n) = tp_0                               ! intial particle temp
            rpart(jtempf,n) = tp_0                               ! intial fluid temp (overwritten)
            rpart(jrho,n)   = param(1)                           ! initial fluid density (overwritten interp)

c           set global particle id (3 part tag)
            ipart(jpid1,n) = nid                                 ! born rank # 
            ipart(jpid2,n) = i_pt_part                           ! born marker #
            ipart(jpid3,n) = icalld                              ! born time step # 

            if( i_pt_part .le. n_IBMpart) then
c           queen_marker
               n_ibm_id          = i_pt_part                     !!! queen_id
               ipart(jrole,n)    = 1                             !!! role = 0- [point]; [1] -queen; 2- worker;
               ipart(jqueen,n)   = n_ibm_id                      !!! ibm_part_id = queen_id
               ipart(jworker1,n) = (n_ibm_id - 1) * n_l(1) + 1   !!! for equal_distribution only
               ipart(jnlm,n)     = n_l(n_ibm_id) 
            else
c           worker_marker
               n_ibm_id          = int((i_pt_part - n_IBMpart-1)
     $              / n_l(1)) + 1
               ipart(jrole,n)    = 2                             !!! x 0- [point] [1] -queen 2- worker
               ipart(jqueen,n)   = n_ibm_id                      !!! ibm_part_id = queen_id
               ipart(jworker1,n) = n_IBMpart + 
     $              (n_ibm_id-1)*n_l(1) + 1                      !!! for equal_distribution only
               ipart(jnlm,n)     = n_l(n_ibm_id) 
            endif

         if(nid.eq.0) write(6,2021) n,ipart(jrole,n),ipart(jqueen,n),
     $        ipart(jworker1,n), ipart(jnlm,n), 
     $        (rpart(jx+j,n),j=0,2)
     $    ,rpart(jdp,n),rpart(jvol,n),(rpart(jv0+j,n),j=0,2)
         enddo

 2021 format('Marker #',5I6,3E12.4,', dp=',E12.4,
     $     ', vol=',E12.4,', vel=',3E12.4)

c--------------------------------------------------------------------------------
         else ! Original code other than IBM 
        
         if(nid.eq.0) write(6,*) 'Place point particles (PP)'

c        correct nwe if discrepancy on rank 0
         nwe         = int(nw/np)                ! num. part per proc
         nw_tmp      = iglsum(nwe,1)
         ndef        = nw - nw_tmp
         if (nid .lt. ndef) nwe = nwe + 1
         
c        main loop to distribute particles
         do i_pt_part = 1,nwe
            n = n + 1
  754 continue
            ! Error checking
            if (n.gt.llpart)then 
               if (nid.eq.0)
     >            write(6,*)'Not enough space to store more particles'
               call exitt
            endif


            ! sphere or cylinder
            if (rxco(1) .gt. -1E7) then

               ! distribute in cylinder 
               if (rxco(6) .gt. -1E7) then
               rnx = rxco(1)
               rny = rxco(2)
               rnz = rxco(3)

               rx0 = rxco(4)
               ry0 = rxco(5)
               rz0 = 1.0
               if (if3d) rz0 = rxco(6)

               rmag = sqrt(rnx**2 + rny**2 + rnz**2)
               rnx = rnx/rmag
               rny = rny/rmag
               rnz = rnz/rmag

               rin  = rxco(7)
               rout = rxco(8)
               rht  = 0.
               if (if3d) rht  = rxco(9)

               rrad = unif_random_cyl(rin,rout)
               rthet= unif_random(0.,2.*pi)
               rxtr = unif_random(-rht/2.,rht/2.)

               do j=0,2
               ! x cylinder
               if ( abs(rnx).gt.abs(rny).and.abs(rnx).gt.abs(rnz)) then
                  if (j.eq. 0) rdum = rx0 + rxtr
                  if (j.eq. 1) rdum = ry0 + rrad*cos(rthet)
                  if (j.eq. 2) rdum = rz0 + rrad*sin(rthet)
               endif
               ! y cylinder
               if ( abs(rny).gt.abs(rnx).and.abs(rny).gt.abs(rnz)) then
                  if (j.eq. 0) rdum = rx0 + rrad*sin(rthet)
                  if (j.eq. 1) rdum = ry0 + rxtr
                  if (j.eq. 2) rdum = rz0 + rrad*cos(rthet)
               endif
               ! z cylinder
               if ( abs(rnz).gt.abs(rnx).and.abs(rnz).gt.abs(rny)) then
                  if (j.eq. 0) rdum = rx0 + rrad*cos(rthet)
                  if (j.eq. 1) rdum = ry0 + rrad*sin(rthet)
                  if (j.eq. 2) rdum = rz0 + rxtr
               endif
                  rpart(jx +j,n) = rdum
                  rpart(jx1+j,n) = rdum
                  rpart(jx2+j,n) = rdum
                  rpart(jx3+j,n) = rdum
               enddo

               ! distribute sphere
               else
               rx0 = rxco(1)
               ry0 = rxco(2)
               rz0 = 1.0
               if (if3d) rz0 = rxco(3)
               rin  = rxco(4)
               rout = rxco(5)

               rrad  = unif_random_sphere(rin,rout)
               rthet1= unif_random(0.,1.)
               rthet1= acos(2.*rthet1 - 1) !tricky tricky
               rthet2= unif_random(0.,2.*pi)

               do j=0,2
                  if (j.eq. 0) rdum = rx0 + rrad*sin(rthet1)*cos(rthet2)
                  if (j.eq. 1) rdum = ry0 + rrad*sin(rthet1)*sin(rthet2)
                  if (j.eq. 2) rdum = rz0 + rrad*cos(rthet1)
                  rpart(jx +j,n) = rdum
                  rpart(jx1+j,n) = rdum
                  rpart(jx2+j,n) = rdum
                  rpart(jx3+j,n) = rdum
               enddo
               endif

               ! check if box defined and if it is, particles must fit 
               ! in that box even if spere or cylinder distribution.
               ! If not, try placing again
               if (rxbo(1,1) .gt. -1E7) then
                  do j=0,ldim-1
                     if ((rpart(jx+j,n) .gt. rxbo(2,j+1)) .or.
     >                   (rpart(jx+j,n) .lt. rxbo(1,j+1))) then
                        goto 754
                     endif
                  enddo
               endif
                  
            ! distribute in box
            else
               do j=0,2
                  rpart(jx +j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
                  rpart(jx1+j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
                  rpart(jx2+j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
                  rpart(jx3+j,n) = unif_random(rxbo(1,j+1),rxbo(2,j+1))
               enddo
            endif

c           set some rpart values for later use
            if (dp_std .gt. 0) then
               rpart(jdp,n) = unif_random_norm(dp(1),dp_std)
            else
               rpart(jdp,n) = unif_random(dp(1),dp(2))
            endif
            rpart(jtaup,n) = rpart(jdp,n)**2*rho_p/18.0d+0/mu_0  ! particle time scale
            rpart(jrhop,n) = rho_p                               ! particle density 
            rpart(jvol,n) = pi*rpart(jdp,n)**3/6.      ! particle volume
            rpart(jspl,n)  = rspl                                ! super particle loading
            rpart(jrpe,n) = rpart(jspl,n)**(1./3.)*rpart(jdp,n)/2.
            rpart(jvol,n) = rpart(jspl,n)*rpart(jvol,n)
         
            rpart(jtemp,n)  = tp_0                               ! intial particle temp
            rpart(jtempf,n) = tp_0                               ! intial fluid temp (overwritten)
            rpart(jrho,n)   = param(1)                           ! initial fluid density (overwritten interp)
         
c           set global particle id (3 part tag)
            ipart(jpid1,n) = nid 
            ipart(jpid2,n) = i_pt_part
            ipart(jpid3,n) = icalld

            write(6,*)"Particle", n, rpart(jvol,n),rpart(jx+0,n) 
     $           ,rpart(jx+1,n) ,rpart(jx+2,n), rpart(jdp,n)

         enddo

         endif ! Marker distribution

      else
         ! read in data
         nread_part = 1
         do j=1,nread_part
            call read_parallel_restart_part
         enddo

      endif


      ! force 2d z to be 1
      if (.not. if3d) then
         do i=1,n
            rpart(jx +2,i) = 1.
            rpart(jx1+2,i) = 1.
            rpart(jx2+2,i) = 1.
            rpart(jx3+2,i) = 1.
         enddo
      endif

      return
      end

      
c----------------------------------------------------------------------      
      subroutine cyl_marker_distribution
c----------------------------------------------------------------------      

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer :: idist_marker
      integer :: n_markers       ! for each circular disc
      integer :: num_of_cylinders
      real cyl_center(3)
      real cyl_diam, cyl_length

c------------------------------------------------------------------------      
      ! assigned data
c      num_of_cylinders = 1      ! for 1 cylinder
c      cyl_center(1) = 0.0
c      cyl_center(2) = 0.02
c      cyl_center(3) = 1.0
c      cyl_diam = 0.02
c      idist_marker = 1
c      call read_cylinder_coornidate

      num_of_cylinders = num_of_IBMpart

      cyl_length = 1.0

      do nn = 1, num_of_cylinders
         cyl_center(1) = ibm_center(nn,1) 
         cyl_center(2) = ibm_center(nn,2) 
         cyl_center(3) = ibm_center(nn,3) 
         cyl_diam      = ibm_diam(nn)
       enddo

       if(if3d) cyl_length = ibm_diam(1)
       if(if3d) num_of_cylinders = int( cyl_length /
     > (ibm_diam(1)/n_dh(1))) + 1

      dh = cyl_length/n_dh(1) 
      write(*,*)"Place" ,num_of_cylinders, " (layer) Cylinders"      
      write(*,*)"Center=",cyl_center(1),cyl_center(2),cyl_center(3)
      write(*,*)"Diameter=",cyl_diam

c------------------------------------------------------------------------
      n_markers = nwe / num_of_cylinders 
      rpi = 4.0*atan(1.0)
      write(*,*)"Markers number = ", n_markers

      do jj = 1, num_of_cylinders

         cyl_rad = cyl_diam/2.0

         if(if3d) cyl_center(3) = (-cyl_length/2.) + (jj-1) * dh 

         do kk = 1, n_markers
            nn = (jj-1) * n_markers + kk

            call cyl_dist(n_markers, nn,cyl_rad, cyl_center,
     &           rpart(jx, nn), rpart(jx+1, nn), rpart(jx+2,nn))
            call cyl_dist(n_markers, nn,cyl_rad, cyl_center,
     &           rpart(jx1,nn), rpart(jx1+1,nn), rpart(jx1+2,nn))
            call cyl_dist(n_markers, nn,cyl_rad, cyl_center, 
     &           rpart(jx2,nn), rpart(jx2+1,nn), rpart(jx2+2,nn))
            call cyl_dist(n_markers, nn,cyl_rad, cyl_center, 
     &           rpart(jx3,nn), rpart(jx3+1,nn), rpart(jx3+2,nn))

            rpart(jdp,nn) =  dp(1)

            if(IBM_marker_shape .eq. 1) then
               rpart(jvol,nn) =  rpart(jdp,nn)**2
            else if (IBM_marker_shape .eq. 2) then
              rpart(jvol,nn) =  rpart(jdp,nn)**3
            else
               rpart(jvol,nn) = rpi*rpart(jdp,nn)**3/6. !
            endif
!     rpart(jdp,nn)**2    ! for 2D ibm test (Uhlmann) only
            !     pi*rpart(jdp,nn)**3/6.      ! particle volume

            rpart(jspl,nn)  = rspl ! super particle loading
            rpart(jrpe,nn) = rpart(jspl,nn)**(1./3.)*rpart(jdp,nn)/2.
            rpart(jvol,nn) = rpart(jspl,nn)*rpart(jvol,nn)

c            write(*,*)"Markers volume = ", nn, rpart(jvol,nn)

         enddo
      enddo
      
c      write(*,*)"Markers width  = ", nn, rpart(jdp,nn)

      return
      end
c----------------------------------------------------------------------
      subroutine cyl_dist(n_markers,nn,rad,center, x,y,z)

c   distribute Lagrange markers for a cylinder geometry
        integer n_markers, nn
        real  rad, center(3)
        real  x,y,z
        
        rpi = 4.*atan(1.0)
        deg = 2. * rpi / n_markers * (nn-1)
        x = center(1) + rad * cos(deg)
        y = center(2) + rad * sin(deg)
        z = center(3)        
      return
      end
        


c----------------------------------------------------------------------
      subroutine set_bounds_box
c
c     set domain and element bounds for a box geometry. Notice that
c     this ONLY works with non curved elements.
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

c     common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
c     integer gfirst, inoassignd, resetFindpts, pload


      if(istep.eq.0.or.istep.eq.1)then
c     if((istep.eq.0) .or. (istep.eq.1).or.(resetFindpts.eq.1)) then 
        call domain_size(xdrange(1,1),xdrange(2,1),xdrange(1,2)
     $                  ,xdrange(2,2),xdrange(1,3),xdrange(2,3))
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine set_part_params
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*132 deathmessage

      if (icalld .lt. 0) then
         rdum   = ran2(-nrandseed*np-nid-1) ! initialize
         icalld = icalld + 1
      endif

c     setup items
      nlxyze = nx1*ny1*nz1*nelt
      call rzero(rpart,lr*llpart)
      call izero(ipart,li*llpart)
      call rzero(ptw,nlxyze*8)
      call rzero(ptdum,iptlen)
      call rzero(pttime,iptlen)

      rtmp = 0.0 ! dummy number, max value it can be

      ! do nothing, no spreading
      if (npro_method .eq. 0) then

      ! gaussian set by user input parameters
      elseif (npro_method .eq. 1) then

         rtmp = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))

      endif

      d2chk(1) = rtmp
      d2chk(2) = d2chk(1)
      d2chk(3) = d2chk(1)

      rsig     = dfilt/(2.*sqrt(2.*log(2.))) ! gaussian filter std. * DP

      mu_0   = abs(param(2))

      return
      end
c----------------------------------------------------------------------
      subroutine set_check_spl_params_ibm
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      character*132 deathmessage

      if (npro_method .eq. 1) then
         rtmp = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))
      endif
      
      rdeff_max=ibm_diam(1)/2.*1.1
      !rdeff_max = dp(2)/2.
      do i = 1,n
         if(ipart(jrole,i) .eq. 1) then
         if (rpart(jrpe,i) .gt. rdeff_max) rdeff_max=rpart(jrpe,i)
         endif
      enddo

      rdeff_max = glmax(rdeff_max,1) 
      if(two_way.gt.2)rdeff_max = rdeff_max*2.0    ! force and motion influence

      if(non_spherical.ne.0) rdeff_max = 3*rdeff_max

      rtmp_col  = rdeff_max        ! collisional zone of influence
      !print*, "rdeff_max, rtmp_col", rdeff_max,rtmp_col
c      d2chk_ibm(1)  = rtmp*rdeff_max
c      d2chk_ibm(2)  = rtmp*rdeff_max
c      d2chk_ibm(3)  = rtmp*rdeff_max

c      d2chk_ibm(1) = max(d2chk_ibm(1),rtmp_col)
c      d2chk_ibm(2) = d2chk_ibm(2)
c      d2chk_ibm(3) = rtmp_col

      d2chk_ibm(1)  = max(rtmp_col, rdeff_max)
      d2chk_ibm(2)  = d2chk_ibm(1)
      d2chk_ibm(3)  = d2chk_ibm(1)

      ! change back the d2chk for the use the create neighbor element
      d2chk(1)  = d2chk_ibm(1)  
      d2chk(2)  = d2chk_ibm(2)  
      d2chk(3)  = d2chk_ibm(3)  
      
      return
      end

c----------------------------------------------------------------------
      subroutine set_check_spl_params
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      character*132 deathmessage

      if (npro_method .eq. 1) then
         rtmp = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))
      endif
      
      rdeff_max = 1e-8 ! dp(2)/2.
      do i = 1,n
         if(ipart(jrole,i) .ne. 1) then 
         if (rpart(jrpe,i) .gt. rdeff_max) rdeff_max=rpart(jrpe,i)
         endif
      enddo
      rdeff_max = glmax(rdeff_max,1)
      rdeff_max = rdeff_max*2. ! to get diameter
      
      rtmp_col = rdeff_max*1.50

      d2chk(1)  = rtmp*rdeff_max
      d2chk(2)  = rtmp*rdeff_max
      d2chk(3)  = rtmp*rdeff_max

      d2chk(1) = max(d2chk(1),rtmp_col)
      d2chk(2) = d2chk(2)
      d2chk(3) = rtmp_col

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_set_dt(rdt_part)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real dt_dum,dt_col,cflp,cflt

      common /save_dt_part/ rdpe_max, rdpe_min

      integer icalld
      save    icalld
      data    icalld  /-1/
      
      if (llpart .eq. 1) goto 1234

      icalld = icalld + 1
      if (icalld .eq. 0) then
         rdpe_max = 0.
         rdpe_min = 100.
         do i=1,n
            if (rpart(jrpe,i).lt.rdpe_min) rdpe_min = rpart(jrpe,i)
            if (rpart(jrpe,i).gt.rdpe_max) rdpe_max = rpart(jrpe,i)
         enddo
         rdpe_min = glmin(rdpe_min,1)
         rdpe_max = glmax(rdpe_max,1)
         rdpe_min = 2.*rdpe_min ! to get diameter
         rdpe_max = 2.*rdpe_max ! to get diameter
      endif

c     ! particle cfl, particles cant move due to velocity
      dt_dum = rdt_part
      dt_part = 1000.
      cflp = 0.10
      rvmag_max = 0.
      do i=1,n
         rvmag  = sqrt(rpart(jv0,i)**2 + rpart(jv0+1,i)**2 
     >                 + rpart(jv0+2,i)**2)

         cflt = dt_dum*rvmag/rdpe_min
         if (cflt .lt. cflp) dt_part = dt_dum
         if (cflt .ge. cflp) dt_part = cflp*rdpe_min/rvmag ! make sure smallest small overlap

         if (rvmag .gt. rvmag_max) rvmag_max = rvmag
      enddo
      rvmag_max = glmax(rvmag_max,1)
      dt_part  = glmin(dt_part,1)

      if(IBM_Particle_shape.ge.1) then
         rdpe_min = ibm_diam(1) 
         rdpe_max = ibm_diam(1)
      endif

      ! resolving collisions
      rm1      = rho_p*pi/6.*rdpe_min**3 ! min
      rm2      = rho_p*pi/6.*rdpe_max**3 ! max
      rm12     = 1./(1./rm1 + 1./rm2)
      n_resolve= 10
      dt_col   = sqrt(rm12/ksp*(log(e_rest)**2+pi**2))/n_resolve ! collision time scale
      
      if (two_way .gt. 2) then
         if(nid.eq.1.and.istep.le.1) 
     $     print*,'dt_part, dt_col', dt_part, dt_col
         rdt_part = min(dt_part,dt_col)
      else
         rdt_part = dt_part ! don't set if no collisions!
      endif

 1234 continue

      return
      end
c----------------------------------------------------------------------
      subroutine set_part_pointers
      include 'SIZE'
      include 'LPM'

      nr   = lr     ! Mandatory for proper striding
      ni   = li     ! Mandatory
      nrgp = lrgp
      nigp = ligp
      nrf  = lrf
      nif  = lif
      n    = 0

      nigp_ibm = ligp_ibm
      nrgp_ibm = lrgp_ibm
      
c     ipart pointers ------------------------------------------------
      jrc   = 1 ! Pointer to findpts return code
      jpt   = 2 ! Pointer to findpts return processor id
      je0   = 3 ! Pointer to findpts return element id
      jps   = 4 ! Pointer to proc id for data swap
      jpid1 = 5 ! initial proc number
      jpid2 = 6 ! initial local particle id, this is overall_id include 3 types of points
      jpid3 = 7 ! initial time step introduced
      jicx  = 8 ! initial time step introduced
      jicy  = 9 ! initial time step introduced
      jicz  = 10 ! initial time step introduced
      jai   = 11 ! Pointer to auxiliary integers

c for ibm particles     
      jrole  = 12   ! IBM: particle role: queen, worker, or regular point
      jqueen = 13   ! IBM: the queen numbering in local 
      jworker1 = 14 ! IBM: first worker id
      jnlm    = 15  ! IBM: number of workers
      
      nai = ni - (jnlm-1)  ! Number of auxiliary integers
      if (nai.le.0) call exitti('Error in nai:$',ni)

c     rpart pointers ------------------------------------------------
      jr  = 1         ! Pointer to findpts return rst variables
      jd  = jr + 3    ! Pointer to findpts return distance
      jx  = jd + 1    ! Pointer to findpts input x value
      jy  = jx + 1    ! Pointer to findpts input y value
      jz  = jy + 1    ! Pointer to findpts input z value
      jv0 = jz + 1    ! particle velocity at this timestep
      ju0 = jv0 + 3   ! fluid velocity at this time step
      jf0 = ju0 + 3   ! particle total force at this timestep
      jq0 = jf0 + 3   ! temperature forcing
      jg0 = jq0 + 1   ! work done by forces
      jquu= jg0 + 1   ! undisturbed unsteady temp. forcing
      jqqs= jquu+ 1   ! quasi-steady temp. forcing

c     forcing
      ii  = jqqs + 1
c     if (part_force(1).ne.0) then ! user specified force
         jfusr = ii
         ii    = ii + 3
c     endif
c     if (part_force(2).ne.0) then ! quasi-steady force
         jfqs  = ii
         ii    = ii + 3
c     endif
c     if (part_force(3).ne.0) then ! undisturbed force
         jfun  = ii
         ii    = ii + 3
c     endif
c     if (part_force(4).ne.0) then ! inviscid unsteady force
         jfiu  = ii
         ii    = ii + 3
c     endif

      jfcol  = ii  ! collisional force

c     other parameters (some may not be used; all at part. location)
      jtaup   = jfcol   + 3 ! particle time scale
      jcd     = jtaup   + 1 ! drag coeff
      jdrhodt = jcd     + 3 ! density material time derivative
      jre     = jdrhodt + 1 ! Relative Reynolds number
      jDuDt   = jre     + 1 ! fluid velocity time derivative
      jtemp   = jDuDt   + 3 ! part. temperature 
      jtempf  = jtemp   + 1 ! fluid temperature at part. loc.
      jrho    = jtempf  + 1 ! fluid denisty 
      jrhop   = jrho    + 1 ! particle material density
      ja      = jrhop   + 1 ! fluid mach number
      jvol    = ja      + 1 ! particle volume 
      jvol1   = jvol    + 1 ! particle volume fraction at part. loc.
      jdp     = jvol1   + 1 ! particle diameter
      jrpe    = jdp     + 1 ! particle effective diameter spl
      jgam    = jrpe    + 1 ! spread to grid correction
      jspl    = jgam    + 1 ! super particle loading
      jcmiu   = jspl    + 1 ! added mass coefficient

c     bdf/ext integration
      jx1 = jcmiu+1 ! Pointer to xyz at t^{n-1}
      jx2 = jx1 +3 ! Pointer to xyz at t^{n-1}
      jx3 = jx2 +3 ! Pointer to xyz at t^{n-1}

      jv1 = jx3+ 3 ! Pointer to particle velocity at t^{n-1}
      jv2 = jv1+ 3 ! Pointer to particle velocity at t^{n-2}
      jv3 = jv2+ 3 ! Pointer to particle velocity at t^{n-3}

      ju1 = jv3+ 3 ! Pointer to fluid velocity at t^{n-1}
      ju2 = ju1+ 3 ! Pointer to fluid velocity at t^{n-2}
      ju3 = ju2+ 3 ! Pointer to fluid velocity at t^{n-3}

      jangle0  = ju3 + 3     ! angular position in 3D 
      jangvel0 = jangle0 + 3 ! angular velocity in 3D

      jangvel1 = jangvel0 + 3 ! rotating velocity at t^{n-1}!
      jangvel2 = jangvel1 + 3 ! rotating velocity at t^{n-2}!
      jangvel3 = jangvel2 + 3 ! rotating velocity at t^{n-3}!

      jtorque0 = jangvel3 + 3 ! torque
      jtq_col  = jtorque0 + 3 ! collisional torque

      jangle1 = jtq_col + 3   ! rotating velocity at t^{n-1}!
      jangle2 = jangle1 + 3 ! rotating velocity at t^{n-2}!
      jangle3 = jangle2 + 3 ! rotating velocity at t^{n-3}!

      jar      = jangle3  + 3   ! Pointer to auxiliary reals
      
      nar = nr - (jar-1)  ! Number of auxiliary reals
      if (nar.le.0) call exitti('Error in nar:$',nr)

c     ghost particle integer pointers -------------------------------
      jgpiic  = 1 ! if gp used in collisions or not
      jgpps   = 2 ! Pointer to proc id for data swap
      jgppt   = 3 ! findpts return processor id
      jgpes   = 4 ! Destination element to be sent to
      jgpicx  = 5 !ii
      jgpicy  = 6 !jj
      jgpicz  = 7 !kk

c     ghost Queen/worker marker communication
      jgp_back     = 8     ! IBM: send back ibm ghost particle, force integration
      jgp_pid1     = 9     ! IBM: ghost particle identifier
      jgp_pid2     = 10    !
      jgp_pid3     = 11    !
      jgp_role     = 12    !
      jgp_queen    = 13    !
      jgp_worker1  = 14    !
      jgp_nlm      = 15    !
      jgp_ax       = 16    

      naigp = ligp_ibm - (jgp_ax-1)     ! Number of auxiliary integers
      if (naigp.le.0) call exitti('Error in naigp:$',ligp_ibm)
      
c     ghost particle real pointers ----------------------------------
      jgpx    = 1 ! ghost particle xloc
      jgpy    = 2 ! ghost particle yloc
      jgpz    = 3 ! ghost particle zloc
      jgpfh   = 4 ! ghost particle hydrodynamic xforce (i+1 > y, i+2 > z)
      jgpvol  = jgpfh+3  ! ghost particle volume
      jgprpe  = jgpvol+1 ! ghost particle effective radius
      jgpspl  = jgprpe+1 ! spreading correction (if used)
      jgpg0   = jgpspl+1 ! 10
      jgpq0   = jgpg0 +1 ! 11
      jgpv0   = jgpq0 +1 ! velocity (3 components) 12-14
c     above are for all ghost particle/marker

      jgp_jx1 = jgpv0 +3        ! previous location for update of worker marker,15-17

c     rotation motion for Ghost Queen marker only
      jgp_angle0  = jgp_jx1  + 3  ! anglular in spherical coordinate,18-20
      jgp_angvel0 = jgp_angle0  + 3  ! angular velocity, 21-23
      jgp_torque0 = jgp_angvel0 + 3  ! torque, 24-26

      return
      end
c----------------------------------------------------------------------
c     particle force routines
c----------------------------------------------------------------------
      subroutine lpm_usr_particles_solver
c
c     call routines in ordered way - main solver structure
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage

      logical ifinject
      integer icalld
      save    icalld
      data    icalld  /-1/

      if (icalld .eq. -1) then
         pttime(1) = 0.
      else
         pttime(1) = pttime(1) + dnekclock() - ptdum(1)
      endif

      icalld = icalld + 1

c     should we inject particles at this time step?
      ifinject = .false.
      if (inject_rate .gt. 0) then
      if ((mod(istep,inject_rate).eq.0)) then 
         ifinject = .true. 
      endif
      endif
      
      if (istep .gt. time_delay) then
      if (stage.eq.1) then
         ! Update coordinates if particle moves outside boundary
         ptdum(2) = dnekclock()
            call update_particle_location  
         pttime(2) = pttime(2) + dnekclock() - ptdum(2)

         ! Inject particles if needed
         if (ifinject) call place_particles

         ! Update where particle is stored at
         ptdum(3) = dnekclock()
            call move_particles_inproc
         pttime(3) = pttime(3) + dnekclock() - ptdum(3)

         if(ibm_debug_bin.eq.1.and.n.ge.1) then
            iqt = 0
            ipt = 0
            do i = 1,n
               if(ipart(jrole,i).eq.1) iqt = iqt + 1 
               if(ipart(jrole,i).eq.2) ipt = ipt + 1 
            enddo
            write(6,'(I4,A,3I6)') nid," contains ",n,ipt,iqt 
         endif


         if (two_way.gt.1) then
            ! Create ghost/wall particles
            ptdum(4) = dnekclock()
            call create_extra_particles
            pttime(4) = pttime(4) + dnekclock() - ptdum(4)
            
            if(ibm_debug_bin.eq.1 .and. nfptsgp_ibm.ge.1) 
     $ write(6,'(I6,A,I6,A)')nid," will send",nfptsgp_ibm," Ghost Queen"
            
            
            ! Send ghost markers
            ptdum(6) = dnekclock()
            call send_ghost_particles
            pttime(6) = pttime(6) + dnekclock() - ptdum(6)

            ptdum(7) = dnekclock()
            call send_ghost_particles_ibm
            pttime(7) = pttime(7) + dnekclock() - ptdum(7)

            !
            ptdum(5) = dnekclock()
            if(ifibm.eq.1) then
               call sort_local_particles_collisions_IBM
            else
               call sort_local_particles_collisions
            endif
            pttime(5) = pttime(5) + dnekclock() - ptdum(5)
            
            ! Projection to Eulerian grid
            ptdum(9) = dnekclock()
            call spread_props_grid
            pttime(9) = pttime(9) + dnekclock() - ptdum(9)  
         endif
      endif

      if(ifIBM .ge. 1) then
         if(stage.ne.1) call update_queen_periodic ! periodic for queen to avoid search outside boundary
         If(stage.ne.1) call create_ghost_particles_full_ibm ! avoid duplicated ghost queen
         if(stage.ne.1) call sort_local_particles_collisions_IBM !
         if(stage.ne.1) call send_ghost_particles_ibm            ! for force collection

         if(ibm_debug_bin.eq.1.and.nfptsgp_ibm.ge.1) 
     $   write(6,'(I6,A,I6,A)') nid," receives "
     $        ,nfptsgp_ibm," Ghost Queen"


         ! update worker location ( follow send ghost queen markers )
         if(ipart_moving.ge.1) then
            ptdum(8) = dnekclock()
            call IBM_local_worker_update
            pttime(8) = pttime(8) + dnekclock() - ptdum(8)
         endif
      endif

      ! Interpolate Eulerian properties to particle location
      ptdum(10) = dnekclock()
         call interp_props_part_location
      pttime(10) = pttime(10) + dnekclock() - ptdum(10)

      ! Evaluate particle force models
      ptdum(11) = dnekclock()
         call usr_particles_forcing
      pttime(11) = pttime(11) + dnekclock() - ptdum(11)

      ! IBM force integration
      if(ifIBM .ge. 1) then
!         call output_queens
!         call output_workers
         ptdum(12) = dnekclock()
         call IBM_part_force_collect
         pttime(12) = pttime(12) + dnekclock() - ptdum(12)

         ptdum(13) = dnekclock()
         call send_back_ghost_particles_ibm
         pttime(13) = pttime(13) + dnekclock() - ptdum(13)

         ptdum(14) = dnekclock()
         call IBM_part_forcing
         pttime(14) = pttime(14) + dnekclock() - ptdum(14)

         if(ioutput_queen.eq.1)call output_queens

         ptdum(15) = dnekclock()
         if(ipart_moving.ge.1) call IBM_part_integrate
         pttime(15) = pttime(15) + dnekclock() - ptdum(15)

      else
         ! Integrate in time
         ptdum(16) = dnekclock()
         if (abs(time_integ) .eq. 1) call rk3_integrate
         if (abs(time_integ) .eq. 2) call bdf_integrate
         pttime(16) = pttime(16) + dnekclock() - ptdum(16)
      endif

      ! Update forces
      ptdum(17) = dnekclock()
      call compute_forcing_post_part
      pttime(17) = pttime(17) + dnekclock() - ptdum(17)

      endif ! time_delay

      ptdum(1) = dnekclock()

      return
      end
c----------------------------------------------------------------------
      subroutine spread_props_grid
c
c     spread particle properties at fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'LPM'

      real               phig(lx1,ly1,lz1,lelt)
      common /otherpvar/ phig

      real    xx,yy,zz,vol,pfx,pfy,pfz,pmass,pmassf,vcell,multfc,multfci
     >       ,qgqf,rvx,rvy,rvz,rcountv(8,nelt),rx2(3)
     >       ,rxyzp(n_walls*2,3)
      integer e

      nlxyze = lx1*ly1*lz1*lelt
      nxyze  = nx1*ny1*nz1*nelt
      call rzero(ptw,nlxyze*8)
      
      ! do nothing
      if (npro_method .eq. 0) then

      ! gaussian spreading
      elseif (abs(npro_method) .eq. 1) then

c       write(*,*) "what is n in spread_prop_grid", n
      ! real particle projection
       do ip=1,n
         rsigp   = rsig*rpart(jrpe,ip)*2.
         multfci = 1./(sqrt(2.*pi)**2 * rsigp**2) ! exponential
         if (if3d) multfci = multfci**(1.5d+0)
         ralph   = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))*
     >           rpart(jrpe,ip)*2.
         ralph2   = ralph**2
         rbexpi   = 1./(-2.*rsigp**2)
         multfc = multfci
         if (.not. if3d) multfc = multfc/(2.*rpart(jrpe,ip))


         pfx = -rpart(jf0,ip)*multfc   ! * rpart(jvol,ip)   !!!!!!!!!!!!!!!! ibm spreading
         pfy = -rpart(jf0+1,ip)*multfc ! * rpart(jvol,ip) !!!!!!!!!!!!!!!!
         pfz = -rpart(jf0+2,ip)*multfc ! * rpart(jvol,ip) !!!!!!!!!!!!!!!!
         vol = rpart(jvol,ip)*multfc

c         if (rpart(jf0+1,ip) .ge. 1d-6 )
c     >   write(*,*) "in spread jf0+1=", rpart(jf0+1,ip)
c         if (abs(pfy) .ge. 1d-6 )
c     >    write(*,*) "in spread pfy=", pfy
c         if (abs(pfy) .ge. 1d-6 )
c     >    write(*,*) multfc,rbexpi,rsigp,ralph2,rpart(jrpe,ip)
         
         qgqf= -(rpart(jg0,ip) + rpart(jq0,ip))*multfc
         rvx = rpart(jv0  ,ip)*vol
         rvy = rpart(jv0+1,ip)*vol
         rvz = rpart(jv0+2,ip)*vol

         ii    = floor((rpart(jx,ip)-xdrange(1,1))/rdxgp) 
         jj    = floor((rpart(jy,ip)-xdrange(1,2))/rdygp) 
         kk    = floor((rpart(jz,ip)-xdrange(1,3))/rdzgp) 

         ! adding wall effects
         rx2(1) = rpart(jx,ip)
         rx2(2) = rpart(jy,ip)
         rx2(3) = rpart(jz,ip)
         call extra_wall_particle_exp(rxyzp,rx2,ic)

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
! ?????????????????????????????? what is these conditions?
      if (     mod_gp_grid(i,j,k,ie,1).ge.ii-1
     >   .and. mod_gp_grid(i,j,k,ie,1).le.ii+1) then
      if (     mod_gp_grid(i,j,k,ie,2).ge.jj-1
     >   .and. mod_gp_grid(i,j,k,ie,2).le.jj+1) then
      if (     mod_gp_grid(i,j,k,ie,3).ge.kk-1
     >   .and. mod_gp_grid(i,j,k,ie,3).le.kk+1) then

         rdist2  = (xm1(i,j,k,ie) - rpart(jx,ip))**2 +
     >           (ym1(i,j,k,ie) - rpart(jy,ip))**2 
         if (if3d) rdist2 = rdist2 + (zm1(i,j,k,ie) - rpart(jz,ip))**2 


         rexp = exp(rdist2*rbexpi)
c         if(abs(rexp) .ge. 1e-6 .and. pfy .ge.1d-6)
c     >   write(*,*) "pfy=",pfy,"rexp",rexp 

         ! add wall effects
         do jjj=1,ic
            rx22 = (xm1(i,j,k,ie) - rxyzp(jjj,1))**2
            ry22 = (ym1(i,j,k,ie) - rxyzp(jjj,2))**2
            rtmp2 = rx22 + ry22
            if (if3d) then
               rz22 = (zm1(i,j,k,ie) - rxyzp(jjj,3))**2
               rtmp2 = rtmp2 + rz22
            endif
            rexp = rexp + exp(rtmp2*rbexpi)
         enddo

         ptw(i,j,k,ie,1) = ptw(i,j,k,ie,1) + pfx*rexp
         ptw(i,j,k,ie,2) = ptw(i,j,k,ie,2) + pfy*rexp
         ptw(i,j,k,ie,3) = ptw(i,j,k,ie,3) + pfz*rexp
         ptw(i,j,k,ie,4) = ptw(i,j,k,ie,4) + vol*rexp
         ptw(i,j,k,ie,5) = ptw(i,j,k,ie,5) + qgqf*rexp
         ptw(i,j,k,ie,6) = ptw(i,j,k,ie,6) + rvx*rexp
         ptw(i,j,k,ie,7) = ptw(i,j,k,ie,7) + rvy*rexp
         ptw(i,j,k,ie,8) = ptw(i,j,k,ie,8) + rvz*rexp
      
      endif
      endif
      endif
      
c       if ( abs(ptw(i,j,k,ie,2)) .ge. 1e-6)
c     $ write(*,*) "in spread ptw = " , ptw(i,j,k,ie,2)

      enddo
      enddo
      enddo
      enddo
      enddo

      ! ghost particle projection
      do ip=1,nfptsgp

         rsigp   = rsig*rptsgp(jgprpe,ip)*2.
         multfci = 1./(sqrt(2.*pi)**2 * rsigp**2) ! exponential
         if (if3d) multfci = multfci**(1.5d+0)
         ralph   = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))*
     >           rptsgp(jgprpe,ip)*2.
         ralph2   = ralph**2
         rbexpi   = 1./(-2.*rsigp**2)
         multfc = multfci
         if (.not. if3d) multfc = multfc/(2.*rptsgp(jgprpe,ip))

         pfx = -rptsgp(jgpfh,ip)*multfc
         pfy = -rptsgp(jgpfh+1,ip)*multfc
         pfz = -rptsgp(jgpfh+2,ip)*multfc
         vol = rptsgp(jgpvol,ip)*multfc
         qgqf= -(rptsgp(jgpg0,ip) + rptsgp(jgpq0,ip))*multfc
         rvx = rptsgp(jgpv0  ,ip)*vol
         rvy = rptsgp(jgpv0+1,ip)*vol
         rvz = rptsgp(jgpv0+2,ip)*vol

         ii    = floor((rptsgp(jgpx,ip)-xdrange(1,1))/rdxgp) 
         jj    = floor((rptsgp(jgpy,ip)-xdrange(1,2))/rdygp) 
         kk    = floor((rptsgp(jgpz,ip)-xdrange(1,3))/rdzgp) 

         ! adding wall effects
         rx2(1) = rptsgp(jgpx,ip)
         rx2(2) = rptsgp(jgpy,ip)
         rx2(3) = rptsgp(jgpz,ip)
         call extra_wall_particle_exp(rxyzp,rx2,ic)

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1

      if (     mod_gp_grid(i,j,k,ie,1).ge.ii-1
     >   .and. mod_gp_grid(i,j,k,ie,1).le.ii+1) then
      if (     mod_gp_grid(i,j,k,ie,2).ge.jj-1
     >   .and. mod_gp_grid(i,j,k,ie,2).le.jj+1) then
      if (     mod_gp_grid(i,j,k,ie,3).ge.kk-1
     >   .and. mod_gp_grid(i,j,k,ie,3).le.kk+1) then

         rdist2  = (xm1(i,j,k,ie) - rptsgp(jgpx,ip))**2 +
     >           (ym1(i,j,k,ie) - rptsgp(jgpy,ip))**2 
         if (if3d) rdist2 = rdist2 + (zm1(i,j,k,ie)-rptsgp(jgpz,ip))**2 

         rexp = exp(rdist2*rbexpi)

         ! add wall effects
         do jjj=1,ic
            rx22 = (xm1(i,j,k,ie) - rxyzp(jjj,1))**2
            ry22 = (ym1(i,j,k,ie) - rxyzp(jjj,2))**2
            rtmp2 = rx22 + ry22
            if (if3d) then
               rz22 = (zm1(i,j,k,ie) - rxyzp(jjj,3))**2
               rtmp2 = rtmp2 + rz22
            endif
            rexp = rexp + exp(rtmp2*rbexpi)
         enddo

         ptw(i,j,k,ie,1) = ptw(i,j,k,ie,1) + pfx*rexp
         ptw(i,j,k,ie,2) = ptw(i,j,k,ie,2) + pfy*rexp
         ptw(i,j,k,ie,3) = ptw(i,j,k,ie,3) + pfz*rexp
         ptw(i,j,k,ie,4) = ptw(i,j,k,ie,4) + vol*rexp
         ptw(i,j,k,ie,5) = ptw(i,j,k,ie,5) + qgqf*rexp
         ptw(i,j,k,ie,6) = ptw(i,j,k,ie,6) + rvx*rexp
         ptw(i,j,k,ie,7) = ptw(i,j,k,ie,7) + rvy*rexp
         ptw(i,j,k,ie,8) = ptw(i,j,k,ie,8) + rvz*rexp

      endif
      endif
      endif

      enddo
      enddo
      enddo
      enddo
      enddo

      endif

#ifdef CMTNEK
      ! filtering makes velocity field smoother for p*grad(phig) in CMT
      wght = 1.0
      ncut = 1
      call filter_s0(ptw(1,1,1,1,4),wght,ncut,'phip') 
#endif

      rvfmax = 0.7405
      rvfmin = 0.0
      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         if (ptw(i,j,k,ie,4) .gt. rvfmax) ptw(i,j,k,ie,4) = rvfmax
         if (ptw(i,j,k,ie,4) .lt. rvfmin) ptw(i,j,k,ie,4) = rvfmin
         phig(i,j,k,ie) = 1. - ptw(i,j,k,ie,4)
      enddo
      enddo
      enddo
      enddo

      ! periodic settling particle conditions
      if(ifIBM.eq.1) then

         if(Periodic_Set.eq. 1) then ! Periodic Settling
         ! Method 1 remove the mean velocity of the whole flow domain

            ntm1  = nx1*ny1*nz1*nelt
            vxbar = glsc2(vx, bm1, ntm1)/ voltm1 ! spatial average velocity
            vybar = glsc2(vy, bm1, ntm1)/ voltm1
            vzbar = glsc2(vz, bm1, ntm1)/ voltm1

            do ie=1,nelt
               do k=1,nz1
                  do j=1,ny1
                     do i=1,nx1
                        vx(i,j,k,ie) = vx(i,j,k,ie) - vxbar
                        vy(i,j,k,ie) = vy(i,j,k,ie) - vybar
                        vz(i,j,k,ie) = vz(i,j,k,ie) - vzbar
                     enddo
                  enddo
               enddo
            enddo


            
         elseif(Periodic_Set.eq.2) then
            ! for periodic BC, subtract the spatial average force term for the compatible
            ! Method 2 remove the mean body force of whole domain

            rfx_sum_loc = 0.0
            rfy_sum_loc = 0.0
            rfz_sum_loc = 0.0

            do ie=1,nelt
               do k=1,nz1
                  do j=1,ny1
                     do i=1,nx1
                        rfx_sum_loc = rfx_sum_loc + ptw(i,j,k,ie,1)
                        rfy_sum_loc = rfy_sum_loc + ptw(i,j,k,ie,2)
                        rfz_sum_loc = rfz_sum_loc + ptw(i,j,k,ie,3)
                     enddo
                  enddo
               enddo
            enddo
            ntot1=lx1*ly1*lz1*nelt
            rfx_avg = glsum(rfx_sum_loc,1) / iglsum(ntot1,1)
            rfy_avg = glsum(rfy_sum_loc,1) / iglsum(ntot1,1)
            rfz_avg = glsum(rfz_sum_loc,1) / iglsum(ntot1,1)
           ! subtract the spatial average of the body force
            do ie=1,nelt
               do k=1,nz1
                  do j=1,ny1
                     do i=1,nx1
                        ptw(i,j,k,ie,1) = ptw(i,j,k,ie,1) - rfx_avg
                        ptw(i,j,k,ie,2) = ptw(i,j,k,ie,2) - rfy_avg
                        ptw(i,j,k,ie,3) = ptw(i,j,k,ie,3) - rfz_avg
                     enddo
                  enddo
               enddo
            enddo           
         endif

      endif
      
      return
      end
c----------------------------------------------------------------------
      subroutine rk3_integrate
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage
      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      common /PARTRK3/ kv_stage_p
      real   kv_stage_p(llpart,13) 

      integer fdim
      real    pmass
      

c      ipart_stationary = 1 ! stationary particles for debug only
      
      jx0 = jx

c     rk3 stage one items ---------------------------------------------
      if (stage.eq.1) then
c        used for time derivative of v in iu force
         call get_bdf_ext_coefs(beta,alpha,times)

c        move data to previous positions
         do j=0,ndim-1
         do i=1,n
            rpart(ju3+j,i)=rpart(ju2+j,i)
            rpart(ju2+j,i)=rpart(ju1+j,i)
            rpart(ju1+j,i)=rpart(ju0+j,i)
            rpart(jv3+j,i)=rpart(jv2+j,i)
            rpart(jv2+j,i)=rpart(jv1+j,i)
            rpart(jv1+j,i)=rpart(jv0+j,i)
            rpart(jx3+j,i)=rpart(jx2+j,i)
            rpart(jx2+j,i)=rpart(jx1+j,i)
            rpart(jx1+j,i)=rpart(jx0+j,i)
         enddo
         enddo

         do i=1,n
            kv_stage_p(i,1) = rpart(jx0  ,i)
            kv_stage_p(i,2) = rpart(jx0+1,i)
            kv_stage_p(i,3) = rpart(jx0+2,i)
            kv_stage_p(i,4) = rpart(jv0  ,i)
            kv_stage_p(i,5) = rpart(jv0+1,i)
            kv_stage_p(i,6) = rpart(jv0+2,i)
            kv_stage_p(i,7) = rpart(jtemp,i)
         enddo
      endif

      if(ipart_moving .eq. 1) then 

c     all rk3 stages items --------------------------------------------
      do i=1,n
         if(ipart(jrole,i).ne.0 ) cycle 
         rpart(jv0  ,i) = tcoef(1,stage)*kv_stage_p(i,4) +
     >                    tcoef(2,stage)*rpart(jv0  ,i)  +
     >                    tcoef(3,stage)*rpart(jf0  ,i)
         rpart(jv0+1,i) = tcoef(1,stage)*kv_stage_p(i,5) +
     >                    tcoef(2,stage)*rpart(jv0+1,i)  +
     >                    tcoef(3,stage)*rpart(jf0+1,i)
         rpart(jv0+2,i) = tcoef(1,stage)*kv_stage_p(i,6) +
     >                    tcoef(2,stage)*rpart(jv0+2,i)  +
     >                    tcoef(3,stage)*rpart(jf0+2,i)
         rpart(jx0  ,i) = tcoef(1,stage)*kv_stage_p(i,1) +
     >                    tcoef(2,stage)*rpart(jx0  ,i)  +
     >                    tcoef(3,stage)*rpart(jv0  ,i)
         rpart(jx0+1,i) = tcoef(1,stage)*kv_stage_p(i,2) +
     >                    tcoef(2,stage)*rpart(jx0+1,i)  +
     >                    tcoef(3,stage)*rpart(jv0+1,i)
         if (if3d)
     >   rpart(jx0+2,i) = tcoef(1,stage)*kv_stage_p(i,3) +
     >                    tcoef(2,stage)*rpart(jx0+2,i)  +
     >                    tcoef(3,stage)*rpart(jv0+2,i)
         rpart(jtemp,i) = tcoef(1,stage)*kv_stage_p(i,7) +
     >                    tcoef(2,stage)*rpart(jtemp,i)  +
     >                    tcoef(3,stage)*rpart(jq0  ,i)
      enddo

      elseif(ipart_moving.eq.2) then
         call ibm_part_force_moving
      else
         if(istep.lt.2 .and.nid.eq.0)
     $        print*, "Stationary Particle"
      endif ! stationary particles for debug only

      return
      end

c-----------------------------------------------------------------------
      subroutine ibm_part_force_moving
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      jx0 = jx
      
      do i=1,n
         rpart(jv0  ,i) = rvpx_force
         rpart(jv0+1,i) = rvpy_force
         if (if3d) rpart(jv0+2,i) = rvpz_force
         rpart(jx  ,i) = rpart(jx  ,i) + rvpx_force*dt 
         rpart(jx+1,i) = rpart(jx+1,i) + rvpy_force*dt
         if (if3d)
     >   rpart(jx+2,i) = rpart(jx+2,i) + rvpz_force*dt
         !rpart(jtemp,i) = 
         if(ipart(jrole,i).eq.1) print*, "new_position",rpart(jx+1,i)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine usr_particles_forcing
c
c     calculate the rhs of particle equation
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real vel_diff,pmass,pmassf

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,13)

      integer local_marker_role
      common /part_force_calc/local_marker_role
      if ((part_force(2) .ne. 0) .or. (part_force(3) .ne. 0))
     >   call calc_substantial_derivative

      ! initialize collisional force and torque
      do i=1,n
         rpart(jfcol  ,i) = 0.
         rpart(jfcol+1,i) = 0.
         rpart(jfcol+2,i) = 0.
      enddo
      do i=1,n
         rpart(jtq_col  ,i) = 0.
         rpart(jtq_col+1,i) = 0.
         rpart(jtq_col+2,i) = 0.
      enddo

      do i=1,n
c        setup values ------------------------------------------------
         pmass = rpart(jvol,i)*rpart(jrhop,i)
         pmassf= rpart(jvol,i)*rpart(jrho,i)
         if(part_force(3).ne.0) pmass = pmass + rpart(jcmiu,i)*pmassf ! am

         call get_part_use_block(i,0)

         local_marker_role = ipart(jrole,i) !!! for IBM Queen/Worker Marker
c        momentum rhs ------------------------------------------------
         lpmforce(1) = 0
         lpmforce(2) = 0
         lpmforce(3) = 0
         lpmforcec(1) = 0 ! coupled to fluid
         lpmforcec(2) = 0 ! coupled to fluid
         lpmforcec(3) = 0 ! coupled to fluid
         call lpm_usr_f  ! user/body force
            rpart(jfusr+0,i) = lpmforce(1) + lpmforcec(1)
            rpart(jfusr+1,i) = lpmforce(2) + lpmforcec(2)
            rpart(jfusr+2,i) = lpmforce(3) + lpmforcec(3)

         if(local_marker_role.eq.0) then

            if (time_integ .gt. 0) then
            call lpm_f_qs
               rpart(jfqs+0,i) = lpmforce(1)
               rpart(jfqs+1,i) = lpmforce(2)
               rpart(jfqs+2,i) = lpmforce(3)
            call lpm_f_un
               rpart(jfun+0,i) = lpmforce(1)
               rpart(jfun+1,i) = lpmforce(2)
               rpart(jfun+2,i) = lpmforce(3)
            call lpm_f_iu
               rpart(jfiu+0,i) = lpmforce(1)
               rpart(jfiu+1,i) = lpmforce(2)
               rpart(jfiu+2,i) = lpmforce(3)
            endif
            call usr_particles_f_col(i) ! colision force 

         elseif(local_marker_role.eq.2) then
            call ibm_f_marker(i)
            rpart(jfqs+0,i) = lpmforce(1)
            rpart(jfqs+1,i) = lpmforce(2)
            rpart(jfqs+2,i) = lpmforce(3)
         elseif(local_marker_role.eq.1) then
            call usr_particles_f_col_IBM(i) ! only for IBM Queen collsion 
         else
            print*, "undefined marker role",nid,i
         endif
      
         rfx = rpart(jfusr+0,i) + 
     >         rpart(jfqs +0,i) +
     >         rpart(jfun +0,i) +
     >         rpart(jfiu +0,i) +
     >         rpart(jfcol+0,i)
         rfy = rpart(jfusr+1,i) + 
     >         rpart(jfqs +1,i) +
     >         rpart(jfun +1,i) +
     >         rpart(jfiu +1,i) +
     >         rpart(jfcol+1,i)
         if (if3d)
     >   rfz = rpart(jfusr+2,i) + 
     >         rpart(jfqs +2,i) +
     >         rpart(jfun +2,i) +
     >         rpart(jfiu +2,i) +
     >         rpart(jfcol+2,i)

         if(ifIBM.eq.0) then
         ! mass weighted total force
            rpart(jf0+0,i) = rfx/pmass
            rpart(jf0+1,i) = rfy/pmass
            rpart(jf0+2,i) = rfz/pmass
         else
            rpart(jf0+0,i) = rfx
            rpart(jf0+1,i) = rfy
            rpart(jf0+2,i) = rfz
         endif
         

c        energy rhs --------------------------------------------------
         lpmheat = 0
         if (time_integ .gt. 0) then
            call lpm_q_uu
               rpart(jquu,i) = lpmheat
            call lpm_q_qs
               rpart(jqqs,i) = lpmheat

            rqq = rpart(jquu,i) +
     >            rpart(jqqs,i)

            pmass = rpart(jvol,i)*rpart(jrhop,i)
            rpart(jq0,i) = rqq/(pmass*cp_p)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_part_use_block(i,ipre)
c
c     set common block for particle-user use
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
#ifdef CMTNEK
      include 'CMTDATA'
#endif

      real MixtPerf_C_GRT_part
      external MixtPerf_C_GRT_part

      if (ipre .eq. 0) then
      lpmx_p(1)    = rpart(jx,i)
      lpmx_p(2)    = rpart(jy,i)
      lpmx_p(3)    = rpart(jz,i)
      lpmvol_p     = rpart(jvol,i)
      lpmvolfrac_p = rpart(jvol1,i)
      lpmvolfrac_f = 1.0 - lpmvolfrac_p
      lpmtau_p     = rpart(jtaup,i)
      lpmv_p(1)    = rpart(jv0+0,i)
      lpmv_p(2)    = rpart(jv0+1,i)
      lpmv_p(3)    = rpart(jv0+2,i)
      lpmv_f(1)    = rpart(ju0+0,i)
      lpmv_f(2)    = rpart(ju0+1,i)
      lpmv_f(3)    = rpart(ju0+2,i)
      lpmvisc_f    = mu_0
      lpmdiam_p    = rpart(jdp,i)
      lpmdens_p    = rpart(jrhop,i)
      lpmdens_f    = rpart(jrho,i)
      lpmtemp_p    = rpart(jtemp,i)
      lpmtemp_f    = rpart(jtempf,i)
      lpmDuDt(1)   = rpart(jDuDt+0,i)
      lpmDuDt(2)   = rpart(jDuDt+1,i)
      lpmDuDt(3)   = rpart(jDuDt+2,i)
      lpmkappa_f   = abs(param(8))
      lpmvdiff_pf  = sqrt((rpart(ju0  ,i)-rpart(jv0  ,i))**2+
     >                    (rpart(ju0+1,i)-rpart(jv0+1,i))**2+
     >                    (rpart(ju0+2,i)-rpart(jv0+2,i))**2)
      if (abs(lpmvdiff_pf) .lt. 1E-6) lpmvdiff_pf=1E-6

      lpmmach_p    = MixtPerf_C_GRT_part(gmaref,rgasref,
     >                                   rpart(jtempf,i),icmtp)
      if (icmtp .eq. 0) then
         lpmach_p  = 0.0
      elseif (icmtp .eq. 1) then
         lpmmach_p = lpmvdiff_pf/lpmmach_p
      endif
      lpmre_p      = rpart(jrho,i)*rpart(jdp,i)*lpmvdiff_pf/mu_0 ! Re
      lpmi         = i

      
      elseif (ipre .eq. 1) then

      rpart(jx,i)        =  lpmx_p(1)    
      rpart(jy,i)        =  lpmx_p(2)    
      rpart(jz,i)        =  lpmx_p(3)    
      rpart(jvol,i)      =  lpmvol_p     
      rpart(jvol1,i)     =  lpmvolfrac_p 
      lpmvolfrac_f       =  1.0 - lpmvolfrac_p
      rpart(jtaup,i)     =  lpmtau_p     
      rpart(jv0+0,i)     =  lpmv_p(1)    
      rpart(jv0+1,i)     =  lpmv_p(2)    
      rpart(jv0+2,i)     =  lpmv_p(3)    
      rpart(ju0+0,i)     =  lpmv_f(1)    
      rpart(ju0+1,i)     =  lpmv_f(2)    
      rpart(ju0+2,i)     =  lpmv_f(3)    
      rpart(jdp,i)       =  lpmdiam_p    
      rpart(jrhop,i)     =  lpmdens_p    
      rpart(jrho,i)      =  lpmdens_f    
      rpart(jtemp,i)     =  lpmtemp_p    
      rpart(jtempf,i)    =  lpmtemp_f    
      rpart(jDuDt+0,i)   =  lpmDuDt(1)   
      rpart(jDuDt+1,i)   =  lpmDuDt(2)   
      rpart(jDuDt+2,i)   =  lpmDuDt(3)   

      lpmvdiff_pf  = sqrt((rpart(ju0  ,i)-rpart(jv0  ,i))**2+
     >                    (rpart(ju0+1,i)-rpart(jv0+1,i))**2+
     >                    (rpart(ju0+2,i)-rpart(jv0+2,i))**2)
      if (abs(lpmvdiff_pf) .lt. 1E-6) lpmvdiff_pf=1E-6

      lpmmach_p    = MixtPerf_C_GRT_part(gmaref,rgasref,
     >                                   rpart(jtempf,i),icmtp)
      lpmmach_p    = lpmvdiff_pf/lpmmach_p
      lpmre_p      = rpart(jrho,i)*rpart(jdp,i)*lpmvdiff_pf/mu_0 ! Re

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_forcing_post_part
c
c     post calculate forces due to factoring of equations
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real uvel(0:2), vvel(0:2), pmass, pmassf,S_qs

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,13)

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      pi  = 4.0d+0*atan(1.0d+0)

      do i=1,n !!nwe ! change from n to nwe ?
         pmass = rpart(jvol,i)*rpart(jrhop,i)
         pmassf= rpart(jvol,i)*rpart(jrho,i)
         if (part_force(3).ne.0) pmass =pmass + rpart(jcmiu,i)*pmassf

         call get_part_use_block(i,0)

         ! get user coupled to fluid forcing
         lpmforcec(1) = 0 ! coupled to fluid
         lpmforcec(2) = 0 ! coupled to fluid
         lpmforcec(3) = 0 ! coupled to fluid
         call lpm_usr_f   ! user/body force

c        momentum forcing to fluid
         do j=0,ndim-1
            rdum = 0.

            rdvdt = rpart(jf0+j,i) ! note already divided by Mp + am
            ram_s = rdvdt*rpart(jcmiu,i)*pmassf
            rpart(jfiu+j,i) = rpart(jfiu+j,i) - ram_s

c           note that no coupled f_un in this formulation
            rdum = rdum + lpmforcec(j+1)
            rdum = rdum + rpart(jfqs+j,i)
            rdum = rdum + rpart(jfiu+j,i)

            rpart(jf0+j,i) = rdum ! now the force to couple with gas

            if(ifIBM.eq.1) rpart(jf0+j,i) = rpart(jfqs+j,i)
         enddo

c        energy forcing to fluid (quasi-steady)
         rpart(jg0,i) = rpart(jv0  ,i)*rpart(jfqs  ,i) + !force work
     >                  rpart(jv0+1,i)*rpart(jfqs+1,i) +
     >                  rpart(jv0+2,i)*rpart(jfqs+2,i)
         rpart(jg0,i) = rpart(jg0,i) + 
     >                  rpart(ju0  ,i)*rpart(jfiu  ,i) + !iu
     >                  rpart(ju0+1,i)*rpart(jfiu+1,i) +
     >                  rpart(ju0+2,i)*rpart(jfiu+2,i)

         rdum = 0.
         rdum = rdum + rpart(jqqs,i)

         rpart(jq0,i) = rdum

      enddo

      return
      end

c---------------------------------------------------------------------------
      subroutine IBM_part_force_collect
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sum up forces from remote worker marker and save to the Ghost Queen marker 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      real dragf(3) , rdragf(3)
      real dragsum0,dragsum1,dragsum2

      real x_o, y_o, z_o
      real x_c, y_c, z_c
      real rl_x, rl_y, rl_z
      real rf_x, rf_y, rf_z
      real rtorque_x, rtorque_y,rtorque_z

      integer isearch_jpid1,isearch_jpid3,isearch_jqueen
      common/search_queen_gl/isearch_jpid1,isearch_jpid3,isearch_jqueen
      
      call rzero(dragf,3)
      dragsum0 = 0.0
      dragsum1 = 0.0
      dragsum2 = 0.0
      do i = 1,n
         do j=0,ndim-1
            if(ipart(jrole,i).ne.2) cycle ! worker
             dragf(j+1) = dragf(j+1) + rpart(jf0+j,i)
             if(j.eq.0) dragsum0 = dragsum0+rpart(jf0+j,i)
             if(j.eq.1) dragsum1 = dragsum1+rpart(jf0+j,i)
             if(j.eq.2) dragsum2 = dragsum2+rpart(jf0+j,i)
         enddo
      enddo
      dragsum_tot0 = glsum(dragsum0,1)
      dragsum_tot1 = glsum(dragsum1,1)
      dragsum_tot2 = glsum(dragsum2,1)
      ntot0 = 0
      do i =1,n
         if(ipart(jrole,i).ne.1) cycle 
         ntot0 = ntot0 + 1
      enddo
      ntot01 = iglsum(ntot0,1)
      if(ibm_debug_worker.eq.1 .and. nid.eq.0) 
     $  write(6,'(a,I6,6e12.4)') "HydroForce worker",
     $  ntot01,dragsum_tot0, dragsum_tot1, dragsum_tot2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

      ! initialize force
      do j=0,ndim-1
         do i_gp = 1,nfptsgp_ibm ! Ghost queen
            rptsgp_ibm(jgpfh       + j, i_gp) = 0.0
            rptsgp_ibm(jgp_torque0 + j, i_gp) = 0.0
         enddo
         do i_qt = 1,n           ! real Queen
            if(ipart(jrole,i_qt).ne.1) cycle
            rpart(jtorque0 + j, i_qt) =
     $      rpart(jtq_col  + j, i_qt)            
            ipart(jai, i_qt)      = 0
         enddo
      enddo

!      do i_qt = 1,n
!         if(ipart(jrole,i_qt).eq.1) then
!             print*,"queen orgin"
!     $       , (rpart(jf0+j,i_qt),j=0,2)
!     $       , (rpart(jtorque0+j,i_qt),j=0,2)
!          endif
!      enddo

      icollect_flag = 1
      ioutput_debug = 0
      
      do ipt = 1,n
         
         if(ipart(jrole,ipt) .ne. 2) cycle ! worker only
         ipt_pid1  = ipart(jpid1,ipt)
         ipt_pid3  = ipart(jpid3,ipt)
         ipt_queen = ipart(jqueen,ipt)
         icollect_flag = 0
         
         ! Local Queen
         do i_qt = 1,n               
            if(ipart(jrole,i_qt) .ne. 1) cycle ! real Queen 
            i_qt_pid1  = ipart(jpid1,i_qt) 
            i_qt_pid3  = ipart(jpid3,i_qt) 
            i_qt_queen = ipart(jqueen,i_qt) 
            if (i_qt_pid1 .eq. ipt_pid1 .and.
     $          i_qt_pid3 .eq. ipt_pid3 .and.
     $          i_qt_queen.eq. ipt_queen ) then               

!               ipart(jai, i_qt) = ipart(jai, i_qt) + 1 ! number of workers

               do j = 0,ndim-1
                  rpart(jf0+j, i_qt) = rpart(jf0+j, i_qt)
     $                               + rpart(jf0+j, ipt)
               enddo
               
               ! rotating torque
               if(ibm_rotation.eq.1) then
               !1 calculate torque on the worker marker 
                  x_o  = rpart(jx,     ipt) ! surface point 
                  y_o  = rpart(jx + 1, ipt)
                  z_o  = rpart(jx + 2, ipt)

                  x_c  = rpart(jx ,    i_qt) ! center
                  y_c  = rpart(jx + 1, i_qt) 
                  z_c  = rpart(jx + 2, i_qt)
                      
                  rl_x = x_o - x_c ! length vecotor r from center
                  rl_y = y_o - y_c
                  rl_z = z_o - z_c

                  rf_x = rpart(jf0,  ipt) ! force vector
                  rf_y = rpart(jf0+1,ipt)
                  rf_z = rpart(jf0+2,ipt)
c     
                  rtorque_x = rl_y * rf_z  - rl_z * rf_y ! calculate rotating velocity
                  rtorque_y = rl_z * rf_x  - rl_x * rf_z
                  rtorque_z = rl_x * rf_y  - rl_y * rf_x

                  rpart(jtorque0,   ipt) = rtorque_x
                  rpart(jtorque0+1, ipt) = rtorque_y
                  rpart(jtorque0+2, ipt) = rtorque_z
                  
                  !2 collect to the queen worker
                  do j = 0,ndim-1
                        rpart(jtorque0+j, i_qt) 
     $               =  rpart(jtorque0+j, i_qt)
     $               +  rpart(jtorque0+j, ipt)
                  enddo

                  if(ibm_debug_worker.eq.1)
     $              print*, "local worker force",nid,ipt
     $              ,  ipart(jpid1,i),ipart(jpid2,i)
     $              ,  (rpart(jf0+j,ipt),j=0,2)
     $              ,  (rpart(jtorque0+j,ipt),j=0,2)

               endif            ! rotate
               icollect_flag = 1
            endif  !match
         enddo

         ! Ghost Queen
         do i_gp = 1,nfptsgp_ibm 
            i_gp_pid1  = iptsgp_ibm(jgp_pid1,i_gp) ! (1) initial proc number  
            i_gp_pid2  = iptsgp_ibm(jgp_pid2,i_gp) ! initial local particle id 
            i_gp_pid3  = iptsgp_ibm(jgp_pid3,i_gp) ! (2) initial time step introduced 
            i_gp_role  = iptsgp_ibm(jgp_role,i_gp) ! 
            i_gp_queen = iptsgp_ibm(jgp_queen,i_gp)! (3)
            if ( i_gp_pid1 .eq. ipt_pid1 .and.
     $           i_gp_pid3 .eq. ipt_pid3 .and.
     $           i_gp_queen.eq. ipt_queen ) then ! match queen & worker

!               iptsgp_ibm(jgp_ax, i_gp) = iptsgp_ibm(jgp_ax, i_gp) + 1 
               do j = 0,ndim-1
                  rptsgp_ibm(jgpfh+j, i_gp) = rptsgp_ibm(jgpfh+j, i_gp)
     $                 + rpart(jf0+j,        ipt) 
               enddo
            
               ! rotating torque
               if(ibm_rotation.eq.1) then
               !1 calculate torque on the worker marker 
                  x_o  = rpart(jx,     ipt) ! surface point, local
                  y_o  = rpart(jx + 1, ipt)
                  z_o  = rpart(jx + 2, ipt)

                  x_c  = rptsgp_ibm(jgpx, i_gp) ! center
                  y_c  = rptsgp_ibm(jgpy, i_gp) 
                  z_c  = rptsgp_ibm(jgpz, i_gp)
                      
                  rl_x = x_o - x_c ! length vecotor r from center
                  rl_y = y_o - y_c
                  rl_z = z_o - z_c

                  rf_x = rpart(jf0,  ipt) ! force vector, local
                  rf_y = rpart(jf0+1,ipt)
                  rf_z = rpart(jf0+2,ipt)
c     
                  rtorque_x = rl_y * rf_z  - rl_z * rf_y ! calculate rotating velocity
                  rtorque_y = rl_z * rf_x  - rl_x * rf_z
                  rtorque_z = rl_x * rf_y  - rl_y * rf_x

                  rpart(jtorque0,   ipt) = rtorque_x
                  rpart(jtorque0+1, ipt) = rtorque_y
                  rpart(jtorque0+2, ipt) = rtorque_z
                !2 collect to the queen worker
                  do j = 0,ndim-1
                      rptsgp_ibm(jgp_torque0+j, i_gp)
     $              = rptsgp_ibm(jgp_torque0+j, i_gp)  
     $              + rpart     (jtorque0+j,     ipt)  ! add worker to Ghost Queen
                  enddo

                  if(ibm_debug_worker.eq.1) then
                  print*, "remote worker force",nid,ipt 
     $               ,  ipart(jpid1,i),ipart(jpid2,i)
     $               ,  (rpart(jf0+j,ipt) ,j=0,2) 
     $               ,  (rpart(jtorque0+j,ipt),j=0,2)
                  endif
                  
               endif            ! rotate
               icollect_flag = 1
            endif ! match
         enddo ! GP

         if(icollect_flag.eq.0) then
            ipart(jtemp,  ipt) = 10000 !flag 
            ipart(jtempf, ipt) = 10000 !flag 
            rxval = rpart(jx,ipt)
            ryval = rpart(jy,ipt)
            rzval = rpart(jz,ipt)

            iip    = floor((rxval-xdrange(1,1))/rdxgp_ibm)
            jjp    = floor((ryval-xdrange(1,2))/rdygp_ibm)
            kkp    = floor((rzval-xdrange(1,3))/rdzgp_ibm)

            print*,"Error in Force collection, ",
     $         nid, ipt
     $        ,ipart(jpid1, ipt),   ipart(jpid2, ipt)
     $        ,ipart(jpid3, ipt),   ipart(jqueen, ipt)
     $        ,ipart(jworker1,ipt),ipart(jnlm, ipt)
     $        , "xyz" ,(rpart(jx+j,ipt), j=0,2)
     $        ,'queen=',ipt_pid1,ipt_pid3,ipt_queen
     $        ," bin=",iip,jjp,kkp
     $        ,(rpart(jf0+j,ipt),j=0,2)                     
     $        ,(rpart(jv0+j,ipt),j=0,2)
            print*,"Error in Force collection, ipart ",
     $         (ipart(j,ipt),j=1,li)
            
c seach where is the queen?
            isearch_jpid1 = ipt_pid1
            isearch_jpid3 = ipt_pid3
            isearch_jqueen = ipt_queen
            call bcast(isearch_jpid1,  isize)
            call bcast(isearch_jpid3,  isize)
            call bcast(isearch_jqueen, isize)
            call seach_queen_global

            ioutput_debug = 1
            call bcast(ioutput_debug, isize) 
!            EXIT
         endif
         
      enddo ! ipt

!      if(ioutput_debug.eq.1) then
!         print*, "Diagnostic ouput", nid
!         call lpm_usr_particles_io
!         print*, "Diagnostic ouput end", nid
!         call exitt
!      endif

c     debug
      if(ibm_debug .eq. 1) then
         do i_qt = 1,n
            if(ipart(jrole,i_qt) .ne. 1) cycle ! local Queen
!            if(ibm_debug_worker .eq. 1)
!     $       write(6,'(A,2I6,3F12.6)')"Local Queen Force",nid
!     $         ,i_qt, (rpart(jf0+j,i_qt),j=0,2)
         enddo         
         if(nfptsgp_ibm.ge.1) then
            do i_gp = 1,nfptsgp_ibm
!               if(ibm_debug_worker .eq. 1)
!     $         write(6,'(A,2I6,3F12.6)')"Ghost Queen Force",nid
!     $         ,i_gp,(rptsgp_ibm(jgpfh+j,i_gp),j=0,2)
            enddo
         endif
      endif

      return
      end subroutine IBM_part_force_collect

c-----------------------------------------------------------------------
      subroutine calc_substantial_derivative
c 
c     calculate rhs of NS, which is just pressure gradient. used for
c     undisturbed and invisicid unsteady force
c
c     no forcing included...should it be?
c
c     rhs_fluidp(i,j,k,e,l)
c   
c        l = 1     >    dP/dx
c        l = 2     >    dP/dy
c        l = 3     >    dP/dz
c        l = 4     >    -P* div(phi_p v), v is eulerian particle vel
c        l = 5     >    P * d phi_g/dx
c        l = 6     >    P * d phi_g/dy
c        l = 7     >    P * d phi_g/dz
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'LPM'

      real               phig(lx1,ly1,lz1,lelt)
      common /otherpvar/ phig

      integer e
      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1),
     >         pm1(lx1,ly1,lz1,lelt,3),udum(lx1,ly1,lz1)

      nxyz=nx1*ny1*nz1
      nlxyze = lx1*ly1*lz1*lelt

      call rzero(rhs_fluidp,nlxyze*7)

      ! if pressure solved on different mesh, map to vel mesh
      if (lx2 .ne. lx1) then
         call mappr(pm1(1,1,1,1,1),pr,pm1(1,1,1,1,2),pm1(1,1,1,1,3))
      else
         call copy(pm1(1,1,1,1,1),pr(1,1,1,1),nlxyze)
      endif

c     compute grad pr
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        pm1(1,1,1,e,1),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,1) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,2) = 1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,3) = 1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
      enddo

      ! div (phi_p * v)
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! x dir
     >                                        ptw(1,1,1,e,6),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! y dir
     >                                        ptw(1,1,1,e,7),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! z dir
     >                                        ptw(1,1,1,e,8),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,4) = -pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,4)
      enddo
      enddo
      enddo
      enddo

c     compute grad phi_g
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        phig(1,1,1,e),lx1,if3d)
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,5) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,6) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,7) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,5) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,5)
         rhs_fluidp(i,j,k,e,6) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,6)
         rhs_fluidp(i,j,k,e,7) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,7)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_col(i)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real mcfac,mcfac_wall, rdum, rdum3(3)

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      real rpx1(3), rpx2(3), rpx0(3),r1(3),r2(3),r3(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      real rv_t12
      common /tangential_vel/ rvx_t12
      
      ptdum(18) = dnekclock()

      mcfac       = 2.*sqrt(ksp)*log(e_rest)/
     >              sqrt(log(e_rest)**2+pi**2)
      mcfac_wall  = 2.*sqrt(ksp_wall)*log(e_rest_wall)/
     >              sqrt(log(e_rest_wall)**2+pi**2)

      if (two_way .gt. 2) then

      rrp1   = rpart(jrpe ,i)
      rvol1  = rpart(jvol ,i)
      rrho1  = rpart(jrhop,i)
      rx1(1) = rpart(jx   ,i)
      rx1(2) = rpart(jx+1 ,i)
      rx1(3) = rpart(jx+2 ,i)
      rv1(1) = rpart(jv0  ,i)
      rv1(2) = rpart(jv0+1,i)
      rv1(3) = rpart(jv0+2,i)

      icx1 = ipart(jicx,i)
      icy1 = ipart(jicy,i)
      
      icz1 = 0
      if (if3d) icz1 = ipart(jicz,i)

      icxm = icx1 -1
      icxp = icx1 +1
      icym = icy1 -1
      icyp = icy1 +1

      iczm = 0
      iczp = 0
      if (if3d) then
         iczm = icz1 -1
         iczp = icz1 +1
      endif

c     let every particle search for itself
c        particles in local elements
         do j = i+1,n
            icx2 = ipart(jicx,j)
            icy2 = ipart(jicy,j)
            icz2 = ipart(jicz,j)
            if ((icx2.ge.icxm).and.(icx2.le.icxp)) then
            if ((icy2.ge.icym).and.(icy2.le.icyp)) then
            if ((icz2.ge.iczm).and.(icz2.le.iczp)) then

            rrp2   = rpart(jrpe ,j)
            rvol2  = rpart(jvol ,j)
            rrho2  = rpart(jrhop,j)
            rx2(1) = rpart(jx   ,j)
            rx2(2) = rpart(jx+1 ,j)
            rx2(3) = rpart(jx+2 ,j)
            rv2(1) = rpart(jv0  ,j)
            rv2(2) = rpart(jv0+1,j)
            rv2(3) = rpart(jv0+2,j)

            idum = 1
            call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rpart(jfcol,j),idum)

            endif
            endif
            endif
         enddo

c        search list of ghost particles
         do j = 1,nfptsgp
            icx2 = iptsgp(jgpicx,j)
            icy2 = iptsgp(jgpicy,j)
            icz2 = iptsgp(jgpicz,j)
            if ((icx2.ge.icxm).and.(icx2.le.icxp)) then
            if ((icy2.ge.icym).and.(icy2.le.icyp)) then
            if ((icz2.ge.iczm).and.(icz2.le.iczp)) then

            rrp2   = rptsgp(jgprpe ,j)
            rvol2  = rptsgp(jgpvol ,j)
            rrho2  = rho_p            ! assume same density. Need2fix
            rx2(1) = rptsgp(jgpx   ,j)
            rx2(2) = rptsgp(jgpx+1 ,j)
            rx2(3) = rptsgp(jgpx+2 ,j)
            rv2(1) = rptsgp(jgpv0  ,j)
            rv2(2) = rptsgp(jgpv0+1,j)
            rv2(3) = rptsgp(jgpv0+2,j)

            idum = 0
            call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rdum3,idum)

            endif
            endif
            endif
         enddo
 1235 continue

         ! plane wall collisions
         do j = 1,np_walls
            rnx = plane_wall_coords(1,j)
            rny = plane_wall_coords(2,j)
            rnz = plane_wall_coords(3,j)
            rpx = plane_wall_coords(4,j)
            rpy = plane_wall_coords(5,j)
            rpz = 1.0
            if (if3d) rpz = plane_wall_coords(6,j)


            rd    = -(rnx*rpx + rny*rpy + rnz*rpz)

            rdist = abs(rnx*rpart(jx,i)+rny*rpart(jy,i)+rnz*rpart(jz,i)
     >                    +rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            rrp2   = 0.
            rvol2  = 1.
            rrho2  = 1E8
            rx2(1) = rpart(jx  ,i) - rdist*rnx
            rx2(2) = rpart(jx+1,i) - rdist*rny
            rx2(3) = rpart(jx+2,i) - rdist*rnz
            rv2(1) = 0.
            rv2(2) = 0.
            rv2(3) = 0.

            idum = 0
            call compute_collide(mcfac_wall,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rdum3,idum)
         enddo

         ! cylinder wall collisions
         do j = 1,nc_walls
            rnx = cyl_wall_coords(1,j)
            rny = cyl_wall_coords(2,j)
            rnz = cyl_wall_coords(3,j)
            rpx = cyl_wall_coords(4,j)
            rpy = cyl_wall_coords(5,j)
            rpz = 1.0
            if (if3d) rpz = cyl_wall_coords(6,j)


            rrad = cyl_wall_coords(7,j)

            rx2(1) = rpart(jx,i)
            rx2(2) = rpart(jy,i)
            rx2(3) = rpart(jz,i)
            ! for now only works with cylinders aligned with axes at
            ! origin
            if (rnz .gt. 0.5) then
               rtheta = atan2(rpart(jy,i)-rpy,rpart(jx,i)-rpx)
               rx2(1) = rpx+rrad*cos(rtheta)
               rx2(2) = rpy+rrad*sin(rtheta)
            elseif (rnx .gt. 0.5) then
               rtheta = atan2(rpart(jz,i)-rpz,rpart(jy,i)-rpy)
               rx2(2) = rpy+rrad*cos(rtheta)
               rx2(3) = rpz+rrad*sin(rtheta)
            elseif (rny .gt. 0.5) then
               rtheta = atan2(rpart(jx,i)-rpx,rpart(jz,i)-rpz)
               rx2(3) = rpz+rrad*cos(rtheta)
               rx2(1) = rpx+rrad*sin(rtheta)
            endif

            rrp2   = 0.
            rvol2  = 1.
            rrho2  = 1E8
            rv2(1) = 0.
            rv2(2) = 0.
            rv2(3) = 0.

            idum = 0
            call compute_collide(mcfac_wall,rrp2,rvol2,rrho2,rx2,rv2,
     >                           rpart(jfcol,i),rdum3,idum)
         enddo

      endif


      pttime(18) = pttime(18) + dnekclock() - ptdum(18)

      return
      end
c----------------------------------------------------------------------
      subroutine compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,fcf1,
     >                                                      fcf2,iflg)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real fcf1(3),fcf2(3),er,eta,ere,mcfac
      integer iflg

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      rthresh  = rrp1 + rrp2
      rthresh2 = rthresh**2

      rxdiff  = rx2(1) - rx1(1)
      rsum2 = rxdiff**2
      if (rsum2 .gt. rthresh2) goto 1511

      rydiff = rx2(2) - rx1(2)
      rydiff2 = rydiff**2
      rsum2 = rsum2 + rydiff2
      if (rsum2 .gt. rthresh2) goto 1511

      if (if3d) then
         rzdiff = rx2(3) - rx1(3)
         rzdiff2 = rzdiff**2
         rsum2 = rsum2 + rzdiff2
         if (rsum2 .gt. rthresh2) goto 1511
      endif

      rdiff = sqrt(rsum2)
      rm1   = rrho1*rvol1
      rm2   = rrho2*rvol2

      ! wall
      if (rm2 .gt. 1E7) then
         rksp_use = ksp_wall ! restitution coefficient
         re_rest_use = e_rest_wall
      ! other particles
      else
         rksp_use = ksp
         re_rest_use = e_rest
      endif

      rmult = 1./sqrt(1./rm1 + 1./rm2)
      eta  = mcfac*rmult

      ! first, handle normal collision part
      rbot     = 1./rdiff
      rn_12x   = rxdiff*rbot
      rn_12y   = rydiff*rbot
      rn_12z   = rzdiff*rbot

      rdelta12 = rthresh - rdiff

      rv12_mag = (rv2(1) - rv1(1))*rn_12x +
     >           (rv2(2) - rv1(2))*rn_12y +
     >           (rv2(3) - rv1(3))*rn_12z

      rv12_mage = rv12_mag*eta

      rksp_max = rksp_use*rdelta12

      rnmag = -rksp_max - rv12_mage

      rfn1 = rnmag*rn_12x
      rfn2 = rnmag*rn_12y
      rfn3 = 0.0
      if (if3d) rfn3 = rnmag*rn_12z

      fcf1(1) = fcf1(1) + rfn1
      fcf1(2) = fcf1(2) + rfn2
      fcf1(3) = fcf1(3) + rfn3

      if (iflg .eq. 1) then
          fcf2(1) = fcf2(1) - rfn1
          fcf2(2) = fcf2(2) - rfn2
          fcf2(3) = fcf2(3) - rfn3
      endif


 1511 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine send_ghost_particles
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl
      logical partl         

      !print*,"transfer",nid,nfptsgp,nigp,nrgp,jgpps
      ! send ghost particles
      call fgslib_crystal_tuple_transfer(i_cr_hndl,nfptsgp,llpart_gp
     $           , iptsgp,nigp,partl,0,rptsgp,nrgp,jgpps) ! jgpps is overwri

      ! sort by element for ghost particles for projection performance
      call fgslib_crystal_tuple_sort    (i_cr_hndl,nfptsgp
     $              , iptsgp,nigp,partl,0,rptsgp,nrgp,jgpes,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine create_extra_particles
c
c     create ghost and wall particles
c
c     bc_part = -1,1  => non-periodic search
c     bc_part = 0  => periodic search
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      character*132 deathmessage

      call create_ghost_particles_full 
      if(ibm_debug_bin.eq.1 .and. nfptsgp.gt.0)
     & print*,"ngp",nid,nfptsgp

      if(IBM_Particle_shape.ge.1) call create_ghost_particles_full_ibm
      if(ibm_debug_bin.eq.1 .and. nfptsgp_ibm.gt.0)
     & print*,"ngp ibm",nid,nfptsgp_ibm

      nmax   = iglmax(n,1)
      ngpmax = iglmax(nfptsgp,1)

      if (nmax .gt. llpart) then
         if (nid.eq.0) write(6,1) nmax, llpart, nid
         call exitt
      elseif (ngpmax .gt. llpart_gp) then
         if (nid.eq.0) write(6,2) ngpmax, llpart_gp, nid
         call exitt
      endif
    1 format('Max number of real particles:',
     >   i9,'. Not moving because llpart =',i9, ' on nid = ', i9)
    2 format('Max number of ghost particles:',
     >   i9,'. Not moving because llpart_gp =',i9, ' on nid = ', i9)

      return
      end
c----------------------------------------------------------------------
      subroutine create_ghost_particles_full
c
c     this routine will create ghost particles by checking if particle
c     is within d2chk of element faces
c
c     ghost particle x,y,z list will be in rptsgp(jgpx,j),rptsgp(jgpy,j),
c     rptsgp(jgpz,j), while processor and local element id are in
c     iptsgp(jgppt,j) and iptsgp(jgpes,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xdlen,ydlen,zdlen,rxdrng(3),rxnew(3)
      integer iadd(3),ntypesl(7)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      xdlen = xdrange(2,1) - xdrange(1,1)
      ydlen = xdrange(2,2) - xdrange(1,2)
      zdlen = -1.
      if (if3d) zdlen = xdrange(2,3) - xdrange(1,3)
      if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) xdlen = -1
      if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) ydlen = -1
      if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      nfptsgp = 0

      do ip=1,n
         rxval = rpart(jx,ip)
         ryval = rpart(jy,ip)
         rzval = 0.
         if(if3d) rzval = rpart(jz,ip)

         iip    = floor((rxval-xdrange(1,1))/rdxgp) 
         jjp    = floor((ryval-xdrange(1,2))/rdygp) 
         kkp    = floor((rzval-xdrange(1,3))/rdzgp) 
         ndump  = iip + ndxgp*jjp + ndxgp*ndygp*kkp

      do i=1,nlist
         ii = ngp_valsp(3,i)
         jj = ngp_valsp(4,i)
         kk = ngp_valsp(5,i)

         ndum = ngp_valsp(2,i)
         nrank= ngp_valsp(1,i)

         ! add this box
         if (nid  .ne. nrank) then
         if (ndum .eq. ndump) then

            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
         
            iadd(1)  = 0
            iadd(2)  = ngp_valsp(1,i)
            iadd(3)  = ipart(je0,ip)
         
            call add_a_ghost_particle(rxnew,iadd,ip)

         endif
         endif

      enddo
      enddo

      do ip=1,n
         rxval = rpart(jx,ip)
         ryval = rpart(jy,ip)
         rzval = 0.
         if(if3d) rzval = rpart(jz,ip)

         iip    = floor((rxval-xdrange(1,1))/rdxgp) 
         jjp    = floor((ryval-xdrange(1,2))/rdygp) 
         kkp    = floor((rzval-xdrange(1,3))/rdzgp) 
         ndump  = iip + ndxgp*jjp + ndxgp*ndygp*kkp

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
            if (iflgsum .ne. 1) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg


            do i=1,nlist
               ndum = ngp_valsp(2,i)
               nrank = ngp_valsp(1,i)
               iin   = ngp_valsp(3,i)
               jjn   = ngp_valsp(4,i)
               kkn   = ngp_valsp(5,i)

               if (ndumn .eq. ndum) then
                  ibctype = iflgx+iflgy+iflgz
                 
                  if (nrank .ne. nid .or. 
     >                    (nrank .eq. nid) .and. (iflgsum.eq.1)) then
                 
                      rxnew(1) = rxval
                      rxnew(2) = ryval
                      rxnew(3) = rzval
                 
                      iadd(1) = ii1
                      iadd(2) = jj1
                      iadd(3) = kk1
                 
                      call check_periodic_gp(rxnew,rxdrng,iadd)
                 
                      iadd(1)  = 0
                      iadd(2)  = nrank
                      iadd(3)  = ipart(je0,ip)
                      
                      call add_a_ghost_particle(rxnew,iadd,ip)
                        
                  endif
               endif
            enddo
         enddo
         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
            if (iflgsum .ne. 2) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg


            do i=1,nlist
               ndum = ngp_valsp(2,i)
               nrank = ngp_valsp(1,i)
               iin   = ngp_valsp(3,i)
               jjn   = ngp_valsp(4,i)
               kkn   = ngp_valsp(5,i)

               if (ndumn .eq. ndum) then
                  ibctype = iflgx+iflgy+iflgz
                 
                  if (nrank .ne. nid .or. 
     >                    (nrank .eq. nid) .and. (iflgsum.eq.2)) then
                 
                      rxnew(1) = rxval
                      rxnew(2) = ryval
                      rxnew(3) = rzval
                 
                      iadd(1) = ii1
                      iadd(2) = jj1
                      iadd(3) = kk1
                 
                      call check_periodic_gp(rxnew,rxdrng,iadd)
                 
                      iadd(1)  = 0
                      iadd(2)  = nrank
                      iadd(3)  = ipart(je0,ip)
                      
                      call add_a_ghost_particle(rxnew,iadd,ip)
                        
                  endif
               endif
            enddo
         enddo
         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif


            iflgsum = iflgx + iflgy + iflgz
            if (iflgsum .ne. 3) cycle

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg


            do i=1,nlist
               ndum = ngp_valsp(2,i)
               nrank = ngp_valsp(1,i)
               iin   = ngp_valsp(3,i)
               jjn   = ngp_valsp(4,i)
               kkn   = ngp_valsp(5,i)

               if (ndumn .eq. ndum) then
                  ibctype = iflgx+iflgy+iflgz
                 
                  if (nrank .ne. nid .or. 
     >                    (nrank .eq. nid) .and. (iflgsum.eq.3)) then
                 
                      rxnew(1) = rxval
                      rxnew(2) = ryval
                      rxnew(3) = rzval
                 
                      iadd(1) = ii1
                      iadd(2) = jj1
                      iadd(3) = kk1
                 
                      call check_periodic_gp(rxnew,rxdrng,iadd)
                 
                      iadd(1)  = 0
                      iadd(2)  = nrank
                      iadd(3)  = ipart(je0,ip)
                      
                      call add_a_ghost_particle(rxnew,iadd,ip)
                        
                  endif
               endif
            enddo
         enddo
      enddo
c      print*, "ngp is ", nid, nfptsgp

      return
      end
c-----------------------------------------------------------------------
      subroutine check_periodic_gp(rxnew,rxdrng,iadd)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
c
      real rxnew(3), rxdrng(3)
      integer iadd(3), irett(3), ntype, ntypel(7)

      xloc = rxnew(1)
      yloc = rxnew(2)
      zloc = rxnew(3)

      xdlen = rxdrng(1)
      ydlen = rxdrng(2)
      zdlen = rxdrng(3)

      ii = iadd(1)
      jj = iadd(2)
      kk = iadd(3)

      irett(1) = 0
      irett(2) = 0
      irett(3) = 0

      if (xdlen .gt. 0 ) then
      if (ii .ge. ndxgp) then
         xloc = xloc - xdlen
         irett(1) = 1
         goto 123
      endif
      endif
      if (xdlen .gt. 0 ) then
      if (ii .lt. 0) then
         xloc = xloc + xdlen
         irett(1) = 1
         goto 123
      endif
      endif

  123 continue    
      if (ydlen .gt. 0 ) then
      if (jj .ge. ndygp) then
         yloc = yloc - ydlen
         irett(2) = 1
         goto 124
      endif
      endif
      if (ydlen .gt. 0 ) then
      if (jj .lt. 0) then
         yloc = yloc + ydlen
         irett(2) = 1
         goto 124
      endif
      endif
  124 continue

      if (if3d) then
         if (zdlen .gt. 0 ) then
         if (kk .ge. ndzgp) then
            zloc = zloc - zdlen
            irett(3) = 1
            goto 125
         endif
         endif
         if (zdlen .gt. 0 ) then
         if (kk .lt. 0) then
            zloc = zloc + zdlen
            irett(3) = 1
            goto 125
         endif
         endif
      endif
  125 continue

      rxnew(1) = xloc
      rxnew(2) = yloc
      rxnew(3) = zloc

      return
      end
c----------------------------------------------------------------------
      subroutine add_a_ghost_particle(rxnew,iadd,i)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real    rxnew(3)
      integer iadd(3)

      nfptsgp = nfptsgp + 1

      rptsgp(jgpx,nfptsgp)    = rxnew(1)       ! x loc
      rptsgp(jgpy,nfptsgp)    = rxnew(2)       ! y loc
      rptsgp(jgpz,nfptsgp)    = rxnew(3)       ! z loc
      rptsgp(jgpfh,nfptsgp)   = rpart(jf0,i)   ! hyd. force x
      rptsgp(jgpfh+1,nfptsgp) = rpart(jf0+1,i) ! hyd. force y
      rptsgp(jgpfh+2,nfptsgp) = rpart(jf0+2,i) ! hyd. force z
      rptsgp(jgpvol,nfptsgp)  = rpart(jvol,i)  ! particle volum
      rptsgp(jgprpe,nfptsgp)  = rpart(jrpe,i)  ! particle rp eff
      rptsgp(jgpspl,nfptsgp)  = rpart(jspl,i)  ! spl
      rptsgp(jgpg0,nfptsgp)   = rpart(jg0,i)   ! work done by forc
      rptsgp(jgpq0,nfptsgp)   = rpart(jq0,i)   ! heating from part 
      rptsgp(jgpv0,nfptsgp)   = rpart(jv0,i)   ! particle velocity
      rptsgp(jgpv0+1,nfptsgp) = rpart(jv0+1,i) ! particle velocity
      rptsgp(jgpv0+2,nfptsgp) = rpart(jv0+2,i) ! particle velocity

      iptsgp(jgpiic,nfptsgp)  = iadd(1)        ! use in collisions
      iptsgp(jgpps,nfptsgp)   = iadd(2)        ! overwritten mpi
      iptsgp(jgppt,nfptsgp)   = iadd(2)        ! dest. mpi rank
      iptsgp(jgpes,nfptsgp)   = iadd(3)        ! dest. elment
      
      
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_neighbor_el_proc
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl

      real pfx,pfy,pfz,vol,qgqf,rvx,rvy,rvz,rexp,multfc,multfci,rx2(3)
     >     ,rxyzp(6,3)
      integer ntypesl(7), ngp_trim(nbox_gp)

      real    rxnew(3), rxdrng(3)
      integer iadd(3)

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (if3d) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

! -------------------------------------------------------
c SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------
      ! how many spacings in each direction
      ndxgp = floor( (xdrange(2,1) - xdrange(1,1))/d2chk(1))+1
      ndygp = floor( (xdrange(2,2) - xdrange(1,2))/d2chk(1))+1
      ndzgp = 1
      if (if3d) ndzgp = floor( (xdrange(2,3) - xdrange(1,3))/d2chk(1))+1

      nreach = floor(real(ndxgp*ndygp*ndzgp)/real(np)) + 1 ! for IBM Yang 

      ! grid spacing for that many spacings
      rdxgp = (xdrange(2,1) - xdrange(1,1))/real(ndxgp)
      rdygp = (xdrange(2,2) - xdrange(1,2))/real(ndygp)
      rdzgp = 1.
      if (if3d) rdzgp = (xdrange(2,3) - xdrange(1,3))/real(ndzgp)

! ------------------------------------------------------------
c Connect boxes to 1D processor map they should be arranged on
! ------------------------------------------------------------
      ! Add my own box and what rank(s) it is in the 1D proc map
      nlist = 0
      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rxval = xm1(i,j,k,ie)
         ryval = ym1(i,j,k,ie)
         rzval = 0.
         if(if3d) rzval = zm1(i,j,k,ie)

c        if (rxval .gt. rxbo(2,1) .or. 
c    >       rxval .lt. rxbo(1,1) .or. 
c    >       ryval .gt. rxbo(2,2) .or. 
c    >       ryval .lt. rxbo(1,2) ) then
c            goto 1234
c        endif
c        if (if3d) then
c        if (rzval .gt. rxbo(2,3) .or. 
c    >       rzval .lt. rxbo(1,3) ) then
c            goto 1234
c        endif
c        endif

         ii    = floor((rxval-xdrange(1,1))/rdxgp) 
         jj    = floor((ryval-xdrange(1,2))/rdygp) 
         kk    = floor((rzval-xdrange(1,3))/rdzgp) 
         if (ii .eq. ndxgp) ii = ndxgp - 1
         if (jj .eq. ndygp) jj = ndygp - 1
         if (kk .eq. ndzgp) kk = ndzgp - 1
         ndum  = ii + ndxgp*jj + ndxgp*ndygp*kk

         mod_gp_grid(i,j,k,ie,1) = ii
         mod_gp_grid(i,j,k,ie,2) = jj
         mod_gp_grid(i,j,k,ie,3) = kk
         mod_gp_grid(i,j,k,ie,4) = ndum

         nlist = nlist + 1
         if (nlist .gt. nbox_gp) then
            write(6,*)'Increase nbox_gp. Need more sub-box storage',
     $                      nbox_gp, nlist, nid, ie
            call exitt
         endif

         ngp_valsp(1,nlist) = nid
         ngp_valsp(2,nlist) = ndum
         ngp_valsp(3,nlist) = ii
         ngp_valsp(4,nlist) = jj
         ngp_valsp(5,nlist) = kk
         ngp_valsp(6,nlist) = floor(real(ndum)/real(nreach))

         if (nlist .gt. 1) then
         do il=1,nlist-1
            if (ngp_valsp(2,il) .eq. ndum) then
               nlist = nlist - 1
               goto 1234
            endif
         enddo
         endif
 1234 continue
      enddo
      enddo
      enddo
      enddo
      
      if (nid.eq.0 .and. ibm_debug_bin.eq.1) then
         do i = 1, nlist
         write(6,'(A,7I6)')"ngp_valsp inner",(ngp_valsp(j,i),j=1,6)
         enddo
      endif

      ! Add connecting boxes and what rank(s) they are in the 1D proc map
      nlist_save = nlist
      do i=1,nlist_save
         ii = ngp_valsp(3,i)
         jj = ngp_valsp(4,i)
         kk = ngp_valsp(5,i)

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = ii + el_face_num(ist+1) 
            jj1 = jj + el_face_num(ist+2)
            kk1 = kk + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            ! periodic if out of domain
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) iig =modulo(iig,ndxgp)
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) jjg =modulo(jjg,ndygp)
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) kkg =modulo(kkg,ndzgp)

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg

            do j=1,nlist
               if (ngp_valsp(2,j) .eq. ndumn) goto 999
            enddo

            nlist = nlist + 1

            ngp_valsp(1,nlist) = nid
            ngp_valsp(2,nlist) = ndumn
            ngp_valsp(6,nlist) = floor(real(ndumn)/real(nreach))
            ! periodic bins
            if (ii1 .gt. ndxgp-1) then
               ii1 = -1
            elseif (ii1 .lt. 0) then
               ii1 = ndxgp
            endif
            if (jj1 .gt. ndygp-1) then
               jj1 = -1
            elseif (jj1 .lt. 0) then
               jj1 = ndygp
            endif
            if (kk1 .gt. ndzgp-1) then
               kk1 = -1
            elseif (kk1 .lt. 0) then
               kk1 = ndzgp
            endif

            ngp_valsp(3,nlist) = ii1
            ngp_valsp(4,nlist) = jj1
            ngp_valsp(5,nlist) = kk1

            if(nid.eq.0 .and.ibm_debug_bin.eq.1)
     $       write(6,'(A,6I5)')"ngp_valsp face",
     $        (ngp_valsp(j,nlist),j=1,6)            

  999 continue
         enddo
         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = ii + el_edge_num(ist+1) 
            jj1 = jj + el_edge_num(ist+2)
            kk1 = kk + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            ! periodic if out of domain
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) iig =modulo(iig,ndxgp)
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) jjg =modulo(jjg,ndygp)
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) kkg =modulo(kkg,ndzgp)

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg

            do j=1,nlist
               if (ngp_valsp(2,j) .eq. ndumn) goto 998
            enddo

            nlist = nlist + 1

            ngp_valsp(1,nlist) = nid
            ngp_valsp(2,nlist) = ndumn
            ngp_valsp(6,nlist) = floor(real(ndumn)/real(nreach))

            if (ii1 .gt. ndxgp-1) then
               ii1 = -1
            elseif (ii1 .lt. 0) then
               ii1 = ndxgp
            endif
            if (jj1 .gt. ndygp-1) then
               jj1 = -1
            elseif (jj1 .lt. 0) then
               jj1 = ndygp
            endif
            if (kk1 .gt. ndzgp-1) then
               kk1 = -1
            elseif (kk1 .lt. 0) then
               kk1 = ndzgp
            endif

            ngp_valsp(3,nlist) = ii1
            ngp_valsp(4,nlist) = jj1
            ngp_valsp(5,nlist) = kk1

            if(nid.eq.0.and.ibm_debug_bin.eq.1)
     $           write(6,'(A,7I6)')"ngp_valsp edge"
     $        ,(ngp_valsp(j,nlist),j=1,6)

  998 continue
         enddo
         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = ii + el_corner_num(ist+1) 
            jj1 = jj + el_corner_num(ist+2)
            kk1 = kk + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            ! periodic if out of domain
            if (iig .lt. 0 .or. iig .gt. ndxgp-1) iig =modulo(iig,ndxgp)
            if (jjg .lt. 0 .or. jjg .gt. ndygp-1) jjg =modulo(jjg,ndygp)
            if (kkg .lt. 0 .or. kkg .gt. ndzgp-1) kkg =modulo(kkg,ndzgp)

            ndumn = iig + ndxgp*jjg + ndxgp*ndygp*kkg

            do j=1,nlist
               if (ngp_valsp(2,j) .eq. ndumn) goto 997
            enddo

            nlist = nlist + 1

            ngp_valsp(1,nlist) = nid
            ngp_valsp(2,nlist) = ndumn
            ngp_valsp(6,nlist) = floor(real(ndumn)/real(nreach))

            if (ii1 .gt. ndxgp-1) then
               ii1 = -1
            elseif (ii1 .lt. 0) then
               ii1 = ndxgp
            endif
            if (jj1 .gt. ndygp-1) then
               jj1 = -1
            elseif (jj1 .lt. 0) then
               jj1 = ndygp
            endif
            if (kk1 .gt. ndzgp-1) then
               kk1 = -1
            elseif (kk1 .lt. 0) then
               kk1 = ndzgp
            endif

            ngp_valsp(3,nlist) = ii1
            ngp_valsp(4,nlist) = jj1
            ngp_valsp(5,nlist) = kk1

            if(nid.eq.0 .and. ibm_debug_bin.eq.1)
     $       write(6,'(A,7I6)')"ngp_valsp corner"
     $        ,(ngp_valsp(j,nlist),j=1,6)            

  997 continue
         enddo
      enddo

! ------------------------
c SEND TO 1D PROCESSOR MAP
! ------------------------
      nps   = 6 ! index of new proc for doing stuff
      nglob = 2 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      call fgslib_crystal_ituple_transfer(i_cr_hndl,ngp_valsp,
     >                 ngpvc,nlist,nbox_gp,nps)
      call fgslib_crystal_ituple_sort(i_cr_hndl,ngp_valsp,
     >                 ngpvc,nlist,nglob,nkey)

      ! trim down so no duplicates
      do i=1,nlist
         ngp_trim(i) = 1
      enddo
      do i=1,nlist
      do j=i+1,nlist
         if (ngp_trim(j) .eq. 1) then
            if (ngp_valsp(1,i) .eq. ngp_valsp(1,j) ) then
            if (ngp_valsp(2,i) .eq. ngp_valsp(2,j) ) then
               ngp_trim(j) = 0
            endif
            endif
         endif
      enddo
      enddo
      ic = 0
      do i=1,nlist
         if (ngp_trim(i) .eq. 1) then
            ic = ic + 1
            call icopy(ngp_valsp(1,ic),ngp_valsp(1,i),ngpvc)
          endif
      enddo

      nlist = ic

      ! create dupicates to send to remote processors
      ! how to understand:  
      ! 1) bin2rank record1 in rank (i) contains (igbl) bin
      ! 2) search for record2 rank(j) contains same (jgbl = igbl) bin
      ! 3) send record1(rank(i)) to rank(j)

      nlist_save = nlist
      do i=1,nlist_save
         irnk = ngp_valsp(1,i)
         igbl = ngp_valsp(2,i)

         do j=1,nlist_save
            jrnk = ngp_valsp(1,j)
            jgbl = ngp_valsp(2,j)
            if (i .eq. j) cycle
            if (irnk .eq. jrnk) cycle  !! in different rank
            if (igbl .ne. jgbl) cycle  !! in the same bin
            nlist = nlist + 1                               ! create new record
            if (nlist .gt. nbox_gp) then
               write(6,*)'Increase nbox_gp. In dup. loop',
     $                         nbox_gp, nlist, nid, nlist_save
               call exitt
            endif
            call icopy(ngp_valsp(1,nlist),ngp_valsp(1,i),ngpvc) !!! create connection between irnk & jrnk in the same bin in the new record
            ngp_valsp(6,nlist) = ngp_valsp(1,j)                 !!! remote rank contains this bin to rank connection
         enddo
      enddo

! ----------------------------------------------
c SEND BACK TO ALL PROCESSORS WITH ADDITIONS NOW
! ----------------------------------------------
      nps   = 6 ! index of new proc for doing stuff
      nglob = 2 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      call fgslib_crystal_ituple_transfer(i_cr_hndl,ngp_valsp,
     >                 ngpvc,nlist,nbox_gp,nps)
      call fgslib_crystal_ituple_sort(i_cr_hndl,ngp_valsp,
     >                 ngpvc,nlist,nglob,nkey)

      return
      end

c-----------------------------------------------------------------------
      subroutine extra_wall_particle_exp(rxyzp,rx2,ic)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
c
      real rxyzp(n_walls*2,3), rx2(3), rx22(3)
      integer ic, ip

      ! plane wall collisions
      ic = 0
      do j = 1,np_walls
         rnx = plane_wall_coords(1,j)
         rny = plane_wall_coords(2,j)
         rnz = plane_wall_coords(3,j)
         rpx = plane_wall_coords(4,j)
         rpy = plane_wall_coords(5,j)
         rpz = 1.0
         if (if3d) rpz = plane_wall_coords(6,j)

         rd    = -(rnx*rpx + rny*rpy + rnz*rpz)

         rdist = abs(rnx*rx2(1)+rny*rx2(2)+rnz*rx2(3)+rd)
         rdist = rdist*2.
         rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

         if (rdist .gt. d2chk(1)) cycle
         ic = ic + 1

         rxyzp(ic,1) = rx2(1) - rdist*rnx
         rxyzp(ic,2) = rx2(2) - rdist*rny
         rxyzp(ic,3) = rx2(3) - rdist*rnz
      enddo

      ! cylinder wall collisions
      do j = 1,nc_walls
         rnx = cyl_wall_coords(1,j)
         rny = cyl_wall_coords(2,j)
         rnz = cyl_wall_coords(3,j)
         rpx = cyl_wall_coords(4,j)
         rpy = cyl_wall_coords(5,j)
         rpz = 1.0
         if (if3d) rpz = cyl_wall_coords(6,j)
         rrad = cyl_wall_coords(7,j)

         rx22(1) = rx2(1)
         rx22(2) = rx2(2)
         rx22(3) = rx2(3)
         ! for now only works with cylinders aligned with axes at
         ! origin
         if (rnz .gt. 0.5) then
            rtheta = atan2(rx22(2),rx22(1))
            rx22(1) = rpx + rrad*cos(rtheta)
            rx22(2) = rpy + rrad*sin(rtheta)
         elseif (rnx .gt. 0.5) then
            rtheta = atan2(rx22(3),rx22(2))
            rx22(2) = rpy + rrad*cos(rtheta)
            rx22(3) = rpz + rrad*sin(rtheta)
         elseif (rny .gt. 0.5) then
            rtheta = atan2(rx22(1),rx22(3))
            rx22(3) = rpz + rrad*cos(rtheta)
            rx22(1) = rpx + rrad*sin(rtheta)
         endif

         rx2d = rx22(1) - rx2(1)
         ry2d = rx22(2) - rx2(2)
         rz2d = rx22(3) - rx2(3)

         rdist = sqrt(rx2d**2 + ry2d**2 + rz2d**2)
         rx2d = rx2d/rdist
         ry2d = ry2d/rdist
         rz2d = rz2d/rdist

         rdist = rdist*2.

         if (rdist .gt. d2chk(1)) cycle
         ic = ic + 1

         rxyzp(ic,1) = rx2(1) + rx2d*rdist
         rxyzp(ic,2) = rx2(2) + ry2d*rdist
         rxyzp(ic,3) = rx2(3) + rz2d*rdist

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_bdf_ext_coefs(beta,alpha,times)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real beta(0:3),alpha(0:3),times(0:3)
      real c(0:8)

      integer ilast,ncoef
      save    ilast,ncoef
      data    ilast,ncoef / -9 , 0 /

      do i=3,1,-1
         times(i)=times(i-1)
      enddo
      times(0) = time

      call rzero(beta ,4)
      call rzero(alpha,4)
      if (istep.ne.ilast) then
         ilast = istep
         ncoef = ncoef + 1
         ncoef = min(ncoef,3) ! Maximum 3rd order in time
      endif
      ncoefm1 = ncoef - 1

      call fd_weights_full(times(0),times(1),ncoefm1,0,alpha(1))
      call fd_weights_full(times(0),times(0),ncoef,1,c)
      do j=0,ncoef
         beta(j) = c(ncoef+1+j)
      enddo
      do j=1,ncoef
         beta(j) = -beta(j)  ! Change sign, for convenience
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine update_particle_location
c     check if particles are outside domain
c     > if bc_part = 0 then it is periodic
c     > if bc_part = -1,1 then particles are killed (outflow)
      include 'SIZE'
      include 'LPM'

      integer in_part(llpart)


      do i=1,n
         in_part(i) = 0
         do j=0,ndim-1
            if (rpart(jx+j,i).lt.xdrange(1,j+1))then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx+j,i) =  xdrange(2,j+1) - 
     &                             abs(xdrange(1,j+1) - rpart(jx+j,i))
                  rpart(jx1+j,i) = xdrange(2,j+1) -
     &                             abs(xdrange(1,j+1) - rpart(jx1+j,i))
                  rpart(jx2+j,i) = xdrange(2,j+1) -
     &                             abs(xdrange(1,j+1) - rpart(jx2+j,i))
                  rpart(jx3+j,i) = xdrange(2,j+1) -
     &                             abs(xdrange(1,j+1) - rpart(jx3+j,i))
                  goto 1512
                endif
            endif
            if (rpart(jx+j,i).gt.xdrange(2,j+1))then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx+j,i) =  xdrange(1,j+1) +
     &                             abs(rpart(jx+j,i)  - xdrange(2,j+1))
                  rpart(jx1+j,i) = xdrange(1,j+1) +
     &                             abs(rpart(jx1+j,i) - xdrange(2,j+1))
                  rpart(jx2+j,i) = xdrange(1,j+1) +
     &                             abs(rpart(jx2+j,i) - xdrange(2,j+1))
                  rpart(jx3+j,i) = xdrange(1,j+1) + 
     &                             abs(rpart(jx3+j,i) - xdrange(2,j+1))
                  goto 1512
                endif
            endif
            if (ipart(jrc,i) .eq. 2) then
               in_part(i) = -1 ! only if periodic check fails it will get here
            endif
 1512 continue
         enddo
      enddo

      ic = 0
      do i=1,n
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               call copy(rpart(1,ic),rpart(1,i),nr)
               call icopy(ipart(1,ic),ipart(1,i),ni)
            endif
         endif
      enddo
      n = ic

      return
      end

c----------------------------------------------------------------------
      subroutine update_queen_periodic
c     check if particles are outside domain
c     > if bc_part = 0 then it is periodic
c     > if bc_part = -1,1 then particles are killed (outflow)
      include 'SIZE'
      include 'LPM'

      integer in_part(llpart)

      common /PARTRK3/ kv_stage_p
      real   kv_stage_p(llpart,13) ! add 6 for rotation

      do i=1,n
         if(ipart(jrole,i).ne.1) cycle
         in_part(i) = 0
         do j=0,ndim-1
            if (rpart(jx+j,i).lt.xdrange(1,j+1))then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx+j,i) =  xdrange(2,j+1) - 
     &                             abs(xdrange(1,j+1) - rpart(jx+j,i))
                  rpart(jx1+j,i) = xdrange(2,j+1) -
     &                             abs(xdrange(1,j+1) - rpart(jx1+j,i))
                  rpart(jx2+j,i) = xdrange(2,j+1) -
     &                             abs(xdrange(1,j+1) - rpart(jx2+j,i))
                  rpart(jx3+j,i) = xdrange(2,j+1) -
     &                             abs(xdrange(1,j+1) - rpart(jx3+j,i))

                  kv_stage_p(i,j+1) = xdrange(2,j+1) -
     &                             abs(xdrange(1,j+1)-kv_stage_p(i,j+1))
                  goto 1515
                endif
            endif
            if (rpart(jx+j,i).gt.xdrange(2,j+1))then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx+j,i) =  xdrange(1,j+1) +
     &                             abs(rpart(jx+j,i)  - xdrange(2,j+1))
                  rpart(jx1+j,i) = xdrange(1,j+1) +
     &                             abs(rpart(jx1+j,i) - xdrange(2,j+1))
                  rpart(jx2+j,i) = xdrange(1,j+1) +
     &                             abs(rpart(jx2+j,i) - xdrange(2,j+1))
                  rpart(jx3+j,i) = xdrange(1,j+1) + 
     &                             abs(rpart(jx3+j,i) - xdrange(2,j+1))

                  kv_stage_p(i,j+1) = xdrange(1,j+1) + 
     &                          abs(kv_stage_p(i,j+1) - xdrange(2,j+1))

                  goto 1515
                endif
            endif
            if (ipart(jrc,i) .eq. 2) then
               in_part(i) = -1 ! only if periodic check fails it will get here
            endif
 1515   continue
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c     interpolation routines
c-----------------------------------------------------------------------
      subroutine tri_interp(ii,jj,kk,x,y,z,rset,xx,yy,zz,fld)
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'LPM'
      include 'GEOM'

      real xx(lx1,ly1,lz1),yy(lx1,ly1,lz1),zz(lx1,ly1,lz1)
     >    ,fld(lx1,ly1,lz1)

      iip = ii +1
      jjp = jj +1
      kkp = kk
      if (if3d) kkp = kk +1

      xd = (x - xx(ii,jj,kk))/(xx(iip,jj,kk)-xx(ii,jj,kk))
      yd = (y - yy(ii,jj,kk))/(yy(ii,jjp,kk)-yy(ii,jj,kk))
      zd = 0.0
      if(if3d) zd = (z - zz(ii,jj,kk))/(zz(ii,jj,kkp)-zz(ii,jj,kk))

      c00=fld(ii,jj,kk)*(1.-xd)+fld(iip,jj,kk)*xd
      c10=fld(ii,jjp,kk)*(1.-xd)+fld(iip,jjp,kk)*xd

      if (if3d) then
      c01=fld(ii,jj,kkp)*(1.-xd)+fld(iip,jj,kkp)*xd
      c11=fld(ii,jjp,kkp)*(1.-xd)+fld(iip,jjp,kkp)*xd
      c1_1 = c01*(1.-yd) + c11*yd
      rset = c1_0*(1.-zd) + c1_1*zd
      else
      rset = c00*(1.-yd) + c10*yd
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_props_part_location
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'LPM'
      include 'GEOM'

      common /myparth/ i_fp_hndl, i_cr_hndl

      if (red_interp .eq. 0) then

      elseif (red_interp .eq. 1) then
      ! note that all data is local now
c         write(*,*) "Velocity interpolation"
         if(interm_vel .ge. 1) then !!! ibm,
!     lagrange interpolation

       call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+0,1),nr,
     &                                   ipart(je0,1)  ,ni,
     &                                   rpart(jr,1)   ,nr,n,
     &                                   vx_tilde)
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+1,1),nr,
     &                                   ipart(je0,1)  ,ni,
     &                                   rpart(jr,1)   ,nr,n,
     &                                   vy_tilde)
      if (if3d)
     &call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+2,1),nr,
     &                                   ipart(je0,1)  ,ni,
     &                                   rpart(jr,1)   ,nr,n,
     &                                   vz_tilde)

         else !!!
            
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+0,1),nr,
     &                                   ipart(je0,1)  ,ni,
     &                                   rpart(jr,1)   ,nr,n,
     &                                   vx)
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+1,1),nr,
     &                                   ipart(je0,1)  ,ni,
     &                                   rpart(jr,1)   ,nr,n,
     &                                   vy)
      if (if3d)
     &call fgslib_findpts_eval_local(i_fp_hndl,rpart(ju0+2,1),nr,
     &                                   ipart(je0,1)  ,ni,
     &                                   rpart(jr,1)   ,nr,n,
     &                                   vz)

      endif !interm_vel_explicit
      
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(jtempf,1),nr,
     &                                   ipart(je0,1)   ,ni,
     &                                   rpart(jr,1)    ,nr,n,
     &                                   t)
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(jrho,1)  ,nr,
     &                                   ipart(je0,1)   ,ni,
     &                                   rpart(jr,1)    ,nr,n,
     &                                   vtrans)
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(jDuDt+0,1) ,nr,
     &                                   ipart(je0,1)     ,ni,
     &                                   rpart(jr,1)      ,nr,n,
     &                                   rhs_fluidp(1,1,1,1,1))
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(jDuDt+1,1) ,nr,
     &                                   ipart(je0,1)     ,ni,
     &                                   rpart(jr,1)      ,nr,n,
     &                                   rhs_fluidp(1,1,1,1,2))
      if (if3d)
     &call fgslib_findpts_eval_local(i_fp_hndl,rpart(jDuDt+2,1) ,nr,
     &                                   ipart(je0,1)     ,ni,
     &                                   rpart(jr,1)      ,nr,n,
     &                                   rhs_fluidp(1,1,1,1,3))
      call fgslib_findpts_eval_local(i_fp_hndl,rpart(jvol1,1)   ,nr,
     &                                   ipart(je0,1)     ,ni,
     &                                   rpart(jr,1)      ,nr,n,
     &                                   ptw(1,1,1,1,4))

      ! trilinear interpolation below
c     do ip=1,n
c        ie = ipart(je0,ip) + 1
c        ii = 1
c        jj = 1
c        kk = 1
c        do k=1,nz1
c        do j=1,ny1
c        do i=1,nx1
c           if (xm1(i,j,k,ie) .lt. rpart(jx,ip)) ii = i
c           if (ym1(i,j,k,ie) .lt. rpart(jy,ip)) jj = j
c           if (zm1(i,j,k,ie) .lt. rpart(jz,ip)) zz = k
c        enddo
c        enddo
c        enddo
c        call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(ju0+0,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),vx(1,1,1,ie))
c        call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(ju0+1,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),vy(1,1,1,ie))
c        if (if3d)
c    >   call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(ju0+2,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),vz(1,1,1,ie))
c        call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(jtempf,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),t(1,1,1,ie,1))
c        call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(jrho,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),vtrans(1,1,1,ie,1))
c        call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(jDuDt+0,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),rhs_fluidp(1,1,1,ie,1))
c        call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(jDuDt+1,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),rhs_fluidp(1,1,1,ie,2))
c        if (if3d)
c    >   call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(jDuDt+2,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),rhs_fluidp(1,1,1,ie,3))
c        call tri_interp(ii,jj,kk,rpart(jx,ip),rpart(jy,ip),
c    >       rpart(jz,ip),rpart(jvol1,ip),xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >       zm1(1,1,1,ie),ptw(1,1,1,ie,4))
c     enddo

      ! check if values outside reasonable bounds
      rvfmax = 0.7405
      rvfmin = 0.0
      do i=1,n
         if (rpart(jvol1,i) .lt. rvfmin) rpart(jvol1,i) = rvfmin
         if (rpart(jvol1,i) .gt. rvfmax) rpart(jvol1,i) = rvfmax
      enddo

      endif

      return
      end
c----------------------------------------------------------------------
c     particle input/output/ restart routines
c----------------------------------------------------------------------
      subroutine lpm_usr_particles_io
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'LPM'

      rdumt = dnekclock()

c     always output fields if 2 or 4 way coupled 
      if (two_way .gt. 1) then
         call output_two_way_io
      endif

c     output diagnostics to logfile
      if (npio_method .lt. 0) then
         call output_particle_timers 
         call output_particle_diagnostics
      endif

      call output_parallel_part
      call output_parallel_restart_part

      pttime(1) = pttime(1) - (dnekclock() - rdumt)

      return
      end
c----------------------------------------------------------------------
      subroutine output_parallel_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      common /myparth/ i_fp_hndl, i_cr_hndl

      character*20 vtufile
      character*14 vtufile1
      character*50 dumstr

      integer icalld
      save    icalld
      data    icalld  /0/

      logical partl         

      integer vtu,vtu1,prevs(2,np)
      integer*4 iint
      real stride_len

      icalld = icalld+1

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,np
         prevs(1,i) = i-1
         prevs(2,i) = n
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call fgslib_crystal_ituple_transfer(i_cr_hndl,prevs,
     >                 ndum,np,np,nps)
      call fgslib_crystal_ituple_sort(i_cr_hndl,prevs,
     >                 ndum,np,nglob,nkey)

      stride_len = 0
      if (nid .ne. 0) then
      do i=1,nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! -----------------------------------------------------
! SEND PARTICLES TO A SMALLER NUMBER OF RANKS TO OUTPUT
! -----------------------------------------------------
      do i = 1,n
         idum = int((stride_len + i)/llpart)
         ipart(jps,i) = idum
      enddo
      nl = 0
      call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     >                  , ipart,ni,partl,nl,rpart,nr,jps)

! -------------------
! GET TOTAL PARTICLES
! -------------------
      nptot = iglsum(n,1)
      npmax_set = int(nptot/llpart) + 1 ! use as few procs as possible
      npmax = np
      if (npmax .gt. npmax_set) npmax = npmax_set

! -----------------------
! WRITE PARALLEL VTU FILE
! -----------------------
      if (nid.eq.0) then
         vtu1 = 107
         write(vtufile1,'(A4,I5.5,A5)')'part',icalld,'.pvtu'
         open(unit=vtu1,file=vtufile1,status='replace',
     >             access='stream',form="unformatted")

         write(vtu1) '<VTKFile '
         write(vtu1) 'type="PUnstructuredGrid" '
         write(vtu1) 'version="1.0" '
         if (lpm_endian .eq. 0) then
            write(vtu1) 'byte_order="LittleEndian"> '
         elseif (lpm_endian .eq. 1) then
            write(vtu1) 'byte_order="BigEndian"> '
         endif
         
         write(vtu1) '<PUnstructuredGrid GhostLevel="0">'

         write(vtu1) '<PPoints> '
         call vtu_write_dataarray(vtu1,"Position    ",3,0,0)
         write(vtu1) '</PPoints> '
         write(vtu1) '<PPointData> '
         call vtu_write_dataarray(vtu1,"VelocityP   ",3,0,0)
         call vtu_write_dataarray(vtu1,"VelocityF   ",3,0,0)
         call vtu_write_dataarray(vtu1,"ForceUSR    ",3,0,0)
         if (part_force(1) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceQS     ",3,0,0)
         if (part_force(2) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceUN     ",3,0,0)
         if (part_force(3) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceIU     ",3,0,0)
         if (two_way .gt. 2)
     >      call vtu_write_dataarray(vtu1,"ForceC      ",3,0,0)
         if (part_force(1) .ne. 0 .or. part_force(3) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"ForceWork   ",1,0,0)
         if (part_force(4) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"HeatQS      ",1,0,0)
         if (part_force(5) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"HeatUN      ",1,0,0)
         if (part_force(4) .ne. 0 .and. part_force(5) .ne. 0)
     >      call vtu_write_dataarray(vtu1,"TemperatureP",1,0,0)
         call vtu_write_dataarray(vtu1,"TemperatureF",1,0,0)
         call vtu_write_dataarray(vtu1,"DensityF    ",1,0,0)
         if (two_way .ge. 2)
     >      call vtu_write_dataarray(vtu1,"ParticleVF  ",1,0,0)
         call vtu_write_dataarray(vtu1,"RadiusCG    ",1,0,0)
         call vtu_write_dataarray(vtu1,"Diameter    ",1,0,0)
         call vtu_write_dataarray(vtu1,"ID_1        ",1,0,0)
         call vtu_write_dataarray(vtu1,"ID_2        ",1,0,0)
         call vtu_write_dataarray(vtu1,"ID_3        ",1,0,0)
         write(vtu1) '</PPointData> '

         write(vtu1) '<PCells> '
         write(vtu1) '<PDataArray '
         write(vtu1) 'type="Int32" '
         write(vtu1) 'Name="connectivity" '
         write(vtu1) 'format="ascii"> '
         write(vtu1) '</PDataArray> '
         write(vtu1) '<PDataArray '
         write(vtu1) 'type="Int32" '
         write(vtu1) 'Name="offsets" '
         write(vtu1) 'format="ascii"> '
         write(vtu1) '</PDataArray> '
         write(vtu1) '<PDataArray '
         write(vtu1) 'type="Int32" '
         write(vtu1) 'Name="types" '
         write(vtu1) 'format="ascii"> '
         write(vtu1) '</PDataArray> '
         write(vtu1) '</PCells> '
         
         do j=0,npmax-1
            write(vtufile,'(A4,I5.5,A2,I5.5,A4)')
     >              'part',icalld,'_p',j,'.vtu'
            write(vtu1) '<Piece Source="',trim(vtufile),'"/>'
         enddo

         write(vtu1) '</PUnstructuredGrid> '
         write(vtu1) '</VTKFile> '

         close(vtu1)

      endif


! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      if (n .gt. 0) then
      write(vtufile,'(A4,I5.5,A2,I5.5,A4)')'part',icalld,'_p',nid,'.vtu'


      vtu=867+nid


      open(unit=vtu,file=vtufile,status='replace',
     >             access='stream',form="unformatted")

! ------------
! FRONT MATTER
! ------------
      write(vtu) '<VTKFile '
      write(vtu) 'type="UnstructuredGrid" '
      write(vtu) 'version="1.0" '
      if (lpm_endian .eq. 0) then
         write(vtu) 'byte_order="LittleEndian"> '
      elseif (lpm_endian .eq. 1) then
         write(vtu) 'byte_order="BigEndian"> '
      endif

      write(vtu) '<UnstructuredGrid> '

      write(vtu) '<FieldData> ' 
      write(vtu) '<DataArray '  ! time
      write(vtu) 'type="Float32" '
      write(vtu) 'Name="TIME" '
      write(vtu) 'NumberOfTuples="1" '
      write(vtu) 'format="ascii"> '
#ifdef CMTNEK
      write(dumstr,'(ES30.16)') time_cmt
#else
      write(dumstr,'(ES30.16)') time
#endif
      write(vtu) trim(dumstr), ' '
      write(vtu) '</DataArray> '

      write(vtu) '<DataArray '  ! cycle
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="CYCLE" '
      write(vtu) 'NumberOfTuples="1" '
      write(vtu) 'format="ascii"> '
      write(dumstr,'(I0)') istep
      write(vtu) trim(dumstr), ' '
      write(vtu) '</DataArray> '
      write(vtu) '</FieldData>'

      write(vtu) '<Piece '
      write(dumstr,'(A16,I0,A2)') 'NumberOfPoints="',n,'" '
      write(vtu) trim(dumstr), ' '
      write(vtu) 'NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      write(vtu) '<Points> '
      iint = 0
      call vtu_write_dataarray(vtu,"Position    ",3,iint,1)
      write(vtu) '</Points> '

! ----
! DATA 
! ----
      idum = 0
      write(vtu) '<PointData> '
      iint = iint + 3*wdsize*n + isize
      call vtu_write_dataarray(vtu,"VelocityP   ",3,iint,1)
      iint = iint + 3*wdsize*n + isize
      call vtu_write_dataarray(vtu,"VelocityF   ",3,iint,1)
      iint = iint + 3*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ForceUSR    ",3,iint,1)
      if (part_force(1) .ne. 0) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceQS     ",3,iint,1)
      endif
      if (part_force(2) .ne. 0) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceUN     ",3,iint,1)
      endif
      if (part_force(3) .ne. 0) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceIU     ",3,iint,1)
      endif
      if (two_way .gt. 2) then
         iint = iint + 3*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceC      ",3,iint,1)
      endif
      idum = 3
      if (part_force(1) .ne. 0 .or. part_force(3) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ForceWork   ",1,iint,1)
         idum = 1
      endif
      if (part_force(4) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"HeatQS      ",1,iint,1)
         idum = 1
      endif
      if (part_force(5) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"HeatUN      ",1,iint,1)
         idum = 1
      endif
      if (part_force(4) .ne. 0 .and. part_force(5) .ne. 0) then
         iint = iint + idum*wdsize*n + isize
         call vtu_write_dataarray(vtu,"TemperatureP",1,iint,1)
         idum = 1
      endif
      iint = iint + idum*wdsize*n + isize
      call vtu_write_dataarray(vtu,"TemperatureF",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"DensityF    ",1,iint,1)
      if (two_way .ge. 2) then
         iint = iint + 1*wdsize*n + isize
         call vtu_write_dataarray(vtu,"ParticleVF  ",1,iint,1)
      endif

      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"RadiusCG    ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"Diameter    ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ID_1        ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ID_2        ",1,iint,1)
      iint = iint + 1*wdsize*n + isize
      call vtu_write_dataarray(vtu,"ID_3        ",1,iint,1)
      write(vtu) '</PointData> '

! ----------
! END MATTER
! ----------
      write(vtu) '<Cells> '
      write(vtu) '<DataArray '
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="connectivity" '
      write(vtu) 'format="ascii"> '
      write(vtu) '</DataArray> '
      write(vtu) '<DataArray '
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="offsets" '
      write(vtu) 'format="ascii"> '
      write(vtu) '</DataArray> '
      write(vtu) '<DataArray '
      write(vtu) 'type="Int32" '
      write(vtu) 'Name="types" '
      write(vtu) 'format="ascii"> '
      write(vtu) '</DataArray> '
      write(vtu) '</Cells> '
      write(vtu) '</Piece> '
      write(vtu) '</UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu) '<AppendedData encoding="raw"> _'

      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jx,i)
         write(vtu) rpart(jy,i)
         write(vtu) rpart(jz,i)
      enddo
      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jv0,i)
         write(vtu) rpart(jv0+1,i)
         write(vtu) rpart(jv0+2,i)
      enddo
      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(ju0,i)
         write(vtu) rpart(ju0+1,i)
         write(vtu) rpart(ju0+2,i)
      enddo
      iint=3*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jfusr,i)
         write(vtu) rpart(jfusr+1,i)
         write(vtu) rpart(jfusr+2,i)
      enddo
      if (part_force(1) .ne. 0) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfqs,i)
            write(vtu) rpart(jfqs+1,i)
            write(vtu) rpart(jfqs+2,i)
         enddo
      endif
      if (part_force(2) .ne. 0) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfun,i)
            write(vtu) rpart(jfun+1,i)
            write(vtu) rpart(jfun+2,i)
         enddo
      endif
      if (part_force(3) .ne. 0) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfiu,i)
            write(vtu) rpart(jfiu+1,i)
            write(vtu) rpart(jfiu+2,i)
         enddo
      endif
      if (two_way .gt. 2) then
         iint=3*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jfcol,i)
            write(vtu) rpart(jfcol+1,i)
            write(vtu) rpart(jfcol+2,i)
         enddo
      endif
      if (part_force(1) .ne. 0 .or. part_force(3) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jg0,i)
         enddo
      endif
      if (part_force(4) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jqqs,i)
         enddo
      endif
      if (part_force(5) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jquu,i)
         enddo
      endif
      if (part_force(4) .ne. 0 .and. part_force(5) .ne. 0) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jtemp,i)
         enddo
      endif
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jtempf,i)
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jrho,i)
      enddo
      if (two_way .ge. 2) then
         iint=1*wdsize*n
         write(vtu) iint
         do i=1,n
            write(vtu) rpart(jvol1,i)
         enddo
      endif
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jrpe,i)
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         write(vtu) rpart(jdp,i)
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         rdum = real(int(ipart(jpid1,i)))
         write(vtu) rdum
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         rdum = real(int(ipart(jpid2,i)))
         write(vtu) rdum
      enddo
      iint=1*wdsize*n
      write(vtu) iint
      do i=1,n
         rdum = real(int(ipart(jpid3,i)))
         write(vtu) rdum
      enddo

      write(vtu) ' </AppendedData> '
      write(vtu) '</VTKFile> '

      close(vtu)
      endif

! ---------------------------------------------------
! SEND PARTICLES BACK TO ORIGINAL PROCESSOR WHEN DONE
! ---------------------------------------------------
      call move_particles_inproc

      return
      end
c----------------------------------------------------------------------
      subroutine vtu_write_dataarray(vtu,dataname,ncomp,idist,ip)
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      integer vtu,ncomp
      integer*4 idist
      character*12 dataname
      character*50 dumstr

      if (ip .eq. 0) then
         write(vtu) '<PDataArray '
         write(vtu) 'type="Float64" '
         write(dumstr,*)'Name="',trim(dataname),'" '
         write(vtu) trim(dumstr), ' '
         write(dumstr,'(A20,I0,A2)')'NumberOfComponents="',ncomp,'" '
         write(vtu) trim(dumstr), ' '
         write(vtu) 'format="append"/> '
      else
         write(vtu) '<DataArray '
         write(vtu) 'type="Float64" '
         write(dumstr,*)'Name="',trim(dataname),'" '
         write(vtu) trim(dumstr), ' '
         write(dumstr,'(A20,I0,A2)')'NumberOfComponents="',ncomp,'" '
         write(vtu) trim(dumstr), ' '
         write(vtu) 'format="append" '
         write(dumstr,'(A8,I0,A3)') 'offset="',idist,'"/>'
         write(vtu) trim(dumstr), ' '
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine output_parallel_restart_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /0/

      character*10 filename

      logical partl         ! This is a dummy placeholder, used in cr()

      integer vtu,vtu1,prevs(3,mp)
      integer*4 iint
      integer stride_len

      integer*8 disp
      integer*4 nptot
      integer pth
      real*4 rnptot, rnptot_new

! ----------------------------------------
! Setup file names to write to mpi
! ----------------------------------------
      icalld = icalld+1
      write(filename,'(A5,I5.5)') 'rpart', icalld
      pth = 174

      if (nid.eq.0) then
         open(unit=pth,iostat=istat,file=filename,status='old')
         if (stat .eq. 0) close(pth, status='delete')
      endif
      nptot = iglsum(n,1)
      rnptot = real(nptot)

      disp   = 0
      icount = 0
      if (nid .eq. 0) icount = 1
      iorank = -1
      call byte_open_mpi(filename,pth,.false.,ierr)
      call byte_set_view(disp,pth)
      call byte_write_mpi(rnptot,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,mp
         prevs(1,i) = i-1
         prevs(2,i) = nid
         prevs(3,i) = n
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 2 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum  = 3
      call fgslib_crystal_ituple_transfer(i_cr_hndl,prevs,
     >                 ndum,mp,mp,nps)
      call fgslib_crystal_ituple_sort(i_cr_hndl,prevs,
     >                 ndum,mp,nglob,nkey)

      stride_len = 0
      if (nid .ne. 0) then
      do i=1,nid
         stride_len = stride_len + prevs(3,i)
      enddo
      endif

! -----------------------------------------------------
! SEND PARTICLES TO A SMALLER NUMBER OF RANKS TO OUTPUT
! -----------------------------------------------------
      do i = 0,n-1
         idum = int((stride_len + i)/llpart)
         ipart(jps,i+1) = idum
      enddo
      nl = 0
      call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     >                  , ipart,ni,partl,nl,rpart,nr,jps)

! ----------------------------------------
! Calculate how many processors to be used
! ----------------------------------------
      npmax = min(nptot/llpart+1,mp)
      stride_len = 0
      if (nid .le. npmax-1 .and. nid. ne. 0) stride_len = nid*llpart

! ----------------------
! COPY DATA FOR EASY I/O
! ----------------------
      do i = 1,n
         rfpts(1,i) = rpart(jx,i)
         rfpts(2,i) = rpart(jy,i)
         rfpts(3,i) = rpart(jz,i)
         rfpts(4,i) = rpart(jx1+0,i)
         rfpts(5,i) = rpart(jx1+1,i)
         rfpts(6,i) = rpart(jx1+2,i)
         rfpts(7,i) = rpart(jx2+0,i)
         rfpts(8,i) = rpart(jx2+1,i)
         rfpts(9,i) = rpart(jx2+2,i)
         rfpts(10,i) = rpart(jx3+0,i)
         rfpts(11,i) = rpart(jx3+1,i)
         rfpts(12,i) = rpart(jx3+2,i)
      
         rfpts(13,i) = rpart(jv0,i)
         rfpts(14,i) = rpart(jv0+1,i)
         rfpts(15,i) = rpart(jv0+2,i)
         rfpts(16,i) = rpart(jv1+0,i)
         rfpts(17,i) = rpart(jv1+1,i)
         rfpts(18,i) = rpart(jv1+2,i)
         rfpts(19,i) = rpart(jv2+0,i)
         rfpts(20,i) = rpart(jv2+1,i)
         rfpts(21,i) = rpart(jv2+2,i)
         rfpts(22,i) = rpart(jv3+0,i)
         rfpts(23,i) = rpart(jv3+1,i)
         rfpts(24,i) = rpart(jv3+2,i)
      
         rfpts(25,i) = rpart(ju0,i)
         rfpts(26,i) = rpart(ju0+1,i)
         rfpts(27,i) = rpart(ju0+2,i)
         rfpts(28,i) = rpart(ju1+0,i)
         rfpts(29,i) = rpart(ju1+1,i)
         rfpts(30,i) = rpart(ju1+2,i)
         rfpts(31,i) = rpart(ju2+0,i)
         rfpts(32,i) = rpart(ju2+1,i)
         rfpts(33,i) = rpart(ju2+2,i)
         rfpts(34,i) = rpart(ju3+0,i)
         rfpts(35,i) = rpart(ju3+1,i)
         rfpts(36,i) = rpart(ju3+2,i)
      
         rfpts(37,i) = rpart(jdp,i)
         rfpts(38,i) = rpart(jspl,i)
         rfpts(39,i) = rpart(jtemp,i)
         rfpts(40,i) = real(ipart(jpid1,i))
         rfpts(41,i) = real(ipart(jpid2,i))
         rfpts(42,i) = real(ipart(jpid3,i))

         rfpts(43,i) = rpart(jangle0+0,i)
         rfpts(44,i) = rpart(jangle0+1,i)
         rfpts(45,i) = rpart(jangle0+2,i)
         rfpts(46,i) = rpart(jangle1+0,i)
         rfpts(47,i) = rpart(jangle1+1,i)
         rfpts(48,i) = rpart(jangle1+2,i)
         rfpts(49,i) = rpart(jangle2+0,i)
         rfpts(50,i) = rpart(jangle2+1,i)
         rfpts(51,i) = rpart(jangle2+2,i)
         rfpts(52,i) = rpart(jangle3+0,i)
         rfpts(53,i) = rpart(jangle3+1,i)
         rfpts(54,i) = rpart(jangle3+2,i)

         rfpts(55,i) = rpart(jangvel0+0,i)
         rfpts(56,i) = rpart(jangvel0+1,i)
         rfpts(57,i) = rpart(jangvel0+2,i)
         rfpts(58,i) = rpart(jangvel1+0,i)
         rfpts(59,i) = rpart(jangvel1+1,i)
         rfpts(60,i) = rpart(jangvel1+2,i)
         rfpts(61,i) = rpart(jangvel2+0,i)
         rfpts(62,i) = rpart(jangvel2+1,i)
         rfpts(63,i) = rpart(jangvel2+2,i)
         rfpts(64,i) = rpart(jangvel3+0,i)
         rfpts(65,i) = rpart(jangvel3+1,i)
         rfpts(66,i) = rpart(jangvel3+2,i)

         rfpts(67,i) = rpart(jtorque0, i)
         rfpts(68,i) = rpart(jtq_col,  i)

         rfpts(69,i) = real(ipart(jgp_role,i))
         rfpts(70,i) = real(ipart(jgp_queen,i))
         rfpts(71,i) = real(ipart(jgp_worker1,i))
         rfpts(72,i) = real(ipart(jgp_nlm,i))
      
      enddo

! -------------------------------
! WRITE PARTICLES TO RESTART FILE
! -------------------------------
      disp   = isize + stride_len*lrf*isize ! is there real*4 var?
      icount = n*lrf
      iorank = -1
      call byte_open_mpi(filename,pth,.false.,ierr)
      call byte_set_view(disp,pth)
      call byte_write_mpi(rfpts,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

! ----------------------------------------
! Move particles back to where they belong
! ----------------------------------------
      call move_particles_inproc


      return
      end
c----------------------------------------------------------------------
      subroutine read_parallel_restart_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      common /myparth/ i_fp_hndl, i_cr_hndl
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*10 filename

      logical partl         ! This is a dummy placeholder, used in cr()

      integer vtu,vtu1,prevs(2,mp)
      integer*4 iint
      integer stride_len

      integer*8 disp
      integer*4 nptot
      integer pth, nread
      real*4 rnptot

      rpi    = 4.0d+0*atan(1.d+0) ! pi

! ------------------------------
! LOAD TOTAL NUMBER OF PARTICLES
! ------------------------------
      write(filename,'(A5,I5.5)') 'rpart', ipart_restartr

      pth = 185

      disp = 0
      icount = 0
      if (nid .eq. 0) icount = 1
      iorank = -1
      call byte_open_mpi(filename,pth,.true.,ierr)
      call byte_set_view(disp,pth)
      call byte_read_mpi(rnptot,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

      nptot = int(rnptot)
      call bcast(nptot, isize)

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      npmax = min(nptot/llpart+1,mp)
      stride_len = 0
      if (nid .le. npmax-1 .and. nid. ne. 0) stride_len = nid*llpart

      nread = llpart
      if (nid .gt. npmax-1) nread = 0

      ndiff = nptot - (npmax-1)*llpart
      if (nid .eq. npmax-1) nread = ndiff

! -------------------------
! Parallel MPI file read in
! -------------------------
      disp   = isize + stride_len*lrf*isize ! is there real*4 var?
      icount = nread*lrf
      iorank = -1
      call byte_open_mpi(filename,pth,.true.,ierr)
      call byte_set_view(disp,pth)
      call byte_read_mpi(rfpts,icount,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

! ----------------------------------------
! Assign values read in to rpart and ipart
! ----------------------------------------
      i = n ! if there are previous particles
      do ii = 1,nread
         i = n + ii
         rpart(jx,i)    =  rfpts(1,ii) 
         rpart(jy,i)    =  rfpts(2,ii)
         rpart(jz,i)    =  rfpts(3,ii)  
         rpart(jx1+0,i) =  rfpts(4,ii)  
         rpart(jx1+1,i) =  rfpts(5,ii)
         rpart(jx1+2,i) =  rfpts(6,ii)  
         rpart(jx2+0,i) =  rfpts(7,ii)  
         rpart(jx2+1,i) =  rfpts(8,ii)
         rpart(jx2+2,i) =  rfpts(9,ii)  
         rpart(jx3+0,i) =  rfpts(10,ii)
         rpart(jx3+1,i) =  rfpts(11,ii)
         rpart(jx3+2,i) =  rfpts(12,ii)
                                        
         rpart(jv0,i)   =  rfpts(13,ii)
         rpart(jv0+1,i) =  rfpts(14,ii)
         rpart(jv0+2,i) =  rfpts(15,ii)
         rpart(jv1+0,i) =  rfpts(16,ii)
         rpart(jv1+1,i) =  rfpts(17,ii)
         rpart(jv1+2,i) =  rfpts(18,ii)
         rpart(jv2+0,i) =  rfpts(19,ii)
         rpart(jv2+1,i) =  rfpts(20,ii)
         rpart(jv2+2,i) =  rfpts(21,ii)
         rpart(jv3+0,i) =  rfpts(22,ii)
         rpart(jv3+1,i) =  rfpts(23,ii)
         rpart(jv3+2,i) =  rfpts(24,ii)
                                        
         rpart(ju0,i)   =  rfpts(25,ii)
         rpart(ju0+1,i) =  rfpts(26,ii)
         rpart(ju0+2,i) =  rfpts(27,ii)
         rpart(ju1+0,i) =  rfpts(28,ii)
         rpart(ju1+1,i) =  rfpts(29,ii)
         rpart(ju1+2,i) =  rfpts(30,ii)
         rpart(ju2+0,i) =  rfpts(31,ii)
         rpart(ju2+1,i) =  rfpts(32,ii)
         rpart(ju2+2,i) =  rfpts(33,ii)
         rpart(ju3+0,i) =  rfpts(34,ii)
         rpart(ju3+1,i) =  rfpts(35,ii)
         rpart(ju3+2,i) =  rfpts(36,ii)
                                        
         rpart(jdp,i)   =  rfpts(37,ii)
         rpart(jspl,i)  =  rfpts(38,ii)
         rpart(jtemp,i) =  rfpts(39,ii)
         ipart(jpid1,i) =  nint(rfpts(40,ii))
         ipart(jpid2,i) =  nint(rfpts(41,ii))
         ipart(jpid3,i) =  nint(rfpts(42,ii))

         ! extra stuff
         rpart(jtaup,i) = rpart(jdp,i)**2*rho_p/18.0d+0/mu_0
         rpart(jrhop,i) = rho_p      ! material density of particle
         rpart(jvol,i)  = rpart(jspl,i)*rpi*rpart(jdp,i)**3/6.d+0! particle volume
         rpart(jgam,i)  = 1.          ! initial integration correction
         rpart(jrpe,i)  = rpart(jspl,i)**(1./3.)*rpart(jdp,i)/2.
         
         rpart(jangle0+0,i) = rfpts(43,ii) 
         rpart(jangle0+1,i) = rfpts(44,ii) 
         rpart(jangle0+2,i) = rfpts(45,ii)
         rpart(jangle1+0,i) = rfpts(46,ii)
         rpart(jangle1+1,i) = rfpts(47,ii)
         rpart(jangle1+2,i) = rfpts(48,ii)
         rpart(jangle2+0,i) = rfpts(49,ii)
         rpart(jangle2+1,i) = rfpts(50,ii)
         rpart(jangle2+2,i) = rfpts(51,ii)
         rpart(jangle3+0,i) = rfpts(52,ii)
         rpart(jangle3+1,i) = rfpts(53,ii)
         rpart(jangle3+2,i) = rfpts(54,ii)

         rpart(jangvel0+0,i) = rfpts(55,ii)
         rpart(jangvel0+1,i) = rfpts(56,ii)
         rpart(jangvel0+2,i) = rfpts(57,ii)
         rpart(jangvel1+0,i) = rfpts(58,ii)
         rpart(jangvel1+1,i) = rfpts(59,ii)
         rpart(jangvel1+2,i) = rfpts(60,ii)
         rpart(jangvel2+0,i) = rfpts(61,ii)
         rpart(jangvel2+1,i) = rfpts(62,ii)
         rpart(jangvel2+2,i) = rfpts(63,ii)
         rpart(jangvel3+0,i) = rfpts(64,ii)
         rpart(jangvel3+1,i) = rfpts(65,ii)
         rpart(jangvel3+2,i) = rfpts(66,ii)

         rpart(jtorque0, i)  = rfpts(67,ii)
         rpart(jtq_col,  i)  = rfpts(68,ii)

         ipart(jgp_role,i)    = nint( rfpts(69,ii) ) 
         ipart(jgp_queen,i)   = nint( rfpts(70,ii) ) 
         ipart(jgp_worker1,i) = nint( rfpts(71,ii) ) 
         ipart(jgp_nlm,i)     = nint( rfpts(72,ii) )

         call get_part_use_block(i,0)
         call lpm_usr_f  ! can overwrite particle properties at init
         call get_part_use_block(i,1)
      enddo
      n = i

      call update_particle_location   
      call move_particles_inproc

      return
      end
c----------------------------------------------------------------------
      subroutine output_two_way_io
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'LPM'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real rfpfluid(3),rfpfluidl(3),msum,msum_tot(3,2)

      
      ! show element map
       do ie=1,nelt
       do iz=1,nz1
       do iy=1,ny1
       do ix=1,nx1
          ptw(ix,iy,iz,ie,5) = real(nid)
       enddo
       enddo
       enddo
       enddo

      
      itmp = 1
      call outpost2(ptw(1,1,1,1,1),         ! fhyd_x
     >              ptw(1,1,1,1,2),         ! fhyd_y
     >              ptw(1,1,1,1,3),         ! fhyd_z
     >              ptw(1,1,1,1,4),         ! phi_p (but not if lx1!=lx2
     >              ptw(1,1,1,1,5),         ! phi_p
     >              itmp          ,        
     >              'ptw')


c     eulerian integrations -----------------------------------------
c     fluid momentum 
      msum_tot(1,1) = glsc3(bm1,vtrans,vx,nx1*ny1*nz1*nelv)
      msum_tot(2,1) = glsc3(bm1,vtrans,vy,nx1*ny1*nz1*nelv)
      msum_tot(3,1) = glsc3(bm1,vtrans,vz,nx1*ny1*nz1*nelv)
c     particle volume fraction
      vf_part_e     = glsc2(bm1,ptw(1,1,1,1,4),nx1*ny1*nz1*nelt)
                                                 ! in z of mono-particle
                                                 ! Dp
c     particle forces on fluid
      rfpfluid(1)   = glsc2(bm1,ptw(1,1,1,1,1),nx1*ny1*nz1*nelt)
      rfpfluid(2)   = glsc2(bm1,ptw(1,1,1,1,2),nx1*ny1*nz1*nelt)
      rfpfluid(3)   = glsc2(bm1,ptw(1,1,1,1,3),nx1*ny1*nz1*nelt)

      if (.not.if3d) vf_part_e   = vf_part_e*dp(1)   ! Here:
      if (.not.if3d) rfpfluid(1) = rfpfluid(1)*dp(1) ! for 2d, assume
      if (.not.if3d) rfpfluid(2) = rfpfluid(2)*dp(1) ! z thicknes of 
      if (.not.if3d) rfpfluid(3) = rfpfluid(3)*dp(1) ! monodisperse Dp


c     lagrangian integrations ---------------------------------------
c     particle momentum
      do ieq=0,2
         msum = 0.0
         rsum = 0.0
         do i=1,n
           msum = msum + 
     >       rpart(jv0+ieq,i)*rpart(jrhop,i)*rpart(jvol,i)
           rsum = rsum + rpart(jf0+ieq,i)
        enddo
         msum_tot(ieq+1,2) = glsum(msum,1)
         rfpfluidl(1+ieq)  = glsum(rsum,1)
      enddo
c     particle volume fraction
      msum = 0.0
      do i=1,n
         msum = msum + rpart(jvol,i)
      enddo
      vf_part_l = glsum(msum,1)

      vf_rel_error = abs(vf_part_l - vf_part_e)/vf_part_l*100.0

c     print to files ------------------------------------------------
c     print properties to logfile
      if (nid.eq.0) write(6,500) "E Fluid Momentum :              ", 
     >                  istep, msum_tot(1,1),msum_tot(2,1),msum_tot(3,1)
      if (nid.eq.0) write(6,500) "E Particle forces:              ", 
     >                  istep, rfpfluid(1),rfpfluid(2),rfpfluid(3)
      if (nid.eq.0) write(6,500) "E Particle Volume:              ", 
     >                  istep, vf_part_e
      if (nid.eq.0) write(6,500) "L Particle Momentum :           ", 
     >                  istep, msum_tot(1,2),msum_tot(2,2),msum_tot(3,2)
      if (nid.eq.0) write(6,500) "L Particle forces:              ", 
     >                  istep, rfpfluidl(1),rfpfluidl(2),rfpfluidl(3)
      if (nid.eq.0) write(6,500) "L Particle Volume:              ", 
     >                  istep, vf_part_l
      if (nid.eq.0) write(6,500) "VF Relative Error %              ", 
     >                  istep, vf_rel_error

  500 FORMAT(A30,I20,3ES20.10)

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_input_defaults
      include 'SIZE'
      include 'LPM'

      ! set some defaults
      nw = 0
      rxbo(1,1) = -1E8
      rxbo(2,1) = -1E8
      rxbo(1,2) = -1E8
      rxbo(2,2) = -1E8
      rxbo(1,3) = -1E8
      rxbo(2,3) = -1E8
      rxco(1) = -1E8
      rxco(2) = -1E8
      rxco(3) = -1E8
      rxco(4) = -1E8
      rxco(5) = -1E8
      rxco(6) = -1E8
      rxco(7) = -1E8
      rxco(8) = -1E8
      rxco(9) = -1E8
      dp(1) = 0.
      dp(2) = 0.
      dp_std = -1.
      tp_0  = 273.
      rho_p = 2500.
      cp_p  = 840.
      do i=1,5
         part_force(i) = 0
      enddo
      time_integ  = 1
      two_way     = 1
      red_interp  = 1
      npio_method = 1
      inject_rate = 0
      time_delay  = 0
      nrandseed   = 1
      lpm_endian  = 0
      npro_method = 1
      rspl        = 1.
      dfilt       = 2.
      ralphdecay  = 1E-3
      do i=1,6
         bc_part(i) = 1
      enddo
      ipart_restartr = 0
      ksp            = 10.
      e_rest         = 0.9
      ksp_wall       = ksp
      e_rest_wall    = e_rest

      np_walls = 0
      nc_walls = 0

      return
      end
c----------------------------------------------------------------------
c     subroutine read_particle_input_par
c     include 'SIZE'
c     include 'INPUT'
c     include 'TSTEP'
c     include 'PARALLEL'
c     include 'LPM'
c
c     only used for load balanced code since there is no par file in 
c     old version of the code. To convert from standard nek5000 input:
c         periodic* = only needs to be included for it to be periodic
c         so yes or no doesn't matter
c
c         remove commas seperating parameter inputs
c
c         must be lowercase key words with one space then an equal sign
c         followed by another space

c     character*72 dum_str

c     character*200 buffer, label
c     integer pos, fh, ios, line, dum
c     parameter(fh = 15)

c     call lpm_input_defaults

c     ios  = 0
c     line = 0

c     open(fh, file='particles.par')
c
c     do while (ios == 0)
c        read(fh, '(A)', iostat=ios) buffer
c        if (ios == 0) then
c           line = line + 1
c     
c           ! Find the first instance of whitespace.  Split label and data.
c           pos = scan(buffer, '=')
c           label = buffer(1:pos)
c           buffer = buffer(pos+1:)
c     
c           select case (label)
c           case ('npart =')
c              read(buffer, *, iostat=ios) nw
c              if(nid.eq.0)write(6,*) 'Read npart: ', nw
c           case ('distributebox =')
c              read(buffer, *, iostat=ios) rxbo(1,1),rxbo(2,1),
c    >                                     rxbo(1,2),rxbo(2,2),
c    >                                     rxbo(1,3),rxbo(2,3)
c              if(nid.eq.0)write(6,*) 'Read distributebox '
c           case ('distributecylinder =')
c              read(buffer, *, iostat=ios) rxco(1),rxco(2),rxco(3),
c    >                                     rxco(4),rxco(5),rxco(6), 
c    >                                     rxco(7),rxco(8),rxco(9)
c              if(nid.eq.0)write(6,*) 'Read distributecylinder '
c           case ('distributesphere =')
c              read(buffer, *, iostat=ios) rxco(1),rxco(2),rxco(3),
c    >                                     rxco(4),rxco(5)
c              if(nid.eq.0)write(6,*) 'Read distributesphere '
c           case ('diameter =')sp
c              read(buffer, *, iostat=ios) dp(1)
c              if(nid.eq.0)write(6,*) 'Read diameter: ', dp(1)
c              dp(2) = dp(1)
c           case ('diameteruniform =')
c              read(buffer, *, iostat=ios) dp(1), dp(2)
c              if(nid.eq.0)write(6,*) 'Read diameteruniform: ', 
c    >                         dp(1), dp(2)
c           case ('diametergaussian =')
c              read(buffer, *, iostat=ios) dp(1), dp_std
c              if(nid.eq.0)write(6,*) 'Read diametergaussian: ', 
c    >                         dp(1), dp_std
c              dp(2) = dp(1)
c           case ('temperature =')
c              read(buffer, *, iostat=ios) tp_0
c              if(nid.eq.0)write(6,*) 'Read temperature: ', tp_0
c           case ('density =')
c              read(buffer, *, iostat=ios) rho_p
c              if(nid.eq.0)write(6,*) 'Read density: ', rho_p
c           case ('specificheat =')
c              read(buffer, *, iostat=ios) cp_p
c              if(nid.eq.0)write(6,*) 'Read specificheat: ', cp_p
c           case ('forceqs =')
c              read(buffer, *, iostat=ios) part_force(1)
c              if(nid.eq.0)write(6,*) 'Read forceqs: ', part_force(1)
c           case ('forceun =')
c              read(buffer, *, iostat=ios) part_force(2)
c              if(nid.eq.0)write(6,*) 'Read forceun: ', part_force(2)
c           case ('forceiu =')
c              read(buffer, *, iostat=ios) part_force(3)
c              if(nid.eq.0)write(6,*) 'Read forceiu: ', part_force(3)
c           case ('heatqs =')
c              read(buffer, *, iostat=ios) part_force(4)
c              if(nid.eq.0)write(6,*) 'Read heatqs: ', part_force(4)
c           case ('heatun =')
c              read(buffer, *, iostat=ios) part_force(5)
c              if(nid.eq.0)write(6,*) 'Read heatun: ', part_force(5)
c           case ('timestepper =') 
c              read(buffer, *, iostat=ios) time_integ
c              if(nid.eq.0)write(6,*) 'Read timestepper: ', time_integ
c           case ('coupling =') 
c              read(buffer, *, iostat=ios) two_way
c              if(nid.eq.0)write(6,*) 'Read coupling: ', two_way
c           case ('interpolation') 
c              read(buffer, *, iostat=ios) red_interp
c              if(nid.eq.0)write(6,*) 'Read interpolation: ', red_interp
c           case ('io =')
c              read(buffer, *, iostat=ios) npio_method
c              if(nid.eq.0)write(6,*) 'Read io: ', npio_method
c           case ('injectionstep =')
c              read(buffer, *, iostat=ios) inject_rate
c              if(nid.eq.0)write(6,*) 'Read injectionstep: ',inject_rate
c           case ('delaystep =')
c              read(buffer, *, iostat=ios) time_delay
c              if(nid.eq.0)write(6,*) 'Read delaystep: ', time_delay
c           case ('seed =')
c              read(buffer, *, iostat=ios) nrandseed
c              if(nid.eq.0)write(6,*) 'Read seed: ', nrandseed
c           case ('projection =')
c              read(buffer, *, iostat=ios) npro_method
c              if(nid.eq.0)write(6,*) 'Read projection: ', npro_method
c           case ('coarsegrain =')
c              read(buffer, *, iostat=ios) rspl
c              if(nid.eq.0)write(6,*) 'Read coarsegrain: ', rspl
c           case ('filter =')
c              read(buffer, *, iostat=ios) dfilt
c              if(nid.eq.0)write(6,*) 'Read filter: ', dfilt
c           case ('alpha =')
c              read(buffer, *, iostat=ios) ralphdecay
c              if(nid.eq.0)write(6,*) 'Read alpha: ', ralphdecay

c           case ('wallp01 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp02 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp03 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp04 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp05 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp06 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp07 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp08 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallp09 =') 
c               call add_wall_list(np_walls,buffer,plane_wall_coords,6)
c           case ('wallc01 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc02 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc03 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc04 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc05 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc06 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc07 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc08 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('wallc09 =') 
c               call add_wall_list(nc_walls,buffer,cyl_wall_coords,7)
c           case ('periodicx =')
c              read(buffer, *, iostat=ios)
c              if(nid.eq.0)write(6,*) 'Read periodicx '
c              bc_part(1) = 0
c              bc_part(2) = 0
c           case ('periodicy =')
c              read(buffer, *, iostat=ios)
c              if(nid.eq.0)write(6,*) 'Read periodicy '
c              bc_part(3) = 0
c              bc_part(4) = 0
c           case ('periodicz =')
c              read(buffer, *, iostat=ios)
c              if(nid.eq.0)write(6,*) 'Read periodicz '
c              bc_part(5) = 0
c              bc_part(6) = 0
c           case ('restartstep =')
c              read(buffer, *, iostat=ios) ipart_restartr
c              if(nid.eq.0)write(6,*)'Read restartstep: ',ipart_restartr
c           case ('spring =')
c              read(buffer, *, iostat=ios) ksp
c              if(nid.eq.0)write(6,*) 'Read spring: ', ksp
c           case ('restitution =')
c              read(buffer, *, iostat=ios) e_rest
c              if(nid.eq.0)write(6,*) 'Read restitution: ', e_rest
c           case default
c              if(nid.eq.0)write(6,*) 'Skipping label at line', line
c           end select
c        end if
c     enddo

c     close(fh)


c     return
c     end
c----------------------------------------------------------------------
c     subroutine add_wall_list(ninc,buf,list,ncomp)
c     include 'SIZE'
c     include 'TOTAL'
c     include 'LPM'

c     integer ninc,ncom
c     real list(ncomp,n_walls)
c     character*200 buf

c     ninc = ninc + 1
c     if (ninc .gt. n_walls) then
c        if (nid.eq.0) then
c            write(6,*) 'Increase max number particle wall',ncomp
c            call exitt
c         endif
c     endif
c     if (ncomp .eq. 7) then
c        read(buf, *, iostat=ios) cyl_wall_coords(1,nc_walls)
c    >                           ,cyl_wall_coords(2,nc_walls)
c    >                           ,cyl_wall_coords(3,nc_walls)
c    >                           ,cyl_wall_coords(4,nc_walls)
c    >                           ,cyl_wall_coords(5,nc_walls)
c    >                           ,cyl_wall_coords(6,nc_walls)
c    >                           ,cyl_wall_coords(7,nc_walls)
c        if(nid.eq.0)write(6,*) 'Read wall_cyl number ', ninc
c     elseif(ncomp .eq. 6) then
c        read(buf, *, iostat=ios) plane_wall_coords(1,np_walls)
c    >                           ,plane_wall_coords(2,np_walls)
c    >                           ,plane_wall_coords(3,np_walls)
c    >                           ,plane_wall_coords(4,np_walls)
c    >                           ,plane_wall_coords(5,np_walls)
c    >                           ,plane_wall_coords(6,np_walls)
c         if(nid.eq.0)write(6,*) 'Read wall_plane number ',ninc
c      endif

c     return
c     end
c----------------------------------------------------------------------
      subroutine output_particle_diagnostics
      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      integer icalld
      save icalld
      data icalld /0/

      real rtpart
      save rtpart

      real rdiags_part(3,2), rvels(4)

      if (icalld .eq. 0) then
         rtpart = time_cmt
         if (icmtp .eq. 0) rtpart = time
         icalld = icalld + 1
      endif
      rtdum  = time_cmt - rtpart
      if (icmtp .eq. 0) rtdum = time - rtpart

      do i=1,3
         rdiags_part(i,1) =  1E8
         rdiags_part(i,2) = -1E8
      enddo
      do i=1,4
         rvels(i) = 0.0
      enddo
      
      do i=1,n
         rvx = rpart(jv0,i)
         rvy = rpart(jv0+1,i)
         rvz = rpart(jv0+2,i)
         rvmag = sqrt(rvx**2 + rvy**2 + rvz**2)

         rxx = rpart(jx,i)
         ryy = rpart(jy,i)
         rzz = rpart(jz,i)
         rrp = rpart(jrpe,i)

         ! maxes
         if (rxx + rrp .gt. rdiags_part(1,2)) rdiags_part(1,2) = rxx+rrp
         if (ryy + rrp .gt. rdiags_part(2,2)) rdiags_part(2,2) = ryy+rrp
         if (rzz + rrp .gt. rdiags_part(3,2)) rdiags_part(3,2) = rzz+rrp

         ! mins
         if (rxx - rrp .lt. rdiags_part(1,1)) rdiags_part(1,1) = rxx-rrp
         if (ryy - rrp .lt. rdiags_part(2,1)) rdiags_part(2,1) = ryy-rrp
         if (rzz - rrp .lt. rdiags_part(3,1)) rdiags_part(3,1) = rzz-rrp

         ! velocities
         if ( abs(rvx)   .gt. abs(rvels(1)) )  rvels(1) = rvx
         if ( abs(rvy)   .gt. abs(rvels(2)) )  rvels(2) = rvy
         if ( abs(rvz)   .gt. abs(rvels(3)) )  rvels(3) = rvz
         if ( abs(rvmag) .gt. abs(rvels(4)) )  rvels(4) = rvmag
      enddo

      ! compute globally now
      do i=1,3
         rdum = rdiags_part(i,1)
         rdiags_part(i,1) =  glmin(rdum,1)
         rdum = rdiags_part(i,2)
         rdiags_part(i,2) =  glmax(rdum,1)
      enddo
      do i=1,4
         rdum  = rvels(i)
         rdum1 = glmin(rdum,1)
         rdum  = rvels(i)
         rdum2 = glmax(rdum,1)
         rvels(i) = rdum1
         if (abs(rdum1) .lt. abs(rdum2)) rvels(i) = rdum2
      enddo

      nptot = iglsum(n,1)

      if (nid .eq. 0) then
      write(6,*)'----- START PARTICLE DIAGNOSTICS: -----'
      write(6,*)'NPART TOTAL      :',nptot
      write(6,*)'XMIN,XMAX        :',rdiags_part(1,1),rdiags_part(1,2)
      write(6,*)'YMIN,YMAX        :',rdiags_part(2,1),rdiags_part(2,2)
      write(6,*)'ZMIN,ZMAX        :',rdiags_part(3,1),rdiags_part(3,2)
      write(6,*)'MAX(VX,VY,VZ,|V|):',rvels(1),rvels(2),rvels(3),rvels(4)
      write(6,*)'PTIME            :',rtdum
      write(6,*)'----- END PARTICLE DIAGNOSTICS: -----'
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine output_particle_timers
      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'LPM'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      rftime_t = 0.
      rptime_t = 0.

      if(nid.eq.0) then
         write(6,*) 'TIMER H: ', istep
      endif

      do i = 1,iptlen
         rdum  = pttime(i)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/mp
         if(nid.eq.0)  write(6,*) 'TIMER #:',i,rtime

         ! fluid and particle total time: note i == iptlen is f_col
         if (i .eq. 1) rftime_t = rtime
         if ((i .gt. 1) .and. (i.ne.iptlen)) rptime_t = rptime_t +
     >                                                  rtime
      enddo

      
      if (nid.eq.0) then
         write (6,*) 'TOTAL F:', rftime_t
         write (6,*) 'TOTAL P:', rptime_t
      endif


      return
      end
c----------------------------------------------------------------------
c     effeciently move particles between processors routines
c----------------------------------------------------------------------
      subroutine move_particles_inproc
c     Interpolate fluid velocity at current xyz points and move
c     data to the processor that owns the points.
c     Input:    n = number of points on this processor
c     Output:   n = number of points on this processor after the move
c     Code checks for n > llpart and will not move data if there
c     is insufficient room.
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LPM'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

c     common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
c     integer gfirst, inoassignd, resetFindpts, pload

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      logical partl         ! This is a dummy placeholder, used in cr()

      nl = 0                ! No logicals exchanged

      if (icalld1.eq.0) then
c     if (icalld1.eq.0 .or. (resetFindpts .eq. 1)) then
         tolin = 1.e-12
         if (wdsize.eq.4) tolin = 1.e-6
         call intpts_setup  (tolin,i_fp_hndl)
         call fgslib_crystal_setup (i_cr_hndl,nekcomm,np)
      endif

      icalld1 = icalld1 + 1

      call fgslib_findpts(i_fp_hndl !  stride     !   call fgslib_findpts( ihndl,
     $        , ipart(jrc,1),li        !   $             rcode,1,
     $        , ipart(jpt,1),li        !   &             proc,1,
     $        , ipart(je0,1),li        !   &             elid,1,
     $        , rpart(jr ,1),lr        !   &             rst,ndim,
     $        , rpart(jd ,1),lr        !   &             dist,1,
     $        , rpart(jx ,1),lr        !   &             pts(    1),1,
     $        , rpart(jy ,1),lr        !   &             pts(  n+1),1,
     $        , rpart(jz ,1),lr ,n)    !   &             pts(2*n+1),1,n
c         if (nid .eq.1) write(*,*) "check 2"
      nmax = iglmax(n,1)
      if (nmax.gt.llpart) then
         if (nid.eq.0) write(6,1) nmax,llpart
    1    format('WARNING: Max number of particles:',
     $   i9,'.  Not moving because llpart =',i9,'.')
      else
c        Move particle info to the processor that owns each particle
c        using crystal router in log P time:

         jps = jpid1-1     ! Pointer to temporary proc id for swapping
         do i=1,n        ! Can't use jpt because it messes up particle info
            ipart(jps,i) = ipart(jpt,i)
         enddo
         call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     $              , ipart,ni,partl,nl,rpart,nr,jps)
c        Sort by element number - for improved local-eval performance
         call fgslib_crystal_tuple_sort    (i_cr_hndl,n 
     $              , ipart,ni,partl,nl,rpart,nr,je0,1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine intpts_setup(tolin,ih)
c
c setup routine for interpolation tool
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c
      include 'SIZE'
      include 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      tol = tolin
      if (tolin.lt.0) tol = 1e-13 ! default tolerance 

      n       = lx1*ly1*lz1*lelt 
      npt_max = 1256   !256 
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.1 ! relative size to expand bounding boxes by
c
      if(nid.eq.0) write(6,*) 'initializing intpts(), tol=', tol
      call fgslib_findpts_setup(ih,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c
      if(nid.eq.0) write(6,*) 'finish intptes'
      return
      end
c----------------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with 
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate 
c between 0.0 and 1.0 (exclusive of the endpoint values). 
c Call with idum a negative integer to initialize; thereafter, do not alter 
c idum between successive deviates in a sequence. RNMX should approximate the 
c largest floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then 
         idum1=max(-idum,1) 
         idum2=idum1
         do j=NTAB+8,1,-1
            k=idum1/IQ1
            idum1=IA1*(idum1-k*IQ1)-k*IR1 
            if (idum1.lt.0) idum1=idum1+IM1 
            if (j.le.NTAB) iv(j)=idum1
         enddo
         iy=iv(1) 
      endif
      k=idum1/IQ1 
      idum1=IA1*(idum1-k*IQ1)-k*IR1
      if (idum1.lt.0) idum1=idum1+IM1 
      k=idum2/IQ2 
      idum2=IA2*(idum2-k*IQ2)-k*IR2 
      if (idum2.lt.0) idum2=idum2+IM2 
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum1 
      if(iy.lt.1)iy=iy+IMM1 
      ran2=min(AM*iy,RNMX)
      return
      END
c----------------------------------------------------------------------
      function unif_random(rxl,rxr)
c
c     must initialize ran2 first
c
      real xl,xr,unif_random

      rdum       = ran2(2)
      rlen       = rxr - rxl
      unif_random= rxl + rdum*rlen

      return
      end
c-----------------------------------------------------------------------
      function unif_random_norm(rmu,rstd)
c
c     must initialize ran2 first
c
      real xl,xr,unif_random_norm,rstd,rxfne(1000),rcdf(1000)

      nxfn  = 1000
      rxlf  = max(0.0,rmu-5.*rstd)
      rxrf  = rmu+5.*rstd
      rdxf  = (rxrf-rxlf)/(nxfn-1.)

      do i=1,nxfn
         rxfne(i) = rxlf + (i-1.)*rdxf
         rcdf(i)  = 0.5*(1. + erf((rxfne(i)-rmu)/(rstd*sqrt(2.))))
      enddo

      rdum = unif_random(0.,1.)

!     find lower min value for inverse sampling
      idum = 0
      rmin = 100.
      do i=1,nxfn
         if (abs(rdum - rcdf(i)) .lt. rmin) then
            rmin = abs(rdum -rcdf(i))
            idum = i
         endif
      enddo
      ml = idum
      if (rdum .lt. rcdf(idum)) ml = ml + 1

      if (rdum .gt. rcdf(nxfn)) then
         unif_random_norm = rxrf 
      elseif (rdum .lt. rcdf(1)) then
         unif_random_norm = rxlf
      else
         rm = (rxfne(ml+1) - rxfne(ml))/(rcdf(ml+1) - rcdf(ml))
         unif_random_norm = rxfne(ml) + rm*(rdum - rcdf(ml))
      endif

      return
      end
c-----------------------------------------------------------------------
      function unif_random_sphere(rxl,rxr)
c
c     must initialize ran2 first
c
      parameter (nxfn = 10000)
      real xl,xr,unif_random_sphere,rxfne(nxfn)

      real    unif_random
      external unif_random

      rxlf  = rxl
      rxrf  = rxr
      rdxf  = (1. - (rxlf/rxrf)**3)/(nxfn-1.)

      do i=1,nxfn
         rxfne(i) = (rxlf/rxrf)**3 + (i-1.)*rdxf
      enddo

      rdum = unif_random((rxl/rxr)**3,rxr/rxr)

      unif_random_sphere = rxr*(rdum)**(1./3.)

      return
      end
c-----------------------------------------------------------------------
      function unif_random_cyl(rxl,rxr)
c
c     must initialize ran2 first
c
      parameter (nxfn = 10000)
      real xl,xr,unif_random_cyl,rxfne(nxfn)

      real    unif_random
      external unif_random

      rxlf  = rxl
      rxrf  = rxr
      rdxf  = (1. - (rxlf/rxrf)**2)/(nxfn-1.)

      do i=1,nxfn
         rxfne(i) = (rxlf/rxrf)**2 + (i-1.)*rdxf
      enddo

      rdum = unif_random((rxl/rxr)**2,rxr/rxr)

      unif_random_cyl = rxr*sqrt(rdum)

      return
      end
c----------------------------------------------------------------------
C> Compute coefficients for Runge-Kutta stages \cite{TVDRK}
      subroutine set_tstep_coef_part(dt_in)

      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      real dt_in

      dt_cmt = dt_in

      tcoef(1,1) = 0.0
      tcoef(2,1) = 1.0 
      tcoef(3,1) = dt_cmt
      tcoef(1,2) = 3.0/4.0
      tcoef(2,2) = 1.0/4.0 
      tcoef(3,2) = dt_cmt/4.0 
      tcoef(1,3) = 1.0/3.0
      tcoef(2,3) = 2.0/3.0 
      tcoef(3,3) = dt_cmt*2.0/3.0 

      return
      end
c----------------------------------------------------------------------
      subroutine compute_phig_qtl(rdt_in,div)
c
c     Computes modified divergence constraint for multiphase dense
c     compressible flow
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real div(lx2,ly2,lz2,lelv)

      common /phig_qtl_blk/ phig_last,phig_qtl
      real phig_last(lx1,ly1,lz1,lelt,4),phig_qtl(lx1,ly1,lz1,lelt)

      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)

      integer icalld
      save    icalld
      data    icalld  /-1/

      icalld = icalld + 1

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         phig_last(ix,iy,iz,ie,1) = 1. - ptw(ix,iy,iz,ie,4)
         phig_qtl(ix,iy,iz,ie) = 0.
      enddo
      enddo
      enddo
      enddo

      if (icalld .eq. 0) then
         do ie=1,nelt
         do iz=1,nz1
         do iy=1,ny1
         do ix=1,nx1
            phig_last(ix,iy,iz,ie,2) = 1. - ptw(ix,iy,iz,ie,4)
            phig_last(ix,iy,iz,ie,3) = 1. - ptw(ix,iy,iz,ie,4)
            phig_last(ix,iy,iz,ie,4) = 1. - ptw(ix,iy,iz,ie,4)
         enddo
         enddo
         enddo
         enddo

         goto 123
      endif
      
      if (icalld .lt. 5) goto 123

      do ie=1,nelt
         call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                   phig_last(1,1,1,ie,1),lx1,if3d)
         
         do iz=1,nz1
         do iy=1,ny1
         do ix=1,nx1
            phig_qtl(ix,iy,iz,ie) = phig_last(ix,iy,iz,ie,1) -
     >                              phig_last(ix,iy,iz,ie,2)
            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie)/rdt_in

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vx(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dx
     >             (ur(ix,iy,iz)*RXM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SXM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TXM1(ix,iy,iz,ie))

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vy(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dy
     >             (ur(ix,iy,iz)*RYM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SYM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TYM1(ix,iy,iz,ie))

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vz(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dz
     >             (ur(ix,iy,iz)*RZM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SZM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TZM1(ix,iy,iz,ie))


            phig_qtl(ix,iy,iz,ie) = -1.0*phig_qtl(ix,iy,iz,ie)/
     >                              phig_last(ix,iy,iz,ie,1)
         enddo
         enddo
         enddo
      enddo

  123 continue

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         phig_last(ix,iy,iz,ie,4) = phig_last(ix,iy,iz,ie,3)
         phig_last(ix,iy,iz,ie,3) = phig_last(ix,iy,iz,ie,2)
         phig_last(ix,iy,iz,ie,2) = phig_last(ix,iy,iz,ie,1)
      enddo
      enddo
      enddo
      enddo

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         div(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie)
      enddo
      enddo
      enddo
      enddo

c     wght = 1.0
c     ncut = 1
c     call filter_s0(div(1,1,1,1),wght,ncut,'div') 

      return
      end
c-----------------------------------------------------------------------
      subroutine pre_sim_collisions
c
c     time stepping routine for pre-simulation collisions/settling
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage
      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      nmax_step = nsteps  ! number of pre-iteration steps
      ninj_step = 3000

      nstage_part = 3
      if (abs(time_integ) .eq. 2) nstage_part = 1

      if(ipart_moving.ne.1) nstage_part = 1

      ! pre simulation iteration for packed bed
      do istep=0,nmax_step

         if (istep.eq.0) then
            time = 0
            pttime(1) = 0.
         else
            pttime(1) = pttime(1) + dnekclock() - ptdum(1)
         endif

         if (nid.eq. 0) write(6,*) 'pre-sim_io time',istep,time,dt_cmt
         if(mod(istep,iostep).eq.0) then
            call lpm_usr_particles_io
         endif

         do stage=1,nstage_part

            if (stage .eq. 1) then
               rdt_part = abs(param(12))
               call lpm_set_dt(rdt_part)
               dt_cmt     = rdt_part
               dt         = dt_cmt
               time       = time + dt_cmt
               time_cmt   = time_cmt + dt_cmt
               if (abs(time_integ).eq.1)
     >            call set_tstep_coef_part(rdt_part)

               ! Update coordinates if particle moves outside boundary
               ptdum(2) = dnekclock()
                  call update_particle_location  
               pttime(2) = pttime(2) + dnekclock() - ptdum(2)

c              ! Update where particle is stored at
               ptdum(3) = dnekclock()
                  call move_particles_inproc
               pttime(3) = pttime(3) + dnekclock() - ptdum(3)

               if (two_way.gt.1) then
                  ! Create ghost/wall particles
                  ptdum(4) = dnekclock()
                     call create_extra_particles
                     call sort_local_particles_collisions
                  pttime(4) = pttime(4) + dnekclock() - ptdum(4)
                  
                  ! Send ghost particles
                  ptdum(5) = dnekclock()
                     call send_ghost_particles
                  pttime(5) = pttime(5) + dnekclock() - ptdum(5)
                  
                  ! Projection to Eulerian grid
                  ptdum(6) = dnekclock()
                     call spread_props_grid
                  pttime(6) = pttime(6) + dnekclock() - ptdum(6)
               endif
            endif
            ! Interpolate Eulerian properties to particle location
            ptdum(7) = dnekclock()
               call interp_props_part_location
            pttime(7) = pttime(7) + dnekclock() - ptdum(7)

            ! Evaluate particle force models
            ptdum(8) = dnekclock()
               call usr_particles_forcing  
            pttime(8) = pttime(8) + dnekclock() - ptdum(8)
   
            ! Integrate in time
            ptdum(9) = dnekclock()
               if (abs(time_integ) .eq. 1) call rk3_integrate
               if (abs(time_integ) .eq. 2) call bdf_integrate
            pttime(9) = pttime(9) + dnekclock() - ptdum(9)
   
            ! Update forces
            ptdum(10) = dnekclock()
               call compute_forcing_post_part
            pttime(10) = pttime(10) + dnekclock() - ptdum(10)

         enddo

         ptdum(1) = dnekclock()
      enddo

      if (nid.eq.0) write(6,*) 'FINISHED PRE COLLISIONS - EXITING NOW'
      call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine sort_local_particles_collisions
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'LPM'

      ! here we can actually shrik rdxgp,rdygp,rdzgp if projection 
      ! distance is larger. We only need collsion distance as dpe_max
      ndxgpc = floor( (xdrange(2,1) - xdrange(1,1))/d2chk(3)) +1
      ndygpc = floor( (xdrange(2,2) - xdrange(1,2))/d2chk(3))+1
      ndzgpc = 1
      if (if3d) ndzgpc =floor( (xdrange(2,3) - xdrange(1,3))/d2chk(3))+1

      ! grid spacing for that many spacings
      rdxgpc = (xdrange(2,1) - xdrange(1,1))/real(ndxgpc)
      rdygpc = (xdrange(2,2) - xdrange(1,2))/real(ndygpc)
      rdzgpc = 1.
      if (if3d) rdzgpc = (xdrange(2,3) - xdrange(1,3))/real(ndzgpc)

      ! set real particles ii,jj,kk
      do i=1,n
         rxval = rpart(jx,i)
         ryval = rpart(jy,i)
         rzval = 0.
         if(if3d) rzval = rpart(jz,i)
  
         ii    = floor((rxval-xdrange(1,1))/rdxgpc) 
         jj    = floor((ryval-xdrange(1,2))/rdygpc) 
         kk    = floor((rzval-xdrange(1,3))/rdzgpc) 
         ndum  = ii + ndxgp*jj + ndxgp*ndygp*kk

         ipart(jicx,i) = ii
         ipart(jicy,i) = jj
         ipart(jicz,i) = kk
      enddo
      ! set ghost particles ii,jj,kk
      do i=1,nfptsgp
         rxval = rptsgp(jgpx,i)
         ryval = rptsgp(jgpy,i)
         rzval = 0.
         if(if3d) rzval = rptsgp(jgpz,i)
  
         ii    = floor((rxval-xdrange(1,1))/rdxgpc) 
         jj    = floor((ryval-xdrange(1,2))/rdygpc) 
         kk    = floor((rzval-xdrange(1,3))/rdzgpc) 
         ndum  = ii + ndxgp*jj + ndxgp*ndygp*kk

         iptsgp(jgpicx,i) = ii
         iptsgp(jgpicy,i) = jj
         iptsgp(jgpicz,i) = kk
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine lpm_userf(ix,iy,iz,e,ffxp,ffyp,ffzp,qvolp)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'LPM'

      integer ix,iy,iz,e

      real ffxp,ffyp,ffzp,qvolp
      
      ! particle forcing
      if (two_way .ge.2) then
         if (istep .gt. time_delay) then

            if(ifibm.eq.1) then
            ffxp =  ptw(ix,iy,iz,e,1)/vtrans(ix,iy,iz,e,1)
            ffyp =  ptw(ix,iy,iz,e,2)/vtrans(ix,iy,iz,e,1)
            ffzp =  ptw(ix,iy,iz,e,3)/vtrans(ix,iy,iz,e,1)
            else
            ffxp =  ptw(ix,iy,iz,e,1)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            ffyp =  ptw(ix,iy,iz,e,2)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            ffzp =  ptw(ix,iy,iz,e,3)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            endif
c     if( abs(ptw(ix,iy,iz,e,2)) .ge. 1.d-5)  
c     >       write(*,*)"ptw in userf", ptw(ix,iy,iz,e,2)
c            if (abs(ffyp) .ge. 1.0d-6)
c     >       write(*,*)"ffyp", ffyp

            ! energy coupling for cmt-nek
            if (icmtp .eq. 1) then
               qvolp= ptw(ix,iy,iz,e,5) + rhs_fluidp(ix,iy,iz,e,4)
            else
               qvolp=0.
            endif
         else
            ffxp = 0.0
            ffyp = 0.0
            ffzp = 0.0
         endif
      else
         ffxp = 0.0
         ffyp = 0.0
         ffzp = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bdf_integrate
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)
      
      real s,pmass

      call get_bdf_ext_coefs(beta,alpha,times)

      jx0 = jx
c     move data to previous positions
      do j=0,ndim-1
      do i=1,n
         rpart(ju3+j,i)=rpart(ju2+j,i)
         rpart(ju2+j,i)=rpart(ju1+j,i)
         rpart(ju1+j,i)=rpart(ju0+j,i)
         rpart(jv3+j,i)=rpart(jv2+j,i)
         rpart(jv2+j,i)=rpart(jv1+j,i)
         rpart(jv1+j,i)=rpart(jv0+j,i)
         rpart(jx3+j,i)=rpart(jx2+j,i)
         rpart(jx2+j,i)=rpart(jx1+j,i)
         rpart(jx1+j,i)=rpart(jx0+j,i)
      enddo
      enddo

c     Solve for velocity at time t^n
      do i=1,n
        do j=0,ndim-1
          rhs = rpart(jf0+j,i)
     $        +     beta (1)*rpart(jv1+j,i)
     $        +     beta (2)*rpart(jv2+j,i)
     $        +     beta (3)*rpart(jv3+j,i)
          rpart(jv0+j,i) = rhs / beta(0)
          rhx = beta (1)*rpart(jx1+j,i)
     $        + beta (2)*rpart(jx2+j,i)
     $        + beta (3)*rpart(jx3+j,i) + rpart(jv0+j,i)
          rpart(jx0+j,i) = rhx / beta(0)     ! Implicit solve for x
        enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      FUNCTION MixtPerf_C_GRT_part(G,R,T,icmt)
      REAL G,R,T! INTENT(IN) 
      REAL MixtPerf_C_GRT_part
      integer icmt
      if (icmt .eq. 0) then
         MixtPerf_C_GRT_part = 1.
      else
         MixtPerf_C_GRT_part = SQRT(G*R*T)
      endif

      END
c-----------------------------------------------------------------------
c     Force models below
c-----------------------------------------------------------------------
      subroutine lpm_f_qs
c
c     quasi-steady force
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (abs(part_force(1)).eq.1) then

         rdum = lpmvol_p*lpmdens_p/lpmtau_p

         if (part_force(1) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum  = rdum*(1. - 2.*rphip)/(1. - rphip)**3
         endif

         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

      elseif (abs(part_force(1)).eq.2) then
         rdum = lpmvol_p*lpmdens_p/lpmtau_p

         rmacr= 0.6
         rcd_std = 1.+0.15*lpmre_p**(0.687) + 
     >               0.42*(1.+42500./lpmre_p**(1.16))**(-1)
         rcd_mcr = 1.+0.15*lpmre_p**(0.684) + 
     >               0.513*(1. + 483./lpmre_p**(0.669))**(-1)
         rcd1 = rcd_std + (rcd_mcr - rcd_std)*lpmmach_p/rmacr
         
         rdum = rdum*rcd1

         if (part_force(1) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum  = rdum*(1. - 2.*rphip)/(1. - rphip)**3
         endif

         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

      elseif (abs(part_force(1)).eq.3) then
         rdum = lpmvol_p*lpmdens_p/lpmtau_p

         rcd_std = 1.+0.15*lpmre_p**(0.687) + 
     >               0.42*(1.+42500./lpmre_p**(1.16))**(-1)
         rdum = rdum*rcd_std

         if (part_force(1) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum  = rdum*(1. - 2.*rphip)/(1. - rphip)**3
         endif

         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

      elseif (abs(part_force(1)).eq.4) then
         ! ergun
         rbeta1 = 150.*lpmvolfrac_p*lpmvolfrac_p*lpmvisc_f/lpmvolfrac_f/
     >            lpmdiam_p**2 + 1.75*lpmvolfrac_p*lpmdens_f*
     >            lpmvdiff_pf/lpmdiam_p

         ! wen-yu
         rrep = lpmvolfrac_f*lpmre_p
         if (rrep .lt. 1000) then
            rcd = 24./rrep*(1. + 0.15*rrep**(0.687))
         else
            rcd = 0.44
         endif

         rbeta2 = 0.75*rcd*lpmvolfrac_f**(-2.65)*
     >        lpmvolfrac_p*lpmvolfrac_f*lpmdens_f*lpmvdiff_pf/lpmdiam_p

         ! stiching
         rs = 0.2
         rpp = atan(150 * 1.75*(lpmvolfrac_p - rs))/pi + 0.5
         rbeta = rpp*rbeta1 + (1.-rpp)*rbeta2


         rpart(jfqs+j,i) = rpart(jvol,i)*rbeta/rphip    
     >              *(rpart(ju0+j,i) - rpart(jv0+j,i))

         rdum = lpmvol_p*rbeta/lpmvolfrac_p
         lpmforce(1) = rdum*(lpmv_f(1) - lpmv_p(1))
         lpmforce(2) = rdum*(lpmv_f(2) - lpmv_p(2))
         lpmforce(3) = rdum*(lpmv_f(3) - lpmv_p(3))

c         write(*,*) "qs drag rd=",rdum,"vel=",lpmv_f(2), lpmv_p(2),
c     > "lpmforce", lpmforce(1), lpmforce(2), lpmforce(3) 
                 
         if (abs(lpmvolfrac_p) .lt. 1E-12) then
             write(6,*) 'Use different drag model w/o volume frac /0'
             call exitt
         endif

      elseif (abs(part_force(1)).eq.5) then !!! ibm force 

         pmassf= lpmvol_p * lpmdens_f

         lpmforce(1) = (lpmv_f(1) - lpmv_p(1)) / dt  * pmassf
         lpmforce(2) = (lpmv_f(2) - lpmv_p(2)) / dt  * pmassf
         lpmforce(3) = (lpmv_f(3) - lpmv_p(3)) / dt  * pmassf

         
c         write(*,*) "IBvel",lpmv_f(2), lpmv_p(2), "dt =", dt
c     >   , "pmassf" , pmassf, lpmvol_p
c     >   , "lpmforce", lpmforce(1), lpmforce(2), lpmforce(3) 

      elseif (part_force(1).eq.0) then
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_f_un
c
c     undisturbed force
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (abs(part_force(2)) .ne. 0) then
         lpmforce(1) = -lpmvol_p*lpmDuDt(1)
         lpmforce(2) = -lpmvol_p*lpmDuDt(2)
         lpmforce(3) = -lpmvol_p*lpmDuDt(3)
      elseif (part_force(2) .eq. 0) then
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_f_iu
c
c     inviscid unsteady force
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (abs(part_force(3)) .eq. 1) then
         ! shape
         rdum = 0.5

         ! volume fraction
         if (part_force(3) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum = rdum*(1.+2.*rphip)/(1.-rphip)
         endif

         ! set coeff for implicit dv/dt solve
         rpart(jcmiu,lpmi) = rdum
         
         ! but don't add a fluid contribution for now, but should later
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0

      elseif (abs(part_force(3)) .eq. 2) then
         ! shape
         rdum = 0.5
         ! mach
         rdum = rdum*(1. + 1.8*lpmmach_p**2 + 7.6*lpmmach_p**4)

         ! volume fraction
         if (part_force(3) .lt. 0) then
            rphip = min(0.3,lpmvolfrac_p)
            rdum = rdum*(1.+2.*rphip)/(1.-rphip)
         endif

         ! set coeff for implicit dv/dt solve
         rpart(jcmiu,lpmi) = rdum
         
         ! but don't add a fluid contribution for now, but should later
         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0

      elseif (part_force(3) .eq. 0) then

         lpmforce(1) = 0.0
         lpmforce(2) = 0.0
         lpmforce(3) = 0.0

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_q_qs
c
c     quasi-steady heat transfer
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real nu_g,rdum,rpra,kappa_g

      nu_g    = lpmvisc_f/lpmdens_f ! kinematic viscosity
      rpra    = nu_g/lpmkappa_f
      rdum    = 2.*pi*lpmkappa_f*lpmdiam_p

      if (abs(part_force(4)) .gt. 0) then
         rdum    = rdum*(1. + 0.3*sqrt(lpmre_p)*rpra**(2./3.)) !re < 500
         lpmheat = rdum*(lpmtemp_f - lpmtemp_p)
      elseif (part_force(4) .eq. 0) then
         lpmheat = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_q_uu
c
c     undisturbed heat transfer
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      if (part_force(5) .eq. 0) then
         lpmheat = 0.0 ! nothing for now
      endif

      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine sphere_marker_distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Purpose: Generate 3D equispaced points of an sphere of the sphere 
!     Description:
!          following the IBM 3D code by Hyungoo Lee
!          ref. E.B Saff and A.B.J. Kuijlaars (1996)
!
!     Input:
!       nwe   : spheres in local MPI rank
!       ibm_center ( nsi+ii, 3 )
!       ibm_diam ( nsi+nn )
!       n_dh  : number of grid points in the diameter
!       n_IBMpart: - number of IBM markers in local MPI rank
!
!     Intermediate variables
!       n_l   : approx. # of Lagrangian points
!       h     : wokrer marker size
!       r     : radius
!       hk (i):
!       tht(i): theta
!     
!     Iterators
!       nn : loop through number of particles: n_IBMpart
!       k  : loop through all markers 
!       l  : numbering of a point on a sphere
!       i  : loop through worker markers on a sphere, i = l, not used for now, for mulit-layers
!       j  : = 0,1,2 - x,y,z component
!       ll : number of layers per particle, not used for now
!     
!     Output:
!       rpart(jx + j,  nn) - Lagrange points locations (j=0,2) 
!       rpart(jx1 + j, nn) - Lagrange points locations (j=0,2) 
!       rpart(jx2 + j, nn) - Lagrange points locations (j=0,2) 
!       rpart(jx3 + j, nn) - Lagrange points locations (j=0,2) 
!       rpart(jdp,     nn) - Diameter 
!       rpart(jvol,    nn) - Volume 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real x_c(num_p),y_c(num_p),z_c(num_p),r_c(num_p)
      real vel_c(num_p,3), omg_c(num_p,3), the_c(num_p,3)
      real hk(n_ll), tht(n_ll), ph(n_ll)

c local
      real h, r
      integer ll,nn,i,j,k,l,mm,i1,i2,j1,j2,k1,k2
      real x_o, y_o, z_o, dv_l
      integer seq_distribute, nsi, random_distribute
      real n_dh_l
      integer        ndef
      common         ndef
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Specifying the particles number, positions and dimensions
c      n_dh=int(2.d0*r_c/deltax)+20
!      h = deltax
      ! num_p = num_of_IBMpart    ! number of particles
      ! num_of_IBMpart  = num_p 
      ! n_l(n) = n_markers number of markers per sphere
      ! n_dh(1)  = 25      !? current understanding, number of points along the diameter,
      ! in other words, diameter_of_sphere/lagrange_point_spacing
      ! n_l_max = 0.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rpi = 4.*atan(1.0)

      write(6,'(A,I6,A,I6)')"Place ",n_IBMpart," spheres in rank",nid
      open(unit=95, file='lagrangian_pts')      

      seq_distribute = 1
      random_distribute = 0
      if(seq_distribute.eq.1) then
            ! ndef: if(nid <  ndef), n_IBMpart+1
            !       if(nid >= ndef), n_IBMpart+0
         if(nid.lt.ndef) nsi = nid * n_IBMpart
         if(nid.ge.ndef) nsi = ndef * (n_IBMpart + 1) + 
     $                        (nid - ndef) * n_IBMpart
         print*,"nid,nsi,ndef",nid,nsi,ndef
      endif
      
      
      ll = 1 ! number of layers per particle, not used for now
!     Step 1 estaimate how many point placed on sphere
      do nn=1, n_IBMpart            ! number of particles in current process
         n_l(nn) = 0  ! number of lagrange points on each surface of each particles

c        get local copy of the global array
         x_c(nn) = ibm_center(nsi + nn, 1)
         y_c(nn) = ibm_center(nsi + nn, 2)
         z_c(nn) = ibm_center(nsi + nn, 3)
         r_c(nn) = ibm_diam  (nsi + nn ) / 2.0d0
         n_dh_l = n_dh(nsi+nn) 

         do j=0,2
            vel_c(nn,j+1) = ibm_vel_init(nsi+nn,j+1) ! translate velocity
            omg_c(nn,j+1) = ibm_omg_init(nsi+nn,j+1) ! rotating velocity
            the_c(nn,j+1) = ibm_the_init(nsi+nn,j+1) ! angle
         enddo
         
         if(random_distribute.eq.1) then
            x_c(nn) = unif_random(rxbo(1,1),rxbo(2,1))
            y_c(nn) = unif_random(rxbo(1,2),rxbo(2,2))
            z_c(nn) = unif_random(rxbo(1,3),rxbo(2,3))
            r_c(nn) = dp(1) / 2.0d0                   ! for uniform only
            n_dh_l  = n_dh(1)                         ! for uniform distribution 
            do j=0,2
               vel_c(nn,j+1) = 0.0 !ibm_vel_init(nsi+nn,j+1) ! translate velocity
               omg_c(nn,j+1) = 0.0 !ibm_omg_init(nsi+nn,j+1) ! rotating velocity
               the_c(nn,j+1) = 0.0
            enddo
         endif

c        intermediate variables         
         r   = r_c(nn)          ! radius
         h   = 2*r/n_dh_l       ! thickness 
         n_l(nn) = int(rpi/3.0*(12.0*(r/h)**2.0+1.0))

         if (n_l(nn) .gt. n_ll) then
            write(26,*) 'error in Lagrangian array size'
            print *, 'error in Lagrangian array size'
            stop
         endif

c     assign data array
         rpart(jx + 0, nn) =  x_c(nn)
         rpart(jx + 1, nn) =  y_c(nn)
         rpart(jx + 2, nn) =  z_c(nn)
         rpart(jdp,  nn)   =  2 * r_c(nn)
         rpart(jvol, nn)   =  4.0/3.0 * rpi * r_c(nn)**3
c         write(6,*)"dp",rpart(jdp,  nn),"vol",rpart(jvol, nn)
         do j =0,2
            rpart(jx1 + j , nn) = rpart(jx + j, nn)
            rpart(jx2 + j , nn) = rpart(jx + j, nn)
            rpart(jx3 + j , nn) = rpart(jx + j, nn) 
         enddo

         do j=0,2
            rpart(jv0 + j,   nn) = vel_c(nn,j+1)    ! intialize velocity 
            rpart(jangvel0 + j,nn) = omg_c(nn,j+1)  ! intialize angular velocity
            rpart(jangle0 + j,nn)  = the_c(nn,j+1)  ! intialize angle 
         enddo

         do j =0,2
            rpart(jv1 + j , nn) = rpart(jv0 + j, nn) ! prev time step
            rpart(jv2 + j , nn) = rpart(jv0 + j, nn)
            rpart(jv3 + j , nn) = rpart(jv0 + j, nn) 
         enddo

         do j =0,2
            rpart(jangvel1 + j , nn) = rpart(jangvel0 + j, nn) ! prev time step
            rpart(jangvel2 + j , nn) = rpart(jangvel0 + j, nn)
            rpart(jangvel3 + j , nn) = rpart(jangvel0 + j, nn)                
         enddo

         do j =0,2
            rpart(jangle1 + j , nn) = rpart(jangle0 + j, nn) ! prev time step
            rpart(jangle2 + j , nn) = rpart(jangle0 + j, nn)
            rpart(jangle3 + j , nn) = rpart(jangle0 + j, nn)                
         enddo

         if(nid.le.1) 
     $    write(6,2022) nid, nn, n_l(nn), 
     $        (rpart(jx+j,nn),j=0,2),rpart(jdp,nn),rpart(jvol,nn)
     $        ,(rpart(jv0+j,nn),j=0,2),(rpart(jangvel0+j,nn),j=0,2)
 2022    format(I6,' Queen #',2I6,' xyz',3F8.3,',dp=',F8.3
     $    ,', vol=',E12.4,', Vel=',3E12.4,', Ang Vel=',3E12.4)

      enddo

!     Step 2 calculate the location for each lagrange point on a sphere
      k = n_IBMpart  ! num_of_IBMpart: total
      do nn = 1, n_IBMpart

         h   = 2 * r_c(nn) / n_dh(nn)  ! thickness             
         r   = r_c(nn)                 ! radius
         i = 0

         do l = 1, n_l(nn)  ! l - is the numbering of a point on a sphere

            i = i + 1
            mm = int(rpi/3.d0*(12.d0*(r/h)**2.d0+1.d0)) ! number of NL in Uhlmann paper
            hk(i) = -1.d0+2.d0*dble(i-1)/dble(mm-1)
            tht(i) = acos(hk(i)) ! theta
c
            if (i .eq. 1 .or. i .eq. mm) then
               ph(i) = 0.0d0
            else
               ph(i) = ph(i-1)+3.809d0/sqrt(dble(mm))/
     &                sqrt(1.d0-hk(i)**2.0d0)
            endif
            if (ph(i) .gt. 2.d0*rpi) then
               ph(i) = ph(i)-2.d0*rpi
            endif

c     position of Lagragian points
            x_o = x_c(nn) + (r-0.3333*h) * sin(tht(i)) * cos(ph(i))!retract marker position based on the Breugem2012
            y_o = y_c(nn) + (r-0.3333*h) * sin(tht(i)) * sin(ph(i))
            z_o = z_c(nn) + (r-0.3333*h) * cos(tht(i))
!            x_o = x_c(nn) + r * sin(tht(i)) * cos(ph(i))
!            y_o = y_c(nn) + r * sin(tht(i)) * sin(ph(i))
!            z_o = z_c(nn) + r * cos(tht(i))

c     volume of each Lagrangian cell (M.Uhlmann 2005)
c            ds(l,n)=4.d0*rpi*r**2.d0/dble(mm)
            dv_l = rpi*h/3.d0/dble(mm)*(12.d0*r**2.d0+h**2.d0)
c
            if (i+1.gt.mm) then
               i=0
               r=r-h
            endif

            ! k - is the total numbering of a point globally
            k = k + 1
            rpart(jx + 0, k) =  x_o
            rpart(jx + 1, k) =  y_o
            rpart(jx + 2, k) =  z_o
            rpart(jvol, k)   =  dv_l
            rpart(jdp,  k)   =  h
                       
            do j = 0,2
               rpart(jx1 + j , k) = rpart(jx + j, k)
               rpart(jx2 + j , k) = rpart(jx + j, k)
               rpart(jx3 + j , k) = rpart(jx + j, k) 
            enddo

            do j=0,2
               rpart(jv0 + j,   k) = vel_c(nn, j+1) ! intialize velocity
               rpart(jangvel0 + j,k) = omg_c(nn, j+1)
            enddo

            do j =0,2
               rpart(jv1 + j , k) = rpart(jv0 + j, k)
               rpart(jv2 + j , k) = rpart(jv0 + j, k)
               rpart(jv3 + j , k) = rpart(jv0 + j, k)                
            enddo
            
         enddo  ! l=n_l(n)
      enddo  !  nn = 1, n_IBMpart
    
      close (unit=95)

      if (k .ne. nwe) write(*,*) "Warning: Particle number must match"
     $ , nwe, k 
      nwe = k

      do nn = 1, n_IBMpart 
         if (n_l_max .lt. n_l(nn)) n_l_max=n_l(nn)
      enddo

      return
      END SUBROUTINE sphere_marker_distribution
c-------------------------------------------------------------------

c-------------------------------------------------------------------
      subroutine save_ibm_neighbor_bin

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
      
      nlist_ibm = nlist
      ndxgp_ibm = ndxgp
      ndygp_ibm = ndygp
      ndzgp_ibm = ndzgp
      
      nlist_total  = nlist_ibm * ngpvc
      call icopy(ngp_valsp_ibm, ngp_valsp, nlist_total)
      nlist_total1 = lx1*ly1*lz1*lelt*4
      call icopy(mod_gp_grid_ibm, mod_gp_grid, nlist_total1)

      rdxgp_ibm = rdxgp
      rdygp_ibm = rdygp
      rdzgp_ibm = rdzgp
      
      end subroutine save_ibm_neighbor_bin
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      subroutine create_ghost_particles_full_ibm
c
c     this routine will create ghost particles by checking if particle
c     is within d2chk of element faces
c
c     ghost particle x,y,z list will be in
c     rptsgp(jgpx,j), rptsgp(jgpy,j),rptsgp(jgpz,j),
c     while processor and local element id are in
c     iptsgp(jgppt,j) and iptsgp(jgpes,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real xdlen,ydlen,zdlen,rxdrng(3),rxnew(3)
      integer iadd(3),ntypesl(7)

      integer ngp_trim(llpart)
! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      xdlen = xdrange(2,1) - xdrange(1,1)
      ydlen = xdrange(2,2) - xdrange(1,2)
      zdlen = -1.
      if (if3d) zdlen = xdrange(2,3) - xdrange(1,3)
      if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) xdlen = -1 ! non-periodic
      if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) ydlen = -1
      if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c      nlist = nlist_ibm
c      
c      rdxgp = rdxgp_ibm
c      rdygp = rdygp_ibm
c      rdzgp = rdzgp_ibm

c      ndxgp = ndxgp_ibm
c      ndygp = ndygp_ibm
c      ndzgp = ndzgp_ibm
      
c      nlist_total  = nlist_ibm * ngpvc
c      call icopy(ngp_valsp, ngp_valsp_ibm, nlist_total)

c      nlist_total1 = lx1*ly1*lz1*lelt*4
c      call icopy(mod_gp_grid, mod_gp_grid_ibm, nlist_total1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      nfptsgp_ibm = 0
      
      !print*,nid, "x-range",vlmin(xm1,lx1*ly1*lz1*lelt),
      !$ vlmax(xm1,lx1*ly1*lz1*lelt)
      !print*,nid, "y-range",vlmin(ym1,lx1*ly1*lz1*lelt),
      !$ vlmax(ym1,lx1*ly1*lz1*lelt)
      !print*,nid, "z-range",vlmin(zm1,lx1*ly1*lz1*lelt),
      !$ vlmax(zm1,lx1*ly1*lz1*lelt)

      do ip=1,n

         if (ipart(jrole,ip) .ne.1) cycle ! only for ibm role is queen

         rxval = rpart(jx,ip)
         ryval = rpart(jy,ip)
         rzval = 0.
         if(if3d) rzval = rpart(jz,ip)

         iip    = floor((rxval-xdrange(1,1))/rdxgp_ibm) 
         jjp    = floor((ryval-xdrange(1,2))/rdygp_ibm) 
         kkp    = floor((rzval-xdrange(1,3))/rdzgp_ibm) 
         ndump  = iip + ndxgp_ibm*jjp + ndxgp_ibm*ndygp_ibm*kkp
c         write(*,*) "bin physical #",nid,rxval,ryval,rzval,iip,jjp,kkp
c         write(*,*) "bin x range",rdxgp_ibm*iip + xdrange(1,1),
c     $        rdxgp_ibm*(iip+1) + xdrange(1,1)
c         write(*,*) "bin y range",rdxgp_ibm*jjp + xdrange(1,2),
c     $        rdxgp_ibm*(jjp+1) + xdrange(1,2)
c         write(*,*) "bin z range",rdxgp_ibm*kkp + xdrange(1,3),
c     $        rdxgp_ibm*(kkp+1) + xdrange(1,3)

         do i=1,nlist_ibm
         ii = ngp_valsp_ibm(3,i)
         jj = ngp_valsp_ibm(4,i)
         kk = ngp_valsp_ibm(5,i)

         ndum = ngp_valsp_ibm(2,i)
         nrank= ngp_valsp_ibm(1,i)
         !write(*,*) "search within bin: ",nid,nrank,ndum,ii,jj,kk
         ! add this box
         if (nid  .ne. nrank) then
         if (ndum .eq. ndump) then

            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
         
            iadd(1)  = 0
            iadd(2)  = ngp_valsp_ibm(1,i) ! where this bin belongs to
            iadd(3)  = ipart(je0,ip)      ! 
c            write(*,*) "ghost IBM send from",nid,"to",nrank,ii,jj,kk
            call add_a_ghost_particle_ibm(rxnew,iadd,ip)

         endif
         endif

      enddo
      enddo

      do ip=1,n

         if (ipart(jrole,ip) .ne.1) cycle ! ibm role is queen

         rxval = rpart(jx,ip)
         ryval = rpart(jy,ip)
         rzval = 0.
         if(if3d) rzval = rpart(jz,ip)

         iip    = floor((rxval-xdrange(1,1))/rdxgp_ibm) 
         jjp    = floor((ryval-xdrange(1,2))/rdygp_ibm) 
         kkp    = floor((rzval-xdrange(1,3))/rdzgp_ibm) 
         ndump  = iip + ndxgp_ibm*jjp + ndxgp_ibm*ndygp_ibm*kkp

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            ! print*,"faces",iig,jjg,kkg,ndygp_ibm-1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp_ibm-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp_ibm)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp_ibm-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp_ibm)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp_ibm-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp_ibm)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif

            !print*,"faces before cycle",iig,jjg,kkg
            iflgsum = iflgx + iflgy + iflgz
            if (iflgsum .ne. 1) cycle

            ndumn = iig + ndxgp_ibm*jjg + ndxgp_ibm*ndygp_ibm*kkg
            
            do i=1,nlist_ibm
               ndum  = ngp_valsp_ibm(2,i)
               nrank = ngp_valsp_ibm(1,i)
               iin   = ngp_valsp_ibm(3,i)
               jjn   = ngp_valsp_ibm(4,i)
               kkn   = ngp_valsp_ibm(5,i)

               if (ndumn .eq. ndum) then
                  ibctype = iflgx+iflgy+iflgz
                 
                  if (nrank .ne. nid .or. 
     >                    (nrank .eq. nid) .and. (iflgsum.eq.1)) then
                 
                      rxnew(1) = rxval
                      rxnew(2) = ryval
                      rxnew(3) = rzval
                 
                      iadd(1) = ii1
                      iadd(2) = jj1
                      iadd(3) = kk1
                 
                      call check_periodic_gp_ibm(rxnew,rxdrng,iadd)
                 
                      iadd(1)  = 0
                      iadd(2)  = nrank
                      iadd(3)  = ipart(je0,ip) 
                      !print*,"face ghost",iin,jjn,kkn
                      call add_a_ghost_particle_ibm(rxnew,iadd,ip)
                        
                  endif
               endif
            enddo
         enddo
         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp_ibm-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp_ibm)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp_ibm-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp_ibm)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp_ibm-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp_ibm)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif

            !print*,"edges",iig,jjg,kkg
            iflgsum = iflgx + iflgy + iflgz
            if (iflgsum .ne. 2) cycle

            ndumn = iig + ndxgp_ibm*jjg + ndxgp_ibm*ndygp_ibm*kkg


            do i=1,nlist_ibm
               ndum  = ngp_valsp_ibm(2,i)
               nrank = ngp_valsp_ibm(1,i)
               iin   = ngp_valsp_ibm(3,i)
               jjn   = ngp_valsp_ibm(4,i)
               kkn   = ngp_valsp_ibm(5,i)

               if (ndumn .eq. ndum) then
                  ibctype = iflgx+iflgy+iflgz
                 
                  if (nrank .ne. nid .or. 
     >                    (nrank .eq. nid) .and. (iflgsum.eq.2)) then
                 
                      rxnew(1) = rxval
                      rxnew(2) = ryval
                      rxnew(3) = rzval
                 
                      iadd(1) = ii1
                      iadd(2) = jj1
                      iadd(3) = kk1
                 
                      call check_periodic_gp_ibm(rxnew,rxdrng,iadd)
                 
                      iadd(1)  = 0
                      iadd(2)  = nrank
                      iadd(3)  = ipart(je0,ip)
                      !print*,"edge ghost",iin,jjn,kkn
                      call add_a_ghost_particle_ibm(rxnew,iadd,ip)
                        
                  endif
               endif
            enddo
         enddo
         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ndxgp_ibm-1) then
               iflgx = 1
               iig =modulo(iig,ndxgp_ibm)
               if (abs(bc_part(1)) + abs(bc_part(2)) .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ndygp_ibm-1) then
               iflgy = 1
               jjg =modulo(jjg,ndygp_ibm)
               if (abs(bc_part(3)) + abs(bc_part(4)) .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ndzgp_ibm-1) then
               iflgz = 1  
               kkg =modulo(kkg,ndzgp_ibm)
               if (abs(bc_part(5)) + abs(bc_part(6)) .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            if (iflgsum .ne. 3) cycle

            ndumn = iig + ndxgp_ibm*jjg + ndxgp_ibm*ndygp_ibm*kkg

            do i=1,nlist_ibm
               ndum  = ngp_valsp_ibm(2,i)
               nrank = ngp_valsp_ibm(1,i)
               iin   = ngp_valsp_ibm(3,i)
               jjn   = ngp_valsp_ibm(4,i)
               kkn   = ngp_valsp_ibm(5,i)

               if (ndumn .eq. ndum) then
                  ibctype = iflgx+iflgy+iflgz
                 
                  if (nrank .ne. nid .or. 
     >                    (nrank .eq. nid) .and. (iflgsum.eq.3)) then
                 
                      rxnew(1) = rxval
                      rxnew(2) = ryval
                      rxnew(3) = rzval
                 
                      iadd(1) = ii1
                      iadd(2) = jj1
                      iadd(3) = kk1
                 
                      call check_periodic_gp_ibm(rxnew,rxdrng,iadd)
                 
                      iadd(1)  = 0
                      iadd(2)  = nrank
                      iadd(3)  = ipart(je0,ip)

                      call add_a_ghost_particle_ibm(rxnew,iadd,ip)
                        
                  endif
               endif
            enddo
         enddo
      enddo

! trim down so no duplicates  in total "nfptsgp_ibm" ghost queen markers
! in list iptsgp_ibm
!     1 jgpps
!     2 jgp_pid1, jgp_pid2, jgp_pid3
      ! trim down so no duplicates

      do i=1,nfptsgp_ibm
         ngp_trim(i) = 1
      enddo

      do i=1,  nfptsgp_ibm
         if(iptsgp_ibm(jgpps,i) .eq. nid)
     &   ngp_trim(i) = 0
      do j=i+1,nfptsgp_ibm
         if (ngp_trim(j) .eq. 1) then
         
         if ( iptsgp_ibm(jgpps,i)    .eq. iptsgp_ibm(jgpps,j)    .and.
     &        iptsgp_ibm(jgp_pid1,i) .eq. iptsgp_ibm(jgp_pid1,j) .and.
     &        iptsgp_ibm(jgp_pid2,i) .eq. iptsgp_ibm(jgp_pid2,j) .and.            
     &        iptsgp_ibm(jgp_pid3,i) .eq. iptsgp_ibm(jgp_pid3,j) )
     &   ngp_trim(j) = 0

         endif
      enddo
      enddo
!      write(*,*) 'ngp_trim', (ngp_trim(j),j=1,nfptsgp_ibm)
      ic = 0
      do i=1,nfptsgp_ibm
         if (ngp_trim(i) .eq. 1) then
            ic = ic + 1
            call icopy(iptsgp_ibm(1,ic),iptsgp_ibm(1,i),ligp) 
            call  copy(rptsgp_ibm(1,ic),rptsgp_ibm(1,i),lrgp) 
          endif
      enddo

      nfptsgp_ibm = ic
      
!      do i = 1, nfptsgp_ibm
!      print*,"GP in",nid,i,rptsgp_ibm(jgpy,i),iptsgp_ibm(jgpps,i)
!      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine check_periodic_gp_ibm(rxnew,rxdrng,iadd)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
c
      real rxnew(3), rxdrng(3)
      integer iadd(3), irett(3), ntype, ntypel(7)

      xloc = rxnew(1)
      yloc = rxnew(2)
      zloc = rxnew(3)

      xdlen = rxdrng(1)
      ydlen = rxdrng(2)
      zdlen = rxdrng(3)

      ii = iadd(1)
      jj = iadd(2)
      kk = iadd(3)

      irett(1) = 0
      irett(2) = 0
      irett(3) = 0

      if (xdlen .gt. 0 ) then
      if (ii .ge. ndxgp_ibm) then
         xloc = xloc - xdlen
         irett(1) = 1
         goto 123
      endif
      endif
      if (xdlen .gt. 0 ) then
      if (ii .lt. 0) then
         xloc = xloc + xdlen
         irett(1) = 1
         goto 123
      endif
      endif

  123 continue    
      if (ydlen .gt. 0 ) then
      if (jj .ge. ndygp_ibm) then
         yloc = yloc - ydlen
         irett(2) = 1
         goto 124
      endif
      endif
      if (ydlen .gt. 0 ) then
      if (jj .lt. 0) then
         yloc = yloc + ydlen
         irett(2) = 1
         goto 124
      endif
      endif
  124 continue

      if (if3d) then
         if (zdlen .gt. 0 ) then
         if (kk .ge. ndzgp_ibm) then
            zloc = zloc - zdlen
            irett(3) = 1
            goto 125
         endif
         endif
         if (zdlen .gt. 0 ) then
         if (kk .lt. 0) then
            zloc = zloc + zdlen
            irett(3) = 1
            goto 125
         endif
         endif
      endif
  125 continue

      rxnew(1) = xloc
      rxnew(2) = yloc
      rxnew(3) = zloc

      return
      end

c----------------------------------------------------------------------
      subroutine add_a_ghost_particle_ibm(rxnew,iadd,i)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer i
      real    rxnew(3)
      integer iadd(3)

      nfptsgp_ibm = nfptsgp_ibm + 1

      rptsgp_ibm(jgpx,    nfptsgp_ibm)   = rxnew(1)       ! x loc
      rptsgp_ibm(jgpy,    nfptsgp_ibm)   = rxnew(2)       ! y loc
      rptsgp_ibm(jgpz,    nfptsgp_ibm)   = rxnew(3)       ! z loc
      rptsgp_ibm(jgpfh,   nfptsgp_ibm)   = rpart(jf0,i)   ! hyd. force x
      rptsgp_ibm(jgpfh+1, nfptsgp_ibm)   = rpart(jf0+1,i) ! hyd. force y
      rptsgp_ibm(jgpfh+2, nfptsgp_ibm)   = rpart(jf0+2,i) ! hyd. force z
      rptsgp_ibm(jgpvol,  nfptsgp_ibm)   = rpart(jvol,i)  ! particle volum
      rptsgp_ibm(jgprpe,  nfptsgp_ibm)   = rpart(jrpe,i)  ! particle rp eff
      rptsgp_ibm(jgpspl,  nfptsgp_ibm)   = rpart(jspl,i)  ! spl
      rptsgp_ibm(jgpg0,   nfptsgp_ibm)   = rpart(jg0,i)   ! work done by forc
      rptsgp_ibm(jgpq0,   nfptsgp_ibm)   = rpart(jq0,i)   ! heating from part 
      rptsgp_ibm(jgpv0,   nfptsgp_ibm)   = rpart(jv0,i)   ! particle velocity
      rptsgp_ibm(jgpv0+1, nfptsgp_ibm)   = rpart(jv0+1,i) ! particle velocity
      rptsgp_ibm(jgpv0+2, nfptsgp_ibm)   = rpart(jv0+2,i) ! particle velocity

      rptsgp_ibm(jgp_jx1,   nfptsgp_ibm)   = rpart(jx1,i)   ! previous x loc
      rptsgp_ibm(jgp_jx1+1, nfptsgp_ibm)   = rpart(jx1+1,i) ! previous y loc
      rptsgp_ibm(jgp_jx1+2, nfptsgp_ibm)   = rpart(jx1+2,i) ! previous z loc

      rptsgp_ibm(jgp_angvel0,   nfptsgp_ibm) = rpart(jangvel0,  i) ! rotating velocity
      rptsgp_ibm(jgp_angvel0+1, nfptsgp_ibm) = rpart(jangvel0+1,i) ! rotating velocity
      rptsgp_ibm(jgp_angvel0+2, nfptsgp_ibm) = rpart(jangvel0+2,i) ! rotating velocity

      rptsgp_ibm(jgp_angle0,   nfptsgp_ibm) = rpart(jangle0,  i) ! rotating velocity
      rptsgp_ibm(jgp_angle0+1, nfptsgp_ibm) = rpart(jangle0+1,i) ! rotating velocity
      rptsgp_ibm(jgp_angle0+2, nfptsgp_ibm) = rpart(jangle0+2,i) ! rotating velocity
 
!      rptsgp_ibm(jgp_f0,  nfptsgp_ibm)   = rpart(jf0,i)   ! forcing for ibm integration
!      rptsgp_ibm(jgp_f0+1,nfptsgp_ibm)   = rpart(jf0+1,i) !
!      rptsgp_ibm(jgp_f0+2,nfptsgp_ibm)   = rpart(jf0+2,i) !

      iptsgp_ibm(jgpiic,  nfptsgp_ibm)   = iadd(1)        ! use in collisions
      iptsgp_ibm(jgpps,   nfptsgp_ibm)   = iadd(2)        ! overwritten mpi
      iptsgp_ibm(jgppt,   nfptsgp_ibm)   = iadd(2)        ! dest. mpi rank
      iptsgp_ibm(jgpes,   nfptsgp_ibm)   = iadd(3)        ! dest. elment
      iptsgp_ibm(jgp_back,nfptsgp_ibm)   = nid            ! record where the physical particle is
      
      iptsgp_ibm(jgp_pid1,    nfptsgp_ibm)    = ipart(jpid1,i)  ! id 1
      iptsgp_ibm(jgp_pid2,    nfptsgp_ibm)    = ipart(jpid2,i)  ! id 2
      iptsgp_ibm(jgp_pid3,    nfptsgp_ibm)    = ipart(jpid3,i)  ! id 3
      iptsgp_ibm(jgp_role,    nfptsgp_ibm)    = ipart(jrole,i)  ! id 4
      iptsgp_ibm(jgp_queen,   nfptsgp_ibm)    = ipart(jqueen,i) ! id 5

      iptsgp_ibm(jgp_worker1, nfptsgp_ibm)    = ipart(jworker1,i) ! first worker id 
      iptsgp_ibm(jgp_nlm,     nfptsgp_ibm)    = ipart(jnlm,i)     ! total workers
      iptsgp_ibm(jgp_ax,      nfptsgp_ibm)    = 0                 ! total workers

!      print*,"Created GP",nfptsgp_ibm, "from,",nid,ipart(jpid1,i),"send"
!     $     ,iptsgp_ibm(jgpps,nfptsgp_ibm), "jgp_back = ",
!     $     jgp_back, jgp_pid1, jgp_pid2, jgp_pid3
!     $     ,iptsgp_ibm(jgp_back,nfptsgp_ibm)
!     $     ,iptsgp_ibm(jgp_pid1,nfptsgp_ibm)
!     $     ,iptsgp_ibm(jgp_pid2,nfptsgp_ibm)
!     $     ,iptsgp_ibm(jgp_pid3,nfptsgp_ibm) 
      return
      end

c-----------------------------------------------------------------------
      subroutine send_ghost_particles_ibm
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
!     -----------
!     nfptsgp_ibm: number of columns
!     iptsgp_ibm
!     nigp_ibm = 8
!     1 - [\         ] jpiic : if gp used in collisions or not
!     2 - [ \        ] jgpps : Pointer to proc id for data swap
!     3 - [  \       ] jgppt : findpts return processor id
!     4 - [   \      ] jgpes : Destination element to be sent to
!     5 - [    \     ] jgpicx : 
!     6 - [     \    ] jgpicy :
!     7 - [      \   ] jgpicz :
!     8 - [       \  ] jgp_back :
!
!     rptsgp_ibm:
!     nrgp_ibm = 14

      common /myparth/ i_fp_hndl, i_cr_hndl
      logical partl         
      ! send ghost particles
      call fgslib_crystal_tuple_transfer(i_cr_hndl,nfptsgp_ibm,llpart_gp
     $        ,iptsgp_ibm,nigp_ibm,partl,0,rptsgp_ibm,nrgp_ibm,jgpps)   ! jgpps is overwri

      ! sort by element for ghost particles for projection performance
      call fgslib_crystal_tuple_sort    (i_cr_hndl,nfptsgp_ibm
     $        ,iptsgp_ibm,nigp_ibm,partl,0,rptsgp_ibm,nrgp_ibm,jgpes,1)

      return
      end

c-----------------------------------------------------------------------
      subroutine send_back_ghost_particles_ibm
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
!     -----------
!     nfptsgp_ibm: number of columns
!     iptsgp_ibm
!     nigp_ibm = 8
!     1 - [\         ] jpiic : if gp used in collisions or not
!     2 - [ \        ] jgpps : Pointer to proc id for data swap
!     3 - [  \       ] jgppt : findpts return processor id
!     4 - [   \      ] jgpes : Destination element to be sent to
!     5 - [    \     ] jgpicx : 
!     6 - [     \    ] jgpicy :
!     7 - [      \   ] jgpicz :
!     8 - [       \  ] jgp_back :
!
!     rptsgp_ibm:
!     nrgp_ibm = 14

      common /myparth/ i_fp_hndl, i_cr_hndl
      logical partl         
      ! send ghost particles
      !print*,"send back,IBM",nid,nfptsgp_ibm,nigp_ibm,nrgp_ibm,jgpps

      call fgslib_crystal_tuple_transfer(i_cr_hndl,nfptsgp_ibm,llpart_gp
     $        ,iptsgp_ibm,nigp_ibm,partl,0,rptsgp_ibm,nrgp_ibm,jgp_back)  ! jgpps is overwri

      ! sort by element for ghost particles for projection performance
      call fgslib_crystal_tuple_sort    (i_cr_hndl,nfptsgp_ibm
     $        ,iptsgp_ibm,nigp_ibm,partl,0,rptsgp_ibm,nrgp_ibm,jgpes,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine IBM_part_forcing
c      
c     sum up the forces from all markers across different processors.
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      ! debug
      if(ibm_debug_bin.eq.1) then
      if(nfptsgp_ibm .gt.0) then
!         print*,"Ghost Queen rank",nid, "receives",nfptsgp_ibm
!         do i = 1, nfptsgp_ibm
!            write(6,'(A,I6,I6,A,3F12.6)')
!     $           "Ghost Queen returned",nid,i," force = ",
!     $          rptsgp_ibm(jgpfh,i)
!     $         ,rptsgp_ibm(jgpfh+1,i),rptsgp_ibm(jgpfh+2,i)
!         enddo       
      endif
      endif

      do i_qt = 1,n
         if(ipart(jrole,i_qt) .ne. 1) cycle ! only for Queens
         i_qt_pid1  = ipart(jpid1, i_qt)
         i_qt_pid3  = ipart(jpid3, i_qt)
         i_qt_queen = ipart(jqueen,i_qt)
         do i_gp = 1, nfptsgp_ibm ! sum ghost queens     
            i_gp_pid1  = iptsgp_ibm(jgp_pid1, i_gp)
            i_gp_pid2  = iptsgp_ibm(jgp_pid2, i_gp)
            i_gp_pid3  = iptsgp_ibm(jgp_pid3, i_gp)
            i_gp_role  = iptsgp_ibm(jgp_role, i_gp)
            i_gp_queen = iptsgp_ibm(jgp_queen,i_gp)
            if (i_gp_pid1 .eq. i_qt_pid1 .and.
     $          i_gp_pid3 .eq. i_qt_pid3 .and.
     $          i_gp_queen.eq. i_qt_queen ) then               
               do j=0,ndim-1
                  rpart(jf0+j,i_qt) = rpart(jf0   + j, i_qt) +
     $                 rptsgp_ibm (jgpfh + j, i_gp)      
!                  ipart(jai, i_qt) =  ipart(jai, i_qt) +
!     $                 iptsgp_ibm (jgp_ax, i_gp)
                  if(ibm_rotation.eq.1) then
                     rpart(jtorque0 + j, i_qt) =
     $               rpart(jtorque0 + j, i_qt) +
     $               rptsgp_ibm(jgp_torque0 + j, i_gp)      
                  endif
               enddo
            endif ! match
          enddo ! i_gp 


!       if( ipart(jai, i_qt) .lt. ipart(jnlm, i_qt) )
!     $ call exitti('Error in Number of workers in :$',nid)

      enddo     ! i_qt

       ! debug
      if(ibm_debug_queen.eq.1) then         
          if(istep.eq.1.and.nid.eq.0) then 
             write(6,'(A)')
     $      "QueenForce variables=istep,nid,i_qt,b_pid,b_jqueen,
     $ Fx,Fy,Fz,Tx,Ty,Tz"
          endif
          do i_qt =1,n
          if(ipart(jrole,i_qt) .ne. 1) cycle
          !          write(6,'(A,5I7,12E12.4)') "QueenForce",
!     $        istep, nid, i_qt, ipart(jpid1,i_qt),ipart(jqueen,i_qt)
!     $     , (rpart(jf0+j,i_qt),j=0,2)
!     $     , (rpart(jtorque0+j,i_qt),j=0,2)
          enddo
       endif
      
       return
       end
c-----------------------------------------------------------------
      subroutine IBM_part_integrate
c      
c     sum up the forces from all markers across different processors.
c      
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage
      real                  tcoef(3,3),dt_cmt,time_cmt
      common /timestepcoef/ tcoef,dt_cmt,time_cmt

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      common /PARTRK3/ kv_stage_p
      real   kv_stage_p(llpart,13) ! add 6 for rotation

      integer fdim
      real    pmass
      real rforce_local(0:2), rvloc_queen(0:2)

      jx0 = jx
      rpi = 4.0*atan(1.0)
            
c     rk3 stage one items ---------------------------------------------
      if (stage.eq.1) then
c        used for time derivative of v in iu force
         call get_bdf_ext_coefs(beta,alpha,times)
         do i=1,n
            if(ipart(jrole,i) .ne. 1) cycle
            kv_stage_p(i,1) = rpart(jx0  ,i)
            kv_stage_p(i,2) = rpart(jx0+1,i)
            kv_stage_p(i,3) = rpart(jx0+2,i)
            kv_stage_p(i,4) = rpart(jv0  ,i)
            kv_stage_p(i,5) = rpart(jv0+1,i)
            kv_stage_p(i,6) = rpart(jv0+2,i)
            kv_stage_p(i,7) = rpart(jtemp,i)
            if(ibm_rotation.eq.1) then
               kv_stage_p(i,8)  = rpart(jangvel0,i)
               kv_stage_p(i,9)  = rpart(jangvel0+1,i)
               kv_stage_p(i,10) = rpart(jangvel0+2,i)
               kv_stage_p(i,11) = rpart(jangle0  ,i)
               kv_stage_p(i,12) = rpart(jangle0+1,i)
               kv_stage_p(i,13) = rpart(jangle0+2,i)
            endif
         enddo
      endif ! if(stage.eq.1)

      rforce_local = 0.0      
      rvloc_queen = 0.0
      nptot0 = 0
c     move data to previous positions, note, differnt reprentation for jv1, jv2, jv3 represent stages, in bdf jv1,jv2 jv3, represent t^{n-1} m t^{n-2}, t^{n-3}
      do i=1,n
         do j=0,ndim-1
            
            if(ipart(jrole,i) .ne. 1) cycle
            rpart(ju3+j,i)=rpart(ju2+j,i)
            rpart(ju2+j,i)=rpart(ju1+j,i)
            rpart(ju1+j,i)=rpart(ju0+j,i)
            rpart(jv3+j,i)=rpart(jv2+j,i)
            rpart(jv2+j,i)=rpart(jv1+j,i)
            rpart(jv1+j,i)=rpart(jv0+j,i)
            rpart(jx3+j,i)=rpart(jx2+j,i)
            rpart(jx2+j,i)=rpart(jx1+j,i)
            rpart(jx1+j,i)=rpart(jx0+j,i)
            if(ibm_rotation.eq.1) then
               rpart(jangvel3+j,i)=rpart(jangvel2+j,i)
               rpart(jangvel2+j,i)=rpart(jangvel1+j,i)
               rpart(jangvel1+j,i)=rpart(jangvel0+j,i)
               rpart(jangle3+j,i)=rpart(jangle2+j,i)  ! 3-4
               rpart(jangle2+j,i)=rpart(jangle1+j,i)
               rpart(jangle1+j,i)=rpart(jangle0+j,i)
            endif
            
             rforce_local(j) = rforce_local(j) + rpart(jf0+j,i)
             rvloc_queen(j)  = rvloc_queen(j)  + rpart(jv0+j,i)
             nptot0 = nptot0 + 1
          enddo
      enddo

      !!! check force on Queens
      nptot01      = iglsum(nptot0,1)
      dragsum_tot0 = glsum(rforce_local(0),1)
      dragsum_tot1 = glsum(rforce_local(1),1)
      dragsum_tot2 = glsum(rforce_local(2),1)
      rvx_avg_queen = glsum(rvloc_queen(0),1)/nptot01
      rvy_avg_queen = glsum(rvloc_queen(1),1)/nptot01
      rvz_avg_queen = glsum(rvloc_queen(2),1)/nptot01
      if(nid.eq.0) write(6,'(a,I6,20e12.4)') "HydroForce Queen",
     $     nptot01,dragsum_tot0, dragsum_tot1, dragsum_tot2
     $     ,rvx_avg_queen, rvy_avg_queen, rvz_avg_queen
      
      ! reassmble lhs and rhs terms in the equations
      do i=1,n
         if(ipart(jrole,i) .ne. 1) cycle
         do j=0,ndim-1
            ! rmp_ibm = rpart(jvol, i) * rpart(jrhop,i) ! mass fraction, wrong
            rvol0  = rpart(jvol, i)
            rho_p0 = rpart(jrhop,i)
            rho_f0 = rpart(jrho, i)

            ! Method 1,! Uhlmann 2005, Eq B.3, fluctuation when rho_p/rho_f < 1.2, for single stage
            rmp_ibm = rvol0 *(rho_p0 - rho_f0) 
            rpart(jf0+j,i) = rpart(jf0+j,i)/rmp_ibm
            
            ! Method2, the added mass effect is assumed from last step
            !rvpc = rho_f0 * rvol0 * (rpart(jv1+j,i) - rpart(jv2+j,i))/dt
            !rpart(jf0+j,i) = rpart(jf0+j,i) + rvpc            
            !rmp_ibm        = rvol0 * rho_p0
            !rpart(jf0+j,i) = rpart(jf0+j,i)/rmp_ibm

            if( ibm_rotation.eq.1  .and. non_spherical.eq.0) then
               r = rpart(jdp,i)/2
               rmom_inertia = (2.0/5.0) * rho_p0*rvol0 * (r**2) ! moment of inertia
               ! Method1 add rates of change of angular momentum, Uhlmann 2005, Eq. (B.4)
               romegac  = (rpart(jangvel1+j,i) - rpart(jangvel2+j,i))/dt
               ! Method2, explicit calculate rate of change of angular momentum
               r_ang_mom = rmom_inertia / rho_p0 * romegac
               rpart(jtorque0+j,i) = rpart(jtorque0+j,i) + r_ang_mom
               rpart(jtorque0+j,i) = rpart(jtorque0+j,i)/rmom_inertia
            endif

         enddo
      enddo

      ! add 6/5/2019 for ellipsoid moment of inertia
      if(non_spherical.eq.1) then ! ellipsoid / Spheroid
         do i=1,n
            if(ipart(jrole,i) .ne. 1) cycle 
            r_1   = rpart(jdp,i)/2      ! a
            eccentricity = rpart(jar,i) ! defined in the auxillary var, 1.5
            r_2   = r_1/eccentricity    ! b
            r_3   = r_1/eccentricity    ! c
            rrho_1 = rpart(jrhop,i)     ! density
            rm_1   = rrho_1 * 4/3 * rpi * r_1 * r_2 * r_3 ! mass

            ! In the principal axis of the body fitted frame 
            rIxx = rm_1 / 5 * (r_2**2 + r_3**2)
            rIxy = 0.0
            rIxz = 0.0
            rIyx = 0.0
            rIyy = rm_1 / 5 * (r_1**2 + r_3**2)
            rIyz = 0.0
            rIzx = 0.0
            rIzy = 0.0
            rIzz = rm_1 / 5 * (r_1**2 + r_2**2)
            
            do j=0,ndim-1
               ! add rates of change of angular momentum, Uhlmann 2005, Eq. (B.4)
               romegac  = (rpart(jangvel1+j,i) - rpart(jangvel2+j,i))/dt ! Method2,explicit calculate rate of change of angular momentum
               r_ang_mom= rIzz / rrho_1 * romegac
               rpart(jtorque0+j,i) = rpart(jtorque0+j,i) + r_ang_mom

               rpart(jtorque0+j,i) = rpart(jtorque0+j,i)/rIzz ! temporally for 2D rotation
            enddo
            if(istep.le.10) print*,"moment of inertia:"
     >           , rIzz, (rpart(jtorque0+j,i),j=0,2)

         enddo
      endif


      if(ipart_moving .eq. 1) then
c     all rk3 stages items --------------------------------------------
      do i=1,n
         if(ipart(jrole,i) .ne. 1) cycle  !!! for IBM Queen only
         rpart(jv0  ,i) = tcoef(1,stage)*kv_stage_p(i,4) +
     >                    tcoef(2,stage)*rpart(jv0  ,i)  +
     >                    tcoef(3,stage)*rpart(jf0  ,i)
         rpart(jv0+1,i) = tcoef(1,stage)*kv_stage_p(i,5) +
     >                    tcoef(2,stage)*rpart(jv0+1,i)  +
     >                    tcoef(3,stage)*rpart(jf0+1,i)
         rpart(jv0+2,i) = tcoef(1,stage)*kv_stage_p(i,6) +
     >                    tcoef(2,stage)*rpart(jv0+2,i)  +
     >                    tcoef(3,stage)*rpart(jf0+2,i)
         rpart(jx0  ,i) = tcoef(1,stage)*kv_stage_p(i,1) +
     >                    tcoef(2,stage)*rpart(jx0  ,i)  +
     >                    tcoef(3,stage)*rpart(jv0  ,i)
         rpart(jx0+1,i) = tcoef(1,stage)*kv_stage_p(i,2) +
     >                    tcoef(2,stage)*rpart(jx0+1,i)  +
     >                    tcoef(3,stage)*rpart(jv0+1,i)
         if (if3d)
     >   rpart(jx0+2,i) = tcoef(1,stage)*kv_stage_p(i,3) +
     >                    tcoef(2,stage)*rpart(jx0+2,i)  +
     >                    tcoef(3,stage)*rpart(jv0+2,i)
         rpart(jtemp,i) = tcoef(1,stage)*kv_stage_p(i,7) +
     >                    tcoef(2,stage)*rpart(jtemp,i)  +
     >                    tcoef(3,stage)*rpart(jq0  ,i)

         if(ibm_rotation.eq.1) then
            rpart(jangvel0  ,i) = tcoef(1,stage) * kv_stage_p(i,  8) +   !!! update ang velocity 
     >                            tcoef(2,stage) * rpart(jangvel0,i) +
     >                            tcoef(3,stage) * rpart(jtorque0,i)
            rpart(jangvel0+1,i) = tcoef(1,stage) * kv_stage_p(i,  9)  +
     >                            tcoef(2,stage) * rpart(jangvel0+1,i)+
     >                            tcoef(3,stage) * rpart(jtorque0+1,i)
            rpart(jangvel0+2,i) = tcoef(1,stage) * kv_stage_p(i, 10)  +
     >                            tcoef(2,stage) * rpart(jangvel0+2,i)+
     >                            tcoef(3,stage) * rpart(jtorque0+2,i)
            rpart(jangle0  ,i)  = tcoef(1,stage) * kv_stage_p(i,  11) +   !!! update angle
     >                            tcoef(2,stage) * rpart(jangle0  ,i) +
     >                            tcoef(3,stage) * rpart(jangvel0  ,i)
            rpart(jangle0+1,i)  = tcoef(1,stage) * kv_stage_p(i,  12) +
     >                            tcoef(2,stage) * rpart(jangle0+1,i) +
     >                            tcoef(3,stage) * rpart(jangvel0+1,i)
            rpart(jangle0+2,i)  = tcoef(1,stage) * kv_stage_p(i,  13) +
     >                            tcoef(2,stage) * rpart(jangle0+2,i) +
     >                            tcoef(3,stage) * rpart(jangvel0+2,i)
         endif

      enddo

      if(ibm_translation .eq. -1) then ! disable translation
         do i=1,n
            if(ipart(jrole,i) .ne. 1) cycle
            rpart(jx0  ,i) = kv_stage_p(i,1) 
            rpart(jx0+1,i) = kv_stage_p(i,2) 
            rpart(jx0+2,i) = kv_stage_p(i,3) 
            rpart(jv0  ,i) = kv_stage_p(i,4) 
            rpart(jv0+1,i) = kv_stage_p(i,5) 
            rpart(jv0+2,i) = kv_stage_p(i,6) 
         enddo
      endif

      elseif(ipart_moving .eq. 2) then
         call ibm_part_force_moving_ibm
      else
         if(istep.lt.2 .and.nid.eq.0) print*, "Stationary Particle"
      endif

      
      ! debug
      if(ibm_debug.eq.1) then
      
c         open(unit=1122,file='Particle_hist.dat')
      !    if(istep.eq.1.and.stage.eq.1.and.nid.eq.0) then
      !     write(1122,*)"istep,nid,stage,i_qt,pid1,pid2
      !,Vx,Vy,Vz,x,y,z,Fx/m,Fy/m,Fz/m
      !    endif
         do i =1,n
            if(ipart(jrole,i).ne.1) cycle

            do j = 0,2
               vintdiff_temp = abs(rpart(jv1+j,i)-rpart(jv0+j,i))
     $          / (abs(rpart(jv1+j,i)) + 1.0e-6)
               if( vintdiff_temp .ge. 0.2 .and.
     $             abs(rpart(jv1+j,i)) .ge. 1.0e-2    )
     $              write(6,'(A,7I7,21E12.4)') "WarnVel"
     $              ,istep, stage,nid, i, j
     $              ,ipart(jpid1,i),ipart(jpid2,i)
     $              ,vintdiff_temp
     $              ,rpart(jv1+j,i), rpart(jv0+j,i), rpart(jf0+j,i)
     $              ,rpart(jx0+0,i), rpart(jx0+1,i), rpart(jx0+2,i)
            enddo
 
            if(ibm_rotation.eq.1) then
               write(6,'(A,6I7,21E12.4)')'Queent ',
     >       istep,nid,stage,i,ipart(jpid1,i),ipart(jpid2,i)
     >      ,(rpart(jv0+j,i),j=0,2)
     >      ,(rpart(jx0+j,i),j=0,2)
     >      ,(rpart(jf0+j,i),j=0,2)
     >      ,(rpart(jangvel0+j,i),j=0,2)
     >      ,(rpart(jangle0+j,i), j=0,2)
     >      ,(rpart(jtorque0+j,i),j=0,2)
            else
               write(6,'(A,6I7,21E12.4)')'Queent ',
     >       istep,nid,stage,i,ipart(jpid1,i),ipart(jpid2,i)
     >      ,(rpart(jv0+j,i),j=0,2)
     >      ,(rpart(jx0+j,i),j=0,2)
     >      ,(rpart(jf0+j,i),j=0,2)
            endif

         enddo
      endif

      return
      end
            
c-----------------------------------------------------------------
      subroutine IBM_local_worker_update

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer iptmove !(llpart)

      real rx_l, ry_l, rz_l
      real omg_x, omg_y, omg_z
      real u_rot, v_rot, w_rot
      real dtheta(3)
      
      jx0 = jx
      do i = 1,n

         iptmove = 0            ! flag if the current marker is updated!

         if(ipart(jrole,i) .ne. 2) cycle ! only for local Worker

         ipt_pid1  = ipart(jpid1, i)
         ipt_pid3  = ipart(jpid3, i)
         ipt_queen = ipart(jqueen,i)

         do i_qt = 1,n
            if(ipart(jrole,i_qt) .ne. 1) cycle ! from local Queen

            ipt_q_pid1  = ipart(jpid1,i_qt) 
            ipt_q_pid3  = ipart(jpid3,i_qt) 
            ipt_q_queen = ipart(jqueen,i_qt) 

            if (ipt_q_pid1 .eq. ipt_pid1 .and.
     $          ipt_q_pid3 .eq. ipt_pid3 .and.
     $          ipt_q_queen.eq. ipt_queen ) then               

               do j=0,ndim-1
                   ! translate 
                  rpart(jv0+j,i) = rpart(jv0+j,i_qt) ! center translate velocity  
                  rpart(jx+j,i)  = rpart(jx+j,i) +
     $                 ( rpart(jx+j,i_qt) - rpart(jx1+j,i_qt) )
               enddo
                   ! rotation
                   if(ibm_rotation.eq.1) then
                      
                      rx_l = rpart(jx,  i) - rpart(jx,   i_qt)  ! length vecotor r from center
                      ry_l = rpart(jx+1,i) - rpart(jx+1, i_qt)
                      rz_l = rpart(jx+2,i) - rpart(jx+2, i_qt)

                      omg_x = rpart(jangvel0,   i_qt) ! angular velocity, center
                      omg_y = rpart(jangvel0+1, i_qt)
                      omg_z = rpart(jangvel0+2, i_qt)                      
     
                      u_rot = omg_y * rz_l - omg_z * ry_l
                      v_rot = omg_z * rx_l - omg_x * rz_l
                      w_rot = omg_x * ry_l - omg_y * rx_l

                      ! (1) update velocity
                      rpart(jv0,  i) = rpart(jv0,  i) + u_rot           
                      rpart(jv0+1,i) = rpart(jv0+1,i) + v_rot
                      rpart(jv0+2,i) = rpart(jv0+2,i) + w_rot
                      
                      ! (2) update angular velocity
                      rpart(jangvel0,   i) = rpart(jangvel0  ,   i_qt)  
                      rpart(jangvel0+1, i) = rpart(jangvel0+1,   i_qt)
                      rpart(jangvel0+2, i) = rpart(jangvel0+2,   i_qt)

                      ! rotating angle
                      dtheta(1)=rpart(jangle0,  i_qt)-rpart(jangle0  ,i)  
                      dtheta(2)=rpart(jangle0+1,i_qt)-rpart(jangle0+1,i)
                      dtheta(3)=rpart(jangle0+2,i_qt)-rpart(jangle0+2,i)

                      ! (3) update sufrace point position
                      call rotate3D(rpart(jx,i),rpart(jx,i_qt),dtheta)

                      ! (4) update angular position
                      rpart(jangle0  , i) = rpart(jangle0  , i_qt) 
                      rpart(jangle0+1, i) = rpart(jangle0+1, i_qt)
                      rpart(jangle0+2, i) = rpart(jangle0+2, i_qt)
                   endif
                                      
                   iptmove = 1  ! update flag

                   if(ibm_debug_worker.eq.1 .and. i.eq.1) then
                      if(ibm_rotation .eq. 1) then
                         write(6,'(A,4I5,4(A,3E12.4))')
     >                       "local Worker t=", istep,nid
     >                     , ipart(jpid1,i),ipart(jpid2,i)
     >                     , " Vel=", (rpart(jv0+j,i),j=0,2)
     >                     , " Pos=", (rpart(jx+j,i),j=0,2)
     >                     , " Angvel=",(rpart(jangvel0+j,i),j=0,2) 
     >                     , " Angle=",(rpart(jangle0+j,i),j=0,2)
                      else
                         write(6,'(A,4I5,4(A,3E12.4))')
     >                       "local Worker t=", istep,nid
     >                     , ipart(jpid1,i),ipart(jpid2,i)
     >                     , " Vel=", (rpart(jv0+j,i),j=0,2)
     >                     , " Pos=", (rpart(jx+j,i),j=0,2)

                      endif                     
                   endif
                   
                   goto 2143
                endif
             enddo

             do i_gp = 1, nfptsgp_ibm          ! from remote Queen 

                i_gp_pid1  = iptsgp_ibm(jgp_pid1, i_gp) ! born rank
                i_gp_pid2  = iptsgp_ibm(jgp_pid2, i_gp) 
                i_gp_pid3  = iptsgp_ibm(jgp_pid3, i_gp) ! born time
                i_gp_role  = iptsgp_ibm(jgp_role, i_gp)
                i_gp_queen = iptsgp_ibm(jgp_queen,i_gp) ! born queen id

                if (i_gp_pid1 .eq. ipt_pid1 .and.
     $              i_gp_pid3 .eq. ipt_pid3 .and.
     $              i_gp_queen.eq. ipt_queen ) then

                   do j=0,ndim-1       
                     ! translate velocity and position 
                     rpart(jv0+j,i) = rptsgp_ibm(jgpv0+j,i_gp) 
                     rpart(jx+j,i) = rpart(jx + j, i) +
     $              (rptsgp_ibm(jgpx+j,i_gp)-rptsgp_ibm(jgp_jx1+j,i_gp)) 
                  enddo
                  
                   ! rotation
                   if(ibm_rotation.eq.1) then
                      
                      rx_l = rpart(jx,  i) - rptsgp_ibm(jgpx, i_gp)  ! length vecotor r from center
                      ry_l = rpart(jx+1,i) - rptsgp_ibm(jgpx+1, i_gp)
                      rz_l = rpart(jx+2,i) - rptsgp_ibm(jgpx+2, i_gp)
                      
                      omg_x = rptsgp_ibm(jgp_angvel0 ,    i_gp) ! obtain angular velocity
                      omg_y = rptsgp_ibm(jgp_angvel0 + 1, i_gp) 
                      omg_z = rptsgp_ibm(jgp_angvel0 + 2, i_gp) 

                      u_rot = omg_y * rz_l - omg_z * ry_l
                      v_rot = omg_z * rx_l - omg_x * rz_l
                      w_rot = omg_x * ry_l - omg_y * rx_l

                      ! (1) update velocity
                      rpart(jv0,  i) = rpart(jv0,  i) + u_rot           
                      rpart(jv0+1,i) = rpart(jv0+1,i) + v_rot
                      rpart(jv0+2,i) = rpart(jv0+2,i) + w_rot
                      
                      ! (2) update angular velocity
                      rpart(jangvel0, i) =rptsgp_ibm(jgp_angvel0  ,i_gp) 
                      rpart(jangvel0+1,i)=rptsgp_ibm(jgp_angvel0+1,i_gp)
                      rpart(jangvel0+2,i)=rptsgp_ibm(jgp_angvel0+2,i_gp)

                      ! rotating angle
                      dtheta(1)=rptsgp_ibm(jgp_angle0,     i_gp)
     $                             - rpart(jangle0   ,  i)  
                      dtheta(2)=rptsgp_ibm(jgp_angle0 + 1, i_gp)
     $                             - rpart(jangle0 + 1 ,i)  
                      dtheta(3)=rptsgp_ibm(jgp_angle0 + 2, i_gp)
     $                             - rpart(jangle0 + 2, i)  
                      
                      ! (3) update sufrace point position
                      call rotate3D(rpart(jx,i),rptsgp_ibm(jgpx, i_gp)
     $                         ,dtheta)

                      ! (4) update angular position
                      rpart(jangle0  ,i)=rptsgp_ibm(jgp_angle0  , i_gp)
                      rpart(jangle0+1,i)=rptsgp_ibm(jgp_angle0+1, i_gp)
                      rpart(jangle0+2,i)=rptsgp_ibm(jgp_angle0+2, i_gp)
                      
                   endif
                  
                   ! update flag
                   iptmove = 1

                   if(ibm_debug_worker.eq.1.and. i.eq.1) then
                     if(ibm_rotation .eq. 1) then
                       write(6,'(A,4I5,4(A,3E12.4))')
     >                       "remote Worker t=", istep,nid
     >                     , ipart(jpid1,i),ipart(jpid2,i)
     >                     , " Vel=", (rpart(jv0+j,i),j=0,2)
     >                     , " Pos=", (rpart(jx+j,i),j=0,2)
     >                     , " Angvel=",(rpart(jangvel0+j,i),j=0,2) 
     >                   , " Angle=",(rpart(jangle0+j,i),j=0,2)
                     else
                       write(6,'(A,4I5,4(A,3E12.4))')
     >                       "remote Worker t=", istep,nid
     >                     , ipart(jpid1,i),ipart(jpid2,i)
     >                     , " Vel=", (rpart(jv0+j,i),j=0,2)
     >                     , " Pos=", (rpart(jx+j,i),j=0,2)                   
                     endif
                  endif
              
                 goto 2143

                endif 

             enddo
             
 2143        if(iptmove.eq.1) continue

         enddo
 
      return
      end

c-----------------------------------------------------------------
      subroutine IBM_ghost_worker_update

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer iptmove !(llpart)

      jx0 = jx

      do j=0,ndim-1
         do i=1,n

            iptmove = 0

            if(ipart(jrole,i) .ne. 2) cycle

             ipt_pid1  = ipart(jpid1, i)
             ipt_pid3  = ipart(jpid3, i)
             ipt_queen = ipart(jqueen,i)

             do i_gp = 1, nfptsgp_ibm

                i_gp_pid1  = iptsgp_ibm(jgp_pid1, i_gp)
                i_gp_pid2  = iptsgp_ibm(jgp_pid2, i_gp)
                i_gp_pid3  = iptsgp_ibm(jgp_pid3, i_gp)
                i_gp_role  = iptsgp_ibm(jgp_role, i_gp)
                i_gp_queen = iptsgp_ibm(jgp_queen,i_gp)

                if (i_gp_pid1 .eq. ipt_pid1 .and.
     $              i_gp_pid3 .eq. ipt_pid3 .and.
     $              i_gp_queen.eq. ipt_queen ) then               
                   
                   rpart(jv0+j,i) = rptsgp_ibm(jgpv0+j,i_gp) ! + w*r
                   rpart(jx+j,i) = rpart(jx + j, i) +
     $            (rptsgp_ibm(jgpx+j,i_gp)-rptsgp_ibm(jgp_jx1+j,i_gp)) ! method1

                   iptmove = 1                   
                   if(i.eq.1 .and. j.eq.1 ) write(6,*)
     >     "Rank=",    pid,"Remote Worker yVel=",rpart(jv0+1,i)
     >    ,"yPos=",   rpart(jx+1,i)
     >    ,"yForce=", rpart(jf0+1,i)
                   goto 2144
                endif 

             enddo

 2144        if(iptmove.eq.1) continue

          enddo
       enddo
          
      return
      end

      
      
c-----------------------------------------------------------------------
      subroutine usr_particles_f_col_IBM(i)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real mcfac,mcfac_wall, rdum, rdum3(3), rdum4(3)

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      real rpx1(3), rpx2(3), rpx0(3),r1(3),r2(3),r3(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      ! friction force and 
      real romega1(3), romega2(3), rft(3)
      common /rel_ang_vel/romega1,romega2,rft

      ! periodic collision check
      integer icx1,icy1,icz1,icx2,icy2,icz2
      common /col_part1/ icx1,icy1,icz1
      common /col_part2/ rx2, icx2,icy2,icz2

      jx0 = jx
      if(ipart(jrole,i).ne.1) return ! for IBM queen only

      ptdum(19) = dnekclock()
      
      mcfac       = 2.*sqrt(ksp)*log(e_rest)/
     >              sqrt(log(e_rest)**2+pi**2)
      mcfac_wall  = 2.*sqrt(ksp_wall)*log(e_rest_wall)/
     >              sqrt(log(e_rest_wall)**2+pi**2)

      if (two_way .gt. 2) then

      rrp1   = rpart(jrpe ,i) ! diameter
      rvol1  = rpart(jvol ,i) ! volume
      rrho1  = rpart(jrhop,i) ! density
      rx1(1) = rpart(jx   ,i) ! x
      rx1(2) = rpart(jx+1 ,i) ! y
      rx1(3) = rpart(jx+2 ,i) ! z
      rv1(1) = rpart(jv0  ,i) ! vx
      rv1(2) = rpart(jv0+1,i) ! vy
      rv1(3) = rpart(jv0+2,i) ! vz

      romega1(1) = rpart(jangvel0+0, i) ! omega1_x
      romega1(2) = rpart(jangvel0+1, i) ! omega1_y
      romega1(3) = rpart(jangvel0+2, i) ! omega1_z
      
      icx1 = ipart(jicx,i)    ! bin ii
      icy1 = ipart(jicy,i)    ! bin jj
      
      icz1 = 0
      if (if3d) icz1 = ipart(jicz,i) ! bin kk

      icxm = icx1 -1 ! neighboring bin ii-1
      icxp = icx1 +1 ! ii+1
      icym = icy1 -1 ! neighboring bin jj-1
      icyp = icy1 +1 ! jj+1

      iczm = 0
      iczp = 0
      if (if3d) then
         iczm = icz1 -1
         iczp = icz1 +1
      endif

c     let every particle search for itself
c        particles in local elements
         do j = i+1,n

            if(ipart(jrole,j).ne.1) cycle ! for queen marker only
            
            icx2 = ipart(jicx,j) ! residing bin i,j,k
            icy2 = ipart(jicy,j) 
            icz2 = ipart(jicz,j)

            rx2(1) = rpart(jx   ,j) ! x
            rx2(2) = rpart(jx+1 ,j) ! y
            rx2(3) = rpart(jx+2 ,j) ! z
            call periodic_collide_check
            !print*, "icx1 icx2",icx1,icy1,icz1,icx2,icy2,icz2,rx1,rx2
            if ((icx2.ge.icxm).and.(icx2.le.icxp)) then ! if the other Queen marker is neighboring
            if ((icy2.ge.icym).and.(icy2.le.icyp)) then 
            if ((icz2.ge.iczm).and.(icz2.le.iczp)) then 

            rrp2   = rpart(jrpe ,j) ! diameter
            rvol2  = rpart(jvol ,j) ! volume
            rrho2  = rpart(jrhop,j) ! density
            rv2(1) = rpart(jv0  ,j) ! vx
            rv2(2) = rpart(jv0+1,j) ! vy
            rv2(3) = rpart(jv0+2,j) ! vz

            romega2(1) = rpart(jangvel0+0, j) ! omega2_x
            romega2(2) = rpart(jangvel0+1, j) ! omega2_y
            romega2(3) = rpart(jangvel0+2, j) ! omega2_z
            
            idum = 1

            !print*, "percol",icx1,icy1,icz1,icx2,icy2,icz2,rx1,rx2
            call compute_collide_IBM(mcfac,rrp2,rvol2,rrho2,rx2,rv2
     >            , rpart(jfcol,i),  rpart(jfcol,j)
     >            , rpart(jtq_col,i),rpart(jtq_col,j)
     >            , idum)

            endif
            endif
            endif
         enddo

c        search list of ghost particles
         do j = 1,nfptsgp_IBM
            
            icx2 = iptsgp_IBM(jgpicx,j)   
            icy2 = iptsgp_IBM(jgpicy,j)   
            icz2 = iptsgp_IBM(jgpicz,j)
            rx2(1) = rptsgp_IBM(jgpx   ,j) ! x
            rx2(2) = rptsgp_IBM(jgpx+1 ,j) ! y
            rx2(3) = rptsgp_IBM(jgpx+2 ,j) ! z
            call periodic_collide_check
            !print*, "icx1 icx2",icx1,icy1,icz1,icx2,icy2,icz2,rx1,rx2
            if ((icx2.ge.icxm).and.(icx2.le.icxp)) then !if the Ghost Queen marker is neighboring
            if ((icy2.ge.icym).and.(icy2.le.icyp)) then
            if ((icz2.ge.iczm).and.(icz2.le.iczp)) then

            rrp2   = rptsgp_IBM(jgprpe ,j) 
            rvol2  = rptsgp_IBM(jgpvol ,j)
            rrho2  = rho_p                 ! assume same density. Need2fix

            rv2(1) = rptsgp_IBM(jgpv0  ,j)
            rv2(2) = rptsgp_IBM(jgpv0+1,j)
            rv2(3) = rptsgp_IBM(jgpv0+2,j)

            romega2(1) = rpart(jgp_angvel0 + 0, j)
            romega2(2) = rpart(jgp_angvel0 + 1, j)
            romega2(3) = rpart(jgp_angvel0 + 2, j)
            
            idum = 0

            !print*, "percol",icx1,icy1,icz1,icx2,icy2,icz2,rx1,rx2
            call compute_collide_IBM(mcfac,rrp2,rvol2,rrho2,rx2,rv2
     >            , rpart(jfcol,i),   rdum3
     >            , rpart(jtq_col,i), rdum4
     >            , idum )
            
            endif
            endif
            endif
         enddo
 1235 continue

         ! plane wall collisions
         do j = 1,np_walls
            rnx = plane_wall_coords(1,j)
            rny = plane_wall_coords(2,j)
            rnz = plane_wall_coords(3,j)
            rpx = plane_wall_coords(4,j)
            rpy = plane_wall_coords(5,j)
            rpz = 1.0
            if (if3d) rpz = plane_wall_coords(6,j)


            rd    = -(rnx*rpx + rny*rpy + rnz*rpz)

            rdist = abs(rnx*rpart(jx,i)+rny*rpart(jy,i)+rnz*rpart(jz,i)
     >                    +rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            rrp2   = 0.
            rvol2  = 1.
            rrho2  = 1E8
            rx2(1) = rpart(jx  ,i) - rdist*rnx ! distance x
            rx2(2) = rpart(jx+1,i) - rdist*rny ! distance y
            rx2(3) = rpart(jx+2,i) - rdist*rnz ! dist z
            rv2(1) = 0.
            rv2(2) = 0.
            rv2(3) = 0.

            romega2(1) = 0.0
            romega2(2) = 0.0
            romega2(3) = 0.0

            idum = 0
            call periodic_collide_check
            call compute_collide_IBM(mcfac_wall,rrp2,rvol2,rrho2,rx2,rv2
     >            , rpart(jfcol,i),   rdum3
     >            , rpart(jtq_col,i), rdum4
     >            , idum )

         enddo

         ! cylinder wall collisions
         do j = 1,nc_walls
            rnx = cyl_wall_coords(1,j)
            rny = cyl_wall_coords(2,j)
            rnz = cyl_wall_coords(3,j)
            rpx = cyl_wall_coords(4,j)
            rpy = cyl_wall_coords(5,j)
            rpz = 1.0
            if (if3d) rpz = cyl_wall_coords(6,j)

            rrad = cyl_wall_coords(7,j)

            rx2(1) = rpart(jx,i)
            rx2(2) = rpart(jy,i)
            rx2(3) = rpart(jz,i)
            ! for now only works with cylinders aligned with axes at
            ! origin
            if (rnz .gt. 0.5) then
               rtheta = atan2(rpart(jy,i)-rpy,rpart(jx,i)-rpx)
               rx2(1) = rpx+rrad*cos(rtheta)
               rx2(2) = rpy+rrad*sin(rtheta)
            elseif (rnx .gt. 0.5) then
               rtheta = atan2(rpart(jz,i)-rpz,rpart(jy,i)-rpy)
               rx2(2) = rpy+rrad*cos(rtheta)
               rx2(3) = rpz+rrad*sin(rtheta)
            elseif (rny .gt. 0.5) then
               rtheta = atan2(rpart(jx,i)-rpx,rpart(jz,i)-rpz)
               rx2(3) = rpz+rrad*cos(rtheta)
               rx2(1) = rpx+rrad*sin(rtheta)
            endif

            rrp2   = 0.
            rvol2  = 1.
            rrho2  = 1E8
            rv2(1) = 0.
            rv2(2) = 0.
            rv2(3) = 0.

            romega2(1) = 0.0
            romega2(2) = 0.0
            romega2(3) = 0.0

            idum = 0
            call periodic_collide_check
            call compute_collide_IBM(mcfac_wall,rrp2,rvol2,rrho2,rx2,rv2
     >            , rpart(jfcol,i), rdum3
     >            , rpart(jtq_col,i),rdum4
     >            , idum)

         enddo

      endif


      pttime(19) = pttime(19) + dnekclock() - ptdum(19)

      return
      end

c-------------------------------------------------------------------------------
      subroutine sort_local_particles_collisions_IBM
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'LPM'

      ! set the residing bin ii,jj,kk for the Queen/Ghost Queen marker 
      ! distance is larger. We only need collsion distance as dpe_max

      ndxgpc_IBM = floor((xdrange(2,1) - xdrange(1,1))/d2chk_IBM(3))+1
      ndygpc_IBM = floor((xdrange(2,2) - xdrange(1,2))/d2chk_IBM(3))+1
      ndzgpc_IBM = 1
      if (if3d)
     $ ndzgpc_IBM =floor((xdrange(2,3) - xdrange(1,3))/d2chk_IBM(3))+1

      ! grid spacing for that many spacings
      rdxgpc_IBM = (xdrange(2,1) - xdrange(1,1))/real(ndxgpc_IBM)
      rdygpc_IBM = (xdrange(2,2) - xdrange(1,2))/real(ndygpc_IBM)
      rdzgpc_IBM = 1.
      if (if3d)
     $rdzgpc_IBM = (xdrange(2,3) - xdrange(1,3))/real(ndzgpc_IBM)

      ! set real particles ii,jj,kk
      do i=1,n
         if(ipart(jrole,i).ne.1) cycle ! for queen only
         rxval = rpart(jx,i)
         ryval = rpart(jy,i)
         rzval = 0.
         if(if3d) rzval = rpart(jz,i)
  
         ii    = floor((rxval-xdrange(1,1))/rdxgpc_IBM) 
         jj    = floor((ryval-xdrange(1,2))/rdygpc_IBM) 
         kk    = floor((rzval-xdrange(1,3))/rdzgpc_IBM) 
         ndum  = ii + ndxgp_IBM*jj + ndxgp_IBM*ndygp_IBM*kk

         ipart(jicx,i) = ii
         ipart(jicy,i) = jj
         ipart(jicz,i) = kk

      enddo

      ! set Ghost Queen particles ii,jj,kk
      do i=1,nfptsgp_IBM
         rxval = rptsgp_IBM(jgpx,i)
         ryval = rptsgp_IBM(jgpy,i)
         rzval = 0.
         if(if3d) rzval = rptsgp_IBM(jgpz,i)
  
         ii    = floor((rxval-xdrange(1,1))/rdxgpc_IBM) 
         jj    = floor((ryval-xdrange(1,2))/rdygpc_IBM) 
         kk    = floor((rzval-xdrange(1,3))/rdzgpc_IBM) 
         ndum  = ii + ndxgp_IBM*jj + ndxgp_IBM*ndygp_IBM*kk

         iptsgp_IBM(jgpicx,i) = ii
         iptsgp_IBM(jgpicy,i) = jj
         iptsgp_IBM(jgpicz,i) = kk

      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine ibm_f_marker(i)
c
c     ibm worker marker force
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
      
      pmassf= lpmvol_p * lpmdens_f

      lpmforce(1) = (lpmv_f(1) - lpmv_p(1)) / dt  * pmassf
      lpmforce(2) = (lpmv_f(2) - lpmv_p(2)) / dt  * pmassf
      lpmforce(3) = (lpmv_f(3) - lpmv_p(3)) / dt  * pmassf
      
      return
      end

c-----------------------------------------------------------------------
      subroutine ibm_part_force_moving_ibm
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      jx0 = jx
      
      do i=1,n
         if(ipart(jrole,i).ne.1) cycle 
         rpart(jv0  ,i) = rvpx_force ! forced velocity
         rpart(jv0+1,i) = rvpy_force 
         if(if3d) rpart(jv0+2,i) = rvpz_force
         rpart(jx  ,i) = rpart(jx  ,i) + rvpx_force*dt 
         rpart(jx+1,i) = rpart(jx+1,i) + rvpy_force*dt
         if (if3d)
     >   rpart(jx+2,i) = rpart(jx+2,i) + rvpz_force*dt         
c         if(ipart(jrole,i).eq.1) print*, "new_position",rpart(jx+1,i)
         if(ibm_rotation.eq.1) then
            rpart(jangvel0  ,i) = ravx_force ! forced angular velocity
            rpart(jangvel0+1,i) = ravy_force
            rpart(jangvel0+2,i) = ravz_force
            rpart(jangle0  ,i)  = rpart(jangle0  ,i) + ravx_force * dt
            rpart(jangle0+1,i)  = rpart(jangle0+1,i) + ravy_force * dt
            rpart(jangle0+2,i)  = rpart(jangle0+2,i) + ravz_force * dt
         endif
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine compute_collide_IBM(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                                     fcf1, fcf2, ftq1, ftq2, iflg)
c
c     extra body forces and torques, fcf1, fcf2, ftq1, ftq2
c
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real fcf1(3),fcf2(3),er,eta,ere,mcfac
      integer iflg ! iflg = 1, two real particle collides, both will be updated, otherwise only first one update

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      ! friction variables
      real rv12(3), rv12o(3), rv12t(3), rv12n(3), rvtr(3) ! total, tangential vectors
      real romega1(3), romega2(3), rt_12(3), rtorq_col(3), rft(3)
      common /rel_ang_vel/romega1,romega2, rft

      real c_fric, eta_t
      real ftq1(3), ftq2(3)

      rDelta_c = (rrp1+rrp2) / 20.0 ! a threshhold
      rthresh  = rrp1 + rrp2 + rDelta_c
      rthresh2 = rthresh**2

      rxdiff  = rx2(1) - rx1(1)
      rsum2 = rxdiff**2
      if (rsum2 .gt. rthresh2) goto 1513

      rydiff = rx2(2) - rx1(2)
      rydiff2 = rydiff**2
      rsum2 = rsum2 + rydiff2
      if (rsum2 .gt. rthresh2) goto 1513

      if (if3d) then
         rzdiff = rx2(3) - rx1(3)
         rzdiff2 = rzdiff**2
         rsum2 = rsum2 + rzdiff2
         if (rsum2 .gt. rthresh2) goto 1513
      endif

      rdiff = sqrt(rsum2)
      rm1   = rrho1*rvol1
      rm2   = rrho2*rvol2

      ! wall
      if (rm2 .gt. 1E7) then
         rksp_use = ksp_wall ! restitution coefficient
         re_rest_use = e_rest_wall
      ! other particles
      else
         rksp_use = ksp
         re_rest_use = e_rest
      endif

      rmult = 1./sqrt(1./rm1 + 1./rm2)
      eta  = mcfac*rmult

      ! first, handle normal collision part
      rbot     = 1./rdiff
      rn_12x   = rxdiff*rbot  ! normal unit vector 
      rn_12y   = rydiff*rbot  
      rn_12z   = rzdiff*rbot

      rdelta12 = rthresh - rdiff ! overlap distance, netative sign

      if(ibm_rotation.eq.1) then

         do j=1,3
            rv12(j) = rv2(j) - rv1(j)
         enddo         
         do j=1,3
            rvtr(j) = rrp1 * romega1(j) + rrp2 * romega2(j)
         enddo

         rv12o(1) = rv12(1) + ( rvtr(2) * rn_12z - rvtr(3) * rn_12y ) ! V = Omega X R
         rv12o(2) = rv12(2) + ( rvtr(3) * rn_12x - rvtr(1) * rn_12z ) ! V = Omega X R
         rv12o(3) = rv12(3) + ( rvtr(1) * rn_12y - rvtr(2) * rn_12x ) ! V = Omega X R

         rv12_mag =  rv12o(1) * rn_12x +
     >               rv12o(2) * rn_12y +
     >               rv12o(3) * rn_12z 
         
      else
        
         rv12_mag = (rv2(1) - rv1(1))*rn_12x +
     >              (rv2(2) - rv1(2))*rn_12y +
     >              (rv2(3) - rv1(3))*rn_12z

      endif

      rv12_mage = rv12_mag*eta      !damping force
      rksp_max = rksp_use*rdelta12  !spring force

      rnmag = -rksp_max - rv12_mage
      
      rfn1 = rnmag*rn_12x
      rfn2 = rnmag*rn_12y
      rfn3 = 0.0
      if (if3d) rfn3 = rnmag*rn_12z

      fcf1(1) = fcf1(1) + rfn1
      fcf1(2) = fcf1(2) + rfn2
      fcf1(3) = fcf1(3) + rfn3

      if (iflg .eq. 1) then
          fcf2(1) = fcf2(1) - rfn1
          fcf2(2) = fcf2(2) - rfn2
          fcf2(3) = fcf2(3) - rfn3
      endif

      ! tangential force

      if(ibm_rotation.eq.1) then
!         if(iflg.eq.0) then
!     $ rDelta_c = ibm_diam(1)/n_dh(1)/10.0     ! 1 compare with Delta_x  (Kidanemariam&Uhlmann 2014)
!            rDelta_c = 0.75*rv12_mag*dt ! 1 compare with Delta_x  (Capecelatro2013)
!         else
         rDelta_c =  0.5 * rv12_mag * dt
!         endif
         rdelta12 = abs(rdelta12) + abs(rDelta_c)    ! overlap width, add Delta_c as the admitted gap by Patankar2001

         rv12_mage = rv12_mag * eta        ! damping force
         rksp_max  = rksp_use * rdelta12   ! spring force         
         rnmag     = abs(rksp_max) + abs(rv12_mage) ! normal force

         rv12n(1) = rv12_mag * rn_12x
         rv12n(2) = rv12_mag * rn_12y
         rv12n(3) = rv12_mag * rn_12z

         do j=1,3
            rv12t(j) = rv12(j)-rv12n(j) ! tangential velocity vector
         enddo

         rv12t_mag = sqrt(rv12t(1)**2+rv12t(2)**2+rv12t(3)**2)
         if (rv12t_mag.le.0.00000001) goto 1513

         do j=1,3
            rv12t(j) = rv12t(j)/rv12t_mag ! unit tangential vector
         enddo

         eta_t         = 0.5*eta  ! tangential damping coefficient,navarro2013  
         c_fric        = 0.10 !0.092 ! friction coefficient 
         rtmag_term1   = c_fric*abs(rnmag)
         kt            = ksp*2.0/7.0 ! tangential spring constant,navarro2013 
         rdelta_t      = 0.0  ! tangential displacement, not known 
         rtmag_term2   = abs(eta_t*rv12t_mag + kt*rdelta_t)
         rtmag    = min(rtmag_term1, rtmag_term2) ! tangential force
         !rtmag     = rtmag_term1

         do j =1,3
            rft(j) = rtmag*rv12t(j)
         enddo

         do j = 1,3
            fcf1(j) = fcf1(j) + rft(j)
         enddo
         
         ! compute torque
         rtorq_col(1) = rrp1 * (rn_12y * rft(3) - rn_12z * rft(2))
         rtorq_col(2) = rrp1 * (rn_12z * rft(1) - rn_12x * rft(3))
         rtorq_col(3) = rrp1 * (rn_12x * rft(2) - rn_12y * rft(1))

         ! update torque
         do j=1,3
            ftq1(j) = ftq1(j) + rtorq_col(j)
         enddo

         if (iflg .eq. 1) then  ! colliding of two real particles
            do j=1,3
               fcf2(j) = fcf2(j) - rft(j)
            enddo
             ! compute torque for particle j
            rtorq_col(1) = rrp2 * (rn_12y * rft(3) - rn_12z * rft(2))
            rtorq_col(2) = rrp2 * (rn_12z * rft(1) - rn_12x * rft(3))
            rtorq_col(3) = rrp2 * (rn_12x * rft(2) - rn_12y * rft(1))
            do j=1,3
               ftq2(j) = ftq2(j) - rtorq_col(j)
            enddo
         endif
      endif

      if(ibm_debug_col.eq.1) then
         print*,"Col: Fn",rfn1,rfn2,rfn3,"Ft",rft        
      endif
!     $        ,"Ft",rft
!     $        ,"Total", fcf1
!     $        ,"Tq",ftq1
!     $        ,"rtorq",rtorq_col(1),rtorq_col(2),rtorq_col(3)
!     $        ,'rn',rn_12x,rn_12y,rn_12z
!     $        ,'rtmag',rtmag
!     $        ,"1 rnmag",c_fric,rnmag
!     $        ,"2 eta_t",eta_t,rv12t_mag
!     $        ,"Delta",rdelta12, rDelta_c
!     $        ,"spring damping",rksp_max,rv12_mage
!     $        ,"eta_t", eta_t,eta,rmult,mcfac,ksp
!     $        ,"rv12t",rv12t,rv12t_mag


 1513 continue

      return
      end

c-----------------------------------------------------------------------
      subroutine rotate3D(xo,xref,theta)
c      3D rotation of xo coordinate theta with reference to xref
c      ! Input xo, xref, theta
c      ! Output: xo 
      real xo(3), theta(3), xref(3)

c     local
      real rotmx(3,3), rotmy(3,3), rotmz(3,3) !3x3 matrix
      real rx1(3) ! relative coordinates

      ! 3 angles
      rtx = theta(1)
      rty = theta(2)
      rtz = theta(3)

      ! define rotation matrix
      rotmx(1,1) = 1.0
      rotmx(1,2) = 0.0
      rotmx(1,3) = 0.0
      rotmx(2,1) = 0.0
      rotmx(2,2) = cos(rtx)
      rotmx(2,3) = - sin(rtx)
      rotmx(3,1) = 0.0
      rotmx(3,2) = sin(rtx)
      rotmx(3,3) = cos(rtx)
      
      rotmy(1,1) = cos(rty)
      rotmy(1,2) = 0.0
      rotmy(1,3) = sin(rty)
      rotmy(2,1) = 0.0
      rotmy(2,2) = 1.0
      rotmy(2,3) = 0.0
      rotmy(3,1) = -sin(rty)
      rotmy(3,2) = 0.0
      rotmy(3,3) = cos(rty)

      rotmz(1,1) = cos(rtz)
      rotmz(1,2) = -sin(rtz)
      rotmz(1,3) = 0.0
      rotmz(2,1) = sin(rtz)
      rotmz(2,2) = cos(rtz)
      rotmz(2,3) = 0.0
      rotmz(3,1) = 0.0
      rotmz(3,2) = 0.0
      rotmz(3,3) = 1.0

      ! new vector
      rx1(1) = xo(1) - xref(1)
      rx1(2) = xo(2) - xref(2)
      rx1(3) = xo(3) - xref(3)
      
      call rot_multi(rotmx,rx1) ! rot x axis
      call rot_multi(rotmy,rx1) ! rot y axis
      call rot_multi(rotmz,rx1) ! rot z axis

      xo(1) = rx1(1) + xref(1)
      xo(2) = rx1(2) + xref(2)
      xo(3) = rx1(3) + xref(3)

      return
      end
      
c-----------------------------------------------------------------------
      subroutine rot_multi(rotm,rx)
c     matrix multiplication with rotate about x,y,z axis

      real rotm(3,3)
      real rx(3)
      real rn(3)

      rn(1) = rx(1)*rotm(1,1) + rx(2)*rotm(1,2) + rx(3)*rotm(1,3)
      rn(2) = rx(1)*rotm(2,1) + rx(2)*rotm(2,2) + rx(3)*rotm(2,3)
      rn(3) = rx(1)*rotm(3,1) + rx(2)*rotm(3,2) + rx(3)*rotm(3,3)

      rx(1) = rn(1)
      rx(2) = rn(2)
      rx(3) = rn(3)

      return
      end
c-----------------------------------------------------------------------

c----------------------------------------------------------------------
      subroutine output_debug_ibm_bins
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real rfpfluid(3),rfpfluidl(3),msum,msum_tot(3,2)

      ! show element map
      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
          ptw(ix,iy,iz,ie,5) = real(nid)
      enddo
      enddo
      enddo
      enddo
      
      itmp = 1
      call outpost2(ptw(1,1,1,1,1),         ! fhyd_x
     >              ptw(1,1,1,1,2),         ! fhyd_y
     >              ptw(1,1,1,1,3),         ! fhyd_z
     >              ptw(1,1,1,1,4),         ! phi_p (but not if lx1!=lx2
     >              ptw(1,1,1,1,4),         ! phi_p
     >              itmp          ,        
     >              'ptw')

      
c     eulerian integrations -----------------------------------------
c     fluid momentum 
      msum_tot(1,1) = glsc3(bm1,vtrans,vx,nx1*ny1*nz1*nelv)
      msum_tot(2,1) = glsc3(bm1,vtrans,vy,nx1*ny1*nz1*nelv)
      msum_tot(3,1) = glsc3(bm1,vtrans,vz,nx1*ny1*nz1*nelv)
c     particle volume fraction
      vf_part_e     = glsc2(bm1,ptw(1,1,1,1,4),nx1*ny1*nz1*nelt)
                                                 ! in z of mono-particle
                                                 ! Dp
c     particle forces on fluid
      rfpfluid(1)   = glsc2(bm1,ptw(1,1,1,1,1),nx1*ny1*nz1*nelt)
      rfpfluid(2)   = glsc2(bm1,ptw(1,1,1,1,2),nx1*ny1*nz1*nelt)
      rfpfluid(3)   = glsc2(bm1,ptw(1,1,1,1,3),nx1*ny1*nz1*nelt)

      if (.not.if3d) vf_part_e   = vf_part_e*dp(1)   ! Here:
      if (.not.if3d) rfpfluid(1) = rfpfluid(1)*dp(1) ! for 2d, assume
      if (.not.if3d) rfpfluid(2) = rfpfluid(2)*dp(1) ! z thicknes of 
      if (.not.if3d) rfpfluid(3) = rfpfluid(3)*dp(1) ! monodisperse Dp


c     lagrangian integrations ---------------------------------------
c     particle momentum
      do ieq=0,2
         msum = 0.0
         rsum = 0.0
         do i=1,n
           msum = msum + 
     >       rpart(jv0+ieq,i)*rpart(jrhop,i)*rpart(jvol,i)
           rsum = rsum + rpart(jf0+ieq,i)
        enddo
         msum_tot(ieq+1,2) = glsum(msum,1)
         rfpfluidl(1+ieq)  = glsum(rsum,1)
      enddo
c     particle volume fraction
      msum = 0.0
      do i=1,n
         msum = msum + rpart(jvol,i)
      enddo
      vf_part_l = glsum(msum,1)

      vf_rel_error = abs(vf_part_l - vf_part_e)/vf_part_l*100.0

c     print to files ------------------------------------------------
c     print properties to logfile
      if (nid.eq.0) write(6,500) "E Fluid Momentum :              ", 
     >                  istep, msum_tot(1,1),msum_tot(2,1),msum_tot(3,1)
      if (nid.eq.0) write(6,500) "E Particle forces:              ", 
     >                  istep, rfpfluid(1),rfpfluid(2),rfpfluid(3)
      if (nid.eq.0) write(6,500) "E Particle Volume:              ", 
     >                  istep, vf_part_e
      if (nid.eq.0) write(6,500) "L Particle Momentum :           ", 
     >                  istep, msum_tot(1,2),msum_tot(2,2),msum_tot(3,2)
      if (nid.eq.0) write(6,500) "L Particle forces:              ", 
     >                  istep, rfpfluidl(1),rfpfluidl(2),rfpfluidl(3)
      if (nid.eq.0) write(6,500) "L Particle Volume:              ", 
     >                  istep, vf_part_l
      if (nid.eq.0) write(6,500) "VF Relative Error %              ", 
     >                  istep, vf_rel_error

  500 FORMAT(A30,I20,3ES20.10)

      return
      end


c----------------------------------------------------------------------
      subroutine seach_queen_global

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer isearch_jpid1,isearch_jpid3,isearch_jqueen
      common/search_queen_gl/isearch_jpid1,isearch_jpid3,isearch_jqueen

      do i_qt = 1,n
         if(ipart(jrole,i_qt) .ne. 1) cycle ! real Queen                                                                                                                                                                                                                                                              
         i_qt_pid1  = ipart(jpid1,i_qt)
         i_qt_pid3  = ipart(jpid3,i_qt)
         i_qt_queen = ipart(jqueen,i_qt)
         if ( i_qt_pid1 .eq. isearch_jpid1 .and.
     $        i_qt_pid3 .eq. isearch_jpid3 .and.
     $        i_qt_queen.eq. isearch_jqueen ) then

            rxval = rpart(jx,i_qt)
            ryval = rpart(jy,i_qt)
            rzval = 0.
            if(if3d) rzval = rpart(jz,i_qt)

            iip    = floor((rxval-xdrange(1,1))/rdxgp_ibm)
            jjp    = floor((ryval-xdrange(1,2))/rdygp_ibm)
            kkp    = floor((rzval-xdrange(1,3))/rdzgp_ibm)
            print*,"Found",isearch_jpid1, isearch_jqueen,"Queen in",nid
     $           ,"Bin in", iip,jjp,kkp

         endif

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine output_queens

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer              stage,nstage
      common /tstepstage/ stage,nstage

      integer pth
      character*20 filename
      write(filename,'(A,I5.5,A4)') 'QP', nid,".dat"
      pth = 2666+nid
      open(unit=pth,file=filename,form="formatted",access='stream')

      do i =1,n
         if(ipart(jrole,i) .ne. 1) cycle ! real Queen
         rxval = rpart(jx,i)
         ryval = rpart(jy,i)
         rzval = rpart(jz,i)

         iip    = floor((rxval-xdrange(1,1))/rdxgp_ibm)
         jjp    = floor((ryval-xdrange(1,2))/rdygp_ibm)
         kkp    = floor((rzval-xdrange(1,3))/rdzgp_ibm)

         if(ibm_rotation.eq.1) then
            write(pth,8347) istep,stage
     $           ,ipart(jpid1,i),ipart(jpid2,i),ipart(jpid3,i)
     $           ,ipart(jqueen,i),iip,jjp,kkp
     $           ,(rpart(jx+j,i), j=0,2)
     $           ,(rpart(jv0+j,i),j=0,2)
     $           ,(rpart(jf0+j,i),j=0,2)
     >           ,(rpart(jangvel0+j,i),j=0,2)
     >           ,(rpart(jangle0+j,i), j=0,2)
     >           ,(rpart(jtorque0+j,i),j=0,2)

         else
            write(pth,8347) istep,stage
     $           ,ipart(jpid1,i),ipart(jpid2,i),ipart(jpid3,i)
     $           ,ipart(jqueen,i),iip,jjp,kkp
     $           ,(rpart(jx+j,i), j=0,2)
     $           ,(rpart(jv0+j,i),j=0,2)
     $           ,(rpart(jf0+j,i),j=0,2)
         endif

      enddo
 8347 FORMAT(9I6,20f12.4)

      return
      end
c----------------------------------------------------------------------
      subroutine output_workers

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'
      integer              stage,nstage
      common /tstepstage/ stage,nstage

      integer pth
      character*20 filename
      write(filename,'(A6,I5.5,A4)') 'WP', nid,".dat"
      pth = 7666+nid
      open(unit=pth,file=filename,form="formatted",access='stream')

      do i =1,n
         if(ipart(jrole,i) .ne. 2) cycle ! real Queen
         rxval = rpart(jx,i)
         ryval = rpart(jy,i)
         rzval = rpart(jz,i)

         iip    = floor((rxval-xdrange(1,1))/rdxgp_ibm)
         jjp    = floor((ryval-xdrange(1,2))/rdygp_ibm)
         kkp    = floor((rzval-xdrange(1,3))/rdzgp_ibm)

         write(pth,8378) istep,stage
     $        ,ipart(jpid1,i),ipart(jpid2,i),ipart(jpid3,i)
     $        ,ipart(jqueen,i),iip,jjp,kkp
     $        ,(rpart(jx+j,i),j=0,2),(rpart(jv0+j,i),j=0,2)
      enddo
 8378 FORMAT(9I6,10f12.6)

      return
      end
c----------------------------------------------------------------------
      subroutine  output_bin_structure

      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      integer pth
      character*20 filename
      write(filename,'(A6,I5.5,A4)') 'Q-bin-', nid,".dat"
      pth = 1666+nid
      open(unit=pth,file=filename,form="formatted",access='stream')
      do i =1,nlist_ibm
         write(pth,8397) nid,(ngp_valsp_ibm(j,i),j=1,6)
      enddo
 8397 FORMAT(10I6)
      close(pth)

      return
      end
c----------------------------------------------------------------------

      subroutine periodic_collide_check
c     For particle collision in periodic boundary conditions
c     update rx2 (icx2, icy2, icz2) based on both particle's bins
      ! (icx1, icy1, icz1) 
      ! (icx2, icy2, icz2)
      include 'SIZE'
      include 'TOTAL'
      include 'LPM'

      real rx2(3)
      integer icx1,icy1,icz1,icx2,icy2,icz2
      common /col_part1/ icx1,icy1,icz1
      common /col_part2/ rx2,icx2,icy2,icz2

      xdlen = xdrange(2,1) - xdrange(1,1)
      ydlen = xdrange(2,2) - xdrange(1,2)
      zdlen = xdrange(2,3) - xdrange(1,3)

      ! x-dir
      if(xdlen .gt. 0 ) then 
         if(icx1 .le. 0 .and. icx2 .ge. ndxgp_ibm-1) then
            rx2(1) = rx2(1) - xdlen
            icx2   = icx2   - ndxgp_ibm 
         endif
         if(icx1 .ge. ndxgp_ibm-1 .and. icx2 .le. 0) then
            rx2(1) = rx2(1) + xdlen
            icx2   = icx2   + ndxgp_ibm 
         endif
      endif

      ! y-dir
      if(ydlen .gt. 0 ) then 
         if(icy1 .le. 0 .and. icy2 .ge. ndygp_ibm-1) then
            rx2(2) = rx2(2) - ydlen
            icy2   = icy2   - ndygp_ibm 
         endif
         if(icy1 .ge. ndygp_ibm-1 .and. icy2 .le. 0) then
            rx2(2) = rx2(2) + ydlen
            icy2   = icy2   + ndygp_ibm 
         endif
      endif

      ! z-dir
      if(zdlen .gt. 0 ) then 
         if(icz1 .le. 0   .and. icz2 .ge. ndzgp_ibm-1) then
            rx2(3) = rx2(3) - zdlen
            icz2   = icz2   - ndzgp_ibm 
         endif
         if(icz1 .ge. ndygp_ibm-1 .and. icz2 .le. 0) then
            rx2(3) = rx2(3) + zdlen
            icz2   = icz2   + ndzgp_ibm
         endif
      endif

         
      return
      end
      
