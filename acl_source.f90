
module actuator_line_source

    use decomp_2d, only: mytype
    use actuator_line_model_utils
    use actuator_line_model

    implicit none
    real(mytype),save :: constant_epsilon, meshFactor, thicknessFactor,chordFactor
    real(mytype),save, allocatable :: Sx(:),Sy(:),Sz(:),Sc(:),Se(:),Sh(:),Su(:),Sv(:),Sw(:),SFX(:),SFY(:),SFZ(:), sum_kernel(:)
    real(mytype),save, allocatable :: Su_part(:),Sv_part(:),Sw_part(:), sum_kernel_part(:)
    real(mytype),save, allocatable :: Snx(:),Sny(:),Snz(:),Stx(:),Sty(:),Stz(:),Ssx(:),Ssy(:),Ssz(:),Ssegm(:)
    real(mytype),save, allocatable :: A(:,:)
    logical, allocatable :: inside_the_domain(:)
    integer :: NSource
    logical, save :: rbf_interpolation=.false.
    logical, save :: pointwise_interpolation=.false.
    logical, save :: anisotropic_projection=.false.
    logical, save :: has_mesh_based_epsilon=.false.
    logical, save :: has_constant_epsilon=.false.
    public get_locations, get_forces, set_vel, initialize_actuator_source

contains

    subroutine initialize_actuator_source

    use param, only: nx, ny, nz

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            ! Blades
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                end do
            end do
            ! Tower
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                end do
            endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
            end do
        end do
    endif
    NSource=counter
    allocate(Sx(NSource),Sy(NSource),Sz(NSource),Sc(Nsource),Su(NSource),Sv(NSource),Sw(NSource),Se(NSource),Sh(NSource),Sfx(NSource),Sfy(NSource),Sfz(NSource), sum_kernel(Nsource))
    allocate(Su_part(NSource),Sv_part(NSource),Sw_part(NSource), sum_kernel_part(NSource))
    allocate(Snx(NSource),Sny(NSource),Snz(NSource),Stx(Nsource),Sty(NSource),Stz(NSource),Ssx(NSource),Ssy(NSource),Ssz(NSource),Ssegm(NSource))
    allocate(A(NSource,NSource))
    allocate(inside_the_domain(NSource))

    end subroutine initialize_actuator_source

    subroutine get_locations

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%Blade(iblade)%PEX(ielem)
                Sy(counter)=Turbine(itur)%Blade(iblade)%PEY(ielem)
                Sz(counter)=Turbine(itur)%Blade(iblade)%PEZ(ielem)
                Sc(counter)=Turbine(itur)%Blade(iblade)%EC(ielem)
                Snx(counter)=Turbine(itur)%Blade(iblade)%nEx(ielem)
                Sny(counter)=Turbine(itur)%Blade(iblade)%nEy(ielem)
                Snz(counter)=Turbine(itur)%Blade(iblade)%nEz(ielem)
                Stx(counter)=Turbine(itur)%Blade(iblade)%tEx(ielem)
                Sty(counter)=Turbine(itur)%Blade(iblade)%tEy(ielem)
                Stz(counter)=Turbine(itur)%Blade(iblade)%tEz(ielem)
                Ssx(counter)=Turbine(itur)%Blade(iblade)%sEx(ielem)
                Ssy(counter)=Turbine(itur)%Blade(iblade)%sEy(ielem)
                Ssz(counter)=Turbine(itur)%Blade(iblade)%sEz(ielem)
                Ssegm(counter)=Turbine(itur)%Blade(iblade)%EDS(ielem)

                end do
            end do
                !Tower
                if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%Tower%PEX(ielem)
                Sy(counter)=Turbine(itur)%Tower%PEY(ielem)
                Sz(counter)=Turbine(itur)%Tower%PEZ(ielem)
                Sc(counter)=Turbine(itur)%Tower%EC(ielem)
                Snx(counter)=Turbine(itur)%Tower%nEx(ielem)
                Sny(counter)=Turbine(itur)%Tower%nEy(ielem)
                Snz(counter)=Turbine(itur)%Tower%nEz(ielem)
                Stx(counter)=Turbine(itur)%Tower%tEx(ielem)
                Sty(counter)=Turbine(itur)%Tower%tEy(ielem)
                Stz(counter)=Turbine(itur)%Tower%tEz(ielem)
                Ssx(counter)=Turbine(itur)%Tower%sEx(ielem)
                Ssy(counter)=Turbine(itur)%Tower%sEy(ielem)
                Ssz(counter)=Turbine(itur)%Tower%sEz(ielem)
                Ssegm(counter)=Turbine(itur)%Tower%EDS(ielem)
                end do
                endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sx(counter)=actuatorline(ial)%PEX(ielem)
                Sy(counter)=actuatorline(ial)%PEY(ielem)
                Sz(counter)=actuatorline(ial)%PEZ(ielem)
                Sc(counter)=actuatorline(ial)%EC(ielem)
                Snx(counter)=actuatorline(ial)%nEx(ielem)
                Sny(counter)=actuatorline(ial)%nEy(ielem)
                Snz(counter)=actuatorline(ial)%nEz(ielem)
                Stx(counter)=actuatorline(ial)%tEx(ielem)
                Sty(counter)=actuatorline(ial)%tEy(ielem)
                Stz(counter)=actuatorline(ial)%tEz(ielem)
                Ssx(counter)=actuatorline(ial)%sEx(ielem)
                Ssy(counter)=actuatorline(ial)%sEy(ielem)
                Ssz(counter)=actuatorline(ial)%sEz(ielem)
                Ssegm(counter)=actuatorline(ial)%EDS(ielem)
            end do
        end do
    endif

    end subroutine get_locations

    subroutine set_vel

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            ! Blades
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Turbine(itur)%Blade(iblade)%EVX(ielem)=Su(counter)
                Turbine(itur)%Blade(iblade)%EVY(ielem)=Sv(counter)
                Turbine(itur)%Blade(iblade)%EVZ(ielem)=Sw(counter)
                Turbine(itur)%Blade(iblade)%Eepsilon(ielem)=Se(counter)
                end do
            end do
            ! Tower
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Turbine(itur)%Tower%EVX(ielem)=Su(counter)
                Turbine(itur)%Tower%EVY(ielem)=Sv(counter)
                Turbine(itur)%Tower%EVZ(ielem)=Sw(counter)
                Turbine(itur)%Tower%Eepsilon(ielem)=Se(counter)
                end do
            endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                actuatorline(ial)%EVX(ielem)=Su(counter)
                actuatorline(ial)%EVY(ielem)=Sv(counter)
                actuatorline(ial)%EVZ(ielem)=Sw(counter)
                actuatorline(ial)%Eepsilon(ielem)=Se(counter)
            end do
        end do
    endif


    end subroutine set_vel

    subroutine get_forces

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            !Blade
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%Blade(iblade)%EFX(ielem)
                Sfy(counter)=Turbine(itur)%Blade(iblade)%EFY(ielem)
                Sfz(counter)=Turbine(itur)%Blade(iblade)%EFZ(ielem)
                end do
            end do

            !Tower
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%Tower%EFX(ielem)
                Sfy(counter)=Turbine(itur)%Tower%EFY(ielem)
                Sfz(counter)=Turbine(itur)%Tower%EFZ(ielem)
                end do
            endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sfx(counter)=actuatorline(ial)%EFX(ielem)
                Sfy(counter)=actuatorline(ial)%EFY(ielem)
                Sfz(counter)=actuatorline(ial)%EFZ(ielem)
            end do
        end do
    endif

    end subroutine get_forces

    subroutine Compute_Momentum_Source_Term_pointwise_new

        use decomp_2d, only: mytype, nproc, xstart, xend, xsize
        use MPI
        use param, only: dx,dy,dz,eps_factor,xnu,yp,istret,xlx,yly,zlz, istret, nx, ny, nz
        use var, only: ux1, uy1, uz1, FTx, FTy, FTz

        implicit none
        real(mytype) :: xmesh, ymesh,zmesh
        real(mytype) :: dist, epsilon, Kernel
        integer :: i,j,k, isource, extended_cells
        integer :: i_source, j_source, k_source, ierr
        integer :: first_i_sample, last_i_sample, first_j_sample,last_j_sample, first_k_sample,last_k_sample

        ! Checks
        real(mytype) :: sum_fx_grid, sum_fy_grid, sum_fz_grid, sum_fx_al, sum_fy_al, sum_fz_al
        real(mytype) :: sum_fx_grid_part, sum_fy_grid_part, sum_fz_grid_part

        if(istret.ne.0) then
            write(*,*) "This part is not ready for refinement, change branch"
            stop
        endif

        ! First we need to compute the locations of the AL points
        call get_locations

        ! Zero the velocities and other variables
        Su(:)=0.0
        Sv(:)=0.0
        Sw(:)=0.0
        Su_part(:)=0.0
        Sv_part(:)=0.0
        Sw_part(:)=0.0
        sum_kernel(:)=0.0
        sum_kernel_part(:)=0.0

        ! Check if the points lie outside the fluid domain
        do isource=1,Nsource
            if((Sx(isource)>xlx).or.(Sx(isource)<0).or.(Sy(isource)>yly).or.(Sy(isource)<0).or.(Sz(isource)>zlz).or.(Sz(isource)<0)) then
                print *, 'Point outside the fluid domain'
                stop
            endif
        enddo

        ! Compute the variables of the Gaussian function
        epsilon = eps_factor*(dx*dy*dz)**(1/3)
        extended_cells = int(5*epsilon/min(dx,dy,dz) + 1)

        ! VELOCITY SAMPLING
        do isource=1,NSource

            i_source = int(Sx(isource)/dx)
            j_source = int(Sy(isource)/dy)
            k_source = int(Sz(isource)/dz)

            ! Check that the cell that I am going to use belongs to the pencil I am in
            if((i_source.ge.xstart(1)).and.(i_source.lt.xend(1))) then
                if((j_source.ge.xstart(2)).and.(j_source.lt.xend(2))) then
                    if((k_source.ge.xstart(3)).and.(k_source.lt.xend(3))) then
                        ! Compute the position of the nodes to be sampled
                        xmesh = (i_source - 1)*dx
                        ymesh = (j_source - 1)*dy
                        zmesh = (k_source - 1)*dz

                        ! Probably I do not need this for a regular mesh
                        ! I dont know how dangerous is to use parallelisation in a place where its not needed
                        Su_part(isource)= trilinear_interpolation(xmesh      ,ymesh      ,zmesh, &
                                                             xmesh+dx   ,ymesh+dy   ,zmesh+dz, &
                                                             Sx(isource),Sy(isource),Sz(isource  ), &
                                                             ux1(i_source  ,j_source  ,k_source  ), &
                                                             ux1(i_source+1,j_source  ,k_source  ), &
                                                             ux1(i_source  ,j_source  ,k_source+1), &
                                                             ux1(i_source+1,j_source  ,k_source+1), &
                                                             ux1(i_source  ,j_source+1,k_source  ), &
                                                             ux1(i_source+1,j_source+1,k_source  ), &
                                                             ux1(i_source  ,j_source+1,k_source+1), &
                                                             ux1(i_source+1,j_source+1,k_source+1))

                        Sv_part(isource)= trilinear_interpolation(xmesh      ,ymesh      ,zmesh, &
                                                             xmesh+dx   ,ymesh+dy   ,zmesh+dz, &
                                                             Sx(isource),Sy(isource),Sz(isource  ), &
                                                             uy1(i_source  ,j_source  ,k_source  ), &
                                                             uy1(i_source+1,j_source  ,k_source  ), &
                                                             uy1(i_source  ,j_source  ,k_source+1), &
                                                             uy1(i_source+1,j_source  ,k_source+1), &
                                                             uy1(i_source  ,j_source+1,k_source  ), &
                                                             uy1(i_source+1,j_source+1,k_source  ), &
                                                             uy1(i_source  ,j_source+1,k_source+1), &
                                                             uy1(i_source+1,j_source+1,k_source+1))

                        Sw_part(isource)= trilinear_interpolation(xmesh      ,ymesh      ,zmesh, &
                                                             xmesh+dx   ,ymesh+dy   ,zmesh+dz, &
                                                             Sx(isource),Sy(isource),Sz(isource  ), &
                                                             uz1(i_source  ,j_source  ,k_source  ), &
                                                             uz1(i_source+1,j_source  ,k_source  ), &
                                                             uz1(i_source  ,j_source  ,k_source+1), &
                                                             uz1(i_source+1,j_source  ,k_source+1), &
                                                             uz1(i_source  ,j_source+1,k_source  ), &
                                                             uz1(i_source+1,j_source+1,k_source  ), &
                                                             uz1(i_source  ,j_source+1,k_source+1), &
                                                             uz1(i_source+1,j_source+1,k_source+1))
                    endif
                endif
            endif
        enddo

        ! Retrieve information from all the processes
        call MPI_ALLREDUCE(Su_part,Su,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(Sv_part,Sv,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(Sw_part,Sw,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        if(nrank==0) then
            write(*,*) "Su", Su
            write(*,*) "Sv", Sv
            write(*,*) "Sw", Sw
        endif

        ! COMPUTATION OF THE INTEGRAL OF THE INTEGRATION KERNEL
        do isource=1,NSource
            ! Indices of the node just before the sampling point
            ! Indices used in the whole domain and the pencils
            i_source = int(Sx(isource)/dx)
            j_source = int(Sy(isource)/dy)
            k_source = int(Sz(isource)/dz)

            ! Indices of the vertices to sample.
            first_i_sample = max(i_source - (extended_cells - 1), xstart(1))
            last_i_sample = min(i_source + extended_cells, xend(1))

            first_j_sample =  max(j_source - (extended_cells - 1), xstart(2))
            last_j_sample = min(j_source + extended_cells, xend(2))

            first_k_sample = max(k_source - (extended_cells - 1), xstart(3))
            last_k_sample = min(k_source + extended_cells, xend(3))

            ! Loop through the points
            do k=first_k_sample, last_k_sample
                do j=first_j_sample, last_j_sample
                    do i=first_i_sample, last_i_sample
                        ! Compute the position of the nodes to be sampled
                        xmesh = (i - 1)*dx
                        ymesh = (j - 1)*dy
                        zmesh = (k - 1)*dz

                        dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
                        Kernel= 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)
                        sum_Kernel_part(isource) = sum_Kernel_part(isource) + Kernel
                    enddo
                enddo
            enddo
        enddo ! loop through the sources

        ! Sum the kernel
        call MPI_ALLREDUCE(sum_kernel_part,sum_kernel,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        ! PROJECT THE FORCES
        FTx(:,:,:)=0.0
        FTy(:,:,:)=0.0
        FTz(:,:,:)=0.0
        ! SFx(:)=0.0
        ! SFy(:)=0.0
        ! SFz(:)=0.0
        Visc=xnu
        !## Send the velocities to the
        call set_vel
        !## Compute forces
        call actuator_line_model_compute_forces
        !## Get Forces
        call get_forces

        ! if(nrank==0) then
        !     write(*,*) "SFx", SFx
        !     write(*,*) "SFy", SFy
        !     write(*,*) "SFz", SFx
        ! endif

        ! Loop through all the nodes of the domain
        ! Each process through theirs
        do k=1,xsize(3)
            do j=1,xsize(2)
                do i=1,xsize(1)
                    xmesh = (i + xstart(1) - 1)*dx
                    ymesh = (j + xstart(2) - 1)*dy
                    zmesh = (k + xstart(3) - 1)*dz

                    do isource=1,Nsource
                        ! Distance from the node to the AL point
                        dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
                        ! Gaussian Kernel
                        ! I see this dangerous anyway
                        if(dist.lt.(extended_cells*min(dx,dy,dz))) then
                            Kernel= 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)
                            FTx(i,j,k)=FTx(i,j,k)-SFx(isource)*Kernel/sum_kernel(isource)
                            FTy(i,j,k)=FTy(i,j,k)-SFy(isource)*Kernel/sum_kernel(isource)
                            FTz(i,j,k)=FTz(i,j,k)-SFz(isource)*Kernel/sum_kernel(isource)
                            ! write(*,*) nrank, "ijk", i, j, k, "FTx", FTx(i,j,k)
                            ! write(*,*) nrank, "ijk", i, j, k, "FTy", FTy(i,j,k)
                            ! write(*,*) nrank, "ijk", i, j, k, "FTz", FTz(i,j,k)
                        endif
                    enddo
                enddo
            enddo
        enddo

        ! checks
        sum_fx_al = 0.0
        sum_fy_al = 0.0
        sum_fz_al = 0.0
        sum_fx_grid_part = 0.0
        sum_fy_grid_part = 0.0
        sum_fz_grid_part = 0.0

        do isource=1,Nsource
            sum_fx_al = sum_fx_al + SFx(isource)
            sum_fy_al = sum_fy_al + SFy(isource)
            sum_fz_al = sum_fz_al + SFz(isource)
        enddo

        do k=1,xsize(3)
            do j=1,xsize(2)
                do i=1,xsize(1)
                    sum_fx_grid_part = sum_fx_grid_part + FTx(i,j,k)
                    sum_fy_grid_part = sum_fy_grid_part + FTy(i,j,k)
                    sum_fz_grid_part = sum_fz_grid_part + FTz(i,j,k)
                enddo
            enddo
        enddo

        call MPI_ALLREDUCE(sum_fx_grid_part,sum_fx_grid,1,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(sum_fy_grid_part,sum_fy_grid,1,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(sum_fz_grid_part,sum_fz_grid,1,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        ! write(*,*) "Check fx", sum_fx_al, sum_fx_grid
        ! write(*,*) "Check fy", sum_fy_al, sum_fy_grid
        ! write(*,*) "Check fz", sum_fz_al, sum_fz_grid

    end subroutine Compute_Momentum_Source_Term_pointwise_new

    subroutine Compute_Momentum_Source_Term_pointwise

        use decomp_2d, only: mytype, nproc, xstart, xend, xsize, update_halo
        use MPI
        use param, only: dx,dy,dz,eps_factor,xnu,yp,istret,xlx,yly,zlz, rho_air
        use var, only: ux1, uy1, uz1, FTx, FTy, FTz

        implicit none
        real(mytype), allocatable, dimension(:,:,:) :: ux1_halo, uy1_halo, uz1_halo
        real(mytype) :: xmesh, ymesh,zmesh
        real(mytype) :: dist, epsilon, Kernel
        real(mytype) :: min_dist, ymax,ymin,zmin,zmax
        real(mytype) :: x0,y0,z0,x1,y1,z1,x,y,z,u000,u100,u001,u101,u010,u110,u011,u111
        real(mytype) :: t1,t2, alm_proj_time, constant
        integer :: min_i,min_j,min_k
        integer :: i_lower, j_lower, k_lower, i_upper, j_upper, k_upper
        integer :: i,j,k, isource, ierr

        real(mytype), dimension(Nsource) :: sum_kernel, sum_kernel_part

        constant = dx*dy*dz*rho_air
        ! First we need to compute the locations
        call get_locations

        ! Zero the velocities
        Su(:)=0.0
        Sv(:)=0.0
        Sw(:)=0.0
        ! This is not optimum but works
        ! Define the domain

        if (istret.eq.0) then
        ymin=(xstart(2)-1)*dy-dy/2. ! Add -dy/2.0 overlap
        ymax=(xend(2)-1)*dy+dy/2.   ! Add +dy/2.0 overlap
        else
        ymin=yp(xstart(2))
        ymax=yp(xend(2))
        endif

        zmin=(xstart(3)-1)*dz-dz/2. ! Add a -dz/2.0 overlap
        zmax=(xend(3)-1)*dz+dz/2.   ! Add a +dz/2.0 overlap

        ! Check if the points lie outside the fluid domain
        do isource=1,Nsource
        if((Sx(isource)>xlx).or.(Sx(isource)<0).or.(Sy(isource)>yly).or.(Sy(isource)<0).or.(Sz(isource)>zlz).or.(Sz(isource)<0)) then
            print *, 'Point outside the fluid domain'
            stop
        endif
        enddo
        !write(*,*) 'Rank=', nrank, 'X index Limits=', xstart(1), xend(1), 'X lims=', (xstart(1)-1)*dx, (xend(1)-1)*dx
        !write(*,*) 'Rank=', nrank, 'Y index Limits=', xstart(2), xend(2), 'Y lims=', ymin, ymax
        !write(*,*) 'Rank=', nrank, 'Z index Limits=', xstart(3), xend(3), 'Z lims=', zmin, zmax
        call update_halo(ux1,ux1_halo,1,opt_global=.true.)
        call update_halo(uy1,uy1_halo,1,opt_global=.true.)
        call update_halo(uz1,uz1_halo,1,opt_global=.true.)
        !print *,  nrank, shape(ux1), shape(ux1_halo)
        !do j=xstart(2),xend(2)
        !print *, nrank, ux1(xstart(1),j,xstart(3)), ux1_halo(xstart(1),j,xstart(3))
        !enddo

        do isource=1,NSource

        min_dist=1e6
        if((Sy(isource)>=ymin).and.(Sy(isource)<ymax).and.(Sz(isource)>=zmin).and.(Sz(isource)<zmax)) then
            !write(*,*) 'nrank= ',nrank, 'owns this node'
            do k=xstart(3),xend(3)
            zmesh=(k-1)*dz
            do j=xstart(2),xend(2)
            if (istret.eq.0) ymesh=(j-1)*dy
            if (istret.ne.0) ymesh=yp(j)
            do i=xstart(1),xend(1)
            xmesh=(i-1)*dx
            dist = sqrt((Sx(isource)-xmesh)**2.+(Sy(isource)-ymesh)**2.+(Sz(isource)-zmesh)**2.)

            if (dist<min_dist) then
                min_dist=dist
                min_i=i
                min_j=j
                min_k=k
            endif

            enddo
            enddo
            enddo

            if(Sy(isource)>ymax.or.Sy(isource)<ymin) then
            write(*,*) 'In processor ', nrank
            write(*,*) 'Sy =', Sy(isource),'is not within the', ymin, ymax, 'limits'
            stop
            endif
            if(Sz(isource)>zmax.or.Sz(isource)<zmin) then
            write(*,*) 'In processor ', nrank
            write(*,*) 'Sz =', Sz(isource),'is not within the', zmin, zmax, 'limits'
            stop
            endif

            if(Sx(isource)>(min_i-1)*dx) then
                i_lower=min_i
                i_upper=min_i+1
            else if(Sx(isource)<(min_i-1)*dx) then
                i_lower=min_i-1
                i_upper=min_i
            else if(Sx(isource)==(min_i-1)*dx) then
                i_lower=min_i
                i_upper=min_i
            endif

            if (istret.eq.0) then
            if(Sy(isource)>(min_j-1)*dy.and.Sy(isource)<(xend(2)-1)*dy) then
                j_lower=min_j
                j_upper=min_j+1
            else if(Sy(isource)>(min_j-1)*dy.and.Sy(isource)>(xend(2)-1)*dy) then
                j_lower=min_j
                j_upper=min_j+1 ! THis is in the Halo domain
            else if(Sy(isource)<(min_j-1)*dy.and.Sy(isource)>(xstart(2)-1)*dy) then
                j_lower=min_j-1
                j_upper=min_j
            else if(Sy(isource)<(min_j-1)*dy.and.Sy(isource)<(xstart(2)-1)*dy) then
                j_lower=min_j-1 ! THis is in the halo domain
                j_upper=min_j
            else if (Sy(isource)==(min_j-1)*dy) then
                j_lower=min_j
                j_upper=min_j
            endif
            else

            if(Sy(isource)>yp(min_j)) then
                j_lower=min_j
                j_upper=min_j+1
            else if(Sy(isource)<yp(min_j)) then
                j_lower=min_j-1
                j_upper=min_j
            else if (Sy(isource)==yp(min_j)) then
                j_lower=min_j
                j_upper=min_j
            endif
            endif

            if(Sz(isource)>(min_k-1)*dz.and.Sz(isource)<(xend(3)-1)*dz) then
                k_lower=min_k
                k_upper=min_k+1
            elseif(Sz(isource)>(min_k-1)*dz.and.Sz(isource)>(xend(3)-1)*dz) then
                k_lower=min_k
                k_upper=min_k+1 ! This in the halo doamin
            else if(Sz(isource)<(min_k-1)*dz.and.Sz(isource)>(xstart(3)-1)*dz) then
                k_lower=min_k-1
                k_upper=min_k
            else if(Sz(isource)<(min_k-1)*dz.and.Sz(isource)<(xstart(3)-1)*dz) then
                k_lower=min_k-1 ! This is in the halo domain
                k_upper=min_k
            else if (Sz(isource)==(min_k-1)*dz) then
                k_lower=min_k
                k_upper=min_k
            endif

            ! Prepare for interpolation
            x0=(i_lower-1)*dx
            x1=(i_upper-1)*dx
            if (istret.eq.0) then
            y0=(j_lower-1)*dy
            y1=(j_upper-1)*dy
            else
            y0=yp(j_lower)
            y1=yp(j_upper)
            endif
            z0=(k_lower-1)*dz
            z1=(k_upper-1)*dz

            x=Sx(isource)
            y=Sy(isource)
            z=Sz(isource)

            !if(x>x1.or.x<x0.or.y>y1.or.y<y0.or.z>z1.or.z<z0) then
            !    write(*,*) 'x0, x1, x', x0, x1, x
            !    write(*,*) 'y0, y1, y', y0, y1, y
            !    write(*,*) 'z0, z1, z', z0, z1, z
            !    write(*,*) 'Problem with the trilinear interpolation';
            !    stop
            !endif
            ! Apply interpolation kernels from 8 neighboring nodes

            Su_part(isource)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  ux1_halo(i_lower,j_lower,k_lower), &
                                                  ux1_halo(i_upper,j_lower,k_lower), &
                                                  ux1_halo(i_lower,j_lower,k_upper), &
                                                  ux1_halo(i_upper,j_lower,k_upper), &
                                                  ux1_halo(i_lower,j_upper,k_lower), &
                                                  ux1_halo(i_upper,j_upper,k_lower), &
                                                  ux1_halo(i_lower,j_upper,k_upper), &
                                                  ux1_halo(i_upper,j_upper,k_upper))

             Sv_part(isource)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  uy1_halo(i_lower,j_lower,k_lower), &
                                                  uy1_halo(i_upper,j_lower,k_lower), &
                                                  uy1_halo(i_lower,j_lower,k_upper), &
                                                  uy1_halo(i_upper,j_lower,k_upper), &
                                                  uy1_halo(i_lower,j_upper,k_lower), &
                                                  uy1_halo(i_upper,j_upper,k_lower), &
                                                  uy1_halo(i_lower,j_upper,k_upper), &
                                                  uy1_halo(i_upper,j_upper,k_upper))

            Sw_part(isource)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  uz1_halo(i_lower,j_lower,k_lower), &
                                                  uz1_halo(i_upper,j_lower,k_lower), &
                                                  uz1_halo(i_lower,j_lower,k_upper), &
                                                  uz1_halo(i_upper,j_lower,k_upper), &
                                                  uz1_halo(i_lower,j_upper,k_lower), &
                                                  uz1_halo(i_upper,j_upper,k_lower), &
                                                  uz1_halo(i_lower,j_upper,k_upper), &
                                                  uz1_halo(i_upper,j_upper,k_upper))

        else
            Su_part(isource)=0.0
            Sv_part(isource)=0.0
            Sw_part(isource)=0.0
            !write(*,*) 'Warning: I do not own this node'
        endif
        enddo
        !$OMP END PARALLEL DO

        call MPI_ALLREDUCE(Su_part,Su,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(Sv_part,Sv,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(Sw_part,Sw,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        ! Zero the Source term at each time step
        FTx(:,:,:)=0.0
        FTy(:,:,:)=0.0
        FTz(:,:,:)=0.0
        Visc=xnu
        !## Send the velocities to the
        call set_vel
        !## Compute forces
        call actuator_line_model_compute_forces
        !## Get Forces
        call get_forces

        if(nrank==0) then
            write(*,*) 'Projecting the AL Momentum Source term ... '
        endif
        t1 = MPI_WTIME()

            !## Compute the integral of the kernel
            sum_kernel_part(:) = 0.0
            sum_kernel(:) = 0.0
            do k=1,xsize(3)
            zmesh=(k+xstart(3)-1-1)*dz
            do j=1,xsize(2)
            if (istret.eq.0) ymesh=(xstart(2)+j-1-1)*dy
            if (istret.ne.0) ymesh=yp(xstart(2)+j-1)
            do i=1,xsize(1)
            xmesh=(i-1)*dx

            do isource=1,NSource

            dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
            epsilon=eps_factor*(dx*dy*dz)**(1.0/3.0)
            if (dist<10.0*epsilon) then
                Kernel = 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)
                sum_Kernel_part(isource) = sum_Kernel_part(isource) + Kernel
            else
                Kernel=0.0
            endif
            enddo
            enddo
            enddo
            enddo
            call MPI_ALLREDUCE(sum_kernel_part,sum_kernel,Nsource,MPI_REAL8,MPI_SUM, &
                MPI_COMM_WORLD,ierr)

            !## Add the source term
            do k=1,xsize(3)
            zmesh=(k+xstart(3)-1-1)*dz
            do j=1,xsize(2)
            if (istret.eq.0) ymesh=(xstart(2)+j-1-1)*dy
            if (istret.ne.0) ymesh=yp(xstart(2)+j-1)
            do i=1,xsize(1)
            xmesh=(i-1)*dx

            do isource=1,NSource

            dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
            epsilon=eps_factor*(dx*dy*dz)**(1.0/3.0)
            if (dist<10.0*epsilon) then
                Kernel= 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)
            else
                Kernel=0.0
            endif
            ! First apply a constant lift to induce the
            FTx(i,j,k)=FTx(i,j,k)-SFx(isource)*Kernel/sum_kernel(isource)/constant
            FTy(i,j,k)=FTy(i,j,k)-SFy(isource)*Kernel/sum_kernel(isource)/constant
            FTz(i,j,k)=FTz(i,j,k)-SFz(isource)*Kernel/sum_kernel(isource)/constant

            enddo

            enddo
            enddo
            enddo

        alm_proj_time=MPI_WTIME()-t1
        call MPI_ALLREDUCE(alm_proj_time,t1,1,MPI_REAL8,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)

        if(nrank==0) then
            alm_proj_time=alm_proj_time/float(nproc)
            write(*,*) 'AL Momentum Source term projection completed in :', alm_proj_time ,'seconds'
        endif

    end subroutine Compute_Momentum_Source_Term_pointwise






    subroutine Compute_Momentum_Source_Term_integral

        use decomp_2d, only: mytype, nproc, xstart, xend, xsize
        use MPI
        use param, only: dx,dy,dz,eps_factor,xnu,yp,xlx,yly,zlz, istret, nx, ny, nz
        use var, only: ux1, uy1, uz1, FTx, FTy, FTz

        implicit none
        real(mytype) :: xmesh, ymesh,zmesh
        real(mytype) :: dist, epsilon, Kernel
        integer :: i,j,k, isource, extended_cells
        integer :: i_source, j_source, k_source, ierr
        integer :: first_i_sample, last_i_sample, first_j_sample,last_j_sample, first_k_sample,last_k_sample

        ! Checks
        real(mytype) :: sum_fx_grid, sum_fy_grid, sum_fz_grid, sum_fx_al, sum_fy_al, sum_fz_al
        real(mytype) :: sum_fx_grid_part, sum_fy_grid_part, sum_fz_grid_part

        if(istret.ne.0) then
            write(*,*) "This part is not ready for refinement, change branch"
            stop
        endif
        ! First we need to compute the locations of the AL points (Sx, Sy, Sz ...)
        call get_locations

        ! Zero the velocities
        Su(:)=0.0
        Sv(:)=0.0
        Sw(:)=0.0
        Su_part(:)=0.0
        Sv_part(:)=0.0
        Sw_part(:)=0.0
        sum_kernel(:)=0.0
        sum_kernel_part(:)=0.0

        ! Check if the points lie outside the fluid domain
        do isource=1,Nsource
          if((Sx(isource)>xlx).or.(Sx(isource)<0).or.(Sy(isource)>yly).or.(Sy(isource)<0).or.(Sz(isource)>zlz).or.(Sz(isource)<0)) then
            print *, 'Point outside the fluid domain'
            stop
          endif
        enddo

        ! Compute the epsilon and the number of cells that the pencils should be
        ! extended to accomodate a minimum of 10 cells around the sources
        ! I will extend all the directions the same number of cells (conservative)
        epsilon = eps_factor*(dx*dy*dz)**(1/3)
        extended_cells = int(5*epsilon/min(dx,dy,dz) + 1)

        ! SAMPLE VELOCITIES AND COMPUTE THE INTEGRATION OF THE KERNEL
        ! Loop through the sources
        do isource=1,NSource

            ! Indices of the node just before the sampling point
            ! Indices used in the whole domain and the pencils
            i_source = int(Sx(isource)/dx)
            j_source = int(Sy(isource)/dy)
            k_source = int(Sz(isource)/dz)

            ! Indices of the vertices to sample.
            first_i_sample = max(i_source - (extended_cells - 1), xstart(1))
            last_i_sample = min(i_source + extended_cells, xend(1))

            first_j_sample =  max(j_source - (extended_cells - 1), xstart(2))
            last_j_sample = min(j_source + extended_cells, xend(2))

            first_k_sample = max(k_source - (extended_cells - 1), xstart(3))
            last_k_sample = min(k_source + extended_cells, xend(3))

            ! Loop through the points
            do k=first_k_sample, last_k_sample
                do j=first_j_sample, last_j_sample
                    do i=first_i_sample, last_i_sample
                        ! Compute the position of the nodes to be sampled
                        xmesh = (i - 1)*dx
                        ymesh = (j - 1)*dy
                        zmesh = (k - 1)*dz

                        ! Distance from the node to the AL point
                        dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
                        ! Gaussian Kernel
                        Kernel= 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)
                        sum_Kernel_part(isource) = sum_Kernel_part(isource) + Kernel
                        ! Integration
                        Su_part(isource) = Su_part(isource) + Kernel*ux1(i, j, k)
                        Sv_part(isource) = Sv_part(isource) + Kernel*uy1(i, j, k)
                        Sw_part(isource) = Sw_part(isource) + Kernel*uz1(i, j, k)
                    enddo
                enddo
            enddo
        enddo ! loop through the sources

        ! Retrieve information from all the processes
        call MPI_ALLREDUCE(Su_part,Su,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(Sv_part,Sv,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(Sw_part,Sw,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        call MPI_ALLREDUCE(sum_kernel_part,sum_kernel,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        ! Apply the integration kernel
        do isource=1,NSource
            Su(isource) = Su(isource)/sum_kernel(isource)
            Sv(isource) = Sv(isource)/sum_kernel(isource)
            Sw(isource) = Sw(isource)/sum_kernel(isource)
        enddo

        ! COMPUTE FORCES
        ! From here on is the same as the poitnwise function except that here the kernel is already computed
        ! Probably it would be a good idea to integrate them.
        ! Zero the Source term at each time step
        FTx(:,:,:)=0.0
        FTy(:,:,:)=0.0
        FTz(:,:,:)=0.0
        SFx(:)=0.0
        SFy(:)=0.0
        SFz(:)=0.0
        Visc=xnu
        !## Send the velocities to the
        call set_vel
        !## Compute forces
        call actuator_line_model_compute_forces
        !## Get Forces
        call get_forces

        ! Loop through all the nodes of the domain
        ! Each process through theirs
        do k=1,xsize(3)
            do j=1,xsize(2)
                do i=1,xsize(1)
                    xmesh = (i + xstart(1) - 1)*dx
                    ymesh = (j + xstart(2) - 1)*dy
                    zmesh = (k + xstart(3) - 1)*dz

                    do isource=1,Nsource
                        ! Distance from the node to the AL point
                        dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
                        ! Gaussian Kernel
                        ! I see this dangerous anyway
                        if(dist.lt.(extended_cells*min(dx,dy,dz))) then
                            Kernel= 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)
                            FTx(i,j,k)=FTx(i,j,k)-SFx(isource)*Kernel/sum_kernel(isource)
                            FTy(i,j,k)=FTy(i,j,k)-SFy(isource)*Kernel/sum_kernel(isource)
                            FTz(i,j,k)=FTz(i,j,k)-SFz(isource)*Kernel/sum_kernel(isource)
                        endif
                    enddo
                enddo
            enddo
        enddo

        ! sum_fx_al = 0.0
        ! sum_fy_al = 0.0
        ! sum_fz_al = 0.0
        ! sum_fx_grid_part = 0.0
        ! sum_fy_grid_part = 0.0
        ! sum_fz_grid_part = 0.0
        !
        ! do isource=1,Nsource
        !     sum_fx_al = sum_fx_al + SFx(isource)
        !     sum_fy_al = sum_fy_al + SFy(isource)
        !     sum_fz_al = sum_fz_al + SFz(isource)
        ! enddo
        !
        ! do k=1,xsize(3)
        !     do j=1,xsize(2)
        !         do i=1,xsize(1)
        !             sum_fx_grid_part = sum_fx_grid_part + FTx(i,j,k)
        !             sum_fy_grid_part = sum_fy_grid_part + FTy(i,j,k)
        !             sum_fz_grid_part = sum_fz_grid_part + FTz(i,j,k)
        !         enddo
        !     enddo
        ! enddo
        !
        ! call MPI_ALLREDUCE(sum_fx_grid_part,sum_fx_grid,1,MPI_REAL8,MPI_SUM, &
        !     MPI_COMM_WORLD,ierr)
        ! call MPI_ALLREDUCE(sum_fy_grid_part,sum_fy_grid,1,MPI_REAL8,MPI_SUM, &
        !     MPI_COMM_WORLD,ierr)
        ! call MPI_ALLREDUCE(sum_fz_grid_part,sum_fz_grid,1,MPI_REAL8,MPI_SUM, &
        !     MPI_COMM_WORLD,ierr)
        !
        ! write(*,*) "Check fx", sum_fx_al, sum_fx_grid
        ! write(*,*) "Check fy", sum_fy_al, sum_fy_grid
        ! write(*,*) "Check fz", sum_fz_al, sum_fz_grid


    end subroutine Compute_Momentum_Source_Term_integral

end module actuator_line_source
