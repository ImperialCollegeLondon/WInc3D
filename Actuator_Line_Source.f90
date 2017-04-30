
module actuator_line_source

    use decomp_2d, only: mytype
    use actuator_line_model_utils
    use actuator_line_model

    implicit none
    real(mytype),save :: constant_epsilon, meshFactor, thicknessFactor,chordFactor
    real(mytype),save, allocatable :: Sx(:),Sy(:),Sz(:),Sc(:),Se(:),Sh(:),Su(:),Sv(:),Sw(:),SFX(:),SFY(:),SFZ(:)
    real(mytype),save, allocatable :: Su_part(:),Sv_part(:),Sw_part(:)
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
    allocate(Sx(NSource),Sy(NSource),Sz(NSource),Sc(Nsource),Su(NSource),Sv(NSource),Sw(NSource),Se(NSource),Sh(NSource),Sfx(NSource),Sfy(NSource),Sfz(NSource))
    allocate(Su_part(NSource),Sv_part(NSource),Sw_part(NSource))
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
    
    subroutine Compute_Momentum_Source_Term_pointwise

        use decomp_2d, only: mytype, xstart, xend, xsize
        use MPI
        use param 
        use variables, only: yp
        use var, only: ux1, uy1, uz1, FTx, FTy, FTz
        
        implicit none
        real(mytype) :: xmesh, ymesh,zmesh
        real(mytype) :: dist, epsilon, Kernel
        real(mytype) :: min_dist, ymax,ymin,zmin,zmax 
        integer :: min_i,min_j,min_k
        integer :: i,j,k, isource, ierr

        ! First we need to compute the locations
        call get_locations 
        
        ! Zero the velocities
        Su(:)=0.0
        Sv(:)=0.0
        Sw(:)=0.0
        ! This is not optimum but works
        ! Define the domain
        ymin=xstart(2)*dy-dy
        ymax=xend(2)*dy
        zmin=xstart(3)*dz-dz
        zmax=xend(3)*dz
        
        do isource=1,NSource
        
        min_dist=1e6
        if((Sy(isource)>=ymin).and.(Sy(isource)<=ymax).and.(Sz(isource)>=zmin).and.(Sz(isource)<=zmax)) then
            !write(*,*) 'Warning: I own this node'
            do k=1,xsize(3)
            zmesh=(k+xstart(3)-1)*dz 
            do j=1,xsize(2)
            if (istret.eq.0) ymesh=(j+xstart(2)-1-1)*dy
            if (istret.ne.0) ymesh=yp(j+xstart(2)-1)
            do i=1,xsize(1)
            xmesh=(i-1)*dx
            dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2) 
            if (dist<min_dist) then
                min_dist=dist
                min_i=i
                min_j=j
                min_k=k
            endif 
            enddo
            enddo
            enddo

            Su_part(isource)=ux1(min_i,min_j,min_k)
            Sv_part(isource)=uy1(min_i,min_j,min_k)
            Sw_part(isource)=uz1(min_i,min_j,min_k)
        else
            Su_part(isource)=0.0
            Sv_part(isource)=0.0
            Sw_part(isource)=0.0
            !write(*,*) 'Warning: I do not own this node' 
        endif
        enddo
           
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

        do isource=1,NSource
            !## Add the source term
            do k=1,xsize(3)
            zmesh=(k+xstart(3)-1)*dz
            do j=1,xsize(2)
            if (istret.eq.0) ymesh=(j+xstart(2)-1-1)*dy
            if (istret.ne.0) ymesh=yp(j+xstart(2)-1) 
            do i=1,xsize(1)
            xmesh=(i-1)*dx
            dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
            epsilon=2.0*(dx*dy*dz)**(1.0/3.0) 
            if (dist<10.0*epsilon) then
                Kernel= 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)            
            else
                Kernel=0.0
            endif
            ! First apply a constant lift to induce the 
            FTx(i,j,k)=FTx(i,j,k)-SFx(isource)*Kernel
            FTy(i,j,k)=FTy(i,j,k)-SFy(isource)*Kernel
            FTz(i,j,k)=FTz(i,j,k)-SFz(isource)*Kernel
            enddo
            enddo
            enddo
        enddo
        end subroutine Compute_Momentum_Source_Term_pointwise

    !subroutine Compute_Momentum_Source_Term_RBF
 
    !    implicit none
    !    integer :: counter,itur,iblade,ielem,jelem,ial
    !    real    :: dx,dy,dz,d
    !    real,allocatable(:) :: Fx,Fy,Fz
    !    
    !    counter=0
    !    if(Ntur>0) then
    !        do itur=1,Ntur
    !            !Blades>
    !            do iblade=1,Turbine(itur)%Nblades
    !                !> Form Matrix A_rbf
    !                allocate(Fx(Turbine(itur)%Blade(iblade)%Nelem),Fy(Turbine(itur)%Blade(iblade)%Nelem),Fx(Turbine(itur)%Blade(iblade)%Nelem)
    !                do jelem=1,Turbine(itur)%Blade(iblade)%Nelem
    !                    do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
    !                    dx=Turbine(itur)%Blade(iblade)%PEx(ielem)-Turbine(itur)%Blade(iblade)%PEy(jelem)
    !                    dy=Turbine(itur)%Blade(iblade)%PEy(ielem)-Turbine(itur)%Blade(iblade)%PEy(jelem)
    !                    dz=Turbine(itur)%Blade(iblade)%PEz(ielem)-Turbine(itur)%Blade(iblade)%PEz(jelem) 
    !                    d=sqrt(dx*dx+dy*dy+dz*dz)
    !                    Turbine(itur)%Blade(iblade)%A_rbf(ielem,jelem)=IsoKernel(d,Turbine(itur)%Blade(iblade)%Eepsilon(ielem),3)
    !                    enddo
    !                enddo
    !
    !                !> Compute the forces by solving 
    !                call solve(Turbine(itur)%Blade(iblade)%A_rfb,fx)
    !                call solve(Turbine(itur)%Blade(iblade)%A_rfb,fy)
    !                call solve(Turbine(itur)%Blade(iblade)%A_rfb,fz)
    !                
    !                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
    !                counter=counter+1
    !                Sfx(counter)=fx
    !                Sfy(counter)=fy
    !                Sfz(counter)=fz    
    !                enddo
    !                deallocate(fx,fy,fz)
    !            enddo
    !        enddo
    !    endif
    !end subroutine Compute_Momentum_Source_Term_RBF
    
end module actuator_line_source
