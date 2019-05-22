module NEQ_SELF_CONSISTENCY
  USE NEQ_VARS_GLOBAL
  implicit none

  private

contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : Initialize the run guessing/reading self-energy
  !+-------------------------------------------------------------------+
  subroutine init_equilibrium_self(self,params)
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    real(8)                 :: wm,res,ims
    logical                 :: bool
    integer                 :: i,j,ik,unit,len

    if(.not.self%status)stop "init_self: self is not allocated"
    if(.not.params%status)stop "init_self: params is not allocated"

    !CHECK IF SIGMA(IW) IS AVAILABLE (START FROM THE EQUILIBRIUM SOLUTION)
    !IF NOT, START FROM NON-INTERACTING (SIGMA=0)
    bool = inquire_keldysh_contour_gf(trim(smats_file))
    if(bool)then
       write(*,"(A)")"Reading initial Sigma(iw) from file "//txt(trim(smats_file))
       unit = free_unit()
       open(unit,file=trim(smats_file),statud='old')
       Lfreq= file_length(unit)
       if(allocated(self%iw))deallocate(self%iw) !if(associated(self%iw))nullify(self%iw)
       allocate(self%iw(Lfreq))
       do i=1,Lfreq
          read(unit,*)wm,res,ims
          self%iw(i) = dcmplx(res,ims)
       enddo
       close(unit)
       !NOW YOU NEED TO PROPAGATE THE SOLUTION to Sigma^{x=M,<,R,\lmix}
       call fftgf_iw2tau(self%iw,self%mats(0:),params%beta)
       self%less(1,1) = -xi*self%mats(params%Ntau)
       self%ret(1,1) =  xi*(self%mats(0)+self%mats(params%Ntau))
       forall(i=0:params%Ntau)self%lmix(1,i)=-xi*self%mats(params%Ntau-i)
    else
       write(*,"(A)")"Start from Non-interacting/Hartree-Fock Sigma(iw)"
       ! if(allocated(self%iw))deallocate(self%iw) !if(associated(self%iw))nullify(self%iw)
       ! allocate(self%iw(Lfreq))
       ! self%iw=zero             !or whatever corresponds to non-interacting hartree-fock condition
       self=zero                !set self-energy to zero or whatever is the initial HF solution
    endif
  end subroutine init_equilibrium_self




end module NEQ_SELF_CONSISTENCY
