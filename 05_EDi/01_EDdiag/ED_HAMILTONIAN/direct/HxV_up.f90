  do idw=1,DimDw
     do iup=1,DimUp
        mup  = Hs(1)%map(iup)
        nup  = bdecomp(mup,Ns)
        i    = iup + (idw-1)*dimUp
        !
        !
        !> H_imp: Off-diagonal elements, i.e. non-local part. 
        !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                   (nup(jorb)==1) .AND. (nup(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
        !
        !> H_Bath: inter-orbital bath hopping contribution.
        if(bath_type=="replica") then   
           do kp=1,Nbath
              do iorb=1,Norb
                 do jorb=1,Norb
                    !
                    ialfa = getBathStride(iorb,kp)
                    ibeta = getBathStride(jorb,kp)
                    Jcondition = &
                         (dmft_bath%h(1,1,iorb,jorb,kp)/=zero) .AND.&
                         (nup(ibeta)==1) .AND. (nup(ialfa)==0)
                    !
                    if (Jcondition)then
                       call c(ibeta,mup,k1,sg1)
                       call cdg(ialfa,k1,k2,sg2)
                       jup = binary_search(Hs(1)%map,k2)
                       j   = jup + (idw-1)*DimUp
                       htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                       !
                       hv(i) = hv(i) + htmp*vin(j)
                       !
                    endif
                 enddo
              enddo
           enddo
        end if
        !
        !
        !>H_hyb: hopping terms for a given spin (imp <--> bath)
        do iorb=1,Norb
           do kp=1,Nbath
              ialfa=getBathStride(iorb,kp)
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. &
                   (nup(iorb)==1) .AND. (nup(ialfa)==0) )then              
                 call c(iorb,mup,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) .AND. &
                   (nup(iorb)==0) .AND. (nup(ialfa)==1) )then
                 call c(ialfa,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 jup = binary_search(Hs(1)%map,k2)
                 j   = jup + (idw-1)*DimUp
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
        !
        !
     enddo
  enddo





