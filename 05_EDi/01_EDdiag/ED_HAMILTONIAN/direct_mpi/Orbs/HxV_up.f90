  do iidw=1,MpiQdw
     do iiup=1,DimUp
        i = iiup + (iidw-1)*DimUp
        call state2indices(i,[DimUps,DimDws],Indices)
        !
        !
        !
        !>H_hyb: hopping terms for a given spin (imp <--> bath)
        do iorb=1,Ns_Ud
           mup = Hs(iorb)%map(Indices(iorb))
           Nups(iorb,:) = Bdecomp(mup,Ns_Orb)
           !
           do kp=1,Nbath
              ialfa=1+kp
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) &
                   .AND. (Nups(iorb,1)==1) .AND. (Nups(iorb,ialfa)==0) )then
                 call c(1,mup,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 Jndices       = Indices
                 Jndices(iorb) = binary_search(Hs(iorb)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*Vin(j)
                 !
              endif
              !
              if( (diag_hybr(1,iorb,kp)/=0d0) &
                   .AND. (Nups(iorb,1)==0) .AND. (Nups(iorb,ialfa)==1) )then
                 call c(ialfa,mup,k1,sg1)
                 call cdg(1,k1,k2,sg2)
                 Jndices       = Indices
                 Jndices(iorb) = binary_search(Hs(iorb)%map,k2)
                 call indices2state(Jndices,[DimUps,DimDws],j)
                 htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                 !
                 Hv(i) = Hv(i) + htmp*Vin(j)
                 !
              endif
           enddo
        enddo
        !
     enddo
  enddo






