  do i=1,Nloc
     iup = iup_index(i+mpiIshift,DimUp)
     idw = idw_index(i+mpiIshift,DimUp)
     !
     mup = Hs(1)%map(iup)
     mdw = Hs(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
     !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
     !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
     if(Norb>1.AND.Jx/=0d0)then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (nup(jorb)==1).AND.&
                   (ndw(iorb)==1).AND.&
                   (ndw(jorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(iorb,mdw,k1,sg1)  !DW
                 call cdg(jorb,k1,k2,sg2) !DW
                 jdw=binary_search(Hs(2)%map,k2)
                 call c(jorb,mup,k1,sg3)  !UP
                 call cdg(iorb,k1,k2,sg4) !UP
                 jup=binary_search(Hs(1)%map,k2)
                 htmp = Jx*sg1*sg2*sg3*sg4
                 j = jup + (jdw-1)*dimup
                 !
                 Hv(i) = Hv(i) + htmp*vt(j)
                 !
              endif
           enddo
        enddo
     endif
     !
     ! PAIR-HOPPING (P-H) TERMS
     !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
     !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
     if(Norb>1.AND.Jp/=0d0)then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (nup(jorb)==1).AND.&
                   (ndw(jorb)==1).AND.&
                   (ndw(iorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(jorb,mdw,k1,sg1)       !c_jorb_dw
                 call cdg(iorb,k1,k2,sg2)      !c^+_iorb_dw
                 jdw = binary_search(Hs(2)%map,k2)
                 call c(jorb,mup,k1,sg1)       !c_jorb_up
                 call cdg(iorb,k1,k2,sg4)      !c^+_iorb_up
                 jup = binary_search(Hs(1)%map,k2)
                 htmp = Jp*sg1*sg2*sg3*sg4
                 j = jup + (jdw-1)*dimup
                 !
                 Hv(i) = Hv(i) + htmp*vt(j)
                 !
              endif
           enddo
        enddo
     endif
     !     
  enddo
