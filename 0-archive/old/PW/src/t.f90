
    mcharge0=0
    mcharge1=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)

             deltakG=norm2(g(1:2,igk_k(ig1,ik0))&
                      -g(1:2,igk_k(ig2,ik))&
                      +xk(1:2,ik0)-xk(1:2,ik))*tpiba
             epsk= splint(eps_data(1,:),eps_data(2,:),eps_data_dy(:),deltakG)
             if (deltak>0.2)      epsk=minval(eps_data(2,:))
             mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
             mcharge1=mcharge1+mcharge0*tpi/deltakG*epsk
         
      Enddo
    Enddo


    mcharge1=mcharge1/dffts%nnr
    write(*,*)  'Mcharge2DnoLFAes ki->kf'   ,ik0,ik,   mcharge1, abs(mcharge1)


    mcharge0=0
    mcharge1=0
    mcharge2=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
    
        mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
        deltakG=norm2(g(:,igk_k(ig1,ik0))&
                   -g(:,igk_k(ig2,ik))&
                   +xk(:,ik0)-xk(:,ik))*tpiba
    
        qxy=norm2(g(1:2,igk_k(ig1,ik0))&
                   -g(1:2,igk_k(ig2,ik))&
                   +xk(1:2,ik0)-xk(1:2,ik))*tpiba
        qz= ((g(3,igk_k(ig1,ik0))-g(3,igk_k(ig2,ik))+ &
             xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
        epsk= splint(eps_data(1,:),eps_data(2,:),eps_data_dy(:),deltakG)
        if (deltakG>0.2) epsk=minval(eps_data(2,:))

        mcharge1=mcharge1+mcharge0*4*pi/(deltakG**2)*epsk
        mcharge2=mcharge2+mcharge0*4*pi/(deltakG**2)*epsk&
          *(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
 
      Enddo
    Enddo
    
    mcharge1=mcharge1/dffts%nnr
    mcharge2=mcharge2/dffts%nnr
    write(*,*)  'Mcharge3DcutnoLFAes ki->kf',ik0,ik,   mcharge1, abs(mcharge1)
    write(*,*)  'Mcharge3DnoLFAes ki->kf'   ,ik0,ik,   mcharge2, abs(mcharge2)


