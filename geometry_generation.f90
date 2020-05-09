    program geometry_generation

    implicit none
    
    integer::i,j,w,imax,jmax,m,n,z,q
    real, allocatable ::TSTAR(:,:,:),xksi(:,:,:),xeta(:,:,:),yksi(:,:,:),yeta(:,:,:),jaco(:,:,:),etax(:,:,:),ksix(:,:,:),etay(:,:,:)
    real, allocatable ::ksiy(:,:,:),xa(:,:,:),ya(:,:,:),coefksi(:,:,:),coefeta(:,:,:),areax(:,:,:),areay(:,:,:)
    real, allocatable :: x(:,:),y(:,:),t(:,:),told(:,:),tzegond(:,:),vol(:,:)
    real, allocatable :: ax(:),bx(:),cx(:),dx(:),sxx(:),ay(:),by(:),cy(:),dy(:),syy(:),ssy(:),teta(:),ssxx(:),ssyy(:),ssxx2(:),ssyy1(:),sy2(:)
    real:: resu,resut,st,sx,sy,sksi,seta,alpha,cprime,czegon,pie,epsilon,epsilont,bc,fig,qdot,steta,r,L1,L2
    print *, "imax="
    read *, imax
    print *, "jmax="
    read *, jmax
    print *, "time_step="
    read *, st
    print *, "fig="
    read*, fig
    w=1
    allocate (TSTAR(IMAX,JMAX,1000),xksi(imax,jmax,5),xeta(imax,jmax,5),yksi(imax,jmax,5),yeta(imax,jmax,5),jaco(imax,jmax,5))
    allocate(etax(imax,jmax,5),ksix(imax,jmax,5),etay(imax,jmax,5),ksiy(imax,jmax,5),xa(imax,jmax,5),ya(imax,jmax,5),coefksi(imax,jmax,5))
    allocate (coefeta(imax,jmax,5),areax(imax,jmax,5),areay(imax,jmax,5))
    allocate (x(imax,jmax),y(imax,jmax),t(imax,jmax),told(imax,jmax),tzegond(imax,jmax),vol(imax,jmax))
    allocate (ax(IMAX),bx(IMAX),cx(IMAX),dx(IMAX),sxx(imax),ay(IMAX),by(IMAX),cy(IMAX),dy(IMAX),syy(imax),ssy(imax),teta(imax),ssyy(imax),ssxx(imax),ssxx2(jmax),ssyy1(imax),sy2(imax)) 
    
    call initial
    call meshgeometrycoef
    
open(21,file="tecplot_result.dat",status="replace",action="write")
write(21,*)'variables = "X" "Y" "T"'
write(21,*) 'Zone T="Time:',(w-1)*st,',Timeset:',st,'"'
write(21,*) 'strandID=1, SolutionTime=',0
write(21,*) 'i=',imax, ', j=',jmax,', ZONETYPE=Ordered'
write(21,*) 'DATAPACKING=Point'

do j=1,jmax
    do i=1,imax
      write(21,*) x(i,j),y(i,j),t(i,j)
    end do 
end do 
 
close(21)
    
    do 
        told(:,:)=t(:,:)
        do i=1,imax
            do j=1,jmax
                TSTAR(i,j,w)=t(i,j)
            end do
        end do
        w=w+1
        do 
            
            tzegond(:,:)=t(:,:)
            m=m+1
            call isweep
            call jsweep

            call resual
            print *, resu,resut,w*st,w
            if(resu<.00001) exit
        end do

        call resualt
        call animate
        if(resut<0.0001) exit
         
    end do
    

    call  ope

    call outp  
    
    print*,"Done"
    pause
    contains
    !*******************************
    subroutine initial
    do i=1,imax
        do j=1,jmax
            t(i,j)=0.1
        end do
    end do  
    end subroutine  initial    
    !**************************************
    subroutine meshgeometrycoef
    pie=3.141592
    alpha=1.0
    r=1.0
   sx=6.0/(imax-2)
   sy=6.0/(jmax-2)
    x(1,1)=0.0
    y(1,1)=0.0
    !fig=1 ! rhombus
    !fig=2 ! square
    !fig=3 ! rectangle
    !fig=4 ! circle
    !fig=5 ! triangle_1
    !fig=6 ! triangle_2
    !fig=7 ! elliptic_1
    !fig=8 ! elliptic_2
   
    !---------------------------------------
    if (fig==1) then
        do i=2,imax
            x(i,1)=x(i-1,1)+sx*cos(45.0*pie/180.0)

            y(i,1)=y(i-1,1)+sx*sin(45.0*pie/180.0)

        end do
        do i=1,imax
        do j=2,jmax
            x(i,j)=x(i,j-1)+sy*sin(-45.0*pie/180.0)
            y(i,j)=y(i,j-1)+sy*cos(-45.0*pie/180.0)
        end do
        end do
    end if 

    !----------------------------------

    if (fig==2) then 
        do i=2,imax 
            x(i,1)=x(i-1,1)+sx
            y(i,1)=0.0
        end do    
       do i=1,imax 
           do j=2,jmax 
              x(i,j)=x(i,j-1)
              y(i,j)=y(i,j-1)+sy
           end do 
       end do 
    end if

    !---------------------------------

    if (fig==3) then
        steta=(pie)/(imax-2)
        teta(1)=0.0
        teta(2)=steta/2.0
        do i=3,imax-1
            teta(i)=teta(i-1)+steta
        end do
        teta(imax)=teta(imax-1)+steta/2.0
        do i=1,imax
            x(i,1)=(0.5-0.5*cos(teta(i)))
            y(i,1)=(0.5*sin(teta(i)))
            x(i,2)=x(i,1)
            ssy(i)=(1.0-0.5*sin(teta(i)))/(jmax-2)
            y(i,2)=y(i,1)+ssy(i)/2.0
            do j=3,jmax-1
                x(i,j)=x(i,j-1)
                y(i,j)=y(i,j-1)+ssy(i)
            end do
            x(i,jmax)=x(i,jmax-1)
            y(i,jmax)=y(i,jmax-1)+ssy(i)/2.0
        end do
    end if
    !----------------------------------

    if (fig==4) then
        y(1,jmax)=sqrt(2.)/2.
        x(1,jmax)=-sqrt(2.)/2.
        y(1,1)=-sqrt(2.)/2.
        x(1,1)=-sqrt(2.)/2.
        !ssyy(i)=(y(1,jmax)-y(1,1))/(jmax-1)
        !ssyy(1)= sqrt(2.)/(jmax-1) 
        x(imax,jmax)=sqrt(2.)/2.  
        ssxx(jmax)=(x(imax,jmax)-x(1,jmax))/(imax-1)
        do i=1,imax
            x(i,jmax)=x(1,jmax)+(i-1)*ssxx(jmax)
            y(i,jmax)=sqrt(r**2-x(i,jmax)**2)
            y(i,1)=-sqrt(r**2-x(i,jmax)**2)
            ssyy(i)=(y(i,jmax)-y(i,1))/(jmax-1) 
        end do 
        do i=1,imax
            do j=1,jmax  
                y(i,j)=y(i,1)+(j-1)*ssyy(i)
            end do 
        end do
        y(imax,jmax)=sqrt(2.)/2.
        y(imax,1)=-sqrt(2.)/2.
        ssyy(imax)=(y(imax,jmax)-y(imax,1))/(jmax-1)
        do j=1,jmax
            y(imax,j)=y(imax,1)+(j-1)*ssyy(imax)
            x(imax,j)=sqrt(1-y(imax,j)**2)
            x(1,j)=-sqrt(1-y(imax,j)**2)
            ssxx(j)=(x(imax,j)-x(1,j))/(imax-1)
        end do
        do i=1,imax
            do j=1,jmax  
                x(i,j)=x(1,j)+(i-1)*ssxx(j)
            end do 
        end do
    end if
    !----------------------------------

    if (fig==5) then
        sy=sin(pie/3.0)/(jmax-1)
        do j=1,jmax
            ssxx(j)=(1.0-(2.0/sqrt(3.0))*(j-1)*sy)/(imax-1)
        end do

        do j=1,jmax
            do i=1,imax 
                x(i,j)=(i-1)*ssxx(j)+((j-1)*sy)/(sqrt(3.0))
                y(i,j)=(j-1)*sy
            end do
        end do
    end if
    !-------------------------------

    if (fig==6) then
        do i=2,imax
           x(i,1)=x(i-1,1)+sx
         y(i,1)=sin(5*pie*i/imax)      !1
          ! y(i,1)=abs(sin(5*pie*i/imax))  !2 
           !if (mod(i,2)==0) then      !3
           !    y(i,1)=1.0  
           ! else
           !    y(i,1)=2.0
           ! end if    
                
        end do
        
        do i=1,imax
            do j=2,jmax
                x(i,j)=x(i,j-1)
                y(i,j)=y(i,j-1)+sy
            end do
        end do
    end if
    !------------------------------

    if (fig==7) then
        !eccencity of ellipse is: e=0.97
        x(1,jmax)=-0.7
        y(1,jmax)=0.7
        x(1,1)=-0.7
        x(1,1)=-0.7
        x(imax,jmax)=0.7
        y(imax,jmax)=0.7
        x(imax,1)=0.7
        y(imax,1)=-0.7
        ssxx(jmax)=(x(imax,jmax)-x(1,jmax))/(imax-1)
        do i=1,imax
            x(i,jmax)=x(1,jmax)+(i-1)*ssxx(jmax)
           ! y(i,jmax)=sqrt(4*r**2-0.16*x(i,jmax)**2)
            y(i,jmax)=sqrt(r**2+6.25*x(i,jmax)**2)
          !  y(i,1)=-sqrt(4*r**2-0.16*x(i,jmax)**2)
            y(i,1)=-sqrt(r**2+6.25*x(i,jmax)**2)
            ssyy(i)=(y(i,jmax)-y(i,1))/(jmax-1) 
        end do 
        do i=1,imax
            do j=1,jmax  
            y(i,j)=y(i,1)+(j-1)*ssyy(i)
            end do 
        end do
          
        ssyy(imax)=(y(imax,jmax)-y(imax,1))/(jmax-1)
        do j=1,jmax
            y(imax,j)=y(imax,1)+(j-1)*ssyy(imax)
            x(imax,j)=sqrt(36*r**2+0.16*y(imax,j)**2)
            !x(1,j)=-sqrt(1-y(imax,j)**2)
          ! x(1,j)=-sqrt(25*r**2-6.25*y(1,j)**2)
            x(1,j)=-sqrt(36*r**2+0.16*y(1,j)**2)
            ssxx(j)=(x(imax,j)-x(1,j))/(imax-1)
        end do
        do i=1,imax
            do j=1,jmax  
                x(i,j)=x(1,j)+(i-1)*ssxx(j)
            end do 
        end do
    end if

    !-----------------------------------------

    if (fig==8) then
        steta=(pie)/(imax-2)
        teta(1)=0.0
        teta(2)=steta/2.0
        do i=3,imax-1
            teta(i)=teta(i-1)+steta
        end do
        teta(imax)=teta(imax-1)+steta/2
        do i=1,imax
            x(i,1)=(0.5-0.5*cos(teta(i)))
            y(i,1)=(0.5*sin(teta(i)))
            x(i,2)=x(i,1)
            ssy(i)=(2.0-sin(teta(i))-0.5*sin(teta(i)))/(jmax-2)
            y(i,2)=y(i,1)+ssy(i)
            do j=3,jmax
                x(i,j)=x(i,j-1)
                y(i,j)=y(i,j-1)+ssy(i)
            end do
        end do
    end if
    
    !----------------
    !------------------------------------
    sksi=1.0
    seta=1.0
    do i=2,imax-1
        do j=2,jmax-1

        !------------in surface 1--------------

        xksi(i,j,1)=(x(i,j)-x(i-1,j))/sksi
        yksi(i,j,1)=(y(i,j)-y(i-1,j))/sksi
        xeta(i,j,1)=.25*(x(i,j+1)-x(i,j-1)+x(i-1,j+1)-x(i-1,j-1))/seta
        yeta(i,j,1)=.25*(y(i,j+1)-y(i,j-1)+y(i-1,j+1)-y(i-1,j-1))/seta

        !-------------in surface 2------------------
        xksi(i,j,2)=(x(i+1,j)-x(i,j))/sksi
        yksi(i,j,2)=(y(i+1,j)-y(i,j))/sksi
        xeta(i,j,2)=.25*(x(i+1,j+1)-x(i+1,j-1)+x(i,j+1)-x(i,j-1))/seta
        yeta(i,j,2)=.25*(y(i+1,j+1)-y(i+1,j-1)+y(i,j+1)-y(i,j-1))/seta

        !-------------in surface 3-------------------
        xksi(i,j,3)=.25*(x(i+1,j-1)-x(i-1,j-1)+x(i+1,j)-x(i-1,j))/sksi
        yksi(i,j,3)=.25*(y(i+1,j-1)-y(i-1,j-1)+y(i+1,j)-y(i-1,j))/sksi
        xeta(i,j,3)=(x(i,j)-x (i,j-1))/seta
        yeta(i,j,3)=(y(i,j)-y(i,j-1))/seta

        !-----------in surface 4---------------------
        xksi(i,j,4)=.25*(x(i+1,j+1)-x(i-1,j+1)+x(i+1,j)-x(i-1,j))/sksi
        yksi(i,j,4)=.25*(y(i+1,j+1)-y(i-1,j+1)+y(i+1,j)-y(i-1,j))/sksi
        xeta(i,j,4)=(x(i,j+1)-x (i,j))/seta
        yeta(i,j,4)=(y(i,j+1)-y(i,j))/seta

        !---------in center of control volume------------
        xksi(i,j,5)=.5*(x(i+1,j)-x (i-1,j))/sksi
        yksi(i,j,5)=.5*(y(i+1,j)-y(i-1,j))/sksi
        xeta(i,j,5)=.5*(x(i,j+1)-x (i,j-1))/seta
        yeta(i,j,5)=.5*(y(i,j+1)-y(i,j-1))/seta

        !---------jacobian---------------------
        jaco(i,j,1)=xksi(i,j,1)*yeta(i,j,1)-xeta(i,j,1)*yksi(i,j,1)
        jaco(i,j,2)=xksi(i,j,2)*yeta(i,j,2)-xeta(i,j,2)*yksi(i,j,2)
        jaco(i,j,3)=xksi(i,j,3)*yeta(i,j,3)-xeta(i,j,3)*yksi(i,j,3)
        jaco(i,j,4)=xksi(i,j,4)*yeta(i,j,4)-xeta(i,j,4)*yksi(i,j,4)
        jaco(i,j,5)=xksi(i,j,5)*yeta(i,j,5)-xeta(i,j,5)*yksi(i,j,5)

        !----------ksi & eta-----------------------
        ksix(i,j,1)=yeta(i,j,1)/jaco(i,j,1)
        etax(i,j,1)=-yksi(i,j,1)/jaco(i,j,1)
        ksiy(i,j,1)=-xeta(i,j,1)/jaco(i,j,1)
        etay(i,j,1)=xksi(i,j,1)/jaco(i,j,1)

        ksix(i,j,2)=yeta(i,j,2)/jaco(i,j,2)
        etax(i,j,2)=-yksi(i,j,2)/jaco(i,j,2)
        ksiy(i,j,2)=-xeta(i,j,2)/jaco(i,j,2)
        etay(i,j,2)=xksi(i,j,2)/jaco(i,j,2)

        ksix(i,j,3)=yeta(i,j,3)/jaco(i,j,3)
        etax(i,j,3)=-yksi(i,j,3)/jaco(i,j,3)
        ksiy(i,j,3)=-xeta(i,j,3)/jaco(i,j,3)
        etay(i,j,3)=xksi(i,j,3)/jaco(i,j,3)

        ksix(i,j,4)=yeta(i,j,4)/jaco(i,j,4)
        etax(i,j,4)=-yksi(i,j,4)/jaco(i,j,4)
        ksiy(i,j,4)=-xeta(i,j,4)/jaco(i,j,4)
        etay(i,j,4)=xksi(i,j,4)/jaco(i,j,4)

        !--------------surface-------------------
        xa(i,j,1)=.25*(x(i,j)+x(i-1,j)+x(i,j-1)+x(i-1,j-1))
        xa(i,j,2)=.25*(x(i,j)+x(i,j+1)+x(i-1,j+1)+x(i-1,j))
        xa(i,j,3)=.25*(x(i,j)+x(i,j-1)+x(i+1,j)+x(i+1,j-1))
        xa(i,j,4)=.25*(x(i,j)+x(i,j+1)+x(i+1,j)+x(i+1,j+1))

        ya(i,j,1)=.25*(y(i,j)+y(i-1,j)+y(i,j-1)+y(i-1,j-1))
        ya(i,j,2)=.25*(y(i,j)+y(i,j+1)+y(i-1,j+1)+y(i-1,j))
        ya(i,j,3)=.25*(y(i,j)+y(i,j-1)+y(i+1,j)+y(i+1,j-1))
        ya(i,j,4)=.25*(y(i,j)+y(i,j+1)+y(i+1,j)+y(i+1,j+1))

        areax(i,j,1)=-(ya(i,j,2)-ya(i,j,1))
        areay(i,j,1)=xa(i,j,2)-xa(i,j,1)
        areax(i,j,2)=ya(i,j,4)-ya(i,j,3)
        areay(i,j,2)=-(xa(i,j,4)-xa(i,j,3))
        areax(i,j,3)=(ya(i,j,3)-ya(i,j,1))
        areay(i,j,3)=-(xa(i,j,3)-xa(i,j,1))
        areax(i,j,4)=-(ya(i,j,4)-ya(i,j,2))
        areay(i,j,4)=(xa(i,j,4)-xa(i,j,2))
        vol(i,j)=jaco(i,j,5)
        !-----coeff---------------
        coefksi(i,j,1)=ksix(i,j,1)*areax(i,j,1)+ksiy(i,j,1)*areay(i,j,1)
        coefeta(i,j,1)=etax(i,j,1)*areax(i,j,1)+etay(i,j,1)*areay(i,j,1)
        coefksi(i,j,2)=ksix(i,j,2)*areax(i,j,2)+ksiy(i,j,2)*areay(i,j,2)
        coefeta(i,j,2)=etax(i,j,2)*areax(i,j,2)+etay(i,j,2)*areay(i,j,2)
        coefksi(i,j,3)=ksix(i,j,3)*areax(i,j,3)+ksiy(i,j,3)*areay(i,j,3)
        coefeta(i,j,3)=etax(i,j,3)*areax(i,j,3)+etay(i,j,3)*areay(i,j,3)
        coefksi(i,j,4)=ksix(i,j,4)*areax(i,j,4)+ksiy(i,j,4)*areay(i,j,4)
        coefeta(i,j,4)=etax(i,j,4)*areax(i,j,4)+etay(i,j,4)*areay(i,j,4)
        end do
    end do

    end subroutine meshgeometrycoef
    !------------------------------------------------------
    subroutine isweep

    do i=2,imax-1
        ay(:)=0.0
        by(:)=0.0
        cy(:)=0.0
        dy(:)=0.0
        do j=2,jmax-1
            ay(j)=-(.25*(coefeta(i,j,1)+coefeta(i,j,2))+coefeta(i,j,3))/seta
            by(j)=-vol(i,j)/(alpha*st)+(coefksi(i,j,1)-coefksi(i,j,2))/sksi+(coefeta(i,j,3)-coefeta(i,j,4))/seta
            cy(j)=(.25*(coefeta(i,j,1)+coefeta(i,j,2))+coefeta(i,j,4))/seta
            dy(j)=-told(i,j)*vol(i,j)/(alpha*st)-(-coefksi(i,j,1)*t(i-1,j)/sksi+.25*coefeta(i,j,1)*(t(i-1,j+1)-t(i-1,j-1))/seta)
            dy(j)=dy(j)-(coefksi(i,j,2)*t(i+1,j)/sksi+.25*coefeta(i,j,2)*(t(i+1,j+1)-t(i+1,j-1))/(seta))
            dy(j)=dy(j)-(.25*coefksi(i,j,3)*(t(i+1,j)-t(i-1,j)+t(i+1,j-1)-t(i-1,j-1))/sksi)
            dy(j)=dy(j)-(.25*coefksi(i,j,4)*(t(i+1,j)-t(i-1,j)+t(i+1,j+1)-t(i-1,j+1))/sksi)
        end do
        bc=fig
        if (bc==1 .or. bc==2 .or. bc==3 .or. bc==6) then
            ay(1)=0.0
            by(1)=1.0
            cy(1)=0.0
            dy(1)=3.0
            ay(jmax)=0.0
            by(jmax)=1.0
            cy(jmax)=0.0
            dy(jmax)=3.0
        end if
        if (bc==4 .or. bc==5 .or. bc==7 .or. bc==8 ) then
            ay(1)=0.0
            by(1)=1.0
            cy(1)=0.0
            dy(1)=0.0
            ay(jmax)=0.0
            by(jmax)=1.0
            cy(jmax)=0.0
            dy(jmax)=3.0
        end if
        call tdmay
        t(i,:)=syy(:)

    end do
    end subroutine isweep
    !---------------------------------------------------
    subroutine jsweep

    do j=2,jmax-1
        ax(:)=0.0
        bx(:)=0.0
        cx(:)=0.0
        dx(:)=0.0
        do i=2,imax-1

        ax(i)=-(.25*(coefksi(i,j,3)+coefksi(i,j,4))+coefksi(i,j,1))/sksi
        bx(i)=-vol(i,j)/(alpha*st)+(coefksi(i,j,1)-coefksi(i,j,2))/sksi+(coefeta(i,j,3)-coefeta(i,j,4))/seta
        cx(i)=(.25*(coefksi(i,j,3)+coefksi(i,j,4))+coefksi(i,j,2))/sksi
        dx(i)=-vol(i,j)*told(i,j)/(alpha*st)-(.25*coefeta(i,j,1)*(t(i,j+1)-t(i,j-1)+t(i-1,j+1)-t(i-1,j-1))/seta)
        dx(i)=dx(i)-(.25*coefeta(i,j,2)*(t(i,j+1)-t(i,j-1)+t(i+1,j+1)-t(i+1,j-1))/seta)
        dx(i)=dx(i)-(.25*coefksi(i,j,3)*(t(i+1,j-1)-t(i-1,j-1))/sksi-coefeta(i,j,3)*t(i,j-1)/seta)
        dx(i)=dx(i)-(.25*coefksi(i,j,4)*(t(i+1,j+1)-t(i-1,j+1))/sksi+coefeta(i,j,4)*t(i,j+1)/seta)
        end do

        bc=fig
        if (bc==1 .or. bc==2 .or. bc==3 .or. bc==6) then
            ax(1)=0.0
            bx(1)=1.0
            cx(1)=0.0
            dx(1)=0.0
            ax(imax)=0.0
            bx(imax)=1.0
            cx(imax)=0.0
            dx(imax)=0.0
        end if
        if (bc==4 .or. bc==5 .or. bc==7 .or. bc==8) then
            ax(1)=0.0
            bx(1)=1.0
            cx(1)=0.0
            dx(1)=0.0
            ax(imax)=0.0
            bx(imax)=1.0
            cx(imax)=0.0
            dx(imax)=0.0
        end if
        call tdmax
        t(:,j)=sxx(:)

    end do
    end subroutine jsweep
    !-----------------------------------------------------
    subroutine tdmay
    real ::f
    n=0
    do z=2,jmax
        f=ay(z)/by(z-1)
        by(z)=by(z)-cy(z-1)*f
        dy(z)=dy(z)-dy(z-1)*f
    end do
    syy(jmax)=dy(jmax)/by(jmax)
    do n=jmax-1,1,-1
        syy(n)=(dy(n)-cy(n)*syy(n+1))/by(n)
    end do
    return
    end subroutine tdmay
    !---------------------------------------------------
    subroutine tdmax
    real ::f
    n=0
    do z=2,imax
        f=ax(z)/bx(z-1)
        bx(z)=bx(z)-cx(z-1)*f
        dx(z)=dx(z)-dx(z-1)*f
    end do
    sxx(imax)=dx(imax)/bx(imax)
    do n=imax-1,1,-1
        sxx(n)=(dx(n)-cx(n)*sxx(n+1))/bx(n)
    end do
    return
    end subroutine tdmax
    !---------------------------------------------------
    subroutine resual
    resu=0.0
    cprime=0.0
    do i=2,imax-1
        do j=2,jmax-1

        resu=resu+abs(t(i,j)-tzegond(i,j))
        cprime=cprime+1
        end do
    end do
    resu=resu/cprime
    end subroutine resual
    !--------------------------------------------------
    subroutine resualt
    resut=0.0

    czegon=0.0
    do i=2,imax-1
        do j=2,jmax-1  

        resut=resut+abs(t(i,j)-told(i,j))
        czegon=czegon+1
        end do
    end do

    resut=resut/czegon
    end subroutine resualt
    !--------------------------------------------------
    subroutine ope
    open(1, file='mesh.plt')
    open(2, file='valid.plt')
    end subroutine ope
    !------------------------------------------------------
    subroutine outp
    write (1,*) 'VARIABLE= "X","Y"'
    write (1,*)  "ZONE ,j=",JMAX ,"i=",IMAX
    do j=1,jmax
        do i=1,imax

        write (1,*)    x(i,j),y(i,j),t(i,j)
        end do
    end do
    write (2,*) 'VARIABLE= "X"'
    do i=1,imax
        write (2,*)   X(i,jmax/2),t(i,jmax/2)
    end do

end subroutine outp
    !--------------------------
subroutine animate 

open(21,file="tecplot_result.dat",status="old",position="append",action="write")
write(21,*)'variables = "X" "Y" "T"'
write(21,*) 'Zone T="Time:',w*st,',Timeset:',st,'"'
write(21,*) 'strandID=1, SolutionTime=',w*st
write(21,*) 'i=',imax, ', j=',jmax,', ZONETYPE=Ordered'
write(21,*) 'DATAPACKING=Point'


do j=1,jmax
    do i=1,imax
      write(21,*) x(i,j),y(i,j),t(i,j)
    end do 
end do 

end subroutine animate
    
    !------------------------------------------
    !------------------------------------------------------
    end program 
