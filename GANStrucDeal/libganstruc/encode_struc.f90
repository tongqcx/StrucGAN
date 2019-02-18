subroutine encode_struc_2D(na, ns, pos, discut, sigma, epos)
real  :: pos(na,3)
real  :: discut
real  :: sigma
real  :: epos(1,ns,ns)
!f2py  intent(in)      :: na, ns, pos,discut,sigma
!f2py  intent(out)     :: epos

integer,allocatable    :: tpos(:,:)
integer                :: i,j,k,na,ns
real                   :: distance
integer                :: mean_x,mean_y,width_x,width_y

mean_x = 0
mean_y = 0
na = size(pos,1)
ns = size(epos,2)
allocate(tpos(na,2))
!do i = 1,na
!print *, pos(i,:)
!enddo
do i = 1,na
    do j = 1,2
        tpos(i,j) = int(pos(i,j)/discut) 
    enddo
enddo

width_x = maxval(tpos(:,1)) - minval(tpos(:,1))
width_y = maxval(tpos(:,2)) - minval(tpos(:,2))
if (width_x > ns.or.width_y > ns) then
    print *, 'Discut is too small!!!'
endif

mean_x = int(sum(tpos(:,1))/na)
mean_y = int(sum(tpos(:,2))/na)
mean_x = ns/2 - mean_x
mean_y = ns/2 - mean_y
do i = 1,na
    tpos(i,1) = tpos(i,1) + mean_x
    tpos(i,2) = tpos(i,2) + mean_y
enddo

epos = 0.d0
do i = 1,ns
    do j = 1,ns
        do k = 1,na
            distance = (i - tpos(k,1))**2 + (j - tpos(k,2))**2
            distance = -1*sqrt(distance)/sigma
            epos(1,i,j) = epos(1,i,j) + exp(distance)
        enddo
    enddo
enddo
!open(1111,file='wwwww')
!do i = 1,ns
!    do j = 1,ns
!write(1111,'(F10.5,$)'), epos(1,i,j)
!enddo
!write(1111,*)
!enddo
!close(1111)
deallocate(tpos)
end subroutine

subroutine encode_struc_3D(na, ns, pos, discut, sigma, epos)
real  :: pos(na,3)
real  :: discut
real  :: sigma
real  :: epos(1,ns,ns,ns)
!f2py  intent(in)      :: na, ns, pos,discut,sigma
!f2py  intent(out)     :: epos

integer,allocatable    :: tpos(:,:)
integer                :: i,j,k,na,ns,i1,i2,i3
real                   :: distance
integer                :: mean_x,mean_y,mean_z,width_x,width_y,width_z

mean_x = 0
mean_y = 0
mean_z = 0
na = size(pos,1)
ns = size(epos,2)
allocate(tpos(na,3))
!do i = 1,na
!print *, pos(i,:)
!enddo
do i = 1,na
    do j = 1,3
        tpos(i,j) = int(pos(i,j)/discut) 
    enddo
enddo

width_x = maxval(tpos(:,1)) - minval(tpos(:,1))
width_y = maxval(tpos(:,2)) - minval(tpos(:,2))
width_z = maxval(tpos(:,3)) - minval(tpos(:,3))
if (width_x > ns .or. width_y > ns .or. width_z > ns) then
    print *, 'Discut is too small!!!'
endif

mean_x = int(sum(tpos(:,1))/na)
mean_y = int(sum(tpos(:,2))/na)
mean_z = int(sum(tpos(:,3))/na)
mean_x = ns/2 - mean_x
mean_y = ns/2 - mean_y
mean_z = ns/2 - mean_z

do i = 1,na
    tpos(i,1) = tpos(i,1) + mean_x
    tpos(i,2) = tpos(i,2) + mean_y
    tpos(i,3) = tpos(i,3) + mean_z
enddo

epos = 0.d0
do i1= 1,ns
    do i2 = 1,ns
        do i3 = 1,ns
            do k = 1,na
                distance = (i1 - tpos(k,1))**2 + (i2 - tpos(k,2))**2 + (i3 - tpos(k,3))**2
                distance = -1*sqrt(distance)/sigma
                epos(1, i1, i2, i3) = epos(1, i1, i2, i3) + exp(distance)
            enddo
        enddo
    enddo
enddo

deallocate(tpos)
end subroutine

subroutine decode_struc_3D(width, discut, na, x, k, xyz)
integer   ::  width
integer   ::  na
real      ::  x(width, width, width)
real      ::  xyz(na,3)
real      ::  discut
!integer,parameter  ::  n_sample = 100
integer            ::  n_sample
!integer            ::  record_index(n_sample,3)
integer,allocatable  ::  record_index(:,:)
integer            ::  i,j,k,n,m,l
logical            ::  lexit
!f2py  intent(in)      :: width, na, x, discut
!f2py  intent(out)     :: xyz, k

n_sample = int(real(na)*1.5)
allocate(record_index(n_sample,3))

do i = 1,n_sample
    record_index(i,:) = maxloc(x)
    n = record_index(i,1)
    m = record_index(i,2)
    l = record_index(i,3)
    x(n,m,l) = 0.0
enddo
!record_index = real(record_index) * discut
xyz = 0.0
xyz(1,:) = record_index(1,:)
k = 1
lexit = .FALSE.
do i = 2,n_sample
     do j = 1,k
         !print*, record_index(i,:),record_index(j,:),length(record_index(i,:),record_index(j,:)),0.5/discut
         if (length(record_index(i,:),record_index(j,:)) < (1.0/discut))  then
             lexit = .True.
             exit
         endif
     enddo
     if (lexit)  then
         cycle
     endif
     k = k + 1
     xyz(k,:) = record_index(i,:)
     if (k == na) then
         exit
     endif
enddo
xyz = xyz * discut
deallocate(record_index)
contains
function length(x,y)
integer,intent(in)  :: x(:),y(:)
real             :: length
integer          :: i,n
n = size(x,1)
length = 0.0
do i = 1,n
    length = length + (x(i) - y(i))**2
enddo
length = real(sqrt(length))
end function
END SUBROUTINE
!-----------------------------------------------------------------
!  decoding structure in 2D
!-----------------------------------------------------------------
subroutine decode_struc_2D(width, discut, na, x, k, xyz)
integer   ::  width
integer   ::  na
real      ::  x(width, width)
real      ::  xyz(na,3)
real      ::  discut
real      ::  record_max
integer,parameter  ::  n_sample = 100
integer            ::  record_index(n_sample,2)
integer            ::  i,j,k,n,m,kk
logical            ::  lexit
!f2py  intent(in)      :: width, na, x, discut
!f2py  intent(out)     :: xyz,k
kk = 0
do i = 1,n_sample
    record_index(i,:) = maxloc(x)
    record_max = maxval(x)
    n = record_index(i,1)
    m = record_index(i,2)
    !if (record_max < 0.75)  then
    !    print*, n,m, 'less than 0.75' 
    !    kk = kk + 1
    !endif
    !l = record_index(i,3)
    x(n,m) = 0.0
enddo
!print*, 'kk',n_sample - kk
!record_index = real(record_index) * discut
xyz = 0.0
xyz(1,1:2) = record_index(1,:)
k = 1
lexit = .FALSE.
do i = 2,n_sample
     do j = 1,k
         !print*, record_index(i,:),record_index(j,:),length(record_index(i,:),record_index(j,:)),0.5/discut
         if (length(record_index(i,:),record_index(j,:)) < (1.0/discut))  then
             lexit = .True.
             exit
         endif
     enddo
     if (lexit)  then
         cycle
     endif
     k = k + 1
     xyz(k,1:2) = record_index(i,:)
     if (k == na) then
         exit
     endif
enddo
xyz = xyz * discut

contains
function length(x,y)
integer,intent(in)  :: x(:),y(:)
real             :: length
integer          :: i,n
n = size(x,1)
length = 0.0
do i = 1,n
    length = length + (x(i) - y(i))**2
enddo
length = real(sqrt(length))
end function

end subroutine
!-----------------------------------------------------------------
!  decoding structure in 2D crystal structure
!-----------------------------------------------------------------
subroutine decode_struc_2D_crystal(width, discut, na, x, k, xyz)
integer   ::  width
integer   ::  na
real      ::  x(width, width)
real      ::  xyz(na,3)
real      ::  discut
real      ::  record_max
integer   ::  n_sample
integer            ::  record_index(300,2)
!integer            ::  i,j,k,n,m,kk
integer            ::  i,k,n,m,kk
logical            ::  lexit
!f2py  intent(in)      :: width, na, x, discut
!f2py  intent(out)     :: xyz,k
kk = 0
i = 0
do while (.TRUE.)
    i = i + 1
    record_index(i,:) = maxloc(x)
    record_max = maxval(x)
    if (record_max < 0.99) then
        exit
    endif
    n = record_index(i,1)
    m = record_index(i,2)
    x(n,m) = 0.0
enddo
n_sample = i
!print*, 'XXXX',i
xyz = 1.0
xyz(1,1:2) = record_index(1,:)
k = 1
lexit = .FALSE.
do i = 2,n_sample
     !do j = 1,k
     !    !print*, record_index(i,:),record_index(j,:),length(record_index(i,:),record_index(j,:)),0.5/discut
     !    if (length(record_index(i,:),record_index(j,:)) < (1.0/discut))  then
     !        lexit = .True.
     !        exit
     !    endif
     !enddo
     !if (lexit)  then
     !    cycle
     !endif
     k = k + 1
     xyz(k,1:2) = record_index(i,:)
     !if (k == na) then
     !    exit
     !endif
enddo
xyz = xyz * discut

contains
function length(x,y)
integer,intent(in)  :: x(:),y(:)
real             :: length
integer          :: i,n
n = size(x,1)
length = 0.0
do i = 1,n
    length = length + (x(i) - y(i))**2
enddo
length = real(sqrt(length))
end function

end subroutine
