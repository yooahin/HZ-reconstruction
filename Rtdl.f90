    ! module for gltd
    ! NOTE: THE ORIGINAL GLTD DATA FOR OBTAINING COSMOLOGICAL PARAMETER ARE USUALLY GIVEN IN FORM OF MCMC
    ! CHANGE!

    module Rtdl
    use CosmologyTypes
    use MatrixUtils
    use likelihood
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: RtdlLikelihood

        character(len=:), allocatable :: gltd_distanceType
        real(mcp):: Rtdl_zd, Rtdl_zs

        !gltd with ddt and dd
        integer :: n_ddt,n_dd
        real(mcp), allocatable, dimension(:) :: ddt_file, dd_file
        real(mcp), allocatable ::  prob_file(:,:)
        !gltd with ddt or dd 
        integer :: n_d
        real(mcp), allocatable, dimension(:) ::  d_file, prob_file_1d
            
    contains
    procedure :: logLikeTheory => Rtdl_LnLike
    procedure :: ReadIni => GLTD_ReadIni
    procedure :: InitProbDist => GLTD_InitProbDist
    end type RtdlLikelihood

!    Type, extends(RtdlLikelihood) :: Ddtdd_Likelihood
!    contains
!    procedure :: logLikeTheory => Ddtdd_loglike
!    end type


    public RtdlLikelihood_Add


    contains

    subroutine RtdlLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(RtdlLikelihood), pointer :: this
    integer i,j
    Type(TTextFile) :: F
    Type(TSettingIni) :: DataSets, OverrideSettings

    if (.not. Ini%Read_Logical('use_Rtdl',.false.)) return

    call Ini%TagValuesForName('gltd_dataset', DataSets, filename=.true.)
    if (DataSets%Count==0) call MpiStop('Use_gltd but no gltd_dataset[NAMETAG] defined')

    do i = 1,DataSets%Count
        call Ini%SettingValuesForTagName('gltd_dataset',DataSets%Name(i),OverrideSettings)
        if (Datasets%Name(i)=='j1206' .or. Datasets%Name(i)=='pg1115' .or. Datasets%Name(i)=='rxj1131'&
           .or. Datasets%Name(i)=='he0435' .or. Datasets%Name(i)=='wfi2033') then
            allocate(RtdlLikelihood::this)
        else
            call MpiStop('unrecgnised GLTD sample : use_Rtdl')
        end if

        call this%ReadDatasetFile(Datasets%Value(i),OverrideSettings)
        this%tag = Datasets%Name(i)
        this%LikelihoodType = 'Rtdl'
        this%needs_background_functions = .true.
        call LikeList%Add(this)
    end do
    end subroutine RtdlLikelihood_Add

    subroutine GLTD_InitProbDist(this,Ini)
    class(RtdlLikelihood) :: this
    class(TSettingIni) :: Ini
    integer :: i,j,ios
    Type(TTExtFile) F
    real(mcp) :: tmp0,tmp1,tmp2
    character(len=:) , allocatable :: gltd_datafile

    this%gltd_distanceType  = Ini%Read_String('distanceType')

    gltd_datafile = Ini%ReadRelativeFileName('datafile')
    call F%Open(gltd_datafile)

    if (this%gltd_distanceType == 'ddtdd') then
        this%n_ddt = Ini%Read_Int('number_of_ddt')
        this%n_dd = Ini%Read_Int('number_of_dd')
        allocate(this%ddt_file(this%n_ddt),this%dd_file(this%n_dd))
        allocate(this%prob_file(this%n_ddt,this%n_dd))!notice!
        do i=1, this%n_ddt
            do j=1, this%n_dd
                read (F%unit,*,iostat=ios) tmp0,tmp1,tmp2
                if (ios /= 0) call MpiStop('Error reading BAO file')
                this%ddt_file(i)  = tmp0
                this%dd_file(j) = tmp1
                this%prob_file(i,j) = tmp2
               ! write(*,*) tmp0,tmp1,tmp2
            end do
        end do
        !Normalize distribution (so that the peak value is 1.0)
        this%prob_file=this%prob_file / maxval(this%prob_file)
    else if (this%gltd_distanceType == 'ddt'.or. this%gltd_distanceType == 'dt' ) then
        this%n_d = Ini%Read_Int('number_of_d')
        allocate(this%d_file(this%n_d))
        allocate(this%prob_file_1d(this%n_d))
        do i=1, this%n_d
            read (F%unit,*,iostat = ios) tmp0,tmp1
                if (ios /= 0) call MpiStop('Error reading BAO file')
                this%d_file(i)  = tmp0
                this%prob_file_1d(i) = tmp1
               ! write(*,*) tmp0,tmp1
        end do
        this%prob_file_1d=this%prob_file_1d / maxval(this%prob_file_1d)
    end if
    end subroutine GLTD_InitProbDist

    subroutine GLTD_ReadIni(this,Ini)
    class(RtdlLikelihood) :: this
    class(TSettingIni) :: Ini

    if (Feedback > 0) write (*,*) 'Reading: Rtdl data : '//trim(this%name)
    this%gltd_distanceType  = Ini%Read_String('distanceType')
        
    this%Rtdl_zd = Ini%Read_Double('zd')
    this%Rtdl_zs = Ini%Read_Double('zs')
!    write(*,*) this%Rtdl_zd,this%Rtdl_zs
    call this%InitProbDist(Ini)

    end subroutine GLTD_ReadIni

    
    real(mcp) function lagerang(x,x0,y0,x1,y1)
    ! give the probality of x using the two points lagrange interpolation
    real(mcp) :: x0,y0,x1,y1,x
    
    lagerang = (x-x1)/(x0-x1)*y0 + (x-x0)/(x1-x0)*y1
    end function

        
    real(mcp) function Rtdl_LnLike(this, CMB)
    Class(CMBParams) CMB
    Class(RtdlLikelihood) :: this
    integer :: i, j,ii,jj
    double precision zd,zs,Ds,Dl,Dds,ddt
    real(mcp) prob
    real(mcp) d_1
!    real(mcp) lagerang

    zd= this%Rtdl_zd
    zs=this%Rtdl_zs
    !Dl =  ComovingDiameterDistance(0.d0,zd)/(1+zd)
    Dl =  this%Calculator%AngularDiameterDistance(zd)
    Ds =  this%Calculator%AngularDiameterDistance(zs)
!    Ds =  ComovingDiameterDistance(0.d0,zs)/(1+zs)
    Dds =  this%Calculator%AngularDiameterDistance2(zd,zs)
    ddt = (1.d0+zd)*Dl*Ds/Dds
    if (this%gltd_distanceType == 'ddtdd') then 
        !write(*,*)'-------------------------'
        !write(*,*) zd,zs,this%Dd_file(1),this%ddt_file(1),this%prob_file(1,1),this%Dd_file(1),this%ddt_file(2),this%prob_file(1,2)
        !write(*,*)'-------------------------'
        if ((Dl < this%Dd_file(1)).or.(Dl > this%Dd_file(this%n_dd-1)).or. &
            &   (ddt < this%ddt_file(1)).or.(ddt > this%ddt_file(this%n_ddt-1))) then
            Rtdl_LnLike = logZero
        else
            do i=1,this%n_dd
                if (Dl - this%Dd_file(i) .le. 0) then
                    ii = i-1
                    exit
                end if
            end do
            do j=1,this%n_ddt
                if (ddt - this%ddt_file(j) .le. 0) then
                    jj = j-1
                    exit
                end if
            end do
            prob=(1./((this%Dd_file(ii+1)-this%Dd_file(ii))*(this%ddt_file(jj+1)-this%ddt_file(jj))))*  &
            &       (this%prob_file(jj,ii)*(this%Dd_file(ii+1)-Dl)*(this%ddt_file(jj+1)-ddt) &
            &       -this%prob_file(jj,ii+1)*(this%Dd_file(ii)-Dl)*(this%ddt_file(jj+1)-ddt) &
            &       -this%prob_file(jj+1,ii)*(this%Dd_file(ii+1)-Dl)*(this%ddt_file(jj)-ddt) &
            &       +this%prob_file(jj+1,ii+1)*(this%Dd_file(ii)-Dl)*(this%ddt_file(jj)-ddt))
            if  (prob > 0) then
                !Rtdl_LnLike = -log( prob )
                Rtdl_LnLike = -log( prob )
            else
                !Rtdl_LnLike = logZero
                Rtdl_LnLike = logZero
            endif
        endif
    else
!------------for 1 distance-------------------------------
        !write(*,*)'-------------------------'
        !write(*,*) zd,zs,this%d_file(1),this%prob_file_1d(1),this%d_file(2),this%prob_file_1d(2)
        !write(*,*)'-------------------------'
        if (this%gltd_distanceType == 'ddt') then
            d_1 = ddt
        elseif (this%gltd_distanceType == 'dd') then
            d_1 = dl
        end if
   
        if ((d_1 < this%d_file(1)).or.(d_1 > this%d_file(this%n_d-1))) then
            Rtdl_LnLike = logZero
        else
            do j=1,this%n_d
                if (d_1 - this%d_file(j) .le. 0) then
                    jj = j-1
                    exit
                end if
            end do
            prob=lagerang(d_1,this%d_file(jj),this%prob_file_1d(jj),this%d_file(jj+1),this%prob_file_1d(jj+1))
            if  (prob > 0) then
                Rtdl_LnLike = -log( prob )
            else
                Rtdl_LnLike = logZero
            endif
        endif
    endif

    if (Feedback > 1) write (*,*) 'Rtdl (dd only) chisq: ', 2.0*Rtdl_LnLike

    end function Rtdl_LnLike



    end module Rtdl
