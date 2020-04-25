    module Rsld
    use CosmologyTypes
    use MatrixUtils
    use likelihood
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: RsldLikelihood
        real(mcp):: Rsld_zd,Rsld_zs
        real(mcp) :: lambdaDdt,mnuDdt,sigmaDdt
        real(mcp) :: lambdaDd,mnuDd,sigmaDd
        character(len=10):: distanceType
        
!double precision :: Rsld_Ninv(Rsld_num,Rsld_Num)
    contains
    procedure :: logLikeTheory => Rsld_LnLike
    end type RsldLikelihood

    public RsldLikelihood,RsldLikelihood_Add
    contains

    subroutine RsldLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(RsldLikelihood), pointer :: this
    Type(TTextFile) :: F
    character (len=20) name

    allocate(this)
    if (.not. Ini%Read_Logical('use_Rsld',.false.)) return
    this%distanceType=Ini%Read_String('distanceType')
    this%LikelihoodType = 'Rsld_1608'
    !this%name= Ini%Read_String('Rsld_type')
    this%name= 'Rsld_1608'
    this%needs_background_functions = .true.
    call LikeList%Add(this)

    if (Feedback > 0) write (*,*) 'Reading: Rsld J1608 Ddt and Dd data...'
    
    call F%Open(trim(DataDir)//'gltdDataset/gltdData/B1608_data.txt')
    read (F%unit,*) name,this%Rsld_zd,this%Rsld_zs,this%lambdaDdt,this%mnuDdt,this%sigmaDdt,&
                             this%lambdaDd,this%mnuDd,this%sigmaDd
    !    write(*,*) this%Rsld_zd(i),this%Rsld_zs(i),this%lambdaDdt(i),this%mnuDdt(i),this%sigmaDdt(i),&
    !                         this%lambdaDd(i),this%mnuDd(i),this%sigmaDd(i)
    call F%Close()

    end subroutine RsldLikelihood_Add


    function Rsld_LnLike(this, CMB)
    !Assume this is called just after CAMB with the correct model  use camb
    Class(CMBParams) CMB
    Class(RsldLikelihood) :: this
    real(mcp) Rsld_LnLike
    integer i
    double precision zd,zs,Ds,Dl,Dds,ddt
    real(mcp) Lnp1,Lnp2

    zd= this%Rsld_zd
    zs=this%Rsld_zs
    Dl =  this%Calculator%AngularDiameterDistance(zd)
    Ds =  this%Calculator%AngularDiameterDistance(zs)
!    Ds =  ComovingDiameterDistance(0.d0,zs)/(1+zs)
    Dds =  this%Calculator%AngularDiameterDistance2(zd,zs)

    ddt = (1.d0+zd)*Dl*Ds/Dds

! minusLnprob
    Lnp1 =&

     !   likelihood from Ddt
               (log(ddt-this%lambdaDdt)-this%mnuDdt)**2/2/this%sigmaDdt**2 + &
                this%sigmaDdt**2/2.0-this%mnuDdt+log(ddt-this%lambdaDdt) !+ &
               
       !   likelihood from dd
    Lnp2 = &
                  (log(Dl-this%lambdaDd)-this%mnuDd)**2/2/this%sigmaDd**2 + &
                      this%sigmaDd**2/2.0-this%mnuDd+log(Dl-this%lambdaDd)
    if (this%distanceType == 'ddt') then
        Rsld_LnLike = Lnp1
    elseif (this%distanceType == 'dd') then
        Rsld_LnLike = Lnp2
    elseif (this%distanceType == 'ddtdd') then
        Rsld_LnLike = Lnp2+Lnp1
    endif

!    write(*,*) '**********this is Rsld_LnLike********',Rsld_LnLike
    if (Feedback > 1) write (*,*) 'Rsld chisq: ', 2.0*Rsld_LnLike
   
    end function Rsld_LnLike

    end module Rsld
