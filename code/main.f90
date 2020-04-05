!******************************************************************************
! main
!******************************************************************************
!
!------------------------------------------------------------------------------
program main

    !modules used-----------------------------

    ! Precision utility.
    use precision, only: &
        wp

    ! Utilities utility.
    use utilities, only: &
        reporterror

    ! Constants utility.
    use constants, only: &
        LENGTH

    ! Setup utility
    use setup

    ! Resume file
    use ResumeFileUtility

    ! Tecplot file
    use TecplotFileUtility

    ! Scatter file
    use ScatterFileUtility
    
    ! scatter file data class
    use ScatterFileDataClass
    
    !topology class
    use TopologicalDataStructures
    
    !TopologyProcess
    use TopologyProcess3D
    
    !checking
    use checking
    
    !RBFprocess
    use RBFprocess
    
    !deformationUtility
    use deformationUtility
    
    !MeshQuality
    use MeshQuality
    
    !qualityDataClass
    use QualityDataClass

    !end modules used-----------------------------------

    implicit none

    ! Local scalars:
    integer :: ierr
    real(wp) :: time

    ! Instantiate data file target:
    type(TecplotFileClass) :: meshFile
    type(ScatterFileClass) :: deformationFile
    
    character(LENGTH) :: deformationFileName
    character(LENGTH) :: meshFileName
    character(LENGTH) :: outputFileName
    
    type(topologymeshclass) :: meshTopology
    
    type(calidadfields) :: calidad
    
    real(wp), allocatable :: PointsIn(:,:) !nodos deformados del scatter
    real(wp), allocatable :: PointsOut(:,:) !nodos a deformar de la malla
    real(wp), allocatable :: ValuesIn(:,:) !deformacion de los nodos del scatter
    real(wp), allocatable :: ValuesOut(:,:) !deformacion de los nodos a deformar de la malla
    integer, allocatable :: PointsIn_id(:) !global id de los nodos deformados del scatter
    integer, allocatable :: PointsOut_id(:) !global id de los nodos a deformar de la malla
    
    real(wp) :: RadiusCoeff 
    
    procedure(RBF), pointer :: RBFfunction
    
    real(wp) :: t1, t2, t3

!- End of header --------------------------------------------------------------

    !-------------------------------------------------------
    ! Tecplot File Processing Application       
    !-------------------------------------------------------
    
    ! open resume file
    ierr = open_resumefile()

    ! Load Setup
    print '(a)', "Loading the setup data file..."
    call loadSetup("settings.ini", deformationFileName, meshFileName, &
                outputFileName, RBFfunction, RadiusCoeff)

    print '(a)', "...done"
    print *

    ! Load Tecplot file into memory.
    print '(a)', "Loading the Tecplot mesh data file..."
    call loadDataFile(meshFile, meshFileName)
    print '(a)', "...done"
    print *

    ! Load deformation file into memory.
    print '(a)', "  Loading the deformation (Scatter) data file..."
    call Scatter_LoadFileData(deformationFile, deformationFilename)
    print '(a)', "  ...done"
    print *
    
    !allocates arrays
    print '(a)', " Loading data into arrays..."
    call PreProcess(meshFile, deformationFile, PointsOut, &
        ValuesIn, ValuesOut, PointsIn, PointsOut_id, PointsIn_id, calidad)
    print '(a)', " ...done"
    print *
    
    !loads topology
    print *, 'saving topology...'
    call TopologyPreProcess3D(meshFile, meshTopology)
    print *, '...done'
    print *
    
    call CPU_TIME(t1)
    
    !deforation process
    print *, 'deformation process'
    call deformationProcess(meshTopology, PointsIn, PointsOut, ValuesIn, ValuesOut, &
        PointsIn_id, PointsOut_id, RBFfunction, radiusCoeff)
    print *, 'done'
    print *
    
    call CPU_TIME(t2)
    
    t3 = t2 - t1
    
    print *, 'tiempo interpolacion', t3
    print *
    
    !deallocates arrays
    print '(a)', " Deallocating arrays"
    call PostProcess(meshFile, deformationFile, PointsIn, PointsOut, &
        ValuesIn, ValuesOut, PointsOut_id, PointsIn_id)
    print '(a)', " ...done"
    print *
    
    !calidad de la malla
    call Quality(MeshFile, calidad)
    
    ! Save Tecplot file.
    print '(a)', "Saving the tecplot data file..."
    call saveDataFile(meshFile, outputFileName, calidad)
    print '(a)', "...done"
    print *
    
    !deallocates calidad
    call qualityDeallocate(calidad)

    ! Finalize the objects.
    call finalize(meshFile)
    call Scatter_Finalize(deformationFile)
    call TopologyPostProcess3D(meshTopology)
    
    !cpu time
    call cpu_time(time)
    print *, time

    ! Termination protocol.
    ierr = close_resumefile()
    
end program main
