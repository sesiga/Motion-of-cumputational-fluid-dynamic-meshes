!******************************************************************************
! TopologicalDataStructures
!******************************************************************************
!
! Description:
!   Tipos derivados para almacenar la topologia de la malla
!------------------------------------------------------------------------------
module TopologicalDataStructures

    use constants, only :   &
        LENGTH

    !---------------------------------------------------------------------------------------------
    ! topologia
    !---------------------------------------------------------------------------------------------
    type meshelement
        integer :: ifield                              ! field al que pertenece el elemento
        integer :: ielement                            ! dentro del bloque iblock, que elemento es
        character(LENGTH) :: elementType               ! topological element type
        integer :: nnodes                              ! no. de nodos que lo componen
        integer :: nedges                              ! no. de aristas que lo componen
        integer :: nfacets                             ! no. de caras que lo componen
        integer, pointer :: nodes(:)                    ! nodos que forman el elemento
        integer, pointer :: edges(:)                    ! aristas que forman el elemento
        integer, pointer :: facets(:)                   ! caras que forman el elemento
    end type meshelement

    type meshfacet
        integer :: nnodes 
        integer :: nedges 
        integer :: nelements 
        integer, pointer :: nodes(:)                    ! vertices que forman la cara
        integer, pointer :: edges(:)                    ! aristas que forman la cara
        integer, pointer :: elements(:)                 ! elementos a los que pertence la cara
    end type meshfacet

    type meshedge
        integer :: nelements
        integer :: nfacets
        integer :: nodes(2)                             ! vertices que forman la arista
        integer, pointer :: facets(:)                   ! caras a las que pertenece la arista
        integer, pointer :: elements(:)                 ! elementos a los que pertenece la arista
    end type meshedge

    type meshnode
        integer :: Nedges
        integer :: Nfacets
        integer :: Nelements
        integer, pointer :: edges(:)                    ! aristas a las que pertenece el nodo
        integer, pointer :: facets(:)                   ! caras a las que pertenece el nodo
        integer, pointer :: elements(:)                 ! elementos a los que pertenece el nodo
    end type meshnode 
 
    type topologymeshclass
        integer :: nelements
        integer :: nfacets
        integer :: nedges
        integer :: nnodes
        type(meshelement), pointer :: Elements(:)
        type(meshfacet), pointer   :: Facets(:)
        type(meshedge), pointer    :: Edges(:)
        type(meshnode), pointer    :: Nodes(:)
    end type topologymeshclass 
 
end module TopologicalDataStructures

