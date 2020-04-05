!TopologyProcess3D
!
!crea la topologia para mallas 3D
    
module TopologyProcess3D

   !modulos usados
   use constants, only : LENGTH
   
   use precision, only : wp
   
   use TopologicalDataStructures
   
   use TopologicalDataStructuresSupport
   
   use TecplotFileDataClass
   
   use TecplotGridDataClass
   
   use TopologyUtility
   !final modulos usados
   
   !variables locales
   !-------------------------------------
   
   !-------------------------------------
      
    contains
    
    !******************************************************************************
    !TopologyPreProcces3D
    !
    !crea la topologia de la malla
    subroutine TopologyPreProcess3D(mesh, meshT)
    
       !subroutine arguments
       type(TecplotFileClass), intent(in) :: mesh
       type(TopologyMeshClass), intent(out) :: meshT
       !end subroutine arguments
       
       !local scalars
       type(EdgeType), pointer :: headEdge
       type(EdgeType), pointer :: currentEdge
       
       type(FacetType), pointer :: headFacet
       type(FacetType), pointer :: currentFacet
       
       type(ElementType), pointer :: headElement(:)
       
       type(NodeType), pointer :: headNode(:)
       type(NodeType_elements), pointer :: currentNode_elements
       
       integer :: ElementNnodes !numero de nodos del elemento
       integer :: ElementNedges !numero de aristas del elemento
       integer :: ElementNfacets !numero de caras del elemento
       character(LENGTH) :: ElementType !tipo de elemento
       integer :: ElementIndex !indice del elemento: meshT%elements(ElementIndex)
       
       integer :: nodeLocal_id !indice local del nodo dentro de cada field
       integer :: nodeGlobal_id !indice global del nodo
       integer, allocatable :: nodes_global_id(:) !vector con el global id de los nodos del elemento
       
       integer :: i, j, k, m
       
       !end of header-----------------------------------------
       
       !allocata los tipos derivados y current los apunta
       allocate(headEdge)
       currentEdge => headEdge
       
       allocate(headFacet)
       currentFacet => headFacet
       
       allocate(headElement(mesh%gridData%no_of_elements))
       
       allocate(headNode(mesh%gridData%no_of_points))
       headNode(:)%Nedges = 0
       headNode(:)%Nfacets = 0
       headNode(:)%Nelements = 0
       do i=1,mesh%gridData%no_of_points
           allocate(headNode(i)%edges)
           allocate(headNode(i)%elements)
       end do
       
       !inicializa los valores
       !----------------------------------------
       ElementIndex = 0 
    
       meshT%Nedges = 0
       meshT%Nfacets = 0
       !---------------------------------------
       
       !meshT: carga numero de elementos y de nodos
       !y allocata elementos y nodos
       !----------------------------------------
       meshT%Nelements = mesh%GridData%no_of_elements
       meshT%Nnodes = mesh%GridData%no_of_points
       
       allocate(meshT%elements(meshT%Nelements))
       allocate(meshT%nodes(meshT%Nnodes))
       
       meshT%nodes(:)%Nedges = 0
       meshT%nodes(:)%Nelements = 0
       !----------------------------------------
       
       !recorre todos los fields y todos los elementos
       !------------------------------------------
       do i=1,size(mesh%grid%fields)
           
           !guarda algunos valores de element para utilizarlos posteriormente
           !--------------------------------
           ElementType = mesh%grid%fields(i)%elementType
           
           select case(ElementType)
           case('triangle')
               ElementNnodes = 3
               ElementNedges = 3
               ElementNfacets = 0
           case('quadrilateral')
               ElementNnodes = 4
               ElementNedges = 4
               ElementNfacets = 0
           end select
           !---------------------------------
           
           do j=1,mesh%grid%fields(i)%fieldHeader%no_of_elements
               
               !siguiente elemento
               ElementIndex = ElementIndex+1
               
               !guarda los datos conocidos del elemento sus nodos
               !--------------------------------------------
               meshT%elements(ElementIndex)%ifield = i
               meshT%elements(ElementIndex)%elementType = ElementType
               meshT%elements(ElementIndex)%Nnodes = ElementNnodes
               meshT%elements(ElementIndex)%Nedges = ElementNedges
               meshT%elements(ElementIndex)%Nfacets = ElementNfacets
               
               allocate(meshT%elements(ElementIndex)%nodes(ElementNnodes))
               allocate(meshT%elements(ElementIndex)%edges(ElementNedges))
               
               allocate(nodes_global_id(ElementNnodes))
               
               allocate(headElement(elementIndex)%edges_state(elementNedges)) !allocata edges_state e inicializa a .true. (todavía no tiene datos guardados)
               headElement(elementIndex)%edges_state = .true.
               
               do k=1,ElementNnodes !guarda los nodos del elemento por su global id en el vector nodes_global_id
                   NodeLocal_id = mesh%grid%fields(i)%connections(j,k)
                   NodeGlobal_id  = mesh%grid%fields(i)%global_id(NodeLocal_id)
                   
                   nodes_global_id(k) = NodeGlobal_id
                   
                   meshT%elements(ElementIndex)%nodes(k) = NodeGlobal_id
                   
                   meshT%nodes(NodeGlobal_id)%Nelements = meshT%nodes(NodeGlobal_id)%Nelements+1
                   
                   currentNode_elements => headNode(NodeGLobal_id)%elements !nodo: guarda el elemento
                   do m=1,headNode(NodeGlobal_id)%Nelements
                       currentNode_elements => currentNode_elements%next
                   end do
                   currentNode_elements%ID = ElementIndex
                   allocate(currentNode_elements%next)
                   headNode(NodeGlobal_id)%Nelements = headNode(NodeGlobal_id)%Nelements+1
                   
               end do
               !-------------------------------------------
               
               !crea lista enlaza con las aristas
               call EdgeCreator(meshT, headEdge, currentEdge, headElement, nodes_global_id, ElementIndex, headNode)
               
               !deallocata nodes_global_id
               deallocate(nodes_global_id)
               
           end do
       end do
       !----------------------------------------------
       
       !carga la topologia en meshT
       call meshTopologyLoad(meshT, headEdge, headNode)
       
       !deallocata headEdge
       do i=1,mesh%gridData%no_of_elements
           deallocate(headElement(i)%edges_state)
       end do
       deallocate(headElement)
       
       !deallocata headNode
       deallocate(headNode)
       
    end subroutine TopologyPreProcess3D
    
    !***************************************************************************************
    !TopologyPostProcess3D
    !
    !deallocata y finaliza meshT
    subroutine TopologyPostProcess3D(meshT)
    
       !subroutine arguments
       type(TopologyMeshClass), intent(inout) :: meshT
       !end subroutine arguments
       
       !local scalars
       integer :: i
       
       !end of header----------------------------------
       
       !deallocata elements
       do i=1,meshT%Nelements
           deallocate(meshT%elements(i)%nodes)
           deallocate(meshT%elements(i)%edges)
       end do
       deallocate(meshT%elements)
       
       !deallocata edges
       do i=1,meshT%Nedges
           deallocate(meshT%edges(i)%elements)
       end do
       deallocate(meshT%edges)
       
       !deallocata nodes
       do i=1,meshT%Nnodes
           deallocate(meshT%nodes(i)%edges)
           deallocate(meshT%nodes(i)%elements)
       end do
       deallocate(meshT%nodes)
       
    end subroutine TopologyPostProcess3D
    
end module TopologyProcess3D