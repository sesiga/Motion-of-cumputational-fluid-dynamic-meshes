!TopologyUtility
!
!utiliza listas enlazadas para crear y guardar la topologia
    
module TopologyUtility

   !modules used
   use precision, only : wp
   
   use TopologicalDataStructures
   
   use TopologicalDataStructuresSupport
   !end modules used
   
   implicit none
   
    contains
    
    !********************************************************************************************
    !EdgeCreator
    !
    !crea las aristas y las guarda en las listas enlazadas
    subroutine EdgeCreator(meshT, headEdge, currentEdge, headElement, &
        figureNodes, ElementIndex, headNode)
    
       !subroutine arguments
       type(TopologyMeshClass), intent(inout) :: meshT
       
       type(EdgeType), pointer :: headEdge !lista enlazada de aristas
       type(EdgeType), pointer :: currentEdge !lista enlazada de aristas
       
       type(NodeType), pointer :: headNode(:) !vector nodos que incluye listas enlazadas para las aristas para cada nodo
       
       type(ElementType), pointer :: headElement(:) !elementos
       
       integer, intent(in) :: figureNodes(:) !nodos del elemento, su global_id
       
       integer, intent(in) :: ElementIndex !indice del elemento en el que esta el bucle principal
       !end subroutine arguments
       
       !local scalars
       type(EdgeType_elements), pointer :: currentEdge_elements => null()
       type(EdgeType), pointer :: auxiliarEdge => null()
       
       type(NodeType_edges), pointer :: currentNode_edges => null()
       
       integer, allocatable :: nodes(:) !nodos del elemento (Nnodes+1), se añade el primer nodo al final para formar todas las aristas
       integer :: Nnodes !numero de nodos que tiene el elemento
       integer :: Nedges !numero de aristas cradas
       
       integer :: k, l, m, n
       
       !end of header-------------------------------------
       
       !crea el vector nodes con los nodos que forman todas las aristas, ordenadas de nodo en nodo
       Nnodes = size(figureNodes)
       allocate(nodes(Nnodes+1))
       nodes = reshape((/figureNodes(:), figureNodes(1)/),(/Nnodes+1/)) 
       
       !crea las aristas y comprueba si ya existen
       bucle_k: do k=1,Nnodes
           
           auxiliarEdge => headEdge
           
           bucle_l: do l=1,meshT%Nedges
               
               !comprueba si existe la arista
               !-------------------------------------------
               if(nodes(k) == auxiliarEdge%nodes(1))then
                   
                   if(nodes(k+1) == auxiliarEdge%nodes(2))then
                       
                       currentEdge_elements => auxiliarEdge%elements !bucle para llegar hasta la ultima arista alocatada
                       do m=1,auxiliarEdge%Nelements-1 
                           currentEdge_elements => currentEdge_elements%next
                       end do 
                       
                       allocate(currentEdge_elements%next) !alocata, apunta, carga los datos
                       currentEdge_elements => currentEdge_elements%next
                       currentEdge_elements%ID = ElementIndex
                       
                       auxiliarEdge%Nelements = auxiliarEdge%Nelements+1 !suma un elemento mas a la arista
                       
                       do n=1,meshT%elements(ElementIndex)%Nedges !añade arista al elemento
                           if(headElement(elementIndex)%edges_state(n))then
                              headElement(elementIndex)%edges_state(n) = .false.
                              meshT%elements(elementIndex)%edges(n) = l  
                              exit
                           end if
                       end do
                                                      
                       cycle bucle_k
                   end if
                   
               elseif(nodes(k) == auxiliarEdge%nodes(2))then
                   
                   if(nodes(k+1) == auxiliarEdge%nodes(1))then
                       
                       currentEdge_elements => auxiliarEdge%elements !bucle para llegar hasta la ultima arista alocatada
                       do m=1,auxiliarEdge%Nelements-1 
                           currentEdge_elements => currentEdge_elements%next
                       end do 
                       
                       allocate(currentEdge_elements%next) !alocata, apunta, carga los datos
                       currentEdge_elements => currentEdge_elements%next
                       currentEdge_elements%ID = ElementIndex
                       
                       auxiliarEdge%Nelements = auxiliarEdge%Nelements+1 !suma un elemento mas a la arista
                       
                       do n=1,meshT%elements(ElementIndex)%Nedges !añade arista al elemento
                           if(headElement(elementIndex)%edges_state(n))then
                              headElement(elementIndex)%edges_state(n) = .false.
                              meshT%elements(elementIndex)%edges(n) = l    
                              exit
                           end if
                       end do
                       
                       cycle bucle_k
                   end if
                   
               end if
               !----------------------------------------------
               
               auxiliarEdge => auxiliarEdge%next !siguiente arista
                              
           end do bucle_l
           
           !si la arista no existe crea una nueva
           currentEdge%Nelements = 1
           currentEdge%nodes(1) = nodes(k)
           currentEdge%nodes(2) = nodes(k+1)
           allocate(currentEdge%elements)
           currentEdge%elements%ID = elementIndex
               
           !suma una arista al numero total
           meshT%Nedges = meshT%Nedges+1
           
           !guarda los datos de la arista en los nodos
           currentNode_edges => headNode(nodes(k))%edges
           do m=1,headNode(nodes(k))%Nedges
               currentNode_edges => currentNode_edges%next
           end do
           currentNode_edges%ID = meshT%Nedges !id de la arista  
           headNode(nodes(k))%Nedges = headNode(nodes(k))%Nedges+1 !numero de aristas que tiene el nodo
           allocate(currentNode_edges%next)
           
           currentNode_edges => headNode(nodes(k+1))%edges
           do m=1,headNode(nodes(k+1))%Nedges
               currentNode_edges => currentNode_edges%next
           end do
           currentNode_edges%ID = meshT%Nedges !id de la arista  
           headNode(nodes(k+1))%Nedges = headNode(nodes(k+1))%Nedges+1 !numero de aristas que tiene el nodo
           allocate(currentNode_edges%next)
           
           !indice de la arista del elemento
           do n=1,meshT%elements(ElementIndex)%Nedges !añade arista al elemento
               if(headElement(elementIndex)%edges_state(n))then
                   headElement(elementIndex)%edges_state(n) = .false.
                   meshT%elements(elementIndex)%edges(n) = meshT%Nedges
                   exit
               end if
           end do
               
           !siguiente arista
           allocate(currentEdge%next)
           currentEdge => currentEdge%next
           
       end do bucle_k
       
       !deallocata el vector de nodos
       deallocate(Nodes)
       
    end subroutine EdgeCreator

    !************************************************************************************
    !FacetCreator
    !
    !utiliza listas enlazadas para crear las caras
    subroutine FacetCreator(meshT, headFacet, currentFacet, headEdge, currentElement, &
        figureNodes, elementIndex, no_nodes_of_facet)
    
       !subroutine arguments
       type(TopologyMeshClass), intent(inout) :: meshT
       
       type(FacetType), pointer :: headFacet !lista enlazada de caras
       type(FacetType), pointer :: currentFacet !lista enlazada de caras
       type(EdgeType),  pointer :: headEdge !lista enlazada de aristas
       type(ElementType), pointer :: currentElement !lista enlazada de elementos
       
       integer, intent(in) :: figureNodes(:) !nodos del elemento
       integer, intent(in) :: elementIndex !indice del elemento
       integer, intent(in) :: no_nodes_of_facet !numero de nodos de la cara
       !end subroutine arguments
       
       !local scalars
       type(facetType), pointer :: auxiliarFacet
       type(facetType_nodes), pointer :: currentFacet_nodes
       type(facetType_elements), pointer :: currentFacet_elements
       type(facetType_edges), pointer :: currentFacet_edges
       type(EdgeType), pointer :: currentEdge
       type(EdgeType_facets), pointer :: currentEdge_facets
       
       integer, allocatable :: nodes(:) !nodos del elemento (figureNodes(:), figureNodes(1)), añadiendo el primero al final
       integer :: sameNode !contador para ver el numero de nodos coincidentes para comprobar si el elemento esta repetido
       integer :: Nnodes !numero de nodos del elemento
       integer :: facetIndex !indice de la cara
       
       logical, allocatable :: facet_state(:) ! .true. si la cara ya esta creada
       
       integer ::  i,  j, k, m
       
       !end of header----------------------------------------------------------------------
       
       !numero de nodos del elemento
       Nnodes = size(figureNodes)
       
       !allocata face_state
       allocate(facet_state(Nnodes))
       facet_state = .false.
       
       !comprueba si la cara esta ya creada
       !--------------------------------------------------
       auxiliarFacet => headFacet
       
       bucle_i: do i=1,meshT%Nfacets
           
           sameNode = 0
           
           currentFacet_nodes => auxiliarFacet%nodes
           
           bucle_j: do j=1,auxiliarFacet%Nnodes             
               
               bucle_k: do k=1,Nnodes
                   
                   if(currentFacet_nodes%ID == figureNodes(k))then
                       sameNode = sameNode+1
                       if(sameNode == auxiliarFacet%Nnodes)then
                           currentFacet%Nelements = currentFacet%Nelements+1
                           currentFacet_elements => currentFacet_elements%next
                           currentFacet_elements%ID = elementIndex
                           facet_state(k-sameNode+1) = .true.
                           cycle bucle_i
                       end if
                   end if
                   
               end do bucle_k
               
               currentFacet_nodes => currentFacet_nodes%next
               
           end do bucle_j
           
           auxiliarFacet => auxiliarFacet%next
           
       end do bucle_i
       !------------------------------------------------------------------------------
       
       !si la cara no existe se crea
       !-------------------------------------------------------------------------------
       allocate(nodes(Nnodes+no_nodes_of_facet-1))
       
       select case(no_nodes_of_facet)
       case(3)
           nodes = reshape((/figureNodes(:), figureNodes(1), figureNodes(2)/),(/Nnodes+no_nodes_of_facet-1/))
       case(4)
           nodes = reshape((/figureNodes(:), figureNodes(1), figureNodes(2), figureNodes(3)/),(/Nnodes+no_nodes_of_facet-1/))
       end select
              
       bucle_i2: do i=1,Nnodes
           
           !si la arista esta creada siguiente iteracion
           if(facet_state(i))then
               cycle bucle_i2
           end if
           
           !suma una cara
           meshT%Nfacets = meshT%Nfacets+1
                      
           !numero de elementos, aristas y nodos
           currentFacet%Nelements = 1
           currentFacet%Nnodes = no_nodes_of_facet
           currentFacet%Nedges = no_nodes_of_facet
           
           !guarda nodos y aristas
           allocate(currentFacet%nodes)
           currentFacet_nodes => currentFacet%nodes
           
           allocate(currentFacet%edges)
           currentFacet_edges => currentFacet%edges
           
           do j=1,currentFacet%Nnodes
               
               !nodos
               !------------------------------------------
               currentFacet_nodes%ID = nodes(i+j-1)
               allocate(currentFacet_nodes%next)
               currentFacet_nodes => currentFacet_nodes%next          
               !--------------------------------------------
               
               !aristas
               !-----------------------------------------------
               currentEdge => headEdge
               
               do k=1,meshT%Nedges
                   
                   if(nodes(i+j-1) == currentEdge%nodes(1))then
                       if(nodes(i+j) == currentEdge%nodes(2))then
                           !cara
                           currentFacet_edges%ID = k
                           allocate(currentFacet_edges%next)
                           currentFacet_edges => currentFacet_edges%next
                           !arista
                           allocate(currentEdge%facets)
                           currentEdge_facets => currentEdge%facets
                           currentEdge_facets%ID = meshT%Nfacets     
                       end if
                   else if(nodes(i+j-1) == currentEdge%nodes(2))then
                       if(nodes(i+j) == currentEdge%nodes(1))then
                           !cara
                           currentFacet_edges%ID = k
                           allocate(currentFacet_edges%next)
                           currentFacet_edges => currentFacet_edges%next
                           !arista
                           allocate(currentEdge%facets)
                           currentEdge_facets => currentEdge%facets
                           currentEdge_facets%ID = meshT%Nfacets
                       end if
                   end if
                   
                   currentEdge => currentEdge%next
                   
               end do
               !----------------------------------------------
               
           end do
                      
           !guarda elemento
           allocate(currentFacet%elements)
           currentFacet%elements%ID = elementIndex
           
           !siguiente cara
           allocate(currentFacet%next)
           currentFacet => currentFacet%next
           
       end do bucle_i2
       !----------------------------------------------------------------------------------
              
       !deallocaata el vector nodes
       deallocate(nodes)
       
    end subroutine FacetCreator

    
    !************************************************************************************
    !meshTopologyLoad
    !
    !carga los datos guardados en la lista enlazada en meshT
    subroutine meshTopologyLoad(meshT, headEdge, headNode)
    
       !subroutine arguments
       type(TopologyMeshClass), intent(inout) :: meshT
       type(EdgeType), pointer :: headEdge !lista enlazada de aristas
       type(NodeType), pointer :: headNode(:)
       !end subroutine arguments
       
       !local scalars
       type(EdgeType), pointer :: currentEdge => null()
       type(EdgeType_elements), pointer :: headEdge_elements => null()
       type(EdgeType_elements), pointer :: currentEdge_elements => null()
       
       type(nodeType_edges), pointer :: headNode_edges => null()
       type(nodeType_edges), pointer :: currentNode_edges => null()
       type(nodeType_elements), pointer :: headNode_elements => null()
       type(nodeType_elements), pointer :: currentNode_elements => null()
       
       integer :: i, j
       
       !end of header--------------------------------------
       
       !carga aristas en meshT
       !----------------------------------------------------
       allocate(meshT%edges(meshT%Nedges)) 
       
       do i=1,meshT%Nedges
           
           currentEdge => headEdge
           
           !aristas: carga elementos
           headEdge_elements => currentEdge%elements
           meshT%edges(i)%Nelements = currentEdge%Nelements
           allocate(meshT%edges(i)%elements(currentEdge%Nelements))
           do j=1,currentEdge%Nelements
               currentEdge_elements => headEdge_elements
               meshT%edges(i)%elements(j) = currentEdge_elements%ID
               headEdge_elements => currentEdge_elements%next
               deallocate(currentEdge_elements)
           end do
           
           !aristas: carga nodos
           meshT%edges(i)%nodes = currentEdge%nodes
           
           headEdge => currentEdge%next
           deallocate(currentEdge)
           
       end do
       !----------------------------------------------------
       
       !carga nodos en meshT
       !----------------------------------------------------
       do i=1,meshT%Nnodes
           
           !nodo: carga aristas
           meshT%nodes(i)%Nedges = headNode(i)%Nedges
           allocate(meshT%nodes(i)%edges(headNode(i)%Nedges))
           headNode_edges => headNode(i)%edges
           do j=1,headNode(i)%Nedges
               currentNode_edges => headNode_edges
               meshT%nodes(i)%edges(j) = currentNode_edges%ID
               headNode_edges => currentNode_edges%next
               deallocate(currentNode_edges)
           end do
           
           !nodo: carga elementos
           meshT%nodes(i)%Nelements = headNode(i)%Nelements
           allocate(meshT%nodes(i)%elements(headNode(i)%Nelements))
           headNode_elements => headNode(i)%elements
           do j=1,headNode(i)%Nelements
               currentNode_elements => headNode_elements
               meshT%nodes(i)%elements(j) = currentNode_elements%ID
               headNode_elements => currentNode_elements%next
               deallocate(currentNode_elements)
           end do
           
       end do
       !----------------------------------------------------
       
    end subroutine meshTopologyLoad

end module TopologyUtility