!TopologyProcess
!
!carga la estructura de la malla usando listas enlazadas
    
module TopologyProcess

   !modules usados
   use TopologicalDataStructure2D
   
   use TecplotFileDataClass
   
   use TecplotGridDataClass
   
   use precision, only : wp
   !end module used
   
   implicit none
   
    contains  
    !end of header-----------------------
    
    
    !*****************************************************************************************
    !TopologyPreProcess
    !
    !utiliza listas enlazadas para cargar la topologia de la malla en un tipo derivado
    subroutine TopologyPreProcess(mesh, meshT)
    
       !subroutine arguments
       type(TecplotFileClass), intent(in) :: mesh !mesh
       type(TopologyType), intent(out) :: meshT !mesh topology
       !end subroutine arguments
       
       !local scalars
       type(ElementType), pointer :: headElement
       type(ElementType), pointer :: currentElement
       
       type(NodeType), pointer :: headNode
       type(NodeType), pointer :: currentNode
       type(NodeType), pointer :: auxiliarNode
       
       type(EdgeType), pointer :: headEdge
       type(EdgeType), pointer :: currentEdge
       type(EdgeType), pointer :: auxiliarEdge
       
       type(NodeTypeElements), pointer :: headNodeElement
       type(NodeTypeElements), pointer :: currentNodeElement
       
       type(EdgeTypeSupport), pointer :: headEdgeSupport
       type(EdgeTypeSupport), pointer :: currentEdgeSupport
       
       integer :: noEdges !numero de aristas
       integer :: noNodes !numero de nodos
       integer :: noElements, noElementsAux !numero de elementos
       
       integer :: edge1, edge2, edge3, edge4
       
       integer :: i, j, k, m, n, p !integer auxiliares 
       
       !end of header------------------------------------------

       !crea lista enlazada para guardar topologia
       allocate(headElement)
       currentElement => headElement
       
       allocate(headNode)
       currentNode => headNode
       
       allocate(headEdge)
       currentEdge => headEdge
       
       allocate(headNodeElement)
       currentNodeElement => headNodeElement
       
       allocate(headEdgeSupport)
       currentEdgeSupport => headEdgeSupport
       
       !calculos
       noNodes = 1
       noEdges = 1
       noElements = 1
       
       print '(a)', "inicia bucle 1"
       
       do i=1,size(mesh%grid%fields)
           do j=1,size(mesh%grid%fields(i)%connections, dim=1)
               
               !Elemento: carga tipo de elemento, numero de nodos y numero de aristas
               currentElement%figure = mesh%grid%fields(i)%elementType
               currentElement%ID = noElements
               currentElement%Nnodes = size(mesh%grid%fields(i)%connections, dim=2)
               currentElement%Nedges = size(mesh%grid%fields(i)%connections, dim=2)
               
               !Elemento: alocata nodos y aristas
               allocate(currentElement%nodes(currentElement%Nnodes))
            
               !Elemento: carga nodos
               currentElement%nodes(:) = mesh%grid%fields(i)%connections(j,:)
               
               !Nodo
               bucle1: do k=1,size(mesh%grid%fields(i)%connections, dim=2)
                   
                   !comprueba que el nodo no este repetido
                   auxiliarNode => headNode
                   bucle2: do m=1,noNodes-1
                       
                       if(auxiliarNode%ID == mesh%grid%fields(i)%connections(j,k))then
                           cycle bucle1
                       end if
                       
                       auxiliarNode => auxiliarNode%next
                       
                   end do bucle2
                   
                   !Nodo: carga el numero de elementos por nodo
                   currentNode%Nelements = 1
                   currentNode%Nedges = 1
                   
                   !lista enlazada auxiliar para saber que elementos estan en cada nodo
                   currentNodeElement%element = noElements
                   allocate(currentNodeElement%next)
                   currentNodeElement => currentNodeElement%next
                   
                   noElementsAux = noElements
                   
                   do m=i,size(mesh%grid%fields)
                       do n=j+1,size(mesh%grid%fields(m)%connections, dim=1)
                           do p=1,size(mesh%grid%fields(m)%connections, dim=2)
                           
                               if(mesh%grid%fields(i)%connections(j,k) == mesh%grid%fields(m)%connections(n,p))then
                                   
                                   currentNode%Nelements = currentNode%Nelements+1
                                   currentNode%Nedges = currentNode%Nedges+1
                                   
                                   currentNodeElement%element = noElementsAux
                                   allocate(currentNodeElement%next)
                                   currentNodeElement => currentNodeElement%next
                                   
                                   noElementsAux = noElementsAux+1
                                   
                               end if
                           
                           end do                           
                       end do
                   end do
                   
                   !Nodo: carga el global_id
                   currentNode%ID = mesh%grid%fields(i)%connections(j,k)
                   
                   !Nodo: allocata siguiente eslabon
                   allocate(currentNode%next)
                   currentNode => currentNode%next
                   
                   !suma 1 unidad al numero de nodos
                   noNodes = noNodes+1
                   
               end do bucle1
               
               !Edge                     
               bucle3: do k=1,size(mesh%grid%fields(i)%connections, dim=2)-1
                   
                   !comprueba que no se repitan aristas
                   auxiliarEdge => headEdge
                   
                   edge1 = mesh%grid%fields(i)%connections(j,k)
                   edge2 = mesh%grid%fields(i)%connections(j,k+1)
                   
                   bucle4: do m=1,noEdges-1
                       
                       edge3 = auxiliarEdge%nodes(1)
                       edge4 = auxiliarEdge%nodes(2)
                       
                       if((edge1 == edge3 .or. edge1 == edge4) .and. (edge2 == edge3 .or. edge2 == edge4))then
                           cycle bucle3
                       end if
                       
                       auxiliarEdge => auxiliarEdge%next
                       
                   end do bucle4
                   
                   !carga los datos de la arista
                   currentEdge%Nnodes = 2
                   currentEdge%ID = noEdges
                   currentEdge%nodes(1) = mesh%grid%fields(i)%connections(j,k)
                   currentEdge%nodes(2) = mesh%grid%fields(i)%connections(j,k+1)
                   
                   !cuenta numero de elementos a los que pertenece cada arista
                   currentEdge%Nelements = 1
                   
                   currentEdgeSupport%element = noElements
                   allocate(currentEdgeSupport%next)
                   currentEdgeSupport => currentEdgeSupport%next
                   
                   noElementsAux = noElements
                   
                   do m=i,size(mesh%grid%fields)
                       do n=j+1,size(mesh%grid%fields(m)%connections, dim=1)
                           
                           noElementsAux = noElementsAux+1
                           
                           do p=1,size(mesh%grid%fields(m)%connections, dim=2)-1
                               
                               edge3 = mesh%grid%fields(m)%connections(n,p)
                               edge4 = mesh%grid%fields(m)%connections(n,p+1)
                           
                               if((edge1 == edge3 .or. edge1 == edge4) .and. (edge2 == edge3 .or. edge2 == edge4))then
                                   
                                   currentEdgeSupport%element = noElementsAux
                                   allocate(currentEdgeSupport%next)
                                   currentEdgeSupport => currentEdgeSupport%next
                                   
                                   currentEdge%Nelements = currentEdge%Nelements+1
                                   
                               end if
                                                          
                           end do
                           
                           edge3 = mesh%grid%fields(m)%connections(n,size(mesh%grid%fields(m)%connections, dim=2))
                           edge4 = mesh%grid%fields(m)%connections(n,1)
                           
                           noElementsAux = noElementsAux+1
                           
                           if((edge1 == edge3 .or. edge1 == edge4) .and. (edge2 == edge3 .or. edge2 == edge4))then
                               
                               currentEdgeSupport%element = noElementsAux
                               allocate(currentEdgeSupport%next)
                               currentEdgeSupport => currentEdgeSupport%next
                               
                               currentEdge%Nelements = currentEdge%Nelements+1
                               
                           end if
                           
                       end do
                   end do
                   
                   !cuenta numero de aristas, asigna ID a la arista
                   noEdges = noEdges+1                   
                   
                   !allocata siguiente arista
                   allocate(currentEdge%next)
                   currentEdge => currentEdge%Next
                                      
               end do bucle3
               
               !misma operacion para el caso de la ultima arista del elemento
               !formada por el ultimo y el primer nodo
               edge1 = mesh%grid%fields(i)%connections(j,size(mesh%grid%fields(i)%connections, dim=2))
               edge2 = mesh%grid%fields(i)%connections(j,1)
               
               !comprueba que la arista no se repita
               auxiliarEdge => headEdge
               
               n = 0
               bucle5: do m=1,noEdges-1
                   
                   edge3 = auxiliarEdge%nodes(1)
                   edge4 = auxiliarEdge%nodes(2)
                   
                   if((edge1 == edge3 .or. edge1 == edge4) .and. (edge2 == edge3 .or. edge2 == edge4))then
                       n = 1
                       exit bucle5
                   end if
                   
                   auxiliarEdge => auxiliarEdge%next
                   
               end do bucle5
               
               !si no se repite carga los datos
               if(n == 0)then
                   
                   !carga los datos de la arista
                   currentEdge%Nnodes = 2
                   currentEdge%ID = noEdges
                   currentEdge%nodes(1) = mesh%grid%fields(i)%connections(j,size(mesh%grid%fields(i)%connections, dim=2))
                   currentEdge%nodes(2) = mesh%grid%fields(i)%connections(j,1)
                   
                   !cuenta numero de elementos a los que pertenece cada arista
                   currentEdge%Nelements = 1
                   
                   noElementsAux = noElements
                   
                   do m=i,size(mesh%grid%fields)
                       do n=j+1,size(mesh%grid%fields(m)%connections, dim=1)
                           
                           noElementsAux = noElementsAux+1
                           
                           do p=1,size(mesh%grid%fields(m)%connections, dim=2)-1
                               
                               edge3 = mesh%grid%fields(m)%connections(n,p)
                               edge4 = mesh%grid%fields(m)%connections(n,p+1)
                           
                               if((edge1 == edge3 .or. edge1 == edge4) .and. (edge2 == edge3 .or. edge2 == edge4))then
                                   
                                   currentEdgeSupport%element = noElementsAux
                                   allocate(currentEdgeSupport%next)
                                   currentEdgeSupport => currentEdgeSupport%next
                                   
                                   currentEdge%Nelements = currentEdge%Nelements+1
                                   
                               end if
                                                          
                           end do
                           
                           edge3 = mesh%grid%fields(m)%connections(n,size(mesh%grid%fields(m)%connections, dim=2))
                           edge4 = mesh%grid%fields(m)%connections(n,1)
                           
                           noElementsAux = noElementsAux+1
                           
                           if((edge1 == edge3 .or. edge1 == edge4) .and. (edge2 == edge3 .or. edge2 == edge4))then
                               
                               currentEdgeSupport%element = noElementsAux
                               allocate(currentEdgeSupport%next)
                               currentEdgeSupport => currentEdgeSupport%next
                               
                               currentEdge%Nelements = currentEdge%Nelements+1
                               
                           end if
                           
                       end do
                   end do
                   
                   !cuenta numero de aristas, asigna ID a la arista
                   noEdges = noEdges+1                   
                   
                   !allocata siguiente arista
                   allocate(currentEdge%next)
                   currentEdge => currentEdge%Next
                   
               end if
               
               !suma una unidad al numero de elementos
               noElements = noElements+1
               
               !allocata siguiente elemento
               allocate(currentElement%next)
               currentElement => currentElement%next
               
           end do
       end do
       
       print '(a)', "bucle 1 terminado"
              
       !Nodo: carga elementos a los que pertenece cada nodo
       !deallocata currentNodeElement
       currentNode => headNode
       
       do i=1,noNodes
           
           allocate(currentNode%elements(currentNode%Nelements))
           
           do j=1,currentNode%Nelements
               
               currentNodeElement => headNodeElement
               
               currentNode%elements(j) = currentNodeElement%element
               
               headNodeElement => currentNodeElement%next
               deallocate(currentNodeElement)
               
           end do
           
           currentNode => currentNode%next
           
       end do
       
       !edge: carga elementos a los que pertenece cara arista
       !deallocata currentEdgeSupport
       currentEdge => headEdge
       
       do i=1,noElements
           
           allocate(currentEdge%elements(currentEdge%Nelements))
        
           do j=1,currentEdge%Nelements
               
               currentEdgeSupport => headEdgeSupport
               
               currentEdge%elements(j) = currentEdgeSupport%element
               
               headEdgeSupport => currentEdgeSupport%next
               deallocate(currentEdgeSupport)
               
           end do
           
           currentEdge => currentEdge%next
           
       end do
       
       !Element: carga aristas
       currentElement => headElement
       
       do i=1,noElements
           
           m = 0
           
           allocate(currentElement%edges(currentElement%Nedges))
           
           currentEdge => headEdge
           
           
           do j=1,noEdges
               do k=1,size(currentEdge%elements)
                   
                   if(currentElement%ID == currentEdge%Elements(k))then
                       m = m+1
                       !currentElement%edges(m) = currentEdge%ID
                   end if
                   
               end do
               
               currentEdge => currentEdge%next
               
           end do
           
           currentElement => currentElement%next
           
       end do
       
     
       
       !carga los datos en el typo meshT-------------------------
       !deallocata las listas enlazadas
       
       !element      
       meshT%Nelements = noElements
       allocate(meshT%elements(noElements))
       
       do i=1,noElements
           
           currentElement => headElement
           
           meshT%elements(i)%figure = currentElement%figure
           meshT%elements(i)%ID = currentElement%ID
           meshT%elements(i)%Nnodes = currentElement%Nnodes
           meshT%elements(i)%Nedges = currentElement%Nedges
           
           allocate(meshT%elements(i)%nodes(currentElement%Nnodes))
           !allocate(meshT%elements(i)%edges(currentElement%Nedges))
           
           meshT%elements(i)%nodes = currentElement%nodes
           !meshT%elements(i)%edges = currentElement%edges
           
           
           headElement => currentElement%next
           
       end do
       
       !node
       meshT%Nnodes = noNodes
       allocate(meshT%nodes(noNodes))
       
       do i=1,noNodes
           
           currentNode => headNode
           
           meshT%nodes(i)%ID = currentNode%ID
           meshT%nodes(i)%Nelements = currentNode%Nelements
           meshT%nodes(i)%Nedges = currentNode%Nedges
           
           allocate(meshT%nodes(i)%elements(currentNode%Nelements))
           !allocate(meshT%nodes(i)%edges(currentNode%Nedges))
           
           meshT%nodes(i)%elements = currentNode%elements
           !meshT%nodes(i)%edges = currentNode%edges
           
           
           headNode => currentNode%next
           
       end do
       
       !edge
       meshT%Nedges = noEdges
       allocate(meshT%edges(noEdges))
       
       do i=1,noEdges
           
           currentEdge => headEdge
           
           meshT%edges(i)%ID = currentEdge%ID
           meshT%edges(i)%Nelements = currentEdge%Nelements
           meshT%edges(i)%Nnodes = currentEdge%Nnodes
           
           allocate(meshT%edges(i)%elements(currentEdge%Nelements))
           
           !meshT%edges(i)%elements = currentEdge%elements
           !meshT%edges(i)%nodes = currentEdge%nodes
           
           
           headEdge => currentEdge%next
           
       end do
       
    end subroutine TopologyPreProcess
    
    !**********************************************************************************
    !TopologyPostProcess
    !
    !deallocata el tipo creado para guardar los datos de la malla
    subroutine TopologyPostProcess(meshT)
    
       !subroutine arguments
       type(TopologyType), intent(inout) :: meshT
       !end subroutine arguments
       
       !local scalars
       integer :: i
       
       !end of header------------------------------------------
       
       !deallocate nodes
       do i=1,meshT%Nnodes
           deallocate(meshT%nodes(i)%elements)
       !    deallocate(meshT%nodes(i)%edges)
       end do
       
       deallocate(meshT%nodes)
       
       !deallocate elements
       do i=1,meshT%Nelements
           deallocate(meshT%elements(i)%nodes)
           !deallocate(meshT%elements(i)%edges)
       end do
       
       deallocate(meshT%elements)
       
       !deallocate edges
       do i=1,meshT%Nedges
           deallocate(meshT%edges(i)%elements)
       end do
       
       deallocate(meshT%edges)      
       
    end subroutine TopologyPostProcess
   
end module TopologyProcess