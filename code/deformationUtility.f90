!deformationProcess
!
!modulo en el que se realiza el proceso de deformacion de la malla

module deformationUtility

!modules used
use precision, only : wp

use TopologicalDataStructures

use RadialBasisFunctions

use RBFinterpolation
!end modules used

implicit none

    type nodes_id
        integer :: id
        type(nodes_id), pointer :: next => null()
    end type nodes_id
              
    type nodes_def
       integer :: Nnodes
       type(nodes_id), pointer :: nodes => null()
    end type nodes_def

    contains
    
    !end of header----------------------------------
    
    !*******************************************************************************
    !deformationProcess
    !
    !subroutina principal en la que se realiza el proceso de deformacion de malla
    
    subroutine deformationProcess(meshT, PointsIn, PointsOut, ValuesIn, ValuesOut, &
        PointsIn_id, PointsOut_id, RBFfunction, radiusCoeff)
    
       !subroutine arguments
       type(TopologyMeshClass), intent(in) :: meshT
       real(wp), intent(in) :: PointsIn(:,:)
       real(wp), intent(in) :: ValuesIn(:,:)
       real(wp), intent(in) :: PointsOut(:,:)
       real(wp), intent(inout) :: ValuesOut(:,:)
       integer, intent(in) :: PointsIn_id(:)
       integer, intent(in) :: PointsOut_id(:)
       procedure(RBF), pointer :: RBFfunction
       real(wp), intent(in) :: radiusCoeff
       !end subroutine arguments
       
       !local scalars
       logical :: nodes_state(meshT%Nnodes) !estado del nodo: .false. si no se ha deformado
       logical :: child_state(meshT%Nnodes) !estado del nodo: .false. si no ha sido child 
       integer :: nodes_layer(meshT%Nnodes) !capa de deformacion a la que pertenece el nodo
       
       type(nodes_def) :: headnodes !global id de los nodos que se utilizan en cada iteracion para deformar los siguientes nodos
       type(nodes_id), pointer :: headNodes_support => null()
       
       type(nodes_def) :: childnodes !global_id de los nodos deformados en cada iteracion
       type(nodes_id), pointer :: childnodes_support => null()
       
       integer :: layer
       integer :: main_node !nodo modificador
       
       type(nodes_def) :: modified_nodes !global id de los nodos a modificar por el nodo modificador
       type(nodes_id), pointer :: modified_nodes_support => null()
       
       type(nodes_def) :: main_nodes !nodos modificadores
       
       real(wp) :: pOut(1,2) !pointOut: coordenadas del nodo a deformar
       real(wp) :: vOut(1,2) !valueOut: deformacion del nodo a deformar
       real(wp), allocatable :: pIn(:,:) !pointsIn: coordenadas de los nodos deformadores
       real(wp), allocatable :: vIn(:,:) !valuesIn : deformacion de los nodos deformadores
       
       real(wp) :: coeficiente_hardy !coeficiente de hardy para deformacion de la capa limite
       real(wp) :: rho 
       
       integer :: i, j
       
       !end of header-------------------------------
       
       !inicializa las variables
       layer = 0
       nodes_state = .false.
       nodes_layer = 0
       child_state = .false.
       
       !crea lista enlaza con los nodos deformados del fichero scatter
       call linkedList_nodes_preProcess(pointsIn_id, headnodes)
       
       !calcula el coeficiente de hardy
       coeficiente_hardy = coef_hardy(PointsIn)
       
       !proceso de deformacion
       do
           
           !suma una capa
           layer = layer + 1
           
           !establece el coeficiente rho
           if(layer < 40)then
               rho = coeficiente_hardy
               RBFfunction => hardy
           else 
               rho = radiusCoeff
               RBFfunction => Spline
           end if
           
           !guarda la capa en la que estan los nodos deformados y cambia el estado
           call layer_of_nodes(headnodes, nodes_layer, nodes_state, layer)
           
           !para cada nodo deformador, busca a nodos que deforma y lo deforma
           headnodes_support => headnodes%nodes
           
           !inicializa childnodes
           childnodes%nnodes = 0
           allocate(childnodes%nodes)
           
           !recorre todos los nodos deformadores
           do i=1,headnodes%nnodes
               
               !global id del nodo modificador
               main_node = headnodes_support%ID
               
               !selecciona los nodos a modificar, guarda esos nodos en childnodes
               call nodes_to_modify_preProcess(main_node, modified_nodes, meshT, nodes_layer, nodes_state, childnodes, child_state)
               
               !apunta a los nodos a deformar
               modified_nodes_support => modified_nodes%nodes
               
               !deforma los nodos
               do j=1,modified_nodes%nnodes
                   
                   !para cada nodo a deformar, selecciona los nodos deformadores                   
                   call modifier_nodes_preProcess(modified_nodes_support%ID, main_node, main_nodes, meshT, nodes_layer)
                   
                   !genera los arrays para el proceso de deformacion
                   call generate_deformation_arrays(pointsIn, pointsOut, valuesIn, valuesOut, &
                   pointsIn_id, pointsOut_id, main_nodes, modified_nodes_support%ID, pOut, pIn, vIn)
                   
                   !interpolacion
                   call RBFinterpolator(RBFfunction, pIn, pOut, vIn, vOut, rho)
                   
                   !guarda vOut en valuesOut
                   call save_deformation(modified_nodes_support%ID, pointsOut_id, vOut, valuesOut)
                   
                   !deallocata la lista enlazada main_nodes
                   call modifier_nodes_postProcess(main_nodes)
                   
                   !siguiente nodo a deformar
                   modified_nodes_support => modified_nodes_support%next
                   
                   !deallocata vectores de nodos deformadores
                   deallocate(pIn, vIn)
                   
               end do
               
               !deallocata la lista enlazada con los nodos a modificar
               call nodes_to_modify_postProcess(modified_nodes)
               
               !apunta al siguiente nodo deformador
               headnodes_support => headnodes_support%next
               
           end do
           
           !pasa la lista enlazada de childnodes a headnodes
           call head_child_process(headnodes, childnodes)
           
           !comprueba si quedan nodos por deformar, si ya estan todos deformados sale del bucle
           if(headnodes%nnodes == 0) exit
           
       end do
              
       !deallocata la lista enlazada headnodes
       call linkedList_nodes_postProcess(headnodes)
       
    end subroutine deformationProcess
    
    !*******************************************************************
    !head_child_process
    !
    !pasa la lista enlazada de childnodes a headnodes, para utilizar los puntos
    !deformados como deformadores en la siguiente iteracion
    subroutine head_child_process(headnodes, childnodes)
    
       !subroutine arguments
       type(nodes_def) :: headnodes
       type(nodes_def) :: childnodes
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: headnodes_head => null()
       type(nodes_id), pointer :: headnodes_current => null()
       type(nodes_id), pointer :: childnodes_head => null()
       type(nodes_id), pointer :: childnodes_current => null()
       
       integer :: i
       
       !end of header------------------------------
       
       !deallocata los nodos guardados en headnodes
       headnodes_head => headnodes%nodes
       do i=1,headnodes%Nnodes
           headnodes_current => headnodes_head
           headnodes_head => headnodes_current%next
           deallocate(headnodes_current)
       end do
       
       !guarda childnodes en headnodes y deallocata childnodes
       headnodes%Nnodes = childnodes%Nnodes
       allocate(headnodes%nodes)
       headnodes_head => headnodes%nodes
       childnodes_head => childnodes%nodes
       do i=1,childnodes%Nnodes
           childnodes_current => childnodes_head
           headnodes_head%id = childnodes_current%id
           allocate(headnodes_head%next)
           headnodes_head => headnodes_head%next
           childnodes_head => childnodes_current%next
           deallocate(childnodes_current)
       end do
           
    end subroutine head_child_process
    
    !*************************************************************************
    !linkedlist_nodes_preprocess
    !
    !crea lista enlaza con los nodos deformados del fichero scatter
    subroutine linkedList_nodes_preProcess(pointsIn_id, headnodes)
    
       !subroutine arguments
       integer, intent(in) :: pointsIn_id(:)
       type(nodes_def) :: headnodes !global id de los nodos que se utilizan en cada iteracion para deformar los siguientes nodos
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: headNodes_support => null()
       
       integer :: i
       
       !end of header-----------------------------
       
       !crea lista enlaza con los nodos deformados del fichero scatter
       headnodes%Nnodes = 0
       allocate(headnodes%nodes)
       headnodes_support => headnodes%nodes
       do i=1,size(pointsIn_id)
           headnodes_support%id = pointsIn_id(i)
           headnodes%Nnodes = headnodes%Nnodes + 1
           allocate(headnodes_support%next)
           headnodes_support => headnodes_support%next
       end do     
       
    end subroutine linkedList_nodes_preProcess
    
    !*************************************************************************
    !linkedList_nodes_postProcess
    !
    !deallocata las listas enlazadas
    subroutine linkedList_nodes_postProcess(headnodes)
    
       !subroutine arguments
       type(nodes_def) :: headnodes
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: headNodes_head => null()
       type(nodes_id), pointer :: headNodes_current => null()
       
       integer :: i
       
       !end of header---------------------------------
       
       !deallocata la lista enlazada
       headNodes_head => headnodes%nodes
       do i=1,headnodes%Nnodes
           headnodes_current => headnodes_head
           headnodes_head => headnodes_current%next
           deallocate(headnodes_current)
       end do
       
    end subroutine linkedList_nodes_postProcess
    
    !************************************************************************
    !layer_of_nodes
    !
    !gaurda la capa en la que esta cada nodo
    subroutine layer_of_nodes(headnodes, nodes_layer, nodes_state, layer)
    
       !subroutine arguments
       type(nodes_def), intent(inout) :: headnodes
       integer, intent(inout) :: nodes_layer(:)
       logical, intent(inout) :: nodes_state(:)
       integer, intent(in) :: layer
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: headNodes_support => null()
       
       integer :: i
       
       !end of header-------------------------------
       
       headnodes_support => headnodes%nodes
       
       do i=1,headnodes%nnodes
           nodes_layer(headnodes_support%ID) =  layer
           nodes_state(headnodes_support%ID) = .true.
           
           headnodes_support => headnodes_support%next
       end do
       
    end subroutine layer_of_nodes
    
    !*************************************************************************
    !nodes_to_modify_preProcess
    !
    !dado el nodo deformador, deforma a los nodos correpondientes
    subroutine nodes_to_modify_preProcess(main_node, modified_nodes, meshT, nodes_layer, nodes_state, childnodes, child_state)
    
       !subroutine arguments
       integer, intent(in) :: main_node
       type(nodes_def), intent(out) :: modified_nodes
       type(topologymeshclass), intent(in) :: meshT
       integer, intent(in) :: nodes_layer(:)
       logical, intent(in) :: nodes_state(:)
       type(nodes_def), intent(inout) :: childnodes
       logical, intent(inout) :: child_state(:)
       !end subroutine arguments
       
       !local scalars
       integer :: edge_id !arista utilizada en cada iteracion del bucle
       integer :: node_id !nodo de la arista
       
       type(nodes_id), pointer :: modified_nodes_support => null()
       type(nodes_id), pointer :: childnodes_support => null()
       
       integer :: i, j
       
       !end of header------------------------------------------
       
       !inicializa la lista enlazada de nodos a modificar
       modified_nodes%nnodes = 0
       allocate(modified_nodes%nodes)
       modified_nodes_support => modified_nodes%nodes
       
       !recorre los childnodes hasta llegar al final de la lista, para guardar los siguientes nodos
       childnodes_support => childnodes%nodes
       do i=1,childnodes%nnodes
           childnodes_support => childnodes_support%next
       end do
       
       !para cada arista del nodo deformador, encuentra los nodos a deformar
       !guarda los nodos a deformar en la lista enlaza modified_nodes
       bucle_i: do i=1,meshT%nodes(main_node)%Nedges
           
           !arista y nodo a deformar de la arista
           edge_id = meshT%nodes(main_node)%edges(i)
           
           if(meshT%edges(edge_id)%nodes(1) == main_node)then
               node_id = meshT%edges(edge_id)%nodes(2)
           else 
               node_id = meshT%edges(edge_id)%nodes(1)
           end if
           
           !si el nodo ya esta deformado se pasa a la siguiente arista
           if(nodes_state(node_id)) cycle bucle_i
           if(child_state(node_id)) cycle bucle_i

           !si el nodo no esta deformado se guarda en las listas enlazadas
           modified_nodes%nnodes = modified_nodes%nnodes + 1
           modified_nodes_support%ID = node_id
           allocate(modified_nodes_support%next)
           modified_nodes_support => modified_nodes_support%next
           
           childnodes%nnodes = childnodes%nnodes + 1
           childnodes_support%ID = node_id
           child_state(node_id) = .true.
           allocate(childnodes_support%next)
           childnodes_support => childnodes_support%next
           
       end do bucle_i
              
    end subroutine nodes_to_modify_preProcess
    
    !*************************************************************************
    !nodes_to_modify_postProcess
    !
    !deallocata las listas enlazadas
    subroutine nodes_to_modify_postProcess(modified_nodes)
    
       !subroutine arguments
       type(nodes_def) :: modified_nodes
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: modified_nodes_head => null()
       type(nodes_id), pointer :: modified_nodes_current => null()
       
       integer ::  i
       
       !end of header----------------------------------
       
       !dealocata la lista enlazada
       modified_nodes_head => modified_nodes%nodes
       
       do i=1,modified_nodes%nnodes
           modified_nodes_current => modified_nodes_head
           modified_nodes_head => modified_nodes_current%next
           deallocate(modified_nodes_current)
       end do
    
    end subroutine nodes_to_modify_postProcess
    
    !***************************************************************************
    !modifier_nodes_preProcess
    !
    !selecciona los nodos modificadores
    subroutine modifier_nodes_preProcess(modified_node, main_node, main_nodes, meshT, nodes_layer)
    
       !subroutine arguments
       integer, intent(in) :: modified_node
       integer, intent(in) :: main_node
       type(nodes_def), intent(out) :: main_nodes
       type(TopologyMeshClass), intent(in) :: meshT
       integer, intent(in) :: nodes_layer(:)
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: main_nodes_support => null()
       
       integer :: element_id
       integer :: node_id
       integer :: edge_id
       
       integer :: i, j, k
       
       !end of header---------------------------
       
       !se busca los elementos en comun de nodo deformador y el nodo a deformar
       !en esos elementos se buscan los demas nodos deformadores
       !tienen que estar en la misma capa
       
       !guarda los nodos
       !allocate(main_nodes%nodes)
       !main_nodes%nnodes = 1
       !main_nodes%nodes%ID = main_node
       !allocate(main_nodes%nodes%next)
       !main_nodes_support => main_nodes%nodes%next
       
       !do i=1,meshT%nodes(main_node)%nedges
           
       !    edge_id = meshT%nodes(main_node)%edges(i)
           
       !    if(meshT%edges(edge_id)%nodes(1) == main_node)then
       !        node_id = meshT%edges(edge_id)%nodes(2)
       !    else
       !        node_id = meshT%edges(edge_id)%nodes(1)
       !    end if
           
       !    if(nodes_layer(node_id) == nodes_layer(main_node))then
               
       !        main_nodes%nnodes = main_nodes%nnodes + 1
       !        main_nodes_support%ID = node_id
       !        allocate(main_nodes_support%next)
       !        main_nodes_support => main_nodes_support%next
       !    end if
           
       !end do
       
       !guarda el nodo principal
       allocate(main_nodes%nodes)
       main_nodes%nodes%ID = main_node
       main_nodes%nnodes = 1
       
       !alocata siguiente nodo
       allocate(main_nodes%nodes%next)
       main_nodes_support => main_nodes%nodes%next
       
       !busca los demas nodos
       do i=1,meshT%nodes(main_node)%nelements
           do j=1,meshT%nodes(modified_node)%nelements
           
               if(meshT%nodes(main_node)%elements(i) == meshT%nodes(modified_node)%elements(j))then
                   
                   !elemento comun
                   element_id = meshT%nodes(main_node)%elements(i)
                   
                   !busca los demas nodos
                   do k=1,meshT%elements(element_id)%nnodes
                       
                       node_id = meshT%elements(element_id)%nodes(k)
                   
                       if(nodes_layer(node_id) == nodes_layer(main_node))then
                           if (node_id /= main_node)then
                               
                               main_nodes%nnodes = main_nodes%nnodes + 1
                               main_nodes_support%ID = node_id
                               allocate(main_nodes_support%next)
                               main_nodes_support => main_nodes_support%next
                               
                           end if
                       end if
                           
                   end do
                   
                   !sale del bucle
                   exit
                   
               end if
               
           end do
       end do
       
    end subroutine modifier_nodes_preProcess
    
    !***************************************************************************
    !modifier_nodes_postProcess
    !
    !deallocata main_nodes
    subroutine modifier_nodes_postProcess(main_nodes)
    
       !subroutine arguments
       type(nodes_def), intent(out) :: main_nodes
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: main_nodes_head => null()
       type(nodes_id), pointer :: main_nodes_current => null()
       
       integer :: i
       
       !end of header-------------------------------------------
       
       main_nodes_head => main_nodes%nodes
       
       do i=1,main_nodes%nnodes
           main_nodes_current => main_nodes_head
           main_nodes_head => main_nodes_current%next
           deallocate(main_nodes_current)
       end do
       
       
    end subroutine modifier_nodes_postProcess
    
    !***************************************************************************
    !generate_deformation_arrays
    !
    !genera los arrays para el proceso de deformacion
    subroutine generate_deformation_arrays(pointsIn, pointsOut, valuesIn, valuesOut, &
        pointsIn_id, pointsOut_id, main_nodes, modified_node, pOut, pIn, vIn)
    
       !subroutine arguments
       real(wp), intent(in) :: PointsIn(:,:)
       real(wp), intent(in) :: ValuesIn(:,:)
       real(wp), intent(in) :: PointsOut(:,:)
       real(wp), intent(in) :: ValuesOut(:,:)
       integer, intent(in) :: PointsIn_id(:)
       integer, intent(in) :: PointsOut_id(:)
       type(nodes_def), intent(in) :: main_nodes
       integer, intent(in) :: modified_node
       real(wp), allocatable, intent(out) :: pIn(:,:)
       real(wp), allocatable, intent(out) :: vIn(:,:)
       real(wp), intent(out) :: pOut(:,:)
       !end subroutine arguments
       
       !local scalars
       type(nodes_id), pointer :: main_nodes_current => null()
       
       integer :: i, j, k
       
       !end of header-----------------------------------
       
       !nodo a deformar
       do i=1,size(pointsOut_id)
           
           if(modified_node == pointsOut_id(i))then
               pOut(1,:) = pointsOut(i,:)
               exit
           end if
           
       end do
       
       !nodos deformadores
       allocate(pIn(main_nodes%nnodes,2))
       allocate(vIn(main_nodes%nnodes,2))
       
       main_nodes_current => main_nodes%nodes
       
       k = 0
       bucle_1: do i=1,main_nodes%nnodes
           
           !comprueba si es nodo guardado en pointsOut y lo guarda
           do j=1,size(pointsOut_id)
               if(main_nodes_current%id == pointsOut_id(j))then
                   k = k+1
                   pIn(k,:) = pointsOut(j,:)
                   vIn(k,:) = valuesOut(j,:)
                   main_nodes_current => main_nodes_current%next
                   cycle bucle_1
               end if
           end do
           
           !comprueba si es nodo guardado en pointsIn y lo guarda
           do j=1,size(pointsIn_id)
               if(main_nodes_current%id == pointsIn_id(j))then
                   k = k+1
                   pIn(k,:) = pointsIn(j,:)
                   vIn(k,:) = valuesIn(j,:)
                   main_nodes_current => main_nodes_current%next
                   cycle bucle_1
               end if
           end do
           
       end do bucle_1
       
        end subroutine generate_deformation_arrays
    
    !********************************************************************************
    !save_deformation
    !
    !guarda la deformacion del nodo en valuesOut
    subroutine save_deformation(modified_node, pointsOut_id, vOut, valuesOut)
    
       !subroutine arguments
       integer, intent(in) :: modified_node
       integer, intent(in) :: pointsOut_id(:)
       real(wp), intent(in) :: vOut(:,:)
       real(wp), intent(out) :: valuesOut(:,:)
       !end subroutine arguments
       
       !local scalars
       integer :: i
       
       !end of header------------------------------------
       
       !guarda la deformacion
       do i=1,size(pointsOut_id)
           
           if(modified_node == pointsOut_id(i))then
               valuesOut(i,:) = vOut(1,:)
               exit
           end if
           
       end do       
       
    end subroutine save_deformation
    
    !***************************************************************************
    !coef_hardy
    !
    !funcion que calcula el coeficiente de hardy: longitud de la arista mas corta del perfil
    function coef_hardy(PointsIn) result(a)
    
       !function arguments
       real(wp), intent(in) :: PointsIn(:,:)
       real(wp) :: a
       !end function arguments
       
       !local scalars
       real(wp) :: d
       
       integer :: i, j, k
       
       !end of header-----------------------------------
       
       !inicializa a 
       a = 1000_wp
       
       !tamaño del vector de puntos deformados de entrada
       k = size(pointsIn, dim=1)
       
       !calcula el coeficiente
       do i=1,k
           do j=i+1,k
               d = norm2(PointsIn(i,:)-PointsIn(j,:))
               if(d < a) a = d
           end do
       end do
       
    end function coef_hardy
    
end module deformationUtility