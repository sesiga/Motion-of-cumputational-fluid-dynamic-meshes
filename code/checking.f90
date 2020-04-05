!checking
!
!comprueba que distintos procesos se ejecutan correctamente

module checking

   !modules used
   use precision, only : wp
   
   use TopologicalDataStructures
   !end moduless used
   
   implicit none
   
    contains
    
    !***************************************************************
    !checkTopology
    !
    !comprueba que se ha cargado bien la topologia
    subroutine checkTopology(meshT)
    
       !subroutine arguments
       type(topologymeshclass), intent(in) :: meshT
       !end subroutine arguments
       
       !local scalars
       integer :: i, j
       
       !end of header-----------------------------
       
       open(1891,file='checkTopology.dat')
       
       write(1891,'(a15,i8)') 'elements = ', meshT%Nelements
       write(1891,'(a)') 'element--- ifield--- ielement ---element type ---no of nodes--- no of edges--- nodes--- edges'
       do i=1,meshT%Nelements
           write(1891,'(3i8,4x,a25,15i8)') i, meshT%elements(i)%ifield, meshT%elements(i)%ielement, meshT%elements(i)%elementType, &
               meshT%elements(i)%Nnodes, meshT%elements(i)%Nedges, meshT%elements(i)%nodes(:), meshT%elements(i)%edges(:)
       end do
       
       write(1891,'(a15,i8)') 'edges = ', meshT%Nedges
       write(1891,'(a)') 'edge--- no of elements--- elements--- nodes'
       do i=1,meshT%Nedges
           write(1891,'(15i8)') i, meshT%edges(i)%Nelements, meshT%edges(i)%elements(:), meshT%edges(i)%nodes(:)
       end do
       
       write(1891,'(a15,i8)') 'nodes = ', meshT%Nnodes
       write(1891,'(a)') 'node--- no of edges--- no of elements--- edges--- elements'
       do i=1,meshT%Nnodes
           write(1891,'(25i8)') i, meshT%nodes(i)%Nedges, meshT%nodes(i)%Nelements, meshT%nodes(i)%edges(:), meshT%nodes(i)%elements(:)
       end do
       
       close(1891)
       
    end subroutine checkTopology
    
end module checking